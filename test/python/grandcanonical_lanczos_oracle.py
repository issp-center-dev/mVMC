"""Independent Fock-space oracle for the grand-canonical first-step Lanczos.

The reference wave function is built from the pair matrix with Fock-space
operations only: the pair wave <x|exp(sum_{I<J} A_IJ c+_I c+_J)|0> = Pf(A_x)
(ascending orbitals, grandcanonical_spin_oracle), translations applied as
T_r|x> with the fermionic/boundary phase (grandcanonical_ap_oracle), finite
Gauss-Legendre S=0 projection as Fock creator rotations, the Re(Proj)
correlation factor counted from the occupation numbers, and the Phase 3b
boundary rule decided by an exact rational sum of the raw translation
weights.  H is applied with fermion operators from the input rows
(Trans: -t c+_out c_in, CoulombIntra U n_up n_dn, CoulombInter V n_i n_j,
Hund -J (n_iu n_ju + n_id n_jd)).  q0 = Psi, q1 = H Psi, q2 = H^2 Psi, the
2x2 S/K/T matrices and the GEVP are evaluated with numpy.

production_qp() additionally writes the production-shaped QP Slater
matrices (translated, 4-block rotated) that the C probe consumes; they are
never used to build the reference.
"""
from __future__ import print_function

import cmath
import math
from fractions import Fraction

import numpy as np

import grandcanonical_ap_oracle as ap
import grandcanonical_exact_oracle as ex
import grandcanonical_spin_oracle as so

RANK_CUTOFF = 2.0 ** -40


def popcount(value):
    return bin(value).count("1")


def domain_basis(nsite, domain):
    masks = [m for m in range(1 << (2 * nsite)) if popcount(m) % 2 == 0]
    if domain == "balanced":
        low = (1 << nsite) - 1
        masks = [m for m in masks
                 if popcount(m & low) == popcount(m >> nsite)]
    return masks


def chain_mappings(nsite, boundary, steps):
    """(r, site) -> (target, phase) for translations by r sites."""
    mappings = {}
    for r in steps:
        for site in range(nsite):
            target = site + r
            wraps = target // nsite
            phase = -1 if (boundary == "apbc" and wraps % 2) else 1
            mappings[(r, site)] = (target % nsite, phase)
    return mappings


def antiparallel_pair_matrix(F):
    return so.antiparallel_matrix(F)


def fixture_pair_matrix(fixture):
    if fixture["pairing"] == "antiparallel":
        return antiparallel_pair_matrix(fixture["F"])
    return [list(row) for row in fixture["A"]]


def exact_sum_is_zero(values):
    real = sum(Fraction(complex(v).real) for v in values)
    imag = sum(Fraction(complex(v).imag) for v in values)
    return real == 0 and imag == 0


def log_projection(fixture, mask):
    nsite = fixture["L"]
    up = [(mask >> s) & 1 for s in range(nsite)]
    down = [(mask >> (s + nsite)) & 1 for s in range(nsite)]
    value = 0.0
    params = fixture.get("projection_parameters", [])
    gutz = fixture.get("gutzwiller_idx")
    ngutz = fixture.get("ngutzwiller", 0)
    if gutz is not None:
        for s in range(nsite):
            value += complex(params[gutz[s]]).real * up[s] * down[s]
    table = fixture.get("jastrow_idx")
    if table is not None:
        for ri in range(nsite):
            xi = up[ri] + down[ri] - 1
            for rj in range(ri + 1, nsite):
                xj = up[rj] + down[rj] - 1
                value += complex(params[ngutz + table[ri][rj]]).real * xi * xj
    return value


def translated_wave(base, mappings, r, nsite):
    """phi_{A_r}(x) = phase(T_r, x) phi_A(T_r x) from Fock translations."""
    result = {}
    for mask in base:
        target, phase = ap.translate_state(mask, mappings, r, nsite)
        result[mask] = phase * base.get(target, 0.0j)
    return result


def pure_triplet_s0(fixture):
    return (fixture["pairing"] == "antiparallel" and
            fixture.get("gauss", 1) > 1 and fixture.get("S", 0) == 0 and
            all(complex(fixture["F"][i][j]) == -complex(fixture["F"][j][i])
                for i in range(fixture["L"]) for j in range(fixture["L"])))


def pair_symmetry(fixture):
    """Exact class of the unrotated F: S=0 reads only the triplet test
    (ZERO -> triplet, SINGLET -> general); S>0 uses the four classes."""
    n, F = fixture["L"], fixture["F"]
    triplet = all(complex(F[i][j]) == -complex(F[j][i])
                  for i in range(n) for j in range(n))
    if fixture.get("S", 0) == 0:
        return "triplet" if triplet else "general"
    singlet = all(complex(F[i][j]) == complex(F[j][i])
                  for i in range(n) for j in range(n))
    if triplet and singlet:
        return "zero"
    return "triplet" if triplet else "singlet" if singlet else "general"


def support_zero(fixture, mask):
    """Structural zero of the base Psi (never applied to H Psi or Phi)."""
    if fixture["pairing"] != "antiparallel" or fixture.get("gauss", 1) <= 1:
        return False
    return so.finite_spin_zero(mask, fixture["L"], fixture.get("S", 0),
                               pair_symmetry(fixture))


def reference_wave(fixture):
    """Psi(x) over the domain basis from Fock operations only."""
    nsite = fixture["L"]
    steps = [r for r, _ in fixture["translations"]]
    mappings = chain_mappings(nsite, fixture["boundary"], steps)
    if fixture["pairing"] == "antiparallel":
        base = so.pair_wave(fixture["F"])
    else:
        values = so.float_pfaffians(fixture_pair_matrix(fixture))
        base = dict((m, values.get(m, 0.0j))
                    for m in range(1 << (2 * nsite)) if popcount(m) % 2 == 0)
    gauss = fixture.get("gauss", 1)
    spin = fixture.get("S", 0)
    nodes = so.independent_spin_nodes(gauss, spin) if gauss > 1 else None
    total = dict((m, 0.0j) for m in base)
    for r, weight in fixture["translations"]:
        wave = translated_wave(base, mappings, r, nsite)
        if nodes is None:
            for mask, value in wave.items():
                total[mask] += complex(weight) * value
        else:
            for node in nodes:
                c, s = so.node_cos_sin(node)
                rotated = so.rotate_fock(wave, nsite, c, s)
                for mask in total:
                    total[mask] += complex(weight) * so.node_weight(node) * \
                        rotated.get(mask, 0.0j)
    boundary_zero = exact_sum_is_zero(
        [w for _, w in fixture.get("raw_translation_weights",
                                   fixture["translations"])])
    basis = domain_basis(nsite, fixture["domain"])
    psi = np.zeros(len(basis), dtype=complex)
    for k, mask in enumerate(basis):
        if support_zero(fixture, mask):
            continue
        if boundary_zero and popcount(mask) in (0, 2 * nsite):
            psi[k] = 0.0
            continue
        psi[k] = total.get(mask, 0.0j) * math.exp(log_projection(fixture, mask))
    return basis, psi


def production_qp(fixture):
    """Production-shaped QP Slater matrices (row-major) and QP weights."""
    nsite = fixture["L"]
    steps = [r for r, _ in fixture["translations"]]
    mappings = chain_mappings(nsite, fixture["boundary"], steps)
    base = fixture_pair_matrix(fixture)
    gauss = fixture.get("gauss", 1)
    nodes = so.independent_spin_nodes(gauss, fixture.get("S", 0))
    slaters, weights = [], []
    for r, weight in fixture["translations"]:
        translated = ap.translated_matrix(base, mappings, r, nsite)
        if fixture["pairing"] == "general":
            if gauss != 1:
                raise ValueError("general pairing has no spin projection")
            slaters.append(translated)
            weights.append(complex(weight))
            continue
        F_r = [row[nsite:] for row in translated[:nsite]]
        for node in nodes:
            slaters.append(so.rotated_block_matrix(
                F_r, node, convert=complex, zero=0.0j,
                mul=lambda a, b: a * b, sub=lambda a, b: a - b,
                add=lambda a, b: a + b, neg=lambda a: -a))
            weights.append(complex(weight) * so.node_weight(node))
    return slaters, weights


def hamiltonian(fixture, basis):
    nsite = fixture["L"]
    index = dict((m, k) for k, m in enumerate(basis))
    h = np.zeros((len(basis), len(basis)), dtype=complex)

    def add(operators, coefficient):
        for col, mask in enumerate(basis):
            out = ex.apply_ops(mask, operators)
            if out is None:
                continue
            if out[0] not in index:
                raise AssertionError("H leaves the configuration domain")
            h[index[out[0]], col] += coefficient * out[1]

    for out_site, out_spin, in_site, in_spin, value in fixture["trans"]:
        add(ex.one_body(out_site + out_spin * nsite, in_site + in_spin * nsite),
            -complex(value))
    for site, value in fixture.get("intra", []):
        add(ex.nbody(((site, site), (site + nsite, site + nsite))), value)
    for i, j, value in fixture.get("inter", []):
        for si in range(2):
            for sj in range(2):
                add(ex.nbody(((i + si * nsite, i + si * nsite),
                              (j + sj * nsite, j + sj * nsite))), value)
    for i, j, value in fixture.get("hund", []):
        for spin in range(2):
            add(ex.nbody(((i + spin * nsite, i + spin * nsite),
                          (j + spin * nsite, j + spin * nsite))), -value)
    if np.max(np.abs(h - h.conj().T), initial=0.0) > 1e-13:
        raise AssertionError("oracle H is not Hermitian")
    return h


def krylov(psi, h):
    q = [psi, h @ psi]
    q.append(h @ q[1])
    s = np.array([[np.vdot(q[i], q[j]) for j in range(2)] for i in range(2)])
    k = np.array([[np.vdot(q[i], q[j + 1]) for j in range(2)]
                  for i in range(2)])
    t = np.array([[np.vdot(q[i + 1], q[j + 1]) for j in range(2)]
                  for i in range(2)])
    return q, s, k, t


def gevp(s, k, t):
    """Independent whitening + eigh; None for rank 0."""
    w, v = np.linalg.eigh((s + s.conj().T) / 2)
    if not np.all(np.isfinite(w)) or w.max(initial=0.0) <= 0.0:
        return None
    keep = w > w.max() * RANK_CUTOFF
    if not np.any(keep):
        return None
    x = v[:, keep] / np.sqrt(w[keep])
    ks = (k + k.conj().T) / 2
    e, y = np.linalg.eigh(x.conj().T @ ks @ x)
    c = x @ y[:, 0]
    norm = np.vdot(c, s @ c).real
    energy = (np.vdot(c, ks @ c) / norm).real
    second = (np.vdot(c, ((t + t.conj().T) / 2) @ c) / norm).real
    gap = e[1] - e[0] if len(e) > 1 else math.inf
    return {"c": c, "energy": energy, "variance": second - energy * energy,
            "rank": int(np.count_nonzero(keep)), "roots": e, "gap": gap}


def packed(matrix):
    return [matrix[0, 0], matrix[0, 1], matrix[1, 1]]


def self_test():
    # Plan anchor: two-state diagonal H with psi = (1, 2).
    q = [np.array([1.0, 2.0])]
    h = np.diag([0.0, 3.0])
    q, s, k, t = krylov(q[0], h)
    np.testing.assert_array_equal(s, [[5.0, 12.0], [12.0, 36.0]])
    np.testing.assert_array_equal(k, [[12.0, 36.0], [36.0, 108.0]])
    np.testing.assert_array_equal(t, [[36.0, 108.0], [108.0, 324.0]])
    # Node fixture: F=[[1,0],[1,1]], mask 0x9 has q0 = 0 and q1 = -2.
    node = {"L": 2, "pairing": "antiparallel", "F": [[1.0, 0.0], [1.0, 1.0]],
            "domain": "balanced", "boundary": "pbc",
            "translations": [(0, 1.0)],
            "trans": [(0, 0, 1, 0, 1.0), (1, 0, 0, 0, 1.0),
                      (0, 1, 1, 1, 1.0), (1, 1, 0, 1, 1.0)]}
    basis, psi = reference_wave(node)
    q, s, k, t = krylov(psi, hamiltonian(node, basis))
    position = basis.index(0x9)
    assert q[0][position] == 0.0 and abs(q[1][position] + 2.0) < 1e-15
    # Fock translation equals the production-shaped translated Pfaffian.
    A = [[0.0j] * 4 for _ in range(4)]
    F = [[1.0, (1 + 2j) / 8], [(3 - 1j) / 16, 0.75]]
    A = so.antiparallel_matrix(F)
    A[0][1], A[1][0] = (1 + 1j) / 4, -(1 + 1j) / 4
    A[2][3], A[3][2] = (2 - 1j) / 8, -(2 - 1j) / 8
    mappings = chain_mappings(2, "apbc", [1])
    base = so.float_pfaffians(A)
    translated = so.float_pfaffians(ap.translated_matrix(A, mappings, 1, 2))
    fock = translated_wave(dict((m, base[m]) for m in base), mappings, 1, 2)
    for mask in base:
        assert abs(fock[mask] - translated[mask]) < 1e-14, mask
    # Exact boundary sum.
    assert exact_sum_is_zero([1.0, -1.0]) and \
        not exact_sum_is_zero([1.0, 2.0 ** -53, -1.0])
    return True


if __name__ == "__main__":
    print("grandcanonical_lanczos_oracle self-test:",
          "PASS" if self_test() else "FAIL")
