"""Independent Fock-space oracle for anti-parallel grand-canonical pairing.

The pair wave function is built by expanding
    |Phi(F)> = exp(sum_ij F_ij c+_{i up} c+_{j down}) |0>
with explicit fermion operators.  Nothing here reuses the mVMC readers,
writers or the OrbitalGeneral 2*f construction, so production results can be
compared against amplitudes and local estimators defined from first
principles.

Orbital numbering: up orbital = site, down orbital = site + nsite.  A mask
holds bit k when orbital k is occupied.  Operator tuples are written left to
right and act on a ket from right to left; "c" creates and "a" annihilates.
"""
from __future__ import print_function

import cmath
import itertools
import math
import sys

import numpy as np


NSITE = 4
GUTZ = 0.08
JASTROW = 0.11
F0 = np.array([
    [.71 + .12j, .34 - .23j, -.29 + .17j, .18 + .31j],
    [-.41 + .26j, .62 - .11j, .27 + .22j, -.33 - .19j],
    [.23 + .37j, -.36 + .14j, .58 + .09j, .31 - .28j],
    [.19 - .24j, .28 + .33j, -.42 + .16j, .67 - .18j],
])
HOPPING = 0.37 + 0.19j
CHEMICAL = 0.17
ONSITE = (0.73, -0.28, 0.41, 0.62)
DENSITY = 0.21
HUND = 0.17
PAIRHOP = 0.13
EXCHANGE = -0.11
INTERALL = -0.16 + 0.09j
NBODY = 0.07
ANOMALOUS = 0.35 - 0.20j

# Reference values printed in the implementation plan (Appendix A).
REFERENCE = {
    False: {"energy": 0.7233489249304274, "energy2": 2.034616794887211,
            "number": 3.1987781458664806, "number2": 14.734215976089576,
            "anomalous": 0.3448501149201095 - 0.05623804890252127j},
    True: {"energy": 0.9450559261818133, "energy2": 2.5389804230511452,
           "number": 3.0950691723876265, "number2": 13.763135351762099,
           "anomalous": 0.3348044990060462 - 0.04937183008980556j},
}
ANOMALOUS_KEY = (1, 0, 0, 0, 1)


def popcount(value):
    return bin(value).count("1")


def finite_close(actual, expected, tolerance):
    try:
        values = (complex(actual), complex(expected))
        tolerance = float(tolerance)
    except (TypeError, ValueError):
        return False
    return (math.isfinite(tolerance) and tolerance >= 0 and
            all(math.isfinite(z.real) and math.isfinite(z.imag)
                for z in values) and
            abs(values[0] - values[1]) <= tolerance)


def apply_ops(mask, operators):
    sign = 1
    for kind, orbital in reversed(operators):
        occupied = (mask >> orbital) & 1
        if occupied == (kind == "c"):
            return None
        sign *= -1 if popcount(mask & ((1 << orbital) - 1)) % 2 else 1
        mask ^= 1 << orbital
    return mask, sign


def one_body(out_orbital, in_orbital):
    return (("c", out_orbital), ("a", in_orbital))


def product(*factors):
    operators = []
    for factor in factors:
        operators.extend(factor)
    return tuple(operators)


def adjoint(operators):
    return tuple(("a" if kind == "c" else "c", orbital)
                 for kind, orbital in reversed(operators))


def ap_signs(ap, nsite=NSITE):
    signs = np.ones((nsite, nsite))
    if ap:
        signs[0, nsite - 1] = signs[nsite - 1, 0] = -1.0
    return signs


def fixture_matrix(ap):
    return F0 * ap_signs(ap)


def pair_wave(F):
    L = len(F)
    result = dict((mask, 0j) for mask in range(1 << (2 * L)))
    term = {0: 1 + 0j}
    result[0] = 1 + 0j
    for order in range(1, L + 1):
        nxt = {}
        for mask, value in term.items():
            for i in range(L):
                for j in range(L):
                    out = apply_ops(mask, (("c", i), ("c", j + L)))
                    if out is not None:
                        target, sign = out
                        nxt[target] = (nxt.get(target, 0j) +
                                       value * F[i, j] * sign / order)
        for mask, value in nxt.items():
            result[mask] += value
        term = nxt
    return result


def occupations(mask, nsite):
    up = [i for i in range(nsite) if (mask >> i) & 1]
    down = [i for i in range(nsite) if (mask >> (i + nsite)) & 1]
    return up, down


def is_balanced(mask, nsite=NSITE):
    up, down = occupations(mask, nsite)
    return len(up) == len(down)


def ring_bonds(nsite):
    if nsite == 2:
        return ((0, 1),)
    return tuple((i, (i + 1) % nsite) for i in range(nsite))


def projection_counts(mask, nsite=NSITE):
    """Gutzwiller, ring-Jastrow and diagonal-Jastrow counts as in MakeProjCnt."""
    up, down = occupations(mask, nsite)
    n = [int(i in up) + int(i in down) for i in range(nsite)]
    doublon = sum(int(i in up and i in down) for i in range(nsite))
    ring = set(tuple(sorted(bond)) for bond in ring_bonds(nsite))
    ring_count = 0
    diagonal_count = 0
    for i, j in itertools.combinations(range(nsite), 2):
        value = (n[i] - 1) * (n[j] - 1)
        if (i, j) in ring:
            ring_count += value
        else:
            diagonal_count += value
    return doublon, ring_count, diagonal_count


def wave(F, gutz=GUTZ, jastrow=JASTROW):
    nsite = len(F)
    result = pair_wave(F)
    for mask in result:
        doublon, ring_count, unused_diagonal = projection_counts(mask, nsite)
        result[mask] *= math.exp(gutz * doublon + jastrow * ring_count)
    return result


def model_terms(ap, nsite=NSITE):
    if nsite != NSITE:
        raise ValueError("the fixed model is defined for four sites")
    L = nsite
    terms = []
    for orbital in range(2 * L):
        terms.append(("transfer", -CHEMICAL, one_body(orbital, orbital)))
    for i, j in ring_bonds(L):
        t = -HOPPING if ap and (i, j) == (L - 1, 0) else HOPPING
        for s in (0, L):
            terms.append(("transfer", -t, one_body(i + s, j + s)))
            terms.append(("transfer", -t.conjugate(), one_body(j + s, i + s)))
    for i, U in enumerate(ONSITE):
        terms.append(("onsite", U, product(one_body(i, i),
                                            one_body(i + L, i + L))))
    for i, j in ring_bonds(L):
        for a, b in itertools.product((0, L), repeat=2):
            terms.append(("density", DENSITY,
                          product(one_body(i + a, i + a),
                                  one_body(j + b, j + b))))
    for s in (0, L):
        terms.append(("hund", -HUND, product(one_body(s, s),
                                             one_body(1 + s, 1 + s))))
    paired = (
        ("pairhop", PAIRHOP, product(one_body(0, 1), one_body(L, L + 1))),
        ("exchange", EXCHANGE, product(one_body(0, 1), one_body(L + 1, L))),
        ("interall", INTERALL, product(one_body(0, L), one_body(L + 1, 1))),
        ("nbody", NBODY, product(one_body(0, L), one_body(L + 1, 1),
                                 one_body(2, 2))),
        ("anomalous", ANOMALOUS, (("c", 0), ("c", L))),
    )
    for label, coefficient, operators in paired:
        terms.append((label, complex(coefficient), operators))
        terms.append((label, complex(coefficient).conjugate(),
                      adjoint(operators)))
    return terms


def operator_matrix(operators, nsite=NSITE):
    dimension = 1 << (2 * nsite)
    matrix = np.zeros((dimension, dimension), dtype=np.complex128)
    for source in range(dimension):
        out = apply_ops(source, operators)
        if out is not None:
            target, sign = out
            matrix[target, source] += sign
    return matrix


def hamiltonian_from_terms(terms, nsite=NSITE):
    dimension = 1 << (2 * nsite)
    matrix = np.zeros((dimension, dimension), dtype=np.complex128)
    for unused_label, coefficient, operators in terms:
        for source in range(dimension):
            out = apply_ops(source, operators)
            if out is not None:
                target, sign = out
                matrix[target, source] += coefficient * sign
    return matrix


def sz_twice(mask, nsite=NSITE):
    up, down = occupations(mask, nsite)
    return len(up) - len(down)


def onebody_operators(nsite=NSITE):
    return dict(((out_orbital, in_orbital), one_body(out_orbital, in_orbital))
                for out_orbital in range(2 * nsite)
                for in_orbital in range(2 * nsite))


def anomalous_operator(key, nsite=NSITE):
    kind, s1, sp1, s2, sp2 = key
    first = s1 + sp1 * nsite
    second = s2 + sp2 * nsite
    letter = "c" if kind == 1 else "a"
    return ((letter, first), (letter, second))


def anomalous_keys(nsite=NSITE):
    keys = []
    for kind in (1, 0):
        for s1, sp1, s2, sp2 in itertools.product(range(nsite), (0, 1),
                                                  range(nsite), (0, 1)):
            if (s1, sp1) != (s2, sp2):
                keys.append((kind, s1, sp1, s2, sp2))
    return keys


def wave_vector(F, gutz=GUTZ, jastrow=JASTROW):
    values = wave(F, gutz, jastrow)
    return np.array([values[mask] for mask in range(len(values))])


def local_estimator(matrix, psi, masks):
    """Return sum_y conj(psi_y) O_yx / conj(psi_x) for each sampled x."""
    numerator = matrix.T @ psi.conj()
    result = {}
    for mask in masks:
        if psi[mask] == 0:
            raise ValueError("local estimator requested on a zero amplitude")
        result[mask] = numerator[mask] / psi[mask].conjugate()
    return result


def orbital_log_derivative(F, signs, mask, nsite):
    """d ln psi / d f for raw parameters f with F_ij = signs_ij f_ij."""
    up, down = occupations(mask, nsite)
    derivative = np.zeros((nsite, nsite), dtype=np.complex128)
    if up:
        inverse = np.linalg.inv(F[np.ix_(up, down)])
        for a, i in enumerate(up):
            for b, j in enumerate(down):
                derivative[i, j] = inverse[b, a] * signs[i, j]
    return derivative


def local_values(F, ap, masks=None, gutz=GUTZ, jastrow=JASTROW,
                 onebody_keys=None, anomalous_list=None):
    nsite = len(F)
    psi = wave_vector(F, gutz, jastrow)
    if masks is None:
        masks = [mask for mask in range(len(psi)) if psi[mask] != 0]
    masks = sorted(set(masks))
    for mask in masks:
        if not (psi[mask] != 0 and np.isfinite(psi[mask])):
            raise ValueError("mask {} has no positive probability".format(mask))
    hamiltonian = hamiltonian_from_terms(model_terms(ap, nsite), nsite)
    energy = local_estimator(hamiltonian, psi, masks)
    signs = ap_signs(ap, nsite)
    if onebody_keys is None:
        onebody_keys = list(onebody_operators(nsite))
    if anomalous_list is None:
        anomalous_list = anomalous_keys(nsite)
    onebody = dict((mask, {}) for mask in masks)
    for key in onebody_keys:
        values = local_estimator(operator_matrix(one_body(*key), nsite),
                                 psi, masks)
        for mask in masks:
            onebody[mask][key] = values[mask]
    anomalous = dict((mask, {}) for mask in masks)
    for key in anomalous_list:
        values = local_estimator(
            operator_matrix(anomalous_operator(key, nsite), nsite), psi, masks)
        for mask in masks:
            anomalous[mask][key] = values[mask]
    derivative = {}
    projection = {}
    for mask in masks:
        base = orbital_log_derivative(F, signs, mask, nsite)
        table = {}
        for i in range(nsite):
            for j in range(nsite):
                table[i, j, False] = base[i, j]
                table[i, j, True] = 1j * base[i, j]
        derivative[mask] = table
        projection[mask] = projection_counts(mask, nsite)
    return {
        "masks": masks,
        "energy": energy,
        "energy2": dict((mask, abs(energy[mask]) ** 2) for mask in masks),
        "number": dict((mask, popcount(mask)) for mask in masks),
        "number2": dict((mask, popcount(mask) ** 2) for mask in masks),
        "onebody": onebody,
        "anomalous": anomalous,
        "derivative": derivative,
        "projection": projection,
    }


def expectation(psi, matrix):
    norm = np.vdot(psi, psi).real
    return np.vdot(psi, matrix @ psi) / norm


def exact(F, ap, gutz=GUTZ, jastrow=JASTROW):
    nsite = len(F)
    psi = wave_vector(F, gutz, jastrow)
    norm = np.vdot(psi, psi).real
    probability = abs(psi) ** 2 / norm
    hamiltonian = hamiltonian_from_terms(model_terms(ap, nsite), nsite)
    action = hamiltonian @ psi
    numbers = np.array([popcount(mask) for mask in range(len(psi))],
                       dtype=float)
    onebody = dict((key, expectation(psi, operator_matrix(ops, nsite)))
                   for key, ops in onebody_operators(nsite).items())
    anomalous = dict(
        (key, expectation(psi, operator_matrix(anomalous_operator(key, nsite),
                                               nsite)))
        for key in anomalous_keys(nsite))
    support = [mask for mask in range(len(psi)) if probability[mask] > 0]
    local = local_values(F, ap, support, gutz, jastrow, onebody_keys=[],
                         anomalous_list=[])
    mean_e = sum(probability[m] * local["energy"][m] for m in support)
    gradient = {}
    for key in local["derivative"][support[0]]:
        mean_o = sum(probability[m] * local["derivative"][m][key]
                     for m in support)
        mean_oh = sum(probability[m] * local["energy"][m] *
                      local["derivative"][m][key] for m in support)
        gradient[key] = 2.0 * (mean_oh.real - mean_e.real * mean_o.real)
    return {
        "wave": psi,
        "probability": probability,
        "energy": np.vdot(psi, action) / norm,
        "energy2": np.vdot(action, action).real / norm,
        "number": float(np.dot(probability, numbers)),
        "number2": float(np.dot(probability, numbers ** 2)),
        "onebody": onebody,
        "anomalous": anomalous,
        "gradient": gradient,
    }


def pfaffian(matrix):
    n = len(matrix)
    if n == 0:
        return 1.0 + 0.0j
    if n % 2:
        return 0.0j
    value = 0.0j
    for column in range(1, n):
        keep = [k for k in range(n) if k not in (0, column)]
        minor = [[matrix[r][c] for c in keep] for r in keep]
        value += (-1) ** (column + 1) * matrix[0][column] * pfaffian(minor)
    return value


def antisymmetric_matrix(F):
    L = len(F)
    A = np.zeros((2 * L, 2 * L), dtype=np.complex128)
    A[:L, L:] = F
    A[L:, :L] = -F.T
    return A


def general_upper_parameters(F):
    """OrbitalGeneral upper-triangle parameters equal to the same state."""
    L = len(F)
    parameters = {}
    for first in range(2 * L):
        for second in range(first + 1, 2 * L):
            if first < L <= second:
                parameters[first, second] = F[first, second - L] / 2.0
            else:
                parameters[first, second] = 0.0j
    return parameters


def general_matrix(parameters, nsite):
    # OrbitalGeneral stores f_IJ=p and f_JI=-p, then uses A_IJ=f_IJ-f_JI.
    f = np.zeros((2 * nsite, 2 * nsite), dtype=np.complex128)
    for (first, second), value in parameters.items():
        f[first, second] = value
        f[second, first] = -value
    return f - f.T


def permutation_sign(order):
    sign = 1
    order = list(order)
    for i in range(len(order)):
        for j in range(i + 1, len(order)):
            if order[i] > order[j]:
                sign = -sign
    return sign


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def self_test():
    F = np.array([[0.7 + 0.2j, 0.3 - 0.1j], [-0.4 + 0.3j, 0.6 + 0.1j]])
    w = pair_wave(F)
    check(finite_close(w[0], 1, 1e-13), "vacuum amplitude")
    check(finite_close(w[(1 << 0) | (1 << 2)], F[0, 0], 1e-13),
          "one-pair amplitude")
    check(finite_close(w[15], -np.linalg.det(F), 1e-13),
          "full-filling amplitude")
    check(w[3] == 0 and w[12] == 0, "same-spin pair amplitude")
    check(not finite_close(float("nan"), 0, 1), "NaN accepted")
    check(not finite_close(0, 0, float("inf")), "infinite tolerance accepted")
    check(not finite_close(complex(0, float("inf")), 0, 1),
          "infinite imaginary part accepted")
    check(not finite_close(0, 0, -1.0), "negative tolerance accepted")
    check(not finite_close(0, float("-inf"), 1), "infinite expected accepted")

    for ap in (False, True):
        Fx = fixture_matrix(ap)
        L = len(Fx)
        amplitudes = pair_wave(Fx)
        A = antisymmetric_matrix(Fx)
        general = general_matrix(general_upper_parameters(Fx), L)
        check(np.max(abs(general - A)) < 1e-15, "General F/2 embedding")
        balanced = [mask for mask in amplitudes if is_balanced(mask, L)]
        check(len(balanced) == 70, "balanced basis size")
        sector_norm = {}
        for mask, amplitude in amplitudes.items():
            up, down = occupations(mask, L)
            if len(up) != len(down):
                check(amplitude == 0, "amplitude outside Sz=0")
                continue
            m = len(up)
            expected = ((-1) ** (m * (m - 1) // 2) *
                        (np.linalg.det(Fx[np.ix_(up, down)]) if m else 1.0))
            check(finite_close(amplitude, expected, 2e-12),
                  "minor determinant mismatch mask={}".format(mask))
            check(abs(amplitude) > 1e-6, "zero amplitude in main fixture")
            orbitals = up + [i + L for i in down]
            pf = pfaffian(A[np.ix_(orbitals, orbitals)].tolist())
            check(finite_close(pf, amplitude, 2e-12), "Pfaffian mismatch")
            pf_general = pfaffian(general[np.ix_(orbitals, orbitals)].tolist())
            check(finite_close(pf_general, amplitude, 2e-12),
                  "General F/2 Pfaffian mismatch")
            if len(orbitals) >= 2:
                order = orbitals[1:] + orbitals[:1]
                pf_order = pfaffian(A[np.ix_(order, order)].tolist())
                check(finite_close(pf_order,
                                   permutation_sign(order) * amplitude,
                                   2e-12), "electron-order sign mismatch")
            sector_norm[2 * m] = sector_norm.get(2 * m, 0.0) + abs(amplitude) ** 2
        check(finite_close(amplitudes[0], 1, 0), "vacuum normalization")
        check(finite_close(amplitudes[(1 << (2 * L)) - 1],
                           (-1) ** (L * (L - 1) // 2) * np.linalg.det(Fx),
                           2e-12), "full filling")
        check(sorted(sector_norm) == [0, 2, 4, 6, 8], "particle sectors")
        check(all(value > 0 for value in sector_norm.values()), "sector norm")

        H = hamiltonian_from_terms(model_terms(ap))
        check(H.shape == (256, 256), "Hamiltonian shape")
        check(np.max(abs(H - H.conj().T)) < 1e-14, "Hamiltonian is not Hermitian")
        for target, source in zip(*np.nonzero(H)):
            check(sz_twice(int(target)) == sz_twice(int(source)),
                  "Hamiltonian changes Sz")
        values = exact(Fx, ap)
        reference = REFERENCE[ap]
        check(finite_close(values["energy"], reference["energy"], 1e-12),
              "energy reference ap={}".format(ap))
        check(finite_close(values["energy2"], reference["energy2"], 1e-12),
              "energy2 reference")
        check(finite_close(values["number"], reference["number"], 1e-12),
              "number reference")
        check(finite_close(values["number2"], reference["number2"], 1e-12),
              "number2 reference")
        check(finite_close(values["anomalous"][ANOMALOUS_KEY],
                           reference["anomalous"], 1e-12),
              "anomalous reference")
        check(finite_close(values["anomalous"][(1, 0, 1, 0, 0)],
                           -reference["anomalous"], 1e-12),
              "reversed anomalous sign")
        check(abs(values["anomalous"][(1, 0, 0, 1, 0)]) < 1e-14,
              "same-spin pair expectation")
        check(abs(values["onebody"][(0, 4)]) < 1e-14, "spin flip expectation")

        probability = values["probability"]
        support = [m for m in range(256) if probability[m] > 0]
        check(len(support) == 70, "support size")
        local = local_values(Fx, ap, support)
        energy = sum(probability[m] * local["energy"][m] for m in support)
        check(finite_close(energy, values["energy"], 1e-12),
              "local energy average")
        energy2 = sum(probability[m] * local["energy2"][m] for m in support)
        check(finite_close(energy2, values["energy2"], 1e-12),
              "local energy2 average (no nodes in Sz=0)")
        for key in ((0, 0), (0, 5), (4, 1), (7, 3)):
            mean = sum(probability[m] * local["onebody"][m][key]
                       for m in support)
            check(finite_close(mean, values["onebody"][key], 1e-12),
                  "one-body local average {}".format(key))
        mean = sum(probability[m] * local["anomalous"][m][ANOMALOUS_KEY]
                   for m in support)
        check(finite_close(mean, values["anomalous"][ANOMALOUS_KEY], 1e-12),
              "anomalous local average")

        epsilon = 1e-6
        signs = ap_signs(ap)
        for (i, j) in ((0, 0), (0, 3), (2, 2), (3, 1)):
            for imaginary in (False, True):
                step = (1j if imaginary else 1.0) * epsilon * signs[i, j]
                plus = Fx.copy()
                minus = Fx.copy()
                plus[i, j] += step
                minus[i, j] -= step
                wave_plus = wave_vector(plus)
                wave_minus = wave_vector(minus)
                psi = values["wave"]
                for mask in support:
                    numeric = ((wave_plus[mask] - wave_minus[mask]) /
                               (2 * epsilon * psi[mask]))
                    analytic = local["derivative"][mask][i, j, imaginary]
                    check(finite_close(analytic, numeric,
                                       2e-8 * (1 + abs(numeric))),
                          "log derivative mismatch {} {}".format((i, j), mask))
        check(all(local["derivative"][0][key] == 0
                  for key in local["derivative"][0]),
              "vacuum orbital derivative")

    rank_one = np.outer([0.5 + 0.1j, -0.3 + 0.2j, 0.4, 0.2 - 0.6j],
                        [0.3 - 0.2j, 0.7, -0.1 + 0.4j, 0.5 + 0.5j])
    low = pair_wave(rank_one)
    for mask, amplitude in low.items():
        up, down = occupations(mask, NSITE)
        if len(up) != len(down):
            check(amplitude == 0, "low-rank amplitude outside Sz=0")
        elif len(up) >= 2:
            check(abs(amplitude) < 1e-15, "low-rank node is not zero")
        else:
            check(abs(amplitude) > 0, "low-rank support missing")
    return exact(fixture_matrix(False), False)


if __name__ == "__main__":
    try:
        result = self_test()
        print("anti-parallel GC oracle passed: E={:.15g} N={:.15g}".format(
            result["energy"].real, result["number"]))
    except Exception as error:  # pragma: no cover - CTest reports the text
        print("ERROR: {}".format(error))
        sys.exit(1)
