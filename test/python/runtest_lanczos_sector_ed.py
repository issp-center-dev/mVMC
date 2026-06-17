from __future__ import print_function

import itertools
import math
import os
import shutil
import subprocess
import sys

import numpy as np


NSITE = 4


def read_out(filename):
    return np.loadtxt(filename, dtype="float").astype("float").flatten()


def occupied(state, orb):
    return (state >> orb) & 1


def popcount(value):
    return bin(value).count("1")


def annihilate(state, orb):
    if not occupied(state, orb):
        return None
    sign = -1 if popcount(state & ((1 << orb) - 1)) % 2 else 1
    return state & ~(1 << orb), sign


def create(state, orb):
    if occupied(state, orb):
        return None
    sign = -1 if popcount(state & ((1 << orb) - 1)) % 2 else 1
    return state | (1 << orb), sign


def apply_ops(state, ops):
    coeff = 1.0
    current = state
    for kind, orb in reversed(ops):
        if kind == "c":
            result = create(current, orb)
        else:
            result = annihilate(current, orb)
        if result is None:
            return None
        current, sign = result
        coeff *= sign
    return current, coeff


def fixed_basis(nup, ndown):
    basis = []
    for ups in itertools.combinations(range(NSITE), nup):
        for downs in itertools.combinations(range(NSITE), ndown):
            state = 0
            for site in ups:
                state |= 1 << site
            for site in downs:
                state |= 1 << (site + NSITE)
            basis.append(state)
    return basis


def doublon_allowed(state):
    for site in range(NSITE):
        if occupied(state, site) != occupied(state, site + NSITE):
            return False
    return True


def tj_allowed(state):
    for site in range(NSITE):
        if occupied(state, site) and occupied(state, site + NSITE):
            return False
    return True


def state_from_sites(ups, downs):
    state = 0
    for site in ups:
        state |= 1 << site
    for site in downs:
        state |= 1 << (site + NSITE)
    return state


def ap_orbital_amplitude(state, orbital):
    ups = [site for site in range(NSITE) if occupied(state, site)]
    downs = [site for site in range(NSITE) if occupied(state, site + NSITE)]
    if len(ups) != len(downs):
        return 0.0
    mat = np.array([[orbital[i, j] for j in downs] for i in ups], dtype=float)
    return np.linalg.det(mat)


def operator_matrix(basis, terms):
    index = dict((state, i) for i, state in enumerate(basis))
    mat = np.zeros((len(basis), len(basis)), dtype=float)
    for col, state in enumerate(basis):
        for coeff, ops in terms:
            result = apply_ops(state, ops)
            if result is None:
                continue
            next_state, sign = result
            if next_state in index:
                mat[index[next_state], col] += coeff * sign
    return mat


def one_body(ri, rj, spin):
    return [("c", ri + spin * NSITE), ("a", rj + spin * NSITE)]


def pair_hop(ri, rj):
    return [
        ("c", ri),
        ("a", rj),
        ("c", ri + NSITE),
        ("a", rj + NSITE),
    ]


def calculate_energy(h1, h2_1, h2_2, h3, h4):
    aa = h2_1 * (h2_1 + h2_2) - 2.0 * h1 * h3
    bb = -h1 * h2_1 + h3
    cc = (
        h2_1 * (h2_1 + h2_2) ** 2
        - h1 ** 2 * h2_1 * (h2_1 + 2.0 * h2_2)
        + 4.0 * h1 ** 3 * h3
        - 2.0 * h1 * (2.0 * h2_1 + h2_2) * h3
        + h3 ** 2
    )
    if cc < -1.0e-12:
        raise AssertionError("negative Lanczos discriminant: {}".format(cc))
    cc = max(cc, 0.0)
    candidates = []
    for alpha in ((bb + math.sqrt(cc)) / aa, (bb - math.sqrt(cc)) / aa):
        numerator = h1 + alpha * (h2_1 + h2_2) + alpha * alpha * h3
        dnorm = 1.0 + 2.0 * alpha * h1 + alpha * alpha * h2_1
        variance_numerator = h2_1 + 2.0 * alpha * h3 + alpha * alpha * h4
        energy = numerator / dnorm
        variance = ((variance_numerator / dnorm) - energy * energy) / (energy * energy)
        candidates.append((energy, variance, alpha))
    return min(candidates, key=lambda item: item[0])


def doublon_pairhop_ed_reference():
    basis = [state for state in fixed_basis(2, 2) if doublon_allowed(state)]
    orbital = np.full((NSITE, NSITE), 0.2, dtype=float)
    np.fill_diagonal(orbital, 1.0)
    psi = np.array([ap_orbital_amplitude(state, orbital) for state in basis], dtype=float)
    norm = float(np.dot(psi, psi))

    # mVMC's PairHop contribution follows the convention exercised by
    # CalculateHamiltonian: the listed pair-hop amplitude contributes twice to
    # the effective doublon hopping term for this AP-orbital sector.
    terms = []
    for ri, rj in ((0, 1), (1, 2), (2, 3), (3, 0)):
        terms.append((-2.0, pair_hop(ri, rj)))
    ham = operator_matrix(basis, terms)

    moments = []
    power = np.identity(len(basis), dtype=float)
    for _ in range(4):
        power = np.dot(ham, power)
        moments.append(float(np.dot(psi, np.dot(power, psi)) / norm))
    energy, variance, alpha = calculate_energy(
        moments[0], moments[1], moments[1], moments[2], moments[3]
    )
    return {
        "h1": moments[0],
        "h2": moments[1],
        "h3": moments[2],
        "h4": moments[3],
        "energy": energy,
        "variance": variance,
        "alpha": alpha,
    }


def assert_close(name, actual, expected, tol):
    diff = abs(actual - expected)
    if diff > tol:
        raise AssertionError(
            "{} mismatch: actual={} expected={} diff={} tol={}".format(
                name, actual, expected, diff, tol
            )
        )


def check_projected_operators():
    doublon_basis = [state for state in fixed_basis(2, 2) if doublon_allowed(state)]

    offdiag_transfer_terms = []
    for ri, rj in ((0, 1), (1, 0), (1, 2), (2, 1), (2, 3), (3, 2), (3, 0), (0, 3)):
        for spin in (0, 1):
            offdiag_transfer_terms.append((-1.0, one_body(ri, rj, spin)))
    offdiag_transfer = operator_matrix(doublon_basis, offdiag_transfer_terms)
    if np.linalg.norm(offdiag_transfer) != 0.0:
        raise AssertionError("doublon-only P K_offdiag P must be zero")

    diagonal_transfer_terms = []
    for site in range(NSITE):
        for spin in (0, 1):
            diagonal_transfer_terms.append((-1.0, one_body(site, site, spin)))
    diagonal_transfer = operator_matrix(doublon_basis, diagonal_transfer_terms)
    if np.linalg.norm(diagonal_transfer) == 0.0:
        raise AssertionError("doublon-only diagonal transfer must survive projection")

    pair_terms = []
    for ri, rj in ((0, 1), (1, 2), (2, 3), (3, 0)):
        pair_terms.append((-2.0, pair_hop(ri, rj)))
    pair_matrix = operator_matrix(doublon_basis, pair_terms)
    if np.linalg.norm(pair_matrix) == 0.0:
        raise AssertionError("doublon-only pair hopping must survive projection")

    start = state_from_sites((1, 3), (1, 3))
    break_doublon_bridge = one_body(2, 0, 0) + pair_hop(0, 1)
    raw = apply_ops(start, break_doublon_bridge)
    if raw is None or doublon_allowed(raw[0]):
        raise AssertionError("expected raw doublon-only bridge sequence to leave the sector")

    preserve_doublon_n3 = one_body(3, 3, 0) + pair_hop(0, 1)
    raw = apply_ops(start, preserve_doublon_n3)
    if raw is None or not doublon_allowed(raw[0]):
        raise AssertionError("expected sector-preserving GreenFuncN(3) sequence")

    preserve_doublon_n4 = pair_hop(2, 3) + pair_hop(0, 1)
    raw = apply_ops(start, preserve_doublon_n4)
    if raw is None or not doublon_allowed(raw[0]):
        raise AssertionError("expected sector-preserving GreenFuncN(4) sequence")

    tj_start = state_from_sites((1,), (3,))
    break_tj_bridge = one_body(0, 3, 1) + one_body(0, 1, 0)
    raw = apply_ops(tj_start, break_tj_bridge)
    if raw is None or tj_allowed(raw[0]):
        raise AssertionError("expected raw t-J bridge sequence to create doublon")

    preserve_tj_bridge = one_body(2, 3, 1) + one_body(0, 1, 0)
    raw = apply_ops(tj_start, preserve_tj_bridge)
    if raw is None or not tj_allowed(raw[0]):
        raise AssertionError("expected sector-preserving t-J GreenFuncN sequence")


def run_mvmc_doublon_only(rootdir):
    refdir = os.path.join(rootdir, "data", "DoublonOnly_Lanczos")
    workdir = os.path.join(rootdir, "work", "LanczosSectorED_DoublonOnly")
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    for name in os.listdir(refdir):
        if name.endswith(".def"):
            shutil.copy2(os.path.join(refdir, name), workdir)

    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    cmd = [bin_to_test, "-e", "namelist.def", "initial.def"]
    result = subprocess.call(cmd, cwd=workdir)
    if result != 0:
        raise AssertionError("vmc.out failed with exit code {}".format(result))

    qqqq = read_out(os.path.join(workdir, "output", "zvo_ls_qqqq_001.dat"))
    ls_out = read_out(os.path.join(workdir, "output", "zvo_ls_out_001.dat"))
    if not np.all(np.isfinite(qqqq)) or not np.all(np.isfinite(ls_out)):
        raise AssertionError("non-finite Lanczos ED verification output")
    return qqqq, ls_out


def check_doublon_only_against_mvmc(rootdir):
    reference = doublon_pairhop_ed_reference()
    qqqq, ls_out = run_mvmc_doublon_only(rootdir)

    assert_close("QQQQ[2] <H>", qqqq[2], reference["h1"], 0.20)
    assert_close("QQQQ[3] local <H^2>", qqqq[3], reference["h2"], 0.60)
    assert_close("QQQQ[10] <H^2>", qqqq[10], reference["h2"], 1.0e-10)
    assert_close("QQQQ[11] <H^3>", qqqq[11], reference["h3"], 0.80)
    assert_close("QQQQ[15] <H^4>", qqqq[15], reference["h4"], 1.0e-10)
    assert_close("Lanczos energy", ls_out[0], reference["energy"], 0.02)
    assert_close("Lanczos alpha", ls_out[2], reference["alpha"], 0.10)


def main():
    rootdir = os.getcwd()
    check_projected_operators()
    check_doublon_only_against_mvmc(rootdir)
    print("Lanczos sector ED verification passed")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AssertionError as exc:
        print("ERROR: {}".format(exc))
        sys.exit(-1)
