from __future__ import print_function

import itertools
import math
import os
import shutil
import subprocess
import sys

import numpy as np

from runtest_bf_fsz import (
    append_namelist_entry,
    prepare_case,
    run_vmc,
    update_modpara,
    write_lanczos_init,
    write_spin_changing_c1_defs,
)


NSITE = 4
COULOMB_INTRA = (0.35, 0.70, 1.05, 1.40)
PROJECTION_BRANCH_GUTZWILLER = 1.0
PROJECTION_BRANCH_RAW_RATIO = 5.0e-13
PROJECTION_BRANCH_BASE_STATE = 75
PROJECTION_BRANCH_TARGET_STATE = 89


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


def operator_matrix_complex(basis, terms):
    index = dict((state, i) for i, state in enumerate(basis))
    mat = np.zeros((len(basis), len(basis)), dtype=complex)
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


def one_body_fsz(ri, spin_i, rj, spin_j):
    return [
        ("c", ri + spin_i * NSITE),
        ("a", rj + spin_j * NSITE),
    ]


def coulomb_intra(site):
    return [
        ("c", site),
        ("a", site),
        ("c", site + NSITE),
        ("a", site + NSITE),
    ]


def general_pairing_values():
    return [
        math.sin(0.37 * float(idx + 1))
        for idx in range((2 * NSITE) * (2 * NSITE - 1) // 2)
    ]


def general_pairing_amplitude(state, values=None):
    occupied_orbitals = [orb for orb in range(2 * NSITE) if occupied(state, orb)]
    if len(occupied_orbitals) != 4:
        return 0.0

    source_values = general_pairing_values() if values is None else values
    pair_values = {}
    idx = 0
    for left in range(2 * NSITE):
        for right in range(left + 1, 2 * NSITE):
            pair_values[(left, right)] = source_values[idx]
            idx += 1

    def pair(left, right):
        if left < right:
            return pair_values[(left, right)]
        return -pair_values[(right, left)]

    a, b, c, d = occupied_orbitals
    return (
        pair(a, b) * pair(c, d)
        - pair(a, c) * pair(b, d)
        + pair(a, d) * pair(b, c)
    )


def doublon_count(state):
    return sum(
        occupied(state, site) * occupied(state, site + NSITE)
        for site in range(NSITE)
    )


def projection_branch_pairing_values():
    values = general_pairing_values()
    pair_index = {}
    idx = 0
    for left in range(2 * NSITE):
        for right in range(left + 1, 2 * NSITE):
            pair_index[(left, right)] = idx
            idx += 1

    base_amplitude = general_pairing_amplitude(
        PROJECTION_BRANCH_BASE_STATE, values)
    desired = PROJECTION_BRANCH_RAW_RATIO * base_amplitude
    p03 = values[pair_index[(0, 3)]]
    p04 = values[pair_index[(0, 4)]]
    p06 = values[pair_index[(0, 6)]]
    p34 = values[pair_index[(3, 4)]]
    p36 = values[pair_index[(3, 6)]]
    values[pair_index[(4, 6)]] = (p04 * p36 - p06 * p34 + desired) / p03
    return values


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


def assert_tight_close(name, actual, expected,
                       abs_tol=1.0e-10, rel_tol=1.0e-9):
    diff = abs(actual - expected)
    tolerance = abs_tol + rel_tol * max(abs(actual), abs(expected))
    if diff > tolerance:
        raise AssertionError(
            "{} mismatch: actual={} expected={} diff={} tolerance={}".format(
                name, actual, expected, diff, tolerance
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


def spin_changing_local_reference():
    basis = [state for state in range(1 << (2 * NSITE)) if popcount(state) == 4]
    psi = np.array([general_pairing_amplitude(state) for state in basis], dtype=float)
    terms = [
        (-0.70 - 0.20j, one_body_fsz(0, 0, 0, 1)),
        (-0.70 + 0.20j, one_body_fsz(0, 1, 0, 0)),
        (0.40 - 0.15j, one_body_fsz(1, 0, 0, 1)),
        (0.40 + 0.15j, one_body_fsz(0, 1, 1, 0)),
    ]
    terms.extend(
        (coupling, coulomb_intra(site))
        for site, coupling in enumerate(COULOMB_INTRA)
    )
    ham = operator_matrix_complex(basis, terms)
    hpsi = np.dot(ham, psi)
    h2psi = np.dot(ham, hpsi)
    result = {}
    for idx, state in enumerate(basis):
        if abs(psi[idx]) <= 1.0e-13:
            continue
        # mVMC stores the bra-oriented estimator <psi|H|x>/<psi|x>.
        # For a Hermitian Hamiltonian this is the complex conjugate of the
        # conventional ket-oriented (H psi)(x)/psi(x) local estimator.
        result[state] = (
            np.conj(hpsi[idx] / psi[idx]),
            np.conj(h2psi[idx] / psi[idx]),
        )
    return result


def write_complex_spin_changing_transfer(workdir):
    rows = [
        (0, 0, 0, 1, 0.70, 0.20),
        (0, 1, 0, 0, 0.70, -0.20),
        (1, 0, 0, 1, -0.40, 0.15),
        (0, 1, 1, 0, -0.40, -0.15),
    ]
    with open(os.path.join(workdir, "trans.def"), "w") as fp:
        fp.write("========================\n")
        fp.write("NTransfer      {}\n".format(len(rows)))
        fp.write("========================\n")
        fp.write("========i_j_s_tijs=======\n")
        fp.write("========================\n")
        for ri, spin_i, rj, spin_j, real, imag in rows:
            fp.write(
                "{:5d} {:5d} {:5d} {:5d} {:25.15f} {:25.15f}\n".format(
                    ri, spin_i, rj, spin_j, real, imag,
                )
            )


def write_coulomb_intra(workdir):
    with open(os.path.join(workdir, "coulombintra.def"), "w") as fp:
        fp.write("=============================================\n")
        fp.write("NCoulombIntra          {}\n".format(len(COULOMB_INTRA)))
        fp.write("=============================================\n")
        fp.write("================== CoulombIntra ================\n")
        fp.write("=============================================\n")
        for site, coupling in enumerate(COULOMB_INTRA):
            fp.write("{:5d} {:25.15f}\n".format(site, coupling))
    append_namelist_entry(workdir, "CoulombIntra", "coulombintra.def")


def write_gutzwiller(workdir):
    with open(os.path.join(workdir, "gutzwilleridx.def"), "w") as fp:
        fp.write("=============================================\n")
        fp.write("NGutzwillerIdx          1\n")
        fp.write("ComplexType             0\n")
        fp.write("=============================================\n")
        fp.write("=============================================\n")
        for site in range(NSITE):
            fp.write("{:5d} {:5d}\n".format(site, 0))
        fp.write("{:5d} {:5d}\n".format(0, 0))
    append_namelist_entry(workdir, "Gutzwiller", "gutzwilleridx.def")


def write_projection_branch_init(workdir, pairing_values):
    values = [0.0] * 6
    values.extend([PROJECTION_BRANCH_GUTZWILLER, 0.0, 0.0])
    for value in pairing_values:
        values.extend([value, 0.0, 0.0])
    filename = "projection_branch_init.dat"
    with open(os.path.join(workdir, filename), "w") as fp:
        fp.write(" ".join("{:.18e}".format(value) for value in values))
        fp.write("\n")
    return filename


def read_lanczos_oracle_dump(path):
    rows = []
    with open(path) as fp:
        for line in fp:
            cols = line.split()
            if len(cols) != 3 + 2 * NSITE + 1 + 8:
                raise AssertionError("invalid FSZ Lanczos oracle row: {}".format(line.rstrip()))
            if cols[0] != "sample" or cols[2] != "occ" or cols[3 + 2 * NSITE] != "lslq":
                raise AssertionError("invalid FSZ Lanczos oracle markers")
            occupation = [int(value) for value in cols[3:3 + 2 * NSITE]]
            state = 0
            for orb, value in enumerate(occupation):
                if value:
                    state |= 1 << orb
            values = [float(value) for value in cols[4 + 2 * NSITE:]]
            lslq = [complex(values[2 * idx], values[2 * idx + 1]) for idx in range(4)]
            rows.append((int(cols[1]), state, lslq))
    if not rows:
        raise AssertionError("FSZ Lanczos oracle dump is empty")
    return rows


def check_spin_changing_fsz_local_oracle(rootdir):
    case_name = "LanczosSectorED_FSZSpinChanging"
    workdir = prepare_case(
        rootdir,
        case_name,
        include_backflow=False,
        orbital_complex_type=1,
        orbital_optimize=True,
    )
    update_modpara(workdir, {
        "NLanczosMode": "1",
        "NVMCCalMode": "1",
        "NVMCSample": "128",
        "NVMCWarmUp": "16",
        "2Sz": "-1",
    })
    write_spin_changing_c1_defs(workdir)
    write_complex_spin_changing_transfer(workdir)
    write_coulomb_intra(workdir)
    init_path = write_lanczos_init(workdir)
    dump_name = "lanczos_sector_ed_oracle.dat"
    proc = run_vmc(
        rootdir,
        workdir,
        init_path=init_path,
        extra_env={
            "MVMC_BF_LANCZOS_STATE_CHECK": "1",
            "MVMC_LANCZOS_ORACLE_DUMP": dump_name,
        },
    )
    if proc.returncode != 0:
        raise AssertionError("FSZ spin-changing vmc.out failed:\n{}".format(proc.stdout))

    reference = spin_changing_local_reference()
    rows = read_lanczos_oracle_dump(os.path.join(workdir, dump_name))
    sampled_nup = set()
    for sample, state, lslq in rows:
        if state not in reference:
            raise AssertionError("sample {} has zero or missing ED amplitude".format(sample))
        expected_h1, expected_h2 = reference[state]
        assert_tight_close(
            "sample {} h1.real".format(sample),
            lslq[1].real, expected_h1.real)
        assert_tight_close(
            "sample {} h1.imag".format(sample),
            lslq[1].imag, expected_h1.imag)
        assert_tight_close(
            "sample {} h2.real".format(sample),
            lslq[3].real, expected_h2.real)
        assert_tight_close(
            "sample {} h2.imag".format(sample),
            lslq[3].imag, expected_h2.imag)
        sampled_nup.add(sum(occupied(state, site) for site in range(NSITE)))
    if len(sampled_nup) < 2:
        raise AssertionError("spin-changing oracle did not sample multiple spin sectors")


def projection_branch_local_reference(pairing_values):
    basis = [state for state in range(1 << (2 * NSITE)) if popcount(state) == 4]
    psi = np.array([
        general_pairing_amplitude(state, pairing_values) * math.exp(
            PROJECTION_BRANCH_GUTZWILLER * doublon_count(state))
        for state in basis
    ], dtype=float)
    terms = [
        (0.40, one_body_fsz(0, 1, 1, 0)),
        (0.40, one_body_fsz(1, 0, 0, 1)),
    ]
    ham = operator_matrix_complex(basis, terms)
    hpsi = np.dot(ham, psi)
    h2psi = np.dot(ham, hpsi)
    result = {}
    for idx, state in enumerate(basis):
        if abs(psi[idx]) <= 1.0e-14:
            continue
        result[state] = (
            np.conj(hpsi[idx] / psi[idx]),
            np.conj(h2psi[idx] / psi[idx]),
        )
    return result


def write_projection_branch_transfer(workdir):
    rows = [
        (0, 1, 1, 0, -0.40, 0.0),
        (1, 0, 0, 1, -0.40, 0.0),
    ]
    with open(os.path.join(workdir, "trans.def"), "w") as fp:
        fp.write("========================\n")
        fp.write("NTransfer      {}\n".format(len(rows)))
        fp.write("========================\n")
        fp.write("========i_j_s_tijs=======\n")
        fp.write("========================\n")
        for ri, spin_i, rj, spin_j, real, imag in rows:
            fp.write(
                "{:5d} {:5d} {:5d} {:5d} {:25.15f} {:25.15f}\n".format(
                    ri, spin_i, rj, spin_j, real, imag,
                )
            )


def check_projection_branch_fsz_local_oracle(rootdir):
    pairing_values = projection_branch_pairing_values()
    base_amplitude = general_pairing_amplitude(
        PROJECTION_BRANCH_BASE_STATE, pairing_values)
    target_amplitude = general_pairing_amplitude(
        PROJECTION_BRANCH_TARGET_STATE, pairing_values)
    raw_ratio = abs(target_amplitude / base_amplitude)
    projection_delta = (
        doublon_count(PROJECTION_BRANCH_TARGET_STATE)
        - doublon_count(PROJECTION_BRANCH_BASE_STATE)
    )
    projected_ratio = raw_ratio * math.exp(
        PROJECTION_BRANCH_GUTZWILLER * projection_delta)
    if not raw_ratio < 1.0e-12:
        raise AssertionError("projection branch fixture raw ratio is not below threshold")
    if not projected_ratio > 1.0e-12:
        raise AssertionError("projection branch fixture full ratio is not above threshold")

    case_name = "LanczosSectorED_FSZProjectionBranch"
    workdir = prepare_case(
        rootdir,
        case_name,
        include_backflow=False,
        orbital_complex_type=1,
        orbital_optimize=True,
    )
    update_modpara(workdir, {
        "NLanczosMode": "1",
        "NVMCCalMode": "1",
        "NVMCSample": "4096",
        "NVMCWarmUp": "32",
        "2Sz": "-1",
    })
    write_projection_branch_transfer(workdir)
    write_gutzwiller(workdir)
    init_path = write_projection_branch_init(workdir, pairing_values)
    dump_name = "lanczos_projection_branch_oracle.dat"
    proc = run_vmc(
        rootdir,
        workdir,
        init_path=init_path,
        extra_env={
            "MVMC_BF_LANCZOS_STATE_CHECK": "1",
            "MVMC_LANCZOS_ORACLE_DUMP": dump_name,
            "MVMC_LANCZOS_TEST_PROJECTION_BRANCH_AUDIT": "1",
        },
    )
    if proc.returncode != 0:
        raise AssertionError("FSZ projection-branch vmc.out failed:\n{}".format(proc.stdout))

    reference = projection_branch_local_reference(pairing_values)
    rows = read_lanczos_oracle_dump(os.path.join(workdir, dump_name))
    sampled_states = set()
    for sample, state, lslq in rows:
        if state not in reference:
            raise AssertionError("sample {} has zero or missing ED amplitude".format(sample))
        expected_h1, expected_h2 = reference[state]
        assert_tight_close(
            "projection sample {} h1".format(sample), lslq[1], expected_h1)
        assert_tight_close(
            "projection sample {} h2".format(sample), lslq[3], expected_h2)
        sampled_states.add(state)
    if PROJECTION_BRANCH_BASE_STATE not in sampled_states:
        raise AssertionError(
            "projection branch fixture did not sample the threshold-crossing base state")


def main():
    rootdir = os.getcwd()
    check_projected_operators()
    check_doublon_only_against_mvmc(rootdir)
    check_spin_changing_fsz_local_oracle(rootdir)
    check_projection_branch_fsz_local_oracle(rootdir)
    print("Lanczos sector ED verification passed")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AssertionError as exc:
        print("ERROR: {}".format(exc))
        sys.exit(-1)
