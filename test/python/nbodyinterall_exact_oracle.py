from __future__ import print_function

import math
import os
import shutil
import subprocess
import sys

import numpy as np

from nbodyg_exact_oracle import (
    apply_ops,
    assert_config_close,
    config_oracle_dump_paths,
    enable_nonidentity_backflow,
    enable_config_oracle_backflow,
    fixed_sz_basis,
    fmt,
    nbody_ops,
    occupied,
    pfaffian,
    parse_config_oracle_dump,
    mpi_command,
    update_modpara,
    validated_work_suffix,
    write,
    write_antiparallel_orbital,
    write_common_defs,
    write_general_orbital,
    write_nbodyg_def,
    validate_config_oracle_rows,
)

NBODYINTERALL_ENERGY_TOL = 1.2e-2


def write_common_nbodyinterall_defs(workdir, nsite, sample_count, seed,
                                    orbital_keyword, orbital_file,
                                    data_qty=1, n_mp_trans=1):
    write_common_defs(
        workdir,
        nsite,
        sample_count,
        seed,
        orbital_keyword,
        orbital_file,
        data_qty=data_qty,
        n_mp_trans=n_mp_trans,
    )
    write(
        os.path.join(workdir, "namelist.def"),
        (
            "         ModPara  modpara.def\n"
            "         LocSpin  locspn.def\n"
            "  NBodyInterAll  nbodyinterall.def\n"
            "  {0:>14s}  {1}\n"
            "        TransSym  qptransidx.def\n"
        ).format(orbital_keyword, orbital_file),
    )


def write_nbodyinterall_def(workdir, terms):
    body = []
    for nbody, factors, coef in terms:
        cols = [str(nbody)] + [str(x) for x in factors]
        cols += [fmt(coef.real), fmt(coef.imag)]
        body.append("{}\n".format(" ".join(cols)))
    write(
        os.path.join(workdir, "nbodyinterall.def"),
        (
            "=============================================\n"
            "NNBodyInterAll   {0}\n"
            "=============================================\n"
            "======== NBodyInterAll interactions =========\n"
            "=============================================\n"
            "{1}"
        ).format(len(terms), "".join(body)),
    )


def exact_nbodyinterall_energy(nsite, nup, ndown, slater_elm, terms):
    basis = fixed_sz_basis(nsite, nup, ndown)
    psi = {}
    for state in basis:
        occ = [orb for orb in range(2 * nsite) if occupied(state, orb)]
        sub = [[slater_elm[i, j] for j in occ] for i in occ]
        psi[state] = pfaffian(sub)

    norm = sum(abs(amp) ** 2 for amp in psi.values())
    if norm == 0.0:
        raise AssertionError("oracle wave function has zero norm")

    energy = 0.0 + 0.0j
    for nbody, factors, coef in terms:
        if nbody < 1:
            raise AssertionError("NBodyInterAll oracle received invalid order")
        ops = nbody_ops(nsite, factors)
        numerator = 0.0 + 0.0j
        for state, amp in psi.items():
            result = apply_ops(state, ops)
            if result is None:
                continue
            next_state, sign = result
            if next_state in psi:
                numerator += np.conj(psi[next_state]) * sign * amp
        energy += coef * numerator / norm
    return energy


def complex_antiparallel_state(nsite):
    f_matrix = np.zeros((nsite, nsite), dtype=complex)
    for i in range(nsite):
        for j in range(nsite):
            real = math.sin(1.3 * (i + 1) * (j + 2))
            real += 0.35 * math.cos(0.4 * (i + 2) * (j + 1))
            imag = math.cos(0.9 * (i + 1) * (j + 1))
            imag -= 0.25 * math.sin((i + 3) * (j + 1))
            f_matrix[i, j] = real + 1j * imag

    slater_elm = np.zeros((2 * nsite, 2 * nsite), dtype=complex)
    for i in range(nsite):
        for j in range(nsite):
            slater_elm[i, nsite + j] = f_matrix[i, j]
            slater_elm[nsite + j, i] = -f_matrix[i, j]
    return f_matrix, slater_elm


def complex_general_state(nsite):
    slater_elm = np.zeros((2 * nsite, 2 * nsite), dtype=complex)
    for i in range(2 * nsite):
        for j in range(i + 1, 2 * nsite):
            real = math.sin(0.7 * (i + 1) * (j + 2))
            real += 0.5 * math.cos(0.3 * (i + 2) * (j + 1))
            imag = math.cos(0.5 * (i + 3) * (j + 1))
            imag += 0.2 * math.sin((i + 1) * (j + 4))
            slater_elm[i, j] = real + 1j * imag
            slater_elm[j, i] = -slater_elm[i, j]
    return slater_elm


def n1_transfer_like_case(workdir):
    nsite = 6
    terms = [
        (1, (0, 0, 2, 0), 0.37 + 0.21j),
        (1, (4, 1, 1, 1), -0.19 + 0.11j),
    ]
    f_matrix, slater_elm = complex_antiparallel_state(nsite)
    write_common_nbodyinterall_defs(workdir, nsite, 12000, 11223,
                                    "Orbital", "orbitalidx.def")
    write_nbodyinterall_def(workdir, terms)
    write_antiparallel_orbital(workdir, f_matrix)
    return nsite, 3, 3, slater_elm, terms


def n2_interall_like_case(workdir):
    nsite = 6
    terms = [
        (2, (0, 0, 2, 0, 5, 1, 1, 1), -0.28 + 0.16j),
        (2, (3, 0, 1, 0, 4, 1, 0, 1), 0.41 - 0.09j),
    ]
    f_matrix, slater_elm = complex_antiparallel_state(nsite)
    write_common_nbodyinterall_defs(workdir, nsite, 12000, 22334,
                                    "Orbital", "orbitalidx.def")
    write_nbodyinterall_def(workdir, terms)
    write_antiparallel_orbital(workdir, f_matrix)
    return nsite, 3, 3, slater_elm, terms


def n3_oracle_case(workdir):
    nsite = 6
    terms = [
        (3, (0, 0, 2, 0, 1, 0, 3, 0, 5, 1, 1, 1), 0.22 - 0.17j),
        (3, (0, 0, 3, 0, 1, 0, 2, 0, 5, 1, 1, 1), -0.31 + 0.13j),
        (3, (2, 0, 0, 0, 3, 0, 1, 0, 1, 1, 5, 1), 0.18 + 0.27j),
    ]
    f_matrix, slater_elm = complex_antiparallel_state(nsite)
    write_common_nbodyinterall_defs(workdir, nsite, 12000, 33445,
                                    "Orbital", "orbitalidx.def")
    write_nbodyinterall_def(workdir, terms)
    write_antiparallel_orbital(workdir, f_matrix)
    return nsite, 3, 3, slater_elm, terms


def fsz_n3_oracle_case(workdir):
    nsite = 6
    terms = [
        (3, (0, 0, 2, 1, 3, 1, 1, 0, 1, 1, 4, 1), 0.23 + 0.19j),
        (3, (0, 0, 4, 1, 3, 1, 1, 0, 1, 1, 2, 1), -0.17 + 0.29j),
        (3, (1, 0, 3, 1, 2, 1, 0, 0, 4, 1, 1, 1), 0.37 - 0.08j),
    ]
    slater_elm = complex_general_state(nsite)
    write_common_nbodyinterall_defs(workdir, nsite, 12000, 44556,
                                    "OrbitalGeneral", "orbitalidxgen.def")
    write_nbodyinterall_def(workdir, terms)
    write_general_orbital(workdir, slater_elm)
    return nsite, 3, 3, slater_elm, terms


def zero_terms_case(workdir):
    nsite = 6
    terms = []
    f_matrix, slater_elm = complex_antiparallel_state(nsite)
    write_common_nbodyinterall_defs(workdir, nsite, 200, 55667,
                                    "Orbital", "orbitalidx.def")
    write_nbodyinterall_def(workdir, terms)
    write_antiparallel_orbital(workdir, f_matrix)
    return nsite, 3, 3, slater_elm, terms


def projector_multiterm_smoke_case(workdir):
    nsite = 6
    terms = [
        (3, (0, 0, 2, 0, 1, 0, 3, 0, 5, 1, 1, 1), 0.22 - 0.17j),
        (3, (0, 0, 3, 0, 1, 0, 2, 0, 5, 1, 1, 1), -0.31 + 0.13j),
        (3, (2, 0, 0, 0, 3, 0, 1, 0, 1, 1, 5, 1), 0.18 + 0.27j),
        (3, (0, 0, 2, 0, 1, 0, 5, 0, 5, 1, 1, 1), -0.11 + 0.09j),
    ]
    f_matrix, slater_elm = complex_antiparallel_state(nsite)
    write_common_nbodyinterall_defs(
        workdir,
        nsite,
        400,
        66778,
        "Orbital",
        "orbitalidx.def",
        data_qty=1,
        n_mp_trans=2,
    )
    write_nbodyinterall_def(workdir, terms)
    write_antiparallel_orbital(workdir, f_matrix)
    return None


def backflow_n1_energy_case(workdir):
    nsite = 6
    terms = [
        (1, (0, 0, 2, 0), 0.37 + 0.21j),
        (1, (4, 1, 1, 1), -0.19 + 0.11j),
    ]
    f_matrix, unused_slater = complex_antiparallel_state(nsite)
    write_common_nbodyinterall_defs(
        workdir, nsite, 400, 71501, "Orbital", "orbitalidx.def"
    )
    write_nbodyinterall_def(workdir, terms)
    slater_values = write_antiparallel_orbital(workdir, f_matrix)
    enable_nonidentity_backflow(workdir, nsite, slater_values)
    return None


def backflow_n2_energy_case(workdir):
    nsite = 6
    terms = [
        (2, (0, 0, 2, 0, 5, 1, 1, 1), -0.28 + 0.16j),
        (2, (3, 0, 1, 0, 4, 1, 0, 1), 0.41 - 0.09j),
    ]
    f_matrix, unused_slater = complex_antiparallel_state(nsite)
    write_common_nbodyinterall_defs(
        workdir, nsite, 400, 71502, "Orbital", "orbitalidx.def"
    )
    write_nbodyinterall_def(workdir, terms)
    slater_values = write_antiparallel_orbital(workdir, f_matrix)
    enable_nonidentity_backflow(workdir, nsite, slater_values)
    return None


def backflow_mixed_energy_case(workdir):
    nsite = 6
    terms = [
        # Effective N=1 dispatch and effective N=2 checked rebuild.
        (1, (0, 0, 2, 0), 0.19 + 0.07j),
        (2, (0, 0, 2, 0, 5, 1, 1, 1), -0.17 + 0.13j),
        # Genuine N=3 and N=4 full rebuilds.
        (3, (0, 0, 2, 0, 1, 0, 3, 0, 5, 1, 1, 1), 0.23 - 0.11j),
        (
            4,
            (
                0, 0, 2, 0,
                1, 0, 3, 0,
                5, 1, 1, 1,
                4, 1, 0, 1,
            ),
            -0.09 + 0.21j,
        ),
        # N=3 contractions to lower effective order.
        (3, (5, 0, 2, 0, 2, 0, 1, 0, 4, 1, 0, 1), 0.08 + 0.06j),
        (3, (2, 0, 1, 0, 4, 0, 2, 0, 5, 0, 5, 0), -0.07 + 0.04j),
        # Complete scalar contraction and exact repeated-annihilation zero.
        (3, (0, 0, 0, 0, 1, 0, 1, 0, 2, 1, 2, 1), 0.03 - 0.02j),
        (3, (3, 0, 0, 0, 4, 0, 0, 0, 5, 1, 1, 1), 0.31 + 0.17j),
    ]
    f_matrix, unused_slater = complex_antiparallel_state(nsite)
    write_common_nbodyinterall_defs(
        workdir, nsite, 600, 71503, "Orbital", "orbitalidx.def"
    )
    write_nbodyinterall_def(workdir, terms)
    slater_values = write_antiparallel_orbital(workdir, f_matrix)
    enable_nonidentity_backflow(workdir, nsite, slater_values)
    return None


def backflow_zero_terms_case(workdir):
    nsite = 6
    terms = []
    f_matrix, unused_slater = complex_antiparallel_state(nsite)
    write_common_nbodyinterall_defs(
        workdir, nsite, 200, 71504, "Orbital", "orbitalidx.def"
    )
    write_nbodyinterall_def(workdir, terms)
    slater_values = write_antiparallel_orbital(workdir, f_matrix)
    enable_nonidentity_backflow(workdir, nsite, slater_values)
    return None


def backflow_fsz_spin_changing_mixed_energy_case(workdir):
    nsite = 6
    terms = [
        (1, (0, 0, 2, 1), 0.19 + 0.07j),
        (2, (0, 0, 2, 1, 3, 1, 1, 1), -0.17 + 0.13j),
        (
            3,
            (0, 0, 2, 1, 3, 1, 1, 0, 1, 1, 4, 1),
            0.23 - 0.11j,
        ),
        (
            4,
            (
                0, 0, 2, 1,
                1, 0, 3, 0,
                5, 1, 1, 1,
                4, 1, 0, 0,
            ),
            -0.09 + 0.21j,
        ),
    ]
    slater_elm = complex_general_state(nsite)
    write_common_nbodyinterall_defs(
        workdir, nsite, 600, 71505, "OrbitalGeneral", "orbitalidxgen.def"
    )
    write_nbodyinterall_def(workdir, terms)
    slater_values = write_general_orbital(workdir, slater_elm)
    enable_nonidentity_backflow(workdir, nsite, slater_values)
    update_modpara(workdir, "2Sz", -1)
    return None


def backflow_fsz_zero_terms_case(workdir):
    nsite = 6
    terms = []
    slater_elm = complex_general_state(nsite)
    write_common_nbodyinterall_defs(
        workdir, nsite, 200, 71506, "OrbitalGeneral", "orbitalidxgen.def"
    )
    write_nbodyinterall_def(workdir, terms)
    write_nbodyg_def(workdir, [])
    with open(os.path.join(workdir, "namelist.def"), "a") as f:
        f.write("          NBodyG  nbodyg.def\n")
    slater_values = write_general_orbital(workdir, slater_elm)
    enable_nonidentity_backflow(workdir, nsite, slater_values)
    update_modpara(workdir, "2Sz", -1)
    return None


def backflow_mixed_config_oracle_case(workdir):
    nsite = 6
    terms = [
        (1, (0, 0, 2, 0), 0.19 + 0.07j),
        (2, (0, 0, 2, 0, 5, 1, 1, 1), -0.17 + 0.13j),
        (3, (0, 0, 2, 0, 1, 0, 3, 0, 5, 1, 1, 1), 0.23 - 0.11j),
        (
            4,
            (
                0, 0, 2, 0,
                1, 0, 3, 0,
                5, 1, 1, 1,
                4, 1, 0, 1,
            ),
            -0.09 + 0.21j,
        ),
        (3, (5, 0, 2, 0, 2, 0, 1, 0, 4, 1, 0, 1), 0.08 + 0.06j),
        (3, (2, 0, 1, 0, 4, 0, 2, 0, 5, 0, 5, 0), -0.07 + 0.04j),
        (3, (0, 0, 0, 0, 1, 0, 1, 0, 2, 1, 2, 1), 0.03 - 0.02j),
        (3, (3, 0, 0, 0, 4, 0, 0, 0, 5, 1, 1, 1), 0.31 + 0.17j),
    ]
    categories = [
        "dispatch1",
        "effective2",
        "full3",
        "full4",
        "contraction",
        "contraction",
        "scalar",
        "zero",
    ]
    f_matrix, unused_slater = complex_antiparallel_state(nsite)
    write_common_nbodyinterall_defs(
        workdir, nsite, 64, 90433, "Orbital", "orbitalidx.def"
    )
    write_nbodyinterall_def(workdir, terms)
    slater_values = write_antiparallel_orbital(workdir, f_matrix)
    enable_config_oracle_backflow(workdir, nsite, slater_values)
    return nsite, terms, categories, {
        "dispatch1",
        "effective2",
        "full3",
        "full4",
        "contraction",
        "scalar",
        "zero",
    }


def backflow_fsz_mixed_config_oracle_case(workdir):
    nsite = 6
    terms = [
        (1, (0, 0, 2, 1), 0.19 + 0.07j),
        (2, (0, 0, 2, 1, 3, 1, 1, 1), -0.17 + 0.13j),
        (
            3,
            (0, 0, 2, 1, 3, 1, 1, 0, 1, 1, 4, 1),
            0.23 - 0.11j,
        ),
        (
            4,
            (
                0, 0, 2, 1,
                1, 0, 3, 0,
                5, 1, 1, 1,
                4, 1, 0, 1,
            ),
            -0.09 + 0.21j,
        ),
    ]
    categories = ["dispatch1", "effective2", "full3", "full4"]
    slater_elm = complex_general_state(nsite)
    write_common_nbodyinterall_defs(
        workdir, nsite, 1024, 90437, "OrbitalGeneral", "orbitalidxgen.def"
    )
    write_nbodyinterall_def(workdir, terms)
    slater_values = write_general_orbital(workdir, slater_elm)
    enable_config_oracle_backflow(workdir, nsite, slater_values)
    update_modpara(workdir, "2Sz", -1)
    return nsite, terms, categories, {
        "dispatch1", "effective2", "full3", "full4"
    }


EXACT_CASES = {
    "NBodyInterAll_N1_TransferLike": n1_transfer_like_case,
    "NBodyInterAll_N2_InterAllLike": n2_interall_like_case,
    "NBodyInterAll_N3_Oracle": n3_oracle_case,
    "NBodyInterAll_fsz_N3_Oracle": fsz_n3_oracle_case,
    "NBodyInterAll_ZeroTerms": zero_terms_case,
}

SMOKE_CASES = {
    "NBodyInterAll_Projector_N3_MultiTerm": projector_multiterm_smoke_case,
    "BackFlow_NBodyInterAll_NonIdentity_N1_Energy": backflow_n1_energy_case,
    "BackFlow_NBodyInterAll_NonIdentity_N2_Energy": backflow_n2_energy_case,
    "BackFlow_NBodyInterAll_NonIdentity_Mixed_Energy": backflow_mixed_energy_case,
    "BackFlow_NBody_ZeroTerms_NonFSZ": backflow_zero_terms_case,
    "BackFlow_FSZ_NBodyInterAll_SpinChanging_Mixed_Energy":
        backflow_fsz_spin_changing_mixed_energy_case,
    "BackFlow_FSZ_NBody_ZeroTerms": backflow_fsz_zero_terms_case,
}

CONFIG_ORACLE_CASES = {
    "BackFlow_NBodyInterAll_NonIdentity_Mixed_ConfigOracle":
        backflow_mixed_config_oracle_case,
    "BackFlow_FSZ_NBodyInterAll_Mixed_ConfigOracle":
        backflow_fsz_mixed_config_oracle_case,
}

NONVACUOUS_ENERGY_CASES = {
    "BackFlow_NBodyInterAll_NonIdentity_N1_Energy",
    "BackFlow_NBodyInterAll_NonIdentity_N2_Energy",
    "BackFlow_NBodyInterAll_NonIdentity_Mixed_Energy",
    "BackFlow_FSZ_NBodyInterAll_SpinChanging_Mixed_Energy",
}


def parse_energy(workdir):
    path = os.path.join(workdir, "output", "zvo_out_001.dat")
    with open(path) as f:
        cols = f.readline().split()
    if len(cols) < 2:
        raise AssertionError("malformed energy output: {}".format(path))
    energy = complex(float(cols[0]), float(cols[1]))
    if not (math.isfinite(energy.real) and math.isfinite(energy.imag)):
        raise AssertionError("energy is not finite: {}".format(energy))
    return energy


def assert_close(
    label, actual, expected, tol=NBODYINTERALL_ENERGY_TOL
):
    diff = abs(actual - expected)
    if diff > tol:
        raise AssertionError(
            "{} mismatch: actual={} expected={} diff={} tol={}".format(
                label, actual, expected, diff, tol
            )
        )


def main():
    all_cases = sorted(
        list(EXACT_CASES) + list(SMOKE_CASES) + list(CONFIG_ORACLE_CASES)
    )
    if len(sys.argv) != 2 or sys.argv[1] not in all_cases:
        print("usage: {} {}".format(sys.argv[0], "|".join(all_cases)))
        return -1

    model = sys.argv[1]
    rootdir = os.getcwd()
    work_suffix = validated_work_suffix()
    workdir = os.path.join(rootdir, "work", model + work_suffix)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)

    if model in EXACT_CASES:
        nsite, nup, ndown, slater_elm, terms = EXACT_CASES[model](workdir)
        expected = exact_nbodyinterall_energy(nsite, nup, ndown, slater_elm, terms)
    elif model in SMOKE_CASES:
        SMOKE_CASES[model](workdir)
        expected = None
    else:
        nsite, terms, categories, required_categories = (
            CONFIG_ORACLE_CASES[model](workdir)
        )
        expected = None

    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    if mpi_procs:
        cmd = mpi_command(
            mpi_procs,
            [bin_to_test, "-e", "namelist.def", "initial.def"])
    else:
        cmd = [bin_to_test, "-e", "namelist.def", "initial.def"]
    result = subprocess.call(cmd, cwd=workdir)
    if result != 0:
        return result

    actual = parse_energy(workdir)
    if model in CONFIG_ORACLE_CASES:
        rows = []
        for rank, dump_path in enumerate(config_oracle_dump_paths(workdir)):
            rows.extend(parse_config_oracle_dump(dump_path, nsite, rank))
        production_by_sample, unused_expected, unused_values = (
            validate_config_oracle_rows(
                rows,
                nsite,
                "h",
                terms,
                categories,
                required_categories,
            )
        )
        average = sum(
            production_by_sample.values(), 0.0 + 0.0j
        ) / float(len(production_by_sample))
        assert_config_close(
            "configuration oracle final energy", actual, average
        )
        print("{} configuration oracle passed".format(model))
        return 0

    if expected is not None:
        assert_close("{} energy".format(model), actual, expected)
    if model in NONVACUOUS_ENERGY_CASES and abs(actual) <= 1.0e-12:
        raise AssertionError("{} energy is vacuous".format(model))
    print("{} NBodyInterAll energy check passed".format(model))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print("ERROR: {}".format(exc))
        sys.exit(-1)
