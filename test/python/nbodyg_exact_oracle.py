from __future__ import print_function

import itertools
import math
import os
import shutil
import subprocess
import sys

import numpy as np

from backflow_def_helper import (
    build_chain_nn_backflow,
    write_chain_nn_backflow,
)


TOL = 1.2e-2


def pfaffian(mat):
    n = len(mat)
    if n == 0:
        return 1.0 + 0.0j
    if n == 2:
        return mat[0][1]

    total = 0.0 + 0.0j
    for j in range(1, n):
        sub = [
            [mat[x][y] for y in range(n) if y not in (0, j)]
            for x in range(n)
            if x not in (0, j)
        ]
        total += ((-1) ** (j + 1)) * mat[0][j] * pfaffian(sub)
    return total


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
    coeff = 1
    current = state
    for kind, orb in reversed(ops):
        result = create(current, orb) if kind == "c" else annihilate(current, orb)
        if result is None:
            return None
        current, sign = result
        coeff *= sign
    return current, coeff


def state_from_sites(nsite, ups, downs):
    state = 0
    for site in ups:
        state |= 1 << site
    for site in downs:
        state |= 1 << (site + nsite)
    return state


def fixed_sz_basis(nsite, nup, ndown):
    basis = []
    for ups in itertools.combinations(range(nsite), nup):
        for downs in itertools.combinations(range(nsite), ndown):
            basis.append(state_from_sites(nsite, ups, downs))
    return basis


def nbody_ops(nsite, factors):
    ops = []
    for k in range(len(factors) // 4):
        site_out, spin_out, site_in, spin_in = factors[4 * k:4 * k + 4]
        ops.append(("c", site_out + spin_out * nsite))
        ops.append(("a", site_in + spin_in * nsite))
    return ops


def exact_nbody_values(nsite, nup, ndown, slater_elm, terms):
    basis = fixed_sz_basis(nsite, nup, ndown)
    psi = {}
    for state in basis:
        occ = [orb for orb in range(2 * nsite) if occupied(state, orb)]
        sub = [[slater_elm[i, j] for j in occ] for i in occ]
        psi[state] = pfaffian(sub)

    norm = sum(abs(amp) ** 2 for amp in psi.values())
    if norm == 0.0:
        raise AssertionError("oracle wave function has zero norm")

    values = []
    for nbody, factors in terms:
        if nbody != 3:
            raise AssertionError("oracle cases are intended to exercise genuine N=3 terms")
        ops = nbody_ops(nsite, factors)
        numerator = 0.0 + 0.0j
        for state, amp in psi.items():
            result = apply_ops(state, ops)
            if result is None:
                continue
            next_state, sign = result
            if next_state in psi:
                numerator += np.conj(psi[next_state]) * sign * amp
        values.append(numerator / norm)
    return values


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def fmt(value):
    return "{:.18e}".format(value)


def write_common_defs(workdir, nsite, sample_count, seed, orbital_keyword, orbital_file,
                      data_qty=1, n_mp_trans=1):
    write(
        os.path.join(workdir, "namelist.def"),
        (
            "         ModPara  modpara.def\n"
            "         LocSpin  locspn.def\n"
            "          NBodyG  nbodyg.def\n"
            "  {0:>14s}  {1}\n"
            "        TransSym  qptransidx.def\n"
        ).format(orbital_keyword, orbital_file),
    )
    write(
        os.path.join(workdir, "modpara.def"),
        (
            "--------------------\n"
            "Model_Parameters   0\n"
            "--------------------\n"
            "VMC_Cal_Parameters\n"
            "--------------------\n"
            "CDataFileHead  zvo\n"
            "CParaFileHead  zqp\n"
            "--------------------\n"
            "NVMCCalMode    1\n"
            "NLanczosMode   0\n"
            "--------------------\n"
            "NDataIdxStart  1\n"
            "NDataQtySmp    {data_qty}\n"
            "--------------------\n"
            "Nsite          {nsite}\n"
            "Ncond          6\n"
            "2Sz            0\n"
            "NSPGaussLeg    1\n"
            "NSPStot        0\n"
            "NMPTrans       {n_mp_trans}\n"
            "NSROptItrStep  20\n"
            "NSROptItrSmp   10\n"
            "DSROptRedCut   0.0000000100\n"
            "DSROptStaDel   0.0100000000\n"
            "DSROptStepDt   0.0100000000\n"
            "NVMCWarmUp     200\n"
            "NVMCInterval   1\n"
            "NVMCSample     {sample_count}\n"
            "NExUpdatePath  0\n"
            "RndSeed        {seed}\n"
            "NSplitSize     1\n"
            "NStore         1\n"
            "NSRCG          0\n"
        ).format(
            nsite=nsite,
            sample_count=sample_count,
            seed=seed,
            data_qty=data_qty,
            n_mp_trans=n_mp_trans,
        ),
    )
    write(
        os.path.join(workdir, "locspn.def"),
        (
            "================================ \n"
            "NlocalSpin     0  \n"
            "================================ \n"
            "========i_1LocSpn_0IteElc ====== \n"
            "================================ \n"
            + "".join("{:5d}      0\n".format(i) for i in range(nsite))
        ),
    )
    write(
        os.path.join(workdir, "qptransidx.def"),
        (
            "=============================================\n"
            "NQPTrans          {0}\n"
            "=============================================\n"
            "======== TrIdx_TrWeight_and_TrIdx_i_xi ======\n"
            "=============================================\n"
            "{1}{2}"
        ).format(
            n_mp_trans,
            "".join("{:d}    1.00000    0.00000\n".format(idx) for idx in range(n_mp_trans)),
            "".join(
                "    {tr}      {site}      {target}      1\n".format(
                    tr=tr, site=site, target=(site + tr) % nsite
                )
                for tr in range(n_mp_trans)
                for site in range(nsite)
            ),
        ),
    )


def write_nbodyg_def(workdir, terms):
    body = []
    for nbody, factors in terms:
        body.append("{} {}\n".format(nbody, " ".join(str(x) for x in factors)))
    write(
        os.path.join(workdir, "nbodyg.def"),
        (
            "=============================================\n"
            "NNBodyG          {0}\n"
            "=============================================\n"
            "======== NBodyG correlation functions =======\n"
            "=============================================\n"
            "{1}"
        ).format(len(terms), "".join(body)),
    )


def write_initial(workdir, slater_values):
    values = [0.0] * 6
    for value in slater_values:
        values.extend([float(np.real(value)), float(np.imag(value)), 0.0])
    write(os.path.join(workdir, "initial.def"), " ".join(fmt(x) for x in values) + "\n")


def write_backflow_initial(workdir, nprojbf, slater_values):
    projbf_template = [
        0.93, 0.08, -0.05, 0.035, -0.025,
        0.015, -0.012, 0.010, -0.007, 0.005,
    ]
    if nprojbf < 2:
        raise AssertionError("non-identity BackFlow case needs at least two ProjBF values")

    projbf_values = []
    for idx in range(nprojbf):
        if idx < len(projbf_template):
            projbf_values.append(projbf_template[idx])
        else:
            projbf_values.append(0.001 / float(idx + 1))
    if projbf_values[0] != 0.93 or not any(value != 0.0 for value in projbf_values[1:]):
        raise AssertionError("BackFlow initial parameters are not non-identity")

    values = [0.0] * 6
    for value in projbf_values:
        values.extend([value, 0.0, 0.0])
    for value in slater_values:
        values.extend([float(np.real(value)), float(np.imag(value)), 0.0])
    write(os.path.join(workdir, "initial.def"), " ".join(fmt(x) for x in values) + "\n")
    return projbf_values


def enable_nonidentity_backflow(workdir, nsite, slater_values):
    definition = build_chain_nn_backflow(length=nsite, optimize=False)
    write_chain_nn_backflow(workdir, length=nsite, optimize=False)
    namelist_path = os.path.join(workdir, "namelist.def")
    with open(namelist_path, "a") as f:
        f.write("              BF  bf.def\n")
        f.write("         BFRange  rangebf.def\n")
    write_backflow_initial(workdir, definition.n_proj_bf, slater_values)
    return definition


def write_antiparallel_orbital(workdir, f_matrix):
    nsite = f_matrix.shape[0]
    lines = [
        "=============================================",
        "NOrbitalIdx         {}".format(nsite * nsite),
        "ComplexType          1",
        "=============================================",
        "=============================================",
    ]
    values = []
    idx = 0
    for i in range(nsite):
        for j in range(nsite):
            lines.append("{:5d} {:6d} {:6d}      1".format(i, j, idx))
            values.append(f_matrix[i, j])
            idx += 1
    for idx in range(nsite * nsite):
        lines.append("{:5d}      1".format(idx))
    write(os.path.join(workdir, "orbitalidx.def"), "\n".join(lines) + "\n")
    write_initial(workdir, values)
    return values


def write_real_antiparallel_orbital(workdir, f_matrix):
    nsite = f_matrix.shape[0]
    lines = [
        "=============================================",
        "NOrbitalIdx         {}".format(nsite * nsite),
        "ComplexType          0",
        "=============================================",
        "=============================================",
    ]
    values = []
    idx = 0
    for i in range(nsite):
        for j in range(nsite):
            lines.append("{:5d} {:6d} {:6d}      1".format(i, j, idx))
            values.append(float(np.real(f_matrix[i, j])))
            idx += 1
    for idx in range(nsite * nsite):
        lines.append("{:5d}      1".format(idx))
    write(os.path.join(workdir, "orbitalidx.def"), "\n".join(lines) + "\n")
    write_initial(workdir, values)
    return values


def write_general_orbital(workdir, slater_elm):
    nsite2 = slater_elm.shape[0]
    nsite = nsite2 // 2
    n_orbital = nsite2 * (nsite2 - 1) // 2
    lines = [
        "=============================================",
        "NOrbitalIdx         {}".format(n_orbital),
        "ComplexType          1",
        "=============================================",
        "=============================================",
    ]
    values = []
    idx = 0
    for orb_i in range(nsite2):
        for orb_j in range(orb_i + 1, nsite2):
            site_i = orb_i % nsite
            spin_i = orb_i // nsite
            site_j = orb_j % nsite
            spin_j = orb_j // nsite
            lines.append(
                "{:5d} {:2d} {:6d} {:2d} {:6d}      1".format(
                    site_i, spin_i, site_j, spin_j, idx
                )
            )
            # OrbitalGeneral constructs slater_elm[I,J] as
            # param[I,J] - param[J,I] = 2 * param[I,J] for I < J.
            values.append(slater_elm[orb_i, orb_j] / 2.0)
            idx += 1
    for idx in range(n_orbital):
        lines.append("{:5d}      1".format(idx))
    write(os.path.join(workdir, "orbitalidxgen.def"), "\n".join(lines) + "\n")
    write_initial(workdir, values)
    return values


def complex_n3_case(workdir):
    nsite = 6
    terms = [
        (3, (0, 0, 2, 0, 1, 0, 3, 0, 5, 1, 1, 1)),
        (3, (0, 0, 3, 0, 1, 0, 2, 0, 5, 1, 1, 1)),
        (3, (2, 0, 0, 0, 3, 0, 1, 0, 1, 1, 5, 1)),
    ]
    f_matrix = np.zeros((nsite, nsite), dtype=complex)
    for i in range(nsite):
        for j in range(nsite):
            real = math.sin(1.7 * (i + 1) * (j + 2))
            real += 0.4 * math.cos(0.2 * (i + 2) * (j + 1))
            imag = math.cos(1.1 * (i + 1) * (j + 1))
            imag -= 0.3 * math.sin((i + 3) * (j + 1))
            f_matrix[i, j] = real + 1j * imag

    slater_elm = np.zeros((2 * nsite, 2 * nsite), dtype=complex)
    for i in range(nsite):
        for j in range(nsite):
            slater_elm[i, nsite + j] = f_matrix[i, j]
            slater_elm[nsite + j, i] = -f_matrix[i, j]

    write_common_defs(workdir, nsite, 8000, 13579, "Orbital", "orbitalidx.def")
    write_nbodyg_def(workdir, terms)
    write_antiparallel_orbital(workdir, f_matrix)
    return nsite, 3, 3, slater_elm, terms


def fsz_n3_case(workdir):
    nsite = 6
    terms = [
        (3, (0, 0, 2, 1, 3, 1, 1, 0, 1, 1, 4, 1)),
        (3, (0, 0, 4, 1, 3, 1, 1, 0, 1, 1, 2, 1)),
        (3, (1, 0, 3, 1, 2, 1, 0, 0, 4, 1, 1, 1)),
    ]
    slater_elm = np.zeros((2 * nsite, 2 * nsite), dtype=complex)
    for i in range(2 * nsite):
        for j in range(i + 1, 2 * nsite):
            real = math.sin(0.7 * (i + 1) * (j + 2))
            real += 0.5 * math.cos(0.3 * (i + 2) * (j + 1))
            imag = math.cos(0.5 * (i + 3) * (j + 1))
            imag += 0.2 * math.sin((i + 1) * (j + 4))
            slater_elm[i, j] = real + 1j * imag
            slater_elm[j, i] = -slater_elm[i, j]

    write_common_defs(workdir, nsite, 8000, 24680, "OrbitalGeneral", "orbitalidxgen.def")
    write_nbodyg_def(workdir, terms)
    write_general_orbital(workdir, slater_elm)
    return nsite, 3, 3, slater_elm, terms


def projector_n3_multiterm_case(workdir):
    nsite = 6
    terms = [
        (3, (0, 0, 2, 0, 1, 0, 3, 0, 5, 1, 1, 1)),
        (3, (0, 0, 3, 0, 1, 0, 2, 0, 5, 1, 1, 1)),
        (3, (2, 0, 0, 0, 3, 0, 1, 0, 1, 1, 5, 1)),
        (3, (0, 0, 2, 0, 1, 0, 5, 0, 5, 1, 1, 1)),
    ]
    f_matrix = np.zeros((nsite, nsite), dtype=complex)
    for i in range(nsite):
        for j in range(nsite):
            real = math.sin(1.7 * (i + 1) * (j + 2))
            real += 0.4 * math.cos(0.2 * (i + 2) * (j + 1))
            imag = math.cos(1.1 * (i + 1) * (j + 1))
            imag -= 0.3 * math.sin((i + 3) * (j + 1))
            f_matrix[i, j] = real + 1j * imag

    write_common_defs(
        workdir,
        nsite,
        400,
        97531,
        "Orbital",
        "orbitalidx.def",
        data_qty=1,
        n_mp_trans=2,
    )
    write_nbodyg_def(workdir, terms)
    write_antiparallel_orbital(workdir, f_matrix)
    return terms, 1


def multibin_case(workdir):
    nsite = 6
    terms = [
        (1, (0, 0, 0, 0)),
        (2, (0, 0, 0, 0, 1, 0, 1, 0)),
        (3, (0, 0, 2, 0, 1, 0, 3, 0, 5, 1, 1, 1)),
    ]
    f_matrix = np.zeros((nsite, nsite), dtype=complex)
    for i in range(nsite):
        for j in range(nsite):
            real = math.sin(0.9 * (i + 1) * (j + 2))
            real += 0.25 * math.cos(0.4 * (i + 2) * (j + 1))
            imag = math.cos(0.8 * (i + 1) * (j + 1))
            imag -= 0.2 * math.sin((i + 3) * (j + 1))
            f_matrix[i, j] = real + 1j * imag

    write_common_defs(
        workdir,
        nsite,
        300,
        86420,
        "Orbital",
        "orbitalidx.def",
        data_qty=2,
        n_mp_trans=1,
    )
    write_nbodyg_def(workdir, terms)
    write_antiparallel_orbital(workdir, f_matrix)
    return terms, 2


def backflow_nonidentity_n3_output_case(workdir):
    nsite = 6
    terms = [
        (3, (0, 0, 2, 0, 1, 0, 3, 0, 5, 1, 1, 1)),
        (3, (0, 0, 3, 0, 1, 0, 2, 0, 5, 1, 1, 1)),
        (3, (2, 0, 0, 0, 3, 0, 1, 0, 1, 1, 5, 1)),
    ]
    f_matrix = np.zeros((nsite, nsite), dtype=complex)
    for i in range(nsite):
        for j in range(nsite):
            real = math.sin(1.7 * (i + 1) * (j + 2))
            real += 0.4 * math.cos(0.2 * (i + 2) * (j + 1))
            imag = math.cos(1.1 * (i + 1) * (j + 1))
            imag -= 0.3 * math.sin((i + 3) * (j + 1))
            f_matrix[i, j] = real + 1j * imag

    write_common_defs(workdir, nsite, 400, 54321, "Orbital", "orbitalidx.def")
    write_nbodyg_def(workdir, terms)
    slater_values = write_antiparallel_orbital(workdir, f_matrix)
    enable_nonidentity_backflow(workdir, nsite, slater_values)
    return terms, 1


def invalid_backflow_state(nsite):
    f_matrix = np.zeros((nsite, nsite), dtype=complex)
    for i in range(nsite):
        for j in range(nsite):
            real = math.sin(1.3 * (i + 1) * (j + 2))
            real += 0.35 * math.cos(0.4 * (i + 2) * (j + 1))
            imag = math.cos(0.9 * (i + 1) * (j + 1))
            imag -= 0.25 * math.sin((i + 3) * (j + 1))
            f_matrix[i, j] = real + 1j * imag

    slater_elm = np.zeros((2 * nsite, 2 * nsite), dtype=complex)
    for i in range(2 * nsite):
        for j in range(i + 1, 2 * nsite):
            value = (
                math.sin(0.7 * (i + 1) * (j + 2))
                + 1j * math.cos(0.5 * (i + 3) * (j + 1))
            )
            slater_elm[i, j] = value
            slater_elm[j, i] = -value
    return f_matrix, slater_elm


def update_modpara(workdir, keyword, value):
    path = os.path.join(workdir, "modpara.def")
    lines = []
    found = False
    with open(path) as f:
        for line in f:
            cols = line.split()
            if cols and cols[0] == keyword:
                lines.append("{:<15} {}\n".format(keyword, value))
                found = True
            else:
                lines.append(line)
    if not found:
        raise AssertionError("{} is missing from modpara.def".format(keyword))
    write(path, "".join(lines))


def invalid_backflow_real_case(workdir):
    nsite = 6
    terms = [(1, (0, 0, 2, 0))]
    f_matrix, unused_slater = invalid_backflow_state(nsite)
    write_common_defs(workdir, nsite, 20, 86420, "Orbital", "orbitalidx.def")
    write_nbodyg_def(workdir, terms)
    slater_values = write_real_antiparallel_orbital(workdir, f_matrix)
    enable_nonidentity_backflow(workdir, nsite, slater_values)


def invalid_backflow_lanczos_case(workdir):
    nsite = 6
    terms = [(1, (0, 0, 2, 0))]
    f_matrix, unused_slater = invalid_backflow_state(nsite)
    write_common_defs(workdir, nsite, 20, 86421, "Orbital", "orbitalidx.def")
    write_nbodyg_def(workdir, terms)
    slater_values = write_antiparallel_orbital(workdir, f_matrix)
    enable_nonidentity_backflow(workdir, nsite, slater_values)
    update_modpara(workdir, "NLanczosMode", 1)


def invalid_backflow_spin_change_nonfsz_case(workdir):
    nsite = 6
    terms = [(1, (0, 0, 2, 1))]
    f_matrix, unused_slater = invalid_backflow_state(nsite)
    write_common_defs(workdir, nsite, 20, 86422, "Orbital", "orbitalidx.def")
    write_nbodyg_def(workdir, terms)
    slater_values = write_antiparallel_orbital(workdir, f_matrix)
    enable_nonidentity_backflow(workdir, nsite, slater_values)


def invalid_backflow_fsz_spin_change_fixedsz_case(workdir):
    nsite = 6
    terms = [(1, (0, 0, 2, 1))]
    unused_f_matrix, slater_elm = invalid_backflow_state(nsite)
    write_common_defs(
        workdir, nsite, 20, 86423, "OrbitalGeneral", "orbitalidxgen.def"
    )
    write_nbodyg_def(workdir, terms)
    slater_values = write_general_orbital(workdir, slater_elm)
    enable_nonidentity_backflow(workdir, nsite, slater_values)


EXACT_CASES = {
    "NBodyG_Complex_N3_Oracle": complex_n3_case,
    "NBodyG_fsz_N3_Oracle": fsz_n3_case,
}

GENERATED_CHECK_CASES = {
    "NBodyG_Projector_N3_MultiTerm": projector_n3_multiterm_case,
    "NBodyG_MultiBin": multibin_case,
    "BackFlow_NBodyG_NonIdentity_N3_Output": backflow_nonidentity_n3_output_case,
}

NONVACUOUS_GENERATED_CASES = {
    "BackFlow_NBodyG_NonIdentity_N3_Output",
}

INVALID_CASES = {
    "NBodyG_InvalidBackFlow": invalid_backflow_real_case,
    "NBodyG_InvalidBackFlow_mpi": invalid_backflow_real_case,
    "NBodyG_InvalidBackFlowLanczos": invalid_backflow_lanczos_case,
    "NBodyG_InvalidBackFlowSpinChangeNonFSZ":
        invalid_backflow_spin_change_nonfsz_case,
    "NBodyG_InvalidBackFlowFszSpinChangeFixedSz":
        invalid_backflow_fsz_spin_change_fixedsz_case,
}


def parse_nbodyg_output(path):
    rows = []
    with open(path) as f:
        for line in f:
            cols = line.split()
            if not cols:
                continue
            nbody = int(cols[0])
            expected_cols = 1 + 4 * nbody + 2
            if len(cols) != expected_cols:
                raise AssertionError(
                    "{}: expected {} columns, got {}".format(path, expected_cols, len(cols))
                )
            factors = tuple(int(x) for x in cols[1:1 + 4 * nbody])
            value = complex(float(cols[-2]), float(cols[-1]))
            rows.append((nbody, factors, value))
    return rows


def assert_close(label, actual, expected):
    diff = abs(actual - expected)
    if diff > TOL:
        raise AssertionError(
            "{} mismatch: actual={} expected={} diff={} tol={}".format(
                label, actual, expected, diff, TOL
            )
        )


def assert_nbodyg_outputs(workdir, terms, bin_count, require_nonzero=False):
    found_n_ge_3 = 0
    max_abs_value = 0.0
    for idx in range(1, bin_count + 1):
        path = os.path.join(workdir, "output", "zvo_NBodyG_{:03d}.dat".format(idx))
        if not os.path.exists(path):
            raise AssertionError("missing NBodyG output: {}".format(path))
        rows = parse_nbodyg_output(path)
        if len(rows) != len(terms):
            raise AssertionError(
                "NBodyG row count mismatch in bin {}: output={} expected={}".format(
                    idx, len(rows), len(terms)
                )
            )
        for row_idx, ((nbody, factors, value), term) in enumerate(zip(rows, terms)):
            if (nbody, factors) != term:
                raise AssertionError(
                    "NBodyG bin {} term {} does not preserve input order".format(idx, row_idx)
                )
            if nbody >= 3:
                found_n_ge_3 += 1
            if not (math.isfinite(value.real) and math.isfinite(value.imag)):
                raise AssertionError("NBodyG bin {} row {} is not finite: {}".format(idx, row_idx, value))
            max_abs_value = max(max_abs_value, abs(value))

    extra_path = os.path.join(workdir, "output", "zvo_NBodyG_{:03d}.dat".format(bin_count + 1))
    if os.path.exists(extra_path):
        raise AssertionError("unexpected extra NBodyG output: {}".format(extra_path))
    if found_n_ge_3 == 0:
        raise AssertionError("generated NBodyG check did not exercise any N>=3 term")
    if require_nonzero and max_abs_value <= 1.0e-12:
        raise AssertionError("generated NBodyG output is vacuous")


def assert_finite_energy_output(workdir):
    path = os.path.join(workdir, "output", "zvo_out_001.dat")
    if not os.path.exists(path):
        raise AssertionError("missing energy output: {}".format(path))
    with open(path) as f:
        cols = f.readline().split()
    if len(cols) < 2:
        raise AssertionError("malformed energy output: {}".format(path))
    energy = complex(float(cols[0]), float(cols[1]))
    if not (math.isfinite(energy.real) and math.isfinite(energy.imag)):
        raise AssertionError("energy is not finite: {}".format(energy))


def main():
    all_cases = sorted(
        list(EXACT_CASES) + list(GENERATED_CHECK_CASES) + list(INVALID_CASES)
    )
    if len(sys.argv) not in (2, 4) or sys.argv[1] not in all_cases:
        print(
            "usage: {} {} [--expect-error <substring>]".format(
                sys.argv[0], "|".join(all_cases)
            )
        )
        return -1

    model = sys.argv[1]
    expected_error = None
    if len(sys.argv) == 4:
        if sys.argv[2] != "--expect-error" or model not in INVALID_CASES:
            print("ERROR: --expect-error is valid only for invalid cases")
            return -1
        expected_error = sys.argv[3]
    elif model in INVALID_CASES:
        print("ERROR: invalid cases require --expect-error")
        return -1

    rootdir = os.getcwd()
    workdir = os.path.join(rootdir, "work", model)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)

    if model in INVALID_CASES:
        INVALID_CASES[model](workdir)
        expected = None
        terms = None
        bin_count = 0
    elif model in EXACT_CASES:
        nsite, nup, ndown, slater_elm, terms = EXACT_CASES[model](workdir)
        expected = exact_nbody_values(nsite, nup, ndown, slater_elm, terms)
        if max(abs(x.real) for x in expected) < 1.0e-3:
            raise AssertionError("oracle terms do not exercise a nonzero real component")
        if max(abs(x.imag) for x in expected) < 1.0e-3:
            raise AssertionError("oracle terms do not exercise a nonzero imaginary component")
        bin_count = 1
    else:
        terms, bin_count = GENERATED_CHECK_CASES[model](workdir)
        expected = None

    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    if model in INVALID_CASES:
        mpi_procs = os.environ.get("MVMC_MPI_PROCS")
        if mpi_procs:
            cmd = [
                "mpirun", "-np", mpi_procs, bin_to_test,
                "-e", "namelist.def",
            ]
        else:
            cmd = [bin_to_test, "-e", "namelist.def"]
        proc = subprocess.run(
            cmd,
            cwd=workdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        with open(os.path.join(workdir, "invalid_test.log"), "w") as f:
            f.write(proc.stdout)
        if proc.returncode == 0:
            raise AssertionError("vmc.out unexpectedly accepted invalid input")
        if expected_error not in proc.stdout:
            raise AssertionError(
                "expected error substring not found: {}\n{}".format(
                    expected_error, proc.stdout
                )
            )
        if "Start: Sampling." in proc.stdout:
            raise AssertionError("invalid input reached sampling")
        if (
            mpi_procs
            and "Error: Definition files(*.def) are incomplete." in proc.stdout
        ):
            raise AssertionError(
                "MPI invalid case stopped at the rank-0 pre-broadcast gate"
            )
        return 0

    result = subprocess.call([bin_to_test, "-e", "namelist.def", "initial.def"], cwd=workdir)
    if result != 0:
        return result

    if expected is None:
        assert_nbodyg_outputs(
            workdir,
            terms,
            bin_count,
            require_nonzero=model in NONVACUOUS_GENERATED_CASES,
        )
        if model in NONVACUOUS_GENERATED_CASES:
            assert_finite_energy_output(workdir)
        print("{} generated NBodyG check passed".format(model))
        return 0

    rows = parse_nbodyg_output(os.path.join(workdir, "output", "zvo_NBodyG_001.dat"))
    if len(rows) != len(terms):
        raise AssertionError("NBodyG row count mismatch: output={} expected={}".format(len(rows), len(terms)))

    for idx, ((nbody, factors, actual), term, ref) in enumerate(zip(rows, terms, expected)):
        if (nbody, factors) != term:
            raise AssertionError("NBodyG term {} does not preserve input order".format(idx))
        if not (math.isfinite(actual.real) and math.isfinite(actual.imag)):
            raise AssertionError("NBodyG oracle row {} is not finite: {}".format(idx, actual))
        assert_close("NBodyG exact oracle row {}".format(idx), actual, ref)
    print("{} exact N=3 oracle passed".format(model))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print("ERROR: {}".format(exc))
        sys.exit(-1)
