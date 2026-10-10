"""Expert-mode checks of anti-parallel grand-canonical pairing.

The fixture is the four-site model of grandcanonical_antiparallel_oracle.py
written as OrbitalAntiParallel (or its Orbital alias) input with 2Sz=0.  The
Hamiltonian rows below are transcribed per keyword; read_input_hamiltonian()
rebuilds the 256x256 matrix from the keyword rules independently of this
writer and of the oracle's physical term list.
"""
from __future__ import print_function

import argparse
import math
import os
import shutil
import subprocess
import sys

import numpy as np

import grandcanonical_antiparallel_oracle as oracle
from runtest_gc import mpi_command, parse_sr_dump, state_dump_records


NSITE = oracle.NSITE
NORBITAL = 2 * NSITE
BAR = "=" * 45 + "\n"
RUN_TIMEOUT = 600
INVALID_TIMEOUT = 30

# Appendix A.1 rows (5-line headers omitted).  APBC flips the real and
# imaginary parts of the four (3,0)/(0,3) hopping rows only.
TRANSFER_ROWS = """0 0 0 0 0.17 0
0 1 0 1 0.17 0
1 0 1 0 0.17 0
1 1 1 1 0.17 0
2 0 2 0 0.17 0
2 1 2 1 0.17 0
3 0 3 0 0.17 0
3 1 3 1 0.17 0
0 0 1 0 0.37 0.19
1 0 0 0 0.37 -0.19
0 1 1 1 0.37 0.19
1 1 0 1 0.37 -0.19
1 0 2 0 0.37 0.19
2 0 1 0 0.37 -0.19
1 1 2 1 0.37 0.19
2 1 1 1 0.37 -0.19
2 0 3 0 0.37 0.19
3 0 2 0 0.37 -0.19
2 1 3 1 0.37 0.19
3 1 2 1 0.37 -0.19
3 0 0 0 0.37 0.19
0 0 3 0 0.37 -0.19
3 1 0 1 0.37 0.19
0 1 3 1 0.37 -0.19
"""
APBC_TRANSFER = {
    "3 0 0 0 0.37 0.19": "3 0 0 0 -0.37 -0.19",
    "0 0 3 0 0.37 -0.19": "0 0 3 0 -0.37 0.19",
    "3 1 0 1 0.37 0.19": "3 1 0 1 -0.37 -0.19",
    "0 1 3 1 0.37 -0.19": "0 1 3 1 -0.37 0.19",
}
COULOMB_INTRA_ROWS = "0 0.73\n1 -0.28\n2 0.41\n3 0.62\n"
COULOMB_INTER_ROWS = "0 1 0.21\n1 2 0.21\n2 3 0.21\n3 0 0.21\n"
HUND_ROWS = "0 1 0.17\n"
PAIRHOP_ROWS = "0 1 0.13\n"
EXCHANGE_ROWS = "0 1 -0.11\n"
INTERALL_ROWS = ("0 0 0 1 1 1 1 0 -0.16 0.09\n"
                 "1 0 1 1 0 1 0 0 -0.16 -0.09\n")
NBODY_ROWS = ("3 0 0 0 1 1 1 1 0 2 0 2 0 0.07 0\n"
              "3 2 0 2 0 1 0 1 1 0 1 0 0 0.07 0\n")
ANOMALOUS_ROWS = "1 0 0 0 1 0.35 -0.20\n0 0 1 0 0 0.35 0.20\n"

# Measurement inventory.  NBodyG and TwoBodyG mix spin-flip factors that
# cancel and include Sz-changing products whose expectation must vanish.
TWO_BODY_PAIRS = (
    ((0, 0), (1, 1)), ((0, 1), (1, 0)), ((0, 4), (4, 0)), ((0, 4), (5, 1)),
    ((1, 5), (6, 2)), ((2, 2), (6, 6)), ((3, 7), (7, 3)), ((4, 5), (5, 4)),
    ((0, 4), (1, 1)),
)
NBODY_G_TERMS = (
    ((0, 4), (5, 1), (2, 2)),
    ((2, 2), (1, 5), (4, 0)),
    ((0, 4), (1, 1), (2, 2)),
    ((1, 1), (6, 6), (3, 7)),
)
ANOMALOUS_G_ROWS = (
    (1, 0, 0, 0, 1), (1, 0, 1, 0, 0), (0, 0, 1, 0, 0), (1, 0, 0, 2, 1),
    (0, 2, 1, 0, 0), (1, 0, 0, 1, 0), (1, 0, 1, 1, 1),
)


def write(path, text):
    with open(path, "w") as stream:
        stream.write(text)


def fmt(value):
    return "{:.18e}".format(value)


def def_text(keyword, count, rows, complex_type=None):
    second = ("ComplexType {}\n".format(complex_type)
              if complex_type is not None else BAR)
    return BAR + "{} {}\n".format(keyword, count) + second + BAR + BAR + \
        "".join(rows)


def fused_to_fields(orbital):
    return orbital % NSITE, orbital // NSITE


def orbital_sign(ap, i, j):
    return -1 if ap and (i, j) in ((0, NSITE - 1), (NSITE - 1, 0)) else 1


def jastrow_class(i, j):
    return 0 if abs(i - j) in (1, NSITE - 1) else 1


def write_modpara(workdir, ap, mode, samples, seed, nstore, nsplit,
                  iterations, init_nelec, data_qty, interval, two_sz,
                  extra=""):
    write(os.path.join(workdir, "modpara.def"), """--------------------
Model_Parameters   0
--------------------
VMC_Cal_Parameters
--------------------
CDataFileHead  zvo
CParaFileHead  zqp
--------------------
NVMCCalMode    {mode}
NLanczosMode   0
--------------------
NDataIdxStart  1
NDataQtySmp    {data_qty}
--------------------
Nsite          {nsite}
Ne             {half}
Ncond          {nsite}
2Sz            {two_sz}
NSPGaussLeg    1
NSPStot        0
NMPTrans       {nmptrans}
NSROptItrStep  {iterations}
NSROptItrSmp   1
DSROptRedCut   0.000000000001
DSROptStaDel   0.01
DSROptStepDt   0.002
NVMCWarmUp     2000
NVMCInterval   {interval}
NVMCSample     {samples}
NExUpdatePath  0
RndSeed        {seed}
NSplitSize     {nsplit}
NStore         {nstore}
NSRCG          0
NGrandCanonical 1
NGCInitNelec   {init_nelec}
{extra}""".format(mode=mode, data_qty=data_qty, nsite=NSITE, half=NSITE // 2,
                  two_sz=two_sz, nmptrans=-1 if ap else 1,
                  iterations=iterations, interval=interval, samples=samples,
                  seed=seed, nsplit=nsplit, nstore=nstore,
                  init_nelec=init_nelec, extra=extra))


def antiparallel_orbital_text(ap, optimize):
    rows = []
    for i in range(NSITE):
        for j in range(NSITE):
            index = NSITE * i + j
            if ap:
                rows.append("{} {} {} {}\n".format(
                    i, j, index, orbital_sign(ap, i, j)))
            else:
                rows.append("{} {} {}\n".format(i, j, index))
    rows.extend("{} {}\n".format(index, 1 if optimize else 0)
                for index in range(NSITE * NSITE))
    return def_text("NOrbitalIdx", NSITE * NSITE, rows, complex_type=1)


def general_pairs():
    return [(first, second) for first in range(NORBITAL)
            for second in range(first + 1, NORBITAL)]


def general_orbital_text(ap, optimize):
    rows = []
    for index, (first, second) in enumerate(general_pairs()):
        site0, spin0 = fused_to_fields(first)
        site1, spin1 = fused_to_fields(second)
        sign = 1
        if spin0 == 0 and spin1 == 1:
            sign = orbital_sign(ap, site0, site1)
        rows.append("{} {} {} {} {} {}\n".format(
            site0, spin0, site1, spin1, index, sign))
    rows.extend("{} {}\n".format(index, 1 if optimize else 0)
                for index in range(len(general_pairs())))
    return def_text("NOrbitalIdx", len(general_pairs()), rows, complex_type=1)


def orbital_parameters(representation):
    """Raw orbital parameters in OptFlag order."""
    if representation == "anti":
        return [oracle.F0[i, j] for i in range(NSITE) for j in range(NSITE)]
    values = []
    for first, second in general_pairs():
        site0, spin0 = fused_to_fields(first)
        site1, spin1 = fused_to_fields(second)
        if spin0 == 0 and spin1 == 1:
            values.append(oracle.F0[site0, site1] / 2.0)
        else:
            values.append(0.0j)
    return values


def write_projection(workdir, optimize):
    write(os.path.join(workdir, "gutzwilleridx.def"),
          def_text("NGutzwillerIdx", 1,
                   ["{} 0\n".format(site) for site in range(NSITE)] +
                   ["0 {}\n".format(1 if optimize else 0)], complex_type=0))
    rows = ["{} {} {}\n".format(i, j, jastrow_class(i, j))
            for i in range(NSITE) for j in range(NSITE) if i != j]
    rows.extend(["0 {}\n".format(1 if optimize else 0), "1 0\n"])
    write(os.path.join(workdir, "jastrowidx.def"),
          def_text("NJastrowIdx", 2, rows, complex_type=0))
    write(os.path.join(workdir, "locspn.def"),
          def_text("NlocalSpin", 0,
                   ["{} 0\n".format(site) for site in range(NSITE)]))
    write(os.path.join(workdir, "qptransidx.def"),
          def_text("NQPTrans", 1,
                   ["0 1.0 0.0\n"] +
                   ["0 {} {} 1\n".format(site, site) for site in range(NSITE)]))


def transfer_text(ap):
    lines = TRANSFER_ROWS.splitlines()
    if ap:
        lines = [APBC_TRANSFER.get(line, line) for line in lines]
    return "".join(line + "\n" for line in lines)


def write_hamiltonian(workdir, ap):
    files = (
        ("trans.def", "NTransfer", transfer_text(ap)),
        ("coulombintra.def", "NCoulombIntra", COULOMB_INTRA_ROWS),
        ("coulombinter.def", "NCoulombInter", COULOMB_INTER_ROWS),
        ("hund.def", "NHund", HUND_ROWS),
        ("pairhop.def", "NPairHop", PAIRHOP_ROWS),
        ("exchange.def", "NExchange", EXCHANGE_ROWS),
        ("interall.def", "NInterAll", INTERALL_ROWS),
        ("nbodyinterall.def", "NNBodyInterAll", NBODY_ROWS),
        ("anomalousterm.def", "NAnomalousTerm", ANOMALOUS_ROWS),
    )
    for name, keyword, body in files:
        rows = [line + "\n" for line in body.splitlines() if line.strip()]
        write(os.path.join(workdir, name), def_text(keyword, len(rows), rows))


def factors_text(pairs):
    fields = []
    for out_orbital, in_orbital in pairs:
        fields.extend(fused_to_fields(out_orbital))
        fields.extend(fused_to_fields(in_orbital))
    return " ".join(str(value) for value in fields)


def write_measurements(workdir, mode):
    one_rows = ["{}\n".format(factors_text(((out_orbital, in_orbital),)))
                for out_orbital in range(NORBITAL)
                for in_orbital in range(NORBITAL)]
    write(os.path.join(workdir, "greenone.def"),
          def_text("NCisAjs", len(one_rows), one_rows))
    two_rows = ["{}\n".format(factors_text(pairs)) for pairs in TWO_BODY_PAIRS]
    write(os.path.join(workdir, "greentwo.def"),
          def_text("NCisAjsCktAltDC", len(two_rows), two_rows))
    nbody_rows = ["{} {}\n".format(len(term), factors_text(term))
                  for term in NBODY_G_TERMS]
    write(os.path.join(workdir, "nbodyg.def"),
          def_text("NNBodyG", len(nbody_rows), nbody_rows))
    names = ["        OneBodyG  greenone.def", "        TwoBodyG  greentwo.def",
             "          NBodyG  nbodyg.def"]
    if mode == 1:
        rows = ["{} {} {} {} {}\n".format(*row) for row in ANOMALOUS_G_ROWS]
        write(os.path.join(workdir, "anomalousg.def"),
              def_text("NAnomalousG", len(rows), rows))
        names.append("      AnomalousG  anomalousg.def")
    return names


def write_initial(workdir, representation, gutz=oracle.GUTZ,
                  jastrow=oracle.JASTROW, parameters=None):
    values = [0.0] * 6
    for value in (gutz, jastrow, 0.0):
        values.extend((value, 0.0, 0.0))
    if parameters is None:
        parameters = orbital_parameters(representation)
    for value in parameters:
        value = complex(value)
        values.extend((value.real, value.imag, 0.0))
    write(os.path.join(workdir, "initial.def"),
          " ".join(fmt(value) for value in values) + "\n")


def write_fixture(workdir, ap=False, mode=1, samples=60000, seed=91807,
                  nstore=1, nsplit=1, init_nelec=4, representation="anti",
                  alias="OrbitalAntiParallel", iterations=1, optimize=False,
                  data_qty=1, interval=4, two_sz=None, measurements=True):
    if two_sz is None:
        two_sz = 0 if representation == "anti" else -1
    write_modpara(workdir, ap, mode, samples, seed, nstore, nsplit,
                  iterations, init_nelec, data_qty, interval, two_sz)
    write_projection(workdir, optimize)
    write_hamiltonian(workdir, ap)
    namelist = [
        "         ModPara  modpara.def",
        "         LocSpin  locspn.def",
        "      Gutzwiller  gutzwilleridx.def",
        "         Jastrow  jastrowidx.def",
        "        TransSym  qptransidx.def",
        "           Trans  trans.def",
        "    CoulombIntra  coulombintra.def",
        "    CoulombInter  coulombinter.def",
        "            Hund  hund.def",
        "         PairHop  pairhop.def",
        "        Exchange  exchange.def",
        "        InterAll  interall.def",
        "   NBodyInterAll  nbodyinterall.def",
        "   AnomalousTerm  anomalousterm.def",
    ]
    if representation == "anti":
        write(os.path.join(workdir, "orbitalidx.def"),
              antiparallel_orbital_text(ap, optimize))
        namelist.append("{:>16}  orbitalidx.def".format(alias))
    else:
        write(os.path.join(workdir, "orbitalidxgen.def"),
              general_orbital_text(ap, optimize))
        namelist.append("  OrbitalGeneral  orbitalidxgen.def")
    if measurements:
        namelist.extend(write_measurements(workdir, mode))
    write(os.path.join(workdir, "namelist.def"), "\n".join(namelist) + "\n")
    write_initial(workdir, representation)


# ---------------------------------------------------------------------------
# Independent read-back of the Hamiltonian keyword rules (calham_gc.c and the
# readdef.c expansions).  Shares no row generation with the writer above.
# ---------------------------------------------------------------------------
def read_namelist(workdir):
    entries = {}
    with open(os.path.join(workdir, "namelist.def")) as stream:
        for line in stream:
            fields = line.split()
            if len(fields) >= 2:
                entries[fields[0]] = os.path.join(workdir, fields[1])
    return entries


def read_body(path):
    with open(path) as stream:
        lines = stream.readlines()
    declared = int(lines[1].split()[1])
    rows = [line.split() for line in lines[5:] if line.split()]
    if len(rows) != declared:
        raise AssertionError("{}: header count {} != body rows {}".format(
            path, declared, len(rows)))
    return rows


def spin_orbital(site, spin):
    return int(site) + int(spin) * NSITE


def one(out_orbital, in_orbital):
    return (("c", out_orbital), ("a", in_orbital))


def read_input_hamiltonian(workdir, nsite=NSITE):
    if nsite != NSITE:
        raise ValueError("read-back is defined for the four-site fixture")
    entries = read_namelist(workdir)
    terms = []
    for row in read_body(entries["Trans"]):
        i, s, j, t = row[:4]
        value = complex(float(row[4]), float(row[5]))
        terms.append((-value, one(spin_orbital(i, s), spin_orbital(j, t))))
    for row in read_body(entries["CoulombIntra"]):
        i = int(row[0])
        terms.append((float(row[1]), one(i, i) + one(i + nsite, i + nsite)))
    for row in read_body(entries["CoulombInter"]):
        i, j, value = int(row[0]), int(row[1]), float(row[2])
        for a in (i, i + nsite):
            for b in (j, j + nsite):
                terms.append((value, one(a, a) + one(b, b)))
    for row in read_body(entries["Hund"]):
        i, j, value = int(row[0]), int(row[1]), float(row[2])
        terms.append((-value, one(i, i) + one(j, j)))
        terms.append((-value, one(i + nsite, i + nsite) +
                      one(j + nsite, j + nsite)))
    for row in read_body(entries["PairHop"]):
        # readdef.c stores (i,j,J) and (j,i,J); each is J C(i,j)up C(i,j)dn.
        i, j, value = int(row[0]), int(row[1]), float(row[2])
        for a, b in ((i, j), (j, i)):
            terms.append((value, one(a, b) + one(a + nsite, b + nsite)))
    for row in read_body(entries["Exchange"]):
        i, j, value = int(row[0]), int(row[1]), float(row[2])
        terms.append((value, one(i, j) + one(j + nsite, i + nsite)))
        terms.append((value, one(i + nsite, j + nsite) + one(j, i)))
    for row in read_body(entries["InterAll"]):
        fields = [int(value) for value in row[:8]]
        value = complex(float(row[8]), float(row[9]))
        terms.append((value, one(spin_orbital(fields[0], fields[1]),
                                 spin_orbital(fields[2], fields[3])) +
                      one(spin_orbital(fields[4], fields[5]),
                          spin_orbital(fields[6], fields[7]))))
    for row in read_body(entries["NBodyInterAll"]):
        count = int(row[0])
        fields = [int(value) for value in row[1:1 + 4 * count]]
        value = complex(float(row[1 + 4 * count]), float(row[2 + 4 * count]))
        operators = ()
        for k in range(count):
            operators += one(spin_orbital(fields[4 * k], fields[4 * k + 1]),
                             spin_orbital(fields[4 * k + 2],
                                          fields[4 * k + 3]))
        terms.append((value, operators))
    for row in read_body(entries["AnomalousTerm"]):
        kind = int(row[0])
        first = spin_orbital(row[1], row[2])
        second = spin_orbital(row[3], row[4])
        value = complex(float(row[5]), float(row[6]))
        letter = "c" if kind == 1 else "a"
        terms.append((value, ((letter, first), (letter, second))))
    dimension = 1 << (2 * nsite)
    matrix = np.zeros((dimension, dimension), dtype=np.complex128)
    for coefficient, operators in terms:
        for source in range(dimension):
            out = oracle.apply_ops(source, operators)
            if out is not None:
                target, sign = out
                matrix[target, source] += coefficient * sign
    return matrix


def check_hamiltonian(workdir, ap):
    h_read = read_input_hamiltonian(workdir, nsite=NSITE)
    h_model = oracle.hamiltonian_from_terms(oracle.model_terms(ap), NSITE)
    if not (h_read.shape == h_model.shape == (256, 256)):
        raise AssertionError("Hamiltonian shape")
    if not (np.isfinite(h_read).all() and np.isfinite(h_model).all()):
        raise AssertionError("nonfinite Hamiltonian entry")
    difference = float(np.max(abs(h_read - h_model)))
    hermitian = float(np.max(abs(h_read - h_read.conj().T)))
    return difference, hermitian


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------
def binary_path(rootdir):
    return os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")


def execute(rootdir, workdir, procs=1, extra_env=None, timeout=RUN_TIMEOUT,
            options=()):
    command = [binary_path(rootdir)] + list(options) + \
        ["-e", "namelist.def", "initial.def"]
    if procs > 1:
        command = mpi_command(procs, command)
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = environment.get("OMP_NUM_THREADS", "1")
    environment["OPENBLAS_NUM_THREADS"] = "1"
    environment["MKL_NUM_THREADS"] = "1"
    if extra_env:
        environment.update(extra_env)
    try:
        process = subprocess.run(command, cwd=workdir, env=environment,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 universal_newlines=True, timeout=timeout)
    except subprocess.TimeoutExpired as error:
        output = error.stdout or ""
        if isinstance(output, bytes):
            output = output.decode("utf-8", "replace")
        write(os.path.join(workdir, "run.log"), output)
        raise AssertionError("vmc.out timed out after {} s in {}".format(
            timeout, workdir))
    write(os.path.join(workdir, "run.log"), process.stdout)
    return process.returncode, process.stdout


def run_binary(rootdir, workdir, procs=1, extra_env=None):
    returncode, output = execute(rootdir, workdir, procs, extra_env)
    if returncode != 0 or "Finish calculation." not in output:
        raise AssertionError("vmc.out failed (rc={}) in {}:\n{}".format(
            returncode, workdir, output[-6000:]))
    return output


def prepare_work(rootdir, name):
    workdir = os.path.join(rootdir, "work_gc_antiparallel", name)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    return workdir


def sample_records(workdir, procs=1, label="SAMPLE"):
    records = []
    for rank in range(procs):
        path = os.path.join(workdir, "state.dat.rank{}".format(rank))
        records.extend((rank,) + record for record in state_dump_records(path)
                       if record[0] == label)
    return records


def check_balanced(records):
    for rank, label, index, ncur, mask, unused in records:
        up = bin(mask & ((1 << NSITE) - 1)).count("1")
        down = bin(mask >> NSITE).count("1")
        if not (up == down == ncur // 2 and ncur == up + down):
            raise AssertionError("{} {} rank {} has mask {} ncur {}".format(
                label, index, rank, mask, ncur))


# ---------------------------------------------------------------------------
# Input mutations (Task 4 corpus on the production fixture, Sz rules, gates)
# ---------------------------------------------------------------------------
def edit_file(workdir, name, function):
    path = os.path.join(workdir, name)
    with open(path) as stream:
        lines = stream.readlines()
    with open(path, "w") as stream:
        stream.writelines(function(lines))


def replace_line(workdir, name, index, text):
    def change(lines):
        lines[index] = text
        return lines
    edit_file(workdir, name, change)


def append_rows(workdir, name, rows, count_delta=None):
    def change(lines):
        if count_delta is not None:
            fields = lines[1].split()
            lines[1] = "{} {}\n".format(fields[0], int(fields[1]) + count_delta)
        return lines + rows
    edit_file(workdir, name, change)


def set_header_count(workdir, name, value):
    def change(lines):
        lines[1] = "{} {}\n".format(lines[1].split()[0], value)
        return lines
    edit_file(workdir, name, change)


def update_modpara(workdir, key, value):
    def change(lines):
        out = []
        found = False
        for line in lines:
            fields = line.split()
            if fields and fields[0] == key:
                found = True
                if value is not None:
                    out.append("{:<15}{}\n".format(key, value))
            else:
                out.append(line)
        if not found and value is not None:
            out.append("{:<15}{}\n".format(key, value))
        return out
    edit_file(workdir, "modpara.def", change)


def add_namelist(workdir, keyword, filename, text=None):
    if text is not None:
        write(os.path.join(workdir, filename), text)
    with open(os.path.join(workdir, "namelist.def"), "a") as stream:
        stream.write("{:>16}  {}\n".format(keyword, filename))


ORBITAL = "orbitalidx.def"
PAIR0 = 5          # first pair row (after the five header lines)
OPT0 = PAIR0 + NSITE * NSITE


def pair_line(ap, i, j, index=None, sign=None):
    index = NSITE * i + j if index is None else index
    if sign is None:
        sign = orbital_sign(ap, i, j)
    return "{} {} {} {}\n".format(i, j, index, sign) if ap else \
        "{} {} {}\n".format(i, j, index)


# name -> (function(workdir, ap), expectation)
# expectation: list of substrings that must appear for a rejected input, or
# None when the input must run to completion.  A dict keys the expectation
# by boundary.
def _m(function, expected):
    return function, expected


def _orbital_line(index, text):
    return lambda workdir, ap: replace_line(workdir, ORBITAL, index, text)


def _orbital_edit(function):
    return lambda workdir, ap: edit_file(workdir, ORBITAL, function)


INPUT_MUTATIONS = {
    # header
    "header_missing": _m(_orbital_edit(lambda lines: lines[:4]),
                         ["GC anti-parallel", "header"]),
    "header_norb_zero": _m(_orbital_line(1, "NOrbitalIdx 0\n"),
                           ["orbital count"]),
    "header_norb_negative": _m(_orbital_line(1, "NOrbitalIdx -16\n"),
                               ["orbital count"]),
    "header_norb_int_overflow": _m(_orbital_line(1, "NOrbitalIdx 2147483648\n"),
                                   ["orbital count"]),
    "header_norb_llong_overflow": _m(
        _orbital_line(1, "NOrbitalIdx 9223372036854775808\n"), ["malformed"]),
    "header_norb_decimal": _m(_orbital_line(1, "NOrbitalIdx 16.0\n"),
                              ["malformed"]),
    "header_complex_type": _m(_orbital_line(2, "ComplexType 2\n"),
                              ["ComplexType"]),
    # pair rows
    "pair_missing": _m(_orbital_edit(lambda lines: lines[:PAIR0] +
                                     lines[PAIR0 + 1:]), ["columns"]),
    "pair_surplus": _m(_orbital_edit(lambda lines: lines[:OPT0] +
                                     [lines[PAIR0]] + lines[OPT0:]),
                       ["columns"]),
    "pair_duplicate": _m(lambda workdir, ap: replace_line(
        workdir, ORBITAL, PAIR0 + 1, pair_line(ap, 0, 0)), ["duplicate pair"]),
    "pair_site_range": _m(lambda workdir, ap: replace_line(
        workdir, ORBITAL, PAIR0, pair_line(ap, 4, 0, 0, 1)), ["site index"]),
    "pair_negative_site": _m(lambda workdir, ap: replace_line(
        workdir, ORBITAL, PAIR0, pair_line(ap, -1, 0, 0, 1)), ["site index"]),
    "pair_param_range": _m(lambda workdir, ap: replace_line(
        workdir, ORBITAL, PAIR0, pair_line(ap, 0, 0, 16)),
        ["parameter index"]),
    "pair_two_columns": _m(_orbital_line(PAIR0, "0 0\n"), ["columns"]),
    "pair_five_columns": _m(_orbital_line(PAIR0, "0 0 0 1 1\n"), ["columns"]),
    "pair_long_line": _m(_orbital_line(PAIR0, "0 0 0 1" + " " * 5000 + "\n"),
                         ["too long"]),
    "pair_blank_middle": _m(_orbital_edit(lambda lines: lines[:PAIR0 + 2] +
                                          ["\n"] + lines[PAIR0 + 2:]),
                            ["blank line"]),
    "pair_comment": _m(_orbital_edit(lambda lines: lines[:PAIR0] +
                                     ["# comment\n"] + lines[PAIR0:]),
                       ["malformed"]),
    # fourth column: required +/-1 under APBC, any int ignored under PBC
    "sign_missing": _m(_orbital_line(PAIR0, "0 0 0\n"),
                       {"apbc": ["4 integer columns"], "pbc": None}),
    "sign_zero": _m(_orbital_line(PAIR0, "0 0 0 0\n"),
                    {"apbc": ["sign must be +1 or -1"], "pbc": None}),
    "sign_two": _m(_orbital_line(PAIR0, "0 0 0 2\n"),
                   {"apbc": ["sign must be +1 or -1"], "pbc": None}),
    "sign_decimal": _m(_orbital_line(PAIR0, "0 0 0 1.5\n"), ["malformed"]),
    "sign_int_overflow": _m(_orbital_line(PAIR0, "0 0 0 2147483648\n"),
                            ["int range"]),
    # OptFlag rows
    "opt_missing": _m(_orbital_edit(lambda lines: lines[:-1]),
                      ["ended after 15 of 16 OptFlag"]),
    "opt_surplus": _m(_orbital_edit(lambda lines: lines + ["0 1\n"]),
                      ["after the OptFlag rows"]),
    "opt_duplicate": _m(_orbital_line(OPT0 + 1, "0 0\n"),
                        ["duplicate OptFlag"]),
    "opt_negative_index": _m(_orbital_line(OPT0, "-1 0\n"), ["OptFlag index"]),
    "opt_upper_index": _m(_orbital_line(OPT0, "16 0\n"), ["OptFlag index"]),
    "opt_flag_minus1": _m(_orbital_line(OPT0, "0 -1\n"),
                          ["flag must be 0 or 1"]),
    "opt_flag_two": _m(_orbital_line(OPT0, "0 2\n"), ["flag must be 0 or 1"]),
    "opt_decimal": _m(_orbital_line(OPT0, "0 0.0\n"), ["malformed"]),
    "opt_three_columns": _m(_orbital_line(OPT0, "0 0 0\n"),
                            ["2 integer columns"]),
    "opt_blank_middle": _m(_orbital_edit(lambda lines: lines[:OPT0 + 1] +
                                         ["\n"] + lines[OPT0 + 1:]),
                           ["blank line"]),
    "orbital_file_missing": _m(lambda workdir, ap: os.remove(
        os.path.join(workdir, ORBITAL)), ["Broken file or Not exist"]),
    # Hamiltonian Sz rules (whole-term)
    "sz_transfer_flip": _m(lambda workdir, ap: append_rows(
        workdir, "trans.def", ["0 0 1 1 0.05 0\n"], 1),
        ["Sz-conserving", "Trans row 25"]),
    "sz_interall_net2": _m(lambda workdir, ap: append_rows(
        workdir, "interall.def", ["0 0 0 1 1 0 1 1 0.05 0\n"], 1),
        ["Sz-conserving", "InterAll row 3"]),
    "sz_nbody_net1": _m(lambda workdir, ap: append_rows(
        workdir, "nbodyinterall.def", ["3 0 0 0 1 1 0 1 0 2 0 2 0 0.05 0\n"],
        1), ["Sz-conserving", "NBodyInterAll row 3"]),
    "sz_anomalous_same_spin": _m(lambda workdir, ap: append_rows(
        workdir, "anomalousterm.def",
        ["1 0 0 1 0 0.10 0.00\n", "0 1 0 0 0 0.10 0.00\n"], 2),
        ["Sz-conserving", "AnomalousTerm row 3"]),
    # accepted inputs
    "zero_nonconserving_rows": _m(lambda workdir, ap: (
        append_rows(workdir, "trans.def", ["0 0 1 1 0 0\n"], 1),
        append_rows(workdir, "interall.def", ["0 0 0 1 1 0 1 1 0 0\n"], 1),
        append_rows(workdir, "anomalousterm.def",
                    ["1 0 0 1 0 0 0\n", "0 1 0 0 0 0 0\n"], 2)), None),
    "alias_orbital": _m(lambda workdir, ap: edit_file(
        workdir, "namelist.def", lambda lines: [
            line.replace("OrbitalAntiParallel", "            Orbital")
            for line in lines]), None),
    "start_vacuum": _m(lambda workdir, ap: update_modpara(
        workdir, "NGCInitNelec", 0), None),
    "start_full": _m(lambda workdir, ap: update_modpara(
        workdir, "NGCInitNelec", 2 * NSITE), None),
    "start_default": _m(lambda workdir, ap: update_modpara(
        workdir, "NGCInitNelec", None), None),
    "start_odd": _m(lambda workdir, ap: update_modpara(
        workdir, "NGCInitNelec", 5), None),
    "ncond_ignored": _m(lambda workdir, ap: update_modpara(
        workdir, "Ncond", 7), None),
    # mode gates
    "twosz_missing": _m(lambda workdir, ap: update_modpara(
        workdir, "2Sz", None), ["requires 2Sz=0", "got -1"]),
    "twosz_minus1": _m(lambda workdir, ap: update_modpara(workdir, "2Sz", -1),
                       ["requires 2Sz=0"]),
    "twosz_two": _m(lambda workdir, ap: update_modpara(workdir, "2Sz", 2),
                    ["requires 2Sz=0", "got 2"]),
    "alias_duplicate": _m(lambda workdir, ap: add_namelist(
        workdir, "Orbital", ORBITAL), ["must not be duplicated"]),
    "all_real": _m(lambda workdir, ap: replace_line(
        workdir, ORBITAL, 2, "ComplexType 0\n"),
        ["requires complex variational parameters"]),
    "gaussleg2": _m(lambda workdir, ap: update_modpara(
        workdir, "NSPGaussLeg", 2), ["NSPGaussLeg=1"]),
    "mptrans2": _m(lambda workdir, ap: update_modpara(
        workdir, "NMPTrans", 2), ["requires NMPTrans=+1 or -1"]),
    "opttrans": _m(lambda workdir, ap: add_namelist(
        workdir, "OptTrans", "opttrans.def",
        BAR + "NQPOptTrans 1\n" + BAR + BAR + BAR + "0 1.0\n" +
        "".join("0 {} {}\n".format(site, site) for site in range(NSITE))),
        ["does not support OptTrans"]),
    "backflow": _m(lambda workdir, ap: (
        add_namelist(workdir, "BF", "bf.def",
                     "====================\nNBackFlowIdx 1\n"
                     "====================\n"),
        add_namelist(workdir, "BFRange", "rangebf.def",
                     "====================\nNrange 1 1\n"
                     "====================\n")),
        ["does not support BackFlow"]),
    "rbm": _m(lambda workdir, ap: add_namelist(
        workdir, "GeneralRBM_HiddenLayer", "rbm.def",
        BAR + "NRBM_HiddenLayerIdx 1\nComplexType          1\n" + BAR),
        ["does not support RBM"]),
    "locspin": _m(lambda workdir, ap: edit_file(
        workdir, "locspn.def", lambda lines: lines[:1] +
        ["NlocalSpin 1\n"] + lines[2:5] + ["0 1\n"] + lines[6:]),
        ["does not support LocSpin"]),
    "lanczos": _m(lambda workdir, ap: update_modpara(
        workdir, "NLanczosMode", 1), ["does not support Lanczos"]),
    "srcg": _m(lambda workdir, ap: update_modpara(workdir, "NSRCG", 1),
               ["supports dense SR only"]),
    "exupdate": _m(lambda workdir, ap: update_modpara(
        workdir, "NExUpdatePath", 1), ["requires NExUpdatePath=0"]),
    "updateweight": _m(lambda workdir, ap: add_namelist(
        workdir, "InUpdateWeight", "updateweight.def",
        "Exchange 1.0\nLocalSpinFlip 1.0\nPairSpinFlip 1.0\n"),
        ["does not support the UpdateWeight input"]),
}
OPTIONS = {"opttrans": ("-o",)}
SIGNAL_MARKERS = ("Segmentation fault", "signal 11", "Signal: ",
                  "AddressSanitizer", "runtime error:", "Abort trap")


def input_expectation(mutation, boundary):
    expected = INPUT_MUTATIONS[mutation][1]
    if isinstance(expected, dict):
        return expected[boundary]
    return expected


def run_invalid(rootdir, mutation, ap, procs):
    boundary = "apbc" if ap else "pbc"
    workdir = prepare_work(rootdir, "input_{}_{}_np{}".format(
        mutation, boundary, procs))
    write_fixture(workdir, ap=ap, mode=1, samples=16, iterations=1)
    INPUT_MUTATIONS[mutation][0](workdir, ap)
    expected = input_expectation(mutation, boundary)
    returncode, output = execute(
        rootdir, workdir, procs,
        extra_env={"MVMC_GC_STATE_DUMP": "state.dat"},
        timeout=INVALID_TIMEOUT if expected is not None else RUN_TIMEOUT,
        options=OPTIONS.get(mutation, ()))
    signals = [marker for marker in SIGNAL_MARKERS if marker in output]
    if signals:
        raise AssertionError("{}: signal/sanitizer report {}:\n{}".format(
            mutation, signals, output[-4000:]))
    if expected is None:
        if returncode != 0 or "Finish calculation." not in output:
            raise AssertionError("{} must be accepted (rc={}):\n{}".format(
                mutation, returncode, output[-4000:]))
        check_balanced(sample_records(workdir, procs))
        return "accepted"
    if returncode == 0 or returncode < 0:
        raise AssertionError("{} must fail with a nonzero exit (rc={}):\n{}"
                             .format(mutation, returncode, output[-4000:]))
    missing = [text for text in expected if text not in output]
    if missing:
        raise AssertionError("{} diagnostic {} missing:\n{}".format(
            mutation, missing, output[-4000:]))
    return "rejected"


def input_case(rootdir, args):
    if args.mutation not in INPUT_MUTATIONS:
        raise SystemExit("unknown input mutation {}".format(args.mutation))
    verdict = run_invalid(rootdir, args.mutation, args.boundary == "apbc",
                          args.np)
    print("GC anti-parallel input {} {} np={}: {} as expected".format(
        args.mutation, args.boundary, args.np, verdict))


def hamiltonian_case(rootdir, args):
    ap = args.boundary == "apbc"
    workdir = prepare_work(rootdir, "hamiltonian_{}".format(args.boundary))
    write_fixture(workdir, ap=ap)
    difference, hermitian = check_hamiltonian(workdir, ap)
    if not (difference < 1e-13 and hermitian < 1e-13):
        raise AssertionError("read-back H differs by {} (hermitian {})".format(
            difference, hermitian))
    print("GC anti-parallel Hamiltonian read-back {}: max diff {:.3g}".format(
        args.boundary, difference))


# ---------------------------------------------------------------------------
# Audit
# ---------------------------------------------------------------------------
def parse_audit(path):
    meta = {}
    orbitals = []
    with open(path) as stream:
        for line in stream:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "ORBITAL":
                orbitals.append(tuple(int(value) for value in fields[1:]))
            elif fields[0] not in ("TRANS", "TRANSWEIGHT"):
                meta[fields[0]] = fields[1]
    return meta, orbitals


def audit_case(rootdir, args):
    ap = args.boundary == "apbc"
    tag = "{}_np{}".format(args.boundary, args.np)
    workdir = prepare_work(rootdir, "audit_" + tag)
    write_fixture(workdir, ap=ap, mode=1, samples=16, measurements=False)
    run_binary(rootdir, workdir, args.np,
               extra_env={"MVMC_GC_INPUT_AUDIT": "audit.txt"})
    meta, orbitals = parse_audit(os.path.join(workdir, "audit.txt"))
    expected_negative = 2 if ap else 0
    if (meta.get("orbital_mode") != "antiparallel" or
            meta.get("nsite") != str(NSITE) or
            meta.get("nslater") != str(NSITE * NSITE) or
            meta.get("ap_flag") != str(int(ap)) or
            meta.get("negative_orbital_input_sign_count") !=
            str(expected_negative)):
        raise AssertionError("audit metadata {}".format(meta))
    expected = [(i, j + NSITE, NSITE * i + j, orbital_sign(ap, i, j))
                for i in range(NSITE) for j in range(NSITE)]
    if orbitals != expected:
        raise AssertionError("audit orbital rows {}".format(orbitals))
    # An empty variable writes nothing; open and write failures stop all ranks.
    workdir = prepare_work(rootdir, "audit_empty_" + tag)
    write_fixture(workdir, ap=ap, mode=1, samples=16, measurements=False)
    run_binary(rootdir, workdir, args.np, extra_env={"MVMC_GC_INPUT_AUDIT": ""})
    if [name for name in os.listdir(workdir) if "audit" in name]:
        raise AssertionError("empty MVMC_GC_INPUT_AUDIT produced a file")
    failures = [("missing_dir", os.path.join("missing", "audit.txt"),
                 "failed to open MVMC_GC_INPUT_AUDIT")]
    if os.path.exists("/dev/full"):
        failures.append(("dev_full", "/dev/full",
                         "failed to write MVMC_GC_INPUT_AUDIT"))
    for name, path, message in failures:
        workdir = prepare_work(rootdir, "audit_{}_{}".format(name, tag))
        write_fixture(workdir, ap=ap, mode=1, samples=16, measurements=False)
        returncode, output = execute(rootdir, workdir, args.np,
                                     extra_env={"MVMC_GC_INPUT_AUDIT": path},
                                     timeout=INVALID_TIMEOUT)
        if returncode <= 0 or message not in output:
            raise AssertionError("audit {} rc={}:\n{}".format(
                name, returncode, output[-3000:]))
    print("GC anti-parallel audit {} passed".format(tag))


# ---------------------------------------------------------------------------
# Smoke
# ---------------------------------------------------------------------------
def read_rows(path):
    with open(path) as stream:
        return [line.split() for line in stream if line.split()]


def check_finite_outputs(workdir, mode, iterations):
    gc_rows = read_rows(os.path.join(workdir, "zvo_gc.dat"))
    out_rows = read_rows(os.path.join(workdir, "output", "zvo_out_001.dat"))
    expected_rows = iterations if mode == 0 else 1
    if len(gc_rows) != expected_rows or len(out_rows) != expected_rows:
        raise AssertionError("output rows gc={} out={} expected {}".format(
            len(gc_rows), len(out_rows), expected_rows))
    for row in gc_rows + out_rows:
        if not all(math.isfinite(float(value)) for value in row):
            raise AssertionError("nonfinite output row {}".format(row))
    for row in out_rows:
        if float(row[4]) != 0.0 or float(row[5]) != 0.0:
            raise AssertionError("Sz columns must vanish: {}".format(row))


def check_burn(workdir, procs):
    for rank in range(procs):
        path = os.path.join(workdir, "state.dat.rank{}".format(rank))
        records = state_dump_records(path)
        stores = [record for record in records if record[0] == "STORE"]
        restores = [record for record in records if record[0] == "RESTORE"]
        if len(stores) < len(restores):
            raise AssertionError("rank {} restores without a store".format(
                rank))
        for store, restore in zip(stores, restores):
            if store[2] != restore[2] or store[3] != restore[3] or \
                    store[4][1:] != restore[4][1:] or \
                    store[4][0].split()[3:] != restore[4][0].split()[3:]:
                raise AssertionError("rank {} burn state changed".format(rank))
        check_balanced([(rank,) + record for record in stores + restores])


def compare_output_files(label, actual_dir, reference_dir, mode):
    for name in ("zvo_gc.dat", os.path.join("output", "zvo_out_001.dat")):
        actual = read_rows(os.path.join(actual_dir, name))
        reference = read_rows(os.path.join(reference_dir, name))
        if len(actual) != len(reference) or not reference:
            raise ComparisonFailure("{} {} rows differ".format(label, name))
        for row_a, row_r in zip(actual, reference):
            if len(row_a) != len(row_r):
                raise ComparisonFailure("{} {} columns".format(label, name))
            for a, r in zip(row_a, row_r):
                compare("{} {}".format(label, name), float(a), float(r),
                        relative_tolerance(float(r)))
    if mode == 1:
        actual = parse_green_outputs(actual_dir, mode)
        reference = parse_green_outputs(reference_dir, mode)
        for group in ("onebody", "twobody", "nbody", "anomalous"):
            compare_keyed("{} {}".format(label, group), actual[group],
                          dict(reference[group]))
    else:
        actual = parse_sr_dump(os.path.join(actual_dir, "sr.dat"))
        reference = parse_sr_dump(os.path.join(reference_dir, "sr.dat"))
        if len(actual) != len(reference) or not reference:
            raise ComparisonFailure("{} SR steps differ".format(label))
        for step_a, step_r in zip(actual, reference):
            if set(step_a["p"]) != set(step_r["p"]):
                raise ComparisonFailure("{} SR P sets differ".format(label))
            for index in step_r["p"]:
                for a, r in zip(step_a["p"][index], step_r["p"][index]):
                    compare("{} SR P{}".format(label, index), a, r,
                            relative_tolerance(r))


def smoke_run(rootdir, name, ap, mode, store, procs, nsplit, iterations):
    workdir = prepare_work(rootdir, name)
    write_fixture(workdir, ap=ap, mode=mode, samples=SMOKE_SAMPLES,
                  iterations=iterations, nstore=store, nsplit=nsplit,
                  optimize=True)
    env = {"MVMC_GC_STATE_DUMP": "state.dat"}
    if mode == 0:
        env["MVMC_GC_SR_DUMP"] = "sr.dat"
    run_binary(rootdir, workdir, procs, extra_env=env)
    return workdir


SMOKE_SAMPLES = 256


def smoke_case(rootdir, args):
    ap = args.boundary == "apbc"
    mode = args.mode
    iterations = 2 if mode == 0 else 1
    name = "smoke_{}_m{}_o{}_np{}_s{}".format(args.boundary, mode, args.store,
                                              args.np, args.nsplit)
    workdir = smoke_run(rootdir, name, ap, mode, args.store, args.np,
                        args.nsplit, iterations)
    records = sample_records(workdir, args.np)
    check_balanced(records)
    check_burn(workdir, args.np)
    check_finite_outputs(workdir, mode, iterations)
    if args.nsplit > 1:
        # One chain split over the projections: every rank holds it.
        by_rank = [[record[1:] for record in records if record[0] == rank]
                   for rank in range(args.np)]
        if any(rank_records != by_rank[0] for rank_records in by_rank):
            raise AssertionError("ranks of one chain disagree on samples")
        chain_records = [record for record in records if record[0] == 0]
    else:
        chain_records = records
    if args.np > 1 and args.nsplit > 1:
        reference = smoke_run(rootdir, name + "_serial", ap, mode, args.store,
                              1, 1, iterations)
        reference_records = sample_records(reference, 1)
        if [record[1:] for record in reference_records] != \
                [record[1:] for record in chain_records]:
            raise AssertionError("MPI occupations differ from serial")
        compare_output_files(name, workdir, reference, mode)
    parameters = [state_parameters()]
    if mode == 0:
        twin = smoke_run(rootdir, name + "_one_iteration", ap, mode,
                         args.store, args.np, args.nsplit, 1)
        twin_records = sample_records(twin, args.np)
        if args.nsplit > 1:
            twin_records = [r for r in twin_records if r[0] == 0]
        if iteration_masks(twin_records, SMOKE_SAMPLES, 0) != \
                iteration_masks(chain_records, SMOKE_SAMPLES, 0):
            raise AssertionError("first iteration is not reproducible")
        parameters.append(read_parameters(os.path.join(
            twin, "output", "zqp_opt.dat")))
    for iteration in range(iterations):
        masks = iteration_masks(chain_records, SMOKE_SAMPLES, iteration)
        expected = replay_values(masks, parameters[iteration], ap, mode)
        label = "{} iteration {}".format(name, iteration + 1)
        compare_physical(label, workdir, expected, row_index=iteration)
        if mode == 1:
            compare_greens(label, workdir, mode, expected)
        else:
            compare_sr(label, os.path.join(workdir, "sr.dat"), expected["sr"],
                       step_index=iteration)
    print("GC anti-parallel smoke {} passed ({} samples)".format(
        name, len(chain_records)))


# ---------------------------------------------------------------------------
# Comparators: actual, expected and tolerance must all be finite.
# ---------------------------------------------------------------------------
class ComparisonFailure(AssertionError):
    pass


def compare(label, actual, expected, tolerance):
    if not oracle.finite_close(actual, expected, tolerance):
        raise ComparisonFailure(
            "{}: actual={!r} expected={!r} tolerance={!r}".format(
                label, actual, expected, tolerance))


def relative_tolerance(expected, scale=2e-9):
    try:
        return scale * (1.0 + abs(complex(expected)))
    except (TypeError, ValueError):
        return float("nan")


def compare_keyed(label, actual_rows, expected, scale=2e-9):
    """actual_rows: list of (key, value).  Keys must match one to one."""
    seen = {}
    for key, value in actual_rows:
        if key in seen:
            raise ComparisonFailure("{}: duplicated row {}".format(label, key))
        seen[key] = value
    if not expected:
        raise ComparisonFailure("{}: empty expectation".format(label))
    if set(seen) != set(expected):
        raise ComparisonFailure("{}: rows differ missing={} extra={}".format(
            label, sorted(set(expected) - set(seen)),
            sorted(set(seen) - set(expected))))
    for key in sorted(expected):
        compare("{} {}".format(label, key), seen[key], expected[key],
                relative_tolerance(expected[key], scale))


# ---------------------------------------------------------------------------
# Output parsing
# ---------------------------------------------------------------------------
def complex_at(row, column):
    return complex(float(row[column]), float(row[column + 1]))


def numeric_row(path, row_index, columns):
    rows = read_rows(path)
    if row_index >= len(rows) or len(rows[row_index]) < columns:
        raise ComparisonFailure("{} lacks row {} with {} columns".format(
            path, row_index, columns))
    return rows[row_index]


def parse_physical_outputs(workdir, data_index=1, row_index=0):
    gc = numeric_row(os.path.join(workdir, "zvo_gc.dat"), row_index, 3)
    out = numeric_row(os.path.join(
        workdir, "output", "zvo_out_{:03d}.dat".format(data_index)),
        row_index, 6)
    return {
        "number": float(gc[0]), "number2": float(gc[1]),
        "variance_number": float(gc[2]),
        "energy": complex_at(out, 0), "energy2": float(out[2]),
        "relative_variance": float(out[3]),
        "sz": float(out[4]), "sz2": float(out[5]),
    }


def parse_green_outputs(workdir, mode, data_index=1):
    output = os.path.join(workdir, "output")
    tag = "{:03d}".format(data_index)
    onebody = []
    for row in read_rows(os.path.join(output, "zvo_cisajs_{}.dat".format(tag))):
        key = (spin_orbital(row[0], row[1]), spin_orbital(row[2], row[3]))
        onebody.append((key, complex_at(row, 4)))
    twobody = []
    for row in read_rows(os.path.join(
            output, "zvo_cisajscktalt_{}.dat".format(tag))):
        key = ((spin_orbital(row[0], row[1]), spin_orbital(row[2], row[3])),
               (spin_orbital(row[4], row[5]), spin_orbital(row[6], row[7])))
        twobody.append((key, complex_at(row, 8)))
    nbody = []
    for row in read_rows(os.path.join(output, "zvo_NBodyG_{}.dat".format(tag))):
        count = int(row[0])
        fields = row[1:1 + 4 * count]
        key = tuple((spin_orbital(fields[4 * k], fields[4 * k + 1]),
                     spin_orbital(fields[4 * k + 2], fields[4 * k + 3]))
                    for k in range(count))
        nbody.append((key, complex_at(row, 1 + 4 * count)))
    anomalous = []
    if mode == 1:
        for row in read_rows(os.path.join(
                output, "zvo_anomalousg_{}.dat".format(tag))):
            anomalous.append((tuple(int(value) for value in row[:5]),
                              complex_at(row, 5)))
    return {"onebody": onebody, "twobody": twobody, "nbody": nbody,
            "anomalous": anomalous}


def parse_sr_headers(path):
    headers = []
    with open(path) as stream:
        for line in stream:
            if line.startswith("STEP "):
                columns = line.split()
                headers.append(dict(zip(columns[::2], columns[1::2])))
    return headers


# ---------------------------------------------------------------------------
# Deterministic replay of saved samples with the independent oracle
# ---------------------------------------------------------------------------
NPROJ = 3
NPARA = NPROJ + NSITE * NSITE


def state_parameters(raw=None, gutz=oracle.GUTZ, jastrow=oracle.JASTROW):
    if raw is None:
        raw = oracle.F0
    return {"raw": np.array(raw, dtype=np.complex128), "gutz": gutz,
            "jastrow": jastrow}


def physical_matrix(params, ap):
    return params["raw"] * oracle.ap_signs(ap)


def operator_tables():
    onebody = dict(((a, b), one(a, b)) for a in range(NORBITAL)
                   for b in range(NORBITAL))
    twobody = dict((pairs, one(*pairs[0]) + one(*pairs[1]))
                   for pairs in TWO_BODY_PAIRS)
    nbody = {}
    for term in NBODY_G_TERMS:
        operators = ()
        for pair in term:
            operators += one(*pair)
        nbody[term] = operators
    anomalous = dict((key, oracle.anomalous_operator(key))
                     for key in ANOMALOUS_G_ROWS)
    return {"onebody": onebody, "twobody": twobody, "nbody": nbody,
            "anomalous": anomalous}


def replay_values(masks, params, ap, mode, greens=True):
    """Averages over the saved masks, duplicates included."""
    if not masks:
        raise AssertionError("no samples to replay")
    counts = {}
    for mask in masks:
        counts[mask] = counts.get(mask, 0) + 1
    total = float(len(masks))
    F = physical_matrix(params, ap)
    local = oracle.local_values(F, ap, list(counts), params["gutz"],
                                params["jastrow"], onebody_keys=[],
                                anomalous_list=[])

    def mean(function):
        return sum(counts[m] * function(m) for m in sorted(counts)) / total

    energy = mean(lambda m: local["energy"][m])
    energy2 = mean(lambda m: local["energy2"][m])
    number = mean(lambda m: local["number"][m])
    number2 = mean(lambda m: local["number2"][m])
    values = {
        "energy": energy, "energy2": energy2, "number": number,
        "number2": number2, "variance_number": number2 - number * number,
        "relative_variance": ((energy2 - energy * energy) /
                              (energy * energy)).real,
    }
    if greens:
        psi = oracle.wave_vector(F, params["gutz"], params["jastrow"])
        for group, table in operator_tables().items():
            if group == "anomalous" and mode != 1:
                continue
            expected = {}
            for key, operators in table.items():
                estimator = oracle.local_estimator(
                    oracle.operator_matrix(operators), psi, list(counts))
                expected[key] = mean(lambda m: estimator[m])
            values[group] = expected
    # SR rows: P(2q+c) for projection q<3 (real counts) and orbitals q>=3.
    sr = {}
    for q in range(NPARA):
        for component in (0, 1):
            if q < NPROJ:
                derivative = (lambda m, q=q, c=component:
                              float(local["projection"][m][q]) if c == 0
                              else 0.0)
            else:
                i, j = divmod(q - NPROJ, NSITE)
                derivative = (lambda m, i=i, j=j, c=component:
                              local["derivative"][m][i, j, bool(c)])
            mean_o = mean(lambda m: complex(derivative(m)))
            mean_oh = mean(lambda m: local["energy"][m] * derivative(m))
            gradient = 2.0 * (mean_oh.real - energy.real * mean_o.real)
            sr[2 * q + component] = (mean_o.real, mean_o.imag, mean_oh.real,
                                     mean_oh.imag, gradient)
    values["sr"] = sr
    return values


def compare_physical(label, workdir, expected, data_index=1, row_index=0,
                     scale=2e-9):
    actual = parse_physical_outputs(workdir, data_index, row_index)
    for key in ("number", "number2", "variance_number", "energy", "energy2",
                "relative_variance"):
        compare("{} {}".format(label, key), actual[key], expected[key],
                relative_tolerance(expected[key], scale))
    compare("{} Sz".format(label), actual["sz"], 0.0, 0.0)
    compare("{} Sz2".format(label), actual["sz2"], 0.0, 0.0)
    return actual


def compare_greens(label, workdir, mode, expected, data_index=1, scale=2e-9):
    actual = parse_green_outputs(workdir, mode, data_index)
    groups = ("onebody", "twobody", "nbody") + (("anomalous",)
                                                if mode == 1 else ())
    for group in groups:
        compare_keyed("{} {}".format(label, group), actual[group],
                      expected[group], scale)


def compare_sr(label, sr_path, expected, step_index=0, scale=2e-9):
    headers = parse_sr_headers(sr_path)
    if not headers:
        raise ComparisonFailure("{}: SR header is missing".format(label))
    for header in headers:
        if int(header["NPARA"]) != NPARA or int(header["NPROJ"]) != NPROJ:
            raise ComparisonFailure("{}: SR header {}".format(label, header))
    steps = parse_sr_dump(sr_path)
    if len(steps) != len(headers):
        raise ComparisonFailure("{}: SR step count".format(label))
    rows = steps[step_index]["p"]
    if set(rows) != set(range(2 * NPARA)):
        raise ComparisonFailure("{}: SR P rows {}".format(label, sorted(rows)))
    for index in range(2 * NPARA):
        actual = rows[index]
        target = expected[index]
        if len(actual) != 5:
            raise ComparisonFailure("{}: SR P {} has {} values".format(
                label, index, len(actual)))
        for column, (a, e) in enumerate(zip(actual, target)):
            compare("{} SR P{} column {}".format(label, index, column), a, e,
                    relative_tolerance(e, scale))
    return rows


def iteration_masks(records, nsample, iteration):
    """SAMPLE records of one iteration (records are in file order per rank)."""
    masks = []
    by_rank = {}
    for record in records:
        by_rank.setdefault(record[0], []).append(record)
    for rank in sorted(by_rank):
        rank_records = by_rank[rank]
        chunk = rank_records[iteration * nsample:(iteration + 1) * nsample]
        if len(chunk) != nsample:
            raise AssertionError("rank {} iteration {} has {} samples".format(
                rank, iteration, len(chunk)))
        if [record[2] for record in chunk] != list(range(nsample)):
            raise AssertionError("sample indices out of order")
        masks.extend(record[4] for record in chunk)
    return masks


def read_parameters(path):
    words = [float(value) for value in open(path).read().split()]
    values = words[6:]
    if len(values) != 3 * NPARA:
        raise AssertionError("{} holds {} parameter words".format(
            path, len(values)))
    params = [complex(values[3 * k], values[3 * k + 1]) for k in range(NPARA)]
    if params[2] != 0:
        raise AssertionError("fixed diagonal Jastrow changed: {}".format(
            params[2]))
    for k in range(NPROJ):
        if params[k].imag != 0:
            raise AssertionError("projection parameter became complex")
    raw = np.array(params[NPROJ:]).reshape(NSITE, NSITE)
    return state_parameters(raw, params[0].real, params[1].real)


def replay_case(rootdir, args):
    ap = args.boundary == "apbc"
    mode = args.mode
    name = "replay_{}_m{}".format(args.boundary, mode)
    workdir = prepare_work(rootdir, name)
    nsample = 256
    write_fixture(workdir, ap=ap, mode=mode, samples=nsample, iterations=1,
                  optimize=(mode == 0))
    env = {"MVMC_GC_STATE_DUMP": "state.dat"}
    if mode == 0:
        env["MVMC_GC_SR_DUMP"] = "sr.dat"
    run_binary(rootdir, workdir, 1, extra_env=env)
    records = sample_records(workdir)
    check_balanced(records)
    for record in records:
        if record[3] != bin(record[4]).count("1"):
            raise AssertionError("Ncur differs from the mask popcount")
    masks = iteration_masks(records, nsample, 0)
    expected = replay_values(masks, state_parameters(), ap, mode)
    compare_physical(name, workdir, expected)
    if mode == 1:
        compare_greens(name, workdir, mode, expected)
    if mode == 0:
        compare_sr(name, os.path.join(workdir, "sr.dat"), expected["sr"])
        # F00, F03, F22 real/imaginary parts are P6/7, P12/13, P26/27.
        for (i, j), packed in (((0, 0), 6), ((0, 3), 12), ((2, 2), 26)):
            if 2 * (NPROJ + NSITE * i + j) != packed:
                raise AssertionError("packed SR index mapping")
    print("GC anti-parallel replay {}: {} samples, {} distinct".format(
        name, len(masks), len(set(masks))))


# ---------------------------------------------------------------------------
# Statistics against exact values (Appendix B thresholds, fixed beforehand)
# ---------------------------------------------------------------------------
THRESHOLDS = {
    "pbc": {"energy": 0.027, "energy2": 0.11, "number": 0.063,
            "number2": 0.47, "anomalous": 0.014,
            6: 0.050, 7: 0.028, 12: 0.069, 13: 0.027, 26: 0.036, 27: 0.042},
    "apbc": {"energy": 0.047, "energy2": 0.12, "number": 0.086,
             "number2": 0.55, "anomalous": 0.012,
             6: 0.042, 7: 0.024, 12: 0.044, 13: 0.052, 26: 0.054, 27: 0.024},
}
VALIDATION_SEED = 91807
VALIDATION_SAMPLES = 60000
SR_COMPONENTS = {6: (0, 0, False), 7: (0, 0, True), 12: (0, 3, False),
                 13: (0, 3, True), 26: (2, 2, False), 27: (2, 2, True)}


def physical_case(rootdir, args):
    ap = args.boundary == "apbc"
    workdir = prepare_work(rootdir, "physical_" + args.boundary)
    write_fixture(workdir, ap=ap, mode=1, samples=VALIDATION_SAMPLES,
                  seed=VALIDATION_SEED, iterations=1)
    run_binary(rootdir, workdir, 1, extra_env={"MVMC_GC_STATE_DUMP": "state.dat"})
    records = sample_records(workdir)
    check_balanced(records)
    masks = iteration_masks(records, VALIDATION_SAMPLES, 0)
    if len(set(masks)) != 70:
        raise AssertionError("chain visited {} of 70 states".format(
            len(set(masks))))
    exact = oracle.exact(oracle.fixture_matrix(ap), ap)
    actual = parse_physical_outputs(workdir)
    limits = THRESHOLDS[args.boundary]
    for key in ("energy", "energy2", "number", "number2"):
        compare("physical {} {}".format(args.boundary, key), actual[key],
                exact[key], limits[key])
    anomalous = dict(parse_green_outputs(workdir, 1)["anomalous"])
    compare("physical {} anomalous".format(args.boundary),
            anomalous[oracle.ANOMALOUS_KEY],
            exact["anomalous"][oracle.ANOMALOUS_KEY], limits["anomalous"])
    print("GC anti-parallel physical {}: E={:.6f} (exact {:.6f}) N={:.5f}".format(
        args.boundary, actual["energy"].real, exact["energy"].real,
        actual["number"]))


def sr_case(rootdir, args):
    ap = args.boundary == "apbc"
    exact = oracle.exact(oracle.fixture_matrix(ap), ap)
    limits = THRESHOLDS[args.boundary]
    runs = {}
    for store in (0, 1):
        workdir = prepare_work(rootdir, "sr_{}_o{}".format(args.boundary,
                                                            store))
        write_fixture(workdir, ap=ap, mode=0, samples=VALIDATION_SAMPLES,
                      seed=VALIDATION_SEED, nstore=store, iterations=1,
                      optimize=True, measurements=False)
        run_binary(rootdir, workdir, 1,
                   extra_env={"MVMC_GC_STATE_DUMP": "state.dat",
                              "MVMC_GC_SR_DUMP": "sr.dat"})
        masks = iteration_masks(sample_records(workdir), VALIDATION_SAMPLES, 0)
        rows = parse_sr_dump(os.path.join(workdir, "sr.dat"))[0]["p"]
        runs[store] = (masks, rows, workdir)
        for packed, key in SR_COMPONENTS.items():
            compare("sr {} NStoreO={} P{}".format(args.boundary, store, packed),
                    rows[packed][4], exact["gradient"][key], limits[packed])
    if runs[0][0] != runs[1][0]:
        raise AssertionError("NStoreO changed the sampled configurations")
    if set(runs[0][1]) != set(range(2 * NPARA)) or \
            set(runs[1][1]) != set(range(2 * NPARA)):
        raise AssertionError("SR P rows are incomplete")
    for index in range(2 * NPARA):
        for column in range(5):
            reference = runs[0][1][index][column]
            compare("sr {} NStoreO P{} column {}".format(
                args.boundary, index, column), runs[1][1][index][column],
                reference, relative_tolerance(reference, 5e-10))
    # Deterministic replay of the same samples.
    expected = replay_values(runs[0][0], state_parameters(), ap, 0,
                             greens=False)
    compare_sr("sr {}".format(args.boundary),
               os.path.join(runs[0][2], "sr.dat"), expected["sr"])
    print("GC anti-parallel SR {} passed".format(args.boundary))


def optimization_case(rootdir, args):
    ap = args.boundary == "apbc"
    results = {}
    for iterations in (1, 2):
        workdir = prepare_work(rootdir, "optimization_{}_i{}".format(
            args.boundary, iterations))
        write_fixture(workdir, ap=ap, mode=0, samples=1024, iterations=iterations,
                      optimize=True, measurements=False)
        run_binary(rootdir, workdir, 1,
                   extra_env={"MVMC_GC_STATE_DUMP": "state.dat"})
        results[iterations] = workdir
    records_1 = sample_records(results[1])
    records_2 = sample_records(results[2])
    check_balanced(records_2)
    first = iteration_masks(records_2, 1024, 0)
    if first != iteration_masks(records_1, 1024, 0):
        raise AssertionError("first iteration differs between the two runs")
    updated = read_parameters(os.path.join(results[1], "output",
                                           "zqp_opt.dat"))
    # No uniform rescale: each orbital moves by a small SR step only.
    change = np.max(abs(updated["raw"] - oracle.F0))
    if not (0.0 < change < 0.05):
        raise AssertionError("orbital update size {}".format(change))
    ratio = abs(updated["raw"]) / abs(oracle.F0)
    if np.max(abs(ratio - 1.0)) > 0.1:
        raise AssertionError("orbital magnitudes were rescaled: {}".format(
            ratio))
    for iteration, params in ((0, state_parameters()), (1, updated)):
        masks = iteration_masks(records_2, 1024, iteration)
        expected = replay_values(masks, params, ap, 0, greens=False)
        compare_physical("optimization {} iteration {}".format(
            args.boundary, iteration + 1), results[2], expected,
            row_index=iteration)
    print("GC anti-parallel optimization {}: max orbital step {:.3g}".format(
        args.boundary, change))


# ---------------------------------------------------------------------------
# Same state as OrbitalGeneral F/2, and mutations that must break agreement
# ---------------------------------------------------------------------------
def expect_replay_failure(label, function):
    try:
        function()
    except ComparisonFailure as error:
        return str(error)
    raise AssertionError("{}: mutated input still matches the oracle".format(
        label))


def replay_against_oracle(workdir, ap, mode, nsample, params=None):
    records = sample_records(workdir)
    check_balanced(records)
    masks = iteration_masks(records, nsample, 0)
    expected = replay_values(masks, params or state_parameters(), ap, mode)
    label = os.path.basename(workdir)
    compare_physical(label, workdir, expected)
    compare_greens(label, workdir, mode, expected)
    return expected


def general_run(rootdir, ap, name, scale=0.5):
    workdir = prepare_work(rootdir, name)
    write_fixture(workdir, ap=ap, mode=1, samples=256, representation="general")
    if scale != 0.5:
        values = [value * (scale / 0.5) for value in
                  orbital_parameters("general")]
        write_initial(workdir, "general", parameters=values)
    run_binary(rootdir, workdir, 1, extra_env={"MVMC_GC_STATE_DUMP": "state.dat"})
    return workdir


def flip_anomalous(workdir, ap):
    def change(lines):
        out = lines[:5]
        for line in lines[5:]:
            fields = line.split()
            if not fields:
                continue
            fields[5] = repr(-float(fields[5]))
            fields[6] = repr(-float(fields[6]))
            out.append(" ".join(fields) + "\n")
        return out
    edit_file(workdir, "anomalousterm.def", change)


HAMILTONIAN_MUTATIONS = {
    # PairHop is one input row with its reverse; adding the reverse doubles it.
    "pairhop_duplicate_reverse": lambda workdir, ap: append_rows(
        workdir, "pairhop.def", ["1 0 0.13\n"], 1),
    # Drop the conjugate partner of the first hopping row.
    "transfer_missing_hc": lambda workdir, ap: edit_file(
        workdir, "trans.def", lambda lines: (
            lines[:1] + ["{} {}\n".format(lines[1].split()[0],
                                          int(lines[1].split()[1]) - 1)] +
            lines[2:5] + [line for line in lines[5:]
                          if line.split() != "1 0 0 0 0.37 -0.19".split()])),
    # Exchange already contains both directions in one row.
    "exchange_double_hc": lambda workdir, ap: append_rows(
        workdir, "exchange.def", ["0 1 -0.11\n"], 1),
}


def run_hamiltonian_mutation(rootdir, mutation, ap):
    boundary = "apbc" if ap else "pbc"
    workdir = prepare_work(rootdir, "mutation_{}_{}".format(mutation,
                                                            boundary))
    write_fixture(workdir, ap=ap, mode=1, samples=256)
    HAMILTONIAN_MUTATIONS[mutation](workdir, ap)
    # 1) The independent read-back must already see the defect.
    difference, unused = check_hamiltonian(workdir, ap)
    if not difference > 1e-3:
        raise AssertionError("{}: read-back H unchanged ({})".format(
            mutation, difference))
    # 2) Bypass that pre-check only here: the production run must reach the
    #    calculation and its energy must disagree with the correct model.
    output = run_binary(rootdir, workdir, 1,
                        extra_env={"MVMC_GC_STATE_DUMP": "state.dat"})
    if "Start: Main calculation" not in output and \
            "Main calculation" not in output:
        raise AssertionError("{}: calculation not reached".format(mutation))
    records = sample_records(workdir)
    check_balanced(records)
    masks = iteration_masks(records, 256, 0)
    expected = replay_values(masks, state_parameters(), ap, 1, greens=False)
    actual = parse_physical_outputs(workdir)
    delta = abs(actual["energy"] - expected["energy"])
    message = expect_replay_failure(mutation, lambda: compare(
        "{} energy".format(mutation), actual["energy"], expected["energy"],
        relative_tolerance(expected["energy"])))
    return difference, delta, message


def mutation_case(rootdir, args):
    ap = args.boundary == "apbc"
    mutation = args.mutation
    if mutation == "general":
        workdir = general_run(rootdir, ap, "general_{}".format(args.boundary))
        replay_against_oracle(workdir, ap, 1, 256)
        print("GC anti-parallel General F/2 {} matches the oracle".format(
            args.boundary))
        return
    if mutation == "general_full_f":
        workdir = general_run(rootdir, ap, "general_full_f_{}".format(
            args.boundary), scale=1.0)
        message = expect_replay_failure(mutation, lambda: replay_against_oracle(
            workdir, ap, 1, 256))
    elif mutation in ("ap_sign_flip", "anomalous_sign"):
        workdir = prepare_work(rootdir, "{}_{}".format(mutation,
                                                       args.boundary))
        write_fixture(workdir, ap=ap, mode=1, samples=256)
        if mutation == "ap_sign_flip":
            if not ap:
                raise SystemExit("ap_sign_flip applies to APBC")
            replace_line(workdir, ORBITAL, PAIR0 + 3, pair_line(ap, 0, 3, 3, 1))
        else:
            flip_anomalous(workdir, ap)
        run_binary(rootdir, workdir, 1,
                   extra_env={"MVMC_GC_STATE_DUMP": "state.dat"})
        message = expect_replay_failure(mutation, lambda: replay_against_oracle(
            workdir, ap, 1, 256))
    elif mutation in HAMILTONIAN_MUTATIONS:
        difference, delta, message = run_hamiltonian_mutation(
            rootdir, mutation, ap)
        print("GC anti-parallel {} {}: read-back |dH|={:.4g}, "
              "production |dE|={:.4g}".format(mutation, args.boundary,
                                              difference, delta))
    else:
        raise SystemExit("unknown mutation {}".format(mutation))
    print("GC anti-parallel mutation {} {} detected: {}".format(
        mutation, args.boundary, message.splitlines()[0][:160]))


# ---------------------------------------------------------------------------
# The comparators themselves must reject nonfinite or missing data.
# ---------------------------------------------------------------------------
def must_fail(label, function):
    try:
        function()
    except ComparisonFailure:
        return
    raise AssertionError("comparator accepted {}".format(label))


def write_physical_files(workdir, values):
    write(os.path.join(workdir, "zvo_gc.dat"),
          " {} {} {}\n".format(*[fmt(values[key]) for key in
                                 ("number", "number2", "variance_number")]))
    os.makedirs(os.path.join(workdir, "output"), exist_ok=True)
    write(os.path.join(workdir, "output", "zvo_out_001.dat"),
          " {} {}  {} {} {} {}\n".format(
              fmt(values["energy"].real), fmt(values["energy"].imag),
              fmt(values["energy2"]), fmt(values["relative_variance"]),
              fmt(0.0), fmt(0.0)))


def write_green_files(workdir, values):
    output = os.path.join(workdir, "output")
    rows = []
    for (a, b), value in sorted(values["onebody"].items()):
        rows.append("{} {} {} {} {} {}\n".format(
            a % NSITE, a // NSITE, b % NSITE, b // NSITE, fmt(value.real),
            fmt(value.imag)))
    write(os.path.join(output, "zvo_cisajs_001.dat"), "".join(rows))
    rows = []
    for pairs, value in sorted(values["twobody"].items()):
        rows.append("{} {} {}\n".format(factors_text(pairs), fmt(value.real),
                                        fmt(value.imag)))
    write(os.path.join(output, "zvo_cisajscktalt_001.dat"), "".join(rows))
    rows = []
    for term, value in sorted(values["nbody"].items()):
        rows.append("{} {} {} {}\n".format(len(term), factors_text(term),
                                           fmt(value.real), fmt(value.imag)))
    write(os.path.join(output, "zvo_NBodyG_001.dat"), "".join(rows))
    rows = []
    for key, value in sorted(values["anomalous"].items()):
        rows.append("{} {} {} {} {} {} {}\n".format(
            *(key + (fmt(value.real), fmt(value.imag)))))
    write(os.path.join(output, "zvo_anomalousg_001.dat"), "".join(rows))


def write_sr_file(path, sr, header=True, drop=None, replace=None):
    lines = []
    if header:
        lines.append("STEP 0 NPARA {} SROPTSIZE {} NPROJ {} STOREO 1\n".format(
            NPARA, NPARA + 1, NPROJ))
    for index in range(2 * NPARA):
        if index == drop:
            continue
        values = list(sr[index])
        if replace is not None and replace[0] == index:
            values[replace[1]] = replace[2]
        lines.append("P {} {}\n".format(index, " ".join(
            "{:.17e}".format(value) for value in values)))
    write(path, "".join(lines))


def comparator_guard_case(rootdir, args):
    workdir = prepare_work(rootdir, "comparator_guard")
    masks = [0, 17, 51, 255, 51]
    expected = replay_values(masks, state_parameters(), False, 1)
    write_physical_files(workdir, expected)
    write_green_files(workdir, expected)
    compare_physical("guard control", workdir, expected)
    compare_greens("guard control", workdir, 1, expected)
    bad_values = (float("nan"), float("inf"), float("-inf"))
    # Nonfinite actual values in every physical column.
    for key in ("number", "number2", "variance_number", "energy", "energy2",
                "relative_variance"):
        for bad in bad_values:
            mutated = dict(expected)
            mutated[key] = complex(bad, 0) if key == "energy" else bad
            write_physical_files(workdir, mutated)
            must_fail("actual {}={}".format(key, bad), lambda: compare_physical(
                "guard", workdir, expected))
            must_fail("expected {}={}".format(key, bad), lambda: compare_physical(
                "guard", workdir, mutated))
    write_physical_files(workdir, expected)
    # Nonfinite Green and anomalous values, duplicated and missing rows.
    for group in ("onebody", "twobody", "nbody", "anomalous"):
        key = sorted(expected[group])[0]
        for bad in bad_values:
            mutated = dict(expected)
            mutated[group] = dict(expected[group])
            mutated[group][key] = complex(bad, 0.0)
            write_green_files(workdir, mutated)
            must_fail("{} actual {}".format(group, bad), lambda: compare_greens(
                "guard", workdir, 1, expected))
            must_fail("{} expected {}".format(group, bad),
                      lambda: compare_greens("guard", workdir, 1, mutated))
        rows = [(k, v) for k, v in expected[group].items()]
        must_fail("{} duplicate".format(group), lambda: compare_keyed(
            "guard", rows + rows[:1], expected[group]))
        must_fail("{} missing".format(group), lambda: compare_keyed(
            "guard", rows[1:], expected[group]))
        must_fail("{} empty".format(group), lambda: compare_keyed(
            "guard", [], {}))
    write_green_files(workdir, expected)
    write(os.path.join(workdir, "output", "zvo_out_001.dat"), "")
    must_fail("removed numeric output", lambda: compare_physical(
        "guard", workdir, expected))
    # Tolerances must be finite and nonnegative.
    for tolerance in (float("nan"), float("inf"), -1.0):
        must_fail("tolerance {}".format(tolerance), lambda: compare(
            "guard", 1.0, 1.0, tolerance))
    # SR dump: control, missing header, missing P, nonfinite entries.
    sr_path = os.path.join(workdir, "sr.dat")
    write_sr_file(sr_path, expected["sr"])
    compare_sr("guard control", sr_path, expected["sr"])
    write_sr_file(sr_path, expected["sr"], header=False)
    must_fail("SR header missing", lambda: compare_sr("guard", sr_path,
                                                      expected["sr"]))
    write_sr_file(sr_path, expected["sr"], drop=7)
    must_fail("SR P missing", lambda: compare_sr("guard", sr_path,
                                                 expected["sr"]))
    for bad in bad_values:
        for column in range(5):
            write_sr_file(sr_path, expected["sr"], replace=(12, column, bad))
            must_fail("SR value {} column {}".format(bad, column),
                      lambda: compare_sr("guard", sr_path, expected["sr"]))
        corrupted = dict(expected["sr"])
        corrupted[6] = (bad,) + tuple(corrupted[6][1:])
        write_sr_file(sr_path, expected["sr"])
        must_fail("SR expected {}".format(bad), lambda: compare_sr(
            "guard", sr_path, corrupted))
    write(sr_path, "")
    must_fail("SR empty", lambda: compare_sr("guard", sr_path,
                                             expected["sr"]))
    print("GC anti-parallel comparator guard passed")


SAMPLE_SAMPLES = 10000


def sample_case(rootdir, args):
    """The published sample equals the PBC fixture and replays exactly."""
    source = args.sample_dir
    if not source or not os.path.isdir(source):
        raise SystemExit("--sample-dir is required")
    reference = prepare_work(rootdir, "sample_reference")
    write_fixture(reference, ap=False, mode=1, samples=SAMPLE_SAMPLES)
    published = sorted(name for name in os.listdir(source)
                       if name != "README.md")
    if published != sorted(os.listdir(reference)):
        raise AssertionError("sample files {} differ from {}".format(
            published, sorted(os.listdir(reference))))
    for name in published:
        with open(os.path.join(source, name), "rb") as stream:
            actual = stream.read()
        with open(os.path.join(reference, name), "rb") as stream:
            expected = stream.read()
        if actual != expected:
            raise AssertionError("sample file {} differs from the fixture"
                                 .format(name))
    workdir = prepare_work(rootdir, "sample_run")
    for name in published:
        shutil.copy(os.path.join(source, name), os.path.join(workdir, name))
    run_binary(rootdir, workdir, 1, extra_env={"MVMC_GC_STATE_DUMP": "state.dat"})
    replay_against_oracle(workdir, False, 1, SAMPLE_SAMPLES)
    print("GC anti-parallel sample reproduces the fixture and replays")


def inventory_case(rootdir, args):
    """The CMake registration must list exactly the runner's mutations."""
    registered = sorted(name for name in (args.mutation or "").split(",")
                        if name)
    known = sorted(INPUT_MUTATIONS)
    if registered != known:
        raise AssertionError(
            "CMake/runner input mutation lists differ: missing={} extra={}"
            .format(sorted(set(known) - set(registered)),
                    sorted(set(registered) - set(known))))
    print("GC anti-parallel input inventory: {} mutations".format(len(known)))


CASES = {
    "sample": sample_case,
    "mutation": mutation_case,
    "comparator_guard": comparator_guard_case,
    "replay": replay_case,
    "physical": physical_case,
    "sr": sr_case,
    "optimization": optimization_case,
    "inventory": inventory_case,
    "hamiltonian": hamiltonian_case,
    "input": input_case,
    "audit": audit_case,
    "smoke": smoke_case,
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case")
    parser.add_argument("--boundary", choices=("pbc", "apbc"), default="pbc")
    parser.add_argument("--np", type=int, default=1)
    parser.add_argument("--nsplit", type=int, default=1)
    parser.add_argument("--mode", type=int, default=1)
    parser.add_argument("--store", type=int, default=1)
    parser.add_argument("--mutation", default=None)
    parser.add_argument("--sample-dir", default=None)
    args = parser.parse_args()
    rootdir = os.getcwd()
    cases = globals().get("CASES", {})
    if args.case not in cases:
        raise SystemExit("unknown case {}".format(args.case))
    cases[args.case](rootdir, args)


if __name__ == "__main__":
    main()
