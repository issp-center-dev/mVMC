"""Grand-canonical S=0 spin projection (anti-parallel pairs, one translation).

The four-site anti-parallel fixture runs with NSPGaussLeg=G>1, NSPStot=0 and
abs(NMPTrans)=1.  Production output is replayed against the finite-G Fock
oracle (grandcanonical_spin_oracle) evaluated on the binary64 quadrature
coefficients dumped by the Testing-only MVMC_GC_SPIN_QP_DUMP hook, and the
quadrature is checked to converge to the exact S=0 projection.
"""
from __future__ import print_function

import argparse
import concurrent.futures
import hashlib
import json
import math
import os
import re
import statistics
import subprocess
import sys
from fractions import Fraction

import grandcanonical_ap_oracle as ap
import grandcanonical_spin_oracle as spin
import runtest_gc_mp as mp
from grandcanonical_exact_oracle import apply_ops, nbody, one_body
from runtest_gc import (
    fmt,
    mpi_command,
    prepare_work,
    run_binary,
    state_dump_records,
    write,
)
from runtest_gc_antiparallel import ComparisonFailure, compare, parse_sr_dump
from runtest_gc_ap import def_text, write_hamiltonian


NSITE = mp.NSITE
NORBITAL = 2 * NSITE
NPROJ = mp.NPROJ
REPLAY_SAMPLES = 256
REPLAY_SEED = 4321
REPLAY_TOLERANCE = 2.0e-9
FORBIDDEN_TOLERANCE = 2.0e-12
NEGATIVE_TIMEOUT = 50
DUMP_NAME = "spinqp.dat"

# S^2 = Sz^2 + (S+S- + S-S+)/2 from the two-body Green functions: S+_i S-_j
# and S-_i S+_j for every site pair.  The density/exchange rows of the
# momentum fixture and two S^z-changing products complete the inventory.
S2_TWO_BODY_KEYS = tuple(
    ((i, i + NSITE), (j + NSITE, j)) for i in range(NSITE)
    for j in range(NSITE)) + tuple(
    ((i + NSITE, i), (j, j + NSITE)) for i in range(NSITE)
    for j in range(NSITE))
TWO_BODY_KEYS = S2_TWO_BODY_KEYS + tuple(mp.TWO_BODY_KEYS) + (
    ((0, 4), (1, 5)), ((0, 4), (5, 5)))
# Spin exchange times a density (total S^z kept) and an S^z-changing one.
NBODY_KEYS = (((0, 4), (5, 1), (2, 2)), ((0, 4), (1, 1), (2, 2)))
# (kind, first, second): up-down pairs in both orders and same-spin pairs.
ANOMALOUS_KEYS = ((1, 3, 4), (0, 4, 3), (1, 4, 3), (0, 3, 4),
                  (1, 0, 1), (0, 1, 0), (1, 5, 6), (0, 6, 5))


def anomalous_operators(kind, first, second):
    letter = "c" if kind == 1 else "a"
    return ((letter, first), (letter, second))


def twice_sz_change(operators):
    """2 Delta S^z of an operator string (creators +, annihilators -)."""
    total = 0
    for kind, orbital in operators:
        twice = 1 if orbital < NSITE else -1
        total += twice if kind == "c" else -twice
    return total


def fixture_pair_amplitudes():
    """Pair amplitudes with S=1 and S=2 content in every allowed sector.

    The antisymmetric part has rank > 2, so N=4 carries S=2 weight.  A rank-2
    antisymmetric part (as in the oracle self-test F) has no S=2 weight; then
    the symmetric quadrature removes S=1 exactly and every G >= 2 gives the
    same state, which would make the G dependence of the replay vacuous.
    """
    F = [[complex((7 + 3 * i - 2 * j) / 16.0, (1 + i + 2 * j) / 32.0)
          for j in range(NSITE)] for i in range(NSITE)]
    F[0][1] += 0.7
    F[2][0] -= 0.3j
    # Avoid the accidental N=6 node at mask 238 in the original probe.
    # A full-support fixture lets the SR rows be checked independently
    # against derivatives of the full Rayleigh quotient.
    F[3][3] += 1.0 / 64.0
    return F


class SpinFixture(object):
    def __init__(self, boundary, gauss, kind="dense"):
        if boundary not in ("pbc", "apbc"):
            raise ValueError("unknown boundary {}".format(boundary))
        if kind != "dense":
            raise ValueError("unknown fixture kind {}".format(kind))
        if gauss < 1:
            raise ValueError("positive NSPGaussLeg required")
        self.boundary = boundary
        self.gauss = gauss
        self.kind = kind
        self.mp = mp.Fixture("anti", boundary, 0)
        self.raw_F = fixture_pair_amplitudes()

    def label(self):
        return "gc_sp_{}_g{}".format(self.boundary, self.gauss)

    def nmptrans(self):
        return -1 if self.boundary == "apbc" else 1

    def parameters(self):
        return [self.raw_F[i][j] for i in range(NSITE) for j in range(NSITE)]

    def log_charge(self, mask):
        return self.mp.log_projection(mask)

    def terms(self):
        return self.mp.terms()


# ---------------------------------------------------------------------------
# Input writers
# ---------------------------------------------------------------------------
def modpara_text(mode_cal, samples, seed, nstore, nsplit, iterations,
                 gauss_token, stot_token, nmptrans, freeze, warmup, order,
                 grand=1, two_sz=0, init_nelec=4):
    data_qty = iterations if mode_cal == 1 else 1
    opt_steps = iterations if mode_cal == 0 else 1
    step = 0.0 if freeze else 0.002
    gc_lines = ""
    if grand is not None:
        gc_lines = "NGrandCanonical {}\n".format(grand)
        if grand == 1:
            gc_lines += "NGCInitNelec   {}\n".format(init_nelec)
    early = gc_lines if order == "early" else ""
    late = gc_lines if order == "late" else ""
    return """--------------------
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
{early}Nsite          {nsite}
Ne             2
Ncond          4
2Sz            {two_sz}
NSPGaussLeg    {gauss}
NSPStot        {stot}
NMPTrans       {nmptrans}
NSROptItrStep  {opt_steps}
NSROptItrSmp   1
DSROptRedCut   0.000000000001
DSROptStaDel   0.01
DSROptStepDt   {step}
NVMCWarmUp     {warmup}
NVMCInterval   4
NVMCSample     {samples}
NExUpdatePath  0
RndSeed        {seed}
NSplitSize     {nsplit}
NStore         {nstore}
NSRCG          0
{late}""".format(mode=mode_cal, data_qty=data_qty, early=early, late=late,
                 nsite=NSITE, two_sz=two_sz, gauss=gauss_token,
                 stot=stot_token, nmptrans=nmptrans, opt_steps=opt_steps,
                 step=fmt(step), warmup=warmup, samples=samples, seed=seed,
                 nsplit=nsplit, nstore=nstore)


def write_identity_translation(workdir, count=1):
    """TransSym with the identity (weight 1), or identity+shift for count 2."""
    if count == 1:
        rows = ["0 {} {}\n".format(fmt(1.0), fmt(0.0))]
        rows.extend("0 {} {} 1\n".format(site, site) for site in range(NSITE))
    else:
        rows = ["{} {} {}\n".format(r, fmt(0.5), fmt(0.0)) for r in range(2)]
        rows.extend("0 {} {} 1\n".format(site, site) for site in range(NSITE))
        rows.extend("1 {} {} 1\n".format(site, (site + 1) % NSITE)
                    for site in range(NSITE))
    write(os.path.join(workdir, "qptransidx.def"),
          def_text("NQPTrans", count, rows))


def factors_text(pairs):
    fields = []
    for out_orbital, in_orbital in pairs:
        fields.extend(mp.fused_fields(out_orbital))
        fields.extend(mp.fused_fields(in_orbital))
    return " ".join(str(value) for value in fields)


def write_measurements(workdir, mode_cal, anomalous_hamiltonian=True):
    if not anomalous_hamiltonian:
        # Canonical compatibility input: spin-diagonal one-body Green only.
        one_rows = ["{}\n".format(factors_text(((first, second),)))
                    for first in range(NORBITAL) for second in range(NORBITAL)
                    if (first < NSITE) == (second < NSITE)]
        write(os.path.join(workdir, "greenone.def"),
              def_text("NCisAjs", len(one_rows), one_rows))
        return ["        OneBodyG  greenone.def"]
    one_rows = ["{}\n".format(factors_text(((first, second),)))
                for first in range(NORBITAL) for second in range(NORBITAL)]
    write(os.path.join(workdir, "greenone.def"),
          def_text("NCisAjs", len(one_rows), one_rows))
    two_rows = ["{}\n".format(factors_text(pairs)) for pairs in TWO_BODY_KEYS]
    write(os.path.join(workdir, "greentwo.def"),
          def_text("NCisAjsCktAltDC", len(two_rows), two_rows))
    nbody_rows = ["{} {}\n".format(len(term), factors_text(term))
                  for term in NBODY_KEYS]
    write(os.path.join(workdir, "nbodyg.def"),
          def_text("NNBodyG", len(nbody_rows), nbody_rows))
    names = ["        OneBodyG  greenone.def",
             "        TwoBodyG  greentwo.def",
             "          NBodyG  nbodyg.def"]
    if anomalous_hamiltonian:
        create = (1,) + mp.fused_fields(3) + mp.fused_fields(4) + (mp.DELTA,)
        remove = ((0,) + mp.fused_fields(4) + mp.fused_fields(3) +
                  (mp.DELTA.conjugate(),))
        rows = ["{} {} {} {} {} {} {}\n".format(
            row[0], row[1], row[2], row[3], row[4],
            fmt(row[5].real), fmt(row[5].imag)) for row in (create, remove)]
        write(os.path.join(workdir, "anomalousterm.def"),
              def_text("NAnomalousTerm", len(rows), rows))
        names.append("   AnomalousTerm  anomalousterm.def")
    if mode_cal == 1 and anomalous_hamiltonian:
        rows = ["{} {} {} {} {}\n".format(
            kind, *(mp.fused_fields(first) + mp.fused_fields(second)))
            for kind, first, second in ANOMALOUS_KEYS]
        write(os.path.join(workdir, "anomalousg.def"),
              def_text("NAnomalousG", len(rows), rows))
        names.append("      AnomalousG  anomalousg.def")
    return names


def write_spin_jastrow(workdir):
    rows = ["{} {} 0\n".format(first, second) for first in range(NSITE)
            for second in range(NSITE) if first != second]
    rows.append("0 1\n")
    write(os.path.join(workdir, "spinjastrow.def"),
          def_text("NSpinJastrowIdx", 1, rows, complex_type=0))
    return "     SpinJastrow  spinjastrow.def"


def write_initial(workdir, parameters, spin_jastrow=False):
    """initial.def: Proj = (Gutzwiller, Jastrow[, SpinJastrow]) + Slater."""
    projection = (mp.GUTZWILLER,) + tuple(mp.JASTROW)
    if spin_jastrow:
        projection += (0.0,)
    values = [0.0] * 6
    for value in projection:
        values.extend((value, 0.0, 0.0))
    for value in parameters:
        values.extend((value.real, value.imag, 0.0))
    write(os.path.join(workdir, "initial.def"),
          " ".join(fmt(value) for value in values) + "\n")


def write_fixture(workdir, fixture, calmode, store, nsplit, samples, seed,
                  iterations, freeze, order="early", warmup=20,
                  gauss_token=None, stot_token="0", nmptrans=None,
                  optimize=None, general=False, grand=1, translations=1,
                  spin_jastrow=False, opttrans=False, momentum=False):
    """Write a complete expert-mode input; returns extra command options."""
    if optimize is None:
        optimize = calmode == 0
    if gauss_token is None:
        gauss_token = str(fixture.gauss)
    if nmptrans is None:
        nmptrans = fixture.nmptrans()
    two_sz = -1 if general else 0
    write(os.path.join(workdir, "modpara.def"), modpara_text(
        calmode, samples, seed, store, nsplit, iterations, gauss_token,
        stot_token, nmptrans, freeze, warmup, order, grand=grand,
        two_sz=two_sz))
    mp.write_projection(workdir, fixture.mp, optimize)
    if momentum:
        mp.write_translations(workdir, fixture.mp.lattice, 0)
    else:
        write_identity_translation(workdir, translations)
    if general:
        general_fixture = mp.Fixture("general", fixture.boundary, 0)
        mp.write_general_orbital(workdir, general_fixture, optimize)
        orbital_name = "  OrbitalGeneral  orbitalidxgen.def"
        parameters = [0.1 + 0.05j] * len(general_fixture.general_pairs)
        model = general_fixture.model(parameters)
    else:
        mp.write_antiparallel_orbital(workdir, fixture.mp, optimize)
        orbital_name = "OrbitalAntiParallel  orbitalidx.def"
        parameters = fixture.parameters()
        model = fixture.mp.model(())
    write_hamiltonian(workdir, model)
    names = [
        "         ModPara  modpara.def",
        "         LocSpin  locspn.def",
        "      Gutzwiller  gutzwilleridx.def",
        "         Jastrow  jastrowidx.def",
        orbital_name,
        "        TransSym  qptransidx.def",
        "           Trans  trans.def",
        "    CoulombIntra  coulombintra.def",
        "    CoulombInter  coulombinter.def",
    ]
    names.extend(write_measurements(workdir, calmode,
                                    anomalous_hamiltonian=(grand == 1)))
    if spin_jastrow:
        names.append(write_spin_jastrow(workdir))
    options = []
    if opttrans:
        write(os.path.join(workdir, "opttrans.def"),
              def_text("NQPOptTrans", 2, []))
        names.append("        OptTrans  opttrans.def")
        options.append("-o")
    write(os.path.join(workdir, "namelist.def"), "\n".join(names) + "\n")
    write_initial(workdir, parameters, spin_jastrow)
    return options


# ---------------------------------------------------------------------------
# Independent readback of the written input
# ---------------------------------------------------------------------------
def read_pair_matrix(workdir):
    """F_ij = sign_ij * f_k from orbitalidx.def and initial.def."""
    path = os.path.join(workdir, "orbitalidx.def")
    with open(path) as stream:
        lines = [line.split() for line in stream if line.strip()]
    rows = [line for line in lines[5:] if len(line) in (3, 4)]
    entries = [line for line in rows[:NSITE * NSITE]]
    if len(entries) != NSITE * NSITE:
        raise AssertionError("orbitalidx.def does not list every pair")
    unused, parameters = mp.read_initial_parameters(
        os.path.join(workdir, "initial.def"), NSITE * NSITE)
    matrix = [[None] * NSITE for unused_row in range(NSITE)]
    for line in entries:
        i, j, index = int(line[0]), int(line[1]), int(line[2])
        sign = int(line[3]) if len(line) == 4 else 1
        if sign not in (1, -1) or matrix[i][j] is not None:
            raise AssertionError("bad orbital row {}".format(line))
        matrix[i][j] = sign * parameters[index]
    return matrix


def read_dump(workdir, gauss):
    dump = spin.read_qp_dump(os.path.join(workdir, DUMP_NAME))
    if dump["gauss"] != gauss or dump["count"] != gauss:
        raise AssertionError("dump sizes G={} Q={} for NSPGaussLeg={}".format(
            dump["gauss"], dump["count"], gauss))
    return dump["nodes"]


# ---------------------------------------------------------------------------
# Expected values
# ---------------------------------------------------------------------------
def local_of(wave, mask, operators):
    if wave[mask] == 0:
        raise AssertionError("local value on a zero amplitude")
    moved = apply_ops(mask, operators)
    return (0.0j if moved is None else
            moved[1] * (wave[moved[0]] / wave[mask]).conjugate())


def mean_over(masks, values):
    return sum(values[mask] for mask in masks) / float(len(masks))


def popcount(value):
    return bin(value).count("1")


def local_tables(wave, terms, unique, calmode):
    energy = {}
    for mask in unique:
        energy[mask] = sum(coefficient * local_of(wave, mask, operators)
                           for unused, coefficient, operators in terms)
    tables = {"energy": energy}
    tables["onebody"] = dict(
        ((first, second), dict((mask, local_of(wave, mask,
                                               one_body(first, second)))
                               for mask in unique))
        for first in range(NORBITAL) for second in range(NORBITAL))
    tables["twobody"] = dict(
        (key, dict((mask, local_of(wave, mask, nbody(key))) for mask in unique))
        for key in TWO_BODY_KEYS)
    tables["nbody"] = dict(
        (key, dict((mask, local_of(wave, mask, nbody(key))) for mask in unique))
        for key in NBODY_KEYS)
    if calmode == 1:
        tables["anomalous"] = dict(
            (key, dict((mask, local_of(wave, mask, anomalous_operators(*key)))
                       for mask in unique)) for key in ANOMALOUS_KEYS)
    # S^2 directly from the fermion-built S-S+ matrix on M=0 states.
    basis = spin.balanced_basis(NSITE)
    matrix = spin.spin2_matrix(NSITE, basis)
    index = dict((mask, k) for k, mask in enumerate(basis))
    s2 = {}
    for mask in unique:
        column = index[mask]
        total = 0.0j
        for row, target in enumerate(basis):
            entry = int(matrix[row][column])
            if entry:
                total += entry * (wave[target] / wave[mask]).conjugate()
        s2[mask] = total
    tables["s2"] = s2
    return tables


def operator_of(group, key):
    if group == "onebody":
        return one_body(*key)
    if group in ("twobody", "nbody"):
        return nbody(key)
    return anomalous_operators(*key)


def check_readback(F, parameters):
    for i in range(NSITE):
        for j in range(NSITE):
            if F[i][j] != parameters[i * NSITE + j]:
                raise AssertionError("input readback changed F")


def replay_expected(fixture, nodes, masks, F, calmode, log_charge=None,
                    wave=None):
    """Expected production output for the sampled masks."""
    if log_charge is None:
        log_charge = fixture.log_charge
    if wave is None:
        reference = spin.independent_nodes(fixture.gauss)
        wave = spin.project_finite(F, nodes, log_charge,
                                   reference_nodes=reference)
    rational = spin.rational_finite(F, nodes)
    for mask in set(masks):
        if not spin.is_balanced(mask, NSITE):
            raise AssertionError("sampled mask {} has M != 0".format(mask))
        if rational[mask] == spin.ZERO:
            raise AssertionError(
                "sampled mask {} has an exactly zero finite-G reference "
                "(rational components {})".format(
                    mask, [spin.component_reference(F, node, mask)
                           for node in nodes]))
    unique = sorted(set(masks))
    tables = local_tables(wave, fixture.terms(), unique, calmode)
    energy = mean_over(masks, tables["energy"])
    energy2 = sum(abs(tables["energy"][mask]) ** 2 for mask in masks) / \
        float(len(masks))
    numbers = dict((mask, popcount(mask)) for mask in unique)
    number = mean_over(masks, numbers)
    number2 = sum(numbers[mask] ** 2 for mask in masks) / float(len(masks))
    expected = {
        "energy": energy,
        "energy2": energy2,
        "relative_variance": ((energy2 - energy * energy) /
                              (energy * energy)).real,
        "number": number,
        "number2": number2,
        "variance_number": number2 - number * number,
        "sz": 0.0,
        "sz2": 0.0,
        "s2": mean_over(masks, tables["s2"]),
    }
    for group in ("onebody", "twobody", "nbody", "anomalous"):
        if group in tables:
            expected[group] = dict(
                (key, mean_over(masks, values))
                for key, values in tables[group].items())
    if calmode == 0:
        expected["sr"] = expected_sr(fixture, F, nodes, masks, unique,
                                     tables["energy"], energy, log_charge)
    return expected


def expected_sr(fixture, F, nodes, masks, unique, local_energy, energy,
                log_charge):
    sr = {}
    derivatives = mp.projection_derivatives(fixture.mp, unique)

    def row(values):
        mean_o = mean_over(masks, values)
        mean_oh = sum(local_energy[mask] * values[mask]
                      for mask in masks) / float(len(masks))
        gradient = 2.0 * (mean_oh.real - energy.real * mean_o.real)
        return (mean_o.real, mean_o.imag, mean_oh.real, mean_oh.imag,
                gradient)

    for index in range(NPROJ):
        for imaginary in (0, 1):
            values = dict((mask, complex(derivatives[index][mask])
                           if not imaginary else 0.0j) for mask in unique)
            sr[2 * index + imaginary] = row(values)
    orbital = spin.log_derivatives(F, nodes, unique, log_charge)
    for parameter in range(NSITE * NSITE):
        for imaginary in (0, 1):
            values = dict((mask, complex(orbital[mask][2 * parameter +
                                                       imaginary]))
                          for mask in unique)
            sr[2 * (NPROJ + parameter) + imaginary] = row(values)
    return sr


# ---------------------------------------------------------------------------
# Production output
# ---------------------------------------------------------------------------
def read_rows(path):
    with open(path) as stream:
        return [line.split() for line in stream if line.split()]


def fused(site_text, spin_text):
    return int(site_text) + NSITE * int(spin_text)


def parse_greens(workdir, data_index, calmode):
    output = os.path.join(workdir, "output")
    tag = "{:03d}".format(data_index)
    groups = {"onebody": [], "twobody": [], "nbody": []}
    for row in read_rows(os.path.join(output,
                                      "zvo_cisajs_{}.dat".format(tag))):
        key = (fused(row[0], row[1]), fused(row[2], row[3]))
        groups["onebody"].append((key, mp.complex_at(row, 4)))
    for row in read_rows(os.path.join(output,
                                      "zvo_cisajscktalt_{}.dat".format(tag))):
        key = ((fused(row[0], row[1]), fused(row[2], row[3])),
               (fused(row[4], row[5]), fused(row[6], row[7])))
        groups["twobody"].append((key, mp.complex_at(row, 8)))
    for row in read_rows(os.path.join(output,
                                      "zvo_NBodyG_{}.dat".format(tag))):
        count = int(row[0])
        fields = row[1:1 + 4 * count]
        key = tuple((fused(fields[4 * k], fields[4 * k + 1]),
                     fused(fields[4 * k + 2], fields[4 * k + 3]))
                    for k in range(count))
        groups["nbody"].append((key, mp.complex_at(row, 1 + 4 * count)))
    if calmode == 1:
        groups["anomalous"] = []
        for row in read_rows(os.path.join(
                output, "zvo_anomalousg_{}.dat".format(tag))):
            key = (int(row[0]), fused(row[1], row[2]), fused(row[3], row[4]))
            groups["anomalous"].append((key, mp.complex_at(row, 5)))
    return groups


def compare_group(label, group, actual_rows, expected):
    seen = {}
    for key, value in actual_rows:
        if key in seen:
            raise ComparisonFailure("{} duplicated {}".format(label, key))
        seen[key] = value
    if set(seen) != set(expected) or not expected:
        raise ComparisonFailure("{} rows differ: missing={} extra={}".format(
            label, sorted(set(expected) - set(seen)),
            sorted(set(seen) - set(expected))))
    for key in sorted(expected):
        if twice_sz_change(operator_of(group, key)) != 0:
            if expected[key] != 0:
                raise AssertionError("oracle forbidden value is nonzero")
            compare("{} {} (S^z-changing)".format(label, key), seen[key],
                    0.0j, FORBIDDEN_TOLERANCE)
        else:
            compare("{} {}".format(label, key), seen[key], expected[key],
                    REPLAY_TOLERANCE * (1.0 + abs(expected[key])))


def production_s2(groups, sz2):
    twobody = dict(groups["twobody"])
    total = 0.0j
    for key in S2_TWO_BODY_KEYS:
        total += 0.5 * twobody[key]
    return total + sz2


def compare_replay(label, workdir, expected, calmode, iteration, nparameters):
    mp.compare_physical(label, workdir, expected, calmode, iteration)
    if calmode == 1:
        groups = parse_greens(workdir, 1 + iteration, calmode)
        for group in ("onebody", "twobody", "nbody", "anomalous"):
            compare_group("{} iteration {} {}".format(label, iteration, group),
                          group, groups[group], expected[group])
        physical = mp.parse_physical(workdir, calmode, iteration)
        compare("{} iteration {} S^2".format(label, iteration),
                production_s2(groups, physical["sz2"]), expected["s2"],
                REPLAY_TOLERANCE * (1.0 + abs(expected["s2"])))
    else:
        mp.compare_sr(label, os.path.join(workdir, "sr.dat"), expected["sr"],
                      nparameters, iteration)


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------
def run_with_dump(rootdir, workdir, procs=1, extra_env=None):
    environment = {"MVMC_GC_SPIN_QP_DUMP": DUMP_NAME}
    if extra_env:
        environment.update(extra_env)
    return run_binary(rootdir, workdir, procs, extra_env=environment)


def replay_case(rootdir, args):
    fixture = SpinFixture(args.boundary, args.gauss)
    threads = args.expect_threads or 1
    workdir = prepare_work(rootdir, "{}_replay_m{}_o{}_t{}".format(
        fixture.label(), args.calmode, args.store, threads))
    parameters = fixture.parameters()
    write_fixture(workdir, fixture, calmode=args.calmode, store=args.store,
                  nsplit=1, samples=REPLAY_SAMPLES, seed=REPLAY_SEED,
                  iterations=2, freeze=True)
    environment = {"MVMC_GC_STATE_DUMP": "state.dat"}
    if args.calmode == 0:
        environment["MVMC_GC_SR_DUMP"] = "sr.dat"
    if args.expect_threads:
        environment["MVMC_GC_OMP_AUDIT"] = "1"
    stdout = run_with_dump(rootdir, workdir, 1, environment)
    nodes = read_dump(workdir, fixture.gauss)
    spin.validate_coefficients(nodes, fixture.gauss)
    if args.calmode == 0:
        mp.assert_frozen_parameters(workdir, parameters)
    records = state_dump_records(os.path.join(workdir, "state.dat.rank0"))
    histogram = {}
    F = read_pair_matrix(workdir)
    check_readback(F, parameters)
    exact = spin.rational_exact(F)
    exact_zero = sorted(mask for mask, value in exact.items()
                        if value == spin.ZERO)
    # Every even-N sector of the exact S=0 state is populated, so no sector
    # of the fixture is a pure triplet (or higher) state.
    for count in range(0, NORBITAL + 1, 2):
        if all(value == spin.ZERO for mask, value in exact.items()
               if popcount(mask) == count):
            raise AssertionError("exact S=0 sector N={} is empty".format(count))
    finite_zero = sorted(mask for mask, value in
                         spin.rational_finite(F, nodes).items()
                         if value == spin.ZERO)
    for iteration, sample_records in enumerate(mp.replay_iterations(records)):
        masks = [record[3] for record in sample_records]
        if len(masks) != REPLAY_SAMPLES:
            raise AssertionError("raw sample count {}".format(len(masks)))
        for mask in masks:
            histogram[popcount(mask)] = histogram.get(popcount(mask), 0) + 1
        expected = replay_expected(fixture, nodes, masks, F, args.calmode)
        compare_replay("replay", workdir, expected, args.calmode, iteration,
                       len(parameters))
        if iteration == 1:
            mp.check_iteration_one_guards(workdir, expected, args.calmode,
                                          len(parameters))
    if args.expect_threads:
        mp.assert_omp_audit(stdout, args.expect_threads,
                            "green" if args.calmode == 1 else "cofactor")
    print("GC SP replay {} m{} o{}: 2 x {} samples passed, N histogram {}, "
          "zero masks exact={} finite={}".format(
              fixture.label(), args.calmode, args.store, REPLAY_SAMPLES,
              dict(sorted(histogram.items())), exact_zero, finite_zero))


INPUT_CASES = {
    # name: (fixture options, expected diagnostic or None for success)
    "g_zero": ({"gauss_token": "0"}, "NSPGaussLeg must be positive"),
    "g_negative": ({"gauss_token": "-1"}, "NSPGaussLeg must be positive"),
    "g_fraction": ({"gauss_token": "2.5"},
                   "NSPGaussLeg must be an exact integer"),
    "g_decimal": ({"gauss_token": "2.0"},
                  "NSPGaussLeg must be an exact integer"),
    "g_nan": ({"gauss_token": "nan"}, "NSPGaussLeg must be an exact integer"),
    "g_inf": ({"gauss_token": "inf"}, "NSPGaussLeg must be an exact integer"),
    "g_minus_inf": ({"gauss_token": "-inf"},
                    "NSPGaussLeg must be an exact integer"),
    "g_int_over": ({"gauss_token": "2147483648"},
                   "NSPGaussLeg must be an exact integer"),
    "g_huge": ({"gauss_token": "1e300"},
               "NSPGaussLeg must be an exact integer"),
    "g_long": ({"gauss_token": "1" + "0" * 299},
               "NSPGaussLeg must be an exact integer"),
    "g_extra": ({"gauss_token": "2 trailing"},
                "NSPGaussLeg must be an exact integer"),
    "g_workspace_over": ({"gauss_token": "1073741824"},
                         "NSPGaussLeg workspace"),
    # INT_MAX/7 < G <= INT_MAX/2: missed by a 2G-only bound.
    "g_workspace_sum_over": ({"gauss_token": "400000000"},
                             "NSPGaussLeg workspace"),
    "s_negative": ({"stot_token": "-1"}, "NSPStot=0"),
    "s_nonzero": ({"stot_token": "1"}, "NSPStot=0"),
    "s_fraction": ({"stot_token": "0.5"}, "NSPStot must be an exact integer"),
    "s_decimal": ({"stot_token": "0.0"}, "NSPStot must be an exact integer"),
    "s_nan": ({"stot_token": "nan"}, "NSPStot must be an exact integer"),
    "s_inf": ({"stot_token": "inf"}, "NSPStot must be an exact integer"),
    "s_int_over": ({"stot_token": "2147483648"},
                   "NSPStot must be an exact integer"),
    "general": ({"general": True},
                "spin projection requires anti-parallel pairing"),
    "combined": ({"nmptrans": 2, "translations": 2},
                 "spin projection currently requires abs(NMPTrans)=1"),
    "spin_jastrow": ({"spin_jastrow": True},
                     "spin projection does not support SpinJastrow"),
    "opttrans": ({"opttrans": True}, "NQPOptTrans<=1"),
    "positive_g2": ({"gauss_token": "2"}, None),
    "positive_g3": ({"gauss_token": "3"}, None),
    "positive_g8": ({"gauss_token": "8"}, None),
    "compat_g1_spin_jastrow": ({"gauss_token": "1", "spin_jastrow": True},
                               None),
    "compat_g1_momentum": ({"gauss_token": "1", "momentum": True,
                            "nmptrans": 4}, None),
    "compat_canonical_g2": ({"gauss_token": "2", "grand": 0}, None),
}
SIGNAL_PATTERN = re.compile(
    r"Segmentation fault|AddressSanitizer|runtime error|signal 11|"
    r"Signal: 11|SIGSEGV|Bus error")


def run_vmc(rootdir, workdir, procs, options, use_initial=True,
            extra_env=None, timeout=NEGATIVE_TIMEOUT):
    binary = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    command = [binary] + list(options) + ["-e", "namelist.def"]
    if use_initial:
        command.append("initial.def")
    if procs > 1:
        command = mpi_command(procs, command)
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    if extra_env:
        environment.update(extra_env)
    try:
        process = subprocess.run(command, cwd=workdir, env=environment,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 universal_newlines=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise AssertionError("vmc.out timed out after {} s".format(timeout))
    write(os.path.join(workdir, "run.log"), process.stdout)
    return process.returncode, process.stdout


def input_case(rootdir, args):
    if args.case not in INPUT_CASES:
        raise AssertionError("unknown input case {}".format(args.case))
    options, message = INPUT_CASES[args.case]
    options = dict(options)
    gauss_token = options.pop("gauss_token", "2")
    gauss = int(gauss_token) if re.match(r"^-?\d{1,9}$", gauss_token) else 2
    fixture = SpinFixture("pbc", max(gauss, 1))
    workdir = prepare_work(rootdir, "gc_sp_input_{}_{}_np{}".format(
        args.case, args.order, args.np))
    command_options = write_fixture(
        workdir, fixture, calmode=1, store=1, nsplit=1, samples=16, seed=11,
        iterations=1, freeze=True, order=args.order, warmup=4,
        gauss_token=gauss_token, **options)
    returncode, output = run_vmc(rootdir, workdir, args.np, command_options)
    if SIGNAL_PATTERN.search(output) or returncode < 0:
        raise AssertionError("input case {} crashed (exit {}):\n{}".format(
            args.case, returncode, output[-4000:]))
    if message is None:
        if returncode != 0 or "Finish calculation." not in output:
            raise AssertionError("positive input case {} failed (exit {}):\n"
                                 "{}".format(args.case, returncode,
                                             output[-4000:]))
    else:
        if returncode == 0:
            raise AssertionError("negative input case {} succeeded".format(
                args.case))
        if message not in output:
            raise AssertionError("input case {} lacks '{}':\n{}".format(
                args.case, message, output[-4000:]))
    print("GC SP input {} order={} np={} passed (exit {})".format(
        args.case, args.order, args.np, returncode))


def weights_case(rootdir, args):
    report = []
    for gauss in (1, 2, 3, 4, 8, 16, 32):
        fixture = SpinFixture("pbc", gauss)
        workdir = prepare_work(rootdir, "gc_sp_weights_g{}".format(gauss))
        write_fixture(workdir, fixture, calmode=1, store=1, nsplit=1,
                      samples=8, seed=5, iterations=1, freeze=True)
        run_with_dump(rootdir, workdir)
        nodes = read_dump(workdir, gauss)
        worst = spin.validate_coefficients(nodes, gauss)
        # A common-scale error in the weights must not survive.
        scaled = [dict(node, weight=2.0 * node["weight"]) for node in nodes]
        try:
            spin.validate_coefficients(scaled, gauss)
        except AssertionError:
            pass
        else:
            raise AssertionError("doubled weights were accepted")
        entry = {"gauss": gauss, "max_abs_difference": worst}
        if gauss >= 4:
            sequential = 0.0
            exact = Fraction(0)
            for node in nodes:
                weight = complex(node["weight"]).real
                sequential += weight * (node["cc"] - node["ss"])
                exact += Fraction(weight) * (Fraction(node["cc"]) -
                                             Fraction(node["ss"]))
            for key, value in (("sequential", sequential),
                               ("rational", float(exact))):
                if abs(value) > spin.COEFFICIENT_TOLERANCE:
                    raise AssertionError("G={} triplet {} sum {}".format(
                        gauss, key, value))
            entry["triplet_sequential"] = sequential
            entry["triplet_rational"] = float(exact)
        report.append(entry)
    for entry in report:
        print("GC SP weights {}".format(entry))


def dumpfail_case(rootdir, args):
    fixture = SpinFixture("pbc", 2)
    workdir = prepare_work(rootdir, "gc_sp_dumpfail_np{}".format(args.np))
    write_fixture(workdir, fixture, calmode=1, store=1, nsplit=1, samples=8,
                  seed=5, iterations=1, freeze=True)
    returncode, output = run_vmc(
        rootdir, workdir, args.np, [],
        extra_env={"MVMC_GC_SPIN_QP_DUMP":
                   os.path.join("no_such_parent", DUMP_NAME)})
    if returncode == 0 or "GC spin qp dump failed" not in output:
        raise AssertionError("dump failure not reported (exit {}):\n{}".format(
            returncode, output[-4000:]))
    if SIGNAL_PATTERN.search(output) or returncode < 0:
        raise AssertionError("dump failure crashed:\n" + output[-4000:])
    if args.np == 1:
        control = prepare_work(rootdir, "gc_sp_dump_unset")
        write_fixture(control, fixture, calmode=1, store=1, nsplit=1,
                      samples=8, seed=5, iterations=1, freeze=True)
        before = set(os.listdir(control))
        run_binary(rootdir, control, 1)
        for name in set(os.listdir(control)) - before:
            path = os.path.join(control, name)
            if os.path.isfile(path):
                with open(path, "rb") as stream:
                    if b"GC_SPIN_QP_V1" in stream.read(64):
                        raise AssertionError("dump written without the env")
    print("GC SP dump failure np={} passed".format(args.np))


def sample_masks(workdir, rank=0):
    records = state_dump_records(os.path.join(
        workdir, "state.dat.rank{}".format(rank)))
    return [record[3] for record in records if record[0] == "SAMPLE"]


def mpi_case(rootdir, args):
    """qp splitting (NSplitSize=np, every rank owns one chain) or sample
    splitting (NSplitSize=1, chains per rank) against the oracle."""
    fixture = SpinFixture(args.boundary, args.gauss)
    workdir = prepare_work(rootdir, "{}_mpi_m{}_o{}_np{}_s{}".format(
        fixture.label(), args.calmode, args.store, args.np, args.nsplit))
    parameters = fixture.parameters()
    write_fixture(workdir, fixture, calmode=args.calmode, store=args.store,
                  nsplit=args.nsplit, samples=REPLAY_SAMPLES, seed=777,
                  iterations=1, freeze=True)
    environment = {"MVMC_GC_STATE_DUMP": "state.dat"}
    if args.calmode == 0:
        environment["MVMC_GC_SR_DUMP"] = "sr.dat"
    run_with_dump(rootdir, workdir, args.np, environment)
    nodes = read_dump(workdir, fixture.gauss)
    ranks = [sample_masks(workdir, rank) for rank in range(args.np)]
    if args.nsplit == args.np:
        for rank in range(1, args.np):
            if ranks[rank] != ranks[0]:
                raise AssertionError(
                    "rank {} chain differs under qp splitting".format(rank))
        masks = ranks[0]
    else:
        if any(not rank_masks for rank_masks in ranks):
            raise AssertionError("a sample-splitting rank drew no samples")
        masks = [mask for rank_masks in ranks for mask in rank_masks]
    if any(not spin.is_balanced(mask, NSITE) for mask in masks):
        raise AssertionError("MPI run sampled M != 0")
    F = read_pair_matrix(workdir)
    check_readback(F, parameters)
    expected = replay_expected(fixture, nodes, masks, F, args.calmode)
    compare_replay("mpi", workdir, expected, args.calmode, 0, len(parameters))
    print("GC SP mpi {} m{} o{} np={} nsplit={}: {} samples ({} per rank) "
          "passed".format(fixture.label(), args.calmode, args.store, args.np,
                          args.nsplit, len(masks),
                          [len(rank_masks) for rank_masks in ranks]))


def store_case(rootdir, args):
    """NStoreO=0/1 sample the same chain and give the same SR rows."""
    fixture = SpinFixture(args.boundary, 8)
    runs = {}
    for store in (0, 1):
        workdir = prepare_work(rootdir, "{}_nstore_o{}".format(
            fixture.label(), store))
        write_fixture(workdir, fixture, calmode=0, store=store, nsplit=1,
                      samples=1024, seed=REPLAY_SEED, iterations=1,
                      freeze=True)
        run_with_dump(rootdir, workdir, 1, {"MVMC_GC_STATE_DUMP": "state.dat",
                                           "MVMC_GC_SR_DUMP": "sr.dat"})
        rows = parse_sr_dump(os.path.join(workdir, "sr.dat"))[0]["p"]
        runs[store] = (rows, sample_masks(workdir))
    if runs[0][1] != runs[1][1]:
        raise AssertionError("NStoreO changed the sampled configurations")
    first, second = runs[0][0], runs[1][0]
    if set(first) != set(second) or len(first) != 2 * (NPROJ + NSITE ** 2):
        raise AssertionError("NStoreO changed the SR row inventory")
    for index in sorted(first):
        for column, (a, b) in enumerate(zip(first[index], second[index])):
            compare("NStoreO {} P{} column {}".format(fixture.label(), index,
                                                      column),
                    a, b, 5.0e-10 * (1.0 + abs(a)))
    print("GC SP NStoreO {}: {} samples, {} SR rows agree".format(
        fixture.label(), len(runs[0][1]), len(first)))


def read_parameters(path):
    projection, orbitals = mp.read_initial_parameters(path, NSITE * NSITE)
    for value in list(projection) + list(orbitals):
        if not (math.isfinite(value.real) and math.isfinite(value.imag)):
            raise AssertionError("nonfinite optimized parameter")
    return projection, orbitals


def optimization_case(rootdir, args):
    """A real dense-SR step: 1-step and 2-step runs share the first
    iteration, and iteration 2 replays against the updated parameters."""
    fixture = SpinFixture(args.boundary, 8)
    samples = 1024
    workdirs = {}
    for steps in (1, 2):
        workdir = prepare_work(rootdir, "{}_optimize_{}step".format(
            fixture.label(), steps))
        write_fixture(workdir, fixture, calmode=0, store=1, nsplit=1,
                      samples=samples, seed=REPLAY_SEED, iterations=steps,
                      freeze=False, optimize=True)
        run_with_dump(rootdir, workdir, 1, {"MVMC_GC_STATE_DUMP": "state.dat",
                                           "MVMC_GC_SR_DUMP": "sr.dat"})
        workdirs[steps] = workdir
    one, two = workdirs[1], workdirs[2]
    first = sample_masks(one)
    groups = mp.replay_iterations(
        state_dump_records(os.path.join(two, "state.dat.rank0")), samples)
    if [record[3] for record in groups[0]] != first:
        raise AssertionError("1-step and 2-step runs differ in iteration 1")
    projection, orbitals = read_parameters(
        os.path.join(one, "output", "zqp_opt.dat"))
    initial = fixture.parameters()
    if max(abs(a - b) for a, b in zip(orbitals, initial)) <= 0.0:
        raise AssertionError("SR step did not change any orbital parameter")
    # Iteration 2 of the 2-step run uses the parameters after one step.
    nodes = read_dump(two, fixture.gauss)
    F = [[orbitals[i * NSITE + j] for j in range(NSITE)]
         for i in range(NSITE)]
    model = fixture.mp.model(()).replace(
        gutzwiller=projection[0].real,
        jastrow=tuple(value.real for value in projection[1:NPROJ]))

    def log_charge(mask):
        return ap.log_projection(model, mask)

    masks = [record[3] for record in groups[1]]
    expected = replay_expected(fixture, nodes, masks, F, 0,
                               log_charge=log_charge)
    mp.compare_physical("optimize", two, expected, 0, 1)
    mp.compare_sr("optimize", os.path.join(two, "sr.dat"), expected["sr"],
                  NSITE * NSITE, 1)
    print("GC SP optimization {}: step changed max |df|={:.3e}, iteration 2 "
          "replay passed".format(fixture.label(), max(
              abs(a - b) for a, b in zip(orbitals, initial))))


def convergence_case(rootdir, args):
    """Finite-G production quadrature converges to the exact S=0 state."""
    report = []
    for gauss in (2, 4, 8, 16, 32):
        current = SpinFixture(args.boundary, gauss)
        workdir = prepare_work(rootdir, "{}_convergence".format(
            current.label()))
        write_fixture(workdir, current, calmode=1, store=1, nsplit=1,
                      samples=8, seed=5, iterations=1, freeze=True)
        run_with_dump(rootdir, workdir)
        nodes = read_dump(workdir, gauss)
        F = read_pair_matrix(workdir)
        finite = spin.project_finite(
            F, nodes, current.log_charge,
            reference_nodes=spin.independent_nodes(gauss))
        exact = spin.project_exact(F, current.log_charge)
        basis = spin.balanced_basis(NSITE)
        error, alpha = spin.wave_error(finite, exact, basis)
        s2 = spin.spin2_expectation(finite, NSITE)
        finite_sectors = spin.sector_norms(finite)
        exact_sectors = spin.sector_norms(exact)
        finite_total = sum(finite_sectors.values())
        exact_total = sum(exact_sectors.values())
        sector_error = max(abs(finite_sectors[count] / finite_total -
                               exact_sectors[count] / exact_total)
                           for count in exact_sectors)
        report.append((gauss, error, s2, sector_error, alpha))
        print("GC SP convergence {} G={}: wave_error={:.3e} S2={:.3e} "
              "sector_error={:.3e} alpha={:.15g}".format(
                  args.boundary, gauss, error, s2, sector_error, alpha))
        if gauss >= 16 and (error > 1.0e-11 or s2 > 1.0e-11 or
                            sector_error > 1.0e-11):
            raise AssertionError("G={} has not converged to S=0".format(gauss))
    if report[0][1] < 1.0e-6:
        raise AssertionError("G=2 already exact: the fixture is not a "
                             "convergence test")


def mutation_case(rootdir, args):
    """The replay comparison rejects wrong oracles for the same samples."""
    fixture = SpinFixture("pbc", 4)
    workdir = prepare_work(rootdir, "gc_sp_mutation_src")
    write_fixture(workdir, fixture, calmode=1, store=1, nsplit=1,
                  samples=REPLAY_SAMPLES, seed=99, iterations=1, freeze=True)
    run_with_dump(rootdir, workdir, 1, {"MVMC_GC_STATE_DUMP": "state.dat"})
    nodes = read_dump(workdir, 4)
    F = read_pair_matrix(workdir)
    masks = sample_masks(workdir)

    def check(label, wave):
        expected = replay_expected(fixture, nodes, masks, F, 1, wave=wave)
        compare_replay(label, workdir, expected, 1, 0, NSITE * NSITE)

    check("control", None)
    no_sine = []
    for node in nodes:
        beta = 2.0 * math.acos(math.sqrt(node["cc"]))
        no_sine.append(dict(node, weight=node["weight"] / math.sin(beta)))
    # cc<->ss (beta -> pi-beta) and cs -> -cs (beta -> -beta) leave the M=0
    # amplitudes of a symmetric Gauss-Legendre rule unchanged, so only
    # changes that break the rule or the rotation are meaningful mutants.
    half = [dict(node, weight=2.0 * node["weight"])
            for node in nodes[:len(nodes) // 2]]
    shifted = []
    for node in nodes:
        beta = 2.0 * math.acos(math.sqrt(node["cc"])) + 0.05
        c, s = math.cos(0.5 * beta), math.sin(0.5 * beta)
        shifted.append(dict(node, cs=c * s, cc=c * c, ss=s * s))
    flipped = [dict(node, ss=-node["ss"]) for node in nodes]
    unrestricted = dict((mask, 0.0j) for mask in spin.even_masks(NSITE))
    base = spin.pair_wave(F)
    for node in nodes:
        c, s = spin.node_cos_sin(node)
        weight = spin.node_weight(node)
        for mask, value in spin.rotate_fock(base, NSITE, c, s).items():
            unrestricted[mask] += weight * value
    unrestricted = dict((mask, value * math.exp(fixture.log_charge(mask)))
                        for mask, value in unrestricted.items())
    mutants = {
        "sin(beta) dropped": spin.project_finite(F, no_sine,
                                                 fixture.log_charge),
        "half of the nodes": spin._finish(NSITE,
                                          spin.rational_finite(F, half),
                                          fixture.log_charge),
        "beta shifted": spin._finish(NSITE, spin.rational_finite(F, shifted),
                                     fixture.log_charge),
        "ss sign flipped": spin._finish(NSITE,
                                        spin.rational_finite(F, flipped),
                                        fixture.log_charge),
        "no M=0 restriction": unrestricted,
        "unprojected (G=1)": spin.project_finite(
            F, spin.independent_nodes(1), fixture.log_charge),
    }
    for label, wave in sorted(mutants.items()):
        try:
            check(label, wave)
        except (AssertionError, ComparisonFailure):
            continue
        raise AssertionError("oracle mutation '{}' was not detected".format(
            label))
    print("GC SP oracle mutations detected: {}".format(sorted(mutants)))


# ---------------------------------------------------------------------------
# Statistical gates (Linux Release pilot, fixed before validation)
# ---------------------------------------------------------------------------
PILOT_SEEDS = (1001, 2003, 3011, 4019, 5021, 6037, 7043, 8053,
               9067, 10079, 11087, 12097, 13109, 14119, 15121, 16127)
VALIDATION_SEED = 91807
ALPHA_FAMILY = 0.002699796063260207  # two-sided normal 3 sigma
PILOT_SAMPLES = 60000
PILOT_WARMUP = 1000
PILOT_GAUSS = 8
PILOT_FIXTURES = (("pbc", 1), ("apbc", 1), ("pbc", 0), ("apbc", 0))
THRESHOLD_FILE = "gc_spin_thresholds.json"
THRESHOLD_SCHEMA = 1
# Readiness-review criterion: twice the largest Phase 2 MP pilot gate of the
# same observable class.  S^2 is derived from Green rows and only recorded.
OVERSIZE_LIMITS = {"energy": 0.096, "energy2": 1.36, "number": 0.094,
                   "number2": 0.50, "green": 0.038, "anomalous": 0.032,
                   "sr": 0.26}


def student_t_cdf(t, dof):
    if dof < 1:
        raise ValueError("positive dof required")

    def integral_cos_power(theta, power):
        base = theta if power % 2 == 0 else math.sin(theta)
        start = 2 if power % 2 == 0 else 3
        for m in range(start, power + 1, 2):
            base = (math.sin(theta) * math.cos(theta) ** (m - 1) / m
                    + (m - 1) * base / m)
        return base

    theta = math.atan(abs(t) / math.sqrt(dof))
    ratio = (integral_cos_power(theta, dof - 1) /
             integral_cos_power(math.pi / 2, dof - 1))
    return 0.5 + math.copysign(0.5 * ratio, t)


def family_critical(nseed, comparisons, alpha):
    if comparisons < 1 or nseed < 2 or not 0.0 < alpha < 1.0:
        raise ValueError("invalid family-wise critical value request")
    target = 1.0 - alpha / (2.0 * comparisons)
    lo, hi = 0.0, 1.0
    while student_t_cdf(hi, nseed - 1) < target:
        hi *= 2.0
    for unused in range(100):
        mid = (lo + hi) / 2.0
        if student_t_cdf(mid, nseed - 1) < target:
            lo = mid
        else:
            hi = mid
    return hi


def summarize_pilot(values, expected, critical):
    if len(values) != len(PILOT_SEEDS) or \
            not all(math.isfinite(v) for v in values):
        raise ValueError("16 finite seed means required")
    if not (math.isfinite(expected) and math.isfinite(critical)):
        raise ValueError("finite expectation and critical value required")
    mean = statistics.mean(values)
    sd = statistics.stdev(values)
    deterministic = 2.e-9 * (1.0 + abs(expected))
    if sd == 0.0:
        if abs(mean - expected) > deterministic:
            raise AssertionError("zero-scatter biased pilot")
        return {"mean": mean, "sd": sd, "gate": deterministic}
    if abs(mean - expected) > critical * sd / 4.0:
        raise AssertionError("family-wise pilot bias gate")
    # Upward two-significant-digit rounding, fixed before validation.
    raw = max(6.0 * sd, deterministic)
    unit = 10.0 ** (math.floor(math.log10(raw)) - 1)
    gate = math.ceil(raw / unit) * unit
    return {"mean": mean, "sd": sd, "gate": gate}


def _expect_raise(label, function, kinds=(AssertionError, ValueError)):
    try:
        function()
    except kinds:
        return
    raise AssertionError("pilot statistics accepted " + label)


def pilotstats_case(rootdir, args):
    """Unit checks of the CDF, critical value and pilot summary."""
    for t in (-3.0, -0.4, 0.0, 0.7, 5.0):
        if abs(student_t_cdf(t, 1) - (0.5 + math.atan(t) / math.pi)) > 2e-15:
            raise AssertionError("dof=1 CDF")
        if abs(student_t_cdf(t, 15) + student_t_cdf(-t, 15) - 1.0) > 2e-15:
            raise AssertionError("CDF symmetry")
    previous = -1.0
    for step in range(-40, 41):
        value = student_t_cdf(step / 8.0, 15)
        if value < previous:
            raise AssertionError("CDF is not monotone")
        previous = value
    for t, reference in ((2.0, 0.9680274963576401),
                         (4.0, 0.9994203415751196)):
        if abs(student_t_cdf(t, 15) - reference) > 2e-13:
            raise AssertionError("dof=15 CDF({})".format(t))
    base = family_critical(16, 100, ALPHA_FAMILY)
    if not family_critical(16, 1000, ALPHA_FAMILY) > base > \
            family_critical(16, 10, ALPHA_FAMILY):
        raise AssertionError("critical value must grow with the family size")
    target = 1.0 - ALPHA_FAMILY / 200.0
    if abs(student_t_cdf(base, 15) - target) > 1e-12:
        raise AssertionError("critical value does not solve the CDF equation")
    scatter = [0.01 * ((-1) ** k) * (1 + k % 3) for k in range(16)]
    centered = [value - statistics.mean(scatter) for value in scatter]
    narrow = summarize_pilot([1.0 + 0.1 * v for v in centered], 1.0, base)
    wide = summarize_pilot([1.0 + v for v in centered], 1.0, base)
    if not wide["gate"] > narrow["gate"]:
        raise AssertionError("6 SD gate does not follow the scatter")
    _expect_raise("a known bias", lambda: summarize_pilot(
        [1.5 + v for v in centered], 1.0, base))
    _expect_raise("NaN", lambda: summarize_pilot(
        [float("nan")] + [1.0] * 15, 1.0, base))
    _expect_raise("a missing seed", lambda: summarize_pilot(
        [1.0] * 15, 1.0, base))
    exact = summarize_pilot([0.0] * 16, 0.0, base)
    if exact["gate"] != 2.0e-9:
        raise AssertionError("SD=0 must use the deterministic tolerance")
    _expect_raise("an SD=0 mismatch", lambda: summarize_pilot(
        [0.1] * 16, 0.0, base))
    print("GC SP pilot statistics: t_crit(15, M=100)={:.6f}".format(base))


def observable_class(key):
    head = key.split(":")[0]
    if head in ("onebody", "twobody", "nbody"):
        return "green"
    if head == "s2":
        return "s2"
    if head in ("sz", "sz2"):
        return "sz"
    return head


def physical_observations(workdir):
    observed = {}
    physical = mp.parse_physical(workdir, 1, 0)
    observed["energy:re"] = physical["energy"].real
    observed["energy:im"] = physical["energy"].imag
    for key in ("energy2", "number", "number2", "sz", "sz2"):
        observed[key] = float(physical[key])
    groups = parse_greens(workdir, 1, 1)
    s2 = production_s2(groups, physical["sz2"])
    observed["s2:re"], observed["s2:im"] = s2.real, s2.imag
    for group in ("onebody", "twobody", "nbody", "anomalous"):
        for key, value in groups[group]:
            name = "{}:{}".format(group, json.dumps(key, separators=(",", ":")))
            observed[name + ":re"] = value.real
            observed[name + ":im"] = value.imag
    return observed


def sr_observations(workdir):
    rows = parse_sr_dump(os.path.join(workdir, "sr.dat"))[0]["p"]
    observed = {}
    for index in sorted(rows):
        for column, value in enumerate(rows[index]):
            observed["sr:{}:{}".format(index, column)] = value
    return observed


def exact_reference(fixture, nodes, F, calmode):
    """|Psi|^2-weighted expectations of every compared scalar."""
    wave = spin.project_finite(F, nodes, fixture.log_charge,
                               reference_nodes=spin.independent_nodes(
                                   fixture.gauss))
    support = sorted(mask for mask, value in wave.items() if value != 0)
    norm = sum(abs(wave[mask]) ** 2 for mask in support)
    probability = dict((mask, abs(wave[mask]) ** 2 / norm) for mask in support)
    tables = local_tables(wave, fixture.terms(), support, calmode)

    def weighted(values):
        return sum(probability[mask] * values[mask] for mask in support)

    energy = weighted(tables["energy"])
    reference = {}
    if calmode == 1:
        reference["energy:re"] = energy.real
        reference["energy:im"] = energy.imag
        reference["energy2"] = weighted(dict(
            (mask, abs(tables["energy"][mask]) ** 2) for mask in support))
        counts = dict((mask, popcount(mask)) for mask in support)
        reference["number"] = weighted(counts)
        reference["number2"] = weighted(dict(
            (mask, counts[mask] ** 2) for mask in support))
        reference["sz"] = 0.0
        reference["sz2"] = 0.0
        s2 = weighted(tables["s2"])
        reference["s2:re"], reference["s2:im"] = s2.real, s2.imag
        for group in ("onebody", "twobody", "nbody", "anomalous"):
            for key, values in tables[group].items():
                name = "{}:{}".format(group, json.dumps(
                    key, separators=(",", ":")))
                value = weighted(values)
                reference[name + ":re"] = value.real
                reference[name + ":im"] = value.imag
        return reference
    derivatives = mp.projection_derivatives(fixture.mp, support)
    orbital = spin.log_derivatives(F, nodes, support, fixture.log_charge)

    def rows_for(index, values):
        mean_o = weighted(values)
        mean_oh = weighted(dict((mask, tables["energy"][mask] * values[mask])
                                for mask in support))
        gradient = 2.0 * (mean_oh.real - energy.real * mean_o.real)
        for column, value in enumerate((mean_o.real, mean_o.imag,
                                        mean_oh.real, mean_oh.imag,
                                        gradient)):
            reference["sr:{}:{}".format(index, column)] = value

    for index in range(NPROJ):
        for imaginary in (0, 1):
            rows_for(2 * index + imaginary, dict(
                (mask, complex(derivatives[index][mask]) if not imaginary
                 else 0.0j) for mask in support))
    for parameter in range(NSITE * NSITE):
        for imaginary in (0, 1):
            rows_for(2 * (NPROJ + parameter) + imaginary, dict(
                (mask, complex(orbital[mask][2 * parameter + imaginary]))
                for mask in support))
    return reference


def rayleigh_energy(fixture, wave):
    """Full-space Rayleigh energy, including Hamiltonian images of nodes."""
    norm = sum(abs(value) ** 2 for value in wave.values())
    if not norm > 0.0:
        raise AssertionError("zero-norm fixture")
    numerator = 0.0j
    for mask, amplitude in wave.items():
        if amplitude == 0.0:
            continue
        for unused, coefficient, operators in fixture.terms():
            moved = apply_ops(mask, operators)
            if moved is not None:
                target, sign = moved
                numerator += (coefficient * sign * amplitude *
                              wave[target].conjugate())
    if abs(numerator.imag) > 2.0e-11 * (1.0 + abs(numerator.real)):
        raise AssertionError("Rayleigh numerator is not real: {}".format(
            numerator))
    return numerator.real / norm


def fixture_contract_case(rootdir, args):
    """The statistical fixture has full support and a non-vacuous SR check."""
    del rootdir, args
    F = fixture_pair_amplitudes()
    nodes = spin.independent_nodes(PILOT_GAUSS)
    basis = spin.balanced_basis(NSITE)
    zeros = [mask for mask, value in spin.rational_finite(F, nodes).items()
             if value == spin.ZERO]
    if zeros or len(basis) != 70:
        raise AssertionError("fixture support is incomplete: {}".format(zeros))

    reports = []
    for boundary in ("pbc", "apbc"):
        fixture = SpinFixture(boundary, PILOT_GAUSS)
        finite = spin.project_finite(F, nodes, fixture.log_charge)
        exact = spin.project_exact(F, fixture.log_charge)
        g2 = spin.project_finite(
            F, spin.independent_nodes(2), fixture.log_charge)
        g2_error = spin.wave_error(g2, exact, basis)[0]
        if g2_error <= 1.0e-6:
            raise AssertionError("{} fixture has no finite-G dependence".format(
                boundary))

        sr = exact_reference(fixture, nodes, F, 0)
        physical = exact_reference(fixture, nodes, F, 1)
        full_energy = rayleigh_energy(fixture, finite)
        if abs(full_energy - physical["energy:re"]) > 2.0e-12 or \
                abs(physical["energy:im"]) > 2.0e-12:
            raise AssertionError("{} full-space energy differs".format(
                boundary))
        worst = 0.0
        for parameter in range(NSITE * NSITE):
            i, j = divmod(parameter, NSITE)
            for imaginary in (0, 1):
                epsilon = 1.0e-5
                delta = (1.0j if imaginary else 1.0) * epsilon
                plus = [list(row) for row in F]
                minus = [list(row) for row in F]
                plus[i][j] += delta
                minus[i][j] -= delta
                plus_wave = spin.project_finite(
                    plus, nodes, fixture.log_charge)
                minus_wave = spin.project_finite(
                    minus, nodes, fixture.log_charge)
                derivative = (rayleigh_energy(fixture, plus_wave) -
                              rayleigh_energy(fixture, minus_wave)) / \
                    (2.0 * epsilon)
                row = 2 * (NPROJ + parameter) + imaginary
                expected = sr["sr:{}:4".format(row)]
                error = abs(derivative - expected)
                worst = max(worst, error)
                if error > 2.0e-8 * (1.0 + abs(expected)):
                    raise AssertionError(
                        "{} full-energy FD row {} differs: {} vs {}".format(
                            boundary, row, derivative, expected))
        reports.append((boundary, g2_error, worst))
    print("GC SP fixture contract: support=70, {}".format(reports))


def expected_comparison_count():
    physical = (9 + 2 * NORBITAL * NORBITAL + 2 * len(TWO_BODY_KEYS) +
                2 * len(NBODY_KEYS) + 2 * len(ANOMALOUS_KEYS))
    sr = 5 * 2 * (NPROJ + NSITE * NSITE)
    return sum(physical if calmode == 1 else sr
               for unused, calmode in PILOT_FIXTURES)


def statistical_run(rootdir, boundary, calmode, seed, samples, label):
    fixture = SpinFixture(boundary, PILOT_GAUSS)
    workdir = prepare_work(rootdir, "{}_{}_m{}_seed{}".format(
        fixture.label(), label, calmode, seed))
    write_fixture(workdir, fixture, calmode=calmode, store=1, nsplit=1,
                  samples=samples, seed=seed, iterations=1, freeze=True,
                  warmup=PILOT_WARMUP)
    environment = {"MVMC_GC_SR_DUMP": "sr.dat"} if calmode == 0 else {}
    run_with_dump(rootdir, workdir, 1, environment)
    observed = (physical_observations(workdir) if calmode == 1
                else sr_observations(workdir))
    return fixture, workdir, observed


def fixture_reference(fixture, workdir, calmode):
    nodes = read_dump(workdir, fixture.gauss)
    spin.validate_coefficients(nodes, fixture.gauss)
    F = read_pair_matrix(workdir)
    check_readback(F, fixture.parameters())
    return exact_reference(fixture, nodes, F, calmode)


def file_sha256(path):
    with open(path, "rb") as stream:
        return hashlib.sha256(stream.read()).hexdigest()


def pilot_case(rootdir, args):
    """Linux Release pilot: 16 seeds x 4 fixtures; writes the gate file."""
    samples = args.samples or PILOT_SAMPLES
    comparisons = expected_comparison_count()
    critical = family_critical(len(PILOT_SEEDS), comparisons, ALPHA_FAMILY)
    jobs = []
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        for boundary, calmode in PILOT_FIXTURES:
            for seed in PILOT_SEEDS:
                jobs.append(((boundary, calmode, seed), pool.submit(
                    statistical_run, rootdir, boundary, calmode, seed,
                    samples, "pilot")))
        results = dict((key, job.result()) for key, job in jobs)
    rows = []
    raw = {}
    oversize = []
    for boundary, calmode in PILOT_FIXTURES:
        first = results[(boundary, calmode, PILOT_SEEDS[0])]
        reference = fixture_reference(first[0], first[1], calmode)
        label = "{}_m{}".format(boundary, calmode)
        keys = sorted(reference)
        for seed in PILOT_SEEDS:
            observed = results[(boundary, calmode, seed)][2]
            if sorted(observed) != keys:
                raise AssertionError("pilot {} seed {} inventory differs".format(
                    label, seed))
            raw["{}:seed{}".format(label, seed)] = observed
        for key in keys:
            values = [results[(boundary, calmode, seed)][2][key]
                      for seed in PILOT_SEEDS]
            summary = summarize_pilot(values, reference[key], critical)
            rows.append({"fixture": label, "key": key,
                         "expected": reference[key], "mean": summary["mean"],
                         "sd": summary["sd"], "gate": summary["gate"]})
            group = "sr" if key.startswith("sr:") else observable_class(key)
            limit = OVERSIZE_LIMITS.get(group)
            if limit is not None and summary["gate"] > limit:
                oversize.append((label, key, summary["gate"], limit))
    if len(rows) != comparisons:
        raise AssertionError("pilot compared {} scalars, declared {}".format(
            len(rows), comparisons))
    binary = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    report = {
        "schema_version": THRESHOLD_SCHEMA,
        "source_commit": args.source_commit,
        "binary_sha256": file_sha256(binary),
        "image_id": args.image_id,
        "toolchain": args.toolchain,
        "thread_env": {"OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
                       "OPENBLAS_NUM_THREADS":
                       os.environ.get("OPENBLAS_NUM_THREADS")},
        "seeds": list(PILOT_SEEDS),
        "validation_seed": VALIDATION_SEED,
        "samples": samples,
        "warmup": PILOT_WARMUP,
        "gauss": PILOT_GAUSS,
        "comparison_count": comparisons,
        "alpha_family": ALPHA_FAMILY,
        "t_critical": critical,
        "oversize_limits": OVERSIZE_LIMITS,
        "oversize": oversize,
        "rows": rows,
    }
    output = args.output or os.path.join(rootdir, "work", "pilot_report.json")
    write(output, json.dumps(report, indent=1, sort_keys=True) + "\n")
    write(output + ".raw.json", json.dumps(raw, sort_keys=True) + "\n")
    print("GC SP pilot: {} comparisons, t_crit={:.6f}, oversize={}".format(
        comparisons, critical, oversize))


THRESHOLD_KEYS = frozenset((
    "alpha_family", "binary_sha256", "comparison_count", "gauss",
    "image_id", "oversize", "oversize_limits", "rows", "samples",
    "schema_version", "seeds", "source_commit", "t_critical", "thread_env",
    "toolchain", "validation_seed", "warmup"))
THRESHOLD_ROW_KEYS = frozenset(
    ("fixture", "key", "expected", "mean", "sd", "gate"))
THREAD_ENV_KEYS = frozenset(("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"))


def load_thresholds(path):
    with open(path) as stream:
        data = json.load(stream)
    validate_thresholds(data)
    return data


def finite_number(value):
    return (isinstance(value, (int, float)) and
            not isinstance(value, bool) and math.isfinite(value))


def same_float(actual, expected):
    return math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-15)


def expected_threshold_keys():
    nodes = spin.independent_nodes(PILOT_GAUSS)
    keys = set()
    for boundary, calmode in PILOT_FIXTURES:
        fixture = SpinFixture(boundary, PILOT_GAUSS)
        reference = exact_reference(
            fixture, nodes, fixture_pair_amplitudes(), calmode)
        label = "{}_m{}".format(boundary, calmode)
        keys.update((label, key) for key in reference)
    return keys


def gate_from_summary(expected, mean, sd, critical):
    deterministic = 2.0e-9 * (1.0 + abs(expected))
    if sd == 0.0:
        if abs(mean - expected) > deterministic:
            raise ValueError("zero-scatter threshold row is biased")
        return deterministic
    if abs(mean - expected) > critical * sd / math.sqrt(len(PILOT_SEEDS)):
        raise ValueError("threshold row violates the pilot bias gate")
    raw = max(6.0 * sd, deterministic)
    unit = 10.0 ** (math.floor(math.log10(raw)) - 1)
    return math.ceil(raw / unit) * unit


def validate_thresholds(data):
    if not isinstance(data, dict) or set(data) != THRESHOLD_KEYS:
        raise ValueError("threshold metadata keys differ")
    if not isinstance(data["schema_version"], int) or \
            isinstance(data["schema_version"], bool) or \
            data["schema_version"] != THRESHOLD_SCHEMA:
        raise ValueError("unsupported threshold schema")
    for key in ("source_commit", "binary_sha256", "image_id", "toolchain"):
        if not isinstance(data[key], str) or not data[key]:
            raise ValueError("threshold metadata {} missing".format(key))
    if not isinstance(data["seeds"], list) or \
            any(not isinstance(seed, int) or isinstance(seed, bool) or seed <= 0
                for seed in data["seeds"]) or \
            len(set(data["seeds"])) != len(data["seeds"]) or \
            data["seeds"] != list(PILOT_SEEDS):
        raise ValueError("threshold pilot seeds differ")
    if not isinstance(data["validation_seed"], int) or \
            isinstance(data["validation_seed"], bool) or \
            data["validation_seed"] != VALIDATION_SEED or \
            data["validation_seed"] in data["seeds"]:
        raise ValueError("threshold validation seed is invalid")
    for key, expected in (("samples", PILOT_SAMPLES),
                          ("warmup", PILOT_WARMUP),
                          ("gauss", PILOT_GAUSS),
                          ("comparison_count", expected_comparison_count())):
        if not isinstance(data[key], int) or isinstance(data[key], bool) or \
                data[key] != expected:
            raise ValueError("threshold {} does not match the driver".format(
                key))
    if not finite_number(data["alpha_family"]) or \
            not 0.0 < data["alpha_family"] < 1.0 or \
            data["alpha_family"] != ALPHA_FAMILY:
        raise ValueError("threshold alpha does not match the driver")
    critical = family_critical(
        len(PILOT_SEEDS), data["comparison_count"], data["alpha_family"])
    if not finite_number(data["t_critical"]) or \
            not same_float(data["t_critical"], critical):
        raise ValueError("threshold critical value is inconsistent")
    if not isinstance(data["thread_env"], dict) or \
            set(data["thread_env"]) != THREAD_ENV_KEYS or \
            any(not isinstance(value, str) or not value
                for value in data["thread_env"].values()):
        raise ValueError("threshold thread environment is incomplete")
    if not isinstance(data["oversize_limits"], dict) or \
            set(data["oversize_limits"]) != set(OVERSIZE_LIMITS):
        raise ValueError("threshold oversize limits differ")
    for key, expected in OVERSIZE_LIMITS.items():
        value = data["oversize_limits"][key]
        if not finite_number(value) or value <= 0.0 or \
                not same_float(value, expected):
            raise ValueError("threshold oversize limit differs for " + key)
    if not isinstance(data["rows"], list) or not isinstance(data["oversize"],
                                                             list):
        raise ValueError("threshold rows/oversize must be lists")
    if not re.fullmatch(r"[0-9a-f]{40}", data["source_commit"]) or \
            not re.fullmatch(r"[0-9a-f]{64}", data["binary_sha256"]):
        raise ValueError("threshold source metadata malformed")
    seen = set()
    computed_oversize = []
    for row in data["rows"]:
        if not isinstance(row, dict) or set(row) != THRESHOLD_ROW_KEYS:
            raise ValueError("threshold row keys differ")
        if not isinstance(row["fixture"], str) or \
                not isinstance(row["key"], str):
            raise ValueError("threshold row name is malformed")
        for key in ("expected", "mean", "sd", "gate"):
            if not finite_number(row[key]):
                raise ValueError("threshold row {} is not finite".format(key))
        if row["gate"] <= 0.0 or row["sd"] < 0.0:
            raise ValueError("threshold gate must be positive")
        name = (row["fixture"], row["key"])
        if name in seen:
            raise ValueError("duplicate threshold row {}".format(name))
        seen.add(name)
        expected_gate = gate_from_summary(
            row["expected"], row["mean"], row["sd"], critical)
        if not same_float(row["gate"], expected_gate):
            raise ValueError("threshold gate is inconsistent for {}".format(
                name))
        group = ("sr" if row["key"].startswith("sr:") else
                 observable_class(row["key"]))
        limit = OVERSIZE_LIMITS.get(group)
        if limit is not None and row["gate"] > limit:
            computed_oversize.append(
                [row["fixture"], row["key"], row["gate"], limit])
    expected_keys = expected_threshold_keys()
    if seen != expected_keys or len(seen) != data["comparison_count"]:
        raise ValueError("threshold rows do not cover the exact family")
    if data["oversize"] != computed_oversize:
        raise ValueError("threshold oversize inventory is inconsistent")


def threshold_path(rootdir):
    return os.path.join(rootdir, THRESHOLD_FILE)


def threshold_contract_case(rootdir, args):
    data = load_thresholds(threshold_path(rootdir))
    mutations = []
    row = data["rows"][0]
    for label, change in (
            ("missing key", lambda d: d["rows"][0].pop("gate")),
            ("extra row key", lambda d: d["rows"][0].update(extra=1)),
            ("extra metadata", lambda d: d.update(extra=1)),
            ("negative gate", lambda d: d["rows"][0].update(gate=-1.0)),
            ("inflated gate", lambda d: d["rows"][0].update(
                gate=2.0 * d["rows"][0]["gate"])),
            ("NaN gate", lambda d: d["rows"][0].update(gate=float("nan"))),
            ("biased mean", lambda d: d["rows"][0].update(
                mean=d["rows"][0]["expected"] + 100.0)),
            ("unknown comparison", lambda d: d["rows"][0].update(
                key="unknown:observable")),
            ("missing source", lambda d: d.pop("source_commit")),
            ("missing binary", lambda d: d.update(binary_sha256="")),
            ("old schema", lambda d: d.update(schema_version=0)),
            ("boolean schema", lambda d: d.update(schema_version=True)),
            ("source suffix", lambda d: d.update(
                source_commit=d["source_commit"] + "\n")),
            ("zero alpha", lambda d: d.update(alpha_family=0.0)),
            ("wrong critical", lambda d: d.update(
                t_critical=d["t_critical"] + 0.1)),
            ("negative samples", lambda d: d.update(samples=-1)),
            ("reused validation seed", lambda d: d.update(
                validation_seed=d["seeds"][0])),
            ("wrong oversize limit", lambda d: d["oversize_limits"].update(
                sr=2.0 * d["oversize_limits"]["sr"])),
            ("false oversize", lambda d: d["oversize"].append(
                ["pbc_m0", "energy:re", 1.0, 0.096])),
            ("missing row", lambda d: d["rows"].pop())):
        mutated = json.loads(json.dumps(data))
        change(mutated)
        try:
            validate_thresholds(mutated)
        except ValueError:
            mutations.append(label)
            continue
        raise AssertionError("threshold contract accepted " + label)
    print("GC SP threshold contract: {} rows, rejected {} (first row {})".format(
        len(data["rows"]), mutations, row["key"]))


def validation_case(rootdir, args, calmode):
    data = load_thresholds(threshold_path(rootdir))
    label = "{}_m{}".format(args.boundary, calmode)
    rows = dict((row["key"], row) for row in data["rows"]
                if row["fixture"] == label)
    if not rows:
        raise AssertionError("no pilot-fixed thresholds for " + label)
    fixture, workdir, observed = statistical_run(
        rootdir, args.boundary, calmode, VALIDATION_SEED,
        data["samples"], "validation")
    reference = fixture_reference(fixture, workdir, calmode)
    if sorted(reference) != sorted(rows) or sorted(observed) != sorted(rows):
        raise AssertionError("validation inventory differs from the gates")
    for key in sorted(rows):
        if abs(reference[key] - rows[key]["expected"]) > \
                1.0e-9 * (1.0 + abs(reference[key])):
            raise AssertionError("oracle expectation changed for " + key)
        compare("{} {}".format(label, key), observed[key], reference[key],
                rows[key]["gate"])
    print("GC SP {} validation {}: {} scalars within pilot gates".format(
        "physical" if calmode == 1 else "SR", label, len(rows)))


def physical_case(rootdir, args):
    validation_case(rootdir, args, 1)


def sr_case(rootdir, args):
    validation_case(rootdir, args, 0)


def sample_case(rootdir, args):
    """The sample generator writes the replay fixture; a short run replays."""
    if not args.sample_dir:
        raise AssertionError("sample case requires --sample-dir")
    generator = os.path.join(args.sample_dir, "generate.py")
    fixture = SpinFixture("pbc", 8)
    generated = prepare_work(rootdir, "gc_sp_sample_m{}".format(args.calmode))
    expected = prepare_work(rootdir, "gc_sp_sample_expected_m{}".format(
        args.calmode))
    command = [sys.executable, generator, "--gauss", "8", "--calmode",
               str(args.calmode), "--samples", str(REPLAY_SAMPLES), "--seed",
               str(REPLAY_SEED), "--boundary", "pbc", "--warmup", "20",
               "--outdir", generated]
    process = subprocess.run(command, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             universal_newlines=True)
    if process.returncode != 0:
        raise AssertionError("sample generator failed:\n" + process.stdout)
    write_fixture(expected, fixture, calmode=args.calmode, store=1, nsplit=1,
                  samples=REPLAY_SAMPLES, seed=REPLAY_SEED, iterations=1,
                  freeze=False, optimize=(args.calmode == 0), warmup=20)
    names = sorted(name for name in os.listdir(expected)
                   if name.endswith(".def"))
    if sorted(name for name in os.listdir(generated)
              if name.endswith(".def")) != names:
        raise AssertionError("sample file inventory differs from the fixture")
    for name in names:
        with open(os.path.join(expected, name), "rb") as stream:
            reference = stream.read()
        with open(os.path.join(generated, name), "rb") as stream:
            if stream.read() != reference:
                raise AssertionError("sample {} differs from the fixture".format(
                    name))
    F = read_pair_matrix(generated)
    check_readback(F, fixture.parameters())
    environment = {"MVMC_GC_STATE_DUMP": "state.dat"}
    if args.calmode == 0:
        environment["MVMC_GC_SR_DUMP"] = "sr.dat"
    run_with_dump(rootdir, generated, 1, environment)
    nodes = read_dump(generated, 8)
    masks = sample_masks(generated)
    if len(masks) != REPLAY_SAMPLES:
        raise AssertionError("sample run drew {} samples".format(len(masks)))
    replay = replay_expected(fixture, nodes, masks, F, args.calmode)
    compare_replay("sample", generated, replay, args.calmode, 0,
                   NSITE * NSITE)
    print("GC SP sample m{}: generator matches the fixture; {} samples "
          "replayed".format(args.calmode, len(masks)))


def oracle_case(rootdir, args):
    spin.self_test()


CASES = {
    "convergence": convergence_case,
    "dumpfail": dumpfail_case,
    "fixture-contract": fixture_contract_case,
    "input": input_case,
    "mpi": mpi_case,
    "mutation": mutation_case,
    "optimization": optimization_case,
    "oracle": oracle_case,
    "physical": physical_case,
    "pilot": pilot_case,
    "pilotstats": pilotstats_case,
    "replay": replay_case,
    "sample": sample_case,
    "sr": sr_case,
    "store": store_case,
    "threshold-contract": threshold_contract_case,
    "weights": weights_case,
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=sorted(CASES))
    parser.add_argument("--boundary", choices=("pbc", "apbc"), default="pbc")
    parser.add_argument("--gauss", type=int, default=2)
    parser.add_argument("--calmode", type=int, choices=(0, 1), default=1)
    parser.add_argument("--store", type=int, choices=(0, 1), default=1)
    parser.add_argument("--np", type=int, default=1)
    parser.add_argument("--nsplit", type=int, default=1)
    parser.add_argument("--case", default=None)
    parser.add_argument("--order", choices=("early", "late"), default="early")
    parser.add_argument("--expect-threads", type=int, default=None)
    parser.add_argument("--fixture", choices=("dense", "pfonly"),
                        default="dense")
    parser.add_argument("--samples", type=int, default=None)
    parser.add_argument("--sample-dir", default=None)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--output", default=None)
    parser.add_argument("--source-commit", default=None)
    parser.add_argument("--image-id", default=None)
    parser.add_argument("--toolchain", default=None)
    args = parser.parse_args()
    CASES[args.command](os.getcwd(), args)


if __name__ == "__main__":
    main()
