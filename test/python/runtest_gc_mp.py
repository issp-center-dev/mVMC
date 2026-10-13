"""Exact-replay checks for grand-canonical momentum projection.

The four-site fixtures deliberately start from a non-translation-invariant
pair function.  The production TransSym sum and an independent Fock-space
oracle therefore both have to construct the requested momentum eigenstate.
"""
from __future__ import print_function

import argparse
import json
import math
import os
import random
import re
import shutil
import statistics

import grandcanonical_ap_oracle as oracle
from gc_component_stats import parse_component_stats
from runtest_gc import (
    fmt,
    prepare_work,
    run_binary,
    state_dump_records,
    validate_burn,
    write,
)
from runtest_gc_ap import def_text, run_expect_failure, write_hamiltonian
from runtest_gc_antiparallel import (
    ComparisonFailure,
    compare,
    compare_keyed,
    parse_sr_dump,
)


NSITE = 4
NORBITAL = 2 * NSITE
GUTZWILLER = -0.31
JASTROW = (0.23, -0.12)
HOPPING = 1.11 - 0.13j
MU = 0.35
COULOMB_INTRA = 0.35
COULOMB_INTER = 0.27
DELTA = 0.33 - 0.21j
NPROJ = 3
REPLAY_SAMPLES = 256
REPLAY_TOLERANCE = 2.0e-9
TWO_BODY_KEYS = (((0, 0), (5, 5)), ((0, 1), (1, 0)))
PILOT_SEEDS = (1001, 2003, 3011, 4019, 5021, 6037, 7043, 8053)
VALIDATION_SEED = 91807
VALIDATION_SAMPLES = 60000
# The pilot inspects at most 141 stochastic scalar components in one mode.
# With only eight seed means, a per-component 3-sigma rule has a Student-t(7)
# false-positive probability of about two percent and is therefore unsuitable
# as an all-components gate.  10.2 is the two-sided Bonferroni critical value,
# rounded upward, for a 3-sigma family-wise error rate over 141 components.
PILOT_SYSTEMATIC_T_CRITICAL = 10.2
PILOT_MAX_STOCHASTIC_COMPONENTS = 141
# Fixed by the Linux Release pilot (8 seeds x 60,000 samples).  Complex
# observables are compared componentwise against the scalar-component gates.
THRESHOLDS = {
    ("anti", "pbc", 2, False): {
        "energy": 0.046,
        "energy2": 0.35,
        "number": 0.038,
        "number2": 0.25,
        "green": 0.017,
        "anomalous": 0.016,
        "sr": {0: 0.015, 2: 0.032, 14: 0.077, 28: 0.13, 29: 0.074},
    },
    ("general", "apbc", 1, False): {
        "energy": 0.035,
        "energy2": 0.68,
        "number": 0.034,
        "number2": 0.25,
        "green": 0.016,
        "anomalous": 0.014,
        "sr": {0: 0.036, 52: 0.041, 56: 0.021, 57: 0.039, 61: 0.025},
    },
    ("general", "pbc", 0, True): {
        "energy": 0.048,
        "energy2": 0.20,
        "number": 0.047,
        "number2": 0.25,
        "green": 0.019,
        "anomalous": 0.015,
        "sr": {0: 0.057, 4: 0.036, 14: 0.052, 28: 0.037, 55: 0.032},
    },
}


def fused_fields(orbital):
    return orbital % NSITE, orbital // NSITE


def popcount(value):
    return bin(value).count("1")


def factors_text(pairs):
    fields = []
    for out_orbital, in_orbital in pairs:
        fields.extend(fused_fields(out_orbital))
        fields.extend(fused_fields(in_orbital))
    return " ".join(str(value) for value in fields)


class Fixture(object):
    def __init__(self, mode, boundary, k_index, structural_zero=False):
        if mode not in ("anti", "general"):
            raise ValueError("unknown pairing mode {}".format(mode))
        if boundary not in ("pbc", "apbc"):
            raise ValueError("unknown boundary {}".format(boundary))
        if k_index not in (0, 1, 2):
            raise ValueError("unsupported momentum index {}".format(k_index))
        self.mode = mode
        self.boundary = boundary
        self.k = k_index
        self.structural_zero = bool(structural_zero)
        edge = oracle.P if boundary == "pbc" else oracle.AP
        self.lattice = oracle.Lattice((NSITE,), (edge,))
        geometry, unused = oracle.orbital_geometry(self.lattice)
        self.general_table, self.general_pairs = oracle.expanded_table(geometry)

    def label(self):
        suffix = "_zero" if self.structural_zero else ""
        return "gc_mp_{}_{}_k{}{}".format(
            self.mode, self.boundary, self.k, suffix)

    def model(self, parameters):
        table = self.general_table
        if self.mode == "general":
            model_parameters = tuple(parameters)
        else:
            model_parameters = tuple(0.0j for unused in self.general_pairs)
        return oracle.Model(
            self.lattice, table, model_parameters,
            gutzwiller=GUTZWILLER, jastrow=JASTROW,
            hopping=HOPPING, mu=MU,
            coulomb_intra=COULOMB_INTRA,
            coulomb_inter=COULOMB_INTER,
            delta=DELTA, anomalous_pair=(3, 4),
            anomalous_green_pairs=((3, 4),))

    def matrix(self, parameters):
        if self.mode == "general":
            return oracle.pair_matrix(
                self.general_table, parameters, NORBITAL)
        rows = [parameters[first * NSITE:(first + 1) * NSITE]
                for first in range(NSITE)]
        return oracle.antiparallel_matrix(rows, NSITE)

    def matrix_builder(self):
        return lambda values: self.matrix(values)

    def log_projection(self, state):
        return oracle.log_projection(self.model(()), state)

    def terms(self):
        return oracle.model_terms(self.model(()))

    def onebody_keys(self):
        return tuple((first, second) for first in range(NORBITAL)
                     for second in range(NORBITAL))

    def twobody_keys(self):
        return TWO_BODY_KEYS

    def anomalous_rows(self):
        return oracle.anomalous_green_operators(self.model(()))

    def nparameters(self):
        return 16 if self.mode == "anti" else 28


def make_parameters(fixture):
    seed = 2741 if fixture.mode == "anti" else 6151
    rng = random.Random(seed)
    values = [complex(rng.uniform(-0.42, 0.42),
                      rng.uniform(-0.42, 0.42))
              for unused in range(fixture.nparameters())]
    if fixture.structural_zero:
        if fixture.mode == "general":
            for index, pair in enumerate(fixture.general_pairs):
                if 0 in pair:
                    values[index] = 0.0j
        else:
            for second in range(NSITE):
                values[second] = 0.0j
    return values


def write_modpara(workdir, fixture, mode_cal, samples, seed, nstore, nsplit,
                  iterations, init_nelec, nmptrans, freeze_parameters,
                  warmup):
    data_qty = iterations if mode_cal == 1 else 1
    opt_steps = iterations if mode_cal == 0 else 1
    step = 0.0 if freeze_parameters else 0.002
    two_sz = 0 if fixture.mode == "anti" else -1
    text = """--------------------
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
Ne             2
Ncond          4
2Sz            {two_sz}
NSPGaussLeg    1
NSPStot        0
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
NGrandCanonical 1
NGCInitNelec   {init_nelec}
""".format(mode=mode_cal, data_qty=data_qty, nsite=NSITE,
           two_sz=two_sz, nmptrans=nmptrans, opt_steps=opt_steps,
           step=fmt(step), samples=samples, seed=seed, nsplit=nsplit,
           nstore=nstore, init_nelec=init_nelec, warmup=warmup)
    write(os.path.join(workdir, "modpara.def"), text)


def write_translations(workdir, lattice, k_index):
    weights, mappings = oracle.translation_set(lattice, lattice.nsite,
                                                k_index)
    rows = ["{} {} {}\n".format(r, fmt(weight.real), fmt(weight.imag))
            for r, weight in enumerate(weights)]
    for r in range(lattice.nsite):
        for site in range(lattice.nsite):
            target, sign = mappings[(r, site)]
            rows.append("{} {} {} {}\n".format(r, site, target, sign))
    write(os.path.join(workdir, "qptransidx.def"),
          def_text("NQPTrans", lattice.nsite, rows))


def write_projection(workdir, fixture, optimize):
    flag = 1 if optimize else 0
    gutz_rows = ["{} 0\n".format(site) for site in range(NSITE)]
    gutz_rows.append("0 {}\n".format(flag))
    write(os.path.join(workdir, "gutzwilleridx.def"),
          def_text("NGutzwillerIdx", 1, gutz_rows, complex_type=0))
    classes = oracle.jastrow_classes(fixture.lattice)[0]
    jastrow_rows = ["{} {} {}\n".format(first, second,
                                         classes[(first, second)])
                    for first in range(NSITE) for second in range(NSITE)
                    if first != second]
    jastrow_rows.extend("{} {}\n".format(index, flag)
                        for index in range(len(JASTROW)))
    write(os.path.join(workdir, "jastrowidx.def"),
          def_text("NJastrowIdx", len(JASTROW), jastrow_rows,
                   complex_type=0))
    write(os.path.join(workdir, "locspn.def"),
          def_text("NlocalSpin", 0,
                   ["{} 0\n".format(site) for site in range(NSITE)]))


def write_general_orbital(workdir, fixture, optimize):
    rows = []
    for first, second in fixture.general_pairs:
        index, sign = fixture.general_table[(first, second)]
        site0, spin0 = fused_fields(first)
        site1, spin1 = fused_fields(second)
        rows.append("{} {} {} {} {} {}\n".format(
            site0, spin0, site1, spin1, index, sign))
    rows.extend("{} {}\n".format(index, 1 if optimize else 0)
                for index in range(len(fixture.general_pairs)))
    write(os.path.join(workdir, "orbitalidxgen.def"),
          def_text("NOrbitalIdx", len(fixture.general_pairs), rows,
                   complex_type=1))


def write_antiparallel_orbital(workdir, fixture, optimize):
    rows = []
    for first in range(NSITE):
        for second in range(NSITE):
            index = first * NSITE + second
            suffix = " 1" if fixture.boundary == "apbc" else ""
            rows.append("{} {} {}{}\n".format(
                first, second, index, suffix))
    rows.extend("{} {}\n".format(index, 1 if optimize else 0)
                for index in range(NSITE * NSITE))
    write(os.path.join(workdir, "orbitalidx.def"),
          def_text("NOrbitalIdx", NSITE * NSITE, rows, complex_type=1))


def write_measurements(workdir, fixture, mode_cal):
    one_rows = ["{}\n".format(factors_text((pair,)))
                for pair in fixture.onebody_keys()]
    write(os.path.join(workdir, "greenone.def"),
          def_text("NCisAjs", len(one_rows), one_rows))
    two_rows = ["{}\n".format(factors_text(pairs))
                for pairs in fixture.twobody_keys()]
    write(os.path.join(workdir, "greentwo.def"),
          def_text("NCisAjsCktAltDC", len(two_rows), two_rows))
    create = (1,) + fused_fields(3) + fused_fields(4) + (DELTA,)
    remove = (0,) + fused_fields(4) + fused_fields(3) + (DELTA.conjugate(),)
    anomalous_hamiltonian = [
        "{} {} {} {} {} {} {}\n".format(
            row[0], row[1], row[2], row[3], row[4],
            fmt(row[5].real), fmt(row[5].imag))
        for row in (create, remove)]
    write(os.path.join(workdir, "anomalousterm.def"),
          def_text("NAnomalousTerm", len(anomalous_hamiltonian),
                   anomalous_hamiltonian))
    names = ["        OneBodyG  greenone.def",
             "        TwoBodyG  greentwo.def"]
    anomalous_terms = []
    for kind, first, second, unused in fixture.anomalous_rows():
        anomalous_terms.append("{} {} {} {} {}\n".format(
            kind, *(fused_fields(first) + fused_fields(second))))
    if mode_cal == 1:
        write(os.path.join(workdir, "anomalousg.def"),
              def_text("NAnomalousG", len(anomalous_terms),
                       anomalous_terms))
        names.append("      AnomalousG  anomalousg.def")
    return names


def write_initial(workdir, parameters):
    values = [0.0] * 6
    for value in (GUTZWILLER,) + JASTROW:
        values.extend((value, 0.0, 0.0))
    for value in parameters:
        values.extend((value.real, value.imag, 0.0))
    write(os.path.join(workdir, "initial.def"),
          " ".join(fmt(value) for value in values) + "\n")


def write_fixture(workdir, fixture, parameters, mode_cal, samples, seed,
                  nstore, nsplit, iterations, init_nelec, optimize,
                  nmptrans=None, freeze_parameters=False, warmup=2000):
    if nmptrans is None:
        nmptrans = -NSITE if fixture.boundary == "apbc" else NSITE
    write_modpara(workdir, fixture, mode_cal, samples, seed, nstore, nsplit,
                  iterations, init_nelec, nmptrans, freeze_parameters,
                  warmup)
    write_projection(workdir, fixture, optimize)
    write_translations(workdir, fixture.lattice, fixture.k)
    if fixture.mode == "anti":
        write_antiparallel_orbital(workdir, fixture, optimize)
        orbital_name = "OrbitalAntiParallel  orbitalidx.def"
    else:
        write_general_orbital(workdir, fixture, optimize)
        orbital_name = "  OrbitalGeneral  orbitalidxgen.def"
    model = fixture.model(parameters)
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
        "   AnomalousTerm  anomalousterm.def",
    ]
    names.extend(write_measurements(workdir, fixture, mode_cal))
    write(os.path.join(workdir, "namelist.def"), "\n".join(names) + "\n")
    write_initial(workdir, parameters)


def modpara_integer(path, key):
    with open(path) as stream:
        rows = [line.split() for line in stream]
    values = [row[1] for row in rows if len(row) >= 2 and row[0] == key]
    if len(values) != 1:
        raise AssertionError("{} occurs {} times in {}".format(
            key, len(values), path))
    return int(values[0])


def read_initial_parameters(path, count):
    words = [float(word) for word in open(path).read().split()]
    expected = 6 + 3 * (NPROJ + count)
    if len(words) != expected:
        raise AssertionError("initial parameter count {} != {}".format(
            len(words), expected))
    values = [complex(words[6 + 3 * index], words[7 + 3 * index])
              for index in range(NPROJ + count)]
    return values[:NPROJ], values[NPROJ:]


def projection_data(fixture, workdir):
    weights, mappings = oracle.parse_trans_sym(
        os.path.join(workdir, "qptransidx.def"), NSITE)
    count = abs(modpara_integer(os.path.join(workdir, "modpara.def"),
                                "NMPTrans"))
    if not 1 <= count <= len(weights):
        raise AssertionError("invalid projection prefix {}".format(count))
    return weights[:count], mappings


def oracle_wave(fixture, parameters, workdir):
    unused_projection, roundtrip = read_initial_parameters(
        os.path.join(workdir, "initial.def"), fixture.nparameters())
    if len(roundtrip) != len(parameters):
        raise AssertionError("initial orbital inventory changed")
    weights, mappings = projection_data(fixture, workdir)
    return oracle.projected_wavefunction(
        fixture.matrix(roundtrip), weights, mappings, NSITE,
        fixture.log_projection)


def replay_iterations(records, samples=REPLAY_SAMPLES):
    groups = []
    current = []
    restored = False
    for record in records:
        label = record[0]
        if label == "RESTORE":
            if len(groups) != 1 or current or restored:
                raise AssertionError("unexpected restore boundary")
            restored = True
        elif label == "SAMPLE":
            if len(groups) >= 2 or (groups and not restored):
                raise AssertionError("missing restore before second iteration")
            if record[2] != popcount(record[3]):
                raise AssertionError("sample Ncur/occupation mismatch")
            current.append(record)
        elif label == "STORE":
            if len(current) != samples:
                raise AssertionError("iteration has {} samples".format(
                    len(current)))
            if [record[1] for record in current] != list(range(samples)):
                raise AssertionError("sample indices are out of order")
            groups.append(current)
            current = []
        else:
            raise AssertionError("unknown state-dump label {}".format(label))
    if len(groups) != 2 or current or not restored:
        raise AssertionError("replay needs two iterations and one restore")
    validate_burn(records)
    return groups


def mean_over_masks(masks, values):
    return sum(values[mask] for mask in masks) / float(len(masks))


def projection_derivatives(fixture, masks):
    result = {0: {}, 1: {}, 2: {}}
    classes = oracle.jastrow_classes(fixture.lattice)[0]
    for mask in masks:
        up = [(mask >> site) & 1 for site in range(NSITE)]
        down = [(mask >> (site + NSITE)) & 1 for site in range(NSITE)]
        result[0][mask] = float(sum(a * b for a, b in zip(up, down)))
        values = [0.0, 0.0]
        charges = [up[site] + down[site] - 1 for site in range(NSITE)]
        for first in range(NSITE):
            for second in range(first + 1, NSITE):
                index = classes[(first, second)]
                values[index] += charges[first] * charges[second]
        result[1][mask], result[2][mask] = values
    return result


def expected_outputs(masks, local, wave, fixture, parameters, workdir,
                     mode_cal):
    unique = sorted(set(masks))
    energy = mean_over_masks(masks, local["energy"])
    energy2 = sum(abs(local["energy"][mask]) ** 2 for mask in masks) / \
        float(len(masks))
    numbers = dict((mask, popcount(mask)) for mask in unique)
    number = mean_over_masks(masks, numbers)
    number2 = sum(numbers[mask] ** 2 for mask in masks) / float(len(masks))
    sz_values = {}
    for mask in unique:
        n_up = popcount(mask & ((1 << NSITE) - 1))
        n_down = popcount(mask >> NSITE)
        sz_values[mask] = 0.5 * (n_up - n_down)
    sz = mean_over_masks(masks, sz_values)
    sz2 = sum(sz_values[mask] ** 2 for mask in masks) / float(len(masks))
    relative_variance = ((energy2 - energy * energy) /
                         (energy * energy)).real
    expected = {
        "energy": energy,
        "energy2": energy2,
        "relative_variance": relative_variance,
        "number": number,
        "number2": number2,
        "variance_number": number2 - number * number,
        "sz": sz,
        "sz2": sz2,
    }
    for group in ("onebody", "twobody", "anomalous"):
        expected[group] = dict(
            (key, mean_over_masks(masks, values))
            for key, values in local[group].items())
    if mode_cal == 0:
        sr = {}
        derivatives = projection_derivatives(fixture, unique)
        for index in range(NPROJ):
            for imaginary in (0, 1):
                values = dict((mask, complex(derivatives[index][mask])
                               if not imaginary else 0.0j)
                              for mask in unique)
                mean_o = mean_over_masks(masks, values)
                mean_oh = sum(local["energy"][mask] * values[mask]
                              for mask in masks) / float(len(masks))
                gradient = 2.0 * (mean_oh.real -
                                  energy.real * mean_o.real)
                sr[2 * index + imaginary] = (
                    mean_o.real, mean_o.imag, mean_oh.real,
                    mean_oh.imag, gradient)
        weights, mappings = projection_data(fixture, workdir)
        builder = fixture.matrix_builder()
        for parameter in range(len(parameters)):
            for imaginary in (0, 1):
                values = oracle.parameter_log_derivative(
                    builder, parameters, weights, mappings, NSITE,
                    fixture.log_projection, unique, parameter,
                    bool(imaginary), base=wave)
                mean_o = mean_over_masks(masks, values)
                mean_oh = sum(local["energy"][mask] * values[mask]
                              for mask in masks) / float(len(masks))
                gradient = 2.0 * (mean_oh.real -
                                  energy.real * mean_o.real)
                packed = 2 * (NPROJ + parameter) + imaginary
                sr[packed] = (mean_o.real, mean_o.imag, mean_oh.real,
                              mean_oh.imag, gradient)
        expected["sr"] = sr
    return expected


def read_rows(path):
    with open(path) as stream:
        return [line.split() for line in stream if line.split()]


def complex_at(row, column):
    return complex(float(row[column]), float(row[column + 1]))


def parse_physical(workdir, mode_cal, iteration):
    gc_rows = read_rows(os.path.join(workdir, "zvo_gc.dat"))
    if iteration >= len(gc_rows):
        raise ComparisonFailure("zvo_gc.dat lacks iteration {}".format(
            iteration))
    gc = gc_rows[iteration]
    data_index = 1 + iteration if mode_cal == 1 else 1
    out_rows = read_rows(os.path.join(
        workdir, "output", "zvo_out_{:03d}.dat".format(data_index)))
    row_index = 0 if mode_cal == 1 else iteration
    if row_index >= len(out_rows):
        raise ComparisonFailure("physical output lacks iteration {}".format(
            iteration))
    out = out_rows[row_index]
    return {
        "number": float(gc[0]),
        "number2": float(gc[1]),
        "variance_number": float(gc[2]),
        "energy": complex_at(out, 0),
        "energy2": float(out[2]),
        "relative_variance": float(out[3]),
        "sz": float(out[4]),
        "sz2": float(out[5]),
    }


def relative_tolerance(expected, scale=REPLAY_TOLERANCE):
    return scale * (1.0 + abs(complex(expected)))


def compare_physical(label, workdir, expected, mode_cal, iteration=0):
    actual = parse_physical(workdir, mode_cal, iteration)
    for key in ("number", "number2", "variance_number", "energy",
                "energy2", "relative_variance", "sz", "sz2"):
        compare("{} iteration {} {}".format(label, iteration, key),
                actual[key], expected[key], relative_tolerance(expected[key]))


def parse_greens(workdir, data_index):
    output = os.path.join(workdir, "output")
    tag = "{:03d}".format(data_index)
    onebody = []
    for row in read_rows(os.path.join(output, "zvo_cisajs_{}.dat".format(tag))):
        first = int(row[0]) + NSITE * int(row[1])
        second = int(row[2]) + NSITE * int(row[3])
        onebody.append(((first, second), complex_at(row, 4)))
    twobody = []
    path = os.path.join(output, "zvo_cisajscktalt_{}.dat".format(tag))
    for row in read_rows(path):
        pairs = []
        for offset in (0, 4):
            first = int(row[offset]) + NSITE * int(row[offset + 1])
            second = int(row[offset + 2]) + NSITE * int(row[offset + 3])
            pairs.append((first, second))
        twobody.append((tuple(pairs), complex_at(row, 8)))
    anomalous = []
    path = os.path.join(output, "zvo_anomalousg_{}.dat".format(tag))
    for row in read_rows(path):
        first = int(row[1]) + NSITE * int(row[2])
        second = int(row[3]) + NSITE * int(row[4])
        anomalous.append(((int(row[0]), first, second), complex_at(row, 5)))
    return {"onebody": onebody, "twobody": twobody,
            "anomalous": anomalous}


def compare_greens(label, workdir, expected, iteration=0):
    actual = parse_greens(workdir, 1 + iteration)
    for group in ("onebody", "twobody", "anomalous"):
        compare_keyed("{} iteration {} {}".format(label, iteration, group),
                      actual[group], expected[group], REPLAY_TOLERANCE)


def parse_sr_headers(path):
    headers = []
    with open(path) as stream:
        for line in stream:
            if line.startswith("STEP "):
                columns = line.split()
                headers.append(dict(zip(columns[::2], columns[1::2])))
    return headers


def compare_sr(label, path, expected, nparameters, iteration=0):
    total = NPROJ + nparameters
    headers = parse_sr_headers(path)
    steps = parse_sr_dump(path)
    if len(headers) != len(steps) or iteration >= len(steps):
        raise ComparisonFailure("{} SR step inventory differs".format(label))
    for header in headers:
        if int(header["NPARA"]) != total or int(header["NPROJ"]) != NPROJ:
            raise ComparisonFailure("{} SR header {}".format(label, header))
    rows = steps[iteration]["p"]
    if set(rows) != set(range(2 * total)):
        raise ComparisonFailure("{} SR P rows differ".format(label))
    for index in range(2 * total):
        for column, (actual, target) in enumerate(
                zip(rows[index], expected[index])):
            compare("{} iteration {} P{} column {}".format(
                label, iteration, index, column), actual, target,
                relative_tolerance(target))


def require_comparison_failure(label, function):
    try:
        function()
    except (ComparisonFailure, AssertionError):
        return
    raise AssertionError("comparator accepted mutated {}".format(label))


def check_iteration_one_guards(workdir, expected, mode_cal, nparameters):
    physical = dict(expected)
    physical["energy"] += 1.0
    require_comparison_failure(
        "iteration-one energy",
        lambda: compare_physical("guard", workdir, physical, mode_cal, 1))
    if mode_cal == 1:
        greens = dict(expected)
        greens["onebody"] = dict(expected["onebody"])
        greens["twobody"] = dict(expected["twobody"])
        key = sorted(greens["twobody"])[0]
        greens["twobody"][key] += 1.0
        require_comparison_failure(
            "iteration-one two-body",
            lambda: compare_greens("guard", workdir, greens, 1))
    else:
        sr = dict(expected["sr"])
        row = list(sr[0])
        row[4] += 1.0
        sr[0] = tuple(row)
        require_comparison_failure(
            "iteration-one SR",
            lambda: compare_sr("guard", os.path.join(workdir, "sr.dat"),
                               sr, nparameters, 1))


def assert_frozen_parameters(workdir, parameters):
    with open(os.path.join(workdir, "modpara.def")) as stream:
        values = [line.split()[1] for line in stream
                  if line.split() and line.split()[0] == "DSROptStepDt"]
    if values != [fmt(0.0)]:
        raise AssertionError("replay did not freeze DSROptStepDt")
    path = os.path.join(workdir, "output", "zqp_opt.dat")
    unused_projection, final = read_initial_parameters(path, len(parameters))
    for index, (actual, expected) in enumerate(zip(final, parameters)):
        if abs(actual - expected) > 1.0e-14:
            raise AssertionError("frozen parameter {} changed: {} vs {}".format(
                index, actual, expected))


def stats_line(stdout):
    return parse_component_stats(stdout)


def assert_omp_audit(stdout, expected, path):
    pattern = (r"GC OMP audit: world_rank=(\d+) path=(green|cofactor) "
               r"team_min=(\d+) team_max=(\d+) workers=(\d+) work=(\d+)")
    rows = [row for row in re.findall(pattern, stdout)
            if row[0] == "0" and row[1] == path]
    if len(rows) != 1:
        raise AssertionError("missing or duplicate work audit for " + path)
    minimum, maximum, workers, work = map(int, rows[0][2:])
    if minimum != expected or maximum != expected or \
            not (2 <= workers <= expected) or work < workers:
        raise AssertionError("OpenMP work audit failed: " + str(rows[0]))


def replay_case(rootdir, args):
    fixture = Fixture(args.mode, args.boundary, args.k, args.zero)
    expected_threads = args.expect_threads or 1
    workdir = prepare_work(
        rootdir, fixture.label() + "_replay_m{}_t{}".format(
            args.calmode, expected_threads))
    parameters = make_parameters(fixture)
    write_fixture(
        workdir, fixture, parameters, mode_cal=args.calmode,
        samples=REPLAY_SAMPLES, seed=4321, nstore=1, nsplit=1,
        iterations=2, init_nelec=4, optimize=(args.calmode == 0),
        freeze_parameters=True)
    environment = {"MVMC_GC_STATE_DUMP": "state.dat"}
    if args.calmode == 0:
        environment["MVMC_GC_SR_DUMP"] = "sr.dat"
    if args.expect_threads:
        environment["MVMC_GC_OMP_AUDIT"] = "1"
    stdout = run_binary(rootdir, workdir, 1, extra_env=environment)
    records = state_dump_records(os.path.join(workdir, "state.dat.rank0"))
    wave = oracle_wave(fixture, parameters, workdir)
    if args.calmode == 0:
        assert_frozen_parameters(workdir, parameters)
    for iteration, sample_records in enumerate(replay_iterations(records)):
        masks = [record[3] for record in sample_records]
        if fixture.mode == "anti" and any(
                popcount(mask) != 2 *
                popcount(mask & ((1 << NSITE) - 1))
                for mask in masks):
            raise AssertionError("anti-parallel sample is not Sz=0")
        if any(wave[mask] == 0 for mask in masks):
            raise AssertionError("production sampled an exact projected node")
        local = oracle.local_values(
            wave, fixture.terms(), sorted(set(masks)),
            fixture.onebody_keys(), fixture.twobody_keys(),
            fixture.anomalous_rows())
        expected = expected_outputs(
            masks, local, wave, fixture, parameters, workdir, args.calmode)
        compare_physical("replay", workdir, expected, args.calmode,
                         iteration)
        if args.calmode == 1:
            compare_greens("replay", workdir, expected, iteration)
        else:
            compare_sr("replay", os.path.join(workdir, "sr.dat"),
                       expected["sr"], len(parameters), iteration)
        if iteration == 1:
            check_iteration_one_guards(
                workdir, expected, args.calmode, len(parameters))
    stats = stats_line(stdout)
    if fixture.structural_zero and stats["to_pfonly"] == 0:
        raise AssertionError("structural-zero fixture never entered PFONLY")
    if args.expect_threads:
        path = "green" if args.calmode == 1 else "cofactor"
        assert_omp_audit(stdout, args.expect_threads, path)
    print("GC MP replay {} m{}: 2 x {} samples passed".format(
        fixture.label(), args.calmode, REPLAY_SAMPLES))


def mpi_case(rootdir, args):
    boundary = "apbc" if args.mode == "general" else "pbc"
    fixture = Fixture(args.mode, boundary, 2, False)
    parameters = make_parameters(fixture)
    workdir = prepare_work(
        rootdir, fixture.label() +
        "_mpi_np{}_s{}_m{}_o{}".format(
            args.np, args.nsplit, args.calmode, args.store))
    if args.np == 3:
        nmptrans = -2 if boundary == "apbc" else 2
    else:
        nmptrans = -4 if boundary == "apbc" else 4
    write_fixture(
        workdir, fixture, parameters, mode_cal=args.calmode,
        samples=REPLAY_SAMPLES, seed=777, nstore=args.store,
        nsplit=args.nsplit, iterations=1, init_nelec=4,
        optimize=(args.calmode == 0), nmptrans=nmptrans)
    environment = {"MVMC_GC_STATE_DUMP": "state.dat"}
    if args.calmode == 0:
        environment["MVMC_GC_SR_DUMP"] = "sr.dat"
    run_binary(rootdir, workdir, args.np, extra_env=environment)
    ranks = []
    for rank in range(args.np):
        records = state_dump_records(os.path.join(
            workdir, "state.dat.rank{}".format(rank)))
        ranks.append([record[3] for record in records
                      if record[0] == "SAMPLE"])
    if args.nsplit == args.np:
        for rank in range(1, args.np):
            if ranks[rank] != ranks[0]:
                raise AssertionError(
                    "rank {} chain differs under qp splitting".format(rank))
    masks = ranks[0]
    if len(masks) != REPLAY_SAMPLES:
        raise AssertionError("MPI replay sample count {}".format(len(masks)))
    wave = oracle_wave(fixture, parameters, workdir)
    if any(wave[mask] == 0 for mask in masks):
        raise AssertionError("MPI run sampled an exact projected node")
    local = oracle.local_values(
        wave, fixture.terms(), sorted(set(masks)), fixture.onebody_keys(),
        fixture.twobody_keys(), fixture.anomalous_rows())
    expected = expected_outputs(
        masks, local, wave, fixture, parameters, workdir, args.calmode)
    compare_physical("mpi", workdir, expected, args.calmode)
    if args.calmode == 1:
        compare_greens("mpi", workdir, expected)
    else:
        compare_sr("mpi", os.path.join(workdir, "sr.dat"),
                   expected["sr"], len(parameters))
    print("GC MP mpi mode={} np={} nsplit={} calmode={} store={} passed".format(
        args.mode, args.np, args.nsplit, args.calmode, args.store))


def vacuum_case(rootdir, args):
    fixture = Fixture("general", "pbc", 2, False)
    workdir = prepare_work(rootdir, "gc_mp_vacuum_k")
    write_fixture(
        workdir, fixture, make_parameters(fixture), mode_cal=1,
        samples=100, seed=1, nstore=1, nsplit=1, iterations=1,
        init_nelec=0, optimize=False)
    returncode, output = run_expect_failure(rootdir, workdir, 1, {})
    if returncode == 0:
        raise AssertionError("vacuum projection unexpectedly succeeded")
    for needle in ("zero projected amplitude", "after 1 attempt", "Ncur=0",
                   "NGCInitNelec"):
        if needle not in output:
            raise AssertionError(
                "vacuum diagnostic lacks '{}':\n{}".format(
                    needle, output[-2000:]))
    print("GC MP vacuum diagnostic passed")


def expect_mutation_failure(label, function):
    try:
        function()
    except (AssertionError, ComparisonFailure):
        return
    raise AssertionError("mutation '{}' was not detected".format(label))


def mutation_case(rootdir, args):
    fixture = Fixture("anti", "apbc", 2, False)
    parameters = make_parameters(fixture)
    workdir = prepare_work(rootdir, "gc_mp_mutation_src")
    write_fixture(
        workdir, fixture, parameters, mode_cal=1, samples=REPLAY_SAMPLES,
        seed=99, nstore=1, nsplit=1, iterations=1, init_nelec=4,
        optimize=False)
    run_binary(rootdir, workdir, 1,
               extra_env={"MVMC_GC_STATE_DUMP": "state.dat"})
    records = state_dump_records(os.path.join(workdir, "state.dat.rank0"))
    masks = [record[3] for record in records if record[0] == "SAMPLE"]
    weights, mappings = projection_data(fixture, workdir)

    def check(label, mutated_weights, mutated_mappings):
        wave = oracle.projected_wavefunction(
            fixture.matrix(parameters), mutated_weights, mutated_mappings,
            NSITE, fixture.log_projection)
        if any(wave[mask] == 0 for mask in masks):
            raise ComparisonFailure(
                "{} makes an observed state an exact node".format(label))
        local = oracle.local_values(
            wave, fixture.terms(), sorted(set(masks)),
            fixture.onebody_keys(), fixture.twobody_keys(),
            fixture.anomalous_rows())
        expected = expected_outputs(
            masks, local, wave, fixture, parameters, workdir, 1)
        compare_physical(label, workdir, expected, 1)
        compare_greens(label, workdir, expected)

    expect_mutation_failure(
        "phase dropped",
        lambda: check("phase dropped", [0.25] * len(weights), mappings))
    no_seam = dict((key, (target, 1))
                   for key, (target, unused_sign) in mappings.items())
    expect_mutation_failure(
        "seam sign dropped",
        lambda: check("seam sign dropped", weights, no_seam))
    print("GC MP mutation detection passed")


def exact_physical(fixture, parameters, workdir):
    wave = oracle_wave(fixture, parameters, workdir)
    observables = oracle.projected_observables(
        fixture.model(parameters), wave)
    return {
        "wave": wave,
        "energy": observables["energy"],
        "energy2": oracle.support_energy2(wave, fixture.terms()),
        "number": observables["number"],
        "number2": observables["number2"],
        "onebody": observables["greens1"],
        "anomalous": observables["anomalous_g"],
    }


def exact_sr(fixture, parameters, workdir):
    wave = oracle_wave(fixture, parameters, workdir)
    support = sorted(state for state, amplitude in wave.items()
                     if amplitude != 0)
    denominator = oracle.norm(wave)
    probabilities = dict(
        (state, abs(wave[state]) ** 2 / denominator) for state in support)
    local = oracle.local_values(
        wave, fixture.terms(), support, (), (), ())

    def weighted(values):
        return sum(probabilities[state] * values[state] for state in support)

    energy = weighted(local["energy"])
    gradients = {}
    projections = projection_derivatives(fixture, support)
    for index in range(NPROJ):
        values = dict((state, complex(projections[index][state]))
                      for state in support)
        mean_o = weighted(values)
        mean_oh = weighted(dict(
            (state, local["energy"][state] * values[state])
            for state in support))
        gradients[2 * index] = 2.0 * (
            mean_oh.real - energy.real * mean_o.real)
        gradients[2 * index + 1] = 0.0
    weights, mappings = projection_data(fixture, workdir)
    builder = fixture.matrix_builder()
    for parameter in range(len(parameters)):
        for imaginary in (0, 1):
            values = oracle.parameter_log_derivative(
                builder, parameters, weights, mappings, NSITE,
                fixture.log_projection, support, parameter, bool(imaginary),
                base=wave)
            mean_o = weighted(values)
            mean_oh = weighted(dict(
                (state, local["energy"][state] * values[state])
                for state in support))
            packed = 2 * (NPROJ + parameter) + imaginary
            gradients[packed] = 2.0 * (
                mean_oh.real - energy.real * mean_o.real)
    return gradients


def scalar_deviations(prefix, actual, expected):
    difference = complex(actual) - complex(expected)
    return {prefix + ":real": difference.real,
            prefix + ":imag": difference.imag}


def compare_complex_components(label, actual, expected, tolerance):
    actual_value = complex(actual)
    expected_value = complex(expected)
    compare(label + " real", actual_value.real, expected_value.real,
            tolerance)
    compare(label + " imag", actual_value.imag, expected_value.imag,
            tolerance)


def pilot_observations(workdir, fixture, parameters, mode_cal,
                       exact_reference=None):
    if mode_cal == 1:
        exact = exact_reference or exact_physical(
            fixture, parameters, workdir)
        actual = parse_physical(workdir, 1, 0)
        observed = {}
        observed.update(scalar_deviations(
            "energy", actual["energy"], exact["energy"]))
        for key in ("energy2", "number", "number2"):
            observed[key] = float(actual[key] - exact[key])
        greens = parse_greens(workdir, 1)
        for group in ("onebody", "anomalous"):
            actual_rows = dict(greens[group])
            exact_rows = exact[group]
            if set(actual_rows) != set(exact_rows):
                raise AssertionError("pilot {} row inventory differs".format(
                    group))
            for number, key in enumerate(sorted(exact_rows)):
                observed.update(scalar_deviations(
                    "{}:{}".format(group, number),
                    actual_rows[key], exact_rows[key]))
        return observed
    exact = exact_reference or exact_sr(fixture, parameters, workdir)
    rows = parse_sr_dump(os.path.join(workdir, "sr.dat"))[0]["p"]
    if set(rows) != set(exact):
        raise AssertionError("pilot SR inventory differs")
    return dict(("sr:{}".format(index), rows[index][4] - exact[index])
                for index in sorted(exact))


def round_up_two_significant(value):
    if value == 0.0:
        return 0.0
    exponent = math.floor(math.log10(abs(value)))
    quantum = 10.0 ** (exponent - 1)
    return math.ceil(value / quantum - 1.0e-12) * quantum


def summarize_pilot(series):
    summary = {}
    stochastic_components = 0
    for key in sorted(series):
        values = series[key]
        if len(values) != len(PILOT_SEEDS) or not all(
                math.isfinite(value) for value in values):
            raise AssertionError("pilot {} has invalid raw deviations".format(
                key))
        mean = statistics.mean(values)
        standard_deviation = statistics.stdev(values)
        rounding = 2.0e-15 * max(1.0, max(abs(value) for value in values))
        if standard_deviation > rounding:
            stochastic_components += 1
        systematic_limit = (PILOT_SYSTEMATIC_T_CRITICAL *
                            standard_deviation / math.sqrt(len(values)))
        if abs(mean) > systematic_limit + rounding:
            raise AssertionError(
                "pilot {} has systematic offset {} > {}".format(
                    key, mean, systematic_limit))
        summary[key] = {
            "mean": mean,
            "sd": standard_deviation,
            "max_abs": max(abs(value) for value in values),
            "gate": round_up_two_significant(6.0 * standard_deviation),
        }
    if stochastic_components > PILOT_MAX_STOCHASTIC_COMPONENTS:
        raise AssertionError(
            "pilot has {} stochastic components; the family-wise critical "
            "value covers at most {}".format(
                stochastic_components, PILOT_MAX_STOCHASTIC_COMPONENTS))
    return summary


def grouped_pilot_gates(summary):
    groups = {}
    for key, values in summary.items():
        if key.startswith("onebody:"):
            group = "green"
        elif key.startswith("anomalous:"):
            group = "anomalous"
        elif key.startswith("energy:"):
            group = "energy"
        elif key.startswith("sr:"):
            group = key
        else:
            group = key
        groups[group] = max(groups.get(group, 0.0), values["gate"])
    return groups


def pilot_case(rootdir, args):
    fixture = Fixture(args.mode, args.boundary, args.k, args.zero)
    seeds = PILOT_SEEDS
    parameters = make_parameters(fixture)
    all_summary = {}
    raw = {}
    for mode_cal in (1, 0):
        series = {}
        exact_reference = None
        for seed in seeds:
            workdir = prepare_work(
                rootdir, "{}_pilot_m{}_seed{}".format(
                    fixture.label(), mode_cal, seed))
            write_fixture(
                workdir, fixture, parameters, mode_cal=mode_cal,
                samples=VALIDATION_SAMPLES, seed=seed, nstore=1, nsplit=1,
                iterations=1, init_nelec=4, optimize=(mode_cal == 0),
                warmup=20000)
            environment = ({"MVMC_GC_SR_DUMP": "sr.dat"}
                           if mode_cal == 0 else {})
            stdout = run_binary(rootdir, workdir, 1, extra_env=environment)
            write(os.path.join(workdir, "pilot_stdout.log"), stdout)
            write(os.path.join(workdir, "pilot_exit_code.txt"), "0\n")
            if exact_reference is None:
                exact_reference = (
                    exact_physical(fixture, parameters, workdir)
                    if mode_cal == 1 else
                    exact_sr(fixture, parameters, workdir))
            observations = pilot_observations(
                workdir, fixture, parameters, mode_cal, exact_reference)
            raw["m{}:seed{}".format(mode_cal, seed)] = observations
            for key, value in observations.items():
                series.setdefault(key, []).append(value)
        summary = summarize_pilot(series)
        all_summary["m{}".format(mode_cal)] = {
            "detail": summary,
            "gates": grouped_pilot_gates(summary),
        }
    report = {
        "fixture": fixture.label(),
        "samples": VALIDATION_SAMPLES,
        "seeds": seeds,
        "warmup": 20000,
        "systematic_t_critical": PILOT_SYSTEMATIC_T_CRITICAL,
        "systematic_max_components": PILOT_MAX_STOCHASTIC_COMPONENTS,
        "raw": raw,
        "summary": all_summary,
    }
    report_path = os.path.join(
        rootdir, "work", fixture.label() + "_pilot_summary.json")
    write(report_path, json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"fixture": fixture.label(),
                      "summary": all_summary}, sort_keys=True))


def threshold_key(fixture):
    return (fixture.mode, fixture.boundary, fixture.k,
            fixture.structural_zero)


def physical_case(rootdir, args):
    fixture = Fixture(args.mode, args.boundary, args.k, args.zero)
    limits = THRESHOLDS.get(threshold_key(fixture))
    if limits is None:
        raise AssertionError("no pilot-fixed thresholds for " + fixture.label())
    parameters = make_parameters(fixture)
    workdir = prepare_work(rootdir, fixture.label() + "_physical")
    write_fixture(
        workdir, fixture, parameters, mode_cal=1,
        samples=VALIDATION_SAMPLES, seed=VALIDATION_SEED, nstore=1,
        nsplit=1, iterations=1, init_nelec=4, optimize=False,
        warmup=20000)
    run_binary(rootdir, workdir, 1)
    exact = exact_physical(fixture, parameters, workdir)
    actual = parse_physical(workdir, 1, 0)
    compare_complex_components(
        "physical {} energy".format(fixture.label()),
        actual["energy"], exact["energy"], limits["energy"])
    for key in ("energy2", "number", "number2"):
        compare("physical {} {}".format(fixture.label(), key),
                actual[key], exact[key], limits[key])
    greens = parse_greens(workdir, 1)
    for group, limit_name in (("onebody", "green"),
                              ("anomalous", "anomalous")):
        actual_rows = dict(greens[group])
        exact_rows = exact[group]
        if set(actual_rows) != set(exact_rows):
            raise AssertionError("physical {} inventory differs".format(
                group))
        for key in sorted(exact_rows):
            compare_complex_components(
                "physical {} {} {}".format(fixture.label(), group, key),
                actual_rows[key], exact_rows[key], limits[limit_name])
    print("GC MP physical {} passed".format(fixture.label()))


def sr_case(rootdir, args):
    fixture = Fixture(args.mode, args.boundary, args.k, args.zero)
    limits = THRESHOLDS.get(threshold_key(fixture))
    if limits is None:
        raise AssertionError("no pilot-fixed thresholds for " + fixture.label())
    if len(limits["sr"]) != 5:
        raise AssertionError("SR gate must select exactly five components")
    parameters = make_parameters(fixture)
    runs = {}
    exact = None
    for store in (0, 1):
        workdir = prepare_work(
            rootdir, fixture.label() + "_sr_o{}".format(store))
        write_fixture(
            workdir, fixture, parameters, mode_cal=0,
            samples=VALIDATION_SAMPLES, seed=VALIDATION_SEED, nstore=store,
            nsplit=1, iterations=1, init_nelec=4, optimize=True,
            warmup=20000)
        run_binary(rootdir, workdir, 1,
                   extra_env={"MVMC_GC_STATE_DUMP": "state.dat",
                              "MVMC_GC_SR_DUMP": "sr.dat"})
        if exact is None:
            exact = exact_sr(fixture, parameters, workdir)
            for packed, tolerance in sorted(limits["sr"].items()):
                if abs(exact[packed]) < 4.0 * tolerance:
                    raise AssertionError(
                        "SR P{} is vacuous: |{}| < 4 * {}".format(
                            packed, exact[packed], tolerance))
        rows = parse_sr_dump(os.path.join(workdir, "sr.dat"))[0]["p"]
        masks = [record[3] for record in state_dump_records(
            os.path.join(workdir, "state.dat.rank0"))
                 if record[0] == "SAMPLE"]
        runs[store] = (rows, masks)
        for packed, tolerance in sorted(limits["sr"].items()):
            compare("SR {} NStoreO={} P{}".format(
                fixture.label(), store, packed), rows[packed][4],
                exact[packed], tolerance)
    if runs[0][1] != runs[1][1]:
        raise AssertionError("NStoreO changed sampled configurations")
    first_rows, second_rows = runs[0][0], runs[1][0]
    if set(first_rows) != set(second_rows):
        raise AssertionError("NStoreO changed SR row inventory")
    for index in sorted(first_rows):
        for column, (first, second) in enumerate(
                zip(first_rows[index], second_rows[index])):
            compare("NStoreO {} P{} column {}".format(
                fixture.label(), index, column), first, second,
                5.0e-10 * (1.0 + abs(first)))
    print("GC MP SR {} passed (components {})".format(
        fixture.label(), sorted(limits["sr"])))


def sample_case(rootdir, args):
    if not args.sample_dir:
        raise AssertionError("sample case requires --sample-dir")
    fixture = Fixture("anti", "pbc", 2, False)
    parameters = make_parameters(fixture)
    generated = prepare_work(rootdir, "gc_mp_sample_expected")
    write_fixture(
        generated, fixture, parameters, mode_cal=1, samples=10000,
        seed=VALIDATION_SEED, nstore=1, nsplit=1, iterations=1,
        init_nelec=4, optimize=False)
    generated_files = sorted(name for name in os.listdir(generated)
                             if name.endswith(".def"))
    sample_files = sorted(name for name in os.listdir(args.sample_dir)
                          if name.endswith(".def"))
    if generated_files != sample_files:
        raise AssertionError("sample file inventory differs")
    for name in generated_files:
        with open(os.path.join(generated, name), "rb") as stream:
            expected = stream.read()
        with open(os.path.join(args.sample_dir, name), "rb") as stream:
            actual = stream.read()
        if actual != expected:
            raise AssertionError("sample {} differs from generator".format(
                name))
    smoke = prepare_work(rootdir, "gc_mp_sample_smoke")
    for name in generated_files:
        shutil.copy2(os.path.join(args.sample_dir, name),
                     os.path.join(smoke, name))
    modpara = os.path.join(smoke, "modpara.def")
    with open(modpara) as stream:
        text = stream.read()
    old = "NVMCSample     10000\n"
    if text.count(old) != 1:
        raise AssertionError("sample NVMCSample line changed")
    write(modpara, text.replace(old, "NVMCSample     500\n"))
    output = run_binary(rootdir, smoke, 1)
    values = parse_physical(smoke, 1, 0)
    for key, value in values.items():
        number = complex(value)
        if not (math.isfinite(number.real) and math.isfinite(number.imag)):
            raise AssertionError("sample produced nonfinite {}".format(key))
    stats_line(output)
    print("GC MP sample matches the generator and smoke test passed")


CASES = {
    "mpi": mpi_case,
    "mutation": mutation_case,
    "pilot": pilot_case,
    "physical": physical_case,
    "replay": replay_case,
    "sample": sample_case,
    "sr": sr_case,
    "vacuum": vacuum_case,
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", choices=sorted(CASES))
    parser.add_argument("--mode", choices=("anti", "general"),
                        default="general")
    parser.add_argument("--boundary", choices=("pbc", "apbc"),
                        default="pbc")
    parser.add_argument("--k", type=int, choices=(0, 1, 2), default=0)
    parser.add_argument("--calmode", type=int, choices=(0, 1), default=1)
    parser.add_argument("--zero", action="store_true")
    parser.add_argument("--np", type=int, default=1)
    parser.add_argument("--nsplit", type=int, default=1)
    parser.add_argument("--store", type=int, default=1)
    parser.add_argument("--sample-dir", default=None)
    parser.add_argument("--expect-threads", type=int, default=None)
    args = parser.parse_args()
    CASES[args.case](os.getcwd(), args)


if __name__ == "__main__":
    main()
