"""End-to-end checks of grand-canonical sampling with NMPTrans=-1 (APBC).

Expected values come from grandcanonical_ap_oracle.py: the lattice geometry
defines the orbital and Hamiltonian tables (cross-checked against StdFace),
an independent parser reads the written files, and exact expectation values
are obtained by full enumeration of the 128 even-parity states of the
four-site ring.  The two-site cases reuse grandcanonical_exact_oracle.py.
"""
from __future__ import print_function

import math
import os
import random
import shutil
import subprocess
import sys

import grandcanonical_ap_oracle as oracle
from grandcanonical_exact_oracle import default_parameters, slater_matrix
from runtest_gc import (
    assert_close,
    fmt,
    parse_complex_rows,
    parse_sr_dump,
    prepare_work,
    read_first_line_bytes,
    read_nonempty_rows,
    run_binary,
    state_dump_records,
    write,
    write_fixture as write_two_site_fixture,
)


BAR = "=============================================\n"


def def_text(keyword, count, rows, complex_type=None):
    second = ("ComplexType {}\n".format(complex_type)
              if complex_type is not None else BAR)
    return BAR + "{} {}\n".format(keyword, count) + second + BAR + BAR + \
        "".join(rows)


def write_modpara(workdir, nsite, mode, samples, seed, nmptrans, nstore,
                  nsplit, iterations, init_nelec, data_qty, interval):
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
2Sz            -1
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
""".format(mode=mode, data_qty=data_qty, nsite=nsite, half=nsite // 2,
           nmptrans=nmptrans, iterations=iterations, interval=interval,
           samples=samples, seed=seed, nsplit=nsplit, nstore=nstore,
           init_nelec=init_nelec))


def orbital_rows(lattice, table):
    rows = []
    nsite = lattice.nsite
    for (first, second) in sorted(table):
        index, sign = table[(first, second)]
        rows.append("{} {} {} {} {} {}\n".format(
            first % nsite, first // nsite, second % nsite, second // nsite,
            index, sign))
    return rows


def write_orbital(workdir, lattice, table, nparameter, optimize):
    flags = ["{} {}\n".format(index, 1 if optimize else 0)
             for index in range(nparameter)]
    write(os.path.join(workdir, "orbitalidxgen.def"),
          def_text("NOrbitalIdx", nparameter,
                   orbital_rows(lattice, table) + flags, complex_type=1))


def write_projection(workdir, model, optimize):
    lattice = model.lattice
    nsite = lattice.nsite
    write(os.path.join(workdir, "gutzwilleridx.def"),
          def_text("NGutzwillerIdx", 1,
                   ["{} 0\n".format(site) for site in range(nsite)] +
                   ["0 {}\n".format(1 if optimize else 0)], complex_type=0))
    rows = []
    for first in range(nsite):
        for second in range(nsite):
            if first != second:
                rows.append("{} {} {}\n".format(
                    first, second, model.jastrow_class[(first, second)]))
    rows.extend("{} {}\n".format(index, 1 if optimize else 0)
                for index in range(len(model.jastrow)))
    write(os.path.join(workdir, "jastrowidx.def"),
          def_text("NJastrowIdx", len(model.jastrow), rows, complex_type=0))
    write(os.path.join(workdir, "locspn.def"),
          def_text("NlocalSpin", 0,
                   ["{} 0\n".format(site) for site in range(nsite)]))
    write(os.path.join(workdir, "qptransidx.def"),
          def_text("NQPTrans", 1,
                   ["0 1.0 0.0\n"] +
                   ["0 {} {} 1\n".format(site, site) for site in range(nsite)]))


def write_hamiltonian(workdir, model, hopping_lattice=None):
    lattice = model.lattice
    nsite = lattice.nsite
    if hopping_lattice is None:
        hopping_lattice = lattice
    rows = []
    for orbital in range(lattice.norbital):
        site, spin = orbital % nsite, orbital // nsite
        rows.append("{} {} {} {} {} {}\n".format(
            site, spin, site, spin, fmt(model.mu), fmt(0.0)))
    hopping = oracle.hopping_coefficients(hopping_lattice, model.hopping)
    for (first, second), value in sorted(hopping.items()):
        for spin in (0, 1):
            rows.append("{} {} {} {} {} {}\n".format(
                first, spin, second, spin, fmt(value.real), fmt(value.imag)))
    write(os.path.join(workdir, "trans.def"),
          def_text("NTransfer", len(rows), rows))
    write(os.path.join(workdir, "coulombintra.def"),
          def_text("NCoulombIntra", nsite,
                   ["{} {}\n".format(site, fmt(model.coulomb_intra))
                    for site in range(nsite)]))
    bonds = sorted(oracle.bond_multiplicity(lattice).items())
    write(os.path.join(workdir, "coulombinter.def"),
          def_text("NCoulombInter", len(bonds),
                   ["{} {} {}\n".format(first, second,
                                        fmt(count * model.coulomb_inter))
                    for (first, second), count in bonds]))


def write_measurements(workdir, model, mode):
    lattice = model.lattice
    nsite = lattice.nsite
    rows = []
    for first in range(lattice.norbital):
        for second in range(lattice.norbital):
            rows.append("{} {} {} {}\n".format(
                first % nsite, first // nsite, second % nsite, second // nsite))
    write(os.path.join(workdir, "greenone.def"),
          def_text("NCisAjs", len(rows), rows))
    names = ["        OneBodyG  greenone.def"]
    if model.delta is not None:
        first, second = model.anomalous_pair
        delta = model.delta
        terms = (
            (1, first % nsite, first // nsite, second % nsite, second // nsite,
             delta),
            (0, second % nsite, second // nsite, first % nsite, first // nsite,
             delta.conjugate()),
        )
        write(os.path.join(workdir, "anomalousterm.def"),
              def_text("NAnomalousTerm", 2, [
                  "{} {} {} {} {} {} {}\n".format(
                      row[0], row[1], row[2], row[3], row[4],
                      fmt(row[5].real), fmt(row[5].imag)) for row in terms]))
        names.append("   AnomalousTerm  anomalousterm.def")
        if mode == 1:
            green = []
            for kind, first, second, unused in \
                    oracle.anomalous_green_operators(model):
                green.append("{} {} {} {} {}\n".format(
                    kind, first % nsite, first // nsite, second % nsite,
                    second // nsite))
            write(os.path.join(workdir, "anomalousg.def"),
                  def_text("NAnomalousG", len(green), green))
            names.append("      AnomalousG  anomalousg.def")
    return names


def write_initial(workdir, model, parameters):
    values = [0.0] * 6
    for value in model.proj_values():
        values.extend((value, 0.0, 0.0))
    for value in parameters:
        values.extend((value.real, value.imag, 0.0))
    write(os.path.join(workdir, "initial.def"),
          " ".join(fmt(value) for value in values) + "\n")


def write_ap_fixture(workdir, model, mode=1, samples=60000, seed=60917,
                     nmptrans=-1, nstore=1, nsplit=1, iterations=2,
                     table=None, parameters=None, orbital_opt=True,
                     projection_opt=True, init_nelec=4, data_qty=1,
                     interval=4, hopping_lattice=None, measure=True):
    if table is None:
        table = model.table
    if parameters is None:
        parameters = model.parameters
    lattice = model.lattice
    write_modpara(workdir, lattice.nsite, mode, samples, seed, nmptrans,
                  nstore, nsplit, iterations, init_nelec, data_qty, interval)
    write_projection(workdir, model, projection_opt)
    write_orbital(workdir, lattice, table, len(parameters), orbital_opt)
    write_hamiltonian(workdir, model, hopping_lattice)
    names = [
        "         ModPara  modpara.def",
        "         LocSpin  locspn.def",
        "      Gutzwiller  gutzwilleridx.def",
        "         Jastrow  jastrowidx.def",
        "  OrbitalGeneral  orbitalidxgen.def",
        "        TransSym  qptransidx.def",
        "           Trans  trans.def",
        "    CoulombIntra  coulombintra.def",
        "    CoulombInter  coulombinter.def",
    ]
    if measure:
        names.extend(write_measurements(workdir, model, mode))
    write(os.path.join(workdir, "namelist.def"), "\n".join(names) + "\n")
    write_initial(workdir, model, parameters)


# ---------------------------------------------------------------------------
# Comparisons against the exact oracle
# ---------------------------------------------------------------------------

def sampled_states(records):
    return [record[3] for record in records if record[0] == "SAMPLE"]


def distribution_chi_square(model, states):
    exact = oracle.exact_observables(model)
    denominator = oracle.norm(exact["wave"])
    counts = {}
    for state in states:
        if state not in exact["wave"]:
            raise AssertionError("chain produced an odd-parity state")
        counts[state] = counts.get(state, 0) + 1
    chi_square = 0.0
    bins = 0
    for state, amplitude in exact["wave"].items():
        expected = len(states) * abs(amplitude) ** 2 / denominator
        if expected >= 5.0:
            chi_square += (counts.get(state, 0) - expected) ** 2 / expected
            bins += 1
    return chi_square, bins


def read_physical(workdir, data_index=1):
    gc_rows = read_nonempty_rows(os.path.join(workdir, "zvo_gc.dat"))
    gc = [float(value) for value in gc_rows[data_index - 1]]
    output_rows = read_nonempty_rows(os.path.join(
        workdir, "output", "zvo_out_{:03d}.dat".format(data_index)))
    output = [float(value) for value in output_rows[0]]
    return gc, output


def read_greens(workdir, nsite, data_index=1):
    rows = parse_complex_rows(os.path.join(
        workdir, "output", "zvo_cisajs_{:03d}.dat".format(data_index)), 4)
    values = {}
    for columns, value in rows:
        key = (int(columns[0]) + int(columns[1]) * nsite,
               int(columns[2]) + int(columns[3]) * nsite)
        values[key] = value
    return values


def read_anomalous(workdir, nsite, data_index=1):
    rows = parse_complex_rows(os.path.join(
        workdir, "output", "zvo_anomalousg_{:03d}.dat".format(data_index)), 5)
    values = {}
    for columns, value in rows:
        key = (int(columns[0]), int(columns[1]) + int(columns[2]) * nsite,
               int(columns[3]) + int(columns[4]) * nsite)
        values[key] = value
    return values


# Gates fixed from a pilot of 8 independent seeds (1001, 2003, 3011, 4019,
# 5021, 6037, 7043, 8053) at the production sample count of 60000.  Each gate
# is max(6 x seed-to-seed SD, 2.5 x largest |deviation|), rounded up; Green
# and AnomalousG gates bound the largest entry of each run.  The pilot SDs
# were N 0.011, N2 0.10, var(N) 0.028, E2 0.10; max |dE| 0.026, max |dG|
# 0.0057, max |dAG| 0.0048, max chi2/bin 1.22.
TOLERANCE = {
    "number": 0.065,
    "number2": 0.62,
    "variance": 0.21,
    "energy": 0.066,
    "energy2": 0.62,
    "green": 0.015,
    "anomalous": 0.012,
    "chi_square_per_bin": 2.1,
}
# First SR step, packed index 2*(NProj+orbital)+imag with NProj=3.  Pilot
# seed-to-seed SDs: P8 0.076, P12 0.034, P13 0.027, P15 0.037, P16 0.032.
GRADIENT_TOLERANCE = {8: 0.46, 12: 0.21, 13: 0.16, 15: 0.22, 16: 0.19}


def compare_physical(model, workdir, label, data_index=1, strict=True):
    """Return a list of (name, actual, expected, tolerance) failures."""
    exact = oracle.exact_observables(model)
    nsite = model.lattice.nsite
    gc, output = read_physical(workdir, data_index)
    checks = [
        ("<N>", gc[0], exact["number"], TOLERANCE["number"]),
        ("<N2>", gc[1], exact["number2"], TOLERANCE["number2"]),
        ("var(N)", gc[2], exact["variance_number"], TOLERANCE["variance"]),
        ("energy", complex(output[0], output[1]), exact["energy"],
         TOLERANCE["energy"]),
        ("energy2", output[2], exact["energy2"].real, TOLERANCE["energy2"]),
    ]
    greens = read_greens(workdir, nsite, data_index)
    if set(greens) != set(exact["greens1"]):
        raise AssertionError("{}: one-body inventory differs".format(label))
    for key in sorted(greens):
        checks.append(("G{}".format(key), greens[key], exact["greens1"][key],
                       TOLERANCE["green"]))
    if model.delta is not None:
        anomalous = read_anomalous(workdir, nsite, data_index)
        if set(anomalous) != set(exact["anomalous_g"]):
            raise AssertionError("{}: AnomalousG inventory differs".format(label))
        for key in sorted(anomalous):
            checks.append(("AG{}".format(key), anomalous[key],
                           exact["anomalous_g"][key], TOLERANCE["anomalous"]))
    failures = [(name, actual, expected, tolerance)
                for name, actual, expected, tolerance in checks
                if abs(actual - expected) > tolerance]
    if strict and failures:
        raise AssertionError("{} mismatches: {}".format(label, failures[:6]))
    worst = max(abs(actual - expected) / tolerance
                for unused, actual, expected, tolerance in checks)
    return failures, worst


def boundary_green_keys(lattice):
    """One-body keys whose bond crosses an anti-periodic seam."""
    keys = []
    hopping = oracle.hopping_coefficients(lattice, 1.0 + 0.0j)
    for (first, second), value in hopping.items():
        if value.real < 0.0:
            for spin in (0, 1):
                keys.append((first + spin * lattice.nsite,
                             second + spin * lattice.nsite))
    return keys


# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------

def audit_env(path="gc_audit.dat"):
    return {"MVMC_GC_INPUT_AUDIT": path}


def check_audit(workdir, lattice, table, nmptrans, nslater):
    header, orbital, trans, weights = oracle.parse_audit(
        os.path.join(workdir, "gc_audit.dat"))
    expected_ap = 1 if nmptrans < 0 else 0
    if header.get("ap_flag") != expected_ap or header.get("nmptrans") != 1:
        raise AssertionError("audit header mismatch: {}".format(header))
    if header.get("nsite") != lattice.nsite or header.get("nslater") != nslater:
        raise AssertionError("audit dimensions mismatch: {}".format(header))
    if expected_ap:
        expected = dict(table)
    else:
        # The periodic reader overwrites the sign column with +1.
        expected = dict((pair, (index, 1))
                        for pair, (index, unused) in table.items())
    if orbital != expected:
        raise AssertionError("audit orbital table differs from expected")
    negative = oracle.negative_sign_count(expected)
    if header.get("negative_orbital_input_sign_count") != negative:
        raise AssertionError("audit negative orbital count {} != {}".format(
            header.get("negative_orbital_input_sign_count"), negative))
    identity = dict(((0, site), (site, 1)) for site in range(lattice.nsite))
    if trans != identity or header.get("negative_qptrans_sign_count") != 0:
        raise AssertionError("audit TransSym is not the identity")
    if weights != {0: 1.0 + 0.0j}:
        raise AssertionError("audit TransSym weight is not 1")
    # The production read-back of the first TransSym pattern equals the file.
    file_weights, mappings = oracle.parse_trans_sym(
        os.path.join(workdir, "qptransidx.def"), lattice.nsite)
    first = dict(((transform, site), value)
                 for (transform, site), value in mappings.items()
                 if transform == 0)
    if first != trans or file_weights[0] != weights[0]:
        raise AssertionError("audit TransSym differs from the input file")
    return header


def legacy_parser_case(rootdir):
    """Two-site fixture: the independent parser reproduces the 2-site oracle
    and the production read-back, for both NMPTrans=+1 and -1."""
    for nmptrans in (1, -1):
        workdir = prepare_work(rootdir, "GC_AP_LegacyParser_{:+d}".format(
            nmptrans))
        write_two_site_fixture(workdir, mode=1, samples=200, seed=11,
                               iterations=1, nmptrans=nmptrans)
        table, nparameter = oracle.parse_orbital_general(
            os.path.join(workdir, "orbitalidxgen.def"), 2)
        unused, parameters = oracle.parse_initial(
            os.path.join(workdir, "initial.def"), 1, nparameter)
        matrix = oracle.pair_matrix(table, parameters, 4)
        reference = slater_matrix(default_parameters())
        worst = max(abs(matrix[row][column] - reference[row][column])
                    for row in range(4) for column in range(4))
        if worst != 0.0:
            raise AssertionError("independent parser differs from 2-site "
                                 "oracle by {}".format(worst))
        run_binary(rootdir, workdir, extra_env=audit_env())
        check_audit(workdir, oracle.Lattice((2,), (oracle.P,)), table,
                    nmptrans, nparameter)
    print("GC AP legacy 2-site parser/audit compatibility passed")


def identity_wiring_case(rootdir):
    """NMPTrans=+1 and -1 with all-positive signs are the same computation."""
    for delta in (None, 0.35 - 0.20j):
        runs = {}
        for nmptrans in (1, -1):
            workdir = prepare_work(rootdir, "GC_AP_Wiring_{}_{:+d}".format(
                "anomalous" if delta is not None else "plain", nmptrans))
            write_two_site_fixture(workdir, mode=1, samples=4000, seed=52711,
                                   iterations=1, delta=delta,
                                   nmptrans=nmptrans)
            run_binary(rootdir, workdir,
                       extra_env={"MVMC_GC_STATE_DUMP": "gc_state.dat"})
            runs[nmptrans] = workdir
        files = ["zvo_gc.dat", os.path.join("output", "zvo_out_001.dat"),
                 os.path.join("output", "zvo_cisajs_001.dat"),
                 os.path.join("output", "zvo_cisajscktalt_001.dat"),
                 os.path.join("output", "zvo_NBodyG_001.dat"),
                 "gc_state.dat.rank0"]
        if delta is not None:
            files.append(os.path.join("output", "zvo_anomalousg_001.dat"))
        for name in files:
            with open(os.path.join(runs[1], name), "rb") as stream:
                periodic = stream.read()
            with open(os.path.join(runs[-1], name), "rb") as stream:
                antiperiodic = stream.read()
            if not periodic or periodic != antiperiodic:
                raise AssertionError("NMPTrans=+1/-1 differ in {}".format(name))
    print("GC AP identity wiring (+1/-1, all signs +1) byte-identical")


def audit_case(rootdir):
    model = oracle.ring_model()
    for nmptrans in (-1, 1):
        workdir = prepare_work(rootdir, "GC_AP_Audit_{:+d}".format(nmptrans))
        write_ap_fixture(workdir, model, samples=50, nmptrans=nmptrans,
                         iterations=1)
        run_binary(rootdir, workdir, extra_env=audit_env())
        parsed, nparameter = oracle.parse_orbital_general(
            os.path.join(workdir, "orbitalidxgen.def"), model.lattice.nsite)
        if parsed != model.table:
            raise AssertionError("written orbital table differs from geometry")
        header = check_audit(workdir, model.lattice, model.table, nmptrans,
                             nparameter)
        print("audit NMPTrans={:+d}: {}".format(nmptrans, header))
    if oracle.negative_sign_count(model.table) <= 0:
        raise AssertionError("AP fixture has no negative upper-triangle sign")
    print("GC AP input audit passed")


def expanded_reference_case(rootdir):
    """Shared AP parameters vs one parameter per pair with the sign absorbed.

    Both runs fix every orbital (OptFlag=0) and read explicit initial values,
    so InitParameter consumes no random numbers for orbitals and the GC chain
    sees the same random sequence.
    """
    model = oracle.ring_model(delta=0.33 - 0.21j)
    expanded, unused = oracle.expanded_table(model.table)
    expanded_values = oracle.expanded_parameters(model.table, model.parameters)
    runs = {}
    for label, nmptrans, table, parameters in (
            ("shared", -1, model.table, model.parameters),
            ("expanded", 1, expanded, expanded_values)):
        workdir = prepare_work(rootdir, "GC_AP_Expanded_{}".format(label))
        write_ap_fixture(workdir, model, samples=6000, seed=31337,
                         nmptrans=nmptrans, table=table, parameters=parameters,
                         orbital_opt=False, iterations=1)
        run_binary(rootdir, workdir,
                   extra_env={"MVMC_GC_STATE_DUMP": "gc_state.dat"})
        runs[label] = workdir
    shared = sampled_states(state_dump_records(
        os.path.join(runs["shared"], "gc_state.dat.rank0")))
    reference = sampled_states(state_dump_records(
        os.path.join(runs["expanded"], "gc_state.dat.rank0")))
    if not shared or shared != reference:
        raise AssertionError("expanded reference occupation sequence differs")
    for name in ("zvo_gc.dat", os.path.join("output", "zvo_out_001.dat"),
                 os.path.join("output", "zvo_cisajs_001.dat"),
                 os.path.join("output", "zvo_anomalousg_001.dat")):
        with open(os.path.join(runs["shared"], name), "rb") as stream:
            first = stream.read()
        with open(os.path.join(runs["expanded"], name), "rb") as stream:
            second = stream.read()
        if not first or first != second:
            raise AssertionError("expanded reference differs in {}".format(name))
    print("GC AP shared/expanded reference byte-identical: samples={}".format(
        len(shared)))


def physical_run(rootdir, name, model, samples, seed, **kwargs):
    workdir = prepare_work(rootdir, name)
    write_ap_fixture(workdir, model, mode=1, samples=samples, seed=seed,
                     iterations=1, **kwargs)
    run_binary(rootdir, workdir,
               extra_env={"MVMC_GC_STATE_DUMP": "gc_state.dat",
                          "MVMC_GC_DEBUG_REBUILD_INTERVAL": "97"})
    return workdir


def physical_case(rootdir, delta=None):
    samples = int(os.environ.get("MVMC_GC_AP_SAMPLES", "60000"))
    seed = int(os.environ.get("MVMC_GC_AP_SEED", "60917"))
    model = oracle.ring_model(delta=delta)
    label = "GC_AP_Physical" if delta is None else "GC_AP_Anomalous"
    workdir = physical_run(rootdir, label, model, samples, seed)
    states = sampled_states(state_dump_records(
        os.path.join(workdir, "gc_state.dat.rank0")))
    chi_square, bins = distribution_chi_square(model, states)
    if chi_square > TOLERANCE["chi_square_per_bin"] * bins:
        raise AssertionError("distribution chi-square {} over {} bins".format(
            chi_square, bins))
    if len(set(states)) != len(oracle.even_basis(8)):
        raise AssertionError("chain did not visit all 128 even states")
    unused, worst = compare_physical(model, workdir, label)
    seam = boundary_green_keys(model.lattice)
    if len(seam) != 4:
        raise AssertionError("seam Green inventory changed: {}".format(seam))
    print("{} passed: samples={} chi2/bins={:.4g}/{} worst/tol={:.3g}".format(
        label, len(states), chi_square, bins, worst))


def sr_case(rootdir):
    samples = int(os.environ.get("MVMC_GC_AP_SR_SAMPLES", "60000"))
    model = oracle.ring_model()
    nproj = len(model.proj_values())
    # Parameter 3 (ud, d=3) is shared by +1 and -1 rows; 4-5 are same-spin.
    selected = ((3, False), (3, True), (1, False), (4, True), (5, False))
    exact = {}
    for orbital, imaginary in selected:
        packed = 2 * (nproj + orbital) + (1 if imaginary else 0)
        gradient = oracle.parameter_gradient(model, orbital, imaginary)
        for epsilon in (1.0e-5, 5.0e-6):
            fd = oracle.finite_difference(model, orbital, imaginary, epsilon)
            assert_close("AP exact FD {} {}".format(packed, epsilon), fd,
                         gradient, 1.0e-7 + 1.0e-5 * abs(gradient))
        if abs(gradient) < 4.0 * GRADIENT_TOLERANCE[packed]:
            raise AssertionError("selected AP gradient P={} is vacuous: {}"
                                 .format(packed, gradient))
        exact[packed] = gradient
    runs = {}
    paths = {}
    for nstore in (0, 1):
        workdir = prepare_work(rootdir, "GC_AP_SR_Store{}".format(nstore))
        write_ap_fixture(workdir, model, mode=0, samples=samples, seed=44821,
                         nstore=nstore, iterations=2, measure=False)
        run_binary(rootdir, workdir, extra_env={"MVMC_GC_SR_DUMP": "gc_sr.dat"})
        runs[nstore] = parse_sr_dump(os.path.join(workdir, "gc_sr.dat"))
        paths[nstore] = workdir
    for packed, gradient in sorted(exact.items()):
        for nstore in (0, 1):
            assert_close("AP sampled SR gradient P={} NStore={}".format(
                packed, nstore), runs[nstore][0]["p"][packed][4], gradient,
                GRADIENT_TOLERANCE[packed])
        # NStoreO=0 accumulates per sample while NStoreO=1 contracts the stored
        # O matrix with BLAS, so the two agree to roundoff (about 1e-12 with
        # OpenBLAS on Linux), not bitwise; GC_SR_Oracle uses the same gate.
        assert_close("NStoreO AP SR gradient P={}".format(packed),
                     runs[0][0]["p"][packed][4], runs[1][0]["p"][packed][4],
                     5.0e-10)
        print("AP SR P={} exact={:.10g} sampled={:.10g}".format(
            packed, gradient, runs[0][0]["p"][packed][4]))
    if read_first_line_bytes(os.path.join(paths[0], "zvo_gc.dat")) != \
            read_first_line_bytes(os.path.join(paths[1], "zvo_gc.dat")):
        raise AssertionError("NStoreO changed AP GC output")
    print("GC AP SR NStoreO/exact-gradient oracle passed")


def mpi_case(rootdir, nsplit, mode):
    # NSplitSize=1 gives each of the two ranks its own chain: 2 x 30000 keeps
    # the combined sample count of the pilot that fixed TOLERANCE.
    samples = 30000 if mode == 1 else 12000
    model = oracle.ring_model(delta=0.33 - 0.21j if mode == 1 else None)
    serial = prepare_work(rootdir, "GC_AP_MPI_serial_s{}_m{}".format(nsplit, mode))
    parallel = prepare_work(rootdir, "GC_AP_MPI_parallel_s{}_m{}".format(
        nsplit, mode))
    for workdir in (serial, parallel):
        write_ap_fixture(workdir, model, mode=mode, samples=samples,
                         seed=67231, nsplit=nsplit, iterations=2,
                         measure=(mode == 1))
    env = {"MVMC_GC_STATE_DUMP": "gc_state.dat"}
    env.update(audit_env())
    run_binary(rootdir, serial, extra_env=env)
    run_binary(rootdir, parallel, procs=2, extra_env=env)
    check_audit(parallel, model.lattice, model.table, -1, len(model.parameters))
    if os.path.exists(os.path.join(parallel, "gc_audit.dat.rank1")):
        raise AssertionError("audit was written by a non-zero rank")
    serial_states = sampled_states(state_dump_records(
        os.path.join(serial, "gc_state.dat.rank0")))
    rank_states = [sampled_states(state_dump_records(os.path.join(
        parallel, "gc_state.dat.rank{}".format(rank)))) for rank in (0, 1)]
    if nsplit == 2:
        for rank, states in enumerate(rank_states):
            if states != serial_states:
                raise AssertionError("NSplitSize=2 rank {} chain differs".format(
                    rank))
        with open(os.path.join(serial, "zvo_gc.dat"), "rb") as stream:
            serial_gc = stream.read()
        with open(os.path.join(parallel, "zvo_gc.dat"), "rb") as stream:
            parallel_gc = stream.read()
        if serial_gc != parallel_gc:
            raise AssertionError("NSplitSize=2 GC output is not byte-identical")
    elif mode == 1:
        chi_square, bins = distribution_chi_square(
            model, rank_states[0] + rank_states[1])
        if chi_square > TOLERANCE["chi_square_per_bin"] * bins:
            raise AssertionError("MPI chi-square {} over {} bins".format(
                chi_square, bins))
        compare_physical(model, parallel, "MPI NSplitSize=1")
    print("GC AP MPI passed: NSplitSize={} CalMode={}".format(nsplit, mode))


def mutation_case(rootdir):
    """Each sign error must be caught by the production-size comparison."""
    samples = int(os.environ.get("MVMC_GC_AP_SAMPLES", "60000"))
    model = oracle.ring_model()
    exact = oracle.exact_observables(model)
    mutated_table = dict((pair, (index, 1))
                         for pair, (index, unused) in model.table.items())
    periodic = oracle.drop_boundary_phase(model.lattice)
    mutations = (
        ("hopping_sign_dropped", {"hopping_lattice": periodic},
         oracle.exact_observables(model, terms=oracle.model_terms(
             model.replace(lattice=periodic)))),
        ("orbital_sign_plus", {"table": mutated_table},
         oracle.exact_observables(model, table=mutated_table)),
    )
    for name, kwargs, mutated in mutations:
        shift = abs(mutated["energy"] - exact["energy"])
        if shift < 4.0 * TOLERANCE["energy"]:
            raise AssertionError("{} energy shift {} is below detection".format(
                name, shift))
        workdir = physical_run(rootdir, "GC_AP_Mutation_{}".format(name),
                               model, samples, 60917, **kwargs)
        failures, worst = compare_physical(model, workdir, name, strict=False)
        if not failures:
            raise AssertionError("{} mutation was not detected".format(name))
        print("mutation {} detected: dE_exact={:.4g} failures={} "
              "worst/tol={:.3g} first={}".format(
                  name, shift, len(failures), worst, failures[0][0]))
    print("GC AP sign mutations detected")


# ---------------------------------------------------------------------------
# Two-dimensional input inventory
# ---------------------------------------------------------------------------

def stdface_inputs(rootdir, workdir, lattice):
    """Generate the same Hubbard lattice with StdFace (vmcdry.out)."""
    if len(lattice.lengths) == 1:
        text = 'model = "Hubbard"\nlattice = "chain"\nL = {}\n'.format(
            lattice.lengths[0])
        phases = "phase0 = {}\n".format(180.0 if lattice.theta[0] < 0 else 0.0)
    else:
        text = 'model = "Hubbard"\nlattice = "square"\nW = {}\nL = {}\n'.format(
            lattice.lengths[0], lattice.lengths[1])
        phases = "phase0 = {}\nphase1 = {}\n".format(
            180.0 if lattice.theta[0] < 0 else 0.0,
            180.0 if lattice.theta[1] < 0 else 0.0)
    text += ("t = 1.0\nU = 4.0\n" + phases +
             "ComplexType = 1\n2Sz = -1\nncond = {}\n".format(lattice.nsite))
    stddir = os.path.join(workdir, "stdface")
    os.makedirs(stddir)
    write(os.path.join(stddir, "stan.in"), text)
    binary = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmcdry.out")
    process = subprocess.run([binary, "stan.in"], cwd=stddir,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True)
    if process.returncode != 0:
        raise AssertionError("vmcdry.out failed:\n{}".format(process.stdout))
    return stddir


def lattice_case(rootdir, lengths, boundaries):
    lattice = oracle.Lattice(lengths, boundaries)
    table, keys = oracle.orbital_geometry(lattice)
    classes, nclass = oracle.jastrow_classes(lattice)
    rng = random.Random(9001 + lattice.nsite)
    parameters = []
    for key in keys:
        if key == ("zero",):
            parameters.append(0.0j)
        else:
            parameters.append(complex(rng.uniform(-0.3, 0.3),
                                      rng.uniform(-0.3, 0.3)))
    model = oracle.Model(lattice, table, parameters, gutzwiller=-0.2,
                         jastrow=tuple(0.05 * (number % 3 - 1)
                                       for number in range(nclass)),
                         hopping=1.0 + 0.0j, mu=0.4, coulomb_intra=4.0,
                         coulomb_inter=0.0)
    name = "GC_AP_Lattice_{}".format(lattice.label())
    workdir = prepare_work(rootdir, name)
    write_ap_fixture(workdir, model, samples=200, seed=271,
                     init_nelec=lattice.nsite, iterations=1, measure=False)

    # Orbital signs: geometry vs file vs StdFace (one parameter per pair).
    parsed, nparameter = oracle.parse_orbital_general(
        os.path.join(workdir, "orbitalidxgen.def"), lattice.nsite)
    if parsed != table or nparameter != len(keys):
        raise AssertionError("{}: written orbital table differs".format(name))
    stddir = stdface_inputs(rootdir, workdir, lattice)
    stdface_table, unused = oracle.parse_orbital_general(
        os.path.join(stddir, "orbitalidxgen.def"), lattice.nsite)
    if set(stdface_table) != set(table):
        raise AssertionError("{}: StdFace pair inventory differs".format(name))
    for pair, (unused, sign) in stdface_table.items():
        if sign != oracle.pair_wrap_sign(lattice, *pair):
            raise AssertionError("{}: StdFace sign differs at {}".format(
                name, pair))
    # Hamiltonian: geometry vs file vs StdFace, coefficient by coefficient.
    expected = {}
    for (first, second), value in oracle.hopping_coefficients(
            lattice, model.hopping).items():
        for spin in (0, 1):
            expected[((first, spin), (second, spin))] = value
    written = oracle.parse_trans(os.path.join(workdir, "trans.def"))
    for site in range(lattice.nsite):
        for spin in (0, 1):
            diagonal = written.pop(((site, spin), (site, spin)))
            assert_close("{} mu".format(name), diagonal, model.mu, 0.0)
    stdface = oracle.parse_trans(os.path.join(stddir, "trans.def"))
    for label, actual in (("written", written), ("StdFace", stdface)):
        keys_union = set(actual) | set(expected)
        worst = max(abs(actual.get(key, 0.0j) - expected.get(key, 0.0j))
                    for key in keys_union)
        if worst > 1.0e-12:
            raise AssertionError("{}: {} Hamiltonian differs by {}".format(
                name, label, worst))
    if oracle.translation_defect(lattice, oracle.pair_matrix(
            table, parameters, lattice.norbital)) > 1.0e-14:
        raise AssertionError("{}: orbitals are not translation covariant"
                             .format(name))

    # Production read-back and a short serial/MPI run.
    env = audit_env()
    env["MVMC_GC_STATE_DUMP"] = "gc_state.dat"
    run_binary(rootdir, workdir, extra_env=env)
    check_audit(workdir, lattice, table, -1, len(keys))
    parallel = prepare_work(rootdir, name + "_mpi")
    write_ap_fixture(parallel, model, samples=200, seed=271,
                     init_nelec=lattice.nsite, iterations=1, measure=False)
    run_binary(rootdir, parallel, procs=2, extra_env=env)
    check_audit(parallel, lattice, table, -1, len(keys))
    for directory in (workdir, parallel):
        gc, output = read_physical(directory)
        if not all(math.isfinite(value) for value in gc + output):
            raise AssertionError("{}: non-finite output".format(name))
    print("{} inventory passed: pairs={} params={} negative={} "
          "seam_bonds={} zero_pairs={}".format(
              name, len(table), len(keys), oracle.negative_sign_count(table),
              sum(1 for value in oracle.hopping_coefficients(
                  lattice, 1.0 + 0.0j).values() if value.real < 0) // 2,
              sum(1 for index, unused in table.values()
                  if keys[index] == ("zero",))))


def stdface_example_case(rootdir):
    """samples/GrandCanonical/APBC_chain is reproducible and runs."""
    sample = os.environ["MVMC_GC_APBC_SAMPLE"]
    workdir = prepare_work(rootdir, "GC_AP_StdFaceExample")
    regenerated = os.path.join(workdir, "expert")
    process = subprocess.run(
        [sys.executable, os.path.join(sample, "make_gc_apbc_input.py"),
         os.path.join(rootdir, "..", "..", "src", "mVMC", "vmcdry.out"),
         regenerated],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True)
    if process.returncode != 0:
        raise AssertionError("sample conversion failed:\n" + process.stdout)
    committed = os.path.join(sample, "expert")
    names = sorted(os.listdir(committed))
    if names != sorted(os.listdir(regenerated)):
        raise AssertionError("sample file inventory differs: {}".format(names))
    for name in names:
        with open(os.path.join(committed, name), "rb") as stream:
            expected = stream.read()
        with open(os.path.join(regenerated, name), "rb") as stream:
            actual = stream.read()
        if actual != expected:
            raise AssertionError("sample {} is not reproduced by StdFace"
                                 .format(name))

    lattice = oracle.RING
    table, unused = oracle.parse_orbital_general(
        os.path.join(committed, "orbitalidxgen.def"), lattice.nsite)
    for pair, (unused, sign) in table.items():
        if sign != oracle.pair_wrap_sign(lattice, *pair):
            raise AssertionError("sample sign differs from geometry at {}"
                                 .format(pair))
    trans = oracle.parse_trans(os.path.join(committed, "trans.def"))
    expected = dict((((site, spin), (site, spin)), 2.0 + 0.0j)
                    for site in range(lattice.nsite) for spin in (0, 1))
    for (first, second), value in oracle.hopping_coefficients(
            lattice, 1.0 + 0.0j).items():
        for spin in (0, 1):
            expected[((first, spin), (second, spin))] = value
    if set(trans) != set(expected) or max(
            abs(trans[key] - expected[key]) for key in trans) > 1.0e-12:
        raise AssertionError("sample Hamiltonian differs from geometry")

    run = os.path.join(workdir, "run")
    shutil.copytree(committed, run)
    lines = []
    with open(os.path.join(run, "modpara.def")) as stream:
        for line in stream:
            words = line.split()
            short = {"NSROptItrStep": "20", "NSROptItrSmp": "10",
                     "NVMCSample": "200"}
            if words and words[0] in short:
                line = "{:<15}{}\n".format(words[0], short[words[0]])
            lines.append(line)
    write(os.path.join(run, "modpara.def"), "".join(lines))
    output = run_binary(rootdir, run, extra_env=audit_env())
    if "Finish calculation." not in output:
        raise AssertionError("sample run did not finish")
    check_audit(run, lattice, table, -1, len(table))
    print("GC AP StdFace example reproduced and ran: pairs={} negative={}"
          .format(len(table), oracle.negative_sign_count(table)))


CASES = {
    "stdface_example": stdface_example_case,
    "legacy_parser": legacy_parser_case,
    "identity_wiring": identity_wiring_case,
    "audit": audit_case,
    "expanded_reference": expanded_reference_case,
    "physical": physical_case,
    "anomalous": lambda root: physical_case(root, delta=0.33 - 0.21j),
    "sr": sr_case,
    "mpi_nsplit2_mode1": lambda root: mpi_case(root, 2, 1),
    "mpi_nsplit1_mode1": lambda root: mpi_case(root, 1, 1),
    "mpi_nsplit2_mode0": lambda root: mpi_case(root, 2, 0),
    "mpi_nsplit1_mode0": lambda root: mpi_case(root, 1, 0),
    "mutation": mutation_case,
    "lattice_4x2": lambda root: lattice_case(root, (4, 2),
                                             (oracle.AP, oracle.P)),
    "lattice_4x4": lambda root: lattice_case(root, (4, 4),
                                             (oracle.AP, oracle.AP)),
}


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in CASES:
        print("usage: {} {}".format(sys.argv[0], "|".join(sorted(CASES))))
        return 2
    CASES[sys.argv[1]](os.getcwd())
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print("ERROR: {}".format(error))
        sys.exit(1)
