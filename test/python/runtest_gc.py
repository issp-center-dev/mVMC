from __future__ import print_function

import math
import os
import shutil
import subprocess
import sys

from grandcanonical_exact_oracle import (
    ANOMALOUS_G_OPERATORS,
    EVEN_BASIS,
    NSITE,
    NORBITAL,
    adjoint,
    default_parameters,
    exact_observables,
    expectation,
    finite_difference,
    model_terms,
    nbody,
    parameter_gradient,
    wavefunction,
)


def write(path, text):
    with open(path, "w") as stream:
        stream.write(text)


def fmt(value):
    return "{:.18e}".format(value)


def fused_to_fields(orbital):
    return orbital % NSITE, orbital // NSITE


def operators_to_factors(operators):
    if len(operators) % 2:
        raise AssertionError("operator list has odd length")
    factors = []
    for index in range(0, len(operators), 2):
        create_op, annihilate_op = operators[index:index + 2]
        if create_op[0] != "c" or annihilate_op[0] != "a":
            raise AssertionError("operator is not in c-dagger/c factor order")
        factors.extend(fused_to_fields(create_op[1]))
        factors.extend(fused_to_fields(annihilate_op[1]))
    return factors


def write_modpara(workdir, mode, sample_count, seed, nstore=1,
                  nsplit=1, reweight=0, iterations=2, init_nelec=2,
                  data_qty=1):
    write(
        os.path.join(workdir, "modpara.def"),
        """--------------------
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
Nsite          2
Ne             1
Ncond          3
2Sz            -1
NSPGaussLeg    1
NSPStot        0
NMPTrans       1
NSROptItrStep  {iterations}
NSROptItrSmp   1
DSROptRedCut   0.000000000001
DSROptStaDel   0.01
DSROptStepDt   0.002
NVMCWarmUp     2000
NVMCInterval   4
NVMCSample     {sample_count}
NExUpdatePath  0
RndSeed        {seed}
NSplitSize     {nsplit}
NStore         {nstore}
NSRCG          0
reweight       {reweight}
NGrandCanonical 1
NGCInitNelec   {init_nelec}
""".format(
            mode=mode,
            iterations=iterations,
            sample_count=sample_count,
            seed=seed,
            nsplit=nsplit,
            nstore=nstore,
            reweight=reweight,
            init_nelec=init_nelec,
            data_qty=data_qty,
        ),
    )


def write_orbital(workdir, parameters, optimize=True):
    rows = []
    index = 0
    for first in range(NORBITAL):
        for second in range(first + 1, NORBITAL):
            site0, spin0 = fused_to_fields(first)
            site1, spin1 = fused_to_fields(second)
            rows.append("{:3d} {:2d} {:3d} {:2d} {:3d} 1".format(
                site0, spin0, site1, spin1, index))
            index += 1
    flags = "\n".join("{:3d} {}".format(index, 1 if optimize else 0)
                      for index in range(len(parameters)))
    write(
        os.path.join(workdir, "orbitalidxgen.def"),
        """=============================================
NOrbitalIdx          6
ComplexType          1
=============================================
=============================================
{rows}
{flags}
""".format(rows="\n".join(rows), flags=flags),
    )


def write_jastrow(workdir, optimize=True):
    write(
        os.path.join(workdir, "jastrowidx.def"),
        """=============================================
NJastrowIdx          1
ComplexType          0
=============================================
=============================================
0 1 0
1 0 0
0 {flag}
""".format(flag=1 if optimize else 0),
    )


def write_initial(workdir, parameters, jastrow=None):
    values = [0.0] * 6
    if jastrow is not None:
        values.extend((jastrow, 0.0, 0.0))
    for value in parameters:
        values.extend((value.real, value.imag, 0.0))
    write(os.path.join(workdir, "initial.def"),
          " ".join(fmt(value) for value in values) + "\n")


def write_identity_projection(workdir):
    write(
        os.path.join(workdir, "qptransidx.def"),
        """=============================================
NQPTrans          1
=============================================
======== TrIdx_TrWeight_and_TrIdx_i_xi ======
=============================================
0 1.0 0.0
0 0 0 1
0 1 1 1
""",
    )
    write(
        os.path.join(workdir, "locspn.def"),
        """================================
NlocalSpin     0
================================
========i_0LocSpn_1IteElc======
================================
0 0
1 0
""",
    )


def write_hamiltonian(workdir, mu):
    trans = []
    for orbital in range(NORBITAL):
        site, spin = fused_to_fields(orbital)
        trans.append((site, spin, site, spin, complex(mu)))
    hopping = 0.37 + 0.19j
    for spin in (0, 1):
        trans.append((0, spin, 1, spin, hopping))
        trans.append((1, spin, 0, spin, hopping.conjugate()))
    write(
        os.path.join(workdir, "trans.def"),
        "========================\nNTransfer {}\n========================\n"
        "========i_j_s_tijs======\n========================\n{}".format(
            len(trans),
            "".join("{} {} {} {} {} {}\n".format(
                a, b, c, d, fmt(value.real), fmt(value.imag))
                    for a, b, c, d, value in trans)),
    )
    write(os.path.join(workdir, "coulombintra.def"),
          "=============================================\nNCoulombIntra 2\n"
          "=============================================\n"
          "============= CoulombIntra ==================\n"
          "=============================================\n"
          "0 0.73\n1 -0.28\n")
    write(os.path.join(workdir, "coulombinter.def"),
          "=============================================\nNCoulombInter 1\n"
          "=============================================\n"
          "============= CoulombInter ==================\n"
          "=============================================\n0 1 0.21\n")
    write(os.path.join(workdir, "hund.def"),
          "=============================================\nNHund 1\n"
          "=============================================\n"
          "================ Hund =======================\n"
          "=============================================\n0 1 0.17\n")
    write(os.path.join(workdir, "pairhopp.def"),
          "=============================================\nNPairHopp 1\n"
          "=============================================\n"
          "================ PairHop ====================\n"
          "=============================================\n0 1 0.13\n")
    write(os.path.join(workdir, "exchange.def"),
          "=============================================\nNExchange 1\n"
          "=============================================\n"
          "================ Exchange ===================\n"
          "=============================================\n0 1 -0.11\n")

    inter_rows = []
    nbody_rows = []
    for label, coefficient, operators in model_terms(mu):
        factors = operators_to_factors(operators)
        if label == "inter_all":
            inter_rows.append("{} {} {}\n".format(
                " ".join(str(value) for value in factors),
                fmt(coefficient.real), fmt(coefficient.imag)))
        elif label == "nbody_inter_all":
            nbody_rows.append("{} {} {} {}\n".format(
                len(factors) // 4,
                " ".join(str(value) for value in factors),
                fmt(coefficient.real), fmt(coefficient.imag)))
    write(os.path.join(workdir, "interall.def"),
          "=============================================\nNInterAll {}\n"
          "=============================================\n"
          "================ InterAll ===================\n"
          "=============================================\n{}".format(len(inter_rows), "".join(inter_rows)))
    write(os.path.join(workdir, "nbodyinterall.def"),
          "=============================================\nNNBodyInterAll {}\n"
          "=============================================\n"
          "============= NBodyInterAll =================\n"
          "=============================================\n{}".format(len(nbody_rows), "".join(nbody_rows)))


def write_measurements(workdir):
    one_rows = []
    for out_orbital in range(NORBITAL):
        for in_orbital in range(NORBITAL):
            one_rows.append("{} {} {} {}\n".format(
                *(fused_to_fields(out_orbital) + fused_to_fields(in_orbital))))
    write(os.path.join(workdir, "greenone.def"),
          "===============================\nNCisAjs {}\n"
          "===============================\n"
          "======== Green functions ======\n"
          "===============================\n{}".format(len(one_rows), "".join(one_rows)))

    two_pairs = []
    for first in range(NORBITAL):
        for second in range(NORBITAL):
            two_pairs.append(((first, first), (second, second)))
            two_pairs.append(((first, second), (second, first)))
    two_rows = []
    for factors in two_pairs:
        fields = []
        for out_orbital, in_orbital in factors:
            fields.extend(fused_to_fields(out_orbital))
            fields.extend(fused_to_fields(in_orbital))
        two_rows.append("{}\n".format(" ".join(str(value) for value in fields)))
    write(os.path.join(workdir, "greentwo.def"),
          "=============================================\nNCisAjsCktAltDC {}\n"
          "=============================================\n"
          "========== Two-body Green ==================\n"
          "=============================================\n{}".format(len(two_rows), "".join(two_rows)))

    up0, up1, down0 = 0, 1, 2
    terms = (
        nbody(((up0, up1), (down0, down0), (up1, up1))),
        adjoint(nbody(((up0, up1), (down0, down0), (up1, up1)))),
        nbody(((up0, up0), (up0, up0), (down0, down0))),
    )
    rows = []
    for operators in terms:
        factors = operators_to_factors(operators)
        rows.append("3 {}\n".format(" ".join(str(value) for value in factors)))
    write(os.path.join(workdir, "nbodyg.def"),
          "=============================================\nNNBodyG {}\n"
          "=============================================\n"
          "============= NBodyG ========================\n"
          "=============================================\n{}".format(len(rows), "".join(rows)))
    return two_pairs, terms


def write_anomalous(workdir, delta):
    term_rows = (
        (1, 0, 0, 1, 1, delta),
        (0, 1, 1, 0, 0, delta.conjugate()),
    )
    green_rows = (
        (1, 0, 0, 1, 1),
        (0, 1, 1, 0, 0),
        (1, 1, 1, 0, 0),
        (0, 0, 0, 1, 1),
    )
    oracle_rows = tuple(
        (operator_type,) + fused_to_fields(first) + fused_to_fields(second)
        for operator_type, first, second, unused_operators
        in ANOMALOUS_G_OPERATORS
    )
    if green_rows != oracle_rows:
        raise AssertionError("AnomalousG fixture and exact oracle order differ")
    write(
        os.path.join(workdir, "anomalousterm.def"),
        "=============================================\n"
        "NAnomalousTerm 2\n"
        "=============================================\n"
        "======== type s1 spin1 s2 spin2 Re Im =======\n"
        "=============================================\n" +
        "".join("{} {} {} {} {} {} {}\n".format(
            row[0], row[1], row[2], row[3], row[4],
            fmt(row[5].real), fmt(row[5].imag)) for row in term_rows),
    )
    write(
        os.path.join(workdir, "anomalousg.def"),
        "=============================================\n"
        "NAnomalousG 4\n"
        "=============================================\n"
        "======== type s1 spin1 s2 spin2 =============\n"
        "=============================================\n" +
        "".join("{} {} {} {} {}\n".format(*row) for row in green_rows),
    )
    return green_rows


def write_fixture(workdir, mode=1, samples=60000, seed=93617, nstore=1,
                  nsplit=1, reweight=0, iterations=2, jastrow=0.23,
                  mu=0.4, vacuum=False, delta=None, init_nelec=None,
                  data_qty=1):
    parameters = tuple(0.0j for unused in range(6)) if vacuum else default_parameters()
    if init_nelec is None:
        init_nelec = 0 if vacuum else 2
    write_modpara(workdir, mode, samples, seed, nstore, nsplit, reweight,
                  iterations, init_nelec, data_qty)
    write_identity_projection(workdir)
    write_orbital(workdir, parameters, optimize=not vacuum)
    namelist = [
        "         ModPara  modpara.def",
        "         LocSpin  locspn.def",
        "  OrbitalGeneral  orbitalidxgen.def",
        "        TransSym  qptransidx.def",
    ]
    if not vacuum:
        write_jastrow(workdir)
        write_hamiltonian(workdir, mu)
        write_measurements(workdir)
        namelist.extend((
            "         Jastrow  jastrowidx.def",
            "           Trans  trans.def",
            "    CoulombIntra  coulombintra.def",
            "    CoulombInter  coulombinter.def",
            "            Hund  hund.def",
            "         PairHop  pairhopp.def",
            "        Exchange  exchange.def",
            "        InterAll  interall.def",
            "  NBodyInterAll  nbodyinterall.def",
            "        OneBodyG  greenone.def",
            "        TwoBodyG  greentwo.def",
            "          NBodyG  nbodyg.def",
        ))
        if delta is not None:
            write_anomalous(workdir, delta)
            namelist.append("   AnomalousTerm  anomalousterm.def")
            if mode == 1:
                namelist.append("      AnomalousG  anomalousg.def")
    write(os.path.join(workdir, "namelist.def"), "\n".join(namelist) + "\n")
    write_initial(workdir, parameters, None if vacuum else jastrow)


def mpi_command(procs, command):
    launcher = os.environ.get("MVMC_MPIEXEC")
    if launcher:
        return launcher.split() + ["-np", str(procs)] + command
    executable = shutil.which("mpiexec") or shutil.which("mpirun")
    if executable is None:
        raise RuntimeError("MPI launcher not found")
    return [executable, "-np", str(procs)] + command


def run_binary(rootdir, workdir, procs=1, extra_env=None):
    binary = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    command = [binary, "-e", "namelist.def", "initial.def"]
    if procs > 1:
        command = mpi_command(procs, command)
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = environment.get("OMP_NUM_THREADS", "1")
    environment["OPENBLAS_NUM_THREADS"] = "1"
    environment["MKL_NUM_THREADS"] = "1"
    if extra_env:
        environment.update(extra_env)
    process = subprocess.run(command, cwd=workdir, env=environment,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True)
    write(os.path.join(workdir, "run.log"), process.stdout)
    if process.returncode != 0:
        raise AssertionError("vmc.out failed:\n{}".format(process.stdout[-8000:]))
    return process.stdout


def parse_complex_rows(path, value_columns):
    rows = []
    with open(path) as stream:
        for line in stream:
            columns = line.split()
            if not columns:
                continue
            rows.append((columns, complex(float(columns[value_columns]),
                                          float(columns[value_columns + 1]))))
    return rows


def state_dump_records(path):
    records = []
    with open(path) as stream:
        for line in stream:
            pieces = line.strip().split("|")
            left = pieces[0].split()
            occupation = tuple(int(value) for value in pieces[2].split())
            mask = sum(value << index for index, value in enumerate(occupation))
            records.append((left[0], int(left[1]), int(left[2]), mask,
                            tuple(piece.strip() for piece in pieces)))
    return records


def assert_close(label, actual, expected, tolerance):
    if abs(actual - expected) > tolerance:
        raise AssertionError("{} mismatch: actual={} expected={} tolerance={}".format(
            label, actual, expected, tolerance))


def validate_burn(records):
    stores = [record for record in records if record[0] == "STORE"]
    restores = [record for record in records if record[0] == "RESTORE"]
    if not restores:
        raise AssertionError("GC Burn restore record is absent")
    if len(stores) < len(restores):
        raise AssertionError("GC Burn store record is absent")
    for index, restore in enumerate(restores):
        # Compare ncur plus active/inactive eleIdx, eleCfg/eleNum and ProjCnt.
        stored = stores[index]
        same = stored[2] == restore[2]
        same = same and stored[4][0].split()[3:] == restore[4][0].split()[3:]
        same = same and stored[4][1:] == restore[4][1:]
        if not same:
            raise AssertionError("GC Burn snapshot/restore byte contract changed")


def validate_distribution(records, jastrow=0.23):
    samples = [record[3] for record in records if record[0] == "SAMPLE"]
    counts = dict((state, 0) for state in EVEN_BASIS)
    for state in samples:
        counts[state] += 1
    wave = wavefunction(jastrow=jastrow)
    denominator = sum(abs(value) ** 2 for value in wave.values())
    chi_square = 0.0
    for state in EVEN_BASIS:
        expected = len(samples) * abs(wave[state]) ** 2 / denominator
        if expected >= 1.0:
            chi_square += (counts[state] - expected) ** 2 / expected
    # The samples are Markov-correlated; this remains a deliberately generous
    # production-chain gate while the independent 32-chain C unit supplies the
    # strict chi-square criterion.
    if chi_square > 120.0:
        raise AssertionError("GC production distribution chi-square={}".format(chi_square))
    if len([count for count in counts.values() if count > 0]) != len(EVEN_BASIS):
        raise AssertionError("GC production chain did not visit all eight even states")
    return len(samples), chi_square


def validate_physical(workdir, samples, jastrow=0.23, mu=0.4):
    expected = exact_observables(jastrow=jastrow, mu=mu)
    with open(os.path.join(workdir, "zvo_gc.dat")) as stream:
        gc = [float(value) for value in stream.readline().split()]
    with open(os.path.join(workdir, "output", "zvo_out_001.dat")) as stream:
        output = [float(value) for value in stream.readline().split()]
    number_tolerance = 10.0 * math.sqrt(expected["variance_number"] / samples) + 0.01
    energy_variance = max(0.0, expected["energy2"].real - abs(expected["energy"]) ** 2)
    energy_tolerance = 10.0 * math.sqrt(energy_variance / samples) + 0.012
    assert_close("<N>", gc[0], expected["number"], number_tolerance)
    assert_close("<N2>", gc[1], expected["number2"], 4.0 * number_tolerance)
    assert_close("var(N)", gc[2], expected["variance_number"], 6.0 * number_tolerance)
    if gc[2] <= 0.0:
        raise AssertionError("production particle-number variance is not positive")
    assert_close("energy real", output[0], expected["energy"].real, energy_tolerance)
    assert_close("energy imag", output[1], expected["energy"].imag, energy_tolerance)
    assert_close("energy squared", output[2], expected["energy2"].real,
                 8.0 * energy_tolerance)

    one_rows = parse_complex_rows(
        os.path.join(workdir, "output", "zvo_cisajs_001.dat"), 4)
    if len(one_rows) != NORBITAL * NORBITAL:
        raise AssertionError("one-body measurement inventory is incomplete")
    for columns, actual in one_rows:
        out_orbital = int(columns[0]) + int(columns[1]) * NSITE
        in_orbital = int(columns[2]) + int(columns[3]) * NSITE
        assert_close("one-body ({},{})".format(out_orbital, in_orbital),
                     actual, expected["greens1"][(out_orbital, in_orbital)], 0.045)

    two_rows = parse_complex_rows(
        os.path.join(workdir, "output", "zvo_cisajscktalt_001.dat"), 8)
    if len(two_rows) != 32:
        raise AssertionError("two-body coincident-index inventory is incomplete")
    wave = expected["wave"]
    for columns, actual in two_rows:
        fused = []
        for offset in (0, 4):
            fused.append((int(columns[offset]) + int(columns[offset + 1]) * NSITE,
                          int(columns[offset + 2]) + int(columns[offset + 3]) * NSITE))
        assert_close("two-body {}".format(fused), actual,
                     expectation(wave, nbody(tuple(fused))), 0.055)

    nbody_rows = parse_complex_rows(
        os.path.join(workdir, "output", "zvo_NBodyG_001.dat"), 13)
    if len(nbody_rows) != 3:
        raise AssertionError("NBodyG n=3 inventory is incomplete")
    for columns, actual in nbody_rows:
        factors = []
        for offset in range(1, 13, 4):
            factors.append((int(columns[offset]) + int(columns[offset + 1]) * NSITE,
                            int(columns[offset + 2]) + int(columns[offset + 3]) * NSITE))
        assert_close("NBodyG {}".format(factors), actual,
                     expectation(wave, nbody(tuple(factors))), 0.055)


def read_nonempty_rows(path):
    with open(path) as stream:
        return [line.split() for line in stream if line.split()]


def read_first_line_bytes(path):
    with open(path, "rb") as stream:
        return stream.readline()


def anomalous_output_path(workdir, data_index=1):
    return os.path.join(workdir, "output",
                        "zvo_anomalousg_{:03d}.dat".format(data_index))


def parse_anomalous_output(workdir, data_index=1):
    path = anomalous_output_path(workdir, data_index)
    rows = parse_complex_rows(path, 5)
    if len(rows) != len(ANOMALOUS_G_OPERATORS):
        raise AssertionError("AnomalousG output inventory is incomplete: {}".format(
            path))
    values = {}
    for columns, value in rows:
        key = (
            int(columns[0]),
            int(columns[1]) + int(columns[2]) * NSITE,
            int(columns[3]) + int(columns[4]) * NSITE,
        )
        if key in values:
            raise AssertionError("duplicate AnomalousG output row: {}".format(key))
        values[key] = value
    return values


def validate_anomalous_file_structure(workdir, data_index):
    path = anomalous_output_path(workdir, data_index)
    with open(path) as stream:
        lines = stream.readlines()
    data_lines = [line for line in lines if line.strip()]
    if len(data_lines) != 4 or not lines or lines[-1].strip():
        raise AssertionError(
            "AnomalousG bin must contain four rows and a trailing blank line: {}"
            .format(path))


def assert_anomalous_files_tightly_equal(first_workdir, second_workdir,
                                          data_index):
    first = parse_complex_rows(
        anomalous_output_path(first_workdir, data_index), 5)
    second = parse_complex_rows(
        anomalous_output_path(second_workdir, data_index), 5)
    if len(first) != len(second):
        raise AssertionError("AnomalousG MPI row counts differ")
    maximum_difference = 0.0
    for (first_columns, first_value), (second_columns, second_value) in zip(
            first, second):
        if first_columns[:5] != second_columns[:5]:
            raise AssertionError("AnomalousG MPI row metadata differ")
        difference = abs(first_value - second_value)
        maximum_difference = max(maximum_difference, difference)
        if difference > 2.0e-14:
            raise AssertionError(
                "AnomalousG MPI reduction differs by {}".format(difference))
    return maximum_difference


def anomalous_tolerances(expected, samples):
    number_tolerance = (
        10.0 * math.sqrt(expected["variance_number"] / samples) + 0.01
    )
    energy_variance = max(
        0.0, expected["energy2"].real - abs(expected["energy"]) ** 2)
    energy_tolerance = 10.0 * math.sqrt(energy_variance / samples) + 0.012
    # All four exact values have magnitude about 0.0845 in this pinned fixture.
    # A 0.025 gate keeps create/remove and operator-order sign mutations visible.
    anomalous_tolerance = 0.025
    return number_tolerance, energy_tolerance, anomalous_tolerance


def validate_anomalous_observables(workdir, samples, delta, data_index=1,
                                   jastrow=0.23, mu=0.4):
    expected = exact_observables(jastrow=jastrow, mu=mu, delta=delta)
    gc_rows = read_nonempty_rows(os.path.join(workdir, "zvo_gc.dat"))
    if data_index > len(gc_rows):
        raise AssertionError("grand-canonical output bin is absent")
    gc = [float(value) for value in gc_rows[data_index - 1]]
    output_rows = read_nonempty_rows(os.path.join(
        workdir, "output", "zvo_out_{:03d}.dat".format(data_index)))
    if len(output_rows) != 1:
        raise AssertionError("physical output must contain exactly one row")
    output = [float(value) for value in output_rows[0]]
    number_tolerance, energy_tolerance, anomalous_tolerance = (
        anomalous_tolerances(expected, samples)
    )
    assert_close("anomalous <N>", gc[0], expected["number"], number_tolerance)
    assert_close("anomalous <N2>", gc[1], expected["number2"],
                 4.0 * number_tolerance)
    assert_close("anomalous var(N)", gc[2], expected["variance_number"],
                 6.0 * number_tolerance)
    assert_close("anomalous energy real", output[0], expected["energy"].real,
                 energy_tolerance)
    assert_close("anomalous energy imag", output[1], expected["energy"].imag,
                 energy_tolerance)
    assert_close("anomalous energy squared", output[2], expected["energy2"].real,
                 8.0 * energy_tolerance)
    actual_anomalous = parse_anomalous_output(workdir, data_index)
    if set(actual_anomalous) != set(expected["anomalous_g"]):
        raise AssertionError("AnomalousG output keys differ from exact oracle")
    for key, actual in actual_anomalous.items():
        assert_close("AnomalousG {}".format(key), actual,
                     expected["anomalous_g"][key], anomalous_tolerance)
    return {
        "expected": expected,
        "gc": gc,
        "output": output,
        "anomalous_g": actual_anomalous,
        "number_tolerance": number_tolerance,
        "energy_tolerance": energy_tolerance,
        "anomalous_tolerance": anomalous_tolerance,
    }


def parse_sr_dump(path):
    steps = []
    current = None
    with open(path) as stream:
        for line in stream:
            columns = line.split()
            if not columns:
                continue
            if columns[0] == "STEP":
                current = {"step": int(columns[1]), "store": int(columns[-1]), "p": {}}
                steps.append(current)
            elif columns[0] == "P":
                current["p"][int(columns[1])] = tuple(float(value) for value in columns[2:])
    return steps


def prepare_work(rootdir, name):
    workdir = os.path.join(rootdir, "work", name)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    return workdir


def physical_case(rootdir):
    samples = int(os.environ.get("MVMC_GC_SAMPLES", "60000"))
    workdir = prepare_work(rootdir, "GC_Exact_Physical")
    write_fixture(workdir, mode=1, samples=samples, iterations=2)
    run_binary(rootdir, workdir,
               extra_env={"MVMC_GC_STATE_DUMP": "gc_state.dat",
                          "MVMC_GC_DEBUG_REBUILD_INTERVAL": "97"})
    records = state_dump_records(os.path.join(workdir, "gc_state.dat.rank0"))
    sample_count, chi_square = validate_distribution(records)
    validate_physical(workdir, sample_count)
    print("GC physical exact oracle passed: samples={} chi2={:.6g}".format(
        sample_count, chi_square))


def sr_case(rootdir):
    runs = {}
    for reweight in (0, 1):
        for nstore in (0, 1):
            workdir = prepare_work(
                rootdir, "GC_SR_Reweight{}_Store{}".format(reweight, nstore))
            write_fixture(workdir, mode=0, samples=50000, seed=44821,
                          nstore=nstore, reweight=reweight, iterations=2)
            run_binary(rootdir, workdir,
                       extra_env={"MVMC_GC_STATE_DUMP": "gc_state.dat",
                                  "MVMC_GC_SR_DUMP": "gc_sr.dat"})
            records = state_dump_records(
                os.path.join(workdir, "gc_state.dat.rank0"))
            validate_burn(records)
            exact = exact_observables(jastrow=0.23, mu=0.4)
            with open(os.path.join(workdir, "zvo_gc.dat")) as stream:
                first_gc = [float(value) for value in stream.readline().split()]
            with open(os.path.join(workdir, "output", "zvo_out_001.dat")) as stream:
                first_energy = [float(value) for value in stream.readline().split()]
            assert_close("SR reweight <N>", first_gc[0], exact["number"], 0.035)
            assert_close("SR reweight <N2>", first_gc[1], exact["number2"], 0.09)
            assert_close("SR reweight energy", first_energy[0],
                         exact["energy"].real, 0.035)
            runs[(reweight, nstore)] = parse_sr_dump(
                os.path.join(workdir, "gc_sr.dat"))
        first_store = runs[(reweight, 0)][0]["p"]
        second_store = runs[(reweight, 1)][0]["p"]
        if set(first_store) != set(second_store):
            raise AssertionError("NStoreO SR parameter maps differ")
        for parameter in first_store:
            assert_close("NStoreO gradient reweight={} parameter={}".format(
                reweight, parameter), first_store[parameter][4],
                second_store[parameter][4], 5.0e-10)

    first = runs[(0, 0)][0]["p"]

    # Parameter 5 real and imaginary use the optimizer's packed
    # mapping pi=2*(NProj+orbital)+component, with NProj=1.
    for orbital, imaginary in ((5, False), (5, True)):
        packed = 2 * (1 + orbital) + (1 if imaginary else 0)
        exact = parameter_gradient(orbital, imaginary, 0.23, 0.4)
        sampled = first[packed][4]
        assert_close("sampled SR covariance {}".format(packed), sampled, exact, 0.012)
        for epsilon in (1.0e-5, 5.0e-6):
            fd = finite_difference(orbital, imaginary, epsilon, 0.23, 0.4)
            assert_close("exact FD {} {}".format(packed, epsilon), fd, exact,
                         1.0e-7 + 1.0e-5 * abs(exact))
    print("GC SR NStoreO/reweight/Burn/gradient oracle passed")


def mpi_case(rootdir, nsplit, mode):
    samples = 24000 if mode == 1 else 20000
    serial = prepare_work(rootdir, "GC_MPI_serial_s{}_m{}".format(nsplit, mode))
    parallel = prepare_work(rootdir, "GC_MPI_parallel_s{}_m{}".format(nsplit, mode))
    for workdir in (serial, parallel):
        write_fixture(workdir, mode=mode, samples=samples, seed=67231,
                      nsplit=nsplit, iterations=2)
    common_env = {"MVMC_GC_STATE_DUMP": "gc_state.dat"}
    run_binary(rootdir, serial, extra_env=common_env)
    run_binary(rootdir, parallel, procs=2, extra_env=common_env)
    serial_records = state_dump_records(os.path.join(serial, "gc_state.dat.rank0"))
    parallel_rank0 = state_dump_records(os.path.join(parallel, "gc_state.dat.rank0"))
    parallel_rank1 = state_dump_records(os.path.join(parallel, "gc_state.dat.rank1"))
    if mode == 0:
        validate_burn(parallel_rank0)
        validate_burn(parallel_rank1)
    if nsplit == 2:
        serial_samples = [record[3] for record in serial_records if record[0] == "SAMPLE"]
        for rank_records in (parallel_rank0, parallel_rank1):
            rank_samples = [record[3] for record in rank_records if record[0] == "SAMPLE"]
            if rank_samples != serial_samples:
                raise AssertionError("NSplitSize=2 saved occupation sequence differs from rank-1")
        with open(os.path.join(serial, "zvo_gc.dat"), "rb") as stream:
            serial_gc = stream.read()
        with open(os.path.join(parallel, "zvo_gc.dat"), "rb") as stream:
            parallel_gc = stream.read()
        if serial_gc != parallel_gc:
            raise AssertionError("NSplitSize=2 normalized GC output is not tightly equal")
    else:
        sample_count, unused_chi = validate_distribution(parallel_rank0)
        sample_count1, unused_chi1 = validate_distribution(parallel_rank1)
        if mode == 1:
            validate_physical(parallel, sample_count + sample_count1)
    print("GC MPI policy passed: NSplitSize={} CalMode={}".format(nsplit, mode))


def vacuum_case(rootdir):
    workdir = prepare_work(rootdir, "GC_Vacuum")
    write_fixture(workdir, mode=1, samples=200, seed=77129, iterations=1,
                  vacuum=True)
    run_binary(rootdir, workdir)
    with open(os.path.join(workdir, "zvo_gc.dat")) as stream:
        values = [float(value) for value in stream.readline().split()]
    if values != [0.0, 0.0, 0.0]:
        raise AssertionError("projector-free vacuum output is not exactly zero: {}".format(values))
    print("GC projector-free vacuum passed")


def mu_response_case(rootdir):
    final_numbers = []
    for mu in (-1.0, 0.0, 1.0):
        workdir = prepare_work(rootdir, "GC_Mu_{:+.0f}".format(mu))
        write_fixture(workdir, mode=0, samples=50000, seed=99123,
                      nstore=1, iterations=4, mu=mu)
        run_binary(rootdir, workdir)
        with open(os.path.join(workdir, "zvo_gc.dat")) as stream:
            rows = [[float(value) for value in line.split()]
                    for line in stream if line.split()]
        if len(rows) != 4:
            raise AssertionError("chemical-potential response iteration count changed")
        final_numbers.append(rows[-1][0])
    if not (final_numbers[0] + 0.01 < final_numbers[1]
            and final_numbers[1] + 0.01 < final_numbers[2]):
        raise AssertionError("<N> does not respond monotonically to mu: {}".format(
            final_numbers))
    print("GC chemical-potential response passed: {}".format(final_numbers))


def anomalous_case(rootdir):
    samples = int(os.environ.get("MVMC_GC_ANOMALOUS_SAMPLES", "60000"))
    delta = 0.35 - 0.20j
    workdir = prepare_work(rootdir, "GC_Anomalous_Exact")
    write_fixture(workdir, mode=1, samples=samples, seed=38611,
                  iterations=2, delta=delta)
    run_binary(rootdir, workdir,
               extra_env={"MVMC_GC_STATE_DUMP": "gc_state.dat",
                          "MVMC_GC_DEBUG_REBUILD_INTERVAL": "97"})
    records = state_dump_records(os.path.join(workdir, "gc_state.dat.rank0"))
    sample_count, chi_square = validate_distribution(records)
    result = validate_anomalous_observables(
        workdir, sample_count, delta)

    # The absolute-energy tolerance includes conservative Markov-chain error.
    # A separate correlated-difference gate pins the smaller isolated pairing
    # contribution without weakening the absolute observable comparison.
    response_tolerance = 0.003
    reference = exact_observables(jastrow=0.23, mu=0.4, delta=0.0)
    exact_response = result["expected"]["energy"] - reference["energy"]
    if abs(exact_response) <= 3.0 * response_tolerance:
        raise AssertionError("anomalous energy-response fixture is vacuous")
    for operator_type in (0, 1):
        magnitudes = [
            abs(value) for key, value in result["expected"]["anomalous_g"].items()
            if key[0] == operator_type
        ]
        if not magnitudes or max(magnitudes) <= 3.0 * result["anomalous_tolerance"]:
            raise AssertionError(
                "AnomalousG type {} fixture is vacuous".format(operator_type))
    print("GC anomalous exact oracle passed: samples={} chi2={:.6g} "
          "dE={:.8g} AGtol={:.4g}".format(
              sample_count, chi_square, exact_response.real,
              result["anomalous_tolerance"]))


def anomalous_energy_case(rootdir):
    samples = int(os.environ.get("MVMC_GC_ANOMALOUS_SAMPLES", "60000"))
    deltas = (0.0, 0.2, 0.4 * complex(math.cos(math.pi / 3.0),
                                      math.sin(math.pi / 3.0)))
    results = []
    gc_bytes = []
    anomalous_bytes = []
    for index, delta in enumerate(deltas):
        workdir = prepare_work(rootdir, "GC_Anomalous_Energy_{}".format(index))
        write_fixture(workdir, mode=1, samples=samples, seed=55291,
                      iterations=2, delta=delta)
        run_binary(rootdir, workdir)
        result = validate_anomalous_observables(workdir, samples, delta)
        results.append(result)
        with open(os.path.join(workdir, "zvo_gc.dat"), "rb") as stream:
            gc_bytes.append(stream.read())
        with open(anomalous_output_path(workdir), "rb") as stream:
            anomalous_bytes.append(stream.read())

    if len(set(gc_bytes)) != 1:
        raise AssertionError("fixed-f GC observables depend on anomalous delta")
    if len(set(anomalous_bytes)) != 1:
        raise AssertionError("fixed-f AnomalousG observables depend on delta")

    response_tolerance = 0.003
    actual_zero = complex(results[0]["output"][0], results[0]["output"][1])
    expected_zero = results[0]["expected"]["energy"]
    for delta, result in zip(deltas[1:], results[1:]):
        actual = complex(result["output"][0], result["output"][1]) - actual_zero
        expected = result["expected"]["energy"] - expected_zero
        if abs(expected) <= 3.0 * response_tolerance:
            raise AssertionError("anomalous delta response is vacuous: {}".format(delta))
        assert_close("isolated anomalous energy response {}".format(delta),
                     actual, expected, response_tolerance)
    print("GC anomalous energy response passed: {}".format([
        result["output"][0] for result in results]))


def sr_parameter_line(path, packed):
    prefix = "P {} ".format(packed)
    with open(path, "rb") as stream:
        for line in stream:
            if line.decode("ascii").startswith(prefix):
                return line
    raise AssertionError("SR parameter {} is absent".format(packed))


def anomalous_sr_case(rootdir):
    delta = 0.35 - 0.20j
    samples = int(os.environ.get("MVMC_GC_ANOMALOUS_SR_SAMPLES", "50000"))
    selected = ((5, False), (3, True))
    exact_gradients = {}
    for orbital in range(6):
        for imaginary in (False, True):
            packed = 2 * (1 + orbital) + (1 if imaginary else 0)
            gradient = parameter_gradient(
                orbital, imaginary, 0.23, 0.4, delta)
            fd_values = tuple(
                finite_difference(orbital, imaginary, epsilon,
                                  0.23, 0.4, delta)
                for epsilon in (1.0e-5, 5.0e-6)
            )
            exact_gradients[packed] = gradient
            for epsilon, fd in zip((1.0e-5, 5.0e-6), fd_values):
                assert_close("anomalous exact FD {} {}".format(packed, epsilon),
                             fd, gradient,
                             1.0e-7 + 1.0e-5 * abs(gradient))
            print("anomalous SR oracle P={} g={:.12g} fd1={:.12g} "
                  "fd2={:.12g}".format(
                      packed, gradient, fd_values[0], fd_values[1]))

    runs = {}
    paths = {}
    for nstore in (0, 1):
        workdir = prepare_work(rootdir, "GC_Anomalous_SR_Store{}".format(nstore))
        write_fixture(workdir, mode=0, samples=samples, seed=44821,
                      nstore=nstore, iterations=2, delta=delta)
        run_binary(rootdir, workdir,
                   extra_env={"MVMC_GC_SR_DUMP": "gc_sr.dat"})
        paths[nstore] = workdir
        runs[nstore] = parse_sr_dump(os.path.join(workdir, "gc_sr.dat"))
        if not runs[nstore]:
            raise AssertionError("anomalous SR dump is empty")

    for orbital, imaginary in selected:
        packed = 2 * (1 + orbital) + (1 if imaginary else 0)
        exact = exact_gradients[packed]
        if abs(exact) < 0.04:
            raise AssertionError("selected anomalous SR gradient is vacuous")
        tolerance = min(0.025, 0.5 * abs(exact))
        se_budget = tolerance / 4.0
        if tolerance < 4.0 * se_budget or tolerance > 0.5 * abs(exact):
            raise AssertionError("anomalous SR detection threshold is inconsistent")
        for nstore in (0, 1):
            sampled = runs[nstore][0]["p"][packed][4]
            assert_close("anomalous sampled SR gradient P={}".format(packed),
                         sampled, exact, tolerance)
        first_line = sr_parameter_line(
            os.path.join(paths[0], "gc_sr.dat"), packed)
        second_line = sr_parameter_line(
            os.path.join(paths[1], "gc_sr.dat"), packed)
        if first_line != second_line:
            raise AssertionError(
                "NStoreO changed selected anomalous SR gradient P={}".format(packed))
        print("anomalous SR selected P={} |g|={:.8g} SE_budget={:.4g} "
              "tolerance={:.4g}".format(
                  packed, abs(exact), se_budget, tolerance))

    gc0 = read_first_line_bytes(os.path.join(paths[0], "zvo_gc.dat"))
    gc1 = read_first_line_bytes(os.path.join(paths[1], "zvo_gc.dat"))
    if gc0 != gc1:
        raise AssertionError("NStoreO changed anomalous SR GC output")
    print("GC anomalous SR NStoreO/exact-gradient oracle passed")


def anomalous_mpi_case(rootdir):
    samples = int(os.environ.get("MVMC_GC_ANOMALOUS_MPI_SAMPLES", "24000"))
    delta = 0.35 - 0.20j
    serial = prepare_work(rootdir, "GC_Anomalous_MPI_serial")
    parallel = prepare_work(rootdir, "GC_Anomalous_MPI_parallel")
    for workdir in (serial, parallel):
        write_fixture(workdir, mode=1, samples=samples, seed=67231,
                      nsplit=2, iterations=2, delta=delta, data_qty=2)
    common_env = {"MVMC_GC_STATE_DUMP": "gc_state.dat"}
    run_binary(rootdir, serial, extra_env=common_env)
    run_binary(rootdir, parallel, procs=2, extra_env=common_env)
    serial_records = state_dump_records(os.path.join(serial, "gc_state.dat.rank0"))
    serial_samples = [record[3] for record in serial_records
                      if record[0] == "SAMPLE"]
    for rank in (0, 1):
        rank_records = state_dump_records(os.path.join(
            parallel, "gc_state.dat.rank{}".format(rank)))
        rank_samples = [record[3] for record in rank_records
                        if record[0] == "SAMPLE"]
        if rank_samples != serial_samples:
            raise AssertionError(
                "anomalous NSplitSize=2 occupation sequence differs on rank {}"
                .format(rank))
    with open(os.path.join(serial, "zvo_gc.dat"), "rb") as stream:
        serial_gc = stream.read()
    with open(os.path.join(parallel, "zvo_gc.dat"), "rb") as stream:
        parallel_gc = stream.read()
    if serial_gc != parallel_gc:
        raise AssertionError("anomalous NSplitSize=2 GC output differs")
    for data_index in (1, 2):
        validate_anomalous_file_structure(serial, data_index)
        validate_anomalous_file_structure(parallel, data_index)
        # The two ranks consume disjoint halves of an identical saved chain.
        # MPI therefore changes floating-point summation grouping, so the
        # 18-digit text is not generally byte-identical.  Pin the five input
        # columns exactly and the complex values to near-roundoff instead.
        maximum_difference = assert_anomalous_files_tightly_equal(
            serial, parallel, data_index)
        print("anomalous MPI bin={} maxdiff={:.3g}".format(
            data_index, maximum_difference))
    print("GC anomalous MPI NSplitSize=2 and two-bin output passed")


def anomalous_vacuum_case(rootdir):
    samples = int(os.environ.get("MVMC_GC_ANOMALOUS_SAMPLES", "60000"))
    delta = 0.35 - 0.20j
    workdir = prepare_work(rootdir, "GC_Anomalous_Vacuum_Smoke")
    # This is an initial-state smoke test.  The variational wave function is
    # unchanged; only the starting occupation is the vacuum.
    write_fixture(workdir, mode=1, samples=samples, seed=77129,
                  iterations=1, delta=delta, init_nelec=0)
    run_binary(rootdir, workdir,
               extra_env={"MVMC_GC_STATE_DUMP": "gc_state.dat"})
    records = state_dump_records(os.path.join(workdir, "gc_state.dat.rank0"))
    sample_count, chi_square = validate_distribution(records)
    validate_anomalous_observables(workdir, sample_count, delta)
    print("GC anomalous initial-vacuum smoke passed: samples={} chi2={:.6g}"
          .format(sample_count, chi_square))


CASES = {
    "anomalous": anomalous_case,
    "anomalous_energy": anomalous_energy_case,
    "anomalous_mpi": anomalous_mpi_case,
    "anomalous_sr": anomalous_sr_case,
    "anomalous_vacuum": anomalous_vacuum_case,
    "physical": physical_case,
    "sr": sr_case,
    "mpi_nsplit2_mode0": lambda root: mpi_case(root, 2, 0),
    "mpi_nsplit2_mode1": lambda root: mpi_case(root, 2, 1),
    "mpi_nsplit1_mode0": lambda root: mpi_case(root, 1, 0),
    "mpi_nsplit1_mode1": lambda root: mpi_case(root, 1, 1),
    "vacuum": vacuum_case,
    "mu_response": mu_response_case,
}


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in CASES:
        print("usage: {} {}".format(sys.argv[0], "|".join(sorted(CASES))))
        return 2
    rootdir = os.getcwd()
    CASES[sys.argv[1]](rootdir)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print("ERROR: {}".format(error))
        sys.exit(1)
