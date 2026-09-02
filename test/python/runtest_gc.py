from __future__ import print_function

import math
import os
import shutil
import subprocess
import sys

from grandcanonical_exact_oracle import (
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
                  nsplit=1, reweight=0, iterations=2, init_nelec=2):
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
NDataQtySmp    1
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


def write_fixture(workdir, mode=1, samples=60000, seed=93617, nstore=1,
                  nsplit=1, reweight=0, iterations=2, jastrow=0.23,
                  mu=0.4, vacuum=False):
    parameters = tuple(0.0j for unused in range(6)) if vacuum else default_parameters()
    write_modpara(workdir, mode, samples, seed, nstore, nsplit, reweight,
                  iterations, 0 if vacuum else 2)
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


CASES = {
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
