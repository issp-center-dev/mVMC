from __future__ import print_function

import os
import shutil
import subprocess
import sys

from backflow_def_helper import write_chain_nn_backflow


MODEL = "HubbardChain_SpinJastrow"


def update_key(path, key, value):
    lines = []
    found = False
    with open(path) as source:
        for line in source:
            words = line.split()
            if words and words[0] == key:
                lines.append("{:<15} {}\n".format(key, value))
                found = True
            else:
                lines.append(line)
    if not found:
        lines.append("{:<15} {}\n".format(key, value))
    with open(path, "w") as destination:
        destination.writelines(lines)


def copy_definitions(rootdir, suite, name):
    source = os.path.join(rootdir, "data", MODEL)
    destination = os.path.join(rootdir, "work", suite, name)
    if os.path.exists(destination):
        shutil.rmtree(destination)
    os.makedirs(destination)
    for filename in os.listdir(source):
        if filename.endswith(".def"):
            shutil.copy(os.path.join(source, filename), destination)
    return destination


def projection_header(keyword, count):
    return [
        "=============================================\n",
        "{:<18} {}\n".format(keyword, count),
        "=============================================\n",
        "======== TrIdx_TrWeight_and_TrIdx_i_xi ======\n",
        "=============================================\n",
    ]


def identity_rows(nsite, transform=0, sign=1):
    return ["{} {} {} {}\n".format(transform, site, site, sign)
            for site in range(nsite)]


def identity_rows_legacy(nsite, transform=0):
    return ["{} {} {}\n".format(transform, site, site)
            for site in range(nsite)]


def shifted_rows(nsite, transform, shift, sign=1):
    return ["{} {} {} {}\n".format(
        transform, site, (site + shift) % nsite, sign)
        for site in range(nsite)]


def write_transsym(path, count, parameters, rows):
    with open(path, "w") as destination:
        destination.writelines(projection_header("NQPTrans", count))
        destination.writelines([line.rstrip("\n") + "\n"
                                for line in parameters])
        destination.writelines([line.rstrip("\n") + "\n"
                                for line in rows])


def write_opttrans(path, count, parameters, rows):
    with open(path, "w") as destination:
        destination.writelines(projection_header("NQPOptTrans", count))
        destination.writelines([line.rstrip("\n") + "\n"
                                for line in parameters])
        destination.writelines([line.rstrip("\n") + "\n"
                                for line in rows])


def add_opttrans_to_namelist(path):
    with open(path, "a") as destination:
        destination.write("        OptTrans  qpopttrans.def\n")


def executable_command(rootdir, workdir, mpi, opttrans):
    binary = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    command = [binary]
    if opttrans:
        command.append("-o")
    command.extend(["-e", "namelist.def"])
    if mpi:
        launcher = os.environ.get("MVMC_MPIEXEC", "mpirun")
        count_flag = os.environ.get("MVMC_MPIEXEC_NUMPROC_FLAG", "-np")
        command = [launcher, count_flag, "2"] + command
    return command


def run_case(rootdir, suite, name, mutate, expected=None, mpi=False,
             opttrans=False):
    workdir = copy_definitions(rootdir, suite, name)
    mutate(workdir)
    process = subprocess.run(
        executable_command(rootdir, workdir, mpi, opttrans),
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=60,
    )
    with open(os.path.join(workdir, "projection_parser.log"), "w") as log:
        log.write(process.stdout)
    if expected is None:
        if process.returncode != 0:
            raise RuntimeError(
                "{} unexpectedly failed (return code {}):\n{}".format(
                    name, process.returncode, process.stdout))
    else:
        if process.returncode == 0:
            raise RuntimeError("{} unexpectedly succeeded".format(name))
        if expected not in process.stdout:
            raise RuntimeError(
                "{} did not report {!r}:\n{}".format(
                    name, expected, process.stdout))
    return workdir


def set_nmp(value):
    def mutate(workdir):
        update_key(os.path.join(workdir, "modpara.def"), "NMPTrans", value)
    return mutate


def replace_transsym(count, parameters, rows, nmp="1"):
    def mutate(workdir):
        update_key(os.path.join(workdir, "modpara.def"), "NMPTrans", nmp)
        write_transsym(os.path.join(workdir, "qptransidx.def"), count,
                       parameters, rows)
    return mutate


def replace_transsym_count(count, nmp="1"):
    def mutate(workdir):
        update_key(os.path.join(workdir, "modpara.def"), "NMPTrans", nmp)
        path = os.path.join(workdir, "qptransidx.def")
        lines = []
        with open(path) as source:
            for line in source:
                if line.split()[:1] == ["NQPTrans"]:
                    lines.append("NQPTrans          {}\n".format(count))
                else:
                    lines.append(line)
        with open(path, "w") as destination:
            destination.writelines(lines)
    return mutate


def invalid_serial_cases(rootdir, suite):
    nsite = 6
    identity = identity_rows(nsite)
    cases = [
        ("nmp_zero", set_nmp("0"), "NMPTrans must be a nonzero integer"),
        ("nmp_int_min", set_nmp("-2147483648"),
         "NMPTrans must be a nonzero integer"),
        ("nmp_fractional", set_nmp("1.5"), "NMPTrans must be an exact integer"),
        ("nmp_too_large", set_nmp("3"), "abs(NMPTrans) <= NQPTrans"),
        ("nmp_negative_too_large", set_nmp("-3"), "abs(NMPTrans) <= NQPTrans"),
        ("nqp_zero", replace_transsym_count("0"),
         "NQPTrans must be an integer in [1"),
        ("nqp_negative", replace_transsym_count("-1"),
         "NQPTrans must be an integer in [1"),
        ("nqp_int_overflow", replace_transsym_count("2147483648"),
         "NQPTrans must be an integer in [1"),
        ("allocation_overflow", replace_transsym_count("400000000"),
         "definition table allocation size exceeds"),
        ("parameter_malformed",
         replace_transsym(1, ["0 1.0 0.0 trailing"], identity),
         "malformed projection parameter row"),
        ("parameter_nonfinite",
         replace_transsym(1, ["0 nan 0.0"], identity),
         "malformed projection parameter row"),
        ("parameter_out_of_range",
         replace_transsym(1, ["1 1.0 0.0"], identity),
         "projection parameter index out of range"),
        ("parameter_duplicate",
         replace_transsym(2, ["0 1.0", "0 1.0"],
                          identity_rows(nsite, 0) + shifted_rows(nsite, 1, 1)),
         "duplicate projection parameter index"),
        ("parameter_missing",
         replace_transsym(2, ["0 1.0"], []),
         "missing projection parameter row"),
        ("mapping_malformed",
         replace_transsym(1, ["0 1.0"], ["0 0"] + identity[1:]),
         "malformed projection mapping row"),
        ("mapping_ap_missing_sign",
         replace_transsym(1, ["0 1.0"], identity_rows_legacy(nsite),
                          nmp="-1"),
         "malformed projection mapping row"),
        ("mapping_transform_range",
         replace_transsym(1, ["0 1.0"], ["1 0 0 1"] + identity[1:]),
         "projection mapping index out of range"),
        ("mapping_source_range",
         replace_transsym(1, ["0 1.0"], ["0 6 0 1"] + identity[1:]),
         "projection mapping index out of range"),
        ("mapping_destination_range",
         replace_transsym(1, ["0 1.0"], ["0 0 6 1"] + identity[1:]),
         "projection mapping index out of range"),
        ("mapping_duplicate_source",
         replace_transsym(1, ["0 1.0"],
                          [identity[0], "0 0 1 1"] + identity[2:]),
         "duplicate projection source"),
        ("mapping_not_permutation",
         replace_transsym(1, ["0 1.0"],
                          [identity[0], "0 1 0 1"] + identity[2:]),
         "projection mapping is not a permutation"),
        ("mapping_bad_sign",
         replace_transsym(1, ["0 1.0"], ["0 0 0 0"] + identity[1:]),
         "projection sign must be +1 or -1"),
        ("mapping_missing",
         replace_transsym(1, ["0 1.0"], identity[:-1]),
         "projection mapping row count mismatch"),
        ("mapping_extra",
         replace_transsym(1, ["0 1.0"], identity + ["0 0 0 1"]),
         "extra projection mapping row"),
    ]
    for name, mutate, expected in cases:
        run_case(rootdir, suite, name, mutate, expected)

    run_case(
        rootdir, suite, "mapping_pbc_legacy_three_fields",
        replace_transsym(1, ["0 1.0"], identity_rows_legacy(nsite)))

    def malformed_opt_parameter(workdir):
        add_opttrans_to_namelist(os.path.join(workdir, "namelist.def"))
        write_opttrans(os.path.join(workdir, "qpopttrans.def"), 1,
                       ["0 1.0 trailing"], identity)

    def malformed_opt_mapping(workdir):
        add_opttrans_to_namelist(os.path.join(workdir, "namelist.def"))
        write_opttrans(os.path.join(workdir, "qpopttrans.def"), 1,
                       ["0 1.0"], ["0 0 0 2"] + identity[1:])

    def invalid_opt_count(workdir):
        add_opttrans_to_namelist(os.path.join(workdir, "namelist.def"))
        write_opttrans(os.path.join(workdir, "qpopttrans.def"), 0, [], [])

    run_case(rootdir, suite, "opt_parameter_malformed",
             malformed_opt_parameter, "malformed OptTrans parameter row",
             opttrans=True)
    run_case(rootdir, suite, "opt_mapping_bad_sign",
             malformed_opt_mapping, "projection sign must be +1 or -1",
             opttrans=True)
    run_case(rootdir, suite, "opt_count_zero",
             invalid_opt_count, "NQPOptTrans must be an integer in [1",
             opttrans=True)

    def add_backflow(workdir):
        update_key(os.path.join(workdir, "modpara.def"), "NMPTrans", "1")
        update_key(os.path.join(workdir, "modpara.def"), "NSPGaussLeg", "1")
        with open(os.path.join(workdir, "namelist.def"), "a") as namelist:
            namelist.write("              BF  bf.def\n")
            namelist.write("         BFRange  rangebf.def\n")
        write_chain_nn_backflow(workdir, length=nsite, optimize=False)

    def legacy_nonidentity(workdir):
        add_backflow(workdir)
        write_transsym(
            os.path.join(workdir, "qptransidx.def"), 2,
            ["0 1.0", "1 1.0"],
            shifted_rows(nsite, 0, 1) + shifted_rows(nsite, 1, 2))

    def nonfsz_opttrans(workdir):
        add_backflow(workdir)
        add_opttrans_to_namelist(os.path.join(workdir, "namelist.def"))
        write_opttrans(os.path.join(workdir, "qpopttrans.def"), 1,
                       ["0 1.0"], identity)

    run_case(rootdir, suite, "nonfsz_nonidentity_first_transform",
             legacy_nonidentity)
    run_case(rootdir, suite, "nonfsz_opttrans_guard",
             nonfsz_opttrans, "non-FSZ BackFlow does not support OptTrans",
             opttrans=True)


def unused_transform_invariance(rootdir, suite):
    nsite = 6

    def variant(shift):
        def mutate(workdir):
            for key, value in (
                    ("NMPTrans", "1"),
                    ("NSPGaussLeg", "1"),
                    ("NSROptItrStep", "1"),
                    ("NSROptItrSmp", "1"),
                    ("NVMCWarmUp", "2"),
                    ("NVMCSample", "8")):
                update_key(os.path.join(workdir, "modpara.def"), key, value)
            write_transsym(
                os.path.join(workdir, "qptransidx.def"), 2,
                ["0 1.0", "1 1.0"],
                identity_rows(nsite, 0) + shifted_rows(nsite, 1, shift))
        return mutate

    first = run_case(rootdir, suite, "unused_shift_1", variant(1))
    second = run_case(rootdir, suite, "unused_shift_2", variant(2))
    outputs = ("zqp_opt.dat", "zvo_out_001.dat", "zvo_var_001.dat")
    for filename in outputs:
        first_path = os.path.join(first, "output", filename)
        second_path = os.path.join(second, "output", filename)
        with open(first_path, "rb") as source:
            first_bytes = source.read()
        with open(second_path, "rb") as source:
            second_bytes = source.read()
        if first_bytes != second_bytes:
            raise RuntimeError(
                "unused QPTrans row changed {}".format(filename))


def mpi_cases(rootdir, suite):
    run_case(rootdir, suite, "cross_field_rank_safe", set_nmp("3"),
             "abs(NMPTrans) <= NQPTrans", mpi=True)
    run_case(
        rootdir, suite, "mapping_rank_safe",
        replace_transsym(1, ["0 1.0"],
                         ["0 0 0 0"] + identity_rows(6)[1:]),
        "projection sign must be +1 or -1", mpi=True)


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in ("serial", "mpi"):
        print("usage: {} serial|mpi".format(sys.argv[0]))
        return 2
    rootdir = os.getcwd()
    suite = "ProjectionParserContract_{}_{}".format(sys.argv[1], os.getpid())
    try:
        if sys.argv[1] == "serial":
            invalid_serial_cases(rootdir, suite)
            unused_transform_invariance(rootdir, suite)
        else:
            mpi_cases(rootdir, suite)
    except (OSError, RuntimeError, subprocess.TimeoutExpired) as error:
        print("ERROR: {}".format(error))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
