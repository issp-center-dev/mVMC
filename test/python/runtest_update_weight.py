from __future__ import print_function

import os
import re
import shutil
import subprocess
import sys


BASE_MODEL = "NBodyG_fsz_Compare"


def replace_setting(path, name, value):
    with open(path) as source:
        text = source.read()
    pattern = re.compile(r"^{}\s+.*$".format(re.escape(name)), re.MULTILINE)
    updated, count = pattern.subn("{:<15s} {}".format(name, value), text)
    if count != 1:
        raise AssertionError("{}: expected one {} setting".format(path, name))
    with open(path, "w") as destination:
        destination.write(updated)


def prepare_case(rootdir, name, pair_weight=None, two_sz=-1, update_path=2):
    suffix = os.environ.get("MVMC_MPI_PROCS")
    work_name = "UpdateWeight_{}{}".format(
        name, "_mpi" if suffix else ""
    )
    workdir = os.path.join(rootdir, "work", work_name)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    source = os.path.join(rootdir, "data", BASE_MODEL)
    for filename in os.listdir(source):
        if filename.endswith(".def"):
            shutil.copy2(
                os.path.join(source, filename),
                os.path.join(workdir, filename),
            )

    modpara = os.path.join(workdir, "modpara.def")
    replace_setting(modpara, "Ncond", "0")
    replace_setting(modpara, "2Sz", str(two_sz))
    replace_setting(modpara, "NExUpdatePath", str(update_path))

    locspin = os.path.join(workdir, "locspn.def")
    with open(locspin) as source_file:
        text = source_file.read()
    text = re.sub(r"^NlocalSpin\s+0$", "NlocalSpin     6", text,
                  flags=re.MULTILINE)
    text = re.sub(r"^(\s*\d+\s+)0$", r"\g<1>1", text,
                  flags=re.MULTILINE)
    with open(locspin, "w") as destination:
        destination.write(text)

    if pair_weight is not None:
        with open(os.path.join(workdir, "updateweight.def"), "w") as output:
            output.write(
                "================================\n"
                "NUpdateWeight 3\n"
                "================================\n"
                "========UpdateType Weight=======\n"
                "================================\n"
                "Exchange 1.0\n"
                "LocalSpinFlip 1.0\n"
                "PairSpinFlip {:.17g}\n".format(pair_weight)
            )
        with open(os.path.join(workdir, "namelist.def"), "a") as namelist:
            namelist.write("  InUpdateWeight  updateweight.def\n")
    return workdir


def run_vmc(rootdir, workdir, expect_success=True, diagnostic=None):
    binary = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    command = [binary, "-e", "namelist.def", "initial.def"]
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    if mpi_procs:
        command = ["mpirun", "-np", mpi_procs] + command
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    process = subprocess.run(
        command,
        cwd=workdir,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )
    if expect_success and process.returncode != 0:
        raise AssertionError("UpdateWeight run failed:\n{}".format(process.stdout))
    if not expect_success and process.returncode == 0:
        raise AssertionError("invalid UpdateWeight run unexpectedly succeeded")
    if diagnostic is not None and diagnostic not in process.stdout:
        raise AssertionError(
            "missing diagnostic {!r}:\n{}".format(diagnostic, process.stdout)
        )
    return process


def read_bytes(path):
    with open(path, "rb") as source:
        return source.read()


def time_lines(workdir):
    path = os.path.join(workdir, "output", "zvo_time_001.dat")
    with open(path) as source:
        return [line.rstrip("\n") for line in source]


def main():
    rootdir = os.getcwd()
    legacy = prepare_case(rootdir, "Legacy")
    explicit_legacy = prepare_case(rootdir, "ExplicitLegacy", pair_weight=0.0)
    pair = prepare_case(rootdir, "Pair", pair_weight=1.0)

    run_vmc(rootdir, legacy)
    run_vmc(rootdir, explicit_legacy)
    run_vmc(rootdir, pair)

    for filename in (
        "zvo_NBodyG_001.dat", "zvo_out_001.dat", "zvo_var_001.dat"
    ):
        legacy_bytes = read_bytes(os.path.join(legacy, "output", filename))
        explicit_bytes = read_bytes(
            os.path.join(explicit_legacy, "output", filename)
        )
        if legacy_bytes != explicit_bytes:
            raise AssertionError(
                "legacy and explicit 1:1:0 outputs differ: {}".format(filename)
            )

    legacy_time = time_lines(legacy)
    explicit_time = time_lines(explicit_legacy)
    pair_time = time_lines(pair)
    if "acc_psf" in legacy_time[0]:
        raise AssertionError("legacy zvo_time format changed")
    if "acc_psf" not in explicit_time[0]:
        raise AssertionError("explicit zvo_time lacks pair counters")
    proposals = sum(
        int(line.split(":", 1)[0].split()[8]) for line in pair_time[1:]
    )
    if proposals <= 0:
        raise AssertionError("PairSpinFlip produced no proposals")

    invalid_path = prepare_case(
        rootdir, "InvalidPath", pair_weight=1.0, update_path=1
    )
    run_vmc(
        rootdir, invalid_path, expect_success=False,
        diagnostic="InUpdateWeight currently requires NExUpdatePath=2",
    )
    invalid_sz = prepare_case(
        rootdir, "InvalidSz", pair_weight=1.0, two_sz=0
    )
    run_vmc(
        rootdir, invalid_sz, expect_success=False,
        diagnostic="LocalSpinFlip with NExUpdatePath=2 requires 2Sz=-1",
    )

    print(
        "UpdateWeight end-to-end{}: PASS (pair proposals={})".format(
            " MPI" if os.environ.get("MVMC_MPI_PROCS") else "", proposals
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
