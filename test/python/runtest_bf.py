from __future__ import print_function

import math
import os
import shutil
import subprocess
import sys

from backflow_def_helper import write_chain_nn_backflow


def parse_nsite(modpara_path):
    with open(modpara_path) as fp:
        for line in fp:
            cols = line.split()
            if len(cols) >= 2 and cols[0] == "Nsite":
                return int(cols[1])
    raise RuntimeError("Nsite was not found in {}".format(modpara_path))


def read_key_value_file(path):
    values = {}
    with open(path) as fp:
        for line in fp:
            cols = line.split()
            if len(cols) == 2:
                values[cols[0]] = cols[1]
    return values


def first_row_is_finite(path):
    with open(path) as fp:
        for line in fp:
            cols = line.split()
            if not cols:
                continue
            return all(math.isfinite(float(x)) for x in cols)
    return False


def main():
    if len(sys.argv) < 2:
        print("usage: {} <model name> [--expect-error <substring>] [--expect-nqp-full <n>]".format(sys.argv[0]))
        return -1

    model = sys.argv[1]
    expected_error = None
    expected_nqp_full = None
    argi = 2
    while argi < len(sys.argv):
        if sys.argv[argi] == "--expect-error" and argi + 1 < len(sys.argv):
            expected_error = sys.argv[argi + 1]
            argi += 2
        elif sys.argv[argi] == "--expect-nqp-full" and argi + 1 < len(sys.argv):
            expected_nqp_full = int(sys.argv[argi + 1])
            argi += 2
        else:
            print("usage: {} <model name> [--expect-error <substring>] [--expect-nqp-full <n>]".format(sys.argv[0]))
            return -1
    rootdir = os.getcwd()
    refdir = os.path.join(rootdir, "data", model)
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    workdir = os.path.join(rootdir, "work", model)
    if mpi_procs:
        workdir += "_mpi{}".format(mpi_procs)

    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)

    for fn in os.listdir(refdir):
        if fn.endswith(".def"):
            shutil.copy(os.path.join(refdir, fn), os.path.join(workdir, fn))

    nsite = parse_nsite(os.path.join(workdir, "modpara.def"))
    write_chain_nn_backflow(workdir, length=nsite, optimize=False)

    os.chdir(workdir)
    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    dump_path = os.path.join(workdir, "bf_identity_dump.dat")
    env = os.environ.copy()
    env["MVMC_BF_IDENTITY_DUMP"] = dump_path

    cmd = [bin_to_test, "-e", "namelist.def"]
    if mpi_procs:
        cmd = ["mpirun", "-np", mpi_procs] + cmd
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env,
    )
    with open("bf_test.log", "w") as fp:
        fp.write(proc.stdout)

    if expected_error is not None:
        if expected_error not in proc.stdout:
            print("ERROR: expected error substring not found: {}".format(expected_error))
            print("---- output begin ----")
            print(proc.stdout)
            print("---- output end ----")
            return -1
        if proc.returncode == 0:
            print("ERROR: vmc.out unexpectedly succeeded for invalid BackFlow input.")
            return -1
        return 0

    if proc.returncode != 0:
        print(proc.stdout)
        return proc.returncode
    if not os.path.exists(dump_path):
        print("ERROR: BackFlow identity dump was not written.")
        return -1

    values = read_key_value_file(dump_path)
    info_no_bf = int(values.get("info_no_bf", "-1"))
    info_bf = int(values.get("info_bf", "-1"))
    nqp_full = int(values.get("nqp_full", "-1"))
    nan_count = int(values.get("nan_count", "-1"))
    max_slater = float(values.get("max_abs_slater_diff", "nan"))
    max_pf = float(values.get("max_abs_pf_diff", "nan"))
    tol = 1.0e-10

    if expected_nqp_full is not None and nqp_full != expected_nqp_full:
        print("ERROR: unexpected nqp_full: got {} expected {}".format(nqp_full, expected_nqp_full))
        return -1
    if info_no_bf != 0 or info_bf != 0:
        print("ERROR: identity dump matrix calculation failed: no_bf={} bf={}".format(info_no_bf, info_bf))
        return -1
    if nan_count != 0:
        print("ERROR: identity dump contains non-finite values: nan_count={}".format(nan_count))
        return -1
    if not math.isfinite(max_slater) or max_slater > tol:
        print("ERROR: SlaterElm identity mismatch: max_abs_slater_diff={}".format(max_slater))
        return -1
    if not math.isfinite(max_pf) or max_pf > tol:
        print("ERROR: PfM identity mismatch: max_abs_pf_diff={}".format(max_pf))
        return -1
    if not first_row_is_finite(os.path.join("output", "zvo_out_001.dat")):
        print("ERROR: zvo_out_001.dat is missing or non-finite.")
        return -1

    return 0


if __name__ == "__main__":
    sys.exit(main())
