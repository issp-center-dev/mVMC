from __future__ import print_function

import math
import os
import shutil
import subprocess
import sys


TOL = 1.0e-10


def copy_def_files(srcdir, dstdir):
    for name in os.listdir(srcdir):
        if name.endswith(".def"):
            shutil.copy(os.path.join(srcdir, name), os.path.join(dstdir, name))


def force_complex_orbital(filename):
    with open(filename) as f:
        lines = f.readlines()

    replaced = False
    for i, line in enumerate(lines):
        if line.split()[:1] == ["ComplexType"]:
            lines[i] = "ComplexType          1\n"
            replaced = True
            break

    if not replaced:
        raise RuntimeError("ComplexType line was not found in {}".format(filename))

    with open(filename, "w") as f:
        f.writelines(lines)


def read_two_body_green(filename):
    rows = []
    with open(filename) as f:
        for line in f:
            cols = line.split()
            if len(cols) != 10:
                continue
            try:
                idx = tuple(int(x) for x in cols[:8])
                value = complex(float(cols[8]), float(cols[9]))
            except ValueError:
                continue
            rows.append((idx, value))
    return rows


def main():
    if len(sys.argv) == 1:
        print("usage: {} <model name> [--complex-orbital]".format(sys.argv[0]))
        return -1

    rootdir = os.getcwd()
    model = sys.argv[1]
    force_complex = "--complex-orbital" in sys.argv[2:]
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    refdir = os.path.join(rootdir, "data", model)
    workdir = os.path.join(rootdir, "work", model)
    if force_complex:
        workdir += "_complex"
    if mpi_procs:
        workdir += "_mpi{}".format(mpi_procs)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    os.chdir(workdir)

    copy_def_files(refdir, workdir)
    if force_complex:
        force_complex_orbital("orbitalidx.def")

    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    cmd = [bin_to_test, "-e", "namelist.def"]
    if os.path.exists("initial.def"):
        cmd.append("initial.def")
    if mpi_procs:
        cmd = ["mpirun", "-np", mpi_procs] + cmd
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )
    with open("doublon_only_test.log", "w") as f:
        f.write(proc.stdout)
    if proc.returncode != 0:
        print(proc.stdout)
        return proc.returncode

    output = os.path.join("output", "zvo_cisajscktalt_001.dat")
    if not os.path.exists(output):
        print("ERROR: {} was not generated".format(output))
        return -1

    rows = read_two_body_green(output)
    if not rows:
        print("ERROR: no two-body Green-function data found in {}".format(output))
        return -1

    result = 0
    for idx, value in rows:
        if math.hypot(value.real, value.imag) > TOL:
            print(
                "ERROR: sector-breaking GreenFunc2 is non-zero for {}: {:.18e} {:.18e}".format(
                    idx, value.real, value.imag
                )
            )
            result = -1

    return result


if __name__ == "__main__":
    sys.exit(main())
