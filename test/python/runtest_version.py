from __future__ import print_function

import os
import re
import subprocess
import sys
import tempfile


def read_version(source_dir):
    header = os.path.join(source_dir, "src", "mVMC", "include", "version.h")
    with open(header) as f:
        text = f.read()
    numbers = [
        re.search(r"^#define\s+MVMC_VERSION_" + part + r"\s+(\d+)", text,
                  re.MULTILINE).group(1)
        for part in ("MAJOR", "MINOR", "PATCH")
    ]
    prerelease = re.search(
        r'^#define\s+MVMC_VERSION_PRERELEASE\s+"([^"]*)"', text,
        re.MULTILINE).group(1)
    version = ".".join(numbers)
    if prerelease:
        version += "-" + prerelease
    return version


def read_git_hash(binary_dir):
    header = os.path.join(binary_dir, "src", "mVMC", "version_git.h")
    with open(header) as f:
        text = f.read()
    return re.search(r'^#define\s+MVMC_GIT_HASH\s+"([^"]*)"', text,
                     re.MULTILINE).group(1)


def main():
    if len(sys.argv) != 3:
        print("usage: {} <source dir> <binary dir>".format(sys.argv[0]))
        return -1
    source_dir, binary_dir = sys.argv[1:3]

    git_hash = read_git_hash(binary_dir)
    if git_hash and not re.match(r"^[0-9a-f]{8}(-dirty)?$", git_hash):
        print("ERROR: unexpected form of the hash: {}".format(git_hash))
        return -1
    expected = "mVMC version " + read_version(source_dir)
    if git_hash:
        expected += " ({})".format(git_hash)

    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    if mpi_procs:
        launcher = [
            os.environ.get("MVMC_MPIEXEC", "mpirun"),
            os.environ.get("MVMC_MPIEXEC_NUMPROC_FLAG", "-np"),
            mpi_procs,
        ]
        programs = ["vmc.out"]
    else:
        launcher = []
        programs = ["vmc.out", "vmcdry.out"]

    status = 0
    scratch = os.path.join(os.path.abspath(binary_dir), "tmp")
    os.makedirs(scratch, exist_ok=True)
    for program in programs:
        for option in ("-v", "--version"):
            command = launcher + [
                os.path.join(binary_dir, "src", "mVMC", program), option]
            # Open MPI 4 singleton startup and teardown can race on the
            # shared ompi.<host>.<uid> directory during parallel CTest runs.
            # Isolate every invocation, including consecutive short runs.
            with tempfile.TemporaryDirectory(prefix="mpi-", dir=scratch) as session:
                env = os.environ.copy()
                env["OMPI_MCA_orte_tmpdir_base"] = session
                proc = subprocess.run(command, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE, text=True, env=env)
            # The version is printed once, also in a run with MPI.
            if proc.returncode != 0 or proc.stdout != expected + "\n":
                print("ERROR: {}".format(" ".join(command)))
                print("  expected: {!r}, exit code 0".format(expected + "\n"))
                print("  obtained: {!r}, exit code {}".format(
                    proc.stdout, proc.returncode))
                print("  stderr: {!r}".format(proc.stderr))
                status = -1
            else:
                print("OK: {} {}: {}".format(program, option, expected))
    return status


if __name__ == "__main__":
    sys.exit(main())
