"""Tests for the OutputDir keyword in modpara.def.

usage: runtest_output_dir.py <case>

cases:
  nested        OutputDir is a relative, not-yet-existing nested directory
  absolute      OutputDir is an absolute path
  mixedcase     the keyword is matched case-insensitively
  data_limit    "<OutputDir>/<CDataFileHead>" is exactly at the length limit
  para_limit    "<OutputDir>/<CParaFileHead>" is exactly at the length limit
  mpi_nested    nested OutputDir with 2 MPI processes
  file          OutputDir names an existing regular file -> error
  spaces        OutputDir value contains whitespace -> error
  data_toolong  "<OutputDir>/<CDataFileHead>" is one byte too long -> error
  para_toolong  "<OutputDir>/<CParaFileHead>" is one byte too long -> error
  mpi_toolong   too long OutputDir with 2 MPI processes -> error, no hang
  greenr2k      greenr2k reads the correlation files from OutputDir

The MPI cases use MVMC_MPIEXEC and MVMC_MPIEXEC_NUMPROC_FLAG
(default: mpirun -np).
"""

from __future__ import print_function

import os
import shutil
import subprocess
import sys
import tempfile

FIXTURE = "HubbardChain_SpinJastrow"
# Distinct heads, so that the tests notice if a head is ignored.
DATA_HEAD = "vmcdata"
PARA_HEAD = "vmcparam"
LONG_HEAD = "head_longer_than_the_other"
# Maximum length in bytes of "<OutputDir>/<head>" accepted by vmc.out.
PREFIX_MAX = 191
TIMEOUT = 300

GREENR2K_STDFACE = """\
L = 4
Lsub = 2
model = "Hubbard"
lattice = "chain"
U = 4.0
t = 1.0
nelec = 4
2Sz = 0
NVMCCalMode = 1
NVMCSample = 10
"""
GREENR2K_KPATH = "2 4\nG 0 0 0\nX 0.5 0 0\n4 1 1\n"


def prepare_workdir(rootdir, case):
    refdir = os.path.join(rootdir, "data", FIXTURE)
    workdir = os.path.join(rootdir, "work", "OutputDir_" + case)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    for name in os.listdir(refdir):
        if name.endswith(".def"):
            shutil.copy(os.path.join(refdir, name), workdir)
    return workdir


def dir_of_length(length):
    """Relative path of exactly `length` bytes with short components."""
    parts = []
    while length > 100:
        parts.append("d" * 99)
        length -= 100
    parts.append("d" * length)
    return "/".join(parts)


def write_modpara(workdir, value, data_head=DATA_HEAD, para_head=PARA_HEAD,
                  keyword="OutputDir"):
    """Set the file heads, shorten the run and append the OutputDir line."""
    path = os.path.join(workdir, "modpara.def")
    with open(path) as stream:
        lines = stream.readlines()
    replaced = {"CDataFileHead": data_head, "CParaFileHead": para_head,
                "NSROptItrStep": "5", "NSROptItrSmp": "1", "NVMCSample": "10"}
    result = []
    for line in lines:
        words = line.split()
        if len(words) == 2 and words[0] in replaced:
            line = "{}  {}\n".format(words[0], replaced[words[0]])
        result.append(line)
    result.append("{}      {}\n".format(keyword, value))
    with open(path, "w") as stream:
        stream.writelines(result)


def run(rootdir, workdir, nproc=None):
    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    cmd = [bin_to_test, "-e", "namelist.def"]
    if nproc is not None:
        mpiexec = os.environ.get("MVMC_MPIEXEC") or "mpirun"
        numproc_flag = os.environ.get("MVMC_MPIEXEC_NUMPROC_FLAG") or "-np"
        cmd = [mpiexec, numproc_flag, str(nproc)] + cmd
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    proc = subprocess.Popen(cmd, cwd=workdir, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    try:
        out, _ = proc.communicate(timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        proc.kill()
        out, _ = proc.communicate()
        print(out.decode("utf-8", "replace"))
        print("FAIL: vmc.out did not finish within {} s".format(TIMEOUT))
        return None, ""
    out = out.decode("utf-8", "replace")
    print(out)
    return proc.returncode, out


def run_greenr2k(rootdir, workdir):
    """vmc.out and greenr2k with an OutputDir containing '/'."""
    build = os.path.join(rootdir, "..", "..")
    with open(os.path.join(workdir, "StdFace.def"), "w") as stream:
        stream.write(GREENR2K_STDFACE)
    code = subprocess.call(
        [os.path.join(build, "src", "mVMC", "vmcdry.out"), "StdFace.def"],
        cwd=workdir)
    if code != 0:
        print("FAIL: vmcdry.out exited with {}".format(code))
        return 1
    with open(os.path.join(workdir, "geometry.dat"), "a") as stream:
        stream.write(GREENR2K_KPATH)
    output_dir = os.path.join("run", "out")
    write_modpara(workdir, output_dir)
    code, _ = run(rootdir, workdir)
    if code != 0:
        print("FAIL: vmc.out exited with {}".format(code))
        return 1
    proc = subprocess.Popen(
        [os.path.join(build, "tool", "greenr2k"), "namelist.def",
         "geometry.dat"],
        cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out, _ = proc.communicate(timeout=TIMEOUT)
    print(out.decode("utf-8", "replace"))
    if proc.returncode != 0:
        print("FAIL: greenr2k exited with {}".format(proc.returncode))
        return 1
    corr = os.path.join(workdir, output_dir, DATA_HEAD + "_corr.dat")
    if not os.path.isfile(corr):
        print("FAIL: greenr2k did not write {}".format(corr))
        return 1
    if os.path.exists(os.path.join(workdir, "output")):
        print("FAIL: default output/ directory was created")
        return 1
    return 0


def expect_success(rootdir, workdir, output_dir, data_head=DATA_HEAD,
                   para_head=PARA_HEAD, nproc=None):
    code, out = run(rootdir, workdir, nproc)
    if code != 0:
        print("FAIL: vmc.out exited with {}".format(code))
        return 1
    if nproc is not None and out.count("Start: Read ModPara File") != 1:
        # Serial binaries started by the launcher each read modpara.def.
        print("FAIL: vmc.out did not run as one {}-process MPI job".format(
            nproc))
        return 1
    target = output_dir if os.path.isabs(output_dir) \
        else os.path.join(workdir, output_dir)
    for name in [para_head + "_opt.dat", data_head + "_out_001.dat"]:
        if not os.path.isfile(os.path.join(target, name)):
            print("FAIL: {} was not written to {}".format(name, target))
            return 1
    if os.path.exists(os.path.join(workdir, "output")):
        print("FAIL: default output/ directory was created")
        return 1
    return 0


def expect_failure(rootdir, workdir, message, nproc=None):
    code, out = run(rootdir, workdir, nproc)
    if code is None:
        return 1
    if code == 0:
        print("FAIL: vmc.out succeeded but an error was expected")
        return 1
    if "is incorrect" in out:
        print("FAIL: OutputDir was rejected as an unknown keyword")
        return 1
    if "OutputDir" not in out or message not in out:
        print("FAIL: expected an OutputDir error containing '{}'".format(
            message))
        return 1
    if os.path.exists(os.path.join(workdir, "output")):
        print("FAIL: default output/ directory was created")
        return 1
    return 0


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    case = sys.argv[1]
    rootdir = os.getcwd()
    workdir = prepare_workdir(rootdir, case)

    if case in ("nested", "mpi_nested"):
        output_dir = os.path.join("run", "out")
        write_modpara(workdir, output_dir)
        nproc = 2 if case == "mpi_nested" else None
        return expect_success(rootdir, workdir, output_dir, nproc=nproc)
    if case == "absolute":
        # The build tree may be deep, so use a short temporary directory.
        tmpdir = tempfile.mkdtemp(prefix="mvmc_od_")
        try:
            output_dir = os.path.join(tmpdir, "abs", "out")
            if len(os.path.join(output_dir, PARA_HEAD)) > PREFIX_MAX:
                print("SKIP: temporary directory path is too long")
                return 77
            write_modpara(workdir, output_dir)
            return expect_success(rootdir, workdir, output_dir)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)
    if case == "mixedcase":
        output_dir = "mixed"
        write_modpara(workdir, output_dir, keyword="outputDIR")
        return expect_success(rootdir, workdir, output_dir)
    if case == "data_limit":
        output_dir = dir_of_length(PREFIX_MAX - 1 - len(LONG_HEAD))
        write_modpara(workdir, output_dir, data_head=LONG_HEAD)
        return expect_success(rootdir, workdir, output_dir,
                              data_head=LONG_HEAD)
    if case == "para_limit":
        output_dir = dir_of_length(PREFIX_MAX - 1 - len(LONG_HEAD))
        write_modpara(workdir, output_dir, para_head=LONG_HEAD)
        return expect_success(rootdir, workdir, output_dir,
                              para_head=LONG_HEAD)
    if case == "file":
        with open(os.path.join(workdir, "blocker"), "w") as stream:
            stream.write("not a directory\n")
        write_modpara(workdir, "blocker")
        return expect_failure(rootdir, workdir, "not a directory")
    if case == "spaces":
        write_modpara(workdir, "run results")
        code = expect_failure(rootdir, workdir, "without whitespace")
        if os.path.exists(os.path.join(workdir, "run")):
            print("FAIL: truncated OutputDir 'run' was created")
            return 1
        return code
    if case == "data_toolong":
        output_dir = dir_of_length(PREFIX_MAX - len(LONG_HEAD))
        write_modpara(workdir, output_dir, data_head=LONG_HEAD)
        return expect_failure(rootdir, workdir, "too long")
    if case == "para_toolong":
        output_dir = dir_of_length(PREFIX_MAX - len(LONG_HEAD))
        write_modpara(workdir, output_dir, para_head=LONG_HEAD)
        return expect_failure(rootdir, workdir, "too long")
    if case == "mpi_toolong":
        output_dir = dir_of_length(PREFIX_MAX - len(PARA_HEAD))
        write_modpara(workdir, output_dir)
        return expect_failure(rootdir, workdir, "too long", nproc=2)

    if case == "greenr2k":
        return run_greenr2k(rootdir, workdir)

    print("unknown case: {}".format(case))
    return 2


if __name__ == "__main__":
    sys.exit(main())
