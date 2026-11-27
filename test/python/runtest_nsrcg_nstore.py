from __future__ import print_function

"""NStore must not change SR-CG results.

For NSRCG!=0 the stored-O path is always taken, so NStore=0 and NStore=1
must give identical zqp_opt.dat.  NSRCG=2 (SR-CG with diagonal scaling) with
NStore=0 used to be accepted while the storage was allocated only for
NSRCG==1 (Julia-mVMC #56); this test pins the accepted combinations:

  NSRCG=1 NStore=0  ==  NSRCG=1 NStore=1
  NSRCG=2 NStore=0  ==  NSRCG=2 NStore=1   (and the DiagScale remark is printed)

Set MVMC_MPI_PROCS to run every case through mpirun.
"""

import math
import os
import shutil
import subprocess
import sys


COMMON_MODPARA_UPDATES = {
    "NSROptItrStep": 6,
    "NSROptItrSmp": 2,
    "NVMCWarmUp": 2,
    "NVMCSample": 100,
    "NSPGaussLeg": 2,
    "NSROptCGMaxIter": 1000,
}

DIAGSCALE_REMARK = "remark: use preconditioned CG (Diag Scale)"


def copy_defs(src, dst):
    for fn in os.listdir(src):
        if fn.endswith(".def"):
            shutil.copy(os.path.join(src, fn), os.path.join(dst, fn))


def update_modpara(filename, updates):
    found = set()
    lines = []
    with open(filename) as f:
        for line in f:
            words = line.split()
            if words and words[0] in updates:
                lines.append("{:<15} {}\n".format(words[0], updates[words[0]]))
                found.add(words[0])
            else:
                lines.append(line)
    for key, value in updates.items():
        if key not in found:
            lines.append("{:<15} {}\n".format(key, value))
    with open(filename, "w") as f:
        f.writelines(lines)


def read_opt(path):
    with open(path) as f:
        values = [float(x) for x in f.read().split()]
    if not values:
        raise RuntimeError("{} is empty".format(path))
    return values


def run_case(rootdir, refdir, workroot, name, updates):
    workdir = os.path.join(workroot, name)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    copy_defs(refdir, workdir)
    update_modpara(os.path.join(workdir, "modpara.def"), updates)

    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    cmd = [bin_to_test, "-e", "namelist.def"]
    if os.path.exists(os.path.join(workdir, "initial.def")):
        cmd.append("initial.def")
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    if mpi_procs:
        cmd = ["mpirun", "-np", mpi_procs] + cmd
    proc = subprocess.run(cmd, cwd=workdir, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True)
    with open(os.path.join(workdir, "run.log"), "w") as f:
        f.write(proc.stdout)
    if proc.returncode != 0:
        print(proc.stdout)
        raise RuntimeError("{} failed with exit code {}".format(name, proc.returncode))
    return read_opt(os.path.join(workdir, "output", "zqp_opt.dat")), proc.stdout


def compare(label, left, right, tol):
    if len(left) != len(right):
        print("{} length mismatch: {} vs {}".format(label, len(left), len(right)))
        return 1
    max_abs_diff = 0.0
    for a, b in zip(left, right):
        if math.isnan(a) or math.isinf(a) or math.isnan(b) or math.isinf(b):
            print("{} non-finite value: {} vs {}".format(label, a, b))
            return 1
        max_abs_diff = max(max_abs_diff, abs(a - b))
    if max_abs_diff > tol:
        print("{} mismatch: max_abs_diff={:.3e}".format(label, max_abs_diff))
        return 1
    print("{} OK: max_abs_diff={:.3e}".format(label, max_abs_diff))
    return 0


def main():
    if len(sys.argv) != 2:
        print("usage: {} <model name>".format(sys.argv[0]))
        return -1

    model = sys.argv[1]
    rootdir = os.getcwd()
    refdir = os.path.join(rootdir, "data", model)
    suffix = "_mpi" if os.environ.get("MVMC_MPI_PROCS") else ""
    workroot = os.path.join(rootdir, "work",
                            "{}_nsrcg_nstore{}".format(model, suffix))
    if os.path.exists(workroot):
        shutil.rmtree(workroot)
    os.makedirs(workroot)

    cases = {}
    for nsrcg in (1, 2):
        for nstore in (0, 1):
            name = "cg{}_store{}".format(nsrcg, nstore)
            cases[name] = dict(COMMON_MODPARA_UPDATES, NSRCG=nsrcg, NStore=nstore)

    results = {}
    logs = {}
    for name in sorted(cases):
        results[name], logs[name] = run_case(rootdir, refdir, workroot, name, cases[name])

    result = 0
    for name in ("cg2_store0", "cg2_store1"):
        if DIAGSCALE_REMARK not in logs[name]:
            print("{}: missing '{}'".format(name, DIAGSCALE_REMARK))
            result = 1
    # Same arithmetic path, so the results are expected to be identical.
    tol = 1.0e-12
    result |= compare("NSRCG=1: NStore=0 vs NStore=1", results["cg1_store0"], results["cg1_store1"], tol)
    result |= compare("NSRCG=2: NStore=0 vs NStore=1", results["cg2_store0"], results["cg2_store1"], tol)
    return -1 if result else 0


if __name__ == "__main__":
    sys.exit(main())
