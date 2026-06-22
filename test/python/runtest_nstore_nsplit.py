from __future__ import print_function

import os
import shutil
import subprocess
import sys


COMMON_MODPARA_UPDATES = {
    "NSROptItrStep": 12,
    "NSROptItrSmp": 4,
    "NVMCWarmUp": 4,
    "NVMCSample": 24,
    "NSRCG": 0,
}


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
    rows = []
    with open(path) as f:
        for line in f:
            cols = line.split()
            if cols:
                rows.append([float(x) for x in cols[2:]])
    if len(rows) > 4:
        rows = rows[4:]
    return rows


def run_case(rootdir, refdir, workroot, model, name, nproc, updates):
    workdir = os.path.join(workroot, name)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    copy_defs(refdir, workdir)
    update_modpara(os.path.join(workdir, "modpara.def"), updates)

    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    cmd = ["mpirun", "-np", str(nproc), bin_to_test, "-e", "namelist.def"]
    if os.path.exists(os.path.join(workdir, "initial.def")):
        cmd.append("initial.def")

    result = subprocess.call(cmd, cwd=workdir)
    if result != 0:
        raise RuntimeError("{} failed with exit code {}".format(name, result))

    return read_opt(os.path.join(workdir, "output", "zqp_opt.dat"))


def compare(label, left, right, tol):
    left_shape = (len(left), len(left[0]) if left else 0)
    right_shape = (len(right), len(right[0]) if right else 0)
    if left_shape != right_shape:
        print("{} shape mismatch: left={} right={}".format(label, left_shape, right_shape))
        return 1

    diffs = []
    for left_row, right_row in zip(left, right):
        diffs.extend(abs(a - b) for a, b in zip(left_row, right_row))
    if not diffs:
        print("{} has no comparable data".format(label))
        return 1

    mean_abs_diff = sum(diffs) / len(diffs)
    max_abs_diff = max(diffs)
    if mean_abs_diff >= tol:
        print(
            "{} mismatch: mean_abs_diff={:.3e}, max_abs_diff={:.3e}".format(
                label, mean_abs_diff, max_abs_diff
            )
        )
        return 1
    return 0


def main():
    if len(sys.argv) != 2:
        print("usage: {} <model name>".format(sys.argv[0]))
        return -1

    model = sys.argv[1]
    rootdir = os.getcwd()
    refdir = os.path.join(rootdir, "data", model)
    workroot = os.path.join(rootdir, "work", "{}_nstore_nsplit".format(model))
    if os.path.exists(workroot):
        shutil.rmtree(workroot)
    os.makedirs(workroot)

    nproc_ref = int(os.environ.get("MVMC_MPI_PROCS_REF", "2"))
    nproc_split = int(os.environ.get("MVMC_MPI_PROCS_SPLIT", "4"))

    cases = {
        "direct_ref": dict(COMMON_MODPARA_UPDATES, NSplitSize=1, NStore=0),
        "store_ref": dict(COMMON_MODPARA_UPDATES, NSplitSize=1, NStore=1),
        "direct_split": dict(COMMON_MODPARA_UPDATES, NSplitSize=2, NStore=0),
        "store_split": dict(COMMON_MODPARA_UPDATES, NSplitSize=2, NStore=1),
    }

    results = {
        "direct_ref": run_case(rootdir, refdir, workroot, model, "direct_ref", nproc_ref, cases["direct_ref"]),
        "store_ref": run_case(rootdir, refdir, workroot, model, "store_ref", nproc_ref, cases["store_ref"]),
        "direct_split": run_case(
            rootdir, refdir, workroot, model, "direct_split", nproc_split, cases["direct_split"]
        ),
        "store_split": run_case(
            rootdir, refdir, workroot, model, "store_split", nproc_split, cases["store_split"]
        ),
    }

    tol = 1.0e-8
    result = 0
    result |= compare("direct_ref vs direct_split", results["direct_ref"], results["direct_split"], tol)
    result |= compare("store_ref vs store_split", results["store_ref"], results["store_split"], tol)
    result |= compare("direct_split vs store_split", results["direct_split"], results["store_split"], tol)
    result |= compare("direct_ref vs store_ref", results["direct_ref"], results["store_ref"], tol)
    return -1 if result else 0


if __name__ == "__main__":
    sys.exit(main())
