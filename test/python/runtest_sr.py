from __future__ import print_function

import os
import shutil
import subprocess
import sys

import numpy as np


CG_MAX_ITER_OVERRIDES = {
    # The default max_iter=nSmat is too tight for this direct-vs-CG comparison.
    "HubbardChain_cmp_DiagCG": 128,
}


def read_out(filename):
    # drop the first two columns
    array = np.loadtxt(filename, dtype="float")
    if array.ndim == 1:
        array = array[2:]
    else:
        array = array[:, 2:]
    return array.astype("float")


def update_modpara(filename, values):
    found = set()
    lines = []
    with open(filename) as f:
        for line in f:
            words = line.split()
            if words and words[0] in values:
                lines.append("{:<15} {}\n".format(words[0], values[words[0]]))
                found.add(words[0])
            else:
                lines.append(line)

    for key, value in values.items():
        if key not in found:
            lines.append("{:<15} {}\n".format(key, value))

    with open(filename, "w") as f:
        f.writelines(lines)


def run_vmc(bin_to_test, dry_bin_to_test, stdef, modpara_updates=None):
    if not modpara_updates:
        return subprocess.call([bin_to_test, "-s", stdef])

    result = subprocess.call([dry_bin_to_test, stdef])
    if result != 0:
        return result
    update_modpara("modpara.def", modpara_updates)
    return subprocess.call([bin_to_test, "-e", "namelist.def"])


if len(sys.argv) == 1:
    print("usage: {} <model name> [--reweight]".format(sys.argv[0]))
    sys.exit(-1)

# with --reweight, run both the direct-SR and the SR-CG side with
# reweight=1 so the w-weighted <O> accumulation of the CG Store path
# is compared against the dposv path
with_reweight = "--reweight" in sys.argv[2:]
common_modpara_updates = {"reweight": 1} if with_reweight else {}

rootdir = os.getcwd()
refdir = os.path.join(rootdir, "data", sys.argv[1])
workname = sys.argv[1] + ("_reweight" if with_reweight else "")
workdir = os.path.join(rootdir, "work", workname)
if os.path.exists(workdir):
    shutil.rmtree(workdir)
os.makedirs(workdir)
os.chdir(workdir)

bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
dry_bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmcdry.out")

result = run_vmc(bin_to_test, dry_bin_to_test, "%s/StdFace.def" % refdir,
                 dict(common_modpara_updates))
if result != 0:
    sys.exit(result)

array_calc_sr = read_out("./output/zqp_opt.dat")[4:]
subprocess.call(["mv", "output", "output_sr"])


cg_modpara_updates = dict(common_modpara_updates)
if sys.argv[1] in CG_MAX_ITER_OVERRIDES:
    cg_modpara_updates["NSROptCGMaxIter"] = CG_MAX_ITER_OVERRIDES[sys.argv[1]]

result_cg = run_vmc(
    bin_to_test,
    dry_bin_to_test,
    "%s/StdFace_CG.def" % refdir,
    cg_modpara_updates,
)
if result_cg != 0:
    sys.exit(result_cg)


array_calc_cg = read_out("./output/zqp_opt.dat")[4:]

result = 0
if array_calc_sr.shape != array_calc_cg.shape:
    print("shape mismatch: sr={} cg={}".format(array_calc_sr.shape, array_calc_cg.shape))
    sys.exit(-1)

abs_diff = np.abs(array_calc_sr - array_calc_cg)
mean_abs_diff = np.mean(abs_diff)
max_abs_diff = np.max(abs_diff)
if mean_abs_diff >= 1e-8:
    print("SR-CG mismatch: mean_abs_diff={:.3e}, max_abs_diff={:.3e}".format(mean_abs_diff, max_abs_diff))
    result = -1

sys.exit(result)
