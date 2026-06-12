from __future__ import print_function

import os
import shutil
import subprocess
import sys

import numpy as np


# Expert-mode smoke models: RescaleSmat=1 together with an RBM block must
# complete normally. Before the sltStart/sltEnd fix in fn_Rescale4SRCG the
# RBM O-columns were rescaled without the compensating parameter division,
# which made the optimization blow up and abort within a few SR steps.
SMOKE_MODELS = {
    "GeneralRBM_cmp": {
        "NSRCG": 1,
        "RescaleSmat": 1,
        "NSROptCGMaxIter": 5000,
        "NSROptItrStep": 200,
        "NSROptItrSmp": 50,
    },
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


def copy_def_files(srcdir, dstdir):
    for name in os.listdir(srcdir):
        if name.endswith(".def"):
            shutil.copy(os.path.join(srcdir, name), os.path.join(dstdir, name))


if len(sys.argv) == 1:
    print("usage: {} <model name>".format(sys.argv[0]))
    sys.exit(-1)

rootdir = os.getcwd()
refdir = os.path.join(rootdir, "data", sys.argv[1])
workdir = os.path.join(rootdir, "work", "Rescale_" + sys.argv[1])
if os.path.exists(workdir):
    shutil.rmtree(workdir)
os.makedirs(workdir)
os.chdir(workdir)

bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
dry_bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmcdry.out")

if sys.argv[1] in SMOKE_MODELS:
    # expert-mode smoke test: the run must finish without SR-CG abort
    copy_def_files(refdir, ".")
    update_modpara("modpara.def", SMOKE_MODELS[sys.argv[1]])
    result = subprocess.call([bin_to_test, "-e", "namelist.def"])
    if result != 0:
        print("RescaleSmat smoke run failed with exit code {}".format(result))
    sys.exit(result)

# standard-mode equivalence test: RescaleSmat is a pure reparametrization of
# the Slater block, so the optimized parameters must agree with the
# RescaleSmat=0 run up to the SR-CG solver tolerance.
result = subprocess.call([bin_to_test, "-s", "%s/StdFace_CG.def" % refdir])
if result != 0:
    sys.exit(result)

array_calc_rs0 = read_out("./output/zqp_opt.dat")[4:]
subprocess.call(["mv", "output", "output_rs0"])

result = subprocess.call([dry_bin_to_test, "%s/StdFace_CG.def" % refdir])
if result != 0:
    sys.exit(result)
update_modpara("modpara.def", {"RescaleSmat": 1})
result = subprocess.call([bin_to_test, "-e", "namelist.def"])
if result != 0:
    sys.exit(result)

array_calc_rs1 = read_out("./output/zqp_opt.dat")[4:]

result = 0
if array_calc_rs0.shape != array_calc_rs1.shape:
    print("shape mismatch: rs0={} rs1={}".format(array_calc_rs0.shape, array_calc_rs1.shape))
    sys.exit(-1)

abs_diff = np.abs(array_calc_rs0 - array_calc_rs1)
mean_abs_diff = np.mean(abs_diff)
max_abs_diff = np.max(abs_diff)
if mean_abs_diff >= 1e-8:
    print("RescaleSmat mismatch: mean_abs_diff={:.3e}, max_abs_diff={:.3e}".format(mean_abs_diff, max_abs_diff))
    result = -1

sys.exit(result)
