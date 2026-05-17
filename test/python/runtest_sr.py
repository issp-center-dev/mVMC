from __future__ import print_function

import os
import shutil
import subprocess
import sys

import numpy as np


def read_out(filename):
    # drop the first two columns
    array = np.loadtxt(filename, dtype="float").astype("float")
    return array


if len(sys.argv) == 1:
    print("usage: {} <model name>".format(sys.argv[0]))
    sys.exit(-1)

rootdir = os.getcwd()
refdir = os.path.join(rootdir, "data", sys.argv[1])
workdir = os.path.join(rootdir, "work", sys.argv[1])
if os.path.exists(workdir):
    shutil.rmtree(workdir)
os.makedirs(workdir)
os.chdir(workdir)

bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")

result = subprocess.call([bin_to_test, "-s", "%s/StdFace.def" % refdir])
if result != 0:
    sys.exit(result)

array_calc_sr = read_out("./output/zqp_opt.dat")[4:]
subprocess.call(["mv", "output", "output_sr"])


result_cg = subprocess.call([bin_to_test, "-s", "%s/StdFace_CG.def" % refdir])
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
