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
if result != 0:
    sys.exit(result)


array_calc_cg = read_out("./output/zqp_opt.dat")[4:]

result = 0
diff = np.mean(array_calc_sr - array_calc_cg)
if diff >= 1e-8:
    result = -1

sys.exit(result)
