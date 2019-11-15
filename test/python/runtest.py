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

array_calc = read_out("./output/zqp_opt.dat")[0:2]
ref_ave = read_out("%s/ref/ref_mean.dat" % refdir)[0:2]
ref_std = read_out("%s/ref/ref_std.dat" % refdir)[0:2]

result = 0
for diff, s in zip(array_calc - ref_ave, ref_std):
    diff = abs(diff)
    if diff >= 3 * s and diff >= 1e-8:
        result = -1

sys.exit(result)
