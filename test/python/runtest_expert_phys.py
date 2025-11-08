from __future__ import print_function

import os
import shutil
import subprocess
import sys

import numpy as np


def read_out(filename):
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

#copy *.def from refdir
os.system("cp %s/*.def ." % refdir)

bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")

#result = subprocess.call([bin_to_test, "-s", "%s/StdFace.def" % refdir])
#result = subprocess.call([bin_to_test, "-e", "%s/namelist.def" % refdir, "%s/initial.def" % refdir])
result = subprocess.call([bin_to_test, "-e", "namelist.def", "initial.def"])
if result != 0:
    sys.exit(result)

if "Twist" in sys.argv[1]:
  output_file = "./output/zvo_twist_001.dat"

array_calc = read_out(output_file)
if "Twist" in sys.argv[1]:
  array_calc[:,3] = np.exp(1j * array_calc[:,3])

num_lines = array_calc.shape[0]
ref_ave = read_out("%s/ref/ref_mean.dat" % refdir)
ref_std = read_out("%s/ref/ref_std.dat" % refdir)

result = 0
for diff, s in zip(array_calc - ref_ave, ref_std):
    diff = np.abs(diff)
    #if diff >= 3 * s and diff >= 1e-8:
    if np.any(diff >= 3 * s):
        result = -1

sys.exit(result)
