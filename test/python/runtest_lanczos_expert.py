from __future__ import print_function

import os
import shutil
import subprocess
import sys

import numpy as np


def read_out(filename):
    array = np.loadtxt(filename, dtype="float").astype("float")
    return array.flatten()


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

os.system("cp %s/*.def ." % refdir)

bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")

cmd = [bin_to_test, "-e", "namelist.def"]
if os.path.exists("initial.def"):
    cmd.append("initial.def")
mpi_procs = os.environ.get("MVMC_MPI_PROCS")
if mpi_procs:
    cmd = ["mpirun", "-np", mpi_procs] + cmd

result = subprocess.call(cmd)
if result != 0:
    sys.exit(result)

array_calc = read_out("./output/zvo_ls_out_001.dat")[0:2]
ref_ave = read_out("%s/ref/ref_mean_Els.dat" % refdir)[0:2]
ref_std = read_out("%s/ref/ref_std_Els.dat" % refdir)[0:2]

if not np.all(np.isfinite(array_calc)):
    print("ERROR: non-finite Lanczos output: {}".format(array_calc))
    sys.exit(-1)
if not np.all(np.isfinite(ref_ave)) or not np.all(np.isfinite(ref_std)):
    print("ERROR: non-finite Lanczos reference")
    sys.exit(-1)

result = 0
for diff, s in zip(array_calc - ref_ave, ref_std):
    diff = abs(diff)
    if diff >= 3 * s and diff >= 1e-8:
        result = -1

# NLanczosMode=2 produces zvo_ls_cisajs_001.dat; check all values are finite.
ls_cisajs = "./output/zvo_ls_cisajs_001.dat"
if os.path.exists(ls_cisajs):
    gf = read_out(ls_cisajs)
    if not np.all(np.isfinite(gf)):
        print("ERROR: non-finite values in {}".format(ls_cisajs))
        result = -1

sys.exit(result)
