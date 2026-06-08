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
mpi_procs = os.environ.get("MVMC_MPI_PROCS")

#result = subprocess.call([bin_to_test, "-s", "%s/StdFace.def" % refdir])
#result = subprocess.call([bin_to_test, "-e", "%s/namelist.def" % refdir, "%s/initial.def" % refdir])
if mpi_procs:
    result = subprocess.call(["mpirun", "-np", mpi_procs, bin_to_test, "-e", "namelist.def", "initial.def"])
else:
    result = subprocess.call([bin_to_test, "-e", "namelist.def", "initial.def"])
if result != 0:
    sys.exit(result)

if "Twist" in sys.argv[1]:
  output_file = "./output/zvo_twist_001.dat"

array_calc = read_out(output_file)

num_lines = array_calc.shape[0]
ref_ave = read_out("%s/ref/ref_mean.dat" % refdir)
ref_std = read_out("%s/ref/ref_std.dat" % refdir)

result = 0
if "Twist" in sys.argv[1]:
    # Columns 0-2 (Re, Im, |z|): per-component 3-sigma check.
    diff_phys = np.abs(array_calc[:, 0:3] - ref_ave[:, 0:3])
    if np.any(diff_phys >= 3 * ref_std[:, 0:3]):
        result = -1
    # Column 3 (arg) vs reference is intentionally not checked: the
    # complex value is already pinned by columns 0 and 1, and a direct
    # phase difference carries a 2*pi ambiguity near the +/-pi branch.
    # Self-consistency: polar form (|z| * exp(i * arg)) must equal the
    # rectangular form (Re + i * Im). Differenced in the complex plane
    # so the test is wrap-free regardless of the carg branch.
    self_polar = array_calc[:, 2] * np.exp(1j * array_calc[:, 3])
    self_rect  = array_calc[:, 0] + 1j * array_calc[:, 1]
    if np.any(np.abs(self_polar - self_rect) >= 1.0e-10):
        result = -1
else:
    for diff, s in zip(array_calc - ref_ave, ref_std):
        diff = np.abs(diff)
        if np.any(diff >= 3 * s):
            result = -1

sys.exit(result)
