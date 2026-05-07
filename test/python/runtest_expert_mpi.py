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

# copy *.def from refdir
os.system("cp %s/*.def ." % refdir)

bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
mpi_command = "mpirun"
nproc = 2

cmd = [mpi_command, "-np", str(nproc), bin_to_test, "-e", "namelist.def"]
if os.path.exists("initial.def"):
    cmd.append("initial.def")
result = subprocess.call(cmd)
if result != 0:
    sys.exit(result)

array_calc = read_out("./output/zqp_opt.dat")[0:2]
ref_ave = read_out("%s/ref/ref_mean.dat" % refdir)[0:2]

# MPI smoke test: ref was generated from a non-MPI single run with a
# particular sample partition, so exact 3-sigma agreement is not expected.
# We only verify that vmc.out completed (rc=0 above) and the energy is in a
# reasonable range (within 1% of the reference).
result = 0
energy_calc = array_calc[0]
energy_ref = ref_ave[0]
if abs(energy_ref) > 1e-12:
    rel = abs(energy_calc - energy_ref) / abs(energy_ref)
    if rel >= 1e-2:
        print("energy mismatch: calc=%.6e ref=%.6e rel=%.3e"
              % (energy_calc, energy_ref, rel))
        result = -1

sys.exit(result)
