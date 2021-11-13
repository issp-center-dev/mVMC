from __future__ import print_function

import glob
import os
import shutil
import subprocess
import sys

import numpy as np


def read_out(filename):
    # drop the first two columns
    array = np.loadtxt(filename, dtype="float").astype("float")
    return array


def read_param(filename):
    res = {}
    with open(filename) as f:
        for line in f:
            words = line.split()
            if len(words) <= 1:
                continue
            res[words[0]] = words[1]
    return res


def read_opt_def(filename, total_site):
    rf = open(filename, "r")
    _arr = np.zeros(total_site, dtype=complex)
    for line in rf.readlines()[5:]:
        line1 = line.split()
        _arr[int(line1[0])] = float(line1[1]) + 1j * float(line1[2])
    rf.close()
    return _arr


def read_trans_def(filename, total_site):
    rf = open(filename, "r")
    _arr = np.zeros(((total_site * 2) ** 2) * 2, dtype=complex).reshape(
        (total_site, 2, total_site, 2, 2)
    )
    for line in rf.readlines()[5:]:
        line1 = line.split()
        _arr[int(line1[0])][int(line1[1])][int(line1[2])][int(line1[3])] = float(
            line1[4]
        ) + 1j * float(line1[5])
    rf.close()
    return _arr


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

bin_to_test = os.path.join(rootdir, "..", "..", "src", "ComplexUHF", "UHF")

for _file in glob.glob("%s/*.def" % refdir):
    shutil.copyfile(_file, os.path.basename(_file))

# run
result = subprocess.call([bin_to_test, "namelist.def"])
if result != 0:
    sys.exit(result)

# get results
array_calc = read_out("./zvo_eigen.dat")
ref_ave = read_out("%s/ref/zvo_eigen.dat" % refdir)
for diff in array_calc - ref_ave:
    diff = abs(diff)
    if (diff >= 1e-8).any():
        sys.exit(-1)

param = read_param("modpara.def")
nsites = int(param["Nsite"])

# get results
array_calc = read_trans_def("./zvo_UHF_cisajs.dat", nsites)
ref_ave = read_trans_def("%s/ref/zvo_UHF_cisajs.dat" % refdir, nsites)
for diff in array_calc - ref_ave:
    diff = np.abs(diff)
    if (diff >= 1e-8).any():
        sys.exit(-1)

sys.exit(0)
