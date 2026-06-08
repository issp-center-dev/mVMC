from __future__ import print_function

import os
import shutil
import subprocess
import sys

import numpy as np


def read_out(filename):
    return np.loadtxt(filename, dtype="float").astype("float").reshape(-1)


def copy_def_files(srcdir, dstdir):
    for name in os.listdir(srcdir):
        if name.endswith(".def"):
            shutil.copy(os.path.join(srcdir, name), os.path.join(dstdir, name))


def enable_reweight(filename):
    with open(filename) as f:
        lines = f.readlines()

    replaced = False
    for i, line in enumerate(lines):
        if line.split()[:1] == ["reweight"]:
            lines[i] = "reweight       1\n"
            replaced = True
            break

    if not replaced:
        lines.append("reweight       1\n")

    with open(filename, "w") as f:
        f.writelines(lines)


def check_finite_output(filename):
    if not os.path.exists(filename):
        print("ERROR: {} was not generated".format(filename))
        return -1

    data = np.loadtxt(filename, dtype="float").astype("float")
    if data.size == 0:
        print("ERROR: {} is empty".format(filename))
        return -1
    if not np.all(np.isfinite(data)):
        print("ERROR: {} contains non-finite values".format(filename))
        return -1
    return 0


def check_reference_output(refdir):
    array_calc = read_out(os.path.join("output", "zqp_opt.dat"))[0:2]
    ref_ave = read_out(os.path.join(refdir, "ref", "ref_mean.dat"))[0:2]
    ref_std = read_out(os.path.join(refdir, "ref", "ref_std.dat"))[0:2]

    result = 0
    for idx, (diff, s) in enumerate(zip(array_calc - ref_ave, ref_std)):
        diff = abs(diff)
        if diff >= 3 * s and diff >= 1e-8:
            print(
                "ERROR: zqp_opt.dat[{}] mismatch: diff={} std={}".format(
                    idx, diff, s
                )
            )
            result = -1
    return result


def main():
    if len(sys.argv) == 1:
        print("usage: {} <model name>".format(sys.argv[0]))
        return -1

    rootdir = os.getcwd()
    model = sys.argv[1]
    refdir = os.path.join(rootdir, "data", model)
    workdir = os.path.join(rootdir, "work", model + "_reweight")
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    os.chdir(workdir)

    copy_def_files(refdir, workdir)
    enable_reweight("modpara.def")

    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    cmd = [bin_to_test, "-e", "namelist.def"]
    if os.path.exists("initial.def"):
        cmd.append("initial.def")
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )
    with open("reweight_test.log", "w") as f:
        f.write(proc.stdout)
    if proc.returncode != 0:
        print(proc.stdout)
        return proc.returncode

    result = check_reference_output(refdir)
    if result != 0:
        return result
    return check_finite_output(os.path.join("output", "zvo_out_001.dat"))


if __name__ == "__main__":
    sys.exit(main())
