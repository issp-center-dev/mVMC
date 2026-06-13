from __future__ import print_function

import os
import shutil
import subprocess
import sys


def main():
    if len(sys.argv) < 3:
        print("usage: {} <model name> <expected_error_substring> [allow_success]".format(sys.argv[0]))
        return -1

    model = sys.argv[1]
    expected = sys.argv[2]
    allow_success = (len(sys.argv) >= 4 and sys.argv[3] == "allow_success")

    rootdir = os.getcwd()
    refdir = os.path.join(rootdir, "data", model)
    workdir = os.path.join(rootdir, "work", model)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    os.chdir(workdir)

    # Copy all definition files required by expert mode.
    for fn in os.listdir(refdir):
        if fn.endswith(".def"):
            shutil.copy(os.path.join(refdir, fn), os.path.join(workdir, fn))

    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    proc = subprocess.run(
        [bin_to_test, "-e", "namelist.def"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    with open("invalid_test.log", "w") as f:
        f.write(proc.stdout)

    if expected not in proc.stdout:
        print("ERROR: expected error substring not found: {}".format(expected))
        print("---- output begin ----")
        print(proc.stdout)
        print("---- output end ----")
        return -1
    # Default: invalid input must terminate with non-zero.
    # Some legacy readers only emit error messages and continue;
    # those cases can opt in to allow_success mode.
    if not allow_success and proc.returncode == 0:
        print("ERROR: vmc.out unexpectedly succeeded for invalid input.")
        return -1
    return 0


if __name__ == "__main__":
    sys.exit(main())
