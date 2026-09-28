from __future__ import print_function

import os
import shutil
import subprocess
import sys

# Kondo chain with 4 local spins and 4 conduction sites, above half filling
# of the conduction band. The itinerant electrons of one spin outnumbered
# the conduction sites with some of the seeds, depending on the directions
# of the local spins, and vmc.out looked for a free site forever.
STDFACE = """L = 4
model = "{model}"
lattice = "chain"
t = 1.0
J = 1.0
ncond = {ncond}
NSROptItrStep = 2
NSROptItrSmp = 1
NVMCSample = 10
RndSeed = {seed}
{extra}
"""

VARIANTS = {
    "real": ("Kondolattice", "2Sz = 0"),
    "cmp": ("Kondolattice", "2Sz = 0\nComplexType = 1"),
    "fsz": ("KondolatticeGC", "ComplexType = 1"),
}

SEEDS = range(21, 41)
TIMEOUT = 20  # seconds for one run, which takes less than 1 second


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in VARIANTS:
        print("usage: {} <{}>".format(sys.argv[0], "|".join(sorted(VARIANTS))))
        return -1
    variant = sys.argv[1]
    model, extra = VARIANTS[variant]

    rootdir = os.getcwd()
    workdir = os.path.join(rootdir, "work", "InitialSample_" + variant)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    os.chdir(workdir)
    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")

    status = 0
    for ncond in (6, 8):
        for seed in SEEDS:
            rundir = os.path.join(workdir,
                                  "ncond{}_seed{}".format(ncond, seed))
            os.makedirs(rundir)
            with open(os.path.join(rundir, "StdFace.def"), "w") as f:
                f.write(STDFACE.format(model=model, ncond=ncond, seed=seed,
                                       extra=extra))
            try:
                proc = subprocess.run(
                    [bin_to_test, "-s", "StdFace.def"], cwd=rundir,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, timeout=TIMEOUT)
            except subprocess.TimeoutExpired:
                print("ERROR: ncond = {}, RndSeed = {}: vmc.out did not end "
                      "in {} s".format(ncond, seed, TIMEOUT))
                status = -1
                continue
            if proc.returncode != 0:
                print("ERROR: ncond = {}, RndSeed = {}: exit code {}".format(
                    ncond, seed, proc.returncode))
                print(proc.stdout)
                status = -1
                continue
            output = os.path.join(rundir, "output", "zvo_out_001.dat")
            with open(output) as f:
                steps = len(f.readlines())
            if steps != 2:
                print("ERROR: ncond = {}, RndSeed = {}: {} steps in "
                      "zvo_out_001.dat".format(ncond, seed, steps))
                status = -1
    return status


if __name__ == "__main__":
    sys.exit(main())
