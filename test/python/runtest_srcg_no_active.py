from __future__ import print_function

"""SR-CG with no active variational parameter must be a clean no-op.

Modes:
  optflag            fixture SRCG_NoActiveParameter (every OptFlag is 0,
                     parameters come from initial.def), NSRCG=1
  optflag_diagscale  the same fixture with NSRCG=2 (diagonal scaling)
  diagcut            HubbardChain_SpinJastrow with NSRCG=1 and
                     DSROptRedCut=2.0, so the diagonal cut removes every
                     parameter (nSmat == 0 although OptFlag == 1)

Checks: the run exits 0 without a crash or sanitizer report, prints the
"no active variational parameter" warning, every zvo_SRinfo.dat row has
Msize == 0, absRmax == 0, imax == -1 and info == 0, and in the OptFlag
modes zqp_opt.dat still equals initial.def (no parameter update).
Set MVMC_MPI_PROCS to run through mpirun.
"""

import os
import re
import shutil
import subprocess
import sys


MODES = {
    "optflag": ("SRCG_NoActiveParameter", {"NSRCG": "1"}, True),
    "optflag_diagscale": ("SRCG_NoActiveParameter", {"NSRCG": "2"}, True),
    "diagcut": ("HubbardChain_SpinJastrow",
                {"NSRCG": "1", "NStore": "1", "DSROptRedCut": "2.0",
                 "NSROptItrStep": "4", "NSROptItrSmp": "2", "NVMCWarmUp": "2",
                 "NVMCSample": "20", "NSPGaussLeg": "2"}, False),
}

WARNING = "SR-CG has no active variational parameter"
CRASH = re.compile(
    r"AddressSanitizer|LeakSanitizer|UndefinedBehaviorSanitizer|"
    r"runtime error:|Segmentation fault|Bus error|core dumped|BLAS error")


def update_modpara(filename, updates):
    found = set()
    lines = []
    with open(filename) as f:
        for line in f:
            words = line.split()
            if words and words[0] in updates:
                lines.append("{:<15} {}\n".format(words[0], updates[words[0]]))
                found.add(words[0])
            else:
                lines.append(line)
    for key, value in updates.items():
        if key not in found:
            lines.append("{:<15} {}\n".format(key, value))
    with open(filename, "w") as f:
        f.writelines(lines)


def read_values(path):
    with open(path) as f:
        return [float(x) for x in f.read().split()]


def check_srinfo(path):
    rows = 0
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split()
            if not cols:
                continue
            rows += 1
            msize = int(cols[1])
            rmax = float(cols[6])
            imax = int(cols[7])
            info = int(cols[8])
            if msize != 0 or rmax != 0.0 or imax != -1 or info != 0:
                raise AssertionError(
                    "unexpected SRinfo row (Msize={}, absRmax={}, imax={}, "
                    "info={}): {}".format(msize, rmax, imax, info, line.rstrip()))
    if rows == 0:
        raise AssertionError("zvo_SRinfo.dat has no data rows")
    return rows


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in MODES:
        print("usage: {} <{}>".format(sys.argv[0], "|".join(sorted(MODES))))
        return -1
    mode = sys.argv[1]
    model, updates, compare_initial = MODES[mode]

    rootdir = os.getcwd()
    refdir = os.path.join(rootdir, "data", model)
    workroot = os.path.join(rootdir, "work")
    os.makedirs(workroot, exist_ok=True)
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    workdir = os.path.join(
        workroot, "srcg_no_active_{}{}_{}".format(
            mode, "_mpi" if mpi_procs else "", os.getpid()))
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    for fn in os.listdir(refdir):
        if fn.endswith(".def"):
            shutil.copy(os.path.join(refdir, fn), os.path.join(workdir, fn))
    os.chdir(workdir)
    update_modpara("modpara.def", updates)

    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    command = [bin_to_test, "-e", "namelist.def"]
    if os.path.exists("initial.def"):
        command.append("initial.def")
    if mpi_procs:
        command = ["mpirun", "-np", mpi_procs] + command
    proc = subprocess.run(command, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True)
    with open("run.log", "w") as f:
        f.write(proc.stdout)

    try:
        if proc.returncode != 0:
            raise AssertionError("exit status {}".format(proc.returncode))
        if CRASH.search(proc.stdout):
            raise AssertionError("crash or sanitizer report")
        if WARNING not in proc.stdout:
            raise AssertionError("missing warning: {}".format(WARNING))
        rows = check_srinfo(os.path.join("output", "zvo_SRinfo.dat"))
        if compare_initial:
            initial = read_values("initial.def")
            final = read_values(os.path.join("output", "zqp_opt.dat"))
            if len(initial) != len(final):
                raise AssertionError("zqp_opt.dat length {} != initial.def {}"
                                     .format(len(final), len(initial)))
            # Triplets (Re, Im, uncertainty).  The first two triplets are
            # <H> and <H^2> (sampling noise, skipped by ReadInitParameter);
            # the parameters start at the third triplet.
            for i in range(6, len(initial), 3):
                if initial[i] != final[i] or initial[i + 1] != final[i + 1]:
                    raise AssertionError(
                        "parameter {} changed: {} {} -> {} {}".format(
                            i // 3, initial[i], initial[i + 1],
                            final[i], final[i + 1]))
    except AssertionError as error:
        print("ERROR: {}".format(error))
        print("---- output begin ----")
        print(proc.stdout)
        print("---- output end ----")
        return -1
    print("OK: mode={} SRinfo rows={}".format(mode, rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
