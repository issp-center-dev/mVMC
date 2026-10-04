from __future__ import print_function

import os
import shutil
import subprocess
import sys

import numpy as np


FIXTURE = "HubbardChainLanczos_Guide"
NBIN = 8
MODES = {
    "legacy": {"DLanczosGuideEps": "0.0"},
    "eps1e-6": {"DLanczosGuideEps": "1.0e-6"},
    "eps1": {"DLanczosGuideEps": "1.0"},
    "invalid_negative": {"DLanczosGuideEps": "-1.0e-6"},
    "invalid_estimator": {
        "DLanczosGuideEps": "1.0e-6",
        "NLanczosEstimatorMode": "1",
    },
}


def set_modpara(path, overrides):
    with open(path) as handle:
        lines = handle.read().splitlines()
    seen = set()
    out = []
    for line in lines:
        fields = line.split()
        key = fields[0] if fields else ""
        if key in overrides:
            out.append("%-24s %s" % (key, overrides[key]))
            seen.add(key)
        else:
            out.append(line)
    for key, value in overrides.items():
        if key not in seen:
            out.append("%-24s %s" % (key, value))
    with open(path, "w") as handle:
        handle.write("\n".join(out) + "\n")


def run_case(rootdir, namespace, name, overrides):
    refdir = os.path.join(rootdir, "data", FIXTURE)
    workdir = os.path.join(
        rootdir, "work", FIXTURE + "_" + namespace + "_" + name)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    shutil.copytree(refdir, workdir)
    set_modpara(os.path.join(workdir, "modpara.def"), overrides)
    binary = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    cmd = [binary, "-e", "namelist.def", "zqp_opt.dat"]
    procs = os.environ.get("MVMC_MPI_PROCS")
    if procs:
        cmd = ["mpirun", "-np", procs] + cmd
    return subprocess.call(cmd, cwd=workdir), workdir


def lanczos_energies(workdir):
    values = []
    for idx in range(1, NBIN + 1):
        row = np.loadtxt(
            os.path.join(workdir, "output", "zvo_ls_out_%03d.dat" % idx),
            dtype="float",
        ).flatten()
        values.append(row[0])
    return np.array(values)


def check_stats(workdir, mode):
    eps = float(MODES[mode]["DLanczosGuideEps"])
    for idx in range(1, NBIN + 1):
        path = os.path.join(
            workdir, "output", "zvo_nzguide_%03d.dat" % idx)
        if not os.path.exists(path):
            print("missing %s" % path)
            return 1
        row = np.loadtxt(path, dtype="float").flatten()
        if row.shape != (9,):
            print("unexpected column count in %s: %s" % (path, row))
            return 1
        (eps_file, n, floored, exact_zero, numeric_skip, sum_w, sum_w2,
         min_w, max_ratio) = row
        ok = (
            np.all(np.isfinite(row))
            and abs(eps_file - eps) <= 1e-12 * eps
            and n > 0
            and 0 <= floored <= n
            and exact_zero == 0
            and numeric_skip == 0
            and sum_w > 0
            and sum_w2 > 0
            and min_w > 0
            and max_ratio >= 1 - 1e-12
        )
        if mode == "eps1":
            ok = ok and floored > 0
        if not ok:
            print("bad guide statistics in %s: %s" % (path, row))
            return 1
    return 0


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in MODES:
        print("usage: %s <%s>" % (sys.argv[0], "|".join(sorted(MODES))))
        return 2
    mode = sys.argv[1]
    rootdir = os.getcwd()
    procs = os.environ.get("MVMC_MPI_PROCS")
    namespace = mode + ("_mpi" + procs if procs else "_serial")
    if mode.startswith("invalid"):
        status, _ = run_case(rootdir, namespace, mode, MODES[mode])
        print("%s: vmc.out exit status %d (nonzero expected)" %
              (mode, status))
        return 0 if status != 0 else 1

    status_legacy, legacy_dir = run_case(
        rootdir, namespace, "legacy", MODES["legacy"])
    if status_legacy != 0:
        return 1
    if mode == "legacy":
        return 0

    status_guide, guide_dir = run_case(
        rootdir, namespace, mode, MODES[mode])
    if status_guide != 0:
        return 1
    e_legacy = lanczos_energies(legacy_dir)
    e_guide = lanczos_energies(guide_dir)
    if not (np.all(np.isfinite(e_legacy)) and np.all(np.isfinite(e_guide))):
        print("non-finite Lanczos energy")
        return 1
    se = np.sqrt(e_legacy.var(ddof=1) / NBIN +
                 e_guide.var(ddof=1) / NBIN)
    diff = abs(e_legacy.mean() - e_guide.mean())
    print("legacy %.6f +- %.6f  guide(%s) %.6f +- %.6f  "
          "diff %.6f  3se %.6f" %
          (e_legacy.mean(), np.sqrt(e_legacy.var(ddof=1) / NBIN), mode,
           e_guide.mean(), np.sqrt(e_guide.var(ddof=1) / NBIN),
           diff, 3 * se))
    if diff >= 3 * se and diff >= 1e-8:
        return 1
    if os.path.exists(
            os.path.join(legacy_dir, "output", "zvo_nzguide_001.dat")):
        print("legacy run must not write guide statistics")
        return 1
    return check_stats(guide_dir, mode)


if __name__ == "__main__":
    sys.exit(main())
