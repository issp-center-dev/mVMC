from __future__ import print_function

import argparse
import filecmp
import os
import shutil
import subprocess
import sys

from runtest_gc import mpi_command, write_fixture


OUTPUTS = (
    "zvo_out_001.dat",
    "zvo_cisajs_001.dat",
    "zvo_cisajscktalt_001.dat",
    "zvo_NBodyG_001.dat",
)


def canonicalize_modpara(path):
    rewritten = []
    with open(path) as stream:
        for line in stream:
            columns = line.split()
            if columns and columns[0] in ("NGrandCanonical", "NGCInitNelec"):
                continue
            if columns and columns[0] == "Ncond":
                line = "Ncond          2\n"
            rewritten.append(line)
    with open(path, "w") as stream:
        stream.writelines(rewritten)


def run(binary, workdir, procs):
    command = [binary, "-e", "namelist.def", "initial.def"]
    if procs > 1:
        command = mpi_command(procs, command)
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment["OPENBLAS_NUM_THREADS"] = "1"
    environment["MKL_NUM_THREADS"] = "1"
    process = subprocess.run(command, cwd=workdir, env=environment,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True)
    with open(os.path.join(workdir, "run.log"), "w") as stream:
        stream.write(process.stdout)
    if process.returncode != 0:
        raise AssertionError("canonical compatibility run failed:\n{}".format(
            process.stdout[-8000:]))


def compare_variant(root, baseline, candidate, procs):
    seed_dir = os.path.join(root, "seed_np{}".format(procs))
    baseline_dir = os.path.join(root, "baseline_np{}".format(procs))
    candidate_dir = os.path.join(root, "candidate_np{}".format(procs))
    for path in (seed_dir, baseline_dir, candidate_dir):
        if os.path.exists(path):
            shutil.rmtree(path)
    os.makedirs(seed_dir)
    write_fixture(seed_dir, mode=1, samples=20000, seed=82351,
                  nsplit=1, iterations=1)
    canonicalize_modpara(os.path.join(seed_dir, "modpara.def"))
    shutil.copytree(seed_dir, baseline_dir)
    shutil.copytree(seed_dir, candidate_dir)
    run(baseline, baseline_dir, procs)
    run(candidate, candidate_dir, procs)
    for filename in OUTPUTS:
        baseline_output = os.path.join(baseline_dir, "output", filename)
        candidate_output = os.path.join(candidate_dir, "output", filename)
        if not filecmp.cmp(baseline_output, candidate_output, shallow=False):
            raise AssertionError("canonical output differs for np={} file={}".format(
                procs, filename))
    print("canonical byte comparison passed: np={} files={}".format(
        procs, len(OUTPUTS)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--work-root", required=True)
    arguments = parser.parse_args()
    baseline = os.path.abspath(arguments.baseline)
    candidate = os.path.abspath(arguments.candidate)
    root = os.path.abspath(arguments.work_root)
    if os.path.exists(root):
        shutil.rmtree(root)
    os.makedirs(root)
    for procs in (1, 2):
        compare_variant(root, baseline, candidate, procs)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print("ERROR: {}".format(error))
        sys.exit(1)
