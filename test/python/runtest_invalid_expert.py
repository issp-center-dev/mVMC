from __future__ import print_function

import os
import shutil
import subprocess
import sys


def safe_path_component(value):
    chars = []
    for ch in value:
        if ch.isalnum() or ch in ("-", "_", "."):
            chars.append(ch)
        else:
            chars.append("_")
    return "".join(chars) or "case"


def make_workdir_name(model, updates):
    pieces = [model]
    for key in sorted(updates):
        pieces.append("{}_{}".format(key, updates[key]))
    return safe_path_component("__".join(pieces))


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


def inject_spin_flip_transfer(filename):
    with open(filename) as source:
        lines = source.readlines()
    for index, line in enumerate(lines):
        words = line.split()
        if len(words) != 6:
            continue
        try:
            int(words[0])
            source_spin = int(words[1])
            int(words[2])
            int(words[3])
            float(words[4])
            float(words[5])
        except ValueError:
            continue
        words[3] = str(1 - source_spin)
        lines[index] = " ".join(words) + "\n"
        with open(filename, "w") as destination:
            destination.writelines(lines)
        return
    raise RuntimeError("no Transfer data row was available for mutation")


def main():
    if len(sys.argv) < 3:
        print(
            "usage: {} <model name> <expected_error_substring> "
            "[allow_success] [KEY=VALUE ...]".format(sys.argv[0])
        )
        return -1

    model = sys.argv[1]
    expected = sys.argv[2]
    allow_success = False
    spin_flip_transfer = False
    reject_before_memory = False
    expected_error_count = None
    modpara_updates = {}
    for arg in sys.argv[3:]:
        if arg == "allow_success":
            allow_success = True
        elif arg == "spin_flip_transfer":
            spin_flip_transfer = True
        elif arg == "reject_before_memory":
            reject_before_memory = True
        elif arg.startswith("expected_error_count="):
            expected_error_count = int(arg.split("=", 1)[1])
        elif "=" in arg:
            key, value = arg.split("=", 1)
            modpara_updates[key] = value
        else:
            print("ERROR: unknown extra argument: {}".format(arg))
            return -1

    # Existing contract/error tests exercise the legacy estimator formulas.
    # Keep that classification explicit in the generated ModPara input; new
    # corrected-selector tests pass NLanczosEstimatorMode=0 themselves.
    try:
        uses_legacy_lanczos = int(
            modpara_updates.get("NLanczosMode", "0")) > 0
    except ValueError:
        uses_legacy_lanczos = False
    if (
        uses_legacy_lanczos
        and "NLanczosEstimatorMode" not in modpara_updates
    ):
        modpara_updates["NLanczosEstimatorMode"] = "1"

    rootdir = os.getcwd()
    refdir = os.path.join(rootdir, "data", model)
    workroot = os.path.join(rootdir, "work")
    if not os.path.exists(workroot):
        os.makedirs(workroot)
    case_name = make_workdir_name(model, modpara_updates)
    workdir = os.path.join(workroot, "{}_{}".format(case_name, os.getpid()))
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    os.chdir(workdir)

    # Copy all definition files required by expert mode.
    for fn in os.listdir(refdir):
        if fn.endswith(".def"):
            shutil.copy(os.path.join(refdir, fn), os.path.join(workdir, fn))

    if modpara_updates:
        update_modpara("modpara.def", modpara_updates)
    if spin_flip_transfer:
        try:
            inject_spin_flip_transfer("trans.def")
        except RuntimeError as error:
            print("ERROR: {}".format(error))
            return -1

    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    command = [bin_to_test, "-e", "namelist.def"]
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    if mpi_procs:
        command = ["mpirun", "-np", mpi_procs] + command
    proc = subprocess.run(
        command,
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
    if expected_error_count is not None and proc.stdout.count(expected) != expected_error_count:
        print(
            "ERROR: expected error substring count {}, observed {}.".format(
                expected_error_count, proc.stdout.count(expected)
            )
        )
        print("---- output begin ----")
        print(proc.stdout)
        print("---- output end ----")
        return -1
    if reject_before_memory:
        forbidden_markers = (
            "End  : Read *def files.",
            "Start: Set memories.",
            "Start: Initialize parameters.",
            "Start: Sampling.",
        )
        observed_markers = [
            marker for marker in forbidden_markers if marker in proc.stdout
        ]
        output_artifacts = sorted(
            name for name in os.listdir(workdir) if name.startswith("zvo")
        )
        if observed_markers or output_artifacts:
            print(
                "ERROR: invalid input was not rejected before memory/output: "
                "markers={} artifacts={}".format(
                    observed_markers, output_artifacts
                )
            )
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
