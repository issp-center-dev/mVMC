from __future__ import print_function

import filecmp
import os
import shutil
import subprocess
import sys


def prepare_definitions(source, destination, dry_binary):
    os.makedirs(destination)
    for filename in os.listdir(source):
        source_path = os.path.join(source, filename)
        if os.path.isfile(source_path) and (
            filename.endswith(".def") or filename.endswith(".dat")
        ):
            shutil.copy(
                source_path,
                os.path.join(destination, filename),
            )
    modpara = os.path.join(destination, "modpara.def")
    if not os.path.exists(modpara):
        standard_input = os.path.join(destination, "StdFace.def")
        proc = subprocess.run(
            [dry_binary, standard_input],
            cwd=destination,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if proc.returncode != 0:
            print(proc.stdout)
            raise RuntimeError("failed to generate expert definition files")


def set_lanczos_step(filename, explicit):
    with open(filename) as source:
        lines = source.readlines()
    lines = [
        line
        for line in lines
        if line.split()[:1]
        not in (["NLanczosStep"], ["NLanczosEstimatorMode"])
    ]
    # This Expert-mode test isolates the NLanczosStep default.  Select the
    # legacy estimator explicitly so an omitted estimator does not enter the
    # corrected route, whose Expert-mode default is intentionally 0.
    lines.append("NLanczosEstimatorMode   1\n")
    if explicit:
        lines.append("NLanczosStep   1\n")
    with open(filename, "w") as destination:
        destination.writelines(lines)


def comparable_outputs(workdir):
    output_dir = os.path.join(workdir, "output")
    result = {}
    for filename in os.listdir(output_dir):
        if "time" in filename.lower() or "calctimer" in filename.lower():
            continue
        path = os.path.join(output_dir, filename)
        if os.path.isfile(path):
            result[filename] = path
    return result


def run_case(binary, workdir):
    proc = subprocess.run(
        [binary, "-e", "namelist.def"],
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    with open(os.path.join(workdir, "run.log"), "w") as destination:
        destination.write(proc.stdout)
    if proc.returncode != 0:
        print(proc.stdout)
    return proc.returncode


def main():
    if len(sys.argv) != 2:
        print("usage: {} <model name>".format(sys.argv[0]))
        return 2

    model = sys.argv[1]
    rootdir = os.getcwd()
    source = os.path.join(rootdir, "data", model)
    workroot = os.path.join(rootdir, "work")
    omitted = os.path.join(workroot, model + "_lanczos_step_omitted")
    explicit = os.path.join(workroot, model + "_lanczos_step_explicit_1")
    binary = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    dry_binary = os.path.join(
        rootdir, "..", "..", "src", "mVMC", "vmcdry.out"
    )
    for workdir in (omitted, explicit):
        if os.path.exists(workdir):
            shutil.rmtree(workdir)
        try:
            prepare_definitions(source, workdir, dry_binary)
        except RuntimeError as error:
            print("ERROR: {}".format(error))
            return 1

    set_lanczos_step(os.path.join(omitted, "modpara.def"), False)
    set_lanczos_step(os.path.join(explicit, "modpara.def"), True)

    if run_case(binary, omitted) != 0 or run_case(binary, explicit) != 0:
        return 1

    omitted_outputs = comparable_outputs(omitted)
    explicit_outputs = comparable_outputs(explicit)
    if not omitted_outputs:
        print("ERROR: compatibility fixture produced no comparable output")
        return 1
    required_lanczos_outputs = ("_ls_out_", "_ls_qqqq_")
    for marker in required_lanczos_outputs:
        if not any(marker in filename for filename in omitted_outputs):
            print(
                "ERROR: compatibility fixture did not produce a {} file".format(
                    marker.strip("_")
                )
            )
            return 1
    if set(omitted_outputs) != set(explicit_outputs):
        print("ERROR: output file sets differ")
        print("omitted:", sorted(omitted_outputs))
        print("explicit:", sorted(explicit_outputs))
        return 1

    for filename in sorted(omitted_outputs):
        if not filecmp.cmp(
            omitted_outputs[filename], explicit_outputs[filename], shallow=False
        ):
            print("ERROR: output differs for {}".format(filename))
            return 1

    print(
        "Lanczos2 step-1 compatibility: PASS ({} byte-identical files)".format(
            len(omitted_outputs)
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
