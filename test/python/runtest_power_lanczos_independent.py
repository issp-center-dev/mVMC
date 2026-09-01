from __future__ import print_function

import os
import pathlib
import re
import shutil
import subprocess
import sys

import numpy as np


def read_result(output):
    text = output.read_text()
    if "# independent_power_lanczos" not in text:
        print(text)
        raise RuntimeError("independent output header is missing: {}".format(
            output))
    rows = [line for line in text.splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    if len(rows) != 1:
        raise RuntimeError("expected one result row in {}, got {}".format(
            output, len(rows)))
    values = np.fromstring(rows[0], sep=" ")
    if values.size != 12 or not np.all(np.isfinite(values)):
        raise RuntimeError("invalid independent result in {}: {}".format(
            output, values))
    if int(round(values[8])) < 1:
        raise RuntimeError("invalid retained rank in {}: {}".format(
            output, values[8]))
    return values


def read_sampling_metadata(output):
    header = output.read_text().splitlines()[0]
    metadata = {}
    for name in ("chains", "samples_per_chain", "coefficient_samples",
                 "final_samples", "blocks"):
        match = re.search(r"\b{}=(\d+)\b".format(name), header)
        if match is None:
            raise RuntimeError("missing {} metadata in {}".format(
                name, output))
        metadata[name] = int(match.group(1))
    return metadata


root = pathlib.Path.cwd()
source = root / "data" / "BackFlow_Identity_InterAll_Real"
extra = root / "data" / "PowerLanczosIndependent"
arguments = set(sys.argv[1:])
unknown_arguments = arguments.difference(
    ("blocked", "interall", "invalid", "projection", "hybrid"))
if unknown_arguments:
    raise RuntimeError("unknown arguments: {}".format(sorted(unknown_arguments)))
processes = os.environ.get("MVMC_MPI_PROCS")
blocked = "blocked" in arguments
invalid = "invalid" in arguments
interall = "interall" in arguments
projection = "projection" in arguments
hybrid = "hybrid" in arguments
if projection and hybrid:
    raise RuntimeError("projection and hybrid modes are mutually exclusive")
if (projection or hybrid) and not processes:
    raise RuntimeError("MPI split test requires MVMC_MPI_PROCS")
suffix_parts = []
if processes:
    suffix_parts.append("mpi")
if projection:
    suffix_parts.append("projection")
if hybrid:
    suffix_parts.append("hybrid")
if invalid:
    suffix_parts.append("invalid")
elif interall:
    suffix_parts.append("interall")
if blocked:
    suffix_parts.append("blocked")
suffix = "_" + "_".join(suffix_parts) if suffix_parts else ""
work = root / "work" / ("PowerLanczosIndependent" + suffix)
if work.exists():
    shutil.rmtree(str(work))
work.mkdir(parents=True)

for name in ("locspn.def", "trans.def", "interall.def", "orbitalidx.def",
             "qptransidx.def"):
    shutil.copy2(str(source / name), str(work / name))
ls_name = "lsinterall.def" if interall else "lstrans.def"
shutil.copy2(str(extra / ls_name), str(work / ls_name))
if invalid:
    ls_text = (work / "lstrans.def").read_text()
    ls_text = ls_text.replace("-0.30        0.0", "-0.30        0.25")
    (work / "lstrans.def").write_text(ls_text)

modpara = (source / "modpara.def").read_text()
replacements = {
    "NLanczosMode   0": "NLanczosMode   1\nNLanczosStep   1\nNLanczosEstimatorMode 1",
    "NVMCWarmUp     4": "NVMCWarmUp     8",
    "NVMCSample     16": "NVMCSample     256",
}
for old, new in replacements.items():
    if old not in modpara:
        raise RuntimeError("missing modpara anchor: {}".format(old))
    modpara = modpara.replace(old, new)
split_size = int(processes) if projection else (2 if hybrid else 1)
if projection or hybrid:
    split_anchor = "NSplitSize     1"
    if split_anchor not in modpara:
        raise RuntimeError("missing NSplitSize anchor")
    modpara = modpara.replace(
        split_anchor, "NSplitSize     {}".format(split_size))
(work / "modpara.def").write_text(modpara)
namelist_lines = [
    " ModPara modpara.def\n",
    " LocSpin locspn.def\n",
    " Trans trans.def\n",
    " InterAll interall.def\n",
    " {} {}\n".format("LsInterAll" if interall else "LsTrans", ls_name),
    " Orbital orbitalidx.def\n",
    " TransSym qptransidx.def\n",
]
(work / "namelist.def").write_text("".join(namelist_lines))

binary = root.parent.parent / "src" / "mVMC" / "vmc.out"
command = [str(binary), "-e", "namelist.def"]
if processes:
    command = ["mpirun", "-np", processes] + command
completed = subprocess.run(command, cwd=str(work), text=True,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)
if invalid:
    if completed.returncode == 0 or "LsTrans row 0 must be finite, real" not in completed.stdout:
        print(completed.stdout)
        raise RuntimeError("invalid LsTrans was not rejected as expected")
    print("PowerLanczosIndependent_invalid: PASS")
    sys.exit(0)
if blocked:
    expected_error = "independent power-Lanczos failed: unsupported model"
    if completed.returncode == 0 or expected_error not in completed.stdout:
        print(completed.stdout)
        raise RuntimeError("blocked-update build was not rejected as expected")
    output = work / "output" / "zvo_pl_out_001.dat"
    if output.exists():
        rows = [line for line in output.read_text().splitlines()
                if line.strip() and not line.lstrip().startswith("#")]
        if rows:
            raise RuntimeError("blocked-update rejection wrote a result row")
    print("PowerLanczosIndependent{}: PASS".format(suffix))
    sys.exit(0)
if completed.returncode != 0:
    print(completed.stdout)
    sys.exit(completed.returncode)

output = work / "output" / "zvo_pl_out_001.dat"
values = read_result(output)
metadata = read_sampling_metadata(output)
if processes:
    expected_chains = (
        int(processes) + split_size - 1) // split_size
    expected_total = 256 * expected_chains
    expected_metadata = {
        "chains": expected_chains,
        "samples_per_chain": 256,
        "coefficient_samples": expected_total,
        "final_samples": expected_total,
    }
    for name, expected in expected_metadata.items():
        if metadata[name] != expected:
            raise RuntimeError("{} metadata mismatch: {} != {}".format(
                name, metadata[name], expected))
    serial_work = work.with_name(work.name + "_serial_reference")
    if serial_work.exists():
        shutil.rmtree(str(serial_work))
    shutil.copytree(str(work), str(serial_work),
                    ignore=shutil.ignore_patterns("output"))
    serial_completed = subprocess.run(
        [str(binary), "-e", "namelist.def"], cwd=str(serial_work), text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if serial_completed.returncode != 0:
        print(serial_completed.stdout)
        raise RuntimeError("serial reference run failed with status {}".format(
            serial_completed.returncode))
    serial_values = read_result(
        serial_work / "output" / "zvo_pl_out_001.dat")
    serial_metadata = read_sampling_metadata(
        serial_work / "output" / "zvo_pl_out_001.dat")
    if serial_metadata["chains"] != 1 or \
            serial_metadata["coefficient_samples"] != 256:
        raise RuntimeError("invalid serial sampling metadata: {}".format(
            serial_metadata))
    if expected_chains == 1:
        if not np.allclose(values, serial_values, rtol=2e-12, atol=2e-12):
            raise RuntimeError("projection-parallel mismatch:\n{}\n{}".format(
                serial_values, values))
    else:
        for value_index, error_index, label in ((0, 1, "energy"),
                                                 (2, 3, "variance")):
            tolerance = max(
                5e-8,
                8.0 * np.hypot(values[error_index],
                               serial_values[error_index]))
            if abs(values[value_index] - serial_values[value_index]) > tolerance:
                raise RuntimeError(
                    "sample-parallel {} mismatch: {} vs {} (tol {})".format(
                        label, values[value_index],
                        serial_values[value_index], tolerance))
print("PowerLanczosIndependent{}: PASS".format(suffix))
