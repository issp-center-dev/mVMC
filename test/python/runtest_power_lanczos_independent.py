from __future__ import print_function

import os
import pathlib
import shutil
import subprocess
import sys

import numpy as np


root = pathlib.Path.cwd()
source = root / "data" / "BackFlow_Identity_InterAll_Real"
extra = root / "data" / "PowerLanczosIndependent"
arguments = set(sys.argv[1:])
unknown_arguments = arguments.difference(("blocked", "interall", "invalid"))
if unknown_arguments:
    raise RuntimeError("unknown arguments: {}".format(sorted(unknown_arguments)))
processes = os.environ.get("MVMC_MPI_PROCS")
blocked = "blocked" in arguments
invalid = "invalid" in arguments
interall = "interall" in arguments
suffix_parts = []
if processes:
    suffix_parts.append("mpi")
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
text = output.read_text()
if "# independent_power_lanczos" not in text:
    print(text)
    raise RuntimeError("independent output header is missing")
rows = [line for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")]
if len(rows) != 1:
    raise RuntimeError("expected one result row, got {}".format(len(rows)))
values = np.fromstring(rows[0], sep=" ")
if values.size != 12 or not np.all(np.isfinite(values)):
    raise RuntimeError("invalid independent result: {}".format(values))
if int(round(values[8])) < 1:
    raise RuntimeError("invalid retained rank: {}".format(values[8]))
if processes:
    serial_output = root / "work" / "PowerLanczosIndependent" / "output" / "zvo_pl_out_001.dat"
    serial_rows = [line for line in serial_output.read_text().splitlines()
                   if line.strip() and not line.lstrip().startswith("#")]
    serial_values = np.fromstring(serial_rows[0], sep=" ")
    if not np.allclose(values, serial_values, rtol=2e-12, atol=2e-12):
        raise RuntimeError("rank-1/rank-2 mismatch:\n{}\n{}".format(
            serial_values, values))
print("PowerLanczosIndependent{}: PASS".format(suffix))
