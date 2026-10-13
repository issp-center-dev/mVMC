"""Verify chain-root composite-zero counting and world reduction."""
import os
import shutil
import subprocess
import sys

from gc_component_stats import parse_component_stats


binary = sys.argv[1]
launcher = os.environ.get("MVMC_MPIEXEC") or shutil.which("mpiexec")
numproc_flag = os.environ.get("MVMC_MPIEXEC_NUMPROC_FLAG", "-np")
if launcher is None:
    raise SystemExit("MPI launcher is unavailable")
command = [launcher, numproc_flag, "2", binary,
           "--component-stats-collective"]
process = subprocess.run(command, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         universal_newlines=True, timeout=30)
print(process.stdout, end="")
if process.returncode != 0:
    raise SystemExit("component stats MPI run failed")
stats = parse_component_stats(process.stdout)
if stats["composite_zero"] != 2:
    raise SystemExit("composite_zero={} instead of 2".format(
        stats["composite_zero"]))
print("GC component stats MPI reduction PASS")
