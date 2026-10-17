"""Boundary sampler/Green integration, including MPI death paths."""
import argparse
import os
import subprocess
from gc_antiparallel_failure import reject_signal_reports, require_expected_failure

parser = argparse.ArgumentParser()
parser.add_argument("sampler")
parser.add_argument("green")
parser.add_argument("--mpiexec")
args = parser.parse_args()
commands = [([args.green, "--boundary"], None)]
for ranks in ([1, 2, 3] if args.mpiexec else [1]):
    prefix = [args.mpiexec, "-np", str(ranks)] if ranks > 1 else []
    for mode, n in [("--boundary", None), ("--boundary-vacuum", 0), ("--boundary-full", 4)]:
        if ranks == 3 and n is not None:
            continue
        commands.append((prefix + [args.sampler, mode], n))
env = dict(os.environ, OMPI_MCA_rmaps_base_oversubscribe="1")
for command, n in commands:
    proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          universal_newlines=True, env=env, timeout=45)
    print(" ".join(command), "exit", proc.returncode)
    print(proc.stdout)
    if n is None:
        reject_signal_reports(proc.stdout)
        if proc.returncode != 0:
            raise AssertionError("boundary unit failed")
    else:
        require_expected_failure(proc.returncode, proc.stdout,
                                 ["zero projected amplitude", "after 1 attempt(s)",
                                  "Ncur={}".format(n)])
print("GC boundary integration: PASS")
