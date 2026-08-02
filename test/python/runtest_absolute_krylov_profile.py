from __future__ import print_function

import math
import os
import shutil
import subprocess
import sys


def command(driver, audit):
    mpi_procs = int(os.environ.get("MVMC_MPI_PROCS", "1"))
    arguments = [driver, "4", "1", "6", str(audit)]
    if mpi_procs == 1:
        return arguments
    mpiexec = os.environ.get("MVMC_MPIEXEC") or shutil.which("mpiexec")
    if not mpiexec:
        raise AssertionError("MVMC_MPIEXEC/mpiexec is unavailable")
    flag = os.environ.get("MVMC_MPIEXEC_NUMPROC_FLAG", "-np")
    return [mpiexec, flag, str(mpi_procs)] + arguments


def run(driver, audit):
    process = subprocess.run(
        command(driver, audit),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        check=False,
    )
    if process.returncode != 0:
        raise AssertionError(
            "profile driver audit={} failed:\n{}".format(
                audit, process.stdout
            )
        )
    return process.stdout.splitlines()


def parse_key_values(line):
    fields = {}
    for token in line.split()[1:]:
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        fields[key] = value
    return fields


def validate_resource(line):
    fields = parse_key_values(line)
    required = (
        "scope", "roots", "recursion", "unique", "memo_hits",
        "memo_misses", "raw_transitions", "terminal_requests",
        "regular_logical", "near_pivot_logical", "singular_logical",
        "total_zero_logical", "local_factorizations",
        "global_factorizations_logical", "frontier_peak", "memo_peak",
        "krylov_workspace_bytes", "model_workspace_bytes",
        "amplitude_workspace_bytes", "rss_bytes", "depth_seconds",
        "amplitude_seconds", "connectivity_seconds", "total_seconds",
    )
    missing = [key for key in required if key not in fields]
    if missing:
        raise AssertionError("resource fields missing: {}".format(missing))
    if fields["scope"] not in ("rank_sum", "rank_max"):
        raise AssertionError("unexpected resource scope")
    for key in (
        "roots", "raw_transitions", "terminal_requests",
        "global_factorizations_logical", "frontier_peak", "memo_peak",
        "krylov_workspace_bytes", "model_workspace_bytes",
        "amplitude_workspace_bytes", "rss_bytes",
    ):
        if int(fields[key]) <= 0:
            raise AssertionError("nonpositive resource {}".format(key))
    for key in (
        "amplitude_seconds", "connectivity_seconds", "total_seconds",
    ):
        value = float(fields[key])
        if not math.isfinite(value) or value < 0.0:
            raise AssertionError("invalid timing {}".format(key))
    depth = [float(value) for value in fields["depth_seconds"].split(",")]
    if len(depth) != 4 or any(
            not math.isfinite(value) or value < 0.0 for value in depth):
        raise AssertionError("invalid depth timings")
    if (float(fields["total_seconds"]) <= 0.0 or sum(depth) <= 0.0 or
            float(fields["amplitude_seconds"]) +
            float(fields["connectivity_seconds"]) <= 0.0):
        raise AssertionError("profile hook recorded no measurable work")
    return fields


def main():
    if len(sys.argv) != 2:
        raise AssertionError("usage: runtest_absolute_krylov_profile.py DRIVER")
    without_audit = run(sys.argv[1], 0)
    with_audit = run(sys.argv[1], 1)
    rows_without = [line for line in without_audit if line.startswith("ROW ")]
    rows_with = [line for line in with_audit if line.startswith("ROW ")]
    if rows_without != rows_with or len(rows_with) != 6:
        raise AssertionError("audit changed deterministic configuration/value trace")
    if any(line.startswith(("STAT ", "RESOURCE ")) for line in without_audit):
        raise AssertionError("audit-off run leaked profile records")
    stats = [line for line in with_audit if line.startswith("STAT ")]
    resources = [line for line in with_audit if line.startswith("RESOURCE ")]
    if len(stats) != 6 or len(resources) != 2:
        raise AssertionError("incomplete audit records")
    parsed = {
        fields["scope"]: fields
        for fields in (validate_resource(line) for line in resources)
    }
    logical = (
        "regular_logical", "near_pivot_logical", "singular_logical",
        "total_zero_logical", "global_factorizations_logical",
    )
    for key in logical:
        if int(parsed["rank_sum"][key]) != int(parsed["rank_max"][key]):
            raise AssertionError(
                "communicator-global counter was rank-multiplied: {}".format(
                    key))
    if (int(parsed["rank_sum"]["local_factorizations"]) !=
            int(parsed["rank_sum"]["global_factorizations_logical"])):
        raise AssertionError("local factorization work count is inconsistent")
    print("absolute Krylov profile audit parity: PASS")


if __name__ == "__main__":
    main()
