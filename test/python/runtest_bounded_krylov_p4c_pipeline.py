from __future__ import print_function

import argparse
import copy
import os
import re
import shutil
import subprocess
import sys
import tempfile

import generate_bounded_krylov_p4c_evidence as evidence


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative gate was accepted: {}".format(label))


def git_head(repository):
    fallback = os.environ.get("MVMC_P4C_SOURCE_COMMIT")
    if fallback and re.match(r"^[0-9a-f]{40}$", fallback):
        return fallback
    process = subprocess.run(
        ["git", "-C", repository, "rev-parse", "HEAD"],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True, check=False,
    )
    if process.returncode != 0:
        raise AssertionError("cannot resolve test repository HEAD")
    return process.stdout.strip()


def write_profile(path, site, qp, ranks, sample_count, cache_bytes,
                  total_scale, constant_rows=True):
    sector_size = evidence.fixed_sector_size(site)
    if sample_count == 0:
        sample_count = sector_size
    if sample_count > sector_size:
        raise AssertionError("synthetic sample count exceeds sector")
    with open(path, "w") as handle:
        handle.write(
            "PROFILE schema=1 fixture={} site_count={} qp_total={} "
            "sample_count={} sector_size={} rank_count={} "
            "cache_requested_bytes={} audit=1\n".format(
                evidence.EXPECTED_PROFILE_FIXTURE, site, qp, sample_count,
                sector_size, ranks, cache_bytes))
        for sample, sector_index, configuration in \
                evidence.expected_profile_rows(site, sample_count):
            values = [1.0, 2.0, 3.0, 4.0]
            if not constant_rows and sample == 0:
                values[3] = 5.0
            handle.write(
                "ROW {} {} {} {} 0 {} 0 {} 0 {} 0\n".format(
                    sample, sector_index, configuration,
                    values[0], values[1], values[2], values[3]))
            handle.write(
                "STAT {} reset_seconds=0 total_seconds={} "
                "depth_seconds={},{},{},{} terminal_calls=4 "
                "raw_transitions=8 cache_hits=1,1,1,1 "
                "cache_misses=1,1,1,1 row_peak=8 "
                "cache_entries_peak=4 workspace_bytes=2048 "
                "cache_allocated_bytes={}\n".format(
                    sample, total_scale, total_scale / 4.0,
                    total_scale / 4.0, total_scale / 4.0,
                    total_scale / 4.0, cache_bytes))
        for scope in ("rank_sum", "rank_max"):
            multiplier = ranks if scope == "rank_sum" else 1
            roots = sample_count * multiplier
            total = total_scale * roots
            aggregated_cache_bytes = cache_bytes * multiplier
            handle.write(
                "RESOURCE scope={} roots={} node_expansions={} "
                "recursion=1,1,1,1 cache_hits=1,1,1,1 "
                "cache_misses=1,1,1,1 cache_insertions={} "
                "cache_evictions=0 cache_epoch_clears=0 "
                "raw_transitions={} merged_duplicates=0 "
                "cancelled_zero=0 terminal_calls={} well_pivoted={} "
                "near_pivot=0 exact_zero_components=0 "
                "numeric_zero_components=0 local_factorizations={} "
                "global_factorizations={} computed_exact_zero=0 "
                "computed_numeric_zero=0 trace_hash=12345 row_peak=8 "
                "cache_entries_peak=4 plan_bytes=1024 model_bytes=1024 "
                "engine_workspace_bytes={} frame_bytes=1024 "
                "scratch_bytes=64 cache_requested_bytes={} "
                "cache_allocated_bytes={} cache_set_count=1 "
                "collective_workspace_bytes=1024 "
                "amplitude_workspace_bytes=1024 "
                "engine_heap_allocations=0 rss_bytes={} "
                "reset_seconds=0 evaluation_seconds={} "
                "depth_seconds={},{},{},{} amplitude_seconds={} "
                "connectivity_seconds={} cache_seconds={} "
                "total_seconds={}\n".format(
                    scope, roots, roots, roots, roots * 8, roots * 4,
                    roots * 4, roots * 4, roots * 4, 2048 + cache_bytes,
                    aggregated_cache_bytes, aggregated_cache_bytes,
                    2000000 + site * 1000,
                    total, total / 4.0, total / 4.0, total / 4.0,
                    total / 4.0, total / 2.0, total / 4.0, total / 16.0,
                    total))


def analyze_args(policy_path, output, resource_logs, exact_logs,
                 allocation_logs, repository, commit, correctness_gate="PASS"):
    return argparse.Namespace(
        source_commit=commit,
        baseline_commit=commit,
        compiler="synthetic-compiler",
        mpi="synthetic-mpi",
        blas="synthetic-blas",
        strict_fp="enabled",
        omp_threads=1,
        blas_threads=1,
        repository=repository,
        policy=policy_path,
        output=output,
        resource_log=resource_logs,
        exact_log=exact_logs,
        allocation_log=allocation_logs,
        correctness_gate=correctness_gate,
        correctness_note="synthetic focused gate",
        allow_smoke_policy=True,
    )


def main():
    if len(sys.argv) != 3:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4c_pipeline.py POLICY REPOSITORY")
    official_policy_path = os.path.abspath(sys.argv[1])
    repository = os.path.abspath(sys.argv[2])
    official = evidence.validate_policy_contract(
        evidence.read_json(official_policy_path))
    commit = git_head(repository)
    validation_repository = (
        repository if os.path.isdir(os.path.join(repository, ".git"))
        else None)
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4c-pipeline-")
    try:
        smoke_policy = copy.deepcopy(official)
        smoke_policy["scope"]["qp_total"] = [1]
        smoke_policy["scope"]["mpi_rank_counts"] = [1, 2]
        smoke_policy["scope"]["repeat_count"] = 2
        smoke_policy["scope"]["cache_capacity_bytes"] = [0, 4096]
        smoke_policy["scope"]["sample_count_by_site"] = {
            "4": 3, "6": 3, "8": 3,
        }
        policy_path = os.path.join(root, "p4c_smoke_policy.json")
        evidence.write_json(policy_path, smoke_policy)
        resource_logs = []
        exact_logs = []
        allocation_logs = []
        inputs = os.path.join(root, "inputs")
        output = os.path.join(root, "evidence")
        os.makedirs(inputs)
        for cache in smoke_policy["scope"]["cache_capacity_bytes"]:
            for site in smoke_policy["scope"]["site_counts"]:
                for rank in smoke_policy["scope"]["mpi_rank_counts"]:
                    for repeat in range(smoke_policy["scope"]["repeat_count"]):
                        path = os.path.join(
                            inputs,
                            "resource_cache{}_L{}_rank{}_rep{}.log".format(
                                cache, site, rank, repeat))
                        write_profile(path, site, 1, rank, 3, cache,
                                      1.0e-7 * (site / 4.0))
                        resource_logs.append(path)
                exact_path = os.path.join(
                    inputs, "exact_cache{}_L{}.log".format(cache, site))
                write_profile(exact_path, site, 1, 1, 0, cache,
                              1.0e-7 * (site / 4.0))
                exact_logs.append(exact_path)
            for rank in smoke_policy["scope"]["mpi_rank_counts"]:
                allocation_path = os.path.join(
                    inputs,
                    "allocation_cache{}_rank{}.log".format(cache, rank))
                write_profile(
                    allocation_path,
                    smoke_policy["scope"]["target_site_count"], 1, rank, 1,
                    cache, 1.0e-6)
                allocation_logs.append(allocation_path)
        evidence.analyze_command(analyze_args(
            policy_path, output, resource_logs, exact_logs,
            allocation_logs, validation_repository, commit))
        decision = evidence.read_json(
            os.path.join(output, "p4c_candidate_decision.json"))
        if decision["p4c_decision"] != "GO" or not decision["promoted_candidates"]:
            raise AssertionError("synthetic P4-C pipeline did not GO")
        checksum_path = os.path.join(output, "checksums.sha256")
        with open(checksum_path, "r") as handle:
            checksum_lines = [line.rstrip("\n") for line in handle]
        if len(checksum_lines) < len(resource_logs) + len(exact_logs):
            raise AssertionError("raw logs are absent from checksum ledger")
        for line in checksum_lines:
            digest, relative = line.split("  ", 1)
            if evidence.sha256_file(os.path.join(output, relative)) != digest:
                raise AssertionError("final checksum ledger mismatch")
        stop_output = os.path.join(root, "stop_evidence")
        evidence.analyze_command(analyze_args(
            policy_path, stop_output, resource_logs, exact_logs,
            allocation_logs, validation_repository, commit,
            correctness_gate="FAIL"))
        stop_decision = evidence.read_json(
            os.path.join(stop_output, "p4c_candidate_decision.json"))
        if stop_decision["p4c_decision"] != "STOP":
            raise AssertionError("correctness gate failure did not STOP")
        bad_exact = os.path.join(root, "bad_exact.log")
        write_profile(bad_exact, 4, 1, 1, 0, 4096, 1.0e-7,
                      constant_rows=False)
        bad_logs = list(exact_logs)
        bad_logs[3] = bad_exact
        expect_failure(
            lambda: evidence.analyze_command(analyze_args(
                policy_path, os.path.join(root, "bad_evidence"),
                resource_logs, bad_logs, allocation_logs,
                validation_repository, commit)),
            "cache invariance")
        print("bounded Krylov P4-C command pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
