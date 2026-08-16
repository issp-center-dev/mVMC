from __future__ import print_function

import argparse
import os
import shutil
import sys
import tempfile

import generate_bounded_krylov_p4s9_hardened_comparison as comparison


COMMIT = "39d48c272998521d78166d38acf60c6b2e8624e5"


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative hardened comparison accepted: {}".format(
        label))


def raw_file(root, relative, text):
    path = os.path.join(root, relative)
    parent = os.path.dirname(path)
    if not os.path.isdir(parent):
        os.makedirs(parent)
    with open(path, "w") as handle:
        handle.write(text)
    return comparison.p4c.sha256_file(path)


def provenance(execution_policy_hash, candidate_id):
    return {
        "schema_version": 1,
        "source_commit": COMMIT,
        "baseline_develop_commit": COMMIT,
        "official_execution_policy_sha256": execution_policy_hash,
        "candidate_id": candidate_id,
    }


def build_statistics(root, policy, execution_policy_hash):
    component = os.path.join(root, "statistics")
    os.makedirs(component)
    cases = []
    for site, qp in policy["official_statistics"]["physical_case_order"]:
        for seed in policy["official_statistics"]["seed_hex_order"]:
            relative = "raw/stats-L{}-QP{}-{}.log".format(site, qp, seed)
            digest = raw_file(
                component, relative,
                "statistics {} {} {}\n".format(site, qp, seed))
            cases.append({
                "site_count": site,
                "qp_total": qp,
                "seed_hex": seed,
                "sample_count": 32768,
                "rank_count": 1,
                "amplitude_generation_hash": 123,
                "raw_log": relative,
                "raw_log_sha256": digest,
                "raw_log_line_count": 1,
                "trace": {"completed": 32768},
                "summary": {"sample_count": 32768},
                "session": {"session_end_pass": 1},
                "manifest": {"rng_stream": 0},
                "decision": {"p4s_decision": "GO"},
                "maximum_tau_int": 2.0,
                "maximum_se_budget_ratio": 0.5,
                "official_statistics_pass": True,
            })
    evidence = provenance(execution_policy_hash, policy["candidate_id"])
    evidence.update({
        "artifact": "p4s9_official_markov_statistics",
        "cases": cases,
        "official_statistics_pass": True,
    })
    decision = provenance(execution_policy_hash, policy["candidate_id"])
    decision.update({
        "artifact": "p4s9_official_statistics_decision",
        "official_statistics_decision": "GO",
        "p4f_authorized": False,
    })
    comparison.p4c.write_json(os.path.join(
        component, "p4s9_official_markov_statistics.json"), evidence)
    comparison.p4c.write_json(os.path.join(
        component, "p4s9_official_statistics_decision.json"), decision)
    comparison.p4c.write_checksums(component)


def build_resource(root, policy, execution_policy_hash):
    component = os.path.join(root, "resource")
    os.makedirs(component)
    resource = policy["official_resource"]
    timing_cases = []
    for rank in resource["session_timing"]["rank_counts"]:
        for site in resource["session_timing"]["site_counts"]:
            for qp in resource["session_timing"]["qp_total"]:
                relative = "raw/timing-R{}-L{}-QP{}.log".format(
                    rank, site, qp)
                digest = raw_file(
                    component, relative,
                    "timing {} {} {}\n".format(rank, site, qp))
                timing_cases.append({
                    "rank_count": rank,
                    "site_count": site,
                    "qp_total": qp,
                    "sample_count": 32768,
                    "median_seconds_per_configuration": 0.01,
                    "maximum_allocated_capacity_bytes_per_rank": 1000,
                    "maximum_peak_rss_bytes_per_rank": 2000,
                    "raw_log": relative,
                    "raw_log_sha256": digest,
                    "raw_log_line_count": 1,
                })
    timing_results = []
    for rank in resource["session_timing"]["rank_counts"]:
        for qp in resource["session_timing"]["qp_total"]:
            timing_results.append({
                "rank_count": rank,
                "qp_total": qp,
                "seconds_per_configuration_l16": {"prediction": 0.01},
                "projected_32768_root_seconds": 327.68,
                "allocated_capacity_bytes_per_rank_l16": 1000,
                "peak_rss_bytes_per_rank_l16": 2000,
                "speedup_over_p4c_cache256_maximum": 16.0,
                "resource_pass": True,
            })
    target_repeats = []
    for rank in resource["direct_target_stress"]["rank_counts"]:
        for qp in resource["direct_target_stress"]["qp_total"]:
            for repeat in range(3):
                relative = "raw/target-R{}-QP{}-N{}.log".format(
                    rank, qp, repeat)
                digest = raw_file(
                    component, relative,
                    "target {} {}\n".format(rank, qp))
                target_repeats.append({
                    "rank_count": rank,
                    "qp_total": qp,
                    "sample_count": 4096,
                    "cache_requested_bytes": 268435456,
                    "row_trace_sha256": "1" * 64,
                    "amplitude_generation_hash": 123,
                    "total_seconds": 40.96,
                    "total_seconds_per_root": 0.01,
                    "allocated_capacity_bytes_per_rank": 1000,
                    "peak_rss_bytes_per_rank": 2000,
                    "cache_hits": [1, 1, 1, 1],
                    "cache_misses": [1, 1, 1, 1],
                    "cache_insertions": 1,
                    "cache_evictions": 1,
                    "cache_entries_peak": 1,
                    "terminal_amplitude_calls": 1,
                    "raw_log": relative,
                    "raw_log_sha256": digest,
                    "raw_log_line_count": 1,
                })
    target_results = []
    for rank in resource["direct_target_stress"]["rank_counts"]:
        for qp in resource["direct_target_stress"]["qp_total"]:
            target_results.append({
                "rank_count": rank,
                "qp_total": qp,
                "repeat_total_seconds_per_root": [0.01, 0.01, 0.01],
                "median_total_seconds_per_root": 0.01,
                "projected_32768_root_seconds": 327.68,
                "maximum_allocated_capacity_bytes_per_rank": 1000,
                "maximum_peak_rss_bytes_per_rank": 2000,
                "speedup_over_p4c_cache256_maximum": 16.0,
                "row_trace_sha256": "1" * 64,
                "trace_invariance_pass": True,
                "resource_pass": True,
            })
    evidence = provenance(execution_policy_hash, policy["candidate_id"])
    evidence.update({
        "artifact": "p4s9_official_resource_evidence",
        "session_timing_cases": timing_cases,
        "session_timing_results": timing_results,
        "direct_target_repeats": target_repeats,
        "direct_target_results": target_results,
        "official_session_timing_cross_rank_trace_pass": True,
        "official_direct_target_cross_rank_trace_pass": True,
        "official_session_timing_pass": True,
        "official_direct_target_pass": True,
        "official_resource_pass": True,
    })
    decision = provenance(execution_policy_hash, policy["candidate_id"])
    decision.update({
        "artifact": "p4s9_official_resource_decision",
        "official_resource_decision": "PASS",
        "p4f_authorized": False,
    })
    comparison.p4c.write_json(os.path.join(
        component, "p4s9_official_resource_evidence.json"), evidence)
    comparison.p4c.write_json(os.path.join(
        component, "p4s9_official_resource_decision.json"), decision)
    comparison.p4c.write_checksums(component)


def build_closure(root, policy, execution_policy_hash):
    component = os.path.join(root, "closure")
    os.makedirs(component)
    parity_runs = []
    for rank in policy["mpi_parity"]["rank_counts"]:
        relative = "raw/parity-R{}.log".format(rank)
        digest = raw_file(component, relative, "parity {}\n".format(rank))
        parity_runs.append({
            "rank_count": rank,
            "raw_log": relative,
            "raw_log_sha256": digest,
            "math_trace_sha256": "2" * 64,
        })
    ctest_relative = "workflow/focused-ctest.log"
    ctest_hash = raw_file(component, ctest_relative, "100% tests passed\n")
    evidence = provenance(execution_policy_hash, policy["candidate_id"])
    evidence.update({
        "artifact": "p4s9_official_closure_evidence",
        "statistics_pass": True,
        "resource_pass": True,
        "mpi_parity_runs": parity_runs,
        "mpi_parity_trace_sha256": "2" * 64,
        "mpi_parity_pass": True,
        "focused_ctest": {
            "raw_log": ctest_relative,
            "sha256": ctest_hash,
            "required_test_count": 13,
            "missing_tests": [],
            "pass_summary_present": True,
            "focused_ctest_pass": True,
        },
        "official_closure_pass": True,
    })
    decision = provenance(execution_policy_hash, policy["candidate_id"])
    decision.update({
        "artifact": "p4s9_official_closure_decision",
        "p4s9_decision": "GO",
        "p4f_authorized": True,
    })
    comparison.p4c.write_json(os.path.join(
        component, "p4s9_official_closure_evidence.json"), evidence)
    comparison.p4c.write_json(os.path.join(
        component, "p4s9_official_closure_decision.json"), decision)
    comparison.p4c.write_checksums(component)


def build_triplet(root, policy, execution_policy_hash):
    os.makedirs(root)
    build_statistics(root, policy, execution_policy_hash)
    build_resource(root, policy, execution_policy_hash)
    build_closure(root, policy, execution_policy_hash)


def make_args(hardened_policy, source_root, primary, reparse, output):
    return argparse.Namespace(
        hardened_policy=hardened_policy,
        source_root=source_root,
        primary_root=primary,
        reparse_root=reparse,
        output=output,
    )


def rewrite_provenance(triplet_root, key, value):
    component_files = {
        "statistics": (
            "p4s9_official_markov_statistics.json",
            "p4s9_official_statistics_decision.json"),
        "resource": (
            "p4s9_official_resource_evidence.json",
            "p4s9_official_resource_decision.json"),
        "closure": (
            "p4s9_official_closure_evidence.json",
            "p4s9_official_closure_decision.json"),
    }
    for name, filenames in component_files.items():
        component_root = os.path.join(triplet_root, name)
        for filename in filenames:
            path = os.path.join(component_root, filename)
            document = comparison.p4c.read_json(path)
            document[key] = value
            comparison.p4c.write_json(path, document)
        comparison.p4c.write_checksums(component_root)


def main():
    if len(sys.argv) != 3:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4s9_hardened_comparison.py "
            "HARDENED_POLICY SOURCE_ROOT")
    hardened_policy = os.path.abspath(sys.argv[1])
    source_root = os.path.abspath(sys.argv[2])
    hardened = comparison.p4c.read_json(hardened_policy)
    parent_path = os.path.join(
        source_root, hardened["parent_official_execution_policy"]["path"])
    policy = comparison.official.validate_policy(
        comparison.p4c.read_json(parent_path))
    policy_hash = comparison.p4c.sha256_file(parent_path)
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4s9-comparison-")
    try:
        primary = os.path.join(root, "primary")
        reparse = os.path.join(root, "reparse")
        build_triplet(primary, policy, policy_hash)
        shutil.copytree(primary, reparse)
        output = os.path.join(root, "comparison-go")
        comparison.analyze_command(make_args(
            hardened_policy, source_root, primary, reparse, output))
        decision = comparison.p4c.read_json(os.path.join(
            output, "p4s9_hardened_reparse_comparison_decision.json"))
        if decision["hardened_v13_reparse_decision"] != "GO" or \
                not decision["p4f_authorized"]:
            raise AssertionError("synthetic hardened comparison GO mismatch")

        resource_path = os.path.join(
            reparse, "resource", "p4s9_official_resource_evidence.json")
        changed = comparison.p4c.read_json(resource_path)
        changed["session_timing_results"][0][
            "projected_32768_root_seconds"] += 1.0
        comparison.p4c.write_json(resource_path, changed)
        comparison.p4c.write_checksums(os.path.join(reparse, "resource"))
        stop_output = os.path.join(root, "comparison-stop")
        comparison.analyze_command(make_args(
            hardened_policy, source_root, primary, reparse, stop_output))
        stop_decision = comparison.p4c.read_json(os.path.join(
            stop_output,
            "p4s9_hardened_reparse_comparison_decision.json"))
        if stop_decision["hardened_v13_reparse_decision"] != "STOP" or \
                stop_decision["p4f_authorized"]:
            raise AssertionError(
                "hardened numeric comparison mismatch did not STOP")

        wrong_primary = os.path.join(root, "wrong-primary")
        wrong_reparse = os.path.join(root, "wrong-reparse")
        shutil.copytree(primary, wrong_primary)
        shutil.copytree(primary, wrong_reparse)
        rewrite_provenance(
            wrong_primary, "official_execution_policy_sha256", "f" * 64)
        rewrite_provenance(
            wrong_reparse, "official_execution_policy_sha256", "f" * 64)
        expect_failure(lambda: comparison.analyze_command(make_args(
            hardened_policy, source_root, wrong_primary, wrong_reparse,
            os.path.join(root, "comparison-wrong-provenance"))),
            "matching but incorrect parent policy provenance")

        wrong_commit_primary = os.path.join(root, "wrong-commit-primary")
        wrong_commit_reparse = os.path.join(root, "wrong-commit-reparse")
        shutil.copytree(primary, wrong_commit_primary)
        shutil.copytree(primary, wrong_commit_reparse)
        rewrite_provenance(
            wrong_commit_primary, "source_commit", "f" * 40)
        rewrite_provenance(
            wrong_commit_reparse, "source_commit", "f" * 40)
        expect_failure(lambda: comparison.analyze_command(make_args(
            hardened_policy, source_root, wrong_commit_primary,
            wrong_commit_reparse,
            os.path.join(root, "comparison-wrong-source-commit"))),
            "matching but incorrect source commit provenance")
        print("bounded Krylov P4-S9 hardened comparison: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
