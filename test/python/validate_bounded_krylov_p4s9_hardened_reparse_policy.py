from __future__ import print_function

import argparse
import os
import re

import generate_bounded_krylov_p4c_evidence as p4c


POLICY_SCHEMA_VERSION = 1
EXPECTED_POLICY_ID = \
    "power-lanczos-zero-support-p4-r0-v13-hardened-reparse"
EXPECTED_CANDIDATE_ID = "direct_neighbor_rho0p01_N32768_session256m"
EXPECTED_PARENT_PATH = \
    "test/python/data/power_lanczos_p4s9_official_execution_policy.json"
EXPECTED_EXECUTION_CONTRACT = {
    "execution_parameter_changes": False,
    "seed_changes": False,
    "sample_count_changes": False,
    "partition_changes": False,
    "rank_changes": False,
    "cache_changes": False,
    "threshold_changes": False,
    "source_commit": "39d48c272998521d78166d38acf60c6b2e8624e5",
    "baseline_develop_commit":
        "39d48c272998521d78166d38acf60c6b2e8624e5",
    "source_archive_sha256":
        "a006c0f3c3b33dc9d794a263818ff5d217116ea17ed41ab1be48458636388495",
    "transfer_package_sha256":
        "ceec699498c885b42ea1f78d9f6ba5c5338002963fe938a326aa33eaa1bf5da8",
    "statistics_log_count": 24,
    "timing_log_count": 18,
    "target_log_count": 18,
    "mpi_parity_log_count": 3,
}
EXPECTED_INVARIANTS = {
    "statistics": [
        "decimal_seed_matches_seed_hex",
        "tau_37_record_exact_census_and_gate_reapplication",
        "entry_18_record_exact_census",
        "denominator_summary_recomputed_from_samples",
        "support_tau_denominator_decision_consistency",
    ],
    "resource": [
        "canonical_input_paths_are_unique",
        "timing_plan_generation_and_stream_trace_match_rank_1_2_4",
        "target_generation_and_row_trace_match_rank_1_2_4",
    ],
    "closure": [
        "component_ledger_matches_complete_regular_file_census",
        "component_evidence_and_decision_are_ledger_bound",
        "symlinks_and_duplicate_ledger_entries_are_rejected",
        "component_schema_artifact_candidate_and_commit_provenance_match",
        "focused_ctest_required_results_each_pass_exactly_once",
        "focused_ctest_has_no_failed_result_and_summary_count_matches",
    ],
}
EXPECTED_PROMOTION = {
    "primary_v12_closure_required": "GO",
    "hardened_v13_reparse_closure_required": "GO",
    "raw_log_sha256_census_must_match": True,
    "case_numeric_results_must_match": True,
    "any_mismatch_action": "STOP_PRESERVE_EVIDENCE_NO_P4F",
}
EXPECTED_AUTHORIZATION = {
    "p4f_authorized": False,
    "remote_mutation_authorized": False,
    "push_or_github_mutation_authorized": False,
}


def checked_source_path(source_root, relative):
    source_root = os.path.realpath(source_root)
    if os.path.isabs(relative) or os.path.normpath(relative) != relative or \
            ".." in relative.split(os.sep):
        raise AssertionError("unsafe hardened parser manifest path")
    lexical = os.path.join(source_root, relative)
    resolved = os.path.realpath(lexical)
    if os.path.commonpath([source_root, resolved]) != source_root or \
            os.path.islink(lexical) or not os.path.isfile(resolved):
        raise AssertionError(
            "missing hardened parser manifest file {}".format(relative))
    return resolved


def validate_policy(policy, source_root):
    if int(policy.get("schema_version", 0)) != POLICY_SCHEMA_VERSION or \
            policy.get("policy_id") != EXPECTED_POLICY_ID or \
            policy.get("candidate_id") != EXPECTED_CANDIDATE_ID or \
            policy.get("registered_before_result_inspection") is not True or \
            re.fullmatch(
                r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\+09:00",
                str(policy.get("registered_at_jst", ""))) is None:
        raise AssertionError("P4-S9 hardened reparse policy identity mismatch")
    if policy.get("execution_contract") != EXPECTED_EXECUTION_CONTRACT or \
            policy.get("hardened_invariants") != EXPECTED_INVARIANTS or \
            policy.get("promotion") != EXPECTED_PROMOTION or \
            policy.get("authorization") != EXPECTED_AUTHORIZATION:
        raise AssertionError("P4-S9 hardened reparse contract mismatch")

    parent = policy.get("parent_official_execution_policy", {})
    if set(parent) != {"path", "sha256"} or \
            parent.get("path") != EXPECTED_PARENT_PATH or \
            re.fullmatch(r"[0-9a-f]{64}", str(parent.get("sha256", ""))) \
            is None:
        raise AssertionError("P4-S9 hardened parent policy mismatch")
    parent_path = checked_source_path(source_root, parent["path"])
    if p4c.sha256_file(parent_path) != parent["sha256"]:
        raise AssertionError("P4-S9 hardened parent policy hash mismatch")

    manifest = policy.get("parser_manifest_sha256", {})
    expected_paths = {
        "test/python/generate_bounded_krylov_p4c_evidence.py",
        "test/python/generate_bounded_krylov_p4s_evidence.py",
        "test/python/generate_bounded_krylov_p4s9_official_statistics_evidence.py",
        "test/python/generate_bounded_krylov_p4s9_session_resource_evidence.py",
        "test/python/generate_bounded_krylov_p4s9_target_session_evidence.py",
        "test/python/generate_bounded_krylov_p4s9_official_resource_evidence.py",
        "test/python/generate_bounded_krylov_p4s9_official_closure.py",
        "test/python/validate_bounded_krylov_p4s9_hardened_reparse_policy.py",
        "test/python/generate_bounded_krylov_p4s9_hardened_comparison.py",
    }
    if set(manifest) != expected_paths:
        raise AssertionError("P4-S9 hardened parser manifest census mismatch")
    for relative, expected_hash in manifest.items():
        if re.fullmatch(r"[0-9a-f]{64}", str(expected_hash)) is None or \
                p4c.sha256_file(checked_source_path(
                    source_root, relative)) != expected_hash:
            raise AssertionError(
                "P4-S9 hardened parser hash mismatch {}".format(relative))
    return policy


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("policy")
    parser.add_argument("source_root")
    args = parser.parse_args(argv)
    validate_policy(p4c.read_json(args.policy), args.source_root)
    print("P4-S9 hardened reparse policy: PASS")


if __name__ == "__main__":
    main()
