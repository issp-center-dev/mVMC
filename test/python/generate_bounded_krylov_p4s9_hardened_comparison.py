from __future__ import print_function

import argparse
import json
import os
import re

import generate_bounded_krylov_p4c_evidence as p4c
import generate_bounded_krylov_p4s9_official_closure as closure
import generate_bounded_krylov_p4s9_official_statistics_evidence as official
import validate_bounded_krylov_p4s9_hardened_reparse_policy as validator


def canonical_value(value):
    if isinstance(value, dict):
        items = [
            [canonical_value(key), canonical_value(item)]
            for key, item in value.items()
        ]
        items.sort(key=lambda item: json.dumps(
            item[0], allow_nan=False, separators=(",", ":"),
            sort_keys=True))
        return {"mapping": items}
    if isinstance(value, tuple):
        return {"tuple": [canonical_value(item) for item in value]}
    if isinstance(value, list):
        return [canonical_value(item) for item in value]
    return value


def canonical_digest(value):
    text = json.dumps(
        canonical_value(value), allow_nan=False, separators=(",", ":"),
        sort_keys=True)
    return p4c.sha256_text(text)


def required_fields(value, names, label):
    missing = [name for name in names if name not in value]
    if missing:
        raise AssertionError("{} missing fields {}".format(label, missing))
    return {name: value[name] for name in names}


def component(root, name, evidence_name, decision_name):
    component_root = os.path.realpath(os.path.join(root, name))
    if os.path.commonpath([os.path.realpath(root), component_root]) != \
            os.path.realpath(root) or not os.path.isdir(component_root):
        raise AssertionError("missing hardened component {}".format(name))
    ledger = closure.verify_ledger(component_root)
    evidence_path = os.path.join(component_root, evidence_name)
    decision_path = os.path.join(component_root, decision_name)
    for path in (evidence_path, decision_path):
        relative = os.path.relpath(path, component_root)
        if ledger["by_path"].get(relative) != p4c.sha256_file(path):
            raise AssertionError("hardened component JSON not ledger-bound")
    return {
        "root": component_root,
        "ledger": ledger,
        "evidence": p4c.read_json(evidence_path),
        "decision": p4c.read_json(decision_path),
    }


def raw_binding(component_value, item, path_field="raw_log",
                hash_field="raw_log_sha256"):
    relative = item.get(path_field)
    expected_hash = item.get(hash_field)
    if not isinstance(relative, str) or not isinstance(expected_hash, str):
        raise AssertionError("missing hardened raw binding")
    lexical = os.path.join(component_value["root"], relative)
    resolved = os.path.realpath(lexical)
    if os.path.isabs(relative) or os.path.normpath(relative) != relative or \
            ".." in relative.split(os.sep) or \
            os.path.commonpath([component_value["root"], resolved]) != \
            component_value["root"] or os.path.islink(lexical) or \
            not os.path.isfile(resolved) or \
            component_value["ledger"]["by_path"].get(relative) != \
            expected_hash or p4c.sha256_file(resolved) != expected_hash:
        raise AssertionError("hardened raw binding mismatch {}".format(
            relative))
    if "raw_log_line_count" in item:
        with open(resolved, "r") as handle:
            actual_line_count = sum(1 for _ in handle)
        if int(item["raw_log_line_count"]) != actual_line_count:
            raise AssertionError(
                "hardened raw line count mismatch {}".format(relative))
    return expected_hash


def load_triplet(root):
    root = os.path.realpath(root)
    return {
        "statistics": component(
            root, "statistics", "p4s9_official_markov_statistics.json",
            "p4s9_official_statistics_decision.json"),
        "resource": component(
            root, "resource", "p4s9_official_resource_evidence.json",
            "p4s9_official_resource_decision.json"),
        "closure": component(
            root, "closure", "p4s9_official_closure_evidence.json",
            "p4s9_official_closure_decision.json"),
    }


def statistics_canonical(component_value, execution_policy):
    fields = (
        "site_count", "qp_total", "seed_hex", "sample_count",
        "rank_count", "amplitude_generation_hash", "raw_log_sha256",
        "raw_log_line_count", "trace", "summary", "session", "manifest",
        "decision", "maximum_tau_int", "maximum_se_budget_ratio",
        "official_statistics_pass",
    )
    records = {}
    hashes = []
    for item in component_value["evidence"].get("cases", []):
        key = (int(item["site_count"]), int(item["qp_total"]),
               str(item["seed_hex"]).lower())
        if key in records:
            raise AssertionError("duplicate hardened statistics case")
        hashes.append(raw_binding(component_value, item))
        records[key] = required_fields(item, fields, "statistics case")
    statistics = execution_policy["official_statistics"]
    expected = {
        (site, qp, seed.lower())
        for site, qp in statistics["physical_case_order"]
        for seed in statistics["seed_hex_order"]
    }
    if set(records) != expected or len(hashes) != 24:
        raise AssertionError("hardened statistics census mismatch")
    return records, sorted(hashes)


def resource_canonical(component_value, execution_policy):
    evidence = component_value["evidence"]
    timing_case_fields = (
        "rank_count", "site_count", "qp_total", "sample_count",
        "median_seconds_per_configuration",
        "maximum_allocated_capacity_bytes_per_rank",
        "maximum_peak_rss_bytes_per_rank", "raw_log_sha256",
        "raw_log_line_count",
    )
    timing_cases = {}
    timing_hashes = []
    for item in evidence.get("session_timing_cases", []):
        key = (int(item["rank_count"]), int(item["site_count"]),
               int(item["qp_total"]))
        if key in timing_cases:
            raise AssertionError("duplicate hardened timing case")
        timing_hashes.append(raw_binding(component_value, item))
        timing_cases[key] = required_fields(
            item, timing_case_fields, "timing case")
    timing = execution_policy["official_resource"]["session_timing"]
    timing_expected = {
        (rank, site, qp)
        for rank in timing["rank_counts"]
        for site in timing["site_counts"]
        for qp in timing["qp_total"]
    }
    if set(timing_cases) != timing_expected or len(timing_hashes) != 18:
        raise AssertionError("hardened timing census mismatch")

    timing_result_fields = (
        "rank_count", "qp_total", "seconds_per_configuration_l16",
        "projected_32768_root_seconds",
        "allocated_capacity_bytes_per_rank_l16",
        "peak_rss_bytes_per_rank_l16",
        "speedup_over_p4c_cache256_maximum", "resource_pass",
    )
    timing_results = {}
    for item in evidence.get("session_timing_results", []):
        key = (int(item["rank_count"]), int(item["qp_total"]))
        if key in timing_results:
            raise AssertionError("duplicate hardened timing result")
        timing_results[key] = required_fields(
            item, timing_result_fields, "timing result")
    expected_results = {
        (rank, qp) for rank in timing["rank_counts"]
        for qp in timing["qp_total"]
    }
    if set(timing_results) != expected_results:
        raise AssertionError("hardened timing result census mismatch")

    target_repeat_fields = (
        "rank_count", "qp_total", "sample_count", "cache_requested_bytes",
        "row_trace_sha256", "amplitude_generation_hash", "total_seconds",
        "total_seconds_per_root", "allocated_capacity_bytes_per_rank",
        "peak_rss_bytes_per_rank", "cache_hits", "cache_misses",
        "cache_insertions", "cache_evictions", "cache_entries_peak",
        "terminal_amplitude_calls", "raw_log_sha256", "raw_log_line_count",
    )
    target_repeats = {}
    target_hashes = []
    for item in evidence.get("direct_target_repeats", []):
        raw_hash = raw_binding(component_value, item)
        key = (int(item["rank_count"]), int(item["qp_total"]))
        target_hashes.append(raw_hash)
        target_repeats.setdefault(key, []).append(required_fields(
            item, target_repeat_fields, "target repeat"))
    target = execution_policy["official_resource"]["direct_target_stress"]
    for rank in target["rank_counts"]:
        for qp in target["qp_total"]:
            key = (rank, qp)
            if len(target_repeats.get(key, [])) != \
                    int(target["repeat_count"]):
                raise AssertionError("hardened target repeat census mismatch")
            target_repeats[key].sort(key=canonical_digest)
    if len(target_hashes) != 18:
        raise AssertionError("hardened target log count mismatch")

    target_result_fields = (
        "rank_count", "qp_total", "repeat_total_seconds_per_root",
        "median_total_seconds_per_root", "projected_32768_root_seconds",
        "maximum_allocated_capacity_bytes_per_rank",
        "maximum_peak_rss_bytes_per_rank",
        "speedup_over_p4c_cache256_maximum", "row_trace_sha256",
        "trace_invariance_pass", "resource_pass",
    )
    target_results = {}
    for item in evidence.get("direct_target_results", []):
        key = (int(item["rank_count"]), int(item["qp_total"]))
        if key in target_results:
            raise AssertionError("duplicate hardened target result")
        target_results[key] = required_fields(
            item, target_result_fields, "target result")
    target_expected = {
        (rank, qp) for rank in target["rank_counts"]
        for qp in target["qp_total"]
    }
    if set(target_results) != target_expected:
        raise AssertionError("hardened target result census mismatch")

    passes = required_fields(evidence, (
        "official_session_timing_pass", "official_direct_target_pass",
        "official_resource_pass"), "resource evidence")
    return {
        "timing_cases": timing_cases,
        "timing_results": timing_results,
        "target_repeats": target_repeats,
        "target_results": target_results,
        "passes": passes,
    }, sorted(timing_hashes), sorted(target_hashes)


def closure_canonical(component_value, execution_policy):
    evidence = component_value["evidence"]
    runs = {}
    hashes = []
    for item in evidence.get("mpi_parity_runs", []):
        key = int(item["rank_count"])
        if key in runs:
            raise AssertionError("duplicate hardened parity run")
        hashes.append(raw_binding(component_value, item))
        runs[key] = required_fields(item, (
            "rank_count", "raw_log_sha256", "math_trace_sha256"),
            "parity run")
    if set(runs) != set(execution_policy["mpi_parity"]["rank_counts"]) or \
            len(hashes) != 3:
        raise AssertionError("hardened parity census mismatch")
    focused = evidence.get("focused_ctest", {})
    ctest_hash = raw_binding(
        component_value, focused, "raw_log", "sha256")
    focused_shared = required_fields(focused, (
        "sha256", "required_test_count", "missing_tests",
        "pass_summary_present", "focused_ctest_pass"), "focused CTest")
    shared = required_fields(evidence, (
        "statistics_pass", "resource_pass", "mpi_parity_trace_sha256",
        "mpi_parity_pass", "official_closure_pass"), "closure evidence")
    return {
        "mpi_parity_runs": runs,
        "focused_ctest": focused_shared,
        "shared": shared,
    }, sorted(hashes), ctest_hash


def component_provenance(component_value):
    evidence = required_fields(component_value["evidence"], (
        "source_commit", "baseline_develop_commit",
        "official_execution_policy_sha256", "candidate_id"),
        "component provenance")
    decision = required_fields(component_value["decision"], (
        "source_commit", "baseline_develop_commit",
        "official_execution_policy_sha256", "candidate_id"),
        "component decision provenance")
    if evidence != decision:
        raise AssertionError("component evidence/decision provenance mismatch")
    return evidence


def validate_triplet_identity(triplet, execution_policy_hash, candidate_id,
                              expected_source_commit,
                              expected_baseline_commit):
    artifacts = {
        "statistics": (
            "p4s9_official_markov_statistics",
            "p4s9_official_statistics_decision"),
        "resource": (
            "p4s9_official_resource_evidence",
            "p4s9_official_resource_decision"),
        "closure": (
            "p4s9_official_closure_evidence",
            "p4s9_official_closure_decision"),
    }
    provenances = []
    for name, expected_artifacts in artifacts.items():
        evidence = triplet[name]["evidence"]
        decision = triplet[name]["decision"]
        if int(evidence.get("schema_version", 0)) != 1 or \
                int(decision.get("schema_version", 0)) != 1 or \
                evidence.get("artifact") != expected_artifacts[0] or \
                decision.get("artifact") != expected_artifacts[1]:
            raise AssertionError("hardened component identity mismatch")
        provenance = component_provenance(triplet[name])
        if provenance["official_execution_policy_sha256"] != \
                execution_policy_hash or \
                provenance["candidate_id"] != candidate_id or \
                provenance["source_commit"] != expected_source_commit or \
                provenance["baseline_develop_commit"] != \
                expected_baseline_commit or \
                re.fullmatch(r"[0-9a-f]{40}", expected_source_commit) is \
                None or re.fullmatch(
                    r"[0-9a-f]{40}", expected_baseline_commit) is None:
            raise AssertionError("hardened component provenance mismatch")
        provenances.append(provenance)
    if any(item != provenances[0] for item in provenances[1:]):
        raise AssertionError("hardened cross-component provenance mismatch")

    statistics_evidence = triplet["statistics"]["evidence"]
    statistics_decision = triplet["statistics"]["decision"]
    statistics_pass = bool(statistics_evidence.get(
        "official_statistics_pass"))
    if (statistics_decision.get("official_statistics_decision") == "GO") \
            != statistics_pass or \
            bool(statistics_decision.get("p4f_authorized", True)):
        raise AssertionError("hardened statistics component mismatch")
    resource_evidence = triplet["resource"]["evidence"]
    resource_decision = triplet["resource"]["decision"]
    resource_pass = bool(resource_evidence.get("official_resource_pass"))
    if (resource_decision.get("official_resource_decision") == "PASS") \
            != resource_pass or \
            bool(resource_decision.get("p4f_authorized", True)):
        raise AssertionError("hardened resource component mismatch")
    closure_evidence = triplet["closure"]["evidence"]
    closure_decision = triplet["closure"]["decision"]
    closure_pass = bool(closure_evidence.get("official_closure_pass"))
    if (closure_decision.get("p4s9_decision") == "GO") != closure_pass or \
            bool(closure_decision.get("p4f_authorized")) != closure_pass:
        raise AssertionError("hardened closure component mismatch")
    return provenances[0]


def analyze_command(args):
    source_root = os.path.realpath(args.source_root)
    hardened_policy = validator.validate_policy(
        p4c.read_json(args.hardened_policy), source_root)
    parent_path = os.path.join(
        source_root,
        hardened_policy["parent_official_execution_policy"]["path"])
    execution_policy = official.validate_policy(p4c.read_json(parent_path))
    execution_policy_hash = p4c.sha256_file(parent_path)
    execution_contract = hardened_policy["execution_contract"]
    primary = load_triplet(args.primary_root)
    reparse = load_triplet(args.reparse_root)
    primary_provenance = validate_triplet_identity(
        primary, execution_policy_hash, execution_policy["candidate_id"],
        execution_contract["source_commit"],
        execution_contract["baseline_develop_commit"])
    reparse_provenance = validate_triplet_identity(
        reparse, execution_policy_hash, execution_policy["candidate_id"],
        execution_contract["source_commit"],
        execution_contract["baseline_develop_commit"])

    primary_statistics, primary_statistics_hashes = statistics_canonical(
        primary["statistics"], execution_policy)
    reparse_statistics, reparse_statistics_hashes = statistics_canonical(
        reparse["statistics"], execution_policy)
    primary_resource, primary_timing_hashes, primary_target_hashes = \
        resource_canonical(primary["resource"], execution_policy)
    reparse_resource, reparse_timing_hashes, reparse_target_hashes = \
        resource_canonical(reparse["resource"], execution_policy)
    primary_closure, primary_parity_hashes, primary_ctest_hash = \
        closure_canonical(primary["closure"], execution_policy)
    reparse_closure, reparse_parity_hashes, reparse_ctest_hash = \
        closure_canonical(reparse["closure"], execution_policy)

    provenance_match = primary_provenance == reparse_provenance
    raw_sha_match = (
        primary_statistics_hashes == reparse_statistics_hashes and
        primary_timing_hashes == reparse_timing_hashes and
        primary_target_hashes == reparse_target_hashes and
        primary_parity_hashes == reparse_parity_hashes and
        primary_ctest_hash == reparse_ctest_hash)
    numeric_match = (
        primary_statistics == reparse_statistics and
        primary_resource == reparse_resource and
        primary_closure == reparse_closure)
    primary_decision = primary["closure"]["decision"]
    reparse_decision = reparse["closure"]["decision"]
    decision_match = required_fields(primary_decision, (
        "p4s9_decision", "p4f_authorized"), "primary closure decision") == \
        required_fields(reparse_decision, (
            "p4s9_decision", "p4f_authorized"),
            "reparse closure decision")
    double_go = (
        primary_decision.get("p4s9_decision") == "GO" and
        bool(primary_decision.get("p4f_authorized")) and
        reparse_decision.get("p4s9_decision") == "GO" and
        bool(reparse_decision.get("p4f_authorized")))
    hardened_pass = (provenance_match and raw_sha_match and numeric_match and
                     decision_match and double_go and
                     bool(reparse["resource"]["evidence"].get(
                         "official_session_timing_cross_rank_trace_pass")) and
                     bool(reparse["resource"]["evidence"].get(
                         "official_direct_target_cross_rank_trace_pass")))

    output = os.path.abspath(args.output)
    if not os.path.isdir(output):
        os.makedirs(output)
    generated_at = p4c.generated_at_jst()
    evidence = {
        "schema_version": 1,
        "artifact": "p4s9_hardened_reparse_comparison_evidence",
        "generated_at_jst": generated_at,
        "candidate_id": hardened_policy["candidate_id"],
        "hardened_reparse_policy_sha256":
            p4c.sha256_file(args.hardened_policy),
        "primary_component_ledger_sha256": {
            name: primary[name]["ledger"]["ledger_sha256"]
            for name in ("statistics", "resource", "closure")
        },
        "reparse_component_ledger_sha256": {
            name: reparse[name]["ledger"]["ledger_sha256"]
            for name in ("statistics", "resource", "closure")
        },
        "raw_log_census": {
            "statistics": len(primary_statistics_hashes),
            "timing": len(primary_timing_hashes),
            "target": len(primary_target_hashes),
            "mpi_parity": len(primary_parity_hashes),
            "focused_ctest": 1,
        },
        "primary_canonical_sha256": {
            "statistics": canonical_digest(primary_statistics),
            "resource": canonical_digest(primary_resource),
            "closure": canonical_digest(primary_closure),
        },
        "reparse_canonical_sha256": {
            "statistics": canonical_digest(reparse_statistics),
            "resource": canonical_digest(reparse_resource),
            "closure": canonical_digest(reparse_closure),
        },
        "provenance_match": provenance_match,
        "raw_log_sha256_census_match": raw_sha_match,
        "case_numeric_results_match": numeric_match,
        "closure_decision_match": decision_match,
        "primary_and_reparse_double_go": double_go,
        "hardened_comparison_pass": hardened_pass,
    }
    decision = {
        "schema_version": 1,
        "artifact": "p4s9_hardened_reparse_comparison_decision",
        "generated_at_jst": generated_at,
        "candidate_id": hardened_policy["candidate_id"],
        "hardened_reparse_policy_sha256":
            p4c.sha256_file(args.hardened_policy),
        "hardened_v13_reparse_decision": "GO" if hardened_pass else "STOP",
        "p4f_authorized": hardened_pass,
        "push_or_github_mutation_authorized": False,
        "failure_action": "begin_local_p4f"
            if hardened_pass else "preserve_evidence_and_stop_before_p4f",
    }
    p4c.write_json(os.path.join(
        output, "p4s9_hardened_reparse_comparison_evidence.json"), evidence)
    p4c.write_json(os.path.join(
        output, "p4s9_hardened_reparse_comparison_decision.json"), decision)
    count = p4c.write_checksums(output)
    print("P4-S9 hardened comparison: {}; {} ledger files".format(
        decision["hardened_v13_reparse_decision"], count))


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hardened-policy", required=True)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--primary-root", required=True)
    parser.add_argument("--reparse-root", required=True)
    parser.add_argument("--output", required=True)
    return parser


def main(argv=None):
    analyze_command(build_parser().parse_args(argv))


if __name__ == "__main__":
    main()
