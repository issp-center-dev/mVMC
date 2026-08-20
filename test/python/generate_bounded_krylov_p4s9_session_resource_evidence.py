from __future__ import print_function

import argparse
import math
import os
import shutil
import statistics

import generate_bounded_krylov_p4c_evidence as p4c


POLICY_SCHEMA_VERSION = 1
DRIVER_SCHEMA_VERSION = 3
EXPECTED_POLICY_ID = "power-lanczos-zero-support-p4-r0-v10"
EXPECTED_FIXTURE = "p4s9_long_direct_session_timing"


def parse_fields(line):
    return {
        key: p4c.parse_number(value)
        for key, value in p4c.parse_key_values(line, 1).items()
    }


def require_close(actual, expected, label, tolerance=1.0e-12):
    actual = float(actual)
    expected = float(expected)
    if not math.isfinite(actual) or not math.isfinite(expected) or \
            abs(actual - expected) > tolerance * max(1.0, abs(expected)):
        raise AssertionError("{} mismatch: {} != {}".format(
            label, actual, expected))


def validate_policy(policy):
    if int(policy.get("schema_version", 0)) != POLICY_SCHEMA_VERSION or \
            policy.get("policy_id") != EXPECTED_POLICY_ID:
        raise AssertionError("P4-S9 resource policy identity mismatch")
    scope = policy.get("scope", {})
    if scope.get("site_counts") != [4, 6, 8] or \
            int(scope.get("target_site_count", 0)) != 16 or \
            scope.get("qp_total") != [1, 4] or \
            scope.get("physical_case_order") != [
                [4, 1], [4, 4], [6, 1], [6, 4], [8, 1], [8, 4]] or \
            int(scope.get("cache_bytes", 0)) != 268435456 or \
            int(scope.get("omp_threads", 0)) != 1 or \
            int(scope.get("blas_threads", 0)) != 1 or \
            scope.get("local_rank_counts") != [1] or \
            scope.get("official_rank_counts") != [1, 2, 4]:
        raise AssertionError("P4-S9 resource scope mismatch")
    driver = policy.get("driver", {})
    if driver.get("mode") != "long-direct-session" or \
            int(driver.get("schema_version", 0)) != DRIVER_SCHEMA_VERSION or \
            driver.get("fixture") != EXPECTED_FIXTURE or \
            int(driver.get("sample_count", 0)) != 32768 or \
            int(driver.get("repeat_count", 0)) != 7 or \
            str(driver.get("seed_hex", "")).lower() != \
            "0x5034533954494d45" or \
            driver.get("stream_order") != list(range(7)) or \
            float(driver.get("rho", 0.0)) != 0.01 or \
            not driver.get("cold_initialization_included_once") or \
            not driver.get("session_end_included") or \
            not driver.get("exact_table_acceptance_audit"):
        raise AssertionError("P4-S9 resource driver mismatch")
    gate = policy.get("resource_gate", {})
    if gate.get("site_extrapolation") != \
            "log(value)=intercept+slope*L_from_L4_L6_L8" or \
            gate.get("rank_reducer") != "maximum" or \
            gate.get("repeat_reducer") != "median" or \
            float(gate.get("seconds_per_configuration_l16_maximum", 0.0)) != \
            0.054931640625 or \
            float(gate.get(
                "combined_session_seconds_per_rank_maximum", 0.0)) != 1800.0 or \
            int(gate.get(
                "allocated_capacity_bytes_per_rank_maximum", 0)) != 4294967296 or \
            int(gate.get("peak_rss_bytes_per_rank_maximum", 0)) != 4294967296 or \
            float(gate.get(
                "p4c_cache256_maximum_seconds_per_configuration", 0.0)) != \
            0.16000264079973991 or \
            float(gate.get("required_speedup_minimum", 0.0)) != \
            2.9127591854032655:
        raise AssertionError("P4-S9 resource gate mismatch")
    for field in ("parent_policy_sha256", "stage_a_eligibility_sha256",
                  "stage_a_evidence_sha256", "p4c_resource_sha256"):
        if len(str(policy.get(field, ""))) != 64:
            raise AssertionError("P4-S9 missing frozen hash {}".format(field))
    return policy


def parse_log(path):
    header = None
    repeats = []
    resources = []
    correctness = []
    with open(path, "r") as handle:
        text = handle.read()
    for line in text.splitlines():
        if line.startswith("TIMING_SESSION_RUN "):
            if header is not None:
                raise AssertionError("duplicate P4-S9 timing header")
            header = parse_fields(line)
        elif line.startswith("TIMING_SESSION_REPEAT "):
            repeats.append(parse_fields(line))
        elif line.startswith("TIMING_SESSION_RESOURCE "):
            resources.append(parse_fields(line))
        elif line.startswith("TIMING_SESSION_CORRECTNESS "):
            correctness.append(parse_fields(line))
    if header is None or len(repeats) != 7 or len(resources) != 7 or \
            len(correctness) != 7:
        raise AssertionError("incomplete P4-S9 timing log {}".format(path))
    return {
        "path": os.path.abspath(path), "sha256": p4c.sha256_file(path),
        "line_count": len(text.splitlines()), "header": header,
        "repeats": repeats, "resources": resources,
        "correctness": correctness,
    }


def validate_run(run, policy, allow_smoke_policy=False):
    header = run["header"]
    scope = policy["scope"]
    driver = policy["driver"]
    site = int(header["site_count"])
    qp = int(header["qp_total"])
    sample_count = int(header["sample_count"])
    cache_bytes = int(header["cache_requested_bytes"])
    if int(header.get("schema", 0)) != DRIVER_SCHEMA_VERSION or \
            header.get("fixture") != EXPECTED_FIXTURE or \
            [site, qp] not in scope["physical_case_order"] or \
            int(header.get("rank_count", 0)) != 1 or \
            int(header.get("repeat_count", 0)) != 7 or \
            int(header.get("sector_size", 0)) != p4c.fixed_sector_size(site) or \
            int(header.get("plan_hash", 0)) == 0 or \
            (not allow_smoke_policy and sample_count !=
             int(driver["sample_count"])) or \
            (not allow_smoke_policy and cache_bytes !=
             int(scope["cache_bytes"])):
        raise AssertionError("P4-S9 timing header mismatch")
    expected_streams = set(driver["stream_order"])
    collections = (run["repeats"], run["resources"], run["correctness"])
    for collection in collections:
        if {int(item["rng_stream"]) for item in collection} != \
                expected_streams:
            raise AssertionError("P4-S9 timing stream census mismatch")
    generation_hashes = set()
    by_stream = {}
    for repeat in run["repeats"]:
        stream = int(repeat["rng_stream"])
        if int(repeat["schema"]) != DRIVER_SCHEMA_VERSION or \
                int(repeat["site_count"]) != site or \
                int(repeat["qp_total"]) != qp or \
                int(repeat["sample_count"]) != sample_count or \
                int(repeat["cache_requested_bytes"]) != cache_bytes or \
                int(repeat["fraction_index"]) != 0 or \
                int(repeat["rho_index"]) != 0 or \
                int(repeat["global_numerator"]) != 0 or \
                int(repeat["global_denominator"]) != 1 or \
                str(repeat["seed_hex"]).lower() != \
                str(driver["seed_hex"]).lower() or \
                int(repeat["accepted"]) + int(repeat["rejected"]) != \
                sample_count or \
                int(repeat["neighbor_attempted"]) != sample_count or \
                int(repeat["global_attempted"]) != 0 or \
                int(repeat["global_self"]) != 0 or \
                int(repeat["proposal_policy_hash"]) == 0 or \
                int(repeat["proposal_model_hash"]) == 0:
            raise AssertionError("P4-S9 timing repeat mismatch")
        require_close(repeat["rho"], 0.01, "P4-S9 timing rho")
        total = float(repeat["total_step_seconds"])
        require_close(repeat["total_seconds_per_step"],
                      total / sample_count, "P4-S9 timing per-step")
        if not math.isfinite(total) or total <= 0.0:
            raise AssertionError("P4-S9 nonpositive timing")
        by_stream[stream] = {"repeat": repeat}
    for resource in run["resources"]:
        stream = int(resource["rng_stream"])
        generation = int(resource["amplitude_generation_hash"])
        combined = float(resource["combined_session_seconds"])
        per_configuration = float(resource["seconds_per_configuration"])
        if int(resource["schema"]) != DRIVER_SCHEMA_VERSION or \
                int(resource["site_count"]) != site or \
                int(resource["qp_total"]) != qp or \
                int(resource["sample_count"]) != sample_count or \
                int(resource["cache_requested_bytes"]) != cache_bytes or \
                generation == 0 or \
                int(resource["session_root_evaluations"]) != \
                sample_count + 1 or \
                int(resource["cache_resets"]) != 1 or \
                int(resource["cache_hits"]) <= 0 or \
                int(resource["cache_misses"]) <= 0 or \
                int(resource["cache_insertions"]) <= 0 or \
                int(resource["cache_evictions"]) > \
                int(resource["cache_insertions"]) or \
                int(resource["cache_resident_peak"]) <= 0 or \
                int(resource["callback_calls"]) <= 0 or \
                int(resource["allocated_capacity_bytes"]) <= 0 or \
                int(resource["peak_rss_bytes"]) <= 0 or \
                float(resource["session_begin_seconds"]) < 0.0 or \
                float(resource["initialization_seconds"]) < 0.0 or \
                not math.isfinite(combined) or combined <= 0.0:
            raise AssertionError("P4-S9 session resource mismatch")
        require_close(per_configuration, combined / sample_count,
                      "P4-S9 session seconds/config")
        generation_hashes.add(generation)
        by_stream[stream]["resource"] = resource
    if len(generation_hashes) != 1:
        raise AssertionError("P4-S9 amplitude generation changed by stream")
    for item in run["correctness"]:
        stream = int(item["rng_stream"])
        if int(item["schema"]) != DRIVER_SCHEMA_VERSION or \
                int(item["site_count"]) != site or \
                int(item["qp_total"]) != qp or \
                int(item["exact_table_match"]) != 1 or \
                int(item["final_sector_index"]) < 0 or \
                int(item["final_sector_index"]) >= int(header["sector_size"]) or \
                int(item["final_configuration"]) <= 0 or \
                int(item["statistical_trace_hash"]) == 0:
            raise AssertionError("P4-S9 correctness mismatch")
        by_stream[stream]["correctness"] = item
    if any(set(item) != {"repeat", "resource", "correctness"}
           for item in by_stream.values()):
        raise AssertionError("P4-S9 timing record join mismatch")
    return {
        "site_count": site, "qp_total": qp, "rank_count": 1,
        "sample_count": sample_count, "cache_requested_bytes": cache_bytes,
        "sector_size": int(header["sector_size"]),
        "plan_hash": int(header["plan_hash"]),
        "amplitude_generation_hash": next(iter(generation_hashes)),
        "records": [by_stream[stream] for stream in driver["stream_order"]],
    }


def materialize_logs(parsed, output):
    root = os.path.join(output, "raw_timing")
    if not os.path.isdir(root):
        os.makedirs(root)
    mapped = {}
    for run in parsed:
        name = "L{}_QP{}-{}-{}".format(
            int(run["header"]["site_count"]),
            int(run["header"]["qp_total"]), run["sha256"][:12],
            os.path.basename(run["path"]))
        destination = os.path.join(root, name)
        shutil.copyfile(run["path"], destination)
        if p4c.sha256_file(destination) != run["sha256"]:
            raise AssertionError("P4-S9 timing log copy mismatch")
        mapped[run["path"]] = os.path.relpath(destination, output)
    return mapped


def p4c_target_allocation(resource, qp):
    key = "cache268435456_qp{}_rank1".format(qp)
    item = resource.get("target_extrapolation", {}).get(key)
    if item is None or int(item.get("cache_requested_bytes", 0)) != 268435456:
        raise AssertionError("missing P4-C allocation reference {}".format(key))
    return int(item["allocated_capacity_bytes_per_rank_l16"])


def analyze_command(args):
    policy = validate_policy(p4c.read_json(args.policy))
    stage_a_eligibility = p4c.read_json(args.stage_a_eligibility)
    p4c_resource = p4c.read_json(args.p4c_resource)
    bindings = {
        "parent_policy_sha256": p4c.sha256_file(args.parent_policy),
        "stage_a_eligibility_sha256": p4c.sha256_file(
            args.stage_a_eligibility),
        "stage_a_evidence_sha256": p4c.sha256_file(args.stage_a_evidence),
        "p4c_resource_sha256": p4c.sha256_file(args.p4c_resource),
    }
    if not args.allow_smoke_policy:
        for field, digest in bindings.items():
            if policy[field] != digest:
                raise AssertionError("P4-S9 frozen input mismatch {}".format(
                    field))
        if stage_a_eligibility.get("stage_a_decision") != "ELIGIBLE" or \
                not stage_a_eligibility.get("stage_b_authorized"):
            raise AssertionError("P4-S9 Stage A did not authorize resource work")
    parsed = [parse_log(path) for path in args.timing_log]
    expected = policy["scope"]["physical_case_order"]
    observed = [[int(run["header"]["site_count"]),
                 int(run["header"]["qp_total"])] for run in parsed]
    if sorted(observed) != sorted(expected) or \
            len(set(tuple(item) for item in observed)) != len(expected):
        raise AssertionError("P4-S9 timing requires six unique cases")
    validated = {
        (int(run["header"]["site_count"]),
         int(run["header"]["qp_total"])):
        validate_run(run, policy, args.allow_smoke_policy) for run in parsed
    }
    output = os.path.abspath(args.output)
    if not os.path.isdir(output):
        os.makedirs(output)
    raw_paths = materialize_logs(parsed, output)
    cases = []
    for site, qp in expected:
        run = next(item for item in parsed
                   if int(item["header"]["site_count"]) == site and
                   int(item["header"]["qp_total"]) == qp)
        case = validated[(site, qp)]
        case["raw_log"] = raw_paths[run["path"]]
        case["raw_log_sha256"] = run["sha256"]
        case["raw_log_line_count"] = run["line_count"]
        case["repeat_seconds_per_configuration"] = [
            float(item["resource"]["seconds_per_configuration"])
            for item in case["records"]]
        case["median_seconds_per_configuration"] = statistics.median(
            case["repeat_seconds_per_configuration"])
        case["maximum_allocated_capacity_bytes"] = max(
            int(item["resource"]["allocated_capacity_bytes"])
            for item in case["records"])
        case["maximum_peak_rss_bytes"] = max(
            int(item["resource"]["peak_rss_bytes"])
            for item in case["records"])
        cases.append(case)
    gate = policy["resource_gate"]
    qp_results = []
    all_pass = True
    for qp in policy["scope"]["qp_total"]:
        runtime_points = [(item["site_count"],
                           item["median_seconds_per_configuration"])
                          for item in cases if item["qp_total"] == qp]
        rss_points = [(item["site_count"], item["maximum_peak_rss_bytes"])
                      for item in cases if item["qp_total"] == qp]
        allocated_points = [
            (item["site_count"], item["maximum_allocated_capacity_bytes"])
            for item in cases if item["qp_total"] == qp]
        runtime_fit = p4c.log_linear_fit(
            runtime_points, policy["scope"]["target_site_count"])
        rss_fit = p4c.log_linear_fit(
            rss_points, policy["scope"]["target_site_count"])
        allocated_fit = p4c.log_linear_fit(
            allocated_points, policy["scope"]["target_site_count"])
        allocated_reference = (max(value for _, value in allocated_points)
                               if args.allow_smoke_policy else
                               p4c_target_allocation(p4c_resource, qp))
        allocated_l16 = max(allocated_reference,
                            allocated_fit["prediction"])
        combined = runtime_fit["prediction"] * int(policy["driver"][
            "sample_count"])
        speedup = float(gate[
            "p4c_cache256_maximum_seconds_per_configuration"]) / \
            runtime_fit["prediction"]
        passed = (
            runtime_fit["prediction"] <= float(gate[
                "seconds_per_configuration_l16_maximum"]) and
            combined <= float(gate[
                "combined_session_seconds_per_rank_maximum"]) and
            allocated_l16 <= int(gate[
                "allocated_capacity_bytes_per_rank_maximum"]) and
            rss_fit["prediction"] <= int(gate[
                "peak_rss_bytes_per_rank_maximum"]) and
            speedup >= float(gate["required_speedup_minimum"]))
        all_pass = all_pass and passed
        qp_results.append({
            "qp_total": qp,
            "seconds_per_configuration_l16": runtime_fit,
            "combined_session_seconds_per_rank_l16": combined,
            "peak_rss_bytes_per_rank_l16": rss_fit,
            "allocated_capacity_bytes_per_rank_l16_fit": allocated_fit,
            "allocated_capacity_bytes_per_rank_l16_reference":
                allocated_reference,
            "allocated_capacity_bytes_per_rank_l16_gate_value":
                allocated_l16,
            "speedup_over_p4c_cache256_maximum": speedup,
            "resource_pass": passed,
        })
    generated_at = p4c.generated_at_jst()
    evidence = {
        "schema_version": POLICY_SCHEMA_VERSION,
        "artifact": "p4s9_local_session_resource_evidence",
        "generated_at_jst": generated_at,
        "source_commit": args.source_commit,
        "producer_binary": os.path.basename(args.producer_binary),
        "producer_binary_sha256": p4c.sha256_file(args.producer_binary),
        "resource_policy_sha256": p4c.sha256_file(args.policy),
        "frozen_input_sha256": bindings,
        "candidate_id": policy["candidate_id"],
        "local_prefilter_pass": all_pass,
        "physical_cases": cases,
        "qp_extrapolations": qp_results,
    }
    eligibility = {
        "schema_version": POLICY_SCHEMA_VERSION,
        "artifact": "p4s9_frozen_local_session_resource_decision",
        "generated_at_jst": generated_at,
        "source_commit": args.source_commit,
        "resource_policy_sha256": evidence["resource_policy_sha256"],
        "candidate_id": policy["candidate_id"],
        "local_resource_decision": "PASS" if all_pass else "STOP",
        "remote_package_preparation_authorized": all_pass,
        "remote_transfer_or_submission_authorized": False,
        "failure_action": (
            "prepare_checksummed_remote_package_and_request_user_approval"
            if all_pass else "formal_local_stop_do_not_prepare_remote_package"),
    }
    p4c.write_json(os.path.join(
        output, "p4s9_local_session_resource_evidence.json"), evidence)
    p4c.write_json(os.path.join(
        output, "p4s9_local_session_resource_decision.json"), eligibility)
    count = p4c.write_checksums(output)
    print("P4-S9 local session resource: {}; {} ledger files".format(
        eligibility["local_resource_decision"], count))


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--policy", required=True)
    parser.add_argument("--parent-policy", required=True)
    parser.add_argument("--stage-a-eligibility", required=True)
    parser.add_argument("--stage-a-evidence", required=True)
    parser.add_argument("--p4c-resource", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--timing-log", nargs="+", required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--producer-binary", required=True)
    parser.add_argument("--allow-smoke-policy", action="store_true")
    return parser


def main(argv=None):
    analyze_command(build_parser().parse_args(argv))


if __name__ == "__main__":
    main()
