from __future__ import print_function

import argparse
import math
import os
import statistics

import generate_bounded_krylov_p4c_evidence as p4c
import generate_bounded_krylov_p4s3_pilot_evidence as p4s3


POLICY_SCHEMA_VERSION = 8
LOG_SCHEMA_VERSION = 2
EXPECTED_FIXTURE = "p4s7_partial_callback_matched_timing"
EXPECTED_CANDIDATES = ((2, 1, 32), (2, 10, 32))
read_json = p4c.read_json
write_json = p4c.write_json
sha256_file = p4c.sha256_file


def parse_key_values(line):
    return {key: p4c.parse_number(value)
            for key, value in p4c.parse_key_values(line, 1).items()}


def candidate_key(candidate):
    return "m{}_alpha{}_K{}".format(*candidate)


def require_close(actual, expected, label, tolerance=1.0e-12):
    actual = float(actual)
    expected = float(expected)
    if not math.isfinite(actual) or not math.isfinite(expected) or \
            abs(actual - expected) > tolerance * max(1.0, abs(expected)):
        raise AssertionError("{} mismatch: {} != {}".format(
            label, actual, expected))


def validate_policy(policy, allow_smoke_policy=False):
    if int(policy.get("schema_version", 0)) != POLICY_SCHEMA_VERSION or \
            policy.get("policy_id") != \
            "power-lanczos-zero-support-p4-r0-v7-stage-c":
        raise AssertionError("P4-S7 timing policy identity mismatch")
    expected_cases = [[4, 1], [4, 4], [6, 1], [6, 4], [8, 1], [8, 4]]
    scope = policy.get("scope", {})
    if scope.get("physical_case_order") != expected_cases or \
            scope.get("site_counts") != [4, 6, 8] or \
            scope.get("qp_total") != [1, 4]:
        raise AssertionError("P4-S7 timing scope mismatch")
    observed = tuple((int(x.get("m", 0)), int(x.get("alpha", 0)),
                      int(x.get("K", 0)))
                     for x in policy.get("candidate_order", []))
    timing = policy.get("timing_census", {})
    if observed != EXPECTED_CANDIDATES or \
            int(timing.get("sample_count", 0)) != 128 or \
            int(timing.get("repeat_count", 0)) != 7 or \
            timing.get("seed_hex", "").lower() != \
            "0x5034533352305631" or \
            timing.get("stream_order") != list(range(7)) or \
            float(timing.get(
                "maximum_total_seconds_per_step_regression_fraction", 0)) != \
            0.10:
        raise AssertionError("P4-S7 timing schedule mismatch")
    baseline = timing.get("baseline", {})
    if int(baseline.get("global_numerator", 0)) != 1 or \
            int(baseline.get("global_denominator", 0)) != 1 or \
            float(baseline.get("rho", 0)) != 0.01 or \
            int(baseline.get("bounded_max_order", 0)) != 3 or \
            int(timing.get("candidate_partial_bounded_max_order", 0)) != 2 or \
            int(timing.get(
                "candidate_callback_evaluations_per_step", 0)) != 32:
        raise AssertionError("P4-S7 timing producer mismatch")
    if not allow_smoke_policy and int(scope.get("cache_bytes", 0)) != \
            268435456:
        raise AssertionError("P4-S7 timing cache mismatch")
    return policy


def parse_log(path):
    header = None
    baseline = []
    candidates = []
    with open(path, "r") as handle:
        text = handle.read()
    for line in text.splitlines():
        if line.startswith("TIMING_PARTIAL_RUN "):
            header = parse_key_values(line)
        elif line.startswith("TIMING_REPEAT "):
            baseline.append(parse_key_values(line))
        elif line.startswith("TIMING_PARTIAL_CALLBACK_REPEAT "):
            candidates.append(parse_key_values(line))
    if header is None or len(baseline) != 7 or len(candidates) != 14:
        raise AssertionError("incomplete P4-S7 timing log {}".format(path))
    return {"path": os.path.abspath(path), "sha256": sha256_file(path),
            "header": header, "run": header, "baseline": baseline,
            "candidates": candidates, "line_count": len(text.splitlines())}


def validate_run(run, policy, allow_smoke_policy=False):
    header = run["header"]
    timing = policy["timing_census"]
    site = int(header["site_count"])
    qp = int(header["qp_total"])
    samples = int(header["sample_count"])
    cache = int(header["cache_requested_bytes"])
    if int(header.get("schema", 0)) != LOG_SCHEMA_VERSION or \
            header.get("fixture") != EXPECTED_FIXTURE or \
            [site, qp] not in policy["scope"]["physical_case_order"] or \
            samples != int(timing["sample_count"]) or \
            int(header["repeat_count"]) != int(timing["repeat_count"]) or \
            int(header["full_plan_hash"]) == 0 or \
            int(header["partial_plan_hash"]) == 0:
        raise AssertionError("P4-S7 timing run mismatch")
    if not allow_smoke_policy and cache != int(policy["scope"]["cache_bytes"]):
        raise AssertionError("P4-S7 timing run cache mismatch")
    streams = set(timing["stream_order"])
    if {int(x["rng_stream"]) for x in run["baseline"]} != streams or \
            {(int(x["candidate_index"]), int(x["rng_stream"]))
             for x in run["candidates"]} != \
            {(candidate, stream) for candidate in range(2)
             for stream in streams}:
        raise AssertionError("P4-S7 timing repeat census mismatch")
    for item in run["baseline"]:
        if int(item["site_count"]) != site or int(item["qp_total"]) != qp or \
                int(item["sample_count"]) != samples or \
                int(item["cache_requested_bytes"]) != cache or \
                int(item["fraction_index"]) != 4 or \
                int(item["rho_index"]) != 0 or \
                int(item["global_numerator"]) != 1 or \
                int(item["global_denominator"]) != 1 or \
                str(item["seed_hex"]).lower() != timing["seed_hex"].lower() or \
                int(item["accepted"]) + int(item["rejected"]) != samples or \
                int(item["global_attempted"]) != samples or \
                int(item["neighbor_attempted"]) != 0 or \
                int(item["proposal_policy_hash"]) == 0 or \
                int(item["proposal_model_hash"]) == 0:
            raise AssertionError("P4-S7 timing baseline mismatch")
        require_close(item["rho"], 0.01, "P4-S7 timing rho")
        require_close(item["total_seconds_per_step"],
                      float(item["total_step_seconds"]) / samples,
                      "P4-S7 baseline seconds per step")
        if float(item["total_step_seconds"]) <= 0.0:
            raise AssertionError("P4-S7 nonpositive baseline")
    for item in run["candidates"]:
        index = int(item["candidate_index"])
        partial_order, alpha, step_count = EXPECTED_CANDIDATES[index]
        attempted = samples * step_count
        if int(item["schema"]) != LOG_SCHEMA_VERSION or \
                int(item["site_count"]) != site or \
                int(item["qp_total"]) != qp or \
                int(item["sample_count"]) != samples or \
                int(item["cache_requested_bytes"]) != cache or \
                int(item["partial_order"]) != partial_order or \
                int(item["floor_multiplier"]) != alpha or \
                int(item["step_count"]) != step_count or \
                str(item["seed_hex"]).lower() != timing["seed_hex"].lower() or \
                int(item["accepted"]) + int(item["rejected"]) != samples or \
                int(item["final_changed"]) + \
                int(item["skipped_full_evaluations"]) != samples or \
                int(item["callback_evaluations"]) != attempted or \
                int(item["callback_terminal_amplitude_calls"]) <= 0 or \
                any(int(item[field]) == 0 for field in
                    ("proposal_policy_hash", "proposal_model_hash",
                     "callback_policy_hash", "partial_plan_hash")):
            raise AssertionError("P4-S7 timing candidate mismatch")
        component_names = (
            "partial_depth0_seconds", "partial_depth1_seconds",
            "partial_depth2_seconds", "partial_evaluation_seconds",
            "inner_overhead_seconds", "final_bounded_seconds")
        total = float(item["total_step_seconds"])
        if not math.isfinite(total) or total <= 0.0:
            raise AssertionError("P4-S7 nonpositive candidate timing")
        for name in component_names:
            value = float(item[name])
            if not math.isfinite(value) or value < 0.0 or \
                    value > total * (1.0 + 1.0e-9):
                raise AssertionError("P4-S7 timing component {}".format(name))
        require_close(item["total_seconds_per_step"], total / samples,
                      "P4-S7 candidate seconds per step")
    return {"site_count": site, "qp_total": qp,
            "rank_count": int(header["rank_count"]),
            "sample_count": samples, "cache_requested_bytes": cache,
            "sector_size": int(header["sector_size"]),
            "full_plan_hash": int(header["full_plan_hash"]),
            "partial_plan_hash": int(header["partial_plan_hash"]),
            "baseline": run["baseline"], "candidates": run["candidates"]}


def values_for(items, candidate_index=None):
    selected = items if candidate_index is None else [
        x for x in items if int(x["candidate_index"]) == candidate_index]
    values = [float(x["total_seconds_per_step"]) for x in selected]
    if len(values) != 7:
        raise AssertionError("P4-S7 timing extraction mismatch")
    return values


def analyze_command(args):
    policy = validate_policy(read_json(args.policy), args.allow_smoke_policy)
    b2 = read_json(args.stage_b2_eligibility)
    if b2.get("artifact") != "p4s7_stage_b2_eligibility" or \
            b2.get("measurement_policy_sha256") != \
            policy["source_stage_b2_policy_sha256"] or \
            b2.get("eligible_candidate_indices") != [0, 1]:
        raise AssertionError("P4-S7 timing B2 input mismatch")
    parsed = [parse_log(path) for path in args.timing_log]
    observed = [[int(x["header"]["site_count"]),
                 int(x["header"]["qp_total"])] for x in parsed]
    if sorted(observed) != sorted(policy["scope"]["physical_case_order"]) or \
            len(set(tuple(x) for x in observed)) != 6:
        raise AssertionError("P4-S7 timing requires six cases")
    validated = {(int(x["header"]["site_count"]),
                  int(x["header"]["qp_total"])):
                 validate_run(x, policy, args.allow_smoke_policy)
                 for x in parsed}
    output = os.path.abspath(args.output)
    if not os.path.isdir(output):
        os.makedirs(output)
    raw_paths = p4s3.materialize_logs(parsed, output)
    cases = []
    for site, qp in policy["scope"]["physical_case_order"]:
        item = next(x for x in parsed if int(x["header"]["site_count"]) == site
                    and int(x["header"]["qp_total"]) == qp)
        case = validated[(site, qp)]
        case["raw_log"] = raw_paths[item["path"]]
        case["raw_log_sha256"] = item["sha256"]
        case["raw_log_line_count"] = item["line_count"]
        cases.append(case)
    threshold = float(policy["timing_census"][
        "maximum_total_seconds_per_step_regression_fraction"])
    summaries = []
    eligible_indices = []
    for index, candidate in enumerate(EXPECTED_CANDIDATES):
        case_results = []
        candidate_pass = True
        for case in cases:
            baseline_values = values_for(case["baseline"])
            candidate_values = values_for(case["candidates"], index)
            baseline_median = statistics.median(baseline_values)
            candidate_median = statistics.median(candidate_values)
            if baseline_median <= 0.0:
                raise AssertionError("P4-S7 timing baseline median")
            regression = candidate_median / baseline_median - 1.0
            timing_pass = regression <= threshold
            candidate_pass = candidate_pass and timing_pass
            case_results.append({
                "site_count": case["site_count"],
                "qp_total": case["qp_total"],
                "baseline_repeat_seconds_per_step": baseline_values,
                "candidate_repeat_seconds_per_step": candidate_values,
                "baseline_median_seconds_per_step": baseline_median,
                "candidate_median_seconds_per_step": candidate_median,
                "regression_fraction": regression,
                "timing_pass": timing_pass})
        summaries.append({
            "candidate_index": index, "candidate": candidate_key(candidate),
            "timing_pass": candidate_pass,
            "maximum_regression_fraction": max(
                x["regression_fraction"] for x in case_results),
            "cases": case_results})
        if candidate_pass:
            eligible_indices.append(index)
    generated = p4c.generated_at_jst()
    policy_sha = sha256_file(args.policy)
    evidence = {"schema_version": POLICY_SCHEMA_VERSION,
                "artifact": "p4s7_partial_callback_timing",
                "generated_at_jst": generated,
                "source_commit": args.source_commit,
                "producer_binary": os.path.basename(args.producer_binary),
                "producer_binary_sha256": sha256_file(args.producer_binary),
                "timing_policy_sha256": policy_sha,
                "timing_policy": policy,
                "stage_b2_eligibility_sha256": sha256_file(
                    args.stage_b2_eligibility),
                "physical_cases": cases,
                "candidate_summaries": summaries}
    eligibility = {
        "schema_version": POLICY_SCHEMA_VERSION,
        "artifact": "p4s7_post_timing_eligibility",
        "generated_at_jst": generated, "source_commit": args.source_commit,
        "timing_policy_sha256": policy_sha,
        "candidate_order": [candidate_key(x) for x in EXPECTED_CANDIDATES],
        "eligible_candidate_indices": eligible_indices,
        "eligible_candidates": [summaries[i]["candidate"]
                                for i in eligible_indices],
        "timing_decision": "GO" if eligible_indices else "STOP",
        "next_action": "prepare_remote_package_and_request_approval"
                       if eligible_indices else
                       "preserve_stop_evidence_and_do_not_prepare_remote_job"}
    write_json(os.path.join(output, "p4s7_timing_evidence.json"), evidence)
    write_json(os.path.join(output, "p4s7_post_timing_eligibility.json"),
               eligibility)
    count = p4c.write_checksums(output)
    print("P4-S7 timing complete: {} candidates eligible; {} ledger files".
          format(len(eligible_indices), count))


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--policy", required=True)
    parser.add_argument("--stage-b2-eligibility", required=True)
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
