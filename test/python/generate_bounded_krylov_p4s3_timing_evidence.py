from __future__ import print_function

import argparse
import math
import os
import shutil
import statistics

import generate_bounded_krylov_p4c_evidence as p4c
import generate_bounded_krylov_p4s3_pilot_evidence as pilot


SCHEMA_VERSION = 1
EXPECTED_FIXTURE = "p4s3_matched_timing"
FRACTIONS = ((0, 1), (1, 8), (1, 4), (1, 2), (1, 1))
RHO_VALUES = (1.0e-2, 1.0e-4, 1.0e-6)


def parse_key_values(line):
    return {
        key: p4c.parse_number(value)
        for key, value in p4c.parse_key_values(line, 1).items()
    }


def parse_timing_log(path):
    header = None
    repeats = []
    with open(path, "r") as handle:
        text = handle.read()
    for line in text.splitlines():
        if line.startswith("TIMING_RUN "):
            if header is not None:
                raise AssertionError("duplicate TIMING_RUN")
            header = parse_key_values(line)
        elif line.startswith("TIMING_REPEAT "):
            repeats.append(parse_key_values(line))
    if header is None or len(repeats) != 5 * 3 * 7:
        raise AssertionError("incomplete P4-S3 timing log {}".format(path))
    return {
        "path": os.path.abspath(path),
        "sha256": p4c.sha256_file(path),
        "header": header,
        "repeats": repeats,
        "line_count": len(text.splitlines()),
    }


def require_close(actual, expected, label, tolerance=1.0e-12):
    actual = float(actual)
    expected = float(expected)
    if not math.isfinite(actual) or not math.isfinite(expected) or \
            abs(actual - expected) > tolerance * max(1.0, abs(expected)):
        raise AssertionError("{} mismatch: {} != {}".format(
            label, actual, expected))


def validate_timing_run(run, policy, allow_smoke_policy=False):
    header = run["header"]
    timing_policy = policy["timing_census"]
    if int(header.get("schema", 0)) != SCHEMA_VERSION or \
            header.get("fixture") != EXPECTED_FIXTURE:
        raise AssertionError("timing run header mismatch")
    site = int(header["site_count"])
    qp = int(header["qp_total"])
    if [site, qp] not in policy["scope"]["physical_case_order"]:
        raise AssertionError("timing physical case outside scope")
    sample_count = int(header["sample_count"])
    cache = int(header["cache_requested_bytes"])
    if not allow_smoke_policy and sample_count != \
            int(timing_policy["sample_count"]):
        raise AssertionError("timing sample count mismatch")
    if int(header["repeat_count"]) != int(timing_policy["repeat_count"]):
        raise AssertionError("timing repeat count mismatch")
    if not allow_smoke_policy and cache != \
            int(policy["scope"]["statistical_cache_bytes"]):
        raise AssertionError("timing cache mismatch")
    expected = {
        (fraction_index, rho_index, stream)
        for fraction_index in range(5)
        for rho_index in range(3)
        for stream in range(7)
    }
    observed = {
        (int(item["fraction_index"]), int(item["rho_index"]),
         int(item["rng_stream"]))
        for item in run["repeats"]
    }
    if observed != expected:
        raise AssertionError("timing schedule mismatch")
    for item in run["repeats"]:
        fraction_index = int(item["fraction_index"])
        rho_index = int(item["rho_index"])
        numerator, denominator = FRACTIONS[fraction_index]
        if int(item["site_count"]) != site or int(item["qp_total"]) != qp:
            raise AssertionError("timing repeat physical-case mismatch")
        if int(item["sample_count"]) != sample_count or \
                int(item["cache_requested_bytes"]) != cache:
            raise AssertionError("timing repeat run-contract mismatch")
        if int(item["global_numerator"]) != numerator or \
                int(item["global_denominator"]) != denominator:
            raise AssertionError("timing proposal fraction mismatch")
        require_close(item["rho"], RHO_VALUES[rho_index], "timing rho")
        if str(item["seed_hex"]).lower() != \
                policy["official_markov"]["official_seed_hex"].lower():
            raise AssertionError("timing seed mismatch")
        if int(item["rng_stream"]) not in \
                timing_policy["stream_order"]:
            raise AssertionError("timing stream mismatch")
        if int(item["accepted"]) + int(item["rejected"]) != sample_count:
            raise AssertionError("timing accept/reject census mismatch")
        if int(item["neighbor_attempted"]) + \
                int(item["global_attempted"]) != sample_count:
            raise AssertionError("timing component census mismatch")
        if int(item["global_self"]) > int(item["global_attempted"]):
            raise AssertionError("timing global-self census mismatch")
        if fraction_index == 0 and int(item["global_attempted"]) != 0:
            raise AssertionError("p=0 timing control used global proposal")
        if fraction_index == 4 and \
                int(item["global_attempted"]) != sample_count:
            raise AssertionError("p=1 timing candidate used neighbor proposal")
        total = float(item["total_step_seconds"])
        per_step = float(item["total_seconds_per_step"])
        require_close(per_step, total / sample_count,
                      "timing seconds per step")
        if total <= 0.0:
            raise AssertionError("nonpositive total timing")
        for name in ("proposal_seconds", "component_selection_seconds",
                     "global_subset_seconds",
                     "bounded_evaluation_seconds"):
            value = float(item[name])
            if not math.isfinite(value) or value < 0.0 or \
                    value > total * (1.0 + 1.0e-9):
                raise AssertionError("invalid timing component {}".format(
                    name))
        if int(item["proposal_policy_hash"]) == 0 or \
                int(item["proposal_model_hash"]) == 0:
            raise AssertionError("timing policy/model hash missing")
    return {
        "site_count": site,
        "qp_total": qp,
        "rank_count": int(header["rank_count"]),
        "sample_count": sample_count,
        "cache_requested_bytes": cache,
        "sector_size": int(header["sector_size"]),
        "plan_hash": int(header["plan_hash"]),
        "repeats": run["repeats"],
    }


def repeat_values(case, fraction_index, rho_index):
    values = [
        float(item["total_seconds_per_step"])
        for item in case["repeats"]
        if int(item["fraction_index"]) == fraction_index and
        int(item["rho_index"]) == rho_index
    ]
    if len(values) != 7:
        raise AssertionError("timing repeat extraction mismatch")
    return values


def materialize_logs(runs, output):
    raw_root = os.path.join(output, "raw_timing")
    if not os.path.isdir(raw_root):
        os.makedirs(raw_root)
    mapped = {}
    for run in runs:
        name = "L{}-QP{}-{}-{}.log".format(
            int(run["header"]["site_count"]),
            int(run["header"]["qp_total"]), run["sha256"][:12],
            os.path.basename(run["path"]))
        destination = os.path.join(raw_root, name)
        shutil.copyfile(run["path"], destination)
        if p4c.sha256_file(destination) != run["sha256"]:
            raise AssertionError("timing raw-log checksum mismatch")
        mapped[run["path"]] = os.path.relpath(destination, output)
    return mapped


def analyze_command(args):
    policy = pilot.validate_policy(
        pilot.read_json(args.policy), args.allow_smoke_policy)
    pilot_eligibility = pilot.read_json(args.pilot_eligibility)
    if pilot_eligibility.get("artifact") != \
            "p4s3_frozen_pilot_eligibility":
        raise AssertionError("timing input is not frozen pilot eligibility")
    if pilot_eligibility.get("measurement_policy_sha256") != \
            p4c.sha256_file(args.policy):
        raise AssertionError("timing/pilot policy checksum mismatch")
    parsed = [parse_timing_log(path) for path in args.timing_log]
    observed_cases = [
        [int(run["header"]["site_count"]),
         int(run["header"]["qp_total"])]
        for run in parsed
    ]
    if sorted(observed_cases) != sorted(policy["scope"]["physical_case_order"]):
        raise AssertionError("timing evidence must contain six unique cases")
    if len(set(tuple(item) for item in observed_cases)) != len(observed_cases):
        raise AssertionError("duplicate timing physical case")
    validated = {
        (int(run["header"]["site_count"]),
         int(run["header"]["qp_total"])):
        validate_timing_run(run, policy, args.allow_smoke_policy)
        for run in parsed
    }
    output = os.path.abspath(args.output)
    if not os.path.isdir(output):
        os.makedirs(output)
    raw_paths = materialize_logs(parsed, output)
    ordered_cases = []
    for site, qp in policy["scope"]["physical_case_order"]:
        run = next(item for item in parsed
                   if int(item["header"]["site_count"]) == site and
                   int(item["header"]["qp_total"]) == qp)
        case = validated[(site, qp)]
        case["raw_log"] = raw_paths[run["path"]]
        case["raw_log_sha256"] = run["sha256"]
        ordered_cases.append(case)
    threshold = float(policy["timing_census"][
        "maximum_total_seconds_per_step_regression_fraction"])
    candidate_results = []
    timing_eligible = []
    for candidate_index in pilot_eligibility["eligible_candidate_indices"]:
        candidate_index = int(candidate_index)
        fraction_index = candidate_index // 3 + 1
        rho_index = candidate_index % 3
        case_results = []
        passes = True
        for case in ordered_cases:
            baseline_values = repeat_values(case, 0, rho_index)
            candidate_values = repeat_values(
                case, fraction_index, rho_index)
            baseline_median = statistics.median(baseline_values)
            candidate_median = statistics.median(candidate_values)
            if baseline_median <= 0.0:
                raise AssertionError("nonpositive timing baseline median")
            regression = candidate_median / baseline_median - 1.0
            case_pass = regression <= threshold
            passes = passes and case_pass
            case_results.append({
                "site_count": case["site_count"],
                "qp_total": case["qp_total"],
                "baseline_repeat_seconds_per_step": baseline_values,
                "candidate_repeat_seconds_per_step": candidate_values,
                "baseline_median_seconds_per_step": baseline_median,
                "candidate_median_seconds_per_step": candidate_median,
                "regression_fraction": regression,
                "timing_pass": case_pass,
            })
        result = {
            "candidate_index": candidate_index,
            "candidate": pilot_eligibility["candidate_order"][candidate_index],
            "timing_pass": passes,
            "maximum_regression_fraction": max(
                item["regression_fraction"] for item in case_results),
            "cases": case_results,
        }
        candidate_results.append(result)
        if passes:
            timing_eligible.append(candidate_index)
    evidence = {
        "schema_version": SCHEMA_VERSION,
        "artifact": "p4s3_matched_timing_evidence",
        "generated_at_jst": p4c.generated_at_jst(),
        "source_commit": args.source_commit,
        "measurement_policy_sha256": p4c.sha256_file(args.policy),
        "pilot_eligibility_sha256": p4c.sha256_file(args.pilot_eligibility),
        "physical_cases": ordered_cases,
        "candidate_results": candidate_results,
    }
    eligibility = {
        "schema_version": SCHEMA_VERSION,
        "artifact": "p4s3_frozen_official_eligibility",
        "generated_at_jst": evidence["generated_at_jst"],
        "source_commit": args.source_commit,
        "measurement_policy_sha256": evidence["measurement_policy_sha256"],
        "pilot_eligibility_sha256": evidence["pilot_eligibility_sha256"],
        "pilot_eligible_candidate_indices":
            pilot_eligibility["eligible_candidate_indices"],
        "timing_eligible_candidate_indices": timing_eligible,
        "eligible_candidates": [
            pilot_eligibility["candidate_order"][index]
            for index in timing_eligible
        ],
        "eligibility_decision": "ELIGIBLE" if timing_eligible else "STOP",
        "failure_action": (
            "freeze_source_package_and_request_remote_submission_approval"
            if timing_eligible else "do_not_submit_remote_job"
        ),
    }
    pilot.write_json(os.path.join(output, "p4s3_timing_evidence.json"),
                     evidence)
    pilot.write_json(os.path.join(
        output, "p4s3_official_eligibility.json"), eligibility)
    count = p4c.write_checksums(output)
    print("P4-S3 timing evidence complete: {} candidates remain; "
          "{} ledger files".format(len(timing_eligible), count))


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--policy", required=True)
    parser.add_argument("--pilot-eligibility", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--timing-log", nargs="+", required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--allow-smoke-policy", action="store_true")
    return parser


def main(argv=None):
    analyze_command(build_parser().parse_args(argv))


if __name__ == "__main__":
    main()
