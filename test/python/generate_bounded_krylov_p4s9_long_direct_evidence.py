from __future__ import print_function

import argparse
import math
import os

import generate_bounded_krylov_p4c_evidence as p4c
import generate_bounded_krylov_p4s3_pilot_evidence as p4s3


POLICY_SCHEMA_VERSION = 1
DRIVER_SCHEMA_VERSION = 8
EXPECTED_FIXTURE = "p4s9_long_direct_exact_table"
EXPECTED_CANDIDATE = "direct_neighbor_rho0p01_N32768_session256m"
FAMILIES = ("S", "K", "B")
COMPONENTS = ("real", "imag")
ENTRY_COUNT = 6


read_json = p4c.read_json
write_json = p4c.write_json
sha256_file = p4c.sha256_file
fixed_sector_size = p4c.fixed_sector_size


def parse_key_values(line, skip=1):
    return {
        key: p4c.parse_number(value)
        for key, value in p4c.parse_key_values(line, skip).items()
    }


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
            "power-lanczos-zero-support-p4-r0-v9":
        raise AssertionError("P4-S9 policy identity mismatch")
    scope = policy.get("scope", {})
    expected_cases = [[4, 1], [4, 4], [6, 1], [6, 4], [8, 1], [8, 4]]
    if scope.get("site_counts") != [4, 6, 8] or \
            scope.get("qp_total") != [1, 4] or \
            scope.get("physical_case_order") != expected_cases or \
            int(scope.get("omp_threads", 0)) != 1 or \
            int(scope.get("blas_threads", 0)) != 1 or \
            int(scope.get("statistical_cache_bytes", -1)) != 268435456:
        raise AssertionError("P4-S9 scope mismatch")
    candidate = policy.get("candidate", {})
    if candidate.get("candidate_id") != EXPECTED_CANDIDATE or \
            int(candidate.get("guide_order", 0)) != 3 or \
            float(candidate.get("rho", 0.0)) != 0.01 or \
            int(candidate.get("global_numerator", -1)) != 0 or \
            int(candidate.get("global_denominator", 0)) != 1 or \
            int(candidate.get("candidate_count", 0)) != 1:
        raise AssertionError("P4-S9 candidate mismatch")
    pilot = policy.get("stage_a_exact_table_pilot", {})
    if pilot.get("driver_mode") != "long-direct" or \
            int(pilot.get("driver_schema_version", 0)) != \
            DRIVER_SCHEMA_VERSION or \
            pilot.get("fixture") != EXPECTED_FIXTURE or \
            (not allow_smoke_policy and
             int(pilot.get("sample_count_per_seed", 0)) != 131072) or \
            pilot.get("seed_hex_order") != [
                "0x5034533950494c31", "0x5034533950494c32",
                "0x5034533950494c33", "0x5034533950494c34"] or \
            int(pilot.get("predicted_official_sample_count", 0)) != 32768 or \
            float(pilot.get("maximum_predicted_se_budget_ratio", 0.0)) != \
            0.80 or float(pilot.get("maximum_tau_int", 0.0)) != 12.0:
        raise AssertionError("P4-S9 pilot gate mismatch")
    resource = policy.get("stage_c_resource", {})
    if int(resource.get("sample_count", 0)) != 32768 or \
            float(resource.get(
                "combined_runtime_seconds_per_rank_maximum", 0.0)) != 1800.0 or \
            float(resource.get(
                "proposal_seconds_per_configuration_maximum", 0.0)) != \
            0.054931640625 or \
            int(resource.get("repeat_count", 0)) != 7:
        raise AssertionError("P4-S9 resource gate mismatch")
    official = policy.get("stage_d_official_markov", {})
    if int(official.get("sample_count", 0)) != 32768 or \
            official.get("official_partition") != {
                "block_count": 16, "block_length": 2048} or \
            official.get("diagnostic_partition") != {
                "block_count": 32, "block_length": 1024} or \
            len(official.get("seed_hex_order", [])) != 4:
        raise AssertionError("P4-S9 official contract mismatch")
    return policy


def parse_pilot_log(path):
    run = None
    candidate = None
    series = []
    entries = []
    seeds = []
    decision = None
    with open(path, "r") as handle:
        text = handle.read()
    for line in text.splitlines():
        if line.startswith("PILOT_RUN "):
            if run is not None:
                raise AssertionError("duplicate P4-S9 run")
            run = parse_key_values(line)
        elif line.startswith("PILOT_CANDIDATE "):
            if candidate is not None:
                raise AssertionError("duplicate P4-S9 candidate")
            candidate = parse_key_values(line)
        elif line.startswith("PILOT_SERIES "):
            series.append(parse_key_values(line))
        elif line.startswith("PILOT_ENTRY "):
            entries.append(parse_key_values(line))
        elif line.startswith("PILOT_SEED "):
            seeds.append(parse_key_values(line))
        elif line.startswith("PILOT_DECISION "):
            if decision is not None:
                raise AssertionError("duplicate P4-S9 decision")
            decision = parse_key_values(line)
    if any(item is None for item in (run, candidate, decision)):
        raise AssertionError("incomplete P4-S9 log {}".format(path))
    return {
        "path": os.path.abspath(path), "sha256": sha256_file(path),
        "line_count": len(text.splitlines()), "run": run,
        "candidate": candidate, "series": series, "entries": entries,
        "seeds": seeds, "decision": decision,
    }


def validate_run(parsed, policy, allow_smoke_policy=False):
    run = parsed["run"]
    candidate = parsed["candidate"]
    decision = parsed["decision"]
    pilot = policy["stage_a_exact_table_pilot"]
    site = int(run["site_count"])
    qp = int(run["qp_total"])
    sample_count = int(run["sample_count_per_seed"])
    sector_size = fixed_sector_size(site)
    if int(run.get("schema", 0)) != DRIVER_SCHEMA_VERSION or \
            run.get("fixture") != EXPECTED_FIXTURE or \
            [site, qp] not in policy["scope"]["physical_case_order"] or \
            int(run["sector_size"]) != sector_size or \
            int(run["anchor_count"]) != sector_size or \
            (not allow_smoke_policy and sample_count !=
             int(pilot["sample_count_per_seed"])) or \
            (not allow_smoke_policy and
             int(run["cache_requested_bytes"]) !=
             int(policy["scope"]["statistical_cache_bytes"])):
        raise AssertionError("P4-S9 run contract mismatch")
    if int(candidate["schema"]) != DRIVER_SCHEMA_VERSION or \
            int(candidate["candidate_index"]) != 0 or \
            int(decision["candidate_index"]) != 0 or \
            int(candidate["site_count"]) != site or \
            int(candidate["qp_total"]) != qp or \
            int(candidate["sample_count_per_seed"]) != sample_count or \
            int(candidate["global_numerator"]) != 0 or \
            int(candidate["global_denominator"]) != 1 or \
            int(candidate["predicted_official_sample_count"]) != 32768:
        raise AssertionError("P4-S9 candidate schedule mismatch")
    require_close(candidate["rho"], 0.01, "P4-S9 rho")
    require_close(candidate["maximum_predicted_ratio_threshold"], 0.80,
                  "P4-S9 ratio threshold")
    require_close(candidate["maximum_tau_threshold"], 12.0,
                  "P4-S9 tau threshold")
    if float(candidate["eta"]) <= 0.0 or \
            int(candidate["anchor_count"]) != sector_size or \
            not bool(int(candidate["exact_balance_pass"])) or \
            not bool(int(candidate["restart_pass"])) or \
            any(int(candidate[field]) == 0 for field in (
                "proposal_policy_hash", "proposal_model_hash",
                "guide_policy_hash", "guide_shape_hash", "table_hash")):
        raise AssertionError("P4-S9 identity/balance mismatch")
    expected_series = {
        (seed, family, entry, component)
        for seed in range(4) for family in FAMILIES
        for entry in range(ENTRY_COUNT) for component in COMPONENTS
    }
    observed_series = {
        (int(x["seed_index"]), x["family"], int(x["entry"]), x["component"])
        for x in parsed["series"]
    }
    expected_entries = {
        (seed, family, entry) for seed in range(4)
        for family in FAMILIES for entry in range(ENTRY_COUNT)
    }
    observed_entries = {
        (int(x["seed_index"]), x["family"], int(x["entry"]))
        for x in parsed["entries"]
    }
    if observed_series != expected_series or observed_entries != expected_entries or \
            len(parsed["seeds"]) != 4:
        raise AssertionError("P4-S9 series census mismatch")
    for seed_index, seed in enumerate(parsed["seeds"]):
        if int(seed["seed_index"]) != seed_index or \
                str(seed["seed_hex"]).lower() != \
                pilot["seed_hex_order"][seed_index].lower() or \
                int(seed["accepted"]) + int(seed["rejected"]) != \
                sample_count or \
                int(seed["neighbor_attempted"]) != sample_count or \
                int(seed["global_attempted"]) != 0 or \
                int(seed["global_self"]) != 0 or \
                int(seed["shell_attempted"]) != 0 or \
                int(seed["trace_hash"]) == 0:
            raise AssertionError("P4-S9 seed census mismatch")
    maximum_tau = max(float(x["tau_int"]) for x in parsed["series"])
    maximum_ratio = max(float(x["predicted_se_budget_ratio"])
                        for x in parsed["entries"])
    require_close(decision["maximum_tau_int"], maximum_tau,
                  "P4-S9 maximum tau")
    require_close(decision["maximum_predicted_se_budget_ratio"],
                  maximum_ratio, "P4-S9 maximum ratio")
    if site == 4:
        tolerance = float(pilot[
            "exact_l4_balance_tolerance_multiplier_dbl_epsilon"
        ]) * 2.220446049250313e-16
        for field in ("proposal_row_residual", "proposal_symmetry_residual",
                      "db_residual", "stationary_residual"):
            if float(candidate[field]) > tolerance:
                raise AssertionError("P4-S9 L4 {}".format(field))
    eligible = maximum_ratio <= float(
        pilot["maximum_predicted_se_budget_ratio"]) and \
        maximum_tau <= float(pilot["maximum_tau_int"]) and \
        bool(int(decision["table_pass"])) and \
        bool(int(decision["balance_pass"])) and \
        bool(int(decision["restart_pass"]))
    if bool(int(decision["physical_case_eligible"])) != eligible:
        raise AssertionError("P4-S9 eligibility mismatch")
    return {
        "site_count": site, "qp_total": qp,
        "rank_count": int(run["rank_count"]),
        "sample_count_per_seed": sample_count,
        "cache_requested_bytes": int(run["cache_requested_bytes"]),
        "sector_size": sector_size, "plan_hash": int(run["plan_hash"]),
        "candidate_id": EXPECTED_CANDIDATE,
        "physical_case_eligible": eligible,
        "maximum_tau_int": maximum_tau,
        "maximum_predicted_se_budget_ratio": maximum_ratio,
        "maximum_iid_se_budget_ratio": float(
            candidate["maximum_iid_se_budget_ratio"]),
        "required_tau_for_0p90": float(candidate["required_tau_for_0p90"]),
        "accepted_by_seed": [int(x["accepted"]) for x in parsed["seeds"]],
        "rejected_by_seed": [int(x["rejected"]) for x in parsed["seeds"]],
        "trace_hash_by_seed": [int(x["trace_hash"]) for x in parsed["seeds"]],
        "proposal_policy_hash": int(candidate["proposal_policy_hash"]),
        "proposal_model_hash": int(candidate["proposal_model_hash"]),
        "guide_policy_hash": int(candidate["guide_policy_hash"]),
        "guide_shape_hash": int(candidate["guide_shape_hash"]),
        "table_hash": int(candidate["table_hash"]),
    }


def analyze_command(args):
    policy = validate_policy(read_json(args.policy), args.allow_smoke_policy)
    parsed = [parse_pilot_log(path) for path in args.pilot_log]
    observed = [[int(x["run"]["site_count"]), int(x["run"]["qp_total"])]
                for x in parsed]
    expected = policy["scope"]["physical_case_order"]
    if sorted(observed) != sorted(expected) or \
            len(set(tuple(x) for x in observed)) != len(expected):
        raise AssertionError("P4-S9 requires six unique physical cases")
    validated = {
        (int(x["run"]["site_count"]), int(x["run"]["qp_total"])):
        validate_run(x, policy, args.allow_smoke_policy) for x in parsed
    }
    output = os.path.abspath(args.output)
    if not os.path.isdir(output):
        os.makedirs(output)
    raw_paths = p4s3.materialize_logs(parsed, output)
    cases = []
    for site, qp in expected:
        item = next(x for x in parsed
                    if int(x["run"]["site_count"]) == site and
                    int(x["run"]["qp_total"]) == qp)
        summary = validated[(site, qp)]
        summary["raw_log"] = raw_paths[item["path"]]
        summary["raw_log_sha256"] = item["sha256"]
        summary["raw_log_line_count"] = item["line_count"]
        cases.append(summary)
    eligible = all(x["physical_case_eligible"] for x in cases)
    generated_at = p4c.generated_at_jst()
    policy_sha = sha256_file(args.policy)
    evidence = {
        "schema_version": POLICY_SCHEMA_VERSION,
        "artifact": "p4s9_long_direct_stage_a_evidence",
        "generated_at_jst": generated_at,
        "source_commit": args.source_commit,
        "producer_binary": os.path.basename(args.producer_binary),
        "producer_binary_sha256": sha256_file(args.producer_binary),
        "measurement_policy_sha256": policy_sha,
        "measurement_policy": policy,
        "candidate_id": EXPECTED_CANDIDATE,
        "stage_a_eligible": eligible,
        "maximum_tau_int": max(x["maximum_tau_int"] for x in cases),
        "maximum_predicted_se_budget_ratio": max(
            x["maximum_predicted_se_budget_ratio"] for x in cases),
        "physical_cases": cases,
    }
    eligibility = {
        "schema_version": POLICY_SCHEMA_VERSION,
        "artifact": "p4s9_frozen_stage_a_eligibility",
        "generated_at_jst": generated_at,
        "source_commit": args.source_commit,
        "measurement_policy_sha256": policy_sha,
        "candidate_id": EXPECTED_CANDIDATE,
        "stage_a_decision": "ELIGIBLE" if eligible else "STOP",
        "stage_b_authorized": eligible,
        "failure_action": "implement_frozen_session_contract"
        if eligible else "do_not_implement_session_or_run_timing",
    }
    write_json(os.path.join(output, "p4s9_stage_a_evidence.json"), evidence)
    write_json(os.path.join(output, "p4s9_stage_a_eligibility.json"),
               eligibility)
    count = p4c.write_checksums(output)
    print("P4-S9 Stage A complete: {}; {} ledger files".format(
        eligibility["stage_a_decision"], count))


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--pilot-log", nargs="+", required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--producer-binary", required=True)
    parser.add_argument("--allow-smoke-policy", action="store_true")
    return parser


def main(argv=None):
    analyze_command(build_parser().parse_args(argv))


if __name__ == "__main__":
    main()
