from __future__ import print_function

import argparse
import math
import os

import generate_bounded_krylov_p4c_evidence as p4c
import generate_bounded_krylov_p4s3_pilot_evidence as p4s3


SCHEMA_VERSION = 2
EXPECTED_FIXTURE = "p4s4_fixed_distance_exact_table"
FAMILIES = ("S", "K", "B")
ENTRY_COUNT = 6
COMPONENTS = ("real", "imag")
CANDIDATE_COUNT = 3


def read_json(path):
    return p4c.read_json(path)


def write_json(path, value):
    p4c.write_json(path, value)


def sha256_file(path):
    return p4c.sha256_file(path)


def parse_number(value):
    return p4c.parse_number(value)


def parse_key_values(line, skip=1):
    return {
        key: parse_number(value)
        for key, value in p4c.parse_key_values(line, skip).items()
    }


def require_close(actual, expected, label, tolerance=1.0e-12):
    actual = float(actual)
    expected = float(expected)
    if not math.isfinite(actual) or not math.isfinite(expected) or \
            abs(actual - expected) > tolerance * max(1.0, abs(expected)):
        raise AssertionError("{} mismatch: {} != {}".format(
            label, actual, expected))


def fixed_sector_size(site_count):
    return p4c.fixed_sector_size(site_count)


def candidate_key(candidate):
    return "d{}/{}_rho_{:.0e}".format(
        int(candidate["distance_numerator"]),
        int(candidate["distance_denominator"]),
        float(candidate["rho"]),
    )


def resolved_distance(site_count, numerator, denominator):
    product = int(site_count) * int(numerator)
    value = product // int(denominator)
    remainder = product % int(denominator)
    if remainder >= (int(denominator) + 1) // 2:
        value += 1
    return max(1, value)


def validate_policy(policy, allow_smoke_policy=False):
    if int(policy.get("schema_version", 0)) != SCHEMA_VERSION:
        raise AssertionError("P4-S4 policy schema mismatch")
    if policy.get("policy_id") != "power-lanczos-zero-support-p4-r0-v4":
        raise AssertionError("P4-S4 policy id mismatch")
    scope = policy.get("scope", {})
    expected_cases = [[4, 1], [4, 4], [6, 1], [6, 4], [8, 1], [8, 4]]
    if scope.get("site_counts") != [4, 6, 8] or \
            scope.get("qp_total") != [1, 4] or \
            scope.get("physical_case_order") != expected_cases or \
            float(scope.get("rho", 0.0)) != 1.0e-2:
        raise AssertionError("P4-S4 physical scope mismatch")
    proposal = policy.get("proposal", {})
    if proposal.get("kernel") != "neighbor_fixed_distance_shell_mixture" or \
            int(proposal.get("policy_version", 0)) != 2 or \
            int(proposal.get("neighbor_numerator", 0)) != 1 or \
            int(proposal.get("neighbor_denominator", 0)) != 8 or \
            int(proposal.get("shell_numerator", 0)) != 7 or \
            int(proposal.get("shell_denominator", 0)) != 8 or \
            proposal.get("distance_rounding") != \
            "round_half_up_v1_minimum_one_no_upper_clamp":
        raise AssertionError("P4-S4 proposal contract mismatch")
    candidates = proposal.get("candidate_order", [])
    observed = [
        (int(item["distance_numerator"]),
         int(item["distance_denominator"]), float(item["rho"]))
        for item in candidates
    ]
    if observed != [(1, 3, 1.0e-2), (1, 2, 1.0e-2),
                    (2, 3, 1.0e-2)]:
        raise AssertionError("P4-S4 candidate order mismatch")
    pilot = policy.get("exact_table_pilot", {})
    if not allow_smoke_policy and \
            int(pilot.get("sample_count_per_seed", 0)) != 131072:
        raise AssertionError("P4-S4 pilot sample count mismatch")
    if pilot.get("seed_hex_order") != [
            "0x5034533450494c31", "0x5034533450494c32",
            "0x5034533450494c33", "0x5034533450494c34"]:
        raise AssertionError("P4-S4 pilot seed order mismatch")
    if float(pilot.get("maximum_predicted_se_budget_ratio", 0.0)) != 0.90 or \
            float(pilot.get("maximum_tau_int", 0.0)) != 12.0 or \
            int(pilot.get("predicted_official_sample_count", 0)) != 4096:
        raise AssertionError("P4-S4 pilot gate mismatch")
    if policy.get("official_markov", {}).get("official_seed_hex") != \
            "0x5034533452305631":
        raise AssertionError("P4-S4 official seed mismatch")
    return policy


def parse_pilot_log(path):
    run = None
    candidates = []
    series = []
    entries = []
    seeds = []
    decisions = []
    with open(path, "r") as handle:
        text = handle.read()
    for line in text.splitlines():
        if line.startswith("PILOT_RUN "):
            if run is not None:
                raise AssertionError("duplicate PILOT_RUN")
            run = parse_key_values(line)
        elif line.startswith("PILOT_CANDIDATE "):
            candidates.append(parse_key_values(line))
        elif line.startswith("PILOT_SERIES "):
            series.append(parse_key_values(line))
        elif line.startswith("PILOT_ENTRY "):
            entries.append(parse_key_values(line))
        elif line.startswith("PILOT_SEED "):
            seeds.append(parse_key_values(line))
        elif line.startswith("PILOT_DECISION "):
            decisions.append(parse_key_values(line))
    if run is None or len(candidates) != CANDIDATE_COUNT or \
            len(decisions) != CANDIDATE_COUNT:
        raise AssertionError("incomplete P4-S4 pilot log {}".format(path))
    return {
        "path": os.path.abspath(path),
        "sha256": sha256_file(path),
        "run": run,
        "candidates": candidates,
        "series": series,
        "entries": entries,
        "seeds": seeds,
        "decisions": decisions,
        "line_count": len(text.splitlines()),
    }


def validate_candidate(run, candidate_index, policy,
                       allow_smoke_policy=False):
    policy_candidate = policy["proposal"]["candidate_order"][candidate_index]
    proposal = policy["proposal"]
    pilot = policy["exact_table_pilot"]
    header = run["candidates"][candidate_index]
    decision = run["decisions"][candidate_index]
    site = int(run["run"]["site_count"])
    qp = int(run["run"]["qp_total"])
    sample_count = int(run["run"]["sample_count_per_seed"])
    if int(header["candidate_index"]) != candidate_index or \
            int(decision["candidate_index"]) != candidate_index or \
            int(header["site_count"]) != site or \
            int(header["qp_total"]) != qp or \
            int(header["sample_count_per_seed"]) != sample_count:
        raise AssertionError("P4-S4 candidate schedule mismatch")
    if int(header["global_numerator"]) != 0 or \
            int(header["global_denominator"]) != 0 or \
            int(header["neighbor_numerator"]) != \
            int(proposal["neighbor_numerator"]) or \
            int(header["neighbor_denominator"]) != \
            int(proposal["neighbor_denominator"]) or \
            int(header["distance_numerator"]) != \
            int(policy_candidate["distance_numerator"]) or \
            int(header["distance_denominator"]) != \
            int(policy_candidate["distance_denominator"]) or \
            int(header["distance_rounding_rule"]) != 1:
        raise AssertionError("P4-S4 proposal schedule mismatch")
    require_close(header["rho"], policy_candidate["rho"], "P4-S4 rho")
    expected_distance = resolved_distance(
        site, policy_candidate["distance_numerator"],
        policy_candidate["distance_denominator"])
    if int(header["shell_max_distance"]) != site or \
            int(header["shell_distance"]) != expected_distance or \
            int(header["shell_count"]) <= 0:
        raise AssertionError("P4-S4 resolved distance mismatch")
    if int(header["anchor_count"]) != fixed_sector_size(site) or \
            int(header["proposal_policy_hash"]) == 0 or \
            int(header["proposal_model_hash"]) == 0 or \
            int(header["table_hash"]) == 0:
        raise AssertionError("P4-S4 anchor/hash contract missing")
    if not allow_smoke_policy and sample_count != \
            int(pilot["sample_count_per_seed"]):
        raise AssertionError("P4-S4 official pilot sample count mismatch")

    candidate_series = [item for item in run["series"]
                        if int(item["candidate_index"]) == candidate_index]
    candidate_entries = [item for item in run["entries"]
                         if int(item["candidate_index"]) == candidate_index]
    candidate_seeds = [item for item in run["seeds"]
                       if int(item["candidate_index"]) == candidate_index]
    if len(candidate_series) != 4 * len(FAMILIES) * ENTRY_COUNT * 2 or \
            len(candidate_entries) != 4 * len(FAMILIES) * ENTRY_COUNT or \
            len(candidate_seeds) != 4:
        raise AssertionError("P4-S4 candidate census mismatch")
    expected_series = {
        (seed_index, family, entry, component)
        for seed_index in range(4)
        for family in FAMILIES
        for entry in range(ENTRY_COUNT)
        for component in COMPONENTS
    }
    observed_series = {
        (int(item["seed_index"]), item["family"], int(item["entry"]),
         item["component"])
        for item in candidate_series
    }
    if observed_series != expected_series:
        raise AssertionError("P4-S4 series identity mismatch")
    for seed_index, seed in enumerate(candidate_seeds):
        if int(seed["seed_index"]) != seed_index or \
                str(seed["seed_hex"]).lower() != \
                pilot["seed_hex_order"][seed_index].lower():
            raise AssertionError("P4-S4 seed order mismatch")
        if int(seed["accepted"]) + int(seed["rejected"]) != sample_count or \
                int(seed["neighbor_attempted"]) + \
                int(seed["shell_attempted"]) != sample_count or \
                int(seed["global_attempted"]) != 0 or \
                int(seed["global_self"]) != 0 or \
                int(seed["trace_hash"]) == 0:
            raise AssertionError("P4-S4 chain census mismatch")
    maximum_tau = max(float(item["tau_int"]) for item in candidate_series)
    maximum_ratio = max(float(item["predicted_se_budget_ratio"])
                        for item in candidate_entries)
    require_close(decision["maximum_tau_int"], maximum_tau,
                  "P4-S4 maximum tau")
    require_close(decision["maximum_predicted_se_budget_ratio"],
                  maximum_ratio, "P4-S4 maximum predicted ratio")
    balance_pass = bool(int(header["exact_balance_pass"])) and \
        bool(int(decision["balance_pass"]))
    restart_pass = bool(int(header["restart_pass"])) and \
        bool(int(decision["restart_pass"]))
    if site == 4:
        tolerance = float(
            pilot["exact_l4_balance_tolerance_multiplier_dbl_epsilon"]
        ) * 2.220446049250313e-16
        for field in ("db_residual", "stationary_residual",
                      "proposal_row_residual",
                      "proposal_symmetry_residual"):
            if float(header[field]) > tolerance:
                raise AssertionError("P4-S4 exact L4 {}".format(field))
        if not bool(int(header["shell_cardinality_pass"])):
            raise AssertionError("P4-S4 shell cardinality oracle failed")
    expected_eligible = (
        maximum_ratio <= float(pilot["maximum_predicted_se_budget_ratio"])
        and maximum_tau <= float(pilot["maximum_tau_int"])
        and bool(int(decision["table_pass"]))
        and balance_pass and restart_pass
    )
    if bool(int(decision["physical_case_eligible"])) != expected_eligible:
        raise AssertionError("P4-S4 physical eligibility mismatch")
    return {
        "candidate_index": candidate_index,
        "candidate": candidate_key(policy_candidate),
        "distance_numerator": int(policy_candidate["distance_numerator"]),
        "distance_denominator": int(policy_candidate["distance_denominator"]),
        "resolved_maximum_distance": int(header["shell_max_distance"]),
        "resolved_distance": int(header["shell_distance"]),
        "shell_count": int(header["shell_count"]),
        "rho": float(policy_candidate["rho"]),
        "physical_case_eligible": expected_eligible,
        "maximum_tau_int": maximum_tau,
        "maximum_predicted_se_budget_ratio": maximum_ratio,
        "proposal_policy_hash": int(header["proposal_policy_hash"]),
        "proposal_model_hash": int(header["proposal_model_hash"]),
        "table_hash": int(header["table_hash"]),
        "anchor_count": int(header["anchor_count"]),
        "exact_balance_pass": balance_pass,
        "shell_cardinality_pass": bool(int(header["shell_cardinality_pass"])),
        "restart_pass": restart_pass,
        "seed_summaries": candidate_seeds,
    }


def validate_run(run, policy, allow_smoke_policy=False):
    header = run["run"]
    if int(header.get("schema", 0)) != SCHEMA_VERSION or \
            header.get("fixture") != EXPECTED_FIXTURE:
        raise AssertionError("P4-S4 run header mismatch")
    site = int(header["site_count"])
    qp = int(header["qp_total"])
    if [site, qp] not in policy["scope"]["physical_case_order"] or \
            int(header["sector_size"]) != fixed_sector_size(site) or \
            int(header["anchor_count"]) != fixed_sector_size(site):
        raise AssertionError("P4-S4 physical case/anchor mismatch")
    if not allow_smoke_policy and \
            int(header["cache_requested_bytes"]) != \
            int(policy["scope"]["statistical_cache_bytes"]):
        raise AssertionError("P4-S4 cache policy mismatch")
    candidates = [
        validate_candidate(run, index, policy, allow_smoke_policy)
        for index in range(CANDIDATE_COUNT)
    ]
    if len(set(item["resolved_distance"] for item in candidates)) != \
            CANDIDATE_COUNT:
        raise AssertionError("P4-S4 electronic candidates collapsed")
    if len(set(item["table_hash"] for item in candidates)) != 1:
        raise AssertionError("P4-S4 exact table changed by proposal distance")
    return {
        "site_count": site,
        "qp_total": qp,
        "rank_count": int(header["rank_count"]),
        "sample_count_per_seed": int(header["sample_count_per_seed"]),
        "cache_requested_bytes": int(header["cache_requested_bytes"]),
        "sector_size": int(header["sector_size"]),
        "plan_hash": int(header["plan_hash"]),
        "candidates": candidates,
    }


def analyze_command(args):
    policy = validate_policy(read_json(args.policy), args.allow_smoke_policy)
    parsed = [parse_pilot_log(path) for path in args.pilot_log]
    observed_cases = [
        [int(run["run"]["site_count"]), int(run["run"]["qp_total"])]
        for run in parsed
    ]
    if sorted(observed_cases) != \
            sorted(policy["scope"]["physical_case_order"]) or \
            len(set(tuple(case) for case in observed_cases)) != len(parsed):
        raise AssertionError("P4-S4 evidence requires six unique cases")
    validated = {
        (int(run["run"]["site_count"]), int(run["run"]["qp_total"])):
        validate_run(run, policy, args.allow_smoke_policy)
        for run in parsed
    }
    output = os.path.abspath(args.output)
    if not os.path.isdir(output):
        os.makedirs(output)
    raw_paths = p4s3.materialize_logs(parsed, output)
    ordered_cases = []
    for site, qp in policy["scope"]["physical_case_order"]:
        run = next(item for item in parsed
                   if int(item["run"]["site_count"]) == site and
                   int(item["run"]["qp_total"]) == qp)
        summary = validated[(site, qp)]
        summary["raw_log"] = raw_paths[run["path"]]
        summary["raw_log_sha256"] = run["sha256"]
        summary["raw_log_line_count"] = run["line_count"]
        ordered_cases.append(summary)

    candidate_summaries = []
    eligible_indices = []
    for candidate_index, policy_candidate in enumerate(
            policy["proposal"]["candidate_order"]):
        cases = [case["candidates"][candidate_index]
                 for case in ordered_cases]
        eligible = all(case["physical_case_eligible"] for case in cases)
        summary = {
            "candidate_index": candidate_index,
            "candidate": candidate_key(policy_candidate),
            "distance_numerator": int(policy_candidate["distance_numerator"]),
            "distance_denominator": int(policy_candidate["distance_denominator"]),
            "rho": float(policy_candidate["rho"]),
            "pilot_eligible": eligible,
            "maximum_tau_int": max(case["maximum_tau_int"]
                                   for case in cases),
            "maximum_predicted_se_budget_ratio": max(
                case["maximum_predicted_se_budget_ratio"]
                for case in cases),
            "cases": [{
                "site_count": ordered_cases[index]["site_count"],
                "qp_total": ordered_cases[index]["qp_total"],
                "resolved_distance": case["resolved_distance"],
                "shell_count": case["shell_count"],
                "physical_case_eligible": case["physical_case_eligible"],
                "maximum_tau_int": case["maximum_tau_int"],
                "maximum_predicted_se_budget_ratio":
                    case["maximum_predicted_se_budget_ratio"],
                "table_hash": case["table_hash"],
                "anchor_count": case["anchor_count"],
            } for index, case in enumerate(cases)],
        }
        candidate_summaries.append(summary)
        if eligible:
            eligible_indices.append(candidate_index)
    generated_at = p4c.generated_at_jst()
    evidence = {
        "schema_version": SCHEMA_VERSION,
        "artifact": "p4s4_fixed_distance_exact_table_pilot",
        "generated_at_jst": generated_at,
        "source_commit": args.source_commit,
        "producer_binary": os.path.basename(args.producer_binary),
        "producer_binary_sha256": sha256_file(args.producer_binary),
        "measurement_policy_sha256": sha256_file(args.policy),
        "measurement_policy": policy,
        "physical_cases": ordered_cases,
        "candidate_summaries": candidate_summaries,
    }
    eligibility = {
        "schema_version": SCHEMA_VERSION,
        "artifact": "p4s4_frozen_pilot_eligibility",
        "generated_at_jst": generated_at,
        "source_commit": args.source_commit,
        "measurement_policy_sha256": evidence["measurement_policy_sha256"],
        "candidate_order": [candidate_key(item) for item in
                            policy["proposal"]["candidate_order"]],
        "eligible_candidate_indices": eligible_indices,
        "eligible_candidates": [candidate_summaries[index]["candidate"]
                                for index in eligible_indices],
        "preferred_candidate_index": eligible_indices[0]
        if eligible_indices else None,
        "pilot_decision": "ELIGIBLE" if eligible_indices else "STOP",
        "failure_action": "freeze_list_then_run_timing_census"
        if eligible_indices else "do_not_submit_remote_job",
    }
    write_json(os.path.join(output, "p4s4_pilot_evidence.json"), evidence)
    write_json(os.path.join(output, "p4s4_pilot_eligibility.json"),
               eligibility)
    count = p4c.write_checksums(output)
    print("P4-S4 pilot evidence complete: {} candidates eligible; "
          "{} ledger files".format(len(eligible_indices), count))


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
