from __future__ import print_function

import argparse
import math
import os

import generate_bounded_krylov_p4c_evidence as p4c
import generate_bounded_krylov_p4s3_pilot_evidence as p4s3
import generate_bounded_krylov_p4s6_pilot_evidence as p4s6


SCHEMA_VERSION = 5
EXPECTED_FIXTURE = "p4s7_partial_guide_exact_q_prefilter"
EXPECTED_CANDIDATES = ((1, 1), (1, 10), (1, 100),
                       (2, 1), (2, 10), (2, 100))
FAMILIES = p4s6.FAMILIES
ENTRY_COUNT = p4s6.ENTRY_COUNT
COMPONENTS = p4s6.COMPONENTS


read_json = p4c.read_json
write_json = p4c.write_json
sha256_file = p4c.sha256_file
parse_number = p4c.parse_number
require_close = p4s6.require_close
fixed_sector_size = p4c.fixed_sector_size


def parse_key_values(line, skip=1):
    return {key: parse_number(value)
            for key, value in p4c.parse_key_values(line, skip).items()}


def candidate_key(partial_order, alpha):
    return "m{}_alpha{}".format(int(partial_order), int(alpha))


def validate_policy(policy, allow_smoke_policy=False):
    if int(policy.get("schema_version", 0)) != SCHEMA_VERSION or \
            policy.get("policy_id") != \
            "power-lanczos-zero-support-p4-r0-v7-stage-a":
        raise AssertionError("P4-S7 policy identity mismatch")
    scope = policy.get("scope", {})
    expected_cases = [[4, 1], [4, 4], [6, 1], [6, 4], [8, 1], [8, 4]]
    if scope.get("site_counts") != [4, 6, 8] or \
            scope.get("qp_total") != [1, 4] or \
            scope.get("physical_case_order") != expected_cases or \
            int(scope.get("omp_threads", 0)) != 1 or \
            int(scope.get("blas_threads", 0)) != 1:
        raise AssertionError("P4-S7 scope mismatch")
    guide = policy.get("measurement_guide", {})
    if int(guide.get("order", 0)) != 3 or \
            float(guide.get("rho", 0.0)) != 0.01 or \
            guide.get("lambda") != [1.0, 1.0, 1.0, 1.0]:
        raise AssertionError("P4-S7 target guide mismatch")
    partial = policy.get("partial_guide", {})
    observed = tuple((int(item.get("m", 0)), int(item.get("alpha", 0)))
                     for item in partial.get("candidate_order", []))
    if observed != EXPECTED_CANDIDATES or \
            partial.get("formula") != \
            "sum_abs_uk_squared_k0_to_m_plus_alpha_eta_v1" or \
            partial.get("proposal") != \
            "exact_q_independence_nonproduction_prefilter" or \
            partial.get("grid_extension") != \
            "forbidden_under_v7_stage_a" or \
            partial.get("stage_b_mapping", {}).get("K_order") != [16, 32]:
        raise AssertionError("P4-S7 frontier mismatch")
    pilot = policy.get("exact_q_prefilter", {})
    if not allow_smoke_policy and \
            int(pilot.get("sample_count_per_seed", 0)) != 131072:
        raise AssertionError("P4-S7 sample count mismatch")
    if pilot.get("seed_hex_order") != [
            "0x5034533750494c31", "0x5034533750494c32",
            "0x5034533750494c33", "0x5034533750494c34"] or \
            int(pilot.get("predicted_official_sample_count", 0)) != 4096 or \
            float(pilot.get("maximum_predicted_se_budget_ratio", 0.0)) != \
            0.90 or float(pilot.get("maximum_tau_int", 0.0)) != 12.0 or \
            int(pilot.get(
                "exact_l4_balance_tolerance_multiplier_dbl_epsilon", 0
            )) != 256 or int(pilot.get(
                "anchor_tolerance_multiplier_dbl_epsilon", 0)) != 1024:
        raise AssertionError("P4-S7 prefilter gate mismatch")
    timing = policy.get("timing_census", {})
    if float(timing.get(
            "maximum_total_seconds_per_step_regression_fraction", -1)) != \
            0.10:
        raise AssertionError("P4-S7 timing gate mismatch")
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
        elif line.startswith("PILOT_PARTIAL_CANDIDATE "):
            candidates.append(parse_key_values(line))
        elif line.startswith("PILOT_SERIES "):
            series.append(parse_key_values(line))
        elif line.startswith("PILOT_ENTRY "):
            entries.append(parse_key_values(line))
        elif line.startswith("PILOT_SEED "):
            seeds.append(parse_key_values(line))
        elif line.startswith("PILOT_DECISION "):
            decisions.append(parse_key_values(line))
    if run is None or len(candidates) != len(EXPECTED_CANDIDATES) or \
            len(decisions) != len(EXPECTED_CANDIDATES):
        raise AssertionError("incomplete P4-S7 log {}".format(path))
    return {"path": os.path.abspath(path), "sha256": sha256_file(path),
            "run": run, "candidates": candidates, "series": series,
            "entries": entries, "seeds": seeds, "decisions": decisions,
            "line_count": len(text.splitlines())}


def validate_diagnostics(item, label):
    minimum = float(item["partial_minimum"])
    maximum = float(item["partial_maximum"])
    mean = float(item["partial_mean"])
    ratio_minimum = float(item["target_ratio_minimum"])
    ratio_maximum = float(item["target_ratio_maximum"])
    ratio_maximum_mean = float(item["target_ratio_maximum_mean_ratio"])
    ess = float(item["target_ratio_ess_fraction"])
    acceptance = float(item["ideal_acceptance"])
    change = float(item["ideal_configuration_change"])
    iid_ratio = float(item["maximum_iid_se_budget_ratio"])
    required_tau = float(item["required_tau_for_0p90"])
    if not (0.0 < minimum <= mean <= maximum) or \
            not (0.0 < ratio_minimum <= ratio_maximum) or \
            ratio_maximum_mean < 1.0 or not (0.0 < ess <= 1.0 + 1e-12) or \
            not (0.0 <= change <= acceptance <= 1.0 + 1e-12) or \
            iid_ratio <= 0.0:
        raise AssertionError("{} diagnostic bounds".format(label))
    require_close(required_tau, (0.90 / iid_ratio) ** 2,
                  label + " required tau")
    if item["iid_limiting_family"] not in FAMILIES or \
            int(item["iid_limiting_entry"]) not in range(ENTRY_COUNT):
        raise AssertionError("{} limiting entry".format(label))
    return {
        "partial_minimum": minimum, "partial_maximum": maximum,
        "partial_mean": mean, "target_ratio_minimum": ratio_minimum,
        "target_ratio_maximum": ratio_maximum,
        "target_ratio_maximum_mean_ratio": ratio_maximum_mean,
        "target_ratio_ess_fraction": ess,
        "ideal_acceptance": acceptance,
        "ideal_configuration_change": change,
        "maximum_iid_se_budget_ratio": iid_ratio,
        "required_tau_for_0p90": required_tau,
    }


def validate_candidate(run, index, policy, allow_smoke_policy=False):
    expected_m, expected_alpha = EXPECTED_CANDIDATES[index]
    pilot = policy["exact_q_prefilter"]
    header = run["candidates"][index]
    decision = run["decisions"][index]
    site = int(run["run"]["site_count"])
    qp = int(run["run"]["qp_total"])
    sector_size = int(run["run"]["sector_size"])
    sample_count = int(run["run"]["sample_count_per_seed"])
    if int(header["schema"]) != SCHEMA_VERSION or \
            int(header["candidate_index"]) != index or \
            int(decision["candidate_index"]) != index or \
            int(header["site_count"]) != site or \
            int(header["qp_total"]) != qp or \
            int(header["sample_count_per_seed"]) != sample_count or \
            int(header["partial_order"]) != expected_m or \
            int(header["floor_multiplier"]) != expected_alpha or \
            header["proposal"] != "exact_q_independence":
        raise AssertionError("P4-S7 candidate schedule mismatch")
    require_close(header["rho"], 0.01, "P4-S7 rho")
    if float(header["eta"]) <= 0.0 or \
            int(header["anchor_count"]) != sector_size or \
            not bool(int(header["exact_balance_pass"])) or \
            not bool(int(header["restart_pass"])) or \
            any(int(header[field]) == 0 for field in (
                "guide_policy_hash", "guide_shape_hash",
                "partial_policy_hash", "target_table_hash",
                "partial_table_hash")):
        raise AssertionError("P4-S7 identity/balance mismatch")
    anchor_tolerance = float(
        pilot["anchor_tolerance_multiplier_dbl_epsilon"]
    ) * 2.220446049250313e-16
    if float(header["maximum_anchor_relative_residual"]) > anchor_tolerance:
        raise AssertionError("P4-S7 raw anchor mismatch")
    diagnostics = validate_diagnostics(header,
                                       "P4-S7 candidate {}".format(index))
    if not allow_smoke_policy and sample_count != \
            int(pilot["sample_count_per_seed"]):
        raise AssertionError("P4-S7 formal sample count mismatch")
    candidate_series = [x for x in run["series"]
                        if int(x["candidate_index"]) == index]
    candidate_entries = [x for x in run["entries"]
                         if int(x["candidate_index"]) == index]
    candidate_seeds = [x for x in run["seeds"]
                       if int(x["candidate_index"]) == index]
    if len(candidate_series) != 4 * len(FAMILIES) * ENTRY_COUNT * 2 or \
            len(candidate_entries) != 4 * len(FAMILIES) * ENTRY_COUNT or \
            len(candidate_seeds) != 4:
        raise AssertionError("P4-S7 census mismatch")
    expected_series = {(seed, family, entry, component)
                       for seed in range(4) for family in FAMILIES
                       for entry in range(ENTRY_COUNT)
                       for component in COMPONENTS}
    observed_series = {(int(x["seed_index"]), x["family"], int(x["entry"]),
                        x["component"]) for x in candidate_series}
    if observed_series != expected_series:
        raise AssertionError("P4-S7 series identity mismatch")
    for seed_index, seed in enumerate(candidate_seeds):
        if int(seed["seed_index"]) != seed_index or \
                str(seed["seed_hex"]).lower() != \
                pilot["seed_hex_order"][seed_index].lower() or \
                int(seed["accepted"]) + int(seed["rejected"]) != \
                sample_count or \
                int(seed["independence_attempted"]) != sample_count or \
                not (0 <= int(seed["independence_self"]) <= sample_count) or \
                int(seed["trace_hash"]) == 0:
            raise AssertionError("P4-S7 seed census mismatch")
    maximum_tau = max(float(x["tau_int"]) for x in candidate_series)
    maximum_ratio = max(float(x["predicted_se_budget_ratio"])
                        for x in candidate_entries)
    require_close(decision["maximum_tau_int"], maximum_tau,
                  "P4-S7 max tau")
    require_close(decision["maximum_predicted_se_budget_ratio"],
                  maximum_ratio, "P4-S7 max ratio")
    if site == 4:
        tolerance = float(pilot[
            "exact_l4_balance_tolerance_multiplier_dbl_epsilon"
        ]) * 2.220446049250313e-16
        for field in ("row_residual", "db_residual",
                      "stationary_residual"):
            if float(header[field]) > tolerance:
                raise AssertionError("P4-S7 exact L4 {}".format(field))
    eligible = maximum_ratio <= float(
        pilot["maximum_predicted_se_budget_ratio"]) and \
        maximum_tau <= float(pilot["maximum_tau_int"]) and \
        bool(int(decision["table_pass"])) and \
        bool(int(decision["balance_pass"])) and \
        bool(int(decision["restart_pass"]))
    if bool(int(decision["physical_case_eligible"])) != eligible:
        raise AssertionError("P4-S7 eligibility mismatch")
    diagnostics.update({
        "candidate_index": index, "candidate": candidate_key(*EXPECTED_CANDIDATES[index]),
        "m": expected_m, "alpha": expected_alpha,
        "physical_case_eligible": eligible,
        "maximum_tau_int": maximum_tau,
        "maximum_predicted_se_budget_ratio": maximum_ratio,
        "guide_policy_hash": int(header["guide_policy_hash"]),
        "guide_shape_hash": int(header["guide_shape_hash"]),
        "partial_policy_hash": int(header["partial_policy_hash"]),
        "target_table_hash": int(header["target_table_hash"]),
        "partial_table_hash": int(header["partial_table_hash"]),
        "seed_summaries": candidate_seeds,
    })
    return diagnostics


def validate_run(run, policy, allow_smoke_policy=False):
    header = run["run"]
    if int(header.get("schema", 0)) != SCHEMA_VERSION or \
            header.get("fixture") != EXPECTED_FIXTURE:
        raise AssertionError("P4-S7 run header mismatch")
    site = int(header["site_count"])
    qp = int(header["qp_total"])
    sector_size = fixed_sector_size(site)
    if [site, qp] not in policy["scope"]["physical_case_order"] or \
            int(header["sector_size"]) != sector_size or \
            int(header["anchor_count"]) != sector_size:
        raise AssertionError("P4-S7 physical case mismatch")
    if not allow_smoke_policy and int(header["cache_requested_bytes"]) != \
            int(policy["scope"]["statistical_cache_bytes"]):
        raise AssertionError("P4-S7 cache mismatch")
    candidates = [validate_candidate(run, i, policy, allow_smoke_policy)
                  for i in range(len(EXPECTED_CANDIDATES))]
    for field in ("guide_policy_hash", "guide_shape_hash",
                  "target_table_hash"):
        if len(set(x[field] for x in candidates)) != 1:
            raise AssertionError("P4-S7 target changed: {}".format(field))
    if len(set(x["partial_policy_hash"] for x in candidates)) != 6 or \
            len(set(x["partial_table_hash"] for x in candidates)) != 6:
        raise AssertionError("P4-S7 partial identity collision")
    return {"site_count": site, "qp_total": qp,
            "rank_count": int(header["rank_count"]),
            "sample_count_per_seed": int(header["sample_count_per_seed"]),
            "cache_requested_bytes": int(header["cache_requested_bytes"]),
            "sector_size": sector_size, "plan_hash": int(header["plan_hash"]),
            "candidates": candidates}


def analyze_command(args):
    policy = validate_policy(read_json(args.policy), args.allow_smoke_policy)
    parsed = [parse_pilot_log(path) for path in args.pilot_log]
    observed = [[int(x["run"]["site_count"]), int(x["run"]["qp_total"])]
                for x in parsed]
    if sorted(observed) != sorted(policy["scope"]["physical_case_order"]) or \
            len(set(tuple(x) for x in observed)) != len(parsed):
        raise AssertionError("P4-S7 requires six unique cases")
    validated = {(int(x["run"]["site_count"]), int(x["run"]["qp_total"])):
                 validate_run(x, policy, args.allow_smoke_policy)
                 for x in parsed}
    output = os.path.abspath(args.output)
    if not os.path.isdir(output):
        os.makedirs(output)
    raw_paths = p4s3.materialize_logs(parsed, output)
    cases = []
    for site, qp in policy["scope"]["physical_case_order"]:
        parsed_run = next(x for x in parsed
                          if int(x["run"]["site_count"]) == site and
                          int(x["run"]["qp_total"]) == qp)
        summary = validated[(site, qp)]
        summary["raw_log"] = raw_paths[parsed_run["path"]]
        summary["raw_log_sha256"] = parsed_run["sha256"]
        summary["raw_log_line_count"] = parsed_run["line_count"]
        cases.append(summary)
    candidate_summaries = []
    eligible_indices = []
    for index, (partial_order, alpha) in enumerate(EXPECTED_CANDIDATES):
        candidate_cases = [case["candidates"][index] for case in cases]
        eligible = all(x["physical_case_eligible"] for x in candidate_cases)
        summary = {
            "candidate_index": index,
            "candidate": candidate_key(partial_order, alpha),
            "m": partial_order, "alpha": alpha,
            "prefilter_eligible": eligible,
            "maximum_tau_int": max(x["maximum_tau_int"]
                                   for x in candidate_cases),
            "maximum_predicted_se_budget_ratio": max(
                x["maximum_predicted_se_budget_ratio"]
                for x in candidate_cases),
            "maximum_target_ratio_maximum_mean_ratio": max(
                x["target_ratio_maximum_mean_ratio"]
                for x in candidate_cases),
            "minimum_target_ratio_ess_fraction": min(
                x["target_ratio_ess_fraction"] for x in candidate_cases),
            "minimum_ideal_acceptance": min(
                x["ideal_acceptance"] for x in candidate_cases),
            "cases": candidate_cases,
        }
        candidate_summaries.append(summary)
        if eligible:
            eligible_indices.append(index)
    generated_at = p4c.generated_at_jst()
    policy_sha = sha256_file(args.policy)
    evidence = {
        "schema_version": SCHEMA_VERSION,
        "artifact": "p4s7_partial_guide_exact_q_prefilter",
        "generated_at_jst": generated_at, "source_commit": args.source_commit,
        "producer_binary": os.path.basename(args.producer_binary),
        "producer_binary_sha256": sha256_file(args.producer_binary),
        "measurement_policy_sha256": policy_sha,
        "measurement_policy": policy, "physical_cases": cases,
        "candidate_summaries": candidate_summaries,
    }
    stage_b = []
    for index in eligible_indices:
        for step_count in policy["partial_guide"]["stage_b_mapping"]["K_order"]:
            stage_b.append({"stage_a_candidate_index": index,
                            "candidate": candidate_summaries[index]["candidate"],
                            "m": candidate_summaries[index]["m"],
                            "alpha": candidate_summaries[index]["alpha"],
                            "K": int(step_count)})
    eligibility = {
        "schema_version": SCHEMA_VERSION,
        "artifact": "p4s7_frozen_stage_a_eligibility",
        "generated_at_jst": generated_at, "source_commit": args.source_commit,
        "measurement_policy_sha256": policy_sha,
        "candidate_order": [candidate_key(*x) for x in EXPECTED_CANDIDATES],
        "eligible_candidate_indices": eligible_indices,
        "eligible_candidates": [candidate_summaries[i]["candidate"]
                                for i in eligible_indices],
        "stage_b_candidate_order": stage_b,
        "prefilter_decision": "ELIGIBLE" if eligible_indices else "STOP",
        "failure_action": "implement_only_frozen_stage_b_mapping"
        if eligible_indices else "do_not_implement_stage_b_or_run_timing",
    }
    write_json(os.path.join(output, "p4s7_prefilter_evidence.json"), evidence)
    write_json(os.path.join(output, "p4s7_prefilter_eligibility.json"),
               eligibility)
    count = p4c.write_checksums(output)
    print("P4-S7 prefilter complete: {} partial guides eligible; {} ledger "
          "files".format(len(eligible_indices), count))


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
