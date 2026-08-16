from __future__ import print_function

import argparse
import math
import os

import generate_bounded_krylov_p4c_evidence as p4c
import generate_bounded_krylov_p4s3_pilot_evidence as p4s3


SCHEMA_VERSION = 4
EXPECTED_FIXTURE = "p4s6_surrogate_exact_table"
FAMILIES = ("S", "K", "B")
ENTRY_COUNT = 6
COMPONENTS = ("real", "imag")
EXPECTED_CANDIDATES = (
    (4, 1), (4, 10), (4, 100),
    (16, 1), (16, 10), (16, 100),
    (32, 1), (32, 10), (32, 100),
)
EXPECTED_ALPHAS = (1, 10, 100)


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


def candidate_key(step_count, alpha):
    return "K{}_alpha{}".format(int(step_count), int(alpha))


def validate_policy(policy, allow_smoke_policy=False):
    if int(policy.get("schema_version", 0)) != SCHEMA_VERSION or \
            policy.get("policy_id") != \
            "power-lanczos-zero-support-p4-r0-v6":
        raise AssertionError("P4-S6 policy identity mismatch")
    scope = policy.get("scope", {})
    expected_cases = [[4, 1], [4, 4], [6, 1], [6, 4], [8, 1], [8, 4]]
    if scope.get("site_counts") != [4, 6, 8] or \
            scope.get("qp_total") != [1, 4] or \
            scope.get("physical_case_order") != expected_cases or \
            int(scope.get("omp_threads", 0)) != 1 or \
            int(scope.get("blas_threads", 0)) != 1:
        raise AssertionError("P4-S6 physical scope mismatch")
    guide = policy.get("measurement_guide", {})
    if int(guide.get("order", 0)) != 3 or \
            float(guide.get("rho", 0.0)) != 0.01 or \
            guide.get("lambda") != [1.0, 1.0, 1.0, 1.0]:
        raise AssertionError("P4-S6 measurement guide mismatch")
    surrogate = policy.get("surrogate", {})
    candidates = tuple(
        (int(item.get("K", 0)), int(item.get("alpha", 0)))
        for item in surrogate.get("candidate_order", []))
    if candidates != EXPECTED_CANDIDATES or \
            surrogate.get("formula") != \
            "Q_alpha_equals_abs_u0_squared_plus_alpha_eta_v1" or \
            surrogate.get("inner_kernel") != \
            "symmetric_fixed_sector_neighbor_mh" or \
            surrogate.get("grid_extension") != "forbidden_under_v6" or \
            int(surrogate.get("policy_version", 0)) != 1 or \
            int(surrogate.get("transaction_version", 0)) != 1:
        raise AssertionError("P4-S6 surrogate frontier mismatch")
    controls = policy.get("deterministic_controls", [])
    if len(controls) != 3 or any(
            item.get("control_id") != "exact_q_independence" or
            int(item.get("K", 0)) != 1 or
            int(item.get("alpha", 0)) != EXPECTED_ALPHAS[index] or
            bool(item.get("promotable"))
            for index, item in enumerate(controls)):
        raise AssertionError("P4-S6 control frontier mismatch")
    pilot = policy.get("exact_table_pilot", {})
    if not allow_smoke_policy and \
            int(pilot.get("sample_count_per_seed", 0)) != 131072:
        raise AssertionError("P4-S6 pilot sample count mismatch")
    if pilot.get("seed_hex_order") != [
            "0x5034533650494c31", "0x5034533650494c32",
            "0x5034533650494c33", "0x5034533650494c34"] or \
            float(pilot.get("maximum_predicted_se_budget_ratio", 0.0)) != \
            0.90 or float(pilot.get("maximum_tau_int", 0.0)) != 12.0 or \
            int(pilot.get("predicted_official_sample_count", 0)) != 4096 or \
            int(pilot.get(
                "exact_l4_balance_tolerance_multiplier_dbl_epsilon", 0
            )) != 256 or int(pilot.get(
                "exact_l4_powered_ratio_tolerance_multiplier_dbl_epsilon",
                0)) != 1024:
        raise AssertionError("P4-S6 pilot gate mismatch")
    if policy.get("official_markov", {}).get("official_seed_hex") != \
            "0x5034533652305631":
        raise AssertionError("P4-S6 official seed mismatch")
    timing = policy.get("timing_census", {})
    if int(timing.get("sample_count", 0)) != 128 or \
            int(timing.get("repeat_count", 0)) != 7 or \
            float(timing.get(
                "maximum_total_seconds_per_step_regression_fraction", -1
            )) != 0.10:
        raise AssertionError("P4-S6 timing gate mismatch")
    return policy


def parse_pilot_log(path):
    run = None
    controls = []
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
        elif line.startswith("PILOT_SURROGATE_CONTROL "):
            controls.append(parse_key_values(line))
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
    if run is None or len(controls) != len(EXPECTED_ALPHAS) or \
            len(candidates) != len(EXPECTED_CANDIDATES) or \
            len(decisions) != len(EXPECTED_CANDIDATES):
        raise AssertionError("incomplete P4-S6 pilot log {}".format(path))
    return {
        "path": os.path.abspath(path),
        "sha256": sha256_file(path),
        "run": run,
        "controls": controls,
        "candidates": candidates,
        "series": series,
        "entries": entries,
        "seeds": seeds,
        "decisions": decisions,
        "line_count": len(text.splitlines()),
    }


def validate_surrogate_diagnostics(item, label):
    minimum = float(item["surrogate_minimum"])
    maximum = float(item["surrogate_maximum"])
    mean = float(item["surrogate_mean"])
    ratio_minimum = float(item["target_ratio_minimum"])
    ratio_maximum = float(item["target_ratio_maximum"])
    ratio_maximum_mean = float(item["target_ratio_maximum_mean_ratio"])
    ess_fraction = float(item["target_ratio_ess_fraction"])
    iid_ratio = float(item["maximum_iid_se_budget_ratio"])
    required_tau = float(item["required_tau_for_0p90"])
    if not (0.0 < minimum <= mean <= maximum) or \
            not (0.0 < ratio_minimum <= ratio_maximum) or \
            ratio_maximum_mean < 1.0 or \
            not (0.0 < ess_fraction <= 1.0 + 1.0e-12) or iid_ratio <= 0.0:
        raise AssertionError("{} diagnostic bounds".format(label))
    require_close(required_tau, (0.90 / iid_ratio) ** 2,
                  label + " required tau")
    family = item.get("iid_limiting_family", item.get("limiting_family"))
    entry = int(item.get("iid_limiting_entry",
                         item.get("limiting_entry", -1)))
    if family not in FAMILIES or entry not in range(ENTRY_COUNT):
        raise AssertionError("{} IID limiting entry".format(label))
    return {
        "surrogate_minimum": minimum,
        "surrogate_maximum": maximum,
        "surrogate_mean": mean,
        "target_ratio_minimum": ratio_minimum,
        "target_ratio_maximum": ratio_maximum,
        "target_ratio_maximum_mean_ratio": ratio_maximum_mean,
        "target_ratio_ess_fraction": ess_fraction,
        "maximum_iid_se_budget_ratio": iid_ratio,
        "required_tau_for_0p90": required_tau,
        "iid_limiting_family": family,
        "iid_limiting_entry": entry,
    }


def validate_controls(run, policy):
    sector_size = int(run["run"]["sector_size"])
    results = []
    for index, item in enumerate(run["controls"]):
        expected = policy["deterministic_controls"][index]
        if int(item["schema"]) != SCHEMA_VERSION or \
                int(item["control_index"]) != index or \
                item["control_id"] != expected["control_id"] or \
                bool(int(item["promotable"])) or \
                int(item["alpha"]) != int(expected["alpha"]) or \
                int(item["anchor_count"]) != sector_size or \
                float(item["eta"]) <= 0.0 or \
                any(int(item[field]) == 0 for field in (
                    "guide_policy_hash", "guide_shape_hash",
                    "surrogate_policy_hash", "target_table_hash",
                    "surrogate_table_hash")):
            raise AssertionError("P4-S6 control identity mismatch")
        require_close(item["rho"], 0.01, "P4-S6 control rho")
        acceptance = float(item["ideal_acceptance"])
        change = float(item["ideal_configuration_change"])
        if not (0.0 <= change <= acceptance <= 1.0 + 1.0e-12):
            raise AssertionError("P4-S6 exact-Q control bounds")
        summary = validate_surrogate_diagnostics(
            item, "P4-S6 control {}".format(index))
        summary.update({
            "control_index": index,
            "control_id": expected["control_id"],
            "alpha": int(expected["alpha"]),
            "ideal_acceptance": acceptance,
            "ideal_configuration_change": change,
            "guide_policy_hash": int(item["guide_policy_hash"]),
            "guide_shape_hash": int(item["guide_shape_hash"]),
            "surrogate_policy_hash": int(item["surrogate_policy_hash"]),
            "target_table_hash": int(item["target_table_hash"]),
            "surrogate_table_hash": int(item["surrogate_table_hash"]),
        })
        results.append(summary)
    return results


def validate_candidate(run, candidate_index, policy,
                       allow_smoke_policy=False):
    expected_k, expected_alpha = EXPECTED_CANDIDATES[candidate_index]
    pilot = policy["exact_table_pilot"]
    header = run["candidates"][candidate_index]
    decision = run["decisions"][candidate_index]
    site = int(run["run"]["site_count"])
    qp = int(run["run"]["qp_total"])
    sector_size = int(run["run"]["sector_size"])
    sample_count = int(run["run"]["sample_count_per_seed"])
    if int(header["schema"]) != SCHEMA_VERSION or \
            int(header["candidate_index"]) != candidate_index or \
            int(decision["candidate_index"]) != candidate_index or \
            int(header["site_count"]) != site or \
            int(header["qp_total"]) != qp or \
            int(header["sample_count_per_seed"]) != sample_count or \
            int(header["surrogate_step_count"]) != expected_k or \
            int(header["floor_multiplier"]) != expected_alpha or \
            header["inner_kernel"] != "neighbor_only":
        raise AssertionError("P4-S6 candidate schedule mismatch")
    require_close(header["rho"], 0.01, "P4-S6 candidate rho")
    if float(header["eta"]) <= 0.0 or \
            int(header["anchor_count"]) != sector_size or \
            not bool(int(header["exact_balance_pass"])) or \
            not bool(int(header["restart_pass"])) or \
            any(int(header[field]) == 0 for field in (
                "proposal_policy_hash", "proposal_model_hash",
                "guide_policy_hash", "guide_shape_hash",
                "surrogate_policy_hash", "target_table_hash",
                "surrogate_table_hash")):
        raise AssertionError("P4-S6 candidate identity/balance mismatch")
    diagnostics = validate_surrogate_diagnostics(
        header, "P4-S6 candidate {}".format(candidate_index))
    if not allow_smoke_policy and sample_count != \
            int(pilot["sample_count_per_seed"]):
        raise AssertionError("P4-S6 official pilot sample count mismatch")

    candidate_series = [item for item in run["series"]
                        if int(item["candidate_index"]) == candidate_index]
    candidate_entries = [item for item in run["entries"]
                         if int(item["candidate_index"]) == candidate_index]
    candidate_seeds = [item for item in run["seeds"]
                       if int(item["candidate_index"]) == candidate_index]
    if len(candidate_series) != 4 * len(FAMILIES) * ENTRY_COUNT * 2 or \
            len(candidate_entries) != 4 * len(FAMILIES) * ENTRY_COUNT or \
            len(candidate_seeds) != 4:
        raise AssertionError("P4-S6 candidate census mismatch")
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
        raise AssertionError("P4-S6 series identity mismatch")
    for seed_index, seed in enumerate(candidate_seeds):
        inner_attempted = sample_count * expected_k
        if int(seed["seed_index"]) != seed_index or \
                str(seed["seed_hex"]).lower() != \
                pilot["seed_hex_order"][seed_index].lower() or \
                int(seed["accepted"]) + int(seed["rejected"]) != \
                sample_count or \
                int(seed["surrogate_inner_attempted"]) != inner_attempted or \
                int(seed["surrogate_inner_accepted"]) + \
                int(seed["surrogate_inner_rejected"]) != inner_attempted or \
                int(seed["surrogate_final_changed"]) + \
                int(seed["surrogate_final_self"]) != sample_count or \
                int(seed["trace_hash"]) == 0:
            raise AssertionError("P4-S6 seed/chain census mismatch")
    maximum_tau = max(float(item["tau_int"]) for item in candidate_series)
    maximum_ratio = max(float(item["predicted_se_budget_ratio"])
                        for item in candidate_entries)
    require_close(decision["maximum_tau_int"], maximum_tau,
                  "P4-S6 maximum tau")
    require_close(decision["maximum_predicted_se_budget_ratio"],
                  maximum_ratio, "P4-S6 maximum ratio")
    if site == 4:
        epsilon = 2.220446049250313e-16
        tolerance = float(pilot[
            "exact_l4_balance_tolerance_multiplier_dbl_epsilon"]) * epsilon
        powered_ratio_tolerance = float(pilot[
            "exact_l4_powered_ratio_tolerance_multiplier_dbl_epsilon"
        ]) * epsilon
        for field in ("inner_row_residual", "inner_db_residual",
                      "powered_row_residual", "powered_db_residual",
                      "outer_db_residual", "outer_stationary_residual"):
            if float(header[field]) > tolerance:
                raise AssertionError("P4-S6 exact L4 {}".format(field))
        if float(header["powered_ratio_residual"]) > powered_ratio_tolerance:
            raise AssertionError("P4-S6 exact L4 powered ratio")
    expected_eligible = (
        maximum_ratio <= float(pilot["maximum_predicted_se_budget_ratio"])
        and maximum_tau <= float(pilot["maximum_tau_int"])
        and bool(int(decision["table_pass"]))
        and bool(int(decision["balance_pass"]))
        and bool(int(decision["restart_pass"]))
    )
    if bool(int(decision["physical_case_eligible"])) != expected_eligible:
        raise AssertionError("P4-S6 physical eligibility mismatch")
    diagnostics.update({
        "candidate_index": candidate_index,
        "candidate": candidate_key(expected_k, expected_alpha),
        "K": expected_k,
        "alpha": expected_alpha,
        "physical_case_eligible": expected_eligible,
        "maximum_tau_int": maximum_tau,
        "maximum_predicted_se_budget_ratio": maximum_ratio,
        "proposal_policy_hash": int(header["proposal_policy_hash"]),
        "proposal_model_hash": int(header["proposal_model_hash"]),
        "guide_policy_hash": int(header["guide_policy_hash"]),
        "guide_shape_hash": int(header["guide_shape_hash"]),
        "surrogate_policy_hash": int(header["surrogate_policy_hash"]),
        "target_table_hash": int(header["target_table_hash"]),
        "surrogate_table_hash": int(header["surrogate_table_hash"]),
        "anchor_count": int(header["anchor_count"]),
        "seed_summaries": candidate_seeds,
    })
    return diagnostics


def validate_run(run, policy, allow_smoke_policy=False):
    header = run["run"]
    if int(header.get("schema", 0)) != SCHEMA_VERSION or \
            header.get("fixture") != EXPECTED_FIXTURE:
        raise AssertionError("P4-S6 run header mismatch")
    site = int(header["site_count"])
    qp = int(header["qp_total"])
    sector_size = fixed_sector_size(site)
    if [site, qp] not in policy["scope"]["physical_case_order"] or \
            int(header["sector_size"]) != sector_size or \
            int(header["anchor_count"]) != sector_size:
        raise AssertionError("P4-S6 physical case/anchor mismatch")
    if not allow_smoke_policy and \
            int(header["cache_requested_bytes"]) != \
            int(policy["scope"]["statistical_cache_bytes"]):
        raise AssertionError("P4-S6 cache policy mismatch")
    controls = validate_controls(run, policy)
    candidates = [
        validate_candidate(run, index, policy, allow_smoke_policy)
        for index in range(len(EXPECTED_CANDIDATES))
    ]
    if len(set(item["proposal_policy_hash"] for item in candidates)) != 1 or \
            len(set(item["proposal_model_hash"] for item in candidates)) != 1:
        raise AssertionError("P4-S6 neighbor proposal identity changed")
    for field in ("guide_policy_hash", "guide_shape_hash",
                  "target_table_hash"):
        if len(set(item[field] for item in candidates)) != 1:
            raise AssertionError("P4-S6 target identity changed: {}".format(
                field))
    if len(set(item["surrogate_policy_hash"] for item in candidates)) != \
            len(EXPECTED_CANDIDATES):
        raise AssertionError("P4-S6 surrogate policy hash collision")
    for alpha_index, alpha in enumerate(EXPECTED_ALPHAS):
        alpha_candidates = [item for item in candidates
                            if item["alpha"] == alpha]
        if len(set(item["surrogate_table_hash"]
                   for item in alpha_candidates)) != 1 or \
                controls[alpha_index]["surrogate_table_hash"] != \
                alpha_candidates[0]["surrogate_table_hash"]:
            raise AssertionError("P4-S6 Q table mismatch for alpha {}".
                                 format(alpha))
    if len(set(item["surrogate_table_hash"] for item in controls)) != 3 or \
            any(control["guide_policy_hash"] !=
                candidates[0]["guide_policy_hash"] or
                control["guide_shape_hash"] !=
                candidates[0]["guide_shape_hash"] or
                control["target_table_hash"] !=
                candidates[0]["target_table_hash"]
                for control in controls):
        raise AssertionError("P4-S6 control/target identity mismatch")
    etas = [float(item["eta"]) for item in run["controls"] +
            run["candidates"]]
    for eta in etas[1:]:
        require_close(eta, etas[0], "P4-S6 eta")
    return {
        "site_count": site,
        "qp_total": qp,
        "rank_count": int(header["rank_count"]),
        "sample_count_per_seed": int(header["sample_count_per_seed"]),
        "cache_requested_bytes": int(header["cache_requested_bytes"]),
        "sector_size": int(header["sector_size"]),
        "plan_hash": int(header["plan_hash"]),
        "controls": controls,
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
        raise AssertionError("P4-S6 evidence requires six unique cases")
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
    for candidate_index, (step_count, alpha) in enumerate(
            EXPECTED_CANDIDATES):
        cases = [case["candidates"][candidate_index]
                 for case in ordered_cases]
        eligible = all(case["physical_case_eligible"] for case in cases)
        summary = {
            "candidate_index": candidate_index,
            "candidate": candidate_key(step_count, alpha),
            "K": step_count,
            "alpha": alpha,
            "pilot_eligible": eligible,
            "maximum_tau_int": max(case["maximum_tau_int"]
                                   for case in cases),
            "maximum_predicted_se_budget_ratio": max(
                case["maximum_predicted_se_budget_ratio"]
                for case in cases),
            "maximum_iid_se_budget_ratio": max(
                case["maximum_iid_se_budget_ratio"] for case in cases),
            "minimum_required_tau_for_0p90": min(
                case["required_tau_for_0p90"] for case in cases),
            "maximum_target_ratio_maximum_mean_ratio": max(
                case["target_ratio_maximum_mean_ratio"] for case in cases),
            "minimum_target_ratio_ess_fraction": min(
                case["target_ratio_ess_fraction"] for case in cases),
            "cases": cases,
        }
        candidate_summaries.append(summary)
        if eligible:
            eligible_indices.append(candidate_index)
    generated_at = p4c.generated_at_jst()
    evidence = {
        "schema_version": SCHEMA_VERSION,
        "artifact": "p4s6_surrogate_exact_table_pilot",
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
        "artifact": "p4s6_frozen_pilot_eligibility",
        "generated_at_jst": generated_at,
        "source_commit": args.source_commit,
        "measurement_policy_sha256": evidence["measurement_policy_sha256"],
        "candidate_order": [candidate_key(*item)
                            for item in EXPECTED_CANDIDATES],
        "eligible_candidate_indices": eligible_indices,
        "eligible_candidates": [candidate_summaries[index]["candidate"]
                                for index in eligible_indices],
        "preferred_candidate_index": eligible_indices[0]
        if eligible_indices else None,
        "pilot_decision": "ELIGIBLE" if eligible_indices else "STOP",
        "failure_action": "freeze_list_then_run_real_callback_timing_census"
        if eligible_indices else "do_not_run_timing_or_submit_remote_job",
    }
    write_json(os.path.join(output, "p4s6_pilot_evidence.json"), evidence)
    write_json(os.path.join(output, "p4s6_pilot_eligibility.json"),
               eligibility)
    count = p4c.write_checksums(output)
    print("P4-S6 pilot evidence complete: {} candidates eligible; "
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
