from __future__ import print_function

import argparse
import math
import os

import generate_bounded_krylov_p4c_evidence as p4c
import generate_bounded_krylov_p4s3_pilot_evidence as p4s3


SCHEMA_VERSION = 3
EXPECTED_FIXTURE = "p4s5_flat_guide_exact_table"
FAMILIES = ("S", "K", "B")
ENTRY_COUNT = 6
COMPONENTS = ("real", "imag")
CANDIDATE_COUNT = 7
CONTROL_IDS = ("rho_0p01", "uniform_guide_limit")


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


def candidate_key(rho):
    return "rho_{:.12g}".format(float(rho))


def validate_policy(policy, allow_smoke_policy=False):
    if int(policy.get("schema_version", 0)) != SCHEMA_VERSION or \
            policy.get("policy_id") != \
            "power-lanczos-zero-support-p4-r0-v5":
        raise AssertionError("P4-S5 policy identity mismatch")
    scope = policy.get("scope", {})
    expected_cases = [[4, 1], [4, 4], [6, 1], [6, 4], [8, 1], [8, 4]]
    if scope.get("site_counts") != [4, 6, 8] or \
            scope.get("qp_total") != [1, 4] or \
            scope.get("physical_case_order") != expected_cases:
        raise AssertionError("P4-S5 physical scope mismatch")
    guide = policy.get("guide", {})
    observed_rho = [float(item) for item in
                    guide.get("candidate_rho_order", [])]
    expected_rho = [0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0]
    if observed_rho != expected_rho or int(guide.get("order", 0)) != 3 or \
            guide.get("lambda") != [1.0, 1.0, 1.0, 1.0] or \
            guide.get("grid_extension") != "forbidden_under_v5":
        raise AssertionError("P4-S5 guide frontier mismatch")
    proposal = policy.get("proposal", {})
    if proposal.get("kernel") != "neighbor_global_mixture" or \
            int(proposal.get("policy_version", 0)) != 2 or \
            int(proposal.get("global_numerator", 0)) != 1 or \
            int(proposal.get("global_denominator", 0)) != 1 or \
            int(proposal.get("neighbor_numerator", -1)) != 0 or \
            int(proposal.get("neighbor_denominator", -1)) != 0:
        raise AssertionError("P4-S5 proposal mismatch")
    controls = policy.get("deterministic_controls", [])
    if len(controls) != 2 or \
            [item.get("control_id") for item in controls] != \
            list(CONTROL_IDS) or \
            float(controls[0].get("rho", 0.0)) != 0.01 or \
            bool(controls[0].get("uniform_guide")) or \
            not bool(controls[1].get("uniform_guide")) or \
            controls[1].get("rho") is not None or \
            any(bool(item.get("promotable")) for item in controls):
        raise AssertionError("P4-S5 deterministic controls mismatch")
    pilot = policy.get("exact_table_pilot", {})
    if not allow_smoke_policy and \
            int(pilot.get("sample_count_per_seed", 0)) != 131072:
        raise AssertionError("P4-S5 pilot sample count mismatch")
    if pilot.get("seed_hex_order") != [
            "0x5034533550494c31", "0x5034533550494c32",
            "0x5034533550494c33", "0x5034533550494c34"]:
        raise AssertionError("P4-S5 seed order mismatch")
    if float(pilot.get("maximum_predicted_se_budget_ratio", 0.0)) != 0.90 or \
            float(pilot.get("maximum_tau_int", 0.0)) != 12.0 or \
            int(pilot.get("predicted_official_sample_count", 0)) != 4096:
        raise AssertionError("P4-S5 pilot gate mismatch")
    if policy.get("official_markov", {}).get("official_seed_hex") != \
            "0x5034533552305631":
        raise AssertionError("P4-S5 official seed mismatch")
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
        elif line.startswith("PILOT_CONTROL "):
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
    if run is None or len(controls) != 2 or \
            len(candidates) != CANDIDATE_COUNT or \
            len(decisions) != CANDIDATE_COUNT:
        raise AssertionError("incomplete P4-S5 pilot log {}".format(path))
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


def validate_diagnostics(item, sector_size, label):
    minimum = float(item["guide_minimum"])
    maximum = float(item["guide_maximum"])
    mean = float(item["guide_mean"])
    concentration = float(item["guide_maximum_minimum_ratio"])
    weight_ratio = float(item["importance_weight_maximum_mean_ratio"])
    ess_fraction = float(item["importance_weight_ess_fraction"])
    iid_ratio = float(item["maximum_iid_se_budget_ratio"])
    required_tau = float(item["required_tau_for_0p90"])
    if not (0.0 < minimum <= mean <= maximum) or \
            weight_ratio < 1.0 or \
            not (0.0 < ess_fraction <= 1.0 + 1.0e-12) or iid_ratio <= 0.0:
        raise AssertionError("{} diagnostic bounds".format(label))
    require_close(concentration, maximum / minimum,
                  label + " guide concentration")
    require_close(item["uniform_sector_self_probability"],
                  1.0 / float(sector_size), label + " uniform self")
    require_close(required_tau, (0.90 / iid_ratio) ** 2,
                  label + " required tau")
    if item.get("limiting_family", item.get("iid_limiting_family")) not in \
            FAMILIES or int(item.get(
                "limiting_entry", item.get("iid_limiting_entry", -1))) not in \
            range(ENTRY_COUNT):
        raise AssertionError("{} IID limiting entry".format(label))
    return {
        "guide_minimum": minimum,
        "guide_maximum": maximum,
        "guide_mean": mean,
        "guide_maximum_minimum_ratio": concentration,
        "importance_weight_maximum_mean_ratio": weight_ratio,
        "importance_weight_ess_fraction": ess_fraction,
        "maximum_iid_se_budget_ratio": iid_ratio,
        "required_tau_for_0p90": required_tau,
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
                bool(int(item["uniform_guide"])) != \
                bool(expected["uniform_guide"]) or \
                int(item["guide_policy_hash"]) == 0 or \
                int(item["guide_shape_hash"]) == 0 or \
                int(item["table_hash"]) == 0 or \
                int(item["anchor_count"]) != sector_size:
            raise AssertionError("P4-S5 control identity mismatch")
        if expected["rho"] is None:
            require_close(item["rho"], 0.0, "uniform control rho")
            for field in ("guide_minimum", "guide_maximum", "guide_mean",
                          "guide_maximum_minimum_ratio"):
                require_close(item[field], 1.0,
                              "uniform control " + field)
        else:
            require_close(item["rho"], expected["rho"], "rho control")
        summary = validate_diagnostics(item, sector_size,
                                       "P4-S5 control {}".format(index))
        summary.update({
            "control_index": index,
            "control_id": expected["control_id"],
            "uniform_guide": bool(expected["uniform_guide"]),
            "rho": expected["rho"],
            "guide_policy_hash": int(item["guide_policy_hash"]),
            "guide_shape_hash": int(item["guide_shape_hash"]),
            "table_hash": int(item["table_hash"]),
        })
        results.append(summary)
    return results


def validate_candidate(run, candidate_index, policy,
                       allow_smoke_policy=False):
    rho = float(policy["guide"]["candidate_rho_order"][candidate_index])
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
            int(header["sample_count_per_seed"]) != sample_count:
        raise AssertionError("P4-S5 candidate schedule mismatch")
    require_close(header["rho"], rho, "P4-S5 candidate rho")
    if int(header["global_numerator"]) != 1 or \
            int(header["global_denominator"]) != 1 or \
            int(header["neighbor_numerator"]) != 0 or \
            int(header["neighbor_denominator"]) != 0 or \
            any(int(header[field]) != 0 for field in (
                "distance_numerator", "distance_denominator",
                "distance_rounding_rule", "shell_max_distance",
                "shell_distance", "shell_count")):
        raise AssertionError("P4-S5 global-only proposal mismatch")
    for field in ("proposal_policy_hash", "proposal_model_hash",
                  "guide_policy_hash", "guide_shape_hash", "table_hash"):
        if int(header[field]) == 0:
            raise AssertionError("P4-S5 zero {}".format(field))
    if int(header["anchor_count"]) != sector_size or \
            not bool(int(header["exact_balance_pass"])) or \
            not bool(int(header["restart_pass"])):
        raise AssertionError("P4-S5 anchor/balance/restart mismatch")
    diagnostics = validate_diagnostics(
        header, sector_size, "P4-S5 candidate {}".format(candidate_index))
    if not allow_smoke_policy and sample_count != \
            int(pilot["sample_count_per_seed"]):
        raise AssertionError("P4-S5 official pilot sample count mismatch")

    candidate_series = [item for item in run["series"]
                        if int(item["candidate_index"]) == candidate_index]
    candidate_entries = [item for item in run["entries"]
                         if int(item["candidate_index"]) == candidate_index]
    candidate_seeds = [item for item in run["seeds"]
                       if int(item["candidate_index"]) == candidate_index]
    if len(candidate_series) != 4 * len(FAMILIES) * ENTRY_COUNT * 2 or \
            len(candidate_entries) != 4 * len(FAMILIES) * ENTRY_COUNT or \
            len(candidate_seeds) != 4:
        raise AssertionError("P4-S5 candidate census mismatch")
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
        raise AssertionError("P4-S5 series identity mismatch")
    for seed_index, seed in enumerate(candidate_seeds):
        if int(seed["seed_index"]) != seed_index or \
                str(seed["seed_hex"]).lower() != \
                pilot["seed_hex_order"][seed_index].lower() or \
                int(seed["accepted"]) + int(seed["rejected"]) != \
                sample_count or \
                int(seed["global_attempted"]) != sample_count or \
                int(seed["neighbor_attempted"]) != 0 or \
                int(seed["shell_attempted"]) != 0 or \
                not (0 <= int(seed["global_self"]) <= sample_count) or \
                int(seed["trace_hash"]) == 0:
            raise AssertionError("P4-S5 seed/chain census mismatch")
    maximum_tau = max(float(item["tau_int"]) for item in candidate_series)
    maximum_ratio = max(float(item["predicted_se_budget_ratio"])
                        for item in candidate_entries)
    require_close(decision["maximum_tau_int"], maximum_tau,
                  "P4-S5 maximum tau")
    require_close(decision["maximum_predicted_se_budget_ratio"],
                  maximum_ratio, "P4-S5 maximum ratio")
    if site == 4:
        tolerance = float(
            pilot["exact_l4_balance_tolerance_multiplier_dbl_epsilon"]
        ) * 2.220446049250313e-16
        for field in ("db_residual", "stationary_residual",
                      "proposal_row_residual",
                      "proposal_symmetry_residual"):
            if float(header[field]) > tolerance:
                raise AssertionError("P4-S5 exact L4 {}".format(field))
    expected_eligible = (
        maximum_ratio <= float(pilot["maximum_predicted_se_budget_ratio"])
        and maximum_tau <= float(pilot["maximum_tau_int"])
        and bool(int(decision["table_pass"]))
        and bool(int(decision["balance_pass"]))
        and bool(int(decision["restart_pass"]))
    )
    if bool(int(decision["physical_case_eligible"])) != expected_eligible:
        raise AssertionError("P4-S5 physical eligibility mismatch")
    diagnostics.update({
        "candidate_index": candidate_index,
        "candidate": candidate_key(rho),
        "rho": rho,
        "physical_case_eligible": expected_eligible,
        "maximum_tau_int": maximum_tau,
        "maximum_predicted_se_budget_ratio": maximum_ratio,
        "proposal_policy_hash": int(header["proposal_policy_hash"]),
        "proposal_model_hash": int(header["proposal_model_hash"]),
        "guide_policy_hash": int(header["guide_policy_hash"]),
        "guide_shape_hash": int(header["guide_shape_hash"]),
        "table_hash": int(header["table_hash"]),
        "anchor_count": int(header["anchor_count"]),
        "seed_summaries": candidate_seeds,
    })
    return diagnostics


def validate_run(run, policy, allow_smoke_policy=False):
    header = run["run"]
    if int(header.get("schema", 0)) != SCHEMA_VERSION or \
            header.get("fixture") != EXPECTED_FIXTURE:
        raise AssertionError("P4-S5 run header mismatch")
    site = int(header["site_count"])
    qp = int(header["qp_total"])
    sector_size = fixed_sector_size(site)
    if [site, qp] not in policy["scope"]["physical_case_order"] or \
            int(header["sector_size"]) != sector_size or \
            int(header["anchor_count"]) != sector_size:
        raise AssertionError("P4-S5 physical case/anchor mismatch")
    if not allow_smoke_policy and \
            int(header["cache_requested_bytes"]) != \
            int(policy["scope"]["statistical_cache_bytes"]):
        raise AssertionError("P4-S5 cache policy mismatch")
    controls = validate_controls(run, policy)
    candidates = [
        validate_candidate(run, index, policy, allow_smoke_policy)
        for index in range(CANDIDATE_COUNT)
    ]
    if len(set(item["proposal_policy_hash"] for item in candidates)) != 1 or \
            len(set(item["proposal_model_hash"] for item in candidates)) != 1:
        raise AssertionError("P4-S5 proposal identity changed across rho")
    if len(set(item["guide_policy_hash"] for item in candidates)) != \
            CANDIDATE_COUNT or \
            len(set(item["guide_shape_hash"] for item in candidates)) != \
            CANDIDATE_COUNT or \
            len(set(item["table_hash"] for item in candidates)) != \
            CANDIDATE_COUNT:
        raise AssertionError("P4-S5 guide/restart/table identity collision")
    eta_over_rho = [float(run["candidates"][index]["eta"]) /
                    float(policy["guide"]["candidate_rho_order"][index])
                    for index in range(CANDIDATE_COUNT)]
    for value in eta_over_rho[1:]:
        require_close(value, eta_over_rho[0], "P4-S5 eta/rho")
    require_close(float(run["controls"][0]["eta"]) / 0.01,
                  eta_over_rho[0], "P4-S5 rho control eta/rho")
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
        raise AssertionError("P4-S5 evidence requires six unique cases")
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
    for candidate_index, rho in enumerate(
            policy["guide"]["candidate_rho_order"]):
        cases = [case["candidates"][candidate_index]
                 for case in ordered_cases]
        eligible = all(case["physical_case_eligible"] for case in cases)
        summary = {
            "candidate_index": candidate_index,
            "candidate": candidate_key(rho),
            "rho": float(rho),
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
            "maximum_importance_weight_maximum_mean_ratio": max(
                case["importance_weight_maximum_mean_ratio"]
                for case in cases),
            "minimum_importance_weight_ess_fraction": min(
                case["importance_weight_ess_fraction"] for case in cases),
            "cases": cases,
        }
        candidate_summaries.append(summary)
        if eligible:
            eligible_indices.append(candidate_index)
    generated_at = p4c.generated_at_jst()
    evidence = {
        "schema_version": SCHEMA_VERSION,
        "artifact": "p4s5_flat_guide_exact_table_pilot",
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
        "artifact": "p4s5_frozen_pilot_eligibility",
        "generated_at_jst": generated_at,
        "source_commit": args.source_commit,
        "measurement_policy_sha256": evidence["measurement_policy_sha256"],
        "candidate_order": [candidate_key(item) for item in
                            policy["guide"]["candidate_rho_order"]],
        "eligible_candidate_indices": eligible_indices,
        "eligible_candidates": [candidate_summaries[index]["candidate"]
                                for index in eligible_indices],
        "preferred_candidate_index": eligible_indices[0]
        if eligible_indices else None,
        "pilot_decision": "ELIGIBLE" if eligible_indices else "STOP",
        "failure_action": "freeze_list_then_run_timing_census"
        if eligible_indices else "do_not_run_timing_or_submit_remote_job",
    }
    write_json(os.path.join(output, "p4s5_pilot_evidence.json"), evidence)
    write_json(os.path.join(output, "p4s5_pilot_eligibility.json"),
               eligibility)
    count = p4c.write_checksums(output)
    print("P4-S5 pilot evidence complete: {} candidates eligible; "
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
