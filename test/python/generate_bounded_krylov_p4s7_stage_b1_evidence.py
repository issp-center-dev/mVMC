from __future__ import print_function

import argparse
import os

import generate_bounded_krylov_p4c_evidence as p4c
import generate_bounded_krylov_p4s3_pilot_evidence as p4s3
import generate_bounded_krylov_p4s6_pilot_evidence as p4s6
import generate_bounded_krylov_p4s7_prefilter_evidence as p4s7


SCHEMA_VERSION = 6
EXPECTED_FIXTURE = "p4s7_partial_guide_neighbor_stage_b1"
EXPECTED_CANDIDATES = ((2, 1, 16), (2, 1, 32),
                       (2, 10, 16), (2, 10, 32))
FAMILIES = p4s6.FAMILIES
ENTRY_COUNT = p4s6.ENTRY_COUNT
COMPONENTS = p4s6.COMPONENTS
read_json = p4c.read_json
write_json = p4c.write_json
sha256_file = p4c.sha256_file
require_close = p4s6.require_close
fixed_sector_size = p4c.fixed_sector_size


def parse_key_values(line):
    return {key: p4c.parse_number(value)
            for key, value in p4c.parse_key_values(line, 1).items()}


def candidate_key(partial_order, alpha, step_count):
    return "m{}_alpha{}_K{}".format(partial_order, alpha, step_count)


def validate_policy(policy, allow_smoke_policy=False):
    if int(policy.get("schema_version", 0)) != SCHEMA_VERSION or \
            policy.get("policy_id") != \
            "power-lanczos-zero-support-p4-r0-v7-stage-b1":
        raise AssertionError("P4-S7 B1 policy identity mismatch")
    expected_cases = [[4, 1], [4, 4], [6, 1], [6, 4], [8, 1], [8, 4]]
    scope = policy.get("scope", {})
    if scope.get("physical_case_order") != expected_cases or \
            scope.get("site_counts") != [4, 6, 8] or \
            scope.get("qp_total") != [1, 4]:
        raise AssertionError("P4-S7 B1 scope mismatch")
    observed = tuple((int(x.get("m", 0)), int(x.get("alpha", 0)),
                      int(x.get("K", 0)))
                     for x in policy.get("candidate_order", []))
    if observed != EXPECTED_CANDIDATES or \
            policy.get("proposal", {}).get("grid_extension") != \
            "forbidden_under_v7_stage_b1":
        raise AssertionError("P4-S7 B1 frontier mismatch")
    pilot = policy.get("exact_table_pilot", {})
    if not allow_smoke_policy and \
            int(pilot.get("sample_count_per_seed", 0)) != 131072:
        raise AssertionError("P4-S7 B1 sample mismatch")
    if pilot.get("seed_hex_order") != [
            "0x5034533850494c31", "0x5034533850494c32",
            "0x5034533850494c33", "0x5034533850494c34"] or \
            int(pilot.get("predicted_official_sample_count", 0)) != 4096 or \
            float(pilot.get("maximum_predicted_se_budget_ratio", 0)) != \
            0.90 or float(pilot.get("maximum_tau_int", 0)) != 12.0:
        raise AssertionError("P4-S7 B1 gate mismatch")
    return policy


def parse_log(path):
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
            run = parse_key_values(line)
        elif line.startswith("PILOT_PARTIAL_NEIGHBOR_CANDIDATE "):
            candidates.append(parse_key_values(line))
        elif line.startswith("PILOT_SERIES "):
            series.append(parse_key_values(line))
        elif line.startswith("PILOT_ENTRY "):
            entries.append(parse_key_values(line))
        elif line.startswith("PILOT_SEED "):
            seeds.append(parse_key_values(line))
        elif line.startswith("PILOT_DECISION "):
            decisions.append(parse_key_values(line))
    if run is None or len(candidates) != 4 or len(decisions) != 4:
        raise AssertionError("incomplete P4-S7 B1 log {}".format(path))
    return {"path": os.path.abspath(path), "sha256": sha256_file(path),
            "run": run, "candidates": candidates, "series": series,
            "entries": entries, "seeds": seeds, "decisions": decisions,
            "line_count": len(text.splitlines())}


def validate_candidate(run, index, policy, allow_smoke_policy=False):
    partial_order, alpha, step_count = EXPECTED_CANDIDATES[index]
    header = run["candidates"][index]
    decision = run["decisions"][index]
    pilot = policy["exact_table_pilot"]
    site = int(run["run"]["site_count"])
    qp = int(run["run"]["qp_total"])
    sector = int(run["run"]["sector_size"])
    samples = int(run["run"]["sample_count_per_seed"])
    if int(header["schema"]) != SCHEMA_VERSION or \
            int(header["candidate_index"]) != index or \
            int(decision["candidate_index"]) != index or \
            int(header["site_count"]) != site or \
            int(header["qp_total"]) != qp or \
            int(header["sample_count_per_seed"]) != samples or \
            int(header["partial_order"]) != partial_order or \
            int(header["floor_multiplier"]) != alpha or \
            int(header["step_count"]) != step_count or \
            header["inner_kernel"] != "neighbor_only":
        raise AssertionError("P4-S7 B1 candidate schedule mismatch")
    require_close(header["rho"], 0.01, "P4-S7 B1 rho")
    for field in ("proposal_policy_hash", "proposal_model_hash",
                  "guide_policy_hash", "guide_shape_hash",
                  "partial_policy_hash", "partial_neighbor_policy_hash",
                  "target_table_hash", "partial_table_hash"):
        if int(header[field]) == 0:
            raise AssertionError("P4-S7 B1 zero {}".format(field))
    if int(header["anchor_count"]) != sector or \
            not bool(int(header["exact_balance_pass"])) or \
            not bool(int(header["restart_pass"])):
        raise AssertionError("P4-S7 B1 balance/restart mismatch")
    anchor_tolerance = float(pilot[
        "anchor_tolerance_multiplier_dbl_epsilon"]) * 2.220446049250313e-16
    if float(header["maximum_anchor_relative_residual"]) > anchor_tolerance:
        raise AssertionError("P4-S7 B1 anchor mismatch")
    diagnostics = p4s7.validate_diagnostics(
        header, "P4-S7 B1 candidate {}".format(index))
    candidate_series = [x for x in run["series"]
                        if int(x["candidate_index"]) == index]
    candidate_entries = [x for x in run["entries"]
                         if int(x["candidate_index"]) == index]
    candidate_seeds = [x for x in run["seeds"]
                       if int(x["candidate_index"]) == index]
    if len(candidate_series) != 4 * len(FAMILIES) * ENTRY_COUNT * 2 or \
            len(candidate_entries) != 4 * len(FAMILIES) * ENTRY_COUNT or \
            len(candidate_seeds) != 4:
        raise AssertionError("P4-S7 B1 census mismatch")
    expected_series = {(seed, family, entry, component)
                       for seed in range(4) for family in FAMILIES
                       for entry in range(ENTRY_COUNT)
                       for component in COMPONENTS}
    if {(int(x["seed_index"]), x["family"], int(x["entry"]), x["component"])
            for x in candidate_series} != expected_series:
        raise AssertionError("P4-S7 B1 series mismatch")
    for seed_index, seed in enumerate(candidate_seeds):
        attempted = samples * step_count
        if int(seed["seed_index"]) != seed_index or \
                str(seed["seed_hex"]).lower() != \
                pilot["seed_hex_order"][seed_index].lower() or \
                int(seed["accepted"]) + int(seed["rejected"]) != samples or \
                int(seed["surrogate_inner_attempted"]) != attempted or \
                int(seed["surrogate_inner_accepted"]) + \
                int(seed["surrogate_inner_rejected"]) != attempted or \
                int(seed["surrogate_final_changed"]) + \
                int(seed["surrogate_final_self"]) != samples or \
                int(seed["trace_hash"]) == 0:
            raise AssertionError("P4-S7 B1 seed mismatch")
    maximum_tau = max(float(x["tau_int"]) for x in candidate_series)
    maximum_ratio = max(float(x["predicted_se_budget_ratio"])
                        for x in candidate_entries)
    require_close(decision["maximum_tau_int"], maximum_tau,
                  "P4-S7 B1 tau")
    require_close(decision["maximum_predicted_se_budget_ratio"],
                  maximum_ratio, "P4-S7 B1 ratio")
    if site == 4:
        eps = 2.220446049250313e-16
        tolerance = float(pilot[
            "exact_l4_balance_tolerance_multiplier_dbl_epsilon"]) * eps
        powered_tolerance = float(pilot[
            "exact_l4_powered_ratio_tolerance_multiplier_dbl_epsilon"]) * eps
        for field in ("inner_row_residual", "inner_db_residual",
                      "powered_row_residual", "powered_db_residual",
                      "outer_db_residual", "outer_stationary_residual"):
            if float(header[field]) > tolerance:
                raise AssertionError("P4-S7 B1 L4 {}".format(field))
        if float(header["powered_ratio_residual"]) > powered_tolerance:
            raise AssertionError("P4-S7 B1 powered ratio")
    eligible = maximum_ratio <= float(
        pilot["maximum_predicted_se_budget_ratio"]) and \
        maximum_tau <= float(pilot["maximum_tau_int"]) and \
        bool(int(decision["table_pass"])) and \
        bool(int(decision["balance_pass"])) and \
        bool(int(decision["restart_pass"]))
    if bool(int(decision["physical_case_eligible"])) != eligible:
        raise AssertionError("P4-S7 B1 eligibility mismatch")
    diagnostics.update({
        "candidate_index": index,
        "candidate": candidate_key(*EXPECTED_CANDIDATES[index]),
        "m": partial_order, "alpha": alpha, "K": step_count,
        "physical_case_eligible": eligible,
        "maximum_tau_int": maximum_tau,
        "maximum_predicted_se_budget_ratio": maximum_ratio,
        "proposal_policy_hash": int(header["proposal_policy_hash"]),
        "proposal_model_hash": int(header["proposal_model_hash"]),
        "guide_policy_hash": int(header["guide_policy_hash"]),
        "guide_shape_hash": int(header["guide_shape_hash"]),
        "partial_policy_hash": int(header["partial_policy_hash"]),
        "partial_neighbor_policy_hash": int(
            header["partial_neighbor_policy_hash"]),
        "target_table_hash": int(header["target_table_hash"]),
        "partial_table_hash": int(header["partial_table_hash"]),
        "seed_summaries": candidate_seeds,
    })
    return diagnostics


def validate_run(run, policy, allow_smoke_policy=False):
    header = run["run"]
    if int(header.get("schema", 0)) != SCHEMA_VERSION or \
            header.get("fixture") != EXPECTED_FIXTURE:
        raise AssertionError("P4-S7 B1 run mismatch")
    site = int(header["site_count"])
    qp = int(header["qp_total"])
    sector = fixed_sector_size(site)
    if [site, qp] not in policy["scope"]["physical_case_order"] or \
            int(header["sector_size"]) != sector or \
            int(header["anchor_count"]) != sector:
        raise AssertionError("P4-S7 B1 physical case mismatch")
    if not allow_smoke_policy and int(header["cache_requested_bytes"]) != \
            int(policy["scope"]["statistical_cache_bytes"]):
        raise AssertionError("P4-S7 B1 cache mismatch")
    candidates = [validate_candidate(run, i, policy, allow_smoke_policy)
                  for i in range(4)]
    for field in ("proposal_policy_hash", "proposal_model_hash",
                  "guide_policy_hash", "guide_shape_hash",
                  "target_table_hash"):
        if len(set(x[field] for x in candidates)) != 1:
            raise AssertionError("P4-S7 B1 shared identity {}".format(field))
    if len(set(x["partial_neighbor_policy_hash"] for x in candidates)) != 4:
        raise AssertionError("P4-S7 B1 neighbor hash collision")
    for alpha in (1, 10):
        subset = [x for x in candidates if x["alpha"] == alpha]
        if len(set(x["partial_policy_hash"] for x in subset)) != 1 or \
                len(set(x["partial_table_hash"] for x in subset)) != 1:
            raise AssertionError("P4-S7 B1 partial table mismatch")
    return {"site_count": site, "qp_total": qp,
            "rank_count": int(header["rank_count"]),
            "sample_count_per_seed": int(header["sample_count_per_seed"]),
            "cache_requested_bytes": int(header["cache_requested_bytes"]),
            "sector_size": sector, "plan_hash": int(header["plan_hash"]),
            "candidates": candidates}


def analyze_command(args):
    policy = validate_policy(read_json(args.policy), args.allow_smoke_policy)
    parsed = [parse_log(path) for path in args.pilot_log]
    observed = [[int(x["run"]["site_count"]), int(x["run"]["qp_total"])]
                for x in parsed]
    if sorted(observed) != sorted(policy["scope"]["physical_case_order"]) or \
            len(set(tuple(x) for x in observed)) != 6:
        raise AssertionError("P4-S7 B1 requires six cases")
    validated = {(int(x["run"]["site_count"]), int(x["run"]["qp_total"])):
                 validate_run(x, policy, args.allow_smoke_policy)
                 for x in parsed}
    output = os.path.abspath(args.output)
    if not os.path.isdir(output):
        os.makedirs(output)
    raw_paths = p4s3.materialize_logs(parsed, output)
    cases = []
    for site, qp in policy["scope"]["physical_case_order"]:
        item = next(x for x in parsed if int(x["run"]["site_count"]) == site
                    and int(x["run"]["qp_total"]) == qp)
        summary = validated[(site, qp)]
        summary["raw_log"] = raw_paths[item["path"]]
        summary["raw_log_sha256"] = item["sha256"]
        summary["raw_log_line_count"] = item["line_count"]
        cases.append(summary)
    summaries = []
    eligible_indices = []
    for index, candidate in enumerate(EXPECTED_CANDIDATES):
        candidate_cases = [case["candidates"][index] for case in cases]
        eligible = all(x["physical_case_eligible"] for x in candidate_cases)
        summary = {"candidate_index": index,
                   "candidate": candidate_key(*candidate),
                   "m": candidate[0], "alpha": candidate[1], "K": candidate[2],
                   "stage_b1_eligible": eligible,
                   "maximum_tau_int": max(x["maximum_tau_int"] for x in candidate_cases),
                   "maximum_predicted_se_budget_ratio": max(
                       x["maximum_predicted_se_budget_ratio"] for x in candidate_cases),
                   "cases": candidate_cases}
        summaries.append(summary)
        if eligible:
            eligible_indices.append(index)
    generated = p4c.generated_at_jst()
    policy_sha = sha256_file(args.policy)
    evidence = {"schema_version": SCHEMA_VERSION,
                "artifact": "p4s7_partial_neighbor_stage_b1",
                "generated_at_jst": generated,
                "source_commit": args.source_commit,
                "producer_binary": os.path.basename(args.producer_binary),
                "producer_binary_sha256": sha256_file(args.producer_binary),
                "measurement_policy_sha256": policy_sha,
                "measurement_policy": policy, "physical_cases": cases,
                "candidate_summaries": summaries}
    eligibility = {"schema_version": SCHEMA_VERSION,
                   "artifact": "p4s7_stage_b1_eligibility",
                   "generated_at_jst": generated,
                   "source_commit": args.source_commit,
                   "measurement_policy_sha256": policy_sha,
                   "candidate_order": [candidate_key(*x) for x in EXPECTED_CANDIDATES],
                   "eligible_candidate_indices": eligible_indices,
                   "eligible_candidates": [summaries[i]["candidate"] for i in eligible_indices],
                   "stage_b1_decision": "ELIGIBLE" if eligible_indices else "STOP",
                   "failure_action": "implement_real_callback_only_for_survivors"
                   if eligible_indices else "do_not_implement_real_callback_or_timing"}
    write_json(os.path.join(output, "p4s7_stage_b1_evidence.json"), evidence)
    write_json(os.path.join(output, "p4s7_stage_b1_eligibility.json"), eligibility)
    count = p4c.write_checksums(output)
    print("P4-S7 B1 complete: {} candidates eligible; {} ledger files".
          format(len(eligible_indices), count))


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
