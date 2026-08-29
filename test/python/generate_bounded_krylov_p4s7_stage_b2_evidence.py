from __future__ import print_function

import argparse
import os

import generate_bounded_krylov_p4c_evidence as p4c
import generate_bounded_krylov_p4s3_pilot_evidence as p4s3


SCHEMA_VERSION = 7
EXPECTED_FIXTURE = "p4s7_partial_guide_real_callback_stage_b2"
EXPECTED_CANDIDATES = ((2, 1, 32), (2, 10, 32))
read_json = p4c.read_json
write_json = p4c.write_json
sha256_file = p4c.sha256_file
fixed_sector_size = p4c.fixed_sector_size


def parse_key_values(line):
    return {key: p4c.parse_number(value)
            for key, value in p4c.parse_key_values(line, 1).items()}


def candidate_key(candidate):
    return "m{}_alpha{}_K{}".format(*candidate)


def validate_policy(policy, allow_smoke_policy=False):
    if int(policy.get("schema_version", 0)) != SCHEMA_VERSION or \
            policy.get("policy_id") != \
            "power-lanczos-zero-support-p4-r0-v7-stage-b2":
        raise AssertionError("P4-S7 B2 policy identity mismatch")
    expected_cases = [[4, 1], [4, 4], [6, 1], [6, 4], [8, 1], [8, 4]]
    scope = policy.get("scope", {})
    if scope.get("physical_case_order") != expected_cases or \
            scope.get("site_counts") != [4, 6, 8] or \
            scope.get("qp_total") != [1, 4]:
        raise AssertionError("P4-S7 B2 scope mismatch")
    observed = tuple((int(x.get("m", 0)), int(x.get("alpha", 0)),
                      int(x.get("K", 0)))
                     for x in policy.get("candidate_order", []))
    if observed != EXPECTED_CANDIDATES:
        raise AssertionError("P4-S7 B2 frontier mismatch")
    callback = policy.get("callback", {})
    if int(callback.get("bounded_max_order", 0)) != 2 or \
            int(callback.get("sample_count_per_seed", 0)) != 256 or \
            callback.get("seed_hex_order") != [
                "0x5034533850494c31", "0x5034533850494c32",
                "0x5034533850494c33", "0x5034533850494c34"] or \
            policy.get("statistics") != \
            "frozen_stage_b1_only_not_reselected_in_stage_b2":
        raise AssertionError("P4-S7 B2 callback contract mismatch")
    if not allow_smoke_policy and int(
            scope.get("statistical_cache_bytes", 0)) != 268435456:
        raise AssertionError("P4-S7 B2 cache mismatch")
    return policy


def parse_log(path):
    run = None
    candidates = []
    seeds = []
    decisions = []
    with open(path, "r") as handle:
        text = handle.read()
    for line in text.splitlines():
        if line.startswith("PILOT_RUN "):
            run = parse_key_values(line)
        elif line.startswith("PILOT_PARTIAL_CALLBACK_CANDIDATE "):
            candidates.append(parse_key_values(line))
        elif line.startswith("PILOT_CALLBACK_SEED "):
            seeds.append(parse_key_values(line))
        elif line.startswith("PILOT_CALLBACK_DECISION "):
            decisions.append(parse_key_values(line))
    if run is None or len(candidates) != 2 or len(seeds) != 8 or \
            len(decisions) != 2:
        raise AssertionError("incomplete P4-S7 B2 log {}".format(path))
    return {"path": os.path.abspath(path), "sha256": sha256_file(path),
            "run": run, "candidates": candidates, "seeds": seeds,
            "decisions": decisions, "line_count": len(text.splitlines())}


def validate_candidate(run, index, policy):
    expected_m, expected_alpha, expected_k = EXPECTED_CANDIDATES[index]
    header = run["candidates"][index]
    decision = run["decisions"][index]
    callback = policy["callback"]
    site = int(run["run"]["site_count"])
    qp = int(run["run"]["qp_total"])
    samples = int(run["run"]["sample_count_per_seed"])
    sector = int(run["run"]["sector_size"])
    if int(header["schema"]) != SCHEMA_VERSION or \
            int(header["candidate_index"]) != index or \
            int(decision["candidate_index"]) != index or \
            int(header["site_count"]) != site or \
            int(header["qp_total"]) != qp or \
            int(header["sample_count_per_seed"]) != samples or \
            int(header["partial_order"]) != expected_m or \
            int(header["floor_multiplier"]) != expected_alpha or \
            int(header["step_count"]) != expected_k or \
            header["callback"] != "max_order_2_bounded" or \
            header["inner_kernel"] != "neighbor_only":
        raise AssertionError("P4-S7 B2 candidate schedule mismatch")
    identity_fields = (
        "proposal_policy_hash", "proposal_model_hash", "guide_policy_hash",
        "guide_shape_hash", "partial_policy_hash", "callback_policy_hash",
        "target_table_hash", "partial_table_hash", "callback_table_hash",
        "partial_plan_hash")
    if any(int(header[field]) == 0 for field in identity_fields) or \
            int(header["anchor_count"]) != sector:
        raise AssertionError("P4-S7 B2 identity mismatch")
    tolerance = float(callback[
        "all_state_anchor_tolerance_multiplier_dbl_epsilon"
    ]) * 2.220446049250313e-16
    if float(header["maximum_anchor_relative_residual"]) > tolerance or \
            float(header["maximum_callback_log_residual"]) > tolerance or \
            float(header[
                "maximum_callback_weight_relative_residual"]) > tolerance or \
            not bool(int(header["anchor_pass"])) or \
            not bool(int(header["exact_balance_pass"])) or \
            not bool(int(header["restart_pass"])):
        raise AssertionError("P4-S7 B2 anchor/restart mismatch")
    if site == 4:
        eps = 2.220446049250313e-16
        balance_tolerance = float(callback[
            "exact_l4_balance_tolerance_multiplier_dbl_epsilon"]) * eps
        ratio_tolerance = float(callback[
            "exact_l4_powered_ratio_tolerance_multiplier_dbl_epsilon"]) * eps
        for field in ("inner_row_residual", "inner_db_residual",
                      "powered_row_residual", "powered_db_residual",
                      "outer_db_residual", "outer_stationary_residual"):
            if float(header[field]) > balance_tolerance:
                raise AssertionError("P4-S7 B2 L4 {}".format(field))
        if float(header["powered_ratio_residual"]) > ratio_tolerance:
            raise AssertionError("P4-S7 B2 powered ratio")
    candidate_seeds = [x for x in run["seeds"]
                       if int(x["candidate_index"]) == index]
    if len(candidate_seeds) != 4:
        raise AssertionError("P4-S7 B2 seed census")
    for seed_index, seed in enumerate(candidate_seeds):
        attempted = samples * expected_k
        if int(seed["seed_index"]) != seed_index or \
                str(seed["seed_hex"]).lower() != \
                callback["seed_hex_order"][seed_index].lower() or \
                not bool(int(seed["replay_pass"])) or \
                int(seed["accepted"]) + int(seed["rejected"]) != samples or \
                int(seed["surrogate_inner_attempted"]) != attempted or \
                int(seed["surrogate_inner_accepted"]) + \
                int(seed["surrogate_inner_rejected"]) != attempted or \
                int(seed["surrogate_final_changed"]) + \
                int(seed["surrogate_final_self"]) != samples or \
                int(seed["callback_evaluations"]) != attempted + 1 or \
                int(seed["callback_terminal_amplitude_calls"]) <= 0 or \
                int(seed["trace_hash"]) == 0:
            raise AssertionError("P4-S7 B2 seed mismatch")
    eligible = all(bool(int(decision[field])) for field in
                   ("table_pass", "balance_pass", "anchor_pass",
                    "restart_pass", "replay_pass"))
    if bool(int(decision["callback_eligible"])) != eligible:
        raise AssertionError("P4-S7 B2 decision mismatch")
    return {"candidate_index": index,
            "candidate": candidate_key(EXPECTED_CANDIDATES[index]),
            "m": expected_m, "alpha": expected_alpha, "K": expected_k,
            "callback_eligible": eligible,
            "callback_policy_hash": int(header["callback_policy_hash"]),
            "partial_plan_hash": int(header["partial_plan_hash"]),
            "maximum_callback_log_residual": float(
                header["maximum_callback_log_residual"]),
            "maximum_callback_weight_relative_residual": float(
                header["maximum_callback_weight_relative_residual"]),
            "seed_summaries": candidate_seeds}


def validate_run(run, policy, allow_smoke_policy=False):
    header = run["run"]
    if int(header.get("schema", 0)) != SCHEMA_VERSION or \
            header.get("fixture") != EXPECTED_FIXTURE:
        raise AssertionError("P4-S7 B2 run mismatch")
    site = int(header["site_count"])
    qp = int(header["qp_total"])
    sector = fixed_sector_size(site)
    if [site, qp] not in policy["scope"]["physical_case_order"] or \
            int(header["sector_size"]) != sector or \
            int(header["anchor_count"]) != sector or \
            int(header["sample_count_per_seed"]) != 256:
        raise AssertionError("P4-S7 B2 physical case mismatch")
    if not allow_smoke_policy and int(header["cache_requested_bytes"]) != \
            int(policy["scope"]["statistical_cache_bytes"]):
        raise AssertionError("P4-S7 B2 cache mismatch")
    candidates = [validate_candidate(run, i, policy) for i in range(2)]
    for field in ("partial_plan_hash",):
        if len(set(x[field] for x in candidates)) != 1:
            raise AssertionError("P4-S7 B2 shared {}".format(field))
    if len(set(x["callback_policy_hash"] for x in candidates)) != 2:
        raise AssertionError("P4-S7 B2 callback hash collision")
    return {"site_count": site, "qp_total": qp,
            "rank_count": int(header["rank_count"]),
            "sample_count_per_seed": 256,
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
        raise AssertionError("P4-S7 B2 requires six cases")
    validated = {(int(x["run"]["site_count"]),
                  int(x["run"]["qp_total"])):
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
        eligible = all(x["callback_eligible"] for x in candidate_cases)
        summaries.append({
            "candidate_index": index, "candidate": candidate_key(candidate),
            "m": candidate[0], "alpha": candidate[1], "K": candidate[2],
            "stage_b2_eligible": eligible,
            "maximum_callback_log_residual": max(
                x["maximum_callback_log_residual"] for x in candidate_cases),
            "maximum_callback_weight_relative_residual": max(
                x["maximum_callback_weight_relative_residual"]
                for x in candidate_cases),
            "cases": candidate_cases})
        if eligible:
            eligible_indices.append(index)
    generated = p4c.generated_at_jst()
    policy_sha = sha256_file(args.policy)
    evidence = {"schema_version": SCHEMA_VERSION,
                "artifact": "p4s7_partial_callback_stage_b2",
                "generated_at_jst": generated,
                "source_commit": args.source_commit,
                "producer_binary": os.path.basename(args.producer_binary),
                "producer_binary_sha256": sha256_file(args.producer_binary),
                "measurement_policy_sha256": policy_sha,
                "measurement_policy": policy, "physical_cases": cases,
                "candidate_summaries": summaries}
    eligibility = {
        "schema_version": SCHEMA_VERSION,
        "artifact": "p4s7_stage_b2_eligibility",
        "generated_at_jst": generated, "source_commit": args.source_commit,
        "measurement_policy_sha256": policy_sha,
        "candidate_order": [candidate_key(x) for x in EXPECTED_CANDIDATES],
        "eligible_candidate_indices": eligible_indices,
        "eligible_candidates": [summaries[i]["candidate"]
                                for i in eligible_indices],
        "stage_b2_decision": "ELIGIBLE" if eligible_indices else "STOP",
        "next_action": "run_real_callback_timing_for_survivors"
                       if eligible_indices else "do_not_run_timing"}
    write_json(os.path.join(output, "p4s7_stage_b2_evidence.json"), evidence)
    write_json(os.path.join(output, "p4s7_stage_b2_eligibility.json"),
               eligibility)
    count = p4c.write_checksums(output)
    print("P4-S7 B2 complete: {} candidates eligible; {} ledger files".
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
