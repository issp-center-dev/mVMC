from __future__ import print_function

import argparse
import math
import os
import shutil

import generate_bounded_krylov_p4c_evidence as p4c


SCHEMA_VERSION = 1
EXPECTED_FIXTURE = "p4s3_exact_table"
FAMILIES = ("S", "K", "B")
ENTRY_COUNT = 6
COMPONENTS = ("real", "imag")


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
    return "p{}/{}_rho_{:.0e}".format(
        int(candidate["global_numerator"]),
        int(candidate["global_denominator"]),
        float(candidate["rho"]),
    )


def validate_policy(policy, allow_smoke_policy=False):
    if int(policy.get("schema_version", 0)) != SCHEMA_VERSION:
        raise AssertionError("P4-S3 policy schema mismatch")
    if policy.get("policy_id") != "power-lanczos-zero-support-p4-r0-v3":
        raise AssertionError("P4-S3 policy id mismatch")
    scope = policy.get("scope", {})
    if scope.get("site_counts") != [4, 6, 8] or \
            scope.get("qp_total") != [1, 4]:
        raise AssertionError("P4-S3 physical scope mismatch")
    expected_cases = [[4, 1], [4, 4], [6, 1], [6, 4], [8, 1], [8, 4]]
    if scope.get("physical_case_order") != expected_cases:
        raise AssertionError("P4-S3 physical case order mismatch")
    candidates = policy.get("proposal", {}).get("candidate_order", [])
    expected_candidates = [
        (1, 8, 1.0e-2), (1, 8, 1.0e-4), (1, 8, 1.0e-6),
        (1, 4, 1.0e-2), (1, 4, 1.0e-4), (1, 4, 1.0e-6),
        (1, 2, 1.0e-2), (1, 2, 1.0e-4), (1, 2, 1.0e-6),
        (1, 1, 1.0e-2), (1, 1, 1.0e-4), (1, 1, 1.0e-6),
    ]
    observed_candidates = [
        (int(item["global_numerator"]),
         int(item["global_denominator"]), float(item["rho"]))
        for item in candidates
    ]
    if observed_candidates != expected_candidates:
        raise AssertionError("P4-S3 candidate order mismatch")
    pilot = policy.get("exact_table_pilot", {})
    if not allow_smoke_policy and \
            int(pilot.get("sample_count_per_seed", 0)) != 131072:
        raise AssertionError("P4-S3 pilot sample count mismatch")
    if pilot.get("seed_hex_order") != [
            "0x5034533350494c31", "0x5034533350494c32",
            "0x5034533350494c33", "0x5034533350494c34"]:
        raise AssertionError("P4-S3 pilot seed order mismatch")
    if float(pilot.get("maximum_predicted_se_budget_ratio", 0.0)) != 0.90:
        raise AssertionError("P4-S3 predicted-SE threshold mismatch")
    if float(pilot.get("maximum_tau_int", 0.0)) != 12.0:
        raise AssertionError("P4-S3 pilot tau threshold mismatch")
    if int(pilot.get("predicted_official_sample_count", 0)) != 4096:
        raise AssertionError("P4-S3 predicted sample count mismatch")
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
    if run is None or len(candidates) != 12 or len(decisions) != 12:
        raise AssertionError("incomplete P4-S3 pilot log {}".format(path))
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
    header = run["candidates"][candidate_index]
    decision = run["decisions"][candidate_index]
    pilot = policy["exact_table_pilot"]
    site = int(run["run"]["site_count"])
    sample_count = int(run["run"]["sample_count_per_seed"])
    cache = int(run["run"]["cache_requested_bytes"])
    if int(header["candidate_index"]) != candidate_index or \
            int(decision["candidate_index"]) != candidate_index:
        raise AssertionError("pilot candidate schedule mismatch")
    if int(header["site_count"]) != site or \
            int(header["qp_total"]) != int(run["run"]["qp_total"]):
        raise AssertionError("pilot candidate physical case mismatch")
    if int(header["sample_count_per_seed"]) != sample_count or \
            int(header["cache_requested_bytes"]) != cache:
        raise AssertionError("pilot candidate run contract mismatch")
    if int(header["global_numerator"]) != \
            int(policy_candidate["global_numerator"]) or \
            int(header["global_denominator"]) != \
            int(policy_candidate["global_denominator"]):
        raise AssertionError("pilot proposal fraction mismatch")
    require_close(header["rho"], policy_candidate["rho"], "pilot rho")
    if int(header["anchor_count"]) != fixed_sector_size(site):
        raise AssertionError("pilot anchor census mismatch")
    if int(header["proposal_policy_hash"]) == 0 or \
            int(header["proposal_model_hash"]) == 0 or \
            int(header["table_hash"]) == 0:
        raise AssertionError("pilot hash contract missing")
    if not allow_smoke_policy and sample_count != \
            int(pilot["sample_count_per_seed"]):
        raise AssertionError("official pilot sample count mismatch")

    candidate_series = [item for item in run["series"]
                        if int(item["candidate_index"]) == candidate_index]
    candidate_entries = [item for item in run["entries"]
                         if int(item["candidate_index"]) == candidate_index]
    candidate_seeds = [item for item in run["seeds"]
                       if int(item["candidate_index"]) == candidate_index]
    if len(candidate_series) != 4 * len(FAMILIES) * ENTRY_COUNT * 2:
        raise AssertionError("pilot series census mismatch")
    if len(candidate_entries) != 4 * len(FAMILIES) * ENTRY_COUNT:
        raise AssertionError("pilot entry census mismatch")
    if len(candidate_seeds) != 4:
        raise AssertionError("pilot seed census mismatch")
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
        raise AssertionError("pilot series identity mismatch")
    expected_entries = {
        (seed_index, family, entry)
        for seed_index in range(4)
        for family in FAMILIES
        for entry in range(ENTRY_COUNT)
    }
    observed_entries = {
        (int(item["seed_index"]), item["family"], int(item["entry"]))
        for item in candidate_entries
    }
    if observed_entries != expected_entries:
        raise AssertionError("pilot entry identity mismatch")
    for seed_index, seed in enumerate(candidate_seeds):
        if int(seed["seed_index"]) != seed_index or \
                str(seed["seed_hex"]).lower() != \
                pilot["seed_hex_order"][seed_index].lower():
            raise AssertionError("pilot seed order mismatch")
        if int(seed["accepted"]) + int(seed["rejected"]) != sample_count:
            raise AssertionError("pilot accept/reject census mismatch")
        if int(seed["neighbor_attempted"]) + \
                int(seed["global_attempted"]) != sample_count:
            raise AssertionError("pilot component census mismatch")
        if int(seed["global_self"]) > int(seed["global_attempted"]):
            raise AssertionError("pilot global-self census mismatch")
        if int(seed["trace_hash"]) == 0:
            raise AssertionError("pilot trace hash missing")
    maximum_tau = max(float(item["tau_int"])
                      for item in candidate_series)
    maximum_ratio = max(float(item["predicted_se_budget_ratio"])
                        for item in candidate_entries)
    require_close(decision["maximum_tau_int"], maximum_tau,
                  "pilot maximum tau")
    require_close(decision["maximum_predicted_se_budget_ratio"],
                  maximum_ratio, "pilot maximum predicted ratio")
    balance_pass = bool(int(header["exact_balance_pass"])) and \
        bool(int(decision["balance_pass"]))
    restart_pass = bool(int(header["restart_pass"])) and \
        bool(int(decision["restart_pass"]))
    if site == 4:
        tolerance = float(
            pilot["exact_l4_balance_tolerance_multiplier_dbl_epsilon"]
        ) * 2.220446049250313e-16
        if float(header["db_residual"]) > tolerance or \
                float(header["stationary_residual"]) > tolerance:
            raise AssertionError("pilot exact L4 balance residual")
    expected_eligible = (
        maximum_ratio <= float(pilot["maximum_predicted_se_budget_ratio"])
        and maximum_tau <= float(pilot["maximum_tau_int"])
        and bool(int(decision["table_pass"]))
        and balance_pass and restart_pass
    )
    if bool(int(decision["physical_case_eligible"])) != expected_eligible:
        raise AssertionError("pilot physical eligibility mismatch")
    return {
        "candidate_index": candidate_index,
        "candidate": candidate_key(policy_candidate),
        "global_numerator": int(policy_candidate["global_numerator"]),
        "global_denominator": int(policy_candidate["global_denominator"]),
        "rho": float(policy_candidate["rho"]),
        "physical_case_eligible": expected_eligible,
        "maximum_tau_int": maximum_tau,
        "maximum_predicted_se_budget_ratio": maximum_ratio,
        "proposal_policy_hash": int(header["proposal_policy_hash"]),
        "proposal_model_hash": int(header["proposal_model_hash"]),
        "table_hash": int(header["table_hash"]),
        "anchor_count": int(header["anchor_count"]),
        "exact_balance_pass": balance_pass,
        "restart_pass": restart_pass,
        "seed_summaries": candidate_seeds,
    }


def validate_run(run, policy, allow_smoke_policy=False):
    header = run["run"]
    if int(header.get("schema", 0)) != SCHEMA_VERSION or \
            header.get("fixture") != EXPECTED_FIXTURE:
        raise AssertionError("pilot run header mismatch")
    site = int(header["site_count"])
    qp = int(header["qp_total"])
    case = [site, qp]
    if case not in policy["scope"]["physical_case_order"]:
        raise AssertionError("pilot physical case outside scope")
    if int(header["sector_size"]) != fixed_sector_size(site) or \
            int(header["anchor_count"]) != fixed_sector_size(site):
        raise AssertionError("pilot sector/anchor size mismatch")
    expected_cache = int(policy["scope"]["statistical_cache_bytes"])
    if not allow_smoke_policy and \
            int(header["cache_requested_bytes"]) != expected_cache:
        raise AssertionError("pilot cache policy mismatch")
    candidates = [
        validate_candidate(run, index, policy, allow_smoke_policy)
        for index in range(12)
    ]
    hashes_by_rho = {}
    for candidate in candidates:
        previous = hashes_by_rho.setdefault(candidate["rho"],
                                             candidate["table_hash"])
        if previous != candidate["table_hash"]:
            raise AssertionError("pilot table changed with mixture fraction")
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


def materialize_logs(runs, output):
    raw_root = os.path.join(output, "raw_pilot")
    if not os.path.isdir(raw_root):
        os.makedirs(raw_root)
    mapped = {}
    for run in runs:
        name = "L{}_QP{}-{}-{}.log".format(
            int(run["run"]["site_count"]), int(run["run"]["qp_total"]),
            run["sha256"][:12], os.path.basename(run["path"]))
        destination = os.path.join(raw_root, name)
        shutil.copyfile(run["path"], destination)
        if sha256_file(destination) != run["sha256"]:
            raise AssertionError("pilot raw-log checksum mismatch")
        mapped[run["path"]] = os.path.relpath(destination, output)
    return mapped


def write_checksums(output):
    return p4c.write_checksums(output)


def analyze_command(args):
    policy = validate_policy(read_json(args.policy), args.allow_smoke_policy)
    parsed = [parse_pilot_log(path) for path in args.pilot_log]
    observed_cases = [
        [int(run["run"]["site_count"]), int(run["run"]["qp_total"])]
        for run in parsed
    ]
    if sorted(observed_cases) != sorted(policy["scope"]["physical_case_order"]):
        raise AssertionError("pilot evidence must contain six unique cases")
    if len(set(tuple(case) for case in observed_cases)) != len(observed_cases):
        raise AssertionError("duplicate pilot physical case")
    validated = {
        (int(run["run"]["site_count"]), int(run["run"]["qp_total"])):
        validate_run(run, policy, args.allow_smoke_policy)
        for run in parsed
    }
    output = os.path.abspath(args.output)
    if not os.path.isdir(output):
        os.makedirs(output)
    raw_paths = materialize_logs(parsed, output)
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
            "global_numerator": int(policy_candidate["global_numerator"]),
            "global_denominator": int(policy_candidate["global_denominator"]),
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
    evidence = {
        "schema_version": SCHEMA_VERSION,
        "artifact": "p4s3_exact_table_pilot",
        "generated_at_jst": p4c.generated_at_jst(),
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
        "artifact": "p4s3_frozen_pilot_eligibility",
        "generated_at_jst": evidence["generated_at_jst"],
        "source_commit": args.source_commit,
        "measurement_policy_sha256": evidence["measurement_policy_sha256"],
        "candidate_order": [candidate_key(item) for item in
                            policy["proposal"]["candidate_order"]],
        "eligible_candidate_indices": eligible_indices,
        "eligible_candidates": [
            candidate_summaries[index]["candidate"]
            for index in eligible_indices
        ],
        "pilot_decision": "ELIGIBLE" if eligible_indices else "STOP",
        "failure_action": (
            "freeze_list_then_run_timing_census"
            if eligible_indices else "do_not_submit_remote_job"
        ),
    }
    write_json(os.path.join(output, "p4s3_pilot_evidence.json"), evidence)
    write_json(os.path.join(output, "p4s3_pilot_eligibility.json"),
               eligibility)
    count = write_checksums(output)
    print("P4-S3 pilot evidence complete: {} candidates eligible; "
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
