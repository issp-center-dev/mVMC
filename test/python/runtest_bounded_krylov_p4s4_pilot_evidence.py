from __future__ import print_function

import argparse
import json
import os
import shutil
import sys
import tempfile

import generate_bounded_krylov_p4s4_pilot_evidence as evidence


SEEDS = [
    "0x5034533450494c31", "0x5034533450494c32",
    "0x5034533450494c33", "0x5034533450494c34",
]


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative P4-S4 evidence accepted: {}".format(label))


def write_synthetic_log(path, site, qp, policy, bad_decision=False,
                        bad_distance=False, bad_table_hash=False):
    sample_count = 256
    sector_size = evidence.fixed_sector_size(site)
    proposal = policy["proposal"]
    with open(path, "w") as handle:
        handle.write(
            "PILOT_RUN schema=2 fixture=p4s4_fixed_distance_exact_table "
            "site_count={} qp_total={} rank_count=1 "
            "sample_count_per_seed={} cache_requested_bytes=0 "
            "sector_size={} anchor_count={} plan_hash=100\n".format(
                site, qp, sample_count, sector_size, sector_size))
        for index, candidate in enumerate(proposal["candidate_order"]):
            distance = evidence.resolved_distance(
                site, candidate["distance_numerator"],
                candidate["distance_denominator"])
            if bad_distance and index == 1:
                distance += 1
            table_hash = 1000 + (99 if bad_table_hash and index == 2 else 0)
            handle.write(
                "PILOT_CANDIDATE schema=2 candidate_index={} "
                "site_count={} qp_total={} sample_count_per_seed={} "
                "cache_requested_bytes=0 rho={} eta=0.1 "
                "proposal_kernel=123 global_numerator=0 "
                "global_denominator=0 neighbor_numerator=1 "
                "neighbor_denominator=8 distance_numerator={} "
                "distance_denominator={} distance_rounding_rule=1 "
                "shell_max_distance={} shell_distance={} shell_count=10 "
                "proposal_policy_hash={} proposal_model_hash={} "
                "table_hash={} anchor_count={} exact_balance_pass=1 "
                "db_residual=0 stationary_residual=0 "
                "proposal_row_residual=0 proposal_symmetry_residual=0 "
                "shell_cardinality_pass=1 restart_pass=1\n".format(
                    index, site, qp, sample_count, candidate["rho"],
                    candidate["distance_numerator"],
                    candidate["distance_denominator"], site, distance,
                    2000 + index, 3000 + site + index, table_hash,
                    sector_size))
            for seed_index, seed_hex in enumerate(SEEDS):
                for family in evidence.FAMILIES:
                    for entry in range(evidence.ENTRY_COUNT):
                        for component in evidence.COMPONENTS:
                            handle.write(
                                "PILOT_SERIES candidate_index={} "
                                "seed_index={} family={} entry={} "
                                "component={} iid_variance=1 tau_int=1\n".
                                format(index, seed_index, family, entry,
                                       component))
                        handle.write(
                            "PILOT_ENTRY candidate_index={} seed_index={} "
                            "family={} entry={} predicted_se=0.01 "
                            "budget=0.02 predicted_se_budget_ratio=0.5\n".
                            format(index, seed_index, family, entry))
                handle.write(
                    "PILOT_SEED candidate_index={} seed_index={} seed={} "
                    "seed_hex={} maximum_tau_int=1 "
                    "maximum_predicted_se_budget_ratio=0.5 accepted=128 "
                    "rejected=128 neighbor_attempted=32 global_attempted=0 "
                    "global_self=0 shell_attempted=224 trace_hash={}\n".
                    format(index, seed_index, int(seed_hex, 16), seed_hex,
                           4000 + index * 4 + seed_index))
            eligible = 0 if bad_decision and index == 0 else 1
            handle.write(
                "PILOT_DECISION candidate_index={} physical_case_eligible={} "
                "maximum_predicted_se_budget_ratio=0.5 maximum_tau_int=1 "
                "table_pass=1 balance_pass=1 restart_pass=1 "
                "resource_limits_hash=5000\n".format(index, eligible))


def analyze_args(policy_path, output, logs):
    return argparse.Namespace(
        policy=policy_path,
        output=output,
        pilot_log=logs,
        source_commit="39d48c272998521d78166d38acf60c6b2e8624e5",
        producer_binary=sys.executable,
        allow_smoke_policy=True,
    )


def main():
    if len(sys.argv) != 2:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4s4_pilot_evidence.py POLICY")
    policy_path = os.path.abspath(sys.argv[1])
    policy = evidence.validate_policy(evidence.read_json(policy_path), True)
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4s4-pilot-")
    try:
        logs = []
        for site, qp in policy["scope"]["physical_case_order"]:
            path = os.path.join(root, "L{}-QP{}.log".format(site, qp))
            write_synthetic_log(path, site, qp, policy)
            logs.append(path)
        output = os.path.join(root, "evidence")
        evidence.analyze_command(analyze_args(policy_path, output, logs))
        with open(os.path.join(output, "p4s4_pilot_eligibility.json"),
                  "r") as handle:
            eligibility = json.load(handle)
        if eligibility["eligible_candidate_indices"] != [0, 1, 2] or \
                eligibility["preferred_candidate_index"] != 0 or \
                eligibility["pilot_decision"] != "ELIGIBLE":
            raise AssertionError("synthetic P4-S4 eligibility mismatch")

        bad_decision = os.path.join(root, "bad-decision.log")
        write_synthetic_log(bad_decision, 4, 1, policy,
                            bad_decision=True)
        expect_failure(
            lambda: evidence.analyze_command(analyze_args(
                policy_path, os.path.join(root, "bad-decision-evidence"),
                [bad_decision] + logs[1:])),
            "decision AND mismatch")

        bad_distance = os.path.join(root, "bad-distance.log")
        write_synthetic_log(bad_distance, 4, 1, policy,
                            bad_distance=True)
        expect_failure(
            lambda: evidence.analyze_command(analyze_args(
                policy_path, os.path.join(root, "bad-distance-evidence"),
                [bad_distance] + logs[1:])),
            "resolved distance mismatch")

        bad_table = os.path.join(root, "bad-table.log")
        write_synthetic_log(bad_table, 4, 1, policy,
                            bad_table_hash=True)
        expect_failure(
            lambda: evidence.analyze_command(analyze_args(
                policy_path, os.path.join(root, "bad-table-evidence"),
                [bad_table] + logs[1:])),
            "proposal-dependent exact table")
        print("bounded Krylov P4-S4 pilot evidence pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
