from __future__ import print_function

import argparse
import json
import os
import shutil
import sys
import tempfile

import generate_bounded_krylov_p4s7_stage_b1_evidence as evidence


SEEDS = ["0x5034533850494c31", "0x5034533850494c32",
         "0x5034533850494c33", "0x5034533850494c34"]


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative P4-S7 B1 evidence accepted: {}".format(
        label))


def write_log(path, site, qp, bad_hash=False, bad_census=False,
              bad_balance=False, bad_decision=False, bad_anchor=False):
    sample_count = 256
    sector_size = evidence.fixed_sector_size(site)
    with open(path, "w") as out:
        out.write(
            "PILOT_RUN schema=6 "
            "fixture=p4s7_partial_guide_neighbor_stage_b1 site_count={} "
            "qp_total={} rank_count=1 sample_count_per_seed={} "
            "cache_requested_bytes=0 sector_size={} anchor_count={} "
            "plan_hash=4000\n".format(site, qp, sample_count, sector_size,
                                      sector_size))
        for index, (partial_order, alpha, step_count) in enumerate(
                evidence.EXPECTED_CANDIDATES):
            neighbor_hash = 7000 + index
            if bad_hash and index == 1:
                neighbor_hash = 7000
            anchor = 1.0 if bad_anchor and index == 0 else 0.0
            balance = 1.0 if bad_balance and index == 0 else 0.0
            partial_hash = 7100 + alpha
            partial_table_hash = 8100 + alpha
            out.write(
                "PILOT_PARTIAL_NEIGHBOR_CANDIDATE schema=6 "
                "candidate_index={} site_count={} qp_total={} "
                "sample_count_per_seed={} cache_requested_bytes=0 "
                "rho=0.01 eta=0.005 partial_order={} "
                "floor_multiplier={} step_count={} "
                "inner_kernel=neighbor_only proposal_policy_hash=5000 "
                "proposal_model_hash=5001 guide_policy_hash=5002 "
                "guide_shape_hash=5003 partial_policy_hash={} "
                "partial_neighbor_policy_hash={} target_table_hash=6000 "
                "partial_table_hash={} anchor_count={} "
                "maximum_anchor_relative_residual={} partial_minimum=1 "
                "partial_maximum=2 partial_mean=1.5 "
                "target_ratio_minimum=0.5 target_ratio_maximum=2 "
                "target_ratio_maximum_mean_ratio=1.2 "
                "target_ratio_ess_fraction=0.8 ideal_acceptance=0.8 "
                "ideal_configuration_change=0.7 "
                "maximum_iid_se_budget_ratio=0.45 "
                "required_tau_for_0p90=4 iid_limiting_family=S "
                "iid_limiting_entry=0 exact_balance_pass=1 "
                "inner_row_residual=0 inner_db_residual={} "
                "powered_row_residual=0 powered_db_residual=0 "
                "powered_ratio_residual=0 outer_db_residual=0 "
                "outer_stationary_residual=0 restart_pass=1\n".format(
                    index, site, qp, sample_count, partial_order, alpha,
                    step_count, partial_hash, neighbor_hash,
                    partial_table_hash, sector_size, anchor, balance))
            for seed_index, seed_hex in enumerate(SEEDS):
                for family in evidence.FAMILIES:
                    for entry in range(evidence.ENTRY_COUNT):
                        for component in evidence.COMPONENTS:
                            out.write(
                                "PILOT_SERIES candidate_index={} "
                                "seed_index={} family={} entry={} "
                                "component={} iid_variance=1 tau_int=1\n".
                                format(index, seed_index, family, entry,
                                       component))
                        out.write(
                            "PILOT_ENTRY candidate_index={} seed_index={} "
                            "family={} entry={} predicted_se=0.01 "
                            "budget=0.02 predicted_se_budget_ratio=0.5\n".
                            format(index, seed_index, family, entry))
                attempted = sample_count * step_count
                if bad_census and index == 0 and seed_index == 0:
                    attempted -= 1
                out.write(
                    "PILOT_SEED candidate_index={} seed_index={} seed={} "
                    "seed_hex={} maximum_tau_int=1 "
                    "maximum_predicted_se_budget_ratio=0.5 accepted=200 "
                    "rejected=56 surrogate_inner_attempted={} "
                    "surrogate_inner_accepted={} surrogate_inner_rejected=0 "
                    "surrogate_final_changed=200 surrogate_final_self=56 "
                    "trace_hash={}\n".format(
                        index, seed_index, int(seed_hex, 16), seed_hex,
                        attempted, attempted,
                        9000 + index * 4 + seed_index))
            eligible = 0 if bad_decision and index == 0 else 1
            out.write(
                "PILOT_DECISION candidate_index={} physical_case_eligible={} "
                "maximum_predicted_se_budget_ratio=0.5 maximum_tau_int=1 "
                "table_pass=1 balance_pass=1 restart_pass=1 "
                "resource_limits_hash=10000\n".format(index, eligible))


def args(policy, output, logs):
    return argparse.Namespace(
        policy=policy, output=output, pilot_log=logs,
        source_commit="39d48c272998521d78166d38acf60c6b2e8624e5",
        producer_binary=sys.executable, allow_smoke_policy=True)


def main():
    if len(sys.argv) != 2:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4s7_stage_b1_evidence.py "
            "POLICY")
    policy_path = os.path.abspath(sys.argv[1])
    policy = evidence.validate_policy(evidence.read_json(policy_path), True)
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4s7-b1-")
    try:
        logs = []
        for site, qp in policy["scope"]["physical_case_order"]:
            path = os.path.join(root, "L{}-QP{}.log".format(site, qp))
            write_log(path, site, qp)
            logs.append(path)
        output = os.path.join(root, "evidence")
        evidence.analyze_command(args(policy_path, output, logs))
        with open(os.path.join(output, "p4s7_stage_b1_eligibility.json"),
                  "r") as handle:
            result = json.load(handle)
        if result["eligible_candidate_indices"] != list(range(4)) or \
                result["stage_b1_decision"] != "ELIGIBLE":
            raise AssertionError("P4-S7 B1 synthetic result mismatch")
        negatives = (("hash", {"bad_hash": True}),
                     ("census", {"bad_census": True}),
                     ("balance", {"bad_balance": True}),
                     ("decision", {"bad_decision": True}),
                     ("anchor", {"bad_anchor": True}))
        for label, options in negatives:
            bad = os.path.join(root, "bad-{}.log".format(label))
            write_log(bad, 4, 1, **options)
            expect_failure(
                lambda bad=bad, label=label: evidence.analyze_command(
                    args(policy_path,
                         os.path.join(root,
                                      "bad-{}-evidence".format(label)),
                         [bad] + logs[1:])), label)
        print("bounded Krylov P4-S7 B1 evidence pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
