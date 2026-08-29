from __future__ import print_function

import argparse
import json
import os
import shutil
import sys
import tempfile

import generate_bounded_krylov_p4s9_long_direct_evidence as evidence


SEEDS = ["0x5034533950494c31", "0x5034533950494c32",
         "0x5034533950494c33", "0x5034533950494c34"]


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative P4-S9 evidence accepted: {}".format(label))


def write_log(path, site, qp, bad_census=False, bad_hash=False,
              bad_decision=False, bad_threshold=False):
    sample_count = 256
    sector_size = evidence.fixed_sector_size(site)
    with open(path, "w") as out:
        out.write(
            "PILOT_RUN schema=8 fixture=p4s9_long_direct_exact_table "
            "site_count={} qp_total={} rank_count=1 "
            "sample_count_per_seed={} cache_requested_bytes=0 "
            "sector_size={} anchor_count={} plan_hash=100\n".format(
                site, qp, sample_count, sector_size, sector_size))
        threshold = 0.79 if bad_threshold else 0.80
        proposal_hash = 0 if bad_hash else 200
        out.write(
            "PILOT_CANDIDATE schema=8 candidate_index=0 site_count={} "
            "qp_total={} sample_count_per_seed={} cache_requested_bytes=0 "
            "rho=0.01 eta=0.1 predicted_official_sample_count=32768 "
            "maximum_predicted_ratio_threshold={} maximum_tau_threshold=12 "
            "proposal_kernel=1 global_numerator=0 global_denominator=1 "
            "neighbor_numerator=1 neighbor_denominator=1 "
            "distance_numerator=0 distance_denominator=1 "
            "distance_rounding_rule=0 shell_max_distance=0 "
            "shell_distance=0 shell_count=0 proposal_policy_hash={} "
            "proposal_model_hash=201 guide_policy_hash=202 "
            "guide_shape_hash=203 table_hash=204 anchor_count={} "
            "guide_minimum=1 guide_maximum=2 guide_mean=1.5 "
            "guide_maximum_minimum_ratio=2 "
            "uniform_sector_self_probability=0.1 "
            "importance_weight_maximum_mean_ratio=1.2 "
            "importance_weight_ess_fraction=0.8 "
            "maximum_iid_se_budget_ratio=0.4 required_tau_for_0p90=5.0625 "
            "iid_limiting_family=S iid_limiting_entry=0 "
            "exact_balance_pass=1 db_residual=0 stationary_residual=0 "
            "proposal_row_residual=0 proposal_symmetry_residual=0 "
            "shell_cardinality_pass=1 restart_pass=1\n".format(
                site, qp, sample_count, threshold, proposal_hash, sector_size))
        for seed_index, seed_hex in enumerate(SEEDS):
            for family in evidence.FAMILIES:
                for entry in range(evidence.ENTRY_COUNT):
                    for component in evidence.COMPONENTS:
                        out.write(
                            "PILOT_SERIES candidate_index=0 seed_index={} "
                            "family={} entry={} component={} iid_variance=1 "
                            "tau_int=1\n".format(seed_index, family, entry,
                                                 component))
                    out.write(
                        "PILOT_ENTRY candidate_index=0 seed_index={} "
                        "family={} entry={} predicted_se=0.01 budget=0.02 "
                        "predicted_se_budget_ratio=0.5\n".format(
                            seed_index, family, entry))
            attempted = sample_count - (1 if bad_census and seed_index == 0
                                        else 0)
            out.write(
                "PILOT_SEED candidate_index=0 seed_index={} seed={} "
                "seed_hex={} maximum_tau_int=1 "
                "maximum_predicted_se_budget_ratio=0.5 accepted=128 "
                "rejected=128 neighbor_attempted={} global_attempted=0 "
                "global_self=0 shell_attempted=0 trace_hash={}\n".format(
                    seed_index, int(seed_hex, 16), seed_hex, attempted,
                    300 + seed_index))
        eligible = 0 if bad_decision else 1
        out.write(
            "PILOT_DECISION candidate_index=0 physical_case_eligible={} "
            "maximum_predicted_se_budget_ratio=0.5 maximum_tau_int=1 "
            "table_pass=1 balance_pass=1 restart_pass=1 "
            "resource_limits_hash=400\n".format(eligible))


def args(policy, output, logs):
    return argparse.Namespace(
        policy=policy, output=output, pilot_log=logs,
        source_commit="39d48c272998521d78166d38acf60c6b2e8624e5",
        producer_binary=sys.executable, allow_smoke_policy=True)


def main():
    if len(sys.argv) != 2:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4s9_long_direct_evidence.py "
            "POLICY")
    policy_path = os.path.abspath(sys.argv[1])
    policy = evidence.validate_policy(evidence.read_json(policy_path), True)
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4s9-")
    try:
        logs = []
        for site, qp in policy["scope"]["physical_case_order"]:
            path = os.path.join(root, "L{}-QP{}.log".format(site, qp))
            write_log(path, site, qp)
            logs.append(path)
        output = os.path.join(root, "evidence")
        evidence.analyze_command(args(policy_path, output, logs))
        with open(os.path.join(output, "p4s9_stage_a_eligibility.json"),
                  "r") as handle:
            result = json.load(handle)
        if result["stage_a_decision"] != "ELIGIBLE" or \
                not result["stage_b_authorized"]:
            raise AssertionError("P4-S9 synthetic eligibility mismatch")
        negatives = (("census", {"bad_census": True}),
                     ("hash", {"bad_hash": True}),
                     ("decision", {"bad_decision": True}),
                     ("threshold", {"bad_threshold": True}))
        for label, options in negatives:
            bad = os.path.join(root, "bad-{}.log".format(label))
            write_log(bad, 4, 1, **options)
            expect_failure(
                lambda bad=bad, label=label: evidence.analyze_command(
                    args(policy_path, os.path.join(
                        root, "bad-{}-evidence".format(label)),
                        [bad] + logs[1:])), label)
        print("bounded Krylov P4-S9 Stage A evidence pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
