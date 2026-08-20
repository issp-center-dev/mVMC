from __future__ import print_function

import argparse
import json
import os
import shutil
import sys
import tempfile

import generate_bounded_krylov_p4s7_stage_b2_evidence as evidence


SEEDS = ["0x5034533850494c31", "0x5034533850494c32",
         "0x5034533850494c33", "0x5034533850494c34"]


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative P4-S7 B2 evidence accepted: {}".format(
        label))


def write_log(path, site, qp, bad_hash=False, bad_census=False,
              bad_anchor=False, bad_replay=False, bad_decision=False):
    samples = 256
    sector = evidence.fixed_sector_size(site)
    with open(path, "w") as out:
        out.write(
            "PILOT_RUN schema=7 "
            "fixture=p4s7_partial_guide_real_callback_stage_b2 "
            "site_count={} qp_total={} rank_count=1 "
            "sample_count_per_seed=256 cache_requested_bytes=0 "
            "sector_size={} anchor_count={} plan_hash=4000\n".format(
                site, qp, sector, sector))
        for index, (partial_order, alpha, step_count) in enumerate(
                evidence.EXPECTED_CANDIDATES):
            callback_hash = 7000 + index
            if bad_hash and index == 1:
                callback_hash = 7000
            residual = 1.0 if bad_anchor and index == 0 else 0.0
            out.write(
                "PILOT_PARTIAL_CALLBACK_CANDIDATE schema=7 "
                "candidate_index={} site_count={} qp_total={} "
                "sample_count_per_seed=256 cache_requested_bytes=0 "
                "rho=0.01 eta=0.005 partial_order={} floor_multiplier={} "
                "step_count={} callback=max_order_2_bounded "
                "inner_kernel=neighbor_only proposal_policy_hash=5000 "
                "proposal_model_hash=5001 guide_policy_hash=5002 "
                "guide_shape_hash=5003 partial_policy_hash={} "
                "callback_policy_hash={} target_table_hash=6000 "
                "partial_table_hash={} callback_table_hash={} "
                "partial_plan_hash=9000 anchor_count={} "
                "maximum_anchor_relative_residual=0 "
                "maximum_callback_log_residual={} "
                "maximum_callback_weight_relative_residual=0 "
                "exact_balance_pass=1 inner_row_residual=0 "
                "inner_db_residual=0 powered_row_residual=0 "
                "powered_db_residual=0 powered_ratio_residual=0 "
                "outer_db_residual=0 outer_stationary_residual=0 "
                "anchor_pass=1 restart_pass=1\n".format(
                    index, site, qp, partial_order, alpha, step_count,
                    7100 + index, callback_hash, 8000 + index,
                    8100 + index, sector, residual))
            for seed_index, seed_hex in enumerate(SEEDS):
                attempted = samples * step_count
                if bad_census and index == 0 and seed_index == 0:
                    attempted -= 1
                replay = 0 if bad_replay and index == 0 and \
                    seed_index == 0 else 1
                out.write(
                    "PILOT_CALLBACK_SEED candidate_index={} seed_index={} "
                    "seed={} seed_hex={} replay_pass={} accepted=200 "
                    "rejected=56 surrogate_inner_attempted={} "
                    "surrogate_inner_accepted={} surrogate_inner_rejected=0 "
                    "surrogate_final_changed=200 surrogate_final_self=56 "
                    "callback_evaluations={} "
                    "callback_terminal_amplitude_calls=100 "
                    "trace_hash={}\n".format(
                        index, seed_index, int(seed_hex, 16), seed_hex,
                        replay, attempted, attempted,
                        samples * step_count + 1,
                        10000 + index * 4 + seed_index))
            eligible = 0 if bad_decision and index == 0 else 1
            out.write(
                "PILOT_CALLBACK_DECISION candidate_index={} "
                "callback_eligible={} table_pass=1 balance_pass=1 "
                "anchor_pass=1 restart_pass=1 replay_pass=1 "
                "resource_limits_hash=12000\n".format(index, eligible))


def args(policy, output, logs):
    return argparse.Namespace(
        policy=policy, output=output, pilot_log=logs,
        source_commit="39d48c272998521d78166d38acf60c6b2e8624e5",
        producer_binary=sys.executable, allow_smoke_policy=True)


def main():
    if len(sys.argv) != 2:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4s7_stage_b2_evidence.py "
            "POLICY")
    policy_path = os.path.abspath(sys.argv[1])
    policy = evidence.validate_policy(evidence.read_json(policy_path), True)
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4s7-b2-")
    try:
        logs = []
        for site, qp in policy["scope"]["physical_case_order"]:
            path = os.path.join(root, "L{}-QP{}.log".format(site, qp))
            write_log(path, site, qp)
            logs.append(path)
        output = os.path.join(root, "evidence")
        evidence.analyze_command(args(policy_path, output, logs))
        with open(os.path.join(output, "p4s7_stage_b2_eligibility.json"),
                  "r") as handle:
            result = json.load(handle)
        if result["eligible_candidate_indices"] != [0, 1] or \
                result["stage_b2_decision"] != "ELIGIBLE":
            raise AssertionError("P4-S7 B2 synthetic result mismatch")
        negatives = (("hash", {"bad_hash": True}),
                     ("census", {"bad_census": True}),
                     ("anchor", {"bad_anchor": True}),
                     ("replay", {"bad_replay": True}),
                     ("decision", {"bad_decision": True}))
        for label, options in negatives:
            bad = os.path.join(root, "bad-{}.log".format(label))
            write_log(bad, 4, 1, **options)
            expect_failure(
                lambda bad=bad, label=label: evidence.analyze_command(
                    args(policy_path,
                         os.path.join(root,
                                      "bad-{}-evidence".format(label)),
                         [bad] + logs[1:])), label)
        print("bounded Krylov P4-S7 B2 evidence pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
