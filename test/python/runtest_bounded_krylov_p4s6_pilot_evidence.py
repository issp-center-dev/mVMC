from __future__ import print_function

import argparse
import json
import os
import shutil
import sys
import tempfile

import generate_bounded_krylov_p4s6_pilot_evidence as evidence


SEEDS = [
    "0x5034533650494c31", "0x5034533650494c32",
    "0x5034533650494c33", "0x5034533650494c34",
]


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative P4-S6 evidence accepted: {}".format(label))


def diagnostics():
    return (
        "surrogate_minimum=1 surrogate_maximum=2 surrogate_mean=1.5 "
        "target_ratio_minimum=0.5 target_ratio_maximum=2 "
        "target_ratio_maximum_mean_ratio=1.2 "
        "target_ratio_ess_fraction=0.8 "
        "maximum_iid_se_budget_ratio=0.45 required_tau_for_0p90=4 "
    )


def write_synthetic_log(path, site, qp, policy, bad_policy_hash=False,
                        bad_table_hash=False, bad_census=False,
                        bad_balance=False, bad_decision=False,
                        bad_control=False):
    sample_count = 256
    sector_size = evidence.fixed_sector_size(site)
    with open(path, "w") as handle:
        handle.write(
            "PILOT_RUN schema=4 fixture=p4s6_surrogate_exact_table "
            "site_count={} qp_total={} rank_count=1 "
            "sample_count_per_seed={} cache_requested_bytes=0 "
            "sector_size={} anchor_count={} plan_hash=4000\n".format(
                site, qp, sample_count, sector_size, sector_size))
        for alpha_index, alpha in enumerate(evidence.EXPECTED_ALPHAS):
            observed_alpha = 7 if bad_control and alpha_index == 0 else alpha
            handle.write(
                "PILOT_SURROGATE_CONTROL schema=4 control_index={} "
                "control_id=exact_q_independence promotable=0 alpha={} "
                "rho=0.01 eta=0.005 guide_policy_hash=5000 "
                "guide_shape_hash=5001 surrogate_policy_hash={} "
                "target_table_hash=6000 surrogate_table_hash={} "
                "anchor_count={} {}ideal_acceptance=0.7 "
                "ideal_configuration_change=0.6 limiting_family=S "
                "limiting_entry=0\n".format(
                    alpha_index, observed_alpha, 7000 + alpha_index,
                    9000 + alpha_index, sector_size, diagnostics()))
        for index, (step_count, alpha) in enumerate(
                evidence.EXPECTED_CANDIDATES):
            alpha_index = evidence.EXPECTED_ALPHAS.index(alpha)
            surrogate_policy_hash = 8000 + index
            if bad_policy_hash and index == 1:
                surrogate_policy_hash = 8000
            surrogate_table_hash = 9000 + alpha_index
            if bad_table_hash and index == 3:
                surrogate_table_hash = 9999
            outer_residual = 1 if bad_balance and index == 0 else 0
            handle.write(
                "PILOT_CANDIDATE schema=4 candidate_index={} "
                "site_count={} qp_total={} sample_count_per_seed={} "
                "cache_requested_bytes=0 rho=0.01 eta=0.005 "
                "surrogate_step_count={} floor_multiplier={} "
                "inner_kernel=neighbor_only proposal_policy_hash=4001 "
                "proposal_model_hash=4002 guide_policy_hash=5000 "
                "guide_shape_hash=5001 surrogate_policy_hash={} "
                "target_table_hash=6000 surrogate_table_hash={} "
                "anchor_count={} {}iid_limiting_family=S "
                "iid_limiting_entry=0 exact_balance_pass=1 "
                "inner_row_residual=0 inner_db_residual=0 "
                "powered_row_residual=0 powered_db_residual=0 "
                "powered_ratio_residual=0 outer_db_residual={} "
                "outer_stationary_residual=0 restart_pass=1\n".format(
                    index, site, qp, sample_count, step_count, alpha,
                    surrogate_policy_hash, surrogate_table_hash, sector_size,
                    diagnostics(), outer_residual))
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
                attempted = sample_count * step_count
                if bad_census and index == 0 and seed_index == 0:
                    attempted -= 1
                accepted_inner = attempted // 2
                rejected_inner = attempted - accepted_inner
                handle.write(
                    "PILOT_SEED candidate_index={} seed_index={} seed={} "
                    "seed_hex={} maximum_tau_int=1 "
                    "maximum_predicted_se_budget_ratio=0.5 accepted=200 "
                    "rejected=56 surrogate_inner_attempted={} "
                    "surrogate_inner_accepted={} surrogate_inner_rejected={} "
                    "surrogate_final_changed=250 surrogate_final_self=6 "
                    "trace_hash={}\n".format(
                        index, seed_index, int(seed_hex, 16), seed_hex,
                        attempted, accepted_inner, rejected_inner,
                        10000 + index * 4 + seed_index))
            eligible = 0 if bad_decision and index == 0 else 1
            handle.write(
                "PILOT_DECISION candidate_index={} physical_case_eligible={} "
                "maximum_predicted_se_budget_ratio=0.5 maximum_tau_int=1 "
                "table_pass=1 balance_pass=1 restart_pass=1 "
                "resource_limits_hash=12000\n".format(index, eligible))


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
            "usage: runtest_bounded_krylov_p4s6_pilot_evidence.py POLICY")
    policy_path = os.path.abspath(sys.argv[1])
    policy = evidence.validate_policy(evidence.read_json(policy_path), True)
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4s6-pilot-")
    try:
        logs = []
        for site, qp in policy["scope"]["physical_case_order"]:
            path = os.path.join(root, "L{}-QP{}.log".format(site, qp))
            write_synthetic_log(path, site, qp, policy)
            logs.append(path)
        output = os.path.join(root, "evidence")
        evidence.analyze_command(analyze_args(policy_path, output, logs))
        with open(os.path.join(output, "p4s6_pilot_eligibility.json"),
                  "r") as handle:
            eligibility = json.load(handle)
        if eligibility["eligible_candidate_indices"] != list(range(9)) or \
                eligibility["preferred_candidate_index"] != 0 or \
                eligibility["pilot_decision"] != "ELIGIBLE":
            raise AssertionError("synthetic P4-S6 eligibility mismatch")

        negative_cases = (
            ("policy-hash", {"bad_policy_hash": True}),
            ("table-hash", {"bad_table_hash": True}),
            ("census", {"bad_census": True}),
            ("balance", {"bad_balance": True}),
            ("decision", {"bad_decision": True}),
            ("control", {"bad_control": True}),
        )
        for label, options in negative_cases:
            path = os.path.join(root, "bad-{}.log".format(label))
            write_synthetic_log(path, 4, 1, policy, **options)
            expect_failure(
                lambda path=path, label=label: evidence.analyze_command(
                    analyze_args(
                        policy_path,
                        os.path.join(root, "bad-{}-evidence".format(label)),
                        [path] + logs[1:])),
                label)
        print("bounded Krylov P4-S6 pilot evidence pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
