from __future__ import print_function

import argparse
import json
import os
import shutil
import sys
import tempfile

import generate_bounded_krylov_p4s5_pilot_evidence as evidence


SEEDS = [
    "0x5034533550494c31", "0x5034533550494c32",
    "0x5034533550494c33", "0x5034533550494c34",
]


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative P4-S5 evidence accepted: {}".format(label))


def diagnostic_fields(scale=1.0):
    iid_ratio = 0.45 * scale
    required_tau = (0.90 / iid_ratio) ** 2
    return (
        "guide_minimum=1 guide_maximum=2 guide_mean=1.5 "
        "guide_maximum_minimum_ratio=2 "
        "uniform_sector_self_probability={self_probability} "
        "importance_weight_maximum_mean_ratio=1.5 "
        "importance_weight_ess_fraction=0.75 "
        "maximum_iid_se_budget_ratio={iid_ratio:.17g} "
        "required_tau_for_0p90={required_tau:.17g}"
    ).format(self_probability="{self_probability}",
             iid_ratio=iid_ratio, required_tau=required_tau)


def write_synthetic_log(path, site, qp, policy, duplicate_guide_hash=False,
                        bad_control=False, bad_required_tau=False,
                        bad_decision=False):
    sample_count = 256
    sector_size = evidence.fixed_sector_size(site)
    eta_per_rho = 2.5
    diagnostic = diagnostic_fields()
    with open(path, "w") as handle:
        handle.write(
            "PILOT_RUN schema=3 fixture=p4s5_flat_guide_exact_table "
            "site_count={} qp_total={} rank_count=1 "
            "sample_count_per_seed={} cache_requested_bytes=0 "
            "sector_size={} anchor_count={} plan_hash=100\n".format(
                site, qp, sample_count, sector_size, sector_size))
        for control_index, control_id in enumerate(evidence.CONTROL_IDS):
            uniform = control_index == 1
            promotable = 1 if bad_control and control_index == 0 else 0
            fields = diagnostic.format(
                self_probability="{:.17g}".format(1.0 / sector_size))
            if uniform:
                fields = fields.replace(
                    "guide_minimum=1 guide_maximum=2 guide_mean=1.5 "
                    "guide_maximum_minimum_ratio=2",
                    "guide_minimum=1 guide_maximum=1 guide_mean=1 "
                    "guide_maximum_minimum_ratio=1")
            handle.write(
                "PILOT_CONTROL schema=3 control_index={} control_id={} "
                "promotable={} uniform_guide={} rho={} eta={} "
                "guide_policy_hash={} guide_shape_hash={} table_hash={} "
                "anchor_count={} {} limiting_family=S limiting_entry=0\n".
                format(control_index, control_id, promotable, int(uniform),
                       0 if uniform else 0.01,
                       1 if uniform else 0.01 * eta_per_rho,
                       1000 + control_index, 2000 + control_index,
                       3000 + control_index, sector_size, fields))
        for index, rho in enumerate(policy["guide"]["candidate_rho_order"]):
            guide_hash = 4000 + (0 if duplicate_guide_hash and index == 1
                                 else index)
            fields = diagnostic.format(
                self_probability="{:.17g}".format(1.0 / sector_size))
            if bad_required_tau and index == 2:
                fields = fields.replace("required_tau_for_0p90=4",
                                        "required_tau_for_0p90=5")
            handle.write(
                "PILOT_CANDIDATE schema=3 candidate_index={} site_count={} "
                "qp_total={} sample_count_per_seed={} cache_requested_bytes=0 "
                "rho={} eta={} proposal_kernel=123 global_numerator=1 "
                "global_denominator=1 neighbor_numerator=0 "
                "neighbor_denominator=0 distance_numerator=0 "
                "distance_denominator=0 distance_rounding_rule=0 "
                "shell_max_distance=0 shell_distance=0 shell_count=0 "
                "proposal_policy_hash=5000 proposal_model_hash=6000 "
                "guide_policy_hash={} guide_shape_hash={} table_hash={} "
                "anchor_count={} {} iid_limiting_family=S "
                "iid_limiting_entry=0 exact_balance_pass=1 db_residual=0 "
                "stationary_residual=0 proposal_row_residual=0 "
                "proposal_symmetry_residual=0 shell_cardinality_pass=1 "
                "restart_pass=1\n".format(
                    index, site, qp, sample_count, rho,
                    float(rho) * eta_per_rho, guide_hash, 7000 + index,
                    8000 + index, sector_size, fields))
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
                    "maximum_predicted_se_budget_ratio=0.5 accepted=200 "
                    "rejected=56 neighbor_attempted=0 global_attempted=256 "
                    "global_self=7 shell_attempted=0 trace_hash={}\n".
                    format(index, seed_index, int(seed_hex, 16), seed_hex,
                           9000 + index * 4 + seed_index))
            eligible = 0 if bad_decision and index == 0 else 1
            handle.write(
                "PILOT_DECISION candidate_index={} physical_case_eligible={} "
                "maximum_predicted_se_budget_ratio=0.5 maximum_tau_int=1 "
                "table_pass=1 balance_pass=1 restart_pass=1 "
                "resource_limits_hash=10000\n".format(index, eligible))


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
            "usage: runtest_bounded_krylov_p4s5_pilot_evidence.py POLICY")
    policy_path = os.path.abspath(sys.argv[1])
    policy = evidence.validate_policy(evidence.read_json(policy_path), True)
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4s5-pilot-")
    try:
        logs = []
        for site, qp in policy["scope"]["physical_case_order"]:
            path = os.path.join(root, "L{}-QP{}.log".format(site, qp))
            write_synthetic_log(path, site, qp, policy)
            logs.append(path)
        output = os.path.join(root, "evidence")
        evidence.analyze_command(analyze_args(policy_path, output, logs))
        with open(os.path.join(output, "p4s5_pilot_eligibility.json"),
                  "r") as handle:
            eligibility = json.load(handle)
        if eligibility["eligible_candidate_indices"] != list(range(7)) or \
                eligibility["preferred_candidate_index"] != 0 or \
                eligibility["pilot_decision"] != "ELIGIBLE":
            raise AssertionError("synthetic P4-S5 eligibility mismatch")

        negative_cases = (
            ("guide-hash", {"duplicate_guide_hash": True}),
            ("control-role", {"bad_control": True}),
            ("required-tau", {"bad_required_tau": True}),
            ("decision", {"bad_decision": True}),
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
        print("bounded Krylov P4-S5 pilot evidence pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
