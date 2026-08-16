from __future__ import print_function

import argparse
import json
import os
import shutil
import sys
import tempfile

import generate_bounded_krylov_p4s7_timing_evidence as evidence


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative P4-S7 timing accepted: {}".format(label))


def write_log(path, site, qp, bad_hash=False, bad_census=False,
              bad_component=False, bad_stream=False):
    samples = 128
    with open(path, "w") as out:
        out.write(
            "TIMING_PARTIAL_RUN schema=2 "
            "fixture=p4s7_partial_callback_matched_timing site_count={} "
            "qp_total={} rank_count=1 sample_count=128 repeat_count=7 "
            "cache_requested_bytes=0 sector_size=36 full_plan_hash=4000 "
            "partial_plan_hash=4001\n".format(site, qp))
        for stream in range(7):
            actual_stream = 0 if bad_stream and stream == 6 else stream
            out.write(
                "TIMING_REPEAT schema=1 site_count={} qp_total={} "
                "sample_count=128 cache_requested_bytes=0 "
                "fraction_index=4 rho_index=0 rho=0.01 "
                "global_numerator=1 global_denominator=1 "
                "seed=5779335701710526001 "
                "seed_hex=0x5034533352305631 rng_stream={} "
                "proposal_policy_hash=5000 proposal_model_hash=5001 "
                "proposal_seconds=1 component_selection_seconds=0 "
                "global_subset_seconds=0.5 bounded_evaluation_seconds=100 "
                "total_step_seconds=128 total_seconds_per_step=1 "
                "accepted=64 rejected=64 neighbor_attempted=0 "
                "global_attempted=128 global_self=1 "
                "final_rng_draws=1000\n".format(
                    site, qp, actual_stream))
            for index, (_, alpha, step_count) in enumerate(
                    evidence.EXPECTED_CANDIDATES):
                total = 134.4 if index == 0 else 139.52
                per_step = total / samples
                callback_hash = 0 if bad_hash and index == 0 and \
                    stream == 0 else 6000 + index
                callback_count = samples * step_count
                if bad_census and index == 0 and stream == 0:
                    callback_count -= 1
                partial_evaluation = total + 1 if bad_component and \
                    index == 0 and stream == 0 else 100
                out.write(
                    "TIMING_PARTIAL_CALLBACK_REPEAT schema=2 "
                    "candidate_index={} site_count={} qp_total={} "
                    "sample_count=128 cache_requested_bytes=0 "
                    "partial_order=2 floor_multiplier={} step_count=32 "
                    "seed=5779335701710526001 "
                    "seed_hex=0x5034533352305631 rng_stream={} "
                    "proposal_policy_hash=5100 proposal_model_hash=5001 "
                    "callback_policy_hash={} partial_plan_hash=4001 "
                    "partial_depth0_seconds=5 partial_depth1_seconds=20 "
                    "partial_depth2_seconds=70 "
                    "partial_evaluation_seconds={} inner_overhead_seconds=5 "
                    "final_bounded_seconds=20 total_step_seconds={} "
                    "total_seconds_per_step={} accepted=120 rejected=8 "
                    "final_changed=124 skipped_full_evaluations=4 "
                    "callback_evaluations={} "
                    "callback_terminal_amplitude_calls=10000 "
                    "final_rng_draws=8321\n".format(
                        index, site, qp, alpha, stream, callback_hash,
                        partial_evaluation, total, per_step,
                        callback_count))


def args(policy, b2, output, logs):
    return argparse.Namespace(
        policy=policy, stage_b2_eligibility=b2, output=output,
        timing_log=logs,
        source_commit="39d48c272998521d78166d38acf60c6b2e8624e5",
        producer_binary=sys.executable, allow_smoke_policy=True)


def main():
    if len(sys.argv) != 2:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4s7_timing_evidence.py POLICY")
    policy_path = os.path.abspath(sys.argv[1])
    policy = evidence.validate_policy(evidence.read_json(policy_path), True)
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4s7-timing-")
    try:
        b2_path = os.path.join(root, "b2.json")
        with open(b2_path, "w") as out:
            json.dump({
                "artifact": "p4s7_stage_b2_eligibility",
                "measurement_policy_sha256":
                    policy["source_stage_b2_policy_sha256"],
                "eligible_candidate_indices": [0, 1]}, out)
        logs = []
        for site, qp in policy["scope"]["physical_case_order"]:
            path = os.path.join(root, "L{}-QP{}.log".format(site, qp))
            write_log(path, site, qp)
            logs.append(path)
        output = os.path.join(root, "evidence")
        evidence.analyze_command(
            args(policy_path, b2_path, output, logs))
        with open(os.path.join(output, "p4s7_post_timing_eligibility.json"),
                  "r") as handle:
            result = json.load(handle)
        if result["eligible_candidate_indices"] != [0, 1] or \
                result["timing_decision"] != "GO":
            raise AssertionError("P4-S7 timing synthetic result mismatch")
        negatives = (("hash", {"bad_hash": True}),
                     ("census", {"bad_census": True}),
                     ("component", {"bad_component": True}),
                     ("stream", {"bad_stream": True}))
        for label, options in negatives:
            bad = os.path.join(root, "bad-{}.log".format(label))
            write_log(bad, 4, 1, **options)
            expect_failure(
                lambda bad=bad, label=label: evidence.analyze_command(
                    args(policy_path, b2_path,
                         os.path.join(root,
                                      "bad-{}-evidence".format(label)),
                         [bad] + logs[1:])), label)
        print("bounded Krylov P4-S7 timing evidence pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
