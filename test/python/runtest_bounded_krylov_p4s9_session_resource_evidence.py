from __future__ import print_function

import argparse
import json
import os
import shutil
import sys
import tempfile

import generate_bounded_krylov_p4s9_session_resource_evidence as evidence


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative P4-S9 resource accepted: {}".format(label))


def write_log(path, policy, site, qp, bad_stream=False, bad_reset=False,
              bad_correctness=False, bad_generation=False,
              slow=False, trace_delta=0):
    sample_count = 16
    cache_bytes = 4194304
    sector_size = evidence.p4c.fixed_sector_size(site)
    base_per_configuration = (0.2 if slow else 1.0e-5) * \
        (1.2 ** (site - 4)) * qp
    with open(path, "w") as out:
        out.write(
            "TIMING_SESSION_RUN schema=3 "
            "fixture=p4s9_long_direct_session_timing site_count={} "
            "qp_total={} rank_count=1 sample_count={} repeat_count=7 "
            "cache_requested_bytes={} sector_size={} plan_hash={}\n".format(
                site, qp, sample_count, cache_bytes, sector_size,
                1000 + site * 10 + qp))
        for stream in range(7):
            output_stream = 0 if bad_stream and stream == 6 else stream
            total_step = base_per_configuration * sample_count * 0.9
            out.write(
                "TIMING_SESSION_REPEAT schema=3 site_count={} qp_total={} "
                "sample_count={} cache_requested_bytes={} fraction_index=0 "
                "rho_index=0 rho=0.01 global_numerator=0 "
                "global_denominator=1 seed={} seed_hex={} rng_stream={} "
                "proposal_policy_hash=101 proposal_model_hash=102 "
                "proposal_seconds=0 component_selection_seconds=0 "
                "global_subset_seconds=0 bounded_evaluation_seconds={} "
                "total_step_seconds={} total_seconds_per_step={} "
                "accepted=8 rejected=8 neighbor_attempted=16 "
                "global_attempted=0 global_self=0 final_rng_draws=33\n".format(
                    site, qp, sample_count, cache_bytes,
                    int(policy["driver"]["seed_hex"], 16),
                    policy["driver"]["seed_hex"], output_stream,
                    total_step * 0.8, total_step,
                    total_step / sample_count))
            combined = base_per_configuration * sample_count
            generation = 999 if not bad_generation or stream != 6 else 998
            resets = 0 if bad_reset and stream == 0 else 1
            out.write(
                "TIMING_SESSION_RESOURCE schema=3 site_count={} qp_total={} "
                "sample_count={} cache_requested_bytes={} rng_stream={} "
                "amplitude_generation_hash={} session_begin_seconds=0 "
                "initialization_seconds={} combined_session_seconds={} "
                "seconds_per_configuration={} session_root_evaluations=17 "
                "cache_resets={} cache_hits=100 cache_misses=50 "
                "cache_insertions=50 cache_evictions=0 "
                "cache_resident_peak=50 callback_calls=36 "
                "allocated_capacity_bytes=270000000 "
                "peak_rss_bytes={}\n".format(
                    site, qp, sample_count, cache_bytes, output_stream,
                    generation, combined * 0.1, combined,
                    base_per_configuration, resets,
                    180000000 + site * 1000000 + qp * 100000))
            exact_match = 0 if bad_correctness and stream == 0 else 1
            out.write(
                "TIMING_SESSION_CORRECTNESS schema=3 site_count={} "
                "qp_total={} rng_stream={} exact_table_match={} "
                "final_sector_index=0 final_configuration=1 "
                "statistical_trace_hash={}\n".format(
                    site, qp, output_stream, exact_match,
                    500 + stream + trace_delta))


def make_args(policy, root, logs, suffix):
    parent = os.path.join(root, "parent.json")
    stage_a_eligibility = os.path.join(root, "stage-a-eligibility.json")
    stage_a_evidence = os.path.join(root, "stage-a-evidence.json")
    p4c_resource = os.path.join(root, "p4c-resource.json")
    return argparse.Namespace(
        policy=policy, parent_policy=parent,
        stage_a_eligibility=stage_a_eligibility,
        stage_a_evidence=stage_a_evidence, p4c_resource=p4c_resource,
        output=os.path.join(root, suffix), timing_log=logs,
        source_commit="39d48c272998521d78166d38acf60c6b2e8624e5",
        producer_binary=sys.executable, allow_smoke_policy=True)


def main():
    if len(sys.argv) != 2:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4s9_session_resource_evidence.py "
            "POLICY")
    policy_path = os.path.abspath(sys.argv[1])
    policy = evidence.validate_policy(evidence.p4c.read_json(policy_path))
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4s9-resource-")
    try:
        for name, value in (
                ("parent.json", {}),
                ("stage-a-eligibility.json",
                 {"stage_a_decision": "ELIGIBLE",
                  "stage_b_authorized": True}),
                ("stage-a-evidence.json", {}),
                ("p4c-resource.json", {})):
            with open(os.path.join(root, name), "w") as handle:
                json.dump(value, handle)
        logs = []
        for site, qp in policy["scope"]["physical_case_order"]:
            path = os.path.join(root, "L{}-QP{}.log".format(site, qp))
            write_log(path, policy, site, qp)
            logs.append(path)
        args = make_args(policy_path, root, logs, "positive")
        evidence.analyze_command(args)
        with open(os.path.join(
                args.output,
                "p4s9_local_session_resource_decision.json")) as handle:
            decision = json.load(handle)
        if decision["local_resource_decision"] != "PASS" or \
                not decision["remote_package_preparation_authorized"] or \
                decision["remote_transfer_or_submission_authorized"]:
            raise AssertionError("P4-S9 synthetic resource PASS mismatch")

        negatives = (
            ("stream", {"bad_stream": True}),
            ("reset", {"bad_reset": True}),
            ("correctness", {"bad_correctness": True}),
            ("generation", {"bad_generation": True}),
        )
        for label, options in negatives:
            bad = os.path.join(root, "bad-{}.log".format(label))
            write_log(bad, policy, 4, 1, **options)
            bad_logs = [bad] + logs[1:]
            expect_failure(
                lambda label=label, bad_logs=bad_logs:
                evidence.analyze_command(make_args(
                    policy_path, root, bad_logs, "bad-{}".format(label))),
                label)

        slow_logs = []
        for site, qp in policy["scope"]["physical_case_order"]:
            path = os.path.join(
                root, "slow-L{}-QP{}.log".format(site, qp))
            write_log(path, policy, site, qp, slow=True)
            slow_logs.append(path)
        slow_args = make_args(policy_path, root, slow_logs, "slow")
        evidence.analyze_command(slow_args)
        with open(os.path.join(
                slow_args.output,
                "p4s9_local_session_resource_decision.json")) as handle:
            slow_decision = json.load(handle)
        if slow_decision["local_resource_decision"] != "STOP" or \
                slow_decision["remote_package_preparation_authorized"]:
            raise AssertionError("P4-S9 resource threshold did not STOP")
        print("bounded Krylov P4-S9 session resource pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
