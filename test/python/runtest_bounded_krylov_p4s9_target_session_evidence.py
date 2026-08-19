from __future__ import print_function

import argparse
import json
import os
import shutil
import sys
import tempfile

import generate_bounded_krylov_p4s9_target_session_evidence as evidence


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative P4-S9 target accepted: {}".format(label))


def write_log(path, qp, bad_schedule=False, bad_reset=False, slow=False,
              row_delta=0.0):
    sample_count = 3
    cache_bytes = 4096
    rows = evidence.p4c.expected_profile_rows(16, sample_count)
    seconds_per_root = 0.1 if slow else 0.001
    with open(path, "w") as out:
        out.write(
            "SESSION_PROFILE schema=2 "
            "fixture=p4s9_target_session_spread_roots site_count=16 "
            "qp_total={} sample_count={} sector_size=165636900 rank_count=1 "
            "cache_requested_bytes={} audit=1\n".format(
                qp, sample_count, cache_bytes))
        resident = 0
        for index, sector_index, configuration in rows:
            if bad_schedule and index == 1:
                sector_index += 1
            values = [float(qp), row_delta, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            out.write("ROW {} {} {} {}\n".format(
                index, sector_index, configuration,
                " ".join(str(value) for value in values)))
            reset = 0 if bad_reset and index == 0 else (1 if index == 0 else 0)
            end = resident + 10
            out.write(
                "STAT {} reset_seconds=0 total_seconds={} "
                "depth_seconds=0,0,0,0 terminal_calls=10 "
                "raw_transitions=10 cache_hits=1,1,1,1 "
                "cache_misses=1,1,1,1 row_peak=4 cache_entries_peak={} "
                "session_active=1 session_generation=999 "
                "session_root_evaluation={} cache_reset_performed={} "
                "cache_resident_start={} cache_resident_end={} "
                "workspace_bytes=10000 cache_allocated_bytes=3000\n".format(
                    index, seconds_per_root, end, index + 1, reset,
                    resident, end))
            resident = end
        resource = (
            "roots=3 node_expansions=30 recursion=3,3,3,3 "
            "cache_hits=3,3,3,3 cache_misses=3,3,3,3 "
            "cache_insertions=30 cache_evictions=0 cache_epoch_clears=0 "
            "raw_transitions=30 merged_duplicates=0 cancelled_zero=0 "
            "terminal_calls=30 well_pivoted=30 near_pivot=0 "
            "exact_zero_components=0 numeric_zero_components=0 "
            "local_factorizations=30 global_factorizations=30 "
            "computed_exact_zero=0 computed_numeric_zero=0 trace_hash=123 "
            "row_peak=4 cache_entries_peak=30 plan_bytes=1000 "
            "model_bytes=500 engine_workspace_bytes=10000 frame_bytes=100 "
            "scratch_bytes=8 cache_requested_bytes=4096 "
            "cache_allocated_bytes=3000 cache_set_count=10 "
            "collective_workspace_bytes=500 amplitude_workspace_bytes=500 "
            "engine_heap_allocations=0 rss_bytes=100000000 reset_seconds=0 "
            "evaluation_seconds={} depth_seconds=0,0,0,0 "
            "amplitude_seconds=0 connectivity_seconds=0 cache_seconds=0 "
            "total_seconds={}")
        total = seconds_per_root * sample_count
        out.write("RESOURCE scope=rank_sum {}\n".format(
            resource.format(total, total)))
        out.write("RESOURCE scope=rank_max {}\n".format(
            resource.format(total, total)))


def make_args(policy, root, logs, suffix):
    return argparse.Namespace(
        policy=policy,
        parent_resource_policy=os.path.join(root, "parent.json"),
        stage_c_evidence=os.path.join(root, "stage-c-evidence.json"),
        stage_c_decision=os.path.join(root, "stage-c-decision.json"),
        output=os.path.join(root, suffix), profile_log=logs,
        source_commit="39d48c272998521d78166d38acf60c6b2e8624e5",
        producer_binary=sys.executable, allow_smoke_policy=True)


def main():
    if len(sys.argv) != 2:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4s9_target_session_evidence.py "
            "POLICY")
    policy_path = os.path.abspath(sys.argv[1])
    evidence.validate_policy(evidence.p4c.read_json(policy_path))
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4s9-target-")
    try:
        for name in ("parent.json", "stage-c-evidence.json",
                     "stage-c-decision.json"):
            with open(os.path.join(root, name), "w") as handle:
                json.dump({}, handle)
        logs = []
        for qp in (1, 4):
            for repeat in range(3):
                path = os.path.join(
                    root, "QP{}-R{}.log".format(qp, repeat))
                write_log(path, qp)
                logs.append(path)
        args = make_args(policy_path, root, logs, "positive")
        evidence.analyze_command(args)
        with open(os.path.join(
                args.output,
                "p4s9_local_target_session_decision.json")) as handle:
            decision = json.load(handle)
        if decision["local_target_stress_decision"] != "PASS" or \
                not decision["remote_package_preparation_authorized"] or \
                decision["remote_transfer_or_submission_authorized"]:
            raise AssertionError("P4-S9 target synthetic PASS mismatch")

        bad_schedule = os.path.join(root, "bad-schedule.log")
        write_log(bad_schedule, 1, bad_schedule=True)
        expect_failure(lambda: evidence.analyze_command(make_args(
            policy_path, root, [bad_schedule] + logs[1:], "bad-schedule")),
                       "schedule")
        bad_reset = os.path.join(root, "bad-reset.log")
        write_log(bad_reset, 1, bad_reset=True)
        expect_failure(lambda: evidence.analyze_command(make_args(
            policy_path, root, [bad_reset] + logs[1:], "bad-reset")),
                       "reset")

        trace_mismatch = os.path.join(root, "trace-mismatch.log")
        write_log(trace_mismatch, 1, row_delta=1.0)
        trace_args = make_args(
            policy_path, root, [trace_mismatch] + logs[1:], "trace-stop")
        evidence.analyze_command(trace_args)
        with open(os.path.join(
                trace_args.output,
                "p4s9_local_target_session_decision.json")) as handle:
            trace_decision = json.load(handle)
        if trace_decision["local_target_stress_decision"] != "STOP":
            raise AssertionError("P4-S9 target trace mismatch did not STOP")

        slow_logs = []
        for qp in (1, 4):
            for repeat in range(3):
                path = os.path.join(
                    root, "slow-QP{}-R{}.log".format(qp, repeat))
                write_log(path, qp, slow=True)
                slow_logs.append(path)
        slow_args = make_args(policy_path, root, slow_logs, "slow-stop")
        evidence.analyze_command(slow_args)
        with open(os.path.join(
                slow_args.output,
                "p4s9_local_target_session_decision.json")) as handle:
            slow_decision = json.load(handle)
        if slow_decision["local_target_stress_decision"] != "STOP":
            raise AssertionError("P4-S9 target runtime did not STOP")
        print("bounded Krylov P4-S9 target session pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
