from __future__ import print_function

import argparse
import os
import shutil
import sys
import tempfile

import generate_bounded_krylov_p4s9_official_resource_evidence as evidence
import runtest_bounded_krylov_p4s9_session_resource_evidence as timing_test
import runtest_bounded_krylov_p4s9_target_session_evidence as target_test


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative official resource accepted: {}".format(
        label))


def make_args(official_policy, timing_policy, target_policy, root,
              timing_logs, target_logs, suffix):
    return argparse.Namespace(
        source_commit="39d48c272998521d78166d38acf60c6b2e8624e5",
        baseline_commit="39d48c272998521d78166d38acf60c6b2e8624e5",
        compiler="synthetic-compiler",
        mpi="synthetic-mpi",
        blas="synthetic-blas",
        strict_fp="synthetic-precise",
        omp_threads=1,
        blas_threads=1,
        repository=None,
        policy=official_policy,
        session_resource_policy=timing_policy,
        target_session_policy=target_policy,
        stage_c_local_decision=official_policy,
        target_local_evidence=official_policy,
        target_local_decision=official_policy,
        producer_timing_binary=sys.executable,
        producer_profile_binary=sys.executable,
        output=os.path.join(root, suffix),
        timing_log=timing_logs,
        profile_log=target_logs,
        allow_smoke_policy=True,
    )


def replace_once(path, old, new):
    with open(path, "r") as handle:
        text = handle.read()
    if text.count(old) != 1:
        raise AssertionError("test mutation target mismatch")
    with open(path, "w") as handle:
        handle.write(text.replace(old, new, 1))


def write_cross_rank_timing_logs(root, policy, label,
                                 mismatch_rank=None):
    logs = []
    for rank in (1, 2, 4):
        for site, qp in policy["scope"]["physical_case_order"]:
            path = os.path.join(
                root, "{}-timing-R{}-L{}-QP{}.log".format(
                    label, rank, site, qp))
            timing_test.write_log(
                path, policy, site, qp,
                trace_delta=1 if rank == mismatch_rank else 0)
            if rank != 1:
                replace_once(path, "rank_count=1", "rank_count={}".format(
                    rank))
            logs.append(path)
    return logs


def write_cross_rank_target_logs(root, label, mismatch_rank=None):
    logs = []
    for rank in (1, 2, 4):
        for qp in (1, 4):
            for repeat in range(3):
                path = os.path.join(
                    root, "{}-target-R{}-QP{}-N{}.log".format(
                        label, rank, qp, repeat))
                target_test.write_log(
                    path, qp,
                    row_delta=1.0 if rank == mismatch_rank else 0.0)
                if rank != 1:
                    replace_once(
                        path, "rank_count=1", "rank_count={}".format(rank))
                logs.append(path)
    return logs


def main():
    if len(sys.argv) != 4:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4s9_official_resource_"
            "evidence.py OFFICIAL_POLICY TIMING_POLICY TARGET_POLICY")
    official_policy_path = os.path.abspath(sys.argv[1])
    timing_policy_path = os.path.abspath(sys.argv[2])
    target_policy_path = os.path.abspath(sys.argv[3])
    official_policy = evidence.official.validate_policy(
        evidence.p4c.read_json(official_policy_path))
    evidence.validate_resource_scope(official_policy)
    timing_policy = evidence.session_resource.validate_policy(
        evidence.p4c.read_json(timing_policy_path))
    evidence.target_session.validate_policy(
        evidence.p4c.read_json(target_policy_path))
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4s9-official-resource-")
    try:
        timing_logs = []
        for site, qp in timing_policy["scope"]["physical_case_order"]:
            path = os.path.join(root, "timing-L{}-QP{}.log".format(site, qp))
            timing_test.write_log(path, timing_policy, site, qp)
            timing_logs.append(path)
        target_logs = []
        for qp in (1, 4):
            for repeat in range(3):
                path = os.path.join(
                    root, "target-QP{}-R{}.log".format(qp, repeat))
                target_test.write_log(path, qp)
                target_logs.append(path)
        args = make_args(
            official_policy_path, timing_policy_path, target_policy_path,
            root, timing_logs, target_logs, "positive")
        evidence.analyze_command(args)
        decision = evidence.p4c.read_json(os.path.join(
            args.output, "p4s9_official_resource_decision.json"))
        if decision["official_resource_decision"] != "PASS" or \
                decision["p4f_authorized"]:
            raise AssertionError("synthetic official resource PASS mismatch")

        slow_logs = []
        for site, qp in timing_policy["scope"]["physical_case_order"]:
            path = os.path.join(
                root, "slow-L{}-QP{}.log".format(site, qp))
            timing_test.write_log(path, timing_policy, site, qp, slow=True)
            slow_logs.append(path)
        slow_args = make_args(
            official_policy_path, timing_policy_path, target_policy_path,
            root, slow_logs, target_logs, "slow")
        evidence.analyze_command(slow_args)
        slow_decision = evidence.p4c.read_json(os.path.join(
            slow_args.output, "p4s9_official_resource_decision.json"))
        if slow_decision["official_resource_decision"] != "STOP":
            raise AssertionError("slow official resource did not STOP")

        trace_mismatch = os.path.join(root, "target-trace-mismatch.log")
        target_test.write_log(trace_mismatch, 1, row_delta=1.0)
        trace_args = make_args(
            official_policy_path, timing_policy_path, target_policy_path,
            root, timing_logs, [trace_mismatch] + target_logs[1:],
            "trace-mismatch")
        evidence.analyze_command(trace_args)
        trace_decision = evidence.p4c.read_json(os.path.join(
            trace_args.output, "p4s9_official_resource_decision.json"))
        if trace_decision["official_resource_decision"] != "STOP":
            raise AssertionError("target trace mismatch did not STOP")

        bad_rank = os.path.join(root, "bad-rank.log")
        shutil.copyfile(timing_logs[0], bad_rank)
        replace_once(bad_rank, "rank_count=1", "rank_count=3")
        bad_logs = [bad_rank] + timing_logs[1:]
        expect_failure(lambda: evidence.analyze_command(make_args(
            official_policy_path, timing_policy_path, target_policy_path,
            root, bad_logs, target_logs, "bad-rank")), "rank scope")

        duplicate_target_logs = (
            [target_logs[0], target_logs[0], target_logs[0]] +
            target_logs[3:])
        expect_failure(lambda: evidence.analyze_command(make_args(
            official_policy_path, timing_policy_path, target_policy_path,
            root, timing_logs, duplicate_target_logs,
            "duplicate-target")), "duplicate target log path")

        bad_provenance = make_args(
            official_policy_path, timing_policy_path, target_policy_path,
            root, timing_logs, target_logs, "bad-provenance")
        bad_provenance.baseline_commit = "39d48c2"
        expect_failure(lambda: evidence.analyze_command(bad_provenance),
                       "baseline commit provenance")

        cross_timing = write_cross_rank_timing_logs(
            root, timing_policy, "cross-good")
        cross_target = write_cross_rank_target_logs(root, "cross-good")
        cross_args = make_args(
            official_policy_path, timing_policy_path, target_policy_path,
            root, cross_timing, cross_target, "cross-rank-positive")
        evidence.analyze_command(cross_args)
        cross_decision = evidence.p4c.read_json(os.path.join(
            cross_args.output, "p4s9_official_resource_decision.json"))
        if cross_decision["official_resource_decision"] != "PASS":
            raise AssertionError("cross-rank resource positive mismatch")

        bad_cross_timing = write_cross_rank_timing_logs(
            root, timing_policy, "cross-timing-bad", mismatch_rank=2)
        bad_cross_timing_args = make_args(
            official_policy_path, timing_policy_path, target_policy_path,
            root, bad_cross_timing, cross_target,
            "cross-rank-timing-trace-stop")
        evidence.analyze_command(bad_cross_timing_args)
        bad_cross_timing_decision = evidence.p4c.read_json(os.path.join(
            bad_cross_timing_args.output,
            "p4s9_official_resource_decision.json"))
        if bad_cross_timing_decision["official_resource_decision"] != "STOP":
            raise AssertionError("timing cross-rank trace mismatch did not STOP")

        bad_cross_target = write_cross_rank_target_logs(
            root, "cross-target-bad", mismatch_rank=2)
        bad_cross_target_args = make_args(
            official_policy_path, timing_policy_path, target_policy_path,
            root, cross_timing, bad_cross_target,
            "cross-rank-target-trace-stop")
        evidence.analyze_command(bad_cross_target_args)
        bad_cross_target_decision = evidence.p4c.read_json(os.path.join(
            bad_cross_target_args.output,
            "p4s9_official_resource_decision.json"))
        if bad_cross_target_decision["official_resource_decision"] != "STOP":
            raise AssertionError("target cross-rank trace mismatch did not STOP")
        print("bounded Krylov P4-S9 official resource pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
