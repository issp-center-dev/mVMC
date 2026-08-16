from __future__ import print_function

import argparse
import os
import shutil
import sys
import tempfile

import generate_bounded_krylov_p4s9_official_closure as closure
import runtest_bounded_krylov_p4s9_official_statistics_evidence as stats_test


COMMIT = "39d48c272998521d78166d38acf60c6b2e8624e5"
CANDIDATE = "direct_neighbor_rho0p01_N32768_session256m"


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative closure accepted: {}".format(label))


def write_component(root, policy_hash, kind, passed=True):
    os.makedirs(root)
    if kind == "statistics":
        evidence_name = "p4s9_official_markov_statistics.json"
        decision_name = "p4s9_official_statistics_decision.json"
        evidence_flag = "official_statistics_pass"
        decision_field = "official_statistics_decision"
        decision_value = "GO" if passed else "STOP"
    else:
        evidence_name = "p4s9_official_resource_evidence.json"
        decision_name = "p4s9_official_resource_decision.json"
        evidence_flag = "official_resource_pass"
        decision_field = "official_resource_decision"
        decision_value = "PASS" if passed else "STOP"
    base = {
        "schema_version": 1,
        "source_commit": COMMIT,
        "baseline_develop_commit": COMMIT,
        "official_execution_policy_sha256": policy_hash,
        "candidate_id": CANDIDATE,
    }
    evidence = dict(base)
    evidence["artifact"] = (
        "p4s9_official_markov_statistics" if kind == "statistics" else
        "p4s9_official_resource_evidence")
    evidence[evidence_flag] = passed
    decision = dict(base)
    decision["artifact"] = (
        "p4s9_official_statistics_decision" if kind == "statistics" else
        "p4s9_official_resource_decision")
    decision[decision_field] = decision_value
    decision["p4f_authorized"] = False
    closure.p4c.write_json(os.path.join(root, evidence_name), evidence)
    closure.p4c.write_json(os.path.join(root, decision_name), decision)
    closure.p4c.write_checksums(root)


def write_parity_log(path, rank):
    stats_test.write_session_log(path, sample_count=128)
    with open(path, "r") as handle:
        text = handle.read()
    text = text.replace("qp_total=1", "qp_total=4")
    text = text.replace("cache_requested_bytes=268435456",
                        "cache_requested_bytes=4194304")
    text = text.replace("rng_stream=0", "rng_stream=7")
    text = text.replace("rank_count=1", "rank_count={}".format(rank), 1)
    text = text.replace("pass=1 required=", "pass=0 required=")
    text = text.replace(
        "DECISION p4s_decision=GO support_pass=1 tau_pass=1",
        "DECISION p4s_decision=STOP support_pass=1 tau_pass=0")
    with open(path, "w") as handle:
        handle.write(text)


def write_ctest(path, policy, omit=None, fail=None):
    names = [name for name in policy["focused_tests"]["required_names"]
             if name != omit]
    with open(path, "w") as handle:
        for index, name in enumerate(names, 1):
            status = "***Failed" if name == fail else "Passed"
            handle.write(
                "{}/{} Test #{}: {} ... {} 0.01 sec\n".format(
                    index, len(names), index, name, status))
        handle.write(
            "100% tests passed, 0 tests failed out of {}\n".format(
                len(names)))


def make_args(policy, stats, resource, parity, ctest, output):
    return argparse.Namespace(
        policy=policy,
        statistics_evidence=stats,
        resource_evidence=resource,
        mpi_parity_log=parity,
        focused_ctest_log=ctest,
        output=output,
        source_commit=COMMIT,
        baseline_commit=COMMIT,
    )


def replace_once(path, old, new):
    with open(path, "r") as handle:
        text = handle.read()
    if text.count(old) < 1:
        raise AssertionError("test mutation target mismatch")
    with open(path, "w") as handle:
        handle.write(text.replace(old, new, 1))


def main():
    if len(sys.argv) != 2:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4s9_official_closure.py POLICY")
    policy_path = os.path.abspath(sys.argv[1])
    policy = closure.official.validate_policy(
        closure.p4c.read_json(policy_path))
    policy_hash = closure.p4c.sha256_file(policy_path)
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4s9-closure-")
    try:
        stats = os.path.join(root, "statistics")
        resource = os.path.join(root, "resource")
        write_component(stats, policy_hash, "statistics")
        write_component(resource, policy_hash, "resource")
        parity = []
        for rank in (1, 2, 4):
            path = os.path.join(root, "parity-rank{}.log".format(rank))
            write_parity_log(path, rank)
            parity.append(path)
        ctest = os.path.join(root, "focused-ctest.log")
        write_ctest(ctest, policy)
        output = os.path.join(root, "closure")
        closure.analyze_command(make_args(
            policy_path, stats, resource, parity, ctest, output))
        decision = closure.p4c.read_json(os.path.join(
            output, "p4s9_official_closure_decision.json"))
        if decision["p4s9_decision"] != "GO" or \
                not decision["p4f_authorized"]:
            raise AssertionError("synthetic P4-S9 closure GO mismatch")

        mismatch = os.path.join(root, "parity-rank2-mismatch.log")
        shutil.copyfile(parity[1], mismatch)
        replace_once(mismatch, "SAMPLE sample=1 configuration=52",
                     "SAMPLE sample=1 configuration=99")
        mismatch_output = os.path.join(root, "closure-mismatch")
        closure.analyze_command(make_args(
            policy_path, stats, resource,
            [parity[0], mismatch, parity[2]], ctest, mismatch_output))
        mismatch_decision = closure.p4c.read_json(os.path.join(
            mismatch_output, "p4s9_official_closure_decision.json"))
        if mismatch_decision["p4s9_decision"] != "STOP" or \
                mismatch_decision["p4f_authorized"]:
            raise AssertionError("MPI parity mismatch did not STOP")

        missing_ctest = os.path.join(root, "missing-ctest.log")
        write_ctest(
            missing_ctest, policy,
            omit=policy["focused_tests"]["required_names"][0])
        missing_output = os.path.join(root, "closure-missing-ctest")
        closure.analyze_command(make_args(
            policy_path, stats, resource, parity, missing_ctest,
            missing_output))
        missing_decision = closure.p4c.read_json(os.path.join(
            missing_output, "p4s9_official_closure_decision.json"))
        if missing_decision["p4s9_decision"] != "STOP":
            raise AssertionError("missing focused test did not STOP")

        failed_ctest = os.path.join(root, "failed-ctest.log")
        write_ctest(
            failed_ctest, policy,
            fail=policy["focused_tests"]["required_names"][0])
        failed_output = os.path.join(root, "closure-failed-ctest")
        closure.analyze_command(make_args(
            policy_path, stats, resource, parity, failed_ctest,
            failed_output))
        failed_decision = closure.p4c.read_json(os.path.join(
            failed_output, "p4s9_official_closure_decision.json"))
        if failed_decision["p4s9_decision"] != "STOP":
            raise AssertionError("failed focused test did not STOP")

        unbound_stats = os.path.join(root, "statistics-unbound")
        shutil.copytree(stats, unbound_stats)
        ledger_path = os.path.join(unbound_stats, "checksums.sha256")
        with open(ledger_path, "r") as handle:
            ledger_lines = handle.readlines()
        evidence_name = "p4s9_official_markov_statistics.json"
        ledger_lines = [line for line in ledger_lines
                        if not line.rstrip("\n").endswith(
                            "  " + evidence_name)]
        with open(ledger_path, "w") as handle:
            handle.writelines(ledger_lines)
        expect_failure(lambda: closure.analyze_command(make_args(
            policy_path, unbound_stats, resource, parity, ctest,
            os.path.join(root, "closure-unbound-component"))),
            "component omitted from ledger")

        bad_provenance = make_args(
            policy_path, stats, resource, parity, ctest,
            os.path.join(root, "closure-bad-provenance"))
        bad_provenance.source_commit = "39d48c2"
        expect_failure(lambda: closure.analyze_command(bad_provenance),
                       "source commit provenance")

        replace_once(
            os.path.join(stats, "p4s9_official_markov_statistics.json"),
            '"official_statistics_pass": true',
            '"official_statistics_pass": false')
        expect_failure(lambda: closure.analyze_command(make_args(
            policy_path, stats, resource, parity, ctest,
            os.path.join(root, "closure-corrupt"))), "ledger corruption")
        print("bounded Krylov P4-S9 official closure pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
