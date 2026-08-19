from __future__ import print_function

import argparse
import os
import shutil
import sys
import tempfile

import generate_bounded_krylov_p4s9_official_statistics_evidence as evidence
import runtest_bounded_krylov_p4s_evidence as p4s_test


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative gate was accepted: {}".format(label))


def write_session_log(path, sample_count=256):
    legacy = path + ".legacy"
    p4s_test.write_synthetic_markov_log(
        legacy, decision="GO", sample_count=sample_count,
        official_se=0.2, diagnostic_se=0.3, budget=1.0)
    with open(legacy, "r") as handle:
        text = handle.read()
    text = text.replace(
        "MARKOV schema=1 fixture=p4s_bounded_markov_real",
        "MARKOV schema=3 fixture=p4s9_long_direct_session_official")
    text = text.replace("cache_requested_bytes=0 rho=0.0001",
                        "cache_requested_bytes=268435456 rho=0.01")
    text = text.replace(
        "seed=5783547389269792305 seed_hex=0x5034535f52305631",
        "seed=5779335727480329777 seed_hex=0x5034533952305631 "
        "rng_stream=0 persistent_session=1 amplitude_generation_hash=123 "
        "global_numerator=0 global_denominator=1")
    text += (
        "SESSION schema=3 amplitude_generation_hash=123 "
        "session_root_evaluations={} cache_reset_count=1 "
        "cache_resident_end=10 cache_entries_peak=10 cache_hits=1 "
        "cache_misses=1 engine_heap_allocations=0 session_end_pass=1\n".
        format(sample_count + 1))
    with open(path, "w") as handle:
        handle.write(text)


def args_for(policy, output, log, producer):
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
        policy=policy,
        long_direct_policy=policy,
        session_resource_policy=policy,
        target_session_policy=policy,
        stage_a_evidence=policy,
        stage_a_decision=policy,
        stage_c_local_decision=policy,
        target_local_evidence=policy,
        target_local_decision=policy,
        producer_binary=producer,
        output=output,
        markov_log=[log],
        allow_smoke_policy=True,
    )


def replace_once(path, old, new):
    with open(path, "r") as handle:
        text = handle.read()
    if text.count(old) != 1:
        raise AssertionError("test mutation target mismatch")
    with open(path, "w") as handle:
        handle.write(text.replace(old, new, 1))


def replace_in_prefixed_line(path, prefix, old, new):
    with open(path, "r") as handle:
        lines = handle.readlines()
    matches = [index for index, line in enumerate(lines)
               if line.startswith(prefix)]
    if len(matches) != 1 or old not in lines[matches[0]]:
        raise AssertionError("test line mutation target mismatch")
    index = matches[0]
    lines[index] = lines[index].replace(old, new, 1)
    with open(path, "w") as handle:
        handle.writelines(lines)


def main():
    if len(sys.argv) != 2:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4s9_official_statistics_"
            "evidence.py POLICY")
    policy = os.path.abspath(sys.argv[1])
    evidence.validate_policy(evidence.p4c.read_json(policy))
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4s9-official-stats-")
    try:
        good = os.path.join(root, "good.log")
        write_session_log(good)
        output = os.path.join(root, "evidence")
        evidence.analyze_command(args_for(policy, output, good, __file__))
        decision = evidence.p4c.read_json(os.path.join(
            output, "p4s9_official_statistics_decision.json"))
        if decision["official_statistics_decision"] != "GO" or \
                decision["p4f_authorized"]:
            raise AssertionError("synthetic P4-S9 official decision mismatch")

        bad_session = os.path.join(root, "bad-session.log")
        shutil.copyfile(good, bad_session)
        replace_once(bad_session, "session_root_evaluations=257",
                     "session_root_evaluations=256")
        expect_failure(lambda: evidence.analyze_command(args_for(
            policy, os.path.join(root, "bad-session-evidence"),
            bad_session, __file__)), "session root census")

        bad_mode = os.path.join(root, "bad-mode.log")
        shutil.copyfile(good, bad_mode)
        replace_once(bad_mode, "persistent_session=1",
                     "persistent_session=0")
        expect_failure(lambda: evidence.analyze_command(args_for(
            policy, os.path.join(root, "bad-mode-evidence"),
            bad_mode, __file__)), "persistent session mode")

        bad_seed = os.path.join(root, "bad-seed.log")
        shutil.copyfile(good, bad_seed)
        replace_once(bad_seed, "seed=5779335727480329777 ",
                     "seed=5779335727480329778 ")
        expect_failure(lambda: evidence.analyze_command(args_for(
            policy, os.path.join(root, "bad-seed-evidence"),
            bad_seed, __file__)), "decimal/hex seed binding")

        bad_tau = os.path.join(root, "bad-tau.log")
        shutil.copyfile(good, bad_tau)
        replace_in_prefixed_line(
            bad_tau, "TAU series=W ", "pass=1", "pass=0")
        expect_failure(lambda: evidence.analyze_command(args_for(
            policy, os.path.join(root, "bad-tau-evidence"),
            bad_tau, __file__)), "raw TAU gate")

        bad_denominator = os.path.join(root, "bad-denominator.log")
        shutil.copyfile(good, bad_denominator)
        replace_in_prefixed_line(
            bad_denominator, "SUMMARY ",
            "denominator_relative_se=0 ",
            "denominator_relative_se=0.5 ")
        expect_failure(lambda: evidence.analyze_command(args_for(
            policy, os.path.join(root, "bad-denominator-evidence"),
            bad_denominator, __file__)), "raw denominator summary")

        bad_provenance = args_for(
            policy, os.path.join(root, "bad-provenance-evidence"),
            good, __file__)
        bad_provenance.source_commit = "39d48c2"
        expect_failure(lambda: evidence.analyze_command(bad_provenance),
                       "source commit provenance")

        duplicate_args = args_for(
            policy, os.path.join(root, "duplicate-evidence"), good, __file__)
        duplicate_args.markov_log = [good, good]
        expect_failure(lambda: evidence.analyze_command(duplicate_args),
                       "duplicate run census")
        print("bounded Krylov P4-S9 official statistics pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
