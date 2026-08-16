from __future__ import print_function

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

import generate_bounded_krylov_p4s_evidence as evidence


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative gate was accepted: {}".format(label))


def git_head(repository):
    fallback = os.environ.get("MVMC_P4S_SOURCE_COMMIT")
    if fallback and re.match(r"^[0-9a-f]{40}$", fallback):
        return fallback
    process = subprocess.run(
        ["git", "-C", repository, "rev-parse", "HEAD"],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True, check=False,
    )
    if process.returncode != 0:
        raise AssertionError("cannot resolve test repository HEAD")
    return process.stdout.strip()


def write_synthetic_markov_log(path, decision="GO", sample_count=512,
                               official_se=0.2, diagnostic_se=0.3,
                               budget=1.0):
    official_block_length = sample_count // 16
    diagnostic_block_length = sample_count // 32
    if official_se < 0.0 or diagnostic_se < 0.0 or budget <= 0.0:
        raise AssertionError("invalid synthetic SE/budget")
    smaller_se = min(official_se, diagnostic_se)
    larger_se = max(official_se, diagnostic_se)
    stability_ratio = 1.0 if larger_se == 0.0 else (
        float("inf") if smaller_se == 0.0 else larger_se / smaller_se)
    stability_pass = int(stability_ratio <= 1.25)
    pathology_pass = int(not (larger_se > 0.0 and smaller_se == 0.0))
    official_se_budget_ratio = official_se / budget
    diagnostic_se_budget_ratio = diagnostic_se / budget
    conservative_se_budget_ratio = larger_se / budget
    budget_pass = int(conservative_se_budget_ratio <= 1.0)
    with open(path, "w") as handle:
        handle.write(
            "MARKOV schema=1 fixture={} site_count=4 qp_total=1 "
            "sample_count={} sector_size=36 rank_count=1 "
            "cache_requested_bytes=0 rho=0.0001 eta=0.125 "
            "seed=5783547389269792305 seed_hex=0x5034535f52305631 "
            "initial_uniform=0.5 initial_sector_index=10 "
            "initial_configuration=51 total_guide_weight=10 "
            "official_block_count=16 official_block_length={} "
            "diagnostic_block_count=32 diagnostic_block_length={} "
            "minimum_leave_one_denominator=1e-12 "
            "minimum_abs_denominator_mean=1e-12 "
            "maximum_denominator_relative_se=1 "
            "maximum_tau_int=16 block_length_multiplier=16 "
            "maximum_block_stability_ratio=1.25 "
            "maximum_conservative_se_budget_ratio=1\n".format(
                evidence.EXPECTED_MARKOV_FIXTURE, sample_count,
                official_block_length, diagnostic_block_length))
        for order in range(4):
            handle.write(
                "SCALE order={} norm={} log_basis_scale=0\n".format(
                    order, 1.0 + order))
        for sample in range(sample_count):
            handle.write(
                "SAMPLE sample={} configuration={} accepted={} "
                "selected_neighbor={} uniform=0.25 "
                "log_acceptance_ratio=0 denominator=1 "
                "zero_target=0 accepted_generation={}\n".format(
                    sample, 51 + sample % 3, sample % 2,
                    sample % 4, sample // 2))
        handle.write(
            "TRACE attempted={0} completed={0} accepted={1} rejected={1} "
            "positive_support={0} support_violation=0 "
            "floor_supported_zero=0 finite_components={2} "
            "exact_zero_components=0 numeric_zero_components=0 "
            "terminal_amplitude_calls={2} min_log_acceptance_ratio=-1 "
            "max_log_acceptance_ratio=1 sum_log_acceptance_ratio=0 "
            "min_proposal_log_guide=0 max_proposal_log_guide=1 "
            "trace_hash=123456\n".format(
                sample_count, sample_count // 2, sample_count * 4))
        handle.write(
            "SUMMARY sample_count={} zero_target_sample_count=0 "
            "denominator_sum={} denominator_mean=1 "
            "denominator_relative_se=0 denominator_stable=1 "
            "effective_sample_fraction=1 zero_target_sample_fraction=0 "
            "minimum_denominator=1 maximum_denominator=1 "
            "denominator_tail_ratio=1 log_contribution_span=2 "
            "hamiltonian_antihermitian_residual=0 hamiltonian_norm=1\n".
            format(sample_count, sample_count))
        handle.write(
            "TAU series=W component=real sample_count={} variance=0.5 "
            "positive_pair_count=1 tau_int=1 block_length={} "
            "block_multiplier={} pass=1 required=1\n".format(
                sample_count, official_block_length, official_block_length))
        for family in evidence.FAMILIES:
            for row, column in evidence.UPPER_ENTRIES:
                name = "{}_{}{}".format(family, row, column)
                handle.write(
                    "TAU series={} component=real sample_count={} "
                    "variance=0.5 positive_pair_count=1 tau_int=1 "
                    "block_length={} block_multiplier={} pass=1 "
                    "required=1\n".format(
                        name, sample_count, official_block_length,
                        official_block_length))
                handle.write(
                    "TAU series={} component=imag sample_count={} "
                    "variance=0 positive_pair_count=0 tau_int=1 "
                    "block_length={} block_multiplier={} pass=1 "
                    "required=0\n".format(
                        name, sample_count, official_block_length,
                        official_block_length))
                handle.write(
                    "ENTRY family={} row={} column={} exact_re=1 exact_im=0 "
                    "estimate_re=1 estimate_im=0 official_se={} "
                    "diagnostic_se={} official_se_real={} "
                    "official_se_imag=0 covariance_real_imag=0 "
                    "stability_ratio={} stability_pass={} "
                    "pathology_pass={} "
                    "denominator_stable=1 unstable_leave_one_blocks=0 "
                    "budget={} se_budget_ratio={} "
                    "diagnostic_se_budget_ratio={} "
                    "conservative_se_budget_ratio={} budget_pass={}\n".format(
                        family, row, column, official_se, diagnostic_se,
                        official_se, stability_ratio, stability_pass,
                        pathology_pass, budget, official_se_budget_ratio,
                        diagnostic_se_budget_ratio,
                        conservative_se_budget_ratio, budget_pass))
        handle.write(
            "MANIFEST rng_algorithm=1 rng_state=2 rng_stream=0 rng_draws={} "
            "policy_hash=3 guide_shape_hash=4 bounded_plan_hash=5 "
            "amplitude_policy_hash=6 current_configuration_hash=7 "
            "current_scale_hash=8 accepted_generation={}\n".format(
                1 + 2 * sample_count, sample_count // 2))
        support_pass = 1
        tau_pass = 1
        denom_pass = 1
        handle.write(
            "DECISION p4s_decision={} support_pass={} tau_pass={} "
            "block_stability_pass={} budget_pass={} denominator_pass={} "
            "block_pathology_pass={} maximum_tau_int=1 "
            "maximum_se_budget_ratio={} "
            "maximum_official_se_budget_ratio={} "
            "maximum_block_stability_ratio={}\n".format(
                decision, support_pass, tau_pass, stability_pass, budget_pass,
                denom_pass, pathology_pass, conservative_se_budget_ratio,
                official_se_budget_ratio, stability_ratio))


def analyze_args(policy_path, output, logs, repository, commit):
    return argparse.Namespace(
        source_commit=commit,
        baseline_commit=commit,
        compiler="synthetic-compiler",
        mpi="synthetic-mpi",
        blas="synthetic-blas",
        strict_fp="enabled",
        omp_threads=1,
        blas_threads=1,
        repository=repository,
        policy=policy_path,
        output=output,
        markov_log=logs,
        allow_smoke_policy=True,
    )


def main():
    if len(sys.argv) != 3:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4s_evidence.py POLICY REPOSITORY")
    policy_path = os.path.abspath(sys.argv[1])
    repository = os.path.abspath(sys.argv[2])
    evidence.validate_p4s_policy(evidence.read_json(policy_path))
    commit = git_head(repository)
    validation_repository = (
        repository if os.path.isdir(os.path.join(repository, ".git"))
        else None)
    root = tempfile.mkdtemp(prefix="bounded-krylov-p4s-evidence-")
    try:
        log_path = os.path.join(root, "markov.log")
        output = os.path.join(root, "evidence")
        write_synthetic_markov_log(log_path)
        evidence.analyze_command(analyze_args(
            policy_path, output, [log_path], validation_repository, commit))
        decision = evidence.read_json(
            os.path.join(output, "p4s_candidate_decision.json"))
        if decision["p4s_decision"] != "GO" or \
                decision["selected_candidate"] != "cache0_r3_rho_1e-04":
            raise AssertionError("synthetic P4-S pipeline did not GO")
        candidate = decision["candidate_summary"]["cache0_r3_rho_1e-04"]
        if candidate["block_stability_pass"] or \
                not candidate["block_pathology_pass"] or \
                not candidate["budget_pass"]:
            raise AssertionError("P4-S2 ratio diagnostic gate mismatch")
        checksum_path = os.path.join(output, "checksums.sha256")
        with open(checksum_path, "r") as handle:
            checksum_lines = [line.rstrip("\n") for line in handle]
        if len(checksum_lines) < 3:
            raise AssertionError("P4-S checksum ledger is incomplete")
        for line in checksum_lines:
            digest, relative = line.split("  ", 1)
            if evidence.sha256_file(os.path.join(output, relative)) != digest:
                raise AssertionError("P4-S checksum ledger mismatch")

        bad_log_path = os.path.join(root, "bad_markov.log")
        bad_output = os.path.join(root, "bad_evidence")
        write_synthetic_markov_log(bad_log_path, decision="STOP")
        expect_failure(
            lambda: evidence.analyze_command(analyze_args(
                policy_path, bad_output, [bad_log_path],
                validation_repository, commit)),
            "decision AND mismatch")
        bad_budget_log_path = os.path.join(root, "bad_budget_markov.log")
        bad_budget_output = os.path.join(root, "bad_budget_evidence")
        write_synthetic_markov_log(
            bad_budget_log_path, decision="GO", official_se=0.2,
            diagnostic_se=1.2, budget=1.0)
        expect_failure(
            lambda: evidence.analyze_command(analyze_args(
                policy_path, bad_budget_output, [bad_budget_log_path],
                validation_repository, commit)),
            "conservative SE budget mismatch")
        print("bounded Krylov P4-S evidence pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
