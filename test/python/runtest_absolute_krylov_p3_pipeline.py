from __future__ import print_function

import argparse
import copy
import os
import shutil
import subprocess
import sys
import tempfile

import generate_absolute_krylov_p3_evidence as evidence


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative gate was accepted: {}".format(label))


def git_head(repository):
    process = subprocess.run(
        ["git", "-C", repository, "rev-parse", "HEAD"],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True, check=False,
    )
    if process.returncode != 0:
        raise AssertionError("cannot resolve test repository HEAD")
    return process.stdout.strip()


def metadata(commit, fixture_hash, policy_hash):
    return {
        "schema_version": evidence.SCHEMA_VERSION,
        "generated_at_jst": evidence.generated_at_jst(),
        "source_commit": commit,
        "baseline_develop_commit": commit,
        "compiler": "synthetic-compiler",
        "mpi": "synthetic-mpi",
        "blas": "synthetic-blas",
        "strict_fp": "enabled",
        "omp_threads": 1,
        "blas_threads": 1,
        "fixture_sha256": fixture_hash,
        "measurement_policy_sha256": policy_hash,
    }


def write_profile(path, site, qp, ranks, sample_count):
    sector_size = evidence.fixed_sector_size(site)
    scale = 1.0e-6 * (2.0 ** (site / 2.0))
    global_factorizations = sample_count * qp
    with open(path, "w") as handle:
        handle.write(
            "PROFILE schema=1 fixture={} site_count={} qp_total={} "
            "sample_count={} sector_size={} rank_count={} audit=1\n".format(
                evidence.EXPECTED_PROFILE_FIXTURE, site, qp, sample_count,
                sector_size, ranks))
        for sample, sector_index, configuration in \
                evidence.expected_profile_rows(site, sample_count):
            handle.write(
                "ROW {} {} {} 1 0 2 0 3 0 4 0\n".format(
                    sample, sector_index, configuration))
        for scope in ("rank_sum", "rank_max"):
            multiplier = ranks if scope == "rank_sum" else 1
            roots = sample_count * multiplier
            elapsed = sample_count * scale * multiplier
            handle.write(
                "RESOURCE scope={} roots={} recursion=1,1,1,1 "
                "unique=1,1,1,1 memo_hits=1,1,1,1 "
                "memo_misses=1,1,1,1 raw_transitions={} "
                "merged_duplicates=0 cancelled_zero=0 "
                "terminal_requests={} terminal_cache_hits=1 "
                "regular_logical={} near_pivot_logical=0 "
                "singular_logical=0 total_zero_logical=0 "
                "local_factorizations={} "
                "global_factorizations_logical={} frontier_peak=1 "
                "memo_peak=1 krylov_workspace_bytes={} "
                "model_workspace_bytes=1024 amplitude_workspace_bytes=1024 "
                "rss_bytes={} depth_seconds={},{},{},{} "
                "amplitude_seconds={} connectivity_seconds={} "
                "total_seconds={}\n".format(
                    scope, roots, roots, roots, global_factorizations,
                    global_factorizations, global_factorizations,
                    1024 * (2 ** site), 2048 * (2 ** site), elapsed,
                    elapsed, elapsed, elapsed, elapsed, elapsed, elapsed))


def profile_args(output, policy_path, logs, commit):
    return argparse.Namespace(
        source_commit=commit, baseline_commit=commit,
        compiler="synthetic-compiler", mpi="synthetic-mpi",
        blas="synthetic-blas", strict_fp="enabled", omp_threads=1,
        blas_threads=1, policy=policy_path,
        model_manifest=os.path.join(output, "p3_model_manifest.json"),
        integrand_statistics=os.path.join(
            output, "p3_integrand_statistics.json"),
        reference_results=os.path.join(output, "p3_reference_results.json"),
        production_anchor=os.path.join(
            output, "p3_production_anchor.json"),
        profile_log=logs, output=output,
    )


def main():
    if len(sys.argv) != 3:
        raise AssertionError(
            "usage: runtest_absolute_krylov_p3_pipeline.py POLICY REPOSITORY")
    policy_path = os.path.abspath(sys.argv[1])
    repository = os.path.abspath(sys.argv[2])
    policy = evidence.validate_policy_contract(evidence.read_json(policy_path))
    policy_hash = evidence.sha256_file(policy_path)
    fixture_hash = "f" * 64
    commit = git_head(repository)
    root = tempfile.mkdtemp(prefix="absolute-krylov-p3-pipeline-")
    try:
        output = os.path.join(root, "evidence")
        inputs = os.path.join(root, "inputs")
        os.makedirs(output)
        os.makedirs(inputs)
        base = metadata(commit, fixture_hash, policy_hash)
        model = dict(base)
        model.update({
            "artifact": "p3_model_manifest",
            "measurement_policy": policy,
        })
        reference = dict(base)
        reference.update({
            "artifact": "p3_reference_results",
            "deterministic_correctness_gate_pass": True,
            "deterministic_correctness_gate": {
                "status": "PASS", "script_sha256": "1" * 64,
                "stdout_sha256": "2" * 64,
            },
        })
        oracle = dict(base)
        oracle.update({"artifact": "p3_exact_oracle"})
        integrands = dict(base)
        integrands.update({
            "artifact": "p3_integrand_statistics",
            "candidate_summary": {
                "r0_rho_1e-06": {
                    "r": 0, "rho": 1.0e-6,
                    "p3_exact_feasibility_gate_pass": True,
                },
            },
            "cases": {
                "electronic_real": {
                    "guides": {
                        "r0_rho_1e-06": {
                            "entries": {
                                "S_00": {"normalized_second_moment": 1.0},
                            },
                        },
                    },
                },
            },
        })
        anchor_directory = os.path.join(output, "raw_anchor")
        os.makedirs(anchor_directory)
        anchor_runs = []
        for ranks in (1, 2):
            relative = os.path.join(
                "raw_anchor", "anchor_rank{}.log".format(ranks))
            path = os.path.join(output, relative)
            with open(path, "w") as handle:
                handle.write("synthetic production anchor rank {} PASS\n".format(
                    ranks))
            anchor_runs.append({
                "file": relative, "sha256": evidence.sha256_file(path),
                "rank_count": ranks, "status": "PASS",
            })
        anchor = dict(base)
        anchor.update({
            "artifact": "p3_production_anchor", "all_runs_pass": True,
            "runs": anchor_runs,
        })
        for filename, value in (
                ("p3_model_manifest.json", model),
                ("p3_reference_results.json", reference),
                ("p3_exact_oracle.json", oracle),
                ("p3_integrand_statistics.json", integrands),
                ("p3_production_anchor.json", anchor)):
            evidence.write_json(os.path.join(output, filename), value)

        logs = []
        grid = policy["official_scaling_grid"]
        for site in grid["site_counts"]:
            requested = grid["sample_count_by_site"][str(site)]
            sample_count = (evidence.fixed_sector_size(site)
                            if requested == 0 else requested)
            for qp in grid["qp_total"]:
                for ranks in grid["mpi_rank_counts"]:
                    for repeat in range(grid["repeat_count"]):
                        path = os.path.join(
                            inputs,
                            "profile_L{}_qp{}_rank{}_rep{}.log".format(
                                site, qp, ranks, repeat))
                        write_profile(path, site, qp, ranks, sample_count)
                        logs.append(path)

        parsed = evidence.parse_profile_log(logs[0])
        invalid_roots = copy.deepcopy(parsed)
        invalid_roots["resources"]["rank_max"]["roots"] += 1
        expect_failure(
            lambda: evidence.validate_profile_run(invalid_roots, policy),
            "root-count invariant")
        invalid_schedule = copy.deepcopy(parsed)
        invalid_schedule["rows"][0] = invalid_schedule["rows"][1]
        expect_failure(
            lambda: evidence.validate_profile_run(invalid_schedule, policy),
            "sample schedule")

        arguments = profile_args(output, policy_path, logs, commit)
        evidence.profile_command(arguments)
        cost = evidence.read_json(
            os.path.join(output, "p3_cost_feasibility.json"))
        if (cost["p4_decision"] != "GO" or
                not cost["correctness_gate"]["all_pass"]):
            raise AssertionError("synthetic command-level AND gate did not GO")

        bad_reference = copy.deepcopy(reference)
        bad_reference["deterministic_correctness_gate_pass"] = False
        bad_reference_path = os.path.join(root, "bad_reference.json")
        evidence.write_json(bad_reference_path, bad_reference)
        bad_arguments = copy.copy(arguments)
        bad_arguments.reference_results = bad_reference_path
        expect_failure(lambda: evidence.profile_command(bad_arguments),
                       "deterministic correctness AND gate")

        bad_anchor = copy.deepcopy(anchor)
        bad_anchor["all_runs_pass"] = False
        bad_anchor_path = os.path.join(root, "bad_anchor.json")
        evidence.write_json(bad_anchor_path, bad_anchor)
        bad_arguments = copy.copy(arguments)
        bad_arguments.production_anchor = bad_anchor_path
        expect_failure(lambda: evidence.profile_command(bad_arguments),
                       "production anchor AND gate")

        changed_log = os.path.join(root, "changed_rank_trace.log")
        shutil.copyfile(logs[7], changed_log)
        with open(changed_log, "r") as handle:
            changed = handle.read()
        changed = changed.replace(" 1 0 2 0 3 0 4 0\n",
                                  " 1 0 2 0 3 0 5 0\n", 1)
        with open(changed_log, "w") as handle:
            handle.write(changed)
        bad_arguments = copy.copy(arguments)
        bad_arguments.profile_log = list(logs)
        bad_arguments.profile_log[7] = changed_log
        expect_failure(lambda: evidence.profile_command(bad_arguments),
                       "cross-rank trace invariance")

        build_config = os.path.join(root, "build_config.txt")
        review = os.path.join(root, "review.md")
        with open(build_config, "w") as handle:
            handle.write("synthetic build configuration\n")
        with open(review, "w") as handle:
            handle.write("synthetic review\n")
        finalize = argparse.Namespace(
            source_commit=commit, baseline_commit=commit,
            repository=repository, build_config=build_config,
            review=review, output=output,
        )
        invalid_commit = copy.copy(finalize)
        invalid_commit.source_commit = "WORKTREE"
        expect_failure(lambda: evidence.finalize_command(invalid_commit),
                       "noncommit provenance")
        missing = os.path.join(output, anchor_runs[0]["file"])
        held = missing + ".held"
        os.rename(missing, held)
        expect_failure(lambda: evidence.finalize_command(finalize),
                       "missing raw anchor")
        os.rename(held, missing)
        evidence.finalize_command(finalize)
        checksum_path = os.path.join(output, "checksums.sha256")
        with open(checksum_path, "r") as handle:
            checksum_lines = [line.rstrip("\n") for line in handle]
        if len(checksum_lines) < len(logs) + len(anchor_runs):
            raise AssertionError("raw evidence is absent from checksum ledger")
        for line in checksum_lines:
            digest, relative = line.split("  ", 1)
            if evidence.sha256_file(os.path.join(output, relative)) != digest:
                raise AssertionError("final checksum ledger mismatch")
        print("absolute Krylov P3 command pipeline: PASS")
    finally:
        shutil.rmtree(root)


if __name__ == "__main__":
    main()
