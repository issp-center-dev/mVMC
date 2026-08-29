#!/usr/bin/env python3

import copy
import hashlib
import io
import json
import math
import os
import pathlib
import tarfile
import tempfile

import validate_power_lanczos_p4f_v3_result as validator


def digest(value):
    return hashlib.sha256(value).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def checksum_tree(root, output, exclude):
    entries = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        if relative == exclude:
            continue
        entries.append(f"{validator.sha256_path(path)}  {relative}\n")
    output.write_text("".join(entries), encoding="ascii")


def dimension(value):
    prefixes = [
        {
            "sample_count": sample,
            "energy_se": error,
            "exact_energy_error": 0.0,
            "exact_energy_pass": True,
            "rank_stability_pass": True,
        }
        for sample, error in zip((8192, 16384, 32768),
                                 (4.0, 2.0, 1.0))
    ]
    return {
        "dimension": value,
        "solver_pass": True,
        "cutoff_scan": {"full_pass": True, "exact_pass": True},
        "prefixes": prefixes,
        "coefficient_gate": {"passed": True},
        "exact": {"energy": -1.0,
                  "antihermitian_residual": 0.0},
    }


def build_fixture(root):
    seeds = [f"0x{index:016x}" for index in range(1, 5)]
    cases = [{"site_count": site, "qp_total": qp}
             for site in (4, 6, 8) for qp in (1, 4)]
    solver_policy = {
        "policy_id": "solver-test",
        "exact_gate": {"energy_tolerance_relative": 1.0e-11},
        "hermitianization": {"exact_maximum_relative_residual": 1.0e-12},
    }
    solver_bytes = (json.dumps(solver_policy, indent=2) + "\n").encode()
    solver_sha = digest(solver_bytes)
    confirmatory_policy = {
        "policy_id": "power_lanczos_zero_support_p4f_confirmatory_v3",
        "solver_policy_sha256": solver_sha,
        "confirmatory_dataset": {
            "trace_count": 24,
            "seed_order": seeds,
            "case_order": cases,
        },
    }
    confirmatory_bytes = (
        json.dumps(confirmatory_policy, indent=2) + "\n").encode()
    confirmatory_sha = digest(confirmatory_bytes)
    source_sha = "3" * 64
    job_bytes = b"#!/bin/sh\n"
    job_sha = digest(job_bytes)
    generator_sha = "4" * 64

    result = root / "result"
    confirmatory = result / "confirmatory"
    raw = result / "raw_traces"
    reparse = confirmatory / "raw_reparse"
    workflow = result / "workflow"
    raw.mkdir(parents=True)
    reparse.mkdir(parents=True)
    workflow.mkdir(parents=True)
    (confirmatory / "solver_policy.json").write_bytes(solver_bytes)
    (confirmatory / "confirmatory_policy.json").write_bytes(
        confirmatory_bytes)

    evidence_cases = []
    index = 0
    for seed in seeds:
        for physical in cases:
            name = f"{index:04d}-trace.log"
            raw_bytes = f"raw {seed} {physical}\n".encode()
            reparse_bytes = f"reparse {seed} {physical}\n".encode()
            (raw / name).write_bytes(raw_bytes)
            (reparse / f"{name}.reparse.log").write_bytes(reparse_bytes)
            evidence_cases.append({
                "trace": name,
                "trace_sha256": digest(raw_bytes),
                "reparse_output": f"raw_reparse/{name}.reparse.log",
                "reparse_output_sha256": digest(reparse_bytes),
                "site_count": physical["site_count"],
                "qp_total": physical["qp_total"],
                "sample_count": 32768,
                "confirmatory_identity": {
                    "seed": seed,
                    "site_count": physical["site_count"],
                    "qp_total": physical["qp_total"],
                    "checks": {
                        key: True for key in (
                            "driver_schema", "session_schema", "fixture",
                            "seed", "case", "sample_count", "rank_count",
                            "cache", "rho", "rng_stream",
                            "persistent_session", "proposal", "partition",
                            "p4s",
                        )
                    },
                    "passed": True,
                },
                "reconstruction": {"passed": True},
                "matrix_diagnostics": [
                    {"family": family, "passed": True}
                    for family in ("K", "S", "B")
                ],
                "dimensions": [dimension(2), dimension(3)],
                "mandatory_gate_pass": True,
                "decision": "GO",
            })
            index += 1

    pooled = []
    for physical in cases:
        group = [item for item in evidence_cases
                 if item["site_count"] == physical["site_count"] and
                 item["qp_total"] == physical["qp_total"]]
        for dimension_value in (2, 3):
            value = validator.recompute_pooled(
                group, dimension_value, solver_policy)
            pooled.append({
                "dimension": dimension_value,
                "seed_count": 4,
                "sample_counts": [8192, 16384, 32768],
                "per_seed_diagnostics": [
                    {
                        "seed": item["seed"],
                        "energy_se": item["errors"],
                        "zero_pattern": item["pattern"],
                        "zero_pattern_pass": item["pattern_pass"],
                        "reported_final_to_initial_ratio": 0.25,
                        "reported_theil_sen_slope": -1.0,
                    }
                    for item in value["per_seed"]
                ],
                "pooled_energy_se": value["pooled"],
                "zero_rule": value["zero_rule"],
                "final_to_initial_ratio": value["ratio"],
                "maximum_ratio": 0.8,
                "theil_sen_log_se_log_n_slope": value["slope"],
                "seed_zero_patterns_pass": value["seed_patterns_pass"],
                "convergence_pass": value["convergence"],
                "passed": value["passed"],
                "site_count": physical["site_count"],
                "qp_total": physical["qp_total"],
            })

    metadata = {
        "solver_policy_sha256": solver_sha,
        "confirmatory_policy_sha256": confirmatory_sha,
        "generator_sha256": generator_sha,
    }
    evidence = {
        "metadata": metadata,
        "policy_id": confirmatory_policy["policy_id"],
        "trace_count": 24,
        "identity_grid_pass": True,
        "cases": evidence_cases,
        "pooled_se_convergence_gates": pooled,
        "stopped_traces": [],
        "stopped_pooled_gates": [],
        "decision": "GO",
    }
    decision_value = {
        "schema": 1,
        "policy_id": confirmatory_policy["policy_id"],
        "solver_policy_sha256": solver_sha,
        "confirmatory_policy_sha256": confirmatory_sha,
        "trace_count": 24,
        "mandatory_trace_go_count": 24,
        "mandatory_trace_stop_count": 0,
        "pooled_gate_go_count": 12,
        "pooled_gate_stop_count": 0,
        "stopped_traces": [],
        "stopped_pooled_gates": [],
        "decision": "GO",
    }
    write_json(confirmatory / "metadata.json", metadata)
    write_json(confirmatory / "p4f_v3_confirmatory_evidence.json", evidence)
    write_json(confirmatory / "p4f_v3_confirmatory_decision.json",
               decision_value)
    checksum_tree(confirmatory, confirmatory / "checksums.sha256",
                  "checksums.sha256")

    package_entries = {
        "frozen_inputs/power_lanczos_zero_support_p4f_solver_policy.json":
            solver_sha,
        "frozen_inputs/power_lanczos_zero_support_p4f_v3_confirmatory_policy.json":
            confirmatory_sha,
        "mVMC-p4f-v3-39d48c2-source.tar.gz": source_sha,
        "workflow/p4f_v3_genkai_job.sh": job_sha,
    }
    package_manifest = "".join(
        f"{value}  {key}\n" for key, value in sorted(package_entries.items()))
    (workflow / "package_manifest.sha256").write_text(
        package_manifest, encoding="ascii")
    (workflow / "source_archive.sha256").write_text(
        f"{source_sha}  mVMC-p4f-v3-39d48c2-source.tar.gz\n",
        encoding="ascii")
    (workflow / "p4f_v3_genkai_job.sh").write_bytes(job_bytes)
    (workflow / "build_config.txt").write_text(
        "\n".join((
            "source_commit=" + "5" * 40,
            "baseline_develop_commit=" + "5" * 40,
            "source_archive_sha256=" + source_sha,
            "solver_policy_sha256=" + solver_sha,
            "confirmatory_policy_sha256=" + confirmatory_sha,
            "PJM_MPI_PROC=4",
            "strict_fp=icc-fp-model-precise-no-fast-math",
        )) + "\n", encoding="ascii")
    (workflow / "genkai_focused_ctest.log").write_text(
        "100% tests passed, 0 tests failed out of 5\n", encoding="ascii")
    (workflow / "remote_scratch_inventory.txt").write_text(
        "raw_trace_count=24\n", encoding="ascii")
    (workflow / "genkai_confirmatory_generator.log").write_text(
        json.dumps(decision_value, sort_keys=True) + "\n", encoding="utf-8")
    checksum_tree(result, result / "artifact-ledger.sha256",
                  "artifact-ledger.sha256")

    expected = {
        "solver_policy_sha256": solver_sha,
        "confirmatory_policy_sha256": confirmatory_sha,
        "package_manifest_sha256": digest(package_manifest.encode()),
        "job_script_sha256": job_sha,
        "generator_sha256": generator_sha,
        "source_archive_sha256": source_sha,
        "source_commit": "5" * 40,
    }
    return result, expected


def pack(root, target, add_symlink=False):
    with tarfile.open(target, "w:gz") as archive:
        archive.add(root, arcname=".", recursive=True)
        if add_symlink:
            info = tarfile.TarInfo("./unsafe-link")
            info.type = tarfile.SYMTYPE
            info.linkname = "raw_traces/0000-trace.log"
            archive.addfile(info)
    checksum = target.with_suffix(target.suffix + ".sha256")
    checksum.write_text(
        f"{validator.sha256_path(target)}  {target.name}\n",
        encoding="ascii")
    return checksum


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError(f"negative result validator fixture accepted: {label}")


def convert_to_valid_stop(root):
    confirmatory = root / "confirmatory"
    evidence_path = confirmatory / "p4f_v3_confirmatory_evidence.json"
    decision_path = confirmatory / "p4f_v3_confirmatory_decision.json"
    solver_policy = json.loads(
        (confirmatory / "solver_policy.json").read_text())
    evidence = json.loads(evidence_path.read_text())
    group = [item for item in evidence["cases"]
             if item["site_count"] == 4 and item["qp_total"] == 1]
    for case in group:
        target = next(item for item in case["dimensions"]
                      if item["dimension"] == 2)
        target["prefixes"][-1]["energy_se"] = 4.0
    recomputed = validator.recompute_pooled(group, 2, solver_policy)
    if recomputed["passed"]:
        raise AssertionError("STOP fixture pooled gate must fail")
    recorded = next(item for item in evidence[
        "pooled_se_convergence_gates"]
                    if item["site_count"] == 4 and item["qp_total"] == 1 and
                    item["dimension"] == 2)
    recorded.update({
        "pooled_energy_se": recomputed["pooled"],
        "zero_rule": recomputed["zero_rule"],
        "final_to_initial_ratio": recomputed["ratio"],
        "theil_sen_log_se_log_n_slope": recomputed["slope"],
        "seed_zero_patterns_pass": recomputed["seed_patterns_pass"],
        "convergence_pass": recomputed["convergence"],
        "passed": recomputed["passed"],
        "per_seed_diagnostics": [
            {
                "seed": item["seed"],
                "energy_se": item["errors"],
                "zero_pattern": item["pattern"],
                "zero_pattern_pass": item["pattern_pass"],
                "reported_final_to_initial_ratio":
                    item["errors"][-1] / item["errors"][0],
                "reported_theil_sen_slope":
                    validator.theil_sen(item["errors"]),
            }
            for item in recomputed["per_seed"]
        ],
    })
    stopped = [{"site_count": 4, "qp_total": 1, "dimension": 2}]
    evidence["stopped_pooled_gates"] = stopped
    evidence["decision"] = "STOP"
    decision = json.loads(decision_path.read_text())
    decision.update({
        "pooled_gate_go_count": 11,
        "pooled_gate_stop_count": 1,
        "stopped_pooled_gates": stopped,
        "decision": "STOP",
    })
    write_json(evidence_path, evidence)
    write_json(decision_path, decision)
    (root / "workflow/genkai_confirmatory_generator.log").write_text(
        json.dumps(decision, sort_keys=True) + "\n", encoding="utf-8")
    checksum_tree(confirmatory, confirmatory / "checksums.sha256",
                  "checksums.sha256")
    checksum_tree(root, root / "artifact-ledger.sha256",
                  "artifact-ledger.sha256")


def convert_to_valid_trace_stop(root):
    confirmatory = root / "confirmatory"
    evidence_path = confirmatory / "p4f_v3_confirmatory_evidence.json"
    decision_path = confirmatory / "p4f_v3_confirmatory_decision.json"
    evidence = json.loads(evidence_path.read_text())
    case = evidence["cases"][0]
    case["confirmatory_identity"]["checks"]["p4s"] = False
    case["confirmatory_identity"]["passed"] = False
    case["mandatory_gate_pass"] = False
    case["decision"] = "STOP"
    stopped = [case["trace"]]
    evidence["stopped_traces"] = stopped
    evidence["decision"] = "STOP"
    decision = json.loads(decision_path.read_text())
    decision.update({
        "mandatory_trace_go_count": 23,
        "mandatory_trace_stop_count": 1,
        "stopped_traces": stopped,
        "decision": "STOP",
    })
    write_json(evidence_path, evidence)
    write_json(decision_path, decision)
    (root / "workflow/genkai_confirmatory_generator.log").write_text(
        json.dumps(decision, sort_keys=True) + "\n", encoding="utf-8")
    checksum_tree(confirmatory, confirmatory / "checksums.sha256",
                  "checksums.sha256")
    checksum_tree(root, root / "artifact-ledger.sha256",
                  "artifact-ledger.sha256")


def main():
    scratch = pathlib.Path.cwd() / "tmp"
    scratch.mkdir(exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(
                prefix="p4f-v3-validator-", dir=scratch) as temp:
            temp = pathlib.Path(temp)
            root, expected = build_fixture(temp)
            archive = temp / "valid.tar.gz"
            checksum = pack(root, archive)
            value = validator.validate_result(archive, checksum, expected)
            if value["decision"] != "GO" or value["trace_count"] != 24:
                raise AssertionError("positive result validator fixture failed")

            stop_root, stop_expected = build_fixture(temp / "stop")
            convert_to_valid_stop(stop_root)
            stop_archive = temp / "valid-stop.tar.gz"
            stop_checksum = pack(stop_root, stop_archive)
            stop_value = validator.validate_result(
                stop_archive, stop_checksum, stop_expected)
            if stop_value["decision"] != "STOP" or \
                    stop_value["pooled_gate_stop_count"] != 1:
                raise AssertionError("valid STOP result fixture failed")

            trace_stop_root, trace_stop_expected = build_fixture(
                temp / "trace-stop")
            convert_to_valid_trace_stop(trace_stop_root)
            trace_stop_archive = temp / "valid-trace-stop.tar.gz"
            trace_stop_checksum = pack(trace_stop_root, trace_stop_archive)
            trace_stop_value = validator.validate_result(
                trace_stop_archive, trace_stop_checksum,
                trace_stop_expected)
            if trace_stop_value["decision"] != "STOP" or \
                    trace_stop_value["mandatory_trace_stop_count"] != 1:
                raise AssertionError(
                    "valid mandatory-trace STOP result fixture failed")

            identity_root, identity_expected = build_fixture(
                temp / "identity-aggregate")
            convert_to_valid_trace_stop(identity_root)
            identity_evidence = identity_root / \
                "confirmatory/p4f_v3_confirmatory_evidence.json"
            identity_value = json.loads(identity_evidence.read_text())
            identity_value["cases"][0]["confirmatory_identity"][
                "passed"] = True
            write_json(identity_evidence, identity_value)
            checksum_tree(identity_root / "confirmatory",
                          identity_root / "confirmatory/checksums.sha256",
                          "checksums.sha256")
            checksum_tree(identity_root,
                          identity_root / "artifact-ledger.sha256",
                          "artifact-ledger.sha256")
            identity_archive = temp / "bad-identity-aggregate.tar.gz"
            identity_checksum = pack(identity_root, identity_archive)
            expect_failure(
                lambda: validator.validate_result(
                    identity_archive, identity_checksum, identity_expected),
                "identity aggregate mismatch")

            identity_type_root, identity_type_expected = build_fixture(
                temp / "identity-type")
            identity_type_evidence = identity_type_root / \
                "confirmatory/p4f_v3_confirmatory_evidence.json"
            identity_type_value = json.loads(identity_type_evidence.read_text())
            identity_type_value["cases"][0]["confirmatory_identity"][
                "checks"]["p4s"] = 0
            identity_type_value["cases"][0]["confirmatory_identity"][
                "passed"] = False
            write_json(identity_type_evidence, identity_type_value)
            checksum_tree(identity_type_root / "confirmatory",
                          identity_type_root /
                          "confirmatory/checksums.sha256",
                          "checksums.sha256")
            checksum_tree(identity_type_root,
                          identity_type_root / "artifact-ledger.sha256",
                          "artifact-ledger.sha256")
            identity_type_archive = temp / "bad-identity-type.tar.gz"
            identity_type_checksum = pack(
                identity_type_root, identity_type_archive)
            expect_failure(
                lambda: validator.validate_result(
                    identity_type_archive, identity_type_checksum,
                    identity_type_expected),
                "non-Boolean identity check")

            unsafe = temp / "unsafe.tar.gz"
            unsafe_checksum = pack(root, unsafe, add_symlink=True)
            expect_failure(
                lambda: validator.validate_result(
                    unsafe, unsafe_checksum, expected), "symlink member")

            changing_root, _ = build_fixture(temp / "changing")
            changing_archive = temp / "changing.tar.gz"
            pack(changing_root, changing_archive)
            changing_view = validator.ArchiveView(changing_archive)
            (changing_root / "raw_traces/0000-trace.log").write_text(
                "archive replaced after initial census\n", encoding="ascii")
            pack(changing_root, changing_archive)
            expect_failure(
                lambda: changing_view.read_bytes(
                    "raw_traces/0000-trace.log"),
                "archive member changed after initial census")

            raw_path = root / "raw_traces/0000-trace.log"
            raw_path.write_text(
                "mutated without ledger update\n", encoding="ascii")
            bad_ledger = temp / "bad-ledger.tar.gz"
            bad_ledger_checksum = pack(root, bad_ledger)
            expect_failure(
                lambda: validator.validate_result(
                    bad_ledger, bad_ledger_checksum, expected),
                "ledger mismatch")

            provenance_root, provenance_expected = build_fixture(
                temp / "provenance")
            build_config = provenance_root / "workflow/build_config.txt"
            build_config.write_text("\n".join(
                line for line in build_config.read_text(
                    encoding="ascii").splitlines()
                if not line.startswith("source_commit=")
            ) + "\n", encoding="ascii")
            checksum_tree(provenance_root,
                          provenance_root / "artifact-ledger.sha256",
                          "artifact-ledger.sha256")
            provenance = temp / "provenance.tar.gz"
            provenance_checksum = pack(provenance_root, provenance)
            expect_failure(
                lambda: validator.validate_result(
                    provenance, provenance_checksum, provenance_expected),
                "missing exact build provenance key")

            duplicate_root, duplicate_expected = build_fixture(
                temp / "duplicate-json")
            duplicate_evidence = duplicate_root / \
                "confirmatory/p4f_v3_confirmatory_evidence.json"
            duplicate_text = duplicate_evidence.read_text(encoding="utf-8")
            duplicate_text = duplicate_text.replace(
                "{\n", '{\n  "decision": "STOP",\n', 1)
            duplicate_evidence.write_text(duplicate_text, encoding="utf-8")
            checksum_tree(duplicate_root / "confirmatory",
                          duplicate_root /
                          "confirmatory/checksums.sha256",
                          "checksums.sha256")
            checksum_tree(duplicate_root,
                          duplicate_root / "artifact-ledger.sha256",
                          "artifact-ledger.sha256")
            duplicate = temp / "duplicate-json.tar.gz"
            duplicate_checksum = pack(duplicate_root, duplicate)
            expect_failure(
                lambda: validator.validate_result(
                    duplicate, duplicate_checksum, duplicate_expected),
                "duplicate JSON decision key")

            root, expected = build_fixture(temp / "semantic")
            evidence_path = root / \
                "confirmatory/p4f_v3_confirmatory_evidence.json"
            evidence = json.loads(evidence_path.read_text())
            evidence["pooled_se_convergence_gates"][0][
                "pooled_energy_se"][0] *= 2.0
            write_json(evidence_path, evidence)
            checksum_tree(root / "confirmatory",
                          root / "confirmatory/checksums.sha256",
                          "checksums.sha256")
            checksum_tree(root, root / "artifact-ledger.sha256",
                          "artifact-ledger.sha256")
            semantic = temp / "semantic.tar.gz"
            semantic_checksum = pack(root, semantic)
            expect_failure(
                lambda: validator.validate_result(
                    semantic, semantic_checksum, expected),
                "pooled semantic mismatch")
    finally:
        try:
            scratch.rmdir()
        except OSError:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
