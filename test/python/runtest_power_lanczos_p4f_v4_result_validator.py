#!/usr/bin/env python3

import json
import pathlib
import tempfile

import runtest_power_lanczos_p4f_v3_result_validator as v3_fixture
import runtest_power_lanczos_p4f_v4_confirmatory_evidence as v4_fixture
import validate_power_lanczos_p4f_v3_result as v3
import validate_power_lanczos_p4f_v4_result as validator


def build_fixture(root, hard_stop=False):
    result, old_expected = v3_fixture.build_fixture(root)
    confirmatory = result / "confirmatory"
    evidence_path = confirmatory / "p4f_v3_confirmatory_evidence.json"
    decision_path = confirmatory / "p4f_v3_confirmatory_decision.json"
    policy_path = confirmatory / "confirmatory_policy.json"
    metadata_path = confirmatory / "metadata.json"
    evidence = json.loads(evidence_path.read_text())
    decision = json.loads(decision_path.read_text())
    policy = json.loads(policy_path.read_text())
    metadata = json.loads(metadata_path.read_text())
    generator_path = pathlib.Path(v4_fixture.evidence.__file__).resolve()
    generator_sha = v3.sha256_path(generator_path)

    policy["policy_id"] = validator.POLICY_ID
    policy_bytes = (json.dumps(policy, indent=2) + "\n").encode()
    policy_path.write_bytes(policy_bytes)
    policy_sha = v3.sha256_bytes(policy_bytes)
    flagged = []
    stopped = []
    raw_dir = result / "raw_traces"
    for index, case in enumerate(evidence["cases"]):
        trace_name = case["trace"]
        diagnostic_flag = index == 0
        trace_hard_pass = not (hard_stop and index == 0)
        raw_bytes = v4_fixture.trace(
            block_diagnostic_pass=not diagnostic_flag,
            hard_pass=trace_hard_pass,
            seed=case["confirmatory_identity"]["seed"],
            site_count=case["site_count"], qp_total=case["qp_total"])
        (raw_dir / trace_name).write_bytes(raw_bytes)
        case["trace_sha256"] = v3.sha256_bytes(raw_bytes)
        parsed = v4_fixture.evidence.parse_raw_trace(raw_bytes)
        identity = case["confirmatory_identity"]
        identity["checks"]["p4s"] = parsed["p4s_hard_pass"]
        identity["passed"] = all(identity["checks"].values())
        identity["block_stability_diagnostic"] = {
            "role": "report_only_not_a_hard_gate",
            "diagnostic_maximum": 1.25,
            "passed": parsed["block_stability_diagnostic_pass"],
            "maximum_ratio": parsed["maximum_block_stability_ratio"],
        }
        case["mandatory_gate_pass"] = identity["passed"]
        case["decision"] = "GO" if identity["passed"] else "STOP"
        if diagnostic_flag:
            flagged.append(trace_name)
        if not identity["passed"]:
            stopped.append(trace_name)

    final = "STOP" if stopped else "GO"
    metadata["confirmatory_policy_sha256"] = policy_sha
    metadata["generator"] = generator_path.name
    metadata["generator_sha256"] = generator_sha
    evidence["schema"] = 1
    evidence["metadata"] = metadata
    evidence["policy_id"] = validator.POLICY_ID
    evidence["block_stability_diagnostic"] = {
        "role": "report_only_not_a_hard_gate",
        "diagnostic_maximum": 1.25,
        "flagged_trace_count": len(flagged),
        "flagged_traces": flagged,
    }
    evidence["stopped_traces"] = stopped
    evidence["decision"] = final
    decision.update({
        "policy_id": validator.POLICY_ID,
        "confirmatory_policy_sha256": policy_sha,
        "mandatory_trace_go_count": 24 - len(stopped),
        "mandatory_trace_stop_count": len(stopped),
        "block_stability_diagnostic_flag_count": len(flagged),
        "flagged_block_stability_diagnostic_traces": flagged,
        "stopped_traces": stopped,
        "decision": final,
    })
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    new_evidence_path = confirmatory / \
        "p4f_v4_confirmatory_evidence.json"
    new_decision_path = confirmatory / \
        "p4f_v4_confirmatory_decision.json"
    v3_fixture.write_json(new_evidence_path, evidence)
    v3_fixture.write_json(new_decision_path, decision)
    evidence_path.unlink()
    decision_path.unlink()
    v3_fixture.checksum_tree(
        confirmatory, confirmatory / "checksums.sha256",
        "checksums.sha256")

    workflow = result / "workflow"
    old_job = workflow / "p4f_v3_genkai_job.sh"
    new_job = workflow / "p4f_v4_genkai_job.sh"
    old_job.rename(new_job)
    source_sha = old_expected["source_archive_sha256"]
    source_name = "mVMC-p4f-v4-source.tar.gz"
    package_entries = {
        "frozen_inputs/power_lanczos_zero_support_p4f_solver_policy.json":
            old_expected["solver_policy_sha256"],
        "frozen_inputs/power_lanczos_zero_support_p4f_v4_confirmatory_policy.json":
            policy_sha,
        source_name: source_sha,
        "workflow/p4f_v4_genkai_job.sh": old_expected["job_script_sha256"],
    }
    package_manifest = "".join(
        f"{value}  {key}\n" for key, value in sorted(package_entries.items()))
    (workflow / "package_manifest.sha256").write_text(
        package_manifest, encoding="ascii")
    (workflow / "source_archive.sha256").write_text(
        f"{source_sha}  {source_name}\n", encoding="ascii")
    build_lines = (workflow / "build_config.txt").read_text().splitlines()
    build_lines = [
        (f"confirmatory_policy_sha256={policy_sha}"
         if line.startswith("confirmatory_policy_sha256=") else line)
        for line in build_lines
    ]
    (workflow / "build_config.txt").write_text(
        "\n".join(build_lines) + "\n", encoding="ascii")
    (workflow / "genkai_focused_ctest.log").write_text(
        "100% tests passed, 0 tests failed out of 6\n", encoding="ascii")
    (workflow / "genkai_confirmatory_generator.log").write_text(
        json.dumps(decision, sort_keys=True) + "\n", encoding="utf-8")
    v3_fixture.checksum_tree(
        result, result / "artifact-ledger.sha256",
        "artifact-ledger.sha256")
    expected = {
        "schema": 1,
        "solver_policy_sha256": old_expected["solver_policy_sha256"],
        "confirmatory_policy_sha256": policy_sha,
        "package_manifest_sha256":
            v3.sha256_bytes(package_manifest.encode("ascii")),
        "job_script_sha256": old_expected["job_script_sha256"],
        "generator_sha256": generator_sha,
        "source_archive_sha256": source_sha,
        "source_commit": old_expected["source_commit"],
    }
    return result, expected


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError(f"negative v4 validator fixture passed: {label}")


def refresh_checksums(root):
    v3_fixture.checksum_tree(
        root / "confirmatory", root / "confirmatory/checksums.sha256",
        "checksums.sha256")
    v3_fixture.checksum_tree(
        root, root / "artifact-ledger.sha256", "artifact-ledger.sha256")


def main():
    scratch = pathlib.Path.cwd() / "tmp"
    scratch.mkdir(exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(
                prefix="p4f-v4-validator-", dir=scratch) as temporary:
            temporary = pathlib.Path(temporary)
            root, expected = build_fixture(temporary / "go")
            archive = temporary / "valid-go.tar.gz"
            checksum = v3_fixture.pack(root, archive)
            value = validator.validate_result(archive, checksum, expected)
            if value["decision"] != "GO" or \
                    value["block_stability_diagnostic_flag_count"] != 1:
                raise AssertionError(
                    "valid v4 diagnostic-flag GO fixture failed")

            wrong_generator = temporary / "wrong-generator.py"
            wrong_generator.write_text("# wrong generator\n", encoding="ascii")
            original_generator_file = validator.generator.__file__
            validator.generator.__file__ = str(wrong_generator)
            try:
                expect_failure(
                    lambda: validator.validate_result(
                        archive, checksum, expected),
                    "active generator provenance mismatch")
            finally:
                validator.generator.__file__ = original_generator_file

            stop_root, stop_expected = build_fixture(
                temporary / "stop", hard_stop=True)
            stop_archive = temporary / "valid-stop.tar.gz"
            stop_checksum = v3_fixture.pack(stop_root, stop_archive)
            stop_value = validator.validate_result(
                stop_archive, stop_checksum, stop_expected)
            if stop_value["decision"] != "STOP" or \
                    stop_value["mandatory_trace_stop_count"] != 1:
                raise AssertionError("valid v4 hard STOP fixture failed")

            mismatch_root, mismatch_expected = build_fixture(
                temporary / "mismatch")
            evidence_path = mismatch_root / \
                "confirmatory/p4f_v4_confirmatory_evidence.json"
            evidence = json.loads(evidence_path.read_text())
            evidence["cases"][0]["confirmatory_identity"][
                "block_stability_diagnostic"]["passed"] = True
            v3_fixture.write_json(evidence_path, evidence)
            refresh_checksums(mismatch_root)
            mismatch_archive = temporary / "mismatch.tar.gz"
            mismatch_checksum = v3_fixture.pack(
                mismatch_root, mismatch_archive)
            expect_failure(
                lambda: validator.validate_result(
                    mismatch_archive, mismatch_checksum,
                    mismatch_expected),
                "block diagnostic/raw mismatch")

            bad_expected = dict(expected)
            bad_expected["unknown"] = "0" * 64
            expect_failure(
                lambda: validator.validate_result(
                    archive, checksum, bad_expected),
                "unexpected expected-provenance key")
    finally:
        try:
            scratch.rmdir()
        except OSError:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
