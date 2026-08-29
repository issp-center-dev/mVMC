#!/usr/bin/env python3

import argparse
import json
import os
import pathlib
import platform
import re
import shutil
import subprocess
import sys
import tempfile

import compare_power_lanczos_p4f_evidence as comparator
import generate_power_lanczos_p4f_v3_confirmatory_evidence as generator
import prepare_power_lanczos_p4f_v3_secondary as preparer
import validate_power_lanczos_p4f_v3_result as validator


GENERATOR_PATH = pathlib.Path(generator.__file__).resolve()
GENERATOR_SHA256 = (
    "789a40d7e1512c9727c84de85417e4f1f3d0bd9508754c5499d2d6b41b4230fe"
)
PREPARER_SHA256 = (
    "fc84c9f02aeae494426295c5cde1125c158d7d48637b7f444d13ef5824328dbd"
)
COMPARATOR_SHA256 = (
    "4b8e1576157eb2ac6ecb163254ecc2fac4c2d50feace954649a0df37ed8b5efd"
)
HEX64 = re.compile(r"[0-9a-f]{64}")
PREPARED_PATHS = {
    "raw_traces": "raw_traces",
    "primary_evidence": "primary/p4f_v3_confirmatory_evidence.json",
    "primary_decision": "primary/p4f_v3_confirmatory_decision.json",
    "primary_metadata": "primary/metadata.json",
    "solver_policy": "frozen_inputs/solver_policy.json",
    "confirmatory_policy": "frozen_inputs/confirmatory_policy.json",
}


def require_regular(path, label, executable=False):
    path = pathlib.Path(path)
    if not path.is_file() or path.is_symlink():
        raise AssertionError(f"{label} must be a regular file")
    if executable and not os.access(path, os.X_OK):
        raise AssertionError(f"{label} must be executable")
    return path


def verify_tool(path, expected_sha, label, executable=False):
    path = require_regular(path, label, executable=executable)
    actual = validator.sha256_path(path)
    if actual != expected_sha:
        raise AssertionError(f"unexpected {label} SHA-256")
    return actual


def tree_file_digests(root, excluded):
    root = pathlib.Path(root)
    digests = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise AssertionError(
                f"symlink in evidence tree: {path.relative_to(root)}")
        if path.is_dir():
            continue
        if not path.is_file():
            raise AssertionError(
                f"non-regular entry in evidence tree: {path.relative_to(root)}")
        relative = path.relative_to(root).as_posix()
        if relative == excluded:
            continue
        digests[relative] = validator.sha256_path(path)
    return digests


def verify_tree_manifest(root, manifest_name):
    root = pathlib.Path(root)
    if not root.is_dir() or root.is_symlink():
        raise AssertionError("evidence root must be a regular directory")
    manifest = require_regular(root / manifest_name, "evidence manifest")
    recorded = validator.parse_checksum_lines(
        manifest.read_bytes(), manifest_name)
    actual = tree_file_digests(root, manifest_name)
    if recorded != actual:
        raise AssertionError("evidence manifest file census/checksum mismatch")
    return {
        "manifest_sha256": validator.sha256_path(manifest),
        "entry_count": len(actual),
    }


def read_prepared(prepared_dir):
    prepared_dir = pathlib.Path(prepared_dir)
    manifest = verify_tree_manifest(
        prepared_dir, "extraction_manifest.sha256")
    report_path = require_regular(
        prepared_dir / "preparation.json", "preparation report")
    report = validator.strict_json(
        report_path.read_bytes(), "secondary preparation report")
    if report.get("schema") != 1 or report.get("valid") is not True or \
            report.get("trace_count") != 24 or \
            report.get("manifest") != "extraction_manifest.sha256" or \
            report.get("manifest_entry_count") != manifest["entry_count"] or \
            report.get("prepared_paths") != PREPARED_PATHS or \
            report.get("validator_sha256") != preparer.VALIDATOR_SHA256:
        raise AssertionError("secondary preparation report contract mismatch")
    for field in ("archive_sha256", "solver_policy_sha256",
                  "confirmatory_policy_sha256"):
        if not isinstance(report.get(field), str) or \
                HEX64.fullmatch(report[field]) is None:
            raise AssertionError(f"invalid preparation hash: {field}")
    raw_dir = prepared_dir / PREPARED_PATHS["raw_traces"]
    raw_paths = sorted(raw_dir.glob("*.log"))
    if not raw_dir.is_dir() or raw_dir.is_symlink() or \
            len(raw_paths) != 24 or any(
                not path.is_file() or path.is_symlink()
                for path in raw_paths):
        raise AssertionError("prepared raw trace census mismatch")
    for key, relative in PREPARED_PATHS.items():
        if key == "raw_traces":
            continue
        require_regular(prepared_dir / relative, f"prepared {key}")
    solver_sha = validator.sha256_path(
        prepared_dir / PREPARED_PATHS["solver_policy"])
    confirmatory_sha = validator.sha256_path(
        prepared_dir / PREPARED_PATHS["confirmatory_policy"])
    if solver_sha != report["solver_policy_sha256"] or \
            confirmatory_sha != report["confirmatory_policy_sha256"]:
        raise AssertionError("prepared policy SHA-256 mismatch")
    return report, manifest


def verify_secondary_tree(root, preparation, expected_driver_sha):
    root = pathlib.Path(root)
    manifest = verify_tree_manifest(root, "checksums.sha256")
    if manifest["entry_count"] != 29:
        raise AssertionError("secondary evidence file census mismatch")
    evidence_path = require_regular(
        root / "p4f_v3_confirmatory_evidence.json", "secondary evidence")
    decision_path = require_regular(
        root / "p4f_v3_confirmatory_decision.json", "secondary decision")
    metadata_path = require_regular(root / "metadata.json", "secondary metadata")
    evidence = validator.strict_json(
        evidence_path.read_bytes(), "secondary evidence")
    decision = validator.strict_json(
        decision_path.read_bytes(), "secondary decision")
    metadata = validator.strict_json(
        metadata_path.read_bytes(), "secondary metadata")
    if evidence.get("metadata") != metadata or \
            evidence.get("decision") not in ("GO", "STOP") or \
            decision.get("decision") != evidence["decision"] or \
            evidence.get("trace_count") != 24 or \
            decision.get("trace_count") != 24:
        raise AssertionError("secondary evidence/decision contract mismatch")
    if metadata.get("generator_sha256") != GENERATOR_SHA256 or \
            metadata.get("reparse_driver_sha256") != expected_driver_sha:
        raise AssertionError("secondary executable provenance mismatch")
    if validator.sha256_path(root / "solver_policy.json") != \
            preparation["solver_policy_sha256"] or \
            validator.sha256_path(root / "confirmatory_policy.json") != \
            preparation["confirmatory_policy_sha256"]:
        raise AssertionError("secondary frozen policy mismatch")
    cases = evidence.get("cases")
    if not isinstance(cases, list) or len(cases) != 24:
        raise AssertionError("secondary evidence case census mismatch")
    reparse_paths = set()
    for case in cases:
        relative = validator.normalize_path(case.get("reparse_output"))
        if not relative.startswith("raw_reparse/") or \
                not relative.endswith(".reparse.log") or \
                relative in reparse_paths:
            raise AssertionError("secondary reparse path mismatch")
        reparse_paths.add(relative)
        path = require_regular(root / relative, "secondary reparse output")
        if validator.sha256_path(path) != case.get("reparse_output_sha256"):
            raise AssertionError("secondary reparse SHA binding mismatch")
    if len(reparse_paths) != 24:
        raise AssertionError("secondary reparse output census mismatch")
    return {
        "evidence": evidence,
        "decision": decision,
        "metadata": metadata,
        "manifest": manifest,
        "evidence_sha256": validator.sha256_path(evidence_path),
        "decision_sha256": validator.sha256_path(decision_path),
    }


def write_closure_manifest(root):
    root = pathlib.Path(root)
    relative = "closure_manifest.sha256"
    digests = tree_file_digests(root, relative)
    manifest = root / relative
    manifest.write_text("".join(
        f"{digest}  {path}\n" for path, digest in sorted(digests.items())
    ), encoding="ascii")
    return validator.sha256_path(manifest), len(digests)


def run_secondary_closure(prepared_dir, reparse_driver, output_dir,
                          expected_reparse_driver_sha256):
    prepared_dir = pathlib.Path(prepared_dir)
    reparse_driver = pathlib.Path(reparse_driver)
    if prepared_dir.is_symlink() or reparse_driver.is_symlink():
        raise AssertionError("secondary closure inputs must not be symlinks")
    prepared_dir = prepared_dir.resolve()
    reparse_driver = reparse_driver.resolve()
    output_dir = pathlib.Path(output_dir)
    if HEX64.fullmatch(expected_reparse_driver_sha256 or "") is None:
        raise AssertionError("invalid expected reparse driver SHA-256")
    if output_dir.exists() or output_dir.is_symlink():
        raise AssertionError("secondary closure output already exists")
    parent = output_dir.parent
    if not parent.is_dir() or parent.is_symlink():
        raise AssertionError("secondary closure parent must be a directory")

    verify_tool(validator.__file__, preparer.VALIDATOR_SHA256,
                "result validator")
    runner_path = require_regular(__file__, "secondary closure runner")
    runner_sha = validator.sha256_path(runner_path)
    verify_tool(preparer.__file__, PREPARER_SHA256, "secondary preparer")
    verify_tool(comparator.__file__, COMPARATOR_SHA256, "evidence comparator")
    verify_tool(GENERATOR_PATH, GENERATOR_SHA256, "secondary generator")
    driver_sha = verify_tool(
        reparse_driver, expected_reparse_driver_sha256,
        "reparse driver", executable=True)
    preparation, prepared_manifest = read_prepared(prepared_dir)
    if preparation["solver_policy_sha256"] != \
            generator.SOLVER_POLICY_SHA256 or \
            preparation["confirmatory_policy_sha256"] != \
            generator.CONFIRMATORY_POLICY_SHA256:
        raise AssertionError("generator/frozen policy binding mismatch")

    staging = pathlib.Path(tempfile.mkdtemp(
        prefix=f".{output_dir.name}.closure.", dir=parent))
    try:
        workflow = staging / "workflow"
        secondary_dir = staging / "secondary"
        workflow.mkdir()
        command = [
            sys.executable, str(GENERATOR_PATH),
            "--trace-dir", str(prepared_dir / PREPARED_PATHS["raw_traces"]),
            "--reparse-driver", str(reparse_driver),
            "--solver-policy",
            str(prepared_dir / PREPARED_PATHS["solver_policy"]),
            "--confirmatory-policy",
            str(prepared_dir / PREPARED_PATHS["confirmatory_policy"]),
            "--output-dir", str(secondary_dir),
        ]
        environment = os.environ.copy()
        environment.update({
            "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
        })
        run = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env=environment, check=False, timeout=4800)
        (workflow / "secondary_generator.stdout.log").write_bytes(run.stdout)
        (workflow / "secondary_generator.stderr.log").write_bytes(run.stderr)
        if run.returncode not in (0, 1):
            raise AssertionError(
                f"secondary generator infrastructure failure: {run.returncode}")

        secondary = verify_secondary_tree(
            secondary_dir, preparation, driver_sha)
        expected_returncode = 0 if secondary["decision"]["decision"] == \
            "GO" else 1
        if run.returncode != expected_returncode:
            raise AssertionError("secondary generator returncode mismatch")

        primary_evidence_path = prepared_dir / \
            PREPARED_PATHS["primary_evidence"]
        primary_decision_path = prepared_dir / \
            PREPARED_PATHS["primary_decision"]
        primary_evidence = comparator.strict_json_text(
            primary_evidence_path.read_text(encoding="utf-8"),
            "primary evidence")
        primary_decision = comparator.strict_json_text(
            primary_decision_path.read_text(encoding="utf-8"),
            "primary decision")
        if primary_decision.get("decision") != \
                preparation["primary_decision"]:
            raise AssertionError("prepared primary decision mismatch")
        evidence_comparison = comparator.compare_documents(
            primary_evidence, secondary["evidence"])
        decision_comparison = comparator.compare_documents(
            primary_decision, secondary["decision"])
        comparison = {
            "schema": 1,
            "evidence": evidence_comparison,
            "decision": decision_comparison,
            "comparison": "PASS" if (
                evidence_comparison["comparison"] == "PASS" and
                decision_comparison["comparison"] == "PASS"
            ) else "FAIL",
        }
        (staging / "primary_secondary_comparison.json").write_text(
            json.dumps(comparison, indent=2, allow_nan=False) + "\n",
            encoding="utf-8")

        primary_value = primary_decision.get("decision")
        secondary_value = secondary["decision"].get("decision")
        if comparison["comparison"] == "PASS" and \
                primary_value == secondary_value == "GO":
            outcome = "GO"
        elif comparison["comparison"] == "PASS" and \
                primary_value == secondary_value == "STOP":
            outcome = "STOP"
        else:
            outcome = "MISMATCH"
        closure = {
            "schema": 1,
            "artifact": "power_lanczos_p4f_v3_secondary_closure",
            "archive_sha256": preparation["archive_sha256"],
            "prepared_manifest_sha256":
                prepared_manifest["manifest_sha256"],
            "validator_sha256": preparer.VALIDATOR_SHA256,
            "runner_sha256": runner_sha,
            "preparer_sha256": PREPARER_SHA256,
            "generator_sha256": GENERATOR_SHA256,
            "comparator_sha256": COMPARATOR_SHA256,
            "reparse_driver_sha256": driver_sha,
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "thread_environment": {
                key: environment[key] for key in (
                    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                )
            },
            "primary_evidence_sha256":
                validator.sha256_path(primary_evidence_path),
            "secondary_evidence_sha256": secondary["evidence_sha256"],
            "secondary_manifest_sha256":
                secondary["manifest"]["manifest_sha256"],
            "primary_decision": primary_value,
            "secondary_decision": secondary_value,
            "comparison": comparison["comparison"],
            "outcome": outcome,
            "p4f_candidate_for_independent_review": outcome == "GO",
            "p5_authorized": False,
        }
        closure_path = staging / "p4f_v3_secondary_closure.json"
        closure_path.write_text(
            json.dumps(closure, indent=2, allow_nan=False) + "\n",
            encoding="utf-8")

        verify_tool(GENERATOR_PATH, GENERATOR_SHA256, "secondary generator")
        verify_tool(runner_path, runner_sha, "secondary closure runner")
        verify_tool(comparator.__file__, COMPARATOR_SHA256,
                    "evidence comparator")
        verify_tool(reparse_driver, driver_sha, "reparse driver",
                    executable=True)
        final_preparation, final_manifest = read_prepared(prepared_dir)
        if final_preparation != preparation or final_manifest != \
                prepared_manifest:
            raise AssertionError("prepared inputs changed during closure")
        closure_entry_count = len(tree_file_digests(
            staging, "closure_manifest.sha256"))
        closure["closure_manifest_entry_count"] = closure_entry_count
        closure_path.write_text(
            json.dumps(closure, indent=2, allow_nan=False) + "\n",
            encoding="utf-8")
        closure_manifest_sha, final_entry_count = write_closure_manifest(
            staging)
        if final_entry_count != closure_entry_count:
            raise AssertionError("closure manifest entry count changed")
        staging.rename(output_dir)
        result = dict(closure)
        result["closure_manifest_sha256"] = closure_manifest_sha
        return result
    except BaseException:
        if staging.is_dir() and not staging.is_symlink():
            shutil.rmtree(staging)
        raise


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared-dir", required=True, type=pathlib.Path)
    parser.add_argument("--reparse-driver", required=True, type=pathlib.Path)
    parser.add_argument("--expected-reparse-driver-sha256", required=True)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()
    result = run_secondary_closure(
        args.prepared_dir, args.reparse_driver, args.output_dir,
        args.expected_reparse_driver_sha256)
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return 0 if result["outcome"] in ("GO", "STOP") else 1


if __name__ == "__main__":
    raise SystemExit(main())
