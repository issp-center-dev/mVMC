#!/usr/bin/env python3
"""Adjudicate the valid P6-C1 v5 runs without rerunning production processes."""

from __future__ import annotations

import argparse
import copy
import gzip
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from typing import Any, Sequence


SCRIPT = Path(__file__).resolve()
WORKTREE = SCRIPT.parents[2]
PROJECT = WORKTREE.parent
POLICIES = PROJECT / "docs" / "policies"
V5_CHECKER_PATH = (
    WORKTREE / "test" / "python" /
    "power_lanczos_p6c1_cartesian_evidence_v5.py")
_SPEC = importlib.util.spec_from_file_location(
    "power_lanczos_p6c1_cartesian_evidence_v5_parent",
    V5_CHECKER_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError("P6-C1 v5 checker module specification")
_V5 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _V5
_SPEC.loader.exec_module(_V5)


class AdjudicationError(RuntimeError):
    """Raised when a frozen adjudication contract does not hold."""


SHA256 = re.compile(r"[0-9a-f]{64}")
PRE_EXECUTION_MANIFEST = (
    POLICIES /
    "power-lanczos-p6c1-cartesian-pre-execution-manifest-v6.json")
CONTROL_VALIDATOR = (
    POLICIES / "validate_power_lanczos_p6c1_control_v6.py")

ATTEMPT_CLAIM = (
    POLICIES / "power-lanczos-p6c1-control-adjudication-attempt-v6.json")
DEFAULT_EVIDENCE = (
    POLICIES / "power-lanczos-p6c1-control-adjudication-evidence-v6.json")
DEFAULT_MANIFEST = (
    POLICIES /
    "power-lanczos-p6c1-control-adjudication-post-result-manifest-v6.json")
FAILURE_STDOUT = (
    POLICIES /
    "power-lanczos-p6c1-control-adjudication-failure-stdout-v6.txt.gz")
FAILURE_STDERR = (
    POLICIES /
    "power-lanczos-p6c1-control-adjudication-failure-stderr-v6.txt.gz")
FAILURE_MANIFEST = (
    POLICIES /
    "power-lanczos-p6c1-control-adjudication-failure-manifest-v6.json")
FAILURE_CANDIDATE_EVIDENCE = (
    POLICIES /
    "power-lanczos-p6c1-control-adjudication-failure-candidate-evidence-v6.json.gz")
FAILURE_CANDIDATE_MANIFEST = (
    POLICIES /
    "power-lanczos-p6c1-control-adjudication-failure-candidate-manifest-v6.json.gz")

PARENT_V5_PRE_EXECUTION_MANIFEST = (
    POLICIES /
    "power-lanczos-p6c1-cartesian-pre-execution-manifest-v5.json")
PARENT_V5_ATTEMPT = (
    POLICIES / "power-lanczos-p6c1-cartesian-attempt-v5.json")
PARENT_V5_INPUT = (
    POLICIES / "power-lanczos-p6c1-cartesian-input-v5.txt.gz")
PARENT_V5_SERIAL = (
    POLICIES / "power-lanczos-p6c1-cartesian-serial-v5.txt.gz")
PARENT_V5_MPI2 = (
    POLICIES / "power-lanczos-p6c1-cartesian-mpi2-v5.txt.gz")
PARENT_V5_MPI4 = (
    POLICIES / "power-lanczos-p6c1-cartesian-mpi4-v5.txt.gz")
PARENT_V5_CANDIDATE_EVIDENCE = (
    POLICIES /
    "power-lanczos-p6c1-cartesian-failure-candidate-evidence-v5.json.gz")
PARENT_V5_CANDIDATE_MANIFEST = (
    POLICIES /
    "power-lanczos-p6c1-cartesian-failure-candidate-manifest-v5.json.gz")
PARENT_V5_FAILURE_STDOUT = (
    POLICIES / "power-lanczos-p6c1-cartesian-failure-stdout-v5.txt.gz")
PARENT_V5_FAILURE_STDERR = (
    POLICIES / "power-lanczos-p6c1-cartesian-failure-stderr-v5.txt.gz")
PARENT_V5_FAILURE_MANIFEST = (
    POLICIES / "power-lanczos-p6c1-cartesian-failure-manifest-v5.json")
PARENT_V5_FAILURE_REVIEW = (
    PROJECT / "docs" /
    "2026-08-21-power-lanczos-p6c1-v5-failure-review.md")
PARENT_V5_PUBLISHED_SUCCESS_ARTIFACTS = (
    POLICIES / "power-lanczos-p6c1-cartesian-evidence-v5.json",
    POLICIES /
    "power-lanczos-p6c1-cartesian-post-result-manifest-v5.json",
)

PARENT_V5_SHA256 = {
    PARENT_V5_PRE_EXECUTION_MANIFEST:
        "04c70b06dabc5bd752d4b874ca9640b294cdaf8a1d98dc2d374e2fa22bbfd886",
    PARENT_V5_ATTEMPT:
        "545a06f7fde01a5bce908f1b6018b78ea9d2215cabfdf9b6a7f80a2104a08f3f",
    PARENT_V5_INPUT:
        "d0b4ced380cbf536ea55dc692eae157682f690ab914708108c77a688d1c34643",
    PARENT_V5_SERIAL:
        "81a8c27523d1eac65c8d0cb4e93d08d8dd7f6a98493350786412354a15ae7e43",
    PARENT_V5_MPI2:
        "53201d4da319a3df8f1585c0a152ffca4b21f78ec3283ef73e3d6d28c021b586",
    PARENT_V5_MPI4:
        "9289ea5a0543d5fec0ca7f47a28a5d701e9bd11078f5cb50e738540701b98145",
    PARENT_V5_CANDIDATE_EVIDENCE:
        "032d42bab910f87884b7a4b7b6f235f0e6ad1d177fb850136f1a0f61fa539299",
    PARENT_V5_CANDIDATE_MANIFEST:
        "899dfe0347a7fb97563a08c868b30461843e55a47f39909a43f01e3ae226afa6",
    PARENT_V5_FAILURE_STDOUT:
        "9ceffb7310338057cfe71a4ae1e2c98d2c485d81cdef906532a801f457a38d64",
    PARENT_V5_FAILURE_STDERR:
        "7dc60fba2ffd52e468ca3b20e356c423b924b2637569036235eb6b8c8e93cf0b",
    PARENT_V5_FAILURE_MANIFEST:
        "9e8fceed02a2f2446c9c5038e5d795438c96fe84aa95f9c49bcddf1394fd622f",
    PARENT_V5_FAILURE_REVIEW:
        "b72855691ee7ee520fcd7fbbe2b9fedbac80861c9e9bad0189c058f4d13704c8",
}
PARENT_V5_CANDIDATE_EVIDENCE_RAW_SHA256 = (
    "6c5a43dc9c343e8ba0d8c990565957b0cf88ad8e53e1fbe92596eb1b8fbfbfea")
PARENT_V5_CANDIDATE_MANIFEST_RAW_SHA256 = (
    "450d7d8cadc8ca6220e47897c4fc17a1e6d5f49b4f4a358c97ec2cf8f0188283")
NORMALIZED_PAYLOAD_SHA256 = (
    "e6800b5dd6b21c4f3887b68f34d71eb19a508c6394db9086be3519b50e49eb4d")
ATTEMPT_SUCCESS_NEXT_PHASE = "p6_c1_v6_independent_post_result_review"
ATTEMPT_FAILURE_NEXT_PHASE = "p6_c1_v6_independent_failure_review"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AdjudicationError(message)


def sha256_file(path: Path) -> str:
    return _V5.sha256_file(path)


def sha256_bytes(payload: bytes) -> str:
    return _V5.sha256_bytes(payload)


def canonical_json(value: Any) -> bytes:
    return _V5.canonical_json(value)


def load_canonical_json(path: Path, label: str) -> dict[str, Any]:
    payload = path.read_bytes()
    value = json.loads(payload)
    require(isinstance(value, dict) and payload == canonical_json(value),
            f"{label}: canonical JSON")
    return value


def load_deterministic_gzip(path: Path, label: str) -> bytes:
    return _V5.load_deterministic_gzip(path, label)


def exclusive_write(path: Path, payload: bytes) -> None:
    _V5.exclusive_write(path, payload)


def exclusive_publish(source: Path, destination: Path) -> None:
    require(not destination.exists(),
            f"official artifact already exists: {destination.name}")
    os.link(source, destination)
    require(destination.read_bytes() == source.read_bytes(),
            f"exclusive publication round-trip: {destination.name}")


def official_artifact_paths() -> tuple[Path, ...]:
    return (
        ATTEMPT_CLAIM, DEFAULT_EVIDENCE, DEFAULT_MANIFEST,
        FAILURE_STDOUT, FAILURE_STDERR, FAILURE_MANIFEST,
        FAILURE_CANDIDATE_EVIDENCE, FAILURE_CANDIDATE_MANIFEST,
    )


def parent_member_records() -> list[dict[str, str]]:
    roles = (
        "v5_pre_execution_manifest", "v5_attempt_claim", "v5_input",
        "v5_serial", "v5_mpi2", "v5_mpi4",
        "v5_candidate_evidence", "v5_candidate_manifest",
        "v5_failure_stdout", "v5_failure_stderr", "v5_failure_manifest",
        "v5_failure_review",
    )
    return [
        {"role": role, "path": path.relative_to(PROJECT).as_posix(),
         "sha256": PARENT_V5_SHA256[path]}
        for role, path in zip(roles, PARENT_V5_SHA256)
    ]


def expected_parent_v5_post_manifest(
        evidence: dict[str, Any], evidence_raw_sha256: str) -> dict[str, Any]:
    members = [
        {"role": "attempt_claim", "path": PARENT_V5_ATTEMPT.name,
         "sha256": sha256_file(PARENT_V5_ATTEMPT)},
        {"role": "input", "path": PARENT_V5_INPUT.name,
         "sha256": sha256_file(PARENT_V5_INPUT)},
        {"role": "serial", "path": PARENT_V5_SERIAL.name,
         "sha256": sha256_file(PARENT_V5_SERIAL)},
        {"role": "mpi2", "path": PARENT_V5_MPI2.name,
         "sha256": sha256_file(PARENT_V5_MPI2)},
        {"role": "mpi4", "path": PARENT_V5_MPI4.name,
         "sha256": sha256_file(PARENT_V5_MPI4)},
        {"role": "evidence",
         "path": "power-lanczos-p6c1-cartesian-evidence-v5.json",
         "sha256": evidence_raw_sha256},
    ]
    return {
        "schema_version": 4,
        "manifest_id": "power_lanczos_p6c1_cartesian_post_result_v5",
        "status": "observed_recovery_pass",
        "production_authorized": False,
        "calibration_or_confirmation_result": False,
        "pre_execution_manifest_sha256":
            evidence["pre_execution_manifest_sha256"],
        "attempt_claim_sha256": evidence["attempt_claim_sha256"],
        "pilot_regression_count": 9,
        "holdout_count": 2007,
        "overall_count": 2016,
        "pilot_regression_verdict":
            evidence["summaries"]["pilot_regression"]["verdict"],
        "holdout_verdict": evidence["summaries"]["holdout"]["verdict"],
        "overall_verdict": evidence["overall_verdict"],
        "hash_algorithm": "sha256",
        "digest_algorithm":
            "sha256(canonical_json_ordered_member_role_path_sha256)",
        "members": members,
        "artifact_set_digest": sha256_bytes(canonical_json(members)),
        "recovery_parent_failure_manifest_sha256":
            _V5.PARENT_V3_FAILURE_MANIFEST_SHA256,
        "serial_process_rerun": False,
        "new_process_world_sizes": [2, 4],
        "next_phase": "p6_c1_v5_independent_recovery_review",
    }


def validate_parent_v5() -> tuple[dict[str, Any], dict[str, Any]]:
    require(all(path.is_file() and sha256_file(path) == expected
                for path, expected in PARENT_V5_SHA256.items()),
            "v5 parent file identities")
    require(not any(path.exists()
                    for path in PARENT_V5_PUBLISHED_SUCCESS_ARTIFACTS),
            "v5 published success artifacts remain absent")
    pre = load_canonical_json(
        PARENT_V5_PRE_EXECUTION_MANIFEST, "v5 pre-execution manifest")
    require(
        pre.get("manifest_id") ==
            "power_lanczos_p6c1_cartesian_pre_execution_v5" and
        pre.get("status") ==
            "accepted_v3_serial_observed_recovery_pending_mpi2_mpi4" and
        pre.get("holdout_observed") is True and
        pre.get("production_authorized") is False and
        pre.get("execution", {}).get("serial_process_rerun") is False and
        pre.get("execution", {}).get("new_process_world_sizes") == [2, 4],
        "v5 pre-execution recovery state")
    attempt = load_canonical_json(PARENT_V5_ATTEMPT, "v5 attempt")
    require(
        attempt.get("attempt_id") ==
            "power_lanczos_p6c1_cartesian_attempt_v5" and
        attempt.get("pre_execution_manifest_sha256") ==
            PARENT_V5_SHA256[PARENT_V5_PRE_EXECUTION_MANIFEST] and
        attempt.get("new_process_run_order") == ["mpi2", "mpi4"] and
        attempt.get("serial_process_rerun") is False and
        attempt.get("rerun_authorized") is False,
        "v5 attempt recovery state")
    failure = load_canonical_json(
        PARENT_V5_FAILURE_MANIFEST, "v5 failure manifest")
    require(
        failure.get("manifest_id") ==
            "power_lanczos_p6c1_cartesian_failure_v5" and
        failure.get("status") == "observed_failed_attempt" and
        failure.get("failed_stage") == "post_result_validation" and
        failure.get("failed_world_size") is None and
        failure.get("process_started") is False and
        failure.get("return_code") is None and
        failure.get("completed_runs") == ["serial", "mpi2", "mpi4"] and
        failure.get("rerun_authorized") is False,
        "v5 failure lifecycle")
    require(load_deterministic_gzip(
                PARENT_V5_FAILURE_STDOUT, "v5 failure stdout") == b"" and
            load_deterministic_gzip(
                PARENT_V5_FAILURE_STDERR, "v5 failure stderr") ==
            b"ValidationError: control validator --post-result failed: "
            b"FAIL post-result parent/lifecycle/census binding\n\n",
            "v5 failure streams")

    run_paths = {
        "serial": (1, PARENT_V5_SERIAL),
        "mpi2": (2, PARENT_V5_MPI2),
        "mpi4": (4, PARENT_V5_MPI4),
    }
    evidence = _V5.build_cartesian_evidence(
        PARENT_V5_INPUT, run_paths, attempt["executable_sha256"],
        PARENT_V5_SHA256[PARENT_V5_PRE_EXECUTION_MANIFEST],
        PARENT_V5_SHA256[PARENT_V5_ATTEMPT])
    evidence_raw = load_deterministic_gzip(
        PARENT_V5_CANDIDATE_EVIDENCE, "v5 candidate evidence")
    require(
        sha256_bytes(evidence_raw) ==
            PARENT_V5_CANDIDATE_EVIDENCE_RAW_SHA256 and
        evidence_raw == canonical_json(evidence),
        "v5 candidate evidence independent regeneration")
    expected_manifest = expected_parent_v5_post_manifest(
        evidence, sha256_bytes(evidence_raw))
    manifest_raw = load_deterministic_gzip(
        PARENT_V5_CANDIDATE_MANIFEST, "v5 candidate manifest")
    require(
        sha256_bytes(manifest_raw) ==
            PARENT_V5_CANDIDATE_MANIFEST_RAW_SHA256 and
        manifest_raw == canonical_json(expected_manifest),
        "v5 candidate manifest independent regeneration")
    require(
        evidence["overall_verdict"] == "PASS" and
        evidence["summaries"]["pilot_regression"]["fixture_count"] == 9 and
        evidence["summaries"]["holdout"]["fixture_count"] == 2007 and
        evidence["summaries"]["overall"]["fixture_count"] == 2016 and
        [item["world_size"] for item in evidence["runs"]] == [1, 2, 4] and
        all(item["case_count"] == 2016 and item["verdict"] == "PASS"
            for item in evidence["runs"]) and
        {item["normalized_payload_sha256"]
         for item in evidence["runs"]} == {NORMALIZED_PAYLOAD_SHA256} and
        evidence["recovery_provenance"]["blind_holdout_claim"] is False,
        "v5 numerical result census")
    return evidence, expected_manifest


def build_attempt_claim(pre_execution_manifest_sha256: str) -> dict[str, Any]:
    require(SHA256.fullmatch(pre_execution_manifest_sha256) is not None,
            "pre-execution manifest SHA grammar")
    return {
        "schema_version": 1,
        "attempt_id": "power_lanczos_p6c1_control_adjudication_attempt_v6",
        "status": "control_adjudication_claimed_publication_pending",
        "production_authorized": False,
        "calibration_or_confirmation_result": False,
        "pre_execution_manifest_sha256": pre_execution_manifest_sha256,
        "parent_v5_failure_manifest_sha256":
            PARENT_V5_SHA256[PARENT_V5_FAILURE_MANIFEST],
        "parent_v5_candidate_evidence_raw_sha256":
            PARENT_V5_CANDIDATE_EVIDENCE_RAW_SHA256,
        "parent_v5_candidate_manifest_raw_sha256":
            PARENT_V5_CANDIDATE_MANIFEST_RAW_SHA256,
        "process_run_order": [],
        "production_process_rerun": False,
        "serial_process_rerun": False,
        "mpi_process_rerun": False,
        "new_process_world_sizes": [],
        "rerun_authorized": False,
        "success_next_phase": ATTEMPT_SUCCESS_NEXT_PHASE,
        "failure_next_phase": ATTEMPT_FAILURE_NEXT_PHASE,
    }


def validate_attempt_claim(claim: dict[str, Any],
                           expected_pre_execution_sha256: str) -> None:
    require(claim == build_attempt_claim(expected_pre_execution_sha256),
            "v6 attempt claim exact object")


def build_adjudication_evidence(
        parent_evidence: dict[str, Any],
        pre_execution_manifest_sha256: str,
        attempt_claim_sha256: str) -> dict[str, Any]:
    require(
        SHA256.fullmatch(pre_execution_manifest_sha256) is not None and
        SHA256.fullmatch(attempt_claim_sha256) is not None,
        "adjudication evidence SHA grammar")
    runs = [{key: item[key] for key in (
                "name", "world_size", "case_count", "raw_sha256",
                "gzip_sha256", "normalized_payload_sha256", "verdict")}
            for item in parent_evidence["runs"]]
    return {
        "schema_version": 1,
        "evidence_id":
            "power_lanczos_p6c1_control_adjudication_evidence_v6",
        "status": "observed_valid_v5_runs_control_adjudication_pass",
        "production_authorized": False,
        "calibration_or_confirmation_result": False,
        "pre_execution_manifest_sha256": pre_execution_manifest_sha256,
        "attempt_claim_sha256": attempt_claim_sha256,
        "current_observation_state": {
            "holdout_observed": True,
            "all_world_sizes_observed": True,
            "observed_world_sizes": [1, 2, 4],
            "blind_holdout_claim": False,
        },
        "process_provenance": {
            "process_run_order": [],
            "production_process_rerun": False,
            "serial_process_rerun": False,
            "mpi_process_rerun": False,
            "new_process_world_sizes": [],
        },
        "parent_v5": {
            "members": parent_member_records(),
            "candidate_evidence_raw_sha256":
                PARENT_V5_CANDIDATE_EVIDENCE_RAW_SHA256,
            "candidate_manifest_raw_sha256":
                PARENT_V5_CANDIDATE_MANIFEST_RAW_SHA256,
            "completed_runs": ["serial", "mpi2", "mpi4"],
            "failure_stage": "post_result_validation",
            "rerun_authorized": False,
        },
        "numerical_result": {
            "pilot_regression_count": 9,
            "holdout_count": 2007,
            "overall_count": 2016,
            "normalized_payload_sha256": NORMALIZED_PAYLOAD_SHA256,
            "runs": runs,
            "overall_verdict": "PASS",
        },
        "control_adjudication": {
            "v5_rejected_binding":
                "historical holdout_observed is false",
            "registered_current_holdout_observed": True,
            "parent_candidate_evidence_independently_regenerated": True,
            "parent_candidate_manifest_independently_regenerated": True,
            "result_publication_only": True,
        },
        "overall_verdict": "PASS",
    }


def validate_adjudication_evidence(
        actual: dict[str, Any], expected: dict[str, Any]) -> None:
    require(set(actual) == {
        "schema_version", "evidence_id", "status",
        "production_authorized", "calibration_or_confirmation_result",
        "pre_execution_manifest_sha256", "attempt_claim_sha256",
        "current_observation_state", "process_provenance", "parent_v5",
        "numerical_result", "control_adjudication", "overall_verdict",
    }, "v6 adjudication evidence top-level census")
    require(actual == expected,
            "v6 adjudication evidence exact object")
    require(
        actual["current_observation_state"] == {
            "holdout_observed": True,
            "all_world_sizes_observed": True,
            "observed_world_sizes": [1, 2, 4],
            "blind_holdout_claim": False,
        } and
        actual["process_provenance"] == {
            "process_run_order": [],
            "production_process_rerun": False,
            "serial_process_rerun": False,
            "mpi_process_rerun": False,
            "new_process_world_sizes": [],
        } and
        actual["numerical_result"]["overall_verdict"] == "PASS" and
        actual["overall_verdict"] == "PASS",
        "v6 adjudication current-state/process binding")


def build_post_result_manifest(
        evidence: dict[str, Any], attempt_path: Path,
        evidence_path: Path) -> dict[str, Any]:
    members = parent_member_records() + [
        {"role": "v6_attempt_claim", "path": attempt_path.name,
         "sha256": sha256_file(attempt_path)},
        {"role": "v6_adjudication_evidence", "path": evidence_path.name,
         "sha256": sha256_file(evidence_path)},
    ]
    return {
        "schema_version": 1,
        "manifest_id":
            "power_lanczos_p6c1_control_adjudication_post_result_v6",
        "status": "observed_control_adjudication_pass",
        "production_authorized": False,
        "calibration_or_confirmation_result": False,
        "p6_c1_exit_candidate": True,
        "pre_execution_manifest_sha256":
            evidence["pre_execution_manifest_sha256"],
        "attempt_claim_sha256": evidence["attempt_claim_sha256"],
        "holdout_observed": True,
        "all_world_sizes_observed": True,
        "observed_world_sizes": [1, 2, 4],
        "blind_holdout_claim": False,
        "overall_count": 2016,
        "overall_verdict": evidence["overall_verdict"],
        "production_process_rerun": False,
        "serial_process_rerun": False,
        "mpi_process_rerun": False,
        "new_process_world_sizes": [],
        "hash_algorithm": "sha256",
        "digest_algorithm":
            "sha256(canonical_json_ordered_member_role_path_sha256)",
        "members": members,
        "artifact_set_digest": sha256_bytes(canonical_json(members)),
        "next_phase": ATTEMPT_SUCCESS_NEXT_PHASE,
    }


def validate_post_result_manifest(
        actual: dict[str, Any], expected: dict[str, Any],
        attempt: dict[str, Any]) -> None:
    require(set(actual) == {
        "schema_version", "manifest_id", "status",
        "production_authorized", "calibration_or_confirmation_result",
        "p6_c1_exit_candidate", "pre_execution_manifest_sha256",
        "attempt_claim_sha256", "holdout_observed",
        "all_world_sizes_observed", "observed_world_sizes",
        "blind_holdout_claim", "overall_count", "overall_verdict",
        "production_process_rerun", "serial_process_rerun",
        "mpi_process_rerun", "new_process_world_sizes", "hash_algorithm",
        "digest_algorithm", "members", "artifact_set_digest",
        "next_phase",
    }, "v6 post-result manifest top-level census")
    require(actual == expected,
            "v6 post-result manifest exact object")
    require(
        actual["status"] == "observed_control_adjudication_pass" and
        actual["p6_c1_exit_candidate"] is True and
        actual["holdout_observed"] is True and
        actual["all_world_sizes_observed"] is True and
        actual["observed_world_sizes"] == [1, 2, 4] and
        actual["blind_holdout_claim"] is False and
        actual["overall_count"] == 2016 and
        actual["overall_verdict"] == "PASS" and
        actual["production_process_rerun"] is False and
        actual["serial_process_rerun"] is False and
        actual["mpi_process_rerun"] is False and
        actual["new_process_world_sizes"] == [] and
        actual["next_phase"] == attempt["success_next_phase"],
        "v6 post-result current-state/lifecycle binding")


def build_failure_manifest(
        pre_execution_manifest_sha256: str, failed_stage: str,
        completed_steps: Sequence[str], error_type: str,
        failure_stdout_sha256: str,
        failure_stderr_sha256: str) -> dict[str, Any]:
    members = parent_member_records()
    if ATTEMPT_CLAIM.exists():
        members.append({
            "role": "v6_attempt_claim", "path": ATTEMPT_CLAIM.name,
            "sha256": sha256_file(ATTEMPT_CLAIM)})
    if FAILURE_CANDIDATE_EVIDENCE.exists():
        members.append({
            "role": "v6_candidate_evidence",
            "path": FAILURE_CANDIDATE_EVIDENCE.name,
            "sha256": sha256_file(FAILURE_CANDIDATE_EVIDENCE)})
    if FAILURE_CANDIDATE_MANIFEST.exists():
        members.append({
            "role": "v6_candidate_manifest",
            "path": FAILURE_CANDIDATE_MANIFEST.name,
            "sha256": sha256_file(FAILURE_CANDIDATE_MANIFEST)})
    if DEFAULT_EVIDENCE.exists():
        members.append({
            "role": "v6_partially_published_evidence",
            "path": DEFAULT_EVIDENCE.name,
            "sha256": sha256_file(DEFAULT_EVIDENCE)})
    if DEFAULT_MANIFEST.exists():
        members.append({
            "role": "v6_partially_published_manifest",
            "path": DEFAULT_MANIFEST.name,
            "sha256": sha256_file(DEFAULT_MANIFEST)})
    members.extend((
        {"role": "v6_failure_stdout", "path": FAILURE_STDOUT.name,
         "sha256": sha256_file(FAILURE_STDOUT)},
        {"role": "v6_failure_stderr", "path": FAILURE_STDERR.name,
         "sha256": sha256_file(FAILURE_STDERR)},
    ))
    return {
        "schema_version": 1,
        "manifest_id":
            "power_lanczos_p6c1_control_adjudication_failure_v6",
        "status": "observed_failed_control_adjudication",
        "production_authorized": False,
        "calibration_or_confirmation_result": False,
        "pre_execution_manifest_sha256": pre_execution_manifest_sha256,
        "attempt_claim_sha256": (
            sha256_file(ATTEMPT_CLAIM) if ATTEMPT_CLAIM.exists() else None),
        "failed_stage": failed_stage,
        "process_started": False,
        "completed_steps": list(completed_steps),
        "error_type": error_type,
        "observation_state":
            "all_world_sizes_observed_control_adjudication_incomplete",
        "production_process_rerun": False,
        "new_process_world_sizes": [],
        "stored_failure_stdout_sha256": failure_stdout_sha256,
        "stored_failure_stderr_sha256": failure_stderr_sha256,
        "stored_failure_streams_sanitized": True,
        "rerun_authorized": False,
        "hash_algorithm": "sha256",
        "digest_algorithm":
            "sha256(canonical_json_ordered_member_role_path_sha256)",
        "members": members,
        "artifact_set_digest": sha256_bytes(canonical_json(members)),
        "next_phase": ATTEMPT_FAILURE_NEXT_PHASE,
    }


def validate_failure_manifest(actual: dict[str, Any],
                              expected: dict[str, Any],
                              attempt: dict[str, Any] | None) -> None:
    require(actual == expected,
            "v6 failure manifest exact object")
    stages = {
        "parent_validation": [],
        "evidence_generation": ["parent_validation"],
        "post_result_validation":
            ["parent_validation", "evidence_generation"],
        "result_publication": [
            "parent_validation", "evidence_generation",
            "post_result_validation"],
    }
    require(
        actual["failed_stage"] in stages and
        actual["completed_steps"] == stages[actual["failed_stage"]] and
        actual["process_started"] is False and
        actual["observation_state"] ==
            "all_world_sizes_observed_control_adjudication_incomplete" and
        actual["production_process_rerun"] is False and
        actual["new_process_world_sizes"] == [] and
        isinstance(actual["error_type"], str) and
        re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*",
                     actual["error_type"]) is not None and
        actual["rerun_authorized"] is False and
        actual["next_phase"] == (
            attempt["failure_next_phase"] if attempt is not None else
            ATTEMPT_FAILURE_NEXT_PHASE),
        "v6 failure lifecycle/next-phase binding")


def validate_post_result_artifacts(
        expected_pre_execution_manifest_sha256: str,
        attempt_path: Path = ATTEMPT_CLAIM,
        evidence_path: Path = DEFAULT_EVIDENCE,
        manifest_path: Path = DEFAULT_MANIFEST) -> None:
    require(
        SHA256.fullmatch(expected_pre_execution_manifest_sha256) is not None and
        PRE_EXECUTION_MANIFEST.is_file() and
        sha256_file(PRE_EXECUTION_MANIFEST) ==
            expected_pre_execution_manifest_sha256,
        "external pre-execution manifest trust anchor")
    parent_evidence, _ = validate_parent_v5()
    attempt = load_canonical_json(attempt_path, "v6 attempt")
    validate_attempt_claim(attempt, expected_pre_execution_manifest_sha256)
    expected_evidence = build_adjudication_evidence(
        parent_evidence, expected_pre_execution_manifest_sha256,
        sha256_file(attempt_path))
    actual_evidence = load_canonical_json(evidence_path, "v6 evidence")
    validate_adjudication_evidence(actual_evidence, expected_evidence)
    expected_manifest = build_post_result_manifest(
        expected_evidence, attempt_path, evidence_path)
    actual_manifest = load_canonical_json(manifest_path, "v6 result manifest")
    validate_post_result_manifest(actual_manifest, expected_manifest, attempt)


def run_control_validator(
        mode: str, expected_pre_execution_manifest_sha256: str,
        cartesian_executable: Path, preflight_executable: Path,
        gevp_unit: Path, asan_gevp_unit: Path,
        production_binary: Path, mpiexec: Path,
        artifact_paths: dict[str, Path] | None = None) -> None:
    command = [
        sys.executable, str(CONTROL_VALIDATOR), mode,
        "--expected-pre-execution-manifest-sha256",
        expected_pre_execution_manifest_sha256,
        "--cartesian-executable", str(cartesian_executable),
        "--preflight-executable", str(preflight_executable),
        "--gevp-unit", str(gevp_unit),
        "--asan-gevp-unit", str(asan_gevp_unit),
        "--production-binary", str(production_binary),
        "--mpiexec", str(mpiexec),
    ]
    for name, path in (artifact_paths or {}).items():
        command.extend((f"--{name.replace('_', '-')}", str(path)))
    completed = subprocess.run(
        command, cwd=PROJECT, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False)
    require(
        completed.returncode == 0,
        f"control validator {mode} failed: " +
        completed.stderr.decode("utf-8", "replace"))


def execute_adjudication(
        expected_pre_execution_manifest_sha256: str,
        cartesian_executable: Path, preflight_executable: Path,
        gevp_unit: Path, asan_gevp_unit: Path,
        production_binary: Path, mpiexec: Path) -> None:
    require(
        PRE_EXECUTION_MANIFEST.is_file() and
        sha256_file(PRE_EXECUTION_MANIFEST) ==
            expected_pre_execution_manifest_sha256,
        "external pre-execution manifest trust anchor")
    run_control_validator(
        "--pre-exec", expected_pre_execution_manifest_sha256,
        cartesian_executable, preflight_executable, gevp_unit,
        asan_gevp_unit, production_binary, mpiexec)
    existing = [path.name for path in official_artifact_paths()
                if path.exists()]
    require(not existing,
            "v6 adjudication already claimed or artifact exists: " +
            ",".join(existing))
    attempt_bytes = canonical_json(build_attempt_claim(
        expected_pre_execution_manifest_sha256))
    exclusive_write(ATTEMPT_CLAIM, attempt_bytes)
    current_stage = "parent_validation"
    completed_steps: list[str] = []
    evidence_bytes: bytes | None = None
    manifest_bytes: bytes | None = None
    try:
        parent_evidence, _ = validate_parent_v5()
        completed_steps.append("parent_validation")
        current_stage = "evidence_generation"
        evidence = build_adjudication_evidence(
            parent_evidence, expected_pre_execution_manifest_sha256,
            sha256_file(ATTEMPT_CLAIM))
        evidence_bytes = canonical_json(evidence)
        completed_steps.append("evidence_generation")
        current_stage = "post_result_validation"
        scratch_root = PROJECT / "tmp"
        scratch_root.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(
                prefix="p6c1-v6-adjudication-", dir=scratch_root) as directory:
            staged_evidence = Path(directory) / DEFAULT_EVIDENCE.name
            staged_manifest = Path(directory) / DEFAULT_MANIFEST.name
            staged_evidence.write_bytes(evidence_bytes)
            manifest = build_post_result_manifest(
                evidence, ATTEMPT_CLAIM, staged_evidence)
            manifest_bytes = canonical_json(manifest)
            staged_manifest.write_bytes(manifest_bytes)
            run_control_validator(
                "--post-result", expected_pre_execution_manifest_sha256,
                cartesian_executable, preflight_executable, gevp_unit,
                asan_gevp_unit, production_binary, mpiexec,
                {"attempt": ATTEMPT_CLAIM,
                 "evidence": staged_evidence,
                 "manifest": staged_manifest})
            completed_steps.append("post_result_validation")
            current_stage = "result_publication"
            exclusive_publish(staged_evidence, DEFAULT_EVIDENCE)
            exclusive_publish(staged_manifest, DEFAULT_MANIFEST)
    except Exception as error:
        if evidence_bytes is not None:
            exclusive_write(
                FAILURE_CANDIDATE_EVIDENCE,
                _V5.deterministic_gzip(evidence_bytes))
        if manifest_bytes is not None:
            exclusive_write(
                FAILURE_CANDIDATE_MANIFEST,
                _V5.deterministic_gzip(manifest_bytes))
        sanitized = _V5.sanitized_failure_stderr(
            (type(error).__name__ + ": " + str(error) + "\n").encode(),
            cartesian_executable, mpiexec,
            (preflight_executable, gevp_unit, asan_gevp_unit,
             production_binary))
        exclusive_write(FAILURE_STDOUT, _V5.deterministic_gzip(b""))
        exclusive_write(FAILURE_STDERR, _V5.deterministic_gzip(sanitized))
        failure = build_failure_manifest(
            expected_pre_execution_manifest_sha256, current_stage,
            completed_steps, type(error).__name__, sha256_bytes(b""),
            sha256_bytes(sanitized))
        failure_bytes = canonical_json(failure)
        exclusive_write(FAILURE_MANIFEST, failure_bytes)
        run_control_validator(
            "--failure-result", expected_pre_execution_manifest_sha256,
            cartesian_executable, preflight_executable, gevp_unit,
            asan_gevp_unit, production_binary, mpiexec)
        raise AdjudicationError(
            "official one-time control adjudication failed; rerun "
            "forbidden; failure=" + sha256_bytes(failure_bytes)) from error
    print(
        "PASS P6-C1 control adjudication evidence=" +
        sha256_bytes(evidence_bytes or b"") + " manifest=" +
        sha256_bytes(manifest_bytes or b""))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("validate-parent")
    post = subparsers.add_parser("validate-post-result-artifacts")
    post.add_argument("--expected-pre-execution-manifest-sha256",
                      required=True)
    post.add_argument("--attempt", type=Path, default=ATTEMPT_CLAIM)
    post.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    post.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    execute = subparsers.add_parser("execute-adjudication")
    execute.add_argument("--expected-pre-execution-manifest-sha256",
                         required=True)
    execute.add_argument("--cartesian-executable", type=Path, required=True)
    execute.add_argument("--preflight-executable", type=Path, required=True)
    execute.add_argument("--gevp-unit", type=Path, required=True)
    execute.add_argument("--asan-gevp-unit", type=Path, required=True)
    execute.add_argument("--production-binary", type=Path, required=True)
    execute.add_argument("--mpiexec", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        if arguments.command == "validate-parent":
            evidence, manifest = validate_parent_v5()
            print(
                "PASS v5 parent candidate evidence=" +
                sha256_bytes(canonical_json(evidence)) + " manifest=" +
                sha256_bytes(canonical_json(manifest)))
        elif arguments.command == "validate-post-result-artifacts":
            validate_post_result_artifacts(
                arguments.expected_pre_execution_manifest_sha256,
                arguments.attempt.resolve(), arguments.evidence.resolve(),
                arguments.manifest.resolve())
            print("PASS v6 adjudication post-result artifacts")
        elif arguments.command == "execute-adjudication":
            execute_adjudication(
                arguments.expected_pre_execution_manifest_sha256,
                arguments.cartesian_executable.resolve(),
                arguments.preflight_executable.resolve(),
                arguments.gevp_unit.resolve(),
                arguments.asan_gevp_unit.resolve(),
                arguments.production_binary.resolve(),
                arguments.mpiexec.resolve())
    except (OSError, ValueError, KeyError, IndexError, TypeError,
            json.JSONDecodeError, AdjudicationError) as error:
        print(f"FAIL {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
