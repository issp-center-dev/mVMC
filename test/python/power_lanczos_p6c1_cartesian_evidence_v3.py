#!/usr/bin/env python3
"""P6-C1 Cartesian v3 lifecycle wrapper for the current-Mac refreeze."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any, Sequence


SCRIPT = Path(__file__).resolve()
BASE_PATH = SCRIPT.with_name("power_lanczos_p6c1_cartesian_evidence.py")
_SPEC = importlib.util.spec_from_file_location(
    "power_lanczos_p6c1_cartesian_evidence_v2_base", BASE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError("P6-C1 v2 checker module specification")
_BASE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _BASE
_SPEC.loader.exec_module(_BASE)

_POLICIES = _BASE.POLICIES
_BASE.SCRIPT = SCRIPT
_BASE.CHECKER_SOURCE = SCRIPT
_BASE.DEFAULT_INPUT = (
    _POLICIES / "power-lanczos-p6c1-cartesian-input-v3.txt.gz")
_BASE.DEFAULT_SERIAL = (
    _POLICIES / "power-lanczos-p6c1-cartesian-serial-v3.txt.gz")
_BASE.DEFAULT_MPI2 = (
    _POLICIES / "power-lanczos-p6c1-cartesian-mpi2-v3.txt.gz")
_BASE.DEFAULT_MPI4 = (
    _POLICIES / "power-lanczos-p6c1-cartesian-mpi4-v3.txt.gz")
_BASE.DEFAULT_EVIDENCE = (
    _POLICIES / "power-lanczos-p6c1-cartesian-evidence-v3.json")
_BASE.ATTEMPT_CLAIM = (
    _POLICIES / "power-lanczos-p6c1-cartesian-attempt-v3.json")
_BASE.PRE_EXECUTION_MANIFEST = (
    _POLICIES /
    "power-lanczos-p6c1-cartesian-pre-execution-manifest-v3.json")
_BASE.DEFAULT_MANIFEST = (
    _POLICIES /
    "power-lanczos-p6c1-cartesian-post-result-manifest-v3.json")
_BASE.FAILURE_STDOUT = (
    _POLICIES /
    "power-lanczos-p6c1-cartesian-failure-stdout-v3.txt.gz")
_BASE.FAILURE_STDERR = (
    _POLICIES /
    "power-lanczos-p6c1-cartesian-failure-stderr-v3.txt.gz")
_BASE.FAILURE_MANIFEST = (
    _POLICIES /
    "power-lanczos-p6c1-cartesian-failure-manifest-v3.json")
_BASE.FAILURE_CANDIDATE_EVIDENCE = (
    _POLICIES /
    "power-lanczos-p6c1-cartesian-failure-candidate-evidence-v3.json.gz")
_BASE.FAILURE_CANDIDATE_MANIFEST = (
    _POLICIES /
    "power-lanczos-p6c1-cartesian-failure-candidate-manifest-v3.json.gz")
_BASE.CONTROL_VALIDATOR = (
    _POLICIES / "validate_power_lanczos_p6c1_control_v3.py")

_BUILD_ATTEMPT_CLAIM_V2 = _BASE.build_attempt_claim
_BUILD_CARTESIAN_EVIDENCE_V2 = _BASE.build_cartesian_evidence
_POST_RESULT_MANIFEST_V2 = _BASE.post_result_manifest
_BUILD_FAILURE_MANIFEST_V2 = _BASE.build_failure_manifest
_OFFICIAL_ARTIFACT_PATHS_V3 = _BASE.official_artifact_paths
_HISTORICAL_V2_OFFICIAL_ARTIFACTS = tuple(
    _POLICIES / name for name in (
        "power-lanczos-p6c1-cartesian-attempt-v2.json",
        "power-lanczos-p6c1-cartesian-input-v2.txt.gz",
        "power-lanczos-p6c1-cartesian-serial-v2.txt.gz",
        "power-lanczos-p6c1-cartesian-mpi2-v2.txt.gz",
        "power-lanczos-p6c1-cartesian-mpi4-v2.txt.gz",
        "power-lanczos-p6c1-cartesian-evidence-v2.json",
        "power-lanczos-p6c1-cartesian-post-result-manifest-v2.json",
        "power-lanczos-p6c1-cartesian-failure-stdout-v2.txt.gz",
        "power-lanczos-p6c1-cartesian-failure-stderr-v2.txt.gz",
        "power-lanczos-p6c1-cartesian-failure-manifest-v2.json",
        "power-lanczos-p6c1-cartesian-failure-candidate-evidence-v2.json.gz",
        "power-lanczos-p6c1-cartesian-failure-candidate-manifest-v2.json.gz",
    ))


def build_attempt_claim(
        executable_sha256: str,
        pre_execution_manifest_sha256: str,
        mpiexec_path: Path) -> dict[str, Any]:
    result = _BUILD_ATTEMPT_CLAIM_V2(
        executable_sha256, pre_execution_manifest_sha256, mpiexec_path)
    result["attempt_id"] = "power_lanczos_p6c1_cartesian_attempt_v3"
    return result


def build_cartesian_evidence(
        input_path: Path,
        run_paths: dict[str, tuple[int, Path]],
        executable_sha256: str,
        pre_execution_manifest_sha256: str,
        attempt_claim_sha256: str) -> dict[str, Any]:
    result = _BUILD_CARTESIAN_EVIDENCE_V2(
        input_path, run_paths, executable_sha256,
        pre_execution_manifest_sha256, attempt_claim_sha256)
    result["evidence_id"] = "power_lanczos_p6c1_cartesian_evidence_v3"
    return result


def post_result_manifest(
        evidence: dict[str, Any],
        attempt_path: Path,
        input_path: Path,
        run_paths: dict[str, tuple[int, Path]],
        evidence_path: Path) -> dict[str, Any]:
    result = _POST_RESULT_MANIFEST_V2(
        evidence, attempt_path, input_path, run_paths, evidence_path)
    result["manifest_id"] = (
        "power_lanczos_p6c1_cartesian_post_result_v3")
    return result


def build_failure_manifest(
        executable_sha256: str,
        failed_stage: str,
        failed_world_size: int | None,
        process_started: bool,
        return_code: int | None,
        error_type: str,
        completed_runs: Sequence[str],
        run_paths: dict[str, tuple[int, Path]],
        stored_failure_stdout_sha256: str,
        stored_failure_stderr_sha256: str) -> dict[str, Any]:
    result = _BUILD_FAILURE_MANIFEST_V2(
        executable_sha256, failed_stage, failed_world_size,
        process_started, return_code, error_type, completed_runs,
        run_paths, stored_failure_stdout_sha256,
        stored_failure_stderr_sha256)
    result["manifest_id"] = (
        "power_lanczos_p6c1_cartesian_failure_v3")
    return result


def official_artifact_paths() -> tuple[Path, ...]:
    """Reject a v3 attempt when either v2 or v3 was already consumed."""
    return _OFFICIAL_ARTIFACT_PATHS_V3() + \
        _HISTORICAL_V2_OFFICIAL_ARTIFACTS


_BASE.build_attempt_claim = build_attempt_claim
_BASE.build_cartesian_evidence = build_cartesian_evidence
_BASE.post_result_manifest = post_result_manifest
_BASE.build_failure_manifest = build_failure_manifest
_BASE.official_artifact_paths = official_artifact_paths

for _NAME in dir(_BASE):
    if not _NAME.startswith("__"):
        globals()[_NAME] = getattr(_BASE, _NAME)


if __name__ == "__main__":
    raise SystemExit(_BASE.main())
