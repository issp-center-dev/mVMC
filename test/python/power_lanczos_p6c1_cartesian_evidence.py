#!/usr/bin/env python3
"""Prepare and independently validate the frozen P6-C1 exact Cartesian run."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import io
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np


SCRIPT = Path(__file__).resolve()
WORKTREE = SCRIPT.parents[2]
PROJECT = WORKTREE.parent
POLICIES = PROJECT / "docs" / "policies"
FIXTURE_MODULE = POLICIES / "power_lanczos_p6a_fixture_v3.py"
REGISTRY = POLICIES / "power-lanczos-p6a-fixture-registry-v3.json"
FROZEN_MANIFEST = POLICIES / "power-lanczos-p6a-artifact-manifest-v9.json"
NUMERICAL_POLICY = (
    POLICIES / "power-lanczos-p6c1-cartesian-numerical-policy-v2.json")
PILOT_POLICY_SNAPSHOT = (
    POLICIES / "power-lanczos-p6c1-focused-pilot-policy-snapshot-v1.json")
PILOT_INPUT_GZIP = (
    POLICIES / "power-lanczos-p6c1-focused-pilot-input-v1.txt.gz")
PILOT_SERIAL_GZIP = (
    POLICIES / "power-lanczos-p6c1-focused-pilot-serial-v1.txt.gz")
PILOT_MPI2_GZIP = (
    POLICIES / "power-lanczos-p6c1-focused-pilot-mpi2-v1.txt.gz")
PILOT_MPI4_GZIP = (
    POLICIES / "power-lanczos-p6c1-focused-pilot-mpi4-v1.txt.gz")
PILOT_EVIDENCE = (
    POLICIES / "power-lanczos-p6c1-focused-pilot-evidence-v1.json")
DESIGN = PROJECT / "docs" / "2026-08-21-power-lanczos-p6c1-cartesian-execution-design.md"
DEFAULT_INPUT = POLICIES / "power-lanczos-p6c1-cartesian-input-v2.txt.gz"
DEFAULT_SERIAL = POLICIES / "power-lanczos-p6c1-cartesian-serial-v2.txt.gz"
DEFAULT_MPI2 = POLICIES / "power-lanczos-p6c1-cartesian-mpi2-v2.txt.gz"
DEFAULT_MPI4 = POLICIES / "power-lanczos-p6c1-cartesian-mpi4-v2.txt.gz"
DEFAULT_EVIDENCE = POLICIES / "power-lanczos-p6c1-cartesian-evidence-v2.json"
ATTEMPT_CLAIM = (
    POLICIES / "power-lanczos-p6c1-cartesian-attempt-v2.json")
PRE_EXECUTION_MANIFEST = (
    POLICIES / "power-lanczos-p6c1-cartesian-pre-execution-manifest-v2.json")
DEFAULT_MANIFEST = (
    POLICIES / "power-lanczos-p6c1-cartesian-post-result-manifest-v2.json")
FAILURE_STDOUT = (
    POLICIES / "power-lanczos-p6c1-cartesian-failure-stdout-v2.txt.gz")
FAILURE_STDERR = (
    POLICIES / "power-lanczos-p6c1-cartesian-failure-stderr-v2.txt.gz")
FAILURE_MANIFEST = (
    POLICIES / "power-lanczos-p6c1-cartesian-failure-manifest-v2.json")
FAILURE_CANDIDATE_EVIDENCE = (
    POLICIES /
    "power-lanczos-p6c1-cartesian-failure-candidate-evidence-v2.json.gz")
FAILURE_CANDIDATE_MANIFEST = (
    POLICIES /
    "power-lanczos-p6c1-cartesian-failure-candidate-manifest-v2.json.gz")
CONTROL_VALIDATOR = (
    POLICIES / "validate_power_lanczos_p6c1_control_v2.py")
DRIVER_SOURCE = WORKTREE / "test" / "power_lanczos_p6c1_cartesian_driver.c"
GEVP_SOURCE = WORKTREE / "src" / "mVMC" / "power_lanczos_gevp.c"
GEVP_HEADER = WORKTREE / "src" / "mVMC" / "include" / "power_lanczos_gevp.h"
CHECKER_SOURCE = SCRIPT

CASE_COUNT = 2016
TRUSTED_FROZEN_MANIFEST_SHA256 = (
    "78e9d30c6603d1529b46da2e2d763fbe52b74be76528f8a6709e7e5f0882d34d")
TRUSTED_FROZEN_ARTIFACT_SET_DIGEST = (
    "f48ca3e29dc5e7f6edff1d73007a33234bc151a960f55a488eaa79df9afb0c84")
TRUSTED_REGISTRY_SHA256 = (
    "ae1e4eb7a1e0766c5c8727c4ae66d0e056f884dd0017cd02ff99d8176d7488fb")
TRUSTED_FIXTURE_MODULE_SHA256 = (
    "2ad8053037e9c974b90d374ca69324ac317aaad3f5641a05265c2cf257345481")
EXPECTED_NUMERICAL_POLICY_SHA256 = (
    "550e488543a88ba50d91294b66262924a1845a6a91fc5b53351d45acfd256389")
PREFLIGHT_FIXTURE_IDS = (
    "F0001", "F1765", "F0254", "F1095", "F0423",
    "F0670", "F2014", "F1007", "F1679",
)
CUTOFFS = {"S48": float.fromhex("0x1p-48"),
           "S40": float.fromhex("0x1p-40"),
           "S32": float.fromhex("0x1p-32"),
           "S24": float.fromhex("0x1p-24")}
EPS64 = 64.0 * sys.float_info.epsilon
BINARY64_ABSOLUTE_TOLERANCE = 5e-14
BINARY64_RELATIVE_TOLERANCE = 5e-15
REPORT_ABSOLUTE_TOLERANCE = 1e-18
REPORT_RELATIVE_TOLERANCE = 2e-10
STRICT_THREAD_ENVIRONMENT = {
    "BLIS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}
PROFILE_IDS = {
    "none": 0, "gutzwiller": 1, "jastrow": 2,
    "spin_jastrow": 3, "dh2": 4, "dh4": 5, "all": 6,
}
SUPPORT_IDS = {"regular": 0, "zero_support": 1, "nqp_singular": 2}
HEX_FLOAT = re.compile(
    r"(?:0x0p\+0|-?0x1(?:\.[0-9a-f]*[1-9a-f])?p[+-](?:0|[1-9][0-9]*)|"
    r"-?0x0\.[0-9a-f]*[1-9a-f]p-1022)\Z")
ORACLE_BINARY64_DECIMAL = re.compile(
    r"-?[0-9]\.[0-9]{17}e[+-][0-9]{3}\Z")
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
INTEGER = re.compile(r"0|-?[1-9][0-9]*\Z")


class ValidationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False) + "\n").encode()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def deterministic_gzip(value: bytes) -> bytes:
    output = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=output,
                       compresslevel=9, mtime=0) as stream:
        stream.write(value)
    return output.getvalue()


def exclusive_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        written = stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    require(written == len(payload) and path.read_bytes() == payload,
            f"exclusive artifact round-trip {path.name}")


def exclusive_publish(source: Path, destination: Path) -> None:
    require(source.is_file() and not destination.exists(),
            f"exclusive publication precondition {destination.name}")
    os.link(source, destination)
    require(destination.read_bytes() == source.read_bytes(),
            f"exclusive publication round-trip {destination.name}")


def sanitized_failure_stderr(payload: bytes,
                             executable_path: Path,
                             mpiexec_path: Path,
                             additional_paths: Sequence[Path] = ()) -> bytes:
    text = payload.decode("utf-8", "replace")
    replacements = (
        (str(executable_path), "<cartesian-driver>"),
        (str(mpiexec_path), "<mpiexec>"),
        (str(WORKTREE), "<worktree>"),
        (str(PROJECT), "<project>"),
        (str(Path.home()), "<home>"),
        *((str(path), f"<artifact-{index}>")
          for index, path in enumerate(additional_paths)),
    )
    for source, target in replacements:
        text = text.replace(source, target)
    text = re.sub(r"(?<![A-Za-z0-9_])/(?:[^\s:]+)",
                  "<absolute-path>", text)
    return text.encode("utf-8")


def load_fixture_module():
    specification = importlib.util.spec_from_file_location(
        "power_lanczos_p6a_fixture_v3", FIXTURE_MODULE)
    require(specification is not None and specification.loader is not None,
            "fixture module import specification")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def validate_numerical_policy(*, accepted_required: bool) -> dict[str, Any]:
    require(sha256_file(NUMERICAL_POLICY) ==
            EXPECTED_NUMERICAL_POLICY_SHA256,
            "numerical policy exact byte identity")
    policy = json.loads(NUMERICAL_POLICY.read_text())
    require(set(policy) == {
        "schema_version", "policy_id", "status", "production_authorized",
        "P6_calibration_result_observed",
        "observed_P6_calibration_artifacts",
        "P6_confirmation_result_observed",
        "observed_P6_confirmation_artifacts", "next_phase", "scope",
        "historical_v3_artifacts_mutated", "precedence",
        "frozen_v3_binding", "focused_pilot", "holdout",
        "oracle_only_numerical_audit",
        "registered_cutoffs", "vector_matrix_comparison", "gevp", "variance",
        "observable_authority", "scaled_zero_evidence", "review_exit_rule",
    }, "numerical policy top-level census")
    require(policy["schema_version"] == 2 and
            policy["policy_id"] ==
                "power_lanczos_p6c1_exact_cartesian_numerical_policy_v2" and
            policy["production_authorized"] is False and
            policy["P6_calibration_result_observed"] is False and
            policy["observed_P6_calibration_artifacts"] == [] and
            policy["P6_confirmation_result_observed"] is False and
            policy["observed_P6_confirmation_artifacts"] == [] and
            policy["historical_v3_artifacts_mutated"] is False,
            "numerical policy calibration/confirmation boundary")
    if accepted_required:
        require(policy["status"] ==
                    "accepted_pilot_observed_holdout_unobserved" and
                policy["next_phase"] ==
                    "p6_c1_one_time_cartesian_holdout_execution",
                "numerical policy accepted holdout boundary")
    else:
        require((policy["status"], policy["next_phase"]) in {
            ("candidate_pilot_observed_holdout_unobserved",
             "p6_c1_preexecution_independent_review"),
            ("accepted_pilot_observed_holdout_unobserved",
             "p6_c1_one_time_cartesian_holdout_execution"),
        }, "numerical policy lifecycle")
    manifest = json.loads(FROZEN_MANIFEST.read_text())
    binding = policy["frozen_v3_binding"]
    require(set(binding) == {
        "artifact_manifest_sha256", "artifact_set_digest",
        "fixture_registry_sha256", "fixture_module_sha256",
    }, "frozen v3 binding census")
    require(binding == {
        "artifact_manifest_sha256": TRUSTED_FROZEN_MANIFEST_SHA256,
        "artifact_set_digest": TRUSTED_FROZEN_ARTIFACT_SET_DIGEST,
        "fixture_registry_sha256": TRUSTED_REGISTRY_SHA256,
        "fixture_module_sha256": TRUSTED_FIXTURE_MODULE_SHA256,
    } and sha256_file(FROZEN_MANIFEST) ==
            TRUSTED_FROZEN_MANIFEST_SHA256 and
            sha256_file(REGISTRY) == TRUSTED_REGISTRY_SHA256 and
            sha256_file(FIXTURE_MODULE) == TRUSTED_FIXTURE_MODULE_SHA256 and
            manifest["artifact_set_digest"] ==
                TRUSTED_FROZEN_ARTIFACT_SET_DIGEST and
            manifest["status"] == "accepted_result_unobserved" and
            manifest["current_status"] == "accepted_result_unobserved" and
            manifest["production_authorized"] is False and
            manifest["observed_P6_calibration_artifacts"] == [] and
            manifest["next_phase"] == "p6_c1_one_time_cartesian",
            "accepted frozen v3 binding")
    focused = policy["focused_pilot"]
    require(focused["observed"] is True and
            focused["fixture_ids"] == list(PREFLIGHT_FIXTURE_IDS) and
            focused["fixture_count"] == len(PREFLIGHT_FIXTURE_IDS) and
            focused["world_sizes"] == [1, 2, 4] and
            focused["execution_policy_snapshot"]["sha256"] ==
                sha256_file(PILOT_POLICY_SNAPSHOT) and
            focused["normalized_payload_sha256"] ==
                "70987998c3748de43625b6d16be6ad72550ec806ca900d8611c0fef20fa05442" and
            focused["policy_influence"]["final_holdout_threshold_changes"] == [],
            "focused pilot exact boundary")
    require(policy["holdout"] == {
        "definition":
            "all_2016_registered_fixtures_minus_focused_pilot_fixture_ids",
        "fixture_count": CASE_COUNT - len(PREFLIGHT_FIXTURE_IDS),
        "observed": False,
        "observed_artifacts": [],
        "reporting_rule": "The nine pilot fixtures are reported as pilot_regression and the remaining 2007 fixtures as holdout; overall also reports all 2016 without treating the pilot as blind confirmation",
    }, "protected holdout exact boundary")
    for artifact in focused["artifacts"]:
        require(set(artifact) == {"role", "path", "sha256"} and
                SHA256.fullmatch(artifact["sha256"]) is not None and
                sha256_file(POLICIES / artifact["path"]) ==
                    artifact["sha256"],
                "focused pilot artifact binding")
    require(policy["registered_cutoffs"] == [
        {"id": "S48", "value": "0x1p-48"},
        {"id": "S40", "value": "0x1p-40"},
        {"id": "S32", "value": "0x1p-32"},
        {"id": "S24", "value": "0x1p-24"},
    ],
        "registered cutoff identity/order")
    comparison = policy["vector_matrix_comparison"]
    require(comparison["maximum_absolute_error"] == "5e-14" and
            comparison["maximum_scale_relative_error"] == "5e-15",
            "registered vector/matrix comparison")
    gevp = policy["gevp"]
    require(gevp["policy_id"] ==
                "power_lanczos_p6c1_exact_cartesian_normwise_gevp_v2" and
            gevp["policy_version"] == 2 and
            gevp["normalization"] == "alpha^H S alpha=1" and
            gevp["acceptance_metric"] ==
                "normwise_backward_error=norm2(K*alpha-E*S*alpha)/max(1,normF(K)*norm2(alpha),abs(E)*normF(S)*norm2(alpha))" and
            float.fromhex(gevp["maximum_normwise_backward_error"]) == EPS64 and
            gevp["diagnostic_report_comparison"] ==
                "absolute error <=1e-18 or relative error <=2e-10 using max(abs(recomputed),DBL_MIN) as the relative scale" and
            gevp["oracle_only_census_failure_count"] == 0,
            "registered GEVP policy")
    require(policy["observable_authority"]["cross_representation_tolerance"] ==
            "absolute error <=1e-10 or error/max(1,abs(state_value)) <=1e-11" and
            policy["variance"]["finite_and_nonnegative_required"] is True,
            "registered observable authority")
    return policy


def parse_fraction(value: Sequence[int]) -> Fraction:
    require(len(value) == 2 and value[1] != 0, "canonical rational")
    return Fraction(int(value[0]), int(value[1]))


def parse_rational_complex(value: dict[str, Sequence[int]]) -> complex:
    return complex(float(parse_fraction(value["re"])),
                   float(parse_fraction(value["im"])))


def hex_float(value: float) -> str:
    require(math.isfinite(value), "finite input value")
    if value == 0.0:
        return "0x0p+0"
    return value.hex()


def input_payload(*, preflight: bool = False,
                  policy_sha_override: str | None = None
                  ) -> tuple[bytes, list[dict[str, Any]]]:
    if policy_sha_override is None:
        validate_numerical_policy(accepted_required=not preflight)
        policy_sha = sha256_file(NUMERICAL_POLICY)
    else:
        require(preflight and SHA256.fullmatch(policy_sha_override) is not None,
                "preflight policy SHA override")
        policy_sha = policy_sha_override
    fixture = load_fixture_module()
    registry = json.loads(REGISTRY.read_text())
    descriptors = fixture.descriptor_census(registry)
    require(len(descriptors) == CASE_COUNT, "fixture descriptor census")
    if preflight:
        by_id = {descriptor["fixture_id"]: descriptor
                 for descriptor in descriptors}
        descriptors = [by_id[fixture_id]
                       for fixture_id in PREFLIGHT_FIXTURE_IDS]
    frozen_shas = (sha256_file(FROZEN_MANIFEST), sha256_file(REGISTRY),
                   sha256_file(FIXTURE_MODULE), policy_sha)
    header = ("P6C1_CARTESIAN_PREFLIGHT_INPUT_V1" if preflight else
              "P6C1_CARTESIAN_INPUT_V1")
    lines = [f"{header} {len(descriptors)} " + " ".join(frozen_shas)]
    logical_rows: list[dict[str, Any]] = []
    for descriptor in descriptors:
        logical = fixture.build_logical_fixture(descriptor)
        rendered = fixture.rendered_manifest(fixture.render_classic_files(logical))
        logical_sha = sha256_bytes(fixture.canonical_bytes(logical))
        model = 0 if descriptor["model"] == "electronic" else 1
        arithmetic = 0 if descriptor["arithmetic"] == "real" else 1
        support = SUPPORT_IDS[str(descriptor["support_case"])]
        profile = PROFILE_IDS[str(descriptor["correlator_profile"])]
        active = logical["correlator"]["active_parameters"]
        projection = [
            float(parse_fraction(value))
            for family in active for value in family["values"]
        ]
        lines.append(
            f"CASE {descriptor['fixture_id']} {model} {descriptor['order']} "
            f"{arithmetic} {descriptor['nqp_full']} {profile} {support} "
            f"{logical_sha} {rendered['aggregate_sha256']} {len(projection)}")
        pair_tokens = []
        for index, item in enumerate(logical["pair_product"]["parameters"]):
            require(item["index"] == index, "pair parameter order")
            value = parse_rational_complex(item["value"])
            pair_tokens.extend((hex_float(value.real), hex_float(value.imag)))
        lines.append("PAIR " + " ".join(pair_tokens))
        lines.append("QP")
        for q, item in enumerate(logical["projection"]["transforms"]):
            require(item["q"] == q, "QP order")
            value = parse_rational_complex(item["weight"])
            mapping = " ".join(str(int(site)) for site in item["site_mapping"])
            lines.append(f"{q} {hex_float(value.real)} {hex_float(value.imag)} {mapping}")
        lines.append("PROJ" + (" " + " ".join(map(hex_float, projection))
                               if projection else ""))
        proof = logical["support_proof"]
        proof_x = int(proof["x_configuration_bits"]) if proof else 0
        proof_y = int(proof["y_configuration_bits"]) if proof else 0
        lines.extend((f"PROOF {proof_x} {proof_y}", "END"))
        logical_rows.append(logical)
    lines.append("END_INPUT")
    return ("\n".join(lines) + "\n").encode(), logical_rows


def preflight_logical_rows() -> list[dict[str, Any]]:
    fixture = load_fixture_module()
    registry = json.loads(REGISTRY.read_text())
    descriptors = fixture.descriptor_census(registry)
    require(len(descriptors) == CASE_COUNT, "fixture descriptor census")
    by_id = {descriptor["fixture_id"]: descriptor
             for descriptor in descriptors}
    require(set(PREFLIGHT_FIXTURE_IDS) <= set(by_id),
            "preflight fixture IDs registered")
    return [fixture.build_logical_fixture(by_id[fixture_id])
            for fixture_id in PREFLIGHT_FIXTURE_IDS]


def prepare_input(path: Path, *, preflight: bool = False) -> None:
    payload, _ = input_payload(preflight=preflight)
    compressed = deterministic_gzip(payload)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(compressed)
    require(gzip.decompress(path.read_bytes()) == payload,
            "input gzip round-trip")
    require(deterministic_gzip(payload) == compressed,
            "input deterministic gzip")
    print(f"prepared {path} "
          f"cases={len(PREFLIGHT_FIXTURE_IDS) if preflight else CASE_COUNT} "
          f"sha256={sha256_bytes(compressed)}")


def prepare_preflight_input(path: Path) -> None:
    payload, _ = input_payload(preflight=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    require(path.read_bytes() == payload, "preflight input round-trip")
    print(f"prepared {path} cases={len(PREFLIGHT_FIXTURE_IDS)} "
          f"sha256={sha256_bytes(payload)}")


def parse_int(token: str, label: str) -> int:
    require(INTEGER.fullmatch(token) is not None, f"{label}: integer grammar")
    return int(token)


def canonical_c_hex(value: float) -> str:
    require(math.isfinite(value), "canonical C hexadecimal float")
    if value == 0.0:
        return "0x0p+0"
    token = value.hex()
    mantissa, exponent = token.split("p", 1)
    if "." in mantissa:
        mantissa = mantissa.rstrip("0").rstrip(".")
    return mantissa + "p" + exponent


def parse_hex(token: str, label: str) -> float:
    require(HEX_FLOAT.fullmatch(token) is not None,
            f"{label}: canonical hexadecimal float grammar")
    value = float.fromhex(token)
    require(math.isfinite(value), f"{label}: finite")
    require(token == canonical_c_hex(value),
            f"{label}: canonical C %a identity")
    return value


def parse_error_bound(token: str, label: str) -> float:
    if token == "NINF":
        return -math.inf
    return parse_hex(token, label)


def parse_complex_tokens(tokens: Sequence[str], offset: int,
                         label: str) -> complex:
    return complex(parse_hex(tokens[offset], label + ".re"),
                   parse_hex(tokens[offset + 1], label + ".im"))


class Lines:
    def __init__(self, payload: bytes):
        require(payload.endswith(b"\n"), "output final newline")
        require(b"\r" not in payload and b"\x00" not in payload,
                "output byte grammar")
        self.lines = payload.decode("ascii").splitlines()
        require(all(line and line == line.strip() and "  " not in line
                    for line in self.lines), "output whitespace grammar")
        self.index = 0

    def take(self, keyword: str | None = None) -> list[str]:
        require(self.index < len(self.lines), "unexpected output EOF")
        tokens = self.lines[self.index].split(" ")
        self.index += 1
        if keyword is not None:
            require(tokens[0] == keyword,
                    f"output line {self.index}: expected {keyword}")
        return tokens

    def peek(self) -> str:
        require(self.index < len(self.lines), "unexpected output EOF")
        return self.lines[self.index].split(" ", 1)[0]

    def finish(self) -> None:
        require(self.index == len(self.lines), "trailing output records")


@dataclass
class ParsedRun:
    world_size: int
    frozen_shas: tuple[str, str, str, str]
    models: dict[tuple[int, int], dict[str, Any]]
    cases: list[dict[str, Any]]
    negative: tuple[int, int]
    raw_sha256: str


def parse_model(lines: Lines) -> tuple[tuple[int, int], dict[str, Any]]:
    header = lines.take("MODEL_H")
    require(len(header) == 6, "MODEL_H field census")
    model, arithmetic, state_count = map(
        lambda pair: parse_int(pair[0], pair[1]),
        zip(header[1:4], ("model", "arithmetic", "state_count")))
    term_count = parse_int(header[4], "term_count")
    operator_count = parse_int(header[5], "operator_count")
    rows = []
    for index in range(state_count * state_count):
        tokens = lines.take("HROW")
        require(len(tokens) == 5, "HROW field census")
        rows.append((parse_int(tokens[1], "HROW target"),
                     parse_int(tokens[2], "HROW source"),
                     parse_complex_tokens(tokens, 3, f"HROW {index}")))
    footer = lines.take("END_MODEL_H")
    require(len(footer) == 3 and parse_int(footer[1], "model footer") == model
            and parse_int(footer[2], "arithmetic footer") == arithmetic,
            "MODEL_H footer")
    return (model, arithmetic), {
        "state_count": state_count, "term_count": term_count,
        "operator_count": operator_count, "rows": rows,
    }


def parse_pack(lines: Lines, expected_name: str, upper_count: int) -> list[complex]:
    tokens = lines.take("PACK")
    require(len(tokens) == 3 + 2 * upper_count and
            tokens[1] == expected_name and
            parse_int(tokens[2], "packed count") == upper_count,
            f"PACK {expected_name} census")
    return [parse_complex_tokens(tokens, 3 + 2 * index,
                                 f"PACK {expected_name} {index}")
            for index in range(upper_count)]


def parse_case(lines: Lines, ordinal: int, world_size: int,
               expected_fixture_id: str | None = None) -> dict[str, Any]:
    tokens = lines.take("BEGIN")
    expected_id = expected_fixture_id or f"F{ordinal:04d}"
    require(len(tokens) == 15 and tokens[1] == expected_id,
            "BEGIN census/order")
    case: dict[str, Any] = {
        "fixture_id": tokens[1],
        "model": parse_int(tokens[2], "case model"),
        "order": parse_int(tokens[3], "case order"),
        "arithmetic": parse_int(tokens[4], "case arithmetic"),
        "nqp": parse_int(tokens[5], "case NQP"),
        "profile": parse_int(tokens[6], "case profile"),
        "support": parse_int(tokens[7], "case support"),
        "logical_sha": tokens[8], "rendered_sha": tokens[9],
        "generation_hash": parse_int(tokens[10], "generation hash"),
        "plan_hash": parse_int(tokens[11], "plan hash"),
        "state_count": parse_int(tokens[12], "state count"),
        "term_count": parse_int(tokens[13], "term count"),
        "operator_count": parse_int(tokens[14], "operator count"),
    }
    require(SHA256.fullmatch(case["logical_sha"]) is not None and
            SHA256.fullmatch(case["rendered_sha"]) is not None and
            case["generation_hash"] > 0 and case["plan_hash"] > 0,
            "BEGIN hash grammar")
    partition = lines.take("PARTITIONS")
    require(len(partition) == 2 + 3 * world_size and
            parse_int(partition[1], "partition world") == world_size,
            "partition census")
    partitions = []
    for rank in range(world_size):
        offset = 2 + 3 * rank
        partitions.append(tuple(parse_int(partition[offset + index],
                                          "partition")
                                for index in range(3)))
    case["partitions"] = partitions
    pivot = lines.take("PIVOT")
    require(len(pivot) == 3, "pivot census")
    case["pivot"] = (parse_hex(pivot[1], "pivot input"),
                     parse_hex(pivot[2], "pivot effective"))
    proof = lines.take("PROOF")
    require(len(proof) == 3, "proof census")
    case["proof"] = (parse_int(proof[1], "proof x"),
                     parse_int(proof[2], "proof y"))
    vectors = case["order"] + 2
    rows = []
    for row_index in range(case["state_count"]):
        row = lines.take("ROW")
        require(len(row) == 13 + 2 * vectors, "ROW field census")
        vector_count = parse_int(row[12], "ROW vector count")
        require(vector_count == vectors, "ROW vector count")
        rows.append({
            "bits": parse_int(row[1], "ROW bits"),
            "scaled_state": parse_int(row[2], "ROW scaled state"),
            "well": parse_int(row[3], "ROW well pivoted"),
            "near": parse_int(row[4], "ROW near"),
            "exact_zero": parse_int(row[5], "ROW exact zero"),
            "numeric_zero": parse_int(row[6], "ROW numeric zero"),
            "local_factorizations": parse_int(row[7], "ROW local factorizations"),
            "global_factorizations": parse_int(row[8], "ROW global factorizations"),
            "amplitude_error_bound": parse_error_bound(
                row[9], f"ROW {row_index} amplitude error bound"),
            "amplitude": parse_complex_tokens(row, 10, f"ROW {row_index} amplitude"),
            "vectors": [parse_complex_tokens(row, 13 + 2 * index,
                                               f"ROW {row_index} v{index}")
                        for index in range(vectors)],
        })
        vector_state = lines.take("VSTATE")
        require(len(vector_state) == 2 + 2 * vectors and
                parse_int(vector_state[1], "VSTATE count") == vectors,
                "VSTATE census")
        rows[-1]["vector_states"] = [
            {"state": parse_int(vector_state[2 + 2 * index], "VSTATE state"),
             "error_bound": parse_error_bound(
                 vector_state[3 + 2 * index], "VSTATE error bound")}
            for index in range(vectors)]
    case["rows"] = rows
    dimension = case["order"] + 1
    upper_count = dimension * (dimension + 1) // 2
    case["packs"] = {name: parse_pack(lines, name, upper_count)
                     for name in ("S", "KF", "KR", "B")}
    solves = {}
    for cutoff_id in CUTOFFS:
        solve = lines.take("SOLVE")
        require(len(solve) >= 4 and solve[1] == cutoff_id,
                "SOLVE order")
        if case["support"] == 2:
            require(solve == ["SOLVE", cutoff_id, "NQP_SINGULAR",
                              "NOT_EVALUATED"], "singular solve grammar")
            solves[cutoff_id] = {"status": "NQP_SINGULAR"}
            continue
        require(len(solve) == 21 and solve[2] == "OK",
                "positive SOLVE census")
        result = {
            "status": "OK", "dimension": parse_int(solve[3], "solve dimension"),
            "rank": parse_int(solve[4], "solve rank"),
            "discarded": parse_int(solve[5], "solve discarded"),
            "multiplicity": parse_int(solve[6], "solve multiplicity"),
            "pivot": parse_int(solve[7], "solve pivot"),
            "normalization": parse_hex(solve[8], "solve normalization"),
            "energy": parse_hex(solve[9], "solve energy"),
            "m2": parse_hex(solve[10], "solve M2"),
            "variance": parse_hex(solve[11], "solve variance"),
            "packed_energy": parse_hex(solve[12], "packed energy"),
            "packed_m2": parse_hex(solve[13], "packed M2"),
            "raw_matrix_variance": parse_hex(solve[14], "raw matrix variance"),
            "backward_error": parse_hex(solve[15], "solve backward error"),
            "raw_action_residual": parse_hex(solve[16], "solve raw action residual"),
            "root_gap": parse_hex(solve[17], "solve root gap"),
            "condition": parse_hex(solve[18], "solve condition"),
            "variance_clamped": parse_int(solve[19], "variance clamp"),
            "matrix_observables_valid": parse_int(solve[20], "matrix observables valid"),
        }
        coefficient = lines.take("COEFF")
        require(len(coefficient) == 3 + 2 * dimension and
                coefficient[1] == cutoff_id and
                parse_int(coefficient[2], "coefficient dimension") == dimension,
                "COEFF census")
        result["coefficient"] = [
            parse_complex_tokens(coefficient, 3 + 2 * index,
                                 f"COEFF {cutoff_id} {index}")
            for index in range(dimension)]
        states = []
        for state_index in range(case["state_count"]):
            state = lines.take("STATE")
            require(len(state) == 7 and state[1] == cutoff_id,
                    "STATE census/order")
            states.append((parse_int(state[2], "STATE bits"),
                           parse_complex_tokens(state, 3, "STATE psi"),
                           parse_complex_tokens(state, 5, "STATE Hpsi")))
        result["states"] = states
        solves[cutoff_id] = result
    case["solves"] = solves
    footer = lines.take("END_CASE")
    require(footer == ["END_CASE", case["fixture_id"]], "END_CASE")
    return case


def parse_output(payload: bytes, *, preflight: bool = False,
                 expected_policy_sha: str | None = None) -> ParsedRun:
    lines = Lines(payload)
    header = lines.take("P6C1_CARTESIAN_PREFLIGHT_OUTPUT_V1" if preflight
                        else "P6C1_CARTESIAN_OUTPUT_V1")
    require(len(header) == 7, "output header census")
    world_size = parse_int(header[1], "world size")
    expected_ids = (PREFLIGHT_FIXTURE_IDS if preflight else
                    tuple(f"F{ordinal:04d}"
                          for ordinal in range(1, CASE_COUNT + 1)))
    require(world_size in (1, 2, 4) and
            parse_int(header[2], "case count") == len(expected_ids),
            "output header identity")
    frozen_shas = tuple(header[3:7])
    require(all(SHA256.fullmatch(value) is not None for value in frozen_shas),
            "output frozen SHA grammar")
    policy_sha = expected_policy_sha or sha256_file(NUMERICAL_POLICY)
    require(SHA256.fullmatch(policy_sha) is not None,
            "expected policy SHA grammar")
    require(frozen_shas == (
        sha256_file(FROZEN_MANIFEST), sha256_file(REGISTRY),
        sha256_file(FIXTURE_MODULE), policy_sha),
        "output frozen/policy SHA binding")
    models: dict[tuple[int, int], dict[str, Any]] = {}
    cases = []
    for ordinal, expected_id in enumerate(expected_ids, 1):
        while lines.peek() == "MODEL_H":
            key, model = parse_model(lines)
            require(key not in models, "duplicate MODEL_H")
            models[key] = model
        cases.append(parse_case(lines, ordinal, world_size, expected_id))
    require(set(models) == {(0, 0), (0, 1), (1, 0), (1, 1)},
            "MODEL_H census")
    negative = lines.take("NEGATIVE")
    require(len(negative) == 4 and negative[:2] == ["NEGATIVE", "ZERO_PACKED_P4"],
            "negative probe grammar")
    negative_result = (parse_int(negative[2], "negative status"),
                       parse_int(negative[3], "negative valid"))
    require(negative_result == (4, 0),
            "zero-packed P4 negative status/valid")
    footer = lines.take("END_OUTPUT")
    require(footer == ["END_OUTPUT", str(len(expected_ids))], "output footer")
    lines.finish()
    return ParsedRun(world_size, frozen_shas, models, cases,
                     negative_result, sha256_bytes(payload))


def configurations(model: int) -> list[int]:
    result = []
    for bits in range(256):
        up = bits & 15
        down = (bits >> 4) & 15
        if up.bit_count() == 2 and down.bit_count() == 2 and (
                model == 0 or down == ((~up) & 15)):
            result.append(bits)
    return result


def fermion_apply(bits: int, annihilate: int,
                   create: int) -> tuple[int, int] | None:
    if not bits & (1 << annihilate) or bits & (1 << create):
        return None
    sign = -1 if (bits & ((1 << annihilate) - 1)).bit_count() % 2 else 1
    after = bits ^ (1 << annihilate)
    if (after & ((1 << create) - 1)).bit_count() % 2:
        sign = -sign
    return after | (1 << create), sign


def independent_hamiltonian(model: int, arithmetic: int) -> tuple[list[int], np.ndarray]:
    states = configurations(model)
    index = {bits: position for position, bits in enumerate(states)}
    matrix = np.zeros((len(states), len(states)), dtype=np.complex128)
    bonds = ((0, 1), (1, 2), (2, 3), (3, 0))
    if model == 0:
        imaginary = (1 / 8, 1 / 16, -1 / 8, -1 / 16)
        intra = (4.0, 17 / 4, 9 / 2, 19 / 4)
        inter = (1.0, 9 / 8, 5 / 4, 11 / 8)
        hund = (1 / 4, 5 / 16, 3 / 8, 7 / 16)
        for source_position, source in enumerate(states):
            for bond_index, (left, right) in enumerate(bonds):
                forward = 1.0 - (imaginary[bond_index] if arithmetic else 0.0) * 1j
                for spin in (0, 1):
                    left_orbital = left + 4 * spin
                    right_orbital = right + 4 * spin
                    for annihilate, create, coefficient in (
                            (right_orbital, left_orbital, forward),
                            (left_orbital, right_orbital, forward.conjugate())):
                        applied = fermion_apply(source, annihilate, create)
                        if applied is not None:
                            matrix[index[applied[0]], source_position] += coefficient * applied[1]
            up = [(source >> site) & 1 for site in range(4)]
            down = [(source >> (site + 4)) & 1 for site in range(4)]
            diagonal = sum(intra[site] * up[site] * down[site]
                           for site in range(4))
            for bond_index, (left, right) in enumerate(bonds):
                diagonal += inter[bond_index] * (up[left] + down[left]) * (up[right] + down[right])
                diagonal -= hund[bond_index] * (up[left] * up[right] + down[left] * down[right])
            matrix[source_position, source_position] += diagonal
    else:
        jxy = (1.0, 9 / 8, 5 / 4, 11 / 8)
        jz = (3 / 4, 13 / 16, 7 / 8, 15 / 16)
        for source_position, source in enumerate(states):
            up = [(source >> site) & 1 for site in range(4)]
            down = [(source >> (site + 4)) & 1 for site in range(4)]
            for bond_index, (left, right) in enumerate(bonds):
                matrix[source_position, source_position] += (
                    jz[bond_index] / 2 *
                    (up[left] * up[right] + down[left] * down[right]))
                if up[left] != up[right]:
                    swapped = source
                    for spin in (0, 1):
                        swapped ^= 1 << (left + 4 * spin)
                        swapped ^= 1 << (right + 4 * spin)
                    matrix[index[swapped], source_position] += jxy[bond_index] / 2
    return states, matrix


def occupied(bits: int, spin: int) -> list[int]:
    return [site for site in range(4) if bits & (1 << (site + 4 * spin))]


def correlator_counts(bits: int) -> dict[str, list[int]]:
    up = [(bits >> site) & 1 for site in range(4)]
    down = [(bits >> (site + 4)) & 1 for site in range(4)]
    occupation = [up[site] + down[site] for site in range(4)]
    charge = [value - 1 for value in occupation]
    spin = [up[site] - down[site] for site in range(4)]
    pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    result = {
        "Gutzwiller": [up[site] * down[site] for site in range(4)],
        "Jastrow": [charge[left] * charge[right] for left, right in pairs],
        "SpinJastrow": [spin[left] * spin[right] for left, right in pairs],
        "DH2": [0] * 6, "DH4": [0] * 10,
    }
    for center in range(4):
        if occupation[center] == 1:
            continue
        xi = occupation[center] // 2
        opposite = 2 if xi == 0 else 0
        near = ((center + 3) % 4, (center + 1) % 4)
        result["DH2"][xi + 2 * sum(occupation[item] == opposite for item in near)] += 1
        around = ((center + 2) % 4, (center + 3) % 4,
                  (center + 1) % 4, (center + 2) % 4)
        result["DH4"][xi + 2 * sum(occupation[item] == opposite for item in around)] += 1
    return result


def independent_psi0(logical: dict[str, Any]) -> np.ndarray:
    descriptor = logical["descriptor"]
    model = 0 if descriptor["model"] == "electronic" else 1
    pair = np.zeros((4, 4), dtype=np.complex128)
    for item in logical["pair_product"]["parameters"]:
        pair[int(item["up_site"]), int(item["down_site"])] = parse_rational_complex(item["value"])
    active = {item["family"]: [float(parse_fraction(value)) for value in item["values"]]
              for item in logical["correlator"]["active_parameters"]}
    values = []
    for bits in configurations(model):
        up = occupied(bits, 0)
        down = occupied(bits, 1)
        total = 0j
        for transform in logical["projection"]["transforms"]:
            mapping = transform["site_mapping"]
            determinant = (
                pair[mapping[up[0]], mapping[down[0]]] *
                pair[mapping[up[1]], mapping[down[1]]] -
                pair[mapping[up[0]], mapping[down[1]]] *
                pair[mapping[up[1]], mapping[down[0]]])
            total -= parse_rational_complex(transform["weight"]) * determinant
        counts = correlator_counts(bits)
        log_factor = sum(parameter * count
                         for family, parameters in active.items()
                         for parameter, count in zip(parameters, counts[family]))
        values.append(total * math.exp(log_factor))
    return np.asarray(values, dtype=np.complex128)


def matrix_multiply(left: list[list[complex]],
                    right: list[list[complex]]) -> list[list[complex]]:
    return [[sum(left[row][inner] * right[inner][column]
                 for inner in range(len(right)))
             for column in range(len(right[0]))]
            for row in range(len(left))]


def conjugate_transpose(matrix: list[list[complex]]) -> list[list[complex]]:
    return [[matrix[row][column].conjugate() for row in range(len(matrix))]
            for column in range(len(matrix[0]))]


def identity(dimension: int) -> list[list[complex]]:
    return [[1.0 + 0j if row == column else 0j
             for column in range(dimension)] for row in range(dimension)]


def hermitian_eigh(matrix: Sequence[Sequence[complex]]) -> tuple[list[float], list[list[complex]]]:
    dimension = len(matrix)
    require(1 <= dimension <= 3 and all(len(row) == dimension for row in matrix),
            "Jacobi dimension")
    values = [[complex(value) for value in row] for row in matrix]
    vectors = identity(dimension)
    for _ in range(128):
        pairs = [(abs(values[row][column]), row, column)
                 for row in range(dimension)
                 for column in range(row + 1, dimension)]
        if not pairs:
            break
        magnitude, p, q = max(pairs, key=lambda item: (item[0], -item[1], -item[2]))
        scale = max(1.0, max(abs(value) for row in values for value in row))
        if magnitude <= 32 * sys.float_info.epsilon * scale:
            break
        phase = values[p][q] / magnitude
        diagonal_phase = identity(dimension)
        diagonal_phase[q][q] = phase.conjugate()
        phased = matrix_multiply(conjugate_transpose(diagonal_phase),
                                 matrix_multiply(values, diagonal_phase))
        app, aqq, apq = phased[p][p].real, phased[q][q].real, phased[p][q].real
        tau = (aqq - app) / (2 * apq)
        tangent = (math.copysign(1.0 / (abs(tau) + math.sqrt(1.0 + tau * tau)), tau)
                   if tau != 0 else 1.0)
        cosine = 1.0 / math.sqrt(1.0 + tangent * tangent)
        sine = tangent * cosine
        rotation = identity(dimension)
        rotation[p][p] = rotation[q][q] = cosine
        rotation[p][q], rotation[q][p] = sine, -sine
        unitary = matrix_multiply(diagonal_phase, rotation)
        values = matrix_multiply(conjugate_transpose(unitary),
                                 matrix_multiply(values, unitary))
        vectors = matrix_multiply(vectors, unitary)
        for row in range(dimension):
            values[row][row] = complex(values[row][row].real, 0.0)
            for column in range(row + 1, dimension):
                average = 0.5 * (values[row][column] + values[column][row].conjugate())
                values[row][column], values[column][row] = average, average.conjugate()
    else:
        raise ValidationError("Jacobi convergence")
    scale = max(1.0, max(abs(value) for row in values for value in row))
    require(max((abs(values[row][column]) for row in range(dimension)
                 for column in range(row + 1, dimension)), default=0.0)
            <= 128 * sys.float_info.epsilon * scale, "Jacobi residual")
    order = sorted(range(dimension), key=lambda index: (values[index][index].real, index))
    return ([values[index][index].real for index in order],
            [[vectors[row][index] for index in order] for row in range(dimension)])


def quadratic(vector: Sequence[complex], matrix: Sequence[Sequence[complex]]) -> complex:
    return sum(vector[row].conjugate() * matrix[row][column] * vector[column]
               for row in range(len(vector)) for column in range(len(vector)))


def matvec(matrix: Sequence[Sequence[complex]],
           vector: Sequence[complex]) -> list[complex]:
    return [sum(matrix[row][column] * vector[column]
                for column in range(len(vector)))
            for row in range(len(matrix))]


def action_quadratic(
        vector: Sequence[complex],
        matrix: Sequence[Sequence[complex]]) -> tuple[complex, list[complex]]:
    """Evaluate c^H M c through the returned matrix action.

    The production API reports diagnostics from a row-wise matrix-vector
    product followed by an inner product.  Keeping that operation order here
    makes the report-consistency check about the returned coefficient, rather
    than about an avoidable reassociation of an ill-conditioned quadratic
    form.  The eigensolver oracle remains independently implemented below.
    """
    action = matvec(matrix, vector)
    value = 0j
    for coefficient, action_value in zip(vector, action):
        value += coefficient.conjugate() * action_value
    return value, action


def dense_from_upper(values: Sequence[complex], dimension: int) -> list[list[complex]]:
    require(len(values) == dimension * (dimension + 1) // 2,
            "upper-packed dimension")
    result = [[0j for _ in range(dimension)] for _ in range(dimension)]
    entry = 0
    for row in range(dimension):
        for column in range(row, dimension):
            value = values[entry]
            result[row][column] = complex(value.real, 0.0) if row == column else value
            result[column][row] = result[row][column].conjugate()
            entry += 1
    return result


def upper_diagonal_indices(dimension: int) -> tuple[int, ...]:
    require(1 <= dimension <= 3, "upper-packed diagonal dimension")
    indices = []
    entry = 0
    for row in range(dimension):
        for column in range(row, dimension):
            if row == column:
                indices.append(entry)
            entry += 1
    return tuple(indices)


def canonical_gevp(overlap: list[list[complex]],
                   hamiltonian: list[list[complex]], cutoff: float) -> dict[str, Any]:
    s_values, s_vectors = hermitian_eigh(overlap)
    largest = max(s_values)
    require(largest > 0.0, "positive overlap spectrum")
    retained = [index for index, value in enumerate(s_values)
                if value > cutoff * largest]
    retained.sort(key=lambda index: (-s_values[index], index))
    require(retained, "retained overlap rank")
    whitening = [[s_vectors[row][index] / math.sqrt(s_values[index])
                  for index in retained] for row in range(len(overlap))]
    reduced = matrix_multiply(conjugate_transpose(whitening),
                              matrix_multiply(hamiltonian, whitening))
    roots, reduced_vectors = hermitian_eigh(reduced)
    generalized = matrix_multiply(whitening, reduced_vectors)
    tolerance = EPS64 * max(1.0, abs(roots[0]))
    degenerate = [index for index, value in enumerate(roots)
                  if value - roots[0] <= tolerance]
    coefficient = None
    for basis_index in range(len(overlap)):
        weights = [sum(generalized[row][root].conjugate() *
                       overlap[row][basis_index]
                       for row in range(len(overlap)))
                   for root in degenerate]
        projected = [sum(generalized[row][root] * weight
                         for root, weight in zip(degenerate, weights))
                     for row in range(len(overlap))]
        norm = quadratic(projected, overlap).real
        if norm > EPS64 * EPS64:
            coefficient = [value / math.sqrt(norm) for value in projected]
            break
    require(coefficient is not None, "canonical root projection")
    pivot = max(range(len(coefficient)), key=lambda index: (abs(coefficient[index]), -index))
    phase = coefficient[pivot].conjugate() / abs(coefficient[pivot])
    coefficient = [value * phase for value in coefficient]
    coefficient[pivot] = complex(abs(coefficient[pivot]), 0.0)
    k_action = matvec(hamiltonian, coefficient)
    s_action = matvec(overlap, coefficient)
    energy = quadratic(coefficient, hamiltonian).real / quadratic(coefficient, overlap).real
    residual = math.sqrt(sum(abs(k_value - energy * s_value) ** 2
                             for k_value, s_value in zip(k_action, s_action)))
    denominator = max(1.0,
                      math.sqrt(sum(abs(value) ** 2 for value in k_action)),
                      abs(energy) * math.sqrt(sum(abs(value) ** 2
                                                  for value in s_action)))
    coefficient_norm = math.sqrt(sum(abs(value) ** 2 for value in coefficient))
    k_frobenius = math.sqrt(sum(abs(value) ** 2 for row in hamiltonian for value in row))
    s_frobenius = math.sqrt(sum(abs(value) ** 2 for row in overlap for value in row))
    backward_denominator = max(1.0, k_frobenius * coefficient_norm,
                               abs(energy) * s_frobenius * coefficient_norm)
    numpy_s_values, numpy_s_vectors = np.linalg.eigh(np.asarray(overlap))
    numpy_retained = [index for index, value in enumerate(numpy_s_values)
                      if value > cutoff * float(numpy_s_values[-1])]
    numpy_w = numpy_s_vectors[:, numpy_retained] / np.sqrt(numpy_s_values[numpy_retained])
    numpy_roots = np.linalg.eigvalsh(numpy_w.conj().T @ np.asarray(hamiltonian) @ numpy_w)
    require(max(abs(left - right) for left, right in zip(s_values, numpy_s_values))
            <= 1e-9 * max(1.0, max(map(abs, s_values))), "Jacobi/NumPy S")
    require(max(abs(left - right) for left, right in zip(roots, numpy_roots))
            <= 1e-9 * max(1.0, max(map(abs, roots))), "Jacobi/NumPy roots")
    return {
        "rank": len(retained), "multiplicity": len(degenerate),
        "coefficient": coefficient, "energy": energy,
        "root_gap": roots[1] - roots[0] if len(roots) > 1 else 0.0,
        "condition": largest / min(s_values[index] for index in retained),
        "raw_action_residual": residual / denominator,
        "backward_error": residual / backward_denominator,
    }


def complex_record(value: complex) -> dict[str, str]:
    return {"re": hex_float(float(value.real)),
            "im": hex_float(float(value.imag))}


def digest_complex(values: Iterable[complex]) -> str:
    return sha256_bytes(canonical_json([complex_record(value) for value in values]))


def close_complex(actual: complex, expected: complex,
                  absolute: float = BINARY64_ABSOLUTE_TOLERANCE,
                  relative: float = BINARY64_RELATIVE_TOLERANCE) -> float:
    error = max(abs(actual.real - expected.real),
                abs(actual.imag - expected.imag))
    scale = max(1.0, abs(expected.real), abs(expected.imag))
    require(error <= absolute or error <= relative * scale,
            f"complex mismatch actual={actual!r} expected={expected!r} error={error}")
    return error


def complex_component_error(actual: complex, expected: complex) -> float:
    return max(abs(actual.real - expected.real),
               abs(actual.imag - expected.imag))


def close_report(actual: float, recomputed: float, label: str) -> None:
    error = abs(actual - recomputed)
    require(error <= REPORT_ABSOLUTE_TOLERANCE or
            error <= REPORT_RELATIVE_TOLERANCE *
                max(abs(recomputed), sys.float_info.min),
            f"{label}: report mismatch actual={actual!r} "
            f"recomputed={recomputed!r} error={error!r}")


def interval_complex_excess(actual: complex,
                            record: dict[str, str]) -> tuple[Decimal, Decimal]:
    radius = Decimal(record["outward_abs_radius"])
    maximum_excess = Decimal(0)
    maximum_scale = Decimal(1)
    for value, key in ((actual.real, "re_mid"), (actual.imag, "im_mid")):
        midpoint = Decimal(record[key])
        excess = abs(Decimal.from_float(float(value)) - midpoint) - radius
        maximum_excess = max(maximum_excess, excess, Decimal(0))
        maximum_scale = max(maximum_scale, abs(midpoint))
    return maximum_excess, maximum_scale


def require_interval_payload(label: str, maximum_excess: Decimal,
                             maximum_scale: Decimal) -> None:
    require(maximum_excess <= Decimal("5e-14") or
            maximum_excess / maximum_scale <=
                Decimal.from_float(BINARY64_RELATIVE_TOLERANCE),
            f"{label}: certified interval excess {maximum_excess} "
            f"at scale {maximum_scale}")


def validate_scaled_state(state: int, error_bound: float, raw: complex,
                          label: str) -> None:
    require(state in (0, 1, 2), f"{label}: scaled state census")
    if state == 0:
        require(raw != 0.0 and
                (error_bound == -math.inf or
                 (math.isfinite(error_bound) and
                  error_bound < math.log(abs(raw)))),
                f"{label}: finite-nonzero state/bound")
    elif state == 1:
        require(raw == 0.0 and error_bound == -math.inf,
                f"{label}: exact-zero state/bound")
    else:
        require(raw == 0.0 and math.isfinite(error_bound),
                f"{label}: numeric-zero state/bound")


def oracle_complex(record: dict[str, str]) -> complex:
    require(set(record) == {"re", "im"} and
            ORACLE_BINARY64_DECIMAL.fullmatch(record["re"]) is not None and
            ORACLE_BINARY64_DECIMAL.fullmatch(record["im"]) is not None,
            "oracle binary64 canonical decimal grammar")
    value = complex(float(record["re"]), float(record["im"]))
    require(math.isfinite(value.real) and math.isfinite(value.imag),
            "oracle binary64 finite")
    return value


def oracle_matrix(record: dict[str, Any], kind: str) -> list[list[complex]]:
    return [[oracle_complex(value) for value in row]
            for row in record["explicit_truth"]["binary64_matrices"][kind]]


def semantic_run_payload(run: ParsedRun) -> dict[str, Any]:
    def convert(value: Any) -> Any:
        if isinstance(value, complex):
            return complex_record(value)
        if isinstance(value, float) and value == -math.inf:
            return "NINF"
        if isinstance(value, dict):
            return {key: convert(item) for key, item in value.items()
                    if key not in ("partitions", "local_factorizations")}
        if isinstance(value, (list, tuple)):
            return [convert(item) for item in value]
        return value
    return {"frozen_shas": list(run.frozen_shas),
            "models": {f"{key[0]}:{key[1]}": convert(value)
                       for key, value in sorted(run.models.items())},
            "cases": convert(run.cases), "negative": list(run.negative)}


def independent_hamiltonian_oracle_sha(states: Sequence[int],
                                        matrix: np.ndarray) -> str:
    payload = []
    for target_index, target in enumerate(states):
        for source_index, source in enumerate(states):
            value = complex(matrix[target_index, source_index])
            if value == 0.0:
                continue
            real = Fraction.from_float(float(value.real))
            imaginary = Fraction.from_float(float(value.imag))
            payload.append({
                "target_configuration_bits": target,
                "source_configuration_bits": source,
                "value": {
                    "re": [real.numerator, real.denominator],
                    "im": [imaginary.numerator, imaginary.denominator],
                },
            })
    return sha256_bytes(canonical_json(payload))


def registered_full_hamiltonian_shas() -> dict[tuple[int, int], str]:
    fixture = load_fixture_module()
    registry = json.loads(REGISTRY.read_text())
    result: dict[tuple[int, int], str] = {}
    for descriptor in fixture.descriptor_census(registry):
        key = (0 if descriptor["model"] == "electronic" else 1,
               0 if descriptor["arithmetic"] == "real" else 1)
        if key in result:
            continue
        logical = fixture.build_logical_fixture(descriptor)
        result[key] = fixture.exact_krylov_oracle(logical)[
            "full_hamiltonian_matrix_sha256"]
        if len(result) == 4:
            break
    require(set(result) == {(0, 0), (0, 1), (1, 0), (1, 1)},
            "registered full-H SHA census")
    return result


def validate_model_blocks(run: ParsedRun) -> dict[str, str]:
    hashes = {}
    registered = registered_full_hamiltonian_shas()
    for key, block in sorted(run.models.items()):
        model, arithmetic = key
        states, expected = independent_hamiltonian(model, arithmetic)
        independent_sha = independent_hamiltonian_oracle_sha(states, expected)
        require(independent_sha == registered[key],
                f"full-H frozen oracle SHA binding {key}")
        require(block["state_count"] == len(states), "model state count")
        require((block["term_count"], block["operator_count"]) ==
                ((44, 144) if model == 0 else (16, 64)),
                f"model term/operator census {key}")
        require(len(block["rows"]) == len(states) ** 2, "model dense census")
        dense = np.zeros_like(expected)
        entry = 0
        for target_index, target in enumerate(states):
            for source_index, source in enumerate(states):
                output_target, output_source, value = block["rows"][entry]
                require((output_target, output_source) == (target, source),
                        "model configuration order")
                close_complex(value, complex(expected[target_index, source_index]),
                              absolute=1e-14, relative=1e-15)
                dense[target_index, source_index] = value
                entry += 1
        require(np.max(np.abs(dense - dense.conj().T)) == 0.0,
                "production full H Hermitian")
        hashes[f"{model}:{arithmetic}"] = digest_complex(dense.flat)
    return hashes


def flatten_upper(matrix: Sequence[Sequence[complex]]) -> list[complex]:
    return [matrix[row][column] for row in range(len(matrix))
            for column in range(row, len(matrix))]


def projective_distance(left: Sequence[complex], right: Sequence[complex],
                        overlap: Sequence[Sequence[complex]]) -> float:
    cross = sum(left[row].conjugate() * overlap[row][column] * right[column]
                for row in range(len(left)) for column in range(len(left)))
    require(abs(cross) > 0.0, "nonzero projective overlap")
    phase = cross.conjugate() / abs(cross)
    difference = [left[index] - phase * right[index]
                  for index in range(len(left))]
    return math.sqrt(max(0.0, quadratic(difference, overlap).real))


def validate_case(case: dict[str, Any], logical: dict[str, Any],
                  oracle: dict[str, Any], world_size: int) -> dict[str, Any]:
    descriptor = logical["descriptor"]
    fixture_id = descriptor["fixture_id"]
    require(case["fixture_id"] == fixture_id, "case fixture identity")
    expected_axes = (
        0 if descriptor["model"] == "electronic" else 1,
        int(descriptor["order"]),
        0 if descriptor["arithmetic"] == "real" else 1,
        int(descriptor["nqp_full"]),
        PROFILE_IDS[str(descriptor["correlator_profile"])],
        SUPPORT_IDS[str(descriptor["support_case"])]
    )
    require(tuple(case[key] for key in
                  ("model", "order", "arithmetic", "nqp", "profile", "support"))
            == expected_axes, f"{fixture_id}: axis binding")
    fixture = load_fixture_module()
    require(case["logical_sha"] == sha256_bytes(fixture.canonical_bytes(logical)),
            f"{fixture_id}: logical SHA")
    rendered = fixture.rendered_manifest(fixture.render_classic_files(logical))
    require(case["rendered_sha"] == rendered["aggregate_sha256"],
            f"{fixture_id}: rendered SHA")
    states = configurations(case["model"])
    require(case["state_count"] == len(states) and
            [row["bits"] for row in case["rows"]] == states,
            f"{fixture_id}: configuration order")
    require(case["pivot"] == (0.0, float.fromhex("0x1p-44")),
            f"{fixture_id}: pivot policy")
    require(case["partitions"] == [
        (rank, case["nqp"] * rank // world_size,
         case["nqp"] * (rank + 1) // world_size)
        for rank in range(world_size)], f"{fixture_id}: QP partition")
    proof = logical["support_proof"]
    expected_proof = ((int(proof["x_configuration_bits"]),
                       int(proof["y_configuration_bits"]))
                      if proof else (0, 0))
    require(case["proof"] == expected_proof, f"{fixture_id}: proof binding")

    truth_vectors = [[oracle_complex(item) for item in vector]
                     for vector in oracle["explicit_truth"]["binary64_vectors"]]
    exact_vectors = oracle["explicit_truth"]["vectors"]
    independent_zero = independent_psi0(logical)
    require(len(truth_vectors) == case["order"] + 2 and
            len(exact_vectors) == len(truth_vectors), "oracle vector count")
    maximum_vector_error = 0.0
    maximum_vector_interval_excess = Decimal(0)
    amplitude_binary_error = 0.0
    amplitude_binary_scale = 1.0
    amplitude_interval_excess = Decimal(0)
    amplitude_interval_scale = Decimal(1)
    vector_binary_errors = [0.0 for _ in truth_vectors]
    vector_interval_excesses = [Decimal(0) for _ in truth_vectors]
    vector_interval_scales = [Decimal(1) for _ in truth_vectors]
    numeric_zero_count = 0
    for state_index, row in enumerate(case["rows"]):
        maximum_vector_error = max(
            maximum_vector_error,
            complex_component_error(
                row["amplitude"], complex(independent_zero[state_index])))
        require(all(row[key] >= 0 for key in (
                    "well", "near", "exact_zero", "numeric_zero",
                    "local_factorizations", "global_factorizations")),
                f"{fixture_id}: nonnegative metadata census")
        require(row["well"] + row["near"] + row["exact_zero"] == case["nqp"],
                f"{fixture_id}: component state census")
        require(row["numeric_zero"] <= row["near"],
                f"{fixture_id}: numeric zero census")
        require(row["global_factorizations"] == case["nqp"] and
                row["local_factorizations"] ==
                    case["partitions"][0][2] - case["partitions"][0][1],
                f"{fixture_id}: factorization census")
        expected_amplitude = truth_vectors[0][state_index]
        amplitude_error = complex_component_error(
            row["amplitude"], expected_amplitude)
        amplitude_binary_error = max(amplitude_binary_error,
                                     amplitude_error)
        amplitude_binary_scale = max(
            amplitude_binary_scale, abs(expected_amplitude.real),
            abs(expected_amplitude.imag))
        maximum_vector_error = max(maximum_vector_error, amplitude_error)
        interval_excess, interval_scale = interval_complex_excess(
            row["amplitude"], exact_vectors[0][state_index])
        amplitude_interval_excess = max(amplitude_interval_excess,
                                        interval_excess)
        amplitude_interval_scale = max(amplitude_interval_scale,
                                       interval_scale)
        maximum_vector_interval_excess = max(
            maximum_vector_interval_excess, interval_excess)
        validate_scaled_state(row["scaled_state"],
                              row["amplitude_error_bound"],
                              row["amplitude"],
                              f"{fixture_id}: amplitude row {state_index}")
        if row["scaled_state"] == 2:
            numeric_zero_count += 1
            require(expected_amplitude == 0.0 and
                    math.isfinite(row["amplitude_error_bound"]),
                    f"{fixture_id}: amplitude numeric-zero proof")
        elif expected_amplitude == 0.0:
            require(row["scaled_state"] == 1,
                    f"{fixture_id}: amplitude exact-zero state")
        else:
            require(row["scaled_state"] == 0,
                    f"{fixture_id}: amplitude finite state")
        for vector_index, actual in enumerate(row["vectors"]):
            expected = truth_vectors[vector_index][state_index]
            binary_error = complex_component_error(actual, expected)
            maximum_vector_error = max(maximum_vector_error, binary_error)
            vector_binary_errors[vector_index] = max(
                vector_binary_errors[vector_index], binary_error)
            exact_record = exact_vectors[vector_index][state_index]
            require(exact_record["configuration_bits"] == row["bits"],
                    f"{fixture_id}: exact vector configuration binding")
            interval_excess, interval_scale = interval_complex_excess(
                actual, exact_record)
            vector_interval_excesses[vector_index] = max(
                vector_interval_excesses[vector_index], interval_excess)
            vector_interval_scales[vector_index] = max(
                vector_interval_scales[vector_index], interval_scale)
            maximum_vector_interval_excess = max(
                maximum_vector_interval_excess, interval_excess)
            metadata = row["vector_states"][vector_index]
            validate_scaled_state(
                metadata["state"], metadata["error_bound"], actual,
                f"{fixture_id}: v{vector_index} row {state_index}")
            if metadata["state"] == 2:
                numeric_zero_count += 1
                require(expected == 0.0 and math.isfinite(metadata["error_bound"]),
                        f"{fixture_id}: vector numeric-zero proof")
            elif expected == 0.0:
                require(metadata["state"] == 1,
                        f"{fixture_id}: vector exact-zero state")
            else:
                require(metadata["state"] == 0,
                        f"{fixture_id}: vector finite state")
    require(amplitude_binary_error <= BINARY64_ABSOLUTE_TOLERANCE or
            amplitude_binary_error <=
                BINARY64_RELATIVE_TOLERANCE * amplitude_binary_scale,
            f"{fixture_id}: amplitude binary64 payload")
    require_interval_payload(
        f"{fixture_id}: amplitude exact interval payload",
        amplitude_interval_excess, amplitude_interval_scale)
    maximum_vector_relative_error = (
        amplitude_binary_error / amplitude_binary_scale)
    maximum_vector_relative_payload = "amplitude"
    maximum_vector_interval_relative_excess = (
        amplitude_interval_excess / amplitude_interval_scale)
    maximum_vector_interval_relative_payload = "amplitude"
    for vector_index in range(len(truth_vectors)):
        vector_scale = max(
            1.0,
            max(abs(value.real) for value in truth_vectors[vector_index]),
            max(abs(value.imag) for value in truth_vectors[vector_index]))
        require(vector_binary_errors[vector_index] <=
                    BINARY64_ABSOLUTE_TOLERANCE or
                vector_binary_errors[vector_index] <=
                    BINARY64_RELATIVE_TOLERANCE * vector_scale,
                f"{fixture_id}: v{vector_index} binary64 payload")
        require_interval_payload(
            f"{fixture_id}: v{vector_index} exact interval payload",
            vector_interval_excesses[vector_index],
            vector_interval_scales[vector_index])
        relative_error = vector_binary_errors[vector_index] / vector_scale
        if relative_error > maximum_vector_relative_error:
            maximum_vector_relative_error = relative_error
            maximum_vector_relative_payload = f"v{vector_index}"
        relative_excess = (vector_interval_excesses[vector_index] /
                           vector_interval_scales[vector_index])
        if relative_excess > maximum_vector_interval_relative_excess:
            maximum_vector_interval_relative_excess = relative_excess
            maximum_vector_interval_relative_payload = f"v{vector_index}"
    zero_support_proof = None
    if case["support"] == 1:
        proof_y = case["proof"][1]
        proof_row = next(row for row in case["rows"]
                         if row["bits"] == proof_y)
        require(proof_row["well"] == 0 and proof_row["near"] == 0 and
                proof_row["exact_zero"] == case["nqp"] and
                proof_row["numeric_zero"] == 0 and
                proof_row["scaled_state"] == 1 and
                proof_row["amplitude_error_bound"] == -math.inf,
                f"{fixture_id}: zero-support proof-y structural metadata")

    dimension = case["order"] + 1
    diagonal_indices = upper_diagonal_indices(dimension)
    for name in ("S", "B"):
        require(all(case["packs"][name][index].imag == 0.0
                    for index in diagonal_indices),
                f"{fixture_id}: {name} diagonal imaginary exact +0")
    require(all((0.5 * (case["packs"]["KF"][index] +
                        case["packs"]["KR"][index].conjugate())).imag == 0.0
                for index in diagonal_indices),
            f"{fixture_id}: Hermitized K diagonal imaginary exact +0")
    truth_matrices = {kind: oracle_matrix(oracle, kind)
                      for kind in ("S", "K", "B")}
    exact_matrices = oracle["explicit_truth"]["matrices"]
    expected_packs = {kind: flatten_upper(matrix)
                      for kind, matrix in truth_matrices.items()}
    maximum_matrix_error = 0.0
    for name, expected_kind in (("S", "S"), ("KF", "K"),
                                ("B", "B")):
        for actual, expected in zip(case["packs"][name],
                                    expected_packs[expected_kind]):
            maximum_matrix_error = max(
                maximum_matrix_error,
                complex_component_error(actual, expected))
    for actual, expected in zip(case["packs"]["KR"], expected_packs["K"]):
        maximum_matrix_error = max(
            maximum_matrix_error,
            complex_component_error(actual.conjugate(), expected))
    overlap = dense_from_upper(case["packs"]["S"], dimension)
    forward = dense_from_upper(case["packs"]["KF"], dimension)
    reverse = dense_from_upper(case["packs"]["KR"], dimension)
    hamiltonian = [[0.5 * forward[row][column] +
                    0.5 * reverse[column][row]
                    for column in range(dimension)] for row in range(dimension)]
    squared_matrix = dense_from_upper(case["packs"]["B"], dimension)
    matrix_actuals = {
        "S": case["packs"]["S"],
        "KF": case["packs"]["KF"],
        "KR": [value.conjugate() for value in case["packs"]["KR"]],
        "K": flatten_upper(hamiltonian),
        "B": case["packs"]["B"],
    }
    matrix_kinds = {"S": "S", "KF": "K", "KR": "K", "K": "K",
                    "B": "B"}
    maximum_matrix_interval_excess = Decimal(0)
    maximum_matrix_relative_error = 0.0
    maximum_matrix_relative_payload = "S"
    maximum_matrix_interval_relative_excess = Decimal(0)
    maximum_matrix_interval_relative_payload = "S"
    for name, actual_values in matrix_actuals.items():
        kind = matrix_kinds[name]
        exact_values = flatten_upper(exact_matrices[kind])
        binary_values = expected_packs[kind]
        require(len(actual_values) == len(exact_values) == len(binary_values),
                f"{fixture_id}: {name} matrix payload census")
        payload_binary_error = 0.0
        payload_interval_excess = Decimal(0)
        payload_interval_scale = Decimal(1)
        binary_scale = 1.0
        for actual_value, exact_record, binary_value in zip(
                actual_values, exact_values, binary_values):
            payload_binary_error = max(
                payload_binary_error,
                max(abs(actual_value.real - binary_value.real),
                    abs(actual_value.imag - binary_value.imag)))
            binary_scale = max(binary_scale, abs(binary_value.real),
                               abs(binary_value.imag))
            interval_excess, interval_scale = interval_complex_excess(
                actual_value, exact_record)
            payload_interval_excess = max(payload_interval_excess,
                                          interval_excess)
            payload_interval_scale = max(payload_interval_scale,
                                         interval_scale)
        require(payload_binary_error <= BINARY64_ABSOLUTE_TOLERANCE or
                payload_binary_error <=
                    BINARY64_RELATIVE_TOLERANCE * binary_scale,
                f"{fixture_id}: {name} binary64 payload")
        require_interval_payload(f"{fixture_id}: {name} exact interval payload",
                                 payload_interval_excess,
                                 payload_interval_scale)
        maximum_matrix_error = max(maximum_matrix_error,
                                   payload_binary_error)
        maximum_matrix_interval_excess = max(
            maximum_matrix_interval_excess, payload_interval_excess)
        relative_error = payload_binary_error / binary_scale
        if relative_error > maximum_matrix_relative_error:
            maximum_matrix_relative_error = relative_error
            maximum_matrix_relative_payload = name
        relative_excess = payload_interval_excess / payload_interval_scale
        if relative_excess > maximum_matrix_interval_relative_excess:
            maximum_matrix_interval_relative_excess = relative_excess
            maximum_matrix_interval_relative_payload = name
    maximum_backward = 0.0
    maximum_raw = 0.0
    maximum_independent_backward = 0.0
    maximum_independent_raw = 0.0
    maximum_observable_reference_abs_error = 0.0
    cutoff_records: dict[str, dict[str, Any]] = {}
    if case["support"] == 2:
        require(all(value == 0.0 for pack in case["packs"].values()
                    for value in pack), f"{fixture_id}: singular exact matrices")
        require(all(result == {"status": "NQP_SINGULAR"}
                    for result in case["solves"].values()),
                f"{fixture_id}: singular skip")
        cutoff_records = {
            cutoff_id: {"status": "NQP_SINGULAR", "verdict": "PASS"}
            for cutoff_id in CUTOFFS
        }
    else:
        reference_coefficient = None
        for cutoff_id, cutoff in CUTOFFS.items():
            actual = case["solves"][cutoff_id]
            expected = canonical_gevp(overlap, hamiltonian, cutoff)
            require(actual["dimension"] == dimension and
                    actual["rank"] == dimension and actual["discarded"] == 0 and
                    actual["multiplicity"] == 1 and
                    expected["rank"] == dimension and expected["multiplicity"] == 1,
                    f"{fixture_id}: GEVP rank/root")
            require(projective_distance(actual["coefficient"],
                                        expected["coefficient"], overlap) <= 1e-9,
                    f"{fixture_id}: coefficient projective distance")
            actual_c = actual["coefficient"]
            actual_pivot = max(range(dimension),
                               key=lambda index: (abs(actual_c[index]), -index))
            expected_pivot = max(
                range(dimension),
                key=lambda index: (abs(expected["coefficient"][index]), -index))
            require(actual["pivot"] == actual_pivot == expected_pivot and
                    actual_c[actual_pivot].real > 0.0 and
                    actual_c[actual_pivot].imag == 0.0,
                    f"{fixture_id}: canonical coefficient phase/pivot")
            coefficient_scale = max(1.0, *(abs(value) for value in
                                           expected["coefficient"]))
            require(max(abs(left - right) for left, right in zip(
                        actual_c, expected["coefficient"])) <=
                    1e-9 * coefficient_scale,
                    f"{fixture_id}: phase-fixed coefficient")
            actual_normalization_value, actual_s_action = action_quadratic(
                actual_c, overlap)
            require(abs(actual_normalization_value.imag) <= 1e-12 and
                    abs(actual_normalization_value.real - 1.0) <= 1e-12,
                    f"{fixture_id}: S normalization")
            close_report(actual["normalization"],
                         actual_normalization_value.real,
                         f"{fixture_id}: normalization")
            actual_energy_value, actual_k_action = action_quadratic(
                actual_c, hamiltonian)
            actual_m2_value, _ = action_quadratic(actual_c, squared_matrix)
            actual_energy = (actual_energy_value.real /
                             actual_normalization_value.real)
            actual_matrix_m2 = (actual_m2_value.real /
                                actual_normalization_value.real)
            recomputed_matrix_variance = (
                actual_matrix_m2 - actual_energy * actual_energy)
            actual_residual_norm = math.sqrt(sum(
                abs(k_value - actual_energy * s_value) ** 2
                for k_value, s_value in zip(actual_k_action, actual_s_action)))
            actual_c_norm = math.sqrt(sum(abs(value) ** 2 for value in actual_c))
            k_frobenius = math.sqrt(sum(abs(value) ** 2
                                        for row in hamiltonian for value in row))
            s_frobenius = math.sqrt(sum(abs(value) ** 2
                                        for row in overlap for value in row))
            recomputed_backward = actual_residual_norm / max(
                1.0, k_frobenius * actual_c_norm,
                abs(actual_energy) * s_frobenius * actual_c_norm)
            recomputed_raw = actual_residual_norm / max(
                1.0,
                math.sqrt(sum(abs(value) ** 2 for value in actual_k_action)),
                abs(actual_energy) * math.sqrt(sum(abs(value) ** 2
                                                   for value in actual_s_action)))
            close_report(actual["backward_error"], recomputed_backward,
                         f"{fixture_id}: backward error")
            close_report(actual["raw_action_residual"], recomputed_raw,
                         f"{fixture_id}: raw action residual")
            close_report(actual["root_gap"], expected["root_gap"],
                         f"{fixture_id}: root gap")
            close_report(actual["condition"], expected["condition"],
                         f"{fixture_id}: condition estimate")
            require(actual["backward_error"] <= EPS64 and
                    recomputed_backward <= EPS64,
                    f"{fixture_id}: backward error gate")
            close_report(actual["packed_energy"], actual_energy,
                         f"{fixture_id}: packed energy")
            close_report(actual["packed_m2"], actual_matrix_m2,
                         f"{fixture_id}: packed M2")
            close_report(actual["raw_matrix_variance"],
                         recomputed_matrix_variance,
                         f"{fixture_id}: raw matrix variance")
            matrix_variance_scale = max(
                1.0, abs(actual_matrix_m2), actual_energy * actual_energy)
            require(actual["matrix_observables_valid"] == int(
                        recomputed_matrix_variance >=
                        -EPS64 * matrix_variance_scale) and
                    actual["variance_clamped"] == 0,
                    f"{fixture_id}: raw matrix diagnostic semantics")
            coefficient = np.asarray(actual["coefficient"])
            vector_array = np.asarray(truth_vectors)
            psi = coefficient @ vector_array[:dimension]
            hpsi = coefficient @ vector_array[1:dimension + 1]
            output_psi = np.asarray([item[1] for item in actual["states"]])
            output_hpsi = np.asarray([item[2] for item in actual["states"]])
            require([item[0] for item in actual["states"]] == states,
                    f"{fixture_id}: optimized state order")
            require(np.linalg.norm(output_psi - psi) <= 1e-9 * max(1.0, np.linalg.norm(psi)) and
                    np.linalg.norm(output_hpsi - hpsi) <= 1e-9 * max(1.0, np.linalg.norm(hpsi)),
                    f"{fixture_id}: optimized state")
            norm = float(np.vdot(output_psi, output_psi).real)
            energy = float(np.vdot(output_psi, output_hpsi).real / norm)
            m2 = float(np.vdot(output_hpsi, output_hpsi).real / norm)
            variance = float(np.vdot(output_hpsi - energy * output_psi,
                                     output_hpsi - energy * output_psi).real / norm)
            reference_c = np.asarray(expected["coefficient"])
            reference_psi = reference_c @ vector_array[:dimension]
            reference_hpsi = reference_c @ vector_array[1:dimension + 1]
            phase_overlap = np.vdot(reference_psi, output_psi)
            require(abs(phase_overlap) > 0.0,
                    f"{fixture_id}: state phase overlap")
            state_phase = phase_overlap.conjugate() / abs(phase_overlap)
            require(np.linalg.norm(state_phase * output_psi - reference_psi)
                    <= 1e-9 * max(1.0, np.linalg.norm(reference_psi)) and
                    np.linalg.norm(state_phase * output_hpsi - reference_hpsi)
                    <= 1e-9 * max(1.0, np.linalg.norm(reference_hpsi)),
                    f"{fixture_id}: independent phase-aligned state")
            reference_norm = float(np.vdot(reference_psi,
                                            reference_psi).real)
            reference_energy = float(
                np.vdot(reference_psi, reference_hpsi).real /
                reference_norm)
            reference_m2 = float(
                np.vdot(reference_hpsi, reference_hpsi).real /
                reference_norm)
            reference_variance = float(np.vdot(
                reference_hpsi - reference_energy * reference_psi,
                reference_hpsi - reference_energy * reference_psi).real /
                reference_norm)
            for label, reported, recomputed in (
                    ("energy", actual["energy"], energy),
                    ("M2", actual["m2"], m2),
                    ("variance", actual["variance"], variance)):
                require(abs(reported - recomputed) <= 1e-10 or
                        abs(reported - recomputed) <= 1e-11 * max(1.0, abs(recomputed)),
                        f"{fixture_id}: state {label}")
            for label, production_value, reference_value in (
                    ("energy", energy, reference_energy),
                    ("M2", m2, reference_m2),
                    ("variance", variance, reference_variance)):
                difference = abs(production_value - reference_value)
                maximum_observable_reference_abs_error = max(
                    maximum_observable_reference_abs_error, difference)
                require(difference <= 1e-10 or
                        difference <=
                            1e-11 * max(1.0, abs(reference_value)),
                        f"{fixture_id}: independent observable {label}")
            for label, packed, state_value in (
                    ("energy", actual_energy, energy),
                    ("M2", actual_matrix_m2, m2)):
                difference = abs(packed - state_value)
                require(difference <= 1e-10 or
                        difference <= 1e-11 * max(1.0, abs(state_value)),
                        f"{fixture_id}: packed/state {label}")
            require(variance >= 0.0 and math.isfinite(variance),
                    f"{fixture_id}: nonnegative variance")
            cutoff_records[cutoff_id] = {
                "status": "OK",
                "dimension": actual["dimension"],
                "retained_rank": actual["rank"],
                "discarded_rank": actual["discarded"],
                "root_multiplicity": actual["multiplicity"],
                "phase_pivot": actual["pivot"],
                "normalization": hex_float(actual["normalization"]),
                "state_energy": hex_float(actual["energy"]),
                "state_second_moment": hex_float(actual["m2"]),
                "state_variance": hex_float(actual["variance"]),
                "packed_energy": hex_float(actual["packed_energy"]),
                "packed_second_moment": hex_float(actual["packed_m2"]),
                "packed_variance": hex_float(
                    actual["raw_matrix_variance"]),
                "normwise_backward_error": hex_float(
                    actual["backward_error"]),
                "independent_normwise_backward_error":
                    hex_float(recomputed_backward),
                "raw_action_relative_residual": hex_float(
                    actual["raw_action_residual"]),
                "independent_raw_action_relative_residual":
                    hex_float(recomputed_raw),
                "root_gap": hex_float(actual["root_gap"]),
                "condition_estimate": hex_float(actual["condition"]),
                "matrix_observables_valid":
                    actual["matrix_observables_valid"],
                "variance_clamped": actual["variance_clamped"],
                "coefficient_sha256": digest_complex(actual_c),
                "independent_coefficient_sha256": digest_complex(
                    expected["coefficient"]),
                "corrected_state_sha256": digest_complex(
                    value for _, psi_value, hpsi_value in actual["states"]
                    for value in (psi_value, hpsi_value)),
                "verdict": "PASS",
            }
            maximum_backward = max(maximum_backward, actual["backward_error"])
            maximum_raw = max(maximum_raw, actual["raw_action_residual"])
            maximum_independent_backward = max(
                maximum_independent_backward, recomputed_backward)
            maximum_independent_raw = max(maximum_independent_raw,
                                          recomputed_raw)
            if reference_coefficient is None:
                reference_coefficient = coefficient
            else:
                require(projective_distance(reference_coefficient, coefficient,
                                            overlap) <= 1e-9,
                        f"{fixture_id}: cutoff invariance")
    if case["support"] == 1:
        x, y = case["proof"]
        state_index = {bits: index for index, bits in enumerate(states)}
        _, full_h = independent_hamiltonian(case["model"], case["arithmetic"])
        require(x != y and independent_zero[state_index[x]] != 0.0 and
                independent_zero[state_index[y]] == 0.0 and
                full_h[state_index[y], state_index[x]] != 0.0 and
                truth_vectors[1][state_index[y]] != 0.0,
                f"{fixture_id}: zero-support bridge")
        zero_support_proof = {
            "x_configuration_bits": x,
            "y_configuration_bits": y,
            "psi0_x_nonzero": True,
            "psi0_y_exact_zero": True,
            "full_h_yx_nonzero": True,
            "h_psi0_y_nonzero": True,
            "all_qp_components_y_structural_zero": True,
            "verdict": "PASS",
        }
    amplitude_scaled_metadata = [{
        "configuration_bits": row["bits"],
        "state": row["scaled_state"],
        "log_abs_error_bound": (
            "NINF" if row["amplitude_error_bound"] == -math.inf
            else hex_float(row["amplitude_error_bound"])),
        "well_pivoted_component_count": row["well"],
        "near_pivot_component_count": row["near"],
        "exact_zero_component_count": row["exact_zero"],
        "numeric_zero_component_count": row["numeric_zero"],
        "local_factorization_count": row["local_factorizations"],
        "global_factorization_count": row["global_factorizations"],
    } for row in case["rows"]]
    return {
        "fixture_id": fixture_id,
        "descriptor": descriptor,
        "logical_sha256": case["logical_sha"],
        "rendered_sha256": case["rendered_sha"],
        "generation_hash": str(case["generation_hash"]),
        "plan_hash": str(case["plan_hash"]),
        "vector_sha256": digest_complex(
            value for row in case["rows"] for value in row["vectors"]),
        "matrix_sha256": sha256_bytes(canonical_json({
            key: [complex_record(value) for value in values]
            for key, values in case["packs"].items()})),
        "registered_oracle_sha256": sha256_bytes(canonical_json(oracle)),
        "amplitude_scaled_metadata": amplitude_scaled_metadata,
        "amplitude_scaled_metadata_sha256":
            sha256_bytes(canonical_json(amplitude_scaled_metadata)),
        "cutoffs": cutoff_records,
        "zero_support_proof": zero_support_proof,
        "numeric_zero_metadata_count": numeric_zero_count,
        "maximum_vector_binary64_abs_error": hex_float(maximum_vector_error),
        "maximum_amplitude_binary64_abs_error":
            hex_float(amplitude_binary_error),
        "maximum_amplitude_exact_interval_excess":
            str(amplitude_interval_excess),
        "maximum_vector_exact_interval_excess":
            str(maximum_vector_interval_excess),
        "maximum_vector_binary64_scale_relative_error":
            hex_float(maximum_vector_relative_error),
        "maximum_vector_binary64_scale_relative_payload":
            maximum_vector_relative_payload,
        "maximum_vector_exact_interval_scale_relative_excess":
            str(maximum_vector_interval_relative_excess),
        "maximum_vector_exact_interval_scale_relative_payload":
            maximum_vector_interval_relative_payload,
        "maximum_matrix_binary64_abs_error": hex_float(maximum_matrix_error),
        "maximum_matrix_exact_interval_excess":
            str(maximum_matrix_interval_excess),
        "maximum_matrix_binary64_scale_relative_error":
            hex_float(maximum_matrix_relative_error),
        "maximum_matrix_binary64_scale_relative_payload":
            maximum_matrix_relative_payload,
        "maximum_matrix_exact_interval_scale_relative_excess":
            str(maximum_matrix_interval_relative_excess),
        "maximum_matrix_exact_interval_scale_relative_payload":
            maximum_matrix_interval_relative_payload,
        "maximum_normwise_backward_error": hex_float(maximum_backward),
        "maximum_independent_normwise_backward_error":
            hex_float(maximum_independent_backward),
        "maximum_raw_action_relative_residual": hex_float(maximum_raw),
        "maximum_independent_raw_action_relative_residual":
            hex_float(maximum_independent_raw),
        "maximum_observable_reference_abs_error":
            hex_float(maximum_observable_reference_abs_error),
        "gate_verdicts": {
            "identity_and_grammar": "PASS",
            "production_full_h_binding": "PASS",
            "scaled_metadata": "PASS",
            "binary64_vectors_and_matrices": "PASS",
            "certified_exact_vectors_and_matrices": "PASS",
            "gevp_and_corrected_state": "PASS",
            "zero_support_or_singular_semantics": "PASS",
        },
        "verdict": "PASS",
    }


def validate_preflight_payload(payload: bytes, expected_world_size: int) -> None:
    validate_numerical_policy(accepted_required=False)
    run = parse_output(payload, preflight=True)
    require(run.world_size == expected_world_size,
            "preflight world-size binding")
    validate_model_blocks(run)
    logical_rows = preflight_logical_rows()
    fixture = load_fixture_module()
    require(len(run.cases) == len(logical_rows) ==
            len(PREFLIGHT_FIXTURE_IDS), "preflight case census")
    for case, logical in zip(run.cases, logical_rows):
        oracle = fixture.exact_krylov_oracle(logical, include_explicit=True)
        validate_case(case, logical, oracle, run.world_size)
    print(f"PASS P6-C1 preflight world_size={run.world_size} "
          f"cases={len(run.cases)} raw_sha256={run.raw_sha256}")


def semantic_json_value(value: Any) -> Any:
    if isinstance(value, complex):
        return complex_record(value)
    if isinstance(value, float):
        if value == -math.inf:
            return "NINF"
        return hex_float(value)
    if isinstance(value, tuple):
        return [semantic_json_value(item) for item in value]
    if isinstance(value, list):
        return [semantic_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): semantic_json_value(item)
                for key, item in value.items()}
    return value


def normalized_run_sha256(run: ParsedRun) -> str:
    models = []
    for (model, arithmetic), payload in sorted(run.models.items()):
        models.append({"model": model, "arithmetic": arithmetic,
                       "payload": semantic_json_value(payload)})
    cases = []
    for payload in run.cases:
        normalized_case = {
            key: value for key, value in payload.items()
            if key != "partitions"
        }
        normalized_case["rows"] = [{
            key: value for key, value in row.items()
            if key != "local_factorizations"
        } for row in payload["rows"]]
        cases.append(semantic_json_value(normalized_case))
    return sha256_bytes(canonical_json({
        "frozen_shas": list(run.frozen_shas),
        "models": models,
        "cases": cases,
        "negative": list(run.negative),
    }))


def validation_record_for_run(payload: bytes, expected_world_size: int,
                              expected_policy_sha: str) -> tuple[ParsedRun, list[dict[str, Any]]]:
    run = parse_output(payload, preflight=True,
                       expected_policy_sha=expected_policy_sha)
    require(run.world_size == expected_world_size,
            "pilot world-size binding")
    validate_model_blocks(run)
    fixture = load_fixture_module()
    logical_rows = preflight_logical_rows()
    require(len(run.cases) == len(logical_rows) == len(PREFLIGHT_FIXTURE_IDS),
            "pilot case census")
    records = []
    for case, logical in zip(run.cases, logical_rows):
        oracle = fixture.exact_krylov_oracle(logical, include_explicit=True)
        records.append(validate_case(case, logical, oracle, run.world_size))
    return run, records


def write_focused_pilot_evidence(
        input_path: Path, serial_path: Path, mpi2_path: Path,
        mpi4_path: Path, execution_driver_sha256: str,
        execution_binary_sha256: str,
        compiler_identity: str) -> None:
    pilot_policy_sha = sha256_file(PILOT_POLICY_SNAPSHOT)
    require(pilot_policy_sha ==
            "190a2053aaff9008ec2a0e8560b5377519558649385afe9d5d67436abc8816a7",
            "pilot policy snapshot identity")
    require(SHA256.fullmatch(execution_driver_sha256) is not None and
            SHA256.fullmatch(execution_binary_sha256) is not None,
            "pilot execution SHA grammar")
    current_driver = DRIVER_SOURCE.read_bytes()
    reconstructed_driver = re.sub(
        rb'(?<=static const char NUMERICAL_POLICY_SHA\[65\] =\n    ")[0-9a-f]{64}(?=";)',
        pilot_policy_sha.encode(), current_driver)
    require(reconstructed_driver != current_driver and
            sha256_bytes(reconstructed_driver) == execution_driver_sha256,
            "pilot driver source reconstruction")
    input_bytes = input_path.read_bytes()
    expected_input, _ = input_payload(
        preflight=True, policy_sha_override=pilot_policy_sha)
    require(input_bytes == expected_input,
            "pilot input canonical byte regeneration")
    run_payloads = {
        "serial": (1, serial_path.read_bytes()),
        "mpi2": (2, mpi2_path.read_bytes()),
        "mpi4": (4, mpi4_path.read_bytes()),
    }
    parsed_runs: dict[str, ParsedRun] = {}
    records_by_run: dict[str, list[dict[str, Any]]] = {}
    normalized_shas = set()
    for name, (world_size, payload) in run_payloads.items():
        run, records = validation_record_for_run(
            payload, world_size, pilot_policy_sha)
        parsed_runs[name] = run
        records_by_run[name] = records
        normalized_shas.add(normalized_run_sha256(run))
    require(len(normalized_shas) == 1,
            "pilot serial/MPI2/MPI4 normalized payload identity")
    generated = {
        PILOT_INPUT_GZIP: deterministic_gzip(input_bytes),
        PILOT_SERIAL_GZIP: deterministic_gzip(run_payloads["serial"][1]),
        PILOT_MPI2_GZIP: deterministic_gzip(run_payloads["mpi2"][1]),
        PILOT_MPI4_GZIP: deterministic_gzip(run_payloads["mpi4"][1]),
    }
    for path, payload in generated.items():
        path.write_bytes(payload)
        require(path.read_bytes() == payload and
                deterministic_gzip(gzip.decompress(payload)) == payload,
                f"deterministic pilot gzip {path.name}")
    evidence = {
        "schema_version": 1,
        "evidence_id": "power_lanczos_p6c1_focused_pilot_evidence_v1",
        "status": "observed_development_pilot",
        "calibration_or_confirmation_result": False,
        "production_authorized": False,
        "selection_rule": "fixed integration coverage selected before the focused run; these nine fixtures are excluded from the protected 2007-fixture holdout claim",
        "fixture_ids": list(PREFLIGHT_FIXTURE_IDS),
        "fixture_count": len(PREFLIGHT_FIXTURE_IDS),
        "holdout_definition": "all 2016 registered fixtures minus fixture_ids",
        "holdout_count": CASE_COUNT - len(PREFLIGHT_FIXTURE_IDS),
        "policy_influence": {
            "interim_post_observation_threshold_change": [
                "vector/matrix scale-relative tolerance was temporarily changed from 5e-15 to 64*DBL_EPSILON",
                "diagnostic report absolute tolerance was temporarily changed from 1e-18 to 4*DBL_EPSILON",
            ],
            "final_holdout_threshold_changes": [],
            "resolution": "Both interim relaxations were reverted after identifying checker-only operation-order and component-local-scale bugs; this pilot is revalidated under the original 5e-15 and 1e-18 limits but is not part of the protected holdout claim.",
            "implementation_or_validation_repairs": [
                "full-H exact-zero delta callback",
                "full-H target-major output grammar",
                "real returned-coefficient diagnostic recomputation",
                "checker row-wise action quadratic-form order",
                "checker payload-wide vector/matrix scale",
            ],
        },
        "comparison_profiles": {
            "execution_header_policy_snapshot": {
                "maximum_absolute_error": "5e-14",
                "maximum_scale_relative_error": "0x1p-46",
                "diagnostic_report_comparison":
                    "absolute error <=4*DBL_EPSILON or relative error <=2e-10",
            },
            "final_holdout_profile_used_for_pilot_revalidation": {
                "maximum_absolute_error": "5e-14",
                "maximum_scale_relative_error": "5e-15",
                "diagnostic_report_comparison":
                    "absolute error <=1e-18 or relative error <=2e-10",
            },
        },
        "frozen_v3_binding": {
            "artifact_manifest_sha256": sha256_file(FROZEN_MANIFEST),
            "artifact_set_digest":
                json.loads(FROZEN_MANIFEST.read_text())["artifact_set_digest"],
            "fixture_registry_sha256": sha256_file(REGISTRY),
            "fixture_module_sha256": sha256_file(FIXTURE_MODULE),
        },
        "execution_binding": {
            "pilot_policy_snapshot_sha256": pilot_policy_sha,
            "driver_source_sha256": execution_driver_sha256,
            "driver_reconstruction": {
                "base_path": "test/power_lanczos_p6c1_cartesian_driver.c",
                "rule": "replace only the embedded NUMERICAL_POLICY_SHA value with pilot_policy_snapshot_sha256",
                "reconstructed_sha256": execution_driver_sha256,
            },
            "gevp_source_sha256": sha256_file(GEVP_SOURCE),
            "gevp_header_sha256": sha256_file(GEVP_HEADER),
            "executable_sha256": execution_binary_sha256,
            "compiler_identity": compiler_identity,
            "strict_fp": ["-fno-fast-math", "-ffp-contract=off"],
        },
        "input": {
            "raw_sha256": sha256_bytes(input_bytes),
            "gzip_path": PILOT_INPUT_GZIP.name,
            "gzip_sha256": sha256_bytes(generated[PILOT_INPUT_GZIP]),
        },
        "runs": [{
            "name": name,
            "world_size": world_size,
            "raw_sha256": parsed_runs[name].raw_sha256,
            "gzip_path": {
                "serial": PILOT_SERIAL_GZIP.name,
                "mpi2": PILOT_MPI2_GZIP.name,
                "mpi4": PILOT_MPI4_GZIP.name,
            }[name],
            "gzip_sha256": sha256_bytes(generated[{
                "serial": PILOT_SERIAL_GZIP,
                "mpi2": PILOT_MPI2_GZIP,
                "mpi4": PILOT_MPI4_GZIP,
            }[name]]),
            "normalized_payload_sha256": normalized_run_sha256(
                parsed_runs[name]),
            "verdict": "PASS",
        } for name, (world_size, _) in run_payloads.items()],
        "case_records": records_by_run["serial"],
        "serial_mpi2_mpi4_semantically_identical": True,
        "overall_verdict": "PASS_DEVELOPMENT_PILOT_ONLY",
    }
    payload = canonical_json(evidence)
    PILOT_EVIDENCE.write_bytes(payload)
    require(PILOT_EVIDENCE.read_bytes() == payload,
            "pilot evidence round-trip")
    print(f"PASS wrote focused pilot evidence sha256={sha256_bytes(payload)} "
          f"normalized_sha256={next(iter(normalized_shas))}")


def run_preflight(command: Sequence[str], expected_world_size: int) -> None:
    require(command, "preflight command")
    completed = subprocess.run(command, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, check=False)
    require(completed.returncode == 0,
            "preflight command failed rc=" + str(completed.returncode) +
            " stderr=" + completed.stderr.decode("utf-8", "replace"))
    validate_preflight_payload(completed.stdout, expected_world_size)


def validate_full_payload(payload: bytes,
                          expected_world_size: int
                          ) -> tuple[ParsedRun, list[dict[str, Any]]]:
    validate_numerical_policy(accepted_required=True)
    run = parse_output(payload, preflight=False)
    require(run.world_size == expected_world_size,
            "full Cartesian world-size binding")
    validate_model_blocks(run)
    _, logical_rows = input_payload(preflight=False)
    require(len(run.cases) == len(logical_rows) == CASE_COUNT,
            "full Cartesian case census")
    fixture = load_fixture_module()
    records = []
    for case, logical in zip(run.cases, logical_rows):
        oracle = fixture.exact_krylov_oracle(logical, include_explicit=True)
        records.append(validate_case(case, logical, oracle, run.world_size))
    return run, records


def metric_value(record: dict[str, Any], key: str) -> float:
    token = str(record[key])
    return (float.fromhex(token) if token.startswith(("0x", "-0x"))
            else float(Decimal(token)))


def maximum_metric(records: Sequence[dict[str, Any]],
                   key: str) -> dict[str, Any]:
    require(records, f"nonempty metric records {key}")
    record = max(records, key=lambda item: metric_value(item, key))
    result = {"value": record[key], "fixture_id": record["fixture_id"]}
    payload_key = key.removesuffix("_error") + "_payload"
    if payload_key in record:
        result["payload"] = record[payload_key]
    return result


def axis_census(records: Sequence[dict[str, Any]]) -> dict[str, int]:
    census: dict[str, int] = {}
    for record in records:
        descriptor = record["descriptor"]
        key = "/".join((
            str(descriptor["model"]), str(descriptor["arithmetic"]),
            f"order{descriptor['order']}", f"NQP{descriptor['nqp_full']}",
            str(descriptor["correlator_profile"]),
            str(descriptor["factor_product"]),
            str(descriptor["projection_weight"]),
            str(descriptor["support_case"]),
        ))
        census[key] = census.get(key, 0) + 1
    return dict(sorted(census.items()))


def marginal_census(records: Sequence[dict[str, Any]]) -> dict[str, dict[str, int]]:
    axes = {
        "model": "model", "arithmetic": "arithmetic", "order": "order",
        "nqp_full": "nqp_full", "correlator_profile": "correlator_profile",
        "factor_product": "factor_product",
        "projection_weight": "projection_weight",
        "support_case": "support_case",
    }
    result: dict[str, dict[str, int]] = {}
    for name, field in axes.items():
        counts: dict[str, int] = {}
        for record in records:
            value = str(record["descriptor"][field])
            counts[value] = counts.get(value, 0) + 1
        result[name] = dict(sorted(counts.items()))
    return result


def load_deterministic_gzip(path: Path, label: str) -> bytes:
    compressed = path.read_bytes()
    payload = gzip.decompress(compressed)
    require(deterministic_gzip(payload) == compressed,
            f"{label} deterministic gzip")
    return payload


def mpi_identity(mpiexec_path: Path) -> str:
    completed = subprocess.run(
        [str(mpiexec_path), "--version"], stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False)
    require(completed.returncode == 0, "MPI launcher identity")
    lines = [line for line in
             completed.stdout.decode("utf-8", "replace").splitlines()
             if line.strip()]
    require(lines, "MPI launcher identity output")
    match = re.fullmatch(
        r"(?:mpiexec|mpirun) \(Open MPI\) ([0-9]+(?:\.[0-9]+)+)",
        lines[0])
    require(match is not None, "MPI launcher identity grammar")
    return "Open MPI " + match.group(1)


def build_attempt_claim(executable_sha256: str,
                        pre_execution_manifest_sha256: str,
                        mpiexec_path: Path) -> dict[str, Any]:
    require(SHA256.fullmatch(executable_sha256) is not None and
            SHA256.fullmatch(pre_execution_manifest_sha256) is not None,
            "attempt claim SHA grammar")
    return {
        "schema_version": 2,
        "attempt_id": "power_lanczos_p6c1_cartesian_attempt_v2",
        "status": "execution_claimed_holdout_observation_pending",
        "production_authorized": False,
        "calibration_or_confirmation_result": False,
        "pre_execution_manifest_sha256": pre_execution_manifest_sha256,
        "numerical_policy_sha256": sha256_file(NUMERICAL_POLICY),
        "executable_sha256": executable_sha256,
        "runtime_binding": {
            "mpiexec_sha256": sha256_file(mpiexec_path),
            "mpi_identity": mpi_identity(mpiexec_path),
            "thread_environment": STRICT_THREAD_ENVIRONMENT,
        },
        "planned_run_order": ["serial", "mpi2", "mpi4"],
        "planned_world_sizes": [1, 2, 4],
        "planned_fixture_count": CASE_COUNT,
        "pilot_regression_count": len(PREFLIGHT_FIXTURE_IDS),
        "protected_holdout_count": CASE_COUNT - len(PREFLIGHT_FIXTURE_IDS),
        "rerun_authorized": False,
        "failure_rule":
            "Any failed or interrupted attempt consumes this version; preserve partial artifacts and require independent review plus a new version before rerun.",
        "success_next_phase": "p6_c1_independent_result_review",
        "failure_next_phase": "p6_c1_independent_failure_review",
    }


def validate_attempt_claim(executable_sha256: str,
                           pre_execution_manifest_sha256: str,
                           mpiexec_path: Path) -> dict[str, Any]:
    payload = ATTEMPT_CLAIM.read_bytes()
    claim = json.loads(payload)
    expected = build_attempt_claim(
        executable_sha256, pre_execution_manifest_sha256, mpiexec_path)
    require(claim == expected and payload == canonical_json(expected),
            "attempt claim canonical exact regeneration")
    return claim


def build_cartesian_evidence(
        input_path: Path, run_paths: dict[str, tuple[int, Path]],
        executable_sha256: str,
        pre_execution_manifest_sha256: str,
        attempt_claim_sha256: str) -> dict[str, Any]:
    require(SHA256.fullmatch(executable_sha256) is not None and
            SHA256.fullmatch(pre_execution_manifest_sha256) is not None and
            SHA256.fullmatch(attempt_claim_sha256) is not None,
            "Cartesian evidence SHA grammar")
    expected_input, _ = input_payload(preflight=False)
    input_payload_bytes = load_deterministic_gzip(input_path, "input")
    require(input_payload_bytes == expected_input,
            "official input canonical byte regeneration")
    parsed_runs: dict[str, ParsedRun] = {}
    records_by_run: dict[str, list[dict[str, Any]]] = {}
    normalized_shas = set()
    run_records = []
    for name, (world_size, path) in run_paths.items():
        raw = load_deterministic_gzip(path, name)
        run, records = validate_full_payload(raw, world_size)
        parsed_runs[name] = run
        records_by_run[name] = records
        normalized = normalized_run_sha256(run)
        normalized_shas.add(normalized)
        partition_census: dict[str, list[dict[str, int]]] = {}
        for case in run.cases:
            nqp_key = str(case["nqp"])
            ranges = [{"rank": rank, "begin": begin, "end": end}
                      for rank, begin, end in case["partitions"]]
            if nqp_key in partition_census:
                require(partition_census[nqp_key] == ranges,
                        f"{name}: stable NQP partition {nqp_key}")
            else:
                partition_census[nqp_key] = ranges
        run_records.append({
            "name": name,
            "world_size": world_size,
            "command_shape": (
                ["<cartesian-driver>", "<canonical-input>"]
                if world_size == 1 else
                ["<mpiexec>", "-n", str(world_size),
                 "<cartesian-driver>", "<canonical-input>"]),
            "case_count": len(records),
            "registered_nqp_full_values": [1, 4, 8],
            "nqp_full_range": [1, 8],
            "qp_partition_rule":
                "rank r owns [floor(NQP*r/P),floor(NQP*(r+1)/P))",
            "qp_partition_census": {
                key: partition_census[key] for key in ("1", "4", "8")},
            "negative_probe": {
                "id": "ZERO_PACKED_P4",
                "status": run.negative[0],
                "valid": run.negative[1],
                "verdict": "PASS",
            },
            "raw_sha256": run.raw_sha256,
            "gzip_path": path.name,
            "gzip_sha256": sha256_file(path),
            "normalized_payload_sha256": normalized,
            "verdict": "PASS",
        })
    require(list(run_paths) == ["serial", "mpi2", "mpi4"] and
            [item[0] for item in run_paths.values()] == [1, 2, 4],
            "official run order/world-size census")
    require(len(normalized_shas) == 1,
            "official serial/MPI2/MPI4 normalized payload identity")
    serial_records = records_by_run["serial"]
    pilot_ids = set(PREFLIGHT_FIXTURE_IDS)
    classified_records = []
    for record in serial_records:
        item = dict(record)
        item["classification"] = (
            "pilot_regression" if record["fixture_id"] in pilot_ids
            else "holdout")
        classified_records.append(item)
    pilot_records = [record for record in classified_records
                     if record["classification"] == "pilot_regression"]
    holdout_records = [record for record in classified_records
                       if record["classification"] == "holdout"]
    require({record["fixture_id"] for record in pilot_records} ==
            pilot_ids and len(pilot_records) == len(PREFLIGHT_FIXTURE_IDS) and
            len(holdout_records) == 2007,
            "pilot/holdout classification")
    metric_keys = (
        "maximum_vector_binary64_abs_error",
        "maximum_vector_binary64_scale_relative_error",
        "maximum_vector_exact_interval_excess",
        "maximum_vector_exact_interval_scale_relative_excess",
        "maximum_matrix_binary64_abs_error",
        "maximum_matrix_binary64_scale_relative_error",
        "maximum_matrix_exact_interval_excess",
        "maximum_matrix_exact_interval_scale_relative_excess",
        "maximum_normwise_backward_error",
        "maximum_independent_normwise_backward_error",
        "maximum_raw_action_relative_residual",
        "maximum_independent_raw_action_relative_residual",
        "maximum_observable_reference_abs_error",
    )
    summaries = {}
    for name, records in (("pilot_regression", pilot_records),
                          ("holdout", holdout_records),
                          ("overall", classified_records)):
        summaries[name] = {
            "fixture_count": len(records),
            "axis_census": axis_census(records),
            "marginal_census": marginal_census(records),
            "maxima": {key: maximum_metric(records, key)
                       for key in metric_keys},
            "verdict": "PASS",
        }
    production_model_hashes = validate_model_blocks(parsed_runs["serial"])
    registered_model_hashes = registered_full_hamiltonian_shas()
    full_hamiltonians = []
    for key, block in sorted(parsed_runs["serial"].models.items()):
        model, arithmetic = key
        dense_rows = [{
            "target_configuration_bits": target,
            "source_configuration_bits": source,
            "value": complex_record(value),
        } for target, source, value in block["rows"]]
        full_hamiltonians.append({
            "model": model,
            "arithmetic": arithmetic,
            "state_count": block["state_count"],
            "term_count": block["term_count"],
            "operator_count": block["operator_count"],
            "production_dense_sha256":
                production_model_hashes[f"{model}:{arithmetic}"],
            "registered_nonzero_oracle_sha256":
                registered_model_hashes[key],
            "production_dense_rows": dense_rows,
            "production_dense_rows_sha256":
                sha256_bytes(canonical_json(dense_rows)),
            "verdict": "PASS",
        })
    return {
        "schema_version": 2,
        "evidence_id": "power_lanczos_p6c1_cartesian_evidence_v2",
        "status": "observed_method_validation_result",
        "calibration_or_confirmation_result": False,
        "production_authorized": False,
        "pre_execution_manifest_sha256":
            pre_execution_manifest_sha256,
        "attempt_claim_sha256": attempt_claim_sha256,
        "frozen_v3_binding": {
            "artifact_manifest_sha256": TRUSTED_FROZEN_MANIFEST_SHA256,
            "artifact_set_digest": TRUSTED_FROZEN_ARTIFACT_SET_DIGEST,
            "fixture_registry_sha256": TRUSTED_REGISTRY_SHA256,
            "fixture_module_sha256": TRUSTED_FIXTURE_MODULE_SHA256,
        },
        "numerical_policy_sha256": sha256_file(NUMERICAL_POLICY),
        "execution_binding": {
            "executable_sha256": executable_sha256,
            "driver_source_sha256": sha256_file(DRIVER_SOURCE),
            "checker_source_sha256": sha256_file(CHECKER_SOURCE),
            "gevp_source_sha256": sha256_file(GEVP_SOURCE),
            "gevp_header_sha256": sha256_file(GEVP_HEADER),
            "strict_fp": ["-fno-fast-math", "-ffp-contract=off"],
        },
        "input": {
            "path": input_path.name,
            "raw_sha256": sha256_bytes(input_payload_bytes),
            "gzip_sha256": sha256_file(input_path),
        },
        "runs": run_records,
        "serial_mpi2_mpi4_semantically_identical": True,
        "normalized_payload_sha256": next(iter(normalized_shas)),
        "full_hamiltonians": full_hamiltonians,
        "summaries": summaries,
        "case_records": classified_records,
        "overall_verdict": "PASS",
    }


def post_result_manifest(evidence: dict[str, Any],
                         attempt_path: Path,
                         input_path: Path,
                         run_paths: dict[str, tuple[int, Path]],
                         evidence_path: Path) -> dict[str, Any]:
    members = [{"role": "attempt_claim", "path": attempt_path.name,
                "sha256": sha256_file(attempt_path)},
               {"role": "input", "path": input_path.name,
                "sha256": sha256_file(input_path)}]
    members.extend({"role": name, "path": path.name,
                    "sha256": sha256_file(path)}
                   for name, (_, path) in run_paths.items())
    members.append({"role": "evidence", "path": evidence_path.name,
                    "sha256": sha256_file(evidence_path)})
    return {
        "schema_version": 2,
        "manifest_id": "power_lanczos_p6c1_cartesian_post_result_v2",
        "status": "observed_pass",
        "production_authorized": False,
        "calibration_or_confirmation_result": False,
        "pre_execution_manifest_sha256":
            evidence["pre_execution_manifest_sha256"],
        "attempt_claim_sha256": evidence["attempt_claim_sha256"],
        "pilot_regression_count": 9,
        "holdout_count": 2007,
        "overall_count": CASE_COUNT,
        "pilot_regression_verdict":
            evidence["summaries"]["pilot_regression"]["verdict"],
        "holdout_verdict": evidence["summaries"]["holdout"]["verdict"],
        "overall_verdict": evidence["overall_verdict"],
        "hash_algorithm": "sha256",
        "digest_algorithm":
            "sha256(canonical_json_ordered_member_role_path_sha256)",
        "members": members,
        "artifact_set_digest": sha256_bytes(canonical_json(members)),
        "next_phase": "p6_c1_independent_result_review",
    }


def run_control_validator(
        mode: str, cartesian_executable_path: Path,
        preflight_executable_path: Path,
        expected_pre_execution_manifest_sha256: str,
        gevp_unit_path: Path, asan_gevp_unit_path: Path,
        production_binary_path: Path, mpiexec_path: Path,
        artifact_paths: dict[str, Path] | None = None) -> None:
    command = [
        sys.executable, str(CONTROL_VALIDATOR), mode,
        "--cartesian-executable", str(cartesian_executable_path),
        "--preflight-executable", str(preflight_executable_path),
        "--expected-pre-execution-manifest-sha256",
        expected_pre_execution_manifest_sha256,
        "--gevp-unit", str(gevp_unit_path),
        "--asan-gevp-unit", str(asan_gevp_unit_path),
        "--production-binary", str(production_binary_path),
        "--mpiexec", str(mpiexec_path)]
    for name, path in (artifact_paths or {}).items():
        command.extend((f"--{name.replace('_', '-')}", str(path)))
    completed = subprocess.run(
        command,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    require(completed.returncode == 0,
            f"control validator {mode} failed: " +
            completed.stderr.decode("utf-8", "replace"))


def validate_post_result_artifacts(
        executable_path: Path,
        expected_pre_execution_manifest_sha256: str,
        mpiexec_path: Path) -> None:
    validate_numerical_policy(accepted_required=True)
    require(PRE_EXECUTION_MANIFEST.is_file() and ATTEMPT_CLAIM.is_file() and
            DEFAULT_INPUT.is_file() and DEFAULT_SERIAL.is_file() and
            DEFAULT_MPI2.is_file() and DEFAULT_MPI4.is_file() and
            DEFAULT_EVIDENCE.is_file() and DEFAULT_MANIFEST.is_file(),
            "complete post-result artifact census")
    require(SHA256.fullmatch(expected_pre_execution_manifest_sha256) is not
            None and sha256_file(PRE_EXECUTION_MANIFEST) ==
                expected_pre_execution_manifest_sha256,
            "external pre-execution manifest trust anchor")
    executable_sha = sha256_file(executable_path)
    validate_attempt_claim(
        executable_sha, expected_pre_execution_manifest_sha256,
        mpiexec_path)
    run_paths = {
        "serial": (1, DEFAULT_SERIAL),
        "mpi2": (2, DEFAULT_MPI2),
        "mpi4": (4, DEFAULT_MPI4),
    }
    expected_evidence = build_cartesian_evidence(
        DEFAULT_INPUT, run_paths, executable_sha,
        sha256_file(PRE_EXECUTION_MANIFEST), sha256_file(ATTEMPT_CLAIM))
    expected_evidence_bytes = canonical_json(expected_evidence)
    require(DEFAULT_EVIDENCE.read_bytes() == expected_evidence_bytes,
            "post-result evidence independent regeneration")
    expected_manifest = post_result_manifest(
        expected_evidence, ATTEMPT_CLAIM, DEFAULT_INPUT, run_paths,
        DEFAULT_EVIDENCE)
    expected_manifest_bytes = canonical_json(expected_manifest)
    require(DEFAULT_MANIFEST.read_bytes() == expected_manifest_bytes,
            "post-result manifest independent regeneration")
    print(f"PASS post-result artifacts evidence="
          f"{sha256_bytes(expected_evidence_bytes)} manifest="
          f"{sha256_bytes(expected_manifest_bytes)}")


def official_artifact_paths() -> tuple[Path, ...]:
    return (
        ATTEMPT_CLAIM, DEFAULT_INPUT, DEFAULT_SERIAL, DEFAULT_MPI2,
        DEFAULT_MPI4, DEFAULT_EVIDENCE, DEFAULT_MANIFEST, FAILURE_STDOUT,
        FAILURE_STDERR, FAILURE_MANIFEST, FAILURE_CANDIDATE_EVIDENCE,
        FAILURE_CANDIDATE_MANIFEST,
    )


def build_failure_manifest(
        executable_sha256: str, failed_stage: str,
        failed_world_size: int | None, process_started: bool,
        return_code: int | None, error_type: str,
        completed_runs: Sequence[str],
        run_paths: dict[str, tuple[int, Path]],
        stored_failure_stdout_sha256: str,
        stored_failure_stderr_sha256: str) -> dict[str, Any]:
    members = [{"role": "attempt_claim", "path": ATTEMPT_CLAIM.name,
                "sha256": sha256_file(ATTEMPT_CLAIM)}]
    if DEFAULT_INPUT.exists():
        members.append({"role": (
                            "unvalidated_input" if
                            failed_stage == "input_generation" else "input"),
                        "path": DEFAULT_INPUT.name,
                        "sha256": sha256_file(DEFAULT_INPUT)})
    for name in ("serial", "mpi2", "mpi4"):
        path = run_paths[name][1]
        if path.exists():
            members.append({
                "role": name if name in completed_runs else
                    f"unvalidated_{name}",
                "path": path.name,
                "sha256": sha256_file(path)})
    if DEFAULT_EVIDENCE.is_file():
        members.append({"role": "unpublished_pass_evidence",
                        "path": DEFAULT_EVIDENCE.name,
                        "sha256": sha256_file(DEFAULT_EVIDENCE)})
    if FAILURE_CANDIDATE_EVIDENCE.is_file():
        members.append({
            "role": "unvalidated_candidate_evidence",
            "path": FAILURE_CANDIDATE_EVIDENCE.name,
            "sha256": sha256_file(FAILURE_CANDIDATE_EVIDENCE)})
    if FAILURE_CANDIDATE_MANIFEST.is_file():
        members.append({
            "role": "unvalidated_candidate_manifest",
            "path": FAILURE_CANDIDATE_MANIFEST.name,
            "sha256": sha256_file(FAILURE_CANDIDATE_MANIFEST)})
    members.extend((
        {"role": "failure_stdout", "path": FAILURE_STDOUT.name,
         "sha256": sha256_file(FAILURE_STDOUT)},
        {"role": "failure_stderr", "path": FAILURE_STDERR.name,
         "sha256": sha256_file(FAILURE_STDERR)},
    ))
    return {
        "schema_version": 2,
        "manifest_id": "power_lanczos_p6c1_cartesian_failure_v2",
        "status": "observed_failed_attempt",
        "production_authorized": False,
        "calibration_or_confirmation_result": False,
        "pre_execution_manifest_sha256":
            sha256_file(PRE_EXECUTION_MANIFEST),
        "attempt_claim_sha256": sha256_file(ATTEMPT_CLAIM),
        "numerical_policy_sha256": sha256_file(NUMERICAL_POLICY),
        "executable_sha256": executable_sha256,
        "failed_stage": failed_stage,
        "failed_world_size": failed_world_size,
        "process_started": process_started,
        "return_code": return_code,
        "error_type": error_type,
        "completed_runs": list(completed_runs),
        "holdout_observation_state": (
            "partial_or_complete_result_observed" if process_started or
            completed_runs else "not_observed_but_attempt_consumed"),
        "stored_failure_stdout_sha256": stored_failure_stdout_sha256,
        "stored_failure_stderr_sha256": stored_failure_stderr_sha256,
        "stored_failure_streams_sanitized": True,
        "rerun_authorized": False,
        "hash_algorithm": "sha256",
        "digest_algorithm":
            "sha256(canonical_json_ordered_member_role_path_sha256)",
        "members": members,
        "artifact_set_digest": sha256_bytes(canonical_json(members)),
        "next_phase": "p6_c1_independent_failure_review",
    }


def execute_cartesian(
        executable_path: Path, preflight_executable_path: Path,
        mpiexec_path: Path,
        gevp_unit_path: Path, asan_gevp_unit_path: Path,
        production_binary_path: Path,
        expected_pre_execution_manifest_sha256: str) -> None:
    validate_numerical_policy(accepted_required=True)
    require(SHA256.fullmatch(expected_pre_execution_manifest_sha256) is not
            None and sha256_file(PRE_EXECUTION_MANIFEST) ==
                expected_pre_execution_manifest_sha256,
            "external pre-execution manifest trust anchor")
    run_control_validator(
        "--pre-exec", executable_path, preflight_executable_path,
        expected_pre_execution_manifest_sha256, gevp_unit_path,
        asan_gevp_unit_path, production_binary_path, mpiexec_path)
    executable_sha = sha256_file(executable_path)
    existing = [path.name for path in official_artifact_paths()
                if path.exists()]
    require(not existing,
            "official attempt already claimed or artifact exists: " +
            ",".join(existing))
    pre_execution_sha = sha256_file(PRE_EXECUTION_MANIFEST)
    claim_bytes = canonical_json(build_attempt_claim(
        executable_sha, pre_execution_sha, mpiexec_path))
    exclusive_write(ATTEMPT_CLAIM, claim_bytes)
    run_paths = {
        "serial": (1, DEFAULT_SERIAL),
        "mpi2": (2, DEFAULT_MPI2),
        "mpi4": (4, DEFAULT_MPI4),
    }
    completed_runs: list[str] = []
    current_stage = "input_generation"
    current_world_size: int | None = None
    process_started = False
    completed_process: subprocess.CompletedProcess[bytes] | None = None
    evidence_bytes: bytes | None = None
    manifest_bytes: bytes | None = None
    try:
        input_bytes, _ = input_payload(preflight=False)
        exclusive_write(DEFAULT_INPUT, deterministic_gzip(input_bytes))
        scratch_root = PROJECT / "tmp"
        scratch_root.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="p6c1-cartesian-",
                                         dir=scratch_root) as directory:
            raw_input = Path(directory) / "input.txt"
            raw_input.write_bytes(input_bytes)
            commands = {
                "serial": [str(executable_path), str(raw_input)],
                "mpi2": [str(mpiexec_path), "-n", "2",
                         str(executable_path), str(raw_input)],
                "mpi4": [str(mpiexec_path), "-n", "4",
                         str(executable_path), str(raw_input)],
            }
            environment = dict(os.environ)
            environment.update(STRICT_THREAD_ENVIRONMENT)
            for name in ("serial", "mpi2", "mpi4"):
                current_stage = name
                current_world_size = run_paths[name][0]
                process_started = True
                completed_process = subprocess.run(
                    commands[name], stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, check=False, env=environment)
                require(completed_process.returncode == 0 and
                        completed_process.stderr == b"",
                        f"official {name} process result")
                validate_full_payload(completed_process.stdout,
                                      current_world_size)
                exclusive_write(
                    run_paths[name][1],
                    deterministic_gzip(completed_process.stdout))
                completed_runs.append(name)
                process_started = False
                completed_process = None
        current_stage = "evidence_generation"
        current_world_size = None
        evidence = build_cartesian_evidence(
            DEFAULT_INPUT, run_paths, executable_sha, pre_execution_sha,
            sha256_file(ATTEMPT_CLAIM))
        evidence_bytes = canonical_json(evidence)
        current_stage = "post_result_validation"
        with tempfile.TemporaryDirectory(prefix="p6c1-post-validation-",
                                         dir=scratch_root) as directory:
            staged_evidence = Path(directory) / DEFAULT_EVIDENCE.name
            staged_manifest = Path(directory) / DEFAULT_MANIFEST.name
            staged_evidence.write_bytes(evidence_bytes)
            manifest = post_result_manifest(
                evidence, ATTEMPT_CLAIM, DEFAULT_INPUT, run_paths,
                staged_evidence)
            manifest_bytes = canonical_json(manifest)
            staged_manifest.write_bytes(manifest_bytes)
            run_control_validator(
                "--post-result", executable_path,
                preflight_executable_path,
                expected_pre_execution_manifest_sha256, gevp_unit_path,
                asan_gevp_unit_path, production_binary_path, mpiexec_path,
                artifact_paths={
                    "attempt": ATTEMPT_CLAIM,
                    "input": DEFAULT_INPUT,
                    "serial": DEFAULT_SERIAL,
                    "mpi2": DEFAULT_MPI2,
                    "mpi4": DEFAULT_MPI4,
                    "evidence": staged_evidence,
                    "manifest": staged_manifest,
                })
            current_stage = "result_publication"
            exclusive_publish(staged_evidence, DEFAULT_EVIDENCE)
            exclusive_publish(staged_manifest, DEFAULT_MANIFEST)
    except Exception as error:
        raw_stdout = (completed_process.stdout
                      if completed_process is not None else b"")
        raw_stderr = (completed_process.stderr
                      if completed_process is not None else
                      (type(error).__name__ + ": " + str(error) + "\n").encode())
        stored_stdout = sanitized_failure_stderr(
            raw_stdout, executable_path, mpiexec_path,
            (preflight_executable_path, gevp_unit_path, asan_gevp_unit_path,
             production_binary_path))
        stored_stderr = sanitized_failure_stderr(
            raw_stderr, executable_path, mpiexec_path,
            (preflight_executable_path, gevp_unit_path, asan_gevp_unit_path,
             production_binary_path))
        if evidence_bytes is not None:
            exclusive_write(FAILURE_CANDIDATE_EVIDENCE,
                            deterministic_gzip(evidence_bytes))
        if manifest_bytes is not None:
            exclusive_write(FAILURE_CANDIDATE_MANIFEST,
                            deterministic_gzip(manifest_bytes))
        exclusive_write(FAILURE_STDOUT, deterministic_gzip(stored_stdout))
        exclusive_write(FAILURE_STDERR, deterministic_gzip(stored_stderr))
        failure = build_failure_manifest(
            executable_sha, current_stage, current_world_size,
            process_started, (completed_process.returncode
                              if completed_process is not None else None),
            type(error).__name__, completed_runs, run_paths,
            sha256_bytes(stored_stdout), sha256_bytes(stored_stderr))
        failure_bytes = canonical_json(failure)
        exclusive_write(FAILURE_MANIFEST, failure_bytes)
        run_control_validator(
            "--failure-result", executable_path,
            preflight_executable_path,
            expected_pre_execution_manifest_sha256, gevp_unit_path,
            asan_gevp_unit_path, production_binary_path, mpiexec_path)
        raise ValidationError(
            "official one-time attempt failed; rerun forbidden; failure=" +
            sha256_bytes(failure_bytes)) from error
    print(f"PASS P6-C1 Cartesian evidence={sha256_bytes(evidence_bytes)} "
          f"manifest={sha256_bytes(manifest_bytes)}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    validate_policy_parser = subparsers.add_parser("validate-policy")
    validate_policy_parser.add_argument(
        "--accepted", action="store_true",
        help="require the accepted_result_unobserved execution boundary")
    prepare_parser = subparsers.add_parser("prepare-input")
    prepare_parser.add_argument("path", type=Path)
    preflight_parser = subparsers.add_parser("prepare-preflight")
    preflight_parser.add_argument("path", type=Path)
    run_preflight_parser = subparsers.add_parser("run-preflight")
    run_preflight_parser.add_argument("--expected-world-size", type=int,
                                      required=True)
    run_preflight_parser.add_argument("command_line", nargs=argparse.REMAINDER)
    pilot_parser = subparsers.add_parser("write-focused-pilot-evidence")
    pilot_parser.add_argument("--input", type=Path, required=True)
    pilot_parser.add_argument("--serial", type=Path, required=True)
    pilot_parser.add_argument("--mpi2", type=Path, required=True)
    pilot_parser.add_argument("--mpi4", type=Path, required=True)
    pilot_parser.add_argument("--execution-driver-sha256", required=True)
    pilot_parser.add_argument("--execution-binary-sha256", required=True)
    pilot_parser.add_argument("--compiler-identity", required=True)
    post_parser = subparsers.add_parser("validate-post-result-artifacts")
    post_parser.add_argument("--executable", type=Path, required=True)
    post_parser.add_argument(
        "--expected-pre-execution-manifest-sha256", required=True)
    post_parser.add_argument("--mpiexec", type=Path, required=True)
    execute_parser = subparsers.add_parser("execute-cartesian")
    execute_parser.add_argument("--executable", type=Path, required=True)
    execute_parser.add_argument("--preflight-executable", type=Path,
                                required=True)
    execute_parser.add_argument("--mpiexec", type=Path, required=True)
    execute_parser.add_argument("--gevp-unit", type=Path, required=True)
    execute_parser.add_argument("--asan-gevp-unit", type=Path, required=True)
    execute_parser.add_argument("--production-binary", type=Path,
                                required=True)
    execute_parser.add_argument(
        "--expected-pre-execution-manifest-sha256", required=True)
    arguments = parser.parse_args()
    try:
        if arguments.command == "validate-policy":
            policy = validate_numerical_policy(
                accepted_required=arguments.accepted)
            print(f"PASS numerical policy {policy['policy_id']} "
                  f"status={policy['status']} sha256={sha256_file(NUMERICAL_POLICY)}")
        elif arguments.command == "prepare-input":
            prepare_input(arguments.path)
        elif arguments.command == "prepare-preflight":
            prepare_preflight_input(arguments.path)
        elif arguments.command == "run-preflight":
            command_line = arguments.command_line
            if command_line and command_line[0] == "--":
                command_line = command_line[1:]
            run_preflight(command_line, arguments.expected_world_size)
        elif arguments.command == "write-focused-pilot-evidence":
            write_focused_pilot_evidence(
                arguments.input, arguments.serial, arguments.mpi2,
                arguments.mpi4, arguments.execution_driver_sha256,
                arguments.execution_binary_sha256,
                arguments.compiler_identity)
        elif arguments.command == "validate-post-result-artifacts":
            validate_post_result_artifacts(
                arguments.executable.resolve(),
                arguments.expected_pre_execution_manifest_sha256,
                arguments.mpiexec.resolve())
        elif arguments.command == "execute-cartesian":
            execute_cartesian(
                arguments.executable.resolve(),
                arguments.preflight_executable.resolve(),
                arguments.mpiexec.resolve(),
                arguments.gevp_unit.resolve(),
                arguments.asan_gevp_unit.resolve(),
                arguments.production_binary.resolve(),
                arguments.expected_pre_execution_manifest_sha256)
        else:  # pragma: no cover - argparse enforces the command census.
            raise ValidationError("unknown command")
    except (OSError, ValueError, ValidationError, json.JSONDecodeError) as error:
        print(f"FAIL {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
