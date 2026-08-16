#!/usr/bin/env python3

import argparse
import hashlib
import json
import math
import pathlib


POLICY_ID = "power_lanczos_zero_support_gate_classification_v1"


class DuplicateKeyError(ValueError):
    pass


def reject_duplicate_keys(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise DuplicateKeyError(f"duplicate JSON key: {key}")
        value[key] = item
    return value


def load_json(path):
    return json.loads(pathlib.Path(path).read_text(encoding="utf-8"),
                      object_pairs_hook=reject_duplicate_keys)


def sha256_path(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite_number(value):
    return (not isinstance(value, bool) and isinstance(value, (int, float))
            and math.isfinite(value))


def item_label(case, dimension=None, suffix=None):
    label = str(case.get("trace", "<unknown-trace>"))
    if dimension is not None:
        label += f":d{dimension}"
    if suffix:
        label += f":{suffix}"
    return label


def validate_policy(policy):
    if not isinstance(policy, dict) or policy.get("schema") != 1 or \
            policy.get("policy_id") != POLICY_ID:
        raise ValueError("invalid gate-classification policy")
    artifact = policy.get("artifact_contract", {})
    exact = policy.get("deterministic_exact_correctness", {})
    required_positive = (
        (artifact, "trace_count"),
        (artifact, "seed_count"),
        (artifact, "sample_count_per_trace"),
        (artifact, "pooled_gate_count"),
        (exact, "maximum_reconstruction_error"),
        (exact, "maximum_relative_residual"),
        (exact, "maximum_relative_antihermitian_residual"),
        (exact, "negative_variance_relative_tolerance"),
        (exact, "minimum_retained_rank"),
    )
    for section, key in required_positive:
        value = section.get(key)
        if not finite_number(value) or value <= 0:
            raise ValueError(f"invalid policy value: {key}")
    return policy


def validate_external_artifact(validation, evidence, policy, findings):
    artifact = policy["artifact_contract"]
    evidence_value = evidence if isinstance(evidence, dict) else {}
    archive = validation.get("archive", {}) if isinstance(validation, dict) \
        else {}
    validator = validation.get("validator", {}) \
        if isinstance(validation, dict) else {}
    required_true = {
        "remote_completed_marker_present": validation.get(
            "remote_completed_marker_present") if isinstance(
                validation, dict) else None,
        "outer_checksum_verified": archive.get("outer_checksum_verified"),
        "artifact_ledger_verified": archive.get("artifact_ledger_verified"),
        "confirmatory_ledger_verified": archive.get(
            "confirmatory_ledger_verified"),
        "validator_valid": validator.get("valid"),
    }
    if not isinstance(validation, dict) or validation.get("schema") != \
            artifact["validation_schema"]:
        findings.append("artifact validation record schema mismatch")
    for name, value in required_true.items():
        if value is not True:
            findings.append(f"artifact validation check failed: {name}")
    for name in ("unsafe_path_count", "duplicate_member_count",
                 "special_file_count"):
        if archive.get(name) != 0:
            findings.append(f"artifact archive check failed: {name}")
    if validator.get("decision") != evidence_value.get("decision"):
        findings.append("validator/evidence frozen decision mismatch")
    if validator.get("trace_count") != artifact["trace_count"]:
        findings.append("validator trace census mismatch")
    stopped = evidence_value.get("stopped_traces")
    if isinstance(stopped, list) and \
            validator.get("mandatory_trace_stop_count") != len(stopped):
        findings.append("validator stopped-trace census mismatch")
    pooled = evidence_value.get("pooled_se_convergence_gates")
    if isinstance(pooled, list):
        pooled_stop_count = sum(
            isinstance(item, dict) and item.get("passed") is not True
            for item in pooled)
        if validator.get("pooled_gate_stop_count") != pooled_stop_count:
            findings.append("validator pooled-gate census mismatch")


def frozen_dimension_mandatory(dimension, policy):
    prefixes = dimension.get("prefixes")
    if not isinstance(prefixes, list):
        return None
    exact = dimension.get("exact", {})
    anti = exact.get("antihermitian_residual")
    if not finite_number(anti):
        return None
    return all((
        dimension.get("solver_pass") is True,
        dimension.get("cutoff_scan", {}).get("full_pass") is True,
        dimension.get("cutoff_scan", {}).get("exact_pass") is True,
        all(item.get("exact_energy_pass") is True for item in prefixes),
        all(item.get("rank_stability_pass") is True for item in prefixes),
        dimension.get("coefficient_gate", {}).get("passed") is True,
        anti <= policy["deterministic_exact_correctness"]
            ["maximum_relative_antihermitian_residual"],
    ))


def frozen_trace_mandatory(case, policy):
    identity = case.get("confirmatory_identity", {})
    reconstruction = case.get("reconstruction", {})
    diagnostics = case.get("matrix_diagnostics")
    dimensions = case.get("dimensions")
    if not isinstance(diagnostics, list) or \
            {item.get("family") for item in diagnostics
             if isinstance(item, dict)} != {"K", "S", "B"} or \
            len(diagnostics) != 3 or not isinstance(dimensions, list):
        return None
    dimension_values = [frozen_dimension_mandatory(item, policy)
                        for item in dimensions]
    if any(value is None for value in dimension_values):
        return None
    return all((
        identity.get("passed") is True,
        reconstruction.get("passed") is True,
        all(item.get("passed") is True for item in diagnostics),
        all(dimension_values),
    ))


def validate_evidence_structure(evidence, policy, findings):
    artifact = policy["artifact_contract"]
    if not isinstance(evidence, dict) or evidence.get("schema") != 1:
        findings.append("evidence schema mismatch")
        return []
    if evidence.get("policy_id") not in policy[
            "applies_to_evidence_policy_ids"]:
        findings.append("evidence policy is outside classification scope")
    if evidence.get("decision") not in ("GO", "STOP"):
        findings.append("invalid frozen evidence decision")
    if evidence.get("identity_grid_pass") is not True:
        findings.append("confirmatory identity grid did not pass")
    if evidence.get("trace_count") != artifact["trace_count"]:
        findings.append("evidence trace_count mismatch")
    cases = evidence.get("cases")
    if not isinstance(cases, list) or len(cases) != artifact["trace_count"]:
        findings.append("evidence case census mismatch")
        return []

    trace_names = set()
    grid_keys = set()
    seeds = set()
    stopped = []
    expected_checks = set(artifact["identity_checks"])
    required_identity = set(artifact["artifact_identity_checks"])
    expected_dimensions = set(artifact["dimensions"])
    expected_prefixes = artifact["prefix_sample_counts"]
    for case in cases:
        if not isinstance(case, dict):
            findings.append("non-object case record")
            continue
        trace = case.get("trace")
        if not isinstance(trace, str) or not trace or "/" in trace or \
                "\\" in trace or trace in trace_names:
            findings.append(f"unsafe or duplicate trace: {trace}")
        else:
            trace_names.add(trace)
        identity = case.get("confirmatory_identity", {})
        checks = identity.get("checks", {}) if isinstance(identity, dict) \
            else {}
        if not isinstance(checks, dict) or set(checks) != expected_checks or \
                any(not isinstance(value, bool)
                    for value in checks.values()):
            findings.append(f"identity check schema mismatch: {trace}")
        else:
            for name in required_identity:
                if checks.get(name) is not True:
                    findings.append(
                        f"artifact identity check failed: {trace}:{name}")
            if identity.get("passed") is not all(checks.values()):
                findings.append(f"identity aggregate mismatch: {trace}")
        seed = str(identity.get("seed", "")).lower()
        site = case.get("site_count")
        qp = case.get("qp_total")
        if not seed or identity.get("site_count") != site or \
                identity.get("qp_total") != qp or \
                site not in artifact["site_counts"] or \
                qp not in artifact["qp_totals"]:
            findings.append(f"case identity mismatch: {trace}")
        key = (seed, site, qp)
        if key in grid_keys:
            findings.append(f"duplicate seed/case grid entry: {trace}")
        grid_keys.add(key)
        seeds.add(seed)
        if case.get("sample_count") != artifact["sample_count_per_trace"]:
            findings.append(f"sample count mismatch: {trace}")

        dimensions = case.get("dimensions")
        if not isinstance(dimensions, list) or \
                {item.get("dimension") for item in dimensions
                 if isinstance(item, dict)} != expected_dimensions or \
                len(dimensions) != len(expected_dimensions):
            findings.append(f"dimension census mismatch: {trace}")
        else:
            for dimension in dimensions:
                prefixes = dimension.get("prefixes")
                if not isinstance(prefixes, list) or \
                        [item.get("sample_count") for item in prefixes
                         if isinstance(item, dict)] != expected_prefixes:
                    findings.append(
                        f"prefix census mismatch: {item_label(case, dimension.get('dimension'))}")
        mandatory = case.get("mandatory_gate_pass")
        decision = case.get("decision")
        if not isinstance(mandatory, bool) or \
                decision != ("GO" if mandatory else "STOP"):
            findings.append(f"recorded trace decision mismatch: {trace}")
        recomputed_mandatory = frozen_trace_mandatory(case, policy)
        if recomputed_mandatory is None:
            findings.append(f"frozen mandatory gate schema mismatch: {trace}")
        elif mandatory is not recomputed_mandatory:
            findings.append(f"frozen mandatory gate value mismatch: {trace}")
        if mandatory is False:
            stopped.append(trace)

    expected_grid_size = (artifact["seed_count"] *
                          len(artifact["site_counts"]) *
                          len(artifact["qp_totals"]))
    if len(seeds) != artifact["seed_count"] or \
            len(grid_keys) != expected_grid_size:
        findings.append("confirmatory seed/case grid mismatch")
    if evidence.get("stopped_traces") != stopped:
        findings.append("frozen stopped-trace list mismatch")
    pooled = evidence.get("pooled_se_convergence_gates")
    if not isinstance(pooled, list) or \
            len(pooled) != artifact["pooled_gate_count"] or \
            any(not isinstance(item, dict) or
                not isinstance(item.get("passed"), bool) for item in pooled):
        findings.append("pooled statistical gate census/schema mismatch")
    else:
        pooled_stops = [item for item in pooled if item["passed"] is not True]
        expected_decision = "STOP" if stopped or pooled_stops else "GO"
        if evidence.get("decision") != expected_decision:
            findings.append("frozen aggregate decision mismatch")
    return cases


def deterministic_exact_findings(cases, policy):
    exact_policy = policy["deterministic_exact_correctness"]
    failures = []
    for case in cases:
        trace = case.get("trace", "<unknown-trace>")
        reconstruction = case.get("reconstruction", {})
        if reconstruction.get("passed") is not True or \
                reconstruction.get("block_sum_pass") is not True or \
                reconstruction.get("corrected_solver_match") is not True:
            failures.append(f"{trace}: exact reconstruction mismatch")
        for name in ("maximum_denominator_error", "maximum_entry_error",
                     "maximum_exact_error", "source_antihermitian_error"):
            value = reconstruction.get(name)
            if not finite_number(value) or value < 0 or \
                    value > exact_policy["maximum_reconstruction_error"]:
                failures.append(f"{trace}: {name} exceeds exact tolerance")
        for dimension in case.get("dimensions", []):
            number = dimension.get("dimension")
            label = item_label(case, number)
            exact = dimension.get("exact", {})
            if exact.get("status") != "ok" or exact.get("valid") != 1:
                failures.append(f"{label}: exact solver invalid")
                continue
            for name in ("energy", "energy_squared", "variance", "residual",
                         "antihermitian_residual"):
                if not finite_number(exact.get(name)):
                    failures.append(f"{label}: nonfinite exact {name}")
            residual = exact.get("residual")
            if finite_number(residual) and (residual < 0 or residual >
                    exact_policy["maximum_relative_residual"]):
                failures.append(f"{label}: exact residual exceeds tolerance")
            anti = exact.get("antihermitian_residual")
            if finite_number(anti) and (anti < 0 or anti >
                    exact_policy["maximum_relative_antihermitian_residual"]):
                failures.append(
                    f"{label}: exact anti-Hermitian residual exceeds tolerance")
            variance = exact.get("variance")
            energy_squared = exact.get("energy_squared")
            if finite_number(variance) and finite_number(energy_squared):
                tolerance = exact_policy[
                    "negative_variance_relative_tolerance"] * max(
                        1.0, abs(energy_squared))
                if variance < -tolerance:
                    failures.append(f"{label}: exact variance is negative")
            rank = exact.get("retained_rank")
            if not isinstance(rank, int) or isinstance(rank, bool) or \
                    rank < exact_policy["minimum_retained_rank"]:
                failures.append(f"{label}: exact retained rank is invalid")
            coefficient = exact.get("coefficient")
            if not isinstance(coefficient, list) or len(coefficient) != number \
                    or any(not isinstance(value, dict) or
                           not finite_number(value.get("real")) or
                           not finite_number(value.get("imag"))
                           for value in coefficient):
                failures.append(f"{label}: exact coefficient is invalid")
            if exact_policy["require_exact_cutoff_scan_pass"] and \
                    dimension.get("cutoff_scan", {}).get("exact_pass") is \
                    not True:
                failures.append(f"{label}: exact cutoff scan failed")
    return failures


def statistical_observations(cases, evidence, policy):
    prefixes = policy["artifact_contract"]["prefix_sample_counts"]
    final_count = prefixes[-1]
    observations = {
        "p4s_combined_failures": [],
        "matrix_diagnostic_failures": [],
        "combined_solver_failures": [],
        "stochastic_cutoff_failures": [],
        "early_prefix_exact_energy_failures": [],
        "final_prefix_exact_energy_failures": [],
        "rank_stability_failures": [],
        "per_trace_se_convergence_failures": [],
        "coefficient_gate_failures": [],
        "pooled_se_convergence_failures": [],
    }
    for case in cases:
        checks = case["confirmatory_identity"]["checks"]
        if checks["p4s"] is not True:
            observations["p4s_combined_failures"].append(case["trace"])
        for diagnostic in case.get("matrix_diagnostics", []):
            if diagnostic.get("passed") is not True:
                observations["matrix_diagnostic_failures"].append(
                    item_label(case, suffix=diagnostic.get("family")))
        for dimension in case["dimensions"]:
            number = dimension["dimension"]
            label = item_label(case, number)
            if dimension.get("solver_pass") is not True:
                observations["combined_solver_failures"].append(label)
            if dimension.get("cutoff_scan", {}).get("full_pass") is not True:
                observations["stochastic_cutoff_failures"].append(label)
            for prefix in dimension["prefixes"]:
                prefix_label = item_label(
                    case, number, f"n{prefix.get('sample_count')}")
                if prefix.get("exact_energy_pass") is not True:
                    key = ("final_prefix_exact_energy_failures"
                           if prefix.get("sample_count") == final_count else
                           "early_prefix_exact_energy_failures")
                    observations[key].append(prefix_label)
                if prefix.get("rank_stability_pass") is not True:
                    observations["rank_stability_failures"].append(
                        prefix_label)
            if dimension.get("se_convergence", {}).get("passed") is not True:
                observations["per_trace_se_convergence_failures"].append(
                    label)
            if dimension.get("coefficient_gate", {}).get("passed") is not True:
                observations["coefficient_gate_failures"].append(label)
    for index, pooled in enumerate(evidence["pooled_se_convergence_gates"]):
        if pooled.get("passed") is not True:
            observations["pooled_se_convergence_failures"].append(index)
    return observations


def classify(evidence, validation, policy):
    validate_policy(policy)
    artifact_findings = []
    validate_external_artifact(validation, evidence, policy,
                               artifact_findings)
    cases = validate_evidence_structure(evidence, policy, artifact_findings)
    artifact_status = "PASS" if not artifact_findings else "ARTIFACT_INVALID"
    correctness_findings = []
    observations = {}
    if artifact_status == "PASS":
        correctness_findings = deterministic_exact_findings(cases, policy)
        observations = statistical_observations(cases, evidence, policy)
    correctness_status = (
        "NOT_EVALUATED" if artifact_status != "PASS" else
        "CORRECTNESS_FAIL" if correctness_findings else "PASS")
    ambiguous = []
    if observations:
        if observations["p4s_combined_failures"]:
            ambiguous.append("p4s combined gate failed without subfield roles")
        if observations["combined_solver_failures"]:
            ambiguous.append(
                "solver_pass combined exact and stochastic scans")
    p5_authorized = (artifact_status == "PASS" and
                     correctness_status == "PASS" and not ambiguous)
    overall = (
        "STOP_ARTIFACT_INVALID" if artifact_status != "PASS" else
        "STOP_CORRECTNESS_FAIL" if correctness_status != "PASS" else
        "BLOCKED_PENDING_GATE_ROLE_DISAGGREGATION" if ambiguous else
        "CONTINUE_P5_TESTING_ONLY_RELEASE_BLOCKED")
    return {
        "schema": 1,
        "classification_policy_id": policy["policy_id"],
        "frozen_policy": {
            "policy_id": evidence.get("policy_id"),
            "decision": evidence.get("decision"),
            "preserved": True,
        },
        "artifact_integrity": {
            "status": artifact_status,
            "findings": artifact_findings,
        },
        "deterministic_exact_correctness": {
            "status": correctness_status,
            "findings": correctness_findings,
        },
        "statistical_adequacy": {
            "status": policy["statistical_adequacy"]["release_status"],
            "bug_implication": False,
            "reason": policy["statistical_adequacy"]["reason"],
            "observations": observations,
        },
        "resource_limit": {
            "status": "NOT_EVALUATED_IN_P4F_EVIDENCE",
            "bug_implication": False,
        },
        "undisaggregated_combined_gate_findings": ambiguous,
        "authorization": {
            "p5_testing_only": {
                "status": "AUTHORIZED" if p5_authorized else "BLOCKED",
                "constraints": policy["authorization"]["p5_testing_only"]
                    ["constraints"],
            },
            "p6_production": {
                "status": policy["authorization"]["p6_production"]
                    ["status"],
            },
        },
        "overall": overall,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", required=True, type=pathlib.Path)
    parser.add_argument("--artifact-validation", required=True,
                        type=pathlib.Path)
    parser.add_argument("--policy", required=True, type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path)
    args = parser.parse_args()
    evidence = load_json(args.evidence)
    validation = load_json(args.artifact_validation)
    policy = load_json(args.policy)
    result = classify(evidence, validation, policy)
    result["inputs"] = {
        "evidence_sha256": sha256_path(args.evidence),
        "artifact_validation_sha256": sha256_path(args.artifact_validation),
        "classification_policy_sha256": sha256_path(args.policy),
        "classifier_sha256": sha256_path(pathlib.Path(__file__).resolve()),
    }
    payload = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return 0 if result["overall"] == \
        "CONTINUE_P5_TESTING_ONLY_RELEASE_BLOCKED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
