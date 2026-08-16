#!/usr/bin/env python3

import copy
import pathlib

import classify_power_lanczos_p4f_v4_evidence as classifier


def build_fixture():
    cases = []
    seeds = [f"0x50344634523056{index:02x}" for index in range(1, 5)]
    for seed_index, seed in enumerate(seeds):
        for site in (4, 6, 8):
            for qp in (1, 4):
                trace = f"seed{seed_index}-L{site}-qp{qp}.log"
                dimensions = []
                for dimension in (2, 3):
                    prefixes = []
                    for sample_count, blocks in ((8192, 4), (16384, 8),
                                                 (32768, 16)):
                        early_stop = ((seed_index, site, qp, dimension,
                                       sample_count) in {
                            (0, 8, 1, 2, 8192),
                            (2, 4, 1, 2, 8192),
                        })
                        prefixes.append({
                            "sample_count": sample_count,
                            "block_count": blocks,
                            "block_length": 2048,
                            "energy_se": 0.01,
                            "exact_energy_pass": not early_stop,
                            "rank_stability_pass": True,
                        })
                    dimension_stop = any(
                        not item["exact_energy_pass"] for item in prefixes)
                    dimensions.append({
                        "dimension": dimension,
                        "exact": {
                            "status": "ok",
                            "valid": 1,
                            "retained_rank": dimension,
                            "energy": -1.0,
                            "energy_squared": 1.25,
                            "variance": 0.25,
                            "residual": 1e-14,
                            "antihermitian_residual": 1e-15,
                            "coefficient": [
                                {"real": 1.0 if index == 0 else 0.0,
                                 "imag": 0.0}
                                for index in range(dimension)
                            ],
                        },
                        "cutoff_scan": {
                            "full_pass": True,
                            "exact_pass": True,
                        },
                        "prefixes": prefixes,
                        "se_convergence": {"passed": True},
                        "coefficient_gate": {"passed": True},
                        "solver_pass": True,
                        "decision": "STOP" if dimension_stop else "GO",
                    })
                mandatory = all(item["decision"] == "GO"
                                for item in dimensions)
                checks = {
                    "driver_schema": True,
                    "session_schema": True,
                    "fixture": True,
                    "seed": True,
                    "case": True,
                    "sample_count": True,
                    "rank_count": True,
                    "cache": True,
                    "rho": True,
                    "rng_stream": True,
                    "persistent_session": True,
                    "proposal": True,
                    "partition": True,
                    "p4s": True,
                }
                cases.append({
                    "trace": trace,
                    "site_count": site,
                    "qp_total": qp,
                    "sample_count": 32768,
                    "confirmatory_identity": {
                        "seed": seed,
                        "site_count": site,
                        "qp_total": qp,
                        "checks": checks,
                        "passed": True,
                    },
                    "reconstruction": {
                        "maximum_denominator_error": 0.0,
                        "maximum_entry_error": 1e-14,
                        "maximum_exact_error": 0.0,
                        "source_antihermitian_error": 0.0,
                        "block_sum_pass": True,
                        "corrected_solver_match": True,
                        "passed": True,
                    },
                    "matrix_diagnostics": [
                        {"family": family, "passed": True}
                        for family in ("K", "S", "B")
                    ],
                    "dimensions": dimensions,
                    "mandatory_gate_pass": mandatory,
                    "decision": "GO" if mandatory else "STOP",
                })
    stopped = [case["trace"] for case in cases
               if not case["mandatory_gate_pass"]]
    evidence = {
        "schema": 1,
        "policy_id": "power_lanczos_zero_support_p4f_confirmatory_v4",
        "trace_count": 24,
        "identity_grid_pass": True,
        "cases": cases,
        "pooled_se_convergence_gates": [
            {"passed": True} for _ in range(12)
        ],
        "stopped_traces": stopped,
        "decision": "STOP" if stopped else "GO",
    }
    validation = {
        "schema": 1,
        "remote_completed_marker_present": True,
        "archive": {
            "outer_checksum_verified": True,
            "artifact_ledger_verified": True,
            "confirmatory_ledger_verified": True,
            "unsafe_path_count": 0,
            "duplicate_member_count": 0,
            "special_file_count": 0,
        },
        "validator": {
            "valid": True,
            "decision": evidence["decision"],
            "trace_count": 24,
            "mandatory_trace_stop_count": len(stopped),
            "pooled_gate_stop_count": 0,
        },
    }
    return evidence, validation


def synchronize_frozen_decision(evidence, validation):
    stopped = [case["trace"] for case in evidence["cases"]
               if not case["mandatory_gate_pass"]]
    evidence["stopped_traces"] = stopped
    evidence["decision"] = "STOP" if stopped else "GO"
    validation["validator"]["decision"] = evidence["decision"]
    validation["validator"]["mandatory_trace_stop_count"] = len(stopped)


def stop_case_for_frozen_policy(evidence, validation, case_index,
                                dimension_index=None):
    case = evidence["cases"][case_index]
    if dimension_index is not None:
        case["dimensions"][dimension_index]["decision"] = "STOP"
    case["mandatory_gate_pass"] = False
    case["decision"] = "STOP"
    synchronize_frozen_decision(evidence, validation)


def require(value, message):
    if not value:
        raise AssertionError(message)


def main():
    policy_path = pathlib.Path(__file__).resolve().parent / \
        "data/power_lanczos_gate_classification_policy.json"
    policy = classifier.load_json(policy_path)
    try:
        classifier.reject_duplicate_keys([("schema", 1), ("schema", 2)])
    except classifier.DuplicateKeyError:
        pass
    else:
        raise AssertionError("duplicate JSON key was accepted")
    evidence, validation = build_fixture()
    result = classifier.classify(evidence, validation, policy)
    require(result["frozen_policy"]["decision"] == "STOP",
            "frozen STOP was not preserved")
    require(result["deterministic_exact_correctness"]["status"] == "PASS",
            "valid exact evidence did not pass")
    require(result["statistical_adequacy"]["status"] ==
            "STATISTICALLY_INCONCLUSIVE", "statistics were overclaimed")
    require(len(result["statistical_adequacy"]["observations"]
                ["early_prefix_exact_energy_failures"]) == 2,
            "early prefix failures were not reported")
    require(result["authorization"]["p5_testing_only"]["status"] ==
            "AUTHORIZED", "diagnostic prefix failure blocked P5 testing")
    require(result["authorization"]["p6_production"]["status"].startswith(
            "BLOCKED"), "uncalibrated statistics authorized production")

    final_failure = copy.deepcopy(evidence)
    final_failure["cases"][0]["dimensions"][0]["prefixes"][2][
        "exact_energy_pass"] = False
    final_validation = copy.deepcopy(validation)
    stop_case_for_frozen_policy(final_failure, final_validation, 0, 0)
    value = classifier.classify(final_failure, final_validation, policy)
    require(value["authorization"]["p5_testing_only"]["status"] ==
            "AUTHORIZED", "final stochastic diagnostic blocked P5 testing")
    require(len(value["statistical_adequacy"]["observations"]
                ["final_prefix_exact_energy_failures"]) == 1,
            "final prefix failure was not reported")

    coefficient_failure = copy.deepcopy(evidence)
    coefficient_failure["cases"][0]["dimensions"][0]["coefficient_gate"][
        "passed"] = False
    coefficient_validation = copy.deepcopy(validation)
    stop_case_for_frozen_policy(
        coefficient_failure, coefficient_validation, 0, 0)
    value = classifier.classify(
        coefficient_failure, coefficient_validation, policy)
    require(value["authorization"]["p5_testing_only"]["status"] ==
            "AUTHORIZED", "uncalibrated coefficient diagnostic blocked P5")

    exact_failure = copy.deepcopy(evidence)
    exact_failure["cases"][0]["dimensions"][0]["exact"]["residual"] = 1e-4
    value = classifier.classify(exact_failure, validation, policy)
    require(value["overall"] == "STOP_CORRECTNESS_FAIL",
            "exact residual failure did not stop")

    reconstruction_failure = copy.deepcopy(evidence)
    reconstruction_validation = copy.deepcopy(validation)
    reconstruction_failure["cases"][0]["reconstruction"]["passed"] = False
    stop_case_for_frozen_policy(
        reconstruction_failure, reconstruction_validation, 0)
    value = classifier.classify(
        reconstruction_failure, reconstruction_validation, policy)
    require(value["overall"] == "STOP_CORRECTNESS_FAIL",
            "reconstruction failure did not stop")

    artifact_failure = copy.deepcopy(validation)
    artifact_failure["archive"]["artifact_ledger_verified"] = False
    value = classifier.classify(evidence, artifact_failure, policy)
    require(value["overall"] == "STOP_ARTIFACT_INVALID",
            "invalid ledger did not stop")

    p4s_ambiguous = copy.deepcopy(evidence)
    p4s_validation = copy.deepcopy(validation)
    case = p4s_ambiguous["cases"][1]
    case["confirmatory_identity"]["checks"]["p4s"] = False
    case["confirmatory_identity"]["passed"] = False
    stop_case_for_frozen_policy(p4s_ambiguous, p4s_validation, 1)
    value = classifier.classify(p4s_ambiguous, p4s_validation, policy)
    require(value["overall"] ==
            "BLOCKED_PENDING_GATE_ROLE_DISAGGREGATION",
            "undisaggregated P4-S failure was misclassified")
    require(value["deterministic_exact_correctness"]["status"] == "PASS",
            "ambiguous P4-S failure was mislabeled as a code bug")

    solver_ambiguous = copy.deepcopy(evidence)
    solver_validation = copy.deepcopy(validation)
    case = solver_ambiguous["cases"][1]
    case["dimensions"][0]["solver_pass"] = False
    stop_case_for_frozen_policy(solver_ambiguous, solver_validation, 1, 0)
    value = classifier.classify(solver_ambiguous, solver_validation, policy)
    require(value["overall"] ==
            "BLOCKED_PENDING_GATE_ROLE_DISAGGREGATION",
            "combined solver failure was misclassified")

    duplicate = copy.deepcopy(evidence)
    duplicate["cases"][1]["trace"] = duplicate["cases"][0]["trace"]
    value = classifier.classify(duplicate, validation, policy)
    require(value["overall"] == "STOP_ARTIFACT_INVALID",
            "duplicate trace was accepted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
