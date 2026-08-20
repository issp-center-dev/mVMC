#!/usr/bin/env python3
"""Generate deterministic P5-B downstream-sensitivity evidence."""

import argparse
import hashlib
import json
import math
import os
import tempfile

import numpy as np


EXPECTED_FIXTURES = [
    "real_dimension_2_support_bridge",
    "complex_dimension_3_krylov",
]
EXPECTED_DIRECTIONS = [
    "coefficient_real",
    "coefficient_imaginary",
    "overlap_hermitian",
    "hamiltonian_hermitian",
    "squared_hamiltonian_hermitian",
    "observable_hermitian",
    "combined_all",
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_policy(path):
    with open(path, "r", encoding="utf-8") as handle:
        policy = json.load(handle)
    require(policy.get("schema") == "power_lanczos_p5b_sensitivity_policy_v1",
            "unexpected policy schema")
    require(policy.get("phase") == "P5-B" and
            policy.get("testing_only") is True and
            policy.get("production_authorized") is False,
            "invalid P5-B authorization fields")
    require(policy.get("fixtures") == EXPECTED_FIXTURES,
            "fixture inventory mismatch")
    require(policy.get("directions") == EXPECTED_DIRECTIONS,
            "direction inventory mismatch")
    scales = policy.get("perturbation_scales")
    require(isinstance(scales, list) and len(scales) >= 3,
            "perturbation scale grid is incomplete")
    require(all(isinstance(value, (int, float)) and
                math.isfinite(value) and value > 0.0 for value in scales),
            "invalid perturbation scale")
    require(scales == sorted(set(scales)),
            "perturbation scales must be strictly increasing")
    targets = policy.get("downstream_absolute_targets", {})
    for key in ("energy", "full_support_variance",
                "quadratic_observable"):
        require(isinstance(targets.get(key), (int, float)) and
                math.isfinite(targets[key]) and targets[key] > 0.0,
                "invalid downstream target {}".format(key))
    factor = policy.get("candidate_envelope_safety_factor")
    require(isinstance(factor, (int, float)) and
            math.isfinite(factor) and 0.0 < factor < 1.0,
            "invalid envelope safety factor")
    required = policy.get("required", {})
    require(required.get("finite_all_cases") is True and
            required.get("positive_overlap") is True,
            "missing hard validity requirements")
    for key in ("solver_residual_maximum",
                "baseline_identity_delta_maximum",
                "global_phase_invariance_maximum"):
        require(isinstance(required.get(key), (int, float)) and
                math.isfinite(required[key]) and required[key] > 0.0,
                "invalid required threshold {}".format(key))
    return policy


def hermitian_frobenius(matrix):
    return float(np.linalg.norm(matrix, ord="fro"))


def normalize_direction(matrix):
    norm = hermitian_frobenius(matrix)
    require(math.isfinite(norm) and norm > 0.0,
            "zero/nonfinite perturbation direction")
    return matrix / norm


def canonical_coefficient(coefficient, overlap):
    norm = float(np.real(np.vdot(coefficient, overlap @ coefficient)))
    require(math.isfinite(norm) and norm > 0.0,
            "non-positive coefficient overlap norm")
    result = coefficient / math.sqrt(norm)
    pivot = int(np.argmax(np.abs(result)))
    phase = result[pivot] / abs(result[pivot])
    return result / phase


def solve_generalized(overlap, hamiltonian):
    overlap_values, overlap_vectors = np.linalg.eigh(overlap)
    require(float(np.min(overlap_values)) > 1.0e-12,
            "candidate overlap is not positive definite")
    inverse_sqrt = ((overlap_vectors *
                     (1.0 / np.sqrt(overlap_values))) @
                    overlap_vectors.conj().T)
    reduced = inverse_sqrt @ hamiltonian @ inverse_sqrt
    reduced = 0.5 * (reduced + reduced.conj().T)
    values, vectors = np.linalg.eigh(reduced)
    coefficient = canonical_coefficient(inverse_sqrt @ vectors[:, 0],
                                        overlap)
    residual = np.linalg.norm(
        hamiltonian @ coefficient - values[0] * overlap @ coefficient)
    scale = max(1.0,
                np.linalg.norm(hamiltonian @ coefficient),
                abs(values[0]) * np.linalg.norm(overlap @ coefficient))
    return coefficient, float(values[0]), float(residual / scale)


def metrics(model, coefficient):
    overlap, hamiltonian, squared, observable = model
    norm = float(np.real(np.vdot(coefficient, overlap @ coefficient)))
    require(math.isfinite(norm) and norm > 0.0,
            "non-positive quadratic norm")
    energy = float(np.real(np.vdot(coefficient,
                                   hamiltonian @ coefficient)) / norm)
    second = float(np.real(np.vdot(coefficient,
                                   squared @ coefficient)) / norm)
    variance = second - energy * energy
    tolerance = 2048.0 * np.finfo(float).eps * max(
        1.0, abs(second) + energy * energy)
    require(variance >= -tolerance,
            "negative full-support variance")
    variance = max(0.0, variance)
    observable_value = float(np.real(np.vdot(
        coefficient, observable @ coefficient)) / norm)
    values = (norm, energy, second, variance, observable_value)
    require(all(math.isfinite(value) for value in values),
            "nonfinite quadratic metric")
    return {
        "norm": norm,
        "energy": energy,
        "hamiltonian_second_moment": second,
        "full_support_variance": variance,
        "quadratic_observable": observable_value,
    }


def projective_distance(baseline, candidate):
    denominator = np.linalg.norm(baseline) * np.linalg.norm(candidate)
    require(float(denominator) > 0.0, "zero coefficient norm")
    cosine = min(1.0, max(0.0,
                         float(abs(np.vdot(baseline, candidate)) /
                               denominator)))
    return math.sqrt(max(0.0, 1.0 - cosine * cosine))


def relative_matrix_error(baseline, candidate):
    return float(np.linalg.norm(candidate - baseline, ord="fro") /
                 max(np.linalg.norm(baseline, ord="fro"),
                     np.finfo(float).tiny))


def downstream_error(baseline, candidate):
    return {
        "energy": abs(candidate["energy"] - baseline["energy"]),
        "full_support_variance": abs(
            candidate["full_support_variance"] -
            baseline["full_support_variance"]),
        "quadratic_observable": abs(
            candidate["quadratic_observable"] -
            baseline["quadratic_observable"]),
    }


def compare(baseline_model, baseline_coefficient,
            candidate_model, candidate_coefficient):
    baseline = metrics(baseline_model, baseline_coefficient)
    coefficient_only_metrics = metrics(baseline_model,
                                       candidate_coefficient)
    matrix_only_metrics = metrics(candidate_model, baseline_coefficient)
    combined_metrics = metrics(candidate_model, candidate_coefficient)
    names = ("overlap", "hamiltonian", "hamiltonian_squared",
             "observable")
    return {
        "coefficient_projective_distance": projective_distance(
            baseline_coefficient, candidate_coefficient),
        "matrix_relative_frobenius_error": {
            name: relative_matrix_error(baseline_matrix, candidate_matrix)
            for name, baseline_matrix, candidate_matrix in zip(
                names, baseline_model, candidate_model)
        },
        "baseline": baseline,
        "coefficient_only": downstream_error(
            baseline, coefficient_only_metrics),
        "matrix_only": downstream_error(baseline, matrix_only_metrics),
        "combined": downstream_error(baseline, combined_metrics),
    }


def fixture_dimension2():
    overlap = np.eye(2, dtype=complex)
    hamiltonian = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    squared = np.eye(2, dtype=complex)
    observable = np.diag([1.0, -1.0]).astype(complex)
    coefficient, energy, residual = solve_generalized(overlap, hamiltonian)
    return (overlap, hamiltonian, squared, observable), coefficient, {
        "exact_energy": energy,
        "solver_residual": residual,
    }


def fixture_dimension3():
    hamiltonian = np.array([
        [-1.1, 0.24 + 0.08j, -0.07j, 0.03],
        [0.24 - 0.08j, -0.2, 0.19 + 0.04j, -0.11j],
        [0.07j, 0.19 - 0.04j, 0.65, 0.16 + 0.05j],
        [0.03, 0.11j, 0.16 - 0.05j, 1.4],
    ], dtype=complex)
    psi = np.array([1.0, 0.45 + 0.2j, -0.3 + 0.25j,
                    0.18 - 0.12j], dtype=complex)
    psi /= np.linalg.norm(psi)
    basis = np.column_stack((psi, hamiltonian @ psi,
                             hamiltonian @ hamiltonian @ psi))
    overlap = basis.conj().T @ basis
    hamiltonian_matrix = basis.conj().T @ hamiltonian @ basis
    h_basis = hamiltonian @ basis
    squared = h_basis.conj().T @ h_basis
    physical_observable = np.diag([1.0, -0.5, 0.25, -0.75])
    observable = basis.conj().T @ physical_observable @ basis
    model = tuple(0.5 * (matrix + matrix.conj().T)
                  for matrix in (overlap, hamiltonian_matrix,
                                 squared, observable))
    coefficient, energy, residual = solve_generalized(model[0], model[1])
    return model, coefficient, {
        "exact_energy": energy,
        "solver_residual": residual,
    }


def matrix_directions(dimension):
    identity = normalize_direction(np.eye(dimension, dtype=complex))
    observable = np.zeros((dimension, dimension), dtype=complex)
    for index in range(dimension):
        observable[index, index] = (-1.0) ** index
    if dimension > 1:
        observable[0, 1] = 0.3 + 0.2j
        observable[1, 0] = observable[0, 1].conjugate()
    return identity, normalize_direction(observable)


def coefficient_direction(dimension, imaginary):
    direction = np.array([complex(index + 1, dimension - index)
                          for index in range(dimension)], dtype=complex)
    if imaginary:
        direction *= 1j
    return direction / np.linalg.norm(direction)


def perturb_case(model, baseline_coefficient, direction, scale):
    overlap, hamiltonian, squared, observable = (
        matrix.copy() for matrix in model)
    dimension = len(baseline_coefficient)
    identity_direction, observable_direction = matrix_directions(dimension)
    coefficient = baseline_coefficient.copy()
    if direction == "coefficient_real":
        coefficient += scale * coefficient_direction(dimension, False)
    elif direction == "coefficient_imaginary":
        coefficient += scale * coefficient_direction(dimension, True)
    elif direction == "overlap_hermitian":
        overlap += scale * hermitian_frobenius(overlap) * identity_direction
        coefficient, _, residual = solve_generalized(overlap, hamiltonian)
    elif direction == "hamiltonian_hermitian":
        hamiltonian += (scale * hermitian_frobenius(hamiltonian) *
                        identity_direction)
        coefficient, _, residual = solve_generalized(overlap, hamiltonian)
    elif direction == "squared_hamiltonian_hermitian":
        squared += scale * hermitian_frobenius(squared) * identity_direction
    elif direction == "observable_hermitian":
        observable += (scale * hermitian_frobenius(observable) *
                       observable_direction)
    elif direction == "combined_all":
        overlap += scale * hermitian_frobenius(overlap) * identity_direction
        hamiltonian += (scale * hermitian_frobenius(hamiltonian) *
                        identity_direction)
        squared += scale * hermitian_frobenius(squared) * identity_direction
        observable += (scale * hermitian_frobenius(observable) *
                       observable_direction)
        coefficient, _, residual = solve_generalized(overlap, hamiltonian)
    else:
        raise ValueError("unknown perturbation direction {}".format(
            direction))
    if direction not in ("overlap_hermitian", "hamiltonian_hermitian",
                         "combined_all"):
        residual = float(np.linalg.norm(
            hamiltonian @ coefficient -
            (np.vdot(coefficient, hamiltonian @ coefficient) /
             np.vdot(coefficient, overlap @ coefficient)) *
            overlap @ coefficient) /
            max(1.0, np.linalg.norm(hamiltonian @ coefficient)))
    candidate = tuple(0.5 * (matrix + matrix.conj().T)
                      for matrix in (overlap, hamiltonian,
                                     squared, observable))
    return candidate, coefficient, residual


def generate(policy, policy_path):
    fixtures = {
        "real_dimension_2_support_bridge": fixture_dimension2(),
        "complex_dimension_3_krylov": fixture_dimension3(),
    }
    required = policy["required"]
    targets = policy["downstream_absolute_targets"]
    cases = []
    identity_maximum = 0.0
    phase_maximum = 0.0
    baseline_records = {}
    for name in EXPECTED_FIXTURES:
        model, coefficient, metadata = fixtures[name]
        identity = compare(model, coefficient, model, coefficient)
        identity_maximum = max(
            identity_maximum,
            max(identity[category][quantity]
                for category in ("coefficient_only", "matrix_only",
                                 "combined")
                for quantity in ("energy", "full_support_variance",
                                 "quadratic_observable")))
        phase_coefficient = coefficient * np.exp(0.731j) * 3.25
        phase_metrics = metrics(model, phase_coefficient)
        baseline_metrics = metrics(model, coefficient)
        phase_maximum = max(
            phase_maximum,
            max(abs(phase_metrics[key] - baseline_metrics[key])
                for key in ("energy", "full_support_variance",
                            "quadratic_observable")))
        baseline_records[name] = {
            **metadata,
            "metrics": baseline_metrics,
        }
        for direction in EXPECTED_DIRECTIONS:
            for scale in policy["perturbation_scales"]:
                candidate_model, candidate_coefficient, residual = (
                    perturb_case(model, coefficient, direction, scale))
                comparison = compare(model, coefficient, candidate_model,
                                     candidate_coefficient)
                require(residual <= required["solver_residual_maximum"] or
                        direction.startswith("coefficient_") or
                        direction in ("squared_hamiltonian_hermitian",
                                      "observable_hermitian"),
                        "solver residual exceeds policy")
                cases.append({
                    "fixture": name,
                    "direction": direction,
                    "scale": scale,
                    "solver_residual": residual,
                    **comparison,
                })
    require(identity_maximum <=
            required["baseline_identity_delta_maximum"],
            "baseline identity gate failed")
    require(phase_maximum <= required["global_phase_invariance_maximum"],
            "global phase gate failed")

    envelopes = []
    largest_passing_scale = None
    observed_failure = False
    for scale in policy["perturbation_scales"]:
        selected = [case for case in cases if case["scale"] == scale]
        maxima = {}
        for quantity in ("energy", "full_support_variance",
                         "quadratic_observable"):
            maxima[quantity] = max(
                case[category][quantity]
                for case in selected
                for category in ("coefficient_only", "matrix_only",
                                 "combined"))
        passed = all(maxima[key] <= targets[key] for key in maxima)
        require(not (observed_failure and passed),
                "perturbation envelope is non-monotonic")
        if passed:
            largest_passing_scale = scale
        else:
            observed_failure = True
        envelopes.append({
            "scale": scale,
            "maximum_downstream_absolute_error": maxima,
            "passes_all_targets": passed,
        })
    require(largest_passing_scale is not None,
            "no perturbation scale satisfies downstream targets")
    with open(policy_path, "rb") as handle:
        policy_bytes = handle.read()
    return {
        "schema": "power_lanczos_p5b_sensitivity_evidence_v1",
        "phase": "P5-B",
        "testing_only": True,
        "production_authorized": False,
        "policy_sha256": hashlib.sha256(policy_bytes).hexdigest(),
        "case_count": len(cases),
        "baseline": baseline_records,
        "identity_maximum": identity_maximum,
        "global_phase_invariance_maximum": phase_maximum,
        "cases": cases,
        "envelopes": envelopes,
        "largest_passing_perturbation_scale": largest_passing_scale,
        "p5c_candidate_input_envelope": (
            largest_passing_scale *
            policy["candidate_envelope_safety_factor"]),
        "decision": "PASS",
    }


def atomic_write_json(path, payload):
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=".p5b-sensitivity-", suffix=".json", dir=directory)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    arguments = parser.parse_args()
    policy = load_policy(arguments.policy)
    payload = generate(policy, arguments.policy)
    atomic_write_json(arguments.output, payload)
    print("P5-B sensitivity evidence: {} ({} cases, envelope {})".format(
        payload["decision"], payload["case_count"],
        payload["p5c_candidate_input_envelope"]))


if __name__ == "__main__":
    main()
