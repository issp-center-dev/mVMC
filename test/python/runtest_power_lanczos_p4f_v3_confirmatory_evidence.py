#!/usr/bin/env python3

import copy
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import generate_power_lanczos_p4f_v3_confirmatory_evidence as evidence


def prefix(energy_se, error=0.0):
    return {"energy_se": energy_se, "exact_energy_error": error}


def case(seed, errors):
    return {
        "confirmatory_identity": {"seed": seed},
        "dimensions": [{
            "dimension": 3,
            "exact": {"energy": -1.0},
            "prefixes": [prefix(value) for value in errors],
            "se_convergence": {
                "theil_sen_log_se_log_n_slope": -0.5,
            },
        }],
    }


def main():
    solver_policy = {
        "exact_gate": {"energy_tolerance_relative": 1.0e-11},
    }
    group = [
        case(f"seed-{index}", errors)
        for index, errors in enumerate((
            (4.0, 2.0, 1.0),
            (1.0, 0.9, 1.1),
            (2.0, 1.0, 0.8),
            (3.0, 1.5, 1.0),
        ))
    ]
    result = evidence.pool_group(group, 3, solver_policy)
    if not result["passed"]:
        raise AssertionError("pooled convergence should pass")
    if result["per_seed_diagnostics"][1][
        "reported_final_to_initial_ratio"
    ] <= 0.8:
        raise AssertionError("fixture must include a noisy per-seed ratio")
    if result["final_to_initial_ratio"] > 0.8 or not (
        result["theil_sen_log_se_log_n_slope"] < 0.0
    ):
        raise AssertionError("pooled formula mismatch")

    mixed = copy.deepcopy(group)
    mixed[0]["dimensions"][0]["prefixes"] = [
        prefix(1.0), prefix(0.0), prefix(0.5)
    ]
    result = evidence.pool_group(mixed, 3, solver_policy)
    if result["passed"] or result["seed_zero_patterns_pass"]:
        raise AssertionError("mixed-zero seed must STOP")

    all_zero = [case(f"zero-{index}", (0.0, 0.0, 0.0))
                for index in range(4)]
    result = evidence.pool_group(all_zero, 3, solver_policy)
    if not result["passed"] or result["zero_rule"] != "all_zero":
        raise AssertionError("exact all-zero pooled fixture should pass")
    all_zero[0]["dimensions"][0]["prefixes"][1][
        "exact_energy_error"
    ] = 1.0e-5
    result = evidence.pool_group(all_zero, 3, solver_policy)
    if result["passed"]:
        raise AssertionError("inexact all-zero pooled fixture must STOP")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
