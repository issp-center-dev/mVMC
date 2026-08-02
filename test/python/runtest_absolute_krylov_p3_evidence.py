from __future__ import print_function

import math
import copy
import sys

import generate_absolute_krylov_p3_evidence as evidence


def close(actual, expected, tolerance=1.0e-13):
    return abs(actual - expected) <= tolerance * (1.0 + abs(expected))


def main():
    if len(sys.argv) != 2:
        raise AssertionError(
            "usage: runtest_absolute_krylov_p3_evidence.py POLICY"
        )
    policy = evidence.read_json(sys.argv[1])
    evidence.validate_policy_contract(policy)
    bad_policy = copy.deepcopy(policy)
    bad_policy["official_scaling_grid"]["sample_count_by_site"].pop("8")
    try:
        evidence.validate_policy_contract(bad_policy)
    except AssertionError:
        pass
    else:
        raise AssertionError("incomplete policy contract was accepted")
    rows = [
        {"configuration": 1,
         "values": [1.0 + 0.0j, 1.0 + 0.0j,
                    1.0 + 1.0j, 2.0 - 1.0j]},
        {"configuration": 2,
         "values": [1.0 + 0.0j, -1.0 + 0.0j,
                    2.0 - 1.0j, -2.0 + 1.0j]},
    ]
    result = evidence.exact_case_statistics(rows, policy)
    if result["sector_size"] != 2 or len(result["guides"]) != 12:
        raise AssertionError("unexpected exact evidence dimensions")
    scales = result["krylov_scales"]
    scaled = [[scales[index] * row["values"][index]
               for index in evidence.ORDERS] for row in rows]
    independent_k00 = sum(row[0].conjugate() * row[1] for row in scaled)
    encoded_k00 = result["matrices_unnormalized"]["K"][0][0]
    if not close(complex(*encoded_k00), independent_k00):
        raise AssertionError("K_00 index convention mismatch")
    independent_b22 = sum(row[3].conjugate() * row[3] for row in scaled)
    encoded_b22 = result["matrices_unnormalized"]["B"][2][2]
    if not close(complex(*encoded_b22), independent_b22):
        raise AssertionError("B_22 index convention mismatch")
    for guide in result["guides"].values():
        if not guide["strict_positive_support"]:
            raise AssertionError("positive floor did not close support")
        if not (0.0 < guide["t_ess_fraction"] <= 1.0 + 1.0e-13):
            raise AssertionError("invalid T-ESS fraction")
        for entry in guide["entries"].values():
            if (entry["normalized_complex_variance"] < 0.0 or
                    entry["absolute_standard_error_at_4096"] < 0.0):
                raise AssertionError("invalid integrand variance/SE")
    fit = evidence.log_linear_fit(
        [(4, math.exp(0.4)), (6, math.exp(0.6)),
         (8, math.exp(0.8))], 16)
    if not close(fit["prediction"], math.exp(1.6), 1.0e-12):
        raise AssertionError("log-linear extrapolation mismatch")
    zero = evidence.extrapolate_nonnegative([(4, 0), (6, 0), (8, 0)], 16)
    if zero["prediction"] != 0.0 or zero["model"] != "structural_zero":
        raise AssertionError("structural zero extrapolation mismatch")
    if evidence.parse_number("fixture-id") != "fixture-id":
        raise AssertionError("profile string parser mismatch")
    malformed = []
    for case in (
            "electronic_real", "electronic_complex", "spin_real",
            "spin_complex"):
        malformed.extend([
            "BEGIN {} rank_count 1 empty_qp_ranks 0".format(case),
            "ROW 0 1 0 2 0 3 0 4 0",
            "END {}".format(case),
        ])
    try:
        evidence.parse_ed_output("\n".join(malformed) + "\n")
    except AssertionError:
        pass
    else:
        raise AssertionError("incomplete exact sector was accepted")
    print("absolute Krylov P3 evidence math: PASS")


if __name__ == "__main__":
    main()
