from __future__ import print_function

import copy
import math
import sys

import generate_bounded_krylov_p4c_evidence as evidence


def close(actual, expected, tolerance=1.0e-13):
    return abs(actual - expected) <= tolerance * (1.0 + abs(expected))


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError("negative gate was accepted: {}".format(label))


def main():
    if len(sys.argv) != 2:
        raise AssertionError(
            "usage: runtest_bounded_krylov_p4c_evidence.py POLICY")
    policy = evidence.validate_policy_contract(
        evidence.read_json(sys.argv[1]))
    bad_policy = copy.deepcopy(policy)
    bad_policy["scope"]["cache_capacity_bytes"] = [0, 4096]
    expect_failure(lambda: evidence.validate_policy_contract(bad_policy),
                   "official cache grid")
    bad_policy = copy.deepcopy(policy)
    bad_policy["guide"]["candidate_family"] = "r2"
    expect_failure(lambda: evidence.validate_policy_contract(bad_policy),
                   "r3-only candidate")
    bad_policy = copy.deepcopy(policy)
    bad_policy["scope"]["sample_count_by_site"]["8"] = 5000
    expect_failure(lambda: evidence.validate_policy_contract(bad_policy),
                   "resource sample count bounds")
    rows = [
        {"configuration": 1,
         "values": [1.0 + 0.0j, 2.0 + 0.0j,
                    3.0 + 0.0j, 4.0 + 0.0j]},
        {"configuration": 2,
         "values": [1.0 + 0.0j, 2.0 + 0.0j,
                    3.0 + 0.0j, 4.0 + 0.0j]},
        {"configuration": 3,
         "values": [1.0 + 0.0j, 2.0 + 0.0j,
                    3.0 + 0.0j, 4.0 + 0.0j]},
    ]
    result = evidence.exact_iid_statistics(rows, policy)
    if result["sector_size"] != 3:
        raise AssertionError("unexpected exact sector size")
    scales = result["krylov_scales"]
    scaled = [[scales[index] * row["values"][index]
               for index in evidence.ORDERS] for row in rows]
    independent_k00 = sum(row[0].conjugate() * row[1] for row in scaled)
    encoded_k00 = result["matrices_unnormalized"]["K"][0][0]
    if not close(complex(*encoded_k00), independent_k00):
        raise AssertionError("K_00 index convention mismatch")
    for guide in result["guides"].values():
        if guide["r"] != 3 or not guide["strict_positive_support"]:
            raise AssertionError("P4-C guide contract mismatch")
        if not guide["finite_weight"]:
            raise AssertionError("finite-weight gate mismatch")
        if not guide["all_entry_budgets_pass"]:
            raise AssertionError("constant exact fixture should have zero SE")
        if guide["maximum_required_entry_se_to_budget_ratio"] != 0.0:
            raise AssertionError("constant exact fixture should have zero ratio")
        for entry in guide["entries"].values():
            if entry["influence_complex_variance"] != 0.0:
                raise AssertionError("constant influence variance mismatch")
            if entry["absolute_standard_error_at_4096"] != 0.0:
                raise AssertionError("constant exact-IID SE mismatch")
    fit = evidence.log_linear_fit(
        [(4, math.exp(0.4)), (6, math.exp(0.6)),
         (8, math.exp(0.8))], 16)
    if not close(fit["prediction"], math.exp(1.6), 1.0e-12):
        raise AssertionError("log-linear extrapolation mismatch")
    print("bounded Krylov P4-C evidence math: PASS")


if __name__ == "__main__":
    main()
