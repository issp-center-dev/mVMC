#!/usr/bin/env python3

import copy
import json
import math

import compare_power_lanczos_p4f_evidence as comparator


def fixture(decision="GO"):
    return {
        "schema": 1,
        "metadata": {
            "generated_at": "2026-08-16 12:00 JST",
            "reparse_driver_sha256": "1" * 64,
            "generator_sha256": "2" * 64,
        },
        "cases": [{
            "trace": "trace.log",
            "reparse_output_sha256": "3" * 64,
            "mandatory_gate_pass": decision == "GO",
            "dimensions": [{
                "dimension": 3,
                "energy": -1.25,
                "residual": 1.0e-14,
                "retained_rank": 3,
            }],
        }],
        "decision": decision,
    }


def require_pass(primary, secondary, label):
    report = comparator.compare_documents(primary, secondary)
    if report["comparison"] != "PASS" or report["mismatch_count"] != 0:
        raise AssertionError(f"comparison should pass ({label}): {report}")
    return report


def require_fail(primary, secondary, label):
    report = comparator.compare_documents(primary, secondary)
    if report["comparison"] != "FAIL" or report["mismatch_count"] == 0:
        raise AssertionError(f"comparison should fail ({label}): {report}")
    json.dumps(report, allow_nan=False)


def main():
    try:
        comparator.strict_json_text(
            '{"decision":"GO","decision":"STOP"}', "duplicate fixture")
    except ValueError:
        pass
    else:
        raise AssertionError("duplicate JSON key was accepted")
    try:
        comparator.strict_json_text('{"value":NaN}', "nonfinite fixture")
    except ValueError:
        pass
    else:
        raise AssertionError("nonfinite JSON token was accepted")

    primary = fixture()
    secondary = copy.deepcopy(primary)
    secondary["metadata"]["generated_at"] = "2026-08-16 13:00 JST"
    secondary["metadata"]["reparse_driver_sha256"] = "4" * 64
    secondary["cases"][0]["reparse_output_sha256"] = "5" * 64
    secondary["cases"][0]["dimensions"][0]["residual"] = math.nextafter(
        primary["cases"][0]["dimensions"][0]["residual"], math.inf)
    report = require_pass(primary, secondary, "allowed provenance and roundoff")
    if report["ignored_value_count"] != 3 or (
        report["float_comparison_count"] != 2
    ):
        raise AssertionError("comparison accounting mismatch")

    primary_stop = fixture("STOP")
    secondary_stop = copy.deepcopy(primary_stop)
    secondary_stop["metadata"]["generated_at"] = "later"
    require_pass(primary_stop, secondary_stop, "matching valid STOP")

    mutated = copy.deepcopy(secondary)
    mutated["decision"] = "STOP"
    require_fail(primary, mutated, "decision mutation")

    mutated = copy.deepcopy(secondary)
    mutated["cases"][0]["dimensions"][0]["energy"] += 1.0e-6
    require_fail(primary, mutated, "numeric mutation")

    mutated = copy.deepcopy(secondary)
    del mutated["cases"][0]["mandatory_gate_pass"]
    require_fail(primary, mutated, "missing gate key")

    mutated = copy.deepcopy(secondary)
    mutated["cases"][0]["dimensions"][0]["residual"] = float("nan")
    require_fail(primary, mutated, "nonfinite float")

    mutated = copy.deepcopy(secondary)
    mutated["metadata"]["reparse_driver_sha256"] = "not-a-sha256"
    require_fail(primary, mutated, "invalid ignored SHA")

    mutated = copy.deepcopy(secondary)
    mutated["metadata"]["generated_at"] = float("nan")
    require_fail(primary, mutated, "invalid ignored timestamp")

    invalid_decision = fixture("UNKNOWN")
    require_fail(invalid_decision, copy.deepcopy(invalid_decision),
                 "invalid matching decision")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
