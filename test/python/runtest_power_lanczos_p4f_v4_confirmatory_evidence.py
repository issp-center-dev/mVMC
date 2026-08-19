#!/usr/bin/env python3

import json
import pathlib
import sys

import generate_power_lanczos_p4f_v4_confirmatory_evidence as evidence


def policy():
    return {
        "confirmatory_dataset": {
            "driver_schema": 4,
            "seed_order": ["0x5034463452305631"],
            "case_order": [{"site_count": 4, "qp_total": 1}],
            "sample_count_per_trace": 32768,
            "cache_requested_bytes": 268435456,
            "rho": 0.01,
            "rng_stream": 0,
        }
    }


def trace(block_diagnostic_pass=False, hard_pass=True, one_zero=False,
          seed="0x5034463452305631", site_count=4, qp_total=1):
    ratio = 1.0 if block_diagnostic_pass else 1.5
    lines = [
        "MARKOV schema=4 fixture=p4s9_long_direct_session_official "
        f"site_count={site_count} qp_total={qp_total} "
        "sample_count=32768 rank_count=1 "
        "cache_requested_bytes=268435456 rho=0.01 "
        f"seed_hex={seed} rng_stream=0 "
        "persistent_session=1 global_numerator=0 global_denominator=1 "
        "official_block_count=16 official_block_length=2048 "
        "diagnostic_block_count=32 diagnostic_block_length=1024",
    ]
    lines.extend(f"SCALE order={index}" for index in range(4))
    lines.extend(f"SAMPLE sample={index}" for index in range(32768))
    lines.append("SESSION schema=4")
    lines.append("TRACE valid=1")
    lines.append("SUMMARY valid=1")
    lines.extend(f"TAU index={index}" for index in range(37))
    entry_index = 0
    for family in ("S", "K", "B"):
        for row, column in ((0, 0), (0, 1), (0, 2),
                            (1, 1), (1, 2), (2, 2)):
            one_zero_entry = one_zero and entry_index == 0
            official_se = 0.0 if one_zero_entry else 1.0
            diagnostic_se = 1.0 if one_zero_entry else ratio
            entry_ratio = (1.7976931348623157e+308
                           if one_zero_entry else ratio)
            pathology = not one_zero_entry
            lines.append(
                f"ENTRY family={family} row={row} column={column} "
                f"official_se={official_se} diagnostic_se={diagnostic_se} "
                f"stability_ratio={entry_ratio} "
                f"stability_pass={int(block_diagnostic_pass)} "
                f"pathology_pass={int(pathology)} budget_pass=1")
            entry_index += 1
    lines.append("MANIFEST valid=1")
    effective_hard = hard_pass and not one_zero
    lines.append(
        "DECISION "
        f"p4s_decision={'GO' if effective_hard else 'STOP'} "
        f"support_pass={int(hard_pass)} tau_pass=1 budget_pass=1 "
        "denominator_pass=1 "
        f"block_pathology_pass={int(not one_zero)} "
        f"block_stability_pass={int(block_diagnostic_pass)} "
        f"maximum_block_stability_ratio="
        f"{1.7976931348623157e+308 if one_zero else ratio}")
    return ("\n".join(lines) + "\n").encode("ascii")


def expect_failure(function, label):
    try:
        function()
    except (ValueError, AssertionError):
        return
    raise AssertionError(f"negative v4 generator fixture passed: {label}")


def main():
    raw = evidence.parse_raw_trace(trace(block_diagnostic_pass=False))
    identity = evidence.validate_trace_identity(raw, policy())
    if not raw["p4s_hard_pass"] or \
            raw["block_stability_diagnostic_pass"] or \
            not identity["passed"] or not identity["checks"]["p4s"] or \
            identity["block_stability_diagnostic"]["passed"]:
        raise AssertionError(
            "finite block-ratio diagnostic flag must not STOP v4 identity")

    hard = evidence.parse_raw_trace(
        trace(block_diagnostic_pass=True, hard_pass=False))
    hard_identity = evidence.validate_trace_identity(hard, policy())
    if hard["p4s_hard_pass"] or hard_identity["passed"] or \
            hard_identity["checks"]["p4s"]:
        raise AssertionError("P4-S hard field failure must STOP v4 identity")

    pathology = evidence.parse_raw_trace(
        trace(block_diagnostic_pass=False, one_zero=True))
    if pathology["p4s_hard_pass"] or \
            evidence.validate_trace_identity(pathology, policy())["passed"]:
        raise AssertionError("one-zero pathology must STOP v4 identity")

    inconsistent = trace(block_diagnostic_pass=False).replace(
        b"block_stability_pass=0", b"block_stability_pass=1")
    expect_failure(lambda: evidence.parse_raw_trace(inconsistent),
                   "diagnostic aggregate mismatch")
    invalid_decision = trace(block_diagnostic_pass=False).replace(
        b"p4s_decision=GO", b"p4s_decision=INVALID")
    expect_failure(lambda: evidence.parse_raw_trace(invalid_decision),
                   "invalid P4-S decision token")
    duplicate_policy = '{"a": 1, "a": 2}'
    expect_failure(
        lambda: evidence.strict_json_text(duplicate_policy, "fixture"),
        "duplicate policy key")

    if len(sys.argv) > 2:
        raise AssertionError("usage: runtest... [V4_POLICY]")
    if len(sys.argv) == 2:
        current_policy = pathlib.Path(sys.argv[1])
        if evidence.base.sha256_path(current_policy) != \
                evidence.CONFIRMATORY_POLICY_SHA256:
            raise AssertionError("v4 generator/policy SHA binding mismatch")
        value = json.loads(current_policy.read_text(encoding="utf-8"))
        if value["block_stability"][
                "finite_nonzero_ratio_above_diagnostic_maximum"] != \
                "report_only_not_a_hard_gate":
            raise AssertionError("v4 block-ratio policy role mismatch")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
