#!/usr/bin/env python3

import math

import analyze_power_lanczos_p4f_block_stability as analysis


def independent_jackknife(blocks):
    denominator = sum(item[1] for item in blocks)
    numerator = sum((item[2] for item in blocks), 0j)
    leave = [(numerator - item[2]) / (denominator - item[1])
             for item in blocks]
    mean = sum(leave) / len(leave)
    factor = (len(leave) - 1.0) / len(leave)
    variance = factor * sum(abs(item - mean) ** 2 for item in leave)
    return math.sqrt(variance)


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError(f"negative block diagnostic fixture passed: {label}")


def main():
    blocks = [
        (8, 7.0 + index, complex((index + 1) ** 2, (-1) ** index))
        for index in range(16)
    ]
    observed = analysis.jackknife(blocks)
    expected = independent_jackknife(blocks)
    if not math.isclose(observed["complex_se"], expected,
                        rel_tol=2.0e-15, abs_tol=1.0e-15):
        raise AssertionError("jackknife oracle mismatch")
    grouped = analysis.aggregate_blocks(blocks, 2)
    if len(grouped) != 8 or any(item[0] != 16 for item in grouped):
        raise AssertionError("contiguous aggregation mismatch")
    if not math.isclose(
            analysis.jackknife(grouped)["complex_se"],
            independent_jackknife(grouped),
            rel_tol=2.0e-15, abs_tol=1.0e-15):
        raise AssertionError("aggregated jackknife oracle mismatch")
    if analysis.symmetric_ratio(0.0, 0.0) != 1.0 or \
            analysis.symmetric_ratio(0.0, 1.0) is not None or \
            analysis.symmetric_ratio(2.0, 3.0) != 1.5:
        raise AssertionError("symmetric ratio boundary mismatch")
    expect_failure(lambda: analysis.jackknife(blocks[:1]),
                   "single jackknife block")
    expect_failure(lambda: analysis.aggregate_blocks(blocks, 3),
                   "nondivisible block aggregation")
    expect_failure(lambda: analysis.fields("ENTRY a=1 a=2"),
                   "duplicate field")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
