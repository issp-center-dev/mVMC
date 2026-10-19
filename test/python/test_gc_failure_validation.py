"""Unit checks for the shared invalid-run exit and diagnostic validator."""

from runtest_invalid_expert import (validate_process_result,
                                    validate_rebuild_fault_reason)


def expect_accept(returncode, output, expected, allow_success=False):
    validate_process_result(returncode, output, expected, allow_success)


def expect_reject(returncode, output, expected, allow_success=False):
    try:
        validate_process_result(returncode, output, expected, allow_success)
    except AssertionError:
        return
    raise AssertionError("validator accepted returncode={} output={!r}".format(
        returncode, output))


def main():
    expected = ["Error: GC measurement", "sample:0"]
    expect_accept(1, "Error: GC measurement rank:0 sample:0 failed.",
                  expected)
    expect_accept(1, "MPI_ABORT was invoked\nError: GC measurement "
                     "rank:0 sample:0 failed.", expected)
    expect_accept(0, "Finish calculation.", ["Finish calculation."], True)
    expect_reject(0, "Error: GC measurement rank:0 sample:0 failed.",
                  expected)
    expect_reject(-11, "Error: GC measurement rank:0 sample:0 failed.",
                  expected)
    expect_reject(139, "Error: GC measurement rank:0 sample:0 failed.",
                  expected)
    for report in ("AddressSanitizer: heap-buffer-overflow",
                   "UndefinedBehaviorSanitizer: undefined-behavior",
                   "runtime error: signed integer overflow",
                   "Segmentation fault", "Bus error", "signal 11"):
        expect_reject(1, "{}\nError: GC measurement rank:0 sample:0 "
                         "failed.".format(report), expected)
    expect_reject(1, "Error: unrelated failure", expected)
    validate_rebuild_fault_reason(
        "GC fault injected: kind=rebuild_nonfinite world_rank=1 sample=0 "
        "ncur=2\nError: GC measurement rank:0 sample:0 failed the matrix "
        "rebuild.")
    validate_rebuild_fault_reason(
        "GC fault injected: kind=rebuild_nonfinite world_rank=1 sample=0 "
        "ncur=0\nError: GC measurement rank:0 sample:0 has a zero or "
        "nonfinite amplitude.")
    for bad_reason in (
            "GC fault injected: kind=rebuild_nonfinite world_rank=1 "
            "sample=0 ncur=2\nError: GC measurement rank:0 sample:0 has a "
            "zero or nonfinite amplitude.",
            "GC fault injected: kind=rebuild_nonfinite world_rank=1 "
            "sample=0 ncur=0\nError: GC measurement rank:0 sample:0 failed "
            "the matrix rebuild."):
        try:
            validate_rebuild_fault_reason(bad_reason)
        except AssertionError:
            pass
        else:
            raise AssertionError("validator accepted mismatched rebuild reason")
    print("GC failure validator: PASS")


if __name__ == "__main__":
    main()
