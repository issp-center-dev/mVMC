"""Reject crashes and sanitizer reports even when an expected error appears."""
import re

SIGNAL_MARKERS = ("Segmentation fault", "Sanitizer", "runtime error:",
                  "Abort trap", "Bus error")


def reject_signal_reports(output):
    reports = [marker for marker in SIGNAL_MARKERS
               if marker.lower() in output.lower()]
    reports.extend(re.findall(r"\bsignal\s*:?\s*\d+\b", output, re.IGNORECASE))
    if reports:
        raise AssertionError("signal/sanitizer report {}:\n{}".format(
            reports, output[-4000:]))


def require_expected_failure(returncode, output, expected):
    reject_signal_reports(output)
    missing = [message for message in expected if message not in output]
    if returncode <= 0 or missing:
        raise AssertionError("unexpected failure: exit={} missing={}:\n{}".format(
            returncode, missing, output[-4000:]))
