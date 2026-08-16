#!/usr/bin/env python3

import argparse
import hashlib
import json
import math
import pathlib
import re
import sys


IGNORED_VALUE_FIELDS = frozenset((
    "generated_at",
    "reparse_driver_sha256",
    "reparse_output_sha256",
))
DEFAULT_EPSILON_MULTIPLIER = 32768.0
MAX_REPORTED_MISMATCHES = 64
HEX64 = re.compile(r"[0-9a-f]{64}")


def sha256_path(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def strict_json_text(value, label):
    def reject_constant(token):
        raise ValueError(f"non-finite JSON token in {label}: {token}")

    def reject_duplicate_keys(pairs):
        result = {}
        for key, item in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON key in {label}: {key}")
            result[key] = item
        return result

    return json.loads(value, parse_constant=reject_constant,
                      object_pairs_hook=reject_duplicate_keys)


def path_text(path):
    if not path:
        return "/"
    return "/" + "/".join(str(value).replace("~", "~0").replace("/", "~1")
                            for value in path)


def compare_documents(primary, secondary,
                      epsilon_multiplier=DEFAULT_EPSILON_MULTIPLIER):
    if not math.isfinite(epsilon_multiplier) or epsilon_multiplier < 0.0:
        raise ValueError("epsilon multiplier must be finite and nonnegative")
    mismatches = []
    mismatch_count = 0
    ignored_value_count = 0
    float_comparison_count = 0
    maximum_scaled_error = 0.0
    maximum_scaled_error_path = None

    def diagnostic_value(value):
        if isinstance(value, float) and not math.isfinite(value):
            return repr(value)
        return value

    def mismatch(path, reason, left=None, right=None):
        nonlocal mismatch_count
        mismatch_count += 1
        if len(mismatches) < MAX_REPORTED_MISMATCHES:
            mismatches.append({
                "path": path_text(path),
                "reason": reason,
                "primary": diagnostic_value(left),
                "secondary": diagnostic_value(right),
            })

    def visit(left, right, path):
        nonlocal ignored_value_count
        nonlocal float_comparison_count
        nonlocal maximum_scaled_error
        nonlocal maximum_scaled_error_path

        if isinstance(left, dict) and isinstance(right, dict):
            left_keys = set(left)
            right_keys = set(right)
            if left_keys != right_keys:
                mismatch(path, "object key set differs",
                         sorted(left_keys), sorted(right_keys))
            for key in sorted(left_keys & right_keys):
                if key in IGNORED_VALUE_FIELDS:
                    ignored_value_count += 1
                    values = (left[key], right[key])
                    if key == "generated_at":
                        valid = all(isinstance(value, str) and value
                                    for value in values)
                    else:
                        valid = all(isinstance(value, str) and
                                    HEX64.fullmatch(value) is not None
                                    for value in values)
                    if not valid:
                        mismatch(path + (key,),
                                 "invalid ignored provenance value",
                                 left[key], right[key])
                else:
                    visit(left[key], right[key], path + (key,))
            return
        if isinstance(left, list) and isinstance(right, list):
            if len(left) != len(right):
                mismatch(path, "array length differs", len(left), len(right))
            for index, (left_item, right_item) in enumerate(zip(left, right)):
                visit(left_item, right_item, path + (index,))
            return
        if type(left) is not type(right):
            mismatch(path, "JSON value type differs",
                     type(left).__name__, type(right).__name__)
            return
        if isinstance(left, float):
            float_comparison_count += 1
            if not math.isfinite(left) or not math.isfinite(right):
                mismatch(path, "nonfinite float", left, right)
                return
            scale = max(1.0, abs(left), abs(right))
            difference = abs(left - right)
            scaled_error = difference / (sys.float_info.epsilon * scale)
            if scaled_error > maximum_scaled_error:
                maximum_scaled_error = scaled_error
                maximum_scaled_error_path = path_text(path)
            if difference > epsilon_multiplier * sys.float_info.epsilon * scale:
                mismatch(path, "float tolerance exceeded", left, right)
            return
        if isinstance(left, (str, int, bool)) or left is None:
            if left != right:
                mismatch(path, "exact value differs", left, right)
            return
        mismatch(path, "unsupported JSON value type",
                 type(left).__name__, type(right).__name__)

    visit(primary, secondary, ())
    for label, document in (("primary", primary), ("secondary", secondary)):
        if not isinstance(document, dict) or document.get("decision") not in (
            "GO", "STOP"
        ):
            mismatch((label, "decision"), "invalid top-level decision",
                     document.get("decision")
                     if isinstance(document, dict) else None, "GO or STOP")
    return {
        "schema": 1,
        "epsilon_multiplier": epsilon_multiplier,
        "ignored_value_fields": sorted(IGNORED_VALUE_FIELDS),
        "ignored_value_count": ignored_value_count,
        "float_comparison_count": float_comparison_count,
        "maximum_scaled_error_in_eps": maximum_scaled_error,
        "maximum_scaled_error_path": maximum_scaled_error_path,
        "mismatch_count": mismatch_count,
        "reported_mismatch_count": len(mismatches),
        "mismatches": mismatches,
        "primary_decision": primary.get("decision")
            if isinstance(primary, dict) else None,
        "secondary_decision": secondary.get("decision")
            if isinstance(secondary, dict) else None,
        "comparison": "PASS" if mismatch_count == 0 else "FAIL",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary", required=True, type=pathlib.Path)
    parser.add_argument("--secondary", required=True, type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--epsilon-multiplier", type=float,
                        default=DEFAULT_EPSILON_MULTIPLIER)
    args = parser.parse_args()

    primary = strict_json_text(
        args.primary.read_text(encoding="utf-8"), "primary evidence")
    secondary = strict_json_text(
        args.secondary.read_text(encoding="utf-8"), "secondary evidence")
    report = compare_documents(primary, secondary, args.epsilon_multiplier)
    report.update({
        "primary": args.primary.name,
        "primary_sha256": sha256_path(args.primary),
        "secondary": args.secondary.name,
        "secondary_sha256": sha256_path(args.secondary),
    })
    serialized = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.output is not None:
        args.output.write_text(serialized, encoding="utf-8")
    print(json.dumps(report, sort_keys=True, allow_nan=False))
    return 0 if report["comparison"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
