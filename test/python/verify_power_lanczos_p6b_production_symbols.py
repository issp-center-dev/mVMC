from __future__ import print_function

import argparse
import json
import os
import re
import subprocess
import sys


DEFAULT_POLICY_ID = "power-lanczos-zero-support-p6b-production-symbols-v2"
DEFAULT_PHASE = "P6-B"


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def read_json(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def normalize_symbol(symbol):
    # Mach-O prepends one underscore to externally visible C symbols.  Policy
    # regexes describe source-level C identifiers and must behave identically
    # on Mach-O and ELF for required and forbidden symbols alike.
    if symbol.startswith("_"):
        return symbol[1:]
    return symbol


def read_global_symbols(binary):
    require(os.path.exists(binary), "binary does not exist: {}".format(binary))
    output = subprocess.check_output(
        ["nm", "-g", binary], stderr=subprocess.STDOUT)
    text = output.decode("utf-8", "replace")
    symbols = set()
    for line in text.splitlines():
        fields = line.split()
        if len(fields) < 2:
            continue
        symbol = normalize_symbol(fields[-1])
        symbols.add(symbol)
    return symbols


def read_binary(binary):
    require(os.path.exists(binary), "binary does not exist: {}".format(binary))
    with open(binary, "rb") as handle:
        return handle.read()


def validate_policy(policy, expected_policy_id, expected_phase):
    require(policy.get("schema_version") == 1, "schema mismatch")
    require(policy.get("policy_id") == expected_policy_id,
            "policy identity mismatch")
    require(policy.get("phase") == expected_phase, "phase mismatch")
    require(policy.get("testing_only") is False,
            "production policy must not be testing-only")
    require(policy.get("production_authorized") is True,
            "production policy authorization mismatch")

    required = policy.get("required_symbols")
    required_strings = policy.get("required_binary_strings", [])
    forbidden = policy.get("forbidden_symbol_regex")
    forbidden_strings = policy.get("forbidden_binary_strings", [])
    require(isinstance(required, list) and required,
            "missing required symbol list")
    require(isinstance(required_strings, list),
            "invalid required binary string list")
    require(isinstance(forbidden, list), "missing forbidden symbol list")
    require(isinstance(forbidden_strings, list),
            "invalid forbidden binary string list")
    for symbol in required:
        require(isinstance(symbol, str) and symbol.startswith("mvmc_"),
                "invalid required symbol {}".format(symbol))
    for pattern in forbidden:
        require(isinstance(pattern, str), "invalid forbidden regex")
        re.compile(pattern)
    for text in required_strings:
        require(isinstance(text, str) and text,
                "invalid required binary string {}".format(text))
        try:
            text.encode("ascii")
        except UnicodeEncodeError:
            raise AssertionError(
                "required binary string is not ASCII: {}".format(text))
    for text in forbidden_strings:
        require(isinstance(text, str) and text.startswith("MVMC_"),
                "invalid forbidden binary string {}".format(text))
        try:
            text.encode("ascii")
        except UnicodeEncodeError:
            raise AssertionError(
                "forbidden binary string is not ASCII: {}".format(text))
    return required, required_strings, forbidden, forbidden_strings


def verify(binary, policy_path, expected_policy_id, expected_phase):
    policy = read_json(policy_path)
    required, required_strings, forbidden_patterns, forbidden_strings = \
        validate_policy(policy, expected_policy_id, expected_phase)
    all_symbols = read_global_symbols(binary)
    binary_payload = read_binary(binary)

    missing = [symbol for symbol in required if symbol not in all_symbols]
    missing_binary_strings = [
        text for text in required_strings
        if text.encode("ascii") not in binary_payload
    ]
    forbidden = []
    for pattern in forbidden_patterns:
        regex = re.compile(pattern)
        forbidden.extend(
            symbol for symbol in all_symbols if regex.search(symbol))
    forbidden = sorted(set(forbidden))
    forbidden_binary_strings = [
        text for text in forbidden_strings
        if text.encode("ascii") in binary_payload
    ]

    if missing or missing_binary_strings or forbidden or forbidden_binary_strings:
        if missing:
            print("missing required production symbols:")
            for symbol in missing:
                print("  {}".format(symbol))
        if missing_binary_strings:
            print("missing required production binary strings:")
            for text in missing_binary_strings:
                print("  {}".format(text))
        if forbidden:
            print("forbidden production symbols:")
            for symbol in forbidden:
                print("  {}".format(symbol))
        if forbidden_binary_strings:
            print("forbidden production binary strings:")
            for text in forbidden_binary_strings:
                print("  {}".format(text))
        return 1

    print("power-Lanczos {} production binary policy: PASS".format(
        expected_phase))
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--expected-policy-id", default=DEFAULT_POLICY_ID)
    parser.add_argument("--expected-phase", default=DEFAULT_PHASE)
    args = parser.parse_args()
    return verify(os.path.abspath(args.binary), os.path.abspath(args.policy),
                  args.expected_policy_id, args.expected_phase)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AssertionError as error:
        print("power-Lanczos production binary policy: FAIL")
        print(error)
        sys.exit(1)
