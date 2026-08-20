from __future__ import print_function

import argparse
import json
import os
import re
import subprocess
import sys


EXPECTED_POLICY_ID = "power-lanczos-zero-support-p6b-production-symbols-v1"


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def read_json(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def normalize_symbol(symbol):
    if symbol.startswith("_mvmc_"):
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
        if symbol.startswith("mvmc_"):
            symbols.add(symbol)
    return symbols


def validate_policy(policy):
    require(policy.get("schema_version") == 1, "schema mismatch")
    require(policy.get("policy_id") == EXPECTED_POLICY_ID,
            "policy identity mismatch")
    require(policy.get("phase") == "P6-B", "phase mismatch")
    require(policy.get("testing_only") is False,
            "P6-B production policy must not be testing-only")
    require(policy.get("production_authorized") is True,
            "P6-B production policy authorization mismatch")

    required = policy.get("required_symbols")
    forbidden = policy.get("forbidden_symbol_regex")
    require(isinstance(required, list) and required,
            "missing required symbol list")
    require(isinstance(forbidden, list), "missing forbidden symbol list")
    for symbol in required:
        require(isinstance(symbol, str) and symbol.startswith("mvmc_"),
                "invalid required symbol {}".format(symbol))
    for pattern in forbidden:
        require(isinstance(pattern, str), "invalid forbidden regex")
        re.compile(pattern)
    return required, forbidden


def verify(binary, policy_path):
    policy = read_json(policy_path)
    required, forbidden_patterns = validate_policy(policy)
    symbols = read_global_symbols(binary)

    missing = [symbol for symbol in required if symbol not in symbols]
    forbidden = []
    for pattern in forbidden_patterns:
        regex = re.compile(pattern)
        forbidden.extend(symbol for symbol in symbols if regex.search(symbol))
    forbidden = sorted(set(forbidden))

    if missing or forbidden:
        if missing:
            print("missing required P6-B production symbols:")
            for symbol in missing:
                print("  {}".format(symbol))
        if forbidden:
            print("forbidden P6-B production symbols:")
            for symbol in forbidden:
                print("  {}".format(symbol))
        return 1

    print("power-Lanczos P6-B production symbol policy: PASS")
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--policy", required=True)
    args = parser.parse_args()
    return verify(os.path.abspath(args.binary), os.path.abspath(args.policy))


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AssertionError as error:
        print("power-Lanczos P6-B production symbol policy: FAIL")
        print(error)
        sys.exit(1)
