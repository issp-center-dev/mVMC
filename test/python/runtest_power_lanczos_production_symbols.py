from __future__ import print_function

import argparse
import subprocess
import sys


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", required=True)
    parser.add_argument("--verifier", required=True)
    parser.add_argument("--binary", required=True)
    parser.add_argument("--policy", required=True)
    args = parser.parse_args()

    command = [
        args.python,
        args.verifier,
        "--binary", args.binary,
        "--policy", args.policy,
        "--expected-policy-id",
        "power-lanczos-symbol-policy-oracle-probe-v1",
        "--expected-phase", "TEST-PROBE",
    ]
    result = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True, check=False)
    require(result.returncode == 1,
            "oracle probe was not rejected:\n{}".format(result.stdout))
    require("forbidden production symbols:" in result.stdout,
            "verifier failed for the wrong reason:\n{}".format(
                result.stdout))
    require("oracle_power_lanczos_forbidden_probe" in result.stdout,
            "oracle symbol was not reported:\n{}".format(result.stdout))
    print("power-Lanczos production symbol oracle probe: PASS")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AssertionError as error:
        print("power-Lanczos production symbol oracle probe: FAIL")
        print(error)
        sys.exit(1)
