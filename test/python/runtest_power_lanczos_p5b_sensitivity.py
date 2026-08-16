#!/usr/bin/env python3
"""Unit tests for deterministic P5-B sensitivity evidence generation."""

import json
import hashlib
import os
import subprocess
import sys
import tempfile


HERE = os.path.dirname(os.path.abspath(__file__))
GENERATOR = os.path.join(HERE, "generate_power_lanczos_p5b_sensitivity.py")
POLICY = os.path.join(HERE, "data",
                      "power_lanczos_p5b_sensitivity_policy.json")


def run(policy, output):
    return subprocess.run(
        [sys.executable, GENERATOR, "--policy", policy,
         "--output", output],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, check=False)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    with tempfile.TemporaryDirectory() as temporary:
        output = os.path.join(temporary, "evidence.json")
        result = run(POLICY, output)
        require(result.returncode == 0,
                "valid policy failed:\n{}".format(result.stdout))
        with open(output, "r", encoding="utf-8") as handle:
            evidence = json.load(handle)
        require(evidence["decision"] == "PASS", "decision is not PASS")
        with open(POLICY, "rb") as handle:
            expected_policy_sha = hashlib.sha256(handle.read()).hexdigest()
        require(evidence["policy_sha256"] == expected_policy_sha,
                "policy SHA mismatch")
        require(evidence["case_count"] == 98, "unexpected case count")
        require(len(evidence["envelopes"]) == 7,
                "unexpected envelope count")
        require(evidence["largest_passing_perturbation_scale"] == 1.0e-5,
                "largest passing scale mismatch")
        require(evidence["p5c_candidate_input_envelope"] == 5.0e-6,
                "P5-C candidate envelope mismatch")
        require([row["scale"] for row in evidence["envelopes"]] ==
                [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2],
                "envelope scale order mismatch")
        require([row["passes_all_targets"]
                 for row in evidence["envelopes"]] ==
                [True, True, True, True, False, False, False],
                "envelope pass boundary mismatch")
        require(evidence["identity_maximum"] <= 1.0e-12 and
                evidence["global_phase_invariance_maximum"] <= 1.0e-12,
                "identity/phase gate failed")

        with open(POLICY, "r", encoding="utf-8") as handle:
            policy = json.load(handle)
        invalid_cases = []
        bad = json.loads(json.dumps(policy))
        bad["perturbation_scales"] = list(reversed(
            bad["perturbation_scales"]))
        invalid_cases.append(bad)
        bad = json.loads(json.dumps(policy))
        bad["downstream_absolute_targets"]["energy"] = 0.0
        invalid_cases.append(bad)
        bad = json.loads(json.dumps(policy))
        bad["directions"] = bad["directions"][:-1]
        invalid_cases.append(bad)
        for index, invalid in enumerate(invalid_cases):
            path = os.path.join(temporary, "bad-{}.json".format(index))
            with open(path, "w", encoding="utf-8") as handle:
                json.dump(invalid, handle)
            bad_output = os.path.join(
                temporary, "bad-output-{}.json".format(index))
            result = run(path, bad_output)
            require(result.returncode != 0,
                    "invalid policy {} unexpectedly passed".format(index))
            require(not os.path.exists(bad_output),
                    "invalid policy {} produced output".format(index))
    print("power-Lanczos P5-B sensitivity unit: PASS")


if __name__ == "__main__":
    main()
