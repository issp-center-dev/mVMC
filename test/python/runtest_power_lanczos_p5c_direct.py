from __future__ import print_function

import argparse
import copy
import json
import os
import tempfile

import generate_power_lanczos_p5c_direct as p5c
import generate_power_lanczos_p5b_sensitivity as p5b


def expect_failure(callable_object, label):
    try:
        callable_object()
    except (AssertionError, ValueError, KeyError, TypeError):
        return
    raise AssertionError("negative case accepted: {}".format(label))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("driver")
    parser.add_argument("policy")
    parser.add_argument("p5b_policy")
    arguments = parser.parse_args()
    policy = p5c.validate_policy(p5c.read_json(arguments.policy))

    bad = copy.deepcopy(policy)
    bad["policy_id"] += "-changed"
    expect_failure(lambda: p5c.validate_policy(bad), "policy identity")
    bad = copy.deepcopy(policy)
    bad["schedule"]["sample_count_per_seed"] = 4096
    expect_failure(lambda: p5c.validate_policy(bad), "schedule mutation")
    bad = copy.deepcopy(policy)
    bad["observable_inventory"]["unsupported_fail_fast"] = []
    expect_failure(lambda: p5c.validate_policy(bad), "observable fallback")
    bad = copy.deepcopy(policy)
    bad["p5b_binding"]["candidate_input_envelope"] = 1e-5
    expect_failure(lambda: p5c.validate_policy(bad), "P5-B envelope")

    with tempfile.TemporaryDirectory(prefix="p5c-direct-unit-") as directory:
        output = os.path.join(directory, "evidence.json")
        p5b_evidence = os.path.join(directory, "p5b-evidence.json")
        logs = os.path.join(directory, "logs")
        p5b.atomic_write_json(
            p5b_evidence,
            p5b.generate(p5b.load_policy(arguments.p5b_policy),
                         arguments.p5b_policy))
        p5c.validate_binding(policy, arguments.p5b_policy, p5b_evidence)
        with open(p5b_evidence, "a", encoding="utf-8") as handle:
            handle.write("\n")
        if p5c.sha256_file(p5b_evidence) == policy["p5b_binding"][
                "reference_evidence_sha256"]:
            raise AssertionError("portable binding fixture did not change SHA")
        p5c.validate_binding(policy, arguments.p5b_policy, p5b_evidence)
        with open(p5b_evidence, "r", encoding="utf-8") as handle:
            p5b_payload = json.load(handle)
        bad_evidence = copy.deepcopy(p5b_payload)
        target_scale = policy["p5b_binding"][
            "largest_passing_perturbation_scale"]
        target_row = next(
            row for row in bad_evidence["envelopes"]
            if row["scale"] == target_scale)
        target_row["maximum_downstream_absolute_error"]["energy"] += 1e-8
        p5b.atomic_write_json(p5b_evidence, bad_evidence)
        expect_failure(
            lambda: p5c.validate_binding(
                policy, arguments.p5b_policy, p5b_evidence),
            "P5-B numeric envelope binding")
        p5b.atomic_write_json(p5b_evidence, p5b_payload)
        p5c.validate_binding(policy, arguments.p5b_policy, p5b_evidence)
        payload = p5c.generate(
            arguments.driver, arguments.policy, arguments.p5b_policy,
            p5b_evidence, logs, allow_smoke=True,
            sample_count=512, burn_in=16, block_length=128,
            seed_count=2, restart_split=123)
        p5c.atomic_write_json(output, payload)
        with open(output, "r", encoding="utf-8") as handle:
            loaded = json.load(handle)
        if loaded["decision"] != "SMOKE_PASS" or \
                loaded["decision_fields"]["artifact"] != "PASS" or \
                loaded["decision_fields"][
                    "deterministic_correctness"] != "PASS" or \
                len(loaded["log_inventory"]) != 5:
            raise AssertionError("P5-C smoke evidence mismatch")
        if loaded["exact"]["two_state_variance_gap"][
                "full_support_variance"] != 1.0 or \
                loaded["exact"]["two_state_variance_gap"][
                    "direct_variance_diagnostic"] != 0.0:
            raise AssertionError("variance support contract mismatch")
        negative_statistics = []
        bad_statistics = copy.deepcopy(loaded["statistics"])
        bad_statistics["metric_results"][0]["prefixes"][-1][
            "absolute_error"] = (
                bad_statistics["metric_results"][0][
                    "simultaneous_half_width"] + 0.01)
        negative_statistics.append(("coverage", bad_statistics))
        bad_statistics = copy.deepcopy(loaded["statistics"])
        bad_statistics["traces"]["complex_ring4_near_zero"][0]["tau"][
            "energy_real"] = 16.5
        negative_statistics.append(("tau", bad_statistics))
        bad_statistics = copy.deepcopy(loaded["statistics"])
        bad_statistics["traces"]["two_state_support_bridge"][1][
            "trace_hash"] = bad_statistics["traces"][
                "two_state_support_bridge"][0]["trace_hash"]
        negative_statistics.append(("trace uniqueness", bad_statistics))
        bad_statistics = copy.deepcopy(loaded["statistics"])
        bad_statistics["traces"]["two_state_support_bridge"][0][
            "acceptance_rate"] = 1.0
        negative_statistics.append(("acceptance", bad_statistics))
        bad_statistics = copy.deepcopy(loaded["statistics"])
        bad_statistics["restart_pass"] = False
        negative_statistics.append(("restart", bad_statistics))
        for label, bad_statistics in negative_statistics:
            if p5c.statistical_gate_passes(policy, bad_statistics):
                raise AssertionError(
                    "statistical negative accepted: {}".format(label))
    print("power-Lanczos P5-C direct unit: PASS")


if __name__ == "__main__":
    main()
