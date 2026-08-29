from __future__ import print_function

import argparse
import hashlib
import json
import math
import os
import statistics
import subprocess
import tempfile
import time


EXPECTED_POLICY_ID = "power-lanczos-zero-support-p5c-testing-v2"
EXPECTED_FIXTURES = [
    "two_state_support_bridge",
    "two_state_variance_gap",
    "complex_ring4_near_zero",
]
STOCHASTIC_FIXTURES = [
    "two_state_support_bridge",
    "complex_ring4_near_zero",
]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def read_json(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_json(path, payload):
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=".p5c-direct-", suffix=".json", dir=directory)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def validate_policy(policy):
    require(policy.get("schema_version") == 1, "P5-C schema mismatch")
    require(policy.get("policy_id") == EXPECTED_POLICY_ID,
            "P5-C policy identity mismatch")
    require(policy.get("phase") == "P5-C" and
            policy.get("testing_only") is True and
            policy.get("production_authorized") is False,
            "P5-C authorization mismatch")
    require(policy.get("fixtures") == EXPECTED_FIXTURES,
            "P5-C fixture order mismatch")
    binding = policy.get("p5b_binding", {})
    require(binding.get("candidate_input_envelope") == 5e-6 and
            binding.get("largest_passing_perturbation_scale") == 1e-5 and
            binding.get("candidate_envelope_safety_factor") == 0.5,
            "P5-B envelope binding mismatch")
    require(policy.get("revision_reason") ==
            "portable_numeric_p5b_binding_and_independent_exact_gate",
            "P5-C revision reason mismatch")
    reference_sha = binding.get("reference_evidence_sha256", "")
    require(isinstance(reference_sha, str) and len(reference_sha) == 64 and
            all(character in "0123456789abcdef"
                for character in reference_sha),
            "invalid P5-B reference evidence SHA")
    tolerance = float(binding.get("binding_numeric_tolerance", -1.0))
    require(math.isfinite(tolerance) and 0.0 < tolerance <= 1e-10,
            "invalid P5-B numeric binding tolerance")
    maxima = binding.get("reference_envelope_maxima", {})
    for name in ("energy", "full_support_variance",
                 "quadratic_observable"):
        require(math.isfinite(float(maxima.get(name, -1.0))) and
                float(maxima[name]) >= 0.0,
                "invalid reference envelope maximum {}".format(name))
    proposal = policy.get("proposal", {})
    require(proposal == {
        "kernel": "canonical_uniform_sector",
        "global_numerator": 1,
        "global_denominator": 1,
        "required_target_support_connectivity": True,
    }, "P5-C proposal mismatch")
    schedule = policy.get("schedule", {})
    require(schedule.get("burn_in") == 1024 and
            schedule.get("sample_count_per_seed") == 65536 and
            schedule.get("prefix_sample_counts") == [4096, 16384, 65536] and
            schedule.get("block_length") == 1024 and
            schedule.get("restart_split_after_completed_steps") == 12345,
            "P5-C schedule mismatch")
    seeds = policy.get("rng", {}).get("seed_hex_order", [])
    require(len(seeds) == 8 and len(set(seeds)) == 8,
            "P5-C seed census mismatch")
    for seed in seeds + [policy.get("rng", {}).get("stream_hex")]:
        require(isinstance(seed, str) and seed.startswith("0x") and
                int(seed, 0) > 0, "invalid P5-C seed/stream")
    coverage = policy.get("coverage", {})
    require(len(coverage.get("metric_census", [])) == 5 and
            coverage.get("family_wise_alpha") == 0.01 and
            coverage.get("multiplicity") ==
            "Bonferroni_two_sided_normal" and
            abs(float(coverage.get("normal_multiplier", 0.0)) -
                3.090232306167813) <= 1e-15 and
            float(coverage.get("maximum_tau_int", 0.0)) == 16.0 and
            float(coverage.get(
                "maximum_full_to_first_conservative_se_ratio", 0.0)) ==
            0.65,
            "P5-C coverage mismatch")
    require(policy.get("near_zero_scan_relative_amplitude_cutoffs") ==
            [0.0, 1e-8, 1e-6, 5e-6, 1e-5, 1e-4],
            "P5-C cutoff schedule mismatch")
    inventory = policy.get("observable_inventory", {})
    require("direct_energy_complex" in inventory.get("supported", []) and
            "full_support_variance_quadratic_reference" in
            inventory.get("supported", []) and
            "arbitrary_off_diagonal" in
            inventory.get("unsupported_fail_fast", []),
            "P5-C observable inventory mismatch")
    require(policy.get("required_decision_fields") == [
        "artifact", "deterministic_correctness",
        "statistical_adequacy", "resource"],
        "P5-C decision-role mismatch")
    return policy


def validate_binding(policy, p5b_policy_path, p5b_evidence_path):
    binding = policy["p5b_binding"]
    policy_sha = sha256_file(p5b_policy_path)
    require(policy_sha == binding["policy_sha256"],
            "P5-B policy SHA mismatch")
    evidence = read_json(p5b_evidence_path)
    tolerance = float(binding["binding_numeric_tolerance"])

    def close_numeric(actual, expected):
        actual = float(actual)
        expected = float(expected)
        return (math.isfinite(actual) and math.isfinite(expected) and
                abs(actual - expected) <=
                tolerance * max(1.0, abs(actual), abs(expected)))

    require(evidence.get("decision") == "PASS" and
            evidence.get("policy_sha256") == policy_sha and
            close_numeric(evidence.get("p5c_candidate_input_envelope"),
                          binding["candidate_input_envelope"]) and
            close_numeric(
                evidence.get("largest_passing_perturbation_scale"),
                binding["largest_passing_perturbation_scale"]),
            "P5-B evidence content mismatch")
    envelope_rows = [
        row for row in evidence.get("envelopes", [])
        if close_numeric(
            row.get("scale"),
            binding["largest_passing_perturbation_scale"])
    ]
    require(len(envelope_rows) == 1 and
            envelope_rows[0].get("passes_all_targets") is True,
            "P5-B envelope row mismatch")
    actual_maxima = envelope_rows[0].get(
        "maximum_downstream_absolute_error", {})
    for name, expected in binding["reference_envelope_maxima"].items():
        require(name in actual_maxima and
                close_numeric(actual_maxima[name], expected),
                "P5-B envelope maximum mismatch: {}".format(name))


def parse_scalar(value):
    if value in ("PASS", "FAIL", "RECORDED", "NA",
                 "global_uniform"):
        return value
    try:
        return int(value, 0)
    except ValueError:
        try:
            parsed = float(value)
        except ValueError:
            return value
        require(math.isfinite(parsed), "nonfinite driver scalar")
        return parsed


def parse_fields(line):
    fields = {}
    for item in line.split()[1:]:
        key, value = item.split("=", 1)
        require(key not in fields, "duplicate driver field {}".format(key))
        fields[key] = parse_scalar(value)
    return fields


def parse_driver_log(path):
    parsed = {"anchors": [], "blocks": []}
    single = {
        "RUN": "run", "CONNECTIVITY": "connectivity", "TRACE": "trace",
        "TAU": "tau", "TAIL": "tail", "MANIFEST": "manifest",
        "DECISION": "decision",
    }
    with open(path, "r", encoding="utf-8") as handle:
        text = handle.read()
    for line in text.splitlines():
        kind = line.split(" ", 1)[0]
        if kind == "ANCHOR":
            parsed["anchors"].append(parse_fields(line))
        elif kind == "BLOCK":
            parsed["blocks"].append(parse_fields(line))
        elif kind in single:
            key = single[kind]
            require(key not in parsed, "duplicate {}".format(kind))
            parsed[key] = parse_fields(line)
        else:
            raise AssertionError("unknown driver line {}".format(line))
    for key in ("run", "connectivity", "decision"):
        require(key in parsed, "missing driver {}".format(key))
    parsed["path"] = os.path.abspath(path)
    parsed["sha256"] = sha256_file(path)
    parsed["line_count"] = len(text.splitlines())
    return parsed


def run_driver(driver, fixture, sample_count, burn_in, block_length,
               seed, stream, restart_split, output_path):
    command = [
        driver, fixture, str(sample_count), str(burn_in), str(block_length),
        seed, stream, str(restart_split),
    ]
    environment = os.environ.copy()
    environment.update({
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
    })
    completed = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, env=environment, check=False, timeout=300)
    if completed.returncode != 0:
        raise AssertionError(
            "driver failed {}\nstdout:\n{}\nstderr:\n{}".format(
                command, completed.stdout, completed.stderr))
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write(completed.stdout)
    return parse_driver_log(output_path)


def complex_field(anchor, prefix):
    return complex(float(anchor[prefix + "_real"]),
                   float(anchor[prefix + "_imag"]))


def anchors_equal(left, right):
    keys = [
        "index", "words", "v0_real", "v0_imag", "v1_real", "v1_imag",
        "v2_real", "v2_imag", "psi_real", "psi_imag", "hpsi_real",
        "hpsi_imag", "sampleable", "psi0_state", "psi_state",
        "hpsi_state",
    ]
    return [{key: item[key] for key in keys} for item in left] == [
        {key: item[key] for key in keys} for item in right]


def close_complex(actual, expected, tolerance=5e-13):
    return abs(actual - expected) <= tolerance * max(1.0, abs(expected))


def ring4_apply(vector):
    require(len(vector) == 4, "ring4 vector dimension mismatch")
    return [
        vector[3] + vector[1],
        vector[0] + vector[2],
        vector[1] + vector[3],
        vector[2] + vector[0],
    ]


def validate_anchor_oracle(fixture, anchors):
    require([int(item["index"]) for item in anchors] ==
            list(range(len(anchors))), "anchor index mismatch")
    if fixture in ("two_state_support_bridge", "two_state_variance_gap"):
        v0 = [1.0 + 0.0j, 0.0 + 0.0j]
        v1 = [0.0 + 0.0j, 1.0 + 0.0j]
        v2 = [1.0 + 0.0j, 0.0 + 0.0j]
        alpha = [1.0, 2.0 if fixture ==
                 "two_state_support_bridge" else 0.0]
        words = [1, 2]
    elif fixture == "complex_ring4_near_zero":
        v0 = [1.0 + 0.0j, 0.0 + 0.0j,
              0.05001 - 0.025j, -0.2 + 0.1j]
        v1 = ring4_apply(v0)
        v2 = ring4_apply(v1)
        alpha = [1.0, 0.25]
        words = [1, 2, 4, 8]
    else:
        raise AssertionError("unknown exact fixture {}".format(fixture))
    require(len(anchors) == len(words), "anchor sector count mismatch")
    for index, anchor in enumerate(anchors):
        expected_psi = alpha[0] * v0[index] + alpha[1] * v1[index]
        expected_hpsi = alpha[0] * v1[index] + alpha[1] * v2[index]
        require(int(anchor["words"]) == words[index] and
                close_complex(complex_field(anchor, "v0"), v0[index]) and
                close_complex(complex_field(anchor, "v1"), v1[index]) and
                close_complex(complex_field(anchor, "v2"), v2[index]) and
                close_complex(
                    complex_field(anchor, "psi"), expected_psi) and
                close_complex(
                    complex_field(anchor, "hpsi"), expected_hpsi) and
                int(anchor["sampleable"]) == int(expected_psi != 0.0),
                "independent anchor oracle mismatch: {} state {}".format(
                    fixture, index))


def exact_summary(anchors, policy):
    psi = [complex_field(item, "psi") for item in anchors]
    hpsi = [complex_field(item, "hpsi") for item in anchors]
    basis = [[complex_field(item, "v{}".format(order))
              for order in range(3)] for item in anchors]
    weights = [abs(value) ** 2 for value in psi]
    norm = sum(weights)
    require(norm > 0.0 and math.isfinite(norm), "invalid exact norm")
    numerator = sum(value.conjugate() * shifted
                    for value, shifted in zip(psi, hpsi))
    energy = numerator / norm
    second = sum(abs(value) ** 2 for value in hpsi) / norm
    full_variance = second - energy.real ** 2
    supported = [index for index, item in enumerate(anchors)
                 if int(item["sampleable"]) == 1]
    direct_norm = sum(weights[index] for index in supported)
    direct_energy = sum(
        psi[index].conjugate() * hpsi[index] for index in supported
    ) / direct_norm
    direct_second = sum(abs(hpsi[index]) ** 2
                        for index in supported) / direct_norm
    direct_variance = direct_second - direct_energy.real ** 2

    mixture = policy["mixture_reference"]
    eta = float(mixture["eta"])
    lambdas = [float(item) for item in mixture["lambda"]]
    guide = [eta + sum(lambdas[order] * abs(row[order]) ** 2
                       for order in range(3)) for row in basis]
    require(all(value > 0.0 and math.isfinite(value) for value in guide),
            "mixture guide is not strictly positive")
    mixture_norm = sum(guide[index] * weights[index] / guide[index]
                       for index in range(len(anchors)))
    mixture_energy_numerator = sum(
        guide[index] * (psi[index].conjugate() * hpsi[index]) /
        guide[index] for index in range(len(anchors)))
    mixture_second_numerator = sum(
        guide[index] * abs(hpsi[index]) ** 2 / guide[index]
        for index in range(len(anchors)))
    mixture_energy = mixture_energy_numerator / mixture_norm
    mixture_second = mixture_second_numerator / mixture_norm
    mixture_variance = mixture_second - mixture_energy.real ** 2
    tolerance = float(mixture["scale_aware_tolerance"])
    require(abs(mixture_energy - energy) <=
            tolerance * max(1.0, abs(energy)) and
            abs(mixture_variance - full_variance) <=
            tolerance * max(1.0, abs(full_variance)),
            "exact mixture reweighting mismatch")

    maximum_amplitude = max(abs(value) for value in psi)
    scans = []
    for cutoff in policy["near_zero_scan_relative_amplitude_cutoffs"]:
        retained = [index for index, value in enumerate(psi)
                    if cutoff == 0.0 or
                    abs(value) / maximum_amplitude >= cutoff]
        cut_norm = sum(weights[index] for index in retained)
        cut_energy = sum(psi[index].conjugate() * hpsi[index]
                         for index in retained) / cut_norm
        cut_second = sum(abs(hpsi[index]) ** 2
                         for index in retained) / cut_norm
        scans.append({
            "relative_amplitude_cutoff": cutoff,
            "excluded_state_count": len(anchors) - len(retained),
            "retained_norm": cut_norm,
            "energy_real": cut_energy.real,
            "energy_imag": cut_energy.imag,
            "direct_second_moment": cut_second,
            "energy_absolute_drift": abs(cut_energy - energy),
            "second_moment_absolute_drift": abs(cut_second - second),
        })
    return {
        "sector_count": len(anchors), "target_support_count": len(supported),
        "norm": norm, "energy_real": energy.real,
        "energy_imag": energy.imag,
        "full_second_moment": second,
        "full_support_variance": full_variance,
        "direct_energy_real": direct_energy.real,
        "direct_energy_imag": direct_energy.imag,
        "direct_second_moment": direct_second,
        "direct_variance_diagnostic": direct_variance,
        "mixture_energy_real": mixture_energy.real,
        "mixture_energy_imag": mixture_energy.imag,
        "mixture_full_support_variance": mixture_variance,
        "mixture_check_classification": "ALGEBRAIC_CONSISTENCY_ONLY",
        "near_zero_scan": scans,
        "minimum_relative_nonzero_amplitude": min(
            abs(value) / maximum_amplitude for value in psi if value != 0.0),
        "psi0_zero_final_nonzero_count": sum(
            int(abs(basis[index][0]) == 0.0 and abs(psi[index]) > 0.0)
            for index in range(len(anchors))),
        "psi_zero_hpsi_nonzero_count": sum(
            int(abs(psi[index]) == 0.0 and abs(hpsi[index]) > 0.0)
            for index in range(len(anchors))),
    }


def validate_exact_fixture(fixture, summary):
    tolerance = 5e-13
    if fixture == "two_state_support_bridge":
        require(abs(summary["energy_real"] - 0.8) <= tolerance and
                abs(summary["energy_imag"]) <= tolerance and
                abs(summary["full_support_variance"] - 0.36) <= tolerance and
                abs(summary["direct_variance_diagnostic"] - 0.36) <=
                tolerance and summary["target_support_count"] == 2 and
                summary["psi0_zero_final_nonzero_count"] == 1,
                "support-bridge exact contract mismatch")
    elif fixture == "two_state_variance_gap":
        require(abs(summary["energy_real"]) <= tolerance and
                abs(summary["full_support_variance"] - 1.0) <= tolerance and
                abs(summary["direct_variance_diagnostic"]) <= tolerance and
                summary["target_support_count"] == 1 and
                summary["psi_zero_hpsi_nonzero_count"] == 1,
                "variance-gap exact contract mismatch")
    elif fixture == "complex_ring4_near_zero":
        require(summary["sector_count"] == 4 and
                summary["target_support_count"] == 4 and
                summary["psi0_zero_final_nonzero_count"] == 1 and
                abs(summary["energy_real"] -
                    0.6315106437881186) <= tolerance and
                abs(summary["energy_imag"]) <= tolerance and
                abs(summary["full_support_variance"] -
                    1.6654963030334775) <= tolerance and
                abs(summary["minimum_relative_nonzero_amplitude"] -
                    1.0522672835285087e-5) <= tolerance,
                "ring4 near-zero exact contract mismatch")


def sample_se(values):
    if len(values) < 2:
        return float("inf")
    return statistics.stdev(values) / math.sqrt(len(values))


def statistical_gate_passes(policy, result):
    try:
        coverage = policy["coverage"]
        metrics = result["metric_results"]
        if len(metrics) != len(coverage["metric_census"]):
            return False
        coverage_pass = True
        scaling_pass = True
        for metric in metrics:
            prefixes = metric["prefixes"]
            if len(prefixes) != 3:
                return False
            full = prefixes[-1]
            covered = (float(full["absolute_error"]) <=
                       float(metric["simultaneous_half_width"]))
            scaled = (float(metric["full_to_first_se_ratio"]) <=
                      float(coverage[
                          "maximum_full_to_first_conservative_se_ratio"]))
            if bool(metric["covered"]) != covered or +                    bool(metric["sample_scaling_pass"]) != scaled:
                return False
            coverage_pass = coverage_pass and covered
            scaling_pass = scaling_pass and scaled
        trace_hashes = []
        maximum_tau = 0.0
        acceptance_pass = True
        for fixture in STOCHASTIC_FIXTURES:
            fixture_traces = result["traces"][fixture]
            if len(fixture_traces) < 2:
                return False
            for trace in fixture_traces:
                trace_hashes.append(int(trace["trace_hash"]))
                rate = float(trace["acceptance_rate"])
                acceptance_pass = (
                    acceptance_pass and
                    float(coverage["minimum_acceptance_rate_exclusive"]) <
                    rate <
                    float(coverage["maximum_acceptance_rate_exclusive"]) and
                    int(trace["changed"]) > 0 and
                    int(trace["exact_zero_proposals"]) == 0)
                maximum_tau = max(
                    maximum_tau,
                    *(float(value) for value in trace["tau"].values()))
        unique_trace_pass = len(trace_hashes) == len(set(trace_hashes))
        tau_pass = maximum_tau <= float(coverage["maximum_tau_int"])
        if (bool(result["all_simultaneous_coverage_pass"]) != coverage_pass or
                bool(result["all_sample_scaling_pass"]) != scaling_pass or
                bool(result["acceptance_pass"]) != acceptance_pass or
                bool(result["independent_trace_hashes_pass"]) !=
                unique_trace_pass or
                bool(result["tau_pass"]) != tau_pass or
                abs(float(result["maximum_tau_int"]) - maximum_tau) >
                1e-13 * max(1.0, maximum_tau)):
            return False
        return (coverage_pass and scaling_pass and tau_pass and
                acceptance_pass and bool(result["restart_pass"]) and
                unique_trace_pass)
    except (KeyError, TypeError, ValueError, OverflowError):
        return False


def metric_definition(fixture, metric, exact):
    if metric == "energy_real":
        return exact["energy_real"]
    if metric == "energy_imag":
        return exact["energy_imag"]
    if metric == "second":
        return exact["direct_second_moment"]
    raise AssertionError("unknown metric {} {}".format(fixture, metric))


def analyze_statistics(runs, exact_by_fixture, policy, prefixes,
                       block_length, enforce_decision):
    coverage = policy["coverage"]
    multiplier = float(coverage["normal_multiplier"])
    metric_map = [
        ("two_state_support_bridge", "energy_real"),
        ("two_state_support_bridge", "energy_imag"),
        ("two_state_support_bridge", "second"),
        ("complex_ring4_near_zero", "energy_real"),
        ("complex_ring4_near_zero", "energy_imag"),
    ]
    require(["{}.{}".format(fixture,
                             "direct_second_moment" if metric == "second"
                             else metric)
             for fixture, metric in metric_map] ==
            coverage["metric_census"], "metric census/code mismatch")
    results = []
    all_coverage = True
    all_scaling = True
    for fixture, metric in metric_map:
        exact_value = metric_definition(
            fixture, metric, exact_by_fixture[fixture])
        prefix_rows = []
        for prefix in prefixes:
            require(prefix % block_length == 0,
                    "prefix is not block aligned")
            blocks_per_seed = prefix // block_length
            block_values = []
            for run in runs[fixture]:
                require(len(run["blocks"]) >= blocks_per_seed,
                        "insufficient blocks")
                block_values.extend(float(item[metric])
                                    for item in run["blocks"][:blocks_per_seed])
            estimate = statistics.fmean(block_values)
            se = sample_se(block_values)
            prefix_rows.append({
                "sample_count_per_seed": prefix,
                "pooled_block_count": len(block_values),
                "estimate": estimate,
                "conservative_se": se,
                "absolute_error": abs(estimate - exact_value),
            })
        first_se = prefix_rows[0]["conservative_se"]
        full_se = prefix_rows[-1]["conservative_se"]
        if first_se == 0.0:
            se_ratio = 0.0 if full_se == 0.0 else float("inf")
        else:
            se_ratio = full_se / first_se
        half_width = multiplier * full_se
        covered = prefix_rows[-1]["absolute_error"] <= half_width
        scaling_pass = se_ratio <= float(coverage[
            "maximum_full_to_first_conservative_se_ratio"])
        all_coverage = all_coverage and covered
        all_scaling = all_scaling and scaling_pass
        results.append({
            "fixture": fixture, "metric": metric,
            "exact": exact_value, "prefixes": prefix_rows,
            "simultaneous_half_width": half_width,
            "covered": covered,
            "full_to_first_se_ratio": se_ratio,
            "sample_scaling_pass": scaling_pass,
        })

    maximum_tau = 0.0
    acceptance_pass = True
    restart_pass = True
    unique_trace_pass = True
    tails = {}
    traces = {}
    direct_variance_diagnostics = {}
    for fixture in STOCHASTIC_FIXTURES:
        fixture_runs = runs[fixture]
        trace_hashes = [int(run["trace"]["trace_hash"])
                        for run in fixture_runs]
        unique_trace_pass = (unique_trace_pass and
                             len(set(trace_hashes)) == len(trace_hashes))
        fixture_tails = []
        fixture_traces = []
        for run in fixture_runs:
            maximum_tau = max(maximum_tau,
                              *(float(value) for value in run["tau"].values()))
            rate = float(run["trace"]["acceptance_rate"])
            acceptance_pass = acceptance_pass and 0.0 < rate < 1.0
            restart_pass = restart_pass and int(
                run["trace"]["restart_pass"]) == 1
            fixture_tails.append(run["tail"])
            fixture_traces.append({
                "seed": run["run"]["seed"],
                "trace_hash": run["trace"]["trace_hash"],
                "acceptance_rate": rate,
                "accepted": run["trace"]["accepted"],
                "changed": run["trace"]["changed"],
                "exact_zero_proposals":
                    run["trace"]["exact_zero_proposals"],
                "tau": run["tau"],
                "manifest": run["manifest"],
            })
        tails[fixture] = fixture_tails
        traces[fixture] = fixture_traces
        full_blocks = {
            name: [float(block[name]) for run in fixture_runs
                   for block in run["blocks"]]
            for name in ("energy_real", "energy_imag", "second")
        }
        energy_real_estimate = statistics.fmean(
            full_blocks["energy_real"])
        energy_imag_estimate = statistics.fmean(
            full_blocks["energy_imag"])
        second_estimate = statistics.fmean(full_blocks["second"])
        direct_variance_estimate = (
            second_estimate - energy_real_estimate ** 2 -
            energy_imag_estimate ** 2)
        exact = exact_by_fixture[fixture]
        reference_eligible = fixture == "two_state_support_bridge"
        direct_variance_diagnostics[fixture] = {
            "energy_real_estimate": energy_real_estimate,
            "energy_imag_estimate": energy_imag_estimate,
            "direct_second_moment_estimate": second_estimate,
            "direct_variance_estimate": direct_variance_estimate,
            "exact_direct_variance": exact["direct_variance_diagnostic"],
            "full_support_variance_reference":
                exact["full_support_variance"],
            "absolute_error_vs_full_reference": abs(
                direct_variance_estimate -
                exact["full_support_variance"]),
            "correctness_reference_eligible": reference_eligible,
            "classification": "SUPPORT_SAFE_COVERAGE_METRIC" if
                reference_eligible else "HEAVY_TAIL_DIAGNOSTIC_ONLY",
            "maximum_observed_local_energy_abs": max(
                float(item["maximum"]) for item in fixture_tails),
        }
    tau_pass = maximum_tau <= float(coverage["maximum_tau_int"])
    payload = {
        "metric_results": results,
        "all_simultaneous_coverage_pass": all_coverage,
        "all_sample_scaling_pass": all_scaling,
        "maximum_tau_int": maximum_tau,
        "tau_pass": tau_pass,
        "acceptance_pass": acceptance_pass,
        "restart_pass": restart_pass,
        "independent_trace_hashes_pass": unique_trace_pass,
        "tails": tails,
        "traces": traces,
        "direct_variance_diagnostics": direct_variance_diagnostics,
    }
    statistical_pass = statistical_gate_passes(policy, payload)
    if enforce_decision:
        require(statistical_pass, "P5-C statistical gate failed")
    payload["decision"] = (
        "PASS" if statistical_pass else
        ("SMOKE_RECORDED" if not enforce_decision else "FAIL"))
    return payload


def generate(driver, policy_path, p5b_policy_path, p5b_evidence_path,
             logs_dir, allow_smoke=False, sample_count=None, burn_in=None,
             block_length=None, seed_count=None, restart_split=None):
    start = time.perf_counter()
    policy = validate_policy(read_json(policy_path))
    validate_binding(policy, p5b_policy_path, p5b_evidence_path)
    schedule = policy["schedule"]
    official_sample_count = int(schedule["sample_count_per_seed"])
    official_burn_in = int(schedule["burn_in"])
    official_block_length = int(schedule["block_length"])
    official_restart = int(schedule["restart_split_after_completed_steps"])
    if allow_smoke:
        sample_count = 512 if sample_count is None else int(sample_count)
        burn_in = 16 if burn_in is None else int(burn_in)
        block_length = 128 if block_length is None else int(block_length)
        seed_count = 2 if seed_count is None else int(seed_count)
        restart_split = 123 if restart_split is None else int(restart_split)
        prefixes = [sample_count // 4, sample_count // 2, sample_count]
    else:
        require(sample_count is None and burn_in is None and
                block_length is None and seed_count is None and
                restart_split is None,
                "official P5-C schedule cannot be overridden")
        sample_count = official_sample_count
        burn_in = official_burn_in
        block_length = official_block_length
        seed_count = len(policy["rng"]["seed_hex_order"])
        restart_split = official_restart
        prefixes = list(schedule["prefix_sample_counts"])
    require(sample_count > 0 and burn_in >= 0 and block_length > 0 and
            sample_count % block_length == 0 and seed_count >= 2 and
            all(prefix > 0 and prefix % block_length == 0 and
                prefix <= sample_count for prefix in prefixes) and
            0 < restart_split < burn_in + sample_count,
            "invalid resolved P5-C schedule")

    seeds = policy["rng"]["seed_hex_order"][:seed_count]
    stream = policy["rng"]["stream_hex"]
    os.makedirs(logs_dir, exist_ok=True)
    runs = {fixture: [] for fixture in STOCHASTIC_FIXTURES}
    exact_by_fixture = {}
    log_inventory = []
    for fixture in STOCHASTIC_FIXTURES:
        reference_anchors = None
        for seed_index, seed in enumerate(seeds):
            path = os.path.join(
                logs_dir, "{}-seed{:02d}.log".format(fixture, seed_index))
            parsed = run_driver(
                driver, fixture, sample_count, burn_in, block_length,
                seed, stream, restart_split, path)
            require(parsed["run"]["schema"] == 2 and
                    parsed["run"]["fixture"] == fixture and
                    parsed["run"]["sample_count"] == sample_count and
                    parsed["run"]["burn_in"] == burn_in and
                    parsed["run"]["block_length"] == block_length and
                    parsed["run"]["seed"] == int(seed, 0),
                    "driver run header mismatch")
            require(parsed["connectivity"]["connected_by_policy"] == 1 and
                    parsed["connectivity"]["proposal"] ==
                    "global_uniform" and
                    parsed["connectivity"]["pairwise_enumerated"] == 0 and
                    parsed["connectivity"]["positive_pair_count"] ==
                    parsed["connectivity"]["target_support_count"] ** 2,
                    "target-support connectivity mismatch")
            require(len(parsed["blocks"]) == sample_count // block_length and
                    [int(block["index"]) for block in parsed["blocks"]] ==
                    list(range(sample_count // block_length)) and
                    parsed["decision"] == {
                        "deterministic": "PASS",
                        "stochastic": "RECORDED"},
                    "driver stochastic census mismatch")
            if reference_anchors is None:
                reference_anchors = parsed["anchors"]
                validate_anchor_oracle(fixture, reference_anchors)
                summary = exact_summary(reference_anchors, policy)
                validate_exact_fixture(fixture, summary)
                exact_by_fixture[fixture] = summary
            else:
                require(anchors_equal(reference_anchors, parsed["anchors"]),
                        "anchor drift across seeds")
            runs[fixture].append(parsed)
            log_inventory.append({
                "fixture": fixture, "seed_index": seed_index,
                "path": os.path.relpath(
                    parsed["path"],
                    os.path.dirname(os.path.abspath(logs_dir))),
                "sha256": parsed["sha256"],
                "line_count": parsed["line_count"],
            })

    variance_path = os.path.join(logs_dir, "two_state_variance_gap.log")
    variance_run = run_driver(
        driver, "two_state_variance_gap", 0, 0, 0, seeds[0], stream, 0,
        variance_path)
    require(variance_run["decision"] == {
        "deterministic": "PASS", "stochastic": "NA"},
        "variance-gap deterministic decision mismatch")
    validate_anchor_oracle(
        "two_state_variance_gap", variance_run["anchors"])
    variance_summary = exact_summary(variance_run["anchors"], policy)
    validate_exact_fixture("two_state_variance_gap", variance_summary)
    exact_by_fixture["two_state_variance_gap"] = variance_summary
    log_inventory.append({
        "fixture": "two_state_variance_gap", "seed_index": None,
        "path": os.path.relpath(
            variance_run["path"],
            os.path.dirname(os.path.abspath(logs_dir))),
        "sha256": variance_run["sha256"],
        "line_count": variance_run["line_count"],
    })

    statistics_result = analyze_statistics(
        runs, exact_by_fixture, policy, prefixes, block_length,
        enforce_decision=not allow_smoke)
    elapsed = time.perf_counter() - start
    decision = "PASS" if statistics_result["decision"] == "PASS" else \
        "SMOKE_PASS"
    return {
        "schema": "power_lanczos_p5c_direct_evidence_v2",
        "phase": "P5-C", "testing_only": True,
        "production_authorized": False,
        "policy_sha256": sha256_file(policy_path),
        "p5b_policy_sha256": sha256_file(p5b_policy_path),
        "p5b_evidence_sha256": sha256_file(p5b_evidence_path),
        "resolved_schedule": {
            "official": not allow_smoke, "sample_count_per_seed": sample_count,
            "burn_in": burn_in, "block_length": block_length,
            "prefix_sample_counts": prefixes, "seed_count": seed_count,
            "restart_split_after_completed_steps": restart_split,
        },
        "exact": exact_by_fixture,
        "statistics": statistics_result,
        "log_inventory": log_inventory,
        "decision_fields": {
            "artifact": "PASS",
            "deterministic_correctness": "PASS",
            "statistical_adequacy": statistics_result["decision"],
            "resource": "RECORDED_LOCAL_NOT_PRODUCTION_GATED",
        },
        "elapsed_seconds": elapsed,
        "decision": decision,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--driver", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--p5b-policy", required=True)
    parser.add_argument("--p5b-evidence", required=True)
    parser.add_argument("--logs-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--allow-smoke", action="store_true")
    parser.add_argument("--sample-count", type=int)
    parser.add_argument("--burn-in", type=int)
    parser.add_argument("--block-length", type=int)
    parser.add_argument("--seed-count", type=int)
    parser.add_argument("--restart-split", type=int)
    arguments = parser.parse_args()
    payload = generate(
        arguments.driver, arguments.policy, arguments.p5b_policy,
        arguments.p5b_evidence, arguments.logs_dir,
        allow_smoke=arguments.allow_smoke,
        sample_count=arguments.sample_count, burn_in=arguments.burn_in,
        block_length=arguments.block_length, seed_count=arguments.seed_count,
        restart_split=arguments.restart_split)
    atomic_write_json(arguments.output, payload)
    print("P5-C direct evidence: {} ({} seeds x {} samples)".format(
        payload["decision"], payload["resolved_schedule"]["seed_count"],
        payload["resolved_schedule"]["sample_count_per_seed"]))


if __name__ == "__main__":
    main()
