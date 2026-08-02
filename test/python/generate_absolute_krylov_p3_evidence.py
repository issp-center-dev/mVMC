from __future__ import print_function

import argparse
import datetime
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys


SCHEMA_VERSION = 1
ORDERS = tuple(range(4))
MATRIX_INDICES = tuple(range(3))
EXPECTED_PROFILE_FIXTURE = "p3_scaling_electronic_real"
EXPECTED_PROFILE_SITES = (4, 6, 8)
EXPECTED_PROFILE_QP = (1, 4)
EXPECTED_PROFILE_RANKS = (1, 2, 4)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def read_json(path):
    with open(path, "r") as handle:
        return json.load(handle)


def write_json(path, value):
    parent = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(parent):
        os.makedirs(parent)
    with open(path, "w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def materialize_raw_files(paths, output, subdirectory):
    destination_root = os.path.join(os.path.abspath(output), subdirectory)
    if not os.path.isdir(destination_root):
        os.makedirs(destination_root)
    result = {}
    names = [os.path.basename(path) for path in paths]
    if len(names) != len(set(names)):
        raise AssertionError("raw evidence basenames are not unique")
    for source, name in zip(paths, names):
        destination = os.path.join(destination_root, name)
        if os.path.abspath(source) != os.path.abspath(destination):
            shutil.copyfile(source, destination)
        if sha256_file(destination) != sha256_file(source):
            raise AssertionError("raw evidence copy checksum mismatch")
        result[os.path.abspath(source)] = os.path.join(subdirectory, name)
    return result


def generated_at_jst():
    zone = datetime.timezone(datetime.timedelta(hours=9))
    return datetime.datetime.now(zone).isoformat(timespec="seconds")


def complex_pair(value):
    return [float(value.real), float(value.imag)]


def finite_complex(value):
    return math.isfinite(value.real) and math.isfinite(value.imag)


def percentile(values, fraction):
    ordered = sorted(float(value) for value in values)
    if not ordered:
        raise AssertionError("percentile of empty sequence")
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def common_metadata(args, fixture_hash):
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at_jst": generated_at_jst(),
        "source_commit": args.source_commit,
        "baseline_develop_commit": args.baseline_commit,
        "compiler": args.compiler,
        "mpi": args.mpi,
        "blas": args.blas,
        "strict_fp": args.strict_fp,
        "omp_threads": args.omp_threads,
        "blas_threads": args.blas_threads,
        "fixture_sha256": fixture_hash,
    }


def validate_policy_contract(policy):
    if policy.get("schema_version") != SCHEMA_VERSION:
        raise AssertionError("measurement policy schema mismatch")
    o2 = policy.get("p3_o2", {})
    if o2.get("lambda") != [1.0] * len(ORDERS):
        raise AssertionError("P3-O2 lambda grid mismatch")
    if o2.get("rho") != [1.0e-6, 1.0e-4, 1.0e-2]:
        raise AssertionError("P3-O2 rho grid mismatch")
    o5 = policy.get("p3_o5", {})
    required = o5.get("required_entries", {})
    if (required.get("families") != ["S", "K", "B"] or
            required.get("i") != list(MATRIX_INDICES) or
            required.get("j") != list(MATRIX_INDICES) or
            int(o5.get("sample_count", 0)) != 4096):
        raise AssertionError("P3-O5 required-entry contract mismatch")
    grid = policy.get("official_scaling_grid", {})
    if (grid.get("site_counts") != list(EXPECTED_PROFILE_SITES) or
            grid.get("qp_total") != list(EXPECTED_PROFILE_QP) or
            grid.get("mpi_rank_counts") != list(EXPECTED_PROFILE_RANKS) or
            int(grid.get("repeat_count", 0)) != 7 or
            int(grid.get("omp_threads", 0)) != 1 or
            int(grid.get("blas_threads", 0)) != 1):
        raise AssertionError("official scaling grid contract mismatch")
    counts = grid.get("sample_count_by_site", {})
    if set(counts) != {str(site) for site in EXPECTED_PROFILE_SITES}:
        raise AssertionError("sample-count grid keys mismatch")
    if any(not isinstance(counts[str(site)], int) or
           counts[str(site)] < 0 for site in EXPECTED_PROFILE_SITES):
        raise AssertionError("invalid sample-count grid")
    return policy


def parse_ed_output(text):
    cases = {}
    current = None
    rank_count = None
    for line in text.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "BEGIN":
            if len(fields) != 6 or fields[2] != "rank_count":
                raise AssertionError("malformed ED BEGIN: {}".format(line))
            current = fields[1]
            rank_count = int(fields[3])
            cases[current] = []
        elif fields[0] == "ROW":
            if current is None or len(fields) != 10:
                raise AssertionError("malformed ED ROW: {}".format(line))
            values = []
            for order in ORDERS:
                value = complex(float(fields[2 + 2 * order]),
                                float(fields[3 + 2 * order]))
                if not finite_complex(value):
                    raise AssertionError("nonfinite ED value")
                values.append(value)
            cases[current].append({
                "configuration": int(fields[1]),
                "values": values,
            })
        elif fields[0] == "END":
            if current is None or fields[1] != current:
                raise AssertionError("malformed ED END")
            current = None
    if current is not None or set(cases) != {
            "electronic_real", "electronic_complex", "spin_real",
            "spin_complex"}:
        raise AssertionError("incomplete ED cases")
    if rank_count != 1:
        raise AssertionError("exact freeze must use the serial ED trace")
    up_masks = [mask for mask in range(1 << 4)
                if bin(mask).count("1") == 2]
    electronic_sector = {
        up | (down << 4) for up in up_masks for down in up_masks
    }
    spin_sector = {
        up | (((~up) & 0x0f) << 4) for up in up_masks
    }
    for name, rows in cases.items():
        observed = [row["configuration"] for row in rows]
        expected = spin_sector if name.startswith("spin_") else electronic_sector
        if len(observed) != len(set(observed)):
            raise AssertionError("duplicate ED configuration in {}".format(name))
        if set(observed) != expected:
            raise AssertionError(
                "ED sector coverage mismatch in {}: missing={} extra={}".format(
                    name, sorted(expected - set(observed)),
                    sorted(set(observed) - expected)))
    return cases


def run_ed_driver(path):
    process = subprocess.run(
        [path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True, check=False,
    )
    if process.returncode != 0:
        raise AssertionError("ED driver failed:\n{}".format(process.stdout))
    return process.stdout


def run_correctness_gate(script, driver):
    environment = os.environ.copy()
    environment.pop("MVMC_MPI_PROCS", None)
    environment.pop("MVMC_MPIEXEC", None)
    environment.pop("MVMC_MPIEXEC_NUMPROC_FLAG", None)
    process = subprocess.run(
        [sys.executable, script, driver],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True, check=False, env=environment,
    )
    expected = (
        "absolute Krylov classic ED passed: 4 families, ranks=1, "
        "orders=0..3, zero-support bridges verified"
    )
    if process.returncode != 0 or expected not in process.stdout:
        raise AssertionError(
            "independent deterministic correctness gate failed:\n{}".format(
                process.stdout))
    return {
        "status": "PASS",
        "script_file": os.path.basename(script),
        "script_sha256": sha256_file(script),
        "stdout_sha256": sha256_text(process.stdout),
        "assertion": expected,
    }


def integrand_indices(family, i, j):
    if family == "S":
        return i, j
    if family == "K":
        return i, j + 1
    if family == "B":
        return i + 1, j + 1
    raise AssertionError("unknown integrand family")


def exact_case_statistics(rows, policy):
    values = [row["values"] for row in rows]
    norms = [sum(abs(row[order]) ** 2 for row in values)
             for order in ORDERS]
    if any(not math.isfinite(value) or value <= 0.0 for value in norms):
        raise AssertionError("zero/nonfinite Krylov norm")
    scales = [1.0] + [math.sqrt(norms[0] / norms[order])
                      for order in ORDERS[1:]]
    scaled = [[scales[order] * row[order] for order in ORDERS]
              for row in values]
    target = [sum(abs(value) ** 2 for value in row) for row in scaled]
    target_sum = sum(target)
    target_mean = target_sum / len(target)
    entry_values = {}
    entry_budget = {}
    matrices = {family: [[None] * 3 for _ in MATRIX_INDICES]
                for family in ("S", "K", "B")}
    for family in ("S", "K", "B"):
        for i in MATRIX_INDICES:
            for j in MATRIX_INDICES:
                left, right = integrand_indices(family, i, j)
                samples = [row[left].conjugate() * row[right]
                           for row in scaled]
                name = "{}_{}{}".format(family, i, j)
                unsigned_scale = sum(abs(value) for value in samples) / target_sum
                budget = max(1.0e-12, 0.02 * unsigned_scale)
                exact_target = sum(samples) / target_sum
                cancellation_denominator = sum(abs(value) for value in samples)
                cancellation = (abs(sum(samples)) / cancellation_denominator
                                if cancellation_denominator != 0.0 else 1.0)
                entry_values[name] = samples
                entry_budget[name] = {
                    "unsigned_scale": unsigned_scale,
                    "absolute_standard_error_budget": budget,
                    "exact_normalized_target": complex_pair(exact_target),
                    "cancellation_ratio": cancellation,
                }
                matrices[family][i][j] = complex_pair(sum(samples))

    guides = {}
    sample_count = int(policy["p3_o5"]["sample_count"])
    for rank in ORDERS:
        for rho in policy["p3_o2"]["rho"]:
            eta = float(rho) * target_mean
            guide = [sum(abs(row[order]) ** 2
                         for order in range(rank + 1)) + eta
                     for row in scaled]
            if any(not math.isfinite(value) or value <= 0.0
                   for value in guide):
                raise AssertionError("guide is not strictly positive")
            guide_sum = sum(guide)
            probabilities = [value / guide_sum for value in guide]
            weights = [target[index] / guide[index]
                       for index in range(len(guide))]
            mean_weight = sum(probabilities[index] * weights[index]
                              for index in range(len(guide)))
            mean_square = sum(probabilities[index] * weights[index] ** 2
                              for index in range(len(guide)))
            ess = mean_weight ** 2 / mean_square
            entry_statistics = {}
            all_entries_pass = True
            cross_term_max = 0.0
            for left in ORDERS:
                for right in ORDERS:
                    cross_term_max = max(
                        cross_term_max,
                        max(abs(row[left].conjugate() * row[right]) /
                            guide[index]
                            for index, row in enumerate(scaled)))
            for name, samples in entry_values.items():
                raw = [samples[index] / guide[index]
                       for index in range(len(guide))]
                raw_mean = sum(probabilities[index] * raw[index]
                               for index in range(len(guide)))
                raw_second = sum(probabilities[index] * abs(raw[index]) ** 2
                                 for index in range(len(guide)))
                normalized = [value * guide_sum / target_sum for value in raw]
                normalized_mean = sum(
                    probabilities[index] * normalized[index]
                    for index in range(len(guide)))
                variance = sum(
                    probabilities[index] *
                    abs(normalized[index] - normalized_mean) ** 2
                    for index in range(len(guide)))
                second_moment = sum(
                    probabilities[index] * abs(normalized[index]) ** 2
                    for index in range(len(guide)))
                standard_error = math.sqrt(max(variance, 0.0) / sample_count)
                budget = entry_budget[name]["absolute_standard_error_budget"]
                passed = standard_error <= budget
                all_entries_pass = all_entries_pass and passed
                tails = [abs(value) for value in raw]
                entry_statistics[name] = {
                    "expectation_f_over_g": complex_pair(raw_mean),
                    "expectation_abs_f_over_g_squared": raw_second,
                    "normalized_expectation": complex_pair(normalized_mean),
                    "normalized_complex_variance": variance,
                    "normalized_second_moment": second_moment,
                    "absolute_standard_error_at_4096": standard_error,
                    "absolute_standard_error_budget": budget,
                    "budget_pass": passed,
                    "cancellation_ratio": entry_budget[name]["cancellation_ratio"],
                    "tail": {
                        "max_abs_f_over_g": max(tails),
                        "median_abs_f_over_g": percentile(tails, 0.5),
                        "p99_abs_f_over_g": percentile(tails, 0.99),
                    },
                    "second_moment_efficiency_per_four_component_evaluations":
                        (1.0 / (4.0 * second_moment)
                         if second_moment > 0.0 else None),
                }
            key = "r{}_rho_{:.0e}".format(rank, float(rho))
            guides[key] = {
                "r": rank,
                "rho": float(rho),
                "eta": eta,
                "guide_min": min(guide),
                "guide_max": max(guide),
                "guide_sum": guide_sum,
                "strict_positive_support": True,
                "t_ess_fraction": ess,
                "t_weight_tail": {
                    "max": max(weights),
                    "median": percentile(weights, 0.5),
                    "p99": percentile(weights, 0.99),
                },
                "cross_term_max_abs_over_g": cross_term_max,
                "entries": entry_statistics,
                "all_entry_budgets_pass": all_entries_pass,
            }
    return {
        "sector_size": len(rows),
        "norms": norms,
        "krylov_scales": scales,
        "sum_t": target_sum,
        "mean_t": target_mean,
        "matrices_unnormalized": matrices,
        "entry_budgets": entry_budget,
        "guides": guides,
    }


def exact_command(args):
    fixture_hash = sha256_file(args.fixture_manifest)
    policy_hash = sha256_file(args.policy)
    fixture = read_json(args.fixture_manifest)
    policy = validate_policy_contract(read_json(args.policy))
    if fixture_hash != policy["fixture_manifest_sha256"]:
        raise AssertionError("policy fixture hash does not match")
    ed_text = run_ed_driver(args.ed_driver)
    cases = parse_ed_output(ed_text)
    correctness_gate = run_correctness_gate(
        args.correctness_script, args.ed_driver)
    metadata = common_metadata(args, fixture_hash)
    metadata["measurement_policy_sha256"] = policy_hash

    reference = dict(metadata)
    reference.update({
        "artifact": "p3_reference_results",
        "fixture_id": fixture["fixture_id"],
        "driver_stdout_sha256": sha256_text(ed_text),
        "deterministic_correctness_gate": correctness_gate,
        "deterministic_correctness_gate_pass": True,
        "cases": {
            name: [{
                "configuration": row["configuration"],
                "v": [complex_pair(value) for value in row["values"]],
            } for row in rows]
            for name, rows in sorted(cases.items())
        },
    })

    analyzed = {name: exact_case_statistics(rows, policy)
                for name, rows in sorted(cases.items())}
    oracle = dict(metadata)
    oracle.update({
        "artifact": "p3_exact_oracle",
        "basis": "u_k=s_k*v_k",
        "matrix_definition": {
            "S_ij": "sum conj(u_i)*u_j",
            "K_ij": "sum conj(u_i)*u_(j+1)",
            "B_ij": "sum conj(u_(i+1))*u_(j+1)",
        },
        "cases": {
            name: {
                key: value for key, value in case.items()
                if key not in ("entry_budgets", "guides")
            } for name, case in analyzed.items()
        },
    })

    model_manifest = dict(metadata)
    model_manifest.update({
        "artifact": "p3_model_manifest",
        "fixture": fixture,
        "measurement_policy": policy,
        "p3_o5_materialized_before_profile_timing": True,
        "entry_budgets": {
            name: case["entry_budgets"] for name, case in analyzed.items()
        },
    })

    candidates = {}
    integrand_cases = {}
    for name, case in analyzed.items():
        integrand_cases[name] = {"guides": case["guides"]}
        for key, guide in case["guides"].items():
            aggregate = candidates.setdefault(key, {
                "r": guide["r"], "rho": guide["rho"],
                "all_cases_entry_budgets_pass": True,
                "all_cases_strict_positive_support": True,
                "all_cases_t_proxy_finite": True,
                "minimum_t_ess_fraction": 1.0,
                "maximum_t_weight_tail": 0.0,
            })
            aggregate["all_cases_entry_budgets_pass"] = (
                aggregate["all_cases_entry_budgets_pass"] and
                guide["all_entry_budgets_pass"])
            aggregate["minimum_t_ess_fraction"] = min(
                aggregate["minimum_t_ess_fraction"],
                guide["t_ess_fraction"])
            aggregate["all_cases_strict_positive_support"] = (
                aggregate["all_cases_strict_positive_support"] and
                guide["strict_positive_support"])
            proxy_values = [
                guide["t_ess_fraction"],
                guide["t_weight_tail"]["max"],
                guide["t_weight_tail"]["median"],
                guide["t_weight_tail"]["p99"],
                guide["cross_term_max_abs_over_g"],
            ]
            aggregate["all_cases_t_proxy_finite"] = (
                aggregate["all_cases_t_proxy_finite"] and
                all(math.isfinite(value) for value in proxy_values))
            aggregate["maximum_t_weight_tail"] = max(
                aggregate["maximum_t_weight_tail"],
                guide["t_weight_tail"]["max"])
    for rank in ORDERS:
        family = [candidate for candidate in candidates.values()
                  if candidate["r"] == rank]
        ess_values = [candidate["minimum_t_ess_fraction"]
                      for candidate in family]
        tail_values = [candidate["maximum_t_weight_tail"]
                       for candidate in family]
        floor_improvement = (
            max(ess_values) > min(ess_values) or
            min(tail_values) < max(tail_values)
        )
        for candidate in family:
            candidate_improves = floor_improvement and (
                candidate["minimum_t_ess_fraction"] > min(ess_values) or
                candidate["maximum_t_weight_tail"] < max(tail_values)
            )
            candidate["floor_scan_improvement_observed"] = \
                candidate_improves
            candidate["p3_exact_feasibility_gate_pass"] = (
                candidate["all_cases_entry_budgets_pass"] and
                candidate["all_cases_strict_positive_support"] and
                candidate["all_cases_t_proxy_finite"] and
                candidate["minimum_t_ess_fraction"] > 0.0 and
                candidate_improves
            )
    integrands = dict(metadata)
    integrands.update({
        "artifact": "p3_integrand_statistics",
        "sample_count_for_budget": policy["p3_o5"]["sample_count"],
        "cases": integrand_cases,
        "candidate_summary": candidates,
        "t_ess_limitation": (
            "aggregate positive-target proxy only; every signed/complex "
            "S/K/B entry is gated separately"
        ),
    })

    write_json(os.path.join(args.output, "p3_reference_results.json"),
               reference)
    write_json(os.path.join(args.output, "p3_exact_oracle.json"), oracle)
    write_json(os.path.join(args.output, "p3_model_manifest.json"),
               model_manifest)
    write_json(os.path.join(args.output, "p3_integrand_statistics.json"),
               integrands)
    print("frozen exact evidence and entry-specific P3-O5 budgets")


def parse_key_values(line):
    result = {}
    for token in line.split()[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            result[key] = value
    return result


def parse_number(value):
    if "," in value:
        return [parse_number(part) for part in value.split(",")]
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def parse_profile_log(path):
    with open(path, "r") as handle:
        text = handle.read()
    header = None
    rows = []
    resources = {}
    for line in text.splitlines():
        if line.startswith("PROFILE "):
            if header is not None:
                raise AssertionError("duplicate PROFILE header")
            header = {key: parse_number(value)
                      for key, value in parse_key_values(line).items()}
        elif line.startswith("ROW "):
            rows.append(line)
        elif line.startswith("RESOURCE "):
            fields = parse_key_values(line)
            scope = fields.pop("scope")
            resources[scope] = {key: parse_number(value)
                                for key, value in fields.items()}
    if header is None or not rows or set(resources) != {"rank_sum", "rank_max"}:
        raise AssertionError("incomplete profile log {}".format(path))
    if int(header["audit"]) != 1 or len(rows) != int(header["sample_count"]):
        raise AssertionError("invalid profile audit trace")
    return {
        "path": os.path.abspath(path),
        "sha256": sha256_file(path),
        "header": header,
        "rows": rows,
        "trace_sha256": sha256_text("\n".join(rows) + "\n"),
        "resources": resources,
    }


def fixed_sector_size(site_count):
    half = site_count // 2
    one_spin = math.factorial(site_count) // (
        math.factorial(half) * math.factorial(site_count - half))
    return one_spin * one_spin


def fixed_masks(site_count):
    half = site_count // 2
    return [mask for mask in range(1 << site_count)
            if bin(mask).count("1") == half]


def expected_profile_rows(site_count, sample_count):
    masks = fixed_masks(site_count)
    sector_size = len(masks) * len(masks)
    rows = []
    for sample in range(sample_count):
        sector_index = ((sector_size - 1) // 2 if sample_count == 1 else
                        sample * (sector_size - 1) // (sample_count - 1))
        up_index = sector_index // len(masks)
        down_index = sector_index % len(masks)
        configuration = (masks[up_index] |
                         (masks[down_index] << site_count))
        rows.append((sample, sector_index, configuration))
    return rows


def validate_profile_run(run, policy):
    header = run["header"]
    required = {
        "schema", "fixture", "site_count", "qp_total", "sample_count",
        "sector_size", "rank_count", "audit",
    }
    missing = required - set(header)
    if missing:
        raise AssertionError("profile header fields missing: {}".format(
            sorted(missing)))
    if (int(header["schema"]) != SCHEMA_VERSION or
            header["fixture"] != EXPECTED_PROFILE_FIXTURE or
            int(header["audit"]) != 1):
        raise AssertionError("profile schema/fixture/audit mismatch")
    site = int(header["site_count"])
    qp = int(header["qp_total"])
    ranks = int(header["rank_count"])
    if (site not in policy["official_scaling_grid"]["site_counts"] or
            qp not in policy["official_scaling_grid"]["qp_total"] or
            ranks not in
            policy["official_scaling_grid"]["mpi_rank_counts"]):
        raise AssertionError("profile header is outside the frozen grid")
    expected_sector = fixed_sector_size(site)
    requested = int(
        policy["official_scaling_grid"]["sample_count_by_site"][str(site)])
    expected_samples = expected_sector if requested == 0 else requested
    if (int(header["sector_size"]) != expected_sector or
            int(header["sample_count"]) != expected_samples):
        raise AssertionError(
            "profile sector/sample contract mismatch for L{}".format(site))
    observed_rows = []
    for line in run["rows"]:
        fields = line.split()
        if len(fields) != 12 or fields[0] != "ROW":
            raise AssertionError("malformed profile ROW")
        try:
            identifiers = tuple(int(value) for value in fields[1:4])
            values = [float(value) for value in fields[4:]]
        except ValueError:
            raise AssertionError("non-numeric profile ROW")
        if any(not math.isfinite(value) for value in values):
            raise AssertionError("nonfinite profile ROW value")
        observed_rows.append(identifiers)
    if observed_rows != expected_profile_rows(site, expected_samples):
        raise AssertionError("profile sample/configuration schedule mismatch")
    rank_sum = run["resources"]["rank_sum"]
    rank_max = run["resources"]["rank_max"]
    if (int(rank_max.get("roots", -1)) != expected_samples or
            int(rank_sum.get("roots", -1)) != expected_samples * ranks):
        raise AssertionError("profile root-count invariant mismatch")
    logical = (
        "regular_logical", "near_pivot_logical", "singular_logical",
        "total_zero_logical", "global_factorizations_logical",
    )
    for field in logical:
        if (field not in run["resources"]["rank_sum"] or
                field not in run["resources"]["rank_max"] or
                run["resources"]["rank_sum"][field] !=
                run["resources"]["rank_max"][field]):
            raise AssertionError(
                "communicator-global logical counter mismatch: {}".format(
                    field))
    if ("local_factorizations" not in rank_sum or
            int(rank_sum["local_factorizations"]) !=
            int(rank_sum["global_factorizations_logical"])):
        raise AssertionError("invalid local/global factorization counters")


def log_linear_fit(points, target_l):
    usable = [(float(site), float(value)) for site, value in points
              if value > 0.0 and math.isfinite(value)]
    if len(usable) != 3:
        raise AssertionError("L=4/6/8 positive values required for fit")
    mean_x = sum(point[0] for point in usable) / len(usable)
    logs = [math.log(point[1]) for point in usable]
    mean_y = sum(logs) / len(logs)
    denominator = sum((point[0] - mean_x) ** 2 for point in usable)
    slope = sum((point[0] - mean_x) * (logs[index] - mean_y)
                for index, point in enumerate(usable)) / denominator
    intercept = mean_y - slope * mean_x
    prediction = math.exp(intercept + slope * float(target_l))
    return {
        "model": "log(value)=intercept+slope*L",
        "intercept": intercept,
        "slope": slope,
        "prediction": prediction,
        "points": [[int(site), value] for site, value in usable],
    }


def extrapolate_nonnegative(points, target_l):
    if all(float(value) == 0.0 for _, value in points):
        return {
            "model": "structural_zero",
            "intercept": None,
            "slope": 0.0,
            "prediction": 0.0,
            "points": [[int(site), float(value)] for site, value in points],
        }
    return log_linear_fit(points, target_l)


def profile_command(args):
    policy = validate_policy_contract(read_json(args.policy))
    model_manifest = read_json(args.model_manifest)
    integrands = read_json(args.integrand_statistics)
    reference = read_json(args.reference_results)
    production_anchor = read_json(args.production_anchor)
    policy_hash = sha256_file(args.policy)
    if model_manifest["measurement_policy_sha256"] != policy_hash:
        raise AssertionError("profile policy differs from frozen manifest")
    runs = [parse_profile_log(path) for path in args.profile_log]
    for run in runs:
        validate_profile_run(run, policy)
    artifacts = (model_manifest, integrands, reference, production_anchor)
    expected_fixture_hash = model_manifest["fixture_sha256"]
    for artifact in artifacts:
        if (artifact.get("source_commit") != args.source_commit or
                artifact.get("baseline_develop_commit") !=
                args.baseline_commit or
                artifact.get("fixture_sha256") != expected_fixture_hash or
                artifact.get("measurement_policy_sha256") != policy_hash):
            raise AssertionError("cross-artifact provenance mismatch")
    if not reference.get("deterministic_correctness_gate_pass", False):
        raise AssertionError("deterministic correctness gate is not PASS")
    anchor_runs = production_anchor.get("runs", [])
    if (not production_anchor.get("all_runs_pass", False) or
            len(anchor_runs) != len(EXPECTED_PROFILE_RANKS[:2]) or
            sorted(run.get("rank_count") for run in anchor_runs) != [1, 2] or
            any(run.get("status") != "PASS" for run in anchor_runs)):
        raise AssertionError("serial/MPI production anchor gate is not PASS")
    groups = {}
    for run in runs:
        header = run["header"]
        key = (int(header["site_count"]), int(header["qp_total"]),
               int(header["rank_count"]))
        groups.setdefault(key, []).append(run)
    expected_repeats = int(policy["official_scaling_grid"]["repeat_count"])
    expected_keys = {(site, qp, rank)
                     for site in policy["official_scaling_grid"]["site_counts"]
                     for qp in policy["official_scaling_grid"]["qp_total"]
                     for rank in policy["official_scaling_grid"]["mpi_rank_counts"]}
    if set(groups) != expected_keys:
        raise AssertionError("profile grid mismatch: missing={} extra={}".format(
            sorted(expected_keys - set(groups)),
            sorted(set(groups) - expected_keys)))

    cross_rank_trace_invariance = {}
    for site in policy["official_scaling_grid"]["site_counts"]:
        for qp in policy["official_scaling_grid"]["qp_total"]:
            traces = {
                member["trace_sha256"]
                for rank in policy["official_scaling_grid"]["mpi_rank_counts"]
                for member in groups[(site, qp, rank)]
            }
            if len(traces) != 1:
                raise AssertionError(
                    "MPI rank-count trace mismatch for L{} qp{}".format(
                        site, qp))
            cross_rank_trace_invariance["L{}_qp{}".format(site, qp)] = {
                "rank_counts":
                    policy["official_scaling_grid"]["mpi_rank_counts"],
                "trace_sha256": next(iter(traces)),
                "status": "PASS",
            }

    summaries = {}
    for key, members in sorted(groups.items()):
        if len(members) != expected_repeats:
            raise AssertionError("repeat count mismatch for {}".format(key))
        traces = {member["trace_sha256"] for member in members}
        if len(traces) != 1:
            raise AssertionError("repeat trace mismatch for {}".format(key))
        roots = [member["resources"]["rank_max"]["roots"]
                 for member in members]
        depth_per_root = []
        for member in members:
            resource = member["resources"]["rank_max"]
            depth_per_root.append([
                float(value) / float(resource["roots"])
                for value in resource["depth_seconds"]
            ])
        median_depth = [statistics.median(row[order]
                                          for row in depth_per_root)
                        for order in ORDERS]
        summaries["L{}_qp{}_rank{}".format(*key)] = {
            "site_count": key[0],
            "qp_total": key[1],
            "rank_count": key[2],
            "repeat_count": len(members),
            "sample_count": members[0]["header"]["sample_count"],
            "sector_size": members[0]["header"]["sector_size"],
            "trace_sha256": members[0]["trace_sha256"],
            "roots_per_rank": statistics.median(roots),
            "median_depth_seconds_per_root": median_depth,
            "median_total_seconds_per_root": statistics.median(
                float(member["resources"]["rank_max"]["total_seconds"]) /
                float(member["resources"]["rank_max"]["roots"])
                for member in members),
            "median_workspace_bytes_per_rank": statistics.median(
                member["resources"]["rank_max"]["krylov_workspace_bytes"] +
                member["resources"]["rank_max"]["model_workspace_bytes"] +
                member["resources"]["rank_max"]["amplitude_workspace_bytes"]
                for member in members),
            "median_rss_bytes_per_rank": statistics.median(
                member["resources"]["rank_max"]["rss_bytes"]
                for member in members),
            "rank_sum_median": {
                key_name: statistics.median(
                    member["resources"]["rank_sum"][key_name]
                    for member in members)
                for key_name in (
                    "roots", "raw_transitions", "merged_duplicates",
                    "cancelled_zero", "terminal_requests",
                    "terminal_cache_hits", "local_factorizations",
                    "global_factorizations_logical")
            },
            "rank_max_peaks": {
                key_name: max(member["resources"]["rank_max"][key_name]
                              for member in members)
                for key_name in ("frontier_peak", "memo_peak")
            },
        }

    fixture_hash = model_manifest["fixture_sha256"]
    metadata = common_metadata(args, fixture_hash)
    metadata["measurement_policy_sha256"] = policy_hash
    raw_profile_files = materialize_raw_files(
        [run["path"] for run in runs], args.output, "raw_profiles")
    resources = dict(metadata)
    resources.update({
        "artifact": "p3_resource_statistics",
        "official_hpc_measurement": True,
        "counter_semantics": {
            "rank_sum": "sum of per-rank work counters",
            "rank_max": "maximum per-rank work/resource counter",
            "logical_component_and_global_factorization_counters": (
                "communicator-global logical counts, recorded once rather "
                "than summing replicated callback results"
            ),
        },
        "raw_logs": [{
            "file": raw_profile_files[run["path"]],
            "sha256": run["sha256"],
            "header": run["header"],
            "trace_sha256": run["trace_sha256"],
            "rank_sum": run["resources"]["rank_sum"],
            "rank_max": run["resources"]["rank_max"],
        } for run in runs],
        "cross_rank_trace_invariance": cross_rank_trace_invariance,
        "group_summary": summaries,
    })

    extrapolations = {}
    o1 = policy["p3_o1"]
    budgets = o1["l16_budgets"]
    workload = o1["representative_workload"]
    guide_passes_all_grid_points = {guide_rank: True for guide_rank in ORDERS}
    for qp in policy["official_scaling_grid"]["qp_total"]:
        for rank_count in policy["official_scaling_grid"]["mpi_rank_counts"]:
            combo = "qp{}_rank{}".format(qp, rank_count)
            combo_result = {
                "qp_total": qp,
                "rank_count": rank_count,
                "guides": {},
            }
            workspace_points = []
            rss_points = []
            for site in (4, 6, 8):
                summary = summaries["L{}_qp{}_rank{}".format(
                    site, qp, rank_count)]
                workspace_points.append(
                    (site, summary["median_workspace_bytes_per_rank"]))
                rss_points.append((site, summary["median_rss_bytes_per_rank"]))
            workspace_l16 = log_linear_fit(workspace_points, 16)
            rss_l16 = log_linear_fit(rss_points, 16)
            workspace_l32 = log_linear_fit(workspace_points, 32)
            rss_l32 = log_linear_fit(rss_points, 32)
            combo_result["rank_memory"] = {
                "l16_workspace": workspace_l16,
                "l16_rss": rss_l16,
                "l16_conservative_bytes": max(
                    workspace_l16["prediction"], rss_l16["prediction"]),
                "l32_workspace_sensitivity": workspace_l32,
                "l32_rss_sensitivity": rss_l32,
            }
            for guide_rank in ORDERS:
                proposal_points = []
                saved_points = []
                full_points = []
                for site in (4, 6, 8):
                    depth = summaries[
                        "L{}_qp{}_rank{}".format(site, qp, rank_count)
                    ]["median_depth_seconds_per_root"]
                    proposal_points.append((site, sum(depth[:guide_rank + 1])))
                    saved_points.append((site, sum(depth[guide_rank + 1:])))
                    full_points.append((site, sum(depth)))
                proposal_l16 = log_linear_fit(proposal_points, 16)
                saved_l16 = extrapolate_nonnegative(saved_points, 16)
                full_l16 = log_linear_fit(full_points, 16)
                full_l32 = log_linear_fit(full_points, 32)
                combined = (
                    float(workload["proposal_count"]) *
                    proposal_l16["prediction"] +
                    float(workload["saved_sample_count"]) *
                    saved_l16["prediction"])
                memory = combo_result["rank_memory"]["l16_conservative_bytes"]
                passed = (
                    proposal_l16["prediction"] <=
                    float(budgets["proposal_seconds_per_configuration"]) and
                    saved_l16["prediction"] <=
                    float(budgets["saved_increment_seconds_per_configuration"]) and
                    combined <=
                    float(budgets["combined_workload_seconds_per_rank"]) and
                    memory <= float(budgets["rank_rss_bytes"])
                )
                guide_passes_all_grid_points[guide_rank] = (
                    guide_passes_all_grid_points[guide_rank] and passed)
                combo_result["guides"]["r{}".format(guide_rank)] = {
                    "proposal_component_evaluation_count": guide_rank + 1,
                    "saved_increment_component_evaluation_count":
                        ORDERS[-1] - guide_rank,
                    "accepted_current_guide_re_evaluation_count": 0,
                    "proposal_l16": proposal_l16,
                    "saved_increment_l16": saved_l16,
                    "full_l16": full_l16,
                    "combined_workload_seconds_l16": combined,
                    "full_l32_sensitivity": full_l32,
                    "o1_budget_pass": passed,
                }
            extrapolations[combo] = combo_result

    viable_r = {guide_rank for guide_rank, passed
                in guide_passes_all_grid_points.items() if passed}
    exact_candidates = integrands["candidate_summary"]
    correctness_gate = {
        "independent_deterministic_ed": {
            "status": "PASS",
            "script_sha256": reference["deterministic_correctness_gate"]
            ["script_sha256"],
            "stdout_sha256": reference["deterministic_correctness_gate"]
            ["stdout_sha256"],
        },
        "production_action_anchor": {
            "status": "PASS",
            "rank_counts": sorted(run["rank_count"] for run in anchor_runs),
            "log_sha256": sorted(run["sha256"] for run in anchor_runs),
        },
        "profile_cross_rank_trace_invariance": {
            "status": "PASS",
            "groups": len(cross_rank_trace_invariance),
        },
        "all_pass": True,
    }
    promoted = []
    for name, candidate in sorted(exact_candidates.items()):
        if (correctness_gate["all_pass"] and
                candidate["p3_exact_feasibility_gate_pass"] and
                int(candidate["r"]) in viable_r):
            promoted.append(name)
    decision = "GO" if promoted else "STOP"
    cost = dict(metadata)
    cost.update({
        "artifact": "p3_cost_feasibility",
        "p3_o1": o1,
        "extrapolation": extrapolations,
        "o1_viable_guide_r_across_all_grid_points": sorted(viable_r),
        "o5_candidate_summary": exact_candidates,
        "correctness_gate": correctness_gate,
        "promoted_candidates": promoted,
        "p4_decision": decision,
        "decision_reason": (
            "at least one candidate passes exact integrand and measured "
            "L=16 cost/RSS gates" if promoted else
            "no candidate passes both exact integrand and conservative "
            "L=16 cost/RSS gates; redesign cache/workspace/API before P4"
        ),
        "p3_o3_api": (
            "promote_testing_value_cache_contract" if promoted else
            "not_fixed_reference_workspace_requires_redesign"
        ),
    })

    l4_cost = summaries["L4_qp1_rank1"]["median_depth_seconds_per_root"]
    for case in integrands["cases"].values():
        for guide in case["guides"].values():
            full_cost = sum(l4_cost)
            for entry in guide["entries"].values():
                moment = entry["normalized_second_moment"]
                entry["hpc_second_moment_efficiency_per_second"] = (
                    1.0 / (moment * full_cost)
                    if moment > 0.0 and full_cost > 0.0 else None)
    integrands["hpc_cost_basis"] = {
        "group": "L4_qp1_rank1",
        "median_full_seconds_per_root": sum(l4_cost),
    }
    integrands["compiler"] = args.compiler
    integrands["mpi"] = args.mpi
    integrands["blas"] = args.blas
    integrands["generated_at_jst"] = generated_at_jst()

    write_json(os.path.join(args.output, "p3_resource_statistics.json"),
               resources)
    write_json(os.path.join(args.output, "p3_cost_feasibility.json"), cost)
    write_json(args.integrand_statistics, integrands)
    print("profile evidence complete: P4 {}".format(decision))


def anchor_command(args):
    fixture_hash = sha256_file(args.fixture_manifest)
    policy_hash = sha256_file(args.policy)
    policy = validate_policy_contract(read_json(args.policy))
    if fixture_hash != policy["fixture_manifest_sha256"]:
        raise AssertionError("anchor policy fixture hash does not match")
    pattern = re.compile(
        r"P3 production action anchor passed: 4 cases, "
        r"([0-9]+) sampled rows, ([0-9]+) case-local distinct states"
    )
    if len(args.anchor_log) != len(args.anchor_rank_count):
        raise AssertionError("anchor log/rank-count cardinality mismatch")
    if sorted(args.anchor_rank_count) != [1, 2]:
        raise AssertionError("production anchor requires rank counts 1 and 2")
    runs = []
    for path, rank_count in zip(args.anchor_log, args.anchor_rank_count):
        with open(path, "r") as handle:
            text = handle.read()
        match = pattern.search(text)
        if match is None:
            raise AssertionError("anchor PASS summary missing: {}".format(path))
        rows = int(match.group(1))
        states = int(match.group(2))
        if rows != 16384 * rank_count or states != 70:
            raise AssertionError("anchor coverage count mismatch")
        runs.append({
            "file": None,
            "source_path": os.path.abspath(path),
            "sha256": sha256_file(path),
            "case_count": 4,
            "rank_count": rank_count,
            "sampled_rows": rows,
            "sampled_rows_per_rank": rows // rank_count,
            "case_local_distinct_states": states,
            "status": "PASS",
        })
    raw_anchor_files = materialize_raw_files(
        [run["source_path"] for run in runs], args.output, "raw_anchor")
    for run in runs:
        run["file"] = raw_anchor_files[run.pop("source_path")]
    metadata = common_metadata(args, fixture_hash)
    metadata["measurement_policy_sha256"] = policy_hash
    artifact = dict(metadata)
    artifact.update({
        "artifact": "p3_production_anchor",
        "comparison": (
            "production F_k versus absolute-reference v_k/v_0 on exact "
            "nonzero support"
        ),
        "runs": runs,
        "all_runs_pass": True,
    })
    write_json(os.path.join(args.output, "p3_production_anchor.json"),
               artifact)
    print("production action anchor evidence complete")


def copy_if_distinct(source, destination):
    if os.path.abspath(source) != os.path.abspath(destination):
        shutil.copyfile(source, destination)


def validate_commit_object(repository, commit):
    if re.match(r"^[0-9a-f]{40}$", commit or "") is None:
        raise AssertionError("commit provenance must be a full lowercase SHA")
    process = subprocess.run(
        ["git", "-C", repository, "rev-parse", "--verify",
         commit + "^{commit}"],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True, check=False,
    )
    if process.returncode != 0 or process.stdout.strip() != commit:
        raise AssertionError("commit provenance is not a repository commit")


def validate_referenced_raw(output, records, expected_count, label):
    if len(records) != expected_count:
        raise AssertionError("{} raw evidence count mismatch".format(label))
    names = [record.get("file") for record in records]
    if len(names) != len(set(names)) or any(not name for name in names):
        raise AssertionError("{} raw evidence names are invalid".format(label))
    output_real = os.path.realpath(output)
    for record, relative in zip(records, names):
        if os.path.isabs(relative):
            raise AssertionError("raw evidence path must be relative")
        path = os.path.realpath(os.path.join(output_real, relative))
        if (os.path.commonpath([output_real, path]) != output_real or
                not os.path.isfile(path) or os.path.islink(path) or
                sha256_file(path) != record.get("sha256")):
            raise AssertionError(
                "{} raw evidence is missing or changed: {}".format(
                    label, relative))


def finalize_command(args):
    required_artifacts = {
        "p3_reference_results.json": "p3_reference_results",
        "p3_model_manifest.json": "p3_model_manifest",
        "p3_resource_statistics.json": "p3_resource_statistics",
        "p3_cost_feasibility.json": "p3_cost_feasibility",
        "p3_exact_oracle.json": "p3_exact_oracle",
        "p3_production_anchor.json": "p3_production_anchor",
        "p3_integrand_statistics.json": "p3_integrand_statistics",
    }
    output = os.path.abspath(args.output)
    if not os.path.isdir(output):
        raise AssertionError("evidence output directory is missing")
    validate_commit_object(args.repository, args.source_commit)
    validate_commit_object(args.repository, args.baseline_commit)
    if (not os.path.isfile(args.build_config) or
            os.path.getsize(args.build_config) == 0 or
            not os.path.isfile(args.review) or
            os.path.getsize(args.review) == 0):
        raise AssertionError("build config and review must be nonempty files")
    fixture_hashes = set()
    policy_hashes = set()
    for filename, artifact_name in sorted(required_artifacts.items()):
        path = os.path.join(output, filename)
        artifact = read_json(path)
        if (artifact.get("schema_version") != SCHEMA_VERSION or
                artifact.get("artifact") != artifact_name or
                artifact.get("source_commit") != args.source_commit or
                artifact.get("baseline_develop_commit") !=
                args.baseline_commit or
                not artifact.get("compiler") or
                not artifact.get("mpi") or
                not artifact.get("blas") or
                int(artifact.get("omp_threads", 0)) != 1 or
                int(artifact.get("blas_threads", 0)) != 1 or
                not artifact.get("strict_fp")):
            raise AssertionError(
                "final evidence metadata mismatch: {}".format(filename))
        fixture_hashes.add(artifact.get("fixture_sha256"))
        policy_hashes.add(artifact.get("measurement_policy_sha256"))
    if len(fixture_hashes) != 1 or None in fixture_hashes or \
            len(policy_hashes) != 1 or None in policy_hashes:
        raise AssertionError("final evidence provenance hashes disagree")
    reference = read_json(os.path.join(output, "p3_reference_results.json"))
    anchor = read_json(os.path.join(output, "p3_production_anchor.json"))
    cost = read_json(os.path.join(output, "p3_cost_feasibility.json"))
    resources = read_json(
        os.path.join(output, "p3_resource_statistics.json"))
    policy = read_json(os.path.join(output, "p3_model_manifest.json"))[
        "measurement_policy"]
    grid = policy["official_scaling_grid"]
    expected_profile_logs = (
        len(grid["site_counts"]) * len(grid["qp_total"]) *
        len(grid["mpi_rank_counts"]) * int(grid["repeat_count"]))
    validate_referenced_raw(
        output, resources.get("raw_logs", []), expected_profile_logs,
        "profile")
    validate_referenced_raw(output, anchor.get("runs", []), 2, "anchor")
    if (not reference.get("deterministic_correctness_gate_pass", False) or
            not anchor.get("all_runs_pass", False) or
            not cost.get("correctness_gate", {}).get("all_pass", False) or
            cost.get("p4_decision") not in ("GO", "STOP")):
        raise AssertionError("final correctness/decision gate is incomplete")
    copy_if_distinct(args.build_config,
                     os.path.join(output, "build_config.txt"))
    copy_if_distinct(args.review, os.path.join(output, "review.md"))
    if (os.path.getsize(os.path.join(output, "build_config.txt")) == 0 or
            os.path.getsize(os.path.join(output, "review.md")) == 0):
        raise AssertionError("materialized build config/review is empty")
    with open(os.path.join(output, "source_commit.txt"), "w") as handle:
        handle.write(args.source_commit + "\n")
    with open(os.path.join(output, "baseline_develop_commit.txt"), "w") as handle:
        handle.write(args.baseline_commit + "\n")
    checksum_path = os.path.join(output, "checksums.sha256")
    entries = []
    for root, directories, files in os.walk(output):
        directories.sort()
        for filename in sorted(files):
            path = os.path.join(root, filename)
            if os.path.abspath(path) == checksum_path:
                continue
            if os.path.islink(path):
                raise AssertionError("symlink is not allowed in evidence ledger")
            relative = os.path.relpath(path, output)
            entries.append((relative, sha256_file(path)))
    with open(checksum_path, "w") as handle:
        for relative, digest in sorted(entries):
            handle.write("{}  {}\n".format(digest, relative))
    print("finalized P3 evidence ledger: {} files".format(len(entries)))


def add_common_arguments(parser):
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--baseline-commit", required=True)
    parser.add_argument("--compiler", required=True)
    parser.add_argument("--mpi", required=True)
    parser.add_argument("--blas", required=True)
    parser.add_argument("--strict-fp", required=True)
    parser.add_argument("--omp-threads", type=int, default=1)
    parser.add_argument("--blas-threads", type=int, default=1)


def build_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    exact = subparsers.add_parser("exact")
    add_common_arguments(exact)
    exact.add_argument("--ed-driver", required=True)
    exact.add_argument("--correctness-script", required=True)
    exact.add_argument("--fixture-manifest", required=True)
    exact.add_argument("--policy", required=True)
    exact.add_argument("--output", required=True)
    exact.set_defaults(function=exact_command)

    profile = subparsers.add_parser("profiles")
    add_common_arguments(profile)
    profile.add_argument("--policy", required=True)
    profile.add_argument("--model-manifest", required=True)
    profile.add_argument("--integrand-statistics", required=True)
    profile.add_argument("--reference-results", required=True)
    profile.add_argument("--production-anchor", required=True)
    profile.add_argument("--profile-log", action="append", required=True)
    profile.add_argument("--output", required=True)
    profile.set_defaults(function=profile_command)

    anchor = subparsers.add_parser("anchor")
    add_common_arguments(anchor)
    anchor.add_argument("--fixture-manifest", required=True)
    anchor.add_argument("--policy", required=True)
    anchor.add_argument("--anchor-log", action="append", required=True)
    anchor.add_argument("--anchor-rank-count", action="append", type=int,
                        required=True)
    anchor.add_argument("--output", required=True)
    anchor.set_defaults(function=anchor_command)

    finalize = subparsers.add_parser("finalize")
    finalize.add_argument("--source-commit", required=True)
    finalize.add_argument("--baseline-commit", required=True)
    finalize.add_argument("--repository", required=True)
    finalize.add_argument("--build-config", required=True)
    finalize.add_argument("--review", required=True)
    finalize.add_argument("--output", required=True)
    finalize.set_defaults(function=finalize_command)
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    if not hasattr(args, "function"):
        parser.error("a subcommand is required")
    args.function(args)


if __name__ == "__main__":
    main()
