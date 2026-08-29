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
EXPECTED_PROFILE_FIXTURE = "p4c_bounded_electronic_real"
OFFICIAL_SITE_COUNTS = [4, 6, 8]
OFFICIAL_TARGET_SITE_COUNT = 16
OFFICIAL_QP_TOTAL = [1, 4]
OFFICIAL_RANK_COUNTS = [1, 2, 4]
OFFICIAL_REPEAT_COUNT = 7
OFFICIAL_CACHE_BYTES = [0, 67108864, 268435456]
OFFICIAL_RESOURCE_SAMPLE_COUNT_BY_SITE = {"4": 0, "6": 128, "8": 128}
OFFICIAL_RHO = [1.0e-6, 1.0e-4, 1.0e-2]


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
    if parent and not os.path.isdir(parent):
        os.makedirs(parent)
    with open(path, "w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


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


def median(values):
    values = list(values)
    if not values:
        raise AssertionError("median of empty sequence")
    return statistics.median(values)


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


def parse_key_values_from_tokens(tokens):
    result = {}
    for token in tokens:
        if "=" in token:
            key, value = token.split("=", 1)
            result[key] = value
    return result


def parse_key_values(line, skip):
    return parse_key_values_from_tokens(line.split()[skip:])


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


def log_linear_fit(points, target_l):
    usable = [(float(site), float(value)) for site, value in points
              if value > 0.0 and math.isfinite(value)]
    if len(usable) != len(points) or len(usable) < 2:
        raise AssertionError("positive finite values required for fit")
    mean_x = sum(point[0] for point in usable) / len(usable)
    logs = [math.log(point[1]) for point in usable]
    mean_y = sum(logs) / len(logs)
    denominator = sum((point[0] - mean_x) ** 2 for point in usable)
    if denominator == 0.0:
        raise AssertionError("fit requires at least two distinct site counts")
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


def validate_cache_values(values, allow_smoke_policy):
    if not isinstance(values, list) or not values:
        raise AssertionError("cache grid must be a nonempty list")
    if any(not isinstance(value, int) or value < 0 for value in values):
        raise AssertionError("invalid cache capacity")
    if allow_smoke_policy:
        return
    if values != OFFICIAL_CACHE_BYTES:
        raise AssertionError("official cache grid mismatch")


def normalize_sample_count_by_site(scope):
    values = scope.get("sample_count_by_site")
    site_counts = scope.get("site_counts")
    if not isinstance(values, dict):
        raise AssertionError("sample_count_by_site must be a map")
    normalized = {}
    for key, value in values.items():
        try:
            site = int(key)
        except (TypeError, ValueError):
            raise AssertionError("invalid sample-count site key")
        if str(site) != str(key):
            raise AssertionError("sample-count site keys must be canonical")
        if not isinstance(value, int) or value < 0:
            raise AssertionError("invalid resource sample count")
        sector_size = fixed_sector_size(site)
        if value > sector_size:
            raise AssertionError("resource sample count exceeds sector")
        normalized[site] = value
    if set(normalized) != set(site_counts):
        raise AssertionError("sample-count map must match site grid")
    return normalized


def expected_resource_sample_count(policy, site_count):
    requested = normalize_sample_count_by_site(policy["scope"])[site_count]
    return fixed_sector_size(site_count) if requested == 0 else requested


def validate_policy_contract(policy, allow_smoke_policy=False):
    if policy.get("schema_version") != SCHEMA_VERSION:
        raise AssertionError("P4 policy schema mismatch")
    scope = policy.get("scope", {})
    guide = policy.get("guide", {})
    gate = policy.get("p4_c_engine_feasibility", {})
    budgets = gate.get("target_budgets", {})
    workload = gate.get("representative_workload", {})
    exact = gate.get("exact_iid_o5", {})
    site_counts = scope.get("site_counts")
    target_site_count = scope.get("target_site_count")
    qp_total = scope.get("qp_total")
    rank_counts = scope.get("mpi_rank_counts")
    repeat_count = scope.get("repeat_count")
    cache_values = scope.get("cache_capacity_bytes")
    sample_count_by_site = normalize_sample_count_by_site(scope)
    if not all(isinstance(value, int) and value > 0 for value in site_counts):
        raise AssertionError("invalid site grid")
    if not isinstance(target_site_count, int) or target_site_count <= 0:
        raise AssertionError("invalid target site count")
    if not all(isinstance(value, int) and value in (1, 4)
               for value in qp_total):
        raise AssertionError("invalid QP grid")
    if not all(isinstance(value, int) and value > 0 for value in rank_counts):
        raise AssertionError("invalid rank grid")
    if not isinstance(repeat_count, int) or repeat_count <= 0:
        raise AssertionError("invalid repeat count")
    validate_cache_values(cache_values, allow_smoke_policy)
    if (int(scope.get("omp_threads", 0)) != 1 or
            int(scope.get("blas_threads", 0)) != 1):
        raise AssertionError("P4-C policy requires OMP/BLAS threads 1")
    if guide.get("candidate_family") != "r3":
        raise AssertionError("P4-C only promotes r=3")
    if guide.get("lambda") != [1.0, 1.0, 1.0, 1.0]:
        raise AssertionError("P4-C lambda contract mismatch")
    if [float(value) for value in guide.get("rho", [])] != OFFICIAL_RHO:
        raise AssertionError("P4-C rho grid mismatch")
    required_budget_keys = (
        "proposal_seconds_per_configuration",
        "saved_increment_seconds_per_configuration",
        "combined_workload_seconds_per_rank",
        "allocated_capacity_bytes_per_rank",
        "peak_rss_bytes_per_rank",
    )
    for key in required_budget_keys:
        if key not in budgets or float(budgets[key]) <= 0.0:
            raise AssertionError("missing/invalid budget {}".format(key))
    for key in ("proposal_count", "saved_sample_count"):
        if key not in workload or int(workload[key]) <= 0:
            raise AssertionError("missing/invalid workload {}".format(key))
    if int(exact.get("sample_count_proxy", 0)) != 4096:
        raise AssertionError("P4-C exact IID sample-count proxy mismatch")
    if allow_smoke_policy:
        if not set(site_counts).issubset(set(OFFICIAL_SITE_COUNTS)):
            raise AssertionError("smoke site grid must be an official subset")
        if not set(qp_total).issubset(set(OFFICIAL_QP_TOTAL)):
            raise AssertionError("smoke QP grid must be an official subset")
        if not set(rank_counts).issubset(set(OFFICIAL_RANK_COUNTS)):
            raise AssertionError("smoke rank grid must be an official subset")
        if repeat_count > OFFICIAL_REPEAT_COUNT:
            raise AssertionError("smoke repeat count exceeds official count")
    else:
        official_samples = {
            int(key): value
            for key, value in OFFICIAL_RESOURCE_SAMPLE_COUNT_BY_SITE.items()
        }
        if (site_counts != OFFICIAL_SITE_COUNTS or
                target_site_count != OFFICIAL_TARGET_SITE_COUNT or
                qp_total != OFFICIAL_QP_TOTAL or
                rank_counts != OFFICIAL_RANK_COUNTS or
                repeat_count != OFFICIAL_REPEAT_COUNT or
                sample_count_by_site != official_samples):
            raise AssertionError("official P4-C grid mismatch")
    return policy


def parse_profile_log(path):
    with open(path, "r") as handle:
        text = handle.read()
    header = None
    rows = []
    stats = []
    resources = {}
    for line in text.splitlines():
        if line.startswith("PROFILE "):
            if header is not None:
                raise AssertionError("duplicate PROFILE header")
            header = {key: parse_number(value)
                      for key, value in parse_key_values(line, 1).items()}
        elif line.startswith("ROW "):
            fields = line.split()
            if len(fields) != 12:
                raise AssertionError("malformed ROW in {}".format(path))
            try:
                values = [complex(float(fields[index]),
                                  float(fields[index + 1]))
                          for index in range(4, 12, 2)]
                row = {
                    "sample": int(fields[1]),
                    "sector_index": int(fields[2]),
                    "configuration": int(fields[3]),
                    "values": values,
                    "line": line,
                }
            except ValueError:
                raise AssertionError("non-numeric ROW in {}".format(path))
            if any(not finite_complex(value) for value in values):
                raise AssertionError("nonfinite ROW value in {}".format(path))
            rows.append(row)
        elif line.startswith("STAT "):
            fields = line.split()
            if len(fields) < 3:
                raise AssertionError("malformed STAT in {}".format(path))
            stat = {key: parse_number(value)
                    for key, value in
                    parse_key_values_from_tokens(fields[2:]).items()}
            stat["sample"] = int(fields[1])
            stats.append(stat)
        elif line.startswith("RESOURCE "):
            fields = parse_key_values(line, 1)
            scope = fields.pop("scope", None)
            if scope not in ("rank_sum", "rank_max"):
                raise AssertionError("invalid RESOURCE scope in {}".format(path))
            resources[scope] = {key: parse_number(value)
                                for key, value in fields.items()}
    if header is None or not rows or set(resources) != {"rank_sum", "rank_max"}:
        raise AssertionError("incomplete profile log {}".format(path))
    if int(header.get("audit", 0)) != 1:
        raise AssertionError("P4-C evidence requires audit logs")
    if len(stats) not in (0, len(rows)):
        raise AssertionError("STAT count mismatch in {}".format(path))
    return {
        "path": os.path.abspath(path),
        "sha256": sha256_file(path),
        "header": header,
        "rows": rows,
        "stats": stats,
        "trace_sha256": sha256_text(
            "\n".join(row["line"] for row in rows) + "\n"),
        "resources": resources,
    }


def validate_profile_run(run, policy, allow_target_site=False):
    header = run["header"]
    required = {
        "schema", "fixture", "site_count", "qp_total", "sample_count",
        "sector_size", "rank_count", "cache_requested_bytes", "audit",
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
    cache_bytes = int(header["cache_requested_bytes"])
    scope = policy["scope"]
    allowed_sites = set(scope["site_counts"])
    if allow_target_site:
        allowed_sites.add(int(scope["target_site_count"]))
    if (site not in allowed_sites or qp not in scope["qp_total"] or
            ranks not in scope["mpi_rank_counts"] or
            cache_bytes not in scope["cache_capacity_bytes"]):
        raise AssertionError("profile header is outside the frozen grid")
    expected_sector = fixed_sector_size(site)
    sample_count = int(header["sample_count"])
    if (int(header["sector_size"]) != expected_sector or
            sample_count <= 0 or sample_count > expected_sector):
        raise AssertionError("profile sector/sample contract mismatch")
    observed_rows = [(row["sample"], row["sector_index"], row["configuration"])
                     for row in run["rows"]]
    if observed_rows != expected_profile_rows(site, sample_count):
        raise AssertionError("profile sample/configuration schedule mismatch")
    if run["stats"]:
        observed_stats = [int(stat["sample"]) for stat in run["stats"]]
        if observed_stats != list(range(sample_count)):
            raise AssertionError("profile STAT schedule mismatch")
    rank_sum = run["resources"]["rank_sum"]
    rank_max = run["resources"]["rank_max"]
    if (int(rank_max.get("roots", -1)) != sample_count or
            int(rank_sum.get("roots", -1)) != sample_count * ranks):
        raise AssertionError("profile root-count invariant mismatch")
    if int(rank_max.get("engine_heap_allocations", -1)) != 0:
        raise AssertionError("bounded engine hot path allocated on heap")
    for resource in (rank_sum, rank_max):
        for field in (
                "node_expansions", "raw_transitions", "terminal_calls",
                "row_peak", "cache_entries_peak", "plan_bytes",
                "model_bytes", "engine_workspace_bytes",
                "cache_requested_bytes", "cache_allocated_bytes",
                "collective_workspace_bytes", "amplitude_workspace_bytes",
                "rss_bytes", "total_seconds", "depth_seconds"):
            if field not in resource:
                raise AssertionError("RESOURCE field missing: {}".format(field))
        if len(resource["depth_seconds"]) != len(ORDERS):
            raise AssertionError("invalid depth timer vector")
        if any(float(value) < 0.0 or not math.isfinite(float(value))
               for value in resource["depth_seconds"]):
            raise AssertionError("invalid depth timer")
        if float(resource["total_seconds"]) <= 0.0:
            raise AssertionError("invalid total timer")
    if int(rank_max["cache_requested_bytes"]) != cache_bytes:
        raise AssertionError("cache request mismatch in rank_max RESOURCE")
    if int(rank_sum["cache_requested_bytes"]) != cache_bytes * ranks:
        raise AssertionError("cache request mismatch in rank_sum RESOURCE")


def validate_resource_run(run, policy):
    validate_profile_run(run, policy, allow_target_site=False)
    header = run["header"]
    site = int(header["site_count"])
    if int(header["sample_count"]) != expected_resource_sample_count(
            policy, site):
        raise AssertionError("resource sample count does not match policy")


def validate_exact_run(run, policy):
    validate_profile_run(run, policy, allow_target_site=False)
    header = run["header"]
    if int(header["rank_count"]) != 1:
        raise AssertionError("exact IID enumeration must use rank_count=1")
    if int(header["sample_count"]) != int(header["sector_size"]):
        raise AssertionError("exact IID enumeration must cover the full sector")


def allocated_capacity_bytes(resource):
    return (int(resource["plan_bytes"]) +
            int(resource["model_bytes"]) +
            int(resource["engine_workspace_bytes"]) +
            int(resource["collective_workspace_bytes"]) +
            int(resource["amplitude_workspace_bytes"]))


def summarize_resource_group(members):
    roots = [int(member["resources"]["rank_max"]["roots"])
             for member in members]
    if len(set(roots)) != 1:
        raise AssertionError("resource group roots mismatch")
    total_per_root = [
        float(member["resources"]["rank_max"]["total_seconds"]) /
        float(member["resources"]["rank_max"]["roots"])
        for member in members
    ]
    depth_per_root = []
    for member in members:
        rank_max = member["resources"]["rank_max"]
        roots_value = float(rank_max["roots"])
        depth_per_root.append([float(value) / roots_value
                               for value in rank_max["depth_seconds"]])
    return {
        "repeat_count": len(members),
        "sample_count": int(members[0]["header"]["sample_count"]),
        "sector_size": int(members[0]["header"]["sector_size"]),
        "trace_sha256": members[0]["trace_sha256"],
        "median_total_seconds_per_root": median(total_per_root),
        "max_total_seconds_per_root": max(total_per_root),
        "median_depth_seconds_per_root": [
            median(row[order] for row in depth_per_root)
            for order in ORDERS
        ],
        "max_allocated_capacity_bytes_per_rank": max(
            allocated_capacity_bytes(member["resources"]["rank_max"])
            for member in members),
        "max_peak_rss_bytes_per_rank": max(
            int(member["resources"]["rank_max"]["rss_bytes"])
            for member in members),
        "rank_max_counter_median": {
            key: median(member["resources"]["rank_max"][key]
                        for member in members)
            for key in (
                "node_expansions", "raw_transitions", "terminal_calls",
                "row_peak", "cache_entries_peak", "cache_insertions",
                "cache_evictions")
        },
        "rank_sum_counter_median": {
            key: median(member["resources"]["rank_sum"][key]
                        for member in members)
            for key in (
                "node_expansions", "raw_transitions", "terminal_calls",
                "cache_insertions", "cache_evictions")
        },
    }


def integrand_indices(family, i, j):
    if family == "S":
        return i, j
    if family == "K":
        return i, j + 1
    if family == "B":
        return i + 1, j + 1
    raise AssertionError("unknown integrand family")


def exact_iid_statistics(rows, policy):
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
    if not math.isfinite(target_sum) or target_sum <= 0.0:
        raise AssertionError("invalid target sum")
    target_mean = target_sum / len(target)
    sample_count = int(
        policy["p4_c_engine_feasibility"]["exact_iid_o5"]
        ["sample_count_proxy"])
    entry_budget = {}
    entry_values = {}
    matrices = {family: [[None] * 3 for _ in MATRIX_INDICES]
                for family in ("S", "K", "B")}
    for family in ("S", "K", "B"):
        for i in MATRIX_INDICES:
            for j in MATRIX_INDICES:
                left, right = integrand_indices(family, i, j)
                samples = [row[left].conjugate() * row[right]
                           for row in scaled]
                unsigned_scale = sum(abs(value) for value in samples) / \
                    target_sum
                budget = max(1.0e-12, 0.02 * unsigned_scale)
                exact_target = sum(samples) / target_sum
                cancellation_denominator = sum(abs(value)
                                               for value in samples)
                cancellation = (abs(sum(samples)) / cancellation_denominator
                                if cancellation_denominator != 0.0 else 1.0)
                name = "{}_{}{}".format(family, i, j)
                entry_values[name] = samples
                entry_budget[name] = {
                    "unsigned_scale": unsigned_scale,
                    "absolute_standard_error_budget": budget,
                    "exact_normalized_target": complex_pair(exact_target),
                    "cancellation_ratio": cancellation,
                }
                matrices[family][i][j] = complex_pair(sum(samples))
    guides = {}
    for rho in policy["guide"]["rho"]:
        eta = float(rho) * target_mean
        guide = [sum(abs(row[order]) ** 2 for order in ORDERS) + eta
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
        if not math.isfinite(mean_weight) or mean_weight <= 0.0:
            raise AssertionError("invalid self-normalized denominator")
        entries = {}
        all_pass = True
        max_ratio = 0.0
        for name, samples in entry_values.items():
            theta = sum(samples) / target_sum
            influence = [
                samples[index] / guide[index] -
                theta * target[index] / guide[index]
                for index in range(len(guide))
            ]
            variance = sum(
                probabilities[index] * abs(influence[index]) ** 2
                for index in range(len(guide))) / (abs(mean_weight) ** 2)
            standard_error = math.sqrt(max(variance, 0.0) / sample_count)
            budget = entry_budget[name]["absolute_standard_error_budget"]
            ratio = standard_error / budget if budget > 0.0 else float("inf")
            passed = standard_error <= budget
            all_pass = all_pass and passed
            max_ratio = max(max_ratio, ratio)
            entries[name] = {
                "theta": complex_pair(theta),
                "influence_complex_variance": variance,
                "absolute_standard_error_at_4096": standard_error,
                "absolute_standard_error_budget": budget,
                "se_to_budget_ratio": ratio,
                "budget_pass": passed,
                "cancellation_ratio":
                    entry_budget[name]["cancellation_ratio"],
            }
        guide_key = "r3_rho_{:.0e}".format(float(rho))
        guides[guide_key] = {
            "r": 3,
            "rho": float(rho),
            "eta": eta,
            "guide_min": min(guide),
            "guide_max": max(guide),
            "guide_sum": guide_sum,
            "strict_positive_support": True,
            "finite_weight": all(math.isfinite(value) for value in weights),
            "t_weight_tail": {
                "max": max(weights),
                "median": percentile(weights, 0.5),
                "p99": percentile(weights, 0.99),
            },
            "entries": entries,
            "all_entry_budgets_pass": all_pass,
            "maximum_required_entry_se_to_budget_ratio": max_ratio,
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


def materialize_raw_files(paths, output, subdirectory):
    destination_root = os.path.join(os.path.abspath(output), subdirectory)
    if not os.path.isdir(destination_root):
        os.makedirs(destination_root)
    result = {}
    used = set()
    for index, source in enumerate(paths):
        digest = sha256_file(source)
        base = os.path.basename(source)
        name = "{:04d}-{}-{}".format(index, digest[:12], base)
        if name in used:
            raise AssertionError("raw evidence name collision")
        used.add(name)
        destination = os.path.join(destination_root, name)
        if os.path.abspath(source) != os.path.abspath(destination):
            shutil.copyfile(source, destination)
        if sha256_file(destination) != digest:
            raise AssertionError("raw evidence copy checksum mismatch")
        result[os.path.abspath(source)] = os.path.join(subdirectory, name)
    return result


def common_metadata(args, policy_hash):
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
        "measurement_policy_sha256": policy_hash,
    }


def group_profile_runs(runs):
    groups = {}
    for run in runs:
        header = run["header"]
        key = (
            int(header["cache_requested_bytes"]),
            int(header["site_count"]),
            int(header["qp_total"]),
            int(header["rank_count"]),
        )
        groups.setdefault(key, []).append(run)
    return groups


def build_resource_statistics(args, policy, runs, allocation_runs, raw_files):
    scope = policy["scope"]
    expected_repeats = int(scope["repeat_count"])
    groups = group_profile_runs(runs)
    expected_keys = {
        (cache, site, qp, rank)
        for cache in scope["cache_capacity_bytes"]
        for site in scope["site_counts"]
        for qp in scope["qp_total"]
        for rank in scope["mpi_rank_counts"]
    }
    if set(groups) != expected_keys:
        raise AssertionError("resource grid mismatch: missing={} extra={}".format(
            sorted(expected_keys - set(groups)),
            sorted(set(groups) - expected_keys)))
    summaries = {}
    cross_rank_trace = {}
    for key, members in sorted(groups.items()):
        if len(members) != expected_repeats:
            raise AssertionError("repeat count mismatch for {}".format(key))
        traces = {member["trace_sha256"] for member in members}
        if len(traces) != 1:
            raise AssertionError("repeat trace mismatch for {}".format(key))
        cache, site, qp, rank = key
        summaries["cache{}_L{}_qp{}_rank{}".format(
            cache, site, qp, rank)] = summarize_resource_group(members)
    for cache in scope["cache_capacity_bytes"]:
        for site in scope["site_counts"]:
            for qp in scope["qp_total"]:
                traces = {
                    member["trace_sha256"]
                    for rank in scope["mpi_rank_counts"]
                    for member in groups[(cache, site, qp, rank)]
                }
                if len(traces) != 1:
                    raise AssertionError(
                        "rank-count trace mismatch for cache{} L{} qp{}".
                        format(cache, site, qp))
                cross_rank_trace["cache{}_L{}_qp{}".format(
                    cache, site, qp)] = {
                    "rank_counts": scope["mpi_rank_counts"],
                    "trace_sha256": next(iter(traces)),
                    "status": "PASS",
                }
    allocation_by_key = {}
    for run in allocation_runs:
        header = run["header"]
        key = (int(header["cache_requested_bytes"]), int(header["qp_total"]),
               int(header["rank_count"]))
        value = allocated_capacity_bytes(run["resources"]["rank_max"])
        allocation_by_key[key] = max(allocation_by_key.get(key, 0), value)
    extrapolation = {}
    budgets = policy["p4_c_engine_feasibility"]["target_budgets"]
    workload = policy["p4_c_engine_feasibility"]["representative_workload"]
    for cache in scope["cache_capacity_bytes"]:
        for qp in scope["qp_total"]:
            for rank in scope["mpi_rank_counts"]:
                runtime_points = []
                rss_points = []
                allocated_points = []
                for site in scope["site_counts"]:
                    summary = summaries["cache{}_L{}_qp{}_rank{}".format(
                        cache, site, qp, rank)]
                    runtime_points.append(
                        (site, summary["median_total_seconds_per_root"]))
                    rss_points.append(
                        (site, summary["max_peak_rss_bytes_per_rank"]))
                    allocated_points.append(
                        (site, summary[
                            "max_allocated_capacity_bytes_per_rank"]))
                proposal_l_target = log_linear_fit(
                    runtime_points, scope["target_site_count"])
                rss_l_target = log_linear_fit(
                    rss_points, scope["target_site_count"])
                allocation_key = (cache, qp, rank)
                allocated_target = allocation_by_key.get(allocation_key)
                allocated_source = "target_allocation_smoke"
                if allocated_target is None:
                    if not args.allow_smoke_policy:
                        raise AssertionError(
                            "missing target allocation smoke for cache{} qp{} "
                            "rank{}".format(cache, qp, rank))
                    allocated_target = log_linear_fit(
                        allocated_points, scope["target_site_count"]
                    )["prediction"]
                    allocated_source = "log_linear_fallback_missing_smoke"
                saved_increment = 0.0
                combined = (
                    float(workload["proposal_count"]) *
                    proposal_l_target["prediction"] +
                    float(workload["saved_sample_count"]) * saved_increment)
                pass_gate = (
                    proposal_l_target["prediction"] <=
                    float(budgets["proposal_seconds_per_configuration"]) and
                    saved_increment <=
                    float(budgets[
                        "saved_increment_seconds_per_configuration"]) and
                    combined <=
                    float(budgets[
                        "combined_workload_seconds_per_rank"]) and
                    float(allocated_target) <=
                    float(budgets[
                        "allocated_capacity_bytes_per_rank"]) and
                    rss_l_target["prediction"] <=
                    float(budgets["peak_rss_bytes_per_rank"])
                )
                extrapolation["cache{}_qp{}_rank{}".format(
                    cache, qp, rank)] = {
                    "cache_requested_bytes": cache,
                    "qp_total": qp,
                    "rank_count": rank,
                    "proposal_seconds_per_configuration_l{}".
                    format(scope["target_site_count"]): proposal_l_target,
                    "saved_increment_seconds_per_configuration_l{}".
                    format(scope["target_site_count"]): saved_increment,
                    "combined_workload_seconds_per_rank_l{}".
                    format(scope["target_site_count"]): combined,
                    "peak_rss_bytes_per_rank_l{}".
                    format(scope["target_site_count"]): rss_l_target,
                    "allocated_capacity_bytes_per_rank_l{}".
                    format(scope["target_site_count"]): allocated_target,
                    "allocated_capacity_source": allocated_source,
                    "resource_gate_pass": pass_gate,
                }
    return {
        "artifact": "p4c_resource_statistics",
        "counter_semantics": {
            "rank_sum": "sum of per-rank work counters",
            "rank_max": "maximum per-rank work/resource counter",
            "allocated_capacity_bytes_per_rank":
                "plan + model + engine workspace + collective + amplitude",
            "runtime_gate_timer": "rank-max total_seconds including reset",
        },
        "raw_logs": [{
            "file": raw_files[run["path"]],
            "sha256": run["sha256"],
            "header": run["header"],
            "trace_sha256": run["trace_sha256"],
            "rank_sum": run["resources"]["rank_sum"],
            "rank_max": run["resources"]["rank_max"],
        } for run in runs],
        "raw_allocation_logs": [{
            "file": raw_files[run["path"]],
            "sha256": run["sha256"],
            "header": run["header"],
            "trace_sha256": run["trace_sha256"],
            "rank_max": run["resources"]["rank_max"],
        } for run in allocation_runs],
        "cross_rank_trace_invariance": cross_rank_trace,
        "group_summary": summaries,
        "target_extrapolation": extrapolation,
    }


def build_exact_statistics(policy, exact_runs, raw_files):
    scope = policy["scope"]
    groups = {}
    for run in exact_runs:
        header = run["header"]
        key = (int(header["site_count"]), int(header["qp_total"]),
               int(header["cache_requested_bytes"]))
        if key in groups:
            raise AssertionError("duplicate exact IID log for {}".format(key))
        groups[key] = run
    expected = {
        (site, qp, cache)
        for site in scope["site_counts"]
        for qp in scope["qp_total"]
        for cache in scope["cache_capacity_bytes"]
    }
    if set(groups) != expected:
        raise AssertionError("exact IID grid mismatch: missing={} extra={}".
                             format(sorted(expected - set(groups)),
                                    sorted(set(groups) - expected)))
    cache_invariance = {}
    cases = {}
    candidate_summary = {}
    for site in scope["site_counts"]:
        for qp in scope["qp_total"]:
            traces = {
                groups[(site, qp, cache)]["trace_sha256"]
                for cache in scope["cache_capacity_bytes"]
            }
            if len(traces) != 1:
                raise AssertionError(
                    "cache changes exact values for L{} qp{}".format(site, qp))
            cache_invariance["L{}_qp{}".format(site, qp)] = {
                "cache_capacity_bytes": scope["cache_capacity_bytes"],
                "trace_sha256": next(iter(traces)),
                "status": "PASS",
            }
            base = groups[(site, qp, scope["cache_capacity_bytes"][0])]
            stats = exact_iid_statistics(base["rows"], policy)
            cases["L{}_qp{}".format(site, qp)] = stats
            for cache in scope["cache_capacity_bytes"]:
                for guide_key, guide in stats["guides"].items():
                    key = "cache{}_{}".format(cache, guide_key)
                    summary = candidate_summary.setdefault(key, {
                        "cache_requested_bytes": cache,
                        "r": 3,
                        "rho": guide["rho"],
                        "all_exact_iid_o5_pass": True,
                        "all_strict_positive_support": True,
                        "all_finite_weight": True,
                        "maximum_required_entry_se_to_budget_ratio": 0.0,
                    })
                    summary["all_exact_iid_o5_pass"] = (
                        summary["all_exact_iid_o5_pass"] and
                        guide["all_entry_budgets_pass"])
                    summary["all_strict_positive_support"] = (
                        summary["all_strict_positive_support"] and
                        guide["strict_positive_support"])
                    summary["all_finite_weight"] = (
                        summary["all_finite_weight"] and
                        guide["finite_weight"])
                    summary[
                        "maximum_required_entry_se_to_budget_ratio"] = max(
                            summary[
                                "maximum_required_entry_se_to_budget_ratio"],
                            guide[
                                "maximum_required_entry_se_to_budget_ratio"])
    return {
        "artifact": "p4c_exact_iid_statistics",
        "sample_count_for_budget":
            policy["p4_c_engine_feasibility"]["exact_iid_o5"]
            ["sample_count_proxy"],
        "raw_logs": [{
            "file": raw_files[run["path"]],
            "sha256": run["sha256"],
            "header": run["header"],
            "trace_sha256": run["trace_sha256"],
        } for run in exact_runs],
        "cache_invariance": cache_invariance,
        "cases": cases,
        "candidate_summary": candidate_summary,
    }


def select_candidates(policy, resource, exact, correctness_pass):
    promoted = []
    all_resource = {}
    for key, value in resource["target_extrapolation"].items():
        candidate_key = "cache{}".format(value["cache_requested_bytes"])
        summary = all_resource.setdefault(candidate_key, {
            "cache_requested_bytes": value["cache_requested_bytes"],
            "resource_gate_pass": True,
            "maximum_combined_runtime_seconds_per_rank": 0.0,
            "maximum_proposal_seconds_per_configuration": 0.0,
            "maximum_peak_rss_bytes_per_rank": 0.0,
            "maximum_allocated_capacity_bytes_per_rank": 0.0,
        })
        summary["resource_gate_pass"] = (
            summary["resource_gate_pass"] and value["resource_gate_pass"])
        combined_keys = [
            name for name in value
            if name.startswith("combined_workload_seconds_per_rank_l")
        ]
        proposal_keys = [
            name for name in value
            if name.startswith("proposal_seconds_per_configuration_l")
        ]
        rss_keys = [
            name for name in value
            if name.startswith("peak_rss_bytes_per_rank_l")
        ]
        alloc_keys = [
            name for name in value
            if name.startswith("allocated_capacity_bytes_per_rank_l")
        ]
        summary["maximum_combined_runtime_seconds_per_rank"] = max(
            summary["maximum_combined_runtime_seconds_per_rank"],
            max(float(value[name]) for name in combined_keys))
        summary["maximum_proposal_seconds_per_configuration"] = max(
            summary["maximum_proposal_seconds_per_configuration"],
            max(float(value[name]["prediction"]) for name in proposal_keys))
        summary["maximum_peak_rss_bytes_per_rank"] = max(
            summary["maximum_peak_rss_bytes_per_rank"],
            max(float(value[name]["prediction"]) for name in rss_keys))
        summary["maximum_allocated_capacity_bytes_per_rank"] = max(
            summary["maximum_allocated_capacity_bytes_per_rank"],
            max(float(value[name]) for name in alloc_keys))
    candidates = {}
    for name, exact_candidate in exact["candidate_summary"].items():
        cache = exact_candidate["cache_requested_bytes"]
        resource_key = "cache{}".format(cache)
        resource_candidate = all_resource.get(resource_key)
        if resource_candidate is None:
            raise AssertionError("missing resource candidate for cache {}".format(
                cache))
        gate_pass = (
            correctness_pass and resource_candidate["resource_gate_pass"] and
            exact_candidate["all_exact_iid_o5_pass"] and
            exact_candidate["all_strict_positive_support"] and
            exact_candidate["all_finite_weight"])
        merged = dict(exact_candidate)
        merged.update({
            "resource_gate": resource_candidate,
            "correctness_gate_pass": correctness_pass,
            "p4c_candidate_gate_pass": gate_pass,
        })
        candidates[name] = merged
        if gate_pass:
            promoted.append(name)
    promoted = sorted(
        promoted,
        key=lambda name: (
            candidates[name]["cache_requested_bytes"],
            candidates[name]["maximum_required_entry_se_to_budget_ratio"],
            candidates[name]["resource_gate"]
            ["maximum_combined_runtime_seconds_per_rank"],
            candidates[name]["rho"],
        ))
    return candidates, promoted


def write_checksums(output):
    checksum_path = os.path.join(output, "checksums.sha256")
    entries = []
    for root, directories, files in os.walk(output):
        directories.sort()
        for filename in sorted(files):
            path = os.path.join(root, filename)
            if os.path.abspath(path) == os.path.abspath(checksum_path):
                continue
            if os.path.islink(path):
                raise AssertionError("symlink is not allowed in evidence ledger")
            entries.append((os.path.relpath(path, output), sha256_file(path)))
    with open(checksum_path, "w") as handle:
        for relative, digest in sorted(entries):
            handle.write("{}  {}\n".format(digest, relative))
    return len(entries)


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


def analyze_command(args):
    policy = validate_policy_contract(
        read_json(args.policy), args.allow_smoke_policy)
    policy_hash = sha256_file(args.policy)
    if args.repository:
        validate_commit_object(args.repository, args.source_commit)
        validate_commit_object(args.repository, args.baseline_commit)
    if int(args.omp_threads) != int(policy["scope"]["omp_threads"]) or \
            int(args.blas_threads) != int(policy["scope"]["blas_threads"]):
        raise AssertionError("thread metadata does not match P4-C policy")
    resource_runs = [parse_profile_log(path) for path in args.resource_log]
    exact_runs = [parse_profile_log(path) for path in args.exact_log]
    allocation_runs = [
        parse_profile_log(path) for path in (args.allocation_log or [])
    ]
    if not args.allow_smoke_policy and not allocation_runs:
        raise AssertionError("official P4-C evidence requires allocation logs")
    for run in resource_runs:
        validate_resource_run(run, policy)
    for run in exact_runs:
        validate_exact_run(run, policy)
    for run in allocation_runs:
        validate_profile_run(run, policy, allow_target_site=True)
        if int(run["header"]["site_count"]) != \
                int(policy["scope"]["target_site_count"]):
            raise AssertionError("allocation smoke must use target site count")
        if int(run["header"]["sample_count"]) != 1:
            raise AssertionError("allocation smoke must use one sample")
    output = os.path.abspath(args.output)
    if not os.path.isdir(output):
        os.makedirs(output)
    all_raw_runs = resource_runs + exact_runs + allocation_runs
    raw_files = materialize_raw_files(
        [run["path"] for run in all_raw_runs], output, "raw_p4c_logs")
    resource = build_resource_statistics(
        args, policy, resource_runs, allocation_runs, raw_files)
    exact = build_exact_statistics(policy, exact_runs, raw_files)
    correctness_pass = args.correctness_gate == "PASS"
    candidates, promoted = select_candidates(
        policy, resource, exact, correctness_pass)
    metadata = common_metadata(args, policy_hash)
    resource.update(metadata)
    exact.update(metadata)
    decision = {
        "artifact": "p4c_candidate_decision",
        "measurement_policy": policy,
        "candidate_summary": candidates,
        "promoted_candidates": promoted,
        "p4c_decision": "GO" if promoted else "STOP",
        "selected_candidate": promoted[0] if promoted else None,
        "correctness_gate": {
            "status": args.correctness_gate,
            "note": args.correctness_note,
            "pass": correctness_pass,
        },
        "decision_rule":
            policy["p4_c_engine_feasibility"]["decision_rule"],
    }
    decision.update(metadata)
    write_json(os.path.join(output, "p4c_resource_statistics.json"), resource)
    write_json(os.path.join(output, "p4c_exact_iid_statistics.json"), exact)
    write_json(os.path.join(output, "p4c_candidate_decision.json"), decision)
    count = write_checksums(output)
    print("P4-C evidence complete: {} ({} ledger files)".format(
        decision["p4c_decision"], count))


def add_common_arguments(parser):
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--baseline-commit", required=True)
    parser.add_argument("--compiler", required=True)
    parser.add_argument("--mpi", required=True)
    parser.add_argument("--blas", required=True)
    parser.add_argument("--strict-fp", required=True)
    parser.add_argument("--omp-threads", type=int, default=1)
    parser.add_argument("--blas-threads", type=int, default=1)
    parser.add_argument("--repository")


def build_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    subparsers.required = True
    analyze = subparsers.add_parser("analyze")
    add_common_arguments(analyze)
    analyze.add_argument("--policy", required=True)
    analyze.add_argument("--output", required=True)
    analyze.add_argument("--resource-log", nargs="+", required=True)
    analyze.add_argument("--exact-log", nargs="+", required=True)
    analyze.add_argument("--allocation-log", nargs="*")
    analyze.add_argument("--correctness-gate", choices=("PASS", "FAIL"),
                         required=True)
    analyze.add_argument("--correctness-note", default="")
    analyze.add_argument("--allow-smoke-policy", action="store_true")
    analyze.set_defaults(func=analyze_command)
    return parser


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
