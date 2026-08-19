from __future__ import print_function

import argparse
import math
import os
import shutil
import statistics

import generate_bounded_krylov_p4c_evidence as p4c


POLICY_SCHEMA_VERSION = 1
DRIVER_SCHEMA_VERSION = 2
EXPECTED_POLICY_ID = "power-lanczos-zero-support-p4-r0-v11"
EXPECTED_FIXTURE = "p4s9_target_session_spread_roots"


def fields(line, skip=1):
    return {
        key: p4c.parse_number(value)
        for key, value in p4c.parse_key_values(line, skip).items()
    }


def validate_policy(policy):
    if int(policy.get("schema_version", 0)) != POLICY_SCHEMA_VERSION or \
            policy.get("policy_id") != EXPECTED_POLICY_ID:
        raise AssertionError("P4-S9 target policy identity mismatch")
    driver = policy.get("driver", {})
    if driver.get("binary") != "bounded_krylov_profile_driver" or \
            driver.get("mode") != "session-profile" or \
            int(driver.get("schema_version", 0)) != DRIVER_SCHEMA_VERSION or \
            driver.get("fixture") != EXPECTED_FIXTURE or \
            int(driver.get("site_count", 0)) != 16 or \
            int(driver.get("sector_size", 0)) != 165636900 or \
            driver.get("qp_total_order") != [1, 4] or \
            int(driver.get("sample_count", 0)) != 4096 or \
            int(driver.get("repeat_count", 0)) != 3 or \
            driver.get("root_schedule") != \
            "evenly_spaced_fixed_sector_indices_including_endpoints" or \
            int(driver.get("cache_bytes", 0)) != 268435456 or \
            int(driver.get("audit", 0)) != 1 or \
            int(driver.get("omp_threads", 0)) != 1 or \
            int(driver.get("blas_threads", 0)) != 1 or \
            int(driver.get("rank_count", 0)) != 1:
        raise AssertionError("P4-S9 target driver mismatch")
    gate = policy.get("resource_gate", {})
    if gate.get("rank_reducer") != "maximum" or \
            gate.get("repeat_reducer") != "median" or \
            float(gate.get("total_seconds_per_root_maximum", 0.0)) != \
            0.054931640625 or \
            float(gate.get("projected_32768_root_seconds_maximum", 0.0)) != \
            1800.0 or \
            int(gate.get("allocated_capacity_bytes_per_rank_maximum", 0)) != \
            4294967296 or \
            int(gate.get("peak_rss_bytes_per_rank_maximum", 0)) != \
            4294967296 or \
            float(gate.get(
                "p4c_cache256_maximum_seconds_per_configuration", 0.0)) != \
            0.16000264079973991 or \
            float(gate.get("required_speedup_minimum", 0.0)) != \
            2.9127591854032655 or \
            int(gate.get("engine_heap_allocations", -1)) != 0 or \
            int(gate.get("cache_reset_count", 0)) != 1:
        raise AssertionError("P4-S9 target resource gate mismatch")
    for name in ("parent_resource_policy_sha256", "stage_c_evidence_sha256",
                 "stage_c_decision_sha256"):
        if len(str(policy.get(name, ""))) != 64:
            raise AssertionError("P4-S9 target frozen hash missing")
    return policy


def parse_log(path):
    header = None
    rows = []
    stats = []
    resources = {}
    with open(path, "r") as handle:
        text = handle.read()
    for line in text.splitlines():
        if line.startswith("SESSION_PROFILE "):
            if header is not None:
                raise AssertionError("duplicate target session header")
            header = fields(line)
        elif line.startswith("ROW "):
            tokens = line.split()
            if len(tokens) != 12:
                raise AssertionError("target ROW shape mismatch")
            values = [float(value) for value in tokens[4:]]
            if any(not math.isfinite(value) for value in values):
                raise AssertionError("target ROW nonfinite")
            rows.append({
                "sample": int(tokens[1]),
                "sector_index": int(tokens[2]),
                "configuration": int(tokens[3]),
                "values": values,
                "line": line,
            })
        elif line.startswith("STAT "):
            item = fields(line, 2)
            item["sample"] = int(line.split()[1])
            stats.append(item)
        elif line.startswith("RESOURCE "):
            item = fields(line)
            scope = item.get("scope")
            if scope in resources:
                raise AssertionError("duplicate target resource scope")
            resources[scope] = item
    if header is None or not rows or len(stats) != len(rows) or \
            set(resources) != {"rank_sum", "rank_max"}:
        raise AssertionError("incomplete target session log {}".format(path))
    return {
        "path": os.path.abspath(path), "sha256": p4c.sha256_file(path),
        "line_count": len(text.splitlines()), "header": header,
        "rows": rows, "stats": stats, "resources": resources,
        "row_trace_sha256": p4c.sha256_text(
            "\n".join(item["line"] for item in rows) + "\n"),
    }


def allocated_capacity(resource):
    return (int(resource["plan_bytes"]) + int(resource["model_bytes"]) +
            int(resource["engine_workspace_bytes"]) +
            int(resource["collective_workspace_bytes"]) +
            int(resource["amplitude_workspace_bytes"]))


def validate_run(run, policy, allow_smoke_policy=False):
    header = run["header"]
    driver = policy["driver"]
    sample_count = int(header["sample_count"])
    cache_bytes = int(header["cache_requested_bytes"])
    qp = int(header["qp_total"])
    if int(header.get("schema", 0)) != DRIVER_SCHEMA_VERSION or \
            header.get("fixture") != EXPECTED_FIXTURE or \
            int(header.get("site_count", 0)) != 16 or \
            qp not in driver["qp_total_order"] or \
            int(header.get("sector_size", 0)) != 165636900 or \
            int(header.get("rank_count", 0)) != 1 or \
            int(header.get("audit", 0)) != 1 or \
            (not allow_smoke_policy and
             sample_count != int(driver["sample_count"])) or \
            (not allow_smoke_policy and
             cache_bytes != int(driver["cache_bytes"])):
        raise AssertionError("P4-S9 target header mismatch")
    expected_rows = p4c.expected_profile_rows(16, sample_count)
    observed_rows = [(item["sample"], item["sector_index"],
                      item["configuration"]) for item in run["rows"]]
    if observed_rows != expected_rows:
        raise AssertionError("P4-S9 target root schedule mismatch")
    generation = None
    previous_resident_end = 0
    for index, stat in enumerate(run["stats"]):
        if int(stat["sample"]) != index or \
                int(stat["session_active"]) != 1 or \
                int(stat["session_root_evaluation"]) != index + 1 or \
                int(stat["cache_reset_performed"]) != (1 if index == 0 else 0) or \
                int(stat["cache_resident_start"]) != previous_resident_end or \
                int(stat["cache_resident_end"]) < \
                int(stat["cache_resident_start"]) or \
                int(stat["workspace_bytes"]) <= 0 or \
                int(stat["cache_allocated_bytes"]) > cache_bytes:
            raise AssertionError("P4-S9 target session lifecycle mismatch")
        current_generation = int(stat["session_generation"])
        if current_generation == 0 or \
                (generation is not None and current_generation != generation):
            raise AssertionError("P4-S9 target generation mismatch")
        generation = current_generation
        previous_resident_end = int(stat["cache_resident_end"])
    resource = run["resources"]["rank_max"]
    if int(resource["roots"]) != sample_count or \
            int(resource["cache_requested_bytes"]) != cache_bytes or \
            int(resource["engine_heap_allocations"]) != 0 or \
            int(resource["engine_workspace_bytes"]) <= 0 or \
            int(resource["cache_allocated_bytes"]) > cache_bytes or \
            int(resource["rss_bytes"]) <= 0 or \
            float(resource["total_seconds"]) <= 0.0 or \
            not math.isfinite(float(resource["total_seconds"])):
        raise AssertionError("P4-S9 target resource mismatch")
    total_per_root = float(resource["total_seconds"]) / sample_count
    return {
        "qp_total": qp, "sample_count": sample_count,
        "cache_requested_bytes": cache_bytes,
        "row_trace_sha256": run["row_trace_sha256"],
        "amplitude_generation_hash": generation,
        "total_seconds": float(resource["total_seconds"]),
        "total_seconds_per_root": total_per_root,
        "allocated_capacity_bytes_per_rank": allocated_capacity(resource),
        "peak_rss_bytes_per_rank": int(resource["rss_bytes"]),
        "cache_hits": resource["cache_hits"],
        "cache_misses": resource["cache_misses"],
        "cache_insertions": int(resource["cache_insertions"]),
        "cache_evictions": int(resource["cache_evictions"]),
        "cache_entries_peak": int(resource["cache_entries_peak"]),
        "terminal_amplitude_calls": int(resource["terminal_calls"]),
    }


def materialize_logs(parsed, output):
    root = os.path.join(output, "raw_target_session")
    if not os.path.isdir(root):
        os.makedirs(root)
    mapped = {}
    for repeat_index, run in enumerate(parsed):
        name = "QP{}-R{}-{}-{}".format(
            int(run["header"]["qp_total"]), repeat_index,
            run["sha256"][:12], os.path.basename(run["path"]))
        destination = os.path.join(root, name)
        shutil.copyfile(run["path"], destination)
        if p4c.sha256_file(destination) != run["sha256"]:
            raise AssertionError("P4-S9 target log copy mismatch")
        mapped[run["path"]] = os.path.relpath(destination, output)
    return mapped


def analyze_command(args):
    policy = validate_policy(p4c.read_json(args.policy))
    bindings = {
        "parent_resource_policy_sha256": p4c.sha256_file(
            args.parent_resource_policy),
        "stage_c_evidence_sha256": p4c.sha256_file(args.stage_c_evidence),
        "stage_c_decision_sha256": p4c.sha256_file(args.stage_c_decision),
    }
    if not args.allow_smoke_policy:
        for name, digest in bindings.items():
            if policy[name] != digest:
                raise AssertionError("P4-S9 target input mismatch {}".format(
                    name))
    parsed = [parse_log(path) for path in args.profile_log]
    expected_count = (len(policy["driver"]["qp_total_order"]) *
                      int(policy["driver"]["repeat_count"]))
    if len(parsed) != expected_count:
        raise AssertionError("P4-S9 target repeat census mismatch")
    validated = [validate_run(run, policy, args.allow_smoke_policy)
                 for run in parsed]
    for qp in policy["driver"]["qp_total_order"]:
        if sum(item["qp_total"] == qp for item in validated) != \
                int(policy["driver"]["repeat_count"]):
            raise AssertionError("P4-S9 target QP repeat census mismatch")
    output = os.path.abspath(args.output)
    if not os.path.isdir(output):
        os.makedirs(output)
    raw_paths = materialize_logs(parsed, output)
    repeats = []
    for run, item in zip(parsed, validated):
        item["raw_log"] = raw_paths[run["path"]]
        item["raw_log_sha256"] = run["sha256"]
        item["raw_log_line_count"] = run["line_count"]
        repeats.append(item)
    gate = policy["resource_gate"]
    qp_results = []
    all_pass = True
    for qp in policy["driver"]["qp_total_order"]:
        members = [item for item in repeats if item["qp_total"] == qp]
        traces = {item["row_trace_sha256"] for item in members}
        times = [item["total_seconds_per_root"] for item in members]
        median_time = statistics.median(times)
        projected = median_time * 32768.0
        maximum_allocated = max(
            item["allocated_capacity_bytes_per_rank"] for item in members)
        maximum_rss = max(item["peak_rss_bytes_per_rank"] for item in members)
        speedup = float(gate[
            "p4c_cache256_maximum_seconds_per_configuration"]) / median_time
        passed = (len(traces) == 1 and
                  median_time <= float(gate[
                      "total_seconds_per_root_maximum"]) and
                  projected <= float(gate[
                      "projected_32768_root_seconds_maximum"]) and
                  maximum_allocated <= int(gate[
                      "allocated_capacity_bytes_per_rank_maximum"]) and
                  maximum_rss <= int(gate[
                      "peak_rss_bytes_per_rank_maximum"]) and
                  speedup >= float(gate["required_speedup_minimum"]))
        all_pass = all_pass and passed
        qp_results.append({
            "qp_total": qp,
            "repeat_total_seconds_per_root": times,
            "median_total_seconds_per_root": median_time,
            "projected_32768_root_seconds": projected,
            "maximum_allocated_capacity_bytes_per_rank": maximum_allocated,
            "maximum_peak_rss_bytes_per_rank": maximum_rss,
            "speedup_over_p4c_cache256_maximum": speedup,
            "row_trace_sha256": next(iter(traces)) if len(traces) == 1 else None,
            "trace_invariance_pass": len(traces) == 1,
            "resource_pass": passed,
        })
    generated_at = p4c.generated_at_jst()
    evidence = {
        "schema_version": POLICY_SCHEMA_VERSION,
        "artifact": "p4s9_local_target_session_stress_evidence",
        "generated_at_jst": generated_at,
        "source_commit": args.source_commit,
        "producer_binary": os.path.basename(args.producer_binary),
        "producer_binary_sha256": p4c.sha256_file(args.producer_binary),
        "target_policy_sha256": p4c.sha256_file(args.policy),
        "frozen_input_sha256": bindings,
        "candidate_id": policy["candidate_id"],
        "local_target_stress_pass": all_pass,
        "repeats": repeats,
        "qp_results": qp_results,
    }
    decision = {
        "schema_version": POLICY_SCHEMA_VERSION,
        "artifact": "p4s9_frozen_local_target_session_decision",
        "generated_at_jst": generated_at,
        "source_commit": args.source_commit,
        "target_policy_sha256": evidence["target_policy_sha256"],
        "candidate_id": policy["candidate_id"],
        "local_target_stress_decision": "PASS" if all_pass else "STOP",
        "remote_package_preparation_authorized": all_pass,
        "remote_transfer_or_submission_authorized": False,
        "failure_action": (
            "prepare_package_and_request_explicit_user_approval"
            if all_pass else "formal_local_stop_no_remote_package"),
    }
    p4c.write_json(os.path.join(
        output, "p4s9_local_target_session_stress_evidence.json"), evidence)
    p4c.write_json(os.path.join(
        output, "p4s9_local_target_session_decision.json"), decision)
    count = p4c.write_checksums(output)
    print("P4-S9 local target session: {}; {} ledger files".format(
        decision["local_target_stress_decision"], count))


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--policy", required=True)
    parser.add_argument("--parent-resource-policy", required=True)
    parser.add_argument("--stage-c-evidence", required=True)
    parser.add_argument("--stage-c-decision", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--profile-log", nargs="+", required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--producer-binary", required=True)
    parser.add_argument("--allow-smoke-policy", action="store_true")
    return parser


def main(argv=None):
    analyze_command(build_parser().parse_args(argv))


if __name__ == "__main__":
    main()
