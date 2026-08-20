from __future__ import print_function

import argparse
import copy
import math
import os
import re
import shutil
import statistics

import generate_bounded_krylov_p4c_evidence as p4c
import generate_bounded_krylov_p4s9_official_statistics_evidence as official
import generate_bounded_krylov_p4s9_session_resource_evidence as session_resource
import generate_bounded_krylov_p4s9_target_session_evidence as target_session


def verify_policy_bindings(policy, args, allow_smoke_policy=False):
    bindings = {
        "session_resource_policy": p4c.sha256_file(
            args.session_resource_policy),
        "target_session_policy": p4c.sha256_file(
            args.target_session_policy),
        "stage_c_local_decision": p4c.sha256_file(
            args.stage_c_local_decision),
        "target_local_evidence": p4c.sha256_file(
            args.target_local_evidence),
        "target_local_decision": p4c.sha256_file(
            args.target_local_decision),
    }
    if not allow_smoke_policy:
        for name, digest in bindings.items():
            if policy["parent_sha256"][name] != digest:
                raise AssertionError(
                    "P4-S9 official resource binding mismatch {}".format(
                        name))
        stage_c = p4c.read_json(args.stage_c_local_decision)
        target = p4c.read_json(args.target_local_decision)
        if stage_c.get("local_resource_decision") != "PASS" or \
                target.get("local_target_stress_decision") != "PASS":
            raise AssertionError("P4-S9 local resource gates not PASS")
    return bindings


def validate_resource_scope(policy):
    resource = policy.get("official_resource", {})
    timing = resource.get("session_timing", {})
    target = resource.get("direct_target_stress", {})
    gate = resource.get("gate", {})
    if timing.get("driver") != "bounded_krylov_markov_timing_driver" or \
            timing.get("mode") != "long-direct-session" or \
            int(timing.get("schema_version", 0)) != 3 or \
            timing.get("rank_counts") != [1, 2, 4] or \
            timing.get("site_counts") != [4, 6, 8] or \
            timing.get("qp_total") != [1, 4] or \
            int(timing.get("sample_count", 0)) != 32768 or \
            int(timing.get("repeat_count", 0)) != 7 or \
            timing.get("rng_stream_order") != list(range(7)) or \
            int(timing.get("cache_bytes", 0)) != 268435456:
        raise AssertionError("P4-S9 official session timing scope mismatch")
    if target.get("driver") != "bounded_krylov_profile_driver" or \
            target.get("mode") != "session-profile" or \
            int(target.get("schema_version", 0)) != 2 or \
            target.get("rank_counts") != [1, 2, 4] or \
            int(target.get("site_count", 0)) != 16 or \
            target.get("qp_total") != [1, 4] or \
            int(target.get("sample_count", 0)) != 4096 or \
            int(target.get("repeat_count", 0)) != 3 or \
            target.get("root_schedule") != \
            "evenly_spaced_fixed_sector_indices_including_endpoints" or \
            int(target.get("cache_bytes", 0)) != 268435456:
        raise AssertionError("P4-S9 official direct target scope mismatch")
    if float(gate.get("seconds_per_configuration_l16_maximum", 0.0)) != \
            0.054931640625 or \
            float(gate.get(
                "projected_32768_root_seconds_maximum", 0.0)) != 1800.0 or \
            int(gate.get(
                "allocated_capacity_bytes_per_rank_maximum", 0)) != \
            4294967296 or \
            int(gate.get("peak_rss_bytes_per_rank_maximum", 0)) != \
            4294967296 or \
            float(gate.get(
                "required_speedup_over_p4c_cache256_maximum", 0.0)) != \
            2.9127591854032655 or \
            gate.get("rank_reducer") != "maximum" or \
            gate.get("repeat_reducer") != "median":
        raise AssertionError("P4-S9 official resource gate mismatch")
    return resource


def validate_timing_run(run, policy, rank_count, allow_smoke_policy=False):
    if rank_count not in [1, 2, 4]:
        raise AssertionError("P4-S9 unsupported official timing rank")
    compatible = copy.deepcopy(run)
    compatible["header"]["rank_count"] = 1
    validated = session_resource.validate_run(
        compatible, policy, allow_smoke_policy)
    validated["rank_count"] = rank_count
    return validated


def validate_target_run(run, policy, rank_count, allow_smoke_policy=False):
    if rank_count not in [1, 2, 4]:
        raise AssertionError("P4-S9 unsupported official target rank")
    compatible = copy.deepcopy(run)
    compatible["header"]["rank_count"] = 1
    validated = target_session.validate_run(
        compatible, policy, allow_smoke_policy)
    validated["rank_count"] = rank_count
    return validated


def copy_log(path, destination_root, prefix):
    if not os.path.isdir(destination_root):
        os.makedirs(destination_root)
    digest = p4c.sha256_file(path)
    name = "{}-{}-{}".format(prefix, digest[:12], os.path.basename(path))
    destination = os.path.join(destination_root, name)
    shutil.copyfile(path, destination)
    if p4c.sha256_file(destination) != digest:
        raise AssertionError("P4-S9 official resource log copy mismatch")
    return destination, digest


def analyze_command(args):
    for label, commit in (("source", args.source_commit),
                          ("baseline", args.baseline_commit)):
        if re.fullmatch(r"[0-9a-f]{40}", commit or "") is None:
            raise AssertionError(
                "P4-S9 {} commit provenance mismatch".format(label))
    policy = official.validate_policy(p4c.read_json(args.policy))
    resource_scope = validate_resource_scope(policy)
    bindings = verify_policy_bindings(
        policy, args, args.allow_smoke_policy)
    timing_policy = session_resource.validate_policy(
        p4c.read_json(args.session_resource_policy))
    target_policy = target_session.validate_policy(
        p4c.read_json(args.target_session_policy))
    if int(args.omp_threads) != 1 or int(args.blas_threads) != 1:
        raise AssertionError("P4-S9 official resource thread mismatch")
    if args.repository:
        p4c.validate_commit_object(args.repository, args.source_commit)
        p4c.validate_commit_object(args.repository, args.baseline_commit)

    timing_parsed = []
    for path in args.timing_log:
        run = session_resource.parse_log(path)
        rank = int(run["header"].get("rank_count", 0))
        validated = validate_timing_run(
            run, timing_policy, rank, args.allow_smoke_policy)
        timing_parsed.append((run, validated))
    if len({os.path.realpath(run["path"])
            for run, _ in timing_parsed}) != \
            len(timing_parsed):
        raise AssertionError("duplicate P4-S9 official timing log path")
    timing_observed = {
        (int(run["header"]["rank_count"]),
         int(run["header"]["site_count"]),
         int(run["header"]["qp_total"])) for run, _ in timing_parsed
    }
    if len(timing_observed) != len(timing_parsed):
        raise AssertionError("duplicate P4-S9 official timing run")
    if not args.allow_smoke_policy:
        timing_expected = {
            (rank, site, qp)
            for rank in resource_scope["session_timing"]["rank_counts"]
            for site in resource_scope["session_timing"]["site_counts"]
            for qp in resource_scope["session_timing"]["qp_total"]
        }
        if timing_observed != timing_expected:
            raise AssertionError("P4-S9 official timing census mismatch")

    target_parsed = []
    for path in args.profile_log:
        run = target_session.parse_log(path)
        rank = int(run["header"].get("rank_count", 0))
        validated = validate_target_run(
            run, target_policy, rank, args.allow_smoke_policy)
        target_parsed.append((run, validated))
    if len({os.path.realpath(run["path"])
            for run, _ in target_parsed}) != \
            len(target_parsed):
        raise AssertionError("duplicate P4-S9 official target log path")
    if not args.allow_smoke_policy:
        for rank in resource_scope["direct_target_stress"]["rank_counts"]:
            for qp in resource_scope["direct_target_stress"]["qp_total"]:
                count = sum(
                    int(run["header"]["rank_count"]) == rank and
                    int(run["header"]["qp_total"]) == qp
                    for run, _ in target_parsed)
                if count != int(resource_scope[
                        "direct_target_stress"]["repeat_count"]):
                    raise AssertionError(
                        "P4-S9 official target repeat census mismatch")

    output = os.path.abspath(args.output)
    frozen = os.path.join(output, "frozen_inputs")
    producer = os.path.join(output, "producer")
    raw_timing = os.path.join(output, "raw_official_timing")
    raw_target = os.path.join(output, "raw_official_target")
    for directory in (output, frozen, producer):
        if not os.path.isdir(directory):
            os.makedirs(directory)
    for path in (args.session_resource_policy, args.target_session_policy,
                 args.stage_c_local_decision, args.target_local_evidence,
                 args.target_local_decision):
        shutil.copyfile(path, os.path.join(frozen, os.path.basename(path)))
    shutil.copyfile(args.policy,
                    os.path.join(producer, os.path.basename(args.policy)))
    for path in (args.producer_timing_binary,
                 args.producer_profile_binary):
        shutil.copyfile(path, os.path.join(producer, os.path.basename(path)))

    timing_cases = []
    for run, validated in timing_parsed:
        rank = int(run["header"]["rank_count"])
        site = int(run["header"]["site_count"])
        qp = int(run["header"]["qp_total"])
        destination, digest = copy_log(
            run["path"], raw_timing,
            "R{}-L{}-QP{}".format(rank, site, qp))
        records = validated["records"]
        math_trace_text = "\n".join(
            "{}:{}:{}:{}".format(
                int(item["repeat"]["rng_stream"]),
                int(item["correctness"]["final_sector_index"]),
                int(item["correctness"]["final_configuration"]),
                int(item["correctness"]["statistical_trace_hash"]))
            for item in records) + "\n"
        timing_cases.append({
            "rank_count": rank,
            "site_count": site,
            "qp_total": qp,
            "sample_count": validated["sample_count"],
            "median_seconds_per_configuration": statistics.median(
                float(item["resource"]["seconds_per_configuration"])
                for item in records),
            "maximum_allocated_capacity_bytes_per_rank": max(
                int(item["resource"]["allocated_capacity_bytes"])
                for item in records),
            "maximum_peak_rss_bytes_per_rank": max(
                int(item["resource"]["peak_rss_bytes"])
                for item in records),
            "plan_hash": validated["plan_hash"],
            "amplitude_generation_hash":
                validated["amplitude_generation_hash"],
            "math_trace_sha256": p4c.sha256_text(math_trace_text),
            "raw_log": os.path.relpath(destination, output),
            "raw_log_sha256": digest,
            "raw_log_line_count": run["line_count"],
        })

    gate = resource_scope["gate"]
    timing_results = []
    timing_pass = True
    for rank in sorted({item["rank_count"] for item in timing_cases}):
        for qp in sorted({item["qp_total"] for item in timing_cases
                          if item["rank_count"] == rank}):
            members = [item for item in timing_cases
                       if item["rank_count"] == rank and
                       item["qp_total"] == qp]
            if len(members) != 3:
                if args.allow_smoke_policy:
                    continue
                raise AssertionError("incomplete official timing fit")
            runtime_fit = p4c.log_linear_fit(
                [(item["site_count"],
                  item["median_seconds_per_configuration"])
                 for item in members], 16)
            allocation_fit = p4c.log_linear_fit(
                [(item["site_count"],
                  item["maximum_allocated_capacity_bytes_per_rank"])
                 for item in members], 16)
            rss_fit = p4c.log_linear_fit(
                [(item["site_count"],
                  item["maximum_peak_rss_bytes_per_rank"])
                 for item in members], 16)
            runtime = runtime_fit["prediction"]
            projected = runtime * 32768.0
            allocation = max(
                allocation_fit["prediction"],
                max(item["maximum_allocated_capacity_bytes_per_rank"]
                    for item in members))
            rss = max(rss_fit["prediction"], max(
                item["maximum_peak_rss_bytes_per_rank"] for item in members))
            speedup = 0.16000264079973991 / runtime
            passed = (
                runtime <= float(gate[
                    "seconds_per_configuration_l16_maximum"]) and
                projected <= float(gate[
                    "projected_32768_root_seconds_maximum"]) and
                allocation <= int(gate[
                    "allocated_capacity_bytes_per_rank_maximum"]) and
                rss <= int(gate["peak_rss_bytes_per_rank_maximum"]) and
                speedup >= float(gate[
                    "required_speedup_over_p4c_cache256_maximum"]))
            timing_pass = timing_pass and passed
            timing_results.append({
                "rank_count": rank, "qp_total": qp,
                "seconds_per_configuration_l16": runtime_fit,
                "projected_32768_root_seconds": projected,
                "allocated_capacity_bytes_per_rank_l16": allocation,
                "peak_rss_bytes_per_rank_l16": rss,
                "speedup_over_p4c_cache256_maximum": speedup,
                "resource_pass": passed,
            })

    timing_trace_results = []
    timing_trace_pass = True
    for site in resource_scope["session_timing"]["site_counts"]:
        for qp in resource_scope["session_timing"]["qp_total"]:
            members = [item for item in timing_cases
                       if item["site_count"] == site and
                       item["qp_total"] == qp]
            if not members:
                if args.allow_smoke_policy:
                    continue
                raise AssertionError("missing official timing trace group")
            signatures = {
                (item["plan_hash"], item["amplitude_generation_hash"],
                 item["math_trace_sha256"])
                for item in members
            }
            ranks = sorted(item["rank_count"] for item in members)
            passed = len(signatures) == 1
            timing_trace_pass = timing_trace_pass and passed
            timing_trace_results.append({
                "site_count": site,
                "qp_total": qp,
                "rank_counts": ranks,
                "math_trace_sha256": next(iter(signatures))[2]
                    if passed else None,
                "cross_rank_trace_pass": passed,
            })
    timing_pass = timing_pass and timing_trace_pass

    target_repeats = []
    for repeat_index, (run, validated) in enumerate(target_parsed):
        rank = int(run["header"]["rank_count"])
        qp = int(run["header"]["qp_total"])
        destination, digest = copy_log(
            run["path"], raw_target,
            "R{}-QP{}-N{}".format(rank, qp, repeat_index))
        item = dict(validated)
        item.update({
            "raw_log": os.path.relpath(destination, output),
            "raw_log_sha256": digest,
            "raw_log_line_count": run["line_count"],
        })
        target_repeats.append(item)
    target_results = []
    target_pass = True
    groups = sorted({(item["rank_count"], item["qp_total"])
                     for item in target_repeats})
    for rank, qp in groups:
        members = [item for item in target_repeats
                   if item["rank_count"] == rank and
                   item["qp_total"] == qp]
        traces = {item["row_trace_sha256"] for item in members}
        generations = {item["amplitude_generation_hash"]
                       for item in members}
        median_time = statistics.median(
            item["total_seconds_per_root"] for item in members)
        projected = median_time * 32768.0
        allocation = max(item["allocated_capacity_bytes_per_rank"]
                         for item in members)
        rss = max(item["peak_rss_bytes_per_rank"] for item in members)
        speedup = 0.16000264079973991 / median_time
        passed = (
            len(traces) == 1 and len(generations) == 1 and
            median_time <= float(gate[
                "seconds_per_configuration_l16_maximum"]) and
            projected <= float(gate[
                "projected_32768_root_seconds_maximum"]) and
            allocation <= int(gate[
                "allocated_capacity_bytes_per_rank_maximum"]) and
            rss <= int(gate["peak_rss_bytes_per_rank_maximum"]) and
            speedup >= float(gate[
                "required_speedup_over_p4c_cache256_maximum"]))
        target_pass = target_pass and passed
        target_results.append({
            "rank_count": rank, "qp_total": qp,
            "repeat_total_seconds_per_root": [
                item["total_seconds_per_root"] for item in members],
            "median_total_seconds_per_root": median_time,
            "projected_32768_root_seconds": projected,
            "maximum_allocated_capacity_bytes_per_rank": allocation,
            "maximum_peak_rss_bytes_per_rank": rss,
            "speedup_over_p4c_cache256_maximum": speedup,
            "row_trace_sha256": next(iter(traces))
                if len(traces) == 1 else None,
            "amplitude_generation_hash": next(iter(generations))
                if len(generations) == 1 else None,
            "trace_invariance_pass":
                len(traces) == 1 and len(generations) == 1,
            "resource_pass": passed,
        })

    target_cross_rank_results = []
    target_cross_rank_pass = True
    for qp in resource_scope["direct_target_stress"]["qp_total"]:
        members = [item for item in target_results
                   if item["qp_total"] == qp]
        if not members:
            if args.allow_smoke_policy:
                continue
            raise AssertionError("missing official target trace group")
        signatures = {
            (item["row_trace_sha256"],
             item["amplitude_generation_hash"])
            for item in members
        }
        ranks = sorted(item["rank_count"] for item in members)
        passed = len(signatures) == 1
        target_cross_rank_pass = target_cross_rank_pass and passed
        target_cross_rank_results.append({
            "qp_total": qp,
            "rank_counts": ranks,
            "row_trace_sha256": next(iter(signatures))[0]
                if passed else None,
            "amplitude_generation_hash": next(iter(signatures))[1]
                if passed else None,
            "cross_rank_trace_pass": passed,
        })
    target_pass = target_pass and target_cross_rank_pass
    all_pass = timing_pass and target_pass
    generated_at = p4c.generated_at_jst()
    metadata = {
        "schema_version": 1,
        "generated_at_jst": generated_at,
        "source_commit": args.source_commit,
        "baseline_develop_commit": args.baseline_commit,
        "compiler": args.compiler,
        "mpi": args.mpi,
        "blas": args.blas,
        "strict_fp": args.strict_fp,
        "omp_threads": args.omp_threads,
        "blas_threads": args.blas_threads,
        "official_execution_policy_sha256": p4c.sha256_file(args.policy),
        "parent_sha256": bindings,
        "producer_timing_binary_sha256": p4c.sha256_file(
            args.producer_timing_binary),
        "producer_profile_binary_sha256": p4c.sha256_file(
            args.producer_profile_binary),
        "candidate_id": policy["candidate_id"],
    }
    evidence = dict(metadata)
    evidence.update({
        "artifact": "p4s9_official_resource_evidence",
        "session_timing_cases": timing_cases,
        "session_timing_results": timing_results,
        "session_timing_cross_rank_trace": timing_trace_results,
        "direct_target_repeats": target_repeats,
        "direct_target_results": target_results,
        "direct_target_cross_rank_trace": target_cross_rank_results,
        "official_session_timing_cross_rank_trace_pass":
            timing_trace_pass,
        "official_direct_target_cross_rank_trace_pass":
            target_cross_rank_pass,
        "official_session_timing_pass": timing_pass,
        "official_direct_target_pass": target_pass,
        "official_resource_pass": all_pass,
    })
    decision = dict(metadata)
    decision.update({
        "artifact": "p4s9_official_resource_decision",
        "official_resource_decision": "PASS" if all_pass else "STOP",
        "p4s9_combined_decision":
            "PENDING_OFFICIAL_STATISTICS_AND_MPI_PARITY"
            if all_pass else "STOP",
        "p4f_authorized": False,
        "failure_action": "continue_to_combined_closure"
            if all_pass else "preserve_evidence_and_stop_before_p4f",
    })
    p4c.write_json(os.path.join(
        output, "p4s9_official_resource_evidence.json"), evidence)
    p4c.write_json(os.path.join(
        output, "p4s9_official_resource_decision.json"), decision)
    count = p4c.write_checksums(output)
    print("P4-S9 official resource: {}; {} ledger files".format(
        decision["official_resource_decision"], count))


def build_parser():
    parser = argparse.ArgumentParser()
    p4c.add_common_arguments(parser)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--session-resource-policy", required=True)
    parser.add_argument("--target-session-policy", required=True)
    parser.add_argument("--stage-c-local-decision", required=True)
    parser.add_argument("--target-local-evidence", required=True)
    parser.add_argument("--target-local-decision", required=True)
    parser.add_argument("--producer-timing-binary", required=True)
    parser.add_argument("--producer-profile-binary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--timing-log", nargs="+", required=True)
    parser.add_argument("--profile-log", nargs="+", required=True)
    parser.add_argument("--allow-smoke-policy", action="store_true")
    return parser


def main(argv=None):
    analyze_command(build_parser().parse_args(argv))


if __name__ == "__main__":
    main()
