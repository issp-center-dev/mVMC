from __future__ import print_function

import argparse
import copy
import math
import os
import re
import shutil

import generate_bounded_krylov_p4c_evidence as p4c
import generate_bounded_krylov_p4s_evidence as p4s


POLICY_SCHEMA_VERSION = 1
EXPECTED_POLICY_ID = "power-lanczos-zero-support-p4-r0-v12"
DRIVER_SCHEMA_VERSION = 3
EXPECTED_FIXTURE = "p4s9_long_direct_session_official"


def validate_policy(policy):
    if int(policy.get("schema_version", 0)) != POLICY_SCHEMA_VERSION or \
            policy.get("policy_id") != EXPECTED_POLICY_ID:
        raise AssertionError("P4-S9 official policy identity mismatch")
    statistics = policy.get("official_statistics", {})
    if statistics.get("driver") != "bounded_krylov_markov_driver" or \
            int(statistics.get("driver_schema_version", 0)) != \
            DRIVER_SCHEMA_VERSION or \
            statistics.get("fixture") != EXPECTED_FIXTURE or \
            statistics.get("mode_argument") != "session" or \
            statistics.get("rank_counts") != [1] or \
            statistics.get("site_counts") != [4, 6, 8] or \
            statistics.get("qp_total") != [1, 4] or \
            statistics.get("physical_case_order") != [
                [4, 1], [4, 4], [6, 1], [6, 4], [8, 1], [8, 4]] or \
            int(statistics.get("cache_bytes", 0)) != 268435456 or \
            float(statistics.get("rho", 0.0)) != 0.01 or \
            int(statistics.get("global_numerator", -1)) != 0 or \
            int(statistics.get("global_denominator", 0)) != 1 or \
            int(statistics.get("rng_stream", -1)) != 0 or \
            int(statistics.get("sample_count", 0)) != 32768 or \
            int(statistics.get("warmup_proposals", -1)) != 0 or \
            statistics.get("seed_hex_order") != [
                "0x5034533952305631", "0x5034533952305632",
                "0x5034533952305633", "0x5034533952305634"] or \
            statistics.get("promotion") != \
            "all_four_seeds_all_six_cases_must_GO":
        raise AssertionError("P4-S9 official statistics scope mismatch")
    official = statistics.get("official_partition", {})
    diagnostic = statistics.get("diagnostic_partition", {})
    session = statistics.get("session_audit", {})
    if official != {"block_count": 16, "block_length": 2048} or \
            diagnostic != {"block_count": 32, "block_length": 1024} or \
            float(statistics.get(
                "maximum_conservative_se_budget_ratio", 0.0)) != 1.0 or \
            float(statistics.get("maximum_tau_int", 0.0)) != 16.0 or \
            float(statistics.get(
                "minimum_block_length_to_tau_ratio", 0.0)) != 16.0 or \
            float(statistics.get("maximum_block_stability_ratio", 0.0)) != \
            1.25 or \
            float(statistics.get("minimum_leave_one_denominator", 0.0)) != \
            1.0e-12 or \
            float(statistics.get("minimum_abs_denominator_mean", 0.0)) != \
            1.0e-12 or \
            float(statistics.get(
                "maximum_denominator_relative_se", 0.0)) != 1.0 or \
            int(session.get("persistent_session", 0)) != 1 or \
            int(session.get("session_root_evaluations", 0)) != 32769 or \
            int(session.get("cache_reset_count", 0)) != 1 or \
            int(session.get("engine_heap_allocations", -1)) != 0 or \
            int(session.get("session_end_pass", 0)) != 1:
        raise AssertionError("P4-S9 official statistical gate mismatch")
    parents = policy.get("parent_sha256", {})
    expected_parent_names = {
        "long_direct_policy", "session_resource_policy",
        "target_session_policy", "stage_a_evidence", "stage_a_decision",
        "stage_c_local_decision", "target_local_evidence",
        "target_local_decision",
    }
    if set(parents) != expected_parent_names or any(
            len(str(parents[name])) != 64 for name in expected_parent_names):
        raise AssertionError("P4-S9 official parent hash set mismatch")
    return policy


def parent_paths(args):
    return {
        "long_direct_policy": args.long_direct_policy,
        "session_resource_policy": args.session_resource_policy,
        "target_session_policy": args.target_session_policy,
        "stage_a_evidence": args.stage_a_evidence,
        "stage_a_decision": args.stage_a_decision,
        "stage_c_local_decision": args.stage_c_local_decision,
        "target_local_evidence": args.target_local_evidence,
        "target_local_decision": args.target_local_decision,
    }


def verify_parent_bindings(policy, args, allow_smoke_policy=False):
    bindings = {
        name: p4c.sha256_file(path) for name, path in parent_paths(args).items()
    }
    if not allow_smoke_policy:
        for name, digest in bindings.items():
            if policy["parent_sha256"][name] != digest:
                raise AssertionError(
                    "P4-S9 official parent mismatch {}".format(name))
        stage_a = p4c.read_json(args.stage_a_decision)
        stage_c = p4c.read_json(args.stage_c_local_decision)
        target = p4c.read_json(args.target_local_decision)
        if stage_a.get("stage_a_decision") != "ELIGIBLE" or \
                stage_c.get("local_resource_decision") != "PASS" or \
                target.get("local_target_stress_decision") != "PASS":
            raise AssertionError("P4-S9 local gates did not authorize HPC")
    return bindings


def parse_session_log(path):
    run = p4s.parse_markov_log(path)
    session = None
    with open(path, "r") as handle:
        for line in handle:
            if line.startswith("SESSION "):
                if session is not None:
                    raise AssertionError("duplicate P4-S9 SESSION record")
                session = p4s.parse_key_values(line, 1)
    if session is None:
        raise AssertionError("missing P4-S9 SESSION record")
    run["session"] = session
    return run


def compatibility_policy(statistics, seed_hex):
    return {
        "scope": {
            "site_counts": statistics["site_counts"],
            "qp_total": statistics["qp_total"],
            "cache_capacity_bytes": [statistics["cache_bytes"]],
        },
        "guide": {"rho": [statistics["rho"]]},
        "p4_s_markov_gate": {
            "sample_count": statistics["sample_count"],
            "official_partition": statistics["official_partition"],
            "diagnostic_partition": statistics["diagnostic_partition"],
            "initialization": {"seed_hex": seed_hex},
            "denominator_stability": {
                "minimum_leave_one_denominator":
                    statistics["minimum_leave_one_denominator"],
                "minimum_abs_denominator_mean":
                    statistics["minimum_abs_denominator_mean"],
                "maximum_denominator_relative_se":
                    statistics["maximum_denominator_relative_se"],
            },
            "autocorrelation": {
                "maximum_tau_int": statistics["maximum_tau_int"],
                "minimum_official_block_length_to_tau_ratio":
                    statistics["minimum_block_length_to_tau_ratio"],
            },
            "block_stability": {
                "diagnostic_maximum":
                    statistics["maximum_block_stability_ratio"],
            },
            "conservative_precision": {
                "maximum":
                    statistics["maximum_conservative_se_budget_ratio"],
            },
        },
    }


def require_close(actual, expected, label, tolerance=1.0e-10):
    actual = float(actual)
    expected = float(expected)
    if not math.isfinite(actual) or not math.isfinite(expected) or \
            abs(actual - expected) > tolerance * max(1.0, abs(expected)):
        raise AssertionError("{} mismatch: {} != {}".format(
            label, actual, expected))


def validate_gate_consistency(run, statistics):
    """Reapply the frozen gate to the producer's raw summary records."""
    header = run["header"]
    sample_count = int(header["sample_count"])
    official_block_length = int(header["official_block_length"])
    maximum_tau = float(statistics["maximum_tau_int"])
    minimum_block_multiplier = float(
        statistics["minimum_block_length_to_tau_ratio"])

    expected_tau = {("W", "real")}
    for family in p4s.FAMILIES:
        for row, column in p4s.UPPER_ENTRIES:
            series = "{}_{}{}".format(family, row, column)
            expected_tau.add((series, "real"))
            expected_tau.add((series, "imag"))
    observed_tau = [(str(item.get("series", "")),
                     str(item.get("component", "")))
                    for item in run["tau"]]
    if len(observed_tau) != len(expected_tau) or \
            set(observed_tau) != expected_tau:
        raise AssertionError("P4-S9 official TAU census mismatch")

    maximum_tau_seen = 0.0
    tau_pass = True
    for item in run["tau"]:
        variance = float(item["variance"])
        tau_int = float(item["tau_int"])
        block_length = float(item["block_length"])
        block_multiplier = float(item["block_multiplier"])
        required = int(item["required"])
        passed = int(item["pass"])
        if int(item["sample_count"]) != sample_count or \
                not math.isfinite(variance) or variance < 0.0 or \
                int(item["positive_pair_count"]) < 0 or \
                not math.isfinite(tau_int) or tau_int < 1.0 or \
                required not in (0, 1) or passed not in (0, 1):
            raise AssertionError("P4-S9 official TAU record mismatch")
        expected_required = int(
            item["series"] == "W" or variance != 0.0)
        expected_pass = int(
            tau_int <= maximum_tau and
            block_length >= minimum_block_multiplier * tau_int)
        require_close(block_length, official_block_length,
                      "P4-S9 official TAU block length")
        require_close(block_multiplier, block_length / tau_int,
                      "P4-S9 official TAU block multiplier")
        if required != expected_required or passed != expected_pass:
            raise AssertionError("P4-S9 official TAU gate mismatch")
        maximum_tau_seen = max(maximum_tau_seen, tau_int)
        tau_pass = tau_pass and (not required or bool(passed))

    entries = run["entries"]
    observed_entries = [
        (item["family"], int(item["row"]), int(item["column"]))
        for item in entries
    ]
    expected_entries = {
        (family, row, column)
        for family in p4s.FAMILIES
        for row, column in p4s.UPPER_ENTRIES
    }
    if len(observed_entries) != len(expected_entries) or \
            set(observed_entries) != expected_entries:
        raise AssertionError("P4-S9 official ENTRY census mismatch")

    denominators = [float(item["denominator"])
                    for item in run["samples"]]
    if len(denominators) != sample_count or any(
            not math.isfinite(value) or value < 0.0
            for value in denominators):
        raise AssertionError("P4-S9 official denominator sample mismatch")
    denominator_sum = math.fsum(denominators)
    denominator_squared_sum = math.fsum(
        value * value for value in denominators)
    denominator_mean = denominator_sum / sample_count
    if denominator_mean <= 0.0 or denominator_squared_sum <= 0.0:
        raise AssertionError("P4-S9 official nonpositive denominator")
    denominator_ss = (
        denominator_squared_sum -
        2.0 * denominator_mean * denominator_sum +
        sample_count * denominator_mean * denominator_mean)
    if denominator_ss < 0.0 and denominator_ss >= \
            -1.0e-10 * max(1.0, denominator_squared_sum):
        denominator_ss = 0.0
    if not math.isfinite(denominator_ss) or denominator_ss < 0.0:
        raise AssertionError("P4-S9 official denominator variance mismatch")
    denominator_relative_se = 0.0
    if sample_count > 1:
        denominator_relative_se = math.sqrt(
            denominator_ss / ((sample_count - 1.0) * sample_count)) / \
            abs(denominator_mean)
    effective_fraction = (
        denominator_sum * denominator_sum /
        (denominator_squared_sum * sample_count))
    zero_target_count = sum(
        int(item["zero_target"]) for item in run["samples"])
    if any(int(item["zero_target"]) not in (0, 1)
           for item in run["samples"]):
        raise AssertionError("P4-S9 official zero-target flag mismatch")

    summary = run["summary"]
    require_close(summary["denominator_sum"], denominator_sum,
                  "P4-S9 official denominator sum")
    require_close(summary["denominator_mean"], denominator_mean,
                  "P4-S9 official denominator mean")
    require_close(summary["denominator_relative_se"],
                  denominator_relative_se,
                  "P4-S9 official denominator relative SE")
    require_close(summary["effective_sample_fraction"], effective_fraction,
                  "P4-S9 official effective sample fraction")
    require_close(summary["zero_target_sample_fraction"],
                  float(zero_target_count) / sample_count,
                  "P4-S9 official zero-target fraction")
    require_close(summary["minimum_denominator"], min(denominators),
                  "P4-S9 official minimum denominator")
    require_close(summary["maximum_denominator"], max(denominators),
                  "P4-S9 official maximum denominator")
    require_close(summary["denominator_tail_ratio"],
                  max(denominators) / denominator_mean,
                  "P4-S9 official denominator tail ratio")
    if int(summary["zero_target_sample_count"]) != zero_target_count:
        raise AssertionError("P4-S9 official zero-target census mismatch")
    expected_summary_stable = (
        abs(denominator_mean) >=
        float(statistics["minimum_abs_denominator_mean"]) and
        denominator_relative_se <=
        float(statistics["maximum_denominator_relative_se"]))
    if int(summary["denominator_stable"]) not in (0, 1) or \
            bool(int(summary["denominator_stable"])) != \
            expected_summary_stable:
        raise AssertionError("P4-S9 official denominator summary gate mismatch")

    entry_denominator_pass = True
    for entry in entries:
        stable = int(entry["denominator_stable"])
        if stable not in (0, 1) or \
                int(entry["unstable_leave_one_blocks"]) < 0:
            raise AssertionError("P4-S9 official ENTRY denominator mismatch")
        entry_denominator_pass = entry_denominator_pass and bool(stable)
    denominator_pass = expected_summary_stable and entry_denominator_pass
    support_pass = int(run["trace"]["support_violation"]) == 0
    decision = run["decision"]
    if bool(int(decision["tau_pass"])) != tau_pass or \
            bool(int(decision["support_pass"])) != support_pass or \
            bool(int(decision["denominator_pass"])) != denominator_pass:
        raise AssertionError("P4-S9 official raw/decision gate mismatch")
    require_close(decision["maximum_tau_int"], maximum_tau_seen,
                  "P4-S9 official maximum tau")
    return run


def validate_run(run, policy, allow_smoke_policy=False):
    statistics = policy["official_statistics"]
    header = run["header"]
    site = int(header.get("site_count", 0))
    qp = int(header.get("qp_total", 0))
    sample_count = int(header.get("sample_count", 0))
    seed_hex = str(header.get("seed_hex", "")).lower()
    allowed_seeds = [value.lower() for value in statistics["seed_hex_order"]]
    if seed_hex not in allowed_seeds:
        raise AssertionError("P4-S9 official seed mismatch")
    try:
        seed = int(seed_hex, 16)
    except ValueError:
        raise AssertionError("P4-S9 official seed encoding mismatch")
    if int(header.get("schema", 0)) != DRIVER_SCHEMA_VERSION or \
            header.get("fixture") != EXPECTED_FIXTURE or \
            [site, qp] not in statistics["physical_case_order"] or \
            int(header.get("rank_count", 0)) != 1 or \
            int(header.get("cache_requested_bytes", 0)) != \
            int(statistics["cache_bytes"]) or \
            float(header.get("rho", 0.0)) != float(statistics["rho"]) or \
            int(header.get("global_numerator", -1)) != 0 or \
            int(header.get("global_denominator", 0)) != 1 or \
            int(header.get("rng_stream", -1)) != \
            int(statistics["rng_stream"]) or \
            int(header.get("persistent_session", 0)) != 1 or \
            int(header.get("amplitude_generation_hash", 0)) == 0 or \
            int(header.get("seed", -1)) != seed or \
            (not allow_smoke_policy and sample_count !=
             int(statistics["sample_count"])):
        raise AssertionError("P4-S9 official Markov header mismatch")
    compatible = copy.deepcopy(run)
    compatible["header"]["schema"] = p4s.SCHEMA_VERSION
    compatible["header"]["fixture"] = p4s.EXPECTED_MARKOV_FIXTURE
    p4s.validate_markov_run(
        compatible, compatibility_policy(statistics, seed_hex),
        allow_smoke_policy)
    validate_gate_consistency(run, statistics)
    session = run["session"]
    if int(session.get("schema", 0)) != DRIVER_SCHEMA_VERSION or \
            int(session.get("amplitude_generation_hash", 0)) != \
            int(header["amplitude_generation_hash"]) or \
            int(session.get("session_root_evaluations", 0)) != \
            sample_count + 1 or \
            int(session.get("cache_reset_count", 0)) != 1 or \
            int(session.get("engine_heap_allocations", -1)) != 0 or \
            int(session.get("session_end_pass", 0)) != 1 or \
            int(session.get("cache_resident_end", 0)) < 0 or \
            int(session.get("cache_entries_peak", 0)) < \
            int(session.get("cache_resident_end", 0)):
        raise AssertionError("P4-S9 official session lifecycle mismatch")
    return run


def materialize_inputs(output, args, runs):
    frozen = os.path.join(output, "frozen_inputs")
    producer = os.path.join(output, "producer")
    if not os.path.isdir(frozen):
        os.makedirs(frozen)
    if not os.path.isdir(producer):
        os.makedirs(producer)
    for path in parent_paths(args).values():
        shutil.copyfile(path, os.path.join(frozen, os.path.basename(path)))
    shutil.copyfile(args.policy,
                    os.path.join(producer, os.path.basename(args.policy)))
    shutil.copyfile(args.producer_binary, os.path.join(
        producer, os.path.basename(args.producer_binary)))
    raw_paths = p4s.materialize_raw_files(
        [run["path"] for run in runs], output)
    return raw_paths


def run_sort_key(run, policy):
    statistics = policy["official_statistics"]
    return (
        [value.lower() for value in statistics["seed_hex_order"]].index(
            str(run["header"]["seed_hex"]).lower()),
        statistics["physical_case_order"].index([
            int(run["header"]["site_count"]),
            int(run["header"]["qp_total"])]),
    )


def analyze_command(args):
    for label, commit in (("source", args.source_commit),
                          ("baseline", args.baseline_commit)):
        if re.fullmatch(r"[0-9a-f]{40}", commit or "") is None:
            raise AssertionError(
                "P4-S9 {} commit provenance mismatch".format(label))
    policy = validate_policy(p4c.read_json(args.policy))
    bindings = verify_parent_bindings(
        policy, args, args.allow_smoke_policy)
    if int(args.omp_threads) != 1 or int(args.blas_threads) != 1:
        raise AssertionError("P4-S9 official thread metadata mismatch")
    if args.repository:
        p4c.validate_commit_object(args.repository, args.source_commit)
        p4c.validate_commit_object(args.repository, args.baseline_commit)
    runs = [validate_run(parse_session_log(path), policy,
                         args.allow_smoke_policy)
            for path in args.markov_log]
    statistics_policy = policy["official_statistics"]
    observed = {
        (int(run["header"]["site_count"]),
         int(run["header"]["qp_total"]),
         str(run["header"]["seed_hex"]).lower()) for run in runs
    }
    if len(observed) != len(runs):
        raise AssertionError("duplicate P4-S9 official statistical run")
    if not args.allow_smoke_policy:
        expected = {
            (site, qp, seed.lower())
            for site, qp in statistics_policy["physical_case_order"]
            for seed in statistics_policy["seed_hex_order"]
        }
        if observed != expected:
            raise AssertionError("P4-S9 official statistical census mismatch")
    output = os.path.abspath(args.output)
    if not os.path.isdir(output):
        os.makedirs(output)
    raw_paths = materialize_inputs(output, args, runs)
    ordered = sorted(runs, key=lambda run: run_sort_key(run, policy))
    cases = []
    all_pass = True
    for run in ordered:
        header = run["header"]
        decision = run["decision"]
        passed = decision["p4s_decision"] == "GO"
        all_pass = all_pass and passed
        cases.append({
            "site_count": int(header["site_count"]),
            "qp_total": int(header["qp_total"]),
            "seed_hex": str(header["seed_hex"]).lower(),
            "sample_count": int(header["sample_count"]),
            "rank_count": int(header["rank_count"]),
            "amplitude_generation_hash":
                int(header["amplitude_generation_hash"]),
            "raw_log": raw_paths[run["path"]],
            "raw_log_sha256": run["sha256"],
            "raw_log_line_count": run["line_count"],
            "trace": run["trace"],
            "summary": run["summary"],
            "session": run["session"],
            "manifest": run["manifest"],
            "decision": decision,
            "maximum_tau_int": float(decision["maximum_tau_int"]),
            "maximum_se_budget_ratio":
                float(decision["maximum_se_budget_ratio"]),
            "official_statistics_pass": passed,
        })
    generated_at = p4c.generated_at_jst()
    metadata = {
        "schema_version": POLICY_SCHEMA_VERSION,
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
        "producer_binary_sha256": p4c.sha256_file(args.producer_binary),
        "parent_sha256": bindings,
        "candidate_id": policy["candidate_id"],
    }
    evidence = dict(metadata)
    evidence.update({
        "artifact": "p4s9_official_markov_statistics",
        "case_count": len(cases),
        "cases": cases,
        "official_statistics_pass": all_pass,
    })
    decision = dict(metadata)
    decision.update({
        "artifact": "p4s9_official_statistics_decision",
        "official_statistics_decision": "GO" if all_pass else "STOP",
        "p4s9_combined_decision": "PENDING_OFFICIAL_RESOURCE_AND_MPI_PARITY"
            if all_pass else "STOP",
        "p4f_authorized": False,
        "failure_action": "continue_to_combined_closure"
            if all_pass else "preserve_evidence_and_stop_before_p4f",
    })
    p4c.write_json(os.path.join(
        output, "p4s9_official_markov_statistics.json"), evidence)
    p4c.write_json(os.path.join(
        output, "p4s9_official_statistics_decision.json"), decision)
    count = p4c.write_checksums(output)
    print("P4-S9 official statistics: {}; {} ledger files".format(
        decision["official_statistics_decision"], count))


def build_parser():
    parser = argparse.ArgumentParser()
    p4c.add_common_arguments(parser)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--long-direct-policy", required=True)
    parser.add_argument("--session-resource-policy", required=True)
    parser.add_argument("--target-session-policy", required=True)
    parser.add_argument("--stage-a-evidence", required=True)
    parser.add_argument("--stage-a-decision", required=True)
    parser.add_argument("--stage-c-local-decision", required=True)
    parser.add_argument("--target-local-evidence", required=True)
    parser.add_argument("--target-local-decision", required=True)
    parser.add_argument("--producer-binary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--markov-log", nargs="+", required=True)
    parser.add_argument("--allow-smoke-policy", action="store_true")
    return parser


def main(argv=None):
    analyze_command(build_parser().parse_args(argv))


if __name__ == "__main__":
    main()
