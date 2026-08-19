from __future__ import print_function

import argparse
import os
import re
import shutil
import subprocess

import generate_bounded_krylov_p4c_evidence as p4c


SCHEMA_VERSION = 1
EXPECTED_MARKOV_FIXTURE = "p4s_bounded_markov_real"
FAMILIES = ("S", "K", "B")
UPPER_ENTRIES = tuple(
    (row, column) for row in range(3) for column in range(row, 3))


def read_json(path):
    return p4c.read_json(path)


def write_json(path, value):
    p4c.write_json(path, value)


def sha256_file(path):
    return p4c.sha256_file(path)


def parse_number(value):
    if value in ("GO", "STOP", "PASS", "FAIL"):
        return value
    return p4c.parse_number(value)


def parse_key_values(line, skip):
    return {
        key: parse_number(value)
        for key, value in p4c.parse_key_values(line, skip).items()
    }


def guide_key(rho):
    return "r3_rho_{:.0e}".format(float(rho))


def candidate_key(cache_bytes, rho):
    return "cache{}_{}".format(int(cache_bytes), guide_key(rho))


def validate_p4s_policy(policy, allow_smoke_policy=False):
    p4c.validate_policy_contract(policy, allow_smoke_policy)
    gate = policy.get("p4_s_markov_gate", {})
    official = gate.get("official_partition", {})
    diagnostic = gate.get("diagnostic_partition", {})
    autocorr = gate.get("autocorrelation", {})
    stability = gate.get("block_stability", {})
    precision = gate.get("conservative_precision", {})
    denominator = gate.get("denominator_stability", {})
    if int(gate.get("sample_count", 0)) != 4096:
        raise AssertionError("P4-S sample count mismatch")
    if (int(official.get("block_count", 0)) != 16 or
            int(official.get("block_length", 0)) != 256):
        raise AssertionError("P4-S official partition mismatch")
    if (int(diagnostic.get("block_count", 0)) != 32 or
            int(diagnostic.get("block_length", 0)) != 128):
        raise AssertionError("P4-S diagnostic partition mismatch")
    if float(stability.get("diagnostic_maximum", 0.0)) < 1.0:
        raise AssertionError("invalid P4-S block-stability diagnostic")
    if float(precision.get("maximum", 0.0)) <= 0.0:
        raise AssertionError("invalid P4-S conservative precision gate")
    if float(autocorr.get("maximum_tau_int", 0.0)) <= 0.0:
        raise AssertionError("invalid P4-S tau threshold")
    if float(autocorr.get(
            "minimum_official_block_length_to_tau_ratio", 0.0)) <= 0.0:
        raise AssertionError("invalid P4-S tau block multiplier")
    for key in ("minimum_abs_denominator_mean",
                "minimum_leave_one_denominator",
                "maximum_denominator_relative_se"):
        if float(denominator.get(key, 0.0)) <= 0.0:
            raise AssertionError("invalid P4-S denominator gate")
    seed_hex = gate.get("initialization", {}).get("seed_hex")
    if not isinstance(seed_hex, str) or \
            re.match(r"^0x[0-9A-Fa-f]+$", seed_hex) is None:
        raise AssertionError("invalid P4-S seed")
    return policy


def parse_markov_log(path):
    header = None
    scales = []
    samples = []
    trace = None
    summary = None
    tau = []
    entries = []
    manifest = None
    decision = None
    with open(path, "r") as handle:
        text = handle.read()
    for line in text.splitlines():
        if line.startswith("MARKOV "):
            if header is not None:
                raise AssertionError("duplicate MARKOV header")
            header = parse_key_values(line, 1)
        elif line.startswith("SCALE "):
            scales.append(parse_key_values(line, 1))
        elif line.startswith("SAMPLE "):
            samples.append(parse_key_values(line, 1))
        elif line.startswith("TRACE "):
            if trace is not None:
                raise AssertionError("duplicate TRACE")
            trace = parse_key_values(line, 1)
        elif line.startswith("SUMMARY "):
            if summary is not None:
                raise AssertionError("duplicate SUMMARY")
            summary = parse_key_values(line, 1)
        elif line.startswith("TAU "):
            tau.append(parse_key_values(line, 1))
        elif line.startswith("ENTRY "):
            entries.append(parse_key_values(line, 1))
        elif line.startswith("MANIFEST "):
            if manifest is not None:
                raise AssertionError("duplicate MANIFEST")
            manifest = parse_key_values(line, 1)
        elif line.startswith("DECISION "):
            if decision is not None:
                raise AssertionError("duplicate DECISION")
            decision = parse_key_values(line, 1)
    if (header is None or trace is None or summary is None or
            manifest is None or decision is None):
        raise AssertionError("incomplete Markov log {}".format(path))
    return {
        "path": os.path.abspath(path),
        "sha256": sha256_file(path),
        "header": header,
        "scales": scales,
        "samples": samples,
        "trace": trace,
        "summary": summary,
        "tau": tau,
        "entries": entries,
        "manifest": manifest,
        "decision": decision,
        "line_count": len(text.splitlines()),
    }


def require_close(actual, expected, label, tolerance=1.0e-12):
    actual = float(actual)
    expected = float(expected)
    if abs(actual - expected) > tolerance * max(1.0, abs(expected)):
        raise AssertionError("{} mismatch: {} != {}".format(
            label, actual, expected))


def validate_markov_run(run, policy, allow_smoke_policy=False):
    header = run["header"]
    gate = policy["p4_s_markov_gate"]
    scope = policy["scope"]
    if int(header.get("schema", 0)) != SCHEMA_VERSION:
        raise AssertionError("Markov schema mismatch")
    if header.get("fixture") != EXPECTED_MARKOV_FIXTURE:
        raise AssertionError("Markov fixture mismatch")
    site = int(header["site_count"])
    qp = int(header["qp_total"])
    cache = int(header["cache_requested_bytes"])
    rho = float(header["rho"])
    sample_count = int(header["sample_count"])
    if site not in scope["site_counts"] or qp not in scope["qp_total"]:
        raise AssertionError("Markov run outside P4-S site/QP scope")
    if cache not in scope["cache_capacity_bytes"]:
        raise AssertionError("Markov run outside cache scope")
    if rho not in [float(value) for value in policy["guide"]["rho"]]:
        raise AssertionError("Markov run outside rho scope")
    if not allow_smoke_policy and sample_count != int(gate["sample_count"]):
        raise AssertionError("official Markov sample count mismatch")
    if len(run["scales"]) != 4:
        raise AssertionError("scale vector must cover orders 0..3")
    observed_orders = [int(item["order"]) for item in run["scales"]]
    if observed_orders != [0, 1, 2, 3]:
        raise AssertionError("scale order mismatch")
    seed_hex = str(header["seed_hex"]).lower()
    if seed_hex != gate["initialization"]["seed_hex"].lower():
        raise AssertionError("P4-S seed mismatch")
    if int(header["official_block_count"]) != \
            int(gate["official_partition"]["block_count"]):
        raise AssertionError("official block count mismatch")
    if int(header["diagnostic_block_count"]) != \
            int(gate["diagnostic_partition"]["block_count"]):
        raise AssertionError("diagnostic block count mismatch")
    if not allow_smoke_policy:
        if int(header["official_block_length"]) != \
                int(gate["official_partition"]["block_length"]):
            raise AssertionError("official block length mismatch")
        if int(header["diagnostic_block_length"]) != \
                int(gate["diagnostic_partition"]["block_length"]):
            raise AssertionError("diagnostic block length mismatch")
    denominator_gate = gate["denominator_stability"]
    require_close(header["minimum_leave_one_denominator"],
                  denominator_gate["minimum_leave_one_denominator"],
                  "minimum leave-one denominator")
    require_close(header["minimum_abs_denominator_mean"],
                  denominator_gate["minimum_abs_denominator_mean"],
                  "minimum denominator mean")
    require_close(header["maximum_denominator_relative_se"],
                  denominator_gate["maximum_denominator_relative_se"],
                  "maximum denominator relative SE")
    require_close(header["maximum_tau_int"],
                  gate["autocorrelation"]["maximum_tau_int"],
                  "maximum tau")
    require_close(header["block_length_multiplier"],
                  gate["autocorrelation"][
                      "minimum_official_block_length_to_tau_ratio"],
                  "tau block multiplier")
    require_close(header["maximum_block_stability_ratio"],
                  gate["block_stability"]["diagnostic_maximum"],
                  "block stability diagnostic maximum")
    require_close(header["maximum_conservative_se_budget_ratio"],
                  gate["conservative_precision"]["maximum"],
                  "conservative SE budget maximum")

    samples = run["samples"]
    if len(samples) != sample_count:
        raise AssertionError("SAMPLE count mismatch")
    for index, sample in enumerate(samples):
        if int(sample["sample"]) != index:
            raise AssertionError("SAMPLE schedule mismatch")
    trace = run["trace"]
    if int(trace["completed"]) != sample_count:
        raise AssertionError("trace completed count mismatch")
    if int(trace["attempted"]) != int(trace["completed"]):
        raise AssertionError("trace attempted/completed mismatch")
    if int(trace["accepted"]) + int(trace["rejected"]) != sample_count:
        raise AssertionError("accept/reject count mismatch")
    if int(trace["positive_support"]) + int(trace["support_violation"]) != \
            sample_count:
        raise AssertionError("support count mismatch")
    summary = run["summary"]
    if int(summary["sample_count"]) != sample_count:
        raise AssertionError("summary sample count mismatch")
    if int(summary["denominator_stable"]) not in (0, 1):
        raise AssertionError("invalid denominator stability flag")

    expected_entries = {
        (family, row, column)
        for family in FAMILIES for row, column in UPPER_ENTRIES
    }
    observed_entries = {
        (entry["family"], int(entry["row"]), int(entry["column"]))
        for entry in run["entries"]
    }
    if observed_entries != expected_entries:
        raise AssertionError("required entry set mismatch")
    tau_series = {(item["series"], item["component"]) for item in run["tau"]}
    if ("W", "real") not in tau_series:
        raise AssertionError("W tau series missing")
    for family, row, column in expected_entries:
        name = "{}_{}{}".format(family, row, column)
        if (name, "real") not in tau_series or (name, "imag") not in tau_series:
            raise AssertionError("influence tau series missing: {}".format(name))
    precision_gate = gate["conservative_precision"]
    stability_gate = gate["block_stability"]
    maximum_official_se_budget_ratio = 0.0
    maximum_conservative_se_budget_ratio = 0.0
    maximum_block_stability_ratio = 0.0
    entry_budget_pass = True
    entry_block_stability_pass = True
    entry_pathology_pass = True
    for entry in run["entries"]:
        official_se = float(entry["official_se"])
        diagnostic_se = float(entry["diagnostic_se"])
        official_ratio = float(entry["se_budget_ratio"])
        diagnostic_ratio = float(entry["diagnostic_se_budget_ratio"])
        conservative_ratio = float(entry["conservative_se_budget_ratio"])
        stability_ratio = float(entry["stability_ratio"])
        if official_se < 0.0 or diagnostic_se < 0.0:
            raise AssertionError("negative entry standard error")
        require_close(
            conservative_ratio, max(official_ratio, diagnostic_ratio),
            "entry conservative SE budget ratio")
        expected_budget_pass = (
            conservative_ratio <= float(precision_gate["maximum"]))
        if bool(int(entry["budget_pass"])) != expected_budget_pass:
            raise AssertionError("entry budget gate mismatch")
        maximum_official_se_budget_ratio = max(
            maximum_official_se_budget_ratio, official_ratio)
        maximum_conservative_se_budget_ratio = max(
            maximum_conservative_se_budget_ratio, conservative_ratio)
        maximum_block_stability_ratio = max(
            maximum_block_stability_ratio, stability_ratio)
        entry_budget_pass = entry_budget_pass and expected_budget_pass
        expected_stability_pass = (
            stability_ratio <= float(stability_gate["diagnostic_maximum"]))
        if bool(int(entry["stability_pass"])) != expected_stability_pass:
            raise AssertionError("entry block-stability diagnostic mismatch")
        entry_block_stability_pass = (
            entry_block_stability_pass and expected_stability_pass)
        one_zero_only = (
            max(official_se, diagnostic_se) > 0.0 and
            min(official_se, diagnostic_se) == 0.0)
        expected_pathology_pass = not one_zero_only
        if bool(int(entry["pathology_pass"])) != expected_pathology_pass:
            raise AssertionError("entry block pathology gate mismatch")
        entry_pathology_pass = (
            entry_pathology_pass and expected_pathology_pass)
    decision = run["decision"]
    if bool(int(decision["budget_pass"])) != entry_budget_pass:
        raise AssertionError("P4-S conservative budget summary mismatch")
    if bool(int(decision["block_stability_pass"])) != \
            entry_block_stability_pass:
        raise AssertionError("P4-S block-stability diagnostic mismatch")
    if bool(int(decision["block_pathology_pass"])) != entry_pathology_pass:
        raise AssertionError("P4-S block pathology summary mismatch")
    require_close(decision["maximum_se_budget_ratio"],
                  maximum_conservative_se_budget_ratio,
                  "maximum conservative SE budget ratio")
    require_close(decision["maximum_official_se_budget_ratio"],
                  maximum_official_se_budget_ratio,
                  "maximum official SE budget ratio")
    require_close(decision["maximum_block_stability_ratio"],
                  maximum_block_stability_ratio,
                  "maximum block stability ratio")
    expected_pass = (
        int(decision["support_pass"]) and int(decision["tau_pass"]) and
        int(decision["budget_pass"]) and
        int(decision["block_pathology_pass"]) and
        int(decision["denominator_pass"]))
    if (decision["p4s_decision"] == "GO") != bool(expected_pass):
        raise AssertionError("P4-S decision does not match AND gate")
    return run


def materialize_raw_files(paths, output):
    destination_root = os.path.join(os.path.abspath(output), "raw_markov")
    if not os.path.isdir(destination_root):
        os.makedirs(destination_root)
    result = {}
    used = set()
    for index, source in enumerate(paths):
        digest = sha256_file(source)
        name = "{:04d}-{}-{}".format(
            index, digest[:12], os.path.basename(source))
        if name in used:
            raise AssertionError("raw Markov evidence name collision")
        used.add(name)
        destination = os.path.join(destination_root, name)
        if os.path.abspath(source) != os.path.abspath(destination):
            shutil.copyfile(source, destination)
        if sha256_file(destination) != digest:
            raise AssertionError("raw Markov checksum mismatch")
        result[os.path.abspath(source)] = os.path.join("raw_markov", name)
    return result


def validate_commit_object(repository, commit):
    return p4c.validate_commit_object(repository, commit)


def common_metadata(args, policy_hash):
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at_jst": p4c.generated_at_jst(),
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


def candidate_sort_key(run, policy):
    header = run["header"]
    scope = policy["scope"]
    return (
        scope["cache_capacity_bytes"].index(int(header["cache_requested_bytes"])),
        [float(value) for value in policy["guide"]["rho"]].index(
            float(header["rho"])),
        int(header["site_count"]),
        int(header["qp_total"]),
    )


def summarize_run(run, raw_files):
    header = run["header"]
    key = candidate_key(header["cache_requested_bytes"], header["rho"])
    return {
        "candidate": key,
        "site_count": int(header["site_count"]),
        "qp_total": int(header["qp_total"]),
        "cache_requested_bytes": int(header["cache_requested_bytes"]),
        "rho": float(header["rho"]),
        "eta": float(header["eta"]),
        "sample_count": int(header["sample_count"]),
        "raw_log": raw_files[run["path"]],
        "raw_log_sha256": run["sha256"],
        "trace": run["trace"],
        "summary": run["summary"],
        "manifest": run["manifest"],
        "decision": run["decision"],
        "maximum_tau_int": float(run["decision"]["maximum_tau_int"]),
        "maximum_se_budget_ratio":
            float(run["decision"]["maximum_se_budget_ratio"]),
        "maximum_official_se_budget_ratio":
            float(run["decision"]["maximum_official_se_budget_ratio"]),
        "maximum_block_stability_ratio":
            float(run["decision"]["maximum_block_stability_ratio"]),
        "entries": run["entries"],
        "tau": run["tau"],
    }


def build_candidate_summary(runs, policy):
    grouped = {}
    for run in runs:
        header = run["header"]
        key = candidate_key(header["cache_requested_bytes"], header["rho"])
        summary = grouped.setdefault(key, {
            "cache_requested_bytes": int(header["cache_requested_bytes"]),
            "r": 3,
            "rho": float(header["rho"]),
            "p4s_markov_gate_pass": True,
            "support_pass": True,
            "tau_pass": True,
            "block_stability_pass": True,
            "block_pathology_pass": True,
            "budget_pass": True,
            "denominator_pass": True,
            "maximum_tau_int": 0.0,
            "maximum_se_budget_ratio": 0.0,
            "maximum_official_se_budget_ratio": 0.0,
            "maximum_block_stability_ratio": 0.0,
            "cases": [],
        })
        decision = run["decision"]
        flags = {
            "support_pass": bool(int(decision["support_pass"])),
            "tau_pass": bool(int(decision["tau_pass"])),
            "block_stability_pass":
                bool(int(decision["block_stability_pass"])),
            "block_pathology_pass":
                bool(int(decision["block_pathology_pass"])),
            "budget_pass": bool(int(decision["budget_pass"])),
            "denominator_pass": bool(int(decision["denominator_pass"])),
        }
        case_pass = decision["p4s_decision"] == "GO"
        summary["p4s_markov_gate_pass"] = (
            summary["p4s_markov_gate_pass"] and case_pass)
        for name, passed in flags.items():
            summary[name] = summary[name] and passed
        summary["maximum_tau_int"] = max(
            summary["maximum_tau_int"], float(decision["maximum_tau_int"]))
        summary["maximum_se_budget_ratio"] = max(
            summary["maximum_se_budget_ratio"],
            float(decision["maximum_se_budget_ratio"]))
        summary["maximum_official_se_budget_ratio"] = max(
            summary["maximum_official_se_budget_ratio"],
            float(decision["maximum_official_se_budget_ratio"]))
        summary["maximum_block_stability_ratio"] = max(
            summary["maximum_block_stability_ratio"],
            float(decision["maximum_block_stability_ratio"]))
        summary["cases"].append({
            "site_count": int(header["site_count"]),
            "qp_total": int(header["qp_total"]),
            "decision": decision["p4s_decision"],
        })
    promoted = [key for key, value in grouped.items()
                if value["p4s_markov_gate_pass"]]
    promoted = sorted(promoted, key=lambda key: (
        policy["scope"]["cache_capacity_bytes"].index(
            grouped[key]["cache_requested_bytes"]),
        [float(value) for value in policy["guide"]["rho"]].index(
            grouped[key]["rho"]),
    ))
    return grouped, promoted


def write_checksums(output):
    return p4c.write_checksums(output)


def analyze_command(args):
    policy = validate_p4s_policy(
        read_json(args.policy), args.allow_smoke_policy)
    policy_hash = sha256_file(args.policy)
    if args.repository:
        validate_commit_object(args.repository, args.source_commit)
        validate_commit_object(args.repository, args.baseline_commit)
    if int(args.omp_threads) != int(policy["scope"]["omp_threads"]) or \
            int(args.blas_threads) != int(policy["scope"]["blas_threads"]):
        raise AssertionError("thread metadata does not match P4 policy")
    runs = [
        validate_markov_run(parse_markov_log(path), policy,
                            args.allow_smoke_policy)
        for path in args.markov_log
    ]
    output = os.path.abspath(args.output)
    if not os.path.isdir(output):
        os.makedirs(output)
    raw_files = materialize_raw_files([run["path"] for run in runs], output)
    ordered_runs = sorted(runs, key=lambda run: candidate_sort_key(run, policy))
    statistics = {
        "artifact": "p4s_markov_statistics",
        "measurement_policy": policy,
        "raw_logs": [summarize_run(run, raw_files) for run in ordered_runs],
    }
    statistics.update(common_metadata(args, policy_hash))
    candidates, promoted = build_candidate_summary(ordered_runs, policy)
    decision = {
        "artifact": "p4s_candidate_decision",
        "measurement_policy": policy,
        "candidate_summary": candidates,
        "promoted_candidates": promoted,
        "p4s_decision": "GO" if promoted else "STOP",
        "selected_candidate": promoted[0] if promoted else None,
        "candidate_application":
            policy["p4_s_markov_gate"]["candidate_application"],
        "p4c_official_decision_note":
            "P4-C official STOP is preserved; this artifact records P4-S "
            "functional Markov evidence for the supplied candidate logs.",
    }
    decision.update(common_metadata(args, policy_hash))
    write_json(os.path.join(output, "p4s_markov_statistics.json"),
               statistics)
    write_json(os.path.join(output, "p4s_candidate_decision.json"),
               decision)
    count = write_checksums(output)
    print("P4-S Markov evidence complete: {} ({} ledger files)".format(
        decision["p4s_decision"], count))


def add_common_arguments(parser):
    p4c.add_common_arguments(parser)


def build_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    subparsers.required = True
    analyze = subparsers.add_parser("analyze")
    add_common_arguments(analyze)
    analyze.add_argument("--policy", required=True)
    analyze.add_argument("--output", required=True)
    analyze.add_argument("--markov-log", nargs="+", required=True)
    analyze.add_argument("--allow-smoke-policy", action="store_true")
    analyze.set_defaults(func=analyze_command)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
