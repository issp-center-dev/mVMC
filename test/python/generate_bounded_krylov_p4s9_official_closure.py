from __future__ import print_function

import argparse
import copy
import os
import re
import shutil

import generate_bounded_krylov_p4c_evidence as p4c
import generate_bounded_krylov_p4s_evidence as p4s
import generate_bounded_krylov_p4s9_official_statistics_evidence as official


def verify_ledger(root):
    root = os.path.realpath(root)
    ledger = os.path.join(root, "checksums.sha256")
    if not os.path.isfile(ledger):
        raise AssertionError("missing P4-S9 evidence ledger")
    entries = []
    by_path = {}
    with open(ledger, "r") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if "  " not in line:
                raise AssertionError("malformed P4-S9 evidence ledger")
            digest, relative = line.split("  ", 1)
            if len(digest) != 64 or not relative or \
                    os.path.isabs(relative) or \
                    ".." in relative.split(os.sep) or \
                    os.path.normpath(relative) != relative:
                raise AssertionError("unsafe P4-S9 evidence ledger entry")
            lexical_path = os.path.join(root, relative)
            path = os.path.realpath(lexical_path)
            if os.path.commonpath([root, path]) != root or \
                    os.path.islink(lexical_path) or \
                    not os.path.isfile(path) or \
                    p4c.sha256_file(path) != digest:
                raise AssertionError(
                    "P4-S9 evidence ledger mismatch {}".format(relative))
            if relative in by_path:
                raise AssertionError(
                    "duplicate P4-S9 evidence ledger entry {}".format(
                        relative))
            entries.append({"path": relative, "sha256": digest})
            by_path[relative] = digest
    if not entries:
        raise AssertionError("empty P4-S9 evidence ledger")
    actual_files = set()
    for directory, directories, files in os.walk(root):
        directories.sort()
        for dirname in directories:
            if os.path.islink(os.path.join(directory, dirname)):
                raise AssertionError("symlink directory in P4-S9 evidence")
        for filename in sorted(files):
            lexical_path = os.path.join(directory, filename)
            if os.path.abspath(lexical_path) == os.path.abspath(ledger):
                continue
            if os.path.islink(lexical_path):
                raise AssertionError("symlink file in P4-S9 evidence")
            actual_files.add(os.path.relpath(lexical_path, root))
    if set(by_path) != actual_files:
        raise AssertionError("P4-S9 evidence ledger file census mismatch")
    return {
        "ledger_sha256": p4c.sha256_file(ledger),
        "entries": entries,
        "by_path": by_path,
    }


def read_component(root, evidence_name, decision_name, policy_hash,
                   source_commit, baseline_commit, candidate_id,
                   decision_field, pass_value, evidence_pass_field,
                   evidence_artifact, decision_artifact):
    ledger = verify_ledger(root)
    evidence_path = os.path.join(root, evidence_name)
    decision_path = os.path.join(root, decision_name)
    for path in (evidence_path, decision_path):
        relative = os.path.relpath(os.path.realpath(path),
                                   os.path.realpath(root))
        digest = p4c.sha256_file(path)
        if ledger["by_path"].get(relative) != digest:
            raise AssertionError(
                "P4-S9 component is not ledger-bound {}".format(relative))
    evidence = p4c.read_json(evidence_path)
    decision = p4c.read_json(decision_path)
    if evidence.get("official_execution_policy_sha256") != policy_hash or \
            decision.get("official_execution_policy_sha256") != policy_hash or \
            evidence.get("source_commit") != source_commit or \
            decision.get("source_commit") != source_commit or \
            evidence.get("baseline_develop_commit") != baseline_commit or \
            decision.get("baseline_develop_commit") != baseline_commit or \
            evidence.get("candidate_id") != candidate_id or \
            decision.get("candidate_id") != candidate_id or \
            int(evidence.get("schema_version", 0)) != 1 or \
            int(decision.get("schema_version", 0)) != 1 or \
            evidence.get("artifact") != evidence_artifact or \
            decision.get("artifact") != decision_artifact or \
            bool(decision.get("p4f_authorized", True)):
        raise AssertionError("P4-S9 component provenance mismatch")
    passed = decision.get(decision_field) == pass_value
    if bool(evidence.get(evidence_pass_field)) != passed:
        raise AssertionError("P4-S9 component decision mismatch")
    return {
        "root": os.path.abspath(root),
        "evidence_path": evidence_path,
        "decision_path": decision_path,
        "evidence": evidence,
        "decision": decision,
        "passed": passed,
        "ledger": ledger,
    }


def selected_math_trace(path):
    prefixes = ("SCALE ", "SAMPLE ", "SUMMARY ", "TAU ", "ENTRY ",
                "MANIFEST ", "DECISION ")
    with open(path, "r") as handle:
        lines = [line.rstrip("\n") for line in handle
                 if line.startswith(prefixes)]
    return "\n".join(lines) + "\n"


def validate_parity_run(path, policy):
    run = official.parse_session_log(path)
    parity = policy["mpi_parity"]
    header = run["header"]
    rank = int(header.get("rank_count", 0))
    if int(header.get("schema", 0)) != official.DRIVER_SCHEMA_VERSION or \
            header.get("fixture") != official.EXPECTED_FIXTURE or \
            rank not in parity["rank_counts"] or \
            int(header.get("site_count", 0)) != int(parity["site_count"]) or \
            int(header.get("qp_total", 0)) != int(parity["qp_total"]) or \
            int(header.get("sample_count", 0)) != \
            int(parity["sample_count"]) or \
            int(header.get("cache_requested_bytes", 0)) != \
            int(parity["cache_bytes"]) or \
            float(header.get("rho", 0.0)) != float(parity["rho"]) or \
            str(header.get("seed_hex", "")).lower() != \
            str(parity["seed_hex"]).lower() or \
            int(header.get("seed", -1)) != int(parity["seed_hex"], 16) or \
            int(header.get("rng_stream", -1)) != int(parity["rng_stream"]) or \
            int(header.get("persistent_session", 0)) != 1 or \
            int(header.get("amplitude_generation_hash", 0)) == 0:
        raise AssertionError("P4-S9 MPI parity header mismatch")
    statistics = copy.deepcopy(policy["official_statistics"])
    statistics["cache_bytes"] = int(parity["cache_bytes"])
    compatible = copy.deepcopy(run)
    compatible["header"]["schema"] = p4s.SCHEMA_VERSION
    compatible["header"]["fixture"] = p4s.EXPECTED_MARKOV_FIXTURE
    p4s.validate_markov_run(
        compatible,
        official.compatibility_policy(statistics, parity["seed_hex"]),
        True)
    official.validate_gate_consistency(run, statistics)
    session = run["session"]
    if int(session.get("schema", 0)) != official.DRIVER_SCHEMA_VERSION or \
            int(session.get("amplitude_generation_hash", 0)) != \
            int(header["amplitude_generation_hash"]) or \
            int(session.get("session_root_evaluations", 0)) != \
            int(parity["sample_count"]) + 1 or \
            int(session.get("cache_reset_count", 0)) != 1 or \
            int(session.get("engine_heap_allocations", -1)) != 0 or \
            int(session.get("session_end_pass", 0)) != 1:
        raise AssertionError("P4-S9 MPI parity session mismatch")
    trace = selected_math_trace(path)
    return {
        "rank_count": rank,
        "path": os.path.abspath(path),
        "raw_log_sha256": p4c.sha256_file(path),
        "math_trace_sha256": p4c.sha256_text(trace),
        "math_trace": trace,
    }


def validate_focused_ctest(path, policy):
    with open(path, "r") as handle:
        text = handle.read()
    records = {}
    duplicates = []
    result_line = re.compile(
        r"^[ \t]*\d+/\d+[ \t]+Test[ \t]+#\d+:[ \t]+(\S+)[ \t]+",
        re.MULTILINE)
    for match in result_line.finditer(text):
        line_end = text.find("\n", match.start())
        line = text[match.start():line_end if line_end >= 0 else len(text)]
        name = match.group(1)
        if name in records:
            duplicates.append(name)
        records.setdefault(name, []).append({
            "passed": re.search(r"\bPassed\b", line) is not None and
            re.search(r"(?:\*\*\*Failed|\bFailed\b|\bNot Run\b)", line)
            is None,
            "line": line.strip(),
        })
    required = policy["focused_tests"]["required_names"]
    missing = [name for name in required if name not in records]
    failed = [name for name in required
              if name in records and
              (len(records[name]) != 1 or not records[name][0]["passed"])]
    summary = policy["focused_tests"]["pass_summary"]
    failed_results = sorted(
        name for name, values in records.items()
        if any(not value["passed"] for value in values))
    summary_match = re.search(
        r"100% tests passed, 0 tests failed out of (\d+)", text)
    parsed_result_count = sum(len(values) for values in records.values())
    summary_count_matches = (
        summary_match is not None and
        int(summary_match.group(1)) == parsed_result_count)
    passed = (not missing and not failed and not duplicates and
              not failed_results and summary in text and
              summary_count_matches)
    return {
        "path": os.path.abspath(path),
        "sha256": p4c.sha256_file(path),
        "required_test_count": len(required),
        "missing_tests": missing,
        "failed_tests": failed,
        "all_failed_test_results": failed_results,
        "duplicate_test_results": sorted(set(duplicates)),
        "pass_summary_present": summary in text,
        "parsed_test_result_count": parsed_result_count,
        "summary_count_matches": summary_count_matches,
        "focused_ctest_pass": passed,
    }


def analyze_command(args):
    for label, commit in (("source", args.source_commit),
                          ("baseline", args.baseline_commit)):
        if re.fullmatch(r"[0-9a-f]{40}", commit or "") is None:
            raise AssertionError(
                "P4-S9 {} commit provenance mismatch".format(label))
    policy = official.validate_policy(p4c.read_json(args.policy))
    policy_hash = p4c.sha256_file(args.policy)
    statistics = read_component(
        args.statistics_evidence,
        "p4s9_official_markov_statistics.json",
        "p4s9_official_statistics_decision.json",
        policy_hash, args.source_commit, args.baseline_commit,
        policy["candidate_id"],
        "official_statistics_decision", "GO",
        "official_statistics_pass",
        "p4s9_official_markov_statistics",
        "p4s9_official_statistics_decision")
    resource = read_component(
        args.resource_evidence,
        "p4s9_official_resource_evidence.json",
        "p4s9_official_resource_decision.json",
        policy_hash, args.source_commit, args.baseline_commit,
        policy["candidate_id"],
        "official_resource_decision", "PASS",
        "official_resource_pass",
        "p4s9_official_resource_evidence",
        "p4s9_official_resource_decision")
    parity_runs = [validate_parity_run(path, policy)
                   for path in args.mpi_parity_log]
    rank_counts = [item["rank_count"] for item in parity_runs]
    if sorted(rank_counts) != policy["mpi_parity"]["rank_counts"] or \
            len(set(rank_counts)) != len(rank_counts):
        raise AssertionError("P4-S9 MPI parity rank census mismatch")
    trace_hashes = {item["math_trace_sha256"] for item in parity_runs}
    parity_pass = len(trace_hashes) == 1
    focused = validate_focused_ctest(args.focused_ctest_log, policy)
    all_pass = (statistics["passed"] and resource["passed"] and
                parity_pass and focused["focused_ctest_pass"])

    output = os.path.abspath(args.output)
    frozen = os.path.join(output, "frozen_inputs")
    raw = os.path.join(output, "raw_mpi_parity")
    workflow = os.path.join(output, "workflow")
    for directory in (output, frozen, raw, workflow):
        if not os.path.isdir(directory):
            os.makedirs(directory)
    shutil.copyfile(args.policy,
                    os.path.join(frozen, os.path.basename(args.policy)))
    component_bindings = {}
    for label, component in (("statistics", statistics),
                             ("resource", resource)):
        for kind, path in (("evidence", component["evidence_path"]),
                           ("decision", component["decision_path"]),
                           ("ledger", os.path.join(
                               component["root"], "checksums.sha256"))):
            destination = os.path.join(
                frozen, "{}-{}-{}".format(label, kind,
                                           os.path.basename(path)))
            shutil.copyfile(path, destination)
            component_bindings["{}_{}_sha256".format(label, kind)] = \
                p4c.sha256_file(path)
    parity_summary = []
    for item in parity_runs:
        destination = os.path.join(
            raw, "rank{}-{}".format(
                item["rank_count"], os.path.basename(item["path"])))
        shutil.copyfile(item["path"], destination)
        parity_summary.append({
            "rank_count": item["rank_count"],
            "raw_log": os.path.relpath(destination, output),
            "raw_log_sha256": item["raw_log_sha256"],
            "math_trace_sha256": item["math_trace_sha256"],
        })
    ctest_destination = os.path.join(
        workflow, os.path.basename(args.focused_ctest_log))
    shutil.copyfile(args.focused_ctest_log, ctest_destination)
    focused["raw_log"] = os.path.relpath(ctest_destination, output)
    focused.pop("path")
    generated_at = p4c.generated_at_jst()
    evidence = {
        "schema_version": 1,
        "artifact": "p4s9_official_closure_evidence",
        "generated_at_jst": generated_at,
        "source_commit": args.source_commit,
        "baseline_develop_commit": args.baseline_commit,
        "official_execution_policy_sha256": policy_hash,
        "candidate_id": policy["candidate_id"],
        "component_sha256": component_bindings,
        "statistics_pass": statistics["passed"],
        "resource_pass": resource["passed"],
        "mpi_parity_runs": parity_summary,
        "mpi_parity_trace_sha256": next(iter(trace_hashes))
            if parity_pass else None,
        "mpi_parity_pass": parity_pass,
        "focused_ctest": focused,
        "official_closure_pass": all_pass,
    }
    decision = {
        "schema_version": 1,
        "artifact": "p4s9_official_closure_decision",
        "generated_at_jst": generated_at,
        "source_commit": args.source_commit,
        "baseline_develop_commit": args.baseline_commit,
        "official_execution_policy_sha256": policy_hash,
        "candidate_id": policy["candidate_id"],
        "p4s9_decision": "GO" if all_pass else "STOP",
        "p4f_authorized": all_pass,
        "push_or_github_mutation_authorized": False,
        "failure_action": "begin_local_p4f"
            if all_pass else "preserve_evidence_and_stop_before_p4f",
    }
    p4c.write_json(os.path.join(
        output, "p4s9_official_closure_evidence.json"), evidence)
    p4c.write_json(os.path.join(
        output, "p4s9_official_closure_decision.json"), decision)
    count = p4c.write_checksums(output)
    print("P4-S9 official closure: {}; {} ledger files".format(
        decision["p4s9_decision"], count))


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--policy", required=True)
    parser.add_argument("--statistics-evidence", required=True)
    parser.add_argument("--resource-evidence", required=True)
    parser.add_argument("--mpi-parity-log", nargs="+", required=True)
    parser.add_argument("--focused-ctest-log", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--baseline-commit", required=True)
    return parser


def main(argv=None):
    analyze_command(build_parser().parse_args(argv))


if __name__ == "__main__":
    main()
