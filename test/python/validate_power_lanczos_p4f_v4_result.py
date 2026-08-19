#!/usr/bin/env python3

import argparse
import json
import pathlib
import re

import generate_power_lanczos_p4f_v4_confirmatory_evidence as generator
import validate_power_lanczos_p4f_v3_result as v3


POLICY_ID = "power_lanczos_zero_support_p4f_confirmatory_v4"
EXPECTED_KEYS = {
    "schema", "solver_policy_sha256", "confirmatory_policy_sha256",
    "package_manifest_sha256", "job_script_sha256", "generator_sha256",
    "source_archive_sha256", "source_commit",
}
HEX40 = re.compile(r"[0-9a-f]{40}")


def validate_expected(expected):
    if not isinstance(expected, dict) or set(expected) != EXPECTED_KEYS or \
            expected.get("schema") != 1:
        raise AssertionError("invalid v4 expected-provenance schema")
    for key in EXPECTED_KEYS - {"schema", "source_commit"}:
        if v3.HEX64.fullmatch(str(expected.get(key, ""))) is None:
            raise AssertionError(f"invalid v4 expected digest: {key}")
    if HEX40.fullmatch(str(expected.get("source_commit", ""))) is None:
        raise AssertionError("invalid v4 expected source commit")
    return dict(expected)


def verify_package_provenance(view, expected):
    package_manifest_path = "workflow/package_manifest.sha256"
    job_path = "workflow/p4f_v4_genkai_job.sh"
    source_name = "mVMC-p4f-v4-source.tar.gz"
    solver_frozen = \
        "frozen_inputs/power_lanczos_zero_support_p4f_solver_policy.json"
    confirmatory_frozen = \
        "frozen_inputs/power_lanczos_zero_support_p4f_v4_confirmatory_policy.json"
    if view.digests.get(package_manifest_path) != \
            expected["package_manifest_sha256"] or \
            view.digests.get(job_path) != expected["job_script_sha256"]:
        raise AssertionError("v4 package workflow provenance mismatch")
    package_entries = v3.parse_checksum_lines(
        view.read_bytes(package_manifest_path), "v4 package manifest")
    required = {
        solver_frozen: expected["solver_policy_sha256"],
        confirmatory_frozen: expected["confirmatory_policy_sha256"],
        source_name: expected["source_archive_sha256"],
        job_path: expected["job_script_sha256"],
    }
    for path, digest in required.items():
        if package_entries.get(path) != digest:
            raise AssertionError(f"v4 package provenance mismatch: {path}")
    source_checksum = v3.parse_checksum_lines(
        view.read_bytes("workflow/source_archive.sha256"),
        "v4 source archive checksum")
    if source_checksum != {source_name: expected["source_archive_sha256"]}:
        raise AssertionError("v4 source archive provenance mismatch")
    build_config = view.read_bytes("workflow/build_config.txt").decode(
        "utf-8", errors="strict")
    required_lines = (
        f"source_commit={expected['source_commit']}",
        f"baseline_develop_commit={expected['source_commit']}",
        f"source_archive_sha256={expected['source_archive_sha256']}",
        f"solver_policy_sha256={expected['solver_policy_sha256']}",
        f"confirmatory_policy_sha256={expected['confirmatory_policy_sha256']}",
        "PJM_MPI_PROC=4",
        "strict_fp=icc-fp-model-precise-no-fast-math",
    )
    lines = build_config.splitlines()
    for required_line in required_lines:
        key = required_line.split("=", 1)[0]
        if [line for line in lines if line.startswith(f"{key}=")] != \
                [required_line]:
            raise AssertionError(f"v4 build provenance mismatch: {key}")
    ctest = view.read_bytes("workflow/genkai_focused_ctest.log").decode(
        "utf-8", errors="strict")
    if "100% tests passed" not in ctest or \
            re.search(r"0 tests failed out of 6\b", ctest) is None:
        raise AssertionError("v4 focused CTest result mismatch")
    inventory = view.read_bytes(
        "workflow/remote_scratch_inventory.txt").decode("ascii")
    if re.search(r"^raw_trace_count=24$", inventory,
                 flags=re.MULTILINE) is None:
        raise AssertionError("v4 remote trace census mismatch")


def validate_block_diagnostic(identity, raw, trace):
    recorded = identity.get("block_stability_diagnostic")
    expected = {
        "role": "report_only_not_a_hard_gate",
        "diagnostic_maximum": generator.RATIO_DIAGNOSTIC_MAXIMUM,
        "passed": raw["block_stability_diagnostic_pass"],
        "maximum_ratio": raw["maximum_block_stability_ratio"],
    }
    if recorded != expected or \
            identity.get("checks", {}).get("p4s") is not \
            raw["p4s_hard_pass"]:
        raise AssertionError(
            f"v4 block diagnostic/hard-gate binding mismatch: {trace}")
    return not raw["block_stability_diagnostic_pass"]


def validate_result(archive_path, checksum_path, expected):
    expected = validate_expected(expected)
    generator_path = pathlib.Path(generator.__file__).resolve()
    if not generator_path.is_file() or generator_path.is_symlink() or \
            v3.sha256_path(generator_path) != expected["generator_sha256"]:
        raise AssertionError("active v4 generator provenance mismatch")
    archive_path = pathlib.Path(archive_path)
    outer_sha = v3.verify_outer_checksum(archive_path, checksum_path)
    view = v3.ArchiveView(archive_path)
    if v3.sha256_path(archive_path) != outer_sha:
        raise AssertionError("v4 result archive changed while validating")
    ledger_path = "artifact-ledger.sha256"
    ledger = v3.parse_checksum_lines(view.read_bytes(ledger_path),
                                     "v4 artifact ledger")
    actual = {path: digest for path, digest in view.digests.items()
              if path != ledger_path}
    if ledger != actual:
        raise AssertionError("v4 artifact ledger census/checksum mismatch")
    nested_path = "confirmatory/checksums.sha256"
    nested = v3.parse_checksum_lines(
        view.read_bytes(nested_path), "v4 confirmatory checksums",
        prefix="confirmatory")
    nested_actual = {path: digest for path, digest in view.digests.items()
                     if path.startswith("confirmatory/") and
                     path != nested_path}
    if nested != nested_actual:
        raise AssertionError("v4 confirmatory checksum census mismatch")

    raw_paths = sorted(path for path in view.digests
                       if path.startswith("raw_traces/") and
                       path.endswith(".log"))
    reparse_paths = sorted(path for path in view.digests
                           if path.startswith("confirmatory/raw_reparse/") and
                           path.endswith(".reparse.log"))
    if len(raw_paths) != 24 or len(reparse_paths) != 24:
        raise AssertionError("v4 raw/reparse trace census mismatch")

    solver_path = "confirmatory/solver_policy.json"
    confirmatory_path = "confirmatory/confirmatory_policy.json"
    if view.digests.get(solver_path) != expected["solver_policy_sha256"] or \
            view.digests.get(confirmatory_path) != \
            expected["confirmatory_policy_sha256"]:
        raise AssertionError("v4 frozen result policy hash mismatch")
    solver_policy = v3.strict_json(view.read_bytes(solver_path), solver_path)
    confirmatory_policy = v3.strict_json(
        view.read_bytes(confirmatory_path), confirmatory_path)
    if confirmatory_policy.get("solver_policy_sha256") != \
            expected["solver_policy_sha256"] or \
            confirmatory_policy.get("policy_id") != POLICY_ID:
        raise AssertionError("v4 frozen policy binding mismatch")
    verify_package_provenance(view, expected)

    evidence_name = "confirmatory/p4f_v4_confirmatory_evidence.json"
    decision_name = "confirmatory/p4f_v4_confirmatory_decision.json"
    evidence = v3.strict_json(view.read_bytes(evidence_name), "v4 evidence")
    decision = v3.strict_json(view.read_bytes(decision_name), "v4 decision")
    metadata = v3.strict_json(
        view.read_bytes("confirmatory/metadata.json"), "v4 metadata")
    if evidence.get("metadata") != metadata or \
            metadata.get("solver_policy_sha256") != \
            expected["solver_policy_sha256"] or \
            metadata.get("confirmatory_policy_sha256") != \
            expected["confirmatory_policy_sha256"] or \
            metadata.get("generator_sha256") != expected["generator_sha256"]:
        raise AssertionError("v4 result metadata provenance mismatch")
    if evidence.get("policy_id") != POLICY_ID or \
            decision.get("policy_id") != POLICY_ID:
        raise AssertionError("v4 result policy ID mismatch")

    dataset = confirmatory_policy.get("confirmatory_dataset", {})
    seeds = [str(value).lower() for value in dataset.get("seed_order", [])]
    case_order = [(item.get("site_count"), item.get("qp_total"))
                  for item in dataset.get("case_order", [])]
    if len(seeds) != 4 or len(set(seeds)) != 4 or \
            len(case_order) != 6 or len(set(case_order)) != 6 or \
            dataset.get("trace_count") != 24:
        raise AssertionError("v4 frozen confirmatory grid mismatch")
    expected_keys = {(seed, site, qp) for seed in seeds
                     for site, qp in case_order}
    cases = evidence.get("cases")
    if not isinstance(cases, list) or len(cases) != 24:
        raise AssertionError("v4 evidence trace census mismatch")
    by_key = {}
    stopped_traces = []
    flagged_traces = []
    for case in cases:
        identity = case.get("confirmatory_identity", {})
        key = (str(identity.get("seed", "")).lower(),
               case.get("site_count"), case.get("qp_total"))
        if identity.get("site_count") != case.get("site_count") or \
                identity.get("qp_total") != case.get("qp_total") or \
                case.get("sample_count") != 32768 or key in by_key:
            raise AssertionError("v4 confirmatory case identity mismatch")
        by_key[key] = case
        trace = case.get("trace")
        if not isinstance(trace, str) or "/" in trace or "\\" in trace:
            raise AssertionError("unsafe v4 evidence trace name")
        raw_path = f"raw_traces/{trace}"
        reparse_relative = case.get("reparse_output")
        if reparse_relative != f"raw_reparse/{trace}.reparse.log":
            raise AssertionError("v4 reparse output binding mismatch")
        reparse_path = f"confirmatory/{reparse_relative}"
        if view.digests.get(raw_path) != case.get("trace_sha256") or \
                view.digests.get(reparse_path) != \
                case.get("reparse_output_sha256"):
            raise AssertionError("v4 raw/reparse SHA binding mismatch")
        raw = generator.parse_raw_trace(view.read_bytes(raw_path))
        if validate_block_diagnostic(identity, raw, trace):
            flagged_traces.append(trace)
        mandatory = v3.trace_mandatory(case, solver_policy)
        if case.get("mandatory_gate_pass") is not mandatory or \
                case.get("decision") != ("GO" if mandatory else "STOP"):
            raise AssertionError("v4 trace mandatory decision mismatch")
        if not mandatory:
            stopped_traces.append(trace)
    if set(by_key) != expected_keys or \
            evidence.get("identity_grid_pass") is not True:
        raise AssertionError("v4 confirmatory identity grid mismatch")

    recorded_pooled = evidence.get("pooled_se_convergence_gates")
    if not isinstance(recorded_pooled, list) or len(recorded_pooled) != 12:
        raise AssertionError("v4 pooled gate census mismatch")
    pooled_by_key = {}
    for item in recorded_pooled:
        key = (item.get("site_count"), item.get("qp_total"),
               item.get("dimension"))
        if key in pooled_by_key:
            raise AssertionError("duplicate v4 pooled gate")
        pooled_by_key[key] = item
    stopped_pooled = []
    for site, qp in case_order:
        group = [by_key[(seed, site, qp)] for seed in seeds]
        for dimension in (2, 3):
            key = (site, qp, dimension)
            if key not in pooled_by_key:
                raise AssertionError("missing v4 pooled gate")
            recomputed = v3.recompute_pooled(group, dimension, solver_policy)
            v3.verify_pooled_record(pooled_by_key[key], recomputed, seeds)
            if not recomputed["passed"]:
                stopped_pooled.append({
                    "site_count": site, "qp_total": qp,
                    "dimension": dimension,
                })
    diagnostic = {
        "role": "report_only_not_a_hard_gate",
        "diagnostic_maximum": generator.RATIO_DIAGNOSTIC_MAXIMUM,
        "flagged_trace_count": len(flagged_traces),
        "flagged_traces": flagged_traces,
    }
    final = "GO" if not stopped_traces and not stopped_pooled else "STOP"
    if evidence.get("trace_count") != 24 or \
            evidence.get("block_stability_diagnostic") != diagnostic or \
            evidence.get("stopped_traces") != stopped_traces or \
            evidence.get("stopped_pooled_gates") != stopped_pooled or \
            evidence.get("decision") != final:
        raise AssertionError("v4 evidence aggregate decision mismatch")
    expected_decision = {
        "schema": 1,
        "policy_id": POLICY_ID,
        "solver_policy_sha256": expected["solver_policy_sha256"],
        "confirmatory_policy_sha256":
            expected["confirmatory_policy_sha256"],
        "trace_count": 24,
        "mandatory_trace_go_count": 24 - len(stopped_traces),
        "mandatory_trace_stop_count": len(stopped_traces),
        "pooled_gate_go_count": 12 - len(stopped_pooled),
        "pooled_gate_stop_count": len(stopped_pooled),
        "block_stability_diagnostic_flag_count": len(flagged_traces),
        "flagged_block_stability_diagnostic_traces": flagged_traces,
        "stopped_traces": stopped_traces,
        "stopped_pooled_gates": stopped_pooled,
        "decision": final,
    }
    if decision != expected_decision:
        raise AssertionError("v4 decision artifact mismatch")
    generator_lines = view.read_bytes(
        "workflow/genkai_confirmatory_generator.log").decode(
            "utf-8", errors="strict").splitlines()
    if not generator_lines or v3.strict_json(
            generator_lines[-1].encode("utf-8"),
            "v4 generator final line") != decision:
        raise AssertionError("v4 generator log decision mismatch")
    return {
        "valid": True,
        "archive_sha256": outer_sha,
        "decision": final,
        "trace_count": 24,
        "mandatory_trace_stop_count": len(stopped_traces),
        "pooled_gate_stop_count": len(stopped_pooled),
        "block_stability_diagnostic_flag_count": len(flagged_traces),
        "artifact_file_count": len(view.digests),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True, type=pathlib.Path)
    parser.add_argument("--checksum", required=True, type=pathlib.Path)
    parser.add_argument("--expected", required=True, type=pathlib.Path)
    args = parser.parse_args()
    expected = v3.strict_json(
        args.expected.read_bytes(), "v4 expected provenance")
    result = validate_result(args.archive, args.checksum, expected)
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
