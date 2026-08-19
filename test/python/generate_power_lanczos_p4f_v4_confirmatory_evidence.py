#!/usr/bin/env python3

import argparse
import datetime as dt
import json
import math
import os
import pathlib
import subprocess
import sys
from zoneinfo import ZoneInfo

import generate_power_lanczos_p4f_gevp_evidence as base
import generate_power_lanczos_p4f_v3_confirmatory_evidence as v3


SOLVER_POLICY_SHA256 = (
    "756310504c90825630a0a86dfbf60903024cc1520c73c4d84eb8554e51b03323"
)
CONFIRMATORY_POLICY_SHA256 = (
    "6056c2944483dd1ff087771fb5e8f1c997ce78c54e77ecfa3fa03722605c00d4"
)
POLICY_ID = "power_lanczos_zero_support_p4f_confirmatory_v4"
RATIO_DIAGNOSTIC_MAXIMUM = 1.25
HARD_P4S_FIELDS = (
    "support_pass", "tau_pass", "budget_pass", "denominator_pass",
    "block_pathology_pass",
)


def strict_json_text(value, label):
    def reject_constant(token):
        raise ValueError(f"nonfinite JSON token in {label}: {token}")

    def reject_duplicate(pairs):
        result = {}
        for key, item in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON key in {label}: {key}")
            result[key] = item
        return result

    return json.loads(value, parse_constant=reject_constant,
                      object_pairs_hook=reject_duplicate)


def finite_field(item, key, minimum=None):
    try:
        value = float(item[key])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"invalid P4-S field: {key}") from error
    if not math.isfinite(value) or (minimum is not None and value < minimum):
        raise ValueError(f"invalid P4-S field: {key}")
    return value


def parse_raw_trace(trace_bytes):
    try:
        lines = trace_bytes.decode("ascii").splitlines()
    except UnicodeDecodeError as error:
        raise ValueError("non-ASCII raw trace") from error
    prefixes = (
        "MARKOV ", "SCALE ", "SAMPLE ", "SESSION ", "TRACE ",
        "SUMMARY ", "TAU ", "ENTRY ", "MANIFEST ", "DECISION ",
    )
    counts = {prefix: sum(line.startswith(prefix) for line in lines)
              for prefix in prefixes}
    expected = {
        "MARKOV ": 1, "SCALE ": 4, "SAMPLE ": 32768,
        "SESSION ": 1, "TRACE ": 1, "SUMMARY ": 1, "TAU ": 37,
        "ENTRY ": 18, "MANIFEST ": 1, "DECISION ": 1,
    }
    if counts != expected:
        raise ValueError(f"raw trace shape mismatch: {counts}")
    header = base.fields(next(line for line in lines
                              if line.startswith("MARKOV ")))
    session = base.fields(next(line for line in lines
                               if line.startswith("SESSION ")))
    decision = base.fields(next(line for line in lines
                                if line.startswith("DECISION ")))
    entries = [base.fields(line) for line in lines
               if line.startswith("ENTRY ")]
    observed_keys = []
    diagnostic_pass = True
    pathology_pass = True
    budget_pass = True
    maximum_ratio = 0.0
    for entry in entries:
        try:
            key = (entry["family"], int(entry["row"]), int(entry["column"]))
        except (KeyError, ValueError) as error:
            raise ValueError("invalid P4-S ENTRY identity") from error
        observed_keys.append(key)
        official_se = finite_field(entry, "official_se", 0.0)
        diagnostic_se = finite_field(entry, "diagnostic_se", 0.0)
        ratio = finite_field(entry, "stability_ratio", 1.0)
        recorded_diagnostic = entry.get("stability_pass")
        recorded_pathology = entry.get("pathology_pass")
        recorded_budget = entry.get("budget_pass")
        expected_diagnostic = ratio <= RATIO_DIAGNOSTIC_MAXIMUM
        one_zero = max(official_se, diagnostic_se) > 0.0 and \
            min(official_se, diagnostic_se) == 0.0
        if max(official_se, diagnostic_se) == 0.0:
            expected_ratio = 1.0
        elif one_zero:
            expected_ratio = sys.float_info.max
        else:
            expected_ratio = max(official_se, diagnostic_se) / \
                min(official_se, diagnostic_se)
        if not math.isclose(ratio, expected_ratio, rel_tol=2.0e-14,
                            abs_tol=1.0e-15):
            raise ValueError("P4-S block ratio value mismatch")
        if recorded_diagnostic not in ("0", "1") or \
                (recorded_diagnostic == "1") is not expected_diagnostic:
            raise ValueError("P4-S block ratio diagnostic mismatch")
        if recorded_pathology not in ("0", "1") or \
                (recorded_pathology == "1") is not (not one_zero):
            raise ValueError("P4-S block pathology mismatch")
        if recorded_budget not in ("0", "1"):
            raise ValueError("P4-S conservative budget flag mismatch")
        diagnostic_pass &= expected_diagnostic
        pathology_pass &= not one_zero
        budget_pass &= recorded_budget == "1"
        maximum_ratio = max(maximum_ratio, ratio)
    expected_keys = {
        (family, row, column)
        for family in ("S", "K", "B")
        for row, column in ((0, 0), (0, 1), (0, 2),
                            (1, 1), (1, 2), (2, 2))
    }
    if len(observed_keys) != len(set(observed_keys)) or \
            set(observed_keys) != expected_keys:
        raise ValueError("P4-S ENTRY census mismatch")
    if any(decision.get(key) not in ("0", "1")
           for key in HARD_P4S_FIELDS + ("block_stability_pass",)):
        raise ValueError("invalid P4-S DECISION flags")
    hard_pass = all(decision[key] == "1" for key in HARD_P4S_FIELDS)
    if decision.get("p4s_decision") not in ("GO", "STOP") or \
            (decision["p4s_decision"] == "GO") is not hard_pass or \
            (decision["budget_pass"] == "1") is not budget_pass or \
            (decision["block_pathology_pass"] == "1") is not \
            pathology_pass or \
            (decision["block_stability_pass"] == "1") is not \
            diagnostic_pass or \
            not math.isclose(
                finite_field(decision, "maximum_block_stability_ratio", 1.0),
                maximum_ratio, rel_tol=2.0e-14, abs_tol=1.0e-15):
        raise ValueError("P4-S DECISION aggregate mismatch")
    return {
        "header": header,
        "session": session,
        "decision": decision,
        "shape_counts": {key.strip(): value for key, value in counts.items()},
        "p4s_hard_pass": hard_pass,
        "block_stability_diagnostic_pass": diagnostic_pass,
        "maximum_block_stability_ratio": maximum_ratio,
    }


def validate_trace_identity(raw, policy):
    header = raw["header"]
    session = raw["session"]
    dataset = policy["confirmatory_dataset"]
    seed = header.get("seed_hex", "").lower()
    expected_seeds = {value.lower() for value in dataset["seed_order"]}
    case = (int(header.get("site_count", -1)),
            int(header.get("qp_total", -1)))
    expected_cases = {(item["site_count"], item["qp_total"])
                      for item in dataset["case_order"]}
    checks = {
        "driver_schema": int(header.get("schema", -1)) ==
            dataset["driver_schema"],
        "session_schema": int(session.get("schema", -1)) ==
            dataset["driver_schema"],
        "fixture": header.get("fixture") ==
            "p4s9_long_direct_session_official",
        "seed": seed in expected_seeds,
        "case": case in expected_cases,
        "sample_count": int(header.get("sample_count", -1)) ==
            dataset["sample_count_per_trace"],
        "rank_count": int(header.get("rank_count", -1)) == 1,
        "cache": int(header.get("cache_requested_bytes", -1)) ==
            dataset["cache_requested_bytes"],
        "rho": math.isclose(float(header.get("rho", "nan")),
                            dataset["rho"], rel_tol=0.0, abs_tol=0.0),
        "rng_stream": int(header.get("rng_stream", -1)) ==
            dataset["rng_stream"],
        "persistent_session": int(header.get("persistent_session", -1)) == 1,
        "proposal": int(header.get("global_numerator", -1)) == 0 and
            int(header.get("global_denominator", -1)) == 1,
        "partition": int(header.get("official_block_count", -1)) == 16 and
            int(header.get("official_block_length", -1)) == 2048 and
            int(header.get("diagnostic_block_count", -1)) == 32 and
            int(header.get("diagnostic_block_length", -1)) == 1024,
        "p4s": raw["p4s_hard_pass"],
    }
    return {
        "seed": seed,
        "site_count": case[0],
        "qp_total": case[1],
        "checks": checks,
        "passed": all(checks.values()),
        "block_stability_diagnostic": {
            "role": "report_only_not_a_hard_gate",
            "diagnostic_maximum": RATIO_DIAGNOSTIC_MAXIMUM,
            "passed": raw["block_stability_diagnostic_pass"],
            "maximum_ratio": raw["maximum_block_stability_ratio"],
        },
    }


dimension_mandatory_pass = v3.dimension_mandatory_pass
pool_group = v3.pool_group


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trace-dir", required=True, type=pathlib.Path)
    parser.add_argument("--reparse-driver", required=True, type=pathlib.Path)
    parser.add_argument("--solver-policy", required=True, type=pathlib.Path)
    parser.add_argument("--confirmatory-policy", required=True,
                        type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()

    solver_policy = strict_json_text(
        args.solver_policy.read_text(encoding="utf-8"), "solver policy")
    confirmatory_policy = strict_json_text(
        args.confirmatory_policy.read_text(encoding="utf-8"),
        "confirmatory policy")
    solver_sha = base.sha256_path(args.solver_policy)
    confirmatory_sha = base.sha256_path(args.confirmatory_policy)
    if solver_sha != SOLVER_POLICY_SHA256:
        raise SystemExit("unexpected solver policy SHA-256")
    if confirmatory_sha != CONFIRMATORY_POLICY_SHA256:
        raise SystemExit("unexpected confirmatory policy SHA-256")
    if confirmatory_policy.get("solver_policy_sha256") != solver_sha or \
            confirmatory_policy.get("policy_id") != POLICY_ID:
        raise SystemExit("confirmatory/solver policy binding mismatch")

    traces = sorted(args.trace_dir.glob("*.log"))
    if len(traces) != confirmatory_policy["confirmatory_dataset"][
            "trace_count"]:
        raise SystemExit(f"expected 24 traces, found {len(traces)}")
    if args.output_dir.exists() or args.output_dir.is_symlink():
        raise SystemExit("confirmatory output already exists")
    args.output_dir.mkdir(parents=False)
    raw_output_dir = args.output_dir / "raw_reparse"
    raw_output_dir.mkdir()
    environment = os.environ.copy()
    environment.update({
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
    })
    cases = []
    observed_pairs = set()
    for index, trace in enumerate(traces):
        trace_bytes = trace.read_bytes()
        raw = parse_raw_trace(trace_bytes)
        identity = validate_trace_identity(raw, confirmatory_policy)
        pair = (identity["seed"], identity["site_count"],
                identity["qp_total"])
        if pair in observed_pairs:
            raise SystemExit(f"duplicate confirmatory trace identity: {pair}")
        observed_pairs.add(pair)
        run = subprocess.run(
            [str(args.reparse_driver), "-"], input=trace_bytes,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env=environment, check=False, timeout=180)
        if run.returncode != 0:
            raise SystemExit(
                f"reparse failed for {trace.name}: " +
                run.stderr.decode("utf-8", errors="replace"))
        output_path = raw_output_dir / f"{trace.name}.reparse.log"
        output_path.write_bytes(run.stdout)
        case = base.evaluate_case(
            trace.name, trace_bytes, run.stdout.decode("utf-8"),
            solver_policy)
        case["v2_per_trace_se_decision_diagnostic"] = case["decision"]
        case["confirmatory_identity"] = identity
        mandatory = (
            identity["passed"] and case["reconstruction"]["passed"] and
            all(item["passed"] for item in case["matrix_diagnostics"]) and
            all(dimension_mandatory_pass(item, solver_policy)
                for item in case["dimensions"])
        )
        case["mandatory_gate_pass"] = mandatory
        case["decision"] = "GO" if mandatory else "STOP"
        case["reparse_output"] = str(
            pathlib.PurePosixPath("raw_reparse") / output_path.name)
        case["reparse_output_sha256"] = base.sha256_path(output_path)
        cases.append(case)
        print(f"[{index + 1:02d}/{len(traces)}] {trace.name}: "
              f"mandatory={case['decision']} "
              f"block_ratio_diagnostic="
              f"{'PASS' if identity['block_stability_diagnostic']['passed'] else 'FLAG'}",
              flush=True)

    expected_pairs = {
        (seed.lower(), item["site_count"], item["qp_total"])
        for seed in confirmatory_policy["confirmatory_dataset"]["seed_order"]
        for item in confirmatory_policy["confirmatory_dataset"]["case_order"]
    }
    identity_grid_pass = observed_pairs == expected_pairs
    pooled_gates = []
    for item in confirmatory_policy["confirmatory_dataset"]["case_order"]:
        group = [case for case in cases
                 if case["site_count"] == item["site_count"] and
                 case["qp_total"] == item["qp_total"]]
        for dimension in (2, 3):
            pooled = pool_group(group, dimension, solver_policy)
            pooled["site_count"] = item["site_count"]
            pooled["qp_total"] = item["qp_total"]
            pooled_gates.append(pooled)

    stopped_traces = [case["trace"] for case in cases
                      if not case["mandatory_gate_pass"]]
    stopped_pooled = [
        {"site_count": item["site_count"], "qp_total": item["qp_total"],
         "dimension": item["dimension"]}
        for item in pooled_gates if not item["passed"]
    ]
    flagged_block_ratio_traces = [
        case["trace"] for case in cases
        if not case["confirmatory_identity"]
        ["block_stability_diagnostic"]["passed"]
    ]
    decision_value = "GO" if (
        identity_grid_pass and not stopped_traces and not stopped_pooled
    ) else "STOP"
    now = dt.datetime.now(ZoneInfo("Asia/Tokyo"))
    metadata = {
        "schema": 1,
        "generated_at": now.strftime("%Y-%m-%d %H:%M JST"),
        "model": "OpenAI GPT-5.6 Codex",
        "solver_policy_sha256": solver_sha,
        "confirmatory_policy_sha256": confirmatory_sha,
        "reparse_driver": args.reparse_driver.name,
        "reparse_driver_sha256": base.sha256_path(args.reparse_driver),
        "generator": pathlib.Path(__file__).name,
        "generator_sha256": base.sha256_path(pathlib.Path(__file__)),
        "thread_environment": {
            key: environment[key] for key in (
                "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
        },
    }
    evidence = {
        "schema": 1,
        "metadata": metadata,
        "policy_id": confirmatory_policy["policy_id"],
        "trace_count": len(cases),
        "identity_grid_pass": identity_grid_pass,
        "cases": cases,
        "pooled_se_convergence_gates": pooled_gates,
        "block_stability_diagnostic": {
            "role": "report_only_not_a_hard_gate",
            "diagnostic_maximum": RATIO_DIAGNOSTIC_MAXIMUM,
            "flagged_trace_count": len(flagged_block_ratio_traces),
            "flagged_traces": flagged_block_ratio_traces,
        },
        "stopped_traces": stopped_traces,
        "stopped_pooled_gates": stopped_pooled,
        "decision": decision_value,
    }
    decision = {
        "schema": 1,
        "policy_id": confirmatory_policy["policy_id"],
        "solver_policy_sha256": solver_sha,
        "confirmatory_policy_sha256": confirmatory_sha,
        "trace_count": len(cases),
        "mandatory_trace_go_count": len(cases) - len(stopped_traces),
        "mandatory_trace_stop_count": len(stopped_traces),
        "pooled_gate_go_count": len(pooled_gates) - len(stopped_pooled),
        "pooled_gate_stop_count": len(stopped_pooled),
        "block_stability_diagnostic_flag_count":
            len(flagged_block_ratio_traces),
        "flagged_block_stability_diagnostic_traces":
            flagged_block_ratio_traces,
        "stopped_traces": stopped_traces,
        "stopped_pooled_gates": stopped_pooled,
        "decision": decision_value,
    }
    evidence_path = args.output_dir / "p4f_v4_confirmatory_evidence.json"
    decision_path = args.output_dir / "p4f_v4_confirmatory_decision.json"
    metadata_path = args.output_dir / "metadata.json"
    solver_copy = args.output_dir / "solver_policy.json"
    confirmatory_copy = args.output_dir / "confirmatory_policy.json"
    evidence_path.write_text(
        json.dumps(evidence, indent=2, allow_nan=False) + "\n",
        encoding="utf-8")
    decision_path.write_text(
        json.dumps(decision, indent=2, allow_nan=False) + "\n",
        encoding="utf-8")
    metadata_path.write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n",
        encoding="utf-8")
    solver_copy.write_bytes(args.solver_policy.read_bytes())
    confirmatory_copy.write_bytes(args.confirmatory_policy.read_bytes())
    checksum_paths = sorted(raw_output_dir.iterdir()) + [
        evidence_path, decision_path, metadata_path, solver_copy,
        confirmatory_copy,
    ]
    (args.output_dir / "checksums.sha256").write_text(
        "".join(
            f"{base.sha256_path(path)}  "
            f"{path.relative_to(args.output_dir)}\n"
            for path in checksum_paths),
        encoding="utf-8")
    print(json.dumps(decision, sort_keys=True, allow_nan=False))
    return 0 if decision_value == "GO" else 1


if __name__ == "__main__":
    raise SystemExit(main())
