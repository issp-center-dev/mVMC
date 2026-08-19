#!/usr/bin/env python3

import argparse
import datetime as dt
import json
import math
import os
import pathlib
import subprocess
from zoneinfo import ZoneInfo

import generate_power_lanczos_p4f_gevp_evidence as base


SOLVER_POLICY_SHA256 = (
    "756310504c90825630a0a86dfbf60903024cc1520c73c4d84eb8554e51b03323"
)
CONFIRMATORY_POLICY_SHA256 = (
    "c1348fee957b4eaefa24025fdfb7c79277c800d377f6b7e08ffbfc21aa890f9f"
)


def parse_raw_trace(trace_bytes):
    lines = trace_bytes.decode("ascii").splitlines()
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
    required_decision = (
        "support_pass", "tau_pass", "block_stability_pass", "budget_pass",
        "denominator_pass", "block_pathology_pass",
    )
    p4s_pass = (decision.get("p4s_decision") == "GO" and
                all(decision.get(key) == "1" for key in required_decision))
    return {
        "header": header,
        "session": session,
        "decision": decision,
        "shape_counts": {key.strip(): value for key, value in counts.items()},
        "p4s_pass": p4s_pass,
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
            int(header.get("official_block_length", -1)) == 2048,
        "p4s": raw["p4s_pass"],
    }
    return {
        "seed": seed,
        "site_count": case[0],
        "qp_total": case[1],
        "checks": checks,
        "passed": all(checks.values()),
    }


def dimension_mandatory_pass(dimension, solver_policy):
    exact_anti_limit = solver_policy["hermitianization"][
        "exact_maximum_relative_residual"
    ]
    return all((
        dimension["solver_pass"],
        dimension["cutoff_scan"]["full_pass"],
        dimension["cutoff_scan"]["exact_pass"],
        all(prefix["exact_energy_pass"] for prefix in
            dimension["prefixes"]),
        all(prefix["rank_stability_pass"] for prefix in
            dimension["prefixes"]),
        dimension["coefficient_gate"]["passed"],
        dimension["exact"]["antihermitian_residual"] <= exact_anti_limit,
    ))


def pool_group(cases, dimension_value, solver_policy):
    if len(cases) != 4:
        raise ValueError("pooled group must contain four seeds")
    dimensions = [next(item for item in case["dimensions"]
                       if item["dimension"] == dimension_value)
                  for case in cases]
    per_seed = []
    seed_patterns_pass = True
    for case, dimension in zip(cases, dimensions):
        errors = [prefix["energy_se"] for prefix in dimension["prefixes"]]
        exact_tolerance = (
            solver_policy["exact_gate"]["energy_tolerance_relative"] *
            max(1.0, abs(dimension["exact"]["energy"]))
        )
        if all(value == 0.0 for value in errors):
            zero_pattern = "all_zero"
            zero_pass = all(
                prefix["exact_energy_error"] <= exact_tolerance
                for prefix in dimension["prefixes"]
            )
        elif any(value == 0.0 for value in errors):
            zero_pattern = "mixed_zero"
            zero_pass = False
        else:
            zero_pattern = "all_nonzero"
            zero_pass = True
        seed_patterns_pass &= zero_pass
        per_seed.append({
            "seed": case["confirmatory_identity"]["seed"],
            "energy_se": errors,
            "zero_pattern": zero_pattern,
            "zero_pattern_pass": zero_pass,
            "reported_final_to_initial_ratio":
                errors[-1] / errors[0] if errors[0] > 0.0 else None,
            "reported_theil_sen_slope":
                dimension["se_convergence"]
                ["theil_sen_log_se_log_n_slope"],
        })
    pooled = [
        math.sqrt(sum(item["energy_se"][prefix] ** 2 for item in per_seed) /
                  len(per_seed))
        for prefix in range(3)
    ]
    if all(value == 0.0 for value in pooled):
        zero_rule = "all_zero"
        all_exact = all(
            prefix["exact_energy_error"] <=
            solver_policy["exact_gate"]["energy_tolerance_relative"] *
            max(1.0, abs(dimension["exact"]["energy"]))
            for dimension in dimensions for prefix in dimension["prefixes"]
        )
        ratio = None
        slope = None
        convergence_pass = all_exact
    elif any(value == 0.0 for value in pooled):
        zero_rule = "mixed_zero"
        ratio = None
        slope = None
        convergence_pass = False
    else:
        zero_rule = "all_nonzero"
        ratio = pooled[-1] / pooled[0]
        slope = base.theil_sen_slope((8192, 16384, 32768), pooled)
        convergence_pass = ratio <= 0.8 and slope < 0.0
    return {
        "dimension": dimension_value,
        "seed_count": len(cases),
        "sample_counts": [8192, 16384, 32768],
        "per_seed_diagnostics": per_seed,
        "pooled_energy_se": pooled,
        "zero_rule": zero_rule,
        "final_to_initial_ratio": ratio,
        "maximum_ratio": 0.8,
        "theil_sen_log_se_log_n_slope": slope,
        "seed_zero_patterns_pass": seed_patterns_pass,
        "convergence_pass": convergence_pass,
        "passed": seed_patterns_pass and convergence_pass,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trace-dir", required=True, type=pathlib.Path)
    parser.add_argument("--reparse-driver", required=True, type=pathlib.Path)
    parser.add_argument("--solver-policy", required=True, type=pathlib.Path)
    parser.add_argument("--confirmatory-policy", required=True,
                        type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()

    solver_policy = json.loads(args.solver_policy.read_text(encoding="utf-8"))
    confirmatory_policy = json.loads(
        args.confirmatory_policy.read_text(encoding="utf-8")
    )
    solver_sha = base.sha256_path(args.solver_policy)
    confirmatory_sha = base.sha256_path(args.confirmatory_policy)
    if solver_sha != SOLVER_POLICY_SHA256:
        raise SystemExit("unexpected solver policy SHA-256")
    if CONFIRMATORY_POLICY_SHA256 is None or (
        confirmatory_sha != CONFIRMATORY_POLICY_SHA256
    ):
        raise SystemExit("unexpected confirmatory policy SHA-256")
    if confirmatory_policy["solver_policy_sha256"] != solver_sha:
        raise SystemExit("confirmatory/solver policy binding mismatch")

    traces = sorted(args.trace_dir.glob("*.log"))
    if len(traces) != confirmatory_policy["confirmatory_dataset"][
        "trace_count"
    ]:
        raise SystemExit(f"expected 24 traces, found {len(traces)}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_output_dir = args.output_dir / "raw_reparse"
    raw_output_dir.mkdir(exist_ok=True)
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
            env=environment, check=False, timeout=180,
        )
        if run.returncode != 0:
            raise SystemExit(
                f"reparse failed for {trace.name}: " +
                run.stderr.decode("utf-8", errors="replace")
            )
        output_path = raw_output_dir / f"{trace.name}.reparse.log"
        output_path.write_bytes(run.stdout)
        case = base.evaluate_case(
            trace.name, trace_bytes, run.stdout.decode("utf-8"),
            solver_policy,
        )
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
            pathlib.PurePosixPath("raw_reparse") / output_path.name
        )
        case["reparse_output_sha256"] = base.sha256_path(output_path)
        cases.append(case)
        print(f"[{index + 1:02d}/{len(traces)}] {trace.name}: "
              f"mandatory={case['decision']}", flush=True)

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
                "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
            )
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
        "stopped_traces": stopped_traces,
        "stopped_pooled_gates": stopped_pooled,
        "decision": decision_value,
    }
    evidence_path = args.output_dir / "p4f_v3_confirmatory_evidence.json"
    decision_path = args.output_dir / "p4f_v3_confirmatory_decision.json"
    metadata_path = args.output_dir / "metadata.json"
    solver_copy = args.output_dir / "solver_policy.json"
    confirmatory_copy = args.output_dir / "confirmatory_policy.json"
    evidence_path.write_text(json.dumps(evidence, indent=2,
                                        allow_nan=False) + "\n",
                             encoding="utf-8")
    decision_path.write_text(json.dumps(decision, indent=2,
                                        allow_nan=False) + "\n",
                             encoding="utf-8")
    metadata_path.write_text(json.dumps(metadata, indent=2,
                                        allow_nan=False) + "\n",
                             encoding="utf-8")
    solver_copy.write_bytes(args.solver_policy.read_bytes())
    confirmatory_copy.write_bytes(args.confirmatory_policy.read_bytes())
    checksum_paths = sorted(raw_output_dir.iterdir()) + [
        evidence_path, decision_path, metadata_path, solver_copy,
        confirmatory_copy,
    ]
    (args.output_dir / "checksums.sha256").write_text(
        "".join(
            f"{base.sha256_path(path)}  {path.relative_to(args.output_dir)}\n"
            for path in checksum_paths
        ),
        encoding="utf-8",
    )
    print(json.dumps(decision, sort_keys=True, allow_nan=False))
    return 0 if decision_value == "GO" else 1


if __name__ == "__main__":
    raise SystemExit(main())
