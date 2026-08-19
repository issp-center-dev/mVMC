#!/usr/bin/env python3

import argparse
import cmath
import datetime as dt
import hashlib
import json
import math
import os
import pathlib
import statistics
import subprocess
import tarfile
from zoneinfo import ZoneInfo


FAMILIES = ("S", "K", "KR", "B")
ENTRIES = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
CUTOFFS = (1.0e-12, 1.0e-10, 1.0e-8)
PREFIX_BLOCKS = (4, 8, 16)
BLOCK_LENGTH = 2048


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_path(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fields(line):
    result = {}
    for token in line.strip().split()[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            result[key] = value
    return result


def close(left, right, multiplier=32768.0):
    scale = max(1.0, abs(left), abs(right))
    return math.isfinite(left) and math.isfinite(right) and abs(left - right) <= (
        multiplier * math.ulp(1.0) * scale
    )


def upper_norm(values):
    total = 0.0
    for entry, value in zip(ENTRIES, values):
        total += (1.0 if entry[0] == entry[1] else 2.0) * abs(value) ** 2
    return math.sqrt(total)


def raw_k_norm(forward, reverse):
    total = 0.0
    for entry, left, right in zip(ENTRIES, forward, reverse):
        total += abs(left) ** 2
        if entry[0] != entry[1]:
            total += abs(right) ** 2
    return math.sqrt(total)


def vector_mean(vectors):
    return [sum(item[index] for item in vectors) / len(vectors)
            for index in range(len(vectors[0]))]


def jackknife_se(values):
    mean = sum(values) / len(values)
    factor = (len(values) - 1.0) / len(values)
    return math.sqrt(factor * sum((value - mean) ** 2 for value in values))


def align(vector, reference):
    overlap = sum(left.conjugate() * right
                  for left, right in zip(reference, vector))
    if abs(overlap) == 0.0:
        return list(vector)
    phase = overlap.conjugate() / abs(overlap)
    return [value * phase for value in vector]


def projective_distance(left, right):
    aligned = align(right, left)
    return math.sqrt(sum(abs(a - b) ** 2 for a, b in zip(left, aligned)))


def theil_sen_slope(sample_counts, errors):
    slopes = []
    for left in range(len(errors)):
        for right in range(left + 1, len(errors)):
            slopes.append((math.log(errors[right]) - math.log(errors[left])) /
                          (math.log(sample_counts[right]) -
                           math.log(sample_counts[left])))
    return statistics.median(slopes)


def sequence_sha256(trace_bytes):
    digest = hashlib.sha256()
    count = 0
    for raw_line in trace_bytes.splitlines():
        if not raw_line.startswith(b"SAMPLE "):
            continue
        item = fields(raw_line.decode("ascii"))
        canonical = (f"{item['sample']}\t{item['configuration']}\t"
                     f"{item['denominator']}\n").encode("ascii")
        digest.update(canonical)
        count += 1
    return digest.hexdigest(), count


def parse_reparse(output):
    summary = None
    decision = None
    entries = {}
    blocks = {}
    gevp = []
    for line in output.splitlines():
        item = fields(line)
        if line.startswith("REPARSE "):
            summary = item
        elif line.startswith("REPARSE_ENTRY "):
            key = (item["family"], int(item["row"]), int(item["column"]))
            entries[key] = {
                "numerator": complex(float(item["numerator_re"]),
                                     float(item["numerator_im"])),
                "denominator": float(item["denominator"]),
            }
        elif line.startswith("REPARSE_BLOCK "):
            block = int(item["block"])
            key = (item["family"], int(item["row"]), int(item["column"]))
            target = blocks.setdefault(block, {"entries": {}})
            denominator = float(item["denominator"])
            sample_count = int(item["sample_count"])
            if "denominator" in target:
                if not close(target["denominator"], denominator) or (
                    target["sample_count"] != sample_count
                ):
                    raise ValueError("inconsistent block denominator/count")
            target["denominator"] = denominator
            target["sample_count"] = sample_count
            target["entries"][key] = complex(
                float(item["numerator_re"]), float(item["numerator_im"])
            )
        elif line.startswith("GEVP "):
            dimension = int(item["dimension"])
            gevp.append({
                "scope": item["scope"],
                "dimension": dimension,
                "cutoff": float(item["cutoff"]),
                "status": item["status"],
                "valid": int(item["valid"]),
                "retained_rank": int(item["retained_rank"]),
                "discarded_rank": int(item["discarded_rank"]),
                "phase_pivot": int(item["phase_pivot"]),
                "root_multiplicity": int(item["root_multiplicity"]),
                "comparison_uses_projector":
                    int(item["comparison_uses_projector"]),
                "energy": float(item["energy"]),
                "energy_squared": float(item["energy_squared"]),
                "variance": float(item["variance"]),
                "residual": float(item["residual"]),
                "antihermitian_residual":
                    float(item["antihermitian_residual"]),
                "condition_estimate": float(item["condition_estimate"]),
                "root_gap": float(item["root_gap"]),
                "coefficient": [
                    complex(float(item[f"coefficient{index}_re"]),
                            float(item[f"coefficient{index}_im"]))
                    for index in range(dimension)
                ],
                "root_subspace_projector": [
                    complex(float(item[f"projector{index}_re"]),
                            float(item[f"projector{index}_im"]))
                    for index in range(dimension * dimension)
                ],
            })
        elif line.startswith("REPARSE_DECISION "):
            decision = item["decision"]
    if summary is None or decision is None:
        raise ValueError("missing reparse summary/decision")
    if len(entries) != len(FAMILIES) * len(ENTRIES):
        raise ValueError("incomplete reconstructed entries")
    if len(blocks) != 16 or any(
        len(block["entries"]) != len(FAMILIES) * len(ENTRIES)
        for block in blocks.values()
    ):
        raise ValueError("incomplete reconstructed blocks")
    return summary, entries, blocks, gevp, decision


def find_gevp(records, scope, dimension, cutoff=1.0e-10):
    candidates = [record for record in records
                  if record["scope"] == scope and
                  record["dimension"] == dimension and
                  math.isclose(record["cutoff"], cutoff,
                               rel_tol=0.0, abs_tol=1.0e-22)]
    if len(candidates) != 1:
        raise ValueError(f"missing/duplicate GEVP {scope} d={dimension} "
                         f"cutoff={cutoff}")
    return candidates[0]


def serializable_gevp(record):
    def finite_or_none(value):
        return value if not isinstance(value, float) or math.isfinite(value) \
            else None

    return {
        key: finite_or_none(value) for key, value in record.items()
        if key not in ("coefficient", "root_subspace_projector")
    } | {
        "coefficient": [
            {"real": value.real, "imag": value.imag}
            for value in record["coefficient"]
        ],
        "root_subspace_projector": [
            {"real": value.real, "imag": value.imag}
            for value in record["root_subspace_projector"]
        ],
    }


def totals_from_entries(entries):
    denominator = entries[("S", 0, 0)]["denominator"]
    matrices = {}
    for family in FAMILIES:
        matrices[family] = []
        for row, column in ENTRIES:
            item = entries[(family, row, column)]
            if not close(item["denominator"], denominator):
                raise ValueError("full entry denominator mismatch")
            matrices[family].append(item["numerator"])
    return denominator, matrices


def block_matrices(block):
    denominator = block["denominator"]
    matrices = {
        family: [block["entries"][(family, row, column)]
                 for row, column in ENTRIES]
        for family in FAMILIES
    }
    return denominator, matrices


def matrix_add(left, right):
    return {family: [a + b for a, b in zip(left[family], right[family])]
            for family in FAMILIES}


def matrix_subtract(left, right):
    return {family: [a - b for a, b in zip(left[family], right[family])]
            for family in FAMILIES}


def matrix_scale(matrix, factor):
    return {family: [value * factor for value in matrix[family]]
            for family in FAMILIES}


def zero_matrices():
    return {family: [0j for _ in ENTRIES] for family in FAMILIES}


def stochastic_matrix_diagnostic(full_denominator, full_matrix,
                                 blocks, family):
    def diagnostic_vector(matrix, denominator):
        normalized = matrix_scale(matrix, 1.0 / denominator)
        if family == "K":
            return [left - right.conjugate()
                    for left, right in zip(normalized["K"],
                                           normalized["KR"])]
        return [complex(0.0, 2.0 * normalized[family][index].imag)
                if row == column else 0j
                for index, (row, column) in enumerate(ENTRIES)]

    full_vector = diagnostic_vector(full_matrix, full_denominator)
    leave_vectors = []
    for block in blocks:
        block_denominator, block_matrix = block_matrices(block)
        leave_vectors.append(diagnostic_vector(
            matrix_subtract(full_matrix, block_matrix),
            full_denominator - block_denominator,
        ))
    mean_leave = vector_mean(leave_vectors)
    covariance_trace = ((len(leave_vectors) - 1.0) / len(leave_vectors)) * sum(
        upper_norm([value - mean for value, mean in zip(vector, mean_leave)]) ** 2
        for vector in leave_vectors
    )
    numerator_norm = upper_norm(full_vector)
    normalized = matrix_scale(full_matrix, 1.0 / full_denominator)
    if family == "K":
        raw_norm = raw_k_norm(normalized["K"], normalized["KR"])
    else:
        raw_norm = upper_norm(normalized[family])
    effect = numerator_norm / raw_norm if raw_norm > 0.0 else numerator_norm
    uncertainty = math.sqrt(max(0.0, covariance_trace))
    numeric_floor = 256.0 * math.ulp(1.0) * raw_norm
    score = numerator_norm / max(uncertainty, numeric_floor)
    passed = ((effect <= 1.0e-10 or score < 4.5) and effect <= 0.25)
    return {
        "family": family,
        "effect_size": effect,
        "jackknife_uncertainty": uncertainty,
        "numeric_floor": numeric_floor,
        "standardized_score": score,
        "passed": passed,
    }


def prefix_result(records, dimension, exact, block_count):
    sample_count = block_count * BLOCK_LENGTH
    full = find_gevp(records, f"prefix{sample_count}", dimension)
    leave = [find_gevp(records,
                       f"prefix{sample_count}_leave{block}", dimension)
             for block in range(block_count)]
    energies = [record["energy"] for record in leave]
    energy_se = jackknife_se(energies)
    energy_gate_limit = max(4.5 * energy_se,
                            1.0e-8 * max(1.0, abs(exact["energy"])))
    aligned = [align(record["coefficient"], full["coefficient"])
               for record in leave]
    coefficient_se = []
    for index in range(dimension):
        coefficient_se.append({
            "real": jackknife_se([value[index].real for value in aligned]),
            "imag": jackknife_se([value[index].imag for value in aligned]),
        })
    projector_se = []
    for index in range(dimension * dimension):
        projector_se.append({
            "real": jackknife_se([
                record["root_subspace_projector"][index].real
                for record in leave
            ]),
            "imag": jackknife_se([
                record["root_subspace_projector"][index].imag
                for record in leave
            ]),
        })
    ranks = [record["retained_rank"] for record in leave]
    multiplicities = [record["root_multiplicity"] for record in leave]
    return {
        "sample_count": sample_count,
        "block_count": block_count,
        "block_length": BLOCK_LENGTH,
        "energy": full["energy"],
        "energy_bias_corrected": (
            block_count * full["energy"] -
            (block_count - 1.0) * statistics.mean(energies)
        ),
        "energy_se": energy_se,
        "exact_energy_error": abs(full["energy"] - exact["energy"]),
        "exact_energy_limit": energy_gate_limit,
        "exact_energy_pass":
            abs(full["energy"] - exact["energy"]) <= energy_gate_limit,
        "full_retained_rank": full["retained_rank"],
        "leave_one_retained_ranks": ranks,
        "full_root_multiplicity": full["root_multiplicity"],
        "leave_one_root_multiplicities": multiplicities,
        "rank_stability_pass":
            all(rank == full["retained_rank"] for rank in ranks) and
            all(value == full["root_multiplicity"]
                for value in multiplicities),
        "coefficient_se": coefficient_se,
        "projector_se": projector_se,
        "full": serializable_gevp(full),
    }


def evaluate_dimension(records, dimension, policy):
    full = find_gevp(records, "full", dimension)
    exact = find_gevp(records, "exact", dimension)
    full_scan = [find_gevp(records, "full", dimension, cutoff)
                 for cutoff in CUTOFFS]
    exact_scan = [find_gevp(records, "exact", dimension, cutoff)
                  for cutoff in CUTOFFS]
    prefixes = [prefix_result(records, dimension, exact, count)
                for count in PREFIX_BLOCKS]
    errors = [item["energy_se"] for item in prefixes]
    exact_tolerance = (policy["exact_gate"]["energy_tolerance_relative"] *
                       max(1.0, abs(exact["energy"])))
    if all(error == 0.0 for error in errors):
        slope = None
        shrinkage_pass = all(item["exact_energy_error"] <= exact_tolerance
                             for item in prefixes)
        zero_rule = "all_zero"
    elif any(error == 0.0 for error in errors):
        slope = None
        shrinkage_pass = False
        zero_rule = "mixed_zero"
    else:
        slope = theil_sen_slope(
            [item["sample_count"] for item in prefixes], errors)
        shrinkage_pass = errors[-1] <= 0.8 * errors[0] and slope < 0.0
        zero_rule = "nonzero"

    final_prefix = prefixes[-1]
    component_checks = []
    uses_projector = (full["root_multiplicity"] > 1 or
                      exact["root_multiplicity"] > 1)
    if uses_projector:
        comparison_kind = "root_subspace_projector"
        for index in range(dimension * dimension):
            for component in ("real", "imag"):
                delta = abs(getattr(
                    full["root_subspace_projector"][index], component) -
                    getattr(exact["root_subspace_projector"][index],
                            component))
                standard_error = final_prefix["projector_se"][index][component]
                limit = max(4.5 * standard_error, 5.0e-8)
                component_checks.append({
                    "index": index,
                    "component": component,
                    "delta": delta,
                    "standard_error": standard_error,
                    "limit": limit,
                    "passed": delta <= limit,
                })
        coefficient_distance = math.sqrt(sum(
            abs(left - right) ** 2 for left, right in zip(
                full["root_subspace_projector"],
                exact["root_subspace_projector"])
        ))
        multiplicity_pass = (
            full["root_multiplicity"] == exact["root_multiplicity"] and
            all(item["rank_stability_pass"] for item in prefixes)
        )
    else:
        comparison_kind = "phase_aligned_coefficient"
        aligned_exact = align(exact["coefficient"], full["coefficient"])
        for index in range(dimension):
            for component in ("real", "imag"):
                delta = abs(getattr(full["coefficient"][index], component) -
                            getattr(aligned_exact[index], component))
                standard_error = final_prefix[
                    "coefficient_se"
                ][index][component]
                limit = max(4.5 * standard_error, 5.0e-8)
                component_checks.append({
                    "index": index,
                    "component": component,
                    "delta": delta,
                    "standard_error": standard_error,
                    "limit": limit,
                    "passed": delta <= limit,
                })
        coefficient_distance = projective_distance(
            full["coefficient"], exact["coefficient"])
        multiplicity_pass = True
    coefficient_pass = (multiplicity_pass and
                        all(item["passed"] for item in component_checks) and
                        coefficient_distance <= 0.1)
    full_span = max(item["energy"] for item in full_scan) - min(
        item["energy"] for item in full_scan)
    full_span_limit = max(
        4.5 * final_prefix["energy_se"],
        1.0e-8 * max(1.0, abs(exact["energy"])),
    )
    exact_span = max(item["energy"] for item in exact_scan) - min(
        item["energy"] for item in exact_scan)
    exact_span_limit = (policy["rank_sensitivity"]
                        ["exact_maximum_energy_span_relative"] *
                        max(1.0, abs(exact["energy"])))
    solver_pass = all(
        item["status"] == "ok" and item["valid"] == 1 and
        item["residual"] <=
        policy["reevaluation"]["maximum_relative_residual"] and
        item["variance"] >= 0.0
        for item in full_scan + exact_scan
    )
    decision = all((
        solver_pass,
        full_span <= full_span_limit,
        exact_span <= exact_span_limit,
        all(item["exact_energy_pass"] and item["rank_stability_pass"]
            for item in prefixes),
        shrinkage_pass,
        coefficient_pass,
        exact["antihermitian_residual"] <=
            policy["hermitianization"]["exact_maximum_relative_residual"],
    ))
    return {
        "dimension": dimension,
        "full": serializable_gevp(full),
        "exact": serializable_gevp(exact),
        "cutoff_scan": {
            "cutoffs": list(CUTOFFS),
            "full_retained_ranks": [item["retained_rank"]
                                    for item in full_scan],
            "full_root_multiplicities": [item["root_multiplicity"]
                                          for item in full_scan],
            "full_energies": [item["energy"] for item in full_scan],
            "full_energy_span": full_span,
            "full_energy_span_limit": full_span_limit,
            "full_pass": full_span <= full_span_limit,
            "exact_retained_ranks": [item["retained_rank"]
                                     for item in exact_scan],
            "exact_root_multiplicities": [item["root_multiplicity"]
                                           for item in exact_scan],
            "exact_energies": [item["energy"] for item in exact_scan],
            "exact_energy_span": exact_span,
            "exact_energy_span_limit": exact_span_limit,
            "exact_pass": exact_span <= exact_span_limit,
        },
        "prefixes": prefixes,
        "se_convergence": {
            "zero_rule": zero_rule,
            "final_to_initial_ratio": (
                errors[-1] / errors[0] if errors[0] > 0.0 else None
            ),
            "theil_sen_log_se_log_n_slope": slope,
            "passed": shrinkage_pass,
        },
        "coefficient_gate": {
            "comparison_kind": comparison_kind,
            "projective_distance": coefficient_distance,
            "maximum_projective_distance": 0.1,
            "multiplicity_pass": multiplicity_pass,
            "components": component_checks,
            "passed": coefficient_pass,
        },
        "solver_pass": solver_pass,
        "decision": "GO" if decision else "STOP",
    }


def evaluate_case(name, trace_bytes, reparse_output, policy):
    summary, entries, block_map, records, reparse_decision = parse_reparse(
        reparse_output
    )
    blocks = [block_map[index] for index in range(16)]
    full_denominator, full_matrix = totals_from_entries(entries)
    summed_denominator = 0.0
    summed_matrix = zero_matrices()
    block_lengths_pass = True
    for block in blocks:
        denominator, matrix = block_matrices(block)
        summed_denominator += denominator
        summed_matrix = matrix_add(summed_matrix, matrix)
        block_lengths_pass &= block["sample_count"] == BLOCK_LENGTH
    block_sum_pass = close(summed_denominator, full_denominator)
    for family in FAMILIES:
        for summed, full in zip(summed_matrix[family], full_matrix[family]):
            block_sum_pass &= close(summed.real, full.real)
            block_sum_pass &= close(summed.imag, full.imag)
    sequence_hash, sequence_count = sequence_sha256(trace_bytes)
    dimensions = [evaluate_dimension(records, dimension, policy)
                  for dimension in (2, 3)]
    matrix_diagnostics = [
        stochastic_matrix_diagnostic(full_denominator, full_matrix,
                                     blocks, family)
        for family in ("K", "S", "B")
    ]
    corrected = float(summary["corrected_antihermitian_residual"])
    solver_corrected = find_gevp(records, "full", 3)[
        "antihermitian_residual"
    ]
    corrected_match = close(corrected, solver_corrected)
    reconstruction_pass = (
        reparse_decision == "GO" and int(summary["all_match"]) == 1 and
        sequence_count == int(summary["sample_count"]) and
        int(summary["sample_count"]) == 32768 and
        int(summary["block_count"]) == 16 and
        int(summary["block_length"]) == BLOCK_LENGTH and
        block_lengths_pass and block_sum_pass and corrected_match
    )
    decision = reconstruction_pass and all(
        item["decision"] == "GO" for item in dimensions
    ) and all(item["passed"] for item in matrix_diagnostics)
    return {
        "trace": name,
        "trace_sha256": sha256_bytes(trace_bytes),
        "configuration_denominator_sequence_sha256": sequence_hash,
        "sequence_canonicalization":
            "ASCII sample<TAB>configuration<TAB>denominator_token<LF>",
        "site_count": int(summary["site_count"]),
        "qp_total": int(summary["qp_total"]),
        "sample_count": int(summary["sample_count"]),
        "sector_size": int(summary["sector_size"]),
        "reconstruction": {
            "sequence_fnv1a64": int(summary["sequence_fnv1a64"]),
            "maximum_denominator_error":
                float(summary["maximum_denominator_error"]),
            "maximum_entry_error": float(summary["maximum_entry_error"]),
            "maximum_exact_error": float(summary["maximum_exact_error"]),
            "source_antihermitian_error":
                float(summary["source_antihermitian_error"]),
            "source_antihermitian_convention":
                summary["source_antihermitian_convention"],
            "legacy_antihermitian_residual":
                float(summary["legacy_antihermitian_residual"]),
            "corrected_antihermitian_residual": corrected,
            "block_sum_pass": block_sum_pass,
            "corrected_solver_match": corrected_match,
            "passed": reconstruction_pass,
        },
        "matrix_diagnostics": matrix_diagnostics,
        "dimensions": dimensions,
        "decision": "GO" if decision else "STOP",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True, type=pathlib.Path)
    parser.add_argument("--reparse-driver", required=True, type=pathlib.Path)
    parser.add_argument("--policy", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()

    policy = json.loads(args.policy.read_text(encoding="utf-8"))
    policy_sha = sha256_path(args.policy)
    if policy["policy_id"] != "power_lanczos_zero_support_p4f_solver_v2":
        raise SystemExit("unexpected P4-F policy ID")
    if policy_sha != "756310504c90825630a0a86dfbf60903024cc1520c73c4d84eb8554e51b03323":
        raise SystemExit("unexpected P4-F policy SHA-256")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = args.output_dir / "raw_reparse"
    raw_dir.mkdir(exist_ok=True)
    cases = []
    environment = os.environ.copy()
    environment.update({
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
    })
    with tarfile.open(args.archive, "r:gz") as archive:
        members = sorted(
            (member for member in archive.getmembers()
             if member.isfile() and
             member.name.startswith("./statistics/raw_markov/") and
             member.name.endswith(".log")),
            key=lambda member: pathlib.PurePosixPath(member.name).name,
        )
        if len(members) != policy["stochastic_gate"]["trace_count"]:
            raise SystemExit(f"expected 24 raw traces, found {len(members)}")
        for index, member in enumerate(members):
            source = archive.extractfile(member)
            if source is None:
                raise SystemExit(f"cannot read {member.name}")
            trace_bytes = source.read()
            run = subprocess.run(
                [str(args.reparse_driver), "-"],
                input=trace_bytes,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=environment,
                check=False,
                timeout=180,
            )
            name = pathlib.PurePosixPath(member.name).name
            output_path = raw_dir / f"{name}.reparse.log"
            output_path.write_bytes(run.stdout)
            if run.returncode != 0:
                raise SystemExit(
                    f"reparse failed for {name}: " +
                    run.stderr.decode("utf-8", errors="replace")
                )
            case = evaluate_case(
                name, trace_bytes, run.stdout.decode("utf-8"), policy
            )
            case["reparse_output"] = str(
                pathlib.PurePosixPath("raw_reparse") / output_path.name
            )
            case["reparse_output_sha256"] = sha256_path(output_path)
            cases.append(case)
            print(f"[{index + 1:02d}/{len(members)}] {name}: "
                  f"{case['decision']}", flush=True)

    stopped = [case["trace"] for case in cases if case["decision"] != "GO"]
    now = dt.datetime.now(ZoneInfo("Asia/Tokyo"))
    evidence = {
        "schema": 1,
        "generated_at": now.strftime("%Y-%m-%d %H:%M JST"),
        "model": "OpenAI GPT-5.6 Codex",
        "scope": "P4-F deterministic reparse and corrected GEVP closure",
        "policy_id": policy["policy_id"],
        "policy_sha256": policy_sha,
        "archive": args.archive.name,
        "archive_sha256": sha256_path(args.archive),
        "reparse_driver": args.reparse_driver.name,
        "reparse_driver_sha256": sha256_path(args.reparse_driver),
        "trace_count": len(cases),
        "cases": cases,
        "aggregate": {
            "go_count": len(cases) - len(stopped),
            "stop_count": len(stopped),
            "stopped_traces": stopped,
            "maximum_denominator_error": max(
                case["reconstruction"]["maximum_denominator_error"]
                for case in cases
            ),
            "maximum_entry_error": max(
                case["reconstruction"]["maximum_entry_error"]
                for case in cases
            ),
            "maximum_corrected_antihermitian_effect": max(
                case["reconstruction"]["corrected_antihermitian_residual"]
                for case in cases
            ),
            "maximum_antihermitian_standardized_score": max(
                next(item["standardized_score"]
                     for item in case["matrix_diagnostics"]
                     if item["family"] == "K")
                for case in cases
            ),
            "maximum_coefficient_projective_distance": max(
                dimension["coefficient_gate"]["projective_distance"]
                for case in cases for dimension in case["dimensions"]
            ),
            "maximum_final_to_initial_energy_se_ratio": max(
                dimension["se_convergence"]["final_to_initial_ratio"] or 0.0
                for case in cases for dimension in case["dimensions"]
            ),
        },
        "decision": "GO" if not stopped else "STOP",
    }
    decision = {
        "schema": 1,
        "policy_id": policy["policy_id"],
        "policy_sha256": policy_sha,
        "trace_count": len(cases),
        "go_count": len(cases) - len(stopped),
        "stop_count": len(stopped),
        "stopped_traces": stopped,
        "decision": evidence["decision"],
    }
    evidence_path = args.output_dir / "p4f_gevp_evidence.json"
    decision_path = args.output_dir / "p4f_gevp_decision.json"
    evidence_path.write_text(json.dumps(evidence, indent=2,
                                        allow_nan=False) + "\n",
                             encoding="utf-8")
    decision_path.write_text(json.dumps(decision, indent=2,
                                        allow_nan=False) + "\n",
                             encoding="utf-8")
    checksums = []
    for path in sorted(raw_dir.iterdir()):
        checksums.append((sha256_path(path), path.relative_to(args.output_dir)))
    checksums.extend((
        (sha256_path(evidence_path), evidence_path.name),
        (sha256_path(decision_path), decision_path.name),
    ))
    (args.output_dir / "checksums.sha256").write_text(
        "".join(f"{digest}  {name}\n" for digest, name in checksums),
        encoding="utf-8",
    )
    print(json.dumps(decision, sort_keys=True, allow_nan=False))
    return 0 if evidence["decision"] == "GO" else 1


if __name__ == "__main__":
    raise SystemExit(main())
