#!/usr/bin/env python3

import argparse
import hashlib
import json
import math
import pathlib
import tempfile

import validate_power_lanczos_p4f_v3_result as validator


FAMILIES = ("S", "K", "B")
ENTRIES = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
ENTRY_KEYS = {(family, row, column) for family in FAMILIES
              for row, column in ENTRIES}
RATIO_DIAGNOSTIC_MAXIMUM = 1.25
BASE_BLOCK_COUNT = 16
BASE_BLOCK_LENGTH = 2048


def fields(line):
    result = {}
    for token in line.split()[1:]:
        if "=" not in token:
            raise AssertionError("malformed evidence field")
        key, value = token.split("=", 1)
        if not key or not value or key in result:
            raise AssertionError("missing/duplicate evidence field")
        result[key] = value
    return result


def finite(value, label, minimum=None):
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise AssertionError(f"invalid numeric {label}") from error
    if not math.isfinite(parsed) or (minimum is not None and
                                     parsed < minimum):
        raise AssertionError(f"invalid numeric {label}")
    return parsed


def close(left, right):
    return math.isclose(left, right, rel_tol=2.0e-12, abs_tol=1.0e-15)


def jackknife(blocks):
    if len(blocks) < 2:
        raise AssertionError("jackknife requires at least two blocks")
    if any(item[0] <= 0 or item[1] <= 0.0 or
           not math.isfinite(item[1]) or
           not math.isfinite(item[2].real) or
           not math.isfinite(item[2].imag) for item in blocks):
        raise AssertionError("invalid jackknife block")
    denominator = math.fsum(item[1] for item in blocks)
    numerator = complex(
        math.fsum(item[2].real for item in blocks),
        math.fsum(item[2].imag for item in blocks),
    )
    leave = []
    for _, block_denominator, block_numerator in blocks:
        leave_denominator = denominator - block_denominator
        if not math.isfinite(leave_denominator) or leave_denominator <= 0.0:
            raise AssertionError("invalid leave-one denominator")
        leave.append((numerator - block_numerator) / leave_denominator)
    mean = complex(
        math.fsum(item.real for item in leave) / len(leave),
        math.fsum(item.imag for item in leave) / len(leave),
    )
    factor = (len(leave) - 1.0) / len(leave)
    variance_real = factor * math.fsum(
        (item.real - mean.real) ** 2 for item in leave)
    variance_imag = factor * math.fsum(
        (item.imag - mean.imag) ** 2 for item in leave)
    complex_se = math.sqrt(max(0.0, variance_real + variance_imag))
    estimate = numerator / denominator
    if not math.isfinite(complex_se) or not math.isfinite(estimate.real) or \
            not math.isfinite(estimate.imag):
        raise AssertionError("nonfinite jackknife result")
    return {"estimate": estimate, "complex_se": complex_se}


def aggregate_blocks(blocks, group_size):
    if group_size <= 0 or len(blocks) % group_size != 0:
        raise AssertionError("invalid contiguous block aggregation")
    result = []
    for offset in range(0, len(blocks), group_size):
        group = blocks[offset:offset + group_size]
        result.append((
            sum(item[0] for item in group),
            math.fsum(item[1] for item in group),
            complex(math.fsum(item[2].real for item in group),
                    math.fsum(item[2].imag for item in group)),
        ))
    return result


def symmetric_ratio(left, right):
    left = finite(left, "left SE", 0.0)
    right = finite(right, "right SE", 0.0)
    larger = max(left, right)
    smaller = min(left, right)
    if larger == 0.0:
        return 1.0
    if smaller == 0.0:
        return None
    return larger / smaller


def parse_raw(value, name):
    try:
        lines = value.decode("ascii").splitlines()
    except UnicodeDecodeError as error:
        raise AssertionError(f"non-ASCII raw trace: {name}") from error
    headers = [fields(line) for line in lines if line.startswith("MARKOV ")]
    entries = [fields(line) for line in lines if line.startswith("ENTRY ")]
    decisions = [fields(line) for line in lines
                 if line.startswith("DECISION ")]
    if len(headers) != 1 or len(decisions) != 1 or len(entries) != 18:
        raise AssertionError(f"invalid raw trace census: {name}")
    header = headers[0]
    if int(header.get("sample_count", 0)) != 32768 or \
            int(header.get("official_block_count", 0)) != 16 or \
            int(header.get("official_block_length", 0)) != 2048 or \
            int(header.get("diagnostic_block_count", 0)) != 32 or \
            int(header.get("diagnostic_block_length", 0)) != 1024:
        raise AssertionError(f"unexpected frozen partition: {name}")
    by_key = {}
    for item in entries:
        key = (item.get("family"), int(item.get("row", -1)),
               int(item.get("column", -1)))
        if key in by_key:
            raise AssertionError(f"duplicate raw entry: {name}")
        by_key[key] = item
    if set(by_key) != ENTRY_KEYS:
        raise AssertionError(f"incomplete raw entry grid: {name}")
    decision = decisions[0]
    hard_fields = ("support_pass", "tau_pass", "budget_pass",
                   "denominator_pass", "block_pathology_pass")
    if any(decision.get(key) not in ("0", "1") for key in hard_fields) or \
            decision.get("block_stability_pass") not in ("0", "1"):
        raise AssertionError(f"invalid raw decision flags: {name}")
    hard_pass = all(decision[key] == "1" for key in hard_fields)
    if decision.get("p4s_decision") not in ("GO", "STOP") or \
            (decision["p4s_decision"] == "GO") is not hard_pass:
        raise AssertionError(f"P4-S hard decision mismatch: {name}")
    return header, by_key, decision


def parse_reparse(path, header):
    if not path.is_file() or path.is_symlink():
        raise AssertionError(f"missing/nonregular reparse output: {path.name}")
    lines = path.read_text(encoding="ascii").splitlines()
    summaries = [fields(line) for line in lines
                 if line.startswith("REPARSE ")]
    decisions = [fields(line) for line in lines
                 if line.startswith("REPARSE_DECISION ")]
    if len(summaries) != 1 or len(decisions) != 1 or \
            summaries[0].get("all_match") != "1" or \
            decisions[0].get("decision") != "GO":
        raise AssertionError(f"invalid reparse closure: {path.name}")
    summary = summaries[0]
    for key in ("site_count", "qp_total", "sample_count"):
        if int(summary.get(key, -1)) != int(header.get(key, -2)):
            raise AssertionError(f"raw/reparse identity mismatch: {path.name}")
    blocks = {}
    for line in lines:
        if not line.startswith("REPARSE_BLOCK "):
            continue
        item = fields(line)
        family = item.get("family")
        if family not in FAMILIES:
            continue
        block = int(item.get("block", -1))
        key = (family, int(item.get("row", -1)),
               int(item.get("column", -1)))
        target = blocks.setdefault(key, {})
        if block in target:
            raise AssertionError(f"duplicate reparse block: {path.name}")
        target[block] = (
            int(item.get("sample_count", 0)),
            finite(item.get("denominator"), "block denominator", 0.0),
            complex(finite(item.get("numerator_re"), "block numerator"),
                    finite(item.get("numerator_im"), "block numerator")),
        )
    if set(blocks) != ENTRY_KEYS or any(
            set(value) != set(range(BASE_BLOCK_COUNT))
            for value in blocks.values()):
        raise AssertionError(f"incomplete reparse block grid: {path.name}")
    return {key: [value[index] for index in range(BASE_BLOCK_COUNT)]
            for key, value in blocks.items()}


def quantile(values, fraction):
    if not values:
        raise AssertionError("empty quantile input")
    ordered = sorted(values)
    return ordered[round((len(ordered) - 1) * fraction)]


def summarize_traces(traces):
    entries = [item for trace in traces for item in trace["entries"]]
    ratios = [item["registered_ratio"] for item in entries]
    prefix_final_to_initial = [
        item["prefix_se"]["32768"] / item["prefix_se"]["8192"]
        for item in entries if item["prefix_se"]["8192"] > 0.0
    ]
    partition_comparisons = {}
    for field in ("coarser_8x4096_to_16x2048_ratio",
                  "coarser_4x8192_to_8x4096_ratio"):
        values = [item[field] for item in entries if item[field] is not None]
        partition_comparisons[field] = {
            "finite_nonzero_entry_count": len(values),
            "one_zero_pathology_entry_count": len(entries) - len(values),
            "ratio_above_1p25_entry_count": sum(
                value > RATIO_DIAGNOSTIC_MAXIMUM for value in values),
            "ratio_quantiles": {
                str(value): quantile(values, value)
                for value in (0.5, 0.9, 0.95, 0.975, 0.99, 1.0)
            },
        }
    case_summary = []
    for site_count in (4, 6, 8):
        for qp_total in (1, 4):
            group = [item for item in traces
                     if item["site_count"] == site_count and
                     item["qp_total"] == qp_total]
            case_summary.append({
                "site_count": site_count,
                "qp_total": qp_total,
                "trace_count": len(group),
                "ratio_diagnostic_stop_count": sum(
                    not item["registered_ratio_diagnostic_pass"]
                    for item in group),
                "maximum_registered_ratio": max(
                    item["maximum_registered_ratio"] for item in group),
            })
    return {
        "trace_count": len(traces),
        "entry_count": len(entries),
        "ratio_diagnostic_stop_trace_count": sum(
            not item["registered_ratio_diagnostic_pass"] for item in traces),
        "ratio_diagnostic_stop_entry_count": sum(
            item["registered_ratio"] > RATIO_DIAGNOSTIC_MAXIMUM
            for item in entries),
        "p4s_hard_stop_trace_count": sum(
            item["p4s_decision"] != "GO" for item in traces),
        "conservative_budget_stop_entry_count": sum(
            not item["conservative_budget_pass"] for item in entries),
        "one_zero_pathology_stop_entry_count": sum(
            not item["pathology_pass"] for item in entries),
        "registered_ratio_quantiles": {
            str(value): quantile(ratios, value)
            for value in (0.5, 0.9, 0.95, 0.975, 0.99, 1.0)
        },
        "prefix_final_to_initial_se_ratio_quantiles": {
            str(value): quantile(prefix_final_to_initial, value)
            for value in (0.5, 0.9, 0.95, 1.0)
        },
        "prefix_final_to_initial_se_ratio_at_most_0p8_entry_count": sum(
            value <= 0.8 for value in prefix_final_to_initial),
        "prefix_final_to_initial_se_ratio_above_1_entry_count": sum(
            value > 1.0 for value in prefix_final_to_initial),
        "partition_comparisons": partition_comparisons,
        "family_ratio_diagnostic_stop_entry_count": {
            family: sum(item["family"] == family and
                        item["registered_ratio"] >
                        RATIO_DIAGNOSTIC_MAXIMUM for item in entries)
            for family in FAMILIES
        },
        "case_summary": case_summary,
    }


def analyze_dataset(label, archive_path, raw_prefix, reparse_dir):
    archive_path = pathlib.Path(archive_path)
    reparse_dir = pathlib.Path(reparse_dir)
    raw_prefix = validator.normalize_path(raw_prefix).rstrip("/") + "/"
    if not reparse_dir.is_dir() or reparse_dir.is_symlink():
        raise AssertionError("reparse directory must be a regular directory")
    view = validator.ArchiveView(archive_path)
    members = sorted(path for path in view.digests
                     if path.startswith(raw_prefix) and
                     path.endswith(".log"))
    if len(members) != 24:
        raise AssertionError(f"dataset {label} requires exactly 24 traces")
    traces = []
    census = {}
    for member in members:
        name = pathlib.PurePosixPath(member).name
        raw_bytes = view.read_bytes(member)
        header, raw_entries, decision = parse_raw(raw_bytes, name)
        reparse_path = reparse_dir / f"{name}.reparse.log"
        blocks = parse_reparse(reparse_path, header)
        census[f"raw/{name}"] = validator.sha256_bytes(raw_bytes)
        census[f"reparse/{reparse_path.name}"] = \
            validator.sha256_path(reparse_path)
        analyzed_entries = []
        for key in sorted(ENTRY_KEYS):
            raw = raw_entries[key]
            base = blocks[key]
            official = jackknife(base)
            official_recorded = finite(raw.get("official_se"),
                                       "official SE", 0.0)
            diagnostic = finite(raw.get("diagnostic_se"),
                                "diagnostic SE", 0.0)
            if not close(official["complex_se"], official_recorded):
                raise AssertionError(
                    f"official SE reconstruction mismatch: {name} {key}")
            ratio = symmetric_ratio(official_recorded, diagnostic)
            recorded_ratio = finite(raw.get("stability_ratio"),
                                    "stability ratio", 1.0)
            one_zero = ratio is None
            if (one_zero and recorded_ratio != float.fromhex(
                    "0x1.fffffffffffffp+1023")) or \
                    (not one_zero and not close(ratio, recorded_ratio)):
                raise AssertionError(
                    f"registered ratio reconstruction mismatch: {name} {key}")
            if raw.get("budget_pass") not in ("0", "1") or \
                    raw.get("pathology_pass") not in ("0", "1") or \
                    raw.get("stability_pass") not in ("0", "1") or \
                    (raw["pathology_pass"] == "1") is one_zero or \
                    (raw["stability_pass"] == "1") is not \
                    (not one_zero and
                     recorded_ratio <= RATIO_DIAGNOSTIC_MAXIMUM):
                raise AssertionError(
                    f"registered ratio/pathology diagnostic mismatch: {name} {key}")
            conservative_pass = raw["budget_pass"] == "1"
            pathology_pass = raw["pathology_pass"] == "1"
            effective_ratio = recorded_ratio if one_zero else ratio
            prefix_se = {}
            for count in (4, 8, 16):
                prefix_se[str(count * BASE_BLOCK_LENGTH)] = \
                    jackknife(base[:count])["complex_se"]
            partition_se = {"16x2048": official["complex_se"],
                            "32x1024": diagnostic}
            for group_size in (2, 4):
                grouped = aggregate_blocks(base, group_size)
                partition_se[
                    f"{len(grouped)}x{BASE_BLOCK_LENGTH * group_size}"
                ] = jackknife(grouped)["complex_se"]
            analyzed_entries.append({
                "family": key[0], "row": key[1], "column": key[2],
                "official_se": official_recorded,
                "diagnostic_se": diagnostic,
                "registered_ratio": effective_ratio,
                "registered_ratio_diagnostic_pass":
                    not one_zero and
                    effective_ratio <= RATIO_DIAGNOSTIC_MAXIMUM,
                "conservative_budget_ratio": finite(
                    raw.get("conservative_se_budget_ratio"),
                    "conservative budget ratio", 0.0),
                "conservative_budget_pass": conservative_pass,
                "pathology_pass": pathology_pass,
                "prefix_se": prefix_se,
                "partition_se": partition_se,
                "coarser_8x4096_to_16x2048_ratio": symmetric_ratio(
                    partition_se["8x4096"], partition_se["16x2048"]),
                "coarser_4x8192_to_8x4096_ratio": symmetric_ratio(
                    partition_se["4x8192"], partition_se["8x4096"]),
            })
        maximum_ratio = max(item["registered_ratio"]
                            for item in analyzed_entries)
        trace_ratio_pass = all(item["registered_ratio_diagnostic_pass"]
                               for item in analyzed_entries)
        if (decision["block_stability_pass"] == "1") is not trace_ratio_pass:
            raise AssertionError(f"trace ratio diagnostic mismatch: {name}")
        traces.append({
            "trace": name,
            "seed": header.get("seed_hex", "").lower(),
            "site_count": int(header["site_count"]),
            "qp_total": int(header["qp_total"]),
            "p4s_decision": decision["p4s_decision"],
            "registered_ratio_diagnostic_pass": trace_ratio_pass,
            "maximum_registered_ratio": maximum_ratio,
            "registered_ratio_stop_entry_count": sum(
                not item["registered_ratio_diagnostic_pass"]
                for item in analyzed_entries),
            "entries": analyzed_entries,
        })
    census_lines = "".join(f"{census[name]}  {name}\n"
                           for name in sorted(census))
    return {
        "label": label,
        "archive": archive_path.name,
        "archive_sha256": validator.sha256_path(archive_path),
        "raw_prefix": raw_prefix,
        "input_file_count": len(census),
        "input_census_sha256": validator.sha256_bytes(
            census_lines.encode("ascii")),
        "summary": summarize_traces(traces),
        "traces": traces,
    }


def write_atomic(path, value):
    path = pathlib.Path(path)
    if path.exists() or path.is_symlink():
        raise AssertionError("diagnostic output already exists")
    if not path.parent.is_dir() or path.parent.is_symlink():
        raise AssertionError("diagnostic output parent must be a directory")
    with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", prefix=f".{path.name}.",
            suffix=".tmp", dir=path.parent, delete=False) as stream:
        temporary = pathlib.Path(stream.name)
        try:
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write("\n")
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    temporary.rename(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset", nargs=4, action="append", required=True,
        metavar=("LABEL", "ARCHIVE", "RAW_PREFIX", "REPARSE_DIR"))
    parser.add_argument("--output", required=True, type=pathlib.Path)
    args = parser.parse_args()
    labels = [item[0] for item in args.dataset]
    if len(labels) != len(set(labels)):
        raise AssertionError("duplicate diagnostic dataset label")
    datasets = [analyze_dataset(*item) for item in args.dataset]
    result = {
        "schema": 1,
        "artifact": "power_lanczos_p4f_block_stability_diagnostic",
        "decision_eligibility": "diagnostic_only_not_eligible_for_GO",
        "registered_ratio_diagnostic_maximum":
            RATIO_DIAGNOSTIC_MAXIMUM,
        "datasets": datasets,
        "combined_summary": summarize_traces([
            trace for dataset in datasets for trace in dataset["traces"]
        ]),
    }
    write_atomic(args.output, result)
    print(json.dumps({
        "output": args.output.name,
        "dataset_count": len(datasets),
        "trace_count": result["combined_summary"]["trace_count"],
        "ratio_diagnostic_stop_trace_count": result[
            "combined_summary"]["ratio_diagnostic_stop_trace_count"],
        "p4s_hard_stop_trace_count": result[
            "combined_summary"]["p4s_hard_stop_trace_count"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
