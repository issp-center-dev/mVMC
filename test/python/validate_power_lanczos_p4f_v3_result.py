#!/usr/bin/env python3

import argparse
import hashlib
import json
import math
import pathlib
import re
import statistics
import tarfile


DEFAULT_EXPECTED = {
    "solver_policy_sha256":
        "756310504c90825630a0a86dfbf60903024cc1520c73c4d84eb8554e51b03323",
    "confirmatory_policy_sha256":
        "c1348fee957b4eaefa24025fdfb7c79277c800d377f6b7e08ffbfc21aa890f9f",
    "package_manifest_sha256":
        "e495d6721d9d028671247b9b30005f339492ad54ee7552c3c886f3c027d8d5bb",
    "job_script_sha256":
        "a1ef9e0a2ee01c72673101ee2117db171ecab702986cf7ba341b6f4b8c52d66c",
    "generator_sha256":
        "789a40d7e1512c9727c84de85417e4f1f3d0bd9508754c5499d2d6b41b4230fe",
    "source_archive_sha256":
        "d85cb506828022d88e0f022ec957ba7ccfb2596355d844bf89b8393f30c904d8",
    "source_commit": "39d48c272998521d78166d38acf60c6b2e8624e5",
}
HEX64 = re.compile(r"[0-9a-f]{64}")


def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def sha256_path(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def strict_json(value, label):
    def reject_constant(token):
        raise AssertionError(f"non-finite JSON token in {label}: {token}")

    def reject_duplicate_keys(pairs):
        result = {}
        for key, item in pairs:
            if key in result:
                raise AssertionError(
                    f"duplicate JSON key in {label}: {key}")
            result[key] = item
        return result

    try:
        return json.loads(value.decode("utf-8"),
                          parse_constant=reject_constant,
                          object_pairs_hook=reject_duplicate_keys)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise AssertionError(f"invalid JSON in {label}: {error}") from error


def normalize_path(value):
    if not isinstance(value, str) or not value or "\\" in value:
        raise AssertionError(f"unsafe archive path: {value!r}")
    while value.startswith("./"):
        value = value[2:]
    value = value.rstrip("/")
    if not value or value.startswith("/"):
        raise AssertionError(f"unsafe archive path: {value!r}")
    parts = value.split("/")
    if any(part in ("", ".", "..") for part in parts):
        raise AssertionError(f"unsafe archive path: {value!r}")
    normalized = pathlib.PurePosixPath(*parts).as_posix()
    if normalized != value:
        raise AssertionError(f"non-canonical archive path: {value!r}")
    return normalized


def parse_checksum_lines(value, label, prefix=""):
    entries = {}
    try:
        lines = value.decode("ascii").splitlines()
    except UnicodeDecodeError as error:
        raise AssertionError(f"non-ASCII checksum file: {label}") from error
    if not lines:
        raise AssertionError(f"empty checksum file: {label}")
    for line in lines:
        if "  " not in line:
            raise AssertionError(f"malformed checksum line in {label}")
        digest, relative = line.split("  ", 1)
        if HEX64.fullmatch(digest) is None:
            raise AssertionError(f"malformed digest in {label}")
        relative = normalize_path(relative)
        full = normalize_path(f"{prefix}/{relative}") if prefix else relative
        if full in entries:
            raise AssertionError(f"duplicate checksum path in {label}: {full}")
        entries[full] = digest
    return entries


class ArchiveView:
    def __init__(self, path):
        self.path = pathlib.Path(path)
        if not self.path.is_file() or self.path.is_symlink():
            raise AssertionError("result archive must be a regular file")
        self.members = {}
        self.digests = {}
        with tarfile.open(self.path, "r:gz") as archive:
            for member in archive.getmembers():
                if member.name in (".", "./") and member.isdir():
                    continue
                relative = normalize_path(member.name)
                if relative in self.members:
                    raise AssertionError(
                        f"duplicate result archive member: {relative}")
                if member.isdir():
                    self.members[relative] = ("directory", member.size)
                    continue
                if not member.isfile():
                    raise AssertionError(
                        f"non-regular result archive member: {relative}")
                source = archive.extractfile(member)
                if source is None:
                    raise AssertionError(f"unreadable archive member: {relative}")
                digest = hashlib.sha256()
                for chunk in iter(lambda: source.read(1024 * 1024), b""):
                    digest.update(chunk)
                self.members[relative] = ("file", member.size)
                self.digests[relative] = digest.hexdigest()

    def read_bytes(self, relative, maximum=64 * 1024 * 1024):
        relative = normalize_path(relative)
        member = self.members.get(relative)
        if member is None or member[0] != "file":
            raise AssertionError(f"missing result file: {relative}")
        if member[1] > maximum:
            raise AssertionError(f"result file too large to parse: {relative}")
        with tarfile.open(self.path, "r:gz") as archive:
            matches = [item for item in archive.getmembers()
                       if item.name not in (".", "./") and
                       normalize_path(item.name) == relative]
            if len(matches) != 1:
                raise AssertionError(f"missing/duplicate result file: {relative}")
            source = archive.extractfile(matches[0])
            if source is None:
                raise AssertionError(f"unreadable result file: {relative}")
            value = source.read(maximum + 1)
            if len(value) != member[1] or len(value) > maximum or \
                    sha256_bytes(value) != self.digests[relative]:
                raise AssertionError(
                    f"result file changed while validating: {relative}")
            return value


def verify_outer_checksum(archive_path, checksum_path):
    archive_path = pathlib.Path(archive_path)
    checksum_path = pathlib.Path(checksum_path)
    if not checksum_path.is_file() or checksum_path.is_symlink():
        raise AssertionError("result checksum must be a regular file")
    entries = parse_checksum_lines(checksum_path.read_bytes(),
                                   "outer checksum")
    if entries != {archive_path.name: sha256_path(archive_path)}:
        raise AssertionError("result archive outer checksum mismatch")
    return entries[archive_path.name]


def finite_number(value, label, minimum=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or \
            not math.isfinite(value):
        raise AssertionError(f"invalid finite number for {label}")
    value = float(value)
    if minimum is not None and value < minimum:
        raise AssertionError(f"out-of-range number for {label}")
    return value


def same_float(left, right):
    left = finite_number(left, "recorded numeric")
    right = finite_number(right, "recomputed numeric")
    return math.isclose(left, right, rel_tol=2.0e-14, abs_tol=1.0e-15)


def dimension_mandatory(dimension, solver_policy):
    if dimension.get("dimension") not in (2, 3):
        raise AssertionError("invalid GEVP dimension")
    prefixes = dimension.get("prefixes")
    if not isinstance(prefixes, list) or len(prefixes) != 3 or \
            [item.get("sample_count") for item in prefixes] != \
            [8192, 16384, 32768]:
        raise AssertionError("invalid prefix census")
    for item in prefixes:
        finite_number(item.get("energy_se"), "prefix energy SE", 0.0)
    exact = dimension.get("exact", {})
    anti = finite_number(exact.get("antihermitian_residual"),
                         "exact anti-Hermitian residual", 0.0)
    anti_limit = finite_number(
        solver_policy["hermitianization"]
        ["exact_maximum_relative_residual"], "anti-Hermitian limit", 0.0)
    return all((
        dimension.get("solver_pass") is True,
        dimension.get("cutoff_scan", {}).get("full_pass") is True,
        dimension.get("cutoff_scan", {}).get("exact_pass") is True,
        all(item.get("exact_energy_pass") is True for item in prefixes),
        all(item.get("rank_stability_pass") is True for item in prefixes),
        dimension.get("coefficient_gate", {}).get("passed") is True,
        anti <= anti_limit,
    ))


def trace_mandatory(case, solver_policy):
    diagnostics = case.get("matrix_diagnostics")
    if not isinstance(diagnostics, list) or \
            {item.get("family") for item in diagnostics} != {"K", "S", "B"}:
        raise AssertionError("invalid matrix diagnostic census")
    dimensions = case.get("dimensions")
    if not isinstance(dimensions, list) or \
            {item.get("dimension") for item in dimensions} != {2, 3}:
        raise AssertionError("invalid dimension census")
    identity = case.get("confirmatory_identity", {})
    checks = identity.get("checks")
    expected_checks = {
        "driver_schema", "session_schema", "fixture", "seed", "case",
        "sample_count", "rank_count", "cache", "rho", "rng_stream",
        "persistent_session", "proposal", "partition", "p4s",
    }
    if not isinstance(checks, dict) or set(checks) != expected_checks or \
            any(not isinstance(value, bool) for value in checks.values()):
        raise AssertionError("confirmatory identity checks are incomplete")
    identity_pass = identity.get("passed")
    if not isinstance(identity_pass, bool) or \
            identity_pass is not all(checks.values()):
        raise AssertionError("confirmatory identity aggregate mismatch")
    return all((
        identity_pass,
        case.get("reconstruction", {}).get("passed") is True,
        all(item.get("passed") is True for item in diagnostics),
        all(dimension_mandatory(item, solver_policy)
            for item in dimensions),
    ))


def theil_sen(errors):
    samples = (8192, 16384, 32768)
    slopes = []
    for left in range(3):
        for right in range(left + 1, 3):
            slopes.append((math.log(errors[right]) - math.log(errors[left])) /
                          (math.log(samples[right]) - math.log(samples[left])))
    return statistics.median(slopes)


def recompute_pooled(group, dimension_value, solver_policy):
    if len(group) != 4:
        raise AssertionError("pooled group must contain four seeds")
    dimensions = []
    per_seed = []
    seed_patterns_pass = True
    tolerance_relative = finite_number(
        solver_policy["exact_gate"]["energy_tolerance_relative"],
        "exact energy tolerance", 0.0)
    for case in group:
        matches = [item for item in case["dimensions"]
                   if item["dimension"] == dimension_value]
        if len(matches) != 1:
            raise AssertionError("missing/duplicate pooled dimension")
        dimension = matches[0]
        dimensions.append(dimension)
        errors = [finite_number(item["energy_se"], "pooled energy SE", 0.0)
                  for item in dimension["prefixes"]]
        tolerance = tolerance_relative * max(
            1.0, abs(finite_number(dimension["exact"]["energy"],
                                   "exact energy")))
        if all(value == 0.0 for value in errors):
            pattern = "all_zero"
            pattern_pass = all(
                finite_number(item["exact_energy_error"],
                              "exact energy error", 0.0) <= tolerance
                for item in dimension["prefixes"])
        elif any(value == 0.0 for value in errors):
            pattern = "mixed_zero"
            pattern_pass = False
        else:
            pattern = "all_nonzero"
            pattern_pass = True
        seed_patterns_pass &= pattern_pass
        per_seed.append({
            "seed": case["confirmatory_identity"]["seed"].lower(),
            "errors": errors,
            "pattern": pattern,
            "pattern_pass": pattern_pass,
        })
    pooled = [math.sqrt(sum(item["errors"][index] ** 2
                            for item in per_seed) / 4.0)
              for index in range(3)]
    if all(value == 0.0 for value in pooled):
        zero_rule = "all_zero"
        ratio = None
        slope = None
        convergence = all(
            finite_number(item["exact_energy_error"],
                          "exact energy error", 0.0) <=
            tolerance_relative * max(
                1.0, abs(finite_number(dimension["exact"]["energy"],
                                       "exact energy")))
            for dimension in dimensions for item in dimension["prefixes"])
    elif any(value == 0.0 for value in pooled):
        zero_rule = "mixed_zero"
        ratio = None
        slope = None
        convergence = False
    else:
        zero_rule = "all_nonzero"
        ratio = pooled[-1] / pooled[0]
        slope = theil_sen(pooled)
        convergence = ratio <= 0.8 and slope < 0.0
    return {
        "per_seed": per_seed,
        "pooled": pooled,
        "zero_rule": zero_rule,
        "ratio": ratio,
        "slope": slope,
        "seed_patterns_pass": seed_patterns_pass,
        "convergence": convergence,
        "passed": seed_patterns_pass and convergence,
    }


def verify_pooled_record(recorded, recomputed, seeds):
    if recorded.get("seed_count") != 4 or \
            recorded.get("sample_counts") != [8192, 16384, 32768] or \
            recorded.get("maximum_ratio") != 0.8 or \
            recorded.get("zero_rule") != recomputed["zero_rule"]:
        raise AssertionError("pooled gate contract mismatch")
    values = recorded.get("pooled_energy_se")
    if not isinstance(values, list) or len(values) != 3 or \
            not all(same_float(left, right) for left, right in
                    zip(values, recomputed["pooled"])):
        raise AssertionError("pooled SE recomputation mismatch")
    for field, expected_value in (
            ("final_to_initial_ratio", recomputed["ratio"]),
            ("theil_sen_log_se_log_n_slope", recomputed["slope"])):
        actual = recorded.get(field)
        if expected_value is None:
            if actual is not None:
                raise AssertionError(f"unexpected pooled {field}")
        elif actual is None or not same_float(actual, expected_value):
            raise AssertionError(f"pooled {field} mismatch")
    diagnostics = recorded.get("per_seed_diagnostics")
    if not isinstance(diagnostics, list) or len(diagnostics) != 4 or \
            [item.get("seed", "").lower() for item in diagnostics] != seeds:
        raise AssertionError("pooled seed diagnostics mismatch")
    for actual, expected in zip(diagnostics, recomputed["per_seed"]):
        errors = actual.get("energy_se")
        if not isinstance(errors, list) or len(errors) != 3 or \
                not all(same_float(left, right) for left, right in
                        zip(errors, expected["errors"])) or \
                actual.get("zero_pattern") != expected["pattern"] or \
                actual.get("zero_pattern_pass") is not \
                expected["pattern_pass"]:
            raise AssertionError("per-seed pooled diagnostic mismatch")
    for field, expected_value in (
            ("seed_zero_patterns_pass", recomputed["seed_patterns_pass"]),
            ("convergence_pass", recomputed["convergence"]),
            ("passed", recomputed["passed"])):
        if recorded.get(field) is not expected_value:
            raise AssertionError(f"pooled boolean mismatch: {field}")


def verify_package_provenance(view, expected):
    package_manifest_path = "workflow/package_manifest.sha256"
    if view.digests.get(package_manifest_path) != \
            expected["package_manifest_sha256"]:
        raise AssertionError("package manifest provenance mismatch")
    if view.digests.get("workflow/p4f_v3_genkai_job.sh") != \
            expected["job_script_sha256"]:
        raise AssertionError("job script provenance mismatch")
    package_entries = parse_checksum_lines(
        view.read_bytes(package_manifest_path), "package manifest")
    required = {
        "frozen_inputs/power_lanczos_zero_support_p4f_solver_policy.json":
            expected["solver_policy_sha256"],
        "frozen_inputs/power_lanczos_zero_support_p4f_v3_confirmatory_policy.json":
            expected["confirmatory_policy_sha256"],
        "mVMC-p4f-v3-39d48c2-source.tar.gz":
            expected["source_archive_sha256"],
        "workflow/p4f_v3_genkai_job.sh": expected["job_script_sha256"],
    }
    for path, digest in required.items():
        if package_entries.get(path) != digest:
            raise AssertionError(f"package provenance entry mismatch: {path}")
    source_checksum = parse_checksum_lines(
        view.read_bytes("workflow/source_archive.sha256"),
        "source archive checksum")
    if source_checksum != {
        "mVMC-p4f-v3-39d48c2-source.tar.gz":
            expected["source_archive_sha256"]
    }:
        raise AssertionError("source archive provenance mismatch")
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
    build_lines = build_config.splitlines()
    for required in required_lines:
        key = required.split("=", 1)[0]
        matches = [line for line in build_lines
                   if line.startswith(f"{key}=")]
        if matches != [required]:
            raise AssertionError(
                f"build configuration provenance mismatch: {key}")
    ctest = view.read_bytes("workflow/genkai_focused_ctest.log").decode(
        "utf-8", errors="strict")
    if "100% tests passed" not in ctest or \
            re.search(r"0 tests failed out of 5\b", ctest) is None:
        raise AssertionError("focused CTest result mismatch")
    inventory = view.read_bytes(
        "workflow/remote_scratch_inventory.txt").decode("ascii")
    if re.search(r"^raw_trace_count=24$", inventory,
                 flags=re.MULTILINE) is None:
        raise AssertionError("remote scratch trace census mismatch")


def validate_result(archive_path, checksum_path, expected=None):
    expected = dict(DEFAULT_EXPECTED if expected is None else expected)
    outer_sha = verify_outer_checksum(archive_path, checksum_path)
    view = ArchiveView(archive_path)
    if sha256_path(pathlib.Path(archive_path)) != outer_sha:
        raise AssertionError("result archive changed while validating")
    ledger_path = "artifact-ledger.sha256"
    ledger = parse_checksum_lines(view.read_bytes(ledger_path),
                                  "artifact ledger")
    actual = {path: digest for path, digest in view.digests.items()
              if path != ledger_path}
    if ledger != actual:
        raise AssertionError("artifact ledger file census/checksum mismatch")
    nested_path = "confirmatory/checksums.sha256"
    nested = parse_checksum_lines(view.read_bytes(nested_path),
                                  "confirmatory checksums",
                                  prefix="confirmatory")
    nested_actual = {path: digest for path, digest in view.digests.items()
                     if path.startswith("confirmatory/") and
                     path != nested_path}
    if nested != nested_actual:
        raise AssertionError("confirmatory checksum census mismatch")

    raw_paths = sorted(path for path in view.digests
                       if path.startswith("raw_traces/") and
                       path.endswith(".log"))
    reparse_paths = sorted(path for path in view.digests
                           if path.startswith("confirmatory/raw_reparse/") and
                           path.endswith(".reparse.log"))
    if len(raw_paths) != 24 or len(reparse_paths) != 24:
        raise AssertionError("raw/reparse trace census mismatch")

    solver_path = "confirmatory/solver_policy.json"
    confirmatory_path = "confirmatory/confirmatory_policy.json"
    if view.digests.get(solver_path) != expected["solver_policy_sha256"] or \
            view.digests.get(confirmatory_path) != \
            expected["confirmatory_policy_sha256"]:
        raise AssertionError("frozen result policy hash mismatch")
    solver_policy = strict_json(view.read_bytes(solver_path), solver_path)
    confirmatory_policy = strict_json(
        view.read_bytes(confirmatory_path), confirmatory_path)
    if confirmatory_policy.get("solver_policy_sha256") != \
            expected["solver_policy_sha256"] or \
            confirmatory_policy.get("policy_id") != \
            "power_lanczos_zero_support_p4f_confirmatory_v3":
        raise AssertionError("frozen policy binding mismatch")
    verify_package_provenance(view, expected)

    evidence = strict_json(view.read_bytes(
        "confirmatory/p4f_v3_confirmatory_evidence.json"), "evidence")
    decision = strict_json(view.read_bytes(
        "confirmatory/p4f_v3_confirmatory_decision.json"), "decision")
    metadata = strict_json(view.read_bytes(
        "confirmatory/metadata.json"), "metadata")
    if evidence.get("metadata") != metadata or \
            metadata.get("solver_policy_sha256") != \
            expected["solver_policy_sha256"] or \
            metadata.get("confirmatory_policy_sha256") != \
            expected["confirmatory_policy_sha256"] or \
            metadata.get("generator_sha256") != expected["generator_sha256"]:
        raise AssertionError("result metadata provenance mismatch")
    if evidence.get("policy_id") != confirmatory_policy["policy_id"] or \
            decision.get("policy_id") != confirmatory_policy["policy_id"]:
        raise AssertionError("result policy ID mismatch")

    dataset = confirmatory_policy["confirmatory_dataset"]
    seeds = [value.lower() for value in dataset["seed_order"]]
    case_order = [(item["site_count"], item["qp_total"])
                  for item in dataset["case_order"]]
    if len(seeds) != 4 or len(set(seeds)) != 4 or len(case_order) != 6 or \
            len(set(case_order)) != 6 or dataset.get("trace_count") != 24:
        raise AssertionError("frozen confirmatory grid mismatch")
    expected_keys = {(seed, site, qp) for seed in seeds
                     for site, qp in case_order}
    cases = evidence.get("cases")
    if not isinstance(cases, list) or len(cases) != 24:
        raise AssertionError("evidence trace census mismatch")
    by_key = {}
    stopped_traces = []
    for case in cases:
        identity = case.get("confirmatory_identity", {})
        key = (str(identity.get("seed", "")).lower(),
               case.get("site_count"), case.get("qp_total"))
        if identity.get("site_count") != case.get("site_count") or \
                identity.get("qp_total") != case.get("qp_total") or \
                case.get("sample_count") != 32768:
            raise AssertionError("confirmatory case/identity mismatch")
        if key in by_key:
            raise AssertionError("duplicate evidence trace identity")
        by_key[key] = case
        trace = case.get("trace")
        if not isinstance(trace, str) or "/" in trace or "\\" in trace:
            raise AssertionError("unsafe evidence trace name")
        raw_path = f"raw_traces/{trace}"
        reparse_relative = case.get("reparse_output")
        if reparse_relative != f"raw_reparse/{trace}.reparse.log":
            raise AssertionError("reparse output binding mismatch")
        reparse_path = f"confirmatory/{reparse_relative}"
        if view.digests.get(raw_path) != case.get("trace_sha256") or \
                view.digests.get(reparse_path) != \
                case.get("reparse_output_sha256"):
            raise AssertionError("raw/reparse SHA binding mismatch")
        mandatory = trace_mandatory(case, solver_policy)
        if case.get("mandatory_gate_pass") is not mandatory or \
                case.get("decision") != ("GO" if mandatory else "STOP"):
            raise AssertionError("trace mandatory decision mismatch")
        if not mandatory:
            stopped_traces.append(trace)
    if set(by_key) != expected_keys or evidence.get("identity_grid_pass") \
            is not True:
        raise AssertionError("confirmatory identity grid mismatch")

    recorded_pooled = evidence.get("pooled_se_convergence_gates")
    if not isinstance(recorded_pooled, list) or len(recorded_pooled) != 12:
        raise AssertionError("pooled gate census mismatch")
    pooled_by_key = {}
    for item in recorded_pooled:
        key = (item.get("site_count"), item.get("qp_total"),
               item.get("dimension"))
        if key in pooled_by_key:
            raise AssertionError("duplicate pooled gate")
        pooled_by_key[key] = item
    stopped_pooled = []
    for site, qp in case_order:
        group = [by_key[(seed, site, qp)] for seed in seeds]
        for dimension in (2, 3):
            key = (site, qp, dimension)
            if key not in pooled_by_key:
                raise AssertionError("missing pooled gate")
            recomputed = recompute_pooled(group, dimension, solver_policy)
            verify_pooled_record(pooled_by_key[key], recomputed, seeds)
            if not recomputed["passed"]:
                stopped_pooled.append({
                    "site_count": site, "qp_total": qp,
                    "dimension": dimension,
                })
    final = "GO" if not stopped_traces and not stopped_pooled else "STOP"
    if evidence.get("trace_count") != 24 or \
            evidence.get("stopped_traces") != stopped_traces or \
            evidence.get("stopped_pooled_gates") != stopped_pooled or \
            evidence.get("decision") != final:
        raise AssertionError("evidence aggregate decision mismatch")
    expected_decision = {
        "schema": 1,
        "policy_id": confirmatory_policy["policy_id"],
        "solver_policy_sha256": expected["solver_policy_sha256"],
        "confirmatory_policy_sha256":
            expected["confirmatory_policy_sha256"],
        "trace_count": 24,
        "mandatory_trace_go_count": 24 - len(stopped_traces),
        "mandatory_trace_stop_count": len(stopped_traces),
        "pooled_gate_go_count": 12 - len(stopped_pooled),
        "pooled_gate_stop_count": len(stopped_pooled),
        "stopped_traces": stopped_traces,
        "stopped_pooled_gates": stopped_pooled,
        "decision": final,
    }
    if decision != expected_decision:
        raise AssertionError("decision artifact mismatch")
    generator_lines = view.read_bytes(
        "workflow/genkai_confirmatory_generator.log").decode(
            "utf-8", errors="strict").splitlines()
    if not generator_lines or strict_json(
            generator_lines[-1].encode("utf-8"), "generator final line") != \
            decision:
        raise AssertionError("generator log decision mismatch")
    return {
        "valid": True,
        "archive_sha256": outer_sha,
        "decision": final,
        "trace_count": 24,
        "mandatory_trace_stop_count": len(stopped_traces),
        "pooled_gate_stop_count": len(stopped_pooled),
        "artifact_file_count": len(view.digests),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True, type=pathlib.Path)
    parser.add_argument("--checksum", required=True, type=pathlib.Path)
    args = parser.parse_args()
    result = validate_result(args.archive, args.checksum)
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
