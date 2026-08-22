#!/usr/bin/env python3
"""Generate, execute, assemble, and validate P6-C2 scaling records.

All process records are append-only.  The runner never replaces a censored
record and never uses a system temporary directory.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import itertools
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from jsonschema import Draft202012Validator


PROTOCOL_ID = "p6c2-operator-count-scaling-v2"
TIMER_ID = "p6c2_observable_source_batch_v2"
ERRATUM_ID = "power-lanczos-p6c2-observable-scaling-erratum-v2"
STRATA = ["SC-ONEBODY", "SC-QUARTIC", "SC-MIXED", "SC-LOW-REUSE"]
COUNTS = [1, 32, 128, 512]
FAMILIES = ["cisajs", "cisajscktaltex", "cisajscktalt"]
PROFILES = {
    "SC-ONEBODY": ("electronic", "real", ["cisajs"]),
    "SC-QUARTIC": ("pure_spin", "real", ["cisajscktaltex"]),
    "SC-MIXED": ("electronic", "complex", FAMILIES),
    "SC-LOW-REUSE": ("electronic", "complex", FAMILIES),
}
OFFICIAL_RESULT = Path(
    "docs/policies/power-lanczos-p6c2-observable-scaling-evidence-v2.json"
)
SHA_KEYS = [
    "source_diff_sha256",
    "binary_sha256",
    "build_sha256",
    "environment_sha256",
]
PROVENANCE_KEYS = ["source_commit", *SHA_KEYS]
DATA_DIR = Path(__file__).resolve().parent / "data"
RAW_SCHEMA_PATH = DATA_DIR / "power_lanczos_observable_scaling_raw_schema_v2.json"
EVIDENCE_SCHEMA_PATH = DATA_DIR / "power_lanczos_observable_scaling_evidence_schema_v2.json"


class CampaignError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CampaignError(message)


def canonical_bytes(value) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def sha_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha_file(path: Path) -> str:
    return sha_bytes(path.read_bytes())


def load_json(path: Path):
    def reject_duplicates(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise CampaignError(f"duplicate JSON key in {path}: {key}")
            result[key] = value
        return result

    return json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=reject_duplicates,
        parse_constant=lambda value: (_ for _ in ()).throw(
            CampaignError(f"nonfinite JSON constant in {path}: {value}")
        ),
    )


def write_new_json(path: Path, value) -> None:
    require(not path.exists(), f"append-only output already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    path.write_text(payload, encoding="utf-8")


def one_body_rows():
    rows = []
    for spin in range(2):
        for output_site in range(16):
            for input_site in range(16):
                rows.append((output_site, spin, input_site, spin))
    return rows


def quartic_rows():
    rows = []
    for first_spin in range(2):
        for first_output in range(16):
            for first_input in range(16):
                for second_spin in range(2):
                    for second_output in range(16):
                        for second_input in range(16):
                            rows.append(
                                (
                                    first_output,
                                    first_spin,
                                    first_input,
                                    first_spin,
                                    second_output,
                                    second_spin,
                                    second_input,
                                    second_spin,
                                )
                            )
                            if len(rows) == 1536:
                                return rows
    raise AssertionError("unreachable")


def pure_spin_quartic_rows():
    rows = []
    for up_site in range(16):
        for down_site in range(16):
            rows.append(
                (up_site, 0, up_site, 0, down_site, 1, down_site, 1)
            )
    for output_site in range(16):
        for input_site in range(16):
            if output_site != input_site:
                rows.append(
                    (
                        output_site,
                        0,
                        input_site,
                        0,
                        input_site,
                        1,
                        output_site,
                        1,
                    )
                )
    for first_site in range(16):
        for second_site in range(first_site + 1, 16):
            rows.append(
                (
                    first_site,
                    0,
                    first_site,
                    0,
                    second_site,
                    0,
                    second_site,
                    0,
                )
            )
    require(len(rows) >= 512, "pure-spin quartic row supply")
    return rows


def low_reuse_rows():
    one = []
    quartic = []
    for spin in range(2):
        for occupied in range(5):
            for empty in range(5, 10):
                one.append((empty, spin, occupied, spin))
    for up_occupied in range(5):
        for up_empty in range(5, 10):
            for down_occupied in range(5):
                for down_empty in range(5, 10):
                    quartic.append(
                        (
                            up_empty,
                            0,
                            up_occupied,
                            0,
                            down_empty,
                            1,
                            down_occupied,
                            1,
                        )
                    )
    queues = {
        "cisajs": list(one),
        "cisajscktaltex": list(quartic[0::2]),
        "cisajscktalt": list(quartic[1::2]),
    }
    ordered = []
    while len(ordered) < 512:
        for family in FAMILIES:
            if queues[family]:
                ordered.append((family, queues[family].pop(0)))
                if len(ordered) == 512:
                    break
    return ordered


def point_rows(stratum: str, count: int):
    require(stratum in STRATA and count in COUNTS, "registered point")
    if stratum == "SC-ONEBODY":
        return [("cisajs", row) for row in one_body_rows()[:count]]
    if stratum == "SC-QUARTIC":
        return [
            ("cisajscktaltex", row) for row in pure_spin_quartic_rows()[:count]
        ]
    if stratum == "SC-LOW-REUSE":
        return low_reuse_rows()[:count]
    queues = {
        "cisajs": one_body_rows(),
        "cisajscktaltex": quartic_rows()[0::2],
        "cisajscktalt": quartic_rows()[1::2],
    }
    ordered = []
    while len(ordered) < count:
        for family in FAMILIES:
            if queues[family]:
                ordered.append((family, queues[family].pop(0)))
                if len(ordered) == count:
                    break
    return ordered


def observable_bytes(family: str, rows) -> bytes:
    labels = {
        "cisajs": "NCisAjs",
        "cisajscktaltex": "NCisAjsCktAlt",
        "cisajscktalt": "NCisAjsCktAltDC",
    }
    lines = [
        "====================",
        f"{labels[family]} {len(rows)}",
        "====================",
        "observable rows",
        "====================",
    ]
    lines.extend(" ".join(str(value) for value in row) for row in rows)
    return ("\n".join(lines) + "\n").encode("ascii")


def generate_inputs(output_dir: Path) -> None:
    require(not output_dir.exists(), f"input output directory already exists: {output_dir}")
    output_dir.mkdir(parents=True)
    top = {
        "schema_version": 2,
        "protocol_id": PROTOCOL_ID,
        "timer_id": TIMER_ID,
        "erratum_id": ERRATUM_ID,
        "points": [],
    }
    for stratum in STRATA:
        model, arithmetic, _ = PROFILES[stratum]
        for count in COUNTS:
            point_dir = output_dir / stratum / f"R{count:03d}"
            point_dir.mkdir(parents=True)
            rows = point_rows(stratum, count)
            by_family = {family: [] for family in FAMILIES}
            for family, row in rows:
                by_family[family].append(row)
            files = {}
            for family in FAMILIES:
                if not by_family[family]:
                    files[family] = None
                    continue
                name = f"{family}.def"
                path = point_dir / name
                payload = observable_bytes(family, by_family[family])
                path.write_bytes(payload)
                files[family] = {
                    "path": name,
                    "sha256": sha_bytes(payload),
                    "bytes": len(payload),
                    "count": len(by_family[family]),
                }
            identity = {
                "stratum_id": stratum,
                "operator_count": count,
                "model": model,
                "arithmetic": arithmetic,
                "Nsite": 16,
                "NQPFull": 8,
                "family_counts": {
                    family: len(by_family[family]) for family in FAMILIES
                },
                "files": files,
            }
            identity["input_sha256"] = sha_bytes(canonical_bytes(identity))
            manifest = {
                "schema_version": 2,
                "protocol_id": PROTOCOL_ID,
                "timer_id": TIMER_ID,
                "erratum_id": ERRATUM_ID,
                **identity,
            }
            manifest_path = point_dir / "input-manifest.json"
            manifest_path.write_text(
                json.dumps(manifest, indent=2, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            top["points"].append(
                {
                    "stratum_id": stratum,
                    "operator_count": count,
                    "manifest": str(manifest_path.relative_to(output_dir)),
                    "manifest_sha256": sha_file(manifest_path),
                    "input_sha256": identity["input_sha256"],
                }
            )
    top["input_set_sha256"] = sha_bytes(canonical_bytes(top["points"]))
    (output_dir / "campaign-input-manifest.json").write_text(
        json.dumps(top, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def git_output(repo: Path, *args: str) -> bytes:
    return subprocess.check_output(["git", "-C", str(repo), *args])


def source_diff_digest(repo: Path) -> str:
    digest = hashlib.sha256()
    digest.update(b"tracked-diff\0")
    digest.update(git_output(repo, "diff", "--binary", "HEAD"))
    untracked = git_output(
        repo, "ls-files", "--others", "--exclude-standard", "-z"
    ).split(b"\0")
    for raw_name in sorted(name for name in untracked if name):
        relative = raw_name.decode("utf-8")
        path = repo / relative
        digest.update(b"untracked\0" + raw_name + b"\0")
        digest.update(path.read_bytes())
    return digest.hexdigest()


def create_snapshot(
    repo: Path | None,
    driver: Path,
    build_dir: Path,
    output: Path,
    source_commit_override: str | None = None,
    source_diff_override: str | None = None,
) -> None:
    if repo is not None:
        repo = repo.resolve()
    driver = driver.resolve()
    build_dir = build_dir.resolve()
    require(driver.is_file(), f"driver does not exist: {driver}")
    cache = build_dir / "CMakeCache.txt"
    link = build_dir / "test/CMakeFiles/power_lanczos_observable_scaling_driver.dir/link.txt"
    require(cache.is_file() and link.is_file(), "build identity inputs")
    require(
        (source_commit_override is None) == (source_diff_override is None),
        "snapshot source overrides must be supplied together",
    )
    if source_commit_override is None:
        require(repo is not None, "snapshot repository is required")
        source_commit = git_output(repo, "rev-parse", "HEAD").decode().strip()
        source_diff_sha256 = source_diff_digest(repo)
    else:
        source_commit = source_commit_override
        source_diff_sha256 = source_diff_override
    require(len(source_commit) == 40, "source commit")
    require(
        source_diff_sha256 is not None and len(source_diff_sha256) == 64,
        "source diff sha256",
    )
    build_identity = {
        "cmake_cache_sha256": sha_file(cache),
        "link_command_sha256": sha_file(link),
        "driver_sha256": sha_file(driver),
    }
    environment = {
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "python": platform.python_version(),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "strict_fp": True,
    }
    snapshot = {
        "schema_version": 2,
        "source_commit": source_commit,
        "source_diff_sha256": source_diff_sha256,
        "binary_sha256": sha_file(driver),
        "build_sha256": sha_bytes(canonical_bytes(build_identity)),
        "environment_sha256": sha_bytes(canonical_bytes(environment)),
        "build_identity": build_identity,
        "environment": environment,
    }
    write_new_json(output, snapshot)


def point_manifest(input_root: Path, stratum: str, count: int):
    path = input_root / stratum / f"R{count:03d}" / "input-manifest.json"
    value = load_json(path)
    require(
        value["protocol_id"] == PROTOCOL_ID
        and value["timer_id"] == TIMER_ID
        and value["erratum_id"] == ERRATUM_ID
        and value["stratum_id"] == stratum
        and value["operator_count"] == count,
        "point manifest identity",
    )
    for family in FAMILIES:
        item = value["files"][family]
        if item is not None:
            file_path = path.parent / item["path"]
            require(
                file_path.is_file()
                and sha_file(file_path) == item["sha256"]
                and file_path.stat().st_size == item["bytes"],
                f"point input bytes: {family}",
            )
    identity = {
        key: value[key]
        for key in [
            "stratum_id",
            "operator_count",
            "model",
            "arithmetic",
            "Nsite",
            "NQPFull",
            "family_counts",
            "files",
        ]
    }
    require(
        sha_bytes(canonical_bytes(identity)) == value["input_sha256"],
        "point input identity digest",
    )
    return path, value


def timing_seed(stratum: str, count: int, kind: str, replicate: int) -> str:
    return sha_bytes(
        canonical_bytes([PROTOCOL_ID, TIMER_ID, stratum, count, kind, replicate])
    )


def censored_record(
    identity,
    manifest_path: Path,
    manifest,
    stratum: str,
    count: int,
    saved_sources: int,
    replicate: int,
    kind: str,
    mpi_world_size: int,
    reason: str,
):
    require(
        reason in {"timeout", "nonzero_exit", "scheduler_eviction", "missing_metric", "identity_mismatch"},
        "registered censor reason",
    )
    seed = timing_seed(stratum, count, kind, replicate)
    record_id = (
        f"{kind}-{stratum}-R{count:03d}-rep{replicate}-mpi{mpi_world_size}"
    )
    return {
        "schema_version": 2,
        "record_id": record_id,
        "protocol_id": PROTOCOL_ID,
        "timer_id": TIMER_ID,
        "erratum_id": ERRATUM_ID,
        "request": {
            "stratum_id": stratum,
            "operator_count": count,
            "saved_source_count": saved_sources,
            "replicate_ordinal": replicate,
            "run_kind": kind,
            "mpi_world_size": mpi_world_size,
            "timing_seed_sha256": seed,
        },
        "input": {
            "manifest_sha256": sha_file(manifest_path),
            "input_sha256": manifest["input_sha256"],
            "family_counts": manifest["family_counts"],
        },
        "provenance": {
            "source_commit": identity["source_commit"],
            **{key: identity[key] for key in SHA_KEYS},
        },
        "process": {
            "driver_basename": "power_lanczos_observable_scaling_driver",
            "exit_code": None,
            "stdout_sha256": sha_bytes(b""),
            "stderr_sha256": sha_bytes(b""),
            "censored": True,
            "censor_reason": reason,
        },
        "metrics": None,
        "artifacts": [],
    }


def seal_missing_records(
    input_root: Path, identity_path: Path, raw_dir: Path, reason: str
) -> int:
    identity = load_json(identity_path)
    for key in SHA_KEYS:
        require(len(identity[key]) == 64, f"snapshot {key}")
    requests = [
        (stratum, count, 0, "warmup", 1)
        for stratum in STRATA
        for count in COUNTS
    ]
    requests += [
        (stratum, count, replicate, "measured", 1)
        for stratum in STRATA
        for count in COUNTS
        for replicate in range(7)
    ]
    requests += [
        (stratum, 512, 0, "parity", mpi_world_size)
        for stratum in STRATA
        for mpi_world_size in [1, 2, 4]
    ]
    created = 0
    raw_dir.mkdir(parents=True, exist_ok=True)
    for stratum, count, replicate, kind, mpi_world_size in requests:
        record_id = (
            f"{kind}-{stratum}-R{count:03d}-rep{replicate}-mpi{mpi_world_size}"
        )
        record_path = raw_dir / f"{record_id}.json"
        if record_path.exists():
            continue
        manifest_path, manifest = point_manifest(input_root, stratum, count)
        record = censored_record(
            identity,
            manifest_path,
            manifest,
            stratum,
            count,
            64,
            replicate,
            kind,
            mpi_world_size,
            reason,
        )
        write_new_json(record_path, record)
        validate_raw_record(record_path, raw_dir)
        created += 1
    return created


def run_one(args) -> None:
    args.driver = args.driver.resolve()
    args.input_root = args.input_root.resolve()
    args.identity = args.identity.resolve()
    args.raw_dir = args.raw_dir.resolve()
    if args.mpiexec is not None:
        args.mpiexec = args.mpiexec.resolve()
    identity = load_json(args.identity)
    for key in SHA_KEYS:
        require(len(identity[key]) == 64, f"snapshot {key}")
    manifest_path, manifest = point_manifest(
        args.input_root, args.stratum, args.operator_count
    )
    record_id = (
        f"{args.kind}-{args.stratum}-R{args.operator_count:03d}-"
        f"rep{args.replicate}-mpi{args.mpi_world_size}"
    )
    record_path = args.raw_dir / f"{record_id}.json"
    require(not record_path.exists(), f"append-only record exists: {record_path}")
    args.raw_dir.mkdir(parents=True, exist_ok=True)
    work_root = args.raw_dir / ".work"
    work_root.mkdir(exist_ok=True)
    work_dir = work_root / record_id
    require(not work_dir.exists(), f"run work directory exists: {work_dir}")
    work_dir.mkdir()
    coefficient_tmp = work_dir / "coefficient.bin"
    final_tmp = work_dir / "final.bin"
    metrics_tmp = work_dir / "metrics.json"
    seed = timing_seed(args.stratum, args.operator_count, args.kind, args.replicate)
    input_paths = []
    for family in FAMILIES:
        item = manifest["files"][family]
        input_paths.append("-" if item is None else str(manifest_path.parent / item["path"]))
    driver_command = [
        str(args.driver),
        args.stratum,
        str(args.operator_count),
        str(args.saved_sources),
        str(args.replicate),
        args.kind,
        seed,
        *input_paths,
        str(coefficient_tmp),
        str(final_tmp),
        str(metrics_tmp),
    ]
    command = driver_command
    if args.mpi_world_size > 1:
        require(args.mpiexec is not None, "--mpiexec required for MPI run")
        command = [str(args.mpiexec), "-np", str(args.mpi_world_size), *driver_command]
    environment = os.environ.copy()
    environment.update(
        {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
        }
    )
    exit_code = None
    stdout = b""
    stderr = b""
    censor_reason = None
    try:
        completed = subprocess.run(
            command,
            cwd=work_dir,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=args.timeout,
            check=False,
        )
        exit_code = completed.returncode
        stdout = completed.stdout
        stderr = completed.stderr
        if exit_code != 0:
            censor_reason = "nonzero_exit"
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout or b""
        stderr = exc.stderr or b""
        censor_reason = "timeout"
    request = {
        "stratum_id": args.stratum,
        "operator_count": args.operator_count,
        "saved_source_count": args.saved_sources,
        "replicate_ordinal": args.replicate,
        "run_kind": args.kind,
        "mpi_world_size": args.mpi_world_size,
        "timing_seed_sha256": seed,
    }
    record = {
        "schema_version": 2,
        "record_id": record_id,
        "protocol_id": PROTOCOL_ID,
        "timer_id": TIMER_ID,
        "erratum_id": ERRATUM_ID,
        "request": request,
        "input": {
            "manifest_sha256": sha_file(manifest_path),
            "input_sha256": manifest["input_sha256"],
            "family_counts": manifest["family_counts"],
        },
        "provenance": {
            "source_commit": identity["source_commit"],
            **{key: identity[key] for key in SHA_KEYS},
        },
        "process": {
            "driver_basename": args.driver.name,
            "exit_code": exit_code,
            "stdout_sha256": sha_bytes(stdout),
            "stderr_sha256": sha_bytes(stderr),
            "censored": censor_reason is not None,
            "censor_reason": censor_reason,
        },
        "metrics": None,
        "artifacts": [],
    }
    if censor_reason is None:
        if not (
            metrics_tmp.is_file()
            and coefficient_tmp.is_file()
            and final_tmp.is_file()
        ):
            censor_reason = "missing_metric"
        else:
            try:
                metrics = load_json(metrics_tmp)
                validate_driver_metrics(metrics, request, manifest)
            except (CampaignError, OSError, json.JSONDecodeError):
                censor_reason = "identity_mismatch"
    if censor_reason is None:
        coefficient_name = f"{record_id}.coefficient.bin"
        final_name = f"{record_id}.final.bin"
        coefficient_path = args.raw_dir / coefficient_name
        final_path = args.raw_dir / final_name
        require(
            not coefficient_path.exists() and not final_path.exists(),
            "payload output exists",
        )
        os.replace(coefficient_tmp, coefficient_path)
        os.replace(final_tmp, final_path)
        artifacts = [
            {
                "role": "coefficient_blocks",
                "path": coefficient_name,
                "bytes": coefficient_path.stat().st_size,
                "sha256": sha_file(coefficient_path),
            },
            {
                "role": "final_blocks",
                "path": final_name,
                "bytes": final_path.stat().st_size,
                "sha256": sha_file(final_path),
            },
        ]
        semantic = hashlib.sha256()
        semantic.update(
            canonical_bytes(
                {
                    "stratum_id": args.stratum,
                    "operator_count": args.operator_count,
                    "operator_census_sha256": metrics[
                        "operator_census_sha256"
                    ],
                    "timing_seed_sha256": seed,
                }
            )
        )
        semantic.update(coefficient_path.read_bytes())
        semantic.update(final_path.read_bytes())
        metrics["normalized_semantic_payload_sha256"] = semantic.hexdigest()
        record["metrics"] = metrics
        record["artifacts"] = artifacts
    else:
        record["process"]["censored"] = True
        record["process"]["censor_reason"] = censor_reason
    write_new_json(record_path, record)
    shutil.rmtree(work_dir)
    if not any(work_root.iterdir()):
        work_root.rmdir()
    validate_raw_record(record_path, args.raw_dir)
    print(
        f"PASS raw scaling record {record_id} censored={record['process']['censored']}"
    )


def validate_driver_metrics(metrics, request, manifest) -> None:
    require(metrics["schema_version"] == 2, "driver metrics schema")
    require(metrics["timer_id"] == TIMER_ID, "driver timer")
    for key in [
        "stratum_id",
        "operator_count",
        "saved_source_count",
        "replicate_ordinal",
        "run_kind",
        "mpi_world_size",
        "timing_seed_sha256",
    ]:
        require(metrics[key] == request[key], f"driver request binding: {key}")
    require(
        (metrics["model"], metrics["arithmetic"])
        == PROFILES[request["stratum_id"]][:2],
        "driver profile",
    )
    require(metrics["Nsite"] == 16 and metrics["NQPFull"] == 8, "driver envelope")
    require(metrics["family_counts"] == manifest["family_counts"], "driver family census")
    require(
        metrics["artifact_payload_bytes"]
        == metrics["coefficient_payload_bytes"] + metrics["final_payload_bytes"]
        == 80 * metrics["block_count"] * metrics["operator_count"],
        "driver payload formula",
    )
    require(
        metrics["wall_seconds"] == metrics["setup_seconds"] + metrics["active_seconds"]
        and metrics["active_seconds"] > 0
        and metrics["setup_seconds"] >= 0,
        "driver timer decomposition",
    )
    require(metrics["status"] == "PASS", "driver status")
    require(
        0 <= metrics["minimum_unique_target_count"]
        <= metrics["maximum_unique_target_count"]
        <= metrics["operator_count"],
        "driver target counts",
    )
    if request["stratum_id"] == "SC-LOW-REUSE":
        require(
            10 * metrics["minimum_unique_target_count"]
            >= 9 * metrics["operator_count"],
            "effective low-reuse gate",
        )


def validate_raw_record(path: Path, raw_root: Path | None = None):
    record = load_json(path)
    raw_schema = load_json(RAW_SCHEMA_PATH)
    Draft202012Validator.check_schema(raw_schema)
    errors = sorted(
        Draft202012Validator(raw_schema).iter_errors(record),
        key=lambda error: list(error.absolute_path),
    )
    require(not errors, f"raw JSON schema: {errors[0].message if errors else ''}")
    require(record["schema_version"] == 2, "raw schema")
    require(
        record["protocol_id"] == PROTOCOL_ID
        and record["timer_id"] == TIMER_ID
        and record["erratum_id"] == ERRATUM_ID,
        "raw protocol identity",
    )
    request = record["request"]
    require(
        request["stratum_id"] in STRATA
        and request["operator_count"] in COUNTS
        and request["run_kind"] in {"warmup", "measured", "parity", "selftest"}
        and request["mpi_world_size"] in {1, 2, 4},
        "raw request envelope",
    )
    require(
        sum(record["input"]["family_counts"].values())
        == request["operator_count"],
        "raw input family/operator count",
    )
    for key in SHA_KEYS:
        require(len(record["provenance"][key]) == 64, f"raw provenance {key}")
    censored = record["process"]["censored"]
    require(censored == (record["process"]["censor_reason"] is not None), "raw censor mapping")
    if censored:
        require(record["metrics"] is None and record["artifacts"] == [], "censored payload absence")
    else:
        metrics = record["metrics"]
        require(metrics is not None and len(record["artifacts"]) == 2, "successful raw artifacts")
        require(metrics["normalized_semantic_payload_sha256"], "semantic payload SHA")
        for key in [
            "stratum_id",
            "operator_count",
            "saved_source_count",
            "replicate_ordinal",
            "run_kind",
            "mpi_world_size",
            "timing_seed_sha256",
        ]:
            require(metrics[key] == request[key], f"raw metrics/request binding: {key}")
        require(
            (metrics["model"], metrics["arithmetic"])
            == PROFILES[request["stratum_id"]][:2],
            "raw effective profile",
        )
        require(
            sum(metrics["family_counts"].values()) == request["operator_count"],
            "raw family/operator count",
        )
        require(
            metrics["family_counts"] == record["input"]["family_counts"],
            "raw metrics/input family binding",
        )
        require(
            metrics["artifact_payload_bytes"]
            == 80 * metrics["block_count"] * request["operator_count"],
            "raw block payload formula",
        )
        require(
            math.isclose(
                metrics["wall_seconds"],
                metrics["setup_seconds"] + metrics["active_seconds"],
                rel_tol=1e-12,
                abs_tol=1e-14,
            ),
            "raw wall decomposition",
        )
        require(
            0 <= metrics["minimum_unique_target_count"]
            <= metrics["maximum_unique_target_count"]
            <= request["operator_count"],
            "raw unique-target bound",
        )
        if request["stratum_id"] == "SC-LOW-REUSE":
            require(
                10 * metrics["minimum_unique_target_count"]
                >= 9 * request["operator_count"],
                "raw low-reuse minimum",
            )
        if raw_root is not None:
            for artifact in record["artifacts"]:
                artifact_path = raw_root / artifact["path"]
                require(
                    artifact_path.is_file()
                    and artifact_path.stat().st_size == artifact["bytes"]
                    and sha_file(artifact_path) == artifact["sha256"],
                    f"raw artifact bytes: {artifact['path']}",
                )
            require(
                sum(item["bytes"] for item in record["artifacts"])
                == metrics["artifact_payload_bytes"],
                "raw payload byte sum",
            )
            semantic = hashlib.sha256()
            semantic.update(
                canonical_bytes(
                    {
                        "stratum_id": request["stratum_id"],
                        "operator_count": request["operator_count"],
                        "operator_census_sha256": metrics[
                            "operator_census_sha256"
                        ],
                        "timing_seed_sha256": request["timing_seed_sha256"],
                    }
                )
            )
            by_role = {item["role"]: item for item in record["artifacts"]}
            semantic.update((raw_root / by_role["coefficient_blocks"]["path"]).read_bytes())
            semantic.update((raw_root / by_role["final_blocks"]["path"]).read_bytes())
            require(
                semantic.hexdigest()
                == metrics["normalized_semantic_payload_sha256"],
                "normalized semantic payload digest",
            )
    return record


def point_summary(records, stratum: str):
    selected = [record for record in records if record["stratum_id"] == stratum]
    points = []
    for count in COUNTS:
        group = [record for record in selected if record["operator_count"] == count]
        require(len(group) == 7, f"summary replicate count: {stratum} R={count}")
        if any(record["censored"] for record in group):
            upper = None
        else:
            upper = max(
                record["active_seconds"] / record["saved_source_count"]
                for record in group
            )
        points.append(
            {
                "operator_count": count,
                "active_seconds_per_saved_source_upper": upper,
            }
        )
    uncensored = all(point["active_seconds_per_saved_source_upper"] is not None for point in points)
    slope = None
    intercept = None
    if uncensored:
        slopes = [0.0]
        for left_index, left in enumerate(points):
            for right in points[left_index + 1 :]:
                slopes.append(
                    (
                        right["active_seconds_per_saved_source_upper"]
                        - left["active_seconds_per_saved_source_upper"]
                    )
                    / (right["operator_count"] - left["operator_count"])
                )
        slope = max(slopes)
        intercept = max(
            0.0,
            max(
                point["active_seconds_per_saved_source_upper"]
                - slope * point["operator_count"]
                for point in points
            ),
        )
    return {
        "stratum_id": stratum,
        "point_uppers": points,
        "slope_upper_seconds_per_source_operator": slope,
        "intercept_upper_seconds_per_source": intercept,
        "setup_upper_seconds": max(record["setup_seconds"] for record in selected if not record["censored"]) if any(not record["censored"] for record in selected) else None,
        "rss_minus_allocation_headroom_upper_bytes": max(
            max(0, record["peak_RSS_bytes"] - record["exact_allocated_bytes"])
            for record in selected
            if not record["censored"]
        ) if any(not record["censored"] for record in selected) else None,
        "passed": uncensored
        and all(
            record["wall_seconds"] <= 3600
            and record["peak_RSS_bytes"] <= 8589934592
            and record["exact_allocated_bytes"] <= 8589934592
            for record in selected
        ),
    }


def assemble(raw_dir: Path, output: Path, allow_official_output: bool) -> None:
    if output.as_posix().endswith(OFFICIAL_RESULT.as_posix()):
        require(allow_official_output, "official result path requires --allow-official-output")
    raw = []
    for path in sorted(raw_dir.glob("*.json")):
        record = validate_raw_record(path, raw_dir)
        record["_record_file_bytes"] = path.stat().st_size
        raw.append(record)
    require(raw, "raw records are absent")
    warmup = [item for item in raw if item["request"]["run_kind"] == "warmup"]
    measured = [item for item in raw if item["request"]["run_kind"] == "measured"]
    parity = [item for item in raw if item["request"]["run_kind"] == "parity"]
    order_key = lambda item: (
        STRATA.index(item["request"]["stratum_id"]),
        COUNTS.index(item["request"]["operator_count"]),
        item["request"]["replicate_ordinal"],
        item["request"]["mpi_world_size"],
    )
    warmup.sort(key=order_key)
    measured.sort(key=order_key)
    parity.sort(key=order_key)
    expected_warmup = [(s, r, 0, 1) for s in STRATA for r in COUNTS]
    expected_measured = [(s, r, rep, 1) for s in STRATA for r in COUNTS for rep in range(7)]
    expected_parity = [(s, 512, 0, size) for s in STRATA for size in [1, 2, 4]]
    identity = lambda item: (
        item["request"]["stratum_id"],
        item["request"]["operator_count"],
        item["request"]["replicate_ordinal"],
        item["request"]["mpi_world_size"],
    )
    require([identity(item) for item in warmup] == expected_warmup, "exact warmup grid/order")
    require([identity(item) for item in measured] == expected_measured, "exact measured grid/order")
    require([identity(item) for item in parity] == expected_parity, "exact parity grid/order")
    provenance_values = {
        key: {item["provenance"][key] for item in measured + parity}
        for key in PROVENANCE_KEYS
    }
    require(all(len(values) == 1 for values in provenance_values.values()), "campaign provenance identity")
    warmup_records = []
    for item in warmup:
        request = item["request"]
        metrics = item["metrics"]
        metadata_bytes = item["_record_file_bytes"]
        payload_bytes = metrics["artifact_payload_bytes"] if metrics else 0
        warmup_records.append(
            {
                "record_id": item["record_id"],
                "stratum_id": request["stratum_id"],
                "operator_count": request["operator_count"],
                "block_count": metrics["block_count"] if metrics else 0,
                "artifact_payload_bytes": payload_bytes,
                "artifact_metadata_bytes": metadata_bytes,
                "artifact_bytes": payload_bytes + metadata_bytes,
                "artifact_files": len(item["artifacts"]) + 1,
                "censored": item["process"]["censored"],
                "censor_reason": item["process"]["censor_reason"],
            }
        )
    records = []
    for item in measured:
        request = item["request"]
        metrics = item["metrics"]
        censored = item["process"]["censored"]
        if not censored:
            require(request["saved_source_count"] == 64, "measured saved-source count")
            if request["stratum_id"] == "SC-LOW-REUSE":
                require(
                    10 * metrics["minimum_unique_target_count"]
                    >= 9 * request["operator_count"],
                    "assembled low reuse",
                )
        records.append(
            {
                "stratum_id": request["stratum_id"],
                "operator_count": request["operator_count"],
                "replicate_ordinal": request["replicate_ordinal"],
                "model": PROFILES[request["stratum_id"]][0],
                "arithmetic": PROFILES[request["stratum_id"]][1],
                "Nsite": 16,
                "NQPFull": 8,
                "family_counts": metrics["family_counts"] if metrics else item["input"]["family_counts"],
                "unique_target_count": metrics["minimum_unique_target_count"] if metrics else None,
                "saved_source_count": request["saved_source_count"],
                "input_sha256": item["input"]["input_sha256"],
                "operator_census_sha256": metrics["operator_census_sha256"] if metrics else None,
                "timing_seed_sha256": request["timing_seed_sha256"],
                **item["provenance"],
                "mpi_world_size": 1,
                "timer_id": TIMER_ID,
                "setup_seconds": metrics["setup_seconds"] if metrics else 0.0,
                "active_seconds": metrics["active_seconds"] if metrics else 0.0,
                "wall_seconds": metrics["wall_seconds"] if metrics else 0.0,
                "peak_RSS_bytes": metrics["peak_RSS_bytes"] if metrics else 0,
                "exact_allocated_bytes": metrics["exact_allocated_bytes"] if metrics else 0,
                "block_count": metrics["block_count"] if metrics else 0,
                "artifact_payload_bytes": metrics["artifact_payload_bytes"] if metrics else 0,
                "artifact_metadata_bytes": item["_record_file_bytes"],
                "artifact_bytes": (metrics["artifact_payload_bytes"] + item["_record_file_bytes"]) if metrics else item["_record_file_bytes"],
                "artifact_files": (len(item["artifacts"]) + 1),
                "normalized_semantic_payload_sha256": metrics["normalized_semantic_payload_sha256"] if metrics else None,
                "censored": censored,
                "censor_reason": item["process"]["censor_reason"],
            }
        )
    summaries = [point_summary(records, stratum) for stratum in STRATA]
    parity_records = []
    for item in parity:
        request = item["request"]
        metrics = item["metrics"]
        parity_records.append(
            {
                "record_id": item["record_id"],
                "stratum_id": request["stratum_id"],
                "operator_count": 512,
                "mpi_world_size": request["mpi_world_size"],
                "normalized_semantic_payload_sha256": metrics["normalized_semantic_payload_sha256"] if metrics else None,
                "block_count": metrics["block_count"] if metrics else 0,
                "artifact_payload_bytes": metrics["artifact_payload_bytes"] if metrics else 0,
                "artifact_metadata_bytes": item["_record_file_bytes"],
                "artifact_bytes": (metrics["artifact_payload_bytes"] + item["_record_file_bytes"]) if metrics else item["_record_file_bytes"],
                "artifact_files": len(item["artifacts"]) + 1,
                "censored": item["process"]["censored"],
                "censor_reason": item["process"]["censor_reason"],
            }
        )
    parity_pass = True
    for stratum in STRATA:
        group = [item for item in parity_records if item["stratum_id"] == stratum]
        parity_pass &= all(not item["censored"] for item in group) and len(
            {item["normalized_semantic_payload_sha256"] for item in group}
        ) == 1
    warmup_pass = all(not item["censored"] for item in warmup_records)
    artifact_records = warmup_records + records + parity_records
    artifact_bytes = sum(item["artifact_bytes"] for item in artifact_records)
    artifact_files = sum(item["artifact_files"] for item in artifact_records)
    passed = (
        warmup_pass
        and all(item["passed"] for item in summaries)
        and parity_pass
        and artifact_bytes <= 107374182400
        and artifact_files <= 25000
    )
    evidence = {
        "schema_version": 2,
        "evidence_id": "power-lanczos-p6c2-observable-scaling-evidence-v2",
        "erratum_id": ERRATUM_ID,
        "result_observed": True,
        "production_authorized": False,
        "status": "PASS" if passed else "FAIL",
        "protocol_id": PROTOCOL_ID,
        "timer_id": TIMER_ID,
        "provenance": {key: next(iter(values)) for key, values in provenance_values.items()},
        "warmup_record_ids": [item["record_id"] for item in warmup],
        "warmup_records": warmup_records,
        "records": records,
        "summaries": summaries,
        "parity_records": parity_records,
        "campaign_artifact_bytes": artifact_bytes,
        "campaign_artifact_files": artifact_files,
        "passed": passed,
    }
    validate_evidence(evidence)
    write_new_json(output, evidence)
    print(f"PASS assembled scaling evidence status={evidence['status']} records=112 parity=12")


def synthetic_raw_campaign(raw_dir: Path) -> None:
    provenance = {
        "source_commit": "1" * 40,
        "source_diff_sha256": "a" * 64,
        "binary_sha256": "b" * 64,
        "build_sha256": "c" * 64,
        "environment_sha256": "d" * 64,
    }
    raw_dir.mkdir(parents=True)
    payload_bases = {}
    for count in COUNTS:
        coefficient = raw_dir / f"base-R{count:03d}.coefficient.bin"
        final = raw_dir / f"base-R{count:03d}.final.bin"
        with coefficient.open("wb") as stream:
            stream.truncate(64 * 32 * count)
        with final.open("wb") as stream:
            stream.truncate(16 * 32 * count)
        payload_bases[count] = (coefficient, final)

    def make_record(stratum: str, count: int, replicate: int, kind: str, mpi_size: int):
        record_id = f"{kind}-{stratum}-R{count:03d}-rep{replicate}-mpi{mpi_size}"
        request = {
            "stratum_id": stratum,
            "operator_count": count,
            "saved_source_count": 64,
            "replicate_ordinal": replicate,
            "run_kind": kind,
            "mpi_world_size": mpi_size,
            "timing_seed_sha256": timing_seed(stratum, count, kind, replicate),
        }
        family_counts = {family: 0 for family in FAMILIES}
        for family, _ in point_rows(stratum, count):
            family_counts[family] += 1
        coefficient_name = f"{record_id}.coefficient.bin"
        final_name = f"{record_id}.final.bin"
        coefficient_path = raw_dir / coefficient_name
        final_path = raw_dir / final_name
        os.link(payload_bases[count][0], coefficient_path)
        os.link(payload_bases[count][1], final_path)
        coefficient_sha = sha_file(coefficient_path)
        final_sha = sha_file(final_path)
        census_sha = sha_bytes(canonical_bytes([stratum, count, "census"]))
        semantic = hashlib.sha256()
        semantic.update(
            canonical_bytes(
                {
                    "stratum_id": stratum,
                    "operator_count": count,
                    "operator_census_sha256": census_sha,
                    "timing_seed_sha256": request["timing_seed_sha256"],
                }
            )
        )
        semantic.update(coefficient_path.read_bytes())
        semantic.update(final_path.read_bytes())
        setup = 0.01 + replicate * 1e-5
        active = 64 * (0.001 + count * 1e-5 + replicate * 1e-8)
        metrics = {
            "schema_version": 2,
            "timer_id": TIMER_ID,
            **request,
            "model": PROFILES[stratum][0],
            "arithmetic": PROFILES[stratum][1],
            "Nsite": 16,
            "NQPFull": 8,
            "family_counts": family_counts,
            "raw_observable_census_sha256": sha_bytes(canonical_bytes([stratum, count, "raw"])),
            "operator_census_sha256": census_sha,
            "setup_seconds": setup,
            "active_seconds": active,
            "wall_seconds": setup + active,
            "peak_RSS_bytes": 200000000,
            "exact_allocated_bytes": 100000000,
            "minimum_unique_target_count": count if stratum == "SC-LOW-REUSE" else min(count, 65),
            "maximum_unique_target_count": count if stratum == "SC-LOW-REUSE" else min(count, 65),
            "block_count": 32,
            "coefficient_payload_bytes": 64 * 32 * count,
            "final_payload_bytes": 16 * 32 * count,
            "artifact_payload_bytes": 80 * 32 * count,
            "artifact_files": 2,
            "amplitude_callback_calls": 1,
            "normalized_semantic_payload_fnv1a64": "1" * 16,
            "status": "PASS",
            "normalized_semantic_payload_sha256": semantic.hexdigest(),
        }
        record = {
            "schema_version": 2,
            "record_id": record_id,
            "protocol_id": PROTOCOL_ID,
            "timer_id": TIMER_ID,
            "erratum_id": ERRATUM_ID,
            "request": request,
            "input": {
                "manifest_sha256": sha_bytes(canonical_bytes([stratum, count, "manifest"])),
                "input_sha256": sha_bytes(canonical_bytes([stratum, count, "input"])),
                "family_counts": family_counts,
            },
            "provenance": provenance,
            "process": {
                "driver_basename": "power_lanczos_observable_scaling_driver",
                "exit_code": 0,
                "stdout_sha256": sha_bytes(b""),
                "stderr_sha256": sha_bytes(b""),
                "censored": False,
                "censor_reason": None,
            },
            "metrics": metrics,
            "artifacts": [
                {"role": "coefficient_blocks", "path": coefficient_name, "bytes": coefficient_path.stat().st_size, "sha256": coefficient_sha},
                {"role": "final_blocks", "path": final_name, "bytes": final_path.stat().st_size, "sha256": final_sha},
            ],
        }
        (raw_dir / f"{record_id}.json").write_text(
            json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )

    for stratum in STRATA:
        for count in COUNTS:
            make_record(stratum, count, 0, "warmup", 1)
            for replicate in range(7):
                make_record(stratum, count, replicate, "measured", 1)
        for mpi_size in [1, 2, 4]:
            make_record(stratum, 512, 0, "parity", mpi_size)


def self_test_assembler() -> None:
    local_tmp = Path.cwd() / "tmp"
    local_tmp.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="p6c2-scaling-assembler.", dir=local_tmp))
    try:
        raw_dir = work / "raw"
        evidence_path = work / "evidence.json"
        synthetic_raw_campaign(raw_dir)
        assemble(raw_dir, evidence_path, False)
        validate_evidence(load_json(evidence_path))
    finally:
        resolved = work.resolve()
        if resolved.parent == local_tmp.resolve() and not resolved.is_symlink():
            shutil.rmtree(resolved)
        if local_tmp.exists() and not any(local_tmp.iterdir()):
            local_tmp.rmdir()


def synthetic_evidence():
    provenance = {
        "source_commit": "1" * 40,
        "source_diff_sha256": "a" * 64,
        "binary_sha256": "b" * 64,
        "build_sha256": "c" * 64,
        "environment_sha256": "d" * 64,
    }
    records = []
    for stratum in STRATA:
        model, arithmetic, _ = PROFILES[stratum]
        for count in COUNTS:
            counts = {family: 0 for family in FAMILIES}
            for family, _ in point_rows(stratum, count):
                counts[family] += 1
            for replicate in range(7):
                setup = 0.01 + replicate * 1e-5
                active = 64 * (0.001 + count * 1e-5 + replicate * 1e-8)
                blocks = 32
                payload = 80 * blocks * count
                metadata = 1024
                records.append(
                    {
                        "stratum_id": stratum,
                        "operator_count": count,
                        "replicate_ordinal": replicate,
                        "model": model,
                        "arithmetic": arithmetic,
                        "Nsite": 16,
                        "NQPFull": 8,
                        "family_counts": counts,
                        "unique_target_count": count if stratum == "SC-LOW-REUSE" else min(count, 65),
                        "saved_source_count": 64,
                        "input_sha256": sha_bytes(canonical_bytes([stratum, count, "input"])),
                        "operator_census_sha256": sha_bytes(canonical_bytes([stratum, count, "census"])),
                        "timing_seed_sha256": timing_seed(stratum, count, "measured", replicate),
                        **provenance,
                        "mpi_world_size": 1,
                        "timer_id": TIMER_ID,
                        "setup_seconds": setup,
                        "active_seconds": active,
                        "wall_seconds": setup + active,
                        "peak_RSS_bytes": 200000000,
                        "exact_allocated_bytes": 100000000,
                        "block_count": blocks,
                        "artifact_payload_bytes": payload,
                        "artifact_metadata_bytes": metadata,
                        "artifact_bytes": payload + metadata,
                        "artifact_files": 3,
                        "normalized_semantic_payload_sha256": sha_bytes(canonical_bytes([stratum, count, replicate, "payload"])),
                        "censored": False,
                        "censor_reason": None,
                    }
                )
    summaries = [point_summary(records, stratum) for stratum in STRATA]
    warmup_records = []
    for stratum in STRATA:
        for count in COUNTS:
            payload = 80 * 32 * count
            warmup_records.append(
                {
                    "record_id": f"warmup-{stratum}-R{count:03d}-rep0-mpi1",
                    "stratum_id": stratum,
                    "operator_count": count,
                    "block_count": 32,
                    "artifact_payload_bytes": payload,
                    "artifact_metadata_bytes": 1024,
                    "artifact_bytes": payload + 1024,
                    "artifact_files": 3,
                    "censored": False,
                    "censor_reason": None,
                }
            )
    parity = []
    for stratum in STRATA:
        payload_sha = sha_bytes(canonical_bytes([stratum, 512, "parity"] ))
        for size in [1, 2, 4]:
            parity.append(
                {
                    "record_id": f"parity-{stratum}-R512-rep0-mpi{size}",
                    "stratum_id": stratum,
                    "operator_count": 512,
                    "mpi_world_size": size,
                    "normalized_semantic_payload_sha256": payload_sha,
                    "block_count": 32,
                    "artifact_payload_bytes": 80 * 32 * 512,
                    "artifact_metadata_bytes": 1024,
                    "artifact_bytes": 80 * 32 * 512 + 1024,
                    "artifact_files": 3,
                    "censored": False,
                    "censor_reason": None,
                }
            )
    artifact_records = warmup_records + records + parity
    artifact_bytes = sum(item["artifact_bytes"] for item in artifact_records)
    artifact_files = sum(item["artifact_files"] for item in artifact_records)
    return {
        "schema_version": 2,
        "evidence_id": "power-lanczos-p6c2-observable-scaling-evidence-v2",
        "erratum_id": ERRATUM_ID,
        "result_observed": True,
        "production_authorized": False,
        "status": "PASS",
        "protocol_id": PROTOCOL_ID,
        "timer_id": TIMER_ID,
        "provenance": provenance,
        "warmup_record_ids": [
            f"warmup-{stratum}-R{count:03d}-rep0-mpi1"
            for stratum in STRATA
            for count in COUNTS
        ],
        "warmup_records": warmup_records,
        "records": records,
        "summaries": summaries,
        "parity_records": parity,
        "campaign_artifact_bytes": artifact_bytes,
        "campaign_artifact_files": artifact_files,
        "passed": True,
    }


def self_test() -> None:
    base = synthetic_evidence()
    validate_evidence(base)
    censored = copy.deepcopy(base)
    censored_record = censored["records"][0]
    censored_record.update(
        {
            "unique_target_count": None,
            "setup_seconds": 0.0,
            "active_seconds": 0.0,
            "wall_seconds": 0.0,
            "peak_RSS_bytes": 0,
            "exact_allocated_bytes": 0,
            "block_count": 0,
            "artifact_payload_bytes": 0,
            "artifact_metadata_bytes": 1024,
            "artifact_bytes": 1024,
            "artifact_files": 1,
            "normalized_semantic_payload_sha256": None,
            "censored": True,
            "censor_reason": "timeout",
        }
    )
    censored["summaries"] = [
        point_summary(censored["records"], stratum) for stratum in STRATA
    ]
    censored["campaign_artifact_bytes"] = sum(
        item["artifact_bytes"]
        for item in censored["warmup_records"]
        + censored["records"]
        + censored["parity_records"]
    )
    censored["campaign_artifact_files"] = sum(
        item["artifact_files"]
        for item in censored["warmup_records"]
        + censored["records"]
        + censored["parity_records"]
    )
    censored["passed"] = False
    censored["status"] = "FAIL"
    validate_evidence(censored)
    mutations = []

    def mutation(mutator):
        candidate = copy.deepcopy(base)
        mutator(candidate)
        mutations.append(candidate)

    mutation(lambda value: value["records"].pop())
    low_index = 3 * len(COUNTS) * 7
    mutation(lambda value: value["records"][low_index].__setitem__("model", "pure_spin"))
    mutation(lambda value: value["records"][low_index + 7].__setitem__("unique_target_count", 1))
    mutation(lambda value: value["records"][1].__setitem__("timing_seed_sha256", value["records"][0]["timing_seed_sha256"]))
    mutation(lambda value: value["summaries"][0]["point_uppers"][0].__setitem__("active_seconds_per_saved_source_upper", 9.0))
    mutation(lambda value: value["parity_records"][2].__setitem__("normalized_semantic_payload_sha256", "f" * 64))
    mutation(lambda value: value.__setitem__("campaign_artifact_bytes", value["campaign_artifact_bytes"] + 1))

    def censor_without_recompute(value):
        value["records"][0]["censored"] = True
        value["records"][0]["censor_reason"] = "timeout"

    mutation(censor_without_recompute)
    mutation(lambda value: value["records"][0].__setitem__("artifact_payload_bytes", value["records"][0]["artifact_payload_bytes"] + 16))
    mutation(lambda value: value.__setitem__("campaign_artifact_files", 25001))
    mutation(lambda value: value["warmup_records"][0].__setitem__("artifact_bytes", 1))
    for index, candidate in enumerate(mutations, 1):
        try:
            validate_evidence(candidate)
        except (CampaignError, KeyError, TypeError, ValueError):
            continue
        raise CampaignError(f"negative evidence mutation accepted: NEG-S{index:02d}")
    self_test_assembler()
    print("PASS scaling campaign self-test positive=3 negative=11 assembler_grid=140")


def validate_evidence(evidence) -> None:
    schema = load_json(EVIDENCE_SCHEMA_PATH)
    Draft202012Validator.check_schema(schema)
    errors = sorted(
        Draft202012Validator(schema).iter_errors(evidence),
        key=lambda error: list(error.absolute_path),
    )
    require(
        not errors,
        "evidence JSON schema at {}: {}".format(
            "/".join(str(item) for item in errors[0].absolute_path)
            if errors
            else "",
            errors[0].message if errors else "",
        ),
    )
    require(evidence["schema_version"] == 2, "evidence schema")
    require(evidence["erratum_id"] == ERRATUM_ID, "evidence erratum")
    require(evidence["protocol_id"] == PROTOCOL_ID and evidence["timer_id"] == TIMER_ID, "evidence protocol")
    require(evidence["production_authorized"] is False, "evidence production boundary")
    require(
        len(evidence["warmup_record_ids"]) == 16
        and len(evidence["warmup_records"]) == 16,
        "evidence warmup count",
    )
    require(len(evidence["records"]) == 112 and len(evidence["parity_records"]) == 12, "evidence record census")
    expected = [(s, r, rep) for s in STRATA for r in COUNTS for rep in range(7)]
    expected_warmup = [(s, r) for s in STRATA for r in COUNTS]
    warmup = evidence["warmup_records"]
    require(
        [(item["stratum_id"], item["operator_count"]) for item in warmup]
        == expected_warmup,
        "evidence warmup grid",
    )
    require(
        [item["record_id"] for item in warmup]
        == evidence["warmup_record_ids"],
        "evidence warmup ID binding",
    )
    for record in warmup:
        require(
            record["censored"] == (record["censor_reason"] is not None),
            "evidence warmup censor mapping",
        )
        require(
            record["artifact_bytes"]
            == record["artifact_payload_bytes"]
            + record["artifact_metadata_bytes"],
            "evidence warmup artifact decomposition",
        )
        if record["censored"]:
            require(
                record["block_count"] == 0
                and record["artifact_payload_bytes"] == 0
                and record["artifact_files"] == 1,
                "evidence censored warmup artifacts",
            )
        else:
            require(
                record["block_count"] >= 32
                and record["artifact_payload_bytes"]
                == 80 * record["block_count"] * record["operator_count"]
                and record["artifact_files"] == 3,
                "evidence warmup artifacts",
            )
    require(
        [(item["stratum_id"], item["operator_count"], item["replicate_ordinal"]) for item in evidence["records"]]
        == expected,
        "evidence measured grid",
    )
    records = evidence["records"]
    provenance = evidence["provenance"]
    seeds = []
    point_identities = {}
    point_census = {}
    for record in records:
        model, arithmetic, families = PROFILES[record["stratum_id"]]
        require(
            (record["model"], record["arithmetic"], record["Nsite"], record["NQPFull"])
            == (model, arithmetic, 16, 8),
            f"evidence frozen profile: {record['stratum_id']}",
        )
        require(record["saved_source_count"] == 64 and record["mpi_world_size"] == 1, "evidence primary execution")
        require(record["timer_id"] == TIMER_ID, "evidence timer identity")
        for key in PROVENANCE_KEYS:
            require(record[key] == provenance[key], f"evidence provenance binding: {key}")
        seeds.append(record["timing_seed_sha256"])
        point_key = (record["stratum_id"], record["operator_count"])
        point_identity = (
            canonical_bytes(record["family_counts"]),
            record["input_sha256"],
        )
        if point_key in point_identities:
            require(point_identities[point_key] == point_identity, f"evidence point identity: {point_key}")
        else:
            point_identities[point_key] = point_identity
        require(sum(record["family_counts"].values()) == record["operator_count"], "evidence family/operator count")
        require(set(family for family, count in record["family_counts"].items() if count) <= set(families), "evidence allowed families")
        if record["operator_census_sha256"] is not None:
            if point_key in point_census:
                require(point_census[point_key] == record["operator_census_sha256"], f"evidence point census: {point_key}")
            else:
                point_census[point_key] = record["operator_census_sha256"]
        if record["censored"]:
            require(
                record["censor_reason"] is not None
                and record["unique_target_count"] is None
                and record["block_count"] == 0
                and record["artifact_payload_bytes"] == 0
                and record["artifact_files"] == 1,
                "evidence censored reason/payload",
            )
        else:
            require(record["censor_reason"] is None, "evidence uncensored reason")
            require(record["operator_census_sha256"] is not None, "evidence uncensored census")
            require(0 <= record["unique_target_count"] <= record["operator_count"], "evidence target bound")
            if record["stratum_id"] == "SC-LOW-REUSE":
                require(10 * record["unique_target_count"] >= 9 * record["operator_count"], "evidence low-reuse minimum")
            require(record["block_count"] >= 32, "evidence block count")
            require(record["artifact_payload_bytes"] == 80 * record["block_count"] * record["operator_count"], "evidence payload formula")
            require(math.isclose(record["wall_seconds"], record["setup_seconds"] + record["active_seconds"], rel_tol=1e-12, abs_tol=1e-14), "evidence wall decomposition")
            require(record["active_seconds"] > 0 and record["wall_seconds"] > 0, "evidence active timing")
            require(record["artifact_files"] == 3, "evidence artifact file count")
        require(
            record["artifact_bytes"]
            == record["artifact_payload_bytes"]
            + record["artifact_metadata_bytes"],
            "evidence artifact decomposition",
        )
    require(len(seeds) == len(set(seeds)), "evidence timing-seed uniqueness")
    require([item["stratum_id"] for item in evidence["summaries"]] == STRATA, "evidence summary order")
    recomputed_summaries = [point_summary(records, stratum) for stratum in STRATA]
    for actual, recomputed in zip(evidence["summaries"], recomputed_summaries):
        require(actual["stratum_id"] == recomputed["stratum_id"] and actual["passed"] == recomputed["passed"], "evidence summary mapping")
        for actual_point, recomputed_point in zip(actual["point_uppers"], recomputed["point_uppers"]):
            require(actual_point["operator_count"] == recomputed_point["operator_count"], "evidence point count")
            left = actual_point["active_seconds_per_saved_source_upper"]
            right = recomputed_point["active_seconds_per_saved_source_upper"]
            require((left is None and right is None) or (left is not None and right is not None and math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-14)), "evidence point upper")
        for key in ["slope_upper_seconds_per_source_operator", "intercept_upper_seconds_per_source", "setup_upper_seconds"]:
            left = actual[key]
            right = recomputed[key]
            require((left is None and right is None) or (left is not None and right is not None and math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-14)), f"evidence summary {key}")
        require(actual["rss_minus_allocation_headroom_upper_bytes"] == recomputed["rss_minus_allocation_headroom_upper_bytes"], "evidence RSS headroom")
    expected_parity = [(stratum, size) for stratum in STRATA for size in [1, 2, 4]]
    parity = evidence["parity_records"]
    require([(item["stratum_id"], item["mpi_world_size"]) for item in parity] == expected_parity, "evidence parity grid")
    parity_pass = True
    expected_parity_ids = [
        f"parity-{stratum}-R512-rep0-mpi{size}"
        for stratum in STRATA
        for size in [1, 2, 4]
    ]
    require(
        [item["record_id"] for item in parity] == expected_parity_ids,
        "evidence parity ID binding",
    )
    for record in parity:
        require(
            record["censored"] == (record["censor_reason"] is not None),
            "evidence parity censor mapping",
        )
        require(
            record["artifact_bytes"]
            == record["artifact_payload_bytes"]
            + record["artifact_metadata_bytes"],
            "evidence parity artifact decomposition",
        )
        if record["censored"]:
            require(
                record["block_count"] == 0
                and record["artifact_payload_bytes"] == 0
                and record["artifact_files"] == 1,
                "evidence censored parity artifacts",
            )
        else:
            require(
                record["block_count"] >= 32
                and record["artifact_payload_bytes"]
                == 80 * record["block_count"] * 512
                and record["artifact_files"] == 3,
                "evidence parity artifacts",
            )
    for stratum in STRATA:
        group = [item for item in parity if item["stratum_id"] == stratum]
        group_pass = all(not item["censored"] and item["censor_reason"] is None for item in group) and len({item["normalized_semantic_payload_sha256"] for item in group}) == 1
        parity_pass &= group_pass
    artifact_records = warmup + records + parity
    artifact_bytes = sum(item["artifact_bytes"] for item in artifact_records)
    artifact_files = sum(item["artifact_files"] for item in artifact_records)
    require(evidence["campaign_artifact_bytes"] == artifact_bytes, "evidence campaign bytes")
    require(evidence["campaign_artifact_files"] == artifact_files, "evidence campaign files")
    warmup_pass = all(not item["censored"] for item in warmup)
    expected_pass = warmup_pass and all(item["passed"] for item in recomputed_summaries) and parity_pass and artifact_bytes <= 107374182400 and artifact_files <= 25000
    require(
        evidence["passed"] == expected_pass
        and evidence["passed"] == (evidence["status"] == "PASS"),
        "evidence status mapping",
    )


def build_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser("generate-inputs")
    generate.add_argument("--output-dir", type=Path, required=True)
    snapshot = subparsers.add_parser("snapshot")
    snapshot.add_argument("--repo", type=Path)
    snapshot.add_argument("--driver", type=Path, required=True)
    snapshot.add_argument("--build-dir", type=Path, required=True)
    snapshot.add_argument("--output", type=Path, required=True)
    snapshot.add_argument("--source-commit")
    snapshot.add_argument("--source-diff-sha256")
    run = subparsers.add_parser("run-one")
    run.add_argument("--driver", type=Path, required=True)
    run.add_argument("--mpiexec", type=Path)
    run.add_argument("--mpi-world-size", type=int, choices=[1, 2, 4], default=1)
    run.add_argument("--input-root", type=Path, required=True)
    run.add_argument("--identity", type=Path, required=True)
    run.add_argument("--raw-dir", type=Path, required=True)
    run.add_argument("--stratum", choices=STRATA, required=True)
    run.add_argument("--operator-count", type=int, choices=COUNTS, required=True)
    run.add_argument("--saved-sources", type=int, default=64)
    run.add_argument("--replicate", type=int, choices=range(7), required=True)
    run.add_argument("--kind", choices=["warmup", "measured", "parity", "selftest"], required=True)
    run.add_argument("--timeout", type=float, default=3600.0)
    validate = subparsers.add_parser("validate-raw")
    validate.add_argument("--record", type=Path, required=True)
    validate.add_argument("--raw-dir", type=Path, required=True)
    assemble_parser = subparsers.add_parser("assemble")
    assemble_parser.add_argument("--raw-dir", type=Path, required=True)
    assemble_parser.add_argument("--output", type=Path, required=True)
    assemble_parser.add_argument("--allow-official-output", action="store_true")
    evidence_parser = subparsers.add_parser("validate-evidence")
    evidence_parser.add_argument("--evidence", type=Path, required=True)
    seal_parser = subparsers.add_parser("seal-missing")
    seal_parser.add_argument("--input-root", type=Path, required=True)
    seal_parser.add_argument("--identity", type=Path, required=True)
    seal_parser.add_argument("--raw-dir", type=Path, required=True)
    seal_parser.add_argument(
        "--reason", choices=["scheduler_eviction"], required=True
    )
    subparsers.add_parser("self-test")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.command == "generate-inputs":
            generate_inputs(args.output_dir)
            print(f"PASS generated 16 scaling input points: {args.output_dir}")
        elif args.command == "snapshot":
            create_snapshot(
                args.repo,
                args.driver,
                args.build_dir,
                args.output,
                args.source_commit,
                args.source_diff_sha256,
            )
            print(f"PASS scaling identity snapshot: {args.output}")
        elif args.command == "run-one":
            run_one(args)
        elif args.command == "validate-raw":
            validate_raw_record(args.record, args.raw_dir)
            print(f"PASS raw scaling record: {args.record}")
        elif args.command == "assemble":
            assemble(args.raw_dir, args.output, args.allow_official_output)
        elif args.command == "validate-evidence":
            validate_evidence(load_json(args.evidence))
            print(f"PASS scaling evidence: {args.evidence}")
        elif args.command == "seal-missing":
            created = seal_missing_records(
                args.input_root, args.identity, args.raw_dir, args.reason
            )
            print(f"PASS sealed missing scaling records: created={created}")
        elif args.command == "self-test":
            self_test()
        return 0
    except (CampaignError, OSError, subprocess.SubprocessError, json.JSONDecodeError) as exc:
        print(f"FAIL P6-C2 scaling campaign: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
