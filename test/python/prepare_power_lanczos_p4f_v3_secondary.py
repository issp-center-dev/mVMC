#!/usr/bin/env python3

import argparse
import json
import pathlib
import shutil
import tempfile

import validate_power_lanczos_p4f_v3_result as validator


VALIDATOR_SHA256 = (
    "13db55a4ad45499e8ced3bece94500c9120b9ace2b501f755195b293ba317fde"
)
PRIMARY_MEMBERS = (
    "confirmatory/p4f_v3_confirmatory_evidence.json",
    "confirmatory/p4f_v3_confirmatory_decision.json",
    "confirmatory/metadata.json",
)
POLICY_MEMBERS = (
    "confirmatory/solver_policy.json",
    "confirmatory/confirmatory_policy.json",
)


def write_manifest(root):
    paths = sorted(path for path in root.rglob("*") if path.is_file() and
                   path.name != "extraction_manifest.sha256")
    manifest = root / "extraction_manifest.sha256"
    manifest.write_text("".join(
        f"{validator.sha256_path(path)}  {path.relative_to(root).as_posix()}\n"
        for path in paths
    ), encoding="ascii")
    return manifest, len(paths)


def prepare_secondary_inputs(archive_path, checksum_path, output_dir,
                             expected=None):
    archive_path = pathlib.Path(archive_path)
    checksum_path = pathlib.Path(checksum_path)
    output_dir = pathlib.Path(output_dir)
    validator_path = pathlib.Path(validator.__file__).resolve()
    validator_sha = validator.sha256_path(validator_path)
    if validator_sha != VALIDATOR_SHA256:
        raise AssertionError("unexpected P4-F v3 result validator SHA-256")
    if output_dir.exists() or output_dir.is_symlink():
        raise AssertionError("secondary preparation output already exists")
    parent = output_dir.parent
    if not parent.is_dir() or parent.is_symlink():
        raise AssertionError("secondary preparation parent must be a directory")

    validation = validator.validate_result(
        archive_path, checksum_path, expected=expected)
    view = validator.ArchiveView(archive_path)
    if validator.sha256_path(archive_path) != validation["archive_sha256"]:
        raise AssertionError("result archive changed before preparation")
    raw_members = sorted(path for path in view.digests
                         if path.startswith("raw_traces/") and
                         path.endswith(".log"))
    if len(raw_members) != 24 or len({pathlib.PurePosixPath(path).name
                                     for path in raw_members}) != 24:
        raise AssertionError("validated raw trace name census mismatch")

    staging = pathlib.Path(tempfile.mkdtemp(
        prefix=f".{output_dir.name}.prepare.", dir=parent))
    try:
        raw_dir = staging / "raw_traces"
        primary_dir = staging / "primary"
        policy_dir = staging / "frozen_inputs"
        raw_dir.mkdir()
        primary_dir.mkdir()
        policy_dir.mkdir()
        for member in raw_members:
            name = pathlib.PurePosixPath(member).name
            (raw_dir / name).write_bytes(view.read_bytes(member))
        for member in PRIMARY_MEMBERS:
            (primary_dir / pathlib.PurePosixPath(member).name).write_bytes(
                view.read_bytes(member))
        for member in POLICY_MEMBERS:
            (policy_dir / pathlib.PurePosixPath(member).name).write_bytes(
                view.read_bytes(member))

        report = {
            "schema": 1,
            "valid": True,
            "archive": archive_path.name,
            "archive_sha256": validation["archive_sha256"],
            "validator": validator_path.name,
            "validator_sha256": validator_sha,
            "solver_policy_sha256": view.digests[POLICY_MEMBERS[0]],
            "confirmatory_policy_sha256": view.digests[POLICY_MEMBERS[1]],
            "primary_decision": validation["decision"],
            "trace_count": validation["trace_count"],
            "mandatory_trace_stop_count":
                validation["mandatory_trace_stop_count"],
            "pooled_gate_stop_count": validation["pooled_gate_stop_count"],
            "prepared_paths": {
                "raw_traces": "raw_traces",
                "primary_evidence":
                    "primary/p4f_v3_confirmatory_evidence.json",
                "primary_decision":
                    "primary/p4f_v3_confirmatory_decision.json",
                "primary_metadata": "primary/metadata.json",
                "solver_policy": "frozen_inputs/solver_policy.json",
                "confirmatory_policy":
                    "frozen_inputs/confirmatory_policy.json",
            },
        }
        (staging / "preparation.json").write_text(
            json.dumps(report, indent=2, allow_nan=False) + "\n",
            encoding="utf-8")
        manifest, file_count = write_manifest(staging)
        report["manifest"] = manifest.name
        report["manifest_entry_count"] = file_count
        (staging / "preparation.json").write_text(
            json.dumps(report, indent=2, allow_nan=False) + "\n",
            encoding="utf-8")
        write_manifest(staging)
        staging.rename(output_dir)
        return report
    except BaseException:
        if staging.is_dir() and not staging.is_symlink():
            shutil.rmtree(staging)
        raise


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True, type=pathlib.Path)
    parser.add_argument("--checksum", required=True, type=pathlib.Path)
    parser.add_argument("--output-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()
    report = prepare_secondary_inputs(
        args.archive, args.checksum, args.output_dir)
    print(json.dumps(report, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
