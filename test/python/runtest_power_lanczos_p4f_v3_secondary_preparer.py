#!/usr/bin/env python3

import json
import pathlib
import tempfile

import prepare_power_lanczos_p4f_v3_secondary as preparer
import runtest_power_lanczos_p4f_v3_result_validator as validator_fixture


def verify_manifest(root):
    recorded = (root / "extraction_manifest.sha256").read_text(
        encoding="ascii").splitlines()
    expected = []
    for path in sorted(item for item in root.rglob("*") if item.is_file() and
                       item.name != "extraction_manifest.sha256"):
        expected.append(
            f"{validator_fixture.validator.sha256_path(path)}  "
            f"{path.relative_to(root).as_posix()}")
    if recorded != expected:
        raise AssertionError("secondary extraction manifest mismatch")


def run_fixture(root, stop):
    result, expected = validator_fixture.build_fixture(root)
    if stop:
        validator_fixture.convert_to_valid_stop(result)
    archive = root / "result.tar.gz"
    checksum = validator_fixture.pack(result, archive)
    output = root / "prepared"
    report = preparer.prepare_secondary_inputs(
        archive, checksum, output, expected=expected)
    expected_decision = "STOP" if stop else "GO"
    if report["primary_decision"] != expected_decision or (
        report["trace_count"] != 24
    ):
        raise AssertionError("secondary preparation decision mismatch")
    if len(list((output / "raw_traces").glob("*.log"))) != 24:
        raise AssertionError("secondary raw trace census mismatch")
    if not all((output / path).is_file() for path in (
        "primary/p4f_v3_confirmatory_evidence.json",
        "primary/p4f_v3_confirmatory_decision.json",
        "primary/metadata.json",
        "frozen_inputs/solver_policy.json",
        "frozen_inputs/confirmatory_policy.json",
        "preparation.json",
        "extraction_manifest.sha256",
    )):
        raise AssertionError("secondary preparation file missing")
    stored = json.loads((output / "preparation.json").read_text())
    if stored != report or stored["manifest_entry_count"] != 30:
        raise AssertionError("secondary preparation report mismatch")
    verify_manifest(output)
    if not stop:
        validator_sha = preparer.VALIDATOR_SHA256
        try:
            preparer.VALIDATOR_SHA256 = "0" * 64
            try:
                preparer.prepare_secondary_inputs(
                    archive, checksum, root / "wrong-validator",
                    expected=expected)
            except AssertionError:
                pass
            else:
                raise AssertionError("validator SHA mismatch was accepted")
        finally:
            preparer.VALIDATOR_SHA256 = validator_sha
    try:
        preparer.prepare_secondary_inputs(
            archive, checksum, output, expected=expected)
    except AssertionError:
        pass
    else:
        raise AssertionError("existing secondary output was overwritten")


def main():
    scratch = pathlib.Path.cwd() / "tmp"
    scratch.mkdir(exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(
                prefix="mvmc-p4f-secondary-go-", dir=scratch) as tmp:
            run_fixture(pathlib.Path(tmp), stop=False)
        with tempfile.TemporaryDirectory(
                prefix="mvmc-p4f-secondary-stop-", dir=scratch) as tmp:
            run_fixture(pathlib.Path(tmp), stop=True)
    finally:
        try:
            scratch.rmdir()
        except OSError:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
