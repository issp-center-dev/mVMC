#!/usr/bin/env python3

import json
import pathlib
import tempfile

import prepare_power_lanczos_p4f_v4_secondary as preparer
import runtest_power_lanczos_p4f_v3_result_validator as archive_fixture
import runtest_power_lanczos_p4f_v4_result_validator as result_fixture
import validate_power_lanczos_p4f_v3_result as common


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError(f"negative v4 preparer fixture passed: {label}")


def main():
    scratch = pathlib.Path.cwd() / "tmp"
    scratch.mkdir(exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(
                prefix="p4f-v4-preparer-", dir=scratch) as temporary:
            temporary = pathlib.Path(temporary)
            root, expected = result_fixture.build_fixture(temporary / "data")
            archive = temporary / "valid.tar.gz"
            checksum = archive_fixture.pack(root, archive)
            expected_path = temporary / "expected.json"
            expected_path.write_text(
                json.dumps(expected, indent=2) + "\n", encoding="utf-8")
            output = temporary / "prepared"
            report = preparer.prepare_secondary_inputs(
                archive, checksum, expected_path, output)
            if report["primary_decision"] != "GO" or \
                    report["block_stability_diagnostic_flag_count"] != 1 or \
                    report["expected_provenance_sha256"] != \
                    common.sha256_path(expected_path) or \
                    not (output / "raw_traces").is_dir():
                raise AssertionError("valid v4 preparation fixture failed")
            expect_failure(
                lambda: preparer.prepare_secondary_inputs(
                    archive, checksum, expected_path, output),
                "existing output")
            bad_expected = temporary / "bad-expected.json"
            bad_expected.write_text("{}\n", encoding="ascii")
            expect_failure(
                lambda: preparer.prepare_secondary_inputs(
                    archive, checksum, bad_expected,
                    temporary / "bad-prepared"),
                "invalid expected provenance")
    finally:
        try:
            scratch.rmdir()
        except OSError:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
