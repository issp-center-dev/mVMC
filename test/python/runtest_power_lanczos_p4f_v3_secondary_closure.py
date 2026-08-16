#!/usr/bin/env python3

import json
import os
import pathlib
import tempfile
import textwrap

import run_power_lanczos_p4f_v3_secondary_closure as closure
import runtest_power_lanczos_p4f_v3_result_validator as fixture


def expect_failure(function, label):
    try:
        function()
    except AssertionError:
        return
    raise AssertionError(f"negative secondary closure accepted: {label}")


def write_fake_generator(path, mismatch=False):
    path.write_text(textwrap.dedent(f"""\
        #!/usr/bin/env python3
        import argparse
        import hashlib
        import json
        import pathlib
        import shutil

        MISMATCH = {mismatch!r}

        def sha256_path(path):
            digest = hashlib.sha256()
            with pathlib.Path(path).open("rb") as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(chunk)
            return digest.hexdigest()

        parser = argparse.ArgumentParser()
        parser.add_argument("--trace-dir", required=True, type=pathlib.Path)
        parser.add_argument("--reparse-driver", required=True,
                            type=pathlib.Path)
        parser.add_argument("--solver-policy", required=True,
                            type=pathlib.Path)
        parser.add_argument("--confirmatory-policy", required=True,
                            type=pathlib.Path)
        parser.add_argument("--output-dir", required=True, type=pathlib.Path)
        args = parser.parse_args()

        prepared = args.trace_dir.parent
        evidence = json.loads((prepared /
            "primary/p4f_v3_confirmatory_evidence.json").read_text())
        decision = json.loads((prepared /
            "primary/p4f_v3_confirmatory_decision.json").read_text())
        metadata = dict(evidence["metadata"])
        metadata["generated_at"] = "2026-08-16 15:00 JST"
        metadata["reparse_driver"] = args.reparse_driver.name
        metadata["reparse_driver_sha256"] = sha256_path(
            args.reparse_driver)
        metadata["generator"] = pathlib.Path(__file__).name
        metadata["generator_sha256"] = sha256_path(pathlib.Path(__file__))
        evidence["metadata"] = metadata
        args.output_dir.mkdir()
        raw_output = args.output_dir / "raw_reparse"
        raw_output.mkdir()
        for case in evidence["cases"]:
            relative = pathlib.PurePosixPath(case["reparse_output"])
            output = args.output_dir / relative
            value = ("secondary " + case["trace"] + "\\n").encode()
            output.write_bytes(value)
            case["reparse_output_sha256"] = hashlib.sha256(value).hexdigest()
        if MISMATCH:
            evidence["cases"][0]["dimensions"][0]["prefixes"][0][
                "energy_se"] *= 2.0
        (args.output_dir / "p4f_v3_confirmatory_evidence.json").write_text(
            json.dumps(evidence, indent=2, allow_nan=False) + "\\n")
        (args.output_dir / "p4f_v3_confirmatory_decision.json").write_text(
            json.dumps(decision, indent=2, allow_nan=False) + "\\n")
        (args.output_dir / "metadata.json").write_text(
            json.dumps(metadata, indent=2, allow_nan=False) + "\\n")
        shutil.copyfile(args.solver_policy,
                        args.output_dir / "solver_policy.json")
        shutil.copyfile(args.confirmatory_policy,
                        args.output_dir / "confirmatory_policy.json")
        paths = sorted(item for item in args.output_dir.rglob("*")
                       if item.is_file() and item.name != "checksums.sha256")
        (args.output_dir / "checksums.sha256").write_text("".join(
            f"{{sha256_path(item)}}  "
            f"{{item.relative_to(args.output_dir).as_posix()}}\\n"
            for item in paths))
        raise SystemExit(0 if decision["decision"] == "GO" else 1)
    """), encoding="utf-8")
    return closure.validator.sha256_path(path)


def patch_primary_provenance(result, expected, generator_path,
                             generator_sha, driver, driver_sha):
    confirmatory = result / "confirmatory"
    metadata_path = confirmatory / "metadata.json"
    evidence_path = confirmatory / "p4f_v3_confirmatory_evidence.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata.update({
        "schema": 1,
        "generated_at": "2026-08-16 14:00 JST",
        "model": "OpenAI GPT-5.6 Codex",
        "reparse_driver": driver.name,
        "reparse_driver_sha256": driver_sha,
        "generator": generator_path.name,
        "generator_sha256": generator_sha,
        "thread_environment": {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
        },
    })
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    evidence["metadata"] = metadata
    fixture.write_json(metadata_path, metadata)
    fixture.write_json(evidence_path, evidence)
    expected["generator_sha256"] = generator_sha
    fixture.checksum_tree(confirmatory,
                          confirmatory / "checksums.sha256",
                          "checksums.sha256")
    fixture.checksum_tree(result, result / "artifact-ledger.sha256",
                          "artifact-ledger.sha256")


def make_prepared(root, stop, generator_path, generator_sha,
                  driver, driver_sha):
    result, expected = fixture.build_fixture(root / "fixture")
    if stop:
        fixture.convert_to_valid_stop(result)
    patch_primary_provenance(
        result, expected, generator_path, generator_sha, driver, driver_sha)
    archive = root / "result.tar.gz"
    checksum = fixture.pack(result, archive)
    prepared = root / "prepared"
    report = closure.preparer.prepare_secondary_inputs(
        archive, checksum, prepared, expected=expected)
    return prepared, report


def configure_runner(generator_path, generator_sha, preparation):
    closure.GENERATOR_PATH = generator_path.resolve()
    closure.GENERATOR_SHA256 = generator_sha
    closure.generator.SOLVER_POLICY_SHA256 = \
        preparation["solver_policy_sha256"]
    closure.generator.CONFIRMATORY_POLICY_SHA256 = \
        preparation["confirmatory_policy_sha256"]


def run_case(root, stop, mismatch=False):
    generator_path = root / (
        "fake_generator_mismatch.py" if mismatch else "fake_generator.py")
    generator_sha = write_fake_generator(generator_path, mismatch=mismatch)
    driver = root / "p4f_matrix_reparse_driver"
    driver.write_text("#!/bin/sh\nexit 0\n", encoding="ascii")
    driver.chmod(0o755)
    driver_sha = closure.validator.sha256_path(driver)
    prepared, preparation = make_prepared(
        root, stop, generator_path, generator_sha, driver, driver_sha)
    configure_runner(generator_path, generator_sha, preparation)
    output = root / "closure"
    report = closure.run_secondary_closure(
        prepared, driver, output, driver_sha)
    expected = "MISMATCH" if mismatch else ("STOP" if stop else "GO")
    if report["outcome"] != expected or report["p5_authorized"]:
        raise AssertionError("secondary closure outcome mismatch")
    if (report["p4f_candidate_for_independent_review"] is not
            (expected == "GO")):
        raise AssertionError("secondary review candidate mismatch")
    stored = json.loads((output / "p4f_v3_secondary_closure.json").read_text())
    if stored["outcome"] != expected or \
            stored["closure_manifest_entry_count"] != 34:
        raise AssertionError("stored secondary closure mismatch")
    closure.verify_tree_manifest(output, "closure_manifest.sha256")
    expect_failure(
        lambda: closure.run_secondary_closure(
            prepared, driver, output, driver_sha),
        "existing closure output")
    expect_failure(
        lambda: closure.run_secondary_closure(
            prepared, driver, root / "bad-driver", "0" * 64),
        "reparse driver SHA mismatch")
    return prepared, driver, driver_sha


def main():
    scratch = pathlib.Path.cwd() / "tmp"
    scratch.mkdir(exist_ok=True)
    saved = (
        closure.GENERATOR_PATH,
        closure.GENERATOR_SHA256,
        closure.generator.SOLVER_POLICY_SHA256,
        closure.generator.CONFIRMATORY_POLICY_SHA256,
    )
    try:
        with tempfile.TemporaryDirectory(
                prefix="mvmc-p4f-closure-go-", dir=scratch) as tmp:
            prepared, driver, driver_sha = run_case(
                pathlib.Path(tmp), stop=False)
            original = closure.COMPARATOR_SHA256
            closure.COMPARATOR_SHA256 = "0" * 64
            try:
                expect_failure(
                    lambda: closure.run_secondary_closure(
                        prepared, driver, pathlib.Path(tmp) / "bad-tool",
                        driver_sha),
                    "comparator SHA mismatch")
            finally:
                closure.COMPARATOR_SHA256 = original
        with tempfile.TemporaryDirectory(
                prefix="mvmc-p4f-closure-stop-", dir=scratch) as tmp:
            run_case(pathlib.Path(tmp), stop=True)
        with tempfile.TemporaryDirectory(
                prefix="mvmc-p4f-closure-mismatch-", dir=scratch) as tmp:
            run_case(pathlib.Path(tmp), stop=False, mismatch=True)
    finally:
        (closure.GENERATOR_PATH,
         closure.GENERATOR_SHA256,
         closure.generator.SOLVER_POLICY_SHA256,
         closure.generator.CONFIRMATORY_POLICY_SHA256) = saved
        try:
            scratch.rmdir()
        except OSError:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
