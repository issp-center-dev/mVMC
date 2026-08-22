#!/usr/bin/env python3
"""Focused serial/MPI/censor/mutation test for the P6-C2 scaling infra."""

from __future__ import annotations

import argparse
import copy
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def run(command, *, expect_success=True):
    completed = subprocess.run(
        [str(item) for item in command],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        env={
            **os.environ,
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
        },
    )
    if (completed.returncode == 0) != expect_success:
        sys.stderr.buffer.write(completed.stdout)
        sys.stderr.buffer.write(completed.stderr)
        raise RuntimeError(
            f"unexpected exit {completed.returncode}: {' '.join(str(item) for item in command)}"
        )
    return completed


def campaign(python, script, *arguments, expect_success=True):
    return run([python, script, *arguments], expect_success=expect_success)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--driver", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--mpiexec", type=Path)
    args = parser.parse_args()
    local_tmp = Path.cwd() / "tmp"
    local_tmp.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="p6c2-scaling-test.", dir=local_tmp))
    try:
        inputs = work / "inputs"
        identity = work / "identity.json"
        raw = work / "raw"
        campaign(args.python, args.campaign, "generate-inputs", "--output-dir", inputs)
        snapshot_source = ["--repo", args.repo]
        if not (args.repo / ".git").exists():
            source_commit = os.environ.get("MVMC_P6C2_SOURCE_COMMIT")
            source_diff = os.environ.get("MVMC_P6C2_SOURCE_DIFF_SHA256")
            if source_commit is None or source_diff is None:
                raise RuntimeError(
                    "archive test requires MVMC_P6C2 source identity"
                )
            snapshot_source = [
                "--source-commit",
                source_commit,
                "--source-diff-sha256",
                source_diff,
            ]
        campaign(
            args.python,
            args.campaign,
            "snapshot",
            *snapshot_source,
            "--driver",
            args.driver,
            "--build-dir",
            args.build_dir,
            "--output",
            identity,
        )
        for stratum in ["SC-ONEBODY", "SC-QUARTIC", "SC-MIXED", "SC-LOW-REUSE"]:
            campaign(
                args.python,
                args.campaign,
                "run-one",
                "--driver",
                args.driver,
                "--input-root",
                inputs,
                "--identity",
                identity,
                "--raw-dir",
                raw,
                "--stratum",
                stratum,
                "--operator-count",
                "32",
                "--saved-sources",
                "2",
                "--replicate",
                "0",
                "--kind",
                "selftest",
                "--mpi-world-size",
                "1",
                "--timeout",
                "120",
            )
        campaign(
            args.python,
            args.campaign,
            "run-one",
            "--driver",
            args.driver,
            "--input-root",
            inputs,
            "--identity",
            identity,
            "--raw-dir",
            raw,
            "--stratum",
            "SC-LOW-REUSE",
            "--operator-count",
            "512",
            "--saved-sources",
            "2",
            "--replicate",
            "0",
            "--kind",
            "selftest",
            "--mpi-world-size",
            "1",
            "--timeout",
            "120",
        )
        if args.mpiexec is not None:
            for size in [2, 4]:
                campaign(
                    args.python,
                    args.campaign,
                    "run-one",
                    "--driver",
                    args.driver,
                    "--mpiexec",
                    args.mpiexec,
                    "--input-root",
                    inputs,
                    "--identity",
                    identity,
                    "--raw-dir",
                    raw,
                    "--stratum",
                    "SC-LOW-REUSE",
                    "--operator-count",
                    "32",
                    "--saved-sources",
                    "2",
                    "--replicate",
                    "0",
                    "--kind",
                    "selftest",
                    "--mpi-world-size",
                    str(size),
                    "--timeout",
                    "120",
                )
            hashes = []
            for size in [1, 2, 4]:
                record = json.loads(
                    (
                        raw
                        / f"selftest-SC-LOW-REUSE-R032-rep0-mpi{size}.json"
                    ).read_text()
                )
                hashes.append(record["metrics"]["normalized_semantic_payload_sha256"])
            if len(set(hashes)) != 1:
                raise RuntimeError("MPI 1/2/4 semantic payload mismatch")

        source_record = json.loads(
            (raw / "selftest-SC-LOW-REUSE-R032-rep0-mpi1.json").read_text()
        )
        mutations = work / "mutations"
        mutations.mkdir()
        bad_model = copy.deepcopy(source_record)
        bad_model["metrics"]["model"] = "pure_spin"
        bad_model_path = mutations / "bad-model.json"
        bad_model_path.write_text(json.dumps(bad_model) + "\n")
        campaign(
            args.python,
            args.campaign,
            "validate-raw",
            "--record",
            bad_model_path,
            "--raw-dir",
            raw,
            expect_success=False,
        )
        bad_reuse = copy.deepcopy(source_record)
        bad_reuse["metrics"]["minimum_unique_target_count"] = 1
        bad_reuse_path = mutations / "bad-reuse.json"
        bad_reuse_path.write_text(json.dumps(bad_reuse) + "\n")
        campaign(
            args.python,
            args.campaign,
            "validate-raw",
            "--record",
            bad_reuse_path,
            "--raw-dir",
            raw,
            expect_success=False,
        )

        censor_raw = work / "censor-raw"
        campaign(
            args.python,
            args.campaign,
            "run-one",
            "--driver",
            args.driver,
            "--input-root",
            inputs,
            "--identity",
            identity,
            "--raw-dir",
            censor_raw,
            "--stratum",
            "SC-LOW-REUSE",
            "--operator-count",
            "512",
            "--saved-sources",
            "64",
            "--replicate",
            "1",
            "--kind",
            "selftest",
            "--mpi-world-size",
            "1",
            "--timeout",
            "0.000001",
        )
        censored = json.loads(
            (censor_raw / "selftest-SC-LOW-REUSE-R512-rep1-mpi1.json").read_text()
        )
        if not censored["process"]["censored"] or censored["metrics"] is not None:
            raise RuntimeError("timeout was not retained as a censored no-payload record")

        empty_driver = work / "power_lanczos_observable_scaling_driver"
        empty_driver.write_text("#!/bin/sh\nexit 0\n", encoding="ascii")
        empty_driver.chmod(0o755)
        missing_raw = work / "missing-raw"
        campaign(
            args.python,
            args.campaign,
            "run-one",
            "--driver",
            empty_driver,
            "--input-root",
            inputs,
            "--identity",
            identity,
            "--raw-dir",
            missing_raw,
            "--stratum",
            "SC-ONEBODY",
            "--operator-count",
            "1",
            "--saved-sources",
            "64",
            "--replicate",
            "2",
            "--kind",
            "selftest",
            "--mpi-world-size",
            "1",
        )
        missing = json.loads(
            (missing_raw / "selftest-SC-ONEBODY-R001-rep2-mpi1.json").read_text()
        )
        if missing["process"]["censor_reason"] != "missing_metric":
            raise RuntimeError("missing metric was not retained as censored")

        sealed_raw = work / "sealed-raw"
        campaign(
            args.python,
            args.campaign,
            "seal-missing",
            "--input-root",
            inputs,
            "--identity",
            identity,
            "--raw-dir",
            sealed_raw,
            "--reason",
            "scheduler_eviction",
        )
        sealed = list(sealed_raw.glob("*.json"))
        if len(sealed) != 140:
            raise RuntimeError(f"scheduler seal grid mismatch: {len(sealed)}")
        if any(
            json.loads(path.read_text())["process"]["censor_reason"]
            != "scheduler_eviction"
            for path in sealed
        ):
            raise RuntimeError("scheduler seal reason mismatch")
        print(
            "PASS P6-C2 scaling driver serial=5 MPI_parity={} raw_mutations=2 censor=2 scheduler_seal=140".format(
                "1/2/4" if args.mpiexec is not None else "SKIP"
            )
        )
        return 0
    except (OSError, RuntimeError, subprocess.SubprocessError, json.JSONDecodeError) as exc:
        print(f"FAIL P6-C2 scaling driver test: {exc}", file=sys.stderr)
        return 1
    finally:
        resolved = work.resolve()
        if resolved.parent == local_tmp.resolve() and not resolved.is_symlink():
            shutil.rmtree(resolved)
        if local_tmp.exists() and not any(local_tmp.iterdir()):
            local_tmp.rmdir()


if __name__ == "__main__":
    raise SystemExit(main())
