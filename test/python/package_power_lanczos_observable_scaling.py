#!/usr/bin/env python3
"""Create and validate a result-unobserved P6-C2 genkai package."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import posixpath
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path, PurePosixPath


SCHEMA_VERSION = 1
PACKAGE_PREFIX = "p6c2-observable-scaling"
SOURCE_DIFF_SHA256 = hashlib.sha256(b"tracked-diff\0").hexdigest()
EXPECTED_GRID = {"warmup": 16, "measured": 112, "parity": 12}
PYTHON_MODULE = "python/3.12.11"
PYTHON_EXECUTABLE = "python3.12"
PYTHON_VERSION = "3.12.11"
PYTHON_WHEEL_TARGET = "CPython 3.12 / manylinux2014 x86_64"
PYTHON_WHEELS = (
    {
        "requirement": "attrs==26.1.0",
        "filename": "attrs-26.1.0-py3-none-any.whl",
        "sha256": "c647aa4a12dfbad9333ca4e71fe62ddc36f4e63b2d260a37a8b83d2f043ac309",
    },
    {
        "requirement": "jsonschema==4.26.0",
        "filename": "jsonschema-4.26.0-py3-none-any.whl",
        "sha256": "d489f15263b8d200f8387e64b4c3a75f06629559fb73deb8fdfb525f2dab50ce",
    },
    {
        "requirement": "jsonschema-specifications==2025.9.1",
        "filename": "jsonschema_specifications-2025.9.1-py3-none-any.whl",
        "sha256": "98802fee3a11ee76ecaca44429fda8a41bff98b00a0f2838151b113f210cc6fe",
    },
    {
        "requirement": "referencing==0.37.0",
        "filename": "referencing-0.37.0-py3-none-any.whl",
        "sha256": "381329a9f99628c9069361716891d34ad94af76e461dcb0335825aecc7692231",
    },
    {
        "requirement": "rpds-py==2026.6.3",
        "filename": "rpds_py-2026.6.3-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl",
        "sha256": "ecabd69db66de867690f9797f2f8fa27ba501bbc24540cbdbdc649cd15888ba6",
    },
    {
        "requirement": "typing-extensions==4.16.0",
        "filename": "typing_extensions-4.16.0-py3-none-any.whl",
        "sha256": "481caa481374e813c1b176ada14e97f1f67a4539ce9cfeb3f350d78d6370c2e8",
    },
)


class PackageError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PackageError(message)


def sha_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def python_requirements_lock() -> str:
    return "".join(
        f"{wheel['requirement']} --hash=sha256:{wheel['sha256']}\n"
        for wheel in PYTHON_WHEELS
    )


def validate_python_wheelhouse(path: Path) -> list[Path]:
    require(path.is_dir() and not path.is_symlink(), "Python wheelhouse directory")
    expected = {wheel["filename"]: wheel for wheel in PYTHON_WHEELS}
    observed = {item.name: item for item in path.iterdir()}
    require(set(observed) == set(expected), "Python wheelhouse exact member set")
    result = []
    for filename, wheel in sorted(expected.items()):
        source = observed[filename]
        require(
            source.is_file()
            and not source.is_symlink()
            and sha_file(source) == wheel["sha256"],
            f"Python wheel: {filename}",
        )
        result.append(source)
    return result


def git(repo: Path, *arguments: str) -> bytes:
    return subprocess.check_output(["git", "-C", str(repo), *arguments])


def safe_name(name: str) -> bool:
    path = PurePosixPath(name)
    return (
        name != ""
        and not path.is_absolute()
        and ".." not in path.parts
        and not any(ord(character) < 32 for character in name)
    )


def git_tree_entries(repo: Path, commit: str, prefix: str):
    raw = git(repo, "ls-tree", "-rz", "--full-tree", commit)
    entries = []
    for row in raw.split(b"\0"):
        if not row:
            continue
        header, raw_name = row.split(b"\t", 1)
        mode, kind, object_id = header.decode("ascii").split()
        name = raw_name.decode("utf-8")
        require(safe_name(name), f"unsafe Git tree path: {name!r}")
        if kind == "commit":
            continue
        require(kind == "blob", f"unsupported Git object type: {kind}")
        payload = git(repo, "cat-file", "blob", object_id)
        archive_name = posixpath.join(prefix, name)
        if mode == "120000":
            linkname = payload.decode("utf-8")
            resolved_link = posixpath.normpath(
                posixpath.join(posixpath.dirname(archive_name), linkname)
            )
            require(
                not linkname.startswith("/")
                and resolved_link.startswith("mVMC/"),
                f"unsafe symlink target: {archive_name}",
            )
            entries.append((archive_name, None, 0o777, linkname))
        else:
            require(mode in {"100644", "100755"}, f"unsupported mode: {mode}")
            entries.append(
                (archive_name, payload, 0o755 if mode == "100755" else 0o644, None)
            )
    return entries


def submodule_commits(repo: Path, commit: str):
    raw = git(repo, "ls-tree", "-rz", "--full-tree", commit)
    result = {}
    for row in raw.split(b"\0"):
        if not row:
            continue
        header, raw_name = row.split(b"\t", 1)
        mode, kind, object_id = header.decode("ascii").split()
        if mode == "160000" and kind == "commit":
            name = raw_name.decode("utf-8")
            require(safe_name(name), f"unsafe submodule path: {name!r}")
            result[name] = object_id
    return result


def write_tar_gz(path: Path, entries) -> None:
    seen = set()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as raw_stream:
        with gzip.GzipFile(
            filename="", mode="wb", compresslevel=9, mtime=0, fileobj=raw_stream
        ) as gzip_stream:
            with tarfile.open(
                fileobj=gzip_stream, mode="w", format=tarfile.GNU_FORMAT
            ) as archive:
                for name, payload, mode, linkname in sorted(entries):
                    require(name not in seen and safe_name(name), f"duplicate path: {name}")
                    seen.add(name)
                    info = tarfile.TarInfo(name)
                    info.uid = 0
                    info.gid = 0
                    info.uname = ""
                    info.gname = ""
                    info.mtime = 0
                    info.mode = mode
                    if linkname is None:
                        require(payload is not None, f"missing payload: {name}")
                        info.size = len(payload)
                        archive.addfile(info, io.BytesIO(payload))
                    else:
                        info.type = tarfile.SYMTYPE
                        info.linkname = linkname
                        info.size = 0
                        archive.addfile(info)


def archive_file_entries(root: Path, paths):
    entries = []
    for relative, mode in paths:
        require(safe_name(relative), f"unsafe package path: {relative}")
        entries.append((relative, (root / relative).read_bytes(), mode, None))
    return entries


def job_script(
    source_name: str,
    source_commit: str,
    source_diff_sha256: str,
    project_code: str,
    resource_group: str,
) -> str:
    short = source_commit[:12]
    return f"""#!/usr/bin/env bash
#PJM -g {project_code}
#PJM -L rscgrp={resource_group}
#PJM -L node=1
#PJM --mpi proc=4
#PJM -L elapse=04:00:00
#PJM -j
#PJM -o p6c2.observable.scaling.{short}.log.%j

set -euo pipefail

RUN_DIR=$(pwd -P)
SOURCE_ARCHIVE=$RUN_DIR/{source_name}
SOURCE_CHECKSUM=$RUN_DIR/source_archive.sha256
PACKAGE_MANIFEST=$RUN_DIR/package_manifest.json
PYTHON_REQUIREMENTS=$RUN_DIR/python-requirements.lock
PYTHON_WHEELHOUSE=$RUN_DIR/python-wheelhouse
PYTHON_WHEELHOUSE_CHECKSUM=$RUN_DIR/python-wheelhouse.sha256
SCRATCH=
EVIDENCE_DIR=
INPUT_DIR=
IDENTITY=
CAMPAIGN=
BOOTSTRAP_PYTHON=
PYTHON=
PYTHON_ENV=
CURRENT_PHASE=bootstrap
FINALIZATION_ACTIVE=0
FINALIZATION_COUNT=0
FINALIZED=0
OUTPUT_ARCHIVE=$RUN_DIR/p6c2-observable-scaling-evidence-{short}.tar.gz
OUTPUT_CHECKSUM=$RUN_DIR/p6c2-observable-scaling-evidence-{short}.sha256
FAILURE_METADATA=$RUN_DIR/p6c2-observable-scaling-failure-{short}.json
FAILURE_CHECKSUM=$RUN_DIR/p6c2-observable-scaling-failure-{short}.sha256
SOURCE_COMMIT={source_commit}
SOURCE_DIFF_SHA256={source_diff_sha256}
export MVMC_P6C2_SOURCE_COMMIT=$SOURCE_COMMIT
export MVMC_P6C2_SOURCE_DIFF_SHA256=$SOURCE_DIFF_SHA256

cleanup() {{
  if [[ -n ${{SCRATCH:-}} && -d $SCRATCH && ! -L $SCRATCH ]]; then
    resolved=$(cd "$SCRATCH" && pwd -P)
    case "$resolved" in
      "$RUN_DIR"/tmp/*) rm -rf -- "$resolved" ;;
      *) echo "refusing unexpected scratch cleanup: $resolved" >&2 ;;
    esac
  fi
  rm -f -- "$OUTPUT_ARCHIVE.partial.$$" "$OUTPUT_CHECKSUM.partial.$$"
  rm -f -- "$FAILURE_METADATA.partial.$$" "$FAILURE_CHECKSUM.partial.$$"
  rmdir "$RUN_DIR/tmp" 2>/dev/null || true
}}

archive_evidence() {{
  local archive_tmp checksum_tmp digest
  [[ -n ${{EVIDENCE_DIR:-}} && -d $EVIDENCE_DIR ]] || return 1
  find "$EVIDENCE_DIR" -type l -print -quit | grep -q . && return 1
  (
    cd "$EVIDENCE_DIR"
    find . -type f ! -name artifact-ledger.sha256 -print0 | sort -z | \
      xargs -0 sha256sum > artifact-ledger.sha256
    sha256sum -c artifact-ledger.sha256
  )
  archive_tmp=$OUTPUT_ARCHIVE.partial.$$
  checksum_tmp=$OUTPUT_CHECKSUM.partial.$$
  rm -f -- "$archive_tmp" "$checksum_tmp"
  COPYFILE_DISABLE=1 tar -czf "$archive_tmp" -C "$EVIDENCE_DIR" .
  tar -tzf "$archive_tmp" >/dev/null
  digest=$(sha256sum "$archive_tmp" | awk '{{print $1}}')
  printf '%s  %s\n' "$digest" "$(basename "$OUTPUT_ARCHIVE")" > "$checksum_tmp"
  mv -f -- "$archive_tmp" "$OUTPUT_ARCHIVE"
  mv -f -- "$checksum_tmp" "$OUTPUT_CHECKSUM"
}}

write_failure_metadata() {{
  local failure_exit=$1 failure_reason=$2 metadata_tmp checksum_tmp digest
  metadata_tmp=$FAILURE_METADATA.partial.$$
  checksum_tmp=$FAILURE_CHECKSUM.partial.$$
  cat > "$metadata_tmp" <<EOF
{{
  "schema_version": 1,
  "status": "failed",
  "exit_code": $failure_exit,
  "phase": "$CURRENT_PHASE",
  "reason": "$failure_reason",
  "finalization_count": $FINALIZATION_COUNT,
  "source_commit": "$SOURCE_COMMIT",
  "source_diff_sha256": "$SOURCE_DIFF_SHA256"
}}
EOF
  digest=$(sha256sum "$metadata_tmp" | awk '{{print $1}}')
  printf '%s  %s\n' "$digest" "$(basename "$FAILURE_METADATA")" > "$checksum_tmp"
  mv -f -- "$metadata_tmp" "$FAILURE_METADATA"
  mv -f -- "$checksum_tmp" "$FAILURE_CHECKSUM"
  if [[ -n ${{EVIDENCE_DIR:-}} && -d $EVIDENCE_DIR ]]; then
    mkdir -p "$EVIDENCE_DIR/workflow"
    cp "$FAILURE_METADATA" "$FAILURE_CHECKSUM" "$EVIDENCE_DIR/workflow/"
  fi
}}

finalize_failure() {{
  local failure_exit=$1 failure_reason=$2 metadata_status archive_status
  if [[ $FINALIZED -eq 1 || $FINALIZATION_ACTIVE -eq 1 ]]; then
    return 0
  fi
  FINALIZATION_ACTIVE=1
  FINALIZATION_COUNT=$((FINALIZATION_COUNT + 1))
  set +e
  write_failure_metadata "$failure_exit" "$failure_reason"
  metadata_status=$?
  archive_evidence
  archive_status=$?
  {{
    echo "status=failed"
    echo "exit_code=$failure_exit"
    echo "phase=$CURRENT_PHASE"
    echo "reason=$failure_reason"
    echo "finalization_count=$FINALIZATION_COUNT"
    echo "failure_metadata_status=$metadata_status"
    echo "evidence_archive_status=$archive_status"
  }} > "$RUN_DIR/FAILED"
  FINALIZED=1
  FINALIZATION_ACTIVE=0
  return 0
}}

on_exit() {{
  local exit_status=$?
  trap - EXIT HUP INT TERM
  set +e
  if [[ $exit_status -ne 0 && $FINALIZED -eq 0 ]]; then
    finalize_failure "$exit_status" nonzero_exit
  fi
  cleanup
  exit "$exit_status"
}}

seal_on_signal() {{
  local signal_name=$1
  trap - HUP INT TERM
  set +e
  echo "scheduler signal received; sealing every absent registered record" >&2
  if [[ -n ${{CAMPAIGN:-}} && -f $CAMPAIGN && -n ${{INPUT_DIR:-}} && \
        -d $INPUT_DIR && -n ${{IDENTITY:-}} && -f $IDENTITY ]]; then
    "$PYTHON" "$CAMPAIGN" seal-missing \
      --input-root "$INPUT_DIR" --identity "$IDENTITY" \
      --raw-dir "$EVIDENCE_DIR/raw" --reason scheduler_eviction
    evidence=$EVIDENCE_DIR/power-lanczos-p6c2-observable-scaling-evidence-v2.json
    if [[ ! -e $evidence ]]; then
      "$PYTHON" "$CAMPAIGN" assemble \
        --raw-dir "$EVIDENCE_DIR/raw" --output "$evidence"
    fi
  fi
  finalize_failure 130 "scheduler_signal_$signal_name"
  exit 130
}}

prepare_scratch() {{
  mkdir -p "$RUN_DIR/tmp"
  SCRATCH=$(mktemp -d "$RUN_DIR/tmp/p6c2-scaling.XXXXXX")
  SOURCE_DIR=$SCRATCH/source
  BUILD_DIR=$SCRATCH/build
  EVIDENCE_DIR=$SCRATCH/evidence
  INPUT_DIR=$EVIDENCE_DIR/inputs
  IDENTITY=$EVIDENCE_DIR/identity.json
  mkdir "$SOURCE_DIR" "$BUILD_DIR" "$EVIDENCE_DIR"
}}

run_lifecycle_selftest() {{
  case "$1" in
    success)
      CURRENT_PHASE=selftest_success
      mkdir -p "$EVIDENCE_DIR/workflow"
      echo selftest-success > "$EVIDENCE_DIR/workflow/selftest.log"
      archive_evidence
      CURRENT_PHASE=selftest_success_marker
      touch "$RUN_DIR/COMPLETED"
      FINALIZED=1
      exit 0
      ;;
    nonzero)
      CURRENT_PHASE=selftest_nonzero
      echo injected-configure-failure > "$EVIDENCE_DIR/configure.log"
      exit 23
      ;;
    dependency)
      CURRENT_PHASE=selftest_python_dependency
      echo injected-python-dependency-failure > "$EVIDENCE_DIR/python-dependency.log"
      exit 24
      ;;
    signal)
      CURRENT_PHASE=selftest_signal
      echo injected-scheduler-signal > "$EVIDENCE_DIR/configure.log"
      kill -TERM "$$"
      exit 99
      ;;
    *)
      CURRENT_PHASE=selftest_unknown
      exit 64
      ;;
  esac
}}

trap on_exit EXIT
trap 'seal_on_signal HUP' HUP
trap 'seal_on_signal INT' INT
trap 'seal_on_signal TERM' TERM

CURRENT_PHASE=prepare_scratch
prepare_scratch
if [[ -n ${{MVMC_P6C2_JOB_SELFTEST_MODE:-}} ]]; then
  run_lifecycle_selftest "$MVMC_P6C2_JOB_SELFTEST_MODE"
fi

CURRENT_PHASE=package_manifest_checksum
test -f "$SOURCE_ARCHIVE" && test ! -L "$SOURCE_ARCHIVE"
(cd "$RUN_DIR" && sha256sum -c package_manifest.sha256)
CURRENT_PHASE=source_archive_checksum
(cd "$RUN_DIR" && sha256sum -c "$(basename "$SOURCE_CHECKSUM")")
CURRENT_PHASE=python_wheelhouse_checksum
: > "$EVIDENCE_DIR/python-dependency.log"
(cd "$RUN_DIR" && sha256sum -c python-wheelhouse.sha256) \
  >> "$EVIDENCE_DIR/python-dependency.log" 2>&1
CURRENT_PHASE=workflow_provenance
mkdir "$EVIDENCE_DIR/workflow"
cp "$PACKAGE_MANIFEST" "$RUN_DIR/package_manifest.sha256" \
  "$SOURCE_CHECKSUM" "$PYTHON_REQUIREMENTS" \
  "$PYTHON_WHEELHOUSE_CHECKSUM" \
  "$RUN_DIR/workflow/p6c2_scaling_genkai_job.sh" \
  "$EVIDENCE_DIR/workflow/"

CURRENT_PHASE=source_extraction
tar -xzf "$SOURCE_ARCHIVE" -C "$SOURCE_DIR"
SOURCE_ROOT=$SOURCE_DIR/mVMC

CURRENT_PHASE=environment_setup
module load intel/2023.2 mvapich/3.0-intel2023.2 {PYTHON_MODULE}
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PIP_CONFIG_FILE=/dev/null

BOOTSTRAP_PYTHON=$(command -v {PYTHON_EXECUTABLE})
MPIEXEC=$(command -v mpiexec)
CAMPAIGN=$SOURCE_ROOT/test/python/power_lanczos_observable_scaling_campaign.py
PACKAGER=$SOURCE_ROOT/test/python/package_power_lanczos_observable_scaling.py
CURRENT_PHASE=python_bootstrap_preflight
{{
  "$BOOTSTRAP_PYTHON" --version
  "$BOOTSTRAP_PYTHON" -I -c \
    'import ensurepip,sys,venv; expected=(3,12,11); assert sys.version_info[:3] == expected, (sys.version_info, expected)'
}} >> "$EVIDENCE_DIR/python-dependency.log" 2>&1
CURRENT_PHASE=package_validation
"$BOOTSTRAP_PYTHON" "$PACKAGER" validate --package-dir "$RUN_DIR" \
  --skip-package-archive
CURRENT_PHASE=python_environment_creation
PYTHON_ENV=$SCRATCH/python-env
"$BOOTSTRAP_PYTHON" -m venv "$PYTHON_ENV" \
  >> "$EVIDENCE_DIR/python-dependency.log" 2>&1
PYTHON=$PYTHON_ENV/bin/python
CURRENT_PHASE=python_dependency_install
"$PYTHON" -m pip --isolated install \
  --disable-pip-version-check --no-cache-dir --no-index --no-deps --no-compile \
  --only-binary=:all: --require-hashes \
  --find-links "$PYTHON_WHEELHOUSE" -r "$PYTHON_REQUIREMENTS" \
  >> "$EVIDENCE_DIR/python-dependency.log" 2>&1
CURRENT_PHASE=python_dependency_preflight
"$PYTHON" -I -c \
  'import importlib.metadata as metadata,sys; from jsonschema import Draft202012Validator; expected=(3,12,11); version=metadata.version("jsonschema"); assert sys.version_info[:3] == expected; assert version == "4.26.0"; Draft202012Validator.check_schema({{"type":"object"}}); print(f"python={{sys.version.split()[0]}} jsonschema={{version}} Draft202012Validator=OK")' \
  >> "$EVIDENCE_DIR/python-dependency.log" 2>&1

CURRENT_PHASE=configure
cmake -S "$SOURCE_ROOT" -B "$BUILD_DIR" \
  -DCONFIG=intel \
  -DTesting=ON \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DCMAKE_C_FLAGS='-fp-model precise' \
  -DCMAKE_CXX_FLAGS='-fp-model precise' \
  -DCMAKE_Fortran_FLAGS='-fp-model precise' \
  -DPYTHON_EXECUTABLE="$PYTHON" \
  > "$EVIDENCE_DIR/configure.log" 2>&1
CURRENT_PHASE=build
cmake --build "$BUILD_DIR" --parallel 1 \
  --target power_lanczos_observable_scaling_driver \
  > "$EVIDENCE_DIR/build.log" 2>&1

CURRENT_PHASE=focused_test
env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  BLIS_NUM_THREADS=1 ctest --test-dir "$BUILD_DIR" --output-on-failure \
  -R '^PowerLanczosObservableScaling(Campaign_Unit|Driver_Focused|Package_Lifecycle)$' \
  > "$EVIDENCE_DIR/focused-ctest.log" 2>&1

DRIVER=$BUILD_DIR/test/power_lanczos_observable_scaling_driver
CURRENT_PHASE=input_generation
"$PYTHON" "$CAMPAIGN" generate-inputs --output-dir "$INPUT_DIR"
CURRENT_PHASE=identity_snapshot
"$PYTHON" "$CAMPAIGN" snapshot \
  --driver "$DRIVER" --build-dir "$BUILD_DIR" --output "$IDENTITY" \
  --source-commit "$SOURCE_COMMIT" \
  --source-diff-sha256 "$SOURCE_DIFF_SHA256"

CURRENT_PHASE=workflow_metadata
mkdir "$EVIDENCE_DIR/raw"

{{
  echo "measurement=p6c2_observable_operator_count_scaling_v2"
  echo "source_commit=$SOURCE_COMMIT"
  echo "source_diff_sha256=$SOURCE_DIFF_SHA256"
  echo "hostname=$(hostname)"
  date '+datetime=%Y-%m-%d %H:%M:%S %Z'
  icc --version | sed -n '1p'
  "$MPIEXEC" --version 2>&1 | sed -n '1p'
  cmake --version | sed -n '1p'
  "$PYTHON" --version
  module -t list 2>&1
  echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
  echo "MKL_NUM_THREADS=$MKL_NUM_THREADS"
  echo "OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
  echo "BLIS_NUM_THREADS=$BLIS_NUM_THREADS"
  echo "PJM_MPI_PROC=4"
  echo "strict_fp=icc-fp-model-precise-no-fast-math"
  cmake -LA -N "$BUILD_DIR"
  ldd "$DRIVER"
}} > "$EVIDENCE_DIR/workflow/build-environment.txt"

run_point() {{
  stratum=$1
  count=$2
  replicate=$3
  kind=$4
  world_size=$5
  mpi_arguments=()
  if [[ $world_size -gt 1 ]]; then
    mpi_arguments=(--mpiexec "$MPIEXEC")
  fi
  "$PYTHON" "$CAMPAIGN" run-one \
    --driver "$DRIVER" "${{mpi_arguments[@]}}" \
    --mpi-world-size "$world_size" \
    --input-root "$INPUT_DIR" --identity "$IDENTITY" \
    --raw-dir "$EVIDENCE_DIR/raw" --stratum "$stratum" \
    --operator-count "$count" --saved-sources 64 \
    --replicate "$replicate" --kind "$kind" --timeout 3600
}}

CURRENT_PHASE=measurement_loop
for stratum in SC-ONEBODY SC-QUARTIC SC-MIXED SC-LOW-REUSE; do
  for count in 1 32 128 512; do
    run_point "$stratum" "$count" 0 warmup 1
    for replicate in 0 1 2 3 4 5 6; do
      run_point "$stratum" "$count" "$replicate" measured 1
    done
  done
  for world_size in 1 2 4; do
    run_point "$stratum" 512 0 parity "$world_size"
  done
  touch "$EVIDENCE_DIR/workflow/batch-${{stratum}}.complete"
done

EVIDENCE_JSON=$EVIDENCE_DIR/power-lanczos-p6c2-observable-scaling-evidence-v2.json
CURRENT_PHASE=evidence_assembly
"$PYTHON" "$CAMPAIGN" assemble \
  --raw-dir "$EVIDENCE_DIR/raw" --output "$EVIDENCE_JSON"
CURRENT_PHASE=evidence_validation
"$PYTHON" "$CAMPAIGN" validate-evidence --evidence "$EVIDENCE_JSON"

CURRENT_PHASE=inventory
{{
  echo "source_entries=$(find "$SOURCE_ROOT" -xdev | wc -l)"
  echo "source_bytes=$(du -sk "$SOURCE_ROOT" | awk '{{print $1 * 1024}}')"
  echo "build_entries=$(find "$BUILD_DIR" -xdev | wc -l)"
  echo "build_bytes=$(du -sk "$BUILD_DIR" | awk '{{print $1 * 1024}}')"
  echo "python_env_entries=$(find "$PYTHON_ENV" -xdev | wc -l)"
  echo "python_env_bytes=$(du -sk "$PYTHON_ENV" | awk '{{print $1 * 1024}}')"
  echo "raw_entries=$(find "$EVIDENCE_DIR/raw" -xdev | wc -l)"
  echo "raw_bytes=$(du -sk "$EVIDENCE_DIR/raw" | awk '{{print $1 * 1024}}')"
}} > "$EVIDENCE_DIR/workflow/remote-scratch-inventory.txt"

CAMPAIGN_PASS=$("$PYTHON" -c \
  'import json,sys; print("1" if json.load(open(sys.argv[1]))["passed"] else "0")' \
  "$EVIDENCE_JSON")
if [[ $CAMPAIGN_PASS != 1 ]]; then
  CURRENT_PHASE=campaign_result
  finalize_failure 2 campaign_evidence_failed
  exit 2
fi
CURRENT_PHASE=evidence_archival
archive_evidence
CURRENT_PHASE=completion_marker
touch "$RUN_DIR/COMPLETED"
FINALIZED=1
echo "P6-C2 scaling campaign completed with PASS evidence"
"""


def member(path: Path, root: Path, role: str):
    return {
        "path": path.relative_to(root).as_posix(),
        "role": role,
        "bytes": path.stat().st_size,
        "sha256": sha_file(path),
    }


def create_package(args) -> None:
    repo = args.repo.resolve()
    output = args.output_dir.resolve()
    python_wheelhouse_source = args.python_wheelhouse.resolve()
    require((repo / ".git").exists(), "repository metadata is required")
    require(not output.exists(), f"package output already exists: {output}")
    python_wheel_sources = validate_python_wheelhouse(python_wheelhouse_source)
    require(
        git(repo, "status", "--porcelain=v1", "--untracked-files=all") == b"",
        "source worktree must be clean",
    )
    source_commit = git(repo, "rev-parse", "HEAD").decode().strip()
    require(len(source_commit) == 40, "source commit")
    submodules = submodule_commits(repo, source_commit)
    require(submodules, "source archive requires initialized submodules")
    for relative, expected_commit in sorted(submodules.items()):
        submodule = repo / relative
        require(submodule.is_dir(), f"missing submodule: {relative}")
        actual_commit = git(submodule, "rev-parse", "HEAD").decode().strip()
        require(actual_commit == expected_commit, f"submodule commit: {relative}")
        require(
            git(submodule, "status", "--porcelain=v1", "--untracked-files=all")
            == b"",
            f"submodule worktree must be clean: {relative}",
        )

    output.mkdir(parents=True)
    workflow = output / "workflow"
    workflow.mkdir()
    short = source_commit[:12]
    source_name = f"mVMC-p6c2-scaling-source-{short}.tar.gz"
    source_path = output / source_name
    source_entries = git_tree_entries(repo, source_commit, "mVMC")
    for relative, commit in sorted(submodules.items()):
        source_entries.extend(
            git_tree_entries(repo / relative, commit, f"mVMC/{relative}")
        )
    write_tar_gz(source_path, source_entries)
    source_checksum = output / "source_archive.sha256"
    source_checksum.write_text(
        f"{sha_file(source_path)}  {source_name}\n", encoding="ascii"
    )

    python_wheelhouse = output / "python-wheelhouse"
    python_wheelhouse.mkdir()
    python_wheel_paths = []
    for source in python_wheel_sources:
        destination = python_wheelhouse / source.name
        shutil.copyfile(source, destination)
        os.chmod(destination, 0o644)
        python_wheel_paths.append(destination)
    python_requirements = output / "python-requirements.lock"
    python_requirements.write_text(python_requirements_lock(), encoding="ascii")
    python_wheelhouse_checksum = output / "python-wheelhouse.sha256"
    python_wheelhouse_checksum.write_text(
        "".join(
            f"{sha_file(path)}  python-wheelhouse/{path.name}\n"
            for path in python_wheel_paths
        ),
        encoding="ascii",
    )

    readme = workflow / "README.md"
    readme.write_text(
        "# P6-C2 observable scaling package\n\n"
        "This package is result-unobserved. Transfer and `pjsub` require "
        "explicit user confirmation. The job uses only `./tmp` below the "
        "submission directory and removes its extracted source/build tree "
        "and its Python virtual environment after checksummed evidence "
        "archival. Python dependencies are installed offline from the "
        "checksummed wheelhouse in this package.\n",
        encoding="utf-8",
    )
    job = workflow / "p6c2_scaling_genkai_job.sh"
    provisional = job_script(
        source_name,
        source_commit,
        SOURCE_DIFF_SHA256,
        args.project_code,
        args.resource_group,
    )
    job.write_text(provisional, encoding="utf-8")
    os.chmod(job, 0o755)

    manifest_path = output / "package_manifest.json"
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "package_id": args.package_id,
        "status": "prepared_result_unobserved",
        "result_observed": False,
        "production_authorized": False,
        "source": {
            "commit": source_commit,
            "source_diff_sha256": SOURCE_DIFF_SHA256,
            "submodule_commits": submodules,
            "archive": member(source_path, output, "source_archive"),
        },
        "execution": {
            "target": "genkai",
            "compiler": "intel/2023.2",
            "mpi": "mvapich/3.0-intel2023.2",
            "omp_num_threads": 1,
            "mpi_process_capacity": 4,
            "exact_grid": EXPECTED_GRID,
            "python_runtime": {
                "module": PYTHON_MODULE,
                "executable": PYTHON_EXECUTABLE,
                "version": PYTHON_VERSION,
                "wheel_target": PYTHON_WHEEL_TARGET,
                "environment": "submission-directory-relative ./tmp virtual environment",
                "installation": "pip --isolated --no-cache-dir --no-index --no-deps --require-hashes --only-binary",
                "requirements_lock": "python-requirements.lock",
                "wheelhouse_checksum": "python-wheelhouse.sha256",
                "wheels": [
                    {
                        "requirement": wheel["requirement"],
                        "path": f"python-wheelhouse/{wheel['filename']}",
                        "sha256": wheel["sha256"],
                    }
                    for wheel in PYTHON_WHEELS
                ],
            },
            "logical_batches": [
                "SC-ONEBODY",
                "SC-QUARTIC",
                "SC-MIXED",
                "SC-LOW-REUSE",
            ],
            "remote_scratch": "submission-directory-relative ./tmp only",
            "cleanup_condition": "after success or failure evidence archive and checksum are written by one-shot EXIT/signal finalization",
            "job_failure_artifacts": {
                "phase_tracking": True,
                "generic_nonzero_exit": True,
                "signal_exit": True,
                "atomic_archive_publication": True,
                "standalone_failure_metadata": True,
                "checksummed_partial_logs": True,
                "python_dependency_log": True,
            },
            "failure_records": [
                "timeout",
                "nonzero_exit",
                "missing_metric",
                "identity_mismatch",
                "scheduler_eviction",
            ],
        },
        "members": [
            member(source_checksum, output, "source_checksum"),
            member(readme, output, "workflow_readme"),
            member(job, output, "scheduler_script"),
            member(python_requirements, output, "python_requirements_lock"),
            member(
                python_wheelhouse_checksum,
                output,
                "python_wheelhouse_checksum",
            ),
            *[
                member(path, output, "python_dependency_wheel")
                for path in python_wheel_paths
            ],
        ],
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    final_manifest_sha = sha_file(manifest_path)
    manifest_checksum = output / "package_manifest.sha256"
    manifest_checksum.write_text(
        f"{final_manifest_sha}  package_manifest.json\n", encoding="ascii"
    )

    package_name = f"{PACKAGE_PREFIX}-package-{short}.tar.gz"
    package_path = output / package_name
    package_members = [
        (source_name, 0o644),
        ("source_archive.sha256", 0o644),
        ("package_manifest.json", 0o644),
        ("package_manifest.sha256", 0o644),
        ("workflow/README.md", 0o644),
        ("workflow/p6c2_scaling_genkai_job.sh", 0o755),
        ("python-requirements.lock", 0o644),
        ("python-wheelhouse.sha256", 0o644),
        *[
            (f"python-wheelhouse/{path.name}", 0o644)
            for path in python_wheel_paths
        ],
    ]
    write_tar_gz(package_path, archive_file_entries(output, package_members))
    (output / f"{package_name}.sha256").write_text(
        f"{sha_file(package_path)}  {package_name}\n", encoding="ascii"
    )
    validate_package(output, False)
    print(
        f"PASS P6-C2 scaling package commit={source_commit} "
        f"submodules={len(submodules)} package={package_path.name}"
    )


def validate_tar(path: Path, required_prefix: str | None = None) -> list[str]:
    names = []
    with tarfile.open(path, "r:gz") as archive:
        for item in archive.getmembers():
            require(safe_name(item.name), f"unsafe tar member: {item.name!r}")
            if required_prefix is not None:
                require(
                    item.name.startswith(required_prefix),
                    f"tar member outside prefix: {item.name}",
                )
            require(not item.isdev() and not item.isfifo(), "special tar member")
            if item.issym() or item.islnk():
                resolved_link = posixpath.normpath(
                    posixpath.join(posixpath.dirname(item.name), item.linkname)
                )
                require(
                    not item.linkname.startswith("/")
                    and resolved_link.startswith("mVMC/"),
                    f"unsafe tar link: {item.name}",
                )
            names.append(item.name)
    require(len(names) == len(set(names)), "duplicate tar member")
    return names


def validate_package(package_dir: Path, skip_package_archive: bool) -> None:
    root = package_dir.resolve()
    manifest_path = root / "package_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    require(manifest["schema_version"] == SCHEMA_VERSION, "package schema")
    require(
        manifest["status"] == "prepared_result_unobserved"
        and manifest["result_observed"] is False
        and manifest["production_authorized"] is False,
        "result-unobserved package boundary",
    )
    require(manifest["execution"]["exact_grid"] == EXPECTED_GRID, "exact grid")
    require(
        manifest["execution"]["python_runtime"]
        == {
            "module": PYTHON_MODULE,
            "executable": PYTHON_EXECUTABLE,
            "version": PYTHON_VERSION,
            "wheel_target": PYTHON_WHEEL_TARGET,
            "environment": "submission-directory-relative ./tmp virtual environment",
            "installation": "pip --isolated --no-cache-dir --no-index --no-deps --require-hashes --only-binary",
            "requirements_lock": "python-requirements.lock",
            "wheelhouse_checksum": "python-wheelhouse.sha256",
            "wheels": [
                {
                    "requirement": wheel["requirement"],
                    "path": f"python-wheelhouse/{wheel['filename']}",
                    "sha256": wheel["sha256"],
                }
                for wheel in PYTHON_WHEELS
            ],
        },
        "Python runtime",
    )
    require(
        manifest["execution"]["job_failure_artifacts"]
        == {
            "phase_tracking": True,
            "generic_nonzero_exit": True,
            "signal_exit": True,
            "atomic_archive_publication": True,
            "standalone_failure_metadata": True,
            "checksummed_partial_logs": True,
            "python_dependency_log": True,
        },
        "job failure artifacts",
    )
    members = [manifest["source"]["archive"], *manifest["members"]]
    for item in members:
        path = root / item["path"]
        require(
            path.is_file()
            and path.stat().st_size == item["bytes"]
            and sha_file(path) == item["sha256"],
            f"package member: {item['path']}",
        )
    validate_python_wheelhouse(root / "python-wheelhouse")
    require(
        (root / "python-requirements.lock").read_text(encoding="ascii")
        == python_requirements_lock(),
        "Python requirements lock",
    )
    expected_wheelhouse_checksum = "".join(
        f"{wheel['sha256']}  python-wheelhouse/{wheel['filename']}\n"
        for wheel in sorted(PYTHON_WHEELS, key=lambda item: item["filename"])
    )
    require(
        (root / "python-wheelhouse.sha256").read_text(encoding="ascii")
        == expected_wheelhouse_checksum,
        "Python wheelhouse checksum",
    )
    source_names = validate_tar(
        root / manifest["source"]["archive"]["path"], "mVMC/"
    )
    for required in [
        "mVMC/test/power_lanczos_observable_scaling_driver.c",
        "mVMC/test/python/power_lanczos_observable_scaling_campaign.py",
        "mVMC/test/python/package_power_lanczos_observable_scaling.py",
    ]:
        require(required in source_names, f"source archive member: {required}")
    for submodule in manifest["source"]["submodule_commits"]:
        require(
            any(name.startswith(f"mVMC/{submodule}/") for name in source_names),
            f"archived submodule: {submodule}",
        )
    job = (root / "workflow/p6c2_scaling_genkai_job.sh").read_text(
        encoding="utf-8"
    )
    for token in [
        "trap on_exit EXIT",
        "trap 'seal_on_signal HUP' HUP",
        "finalize_failure",
        "write_failure_metadata",
        "FINALIZATION_COUNT",
        "CURRENT_PHASE=configure",
        "CURRENT_PHASE=python_dependency_install",
        "CURRENT_PHASE=workflow_provenance",
        "python-dependency.log",
        "python_env_entries=",
        "module load intel/2023.2 mvapich/3.0-intel2023.2 python/3.12.11",
        "--no-cache-dir --no-index --no-deps --no-compile",
        "--only-binary=:all: --require-hashes",
        "MVMC_P6C2_JOB_SELFTEST_MODE",
        'mktemp -d "$RUN_DIR/tmp/',
        "seal-missing",
        "scheduler_eviction",
        "SC-LOW-REUSE",
        "for world_size in 1 2 4",
        "archive_evidence",
    ]:
        require(token in job, f"scheduler lifecycle token: {token}")
    provenance_order = [
        "CURRENT_PHASE=python_wheelhouse_checksum",
        "CURRENT_PHASE=workflow_provenance",
        "CURRENT_PHASE=source_extraction",
        "CURRENT_PHASE=python_dependency_install",
    ]
    require(
        [job.index(token) for token in provenance_order]
        == sorted(job.index(token) for token in provenance_order),
        "scheduler dependency provenance order",
    )
    require("$TMPDIR" not in job and "/var/tmp" not in job, "system temporary path")
    subprocess.run(["bash", "-n", str(root / "workflow/p6c2_scaling_genkai_job.sh")], check=True)
    if not skip_package_archive:
        packages = list(root.glob(f"{PACKAGE_PREFIX}-package-*.tar.gz"))
        require(len(packages) == 1, "single transfer package")
        names = validate_tar(packages[0])
        expected_names = sorted(
            [
                manifest["source"]["archive"]["path"],
                "package_manifest.json",
                "package_manifest.sha256",
                "source_archive.sha256",
                "workflow/README.md",
                "workflow/p6c2_scaling_genkai_job.sh",
                "python-requirements.lock",
                "python-wheelhouse.sha256",
                *[
                    f"python-wheelhouse/{wheel['filename']}"
                    for wheel in PYTHON_WHEELS
                ],
            ]
        )
        require(
            names == expected_names,
            "transfer package members",
        )
        sidecar = root / f"{packages[0].name}.sha256"
        expected = f"{sha_file(packages[0])}  {packages[0].name}\n"
        require(sidecar.read_text(encoding="ascii") == expected, "package checksum")


def validate_checksum_sidecar(path: Path, sidecar: Path) -> None:
    fields = sidecar.read_text(encoding="ascii").strip().split()
    require(
        len(fields) == 2
        and fields[0] == sha_file(path)
        and fields[1] == path.name,
        f"checksum sidecar: {sidecar.name}",
    )


def validate_evidence_archive(path: Path, required_suffixes: set[str]) -> None:
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        require(
            all(item.isfile() or item.isdir() for item in members),
            "self-test archive special member",
        )
        payloads = {
            item.name: archive.extractfile(item).read()
            for item in members
            if item.isfile()
        }
    ledger_names = [name for name in payloads if name.endswith("artifact-ledger.sha256")]
    require(
        len(ledger_names) == 1,
        f"self-test artifact ledger: ledgers={ledger_names} members={sorted(payloads)}",
    )
    ledger_name = ledger_names[0]
    observed = {}
    for raw_line in payloads[ledger_name].decode("ascii").splitlines():
        digest, relative = raw_line.split(maxsplit=1)
        observed[relative] = digest
    expected = {
        name: sha_bytes(payload)
        for name, payload in payloads.items()
        if name != ledger_name
    }
    require(observed == expected, "self-test artifact ledger contents")
    for suffix in required_suffixes:
        require(
            any(name.endswith(suffix) for name in payloads),
            f"self-test archive member: {suffix}",
        )


def self_test_lifecycle() -> None:
    source_commit = "1" * 40
    short = source_commit[:12]
    script_text = job_script(
        f"mVMC-p6c2-scaling-source-{short}.tar.gz",
        source_commit,
        SOURCE_DIFF_SHA256,
        "selftest-project",
        "selftest-resource",
    )
    local_tmp = Path.cwd() / "tmp"
    local_tmp.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="p6c2-package-lifecycle.", dir=local_tmp))
    try:
        cases = {
            "success": 0,
            "nonzero": 23,
            "dependency": 24,
            "signal": 130,
        }
        for mode, expected_exit in cases.items():
            case_root = work / mode
            case_root.mkdir()
            script = case_root / "p6c2_scaling_genkai_job.sh"
            script.write_text(script_text, encoding="utf-8")
            script.chmod(0o755)
            completed = subprocess.run(
                ["bash", str(script)],
                cwd=case_root,
                env={**os.environ, "MVMC_P6C2_JOB_SELFTEST_MODE": mode},
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            require(
                completed.returncode == expected_exit,
                f"self-test {mode} exit: {completed.returncode}",
            )
            require(not (case_root / "tmp").exists(), f"self-test {mode} cleanup")
            require(
                not list(case_root.glob("*.partial.*")),
                f"self-test {mode} partial artifact",
            )
            archive = (
                case_root
                / f"p6c2-observable-scaling-evidence-{short}.tar.gz"
            )
            checksum = (
                case_root / f"p6c2-observable-scaling-evidence-{short}.sha256"
            )
            require(archive.is_file() and checksum.is_file(), f"self-test {mode} archive")
            validate_checksum_sidecar(archive, checksum)
            if mode == "success":
                require((case_root / "COMPLETED").is_file(), "self-test success marker")
                require(not (case_root / "FAILED").exists(), "self-test success failure marker")
                validate_evidence_archive(
                    archive, {"workflow/selftest.log", "artifact-ledger.sha256"}
                )
                continue

            require(not (case_root / "COMPLETED").exists(), f"self-test {mode} completed")
            require((case_root / "FAILED").is_file(), f"self-test {mode} marker")
            metadata = (
                case_root / f"p6c2-observable-scaling-failure-{short}.json"
            )
            metadata_checksum = (
                case_root / f"p6c2-observable-scaling-failure-{short}.sha256"
            )
            require(
                metadata.is_file() and metadata_checksum.is_file(),
                f"self-test {mode} metadata",
            )
            validate_checksum_sidecar(metadata, metadata_checksum)
            record = json.loads(metadata.read_text(encoding="utf-8"))
            expected_phase = (
                "selftest_python_dependency"
                if mode == "dependency"
                else f"selftest_{mode}"
            )
            expected_reason = (
                "scheduler_signal_TERM" if mode == "signal" else "nonzero_exit"
            )
            require(
                record["status"] == "failed"
                and record["exit_code"] == expected_exit
                and record["phase"] == expected_phase
                and record["reason"] == expected_reason
                and record["finalization_count"] == 1,
                f"self-test {mode} failure metadata",
            )
            marker = dict(
                line.split("=", 1)
                for line in (case_root / "FAILED").read_text(encoding="ascii").splitlines()
            )
            require(
                marker["exit_code"] == str(expected_exit)
                and marker["phase"] == expected_phase
                and marker["reason"] == expected_reason
                and marker["finalization_count"] == "1"
                and marker["failure_metadata_status"] == "0"
                and marker["evidence_archive_status"] == "0",
                f"self-test {mode} failure marker",
            )
            validate_evidence_archive(
                archive,
                {
                    (
                        "python-dependency.log"
                        if mode == "dependency"
                        else "configure.log"
                    ),
                    f"workflow/{metadata.name}",
                    f"workflow/{metadata_checksum.name}",
                    "artifact-ledger.sha256",
                },
            )
        print(
            "PASS P6-C2 package lifecycle normal=1 nonzero=1 dependency=1 "
            "signal=1 finalize_once=3 cleanup=4 archives=4"
        )
    finally:
        resolved = work.resolve()
        if resolved.parent == local_tmp.resolve() and not resolved.is_symlink():
            shutil.rmtree(resolved)
        if local_tmp.exists() and not any(local_tmp.iterdir()):
            local_tmp.rmdir()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    create = subparsers.add_parser("create")
    create.add_argument("--repo", type=Path, required=True)
    create.add_argument("--output-dir", type=Path, required=True)
    create.add_argument("--python-wheelhouse", type=Path, required=True)
    create.add_argument("--package-id", required=True)
    create.add_argument("--project-code", default="pj25000164")
    create.add_argument("--resource-group", default="a-pj25000164")
    validate = subparsers.add_parser("validate")
    validate.add_argument("--package-dir", type=Path, required=True)
    validate.add_argument("--skip-package-archive", action="store_true")
    subparsers.add_parser("self-test-lifecycle")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.command == "create":
            create_package(args)
        elif args.command == "validate":
            validate_package(args.package_dir, args.skip_package_archive)
            print(f"PASS P6-C2 scaling package validation: {args.package_dir}")
        else:
            self_test_lifecycle()
        return 0
    except (
        PackageError,
        OSError,
        subprocess.CalledProcessError,
        json.JSONDecodeError,
        tarfile.TarError,
    ) as exc:
        print(f"FAIL P6-C2 scaling package: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
