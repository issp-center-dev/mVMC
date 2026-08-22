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
import subprocess
import sys
import tarfile
from pathlib import Path, PurePosixPath


SCHEMA_VERSION = 1
PACKAGE_PREFIX = "p6c2-observable-scaling"
SOURCE_DIFF_SHA256 = hashlib.sha256(b"tracked-diff\0").hexdigest()
EXPECTED_GRID = {"warmup": 16, "measured": 112, "parity": 12}


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
SCRATCH=
EVIDENCE_DIR=
INPUT_DIR=
IDENTITY=
CAMPAIGN=
OUTPUT_ARCHIVE=$RUN_DIR/p6c2-observable-scaling-evidence-{short}.tar.gz
OUTPUT_CHECKSUM=$RUN_DIR/p6c2-observable-scaling-evidence-{short}.sha256
SOURCE_COMMIT={source_commit}
SOURCE_DIFF_SHA256={source_diff_sha256}
export MVMC_P6C2_SOURCE_COMMIT=$SOURCE_COMMIT
export MVMC_P6C2_SOURCE_DIFF_SHA256=$SOURCE_DIFF_SHA256

cleanup() {{
  if [[ -n ${{SCRATCH:-}} && -d $SCRATCH && ! -L $SCRATCH ]]; then
    resolved=$(readlink -f "$SCRATCH")
    case "$resolved" in
      "$RUN_DIR"/tmp/*) rm -rf -- "$resolved" ;;
      *) echo "refusing unexpected scratch cleanup: $resolved" >&2 ;;
    esac
  fi
  rmdir "$RUN_DIR/tmp" 2>/dev/null || true
}}

archive_evidence() {{
  [[ -n ${{EVIDENCE_DIR:-}} && -d $EVIDENCE_DIR ]] || return 0
  find "$EVIDENCE_DIR" -type l -print -quit | grep -q . && return 1
  (
    cd "$EVIDENCE_DIR"
    find . -type f ! -name artifact-ledger.sha256 -print0 | sort -z | \
      xargs -0 sha256sum > artifact-ledger.sha256
    sha256sum -c artifact-ledger.sha256
  )
  tar -czf "$OUTPUT_ARCHIVE" -C "$EVIDENCE_DIR" .
  tar -tzf "$OUTPUT_ARCHIVE" >/dev/null
  sha256sum "$OUTPUT_ARCHIVE" > "$OUTPUT_CHECKSUM"
}}

seal_on_signal() {{
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
  archive_evidence
  exit 130
}}

trap cleanup EXIT
trap seal_on_signal HUP INT TERM

test -f "$SOURCE_ARCHIVE" && test ! -L "$SOURCE_ARCHIVE"
(cd "$RUN_DIR" && sha256sum -c package_manifest.sha256)
(cd "$RUN_DIR" && sha256sum -c "$(basename "$SOURCE_CHECKSUM")")

mkdir -p "$RUN_DIR/tmp"
SCRATCH=$(mktemp -d "$RUN_DIR/tmp/p6c2-scaling.XXXXXX")
SOURCE_DIR=$SCRATCH/source
BUILD_DIR=$SCRATCH/build
EVIDENCE_DIR=$SCRATCH/evidence
INPUT_DIR=$EVIDENCE_DIR/inputs
IDENTITY=$EVIDENCE_DIR/identity.json
mkdir "$SOURCE_DIR" "$BUILD_DIR" "$EVIDENCE_DIR"
tar -xzf "$SOURCE_ARCHIVE" -C "$SOURCE_DIR"
SOURCE_ROOT=$SOURCE_DIR/mVMC

module load intel/2023.2 mvapich/3.0-intel2023.2
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1

PYTHON=$(command -v python3.11 2>/dev/null || command -v python3)
MPIEXEC=$(command -v mpiexec)
CAMPAIGN=$SOURCE_ROOT/test/python/power_lanczos_observable_scaling_campaign.py
PACKAGER=$SOURCE_ROOT/test/python/package_power_lanczos_observable_scaling.py
"$PYTHON" "$PACKAGER" validate --package-dir "$RUN_DIR" \
  --skip-package-archive

cmake -S "$SOURCE_ROOT" -B "$BUILD_DIR" \
  -DCONFIG=intel \
  -DTesting=ON \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DCMAKE_C_FLAGS='-fp-model precise' \
  -DCMAKE_CXX_FLAGS='-fp-model precise' \
  -DCMAKE_Fortran_FLAGS='-fp-model precise' \
  -DPYTHON_EXECUTABLE="$PYTHON" \
  > "$EVIDENCE_DIR/configure.log" 2>&1
cmake --build "$BUILD_DIR" --parallel 1 \
  --target power_lanczos_observable_scaling_driver \
  > "$EVIDENCE_DIR/build.log" 2>&1

env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  BLIS_NUM_THREADS=1 ctest --test-dir "$BUILD_DIR" --output-on-failure \
  -R '^PowerLanczosObservableScaling(Campaign_Unit|Driver_Focused)$' \
  > "$EVIDENCE_DIR/focused-ctest.log" 2>&1

DRIVER=$BUILD_DIR/test/power_lanczos_observable_scaling_driver
"$PYTHON" "$CAMPAIGN" generate-inputs --output-dir "$INPUT_DIR"
"$PYTHON" "$CAMPAIGN" snapshot \
  --driver "$DRIVER" --build-dir "$BUILD_DIR" --output "$IDENTITY" \
  --source-commit "$SOURCE_COMMIT" \
  --source-diff-sha256 "$SOURCE_DIFF_SHA256"

mkdir "$EVIDENCE_DIR/raw" "$EVIDENCE_DIR/workflow"
cp "$PACKAGE_MANIFEST" "$EVIDENCE_DIR/workflow/"
cp "$RUN_DIR/package_manifest.sha256" "$EVIDENCE_DIR/workflow/"
cp "$SOURCE_CHECKSUM" "$EVIDENCE_DIR/workflow/"
cp "$RUN_DIR/workflow/p6c2_scaling_genkai_job.sh" \
  "$EVIDENCE_DIR/workflow/"

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
"$PYTHON" "$CAMPAIGN" assemble \
  --raw-dir "$EVIDENCE_DIR/raw" --output "$EVIDENCE_JSON"
"$PYTHON" "$CAMPAIGN" validate-evidence --evidence "$EVIDENCE_JSON"

{{
  echo "source_entries=$(find "$SOURCE_ROOT" -xdev | wc -l)"
  echo "source_bytes=$(du -sk "$SOURCE_ROOT" | awk '{{print $1 * 1024}}')"
  echo "build_entries=$(find "$BUILD_DIR" -xdev | wc -l)"
  echo "build_bytes=$(du -sk "$BUILD_DIR" | awk '{{print $1 * 1024}}')"
  echo "raw_entries=$(find "$EVIDENCE_DIR/raw" -xdev | wc -l)"
  echo "raw_bytes=$(du -sk "$EVIDENCE_DIR/raw" | awk '{{print $1 * 1024}}')"
}} > "$EVIDENCE_DIR/workflow/remote-scratch-inventory.txt"

CAMPAIGN_PASS=$("$PYTHON" -c \
  'import json,sys; print("1" if json.load(open(sys.argv[1]))["passed"] else "0")' \
  "$EVIDENCE_JSON")
archive_evidence
touch "$RUN_DIR/COMPLETED"
if [[ $CAMPAIGN_PASS != 1 ]]; then
  echo "P6-C2 scaling campaign completed with FAIL evidence" >&2
  exit 2
fi
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
    require((repo / ".git").exists(), "repository metadata is required")
    require(not output.exists(), f"package output already exists: {output}")
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

    readme = workflow / "README.md"
    readme.write_text(
        "# P6-C2 observable scaling package\n\n"
        "This package is result-unobserved. Transfer and `pjsub` require "
        "explicit user confirmation. The job uses only `./tmp` below the "
        "submission directory and removes its extracted source/build tree "
        "after checksummed evidence archival.\n",
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
            "logical_batches": [
                "SC-ONEBODY",
                "SC-QUARTIC",
                "SC-MIXED",
                "SC-LOW-REUSE",
            ],
            "remote_scratch": "submission-directory-relative ./tmp only",
            "cleanup_condition": "after evidence archive and checksum are written; EXIT/signal trap otherwise",
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
    members = [manifest["source"]["archive"], *manifest["members"]]
    for item in members:
        path = root / item["path"]
        require(
            path.is_file()
            and path.stat().st_size == item["bytes"]
            and sha_file(path) == item["sha256"],
            f"package member: {item['path']}",
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
        "trap cleanup EXIT",
        "trap seal_on_signal HUP INT TERM",
        'mktemp -d "$RUN_DIR/tmp/',
        "seal-missing",
        "scheduler_eviction",
        "SC-LOW-REUSE",
        "for world_size in 1 2 4",
        "archive_evidence",
    ]:
        require(token in job, f"scheduler lifecycle token: {token}")
    require("$TMPDIR" not in job and "/var/tmp" not in job, "system temporary path")
    subprocess.run(["bash", "-n", str(root / "workflow/p6c2_scaling_genkai_job.sh")], check=True)
    if not skip_package_archive:
        packages = list(root.glob(f"{PACKAGE_PREFIX}-package-*.tar.gz"))
        require(len(packages) == 1, "single transfer package")
        names = validate_tar(packages[0])
        require(
            names
            == [
                manifest["source"]["archive"]["path"],
                "package_manifest.json",
                "package_manifest.sha256",
                "source_archive.sha256",
                "workflow/README.md",
                "workflow/p6c2_scaling_genkai_job.sh",
            ],
            "transfer package members",
        )
        sidecar = root / f"{packages[0].name}.sha256"
        expected = f"{sha_file(packages[0])}  {packages[0].name}\n"
        require(sidecar.read_text(encoding="ascii") == expected, "package checksum")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    create = subparsers.add_parser("create")
    create.add_argument("--repo", type=Path, required=True)
    create.add_argument("--output-dir", type=Path, required=True)
    create.add_argument("--package-id", required=True)
    create.add_argument("--project-code", default="pj25000164")
    create.add_argument("--resource-group", default="a-pj25000164")
    validate = subparsers.add_parser("validate")
    validate.add_argument("--package-dir", type=Path, required=True)
    validate.add_argument("--skip-package-archive", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.command == "create":
            create_package(args)
        else:
            validate_package(args.package_dir, args.skip_package_archive)
            print(f"PASS P6-C2 scaling package validation: {args.package_dir}")
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
