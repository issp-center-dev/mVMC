"""Exercise MultiDef initialization and compare groups with independent runs."""

import argparse
import json
import math
import os
from pathlib import Path
import shlex
import signal
import subprocess
import tempfile

from runtest_varbin import prepare_case


def run(binary, directory, arguments, ranks, launcher, log_name):
    command = [str(binary)] + arguments
    if launcher:
        command = (shlex.split(os.environ.get("MVMC_MPIEXEC", "mpiexec")) +
                   [os.environ.get("MVMC_MPIEXEC_NUMPROC_FLAG", "-n"),
                    str(ranks)] + command)
    environment = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        environment[name] = "1"
    process = subprocess.Popen(command, cwd=directory, env=environment,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, errors="replace", start_new_session=True)
    try:
        output, _ = process.communicate(timeout=45)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        output, _ = process.communicate()
        (directory / log_name).write_text(output)
        raise AssertionError("MultiDef did not terminate within 45 seconds")
    (directory / log_name).write_text(output)
    return process.returncode, output


def check_failure(result, diagnostic):
    code, output = result
    diagnostics = (diagnostic,) if isinstance(diagnostic, str) else diagnostic
    if code == 0 or not any(text in output for text in diagnostics):
        raise AssertionError("expected failure containing {!r}, got {}:\n{}".format(
            diagnostic, code, output))
    bad = ("Segmentation fault", "Floating point exception", "AddressSanitizer",
           "runtime error:", "SIGSEGV", "SIGFPE")
    if code < 0 or any(marker in output for marker in bad):
        raise AssertionError("input error caused a signal/sanitizer failure:\n" + output)


def fixture(path, seed=2468, split=1):
    path.parent.mkdir(parents=True, exist_ok=True)
    prepare_case(path, complex_type=0, gutz_count=1, mode=1)
    modpara = path / "modpara.def"
    modpara.write_text(modpara.read_text().replace(
        "RndSeed 2468", "RndSeed {}".format(seed)).replace(
        "NSplitSize 1", "NSplitSize {}".format(split)))


def check_same_output(reference, actual):
    for name in ("zvo_out_007.dat", "zvo_out_008.dat",
                 "zvo_var_007.dat", "zvo_var_008.dat"):
        expected = (reference / "output" / name).read_bytes()
        observed = (actual / "output" / name).read_bytes()
        if expected != observed:
            raise AssertionError("MultiDef changed " + name)
        if not expected or not all(math.isfinite(float(x)) for x in expected.split()):
            raise AssertionError("empty/non-finite output: " + name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--mpi", action="store_true")
    parser.add_argument("--only", choices=("invalid", "normal"))
    args = parser.parse_args()
    binary = args.binary.resolve()
    if args.ranks < 1 or (args.ranks != 1 and not args.mpi):
        parser.error("multiple ranks require --mpi")
    work = Path.cwd() / "work"
    work.mkdir(exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="MultiDef_{}_".format(args.ranks), dir=work))
    print("MultiDef evidence: {}".format(root), flush=True)
    results = []

    def case(name, function):
        directory = root / name
        directory.mkdir()
        try:
            function(directory)
            results.append({"case": name, "status": "PASS"})
            print(name + ": PASS", flush=True)
        except (AssertionError, OSError, ValueError) as error:
            results.append({"case": name, "status": "FAIL", "error": str(error)})
            print("{}: FAIL: {}".format(name, error), flush=True)

    def invalid(directory, number, content, diagnostic, list_name="dirs"):
        if content is not None:
            (directory / list_name).write_text(content)
        result = run(binary, directory,
                     ["-m", str(number), list_name, "namelist.def"],
                     args.ranks, args.mpi, "run.log")
        check_failure(result, diagnostic)

    if args.only != "normal":
        for number in (0, -1, -2147483648, args.ranks + 1, 2147483647):
            case("count_{}".format(number), lambda d, n=number: invalid(
                d, n, None, "N must be between 1 and MPI size"))
        # macOS strtol() reports EINVAL for "x"; glibc reaches the no-digit check.
        for number, message in (("2147483648", "Numerical result out of range"),
                                ("9999999999999999999999999999", "strtol()"),
                                ("x", ("No digits were found",
                                       "strtol(): Invalid argument"))):
            case("parse_" + number, lambda d, n=number, m=message: invalid(d, n, None, m))
        case("missing_list", lambda d: invalid(d, 1, None, "DirListFile"))
        case("empty_list", lambda d: invalid(d, 1, " \n\t", "incomplete"))
        case("missing_directory", lambda d: invalid(d, 1, "absent\n", "chdir()"))
        case("not_a_directory", lambda d: invalid(d, 1, "dirs\n", "chdir()"))
        for length in (256, 4096):
            case("long_name_{}".format(length), lambda d, n=length: invalid(
                d, 1, "a" * n + "\n", "directory name is too long"))
        # Directories can be opened as streams on Linux; fgetc reports EISDIR.
        case("unreadable_list", lambda d: invalid(
            d, 1, None, "incomplete or unreadable", list_name="."))
        if args.ranks > 1:
            case("short_list", lambda d: invalid(d, args.ranks, ".\n", "incomplete"))
            # A nonzero world rank owns the failing group. Other roots must not
            # proceed into ReadDefFile (the valid group has no namelist).
            def late_chdir(directory):
                invalid(directory, 2, ".\nabsent\n", "chdir()")
                text = (directory / "run.log").read_text()
                if "namelist.def" in text:
                    raise AssertionError("another group read input after chdir failure")
            case("other_group_chdir", late_chdir)

    def normal(directory, groups, split=1, long_name=False):
        group_sizes = [args.ranks // groups + (i < args.ranks % groups)
                       for i in range(groups)]
        names = []
        for index, size in enumerate(group_sizes):
            name = ("a" * 120 + "/" + "b" * 134 if long_name
                    else "group_{}".format(index))
            names.append(name)
            # Give each directory distinct input for its independent-run comparison.
            seed = 2468 + 17 * index
            fixture(directory / name, seed, split)
            reference = directory / "reference_{}".format(index)
            fixture(reference, seed, split)
            code, output = run(binary, reference,
                               ["-e", "namelist.def", "initial.def"],
                               size, args.mpi, "run.log")
            if code:
                raise AssertionError("independent run failed:\n" + output)
        # Exercise all fscanf-style whitespace, and EOF without a final newline.
        (directory / "dirs").write_text(" \n\t" + "\t \n".join(names))
        code, output = run(binary, directory,
                           ["-m", str(groups), "dirs", "namelist.def", "initial.def"],
                           args.ranks, args.mpi, "run.log")
        if code:
            raise AssertionError("MultiDef run failed:\n" + output)
        if args.ranks % groups and "load imbalance" not in output:
            raise AssertionError("missing load imbalance diagnostic")
        for index, name in enumerate(names):
            check_same_output(directory / "reference_{}".format(index), directory / name)

    if args.only != "invalid":
        for groups in sorted({1, args.ranks, max(1, args.ranks - 1), min(2, args.ranks)}):
            case("groups_{}".format(groups), lambda d, g=groups: normal(d, g))
        if args.ranks >= 2:
            case("inner_split_2", lambda d: normal(d, max(1, args.ranks // 2), split=2))
        case("name_boundary_255", lambda d: normal(d, 1, long_name=True))
        def overridden_mode(directory):
            fixture(directory / "expert")
            code, output = run(binary, directory / "expert",
                               ["-m", "0", "-e", "namelist.def", "initial.def"],
                               args.ranks, args.mpi, "run.log")
            if code:
                raise AssertionError("-e no longer overrides -m:\n" + output)
        case("expert_override", overridden_mode)

    (root / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    failures = sum(result["status"] != "PASS" for result in results)
    print("MultiDef: {}/{} passed".format(len(results) - failures, len(results)), flush=True)
    return int(failures != 0)


if __name__ == "__main__":
    raise SystemExit(main())
