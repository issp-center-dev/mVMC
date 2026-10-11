"""Fault-inject child output to test the actual death and audit checkers."""
import argparse
import os
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

import runtest_gc_antiparallel as runner


REPORTS = ("AddressSanitizer: heap-buffer-overflow",
           "UndefinedBehaviorSanitizer: undefined-behavior",
           "ERROR: LeakSanitizer: detected memory leaks",
           "WARNING: ThreadSanitizer: data race",
           "runtime error: signed integer overflow",
           "Segmentation fault", "Signal: 11", "signal 11", "Abort trap",
           "process exited on signal 6 (Aborted)", "Bus error")


def executable(path, body):
    path.write_text("#!{}\n{}".format(sys.executable, body))
    path.chmod(0o755)
    return str(path)


def death_guard(work):
    binary = executable(work / "death_child.py", """
import os, sys
print('GC anti-parallel sampler nonfinite', file=sys.stderr)
print(os.environ['GUARD_REPORT'], file=getattr(sys, os.environ['GUARD_STREAM']))
sys.exit(1)
""")
    wrapper = str(Path(__file__).with_name("runtest_gc_antiparallel_death.py"))
    missed = []
    for stream in ("stdout", "stderr"):
        for report in ("",) + REPORTS:
            environment = dict(os.environ, GUARD_REPORT=report,
                               GUARD_STREAM=stream)
            result = subprocess.run(
                [sys.executable, wrapper, binary,
                 "--antiparallel-death-nonfinite"], env=environment,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                universal_newlines=True, timeout=40)
            if (result.returncode == 0) != (report == ""):
                missed.append((stream, report, result.stdout))
    if missed:
        raise AssertionError("death checker accepted fault output: {}".format(missed))


def audit_guard(rootdir, work):
    # Keep fixture writing, real vmc execution, and audit parsing intact. The
    # child adds a report only after the real binary's expected audit failure.
    real_binary = os.path.abspath(runner.binary_path(rootdir))
    binary = executable(work / "audit_child.py", """
import os, subprocess, sys
code = subprocess.call([%r] + sys.argv[1:])
if code > 0:
    print(os.environ['GUARD_REPORT'], file=getattr(sys, os.environ['GUARD_STREAM']))
sys.exit(code)
""" % real_binary)
    args = argparse.Namespace(boundary="pbc", np=1)
    missed = []
    for stream in ("stdout", "stderr"):
        for report in ("",) + REPORTS:
            with patch.dict(os.environ, GUARD_REPORT=report, GUARD_STREAM=stream), \
                    patch.object(runner, "binary_path", return_value=binary):
                try:
                    runner.audit_case(str(work), args)
                except AssertionError as error:
                    if not report or "signal/sanitizer" not in str(error):
                        raise
                else:
                    if report:
                        missed.append((stream, report))
    if missed:
        raise AssertionError("audit checker accepted fault output: {}".format(missed))


if __name__ == "__main__":
    case = sys.argv[1]
    rootdir = os.getcwd()
    work = Path(runner.prepare_work(rootdir, "failure_guard_" + case))
    if case == "death":
        death_guard(work)
    elif case == "audit":
        audit_guard(rootdir, work)
    else:
        raise SystemExit("unknown guard case " + case)
    print("GC anti-parallel {} failure guard passed".format(case))
