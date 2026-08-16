#!/usr/bin/env python3

import pathlib
import subprocess
import sys
import tempfile


def run_case(markov, reparse, trace, arguments, required):
    generated = subprocess.run(
        [str(markov)] + arguments,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if generated.returncode != 0:
        print(generated.stderr, file=sys.stderr)
        return False
    trace.write_text(generated.stdout, encoding="utf-8")
    replayed = subprocess.run(
        [str(reparse), str(trace)],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if replayed.returncode != 0:
        print(replayed.stderr, file=sys.stderr)
        return False
    missing = [needle for needle in required if needle not in replayed.stdout]
    if missing:
        print("missing reparse markers: " + repr(missing), file=sys.stderr)
        return False
    return True


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: runtest_p4f_matrix_reparse.py MARKOV REPARSE",
              file=sys.stderr)
        return 2

    markov = pathlib.Path(sys.argv[1]).resolve()
    reparse = pathlib.Path(sys.argv[2]).resolve()
    with tempfile.TemporaryDirectory(prefix="mvmc-p4f-reparse-") as tmp:
        trace = pathlib.Path(tmp) / "trace.log"
        common_required = (
            "REPARSE schema=1",
            " sample_count=128 ",
            " block_count=16 block_length=8 ",
            " all_match=1",
            "GEVP scope=full dimension=2 cutoff=1e-10 status=ok valid=1",
            "GEVP scope=full dimension=3 cutoff=1e-10 status=ok valid=1",
            "GEVP scope=prefix128_leave15 dimension=3 cutoff=1e-10 "
            "status=ok valid=1",
            "REPARSE_DECISION decision=GO",
        )
        if not run_case(
            markov, reparse, trace,
            ["4", "1", "128", "0", "1e-4", "0x5034535F52305631"],
            common_required + (
                "source_schema=2 fixture=p4s_bounded_markov_real",
                "source_antihermitian_convention=corrected_full_frobenius",
            ),
        ):
            return 1
        if not run_case(
            markov, reparse, trace,
            ["4", "1", "128", "4194304", "1e-2",
             "0x5034535F52305631", "0", "1", "7", "session"],
            common_required + (
                "source_schema=4 "
                "fixture=p4s9_long_direct_session_official",
                "source_antihermitian_convention=corrected_full_frobenius",
            ),
        ):
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
