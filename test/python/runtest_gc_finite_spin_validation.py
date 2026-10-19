"""Reject damaged finite-spin probe output and nonfinite oracle values.

Use real production output as the control, then corrupt only its text before
the normal parser/comparator reads it. No production source is mutated.
"""
import argparse
import contextlib
import io
import os
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import runtest_gc_finite_spin as runner


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary")
    parser.add_argument("--work", required=True)
    args = parser.parse_args()
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    controls = {}
    for name in ("L4_S1_general", "L4_S1_triplet", "L4_S1_near"):
        fixture, output = work / (name + ".txt"), work / (name + ".probe")
        runner.write_fixture(str(fixture), name, 4)
        runner.run([os.path.abspath(args.binary), "--probe", str(fixture),
                    str(output)])
        controls[name] = output.read_text().splitlines()

    checked = 0

    def check(label, lines, fixture="L4_S1_general", reject=None):
        nonlocal checked
        case = SimpleNamespace(binary=args.binary, fixture=fixture, gauss=4,
                               work=str(work / "comparison"))

        def write_probe(command):
            Path(command[-1]).write_text("\n".join(lines) + "\n")
            return ""

        try:
            with mock.patch.object(runner, "run", write_probe), \
                    contextlib.redirect_stdout(io.StringIO()):
                runner.probe_case(case)
        except AssertionError as exc:
            if reject is None or reject not in str(exc):
                raise AssertionError("{}: unexpected failure: {}".format(
                    label, exc)) from exc
        else:
            if reject is not None:
                raise AssertionError("{}: corrupted probe was accepted".format(label))
        checked += 1

    for fixture, lines in controls.items():
        check("control " + fixture, lines, fixture)

    # Every numeric component, including imaginary components and both
    # derivative directions, must reject NaN and either sign of infinity.
    numeric_columns = {"Q": range(2, 7), "A": range(4, 6),
                       "D": range(3, 7), "Z": range(3, 5),
                       "G1": range(4, 6), "PA": range(4, 6),
                       "PR": range(4, 6), "G2": range(6, 8)}
    for group, columns in numeric_columns.items():
        fixture = "L4_S1_triplet" if group == "Z" else "L4_S1_general"
        lines = controls[fixture]
        # For Green rows choose a nonzero value: zeros already had an exact
        # equality check, whereas the ratio comparison used to accept NaN.
        index = next(i for i, line in enumerate(lines)
                     if line.split()[0] == group and
                     (group not in ("G1", "PA", "PR", "G2") or
                      any(float.fromhex(v) != 0 for v in line.split()[-2:])))
        for column in columns:
            for token in ("nan", "inf", "-inf"):
                damaged = list(lines)
                words = damaged[index].split()
                words[column] = token
                damaged[index] = " ".join(words)
                check("{} column {} {}".format(group, column, token),
                      damaged, fixture, reject="nonfinite")

    # The ratio-conditioning cut must not hide damaged output. The near
    # fixture has tiny nonzero amplitudes that are excluded from ratios.
    lines = controls["L4_S1_near"]
    amplitudes = {int(w[1]): abs(complex(float.fromhex(w[4]), float.fromhex(w[5])))
                  for w in (line.split() for line in lines) if w[0] == "A"}
    small = {mask for mask, value in amplitudes.items()
             if 0 < value < runner.CONDITION_CUT * max(amplitudes.values())}
    assert small, "near fixture no longer exercises the ratio cut"
    index = next(i for i, line in enumerate(lines)
                 if line.startswith("D ") and int(line.split()[1]) in small)
    damaged = list(lines)
    words = damaged[index].split()
    words[3] = "nan"
    damaged[index] = " ".join(words)
    check("NaN below ratio cut", damaged, "L4_S1_near", reject="nonfinite")

    # One missing row, an entirely absent group, duplicate rows and extra
    # keys all fail. A same-length file with the wrong operator also fails.
    for group in ("D", "Z", "G1", "PA", "PR", "G2"):
        fixture = "L4_S1_triplet" if group == "Z" else "L4_S1_general"
        lines = controls[fixture]
        index = next(i for i, line in enumerate(lines) if line.startswith(group + " "))
        check("missing " + group, lines[:index] + lines[index + 1:], fixture,
              reject=group + " keys differ")
        check("absent " + group,
              [line for line in lines if not line.startswith(group + " ")],
              fixture, reject=group + " keys differ")
        check("duplicate " + group, lines + [lines[index]], fixture,
              reject="duplicate " + group)
        words = lines[index].split()
        words[2] = "999"  # invalid orbital/parameter with unchanged row count
        damaged = list(lines)
        damaged[index] = " ".join(words)
        check("wrong key " + group, damaged, fixture, reject=group + " keys differ")
        check("extra " + group, lines + [" ".join(words)], fixture,
              reject=group + " keys differ")
    lines = controls["L4_S1_general"]
    check("all Green groups absent",
          [line for line in lines if line.split()[0] not in ("G1", "PA", "PR", "G2")],
          reject="G1 keys differ")
    for group in ("Q", "A"):
        row = next(line for line in lines if line.startswith(group + " "))
        check("duplicate " + group, lines + [row], reject="duplicate " + group)

    # Oracle failures must fail as well, before comparisons/conditioning
    # could turn NaN into an apparently zero maximum error.
    nodes = runner.oracle.independent_spin_nodes
    derivative = runner.oracle.affine_spin_derivative
    wave = runner.oracle.masked_finite_spin
    for bad in (float("nan"), complex(0, float("inf")), -float("inf")):
        def bad_nodes(*a, **kw):
            result = nodes(*a, **kw)
            result[0]["weight"] = bad
            return result

        def bad_derivative(*a, **kw):
            result = derivative(*a, **kw)
            return {key: bad for key in result}

        def bad_zero_derivative(*a, **kw):
            result = derivative(*a, **kw)
            return {key: bad if runner.oracle.popcount(key) == 4 else value
                    for key, value in result.items()}

        def bad_wave(*a, **kw):
            result = wave(*a, **kw)
            if not kw.get("exact_arithmetic", False):
                result = {key: bad for key in result}
            return result

        for function, replacement, fixture in (
                ("independent_spin_nodes", bad_nodes, "L4_S1_general"),
                ("qfloat", lambda value: bad, "L4_S1_general"),
                ("affine_spin_derivative", bad_derivative, "L4_S1_general"),
                ("affine_spin_derivative", bad_zero_derivative, "L4_S1_triplet"),
                ("masked_finite_spin", bad_wave, "L4_S1_general")):
            with mock.patch.object(runner.oracle, function, replacement):
                check("nonfinite reference " + function, controls[fixture],
                      fixture, reject="nonfinite")

    print("finite-spin probe validation: {} checks PASS".format(checked))


if __name__ == "__main__":
    main()
