#!/usr/bin/env python3

import re
import sys


FUNCTIONS = (
    "CalculateHamiltonianBF_fcmp",
    "CalculateGreenFuncBF",
)


def function_segment(source, name):
    match = re.search(r"\b{}\s*\(".format(re.escape(name)), source)
    if match is None:
        raise AssertionError("function not found: {}".format(name))
    next_function = re.search(
        r"\n(?:double complex|double|int|void)\s+[A-Za-z_][A-Za-z0-9_]*\s*\(",
        source[match.end():],
    )
    end = (
        match.end() + next_function.start()
        if next_function is not None
        else len(source)
    )
    return source[match.start():end]


def parallel_clause(segment, name):
    marker = "#pragma omp parallel default(shared)"
    start = segment.find(marker)
    if start < 0:
        raise AssertionError("outer OpenMP parallel clause not found: {}".format(name))
    end = segment.find("\n  {", start)
    if end < 0:
        raise AssertionError("outer OpenMP parallel clause is truncated: {}".format(name))
    return segment[start:end]


def main():
    if len(sys.argv) != 2:
        print("usage: {} <backflow_measure.c>".format(sys.argv[0]))
        return 2

    with open(sys.argv[1]) as source_file:
        source = source_file.read()

    for name in FUNCTIONS:
        clause = parallel_clause(function_segment(source, name), name)
        private_match = re.search(r"private\s*\((.*?)\)", clause, re.DOTALL)
        if private_match is None:
            raise AssertionError("private clause not found: {}".format(name))
        private_names = set(
            token.strip() for token in private_match.group(1).split(",")
        )
        if "idx" not in private_names:
            raise AssertionError(
                "{} outer OpenMP private clause is missing idx".format(name)
            )

    print("BackFlow complex OpenMP private-clause contract PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
