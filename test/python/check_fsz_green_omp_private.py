#!/usr/bin/env python3

import re
import sys


FUNCTIONS = ("CalculateGreenFunc_fsz", "CalculateGreenFuncBF_fsz")
SCALAR_TYPE = re.compile(
    r"(?m)^\s*(?:const\s+)?"
    r"(?:unsigned\s+|signed\s+)?"
    r"(?:int|long(?:\s+long)?|size_t|double(?:\s+complex)?|"
    r"BFFSZ[A-Za-z0-9_]*|BFNBody[A-Za-z0-9_]*)\s+([^;]+);"
)
ASSIGNMENT = re.compile(
    r"(?<![.>A-Za-z0-9_])(?:"
    r"(?:\+\+|--)\s*([A-Za-z_][A-Za-z0-9_]*)|"
    r"([A-Za-z_][A-Za-z0-9_]*)\s*"
    r"(?:\+\+|--|[+\-*/%&|^]?=(?!=)))"
)


def strip_comments(source):
    def blank(match):
        return "".join("\n" if char == "\n" else " " for char in match.group(0))

    source = re.sub(r"/\*.*?\*/", blank, source, flags=re.DOTALL)
    return re.sub(r"//[^\n]*", blank, source)


def matching_delimiter(source, opening, left, right):
    depth = 0
    for index in range(opening, len(source)):
        if source[index] == left:
            depth += 1
        elif source[index] == right:
            depth -= 1
            if depth == 0:
                return index
    raise AssertionError("unmatched {} in source".format(left))


def function_segment(source, name):
    match = re.search(r"\b{}\s*\([^;]*?\)\s*\{{".format(re.escape(name)), source, re.DOTALL)
    if match is None:
        raise AssertionError("function not found: {}".format(name))
    opening = source.find("{", match.start())
    closing = matching_delimiter(source, opening, "{", "}")
    return source[match.start():closing + 1]


def scalar_locals(segment):
    names = set()
    for match in SCALAR_TYPE.finditer(segment):
        for declarator in match.group(1).split(","):
            declarator = declarator.strip()
            if not declarator or declarator.startswith("*") or "[" in declarator:
                continue
            name = re.match(r"([A-Za-z_][A-Za-z0-9_]*)", declarator)
            if name is not None:
                names.add(name.group(1))
    return names


def clause_names(clause, kind):
    names = set()
    pattern = r"\b{}\s*\((.*?)\)".format(kind)
    for match in re.finditer(pattern, clause, re.DOTALL):
        payload = match.group(1)
        if kind == "reduction":
            if ":" not in payload:
                continue
            payload = payload.split(":", 1)[1]
        names.update(
            token.strip() for token in payload.split(",") if token.strip()
        )
    return names


def outer_parallel_private(segment):
    match = re.search(r"#pragma\s+omp\s+parallel(?P<clause>.*?)\{", segment, re.DOTALL)
    if match is None:
        raise AssertionError("outer OpenMP parallel region not found")
    return clause_names(match.group("clause"), "private")


def omp_for_loops(segment):
    pattern = re.compile(
        r"#pragma[ \t]+omp[ \t]+for\b(?P<clause>[^\n]*)\n"
        r"[ \t]*for[ \t]*\(",
    )
    for match in pattern.finditer(segment):
        header_open = segment.find("(", match.start())
        header_close = matching_delimiter(segment, header_open, "(", ")")
        body_open = segment.find("{", header_close)
        if body_open < 0:
            raise AssertionError("OpenMP for loop has no braced body")
        body_close = matching_delimiter(segment, body_open, "{", "}")
        yield (
            match.group("clause"),
            segment[header_open + 1:header_close],
            segment[body_open + 1:body_close],
        )


def assigned_scalars(text, locals_):
    assigned = set()
    for match in ASSIGNMENT.finditer(text):
        name = match.group(1) or match.group(2)
        if name in locals_:
            assigned.add(name)
    return assigned


def validate_function(source, name):
    segment = function_segment(source, name)
    locals_ = scalar_locals(segment)
    parallel_private = outer_parallel_private(segment)
    loop_count = 0
    for clause, header, body in omp_for_loops(segment):
        loop_count += 1
        assigned = assigned_scalars(header + "\n" + body, locals_)
        body_locals = scalar_locals(body)
        allowed = (
            parallel_private
            | clause_names(clause, "private")
            | clause_names(clause, "firstprivate")
            | clause_names(clause, "lastprivate")
            | clause_names(clause, "reduction")
            | body_locals
        )
        missing = sorted(assigned - allowed)
        if missing:
            label = " ".join(header.split())
            raise AssertionError(
                "{} OpenMP loop '{}' has assigned scalar(s) missing from "
                "private/reduction ownership: {}".format(
                    name, label, ", ".join(missing)
                )
            )
    if loop_count == 0:
        raise AssertionError("no OpenMP for loops found: {}".format(name))
    return loop_count


def validate_source(source):
    return sum(validate_function(source, name) for name in FUNCTIONS)


def mutation_self_check(source):
    mutated, count = re.subn(
        r"(#pragma\s+omp\s+for\s+private\([^)]*)\btmp\b\s*,?",
        r"\1",
        source,
        count=1,
    )
    if count != 1:
        raise AssertionError("could not create missing-private mutation")
    try:
        validate_source(mutated)
    except AssertionError as error:
        if "missing from private/reduction ownership: tmp" not in str(error):
            raise
        return
    raise AssertionError("missing-private mutation was not rejected")


def main():
    if len(sys.argv) != 2:
        print("usage: {} <calgrn_fsz.c>".format(sys.argv[0]))
        return 2

    with open(sys.argv[1]) as source_file:
        source = strip_comments(source_file.read())

    loop_count = validate_source(source)
    mutation_self_check(source)
    print(
        "FSZ Green-function OpenMP loop-private contract PASS "
        "(loops={}, mutation rejected)".format(loop_count)
    )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AssertionError as error:
        print("FSZ Green-function OpenMP loop-private contract FAIL")
        print(error)
        sys.exit(1)
