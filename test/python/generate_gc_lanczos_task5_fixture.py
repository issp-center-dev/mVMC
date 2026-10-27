#!/usr/bin/env python3
"""Frozen Task 5 GC Lanczos fixtures as a C header.

The JSON files in data/ are byte copies of the frozen design fixtures
(schema v2: R1, E1, R1_G2, S2_G2).  Every binary64 input is copied from its
saved hex string; nothing is recomputed from F or the quadrature.

  generate_gc_lanczos_task5_fixture.py --output <header>
  generate_gc_lanczos_task5_fixture.py --check <header>

--check verifies both SHA-256 digests and that the checked-in header is
byte-identical to the regenerated one.
"""
import argparse
import hashlib
import json
import pathlib
import sys

sys.dont_write_bytecode = True

HERE = pathlib.Path(__file__).resolve().parent
INPUTS = HERE / "data" / "gc_lanczos_task5_inputs.json"
EVIDENCE = HERE / "data" / "gc_lanczos_task5_evidence.json"
INPUTS_SHA256 = (
    "d7b840e97c73df0a1bd429d2e8dca37814608dd9046c11745f2195e7bed7b0a3")
EVIDENCE_SHA256 = (
    "d4278eb43d7c2cc8d5823ee4127e6f8bb0dc3767d3fd2dbf9ad14d33b101996b")
SYMMETRY = {"general": 0, "singlet": 1, "triplet": 2, "zero": 3}
VALUES = ("Psi", "HPsi", "H2Psi", "Phi", "HPhi")
CASES = (("R1_base_zero_rescue", "r1", "R1"),
         ("E1_finite_spin_exponential", "e1", "E1"),
         ("R1_G2_empty_rank", "r1_g2", "R1_G2"),
         ("S2_G2_mask_mutation", "s2_g2", "S2_G2"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def literal(text):
    """A C hex-float literal with the exact value of a saved hex string."""
    value = float.fromhex(text)
    if value == 0.0:
        return "-0.0" if text.startswith("-") else "0.0"
    return value.hex()


def neg(text):
    return text[1:] if text.startswith("-") else "-" + text


def pair(value):
    return "{%s, %s}" % (literal(value[0]), literal(value[1]))


def scaled_rows(reference):
    rows = []
    for mask in sorted(reference, key=int):
        for index, name in enumerate(VALUES):
            entry = reference[mask][name]
            if entry["state"] == "EXACT_ZERO":
                rows.append("  {%s, %d, 0, {0.0, 0.0}, 0.0}" %
                            (mask, index))
            else:
                if entry["state"] != "FINITE_NONZERO":
                    raise SystemExit("unexpected state %r" % entry["state"])
                rows.append("  {%s, %d, 1, %s, %s}" % (
                    mask, index, pair(entry["phase_hex"]),
                    literal(entry["log_abs_hex"])))
    return rows


def emit_case(case, tag, evidence, lines):
    L = case["L"]
    norbital = 2 * L
    slater = []
    for qp in case["qp"]:
        matrix = [["0.0", "0.0"] for _ in range(norbital * norbital)]
        for i, j, value in qp["upper"]:
            if not 0 <= i < j < norbital:
                raise SystemExit("%s: entry outside the upper triangle" % tag)
            matrix[i * norbital + j] = list(value)
            matrix[j * norbital + i] = [neg(value[0]), neg(value[1])]
        slater.extend(matrix)
    lines.append("TASK5_DATA double task5_%s_slater[%d][2] = {" %
                 (tag, len(slater)))
    lines.append(",\n".join("  " + pair(v) for v in slater))
    lines.append("};")
    lines.append("TASK5_DATA double task5_%s_weight[%d][2] = {" %
                 (tag, len(case["qp"])))
    lines.append(",\n".join("  " + pair(qp["weight"]) for qp in case["qp"]))
    lines.append("};")
    terms = case["hamiltonian"]
    for term in terms:
        if len(term["operators"]) != 2:
            raise SystemExit("%s: only two-operator terms" % tag)
    lines.append("TASK5_DATA Task5Term task5_%s_terms[%d] = {" %
                 (tag, len(terms)))
    lines.append(",\n".join(
        "  {%s, {{%d, %d}, {%d, %d}}}" % (
            pair(t["coefficient"]), t["operators"][0][0],
            int(t["operators"][0][1]), t["operators"][1][0],
            int(t["operators"][1][1])) for t in terms))
    lines.append("};")
    reference = evidence.get("scaled_reference")
    if reference:
        rows = scaled_rows(reference)
        lines.append("TASK5_DATA Task5Scaled task5_%s_scaled[%d] = {" %
                     (tag, len(rows)))
        lines.append(",\n".join(rows))
        lines.append("};")
    charge = case["charge"]
    if (charge["kind"] != "uniform_gutzwiller" or
            charge["exponent_count"] != "number_of_doubly_occupied_sites"):
        raise SystemExit("%s: unsupported charge" % tag)
    lines.append("TASK5_DATA Task5Case task5_%s = {" % tag)
    lines.append('  "%s", %d, %d, %d, %d, %d,' % (
        case["id"], L, case["S"], case["G"], SYMMETRY[case["symmetry"]],
        len(case["qp"])))
    lines.append("  task5_%s_slater, task5_%s_weight," % (tag, tag))
    lines.append("  {%s, %s}," % (pair(case["c"][0]), pair(case["c"][1])))
    lines.append("  %s," % pair(charge["parameter"]))
    lines.append("  %d, task5_%s_terms," % (len(terms), tag))
    if reference:
        lines.append("  %d, task5_%s_scaled" % (len(scaled_rows(reference)),
                                                tag))
    else:
        lines.append("  0, NULL")
    lines.append("};")


def build():
    inputs = json.loads(INPUTS.read_text())
    evidence = json.loads(EVIDENCE.read_text())
    if inputs.get("schema_version") != 2 or \
            evidence.get("inputs_sha256") != INPUTS_SHA256:
        raise SystemExit("fixture schema or evidence binding changed")
    by_id = {case["id"]: case for case in inputs["cases"]}
    lines = [
        "/* Generated by test/python/generate_gc_lanczos_task5_fixture.py; "
        "do not edit.",
        " * inputs   SHA-256 %s" % INPUTS_SHA256,
        " * evidence SHA-256 %s" % EVIDENCE_SHA256,
        " * Frozen binary64 Task 5 fixtures (source commit %s)." %
        inputs["source_commit"][:12],
        " * Slater: QP-major, full row-major (2L)^2, lower = -upper, diagonal 0;",
        " * orbitals up0..up(L-1), down0..down(L-1).  Operators are listed",
        " * left to right and act right to left.  Not a production QP dump. */",
        "#ifndef GC_LANCZOS_TASK5_FIXTURE_H",
        "#define GC_LANCZOS_TASK5_FIXTURE_H",
        "",
        "#include <stddef.h>",
        "",
        "/* Each test uses a subset of the tables. */",
        "#if defined(__GNUC__) || defined(__clang__)",
        "#define TASK5_DATA static const __attribute__((unused))",
        "#else",
        "#define TASK5_DATA static const",
        "#endif",
        "",
        "typedef struct {",
        "  double coefficient[2];",
        "  int op[2][2]; /* {orbital, create} */",
        "} Task5Term;",
        "",
        "/* value: 0 Psi, 1 HPsi, 2 H2Psi, 3 Phi, 4 HPhi; finite 0 = EXACT_ZERO. */",
        "typedef struct {",
        "  unsigned configuration;",
        "  int value;",
        "  int finite;",
        "  double phase[2];",
        "  double log_abs;",
        "} Task5Scaled;",
        "",
        "typedef struct {",
        "  const char *id;",
        "  int site_count, total_spin, gauss_count, symmetry, qp_total;",
        "  const double (*slater)[2];",
        "  const double (*weight)[2];",
        "  double coefficient[2][2];",
        "  double gutzwiller[2];",
        "  size_t term_count;",
        "  const Task5Term *terms;",
        "  size_t scaled_count;",
        "  const Task5Scaled *scaled;",
        "} Task5Case;",
        ""]
    for case_id, tag, key in CASES:
        emit_case(by_id[case_id], tag, evidence[key], lines)
        lines.append("")
    complex_control = evidence["R1"]["complex_c1_control"]
    rows = scaled_rows(complex_control["scaled_reference"])
    lines.append("/* R1C: R1 inputs with c1 = %s. */" %
                 " ".join(complex_control["c1"]))
    lines.append("TASK5_DATA double task5_r1c_c1[2] = %s;" %
                 pair(complex_control["c1"]))
    lines.append("TASK5_DATA Task5Scaled task5_r1c_scaled[%d] = {" %
                 len(rows))
    lines.append(",\n".join(rows))
    lines.append("};")
    lines.append("")
    draws = evidence["draw_enumeration"]
    lines.append("/* Structured initial draw: configuration probabilities "
                 "num/den. */")
    lines.append("typedef struct { unsigned configuration; unsigned "
                 "numerator, denominator; } Task5DrawProbability;")
    lines.append("typedef struct { int site_count, total_spin, "
                 "electron_count; unsigned leaves; size_t count; "
                 "const Task5DrawProbability *probability; } Task5Draw;")
    for index, item in enumerate(draws):
        entries = []
        for mask in sorted(item["probabilities"], key=int):
            num, den = item["probabilities"][mask].split("/")
            entries.append("{%s, %s, %s}" % (mask, num, den))
        lines.append("TASK5_DATA Task5DrawProbability task5_draw%d[%d] = {"
                     % (index, len(entries)))
        lines.append("  " + ",\n  ".join(entries))
        lines.append("};")
    lines.append("TASK5_DATA Task5Draw task5_draws[%d] = {" % len(draws))
    lines.append(",\n".join(
        "  {%d, %d, %d, %d, %d, task5_draw%d}" % (
            item["L"], item["S"], item["N"], item["rng_leaves"],
            len(item["probabilities"]), index)
        for index, item in enumerate(draws)))
    lines.append("};")
    lines.append("")
    lines.append("#endif /* GC_LANCZOS_TASK5_FIXTURE_H */")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--output")
    group.add_argument("--check")
    args = parser.parse_args()
    if digest(INPUTS) != INPUTS_SHA256 or digest(EVIDENCE) != EVIDENCE_SHA256:
        print("frozen fixture JSON digest mismatch")
        return 1
    text = build()
    if args.output:
        pathlib.Path(args.output).write_text(text)
        return 0
    if pathlib.Path(args.check).read_text() != text:
        print("%s differs from the regenerated header" % args.check)
        return 1
    print("GC Lanczos Task 5 fixture header: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
