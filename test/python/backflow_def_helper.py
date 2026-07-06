from __future__ import print_function

import argparse
from collections import namedtuple
import os
import shutil
import sys
import tempfile


SEPARATOR = "===================="

BackflowDefinition = namedtuple(
    "BackflowDefinition",
    [
        "length",
        "nrange",
        "nzbf",
        "n_backflow_idx",
        "nrange_idx",
        "n_bf_idx_total",
        "n_proj_bf",
        "positions",
        "range_rows",
        "bf_rows",
        "opt_rows",
        "rangebf_text",
        "bf_text",
        "bf_compact_text",
    ],
)


def _require(condition, message):
    if not condition:
        raise AssertionError(message)


def _body_header(keyword, values):
    count_line = "{} {}".format(keyword, " ".join(str(v) for v in values))
    return [
        SEPARATOR,
        count_line,
        SEPARATOR,
        SEPARATOR,
        SEPARATOR,
        SEPARATOR,
        count_line,
        SEPARATOR,
        SEPARATOR,
        SEPARATOR,
    ]


def _render_rangebf(nrange, nzbf, rows):
    lines = _body_header("Nrange", [nrange, nzbf])
    lines.extend("{} {} {}".format(i, j, shell) for i, j, shell in rows)
    return "\n".join(lines) + "\n"


def _render_bf(n_backflow_idx, bf_rows, opt_rows):
    lines = _body_header("NBackFlowIdx", [n_backflow_idx])
    lines.extend("{} {} {} {} {}".format(i, j, x0, x1, group) for i, j, x0, x1, group in bf_rows)
    lines.extend("{} {}".format(idx, flag) for idx, flag in opt_rows)
    return "\n".join(lines) + "\n"


def _render_bf_compact(n_backflow_idx, opt_rows):
    lines = _body_header("NBackFlowIdx", [n_backflow_idx])
    lines.append("BFCompact 1")
    lines.extend("{} {}".format(idx, flag) for idx, flag in opt_rows)
    return "\n".join(lines) + "\n"


def build_chain_nn_backflow(length=4, optimize=False):
    """Build minimal 1D PBC nearest-neighbor BackFlow definition text."""
    if length < 4:
        raise ValueError("length must be at least 4 for the NN helper")

    nrange = 3
    nzbf = 2
    n_backflow_idx = 1
    nrange_idx = 3 * (nrange - 1) // nzbf + 1
    n_bf_idx_total = nrange_idx * (nrange_idx + 1) // 2
    n_proj_bf = n_bf_idx_total * n_backflow_idx

    positions = {}
    range_rows = []
    for site in range(length):
        positions[site] = [site, (site - 1) % length, (site + 1) % length]
        for pos in positions[site]:
            shell = 0 if pos == site else 1
            range_rows.append((site, pos, shell))

    bf_rows = []
    for i in range(length):
        for j in range(length):
            for x0 in positions[i]:
                for x1 in positions[j]:
                    bf_rows.append((i, j, x0, x1, 0))

    flag = 1 if optimize else 0
    opt_rows = [(idx, flag) for idx in range(n_proj_bf)]

    definition = BackflowDefinition(
        length=length,
        nrange=nrange,
        nzbf=nzbf,
        n_backflow_idx=n_backflow_idx,
        nrange_idx=nrange_idx,
        n_bf_idx_total=n_bf_idx_total,
        n_proj_bf=n_proj_bf,
        positions=positions,
        range_rows=range_rows,
        bf_rows=bf_rows,
        opt_rows=opt_rows,
        rangebf_text=_render_rangebf(nrange, nzbf, range_rows),
        bf_text=_render_bf(n_backflow_idx, bf_rows, opt_rows),
        bf_compact_text=_render_bf_compact(n_backflow_idx, opt_rows),
    )
    verify_definition(definition, expected_opt_flag=flag)
    return definition


def write_chain_nn_backflow(output_dir, length=4, optimize=False, compact=False):
    definition = build_chain_nn_backflow(length=length, optimize=optimize)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    range_path = os.path.join(output_dir, "rangebf.def")
    bf_path = os.path.join(output_dir, "bf.def")
    with open(range_path, "w") as fp:
        fp.write(definition.rangebf_text)
    with open(bf_path, "w") as fp:
        fp.write(definition.bf_compact_text if compact else definition.bf_text)
    return range_path, bf_path


def _verify_header(lines, keyword, values):
    count_line = "{} {}".format(keyword, " ".join(str(v) for v in values))
    _require(len(lines) >= 10, "{} file is shorter than the BF header".format(keyword))
    _require(lines[1] == count_line, "{} line 2 count mismatch".format(keyword))
    _require(lines[6] == count_line, "{} line 7 count mismatch".format(keyword))
    for idx in (0, 2, 3, 4, 5, 7, 8, 9):
        _require(lines[idx] == SEPARATOR, "{} header line {} mismatch".format(keyword, idx + 1))


def verify_definition(definition, expected_opt_flag=0):
    length = definition.length
    nrange = definition.nrange
    positions = definition.positions

    _require(definition.nrange == 3, "minimal helper must use Nrange=3")
    _require(definition.nzbf == 2, "minimal helper must use NzBF=2")
    _require(definition.n_backflow_idx == 1, "minimal helper must use one BF group")
    _require(definition.nrange_idx == 4, "unexpected NrangeIdx")
    _require(definition.n_bf_idx_total == 10, "unexpected NBFIdxTotal")
    _require(definition.n_proj_bf == 10, "unexpected NProjBF")

    range_lines = definition.rangebf_text.splitlines()
    bf_lines = definition.bf_text.splitlines()
    compact_lines = definition.bf_compact_text.splitlines()
    _verify_header(range_lines, "Nrange", [definition.nrange, definition.nzbf])
    _verify_header(bf_lines, "NBackFlowIdx", [definition.n_backflow_idx])
    _verify_header(compact_lines, "NBackFlowIdx", [definition.n_backflow_idx])
    _require(len(range_lines) == 10 + length * nrange, "rangebf row count mismatch")
    _require(len(bf_lines) == 10 + length * length * nrange * nrange + definition.n_proj_bf,
             "bf row count mismatch")
    _require(compact_lines[10] == "BFCompact 1", "compact bf.def marker mismatch")
    _require(len(compact_lines) == 11 + definition.n_proj_bf, "compact bf.def row count mismatch")

    range_seen = set()
    per_center = dict((site, 0) for site in range(length))
    for i, j, shell in definition.range_rows:
        _require(0 <= i < length and 0 <= j < length, "BFRange site index out of range")
        _require(shell in (0, 1), "BFRange shell out of MVP range")
        _require((i == j and shell == 0) or (i != j and shell == 1),
                 "BFRange self/shell contract mismatch")
        _require((i, j) not in range_seen, "duplicated BFRange pair")
        range_seen.add((i, j))
        per_center[i] += 1
    for site in range(length):
        _require(per_center[site] == nrange, "BFRange center row count mismatch")
        _require(positions[site] == [site, (site - 1) % length, (site + 1) % length],
                 "unexpected PosBF order for site {}".format(site))

    bf_seen = set()
    for i, j, x0, x1, group in definition.bf_rows:
        _require(0 <= i < length and 0 <= j < length, "BF row base site out of range")
        _require(0 <= x0 < length and 0 <= x1 < length, "BF row range site out of range")
        _require(group == 0, "MVP BF group must be 0")
        _require(x0 in positions[i] and x1 in positions[j], "BF row uses site outside PosBF")
        key = (i, j, x0, x1)
        _require(key not in bf_seen, "duplicated BF row")
        bf_seen.add(key)
    _require(len(bf_seen) == length * length * nrange * nrange, "BF coverage mismatch")

    _require(len(definition.opt_rows) == definition.n_proj_bf, "OptFlag row count mismatch")
    for expected_idx, (idx, flag) in enumerate(definition.opt_rows):
        _require(idx == expected_idx, "OptFlag idx must be sequential")
        _require(flag == expected_opt_flag, "OptFlag flag mismatch")


def self_test():
    definition = build_chain_nn_backflow(length=4, optimize=False)
    verify_definition(definition, expected_opt_flag=0)
    optimized = build_chain_nn_backflow(length=4, optimize=True)
    verify_definition(optimized, expected_opt_flag=1)

    tmpdir = tempfile.mkdtemp(prefix="backflow_def_helper_")
    try:
        range_path, bf_path = write_chain_nn_backflow(tmpdir, length=4, optimize=False)
        with open(range_path) as fp:
            _require(fp.read() == definition.rangebf_text, "written rangebf.def mismatch")
        with open(bf_path) as fp:
            _require(fp.read() == definition.bf_text, "written bf.def mismatch")
        range_path, bf_path = write_chain_nn_backflow(tmpdir, length=4, optimize=True, compact=True)
        with open(range_path) as fp:
            _require(fp.read() == optimized.rangebf_text, "written compact rangebf.def mismatch")
        with open(bf_path) as fp:
            _require(fp.read() == optimized.bf_compact_text, "written compact bf.def mismatch")
    finally:
        shutil.rmtree(tmpdir)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Generate minimal BackFlow definition files for tests.")
    parser.add_argument("--length", type=int, default=4,
                        help="1D PBC chain length for the NN helper (default: 4)")
    parser.add_argument("--output", help="directory where rangebf.def and bf.def are written")
    parser.add_argument("--optimize", action="store_true",
                        help="write BF OptFlag=1 instead of the identity-test default OptFlag=0")
    parser.add_argument("--compact", action="store_true",
                        help="write compact bf.def with BFCompact marker and OptFlag rows only")
    parser.add_argument("--self-test", action="store_true",
                        help="run helper contract checks")
    args = parser.parse_args(argv)

    if args.self_test:
        self_test()
        print("backflow_def_helper: self-test passed")

    if args.output:
        range_path, bf_path = write_chain_nn_backflow(
            args.output,
            length=args.length,
            optimize=args.optimize,
            compact=args.compact,
        )
        print("wrote {}".format(range_path))
        print("wrote {}".format(bf_path))

    if not args.self_test and not args.output:
        parser.error("specify --self-test or --output DIR")

    return 0


if __name__ == "__main__":
    sys.exit(main())
