#!/usr/bin/env python3
"""Unit checks for rank-local BackFlow sampler dump layout validation."""

import contextlib
import io
import os
import tempfile

from runtest_bf import check_bf_sampler_rebuild_dump


def dump_line(world_rank, comm_rank, comm_size, qp_start, qp_end, nqp_full,
              local_value=0.0):
    return (
        "path complex world_rank {world_rank} comm_rank {comm_rank} "
        "comm_size {comm_size} qp_start {qp_start} qp_end {qp_end} "
        "qp_count {qp_count} nqp_full {nqp_full} "
        "proposals 1 accepts 1 rejects 1 "
        "max_candidate_slater_diff 0 "
        "max_candidate_pf_reldiff {local_value:.17e} "
        "max_accept_inv_reldiff {local_value:.17e} "
        "max_accept_pf_reldiff {local_value:.17e} "
        "max_reject_slater_diff 0 "
        "max_reject_inv_reldiff {local_value:.17e} "
        "row_select_adds 1 row_select_low_count_adds 1\n"
    ).format(
        world_rank=world_rank,
        comm_rank=comm_rank,
        comm_size=comm_size,
        qp_start=qp_start,
        qp_end=qp_end,
        qp_count=qp_end - qp_start,
        nqp_full=nqp_full,
        local_value=local_value,
    )


def run_case(name, rows, expect_pass, expect_zero_qp_rank=False):
    with tempfile.TemporaryDirectory(prefix="mvmc_bf_dump_parser_") as workdir:
        base = os.path.join(workdir, "sampler.dat")
        for world_rank, line in rows:
            with open(base + ".rank{}".format(world_rank), "w") as fp:
                fp.write(line)
        captured = io.StringIO()
        with contextlib.redirect_stdout(captured):
            result = check_bf_sampler_rebuild_dump(
                base, 1.0e-10, True, len(rows), expect_zero_qp_rank)
        passed = result == 0
        if passed != expect_pass:
            print("{}: expected pass={} but result was {}".format(name, expect_pass, result))
            print(captured.getvalue())
            return False
    return True


def main():
    cases = [
        (
            "valid two-QP split",
            [
                (0, dump_line(0, 0, 2, 0, 1, 2)),
                (1, dump_line(1, 1, 2, 1, 2, 2)),
            ],
            True,
            False,
        ),
        (
            "duplicate interval with missing tail",
            [
                (0, dump_line(0, 0, 2, 0, 1, 2)),
                (1, dump_line(1, 1, 2, 0, 1, 2)),
            ],
            False,
            False,
        ),
        (
            "valid zero-QP rank",
            [
                (0, dump_line(0, 0, 2, 0, 1, 1)),
                (1, dump_line(1, 1, 2, 1, 1, 1)),
            ],
            True,
            True,
        ),
        (
            "zero-QP rank with stale local value",
            [
                (0, dump_line(0, 0, 2, 0, 1, 1)),
                (1, dump_line(1, 1, 2, 1, 1, 1, 3.2e-315)),
            ],
            False,
            True,
        ),
    ]
    if not all(run_case(*case) for case in cases):
        return 1
    print("BackFlow sampler dump parser unit checks passed ({} cases).".format(len(cases)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
