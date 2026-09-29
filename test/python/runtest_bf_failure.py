"""D-10: identical proposals/random draws and atomic recovery, including empty-QP ranks."""
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from backflow_def_helper import write_chain_nn_backflow
from runtest_bf import update_modpara, parse_norbitalidx, write_nonidentity_init_parameter


def compare_trace(reference, actual):
    assert len(reference) == len(actual), (len(reference), len(actual))
    for step, (left, right) in enumerate(zip(reference, actual)):
        for key in ("kind", "accept", "u", "ele_idx", "candidate_ele_idx", "count", "bf_count", "qp_start", "qp_end"):
            assert left[key] == right[key], (step, key, left[key], right[key])
        # Regular four-site fixture: absolute + scaled tolerance fixed before
        # comparing recovery results. Configuration/count/decision comparisons
        # above remain exact, including rejected proposals after recovery.
        for key in ("pf", "inverse", "slater"):
            assert len(left[key]) == len(right[key]), (step, key)
            for a, b in zip(left[key], right[key]):
                a, b = complex(*a), complex(*b)
                assert abs(a-b) <= 2e-11 + 1e-10*max(abs(a), abs(b)), (step, key, a, b)
        for key in ("log_committed", "log_candidate"):
            a, b = float(left[key]), float(right[key])
            assert a == b or abs(a-b) <= 2e-11 + 1e-10*max(abs(a), abs(b)), (step, key, a, b)


def main():
    mode, boundary = sys.argv[1:]
    root = Path.cwd()
    ranks = int(os.environ.get("MVMC_MPI_PROCS", "1"))
    parent = root / "work" / ("bf_failure_{}_{}_np{}".format(mode, boundary, ranks))
    if parent.exists():
        shutil.rmtree(parent)
    parent.mkdir(parents=True)
    binary = (root / "../../src/mVMC/vmc.out").resolve()
    prefix = [] if ranks == 1 else [os.environ.get("MPIRUN", "mpiexec"), "-np", str(ranks)] + os.environ.get("MVMC_MPI_ARGS", "").split()
    model = "BackFlow_Identity_" + mode.capitalize()

    def run(failure, node=False):
        work = parent / ("physical_zero" if node else (failure.replace(",", "_") if failure else "baseline"))
        shutil.copytree(root / "data" / model, work)
        write_chain_nn_backflow(str(work), antiperiodic=boundary == "ap")
        update_modpara(str(work / "modpara.def"), {
            "NMPTrans": "-1" if boundary == "ap" else "1", "NSplitSize": str(ranks),
            "NVMCWarmUp": "0", "NVMCInterval": "1", "NVMCSample": "32",
            "NExUpdatePath": "1", "RndSeed": "8271", "NDataQtySmp": "1"})
        write_nonidentity_init_parameter(str(work / "initial.def"), 10,
            parse_norbitalidx(str(work / "orbitalidx.def")), complex_orbitals=mode == "complex")
        initial = work / "initial.def"
        values = initial.read_text().split()
        for i in range(parse_norbitalidx(str(work / "orbitalidx.def"))):
            values[6 + 30 + 3*i] = "{:.18e}".format(0.4*math.sin((i+1)**2*0.731) + 0.3*math.cos((i+2)*0.331))
        if node:
            # Identity BF and diagonal pair orbitals: any one-electron hop
            # from a paired configuration has an exactly zero determinant.
            for i in range(10):
                values[6+3*i] = "1" if i == 0 else "0"
                values[7+3*i] = "0"
            for i in range(16):
                values[36+3*i] = "1" if i//4 == i%4 else "0"
                values[37+3*i] = "0"
        initial.write_text(" ".join(values) + "\n")
        env = dict(os.environ, MVMC_BF_TEST_TRACE="1", MVMC_BF_TEST_FAILURE=failure)
        command = prefix + [str(binary), "-e", "namelist.def", "initial.def"]
        proc = subprocess.run(command, cwd=work, env=env, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True, timeout=30)
        (work / "run.log").write_text(proc.stdout)
        return work, proc

    baseline, proc = run("")
    assert proc.returncode == 0, proc.stdout
    read = lambda p: [json.loads(line) for line in p.read_text().splitlines()]
    references = {p.name: read(p) for p in baseline.glob("bf_transaction_trace.rank*.jsonl")}
    assert len(references) == ranks, references.keys()
    first = next(iter(references.values()))
    assert any(r["kind"] >= 0 and r["accept"] for r in first)
    assert any(r["kind"] >= 0 and not r["accept"] for r in first)
    if ranks == 2:
        assert sorted(rows[0]["qp_end"]-rows[0]["qp_start"] for rows in references.values()) == [0, 1]

    for failure, counter in (
            ("proposal", "proposal"), ("proposal_nan", "proposal"),
            ("accept_lu", "accept"), ("accept_inverse", "accept"),
            ("accept_pf", "accept"), ("periodic", "periodic"),
            ("accept_reject", "accept")):
        work, proc = run(failure)
        assert proc.returncode == 0, (failure, proc.stdout)
        counts = re.findall(r"BackFlow recovery: proposal=(\d+) accept=(\d+) periodic=(\d+) zero=(\d+)", proc.stdout)
        assert counts and int(counts[-1][("proposal", "accept", "periodic").index(counter)]) > 0, (failure, proc.stdout)
        for name, reference in references.items():
            actual = read(work / name)
            compare_trace(reference, actual)
            if failure == "accept_reject":
                assert any(r["accept_recovery"] and not r["accept"] for r in actual)
        print("recovered", mode, boundary, ranks, failure, flush=True)

    for failure in ("proposal_negative", "proposal,full_pf", "accept_lu,full_inverse", "periodic,full_inverse"):
        work, proc = run(failure)
        assert proc.returncode != 0, (failure, "unexpected success")
        assert "Error: BackFlow" in proc.stdout, (failure, proc.stdout)
        assert "invalid numerical arguments" in proc.stdout or "checked recovery failed" in proc.stdout, proc.stdout
        print("stopped", mode, boundary, ranks, failure, flush=True)
    work, proc = run("", node=True)
    assert proc.returncode == 0, proc.stdout
    counts = re.findall(r"BackFlow recovery: proposal=(\d+) accept=(\d+) periodic=(\d+) zero=(\d+)", proc.stdout)
    assert counts and list(map(int, counts[-1][:3])) == [0, 0, 0] and int(counts[-1][3]) > 0, proc.stdout
    for path in work.glob("bf_transaction_trace.rank*.jsonl"):
        trace = read(path)
        assert len(trace) > 1
        for row in trace[1:]:
            assert not row["accept"] and float(row["log_candidate"]) == -math.inf
            for key in ("ele_idx", "count", "bf_count", "pf", "inverse", "slater", "log_committed"):
                assert row[key] == trace[0][key], key
    print("physical zero rejected with unchanged state", mode, boundary, ranks, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
