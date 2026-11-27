"""Check sampler invariance under uneven and empty MPI QP partitions.

Three distinct translations are split over 1, 2 and 4 ranks in one sampler
group. The seed and sample count therefore stay fixed. Compare scientific
outputs, including every optimization iteration, with and without NStore.
"""

import argparse
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from runtest_nstore_nsplit import update_modpara


def numbers(path):
    rows = [[float(x) for x in line.split()] for line in path.read_text().splitlines()
            if line.strip()]
    if not rows or any(not math.isfinite(x) for row in rows for x in row):
        raise AssertionError("empty or non-finite output: {}".format(path))
    return rows


def compare(left, right, label):
    if [len(row) for row in left] != [len(row) for row in right]:
        raise AssertionError("{}: output shape differs".format(label))
    largest = 0.0
    for i, (a_row, b_row) in enumerate(zip(left, right)):
        for j, (a, b) in enumerate(zip(a_row, b_row)):
            error = abs(a - b) / max(1.0, abs(a), abs(b))
            largest = max(largest, error)
            if error > 1.0e-8:
                raise AssertionError("{} row {} col {}: {} != {} (scaled error {})"
                                     .format(label, i, j, a, b, error))
    return largest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complex", action="store_true")
    parser.add_argument("--fsz", action="store_true")
    parser.add_argument("--binary", type=Path,
                        default=Path("../../src/mVMC/vmc.out"))
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--np-flag", default="-np")
    args = parser.parse_args()
    binary = args.binary.resolve()
    fixture = Path("data") / ("HubbardChain_fsz_SpinJastrow" if args.fsz
                              else "HubbardChain_SpinJastrow")
    Path("work").mkdir(exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="qp_split_{}{}_".format(
        "complex" if args.complex else "real", "_fsz" if args.fsz else ""),
        dir="work")).resolve()
    print("QP split outputs: {}".format(root), flush=True)
    max_error = 0.0
    for exchange in (0, 1):
        reference = None
        for nproc in (1, 2, 4):
            for store in (0, 1):
                case = root / "exchange{}_np{}_store{}".format(exchange, nproc, store)
                case.mkdir()
                for path in fixture.glob("*.def"):
                    shutil.copy2(path, case / path.name)
                for path in case.glob("*.def"):
                    text = re.sub(r"(?m)^(ComplexType\s+)\d+",
                                  lambda m: m[1] + str(int(args.complex)),
                                  path.read_text())
                    path.write_text(text)
                update_modpara(str(case / "modpara.def"), {
                    "NMPTrans": 3, "NSPGaussLeg": 1, "NSplitSize": nproc,
                    "NStore": store, "NExUpdatePath": exchange,
                    "NSROptItrStep": 4, "NSROptItrSmp": 2, "NSRCG": 0,
                    "NVMCWarmUp": 4, "NVMCInterval": 1, "NVMCSample": 24,
                })
                lines = ["====================", "NQPTrans 3"] + ["===================="] * 3
                lines += ["{} {}".format(q, weight)
                          for q, weight in enumerate((0.6, 0.3, 0.1))]
                lines += ["{} {} {} 1".format(q, i, (i + q) % 6)
                          for q in range(3) for i in range(6)]
                (case / "qptransidx.def").write_text("\n".join(lines) + "\n")
                command = [args.mpiexec, args.np_flag, str(nproc), str(binary),
                           "-e", "namelist.def"]
                with (case / "run.log").open("w") as log:
                    result = subprocess.run(command, cwd=str(case), stdout=log,
                                            stderr=subprocess.STDOUT, timeout=120)
                if result.returncode:
                    raise AssertionError("{} failed: see {}".format(command, case / "run.log"))
                output = {name: numbers(case / "output" / name)
                          for name in ("zvo_out_001.dat", "zvo_var_001.dat", "zqp_opt.dat")}
                if reference is None:
                    reference = output
                else:
                    for name in output:
                        max_error = max(max_error, compare(reference[name], output[name],
                                                          "{}/{}".format(case.name, name)))
    print("PASS: 12 runs, uneven/empty QP slices, max scaled error {:.3e}".format(max_error))


if __name__ == "__main__":
    main()
