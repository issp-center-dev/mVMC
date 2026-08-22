"""Manual, non-CI structural scaling probe for second power-Lanczos.

The qualitative checks below verify that the expected H^3 work is exercised
as NTransfer grows.  The projection grid records both the identity fast path
and the NQPFull>1 direct-contraction path.  These measurements are not a
portable performance-regression gate: absolute ratios depend on the compiler,
architecture, and linear-algebra stack.  Establish quantitative baselines on
the same Linux/HPC environment.
"""

from __future__ import print_function

import argparse
import hashlib
import json
import math
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile


TIMER_PATTERN = re.compile(r"\[(41|95)\]\s+([0-9.eE+-]+)\s*$")


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_standard_input(path, size, samples, warmup):
    with open(path, "w") as destination:
        destination.write(
            """L = {size}
Lsub = 2
model = "Hubbard"
lattice = "chain"
U = 4
t = 1
ncond = {size}
2Sz = 0
NSROptItrStep = 1
NVMCSample = {samples}
NVMCWarmUp = {warmup}
NVMCCalMode = 1
NLanczosMode = 1
NLanczosEstimatorMode = 1
RndSeed = 2468
""".format(size=size, samples=samples, warmup=warmup)
        )


def write_momentum_projection(path, size, nmptrans):
    """Write a k=0 projector for a cyclic translation subgroup."""
    if size % nmptrans != 0:
        raise RuntimeError(
            "NMPTrans={} must divide L={} to form a cyclic subgroup".format(
                nmptrans, size
            )
        )
    stride = size // nmptrans
    with open(path, "w") as destination:
        destination.write("=============================================\n")
        destination.write("NQPTrans          {}\n".format(nmptrans))
        destination.write("=============================================\n")
        destination.write(
            "======== TrIdx_TrWeight_and_TrIdx_i_xi ======\n"
        )
        destination.write("=============================================\n")
        for index in range(nmptrans):
            destination.write(
                "{:5d} {:25.15f} {:25.15f}\n".format(index, 1.0, 0.0)
            )
        for index in range(nmptrans):
            displacement = index * stride
            for origin in range(size):
                destination.write(
                    "{:5d} {:5d} {:5d} {:5d}\n".format(
                        index,
                        origin,
                        (origin + displacement) % size,
                        1,
                    )
                )


def replace_parameters(path, updates):
    with open(path) as source:
        lines = source.readlines()
    remaining = dict(updates)
    result = []
    for line in lines:
        columns = line.split()
        if columns and columns[0] in remaining:
            key = columns[0]
            result.append("{:<16s} {}\n".format(key, remaining.pop(key)))
        else:
            result.append(line)
    for key in sorted(remaining):
        result.append("{:<16s} {}\n".format(key, remaining[key]))
    with open(path, "w") as destination:
        destination.writelines(result)


def read_definition_count(path, keyword):
    with open(path) as source:
        for line in source:
            columns = line.split()
            if len(columns) >= 2 and columns[0] == keyword:
                return int(columns[1])
    raise RuntimeError("{} is missing from {}".format(keyword, path))


def write_initial_parameter(path, workdir):
    gutzwiller = read_definition_count(
        os.path.join(workdir, "gutzwilleridx.def"), "NGutzwillerIdx"
    )
    jastrow = read_definition_count(
        os.path.join(workdir, "jastrowidx.def"), "NJastrowIdx"
    )
    orbital = read_definition_count(
        os.path.join(workdir, "orbitalidx.def"), "NOrbitalIdx"
    )
    values = [0.0] * 6
    for index in range(gutzwiller + jastrow):
        values.extend(
            [
                0.08 * math.sin(0.31 * (index + 1)),
                0.0,
                0.0,
            ]
        )
    for index in range(orbital):
        values.extend(
            [
                math.sin(0.37 * (index + 1))
                + 0.23 * math.cos(0.19 * (index + 2)),
                0.0,
                0.0,
            ]
        )
    with open(path, "w") as destination:
        destination.write(
            " ".join("{:.18e}".format(value) for value in values)
        )
        destination.write("\n")


def read_transfer_count(path):
    return read_definition_count(path, "NTransfer")


def read_timers(path):
    timers = {}
    with open(path) as source:
        for line in source:
            match = TIMER_PATTERN.search(line)
            if match is not None:
                timers[int(match.group(1))] = float(match.group(2))
    missing = sorted(set((41, 95)) - set(timers))
    if missing:
        raise RuntimeError(
            "{} is missing timer(s) {}".format(path, ", ".join(map(str, missing)))
        )
    if timers[41] <= 0.0 or timers[95] <= 0.0:
        raise RuntimeError("timers must be positive: {}".format(timers))
    return timers


def run_checked(command, workdir, environment):
    process = subprocess.run(
        command,
        cwd=workdir,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if process.returncode != 0:
        raise RuntimeError(
            "command failed ({}):\n{}".format(
                " ".join(command), process.stdout
            )
        )
    return process.stdout


def projection_label(nspgaussleg, nmptrans):
    return "sp{}-mp{}".format(nspgaussleg, nmptrans)


def prepare_case(
    root,
    vmcdry,
    size,
    samples,
    warmup,
    nspgaussleg,
    nmptrans,
    environment,
):
    case = os.path.join(
        root,
        projection_label(nspgaussleg, nmptrans),
        "L{}".format(size),
        "template",
    )
    os.makedirs(case)
    standard_input = os.path.join(case, "StdFace.def")
    write_standard_input(standard_input, size, samples, warmup)
    run_checked([vmcdry, standard_input], case, environment)
    replace_parameters(
        os.path.join(case, "modpara.def"),
        {
            "NLanczosMode": "1",
            "NLanczosEstimatorMode": "1",
            "NLanczosStep": "2",
            "NVMCCalMode": "1",
            "NVMCSample": str(samples),
            "NVMCWarmUp": str(warmup),
            "NExUpdatePath": "0",
            "NSPGaussLeg": str(nspgaussleg),
            "NMPTrans": str(nmptrans),
            "NDataQtySmp": "1",
        },
    )
    write_momentum_projection(
        os.path.join(case, "qptransidx.def"), size, nmptrans
    )
    write_initial_parameter(os.path.join(case, "initial.dat"), case)
    transfer_count = read_transfer_count(os.path.join(case, "trans.def"))
    expected_transfer_count = 4 * size
    if transfer_count != expected_transfer_count:
        raise RuntimeError(
            "L={} generated NTransfer={}, expected {}".format(
                size, transfer_count, expected_transfer_count
            )
        )
    return case, transfer_count


def measure_case(
    root,
    template,
    vmc,
    size,
    repeats,
    nspgaussleg,
    nmptrans,
    environment,
):
    measurements = []
    for repeat in range(repeats):
        workdir = os.path.join(
            root,
            projection_label(nspgaussleg, nmptrans),
            "L{}".format(size),
            "repeat{}".format(repeat + 1),
        )
        shutil.copytree(template, workdir)
        replace_parameters(
            os.path.join(workdir, "modpara.def"),
            {"RndSeed": str(2468 + repeat)},
        )
        run_checked(
            [vmc, "-e", "namelist.def", "initial.dat"],
            workdir,
            environment,
        )
        timers = read_timers(
            os.path.join(workdir, "output", "zvo_CalcTimer.dat")
        )
        measurements.append(
            {
                "timer41_seconds": timers[41],
                "timer95_seconds": timers[95],
                "timer95_over_timer41": timers[95] / timers[41],
            }
        )
    return measurements


def summarize(
    size, transfer_count, nspgaussleg, nmptrans, measurements
):
    timer41 = statistics.median(
        row["timer41_seconds"] for row in measurements
    )
    timer95 = statistics.median(
        row["timer95_seconds"] for row in measurements
    )
    ratio = statistics.median(
        row["timer95_over_timer41"] for row in measurements
    )
    return {
        "L": size,
        "NTransfer": transfer_count,
        "NSPGaussLeg": nspgaussleg,
        "NMPTrans": nmptrans,
        "NQPFull": nspgaussleg * nmptrans,
        "timer41_seconds_median": timer41,
        "timer95_seconds_median": timer95,
        "timer95_over_timer41_median": ratio,
        "ratio_per_transfer": ratio / transfer_count,
        "measurements": measurements,
    }


def projection_groups(rows):
    groups = {}
    for row in rows:
        key = (row["NSPGaussLeg"], row["NMPTrans"])
        groups.setdefault(key, []).append(row)
    return [
        (key, sorted(group, key=lambda row: row["L"]))
        for key, group in sorted(groups.items())
    ]


def check_scaling(rows):
    # Deliberately check only qualitative growth.  Cross-platform absolute
    # timing bands belong to a separately recorded Linux/HPC baseline.
    for (nspgaussleg, nmptrans), group in projection_groups(rows):
        if len(group) < 2:
            raise RuntimeError(
                "--check-scaling needs at least two sizes for "
                "NSPGaussLeg={}, NMPTrans={}".format(
                    nspgaussleg, nmptrans
                )
            )
        ratios = [
            row["timer95_over_timer41_median"] for row in group
        ]
        for previous, current in zip(ratios, ratios[1:]):
            if current < 0.9 * previous:
                raise RuntimeError(
                    "Timer[95]/Timer[41] is not monotone within 10% for "
                    "NSPGaussLeg={}, NMPTrans={}: {}".format(
                        nspgaussleg, nmptrans, ratios
                    )
                )
        if ratios[-1] < 1.3 * ratios[0]:
            raise RuntimeError(
                "endpoint Timer[95]/Timer[41] growth is below 1.3x for "
                "NSPGaussLeg={}, NMPTrans={}: {}".format(
                    nspgaussleg,
                    nmptrans,
                    ratios[-1] / ratios[0],
                )
            )


def print_markdown(result):
    print("# Second-Lanczos cost probe")
    print()
    print(
        "platform: `{}`; OMP/BLAS threads: 1; samples: {}; warmup: {}; "
        "repeats: {}".format(
            result["platform"],
            result["samples"],
            result["warmup"],
            result["repeats"],
        )
    )
    print()
    print(
        "| L | NSPGaussLeg | NMPTrans | NQPFull | NTransfer "
        "| Timer[41] median (s) | Timer[95] median (s) "
        "| [95]/[41] | ratio/NTransfer |"
    )
    print("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in result["rows"]:
        print(
            "| {L} | {NSPGaussLeg} | {NMPTrans} | {NQPFull} "
            "| {NTransfer} | {timer41_seconds_median:.6g} "
            "| {timer95_seconds_median:.6g} "
            "| {timer95_over_timer41_median:.6g} "
            "| {ratio_per_transfer:.6g} |".format(**row)
        )
    for (nspgaussleg, nmptrans), group in projection_groups(
        result["rows"]
    ):
        if len(group) >= 2:
            endpoint = (
                group[-1]["timer95_over_timer41_median"]
                / group[0]["timer95_over_timer41_median"]
            )
            print()
            print(
                "Endpoint growth of Timer[95]/Timer[41] for "
                "NSPGaussLeg={}, NMPTrans={} (NQPFull={}): {:.6g}x.".format(
                    nspgaussleg,
                    nmptrans,
                    nspgaussleg * nmptrans,
                    endpoint,
                )
            )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Manual, non-CI Timer[95]/Timer[41] scaling and projection "
            "cost probe."
        )
    )
    parser.add_argument("--vmc", required=True, help="path to vmc.out")
    parser.add_argument("--vmcdry", required=True, help="path to vmcdry.out")
    parser.add_argument("--sizes", type=int, nargs="+", default=[4, 6, 8])
    parser.add_argument(
        "--nspgaussleg",
        type=int,
        nargs="+",
        default=[1],
        help=(
            "spin-projection quadrature counts; all values are crossed "
            "with --nmptrans (default: 1)"
        ),
    )
    parser.add_argument(
        "--nmptrans",
        type=int,
        nargs="+",
        default=[1],
        help=(
            "cyclic momentum-projection subgroup sizes; each value must "
            "divide every L (default: 1)"
        ),
    )
    parser.add_argument("--samples", type=int, default=4096)
    parser.add_argument("--warmup", type=int, default=64)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--work-parent", default=None)
    parser.add_argument("--keep-work", action="store_true")
    parser.add_argument(
        "--check-scaling",
        action="store_true",
        help=(
            "check qualitative growth only; this is not a portable "
            "performance-regression baseline"
        ),
    )
    parser.add_argument("--json-output", default=None)
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    if (
        arguments.samples <= 0
        or arguments.warmup < 0
        or arguments.repeats <= 0
        or not arguments.sizes
        or any(size <= 0 or size % 2 != 0 for size in arguments.sizes)
        or any(value <= 0 for value in arguments.nspgaussleg)
        or any(value <= 0 for value in arguments.nmptrans)
    ):
        raise RuntimeError(
            "sizes must be positive even integers; projection counts and "
            "samples/repeats must be positive; warmup must be non-negative"
        )
    sizes = sorted(set(arguments.sizes))
    nspgaussleg_values = sorted(set(arguments.nspgaussleg))
    nmptrans_values = sorted(set(arguments.nmptrans))
    for size in sizes:
        for nmptrans in nmptrans_values:
            if size % nmptrans != 0:
                raise RuntimeError(
                    "NMPTrans={} must divide L={}; choose a cyclic "
                    "translation subgroup for every measured size".format(
                        nmptrans, size
                    )
                )
    if arguments.check_scaling and len(sizes) < 2:
        raise RuntimeError("--check-scaling requires at least two sizes")
    vmc = os.path.abspath(arguments.vmc)
    vmcdry = os.path.abspath(arguments.vmcdry)
    for binary in (vmc, vmcdry):
        if not os.path.isfile(binary) or not os.access(binary, os.X_OK):
            raise RuntimeError("executable not found: {}".format(binary))

    environment = os.environ.copy()
    environment.update(
        {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
        }
    )
    workroot = tempfile.mkdtemp(
        prefix="lanczos2-cost-probe-", dir=arguments.work_parent
    )
    try:
        rows = []
        for nspgaussleg in nspgaussleg_values:
            for nmptrans in nmptrans_values:
                for size in sizes:
                    template, transfer_count = prepare_case(
                        workroot,
                        vmcdry,
                        size,
                        arguments.samples,
                        arguments.warmup,
                        nspgaussleg,
                        nmptrans,
                        environment,
                    )
                    measurements = measure_case(
                        workroot,
                        template,
                        vmc,
                        size,
                        arguments.repeats,
                        nspgaussleg,
                        nmptrans,
                        environment,
                    )
                    rows.append(
                        summarize(
                            size,
                            transfer_count,
                            nspgaussleg,
                            nmptrans,
                            measurements,
                        )
                    )
        if arguments.check_scaling:
            check_scaling(rows)
        result = {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "probe_sha256": sha256_file(os.path.abspath(__file__)),
            "vmc_sha256": sha256_file(vmc),
            "vmcdry_sha256": sha256_file(vmcdry),
            "threads": {
                "OMP_NUM_THREADS": 1,
                "OPENBLAS_NUM_THREADS": 1,
                "MKL_NUM_THREADS": 1,
                "VECLIB_MAXIMUM_THREADS": 1,
            },
            "projection_definition": "cyclic k=0 translation subgroup",
            "NQPOptTrans": 1,
            "samples": arguments.samples,
            "warmup": arguments.warmup,
            "repeats": arguments.repeats,
            "nspgaussleg": nspgaussleg_values,
            "nmptrans": nmptrans_values,
            "rows": rows,
        }
        print_markdown(result)
        if arguments.json_output is not None:
            with open(arguments.json_output, "w") as destination:
                json.dump(result, destination, indent=2, sort_keys=True)
                destination.write("\n")
    finally:
        if arguments.keep_work:
            print("work directory kept at {}".format(workroot), file=sys.stderr)
        else:
            shutil.rmtree(workroot)
    return 0


if __name__ == "__main__":
    sys.exit(main())
