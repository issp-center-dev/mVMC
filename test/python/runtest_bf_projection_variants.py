from __future__ import print_function

import math
import os
import subprocess
import sys


MODEL = "BackFlow_Optimization_Complex_MultiQP"
TOL = 5.0e-10
MUTATION_MIN_DELTA = 1.0e-8


def safe_suffix(value):
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value)


def workdir(rootdir, mode, variant):
    suffix = "_projbf_fd_nonidentity"
    if variant == "reverse":
        suffix += "_reverse_projection"
    if mode == "ap":
        suffix += "_ap"
    if variant == "mutation":
        suffix += "_all_signs_positive"
    suffix += safe_suffix(os.environ.get("MVMC_BF_TEST_WORK_SUFFIX", ""))
    return os.path.join(rootdir, "work", MODEL + suffix)


def run_variant(rootdir, mode, variant):
    command = [
        sys.executable,
        os.path.join(rootdir, "runtest_bf.py"),
        MODEL,
        "--compare-proj-bf-finite-diff",
        "--use-nonidentity-init",
        "--check-multiqp-full-rebuild",
    ]
    if variant == "reverse":
        command.append("--reverse-projection-order")
    if mode == "ap":
        command.append("--use-ap-projection")
    if variant == "mutation":
        command.append("--mutate-ap-signs-positive")
    proc = subprocess.run(
        command,
        cwd=rootdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=os.environ.copy(),
    )
    if proc.returncode != 0:
        print(proc.stdout)
        raise RuntimeError("{} {} variant failed".format(mode, variant))
    path = os.path.join(
        workdir(rootdir, mode, variant), "bf_projbf_fd_dump.dat")
    if not os.path.exists(path):
        raise RuntimeError("finite-difference dump was not written: {}".format(path))
    return read_blocks(path)


def read_blocks(path):
    blocks = []
    current = {}
    with open(path) as source:
        for line in source:
            fields = line.split()
            if not fields:
                if current:
                    blocks.append(current)
                    current = {}
                continue
            current[fields[0]] = fields[1:]
    if current:
        blocks.append(current)
    if not blocks:
        raise RuntimeError("empty finite-difference dump: {}".format(path))
    return blocks


def integer(block, key):
    if key not in block or len(block[key]) != 1:
        raise RuntimeError("missing integer field {}".format(key))
    return int(block[key][0])


def floats(block, key, count):
    if key not in block or len(block[key]) != count:
        raise RuntimeError("invalid numeric field {}".format(key))
    values = [float(value) for value in block[key]]
    if not all(math.isfinite(value) for value in values):
        raise RuntimeError("non-finite numeric field {}".format(key))
    return values


def payload(block):
    nprojbf = integer(block, "nprojbf")
    nslater = integer(block, "nslater")
    values = floats(block, "ip_base", 2)
    values.extend(floats(block, "energy_base", 2))
    for idx in range(nprojbf):
        values.extend(floats(block, "projbf_derivative_{}".format(idx), 4))
    for idx in range(nslater):
        values.extend(floats(block, "orbital_derivative_{}".format(idx), 4))
    return values


def assert_fixture(blocks, mode, mutated=False):
    expected_ap = 1 if mode == "ap" else 0
    for block in blocks:
        if integer(block, "ap_flag") != expected_ap:
            raise RuntimeError("APFlag fixture assertion failed")
        if integer(block, "nmptrans") != 2:
            raise RuntimeError("NMPTrans fixture assertion failed")
        if integer(block, "nonzero_theta_count") <= 0:
            raise RuntimeError("Theta fixture assertion failed")
        negative = integer(block, "negative_qptrans_sign_count")
        if mode == "ap" and not mutated and negative <= 0:
            raise RuntimeError("AP sign census fixture assertion failed")
        if (mode == "pbc" or mutated) and negative != 0:
            raise RuntimeError("positive-sign fixture assertion failed")


def max_payload_difference(left, right):
    if len(left) != len(right):
        raise RuntimeError("sample count changed between projection variants")
    maximum = 0.0
    for left_block, right_block in zip(left, right):
        if integer(left_block, "sample") != integer(right_block, "sample"):
            raise RuntimeError("sample ordering changed between projection variants")
        left_values = payload(left_block)
        right_values = payload(right_block)
        if len(left_values) != len(right_values):
            raise RuntimeError("payload shape changed between projection variants")
        for left_value, right_value in zip(left_values, right_values):
            maximum = max(maximum, abs(left_value - right_value))
    return maximum


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in ("pbc", "ap"):
        print("usage: {} <pbc|ap>".format(sys.argv[0]))
        return 2
    mode = sys.argv[1]
    rootdir = os.getcwd()
    base = run_variant(rootdir, mode, "base")
    reverse = run_variant(rootdir, mode, "reverse")
    assert_fixture(base, mode)
    assert_fixture(reverse, mode)
    order_delta = max_payload_difference(base, reverse)
    if not math.isfinite(order_delta) or order_delta > TOL:
        print("ERROR: projection-order invariance mismatch: {:.3e}".format(order_delta))
        return 1

    if mode == "ap":
        mutation = run_variant(rootdir, mode, "mutation")
        assert_fixture(mutation, mode, mutated=True)
        mutation_delta = max_payload_difference(base, mutation)
        if (not math.isfinite(mutation_delta) or
                mutation_delta <= MUTATION_MIN_DELTA):
            print("ERROR: AP all-positive-sign mutation was not detected: "
                  "max_abs_delta={:.3e}".format(mutation_delta))
            return 1
        print("AP projection variants: order_delta={:.3e} "
              "mutation_delta={:.3e}".format(order_delta, mutation_delta))
    else:
        print("PBC projection variants: order_delta={:.3e}".format(order_delta))
    return 0


if __name__ == "__main__":
    sys.exit(main())
