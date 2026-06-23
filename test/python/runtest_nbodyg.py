from __future__ import print_function

import math
import os
import shutil
import subprocess
import sys


TOL = 1.0e-8


def copy_defs(src, dst):
    for fn in os.listdir(src):
        if fn.endswith(".def"):
            shutil.copy(os.path.join(src, fn), os.path.join(dst, fn))


def parse_complex_rows(path, nidx):
    rows = {}
    with open(path) as f:
        for line in f:
            cols = line.split()
            if not cols:
                continue
            if len(cols) != nidx + 2:
                raise ValueError("{}: expected {} columns, got {}".format(path, nidx + 2, len(cols)))
            key = tuple(int(cols[i]) for i in range(nidx))
            rows[key] = complex(float(cols[nidx]), float(cols[nidx + 1]))
    return rows


def parse_nbodyg_def(path):
    terms = []
    with open(path) as f:
        lines = f.readlines()[5:]
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        cols = stripped.split()
        n = int(cols[0])
        expected = 1 + 4 * n
        if len(cols) != expected:
            raise ValueError("{}: expected {} fields, got {}".format(path, expected, len(cols)))
        terms.append((n, tuple(int(c) for c in cols[1:])))
    return terms


def parse_nbodyg_output(path):
    rows = []
    with open(path) as f:
        for line in f:
            cols = line.split()
            if not cols:
                continue
            n = int(cols[0])
            expected = 1 + 4 * n + 2
            if len(cols) != expected:
                raise ValueError("{}: expected {} fields, got {}".format(path, expected, len(cols)))
            factors = tuple(int(c) for c in cols[1:1 + 4 * n])
            value = complex(float(cols[-2]), float(cols[-1]))
            rows.append((n, factors, value))
    return rows


def assert_close(label, actual, expected):
    diff = abs(actual - expected)
    if diff > TOL:
        raise AssertionError("{} mismatch: actual={} expected={} diff={}".format(label, actual, expected, diff))


def assert_finite(label, value):
    if not (math.isfinite(value.real) and math.isfinite(value.imag)):
        raise AssertionError("{} is not finite: {}".format(label, value))


def reduce_number_operator_factors(n, factors):
    reduced = []
    seen = set()
    for k in range(n):
        factor = factors[4 * k:4 * k + 4]
        if factor[0] != factor[2] or factor[1] != factor[3]:
            return None
        if factor not in seen:
            reduced.extend(factor)
            seen.add(factor)
    return tuple(reduced)


def main():
    if len(sys.argv) != 2:
        print("usage: {} <model name>".format(sys.argv[0]))
        return -1

    model = sys.argv[1]
    rootdir = os.getcwd()
    refdir = os.path.join(rootdir, "data", model)
    workdir = os.path.join(rootdir, "work", model)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    os.chdir(workdir)

    copy_defs(refdir, workdir)

    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    cmd = [bin_to_test, "-e", "namelist.def"]
    if os.path.exists("initial.def"):
        cmd.append("initial.def")
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    if mpi_procs:
        cmd = ["mpirun", "-np", mpi_procs] + cmd
    result = subprocess.call(cmd)
    if result != 0:
        return result

    terms_def = parse_nbodyg_def("nbodyg.def")
    nbodyg_output = os.path.join("output", "zvo_NBodyG_001.dat")
    if not terms_def:
        if os.path.exists(nbodyg_output):
            raise AssertionError("NBodyG output exists for NNBodyG=0")
        return 0

    rows = parse_nbodyg_output(nbodyg_output)
    if len(rows) != len(terms_def):
        raise AssertionError("NBodyG row count mismatch: output={} def={}".format(len(rows), len(terms_def)))

    one = parse_complex_rows(os.path.join("output", "zvo_cisajs_001.dat"), 4)
    two = parse_complex_rows(os.path.join("output", "zvo_cisajscktalt_001.dat"), 8)

    matched_one = 0
    matched_two = 0
    n_ge_3 = 0
    matched_reduced_n_ge_3 = 0
    comparable_reduced_n_ge_3 = 0
    for idx, (row, term_def) in enumerate(zip(rows, terms_def)):
        n, factors, value = row
        if (n, factors) != term_def:
            raise AssertionError("NBodyG term {} does not preserve input order".format(idx))
        assert_finite("NBodyG row {}".format(idx), value)
        if n == 1 and factors in one:
            assert_close("NBodyG N=1 row {}".format(idx), value, one[factors])
            matched_one += 1
        elif n == 2 and factors in two:
            assert_close("NBodyG N=2 row {}".format(idx), value, two[factors])
            matched_two += 1
        elif n >= 3:
            n_ge_3 += 1
            reduced = reduce_number_operator_factors(n, factors)
            if reduced is not None:
                if len(reduced) == 4 and reduced in one:
                    comparable_reduced_n_ge_3 += 1
                    assert_close("NBodyG reduced N>=3 row {}".format(idx), value, one[reduced])
                    matched_reduced_n_ge_3 += 1
                elif len(reduced) == 8 and reduced in two:
                    comparable_reduced_n_ge_3 += 1
                    assert_close("NBodyG reduced N>=3 row {}".format(idx), value, two[reduced])
                    matched_reduced_n_ge_3 += 1

    if matched_one == 0:
        raise AssertionError("no NBodyG N=1 row matched OneBodyG output")
    if matched_two == 0:
        raise AssertionError("no NBodyG N=2 row matched TwoBodyG output")
    if n_ge_3 == 0:
        raise AssertionError("no NBodyG N>=3 row was exercised")
    if model == "NBodyG_Compare" and matched_reduced_n_ge_3 == 0:
        raise AssertionError("no NBodyG N>=3 reduced row matched existing outputs")
    if matched_reduced_n_ge_3 != comparable_reduced_n_ge_3:
        raise AssertionError("not all comparable NBodyG N>=3 reduced rows were checked")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print("ERROR: {}".format(e))
        sys.exit(-1)
