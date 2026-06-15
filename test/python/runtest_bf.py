from __future__ import print_function

import math
import os
import shutil
import subprocess
import sys

from backflow_def_helper import write_chain_nn_backflow


def parse_nsite(modpara_path):
    with open(modpara_path) as fp:
        for line in fp:
            cols = line.split()
            if len(cols) >= 2 and cols[0] == "Nsite":
                return int(cols[1])
    raise RuntimeError("Nsite was not found in {}".format(modpara_path))


def read_key_value_file(path):
    values = {}
    with open(path) as fp:
        for line in fp:
            cols = line.split()
            if len(cols) == 2:
                values[cols[0]] = cols[1]
    return values


def first_row_is_finite(path):
    with open(path) as fp:
        for line in fp:
            cols = line.split()
            if not cols:
                continue
            return all(math.isfinite(float(x)) for x in cols)
    return False


def read_first_float_row(path):
    with open(path) as fp:
        for line in fp:
            cols = line.split()
            if cols:
                return [float(x) for x in cols]
    raise RuntimeError("no numeric rows in {}".format(path))


def read_key_value_blocks(path):
    blocks = []
    current = {}
    with open(path) as fp:
        for line in fp:
            cols = line.split()
            if not cols:
                if current:
                    blocks.append(current)
                    current = {}
                continue
            current[cols[0]] = cols[1:]
    if current:
        blocks.append(current)
    return blocks


def copy_def_files(refdir, workdir, include_backflow):
    for fn in os.listdir(refdir):
        if not fn.endswith(".def"):
            continue
        src_path = os.path.join(refdir, fn)
        dst_path = os.path.join(workdir, fn)
        if fn == "namelist.def" and not include_backflow:
            lines = []
            with open(src_path) as fp:
                for line in fp:
                    cols = line.split()
                    if cols and cols[0] in ("BF", "BFRange"):
                        continue
                    lines.append(line)
            with open(dst_path, "w") as fp:
                fp.writelines(lines)
        else:
            shutil.copy(src_path, dst_path)


def run_vmc(rootdir, workdir, mpi_procs, dump_path=None, diff_dump_path=None, fd_dump_path=None,
            log_name="bf_test.log"):
    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    env = os.environ.copy()
    if dump_path is not None:
        env["MVMC_BF_IDENTITY_DUMP"] = dump_path
    if diff_dump_path is not None:
        env["MVMC_BF_DIFF_DUMP"] = diff_dump_path
    if fd_dump_path is not None:
        env["MVMC_BF_FD_DUMP"] = fd_dump_path

    cmd = [bin_to_test, "-e", "namelist.def"]
    if mpi_procs:
        cmd = ["mpirun", "-np", mpi_procs] + cmd
    proc = subprocess.run(
        cmd,
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env,
    )
    with open(os.path.join(workdir, log_name), "w") as fp:
        fp.write(proc.stdout)
    return proc


def contains_rejected_output(output, rejected_outputs):
    for substring in rejected_outputs:
        if substring in output:
            return substring
    return None


def compare_no_bf_energy(rootdir, refdir, bf_workdir, mpi_procs, tol, rejected_outputs):
    no_bf_workdir = bf_workdir + "_nobf_energy"

    if os.path.exists(no_bf_workdir):
        shutil.rmtree(no_bf_workdir)
    os.makedirs(no_bf_workdir)
    copy_def_files(refdir, no_bf_workdir, include_backflow=False)

    proc = run_vmc(rootdir, no_bf_workdir, mpi_procs, log_name="bf_test_nobf.log")
    if proc.returncode != 0:
        print(proc.stdout)
        return proc.returncode
    rejected = contains_rejected_output(proc.stdout, rejected_outputs)
    if rejected is not None:
        print("ERROR: rejected no-BF output substring found: {}".format(rejected))
        print("---- output begin ----")
        print(proc.stdout)
        print("---- output end ----")
        return -1

    bf_out_path = os.path.join(bf_workdir, "output", "zvo_out_001.dat")
    no_bf_out_path = os.path.join(no_bf_workdir, "output", "zvo_out_001.dat")
    bf_row = read_first_float_row(bf_out_path)
    no_bf_row = read_first_float_row(no_bf_out_path)
    if len(bf_row) != len(no_bf_row):
        print("ERROR: zvo_out column mismatch: bf={} no_bf={}".format(len(bf_row), len(no_bf_row)))
        return -1

    diffs = [abs(x - y) for x, y in zip(bf_row, no_bf_row)]
    max_diff = max(diffs) if diffs else 0.0
    if not math.isfinite(max_diff) or max_diff > tol:
        print("ERROR: identity BackFlow energy mismatch: max_abs_diff={:.3e}".format(max_diff))
        print("BF    zvo_out_001.dat: {}".format(" ".join("{:.18e}".format(x) for x in bf_row)))
        print("no-BF zvo_out_001.dat: {}".format(" ".join("{:.18e}".format(x) for x in no_bf_row)))
        return -1
    return 0


def check_no_bf_gradient_dump(path, tol):
    if not os.path.exists(path):
        print("ERROR: BackFlow diff dump was not written.")
        return -1

    blocks = read_key_value_blocks(path)
    if not blocks:
        print("ERROR: BackFlow diff dump is empty.")
        return -1

    maxima = {
        "max_abs_slater_diff_o": 0.0,
        "max_abs_main_order_slater_diff_o": 0.0,
        "abs_energy_diff": 0.0,
    }
    for block in blocks:
        info_no_bf = int(block.get("info_no_bf", ["-1"])[0])
        info_bf = int(block.get("info_bf", ["-1"])[0])
        nan_count = int(block.get("nan_count", ["-1"])[0])
        if info_no_bf != 0 or info_bf != 0:
            print("ERROR: diff dump matrix calculation failed: no_bf={} bf={}".format(info_no_bf, info_bf))
            return -1
        if nan_count != 0:
            print("ERROR: diff dump contains non-finite values: nan_count={}".format(nan_count))
            return -1
        for key in maxima:
            value = float(block.get(key, ["nan"])[0])
            if not math.isfinite(value):
                print("ERROR: diff dump contains non-finite {}.".format(key))
                return -1
            if value > maxima[key]:
                maxima[key] = value

    for key, value in maxima.items():
        if value > tol:
            print("ERROR: identity BackFlow gradient mismatch: {}={:.3e}".format(key, value))
            return -1
    return 0


def check_proj_bf_finite_diff_dump(path, tol):
    if not os.path.exists(path):
        print("ERROR: BackFlow finite-difference dump was not written.")
        return -1

    blocks = read_key_value_blocks(path)
    if not blocks:
        print("ERROR: BackFlow finite-difference dump is empty.")
        return -1

    total_nonzero = 0
    maxima = {
        "max_abs_projbf_fd_real": 0.0,
        "max_abs_projbf_fd_imag": 0.0,
        "max_abs_projbf_fd_diff": 0.0,
        "max_abs_fd_value": 0.0,
    }
    for block in blocks:
        info_base = int(block.get("info_base", ["-1"])[0])
        fd_fail_count = int(block.get("fd_fail_count", ["-1"])[0])
        nan_count = int(block.get("nan_count", ["-1"])[0])
        nprojbf = int(block.get("nprojbf", ["0"])[0])
        total_nonzero += int(block.get("nonzero_fd_count", ["0"])[0])
        if info_base != 0:
            print("ERROR: finite-difference base matrix calculation failed: info_base={}".format(info_base))
            return -1
        if fd_fail_count != 0:
            print("ERROR: finite-difference perturbation failed: fd_fail_count={}".format(fd_fail_count))
            return -1
        if nan_count != 0:
            print("ERROR: finite-difference dump contains non-finite values: nan_count={}".format(nan_count))
            return -1
        if nprojbf <= 0:
            print("ERROR: finite-difference dump has invalid nprojbf={}".format(nprojbf))
            return -1
        for key in maxima:
            value = float(block.get(key, ["nan"])[0])
            if not math.isfinite(value):
                print("ERROR: finite-difference dump contains non-finite {}.".format(key))
                return -1
            if value > maxima[key]:
                maxima[key] = value

    if total_nonzero == 0 or maxima["max_abs_fd_value"] <= 1.0e-12:
        print("ERROR: finite-difference check was vacuous: nonzero_fd_count={}".format(total_nonzero))
        return -1
    if maxima["max_abs_projbf_fd_diff"] > tol:
        print("ERROR: ProjBF finite-difference mismatch: max_abs_diff={:.3e}".format(
            maxima["max_abs_projbf_fd_diff"]))
        return -1
    return 0


def main():
    if len(sys.argv) < 2:
        print("usage: {} <model name> [--expect-error <substring>] [--expect-nqp-full <n>] [--compare-no-bf-energy] [--compare-no-bf-gradient] [--compare-proj-bf-finite-diff] [--reject-output <substring>]".format(sys.argv[0]))
        return -1

    model = sys.argv[1]
    expected_error = None
    expected_nqp_full = None
    compare_energy = False
    compare_gradient = False
    compare_proj_bf_fd = False
    rejected_outputs = []
    argi = 2
    while argi < len(sys.argv):
        if sys.argv[argi] == "--expect-error" and argi + 1 < len(sys.argv):
            expected_error = sys.argv[argi + 1]
            argi += 2
        elif sys.argv[argi] == "--expect-nqp-full" and argi + 1 < len(sys.argv):
            expected_nqp_full = int(sys.argv[argi + 1])
            argi += 2
        elif sys.argv[argi] == "--compare-no-bf-energy":
            compare_energy = True
            argi += 1
        elif sys.argv[argi] == "--compare-no-bf-gradient":
            compare_gradient = True
            argi += 1
        elif sys.argv[argi] == "--compare-proj-bf-finite-diff":
            compare_proj_bf_fd = True
            argi += 1
        elif sys.argv[argi] == "--reject-output" and argi + 1 < len(sys.argv):
            rejected_outputs.append(sys.argv[argi + 1])
            argi += 2
        else:
            print("usage: {} <model name> [--expect-error <substring>] [--expect-nqp-full <n>] [--compare-no-bf-energy] [--compare-no-bf-gradient] [--compare-proj-bf-finite-diff] [--reject-output <substring>]".format(sys.argv[0]))
            return -1
    rootdir = os.getcwd()
    refdir = os.path.join(rootdir, "data", model)
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    work_suffix = ""
    if compare_energy:
        work_suffix += "_energy"
    if compare_gradient:
        work_suffix += "_gradient"
    if compare_proj_bf_fd:
        work_suffix += "_projbf_fd"
    workdir = os.path.join(rootdir, "work", model + work_suffix)
    if mpi_procs:
        workdir += "_mpi{}".format(mpi_procs)

    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)

    copy_def_files(refdir, workdir, include_backflow=True)

    nsite = parse_nsite(os.path.join(workdir, "modpara.def"))
    write_chain_nn_backflow(workdir, length=nsite, optimize=compare_proj_bf_fd)

    dump_path = os.path.join(workdir, "bf_identity_dump.dat")
    diff_dump_path = os.path.join(workdir, "bf_diff_dump.dat") if compare_gradient else None
    fd_dump_path = os.path.join(workdir, "bf_projbf_fd_dump.dat") if compare_proj_bf_fd else None
    proc = run_vmc(rootdir, workdir, mpi_procs, dump_path=dump_path, diff_dump_path=diff_dump_path,
                   fd_dump_path=fd_dump_path)

    if expected_error is not None:
        if expected_error not in proc.stdout:
            print("ERROR: expected error substring not found: {}".format(expected_error))
            print("---- output begin ----")
            print(proc.stdout)
            print("---- output end ----")
            return -1
        if proc.returncode == 0:
            print("ERROR: vmc.out unexpectedly succeeded for invalid BackFlow input.")
            return -1
        return 0

    if proc.returncode != 0:
        print(proc.stdout)
        return proc.returncode
    rejected = contains_rejected_output(proc.stdout, rejected_outputs)
    if rejected is not None:
        print("ERROR: rejected BackFlow output substring found: {}".format(rejected))
        print("---- output begin ----")
        print(proc.stdout)
        print("---- output end ----")
        return -1
    if not os.path.exists(dump_path):
        print("ERROR: BackFlow identity dump was not written.")
        return -1

    values = read_key_value_file(dump_path)
    info_no_bf = int(values.get("info_no_bf", "-1"))
    info_bf = int(values.get("info_bf", "-1"))
    nqp_full = int(values.get("nqp_full", "-1"))
    nan_count = int(values.get("nan_count", "-1"))
    max_slater = float(values.get("max_abs_slater_diff", "nan"))
    max_pf = float(values.get("max_abs_pf_diff", "nan"))
    tol = 1.0e-10

    if expected_nqp_full is not None and nqp_full != expected_nqp_full:
        print("ERROR: unexpected nqp_full: got {} expected {}".format(nqp_full, expected_nqp_full))
        return -1
    if info_no_bf != 0 or info_bf != 0:
        print("ERROR: identity dump matrix calculation failed: no_bf={} bf={}".format(info_no_bf, info_bf))
        return -1
    if nan_count != 0:
        print("ERROR: identity dump contains non-finite values: nan_count={}".format(nan_count))
        return -1
    if not math.isfinite(max_slater) or max_slater > tol:
        print("ERROR: SlaterElm identity mismatch: max_abs_slater_diff={}".format(max_slater))
        return -1
    if not math.isfinite(max_pf) or max_pf > tol:
        print("ERROR: PfM identity mismatch: max_abs_pf_diff={}".format(max_pf))
        return -1
    if not first_row_is_finite(os.path.join(workdir, "output", "zvo_out_001.dat")):
        print("ERROR: zvo_out_001.dat is missing or non-finite.")
        return -1
    if compare_energy:
        result = compare_no_bf_energy(rootdir, refdir, workdir, mpi_procs, tol, rejected_outputs)
        if result != 0:
            return result
    if compare_gradient:
        result = check_no_bf_gradient_dump(diff_dump_path, tol)
        if result != 0:
            return result
    if compare_proj_bf_fd:
        result = check_proj_bf_finite_diff_dump(fd_dump_path, 1.0e-6)
        if result != 0:
            return result

    return 0


if __name__ == "__main__":
    sys.exit(main())
