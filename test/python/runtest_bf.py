from __future__ import print_function

import math
import os
import shutil
import subprocess
import sys

from backflow_def_helper import build_chain_nn_backflow, write_chain_nn_backflow


def parse_nsite(modpara_path):
    with open(modpara_path) as fp:
        for line in fp:
            cols = line.split()
            if len(cols) >= 2 and cols[0] == "Nsite":
                return int(cols[1])
    raise RuntimeError("Nsite was not found in {}".format(modpara_path))


def parse_norbitalidx(orbitalidx_path):
    with open(orbitalidx_path) as fp:
        for line in fp:
            cols = line.split()
            if len(cols) >= 2 and cols[0] == "NOrbitalIdx":
                return int(cols[1])
    raise RuntimeError("NOrbitalIdx was not found in {}".format(orbitalidx_path))


def read_namelist_keywords(namelist_path):
    keywords = set()
    with open(namelist_path) as fp:
        for line in fp:
            cols = line.split()
            if len(cols) >= 2:
                keywords.add(cols[0])
    return keywords


def assert_nonidentity_init_layout(workdir):
    keywords = read_namelist_keywords(os.path.join(workdir, "namelist.def"))
    nproj_keywords = set([
        "Gutzwiller",
        "Jastrow",
        "SpinJastrow",
        "DH2",
        "DH4",
    ])
    rbm_keywords = set([
        "ChargeRBM_HiddenLayer",
        "ChargeRBM_PhysLayer",
        "ChargeRBM_PhysHidden",
        "SpinRBM_HiddenLayer",
        "SpinRBM_PhysLayer",
        "SpinRBM_PhysHidden",
        "GeneralRBM_HiddenLayer",
        "GeneralRBM_PhysLayer",
        "GeneralRBM_PhysHidden",
    ])
    unsupported = (keywords & nproj_keywords) | (keywords & rbm_keywords)
    if "OptTrans" in keywords:
        unsupported.add("OptTrans")
    if unsupported:
        raise RuntimeError(
            "non-identity init writer assumes NProj=0, NRBM=0, and NOptTrans=0; unsupported entries: {}".format(
                ", ".join(sorted(unsupported))))


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


def max_abs_row_diff(left, right):
    if len(left) != len(right):
        raise RuntimeError("row length mismatch: left={} right={}".format(len(left), len(right)))
    diffs = [abs(x - y) for x, y in zip(left, right)]
    return max(diffs) if diffs else 0.0


def max_scaled_row_diff(left, right, abs_tol, rel_tol):
    if len(left) != len(right):
        raise RuntimeError("row length mismatch: left={} right={}".format(len(left), len(right)))
    scaled = []
    for x, y in zip(left, right):
        scale = max(abs_tol, rel_tol * max(abs(x), abs(y), 1.0))
        scaled.append(abs(x - y) / scale)
    return max(scaled) if scaled else 0.0


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


def write_nonidentity_init_parameter(path, nprojbf, nslater):
    projbf_values = [
        0.93, 0.08, -0.05, 0.035, -0.025,
        0.015, -0.012, 0.010, -0.007, 0.005,
    ]
    if nprojbf < 2:
        raise RuntimeError("non-identity BackFlow test requires at least two ProjBF parameters")

    values = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for idx in range(nprojbf):
        real = projbf_values[idx] if idx < len(projbf_values) else 0.001 / float(idx + 1)
        values.extend([real, 0.0, 0.0])
    matrix_width = int(round(math.sqrt(nslater)))
    for idx in range(nslater):
        if matrix_width * matrix_width == nslater:
            row = idx // matrix_width
            col = idx % matrix_width
            if row == col:
                real = 1.0 + 0.1 * row
            else:
                real = 0.03 * float(row - col)
        else:
            real = math.sin(0.37 * (idx + 1))
        values.extend([real, 0.0, 0.0])

    with open(path, "w") as fp:
        fp.write(" ".join("{:.18e}".format(value) for value in values))
        fp.write("\n")
    return projbf_values[:nprojbf]


def write_orbital_opt_flags(path, nsite, nslater, opt_flag):
    with open(path) as fp:
        lines = fp.readlines()

    mapping_left = nsite * nsite
    opt_left = nslater
    updated = []
    for line in lines:
        cols = line.split()
        if mapping_left > 0 and len(cols) >= 4:
            mapping_left -= 1
            updated.append(line)
        elif mapping_left == 0 and opt_left > 0 and len(cols) >= 2 and cols[0].lstrip("-").isdigit():
            idx = int(cols[0])
            updated.append("{:5d} {:6d}\n".format(idx, opt_flag))
            opt_left -= 1
        else:
            updated.append(line)

    if mapping_left != 0 or opt_left != 0:
        raise RuntimeError("failed to rewrite orbital OptFlag rows in {}".format(path))

    with open(path, "w") as fp:
        fp.writelines(updated)


def run_vmc(rootdir, workdir, mpi_procs, dump_path=None, diff_dump_path=None, fd_dump_path=None,
            init_path=None, log_name="bf_test.log"):
    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    env = os.environ.copy()
    if dump_path is not None:
        env["MVMC_BF_IDENTITY_DUMP"] = dump_path
    if diff_dump_path is not None:
        env["MVMC_BF_DIFF_DUMP"] = diff_dump_path
    if fd_dump_path is not None:
        env["MVMC_BF_FD_DUMP"] = fd_dump_path

    cmd = [bin_to_test, "-e", "namelist.def"]
    if init_path is not None:
        cmd.append(init_path)
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


def compare_real_complex_nonidentity(rootdir, real_model, complex_model, mpi_procs, tol, rejected_outputs):
    if mpi_procs:
        print("ERROR: non-identity real/complex comparison is a single-rank smoke test.")
        return -1

    real_refdir = os.path.join(rootdir, "data", real_model)
    complex_refdir = os.path.join(rootdir, "data", complex_model)
    real_workdir = os.path.join(rootdir, "work", real_model + "_nonidentity_real")
    complex_workdir = os.path.join(rootdir, "work", complex_model + "_nonidentity_complex")

    for refdir, workdir in ((real_refdir, real_workdir), (complex_refdir, complex_workdir)):
        if os.path.exists(workdir):
            shutil.rmtree(workdir)
        os.makedirs(workdir)
        copy_def_files(refdir, workdir, include_backflow=True)
        assert_nonidentity_init_layout(workdir)

    real_nsite = parse_nsite(os.path.join(real_workdir, "modpara.def"))
    complex_nsite = parse_nsite(os.path.join(complex_workdir, "modpara.def"))
    if real_nsite != complex_nsite:
        print("ERROR: Nsite mismatch: real={} complex={}".format(real_nsite, complex_nsite))
        return -1

    real_definition = build_chain_nn_backflow(length=real_nsite, optimize=False)
    complex_definition = build_chain_nn_backflow(length=complex_nsite, optimize=False)
    if real_definition.n_proj_bf != complex_definition.n_proj_bf:
        print("ERROR: NProjBF mismatch: real={} complex={}".format(
            real_definition.n_proj_bf, complex_definition.n_proj_bf))
        return -1
    write_chain_nn_backflow(real_workdir, length=real_nsite, optimize=False)
    write_chain_nn_backflow(complex_workdir, length=complex_nsite, optimize=False)

    real_nslater = parse_norbitalidx(os.path.join(real_workdir, "orbitalidx.def"))
    complex_nslater = parse_norbitalidx(os.path.join(complex_workdir, "orbitalidx.def"))
    if real_nslater != complex_nslater:
        print("ERROR: NSlater mismatch: real={} complex={}".format(real_nslater, complex_nslater))
        return -1
    write_orbital_opt_flags(os.path.join(real_workdir, "orbitalidx.def"), real_nsite, real_nslater, 0)
    write_orbital_opt_flags(os.path.join(complex_workdir, "orbitalidx.def"), complex_nsite, complex_nslater, 0)

    init_name = "nonidentity_init.dat"
    real_projbf = write_nonidentity_init_parameter(
        os.path.join(real_workdir, init_name),
        real_definition.n_proj_bf,
        real_nslater,
    )
    complex_projbf = write_nonidentity_init_parameter(
        os.path.join(complex_workdir, init_name),
        complex_definition.n_proj_bf,
        complex_nslater,
    )
    if real_projbf != complex_projbf or real_projbf[0] == 1.0 or all(value == 0.0 for value in real_projbf[1:]):
        print("ERROR: non-identity ProjBF initialization is invalid.")
        return -1

    real_proc = run_vmc(rootdir, real_workdir, mpi_procs, init_path=init_name,
                        log_name="bf_test_nonidentity_real.log")
    complex_proc = run_vmc(rootdir, complex_workdir, mpi_procs, init_path=init_name,
                           log_name="bf_test_nonidentity_complex.log")
    for label, proc in (("real", real_proc), ("complex", complex_proc)):
        if proc.returncode != 0:
            print("ERROR: {} BackFlow run failed.".format(label))
            print(proc.stdout)
            return proc.returncode
        rejected = contains_rejected_output(proc.stdout, rejected_outputs)
        if rejected is not None:
            print("ERROR: rejected {} BackFlow output substring found: {}".format(label, rejected))
            print("---- output begin ----")
            print(proc.stdout)
            print("---- output end ----")
            return -1

    real_row = read_first_float_row(os.path.join(real_workdir, "output", "zvo_out_001.dat"))
    complex_row = read_first_float_row(os.path.join(complex_workdir, "output", "zvo_out_001.dat"))
    if not all(math.isfinite(value) for value in real_row + complex_row):
        print("ERROR: non-identity real/complex output contains non-finite values.")
        return -1
    max_diff = max_abs_row_diff(real_row, complex_row)
    rel_tol = 1.0e-12
    max_scaled_diff = max_scaled_row_diff(real_row, complex_row, tol, rel_tol)
    if not math.isfinite(max_diff) or not math.isfinite(max_scaled_diff) or max_scaled_diff > 1.0:
        print("ERROR: non-identity real/complex mismatch: max_abs_diff={:.3e} max_scaled_diff={:.3e}".format(
            max_diff, max_scaled_diff))
        print("real    zvo_out_001.dat: {}".format(" ".join("{:.18e}".format(x) for x in real_row)))
        print("complex zvo_out_001.dat: {}".format(" ".join("{:.18e}".format(x) for x in complex_row)))
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
        print("usage: {} <model name> [--expect-error <substring>] [--expect-nqp-full <n>] [--compare-no-bf-energy] [--compare-no-bf-gradient] [--compare-proj-bf-finite-diff] [--compare-real-complex-nonidentity <complex model>] [--reject-output <substring>]".format(sys.argv[0]))
        return -1

    model = sys.argv[1]
    expected_error = None
    expected_nqp_full = None
    compare_energy = False
    compare_gradient = False
    compare_proj_bf_fd = False
    compare_real_complex_model = None
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
        elif sys.argv[argi] == "--compare-real-complex-nonidentity" and argi + 1 < len(sys.argv):
            compare_real_complex_model = sys.argv[argi + 1]
            argi += 2
        elif sys.argv[argi] == "--reject-output" and argi + 1 < len(sys.argv):
            rejected_outputs.append(sys.argv[argi + 1])
            argi += 2
        else:
            print("usage: {} <model name> [--expect-error <substring>] [--expect-nqp-full <n>] [--compare-no-bf-energy] [--compare-no-bf-gradient] [--compare-proj-bf-finite-diff] [--compare-real-complex-nonidentity <complex model>] [--reject-output <substring>]".format(sys.argv[0]))
            return -1
    rootdir = os.getcwd()
    refdir = os.path.join(rootdir, "data", model)
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    if compare_real_complex_model is not None:
        return compare_real_complex_nonidentity(
            rootdir,
            model,
            compare_real_complex_model,
            mpi_procs,
            1.0e-10,
            rejected_outputs,
        )

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
