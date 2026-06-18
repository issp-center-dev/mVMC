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


def update_modpara(path, values):
    found = set()
    lines = []
    with open(path) as fp:
        for line in fp:
            cols = line.split()
            if cols and cols[0] in values:
                lines.append("{:<15} {}\n".format(cols[0], values[cols[0]]))
                found.add(cols[0])
            else:
                lines.append(line)

    for key, value in values.items():
        if key not in found:
            lines.append("{:<15} {}\n".format(key, value))

    with open(path, "w") as fp:
        fp.writelines(lines)


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


def read_float_rows(path):
    rows = []
    with open(path) as fp:
        for line in fp:
            cols = line.split()
            if cols:
                rows.append([float(x) for x in cols])
    return rows


def read_child_opt_values(path):
    values = []
    with open(path) as fp:
        for line in fp:
            cols = line.split()
            if len(cols) != 3:
                continue
            try:
                idx = int(cols[0])
                real = float(cols[1])
                imag = float(cols[2])
            except ValueError:
                continue
            values.append((idx, real, imag))
    return values


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


def write_minimal_twobodyg(workdir, nsite):
    if nsite < 4:
        raise RuntimeError("minimal TwoBodyG BackFlow smoke requires Nsite >= 4")
    rows = [
        (0, 0, 0, 0, 0, 1, 0, 1),
        (0, 0, 0, 0, 1, 1, 1, 1),
        (0, 0, 1, 0, 2, 1, 3, 1),
        (0, 0, 1, 0, 2, 0, 3, 0),
        (1, 1, 0, 1, 3, 0, 2, 0),
        (2, 0, 3, 0, 1, 1, 0, 1),
    ]
    with open(os.path.join(workdir, "greentwo.def"), "w") as fp:
        fp.write("=============================================\n")
        fp.write("NCisAjsCktAltDC         {}\n".format(len(rows)))
        fp.write("=============================================\n")
        fp.write("======== Green functions for BF smoke =======\n")
        fp.write("=============================================\n")
        for row in rows:
            fp.write(" ".join("{:5d}".format(value) for value in row))
            fp.write("\n")

    namelist_path = os.path.join(workdir, "namelist.def")
    with open(namelist_path) as fp:
        lines = fp.readlines()
    if not any(line.split() and line.split()[0] == "TwoBodyG" for line in lines):
        lines.append("        TwoBodyG  greentwo.def\n")
        with open(namelist_path, "w") as fp:
            fp.writelines(lines)


def write_minimal_twobodygex(workdir, nsite):
    if nsite < 4:
        raise RuntimeError("minimal TwoBodyGEx BackFlow smoke requires Nsite >= 4")
    rows = [
        (0, 0, 1, 0, 2, 1, 3, 1),
        (0, 0, 1, 0, 2, 0, 3, 0),
        (1, 1, 0, 1, 3, 0, 2, 0),
        (2, 0, 3, 0, 1, 1, 0, 1),
    ]
    with open(os.path.join(workdir, "greentwoex.def"), "w") as fp:
        fp.write("=============================================\n")
        fp.write("NCisAjsCktAlt         {}\n".format(len(rows)))
        fp.write("=============================================\n")
        fp.write("======== Green functions for BF smoke =======\n")
        fp.write("=============================================\n")
        for row in rows:
            fp.write(" ".join("{:5d}".format(value) for value in row))
            fp.write("\n")

    namelist_path = os.path.join(workdir, "namelist.def")
    with open(namelist_path) as fp:
        lines = fp.readlines()
    if not any(line.split() and line.split()[0] == "TwoBodyGEx" for line in lines):
        lines.append("      TwoBodyGEx  greentwoex.def\n")
        with open(namelist_path, "w") as fp:
            fp.writelines(lines)


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
            green2_dump_path=None, init_path=None, log_name="bf_test.log"):
    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    env = os.environ.copy()
    if dump_path is not None:
        env["MVMC_BF_IDENTITY_DUMP"] = dump_path
    if diff_dump_path is not None:
        env["MVMC_BF_DIFF_DUMP"] = diff_dump_path
    if fd_dump_path is not None:
        env["MVMC_BF_FD_DUMP"] = fd_dump_path
    if green2_dump_path is not None:
        env["MVMC_BF_GREEN2_DUMP"] = green2_dump_path

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


def compare_no_bf_twobodyg(rootdir, refdir, bf_workdir, mpi_procs, tol, rejected_outputs):
    no_bf_workdir = bf_workdir + "_nobf_twobodyg"

    if os.path.exists(no_bf_workdir):
        shutil.rmtree(no_bf_workdir)
    os.makedirs(no_bf_workdir)
    copy_def_files(refdir, no_bf_workdir, include_backflow=False)
    nsite = parse_nsite(os.path.join(no_bf_workdir, "modpara.def"))
    write_minimal_twobodyg(no_bf_workdir, nsite)

    proc = run_vmc(rootdir, no_bf_workdir, mpi_procs, log_name="bf_test_nobf_twobodyg.log")
    if proc.returncode != 0:
        print(proc.stdout)
        return proc.returncode
    rejected = contains_rejected_output(proc.stdout, rejected_outputs)
    if rejected is not None:
        print("ERROR: rejected no-BF TwoBodyG output substring found: {}".format(rejected))
        print("---- output begin ----")
        print(proc.stdout)
        print("---- output end ----")
        return -1

    bf_path = os.path.join(bf_workdir, "output", "zvo_cisajscktalt_001.dat")
    no_bf_path = os.path.join(no_bf_workdir, "output", "zvo_cisajscktalt_001.dat")
    bf_rows = read_float_rows(bf_path)
    no_bf_rows = read_float_rows(no_bf_path)
    if len(bf_rows) != len(no_bf_rows):
        print("ERROR: TwoBodyG row count mismatch: bf={} no_bf={}".format(len(bf_rows), len(no_bf_rows)))
        return -1
    if not bf_rows:
        print("ERROR: TwoBodyG output is empty.")
        return -1

    max_diff = 0.0
    max_no_bf_observable = 0.0
    for row_idx, (bf_row, no_bf_row) in enumerate(zip(bf_rows, no_bf_rows)):
        if len(bf_row) != len(no_bf_row):
            print("ERROR: TwoBodyG column mismatch at row {}: bf={} no_bf={}".format(
                row_idx, len(bf_row), len(no_bf_row)))
            return -1
        max_diff = max(max_diff, max_abs_row_diff(bf_row, no_bf_row))
        if len(no_bf_row) >= 2:
            max_no_bf_observable = max(max_no_bf_observable, abs(no_bf_row[-2]), abs(no_bf_row[-1]))
    if max_no_bf_observable <= 1.0e-12:
        print("ERROR: TwoBodyG comparison is vacuous: no-BF observables are all zero.")
        return -1
    if not math.isfinite(max_diff) or max_diff > tol:
        print("ERROR: identity BackFlow TwoBodyG mismatch: max_abs_diff={:.3e}".format(max_diff))
        return -1
    return 0


def compare_no_bf_twobodygex(rootdir, refdir, bf_workdir, mpi_procs, tol, rejected_outputs):
    no_bf_workdir = bf_workdir + "_nobf_twobodygex"

    if os.path.exists(no_bf_workdir):
        shutil.rmtree(no_bf_workdir)
    os.makedirs(no_bf_workdir)
    copy_def_files(refdir, no_bf_workdir, include_backflow=False)
    nsite = parse_nsite(os.path.join(no_bf_workdir, "modpara.def"))
    write_minimal_twobodygex(no_bf_workdir, nsite)

    proc = run_vmc(rootdir, no_bf_workdir, mpi_procs, log_name="bf_test_nobf_twobodygex.log")
    if proc.returncode != 0:
        print(proc.stdout)
        return proc.returncode
    rejected = contains_rejected_output(proc.stdout, rejected_outputs)
    if rejected is not None:
        print("ERROR: rejected no-BF TwoBodyGEx output substring found: {}".format(rejected))
        print("---- output begin ----")
        print(proc.stdout)
        print("---- output end ----")
        return -1

    bf_path = os.path.join(bf_workdir, "output", "zvo_cisajscktaltex_001.dat")
    no_bf_path = os.path.join(no_bf_workdir, "output", "zvo_cisajscktaltex_001.dat")
    bf_row = read_first_float_row(bf_path)
    no_bf_row = read_first_float_row(no_bf_path)
    max_observable = max(abs(value) for value in no_bf_row) if no_bf_row else 0.0
    if max_observable <= 1.0e-12:
        print("ERROR: TwoBodyGEx comparison is vacuous: no-BF observables are all zero.")
        return -1
    max_diff = max_abs_row_diff(bf_row, no_bf_row)
    if not math.isfinite(max_diff) or max_diff > tol:
        print("ERROR: identity BackFlow TwoBodyGEx mismatch: max_abs_diff={:.3e}".format(max_diff))
        print("BF    zvo_cisajscktaltex_001.dat: {}".format(" ".join("{:.18e}".format(x) for x in bf_row)))
        print("no-BF zvo_cisajscktaltex_001.dat: {}".format(" ".join("{:.18e}".format(x) for x in no_bf_row)))
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


def check_child_opt_block(all_values, child_values, start, label, tol):
    if len(child_values) == 0:
        print("ERROR: {} opt file contains no parameter rows.".format(label))
        return -1
    if start + 3 * len(child_values) > len(all_values):
        print("ERROR: {} block overruns zqp_opt.dat: start={} rows={} len={}".format(
            label, start, len(child_values), len(all_values)))
        return -1

    max_diff = 0.0
    for row, (idx, real, imag) in enumerate(child_values):
        if idx != row:
            print("ERROR: {} opt row index mismatch: row={} idx={}".format(label, row, idx))
            return -1
        all_real = all_values[start + 3 * row]
        all_imag = all_values[start + 3 * row + 1]
        max_diff = max(max_diff, abs(all_real - real), abs(all_imag - imag))
    if not math.isfinite(max_diff) or max_diff > tol:
        print("ERROR: {} opt block mismatch: max_abs_diff={:.3e}".format(label, max_diff))
        return -1
    return 0


def check_opt_output_restart(rootdir, model, mpi_procs, rejected_outputs):
    if mpi_procs:
        print("ERROR: BackFlow opt-output restart smoke is a single-rank test.")
        return -1

    refdir = os.path.join(rootdir, "data", model)
    workdir = os.path.join(rootdir, "work", model + "_opt_output_restart")
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    copy_def_files(refdir, workdir, include_backflow=True)
    assert_nonidentity_init_layout(workdir)

    update_modpara(
        os.path.join(workdir, "modpara.def"),
        {
            "NSROptItrStep": 2,
            "NSROptItrSmp": 2,
            "NVMCSample": 16,
            "NVMCWarmUp": 4,
            "NStore": 1,
        },
    )

    nsite = parse_nsite(os.path.join(workdir, "modpara.def"))
    definition = build_chain_nn_backflow(length=nsite, optimize=True)
    write_chain_nn_backflow(workdir, length=nsite, optimize=True)
    nslater = parse_norbitalidx(os.path.join(workdir, "orbitalidx.def"))

    proc = run_vmc(rootdir, workdir, mpi_procs, log_name="bf_test_opt_output.log")
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

    zqp_opt = os.path.join(workdir, "output", "zqp_opt.dat")
    zqp_bf_opt = os.path.join(workdir, "output", "zqp_bf_opt.dat")
    zqp_orbital_opt = os.path.join(workdir, "output", "zqp_orbital_opt.dat")
    for path in (zqp_opt, zqp_bf_opt, zqp_orbital_opt):
        if not os.path.exists(path):
            print("ERROR: expected opt output was not written: {}".format(path))
            return -1

    all_values = read_first_float_row(zqp_opt)
    expected_len = 6 + 3 * (definition.n_proj_bf + nslater)
    if len(all_values) != expected_len:
        print("ERROR: zqp_opt.dat length mismatch: got {} expected {}".format(
            len(all_values), expected_len))
        return -1
    if not all(math.isfinite(value) for value in all_values):
        print("ERROR: zqp_opt.dat contains non-finite values.")
        return -1

    bf_values = read_child_opt_values(zqp_bf_opt)
    if len(bf_values) != definition.n_proj_bf:
        print("ERROR: zqp_bf_opt.dat row count mismatch: got {} expected {}".format(
            len(bf_values), definition.n_proj_bf))
        return -1
    result = check_child_opt_block(all_values, bf_values, 6, "BackFlow", 1.0e-15)
    if result != 0:
        return result

    orbital_values = read_child_opt_values(zqp_orbital_opt)
    if len(orbital_values) != nslater:
        print("ERROR: zqp_orbital_opt.dat row count mismatch: got {} expected {}".format(
            len(orbital_values), nslater))
        return -1
    result = check_child_opt_block(
        all_values,
        orbital_values,
        6 + 3 * definition.n_proj_bf,
        "Slater",
        1.0e-15,
    )
    if result != 0:
        return result

    restart_proc = run_vmc(
        rootdir,
        workdir,
        mpi_procs,
        init_path=os.path.join("output", "zqp_opt.dat"),
        log_name="bf_test_opt_output_restart.log",
    )
    if restart_proc.returncode != 0:
        print(restart_proc.stdout)
        return restart_proc.returncode
    rejected = contains_rejected_output(restart_proc.stdout, rejected_outputs)
    if rejected is not None:
        print("ERROR: rejected BackFlow restart output substring found: {}".format(rejected))
        print("---- output begin ----")
        print(restart_proc.stdout)
        print("---- output end ----")
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


def check_bf_green2_bruteforce_dump(path, tol, expected_all_complex_flag):
    if not os.path.exists(path):
        print("ERROR: BackFlow GreenFunc2 brute-force dump was not written.")
        return -1

    blocks = read_key_value_blocks(path)
    if not blocks:
        print("ERROR: BackFlow GreenFunc2 brute-force dump is empty.")
        return -1

    total_compared = 0
    total_nonzero = 0
    max_diff = 0.0
    max_bruteforce = 0.0
    for block in blocks:
        if expected_all_complex_flag is not None:
            all_complex_flag = int(block.get("all_complex_flag", ["-1"])[0])
            if all_complex_flag != expected_all_complex_flag:
                print("ERROR: unexpected all_complex_flag in GreenFunc2 dump: got {} expected {}".format(
                    all_complex_flag, expected_all_complex_flag))
                return -1
        info_fail_count = int(block.get("info_fail_count", ["-1"])[0])
        nan_count = int(block.get("nan_count", ["-1"])[0])
        if info_fail_count != 0:
            print("ERROR: GreenFunc2 brute-force dump has info_fail_count={}".format(info_fail_count))
            return -1
        if nan_count != 0:
            print("ERROR: GreenFunc2 brute-force dump has nan_count={}".format(nan_count))
            return -1
        compared = int(block.get("compared_count", ["0"])[0])
        nonzero = int(block.get("nonzero_bruteforce_count", ["0"])[0])
        block_diff = float(block.get("max_abs_green2_diff", ["nan"])[0])
        block_bruteforce = float(block.get("max_abs_bruteforce", ["nan"])[0])
        if not math.isfinite(block_diff) or not math.isfinite(block_bruteforce):
            print("ERROR: GreenFunc2 brute-force dump contains non-finite maxima.")
            return -1
        total_compared += compared
        total_nonzero += nonzero
        max_diff = max(max_diff, block_diff)
        max_bruteforce = max(max_bruteforce, block_bruteforce)

    if total_compared == 0:
        print("ERROR: GreenFunc2 brute-force check was vacuous: compared_count=0.")
        return -1
    if total_nonzero == 0 or max_bruteforce <= 1.0e-12:
        print("ERROR: GreenFunc2 brute-force check was vacuous: nonzero_bruteforce_count={}".format(
            total_nonzero))
        return -1
    if max_diff > tol:
        print("ERROR: BackFlow GreenFunc2 brute-force mismatch: max_abs_diff={:.3e}".format(max_diff))
        return -1
    return 0


def main():
    if len(sys.argv) < 2:
        print("usage: {} <model name> [--expect-error <substring>] [--expect-nqp-full <n>] [--compare-no-bf-energy] [--compare-no-bf-twobodyg] [--compare-no-bf-twobodygex] [--compare-no-bf-gradient] [--compare-proj-bf-finite-diff] [--check-bf-green2-bruteforce] [--use-nonidentity-init] [--set-ncond <n>] [--expect-all-complex-flag <0|1>] [--compare-real-complex-nonidentity <complex model>] [--check-opt-output-restart] [--reject-output <substring>]".format(sys.argv[0]))
        return -1

    model = sys.argv[1]
    expected_error = None
    expected_nqp_full = None
    compare_energy = False
    compare_twobodyg = False
    compare_twobodygex = False
    compare_gradient = False
    compare_proj_bf_fd = False
    check_bf_green2_bruteforce = False
    use_nonidentity_init = False
    expected_all_complex_flag = None
    compare_real_complex_model = None
    check_opt_restart = False
    ncond_override = None
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
        elif sys.argv[argi] == "--compare-no-bf-twobodyg":
            compare_twobodyg = True
            argi += 1
        elif sys.argv[argi] == "--compare-no-bf-twobodygex":
            compare_twobodygex = True
            argi += 1
        elif sys.argv[argi] == "--compare-no-bf-gradient":
            compare_gradient = True
            argi += 1
        elif sys.argv[argi] == "--compare-proj-bf-finite-diff":
            compare_proj_bf_fd = True
            argi += 1
        elif sys.argv[argi] == "--check-bf-green2-bruteforce":
            check_bf_green2_bruteforce = True
            argi += 1
        elif sys.argv[argi] == "--use-nonidentity-init":
            use_nonidentity_init = True
            argi += 1
        elif sys.argv[argi] == "--set-ncond" and argi + 1 < len(sys.argv):
            ncond_override = str(int(sys.argv[argi + 1]))
            argi += 2
        elif sys.argv[argi] == "--expect-all-complex-flag" and argi + 1 < len(sys.argv):
            expected_all_complex_flag = int(sys.argv[argi + 1])
            argi += 2
        elif sys.argv[argi] == "--compare-real-complex-nonidentity" and argi + 1 < len(sys.argv):
            compare_real_complex_model = sys.argv[argi + 1]
            argi += 2
        elif sys.argv[argi] == "--check-opt-output-restart":
            check_opt_restart = True
            argi += 1
        elif sys.argv[argi] == "--reject-output" and argi + 1 < len(sys.argv):
            rejected_outputs.append(sys.argv[argi + 1])
            argi += 2
        else:
            print("usage: {} <model name> [--expect-error <substring>] [--expect-nqp-full <n>] [--compare-no-bf-energy] [--compare-no-bf-twobodyg] [--compare-no-bf-twobodygex] [--compare-no-bf-gradient] [--compare-proj-bf-finite-diff] [--check-bf-green2-bruteforce] [--use-nonidentity-init] [--set-ncond <n>] [--expect-all-complex-flag <0|1>] [--compare-real-complex-nonidentity <complex model>] [--check-opt-output-restart] [--reject-output <substring>]".format(sys.argv[0]))
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
    if check_opt_restart:
        return check_opt_output_restart(rootdir, model, mpi_procs, rejected_outputs)

    work_suffix = ""
    if compare_energy:
        work_suffix += "_energy"
    if compare_twobodyg:
        work_suffix += "_twobodyg"
    if compare_twobodygex:
        work_suffix += "_twobodygex"
    if compare_gradient:
        work_suffix += "_gradient"
    if compare_proj_bf_fd:
        work_suffix += "_projbf_fd"
    if check_bf_green2_bruteforce:
        work_suffix += "_green2_bruteforce"
    if use_nonidentity_init:
        work_suffix += "_nonidentity"
    if ncond_override is not None:
        work_suffix += "_ncond{}".format(ncond_override)
    workdir = os.path.join(rootdir, "work", model + work_suffix)
    if mpi_procs:
        workdir += "_mpi{}".format(mpi_procs)

    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)

    copy_def_files(refdir, workdir, include_backflow=True)
    if ncond_override is not None:
        update_modpara(os.path.join(workdir, "modpara.def"), {"Ncond": ncond_override})

    nsite = parse_nsite(os.path.join(workdir, "modpara.def"))
    definition = build_chain_nn_backflow(length=nsite, optimize=compare_proj_bf_fd)
    write_chain_nn_backflow(workdir, length=nsite, optimize=compare_proj_bf_fd)
    if compare_twobodyg or check_bf_green2_bruteforce:
        write_minimal_twobodyg(workdir, nsite)
    if compare_twobodygex:
        write_minimal_twobodygex(workdir, nsite)

    init_path = None
    if use_nonidentity_init:
        assert_nonidentity_init_layout(workdir)
        nslater = parse_norbitalidx(os.path.join(workdir, "orbitalidx.def"))
        write_orbital_opt_flags(os.path.join(workdir, "orbitalidx.def"), nsite, nslater, 0)
        init_name = "nonidentity_init.dat"
        projbf = write_nonidentity_init_parameter(
            os.path.join(workdir, init_name),
            definition.n_proj_bf,
            nslater,
        )
        if projbf[0] == 1.0 or all(value == 0.0 for value in projbf[1:]):
            print("ERROR: non-identity ProjBF initialization is invalid.")
            return -1
        init_path = init_name

    dump_path = None if use_nonidentity_init else os.path.join(workdir, "bf_identity_dump.dat")
    diff_dump_path = os.path.join(workdir, "bf_diff_dump.dat") if compare_gradient else None
    fd_dump_path = os.path.join(workdir, "bf_projbf_fd_dump.dat") if compare_proj_bf_fd else None
    green2_dump_path = os.path.join(workdir, "bf_green2_bruteforce_dump.dat") if check_bf_green2_bruteforce else None
    proc = run_vmc(rootdir, workdir, mpi_procs, dump_path=dump_path, diff_dump_path=diff_dump_path,
                   fd_dump_path=fd_dump_path, green2_dump_path=green2_dump_path, init_path=init_path)

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
    tol = 1.0e-10
    if dump_path is not None:
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
    elif expected_nqp_full is not None:
        print("ERROR: --expect-nqp-full requires an identity BackFlow dump.")
        return -1
    if not first_row_is_finite(os.path.join(workdir, "output", "zvo_out_001.dat")):
        print("ERROR: zvo_out_001.dat is missing or non-finite.")
        return -1
    if compare_energy:
        result = compare_no_bf_energy(rootdir, refdir, workdir, mpi_procs, tol, rejected_outputs)
        if result != 0:
            return result
    if compare_twobodyg:
        result = compare_no_bf_twobodyg(rootdir, refdir, workdir, mpi_procs, tol, rejected_outputs)
        if result != 0:
            return result
    if compare_twobodygex:
        result = compare_no_bf_twobodygex(rootdir, refdir, workdir, mpi_procs, tol, rejected_outputs)
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
    if check_bf_green2_bruteforce:
        result = check_bf_green2_bruteforce_dump(green2_dump_path, 1.0e-10, expected_all_complex_flag)
        if result != 0:
            return result

    return 0


if __name__ == "__main__":
    sys.exit(main())
