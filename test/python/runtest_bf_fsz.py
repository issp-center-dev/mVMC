from __future__ import print_function

import math
import os
import shutil
import subprocess
import sys

from backflow_def_helper import build_chain_nn_backflow, write_chain_nn_backflow


def copy_base_defs(refdir, workdir):
    for fn in ("modpara.def", "locspn.def", "trans.def", "qptransidx.def"):
        shutil.copy(os.path.join(refdir, fn), os.path.join(workdir, fn))


def write_namelist(workdir, include_backflow):
    lines = [
        "         ModPara  modpara.def",
        "         LocSpin  locspn.def",
        "           Trans  trans.def",
        "        OneBodyG  greenone.def",
        "        TwoBodyG  greentwo.def",
        " OrbitalGeneral  orbitalidxgen.def",
        "        TransSym  qptransidx.def",
    ]
    if include_backflow:
        lines.extend([
            "              BF  bf.def",
            "         BFRange  rangebf.def",
        ])
    with open(os.path.join(workdir, "namelist.def"), "w") as fp:
        fp.write("\n".join(lines) + "\n")


def write_orbital_general(workdir, nsite):
    nsite2 = 2 * nsite
    norbital = nsite2 * (nsite2 - 1) // 2
    lines = [
        "=============================================",
        "NOrbitalIdx         {}".format(norbital),
        "ComplexType          1",
        "=============================================",
        "=============================================",
    ]
    idx = 0
    for orb_i in range(nsite2):
        for orb_j in range(orb_i + 1, nsite2):
            lines.append(
                "{:5d} {:2d} {:6d} {:2d} {:6d}      1".format(
                    orb_i % nsite, orb_i // nsite,
                    orb_j % nsite, orb_j // nsite,
                    idx,
                )
            )
            idx += 1
    for idx in range(norbital):
        lines.append("{:5d}      1".format(idx))
    with open(os.path.join(workdir, "orbitalidxgen.def"), "w") as fp:
        fp.write("\n".join(lines) + "\n")


def write_green_defs(workdir):
    one_body_rows = [
        (0, 0, 0, 0),
        (0, 0, 1, 0),
        (0, 1, 0, 1),
        (0, 1, 1, 1),
    ]
    with open(os.path.join(workdir, "greenone.def"), "w") as fp:
        fp.write("===============================\n")
        fp.write("NCisAjs         {}\n".format(len(one_body_rows)))
        fp.write("===============================\n")
        fp.write("======== Green functions ======\n")
        fp.write("===============================\n")
        for row in one_body_rows:
            fp.write(" ".join("{:5d}".format(value) for value in row) + "\n")

    two_body_rows = [
        (0, 0, 0, 0, 0, 1, 0, 1),
        (0, 0, 1, 0, 2, 1, 3, 1),
        (0, 1, 1, 1, 2, 0, 3, 0),
    ]
    with open(os.path.join(workdir, "greentwo.def"), "w") as fp:
        fp.write("=============================================\n")
        fp.write("NCisAjsCktAltDC         {}\n".format(len(two_body_rows)))
        fp.write("=============================================\n")
        fp.write("======== Green functions for BF-FSZ ========\n")
        fp.write("=============================================\n")
        for row in two_body_rows:
            fp.write(" ".join("{:5d}".format(value) for value in row) + "\n")


def replace_data_line(workdir, filename, row):
    path = os.path.join(workdir, filename)
    with open(path) as fp:
        lines = fp.readlines()
    if len(lines) < 6:
        raise RuntimeError("{} has no data row".format(filename))
    lines[5] = row
    with open(path, "w") as fp:
        fp.writelines(lines)


def make_spin_changing_transfer(workdir):
    replace_data_line(
        workdir,
        "trans.def",
        "    1     0     0     1         1.000000000000000         0.000000000000000\n",
    )


def make_spin_changing_one_body_g(workdir):
    replace_data_line(
        workdir,
        "greenone.def",
        "    0     0     1     1\n",
    )


def make_spin_changing_two_body_g(workdir):
    replace_data_line(
        workdir,
        "greentwo.def",
        "    0     0     0     1     0     1     0     1\n",
    )


def make_spin_changing_two_body_g_ex(workdir):
    rows = [
        (0, 0, 1, 1, 2, 0, 3, 0),
    ]
    with open(os.path.join(workdir, "greentwoex.def"), "w") as fp:
        fp.write("=============================================\n")
        fp.write("NCisAjsCktAlt         {}\n".format(len(rows)))
        fp.write("=============================================\n")
        fp.write("======== Green functions for BF-FSZ ========\n")
        fp.write("=============================================\n")
        for row in rows:
            fp.write(" ".join("{:5d}".format(value) for value in row) + "\n")

    namelist_path = os.path.join(workdir, "namelist.def")
    with open(namelist_path) as fp:
        lines = fp.readlines()
    if not any(line.split() and line.split()[0] == "TwoBodyGEx" for line in lines):
        lines.append("      TwoBodyGEx  greentwoex.def\n")
    with open(namelist_path, "w") as fp:
        fp.writelines(lines)


def make_momentum_projection(workdir):
    with open(os.path.join(workdir, "qptransidx.def"), "w") as fp:
        fp.write("=============================================\n")
        fp.write("NQPTrans          2\n")
        fp.write("=============================================\n")
        fp.write("======== TrIdx_TrWeight_and_TrIdx_i_xi ======\n")
        fp.write("=============================================\n")
        fp.write("0    1.00000    0.00000\n")
        fp.write("1    1.00000    0.00000\n")
        fp.write("    0      0      0      1\n")
        fp.write("    0      1      1      1\n")
        fp.write("    0      2      2      1\n")
        fp.write("    0      3      3      1\n")
        fp.write("    1      0      1      1\n")
        fp.write("    1      1      2      1\n")
        fp.write("    1      2      3      1\n")
        fp.write("    1      3      0      1\n")


def write_locspn(workdir, local_sites, nsite=4):
    local_sites = set(local_sites)
    with open(os.path.join(workdir, "locspn.def"), "w") as fp:
        fp.write("================================\n")
        fp.write("NlocalSpin     {}\n".format(len(local_sites)))
        fp.write("================================\n")
        fp.write("========i_1LocSpn_0IteElc ======\n")
        fp.write("================================\n")
        for site in range(nsite):
            fp.write("{:5d} {:6d}\n".format(site, 1 if site in local_sites else 0))


def update_modpara(workdir, updates):
    path = os.path.join(workdir, "modpara.def")
    found = set()
    lines = []
    with open(path) as fp:
        for line in fp:
            cols = line.split()
            if cols and cols[0] in updates:
                lines.append("{:<15} {}\n".format(cols[0], updates[cols[0]]))
                found.add(cols[0])
            else:
                lines.append(line)
    for key, value in updates.items():
        if key not in found:
            lines.append("{:<15} {}\n".format(key, value))
    with open(path, "w") as fp:
        fp.writelines(lines)


def parse_norbitalidx(path):
    with open(path) as fp:
        for line in fp:
            cols = line.split()
            if len(cols) >= 2 and cols[0] == "NOrbitalIdx":
                return int(cols[1])
    raise RuntimeError("NOrbitalIdx was not found in {}".format(path))


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
            real = 1.0 + 0.1 * row if row == col else 0.03 * float(row - col)
        else:
            real = math.sin(0.37 * (idx + 1))
        values.extend([real, 0.0, 0.0])

    with open(path, "w") as fp:
        fp.write(" ".join("{:.18e}".format(value) for value in values))
        fp.write("\n")


def write_nonidentity_init(workdir, nsite=4):
    definition = build_chain_nn_backflow(length=nsite, optimize=False)
    nslater = parse_norbitalidx(os.path.join(workdir, "orbitalidxgen.def"))
    init_name = "nonidentity_init.dat"
    write_nonidentity_init_parameter(
        os.path.join(workdir, init_name),
        definition.n_proj_bf,
        nslater,
    )
    return init_name


def run_vmc(rootdir, workdir, mpi_procs=None, init_path=None):
    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
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
    )
    with open(os.path.join(workdir, "vmc.log"), "w") as fp:
        fp.write(proc.stdout)
    return proc


def read_first_row(path):
    with open(path) as fp:
        for line in fp:
            cols = line.split()
            if cols:
                return [float(x) for x in cols]
    raise RuntimeError("no numeric row in {}".format(path))


def read_float_rows(path):
    rows = []
    with open(path) as fp:
        for line in fp:
            cols = line.split()
            if cols:
                rows.append([float(x) for x in cols])
    return rows


def assert_finite_nonzero_rows(path, label):
    rows = read_float_rows(path)
    if not rows:
        print("ERROR: {} output is empty.".format(label))
        return -1
    max_abs = 0.0
    for row in rows:
        for value in row:
            if not math.isfinite(value):
                print("ERROR: {} output contains non-finite values.".format(label))
                return -1
            max_abs = max(max_abs, abs(value))
    if max_abs <= 1.0e-12:
        print("ERROR: {} output is vacuous: all values are zero.".format(label))
        return -1
    return 0


def max_abs_diff(left_rows, right_rows):
    if len(left_rows) != len(right_rows):
        raise RuntimeError("row count mismatch: left={} right={}".format(len(left_rows), len(right_rows)))
    max_diff = 0.0
    for left, right in zip(left_rows, right_rows):
        if len(left) != len(right):
            raise RuntimeError("column count mismatch: left={} right={}".format(len(left), len(right)))
        for x, y in zip(left, right):
            max_diff = max(max_diff, abs(x-y))
    return max_diff


def prepare_case(rootdir, name, include_backflow):
    refdir = os.path.join(rootdir, "data", "BackFlow_Identity_Complex")
    workdir = os.path.join(rootdir, "work", name)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)

    copy_base_defs(refdir, workdir)
    write_namelist(workdir, include_backflow)
    write_orbital_general(workdir, nsite=4)
    write_green_defs(workdir)
    if include_backflow:
        write_chain_nn_backflow(workdir, length=4, optimize=False)
    return workdir


def expect_invalid(rootdir, name, updates, expected, mpi_procs=None, local_sites=None, mutate_defs=None):
    workdir = prepare_case(rootdir, name, True)
    update_modpara(workdir, updates)
    if local_sites is not None:
        write_locspn(workdir, local_sites)
    if mutate_defs is not None:
        mutate_defs(workdir)
    proc = run_vmc(rootdir, workdir, mpi_procs)
    if proc.returncode == 0:
        print("ERROR: {} unexpectedly succeeded".format(name))
        return -1
    if expected not in proc.stdout:
        print("ERROR: {} did not report expected error: {}".format(name, expected))
        print(proc.stdout)
        return -1
    if "Start: Sampling." in proc.stdout:
        print("ERROR: {} reached sampling before failing".format(name))
        print(proc.stdout)
        return -1
    return 0


def get_named_invalid_case(case_name):
    invalid_cases = {
        "BackFlow_FSZ_InvalidSpinChangingTransfer": (
            {},
            "BackFlow FSZ supports only spin-conserving Transfer terms",
            None,
            make_spin_changing_transfer,
        ),
        "BackFlow_FSZ_InvalidSpinChangingOneBodyG": (
            {},
            "BackFlow FSZ supports only spin-conserving OneBodyG terms",
            None,
            make_spin_changing_one_body_g,
        ),
        "BackFlow_FSZ_InvalidSpinChangingOneBodyG_mpi": (
            {},
            "BackFlow FSZ supports only spin-conserving OneBodyG terms",
            None,
            make_spin_changing_one_body_g,
        ),
        "BackFlow_FSZ_InvalidSpinChangingTwoBodyG": (
            {},
            "BackFlow FSZ supports only spin-conserving TwoBodyG terms",
            None,
            make_spin_changing_two_body_g,
        ),
        "BackFlow_FSZ_InvalidSpinChangingTwoBodyGEx": (
            {},
            "including entries derived from TwoBodyGEx",
            None,
            make_spin_changing_two_body_g_ex,
        ),
    }
    return invalid_cases.get(case_name)


def main():
    case_name = sys.argv[1] if len(sys.argv) > 1 else "BackFlow_FSZ_IdentityEnergy_Complex"
    rootdir = os.getcwd()
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    momentum_projection_cases = (
        "BackFlow_FSZ_MomentumProjection_Complex",
        "BackFlow_FSZ_MomentumProjection_Complex_mpi",
        "BackFlow_FSZ_MomentumProjection_NonIdentity_Complex",
    )
    nonidentity_momentum_case = case_name == "BackFlow_FSZ_MomentumProjection_NonIdentity_Complex"
    named_invalid = get_named_invalid_case(case_name)
    if named_invalid is not None:
        updates, expected, local_sites, mutate_defs = named_invalid
        return expect_invalid(rootdir, case_name, updates, expected, mpi_procs, local_sites, mutate_defs)

    bf_workdir = prepare_case(rootdir, case_name, True)
    no_bf_workdir = prepare_case(rootdir, case_name + "_nobf", False)
    if case_name in momentum_projection_cases:
        for workdir in (bf_workdir, no_bf_workdir):
            update_modpara(workdir, {"NMPTrans": "2"})
            make_momentum_projection(workdir)
    init_path = write_nonidentity_init(bf_workdir) if nonidentity_momentum_case else None

    bf_proc = run_vmc(rootdir, bf_workdir, mpi_procs, init_path)
    if bf_proc.returncode != 0:
        print(bf_proc.stdout)
        return bf_proc.returncode
    if nonidentity_momentum_case:
        for filename, label in (
            ("zvo_out_001.dat", "zvo_out"),
            ("zvo_cisajs_001.dat", "OneBodyG"),
            ("zvo_cisajscktalt_001.dat", "TwoBodyG"),
        ):
            status = assert_finite_nonzero_rows(os.path.join(bf_workdir, "output", filename), label)
            if status != 0:
                return status
        return 0

    no_bf_proc = run_vmc(rootdir, no_bf_workdir, mpi_procs)
    if no_bf_proc.returncode != 0:
        print(no_bf_proc.stdout)
        return no_bf_proc.returncode

    bf_row = read_first_row(os.path.join(bf_workdir, "output", "zvo_out_001.dat"))
    no_bf_row = read_first_row(os.path.join(no_bf_workdir, "output", "zvo_out_001.dat"))
    if len(bf_row) != len(no_bf_row):
        print("ERROR: zvo_out column mismatch: bf={} no_bf={}".format(len(bf_row), len(no_bf_row)))
        return -1

    max_diff = max(abs(x - y) for x, y in zip(bf_row, no_bf_row))
    if max_diff > 1.0e-10:
        print("ERROR: BackFlow FSZ identity energy mismatch: max_abs_diff={:.3e}".format(max_diff))
        print("bf    {}".format(bf_row))
        print("no-bf {}".format(no_bf_row))
        return -1

    comparisons = [
        ("zvo_cisajs_001.dat", "OneBodyG"),
        ("zvo_cisajscktalt_001.dat", "TwoBodyG"),
    ]
    for filename, label in comparisons:
        bf_rows = read_float_rows(os.path.join(bf_workdir, "output", filename))
        no_bf_rows = read_float_rows(os.path.join(no_bf_workdir, "output", filename))
        diff = max_abs_diff(bf_rows, no_bf_rows)
        if diff > 1.0e-10:
            print("ERROR: BackFlow FSZ identity {} mismatch: max_abs_diff={:.3e}".format(label, diff))
            return -1

    invalid_cases = [
        (
            case_name + "_InvalidOptimization",
            {"NVMCCalMode": "0"},
            "BackFlow FSZ optimization is not implemented yet",
            None,
        ),
        (
            case_name + "_InvalidLocSpin",
            {"NExUpdatePath": "1"},
            "BackFlow MVP does not support NLocalSpin > 0",
            [0],
        ),
    ]
    for name, updates, expected, local_sites in invalid_cases:
        status = expect_invalid(rootdir, name, updates, expected, mpi_procs, local_sites)
        if status != 0:
            return status

    return 0


if __name__ == "__main__":
    sys.exit(main())
