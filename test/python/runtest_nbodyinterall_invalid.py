from __future__ import print_function

import os
import shutil
import subprocess
import sys

from nbodyinterall_exact_oracle import (
    complex_antiparallel_state,
    write_common_nbodyinterall_defs,
    write_nbodyinterall_def,
)
from nbodyg_exact_oracle import write, write_antiparallel_orbital


def write_real_antiparallel_orbital(workdir, nsite):
    lines = [
        "=============================================",
        "NOrbitalIdx         {}".format(nsite * nsite),
        "ComplexType          0",
        "=============================================",
        "=============================================",
    ]
    idx = 0
    for i in range(nsite):
        for j in range(nsite):
            lines.append("{:5d} {:6d} {:6d}      1".format(i, j, idx))
            idx += 1
    for idx in range(nsite * nsite):
        lines.append("{:5d}      1".format(idx))
    write(os.path.join(workdir, "orbitalidx.def"), "\n".join(lines) + "\n")


def write_base_case(workdir, orbital_complex=True, lanczos_mode=0, backflow=False):
    nsite = 6
    f_matrix, _ = complex_antiparallel_state(nsite)
    write_common_nbodyinterall_defs(workdir, nsite, 20, 77889,
                                    "Orbital", "orbitalidx.def")
    if orbital_complex:
        write_antiparallel_orbital(workdir, f_matrix)
    else:
        write_real_antiparallel_orbital(workdir, nsite)

    modpara = os.path.join(workdir, "modpara.def")
    text = open(modpara).read()
    text = text.replace("NLanczosMode   0", "NLanczosMode   {}".format(lanczos_mode))
    write(modpara, text)

    if backflow:
        namelist = os.path.join(workdir, "namelist.def")
        text = open(namelist).read() + "              BF  backflow.def\n"
        write(namelist, text)
        write(
            os.path.join(workdir, "backflow.def"),
            (
                "=============================================\n"
                "NBackFlowIdx     1\n"
                "=============================================\n"
                "======== BackFlow parameters ================\n"
                "=============================================\n"
            ),
        )


def write_raw_nbodyinterall(workdir, count, body):
    write(
        os.path.join(workdir, "nbodyinterall.def"),
        (
            "=============================================\n"
            "NNBodyInterAll   {0}\n"
            "=============================================\n"
            "======== NBodyInterAll interactions =========\n"
            "=============================================\n"
            "{1}"
        ).format(count, body),
    )


def case_negative_count(workdir):
    write_base_case(workdir)
    write_raw_nbodyinterall(workdir, -1, "")


def case_zero_order(workdir):
    write_base_case(workdir)
    write_raw_nbodyinterall(workdir, 1, "0 1.0 0.0\n")


def case_field_count(workdir):
    write_base_case(workdir)
    write_raw_nbodyinterall(workdir, 1, "1 0 0 0 0\n")


def case_extra_token(workdir):
    write_base_case(workdir)
    write_raw_nbodyinterall(workdir, 1, "1 0 0 0 0 1.0 0.0 extra\n")


def case_site(workdir):
    write_base_case(workdir)
    write_raw_nbodyinterall(workdir, 1, "1 99 0 0 0 1.0 0.0\n")


def case_spin(workdir):
    write_base_case(workdir)
    write_raw_nbodyinterall(workdir, 1, "1 0 2 0 0 1.0 0.0\n")


def case_coeff(workdir):
    write_base_case(workdir)
    write_raw_nbodyinterall(workdir, 1, "1 0 0 0 0 nan 0.0\n")


def case_spin_change_non_fsz(workdir):
    write_base_case(workdir)
    write_raw_nbodyinterall(workdir, 1, "1 0 0 0 1 1.0 0.0\n")


def case_real_path(workdir):
    write_base_case(workdir, orbital_complex=False)
    write_nbodyinterall_def(workdir, [(1, (0, 0, 0, 0), 1.0 + 0.0j)])


def case_backflow(workdir):
    write_base_case(workdir, backflow=True)
    write_nbodyinterall_def(workdir, [(1, (0, 0, 0, 0), 1.0 + 0.0j)])


def case_lanczos(workdir):
    write_base_case(workdir, lanczos_mode=1)
    write_nbodyinterall_def(workdir, [(1, (0, 0, 0, 0), 1.0 + 0.0j)])


CASES = {
    "NBodyInterAll_InvalidNegativeCount": case_negative_count,
    "NBodyInterAll_InvalidZeroOrder": case_zero_order,
    "NBodyInterAll_InvalidFieldCount": case_field_count,
    "NBodyInterAll_InvalidExtraToken": case_extra_token,
    "NBodyInterAll_InvalidSite": case_site,
    "NBodyInterAll_InvalidSpin": case_spin,
    "NBodyInterAll_InvalidCoeff": case_coeff,
    "NBodyInterAll_InvalidSpinChangeNonFsz": case_spin_change_non_fsz,
    "NBodyInterAll_InvalidRealPath": case_real_path,
    "NBodyInterAll_InvalidBackFlow": case_backflow,
    "NBodyInterAll_InvalidLanczos": case_lanczos,
}


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in CASES:
        print("usage: {} <model name> <expected_error_substring>".format(sys.argv[0]))
        return -1

    model = sys.argv[1]
    expected = sys.argv[2]
    rootdir = os.getcwd()
    workdir = os.path.join(rootdir, "work", model)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    CASES[model](workdir)

    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    proc = subprocess.run(
        [bin_to_test, "-e", "namelist.def"],
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    with open(os.path.join(workdir, "invalid_test.log"), "w") as f:
        f.write(proc.stdout)

    if expected not in proc.stdout:
        print("ERROR: expected error substring not found: {}".format(expected))
        print("---- output begin ----")
        print(proc.stdout)
        print("---- output end ----")
        return -1
    if proc.returncode == 0:
        print("ERROR: vmc.out unexpectedly succeeded for invalid input.")
        return -1
    return 0


if __name__ == "__main__":
    sys.exit(main())
