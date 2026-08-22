from __future__ import print_function

import os
import shutil
import sys

import numpy as np

from runtest_bf_fsz import run_vmc, update_modpara
from runtest_lanczos2_real_ed import stochastic_hankel_standard_error


MODEL = "HubbardTetragonal_MomentumProjection"


def modpara_values(path):
    values = {}
    with open(path) as source:
        for line in source:
            columns = line.split()
            if len(columns) >= 2:
                values[columns[0]] = columns[1]
    return values


def read_oracle(path, nsite):
    rows = []
    with open(path) as source:
        for line in source:
            columns = line.split()
            marker = 3 + 2 * nsite
            if (
                len(columns) != marker + 1 + 8
                or columns[0] != "sample"
                or columns[2] != "occ"
                or columns[marker] != "ls2power"
            ):
                raise AssertionError(
                    "invalid projection oracle row: {}".format(line.rstrip())
                )
            occupation = tuple(
                int(value) for value in columns[3:marker]
            )
            values = [float(value) for value in columns[marker + 1:]]
            power = np.array(
                [
                    complex(values[2 * index], values[2 * index + 1])
                    for index in range(4)
                ],
                dtype=complex,
            )
            if not np.all(np.isfinite(power)):
                raise AssertionError("projection oracle has non-finite power")
            rows.append((occupation, power))
    if not rows:
        raise AssertionError("projection oracle dump is empty")
    return rows


def read_moment(path):
    with open(path) as source:
        lines = source.readlines()
    prefix = "# mVMC ls2_moment v2 basis_shift="
    if len(lines) != 2 or not lines[0].startswith(prefix):
        raise AssertionError("projection moment output has invalid layout")
    shift = float(lines[0][len(prefix):].split()[0])
    values = [float(value) for value in lines[1].split()]
    if len(values) != 32:
        raise AssertionError("projection moment output has invalid width")
    moment = np.array(
        [
            complex(values[2 * index], values[2 * index + 1])
            for index in range(16)
        ],
        dtype=complex,
    ).reshape((4, 4))
    return shift, moment


def shift_local_power(power, shift):
    return np.array(
        [
            power[0],
            power[1] - shift * power[0],
            power[2] - 2.0 * shift * power[1]
            + shift * shift * power[0],
            power[3] - 3.0 * shift * power[2]
            + 3.0 * shift * shift * power[1]
            - shift * shift * shift * power[0],
        ],
        dtype=complex,
    )


def main():
    rootdir = os.getcwd()
    source_workdir = os.path.join(rootdir, "work", MODEL)
    if not os.path.isdir(source_workdir):
        raise AssertionError(
            "{} setup fixture did not create its work directory".format(MODEL)
        )

    workdir = os.path.join(rootdir, "work", "Lanczos2_Projection_Smoke")
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    shutil.copytree(source_workdir, workdir)
    output_dir = os.path.join(workdir, "output")
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    update_modpara(
        workdir,
        {
            "NLanczosMode": "1",
            "NLanczosEstimatorMode": "1",
            "NLanczosStep": "2",
            "NVMCCalMode": "1",
            "NDataQtySmp": "1",
            "NVMCSample": "128",
            "NVMCWarmUp": "32",
            "NVMCInterval": "2",
            "NSplitSize": "1",
        },
    )
    parameters = modpara_values(os.path.join(workdir, "modpara.def"))
    nsite = int(parameters["Nsite"])
    if int(parameters["NMPTrans"]) <= 1:
        raise AssertionError("projection smoke does not use momentum projection")
    if int(parameters["NSPGaussLeg"]) <= 1:
        raise AssertionError("projection smoke does not use spin projection")

    initial_path = os.path.join(rootdir, "data", MODEL, "initial.def")
    dump_name = "lanczos2_projection_oracle.dat"
    process = run_vmc(
        rootdir,
        workdir,
        init_path=initial_path,
        extra_env={"MVMC_LANCZOS_ORACLE_DUMP": dump_name},
    )
    if process.returncode != 0:
        raise AssertionError(
            "Lanczos2 projection vmc.out failed:\n{}".format(process.stdout)
        )

    rows = read_oracle(os.path.join(workdir, dump_name), nsite)
    occupation_count = len(set(occupation for occupation, _ in rows))
    if occupation_count < 8:
        raise AssertionError("projection smoke sampled too few occupations")
    basis_shift, observed_moment = read_moment(
        os.path.join(output_dir, "zvo_ls2_moment_001.dat")
    )
    expected_moment = sum(
        np.outer(
            np.conj(shift_local_power(power, basis_shift)),
            shift_local_power(power, basis_shift),
        )
        for _, power in rows
    ) / float(len(rows))
    if not np.allclose(
        observed_moment, expected_moment, rtol=3.0e-12, atol=3.0e-12
    ):
        raise AssertionError(
            "projection shifted moment does not match sample-local powers"
        )
    if not np.allclose(
        observed_moment, np.conj(observed_moment.T),
        rtol=3.0e-12, atol=3.0e-12,
    ):
        raise AssertionError("projection shifted moment is not Hermitian")
    hankel_standard_error = stochastic_hankel_standard_error(
        [power for _, power in rows], basis_shift
    )

    print(
        "Lanczos2 projection smoke: PASS "
        "({} samples, {} occupations, NMPTrans={}, NSPGaussLeg={}, "
        "max Hankel SE={:.3f})".format(
            len(rows),
            occupation_count,
            parameters["NMPTrans"],
            parameters["NSPGaussLeg"],
            hankel_standard_error,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AssertionError as error:
        print("ERROR: {}".format(error))
        sys.exit(1)
