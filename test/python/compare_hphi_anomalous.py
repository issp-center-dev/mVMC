from __future__ import print_function

import argparse
import os
import shutil
import subprocess
import sys

import numpy as np

from grandcanonical_exact_oracle import (
    ANOMALOUS_G_OPERATORS,
    NORBITAL,
    NSITE,
    ed_ground_state,
    even_sector_eigenvalues,
    sector_hamiltonian,
)
from runtest_gc import write, write_anomalous, write_hamiltonian


HPHI_REFERENCE_COMMIT = "481ac1075c7caebbf0173a77723f611a3252bcae"
TOLERANCE = 1.0e-8


def match_spectrum(expected, actual, tolerance=TOLERANCE):
    unused = list(float(value) for value in actual)
    for expected_value in sorted(float(value) for value in expected):
        if not unused:
            raise AssertionError("HPhi spectrum ran out of unmatched values")
        index = min(range(len(unused)),
                    key=lambda candidate: abs(unused[candidate] - expected_value))
        difference = abs(unused[index] - expected_value)
        if difference > tolerance:
            raise AssertionError(
                "HPhi spectrum has no one-to-one match for {} (nearest {}, diff {})"
                .format(expected_value, unused[index], difference))
        unused.pop(index)


def assert_complex_maps_close(expected, actual, tolerance=TOLERANCE):
    if set(expected) != set(actual):
        raise AssertionError("AnomalousG operator inventories differ")
    for key, expected_value in expected.items():
        difference = abs(actual[key] - expected_value)
        if difference > tolerance:
            raise AssertionError(
                "AnomalousG {} mismatch: actual={} expected={} diff={}"
                .format(key, actual[key], expected_value, difference))


def self_test():
    match_spectrum((0.0, 1.0), (1.0, 0.0, 2.0), 1.0e-12)
    try:
        match_spectrum((0.0, 1.0), (0.0, 0.0), 1.0e-12)
    except AssertionError:
        pass
    else:
        raise AssertionError(
            "one-to-one spectrum matcher accepted a duplicated eigenvalue")

    expected = {(1, 0, 3): 0.2 + 0.3j}
    try:
        assert_complex_maps_close(
            expected, {(1, 0, 3): expected[(1, 0, 3)].conjugate()},
            1.0e-12)
    except AssertionError:
        pass
    else:
        raise AssertionError("AnomalousG matcher accepted complex conjugation")


def write_hphi_control_files(workdir):
    write(
        os.path.join(workdir, "calcmod.def"),
        "# CalcType: 2=FullDiag; CalcModel: 3=HubbardGC\n"
        "CalcType   2\n"
        "CalcModel  3\n"
        "ReStart    0\n"
        "CalcSpec   0\n"
        "CalcEigenVec 0\n"
        "InitialVecType 0\n"
        "InputEigenVec 0\n"
        "OutputEigenVec 0\n"
        "InputHam 0\n"
        "OutputHam 0\n"
        "OutputExVec 0\n",
    )
    write(
        os.path.join(workdir, "modpara.def"),
        "--------------------\n"
        "Model_Parameters   0\n"
        "--------------------\n"
        "HPhi_Cal_Parameters\n"
        "--------------------\n"
        "CDataFileHead  zvo\n"
        "CParaFileHead  zqp\n"
        "--------------------\n"
        "Nsite          2\n"
        "Lanczos_max    64\n"
        "initial_iv     1\n"
        "exct           1\n"
        "LanczosEps     14\n"
        "LanczosTarget  2\n"
        "LargeValue     0.0\n"
        "NumAve         1\n"
        "ExpecInterval  1\n"
        "NOmega         8\n"
        "OmegaMax       0.0 0.0\n"
        "OmegaMin       0.0 0.0\n"
        "OmegaOrg       0.0 0.0\n"
        "PreCG          1\n",
    )
    write(
        os.path.join(workdir, "locspn.def"),
        "================================\n"
        "NlocalSpin     0\n"
        "================================\n"
        "========i_1LocSpn_0IteElc======\n"
        "================================\n"
        "0 0\n"
        "1 0\n",
    )
    namelist = (
        "         ModPara  modpara.def",
        "         LocSpin  locspn.def",
        "           Trans  trans.def",
        "    CoulombIntra  coulombintra.def",
        "    CoulombInter  coulombinter.def",
        "            Hund  hund.def",
        "         PairHop  pairhopp.def",
        "        Exchange  exchange.def",
        "        InterAll  interall.def",
        "  NBodyInterAll  nbodyinterall.def",
        "   AnomalousTerm  anomalousterm.def",
        "      AnomalousG  anomalousg.def",
        "         CalcMod  calcmod.def",
    )
    write(os.path.join(workdir, "namelist.def"), "\n".join(namelist) + "\n")


def prepare_fixture(workdir, mu, delta):
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    write_hphi_control_files(workdir)
    write_hamiltonian(workdir, mu)
    write_anomalous(workdir, delta)


def parse_hphi_spectrum(path):
    eigenvalues = []
    with open(path) as stream:
        for line_number, line in enumerate(stream, 1):
            columns = line.split()
            if line_number == 1 or not columns:
                continue
            try:
                eigenvalues.append(float(columns[0]))
            except (ValueError, IndexError):
                raise AssertionError("malformed HPhi physical row: {}".format(line))
    if not eigenvalues:
        raise AssertionError("HPhi FullDiag spectrum is empty")
    return eigenvalues


def parse_hphi_eigenvalues(path):
    eigenvalues = []
    with open(path) as stream:
        for line in stream:
            columns = line.split()
            if not columns:
                continue
            if len(columns) != 2:
                raise AssertionError("malformed HPhi eigenvalue row: {}".format(line))
            eigenvalues.append(float(columns[1]))
    if not eigenvalues:
        raise AssertionError("HPhi high-precision spectrum is empty")
    return eigenvalues


def parse_hphi_anomalous(path):
    values = {}
    with open(path) as stream:
        for line in stream:
            columns = line.split()
            if not columns:
                continue
            if len(columns) != 7:
                raise AssertionError("malformed HPhi AnomalousG row: {}".format(line))
            key = (
                int(columns[0]),
                int(columns[1]) + int(columns[2]) * NSITE,
                int(columns[3]) + int(columns[4]) * NSITE,
            )
            if key in values:
                raise AssertionError("duplicate HPhi AnomalousG row: {}".format(key))
            values[key] = complex(float(columns[5]), float(columns[6]))
    return values


def source_directory_from_cmake_cache(binary):
    build_root = os.path.dirname(os.path.dirname(os.path.realpath(binary)))
    cache = os.path.join(build_root, "CMakeCache.txt")
    if not os.path.isfile(cache):
        return None
    with open(cache) as stream:
        for line in stream:
            prefix = "CMAKE_HOME_DIRECTORY:INTERNAL="
            if line.startswith(prefix):
                return line[len(prefix):].strip()
    return None


def hphi_commit(binary):
    explicit = os.environ.get("MVMC_HPHI_COMMIT")
    if explicit:
        return explicit
    source = source_directory_from_cmake_cache(binary)
    if source is None:
        return "unknown"
    process = subprocess.run(
        ["git", "-C", source, "rev-parse", "HEAD"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True)
    return process.stdout.strip() if process.returncode == 0 else "unknown"


def run_cross_check(binary, workdir):
    mu = 0.4
    delta = 0.35 - 0.20j
    odd_basis = tuple(
        state for state in range(1 << NORBITAL)
        if bin(state).count("1") % 2 == 1
    )
    even_spectrum = even_sector_eigenvalues(mu, delta)
    odd_spectrum = np.linalg.eigvalsh(sector_hamiltonian(
        odd_basis, mu, delta))
    if even_spectrum[0] + 1.0e-6 >= odd_spectrum[0]:
        raise AssertionError(
            "HPhi fixture ground state is not robustly even parity")
    # AnomalousG expectation values are basis dependent inside a degenerate
    # ground-state manifold, so require a clear gap in the even sector.
    if even_spectrum[1] - even_spectrum[0] <= 1.0e-6:
        raise AssertionError(
            "HPhi fixture even-sector ground state is degenerate")

    prepare_fixture(workdir, mu, delta)
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment["OPENBLAS_NUM_THREADS"] = "1"
    environment["MKL_NUM_THREADS"] = "1"
    process = subprocess.run(
        [binary, "-e", "namelist.def"], cwd=workdir, env=environment,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True)
    write(os.path.join(workdir, "hphi.log"), process.stdout)
    if process.returncode != 0:
        raise AssertionError("HPhi failed:\n{}".format(process.stdout[-8000:]))

    output = os.path.join(workdir, "output")
    # zvo_phys.dat is the FullDiag GC contract and pins row order/inventory,
    # but HPhi formats its energy column with only six digits after the
    # decimal.  Eigenvalue.dat is produced by the same diagonalization with
    # ten digits, so use it for the requested 1e-8 numerical gate after
    # confirming that both files describe the same ordered rows.
    hphi_phys_spectrum = parse_hphi_spectrum(
        os.path.join(output, "zvo_phys.dat"))
    hphi_spectrum = parse_hphi_eigenvalues(
        os.path.join(output, "Eigenvalue.dat"))
    if len(hphi_phys_spectrum) != len(hphi_spectrum):
        raise AssertionError("HPhi physical/eigenvalue inventories differ")
    for physical, precise in zip(hphi_phys_spectrum, hphi_spectrum):
        if abs(physical - precise) > 5.1e-7:
            raise AssertionError("HPhi physical/eigenvalue row order differs")
    match_spectrum(even_spectrum, hphi_spectrum)
    hphi_anomalous = parse_hphi_anomalous(os.path.join(
        output, "zvo_AnomalousG_eigen0.dat"))
    expected_ground = ed_ground_state(mu, delta)
    expected_keys = set(
        (operator_type, first, second)
        for operator_type, first, second, unused_operators
        in ANOMALOUS_G_OPERATORS
    )
    if set(expected_ground["anomalous_g"]) != expected_keys:
        raise AssertionError("ED AnomalousG oracle inventory changed")
    assert_complex_maps_close(expected_ground["anomalous_g"], hphi_anomalous)
    if abs(hphi_spectrum[0] - expected_ground["energy"]) > TOLERANCE:
        raise AssertionError("HPhi and ED ground-state energies differ")
    print("HPhi cross-check passed: full={} even_matched={} E0={:.12g}"
          .format(len(hphi_spectrum), len(even_spectrum), hphi_spectrum[0]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    arguments = parser.parse_args()
    self_test()
    if arguments.self_test:
        print("HPhi cross-check negative self-tests passed")
        return 0

    binary = os.environ.get("MVMC_HPHI_BINARY")
    if not binary:
        print("SKIP: MVMC_HPHI_BINARY is not set")
        return 77
    binary = os.path.realpath(binary)
    if not os.path.isfile(binary) or not os.access(binary, os.X_OK):
        raise AssertionError("MVMC_HPHI_BINARY is not executable: {}".format(binary))
    commit = hphi_commit(binary)
    print("HPhi binary: {}".format(binary))
    print("HPhi commit: {} (reference {})".format(
        commit, HPHI_REFERENCE_COMMIT))
    if commit != "unknown" and not HPHI_REFERENCE_COMMIT.startswith(commit) \
            and not commit.startswith(HPHI_REFERENCE_COMMIT):
        raise AssertionError("HPhi binary is not from the pinned reference commit")
    workdir = os.environ.get(
        "MVMC_HPHI_WORKDIR",
        os.path.join(os.getcwd(), "work", "GC_HPhi_CrossCheck"))
    run_cross_check(binary, os.path.realpath(workdir))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print("ERROR: {}".format(error))
        sys.exit(1)
