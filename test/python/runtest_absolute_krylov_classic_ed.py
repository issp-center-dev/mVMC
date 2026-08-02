from __future__ import print_function

import os
import shutil
import subprocess
import sys

import numpy as np


SITE_COUNT = 4
ORBITAL_COUNT = 2 * SITE_COUNT
GUTZWILLER_PARAMETER = -0.27
BONDS = ((0, 1), (1, 2), (2, 3), (3, 0))
INTRA = (0.25, -0.125, 0.375, -0.3125)
INTER = ((0, 2, 0.1875), (1, 3, -0.15625))
ELECTRONIC_HUND = ((0, 1, 0.21875), (2, 3, -0.09375))
SPIN_HUND = (0.375, -0.25, 0.1875, -0.125)
SPIN_EXCHANGE = (0.625, -0.3125, 0.4375, -0.5625)


def popcount(value):
    return bin(value).count("1")


def fixed_basis(pure_spin):
    basis = []
    for up_mask in range(1 << SITE_COUNT):
        if popcount(up_mask) != 2:
            continue
        for down_mask in range(1 << SITE_COUNT):
            if popcount(down_mask) != 2:
                continue
            if pure_spin and down_mask != ((~up_mask) & 0x0f):
                continue
            basis.append(up_mask | (down_mask << SITE_COUNT))
    return basis


def apply_operators(state, operators):
    sign = 1
    updated = state
    for kind, orbital in reversed(operators):
        occupied = (updated >> orbital) & 1
        parity = popcount(updated & ((1 << orbital) - 1))
        if kind == "annihilate":
            if not occupied:
                return None
            if parity & 1:
                sign = -sign
            updated &= ~(1 << orbital)
        else:
            if occupied:
                return None
            if parity & 1:
                sign = -sign
            updated |= 1 << orbital
    return updated, sign


def density_term(first, second, coefficient):
    return coefficient, (
        ("create", first), ("annihilate", first),
        ("create", second), ("annihilate", second),
    )


def electronic_terms(complex_case):
    terms = []
    for bond_index, (left, right) in enumerate(BONDS):
        hopping = 0.75
        if complex_case:
            hopping += 1.0j * (0.125 + 0.0625 * bond_index)
        for spin in (0, 1):
            left_orbital = left + spin * SITE_COUNT
            right_orbital = right + spin * SITE_COUNT
            terms.append((
                -hopping,
                (("create", left_orbital),
                 ("annihilate", right_orbital)),
            ))
            terms.append((
                -np.conj(hopping),
                (("create", right_orbital),
                 ("annihilate", left_orbital)),
            ))
    terms.append((
        -0.4375,
        (("create", 1), ("annihilate", 1)),
    ))
    for site, coefficient in enumerate(INTRA):
        terms.append(density_term(
            site, site + SITE_COUNT, coefficient
        ))
    for first, second, coefficient in INTER:
        for first_spin in (0, 1):
            for second_spin in (0, 1):
                terms.append(density_term(
                    first + first_spin * SITE_COUNT,
                    second + second_spin * SITE_COUNT,
                    coefficient,
                ))
    for first, second, coefficient in ELECTRONIC_HUND:
        for spin in (0, 1):
            terms.append(density_term(
                first + spin * SITE_COUNT,
                second + spin * SITE_COUNT,
                -coefficient,
            ))
    return terms


def spin_terms():
    terms = []
    for bond_index, (first, second) in enumerate(BONDS):
        for spin in (0, 1):
            terms.append(density_term(
                first + spin * SITE_COUNT,
                second + spin * SITE_COUNT,
                -SPIN_HUND[bond_index],
            ))
        up_first = first
        up_second = second
        down_first = first + SITE_COUNT
        down_second = second + SITE_COUNT
        terms.append((
            SPIN_EXCHANGE[bond_index],
            (("create", up_first), ("annihilate", up_second),
             ("create", down_second), ("annihilate", down_first)),
        ))
        terms.append((
            SPIN_EXCHANGE[bond_index],
            (("create", up_second), ("annihilate", up_first),
             ("create", down_first), ("annihilate", down_second)),
        ))
    return terms


def dense_hamiltonian(basis, terms):
    state_index = {state: index for index, state in enumerate(basis)}
    matrix = np.zeros((len(basis), len(basis)), dtype=complex)
    for column, state in enumerate(basis):
        for coefficient, operators in terms:
            applied = apply_operators(state, operators)
            if applied is None:
                continue
            destination, sign = applied
            if destination not in state_index:
                raise AssertionError(
                    "term escaped the fixed sector: {} -> {}".format(
                        state, destination
                    )
                )
            matrix[state_index[destination], column] += sign * coefficient
    residual = np.max(np.abs(matrix - np.conj(matrix.T)))
    if residual > 2.0e-13:
        raise AssertionError(
            "independent Hamiltonian is not Hermitian: {:.3e}".format(
                residual
            )
        )
    return matrix


def pair_matrix(complex_case):
    if complex_case:
        pairs = np.array([
            [2.0 + 1.0j, 1.0 + 1.5j, 0.5 - 0.5j, 0.5 - 1.5j],
            [-2.0 - 1.0j, -1.0 - 1.5j, -0.5 + 0.5j, -0.5 + 1.5j],
            [0.5 + 1.0j, 0.5 + 2.0j, -0.5 - 2.0j, -0.5 - 1.5j],
            [1.5 + 1.5j, 0.5 + 0.5j, 1.5 - 1.5j, 1.5 + 1.5j],
        ], dtype=complex)
        return pairs
    return np.array([
        [1.5, -0.5, 1.5, 1.0],
        [-1.5, 0.5, -1.5, -1.0],
        [-0.5, 1.5, 1.5, 2.0],
        [-2.0, -2.0, -1.0, 0.5],
    ], dtype=complex)


def pfaffian(matrix):
    size = matrix.shape[0]
    if size == 0:
        return 1.0 + 0.0j
    value = 0.0 + 0.0j
    for column in range(1, size):
        keep = [index for index in range(1, size) if index != column]
        minor = matrix[np.ix_(keep, keep)]
        sign = 1.0 if column % 2 == 1 else -1.0
        value += sign * matrix[0, column] * pfaffian(minor)
    return value


def amplitude(state, complex_case, pure_spin):
    up_sites = [
        site for site in range(SITE_COUNT) if (state >> site) & 1
    ]
    down_sites = [
        site for site in range(SITE_COUNT)
        if (state >> (site + SITE_COUNT)) & 1
    ]
    pairs = pair_matrix(complex_case)
    matrix = np.zeros((4, 4), dtype=complex)
    for up_index, up_site in enumerate(up_sites):
        for down_index, down_site in enumerate(down_sites):
            value = pairs[up_site, down_site]
            matrix[up_index, 2 + down_index] = value
            matrix[2 + down_index, up_index] = -value
    value = pfaffian(matrix)
    if not pure_spin:
        doublons = sum(
            ((state >> site) & 1)
            and ((state >> (site + SITE_COUNT)) & 1)
            for site in range(SITE_COUNT)
        )
        value *= np.exp(GUTZWILLER_PARAMETER * doublons)
    return value


def independent_reference(name):
    pure_spin = name.startswith("spin_")
    complex_case = name.endswith("_complex")
    basis = fixed_basis(pure_spin)
    terms = spin_terms() if pure_spin else electronic_terms(complex_case)
    hamiltonian = dense_hamiltonian(basis, terms)
    powers = [np.array([
        amplitude(state, complex_case, pure_spin) for state in basis
    ], dtype=complex)]
    for _ in range(3):
        powers.append(np.dot(hamiltonian, powers[-1]))
    zero_indices = [
        index for index, value in enumerate(powers[0])
        if abs(value) <= 2.0e-14
    ]
    if not zero_indices:
        raise AssertionError("{} has no exact zero-support state".format(name))
    if not any(
            max(abs(powers[order][index]) for order in range(1, 4))
            > 1.0e-8 for index in zero_indices):
        raise AssertionError(
            "{} has no zero-support bridge into H^k psi".format(name)
        )
    return {
        state: np.array([powers[order][index] for order in range(4)])
        for index, state in enumerate(basis)
    }


def parse_driver_output(output):
    cases = {}
    metadata = {}
    active = None
    for line in output.splitlines():
        columns = line.split()
        if not columns:
            continue
        if columns[0] == "BEGIN":
            active = columns[1]
            cases[active] = {}
            metadata[active] = {
                "rank_count": int(columns[3]),
                "empty_qp_ranks": int(columns[5]),
            }
        elif columns[0] == "ROW":
            if active is None or len(columns) != 10:
                raise AssertionError("malformed ROW: {}".format(line))
            values = np.array([
                float(columns[2 + 2 * order])
                + 1.0j * float(columns[3 + 2 * order])
                for order in range(4)
            ])
            assert_finite_power("driver ROW", values)
            cases[active][int(columns[1])] = values
        elif columns[0] == "END":
            if active != columns[1]:
                raise AssertionError("mismatched END: {}".format(line))
            active = None
    if active is not None:
        raise AssertionError("driver output ended inside {}".format(active))
    return cases, metadata


def assert_finite_power(label, values):
    if not np.all(np.isfinite(np.asarray(values, dtype=complex))):
        raise AssertionError("{} contains NaN or Inf: {}".format(label, values))


def check_nonfinite_guard():
    for value in (np.nan + 0.0j, np.inf + 0.0j, 0.0 + np.inf * 1.0j):
        try:
            assert_finite_power("negative-control", np.array([value]))
        except AssertionError:
            continue
        raise AssertionError("nonfinite comparison guard accepted {}".format(value))


def driver_command(driver):
    mpi_procs = int(os.environ.get("MVMC_MPI_PROCS", "1"))
    if mpi_procs == 1:
        return [driver], mpi_procs
    mpiexec = os.environ.get("MVMC_MPIEXEC") or shutil.which("mpiexec")
    if not mpiexec:
        raise AssertionError("MVMC_MPIEXEC/mpiexec is unavailable")
    numproc_flag = os.environ.get("MVMC_MPIEXEC_NUMPROC_FLAG", "-np")
    return [mpiexec, numproc_flag, str(mpi_procs), driver], mpi_procs


def main():
    if len(sys.argv) != 2:
        raise AssertionError("usage: runtest_absolute_krylov_classic_ed.py DRIVER")
    check_nonfinite_guard()
    command, expected_ranks = driver_command(sys.argv[1])
    process = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True, check=False,
    )
    if process.returncode != 0:
        raise AssertionError(
            "absolute Krylov classic driver failed:\n{}".format(
                process.stdout
            )
        )
    cases, metadata = parse_driver_output(process.stdout)
    expected_cases = (
        "electronic_real", "electronic_complex", "spin_real", "spin_complex"
    )
    if set(cases) != set(expected_cases):
        raise AssertionError("driver cases mismatch: {}".format(sorted(cases)))
    if expected_ranks > 1:
        serial = subprocess.run(
            [sys.argv[1]], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            universal_newlines=True, check=False,
        )
        if serial.returncode != 0:
            raise AssertionError(
                "serial rank-invariance reference failed:\n{}".format(
                    serial.stdout
                )
            )
        serial_cases, _ = parse_driver_output(serial.stdout)
        for name in expected_cases:
            if set(serial_cases[name]) != set(cases[name]):
                raise AssertionError("{} MPI/serial basis drift".format(name))
            for state in cases[name]:
                if not np.array_equal(
                        cases[name][state], serial_cases[name][state]):
                    raise AssertionError(
                        "{} state={} is not bitwise rank invariant".format(
                            name, state
                        )
                    )
    for name in expected_cases:
        reference = independent_reference(name)
        if set(cases[name]) != set(reference):
            raise AssertionError("{} basis mismatch".format(name))
        if metadata[name]["rank_count"] != expected_ranks:
            raise AssertionError("{} rank-count mismatch".format(name))
        expected_empty = max(0, expected_ranks - 1)
        if metadata[name]["empty_qp_ranks"] != expected_empty:
            raise AssertionError("{} empty-QP audit mismatch".format(name))
        for state, expected in reference.items():
            actual = cases[name][state]
            assert_finite_power("{} state={} expected".format(name, state),
                                expected)
            error = np.max(np.abs(actual - expected))
            scale = max(1.0, float(np.max(np.abs(expected))))
            absolute_tolerance = 2.0e-13 if name.endswith("_real") else 2.0e-12
            relative_tolerance = 2.0e-12 if name.endswith("_real") else 2.0e-11
            allowed = absolute_tolerance + relative_tolerance * scale
            if error > allowed:
                raise AssertionError(
                    "{} state={} mismatch: error={:.3e} allowed={:.3e}\n"
                    "actual={}\nexpected={}".format(
                        name, state, error, allowed, actual, expected
                    )
                )
        if name.endswith("_real"):
            imaginary = max(
                float(np.max(np.abs(values.imag)))
                for values in cases[name].values()
            )
            if imaginary > 2.0e-12:
                raise AssertionError(
                    "{} unexpectedly complex: {:.3e}".format(name, imaginary)
                )
    print(
        "absolute Krylov classic ED passed: 4 families, ranks={}, "
        "orders=0..3, zero-support bridges verified".format(expected_ranks)
    )


if __name__ == "__main__":
    main()
