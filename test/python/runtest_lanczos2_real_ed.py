from __future__ import print_function

import os
import shutil
import sys

import numpy as np

from runtest_bf_fsz import (
    run_vmc,
    update_modpara,
)
from runtest_lanczos_sector_ed import (
    COULOMB_INTRA,
    NSITE,
    ap_orbital_amplitude,
    coulomb_intra,
    doublon_count,
    fixed_basis,
    one_body,
    operator_matrix,
    write_coulomb_intra,
)

RING_HOPS = (
    (1, 0), (0, 1), (2, 1), (1, 2),
    (3, 2), (2, 3), (0, 3), (3, 0),
)
DIAGONAL_TRANSFER = (0, 0, 0, 0.17)
GUTZWILLER_PARAMETER = -0.43
# This is a statistical sanity gate rather than a precision target.  Keep it
# below unity; revise it only with Linux/CI evidence, not a single local run.
HANKEL_RESIDUAL_LIMIT = 0.5
# The stochastic Hankel identities are checked with contiguous-block standard
# errors, so the gate scales with the observed Monte Carlo uncertainty instead
# of relying only on the absolute residual above.  Eight standard errors is a
# conservative CI-stable outlier threshold across the ten matrix comparisons.
HANKEL_STANDARD_ERROR_LIMIT = 8.0


def real_orbitals():
    orbital = np.array(
        [
            [
                np.sin(0.37 * (row + 1) * (column + 2))
                + 0.3 * np.cos(0.19 * (row + 2) * (column + 1))
                for column in range(NSITE)
            ]
            for row in range(NSITE)
        ],
        dtype=float,
    )
    orbital[1, 1] = orbital[1, 0] * orbital[0, 1] / orbital[0, 0]
    return orbital


def write_real_init(workdir):
    values = [0.0] * 6
    values.extend([GUTZWILLER_PARAMETER, 0.0, 0.0])
    for value in real_orbitals().reshape(NSITE * NSITE):
        values.extend([value, 0.0, 0.0])
    filename = "lanczos2_real_init.dat"
    with open(os.path.join(workdir, filename), "w") as destination:
        destination.write(
            " ".join("{:.18e}".format(value) for value in values)
        )
        destination.write("\n")
    return filename


def write_transfer_with_diagonal(workdir):
    rows = []
    for left, right in RING_HOPS:
        for spin in (0, 1):
            rows.append((left, spin, right, spin, 1.0, 0.0))
    site, spin, _, coupling = DIAGONAL_TRANSFER
    rows.append((site, spin, site, spin, coupling, 0.0))
    with open(os.path.join(workdir, "trans.def"), "w") as destination:
        destination.write("========================\n")
        destination.write("NTransfer      {}\n".format(len(rows)))
        destination.write("========================\n")
        destination.write("========i_j_s_tijs======\n")
        destination.write("========================\n")
        for row in rows:
            destination.write(
                "{:5d} {:5d} {:5d} {:5d} {:25.15f} {:25.15f}\n".format(
                    *row
                )
            )


def write_gutzwiller(workdir):
    with open(os.path.join(workdir, "gutzwilleridx.def"), "w") as destination:
        destination.write("=============================================\n")
        destination.write("NGutzwillerIdx          1\n")
        destination.write("ComplexType             0\n")
        destination.write("=============================================\n")
        destination.write("=============================================\n")
        for site in range(NSITE):
            destination.write("{:5d} {:5d}\n".format(site, 0))
        destination.write("{:5d} {:5d}\n".format(0, 0))


def read_ls2_oracle(path):
    rows = []
    with open(path) as source:
        for line in source:
            columns = line.split()
            if len(columns) != 4 + 2 * NSITE + 8:
                raise AssertionError(
                    "invalid Lanczos2 oracle width: {}".format(line.rstrip())
                )
            if (
                columns[0] != "sample"
                or columns[2] != "occ"
                or columns[3 + 2 * NSITE] != "ls2power"
            ):
                raise AssertionError("invalid Lanczos2 oracle markers")
            occupation = tuple(
                int(value) for value in columns[3 : 3 + 2 * NSITE]
            )
            values = [
                float(value) for value in columns[4 + 2 * NSITE :]
            ]
            power = np.array(
                [
                    complex(values[2 * index], values[2 * index + 1])
                    for index in range(4)
                ],
                dtype=complex,
            )
            rows.append((int(columns[1]), occupation, power))
    if not rows:
        raise AssertionError("Lanczos2 oracle dump is empty")
    return rows


def write_lanczos2_namelist(workdir):
    entries = [
        ("ModPara", "modpara.def"),
        ("LocSpin", "locspn.def"),
        ("Trans", "trans.def"),
        ("CoulombIntra", "coulombintra.def"),
        ("Gutzwiller", "gutzwilleridx.def"),
        ("Orbital", "orbitalidx.def"),
        ("TransSym", "qptransidx.def"),
    ]
    with open(os.path.join(workdir, "namelist.def"), "w") as destination:
        for keyword, filename in entries:
            destination.write("{:>16}  {}\n".format(keyword, filename))


def occupation_state(occupation):
    state = 0
    for orbital, value in enumerate(occupation):
        if value:
            state |= 1 << orbital
    return state


def one_hop_states(state):
    for left, right in RING_HOPS:
        for spin in (0, 1):
            left_orbital = left + spin * NSITE
            right_orbital = right + spin * NSITE
            if ((state >> right_orbital) & 1) and not (
                (state >> left_orbital) & 1
            ):
                yield state ^ (1 << left_orbital) ^ (1 << right_orbital)


def read_numeric_record(path):
    numeric = []
    with open(path) as source:
        lines = source.readlines()
    if len(lines) != 2 or not lines[0].startswith("# mVMC "):
        raise AssertionError("{} does not have the v1 two-line layout".format(path))
    if not lines[1].endswith("\n"):
        raise AssertionError("{} data record lacks final newline".format(path))
    numeric.extend(float(value) for value in lines[1].split())
    return numeric


def read_shifted_moment(path):
    with open(path) as source:
        lines = source.readlines()
    prefix = "# mVMC ls2_moment v2 basis_shift="
    if len(lines) != 2 or not lines[0].startswith(prefix):
        raise AssertionError(
            "{} does not have the v2 shifted-moment layout".format(path)
        )
    shift = float(lines[0][len(prefix):].split()[0])
    values = [float(value) for value in lines[1].split()]
    if len(values) != 32:
        raise AssertionError(
            "ls2_moment column count is {}".format(len(values))
        )
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


def stochastic_hankel_standard_error(powers, shift):
    shifted = np.array(
        [shift_local_power(power, shift) for power in powers],
        dtype=complex,
    )
    if len(shifted) < 32:
        raise AssertionError(
            "stochastic Hankel diagnostic needs at least 32 samples"
        )

    block_count = min(16, len(shifted) // 8)
    maximum_score = 0.0
    comparison_count = 0
    epsilon = np.finfo(float).eps
    for order in range(1, 6):
        pairs = [
            (left, order - left)
            for left in range(4)
            if 0 <= order - left < 4
        ]
        reference_pair = pairs[0]
        reference = (
            np.conj(shifted[:, reference_pair[0]])
            * shifted[:, reference_pair[1]]
        )
        for pair in pairs[1:]:
            estimate = np.conj(shifted[:, pair[0]]) * shifted[:, pair[1]]
            difference = estimate - reference
            sample_scale = max(
                np.max(np.abs(estimate)),
                np.max(np.abs(reference)),
                1.0,
            )
            if np.max(np.abs(difference)) <= 128.0 * epsilon * sample_scale:
                continue

            block_means = np.array(
                [
                    np.mean(block)
                    for block in np.array_split(difference, block_count)
                ],
                dtype=complex,
            )
            standard_error = np.sqrt(
                np.sum(
                    np.abs(block_means - np.mean(block_means)) ** 2
                )
                / float(block_count * (block_count - 1))
            )
            roundoff_floor = 128.0 * epsilon * max(
                np.mean(np.abs(estimate)),
                np.mean(np.abs(reference)),
                1.0,
            )
            score = abs(np.mean(difference)) / max(
                standard_error, roundoff_floor
            )
            if not np.isfinite(score):
                raise AssertionError(
                    "stochastic Hankel standard-error score is not finite"
                )
            maximum_score = max(maximum_score, score)
            comparison_count += 1

    if comparison_count == 0:
        raise AssertionError(
            "stochastic Hankel diagnostic found no nontrivial comparison"
        )
    if maximum_score >= HANKEL_STANDARD_ERROR_LIMIT:
        raise AssertionError(
            "stochastic Hankel identity exceeds its block standard error: "
            "max score={} (limit={})".format(
                maximum_score, HANKEL_STANDARD_ERROR_LIMIT
            )
        )
    return maximum_score


def build_reference():
    basis = fixed_basis(2, 2)
    orbital = real_orbitals()
    unprojected_psi = np.array(
        [ap_orbital_amplitude(state, orbital) for state in basis],
        dtype=float,
    )
    psi = np.array(
        [
            amplitude * np.exp(
                GUTZWILLER_PARAMETER * doublon_count(state)
            )
            for state, amplitude in zip(basis, unprojected_psi)
        ],
        dtype=float,
    )
    terms = []
    for left, right in RING_HOPS:
        for spin in (0, 1):
            terms.append((-1.0, one_body(left, right, spin)))
    site, spin, _, coupling = DIAGONAL_TRANSFER
    terms.append((-coupling, one_body(site, site, spin)))
    terms.extend(
        (coupling, coulomb_intra(site))
        for site, coupling in enumerate(COULOMB_INTRA)
    )
    hamiltonian = operator_matrix(basis, terms)
    powers = [psi]
    for _ in range(3):
        powers.append(np.dot(hamiltonian, powers[-1]))
    unprojected_powers = [unprojected_psi]
    for _ in range(3):
        unprojected_powers.append(
            np.dot(hamiltonian, unprojected_powers[-1])
        )
    projection_sensitivity = 0.0
    for index in range(len(basis)):
        if abs(psi[index]) <= 1.0e-13 or \
                abs(unprojected_psi[index]) <= 1.0e-13:
            continue
        projection_sensitivity = max(
            projection_sensitivity,
            abs(
                powers[3][index] / psi[index]
                - unprojected_powers[3][index] / unprojected_psi[index]
            ),
        )
    if projection_sensitivity <= 1.0e-6:
        raise AssertionError(
            "real ED fixture does not exercise ProjRatio in F3"
        )
    validate_full_enumeration(powers)
    reference = {}
    zero_overlap_states = set()
    for index, state in enumerate(basis):
        if abs(psi[index]) <= 1.0e-13:
            zero_overlap_states.add(state)
            continue
        reference[state] = np.array(
            [powers[order][index] / psi[index] for order in range(4)],
            dtype=complex,
        )
    return reference, zero_overlap_states


def validate_full_enumeration(powers):
    normalization = np.vdot(powers[0], powers[0])
    moment = np.array(
        [
            [
                np.vdot(powers[left], powers[right]) / normalization
                for right in range(4)
            ]
            for left in range(4)
        ],
        dtype=complex,
    )
    scale = max(np.max(np.abs(moment)), 1.0)
    for left in range(4):
        for right in range(4):
            for other_left in range(4):
                other_right = left + right - other_left
                if 0 <= other_right < 4:
                    if (
                        abs(
                            moment[left, right]
                            - moment[other_left, other_right]
                        )
                        > 2.0e-12 * scale
                    ):
                        raise AssertionError(
                            "full-enumeration Hankel identity failed"
                        )

    energies = []
    for dimension in (1, 2, 3):
        overlap = moment[:dimension, :dimension]
        reduced_hamiltonian = moment[:dimension, 1 : dimension + 1]
        eigenvalues = np.linalg.eigvals(
            np.linalg.solve(overlap, reduced_hamiltonian)
        )
        if np.max(np.abs(np.imag(eigenvalues))) > 2.0e-11:
            raise AssertionError(
                "full-enumeration Krylov eigenvalue is not real"
            )
        energies.append(np.min(np.real(eigenvalues)))
    tolerance = 2.0e-11 * max(max(abs(value) for value in energies), 1.0)
    if energies[2] > energies[1] + tolerance:
        raise AssertionError("E2 is above E1 in full enumeration")
    if energies[1] > energies[0] + tolerance:
        raise AssertionError("E1 is above Evar in full enumeration")


def assert_scaled_close(label, actual, expected, absolute=2.0e-10,
                        relative=2.0e-9):
    difference = abs(actual - expected)
    allowed = absolute + relative * max(abs(actual), abs(expected))
    if difference > allowed:
        raise AssertionError(
            "{} mismatch: actual={} expected={} difference={} allowed={}".format(
                label, actual, expected, difference, allowed
            )
        )


def check_solver_output(workdir, rows):
    output_path = os.path.join(workdir, "output", "zvo_ls2_out_001.dat")
    moment_path = os.path.join(workdir, "output", "zvo_ls2_moment_001.dat")
    output = read_numeric_record(output_path)
    basis_shift, observed_moment = read_shifted_moment(moment_path)
    if len(output) != 17:
        raise AssertionError("ls2_out column count is {}".format(len(output)))
    if abs(output[15]) > 1.0e-12:
        raise AssertionError(
            "antihermitian_residual is unexpectedly large: {}".format(
                output[15]
            )
        )
    if not np.isfinite(output[16]) or \
            output[16] >= HANKEL_RESIDUAL_LIMIT:
        raise AssertionError(
            "hankel_residual is not a finite sub-unity diagnostic: "
            "{} (limit={})".format(
                output[16], HANKEL_RESIDUAL_LIMIT
            )
        )
    if os.path.exists(os.path.join(workdir, "output", "zvo_ls_out_001.dat")):
        raise AssertionError("NLanczosStep=2 wrote a legacy ls_out file")
    if os.path.exists(os.path.join(workdir, "output", "zvo_ls_qqqq_001.dat")):
        raise AssertionError("NLanczosStep=2 wrote a legacy ls_qqqq file")

    expected_raw_moment = np.zeros((4, 4), dtype=complex)
    expected_shifted_moment = np.zeros((4, 4), dtype=complex)
    for _, _, power in rows:
        shifted_power = shift_local_power(power, basis_shift)
        expected_raw_moment += np.outer(np.conj(power), power)
        expected_shifted_moment += np.outer(
            np.conj(shifted_power), shifted_power
        )
    expected_raw_moment /= float(len(rows))
    expected_shifted_moment /= float(len(rows))
    if not np.allclose(
        observed_moment, expected_shifted_moment,
        rtol=2.0e-12, atol=2.0e-12
    ):
        difference = np.max(
            np.abs(observed_moment - expected_shifted_moment)
        )
        raise AssertionError(
            "shifted moment does not match centered sample-local powers: "
            "max diff={}".format(difference)
        )
    hankel_standard_error = stochastic_hankel_standard_error(
        [power for _, _, power in rows], basis_shift
    )

    overlap = expected_raw_moment[:3, :3]
    reduced_hamiltonian = 0.5 * (
        expected_raw_moment[:3, 1:4] + expected_raw_moment[1:4, :3]
    )
    hamiltonian_squared = expected_raw_moment[1:4, 1:4]
    coefficient = np.array(
        [
            complex(output[2], output[3]),
            complex(output[4], output[5]),
            complex(output[6], output[7]),
        ],
        dtype=complex,
    )
    norm = np.vdot(coefficient, np.dot(overlap, coefficient)).real
    energy = np.vdot(
        coefficient, np.dot(reduced_hamiltonian, coefficient)
    ).real
    h2 = np.vdot(
        coefficient, np.dot(hamiltonian_squared, coefficient)
    ).real
    variance_ratio = (h2 - energy * energy) / (energy * energy)
    assert_scaled_close("coefficient norm", norm, 1.0, 2.0e-11, 2.0e-11)
    assert_scaled_close("reported energy", output[0], energy, 2.0e-11, 2.0e-11)
    assert_scaled_close(
        "reported variance/E2", output[1], variance_ratio, 2.0e-11, 2.0e-11
    )
    residual = np.linalg.norm(
        np.dot(reduced_hamiltonian, coefficient)
        - output[0] * np.dot(overlap, coefficient)
    )
    scale = max(
        np.linalg.norm(np.dot(reduced_hamiltonian, coefficient)), 1.0
    )
    if residual / scale > 2.0e-10:
        raise AssertionError(
            "reported coefficient GEV residual={}".format(residual / scale)
        )

    timer_path = os.path.join(
        workdir, "output", "zvo_CalcTimer.dat"
    )
    with open(timer_path) as source:
        timer_rows = [line for line in source if "[95]" in line]
    if len(timer_rows) != 1:
        raise AssertionError("Lanczos2 timer [95] is missing or duplicated")
    if float(timer_rows[0].split()[-1]) <= 0.0:
        raise AssertionError("Lanczos2 timer [95] is not positive")
    return hankel_standard_error


def main():
    rootdir = os.getcwd()
    expected_error = sys.argv[1] if len(sys.argv) == 2 else None
    if len(sys.argv) > 2:
        raise AssertionError("usage: runtest_lanczos2_real_ed.py [expected-error]")
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    split_size = os.environ.get("MVMC_NSPLIT_SIZE", "1")
    reference_dir = os.path.join(rootdir, "data", "BackFlow_Identity_Real")
    suffix = "_split{}".format(split_size) if split_size != "1" else ""
    if mpi_procs:
        suffix += "_mpi"
    workdir = os.path.join(rootdir, "work", "Lanczos2_Real_ED" + suffix)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    for filename in os.listdir(reference_dir):
        source = os.path.join(reference_dir, filename)
        if os.path.isfile(source):
            shutil.copy2(source, os.path.join(workdir, filename))

    write_coulomb_intra(workdir)
    write_lanczos2_namelist(workdir)
    write_transfer_with_diagonal(workdir)
    write_gutzwiller(workdir)
    update_modpara(
        workdir,
        {
            "NLanczosMode": "1",
            "NLanczosStep": "2",
            "NVMCCalMode": "1",
            "NVMCSample": "256",
            "NVMCWarmUp": "32",
            "2Sz": "0",
            "NMPTrans": "1",
            "NSplitSize": split_size,
        },
    )
    init_path = write_real_init(workdir)
    dump_name = "lanczos2_real_oracle.dat"
    proc = run_vmc(
        rootdir,
        workdir,
        mpi_procs=mpi_procs,
        init_path=init_path,
        extra_env={"MVMC_LANCZOS_ORACLE_DUMP": dump_name},
    )
    if expected_error is not None:
        if proc.returncode == 0:
            raise AssertionError(
                "Lanczos2 fault injection unexpectedly succeeded"
            )
        if expected_error not in proc.stdout:
            raise AssertionError(
                "expected error {!r} was not found:\n{}".format(
                    expected_error, proc.stdout
                )
            )
        print("Lanczos2 non-finite fault injection: PASS")
        return 0
    if proc.returncode != 0:
        raise AssertionError("Lanczos2 real vmc.out failed:\n{}".format(proc.stdout))

    if mpi_procs:
        rows = []
        for rank in range(int(mpi_procs)):
            rows.extend(
                read_ls2_oracle(
                    os.path.join(
                        workdir, dump_name + ".rank{:04d}".format(rank)
                    )
                )
            )
    else:
        rows = read_ls2_oracle(os.path.join(workdir, dump_name))
    reference, zero_overlap_states = build_reference()
    if not zero_overlap_states:
        raise AssertionError("real ED fixture has no zero-overlap state")
    seen = {}
    exercised_zero_overlap = False
    for sample, occupation, power in rows:
        state = occupation_state(occupation)
        if state not in reference:
            raise AssertionError(
                "sample {} has zero or missing ED amplitude".format(sample)
            )
        for order in range(4):
            assert_scaled_close(
                "sample {} F{}".format(sample, order),
                power[order],
                reference[state][order],
            )
        if occupation in seen:
            if not np.allclose(power, seen[occupation], rtol=2.0e-12,
                               atol=2.0e-12):
                raise AssertionError(
                    "local power changed for repeated occupation"
                )
        else:
            seen[occupation] = power
        if any(
            neighbor in zero_overlap_states
            for neighbor in one_hop_states(state)
        ):
            exercised_zero_overlap = True
    if len(seen) < 8:
        raise AssertionError(
            "Lanczos2 ED fixture sampled only {} occupations".format(len(seen))
        )
    if not exercised_zero_overlap:
        raise AssertionError(
            "real ED fixture did not exercise a zero-overlap intermediate"
        )
    hankel_standard_error = check_solver_output(workdir, rows)
    print(
        "Lanczos2 real ED: PASS "
        "({} samples, {} occupations, max Hankel SE={:.3f})".format(
            len(rows), len(seen), hankel_standard_error
        )
    )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AssertionError as error:
        print("ERROR: {}".format(error))
        sys.exit(1)
