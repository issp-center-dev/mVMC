from __future__ import print_function

import os
import shutil
import sys

import numpy as np

from runtest_bf_fsz import run_vmc, update_modpara
from runtest_lanczos2_real_ed import (
    assert_scaled_close,
    check_solver_output,
    DIAGONAL_TRANSFER,
    occupation_state,
    one_hop_states,
    RING_HOPS,
    read_ls2_oracle,
    GUTZWILLER_PARAMETER,
    write_gutzwiller,
    write_lanczos2_namelist,
)
from runtest_lanczos_sector_ed import (
    COULOMB_INTRA,
    NSITE,
    coulomb_intra,
    doublon_count,
    fixed_basis,
    one_body,
    operator_matrix_complex,
    write_coulomb_intra,
)

PEIERLS_PHASE = 0.37
FORWARD_RING_HOPS = frozenset(RING_HOPS[::2])


def hopping_amplitude(left, right):
    phase = PEIERLS_PHASE if (left, right) in FORWARD_RING_HOPS \
        else -PEIERLS_PHASE
    return np.exp(1.0j * phase)


def write_complex_transfer(workdir):
    rows = []
    for left, right in RING_HOPS:
        amplitude = hopping_amplitude(left, right)
        for spin in (0, 1):
            rows.append(
                (left, spin, right, spin, amplitude.real, amplitude.imag)
            )
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


def complex_orbitals():
    orbital = np.array(
        [
            [
                np.sin(0.37 * (row + 1) * (column + 2))
                + 0.3 * np.cos(0.19 * (row + 2) * (column + 1))
                + 0.2j * np.cos(0.23 * (row + 2) * (column + 1))
                for column in range(NSITE)
            ]
            for row in range(NSITE)
        ],
        dtype=complex,
    )
    orbital[1, 1] = orbital[1, 0] * orbital[0, 1] / orbital[0, 0]
    return orbital.reshape(NSITE * NSITE)


def complex_ap_orbital_amplitude(state, orbital):
    up_sites = [
        site for site in range(NSITE) if (state >> site) & 1
    ]
    down_sites = [
        site for site in range(NSITE) if (state >> (site + NSITE)) & 1
    ]
    matrix = np.array(
        [[orbital[up, down] for down in down_sites] for up in up_sites],
        dtype=complex,
    )
    return np.linalg.det(matrix)


def write_complex_init(workdir):
    values = [0.0] * 6
    values.extend([GUTZWILLER_PARAMETER, 0.0, 0.0])
    for value in complex_orbitals():
        values.extend([value.real, value.imag, 0.0])
    filename = "lanczos2_complex_init.dat"
    with open(os.path.join(workdir, filename), "w") as destination:
        destination.write(
            " ".join("{:.18e}".format(value) for value in values)
        )
        destination.write("\n")
    return filename


def build_reference():
    basis = fixed_basis(2, 2)
    orbital = complex_orbitals().reshape((NSITE, NSITE))
    psi = np.array(
        [
            complex_ap_orbital_amplitude(state, orbital)
            * np.exp(GUTZWILLER_PARAMETER * doublon_count(state))
            for state in basis
        ],
        dtype=complex,
    )
    terms = []
    for left, right in RING_HOPS:
        for spin in (0, 1):
            terms.append(
                (-hopping_amplitude(left, right),
                 one_body(left, right, spin))
            )
    site, spin, _, coupling = DIAGONAL_TRANSFER
    terms.append((-coupling, one_body(site, site, spin)))
    terms.extend(
        (coupling, coulomb_intra(site))
        for site, coupling in enumerate(COULOMB_INTRA)
    )
    hamiltonian = operator_matrix_complex(basis, terms)
    hermiticity_residual = np.max(
        np.abs(hamiltonian - np.conj(hamiltonian.T))
    )
    if hermiticity_residual > 2.0e-13:
        raise AssertionError(
            "complex-hopping ED Hamiltonian is not Hermitian: {}".format(
                hermiticity_residual
            )
        )
    powers = [psi]
    for _ in range(3):
        powers.append(np.dot(hamiltonian, powers[-1]))
    reference = {}
    unconjugated_reference = {}
    zero_overlap_states = set()
    for index, state in enumerate(basis):
        if abs(psi[index]) <= 1.0e-13:
            zero_overlap_states.add(state)
            continue
        # mVMC evaluates the bra-oriented local estimator.
        reference[state] = np.array(
            [
                np.conj(powers[order][index] / psi[index])
                for order in range(4)
            ],
            dtype=complex,
        )
        unconjugated_reference[state] = np.array(
            [
                powers[order][index] / psi[index]
                for order in range(4)
            ],
            dtype=complex,
        )
    return reference, unconjugated_reference, zero_overlap_states


def main():
    rootdir = os.getcwd()
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    reference_dir = os.path.join(
        rootdir, "data", "BackFlow_Identity_Complex"
    )
    suffix = "_mpi" if mpi_procs else ""
    workdir = os.path.join(rootdir, "work", "Lanczos2_Complex_ED" + suffix)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    for filename in os.listdir(reference_dir):
        source = os.path.join(reference_dir, filename)
        if os.path.isfile(source):
            shutil.copy2(source, os.path.join(workdir, filename))

    write_coulomb_intra(workdir)
    write_lanczos2_namelist(workdir)
    write_complex_transfer(workdir)
    write_gutzwiller(workdir)
    update_modpara(
        workdir,
        {
            "NLanczosMode": "1",
            "NLanczosEstimatorMode": "1",
            "NLanczosStep": "2",
            "NVMCCalMode": "1",
            "NVMCSample": "256",
            "NVMCWarmUp": "32",
            "2Sz": "0",
            "NMPTrans": "1",
        },
    )
    init_path = write_complex_init(workdir)
    dump_name = "lanczos2_complex_oracle.dat"
    proc = run_vmc(
        rootdir,
        workdir,
        mpi_procs=mpi_procs,
        init_path=init_path,
        extra_env={"MVMC_LANCZOS_ORACLE_DUMP": dump_name},
    )
    if proc.returncode != 0:
        raise AssertionError(
            "Lanczos2 complex vmc.out failed:\n{}".format(proc.stdout)
        )

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
    reference, unconjugated_reference, zero_overlap_states = build_reference()
    if not zero_overlap_states:
        raise AssertionError("complex ED fixture has no zero-overlap state")
    seen = {}
    saw_imaginary = False
    exercised_zero_overlap = False
    conjugation_sensitive = False
    for sample, occupation, power in rows:
        state = occupation_state(occupation)
        if state not in reference:
            raise AssertionError(
                "sample {} has zero or missing ED amplitude".format(sample)
            )
        for order in range(4):
            assert_scaled_close(
                "complex sample {} F{}".format(sample, order),
                power[order],
                reference[state][order],
            )
            if abs(power[order].imag) > 1.0e-8:
                saw_imaginary = True
            if abs(
                reference[state][order]
                - unconjugated_reference[state][order]
            ) > 1.0e-6:
                conjugation_sensitive = True
        if occupation in seen:
            if not np.allclose(power, seen[occupation], rtol=2.0e-12,
                               atol=2.0e-12):
                raise AssertionError(
                    "complex local power changed for repeated occupation"
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
            "Lanczos2 complex ED fixture sampled only {} occupations".format(
                len(seen)
            )
        )
    if not saw_imaginary:
        raise AssertionError("complex Lanczos2 fixture did not exercise Im(F)")
    if not conjugation_sensitive:
        raise AssertionError(
            "complex-hopping fixture is insensitive to the local-power "
            "conjugation convention"
        )
    if not exercised_zero_overlap:
        raise AssertionError(
            "complex ED fixture did not exercise a zero-overlap intermediate"
        )
    hankel_standard_error = check_solver_output(workdir, rows)
    print(
        "Lanczos2 complex ED: PASS "
        "({} samples, {} occupations, Peierls phase={}, "
        "max Hankel SE={:.3f})".format(
            len(rows), len(seen), PEIERLS_PHASE, hankel_standard_error
        )
    )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AssertionError as error:
        print("ERROR: {}".format(error))
        sys.exit(1)
