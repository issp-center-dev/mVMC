from __future__ import print_function

import os
import shutil
import sys

import numpy as np

from runtest_bf_fsz import run_vmc, update_modpara
from runtest_lanczos2_real_ed import (
    DIAGONAL_TRANSFER,
    GUTZWILLER_PARAMETER,
    RING_HOPS,
    assert_scaled_close,
    check_solver_output,
    occupation_state,
    read_ls2_oracle,
    real_orbitals,
    write_gutzwiller,
    write_lanczos2_namelist,
    write_real_init,
    write_transfer_with_diagonal,
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


TRANSLATIONS = tuple(
    tuple((site + displacement) % NSITE for site in range(NSITE))
    for displacement in range(NSITE)
)
MOMENTUM_PI_WEIGHTS = (1.0, -1.0, 1.0, -1.0)


def write_translation_projection(workdir):
    with open(os.path.join(workdir, "qptransidx.def"), "w") as destination:
        destination.write("=============================================\n")
        destination.write("NQPTrans          {}\n".format(len(TRANSLATIONS)))
        destination.write("=============================================\n")
        destination.write(
            "======== TrIdx_TrWeight_and_TrIdx_i_xi ======\n"
        )
        destination.write("=============================================\n")
        for index, weight in enumerate(MOMENTUM_PI_WEIGHTS):
            destination.write(
                "{:5d} {:25.15f} {:25.15f}\n".format(index, weight, 0.0)
            )
        for index, translation in enumerate(TRANSLATIONS):
            for origin, translated in enumerate(translation):
                destination.write(
                    "{:5d} {:5d} {:5d} {:5d}\n".format(
                        index, origin, translated, 1
                    )
                )


def translated_ap_amplitude(state, orbital, translation):
    up_sites = [
        site for site in range(NSITE) if (state >> site) & 1
    ]
    down_sites = [
        site for site in range(NSITE)
        if (state >> (site + NSITE)) & 1
    ]
    matrix = np.array(
        [
            [
                orbital[translation[up], translation[down]]
                for down in down_sites
            ]
            for up in up_sites
        ],
        dtype=float,
    )
    return np.linalg.det(matrix)


def projected_amplitude(state, orbital):
    return sum(
        weight * translated_ap_amplitude(state, orbital, translation)
        for weight, translation in zip(
            MOMENTUM_PI_WEIGHTS, TRANSLATIONS
        )
    )


def build_projected_reference():
    basis = fixed_basis(2, 2)
    orbital = real_orbitals()
    unprojected = np.array(
        [
            ap_orbital_amplitude(state, orbital)
            * np.exp(GUTZWILLER_PARAMETER * doublon_count(state))
            for state in basis
        ],
        dtype=float,
    )
    projected = np.array(
        [
            projected_amplitude(state, orbital)
            * np.exp(GUTZWILLER_PARAMETER * doublon_count(state))
            for state in basis
        ],
        dtype=float,
    )
    if np.min(np.abs(projected)) <= 1.0e-13:
        raise AssertionError(
            "projected ED fixture unexpectedly has a zero-overlap state"
        )

    fidelity = (
        abs(np.vdot(unprojected, projected)) ** 2
        / (
            np.vdot(unprojected, unprojected).real
            * np.vdot(projected, projected).real
        )
    )
    if fidelity >= 0.95:
        raise AssertionError(
            "momentum projection is too close to the identity: "
            "fidelity={}".format(fidelity)
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
    powers = [projected]
    for _ in range(3):
        powers.append(np.dot(hamiltonian, powers[-1]))

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
                if 0 <= other_right < 4 and abs(
                    moment[left, right]
                    - moment[other_left, other_right]
                ) > 2.0e-12 * scale:
                    raise AssertionError(
                        "projected full-enumeration Hankel identity failed"
                    )

    reference = {}
    for index, state in enumerate(basis):
        reference[state] = np.array(
            [
                powers[order][index] / projected[index]
                for order in range(4)
            ],
            dtype=complex,
        )
    return reference, fidelity


def main():
    rootdir = os.getcwd()
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    reference_dir = os.path.join(rootdir, "data", "BackFlow_Identity_Real")
    suffix = "_mpi" if mpi_procs else ""
    workdir = os.path.join(rootdir, "work", "Lanczos2_Projected_ED" + suffix)
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
    write_translation_projection(workdir)
    update_modpara(
        workdir,
        {
            "NLanczosMode": "1",
            "NLanczosStep": "2",
            "NVMCCalMode": "1",
            "NVMCSample": "256",
            "NVMCWarmUp": "32",
            "2Sz": "0",
            "NSPGaussLeg": "1",
            "NSPStot": "0",
            "NMPTrans": str(len(TRANSLATIONS)),
        },
    )
    init_path = write_real_init(workdir)
    dump_name = "lanczos2_projected_oracle.dat"
    process = run_vmc(
        rootdir,
        workdir,
        mpi_procs=mpi_procs,
        init_path=init_path,
        extra_env={"MVMC_LANCZOS_ORACLE_DUMP": dump_name},
    )
    if process.returncode != 0:
        raise AssertionError(
            "Lanczos2 projected vmc.out failed:\n{}".format(process.stdout)
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

    reference, identity_fidelity = build_projected_reference()
    seen = {}
    for sample, occupation, power in rows:
        state = occupation_state(occupation)
        if state not in reference:
            raise AssertionError(
                "sample {} is outside the projected ED basis".format(sample)
            )
        for order in range(4):
            assert_scaled_close(
                "projected sample {} F{}".format(sample, order),
                power[order],
                reference[state][order],
                absolute=3.0e-9,
                relative=3.0e-9,
            )
        if occupation in seen and not np.allclose(
            power, seen[occupation], rtol=3.0e-12, atol=3.0e-12
        ):
            raise AssertionError(
                "projected local power changed for repeated occupation"
            )
        seen[occupation] = power
    if len(seen) < 8:
        raise AssertionError(
            "projected ED fixture sampled only {} occupations".format(
                len(seen)
            )
        )

    hankel_standard_error = check_solver_output(workdir, rows)
    print(
        "Lanczos2 projected ED: PASS "
        "({} samples, {} occupations, NMPTrans={}, k=pi, "
        "identity fidelity={:.3f}, max Hankel SE={:.3f})".format(
            len(rows),
            len(seen),
            len(TRANSLATIONS),
            identity_fidelity,
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
