from __future__ import print_function

import os
import re
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
        dtype=orbital.dtype,
    )
    return np.linalg.det(matrix)


def projected_amplitude(state, orbital):
    return sum(
        weight * translated_ap_amplitude(state, orbital, translation)
        for weight, translation in zip(
            MOMENTUM_PI_WEIGHTS, TRANSLATIONS
        )
    )


def build_projected_reference(complex_mode):
    basis = fixed_basis(2, 2)
    if complex_mode:
        from runtest_lanczos2_complex_ed import (
            complex_ap_orbital_amplitude,
            complex_orbitals,
            hopping_amplitude,
        )
        orbital = complex_orbitals().reshape((NSITE, NSITE))
        amplitude = complex_ap_orbital_amplitude
    else:
        orbital = real_orbitals()
        amplitude = ap_orbital_amplitude
    unprojected = np.array(
        [
            amplitude(state, orbital)
            * np.exp(GUTZWILLER_PARAMETER * doublon_count(state))
            for state in basis
        ],
        dtype=complex if complex_mode else float,
    )
    projected = np.array(
        [
            projected_amplitude(state, orbital)
            * np.exp(GUTZWILLER_PARAMETER * doublon_count(state))
            for state in basis
        ],
        dtype=complex if complex_mode else float,
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
            coefficient = (
                -hopping_amplitude(left, right)
                if complex_mode else -1.0
            )
            terms.append((coefficient, one_body(left, right, spin)))
    site, spin, _, coupling = DIAGONAL_TRANSFER
    terms.append((-coupling, one_body(site, site, spin)))
    terms.extend(
        (coupling, coulomb_intra(site))
        for site, coupling in enumerate(COULOMB_INTRA)
    )
    if complex_mode:
        from runtest_lanczos_sector_ed import operator_matrix_complex
        hamiltonian = operator_matrix_complex(basis, terms)
    else:
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
                (
                    np.conj(powers[order][index] / projected[index])
                    if complex_mode
                    else powers[order][index] / projected[index]
                )
                for order in range(4)
            ],
            dtype=complex,
        )
    return reference, fidelity, orbital


def moved_state(state, left, right, spin):
    source = right + spin * NSITE
    target = left + spin * NSITE
    if not ((state >> source) & 1) or ((state >> target) & 1):
        return None
    return state ^ (1 << source) ^ (1 << target)


def count_zero_component_transitions(sampled_states, orbital):
    count = 0
    for state in sampled_states:
        if abs(projected_amplitude(state, orbital)) <= 1.0e-12:
            continue
        for left, right in RING_HOPS:
            for spin in (0, 1):
                candidate = moved_state(state, left, right, spin)
                if candidate is None:
                    continue
                components = [
                    translated_ap_amplitude(
                        candidate, orbital, translation
                    )
                    for translation in TRANSLATIONS
                ]
                projected = sum(
                    weight * component
                    for weight, component in zip(
                        MOMENTUM_PI_WEIGHTS, components
                    )
                )
                if (
                    abs(projected) > 1.0e-12
                    and any(
                        abs(component) <= 1.0e-12
                        for component in components
                    )
                ):
                    count += 1
    return count


def main():
    arguments = set(sys.argv[1:])
    if arguments - {"--complex"}:
        raise AssertionError(
            "unknown projected ED argument: {}".format(
                sorted(arguments - {"--complex"})
            )
        )
    complex_mode = "--complex" in arguments
    guard_audit = os.environ.get(
        "MVMC_LANCZOS2_TEST_CALHCA_GUARD_AUDIT", "0"
    ) not in ("", "0")
    rootdir = os.getcwd()
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    reference_dir = os.path.join(
        rootdir,
        "data",
        (
            "BackFlow_Identity_Complex"
            if complex_mode else "BackFlow_Identity_Real"
        ),
    )
    suffix = "_mpi" if mpi_procs else ""
    case_name = (
        "Lanczos2_Projected_Complex_ED"
        if complex_mode else "Lanczos2_Projected_ED"
    )
    workdir = os.path.join(rootdir, "work", case_name + suffix)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    for filename in os.listdir(reference_dir):
        source = os.path.join(reference_dir, filename)
        if os.path.isfile(source):
            shutil.copy2(source, os.path.join(workdir, filename))

    write_coulomb_intra(workdir)
    write_lanczos2_namelist(workdir)
    if complex_mode:
        from runtest_lanczos2_complex_ed import (
            write_complex_init,
            write_complex_transfer,
        )
        write_complex_transfer(workdir)
    else:
        write_transfer_with_diagonal(workdir)
    write_gutzwiller(workdir)
    write_translation_projection(workdir)
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
            "NSPGaussLeg": "1",
            "NSPStot": "0",
            "NMPTrans": str(len(TRANSLATIONS)),
        },
    )
    init_path = (
        write_complex_init(workdir)
        if complex_mode else write_real_init(workdir)
    )
    dump_name = (
        "lanczos2_projected_complex_oracle.dat"
        if complex_mode else "lanczos2_projected_oracle.dat"
    )
    process = run_vmc(
        rootdir,
        workdir,
        mpi_procs=mpi_procs,
        init_path=init_path,
        extra_env={
            "MVMC_LANCZOS_ORACLE_DUMP": dump_name,
        },
    )
    if process.returncode != 0:
        raise AssertionError(
            "Lanczos2 projected vmc.out failed:\n{}".format(process.stdout)
        )
    if guard_audit:
        expected_kind = "complex" if complex_mode else "real"
        match = re.search(
            r"Lanczos2 calHCA guard audit: {} direct=(\d+) "
            r"zero_component=(\d+)".format(expected_kind),
            process.stdout,
        )
        if match is None:
            raise AssertionError(
                "Lanczos2 {} calHCA guard audit marker is missing:\n{}".format(
                    expected_kind, process.stdout
                )
            )
        if int(match.group(1)) <= 0 or int(match.group(2)) <= 0:
            raise AssertionError(
                "Lanczos2 {} calHCA guard audit did not exercise the "
                "zero-component branch".format(expected_kind)
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

    reference, identity_fidelity, orbital = build_projected_reference(
        complex_mode
    )
    seen = {}
    sampled_states = set()
    for sample, occupation, power in rows:
        state = occupation_state(occupation)
        sampled_states.add(state)
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
    zero_component_transitions = count_zero_component_transitions(
        sampled_states, orbital
    )
    if zero_component_transitions <= 0:
        raise AssertionError(
            "projected {} ED fixture did not sample a nonzero-overlap "
            "transition with a zero projection component".format(
                "complex" if complex_mode else "real"
            )
        )

    hankel_standard_error = check_solver_output(workdir, rows)
    print(
        "Lanczos2 projected {} ED: PASS "
        "({} samples, {} occupations, NMPTrans={}, k=pi, "
        "identity fidelity={:.3f}, zero-component transitions={}, "
        "max Hankel SE={:.3f})".format(
            "complex" if complex_mode else "real",
            len(rows),
            len(seen),
            len(TRANSLATIONS),
            identity_fidelity,
            zero_component_transitions,
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
