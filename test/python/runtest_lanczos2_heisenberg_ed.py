from __future__ import print_function

import argparse
import math
import os
import shutil
import subprocess
import sys

import numpy as np

from runtest_bf_fsz import run_vmc, update_modpara
from runtest_lanczos2_real_ed import (
    assert_scaled_close,
    check_solver_output,
    occupation_state,
    read_ls2_oracle,
    validate_full_enumeration,
)
from runtest_lanczos_sector_ed import (
    NSITE,
    ap_orbital_amplitude,
    fixed_basis,
    occupied,
    operator_matrix_complex,
    tj_allowed,
)


HUND_BONDS = (
    (0, 1, 0.31),
    (1, 2, -0.27),
    (2, 3, 0.43),
    (3, 0, 0.19),
    (0, 2, -0.11),
)
EXCHANGE_BONDS = (
    (0, 1, -0.47),
    (1, 2, -0.63),
    (2, 3, -0.38),
    (3, 0, -0.56),
    (0, 2, 0.17),
)
HEISENBERG_HANKEL_RESIDUAL_LIMIT = 0.8
HEISENBERG_HANKEL_STANDARD_ERROR_LIMIT = 12.0
ZERO_AMPLITUDE_UPS = (0, 1)
ZERO_AMPLITUDE_DOWNS = (2, 3)
EXACT_SUPPORT_REFERENCE = {
    False: {
        "mu2": 1.621983259441473,
        "restricted_m11": 1.266138735728027,
        "delta2": 0.355844523713446,
        "delta3_abs": 0.583745170414107,
        "p1_full": -1.723240632552358,
        "p1_restricted": -1.703521587320468,
        "p2_full": -1.795545994455297,
        "p2_restricted": -1.777446224885233,
    },
    True: {
        "mu2": 1.642199795817063,
        "restricted_m11": 1.282974521245426,
        "delta2": 0.359225274571637,
        "delta3_abs": 0.589315487817609,
        "p1_full": -1.722520204414543,
        "p1_restricted": -1.702048688627278,
        "p2_full": -1.795455642029881,
        "p2_restricted": -1.777603929611767,
    },
}


def read_support_diagnostic(path):
    with open(path) as source:
        lines = [line.strip() for line in source if line.strip()]
    if len(lines) != 3:
        raise AssertionError(
            "support diagnostic has {} lines".format(len(lines))
        )
    metadata = {}
    for token in lines[0].split()[4:]:
        if "=" in token:
            key, value = token.split("=", 1)
            metadata[key] = value
    columns = lines[1].lstrip("# ").split()
    values = lines[2].split()
    if len(columns) != len(values):
        raise AssertionError(
            "support diagnostic columns={} values={}".format(
                len(columns), len(values)
            )
        )
    return metadata, {
        name: float(value) for name, value in zip(columns, values)
    }


def check_support_diagnostic(
        workdir, rows, expected_mode, expected_result, expected_step=2):
    path = os.path.join(
        workdir, "output", "zvo_ls_support_001.dat"
    )
    metadata, values = read_support_diagnostic(path)
    if expected_result == "pass":
        expected_quality = "support-check-passed-not-proof"
    elif expected_mode == "experimental":
        expected_quality = "biased-diagnostic-only"
    else:
        expected_quality = "invalid-biased-estimator"
    for key, expected in (
        ("step", str(expected_step)),
        ("mode", expected_mode),
        ("result", expected_result),
        ("quality", expected_quality),
        ("scope", "necessary-not-sufficient"),
    ):
        if metadata.get(key) != expected:
            raise AssertionError(
                "support diagnostic {}={} expected {}".format(
                    key, metadata.get(key), expected
                )
            )
    if rows is None:
        return values
    powers = np.array([power for _, _, power in rows], dtype=complex)
    expected = {
        "M02": np.mean(np.conj(powers[:, 0]) * powers[:, 2]),
        "M11": np.mean(np.conj(powers[:, 1]) * powers[:, 1]),
        "M03": np.mean(np.conj(powers[:, 0]) * powers[:, 3]),
        "M12": np.mean(np.conj(powers[:, 1]) * powers[:, 2]),
    }
    for name, target in expected.items():
        observed = complex(values[name + "_re"], values[name + "_im"])
        assert_scaled_close(
            "support diagnostic " + name,
            observed,
            target,
            2.0e-12,
            2.0e-12,
        )
    return values


def generalized_minimum(gram, dimension):
    overlap = gram[:dimension, :dimension]
    hamiltonian = 0.5 * (
        gram[:dimension, 1:dimension + 1]
        + gram[1:dimension + 1, :dimension]
    )
    eigenvalue, eigenvector = np.linalg.eigh(overlap)
    keep = eigenvalue > 1.0e-12 * np.max(eigenvalue)
    transform = eigenvector[:, keep] / np.sqrt(eigenvalue[keep])
    projected = np.dot(
        np.conj(transform.T), np.dot(hamiltonian, transform)
    )
    projected = 0.5 * (projected + np.conj(projected.T))
    return np.linalg.eigvalsh(projected)[0]


def exact_support_negative_control(complex_mode):
    basis = pure_spin_basis()
    orbital = pair_orbitals(complex_mode)
    psi = np.array(
        [pair_amplitude(state, orbital) for state in basis], dtype=complex
    )
    target_index = basis.index(zero_amplitude_state())
    psi[target_index] = 0.0
    hamiltonian = build_hamiltonian(basis, False)
    powers = [psi]
    for _ in range(3):
        powers.append(np.dot(hamiltonian, powers[-1]))
    support = np.abs(psi) > 1.0e-13
    normalization = np.vdot(psi, psi).real
    full = np.array(
        [
            [np.vdot(powers[a], powers[b]) / normalization
             for b in range(4)]
            for a in range(4)
        ]
    )
    restricted = np.array(
        [
            [np.vdot(powers[a][support], powers[b][support])
             / normalization for b in range(4)]
            for a in range(4)
        ]
    )
    observed = {
        "mu2": full[0, 2].real,
        "restricted_m11": restricted[1, 1].real,
        "delta2": (full[0, 2] - restricted[1, 1]).real,
        "delta3_abs": abs(full[0, 3] - restricted[1, 2]),
        "p1_full": generalized_minimum(full, 2),
        "p1_restricted": generalized_minimum(restricted, 2),
        "p2_full": generalized_minimum(full, 3),
        "p2_restricted": generalized_minimum(restricted, 3),
    }
    if len(basis) != 6 or np.count_nonzero(support) != 5:
        raise AssertionError(
            "exact support fixture is {}/{} instead of 5/6".format(
                np.count_nonzero(support), len(basis)
            )
        )
    if abs(full[0, 2] - full[1, 1]) > 5.0e-13 or \
            abs(full[0, 3] - full[1, 2]) > 5.0e-13:
        raise AssertionError("full-basis power identities do not close")
    if abs(restricted[0, 2] - full[0, 2]) > 5.0e-13:
        raise AssertionError("restricted M02 changed despite its v0 factor")
    for name, expected in EXACT_SUPPORT_REFERENCE[complex_mode].items():
        assert_scaled_close(
            "exact support " + name,
            observed[name],
            expected,
            5.0e-13,
            5.0e-13,
        )
    if not observed["p1_restricted"] > observed["p1_full"] or \
            not observed["p2_restricted"] > observed["p2_full"]:
        raise AssertionError("restricted energies did not retain known bias")
    return observed


def write_standard_input(path):
    with open(path, "w") as destination:
        destination.write(
            """L = 4
Lsub = 2
model = "Spin"
lattice = "chain"
J = 1.0
NSROptItrStep = 1
NVMCSample = 128
NVMCWarmUp = 16
2Sz = 0
RndSeed = 2468
"""
        )


def read_definition_count(path, keyword):
    with open(path) as source:
        for line in source:
            columns = line.split()
            if len(columns) >= 2 and columns[0] == keyword:
                return int(columns[1])
    raise AssertionError("{} is missing from {}".format(keyword, path))


def write_standard_initial(path, workdir):
    projector_count = sum(
        read_definition_count(
            os.path.join(workdir, filename), keyword
        )
        for filename, keyword in (
            ("gutzwilleridx.def", "NGutzwillerIdx"),
            ("jastrowidx.def", "NJastrowIdx"),
        )
    )
    orbital_count = read_definition_count(
        os.path.join(workdir, "orbitalidx.def"), "NOrbitalIdx"
    )
    values = [0.0] * 6
    for index in range(projector_count):
        values.extend(
            [0.07 * math.sin(0.31 * (index + 1)), 0.0, 0.0]
        )
    for index in range(orbital_count):
        values.extend(
            [
                math.sin(0.37 * (index + 1))
                + 0.23 * math.cos(0.19 * (index + 2)),
                0.0,
                0.0,
            ]
        )
    with open(path, "w") as destination:
        destination.write(
            " ".join("{:.18e}".format(value) for value in values)
        )
        destination.write("\n")


def run_standard_smoke(rootdir):
    workdir = os.path.join(
        rootdir, "work", "Lanczos2_Heisenberg_Standard_Smoke"
    )
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    write_standard_input(os.path.join(workdir, "StdFace.def"))

    dry_binary = os.path.join(
        rootdir, "..", "..", "src", "mVMC", "vmcdry.out"
    )
    dry_process = subprocess.run(
        [dry_binary, "StdFace.def"],
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if dry_process.returncode != 0:
        raise AssertionError(
            "standard-mode generation failed:\n{}".format(
                dry_process.stdout
            )
        )
    if read_definition_count(
        os.path.join(workdir, "locspn.def"), "NlocalSpin"
    ) != NSITE:
        raise AssertionError(
            "standard-mode input did not generate a pure-spin model"
        )
    if read_definition_count(
        os.path.join(workdir, "trans.def"), "NTransfer"
    ) != 0:
        raise AssertionError(
            "standard-mode pure-spin input generated Transfer terms"
        )
    if read_definition_count(
        os.path.join(workdir, "exchange.def"), "NExchange"
    ) <= 0:
        raise AssertionError(
            "standard-mode pure-spin input generated no Exchange terms"
        )

    update_modpara(
        workdir,
        {
            "NLanczosMode": "1",
            "NLanczosEstimatorMode": "1",
            "NLanczosStep": "2",
            "NVMCCalMode": "1",
            "NVMCSample": "128",
            "NVMCWarmUp": "16",
            "NMPTrans": "1",
            "NSPGaussLeg": "1",
            "NExUpdatePath": "2",
        },
    )
    initial_name = "lanczos2_heisenberg_standard_init.dat"
    write_standard_initial(
        os.path.join(workdir, initial_name), workdir
    )
    process = run_vmc(
        rootdir,
        workdir,
        init_path=initial_name,
        extra_env={
            "MVMC_LANCZOS_ORACLE_DUMP":
                "lanczos2_heisenberg_standard_oracle.dat"
        },
    )
    if process.returncode != 0:
        raise AssertionError(
            "standard-generated Heisenberg step=2 run failed:\n{}".format(
                process.stdout
            )
        )
    required_outputs = (
        "zvo_ls2_out_001.dat",
        "zvo_ls2_moment_001.dat",
    )
    for filename in required_outputs:
        path = os.path.join(workdir, "output", filename)
        if not os.path.exists(path):
            raise AssertionError(
                "standard-generated run did not write {}".format(filename)
            )
        values = np.loadtxt(path)
        if not np.all(np.isfinite(values)):
            raise AssertionError(
                "standard-generated {} contains non-finite data".format(
                    filename
                )
            )
    if os.path.exists(
        os.path.join(workdir, "output", "zvo_ls_out_001.dat")
    ):
        raise AssertionError(
            "standard-generated step=2 run wrote a legacy ls_out file"
        )
    print("Lanczos2 Heisenberg standard-mode smoke: PASS")
    return 0


def pure_spin_basis():
    return [state for state in fixed_basis(2, 2) if tj_allowed(state)]


def pair_orbitals(complex_mode):
    orbital = np.array(
        [
            [
                np.sin(0.37 * (row + 1) * (column + 2))
                + 0.3 * np.cos(0.19 * (row + 2) * (column + 1))
                + (
                    0.2j * np.cos(0.23 * (row + 2) * (column + 1))
                    if complex_mode
                    else 0.0
                )
                for column in range(NSITE)
            ]
            for row in range(NSITE)
        ],
        dtype=complex if complex_mode else float,
    )
    # det(F[ups=(0,1), downs=(2,3)]) = 0.  The corresponding spin
    # configuration is an intermediate state for several Exchange paths.
    orbital[1, 3] = (
        orbital[0, 3] * orbital[1, 2] / orbital[0, 2]
    )
    return orbital


def pair_amplitude(state, orbital):
    if not np.iscomplexobj(orbital):
        return ap_orbital_amplitude(state, orbital)
    ups = [site for site in range(NSITE) if occupied(state, site)]
    downs = [
        site for site in range(NSITE)
        if occupied(state, site + NSITE)
    ]
    matrix = np.array(
        [[orbital[up, down] for down in downs] for up in ups],
        dtype=complex,
    )
    return np.linalg.det(matrix)


def exchange_operators(left, right):
    return (
        [
            ("c", left),
            ("a", right),
            ("c", right + NSITE),
            ("a", left + NSITE),
        ],
        [
            ("c", left + NSITE),
            ("a", right + NSITE),
            ("c", right),
            ("a", left),
        ],
    )


def diagonal_energy(state):
    value = 0.0
    for left, right, coupling in HUND_BONDS:
        parallel = (
            occupied(state, left) * occupied(state, right)
            + occupied(state, left + NSITE)
            * occupied(state, right + NSITE)
        )
        value -= coupling * parallel
    return value


def build_hamiltonian(basis, ising):
    terms = []
    if not ising:
        for left, right, coupling in EXCHANGE_BONDS:
            for operators in exchange_operators(left, right):
                terms.append((coupling, operators))
    hamiltonian = operator_matrix_complex(basis, terms)
    for index, state in enumerate(basis):
        hamiltonian[index, index] += diagonal_energy(state)
    residual = np.max(np.abs(hamiltonian - np.conj(hamiltonian.T)))
    if residual > 2.0e-13:
        raise AssertionError(
            "Heisenberg ED Hamiltonian is not Hermitian: {}".format(
                residual
            )
        )
    return hamiltonian


def zero_amplitude_state():
    state = 0
    for site in ZERO_AMPLITUDE_UPS:
        state |= 1 << site
    for site in ZERO_AMPLITUDE_DOWNS:
        state |= 1 << (site + NSITE)
    return state


def build_reference(complex_mode, ising):
    basis = pure_spin_basis()
    orbital = pair_orbitals(complex_mode)
    psi = np.array(
        [pair_amplitude(state, orbital) for state in basis],
        dtype=complex,
    )
    target = zero_amplitude_state()
    target_index = basis.index(target)
    if abs(psi[target_index]) > 1.0e-13:
        raise AssertionError(
            "zero-overlap fixture amplitude is {}".format(
                psi[target_index]
            )
        )
    psi[target_index] = 0.0

    hamiltonian = build_hamiltonian(basis, ising)
    powers = [psi]
    for _ in range(3):
        powers.append(np.dot(hamiltonian, powers[-1]))
    validate_full_enumeration(powers)

    reference = {}
    unconjugated = {}
    zero_states = set()
    for index, state in enumerate(basis):
        if abs(psi[index]) <= 1.0e-13:
            zero_states.add(state)
            continue
        values = np.array(
            [powers[order][index] / psi[index] for order in range(4)],
            dtype=complex,
        )
        unconjugated[state] = values
        reference[state] = np.conj(values) if complex_mode else values

    bridge_sources = set()
    if not ising:
        first_power = powers[1]
        for source_index, source in enumerate(basis):
            if source not in reference:
                continue
            coupling_to_zero = hamiltonian[target_index, source_index]
            bridge = coupling_to_zero * first_power[target_index]
            if abs(coupling_to_zero) > 1.0e-13 and abs(bridge) > 1.0e-10:
                bridge_sources.add(source)
        if not bridge_sources:
            raise AssertionError(
                "fixture has no nonzero bridge through its zero-overlap state"
            )

    if ising:
        for state, values in reference.items():
            diagonal = diagonal_energy(state)
            expected = np.array(
                [1.0, diagonal, diagonal ** 2, diagonal ** 3],
                dtype=complex,
            )
            if not np.allclose(values, expected, rtol=2.0e-13, atol=2.0e-13):
                raise AssertionError(
                    "Ising local-power reduction failed for state {}".format(
                        state
                    )
                )

    return reference, unconjugated, zero_states, bridge_sources


def write_locspin_rows(workdir, rows, header_count=NSITE):
    with open(os.path.join(workdir, "locspn.def"), "w") as destination:
        destination.write("================================\n")
        destination.write("NlocalSpin     {}\n".format(header_count))
        destination.write("================================\n")
        destination.write("========i_1LocSpn_0IteElc ======\n")
        destination.write("================================\n")
        for site, flag in rows:
            if flag is None:
                destination.write("{}\n".format(site))
            else:
                destination.write("{:5d} {:5d}\n".format(site, flag))


def write_locspin(workdir):
    write_locspin_rows(
        workdir, tuple((site, 1) for site in range(NSITE))
    )


def write_zero_transfer(workdir):
    with open(os.path.join(workdir, "trans.def"), "w") as destination:
        destination.write("========================\n")
        destination.write("NTransfer       0\n")
        destination.write("========================\n")
        destination.write("========i_j_s_tijs======\n")
        destination.write("========================\n")


def write_bond_definition(workdir, filename, count_name, rows):
    with open(os.path.join(workdir, filename), "w") as destination:
        destination.write("=============================================\n")
        destination.write("{}          {}\n".format(count_name, len(rows)))
        destination.write("=============================================\n")
        destination.write("============ i j coupling ===================\n")
        destination.write("=============================================\n")
        for left, right, coupling in rows:
            destination.write(
                "{:5d} {:5d} {:25.15f}\n".format(
                    left, right, coupling
                )
            )


def write_namelist(workdir):
    entries = [
        ("ModPara", "modpara.def"),
        ("LocSpin", "locspn.def"),
        ("Trans", "trans.def"),
        ("Hund", "hund.def"),
        ("Exchange", "exchange.def"),
        ("Orbital", "orbitalidx.def"),
        ("TransSym", "qptransidx.def"),
    ]
    with open(os.path.join(workdir, "namelist.def"), "w") as destination:
        for keyword, filename in entries:
            destination.write(
                "{:>16}  {}\n".format(keyword, filename)
            )


def write_init(workdir, complex_mode):
    values = [0.0] * 6
    for value in pair_orbitals(complex_mode).reshape(NSITE * NSITE):
        values.extend([value.real, value.imag, 0.0])
    filename = (
        "lanczos2_heisenberg_complex_init.dat"
        if complex_mode
        else "lanczos2_heisenberg_real_init.dat"
    )
    with open(os.path.join(workdir, filename), "w") as destination:
        destination.write(
            " ".join("{:.18e}".format(value) for value in values)
        )
        destination.write("\n")
    return filename


def prepare_case(rootdir, complex_mode, ising, workdir_stem=None):
    reference_name = (
        "BackFlow_Identity_Complex"
        if complex_mode
        else "BackFlow_Identity_Real"
    )
    reference_dir = os.path.join(rootdir, "data", reference_name)
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    split_size = os.environ.get("MVMC_NSPLIT_SIZE", "1")
    suffix = "_Complex" if complex_mode else "_Real"
    if ising:
        suffix += "_Ising"
    if split_size != "1":
        suffix += "_Split{}".format(split_size)
    if mpi_procs:
        suffix += "_mpi"
    if workdir_stem is None:
        workdir_name = "Lanczos2_Heisenberg_ED" + suffix
    else:
        workdir_name = workdir_stem
        if mpi_procs:
            workdir_name += "_mpi"
    workdir = os.path.join(rootdir, "work", workdir_name)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    for filename in ("modpara.def", "orbitalidx.def", "qptransidx.def"):
        shutil.copy2(
            os.path.join(reference_dir, filename),
            os.path.join(workdir, filename),
        )

    write_locspin(workdir)
    write_zero_transfer(workdir)
    write_bond_definition(workdir, "hund.def", "NHund", HUND_BONDS)
    write_bond_definition(
        workdir,
        "exchange.def",
        "NExchange",
        () if ising else EXCHANGE_BONDS,
    )
    write_namelist(workdir)
    update_modpara(
        workdir,
        {
            "NLanczosMode": "1",
            "NLanczosEstimatorMode": "1",
            "NLanczosStep": "2",
            "NVMCCalMode": "1",
            "NVMCSample": "512",
            "NVMCWarmUp": "32",
            "NVMCInterval": "1",
            "Ncond": "0",
            "2Sz": "0",
            "NMPTrans": "1",
            "NExUpdatePath": "2",
            "NSplitSize": split_size,
        },
    )
    return workdir, write_init(workdir, complex_mode)


def collect_oracle_rows(workdir, dump_name, mpi_procs):
    if not mpi_procs:
        return read_ls2_oracle(os.path.join(workdir, dump_name))
    rows = []
    for rank in range(int(mpi_procs)):
        rows.extend(
            read_ls2_oracle(
                os.path.join(
                    workdir,
                    dump_name + ".rank{:04d}".format(rank),
                )
            )
        )
    return rows


def run_invalid_locspin(rootdir):
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    cases = (
        (
            "count_mismatch",
            ((0, 1), (1, 1), (2, 0), (3, 0)),
            "NLocalSpin header (4) does not match LocSpin flag count (2)",
        ),
        (
            "invalid_site",
            ((0, 1), (1, 1), (2, 1), (NSITE, 1)),
            "Site index is incorrect",
        ),
        (
            "invalid_flag",
            ((0, 1), (1, 1), (2, 1), (3, 2)),
            "LocSpin flag must be 0 or 1",
        ),
        (
            "duplicate_site",
            ((0, 1), (1, 1), (2, 1), (2, 1)),
            "Duplicate LocSpin site index",
        ),
        (
            "missing_site",
            ((0, 1), (1, 1), (2, 1)),
            "Missing LocSpin site index",
        ),
        (
            "malformed_row",
            ((0, 1), (1, 1), (2, 1), ("not-an-index", None)),
            "Malformed LocSpin definition",
        ),
    )
    for name, rows, expected in cases:
        workdir, init_path = prepare_case(
            rootdir,
            False,
            False,
            workdir_stem="Lanczos2_Heisenberg_InvalidLocSpin",
        )
        write_locspin_rows(workdir, rows)
        process = run_vmc(
            rootdir,
            workdir,
            mpi_procs=mpi_procs,
            init_path=init_path,
        )
        if process.returncode == 0:
            raise AssertionError(
                "invalid LocSpin case {} unexpectedly succeeded".format(name)
            )
        if expected not in process.stdout:
            raise AssertionError(
                "invalid LocSpin case {} missed diagnostic {!r}:\n{}".format(
                    name, expected, process.stdout
                )
            )
    print(
        "Lanczos2 Heisenberg invalid LocSpin{}: PASS ({} cases)".format(
            " MPI" if mpi_procs else "", len(cases)
        )
    )
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--complex", action="store_true")
    parser.add_argument("--ising", action="store_true")
    parser.add_argument("--standard-smoke", action="store_true")
    parser.add_argument("--invalid-locspin", action="store_true")
    args = parser.parse_args()

    rootdir = os.getcwd()
    if args.invalid_locspin:
        if args.complex or args.ising or args.standard_smoke:
            raise AssertionError(
                "--invalid-locspin cannot be combined with other modes"
            )
        return run_invalid_locspin(rootdir)
    if args.standard_smoke:
        if args.complex or args.ising:
            raise AssertionError(
                "--standard-smoke cannot be combined with other modes"
            )
        return run_standard_smoke(rootdir)
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    workdir, init_path = prepare_case(
        rootdir, args.complex, args.ising
    )
    dump_name = "lanczos2_heisenberg_oracle.dat"
    extra_env = {"MVMC_LANCZOS_ORACLE_DUMP": dump_name}
    if not args.ising and not mpi_procs:
        extra_env["MVMC_LANCZOS2_HEISENBERG_AUDIT"] = "1"
    process = run_vmc(
        rootdir,
        workdir,
        mpi_procs=mpi_procs,
        init_path=init_path,
        extra_env=extra_env,
    )
    if not args.ising:
        if process.returncode == 0:
            raise AssertionError(
                "strict zero-support negative control unexpectedly succeeded"
            )
        if "power-Lanczos support mismatch" not in process.stdout:
            raise AssertionError(
                "strict zero-support run missed support mismatch:\n{}".format(
                    process.stdout
                )
            )
        check_support_diagnostic(workdir, None, "strict", "mismatch")
        update_modpara(workdir, {"NLanczosStep": "1"})
        step1_process = run_vmc(
            rootdir,
            workdir,
            mpi_procs=mpi_procs,
            init_path=init_path,
            extra_env=extra_env,
        )
        if step1_process.returncode == 0 or \
                "power-Lanczos support mismatch" not in step1_process.stdout:
            raise AssertionError(
                "strict first-step zero-support negative control failed:\n{}"
                .format(step1_process.stdout)
            )
        check_support_diagnostic(
            workdir, None, "strict", "mismatch", expected_step=1
        )
        for filename in ("zvo_ls_out_001.dat", "zvo_ls_qqqq_001.dat"):
            path = os.path.join(workdir, "output", filename)
            if os.path.exists(path):
                os.remove(path)
        update_modpara(
            workdir,
            {"NLanczosStep": "2", "NLanczosSupportMode": "1"},
        )
        process = run_vmc(
            rootdir,
            workdir,
            mpi_procs=mpi_procs,
            init_path=init_path,
            extra_env=extra_env,
        )
    if process.returncode != 0:
        raise AssertionError(
            "Lanczos2 Heisenberg vmc.out failed:\n{}".format(
                process.stdout
            )
        )
    if not args.ising and not mpi_procs and \
            "Lanczos2 Heisenberg audit:" not in process.stdout:
        raise AssertionError("Heisenberg audit completion marker is missing")

    rows = collect_oracle_rows(workdir, dump_name, mpi_procs)
    reference, unconjugated, zero_states, bridge_sources = build_reference(
        args.complex, args.ising
    )
    if not args.ising:
        exact_support_negative_control(args.complex)
    if len(zero_states) != 1:
        raise AssertionError(
            "expected one zero-overlap state, got {}".format(
                len(zero_states)
            )
        )

    seen = {}
    saw_imaginary = False
    conjugation_sensitive = False
    bridge_exercised = args.ising
    for sample, occupation, power in rows:
        state = occupation_state(occupation)
        if state not in reference:
            raise AssertionError(
                "sample {} has zero or missing ED amplitude".format(sample)
            )
        for order in range(4):
            assert_scaled_close(
                "Heisenberg sample {} F{}".format(sample, order),
                power[order],
                reference[state][order],
            )
            if abs(power[order].imag) > 1.0e-8:
                saw_imaginary = True
            if abs(
                reference[state][order] - unconjugated[state][order]
            ) > 1.0e-6:
                conjugation_sensitive = True
        if occupation in seen and not np.allclose(
            power, seen[occupation], rtol=2.0e-12, atol=2.0e-12
        ):
            raise AssertionError(
                "Heisenberg local power changed for repeated occupation"
            )
        seen[occupation] = power
        if state in bridge_sources:
            bridge_exercised = True

    if len(seen) != len(reference):
        raise AssertionError(
            "Heisenberg fixture sampled {}/{} nonzero occupations".format(
                len(seen), len(reference)
            )
        )
    if not bridge_exercised:
        raise AssertionError(
            "Heisenberg sampling did not exercise a nonzero zero-overlap "
            "bridge source"
        )
    if args.complex and not args.ising:
        if not saw_imaginary:
            raise AssertionError(
                "complex Heisenberg fixture did not exercise Im(F)"
            )
        if not conjugation_sensitive:
            raise AssertionError(
                "complex Heisenberg fixture is insensitive to conjugation"
            )

    check_support_diagnostic(
        workdir,
        rows,
        "strict" if args.ising else "experimental",
        "pass" if args.ising else "mismatch",
    )

    hankel_standard_error = check_solver_output(
        workdir,
        rows,
        HEISENBERG_HANKEL_RESIDUAL_LIMIT,
        require_nontrivial_hankel=not args.ising,
        hankel_standard_error_limit=(
            HEISENBERG_HANKEL_STANDARD_ERROR_LIMIT
        ),
    )
    print(
        "Lanczos2 Heisenberg {}{} ED: PASS "
        "({} samples, {} occupations, max Hankel SE={:.3f})".format(
            "complex" if args.complex else "real",
            " Ising" if args.ising else "",
            len(rows),
            len(seen),
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
