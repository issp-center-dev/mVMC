from __future__ import print_function

import json
import os
import shutil
import subprocess
import sys

import numpy as np

from runtest_absolute_krylov_classic_ed import (
    BONDS,
    ELECTRONIC_HUND,
    GUTZWILLER_PARAMETER,
    INTER,
    INTRA,
    SITE_COUNT,
    SPIN_EXCHANGE,
    SPIN_HUND,
    assert_finite_power,
    driver_command,
    pair_matrix,
    parse_driver_output,
)
from runtest_bf_fsz import run_vmc, update_modpara
from runtest_lanczos2_real_ed import (
    assert_scaled_close,
    occupation_state,
    read_ls2_oracle,
)


def write_site_definition(workdir, filename, count_name, values):
    with open(os.path.join(workdir, filename), "w") as destination:
        destination.write("=============================================\n")
        destination.write("{}          {}\n".format(count_name, len(values)))
        destination.write("=============================================\n")
        destination.write("============ i coupling =====================\n")
        destination.write("=============================================\n")
        for site, coefficient in enumerate(values):
            destination.write(
                "{:5d} {:25.15f}\n".format(site, coefficient)
            )


def write_bond_definition(workdir, filename, count_name, rows):
    with open(os.path.join(workdir, filename), "w") as destination:
        destination.write("=============================================\n")
        destination.write("{}          {}\n".format(count_name, len(rows)))
        destination.write("=============================================\n")
        destination.write("============ i j coupling ===================\n")
        destination.write("=============================================\n")
        for first, second, coefficient in rows:
            destination.write(
                "{:5d} {:5d} {:25.15f}\n".format(
                    first, second, coefficient
                )
            )


def write_transfer(workdir, complex_case):
    rows = []
    for bond_index, (left, right) in enumerate(BONDS):
        hopping = 0.75
        if complex_case:
            hopping += 1.0j * (0.125 + 0.0625 * bond_index)
        for spin in (0, 1):
            rows.append((left, spin, right, spin, hopping.real, hopping.imag))
            reverse = np.conj(hopping)
            rows.append((right, spin, left, spin,
                         reverse.real, reverse.imag))
    rows.append((1, 0, 1, 0, 0.4375, 0.0))
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


def write_zero_transfer(workdir):
    with open(os.path.join(workdir, "trans.def"), "w") as destination:
        destination.write("========================\n")
        destination.write("NTransfer       0\n")
        destination.write("========================\n")
        destination.write("========i_j_s_tijs======\n")
        destination.write("========================\n")


def write_locspin(workdir, pure_spin):
    with open(os.path.join(workdir, "locspn.def"), "w") as destination:
        destination.write("================================\n")
        destination.write("NlocalSpin     {}\n".format(
            SITE_COUNT if pure_spin else 0
        ))
        destination.write("================================\n")
        destination.write("========i_1LocSpn_0IteElc ======\n")
        destination.write("================================\n")
        for site in range(SITE_COUNT):
            destination.write(
                "{:5d} {:5d}\n".format(site, 1 if pure_spin else 0)
            )


def write_gutzwiller(workdir):
    with open(os.path.join(workdir, "gutzwilleridx.def"), "w") as destination:
        destination.write("=============================================\n")
        destination.write("NGutzwillerIdx          1\n")
        destination.write("ComplexType             0\n")
        destination.write("=============================================\n")
        destination.write("=============================================\n")
        for site in range(SITE_COUNT):
            destination.write("{:5d} {:5d}\n".format(site, 0))
        destination.write("{:5d} {:5d}\n".format(0, 0))


def write_namelist(workdir, pure_spin):
    entries = [
        ("ModPara", "modpara.def"),
        ("LocSpin", "locspn.def"),
        ("Trans", "trans.def"),
    ]
    if pure_spin:
        entries.extend([
            ("Hund", "hund.def"),
            ("Exchange", "exchange.def"),
        ])
    else:
        entries.extend([
            ("CoulombIntra", "coulombintra.def"),
            ("CoulombInter", "coulombinter.def"),
            ("Hund", "hund.def"),
            ("Gutzwiller", "gutzwilleridx.def"),
        ])
    entries.extend([
        ("Orbital", "orbitalidx.def"),
        ("TransSym", "qptransidx.def"),
    ])
    with open(os.path.join(workdir, "namelist.def"), "w") as destination:
        for keyword, filename in entries:
            destination.write("{:>16}  {}\n".format(keyword, filename))


def write_initial(workdir, complex_case, pure_spin):
    values = [0.0] * 6
    if not pure_spin:
        values.extend([GUTZWILLER_PARAMETER, 0.0, 0.0])
    for value in pair_matrix(complex_case).reshape(SITE_COUNT * SITE_COUNT):
        values.extend([value.real, value.imag, 0.0])
    filename = "p3_anchor_{}_{}.dat".format(
        "spin" if pure_spin else "electronic",
        "complex" if complex_case else "real",
    )
    with open(os.path.join(workdir, filename), "w") as destination:
        destination.write(
            " ".join("{:.18e}".format(value) for value in values)
        )
        destination.write("\n")
    return filename


def prepare_case(rootdir, name):
    pure_spin = name.startswith("spin_")
    complex_case = name.endswith("_complex")
    reference_name = (
        "BackFlow_Identity_Complex"
        if complex_case else "BackFlow_Identity_Real"
    )
    reference_dir = os.path.join(rootdir, "data", reference_name)
    suffix = "_mpi" if os.environ.get("MVMC_MPI_PROCS") else ""
    workdir = os.path.join(rootdir, "work", "P3ProductionAnchor_" + name + suffix)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    for filename in ("modpara.def", "orbitalidx.def", "qptransidx.def"):
        shutil.copy2(
            os.path.join(reference_dir, filename),
            os.path.join(workdir, filename),
        )
    write_locspin(workdir, pure_spin)
    if pure_spin:
        write_zero_transfer(workdir)
        write_bond_definition(
            workdir, "hund.def", "NHund",
            tuple((left, right, SPIN_HUND[index])
                  for index, (left, right) in enumerate(BONDS)),
        )
        write_bond_definition(
            workdir, "exchange.def", "NExchange",
            tuple((left, right, SPIN_EXCHANGE[index])
                  for index, (left, right) in enumerate(BONDS)),
        )
    else:
        write_transfer(workdir, complex_case)
        write_site_definition(
            workdir, "coulombintra.def", "NCoulombIntra", INTRA
        )
        write_bond_definition(
            workdir, "coulombinter.def", "NCoulombInter", INTER
        )
        write_bond_definition(
            workdir, "hund.def", "NHund", ELECTRONIC_HUND
        )
        write_gutzwiller(workdir)
    write_namelist(workdir, pure_spin)
    update_modpara(
        workdir,
        {
            "NLanczosMode": "1",
            "NLanczosEstimatorMode": "1",
            "NLanczosStep": "2",
            "NLanczosSupportMode": "1",
            "NVMCCalMode": "1",
            "NVMCSample": "4096",
            "NVMCWarmUp": "64",
            "NVMCInterval": "1",
            "Ncond": "0" if pure_spin else "4",
            "2Sz": "0",
            "NMPTrans": "1",
            "NExUpdatePath": "2" if pure_spin else "0",
            "NSplitSize": "1",
        },
    )
    return workdir, write_initial(workdir, complex_case, pure_spin)


def collect_rows(workdir, dump_name):
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    if not mpi_procs:
        return read_ls2_oracle(os.path.join(workdir, dump_name))
    rows = []
    for rank in range(int(mpi_procs)):
        rows.extend(read_ls2_oracle(os.path.join(
            workdir, dump_name + ".rank{:04d}".format(rank)
        )))
    return rows


def run_reference_driver(driver):
    command, _ = driver_command(driver)
    process = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True, check=False,
    )
    if process.returncode != 0:
        raise AssertionError(
            "P3 reference driver failed:\n{}".format(process.stdout)
        )
    return parse_driver_output(process.stdout)[0]


def validate_manifest(rootdir):
    path = os.path.join(rootdir, "data", "absolute_krylov_p3_manifest.json")
    with open(path) as source:
        manifest = json.load(source)
    if manifest.get("schema_version") != 1:
        raise AssertionError("P3 fixture manifest schema mismatch")
    if manifest.get("fixture_id") != "absolute-krylov-p3-l4-v1" or \
            manifest.get("site_count") != SITE_COUNT or \
            manifest.get("orders") != [0, 1, 2, 3]:
        raise AssertionError("P3 fixture manifest core metadata mismatch")
    expected_amplitude = {
        "qp_total": 1,
        "electronic_projection": "GUTZWILLER_ONLY",
        "spin_projection": "NONE",
        "gutzwiller_parameter": GUTZWILLER_PARAMETER,
        "zero_support_bridge_required": True,
    }
    if manifest.get("amplitude") != expected_amplitude:
        raise AssertionError("P3 fixture amplitude metadata mismatch")
    expected_cases = [
        {
            "name": "electronic_real", "scalar": "real",
            "term_families": [
                "Transfer", "DiagonalTransfer", "CoulombIntra",
                "CoulombInter", "Hund",
            ],
        },
        {
            "name": "electronic_complex", "scalar": "complex",
            "term_families": [
                "Transfer", "DiagonalTransfer", "CoulombIntra",
                "CoulombInter", "Hund",
            ],
        },
        {
            "name": "spin_real", "scalar": "real",
            "term_families": ["Hund", "Exchange"],
        },
        {
            "name": "spin_complex", "scalar": "complex",
            "term_families": ["Hund", "Exchange"],
        },
    ]
    if manifest.get("cases") != expected_cases:
        raise AssertionError("P3 fixture per-case coverage metadata mismatch")
    required = {
        "Transfer", "DiagonalTransfer", "CoulombIntra", "CoulombInter",
        "Hund", "Exchange",
    }
    observed = set()
    for case in manifest.get("cases", []):
        observed.update(case.get("term_families", []))
    if observed != required:
        raise AssertionError(
            "P3 manifest term coverage mismatch: {}".format(sorted(observed))
        )
    if manifest.get("pilot_scan") != {
            "lambda": [1.0, 1.0, 1.0, 1.0],
            "rho": [1e-06, 0.0001, 0.01]}:
        raise AssertionError("P3-O2 pilot scan is not frozen in manifest")
    parameters = manifest.get("model_parameters", {})
    electronic = parameters.get("electronic", {})
    spin = parameters.get("spin", {})
    expected_parameters = {
        "ring_bonds": [list(bond) for bond in BONDS],
        "electronic": {
            "transfer_real": 0.75,
            "transfer_imaginary_by_bond": [
                0.125 + 0.0625 * index for index in range(len(BONDS))
            ],
            "diagonal_transfer": [1, 0, 0.4375],
            "coulomb_intra": list(INTRA),
            "coulomb_inter": [list(row) for row in INTER],
            "hund": [list(row) for row in ELECTRONIC_HUND],
        },
        "spin": {
            "hund_by_bond": list(SPIN_HUND),
            "exchange_by_bond": list(SPIN_EXCHANGE),
        },
    }
    if parameters != expected_parameters or not electronic or not spin:
        raise AssertionError("P3 fixture model parameters drifted from manifest")
    pair_data = manifest.get("pair_matrix", {})
    manifest_real = np.array(pair_data.get("real", []), dtype=complex)
    manifest_complex = np.array(
        pair_data.get("complex_real", []), dtype=float
    ) + 1.0j * np.array(
        pair_data.get("complex_imaginary", []), dtype=float
    )
    if not np.array_equal(manifest_real, pair_matrix(False)) or \
            not np.array_equal(manifest_complex, pair_matrix(True)):
        raise AssertionError("P3 fixture pair matrix drifted from manifest")


def main():
    if len(sys.argv) != 2:
        raise AssertionError(
            "usage: runtest_absolute_krylov_production_anchor.py DRIVER"
        )
    rootdir = os.getcwd()
    validate_manifest(rootdir)
    p3_results = run_reference_driver(sys.argv[1])
    cases = (
        "electronic_real", "electronic_complex", "spin_real", "spin_complex"
    )
    total_rows = 0
    total_states = 0
    for name in cases:
        complex_case = name.endswith("_complex")
        support = {
            state for state, absolute in p3_results[name].items()
            if abs(absolute[0]) > 2.0e-13
        }
        workdir, initial = prepare_case(rootdir, name)
        dump_name = "p3_production_anchor_{}.dat".format(name)
        process = run_vmc(
            rootdir, workdir,
            mpi_procs=os.environ.get("MVMC_MPI_PROCS"),
            init_path=initial,
            extra_env={"MVMC_LANCZOS_ORACLE_DUMP": dump_name},
        )
        if process.returncode != 0:
            raise AssertionError(
                "{} production run failed:\n{}".format(name, process.stdout)
            )
        rows = collect_rows(workdir, dump_name)
        seen = set()
        for sample, occupation, production_power in rows:
            state = occupation_state(occupation)
            if state not in p3_results[name]:
                raise AssertionError(
                    "{} sample {} state {} is outside P3 basis".format(
                        name, sample, state
                    )
                )
            absolute = p3_results[name][state]
            if abs(absolute[0]) <= 2.0e-13:
                raise AssertionError(
                    "{} production sampled zero-support state {}".format(
                        name, state
                    )
                )
            expected = absolute / absolute[0]
            if complex_case:
                expected = np.conj(expected)
            assert_finite_power(
                "{} sample {} production".format(name, sample),
                production_power,
            )
            assert_finite_power(
                "{} sample {} expected".format(name, sample), expected
            )
            for order in range(4):
                assert_scaled_close(
                    "{} sample {} F{}".format(name, sample, order),
                    production_power[order], expected[order],
                    3.0e-10, 3.0e-9,
                )
            seen.add(state)
        if seen != support:
            raise AssertionError(
                "{} production anchor did not enumerate nonzero support: "
                "seen={}/{} missing={} extra={}".format(
                    name, len(seen), len(support),
                    sorted(support - seen), sorted(seen - support)
                )
            )
        total_rows += len(rows)
        total_states += len(seen)
    print(
        "P3 production action anchor passed: 4 cases, {} sampled rows, "
        "{} case-local distinct states".format(total_rows, total_states)
    )


if __name__ == "__main__":
    main()
