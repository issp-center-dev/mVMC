from __future__ import print_function

import hashlib
import json
import math
import os
import shutil
import statistics
import sys

from backflow_def_helper import write_chain_nn_backflow
from runtest_bf import (
    check_multiqp_full_rebuild_profile,
    copy_def_files,
    parse_norbitalidx,
    run_vmc,
    update_modpara,
    write_nonidentity_init_parameter,
    write_orbital_opt_flags,
    write_uniform_gutzwiller,
)


FIXTURE_NAME = "BackFlow_HPhi_Chain_L4"
BASE_MODEL = "BackFlow_Optimization_Complex_MultiQP"


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest(rootdir):
    fixture_root = os.path.join(rootdir, "data", FIXTURE_NAME)
    manifest_path = os.path.join(fixture_root, "manifest.json")
    with open(manifest_path) as source:
        manifest = json.load(source)
    return fixture_root, manifest


def read_hphi_energy(path):
    with open(path) as source:
        for line in source:
            fields = line.split()
            if len(fields) == 2 and fields[0] == "Energy":
                return float(fields[1])
    raise ValueError("missing HPhi Energy in {}".format(path))


def read_sector_fields(path):
    fields = {}
    with open(path) as source:
        for line in source:
            words = line.split()
            if len(words) == 2:
                fields[words[0]] = words[1]
    return fields


def verify_manifest(rootdir, selected_case=None):
    fixture_root, manifest = load_manifest(rootdir)
    if manifest.get("schema_version") != 1:
        raise ValueError("unsupported HPhi fixture manifest schema")
    case_names = [selected_case] if selected_case else sorted(manifest["cases"])
    for case_name in case_names:
        case = manifest["cases"][case_name]
        case_root = os.path.join(fixture_root, case_name)
        for relative, expected in sorted(case["files"].items()):
            path = os.path.join(case_root, relative)
            actual = sha256_file(path)
            if actual != expected:
                raise ValueError(
                    "HPhi fixture checksum mismatch for {}: got {} expected {}".format(
                        path, actual, expected))
        energy = read_hphi_energy(os.path.join(case_root, "zvo_energy.dat"))
        if energy != float(case["exact_ground_energy"]):
            raise ValueError("HPhi fixture energy/manifest mismatch for {}".format(
                case_name))
        sector = read_sector_fields(os.path.join(case_root, "sector_check.txt"))
        residual = float(sector["translation_residual_l2"])
        if residual > float(case["sector_residual_ceiling"]):
            raise ValueError("HPhi sector residual exceeds ceiling for {}".format(
                case_name))
        eigenvalue = complex(
            float(sector["translation_eigenvalue_real"]),
            float(sector["translation_eigenvalue_imag"]))
        expected_eigenvalue = complex(*case["translation_eigenvalue"])
        if abs(eigenvalue - expected_eigenvalue) > float(
                case["sector_residual_ceiling"]):
            raise ValueError("HPhi sector eigenvalue mismatch for {}".format(
                case_name))
        if sector["eigenvector_sha256"] != case["eigenvector_sha256"]:
            raise ValueError("HPhi sector eigenvector checksum mismatch for {}".format(
                case_name))
    return fixture_root, manifest


def write_full_translation(path, case, nsite):
    weights = case["projection_weights"]
    if len(weights) != nsite or int(case["nqptrans"]) != nsite:
        raise ValueError("HPhi fixture must use the complete translation group")
    antiperiodic = int(case["ap_flag"]) != 0
    lines = [
        "=============================================\n",
        "NQPTrans          {}\n".format(nsite),
        "=============================================\n",
        "======== TrIdx_TrWeight_and_TrIdx_i_xi ======\n",
        "=============================================\n",
    ]
    for transform, weight in enumerate(weights):
        lines.append("{:5d} {:.17e} {:.17e}\n".format(
            transform, float(weight[0]), float(weight[1])))
    negative_signs = 0
    for transform in range(nsite):
        for source in range(nsite):
            destination = (source + transform) % nsite
            sign = -1 if antiperiodic and source + transform >= nsite else 1
            negative_signs += sign < 0
            lines.append("{:5d} {:6d} {:6d} {:6d}\n".format(
                transform, source, destination, sign))
    if antiperiodic and negative_signs == 0:
        raise ValueError("AP HPhi fixture has no negative translation signs")
    with open(path, "w") as destination:
        destination.writelines(lines)


def add_namelist_entry(path, keyword, filename):
    with open(path) as source:
        lines = source.readlines()
    if not any(line.split() and line.split()[0] == keyword for line in lines):
        lines.append("{:>16}  {}\n".format(keyword, filename))
    with open(path, "w") as destination:
        destination.writelines(lines)


def prepare_workdir(rootdir, fixture_root, manifest, case_name, seed, stage,
                    optimized_parameter=None):
    case = manifest["cases"][case_name]
    acceptance = manifest["mvmc_acceptance"]
    workdir = os.path.join(
        rootdir, "work", "{}_{}_seed{}_{}".format(
            FIXTURE_NAME, case_name, seed, stage))
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    base = os.path.join(rootdir, "data", BASE_MODEL)
    copy_def_files(base, workdir, include_backflow=True)
    shutil.copy(os.path.join(fixture_root, case_name, "trans.def"),
                os.path.join(workdir, "trans.def"))
    shutil.copy(os.path.join(fixture_root, case_name, "coulombintra.def"),
                os.path.join(workdir, "coulombintra.def"))
    add_namelist_entry(os.path.join(workdir, "namelist.def"),
                       "CoulombIntra", "coulombintra.def")
    write_chain_nn_backflow(workdir, length=4, optimize=True)
    write_uniform_gutzwiller(workdir, 4)
    add_namelist_entry(os.path.join(workdir, "namelist.def"),
                       "Gutzwiller", "gutzwilleridx.def")
    write_full_translation(os.path.join(workdir, "qptransidx.def"), case, 4)
    nslater = parse_norbitalidx(os.path.join(workdir, "orbitalidx.def"))
    write_orbital_opt_flags(
        os.path.join(workdir, "orbitalidx.def"), 4, nslater, 1)

    updates = {
        "NMPTrans": str(int(case["nmptrans"])),
        "RndSeed": str(seed if stage == "opt" else seed + 1000000),
        "NSplitSize": "1",
        "NStore": "1",
    }
    if stage == "opt":
        updates.update(dict((key, str(value))
                            for key, value in acceptance["optimization"].items()))
        updates.update({"NVMCCalMode": "0", "NDataQtySmp": "1"})
        init_name = "nonidentity_init.dat"
        projbf = write_nonidentity_init_parameter(
            os.path.join(workdir, init_name), 10, nslater,
            proj_values=(0.2,))
        if projbf[0] == 1.0 or all(value == 0.0 for value in projbf[1:]):
            raise ValueError("HPhi optimization initial ProjBF is identity")
    elif stage == "measure":
        updates.update(dict((key, str(value))
                            for key, value in acceptance["measurement"].items()))
        updates.update({"NVMCCalMode": "1"})
        init_name = "optimized_init.dat"
        shutil.copy(optimized_parameter, os.path.join(workdir, init_name))
    else:
        raise ValueError("unknown HPhi mVMC stage {}".format(stage))
    update_modpara(os.path.join(workdir, "modpara.def"), updates)
    return workdir, init_name


def output_energies(workdir, expected_count):
    energies = []
    for index in range(1, expected_count + 1):
        path = os.path.join(workdir, "output", "zvo_out_{:03d}.dat".format(index))
        with open(path) as source:
            rows = [line.split() for line in source if line.split()]
        if len(rows) != 1:
            raise ValueError("unexpected mVMC energy row count in {}".format(path))
        value = float(rows[0][0])
        if not math.isfinite(value):
            raise ValueError("non-finite mVMC energy in {}".format(path))
        energies.append(value)
    return energies


def run_case(rootdir, case_name):
    fixture_root, manifest = verify_manifest(rootdir, case_name)
    acceptance = manifest["mvmc_acceptance"]
    measurements = []
    for seed in acceptance["seeds"]:
        opt_workdir, opt_init = prepare_workdir(
            rootdir, fixture_root, manifest, case_name, seed, "opt")
        proc = run_vmc(
            rootdir, opt_workdir, None, init_path=opt_init,
            log_name="hphi_opt.log", extra_env={"MVMC_BF_PROFILE": "1"})
        if proc.returncode != 0:
            print(proc.stdout)
            raise RuntimeError("mVMC HPhi optimization failed for seed {}".format(seed))
        optimized_parameter = os.path.join(opt_workdir, "output", "zqp_opt.dat")
        if not os.path.exists(optimized_parameter):
            raise RuntimeError("mVMC HPhi optimization did not write zqp_opt.dat")
        if check_multiqp_full_rebuild_profile(opt_workdir) != 0:
            raise RuntimeError("mVMC HPhi optimization used an invalid BF route")

        measure_workdir, measure_init = prepare_workdir(
            rootdir, fixture_root, manifest, case_name, seed, "measure",
            optimized_parameter=optimized_parameter)
        proc = run_vmc(
            rootdir, measure_workdir, None, init_path=measure_init,
            log_name="hphi_measure.log", extra_env={"MVMC_BF_PROFILE": "1"})
        if proc.returncode != 0:
            print(proc.stdout)
            raise RuntimeError("mVMC HPhi measurement failed for seed {}".format(seed))
        if check_multiqp_full_rebuild_profile(measure_workdir) != 0:
            raise RuntimeError("mVMC HPhi measurement used an invalid BF route")
        measurements.extend(output_energies(
            measure_workdir, int(acceptance["measurement"]["NDataQtySmp"])))

    if len(measurements) < 2:
        raise RuntimeError("HPhi gate requires at least two independent measurements")
    mean = statistics.mean(measurements)
    stderr = statistics.stdev(measurements) / math.sqrt(len(measurements))
    exact = float(manifest["cases"][case_name]["exact_ground_energy"])
    delta = mean - exact
    lower_sigma = float(acceptance["lower_bound_sigma"])
    lower_slack = float(acceptance["lower_bound_numerical_slack"])
    if mean + lower_sigma * stderr < exact - lower_slack:
        raise RuntimeError(
            "variational lower-bound sanity failed: mean={} stderr={} exact={}".format(
                mean, stderr, exact))
    if delta > float(acceptance["accuracy_tolerance"]):
        raise RuntimeError(
            "mVMC did not reach HPhi accuracy tolerance: delta={}".format(delta))
    if stderr > float(acceptance["uncertainty_ceiling"]):
        raise RuntimeError(
            "mVMC HPhi uncertainty ceiling failed: stderr={}".format(stderr))
    print("case {} exact {:.16e} mean {:.16e} stderr {:.6e} delta {:.6e} samples {}".format(
        case_name, exact, mean, stderr, delta, len(measurements)))
    return 0


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in ("fixture", "pbc", "apbc"):
        print("usage: {} <fixture|pbc|apbc>".format(sys.argv[0]))
        return 1
    rootdir = os.getcwd()
    try:
        if sys.argv[1] == "fixture":
            verify_manifest(rootdir)
            print("HPhi fixture manifest and sector records verified")
            return 0
        return run_case(rootdir, sys.argv[1])
    except (IOError, OSError, RuntimeError, ValueError) as error:
        print("ERROR: {}".format(error))
        return 1


if __name__ == "__main__":
    sys.exit(main())
