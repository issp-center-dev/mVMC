from __future__ import print_function

import glob
import os
import shutil
import subprocess
import sys

import numpy as np


RBM_NAMELIST_KEYS = (
    "ChargeRBM_HiddenLayer",
    "ChargeRBM_PhysLayer",
    "ChargeRBM_PhysHidden",
    "SpinRBM_HiddenLayer",
    "SpinRBM_PhysLayer",
    "SpinRBM_PhysHidden",
    "GeneralRBM_HiddenLayer",
    "GeneralRBM_PhysLayer",
    "GeneralRBM_PhysHidden",
)


def copy_def_files(source, destination):
    for name in os.listdir(source):
        if name.endswith(".def") and name != "initial.def":
            shutil.copy(os.path.join(source, name), destination)


def update_modpara(filename, updates, remove_prefixes=()):
    found = set()
    result = []
    with open(filename) as source:
        for line in source:
            words = line.split()
            if words and any(words[0].startswith(prefix) for prefix in remove_prefixes):
                continue
            if words and words[0] in updates:
                result.append("{:<18} {}\n".format(words[0], updates[words[0]]))
                found.add(words[0])
            else:
                result.append(line)
    for key, value in updates.items():
        if key not in found:
            result.append("{:<18} {}\n".format(key, value))
    with open(filename, "w") as destination:
        destination.writelines(result)


def read_modpara_int(filename, key):
    with open(filename) as source:
        for line in source:
            words = line.split()
            if len(words) >= 2 and words[0] == key:
                return int(words[1])
    raise RuntimeError("{} not found in {}".format(key, filename))


def set_all_complex_types(workdir, value):
    for filename in glob.glob(os.path.join(workdir, "*.def")):
        result = []
        changed = False
        with open(filename) as source:
            for line in source:
                words = line.split()
                if words and words[0] == "ComplexType":
                    result.append("ComplexType       {}\n".format(value))
                    changed = True
                else:
                    result.append(line)
        if changed:
            with open(filename, "w") as destination:
                destination.writelines(result)


def disable_all_variational_initialization(workdir):
    """Set definition OptFlags to zero so real/complex InitParameter consumes
    the same number of random values before a restart file is loaded."""
    for filename in glob.glob(os.path.join(workdir, "*.def")):
        with open(filename) as source:
            lines = source.readlines()
        if not any(line.split()[:1] == ["ComplexType"] for line in lines):
            continue
        count = None
        for line in lines:
            words = line.split()
            if len(words) >= 2 and words[0].startswith("N"):
                try:
                    count = int(words[1])
                    break
                except ValueError:
                    pass
        if not count:
            continue
        nonempty = [index for index, line in enumerate(lines) if line.split()]
        flag_rows = nonempty[-count:]
        for index in flag_rows:
            words = lines[index].split()
            if len(words) != 2:
                raise AssertionError(
                    "unexpected OptFlag row in {}: {}".format(
                        filename, lines[index].rstrip()
                    )
                )
            lines[index] = "{}\t0\n".format(words[0])
        with open(filename, "w") as destination:
            destination.writelines(lines)


def replace_rbm_namelist(filename, entries):
    result = []
    with open(filename) as source:
        for line in source:
            words = line.split()
            if words and words[0] in RBM_NAMELIST_KEYS:
                continue
            result.append(line)
    for key, value in entries:
        result.append("{:<30} {}\n".format(key, value))
    with open(filename, "w") as destination:
        destination.writelines(result)


def write_layer(filename, count_name, count, complex_type, rows, opt_flags,
                columns):
    with open(filename, "w") as output:
        output.write("--------------------\n")
        output.write("{}\t{}\n".format(count_name, count))
        output.write("ComplexType\t{}\n".format(complex_type))
        output.write("{}\n".format(columns))
        output.write("--------------------\n")
        for row in rows:
            output.write("\t".join(str(value) for value in row) + "\n")
        for index, flag in enumerate(opt_flags):
            output.write("{}\t{}\n".format(index, flag))


def configure_rbm_family(workdir, family, opt_flag=1, mixed_complex=False):
    nsite = read_modpara_int(os.path.join(workdir, "modpara.def"), "Nsite")
    prefix = family.lower() + "_rbm"
    hidden_name = prefix + "_hidden.def"
    phys_name = prefix + "_phys.def"
    coupling_name = prefix + "_coupling.def"
    hidden_path = os.path.join(workdir, hidden_name)
    phys_path = os.path.join(workdir, phys_name)
    coupling_path = os.path.join(workdir, coupling_name)
    family_key = family + "RBM"

    hidden_rows = [(0, 0)]
    if family == "General":
        phys_count = 4
        phys_rows = []
        coupling_rows = []
        for spin in range(2):
            for site in range(nsite):
                class_index = (site + spin * nsite) % phys_count
                phys_rows.append((site, spin, class_index))
                coupling_rows.append((site, spin, 0, class_index))
        phys_columns = "i s RBM_PhysLayer_Idx"
        coupling_columns = "i s k RBM_PhysHidden_Idx"
    else:
        phys_count = 2
        phys_rows = [(site, site % phys_count) for site in range(nsite)]
        coupling_rows = [
            (site, 0, site % phys_count) for site in range(nsite)
        ]
        phys_columns = "i RBM_PhysLayer_Idx"
        coupling_columns = "i k RBM_PhysHidden_Idx"

    write_layer(
        hidden_path,
        "N{}RBM_HiddenLayerIdx".format(family),
        1,
        0,
        hidden_rows,
        [opt_flag],
        "k RBM_HiddenLayer_Idx",
    )
    write_layer(
        phys_path,
        "N{}RBM_PhysLayerIdx".format(family),
        phys_count,
        1 if mixed_complex else 0,
        phys_rows,
        [opt_flag] * phys_count,
        phys_columns,
    )
    write_layer(
        coupling_path,
        "N{}RBM_PhysHiddenIdx".format(family),
        phys_count,
        0,
        coupling_rows,
        [opt_flag] * phys_count,
        coupling_columns,
    )
    replace_rbm_namelist(
        os.path.join(workdir, "namelist.def"),
        [
            (family_key + "_HiddenLayer", hidden_name),
            (family_key + "_PhysLayer", phys_name),
            (family_key + "_PhysHidden", coupling_name),
        ],
    )
    update_modpara(
        os.path.join(workdir, "modpara.def"),
        {"Nneuron" + family: 1},
        remove_prefixes=("Nneuron",),
    )


def remove_rbm(workdir):
    replace_rbm_namelist(os.path.join(workdir, "namelist.def"), [])
    update_modpara(
        os.path.join(workdir, "modpara.def"), {},
        remove_prefixes=("Nneuron",),
    )


def configure_short_run(workdir, extra=None):
    values = {
        "NVMCCalMode": 0,
        "NSROptItrStep": 4,
        "NSROptItrSmp": 2,
        "NVMCWarmUp": 4,
        "NVMCInterval": 1,
        "NVMCSample": 40,
        "NExUpdatePath": 0,
        "NSplitSize": 1,
        "NStore": 1,
        "NSRCG": 0,
        "RndSeed": 39157,
    }
    if extra:
        values.update(extra)
    update_modpara(os.path.join(workdir, "modpara.def"), values)


def make_case(rootdir, case_name, base="GeneralRBM_cmp"):
    suffix = os.environ.get("MVMC_REAL_RBM_TEST_SUFFIX", "")
    workdir = os.path.join(
        rootdir,
        "work",
        "RealRBM_{}_{}{}_{}".format(case_name, os.getpid(), suffix, base),
    )
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    copy_def_files(os.path.join(rootdir, "data", base), workdir)
    set_all_complex_types(workdir, 0)
    return workdir


def vmc_command(rootdir, initial=None):
    binary = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    command = [binary, "-e", "namelist.def"]
    if initial is not None:
        command.append(initial)
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    if mpi_procs:
        command = ["mpirun", "-np", mpi_procs] + command
    return command


def run_vmc(rootdir, workdir, initial=None, expected_error=None):
    proc = subprocess.run(
        vmc_command(rootdir, initial),
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )
    with open(os.path.join(workdir, "real_rbm_test.log"), "w") as output:
        output.write(proc.stdout)
    if expected_error is not None:
        if proc.returncode == 0:
            raise AssertionError("invalid input unexpectedly succeeded")
        if expected_error not in proc.stdout:
            raise AssertionError(
                "expected error {!r} was absent:\n{}".format(
                    expected_error, proc.stdout
                )
            )
        return
    if proc.returncode != 0:
        raise AssertionError(
            "vmc.out failed with code {}:\n{}".format(proc.returncode, proc.stdout)
        )


def load_finite(filename):
    if not os.path.exists(filename):
        raise AssertionError("missing output: {}".format(filename))
    values = np.loadtxt(filename, dtype=float)
    values = np.asarray(values, dtype=float)
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise AssertionError("empty or non-finite output: {}".format(filename))
    return values


def load_numeric_rows(filename):
    if not os.path.exists(filename):
        raise AssertionError("missing output: {}".format(filename))
    rows = []
    with open(filename) as source:
        for line in source:
            words = line.split()
            if not words:
                continue
            try:
                rows.append([float(word) for word in words])
            except ValueError:
                continue
    if not rows:
        raise AssertionError("no numeric rows in output: {}".format(filename))
    values = np.asarray(rows, dtype=float)
    if not np.all(np.isfinite(values)):
        raise AssertionError("non-finite output: {}".format(filename))
    return values


def load_time_counters(filename):
    if not os.path.exists(filename):
        raise AssertionError("missing output: {}".format(filename))
    rows = []
    with open(filename) as source:
        for line in source:
            words = line.split()
            if len(words) < 7:
                continue
            try:
                rows.append([float(word) for word in words[:7]])
            except ValueError:
                continue
    if not rows:
        raise AssertionError("missing acceptance counters: {}".format(filename))
    return np.asarray(rows, dtype=float)


def check_real_parameter_output(workdir, family):
    output = os.path.join(workdir, "output")
    combined = load_finite(os.path.join(output, "zqp_opt.dat")).reshape(-1)
    if combined.size % 3 != 0:
        raise AssertionError("zqp_opt.dat does not contain triplets")
    imaginary = combined.reshape((-1, 3))[:, 1]
    if np.max(np.abs(imaginary)) > 1.0e-13:
        raise AssertionError("all-real zqp_opt.dat has a nonzero imaginary part")
    family_lower = family.lower()
    expected = (
        "zqp_{}RBM_physlayer_opt.dat".format(family_lower),
        "zqp_{}RBM_hiddenlayer_opt.dat".format(family_lower),
        "zqp_{}RBM_physhidden_opt.dat".format(family_lower),
    )
    rbm_values = []
    for name in expected:
        rows = load_numeric_rows(os.path.join(output, name))
        rbm_values.extend(rows[:, 1].tolist())
    if not np.any(np.abs(np.asarray(rbm_values)) > 1.0e-14):
        raise AssertionError("all RBM parameters remained zero")
    load_finite(os.path.join(output, "zvo_out_001.dat"))


def run_family(rootdir, family, extra=None, mixed_complex=False):
    workdir = make_case(rootdir, "family_{}".format(family.lower()))
    configure_rbm_family(workdir, family, mixed_complex=mixed_complex)
    configure_short_run(workdir, extra)
    run_vmc(rootdir, workdir)
    if mixed_complex:
        load_finite(os.path.join(workdir, "output", "zqp_opt.dat"))
        load_finite(os.path.join(workdir, "output", "zvo_out_001.dat"))
    else:
        check_real_parameter_output(workdir, family)
        counters = load_time_counters(
            os.path.join(workdir, "output", "zvo_time_001.dat")
        )
        if counters[-1, 4] <= 0:
            raise AssertionError("the real-RBM fixture made no hopping proposal")
        if extra and int(extra.get("NExUpdatePath", 0)) == 1:
            if counters[-1, 5] <= 0:
                raise AssertionError(
                    "the NExUpdatePath=1 real-RBM fixture made no exchange proposal"
                )


def compare_arrays(label, left, right, atol=2.0e-11, rtol=2.0e-11):
    if left.shape != right.shape or not np.allclose(
        left, right, atol=atol, rtol=rtol
    ):
        max_error = float("inf")
        if left.shape == right.shape:
            max_error = np.max(np.abs(left - right))
        raise AssertionError(
            "{} mismatch: left shape={} right shape={} max_error={}".format(
                label, left.shape, right.shape, max_error
            )
        )


def run_zero_point(rootdir):
    no_rbm = make_case(rootdir, "zero_no_rbm")
    remove_rbm(no_rbm)
    configure_short_run(
        no_rbm,
        {"NSROptItrStep": 1, "NSROptItrSmp": 1, "NVMCSample": 60},
    )

    zero_rbm = make_case(rootdir, "zero_general_rbm")
    configure_rbm_family(zero_rbm, "General", opt_flag=0)
    configure_short_run(
        zero_rbm,
        {"NSROptItrStep": 1, "NSROptItrSmp": 1, "NVMCSample": 60},
    )
    run_vmc(rootdir, no_rbm)
    run_vmc(rootdir, zero_rbm)
    name = "zvo_out_001.dat"
    left = load_finite(os.path.join(no_rbm, "output", name))
    right = load_finite(os.path.join(zero_rbm, "output", name))
    compare_arrays("zero-point " + name, left, right, atol=1.0e-13,
                   rtol=1.0e-13)
    name = "zvo_time_001.dat"
    left = load_time_counters(os.path.join(no_rbm, "output", name))
    right = load_time_counters(os.path.join(zero_rbm, "output", name))
    compare_arrays("zero-point counters", left, right, atol=0.0, rtol=0.0)


def configure_aft(workdir):
    configure_short_run(
        workdir,
        {
            "NVMCCalMode": 1,
            "NVMCWarmUp": 4,
            "NVMCSample": 80,
            "NDataQtySmp": 1,
            "RndSeed": 7193,
        },
    )


def run_equivalence(rootdir):
    optimizer = make_case(rootdir, "equivalence_optimizer")
    configure_rbm_family(optimizer, "General")
    configure_short_run(
        optimizer,
        {"NSROptItrStep": 2, "NSROptItrSmp": 1, "NVMCSample": 60},
    )
    run_vmc(rootdir, optimizer)
    initial_source = os.path.join(optimizer, "output", "zqp_opt.dat")
    load_finite(initial_source)

    real_work = make_case(rootdir, "equivalence_real")
    configure_rbm_family(real_work, "General", opt_flag=0)
    configure_aft(real_work)
    disable_all_variational_initialization(real_work)
    shutil.copy(initial_source, os.path.join(real_work, "restart.def"))
    run_vmc(rootdir, real_work, "restart.def")

    repeat_work = make_case(rootdir, "equivalence_real_repeat")
    configure_rbm_family(repeat_work, "General", opt_flag=0)
    configure_aft(repeat_work)
    disable_all_variational_initialization(repeat_work)
    shutil.copy(initial_source, os.path.join(repeat_work, "restart.def"))
    run_vmc(rootdir, repeat_work, "restart.def")

    complex_work = make_case(rootdir, "equivalence_complex")
    configure_rbm_family(complex_work, "General", opt_flag=0)
    configure_aft(complex_work)
    disable_all_variational_initialization(complex_work)
    set_all_complex_types(complex_work, 1)
    shutil.copy(initial_source, os.path.join(complex_work, "restart.def"))
    run_vmc(rootdir, complex_work, "restart.def")

    outputs = (
        "zvo_out_001.dat",
        "zvo_var_001.dat",
        "zvo_cisajs_001.dat",
        "zvo_cisajscktalt_001.dat",
    )
    for name in outputs:
        real_values = load_finite(os.path.join(real_work, "output", name))
        repeat_values = load_finite(os.path.join(repeat_work, "output", name))
        complex_values = load_finite(os.path.join(complex_work, "output", name))
        compare_arrays("restart " + name, real_values, repeat_values,
                       atol=0.0, rtol=0.0)
        compare_arrays("real/complex " + name, real_values, complex_values)


def run_sr_variant(rootdir, variant):
    settings = {
        "direct": {"NStore": 0, "NSRCG": 0},
        "stored": {"NStore": 1, "NSRCG": 0},
        "cg": {
            "NStore": 1,
            "NSRCG": 1,
            "NSROptCGMaxIter": 1000,
            "DSROptStaDel": 0.1,
            "DSROptCGTol": 1.0e-6,
        },
        "reweight": {"NStore": 1, "NSRCG": 0, "reweight": 1},
        "rescale": {
            "NStore": 1,
            "NSRCG": 1,
            "RescaleSmat": 1,
            "NSROptCGMaxIter": 1000,
            "DSROptStaDel": 0.1,
            "DSROptCGTol": 1.0e-6,
        },
    }[variant]
    run_family(rootdir, "General", settings)


def run_negative(rootdir, case):
    if case == "fsz":
        workdir = make_case(
            rootdir, "negative_fsz", base="HubbardChain_fsz_SpinJastrow"
        )
        expected = "real-valued RBM does not support orbital-general (FSZ) mode"
    else:
        workdir = make_case(rootdir, "negative_" + case)
        expected = {
            "lanczos": "real-valued RBM does not support Lanczos mode",
            "path2": "real-valued RBM does not support NExUpdatePath=2",
        }[case]
    configure_rbm_family(workdir, "General")
    extra = {"NSROptItrStep": 1, "NSROptItrSmp": 1, "NVMCSample": 8}
    if case == "lanczos":
        extra.update({"NVMCCalMode": 1, "NLanczosMode": 1})
    if case == "path2":
        extra["NExUpdatePath"] = 2
    configure_short_run(workdir, extra)
    run_vmc(rootdir, workdir, expected_error=expected)


def main():
    if len(sys.argv) != 2:
        print("usage: {} <case>".format(sys.argv[0]))
        return 2
    case = sys.argv[1]
    rootdir = os.getcwd()
    try:
        if case.startswith("family-"):
            family = case.split("-", 1)[1].capitalize()
            run_family(rootdir, family, {"NExUpdatePath": 1 if family == "General" else 0})
        elif case == "mixed-dispatch":
            run_family(rootdir, "General", mixed_complex=True)
        elif case == "zero-point":
            run_zero_point(rootdir)
        elif case == "equivalence":
            run_equivalence(rootdir)
        elif case.startswith("sr-"):
            run_sr_variant(rootdir, case.split("-", 1)[1])
        elif case.startswith("negative-"):
            run_negative(rootdir, case.split("-", 1)[1])
        else:
            raise AssertionError("unknown real-RBM test case: {}".format(case))
    except (AssertionError, OSError, RuntimeError, ValueError) as error:
        print("ERROR: {}".format(error))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
