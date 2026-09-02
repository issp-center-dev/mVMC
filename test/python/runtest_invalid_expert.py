from __future__ import print_function

import os
import shutil
import subprocess
import sys


def safe_path_component(value):
    chars = []
    for ch in value:
        if ch.isalnum() or ch in ("-", "_", "."):
            chars.append(ch)
        else:
            chars.append("_")
    return "".join(chars) or "case"


def make_workdir_name(model, updates):
    pieces = [model]
    for key in sorted(updates):
        pieces.append("{}_{}".format(key, updates[key]))
    return safe_path_component("__".join(pieces))


def update_modpara(filename, updates):
    found = set()
    lines = []
    with open(filename) as f:
        for line in f:
            words = line.split()
            if words and words[0] in updates:
                lines.append("{:<15} {}\n".format(words[0], updates[words[0]]))
                found.add(words[0])
            else:
                lines.append(line)
    for key, value in updates.items():
        if key not in found:
            lines.append("{:<15} {}\n".format(key, value))
    with open(filename, "w") as f:
        f.writelines(lines)


def inject_spin_flip_transfer(filename):
    with open(filename) as source:
        lines = source.readlines()
    for index, line in enumerate(lines):
        words = line.split()
        if len(words) != 6:
            continue
        try:
            int(words[0])
            source_spin = int(words[1])
            int(words[2])
            int(words[3])
            float(words[4])
            float(words[5])
        except ValueError:
            continue
        words[3] = str(1 - source_spin)
        lines[index] = " ".join(words) + "\n"
        with open(filename, "w") as destination:
            destination.writelines(lines)
        return
    raise RuntimeError("no Transfer data row was available for mutation")


def inject_two_body_g_ex(namelist):
    filename = "greentwoex.def"
    with open(filename, "w") as destination:
        destination.write("===============================\n")
        destination.write("NCisAjsCktAlt  1\n")
    with open(namelist) as source:
        lines = source.readlines()
    if not any(line.split() and line.split()[0] == "TwoBodyGEx"
               for line in lines):
        lines.append("      TwoBodyGEx  {}\n".format(filename))
    with open(namelist, "w") as destination:
        destination.writelines(lines)


def append_namelist_entry(filename, keyword, definition):
    with open(filename) as source:
        lines = source.readlines()
    lines.append("{:>18}  {}\n".format(keyword, definition))
    with open(filename, "w") as destination:
        destination.writelines(lines)


def write_anomalous_g(filename, count=2, rows=None):
    if rows is None:
        rows = ["1 0 0 1 1\n", "0 1 1 0 0\n"]
    with open(filename, "w") as destination:
        destination.write("=============================================\n")
        destination.write("NAnomalousG {}\n".format(count))
        destination.write("=============================================\n")
        destination.write("============ type s1 spin1 s2 spin2 =========\n")
        destination.write("=============================================\n")
        destination.writelines(rows)


def replace_anomalous_term(old, new):
    with open("anomalousterm.def") as source:
        text = source.read()
    if old not in text:
        raise RuntimeError("anomalous term mutation target not found: {}".format(old))
    with open("anomalousterm.def", "w") as destination:
        destination.write(text.replace(old, new, 1))


def set_anomalous_count(filename, keyword, value):
    with open(filename) as source:
        lines = source.readlines()
    if len(lines) < 2:
        raise RuntimeError("incomplete anomalous definition: {}".format(filename))
    lines[1] = "{} {}\n".format(keyword, value)
    with open(filename, "w") as destination:
        destination.writelines(lines)


def apply_gc_fixture_mutation(action):
    if action == "gc_locspin":
        with open("locspn.def") as source:
            text = source.read()
        text = text.replace("NlocalSpin     0", "NlocalSpin     1")
        text = text.replace("    0      0", "    0      1", 1)
        with open("locspn.def", "w") as destination:
            destination.write(text)
    elif action == "gc_noorbgen":
        with open("namelist.def") as source:
            lines = [line for line in source
                     if not line.split() or line.split()[0] != "OrbitalGeneral"]
        lines.append("           Orbital  orbitalidx.def\n")
        with open("namelist.def", "w") as destination:
            destination.writelines(lines)
    elif action == "gc_realparam":
        with open("orbitalidxgen.def") as source:
            text = source.read()
        text = text.replace("ComplexType          1", "ComplexType          0")
        with open("orbitalidxgen.def", "w") as destination:
            destination.write(text)
    elif action == "gc_updateweight":
        append_namelist_entry("namelist.def", "InUpdateWeight", "updateweight.def")
    elif action == "gc_rbm":
        append_namelist_entry("namelist.def", "GeneralRBM_HiddenLayer", "rbm.def")
    elif action == "gc_backflow":
        append_namelist_entry("namelist.def", "BF", "bf.def")
        append_namelist_entry("namelist.def", "BFRange", "rangebf.def")
    elif action == "gc_opttrans1":
        append_namelist_entry("namelist.def", "OptTrans", "opttrans1.def")
    elif action == "gc_opttrans2":
        append_namelist_entry("namelist.def", "OptTrans", "opttrans2.def")
    elif action == "anomalous_on":
        append_namelist_entry("namelist.def", "AnomalousTerm", "anomalousterm.def")
    elif action == "anomalous_g_on":
        write_anomalous_g("anomalousg.def")
        append_namelist_entry("namelist.def", "AnomalousG", "anomalousg.def")
    elif action == "anomalous_broken_pair":
        replace_anomalous_term("0.30 -0.10", "0.30 0.10")
        append_namelist_entry("namelist.def", "AnomalousTerm", "anomalousterm.def")
    elif action == "anomalous_odd":
        set_anomalous_count("anomalousterm.def", "NAnomalousTerm", 1)
        append_namelist_entry("namelist.def", "AnomalousTerm", "anomalousterm.def")
    elif action == "anomalous_bad_site":
        replace_anomalous_term("1 0 0 1 1", "1 9 0 1 1")
        append_namelist_entry("namelist.def", "AnomalousTerm", "anomalousterm.def")
    elif action == "anomalous_bad_spin":
        replace_anomalous_term("1 0 0 1 1", "1 0 2 1 1")
        append_namelist_entry("namelist.def", "AnomalousTerm", "anomalousterm.def")
    elif action == "anomalous_bad_type":
        replace_anomalous_term("1 0 0 1 1", "2 0 0 1 1")
        append_namelist_entry("namelist.def", "AnomalousTerm", "anomalousterm.def")
    elif action == "anomalous_same_op":
        replace_anomalous_term("1 0 0 1 1", "1 0 0 0 0")
        append_namelist_entry("namelist.def", "AnomalousTerm", "anomalousterm.def")
    elif action == "anomalous_nonfinite":
        replace_anomalous_term("0.30 0.10", "nan 0.10")
        append_namelist_entry("namelist.def", "AnomalousTerm", "anomalousterm.def")
    elif action == "anomalous_extra_token":
        replace_anomalous_term("0.30 0.10\n", "0.30 0.10 9\n")
        append_namelist_entry("namelist.def", "AnomalousTerm", "anomalousterm.def")
    elif action == "anomalous_extra_row":
        with open("anomalousterm.def", "a") as destination:
            destination.write("1 0 1 1 0 0.20 0.00\n")
        append_namelist_entry("namelist.def", "AnomalousTerm", "anomalousterm.def")
    elif action == "anomalous_g_malformed":
        write_anomalous_g("anomalousg.def", count=1, rows=["1 0 0 1\n"])
        append_namelist_entry("namelist.def", "AnomalousG", "anomalousg.def")
    elif action == "anomalous_negative_term_count":
        set_anomalous_count("anomalousterm.def", "NAnomalousTerm", -1)
        append_namelist_entry("namelist.def", "AnomalousTerm", "anomalousterm.def")
    elif action == "anomalous_negative_g_count":
        write_anomalous_g("anomalousg.def", count=-1, rows=[])
        append_namelist_entry("namelist.def", "AnomalousG", "anomalousg.def")
    elif action == "anomalous_oversized_term_count":
        set_anomalous_count("anomalousterm.def", "NAnomalousTerm", 2147483648)
        append_namelist_entry("namelist.def", "AnomalousTerm", "anomalousterm.def")
    elif action == "anomalous_oversized_g_count":
        write_anomalous_g("anomalousg.def", count=2147483648, rows=[])
        append_namelist_entry("namelist.def", "AnomalousG", "anomalousg.def")
    elif action == "anomalous_blank_comment":
        # HPhi-compatible layout: comment lines and blank lines between and
        # after the data rows must be skipped by both parsers.
        with open("anomalousterm.def") as source:
            lines = source.readlines()
        header, rows = lines[:5], lines[5:]
        with open("anomalousterm.def", "w") as destination:
            destination.writelines(header)
            destination.write("# Hermitian pair\n")
            destination.write(rows[0])
            destination.write("\n")
            destination.writelines(rows[1:])
            destination.write("   \n\n")
        write_anomalous_g("anomalousg.def",
                          rows=["# create\n", "1 0 0 1 1\n", "\n",
                                "0 1 1 0 0\n", "\n"])
        append_namelist_entry("namelist.def", "AnomalousTerm", "anomalousterm.def")
        append_namelist_entry("namelist.def", "AnomalousG", "anomalousg.def")
    elif action == "anomalous_counts_zero":
        set_anomalous_count("anomalousterm.def", "NAnomalousTerm", 0)
        with open("anomalousterm.def") as source:
            lines = source.readlines()[:5]
        with open("anomalousterm.def", "w") as destination:
            destination.writelines(lines)
        write_anomalous_g("anomalousg.def", count=0, rows=[])
        append_namelist_entry("namelist.def", "AnomalousTerm", "anomalousterm.def")
        append_namelist_entry("namelist.def", "AnomalousG", "anomalousg.def")
    else:
        raise RuntimeError("unknown GC fixture mutation: {}".format(action))


def main():
    if len(sys.argv) < 3:
        print(
            "usage: {} <model name> <expected_error_substring> "
            "[allow_success] [KEY=VALUE ...]".format(sys.argv[0])
        )
        return -1

    model = sys.argv[1]
    expected = sys.argv[2]
    expected_substrings = [expected]
    allow_success = False
    spin_flip_transfer = False
    two_body_g_ex = False
    opttrans_mode = False
    gc_fixture_mutations = []
    modpara_updates = {}
    for arg in sys.argv[3:]:
        if arg == "allow_success":
            allow_success = True
        elif arg == "spin_flip_transfer":
            spin_flip_transfer = True
        elif arg == "two_body_g_ex":
            two_body_g_ex = True
        elif arg == "opttrans_mode":
            opttrans_mode = True
        elif arg.startswith("gc_") or arg.startswith("anomalous_"):
            gc_fixture_mutations.append(arg)
        elif arg.startswith("expect:"):
            expected_substrings.append(arg[len("expect:"):])
        elif "=" in arg:
            key, value = arg.split("=", 1)
            modpara_updates[key] = value
        else:
            print("ERROR: unknown extra argument: {}".format(arg))
            return -1

    rootdir = os.getcwd()
    refdir = os.path.join(rootdir, "data", model)
    workroot = os.path.join(rootdir, "work")
    if not os.path.exists(workroot):
        os.makedirs(workroot)
    case_name = make_workdir_name(model, modpara_updates)
    workdir = os.path.join(workroot, "{}_{}".format(case_name, os.getpid()))
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    os.chdir(workdir)

    # Copy all definition files required by expert mode.
    for fn in os.listdir(refdir):
        if fn.endswith(".def"):
            shutil.copy(os.path.join(refdir, fn), os.path.join(workdir, fn))

    if modpara_updates:
        update_modpara("modpara.def", modpara_updates)
    if spin_flip_transfer:
        try:
            inject_spin_flip_transfer("trans.def")
        except RuntimeError as error:
            print("ERROR: {}".format(error))
            return -1
    if two_body_g_ex:
        inject_two_body_g_ex("namelist.def")
    for mutation in gc_fixture_mutations:
        try:
            apply_gc_fixture_mutation(mutation)
        except RuntimeError as error:
            print("ERROR: {}".format(error))
            return -1

    bin_to_test = os.path.join(rootdir, "..", "..", "src", "mVMC", "vmc.out")
    command = [bin_to_test]
    if opttrans_mode:
        command.append("-o")
    command.extend(["-e", "namelist.def"])
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    if mpi_procs:
        command = ["mpirun", "-np", mpi_procs] + command
    proc = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    with open("invalid_test.log", "w") as f:
        f.write(proc.stdout)

    for expected_substring in expected_substrings:
        if expected_substring not in proc.stdout:
            print("ERROR: expected error substring not found: {}".format(
                expected_substring))
            print("---- output begin ----")
            print(proc.stdout)
            print("---- output end ----")
            return -1
    # Default: invalid input must terminate with non-zero.
    # Some legacy readers only emit error messages and continue;
    # those cases can opt in to allow_success mode.
    if not allow_success and proc.returncode == 0:
        print("ERROR: vmc.out unexpectedly succeeded for invalid input.")
        return -1
    return 0


if __name__ == "__main__":
    sys.exit(main())
