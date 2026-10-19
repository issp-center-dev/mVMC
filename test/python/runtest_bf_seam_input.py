"""Exercise the real reader and its MPI failure agreement, including FSZ."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from backflow_def_helper import write_chain_nn_backflow
from runtest_bf import update_modpara
from runtest_bf_fsz import prepare_case, make_ap_momentum_projection, make_ap_opt_projection


def main():
    root = Path.cwd()
    ranks = int(os.environ.get("MVMC_MPI_PROCS", "1"))
    workspace = root / "work" / ("bf_seam_input_np" + str(ranks))
    if workspace.exists():
        shutil.rmtree(workspace)
    workspace.mkdir(parents=True)
    binary = (root / "../../src/mVMC/vmc.out").resolve()
    prefix = [] if ranks == 1 else [os.environ.get("MPIRUN", "mpiexec"), "-np", str(ranks)]
    prefix += os.environ.get("MVMC_MPI_ARGS", "").split() if ranks > 1 else []

    def case(name, ap=False, phase=True):
        path = workspace / name
        shutil.copytree(root / "data/BackFlow_Identity_Real", path)
        write_chain_nn_backflow(str(path), antiperiodic=ap, phase_column=phase)
        update_modpara(str(path / "modpara.def"), {
            "NMPTrans": "-1" if ap else "1", "NVMCSample": "4",
            "NVMCWarmUp": "0", "NVMCInterval": "1"})
        for keyword, filename, count_key, row in (
                ("OneBodyG", "greenone.def", "NCisAjs", "0 0 1 0"),
                ("TwoBodyG", "greentwo.def", "NCisAjsCktAltDC", "0 0 1 0 2 1 3 1")):
            (path / filename).write_text("====\n{} 1\n====\n====\n====\n{}\n".format(count_key, row))
            with (path / "namelist.def").open("a") as fp:
                fp.write("{} {}\n".format(keyword, filename))
        return path

    def run(path, error=None, opt=False, force=False):
        cmd = prefix + [str(binary), "-e"] + (["-o"] if opt else []) + ["namelist.def"]
        env = dict(os.environ, MVMC_BF_FORCE_CANONICAL_NONFSZ="1" if force else "0",
                   MVMC_BF_PROFILE="1")
        proc = subprocess.run(cmd, cwd=path, env=env, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True, timeout=30)
        (path / "reader.log").write_text(proc.stdout)
        if error is None:
            if proc.returncode:
                raise AssertionError(str(path) + "\n" + proc.stdout)
        elif proc.returncode == 0 or error not in proc.stdout:
            raise AssertionError("Expected rejection: " + error + "\n" + proc.stdout)
        return proc

    p3, p4 = case("pbc3", phase=False), case("pbc4")
    run(p3); run(p4)
    for filename in ("zvo_out_001.dat", "zvo_cisajs_001.dat", "zvo_cisajscktalt_001.dat"):
        assert (p3 / "output" / filename).read_bytes() == (p4 / "output" / filename).read_bytes(), filename
    run(case("ap4", ap=True))
    # Classify the actual effective transformation, independently of APFlag.
    # A global minus sign is covariant but deliberately ineligible for rows.
    reasons = {"identity": "identity_all_positive", "minus": "nonpositive_transformation_sign",
               "shift": "nonidentity_transformation", "multi": "multiple_transformations",
               "forced": "forced_canonical"}
    for ap in (False, True):
        for complex_mode in (False, True):
            for kind, reason in reasons.items():
                path = case("dispatch_{}_{}_{}".format(ap, complex_mode, kind), ap=ap)
                if complex_mode:
                    orbital = path / "orbitalidx.def"
                    orbital.write_text(re.sub(r"(ComplexType\s+)0", r"\g<1>1", orbital.read_text()))
                count = 2 if kind == "multi" else 1
                update_modpara(str(path / "modpara.def"), {
                    "NMPTrans": str(-count if ap else count), "NSplitSize": str(ranks),
                    "NVMCSample": "16", "RndSeed": "8271"})
                qp = path / "qptransidx.def"
                header = qp.read_text().splitlines()[:5]
                header[1] = "NQPTrans " + str(count)
                rows = ["{} 1.0 0.0".format(q) for q in range(count)]
                for q in range(count):
                    shift = 1 if kind == "shift" else q
                    rows += ["{} {} {} {}".format(q, i, (i+shift)%4,
                        -1 if kind == "minus" or (ap and i+shift >= 4) else 1) for i in range(4)]
                qp.write_text("\n".join(header+rows)+"\n")
                proc = run(path, force=kind == "forced")
                # PBC's projection reader normalizes transformation signs to +1.
                effective_identity = kind == "identity" or (kind == "minus" and not ap)
                expected = "legacy" if effective_identity else "canonical"
                if effective_identity:
                    reason = "identity_all_positive"
                assert "BackFlow non-FSZ path: {} (reason={},".format(expected, reason) in proc.stdout, proc.stdout
                timer, = (path / "output").glob("*CalcTimer.dat")
                counters = {label.strip(): int(value) for label, value in
                            re.findall(r"(BF[^\n]+?)\s+\[\d+\]\s+(\d+)\s*$", timer.read_text(), re.M)}
                assert counters["BF multi-QP legacy incremental"] == 0, counters
                for label in ("BF legacy non-FSZ proposal", "BF legacy non-FSZ accept prep", "BF legacy non-FSZ Green rows"):
                    assert (counters[label] > 0) if expected == "legacy" else (counters[label] == 0), counters
                if expected == "canonical":
                    assert counters["BF canonical full rebuild"] > 0, counters
    mutations = {
        "ap_missing_phase": (True, 0, "0 0 0\n", "requires a fourth"),
        "mixed_columns": (False, 1, "0 1 1\n", "mixed BFRange"),
        "pbc_negative": (False, 1, "0 3 1 -1\n", "invalid BFRange seam phase"),
        "self_negative": (True, 0, "0 0 0 -1\n", "invalid BFRange seam phase"),
        "zero": (True, 1, "0 3 1 0\n", "invalid BFRange seam phase"),
        "two": (True, 1, "0 3 1 2\n", "invalid BFRange seam phase"),
        "inline_comment": (True, 0, "0 0 0 1 #comment\n", "failed to read BFRange"),
        "asymmetric": (True, 1, "0 3 1 1\n", "must be symmetric"),
        "decimal": (True, 0, "0 0 0 1.0\n", "failed to read BFRange"),
        "overflow": (True, 0, "0 0 0 2147483648\n", "failed to read BFRange"),
        "suffix": (True, 0, "0 0 0 1junk\n", "failed to read BFRange"),
        "extra_column": (True, 0, "0 0 0 1 1\n", "failed to read BFRange"),
        "empty": (True, 0, "\n", "failed to read BFRange"),
        "comment": (True, 0, "# comment\n", "failed to read BFRange"),
        "embedded_nul": (True, 0, "0 0 0 1\x00junk\n", "failed to read BFRange"),
        "long_line": (True, 0, "0 0 0 1" + " " * 300 + "\n", "failed to read BFRange"),
        "negative_site": (True, 0, "-1 0 0 1\n", "site index out of range"),
        "trailing": (True, 12, "\n", "extra data"),
    }
    for name, (ap, row, content, error) in mutations.items():
        path = case(name, ap=ap)
        lines = (path / "rangebf.def").read_text().splitlines(keepends=True)
        if row == 12:
            lines.append(content)
        else:
            lines[10 + row] = content
        (path / "rangebf.def").write_text("".join(lines))
        run(path, error)

    for fsz in (False, True):
        path = Path(prepare_case(str(root), "bf_seam_transform_" + str(fsz) + "_np" + str(ranks), True)) if fsz else case("qp_transform", ap=True)
        update_modpara(str(path / "modpara.def"), {"NMPTrans": "-2", "NVMCSample": "1", "NVMCWarmUp": "0"})
        make_ap_momentum_projection(str(path))
        if fsz:
            make_ap_opt_projection(str(path))
        run(path, opt=fsz)
        if fsz:
            make_ap_opt_projection(str(path), all_positive=True)
        else:
            make_ap_momentum_projection(str(path), all_positive=True)
        run(path, "BackFlow seam transform mismatch", opt=fsz)
        if not fsz:
            # An unused, structurally valid row is outside V-7's contract.
            update_modpara(str(path / "modpara.def"), {"NMPTrans": "-1"})
            run(path)
            # A used row is checked even when its projection weight is zero.
            update_modpara(str(path / "modpara.def"), {"NMPTrans": "-2"})
            qp = path / "qptransidx.def"
            lines = qp.read_text().splitlines(keepends=True)
            lines[6] = "1 0.0 0.0\n"
            qp.write_text("".join(lines))
            run(path, "BackFlow seam transform mismatch")
        # Existing structural failures must abort before V-7 sees the table.
        make_ap_momentum_projection(str(path))
        if fsz:
            make_ap_opt_projection(str(path))
        for filename in (("qptransidx.def", "qpopttrans.def") if fsz else ("qptransidx.def",)):
            table = path / filename
            original = table.read_text().splitlines(keepends=True)
            for mutation, error in (("range", "index out of range"),
                                    ("duplicate", "duplicate projection source"),
                                    ("missing", "mapping row count mismatch")):
                lines = list(original)
                if mutation == "range":
                    fields = lines[-1].split()
                    fields[2] = "999"
                    lines[-1] = " ".join(fields)+"\n"
                elif mutation == "duplicate":
                    lines[-1] = lines[-2]
                else:
                    lines.pop()
                table.write_text("".join(lines))
                proc = run(path, error, opt=fsz)
                assert "BackFlow seam transform mismatch" not in proc.stdout
                (path / (filename+"_"+mutation+".log")).write_text(proc.stdout)
            table.write_text("".join(original))
        run(path, opt=fsz)
    for ending in ("no_final_newline", "premature_eof"):
        path = case(ending, ap=True)
        table = path / "rangebf.def"
        if ending == "no_final_newline":
            table.write_text(table.read_text().rstrip("\n"))
            run(path)
        else:
            table.write_text("\n".join(table.read_text().splitlines()[:-1])+"\n")
            run(path,"failed to read BFRange")
    if ranks > 1:
        # Only world rank 0 owns dispatch flag parsing. Deliberately invalid
        # non-root values must neither change the route nor cause a deadlock.
        wrapper = workspace / "rank_flags.py"
        wrapper.write_text('''import os, sys
rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", os.environ.get("PMI_RANK", "-1")))
assert rank >= 0
mode = sys.argv[1]
os.environ["MVMC_BF_FORCE_CANONICAL_NONFSZ"] = ("1" if mode == "canonical" else "0") if rank == 0 else "invalid_nonroot"
if mode == "invalid" and rank == 0:
    os.environ["MVMC_BF_FORCE_CANONICAL_NONFSZ"] = "invalid_root"
os.execv(sys.argv[2], sys.argv[2:])
''')
        for mode in ("canonical", "legacy", "invalid"):
            path = case("rank_flags_"+mode, ap=True)
            proc = subprocess.run(prefix+[sys.executable,str(wrapper),mode,str(binary),"-e","namelist.def"],
                cwd=path,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=30)
            (path / "reader.log").write_text(proc.stdout)
            if mode == "invalid":
                assert proc.returncode != 0 and "path flags require 0 or 1" in proc.stdout, proc.stdout
            else:
                assert proc.returncode == 0 and "BackFlow non-FSZ path: "+mode in proc.stdout, proc.stdout
    print("BF seam reader passed; ranks={}".format(ranks))
    return 0


if __name__ == "__main__":
    sys.exit(main())
