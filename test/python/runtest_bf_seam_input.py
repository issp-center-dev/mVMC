"""Exercise the real reader and its MPI failure agreement, including FSZ."""
import os
from pathlib import Path
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

    def run(path, error=None, opt=False):
        cmd = prefix + [str(binary), "-e"] + (["-o"] if opt else []) + ["namelist.def"]
        proc = subprocess.run(cmd, cwd=path, stdout=subprocess.PIPE,
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
    mutations = {
        "ap_missing_phase": (True, 0, "0 0 0\n", "requires a fourth"),
        "mixed_columns": (False, 1, "0 1 1\n", "mixed BFRange"),
        "pbc_negative": (False, 1, "0 3 1 -1\n", "invalid BFRange seam phase"),
        "self_negative": (True, 0, "0 0 0 -1\n", "invalid BFRange seam phase"),
        "zero": (True, 1, "0 3 1 0\n", "invalid BFRange seam phase"),
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
        # Existing reader failures must abort before the seam validator runs.
        qp = path / "qptransidx.def"
        lines = qp.read_text().splitlines(keepends=True)
        lines[-1] = "1 3 999 1\n"
        qp.write_text("".join(lines))
        proc = run(path, "Error", opt=fsz)
        assert "BackFlow seam transform mismatch" not in proc.stdout
    print("BF seam reader passed; ranks={}".format(ranks))
    return 0


if __name__ == "__main__":
    sys.exit(main())
