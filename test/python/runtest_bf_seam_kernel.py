"""Small bounded AP fixtures: coefficients, values, amplitudes and FD h scan."""
import json
import math
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys

from backflow_def_helper import write_chain_nn_backflow
from bf_canonical_model import CanonicalModel, _read_blocks
from bf_seam_oracle import check_coefficients, scaled_error
from runtest_bf import update_modpara, write_nonidentity_init_parameter


def main():
    mode = sys.argv[1]
    real = mode.endswith("_real")
    legacy = mode.startswith("legacy")
    root = Path.cwd()
    binary = (root / "../../src/mVMC/vmc.out").resolve()
    parent = root / "work" / ("bf_seam_kernel_"+mode)
    if parent.exists():
        shutil.rmtree(parent)
    parent.mkdir(parents=True)
    results = []
    for h in (1e-4, 3e-5, 1e-5, 3e-6, 1e-6):
        work = parent / str(h)
        shutil.copytree(root / ("data/BackFlow_Identity_Real" if real else "data/BackFlow_Identity_Complex"), work)
        write_chain_nn_backflow(str(work), antiperiodic=True)
        update_modpara(str(work / "modpara.def"), {
            "NMPTrans": "-2" if mode == "projected" else "-1",
            "NVMCCalMode":"0", "NSROptItrStep":"1", "NSROptItrSmp":"1",
            "NVMCSample":"4", "NVMCWarmUp":"0", "NVMCInterval":"1",
            "RndSeed":"8271", "NExUpdatePath":"1"})
        if mode == "projected":
            path = work / "qptransidx.def"
            header = path.read_text().splitlines()[:5]
            header[1] = "NQPTrans 2"
            rows = ["0 1.0 0.0", "1 0.3 0.2"]
            rows += ["{} {} {} {}".format(q, i, (i+q)%4, -1 if i+q >= 4 else 1)
                     for q in range(2) for i in range(4)]
            path.write_text("\n".join(header+rows)+"\n")
        write_nonidentity_init_parameter(str(work / "initial.def"),10,16,complex_orbitals=not real)
        path = work / "initial.def"
        values = path.read_text().split()
        for i in range(16):
            values[36+3*i] = str(0.4*math.sin((i+1)**2*0.731)+0.3*math.cos((i+2)*0.331))
        path.write_text(" ".join(values)+"\n")
        env = dict(os.environ, MVMC_BF_FD_DUMP="fd.dat",
                   MVMC_BF_TEST_FD_STEP=str(h), MVMC_BF_TEST_COEFFICIENT_DUMP="1",
                   MVMC_BF_PROFILE="1",
                   MVMC_BF_FORCE_CANONICAL_NONFSZ="0" if legacy else "1",
                   MVMC_BF_TEST_FORCE_LEGACY_AP="1" if legacy else "0")
        proc = subprocess.run([str(binary),"-e","namelist.def","initial.def"],cwd=work,
                              env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                              text=True,timeout=60)
        (work / "run.log").write_text(proc.stdout)
        assert proc.returncode == 0, proc.stdout
        expected = "legacy" if legacy else "canonical"
        assert "BackFlow non-FSZ path: "+expected in proc.stdout
        timer = list((work / "output").glob("*CalcTimer.dat"))
        assert len(timer) == 1, timer
        counters = {label.strip(): int(count) for label, count in
                    re.findall(r"(BF[^\n]+?)\s+\[\d+\]\s+(\d+)\s*$", timer[0].read_text(), re.M)}
        assert counters["BF multi-QP legacy incremental"] == 0, counters
        for label in ("BF legacy non-FSZ proposal", "BF legacy non-FSZ accept prep", "BF legacy non-FSZ Green rows"):
            assert (counters[label] > 0) if legacy else (counters[label] == 0), (label, counters)
        if not legacy:
            assert counters["BF canonical full rebuild"] > 0, counters
        dump = work / "fd.dat"
        result = check_coefficients(dump)
        blocks = _read_blocks(dump)
        assert len(blocks) == 4, len(blocks)
        result.update(h=h,blocks=len(blocks),fd_scaled=0.0,amplitude=0.0,slater=0.0,
                      route_counters={k:v for k,v in counters.items() if "legacy" in k or "full rebuild" in k})
        for block in blocks:
            assert int(block["fd_fail_count"][0]) == int(block["nan_count"][0]) == 0
            model = CanonicalModel(block)
            result["slater"] = max(result["slater"],scaled_error(model.build(),model.c_value))
            result["amplitude"] = max(result["amplitude"],scaled_error(model.projected_ip(),model.c_ip))
            for key,value_key in (("max_abs_projbf_fd_diff","max_abs_fd_value"),
                                  ("max_abs_orbital_fd_diff","max_abs_orbital_fd_value")):
                diff,value = float(block[key][0]),float(block[value_key][0])
                result["fd_scaled"] = max(result["fd_scaled"],diff/(2e-7+1e-8*value))
        assert result["slater"] <= 1e-13 and result["amplitude"] <= 1e-13, result
        # Larger h values characterize truncation; the last three are the fixed
        # correctness window. Every block and both parameter families enter max.
        if h <= 1e-5:
            assert result["fd_scaled"] <= 1, result
        results.append(result)
    (parent / "results.json").write_text(json.dumps(results,indent=2)+"\n")
    print(json.dumps(results,indent=2))


if __name__ == "__main__":
    main()
