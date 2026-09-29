"""Compare measurement routes on saved configurations, independent of sampling."""
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

from backflow_def_helper import write_chain_nn_backflow
from runtest_bf import (update_modpara, write_nonidentity_init_parameter,
    write_minimal_twobodyg, write_nbody_failure_def, write_nbodyg_failure_def)


def main():
    mode = sys.argv[1]
    ranks = int(os.environ.get("MVMC_MPI_PROCS", "1"))
    root = Path.cwd()
    parent = root / "work" / ("bf_replay_{}_np{}".format(mode, ranks))
    if parent.exists():
        shutil.rmtree(parent)
    parent.mkdir(parents=True)
    binary = (root / "../../src/mVMC/vmc.out").resolve()
    prefix = [] if ranks == 1 else [os.environ.get("MPIRUN", "mpiexec"), "-np", str(ranks)] + os.environ.get("MVMC_MPI_ARGS", "").split()
    source = parent / "source"
    shutil.copytree(root / ("data/BackFlow_Identity_" + mode.capitalize()), source)
    write_chain_nn_backflow(str(source), antiperiodic=True)
    write_nonidentity_init_parameter(str(source / "initial.def"), 10, 16,
                                    complex_orbitals=mode == "complex")
    path = source / "initial.def"
    values = path.read_text().split()
    for i in range(16):
        values[36+3*i] = str(0.4*math.sin((i+1)**2*0.731)+0.3*math.cos((i+2)*0.331))
    path.write_text(" ".join(values)+"\n")
    update_modpara(str(source / "modpara.def"), {
        "NVMCCalMode":"1", "NLanczosMode":"0", "NMPTrans":"-1",
        "NDataQtySmp":"1", "NVMCSample":"16", "NVMCWarmUp":"0",
        "NVMCInterval":"1", "NExUpdatePath":"1", "RndSeed":"8271",
        "NSplitSize":str(ranks)})
    write_minimal_twobodyg(str(source), 4, all_general=True)
    rows = [(i,s,j,s) for i in range(4) for j in range(4) for s in (0,1)]
    (source / "greenone.def").write_text("====\nNCisAjs 32\n====\n====\n====\n"+
        "\n".join(" ".join(map(str,row)) for row in rows)+"\n")
    with (source / "namelist.def").open("a") as fp:
        fp.write("OneBodyG greenone.def\n")
    if mode == "complex":
        write_nbody_failure_def(str(source))
        write_nbodyg_failure_def(str(source))

    def run(name, record=False, bad=None):
        work = parent / name
        shutil.copytree(source, work)
        legacy = name == "legacy"
        env = dict(os.environ, MVMC_BF_FORCE_CANONICAL_NONFSZ="0" if legacy else "1",
                   MVMC_BF_TEST_FORCE_LEGACY_AP="1" if legacy else "0")
        if record:
            env["MVMC_BF_TEST_RECORD_CONFIG"] = str(parent / "configurations.txt")
        else:
            # Different RNG seed: only the replayed measurement configuration
            # may determine these observables, never the new sampler trajectory.
            update_modpara(str(work / "modpara.def"), {"RndSeed":"98317"})
            env["MVMC_BF_TEST_REPLAY_CONFIG"] = str(bad or parent / "configurations.txt")
        proc = subprocess.run(prefix+[str(binary),"-e","namelist.def","initial.def"],
            cwd=work,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=60)
        (work / "run.log").write_text(proc.stdout)
        if bad is not None:
            assert proc.returncode != 0 and "Error: invalid BackFlow" in proc.stdout, proc.stdout
        else:
            assert proc.returncode == 0, proc.stdout
            if not record:
                assert "measurement replay: samples=16 particles=4" in proc.stdout
        return work

    record = run("record", record=True)
    canonical = run("canonical")
    legacy = run("legacy")
    results = {}
    filenames = ["zvo_out_001.dat", "zvo_cisajs_001.dat", "zvo_cisajscktalt_001.dat"]
    if mode == "complex":
        nbody = list((canonical / "output").glob("*_NBodyG_*.dat"))
        assert len(nbody) == 1, nbody
        assert any(abs(float(x)) > 1e-12 for line in nbody[0].read_text().splitlines()
                   if line.strip() and not line.lstrip().startswith("#") for x in line.split()[-2:])
        filenames.append(nbody[0].name)
    for filename in filenames:
        def numbers(work):
            return [float(x) for line in (work / "output" / filename).read_text().splitlines()
                    if line.strip() and not line.lstrip().startswith("#") for x in line.split()]
        reference = numbers(record)
        assert reference and all(math.isfinite(x) for x in reference)
        for label, work in (("canonical",canonical),("legacy",legacy)):
            actual = numbers(work)
            assert len(actual) == len(reference) and all(math.isfinite(x) for x in actual)
            worst = max(abs(a-b)/(2e-11+1e-10*max(abs(a),abs(b))) for a,b in zip(reference,actual))
            assert worst <= 1, (filename,label,worst)
            results[filename+":"+label] = worst
    original = (parent / "configurations.txt").read_text()
    invalid_site = original.split()
    invalid_site[2] = "999"
    duplicate_site = original.split()
    duplicate_site[3] = duplicate_site[2]
    for name, content in (("truncated",original.rsplit(" ",1)[0]),
                          ("out_of_range"," ".join(invalid_site)),
                          ("duplicate_site"," ".join(duplicate_site)),
                          ("extra",original+" 1\n")):
        path = parent / (name+".txt")
        path.write_text(content)
        run(name,bad=path)
    (parent / "results.json").write_text(json.dumps(results,indent=2)+"\n")
    print(json.dumps(results,indent=2))


if __name__ == "__main__":
    main()
