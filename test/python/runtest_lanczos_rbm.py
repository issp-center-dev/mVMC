"""First-step Lanczos output and nonzero RBM weights in the expert interface."""

import math
import os
import sys
from pathlib import Path

import numpy as np

from runtest_real_rbm import (
    configure_rbm_family, configure_short_run, make_case, run_vmc,
    set_all_complex_types, update_modpara,
)
from runtest_varbin import write_def


def make_dimer(root, name, complex_type, mode, lanczos, fixed=False, singular=False):
    # One up and one down electron. The up electron moves within (0,1) or
    # (2,3); the down electron is a spectator. All pair amplitudes are one.
    work = Path(make_case(root, name, base="TJChainL4_Ncond2_J_path4"))
    (work / "namelist.def").write_text(
        "ModPara modpara.def\nLocSpin locspn.def\nTrans trans.def\n"
        "Orbital orbitalidx.def\nInOrbital orbital.def\nTransSym qptransidx.def\n"
        "OneBodyG greenone.def\n")
    write_def(work, "orbitalidx.def", "NOrbitalIdx", 1,
              ["{} {} 0".format(i, j) for i in range(4) for j in range(4)]
              + ["0 1"], complex_type)
    write_def(work, "orbital.def", "NParameter", 1, ["0 1 0"], complex_type)
    write_def(work, "greenone.def", "NCisAjs", 4,
              ["{0} 0 {0} 0".format(i) for i in range(4)])
    trans = []
    potentials = (3.3, 3.3) if singular else (1, 3)
    for i in range(4):
        trans.append("{0} 0 {0} 0 {1} 0".format(i, -potentials[i % 2]))
        if not singular:
            trans.append("{} 0 {} 0 1 0".format(i, i ^ 1))
    write_def(work, "trans.def", "NTransfer", len(trans), trans)
    configure_short_run(work, {"NVMCCalMode": mode, "NLanczosMode": lanczos,
                               "NSPGaussLeg": 1, "NSROptItrStep": 1,
                               "NSROptItrSmp": 1, "NVMCSample": 128,
                               "NVMCInterval": 0 if fixed else 1})
    return work


def check_singular(root):
    ranks = int(os.environ.get("MVMC_MPI_PROCS", "1"))
    for complex_type in (0, 1):
        for lanczos in (1, 2):
            work = make_dimer(root, "singular_{}_{}".format(complex_type, lanczos),
                              complex_type, 1, lanczos, singular=True)
            # H = 3.3 on every sample. Summing this value in binary64
            # gives a positive apparent variance at 10,000 samples.
            split_size = 1 if lanczos == 1 else ranks
            update_modpara(work / "modpara.def", {"NVMCSample": 10000, "NDataQtySmp": 2,
                                                  "NSplitSize": split_size})
            run_vmc(root, work, expected_error="singular or non-finite first-step Lanczos")
            out = work / "output"
            for name in ("out", "var"):
                energy = np.loadtxt(out / ("zvo_{}_001.dat".format(name))).flat[0]
                np.testing.assert_allclose(energy, 3.3, rtol=0, atol=2e-11)
            occupations = np.loadtxt(out / "zvo_cisajs_001.dat")[:, 4]
            np.testing.assert_allclose(sum(occupations), 1, rtol=0, atol=1e-12)
            sample_count = np.loadtxt(out / "zvo_ls_support_001.dat").flat[0]
            if sample_count != 10000 * (ranks // split_size):
                raise AssertionError("unexpected support sample count: {}".format(sample_count))
            if (out / "zvo_out_002.dat").exists():
                raise AssertionError("started a second bin after the singular first bin")
            names = ["out", "qqqq"]
            if lanczos == 2:
                names.extend(["cisajs", "cisajscktalt", "cisajscktaltex"])
            for name in names:
                path = out / ("zvo_ls_{}_001.dat".format(name))
                if path.stat().st_size:
                    raise AssertionError("singular bin wrote a Lanczos result to " + str(path))


def check_rbm(root, complex_type, mode, lanczos):
    # H = [[1,-1],[-1,3]]. A visible bias on even sites changes the
    # even/odd wave-function amplitude ratio from 1 to 2.
    even_occupations = []
    for ratio in (1, 2):
        work = make_dimer(root, "weights_{}_{}_{}_{}".format(complex_type, mode, lanczos, ratio),
                          complex_type, mode, lanczos, fixed=not lanczos)
        configure_rbm_family(work, "General", opt_flag=0)
        write_def(work, "bias.def", "NParameter", 4,
                  ["{} {} 0".format(i, math.log(ratio) / 2 if i % 2 == 0 else 0)
                   for i in range(4)], 0)
        with (work / "namelist.def").open("a") as stream:
            stream.write("InGeneralRBM_PhysLayer bias.def\n")
        if lanczos and not complex_type:
            run_vmc(root, work, expected_error="real-valued RBM does not support Lanczos mode")
            continue
        run_vmc(root, work)
        out = work / "output"
        energy = np.loadtxt(out / "zvo_out_001.dat").flat[0]
        local_e = np.array([1 - 1 / ratio, 3 - ratio])
        if mode == 0:
            # With updates frozen the measured energy is one of the two local
            # values. The sets for zero and nonzero weights are disjoint.
            if min(abs(local_e - energy)) > 1e-11:
                raise AssertionError("unexpected optimization energy {}".format(energy))
            continue
        occupations = np.loadtxt(out / "zvo_cisajs_001.dat")[:, 4]
        weights = [sum(occupations[::2]), sum(occupations[1::2])]
        np.testing.assert_allclose(energy, np.dot(weights, local_e), atol=1e-11, rtol=0)
        if lanczos:
            even_occupations.append(weights[0])
            local_h2 = np.array([2 - 4 / ratio, 10 - 4 * ratio])
            moments = np.loadtxt(out / "zvo_ls_qqqq_001.dat")
            np.testing.assert_allclose(moments[[2, 3, 10, 11, 15]],
                                       np.dot([local_e, local_e**2, local_h2,
                                               local_e*local_h2, local_h2**2], weights),
                                       atol=1e-11, rtol=0)
            if not np.isfinite(np.loadtxt(out / "zvo_ls_out_001.dat")).all():
                raise AssertionError("non-finite Lanczos result")
    # With a fixed seed, the even-site probability rises from 1/2 to 4/5.
    if even_occupations and even_occupations[1] <= even_occupations[0] + 0.1:
        raise AssertionError("RBM sampling did not shift toward even sites: {}".format(even_occupations))


def check_fsz(root):
    for complex_type in (0, 1):
        for mode, lanczos in ((0, 0), (1, 0), (1, 1)):
            work = make_case(root, "fsz_{}_{}_{}".format(complex_type, mode, lanczos),
                             base="HubbardChain_fsz_SpinJastrow")
            configure_rbm_family(work, "General")
            set_all_complex_types(work, complex_type)
            configure_short_run(work, {"NVMCCalMode": mode, "NLanczosMode": lanczos})
            run_vmc(root, work, expected_error="RBM does not support orbital-general (FSZ) mode")


if __name__ == "__main__":
    root = Path.cwd()
    if sys.argv[1] == "singular":
        check_singular(root)
    elif sys.argv[1] == "fsz":
        check_fsz(root)
    elif sys.argv[1] == "weights":
        for complex_type in (0, 1):
            for mode, lanczos in ((0, 0), (1, 0), (1, 1), (1, 2)):
                check_rbm(root, complex_type, mode, lanczos)
    else:
        raise ValueError(sys.argv[1])
