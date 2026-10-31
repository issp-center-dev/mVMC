"""Check the BackFlow measurement drift summary and reweighting on a short run.

usage: runtest_bf_drift_log.py <model> --mode <drift|reweight>

The work directory is prepared like ``runtest_bf.py <model> --lanczos-mode 1
--use-nonidentity-init``: non-FSZ BackFlow, first-Lanczos measurement, one QP
projection (the production BackFlow setting). MVMC_MPI_PROCS selects the rank
count as in the other BackFlow tests. MVMC_SAMPLER_* settings inherited from
the environment are discarded so each run sees only its own switches.
"""
import collections
import hashlib
import math
import os
import shutil
import sys

import runtest_bf as bf

PHYSICAL_PREFIXES = ("zvo_out_", "zvo_var_", "zvo_ls_")
LANCZOS_SAMPLES = 128
BINS = 2


def physical(workdir):
    out = os.path.join(workdir, "output")
    result = {}
    for name in sorted(os.listdir(out)):
        if name.startswith(PHYSICAL_PREFIXES) and name.endswith(".dat"):
            with open(os.path.join(out, name), "rb") as fp:
                result[name] = hashlib.sha256(fp.read()).hexdigest()
    return result


def drift_files(workdir):
    return sorted(name for name in os.listdir(workdir)
                  if name.startswith("sampler_drift_r") and name.endswith(".dat"))


def read_rows(path):
    with open(path) as fp:
        return [line.split() for line in fp.read().splitlines()
                if line and line[0] != "#"]


def prepare(rootdir, model, workdir, modpara_updates):
    refdir = os.path.join(rootdir, "data", model)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    bf.copy_def_files(refdir, workdir, include_backflow=True)
    modpara = os.path.join(workdir, "modpara.def")
    bf.update_modpara(modpara, modpara_updates)
    for key, value in modpara_updates.items():
        bf.assert_modpara_value(modpara, key, value)
    nsite = bf.parse_nsite(modpara)
    definition = bf.build_chain_nn_backflow(length=nsite, optimize=False,
                                            antiperiodic=False)
    bf.write_chain_nn_backflow(workdir, length=nsite, optimize=False,
                               compact=False, antiperiodic=False)
    bf.assert_nonidentity_init_layout(workdir)
    orbitalidx = os.path.join(workdir, "orbitalidx.def")
    nslater = bf.parse_norbitalidx(orbitalidx)
    bf.write_orbital_opt_flags(orbitalidx, nsite, nslater, 0)
    init_name = "nonidentity_init.dat"
    projbf = bf.write_nonidentity_init_parameter(
        os.path.join(workdir, init_name), definition.n_proj_bf, nslater)
    if projbf[0] == 1.0 or all(value == 0.0 for value in projbf[1:]):
        raise RuntimeError("non-identity ProjBF initialization is invalid")
    return init_name


class Case(object):
    def __init__(self, rootdir, model, base, mpi_procs):
        self.rootdir = rootdir
        self.model = model
        self.base = base
        self.mpi_procs = mpi_procs
        self.rank_count = int(mpi_procs) if mpi_procs else 1

    def run(self, name, modpara_updates, settings, expect_error=None):
        workdir = self.base + "_" + name
        init_name = prepare(self.rootdir, self.model, workdir, modpara_updates)
        proc = bf.run_vmc(self.rootdir, workdir, self.mpi_procs,
                          init_path=init_name, extra_env=dict(settings),
                          log_name="bf_drift_log_{}.log".format(name))
        if expect_error is None:
            if proc.returncode != 0:
                print(proc.stdout)
                raise RuntimeError("{} failed".format(name))
        else:
            if proc.returncode == 0 or expect_error not in proc.stdout:
                print(proc.stdout)
                raise RuntimeError("{} did not stop with: {}".format(name, expect_error))
            if drift_files(workdir):
                raise RuntimeError("{} left drift files after a rejected start".format(name))
        return workdir, proc.stdout


def check_drift_file(path, rank, rank_count, nsample, bins, max_drift):
    rows = read_rows(path)
    kinds = collections.defaultdict(list)
    for row in rows:
        kinds[row[0]].append(row)
    q, z, r, b, e = (kinds[k] for k in ("Q", "Z", "R", "B", "E"))
    assert len(q) == len(z) == len(r) == len(b) == bins, path
    assert [int(row[1]) for row in q] == list(range(bins)), path
    # With NSplitSize == 1 every rank measures its own NVMCSample configurations.
    for row in q:
        assert len(row) == 11, row
        assert int(row[2]) == nsample, (row, nsample)
        assert int(row[3]) == 0, row
        assert float(row[4]) <= max_drift, (row, max_drift)
        # No sampler hooks exist for BackFlow: accept/refresh/component columns stay 0.
        assert all(int(v) == 0 for v in row[7:11]), row
    for row in z:
        assert len(row) == 16, row
    for row in r:
        assert len(row) == 6 and all(int(v) == 0 for v in row[2:]), row
    for row in b:
        assert len(row) == 4 and all(int(v) == 0 for v in row[2:]), row
    assert len(e) <= 16, path
    assert all(n <= 4 for n in collections.Counter(row[5] for row in e).values()), path
    assert rows[-1][0] == "F", rows[-1]
    assert os.path.getsize(path) < 16384, path
    return q


def drift_mode(case):
    lanczos = {"NVMCCalMode": "1", "NLanczosMode": "1", "NLanczosStep": "1",
               "NVMCSample": str(LANCZOS_SAMPLES), "NVMCWarmUp": "8",
               "NDataQtySmp": str(BINS)}
    work = {}
    for name, settings, flag in [("default", {}, 1),
                                 ("quiet", {"MVMC_SAMPLER_DRIFT_LOG": "0"}, 0),
                                 ("explicit", {"MVMC_SAMPLER_DRIFT_LOG": "1"}, 1)]:
        workdir, out = case.run(name, lanczos, settings)
        expected = "Sampler component repair: off; growth=1000; drift_log={}".format(flag)
        assert expected in out, (name, expected)
        work[name] = workdir
    assert physical(work["default"]) == physical(work["quiet"]) == physical(work["explicit"])
    assert len(physical(work["default"])) >= 2 * BINS
    assert not drift_files(work["quiet"])
    for name in ("default", "explicit"):
        files = drift_files(work[name])
        assert len(files) == case.rank_count, (name, files)
        for rank, fname in enumerate(files):
            check_drift_file(os.path.join(work[name], fname), rank, case.rank_count,
                             LANCZOS_SAMPLES, BINS, 1.0e-9)
    # Repair is not implemented for the BackFlow update path.
    case.run("force-repair", lanczos, {"MVMC_SAMPLER_REPAIR": "1"},
             expect_error="enabled outside")
    case.run("invalid-switch", lanczos, {"MVMC_SAMPLER_DRIFT_LOG": "2"},
             expect_error="accept only 0 or 1")
    # Plain measurement: off unless requested.
    plain = dict(lanczos, NLanczosMode="0")
    workdir, out = case.run("measure-default", plain, {})
    assert "drift_log=0" in out and not drift_files(workdir)
    workdir, out = case.run("measure-explicit", plain, {"MVMC_SAMPLER_DRIFT_LOG": "1"})
    assert "drift_log=1" in out
    files = drift_files(workdir)
    assert len(files) == case.rank_count
    for rank, fname in enumerate(files):
        check_drift_file(os.path.join(workdir, fname), rank, case.rank_count,
                         LANCZOS_SAMPLES, BINS, 1.0e-9)
    # Optimization never initializes the summary.
    optimization = {"NVMCCalMode": "0", "NLanczosMode": "0",
                    "NSROptItrStep": "1", "NSROptItrSmp": "1"}
    workdir, out = case.run("optimization", optimization, {})
    assert "Sampler component repair" not in out and not drift_files(workdir)
    print("BackFlow drift summary PASS ({} rank)".format(case.rank_count))


def first_energy(workdir, bin_index):
    path = os.path.join(workdir, "output", "zvo_out_{:03d}.dat".format(bin_index + 1))
    with open(path) as fp:
        return float(fp.readline().split()[0])


def read_weight_dump(workdir, rank_count):
    rows = {}
    for rank in range(rank_count):
        path = os.path.join(workdir, "weight_dump.dat.rank{:04d}".format(rank))
        with open(path) as fp:
            rows[rank] = [[float(v) for v in line.split()] for line in fp if line.strip()]
    return rows


def reweight_mode(case):
    lanczos = {"NVMCCalMode": "1", "NLanczosMode": "1", "NLanczosStep": "1",
               "NVMCSample": str(LANCZOS_SAMPLES), "NVMCWarmUp": "8",
               "NDataQtySmp": str(BINS)}
    shift = math.log(4.0)
    shift_env = {"MVMC_BF_TEST_STORED_WEIGHT_SHIFT": repr(shift)}
    # reweight=0 is the previous fixed weight: identical to an unspecified key.
    unset, _ = case.run("reweight-unset", lanczos, {})
    off, _ = case.run("reweight-off", dict(lanczos, reweight="0"), {})
    assert physical(unset) == physical(off)
    # The stable sampler stores the amplitude the measurement recomputes:
    # reweighting changes nothing beyond round-off.
    on, _ = case.run("reweight-on", dict(lanczos, reweight="1"), {})
    for rank, fname in enumerate(drift_files(on)):
        check_drift_file(os.path.join(on, fname), rank, case.rank_count,
                         LANCZOS_SAMPLES, BINS, 1.0e-9)
    for b in range(BINS):
        e_off, e_on = first_energy(off, b), first_energy(on, b)
        assert abs(e_off - e_on) <= 1.0e-9 * max(1.0, abs(e_off)), (b, e_off, e_on)
    # A shifted stored weight is ignored with reweight=0 ...
    shifted_off, _ = case.run("shift-off", dict(lanczos, reweight="0"), shift_env)
    assert physical(shifted_off) == physical(off)
    # ... and applied with reweight=1: w = exp(2(direct + projection) - stored).
    shifted_on, _ = case.run("shift-on", dict(lanczos, reweight="1"),
                             dict(shift_env, MVMC_BF_TEST_WEIGHT_DUMP="weight_dump.dat"))
    assert physical(shifted_on) != physical(on)
    dump = read_weight_dump(shifted_on, case.rank_count)
    for rank in range(case.rank_count):
        rows = dump[rank]
        assert len(rows) == BINS * LANCZOS_SAMPLES, (rank, len(rows))
        for sample, stored, direct, projection, w, _energy in rows:
            expected = math.exp(2.0 * (direct + projection) - stored)
            assert abs(w - expected) <= 1.0e-12 * max(1.0, expected), (sample, w, expected)
            target = 0.25 if int(sample) % 2 == 0 else 1.0
            assert abs(w - target) <= 1.0e-9, (sample, w, target)
    for b in range(BINS):
        num = den = 0.0
        for rank in range(case.rank_count):
            for _s, _st, _d, _p, w, energy in dump[rank][b * LANCZOS_SAMPLES:(b + 1) * LANCZOS_SAMPLES]:
                num += w * energy
                den += w
        e_out = first_energy(shifted_on, b)
        assert abs(e_out - num / den) <= 1.0e-10 * max(1.0, abs(e_out)), (b, e_out, num / den)
    for rank, fname in enumerate(drift_files(shifted_on)):
        rows = read_rows(os.path.join(shifted_on, fname))
        q = [r for r in rows if r[0] == "Q"]
        assert len(q) == BINS
        for row in q:
            assert abs(float(row[4]) - shift) <= 1.0e-9, row
            assert int(row[7]) == LANCZOS_SAMPLES // 2, row
    print("BackFlow reweight PASS ({} rank)".format(case.rank_count))


def main():
    if len(sys.argv) < 4 or sys.argv[2] != "--mode" or sys.argv[3] not in ("drift", "reweight"):
        print(__doc__)
        return 2
    model, mode = sys.argv[1], sys.argv[3]
    for key in list(os.environ):
        if key.startswith("MVMC_SAMPLER_"):
            del os.environ[key]
    rootdir = os.getcwd()
    mpi_procs = os.environ.get("MVMC_MPI_PROCS")
    suffix = os.environ.get("MVMC_BF_TEST_WORK_SUFFIX", "")
    suffix = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in suffix)
    base = os.path.join(rootdir, "work", "{}_{}{}".format(model, mode, suffix))
    if mpi_procs:
        base += "_mpi{}".format(mpi_procs)
    case = Case(rootdir, model, base, mpi_procs)
    if mode == "drift":
        drift_mode(case)
    else:
        reweight_mode(case)
    return 0


if __name__ == "__main__":
    sys.exit(main())
