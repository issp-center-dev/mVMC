#!/usr/bin/env python3
"""Finite-S first Lanczos: full-run snapshot replay and independent statistics.

The Fock-space oracle never uses the QP Slater matrices from production.
All sampled snapshots are checked before any statistical acceptance test.
"""
import argparse
import gzip
import json
import math
import pathlib
import re
import shutil
import tempfile

import numpy as np

import grandcanonical_lanczos_oracle as gl
import grandcanonical_spin_oracle as so
import runtest_gc_lanczos as base
from runtest_gc_finite_spin import student_t_quantile

PILOT_SEEDS = (11, 23)
CONFIRM_SEEDS = (101, 211, 307, 401, 503, 601, 701, 809)
STAT_CASES = ("general", "triplet", "maximal", "multi", "pbc_k1", "apbc_k1")


def read_trace(path):
    if path.exists():
        return path.read_text()
    with gzip.open(str(path) + ".gz", "rt") as stream:
        return stream.read()


def compress_passed_traces(workdir, previous_traces=()):
    """Retain exact evidence compactly, only after the entire mode passes."""
    for path in sorted(pathlib.Path(workdir).glob("*/gc-b0-lanczos.chain*")):
        if path.suffix == ".gz" or path in previous_traces:
            continue
        with path.open("rb") as source, gzip.open(str(path) + ".gz", "xb") as target:
            shutil.copyfileobj(source, target)
        path.unlink()


def fixture(name):
    source = {"minimal": "L2_S1", "general": "L4_S1_general",
              "triplet": "L4_S1_triplet_stat", "maximal": "L4_S2",
              "multi": "L5_S2", "pbc_k1": "L4_S1_general",
              "apbc_k1": "L4_S1_general", "non_su2": "L4_S2",
              "upper": "L6_S3"}[name]
    F, spin = so.finite_spin_fixture(source)
    L = len(F)
    boundary = "apbc" if name == "apbc_k1" else "pbc"
    translations = base.momentum(L, 1) if name.endswith("k1") else None
    rows = base.model_rows(L, base.ring_bonds(L), boundary)
    if name == "non_su2":
        rows["trans"] = [(i, s, j, t, v * (0.5 if s else 1.0) if i != j else v)
                         for i, s, j, t, v in rows["trans"]]
    case = base.fixture(name, L, "antiparallel", F, boundary, translations,
                        gauss=4, rows=rows)
    case["S"] = spin
    return case


def close(label, value, reference, tol=2e-8):
    value, reference = np.asarray(value), np.asarray(reference)
    base.require(np.all(np.isfinite(value)) and np.all(np.isfinite(reference)),
                 label + ": nonfinite")
    delta = np.abs(value - reference) / (1 + np.abs(reference))
    base.require(np.all(delta <= tol), "%s: error %.4g" % (label, np.max(delta)))


def reference(case):
    basis, psi = gl.reference_wave(case)
    h = gl.hamiltonian(case, basis)
    q = np.array([psi, h @ psi, h @ (h @ psi)])
    base.require(np.vdot(psi, psi).real > 0, "zero reference norm")
    return basis, h, q


def decode_scaled(fields):
    state = int(fields[0])
    phase = complex(float(fields[1]), float(fields[2]))
    logabs = float(fields[3])
    base.require(state in (0, 1, 2), "invalid scaled state")
    if state == 0:
        base.require(math.isfinite(logabs) and abs(abs(phase) - 1) < 1e-10,
                     "invalid finite phase/log")
        return state, phase * math.exp(logabs)
    return state, 0j


def require_run_metadata(result, summary):
    base.require(np.array_equal(result["number"]["c"], summary["c"]),
                 "public coefficients differ from sampled state")
    for entry in summary["bootstrap"]:
        base.require(len(entry) == 5 and entry[0] >= 1 and entry[4] == 0,
                     "invalid successful bootstrap diagnostics")


def replay(case, work, result, summary, samples, chains, ref):
    require_run_metadata(result, summary)
    basis, h, q = ref
    index = {m: k for k, m in enumerate(basis)}
    c = summary["c"]
    phi = c[0] * q[0] + c[1] * q[1]
    hphi = c[0] * q[1] + c[1] * q[2]
    records, energies, numbers = [], [], []
    node_samples = 0
    for chain in range(chains):
        lines = read_trace(work / ("gc-b0-lanczos.chain%d" % chain)).splitlines()
        base.require(lines[-1] == "END 0", "incomplete trace")
        final = [line.split() for line in lines if line.startswith("F1 ")]
        base.require(len(final) == samples, "wrong final snapshot inventory")
        number_trace = [int(row.split()[3]) for row in lines
                        if row.startswith("S final-sample ")]
        previous = None
        for sample, row in enumerate(final):
            word_count = int(row[4])
            base.require(word_count == (2 * case["L"] + 63) // 64 and
                         len(row) == 27 + word_count and int(row[1]) == sample,
                         "invalid trace shape/index")
            mask = sum(int(row[5 + w], 16) << (64 * w) for w in range(word_count))
            base.require(mask in index, "sample outside domain")
            generation, integrity = int(row[2]), int(row[3])
            start = 5 + word_count
            energy = complex(float(row[start]), float(row[start + 1]))
            scaled = [decode_scaled(row[start + 2 + 4 * k:start + 6 + 4 * k])
                      for k in range(5)]
            k = index[mask]
            expected = list(q[:, k]) + [phi[k], hphi[k]]
            close("snapshot q/Phi/HPhi", [v for _, v in scaled], expected)
            base.require(scaled[3][0] == 0 and abs(phi[k]) > 0,
                         "sampled zero final amplitude")
            close("sample local energy", energy, hphi[k] / phi[k])
            if gl.support_zero(case, mask):
                base.require(scaled[0][0] == 1, "base structural zero not exact")
                node_samples += 1
            n = gl.popcount(mask)
            base.require(n == number_trace[sample], "N trace mismatch")
            if previous is not None:
                pg, pi, pm = previous
                base.require(pg <= generation <= pg + case["L"], "invalid generation step")
                if pg == generation:
                    base.require((pi, pm) == (integrity, mask), "rejected state changed")
            previous = generation, integrity, mask
            records.append((chain, sample, mask, generation))
            energies.append(energy)
            numbers.append(n)
    energies = np.asarray(energies)
    numbers = np.asarray(numbers)
    close("final energy replay", np.mean(energies),
          complex(result["energy"], result["energy_imag"]))
    counts = np.bincount(numbers, minlength=2 * case["L"] + 1).tolist()
    base.require(counts == summary["sectors"] == result["number"]["sectors"],
                 "sector inventory mismatch")
    blocks = base.block_count(samples * chains)
    block_n = numbers.reshape(blocks, -1)
    sums = np.array([block_n.sum(axis=1), (block_n ** 2).sum(axis=1)]).T
    base.require(np.array_equal(sums, result["number"]["block_sums"]),
                 "all-sample N/N2 block mismatch")
    weight = abs(phi) ** 2
    weight /= weight.sum()
    ns = np.asarray([gl.popcount(m) for m in basis])
    exact_e = np.vdot(phi, hphi) / np.vdot(phi, phi)
    exact_e2 = np.vdot(hphi, hphi).real / np.vdot(phi, phi).real
    observables = {"E": (result["energy"], exact_e.real),
                   "ImE": (result["energy_imag"], exact_e.imag),
                   # Physical H^2 comes from coefficient matrices, including
                   # H Phi on nodes; mean |E_local|^2 is not a substitute.
                   "H2": (result["variance"] +
                          (complex(result["energy"], result["energy_imag"]) ** 2).real,
                          exact_e2),
                   "N": (float(numbers.mean()), float(weight @ ns)),
                   "N2": (float((numbers ** 2).mean()), float(weight @ ns ** 2))}
    for n in range(0, 2 * case["L"] + 1, 2):
        observables["P%d" % n] = (float(np.mean(numbers == n)),
                                  float(weight[ns == n].sum()))
    curves = {}
    for block in (64, 128, 256, 512):
        curves[str(block)] = {"E_se": float(np.std(energies.real.reshape(-1, block).mean(axis=1), ddof=1) /
                                             math.sqrt(len(energies) // block)),
                              "N_se": float(np.std(numbers.reshape(-1, block).mean(axis=1), ddof=1) /
                                             math.sqrt(len(numbers) // block))}
    return dict(records=records, node_samples=node_samples, observables=observables,
                block_curves=curves, ess=result["effective_samples"], tau=result["tau_int"])


def run(args, case, seed, samples, warmup, label, procs=1, split=1, initial=None, ref=None, force=False):
    work = pathlib.Path(args.workdir) / label
    base.require(not work.exists(), "refusing to overwrite " + str(work))
    base.write_e2e_input(case, work, seed=seed, samples=samples, warmup=warmup,
                         nsplit=split, init_nelec=initial)
    env = {"MVMC_GC_LANCZOS_TEST_SUMMARY": "1", "MVMC_GC_B0_LANCZOS_TRACE": "1"}
    if force:
        env["MVMC_GC_LANCZOS_TEST_FAULT"] = "refine_final"
    done = base.run_vmc(args.vmc, work, procs, env=env, timeout=600)
    chains = base.chain_count(procs, split)
    result = base.require_success(label, done, work, samples, chains)
    summary = base.parse_summary(label, done.stdout, chains, case["L"])
    ref = reference(case) if ref is None else ref
    evidence = replay(case, work, result, summary, samples, chains, ref)
    evidence["refinement"] = re.findall(r"GC Lanczos refinement:.*", done.stdout)
    evidence["coefficients"] = [[z.real, z.imag] for z in summary["c"]]
    (work / "replay.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print("PASS %s samples=%d nodes=%d ESS=%.0f" %
          (label, samples * chains, evidence["node_samples"], evidence["ess"]), flush=True)
    return evidence


def replay_mode(args):
    for name in (args.case,) if args.case else ("minimal",) + STAT_CASES + ("non_su2",):
        case = fixture(name)
        row = run(args, case, 17, 8192, 512, name)
        if name == "non_su2":
            base.require(row["node_samples"] > 0, "no samples on base nodes")


def sample_mode(args):
    case = fixture("multi")
    case["gauss"] = 16
    work = pathlib.Path(args.workdir) / "sample"
    shutil.copytree(args.sample, work)
    done = base.run_vmc(args.vmc, work, env={
        "MVMC_GC_LANCZOS_TEST_SUMMARY": "1", "MVMC_GC_B0_LANCZOS_TRACE": "1"}, timeout=600)
    result = base.require_success("sample", done, work, 8192, 1)
    summary = base.parse_summary("sample", done.stdout, 1, case["L"])
    evidence = replay(case, work, result, summary, 8192, 1, reference(case))
    base.require(all(result["number"]["sectors"][n] > 0 for n in (4, 6)),
                 "sample did not visit both N sectors")
    (work / "replay.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print("PASS public S=2 sample: all 8192 snapshots and N=4,6", flush=True)


def statistical_comparisons(name, group):
    # This shared helper takes an upper-tail probability, not a CDF value.
    # Cauchy has t_{tail=1/4,df=1}=1: catch probability-convention mistakes.
    close("Student t upper-tail convention", student_t_quantile(.25, 1), 1.0, 1e-12)
    critical = student_t_quantile(.01 / (2 * 66), len(CONFIRM_SEEDS) - 1)
    base.require(2 < critical < 100, "invalid simultaneous critical value")
    comparisons = []
    base.require(len(group) == len(CONFIRM_SEEDS), "missing confirmation seeds")
    for key in group[0]["observables"]:
        values = np.array([r["observables"][key] for r in group])
        residual = values[:, 0] - values[:, 1]
        error = abs(residual.mean())
        se = residual.std(ddof=1) / math.sqrt(len(residual))
        bound = critical * se + 1e-8 * (1 + np.max(abs(values[:, 1])))
        comparisons.append(dict(case=name, observable=key, error=float(error),
                                se=float(se), bound=float(bound), pass_=bool(error <= bound)))
    return critical, comparisons


def statistics_mode(args):
    phase = args.mode
    seeds = PILOT_SEEDS if phase == "pilot" else CONFIRM_SEEDS
    samples, warmup = (4096, 1024) if phase == "pilot" else (16384, 2048)
    rows, comparisons = [], []
    critical = student_t_quantile(.01 / (2 * 66), len(CONFIRM_SEEDS) - 1)
    for name in (args.case,) if args.case else STAT_CASES:
        case = fixture(name)
        ref = reference(case)
        group = []
        for seed in seeds:
            row = run(args, case, seed, samples, warmup, "%s-%d" % (name, seed), ref=ref)
            base.require(row["ess"] >= 128, "ESS below 128")
            row.pop("records")
            group.append(row)
            rows.append(dict(case=name, seed=seed, **row))
            (pathlib.Path(args.workdir) / "statistics.json").write_text(
                json.dumps(dict(phase=phase, rows=rows, comparisons=comparisons), indent=2) + "\n")
        if phase == "confirm":
            critical, checks = statistical_comparisons(name, group)
            comparisons.extend(checks)
    (pathlib.Path(args.workdir) / "statistics.json").write_text(
        json.dumps(dict(phase=phase, rows=rows, comparisons=comparisons,
                        critical=critical, alpha=.01, family_comparisons=66), indent=2) + "\n")
    base.require(all(c["pass_"] for c in comparisons), "statistical comparison failed")


def mpi_mode(args):
    # Same chain seeds and QP partitions, including G=2 / np3 empty QP rank.
    case = fixture("general")
    case["gauss"] = 2
    ref = reference(case)
    first = run(args, case, 17, 8192, 512, "serial", ref=ref)
    for np_ in (2, 3):
        row = run(args, case, 17, 8192, 512, "np%d" % np_, procs=np_, split=np_, ref=ref)
        base.require(row["records"] == first["records"], "QP partition changed sample trace")
        close("QP coefficients", row["coefficients"], first["coefficients"], 2e-10)
    two = run(args, case, 29, 8192, 512, "two-chains", procs=2, split=1, ref=ref)
    split = run(args, case, 29, 8192, 512, "two-by-two", procs=4, split=2, ref=ref)
    base.require(two["records"] == split["records"], "two-chain QP partition changed trace")


def oracle_mode(args):
    for name in ("maximal", "multi", "upper", "pbc_k1", "apbc_k1", "non_su2"):
        case = fixture(name)
        work = pathlib.Path(args.workdir) / name
        work.mkdir(parents=True, exist_ok=False)
        code, data = base.run_probe(args.probe, case, work, 1)
        base.require(code == 0, "probe exit %d" % code)
        base.compare_case(case, data)
        print("PASS full-space probe " + name, flush=True)


def refined_mode(args):
    for name in ("general", "multi"):
        case = fixture(name)
        if name == "general":
            case["gauss"] = 2
        ref = reference(case)
        control = run(args,case,17,8192,512,name+"-ordinary",ref=ref)
        layouts = ((1,1),) if args.no_mpi or name == "multi" else ((1,1),(2,2),(3,3))
        for procs,split in layouts:
            row=run(args,case,17,8192,512,name+"-refined-np%d"%procs,
                    procs=procs,split=split,ref=ref,force=True)
            accepted=sum(int(re.search(r" accepted=(\d+)",line).group(1)) for line in row["refinement"])
            base.require(accepted>0,"full run did not accept a refined step")
            base.require(row["records"]==control["records"],"refinement changed same-uniform trajectory")
            print("PASS refined accepts %s np%d: %d"%(name,procs,accepted),flush=True)
        if name == "general" and not args.no_mpi:
            control = run(args,case,29,8192,512,"two-chain-ordinary",procs=2,split=1,ref=ref)
            row = run(args,case,29,8192,512,"two-chain-refined",procs=4,split=2,ref=ref,force=True)
            accepts = [int(re.search(r" accepted=(\d+)",line).group(1)) for line in row["refinement"]]
            base.require(len(accepts)==2 and min(accepts)>0, "missing refined accept on a chain")
            base.require(row["records"]==control["records"], "two-chain refined trajectory changed")
            print("PASS two-chain refined accepts: %s" % accepts,flush=True)


def convergence_mode(args):
    rows = []
    for name in ("minimal", "general", "triplet", "maximal", "multi", "upper", "pbc_k1", "apbc_k1", "non_su2"):
        case = fixture(name)
        L, spin = case["L"], case["S"]
        raw = {m: so.qfloat(v) for m, v in so.raw_exact_spin(case["F"], spin).items()}
        steps = [r for r, _ in case["translations"]]
        mappings = gl.chain_mappings(L, case["boundary"], steps)
        exact = {m: 0j for m in raw}
        for r, w in case["translations"]:
            part = gl.translated_wave(raw, mappings, r, L)
            for m in exact:
                exact[m] += w * part[m]
        basis = gl.domain_basis(L, "balanced")
        target = np.array([exact[m] for m in basis])
        norm = np.linalg.norm(target)
        base.require(norm > 1e-12, "zero exact-spin reference")
        h = gl.hamiltonian(case, basis)
        htarget = h @ target
        h_residual = so.spin_residual(dict(zip(basis, htarget)), L, spin)
        exact_gevp = gl.gevp(*gl.krylov(target,h)[1:])
        base.require(exact_gevp is not None, "missing exact GEVP")
        phi = exact_gevp["c"][0]*target + exact_gevp["c"][1]*htarget
        phi_residual = so.spin_residual(dict(zip(basis,phi)),L,spin)
        if name != "non_su2":
            base.require(max(h_residual,phi_residual) < 1e-10, "SU(2) H/Lanczos lost exact spin")
        else:
            base.require(min(h_residual,phi_residual) > .1, "non-SU(2) counterexample vanished")
            base.require(any(abs(htarget[i]) > 1e-8 and gl.support_zero(case, m)
                             for i, m in enumerate(basis)), "missing H Psi on base nodes")
        for gauss in (2, 3, 4, 8, 16, 32, 64):
            _, wave = gl.reference_wave(dict(case, gauss=gauss))
            error = np.linalg.norm(wave - target) / norm
            residual = so.spin_residual(dict(zip(basis, wave)), L, spin)
            p, t = abs(wave)**2, abs(target)**2
            p, t = p / p.sum(), t / t.sum()
            sector_error = max(abs(sum(p[i]-t[i] for i,m in enumerate(basis)
                                       if gl.popcount(m)==n)) for n in range(2*L+1))
            rows.append(dict(case=name, gauss=gauss, wave_error=float(error),
                             spin_residual=float(residual), sector_error=float(sector_error),
                             h_spin_residual=float(h_residual), final_spin_residual=float(phi_residual)))
            if gauss >= 32:
                base.require(max(error, residual, sector_error) <= 1e-10,
                             "exact-spin convergence failed: %s G%d" % (name,gauss))
        print("PASS spin convergence " + name, flush=True)
    (pathlib.Path(args.workdir) / "convergence.json").write_text(json.dumps(rows,indent=2)+"\n")


def initial_mode(args):
    case = fixture("multi")
    ref = reference(case)
    for n in (4, 6):
        row = run(args, case, 29, 8192, 512, "initial-N%d" % n, initial=n, ref=ref)
        lines = base.initialization_lines((pathlib.Path(args.workdir)/("initial-N%d"%n)/"stdout.txt").read_text())
        base.require((0,n) in lines and lines[(0,n)]["selected"] == "1", "explicit N not kept")
        base.require({4,6} <= {gl.popcount(record[2]) for record in row["records"]},
                     "multi-N chain did not visit both sectors")


def fault_mode(args):
    case = fixture("general")
    case["gauss"] = 2
    for procs,split in (((1,1),) if args.no_mpi else ((1,1),(2,2),(3,3),(4,2))):
        faults = ("bootstrap_exhausted", "final_numeric")
        if procs > 1:
            faults += ("allocation:rank=1", "nonfinite:rank=1")
        for fault in faults:
            work = pathlib.Path(args.workdir)/("np%d-%s" % (procs,fault.replace(':','-')))
            base.require(not work.exists(), "existing fault directory")
            base.write_e2e_input(case, work, seed=17, samples=8192, warmup=8,
                                 nsplit=split, init_nelec=None)
            env = {"MVMC_GC_LANCZOS_TEST_FAULT":fault, "MVMC_GC_LANCZOS_TEST_SUMMARY":"1"}
            done=base.run_vmc(args.vmc,work,procs,env=env,timeout=120)
            needles={"bootstrap_exhausted":["stage=final-bootstrap","exhausted=1"],
                     "final_numeric":["stage=final-sample","amplitude_state=numeric_zero"],
                     "allocation:rank=1":["allocation failure"],
                     "nonfinite:rank=1":["stage=coefficient-initial","nonfinite value"]}[fault]
            base.require(done.returncode == 1, "unexpected failure exit")
            base.require_failure(str(work.name),done,work,needles)
            base.require("GC_LANCZOS_SUMMARY" not in done.stdout, "failure published summary")
            base.require(not (work/'output/zvo_gc_pl_number_001.dat').exists(), "failure published N")
            base.require(not re.search(r"AddressSanitizer|UndefinedBehaviorSanitizer|runtime error:",done.stdout+done.stderr), "sanitizer report")
            print("PASS finite-S fault " + work.name,flush=True)


def contract_mode(args):
    from unittest.mock import patch
    import copy
    import subprocess
    case = fixture("general")
    run(args,case,17,8192,512,"control")
    work=pathlib.Path(args.workdir)/"control"
    done=subprocess.CompletedProcess([],0,(work/"stdout.txt").read_text(),(work/"stderr.txt").read_text())
    result=base.require_success("contract",done,work,8192,1)
    summary=base.parse_summary("contract",done.stdout,1,case["L"])
    ref=reference(case)
    for name in ("public_coefficients", "bootstrap"):
        altered_result, altered_summary = copy.deepcopy(result), copy.deepcopy(summary)
        if name == "public_coefficients":
            altered_result["number"]["c"][0] += 1
        else:
            altered_summary["bootstrap"][0] = (1, 0, 0, 0, 1)
        try:
            require_run_metadata(altered_result, altered_summary)
        except base.CheckFailure:
            pass
        else:
            raise base.CheckFailure("undetected metadata mutation: " + name)
    trace=work/"gc-b0-lanczos.chain0"
    lines=trace.read_text().splitlines()
    index=next(i for i,line in enumerate(lines) if line.startswith("F1 "))
    original=pathlib.Path.read_text
    for name in ("missing", "duplicate", "configuration", "energy", "q0", "zero_phi", "generation"):
        altered=lines.copy()
        fields=altered[index].split()
        if name=="missing":
            altered.pop(index)
        elif name=="duplicate":
            fields[1]="1"
        elif name=="configuration":
            fields[5]="ffffffffffffffff"
        elif name=="energy":
            fields[6]=str(float(fields[6])+1)
        elif name=="q0":
            fields[9]=str(float(fields[9])+1)
        elif name=="zero_phi":
            fields[20]="1"
        else:
            fields[2]="99999999"
        if name!="missing":
            altered[index]=" ".join(fields)
        content="\n".join(altered)+"\n"
        def read(path,*a,**kw):
            return content if path==trace else original(path,*a,**kw)
        try:
            with patch.object(pathlib.Path,"read_text",read):
                replay(case,work,result,summary,8192,1,ref)
        except base.CheckFailure:
            pass
        else:
            raise base.CheckFailure("undetected trace mutation: "+name)
    group=[dict(observables={"E":[1+(.01 if i%2 else -.01),1.]}) for i in range(8)]
    _, checks=statistical_comparisons("synthetic",group)
    base.require(all(c["pass_"] for c in checks),"centered control rejected")
    shifted=copy.deepcopy(group)
    for row in shifted:
        row["observables"]["E"][0]+=1
    _,checks=statistical_comparisons("biased",shifted)
    base.require(not checks[0]["pass_"],"statistical bias mutation not rejected")
    print("PASS seven trace mutations, two metadata mutations and statistical bias mutation",flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('replay', 'pilot', 'confirm', 'mpi', 'oracle', 'convergence', 'initial', 'fault', 'refined', 'contract', 'sample'), required=True)
    parser.add_argument('--vmc')
    parser.add_argument('--probe')
    parser.add_argument('--sample')
    parser.add_argument('--case')
    parser.add_argument('--no-mpi', action='store_true')
    parser.add_argument('--unique-run', action='store_true')
    parser.add_argument('--workdir', required=True)
    args = parser.parse_args()
    pathlib.Path(args.workdir).mkdir(parents=True, exist_ok=True)
    if args.unique_run:
        args.workdir = tempfile.mkdtemp(prefix='run-', dir=args.workdir)
    previous_traces = set(pathlib.Path(args.workdir).glob("*/gc-b0-lanczos.chain*"))
    if args.mode in ('pilot', 'confirm'):
        statistics_mode(args)
    else:
        globals()[args.mode + '_mode'](args)
    compress_passed_traces(args.workdir, previous_traces)
    print('Finite-S Lanczos ' + args.mode + ': PASS')


if __name__ == '__main__':
    main()
