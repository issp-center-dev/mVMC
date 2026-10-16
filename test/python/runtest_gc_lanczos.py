#!/usr/bin/env python3
"""GC first-step stabilized Lanczos checks.

modes:
  proposal  : exact Fraction cross-check of the symmetric GC proposals
              enumerated by gc_krylov_proposal_unit --dump.
  oracle    : enumeration probe (production GC bridge, bounded engine,
              measurement, GEVP) against the independent Fock oracle
              (grandcanonical_lanczos_oracle) for F2/G2/F4, the six
              pairing/projection combinations x PBC/APBC, boundary, node,
              singular, rank and 2^+-300 scale fixtures.
  input     : expert-input mutations (non-zero exit, named diagnostic, no
              result row) and accepted inputs, serial or --np 2.
  production: F2/G2/F2+spin+momentum expert inputs through vmc.out, seeds
              17/29/43/71, compared with the oracle for the frozen c.
  mpi       : --np 2 (QP partition == serial, two chains) or --np 3
              (chains 2+1 with an empty QP rank == two single-rank chains).
  failure   : MVMC_GC_LANCZOS_TEST_FAULT contracts and the real singular-F
              classification.
  sample    : samples/GrandCanonical/FirstLanczos, copied to --workdir.

Tolerances are the frozen plan values (fixture-policy):
  q        : |a - r| <= 3e-12 + 3e-11 |r| on one common scale
  scale    : phase <= 3e-12, log_abs <= 3e-11
  S/K/T    : 2e-11 + 2e-10 |r| after the oracle's common normalization
  GEVP     : E/variance <= 2e-10 (1 + |r|), Phi norm diff <= 1e-9
"""
from __future__ import print_function

import argparse
import cmath
import math
import os
import pathlib
import subprocess
import sys

import numpy as np

import grandcanonical_lanczos_oracle as gl
import grandcanonical_spin_oracle as so

UNIT_PHASES = (1.0 + 0.0j, -1.0j, -1.0 + 0.0j, 1.0j)
F2 = [[1.0, (1 + 2j) / 8], [(3 - 1j) / 16, 0.75]]
NODE = [[1.0, 0.0], [1.0, 1.0]]
SINGULAR = [[1.0, 1.0], [1.0, 1.0]]
NEAR_SINGULAR = [[1.0, 1.0], [1.0, 1.0 + 2.0 ** -40]]


class CheckFailure(Exception):
    pass


def require(condition, message):
    if not condition:
        raise CheckFailure(message)


def f4_matrix():
    return [[(1.0 if i == j else 0.0) + (i + 1) / 32.0 + (j + 1) * 1j / 64.0
             for j in range(4)] for i in range(4)]


def general_matrix(F):
    nsite = len(F)
    A = so.antiparallel_matrix(F)
    A[0][1], A[1][0] = (1 + 1j) / 4, -(1 + 1j) / 4
    A[nsite][nsite + 1] = (2 - 1j) / 8
    A[nsite + 1][nsite] = -(2 - 1j) / 8
    return A


def momentum(nsite, k):
    """Exact exp(-2 pi i k r / L) / L for L in {1, 2, 4}."""
    weights = []
    for r in range(nsite):
        numerator = 4 * k * r
        if numerator % nsite:
            raise ValueError("phase not tabulated")
        weights.append((r, UNIT_PHASES[(numerator // nsite) % 4] / nsite))
    return weights


def model_rows(nsite, bonds, boundary, t=1.0, phase=1.0, mu=0.25, U=1.5,
               V=None, J=None, scale=1.0):
    trans = []
    for i, j in bonds:
        sign = -1.0 if (boundary == "apbc" and j < i) else 1.0
        value = scale * t * phase * sign
        for spin in range(2):
            trans.append((i, spin, j, spin, value))
            trans.append((j, spin, i, spin, complex(value).conjugate()))
    for site in range(nsite):
        if mu != 0.0:
            for spin in range(2):
                trans.append((site, spin, site, spin, scale * mu))
    rows = {"trans": trans,
            "intra": [(s, scale * U) for s in range(nsite)] if U else [],
            "inter": [(i, j, scale * V) for i, j in bonds] if V else [],
            "hund": [(i, j, scale * J) for i, j in bonds] if J else []}
    return rows


def ring_bonds(nsite):
    if nsite == 2:
        return [(0, 1)]
    return [(i, (i + 1) % nsite) for i in range(nsite)]


def fixture(name, nsite, pairing, matrix, boundary="pbc", translations=None,
            gauss=1, rows=None, correlation=False, raw_weights=None):
    item = {"name": name, "L": nsite, "pairing": pairing,
            "domain": "balanced" if pairing == "antiparallel" else "even",
            "boundary": boundary,
            "translations": translations or [(0, 1.0)], "gauss": gauss}
    if pairing == "antiparallel":
        item["F"] = matrix
    else:
        item["A"] = matrix
    if raw_weights is not None:
        item["raw_translation_weights"] = [(i, w) for i, w in
                                           enumerate(raw_weights)]
    item.update(rows if rows is not None else
                model_rows(nsite, ring_bonds(nsite), boundary))
    if correlation:
        item["ngutzwiller"] = 1
        item["gutzwiller_idx"] = [0] * nsite
        item["jastrow_idx"] = [[-1 if i == j else 0 for j in range(nsite)]
                               for i in range(nsite)]
        item["projection_parameters"] = [-0.375 + 0.75j, 0.1875 - 0.5j]
    return item


def oracle_cases():
    cases = []
    f2_rows = lambda boundary, phase=1.0: model_rows(
        2, [(0, 1)], boundary, phase=phase, V=0.125, J=0.0625)
    g2 = general_matrix(F2)
    for boundary in ("pbc", "apbc"):
        tag = boundary
        cases.append(fixture("F2_none_" + tag, 2, "antiparallel", F2, boundary,
                             rows=f2_rows(boundary)))
        for k in (0, 1):
            cases.append(fixture("F2_MP_k%d_%s" % (k, tag), 2, "antiparallel",
                                 F2, boundary, momentum(2, k),
                                 rows=f2_rows(boundary)))
        cases.append(fixture("F2_SP_G2_" + tag, 2, "antiparallel", F2,
                             boundary, gauss=2, rows=f2_rows(boundary)))
        cases.append(fixture("F2_SMP_G2_k0_" + tag, 2, "antiparallel", F2,
                             boundary, momentum(2, 0), gauss=2,
                             rows=f2_rows(boundary)))
        cases.append(fixture("G2_none_" + tag, 2, "general", g2, boundary,
                             rows=f2_rows(boundary)))
        for k in (0, 1):
            cases.append(fixture("G2_MP_k%d_%s" % (k, tag), 2, "general", g2,
                                 boundary, momentum(2, k),
                                 rows=f2_rows(boundary)))
        F4 = f4_matrix()
        G4 = general_matrix(F4)
        cases.append(fixture("F4_none_" + tag, 4, "antiparallel", F4,
                             boundary))
        for k in (0, 1):
            cases.append(fixture("F4_MP_k%d_%s" % (k, tag), 4,
                                 "antiparallel", F4, boundary,
                                 momentum(4, k)))
        for gauss in (2, 4):
            cases.append(fixture("F4_SP_G%d_%s" % (gauss, tag), 4,
                                 "antiparallel", F4, boundary, gauss=gauss))
        cases.append(fixture("F4_SMP_G2_k0_" + tag, 4, "antiparallel", F4,
                             boundary, momentum(4, 0), gauss=2))
        cases.append(fixture("F4_SMP_G4_k1_" + tag, 4, "antiparallel", F4,
                             boundary, momentum(4, 1), gauss=4))
        cases.append(fixture("G4_none_" + tag, 4, "general", G4, boundary))
        for k in (0, 1):
            cases.append(fixture("G4_MP_k%d_%s" % (k, tag), 4, "general", G4,
                                 boundary, momentum(4, k)))
    cases.append(fixture("F2_complex_hopping", 2, "antiparallel", F2,
                         rows=f2_rows("pbc", phase=(3 + 4j) / 5)))
    cases.append(fixture("F2_correlation", 2, "antiparallel", F2,
                         rows=f2_rows("pbc"), correlation=True))
    cases.append(fixture("G2_correlation_MP_k1", 2, "general", g2, "apbc",
                         momentum(2, 1), rows=f2_rows("apbc"),
                         correlation=True))
    boundary_case = fixture("F2_boundary_exact", 2, "antiparallel", F2,
                            "pbc", [(0, 1.0), (1, -1.0)],
                            rows=f2_rows("pbc"))
    boundary_case["expect_boundary_zero"] = True
    cases.append(boundary_case)
    near = fixture("F2_boundary_near", 2, "antiparallel", F2, "pbc",
                   [(0, 1.0), (1, 2.0 ** -53), (0, -1.0)],
                   rows=f2_rows("pbc"))
    near["expect_boundary_zero"] = False
    # 1 + 2^-53 - 1 cancels to rounding level: only the boundary state is
    # compared, values are rounding dominated in both implementations.
    near["compare_values"] = False
    near["skip_gevp"] = True
    cases.append(near)
    node = fixture("Node", 2, "antiparallel", NODE,
                   rows=model_rows(2, [(0, 1)], "pbc", mu=0.0, U=0.0))
    node["expect_node"] = True
    cases.append(node)
    zero = fixture("QP_zero", 2, "antiparallel", [[0.0, 0.0], [0.0, 0.0]],
                   rows=f2_rows("pbc"))
    zero["expect_gevp_failure"] = True  # vacuum-only Psi has H Psi == 0
    cases.append(zero)
    singular = fixture("Singular", 2, "antiparallel", SINGULAR,
                       rows=f2_rows("pbc"))
    singular["classify_full"] = True
    cases.append(singular)
    near_singular = fixture("NearSingular", 2, "antiparallel", NEAR_SINGULAR,
                            rows=f2_rows("pbc"))
    near_singular["classify_full"] = True
    near_singular["near_singular"] = True
    cases.append(near_singular)
    eigen = fixture("Rank1_eigenstate", 2, "antiparallel", F2, "pbc",
                    momentum(2, 1),
                    rows=model_rows(2, [(0, 1)], "pbc", t=0.0, mu=0.25, U=0.0))
    eigen["expect_rank"] = 1
    cases.append(eigen)
    # Decision A: q1 == 0 (zero H, H Psi == 0) is a GEVP failure, not rank 1.
    rank1 = fixture("ZeroH_q1zero", 2, "antiparallel", F2,
                    rows={"trans": [], "intra": [], "inter": [], "hund": []})
    rank1["expect_gevp_failure"] = True
    cases.append(rank1)
    rank0 = fixture("Rank0_allzero", 2, "antiparallel",
                    [[0.0, 0.0], [0.0, 0.0]], "pbc", [(0, 1.0), (1, -1.0)],
                    rows=f2_rows("pbc"))
    rank0["expect_rank"] = 0
    cases.append(rank0)
    for sign in (1, -1):
        F4 = f4_matrix()
        scaled = fixture("F4_scale_%s300" % ("+" if sign > 0 else "-"), 4,
                         "antiparallel",
                         [[v * 2.0 ** (300 * sign) for v in row] for row in F4],
                         "pbc", momentum(4, 1), gauss=2)
        scaled["scale_shift"] = 300.0 * sign * math.log(2.0)
        scaled["unscaled_F"] = F4
        cases.append(scaled)
    return cases


def fmt(value):
    return "%.17e" % value


def write_probe_input(case, path):
    nsite = case["L"]
    slaters, weights = gl.production_qp(case)
    raw = [w for _, w in case.get("raw_translation_weights",
                                   case["translations"])]
    lines = ["GC_LANCZOS_PROBE_V1",
             "%d %d %d %d" % (nsite, 1 if case["domain"] == "balanced" else 2,
                              len(slaters), len(raw))]
    lines += ["%s %s" % (fmt(complex(w).real), fmt(complex(w).imag))
              for w in raw]
    lines += ["%s %s" % (fmt(w.real), fmt(w.imag)) for w in weights]
    for slater in slaters:
        for row in slater:
            lines.append(" ".join("%s %s" % (fmt(complex(v).real),
                                             fmt(complex(v).imag))
                                  for v in row))
    lines.append(str(len(case["trans"])))
    for out_site, out_spin, in_site, in_spin, value in case["trans"]:
        value = complex(value)
        lines.append("%d %d %d %d %s %s" % (out_site, out_spin, in_site,
                                            in_spin, fmt(value.real),
                                            fmt(value.imag)))
    for key in ("intra",):
        lines.append(str(len(case.get(key, []))))
        lines += ["%d %s" % (s, fmt(v)) for s, v in case.get(key, [])]
    for key in ("inter", "hund"):
        lines.append(str(len(case.get(key, []))))
        lines += ["%d %d %s" % (i, j, fmt(v)) for i, j, v in case.get(key, [])]
    params = case.get("projection_parameters", [])
    gutz = case.get("gutzwiller_idx")
    table = case.get("jastrow_idx")
    lines.append("%d %d" % (case.get("ngutzwiller", 0) if gutz else 0,
                            1 if table else 0))
    if gutz:
        lines.append(" ".join(str(v) for v in gutz))
    if table:
        lines.append(" ".join(str(v) for row in table for v in row))
    lines.append(str(len(params)))
    lines += ["%s %s" % (fmt(complex(p).real), fmt(complex(p).imag))
              for p in params]
    lines.append("END")
    pathlib.Path(path).write_text("\n".join(lines) + "\n")


def run_probe(probe, case, workdir, nproc):
    work = pathlib.Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    source = work / (case["name"] + ".in")
    target = work / (case["name"] + ".tsv")
    write_probe_input(case, source)
    command = [probe, str(source), str(target)]
    if nproc > 1:
        command = ["mpirun", "-np", str(nproc)] + command
    done = subprocess.run(command, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, timeout=180)
    if not target.exists():
        raise CheckFailure("%s: probe produced no output (exit %d)\n%s" %
                           (case["name"], done.returncode, done.stdout))
    return done.returncode, parse_probe(target)


def parse_probe(path):
    result = {"Q": {}, "M": {}, "GEVP": None, "STATUS": None, "BRIDGE": None,
              "MEASUREMENT": None}
    for line in pathlib.Path(path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0] == "Q":
            mask = int(fields[1])
            values = []
            for order in range(3):
                base = 4 + 4 * order
                if len(fields) <= base:
                    break
                values.append((fields[base], complex(float(fields[base + 1]),
                                                     float(fields[base + 2])),
                               float(fields[base + 3])))
            result["Q"][mask] = (fields[3], values)
        elif fields[0] == "M":
            numbers = [float(v) for v in fields[2:]]
            result["M"][int(fields[1])] = [
                complex(numbers[2 * i], numbers[2 * i + 1]) for i in range(4)]
        elif fields[0] == "GEVP":
            result["GEVP"] = {
                "status": fields[1], "rank": int(fields[2]),
                "energy": float(fields[3]), "variance": float(fields[4]),
                "c": np.array([complex(float(fields[5]), float(fields[6])),
                               complex(float(fields[7]), float(fields[8]))]),
                "residual": float(fields[9])}
        elif fields[0] in ("STATUS", "BRIDGE", "MEASUREMENT"):
            result[fields[0]] = fields[1]
    return result


def scaled_value(entry, common):
    state, phase, log_abs = entry
    if state in ("EXACT_ZERO", "NUMERIC_ZERO"):
        return 0.0j
    require(state == "FINITE", "nonfinite probe value")
    return phase * math.exp(log_abs - common)


def compare_case(case, probe_result):
    name = case["name"]
    nsite = case["L"]
    require(probe_result["BRIDGE"] == "ok" and
            probe_result["STATUS"] == "ok",
            "%s: probe status %s/%s" % (name, probe_result["BRIDGE"],
                                         probe_result["STATUS"]))
    if "scale_shift" not in case:
        require(probe_result["MEASUREMENT"] == "ok",
                "%s: probe measurement failed" % name)
    if "scale_shift" in case:
        unscaled = dict(case, F=case["unscaled_F"])
        basis, psi = gl.reference_wave(unscaled)
        q, s, k, t = gl.krylov(psi, gl.hamiltonian(unscaled, basis))
        sector_scale = {}
        for order in range(3):
            for position, mask in enumerate(basis):
                key = (order, gl.popcount(mask))
                sector_scale[key] = max(sector_scale.get(key, 0.0),
                                        abs(q[order][position]))
        for position, mask in enumerate(basis):
            status, values = probe_result["Q"][mask]
            require(status == "ok", "%s: mask %d status" % (name, mask))
            shift = case["scale_shift"] * gl.popcount(mask) / 2.0
            for order in range(3):
                ref = q[order][position]
                # Symmetry-cancelled entries are rounding noise in the oracle.
                if abs(ref) <= 1e-10 * sector_scale[(order,
                                                     gl.popcount(mask))]:
                    continue
                state, phase, log_abs = values[order]
                require(state == "FINITE", "%s: scaled state" % name)
                require(abs(phase - ref / abs(ref)) <= 3e-12 and
                        abs(log_abs - (math.log(abs(ref)) + shift)) <= 3e-11,
                        "%s: scaled phase/log mask %d q%d" %
                        (name, mask, order))
        return "scale"
    basis, psi = gl.reference_wave(case)
    q, s, k, t = gl.krylov(psi, gl.hamiltonian(case, basis))
    require(sorted(probe_result["Q"]) == basis,
            "%s: configuration set differs" % name)
    logs = [math.log(abs(v)) for order in range(3) for v in q[order]
            if abs(v) > 0.0]
    common = max(logs) if logs else 0.0
    worst = 0.0
    for position, mask in enumerate(basis):
        status, values = probe_result["Q"][mask]
        require(status == "ok", "%s: mask %d status %s" % (name, mask, status))
        if not case.get("compare_values", True):
            continue
        for order in range(3):
            actual = scaled_value(values[order], common)
            ref = q[order][position] * math.exp(-common)
            error = abs(actual - ref)
            worst = max(worst, error / (3e-12 + 3e-11 * abs(ref)))
            require(error <= 3e-12 + 3e-11 * abs(ref),
                    "%s: mask %d q%d actual %r ref %r" %
                    (name, mask, order, actual, ref))
    if case.get("expect_boundary_zero") is not None:
        for mask in (0, (1 << (2 * nsite)) - 1):
            state = probe_result["Q"][mask][1][0][0]
            if case["expect_boundary_zero"]:
                require(state == "EXACT_ZERO",
                        "%s: boundary mask %d is %s" % (name, mask, state))
            else:
                require(state != "EXACT_ZERO",
                        "%s: near boundary overridden" % name)
    if case.get("expect_node"):
        node = probe_result["Q"][0x9][1]
        require(node[0][0] == "EXACT_ZERO", "%s: node q0 state" % name)
        q1 = scaled_value(node[1], 0.0)
        require(abs(q1 - (-2.0)) <= 3e-12, "%s: node q1 %r" % (name, q1))
    classification = None
    if case.get("classify_full"):
        full = (1 << (2 * nsite)) - 1
        classification = probe_result["Q"][full][1][0][0]
        if case.get("near_singular"):
            ref = q[0][basis.index(full)]
            actual = scaled_value(probe_result["Q"][full][1][0], 0.0)
            require(classification == "FINITE" and
                    abs(actual - ref) <= 1e-3 * abs(ref),
                    "%s: near-singular full filling %s %r vs %r" %
                    (name, classification, actual, ref))
    if case.get("skip_gevp"):
        return "values-skipped"
    # S/K/T after the oracle's common normalization.
    gevp_actual = probe_result["GEVP"]
    oracle = gl.gevp(s, k, t)
    expected_rank = case.get("expect_rank")
    if case.get("expect_gevp_failure"):
        require(s[1, 1] == 0.0, "%s: oracle q1 must vanish" % name)
        require(gevp_actual is not None and gevp_actual["status"] != "ok" and
                gevp_actual["rank"] == 0,
                "%s: q1 == 0 must fail in the GEVP (%r)" % (name, gevp_actual))
        return "gevp-fails(q1=0)"
    if expected_rank == 0:
        require(oracle is None, "%s: oracle rank" % name)
        require(gevp_actual is not None and gevp_actual["status"] != "ok" and
                gevp_actual["rank"] == 0,
                "%s: rank-0 GEVP must fail (%r)" % (name, gevp_actual))
        return "rank0"
    norm = s[0, 0].real
    require(norm > 0.0, "%s: oracle normalization" % name)
    for index, (rs, rk, rka, rt) in enumerate(zip(gl.packed(s), gl.packed(k),
                                                   gl.packed(k.T), gl.packed(t))):
        actual = probe_result["M"][index]
        for got, ref in zip(actual, (rs, rk, rka, rt)):
            got, ref = got / norm, ref / norm
            require(abs(got - ref) <= 2e-11 + 2e-10 * abs(ref),
                    "%s: matrix entry %d %r vs %r" % (name, index, got, ref))
    require(gevp_actual is not None and gevp_actual["status"] == "ok",
            "%s: GEVP status %r" % (name, gevp_actual))
    require(abs(gevp_actual["energy"] - oracle["energy"]) <=
            2e-10 * (1 + abs(oracle["energy"])),
            "%s: GEVP energy %r vs %r" % (name, gevp_actual["energy"],
                                          oracle["energy"]))
    require(abs(gevp_actual["variance"] - oracle["variance"]) <=
            2e-10 * (1 + abs(oracle["variance"])),
            "%s: GEVP variance %r vs %r" % (name, gevp_actual["variance"],
                                            oracle["variance"]))
    if expected_rank is not None:
        require(gevp_actual["rank"] == expected_rank and
                oracle["rank"] == expected_rank,
                "%s: rank %d/%d" % (name, gevp_actual["rank"], oracle["rank"]))
    phi_actual = gevp_actual["c"][0] * q[0] + gevp_actual["c"][1] * q[1]
    phi_oracle = oracle["c"][0] * q[0] + oracle["c"][1] * q[1]
    na, no = np.linalg.norm(phi_actual), np.linalg.norm(phi_oracle)
    require(na > 0 and no > 0, "%s: zero corrected state" % name)
    phi_actual, phi_oracle = phi_actual / na, phi_oracle / no
    overlap = np.vdot(phi_oracle, phi_actual)
    phi_actual *= abs(overlap) / overlap if overlap != 0 else 1.0
    if oracle["gap"] > 1e-8 * (1 + abs(oracle["energy"])):
        require(np.linalg.norm(phi_actual - phi_oracle) <= 1e-9,
                "%s: corrected state differs" % name)
    return "%.2f%s" % (worst, "" if classification is None
                       else " full=" + classification)


def oracle_mode(args):
    gl.self_test()
    failures = 0
    for case in oracle_cases():
        if args.case and args.case != case["name"]:
            continue
        try:
            code, result = run_probe(args.probe, case, args.workdir, args.np)
            summary = compare_case(case, result)
            require(code == 0 or case.get("expect_rank") == 0,
                    "%s: probe exit %d" % (case["name"], code))
            print("PASS %-26s %s" % (case["name"], summary))
        except CheckFailure as error:
            failures += 1
            print("FAIL %s" % error)
    if failures:
        print("GCLanczos oracle: %d failure(s)" % failures)
        return 1
    print("GCLanczos oracle: PASS")
    return 0


# --- proposal: exact Fraction oracle of the symmetric GC proposals --------
from collections import defaultdict
from fractions import Fraction


def enumerate_draws(draw_function):
    """All bounded draw sequences of draw_function(choose) -> mask."""
    result = defaultdict(Fraction)
    path = []
    while True:
        used = []

        def choose(bound):
            position = len(used)
            used.append(bound)
            if position >= len(path):
                path.append(0)
            return path[position]

        mask = draw_function(choose)
        del path[len(used):]
        probability = Fraction(1)
        for bound in used:
            probability /= bound
        result[mask] += probability
        position = len(used)
        while position > 0:
            path[position - 1] += 1
            if path[position - 1] < used[position - 1]:
                del path[position:]
                break
            position -= 1
        if position == 0:
            return result


def even_uniform(nsite, choose):
    mask, parity = 0, 0
    for orbital in range(2 * nsite - 1):
        if choose(2):
            mask |= 1 << orbital
            parity ^= 1
    return mask | (parity << (2 * nsite - 1))


def balanced_uniform(nsite, choose):
    """Plan's independent oracle: L of the labels U_i, Dbar_i."""
    remaining = nsite
    selected = []
    for label in range(2 * nsite):
        take = choose(2 * nsite - label) < remaining
        selected.append(take)
        remaining -= int(take)
    up = sum(int(selected[i]) << i for i in range(nsite))
    down = sum(int(not selected[nsite + i]) << (nsite + i)
               for i in range(nsite))
    return up | down


def subset(count, k, choose, offset=0):
    mask, remaining = 0, k
    for item in range(count):
        if choose(count - item) < remaining:
            mask |= 1 << (item + offset)
            remaining -= 1
    return mask


def initial(domain, nsite, n, choose):
    if domain == 2:
        return subset(2 * nsite, n, choose)
    return subset(nsite, n // 2, choose) | subset(nsite, n // 2, choose, nsite)


def local(domain, nsite, x, choose):
    def ordered_pair(count):
        a = choose(count)
        b = choose(count - 1)
        return a, b if b < a else b + 1
    if domain == 2:
        a, b = ordered_pair(2 * nsite)
        return x ^ (1 << a) ^ (1 << b)
    if choose(2) == 0:
        i, j = choose(nsite), choose(nsite)
        if ((x >> i) & 1) == ((x >> (j + nsite)) & 1):
            return x ^ (1 << i) ^ (1 << (j + nsite))
        return x
    spin = choose(2)
    if nsite < 2:
        return x
    a, b = ordered_pair(nsite)
    a, b = a + spin * nsite, b + spin * nsite
    if ((x >> a) & 1) != ((x >> b) & 1):
        return x ^ (1 << a) ^ (1 << b)
    return x


def proposal_mode(args):
    dump = pathlib.Path(args.workdir) / "proposal-dump.txt"
    dump.parent.mkdir(parents=True, exist_ok=True)
    done = subprocess.run([args.unit, "--dump", str(dump)],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, timeout=120)
    if done.returncode != 0:
        print(done.stdout)
        print("GCLanczos proposal: unit failed")
        return 1
    actual = defaultdict(dict)
    for line in dump.read_text().splitlines():
        kind, domain, nsite, extra, x, y, num, den = line.split()
        actual[(kind, int(domain), int(nsite), int(extra), int(x))][int(y)] = \
            Fraction(int(num), int(den))
    failures = 0
    rows = 0
    for (kind, domain, nsite, extra, x), row in sorted(actual.items()):
        if kind == "UNIFORM":
            expected = enumerate_draws(
                lambda c: (even_uniform if domain == 2 else balanced_uniform)(
                    nsite, c))
        elif kind == "INITIAL":
            expected = enumerate_draws(
                lambda c: initial(domain, nsite, extra, c))
        else:
            expected = enumerate_draws(lambda c: local(domain, nsite, x, c))
        expected = dict((k, v) for k, v in expected.items() if v != 0)
        if expected != row:
            failures += 1
            print("FAIL %s domain=%d L=%d extra=%d x=%d" %
                  (kind, domain, nsite, extra, x))
        if sum(row.values()) != 1:
            failures += 1
            print("FAIL row sum %s domain=%d L=%d x=%d" %
                  (kind, domain, nsite, x))
        rows += 1
    # Exact symmetry of the local kernels from the dump itself.
    for (kind, domain, nsite, extra, x), row in actual.items():
        if kind != "LOCAL":
            continue
        for y, value in row.items():
            reverse = actual.get(("LOCAL", domain, nsite, 0, y), {}).get(x)
            if reverse != value:
                failures += 1
                print("FAIL symmetry domain=%d L=%d %d<->%d" %
                      (domain, nsite, x, y))
    if failures or rows == 0:
        print("GCLanczos proposal: %d failure(s), %d rows" % (failures, rows))
        return 1
    print("GCLanczos proposal: PASS (%d kernels, Fraction exact)" % rows)
    return 0


# --- E2E: expert input -> vmc.out -> xxx_pl_out (11 columns) -------------
import re
import shutil

SERIAL_TIMEOUT = 120
MPI_TIMEOUT = 180
STAT_TIMEOUT = 600
STAT_SEEDS = (17, 29, 43, 71)
STAT_SAMPLES = 8192
STAT_WARMUP = 512
PL_COLUMNS = ("energy energy_se variance variance_se energy_imag "
              "variance_imag tau_int effective_samples retained_rank "
              "condition gevp_residual").split()
PL_HEADER = re.compile(
    r"^# stabilized_power_lanczos version=2 order=1 chains=(\d+) "
    r"samples_per_chain=(\d+) coefficient_samples=(\d+) "
    r"final_samples=(\d+) blocks=(\d+)$")
SUMMARY = re.compile(
    r"^GC_LANCZOS_SUMMARY v1 c0=(\S+),(\S+) c1=(\S+),(\S+) "
    r"ecoeff=(\S+),(\S+) chains=(\d+) bootstrap=(\S+) sectors=(\S+)$",
    re.M)
BAR = "=" * 45 + "\n"
STAT_CASES = ("f2", "g2", "f2_combined")


def e2e_case(name):
    rows = model_rows(2, [(0, 1)], "pbc", V=0.125, J=0.0625)
    if name == "f2":
        return fixture("F2", 2, "antiparallel", F2, rows=rows)
    if name == "g2":
        return fixture("G2", 2, "general", general_matrix(F2), rows=rows)
    if name == "f2_combined":
        return fixture("F2_combined", 2, "antiparallel", F2, "pbc",
                       momentum(2, 0), gauss=2, rows=rows)
    if name == "singular":
        return fixture("Singular", 2, "antiparallel", SINGULAR, rows=rows)
    if name == "sample":
        # samples/GrandCanonical/FirstLanczos: F2 without V and Hund.
        return fixture("Sample", 2, "antiparallel", F2,
                       rows=model_rows(2, [(0, 1)], "pbc"))
    if name == "f2_diagonal":
        # t=0: q_k(x) = E(x)^k Psi(x), so every sampled S/K/T is real
        # symmetric and the GEVP input gate holds for any sample count.
        return fixture("F2_diagonal", 2, "antiparallel", F2,
                       rows=model_rows(2, [(0, 1)], "pbc", t=0.0, V=0.125,
                                       J=0.0625))
    raise CheckFailure("unknown E2E case %s" % name)


def chain_count(nproc, nsplit):
    return max(1, (nproc + nsplit - 1) // nsplit)


def block_count(samples):
    count = 16 if samples >= 32 else samples // 2
    while count >= 4 and samples % count:
        count -= 1
    return count if count >= 4 else 0


def def_block(keyword, count, rows, complex_type=None):
    second = ("ComplexType %d\n" % complex_type
              if complex_type is not None else BAR)
    return BAR + "%s %d\n" % (keyword, count) + second + BAR + BAR + \
        "".join(rows)


def e2e_modpara(case, seed, samples, warmup, nsplit, init_nelec):
    general = case["pairing"] == "general"
    return [
        ("NVMCCalMode", 1), ("NLanczosMode", 1), ("NLanczosStep", 1),
        ("NLanczosEstimatorMode", 1), None,
        ("NDataIdxStart", 1), ("NDataQtySmp", 1), None,
        ("Nsite", case["L"]), ("Ne", case["L"] // 2),
        ("2Sz", -1 if general else 0),
        ("NSPGaussLeg", case["gauss"]), ("NSPStot", 0),
        ("NMPTrans", len(case["translations"])),
        ("NSROptItrStep", 1), ("NSROptItrSmp", 1),
        ("DSROptRedCut", "0.000000000001"), ("DSROptStaDel", "0.01"),
        ("DSROptStepDt", "0.002"),
        ("NVMCWarmUp", warmup), ("NVMCInterval", 1),
        ("NVMCSample", samples), ("NExUpdatePath", 0), ("RndSeed", seed),
        ("NSplitSize", nsplit), ("NStore", 1), ("NSRCG", 0),
        ("NGrandCanonical", 1), ("NGCInitNelec", init_nelec)]


def write_modpara(path, entries, overrides):
    keys = set(entry[0] for entry in entries if entry is not None)
    lines = ["--------------------", "Model_Parameters   0",
             "--------------------", "VMC_Cal_Parameters",
             "--------------------", "CDataFileHead  zvo",
             "CParaFileHead  zqp", "--------------------"]
    for entry in entries:
        if entry is None:
            lines.append("--------------------")
            continue
        key, value = entry
        value = overrides.get(key, value)
        if value is not None:
            lines.append("%-22s %s" % (key, value))
    for key, value in overrides.items():
        if key not in keys and value is not None:
            lines.append("%-22s %s" % (key, value))
    pathlib.Path(path).write_text("\n".join(lines) + "\n")


def write_e2e_input(case, workdir, seed=17, samples=STAT_SAMPLES,
                    warmup=STAT_WARMUP, nsplit=1, init_nelec=2,
                    overrides=None):
    """Expert input for one fixture; returns the namelist entries."""
    work = pathlib.Path(workdir)
    if work.exists():
        shutil.rmtree(str(work))
    work.mkdir(parents=True)
    nsite = case["L"]
    write_modpara(work / "modpara.def",
                  e2e_modpara(case, seed, samples, warmup, nsplit,
                              init_nelec), overrides or {})
    names = [("ModPara", "modpara.def"), ("LocSpin", "locspn.def")]
    (work / "locspn.def").write_text(def_block(
        "NlocalSpin", 0, ["%d 0\n" % s for s in range(nsite)]))
    rows = ["%d %d %d %d %s %s\n" % (i, si, j, sj, fmt(complex(v).real),
                                      fmt(complex(v).imag))
            for i, si, j, sj, v in case["trans"]]
    (work / "trans.def").write_text(def_block("NTransfer", len(rows), rows))
    names.append(("Trans", "trans.def"))
    if case.get("intra"):
        (work / "coulombintra.def").write_text(def_block(
            "NCoulombIntra", len(case["intra"]),
            ["%d %s\n" % (s, fmt(v)) for s, v in case["intra"]]))
        names.append(("CoulombIntra", "coulombintra.def"))
    if case.get("inter"):
        (work / "coulombinter.def").write_text(def_block(
            "NCoulombInter", len(case["inter"]),
            ["%d %d %s\n" % (i, j, fmt(v)) for i, j, v in case["inter"]]))
        names.append(("CoulombInter", "coulombinter.def"))
    if case.get("hund"):
        (work / "hund.def").write_text(def_block(
            "NHund", len(case["hund"]),
            ["%d %d %s\n" % (i, j, fmt(v)) for i, j, v in case["hund"]]))
        names.append(("Hund", "hund.def"))
    if case["pairing"] == "antiparallel":
        pairs = [(i, j) for i in range(nsite) for j in range(nsite)]
        rows = ["%d %d %d\n" % (i, j, k) for k, (i, j) in enumerate(pairs)]
        values = [complex(case["F"][i][j]) for i, j in pairs]
        name = ("OrbitalAntiParallel", "orbitalidx.def")
    else:
        pairs = [(a, b) for a in range(2 * nsite)
                 for b in range(a + 1, 2 * nsite)]
        rows = ["%d %d %d %d %d 1\n" % (a % nsite, a // nsite, b % nsite,
                                         b // nsite, k)
                for k, (a, b) in enumerate(pairs)]
        # OrbitalGeneral builds SlaterElm_IJ = f_IJ - f_JI = 2 f_IJ.
        values = [complex(case["A"][a][b]) / 2 for a, b in pairs]
        name = ("OrbitalGeneral", "orbitalidxgen.def")
    rows += ["%d 1\n" % k for k in range(len(pairs))]
    (work / name[1]).write_text(def_block("NOrbitalIdx", len(pairs), rows,
                                          complex_type=1))
    names.append(name)
    steps = [r for r, _ in case["translations"]]
    if steps != list(range(len(steps))):
        raise CheckFailure("E2E translations must be 0..R-1")
    mappings = gl.chain_mappings(nsite, case["boundary"], steps)
    rows = ["%d %s %s\n" % (r, fmt(complex(w).real), fmt(complex(w).imag))
            for r, w in case["translations"]]
    rows += ["%d %d %d %d\n" % ((r, s) + mappings[(r, s)])
             for r in steps for s in range(nsite)]
    (work / "qptransidx.def").write_text(def_block("NQPTrans", len(steps),
                                                   rows))
    names.append(("TransSym", "qptransidx.def"))
    initial = [0.0] * 6
    for value in values:
        initial += [value.real, value.imag, 0.0]
    (work / "initial.def").write_text(
        " ".join(fmt(v) for v in initial) + "\n")
    write_namelist(work, names)
    return names


def write_namelist(work, names):
    (pathlib.Path(work) / "namelist.def").write_text(
        "".join("%16s  %s\n" % (key, path) for key, path in names))


def mpi_prefix(nproc):
    launcher = os.environ.get("MVMC_MPIEXEC") or shutil.which("mpiexec") \
        or shutil.which("mpirun")
    if not launcher:
        raise CheckFailure("MPI launcher not found")
    flag = os.environ.get("MVMC_MPIEXEC_NUMPROC_FLAG") or "-np"
    return [launcher, flag, str(nproc)]


def run_vmc(vmc, workdir, nproc=1, env=None, timeout=None, options=()):
    command = [os.path.abspath(vmc)] + list(options) + [
        "-e", "namelist.def", "initial.def"]
    if nproc > 1:
        command = mpi_prefix(nproc) + command
    if timeout is None:
        timeout = MPI_TIMEOUT if nproc > 1 else SERIAL_TIMEOUT
    environment = dict(os.environ, OMP_NUM_THREADS="1",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       VECLIB_MAXIMUM_THREADS="1")
    environment.pop("MVMC_GC_LANCZOS_TEST_SUMMARY", None)
    environment.pop("MVMC_GC_LANCZOS_TEST_FAULT", None)
    environment.update(env or {})
    try:
        done = subprocess.run(command, cwd=str(workdir), env=environment,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise CheckFailure("%s: timeout after %d s (never an expected error)"
                           % (workdir, timeout))
    work = pathlib.Path(workdir)
    (work / "stdout.txt").write_text(done.stdout)
    (work / "stderr.txt").write_text(done.stderr)
    return done


def pl_out_rows(workdir):
    path = pathlib.Path(workdir) / "output" / "zvo_pl_out_001.dat"
    if not path.exists():
        return None, []
    lines = path.read_text().splitlines()
    headers = [line for line in lines if line.startswith("#")]
    rows = [line.split() for line in lines
            if line.strip() and not line.startswith("#")]
    return headers, rows


def require_success(label, done, workdir, samples, chains):
    """Exit 0 and exactly one finite 11-column row with a valid header."""
    require(done.returncode == 0, "%s: exit %d\n%s%s" % (
        label, done.returncode, done.stdout[-3000:], done.stderr[-3000:]))
    headers, rows = pl_out_rows(workdir)
    require(headers is not None, "%s: no pl_out file" % label)
    require(len(headers) == 2, "%s: header lines %r" % (label, headers))
    match = PL_HEADER.match(headers[0])
    require(match is not None, "%s: header %r" % (label, headers[0]))
    total = samples * chains
    expected = (chains, samples, total, total, block_count(total))
    require(tuple(int(v) for v in match.groups()) == expected,
            "%s: header %r, expected %r" % (label, headers[0], expected))
    require(headers[1] == "# " + " ".join(PL_COLUMNS),
            "%s: column order %r" % (label, headers[1]))
    require(len(rows) == 1 and len(rows[0]) == 11,
            "%s: data rows %r" % (label, rows))
    values = [float(v) for v in rows[0]]
    require(all(math.isfinite(v) for v in values),
            "%s: nonfinite row %r" % (label, rows[0]))
    result = dict(zip(PL_COLUMNS, values))
    require(result["retained_rank"] in (1.0, 2.0),
            "%s: retained_rank %r" % (label, result["retained_rank"]))
    require(result["energy_se"] >= 0 and result["variance_se"] >= 0 and
            result["condition"] >= 1 and result["gevp_residual"] >= 0,
            "%s: invalid diagnostics %r" % (label, result))
    require(0 < result["effective_samples"] <= total and
            result["tau_int"] >= 1,
            "%s: ESS/tau %r" % (label, result))
    return result


def require_failure(label, done, workdir, needles):
    output = done.stdout + done.stderr
    require(done.returncode != 0, "%s: unexpectedly succeeded\n%s" % (
        label, output[-3000:]))
    for needle in needles:
        require(needle in output, "%s: missing diagnostic %r\n%s" % (
            label, needle, output[-3000:]))
    unused, rows = pl_out_rows(workdir)
    require(not rows, "%s: %d result rows written on failure" % (
        label, len(rows)))


def parse_summary(label, stdout, chains, nsite):
    found = SUMMARY.findall(stdout)
    require(len(found) == 1, "%s: %d summary records" % (label, len(found)))
    fields = found[0]
    c = np.array([complex(float(fields[0]), float(fields[1])),
                  complex(float(fields[2]), float(fields[3]))])
    ecoeff = complex(float(fields[4]), float(fields[5]))
    require(int(fields[6]) == chains, "%s: summary chains %s" % (
        label, fields[6]))
    bootstrap = [tuple(int(v) for v in item.split(":"))
                 for item in fields[7].split(";")]
    require(len(bootstrap) == min(chains, 16), "%s: bootstrap %r" % (
        label, fields[7]))
    sectors = [int(v) for v in fields[8].split(",")]
    require(len(sectors) == 2 * nsite + 1, "%s: sectors %r" % (
        label, sectors))
    require(np.all(np.isfinite(c)) and np.any(c != 0),
            "%s: coefficients %r" % (label, c))
    return {"c": c, "ecoeff": ecoeff, "bootstrap": bootstrap,
            "sectors": sectors}


def oracle_reference(case):
    basis, psi = gl.reference_wave(case)
    h = gl.hamiltonian(case, basis)
    norm = float(np.max(np.abs(np.linalg.eigvalsh(h))))
    return basis, psi, h, norm


def oracle_observables(reference, c):
    """E, variance and N probabilities of Phi = c0 Psi + c1 H Psi."""
    basis, psi, h, unused = reference
    phi = c[0] * psi + c[1] * (h @ psi)
    weight = np.abs(phi) ** 2
    total = float(np.sum(weight))
    require(total > 0, "oracle Phi vanishes")
    hphi = h @ phi
    energy = (np.vdot(phi, hphi) / total).real
    second = (np.vdot(hphi, hphi) / total).real
    probability = {}
    for mask, value in zip(basis, weight):
        n = gl.popcount(mask)
        probability[n] = probability.get(n, 0.0) + value / total
    return energy, second - energy * energy, probability


def statistical_gates(label, result, summary, reference, samples, chains):
    """Per-run gates of fixture-policy; returns the oracle probabilities."""
    norm = reference[3]
    energy, variance, probability = oracle_observables(reference,
                                                       summary["c"])
    scale = max(1.0, norm)
    checks = [
        ("energy", abs(result["energy"] - energy),
         8 * result["energy_se"] + 2e-10 * (1 + abs(energy))),
        ("energy_se", result["energy_se"], 0.03 * scale),
        ("variance", abs(result["variance"] - variance),
         8 * result["variance_se"] + 0.02 * scale ** 2),
        ("variance_se", result["variance_se"], 0.05 * scale ** 2),
    ]
    for name, value, bound in checks:
        require(value <= bound, "%s: %s %.6g > %.6g (oracle E=%.12g var=%.12g)"
                % (label, name, value, bound, energy, variance))
    require(result["effective_samples"] >= 128,
            "%s: ESS %.6g < 128" % (label, result["effective_samples"]))
    require(sum(summary["sectors"]) == samples * chains,
            "%s: sector counts %r" % (label, summary["sectors"]))
    for n, count in enumerate(summary["sectors"]):
        if count:
            require(n % 2 == 0 and probability.get(n, 0.0) > 0,
                    "%s: visited forbidden sector N=%d" % (label, n))
        if probability.get(n, 0.0) >= 0.05:
            require(count > 0, "%s: sector N=%d (p=%.3f) not visited" % (
                label, n, probability[n]))
    for entry in summary["bootstrap"]:
        require(entry[0] >= 1 and entry[4] == 0,
                "%s: bootstrap %r" % (label, entry))
    print("  %-28s E=%.10f oracle=%.10f |d|/se=%.2f se=%.4g var=%.6f "
          "oracle=%.6f tau=%.2f ESS=%.0f E_coeff-E=%.3g" % (
              label, result["energy"], energy,
              abs(result["energy"] - energy) / max(result["energy_se"],
                                                   1e-300),
              result["energy_se"], result["variance"], variance,
              result["tau_int"], result["effective_samples"],
              summary["ecoeff"].real - result["energy"]))
    return probability


def stat_run(args, case_name, seed, nproc=1, nsplit=1, tag=""):
    case = e2e_case(case_name)
    label = "%s_seed%d_np%d_split%d%s" % (case_name, seed, nproc, nsplit, tag)
    work = pathlib.Path(args.workdir) / label
    write_e2e_input(case, work, seed=seed, nsplit=nsplit)
    done = run_vmc(args.vmc, work, nproc,
                   env={"MVMC_GC_LANCZOS_TEST_SUMMARY": "1"},
                   timeout=STAT_TIMEOUT)
    chains = chain_count(nproc, nsplit)
    result = require_success(label, done, work, STAT_SAMPLES, chains)
    summary = parse_summary(label, done.stdout, chains, case["L"])
    return label, case, result, summary, chains


def production_mode(args):
    failures = 0
    for case_name in STAT_CASES:
        if args.case and args.case != case_name:
            continue
        reference = oracle_reference(e2e_case(case_name))
        pooled = None
        oracle_mean = {}
        try:
            for seed in STAT_SEEDS:
                label, case, result, summary, chains = stat_run(
                    args, case_name, seed)
                probability = statistical_gates(
                    label, result, summary, reference, STAT_SAMPLES, chains)
                counts = np.array(summary["sectors"], dtype=float)
                pooled = counts if pooled is None else pooled + counts
                for n, p in probability.items():
                    oracle_mean[n] = oracle_mean.get(n, 0.0) + \
                        p / len(STAT_SEEDS)
            frequency = pooled / np.sum(pooled)
            for n in range(len(frequency)):
                diff = abs(frequency[n] - oracle_mean.get(n, 0.0))
                require(diff <= 0.06, "%s: pooled N=%d frequency %.4f vs "
                        "oracle %.4f" % (case_name, n, frequency[n],
                                         oracle_mean.get(n, 0.0)))
            print("PASS %-12s ||H||=%.5f pooled N freq=%s oracle=%s" % (
                case_name, reference[3],
                ",".join("%.4f" % f for f in frequency),
                ",".join("%.4f" % oracle_mean.get(n, 0.0)
                         for n in range(len(frequency)))))
        except CheckFailure as error:
            failures += 1
            print("FAIL %s" % error)
    print("GCLanczos production: %s" % ("PASS" if not failures else
                                        "%d failure(s)" % failures))
    return 1 if failures else 0


def deterministic_match(label, first, second):
    """Same chain seeds, only the QP partition changed."""
    (r1, s1), (r2, s2) = first, second
    for key in PL_COLUMNS:
        a, b = r1[key], r2[key]
        require(abs(a - b) <= 2e-10 * (1 + abs(b)),
                "%s: %s %.17g vs %.17g" % (label, key, a, b))
    for index in range(2):
        a, b = s1["c"][index], s2["c"][index]
        require(abs(a - b) <= 2e-10 * (1 + abs(b)),
                "%s: c%d %r vs %r" % (label, index, a, b))
    require(s1["sectors"] == s2["sectors"] and
            s1["bootstrap"] == s2["bootstrap"],
            "%s: counts differ %r/%r vs %r/%r" % (
                label, s1["sectors"], s1["bootstrap"], s2["sectors"],
                s2["bootstrap"]))


def mpi_mode(args):
    require(args.np in (2, 3), "mpi mode requires --np 2 or 3")
    failures = 0
    seed = STAT_SEEDS[0]
    try:
        if args.np == 2:
            # QP partition of one chain against the serial run.
            serial = stat_run(args, "f2_combined", seed)
            split = stat_run(args, "f2_combined", seed, nproc=2, nsplit=2)
            require(split[4] == 1, "np2 split: chains %d" % split[4])
            deterministic_match("f2_combined QP split np2",
                                (split[2], split[3]), (serial[2], serial[3]))
            print("PASS f2_combined QP partition (np2, 1 chain) == serial")
            # Two independent chains share one coefficient.
            for case_name in ("f2", "g2"):
                label, case, result, summary, chains = stat_run(
                    args, case_name, seed, nproc=2, nsplit=1)
                require(chains == 2, "%s: chains %d" % (label, chains))
                statistical_gates(label, result, summary,
                                  oracle_reference(case), STAT_SAMPLES,
                                  chains)
                print("PASS %s two chains (np2)" % case_name)
        else:
            # Chains of size 2 + 1; F2 has one QP, so one rank is empty.
            reference = stat_run(args, "f2", seed, nproc=2, nsplit=1,
                                 tag="_ref")
            uneven = stat_run(args, "f2", seed, nproc=3, nsplit=2)
            require(uneven[4] == 2, "np3: chains %d" % uneven[4])
            statistical_gates(uneven[0], uneven[2], uneven[3],
                              oracle_reference(uneven[1]), STAT_SAMPLES, 2)
            deterministic_match("f2 np3 (2+1, empty QP) vs np2 (1+1)",
                                (uneven[2], uneven[3]),
                                (reference[2], reference[3]))
            print("PASS f2 uneven chains 2+1 with an empty QP rank (np3)")
            combined_ref = stat_run(args, "f2_combined", seed, nproc=2,
                                    nsplit=1, tag="_ref")
            combined = stat_run(args, "f2_combined", seed, nproc=3,
                                nsplit=2)
            deterministic_match("f2_combined np3 (2+1) vs np2 (1+1)",
                                (combined[2], combined[3]),
                                (combined_ref[2], combined_ref[3]))
            print("PASS f2_combined uneven chains 2+1 (np3)")
    except CheckFailure as error:
        failures += 1
        print("FAIL %s" % error)
    print("GCLanczos mpi np%d: %s" % (args.np, "PASS" if not failures else
                                      "%d failure(s)" % failures))
    return 1 if failures else 0


# --- input mutations -------------------------------------------------------
def add_file(work, keyword, name, text):
    work = pathlib.Path(work)
    (work / name).write_text(text)
    with open(str(work / "namelist.def"), "a") as stream:
        stream.write("%16s  %s\n" % (keyword, name))


def edit_rows(work, name, function):
    path = pathlib.Path(work) / name
    lines = path.read_text().splitlines(True)
    path.write_text("".join(function(lines)))


def set_count(lines, count):
    keyword = lines[1].split()[0]
    return lines[:1] + ["%s %d\n" % (keyword, count)] + lines[2:]


def append_trans(work, rows):
    edit_rows(work, "trans.def", lambda lines: set_count(
        lines, len(lines) - 5 + len(rows)) + rows)


def nonhermitian(work):
    def change(lines):
        body = lines[5:]
        fields = body[1].split()
        fields[4] = fmt(0.5 * float(fields[4]))
        body[1] = " ".join(fields) + "\n"
        return lines[:5] + body
    edit_rows(work, "trans.def", change)


SPIN_JASTROW = def_block("NSpinJastrowIdx", 1,
                         ["0 1 0\n", "1 0 0\n", "0 1\n"], complex_type=0)
LATTICE = BAR.replace("=", "-") + "NLattice 2 1 1 1\n" + \
    BAR.replace("=", "-") + "i_x_y_z_orb\n" + BAR.replace("=", "-") + \
    "0 0 0 0 0\n1 1 0 0 0\n"
TWIST = BAR.replace("=", "-") + "NTwist 1\n" + BAR.replace("=", "-") + \
    "idx_site_s_x_y_z\n" + BAR.replace("=", "-") + \
    "".join("0 %d %d 0.0 0.0 0.0\n" % (site, spin)
            for site in range(2) for spin in range(2))
OPTTRANS = BAR + "NQPOptTrans 1\n" + BAR + BAR + BAR + "0 1.0\n" + \
    "".join("0 %d %d\n" % (site, site) for site in range(2))
FIRST_STEP = "supports Lanczos only as the first-step stabilized estimator"
TERMS = "Lanczos supports only Trans, CoulombIntra, CoulombInter and Hund"
GREEN = "Lanczos does not measure Green functions"

# name: (base case, mutation(work), diagnostics)
INVALID_INPUTS = {
    "legacy": ("f2", {"NLanczosEstimatorMode": 0}, [FIRST_STEP]),
    "step2": ("f2", {"NLanczosStep": 2}, [FIRST_STEP]),
    "mode2": ("f2", {"NLanczosMode": 2}, [FIRST_STEP]),
    "calmode0": ("f2", {"NVMCCalMode": 0}, [FIRST_STEP]),
    "reweight": ("f2", {"reweight": 1}, ["Lanczos requires reweight=0"]),
    "spinflip_transfer": ("g2", lambda w: append_trans(w, [
        "0 0 1 1 %s 0.0\n" % fmt(0.1), "1 1 0 0 %s 0.0\n" % fmt(0.1)]),
        ["Lanczos supports only spin-conserving Trans"]),
    "pair_hopping": ("f2", lambda w: add_file(
        w, "PairHop", "pairhop.def", def_block("NPairHop", 1,
                                               ["0 1 0.1\n"])), [TERMS]),
    "exchange": ("f2", lambda w: add_file(
        w, "Exchange", "exchange.def", def_block("NExchange", 1,
                                                 ["0 1 0.1\n"])), [TERMS]),
    "interall": ("f2", lambda w: add_file(
        w, "InterAll", "interall.def", def_block(
            "NInterAll", 1, ["0 0 0 0 1 0 1 0 0.1 0.0\n"])), [TERMS]),
    "nbody_interall": ("f2", lambda w: add_file(
        w, "NBodyInterAll", "nbodyinterall.def", def_block(
            "NNBodyInterAll", 1, ["1 0 0 0 0 0.1 0\n"])), [TERMS]),
    "anomalous_term": ("f2", lambda w: add_file(
        w, "AnomalousTerm", "anomalousterm.def", def_block(
            "NAnomalousTerm", 2, ["1 0 0 0 1 0.35 -0.20\n",
                                  "0 0 1 0 0 0.35 0.20\n"])), [TERMS]),
    "onebodyg": ("f2", lambda w: add_file(
        w, "OneBodyG", "greenone.def", def_block(
            "NCisAjs", 1, ["0 0 0 0\n"])), [GREEN]),
    "twobodyg": ("f2", lambda w: add_file(
        w, "TwoBodyG", "greentwo.def", def_block(
            "NCisAjsCktAltDC", 1, ["0 0 0 0 1 0 1 0\n"])), [GREEN]),
    "twobodygex": ("f2", lambda w: add_file(
        w, "TwoBodyGEx", "greentwoex.def", def_block(
            "NCisAjsCktAlt", 1, ["0 0 1 0 1 1 0 1\n"])), [GREEN]),
    "nbodyg": ("f2", lambda w: add_file(
        w, "NBodyG", "nbodyg.def", def_block(
            "NNBodyG", 1, ["1 0 0 0 0\n"])), [GREEN]),
    "anomalousg": ("f2", lambda w: add_file(
        w, "AnomalousG", "anomalousg.def", def_block(
            "NAnomalousG", 1, ["1 0 0 0 1\n"])), [GREEN]),
    "independent_h": ("f2", lambda w: add_file(
        w, "LsTrans", "lstrans.def", def_block(
            "NLsTransfer", 1, ["0 0 0 0 -0.30 0.0\n"])),
        ["Lanczos does not support an independent operator"]),
    "general_spin_projection": ("g2", {"NSPGaussLeg": 2},
                                ["GC spin projection requires anti-parallel"]),
    "spin_jastrow_projection": ("f2", ({"NSPGaussLeg": 2}, lambda w: add_file(
        w, "SpinJastrow", "spinjastrow.def", SPIN_JASTROW)),
        ["spin projection does not support SpinJastrow"]),
    "opttrans": ("f2", lambda w: add_file(w, "OptTrans", "opttrans.def",
                                          OPTTRANS),
                 ["does not support OptTrans"]),
    "twist": ("f2", lambda w: (
        add_file(w, "Lattice", "lattice.def", LATTICE),
        add_file(w, "Twist", "twist.def", TWIST)),
        ["Lanczos does not support Twist"]),
    "backflow": ("f2", lambda w: (
        add_file(w, "BF", "bf.def", "====================\n"
                 "NBackFlowIdx 1\n====================\n"),
        add_file(w, "BFRange", "rangebf.def", "====================\n"
                 "Nrange 1 1\n====================\n")),
        ["does not support BackFlow"]),
    "rbm": ("f2", lambda w: add_file(
        w, "GeneralRBM_HiddenLayer", "rbm.def",
        BAR + "NRBM_HiddenLayerIdx 1\nComplexType          1\n" + BAR),
        ["does not support RBM"]),
    "locspin": ("f2", lambda w: edit_rows(w, "locspn.def", lambda lines: (
        lines[:1] + ["NlocalSpin 1\n"] + lines[2:5] + ["0 1\n"] +
        lines[6:])), ["does not support LocSpin"]),
    "real_only": ("f2", lambda w: edit_rows(
        w, "orbitalidx.def", lambda lines: lines[:2] +
        ["ComplexType 0\n"] + lines[3:]),
        ["requires complex variational parameters"]),
    "update_path": ("f2", {"NExUpdatePath": 1}, ["requires NExUpdatePath=0"]),
    "negative_site": ("f2", {"Nsite": -2}, ["Nsite must be positive"]),
    "negative_qp": ("f2", {"NSPGaussLeg": -1},
                    ["GC NSPGaussLeg must be positive"]),
    "negative_sample": ("f2", {"NVMCSample": -16},
                        ["NVMCSample must be non-negative"]),
    "sample9_chain2": ("f2", {"NVMCSample": 9},
                       ["requires NVMCSample >= 8 and a divisor"]),
    "initial_n_negative": ("f2", {"NGCInitNelec": -2},
                           ["0 <= NGCInitNelec <= 2*Nsite"]),
    "initial_n_too_large": ("f2", {"NGCInitNelec": 5},
                            ["0 <= NGCInitNelec <= 2*Nsite"]),
    "size_overflow": ("f2", {"NSPGaussLeg": 2147483647},
                      ["GC NSPGaussLeg workspace exceeds int range"]),
    "nonhermitian": ("f2", nonhermitian,
                     ["Lanczos requires a Hermitian Trans input"]),
}


# OptTrans is enabled by the -o option, as in the existing GC tests.
OPTIONS = {"opttrans": ("-o",)}


def split_mutation(mutation):
    """A mutation is modpara overrides, a file edit, or (overrides, edit)."""
    if isinstance(mutation, tuple):
        return mutation
    if isinstance(mutation, dict):
        return mutation, None
    return {}, mutation


def input_case(args, name, base, mutation, samples=64, warmup=8):
    case = e2e_case(base)
    work = pathlib.Path(args.workdir) / ("%s_np%d" % (name, args.np))
    overrides, edit = split_mutation(mutation)
    write_e2e_input(case, work, samples=samples, warmup=warmup,
                    overrides=overrides)
    if edit is not None:
        edit(work)
    return case, work


def positive_inputs(args):
    """Accepted inputs: (name, base, overrides/mutation, samples, check).

    A sampled K is Hermitian only in expectation, and the existing GEVP
    input gate rejects a large anti-Hermitian part of a few-sample estimate
    (canonical runs behave the same).  Small sample counts therefore use the
    diagonal-H fixture; the other accepted inputs use 8192 samples.
    """
    chains = args.np
    items = []
    for samples in (8, 16, 32):
        items.append(("samples%d" % samples, "f2_diagonal", None, samples,
                      None))
    items.append(("init_odd_rounded", "f2", {"NGCInitNelec": 3},
                  STAT_SAMPLES, "initial GC particle number is 2"))
    items.append(("zero_element_green", "f2", lambda w: (
        add_file(w, "OneBodyG", "greenone.def",
                 def_block("NCisAjs", 0, [])),
        add_file(w, "TwoBodyG", "greentwo.def",
                 def_block("NCisAjsCktAltDC", 0, []))), STAT_SAMPLES, None))
    items.append(("general_g2", "g2", None, STAT_SAMPLES, None))
    return chains, items


def input_mode(args):
    names = tuple(INVALID_INPUTS)
    require(len(set(names)) == len(names) == 34, "mutation inventory")
    failures = 0
    for name, (base, mutation, needles) in INVALID_INPUTS.items():
        if args.case and args.case != name:
            continue
        try:
            unused, work = input_case(args, name, base, mutation)
            done = run_vmc(args.vmc, work, args.np,
                           options=OPTIONS.get(name, ()))
            require_failure(name, done, work, needles)
            print("PASS invalid %-24s %s" % (name, needles[0]))
        except CheckFailure as error:
            failures += 1
            print("FAIL %s" % error)
    chains, items = positive_inputs(args)
    for name, base, mutation, samples, needle in items:
        if args.case and args.case != name:
            continue
        try:
            unused, work = input_case(args, name, base, mutation or {},
                                      samples=samples)
            done = run_vmc(args.vmc, work, args.np)
            require_success(name, done, work, samples, chains)
            if needle:
                require(needle in done.stdout, "%s: missing %r" % (
                    name, needle))
            print("PASS accepted %-23s chains=%d samples=%d" % (
                name, chains, samples))
        except CheckFailure as error:
            failures += 1
            print("FAIL %s" % error)
    for name, function in (("normal_gc", check_normal_gc),
                           ("canonical_green", check_canonical_green)):
        if args.case and args.case != name:
            continue
        try:
            function(args, name)
        except CheckFailure as error:
            failures += 1
            print("FAIL %s" % error)
    print("GCLanczos input np%d: %s" % (args.np, "PASS" if not failures else
                                        "%d failure(s)" % failures))
    return 1 if failures else 0


def check_normal_gc(args, name):
    """NLanczosMode=0 keeps the ordinary GC sampler (no pl_out)."""
    unused, work = input_case(args, name, "f2", {"NLanczosMode": 0,
                                                 "NLanczosStep": None,
                                                 "NLanczosEstimatorMode":
                                                 None})
    done = run_vmc(args.vmc, work, args.np)
    require(done.returncode == 0, "%s: exit %d\n%s" % (
        name, done.returncode, done.stderr[-3000:]))
    headers, rows = pl_out_rows(work)
    require(headers is None, "%s: pl_out written by normal GC" % name)
    out = work / "output" / "zvo_out_001.dat"
    require(out.exists() and out.read_text().split(),
            "%s: no physical output" % name)
    gc = work / "zvo_gc.dat"
    require(gc.exists() and gc.read_text().split(),
            "%s: no zvo_gc.dat record" % name)
    print("PASS accepted %-23s ordinary GC sampler" % name)


def check_canonical_green(args, name):
    """Canonical stabilized still accepts (and skips) OneBodyG/TwoBodyG."""
    unused, work = input_case(args, name, "f2", lambda w: (
        add_file(w, "OneBodyG", "greenone.def",
                 def_block("NCisAjs", 1, ["0 0 0 0\n"])),
        add_file(w, "TwoBodyG", "greentwo.def",
                 def_block("NCisAjsCktAltDC", 1, ["0 0 0 0 1 0 1 0\n"]))))
    write_modpara(work / "modpara.def",
                  e2e_modpara(e2e_case("f2"), 17, STAT_SAMPLES, 8, 1, 2),
                  {"NGrandCanonical": None, "NGCInitNelec": None,
                   "Ncond": 2})
    done = run_vmc(args.vmc, work, args.np)
    require_success(name, done, work, STAT_SAMPLES, args.np)
    print("PASS accepted %-23s canonical Green input unchanged" % name)


# --- failure contracts ------------------------------------------------------
def failure_mode(args):
    failures = 0
    items = [
        ("bootstrap_numeric", "f2", "bootstrap_numeric", None),
        ("bootstrap_exhausted", "f2", "bootstrap_exhausted",
         ["GC Lanczos failure: stage=final-bootstrap",
          "attempts=1024 exact_zero=0 numeric_zero=1024 local_energy=0 "
          "exhausted=1", "amplitude_state=numeric_zero"]),
        ("final_numeric", "f2", "final_numeric",
         ["GC Lanczos failure: stage=final-sample",
          "amplitude_state=numeric_zero"]),
        ("unknown_fault", "f2", "no_such_fault",
         ["unknown MVMC_GC_LANCZOS_TEST_FAULT value"]),
    ]
    if args.np > 1:
        items += [
            ("allocation_rank1", "f2_combined", "allocation:rank=1",
             ["stabilized power-Lanczos failed: allocation failure"]),
            # The first evaluation (coefficient start) fails on rank 1;
            # candidate diagnostics cover the final-chain stages only.
            ("nonfinite_rank1", "f2_combined", "nonfinite:rank=1",
             ["stabilized power-Lanczos failed: nonfinite value"]),
        ]
    nsplit = 1 if args.np == 1 else 2
    chains = chain_count(args.np, nsplit)
    for name, base, fault, needles in items:
        if args.case and args.case != name:
            continue
        try:
            case = e2e_case(base)
            work = pathlib.Path(args.workdir) / ("%s_np%d" % (name, args.np))
            write_e2e_input(case, work, samples=STAT_SAMPLES, warmup=8,
                            nsplit=nsplit)
            done = run_vmc(args.vmc, work, args.np, env={
                "MVMC_GC_LANCZOS_TEST_FAULT": fault,
                "MVMC_GC_LANCZOS_TEST_SUMMARY": "1"})
            if needles is None:
                require_success(name, done, work, STAT_SAMPLES, chains)
                summary = parse_summary(name, done.stdout, chains, case["L"])
                retried = summary["bootstrap"][0]
                require(retried[2] >= 1 and retried[4] == 0,
                        "%s: bootstrap %r" % (name, retried))
                print("PASS %-20s numeric-zero candidate redrawn %r" % (
                    name, retried))
            else:
                require_failure(name, done, work, needles)
                require("GC_LANCZOS_SUMMARY" not in done.stdout,
                        "%s: summary printed on failure" % name)
                print("PASS %-20s %s" % (name, needles[0]))
        except CheckFailure as error:
            failures += 1
            print("FAIL %s" % error)
    if not args.case or args.case == "singular":
        try:
            classify_singular(args)
        except CheckFailure as error:
            failures += 1
            print("FAIL %s" % error)
    print("GCLanczos failure np%d: %s" % (args.np, "PASS" if not failures else
                                          "%d failure(s)" % failures))
    return 1 if failures else 0


def classify_singular(args):
    """Real classification of the singular F fixture (no forced fault)."""
    case = e2e_case("singular")
    work = pathlib.Path(args.workdir) / ("singular_np%d" % args.np)
    nsplit = 1 if args.np == 1 else 2
    write_e2e_input(case, work, samples=STAT_SAMPLES, warmup=8,
                    nsplit=nsplit)
    done = run_vmc(args.vmc, work, args.np,
                   env={"MVMC_GC_LANCZOS_TEST_SUMMARY": "1"})
    chains = chain_count(args.np, nsplit)
    bootstrap = re.search(r"^GC Lanczos bootstrap: .*$", done.stdout, re.M)
    if done.returncode == 0:
        result = require_success("singular", done, work, STAT_SAMPLES,
                                 chains)
        summary = parse_summary("singular", done.stdout, chains, case["L"])
        print("PASS singular F: completed, E=%.10f se=%.3g sectors=%s; %s" % (
            result["energy"], result["energy_se"], summary["sectors"],
            bootstrap.group(0) if bootstrap else "no bootstrap line"))
    else:
        require_failure("singular", done, work,
                        ["GC Lanczos failure: stage="])
        failure = re.search(r"^GC Lanczos failure: .*$",
                            done.stdout + done.stderr, re.M)
        print("PASS singular F: terminated, %s" % failure.group(0))


# --- public sample ----------------------------------------------------------
def sample_mode(args):
    source = pathlib.Path(args.sample)
    work = pathlib.Path(args.workdir)
    if work.exists():
        shutil.rmtree(str(work))
    shutil.copytree(str(source), str(work))
    try:
        done = run_vmc(args.vmc, work, 1, timeout=STAT_TIMEOUT,
                       env={"MVMC_GC_LANCZOS_TEST_SUMMARY": "1"})
        result = require_success("sample", done, work, STAT_SAMPLES, 1)
        require("Stabilized power-Lanczos (grand canonical): E=" in
                done.stdout, "sample: missing stdout energy line")
        require(not (work / "zvo_gc.dat").exists(),
                "sample: GC Lanczos must not create zvo_gc.dat")
        # The README describes the same H and pair function as the oracle.
        case = e2e_case("sample")
        summary = parse_summary("sample", done.stdout, 1, case["L"])
        reference = oracle_reference(case)
        statistical_gates("sample", result, summary, reference,
                          STAT_SAMPLES, 1)
        basis, psi, h, unused = reference
        gevp = gl.gevp(*gl.krylov(psi, h)[1:])
        print("PASS sample E=%.10f +/- %.3g variance=%.6f (oracle GEVP "
              "E=%.10f, exact ground state %.10f)" % (
                  result["energy"], result["energy_se"], result["variance"],
                  gevp["energy"], float(np.min(np.linalg.eigvalsh(h)))))
    except CheckFailure as error:
        print("FAIL %s" % error)
        return 1
    print("GCLanczos sample: PASS")
    return 0

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("oracle", "proposal", "input",
                                         "production", "mpi", "failure",
                                         "sample"))
    parser.add_argument("--unit")
    parser.add_argument("--probe")
    parser.add_argument("--vmc")
    parser.add_argument("--case")
    parser.add_argument("--np", type=int, default=1)
    parser.add_argument("--workdir", default="work/GCLanczos")
    parser.add_argument("--sample")
    args = parser.parse_args()
    if args.mode == "oracle":
        if not args.probe:
            parser.error("oracle mode requires --probe")
        return oracle_mode(args)
    if args.mode == "proposal":
        if not args.unit:
            parser.error("proposal mode requires --unit")
        return proposal_mode(args)
    if not args.vmc:
        parser.error("%s mode requires --vmc" % args.mode)
    if args.mode == "sample":
        if not args.sample:
            parser.error("sample mode requires --sample")
        return sample_mode(args)
    modes = {"input": input_mode, "production": production_mode,
             "mpi": mpi_mode, "failure": failure_mode}
    try:
        return modes[args.mode](args)
    except CheckFailure as error:
        print("FAIL %s" % error)
        return 1


if __name__ == "__main__":
    sys.exit(main())
