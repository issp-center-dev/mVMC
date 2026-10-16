#!/usr/bin/env python3
"""GC first-step stabilized Lanczos checks.

modes:
  proposal: exact Fraction cross-check of the symmetric GC proposals
            enumerated by gc_krylov_proposal_unit --dump.
  oracle  : enumeration probe (production GC bridge, bounded engine,
            measurement, GEVP) against the independent Fock oracle
            (grandcanonical_lanczos_oracle) for F2/G2/F4, the six
            pairing/projection combinations x PBC/APBC, boundary, node,
            singular, rank and 2^+-300 scale fixtures.

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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("oracle", "proposal"))
    parser.add_argument("--unit")
    parser.add_argument("--probe")
    parser.add_argument("--vmc")
    parser.add_argument("--case")
    parser.add_argument("--np", type=int, default=1)
    parser.add_argument("--workdir", default="work/GCLanczos")
    args = parser.parse_args()
    if args.mode == "oracle":
        if not args.probe:
            parser.error("oracle mode requires --probe")
        return oracle_mode(args)
    if args.mode == "proposal":
        if not args.unit:
            parser.error("proposal mode requires --unit")
        return proposal_mode(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())
