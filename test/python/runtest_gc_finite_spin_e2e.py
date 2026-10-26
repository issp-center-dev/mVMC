"""Finite-S public-input tests, with independent Fock-space output replay.

All generated inputs/results stay below the CTest working directory.  The
reference rotates creation operators in Fock space, checks dumped quadrature
coefficients independently, and never calls the production Pfaffian evaluator.
"""
import argparse
import cmath
import functools
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import numpy as np

import grandcanonical_spin_oracle as oracle
from grandcanonical_exact_oracle import apply_ops, one_body
from runtest_gc import mpi_command, state_dump_records
from runtest_gc_spin import local_of
from runtest_gc_spin_sector import table, write_input
from runtest_gc_mp import parse_physical, replay_iterations
from runtest_gc_finite_spin import student_t_quantile


def close(label, actual, expected, tolerance=2e-8):
    if not np.isfinite(actual) or not np.isfinite(expected):
        raise AssertionError(label + ": nonfinite comparison")
    if abs(actual - expected) > tolerance * (1 + abs(expected)):
        raise AssertionError("{}: {} != {}".format(label, actual, expected))


def option(work, key, value):
    path = work / "modpara.def"
    text, count = re.subn(r"(?m)^" + re.escape(key) + r"\s+[^\n]+$",
                          "{} {}".format(key, value), path.read_text())
    if count > 1:
        raise AssertionError("duplicate input option " + key)
    path.write_text(text if count else text + "{} {}\n".format(key, value))


class Fixture:
    def __init__(self, name, gauss=8, boundary="pbc", momentum=None, shared=False):
        self.F, self.spin = oracle.finite_spin_fixture(name)
        self.F = np.asarray(self.F)
        self.name, self.L, self.gauss = name, len(self.F), gauss
        self.apbc, self.momentum = boundary == "apbc", momentum
        self.signs = np.ones((self.L, self.L), dtype=int)
        if self.apbc:
            self.signs[np.tril_indices(self.L, -1)] = -1
        self.indices = np.arange(self.L ** 2).reshape(self.L, self.L)
        if shared:
            # Two entries share one variational parameter with opposite signs.
            # This reference does not assume affine dependence on that parameter.
            self.F[1, 0] = -self.F[0, 1]
            self.signs[1, 0] = -self.signs[0, 1]
            self.indices[1, 0] = self.indices[0, 1]
            ids = sorted(set(self.indices.flat))
            self.indices = np.asarray([[ids.index(int(k)) for k in row]
                                       for row in self.indices])
        self.parameters = np.empty(int(self.indices.max()) + 1, dtype=complex)
        for i in range(self.L):
            for j in range(self.L):
                self.parameters[self.indices[i, j]] = self.F[i, j] / self.signs[i, j]
        self.gutz = -0.19
        self.translations = list(range(self.L)) if momentum is not None else [0]
        self.weights = [cmath.exp(-2j * math.pi * (momentum or 0) * r / self.L)
                        / len(self.translations) for r in self.translations]

    def matrix(self, parameters):
        return np.asarray(parameters)[self.indices] * self.signs

    def doublons(self, mask):
        return sum((mask >> i & 1) * (mask >> (i + self.L) & 1)
                   for i in range(self.L))

    def translated(self, F, r):
        return [[F[(i + r) % self.L, (j + r) % self.L] *
                 (-1 if self.apbc and i + r >= self.L else 1) *
                 (-1 if self.apbc and j + r >= self.L else 1)
                 for j in range(self.L)] for i in range(self.L)]

    def write(self, work, mode=1, store=1, reweight=0, samples=128, seed=4321,
              initial=-1, split=1, iterations=2, freeze=True, warmup=128):
        hopping, _ = write_input(work, self.F, self.gauss, initial, seed,
                                  samples, split, self.momentum is not None,
                                  self.apbc)
        for key, value in dict(NSPStot=self.spin, NVMCCalMode=mode, NStore=store,
                               reweight=reweight, NVMCWarmUp=warmup,
                               NMPTrans=(-len(self.translations)
                                         if np.any(self.signs < 0) else len(self.translations)),
                               NDataQtySmp=iterations if mode else 1,
                               NSROptItrStep=iterations if not mode else 1,
                               DSROptStepDt=0 if freeze else 0.002).items():
            option(work, key, value)
        rows = ["{} {} {} {}".format(i, j, self.indices[i, j], self.signs[i, j])
                for i in range(self.L) for j in range(self.L)]
        rows += ["{} 1".format(k) for k in range(len(self.parameters))]
        (work / "orbitalidx.def").write_text(
            table("NOrbitalIdx", len(self.parameters), rows, True))
        values = [0.0] * 6 + [self.gutz, 0.0, 0.0]
        for value in self.parameters:
            values += [value.real, value.imag, 0.0]
        (work / "initial.def").write_text(" ".join(format(x, ".17g") for x in values) + "\n")
        (work / "gutz.def").write_text(table("NGutzwillerIdx", 1,
            ["{} 0".format(i) for i in range(self.L)] + ["0 1"]))
        rows = ["{} {:.17g} {:.17g}".format(r, w.real, w.imag)
                for r, w in zip(self.translations, self.weights)]
        rows += ["{} {} {} {}".format(r, i, (i + r) % self.L,
                  -1 if self.apbc and i + r >= self.L else 1)
                 for r in self.translations for i in range(self.L)]
        (work / "qptransidx.def").write_text(table("NQPTrans", len(self.translations), rows))
        # A number-changing, Sz-preserving Hamiltonian probes relative N weights.
        self.pair_ops = ((("c", 0), ("c", self.L + 1)),
                         (("a", self.L + 1), ("a", 0)))
        self.delta = 0.13 + 0.07j
        (work / "pair.def").write_text(table("NAnomalousTerm", 2, [
            "1 0 0 1 1 {:.17g} {:.17g}".format(self.delta.real, self.delta.imag),
            "0 1 1 0 0 {:.17g} {:.17g}".format(self.delta.real, -self.delta.imag)]))
        # Every term in (S+ S- + S- S+)/2, plus a forbidden Delta Sz term.
        self.two_keys = [
            ((i, i + self.L), (j + self.L, j)) for i in range(self.L)
            for j in range(self.L)] + [
            ((i + self.L, i), (j, j + self.L)) for i in range(self.L)
            for j in range(self.L)]
        self.two_keys.append(((0, self.L), (1, self.L + 1)))
        fields = lambda pairs: " ".join(str(x) for pair in pairs for orb in pair
                                        for x in (orb % self.L, orb // self.L))
        (work / "greentwo.def").write_text(table("NCisAjsCktAltDC",
            len(self.two_keys), [fields(key) for key in self.two_keys]))
        self.nbody_key = ((0, self.L), (self.L + 1, 1), (0, 0))
        (work / "nbodyg.def").write_text(table("NNBodyG", 1,
            ["3 " + fields(self.nbody_key)]))
        # Match the pair observables to the number-changing Hamiltonian.
        (work / "anomalousg.def").write_text(table("NAnomalousG", 2,
            ["1 0 0 1 1", "0 1 1 0 0"]))
        path = work / "namelist.def"
        names = path.read_text() + "Gutzwiller gutz.def\nAnomalousTerm pair.def\n"
        if mode:
            names += "TwoBodyG greentwo.def\nNBodyG nbodyg.def\n"
        else:
            names = "\n".join(line for line in names.splitlines()
                              if not line.startswith(("OneBodyG ", "AnomalousG "))) + "\n"
        path.write_text(names)
        self.hopping = hopping


@functools.lru_cache(maxsize=32)
def rotation(L, gauss, spin):
    return oracle.rotation_operator(L, oracle.independent_spin_nodes(gauss, spin))


class Reference:
    def __init__(self, fixture):
        self.f = fixture
        self.basis, self.rotation = rotation(fixture.L, fixture.gauss, fixture.spin)

    def raw(self, parameters=None, exact=False):
        f = self.f
        F = f.matrix(f.parameters if parameters is None else parameters)
        symmetry = oracle.pair_symmetry(F)
        total = np.zeros(len(self.basis), dtype=complex)
        for r, w in zip(f.translations, f.weights):
            translated = f.translated(F, r)
            if exact:
                values = oracle.raw_exact_spin(translated, f.spin)
                total += w * np.array([oracle.qfloat(values[m]) for m in self.basis])
            else:
                pair = oracle.pair_wave(translated)
                total += w * (self.rotation @ np.array([pair[m] for m in self.basis]))
        # These proven masks are independent of the rotation implementation.
        for i, mask in enumerate(self.basis):
            if oracle.finite_spin_zero(mask, f.L, f.spin, symmetry):
                total[i] = 0
        return dict(zip(self.basis, total))

    def wave(self, parameters=None, gutz=None, exact=False):
        raw = self.raw(parameters, exact)
        g = self.f.gutz if gutz is None else gutz
        return {m: raw.get(m, 0j) * math.exp(g * self.f.doublons(m))
                for m in oracle.even_masks(self.f.L)}

    def derivatives(self, masks, parameters=None):
        parameters = self.f.parameters if parameters is None else np.asarray(parameters)
        base = self.raw(parameters)
        result = {0: {m: complex(self.f.doublons(m)) for m in masks},
                  1: {m: 0j for m in masks}}
        for k in range(len(parameters)):
            for imaginary in (0, 1):
                estimates = []
                for epsilon in (1e-5, 5e-6):
                    plus, minus = parameters.copy(), parameters.copy()
                    delta = epsilon * (1j if imaginary else 1)
                    plus[k] += delta
                    minus[k] -= delta
                    a, b = self.raw(plus), self.raw(minus)
                    estimates.append({m: (a[m] - b[m]) / (2 * epsilon)
                                      for m in masks})
                for m in masks:
                    close("finite difference", estimates[0][m], estimates[1][m], 2e-7)
                result[2 * (k + 1) + imaginary] = {
                    m: (4 * estimates[1][m] - estimates[0][m]) / (3 * base[m])
                    for m in masks}
        return result


def run(binary, work, procs=1, trace=True, extra=None):
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "1"),
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    if trace:
        env.update(MVMC_GC_STATE_DUMP="state.dat", MVMC_GC_STATE_AMPLITUDE="1",
                   MVMC_GC_SPIN_QP_DUMP="qp.dat", MVMC_GC_SR_DUMP="sr.dat")
    if extra:
        env.update(extra)
    command = [str(binary), "-e", "namelist.def", "initial.def"]
    if procs > 1:
        command = mpi_command(procs, command)
    process = subprocess.run(command, cwd=work, env=env, text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             timeout=180)
    (work / "run.log").write_text(process.stdout)
    if process.returncode or "Finish calculation." not in process.stdout:
        raise AssertionError("{} exit {}:\n{}".format(work, process.returncode,
                                                      process.stdout[-6000:]))
    return process.stdout


def check_qp(f, work):
    dump = oracle.read_qp_dump(str(work / "qp.dat"), allow_complex_weights=True)
    if dump["count"] != f.gauss * len(f.translations) or dump["gauss"] != f.gauss:
        raise AssertionError("QP inventory")
    nodes = oracle.independent_spin_nodes(f.gauss, f.spin)
    for q, actual in enumerate(dump["nodes"]):
        if actual["q"] != q:
            raise AssertionError("QP order")
        reference = nodes[q % f.gauss]
        for key in ("cs", "cc", "ss", "weight"):
            target = reference[key] * f.weights[q // f.gauss] if key == "weight" else reference[key]
            close("QP " + key, actual[key], target, 2e-13)


def energy_table(f, wave, masks):
    return {m: 1.5 * f.doublons(m)
            - sum(t * local_of(wave, m, one_body(i, j)) for i, j, t in f.hopping)
            + f.delta * local_of(wave, m, f.pair_ops[0])
            + f.delta.conjugate() * local_of(wave, m, f.pair_ops[1]) for m in masks}


def explicit_hamiltonian(f, basis):
    """Dense Fock-space H, including its action on nodes of the trial state."""
    index = {m: i for i, m in enumerate(basis)}
    h = np.diag([1.5*f.doublons(m) for m in basis]).astype(complex)
    terms = [(-t, one_body(i, j)) for i, j, t in f.hopping]
    terms += [(f.delta, f.pair_ops[0]),
              (f.delta.conjugate(), f.pair_ops[1])]
    for column, mask in enumerate(basis):
        for coefficient, operators in terms:
            moved = apply_ops(mask, operators)
            if moved is not None:
                h[index[moved[0]], column] += coefficient*moved[1]
    if np.max(abs(h-h.conjugate().T)) > 1e-14:
        raise AssertionError("reference Hamiltonian is not Hermitian")
    return h


def energy_gradient_case(args, root):
    """Bridge same-sample SR replay to the full Rayleigh-energy derivative.

    These general-pair fixtures have no parameter-dependent exact nodes.
    At such nodes the usual |Psi|^2 SR estimator can omit both gradient and
    tangent-norm contributions; same-sample replay alone cannot detect that.
    """
    f = Fixture(args.fixture, args.gauss, args.boundary, args.k, args.shared)
    work = root / "binary"
    f.write(work, mode=0, samples=128)
    ref = Reference(f)
    wave = ref.wave()
    basis = sorted(wave)
    vector = lambda values: np.array([values[m] for m in basis])
    psi = vector(wave)
    support = psi != 0
    masks = [m for m in basis if wave[m] != 0]
    norm = np.vdot(psi, psi).real
    h = explicit_hamiltonian(f, basis)
    rayleigh = lambda v: np.vdot(v, h @ v).real / np.vdot(v, v).real
    energy = rayleigh(psi)
    local = energy_table(f, wave, masks)
    average = lambda values: sum(abs(wave[m])**2*values[m] for m in masks)/norm
    close("full Rayleigh energy", average(local), energy, 1e-12)
    derivatives = ref.derivatives(masks)
    records = []
    # The imaginary Gutzwiller coordinate is disabled by this input.
    for pi in [0] + list(range(2, 2*(len(f.parameters)+1))):
        o = derivatives[pi]
        sr = 2*(average({m: o[m]*local[m] for m in masks}).real -
                energy*average(o).real)
        differences, tangents = [], []
        for epsilon in (1e-5, 5e-6):
            if pi == 0:
                a, b = vector(ref.wave(gutz=f.gutz+epsilon)), \
                       vector(ref.wave(gutz=f.gutz-epsilon))
            else:
                k, imaginary = pi//2-1, pi % 2
                plus, minus = f.parameters.copy(), f.parameters.copy()
                delta = epsilon*(1j if imaginary else 1)
                plus[k] += delta
                minus[k] -= delta
                a, b = vector(ref.wave(plus)), vector(ref.wave(minus))
            differences.append((rayleigh(a)-rayleigh(b))/(2*epsilon))
            tangents.append((a-b)/(2*epsilon))
        tangent = (4*tangents[1]-tangents[0])/3
        node_norm = np.vdot(tangent[~support], tangent[~support]).real/norm
        if node_norm > 1e-16:
            raise AssertionError("gradient fixture has parameter-dependent nodes")
        for difference in differences:
            close("full energy gradient {}".format(pi), sr, difference, 2e-7)
        records.append(dict(parameter=pi, sr_gradient=sr,
                            full_energy_fd=differences,
                            node_tangent_norm=float(node_norm)))
    run(args.binary, work)
    replay(f, work, 0, 128)
    (root / "energy-gradient.json").write_text(
        json.dumps(records, indent=2, allow_nan=False) + "\n")


def read_parameters(work, nparameters):
    values = [float(x) for x in (work / "output/zqp_opt.dat").read_text().split()]
    if len(values) != 6 + 3 * (1 + nparameters) or not np.isfinite(values).all():
        raise AssertionError("optimized parameter inventory")
    return values[6], np.array([complex(values[9 + 3*k], values[10 + 3*k])
                               for k in range(nparameters)])


def strict_sr(path, expected_parameters):
    steps = []
    for line in path.read_text().splitlines():
        words = line.split()
        if not words:
            continue
        if words[0] == "STEP":
            header = dict(zip(words[::2], words[1::2]))
            if len(words) != 10 or set(header) != {
                    "STEP", "NPARA", "SROPTSIZE", "NPROJ", "STOREO"} or \
                    int(header["STEP"]) != len(steps) or \
                    int(header["NPARA"]) != expected_parameters or \
                    int(header["SROPTSIZE"]) != expected_parameters+1 or \
                    int(header["NPROJ"]) != 1:
                raise AssertionError("SR header inventory")
            steps.append({})
        elif words[0] == "P" and steps:
            key = int(words[1])
            if len(words) != 7 or key in steps[-1]:
                raise AssertionError("duplicate/malformed SR row")
            values = tuple(float(x) for x in words[2:])
            if not np.isfinite(values).all():
                raise AssertionError("nonfinite SR row")
            steps[-1][key] = values
        else:
            raise AssertionError("unknown SR row")
    for rows in steps:
        if set(rows) != set(range(2 * expected_parameters)):
            raise AssertionError("missing SR rows")
    return steps


def check_physical_inventory(work, mode, iterations):
    paths = [(work / "zvo_gc.dat", iterations, 3)]
    paths += [(work / ("output/zvo_out_{:03d}.dat".format(i+1)),
               1 if mode else iterations, 6)
              for i in range(iterations if mode else 1)]
    for path, count, columns in paths:
        rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
        if len(rows) != count or any(len(row) != columns for row in rows):
            raise AssertionError("physical output inventory: " + path.name)
        if not np.isfinite([[float(x) for x in row] for row in rows]).all():
            raise AssertionError("nonfinite physical output: " + path.name)


def output_orbital(row, offset, sites):
    site, spin = int(row[offset]), int(row[offset+1])
    if not 0 <= site < sites or spin not in (0, 1):
        raise AssertionError("Green orbital index outside range")
    return site + sites*spin


def replay(f, work, mode, samples, procs=1, split=None, parameters_by_step=None):
    check_qp(f, work)
    check_physical_inventory(work, mode, 2)
    split = procs if split is None else split
    groups = []
    rank_records = []
    for rank in range(procs):
        records = state_dump_records(work / ("state.dat.rank{}".format(rank)))
        rank_records.append(records)
        iterations = replay_iterations(records, samples)
        if rank % split == 0:
            groups.append(iterations)
        elif records != rank_records[rank - rank % split]:
            raise AssertionError("QP ranks disagree on sample/burn state")
    ref = Reference(f)
    all_masks = []
    for step in range(2):
        gutz, parameters = ((f.gutz, f.parameters) if parameters_by_step is None
                            else parameters_by_step[step])
        raw, wave = ref.raw(parameters), ref.wave(parameters, gutz)
        records = [r for chain in groups for r in chain[step]]
        masks = [r[3] for r in records]
        all_masks.append(masks)
        for record in records:
            if not np.isfinite(wave[record[3]]) or wave[record[3]] == 0:
                raise AssertionError("sampled reference node")
            pieces = record[4]
            words = pieces[0].split()
            order = [int(x) for x in words[3:3 + record[2]]]
            parity = sum(order[i] > order[j] for i in range(len(order))
                         for j in range(i + 1, len(order)))
            logip = complex(*(float.fromhex(x) for x in pieces[4].split()[1:]))
            close("sample amplitude", cmath.exp(logip), (-1)**parity * raw[record[3]], 1e-11)
        unique = sorted(set(masks))
        energies = energy_table(f, wave, unique)
        mean = lambda values: sum(values[m] for m in masks) / len(masks)
        e = mean(energies)
        e2 = mean({m: abs(energies[m])**2 for m in unique})
        n, n2 = np.mean([m.bit_count() for m in masks]), np.mean([m.bit_count()**2 for m in masks])
        expected = dict(energy=e, energy2=e2, relative_variance=((e2-e*e)/(e*e)).real,
                        number=n, number2=n2, variance_number=n2-n*n, sz=0, sz2=0)
        actual = parse_physical(str(work), mode, step)
        for key in expected:
            close("step {} {}".format(step, key), actual[key], expected[key])
        if mode:
            tag = "{:03d}".format(step + 1)
            key_specs = [
                ("cisajs", [(i, j) for i in range(2*f.L) for j in range(2*f.L)], 4),
                ("cisajscktalt", f.two_keys, 8),
                ("NBodyG", [f.nbody_key], 13),
                ("anomalousg", [0, 1], 5)]
            s2 = 0j
            for group, keys, offset in key_specs:
                rows = [line.split() for line in
                        (work / ("output/zvo_{}_{}.dat".format(group, tag))).read_text().splitlines()
                        if line.strip()]
                if len(rows) != len(keys):
                    raise AssertionError(group + " inventory")
                for row, key in zip(rows, keys):
                    if len(row) != offset + 2:
                        raise AssertionError(group + " column inventory")
                    if group == "cisajs":
                        parsed = tuple(output_orbital(row, k, f.L) for k in (0, 2))
                        ops = one_body(*key)
                    elif group in ("cisajscktalt", "NBodyG"):
                        start = 1 if group == "NBodyG" else 0
                        if group == "NBodyG" and int(row[0]) != len(key):
                            raise AssertionError("NBodyG operator count")
                        parsed = tuple((output_orbital(row, k, f.L),
                                        output_orbital(row, k+2, f.L))
                                       for k in range(start, offset, 4))
                        ops = tuple(op for a, b in key for op in one_body(a, b))
                    else:
                        parsed = key
                        expected_prefix = [1, 0, 0, 1, 1] if key == 0 else [0, 1, 1, 0, 0]
                        if list(map(int, row[:5])) != expected_prefix:
                            raise AssertionError("anomalous keys")
                        ops = f.pair_ops[key]
                    if parsed != key:
                        raise AssertionError(group + " keys")
                    value = complex(float(row[offset]), float(row[offset+1]))
                    target = mean({m: local_of(wave, m, ops) for m in unique})
                    close("step {} {} {}".format(step, group, key), value, target)
                    if group == "cisajscktalt" and key != f.two_keys[-1]:
                        s2 += value / 2
            expected["s2"] = s2
        else:
            derivatives = ref.derivatives(unique, parameters)
            steps = strict_sr(work / "sr.dat", 1 + len(parameters))
            if len(steps) != 2:
                raise AssertionError("SR row inventory")
            rows = steps[step]
            if set(rows) != set(derivatives):
                raise AssertionError("SR parameter inventory")
            for pi, values in derivatives.items():
                o = mean(values)
                oh = mean({m: values[m] * energies[m] for m in unique})
                target = (o.real, o.imag, oh.real, oh.imag, 2*(oh.real-e.real*o.real))
                if len(rows[pi]) != 5:
                    raise AssertionError("SR column inventory")
                for a, b in zip(rows[pi], target):
                    close("SR row {}".format(pi), a, b, 2e-7)
    return all_masks


def validation_case(args, root):
    """Corrupt actual binary output, ensuring the replay cannot pass vacuously."""
    rejected = []
    for mode in (1, 0):
        f = Fixture("L2_S1")
        work = root / ("mode" + str(mode))
        f.write(work, mode=mode, samples=8)
        run(args.binary, work)
        replay(f, work, mode, 8)
        paths = [work / "zvo_gc.dat", work / "output/zvo_out_001.dat"]
        paths += ([work / ("output/zvo_{}_001.dat".format(group)) for group in
                   ("cisajs", "cisajscktalt", "NBodyG", "anomalousg")]
                  if mode else [work / "sr.dat"])
        for path in paths:
            original = path.read_text()
            lines = original.splitlines()
            index = 1 if path.name == "sr.dat" else 0
            words = lines[index].split()
            mutations = {
                "missing": lines[:index] + lines[index+1:],
                "duplicate": lines[:index] + [lines[index]] + lines[index:],
                "extra-column": lines[:index] + [lines[index] + " 0"] + lines[index+1:]}
            for value in ("nan", "inf", "-inf"):
                mutations[value] = lines[:index] + [" ".join(words[:-1] + [value])] + lines[index+1:]
            if path.name == "sr.dat" or any(group in path.name for group in
                    ("cisajs", "NBodyG", "anomalousg")):
                wrong = list(words)
                wrong[1 if path.name == "sr.dat" else 0] = "999"
                mutations["wrong-key"] = lines[:index] + [" ".join(wrong)] + lines[index+1:]
            if "cisajs" in path.name:
                alias = list(words)
                alias[0], alias[1] = str(int(alias[0])+f.L), str(int(alias[1])-1)
                mutations["aliased-invalid-orbital"] = (
                    lines[:index] + [" ".join(alias)] + lines[index+1:])
            for label, changed in mutations.items():
                try:
                    path.write_text("\n".join(changed) + "\n")
                    try:
                        replay(f, work, mode, 8)
                    except AssertionError:
                        rejected.append("{}:{}:{}".format(mode, path.name, label))
                    else:
                        raise RuntimeError("invalid replay passed: " + str(path) + " " + label)
                finally:
                    path.write_text(original)
        replay(f, work, mode, 8)
    (root / "validation.json").write_text(json.dumps(rejected, indent=2) + "\n")


def replay_case(args, root):
    f = Fixture(args.fixture, args.gauss, args.boundary, args.k, args.shared)
    results = {}
    for mode, store, weight in [(1, 0, 0), (0, 0, 0), (0, 1, 0), (0, 1, 1)]:
        work = root / "m{}-o{}-w{}".format(mode, store, weight)
        f.write(work, mode=mode, store=store, reweight=weight, split=args.split or args.np)
        run(args.binary, work, args.np)
        results[(mode, store, weight)] = replay(f, work, mode, 128, args.np, args.split)
    if results[(0, 0, 0)] != results[(0, 1, 0)] or results[(0, 1, 0)] != results[(0, 1, 1)]:
        raise AssertionError("NStore/reweight changed the fixed-parameter chain")


def optimization_case(args, root):
    f = Fixture(args.fixture, args.gauss, args.boundary, args.k, args.shared)
    workdirs = []
    for steps in (1, 2):
        work = root / ("steps" + str(steps))
        f.write(work, mode=0, store=1, samples=1024, freeze=False, iterations=steps)
        run(args.binary, work)
        workdirs.append(work)
    gutz, updated = read_parameters(workdirs[0], len(f.parameters))
    if np.max(abs(updated-f.parameters)) == 0:
        raise AssertionError("SR did not update parameters")
    first = [r for r in state_dump_records(workdirs[0] / "state.dat.rank0") if r[0] == "SAMPLE"]
    second = replay_iterations(state_dump_records(workdirs[1] / "state.dat.rank0"), 1024)
    if first != second[0]:
        raise AssertionError("first SR iteration changed with iteration count")
    replay(f, workdirs[1], 0, 1024,
           parameters_by_step=[(f.gutz, f.parameters), (gutz, updated)])


def convergence_case(args, root):
    records = []
    for gauss in (2, 3, 4, 8, 16, 32, 64):
        f = Fixture(args.fixture, gauss, args.boundary, args.k)
        work = root / ("g" + str(gauss))
        f.write(work, samples=16)
        run(args.binary, work)
        replay(f, work, 1, 16)
        ref = Reference(f)
        raw, exact = ref.raw(), ref.raw(exact=True)
        basis = ref.basis
        a = np.array([raw[m] for m in basis])
        b = np.array([exact[m] for m in basis])
        norm = np.linalg.norm(b)
        if norm < 1e-10:
            raise AssertionError("convergence fixture has no target-S component")
        error = np.linalg.norm(a-b)/norm
        residual = oracle.spin_residual(raw, f.L, f.spin)
        # Explicit weight outside the target S eigenspace, not only <S^2>.
        eigenvalues, eigenvectors = np.linalg.eigh(oracle.spin2_matrix(f.L, basis))
        other = abs(eigenvalues-f.spin*(f.spin+1)) > 1e-8
        leakage = float(np.sum(abs((eigenvectors.conj().T @ a)[other])**2)/np.vdot(a, a).real)
        pa, pb = abs(a)**2/np.vdot(a, a).real, abs(b)**2/np.vdot(b, b).real
        sector_error = max(abs(sum(pa[i]-pb[i] for i, m in enumerate(basis)
                                   if m.bit_count() == n)) for n in range(0, 2*f.L+1, 2))
        if gauss >= 32 and max(error, residual, leakage, sector_error) > 1e-10:
            raise AssertionError("finite-S quadrature failed to converge")
        records.append(dict(gauss=gauss, relative_error=float(error), residual=residual,
                            leakage=leakage, sector_error=float(sector_error)))
    (root / "convergence.json").write_text(json.dumps(records, indent=2) + "\n")


def statistics_case(args, root):
    """Fixed design: four pilot runs, then eight fresh seeds per initial N.

    The pilot is diagnostic; the confirmation length, seeds and simultaneous
    Student-t gate below are fixed before observing either phase.  Chain means
    are independent between initial groups. Constant observables use no ESS.
    """
    close("Student-t tail convention", student_t_quantile(0.025, 7),
          2.36462425, 1e-8)
    f = Fixture(args.fixture, args.gauss, args.boundary, args.k)
    ref = Reference(f)
    wave = ref.wave()
    norm = sum(abs(z)**2 for z in wave.values())
    if norm <= 0 or not np.isfinite(norm):
        raise AssertionError("zero statistical reference")
    probability = {m: abs(z)**2/norm for m, z in wave.items() if z != 0}
    sectors = sorted(n for n in {m.bit_count() for m in probability}
                     if sum(p for m, p in probability.items() if m.bit_count() == n) > 1e-10)
    initials = sorted(set([min(sectors), max(sectors)]))
    design = dict(pilot_seeds=[11, 23, 37, 53], confirm_seeds=[
        [101 + 37*i + 1009*j for i in range(8)] for j in range(len(initials))],
        pilot_samples=4096, confirm_samples=16384, warmup=2048,
        initials=initials, family_alpha=0.001, block_sizes=[128, 256, 512])
    (root / "design.json").write_text(json.dumps(design, indent=2) + "\n")
    rows = []
    reference = None
    for phase in ("pilot", "confirm"):
        for j, initial in enumerate(initials):
            seeds = (design["pilot_seeds"][j::len(initials)] if phase == "pilot"
                     else design["confirm_seeds"][j])
            for seed in seeds:
                samples = design[phase + "_samples"]
                work = root / "{}-N{}-seed{}".format(phase, initial, seed)
                f.write(work, samples=samples, iterations=1, seed=seed,
                        initial=initial, warmup=design["warmup"])
                run(args.binary, work)
                check_qp(f, work)
                masks = [r[3] for r in state_dump_records(work / "state.dat.rank0")
                         if r[0] == "SAMPLE"]
                if len(masks) != samples or any(probability.get(m, 0) == 0 for m in masks):
                    raise AssertionError("statistical sample inventory/support")
                energy = energy_table(f, wave, probability)
                series = {
                    "N": np.array([m.bit_count() for m in masks], dtype=float),
                    "N2": np.array([m.bit_count()**2 for m in masks], dtype=float),
                    "energy_real": np.array([energy[m].real for m in masks]),
                    "energy_imag": np.array([energy[m].imag for m in masks])}
                series.update({"P{}".format(n): np.array([m.bit_count() == n for m in masks],
                                                          dtype=float) for n in sectors})
                reference = {
                    "N": sum(m.bit_count()*p for m, p in probability.items()),
                    "N2": sum(m.bit_count()**2*p for m, p in probability.items()),
                    "energy_real": sum(energy[m].real*p for m, p in probability.items()),
                    "energy_imag": sum(energy[m].imag*p for m, p in probability.items())}
                reference.update({"P{}".format(n): sum(p for m, p in probability.items()
                    if m.bit_count() == n) for n in sectors})
                output = parse_physical(str(work), 1, 0)
                close("stat N replay", output["number"], series["N"].mean())
                close("stat N2 replay", output["number2"], series["N2"].mean())
                close("stat energy replay", output["energy"],
                      series["energy_real"].mean() + 1j*series["energy_imag"].mean())
                row = dict(phase=phase, initial=initial, seed=seed, observables={})
                for key, values in series.items():
                    var = float(values.var())
                    blocks = []
                    for size in design["block_sizes"]:
                        means = values.reshape(-1, size).mean(axis=1)
                        se = float(means.std(ddof=1)/math.sqrt(len(means)))
                        g = max(1., samples*se*se/var) if var > 1e-28 else None
                        blocks.append(dict(size=size, se=se, g=g,
                                           ESS=samples/g if g is not None else None))
                    row["observables"][key] = dict(mean=float(values.mean()),
                                                   blocks=blocks)
                rows.append(row)
        # Preserve the pilot and its block curves before starting confirmation.
        (root / (phase + ".json")).write_text(json.dumps(rows, indent=2, allow_nan=False) + "\n")
    family = len(initials) * len(reference)
    # The helper accepts the upper-tail probability, not the CDF value.
    critical = student_t_quantile(design["family_alpha"]/(2*family), 7)
    comparisons = []
    for initial in initials:
        group = [row for row in rows if row["phase"] == "confirm" and row["initial"] == initial]
        for key, exact in reference.items():
            values = np.array([row["observables"][key]["mean"] for row in group])
            se = float(values.std(ddof=1)/math.sqrt(len(values)))
            error = abs(float(values.mean()) - exact)
            comparisons.append(dict(initial=initial, observable=key, exact=exact,
                mean=float(values.mean()), se=se, critical=critical, error=error))
            if error > critical*se + 1e-10*(1+abs(exact)):
                raise AssertionError("confirmation bias: " + str(comparisons[-1]))
    (root / "statistics.json").write_text(json.dumps(dict(design=design,
        comparisons=comparisons, chains=rows), indent=2, allow_nan=False) + "\n")


def negative_case(args, root):
    cases = [
        ("negative", "L4_S1_general", {"NSPStot": -1}, "0 <= NSPStot <= Nsite/2"),
        ("above_max", "L4_S1_general", {"NSPStot": 3}, "0 <= NSPStot <= Nsite/2"),
        ("int_max", "L4_S1_general", {"NSPStot": 2147483647}, "0 <= NSPStot <= Nsite/2"),
        ("fraction", "L4_S1_general", {"NSPStot": "1.5"}, "must be an exact integer"),
        ("nan", "L4_S1_general", {"NSPStot": "nan"}, "must be an exact integer"),
        ("inf", "L4_S1_general", {"NSPStot": "inf"}, "must be an exact integer"),
        ("initial_low", "L4_S1_general", {"NGCInitNelec": 0}, "cannot carry total spin"),
        ("initial_high", "L4_S1_general", {"NGCInitNelec": 8}, "cannot carry total spin"),
        ("initial_triplet_parity", "L4_S1_triplet_stat", {"NGCInitNelec": 4}, "cannot carry total spin"),
        ("singlet", "L4_S1_singlet", {}, "has no total spin"),
        ("lanczos", "L4_S1_general", {"NLanczosMode": 1, "NLanczosStep": 2,
             "NLanczosEstimatorMode": 1}, "NLanczosStep=1")]
    for order in ("early", "late"):
        for name, fixture, options, message in cases:
            f = Fixture(fixture)
            work = root / (name + "-" + order)
            f.write(work, samples=8, iterations=1, split=args.np)
            for key, value in options.items():
                option(work, key, value)
            path = work / "modpara.def"
            text = path.read_text()
            text = re.sub(r"(?m)^NGrandCanonical\s+[^\n]+\n", "", text)
            marker = "NVMCCalMode "
            path.write_text(text.replace(marker, "NGrandCanonical 1\n" + marker)
                            if order == "early" else text + "NGrandCanonical 1\n")
            if name == "lanczos":
                # Isolate the unsupported Lanczos order from other features.
                names = work / "namelist.def"
                names.write_text("\n".join(line for line in names.read_text().splitlines()
                    if not line.startswith(("OneBodyG ", "TwoBodyG ", "NBodyG ",
                                            "AnomalousG ", "AnomalousTerm "))) + "\n")
            command = [str(args.binary), "-e", "namelist.def", "initial.def"]
            if args.np > 1:
                command = mpi_command(args.np, command)
            process = subprocess.run(command, cwd=work, text=True, timeout=60,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                env=dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1"))
            (work / "run.log").write_text(process.stdout)
            if process.returncode <= 0 or message not in process.stdout or re.search(
                    r"AddressSanitizer|runtime error:|SIGSEGV|Segmentation fault", process.stdout):
                raise AssertionError("negative {}:\n{}".format(name, process.stdout[-5000:]))


def sample_case(args, root):
    for name, total_spin in (("S1", 1), ("S2", 2)):
        for mode in ("measurement", "optimization"):
            work = root / (name + "-" + mode)
            shutil.copytree(args.samples_dir / name, work)
            if mode == "optimization":
                shutil.copyfile(work / "modpara_opt.def", work / "modpara.def")
                shutil.copyfile(work / "namelist_opt.def", work / "namelist.def")
                # Keep the supplied 20-step sample, using its original filenames.
                text = (work / "namelist.def").read_text().replace(
                    "modpara_opt.def", "modpara.def")
                (work / "namelist.def").write_text(text)
            run(args.binary, work, trace=False)
            numbers = np.loadtxt(work / "zvo_gc.dat", ndmin=2)
            output = np.loadtxt(work / "output/zvo_out_001.dat", ndmin=2)
            if not np.isfinite(numbers).all() or not np.isfinite(output).all():
                raise AssertionError("nonfinite public sample")
            L = 4 if name == "S1" else 5
            if np.any(numbers[:, 0] < 2*total_spin-1e-10) or \
                    np.any(numbers[:, 0] > 2*L-2*total_spin+1e-10):
                raise AssertionError("sample N outside spin support")
            if mode == "measurement":
                rows = np.loadtxt(work / "output/zvo_cisajscktalt_001.dat", ndmin=2)
                if len(rows) != 2*L*L+1:
                    raise AssertionError("sample S2 row inventory")
                s2 = np.sum(rows[:-1, 8] + 1j*rows[:-1, 9])/2
                close("public sample S2", s2, total_spin*(total_spin+1), 1e-6)
            else:
                if len(output) != 20:
                    raise AssertionError("sample optimization step count")
                parameters = np.loadtxt(work / "output/zqp_opt.dat")
                if not np.isfinite(parameters).all():
                    raise AssertionError("nonfinite optimized sample")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("replay", "optimization", "convergence",
                                        "statistics", "negative", "sample", "validation",
                                        "energy_gradient"))
    parser.add_argument("--binary", type=Path,
                        default=Path("../../src/mVMC/vmc.out"))
    parser.add_argument("--fixture", default="L4_S1_general")
    parser.add_argument("--gauss", type=int, default=8)
    parser.add_argument("--boundary", choices=("pbc", "apbc"), default="pbc")
    parser.add_argument("--k", type=int)
    parser.add_argument("--shared", action="store_true")
    parser.add_argument("--np", type=int, default=1)
    parser.add_argument("--split", type=int)
    parser.add_argument("--samples-dir", type=Path,
                        default=Path(__file__).resolve().parents[2] /
                        "samples/GrandCanonical/FiniteSpin")
    args = parser.parse_args()
    args.binary = args.binary.resolve()
    Path("work").mkdir(exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="gc-fs-" + args.mode + "-", dir="work")).resolve()
    globals()[args.mode + "_case"](args, root)
    print("finite-S {} {} PASS: {}".format(args.mode, args.fixture, root))


if __name__ == "__main__":
    main()
