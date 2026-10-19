"""Finite total spin S>0: production probe and chains versus the oracle.

The public S>0 input gate is closed (Task 1/2 of the finite-S plan), so the
production code is driven by the Testing-only gc_finite_spin_unit binary.
Reference values come only from grandcanonical_spin_oracle (independent
Fock-space and rational evaluation); the binary never builds them.

  probe --fixture NAME --gauss G   amplitudes, QP coefficients, log and
                                   zero-point derivatives, Green functions
  chain --fixture NAME --gauss G   production sampler histogram (many
                                   independent chains) against |Psi|^2
"""
from __future__ import print_function

import argparse
import math
import os
import subprocess
import sys

import grandcanonical_spin_oracle as oracle
from grandcanonical_exact_oracle import apply_ops, one_body

# Raw amplitude gate of plan section 6.2 (initial value).
AMPLITUDE_TOLERANCE = 1.0e-11
# Coefficients against an independent leggauss/legval.
COEFFICIENT_TOLERANCE = 2.0e-13
# Derivative gate of plan section 6.2. The oracle uses the exact affine
# difference (cross-checked against eps 1e-5/5e-6 finite differences in the
# oracle self-test), so the gate is not limited by 1/eps roundoff.
DERIVATIVE_TOLERANCE = 2.0e-7
# Ratios of production amplitudes against the oracle wave.
GREEN_TOLERANCE = 1.0e-9
# Ratio checks (log derivatives, Green functions) use only configurations
# with |Psi| >= CONDITION_CUT * max|Psi|. Smaller values are quadrature or
# cancellation residue (e.g. 1e-16 of a near-triplet N=4 sector) whose
# ratios are ill-conditioned; their amplitudes are still checked absolutely.
CONDITION_CUT = 1.0e-8

GUTZWILLER = {"L2_S1": 0.0, "L4_S2": -0.3, "L4_S1_triplet": 0.0,
              "L4_S1_triplet_stat": 0.0, "L4_S1_triplet_probe": 0.0,
              "L4_S1_general": 0.4, "L5_S2": 0.25, "L6_S3": 0.0,
              "L4_S1_near": -0.2}


def write_fixture(path, name, gauss):
    F, spin = oracle.finite_spin_fixture(name)
    nsite = len(F)
    with open(path, "w") as stream:
        stream.write("GC_FINITE_SPIN_FIXTURE_V1\n")
        stream.write("{} {} {} {}\n".format(nsite, spin, gauss,
                                            float(GUTZWILLER[name]).hex()))
        for row in F:
            for value in row:
                stream.write("{} {}\n".format(float(value.real).hex(),
                                              float(value.imag).hex()))
    return F, spin


def run(command):
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "1")
    process = subprocess.run(command, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             universal_newlines=True, env=environment)
    if process.returncode != 0:
        raise AssertionError("{} failed:\n{}".format(
            " ".join(command), process.stdout[-4000:]))
    return process.stdout


def hexc(re_text, im_text):
    return complex(float.fromhex(re_text), float.fromhex(im_text))


def parse_probe(path):
    data = {"Q": [], "A": {}, "D": {}, "Z": {}, "G1": {}, "PA": {},
            "PR": {}, "G2": {}}
    with open(path) as stream:
        header = stream.readline().split()
        if header[0] != "GC_FINITE_SPIN_PROBE_V1":
            raise AssertionError("bad probe header {}".format(header))
        data["header"] = [int(word) for word in header[1:]]
        for line in stream:
            w = line.split()
            kind = w[0]
            if kind == "Q":
                data["Q"].append({"q": int(w[1]),
                                  "cs": float.fromhex(w[2]),
                                  "cc": float.fromhex(w[3]),
                                  "ss": float.fromhex(w[4]),
                                  "weight": hexc(w[5], w[6])})
            elif kind == "A":
                data["A"][int(w[1])] = (int(w[2]), int(w[3]),
                                        hexc(w[4], w[5]))
            elif kind == "D":
                data["D"][(int(w[1]), int(w[2]))] = (hexc(w[3], w[4]),
                                                     hexc(w[5], w[6]))
            elif kind == "Z":
                data["Z"][(int(w[1]), int(w[2]))] = hexc(w[3], w[4])
            elif kind in ("G1", "PA", "PR"):
                data[kind][(int(w[1]), int(w[2]), int(w[3]))] = \
                    hexc(w[4], w[5])
            elif kind == "G2":
                data[kind][tuple(int(x) for x in w[1:6])] = hexc(w[6], w[7])
            else:
                raise AssertionError("unknown probe row {}".format(kind))
    return data


def doublons(mask, nsite):
    low = (1 << nsite) - 1
    return oracle.popcount((mask & low) & (mask >> nsite))


def production_nodes(data, gauss):
    """Production coefficients as oracle nodes (validated separately)."""
    nodes = []
    for row in data["Q"]:
        if row["weight"].imag != 0.0:
            raise AssertionError("complex QP weight with identity "
                                 "translation")
        cc, ss = row["cc"], row["ss"]
        nodes.append({"spidx": row["q"] % gauss, "cs": row["cs"],
                      "cc": cc, "ss": ss, "c": math.sqrt(cc),
                      "s": math.sqrt(ss), "weight": row["weight"].real})
    return nodes


def check_coefficients(data, gauss, spin):
    reference = oracle.independent_spin_nodes(gauss, spin)
    worst = 0.0
    if len(data["Q"]) != gauss:
        raise AssertionError("QP count {} != G {}".format(len(data["Q"]),
                                                          gauss))
    for row, ref in zip(data["Q"], reference):
        for key in ("cs", "cc", "ss"):
            worst = max(worst, abs(row[key] - ref[key]))
        worst = max(worst, abs(row["weight"] - ref["weight"]))
    if worst > COEFFICIENT_TOLERANCE:
        raise AssertionError("production a_g(S) differ from leggauss/legval "
                             "by {}".format(worst))
    return worst


def oracle_wave(F, nodes, spin, gutzwiller, exact=False):
    nsite = len(F)
    masked = oracle.masked_finite_spin(F, nodes, spin,
                                       exact_arithmetic=exact)
    if exact:
        return masked
    return dict((mask, value * math.exp(gutzwiller * doublons(mask, nsite)))
                for mask, value in masked.items())


def probe_case(args):
    workdir = os.path.abspath(args.work)
    if not os.path.isdir(workdir):
        os.makedirs(workdir)
    fixture = os.path.join(workdir, "fixture.txt")
    output = os.path.join(workdir, "probe.txt")
    F, spin = write_fixture(fixture, args.fixture, args.gauss)
    nsite = len(F)
    gutzwiller = GUTZWILLER[args.fixture]
    run([args.binary, "--probe", fixture, output])
    data = parse_probe(output)
    symmetry = oracle.pair_symmetry(F)
    header = data["header"]
    expected = [nsite, spin, args.gauss, 1, int(symmetry == "triplet"),
                int(symmetry == "singlet")]
    if header != expected:
        raise AssertionError("probe header {} != {}".format(header, expected))
    worst_coefficient = check_coefficients(data, args.gauss, spin)

    # Amplitudes: exact rational reference with the production coefficients,
    # proven zeros exact in production, and no legitimate nonzero dropped.
    nodes = production_nodes(data, args.gauss)
    exact = oracle.masked_finite_spin(F, nodes, spin, exact_arithmetic=True)
    reference = dict((mask, oracle.qfloat(value))
                     for mask, value in exact.items())
    scale = max(abs(value) for value in reference.values())
    if set(data["A"]) != set(reference):
        raise AssertionError("probe masks differ from the M=0 basis")
    worst_amplitude = 0.0
    zeros = nonzeros = 0
    for mask, (count, singles, actual) in data["A"].items():
        if count != oracle.popcount(mask) or \
                singles != oracle.singles(mask, nsite):
            raise AssertionError("probe N/singles mismatch at {}".format(mask))
        proven = oracle.finite_spin_zero(mask, nsite, spin, symmetry)
        if proven:
            zeros += 1
            if actual != 0:
                raise AssertionError("proven zero {} is {} in production"
                                     .format(mask, actual))
            continue
        error = abs(actual - reference[mask])
        worst_amplitude = max(worst_amplitude, error)
        if error > AMPLITUDE_TOLERANCE * (1.0 + scale):
            raise AssertionError("amplitude {} differs: {} vs {}".format(
                mask, actual, reference[mask]))
        if exact[mask] != oracle.ZERO:
            nonzeros += 1
            if actual == 0:
                raise AssertionError("nonzero amplitude {} dropped".format(
                    mask))

    # Log derivatives at every nonzero configuration (finite differences of
    # the masked state, reclassifying F at each displaced point).
    wave = oracle_wave(F, nodes, spin, gutzwiller)
    nonzero = sorted(mask for mask, row in data["A"].items()
                     if row[2] != 0 and
                     abs(reference[mask]) >= CONDITION_CUT * scale)
    ill_conditioned = sum(1 for row in data["A"].values() if row[2] != 0) - \
        len(nonzero)
    worst_derivative = 0.0
    derivative_cache = {}

    def derivative(k, imaginary):
        key = (k, imaginary)
        if key not in derivative_cache:
            derivative_cache[key] = oracle.affine_spin_derivative(
                F, nodes, spin, k // nsite, k % nsite, imaginary)
        return derivative_cache[key]

    for k in range(nsite * nsite):
        for imaginary in (0, 1):
            values = derivative(k, imaginary)
            for mask in nonzero:
                actual = data["D"][(mask, k)][imaginary]
                expect = values[mask] / reference[mask]
                error = abs(actual - expect)
                worst_derivative = max(worst_derivative,
                                       error / (1.0 + abs(expect)))
                if error > DERIVATIVE_TOLERANCE * (1.0 + abs(expect)):
                    raise AssertionError(
                        "log derivative mask={} k={} im={}: {} vs {}".format(
                            mask, k, imaginary, actual, expect))
    # Zero-point derivatives of pure-triplet parity zeros.
    worst_zero_point = 0.0
    largest_zero_point = 0.0
    for (mask, k), actual in data["Z"].items():
        expect = derivative(k, 0)[mask]
        largest_zero_point = max(largest_zero_point, abs(expect))
        error = abs(actual - expect)
        worst_zero_point = max(worst_zero_point, error / (1.0 + abs(expect)))
        if error > DERIVATIVE_TOLERANCE * (1.0 + abs(expect)):
            raise AssertionError("zero-point derivative mask={} k={}: {} vs {}"
                                 .format(mask, k, actual, expect))
    expected_z = set(mask for mask in reference if symmetry == "triplet" and
                     oracle.configuration_allowed(mask, nsite, spin) and
                     oracle.triplet_sector_forbidden(oracle.popcount(mask),
                                                     spin))
    if set(mask for mask, unused in data["Z"]) != expected_z or \
            len(data["Z"]) != len(expected_z) * nsite * nsite:
        raise AssertionError("zero-point derivative rows differ from the "
                             "triplet parity zeros")
    if data["Z"] and largest_zero_point < 1.0e-3:
        raise AssertionError("zero-point derivatives vanish: the fixture "
                             "does not test the symmetry-breaking direction")

    # Green functions: ratios of the oracle wave including virtual
    # configurations with fewer than 2S singles (exact zeros).
    worst_green = 0.0
    exact_virtual_zeros = 0
    for group, build in (("G1", lambda key: one_body(key[1], key[2])),
                         ("PA", lambda key: (("c", key[1]), ("c", key[2]))),
                         ("PR", lambda key: (("a", key[1]), ("a", key[2]))),
                         ("G2", lambda key: one_body(key[1], key[2]) +
                          one_body(key[3], key[4]))):
        for key, actual in data[group].items():
            mask = key[0]
            if mask not in nonzero:
                continue
            moved = apply_ops(mask, build(key))
            if moved is None or moved[0] not in wave:
                expect = 0.0j
                target_zero = True
            else:
                target = moved[0]
                target_zero = exact.get(target, oracle.ZERO) == oracle.ZERO \
                    and oracle.finite_spin_zero(target, nsite, spin, symmetry)
                expect = moved[1] * (wave[target] / wave[mask]).conjugate()
            if target_zero:
                if actual != 0:
                    raise AssertionError("{} {} reaches a proven zero but is "
                                         "{}".format(group, key, actual))
                exact_virtual_zeros += 1
                continue
            error = abs(actual - expect)
            worst_green = max(worst_green, error / (1.0 + abs(expect)))
            if error > GREEN_TOLERANCE * (1.0 + abs(expect)):
                raise AssertionError("{} {}: {} vs {}".format(group, key,
                                                              actual, expect))
    print("probe {} G={} S={} ({}): {} proven zeros, {} nonzero; worst "
          "coefficient {:.2e}, amplitude {:.2e}, log derivative {:.2e}, "
          "zero-point derivative {:.2e} (largest {:.3g}, {} rows), Green "
          "{:.2e} ({} exact virtual zeros); {} ill-conditioned masks only "
          "checked absolutely PASS".format(
              args.fixture, args.gauss, spin, symmetry, zeros, nonzeros,
              worst_coefficient, worst_amplitude, worst_derivative,
              worst_zero_point, largest_zero_point, len(data["Z"]),
              worst_green, exact_virtual_zeros, ill_conditioned))


# ---------------------------------------------------------------------------
# Production chains
# ---------------------------------------------------------------------------
# Fixed in advance (not tuned on results): seeds, chain count and length.
CHAIN_SEEDS = tuple(1009 + 7919 * k for k in range(16))


def student_t_quantile(probability, dof):
    """Upper quantile of Student t by bisection on the regularized beta."""
    def cdf(t):
        x = dof / (dof + t * t)
        return 1.0 - 0.5 * _betainc(0.5 * dof, 0.5, x)
    low, high = 0.0, 100.0
    for unused in range(200):
        middle = 0.5 * (low + high)
        if cdf(middle) < 1.0 - probability:
            low = middle
        else:
            high = middle
    return high


def _betainc(a, b, x):
    """Regularized incomplete beta I_x(a,b) by Lentz continued fraction."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_front = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) +
                 a * math.log(x) + b * math.log(1.0 - x))
    if x > (a + 1.0) / (a + b + 2.0):
        return 1.0 - _betainc(b, a, 1.0 - x)
    tiny = 1.0e-300
    c, d = 1.0, 1.0 - (a + b) * x / (a + 1.0)
    d = 1.0 / (d if abs(d) > tiny else tiny)
    result = d
    for m in range(1, 400):
        numerator = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
        d = 1.0 + numerator * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + numerator / c
        c = c if abs(c) > tiny else tiny
        result *= d * c
        numerator = -(a + m) * (a + b + m) * x / \
            ((a + 2 * m) * (a + 2 * m + 1))
        d = 1.0 + numerator * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + numerator / c
        c = c if abs(c) > tiny else tiny
        delta = d * c
        result *= delta
        if abs(delta - 1.0) < 1.0e-15:
            break
    return math.exp(log_front) * result / a


def chain_case(args):
    workdir = os.path.abspath(args.work)
    if not os.path.isdir(workdir):
        os.makedirs(workdir)
    fixture = os.path.join(workdir, "fixture.txt")
    F, spin = write_fixture(fixture, args.fixture, args.gauss)
    nsite = len(F)
    gutzwiller = GUTZWILLER[args.fixture]
    probe = os.path.join(workdir, "probe.txt")
    run([args.binary, "--probe", fixture, probe])
    data = parse_probe(probe)
    check_coefficients(data, args.gauss, spin)
    nodes = production_nodes(data, args.gauss)
    wave = oracle_wave(F, nodes, spin, gutzwiller)
    weights = dict((mask, abs(value) ** 2) for mask, value in wave.items()
                   if value != 0)
    total = sum(weights.values())
    target = dict((mask, value / total) for mask, value in weights.items())
    seeds = CHAIN_SEEDS[:args.chains]
    frequencies = []
    for seed in seeds:
        output = os.path.join(workdir, "chain-{}.txt".format(seed))
        run([args.binary, "--chain", fixture, str(seed), str(args.sweeps),
             output])
        counts = {}
        with open(output) as stream:
            header = stream.readline().split()
            if header[0] != "GC_FINITE_SPIN_CHAIN_V1" or \
                    int(header[2]) != args.sweeps:
                raise AssertionError("bad chain header {}".format(header))
            for line in stream:
                word = line.split()
                counts[int(word[1])] = int(word[2])
        outside = [mask for mask in counts if mask not in target]
        if outside:
            raise AssertionError("chain {} visited zero-amplitude masks {}"
                                 .format(seed, outside[:8]))
        frequencies.append(dict((mask, counts.get(mask, 0) /
                                 float(args.sweeps)) for mask in target))
    # Between-chain standard errors; Bonferroni over the support masks and
    # the N sectors, Student t with chains-1 degrees of freedom.
    groups = dict((("mask", mask), [mask]) for mask in target)
    for mask in target:
        groups.setdefault(("N", oracle.popcount(mask)), []).append(mask)
    alpha = 1.0e-3
    critical = student_t_quantile(alpha / (2.0 * len(groups)),
                                  len(seeds) - 1)
    worst = 0.0
    for key, masks in sorted(groups.items()):
        expect = sum(target[mask] for mask in masks)
        values = [sum(freq[mask] for mask in masks) for freq in frequencies]
        mean = sum(values) / len(values)
        variance = sum((value - mean) ** 2 for value in values) / \
            (len(values) - 1)
        error = math.sqrt(variance / len(values))
        if error == 0.0:
            if abs(mean - expect) > 1.0e-12:
                raise AssertionError("{} frequency {} without spread, "
                                     "expected {}".format(key, mean, expect))
            continue
        score = abs(mean - expect) / error
        worst = max(worst, score)
        if score > critical:
            raise AssertionError("{} frequency {:.6f} +- {:.6f}, expected "
                                 "{:.6f} (|t|={:.2f} > {:.2f})".format(
                                     key, mean, error, expect, score,
                                     critical))
        if key[0] == "mask" and min(values) == 0.0 and expect > 0.01:
            raise AssertionError("support mask {} unvisited by a chain"
                                 .format(key[1]))
    sectors = sorted(set(oracle.popcount(mask) for mask in target))
    print("chain {} G={} S={}: {} chains x {} sweeps, {} support masks, N "
          "sectors {}, worst |t| {:.2f} <= {:.2f} PASS".format(
              args.fixture, args.gauss, spin, len(seeds), args.sweeps,
              len(target), sectors, worst, critical))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("probe", "chain"))
    parser.add_argument("binary")
    parser.add_argument("--fixture", required=True,
                        choices=sorted(GUTZWILLER))
    parser.add_argument("--gauss", type=int, required=True)
    parser.add_argument("--work", required=True)
    parser.add_argument("--chains", type=int, default=16)
    parser.add_argument("--sweeps", type=int, default=20000)
    args = parser.parse_args()
    if args.mode == "probe":
        probe_case(args)
    else:
        if args.chains < 4 or args.chains > len(CHAIN_SEEDS):
            parser.error("4 <= --chains <= {}".format(len(CHAIN_SEEDS)))
        chain_case(args)


if __name__ == "__main__":
    sys.exit(main())
