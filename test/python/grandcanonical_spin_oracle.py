"""Independent Fock-space oracle for grand-canonical S=0 spin projection.

The anti-parallel pair state |phi_F> = exp(sum_ij F_ij c+_{i up} c+_{j dn})|0>
is projected in two independent ways:

* finite quadrature: sum_g a_g <x|exp(-i beta_g S^y)|phi_F> for M=0 states x,
  where the rotation is applied in Fock space (c+_{i up} -> c c+_{i up} +
  s c+_{i dn}, c+_{i dn} -> -s c+_{i up} + c c+_{i dn}, c=cos(beta/2),
  s=sin(beta/2)).  A rational reference evaluates the production 4-block
  matrices with the dumped binary64 coefficients and is cross-checked
  against the Fock rotation.
* exact projection: P0 = prod_s [1 - S^2/(s(s+1))] applied with rational
  arithmetic in each M=0 particle-number sector.

Only the fermion conventions of grandcanonical_exact_oracle (apply_ops,
one_body, creation order) are shared with the other oracles.  Orbital o is
site o % L with spin o // L, and |x> = c+_{o1} ... c+_{on}|0> with o1 < ... <
on, so the amplitude of x is Pf(A_x) over ascending orbitals.
"""
from __future__ import print_function

import argparse
import math
import sys
from collections import deque
from fractions import Fraction

import numpy

from grandcanonical_exact_oracle import apply_ops, one_body, pfaffian


ZERO = (Fraction(0), Fraction(0))
ONE = (Fraction(1), Fraction(0))
# Raw amplitude agreement between the rational 4-block reference and the
# independent floating-point Fock rotation.
FOCK_TOLERANCE = 2.0e-12
# Independent Gauss-Legendre coefficients versus the production dump.
COEFFICIENT_TOLERANCE = 2.0e-13


# ---------------------------------------------------------------------------
# Small exact complex arithmetic over Q + iQ
# ---------------------------------------------------------------------------
def rational(value):
    """Exact Q+iQ value of a finite binary64-compatible complex number."""
    value = complex(value)
    if not (math.isfinite(value.real) and math.isfinite(value.imag)):
        raise ValueError("nonfinite oracle input")
    return Fraction(value.real), Fraction(value.imag)


def qadd(first, second):
    return first[0] + second[0], first[1] + second[1]


def qsub(first, second):
    return first[0] - second[0], first[1] - second[1]


def qmul(first, second):
    return (first[0] * second[0] - first[1] * second[1],
            first[0] * second[1] + first[1] * second[0])


def qscale(value, scalar):
    return value[0] * scalar, value[1] * scalar


def qfloat(value):
    return complex(float(value[0]), float(value[1]))


# ---------------------------------------------------------------------------
# Basis, masks and the S^2 operator
# ---------------------------------------------------------------------------
def popcount(value):
    return bin(value).count("1")


def spin_counts(mask, nsite):
    low = (1 << nsite) - 1
    return popcount(mask & low), popcount(mask >> nsite)


def is_balanced(mask, nsite):
    up, down = spin_counts(mask, nsite)
    return up == down


def even_masks(nsite):
    return tuple(mask for mask in range(1 << (2 * nsite))
                 if popcount(mask) % 2 == 0)


def balanced_basis(nsite):
    """All M=0 masks (every even particle number) in ascending order."""
    if nsite <= 0:
        raise ValueError("positive site count required")
    return tuple(mask for mask in range(1 << (2 * nsite))
                 if is_balanced(mask, nsite))


def spin2_matrix(nsite, basis):
    """S^2 = S-S+ on M=0 states, built from fermion operators only."""
    index = dict((mask, k) for k, mask in enumerate(basis))
    out = numpy.zeros((len(basis), len(basis)), dtype=int)
    for col, mask in enumerate(basis):
        if not is_balanced(mask, nsite):
            raise ValueError("S^2 basis must have M=0")
        for i in range(nsite):
            middle = apply_ops(mask, one_body(i, i + nsite))
            if middle is None:
                continue
            for j in range(nsite):
                final = apply_ops(middle[0], one_body(j + nsite, j))
                if final is not None:
                    out[index[final[0]], col] += middle[1] * final[1]
    return out


# ---------------------------------------------------------------------------
# Pfaffians of every principal submatrix by subset recursion
# ---------------------------------------------------------------------------
def _subset_pfaffians(matrix, norbital, add, mul, neg, zero, one):
    """Pf(A_x) for every even mask x, expanded along the lowest orbital.

    Pf(S) = sum_k (-1)^(k+1) A[s0][sk] Pf(S without s0, sk) for the sorted
    members s0 < s1 < ...; this is the expansion of grandcanonical
    exact_oracle.pfaffian applied to the ascending principal submatrix.
    """
    values = {0: one}
    for mask in range(1, 1 << norbital):
        if popcount(mask) % 2:
            continue
        members = [orbital for orbital in range(norbital)
                   if (mask >> orbital) & 1]
        first = members[0]
        rest = mask & ~(1 << first)
        total = zero
        for k in range(1, len(members)):
            other = members[k]
            term = mul(matrix[first][other], values[rest & ~(1 << other)])
            total = add(total, term if k % 2 else neg(term))
        values[mask] = total
    return values


def float_pfaffians(matrix):
    return _subset_pfaffians(
        matrix, len(matrix), lambda a, b: a + b, lambda a, b: a * b,
        lambda a: -a, 0.0j, 1.0 + 0.0j)


def rational_pfaffians(matrix):
    return _subset_pfaffians(
        matrix, len(matrix), qadd, qmul, lambda a: (-a[0], -a[1]),
        ZERO, ONE)


def _as_rows(F):
    rows = [[complex(value) for value in row] for row in F]
    nsite = len(rows)
    if nsite == 0 or any(len(row) != nsite for row in rows):
        raise ValueError("F must be a nonempty square matrix")
    for row in rows:
        for value in row:
            if not (math.isfinite(value.real) and math.isfinite(value.imag)):
                raise ValueError("nonfinite pair amplitude")
    return rows


def antiparallel_matrix(F, convert=complex, zero=0.0j):
    rows = _as_rows(F)
    nsite = len(rows)
    matrix = [[zero] * (2 * nsite) for unused in range(2 * nsite)]
    for i in range(nsite):
        for j in range(nsite):
            value = convert(rows[i][j])
            matrix[i][j + nsite] = value
            matrix[j + nsite][i] = ((-value[0], -value[1])
                                    if isinstance(value, tuple) else -value)
    return matrix


def pair_wave(F):
    """<x|phi_F> for every even mask; M != 0 entries are exactly zero."""
    rows = _as_rows(F)
    nsite = len(rows)
    values = float_pfaffians(antiparallel_matrix(rows))
    return dict((mask, values[mask] if is_balanced(mask, nsite) else 0.0j)
                for mask in even_masks(nsite))


def _rational_pair_wave(F):
    rows = _as_rows(F)
    nsite = len(rows)
    values = rational_pfaffians(antiparallel_matrix(rows, rational, ZERO))
    return dict((mask, values[mask] if is_balanced(mask, nsite) else ZERO)
                for mask in even_masks(nsite))


# ---------------------------------------------------------------------------
# Fock-space spin rotation (independent of the production 4-block form)
# ---------------------------------------------------------------------------
def _rotated_creators(orbital, nsite, c, s):
    site, spin = orbital % nsite, orbital // nsite
    up, down = site, site + nsite
    if spin == 0:
        return ((up, c), (down, s))
    return ((up, -s), (down, c))


def rotate_state(mask, nsite, c, s):
    """R(beta)|mask> as {mask: amplitude} using creation-operator images."""
    occupied = [orbital for orbital in range(2 * nsite)
                if (mask >> orbital) & 1]
    result = {}
    choices = [_rotated_creators(orbital, nsite, c, s) for orbital in occupied]

    def expand(position, operators, amplitude):
        if amplitude == 0:
            return
        if position == len(choices):
            created = apply_ops(0, tuple(("c", target)
                                         for target in operators))
            if created is None:
                return
            target, sign = created
            result[target] = result.get(target, 0.0) + sign * amplitude
            return
        for target, factor in choices[position]:
            if target in operators:
                continue
            expand(position + 1, operators + (target,), amplitude * factor)

    expand(0, (), 1.0)
    return result


def rotate_fock(wave, nsite, c, s):
    """Apply R(beta) to a Fock-space wave function (all even masks)."""
    for value in (c, s):
        if not math.isfinite(value):
            raise ValueError("nonfinite rotation")
    result = dict((mask, 0.0j) for mask in even_masks(nsite))
    for mask, amplitude in wave.items():
        if amplitude == 0:
            continue
        for target, factor in rotate_state(mask, nsite, c, s).items():
            result[target] += factor * amplitude
    return result


def rotation_operator(nsite, nodes):
    """K = sum_g a_g R(beta_g) restricted to M=0 on both sides (numpy)."""
    basis = balanced_basis(nsite)
    index = dict((mask, k) for k, mask in enumerate(basis))
    operator = numpy.zeros((len(basis), len(basis)), dtype=complex)
    for node in nodes:
        c, s = node_cos_sin(node)
        weight = node_weight(node)
        for col, mask in enumerate(basis):
            for target, factor in rotate_state(mask, nsite, c, s).items():
                if target in index:
                    operator[index[target], col] += weight * factor
    return basis, operator


# ---------------------------------------------------------------------------
# Quadrature nodes and the production coefficient dump
# ---------------------------------------------------------------------------
def node_weight(node):
    weight = complex(node["weight"])
    if weight.imag != 0.0:
        raise ValueError("identity-translation fixture requires real weights")
    if not math.isfinite(weight.real):
        raise ValueError("nonfinite quadrature weight")
    return weight.real


def node_cos_sin(node):
    if "c" in node and "s" in node:
        return float(node["c"]), float(node["s"])
    cc, ss = float(node["cc"]), float(node["ss"])
    if not (math.isfinite(cc) and math.isfinite(ss)) or cc < 0 or ss < 0:
        raise ValueError("invalid rotation coefficients")
    return math.sqrt(cc), math.sqrt(ss)


def independent_nodes(gauss):
    """Gauss-Legendre nodes on beta in [0, pi] from numpy, a_g included."""
    if gauss < 1:
        raise ValueError("positive Gauss-Legendre order required")
    if gauss == 1:
        return [{"spidx": 0, "c": 1.0, "s": 0.0, "cs": 0.0, "cc": 1.0,
                 "ss": 0.0, "weight": 1.0}]
    points, weights = numpy.polynomial.legendre.leggauss(gauss)
    nodes = []
    for index, (point, weight) in enumerate(zip(points, weights)):
        beta = 0.5 * math.pi * (float(point) + 1.0)
        c, s = math.cos(0.5 * beta), math.sin(0.5 * beta)
        nodes.append({
            "spidx": index, "beta": beta, "c": c, "s": s,
            "cs": c * s, "cc": c * c, "ss": s * s,
            "weight": 0.5 * math.sin(beta) * 0.5 * math.pi * float(weight),
        })
    return nodes


def fixture_nodes(betas, weights):
    nodes = []
    for index, (beta, weight) in enumerate(zip(betas, weights)):
        c, s = math.cos(0.5 * beta), math.sin(0.5 * beta)
        nodes.append({"spidx": index, "beta": beta, "c": c, "s": s,
                      "cs": c * s, "cc": c * c, "ss": s * s,
                      "weight": weight})
    return nodes


def read_qp_dump(path, allow_complex_weights=False):
    """Parse GC_SPIN_QP_V1 rows written with C %a (lossless hex floats)."""
    with open(path) as stream:
        lines = [line.split() for line in stream if line.strip()]
    if not lines:
        raise ValueError("empty spin qp dump")
    header = lines[0]
    if len(header) != 5 or header[0] != "GC_SPIN_QP_V1" or \
            header[1] != "G" or header[3] != "Q":
        raise ValueError("bad spin qp dump header: {}".format(header))
    gauss, count = int(header[2]), int(header[4])
    if gauss < 1 or count < gauss or count % gauss:
        raise ValueError("inconsistent spin qp dump sizes")
    rows = lines[1:]
    if len(rows) != count:
        raise ValueError("spin qp dump has {} rows, expected {}".format(
            len(rows), count))
    nodes = []
    for expected_q, row in enumerate(rows):
        if len(row) != 7:
            raise ValueError("spin qp dump row has {} fields".format(len(row)))
        q, spidx = int(row[0]), int(row[1])
        if q != expected_q:
            raise ValueError("spin qp dump index {} out of order".format(q))
        if spidx != q % gauss:
            raise ValueError("spin qp dump spidx {} for q {}".format(spidx, q))
        values = [float.fromhex(word) for word in row[2:]]
        if not all(math.isfinite(value) for value in values):
            raise ValueError("nonfinite spin qp dump value")
        if not allow_complex_weights and values[4] != 0.0:
            raise ValueError("spin qp dump weight has an imaginary part")
        nodes.append({"q": q, "spidx": spidx, "cs": values[0],
                      "cc": values[1], "ss": values[2],
                      "weight": complex(values[3], values[4])})
    return {"gauss": gauss, "count": count, "nodes": nodes}


def validate_coefficients(nodes, gauss, tolerance=COEFFICIENT_TOLERANCE):
    """Compare dumped coefficients with an independent leggauss(G)."""
    reference = independent_nodes(gauss)
    if len(nodes) != len(reference):
        raise AssertionError("coefficient count {} != {}".format(
            len(nodes), len(reference)))
    worst = 0.0
    for node, ref in zip(nodes, reference):
        if node["spidx"] != ref["spidx"]:
            raise AssertionError("spidx order differs")
        for key in ("cs", "cc", "ss"):
            worst = max(worst, abs(float(node[key]) - ref[key]))
        worst = max(worst, abs(complex(node["weight"]) - ref["weight"]))
        if abs(float(node["cc"]) + float(node["ss"]) - 1.0) > tolerance:
            raise AssertionError("cc+ss != 1 at spidx {}".format(
                node["spidx"]))
        if abs(float(node["cs"]) ** 2 -
               float(node["cc"]) * float(node["ss"])) > tolerance:
            raise AssertionError("cs^2 != cc*ss at spidx {}".format(
                node["spidx"]))
    if worst > tolerance:
        raise AssertionError("coefficients differ from leggauss by {}".format(
            worst))
    return worst


# ---------------------------------------------------------------------------
# Projected wave functions
# ---------------------------------------------------------------------------
def rotated_block_matrix(F, node, convert=rational, zero=ZERO,
                         mul=qmul, sub=qsub, add=qadd, neg=None):
    """Production 4-block matrix with the node's coefficients (Q+iQ)."""
    rows = _as_rows(F)
    nsite = len(rows)
    cs, cc, ss = (convert(node[key]) for key in ("cs", "cc", "ss"))
    if neg is None:
        def neg(value):
            return (-value[0], -value[1])
    matrix = [[zero] * (2 * nsite) for unused in range(2 * nsite)]
    for i in range(nsite):
        for j in range(nsite):
            fij, fji = convert(rows[i][j]), convert(rows[j][i])
            antisym = sub(fij, fji)
            matrix[i][j] = neg(mul(antisym, cs))
            matrix[i][j + nsite] = add(mul(fij, cc), mul(fji, ss))
            matrix[i + nsite][j] = neg(add(mul(fij, ss), mul(fji, cc)))
            matrix[i + nsite][j + nsite] = mul(antisym, cs)
    return matrix


def component_reference(F, node, mask):
    """Exact Pf of one rotated component on mask (no weight, no charge)."""
    values = rational_pfaffians(rotated_block_matrix(F, node))
    if popcount(mask) % 2:
        return ZERO
    return values[mask]


def _rational_finite(F, nodes):
    rows = _as_rows(F)
    nsite = len(rows)
    total = dict((mask, ZERO) for mask in balanced_basis(nsite))
    for node in nodes:
        weight = rational(node_weight(node))
        values = rational_pfaffians(rotated_block_matrix(rows, node))
        for mask in total:
            total[mask] = qadd(total[mask], qmul(weight, values[mask]))
    return total


def _finish(nsite, amplitudes, log_charge):
    result = dict((mask, 0.0j) for mask in even_masks(nsite))
    for mask, value in amplitudes.items():
        if value == ZERO:
            continue
        logc = log_charge(mask)
        if not math.isfinite(logc):
            raise ValueError("nonfinite log charge factor")
        amplitude = qfloat(value) * math.exp(logc)
        if amplitude == 0 or not (math.isfinite(amplitude.real) and
                                  math.isfinite(amplitude.imag)):
            raise ValueError("nonzero amplitude is not representable")
        result[mask] = amplitude
    return result


def finite_fock(F, nodes):
    """Floating Fock-rotation sum restricted to M=0 (no charge factor)."""
    rows = _as_rows(F)
    nsite = len(rows)
    base = pair_wave(rows)
    total = dict((mask, 0.0j) for mask in balanced_basis(nsite))
    for node in nodes:
        c, s = node_cos_sin(node)
        weight = node_weight(node)
        rotated = rotate_fock(base, nsite, c, s)
        for mask in total:
            total[mask] += weight * rotated[mask]
    return total


def project_finite(F, nodes, log_charge, reference_nodes=None):
    """Finite-G wave from the rational reference of the given coefficients.

    The rational 4-block sum is the main value.  It must agree with the
    independent floating Fock rotation evaluated at reference_nodes (the
    same nodes when omitted) within FOCK_TOLERANCE in raw amplitude.
    """
    rows = _as_rows(F)
    nsite = len(rows)
    exact = _rational_finite(rows, nodes)
    fock = finite_fock(rows, reference_nodes or nodes)
    for mask, value in exact.items():
        reference = qfloat(value)
        if abs(reference - fock[mask]) > FOCK_TOLERANCE * (1.0 + abs(reference)):
            raise AssertionError(
                "rational and Fock rotation differ at mask {}: {} vs {}".format(
                    mask, reference, fock[mask]))
    return _finish(nsite, exact, log_charge)


def rational_finite(F, nodes):
    """Q+iQ finite-G amplitudes on M=0 masks (for exact zero classification)."""
    return _rational_finite(F, nodes)


def _sector_projector_apply(vector, matrix, smax):
    """prod_{s=1..smax} (1 - S^2/(s(s+1))) applied to a Q+iQ vector."""
    current = list(vector)
    size = len(current)
    for spin in range(1, smax + 1):
        scale = Fraction(1, spin * (spin + 1))
        updated = []
        for row in range(size):
            acc = ZERO
            for col in range(size):
                entry = int(matrix[row][col])
                if entry:
                    acc = qadd(acc, qscale(current[col], entry))
            updated.append(qsub(current[row], qscale(acc, scale)))
        current = updated
    return current


def rational_exact(F):
    """Exact S=0 projection P0|phi_F> on M=0 masks over Q+iQ."""
    rows = _as_rows(F)
    nsite = len(rows)
    pair = _rational_pair_wave(rows)
    result = {}
    for count in range(0, 2 * nsite + 1, 2):
        sector = tuple(mask for mask in balanced_basis(nsite)
                       if popcount(mask) == count)
        if not sector:
            continue
        smax = min(count, 2 * nsite - count) // 2
        matrix = spin2_matrix(nsite, sector)
        projected = _sector_projector_apply(
            [pair[mask] for mask in sector], matrix, smax)
        for mask, value in zip(sector, projected):
            result[mask] = value
    return result


def project_exact(F, log_charge):
    rows = _as_rows(F)
    return _finish(len(rows), rational_exact(rows), log_charge)


# ---------------------------------------------------------------------------
# Derivatives, support graph and observables
# ---------------------------------------------------------------------------
def _finite_difference(evaluate, F, i, j, imaginary, epsilon):
    step = 1j * epsilon if imaginary else epsilon
    plus = [list(row) for row in F]
    minus = [list(row) for row in F]
    plus[i][j] += step
    minus[i][j] -= step
    first, second = evaluate(plus), evaluate(minus)
    return dict((mask, (first[mask] - second[mask]) / (2.0 * epsilon))
                for mask in first)


def log_derivatives(F, nodes, masks, log_charge, epsilons=(1.0e-5, 5.0e-6)):
    """d log Psi / d Re F_ij and d Im F_ij on masks via the Fock rotation.

    The charge factor does not depend on F and cancels.  Two central
    differences are combined by Richardson extrapolation and must agree.
    Returns {mask: array [dReF00, dImF00, dReF01, ...]}.
    """
    rows = _as_rows(F)
    nsite = len(rows)
    basis, operator = rotation_operator(nsite, nodes)
    index = dict((mask, k) for k, mask in enumerate(basis))

    def evaluate(matrix):
        pair = pair_wave(matrix)
        vector = numpy.array([pair[mask] for mask in basis])
        projected = operator @ vector
        return dict((mask, projected[index[mask]]) for mask in masks)

    base = evaluate(rows)
    for mask in masks:
        if base[mask] == 0:
            raise ValueError("log derivative requested on a zero amplitude")
    result = dict((mask, numpy.zeros(2 * nsite * nsite, dtype=complex))
                  for mask in masks)
    for i in range(nsite):
        for j in range(nsite):
            for imaginary in (0, 1):
                coarse = _finite_difference(evaluate, rows, i, j, imaginary,
                                            epsilons[0])
                fine = _finite_difference(evaluate, rows, i, j, imaginary,
                                          epsilons[1])
                slot = 2 * (i * nsite + j) + imaginary
                for mask in masks:
                    if abs(coarse[mask] - fine[mask]) > \
                            1.0e-6 * (1.0 + abs(fine[mask])):
                        raise AssertionError("unstable finite difference")
                    ratio = (epsilons[0] / epsilons[1]) ** 2
                    value = (ratio * fine[mask] - coarse[mask]) / (ratio - 1.0)
                    result[mask][slot] = value / base[mask]
    return result


def support_reachable(support, nsite, start):
    """BFS over same-spin hops and up-down pair add/remove inside support."""
    support = set(support)
    if start not in support:
        raise AssertionError("BFS start {} is outside the support".format(start))
    seen = set([start])
    queue = deque([start])
    while queue:
        mask = queue.popleft()
        neighbours = []
        for spin in (0, 1):
            for i in range(nsite):
                for j in range(nsite):
                    if i == j:
                        continue
                    moved = apply_ops(mask, one_body(j + spin * nsite,
                                                     i + spin * nsite))
                    if moved is not None:
                        neighbours.append(moved[0])
        for i in range(nsite):
            for j in range(nsite):
                up, down = 1 << i, 1 << (j + nsite)
                if not mask & up and not mask & down:
                    neighbours.append(mask | up | down)
                if mask & up and mask & down:
                    neighbours.append(mask & ~up & ~down)
        for target in neighbours:
            if target in support and target not in seen:
                seen.add(target)
                queue.append(target)
    missing = support - seen
    if missing:
        raise AssertionError("support masks unreachable from {}: {}".format(
            start, sorted(missing)))
    return len(seen)


def norm(wave):
    return sum(abs(value) ** 2 for value in wave.values())


def sector_norms(wave):
    sectors = {}
    for mask, value in wave.items():
        count = popcount(mask)
        sectors[count] = sectors.get(count, 0.0) + abs(value) ** 2
    return sectors


def spin2_expectation(wave, nsite):
    basis = balanced_basis(nsite)
    vector = numpy.array([wave.get(mask, 0.0j) for mask in basis])
    matrix = spin2_matrix(nsite, basis)
    return (numpy.vdot(vector, matrix @ vector) /
            numpy.vdot(vector, vector)).real


def wave_error(finite, exact, basis):
    """Relative error after removing the common complex normalization."""
    f = numpy.array([finite.get(mask, 0.0j) for mask in basis])
    e = numpy.array([exact.get(mask, 0.0j) for mask in basis])
    alpha = numpy.vdot(e, f) / numpy.vdot(e, e)
    return float(numpy.linalg.norm(f - alpha * e) / numpy.linalg.norm(f)), alpha


def triplet_probe(gauss):
    """Pure-triplet weights sum_g a_g d^1_00(beta_g) three ways."""
    nodes = independent_nodes(gauss)
    by_cos = 0.0
    by_cc_ss = 0.0
    exact = Fraction(0)
    for node in nodes:
        by_cos += node["weight"] * math.cos(node["beta"])
        by_cc_ss += node["weight"] * (node["cc"] - node["ss"])
        exact += Fraction(node["weight"]) * (Fraction(node["cc"]) -
                                             Fraction(node["ss"]))
    return {"gauss": gauss, "sum_cos": by_cos, "sum_cc_minus_ss": by_cc_ss,
            "rational_cc_minus_ss": float(exact)}


# ---------------------------------------------------------------------------
# Self test
# ---------------------------------------------------------------------------
def _check(condition, message):
    if not condition:
        raise AssertionError(message)


def _fixture_F(nsite):
    return numpy.array([[complex((11 + 3 * i - 2 * j + 17 * (i == j)) / 32.0,
                                 (1 + 2 * i + j) / 64.0)
                         for j in range(nsite)] for i in range(nsite)])


def _probe_F(nsite):
    """Design-probe state whose antisymmetric part has rank > 2, so its
    N=4 sector carries S=2 weight (the fixture F above has none)."""
    F = numpy.array([[complex((7 + 3 * i - 2 * j) / 16.0,
                              (1 + i + 2 * j) / 32.0)
                      for j in range(nsite)] for i in range(nsite)])
    F[0, 1] += 0.7
    F[2, 0] -= 0.3j
    F[3, 3] += 1.0 / 64.0
    return F


def _no_charge(unused_mask):
    return 0.0


def _expect_failure(label, function):
    try:
        function()
    except (AssertionError, ValueError):
        return
    raise AssertionError("mutation/negative '{}' was not detected".format(label))


def self_test():
    L = 4
    basis = balanced_basis(L)
    _check(len(basis) == 70, "M=0 basis size")
    F = _fixture_F(L)
    exact = project_exact(F, _no_charge)
    s2 = spin2_matrix(L, basis)
    v = numpy.array([exact[x] for x in basis])
    _check(numpy.max(numpy.abs(s2 @ v)) < 1.e-11, "S^2 P0 != 0")
    _check(exact[0] == 1.0, "vacuum amplitude of P0")

    # Rational S^2 P0 = 0 and P0^2 = P0 in every sector.
    rational = rational_exact(F)
    for count in range(0, 2 * L + 1, 2):
        sector = tuple(mask for mask in basis if popcount(mask) == count)
        smax = min(count, 2 * L - count) // 2
        matrix = spin2_matrix(L, sector)
        vector = [rational[mask] for mask in sector]
        again = _sector_projector_apply(vector, matrix, smax)
        _check(again == vector, "P0^2 != P0 in N={}".format(count))
        applied = _sector_projector_apply(vector, matrix, 0)
        for row in range(len(sector)):
            acc = ZERO
            for col in range(len(sector)):
                if matrix[row][col]:
                    acc = qadd(acc, qscale(vector[col], int(matrix[row][col])))
            _check(acc == ZERO, "rational S^2 P0 != 0 in N={}".format(count))
        _check(applied == vector, "empty projector product is not identity")
    full = (1 << (2 * L)) - 1
    _check(abs(exact[full]) > 0, "full filling is a singlet")

    # N=2: a symmetric F is a pure singlet and an antisymmetric F a pure
    # triplet whose exact S=0 component vanishes identically.
    symmetric = numpy.array([[0.3, 0.7 - 0.2j], [0.7 - 0.2j, -0.4]])
    antisymmetric = numpy.array([[0.0, 0.6 + 0.1j], [-0.6 - 0.1j, 0.0]])
    pair_sym = pair_wave(symmetric)
    exact_sym = project_exact(symmetric, _no_charge)
    for mask in balanced_basis(2):
        if popcount(mask) == 2:
            _check(abs(exact_sym[mask] - pair_sym[mask]) < 1e-15,
                   "symmetric F is not a pure singlet")
    triplet_exact = rational_exact(antisymmetric)
    for mask, value in triplet_exact.items():
        if popcount(mask) == 2:
            _check(value == ZERO, "triplet has an S=0 component")

    # A non-scalar SpinJastrow factor breaks S=0 (design probe state, where
    # exp(0.6 m0 m1) P0|phi> has S^2 ~ 1.0e-2).
    def spin_jastrow(mask):
        up0, up1 = mask & 1, (mask >> 1) & 1
        down0, down1 = (mask >> L) & 1, (mask >> (L + 1)) & 1
        return 0.6 * (up0 - down0) * (up1 - down1)
    probe_F = _probe_F(L)
    _check(spin2_expectation(project_exact(probe_F, _no_charge), L) < 1e-14,
           "probe state is not a singlet")
    with_spin_jastrow = project_exact(probe_F, spin_jastrow)
    _check(spin2_expectation(with_spin_jastrow, L) > 1e-4,
           "SpinJastrow should break the singlet")

    # Fock rotation versus the rational 4-block reference, several G and L.
    for nsite in (2, 4):
        Fn = _fixture_F(nsite)
        exact_n = project_exact(Fn, _no_charge)
        basis_n = balanced_basis(nsite)
        for gauss in (1, 2, 3, 4, 8, 16, 32):
            nodes = independent_nodes(gauss)
            finite = project_finite(Fn, nodes, _no_charge)
            vacuum = sum(node["weight"] for node in nodes)
            _check(abs(finite[0] - vacuum) < 1e-15,
                   "raw vacuum amplitude is not sum a_g")
            for mask in finite:
                if not is_balanced(mask, nsite):
                    _check(finite[mask] == 0, "M != 0 amplitude survived")
            if gauss >= 16:
                error, alpha = wave_error(finite, exact_n, basis_n)
                _check(error <= 1e-11, "G={} L={} wave error {}".format(
                    gauss, nsite, error))
                _check(abs(alpha - 1.0) < 1e-11,
                       "finite G normalization differs from P0")
                _check(spin2_expectation(finite, nsite) <= 1e-11,
                       "G={} S^2 residual".format(gauss))
                finite_sectors = sector_norms(finite)
                exact_sectors = sector_norms(exact_n)
                for count, value in exact_sectors.items():
                    _check(abs(finite_sectors[count] - value) <=
                           1e-11 * (1.0 + value),
                           "relative sector weight N={}".format(count))
    # Mutations: weight without 1/2 (raw vacuum) and without sin(beta).
    nodes = independent_nodes(8)
    doubled = [dict(node, weight=2.0 * node["weight"]) for node in nodes]
    reference = project_finite(F, nodes, _no_charge)
    _expect_failure("weight 1/2 dropped", lambda: _check(
        abs(project_finite(F, doubled, _no_charge)[0] - reference[0]) < 1e-12,
        "raw vacuum"))
    no_sine = [dict(node, weight=node["weight"] / math.sin(node["beta"]))
               for node in nodes]
    _check(abs(project_finite(F, no_sine, _no_charge)[0] - 1.0) > 0.1,
           "missing sin(beta) is not detected in the raw vacuum")
    probe_F = _probe_F(L)
    mutated = project_finite(probe_F, no_sine, _no_charge)
    error, unused_alpha = wave_error(
        mutated, project_exact(probe_F, _no_charge), basis)
    _check(error > 1e-3 and spin2_expectation(mutated, L) > 1e-5,
           "missing sin(beta) is not detected in a nontrivial spin state")
    _expect_failure("nan F", lambda: pair_wave(
        [[float("nan"), 0.0], [0.0, 0.0]]))
    _expect_failure("inf weight", lambda: project_finite(
        F, [dict(nodes[0], weight=float("inf"))], _no_charge))

    # Derivatives agree with a direct central difference of the rational
    # reference on two masks.
    nodes = fixture_nodes((0.37, 1.21), (0.7, 0.3))
    masks = [m for m in basis if popcount(m) == 4][:2]
    derivatives = log_derivatives(F, nodes, masks, _no_charge)
    for mask in masks:
        base = qfloat(_rational_finite(F, nodes)[mask])
        for slot in (0, 1, 7):
            i, j, imaginary = slot // 2 // L, slot // 2 % L, slot % 2
            step = 1e-6 * (1j if imaginary else 1.0)
            plus = [list(row) for row in F]
            minus = [list(row) for row in F]
            plus[i][j] += step
            minus[i][j] -= step
            direct = (qfloat(_rational_finite(plus, nodes)[mask]) -
                      qfloat(_rational_finite(minus, nodes)[mask])) / 2e-6
            _check(abs(direct / base - derivatives[mask][slot]) <
                   1e-7 * (1.0 + abs(derivatives[mask][slot])),
                   "log derivative mismatch")

    # Support graph of a finite-G reference is connected from the vacuum.
    finite = rational_finite(F, independent_nodes(4))
    support = [mask for mask, value in finite.items() if value != ZERO]
    support_reachable(support, L, 0)

    # Pure-triplet probe: tiny but not necessarily zero in binary64.
    for gauss in (4, 8, 16, 32):
        probe = triplet_probe(gauss)
        for key in ("sum_cos", "sum_cc_minus_ss", "rational_cc_minus_ss"):
            _check(abs(probe[key]) <= COEFFICIENT_TOLERANCE,
                   "triplet probe {} G={}".format(key, gauss))
    print("grandcanonical_spin_oracle self-test passed")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if not args.self_test:
        parser.error("only --self-test is available")
    self_test()


if __name__ == "__main__":
    sys.exit(main())
