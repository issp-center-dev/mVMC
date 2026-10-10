"""Independent exact oracle for grand-canonical anti-periodic fixtures.

The lattice geometry, not the input files, defines the expected orbital and
Hamiltonian tables: every pair orbital F_IJ = phi_{s_I s_J}(r_J - r_I) of a
translation-invariant state obeys phi(r + L e_k) = theta_k phi(r) with
theta_k = -1 for an anti-periodic direction.  The independent file parsers
below read what a fixture actually wrote, so a fixture or production table
that drifts from the geometry is detected instead of being re-read as truth.

Only the site-count independent Pfaffian and Fock-operator helpers are shared
with the two-site oracle in grandcanonical_exact_oracle.py.
"""
from __future__ import print_function

import math
import random
import sys

from grandcanonical_exact_oracle import (
    adjoint,
    apply_ops,
    nbody,
    one_body,
    pair_create,
    pfaffian,
)


AP = "AP"
P = "P"


class Lattice(object):
    """Hypercubic lattice; site = x + Lx*y + ..., orbital = site + spin*N."""

    def __init__(self, lengths, boundaries):
        if len(lengths) != len(boundaries):
            raise ValueError("one boundary per direction is required")
        self.lengths = tuple(int(length) for length in lengths)
        self.boundaries = tuple(boundaries)
        for boundary in self.boundaries:
            if boundary not in (AP, P):
                raise ValueError("boundary must be AP or P")
        self.theta = tuple(-1 if boundary == AP else 1
                           for boundary in self.boundaries)
        self.nsite = 1
        for length in self.lengths:
            self.nsite *= length
        self.norbital = 2 * self.nsite

    def coordinates(self, site):
        coords = []
        rest = site
        for length in self.lengths:
            coords.append(rest % length)
            rest //= length
        return tuple(coords)

    def site(self, coords):
        site = 0
        stride = 1
        for coordinate, length in zip(coords, self.lengths):
            site += (coordinate % length) * stride
            stride *= length
        return site

    def displacement(self, first, second):
        """Canonical displacement d in [0,L) and wrap count of r_j - r_i."""
        canonical = []
        wraps = []
        for a, b, length in zip(self.coordinates(first),
                                self.coordinates(second), self.lengths):
            r = b - a
            d = r % length
            canonical.append(d)
            wraps.append((r - d) // length)
        return tuple(canonical), tuple(wraps)

    def wrap_phase(self, wraps):
        phase = 1
        for wrap, theta in zip(wraps, self.theta):
            if wrap % 2:
                phase *= theta
        return phase

    def neighbor(self, site, axis):
        """Site reached by +e_axis and the boundary phase picked up."""
        coords = list(self.coordinates(site))
        coords[axis] += 1
        wrapped = coords[axis] == self.lengths[axis]
        return self.site(coords), (self.theta[axis] if wrapped else 1)

    def translate_orbital(self, orbital, axis):
        site, spin = orbital % self.nsite, orbital // self.nsite
        target, phase = self.neighbor(site, axis)
        return target + spin * self.nsite, phase

    def label(self):
        return "x".join(str(length) for length in self.lengths) + "_" + \
            "_".join(self.boundaries)


# ---------------------------------------------------------------------------
# Geometry-defined expected tables
# ---------------------------------------------------------------------------

def orbital_key_and_sign(lattice, first, second):
    """Shared parameter key and sign for the upper-triangle entry I<J.

    Opposite spins use phi_ud(r) for all r.  For equal spins antisymmetry
    gives phi(-r) = -phi(r), hence phi(d') = c phi(d) with d' = -d mod L and
    c = -prod_{k: d_k != 0} theta_k; a self-conjugate d with c = -1 is
    forced to vanish.
    """
    nsite = lattice.nsite
    site_i, spin_i = first % nsite, first // nsite
    site_j, spin_j = second % nsite, second // nsite
    canonical, wraps = lattice.displacement(site_i, site_j)
    sign = lattice.wrap_phase(wraps)
    if spin_i != spin_j:
        return ("ud", canonical), sign
    conjugate = tuple((-component) % length
                      for component, length in zip(canonical, lattice.lengths))
    factor = -1
    for component, theta in zip(canonical, lattice.theta):
        if component != 0:
            factor *= theta
    if canonical == conjugate:
        if factor == -1:
            return ("zero",), 1
        return ("ss", spin_i, canonical), sign
    if canonical < conjugate:
        return ("ss", spin_i, canonical), sign
    return ("ss", spin_i, conjugate), sign * factor


def key_order(key):
    if key[0] == "ud":
        return (0, key[1])
    if key[0] == "ss":
        return (1 + key[1], key[2])
    return (3, ())


def orbital_geometry(lattice):
    """Return ({(I,J): (index, sign)} for I<J, ordered parameter keys)."""
    keyed = {}
    for first in range(lattice.norbital):
        for second in range(first + 1, lattice.norbital):
            keyed[(first, second)] = orbital_key_and_sign(lattice, first, second)
    keys = sorted(set(key for key, unused in keyed.values()), key=key_order)
    index = dict((key, number) for number, key in enumerate(keys))
    table = dict((pair, (index[key], sign))
                 for pair, (key, sign) in keyed.items())
    return table, keys


def expanded_table(table):
    """One independent parameter per pair, all signs +1 (parameter order)."""
    pairs = sorted(table)
    return dict((pair, (number, 1)) for number, pair in enumerate(pairs)), pairs


def expanded_parameters(table, parameters):
    unused, pairs = expanded_table(table)
    return tuple(table[pair][1] * parameters[table[pair][0]] for pair in pairs)


def pair_wrap_sign(lattice, first, second):
    """Boundary phase of r_J - r_I alone (one parameter per pair, no sharing)."""
    nsite = lattice.nsite
    unused, wraps = lattice.displacement(first % nsite, second % nsite)
    return lattice.wrap_phase(wraps)


def negative_sign_count(table):
    return sum(1 for unused, sign in table.values() if sign < 0)


def full_sign_tables(table, norbital):
    """Expand I<J rows with the reader's antisymmetric completion."""
    index = {}
    sign = {}
    for (first, second), (parameter, value) in table.items():
        index[(first, second)] = parameter
        index[(second, first)] = parameter
        sign[(first, second)] = value
        sign[(second, first)] = -value
    if len(index) != norbital * (norbital - 1):
        raise AssertionError("pair table does not cover all I != J")
    return index, sign


def pair_matrix(table, parameters, norbital):
    index, sign = full_sign_tables(table, norbital)
    matrix = [[0.0j for unused in range(norbital)] for unused in range(norbital)]
    for first in range(norbital):
        for second in range(norbital):
            if first == second:
                continue
            # F_IJ = sign_IJ f[index_IJ] - sign_JI f[index_JI]
            matrix[first][second] = (
                sign[(first, second)] * parameters[index[(first, second)]] -
                sign[(second, first)] * parameters[index[(second, first)]])
    return matrix


def translation_defect(lattice, matrix):
    """max |F_{TI,TJ} - s_I s_J F_IJ| over unit translations T."""
    worst = 0.0
    for axis in range(len(lattice.lengths)):
        if lattice.lengths[axis] == 1:
            continue
        for first in range(lattice.norbital):
            target_i, phase_i = lattice.translate_orbital(first, axis)
            for second in range(lattice.norbital):
                if first == second:
                    continue
                target_j, phase_j = lattice.translate_orbital(second, axis)
                worst = max(worst, abs(
                    matrix[target_i][target_j] -
                    phase_i * phase_j * matrix[first][second]))
    return worst


def hopping_coefficients(lattice, hopping):
    """{(i,j): t_ij} with H_hop = -sum_{ij,s} t_ij c+_is c_js.

    Each +e_axis bond contributes t*phase to (i,j) and conj(t)*phase to
    (j,i); a length-2 periodic direction therefore keeps both bonds.
    """
    coefficients = {}
    for site in range(lattice.nsite):
        for axis, length in enumerate(lattice.lengths):
            if length == 1:
                continue
            target, phase = lattice.neighbor(site, axis)
            for key, value in (((site, target), hopping * phase),
                               ((target, site), hopping.conjugate() * phase)):
                coefficients[key] = coefficients.get(key, 0.0j) + value
    return coefficients


def bond_multiplicity(lattice):
    """{(i,j) with i<j: number of +e_axis bonds joining i and j}."""
    bonds = {}
    for site in range(lattice.nsite):
        for axis, length in enumerate(lattice.lengths):
            if length == 1:
                continue
            target, unused = lattice.neighbor(site, axis)
            key = (min(site, target), max(site, target))
            bonds[key] = bonds.get(key, 0) + 1
    return bonds


def hopping_translation_defect(lattice, coefficients):
    worst = 0.0
    for axis in range(len(lattice.lengths)):
        if lattice.lengths[axis] == 1:
            continue
        for (first, second), value in coefficients.items():
            target_i, phase_i = lattice.neighbor(first, axis)
            target_j, phase_j = lattice.neighbor(second, axis)
            worst = max(worst, abs(
                coefficients.get((target_i, target_j), 0.0j) -
                phase_i * phase_j * value))
    return worst


def jastrow_classes(lattice):
    """{(i,j): class} for i != j by minimum-image distance per direction."""
    raw = {}
    for first in range(lattice.nsite):
        for second in range(lattice.nsite):
            if first == second:
                continue
            canonical, unused = lattice.displacement(first, second)
            raw[(first, second)] = tuple(
                min(component, length - component)
                for component, length in zip(canonical, lattice.lengths))
    keys = sorted(set(raw.values()))
    index = dict((key, number) for number, key in enumerate(keys))
    return dict((pair, index[key]) for pair, key in raw.items()), len(keys)


# ---------------------------------------------------------------------------
# Model and exact expectation values
# ---------------------------------------------------------------------------

class Model(object):
    def __init__(self, lattice, table, parameters, gutzwiller, jastrow,
                 hopping, mu, coulomb_intra, coulomb_inter, delta=None,
                 anomalous_pair=None, anomalous_green_pairs=()):
        self.lattice = lattice
        self.table = table
        self.parameters = tuple(parameters)
        self.gutzwiller = gutzwiller
        self.jastrow = tuple(jastrow)
        self.hopping = hopping
        self.mu = mu
        self.coulomb_intra = coulomb_intra
        self.coulomb_inter = coulomb_inter
        self.delta = delta
        self.anomalous_pair = anomalous_pair
        self.anomalous_green_pairs = tuple(anomalous_green_pairs)
        self.jastrow_class, count = jastrow_classes(lattice)
        if count != len(self.jastrow):
            raise ValueError("Jastrow class count mismatch")

    def replace(self, **changes):
        values = dict(self.__dict__)
        values.pop("jastrow_class")
        values.update(changes)
        return Model(**values)

    def proj_values(self):
        return (self.gutzwiller,) + self.jastrow


def even_basis(norbital):
    return tuple(state for state in range(1 << norbital)
                 if bin(state).count("1") % 2 == 0)


def log_projection(model, state):
    nsite = model.lattice.nsite
    up = [(state >> site) & 1 for site in range(nsite)]
    down = [(state >> (site + nsite)) & 1 for site in range(nsite)]
    value = model.gutzwiller * sum(a * b for a, b in zip(up, down))
    for first in range(nsite):
        x_first = up[first] + down[first] - 1
        for second in range(first + 1, nsite):
            x_second = up[second] + down[second] - 1
            value += (model.jastrow[model.jastrow_class[(first, second)]] *
                      x_first * x_second)
    return value


def wavefunction(model, parameters=None, table=None):
    if parameters is None:
        parameters = model.parameters
    if table is None:
        table = model.table
    norbital = model.lattice.norbital
    matrix = pair_matrix(table, parameters, norbital)
    result = {}
    for state in even_basis(norbital):
        orbitals = [orbital for orbital in range(norbital)
                    if (state >> orbital) & 1]
        restricted = [[matrix[row][column] for column in orbitals]
                      for row in orbitals]
        result[state] = pfaffian(restricted) * math.exp(
            log_projection(model, state))
    return result


def model_terms(model):
    lattice = model.lattice
    nsite = lattice.nsite
    terms = []
    for orbital in range(lattice.norbital):
        terms.append(("trans_mu", -model.mu, one_body(orbital, orbital)))
    for (first, second), value in sorted(
            hopping_coefficients(lattice, model.hopping).items()):
        for spin in (0, 1):
            terms.append(("trans", -value,
                          one_body(first + spin * nsite,
                                   second + spin * nsite)))
    for site in range(nsite):
        terms.append(("coulomb_intra", model.coulomb_intra,
                      nbody(((site, site), (site + nsite, site + nsite)))))
    for (first, second), count in sorted(bond_multiplicity(lattice).items()):
        for a in (first, first + nsite):
            for b in (second, second + nsite):
                terms.append(("coulomb_inter", count * model.coulomb_inter,
                              nbody(((a, a), (b, b)))))
    if model.delta is not None:
        operators = pair_create(*model.anomalous_pair)
        terms.append(("anomalous", model.delta, operators))
        terms.append(("anomalous", model.delta.conjugate(), adjoint(operators)))
    return tuple(terms)


def anomalous_green_operators(model):
    """Rows (type, orbital1, orbital2, operators) as in AnomalousG files."""
    rows = []
    for first, second in model.anomalous_green_pairs:
        rows.append((1, first, second, (("c", first), ("c", second))))
        rows.append((0, second, first, (("a", second), ("a", first))))
        rows.append((1, second, first, (("c", second), ("c", first))))
        rows.append((0, first, second, (("a", first), ("a", second))))
    return tuple(rows)


def norm(wave):
    return sum(abs(amplitude) ** 2 for amplitude in wave.values())


def expectation(wave, operators, denominator=None):
    if denominator is None:
        denominator = norm(wave)
    numerator = 0.0j
    for state, amplitude in wave.items():
        transformed = apply_ops(state, operators)
        if transformed is None:
            continue
        target, sign = transformed
        numerator += wave.get(target, 0.0j).conjugate() * sign * amplitude
    return numerator / denominator


def hamiltonian_action(wave, terms):
    result = dict((state, 0.0j) for state in wave)
    for unused_label, coefficient, operators in terms:
        for state, amplitude in wave.items():
            transformed = apply_ops(state, operators)
            if transformed is None:
                continue
            target, sign = transformed
            if target in result:
                result[target] += coefficient * sign * amplitude
    return result


def exact_observables(model, parameters=None, table=None, terms=None):
    wave = wavefunction(model, parameters, table)
    denominator = norm(wave)
    if terms is None:
        terms = model_terms(model)
    action = hamiltonian_action(wave, terms)
    energy = sum(wave[state].conjugate() * action[state]
                 for state in wave) / denominator
    energy2 = sum(abs(action[state]) ** 2 for state in wave) / denominator
    number = sum(abs(amplitude) ** 2 * bin(state).count("1")
                 for state, amplitude in wave.items()) / denominator
    number2 = sum(abs(amplitude) ** 2 * bin(state).count("1") ** 2
                  for state, amplitude in wave.items()) / denominator
    norbital = model.lattice.norbital
    greens1 = dict(
        ((first, second), expectation(wave, one_body(first, second),
                                      denominator))
        for first in range(norbital) for second in range(norbital))
    anomalous_g = dict(
        ((kind, first, second), expectation(wave, operators, denominator))
        for kind, first, second, operators in anomalous_green_operators(model))
    sectors = {}
    for state, amplitude in wave.items():
        count = bin(state).count("1")
        sectors[count] = sectors.get(count, 0.0) + abs(amplitude) ** 2 / denominator
    return {
        "wave": wave,
        "energy": energy,
        "energy2": energy2,
        "number": number,
        "number2": number2,
        "variance_number": number2 - number * number,
        "greens1": greens1,
        "anomalous_g": anomalous_g,
        "sectors": sectors,
    }


def energy(model, parameters):
    return exact_observables(model, parameters)["energy"].real


def parameter_gradient(model, parameter_index, imaginary):
    """Exact SR force 2 Re(<O* H> - <O>* <H>) for one orbital parameter."""
    parameters = list(model.parameters)
    wave = wavefunction(model, parameters)
    action = hamiltonian_action(wave, model_terms(model))
    denominator = norm(wave)
    mean_h = sum(wave[state].conjugate() * action[state]
                 for state in wave) / denominator
    epsilon = 1.0e-7
    step = 1j * epsilon if imaginary else epsilon
    plus = list(parameters)
    minus = list(parameters)
    plus[parameter_index] += step
    minus[parameter_index] -= step
    wave_plus = wavefunction(model, plus)
    wave_minus = wavefunction(model, minus)
    mean_o = 0.0j
    mean_oh = 0.0j
    for state, amplitude in wave.items():
        if amplitude == 0.0:
            continue
        observable = (wave_plus[state] - wave_minus[state]) / (
            2.0 * epsilon * amplitude)
        probability = abs(amplitude) ** 2 / denominator
        local_energy = action[state] / amplitude
        mean_o += probability * observable
        mean_oh += probability * observable.conjugate() * local_energy
    return 2.0 * (mean_oh - mean_o.conjugate() * mean_h).real


def finite_difference(model, parameter_index, imaginary, epsilon):
    plus = list(model.parameters)
    minus = list(model.parameters)
    step = 1j * epsilon if imaginary else epsilon
    plus[parameter_index] += step
    minus[parameter_index] -= step
    return (energy(model, plus) - energy(model, minus)) / (2.0 * epsilon)


# ---------------------------------------------------------------------------
# Independent input-file parsers
# ---------------------------------------------------------------------------

def _data_lines(path, skip=5):
    with open(path) as stream:
        lines = stream.readlines()
    return [line for line in lines[skip:] if line.strip()]


def _header_count(path, keyword):
    with open(path) as stream:
        lines = stream.readlines()
    words = lines[1].split()
    if len(words) != 2 or words[0] != keyword:
        raise AssertionError("{}: header {} expected".format(path, keyword))
    return int(words[1])


def parse_orbital_general(path, nsite):
    """Return ({(I,J): (index, sign)}, nparameter) from orbitalidxgen.def."""
    nparameter = _header_count(path, "NOrbitalIdx")
    norbital = 2 * nsite
    expected = norbital * (norbital - 1) // 2
    lines = _data_lines(path)
    if len(lines) != expected + nparameter:
        raise AssertionError("{}: unexpected row count".format(path))
    table = {}
    for line in lines[:expected]:
        words = line.split()
        if len(words) != 6:
            raise AssertionError("{}: pair row is not 6 columns".format(path))
        site_i, spin_i, site_j, spin_j, index, sign = (int(w) for w in words)
        first = site_i + spin_i * nsite
        second = site_j + spin_j * nsite
        if not first < second or (first, second) in table:
            raise AssertionError("{}: bad pair {}".format(path, words))
        if sign not in (1, -1) or not 0 <= index < nparameter:
            raise AssertionError("{}: bad index/sign {}".format(path, words))
        table[(first, second)] = (index, sign)
    flags = {}
    for line in lines[expected:]:
        words = line.split()
        flags[int(words[0])] = int(words[1])
    if sorted(flags) != list(range(nparameter)):
        raise AssertionError("{}: OptFlag rows incomplete".format(path))
    return table, nparameter


def parse_trans_sym(path, nsite):
    count = _header_count(path, "NQPTrans")
    lines = _data_lines(path)
    weights = []
    for line in lines[:count]:
        words = line.split()
        # StdFace writes the weight without an imaginary column.
        imaginary = float(words[2]) if len(words) > 2 else 0.0
        weights.append(complex(float(words[1]), imaginary))
    mappings = {}
    for line in lines[count:]:
        words = [int(word) for word in line.split()]
        if len(words) != 4:
            raise AssertionError("{}: TransSym row is not 4 columns".format(path))
        mappings[(words[0], words[1])] = (words[2], words[3])
    if len(mappings) != count * nsite:
        raise AssertionError("{}: TransSym rows incomplete".format(path))
    return weights, mappings


def parse_trans(path):
    """{((i,s),(j,t)): summed coefficient} from a Trans file."""
    coefficients = {}
    for line in _data_lines(path):
        words = line.split()
        key = ((int(words[0]), int(words[1])), (int(words[2]), int(words[3])))
        value = complex(float(words[4]), float(words[5]))
        coefficients[key] = coefficients.get(key, 0.0j) + value
    return coefficients


def parse_initial(path, nproj, nslater):
    with open(path) as stream:
        values = [float(word) for word in stream.read().split()]
    if len(values) != 6 + 3 * (nproj + nslater):
        raise AssertionError("{}: unexpected initial parameter count".format(path))
    triples = [complex(values[6 + 3 * n], values[7 + 3 * n])
               for n in range(nproj + nslater)]
    return triples[:nproj], triples[nproj:]


def parse_audit(path):
    header = {}
    orbital = {}
    trans = {}
    weights = {}
    with open(path) as stream:
        for line in stream:
            words = line.split()
            if not words:
                continue
            if words[0] == "ORBITAL":
                orbital[(int(words[1]), int(words[2]))] = (int(words[3]),
                                                           int(words[4]))
            elif words[0] == "TRANS":
                trans[(int(words[1]), int(words[2]))] = (int(words[3]),
                                                         int(words[4]))
            elif words[0] == "TRANSWEIGHT":
                weights[int(words[1])] = complex(float(words[2]),
                                                 float(words[3]))
            else:
                header[words[0]] = int(words[1])
    return header, orbital, trans, weights


# ---------------------------------------------------------------------------
# Pinned four-site fixture
# ---------------------------------------------------------------------------

RING = Lattice((4,), (AP,))


def ring_parameters(keys):
    values = {
        # Every even sector N=0..8 carries >= 5% weight, and the seam bond
        # carries enough coherence that dropping its sign moves the energy by
        # many standard errors (selected by a seeded search over fixtures).
        ("ud", (0,)): -0.113 - 0.016j,
        ("ud", (1,)): -0.247 + 0.169j,
        ("ud", (2,)): -0.278 - 0.182j,
        ("ud", (3,)): -0.023 + 0.212j,
        ("ss", 0, (1,)): 0.27 - 0.128j,
        ("ss", 0, (2,)): -0.053 + 0.428j,
        ("ss", 1, (1,)): 0.434 + 0.119j,
        ("ss", 1, (2,)): -0.395 - 0.044j,
    }
    return tuple(values[key] for key in keys)


def ring_model(delta=None):
    table, keys = orbital_geometry(RING)
    # The pair source straddles the anti-periodic seam: c+_{3,up} c+_{0,dn}.
    return Model(
        RING, table, ring_parameters(keys),
        gutzwiller=-0.31, jastrow=(0.23, -0.12),
        hopping=1.11 - 0.13j, mu=0.35,
        coulomb_intra=0.35, coulomb_inter=0.27,
        delta=delta, anomalous_pair=(3, 4) if delta is not None else None,
        anomalous_green_pairs=((3, 4), (0, 5)))


def drop_boundary_phase(lattice):
    """Same geometry with every anti-periodic direction made periodic."""
    return Lattice(lattice.lengths, tuple(P for unused in lattice.lengths))


def self_test():
    model = ring_model()
    table, keys = orbital_geometry(RING)
    if len(keys) != 8 or len(table) != 28:
        raise AssertionError("4-site AP geometry inventory changed")
    negative = negative_sign_count(table)
    if negative != 6:
        raise AssertionError("4-site AP upper-triangle negative count is {}"
                             .format(negative))
    rng = random.Random(1729)
    for lattice in (RING, Lattice((4, 2), (AP, P)), Lattice((4, 4), (AP, AP)),
                    Lattice((4,), (P,))):
        lattice_table, lattice_keys = orbital_geometry(lattice)
        values = [complex(rng.uniform(-1, 1), rng.uniform(-1, 1))
                  for unused in lattice_keys]
        for number, key in enumerate(lattice_keys):
            if key == ("zero",):
                values[number] = 0.0j
        matrix = pair_matrix(lattice_table, values, lattice.norbital)
        defect = translation_defect(lattice, matrix)
        if defect > 1.0e-14:
            raise AssertionError("{} orbitals are not translation covariant: {}"
                                 .format(lattice.label(), defect))
        hopping = hopping_coefficients(lattice, 0.7 - 0.3j)
        if hopping_translation_defect(lattice, hopping) > 1.0e-14:
            raise AssertionError("{} hopping is not translation covariant"
                                 .format(lattice.label()))
    # Dropping the seam phase must break covariance of the AP table.
    periodic = drop_boundary_phase(RING)
    wrong = pair_matrix(orbital_geometry(RING)[0], model.parameters, 8)
    if translation_defect(periodic, wrong) < 1.0e-3:
        raise AssertionError("AP/PBC covariance check is vacuous")

    observable = exact_observables(model)
    if abs(observable["energy"].imag) > 1.0e-12:
        raise AssertionError("Hamiltonian inventory is not Hermitian")
    for count in range(0, 9, 2):
        if observable["sectors"].get(count, 0.0) < 0.04:
            raise AssertionError("sector N={} is not sampled: {}".format(
                count, observable["sectors"]))
    delta_model = ring_model(delta=0.33 - 0.21j)
    observable_delta = exact_observables(delta_model)
    if abs(observable_delta["energy"].imag) > 1.0e-12:
        raise AssertionError("anomalous inventory is not Hermitian")
    if abs(observable_delta["number"] - observable["number"]) > 1.0e-12:
        raise AssertionError("fixed-f particle number depends on delta")
    # The expanded reference parameterization is the same physical state.
    expanded, unused = expanded_table(model.table)
    wave_ap = wavefunction(model)
    wave_expanded = wavefunction(model, expanded_parameters(model.table,
                                                            model.parameters),
                                 expanded)
    worst = max(abs(wave_ap[state] - wave_expanded[state]) for state in wave_ap)
    if worst > 1.0e-14:
        raise AssertionError("expanded reference state differs: {}".format(worst))
    for parameter_index, imaginary in ((1, False), (1, True), (5, False)):
        exact = parameter_gradient(model, parameter_index, imaginary)
        fd1 = finite_difference(model, parameter_index, imaginary, 1.0e-5)
        fd2 = finite_difference(model, parameter_index, imaginary, 5.0e-6)
        tolerance = 1.0e-7 + 1.0e-5 * abs(exact)
        if abs(fd1 - exact) > tolerance or abs(fd2 - exact) > tolerance:
            raise AssertionError("finite-difference gradient mismatch")
    return observable


if __name__ == "__main__":
    try:
        values = self_test()
        print("grand-canonical AP oracle passed: E={:.12g} N={:.12g} "
              "varN={:.12g} sectors={}".format(
                  values["energy"].real, values["number"],
                  values["variance_number"],
                  dict((k, round(v, 4)) for k, v in
                       sorted(values["sectors"].items()))))
    except Exception as error:
        print("ERROR: {}".format(error))
        sys.exit(1)
