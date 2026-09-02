from __future__ import print_function

import cmath
import math
import sys

import numpy as np


NSITE = 2
NORBITAL = 2 * NSITE
EVEN_BASIS = tuple(
    state for state in range(1 << NORBITAL)
    if bin(state).count("1") % 2 == 0
)


def pfaffian(matrix):
    n = len(matrix)
    if n == 0:
        return 1.0 + 0.0j
    if n % 2:
        return 0.0 + 0.0j
    if n == 2:
        return matrix[0][1]
    value = 0.0 + 0.0j
    for column in range(1, n):
        minor = [
            [matrix[row][col] for col in range(n) if col not in (0, column)]
            for row in range(n) if row not in (0, column)
        ]
        value += ((-1) ** (column + 1)) * matrix[0][column] * pfaffian(minor)
    return value


def occupied(state, orbital):
    return (state >> orbital) & 1


def annihilate(state, orbital):
    if not occupied(state, orbital):
        return None
    lower = state & ((1 << orbital) - 1)
    sign = -1 if bin(lower).count("1") % 2 else 1
    return state & ~(1 << orbital), sign


def create(state, orbital):
    if occupied(state, orbital):
        return None
    lower = state & ((1 << orbital) - 1)
    sign = -1 if bin(lower).count("1") % 2 else 1
    return state | (1 << orbital), sign


def apply_ops(state, operators):
    sign = 1
    current = state
    for kind, orbital in reversed(operators):
        result = create(current, orbital) if kind == "c" else annihilate(current, orbital)
        if result is None:
            return None
        current, step_sign = result
        sign *= step_sign
    return current, sign


def one_body(out_orbital, in_orbital):
    return (("c", out_orbital), ("a", in_orbital))


def nbody(fused_pairs):
    operators = []
    for out_orbital, in_orbital in fused_pairs:
        operators.extend(one_body(out_orbital, in_orbital))
    return tuple(operators)


def pair_create(first, second):
    return (("c", first), ("c", second))


def pair_remove(first, second):
    return (("a", first), ("a", second))


def adjoint(operators):
    return tuple(("a" if kind == "c" else "c", orbital)
                 for kind, orbital in reversed(operators))


ANOMALOUS_G_OPERATORS = (
    (1, 0, 3, pair_create(0, 3)),
    (0, 3, 0, pair_remove(3, 0)),
    (1, 3, 0, pair_create(3, 0)),
    (0, 0, 3, pair_remove(0, 3)),
)


def default_parameters():
    # Raw OrbitalGeneral parameters.  The largest magnitude is deliberately
    # below D_AmpMax=4 so the GC no-rescale guard is observable end to end.
    return (
        2.0 + 0.0j,
        0.42 - 0.31j,
        -0.27 + 0.36j,
        0.33 + 0.18j,
        -0.29 - 0.22j,
        0.21 + 0.41j,
    )


def slater_matrix(parameters):
    matrix = [[0.0j for unused in range(NORBITAL)] for unused in range(NORBITAL)]
    index = 0
    for row in range(NORBITAL):
        for column in range(row + 1, NORBITAL):
            # OrbitalGeneral expands one upper-triangle parameter as
            # f_ij-f_ji = 2*f_ij.  CalculateMAllGC forms its Pfaffian matrix
            # from -SlaterElm, so the minus sign is part of the mVMC BCS
            # input convention.  It is invisible to number-conserving
            # observables but fixes the relative phase between N sectors.
            matrix[row][column] = -2.0 * parameters[index]
            matrix[column][row] = -matrix[row][column]
            index += 1
    return matrix


def jastrow_count(state):
    n0 = occupied(state, 0) + occupied(state, NSITE)
    n1 = occupied(state, 1) + occupied(state, NSITE + 1)
    return (n0 - 1) * (n1 - 1)


def wavefunction(parameters=None, jastrow=0.0):
    if parameters is None:
        parameters = default_parameters()
    matrix = slater_matrix(parameters)
    result = {}
    for state in EVEN_BASIS:
        orbitals = [orbital for orbital in range(NORBITAL) if occupied(state, orbital)]
        restricted = [[matrix[row][column] for column in orbitals] for row in orbitals]
        result[state] = pfaffian(restricted) * cmath.exp(jastrow * jastrow_count(state))
    return result


def norm(wave):
    return sum(abs(amplitude) ** 2 for amplitude in wave.values())


def expectation(wave, operators):
    denominator = norm(wave)
    numerator = 0.0 + 0.0j
    for state, amplitude in wave.items():
        transformed = apply_ops(state, operators)
        if transformed is None:
            continue
        target, sign = transformed
        numerator += wave.get(target, 0.0j).conjugate() * sign * amplitude
    return numerator / denominator


def model_terms(mu=0.0, delta=None):
    up0, up1, down0, down1 = 0, 1, 2, 3
    terms = []

    # Trans, including diagonal chemical potential and complex hopping.
    for orbital in range(NORBITAL):
        terms.append(("trans_mu", -mu, one_body(orbital, orbital)))
    hopping = 0.37 + 0.19j
    for left, right in ((up0, up1), (down0, down1)):
        terms.append(("trans", -hopping, one_body(left, right)))
        terms.append(("trans", -hopping.conjugate(), one_body(right, left)))

    # Density interactions and Hund's diagonal same-spin contribution.
    terms.append(("coulomb_intra", 0.73, nbody(((up0, up0), (down0, down0)))))
    terms.append(("coulomb_intra", -0.28, nbody(((up1, up1), (down1, down1)))))
    for first in (up0, down0):
        for second in (up1, down1):
            terms.append(("coulomb_inter", 0.21, nbody(((first, first), (second, second)))))
    terms.append(("hund", -0.17, nbody(((up0, up0), (up1, up1)))))
    terms.append(("hund", -0.17, nbody(((down0, down0), (down1, down1)))))

    # Pair hopping, exchange, and a general complex spin-flip interaction.
    terms.append(("pair_hop", 0.13, nbody(((up0, up1), (down0, down1)))))
    terms.append(("pair_hop", 0.13, nbody(((up1, up0), (down1, down0)))))
    terms.append(("exchange", -0.11, nbody(((up0, up1), (down1, down0)))))
    terms.append(("exchange", -0.11, nbody(((down0, down1), (up1, up0)))))
    general = -0.16 + 0.09j
    general_ops = nbody(((up0, down1), (up1, down0)))
    terms.append(("inter_all", general, general_ops))
    terms.append(("inter_all", general.conjugate(), adjoint(general_ops)))

    # A genuine input-order-three expression which reduces through coincident
    # factors on this four-orbital Hilbert space.
    nbody_ops = nbody(((up0, up1), (down0, down0), (up1, up1)))
    terms.append(("nbody_inter_all", 0.07, nbody_ops))
    terms.append(("nbody_inter_all", 0.07, adjoint(nbody_ops)))
    if delta is not None:
        create_ops = pair_create(0, 3)
        terms.append(("anomalous", delta, create_ops))
        terms.append(("anomalous", delta.conjugate(),
                      adjoint(create_ops)))
    return tuple(terms)


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


def exact_observables(parameters=None, jastrow=0.0, mu=0.0, delta=None):
    wave = wavefunction(parameters, jastrow)
    denominator = norm(wave)
    action = hamiltonian_action(wave, model_terms(mu, delta))
    energy = sum(wave[state].conjugate() * action[state] for state in wave) / denominator
    energy2 = sum(abs(action[state]) ** 2 for state in wave) / denominator
    number = sum(abs(amplitude) ** 2 * bin(state).count("1")
                 for state, amplitude in wave.items()) / denominator
    number2 = sum(abs(amplitude) ** 2 * bin(state).count("1") ** 2
                  for state, amplitude in wave.items()) / denominator
    greens1 = dict(
        ((out_orbital, in_orbital), expectation(wave, one_body(out_orbital, in_orbital)))
        for out_orbital in range(NORBITAL) for in_orbital in range(NORBITAL)
    )
    anomalous_g = dict(
        ((operator_type, first, second), expectation(wave, operators))
        for operator_type, first, second, operators in ANOMALOUS_G_OPERATORS
    )
    return {
        "wave": wave,
        "energy": energy,
        "energy2": energy2,
        "number": number,
        "number2": number2,
        "variance_number": number2 - number * number,
        "greens1": greens1,
        "anomalous_g": anomalous_g,
    }


def sector_hamiltonian(basis, mu=0.0, delta=None):
    basis_index = dict((state, index) for index, state in enumerate(basis))
    matrix = np.zeros((len(basis), len(basis)), dtype=np.complex128)
    for column, state in enumerate(basis):
        for unused_label, coefficient, operators in model_terms(mu, delta):
            transformed = apply_ops(state, operators)
            if transformed is None:
                continue
            target, sign = transformed
            if target in basis_index:
                matrix[basis_index[target], column] += coefficient * sign
    return matrix


def even_sector_eigenvalues(mu=0.0, delta=None):
    return np.linalg.eigvalsh(sector_hamiltonian(EVEN_BASIS, mu, delta))


def ed_ground_state(mu=0.0, delta=None):
    eigenvalues, eigenvectors = np.linalg.eigh(
        sector_hamiltonian(EVEN_BASIS, mu, delta))
    vector = eigenvectors[:, 0]
    wave = dict((state, vector[index])
                for index, state in enumerate(EVEN_BASIS))
    number = sum(abs(vector[index]) ** 2 * bin(state).count("1")
                 for index, state in enumerate(EVEN_BASIS))
    anomalous_g = dict(
        ((operator_type, first, second), expectation(wave, operators))
        for operator_type, first, second, operators in ANOMALOUS_G_OPERATORS
    )
    return {
        "energy": eigenvalues[0],
        "number": number,
        "anomalous_g": anomalous_g,
        "vector": vector,
    }


def spectra_match(first, second, tolerance):
    return len(first) == len(second) and all(
        abs(left - right) <= tolerance
        for left, right in zip(first, second)
    )


def energy(parameters, component_jastrow=0.0, mu=0.0, delta=None):
    return exact_observables(parameters, component_jastrow, mu,
                             delta)["energy"].real


def parameter_gradient(parameter_index, imaginary=False, jastrow=0.0,
                       mu=0.0, delta=None):
    parameters = list(default_parameters())
    wave = wavefunction(parameters, jastrow)
    action = hamiltonian_action(wave, model_terms(mu, delta))
    denominator = norm(wave)
    mean_h = sum(wave[state].conjugate() * action[state] for state in wave) / denominator
    epsilon = 1.0e-7
    plus = list(parameters)
    minus = list(parameters)
    parameter_step = 1j * epsilon if imaginary else epsilon
    plus[parameter_index] += parameter_step
    minus[parameter_index] -= parameter_step
    wave_plus = wavefunction(plus, jastrow)
    wave_minus = wavefunction(minus, jastrow)
    mean_o = 0.0j
    mean_oh = 0.0j
    for state, amplitude in wave.items():
        derivative = (wave_plus[state] - wave_minus[state]) / (2.0 * epsilon)
        observable = derivative / amplitude if amplitude != 0.0 else 0.0j
        probability = abs(amplitude) ** 2 / denominator
        local_energy = action[state] / amplitude if amplitude != 0.0 else 0.0j
        mean_o += probability * observable
        mean_oh += probability * observable.conjugate() * local_energy
    covariance = 2.0 * (mean_oh - mean_o.conjugate() * mean_h).real
    return covariance


def finite_difference(parameter_index, imaginary, epsilon, jastrow=0.0,
                      mu=0.0, delta=None):
    plus = list(default_parameters())
    minus = list(default_parameters())
    parameter_step = 1j * epsilon if imaginary else epsilon
    plus[parameter_index] += parameter_step
    minus[parameter_index] -= parameter_step
    return (energy(plus, jastrow, mu, delta) -
            energy(minus, jastrow, mu, delta)) / (2.0 * epsilon)


def self_test():
    if len(EVEN_BASIS) != 8:
        raise AssertionError("Nsite=2 even basis must contain eight states")
    observable = exact_observables(jastrow=0.23, mu=0.4)
    if observable["variance_number"] <= 0.0:
        raise AssertionError("fixture must have particle-number fluctuations")
    # Hermiticity of the full term inventory is exposed by a real energy.
    if abs(observable["energy"].imag) > 1.0e-12:
        raise AssertionError("Hamiltonian inventory is not Hermitian")
    observable_delta = exact_observables(jastrow=0.23, mu=0.4,
                                         delta=0.4)
    if abs(observable_delta["energy"].imag) > 1.0e-12:
        raise AssertionError("anomalous Hamiltonian inventory is not Hermitian")
    if not spectra_match(even_sector_eigenvalues(0.4, None),
                         even_sector_eigenvalues(0.4, 0.0), 1.0e-12):
        raise AssertionError("delta=0 spectrum changed")
    if abs(ed_ground_state(0.4, 0.4)["number"] -
           ed_ground_state(0.4, None)["number"]) <= 1.0e-3:
        raise AssertionError("ED ground state does not respond to pairing")
    if abs(observable_delta["number"] - observable["number"]) > 1.0e-12:
        raise AssertionError("fixed-f particle number depends on delta")
    if abs(observable_delta["energy"] - observable["energy"]) <= 1.0e-3:
        raise AssertionError("fixed-f energy does not respond to delta")
    for parameter_index, imaginary in ((5, False), (5, True)):
        exact = parameter_gradient(parameter_index, imaginary, 0.23, 0.4)
        fd1 = finite_difference(parameter_index, imaginary, 1.0e-5, 0.23, 0.4)
        fd2 = finite_difference(parameter_index, imaginary, 5.0e-6, 0.23, 0.4)
        tolerance = 1.0e-7 + 1.0e-5 * abs(exact)
        if abs(fd1 - exact) > tolerance or abs(fd2 - exact) > tolerance:
            raise AssertionError("finite-difference gradient mismatch")
        if abs(fd1 - fd2) > 2.0e-7 + 1.0e-5 * abs(exact):
            raise AssertionError("finite-difference epsilon mismatch")
    return observable


if __name__ == "__main__":
    try:
        values = self_test()
        print("grand-canonical exact oracle passed: E={:.12g} N={:.12g} varN={:.12g}".format(
            values["energy"].real, values["number"], values["variance_number"]))
    except Exception as error:
        print("ERROR: {}".format(error))
        sys.exit(1)
