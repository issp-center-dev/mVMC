from __future__ import print_function

import argparse
import cmath
import hashlib
import math
import os
import struct
import sys


def fixed_sz_basis(nsite, nup, ndown):
    states = []
    for state in range(1 << (2 * nsite)):
        count_up = sum((state >> (2 * site)) & 1 for site in range(nsite))
        count_down = sum((state >> (2 * site + 1)) & 1 for site in range(nsite))
        if count_up == nup and count_down == ndown:
            states.append(state)
    return states


def permutation_sign(values):
    inversions = 0
    for left in range(len(values)):
        for right in range(left + 1, len(values)):
            inversions += values[left] > values[right]
    return -1 if inversions % 2 else 1


def translate_state(state, nsite, antiperiodic):
    occupied = [orbital for orbital in range(2 * nsite)
                if (state >> orbital) & 1]
    transformed = []
    sign = 1
    for orbital in occupied:
        site, spin = divmod(orbital, 2)
        if antiperiodic and site == nsite - 1:
            sign = -sign
        transformed.append(2 * ((site + 1) % nsite) + spin)
    sign *= permutation_sign(transformed)
    transformed_state = sum(1 << orbital for orbital in transformed)
    return transformed_state, sign


def read_hphi_eigenvector(path):
    with open(path, "rb") as source:
        data = source.read()
    if len(data) < 28:
        raise ValueError("HPhi eigenvector file is too short")
    iteration = struct.unpack_from("<i", data, 0)[0]
    dimension = struct.unpack_from("<Q", data, 4)[0]
    expected_size = 12 + 16 * (dimension + 1)
    if len(data) != expected_size:
        raise ValueError(
            "unexpected HPhi eigenvector size: got {} expected {}".format(
                len(data), expected_size))
    dummy = complex(*struct.unpack_from("<dd", data, 12))
    if abs(dummy) > 1.0e-14:
        raise ValueError("HPhi one-based eigenvector dummy entry is nonzero")
    vector = [complex(*struct.unpack_from("<dd", data, 28 + 16 * idx))
              for idx in range(dimension)]
    return iteration, vector


def check_sector(path, nsite, nup, ndown, antiperiodic,
                 expected_eigenvalue, residual_tolerance):
    iteration, vector = read_hphi_eigenvector(path)
    basis = fixed_sz_basis(nsite, nup, ndown)
    if len(vector) != len(basis):
        raise ValueError("HPhi vector dimension does not match fixed-Sz basis")
    index = dict((state, idx) for idx, state in enumerate(basis))
    norm = sum(abs(value) ** 2 for value in vector)
    if abs(norm - 1.0) > 1.0e-10:
        raise ValueError("HPhi eigenvector norm mismatch: {:.17e}".format(norm))

    transformed = [0.0j for _ in vector]
    for idx, state in enumerate(basis):
        target, sign = translate_state(state, nsite, antiperiodic)
        transformed[index[target]] += sign * vector[idx]
    eigenvalue = sum(value.conjugate() * moved
                     for value, moved in zip(vector, transformed))
    residual = math.sqrt(sum(abs(moved - eigenvalue * value) ** 2
                             for value, moved in zip(vector, transformed)))
    if abs(eigenvalue - expected_eigenvalue) > residual_tolerance:
        raise ValueError(
            "translation eigenvalue mismatch: got {} expected {}".format(
                eigenvalue, expected_eigenvalue))
    if residual > residual_tolerance:
        raise ValueError("translation residual exceeds tolerance: {:.3e}".format(
            residual))
    with open(path, "rb") as source:
        digest = hashlib.sha256(source.read()).hexdigest()
    return {
        "iteration": iteration,
        "dimension": len(vector),
        "norm": norm,
        "eigenvalue": eigenvalue,
        "residual": residual,
        "sha256": digest,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Verify the translation sector of an HPhi fixed-Sz eigenvector")
    parser.add_argument("vector")
    parser.add_argument("--nsite", type=int, required=True)
    parser.add_argument("--nup", type=int, required=True)
    parser.add_argument("--ndown", type=int, required=True)
    parser.add_argument("--antiperiodic", action="store_true")
    parser.add_argument("--expected-real", type=float, required=True)
    parser.add_argument("--expected-imag", type=float, default=0.0)
    parser.add_argument("--residual-tolerance", type=float, default=2.0e-6)
    args = parser.parse_args()
    try:
        result = check_sector(
            args.vector, args.nsite, args.nup, args.ndown,
            args.antiperiodic,
            complex(args.expected_real, args.expected_imag),
            args.residual_tolerance)
    except (IOError, OSError, ValueError) as error:
        print("ERROR: {}".format(error))
        return 1
    print("iteration {}".format(result["iteration"]))
    print("dimension {}".format(result["dimension"]))
    print("norm {:.17e}".format(result["norm"]))
    print("translation_eigenvalue_real {:.17e}".format(
        result["eigenvalue"].real))
    print("translation_eigenvalue_imag {:.17e}".format(
        result["eigenvalue"].imag))
    print("translation_eigenvalue_abs {:.17e}".format(
        abs(result["eigenvalue"])))
    print("translation_residual_l2 {:.17e}".format(result["residual"]))
    print("eigenvector_sha256 {}".format(result["sha256"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
