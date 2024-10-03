/**
 * \file skslc.tcc
 * Get a slice from antisymmetric matrix.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "lapack2eigen.hh"

template<typename Mat>
inline Eigen::Matrix<typename Mat::Scalar, Eigen::Dynamic, 1>
  skslc(const char &uploA, unsigned i, const Mat &A)
{
    Eigen::Matrix<typename Mat::Scalar, Eigen::Dynamic, 1> x(A.rows());

    x[i] = 0.0;
    switch (uploA) {
    case 'U':
    case 'u':
        for (unsigned j = 0; j < i; ++j)
            x[j] =  A(j, i);
        for (unsigned j = i+1; j < A.rows(); ++j)
            x[j] = -A(i, j);
        break;

    default:
        std::cerr << "SKSLC: Lower triangular storage not implemented. Sorry." << std::endl;
        break;;
    }
    return x;
}

