/**
 * \file skr2k.tcc
 * Dispatcher and minimal implementation for SKR2k.
 * Minimal implementation is used when M is extremely small.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "lapack2eigen.hh"
#include <iostream>


template <typename Mat0, typename Mat1, typename Mat2>
inline void skr2k(const char &uploC,
                  const char &transAB,
                  const typename Mat2::Scalar &alpha,
                  const Mat0 &A,
                  const Mat1 &B,
                  const typename Mat2::Scalar &beta,
                        Mat2 &C)
{
    using T = typename Mat0::Scalar;

    l2e::le_mat_t<T> C_(C);
    if (transAB == 'N' || transAB == 'n') {
        gemmt<T>(uploC, 'N', 'T',  alpha, l2e::le_mat_t<T>(A), l2e::le_mat_t<T>(B), beta, C_);
        gemmt<T>(uploC, 'N', 'T', -alpha, l2e::le_mat_t<T>(B), l2e::le_mat_t<T>(A), T(1), C_);
    } else {
        gemmt<T>(uploC, 'T', 'N',  alpha, l2e::le_mat_t<T>(A), l2e::le_mat_t<T>(B), beta, C_);
        gemmt<T>(uploC, 'T', 'N', -alpha, l2e::le_mat_t<T>(B), l2e::le_mat_t<T>(A), T(1), C_);
    }
}

