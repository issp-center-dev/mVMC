/**
 * \file skr2k.tcc
 * Dispatcher and minimal implementation for SKR2k.
 * Minimal implementation is used when M is extremely small.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#include "blalink.hh"
#include "blalink_gemmt.hh"
#include "colmaj.hh"
#include <iostream>


template <typename T>
inline void skr2k(uplo_t uploc,
                  trans_t transab,
                  dim_t m, dim_t k,
                  T alpha,
                  T *_A, inc_t ldA,
                  T *_B, inc_t ldB,
                  T beta,
                  T *_C, inc_t ldC)
{
    T one = 1.0;
    if (transab == BLIS_NO_TRANSPOSE) {
        gemmt<T>(uploc, BLIS_NO_TRANSPOSE, BLIS_TRANSPOSE, m, k, alpha, _A, ldA, _B, ldB, beta, _C, ldC);
        gemmt<T>(uploc, BLIS_NO_TRANSPOSE, BLIS_TRANSPOSE, m, k,-alpha, _B, ldB, _A, ldA, one , _C, ldC);
    } else {
        gemmt<T>(uploc, BLIS_TRANSPOSE, BLIS_NO_TRANSPOSE, m, k, alpha, _A, ldA, _B, ldB, beta, _C, ldC);
        gemmt<T>(uploc, BLIS_TRANSPOSE, BLIS_NO_TRANSPOSE, m, k,-alpha, _B, ldB, _A, ldA, one , _C, ldC);
    }
}

