/**
 * \copyright Copyright (c) Dept. Phys., Univ. Tokyo
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "colmaj.tcc"
#include "blalink.hh"
#include <vector>

template <typename T>
void skmv(uplo_t uploA,
          dim_t m,
          T alpha,
          T *_A, inc_t ldA,
          T *X,
          T *Y) {
  using namespace std;
  colmaj<T> A(_A, ldA);

  // Way 1:
  // skmm(BLIS_LEFT, uploA, BLIS_NO_CONJUGATE, BLIS_NO_TRANSPOSE, m, 1,
  //      T(1.0), &A(0, 0), A.ld, &X[0], ldA, T(0.0), &Y[0], ldA);

  // Way 2:
  vector<T> Z(m);
  for (dim_t i = 0; i < m; ++i) {
    Z[i] = X[i];
    Y[i] = X[i];
  }
  trmv(uploA, BLIS_NO_TRANSPOSE,
       m, alpha, &A(0, 0), A.ld, &Y[0], 1);
  trmv(uploA, BLIS_TRANSPOSE,
       m, alpha, &A(0, 0), A.ld, &Z[0], 1);
  for (dim_t i = 0; i < m; ++i)
    Y[i] -= Z[i];

  // Way 3:
  // for (dim_t i = 0; i < m; ++i)
  //   Y[i] = 0.0;
  // for (dim_t j = 0; j < m; ++j) {
  //   T x_j = X[j];
  //   for (dim_t i = 0; i < j; ++i) {
  //     T A_ij = A(i, j);
  //     Y[i] += A_ij * x_j;
  //     Y[j] -= A_ij * X[i];
  //   }
  // }
  // for (dim_t i = 0; i < m; ++i)
  //   Y[i] *= alpha;
  
}

