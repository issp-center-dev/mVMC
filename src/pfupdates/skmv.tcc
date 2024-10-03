/**
 * \copyright Copyright (c) Dept. Phys., Univ. Tokyo
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "lapack2eigen.hh"

template <typename Mat, typename Vec>
Eigen::Matrix<typename Mat::Scalar, Eigen::Dynamic, 1>
  skmv(const char &uploA,
       const typename Mat::Scalar &alpha,
       const Mat &A,
       const Vec &X) {
  using namespace Eigen;

  Eigen::Matrix<typename Mat::Scalar, Eigen::Dynamic, 1> Y;
  // Way 1:
  if (uploA == 'U' || uploA == 'u') {
    Y = A.template triangularView<StrictlyUpper>() * X - A.template triangularView<StrictlyUpper>().transpose() * X;
  } else {
    Y = A.template triangularView<StrictlyLower>() * X - A.template triangularView<StrictlyLower>().transpose() * X;
  }
  Y *= alpha;

  // Way 2:
  // for (int i = 0; i < m; ++i)
  //   Y[i] = 0.0;
  // for (int j = 0; j < m; ++j) {
  //   T x_j = X[j];
  //   for (int i = 0; i < j; ++i) {
  //     T A_ij = A(i, j);
  //     Y[i] += A_ij * x_j;
  //     Y[j] -= A_ij * X[i];
  //   }
  // }
  // for (int i = 0; i < m; ++i)
  //   Y[i] *= alpha;
  
  return Y;
}

