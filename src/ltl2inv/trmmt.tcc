/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "lapack2eigen.hh"
#include "ilaenv_lauum.hh"


template <typename Mat, typename Mat_>
inline void trmmt(const char &uploAB, const char &diagA, const typename Mat::Scalar &alpha, Mat &A, Mat_ &B)
{
  using namespace Eigen;
  using T = typename Mat::Scalar;
  int n = A.rows();
  int npanel = ilaenv_lauum<T>(uploAB, n);

  if (uploAB == 'U' || uploAB == 'u') {
    // A upper => A' lower -> write to UPPER part of B.
    for (int j = 0; j < n; j += npanel) {
      int nloc = std::min<int>(npanel, n - j);

      auto Bslice = B(seq(0, j+nloc-1), seq(j, j+nloc-1));
      if (diagA == 'U' || diagA == 'u')
        Bslice = A(seq(0, j+nloc-1), seq(0, j+nloc-1)).template triangularView<UnitUpper>().transpose() * Bslice * alpha;
      else
        Bslice = A(seq(0, j+nloc-1), seq(0, j+nloc-1)).template triangularView<Upper>().transpose() * Bslice * alpha;
    }
  } else {
    // A lower => A' upper -> write to LOWER part of B.
    for (int j = 0; j < n; j += npanel) {
      int nloc = std::min<int>(npanel, n - j);

      auto Bslice = B(seq(j, n-1), seq(j, j+nloc-1));
      if (diagA == 'U' || diagA == 'u')
        Bslice = A(seq(j, n-1), seq(j, n-1)).template triangularView<UnitLower>().transpose() * Bslice * alpha;
      else
        Bslice = A(seq(j, n-1), seq(j, n-1)).template triangularView<Lower>().transpose() * Bslice * alpha;
    }
  }
}

