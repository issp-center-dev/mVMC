/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "colmaj.hh"
#include "blalink.hh"
#include "ilaenv_lauum.hh"


template <typename T>
inline void trmmt(uplo_t uploab, diag_t diaga, int n, T alpha, colmaj<T> &A, colmaj<T> &B)
{
  int npanel = ilaenv_lauum<T>(uploab, n);

  if (uploab == BLIS_UPPER) {
    // A upper => A' lower -> write to UPPER part of B.
    for (int j = 0; j < n; j += npanel) {
      int nloc = std::min<int>(npanel, n - j);

      trmm<T>(BLIS_LEFT, BLIS_UPPER, BLIS_TRANSPOSE, diaga,
              j+nloc, nloc, alpha, &A(0, 0), A.ld, &B(0, j), B.ld);
    }
  } else {
    // A lower => A' upper -> write to LOWER part of B.
    for (int j = 0; j < n; j += npanel) {
      int nloc = std::min<int>(npanel, n - j);

      trmm<T>(BLIS_LEFT, BLIS_LOWER, BLIS_TRANSPOSE, diaga,
              n-j, nloc, alpha, &A(j, j), A.ld, &B(j, j), B.ld);
    }
  }
}

