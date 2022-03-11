/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "colmaj.hh"
#include "blalink.hh"

template <typename T>
T ltl2pfa(int n, T *A_, int ldA, int *iPiv)
{
  T pfa = 1.0;
  colmaj<T> A(A_, ldA);

  for (int i = 0; i < n; i += 2) {
    pfa *= -A(i+1, i);

    if (iPiv[i+1]-1 != i+1) pfa = -pfa;
  }
  return pfa;
}

