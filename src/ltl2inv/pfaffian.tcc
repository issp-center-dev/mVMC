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
    if (iPiv[i+2]-1 != i+2) pfa = -pfa;
  }
  return pfa;
}

template <typename T>
T utu2pfa(int n, T *A_, int ldA, int *iPiv)
{
  T pfa = 1.0;
  colmaj<T> A(A_, ldA);

  for (int i = 0; i < n; i += 2) {
    pfa *= A(i, i+1);

    if (iPiv[i+1]-1 != i+1) pfa = -pfa;
    if (iPiv[i+2]-1 != i+2) pfa = -pfa;
  }
  return pfa;
}

template <typename T>
T ltl2pfa(uplo_t uplo, int n, T *A_, int ldA, int *iPiv)
{
  switch (uplo) {
  case BLIS_LOWER:
    return ltl2pfa(n, A_, ldA, iPiv);
  case BLIS_UPPER:
    return utu2pfa(n, A_, ldA, iPiv);
  default:
    return 0.0;
  }
}

