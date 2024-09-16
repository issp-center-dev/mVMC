/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "colmaj.hh"
#include "blalink.hh"
#include "trmmt.tcc"

template <typename T>
void sktdsmx(int n, T *vT, T *B_, int ldB, T *C_, int ldC)
{
  colmaj<T> B(B_, ldB);
  colmaj<T> C(C_, ldC);

  for (int j = 0; j < n; ++j)
    C(1, j) = B(0, j) / -vT[0];
  for (int i = 2; i < n; i += 2)
    for (int j = 0; j < n; ++j)
      C(i+1, j) = (B(i, j) - C(i-1, j) * vT[i-1]) / -vT[i];

  for (int j = 0; j < n; ++j)
    C(n-2, j) = B(n-1, j) / vT[n-2];
  for (int i = n-3; i >= 0; i -= 2)
    for (int j = 0; j < n; ++j)
      C(i-1, j) = (B(i, j) + C(i+1, j) * vT[i]) / vT[i-1];
}

template <typename T>
void ltl2inv(int n, T *A_, int ldA, int *iPiv, T *vT, T *M_, int ldM)
{
  colmaj<T> A(A_, ldA);
  colmaj<T> M(M_, ldM);

  int npanel_trmm = ilaenv_lauum<T>(BLIS_LOWER, n);
  int full = npanel_trmm <= 1 || npanel_trmm >= n;

  // Set M.
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      M(i, j) = T(!(i - j));

  trtri(BLIS_LOWER, BLIS_UNIT_DIAG, n-1, &A(1, 0), ldA);
  lacpy(BLIS_LOWER, n-2, n-2, &A(2, 0), ldA, &M(2, 1), ldM);

  for (int i = 0; i < n-1; ++i)
    vT[i] = A(i+1, i);
  sktdsmx<T>(n, vT, &M(0, 0), ldM, &A(0, 0), ldA);

  if (!full) {
    trmmt(BLIS_LOWER, BLIS_UNIT_DIAG, n, T(1.0), M, A);
    for (int j = 0; j < n; ++j) {
      A(j, j) = 0.0;
      for (int i = 0; i < j; ++i)
        A(i, j) = -A(j, i);
    }
  }

  // In-place permute columns.
  for (int j = n-1; j >= 0; --j)
    if (iPiv[j]-1 != j)
      swap(n, &A(0, j), 1, &A(0, iPiv[j]-1), 1);

  if (full)
    trmm(BLIS_LEFT, BLIS_LOWER, BLIS_TRANSPOSE, BLIS_UNIT_DIAG,
         n, n, T(1.0), &M(0, 0), ldM, &A(0, 0), ldA);

  // In-place permute rows.
  for (int i = n-1; i >= 0; --i)
    if (iPiv[i]-1 != i)
      swap(n, &A(i, 0), ldA, &A(iPiv[i]-1, 0), ldA);
}

template <typename T>
void utu2inv(int n, T *A_, int ldA, int *iPiv, T *vT, T *M_, int ldM)
{
  colmaj<T> A(A_, ldA);
  colmaj<T> M(M_, ldM);

  int npanel_trmm = ilaenv_lauum<T>(BLIS_UPPER, n);
  int full = npanel_trmm <= 1 || npanel_trmm >= n;

  // Set M.
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      M(i, j) = T(!(i - j));

  trtri(BLIS_UPPER, BLIS_UNIT_DIAG, n-1, &A(0, 1), ldA);
  lacpy(BLIS_UPPER, n-2, n-2, &A(0, 2), ldA, &M(0, 1), ldM);

  for (int i = 0; i < n-1; ++i)
    vT[i] = -A(i, i+1);
  sktdsmx<T>(n, vT, &M(0, 0), ldM, &A(0, 0), ldA);

  if (!full) {
    trmmt(BLIS_UPPER, BLIS_UNIT_DIAG, n, T(1.0), M, A);
      for (int j = 0; j < n; ++j) {
        A(j, j) = 0.0;
        for (int i = j + 1; i < n; ++i)
          A(i, j) = -A(j, i);
      }
  }

  // In-place permute columns.
  for (int j = 0; j < n; ++j)
    if (iPiv[j]-1 != j)
      swap(n, &A(0, j), 1, &A(0, iPiv[j]-1), 1);

  if (full)
    trmm(BLIS_LEFT, BLIS_UPPER, BLIS_TRANSPOSE, BLIS_UNIT_DIAG,
         n, n, T(1.0), &M(0, 0), ldM, &A(0, 0), ldA);

  // In-place permute rows.
  for (int i = 0; i < n; ++i)
    if (iPiv[i]-1 != i)
      swap(n, &A(i, 0), ldA, &A(iPiv[i]-1, 0), ldA);
}

template <typename T>
void ltl2inv(uplo_t uplo, int n, T *A_, int ldA, int *iPiv, T *vT, T *M_, int ldM)
{
  switch (uplo) {
  case BLIS_LOWER:
    ltl2inv(n, A_, ldA, iPiv, vT, M_, ldM); break;
  case BLIS_UPPER:
    utu2inv(n, A_, ldA, iPiv, vT, M_, ldM); break;
  default:
    break;
  }
}

