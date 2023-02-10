/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "lapack2eigen.hh"
#include "trmmt.tcc"
#include "error.hh"

template <typename Vec, typename Mat>
void sktdsmx(const Vec &vT, const Mat &B, Mat &C)
{
  int n = B.rows();
  assert_(B.cols() == n && C.cols() == n && C.rows() == n, "sktdsmx: Dimension mismatch.");

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

template <typename iVec, typename Vec, typename Mat>
void ltl2inv(Mat &A, const iVec &iPiv, Vec &vT, Mat &M)
{
  using namespace Eigen;
  using T = typename Mat::Scalar;
  int n = A.rows();
  int npanel_trmm = ilaenv_lauum<T>('L', n);
  int full = npanel_trmm <= 1 || npanel_trmm >= n;

  // Set M.
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      M(i, j) = T(!(i - j));

  l2e::le_mat_t<T> A_(A);
  l2e::le_mat_t<T> M_(M);
  l2e::le_mat_t<T> A_Lx(A(seq(1, n-1), seq(0, n-2)));
  l2e::le_mat_t<T> A_L(A(seq(2, n-1), seq(0, n-3)));
  l2e::le_mat_t<T> M_L(M(seq(2, n-1), seq(1, n-2)));

  trtri('L', 'U', A_Lx);
  lacpy('L', A_L, M_L);

  for (int i = 0; i < n-1; ++i)
    vT[i] = A(i+1, i);
  sktdsmx(vT, M, A);

  if (!full) {
    trmmt('L', 'U', T(1.0), M, A);
    for (int j = 0; j < n; ++j) {
      A(j, j) = 0.0;
      for (int i = 0; i < j; ++i)
        A(i, j) = -A(j, i);
    }
  }

  // In-place permute columns.
  for (int j = n-1; j >= 0; --j)
    if (iPiv[j]-1 != j)
      A(Eigen::all, j).swap(A(Eigen::all, iPiv[j] - 1));

  if (full)
    A = M.template triangularView<UnitLower>().transpose() * A;

  // In-place permute rows.
  for (int i = n-1; i >= 0; --i)
    if (iPiv[i]-1 != i)
      A(i, Eigen::all).swap(A(iPiv[i] - 1, Eigen::all));
}

template <typename iVec, typename Vec, typename Mat>
void utu2inv(Mat &A, const iVec &iPiv, Vec &vT, Mat &M)
{
  using namespace Eigen;
  using T = typename Mat::Scalar;
  int n = A.rows();
  int npanel_trmm = ilaenv_lauum<T>('U', n);
  int full = npanel_trmm <= 1 || npanel_trmm >= n;

  // Set M.
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      M(i, j) = T(!(i - j));

  l2e::le_mat_t<T> A_(A);
  l2e::le_mat_t<T> M_(M);
  l2e::le_mat_t<T> A_Ux(A(seq(0, n-2), seq(1, n-1)));
  l2e::le_mat_t<T> A_U(A(seq(0, n-3), seq(2, n-1)));
  l2e::le_mat_t<T> M_U(M(seq(0, n-3), seq(1, n-2)));

  trtri('U', 'U', A_Ux);
  lacpy('U', A_U, M_U);

  for (int i = 0; i < n-1; ++i)
    vT[i] = -A(i, i+1);
  sktdsmx(vT, M, A);

  if (!full) {
    trmmt('U', 'U', T(1.0), M, A);
      for (int j = 0; j < n; ++j) {
        A(j, j) = 0.0;
        for (int i = j + 1; i < n; ++i)
          A(i, j) = -A(j, i);
      }
  }

  // In-place permute columns.
  for (int j = 0; j < n; ++j)
    if (iPiv[j]-1 != j)
      A(Eigen::all, j).swap(A(Eigen::all, iPiv[j]-1));

  if (full)
    A = M.template triangularView<UnitUpper>().transpose() * A;

  // In-place permute rows.
  for (int i = 0; i < n; ++i)
    if (iPiv[i]-1 != i)
      A(i, Eigen::all).swap(A(iPiv[i]-1, Eigen::all));
}

template <typename iVec, typename Vec, typename Mat>
void ltl2inv(const char uplo, Mat &A, const iVec &iPiv, Vec &vT, Mat &M)
{
  switch (uplo) {
  case 'l':
  case 'L':
    ltl2inv(A, iPiv, vT, M); break;
  case 'u':
  case 'U':
    utu2inv(A, iPiv, vT, M); break;
  default:
    break;
  }
}

