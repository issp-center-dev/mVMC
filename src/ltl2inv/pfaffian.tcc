/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "lapack2eigen.hh"


template <typename T_>
struct SignedLogScale
{
  using T = T_;

  T_  val;
  int sgn;
};


template <typename Mat, typename Vec,
          typename Amp = typename Mat::Scalar>
SignedLogScale<Amp> utu2logpfa(const Mat &A, const Vec &iPiv)
{
  SignedLogScale<Amp> logpfa{ 0.0, 1 };

  for (int i = 0; i < A.rows(); i += 2) {
    logpfa.val += std::log(std::abs(Amp( A(i, i+1) )));
    logpfa.sgn *= A(i, i+1) > 0 ? 1 : -1;

    if (iPiv(i)  -1 != i  ) logpfa.sgn *= -1;
    if (iPiv(i+1)-1 != i+1) logpfa.sgn *= -1;
  }
  return logpfa;
}

// Additional parameter for return type deduction.
template <typename Mat, typename Vec, typename Ret>
Ret utu2logpfa(const Mat &A, const Vec &iPiv, Ret &ret)
{ return utu2logpfa<Mat, Vec, typename Ret::T>(A, iPiv); }



template <typename Mat, typename Vec,
          typename Amp = typename Mat::Scalar>
Amp ltl2pfa(const Mat &A, const Vec &iPiv)
{
  Amp pfa = 1.0;

  for (int i = 0; i < A.rows(); i += 2) {
    pfa *= -A(i+1, i);

    if (iPiv(i)  -1 != i  ) pfa = -pfa;
    if (iPiv(i+1)-1 != i+1) pfa = -pfa;
  }
  return pfa;
}

template <typename Mat, typename Vec,
          typename Amp = typename Mat::Scalar>
Amp utu2pfa(const Mat &A, const Vec &iPiv)
{
  Amp pfa = 1.0;

  for (int i = 0; i < A.rows(); i += 2) {
    pfa *= A(i, i+1);

    if (iPiv(i)  -1 != i  ) pfa = -pfa;
    if (iPiv(i+1)-1 != i+1) pfa = -pfa;
  }
  return pfa;
}

template <typename Mat, typename Vec,
          typename Amp = typename Mat::Scalar>
Amp ltl2pfa(const char &uplo, const Mat &A, const Vec &iPiv)
{
  switch (uplo) {
  case 'l':
  case 'L':
    return ltl2pfa(A, iPiv);
  case 'u':
  case 'U':
    return utu2pfa(A, iPiv);
  default:
    return 0.0;
  }
}

