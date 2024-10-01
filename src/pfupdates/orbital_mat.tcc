/**
 * \copyright Copyright (c) Dept. Phys., Univ. Tokyo
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include "blis.h"
#include "colmaj.hh"
#include <random>
template <typename T>
using matrix_t = colmaj<T>;

template <typename T>
struct orbital_mat
{
  const uplo_t uplo;
  const dim_t nsite;
  matrix_t<T> X; ///< nsite*nsite.

  orbital_mat(uplo_t uplo_, dim_t nsite_, T *X_, inc_t ldX)
      : uplo(uplo_), nsite(nsite_),
        X(X_, ldX) { }

  orbital_mat(uplo_t uplo_, dim_t nsite_, matrix_t<T> &X_)
  : uplo(uplo_), nsite(nsite_), X(X_) { }

  void randomize(double amplitude, unsigned seed) { 
    using namespace std;
    mt19937_64 rng(seed);
    uniform_real_distribution<double> dist(-0.1, 1.0);

    for (dim_t j = 0; j < nsite; ++j) {
      for (dim_t i = 0; i < j; ++i) {
        X(i, j) = T(dist(rng)) * amplitude;
        X(j, i) = -X(i, j);
      }
      X(j, j) = T(0.0);
    }
  }

  void randomize(double amplitude) { randomize(amplitude, 511); }

  T operator()(dim_t osi, dim_t osj) {
    if (osi == osj)
      return T(0.0);

    switch (uplo) {
    case BLIS_UPPER:
      if (osi < osj)
        return X(osi, osj);
      else
        return -X(osj, osi);
    
    case BLIS_LOWER:
      if (osi > osj)
        return X(osi, osj);
      else
        return -X(osj, osi);

    default:
      return X(osi, osj);
    }
  }
};
