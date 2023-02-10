/**
 * \copyright Copyright (c) Dept. Phys., Univ. Tokyo
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#pragma once
#include <Eigen/Dense>
#include <random>

namespace vmc
{
namespace orbital
{
template <typename matrix_t_>
struct orbital_mat
{
  using matrix_t = matrix_t_;
  using T = typename matrix_t::Scalar;
  using index_t = int32_t;

  const char uplo;
  const index_t norb_;
  matrix_t &X; ///< norb*norb.

  orbital_mat(char uplo_, index_t norb__, matrix_t &X_)
  : uplo(uplo_), norb_(norb__), X(X_) { }

  virtual index_t norb() { return norb_; }

  virtual void randomize(double amplitude, unsigned seed) {
    using namespace std;
    mt19937_64 rng(seed);
    uniform_real_distribution<double> dist(-0.1, 1.0);

    for (index_t j = 0; j < norb_; ++j) {
      for (index_t i = 0; i < j; ++i) {
        X(i, j) = T(dist(rng)) * amplitude;
        X(j, i) = -X(i, j);
      }
      X(j, j) = T(0.0);
    }
  }

  void randomize(double amplitude) { randomize(amplitude, 511); }

  virtual index_t pack_idx(index_t osi, index_t osj) {
    abort();
  }

  virtual void unpack_idx(index_t &osi, index_t &osj, index_t idx) {
    abort();
  }

  virtual T operator()(index_t osi, index_t osj) {
    if (osi == osj)
      return T(0.0);

    switch (uplo) {
    case 'U':
    case 'u':
      if (osi < osj)
        return X(osi, osj);
      else
        return -X(osj, osi);
    
    case 'L':
    case 'l':
      if (osi > osj)
        return X(osi, osj);
      else
        return -X(osj, osi);

    default:
      return X(osi, osj);
    }
  }
};
}
}
