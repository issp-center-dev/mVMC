/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "classic_pfaffian_matrix.h"

#include <limits.h>
#include <stddef.h>

static int valid_layout(int nsite, int nsize, const int *ele_idx) {
  int index;

  if (nsite <= 0 || nsite > INT_MAX / 2 || nsize <= 0 ||
      (nsize % 2) != 0 || nsize / 2 > nsite || ele_idx == NULL) {
    return 0;
  }
  for (index = 0; index < nsize; ++index) {
    if (ele_idx[index] < 0 || ele_idx[index] >= nsite) return 0;
  }
  return 1;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_build_real_matrix(
    const double *slater, int nsite, int nsize, const int *ele_idx,
    double *matrix) {
  int ne;
  int nsite2;
  int column, row;

  if (slater == NULL || matrix == NULL ||
      !valid_layout(nsite, nsize, ele_idx)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  ne = nsize / 2;
  nsite2 = 2 * nsite;
  for (column = 0; column < nsize; ++column) {
    const int slater_column = ele_idx[column] + (column / ne) * nsite;
    for (row = 0; row < nsize; ++row) {
      const int slater_row = ele_idx[row] + (row / ne) * nsite;
      matrix[row + (size_t)column * (size_t)nsize] =
          -slater[(size_t)slater_column * (size_t)nsite2 +
                  (size_t)slater_row];
    }
  }
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_build_complex_matrix(
    const double complex *slater, int nsite, int nsize,
    const int *ele_idx, double complex *matrix) {
  int ne;
  int nsite2;
  int column, row;

  if (slater == NULL || matrix == NULL ||
      !valid_layout(nsite, nsize, ele_idx)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  ne = nsize / 2;
  nsite2 = 2 * nsite;
  for (column = 0; column < nsize; ++column) {
    const int slater_column = ele_idx[column] + (column / ne) * nsite;
    for (row = 0; row < nsize; ++row) {
      const int slater_row = ele_idx[row] + (row / ne) * nsite;
      matrix[row + (size_t)column * (size_t)nsize] =
          -slater[(size_t)slater_column * (size_t)nsite2 +
                  (size_t)slater_row];
    }
  }
  return MVMC_PFAFFIAN_STATUS_OK;
}
