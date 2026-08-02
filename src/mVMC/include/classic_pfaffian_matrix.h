/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_CLASSIC_PFAFFIAN_MATRIX_H
#define MVMC_CLASSIC_PFAFFIAN_MATRIX_H

#include "absolute_pfaffian.h"

#include <complex.h>

/*
 * Build the classic non-FSZ nsize-by-nsize skew matrix for one QP component.
 * Slater storage is a column-major (2*nsite)-square matrix.  The first and
 * second nsize/2 electron indices represent up and down spin respectively.
 * No input is modified and no allocation is performed.
 */
MVMCPfaffianStatus mvmc_classic_pfaffian_build_real_matrix(
    const double *slater, int nsite, int nsize, const int *ele_idx,
    double *matrix);

MVMCPfaffianStatus mvmc_classic_pfaffian_build_complex_matrix(
    const double complex *slater, int nsite, int nsize, const int *ele_idx,
    double complex *matrix);

#endif /* MVMC_CLASSIC_PFAFFIAN_MATRIX_H */
