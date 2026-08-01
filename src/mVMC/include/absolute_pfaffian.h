/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_ABSOLUTE_PFAFFIAN_H
#define MVMC_ABSOLUTE_PFAFFIAN_H

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

typedef enum {
  MVMC_PFAFFIAN_REGULAR = 0,
  MVMC_PFAFFIAN_NEAR_SINGULAR,
  MVMC_PFAFFIAN_SINGULAR,
  MVMC_PFAFFIAN_NONFINITE,
  MVMC_PFAFFIAN_INVALID
} MVMCPfaffianState;

typedef enum {
  MVMC_PFAFFIAN_STATUS_OK = 0,
  MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
  MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE,
  MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE,
  MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE
} MVMCPfaffianStatus;

typedef struct {
  MVMCPfaffianState state;
  int inverse_valid;
  int factor_info;
  int inverse_info;
  uint64_t rebuild_generation;
  double matrix_scale;
  double scaled_min_pivot;
  double inverse_residual;
  double complex pfaffian;
} MVMCAbsolutePfaffianResult;

typedef struct {
  int valid;
  double complex total;
  double sum_abs;
  double cancellation_ratio;
  size_t regular_count;
  size_t near_singular_count;
  size_t singular_count;
} MVMCProjectedAmplitudeResult;

/*
 * Evaluate an even-dimensional skew-symmetric matrix without modifying it.
 * Matrices use LAPACK column-major layout.  A non-positive tolerance selects
 * the scale-aware default.  inverse_out is modified only when inverse_valid
 * is 1; singular, near-singular, nonfinite, and error paths leave it intact.
 */
MVMCPfaffianStatus mvmc_absolute_pfaffian_real(
    const double *matrix, int n, int lda,
    double *inverse_out, int inverse_lda,
    uint64_t rebuild_generation,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCAbsolutePfaffianResult *result);

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex(
    const double complex *matrix, int n, int lda,
    double complex *inverse_out, int inverse_lda,
    uint64_t rebuild_generation,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCAbsolutePfaffianResult *result);

/*
 * Sum QPFullWeight[q] * PfM[q] without consulting inverse availability.
 * The result is usable only when the returned status is OK and result->valid
 * is 1.  Error paths leave a zeroed result with valid=0.
 */
MVMCPfaffianStatus mvmc_projected_amplitude(
    const MVMCAbsolutePfaffianResult *components,
    const double complex *weights, size_t component_count,
    MVMCProjectedAmplitudeResult *result);

/*
 * Evaluate a rank-local contiguous QP slice against a replicated global
 * weight array.  MPI reduction, if needed, is deliberately owned by the
 * caller; this function never reads rank-0-only globals.
 */
MVMCPfaffianStatus mvmc_projected_amplitude_slice(
    const MVMCAbsolutePfaffianResult *local_components,
    size_t local_component_count, const double complex *global_weights,
    int qp_total, int qp_start, int qp_end,
    MVMCProjectedAmplitudeResult *result);

/* False when this translation unit was compiled with fast-math semantics. */
int mvmc_absolute_pfaffian_strict_fp_enabled(void);

const char *mvmc_pfaffian_state_name(MVMCPfaffianState state);

#endif /* MVMC_ABSOLUTE_PFAFFIAN_H */
