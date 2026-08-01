/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "absolute_pfaffian.h"

#ifndef _BLAS_EXTERNS_H
#ifdef _mpi_use
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif
#include "blas_externs.h"
#endif

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* GCC 12 ICEs in vect_transform_reduction on this file at -O3. */
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ == 12
/* GCC 12 ICEs while vectorizing the compensated aggregation loop at -O3. */
#pragma GCC optimize("no-tree-loop-vectorize")
#endif

static void initialize_result(MVMCAbsolutePfaffianResult *result,
                              uint64_t generation) {
  result->state = MVMC_PFAFFIAN_INVALID;
  result->inverse_valid = 0;
  result->factor_info = 0;
  result->inverse_info = 0;
  result->rebuild_generation = generation;
  result->matrix_scale = NAN;
  result->scaled_min_pivot = NAN;
  result->inverse_residual = NAN;
  result->pfaffian = 0.0;
}

static int checked_square_size(int n, size_t element_size, size_t *count) {
  size_t dimension;

  if (n <= 0) return 0;
  dimension = (size_t)n;
  if (dimension > SIZE_MAX / dimension) return 0;
  *count = dimension * dimension;
  if (*count > SIZE_MAX / element_size) return 0;
  return 1;
}

static int real_matrix_properties(const double *matrix, int n, int lda,
                                  double *scale) {
  int column, row;
  double matrix_scale = 0.0;
  double skew_error = 0.0;

  for (column = 0; column < n; ++column) {
    for (row = 0; row < n; ++row) {
      const double value = matrix[row + (size_t)column * (size_t)lda];
      if (!isfinite(value)) return -1;
      matrix_scale = fmax(matrix_scale, fabs(value));
    }
  }
  for (column = 0; column < n; ++column) {
    skew_error = fmax(skew_error,
                      fabs(matrix[column + (size_t)column * (size_t)lda]));
    for (row = 0; row < column; ++row) {
      const double pair_sum =
          matrix[row + (size_t)column * (size_t)lda] +
          matrix[column + (size_t)row * (size_t)lda];
      skew_error = fmax(skew_error, fabs(pair_sum));
    }
  }
  *scale = matrix_scale;
  return skew_error <= 64.0 * DBL_EPSILON * fmax(DBL_MIN, matrix_scale);
}

static int complex_matrix_properties(const double complex *matrix, int n,
                                     int lda, double *scale) {
  int column, row;
  double matrix_scale = 0.0;
  double skew_error = 0.0;

  for (column = 0; column < n; ++column) {
    for (row = 0; row < n; ++row) {
      const double complex value =
          matrix[row + (size_t)column * (size_t)lda];
      if (!isfinite(creal(value)) || !isfinite(cimag(value))) return -1;
      matrix_scale = fmax(matrix_scale, cabs(value));
    }
  }
  for (column = 0; column < n; ++column) {
    skew_error = fmax(
        skew_error,
        cabs(matrix[column + (size_t)column * (size_t)lda]));
    for (row = 0; row < column; ++row) {
      const double complex pair_sum =
          matrix[row + (size_t)column * (size_t)lda] +
          matrix[column + (size_t)row * (size_t)lda];
      skew_error = fmax(skew_error, cabs(pair_sum));
    }
  }
  *scale = matrix_scale;
  return skew_error <= 64.0 * DBL_EPSILON * fmax(DBL_MIN, matrix_scale);
}

static double real_inverse_residual(const double *matrix, int lda,
                                    const double *inverse, int n) {
  int column, inner, row;
  double matrix_scale = 0.0;
  double inverse_scale = 0.0;
  double residual = 0.0;

  for (column = 0; column < n; ++column) {
    for (row = 0; row < n; ++row) {
      matrix_scale = fmax(
          matrix_scale,
          fabs(matrix[row + (size_t)column * (size_t)lda]));
      inverse_scale = fmax(
          inverse_scale, fabs(inverse[row + (size_t)column * (size_t)n]));
    }
  }
  for (column = 0; column < n; ++column) {
    for (row = 0; row < n; ++row) {
      double product = 0.0;
      for (inner = 0; inner < n; ++inner) {
        product += matrix[row + (size_t)inner * (size_t)lda] *
                   inverse[inner + (size_t)column * (size_t)n];
      }
      product -= (row == column) ? 1.0 : 0.0;
      residual = fmax(residual, fabs(product));
    }
  }
  return residual /
         fmax(1.0, (double)n * matrix_scale * inverse_scale);
}

static double complex_inverse_residual(const double complex *matrix, int lda,
                                       const double complex *inverse, int n) {
  int column, inner, row;
  double matrix_scale = 0.0;
  double inverse_scale = 0.0;
  double residual = 0.0;

  for (column = 0; column < n; ++column) {
    for (row = 0; row < n; ++row) {
      matrix_scale = fmax(
          matrix_scale,
          cabs(matrix[row + (size_t)column * (size_t)lda]));
      inverse_scale = fmax(
          inverse_scale, cabs(inverse[row + (size_t)column * (size_t)n]));
    }
  }
  for (column = 0; column < n; ++column) {
    for (row = 0; row < n; ++row) {
      double complex product = 0.0;
      for (inner = 0; inner < n; ++inner) {
        product += matrix[row + (size_t)inner * (size_t)lda] *
                   inverse[inner + (size_t)column * (size_t)n];
      }
      product -= (row == column) ? 1.0 : 0.0;
      residual = fmax(residual, cabs(product));
    }
  }
  return residual /
         fmax(1.0, (double)n * matrix_scale * inverse_scale);
}

static double effective_pivot_tolerance(double tolerance, int n) {
  if (tolerance > 0.0 && isfinite(tolerance)) return tolerance;
  return 64.0 * DBL_EPSILON * (double)n;
}

static double effective_residual_tolerance(double tolerance, int n) {
  if (tolerance > 0.0 && isfinite(tolerance)) return tolerance;
  return 256.0 * DBL_EPSILON * (double)n;
}

int mvmc_absolute_pfaffian_strict_fp_enabled(void) {
#ifdef __FAST_MATH__
  return 0;
#else
  return 1;
#endif
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_real(
    const double *matrix, int n, int lda,
    double *inverse_out, int inverse_lda,
    uint64_t rebuild_generation,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCAbsolutePfaffianResult *result) {
  double *factor = NULL;
  double *inverse_factor = NULL;
  double *inverse_work = NULL;
  double *factor_work = NULL;
  int *pivots = NULL;
  size_t matrix_count;
  int column, row, info = 0, lwork;
  double work_query = 0.0;
  double pfaffian = 0.0;
  double min_pivot = INFINITY;
  MVMCPfaffianStatus status = MVMC_PFAFFIAN_STATUS_OK;

  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  initialize_result(result, rebuild_generation);
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (matrix == NULL || inverse_out == NULL || n <= 0 || (n % 2) != 0 ||
      lda < n || inverse_lda < n ||
      !isfinite(scaled_pivot_tolerance) ||
      !isfinite(residual_tolerance) ||
      !checked_square_size(n, sizeof(*factor), &matrix_count)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  info = real_matrix_properties(matrix, n, lda, &result->matrix_scale);
  if (info < 0) {
    result->state = MVMC_PFAFFIAN_NONFINITE;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (info == 0) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;

  factor = (double *)malloc(matrix_count * sizeof(*factor));
  pivots = (int *)malloc((size_t)n * sizeof(*pivots));
  if (factor == NULL || pivots == NULL) {
    status = MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
    goto cleanup;
  }
  for (column = 0; column < n; ++column) {
    for (row = 0; row < n; ++row) {
      factor[row + (size_t)column * (size_t)n] =
          matrix[row + (size_t)column * (size_t)lda];
    }
  }

  lwork = -1;
  M_DSKTRF("U", "P", &n, factor, &n, pivots, &work_query, &lwork, &info);
  if (info != 0 || !isfinite(work_query) || work_query < 1.0 ||
      work_query > (double)INT_MAX) {
    status = MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
    goto cleanup;
  }
  lwork = (int)ceil(work_query);
  factor_work = (double *)malloc((size_t)lwork * sizeof(*factor_work));
  if (factor_work == NULL) {
    status = MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
    goto cleanup;
  }

  M_DSKTRF("U", "P", &n, factor, &n, pivots, factor_work, &lwork, &info);
  result->factor_info = info;
  if (info < 0) {
    status = MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
    goto cleanup;
  }
  if (info > 0) {
    result->state = MVMC_PFAFFIAN_SINGULAR;
    result->scaled_min_pivot = 0.0;
    result->pfaffian = 0.0;
    goto cleanup;
  }

  utu2pfa_d(n, factor, n, pivots, &pfaffian);
  if (!isfinite(pfaffian)) {
    result->state = MVMC_PFAFFIAN_NONFINITE;
    goto cleanup;
  }
  result->pfaffian = pfaffian;
  for (row = 0; row < n; row += 2) {
    min_pivot = fmin(min_pivot,
                     fabs(factor[row + (size_t)(row + 1) * (size_t)n]));
  }
  result->scaled_min_pivot =
      min_pivot / fmax(result->matrix_scale, DBL_MIN);

  if (!isfinite(result->scaled_min_pivot)) {
    result->state = MVMC_PFAFFIAN_NONFINITE;
    goto cleanup;
  }
  if (result->scaled_min_pivot <
      effective_pivot_tolerance(scaled_pivot_tolerance, n)) {
    result->state = MVMC_PFAFFIAN_NEAR_SINGULAR;
    goto cleanup;
  }

  free(factor_work);
  factor_work = NULL;
  inverse_factor =
      (double *)malloc(matrix_count * sizeof(*inverse_factor));
  if (inverse_factor == NULL) {
    status = MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
    goto cleanup;
  }
  for (column = 0; column < n; ++column) {
    for (row = 0; row < n; ++row) {
      inverse_factor[row + (size_t)column * (size_t)n] =
          matrix[row + (size_t)column * (size_t)lda];
    }
  }
  M_DGETRF(&n, &n, inverse_factor, &n, pivots, &info);
  result->inverse_info = info;
  if (info < 0) {
    status = MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
    goto cleanup;
  }
  if (info > 0) {
    result->state = MVMC_PFAFFIAN_NEAR_SINGULAR;
    goto cleanup;
  }
  lwork = -1;
  M_DGETRI(&n, inverse_factor, &n, pivots, &work_query, &lwork, &info);
  result->inverse_info = info;
  if (info < 0 || !isfinite(work_query) || work_query < 1.0 ||
      work_query > (double)INT_MAX) {
    status = MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
    goto cleanup;
  }
  if (info > 0) {
    result->state = MVMC_PFAFFIAN_NEAR_SINGULAR;
    goto cleanup;
  }
  lwork = (int)ceil(work_query);
  inverse_work = (double *)malloc((size_t)lwork * sizeof(*inverse_work));
  if (inverse_work == NULL) {
    status = MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
    goto cleanup;
  }
  M_DGETRI(&n, inverse_factor, &n, pivots, inverse_work, &lwork, &info);
  result->inverse_info = info;
  if (info < 0) {
    status = MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
    goto cleanup;
  }
  if (info > 0) {
    result->state = MVMC_PFAFFIAN_NEAR_SINGULAR;
    goto cleanup;
  }
  result->inverse_residual =
      real_inverse_residual(matrix, lda, inverse_factor, n);
  if (!isfinite(result->inverse_residual)) {
    result->state = MVMC_PFAFFIAN_NONFINITE;
    goto cleanup;
  }
  if (result->inverse_residual >
      effective_residual_tolerance(residual_tolerance, n)) {
    result->state = MVMC_PFAFFIAN_NEAR_SINGULAR;
    goto cleanup;
  }

  for (column = 0; column < n; ++column) {
    for (row = 0; row < n; ++row) {
      inverse_out[row + (size_t)column * (size_t)inverse_lda] =
          inverse_factor[row + (size_t)column * (size_t)n];
    }
  }
  result->state = MVMC_PFAFFIAN_REGULAR;
  result->inverse_valid = 1;

cleanup:
  free(factor_work);
  free(pivots);
  free(inverse_work);
  free(inverse_factor);
  free(factor);
  return status;
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex(
    const double complex *matrix, int n, int lda,
    double complex *inverse_out, int inverse_lda,
    uint64_t rebuild_generation,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCAbsolutePfaffianResult *result) {
  double complex *factor = NULL;
  double complex *inverse_factor = NULL;
  double complex *inverse_work = NULL;
  double complex *factor_work = NULL;
  int *pivots = NULL;
  size_t matrix_count;
  int column, row, info = 0, lwork;
  double complex work_query = 0.0;
  double complex pfaffian = 0.0;
  double min_pivot = INFINITY;
  MVMCPfaffianStatus status = MVMC_PFAFFIAN_STATUS_OK;

  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  initialize_result(result, rebuild_generation);
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (matrix == NULL || inverse_out == NULL || n <= 0 || (n % 2) != 0 ||
      lda < n || inverse_lda < n ||
      !isfinite(scaled_pivot_tolerance) ||
      !isfinite(residual_tolerance) ||
      !checked_square_size(n, sizeof(*factor), &matrix_count)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  info = complex_matrix_properties(matrix, n, lda, &result->matrix_scale);
  if (info < 0) {
    result->state = MVMC_PFAFFIAN_NONFINITE;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (info == 0) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;

  factor = (double complex *)malloc(matrix_count * sizeof(*factor));
  pivots = (int *)malloc((size_t)n * sizeof(*pivots));
  if (factor == NULL || pivots == NULL) {
    status = MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
    goto cleanup;
  }
  for (column = 0; column < n; ++column) {
    for (row = 0; row < n; ++row) {
      factor[row + (size_t)column * (size_t)n] =
          matrix[row + (size_t)column * (size_t)lda];
    }
  }

  lwork = -1;
  M_ZSKTRF("U", "P", &n, factor, &n, pivots, &work_query, &lwork, &info);
  if (info != 0 || !isfinite(creal(work_query)) || cimag(work_query) != 0.0 ||
      creal(work_query) < 1.0 || creal(work_query) > (double)INT_MAX) {
    status = MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
    goto cleanup;
  }
  lwork = (int)ceil(creal(work_query));
  factor_work =
      (double complex *)malloc((size_t)lwork * sizeof(*factor_work));
  if (factor_work == NULL) {
    status = MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
    goto cleanup;
  }

  M_ZSKTRF("U", "P", &n, factor, &n, pivots, factor_work, &lwork, &info);
  result->factor_info = info;
  if (info < 0) {
    status = MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
    goto cleanup;
  }
  if (info > 0) {
    result->state = MVMC_PFAFFIAN_SINGULAR;
    result->scaled_min_pivot = 0.0;
    result->pfaffian = 0.0;
    goto cleanup;
  }

  utu2pfa_z(n, factor, n, pivots, &pfaffian);
  if (!isfinite(creal(pfaffian)) || !isfinite(cimag(pfaffian))) {
    result->state = MVMC_PFAFFIAN_NONFINITE;
    goto cleanup;
  }
  result->pfaffian = pfaffian;
  for (row = 0; row < n; row += 2) {
    min_pivot = fmin(
        min_pivot,
        cabs(factor[row + (size_t)(row + 1) * (size_t)n]));
  }
  result->scaled_min_pivot =
      min_pivot / fmax(result->matrix_scale, DBL_MIN);

  if (!isfinite(result->scaled_min_pivot)) {
    result->state = MVMC_PFAFFIAN_NONFINITE;
    goto cleanup;
  }
  if (result->scaled_min_pivot <
      effective_pivot_tolerance(scaled_pivot_tolerance, n)) {
    result->state = MVMC_PFAFFIAN_NEAR_SINGULAR;
    goto cleanup;
  }

  free(factor_work);
  factor_work = NULL;
  inverse_factor =
      (double complex *)malloc(matrix_count * sizeof(*inverse_factor));
  if (inverse_factor == NULL) {
    status = MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
    goto cleanup;
  }
  for (column = 0; column < n; ++column) {
    for (row = 0; row < n; ++row) {
      inverse_factor[row + (size_t)column * (size_t)n] =
          matrix[row + (size_t)column * (size_t)lda];
    }
  }
  M_ZGETRF(&n, &n, inverse_factor, &n, pivots, &info);
  result->inverse_info = info;
  if (info < 0) {
    status = MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
    goto cleanup;
  }
  if (info > 0) {
    result->state = MVMC_PFAFFIAN_NEAR_SINGULAR;
    goto cleanup;
  }
  lwork = -1;
  M_ZGETRI(&n, inverse_factor, &n, pivots, &work_query, &lwork, &info);
  result->inverse_info = info;
  if (info < 0 || !isfinite(creal(work_query)) ||
      cimag(work_query) != 0.0 || creal(work_query) < 1.0 ||
      creal(work_query) > (double)INT_MAX) {
    status = MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
    goto cleanup;
  }
  if (info > 0) {
    result->state = MVMC_PFAFFIAN_NEAR_SINGULAR;
    goto cleanup;
  }
  lwork = (int)ceil(creal(work_query));
  inverse_work =
      (double complex *)malloc((size_t)lwork * sizeof(*inverse_work));
  if (inverse_work == NULL) {
    status = MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
    goto cleanup;
  }
  M_ZGETRI(&n, inverse_factor, &n, pivots, inverse_work, &lwork, &info);
  result->inverse_info = info;
  if (info < 0) {
    status = MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
    goto cleanup;
  }
  if (info > 0) {
    result->state = MVMC_PFAFFIAN_NEAR_SINGULAR;
    goto cleanup;
  }
  result->inverse_residual =
      complex_inverse_residual(matrix, lda, inverse_factor, n);
  if (!isfinite(result->inverse_residual)) {
    result->state = MVMC_PFAFFIAN_NONFINITE;
    goto cleanup;
  }
  if (result->inverse_residual >
      effective_residual_tolerance(residual_tolerance, n)) {
    result->state = MVMC_PFAFFIAN_NEAR_SINGULAR;
    goto cleanup;
  }

  for (column = 0; column < n; ++column) {
    for (row = 0; row < n; ++row) {
      inverse_out[row + (size_t)column * (size_t)inverse_lda] =
          inverse_factor[row + (size_t)column * (size_t)n];
    }
  }
  result->state = MVMC_PFAFFIAN_REGULAR;
  result->inverse_valid = 1;

cleanup:
  free(factor_work);
  free(pivots);
  free(inverse_work);
  free(inverse_factor);
  free(factor);
  return status;
}

MVMCPfaffianStatus mvmc_projected_amplitude(
    const MVMCAbsolutePfaffianResult *components,
    const double complex *weights, size_t component_count,
    MVMCProjectedAmplitudeResult *result) {
  size_t index;
  double real_sum = 0.0, real_compensation = 0.0;
  double imag_sum = 0.0, imag_compensation = 0.0;
  double abs_sum = 0.0, abs_compensation = 0.0;
  MVMCProjectedAmplitudeResult accumulated;

  if (result == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  memset(result, 0, sizeof(*result));
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (components == NULL || weights == NULL || component_count == 0) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  memset(&accumulated, 0, sizeof(accumulated));
  for (index = 0; index < component_count; ++index) {
    double complex term;
    double value, corrected, updated;

    if (!isfinite(creal(weights[index])) ||
        !isfinite(cimag(weights[index])) ||
        !isfinite(creal(components[index].pfaffian)) ||
        !isfinite(cimag(components[index].pfaffian)) ||
        components[index].state == MVMC_PFAFFIAN_NONFINITE ||
        components[index].state == MVMC_PFAFFIAN_INVALID) {
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    switch (components[index].state) {
      case MVMC_PFAFFIAN_REGULAR:
        ++accumulated.regular_count;
        break;
      case MVMC_PFAFFIAN_NEAR_SINGULAR:
        ++accumulated.near_singular_count;
        break;
      case MVMC_PFAFFIAN_SINGULAR:
        if (components[index].pfaffian != 0.0) {
          return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
        }
        ++accumulated.singular_count;
        break;
      default:
        return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }

    term = weights[index] * components[index].pfaffian;
    if (!isfinite(creal(term)) || !isfinite(cimag(term))) {
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    value = cabs(term);
    corrected = value - abs_compensation;
    updated = abs_sum + corrected;
    abs_compensation = (updated - abs_sum) - corrected;
    abs_sum = updated;

    value = creal(term);
    corrected = value - real_compensation;
    updated = real_sum + corrected;
    real_compensation = (updated - real_sum) - corrected;
    real_sum = updated;

    value = cimag(term);
    corrected = value - imag_compensation;
    updated = imag_sum + corrected;
    imag_compensation = (updated - imag_sum) - corrected;
    imag_sum = updated;
  }
  if (!isfinite(abs_sum)) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  accumulated.total = real_sum + I * imag_sum;
  accumulated.sum_abs = abs_sum;
  accumulated.cancellation_ratio =
      accumulated.sum_abs == 0.0
          ? 0.0
          : fmin(1.0, cabs(accumulated.total) / accumulated.sum_abs);
  accumulated.valid = 1;
  *result = accumulated;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_projected_amplitude_slice(
    const MVMCAbsolutePfaffianResult *local_components,
    size_t local_component_count, const double complex *global_weights,
    int qp_total, int qp_start, int qp_end,
    MVMCProjectedAmplitudeResult *result) {
  if (result != NULL) memset(result, 0, sizeof(*result));
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (qp_total <= 0 || qp_start < 0 || qp_end <= qp_start ||
      qp_end > qp_total ||
      (size_t)(qp_end - qp_start) != local_component_count) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  return mvmc_projected_amplitude(
      local_components, global_weights == NULL ? NULL : global_weights + qp_start,
      local_component_count, result);
}

const char *mvmc_pfaffian_state_name(MVMCPfaffianState state) {
  switch (state) {
    case MVMC_PFAFFIAN_REGULAR:
      return "REGULAR";
    case MVMC_PFAFFIAN_NEAR_SINGULAR:
      return "NEAR_SINGULAR";
    case MVMC_PFAFFIAN_SINGULAR:
      return "SINGULAR";
    case MVMC_PFAFFIAN_NONFINITE:
      return "NONFINITE";
    case MVMC_PFAFFIAN_INVALID:
      return "INVALID";
  }
  return "INVALID";
}
