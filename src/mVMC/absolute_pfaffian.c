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

struct MVMCAbsolutePfaffianRealWorkspace {
  int n;
  size_t matrix_count;
  int factor_lwork;
  int inverse_lwork;
  int *pivots;
  double *factor;
  double *inverse_factor;
  double *factor_work;
  double *inverse_work;
};

struct MVMCAbsolutePfaffianComplexWorkspace {
  int n;
  size_t matrix_count;
  int factor_lwork;
  int inverse_lwork;
  int *pivots;
  double complex *factor;
  double complex *inverse_factor;
  double complex *factor_work;
  double complex *inverse_work;
};

struct MVMCAbsolutePfaffianRealValueWorkspace {
  int n;
  size_t matrix_count;
  int factor_lwork;
  int *pivots;
  double *factor;
  double *factor_work;
};

struct MVMCAbsolutePfaffianComplexValueWorkspace {
  int n;
  size_t matrix_count;
  int factor_lwork;
  int *pivots;
  double complex *factor;
  double complex *factor_work;
};

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

static void initialize_value_result(MVMCAbsolutePfaffianValueResult *result) {
  result->state = MVMC_PFAFFIAN_VALUE_INVALID;
  result->factor_info = 0;
  result->matrix_scale = NAN;
  result->scaled_min_pivot = NAN;
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

static int valid_real_lwork_query(double query, int *lwork) {
  if (!isfinite(query) || query < 1.0 || query > (double)INT_MAX) return 0;
  *lwork = (int)ceil(query);
  return 1;
}

static int valid_complex_lwork_query(double complex query, int *lwork) {
  if (!isfinite(creal(query)) || cimag(query) != 0.0 ||
      creal(query) < 1.0 || creal(query) > (double)INT_MAX) {
    return 0;
  }
  *lwork = (int)ceil(creal(query));
  return 1;
}

void mvmc_absolute_pfaffian_real_workspace_destroy(
    MVMCAbsolutePfaffianRealWorkspace *workspace) {
  if (workspace == NULL) return;
  free(workspace->inverse_work);
  free(workspace->factor_work);
  free(workspace->inverse_factor);
  free(workspace->factor);
  free(workspace->pivots);
  free(workspace);
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_real_workspace_create(
    int n, MVMCAbsolutePfaffianRealWorkspace **workspace) {
  MVMCAbsolutePfaffianRealWorkspace *created = NULL;
  double factor_query = 0.0;
  double inverse_query = 0.0;
  int info = 0;
  int query_lwork = -1;
  int index;
  size_t matrix_count;

  if (workspace == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (n <= 0 || (n % 2) != 0 ||
      !checked_square_size(n, sizeof(double), &matrix_count) ||
      (size_t)n > SIZE_MAX / sizeof(int)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  created = (MVMCAbsolutePfaffianRealWorkspace *)calloc(1, sizeof(*created));
  if (created == NULL) return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  created->n = n;
  created->matrix_count = matrix_count;
  created->pivots = (int *)malloc((size_t)n * sizeof(*created->pivots));
  created->factor =
      (double *)calloc(created->matrix_count, sizeof(*created->factor));
  created->inverse_factor = (double *)calloc(
      created->matrix_count, sizeof(*created->inverse_factor));
  if (created->pivots == NULL || created->factor == NULL ||
      created->inverse_factor == NULL) {
    mvmc_absolute_pfaffian_real_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  }
  for (index = 0; index < n; ++index) created->pivots[index] = index + 1;

  M_DSKTRF("U", "P", &n, created->factor, &n, created->pivots,
           &factor_query, &query_lwork, &info);
  if (info != 0 ||
      !valid_real_lwork_query(factor_query, &created->factor_lwork) ||
      (size_t)created->factor_lwork >
          SIZE_MAX / sizeof(*created->factor_work)) {
    mvmc_absolute_pfaffian_real_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
  }
  M_DGETRI(&n, created->inverse_factor, &n, created->pivots,
           &inverse_query, &query_lwork, &info);
  if (info != 0 ||
      !valid_real_lwork_query(inverse_query, &created->inverse_lwork) ||
      (size_t)created->inverse_lwork >
          SIZE_MAX / sizeof(*created->inverse_work)) {
    mvmc_absolute_pfaffian_real_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
  }
  created->factor_work = (double *)malloc(
      (size_t)created->factor_lwork * sizeof(*created->factor_work));
  created->inverse_work = (double *)malloc(
      (size_t)created->inverse_lwork * sizeof(*created->inverse_work));
  if (created->factor_work == NULL || created->inverse_work == NULL) {
    mvmc_absolute_pfaffian_real_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  }
  *workspace = created;
  return MVMC_PFAFFIAN_STATUS_OK;
}

void mvmc_absolute_pfaffian_complex_workspace_destroy(
    MVMCAbsolutePfaffianComplexWorkspace *workspace) {
  if (workspace == NULL) return;
  free(workspace->inverse_work);
  free(workspace->factor_work);
  free(workspace->inverse_factor);
  free(workspace->factor);
  free(workspace->pivots);
  free(workspace);
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex_workspace_create(
    int n, MVMCAbsolutePfaffianComplexWorkspace **workspace) {
  MVMCAbsolutePfaffianComplexWorkspace *created = NULL;
  double complex factor_query = 0.0;
  double complex inverse_query = 0.0;
  int info = 0;
  int query_lwork = -1;
  int index;
  size_t matrix_count;

  if (workspace == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (n <= 0 || (n % 2) != 0 ||
      !checked_square_size(n, sizeof(double complex), &matrix_count) ||
      (size_t)n > SIZE_MAX / sizeof(int)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  created =
      (MVMCAbsolutePfaffianComplexWorkspace *)calloc(1, sizeof(*created));
  if (created == NULL) return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  created->n = n;
  created->matrix_count = matrix_count;
  created->pivots = (int *)malloc((size_t)n * sizeof(*created->pivots));
  created->factor = (double complex *)calloc(
      created->matrix_count, sizeof(*created->factor));
  created->inverse_factor = (double complex *)calloc(
      created->matrix_count, sizeof(*created->inverse_factor));
  if (created->pivots == NULL || created->factor == NULL ||
      created->inverse_factor == NULL) {
    mvmc_absolute_pfaffian_complex_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  }
  for (index = 0; index < n; ++index) created->pivots[index] = index + 1;

  M_ZSKTRF("U", "P", &n, created->factor, &n, created->pivots,
           &factor_query, &query_lwork, &info);
  if (info != 0 ||
      !valid_complex_lwork_query(factor_query, &created->factor_lwork) ||
      (size_t)created->factor_lwork >
          SIZE_MAX / sizeof(*created->factor_work)) {
    mvmc_absolute_pfaffian_complex_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
  }
  M_ZGETRI(&n, created->inverse_factor, &n, created->pivots,
           &inverse_query, &query_lwork, &info);
  if (info != 0 ||
      !valid_complex_lwork_query(inverse_query, &created->inverse_lwork) ||
      (size_t)created->inverse_lwork >
          SIZE_MAX / sizeof(*created->inverse_work)) {
    mvmc_absolute_pfaffian_complex_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
  }
  created->factor_work = (double complex *)malloc(
      (size_t)created->factor_lwork * sizeof(*created->factor_work));
  created->inverse_work = (double complex *)malloc(
      (size_t)created->inverse_lwork * sizeof(*created->inverse_work));
  if (created->factor_work == NULL || created->inverse_work == NULL) {
    mvmc_absolute_pfaffian_complex_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  }
  *workspace = created;
  return MVMC_PFAFFIAN_STATUS_OK;
}

void mvmc_absolute_pfaffian_real_value_workspace_destroy(
    MVMCAbsolutePfaffianRealValueWorkspace *workspace) {
  if (workspace == NULL) return;
  free(workspace->factor_work);
  free(workspace->factor);
  free(workspace->pivots);
  free(workspace);
}

size_t mvmc_absolute_pfaffian_real_value_workspace_bytes(
    const MVMCAbsolutePfaffianRealValueWorkspace *workspace) {
  if (workspace == NULL) return 0;
  return sizeof(*workspace) +
         (size_t)workspace->n * sizeof(*workspace->pivots) +
         workspace->matrix_count * sizeof(*workspace->factor) +
         (size_t)workspace->factor_lwork * sizeof(*workspace->factor_work);
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_real_value_workspace_create(
    int n, MVMCAbsolutePfaffianRealValueWorkspace **workspace) {
  MVMCAbsolutePfaffianRealValueWorkspace *created = NULL;
  double factor_query = 0.0;
  int info = 0;
  int query_lwork = -1;
  int index;
  size_t matrix_count;

  if (workspace == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (n <= 0 || (n % 2) != 0 ||
      !checked_square_size(n, sizeof(double), &matrix_count) ||
      (size_t)n > SIZE_MAX / sizeof(int)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  created = (MVMCAbsolutePfaffianRealValueWorkspace *)calloc(
      1, sizeof(*created));
  if (created == NULL) return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  created->n = n;
  created->matrix_count = matrix_count;
  created->pivots = (int *)malloc((size_t)n * sizeof(*created->pivots));
  created->factor =
      (double *)calloc(created->matrix_count, sizeof(*created->factor));
  if (created->pivots == NULL || created->factor == NULL) {
    mvmc_absolute_pfaffian_real_value_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  }
  for (index = 0; index < n; ++index) created->pivots[index] = index + 1;
  M_DSKTRF("U", "P", &n, created->factor, &n, created->pivots,
           &factor_query, &query_lwork, &info);
  if (info != 0 ||
      !valid_real_lwork_query(factor_query, &created->factor_lwork) ||
      (size_t)created->factor_lwork >
          SIZE_MAX / sizeof(*created->factor_work)) {
    mvmc_absolute_pfaffian_real_value_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
  }
  created->factor_work = (double *)malloc(
      (size_t)created->factor_lwork * sizeof(*created->factor_work));
  if (created->factor_work == NULL) {
    mvmc_absolute_pfaffian_real_value_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  }
  *workspace = created;
  return MVMC_PFAFFIAN_STATUS_OK;
}

void mvmc_absolute_pfaffian_complex_value_workspace_destroy(
    MVMCAbsolutePfaffianComplexValueWorkspace *workspace) {
  if (workspace == NULL) return;
  free(workspace->factor_work);
  free(workspace->factor);
  free(workspace->pivots);
  free(workspace);
}

size_t mvmc_absolute_pfaffian_complex_value_workspace_bytes(
    const MVMCAbsolutePfaffianComplexValueWorkspace *workspace) {
  if (workspace == NULL) return 0;
  return sizeof(*workspace) +
         (size_t)workspace->n * sizeof(*workspace->pivots) +
         workspace->matrix_count * sizeof(*workspace->factor) +
         (size_t)workspace->factor_lwork * sizeof(*workspace->factor_work);
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex_value_workspace_create(
    int n, MVMCAbsolutePfaffianComplexValueWorkspace **workspace) {
  MVMCAbsolutePfaffianComplexValueWorkspace *created = NULL;
  double complex factor_query = 0.0;
  int info = 0;
  int query_lwork = -1;
  int index;
  size_t matrix_count;

  if (workspace == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (n <= 0 || (n % 2) != 0 ||
      !checked_square_size(n, sizeof(double complex), &matrix_count) ||
      (size_t)n > SIZE_MAX / sizeof(int)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  created = (MVMCAbsolutePfaffianComplexValueWorkspace *)calloc(
      1, sizeof(*created));
  if (created == NULL) return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  created->n = n;
  created->matrix_count = matrix_count;
  created->pivots = (int *)malloc((size_t)n * sizeof(*created->pivots));
  created->factor = (double complex *)calloc(
      created->matrix_count, sizeof(*created->factor));
  if (created->pivots == NULL || created->factor == NULL) {
    mvmc_absolute_pfaffian_complex_value_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  }
  for (index = 0; index < n; ++index) created->pivots[index] = index + 1;
  M_ZSKTRF("U", "P", &n, created->factor, &n, created->pivots,
           &factor_query, &query_lwork, &info);
  if (info != 0 ||
      !valid_complex_lwork_query(factor_query, &created->factor_lwork) ||
      (size_t)created->factor_lwork >
          SIZE_MAX / sizeof(*created->factor_work)) {
    mvmc_absolute_pfaffian_complex_value_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
  }
  created->factor_work = (double complex *)malloc(
      (size_t)created->factor_lwork * sizeof(*created->factor_work));
  if (created->factor_work == NULL) {
    mvmc_absolute_pfaffian_complex_value_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  }
  *workspace = created;
  return MVMC_PFAFFIAN_STATUS_OK;
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

static MVMCPfaffianStatus factorize_real_value(
    const double *matrix, int n, int lda, double *factor, int *pivots,
    double *factor_work, int factor_lwork, double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianValueResult *result) {
  double pfaffian = 0.0;
  double min_pivot = INFINITY;
  int column, row, info = 0;

  info = real_matrix_properties(matrix, n, lda, &result->matrix_scale);
  if (info < 0) {
    result->state = MVMC_PFAFFIAN_VALUE_NONFINITE;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (info == 0) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  for (column = 0; column < n; ++column) {
    for (row = 0; row < n; ++row) {
      factor[row + (size_t)column * (size_t)n] =
          matrix[row + (size_t)column * (size_t)lda];
    }
  }
  M_DSKTRF("U", "P", &n, factor, &n, pivots, factor_work,
           &factor_lwork, &info);
  result->factor_info = info;
  if (info < 0) return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
  if (info > 0) {
    result->state = MVMC_PFAFFIAN_VALUE_SINGULAR;
    result->scaled_min_pivot = 0.0;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  utu2pfa_d(n, factor, n, pivots, &pfaffian);
  if (!isfinite(pfaffian)) {
    result->state = MVMC_PFAFFIAN_VALUE_NONFINITE;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  result->pfaffian = pfaffian;
  for (row = 0; row < n; row += 2) {
    min_pivot = fmin(
        min_pivot, fabs(factor[row + (size_t)(row + 1) * (size_t)n]));
  }
  result->scaled_min_pivot =
      min_pivot / fmax(result->matrix_scale, DBL_MIN);
  if (!isfinite(result->scaled_min_pivot)) {
    result->state = MVMC_PFAFFIAN_VALUE_NONFINITE;
  } else if (result->scaled_min_pivot <
             effective_pivot_tolerance(scaled_pivot_tolerance, n)) {
    result->state = MVMC_PFAFFIAN_VALUE_NEAR_PIVOT;
  } else {
    result->state = MVMC_PFAFFIAN_VALUE_WELL_PIVOTED;
  }
  return MVMC_PFAFFIAN_STATUS_OK;
}

static MVMCPfaffianStatus factorize_complex_value(
    const double complex *matrix, int n, int lda, double complex *factor,
    int *pivots, double complex *factor_work, int factor_lwork,
    double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianValueResult *result) {
  double complex pfaffian = 0.0;
  double min_pivot = INFINITY;
  int column, row, info = 0;

  info = complex_matrix_properties(matrix, n, lda, &result->matrix_scale);
  if (info < 0) {
    result->state = MVMC_PFAFFIAN_VALUE_NONFINITE;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (info == 0) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  for (column = 0; column < n; ++column) {
    for (row = 0; row < n; ++row) {
      factor[row + (size_t)column * (size_t)n] =
          matrix[row + (size_t)column * (size_t)lda];
    }
  }
  M_ZSKTRF("U", "P", &n, factor, &n, pivots, factor_work,
           &factor_lwork, &info);
  result->factor_info = info;
  if (info < 0) return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
  if (info > 0) {
    result->state = MVMC_PFAFFIAN_VALUE_SINGULAR;
    result->scaled_min_pivot = 0.0;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  utu2pfa_z(n, factor, n, pivots, &pfaffian);
  if (!isfinite(creal(pfaffian)) || !isfinite(cimag(pfaffian))) {
    result->state = MVMC_PFAFFIAN_VALUE_NONFINITE;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  result->pfaffian = pfaffian;
  for (row = 0; row < n; row += 2) {
    min_pivot = fmin(
        min_pivot, cabs(factor[row + (size_t)(row + 1) * (size_t)n]));
  }
  result->scaled_min_pivot =
      min_pivot / fmax(result->matrix_scale, DBL_MIN);
  if (!isfinite(result->scaled_min_pivot)) {
    result->state = MVMC_PFAFFIAN_VALUE_NONFINITE;
  } else if (result->scaled_min_pivot <
             effective_pivot_tolerance(scaled_pivot_tolerance, n)) {
    result->state = MVMC_PFAFFIAN_VALUE_NEAR_PIVOT;
  } else {
    result->state = MVMC_PFAFFIAN_VALUE_WELL_PIVOTED;
  }
  return MVMC_PFAFFIAN_STATUS_OK;
}

int mvmc_absolute_pfaffian_strict_fp_enabled(void) {
#ifdef __FAST_MATH__
  return 0;
#else
  return 1;
#endif
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_real_with_workspace(
    MVMCAbsolutePfaffianRealWorkspace *workspace,
    const double *matrix, int n, int lda,
    double *inverse_out, int inverse_lda,
    uint64_t rebuild_generation,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCAbsolutePfaffianResult *result) {
  double *factor = workspace == NULL ? NULL : workspace->factor;
  double *inverse_factor =
      workspace == NULL ? NULL : workspace->inverse_factor;
  double *inverse_work = workspace == NULL ? NULL : workspace->inverse_work;
  double *factor_work = workspace == NULL ? NULL : workspace->factor_work;
  int *pivots = workspace == NULL ? NULL : workspace->pivots;
  size_t matrix_count;
  int column, row, info = 0, lwork;
  MVMCAbsolutePfaffianValueResult value_result;
  MVMCPfaffianStatus status = MVMC_PFAFFIAN_STATUS_OK;

  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  initialize_result(result, rebuild_generation);
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (workspace == NULL || workspace->n != n || matrix == NULL ||
      inverse_out == NULL || n <= 0 || (n % 2) != 0 ||
      lda < n || inverse_lda < n ||
      !isfinite(scaled_pivot_tolerance) ||
      !isfinite(residual_tolerance) ||
      !checked_square_size(n, sizeof(*factor), &matrix_count) ||
      workspace->matrix_count != matrix_count) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  initialize_value_result(&value_result);
  status = factorize_real_value(
      matrix, n, lda, factor, pivots, factor_work,
      workspace->factor_lwork, scaled_pivot_tolerance, &value_result);
  result->factor_info = value_result.factor_info;
  result->matrix_scale = value_result.matrix_scale;
  result->scaled_min_pivot = value_result.scaled_min_pivot;
  result->pfaffian = value_result.pfaffian;
  if (status != MVMC_PFAFFIAN_STATUS_OK) goto cleanup;
  if (value_result.state == MVMC_PFAFFIAN_VALUE_NONFINITE) {
    result->state = MVMC_PFAFFIAN_NONFINITE;
    goto cleanup;
  }
  if (value_result.state == MVMC_PFAFFIAN_VALUE_SINGULAR) {
    result->state = MVMC_PFAFFIAN_SINGULAR;
    goto cleanup;
  }
  if (value_result.state == MVMC_PFAFFIAN_VALUE_NEAR_PIVOT) {
    result->state = MVMC_PFAFFIAN_NEAR_SINGULAR;
    goto cleanup;
  }
  if (value_result.state != MVMC_PFAFFIAN_VALUE_WELL_PIVOTED) {
    status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
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
  lwork = workspace->inverse_lwork;
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
  return status;
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_real(
    const double *matrix, int n, int lda,
    double *inverse_out, int inverse_lda,
    uint64_t rebuild_generation,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCAbsolutePfaffianResult *result) {
  MVMCAbsolutePfaffianRealWorkspace *workspace = NULL;
  MVMCPfaffianStatus status;
  size_t matrix_count;

  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  initialize_result(result, rebuild_generation);
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (matrix == NULL || inverse_out == NULL || n <= 0 || (n % 2) != 0 ||
      lda < n || inverse_lda < n ||
      !isfinite(scaled_pivot_tolerance) ||
      !isfinite(residual_tolerance) ||
      !checked_square_size(n, sizeof(*matrix), &matrix_count)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_absolute_pfaffian_real_workspace_create(n, &workspace);
  if (status == MVMC_PFAFFIAN_STATUS_OK) {
    status = mvmc_absolute_pfaffian_real_with_workspace(
        workspace, matrix, n, lda, inverse_out, inverse_lda,
        rebuild_generation, scaled_pivot_tolerance, residual_tolerance,
        result);
  }
  mvmc_absolute_pfaffian_real_workspace_destroy(workspace);
  return status;
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex_with_workspace(
    MVMCAbsolutePfaffianComplexWorkspace *workspace,
    const double complex *matrix, int n, int lda,
    double complex *inverse_out, int inverse_lda,
    uint64_t rebuild_generation,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCAbsolutePfaffianResult *result) {
  double complex *factor = workspace == NULL ? NULL : workspace->factor;
  double complex *inverse_factor =
      workspace == NULL ? NULL : workspace->inverse_factor;
  double complex *inverse_work =
      workspace == NULL ? NULL : workspace->inverse_work;
  double complex *factor_work =
      workspace == NULL ? NULL : workspace->factor_work;
  int *pivots = workspace == NULL ? NULL : workspace->pivots;
  size_t matrix_count;
  int column, row, info = 0, lwork;
  MVMCAbsolutePfaffianValueResult value_result;
  MVMCPfaffianStatus status = MVMC_PFAFFIAN_STATUS_OK;

  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  initialize_result(result, rebuild_generation);
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (workspace == NULL || workspace->n != n || matrix == NULL ||
      inverse_out == NULL || n <= 0 || (n % 2) != 0 ||
      lda < n || inverse_lda < n ||
      !isfinite(scaled_pivot_tolerance) ||
      !isfinite(residual_tolerance) ||
      !checked_square_size(n, sizeof(*factor), &matrix_count) ||
      workspace->matrix_count != matrix_count) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  initialize_value_result(&value_result);
  status = factorize_complex_value(
      matrix, n, lda, factor, pivots, factor_work,
      workspace->factor_lwork, scaled_pivot_tolerance, &value_result);
  result->factor_info = value_result.factor_info;
  result->matrix_scale = value_result.matrix_scale;
  result->scaled_min_pivot = value_result.scaled_min_pivot;
  result->pfaffian = value_result.pfaffian;
  if (status != MVMC_PFAFFIAN_STATUS_OK) goto cleanup;
  if (value_result.state == MVMC_PFAFFIAN_VALUE_NONFINITE) {
    result->state = MVMC_PFAFFIAN_NONFINITE;
    goto cleanup;
  }
  if (value_result.state == MVMC_PFAFFIAN_VALUE_SINGULAR) {
    result->state = MVMC_PFAFFIAN_SINGULAR;
    goto cleanup;
  }
  if (value_result.state == MVMC_PFAFFIAN_VALUE_NEAR_PIVOT) {
    result->state = MVMC_PFAFFIAN_NEAR_SINGULAR;
    goto cleanup;
  }
  if (value_result.state != MVMC_PFAFFIAN_VALUE_WELL_PIVOTED) {
    status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
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
  lwork = workspace->inverse_lwork;
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
  return status;
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex(
    const double complex *matrix, int n, int lda,
    double complex *inverse_out, int inverse_lda,
    uint64_t rebuild_generation,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCAbsolutePfaffianResult *result) {
  MVMCAbsolutePfaffianComplexWorkspace *workspace = NULL;
  MVMCPfaffianStatus status;
  size_t matrix_count;

  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  initialize_result(result, rebuild_generation);
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (matrix == NULL || inverse_out == NULL || n <= 0 || (n % 2) != 0 ||
      lda < n || inverse_lda < n ||
      !isfinite(scaled_pivot_tolerance) ||
      !isfinite(residual_tolerance) ||
      !checked_square_size(n, sizeof(*matrix), &matrix_count)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_absolute_pfaffian_complex_workspace_create(n, &workspace);
  if (status == MVMC_PFAFFIAN_STATUS_OK) {
    status = mvmc_absolute_pfaffian_complex_with_workspace(
        workspace, matrix, n, lda, inverse_out, inverse_lda,
        rebuild_generation, scaled_pivot_tolerance, residual_tolerance,
        result);
  }
  mvmc_absolute_pfaffian_complex_workspace_destroy(workspace);
  return status;
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_real_value_with_workspace(
    MVMCAbsolutePfaffianRealValueWorkspace *workspace,
    const double *matrix, int n, int lda, double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianValueResult *result) {
  size_t matrix_count;

  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  initialize_value_result(result);
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (workspace == NULL || workspace->n != n || matrix == NULL ||
      n <= 0 || (n % 2) != 0 || lda < n ||
      !isfinite(scaled_pivot_tolerance) ||
      !checked_square_size(n, sizeof(*matrix), &matrix_count) ||
      workspace->matrix_count != matrix_count) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  return factorize_real_value(
      matrix, n, lda, workspace->factor, workspace->pivots,
      workspace->factor_work, workspace->factor_lwork,
      scaled_pivot_tolerance, result);
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex_value_with_workspace(
    MVMCAbsolutePfaffianComplexValueWorkspace *workspace,
    const double complex *matrix, int n, int lda,
    double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianValueResult *result) {
  size_t matrix_count;

  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  initialize_value_result(result);
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (workspace == NULL || workspace->n != n || matrix == NULL ||
      n <= 0 || (n % 2) != 0 || lda < n ||
      !isfinite(scaled_pivot_tolerance) ||
      !checked_square_size(n, sizeof(*matrix), &matrix_count) ||
      workspace->matrix_count != matrix_count) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  return factorize_complex_value(
      matrix, n, lda, workspace->factor, workspace->pivots,
      workspace->factor_work, workspace->factor_lwork,
      scaled_pivot_tolerance, result);
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_real_value(
    const double *matrix, int n, int lda, double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianValueResult *result) {
  MVMCAbsolutePfaffianRealValueWorkspace *workspace = NULL;
  MVMCPfaffianStatus status;
  size_t matrix_count;

  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  initialize_value_result(result);
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (matrix == NULL || n <= 0 || (n % 2) != 0 || lda < n ||
      !isfinite(scaled_pivot_tolerance) ||
      !checked_square_size(n, sizeof(*matrix), &matrix_count)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_absolute_pfaffian_real_value_workspace_create(n, &workspace);
  if (status == MVMC_PFAFFIAN_STATUS_OK) {
    status = mvmc_absolute_pfaffian_real_value_with_workspace(
        workspace, matrix, n, lda, scaled_pivot_tolerance, result);
  }
  mvmc_absolute_pfaffian_real_value_workspace_destroy(workspace);
  return status;
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex_value(
    const double complex *matrix, int n, int lda,
    double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianValueResult *result) {
  MVMCAbsolutePfaffianComplexValueWorkspace *workspace = NULL;
  MVMCPfaffianStatus status;
  size_t matrix_count;

  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  initialize_value_result(result);
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (matrix == NULL || n <= 0 || (n % 2) != 0 || lda < n ||
      !isfinite(scaled_pivot_tolerance) ||
      !checked_square_size(n, sizeof(*matrix), &matrix_count)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  status =
      mvmc_absolute_pfaffian_complex_value_workspace_create(n, &workspace);
  if (status == MVMC_PFAFFIAN_STATUS_OK) {
    status = mvmc_absolute_pfaffian_complex_value_with_workspace(
        workspace, matrix, n, lda, scaled_pivot_tolerance, result);
  }
  mvmc_absolute_pfaffian_complex_value_workspace_destroy(workspace);
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
    if (!isfinite(updated)) {
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    abs_compensation = (updated - abs_sum) - corrected;
    abs_sum = updated;

    value = creal(term);
    corrected = value - real_compensation;
    updated = real_sum + corrected;
    if (!isfinite(updated)) {
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    real_compensation = (updated - real_sum) - corrected;
    real_sum = updated;

    value = cimag(term);
    corrected = value - imag_compensation;
    updated = imag_sum + corrected;
    if (!isfinite(updated)) {
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
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

MVMCPfaffianStatus mvmc_projected_amplitude_values(
    const MVMCAbsolutePfaffianValueResult *components,
    const double complex *weights, size_t component_count,
    MVMCProjectedAmplitudeResult *result) {
  size_t index;
  double real_sum = 0.0, real_compensation = 0.0;
  double imag_sum = 0.0, imag_compensation = 0.0;
  double abs_sum = 0.0, abs_compensation = 0.0;
  MVMCProjectedAmplitudeResult accumulated;

  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
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
        components[index].state == MVMC_PFAFFIAN_VALUE_NONFINITE ||
        components[index].state == MVMC_PFAFFIAN_VALUE_INVALID) {
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    switch (components[index].state) {
      case MVMC_PFAFFIAN_VALUE_WELL_PIVOTED:
        ++accumulated.regular_count;
        break;
      case MVMC_PFAFFIAN_VALUE_NEAR_PIVOT:
        ++accumulated.near_singular_count;
        break;
      case MVMC_PFAFFIAN_VALUE_SINGULAR:
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
    if (!isfinite(updated)) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    abs_compensation = (updated - abs_sum) - corrected;
    abs_sum = updated;

    value = creal(term);
    corrected = value - real_compensation;
    updated = real_sum + corrected;
    if (!isfinite(updated)) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    real_compensation = (updated - real_sum) - corrected;
    real_sum = updated;

    value = cimag(term);
    corrected = value - imag_compensation;
    updated = imag_sum + corrected;
    if (!isfinite(updated)) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    imag_compensation = (updated - imag_sum) - corrected;
    imag_sum = updated;
  }
  accumulated.total = real_sum + I * imag_sum;
  accumulated.sum_abs = abs_sum;
  accumulated.cancellation_ratio =
      abs_sum == 0.0 ? 0.0 : fmin(1.0, cabs(accumulated.total) / abs_sum);
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

MVMCPfaffianStatus mvmc_projected_amplitude_value_slice(
    const MVMCAbsolutePfaffianValueResult *local_components,
    size_t local_component_count, const double complex *global_weights,
    int qp_total, int qp_start, int qp_end,
    MVMCProjectedAmplitudeResult *result) {
  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  memset(result, 0, sizeof(*result));
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (qp_total <= 0 || qp_start < 0 || qp_end < qp_start ||
      qp_end > qp_total ||
      (size_t)(qp_end - qp_start) != local_component_count ||
      global_weights == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  if (local_component_count == 0) {
    result->valid = 1;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  return mvmc_projected_amplitude_values(
      local_components, global_weights + qp_start,
      local_component_count, result);
}

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
static int finite_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static int negative_infinity(double value) {
  return isinf(value) && value < 0.0;
}

static double log_add_bounds(double left, double right) {
  double maximum;
  if (negative_infinity(left)) return right;
  if (negative_infinity(right)) return left;
  if (!isfinite(left) || !isfinite(right)) return INFINITY;
  maximum = fmax(left, right);
  return maximum + log1p(exp(fmin(left, right) - maximum));
}

static MVMCScaledComplex scaled_nonfinite(void) {
  MVMCScaledComplex value;
  value.state = MVMC_SCALED_COMPLEX_NONFINITE;
  value.phase = NAN + I * NAN;
  value.log_abs = NAN;
  value.log_abs_error_bound = NAN;
  value.max_input_log_abs = NAN;
  value.cancellation_log_abs = NAN;
  value.cancellation_ratio = NAN;
  return value;
}

static MVMCScaledComplex scaled_exact_zero(void) {
  MVMCScaledComplex value;
  value.state = MVMC_SCALED_COMPLEX_EXACT_ZERO;
  value.phase = 1.0;
  value.log_abs = -INFINITY;
  value.log_abs_error_bound = -INFINITY;
  value.max_input_log_abs = -INFINITY;
  value.cancellation_log_abs = -INFINITY;
  value.cancellation_ratio = 0.0;
  return value;
}

static MVMCScaledComplex scaled_numeric_zero(
    double log_abs_error_bound, double max_input_log_abs,
    double cancellation_log_abs, double cancellation_ratio) {
  MVMCScaledComplex value;
  value.state = MVMC_SCALED_COMPLEX_NUMERIC_ZERO;
  value.phase = 1.0;
  value.log_abs = -INFINITY;
  value.log_abs_error_bound = log_abs_error_bound;
  value.max_input_log_abs = max_input_log_abs;
  value.cancellation_log_abs = cancellation_log_abs;
  value.cancellation_ratio = cancellation_ratio;
  return value;
}

static int valid_error_log(double value) {
  return isfinite(value) || negative_infinity(value);
}

static int valid_unit_phase(double complex phase) {
  double magnitude;
  if (!finite_complex(phase)) return 0;
  magnitude = cabs(phase);
  return isfinite(magnitude) && magnitude > 0.0 &&
         fabs(magnitude - 1.0) <= 32.0 * DBL_EPSILON;
}

static int decompose_finite_complex(
    double complex value, double complex *phase, double *log_abs) {
  const double real_abs = fabs(creal(value));
  const double imag_abs = fabs(cimag(value));
  const double scale = fmax(real_abs, imag_abs);
  double complex scaled;
  double scaled_abs;
  if (phase == NULL || log_abs == NULL || !finite_complex(value) ||
      scale == 0.0 || !isfinite(scale)) {
    return 0;
  }
  scaled = (creal(value) / scale) + I * (cimag(value) / scale);
  scaled_abs = hypot(creal(scaled), cimag(scaled));
  if (!isfinite(scaled_abs) || scaled_abs == 0.0) return 0;
  *phase = scaled / scaled_abs;
  *log_abs = log(scale) + log(scaled_abs);
  return valid_unit_phase(*phase) && isfinite(*log_abs);
}

int mvmc_scaled_complex_is_valid(const MVMCScaledComplex *value) {
  if (value == NULL) return 0;
  switch (value->state) {
    case MVMC_SCALED_COMPLEX_FINITE_NONZERO:
      return valid_unit_phase(value->phase) &&
             isfinite(value->log_abs) &&
             valid_error_log(value->log_abs_error_bound) &&
             (negative_infinity(value->log_abs_error_bound) ||
              value->log_abs_error_bound < value->log_abs) &&
             isfinite(value->max_input_log_abs) &&
             isfinite(value->cancellation_log_abs) &&
             isfinite(value->cancellation_ratio) &&
             value->cancellation_ratio >= 0.0 &&
             value->cancellation_ratio <= 1.0;
    case MVMC_SCALED_COMPLEX_EXACT_ZERO:
      return valid_unit_phase(value->phase) &&
             negative_infinity(value->log_abs) &&
             negative_infinity(value->log_abs_error_bound) &&
             negative_infinity(value->max_input_log_abs) &&
             negative_infinity(value->cancellation_log_abs) &&
             value->cancellation_ratio == 0.0;
    case MVMC_SCALED_COMPLEX_NUMERIC_ZERO:
      return valid_unit_phase(value->phase) &&
             negative_infinity(value->log_abs) &&
             isfinite(value->log_abs_error_bound) &&
             valid_error_log(value->max_input_log_abs) &&
             valid_error_log(value->cancellation_log_abs) &&
             isfinite(value->cancellation_ratio) &&
             value->cancellation_ratio >= 0.0 &&
             value->cancellation_ratio <= 1.0;
    case MVMC_SCALED_COMPLEX_NONFINITE:
      return 1;
  }
  return 0;
}

static double phase_rounding_error_log(double log_abs) {
  const double relative_bound = 32.0 * DBL_EPSILON;
  return log_abs + log(relative_bound);
}

static int normalize_phase(double complex phase, double complex *normalized) {
  double ignored_log_abs;
  return decompose_finite_complex(phase, normalized, &ignored_log_abs);
}

MVMCPfaffianStatus mvmc_scaled_complex_make_finite(
    double complex phase, double log_abs, double log_abs_error_bound,
    MVMCScaledComplex *result) {
  MVMCScaledComplex candidate;
  double complex normalized;
  if (result == NULL || !isfinite(log_abs) ||
      !valid_error_log(log_abs_error_bound)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  if (!normalize_phase(phase, &normalized)) {
    candidate = scaled_nonfinite();
    *result = candidate;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (isfinite(log_abs_error_bound) &&
      log_abs <= log_abs_error_bound) {
    candidate = scaled_numeric_zero(log_abs_error_bound, log_abs, log_abs,
                                    1.0);
  } else {
    candidate.state = MVMC_SCALED_COMPLEX_FINITE_NONZERO;
    candidate.phase = normalized;
    candidate.log_abs = log_abs;
    candidate.log_abs_error_bound = log_abs_error_bound;
    candidate.max_input_log_abs = log_abs;
    candidate.cancellation_log_abs = log_abs;
    candidate.cancellation_ratio = 1.0;
  }
  *result = candidate;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_scaled_complex_make_exact_zero(
    MVMCScaledComplex *result) {
  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  *result = scaled_exact_zero();
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_scaled_complex_make_numeric_zero(
    double log_abs_error_bound, double max_input_log_abs,
    double cancellation_log_abs, double cancellation_ratio,
    MVMCScaledComplex *result) {
  if (result == NULL || !isfinite(log_abs_error_bound) ||
      !valid_error_log(max_input_log_abs) ||
      !valid_error_log(cancellation_log_abs) ||
      !isfinite(cancellation_ratio) || cancellation_ratio < 0.0 ||
      cancellation_ratio > 1.0) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  *result = scaled_numeric_zero(log_abs_error_bound, max_input_log_abs,
                                cancellation_log_abs, cancellation_ratio);
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_scaled_complex_from_raw_testing(
    double complex value, MVMCScaledComplex *result) {
  MVMCScaledComplex candidate;
  double complex phase;
  double log_abs;
  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  if (!finite_complex(value)) {
    *result = scaled_nonfinite();
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (creal(value) == 0.0 && cimag(value) == 0.0) {
    /* A raw zero has lost its producer scale, so its error is unknown. */
    *result = scaled_numeric_zero(log(DBL_MAX), -INFINITY, -INFINITY,
                                  0.0);
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (!decompose_finite_complex(value, &phase, &log_abs)) {
    *result = scaled_nonfinite();
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (mvmc_scaled_complex_make_finite(phase, log_abs, -INFINITY,
                                      &candidate) !=
      MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  *result = candidate;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_scaled_complex_multiply(
    const MVMCScaledComplex *left, const MVMCScaledComplex *right,
    MVMCScaledComplex *result) {
  MVMCScaledComplex candidate;
  double error_bound = -INFINITY;
  double output_log_abs;
  double max_input;
  if (result == NULL || !mvmc_scaled_complex_is_valid(left) ||
      !mvmc_scaled_complex_is_valid(right)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  if (left->state == MVMC_SCALED_COMPLEX_NONFINITE ||
      right->state == MVMC_SCALED_COMPLEX_NONFINITE) {
    *result = scaled_nonfinite();
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (left->state == MVMC_SCALED_COMPLEX_EXACT_ZERO ||
      right->state == MVMC_SCALED_COMPLEX_EXACT_ZERO) {
    *result = scaled_exact_zero();
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (left->state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO ||
      right->state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO) {
    if (left->state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO &&
        right->state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO) {
      error_bound = left->log_abs_error_bound +
                    right->log_abs_error_bound;
    } else if (left->state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO) {
      error_bound = left->log_abs_error_bound + right->log_abs;
      if (isfinite(right->log_abs_error_bound)) {
        error_bound = log_add_bounds(
            error_bound, left->log_abs_error_bound +
                             right->log_abs_error_bound);
      }
    } else {
      error_bound = right->log_abs_error_bound + left->log_abs;
      if (isfinite(left->log_abs_error_bound)) {
        error_bound = log_add_bounds(
            error_bound, right->log_abs_error_bound +
                             left->log_abs_error_bound);
      }
    }
    if (!isfinite(error_bound)) {
      *result = scaled_nonfinite();
      return MVMC_PFAFFIAN_STATUS_OK;
    }
    max_input = fmax(left->max_input_log_abs, right->max_input_log_abs);
    *result = scaled_numeric_zero(error_bound, max_input, -INFINITY, 0.0);
    return MVMC_PFAFFIAN_STATUS_OK;
  }

  output_log_abs = left->log_abs + right->log_abs;
  if (!isfinite(output_log_abs)) {
    *result = scaled_nonfinite();
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (isfinite(left->log_abs_error_bound)) {
    error_bound = log_add_bounds(
        error_bound, left->log_abs_error_bound + right->log_abs);
  }
  if (isfinite(right->log_abs_error_bound)) {
    error_bound = log_add_bounds(
        error_bound, right->log_abs_error_bound + left->log_abs);
  }
  if (isfinite(left->log_abs_error_bound) &&
      isfinite(right->log_abs_error_bound)) {
    error_bound = log_add_bounds(
        error_bound, left->log_abs_error_bound +
                         right->log_abs_error_bound);
  }
  error_bound = log_add_bounds(error_bound,
                               phase_rounding_error_log(output_log_abs));
  if (!isfinite(error_bound)) {
    *result = scaled_nonfinite();
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (mvmc_scaled_complex_make_finite(
          left->phase * right->phase, output_log_abs, error_bound,
          &candidate) != MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  candidate.max_input_log_abs =
      fmax(left->max_input_log_abs, right->max_input_log_abs);
  candidate.cancellation_log_abs = output_log_abs;
  candidate.cancellation_ratio = 1.0;
  *result = candidate;
  return MVMC_PFAFFIAN_STATUS_OK;
}

static int neumaier_add(double value, double *sum, double *compensation) {
  double updated;
  if (sum == NULL || compensation == NULL || !isfinite(value)) return 0;
  updated = *sum + value;
  if (!isfinite(updated)) return 0;
  if (fabs(*sum) >= fabs(value)) {
    *compensation += (*sum - updated) + value;
  } else {
    *compensation += (value - updated) + *sum;
  }
  if (!isfinite(*compensation)) return 0;
  *sum = updated;
  return 1;
}

MVMCPfaffianStatus mvmc_scaled_complex_sum_ordered(
    const MVMCScaledComplex *values, size_t value_count,
    MVMCScaledComplex *result) {
  MVMCScaledComplex candidate;
  double scale = -INFINITY;
  double error_bound = -INFINITY;
  double max_input = -INFINITY;
  double real_sum = 0.0, real_compensation = 0.0;
  double imag_sum = 0.0, imag_compensation = 0.0;
  double abs_sum = 0.0, abs_compensation = 0.0;
  double complex central;
  double central_abs;
  double central_log;
  double cancellation_ratio;
  size_t index;
  int finite_count = 0;
  int numeric_count = 0;

  if (result == NULL || values == NULL || value_count == 0) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  for (index = 0; index < value_count; ++index) {
    if (!mvmc_scaled_complex_is_valid(values + index)) {
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    if (values[index].state == MVMC_SCALED_COMPLEX_NONFINITE) {
      *result = scaled_nonfinite();
      return MVMC_PFAFFIAN_STATUS_OK;
    }
    max_input = fmax(max_input, values[index].max_input_log_abs);
    if (values[index].state == MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
      scale = fmax(scale, values[index].log_abs);
      ++finite_count;
    } else if (values[index].state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO) {
      error_bound = log_add_bounds(error_bound,
                                   values[index].log_abs_error_bound);
      ++numeric_count;
    }
  }
  if (finite_count == 0) {
    if (numeric_count == 0) {
      *result = scaled_exact_zero();
    } else if (isfinite(error_bound)) {
      *result = scaled_numeric_zero(error_bound, max_input, -INFINITY, 0.0);
    } else {
      *result = scaled_nonfinite();
    }
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  for (index = 0; index < value_count; ++index) {
    double magnitude;
    double complex term;
    if (values[index].state != MVMC_SCALED_COMPLEX_FINITE_NONZERO) continue;
    magnitude = exp(values[index].log_abs - scale);
    term = magnitude * values[index].phase;
    if (!finite_complex(term) ||
        !neumaier_add(creal(term), &real_sum, &real_compensation) ||
        !neumaier_add(cimag(term), &imag_sum, &imag_compensation) ||
        !neumaier_add(magnitude, &abs_sum, &abs_compensation)) {
      *result = scaled_nonfinite();
      return MVMC_PFAFFIAN_STATUS_OK;
    }
    error_bound = log_add_bounds(error_bound,
                                 values[index].log_abs_error_bound);
  }
  central = (real_sum + real_compensation) +
            I * (imag_sum + imag_compensation);
  central_abs = cabs(central);
  abs_sum += abs_compensation;
  if (!finite_complex(central) || !isfinite(abs_sum) || abs_sum <= 0.0) {
    *result = scaled_nonfinite();
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  cancellation_ratio = fmin(1.0, central_abs / abs_sum);
  if (central_abs == 0.0) {
    error_bound = log_add_bounds(
        error_bound, scale + log(32.0 * DBL_EPSILON * abs_sum));
    *result = isfinite(error_bound)
                  ? scaled_numeric_zero(error_bound, max_input, -INFINITY,
                                        cancellation_ratio)
                  : scaled_nonfinite();
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  central_log = scale + log(central_abs);
  error_bound = log_add_bounds(
      error_bound, scale + log(32.0 * DBL_EPSILON * abs_sum));
  if (!isfinite(central_log) || !isfinite(error_bound)) {
    *result = scaled_nonfinite();
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (central_log <= error_bound) {
    *result = scaled_numeric_zero(error_bound, max_input, central_log,
                                  cancellation_ratio);
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (mvmc_scaled_complex_make_finite(central / central_abs, central_log,
                                      error_bound, &candidate) !=
      MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  candidate.max_input_log_abs = max_input;
  candidate.cancellation_log_abs = central_log;
  candidate.cancellation_ratio = cancellation_ratio;
  *result = candidate;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCScaledComplexExportStatus mvmc_scaled_complex_export_common_scale(
    const MVMCScaledComplex *value, double common_log_scale,
    double complex *result) {
  double delta;
  double magnitude;
  if (result == NULL || !mvmc_scaled_complex_is_valid(value) ||
      !isfinite(common_log_scale)) {
    return MVMC_SCALED_EXPORT_INVALID;
  }
  switch (value->state) {
    case MVMC_SCALED_COMPLEX_EXACT_ZERO:
      *result = 0.0;
      return MVMC_SCALED_EXPORT_EXACT_ZERO;
    case MVMC_SCALED_COMPLEX_NUMERIC_ZERO:
      *result = 0.0;
      return MVMC_SCALED_EXPORT_NUMERIC_ZERO;
    case MVMC_SCALED_COMPLEX_NONFINITE:
      *result = NAN + I * NAN;
      return MVMC_SCALED_EXPORT_NONFINITE;
    case MVMC_SCALED_COMPLEX_FINITE_NONZERO:
      break;
  }
  delta = value->log_abs - common_log_scale;
  if (!isfinite(delta)) return MVMC_SCALED_EXPORT_INVALID;
  if (delta > log(DBL_MAX)) {
    *result = NAN + I * NAN;
    return MVMC_SCALED_EXPORT_OVERFLOW;
  }
  if (delta < log(nextafter(0.0, 1.0))) {
    *result = 0.0;
    return MVMC_SCALED_EXPORT_UNDERFLOW;
  }
  magnitude = exp(delta);
  if (magnitude == 0.0) {
    *result = 0.0;
    return MVMC_SCALED_EXPORT_UNDERFLOW;
  }
  *result = magnitude * value->phase;
  if (!finite_complex(*result)) {
    *result = NAN + I * NAN;
    return MVMC_SCALED_EXPORT_OVERFLOW;
  }
  return MVMC_SCALED_EXPORT_OK;
}

static MVMCScaledComplex scaled_exact_weight(double complex weight) {
  MVMCScaledComplex value;
  double complex phase;
  double log_abs;
  if (!finite_complex(weight)) return scaled_nonfinite();
  if (creal(weight) == 0.0 && cimag(weight) == 0.0) {
    return scaled_exact_zero();
  }
  if (!decompose_finite_complex(weight, &phase, &log_abs)) {
    return scaled_nonfinite();
  }
  if (mvmc_scaled_complex_make_finite(phase, log_abs, -INFINITY, &value) !=
      MVMC_PFAFFIAN_STATUS_OK) {
    return scaled_nonfinite();
  }
  return value;
}

static double zero_pivot_error_log(double matrix_scale, int n) {
  double bound = 64.0 * DBL_EPSILON * (double)n;
  if (!isfinite(matrix_scale) || matrix_scale < 0.0) return NAN;
  if (matrix_scale > DBL_MAX / bound) return log(DBL_MAX);
  bound *= fmax(matrix_scale, DBL_MIN);
  if (!isfinite(bound) || bound <= 0.0) return log(DBL_MAX);
  return log(bound);
}

static void initialize_scaled_pfaffian_result(
    MVMCAbsolutePfaffianScaledValueResult *result) {
  result->factor_state = MVMC_PFAFFIAN_VALUE_INVALID;
  result->factor_info = 0;
  result->matrix_scale = NAN;
  result->scaled_min_pivot = NAN;
  result->value = scaled_nonfinite();
}

static MVMCPfaffianStatus factorize_real_scaled_value(
    const double *matrix, int n, int lda, double *factor, int *pivots,
    double *factor_work, int factor_lwork, double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianScaledValueResult *result) {
  MVMCAbsolutePfaffianScaledValueResult candidate;
  MVMCScaledComplex product;
  double min_pivot = INFINITY;
  int column, row, info = 0;
  initialize_scaled_pfaffian_result(&candidate);
  info = real_matrix_properties(matrix, n, lda, &candidate.matrix_scale);
  if (info < 0) {
    candidate.factor_state = MVMC_PFAFFIAN_VALUE_NONFINITE;
    candidate.value = scaled_nonfinite();
    *result = candidate;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (info == 0) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  for (column = 0; column < n; ++column) {
    pivots[column] = column + 1;
    for (row = 0; row < n; ++row) {
      factor[row + (size_t)column * (size_t)n] =
          matrix[row + (size_t)column * (size_t)lda];
    }
  }
  M_DSKTRF("U", "P", &n, factor, &n, pivots, factor_work,
           &factor_lwork, &info);
  candidate.factor_info = info;
  if (info < 0) return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
  if (info > 0) {
    candidate.factor_state = MVMC_PFAFFIAN_VALUE_SINGULAR;
    candidate.scaled_min_pivot = 0.0;
    candidate.value = scaled_exact_zero();
    *result = candidate;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (mvmc_scaled_complex_make_finite(1.0, 0.0, -INFINITY, &product) !=
      MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
  }
  for (row = 0; row < n; row += 2) {
    MVMCScaledComplex pivot_value;
    MVMCScaledComplex updated;
    double pivot = factor[row + (size_t)(row + 1) * (size_t)n];
    const double pivot_abs = fabs(pivot);
    if (!isfinite(pivot)) {
      candidate.factor_state = MVMC_PFAFFIAN_VALUE_NONFINITE;
      candidate.value = scaled_nonfinite();
      *result = candidate;
      return MVMC_PFAFFIAN_STATUS_OK;
    }
    min_pivot = fmin(min_pivot, pivot_abs);
    if (pivot == 0.0) {
      pivot_value = scaled_numeric_zero(
          zero_pivot_error_log(candidate.matrix_scale, n), -INFINITY,
          -INFINITY, 0.0);
    } else if (mvmc_scaled_complex_from_raw_testing(pivot, &pivot_value) !=
               MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
    }
    if (pivots[row] - 1 != row) pivot_value.phase = -pivot_value.phase;
    if (pivots[row + 1] - 1 != row + 1) {
      pivot_value.phase = -pivot_value.phase;
    }
    if (mvmc_scaled_complex_multiply(&product, &pivot_value, &updated) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
    }
    product = updated;
  }
  candidate.scaled_min_pivot =
      min_pivot / fmax(candidate.matrix_scale, DBL_MIN);
  if (!isfinite(candidate.scaled_min_pivot) ||
      product.state == MVMC_SCALED_COMPLEX_NONFINITE) {
    candidate.factor_state = MVMC_PFAFFIAN_VALUE_NONFINITE;
    candidate.value = scaled_nonfinite();
  } else if (candidate.scaled_min_pivot <
                 effective_pivot_tolerance(scaled_pivot_tolerance, n) ||
             product.state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO) {
    candidate.factor_state = MVMC_PFAFFIAN_VALUE_NEAR_PIVOT;
    candidate.value = product;
  } else {
    candidate.factor_state = MVMC_PFAFFIAN_VALUE_WELL_PIVOTED;
    candidate.value = product;
  }
  *result = candidate;
  return MVMC_PFAFFIAN_STATUS_OK;
}

static MVMCPfaffianStatus factorize_complex_scaled_value(
    const double complex *matrix, int n, int lda, double complex *factor,
    int *pivots, double complex *factor_work, int factor_lwork,
    double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianScaledValueResult *result) {
  MVMCAbsolutePfaffianScaledValueResult candidate;
  MVMCScaledComplex product;
  double min_pivot = INFINITY;
  int column, row, info = 0;
  initialize_scaled_pfaffian_result(&candidate);
  info = complex_matrix_properties(matrix, n, lda, &candidate.matrix_scale);
  if (info < 0) {
    candidate.factor_state = MVMC_PFAFFIAN_VALUE_NONFINITE;
    candidate.value = scaled_nonfinite();
    *result = candidate;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (info == 0) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  for (column = 0; column < n; ++column) {
    pivots[column] = column + 1;
    for (row = 0; row < n; ++row) {
      factor[row + (size_t)column * (size_t)n] =
          matrix[row + (size_t)column * (size_t)lda];
    }
  }
  M_ZSKTRF("U", "P", &n, factor, &n, pivots, factor_work,
           &factor_lwork, &info);
  candidate.factor_info = info;
  if (info < 0) return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
  if (info > 0) {
    candidate.factor_state = MVMC_PFAFFIAN_VALUE_SINGULAR;
    candidate.scaled_min_pivot = 0.0;
    candidate.value = scaled_exact_zero();
    *result = candidate;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (mvmc_scaled_complex_make_finite(1.0, 0.0, -INFINITY, &product) !=
      MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
  }
  for (row = 0; row < n; row += 2) {
    MVMCScaledComplex pivot_value;
    MVMCScaledComplex updated;
    double complex pivot = factor[row + (size_t)(row + 1) * (size_t)n];
    const double pivot_abs = cabs(pivot);
    if (!finite_complex(pivot) || !isfinite(pivot_abs)) {
      candidate.factor_state = MVMC_PFAFFIAN_VALUE_NONFINITE;
      candidate.value = scaled_nonfinite();
      *result = candidate;
      return MVMC_PFAFFIAN_STATUS_OK;
    }
    min_pivot = fmin(min_pivot, pivot_abs);
    if (pivot_abs == 0.0) {
      pivot_value = scaled_numeric_zero(
          zero_pivot_error_log(candidate.matrix_scale, n), -INFINITY,
          -INFINITY, 0.0);
    } else if (mvmc_scaled_complex_from_raw_testing(pivot, &pivot_value) !=
               MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
    }
    if (pivots[row] - 1 != row) pivot_value.phase = -pivot_value.phase;
    if (pivots[row + 1] - 1 != row + 1) {
      pivot_value.phase = -pivot_value.phase;
    }
    if (mvmc_scaled_complex_multiply(&product, &pivot_value, &updated) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
    }
    product = updated;
  }
  candidate.scaled_min_pivot =
      min_pivot / fmax(candidate.matrix_scale, DBL_MIN);
  if (!isfinite(candidate.scaled_min_pivot) ||
      product.state == MVMC_SCALED_COMPLEX_NONFINITE) {
    candidate.factor_state = MVMC_PFAFFIAN_VALUE_NONFINITE;
    candidate.value = scaled_nonfinite();
  } else if (candidate.scaled_min_pivot <
                 effective_pivot_tolerance(scaled_pivot_tolerance, n) ||
             product.state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO) {
    candidate.factor_state = MVMC_PFAFFIAN_VALUE_NEAR_PIVOT;
    candidate.value = product;
  } else {
    candidate.factor_state = MVMC_PFAFFIAN_VALUE_WELL_PIVOTED;
    candidate.value = product;
  }
  *result = candidate;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_real_scaled_value_with_workspace(
    MVMCAbsolutePfaffianRealValueWorkspace *workspace,
    const double *matrix, int n, int lda, double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianScaledValueResult *result) {
  size_t matrix_count;
  if (result == NULL || !mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return result == NULL ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT
                          : MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (workspace == NULL || workspace->n != n || matrix == NULL || n <= 0 ||
      (n % 2) != 0 || lda < n || !isfinite(scaled_pivot_tolerance) ||
      !checked_square_size(n, sizeof(*matrix), &matrix_count) ||
      workspace->matrix_count != matrix_count) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  return factorize_real_scaled_value(
      matrix, n, lda, workspace->factor, workspace->pivots,
      workspace->factor_work, workspace->factor_lwork,
      scaled_pivot_tolerance, result);
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex_scaled_value_with_workspace(
    MVMCAbsolutePfaffianComplexValueWorkspace *workspace,
    const double complex *matrix, int n, int lda,
    double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianScaledValueResult *result) {
  size_t matrix_count;
  if (result == NULL || !mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return result == NULL ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT
                          : MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (workspace == NULL || workspace->n != n || matrix == NULL || n <= 0 ||
      (n % 2) != 0 || lda < n || !isfinite(scaled_pivot_tolerance) ||
      !checked_square_size(n, sizeof(*matrix), &matrix_count) ||
      workspace->matrix_count != matrix_count) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  return factorize_complex_scaled_value(
      matrix, n, lda, workspace->factor, workspace->pivots,
      workspace->factor_work, workspace->factor_lwork,
      scaled_pivot_tolerance, result);
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_real_scaled_value(
    const double *matrix, int n, int lda, double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianScaledValueResult *result) {
  MVMCAbsolutePfaffianRealValueWorkspace *workspace = NULL;
  MVMCPfaffianStatus status;
  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  status = mvmc_absolute_pfaffian_real_value_workspace_create(n, &workspace);
  if (status == MVMC_PFAFFIAN_STATUS_OK) {
    status = mvmc_absolute_pfaffian_real_scaled_value_with_workspace(
        workspace, matrix, n, lda, scaled_pivot_tolerance, result);
  }
  mvmc_absolute_pfaffian_real_value_workspace_destroy(workspace);
  return status;
}

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex_scaled_value(
    const double complex *matrix, int n, int lda,
    double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianScaledValueResult *result) {
  MVMCAbsolutePfaffianComplexValueWorkspace *workspace = NULL;
  MVMCPfaffianStatus status;
  if (result == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  status =
      mvmc_absolute_pfaffian_complex_value_workspace_create(n, &workspace);
  if (status == MVMC_PFAFFIAN_STATUS_OK) {
    status = mvmc_absolute_pfaffian_complex_scaled_value_with_workspace(
        workspace, matrix, n, lda, scaled_pivot_tolerance, result);
  }
  mvmc_absolute_pfaffian_complex_value_workspace_destroy(workspace);
  return status;
}

MVMCPfaffianStatus mvmc_projected_scaled_amplitude_values(
    const MVMCAbsolutePfaffianScaledValueResult *components,
    const double complex *weights, size_t component_count,
    MVMCProjectedScaledAmplitudeResult *result) {
  MVMCProjectedScaledAmplitudeResult candidate;
  MVMCScaledComplex *terms;
  double log_sum_abs = -INFINITY;
  size_t index;
  MVMCPfaffianStatus status;
  if (result == NULL || components == NULL || weights == NULL ||
      component_count == 0 ||
      component_count > SIZE_MAX / sizeof(*terms)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  terms = (MVMCScaledComplex *)malloc(component_count * sizeof(*terms));
  if (terms == NULL) return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  memset(&candidate, 0, sizeof(candidate));
  candidate.total = scaled_nonfinite();
  candidate.log_sum_abs = -INFINITY;
  for (index = 0; index < component_count; ++index) {
    MVMCScaledComplex scaled_weight;
    if (!finite_complex(weights[index]) ||
        !mvmc_scaled_complex_is_valid(&components[index].value)) {
      free(terms);
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    switch (components[index].factor_state) {
      case MVMC_PFAFFIAN_VALUE_WELL_PIVOTED:
        if (components[index].value.state !=
            MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
          free(terms);
          return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
        }
        ++candidate.well_pivoted_count;
        break;
      case MVMC_PFAFFIAN_VALUE_NEAR_PIVOT:
        if (components[index].value.state !=
                MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
            components[index].value.state !=
                MVMC_SCALED_COMPLEX_NUMERIC_ZERO) {
          free(terms);
          return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
        }
        ++candidate.near_pivot_count;
        break;
      case MVMC_PFAFFIAN_VALUE_SINGULAR:
        if (components[index].value.state !=
            MVMC_SCALED_COMPLEX_EXACT_ZERO) {
          free(terms);
          return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
        }
        break;
      case MVMC_PFAFFIAN_VALUE_NONFINITE:
      case MVMC_PFAFFIAN_VALUE_INVALID:
        free(terms);
        return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    if (components[index].value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO) {
      ++candidate.exact_zero_count;
    } else if (components[index].value.state ==
               MVMC_SCALED_COMPLEX_NUMERIC_ZERO) {
      ++candidate.numeric_zero_count;
    } else if (components[index].value.state ==
               MVMC_SCALED_COMPLEX_NONFINITE) {
      free(terms);
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    scaled_weight = scaled_exact_weight(weights[index]);
    if (scaled_weight.state == MVMC_SCALED_COMPLEX_NONFINITE) {
      free(terms);
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    status = mvmc_scaled_complex_multiply(
        &scaled_weight, &components[index].value, terms + index);
    if (status != MVMC_PFAFFIAN_STATUS_OK ||
        terms[index].state == MVMC_SCALED_COMPLEX_NONFINITE) {
      free(terms);
      return status == MVMC_PFAFFIAN_STATUS_OK
                 ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT
                 : status;
    }
    if (terms[index].state == MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
      log_sum_abs = log_add_bounds(log_sum_abs, terms[index].log_abs);
    } else if (terms[index].state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO) {
      log_sum_abs = log_add_bounds(log_sum_abs,
                                   terms[index].log_abs_error_bound);
    }
  }
  status = mvmc_scaled_complex_sum_ordered(terms, component_count,
                                           &candidate.total);
  free(terms);
  if (status != MVMC_PFAFFIAN_STATUS_OK ||
      candidate.total.state == MVMC_SCALED_COMPLEX_NONFINITE ||
      !valid_error_log(log_sum_abs)) {
    return status == MVMC_PFAFFIAN_STATUS_OK
               ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT
               : status;
  }
  candidate.log_sum_abs = log_sum_abs;
  candidate.valid = 1;
  *result = candidate;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_projected_scaled_amplitude_value_slice(
    const MVMCAbsolutePfaffianScaledValueResult *local_components,
    size_t local_component_count, const double complex *global_weights,
    int qp_total, int qp_start, int qp_end,
    MVMCProjectedScaledAmplitudeResult *result) {
  MVMCProjectedScaledAmplitudeResult candidate;
  if (result == NULL || qp_total <= 0 || qp_start < 0 || qp_end < qp_start ||
      qp_end > qp_total ||
      (size_t)(qp_end - qp_start) != local_component_count ||
      global_weights == NULL ||
      (local_component_count != 0 && local_components == NULL)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  if (local_component_count == 0) {
    memset(&candidate, 0, sizeof(candidate));
    candidate.valid = 1;
    candidate.total = scaled_exact_zero();
    candidate.log_sum_abs = -INFINITY;
    *result = candidate;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  return mvmc_projected_scaled_amplitude_values(
      local_components, global_weights + qp_start, local_component_count,
      result);
}
#endif /* MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE */

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

const char *mvmc_pfaffian_value_state_name(MVMCPfaffianValueState state) {
  switch (state) {
    case MVMC_PFAFFIAN_VALUE_WELL_PIVOTED:
      return "WELL_PIVOTED";
    case MVMC_PFAFFIAN_VALUE_NEAR_PIVOT:
      return "NEAR_PIVOT";
    case MVMC_PFAFFIAN_VALUE_SINGULAR:
      return "SINGULAR";
    case MVMC_PFAFFIAN_VALUE_NONFINITE:
      return "NONFINITE";
    case MVMC_PFAFFIAN_VALUE_INVALID:
      return "INVALID";
  }
  return "INVALID";
}
