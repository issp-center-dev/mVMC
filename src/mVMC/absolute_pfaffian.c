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
