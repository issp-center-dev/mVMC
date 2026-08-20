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
  MVMC_PFAFFIAN_VALUE_WELL_PIVOTED = 0,
  MVMC_PFAFFIAN_VALUE_NEAR_PIVOT,
  MVMC_PFAFFIAN_VALUE_SINGULAR,
  MVMC_PFAFFIAN_VALUE_NONFINITE,
  MVMC_PFAFFIAN_VALUE_INVALID
} MVMCPfaffianValueState;

typedef enum {
  /* Keep this list contiguous and synchronize collective status severity. */
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
  MVMCPfaffianValueState state;
  int factor_info;
  double matrix_scale;
  double scaled_min_pivot;
  double complex pfaffian;
} MVMCAbsolutePfaffianValueResult;

typedef struct {
  int valid;
  double complex total;
  double sum_abs;
  double cancellation_ratio;
  size_t regular_count;
  size_t near_singular_count;
  size_t singular_count;
} MVMCProjectedAmplitudeResult;

typedef struct MVMCAbsolutePfaffianRealWorkspace
    MVMCAbsolutePfaffianRealWorkspace;
typedef struct MVMCAbsolutePfaffianComplexWorkspace
    MVMCAbsolutePfaffianComplexWorkspace;
typedef struct MVMCAbsolutePfaffianRealValueWorkspace
    MVMCAbsolutePfaffianRealValueWorkspace;
typedef struct MVMCAbsolutePfaffianComplexValueWorkspace
    MVMCAbsolutePfaffianComplexValueWorkspace;

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
/*
 * Testing-gated P4 numeric core.  Production producers must retain the
 * phase/log representation until the P4-C promotion gate has passed.
 */
typedef enum {
  MVMC_SCALED_COMPLEX_FINITE_NONZERO = 0,
  MVMC_SCALED_COMPLEX_EXACT_ZERO,
  MVMC_SCALED_COMPLEX_NUMERIC_ZERO,
  MVMC_SCALED_COMPLEX_NONFINITE
} MVMCScaledComplexState;

typedef enum {
  MVMC_SCALED_EXPORT_OK = 0,
  MVMC_SCALED_EXPORT_EXACT_ZERO,
  MVMC_SCALED_EXPORT_NUMERIC_ZERO,
  MVMC_SCALED_EXPORT_UNDERFLOW,
  MVMC_SCALED_EXPORT_OVERFLOW,
  MVMC_SCALED_EXPORT_NONFINITE,
  MVMC_SCALED_EXPORT_INVALID
} MVMCScaledComplexExportStatus;

typedef struct {
  MVMCScaledComplexState state;
  double complex phase;
  double log_abs;
  double log_abs_error_bound;
  double max_input_log_abs;
  double cancellation_log_abs;
  double cancellation_ratio;
} MVMCScaledComplex;

typedef struct {
  MVMCPfaffianValueState factor_state;
  int factor_info;
  double matrix_scale;
  double scaled_min_pivot;
  MVMCScaledComplex value;
} MVMCAbsolutePfaffianScaledValueResult;

typedef struct {
  int valid;
  MVMCScaledComplex total;
  double log_sum_abs;
  size_t well_pivoted_count;
  size_t near_pivot_count;
  size_t exact_zero_count;
  size_t numeric_zero_count;
} MVMCProjectedScaledAmplitudeResult;
#endif

/*
 * Reusable factorization workspaces for sampler hot paths.  Creation may
 * allocate and query LAPACK/PFAPACK workspace sizes; evaluation with an
 * existing workspace performs no allocation.  A workspace is bound to n and
 * requires exclusive use while an evaluation is in progress.
 */
MVMCPfaffianStatus mvmc_absolute_pfaffian_real_workspace_create(
    int n, MVMCAbsolutePfaffianRealWorkspace **workspace);
void mvmc_absolute_pfaffian_real_workspace_destroy(
    MVMCAbsolutePfaffianRealWorkspace *workspace);

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex_workspace_create(
    int n, MVMCAbsolutePfaffianComplexWorkspace **workspace);
void mvmc_absolute_pfaffian_complex_workspace_destroy(
    MVMCAbsolutePfaffianComplexWorkspace *workspace);

/*
 * Value-only workspaces contain the PFAPACK factorization storage but no
 * inverse factor, inverse workspace, or inverse-residual machinery.  They
 * are independent from the full P1 workspaces and are intended for absolute
 * terminal-amplitude evaluation.
 */
MVMCPfaffianStatus mvmc_absolute_pfaffian_real_value_workspace_create(
    int n, MVMCAbsolutePfaffianRealValueWorkspace **workspace);
void mvmc_absolute_pfaffian_real_value_workspace_destroy(
    MVMCAbsolutePfaffianRealValueWorkspace *workspace);
size_t mvmc_absolute_pfaffian_real_value_workspace_bytes(
    const MVMCAbsolutePfaffianRealValueWorkspace *workspace);

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex_value_workspace_create(
    int n, MVMCAbsolutePfaffianComplexValueWorkspace **workspace);
void mvmc_absolute_pfaffian_complex_value_workspace_destroy(
    MVMCAbsolutePfaffianComplexValueWorkspace *workspace);
size_t mvmc_absolute_pfaffian_complex_value_workspace_bytes(
    const MVMCAbsolutePfaffianComplexValueWorkspace *workspace);

MVMCPfaffianStatus mvmc_absolute_pfaffian_real_with_workspace(
    MVMCAbsolutePfaffianRealWorkspace *workspace,
    const double *matrix, int n, int lda,
    double *inverse_out, int inverse_lda,
    uint64_t rebuild_generation,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCAbsolutePfaffianResult *result);

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex_with_workspace(
    MVMCAbsolutePfaffianComplexWorkspace *workspace,
    const double complex *matrix, int n, int lda,
    double complex *inverse_out, int inverse_lda,
    uint64_t rebuild_generation,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCAbsolutePfaffianResult *result);

MVMCPfaffianStatus mvmc_absolute_pfaffian_real_value_with_workspace(
    MVMCAbsolutePfaffianRealValueWorkspace *workspace,
    const double *matrix, int n, int lda,
    double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianValueResult *result);

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex_value_with_workspace(
    MVMCAbsolutePfaffianComplexValueWorkspace *workspace,
    const double complex *matrix, int n, int lda,
    double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianValueResult *result);

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

MVMCPfaffianStatus mvmc_absolute_pfaffian_real_value(
    const double *matrix, int n, int lda,
    double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianValueResult *result);

MVMCPfaffianStatus mvmc_absolute_pfaffian_complex_value(
    const double complex *matrix, int n, int lda,
    double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianValueResult *result);

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
/* Pure scaled-complex operations.  Error returns preserve *result. */
MVMCPfaffianStatus mvmc_scaled_complex_make_finite(
    double complex phase, double log_abs, double log_abs_error_bound,
    MVMCScaledComplex *result);
MVMCPfaffianStatus mvmc_scaled_complex_make_exact_zero(
    MVMCScaledComplex *result);
MVMCPfaffianStatus mvmc_scaled_complex_make_numeric_zero(
    double log_abs_error_bound, double max_input_log_abs,
    double cancellation_log_abs, double cancellation_ratio,
    MVMCScaledComplex *result);
int mvmc_scaled_complex_is_valid(const MVMCScaledComplex *value);
MVMCPfaffianStatus mvmc_scaled_complex_from_raw(
    double complex value, MVMCScaledComplex *result);
MVMCPfaffianStatus mvmc_scaled_complex_multiply(
    const MVMCScaledComplex *left, const MVMCScaledComplex *right,
    MVMCScaledComplex *result);
MVMCPfaffianStatus mvmc_scaled_complex_sum_ordered(
    const MVMCScaledComplex *values, size_t value_count,
    MVMCScaledComplex *result);
MVMCScaledComplexExportStatus mvmc_scaled_complex_export_common_scale(
    const MVMCScaledComplex *value, double common_log_scale,
    double complex *result);

/* Value-only Pfaffian path that never forms the raw pivot product. */
MVMCPfaffianStatus mvmc_absolute_pfaffian_real_scaled_value_with_workspace(
    MVMCAbsolutePfaffianRealValueWorkspace *workspace,
    const double *matrix, int n, int lda, double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianScaledValueResult *result);
MVMCPfaffianStatus mvmc_absolute_pfaffian_complex_scaled_value_with_workspace(
    MVMCAbsolutePfaffianComplexValueWorkspace *workspace,
    const double complex *matrix, int n, int lda,
    double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianScaledValueResult *result);
MVMCPfaffianStatus mvmc_absolute_pfaffian_real_scaled_value(
    const double *matrix, int n, int lda, double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianScaledValueResult *result);
MVMCPfaffianStatus mvmc_absolute_pfaffian_complex_scaled_value(
    const double complex *matrix, int n, int lda,
    double scaled_pivot_tolerance,
    MVMCAbsolutePfaffianScaledValueResult *result);

/* Canonical global-QP ordered aggregation and its rank-local slice form. */
MVMCPfaffianStatus mvmc_projected_scaled_amplitude_values(
    const MVMCAbsolutePfaffianScaledValueResult *components,
    const double complex *weights, size_t component_count,
    MVMCProjectedScaledAmplitudeResult *result);
MVMCPfaffianStatus mvmc_projected_scaled_amplitude_value_slice(
    const MVMCAbsolutePfaffianScaledValueResult *local_components,
    size_t local_component_count, const double complex *global_weights,
    int qp_total, int qp_start, int qp_end,
    MVMCProjectedScaledAmplitudeResult *result);
#endif

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

MVMCPfaffianStatus mvmc_projected_amplitude_values(
    const MVMCAbsolutePfaffianValueResult *components,
    const double complex *weights, size_t component_count,
    MVMCProjectedAmplitudeResult *result);

MVMCPfaffianStatus mvmc_projected_amplitude_value_slice(
    const MVMCAbsolutePfaffianValueResult *local_components,
    size_t local_component_count, const double complex *global_weights,
    int qp_total, int qp_start, int qp_end,
    MVMCProjectedAmplitudeResult *result);

/* False when this translation unit was compiled with fast-math semantics. */
int mvmc_absolute_pfaffian_strict_fp_enabled(void);

const char *mvmc_pfaffian_state_name(MVMCPfaffianState state);
const char *mvmc_pfaffian_value_state_name(MVMCPfaffianValueState state);

#endif /* MVMC_ABSOLUTE_PFAFFIAN_H */
