/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_KRYLOV_FINAL_STATE_SENSITIVITY_H
#define MVMC_KRYLOV_FINAL_STATE_SENSITIVITY_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_P5_TESTING)

#include "krylov_final_state.h"

#define MVMC_KRYLOV_FINAL_SENSITIVITY_MODEL_VERSION UINT64_C(1)
#define MVMC_KRYLOV_FINAL_SENSITIVITY_MAX_UPPER                         \
  (((MVMC_KRYLOV_MAX_ORDER + 1) * (MVMC_KRYLOV_MAX_ORDER + 2)) / 2)

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int dimension;
  size_t upper_count;
  int has_observable;
  uint64_t model_hash;
  double complex overlap_upper[MVMC_KRYLOV_FINAL_SENSITIVITY_MAX_UPPER];
  double complex
      hamiltonian_upper[MVMC_KRYLOV_FINAL_SENSITIVITY_MAX_UPPER];
  double complex
      hamiltonian_squared_upper[MVMC_KRYLOV_FINAL_SENSITIVITY_MAX_UPPER];
  double complex observable_upper[MVMC_KRYLOV_FINAL_SENSITIVITY_MAX_UPPER];
} MVMCKrylovFinalStateQuadraticModel;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int dimension;
  int has_observable;
  double norm;
  double energy;
  double hamiltonian_second_moment;
  double full_support_variance;
  double observable;
} MVMCKrylovFinalStateQuadraticMetrics;

typedef struct {
  double energy_absolute;
  double energy_scaled;
  double variance_absolute;
  double variance_scaled;
  double observable_absolute;
  double observable_scaled;
} MVMCKrylovFinalStateDownstreamError;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int dimension;
  int has_observable;
  double coefficient_projective_distance;
  double overlap_relative_frobenius_error;
  double hamiltonian_relative_frobenius_error;
  double hamiltonian_squared_relative_frobenius_error;
  double observable_relative_frobenius_error;
  MVMCKrylovFinalStateQuadraticMetrics baseline;
  MVMCKrylovFinalStateQuadraticMetrics coefficient_only_metrics;
  MVMCKrylovFinalStateQuadraticMetrics matrix_only_metrics;
  MVMCKrylovFinalStateQuadraticMetrics combined_metrics;
  MVMCKrylovFinalStateDownstreamError coefficient_only;
  MVMCKrylovFinalStateDownstreamError matrix_only;
  MVMCKrylovFinalStateDownstreamError combined;
} MVMCKrylovFinalStateSensitivityComparison;

MVMCKrylovStatus mvmc_krylov_final_state_quadratic_model_create(
    int dimension, const double complex *overlap_upper,
    const double complex *hamiltonian_upper,
    const double complex *hamiltonian_squared_upper,
    const double complex *observable_upper, size_t upper_count,
    MVMCKrylovFinalStateQuadraticModel *model);

MVMCKrylovStatus mvmc_krylov_final_state_quadratic_metrics(
    const MVMCKrylovFinalStateQuadraticModel *model,
    const double complex *coefficient, size_t coefficient_count,
    MVMCKrylovFinalStateQuadraticMetrics *metrics);

MVMCKrylovStatus mvmc_krylov_final_state_sensitivity_compare(
    const MVMCKrylovFinalStateQuadraticModel *baseline_model,
    const double complex *baseline_coefficient,
    const MVMCKrylovFinalStateQuadraticModel *candidate_model,
    const double complex *candidate_coefficient,
    size_t coefficient_count,
    MVMCKrylovFinalStateSensitivityComparison *comparison);

#endif /* MVMC_ENABLE_POWER_LANCZOS_P5_TESTING */

#endif /* MVMC_KRYLOV_FINAL_STATE_SENSITIVITY_H */
