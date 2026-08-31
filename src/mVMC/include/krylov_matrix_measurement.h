/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_KRYLOV_MATRIX_MEASUREMENT_H
#define MVMC_KRYLOV_MATRIX_MEASUREMENT_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "absolute_pfaffian.h"
#include "krylov_streaming_statistics.h"

typedef enum {
  MVMC_KRYLOV_MATRIX_OVERLAP = 0,
  MVMC_KRYLOV_MATRIX_HAMILTONIAN,
  MVMC_KRYLOV_MATRIX_HAMILTONIAN_SQUARED,
  MVMC_KRYLOV_MATRIX_HAMILTONIAN_ADJOINT
} MVMCKrylovMatrixKind;

typedef struct {
  int order;
  double eta;
  double guide_lambda[MVMC_KRYLOV_MAX_ORDER + 1];
  double target_weight[MVMC_KRYLOV_MAX_ORDER + 1];
  double log_basis_scale[MVMC_KRYLOV_MAX_ORDER + 1];
} MVMCKrylovMatrixMeasurementPolicy;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int order;
  size_t dimension;
  size_t upper_count;
  double log_guide;
  double denominator;
  double target_over_guide;
  double guide_signal_over_guide;
  double minimum_log_contribution;
  double maximum_log_contribution;
  int finite_component_count;
  int zero_component_count;
  int zero_target_sample;
} MVMCKrylovMatrixMeasurementSampleDiagnostics;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int order;
  size_t dimension;
  size_t upper_count;
  uint64_t sample_count;
  uint64_t zero_target_sample_count;
  double denominator_sum;
  double denominator_mean;
  double denominator_sample_variance;
  double denominator_relative_se;
  double minimum_abs_denominator_mean;
  double maximum_denominator_relative_se;
  int denominator_stable;
  double effective_sample_count;
  double effective_sample_fraction;
  double zero_target_sample_fraction;
  double minimum_denominator;
  double maximum_denominator;
  double denominator_tail_ratio;
  double minimum_log_contribution;
  double maximum_log_contribution;
  double log_contribution_span;
  /* ||K - K^dagger||_F / ||K||_F, or the absolute numerator for zero norm. */
  double hamiltonian_antihermitian_residual;
  /* Full-matrix Frobenius norm ||K||_F. */
  double hamiltonian_norm;
  int hamiltonian_residual_available;
} MVMCKrylovMatrixMeasurementDiagnosticsSummary;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int order;
  size_t dimension;
  size_t upper_count;
  uint64_t sample_count;
  uint64_t zero_target_sample_count;
  MVMCKrylovStreamingMatrixAccumulator overlap;
  MVMCKrylovStreamingMatrixAccumulator hamiltonian;
  MVMCKrylovStreamingMatrixAccumulator hamiltonian_squared;
  int has_hamiltonian_adjoint;
  MVMCKrylovStreamingMatrixAccumulator hamiltonian_adjoint;
  MVMCKrylovStreamingRealSum denominator;
  MVMCKrylovStreamingRealSum denominator_squared;
  double minimum_denominator;
  double maximum_denominator;
  double minimum_log_contribution;
  double maximum_log_contribution;
} MVMCKrylovMatrixMeasurementAccumulator;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int order;
  size_t dimension;
  size_t upper_count;
  size_t block_count;
  uint64_t block_length;
  uint64_t sample_count;
  int has_hamiltonian_adjoint;
  MVMCKrylovMatrixMeasurementAccumulator *blocks;
} MVMCKrylovMatrixMeasurementBlockAccumulator;

MVMCKrylovStatus mvmc_krylov_matrix_measurement_dimension(
    int order, size_t *dimension, size_t *upper_count);

MVMCKrylovStatus mvmc_krylov_matrix_measurement_sample(
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    const MVMCScaledComplex *values, size_t value_count,
    double complex *overlap_upper, double complex *hamiltonian_upper,
    double complex *hamiltonian_squared_upper, size_t upper_count,
    MVMCKrylovMatrixMeasurementSampleDiagnostics *diagnostics);

MVMCKrylovStatus mvmc_krylov_matrix_measurement_sample_with_adjoint(
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    const MVMCScaledComplex *values, size_t value_count,
    double complex *overlap_upper, double complex *hamiltonian_upper,
    double complex *hamiltonian_adjoint_upper,
    double complex *hamiltonian_squared_upper, size_t upper_count,
    MVMCKrylovMatrixMeasurementSampleDiagnostics *diagnostics);

MVMCKrylovStatus mvmc_krylov_matrix_measurement_accumulator_init(
    int order, MVMCKrylovStreamingComplexSum *overlap_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_squared_entries,
    size_t upper_count, MVMCKrylovMatrixMeasurementAccumulator *accumulator);

MVMCKrylovStatus mvmc_krylov_matrix_measurement_accumulator_init_with_diagnostics(
    int order, MVMCKrylovStreamingComplexSum *overlap_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_adjoint_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_squared_entries,
    size_t upper_count, MVMCKrylovMatrixMeasurementAccumulator *accumulator);

MVMCKrylovStatus mvmc_krylov_matrix_measurement_accumulator_add_sample(
    MVMCKrylovMatrixMeasurementAccumulator *accumulator,
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    const MVMCScaledComplex *values, size_t value_count,
    MVMCKrylovMatrixMeasurementSampleDiagnostics *diagnostics);

MVMCKrylovStatus mvmc_krylov_matrix_measurement_accumulator_mean(
    const MVMCKrylovMatrixMeasurementAccumulator *accumulator,
    double complex *overlap_upper, double complex *hamiltonian_upper,
    double complex *hamiltonian_squared_upper, size_t upper_count,
    double *denominator_mean);

MVMCKrylovStatus mvmc_krylov_matrix_measurement_extract_block(
    const MVMCKrylovMatrixMeasurementAccumulator *accumulator,
    MVMCKrylovMatrixKind matrix_kind, size_t row, size_t column,
    MVMCKrylovJackknifeBlock *block);

MVMCKrylovStatus mvmc_krylov_matrix_measurement_block_accumulator_init(
    int order, size_t block_count, uint64_t block_length,
    MVMCKrylovMatrixMeasurementAccumulator *blocks,
    MVMCKrylovStreamingComplexSum *overlap_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_adjoint_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_squared_entries,
    size_t storage_entry_count,
    MVMCKrylovMatrixMeasurementBlockAccumulator *accumulator);

MVMCKrylovStatus mvmc_krylov_matrix_measurement_block_accumulator_add_sample(
    MVMCKrylovMatrixMeasurementBlockAccumulator *accumulator,
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    const MVMCScaledComplex *values, size_t value_count,
    MVMCKrylovMatrixMeasurementSampleDiagnostics *diagnostics);

MVMCKrylovStatus mvmc_krylov_matrix_measurement_block_accumulator_extract_blocks(
    const MVMCKrylovMatrixMeasurementBlockAccumulator *accumulator,
    MVMCKrylovMatrixKind matrix_kind, size_t row, size_t column,
    MVMCKrylovJackknifeBlock *blocks, size_t block_count);

MVMCKrylovStatus mvmc_krylov_matrix_measurement_diagnostics_summary(
    const MVMCKrylovMatrixMeasurementAccumulator *accumulator,
    double minimum_abs_denominator_mean,
    double maximum_denominator_relative_se,
    MVMCKrylovMatrixMeasurementDiagnosticsSummary *summary);

#endif /* MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE */

#endif /* MVMC_KRYLOV_MATRIX_MEASUREMENT_H */
