/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_KRYLOV_STREAMING_STATISTICS_H
#define MVMC_KRYLOV_STREAMING_STATISTICS_H

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)

#include "krylov_fock_reference.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS 16
#define MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCK_LENGTH 256
#define MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS 32
#define MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCK_LENGTH 128

typedef struct {
  double sum;
  double compensation;
} MVMCKrylovStreamingRealSum;

typedef struct {
  MVMCKrylovStreamingRealSum real;
  MVMCKrylovStreamingRealSum imag;
} MVMCKrylovStreamingComplexSum;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  size_t dimension;
  size_t upper_count;
  uint64_t sample_count;
  MVMCKrylovStreamingComplexSum *entries;
} MVMCKrylovStreamingMatrixAccumulator;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  size_t dimension;
  size_t upper_count;
  uint64_t sample_count;
  size_t storage_bytes;
} MVMCKrylovStreamingMatrixSummary;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  uint64_t sample_count;
  MVMCKrylovStreamingRealSum denominator;
  MVMCKrylovStreamingComplexSum numerator;
  MVMCKrylovStreamingRealSum denominator_squared;
  MVMCKrylovStreamingRealSum numerator_real_squared;
  MVMCKrylovStreamingRealSum numerator_imag_squared;
  MVMCKrylovStreamingRealSum numerator_real_imag;
  MVMCKrylovStreamingRealSum denominator_numerator_real;
  MVMCKrylovStreamingRealSum denominator_numerator_imag;
} MVMCKrylovStreamingRatioAccumulator;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  uint64_t sample_count;
  double complex estimate;
  double complex numerator_sum;
  double denominator_sum;
  double denominator_mean;
  double denominator_sample_variance;
  double denominator_relative_se;
  int denominator_stable;
  double variance_real;
  double variance_imag;
  double covariance_real_imag;
  double se_real;
  double se_imag;
  double complex_se;
} MVMCKrylovStreamingRatioResult;

typedef struct {
  double complex numerator;
  double denominator;
  uint64_t sample_count;
} MVMCKrylovJackknifeBlock;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  size_t block_count;
  uint64_t sample_count;
  double complex estimate;
  double complex jackknife_mean;
  double denominator_total;
  double minimum_leave_one_denominator;
  size_t unstable_leave_one_blocks;
  int denominator_stable;
  double variance_real;
  double variance_imag;
  double covariance_real_imag;
  double se_real;
  double se_imag;
  double complex_se;
} MVMCKrylovJackknifeResult;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  size_t official_block_count;
  size_t diagnostic_block_count;
  double official_complex_se;
  double diagnostic_complex_se;
  double symmetric_se_ratio;
  double maximum_symmetric_se_ratio;
  int passed;
} MVMCKrylovBlockStabilityResult;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  size_t sample_count;
  size_t positive_pair_count;
  double mean;
  double variance;
  double tau_int;
} MVMCKrylovTauIntResult;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  double tau_int;
  double block_length;
  double maximum_tau_int;
  double block_length_multiplier;
  double minimum_block_length;
  double actual_block_multiplier;
  int passed;
} MVMCKrylovTauIntGateResult;

MVMCKrylovStatus mvmc_krylov_streaming_upper_count(
    size_t dimension, size_t *upper_count);

MVMCKrylovStatus mvmc_krylov_streaming_upper_index(
    size_t dimension, size_t row, size_t column, size_t *index);

MVMCKrylovStatus mvmc_krylov_streaming_matrix_accumulator_init(
    size_t dimension, MVMCKrylovStreamingComplexSum *entries,
    size_t entry_count, MVMCKrylovStreamingMatrixAccumulator *accumulator);

MVMCKrylovStatus mvmc_krylov_streaming_matrix_accumulator_add_sample(
    MVMCKrylovStreamingMatrixAccumulator *accumulator,
    const double complex *upper_values, size_t upper_count, double weight);

MVMCKrylovStatus mvmc_krylov_streaming_matrix_accumulator_merge(
    MVMCKrylovStreamingMatrixAccumulator *destination,
    const MVMCKrylovStreamingMatrixAccumulator *source);

MVMCKrylovStatus mvmc_krylov_streaming_matrix_accumulator_mean(
    const MVMCKrylovStreamingMatrixAccumulator *accumulator,
    double complex *upper_mean, size_t upper_count);

MVMCKrylovStatus mvmc_krylov_streaming_matrix_accumulator_summary(
    const MVMCKrylovStreamingMatrixAccumulator *accumulator,
    MVMCKrylovStreamingMatrixSummary *summary);

MVMCKrylovStatus mvmc_krylov_streaming_ratio_accumulator_reset(
    MVMCKrylovStreamingRatioAccumulator *accumulator);

MVMCKrylovStatus mvmc_krylov_streaming_ratio_accumulator_add_sample(
    MVMCKrylovStreamingRatioAccumulator *accumulator,
    double complex numerator, double denominator);

MVMCKrylovStatus mvmc_krylov_streaming_ratio_accumulator_finalize(
    const MVMCKrylovStreamingRatioAccumulator *accumulator,
    double minimum_abs_denominator_mean,
    double maximum_denominator_relative_se,
    MVMCKrylovStreamingRatioResult *result);

MVMCKrylovStatus mvmc_krylov_jackknife_ratio(
    const MVMCKrylovJackknifeBlock *blocks, size_t block_count,
    double minimum_leave_one_denominator,
    MVMCKrylovJackknifeResult *result);

MVMCKrylovStatus mvmc_krylov_block_stability_check(
    const MVMCKrylovJackknifeResult *official,
    const MVMCKrylovJackknifeResult *diagnostic,
    double maximum_symmetric_se_ratio,
    MVMCKrylovBlockStabilityResult *result);

MVMCKrylovStatus mvmc_krylov_tau_int_geyer_initial_positive(
    const double *values, size_t sample_count,
    MVMCKrylovTauIntResult *result);

MVMCKrylovStatus mvmc_krylov_tau_int_gate_check(
    double tau_int, double block_length, double maximum_tau_int,
    double block_length_multiplier,
    MVMCKrylovTauIntGateResult *result);

#endif /* MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE */

#endif /* MVMC_KRYLOV_STREAMING_STATISTICS_H */
