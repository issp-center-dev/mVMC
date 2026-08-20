/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "krylov_streaming_statistics.h"

#if !defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) &&                        \
    !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "krylov_streaming_statistics.c requires the power-Lanczos core"
#endif

#include <float.h>
#include <math.h>
#include <string.h>

static void real_sum_reset(MVMCKrylovStreamingRealSum *sum) {
  sum->sum = 0.0;
  sum->compensation = 0.0;
}

static MVMCKrylovStatus real_sum_add(
    MVMCKrylovStreamingRealSum *sum, double value) {
  double next;
  if (sum == NULL || !isfinite(value)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  next = sum->sum + value;
  if (!isfinite(next)) return MVMC_KRYLOV_STATUS_NONFINITE;
  if (fabs(sum->sum) >= fabs(value)) {
    sum->compensation += (sum->sum - next) + value;
  } else {
    sum->compensation += (value - next) + sum->sum;
  }
  if (!isfinite(sum->compensation)) return MVMC_KRYLOV_STATUS_NONFINITE;
  sum->sum = next;
  return MVMC_KRYLOV_STATUS_OK;
}

static double real_sum_value(const MVMCKrylovStreamingRealSum *sum) {
  return sum->sum + sum->compensation;
}

static void complex_sum_reset(MVMCKrylovStreamingComplexSum *sum) {
  real_sum_reset(&sum->real);
  real_sum_reset(&sum->imag);
}

static MVMCKrylovStatus complex_sum_add(
    MVMCKrylovStreamingComplexSum *sum, double complex value) {
  MVMCKrylovStatus status;
  if (sum == NULL || !isfinite(creal(value)) || !isfinite(cimag(value))) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = real_sum_add(&sum->real, creal(value));
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  return real_sum_add(&sum->imag, cimag(value));
}

static double complex complex_sum_value(
    const MVMCKrylovStreamingComplexSum *sum) {
  return real_sum_value(&sum->real) + I * real_sum_value(&sum->imag);
}

static int finite_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static void invalidate_matrix_accumulator(
    MVMCKrylovStatus status,
    MVMCKrylovStreamingMatrixAccumulator *accumulator) {
  if (accumulator == NULL) return;
  memset(accumulator, 0, sizeof(*accumulator));
  accumulator->status = status;
}

static void invalidate_matrix_summary(
    MVMCKrylovStatus status, MVMCKrylovStreamingMatrixSummary *summary) {
  if (summary == NULL) return;
  memset(summary, 0, sizeof(*summary));
  summary->status = status;
}

static void invalidate_ratio_result(
    MVMCKrylovStatus status, MVMCKrylovStreamingRatioResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
}

static void invalidate_jackknife_result(
    MVMCKrylovStatus status, MVMCKrylovJackknifeResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
}

static void invalidate_block_stability_result(
    MVMCKrylovStatus status, MVMCKrylovBlockStabilityResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
}

static void invalidate_tau_int_result(
    MVMCKrylovStatus status, MVMCKrylovTauIntResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
}

static void invalidate_tau_int_gate_result(
    MVMCKrylovStatus status, MVMCKrylovTauIntGateResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
}

static int sample_count_to_double(uint64_t sample_count, double *value) {
  if (value == NULL || sample_count == 0) return 0;
  *value = (double)sample_count;
  return isfinite(*value) && *value > 0.0;
}

static int clamp_nonnegative(double *value, double scale) {
  const double tolerance = 1.0e-10 * fmax(1.0, scale);
  if (value == NULL || !isfinite(*value)) return 0;
  if (*value >= 0.0) return 1;
  if (*value >= -tolerance) {
    *value = 0.0;
    return 1;
  }
  return 0;
}

MVMCKrylovStatus mvmc_krylov_streaming_upper_count(
    size_t dimension, size_t *upper_count) {
  size_t left;
  size_t right;
  if (upper_count == NULL || dimension == 0 || dimension == SIZE_MAX) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  left = dimension;
  right = dimension + 1;
  if ((left & (size_t)1) == 0) {
    left /= 2;
  } else {
    right /= 2;
  }
  if (left != 0 && right > SIZE_MAX / left) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  *upper_count = left * right;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_streaming_upper_index(
    size_t dimension, size_t row, size_t column, size_t *index) {
  size_t total = 0;
  size_t remaining = 0;
  size_t upper_row;
  size_t upper_column;
  MVMCKrylovStatus status;

  if (index == NULL || row >= dimension || column >= dimension) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  upper_row = row <= column ? row : column;
  upper_column = row <= column ? column : row;
  status = mvmc_krylov_streaming_upper_count(dimension, &total);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  status = mvmc_krylov_streaming_upper_count(
      dimension - upper_row, &remaining);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  *index = total - remaining + (upper_column - upper_row);
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_streaming_matrix_accumulator_init(
    size_t dimension, MVMCKrylovStreamingComplexSum *entries,
    size_t entry_count, MVMCKrylovStreamingMatrixAccumulator *accumulator) {
  size_t upper_count = 0;
  size_t entry;
  MVMCKrylovStatus status;

  if (accumulator == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_matrix_accumulator(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                                accumulator);
  status = mvmc_krylov_streaming_upper_count(dimension, &upper_count);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_matrix_accumulator(status, accumulator);
    return status;
  }
  if (entries == NULL || entry_count != upper_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (entry = 0; entry < entry_count; ++entry) {
    complex_sum_reset(&entries[entry]);
  }
  accumulator->valid = 1;
  accumulator->status = MVMC_KRYLOV_STATUS_OK;
  accumulator->dimension = dimension;
  accumulator->upper_count = upper_count;
  accumulator->sample_count = 0;
  accumulator->entries = entries;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_streaming_matrix_accumulator_add_sample(
    MVMCKrylovStreamingMatrixAccumulator *accumulator,
    const double complex *upper_values, size_t upper_count, double weight) {
  size_t entry;
  MVMCKrylovStatus status;

  if (accumulator == NULL || !accumulator->valid ||
      accumulator->status != MVMC_KRYLOV_STATUS_OK ||
      accumulator->entries == NULL || upper_values == NULL ||
      upper_count != accumulator->upper_count || !isfinite(weight)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (accumulator->sample_count == UINT64_MAX) {
    accumulator->status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  for (entry = 0; entry < upper_count; ++entry) {
    if (!finite_complex(upper_values[entry]) ||
        !finite_complex(upper_values[entry] * weight)) {
      accumulator->status = MVMC_KRYLOV_STATUS_NONFINITE;
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  for (entry = 0; entry < upper_count; ++entry) {
    status = complex_sum_add(&accumulator->entries[entry],
                             upper_values[entry] * weight);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      accumulator->status = status;
      return status;
    }
  }
  ++accumulator->sample_count;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_streaming_matrix_accumulator_merge(
    MVMCKrylovStreamingMatrixAccumulator *destination,
    const MVMCKrylovStreamingMatrixAccumulator *source) {
  size_t entry;
  MVMCKrylovStatus status;

  if (destination == NULL || source == NULL || !destination->valid ||
      !source->valid || destination->status != MVMC_KRYLOV_STATUS_OK ||
      source->status != MVMC_KRYLOV_STATUS_OK ||
      destination->dimension != source->dimension ||
      destination->upper_count != source->upper_count ||
      destination->entries == NULL || source->entries == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (source->sample_count > UINT64_MAX - destination->sample_count) {
    destination->status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  for (entry = 0; entry < destination->upper_count; ++entry) {
    status = real_sum_add(&destination->entries[entry].real,
                          source->entries[entry].real.sum);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = real_sum_add(&destination->entries[entry].real,
                            source->entries[entry].real.compensation);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = real_sum_add(&destination->entries[entry].imag,
                            source->entries[entry].imag.sum);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = real_sum_add(&destination->entries[entry].imag,
                            source->entries[entry].imag.compensation);
    }
    if (status != MVMC_KRYLOV_STATUS_OK) {
      destination->status = status;
      return status;
    }
  }
  destination->sample_count += source->sample_count;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_streaming_matrix_accumulator_mean(
    const MVMCKrylovStreamingMatrixAccumulator *accumulator,
    double complex *upper_mean, size_t upper_count) {
  double sample_count;
  size_t entry;
  if (accumulator == NULL || !accumulator->valid ||
      accumulator->status != MVMC_KRYLOV_STATUS_OK ||
      accumulator->entries == NULL || upper_mean == NULL ||
      upper_count != accumulator->upper_count ||
      !sample_count_to_double(accumulator->sample_count, &sample_count)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (entry = 0; entry < upper_count; ++entry) {
    upper_mean[entry] = complex_sum_value(&accumulator->entries[entry]) /
                        sample_count;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_streaming_matrix_accumulator_summary(
    const MVMCKrylovStreamingMatrixAccumulator *accumulator,
    MVMCKrylovStreamingMatrixSummary *summary) {
  MVMCKrylovStreamingMatrixSummary candidate;
  if (summary == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_matrix_summary(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, summary);
  if (accumulator == NULL || !accumulator->valid ||
      accumulator->status != MVMC_KRYLOV_STATUS_OK ||
      accumulator->entries == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (accumulator->upper_count > SIZE_MAX /
                                      sizeof(MVMCKrylovStreamingComplexSum)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.dimension = accumulator->dimension;
  candidate.upper_count = accumulator->upper_count;
  candidate.sample_count = accumulator->sample_count;
  candidate.storage_bytes = accumulator->upper_count *
                            sizeof(MVMCKrylovStreamingComplexSum);
  *summary = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_streaming_ratio_accumulator_reset(
    MVMCKrylovStreamingRatioAccumulator *accumulator) {
  if (accumulator == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(accumulator, 0, sizeof(*accumulator));
  accumulator->valid = 1;
  accumulator->status = MVMC_KRYLOV_STATUS_OK;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_streaming_ratio_accumulator_add_sample(
    MVMCKrylovStreamingRatioAccumulator *accumulator,
    double complex numerator, double denominator) {
  const double numerator_real = creal(numerator);
  const double numerator_imag = cimag(numerator);
  const double denominator_squared = denominator * denominator;
  const double numerator_real_squared = numerator_real * numerator_real;
  const double numerator_imag_squared = numerator_imag * numerator_imag;
  const double numerator_real_imag = numerator_real * numerator_imag;
  const double denominator_numerator_real = denominator * numerator_real;
  const double denominator_numerator_imag = denominator * numerator_imag;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;

  if (accumulator == NULL || !accumulator->valid ||
      accumulator->status != MVMC_KRYLOV_STATUS_OK ||
      !finite_complex(numerator) || !isfinite(denominator) ||
      denominator <= 0.0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (accumulator->sample_count == UINT64_MAX) {
    accumulator->status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  if (!isfinite(denominator_squared) || !isfinite(numerator_real_squared) ||
      !isfinite(numerator_imag_squared) || !isfinite(numerator_real_imag) ||
      !isfinite(denominator_numerator_real) ||
      !isfinite(denominator_numerator_imag)) {
    accumulator->status = MVMC_KRYLOV_STATUS_NONFINITE;
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }

  status = real_sum_add(&accumulator->denominator, denominator);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = complex_sum_add(&accumulator->numerator, numerator);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = real_sum_add(&accumulator->denominator_squared,
                          denominator_squared);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = real_sum_add(&accumulator->numerator_real_squared,
                          numerator_real_squared);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = real_sum_add(&accumulator->numerator_imag_squared,
                          numerator_imag_squared);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = real_sum_add(&accumulator->numerator_real_imag,
                          numerator_real_imag);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = real_sum_add(&accumulator->denominator_numerator_real,
                          denominator_numerator_real);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = real_sum_add(&accumulator->denominator_numerator_imag,
                          denominator_numerator_imag);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    accumulator->status = status;
    return status;
  }
  ++accumulator->sample_count;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_streaming_ratio_accumulator_finalize(
    const MVMCKrylovStreamingRatioAccumulator *accumulator,
    double minimum_abs_denominator_mean,
    double maximum_denominator_relative_se,
    MVMCKrylovStreamingRatioResult *result) {
  MVMCKrylovStreamingRatioResult candidate;
  double sample_count;
  double denominator_sum;
  double denominator_mean;
  double complex numerator_sum;
  double complex estimate;
  double denominator_ss;
  double denominator_scale;
  double denominator_sample_variance = 0.0;
  double denominator_relative_se = 0.0;
  double residual_real_ss = 0.0;
  double residual_imag_ss = 0.0;
  double residual_real_imag_ss = 0.0;
  double residual_real_scale = 1.0;
  double residual_imag_scale = 1.0;
  double denominator_squared;
  double numerator_real_squared;
  double numerator_imag_squared;
  double numerator_real_imag;
  double denominator_numerator_real;
  double denominator_numerator_imag;
  double denominator_mean_squared;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_ratio_result(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (accumulator == NULL || !accumulator->valid ||
      accumulator->status != MVMC_KRYLOV_STATUS_OK ||
      !sample_count_to_double(accumulator->sample_count, &sample_count) ||
      !isfinite(minimum_abs_denominator_mean) ||
      minimum_abs_denominator_mean < 0.0 ||
      !isfinite(maximum_denominator_relative_se) ||
      maximum_denominator_relative_se < 0.0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }

  denominator_sum = real_sum_value(&accumulator->denominator);
  numerator_sum = complex_sum_value(&accumulator->numerator);
  denominator_squared = real_sum_value(&accumulator->denominator_squared);
  numerator_real_squared =
      real_sum_value(&accumulator->numerator_real_squared);
  numerator_imag_squared =
      real_sum_value(&accumulator->numerator_imag_squared);
  numerator_real_imag = real_sum_value(&accumulator->numerator_real_imag);
  denominator_numerator_real =
      real_sum_value(&accumulator->denominator_numerator_real);
  denominator_numerator_imag =
      real_sum_value(&accumulator->denominator_numerator_imag);
  if (!isfinite(denominator_sum) || denominator_sum <= 0.0 ||
      !finite_complex(numerator_sum) || !isfinite(denominator_squared) ||
      !isfinite(numerator_real_squared) ||
      !isfinite(numerator_imag_squared) ||
      !isfinite(numerator_real_imag) ||
      !isfinite(denominator_numerator_real) ||
      !isfinite(denominator_numerator_imag)) {
    invalidate_ratio_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }

  denominator_mean = denominator_sum / sample_count;
  denominator_mean_squared = denominator_mean * denominator_mean;
  estimate = numerator_sum / denominator_sum;
  if (!isfinite(denominator_mean) || denominator_mean <= 0.0 ||
      !finite_complex(estimate) || !isfinite(denominator_mean_squared) ||
      denominator_mean_squared <= 0.0) {
    invalidate_ratio_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }

  denominator_ss = denominator_squared - 2.0 * denominator_mean *
                                             denominator_sum +
                   sample_count * denominator_mean_squared;
  denominator_scale = fabs(denominator_squared) +
                      fabs(sample_count * denominator_mean_squared);
  if (!clamp_nonnegative(&denominator_ss, denominator_scale)) {
    invalidate_ratio_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }

  if (accumulator->sample_count > 1) {
    const double denominator_estimator_variance_denominator =
        (sample_count - 1.0) * sample_count;
    denominator_sample_variance = denominator_ss / (sample_count - 1.0);
    denominator_relative_se =
        sqrt(denominator_ss / denominator_estimator_variance_denominator) /
        fabs(denominator_mean);

    residual_real_ss =
        numerator_real_squared -
        2.0 * creal(estimate) * denominator_numerator_real +
        creal(estimate) * creal(estimate) * denominator_squared;
    residual_imag_ss =
        numerator_imag_squared -
        2.0 * cimag(estimate) * denominator_numerator_imag +
        cimag(estimate) * cimag(estimate) * denominator_squared;
    residual_real_imag_ss =
        numerator_real_imag -
        creal(estimate) * denominator_numerator_imag -
        cimag(estimate) * denominator_numerator_real +
        creal(estimate) * cimag(estimate) * denominator_squared;
    residual_real_scale =
        fabs(numerator_real_squared) +
        fabs(creal(estimate) * creal(estimate) * denominator_squared) +
        fabs(2.0 * creal(estimate) * denominator_numerator_real);
    residual_imag_scale =
        fabs(numerator_imag_squared) +
        fabs(cimag(estimate) * cimag(estimate) * denominator_squared) +
        fabs(2.0 * cimag(estimate) * denominator_numerator_imag);
    if (!clamp_nonnegative(&residual_real_ss, residual_real_scale) ||
        !clamp_nonnegative(&residual_imag_ss, residual_imag_scale) ||
        !isfinite(residual_real_imag_ss)) {
      invalidate_ratio_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }

  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.sample_count = accumulator->sample_count;
  candidate.estimate = estimate;
  candidate.numerator_sum = numerator_sum;
  candidate.denominator_sum = denominator_sum;
  candidate.denominator_mean = denominator_mean;
  candidate.denominator_sample_variance = denominator_sample_variance;
  candidate.denominator_relative_se = denominator_relative_se;
  candidate.denominator_stable =
      fabs(denominator_mean) >= minimum_abs_denominator_mean &&
      denominator_relative_se <= maximum_denominator_relative_se;
  if (accumulator->sample_count > 1) {
    const double estimator_denominator =
        (sample_count - 1.0) * sample_count * denominator_mean_squared;
    candidate.variance_real = residual_real_ss / estimator_denominator;
    candidate.variance_imag = residual_imag_ss / estimator_denominator;
    candidate.covariance_real_imag =
        residual_real_imag_ss / estimator_denominator;
    candidate.se_real = sqrt(candidate.variance_real);
    candidate.se_imag = sqrt(candidate.variance_imag);
    candidate.complex_se =
        sqrt(fmax(0.0, candidate.variance_real + candidate.variance_imag));
  }
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_jackknife_ratio(
    const MVMCKrylovJackknifeBlock *blocks, size_t block_count,
    double minimum_leave_one_denominator,
    MVMCKrylovJackknifeResult *result) {
  MVMCKrylovJackknifeResult candidate;
  double complex total_numerator = 0.0;
  double denominator_total = 0.0;
  uint64_t sample_count = 0;
  double complex leave_one_sum = 0.0;
  double complex jackknife_mean;
  double variance_real = 0.0;
  double variance_imag = 0.0;
  double covariance_real_imag = 0.0;
  size_t block;
  size_t unstable_leave_one_blocks = 0;
  double block_count_double;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_jackknife_result(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (blocks == NULL || block_count < 2 ||
      !isfinite(minimum_leave_one_denominator) ||
      minimum_leave_one_denominator < 0.0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  block_count_double = (double)block_count;
  if (!isfinite(block_count_double) || block_count_double <= 0.0) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }

  for (block = 0; block < block_count; ++block) {
    if (!finite_complex(blocks[block].numerator) ||
        !isfinite(blocks[block].denominator) ||
        blocks[block].denominator < 0.0 ||
        (blocks[block].denominator == 0.0 &&
         cabs(blocks[block].numerator) != 0.0) ||
        blocks[block].sample_count == 0 ||
        blocks[block].sample_count > UINT64_MAX - sample_count) {
      invalidate_jackknife_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    total_numerator += blocks[block].numerator;
    denominator_total += blocks[block].denominator;
    sample_count += blocks[block].sample_count;
  }
  if (!finite_complex(total_numerator) || !isfinite(denominator_total) ||
      denominator_total <= 0.0) {
    invalidate_jackknife_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }

  for (block = 0; block < block_count; ++block) {
    const double complex leave_numerator =
        total_numerator - blocks[block].numerator;
    const double leave_denominator =
        denominator_total - blocks[block].denominator;
    if (!isfinite(leave_denominator) || leave_denominator <= 0.0 ||
        !finite_complex(leave_numerator)) {
      invalidate_jackknife_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    if (leave_denominator < minimum_leave_one_denominator) {
      ++unstable_leave_one_blocks;
    }
    leave_one_sum += leave_numerator / leave_denominator;
  }
  jackknife_mean = leave_one_sum / block_count_double;
  if (!finite_complex(jackknife_mean)) {
    invalidate_jackknife_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }

  for (block = 0; block < block_count; ++block) {
    const double complex leave_numerator =
        total_numerator - blocks[block].numerator;
    const double leave_denominator =
        denominator_total - blocks[block].denominator;
    const double complex delta =
        leave_numerator / leave_denominator - jackknife_mean;
    variance_real += creal(delta) * creal(delta);
    variance_imag += cimag(delta) * cimag(delta);
    covariance_real_imag += creal(delta) * cimag(delta);
  }
  variance_real *= (block_count_double - 1.0) / block_count_double;
  variance_imag *= (block_count_double - 1.0) / block_count_double;
  covariance_real_imag *= (block_count_double - 1.0) / block_count_double;
  if (!isfinite(variance_real) || !isfinite(variance_imag) ||
      !isfinite(covariance_real_imag) ||
      !clamp_nonnegative(&variance_real, variance_real) ||
      !clamp_nonnegative(&variance_imag, variance_imag)) {
    invalidate_jackknife_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }

  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.block_count = block_count;
  candidate.sample_count = sample_count;
  candidate.estimate = total_numerator / denominator_total;
  candidate.jackknife_mean = jackknife_mean;
  candidate.denominator_total = denominator_total;
  candidate.minimum_leave_one_denominator = minimum_leave_one_denominator;
  candidate.unstable_leave_one_blocks = unstable_leave_one_blocks;
  candidate.denominator_stable = unstable_leave_one_blocks == 0;
  candidate.variance_real = variance_real;
  candidate.variance_imag = variance_imag;
  candidate.covariance_real_imag = covariance_real_imag;
  candidate.se_real = sqrt(variance_real);
  candidate.se_imag = sqrt(variance_imag);
  candidate.complex_se = sqrt(fmax(0.0, variance_real + variance_imag));
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_block_stability_check(
    const MVMCKrylovJackknifeResult *official,
    const MVMCKrylovJackknifeResult *diagnostic,
    double maximum_symmetric_se_ratio,
    MVMCKrylovBlockStabilityResult *result) {
  MVMCKrylovBlockStabilityResult candidate;
  double smaller;
  double larger;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_block_stability_result(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                                    result);
  if (official == NULL || diagnostic == NULL || !official->valid ||
      !diagnostic->valid || official->status != MVMC_KRYLOV_STATUS_OK ||
      diagnostic->status != MVMC_KRYLOV_STATUS_OK ||
      !isfinite(official->complex_se) || official->complex_se < 0.0 ||
      !isfinite(diagnostic->complex_se) || diagnostic->complex_se < 0.0 ||
      !isfinite(maximum_symmetric_se_ratio) ||
      maximum_symmetric_se_ratio < 1.0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  smaller = fmin(official->complex_se, diagnostic->complex_se);
  larger = fmax(official->complex_se, diagnostic->complex_se);
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.official_block_count = official->block_count;
  candidate.diagnostic_block_count = diagnostic->block_count;
  candidate.official_complex_se = official->complex_se;
  candidate.diagnostic_complex_se = diagnostic->complex_se;
  candidate.maximum_symmetric_se_ratio = maximum_symmetric_se_ratio;
  if (larger == 0.0) {
    candidate.symmetric_se_ratio = 1.0;
  } else if (smaller == 0.0) {
    candidate.symmetric_se_ratio = DBL_MAX;
  } else {
    candidate.symmetric_se_ratio = larger / smaller;
  }
  candidate.passed =
      candidate.symmetric_se_ratio <= maximum_symmetric_se_ratio;
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

static double autocovariance_biased(
    const double *values, size_t sample_count, double mean, size_t lag) {
  double sum = 0.0;
  size_t index;
  for (index = 0; index + lag < sample_count; ++index) {
    sum += (values[index] - mean) * (values[index + lag] - mean);
  }
  return sum / (double)sample_count;
}

MVMCKrylovStatus mvmc_krylov_tau_int_geyer_initial_positive(
    const double *values, size_t sample_count,
    MVMCKrylovTauIntResult *result) {
  MVMCKrylovTauIntResult candidate;
  double mean = 0.0;
  double variance;
  double pair_sum = 0.0;
  double sample_count_double;
  size_t index;
  size_t lag;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_tau_int_result(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (values == NULL || sample_count == 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  sample_count_double = (double)sample_count;
  if (!isfinite(sample_count_double) || sample_count_double <= 0.0) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  for (index = 0; index < sample_count; ++index) {
    if (!isfinite(values[index])) {
      invalidate_tau_int_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    mean += values[index];
  }
  mean /= sample_count_double;
  variance = autocovariance_biased(values, sample_count, mean, 0);
  if (!isfinite(mean) || !isfinite(variance) || variance < 0.0) {
    invalidate_tau_int_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }

  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.sample_count = sample_count;
  candidate.mean = mean;
  candidate.variance = variance;
  candidate.tau_int = 1.0;
  if (variance == 0.0) {
    *result = candidate;
    return MVMC_KRYLOV_STATUS_OK;
  }

  for (lag = 0; lag < sample_count; lag += 2) {
    const double gamma_even =
        autocovariance_biased(values, sample_count, mean, lag);
    const double gamma_odd =
        lag + 1 < sample_count
            ? autocovariance_biased(values, sample_count, mean, lag + 1)
            : 0.0;
    const double pair = (gamma_even + gamma_odd) / variance;
    if (!isfinite(gamma_even) || !isfinite(gamma_odd) ||
        !isfinite(pair)) {
      invalidate_tau_int_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    if (pair <= 0.0) break;
    pair_sum += pair;
    ++candidate.positive_pair_count;
  }
  candidate.tau_int = -1.0 + 2.0 * pair_sum;
  if (!isfinite(candidate.tau_int)) {
    invalidate_tau_int_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  if (candidate.tau_int < 1.0) candidate.tau_int = 1.0;
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_tau_int_gate_check(
    double tau_int, double block_length, double maximum_tau_int,
    double block_length_multiplier,
    MVMCKrylovTauIntGateResult *result) {
  MVMCKrylovTauIntGateResult candidate;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_tau_int_gate_result(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                                 result);
  if (!isfinite(tau_int) || tau_int < 0.0 || !isfinite(block_length) ||
      block_length <= 0.0 || !isfinite(maximum_tau_int) ||
      maximum_tau_int <= 0.0 || !isfinite(block_length_multiplier) ||
      block_length_multiplier <= 0.0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.tau_int = tau_int;
  candidate.block_length = block_length;
  candidate.maximum_tau_int = maximum_tau_int;
  candidate.block_length_multiplier = block_length_multiplier;
  candidate.minimum_block_length = block_length_multiplier * tau_int;
  candidate.actual_block_multiplier =
      tau_int > 0.0 ? block_length / tau_int : DBL_MAX;
  candidate.passed = tau_int <= maximum_tau_int &&
                     block_length >= candidate.minimum_block_length;
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}
