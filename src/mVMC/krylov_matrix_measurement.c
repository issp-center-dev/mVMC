/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "krylov_matrix_measurement.h"

#if !defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                         \
    !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "krylov_matrix_measurement.c is Testing-only"
#endif

#include <float.h>
#include <math.h>
#include <string.h>

typedef struct {
  int zero;
  double complex phase;
  double log_abs;
} DimensionlessValue;

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

static int finite_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static void invalidate_sample_diagnostics(
    MVMCKrylovStatus status,
    MVMCKrylovMatrixMeasurementSampleDiagnostics *diagnostics) {
  if (diagnostics == NULL) return;
  memset(diagnostics, 0, sizeof(*diagnostics));
  diagnostics->status = status;
  diagnostics->log_guide = NAN;
  diagnostics->denominator = NAN;
  diagnostics->target_over_guide = NAN;
  diagnostics->guide_signal_over_guide = NAN;
  diagnostics->minimum_log_contribution = NAN;
  diagnostics->maximum_log_contribution = NAN;
}

static void invalidate_accumulator(
    MVMCKrylovStatus status,
    MVMCKrylovMatrixMeasurementAccumulator *accumulator) {
  if (accumulator == NULL) return;
  memset(accumulator, 0, sizeof(*accumulator));
  accumulator->status = status;
  accumulator->minimum_denominator = NAN;
  accumulator->maximum_denominator = NAN;
  accumulator->minimum_log_contribution = NAN;
  accumulator->maximum_log_contribution = NAN;
}

static void invalidate_block_accumulator(
    MVMCKrylovStatus status,
    MVMCKrylovMatrixMeasurementBlockAccumulator *accumulator) {
  if (accumulator == NULL) return;
  memset(accumulator, 0, sizeof(*accumulator));
  accumulator->status = status;
}

static void invalidate_diagnostics_summary(
    MVMCKrylovStatus status,
    MVMCKrylovMatrixMeasurementDiagnosticsSummary *summary) {
  if (summary == NULL) return;
  memset(summary, 0, sizeof(*summary));
  summary->status = status;
  summary->denominator_sum = NAN;
  summary->denominator_mean = NAN;
  summary->denominator_sample_variance = NAN;
  summary->denominator_relative_se = NAN;
  summary->minimum_abs_denominator_mean = NAN;
  summary->maximum_denominator_relative_se = NAN;
  summary->effective_sample_count = NAN;
  summary->effective_sample_fraction = NAN;
  summary->zero_target_sample_fraction = NAN;
  summary->minimum_denominator = NAN;
  summary->maximum_denominator = NAN;
  summary->denominator_tail_ratio = NAN;
  summary->minimum_log_contribution = NAN;
  summary->maximum_log_contribution = NAN;
  summary->log_contribution_span = NAN;
  summary->hamiltonian_antihermitian_residual = NAN;
  summary->hamiltonian_norm = NAN;
}

static MVMCKrylovStatus policy_validate(
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    size_t *dimension, size_t *upper_count) {
  int index;
  int has_target = 0;
  MVMCKrylovStatus status;
  if (policy == NULL || policy->order <= 0 ||
      policy->order > MVMC_KRYLOV_MAX_ORDER || !isfinite(policy->eta) ||
      policy->eta <= 0.0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_matrix_measurement_dimension(
      policy->order, dimension, upper_count);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  for (index = 0; index <= policy->order; ++index) {
    if (!isfinite(policy->guide_lambda[index]) ||
        policy->guide_lambda[index] < 0.0 ||
        !isfinite(policy->target_weight[index]) ||
        policy->target_weight[index] < 0.0 ||
        !isfinite(policy->log_basis_scale[index])) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    has_target = has_target || policy->target_weight[index] > 0.0;
  }
  /* eta is strictly positive, so an all-zero guide signal is the exact
     finite-sector uniform-guide boundary and remains a valid guide. */
  return has_target ? MVMC_KRYLOV_STATUS_OK
                    : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
}

static double log_add(double left, double right) {
  if (left == -INFINITY) return right;
  if (right == -INFINITY) return left;
  if (left < right) {
    const double temporary = left;
    left = right;
    right = temporary;
  }
  return left + log1p(exp(right - left));
}

static MVMCKrylovStatus safe_exp_value(double log_value, double *value) {
  if (value == NULL || !isfinite(log_value)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (log_value > log(DBL_MAX)) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  if (log_value < log(nextafter(0.0, 1.0))) {
    *value = 0.0;
    return MVMC_KRYLOV_STATUS_OK;
  }
  *value = exp(log_value);
  return isfinite(*value) ? MVMC_KRYLOV_STATUS_OK
                          : MVMC_KRYLOV_STATUS_NONFINITE;
}

static MVMCKrylovStatus scaled_to_dimensionless(
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    const MVMCScaledComplex *values, size_t value_count,
    DimensionlessValue *dimensionless, int *finite_count, int *zero_count) {
  int index;
  if (values == NULL || dimensionless == NULL ||
      finite_count == NULL || zero_count == NULL ||
      value_count < (size_t)policy->order + 1) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *finite_count = 0;
  *zero_count = 0;
  for (index = 0; index <= policy->order; ++index) {
    const MVMCScaledComplex *value = &values[index];
    memset(&dimensionless[index], 0, sizeof(dimensionless[index]));
    if (!mvmc_scaled_complex_is_valid(value)) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    switch (value->state) {
      case MVMC_SCALED_COMPLEX_FINITE_NONZERO:
        if (!finite_complex(value->phase) || !isfinite(value->log_abs) ||
            !isfinite(value->log_abs + policy->log_basis_scale[index])) {
          return MVMC_KRYLOV_STATUS_NONFINITE;
        }
        dimensionless[index].zero = 0;
        dimensionless[index].phase = value->phase;
        dimensionless[index].log_abs =
            value->log_abs + policy->log_basis_scale[index];
        ++(*finite_count);
        break;
      case MVMC_SCALED_COMPLEX_EXACT_ZERO:
      case MVMC_SCALED_COMPLEX_NUMERIC_ZERO:
        dimensionless[index].zero = 1;
        ++(*zero_count);
        break;
      case MVMC_SCALED_COMPLEX_NONFINITE:
        return MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus compute_log_guide(
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    const DimensionlessValue *dimensionless, double *log_guide) {
  int index;
  double result;
  if (log_guide == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  result = log(policy->eta);
  for (index = 0; index <= policy->order; ++index) {
    if (!dimensionless[index].zero && policy->guide_lambda[index] > 0.0) {
      const double term =
          log(policy->guide_lambda[index]) +
          2.0 * dimensionless[index].log_abs;
      if (!isfinite(term)) return MVMC_KRYLOV_STATUS_NONFINITE;
      result = log_add(result, term);
    }
  }
  if (!isfinite(result)) return MVMC_KRYLOV_STATUS_NONFINITE;
  *log_guide = result;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus weighted_norm_over_guide(
    const DimensionlessValue *dimensionless, int order,
    const double *weight, double log_guide, double *result) {
  int index;
  double sum = 0.0;
  if (result == NULL || dimensionless == NULL || weight == NULL ||
      !isfinite(log_guide)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (index = 0; index <= order; ++index) {
    double contribution = 0.0;
    MVMCKrylovStatus status;
    if (dimensionless[index].zero || weight[index] == 0.0) continue;
    status = safe_exp_value(
        log(weight[index]) + 2.0 * dimensionless[index].log_abs -
            log_guide,
        &contribution);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    sum += contribution;
    if (!isfinite(sum)) return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  *result = sum;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus guide_normalized_product(
    const DimensionlessValue *left, const DimensionlessValue *right,
    double log_positive_factor, double log_guide,
    double complex *result, double *log_contribution) {
  double log_magnitude;
  double magnitude = 0.0;
  MVMCKrylovStatus status;
  if (result == NULL || log_contribution == NULL || left == NULL ||
      right == NULL || !isfinite(log_positive_factor) ||
      !isfinite(log_guide)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (left->zero || right->zero) {
    *result = 0.0;
    *log_contribution = -INFINITY;
    return MVMC_KRYLOV_STATUS_OK;
  }
  log_magnitude =
      left->log_abs + right->log_abs + log_positive_factor - log_guide;
  if (!isfinite(log_magnitude)) return MVMC_KRYLOV_STATUS_NONFINITE;
  status = safe_exp_value(log_magnitude, &magnitude);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  *result = magnitude * conj(left->phase) * right->phase;
  if (!finite_complex(*result)) return MVMC_KRYLOV_STATUS_NONFINITE;
  *log_contribution = log_magnitude;
  return MVMC_KRYLOV_STATUS_OK;
}

static void note_log_contribution(
    double log_contribution,
    MVMCKrylovMatrixMeasurementSampleDiagnostics *diagnostics) {
  if (diagnostics == NULL || log_contribution == -INFINITY) return;
  if (!isfinite(diagnostics->minimum_log_contribution) ||
      log_contribution < diagnostics->minimum_log_contribution) {
    diagnostics->minimum_log_contribution = log_contribution;
  }
  if (!isfinite(diagnostics->maximum_log_contribution) ||
      log_contribution > diagnostics->maximum_log_contribution) {
    diagnostics->maximum_log_contribution = log_contribution;
  }
}

static int all_upper_entries_zero(
    const double complex *values, size_t upper_count) {
  size_t entry;
  if (values == NULL) return 1;
  for (entry = 0; entry < upper_count; ++entry) {
    if (creal(values[entry]) != 0.0 || cimag(values[entry]) != 0.0) {
      return 0;
    }
  }
  return 1;
}

MVMCKrylovStatus mvmc_krylov_matrix_measurement_dimension(
    int order, size_t *dimension, size_t *upper_count) {
  MVMCKrylovStatus status;
  if (dimension == NULL || upper_count == NULL || order <= 0 ||
      order > MVMC_KRYLOV_MAX_ORDER) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_streaming_upper_count((size_t)order, upper_count);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  *dimension = (size_t)order;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_matrix_measurement_sample(
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    const MVMCScaledComplex *values, size_t value_count,
    double complex *overlap_upper, double complex *hamiltonian_upper,
    double complex *hamiltonian_squared_upper, size_t upper_count,
    MVMCKrylovMatrixMeasurementSampleDiagnostics *diagnostics) {
  return mvmc_krylov_matrix_measurement_sample_with_adjoint(
      policy, values, value_count, overlap_upper, hamiltonian_upper, NULL,
      hamiltonian_squared_upper, upper_count, diagnostics);
}

MVMCKrylovStatus mvmc_krylov_matrix_measurement_sample_with_adjoint(
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    const MVMCScaledComplex *values, size_t value_count,
    double complex *overlap_upper, double complex *hamiltonian_upper,
    double complex *hamiltonian_adjoint_upper,
    double complex *hamiltonian_squared_upper, size_t upper_count,
    MVMCKrylovMatrixMeasurementSampleDiagnostics *diagnostics) {
  DimensionlessValue dimensionless[MVMC_KRYLOV_MAX_ORDER + 1];
  MVMCKrylovMatrixMeasurementSampleDiagnostics candidate;
  size_t dimension = 0;
  size_t required_upper_count = 0;
  size_t row;
  size_t column;
  double log_guide = NAN;
  MVMCKrylovStatus status;

  if (diagnostics == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_sample_diagnostics(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                                diagnostics);
  status = policy_validate(policy, &dimension, &required_upper_count);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_sample_diagnostics(status, diagnostics);
    return status;
  }
  if (overlap_upper == NULL || hamiltonian_upper == NULL ||
      hamiltonian_squared_upper == NULL ||
      upper_count != required_upper_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }

  memset(&candidate, 0, sizeof(candidate));
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.order = policy->order;
  candidate.dimension = dimension;
  candidate.upper_count = upper_count;
  candidate.minimum_log_contribution = NAN;
  candidate.maximum_log_contribution = NAN;

  status = scaled_to_dimensionless(
      policy, values, value_count, dimensionless,
      &candidate.finite_component_count, &candidate.zero_component_count);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_sample_diagnostics(status, diagnostics);
    return status;
  }
  status = compute_log_guide(policy, dimensionless, &log_guide);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_sample_diagnostics(status, diagnostics);
    return status;
  }
  candidate.log_guide = log_guide;
  status = weighted_norm_over_guide(
      dimensionless, policy->order, policy->target_weight, log_guide,
      &candidate.target_over_guide);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_sample_diagnostics(status, diagnostics);
    return status;
  }
  status = weighted_norm_over_guide(
      dimensionless, policy->order, policy->guide_lambda, log_guide,
      &candidate.guide_signal_over_guide);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_sample_diagnostics(status, diagnostics);
    return status;
  }
  candidate.denominator = candidate.target_over_guide;
  if (!isfinite(candidate.denominator) || candidate.denominator < 0.0) {
    invalidate_sample_diagnostics(MVMC_KRYLOV_STATUS_NONFINITE,
                                  diagnostics);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  candidate.zero_target_sample = candidate.denominator == 0.0;

  for (row = 0; row < dimension; ++row) {
    for (column = row; column < dimension; ++column) {
      size_t index = 0;
      double log_contribution = -INFINITY;
      status = mvmc_krylov_streaming_upper_index(
          dimension, row, column, &index);
      if (status != MVMC_KRYLOV_STATUS_OK) {
        invalidate_sample_diagnostics(status, diagnostics);
        return status;
      }
      status = guide_normalized_product(
          &dimensionless[row], &dimensionless[column], 0.0, log_guide,
          &overlap_upper[index], &log_contribution);
      if (status != MVMC_KRYLOV_STATUS_OK) {
        invalidate_sample_diagnostics(status, diagnostics);
        return status;
      }
      note_log_contribution(log_contribution, &candidate);

      status = guide_normalized_product(
          &dimensionless[row], &dimensionless[column + 1],
          policy->log_basis_scale[column] -
              policy->log_basis_scale[column + 1],
          log_guide, &hamiltonian_upper[index], &log_contribution);
      if (status != MVMC_KRYLOV_STATUS_OK) {
        invalidate_sample_diagnostics(status, diagnostics);
        return status;
      }
      note_log_contribution(log_contribution, &candidate);

      if (hamiltonian_adjoint_upper != NULL) {
        status = guide_normalized_product(
            &dimensionless[column], &dimensionless[row + 1],
            policy->log_basis_scale[row] -
                policy->log_basis_scale[row + 1],
            log_guide, &hamiltonian_adjoint_upper[index],
            &log_contribution);
        if (status != MVMC_KRYLOV_STATUS_OK) {
          invalidate_sample_diagnostics(status, diagnostics);
          return status;
        }
        note_log_contribution(log_contribution, &candidate);
      }

      status = guide_normalized_product(
          &dimensionless[row + 1], &dimensionless[column + 1],
          policy->log_basis_scale[row] - policy->log_basis_scale[row + 1] +
              policy->log_basis_scale[column] -
              policy->log_basis_scale[column + 1],
          log_guide, &hamiltonian_squared_upper[index],
          &log_contribution);
      if (status != MVMC_KRYLOV_STATUS_OK) {
        invalidate_sample_diagnostics(status, diagnostics);
        return status;
      }
      note_log_contribution(log_contribution, &candidate);
    }
  }
  if (candidate.zero_target_sample &&
      (!all_upper_entries_zero(overlap_upper, upper_count) ||
       !all_upper_entries_zero(hamiltonian_upper, upper_count) ||
       !all_upper_entries_zero(hamiltonian_adjoint_upper, upper_count) ||
       !all_upper_entries_zero(hamiltonian_squared_upper, upper_count))) {
    invalidate_sample_diagnostics(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                                  diagnostics);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }

  candidate.valid = 1;
  *diagnostics = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_matrix_measurement_accumulator_init(
    int order, MVMCKrylovStreamingComplexSum *overlap_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_squared_entries,
    size_t upper_count, MVMCKrylovMatrixMeasurementAccumulator *accumulator) {
  return mvmc_krylov_matrix_measurement_accumulator_init_with_diagnostics(
      order, overlap_entries, hamiltonian_entries, NULL,
      hamiltonian_squared_entries, upper_count, accumulator);
}

MVMCKrylovStatus mvmc_krylov_matrix_measurement_accumulator_init_with_diagnostics(
    int order, MVMCKrylovStreamingComplexSum *overlap_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_adjoint_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_squared_entries,
    size_t upper_count, MVMCKrylovMatrixMeasurementAccumulator *accumulator) {
  MVMCKrylovMatrixMeasurementAccumulator candidate;
  size_t dimension = 0;
  size_t required_upper_count = 0;
  MVMCKrylovStatus status;
  if (accumulator == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_accumulator(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, accumulator);
  status = mvmc_krylov_matrix_measurement_dimension(
      order, &dimension, &required_upper_count);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_accumulator(status, accumulator);
    return status;
  }
  if (upper_count != required_upper_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (overlap_entries == hamiltonian_entries ||
      overlap_entries == hamiltonian_squared_entries ||
      hamiltonian_entries == hamiltonian_squared_entries ||
      (hamiltonian_adjoint_entries != NULL &&
       (hamiltonian_adjoint_entries == overlap_entries ||
        hamiltonian_adjoint_entries == hamiltonian_entries ||
        hamiltonian_adjoint_entries == hamiltonian_squared_entries))) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  status = mvmc_krylov_streaming_matrix_accumulator_init(
      dimension, overlap_entries, upper_count, &candidate.overlap);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_streaming_matrix_accumulator_init(
        dimension, hamiltonian_entries, upper_count, &candidate.hamiltonian);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_streaming_matrix_accumulator_init(
        dimension, hamiltonian_squared_entries, upper_count,
        &candidate.hamiltonian_squared);
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      hamiltonian_adjoint_entries != NULL) {
    status = mvmc_krylov_streaming_matrix_accumulator_init(
        dimension, hamiltonian_adjoint_entries, upper_count,
        &candidate.hamiltonian_adjoint);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      candidate.has_hamiltonian_adjoint = 1;
    }
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_accumulator(status, accumulator);
    return status;
  }
  real_sum_reset(&candidate.denominator);
  real_sum_reset(&candidate.denominator_squared);
  candidate.minimum_denominator = NAN;
  candidate.maximum_denominator = NAN;
  candidate.minimum_log_contribution = NAN;
  candidate.maximum_log_contribution = NAN;
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.order = order;
  candidate.dimension = dimension;
  candidate.upper_count = upper_count;
  *accumulator = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_matrix_measurement_accumulator_add_sample(
    MVMCKrylovMatrixMeasurementAccumulator *accumulator,
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    const MVMCScaledComplex *values, size_t value_count,
    MVMCKrylovMatrixMeasurementSampleDiagnostics *diagnostics) {
  double complex overlap[MVMC_KRYLOV_MAX_ORDER *
                         (MVMC_KRYLOV_MAX_ORDER + 1) / 2];
  double complex hamiltonian[MVMC_KRYLOV_MAX_ORDER *
                             (MVMC_KRYLOV_MAX_ORDER + 1) / 2];
  double complex hamiltonian_adjoint[MVMC_KRYLOV_MAX_ORDER *
                                     (MVMC_KRYLOV_MAX_ORDER + 1) / 2];
  double complex hamiltonian_squared[MVMC_KRYLOV_MAX_ORDER *
                                     (MVMC_KRYLOV_MAX_ORDER + 1) / 2];
  MVMCKrylovMatrixMeasurementSampleDiagnostics local_diagnostics;
  MVMCKrylovStatus status;

  if (accumulator == NULL || !accumulator->valid ||
      accumulator->status != MVMC_KRYLOV_STATUS_OK ||
      accumulator->sample_count == UINT64_MAX ||
      accumulator->overlap.sample_count != accumulator->sample_count ||
      accumulator->hamiltonian.sample_count != accumulator->sample_count ||
      accumulator->hamiltonian_squared.sample_count !=
          accumulator->sample_count ||
      (accumulator->has_hamiltonian_adjoint &&
       accumulator->hamiltonian_adjoint.sample_count !=
           accumulator->sample_count)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_matrix_measurement_sample_with_adjoint(
      policy, values, value_count, overlap, hamiltonian,
      accumulator->has_hamiltonian_adjoint ? hamiltonian_adjoint : NULL,
      hamiltonian_squared, accumulator->upper_count, &local_diagnostics);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    accumulator->status = status;
    if (diagnostics != NULL) *diagnostics = local_diagnostics;
    return status;
  }
  status = mvmc_krylov_streaming_matrix_accumulator_add_sample(
      &accumulator->overlap, overlap, accumulator->upper_count, 1.0);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_streaming_matrix_accumulator_add_sample(
        &accumulator->hamiltonian, hamiltonian, accumulator->upper_count,
        1.0);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    if (accumulator->has_hamiltonian_adjoint) {
      status = mvmc_krylov_streaming_matrix_accumulator_add_sample(
          &accumulator->hamiltonian_adjoint, hamiltonian_adjoint,
          accumulator->upper_count, 1.0);
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_streaming_matrix_accumulator_add_sample(
        &accumulator->hamiltonian_squared, hamiltonian_squared,
        accumulator->upper_count, 1.0);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = real_sum_add(&accumulator->denominator,
                          local_diagnostics.denominator);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = real_sum_add(&accumulator->denominator_squared,
                          local_diagnostics.denominator *
                              local_diagnostics.denominator);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    accumulator->status = status;
    if (diagnostics != NULL) *diagnostics = local_diagnostics;
    return status;
  }
  if (local_diagnostics.zero_target_sample) {
    ++accumulator->zero_target_sample_count;
  }
  if (!isfinite(accumulator->minimum_denominator) ||
      local_diagnostics.denominator < accumulator->minimum_denominator) {
    accumulator->minimum_denominator = local_diagnostics.denominator;
  }
  if (!isfinite(accumulator->maximum_denominator) ||
      local_diagnostics.denominator > accumulator->maximum_denominator) {
    accumulator->maximum_denominator = local_diagnostics.denominator;
  }
  if (isfinite(local_diagnostics.minimum_log_contribution) &&
      (!isfinite(accumulator->minimum_log_contribution) ||
       local_diagnostics.minimum_log_contribution <
           accumulator->minimum_log_contribution)) {
    accumulator->minimum_log_contribution =
        local_diagnostics.minimum_log_contribution;
  }
  if (isfinite(local_diagnostics.maximum_log_contribution) &&
      (!isfinite(accumulator->maximum_log_contribution) ||
       local_diagnostics.maximum_log_contribution >
           accumulator->maximum_log_contribution)) {
    accumulator->maximum_log_contribution =
        local_diagnostics.maximum_log_contribution;
  }
  ++accumulator->sample_count;
  if (diagnostics != NULL) *diagnostics = local_diagnostics;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_matrix_measurement_accumulator_mean(
    const MVMCKrylovMatrixMeasurementAccumulator *accumulator,
    double complex *overlap_upper, double complex *hamiltonian_upper,
    double complex *hamiltonian_squared_upper, size_t upper_count,
    double *denominator_mean) {
  double sample_count;
  MVMCKrylovStatus status;
  if (accumulator == NULL || !accumulator->valid ||
      accumulator->status != MVMC_KRYLOV_STATUS_OK ||
      denominator_mean == NULL || accumulator->sample_count == 0 ||
      upper_count != accumulator->upper_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  sample_count = (double)accumulator->sample_count;
  if (!isfinite(sample_count) || sample_count <= 0.0) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  status = mvmc_krylov_streaming_matrix_accumulator_mean(
      &accumulator->overlap, overlap_upper, upper_count);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_streaming_matrix_accumulator_mean(
        &accumulator->hamiltonian, hamiltonian_upper, upper_count);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_streaming_matrix_accumulator_mean(
        &accumulator->hamiltonian_squared, hamiltonian_squared_upper,
        upper_count);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  *denominator_mean = real_sum_value(&accumulator->denominator) /
                      sample_count;
  return isfinite(*denominator_mean)
             ? MVMC_KRYLOV_STATUS_OK
             : MVMC_KRYLOV_STATUS_NONFINITE;
}

MVMCKrylovStatus mvmc_krylov_matrix_measurement_extract_block(
    const MVMCKrylovMatrixMeasurementAccumulator *accumulator,
    MVMCKrylovMatrixKind matrix_kind, size_t row, size_t column,
    MVMCKrylovJackknifeBlock *block) {
  const MVMCKrylovStreamingMatrixAccumulator *matrix = NULL;
  size_t index = 0;
  MVMCKrylovStatus status;
  if (block == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(block, 0, sizeof(*block));
  if (accumulator == NULL || !accumulator->valid ||
      accumulator->status != MVMC_KRYLOV_STATUS_OK ||
      row >= accumulator->dimension || column >= accumulator->dimension) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  switch (matrix_kind) {
    case MVMC_KRYLOV_MATRIX_OVERLAP:
      matrix = &accumulator->overlap;
      break;
    case MVMC_KRYLOV_MATRIX_HAMILTONIAN:
      matrix = &accumulator->hamiltonian;
      break;
    case MVMC_KRYLOV_MATRIX_HAMILTONIAN_ADJOINT:
      if (!accumulator->has_hamiltonian_adjoint) {
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      }
      matrix = &accumulator->hamiltonian_adjoint;
      break;
    case MVMC_KRYLOV_MATRIX_HAMILTONIAN_SQUARED:
      matrix = &accumulator->hamiltonian_squared;
      break;
    default:
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_streaming_upper_index(
      accumulator->dimension, row, column, &index);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  block->numerator =
      matrix->entries[index].real.sum + matrix->entries[index].real.compensation +
      I * (matrix->entries[index].imag.sum +
           matrix->entries[index].imag.compensation);
  block->denominator = real_sum_value(&accumulator->denominator);
  block->sample_count = accumulator->sample_count;
  return isfinite(block->denominator) && finite_complex(block->numerator)
             ? MVMC_KRYLOV_STATUS_OK
             : MVMC_KRYLOV_STATUS_NONFINITE;
}

static double complex streaming_entry_value(
    const MVMCKrylovStreamingMatrixAccumulator *matrix, size_t index) {
  return matrix->entries[index].real.sum +
         matrix->entries[index].real.compensation +
         I * (matrix->entries[index].imag.sum +
              matrix->entries[index].imag.compensation);
}

static MVMCKrylovStatus block_accumulator_capacity(
    const MVMCKrylovMatrixMeasurementBlockAccumulator *accumulator,
    uint64_t *capacity) {
  if (capacity == NULL || accumulator == NULL ||
      accumulator->block_count == 0 || accumulator->block_length == 0 ||
      accumulator->block_count >
          (size_t)(UINT64_MAX / accumulator->block_length)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  *capacity = (uint64_t)accumulator->block_count *
              accumulator->block_length;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_matrix_measurement_block_accumulator_init(
    int order, size_t block_count, uint64_t block_length,
    MVMCKrylovMatrixMeasurementAccumulator *blocks,
    MVMCKrylovStreamingComplexSum *overlap_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_adjoint_entries,
    MVMCKrylovStreamingComplexSum *hamiltonian_squared_entries,
    size_t storage_entry_count,
    MVMCKrylovMatrixMeasurementBlockAccumulator *accumulator) {
  MVMCKrylovMatrixMeasurementBlockAccumulator candidate;
  size_t dimension = 0;
  size_t upper_count = 0;
  size_t required_storage_count;
  size_t block;
  MVMCKrylovStatus status;

  if (accumulator == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_block_accumulator(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                               accumulator);
  status = mvmc_krylov_matrix_measurement_dimension(
      order, &dimension, &upper_count);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_block_accumulator(status, accumulator);
    return status;
  }
  if (block_count == 0 || block_length == 0 ||
      block_count > (size_t)(UINT64_MAX / block_length) ||
      upper_count > SIZE_MAX / block_count) {
    invalidate_block_accumulator(MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                                 accumulator);
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  required_storage_count = upper_count * block_count;
  if (blocks == NULL || overlap_entries == NULL ||
      hamiltonian_entries == NULL || hamiltonian_squared_entries == NULL ||
      storage_entry_count != required_storage_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }

  memset(&candidate, 0, sizeof(candidate));
  for (block = 0; block < block_count; ++block) {
    const size_t offset = block * upper_count;
    status =
        mvmc_krylov_matrix_measurement_accumulator_init_with_diagnostics(
            order, overlap_entries + offset, hamiltonian_entries + offset,
            hamiltonian_adjoint_entries == NULL
                ? NULL
                : hamiltonian_adjoint_entries + offset,
            hamiltonian_squared_entries + offset, upper_count,
            &blocks[block]);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      invalidate_block_accumulator(status, accumulator);
      return status;
    }
  }
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.order = order;
  candidate.dimension = dimension;
  candidate.upper_count = upper_count;
  candidate.block_count = block_count;
  candidate.block_length = block_length;
  candidate.sample_count = 0;
  candidate.has_hamiltonian_adjoint =
      hamiltonian_adjoint_entries != NULL ? 1 : 0;
  candidate.blocks = blocks;
  *accumulator = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_matrix_measurement_block_accumulator_add_sample(
    MVMCKrylovMatrixMeasurementBlockAccumulator *accumulator,
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    const MVMCScaledComplex *values, size_t value_count,
    MVMCKrylovMatrixMeasurementSampleDiagnostics *diagnostics) {
  uint64_t capacity = 0;
  size_t block_index;
  MVMCKrylovStatus status;

  if (accumulator == NULL || !accumulator->valid ||
      accumulator->status != MVMC_KRYLOV_STATUS_OK ||
      accumulator->blocks == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = block_accumulator_capacity(accumulator, &capacity);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    accumulator->status = status;
    return status;
  }
  if (accumulator->sample_count >= capacity) {
    accumulator->status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  block_index =
      (size_t)(accumulator->sample_count / accumulator->block_length);
  if (block_index >= accumulator->block_count) {
    accumulator->status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  status = mvmc_krylov_matrix_measurement_accumulator_add_sample(
      &accumulator->blocks[block_index], policy, values, value_count,
      diagnostics);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    accumulator->status = status;
    return status;
  }
  ++accumulator->sample_count;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_matrix_measurement_block_accumulator_extract_blocks(
    const MVMCKrylovMatrixMeasurementBlockAccumulator *accumulator,
    MVMCKrylovMatrixKind matrix_kind, size_t row, size_t column,
    MVMCKrylovJackknifeBlock *blocks, size_t block_count) {
  uint64_t capacity = 0;
  size_t block;
  MVMCKrylovStatus status;

  if (blocks == NULL || accumulator == NULL || !accumulator->valid ||
      accumulator->status != MVMC_KRYLOV_STATUS_OK ||
      accumulator->blocks == NULL || block_count != accumulator->block_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = block_accumulator_capacity(accumulator, &capacity);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (accumulator->sample_count != capacity) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (block = 0; block < block_count; ++block) {
    if (accumulator->blocks[block].sample_count !=
        accumulator->block_length) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    status = mvmc_krylov_matrix_measurement_extract_block(
        &accumulator->blocks[block], matrix_kind, row, column,
        &blocks[block]);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus hamiltonian_antihermitian_residual(
    const MVMCKrylovMatrixMeasurementAccumulator *accumulator,
    double *residual, double *norm) {
  double residual_squared = 0.0;
  double norm_squared = 0.0;
  size_t row;
  size_t column;
  double sample_count;
  if (residual == NULL || norm == NULL || accumulator == NULL ||
      !accumulator->has_hamiltonian_adjoint ||
      accumulator->sample_count == 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  sample_count = (double)accumulator->sample_count;
  if (!isfinite(sample_count) || sample_count <= 0.0) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  for (row = 0; row < accumulator->dimension; ++row) {
    for (column = row; column < accumulator->dimension; ++column) {
      size_t index = 0;
      const MVMCKrylovStatus status = mvmc_krylov_streaming_upper_index(
          accumulator->dimension, row, column, &index);
      double complex forward;
      double complex reverse;
      double complex delta;
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      forward =
          streaming_entry_value(&accumulator->hamiltonian, index) /
          sample_count;
      reverse =
          streaming_entry_value(&accumulator->hamiltonian_adjoint, index) /
          sample_count;
      delta = forward - conj(reverse);
      residual_squared +=
          (row == column ? 1.0 : 2.0) *
          (creal(delta) * creal(delta) + cimag(delta) * cimag(delta));
      norm_squared += creal(forward) * creal(forward) +
                      cimag(forward) * cimag(forward);
      if (row != column) {
        norm_squared += creal(reverse) * creal(reverse) +
                        cimag(reverse) * cimag(reverse);
      }
      if (!isfinite(residual_squared) || !isfinite(norm_squared)) {
        return MVMC_KRYLOV_STATUS_NONFINITE;
      }
    }
  }
  *norm = sqrt(norm_squared);
  *residual = *norm > 0.0 ? sqrt(residual_squared) / *norm
                          : sqrt(residual_squared);
  return isfinite(*norm) && isfinite(*residual)
             ? MVMC_KRYLOV_STATUS_OK
             : MVMC_KRYLOV_STATUS_NONFINITE;
}

MVMCKrylovStatus mvmc_krylov_matrix_measurement_diagnostics_summary(
    const MVMCKrylovMatrixMeasurementAccumulator *accumulator,
    double minimum_abs_denominator_mean,
    double maximum_denominator_relative_se,
    MVMCKrylovMatrixMeasurementDiagnosticsSummary *summary) {
  MVMCKrylovMatrixMeasurementDiagnosticsSummary candidate;
  double sample_count;
  double denominator_sum;
  double denominator_squared_sum;
  double denominator_mean;
  double denominator_ss;
  double denominator_variance = 0.0;
  double denominator_relative_se = 0.0;
  MVMCKrylovStatus status;

  if (summary == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_diagnostics_summary(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                                 summary);
  if (accumulator == NULL || !accumulator->valid ||
      accumulator->status != MVMC_KRYLOV_STATUS_OK ||
      accumulator->sample_count == 0 ||
      accumulator->zero_target_sample_count > accumulator->sample_count ||
      !isfinite(minimum_abs_denominator_mean) ||
      minimum_abs_denominator_mean < 0.0 ||
      !isfinite(maximum_denominator_relative_se) ||
      maximum_denominator_relative_se < 0.0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  sample_count = (double)accumulator->sample_count;
  if (!isfinite(sample_count) || sample_count <= 0.0) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  denominator_sum = real_sum_value(&accumulator->denominator);
  denominator_squared_sum = real_sum_value(&accumulator->denominator_squared);
  if (!isfinite(denominator_sum) || denominator_sum <= 0.0 ||
      !isfinite(denominator_squared_sum) || denominator_squared_sum <= 0.0) {
    invalidate_diagnostics_summary(MVMC_KRYLOV_STATUS_NONFINITE, summary);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  denominator_mean = denominator_sum / sample_count;
  denominator_ss =
      denominator_squared_sum - 2.0 * denominator_mean * denominator_sum +
      sample_count * denominator_mean * denominator_mean;
  if (denominator_ss < 0.0 &&
      denominator_ss >= -1.0e-10 * fmax(1.0, denominator_squared_sum)) {
    denominator_ss = 0.0;
  }
  if (!isfinite(denominator_mean) || denominator_mean <= 0.0 ||
      !isfinite(denominator_ss) || denominator_ss < 0.0) {
    invalidate_diagnostics_summary(MVMC_KRYLOV_STATUS_NONFINITE, summary);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  if (accumulator->sample_count > 1) {
    denominator_variance = denominator_ss / (sample_count - 1.0);
    denominator_relative_se =
        sqrt(denominator_ss / ((sample_count - 1.0) * sample_count)) /
        fabs(denominator_mean);
  }

  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.order = accumulator->order;
  candidate.dimension = accumulator->dimension;
  candidate.upper_count = accumulator->upper_count;
  candidate.sample_count = accumulator->sample_count;
  candidate.zero_target_sample_count =
      accumulator->zero_target_sample_count;
  candidate.denominator_sum = denominator_sum;
  candidate.denominator_mean = denominator_mean;
  candidate.denominator_sample_variance = denominator_variance;
  candidate.denominator_relative_se = denominator_relative_se;
  candidate.minimum_abs_denominator_mean = minimum_abs_denominator_mean;
  candidate.maximum_denominator_relative_se =
      maximum_denominator_relative_se;
  candidate.denominator_stable =
      accumulator->sample_count > 1 &&
      fabs(denominator_mean) >= minimum_abs_denominator_mean &&
      denominator_relative_se <= maximum_denominator_relative_se;
  candidate.effective_sample_count =
      denominator_sum * denominator_sum / denominator_squared_sum;
  candidate.effective_sample_fraction =
      candidate.effective_sample_count / sample_count;
  candidate.zero_target_sample_fraction =
      (double)accumulator->zero_target_sample_count / sample_count;
  candidate.minimum_denominator = accumulator->minimum_denominator;
  candidate.maximum_denominator = accumulator->maximum_denominator;
  candidate.denominator_tail_ratio =
      accumulator->maximum_denominator / denominator_mean;
  candidate.minimum_log_contribution =
      accumulator->minimum_log_contribution;
  candidate.maximum_log_contribution =
      accumulator->maximum_log_contribution;
  candidate.log_contribution_span =
      isfinite(candidate.minimum_log_contribution) &&
              isfinite(candidate.maximum_log_contribution)
          ? candidate.maximum_log_contribution -
                candidate.minimum_log_contribution
          : NAN;
  if (!isfinite(candidate.effective_sample_count) ||
      !isfinite(candidate.effective_sample_fraction) ||
      !isfinite(candidate.zero_target_sample_fraction) ||
      !isfinite(candidate.denominator_tail_ratio)) {
    invalidate_diagnostics_summary(MVMC_KRYLOV_STATUS_NONFINITE, summary);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  if (accumulator->has_hamiltonian_adjoint) {
    status = hamiltonian_antihermitian_residual(
        accumulator, &candidate.hamiltonian_antihermitian_residual,
        &candidate.hamiltonian_norm);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      invalidate_diagnostics_summary(status, summary);
      return status;
    }
    candidate.hamiltonian_residual_available = 1;
  } else {
    candidate.hamiltonian_antihermitian_residual = NAN;
    candidate.hamiltonian_norm = NAN;
    candidate.hamiltonian_residual_available = 0;
  }
  *summary = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}
