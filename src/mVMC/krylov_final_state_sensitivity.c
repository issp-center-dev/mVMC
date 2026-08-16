/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "krylov_final_state_sensitivity.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_P5_TESTING) ||                          \
    !defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                         \
    !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "krylov_final_state_sensitivity.c is Testing-only"
#endif

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

static int finite_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static uint64_t hash_u64(uint64_t hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    hash ^= value & UINT64_C(0xff);
    hash *= UINT64_C(1099511628211);
    value >>= 8;
  }
  return hash;
}

static uint64_t hash_double(uint64_t hash, double value) {
  uint64_t bits;
  memcpy(&bits, &value, sizeof(bits));
  return hash_u64(hash, bits);
}

static size_t required_upper_count(int dimension) {
  return (size_t)dimension * ((size_t)dimension + 1) / 2;
}

static uint64_t compute_model_hash(
    const MVMCKrylovFinalStateQuadraticModel *model) {
  const double complex *matrices[] = {
      model->overlap_upper, model->hamiltonian_upper,
      model->hamiltonian_squared_upper, model->observable_upper};
  uint64_t hash = UINT64_C(1469598103934665603);
  size_t matrix;
  size_t entry;
  hash = hash_u64(hash, MVMC_KRYLOV_FINAL_SENSITIVITY_MODEL_VERSION);
  hash = hash_u64(hash, (uint64_t)(unsigned int)model->dimension);
  hash = hash_u64(hash, (uint64_t)(unsigned int)model->has_observable);
  for (matrix = 0; matrix < 4; ++matrix) {
    for (entry = 0; entry < model->upper_count; ++entry) {
      hash = hash_double(hash, creal(matrices[matrix][entry]));
      hash = hash_double(hash, cimag(matrices[matrix][entry]));
    }
  }
  return hash == 0 ? UINT64_C(1) : hash;
}

static int model_valid(const MVMCKrylovFinalStateQuadraticModel *model) {
  size_t entry;
  int row;
  int column;
  if (model == NULL || !model->valid ||
      model->status != MVMC_KRYLOV_STATUS_OK || model->dimension < 2 ||
      model->dimension > MVMC_KRYLOV_MAX_ORDER ||
      model->upper_count != required_upper_count(model->dimension) ||
      (model->has_observable != 0 && model->has_observable != 1) ||
      model->model_hash == 0) {
    return 0;
  }
  entry = 0;
  for (row = 0; row < model->dimension; ++row) {
    for (column = row; column < model->dimension; ++column, ++entry) {
      if (!finite_complex(model->overlap_upper[entry]) ||
          !finite_complex(model->hamiltonian_upper[entry]) ||
          !finite_complex(model->hamiltonian_squared_upper[entry]) ||
          !finite_complex(model->observable_upper[entry])) {
        return 0;
      }
      if (row == column &&
          (cimag(model->overlap_upper[entry]) != 0.0 ||
           cimag(model->hamiltonian_upper[entry]) != 0.0 ||
           cimag(model->hamiltonian_squared_upper[entry]) != 0.0 ||
           cimag(model->observable_upper[entry]) != 0.0)) {
        return 0;
      }
      if (!model->has_observable &&
          (creal(model->observable_upper[entry]) != 0.0 ||
           cimag(model->observable_upper[entry]) != 0.0)) {
        return 0;
      }
    }
  }
  return model->model_hash == compute_model_hash(model);
}

static void invalidate_model(MVMCKrylovStatus status,
                             MVMCKrylovFinalStateQuadraticModel *model) {
  if (model == NULL) return;
  memset(model, 0, sizeof(*model));
  model->status = status;
}

static void invalidate_metrics(MVMCKrylovStatus status,
                               MVMCKrylovFinalStateQuadraticMetrics *metrics) {
  if (metrics == NULL) return;
  memset(metrics, 0, sizeof(*metrics));
  metrics->status = status;
  metrics->norm = NAN;
  metrics->energy = NAN;
  metrics->hamiltonian_second_moment = NAN;
  metrics->full_support_variance = NAN;
  metrics->observable = NAN;
}

static void invalidate_comparison(
    MVMCKrylovStatus status,
    MVMCKrylovFinalStateSensitivityComparison *comparison) {
  if (comparison == NULL) return;
  memset(comparison, 0, sizeof(*comparison));
  comparison->status = status;
  comparison->coefficient_projective_distance = NAN;
  comparison->overlap_relative_frobenius_error = NAN;
  comparison->hamiltonian_relative_frobenius_error = NAN;
  comparison->hamiltonian_squared_relative_frobenius_error = NAN;
  comparison->observable_relative_frobenius_error = NAN;
}

static double complex quadratic_form(
    const double complex *upper, int dimension,
    const double complex *coefficient, double *absolute_sum) {
  double complex sum = 0.0;
  double magnitude_sum = 0.0;
  size_t entry = 0;
  int row;
  int column;
  for (row = 0; row < dimension; ++row) {
    for (column = row; column < dimension; ++column, ++entry) {
      const double complex forward =
          conj(coefficient[row]) * upper[entry] * coefficient[column];
      if (row == column) {
        sum += forward;
        magnitude_sum += cabs(forward);
      } else {
        const double complex reverse = conj(forward);
        sum += forward + reverse;
        magnitude_sum += cabs(forward) + cabs(reverse);
      }
    }
  }
  if (absolute_sum != NULL) *absolute_sum = magnitude_sum;
  return sum;
}

static int real_quadratic_value(double complex value, double absolute_sum,
                                double *real_value) {
  const double tolerance =
      256.0 * DBL_EPSILON * fmax(1.0, absolute_sum);
  if (real_value == NULL || !finite_complex(value) ||
      !isfinite(absolute_sum) || absolute_sum < 0.0 ||
      fabs(cimag(value)) > tolerance) {
    return 0;
  }
  *real_value = creal(value);
  return isfinite(*real_value);
}

MVMCKrylovStatus mvmc_krylov_final_state_quadratic_model_create(
    int dimension, const double complex *overlap_upper,
    const double complex *hamiltonian_upper,
    const double complex *hamiltonian_squared_upper,
    const double complex *observable_upper, size_t upper_count,
    MVMCKrylovFinalStateQuadraticModel *model) {
  MVMCKrylovFinalStateQuadraticModel candidate;
  size_t entry = 0;
  int row;
  int column;
  if (model == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_model(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, model);
  if (dimension < 2 || dimension > MVMC_KRYLOV_MAX_ORDER ||
      upper_count != required_upper_count(dimension) ||
      overlap_upper == NULL || hamiltonian_upper == NULL ||
      hamiltonian_squared_upper == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.dimension = dimension;
  candidate.upper_count = upper_count;
  candidate.has_observable = observable_upper != NULL;
  for (row = 0; row < dimension; ++row) {
    for (column = row; column < dimension; ++column, ++entry) {
      const double complex observable =
          observable_upper == NULL ? 0.0 : observable_upper[entry];
      if (!finite_complex(overlap_upper[entry]) ||
          !finite_complex(hamiltonian_upper[entry]) ||
          !finite_complex(hamiltonian_squared_upper[entry]) ||
          !finite_complex(observable) ||
          (row == column &&
           (cimag(overlap_upper[entry]) != 0.0 ||
            cimag(hamiltonian_upper[entry]) != 0.0 ||
            cimag(hamiltonian_squared_upper[entry]) != 0.0 ||
            cimag(observable) != 0.0))) {
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      }
      candidate.overlap_upper[entry] = overlap_upper[entry];
      candidate.hamiltonian_upper[entry] = hamiltonian_upper[entry];
      candidate.hamiltonian_squared_upper[entry] =
          hamiltonian_squared_upper[entry];
      candidate.observable_upper[entry] = observable;
    }
  }
  candidate.model_hash = compute_model_hash(&candidate);
  candidate.valid = 1;
  *model = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_final_state_quadratic_metrics(
    const MVMCKrylovFinalStateQuadraticModel *model,
    const double complex *coefficient, size_t coefficient_count,
    MVMCKrylovFinalStateQuadraticMetrics *metrics) {
  MVMCKrylovFinalStateQuadraticMetrics candidate;
  double norm_abs;
  double energy_abs;
  double second_abs;
  double observable_abs = 0.0;
  double norm;
  double energy_numerator;
  double second_numerator;
  double observable_numerator = 0.0;
  double complex norm_value;
  double complex energy_value;
  double complex second_value;
  double complex observable_value = 0.0;
  double negative_tolerance;
  size_t index;
  if (metrics == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_metrics(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, metrics);
  if (!model_valid(model) || coefficient == NULL ||
      coefficient_count != (size_t)model->dimension) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (index = 0; index < coefficient_count; ++index) {
    if (!finite_complex(coefficient[index])) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  norm_value = quadratic_form(model->overlap_upper, model->dimension,
                              coefficient, &norm_abs);
  energy_value = quadratic_form(model->hamiltonian_upper, model->dimension,
                                coefficient, &energy_abs);
  second_value = quadratic_form(model->hamiltonian_squared_upper,
                                model->dimension, coefficient, &second_abs);
  if (model->has_observable) {
    observable_value = quadratic_form(model->observable_upper,
                                      model->dimension, coefficient,
                                      &observable_abs);
  }
  if (!real_quadratic_value(norm_value, norm_abs, &norm) ||
      !real_quadratic_value(energy_value, energy_abs, &energy_numerator) ||
      !real_quadratic_value(second_value, second_abs, &second_numerator) ||
      (model->has_observable &&
       !real_quadratic_value(observable_value, observable_abs,
                             &observable_numerator)) ||
      !isfinite(norm) || norm <= 0.0) {
    invalidate_metrics(MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE,
                       metrics);
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.dimension = model->dimension;
  candidate.has_observable = model->has_observable;
  candidate.norm = norm;
  candidate.energy = energy_numerator / norm;
  candidate.hamiltonian_second_moment = second_numerator / norm;
  candidate.full_support_variance =
      candidate.hamiltonian_second_moment -
      candidate.energy * candidate.energy;
  candidate.observable = model->has_observable
                             ? observable_numerator / norm
                             : NAN;
  negative_tolerance = 512.0 * DBL_EPSILON *
                       fmax(1.0, fabs(candidate.hamiltonian_second_moment) +
                                     candidate.energy * candidate.energy);
  if (!isfinite(candidate.energy) ||
      !isfinite(candidate.hamiltonian_second_moment) ||
      (model->has_observable && !isfinite(candidate.observable)) ||
      !isfinite(candidate.full_support_variance) ||
      candidate.full_support_variance < -negative_tolerance) {
    invalidate_metrics(MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE,
                       metrics);
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (candidate.full_support_variance < 0.0) {
    candidate.full_support_variance = 0.0;
  }
  candidate.valid = 1;
  *metrics = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

static double relative_frobenius_error(
    const double complex *baseline, const double complex *candidate,
    int dimension) {
  double baseline_squared = 0.0;
  double difference_squared = 0.0;
  size_t entry = 0;
  int row;
  int column;
  for (row = 0; row < dimension; ++row) {
    for (column = row; column < dimension; ++column, ++entry) {
      const double weight = row == column ? 1.0 : 2.0;
      baseline_squared += weight * pow(cabs(baseline[entry]), 2.0);
      difference_squared +=
          weight * pow(cabs(candidate[entry] - baseline[entry]), 2.0);
    }
  }
  if (!isfinite(baseline_squared) || !isfinite(difference_squared)) {
    return NAN;
  }
  return sqrt(difference_squared) /
         fmax(sqrt(baseline_squared), DBL_MIN);
}

static double coefficient_projective_distance(
    const double complex *baseline, const double complex *candidate,
    size_t count) {
  double baseline_norm = 0.0;
  double candidate_norm = 0.0;
  double complex inner = 0.0;
  double cosine;
  size_t index;
  for (index = 0; index < count; ++index) {
    baseline_norm += pow(cabs(baseline[index]), 2.0);
    candidate_norm += pow(cabs(candidate[index]), 2.0);
    inner += conj(baseline[index]) * candidate[index];
  }
  if (!isfinite(baseline_norm) || !isfinite(candidate_norm) ||
      !finite_complex(inner) || baseline_norm <= 0.0 ||
      candidate_norm <= 0.0) {
    return NAN;
  }
  cosine = cabs(inner) / sqrt(baseline_norm * candidate_norm);
  cosine = fmax(0.0, fmin(1.0, cosine));
  return sqrt(fmax(0.0, 1.0 - cosine * cosine));
}

static MVMCKrylovFinalStateDownstreamError downstream_error(
    const MVMCKrylovFinalStateQuadraticMetrics *baseline,
    const MVMCKrylovFinalStateQuadraticMetrics *candidate) {
  MVMCKrylovFinalStateDownstreamError error;
  memset(&error, 0, sizeof(error));
  error.energy_absolute = fabs(candidate->energy - baseline->energy);
  error.energy_scaled = error.energy_absolute /
                        fmax(1.0, fabs(baseline->energy));
  error.variance_absolute =
      fabs(candidate->full_support_variance -
           baseline->full_support_variance);
  error.variance_scaled =
      error.variance_absolute /
      fmax(1.0, fabs(baseline->full_support_variance));
  if (baseline->has_observable) {
    error.observable_absolute =
        fabs(candidate->observable - baseline->observable);
    error.observable_scaled =
        error.observable_absolute /
        fmax(1.0, fabs(baseline->observable));
  } else {
    error.observable_absolute = NAN;
    error.observable_scaled = NAN;
  }
  return error;
}

MVMCKrylovStatus mvmc_krylov_final_state_sensitivity_compare(
    const MVMCKrylovFinalStateQuadraticModel *baseline_model,
    const double complex *baseline_coefficient,
    const MVMCKrylovFinalStateQuadraticModel *candidate_model,
    const double complex *candidate_coefficient,
    size_t coefficient_count,
    MVMCKrylovFinalStateSensitivityComparison *comparison) {
  MVMCKrylovFinalStateSensitivityComparison candidate;
  MVMCKrylovStatus status;
  if (comparison == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_comparison(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, comparison);
  if (!model_valid(baseline_model) || !model_valid(candidate_model) ||
      baseline_model->dimension != candidate_model->dimension ||
      baseline_model->has_observable != candidate_model->has_observable ||
      baseline_coefficient == NULL || candidate_coefficient == NULL ||
      coefficient_count != (size_t)baseline_model->dimension) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.dimension = baseline_model->dimension;
  candidate.has_observable = baseline_model->has_observable;
  status = mvmc_krylov_final_state_quadratic_metrics(
      baseline_model, baseline_coefficient, coefficient_count,
      &candidate.baseline);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_final_state_quadratic_metrics(
        baseline_model, candidate_coefficient, coefficient_count,
        &candidate.coefficient_only_metrics);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_final_state_quadratic_metrics(
        candidate_model, baseline_coefficient, coefficient_count,
        &candidate.matrix_only_metrics);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_final_state_quadratic_metrics(
        candidate_model, candidate_coefficient, coefficient_count,
        &candidate.combined_metrics);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_comparison(status, comparison);
    return status;
  }
  candidate.coefficient_projective_distance =
      coefficient_projective_distance(baseline_coefficient,
                                      candidate_coefficient,
                                      coefficient_count);
  candidate.overlap_relative_frobenius_error = relative_frobenius_error(
      baseline_model->overlap_upper, candidate_model->overlap_upper,
      baseline_model->dimension);
  candidate.hamiltonian_relative_frobenius_error =
      relative_frobenius_error(
          baseline_model->hamiltonian_upper,
          candidate_model->hamiltonian_upper, baseline_model->dimension);
  candidate.hamiltonian_squared_relative_frobenius_error =
      relative_frobenius_error(
          baseline_model->hamiltonian_squared_upper,
          candidate_model->hamiltonian_squared_upper,
          baseline_model->dimension);
  candidate.observable_relative_frobenius_error =
      baseline_model->has_observable
          ? relative_frobenius_error(
                baseline_model->observable_upper,
                candidate_model->observable_upper,
                baseline_model->dimension)
          : NAN;
  if (!isfinite(candidate.coefficient_projective_distance) ||
      !isfinite(candidate.overlap_relative_frobenius_error) ||
      !isfinite(candidate.hamiltonian_relative_frobenius_error) ||
      !isfinite(candidate.hamiltonian_squared_relative_frobenius_error) ||
      (candidate.has_observable &&
       !isfinite(candidate.observable_relative_frobenius_error))) {
    invalidate_comparison(MVMC_KRYLOV_STATUS_NONFINITE, comparison);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  candidate.coefficient_only = downstream_error(
      &candidate.baseline, &candidate.coefficient_only_metrics);
  candidate.matrix_only =
      downstream_error(&candidate.baseline, &candidate.matrix_only_metrics);
  candidate.combined =
      downstream_error(&candidate.baseline, &candidate.combined_metrics);
  candidate.valid = 1;
  *comparison = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}
