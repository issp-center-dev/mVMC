/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "krylov_final_state.h"

#if (!defined(MVMC_ENABLE_POWER_LANCZOS_P5_CORE) &&                            \
     !defined(MVMC_ENABLE_POWER_LANCZOS_P5_TESTING)) ||                        \
    !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "krylov_final_state.c requires the P5 final-state core"
#endif

#include <complex.h>
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

static uint64_t configuration_hash(const uint64_t *words,
                                   size_t word_count) {
  uint64_t hash = UINT64_C(1469598103934665603);
  size_t index;
  hash = hash_u64(hash, (uint64_t)word_count);
  for (index = 0; index < word_count; ++index) {
    hash = hash_u64(hash, words[index]);
  }
  return hash == 0 ? UINT64_C(1) : hash;
}

static uint64_t hash_scaled_complex(uint64_t hash,
                                    const MVMCScaledComplex *value) {
  hash = hash_u64(hash, (uint64_t)(unsigned int)value->state);
  hash = hash_double(hash, creal(value->phase));
  hash = hash_double(hash, cimag(value->phase));
  hash = hash_double(hash, value->log_abs);
  hash = hash_double(hash, value->log_abs_error_bound);
  hash = hash_double(hash, value->max_input_log_abs);
  hash = hash_double(hash, value->cancellation_log_abs);
  hash = hash_double(hash, value->cancellation_ratio);
  return hash;
}

uint64_t mvmc_krylov_final_state_snapshot_integrity_hash(
    const MVMCKrylovFinalStateSnapshot *snapshot) {
  uint64_t hash = UINT64_C(1469598103934665603);
  int order;
  if (snapshot == NULL || !snapshot->valid ||
      snapshot->status != MVMC_KRYLOV_STATUS_OK ||
      snapshot->krylov.evaluated_order < 0 ||
      snapshot->krylov.evaluated_order > MVMC_KRYLOV_MAX_ORDER ||
      !snapshot->final_state.valid ||
      snapshot->final_state.status != MVMC_KRYLOV_STATUS_OK) {
    return 0;
  }
  hash = hash_u64(hash, snapshot->policy_hash);
  hash = hash_u64(hash, snapshot->configuration_hash);
  hash = hash_u64(hash, snapshot->amplitude_generation_hash);
  hash = hash_u64(
      hash, (uint64_t)(unsigned int)snapshot->krylov.evaluated_order);
  for (order = 0; order <= snapshot->krylov.evaluated_order; ++order) {
    hash = hash_scaled_complex(hash, &snapshot->krylov.value[order]);
  }
  hash = hash_u64(
      hash, (uint64_t)(unsigned int)snapshot->final_state.order);
  hash = hash_u64(hash, snapshot->final_state.policy_hash);
  hash = hash_u64(
      hash, (uint64_t)(unsigned int)
                snapshot->final_state.coefficient_zero_count);
  hash = hash_u64(
      hash, (uint64_t)(unsigned int)
                snapshot->final_state.basis_finite_count);
  hash = hash_u64(
      hash, (uint64_t)(unsigned int)
                snapshot->final_state.basis_exact_zero_count);
  hash = hash_u64(
      hash, (uint64_t)(unsigned int)
                snapshot->final_state.basis_numeric_zero_count);
  hash = hash_u64(
      hash, (uint64_t)(unsigned int)snapshot->final_state.sampleable);
  hash = hash_u64(
      hash, (uint64_t)(unsigned int)
                snapshot->final_state.local_energy_available);
  hash = hash_double(hash, snapshot->final_state.log_weight);
  hash = hash_scaled_complex(hash, &snapshot->final_state.amplitude);
  hash = hash_scaled_complex(
      hash, &snapshot->final_state.hamiltonian_amplitude);
  hash = hash_scaled_complex(hash, &snapshot->final_state.local_energy);
  return hash == 0 ? UINT64_C(1) : hash;
}

static int word_ranges_overlap(const uint64_t *left, size_t left_count,
                               const uint64_t *right,
                               size_t right_count) {
  size_t left_bytes;
  size_t right_bytes;
  uintptr_t left_begin;
  uintptr_t right_begin;
  uintptr_t left_end;
  uintptr_t right_end;
  if (left == NULL || right == NULL ||
      left_count > SIZE_MAX / sizeof(left[0]) ||
      right_count > SIZE_MAX / sizeof(right[0])) {
    return 1;
  }
  left_bytes = left_count * sizeof(left[0]);
  right_bytes = right_count * sizeof(right[0]);
  left_begin = (uintptr_t)left;
  right_begin = (uintptr_t)right;
  if (left_begin > UINTPTR_MAX - left_bytes ||
      right_begin > UINTPTR_MAX - right_bytes) return 1;
  left_end = left_begin + left_bytes;
  right_end = right_begin + right_bytes;
  return left_begin < right_end && right_begin < left_end;
}

static int source_valid(MVMCKrylovFinalCoefficientSource source) {
  return source == MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT ||
         source == MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC ||
         source == MVMC_KRYLOV_FINAL_COEFFICIENT_PERTURBED_TESTING ||
         source == MVMC_KRYLOV_FINAL_COEFFICIENT_PRODUCTION_NOISY;
}

static int order_valid(int order) {
  return order > 0 && order < MVMC_KRYLOV_MAX_ORDER;
}

static uint64_t compute_policy_hash(
    const MVMCKrylovFinalStatePolicy *policy) {
  uint64_t hash = UINT64_C(1469598103934665603);
  int index;
  hash = hash_u64(hash, MVMC_KRYLOV_FINAL_STATE_POLICY_VERSION);
  hash = hash_u64(hash, (uint64_t)(unsigned int)policy->order);
  hash = hash_u64(hash, (uint64_t)(unsigned int)policy->source);
  hash = hash_u64(hash, policy->provenance_hash);
  for (index = 0; index <= MVMC_KRYLOV_MAX_ORDER; ++index) {
    hash = hash_double(hash, creal(policy->coefficient[index]));
    hash = hash_double(hash, cimag(policy->coefficient[index]));
  }
  return hash == 0 ? UINT64_C(1) : hash;
}

static int policy_valid(const MVMCKrylovFinalStatePolicy *policy) {
  int index;
  int nonzero = 0;
  if (policy == NULL || !order_valid(policy->order) ||
      !source_valid(policy->source) || policy->provenance_hash == 0 ||
      policy->policy_hash == 0) {
    return 0;
  }
  for (index = 0; index <= MVMC_KRYLOV_MAX_ORDER; ++index) {
    const double complex coefficient = policy->coefficient[index];
    if (!finite_complex(coefficient)) return 0;
    if (index <= policy->order) {
      nonzero = nonzero || creal(coefficient) != 0.0 ||
                    cimag(coefficient) != 0.0;
    } else if (creal(coefficient) != 0.0 || cimag(coefficient) != 0.0) {
      return 0;
    }
  }
  return nonzero && policy->policy_hash == compute_policy_hash(policy);
}

static void invalidate_evaluation(
    MVMCKrylovStatus status, MVMCKrylovFinalStateEvaluation *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
  result->log_weight = NAN;
}

static void invalidate_acceptance(
    MVMCKrylovStatus status, MVMCKrylovFinalStateAcceptance *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
  result->uniform = NAN;
  result->log_target_ratio = NAN;
  result->log_proposal_ratio = NAN;
  result->log_acceptance_ratio = NAN;
}

static void invalidate_snapshot(
    MVMCKrylovStatus status, MVMCKrylovFinalStateSnapshot *snapshot) {
  if (snapshot == NULL) return;
  memset(snapshot, 0, sizeof(*snapshot));
  snapshot->status = status;
}

static void invalidate_step(
    MVMCKrylovStatus status, MVMCKrylovFinalStateStepResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
}

static void invalidate_proposal_step(
    MVMCKrylovStatus status,
    MVMCKrylovFinalStateProposalStepResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
  result->log_proposal_ratio = NAN;
}

static int evaluation_valid(
    const MVMCKrylovFinalStateEvaluation *evaluation) {
  int basis_count;
  if (evaluation == NULL || !evaluation->valid ||
      evaluation->status != MVMC_KRYLOV_STATUS_OK ||
      !order_valid(evaluation->order) || evaluation->policy_hash == 0 ||
      evaluation->coefficient_zero_count < 0 ||
      evaluation->coefficient_zero_count > evaluation->order + 1 ||
      evaluation->basis_finite_count < 0 ||
      evaluation->basis_exact_zero_count < 0 ||
      evaluation->basis_numeric_zero_count < 0 ||
      !mvmc_scaled_complex_is_valid(&evaluation->amplitude) ||
      !mvmc_scaled_complex_is_valid(
          &evaluation->hamiltonian_amplitude)) {
    return 0;
  }
  basis_count = evaluation->basis_finite_count +
                evaluation->basis_exact_zero_count +
                evaluation->basis_numeric_zero_count;
  if (basis_count != evaluation->order + 2) return 0;
  if (evaluation->amplitude.state ==
      MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
    if (!isfinite(evaluation->log_weight) ||
        !mvmc_scaled_complex_is_valid(&evaluation->local_energy) ||
        evaluation->local_energy_available !=
            (evaluation->local_energy.state ==
                 MVMC_SCALED_COMPLEX_FINITE_NONZERO ||
             evaluation->local_energy.state ==
                 MVMC_SCALED_COMPLEX_EXACT_ZERO) ||
        evaluation->sampleable != evaluation->local_energy_available) {
      return 0;
    }
  } else {
    if (evaluation->sampleable || evaluation->local_energy_available) {
      return 0;
    }
    if (evaluation->amplitude.state == MVMC_SCALED_COMPLEX_EXACT_ZERO) {
      if (!(isinf(evaluation->log_weight) &&
            evaluation->log_weight < 0.0)) {
        return 0;
      }
    } else if (!isnan(evaluation->log_weight)) {
      return 0;
    }
  }
  return evaluation->amplitude.state != MVMC_SCALED_COMPLEX_NONFINITE &&
         evaluation->hamiltonian_amplitude.state !=
             MVMC_SCALED_COMPLEX_NONFINITE;
}

static MVMCKrylovStatus scaled_coefficient(
    double complex coefficient, MVMCScaledComplex *result) {
  MVMCPfaffianStatus status;
  if (creal(coefficient) == 0.0 && cimag(coefficient) == 0.0) {
    status = mvmc_scaled_complex_make_exact_zero(result);
  } else {
    status = mvmc_scaled_complex_from_raw(coefficient, result);
  }
  return status == MVMC_PFAFFIAN_STATUS_OK
             ? MVMC_KRYLOV_STATUS_OK
             : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

MVMCKrylovStatus mvmc_krylov_final_state_policy_create(
    int order, MVMCKrylovFinalCoefficientSource source,
    uint64_t provenance_hash, const double complex *coefficient,
    size_t coefficient_count, MVMCKrylovFinalStatePolicy *policy) {
  MVMCKrylovFinalStatePolicy candidate;
  size_t index;
  int nonzero = 0;
  if (policy == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(policy, 0, sizeof(*policy));
  if (!order_valid(order) || !source_valid(source) ||
      provenance_hash == 0 || coefficient == NULL ||
      coefficient_count != (size_t)order + 1) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.order = order;
  candidate.source = source;
  candidate.provenance_hash = provenance_hash;
  for (index = 0; index < coefficient_count; ++index) {
    double real_part;
    double imaginary_part;
    if (!finite_complex(coefficient[index])) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    real_part = creal(coefficient[index]);
    imaginary_part = cimag(coefficient[index]);
    if (real_part == 0.0) real_part = 0.0;
    if (imaginary_part == 0.0) imaginary_part = 0.0;
    candidate.coefficient[index] = real_part + I * imaginary_part;
    nonzero = nonzero || real_part != 0.0 || imaginary_part != 0.0;
  }
  if (!nonzero) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  candidate.policy_hash = compute_policy_hash(&candidate);
  *policy = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_final_state_policy_create_scaled_basis(
    int order, MVMCKrylovFinalCoefficientSource source,
    uint64_t provenance_hash, const double complex *scaled_coefficient,
    const double *log_basis_scale, size_t coefficient_count,
    MVMCKrylovFinalStatePolicy *policy) {
  double log_raw_magnitude[MVMC_KRYLOV_MAX_ORDER + 1];
  double complex raw_coefficient[MVMC_KRYLOV_MAX_ORDER + 1];
  double maximum_log_raw_magnitude = -INFINITY;
  size_t index;
  int nonzero = 0;
  if (policy == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(policy, 0, sizeof(*policy));
  if (!order_valid(order) || !source_valid(source) ||
      provenance_hash == 0 || scaled_coefficient == NULL ||
      log_basis_scale == NULL ||
      coefficient_count != (size_t)order + 1) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(log_raw_magnitude, 0, sizeof(log_raw_magnitude));
  memset(raw_coefficient, 0, sizeof(raw_coefficient));
  for (index = 0; index < coefficient_count; ++index) {
    const double complex coefficient = scaled_coefficient[index];
    double magnitude;
    double log_magnitude;
    if (!finite_complex(coefficient) || !isfinite(log_basis_scale[index])) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    magnitude = hypot(creal(coefficient), cimag(coefficient));
    if (!isfinite(magnitude)) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    if (magnitude == 0.0) {
      log_raw_magnitude[index] = -INFINITY;
      continue;
    }
    log_magnitude = log(magnitude) + log_basis_scale[index];
    if (!isfinite(log_magnitude)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    log_raw_magnitude[index] = log_magnitude;
    if (log_magnitude > maximum_log_raw_magnitude) {
      maximum_log_raw_magnitude = log_magnitude;
    }
    nonzero = 1;
  }
  if (!nonzero) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  for (index = 0; index < coefficient_count; ++index) {
    const double complex coefficient = scaled_coefficient[index];
    const double magnitude = hypot(creal(coefficient), cimag(coefficient));
    double relative_magnitude;
    double phase_real;
    double phase_imaginary;
    if (magnitude == 0.0) continue;
    relative_magnitude =
        exp(log_raw_magnitude[index] - maximum_log_raw_magnitude);
    if (!isfinite(relative_magnitude)) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    if (relative_magnitude == 0.0) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    phase_real = creal(coefficient) / magnitude;
    phase_imaginary = cimag(coefficient) / magnitude;
    raw_coefficient[index] =
        phase_real * relative_magnitude +
        I * (phase_imaginary * relative_magnitude);
    if (!finite_complex(raw_coefficient[index]) ||
        (creal(raw_coefficient[index]) == 0.0 &&
         cimag(raw_coefficient[index]) == 0.0)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  }
  return mvmc_krylov_final_state_policy_create(
      order, source, provenance_hash, raw_coefficient, coefficient_count,
      policy);
}

uint64_t mvmc_krylov_final_state_policy_hash(
    const MVMCKrylovFinalStatePolicy *policy) {
  return policy_valid(policy) ? policy->policy_hash : UINT64_C(0);
}

MVMCKrylovStatus mvmc_krylov_final_state_evaluate(
    const MVMCKrylovFinalStatePolicy *policy,
    const MVMCKrylovBoundedResult *krylov,
    MVMCKrylovFinalStateEvaluation *result) {
  MVMCKrylovFinalStateEvaluation candidate;
  MVMCScaledComplex amplitude_terms[MVMC_KRYLOV_MAX_ORDER + 1];
  MVMCScaledComplex hamiltonian_terms[MVMC_KRYLOV_MAX_ORDER + 1];
  int index;
  MVMCKrylovStatus status;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_evaluation(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (!policy_valid(policy) || krylov == NULL || !krylov->valid ||
      krylov->status != MVMC_KRYLOV_STATUS_OK ||
      krylov->evaluated_order < policy->order + 1) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.order = policy->order;
  candidate.policy_hash = policy->policy_hash;
  candidate.log_weight = NAN;

  for (index = 0; index <= policy->order + 1; ++index) {
    const MVMCScaledComplex *value = &krylov->value[index];
    if (!mvmc_scaled_complex_is_valid(value)) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    switch (value->state) {
      case MVMC_SCALED_COMPLEX_FINITE_NONZERO:
        ++candidate.basis_finite_count;
        break;
      case MVMC_SCALED_COMPLEX_EXACT_ZERO:
        ++candidate.basis_exact_zero_count;
        break;
      case MVMC_SCALED_COMPLEX_NUMERIC_ZERO:
        ++candidate.basis_numeric_zero_count;
        break;
      case MVMC_SCALED_COMPLEX_NONFINITE:
        invalidate_evaluation(MVMC_KRYLOV_STATUS_NONFINITE, result);
        return MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }

  for (index = 0; index <= policy->order; ++index) {
    MVMCScaledComplex coefficient_value;
    status = scaled_coefficient(policy->coefficient[index],
                                &coefficient_value);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      invalidate_evaluation(status, result);
      return status;
    }
    if (coefficient_value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO) {
      ++candidate.coefficient_zero_count;
    }
    if (mvmc_scaled_complex_multiply(
            &coefficient_value, &krylov->value[index],
            &amplitude_terms[index]) != MVMC_PFAFFIAN_STATUS_OK ||
        mvmc_scaled_complex_multiply(
            &coefficient_value, &krylov->value[index + 1],
            &hamiltonian_terms[index]) != MVMC_PFAFFIAN_STATUS_OK) {
      invalidate_evaluation(
          MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE, result);
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
  }
  if (mvmc_scaled_complex_sum_ordered(
          amplitude_terms, (size_t)policy->order + 1,
          &candidate.amplitude) != MVMC_PFAFFIAN_STATUS_OK ||
      mvmc_scaled_complex_sum_ordered(
          hamiltonian_terms, (size_t)policy->order + 1,
          &candidate.hamiltonian_amplitude) !=
          MVMC_PFAFFIAN_STATUS_OK) {
    invalidate_evaluation(
        MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE, result);
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (candidate.amplitude.state == MVMC_SCALED_COMPLEX_NONFINITE ||
      candidate.hamiltonian_amplitude.state ==
          MVMC_SCALED_COMPLEX_NONFINITE) {
    invalidate_evaluation(MVMC_KRYLOV_STATUS_NONFINITE, result);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  if (candidate.amplitude.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
    candidate.log_weight = 2.0 * candidate.amplitude.log_abs;
    if (!isfinite(candidate.log_weight)) {
      invalidate_evaluation(MVMC_KRYLOV_STATUS_NONFINITE, result);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    status = mvmc_power_lanczos_scaled_divide(
        &candidate.hamiltonian_amplitude, &candidate.amplitude,
        &candidate.local_energy);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      invalidate_evaluation(status, result);
      return status;
    }
    candidate.local_energy_available =
        candidate.local_energy.state ==
            MVMC_SCALED_COMPLEX_FINITE_NONZERO ||
        candidate.local_energy.state == MVMC_SCALED_COMPLEX_EXACT_ZERO;
    candidate.sampleable = candidate.local_energy_available;
  } else if (candidate.amplitude.state ==
             MVMC_SCALED_COMPLEX_EXACT_ZERO) {
    candidate.log_weight = -INFINITY;
  }
  candidate.valid = 1;
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCScaledComplexExportStatus mvmc_krylov_final_state_local_energy_export(
    const MVMCKrylovFinalStateEvaluation *evaluation,
    double complex *local_energy) {
  if (!evaluation_valid(evaluation) ||
      !evaluation->local_energy_available || local_energy == NULL) {
    return MVMC_SCALED_EXPORT_INVALID;
  }
  return mvmc_scaled_complex_export_common_scale(
      &evaluation->local_energy, 0.0, local_energy);
}

MVMCKrylovStatus mvmc_krylov_final_state_local_energy_evidence(
    const MVMCKrylovFinalStateEvaluation *evaluation,
    MVMCPowerLanczosNumericEvidence *evidence) {
  if (!evaluation_valid(evaluation) ||
      !evaluation->local_energy_available || evidence == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return mvmc_power_lanczos_scaled_export_evidence(
      &evaluation->local_energy, 0.0, evidence);
}

MVMCKrylovStatus mvmc_krylov_final_state_acceptance(
    const MVMCKrylovFinalStateEvaluation *current,
    const MVMCKrylovFinalStateEvaluation *proposal,
    double log_proposal_ratio, double uniform,
    MVMCKrylovFinalStateAcceptance *result) {
  MVMCKrylovFinalStateAcceptance candidate;
  double log_target_ratio;
  double log_acceptance_ratio;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_acceptance(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (!evaluation_valid(current) || !evaluation_valid(proposal) ||
      current->order != proposal->order ||
      current->policy_hash != proposal->policy_hash ||
      current->amplitude.state != MVMC_SCALED_COMPLEX_FINITE_NONZERO ||
      !current->sampleable || !current->local_energy_available ||
      !isfinite(current->log_weight) || !isfinite(log_proposal_ratio) ||
      !isfinite(uniform) || uniform < 0.0 || uniform >= 1.0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (proposal->amplitude.state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO ||
      (proposal->amplitude.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
       (!proposal->sampleable || !proposal->local_energy_available))) {
    invalidate_acceptance(MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE, result);
    return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  if (proposal->amplitude.state == MVMC_SCALED_COMPLEX_NONFINITE) {
    invalidate_acceptance(MVMC_KRYLOV_STATUS_NONFINITE, result);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.uniform = uniform;
  candidate.log_proposal_ratio = log_proposal_ratio;
  if (proposal->amplitude.state == MVMC_SCALED_COMPLEX_EXACT_ZERO) {
    candidate.valid = 1;
    candidate.exact_zero_proposal = 1;
    candidate.log_target_ratio = -INFINITY;
    candidate.log_acceptance_ratio = -INFINITY;
    *result = candidate;
    return MVMC_KRYLOV_STATUS_OK;
  }
  if (!isfinite(proposal->log_weight)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  log_target_ratio = proposal->log_weight - current->log_weight;
  log_acceptance_ratio = log_target_ratio + log_proposal_ratio;
  if (!isfinite(log_target_ratio) || !isfinite(log_acceptance_ratio)) {
    invalidate_acceptance(MVMC_KRYLOV_STATUS_NONFINITE, result);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  candidate.valid = 1;
  candidate.log_target_ratio = log_target_ratio;
  candidate.log_acceptance_ratio = log_acceptance_ratio;
  candidate.accepted = log_acceptance_ratio >= 0.0 || uniform == 0.0 ||
                       log(uniform) < log_acceptance_ratio;
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus evaluate_krylov(
    MVMCKrylovBoundedWorkspace *workspace,
    const uint64_t *words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovBoundedResult *result) {
  if (mvmc_bounded_krylov_session_is_active(workspace)) {
    return mvmc_bounded_krylov_session_evaluate_bound(
        workspace, words, word_count, amplitude, amplitude_context, result);
  }
  return mvmc_bounded_krylov_evaluate(
      workspace, words, word_count, amplitude, amplitude_context, result);
}

static MVMCKrylovStatus evaluate_snapshot(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy,
    const uint64_t *words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovFinalStateSnapshot *snapshot) {
  MVMCKrylovFinalStateSnapshot candidate;
  MVMCKrylovStatus status;
  if (snapshot == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_snapshot(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, snapshot);
  if (workspace == NULL || !policy_valid(policy) || words == NULL ||
      word_count == 0 || amplitude == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  status = evaluate_krylov(workspace, words, word_count, amplitude,
                           amplitude_context, &candidate.krylov);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_snapshot(status, snapshot);
    return status;
  }
  status = mvmc_krylov_final_state_evaluate(
      policy, &candidate.krylov, &candidate.final_state);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_snapshot(status, snapshot);
    return status;
  }
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.word_count = word_count;
  candidate.policy_hash = policy->policy_hash;
  candidate.configuration_hash = configuration_hash(words, word_count);
  candidate.amplitude_generation_hash =
      candidate.krylov.statistics.amplitude_generation_hash;
  candidate.integrity_hash =
      mvmc_krylov_final_state_snapshot_integrity_hash(&candidate);
  if (candidate.integrity_hash == 0) {
    invalidate_snapshot(
        MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE, snapshot);
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  *snapshot = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_final_state_sampler_initialize(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy,
    const uint64_t *current_words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovFinalStateSnapshot *snapshot) {
  MVMCKrylovStatus status = evaluate_snapshot(
      workspace, policy, current_words, word_count, amplitude,
      amplitude_context, snapshot);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (!snapshot->final_state.sampleable ||
      !snapshot->final_state.local_energy_available) {
    invalidate_snapshot(MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE, snapshot);
    return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_final_state_sampler_step(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy,
    uint64_t *current_words, size_t current_word_count,
    MVMCKrylovFinalStateSnapshot *current,
    const uint64_t *proposal_words, size_t proposal_word_count,
    double log_proposal_ratio, double uniform,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovFinalStateStepResult *result) {
  MVMCKrylovFinalStateSnapshot proposal;
  MVMCKrylovFinalStateSnapshot updated;
  MVMCKrylovFinalStateStepResult candidate;
  MVMCKrylovFinalStateAcceptance acceptance;
  MVMCKrylovStatus status;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_step(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (workspace == NULL || !policy_valid(policy) || current_words == NULL ||
      current == NULL || !current->valid ||
      current->status != MVMC_KRYLOV_STATUS_OK ||
      current->word_count != current_word_count ||
      current->policy_hash != policy->policy_hash ||
      current->integrity_hash == 0 ||
      current->integrity_hash !=
          mvmc_krylov_final_state_snapshot_integrity_hash(current) ||
      current->final_state.policy_hash != policy->policy_hash ||
      current->final_state.order != policy->order ||
      current->configuration_hash !=
          configuration_hash(current_words, current_word_count) ||
      !current->final_state.sampleable ||
      !current->final_state.local_energy_available ||
      proposal_words == NULL ||
      word_ranges_overlap(current_words, current_word_count,
                          proposal_words, proposal_word_count) ||
      proposal_word_count != current_word_count || current_word_count == 0 ||
      current_word_count > SIZE_MAX / sizeof(current_words[0]) ||
      !isfinite(log_proposal_ratio) || !isfinite(uniform) ||
      uniform < 0.0 || uniform >= 1.0 || amplitude == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = evaluate_snapshot(workspace, policy, proposal_words,
                             proposal_word_count, amplitude,
                             amplitude_context, &proposal);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_step(status, result);
    return status;
  }
  if (proposal.amplitude_generation_hash !=
      current->amplitude_generation_hash) {
    invalidate_step(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_final_state_acceptance(
      &current->final_state, &proposal.final_state, log_proposal_ratio,
      uniform, &acceptance);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_step(status, result);
    return status;
  }
  if (acceptance.accepted && current->accepted_generation == UINT64_MAX) {
    invalidate_step(MVMC_KRYLOV_STATUS_RESOURCE_LIMIT, result);
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }

  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.accepted = acceptance.accepted;
  candidate.proposal_krylov = proposal.krylov;
  candidate.proposal_final_state = proposal.final_state;
  candidate.acceptance = acceptance;
  updated = *current;
  if (acceptance.accepted) {
    updated.krylov = proposal.krylov;
    updated.final_state = proposal.final_state;
    updated.configuration_hash = proposal.configuration_hash;
    updated.amplitude_generation_hash =
        proposal.amplitude_generation_hash;
    updated.integrity_hash = proposal.integrity_hash;
    ++updated.accepted_generation;
    memcpy(current_words, proposal_words,
           current_word_count * sizeof(current_words[0]));
    *current = updated;
  }
  candidate.accepted_generation = current->accepted_generation;
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_final_state_sampler_step_selected_neighbor(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy,
    const MVMCKrylovFockModel *proposal_model,
    uint64_t *current_words, size_t current_word_count,
    MVMCKrylovFinalStateSnapshot *current,
    size_t selected_neighbor_index, double uniform,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words, size_t proposal_word_count,
    MVMCKrylovFinalStateProposalStepResult *result) {
  MVMCKrylovStatus status;
  MVMCKrylovFinalStateStepResult step;
  size_t neighbor_count = 0;
  double log_proposal_ratio = NAN;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_proposal_step(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (current_words == NULL || current_word_count == 0 ||
      proposal_words == NULL ||
      word_ranges_overlap(current_words, current_word_count,
                          proposal_words, proposal_word_count) ||
      proposal_word_count != current_word_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_fock_proposal_select_neighbor(
      proposal_model, current_words, current_word_count,
      selected_neighbor_index, proposal_words, proposal_word_count,
      &neighbor_count);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_proposal_step(status, result);
    return status;
  }
  status = mvmc_krylov_fock_proposal_log_ratio(
      proposal_model, current_words, current_word_count,
      proposal_words, proposal_word_count, &log_proposal_ratio);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_proposal_step(status, result);
    return status;
  }
  status = mvmc_krylov_final_state_sampler_step(
      workspace, policy, current_words, current_word_count, current,
      proposal_words, proposal_word_count, log_proposal_ratio, uniform,
      amplitude, amplitude_context, &step);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_proposal_step(status, result);
    return status;
  }
  result->valid = 1;
  result->status = MVMC_KRYLOV_STATUS_OK;
  result->configuration_changed = step.accepted;
  result->neighbor_count = neighbor_count;
  result->selected_neighbor_index = selected_neighbor_index;
  result->log_proposal_ratio = log_proposal_ratio;
  result->step = step;
  return MVMC_KRYLOV_STATUS_OK;
}
