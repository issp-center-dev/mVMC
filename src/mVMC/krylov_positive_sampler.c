/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "krylov_positive_sampler.h"

#if !defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                         \
    !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "krylov_positive_sampler.c is Testing-only"
#endif

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double sampler_wall_seconds(void) {
  struct timespec now;
  if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return -1.0;
  return (double)now.tv_sec + 1.0e-9 * (double)now.tv_nsec;
}

static double sampler_elapsed_seconds(double start) {
  const double end = sampler_wall_seconds();
  return start >= 0.0 && end >= start ? end - start : -1.0;
}

static void hash_byte(uint64_t *hash, unsigned char value) {
  *hash ^= (uint64_t)value;
  *hash *= UINT64_C(1099511628211);
}

static void hash_u64(uint64_t *hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    hash_byte(hash, (unsigned char)(value & UINT64_C(0xff)));
    value >>= 8;
  }
}

static uint64_t double_bits(double value) {
  uint64_t bits = 0;
  if (sizeof(bits) == sizeof(value)) {
    memcpy(&bits, &value, sizeof(bits));
  }
  return bits;
}

static int increment_u64(uint64_t *value) {
  if (value == NULL || *value == UINT64_MAX) return 0;
  ++(*value);
  return 1;
}

static int add_u64(uint64_t *value, uint64_t increment) {
  if (value == NULL || increment > UINT64_MAX - *value) return 0;
  *value += increment;
  return 1;
}

static uint64_t hash_configuration_words(const uint64_t *words,
                                         size_t word_count) {
  uint64_t hash = UINT64_C(1469598103934665603);
  size_t word;
  hash_u64(&hash, (uint64_t)word_count);
  for (word = 0; word < word_count; ++word) {
    hash_u64(&hash, words[word]);
  }
  return hash;
}

static uint64_t hash_guide_policy(
    const MVMCKrylovPositiveGuidePolicy *policy) {
  uint64_t hash = UINT64_C(1469598103934665603);
  int order;
  hash_u64(&hash, policy->policy_hash);
  hash_u64(&hash, (uint64_t)(unsigned int)policy->order);
  hash_u64(&hash, double_bits(policy->eta));
  for (order = 0; order <= policy->order; ++order) {
    hash_u64(&hash, double_bits(policy->lambda[order]));
    hash_u64(&hash, double_bits(policy->log_basis_scale[order]));
  }
  return hash;
}

static uint64_t hash_surrogate_policy_values(
    const MVMCKrylovPositiveSamplerSurrogatePolicy *policy) {
  uint64_t hash = UINT64_C(1469598103934665603);
  hash_u64(&hash, policy->version);
  hash_u64(&hash, policy->guide_shape_hash);
  hash_u64(&hash, policy->formula_version);
  hash_u64(&hash, policy->transaction_version);
  hash_u64(&hash, (uint64_t)policy->step_count);
  hash_u64(&hash, (uint64_t)(unsigned int)policy->partial_order);
  hash_u64(&hash, (uint64_t)policy->floor_multiplier);
  hash_u64(&hash, policy->rng_algorithm);
  return hash;
}

static uint64_t hash_scaled_values(
    const MVMCScaledComplex *values, int evaluated_order) {
  uint64_t hash = UINT64_C(1469598103934665603);
  int order;
  hash_u64(&hash, (uint64_t)(unsigned int)evaluated_order);
  for (order = 0; order <= evaluated_order; ++order) {
    const MVMCScaledComplex *value = &values[order];
    hash_u64(&hash, (uint64_t)(unsigned int)value->state);
    hash_u64(&hash, double_bits(creal(value->phase)));
    hash_u64(&hash, double_bits(cimag(value->phase)));
    hash_u64(&hash, double_bits(value->log_abs));
    hash_u64(&hash, double_bits(value->log_abs_error_bound));
    hash_u64(&hash, double_bits(value->max_input_log_abs));
    hash_u64(&hash, double_bits(value->cancellation_log_abs));
    hash_u64(&hash, double_bits(value->cancellation_ratio));
  }
  return hash;
}

static int rng_valid(const MVMCKrylovPositiveSamplerRng *rng) {
  return rng != NULL && rng->valid &&
         rng->state_version == MVMC_KRYLOV_POSITIVE_SAMPLER_STATE_VERSION &&
         rng->algorithm == MVMC_KRYLOV_POSITIVE_SAMPLER_RNG_SPLITMIX64;
}

static uint64_t splitmix64_scramble(uint64_t value) {
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

static void invalidate_snapshot(MVMCKrylovStatus status,
                                MVMCKrylovPositiveSamplerSnapshot *snapshot) {
  if (snapshot == NULL) return;
  memset(snapshot, 0, sizeof(*snapshot));
  snapshot->status = status;
}

static void invalidate_step(MVMCKrylovStatus status,
                            MVMCKrylovPositiveSamplerStepResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
}

static void invalidate_proposal_step(
    MVMCKrylovStatus status,
    MVMCKrylovPositiveSamplerProposalStepResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
}

static void invalidate_proposal_draw(
    MVMCKrylovStatus status,
    MVMCKrylovPositiveSamplerProposalDrawResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
}

static void invalidate_surrogate_weight(
    MVMCKrylovStatus status,
    MVMCKrylovPositiveSamplerSurrogateWeightResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
  result->log_weight = NAN;
  result->log_floor = NAN;
  result->log_signal = NAN;
}

static void invalidate_surrogate_draw(
    MVMCKrylovStatus status,
    MVMCKrylovPositiveSamplerSurrogateDrawResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
  result->initial_log_weight = NAN;
  result->final_log_weight = NAN;
  result->log_proposal_ratio = NAN;
  result->proposal_seconds = -1.0;
}

static void invalidate_manifest(
    MVMCKrylovStatus status,
    MVMCKrylovPositiveSamplerManifest *manifest) {
  if (manifest == NULL) return;
  memset(manifest, 0, sizeof(*manifest));
  manifest->status = status;
}

static void invalidate_trace_statistics(
    MVMCKrylovStatus status,
    MVMCKrylovPositiveSamplerTraceStatistics *statistics) {
  if (statistics == NULL) return;
  memset(statistics, 0, sizeof(*statistics));
  statistics->status = status;
  statistics->min_log_acceptance_ratio = NAN;
  statistics->max_log_acceptance_ratio = NAN;
  statistics->sum_log_acceptance_ratio = NAN;
  statistics->min_proposal_log_guide = NAN;
  statistics->max_proposal_log_guide = NAN;
}

static int policy_shape_valid(const MVMCKrylovPositiveGuidePolicy *policy) {
  return policy != NULL && policy->order >= 0 &&
         policy->order <= MVMC_KRYLOV_MAX_ORDER;
}

static int policy_value_valid(const MVMCKrylovPositiveGuidePolicy *policy) {
  int order;
  if (!policy_shape_valid(policy) || !isfinite(policy->eta) ||
      policy->eta <= 0.0) {
    return 0;
  }
  for (order = 0; order <= policy->order; ++order) {
    if (!isfinite(policy->lambda[order]) || policy->lambda[order] < 0.0 ||
        !isfinite(policy->log_basis_scale[order])) {
      return 0;
    }
  }
  return 1;
}

static uint64_t hash_proposal_policy_values(
    const MVMCKrylovPositiveSamplerProposalPolicy *policy) {
  uint64_t hash = UINT64_C(1469598103934665603);
  hash_u64(&hash, policy->version);
  hash_u64(&hash, policy->kernel_id);
  hash_u64(&hash, (uint64_t)policy->global_numerator);
  hash_u64(&hash, (uint64_t)policy->global_denominator);
  hash_u64(&hash, (uint64_t)policy->neighbor_numerator);
  hash_u64(&hash, (uint64_t)policy->neighbor_denominator);
  hash_u64(&hash, (uint64_t)policy->distance_numerator);
  hash_u64(&hash, (uint64_t)policy->distance_denominator);
  hash_u64(&hash, policy->distance_rounding_rule);
  hash_u64(&hash, policy->subset_order);
  hash_u64(&hash, policy->self_report_rule);
  hash_u64(&hash, policy->rng_algorithm);
  return hash;
}

static uint64_t hash_proposal_model_values(
    const MVMCKrylovFockModel *proposal_model) {
  uint64_t hash = UINT64_C(1469598103934665603);
  hash_u64(&hash, MVMC_KRYLOV_PROPOSAL_MODEL_CONTRACT_VERSION);
  hash_u64(&hash, (uint64_t)proposal_model->site_count);
  hash_u64(&hash, (uint64_t)proposal_model->up_electron_count);
  hash_u64(&hash, (uint64_t)proposal_model->down_electron_count);
  hash_u64(&hash, (uint64_t)(unsigned int)proposal_model->pure_spin);
  return hash;
}

static int proposal_policy_valid(
    const MVMCKrylovPositiveSamplerProposalPolicy *policy) {
  int kernel_valid;
  if (policy == NULL || !policy->valid ||
      policy->status != MVMC_KRYLOV_STATUS_OK ||
      policy->version != MVMC_KRYLOV_PROPOSAL_POLICY_VERSION ||
      policy->subset_order != MVMC_KRYLOV_PROPOSAL_SUBSET_CANONICAL_V1 ||
      policy->self_report_rule != MVMC_KRYLOV_PROPOSAL_SELF_REPORT_V1 ||
      policy->rng_algorithm != MVMC_KRYLOV_POSITIVE_SAMPLER_RNG_SPLITMIX64 ||
      policy->policy_hash != hash_proposal_policy_values(policy)) {
    return 0;
  }
  if (policy->kernel_id ==
      MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_GLOBAL_MIXTURE) {
    kernel_valid =
        policy->global_denominator != 0 &&
        policy->global_numerator <= policy->global_denominator &&
        (size_t)((uint64_t)policy->global_numerator) ==
            policy->global_numerator &&
        (size_t)((uint64_t)policy->global_denominator) ==
            policy->global_denominator &&
        policy->neighbor_numerator == 0 &&
        policy->neighbor_denominator == 0 &&
        policy->distance_numerator == 0 &&
        policy->distance_denominator == 0 &&
        policy->distance_rounding_rule == 0;
  } else if (policy->kernel_id ==
             MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_SHELL_MIXTURE) {
    kernel_valid =
        policy->global_numerator == 0 &&
        policy->global_denominator == 0 &&
        policy->neighbor_denominator != 0 &&
        policy->neighbor_numerator != 0 &&
        policy->neighbor_numerator < policy->neighbor_denominator &&
        policy->distance_denominator != 0 &&
        policy->distance_numerator != 0 &&
        policy->distance_numerator <= policy->distance_denominator &&
        (size_t)((uint64_t)policy->neighbor_numerator) ==
            policy->neighbor_numerator &&
        (size_t)((uint64_t)policy->neighbor_denominator) ==
            policy->neighbor_denominator &&
        (size_t)((uint64_t)policy->distance_numerator) ==
            policy->distance_numerator &&
        (size_t)((uint64_t)policy->distance_denominator) ==
            policy->distance_denominator &&
        policy->distance_rounding_rule ==
            MVMC_KRYLOV_PROPOSAL_DISTANCE_ROUND_HALF_UP_V1;
  } else {
    kernel_valid = 0;
  }
  return kernel_valid;
}

static int surrogate_policy_valid(
    const MVMCKrylovPositiveSamplerSurrogatePolicy *policy) {
  return policy != NULL && policy->valid &&
         policy->status == MVMC_KRYLOV_STATUS_OK &&
         policy->version == MVMC_KRYLOV_SURROGATE_POLICY_VERSION &&
         policy->policy_hash != 0 && policy->guide_shape_hash != 0 &&
         policy->formula_version == MVMC_KRYLOV_SURROGATE_Q_PARTIAL_FLOOR_V1 &&
         policy->transaction_version == MVMC_KRYLOV_SURROGATE_TRANSACTION_V1 &&
         policy->step_count > 0 && policy->partial_order >= 0 &&
         policy->partial_order <= MVMC_KRYLOV_MAX_ORDER &&
         policy->floor_multiplier > 0 &&
         (size_t)((uint64_t)policy->step_count) == policy->step_count &&
         (size_t)((uint64_t)policy->floor_multiplier) ==
             policy->floor_multiplier &&
         policy->rng_algorithm == MVMC_KRYLOV_POSITIVE_SAMPLER_RNG_SPLITMIX64 &&
         policy->policy_hash == hash_surrogate_policy_values(policy);
}

static int bounded_limits_valid(const MVMCKrylovBoundedLimits *limits) {
  return limits != NULL && limits->max_order >= 0 &&
         limits->max_order <= MVMC_KRYLOV_MAX_ORDER &&
         limits->max_workspace_bytes != 0 &&
         limits->max_node_expansions != 0 &&
         limits->max_terminal_amplitude_calls != 0 &&
         limits->max_total_row_transitions != 0;
}

static int manifest_equal(
    const MVMCKrylovPositiveSamplerManifest *left,
    const MVMCKrylovPositiveSamplerManifest *right) {
  return left->valid == right->valid &&
         left->status == right->status &&
         left->state_version == right->state_version &&
         left->rng_algorithm == right->rng_algorithm &&
         left->rng_state == right->rng_state &&
         left->rng_stream == right->rng_stream &&
         left->rng_draws == right->rng_draws &&
         left->policy_hash == right->policy_hash &&
         left->proposal_policy_hash == right->proposal_policy_hash &&
         left->proposal_model_hash == right->proposal_model_hash &&
         left->surrogate_policy_hash == right->surrogate_policy_hash &&
         left->guide_shape_hash == right->guide_shape_hash &&
         left->bounded_plan_hash == right->bounded_plan_hash &&
         left->amplitude_policy_hash == right->amplitude_policy_hash &&
         left->current_configuration_hash ==
             right->current_configuration_hash &&
         left->current_scale_hash == right->current_scale_hash &&
         left->accepted_generation == right->accepted_generation &&
         left->word_count == right->word_count &&
         left->cache_bytes == right->cache_bytes &&
         left->guide_order == right->guide_order &&
         left->bounded_max_order == right->bounded_max_order;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_rng_seed(
    uint64_t seed, uint64_t stream,
    MVMCKrylovPositiveSamplerRng *rng) {
  MVMCKrylovPositiveSamplerRng candidate;
  if (rng == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.state_version = MVMC_KRYLOV_POSITIVE_SAMPLER_STATE_VERSION;
  candidate.algorithm = MVMC_KRYLOV_POSITIVE_SAMPLER_RNG_SPLITMIX64;
  candidate.stream = stream;
  candidate.state = seed ^ splitmix64_scramble(stream + UINT64_C(0x9e3779b97f4a7c15));
  candidate.draws = 0;
  *rng = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_rng_draw_uint64(
    MVMCKrylovPositiveSamplerRng *rng, uint64_t *value) {
  MVMCKrylovPositiveSamplerRng candidate;
  if (!rng_valid(rng) || value == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (rng->draws == UINT64_MAX) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  candidate = *rng;
  candidate.state += UINT64_C(0x9e3779b97f4a7c15);
  *value = splitmix64_scramble(candidate.state);
  ++candidate.draws;
  *rng = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_rng_draw_bounded(
    MVMCKrylovPositiveSamplerRng *rng, size_t bound, size_t *value) {
  uint64_t uint_bound;
  uint64_t threshold;
  uint64_t drawn = 0;
  if (!rng_valid(rng) || value == NULL || bound == 0 ||
      (size_t)((uint64_t)bound) != bound) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  uint_bound = (uint64_t)bound;
  threshold = (uint64_t)(-uint_bound) % uint_bound;
  for (;;) {
    MVMCKrylovStatus status =
        mvmc_krylov_positive_sampler_rng_draw_uint64(rng, &drawn);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    if (drawn >= threshold) break;
  }
  *value = (size_t)(drawn % uint_bound);
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_rng_draw_uniform(
    MVMCKrylovPositiveSamplerRng *rng, double *uniform) {
  uint64_t drawn = 0;
  MVMCKrylovStatus status;
  if (!rng_valid(rng) || uniform == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_positive_sampler_rng_draw_uint64(rng, &drawn);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  *uniform = (double)(drawn >> 11) * (1.0 / 9007199254740992.0);
  return isfinite(*uniform) && *uniform >= 0.0 && *uniform < 1.0
             ? MVMC_KRYLOV_STATUS_OK
             : MVMC_KRYLOV_STATUS_NONFINITE;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_proposal_policy_create(
    size_t global_numerator, size_t global_denominator,
    MVMCKrylovPositiveSamplerProposalPolicy *policy) {
  MVMCKrylovPositiveSamplerProposalPolicy candidate;
  if (policy == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(policy, 0, sizeof(*policy));
  policy->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (global_denominator == 0 || global_numerator > global_denominator ||
      (size_t)((uint64_t)global_numerator) != global_numerator ||
      (size_t)((uint64_t)global_denominator) != global_denominator) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.version = MVMC_KRYLOV_PROPOSAL_POLICY_VERSION;
  candidate.kernel_id =
      MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_GLOBAL_MIXTURE;
  candidate.global_numerator = global_numerator;
  candidate.global_denominator = global_denominator;
  candidate.subset_order = MVMC_KRYLOV_PROPOSAL_SUBSET_CANONICAL_V1;
  candidate.self_report_rule = MVMC_KRYLOV_PROPOSAL_SELF_REPORT_V1;
  candidate.rng_algorithm = MVMC_KRYLOV_POSITIVE_SAMPLER_RNG_SPLITMIX64;
  candidate.policy_hash = hash_proposal_policy_values(&candidate);
  *policy = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_shell_policy_create(
    size_t neighbor_numerator, size_t neighbor_denominator,
    size_t distance_numerator, size_t distance_denominator,
    MVMCKrylovPositiveSamplerProposalPolicy *policy) {
  MVMCKrylovPositiveSamplerProposalPolicy candidate;
  if (policy == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(policy, 0, sizeof(*policy));
  policy->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (neighbor_denominator == 0 || neighbor_numerator == 0 ||
      neighbor_numerator >= neighbor_denominator ||
      distance_denominator == 0 || distance_numerator == 0 ||
      distance_numerator > distance_denominator ||
      (size_t)((uint64_t)neighbor_numerator) != neighbor_numerator ||
      (size_t)((uint64_t)neighbor_denominator) != neighbor_denominator ||
      (size_t)((uint64_t)distance_numerator) != distance_numerator ||
      (size_t)((uint64_t)distance_denominator) != distance_denominator) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.version = MVMC_KRYLOV_PROPOSAL_POLICY_VERSION;
  candidate.kernel_id =
      MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_SHELL_MIXTURE;
  candidate.neighbor_numerator = neighbor_numerator;
  candidate.neighbor_denominator = neighbor_denominator;
  candidate.distance_numerator = distance_numerator;
  candidate.distance_denominator = distance_denominator;
  candidate.distance_rounding_rule =
      MVMC_KRYLOV_PROPOSAL_DISTANCE_ROUND_HALF_UP_V1;
  candidate.subset_order = MVMC_KRYLOV_PROPOSAL_SUBSET_CANONICAL_V1;
  candidate.self_report_rule = MVMC_KRYLOV_PROPOSAL_SELF_REPORT_V1;
  candidate.rng_algorithm = MVMC_KRYLOV_POSITIVE_SAMPLER_RNG_SPLITMIX64;
  candidate.policy_hash = hash_proposal_policy_values(&candidate);
  *policy = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_surrogate_policy_create(
    size_t step_count, size_t floor_multiplier,
    const MVMCKrylovPositiveGuidePolicy *guide_policy,
    MVMCKrylovPositiveSamplerSurrogatePolicy *policy) {
  return mvmc_krylov_positive_sampler_surrogate_partial_policy_create(
      step_count, 0, floor_multiplier, guide_policy, policy);
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_surrogate_partial_policy_create(
    size_t step_count, int partial_order, size_t floor_multiplier,
    const MVMCKrylovPositiveGuidePolicy *guide_policy,
    MVMCKrylovPositiveSamplerSurrogatePolicy *policy) {
  MVMCKrylovPositiveSamplerSurrogatePolicy candidate;
  if (policy == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(policy, 0, sizeof(*policy));
  policy->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (step_count == 0 || step_count == SIZE_MAX || partial_order < 0 ||
      partial_order > MVMC_KRYLOV_MAX_ORDER || floor_multiplier == 0 ||
      (size_t)((uint64_t)step_count) != step_count ||
      (size_t)((uint64_t)floor_multiplier) != floor_multiplier ||
      (uint64_t)floor_multiplier > UINT64_C(9007199254740992) ||
      !policy_value_valid(guide_policy) ||
      partial_order > guide_policy->order) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.version = MVMC_KRYLOV_SURROGATE_POLICY_VERSION;
  candidate.guide_shape_hash = hash_guide_policy(guide_policy);
  candidate.formula_version = MVMC_KRYLOV_SURROGATE_Q_PARTIAL_FLOOR_V1;
  candidate.transaction_version = MVMC_KRYLOV_SURROGATE_TRANSACTION_V1;
  candidate.step_count = step_count;
  candidate.partial_order = partial_order;
  candidate.floor_multiplier = floor_multiplier;
  candidate.rng_algorithm = MVMC_KRYLOV_POSITIVE_SAMPLER_RNG_SPLITMIX64;
  candidate.policy_hash = hash_surrogate_policy_values(&candidate);
  *policy = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_surrogate_weight_zeroth(
    const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy,
    const MVMCKrylovPositiveGuidePolicy *guide_policy,
    const MVMCScaledComplex *zeroth_value,
    MVMCKrylovPositiveSamplerSurrogateWeightResult *result) {
  if (surrogate_policy == NULL || surrogate_policy->partial_order != 0) {
    if (result != NULL) {
      invalidate_surrogate_weight(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                                  result);
    }
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return mvmc_krylov_positive_sampler_surrogate_weight_partial(
      surrogate_policy, guide_policy, zeroth_value, 1, result);
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_surrogate_weight_partial(
    const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy,
    const MVMCKrylovPositiveGuidePolicy *guide_policy,
    const MVMCScaledComplex *values, size_t value_count,
    MVMCKrylovPositiveSamplerSurrogateWeightResult *result) {
  MVMCKrylovPositiveSamplerSurrogateWeightResult candidate;
  int order;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_surrogate_weight(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (!surrogate_policy_valid(surrogate_policy) ||
      !policy_value_valid(guide_policy) || values == NULL ||
      value_count <= (size_t)surrogate_policy->partial_order ||
      surrogate_policy->guide_shape_hash != hash_guide_policy(guide_policy) ||
      surrogate_policy->partial_order > guide_policy->order) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.log_floor =
      log((double)surrogate_policy->floor_multiplier) +
      log(guide_policy->eta);
  candidate.log_signal = -INFINITY;
  candidate.floor_only = 1;
  if (!isfinite(candidate.log_floor)) {
    invalidate_surrogate_weight(MVMC_KRYLOV_STATUS_NONFINITE, result);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  for (order = 0; order <= surrogate_policy->partial_order; ++order) {
    const MVMCScaledComplex *value = &values[order];
    double signal;
    double maximum;
    if (!mvmc_scaled_complex_is_valid(value)) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    if (value->state == MVMC_SCALED_COMPLEX_NONFINITE) {
      invalidate_surrogate_weight(MVMC_KRYLOV_STATUS_NONFINITE, result);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    if (value->state == MVMC_SCALED_COMPLEX_EXACT_ZERO ||
        value->state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO) {
      continue;
    }
    signal = 2.0 * (value->log_abs +
                    guide_policy->log_basis_scale[order]);
    if (!isfinite(signal)) {
      invalidate_surrogate_weight(MVMC_KRYLOV_STATUS_NONFINITE, result);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    maximum = fmax(candidate.log_signal, signal);
    candidate.log_signal = isinf(candidate.log_signal) &&
                                   candidate.log_signal < 0.0
                               ? signal
                               : maximum + log1p(exp(
                                     fmin(candidate.log_signal, signal) -
                                     maximum));
    candidate.floor_only = 0;
  }
  if (candidate.floor_only) {
    candidate.log_weight = candidate.log_floor;
  } else {
    const double maximum = fmax(candidate.log_floor, candidate.log_signal);
    candidate.log_weight = maximum + log1p(exp(
        fmin(candidate.log_floor, candidate.log_signal) - maximum));
  }
  if (!isfinite(candidate.log_weight)) {
    invalidate_surrogate_weight(MVMC_KRYLOV_STATUS_NONFINITE, result);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  candidate.valid = 1;
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_proposal_model_hash(
    const MVMCKrylovFockModel *proposal_model,
    const uint64_t *configuration_words, size_t word_count,
    uint64_t *proposal_model_hash) {
  MVMCKrylovStatus status;
  size_t neighbor_count = 0;
  if (proposal_model_hash == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *proposal_model_hash = 0;
  status = mvmc_krylov_fock_proposal_count_neighbors(
      proposal_model, configuration_words, word_count, &neighbor_count);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  (void)neighbor_count;
  *proposal_model_hash = hash_proposal_model_values(proposal_model);
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_proposal_model_policy_hash(
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const MVMCKrylovFockModel *proposal_model,
    const uint64_t *configuration_words, size_t word_count,
    uint64_t *proposal_model_hash) {
  MVMCKrylovStatus status;
  uint64_t hash;
  size_t maximum = 0;
  size_t distance = 0;
  if (proposal_model_hash == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *proposal_model_hash = 0;
  if (!proposal_policy_valid(proposal_policy)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_positive_sampler_proposal_model_hash(
      proposal_model, configuration_words, word_count, &hash);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (proposal_policy->kernel_id ==
      MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_SHELL_MIXTURE) {
    status = mvmc_krylov_fock_proposal_resolve_shell_distance(
        proposal_model, proposal_policy->distance_numerator,
        proposal_policy->distance_denominator, &maximum, &distance);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    hash_u64(&hash, proposal_policy->kernel_id);
    hash_u64(&hash, (uint64_t)proposal_policy->neighbor_numerator);
    hash_u64(&hash, (uint64_t)proposal_policy->neighbor_denominator);
    hash_u64(&hash, (uint64_t)proposal_policy->distance_numerator);
    hash_u64(&hash, (uint64_t)proposal_policy->distance_denominator);
    hash_u64(&hash, proposal_policy->distance_rounding_rule);
    hash_u64(&hash, (uint64_t)maximum);
    hash_u64(&hash, (uint64_t)distance);
  }
  *proposal_model_hash = hash;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus evaluate_snapshot(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *policy,
    const uint64_t *words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovPositiveSamplerSnapshot *snapshot) {
  MVMCKrylovStatus status;
  MVMCKrylovPositiveSamplerSnapshot candidate;

  if (snapshot == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_snapshot(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, snapshot);
  if (workspace == NULL || !policy_shape_valid(policy) || words == NULL ||
      word_count == 0 || amplitude == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }

  memset(&candidate, 0, sizeof(candidate));
  status = mvmc_bounded_krylov_session_is_active(workspace)
               ? mvmc_bounded_krylov_session_evaluate_bound(
                     workspace, words, word_count, amplitude,
                     amplitude_context, &candidate.krylov)
               : mvmc_bounded_krylov_evaluate(
                     workspace, words, word_count, amplitude,
                     amplitude_context, &candidate.krylov);
  if (status != MVMC_KRYLOV_STATUS_OK || !candidate.krylov.valid ||
      candidate.krylov.evaluated_order < policy->order) {
    invalidate_snapshot(status == MVMC_KRYLOV_STATUS_OK
                            ? MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE
                            : status,
                        snapshot);
    return snapshot->status;
  }

  status = mvmc_krylov_positive_guide_evaluate(
      policy, candidate.krylov.value, (size_t)policy->order + 1,
      &candidate.guide);
  if (status != MVMC_KRYLOV_STATUS_OK || !candidate.guide.valid) {
    invalidate_snapshot(status, snapshot);
    return status;
  }

  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.word_count = word_count;
  candidate.policy_hash = policy->policy_hash;
  candidate.accepted_generation = 0;
  *snapshot = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_initialize(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *policy,
    const uint64_t *current_words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovPositiveSamplerSnapshot *snapshot) {
  return evaluate_snapshot(workspace, policy, current_words, word_count,
                           amplitude, amplitude_context, snapshot);
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_step(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *policy,
    uint64_t *current_words, size_t current_word_count,
    MVMCKrylovPositiveSamplerSnapshot *current,
    const uint64_t *proposal_words, size_t proposal_word_count,
    double log_proposal_ratio, double uniform,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovPositiveSamplerStepResult *result) {
  MVMCKrylovStatus status;
  MVMCKrylovPositiveSamplerSnapshot proposal;
  MVMCKrylovPositiveGuideAcceptance acceptance;
  MVMCKrylovPositiveSamplerStepResult candidate;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_step(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (workspace == NULL || !policy_shape_valid(policy) ||
      current_words == NULL || current == NULL || !current->valid ||
      current->status != MVMC_KRYLOV_STATUS_OK ||
      current->word_count != current_word_count ||
      current->policy_hash != policy->policy_hash ||
      proposal_words == NULL || proposal_word_count != current_word_count ||
      current_word_count == 0 ||
      current_word_count > SIZE_MAX / sizeof(uint64_t) ||
      !isfinite(log_proposal_ratio) || !isfinite(uniform) ||
      amplitude == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }

  status = evaluate_snapshot(workspace, policy, proposal_words,
                             proposal_word_count, amplitude,
                             amplitude_context, &proposal);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_step(status, result);
    return status;
  }
  status = mvmc_krylov_positive_guide_acceptance(
      &current->guide, &proposal.guide, log_proposal_ratio, uniform,
      &acceptance);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_step(status, result);
    return status;
  }
  if (acceptance.accepted &&
      current->accepted_generation == UINT64_MAX) {
    invalidate_step(MVMC_KRYLOV_STATUS_RESOURCE_LIMIT, result);
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }

  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.accepted = acceptance.accepted;
  candidate.proposal_krylov = proposal.krylov;
  candidate.proposal_guide = proposal.guide;
  candidate.acceptance = acceptance;

  if (acceptance.accepted) {
    memcpy(current_words, proposal_words,
           current_word_count * sizeof(current_words[0]));
    current->krylov = proposal.krylov;
    current->guide = proposal.guide;
    ++current->accepted_generation;
  }
  candidate.accepted_generation = current->accepted_generation;
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_step_selected_neighbor(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *policy,
    const MVMCKrylovFockModel *proposal_model,
    uint64_t *current_words, size_t current_word_count,
    MVMCKrylovPositiveSamplerSnapshot *current,
    size_t selected_neighbor_index, double uniform,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words, size_t proposal_word_count,
    MVMCKrylovPositiveSamplerProposalStepResult *result) {
  MVMCKrylovStatus status;
  size_t neighbor_count = 0;
  double log_proposal_ratio = 0.0;
  MVMCKrylovPositiveSamplerStepResult step;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_proposal_step(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (current_words == NULL || current_word_count == 0 ||
      proposal_words == NULL || proposal_word_count != current_word_count) {
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
      proposal_model, current_words, current_word_count, proposal_words,
      proposal_word_count, &log_proposal_ratio);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_proposal_step(status, result);
    return status;
  }
  status = mvmc_krylov_positive_sampler_step(
      workspace, policy, current_words, current_word_count, current,
      proposal_words, proposal_word_count, log_proposal_ratio, uniform,
      amplitude, amplitude_context, &step);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_proposal_step(status, result);
    return status;
  }

  result->valid = 1;
  result->status = MVMC_KRYLOV_STATUS_OK;
  result->component = MVMC_KRYLOV_PROPOSAL_COMPONENT_NEIGHBOR;
  result->self_proposal = 0;
  result->configuration_changed = step.accepted;
  result->neighbor_count = neighbor_count;
  result->selected_neighbor_index = selected_neighbor_index;
  result->uniform = uniform;
  result->log_proposal_ratio = log_proposal_ratio;
  result->step = step;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_step_rng(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *policy,
    const MVMCKrylovFockModel *proposal_model,
    uint64_t *current_words, size_t current_word_count,
    MVMCKrylovPositiveSamplerSnapshot *current,
    MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words, size_t proposal_word_count,
    MVMCKrylovPositiveSamplerProposalStepResult *result) {
  MVMCKrylovStatus status;
  MVMCKrylovPositiveSamplerRng candidate_rng;
  size_t neighbor_count = 0;
  size_t selected_neighbor_index = 0;
  double uniform = NAN;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_proposal_step(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (!rng_valid(rng) || current_words == NULL || current_word_count == 0 ||
      proposal_words == NULL || proposal_word_count != current_word_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_fock_proposal_count_neighbors(
      proposal_model, current_words, current_word_count, &neighbor_count);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_proposal_step(status, result);
    return status;
  }
  if (neighbor_count == 0) {
    invalidate_proposal_step(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }

  candidate_rng = *rng;
  status = mvmc_krylov_positive_sampler_rng_draw_bounded(
      &candidate_rng, neighbor_count, &selected_neighbor_index);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_proposal_step(status, result);
    return status;
  }
  status = mvmc_krylov_positive_sampler_rng_draw_uniform(
      &candidate_rng, &uniform);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_proposal_step(status, result);
    return status;
  }
  status = mvmc_krylov_positive_sampler_step_selected_neighbor(
      workspace, policy, proposal_model, current_words, current_word_count,
      current, selected_neighbor_index, uniform, amplitude, amplitude_context,
      proposal_words, proposal_word_count, result);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_proposal_step(status, result);
    return status;
  }
  *rng = candidate_rng;
  result->rng_after = candidate_rng;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus sampler_rng_bounded_draw(
    void *context, size_t bound, size_t *value) {
  return mvmc_krylov_positive_sampler_rng_draw_bounded(
      (MVMCKrylovPositiveSamplerRng *)context, bound, value);
}

static int configurations_equal(const uint64_t *left,
                                const uint64_t *right,
                                size_t word_count) {
  return memcmp(left, right, word_count * sizeof(left[0])) == 0;
}

static void zero_proposal_words(uint64_t *proposal_words,
                                size_t proposal_word_count) {
  if (proposal_words != NULL &&
      proposal_word_count <= SIZE_MAX / sizeof(proposal_words[0])) {
    memset(proposal_words, 0,
           proposal_word_count * sizeof(proposal_words[0]));
  }
}

static int sampler_word_buffers_overlap(
    const uint64_t *first, size_t first_count,
    const uint64_t *second, size_t second_count) {
  uintptr_t first_start;
  uintptr_t second_start;
  size_t first_bytes;
  size_t second_bytes;
  if (first == NULL || second == NULL ||
      first_count > SIZE_MAX / sizeof(first[0]) ||
      second_count > SIZE_MAX / sizeof(second[0])) {
    return 1;
  }
  first_bytes = first_count * sizeof(first[0]);
  second_bytes = second_count * sizeof(second[0]);
  first_start = (uintptr_t)first;
  second_start = (uintptr_t)second;
  if (first_start > UINTPTR_MAX - first_bytes ||
      second_start > UINTPTR_MAX - second_bytes) {
    return 1;
  }
  return first_start < second_start + second_bytes &&
         second_start < first_start + first_bytes;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_draw_mixture_rng(
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const uint64_t *current_words, size_t current_word_count,
    const MVMCKrylovPositiveSamplerRng *rng,
    uint64_t *proposal_words, size_t proposal_word_count,
    MVMCKrylovPositiveSamplerProposalDrawResult *result) {
  MVMCKrylovStatus status;
  MVMCKrylovPositiveSamplerRng candidate_rng;
  MVMCKrylovFockUniformProposalResult global_result;
  MVMCKrylovFockShellProposalResult shell_result;
  size_t neighbor_count = 0;
  size_t selected_neighbor_index = SIZE_MAX;
  size_t component_draw = 0;
  uint64_t proposal_model_hash = 0;
  double log_proposal_ratio = 0.0;
  double uniform = NAN;
  const double proposal_start = sampler_wall_seconds();
  double component_start;
  double component_selection_seconds;
  double global_subset_seconds = 0.0;
  double shell_generation_seconds = 0.0;
  size_t shell_max_distance = 0;
  size_t shell_distance = 0;
  int use_global = 0;
  int use_shell = 0;
  int self_proposal = 0;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_proposal_draw(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (current_words == NULL || current_word_count == 0 ||
      proposal_words == NULL || proposal_word_count != current_word_count) {
    zero_proposal_words(proposal_words, proposal_word_count);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (sampler_word_buffers_overlap(current_words, current_word_count,
                                   proposal_words, proposal_word_count)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (!proposal_policy_valid(proposal_policy) || !rng_valid(rng) ||
      rng->algorithm != proposal_policy->rng_algorithm) {
    zero_proposal_words(proposal_words, proposal_word_count);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_fock_proposal_count_neighbors(
      proposal_model, current_words, current_word_count, &neighbor_count);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    zero_proposal_words(proposal_words, proposal_word_count);
    invalidate_proposal_draw(status, result);
    return status;
  }
  status = mvmc_krylov_positive_sampler_proposal_model_policy_hash(
      proposal_policy, proposal_model, current_words, current_word_count,
      &proposal_model_hash);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    zero_proposal_words(proposal_words, proposal_word_count);
    invalidate_proposal_draw(status, result);
    return status;
  }
  if (proposal_policy->kernel_id ==
      MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_SHELL_MIXTURE) {
    status = mvmc_krylov_fock_proposal_resolve_shell_distance(
        proposal_model, proposal_policy->distance_numerator,
        proposal_policy->distance_denominator, &shell_max_distance,
        &shell_distance);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      zero_proposal_words(proposal_words, proposal_word_count);
      invalidate_proposal_draw(status, result);
      return status;
    }
  }
  if (neighbor_count == 0 &&
      (proposal_policy->kernel_id ==
           MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_SHELL_MIXTURE ||
       proposal_policy->global_numerator <
           proposal_policy->global_denominator)) {
    zero_proposal_words(proposal_words, proposal_word_count);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }

  candidate_rng = *rng;
  component_start = sampler_wall_seconds();
  if (proposal_policy->kernel_id ==
      MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_GLOBAL_MIXTURE) {
    if (proposal_policy->global_numerator == 0) {
      use_global = 0;
    } else if (proposal_policy->global_numerator ==
               proposal_policy->global_denominator) {
      use_global = 1;
    } else {
      status = mvmc_krylov_positive_sampler_rng_draw_bounded(
          &candidate_rng, proposal_policy->global_denominator,
          &component_draw);
      if (status != MVMC_KRYLOV_STATUS_OK) goto fail;
      use_global = component_draw < proposal_policy->global_numerator;
    }
  } else {
    status = mvmc_krylov_positive_sampler_rng_draw_bounded(
        &candidate_rng, proposal_policy->neighbor_denominator,
        &component_draw);
    if (status != MVMC_KRYLOV_STATUS_OK) goto fail;
    use_shell = component_draw >= proposal_policy->neighbor_numerator;
  }
  component_selection_seconds = sampler_elapsed_seconds(component_start);

  memset(&global_result, 0, sizeof(global_result));
  memset(&shell_result, 0, sizeof(shell_result));
  if (use_global) {
    const double subset_start = sampler_wall_seconds();
    status = mvmc_krylov_fock_proposal_draw_uniform_sector(
        proposal_model, sampler_rng_bounded_draw, &candidate_rng,
        proposal_words, proposal_word_count, &global_result);
    if (status != MVMC_KRYLOV_STATUS_OK || !global_result.valid) {
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      }
      goto fail;
    }
    global_subset_seconds = sampler_elapsed_seconds(subset_start);
    self_proposal = configurations_equal(
        current_words, proposal_words, current_word_count);
  } else if (use_shell) {
    const double shell_start = sampler_wall_seconds();
    status = mvmc_krylov_fock_proposal_draw_shell(
        proposal_model, current_words, current_word_count, shell_distance,
        sampler_rng_bounded_draw, &candidate_rng, proposal_words,
        proposal_word_count, &shell_result);
    if (status != MVMC_KRYLOV_STATUS_OK || !shell_result.valid) {
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      }
      goto fail;
    }
    if (shell_result.max_distance != shell_max_distance ||
        shell_result.distance != shell_distance) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      goto fail;
    }
    shell_generation_seconds = sampler_elapsed_seconds(shell_start);
  } else {
    status = mvmc_krylov_positive_sampler_rng_draw_bounded(
        &candidate_rng, neighbor_count, &selected_neighbor_index);
    if (status != MVMC_KRYLOV_STATUS_OK) goto fail;
    status = mvmc_krylov_fock_proposal_select_neighbor(
        proposal_model, current_words, current_word_count,
        selected_neighbor_index, proposal_words, proposal_word_count,
        &neighbor_count);
    if (status != MVMC_KRYLOV_STATUS_OK) goto fail;
    status = mvmc_krylov_fock_proposal_log_ratio(
        proposal_model, current_words, current_word_count, proposal_words,
        proposal_word_count, &log_proposal_ratio);
    if (status != MVMC_KRYLOV_STATUS_OK) goto fail;
    if (log_proposal_ratio != 0.0) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      goto fail;
    }
  }
  status = mvmc_krylov_positive_sampler_rng_draw_uniform(
      &candidate_rng, &uniform);
  if (status != MVMC_KRYLOV_STATUS_OK) goto fail;

  result->valid = 1;
  result->status = MVMC_KRYLOV_STATUS_OK;
  result->component =
      use_global ? MVMC_KRYLOV_PROPOSAL_COMPONENT_GLOBAL
                 : (use_shell ? MVMC_KRYLOV_PROPOSAL_COMPONENT_SHELL
                              : MVMC_KRYLOV_PROPOSAL_COMPONENT_NEIGHBOR);
  result->self_proposal = self_proposal;
  result->neighbor_count = neighbor_count;
  result->selected_neighbor_index = selected_neighbor_index;
  result->component_draw_count =
      proposal_policy->kernel_id ==
              MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_SHELL_MIXTURE ||
              (proposal_policy->global_numerator != 0 &&
               proposal_policy->global_numerator !=
                   proposal_policy->global_denominator)
          ? 1 : 0;
  result->global_subset_draw_count =
      use_global ? global_result.draw_count : 0;
  result->shell_draw_count = use_shell ? shell_result.draw_count : 0;
  result->shell_max_distance = use_shell ? shell_result.max_distance : 0;
  result->shell_distance = use_shell ? shell_result.distance : 0;
  result->shell_up_distance = use_shell ? shell_result.up_distance : 0;
  result->shell_down_distance = use_shell ? shell_result.down_distance : 0;
  result->shell_count = use_shell ? shell_result.shell_count : 0;
  result->uniform = uniform;
  result->log_proposal_ratio = 0.0;
  result->proposal_seconds = sampler_elapsed_seconds(proposal_start);
  result->component_selection_seconds = component_selection_seconds;
  result->global_subset_seconds = global_subset_seconds;
  result->shell_generation_seconds = shell_generation_seconds;
  result->proposal_policy_hash = proposal_policy->policy_hash;
  result->proposal_model_hash = proposal_model_hash;
  result->rng_after = candidate_rng;
  return MVMC_KRYLOV_STATUS_OK;

fail:
  zero_proposal_words(proposal_words, proposal_word_count);
  invalidate_proposal_draw(status, result);
  return status;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_draw_surrogate_rng(
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy,
    const uint64_t *current_words, size_t current_word_count,
    double current_log_weight,
    const MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovSurrogateLogWeightCallback log_weight_callback,
    void *log_weight_context,
    uint64_t *proposal_words, size_t proposal_word_count,
    MVMCKrylovPositiveSamplerSurrogateDrawResult *result) {
  MVMCKrylovPositiveSamplerSurrogateDrawResult candidate_result;
  MVMCKrylovPositiveSamplerRng candidate_rng;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  uint64_t *scratch_allocation = NULL;
  uint64_t *candidate_allocation = NULL;
  uint64_t *scratch_words;
  uint64_t *candidate_words;
  uint64_t proposal_model_hash = 0;
  double scratch_log_weight = current_log_weight;
  const double start = sampler_wall_seconds();
  size_t step;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_surrogate_draw(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (!surrogate_policy_valid(surrogate_policy) || !rng_valid(rng) ||
      rng->algorithm != surrogate_policy->rng_algorithm ||
      current_words == NULL || current_word_count == 0 ||
      proposal_words == NULL || proposal_word_count != current_word_count ||
      sampler_word_buffers_overlap(current_words, current_word_count,
                                   proposal_words, proposal_word_count) ||
      current_word_count > SIZE_MAX / sizeof(current_words[0]) ||
      !isfinite(current_log_weight) || log_weight_callback == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_positive_sampler_proposal_model_hash(
      proposal_model, current_words, current_word_count,
      &proposal_model_hash);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_surrogate_draw(status, result);
    return status;
  }
  scratch_allocation =
      (uint64_t *)malloc(current_word_count * sizeof(current_words[0]));
  candidate_allocation =
      (uint64_t *)malloc(current_word_count * sizeof(current_words[0]));
  if (scratch_allocation == NULL || candidate_allocation == NULL) {
    free(candidate_allocation);
    free(scratch_allocation);
    invalidate_surrogate_draw(MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE, result);
    return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
  scratch_words = scratch_allocation;
  candidate_words = candidate_allocation;
  memcpy(scratch_words, current_words,
         current_word_count * sizeof(current_words[0]));
  candidate_rng = *rng;
  memset(&candidate_result, 0, sizeof(candidate_result));
  candidate_result.status = MVMC_KRYLOV_STATUS_OK;
  candidate_result.initial_log_weight = current_log_weight;

  for (step = 0; step < surrogate_policy->step_count; ++step) {
    size_t neighbor_count = 0;
    size_t selected_neighbor_index = 0;
    double neighbor_log_proposal_ratio = NAN;
    double proposal_log_weight = NAN;
    double uniform = NAN;
    int accepted;
    status = mvmc_krylov_fock_proposal_count_neighbors(
        proposal_model, scratch_words, current_word_count, &neighbor_count);
    if (status != MVMC_KRYLOV_STATUS_OK || neighbor_count == 0) {
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      }
      goto surrogate_fail;
    }
    status = mvmc_krylov_positive_sampler_rng_draw_bounded(
        &candidate_rng, neighbor_count, &selected_neighbor_index);
    if (status != MVMC_KRYLOV_STATUS_OK) goto surrogate_fail;
    status = mvmc_krylov_fock_proposal_select_neighbor(
        proposal_model, scratch_words, current_word_count,
        selected_neighbor_index, candidate_words, current_word_count,
        &neighbor_count);
    if (status != MVMC_KRYLOV_STATUS_OK) goto surrogate_fail;
    status = mvmc_krylov_fock_proposal_log_ratio(
        proposal_model, scratch_words, current_word_count, candidate_words,
        current_word_count, &neighbor_log_proposal_ratio);
    if (status != MVMC_KRYLOV_STATUS_OK ||
        neighbor_log_proposal_ratio != 0.0) {
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      }
      goto surrogate_fail;
    }
    status = log_weight_callback(candidate_words, current_word_count,
                                 log_weight_context,
                                 &proposal_log_weight);
    if (status != MVMC_KRYLOV_STATUS_OK ||
        !isfinite(proposal_log_weight)) {
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = MVMC_KRYLOV_STATUS_NONFINITE;
      }
      goto surrogate_fail;
    }
    status = mvmc_krylov_positive_sampler_rng_draw_uniform(
        &candidate_rng, &uniform);
    if (status != MVMC_KRYLOV_STATUS_OK) goto surrogate_fail;
    accepted = proposal_log_weight >= scratch_log_weight || uniform == 0.0 ||
               log(uniform) < proposal_log_weight - scratch_log_weight;
    if (candidate_result.surrogate_evaluation_count == SIZE_MAX ||
        (accepted &&
         (candidate_result.inner_accepted == SIZE_MAX ||
          candidate_result.inner_configuration_changes == SIZE_MAX)) ||
        (!accepted && candidate_result.inner_rejected == SIZE_MAX)) {
      status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      goto surrogate_fail;
    }
    if (accepted) {
      uint64_t *temporary = scratch_words;
      scratch_words = candidate_words;
      candidate_words = temporary;
      scratch_log_weight = proposal_log_weight;
      ++candidate_result.inner_accepted;
      ++candidate_result.inner_configuration_changes;
    } else {
      ++candidate_result.inner_rejected;
    }
    ++candidate_result.surrogate_evaluation_count;
  }
  candidate_result.log_proposal_ratio =
      current_log_weight - scratch_log_weight;
  if (!isfinite(candidate_result.log_proposal_ratio)) {
    status = MVMC_KRYLOV_STATUS_NONFINITE;
    goto surrogate_fail;
  }
  memcpy(proposal_words, scratch_words,
         current_word_count * sizeof(current_words[0]));
  candidate_result.valid = 1;
  candidate_result.step_count = surrogate_policy->step_count;
  candidate_result.final_configuration_changed =
      !configurations_equal(current_words, scratch_words,
                            current_word_count);
  candidate_result.final_log_weight = scratch_log_weight;
  candidate_result.proposal_seconds = sampler_elapsed_seconds(start);
  candidate_result.surrogate_policy_hash = surrogate_policy->policy_hash;
  candidate_result.proposal_model_hash = proposal_model_hash;
  candidate_result.final_configuration_hash =
      hash_configuration_words(scratch_words, current_word_count);
  candidate_result.rng_after = candidate_rng;
  free(candidate_allocation);
  free(scratch_allocation);
  *result = candidate_result;
  return MVMC_KRYLOV_STATUS_OK;

surrogate_fail:
  free(candidate_allocation);
  free(scratch_allocation);
  invalidate_surrogate_draw(status, result);
  return status;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_step_mixture_rng(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *policy,
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    uint64_t *current_words, size_t current_word_count,
    MVMCKrylovPositiveSamplerSnapshot *current,
    MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words, size_t proposal_word_count,
    MVMCKrylovPositiveSamplerProposalStepResult *result) {
  MVMCKrylovStatus status;
  MVMCKrylovPositiveSamplerProposalDrawResult draw;
  MVMCKrylovPositiveSamplerStepResult step;
  const double total_start = sampler_wall_seconds();

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_proposal_step(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  status = mvmc_krylov_positive_sampler_draw_mixture_rng(
      proposal_model, proposal_policy, current_words, current_word_count,
      rng, proposal_words, proposal_word_count, &draw);
  if (status != MVMC_KRYLOV_STATUS_OK || !draw.valid) {
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    invalidate_proposal_step(status, result);
    return status;
  }
  status = mvmc_krylov_positive_sampler_step(
      workspace, policy, current_words, current_word_count, current,
      proposal_words, proposal_word_count, draw.log_proposal_ratio,
      draw.uniform, amplitude, amplitude_context, &step);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    zero_proposal_words(proposal_words, proposal_word_count);
    invalidate_proposal_step(status, result);
    return status;
  }

  result->valid = 1;
  result->status = MVMC_KRYLOV_STATUS_OK;
  result->component = draw.component;
  result->self_proposal = draw.self_proposal;
  result->configuration_changed = step.accepted && !draw.self_proposal;
  result->neighbor_count = draw.neighbor_count;
  result->selected_neighbor_index = draw.selected_neighbor_index;
  result->component_draw_count = draw.component_draw_count;
  result->global_subset_draw_count = draw.global_subset_draw_count;
  result->shell_draw_count = draw.shell_draw_count;
  result->shell_max_distance = draw.shell_max_distance;
  result->shell_distance = draw.shell_distance;
  result->shell_up_distance = draw.shell_up_distance;
  result->shell_down_distance = draw.shell_down_distance;
  result->shell_count = draw.shell_count;
  result->uniform = draw.uniform;
  result->log_proposal_ratio = draw.log_proposal_ratio;
  result->proposal_seconds = draw.proposal_seconds;
  result->component_selection_seconds =
      draw.component_selection_seconds;
  result->global_subset_seconds = draw.global_subset_seconds;
  result->shell_generation_seconds = draw.shell_generation_seconds;
  result->bounded_evaluation_seconds =
      step.proposal_krylov.statistics.total_seconds;
  result->total_step_seconds = sampler_elapsed_seconds(total_start);
  result->proposal_policy_hash = draw.proposal_policy_hash;
  result->proposal_model_hash = draw.proposal_model_hash;
  result->rng_after = draw.rng_after;
  result->step = step;
  *rng = draw.rng_after;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_manifest_create(
    const MVMCKrylovPositiveGuidePolicy *policy,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovBoundedLimits *limits,
    uint64_t bounded_plan_hash,
    const uint64_t *current_words, size_t current_word_count,
    const MVMCKrylovPositiveSamplerSnapshot *current,
    const MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovPositiveSamplerManifest *manifest) {
  MVMCKrylovPositiveSamplerManifest candidate;
  MVMCKrylovStatus model_status;
  uint64_t proposal_model_hash = 0;

  if (manifest == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_manifest(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, manifest);
  if (sizeof(double) != sizeof(uint64_t)) {
    invalidate_manifest(MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE,
                        manifest);
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (!policy_value_valid(policy) ||
      !proposal_policy_valid(proposal_policy) ||
      !bounded_limits_valid(limits) ||
      current_words == NULL || current_word_count == 0 ||
      current == NULL || !current->valid ||
      current->status != MVMC_KRYLOV_STATUS_OK ||
      current->word_count != current_word_count ||
      current->policy_hash != policy->policy_hash ||
      !current->krylov.valid ||
      current->krylov.status != MVMC_KRYLOV_STATUS_OK ||
      current->krylov.evaluated_order < policy->order ||
      !current->guide.valid ||
      current->guide.status != MVMC_KRYLOV_STATUS_OK ||
      !rng_valid(rng)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  model_status = mvmc_krylov_positive_sampler_proposal_model_policy_hash(
      proposal_policy, proposal_model, current_words, current_word_count,
      &proposal_model_hash);
  if (model_status != MVMC_KRYLOV_STATUS_OK) return model_status;

  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.state_version = MVMC_KRYLOV_POSITIVE_SAMPLER_STATE_VERSION;
  candidate.rng_algorithm = rng->algorithm;
  candidate.rng_state = rng->state;
  candidate.rng_stream = rng->stream;
  candidate.rng_draws = rng->draws;
  candidate.policy_hash = policy->policy_hash;
  candidate.proposal_policy_hash = proposal_policy->policy_hash;
  candidate.proposal_model_hash = proposal_model_hash;
  candidate.guide_shape_hash = hash_guide_policy(policy);
  candidate.bounded_plan_hash = bounded_plan_hash;
  candidate.amplitude_policy_hash = limits->amplitude_policy_hash;
  candidate.current_configuration_hash =
      hash_configuration_words(current_words, current_word_count);
  candidate.current_scale_hash = hash_scaled_values(
      current->krylov.value, current->krylov.evaluated_order);
  candidate.accepted_generation = current->accepted_generation;
  candidate.word_count = current_word_count;
  candidate.cache_bytes = limits->cache_bytes;
  candidate.guide_order = policy->order;
  candidate.bounded_max_order = limits->max_order;
  *manifest = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_manifest_matches(
    const MVMCKrylovPositiveSamplerManifest *expected,
    const MVMCKrylovPositiveGuidePolicy *policy,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovBoundedLimits *limits,
    uint64_t bounded_plan_hash,
    const uint64_t *current_words, size_t current_word_count,
    const MVMCKrylovPositiveSamplerSnapshot *current,
    const MVMCKrylovPositiveSamplerRng *rng,
    int *matches) {
  MVMCKrylovPositiveSamplerManifest candidate;
  MVMCKrylovStatus status;
  if (matches == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *matches = 0;
  if (expected == NULL || !expected->valid ||
      expected->status != MVMC_KRYLOV_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_positive_sampler_manifest_create(
      policy, proposal_policy, proposal_model, limits, bounded_plan_hash,
      current_words, current_word_count, current, rng, &candidate);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  *matches = manifest_equal(expected, &candidate);
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_surrogate_manifest_create(
    const MVMCKrylovPositiveGuidePolicy *policy,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy,
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovBoundedLimits *limits,
    uint64_t bounded_plan_hash,
    const uint64_t *current_words, size_t current_word_count,
    const MVMCKrylovPositiveSamplerSnapshot *current,
    const MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovPositiveSamplerManifest *manifest) {
  MVMCKrylovStatus status;
  if (manifest == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_manifest(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, manifest);
  if (!surrogate_policy_valid(surrogate_policy) ||
      !policy_value_valid(policy) ||
      surrogate_policy->guide_shape_hash != hash_guide_policy(policy)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_positive_sampler_manifest_create(
      policy, proposal_policy, proposal_model, limits, bounded_plan_hash,
      current_words, current_word_count, current, rng, manifest);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  manifest->surrogate_policy_hash = surrogate_policy->policy_hash;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_surrogate_manifest_matches(
    const MVMCKrylovPositiveSamplerManifest *expected,
    const MVMCKrylovPositiveGuidePolicy *policy,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy,
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovBoundedLimits *limits,
    uint64_t bounded_plan_hash,
    const uint64_t *current_words, size_t current_word_count,
    const MVMCKrylovPositiveSamplerSnapshot *current,
    const MVMCKrylovPositiveSamplerRng *rng,
    int *matches) {
  MVMCKrylovPositiveSamplerManifest candidate;
  MVMCKrylovStatus status;
  if (matches == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *matches = 0;
  if (expected == NULL || !expected->valid ||
      expected->status != MVMC_KRYLOV_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_positive_sampler_surrogate_manifest_create(
      policy, proposal_policy, surrogate_policy, proposal_model, limits,
      bounded_plan_hash, current_words, current_word_count, current, rng,
      &candidate);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  *matches = manifest_equal(expected, &candidate);
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_trace_statistics_reset(
    MVMCKrylovPositiveSamplerTraceStatistics *statistics) {
  if (statistics == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(statistics, 0, sizeof(*statistics));
  statistics->valid = 1;
  statistics->status = MVMC_KRYLOV_STATUS_OK;
  statistics->trace_hash = UINT64_C(1469598103934665603);
  statistics->min_log_acceptance_ratio = NAN;
  statistics->max_log_acceptance_ratio = NAN;
  statistics->min_proposal_log_guide = NAN;
  statistics->max_proposal_log_guide = NAN;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_sampler_trace_statistics_record_step(
    const MVMCKrylovPositiveSamplerProposalStepResult *step,
    const uint64_t *current_words, size_t current_word_count,
    MVMCKrylovPositiveSamplerTraceStatistics *statistics) {
  MVMCKrylovPositiveSamplerTraceStatistics candidate;
  const MVMCKrylovPositiveGuideResult *guide;
  const MVMCKrylovPositiveGuideAcceptance *acceptance;
  const uint64_t exact_zero_count =
      step != NULL
          ? (uint64_t)step->step.proposal_guide.exact_zero_component_count
          : 0;
  const uint64_t numeric_zero_count =
      step != NULL
          ? (uint64_t)step->step.proposal_guide.numeric_zero_component_count
          : 0;
  const uint64_t finite_count =
      step != NULL
          ? (uint64_t)step->step.proposal_guide.finite_component_count
          : 0;
  const uint64_t zero_count = exact_zero_count + numeric_zero_count;

  if (statistics == NULL || !statistics->valid ||
      statistics->status != MVMC_KRYLOV_STATUS_OK || step == NULL ||
      !step->valid || step->status != MVMC_KRYLOV_STATUS_OK ||
      !step->step.valid || step->step.status != MVMC_KRYLOV_STATUS_OK ||
      current_words == NULL || current_word_count == 0 ||
      step->step.proposal_guide.exact_zero_component_count < 0 ||
      step->step.proposal_guide.numeric_zero_component_count < 0 ||
      step->step.proposal_guide.finite_component_count < 0 ||
      zero_count < exact_zero_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }

  guide = &step->step.proposal_guide;
  acceptance = &step->step.acceptance;
  if (!guide->valid || guide->status != MVMC_KRYLOV_STATUS_OK ||
      !acceptance->valid || acceptance->status != MVMC_KRYLOV_STATUS_OK ||
      !isfinite(guide->log_guide) ||
      !isfinite(acceptance->log_acceptance_ratio)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }

  candidate = *statistics;
  if (!increment_u64(&candidate.attempted_steps) ||
      !increment_u64(&candidate.completed_steps)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  if (step->step.accepted) {
    if (!increment_u64(&candidate.accepted_steps)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  } else if (!increment_u64(&candidate.rejected_steps)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  if (step->component == MVMC_KRYLOV_PROPOSAL_COMPONENT_NEIGHBOR) {
    if (!increment_u64(&candidate.neighbor_attempted_steps) ||
        (step->step.accepted
             ? !increment_u64(&candidate.neighbor_accepted_steps)
             : !increment_u64(&candidate.neighbor_rejected_steps))) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  } else if (step->component == MVMC_KRYLOV_PROPOSAL_COMPONENT_GLOBAL) {
    if (!increment_u64(&candidate.global_attempted_steps) ||
        (step->step.accepted
             ? !increment_u64(&candidate.global_accepted_steps)
             : !increment_u64(&candidate.global_rejected_steps)) ||
        (step->self_proposal &&
         !increment_u64(&candidate.global_self_proposals))) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  } else if (step->component == MVMC_KRYLOV_PROPOSAL_COMPONENT_SHELL) {
    if (step->self_proposal || step->shell_distance == 0 ||
        step->shell_distance > step->shell_max_distance ||
        step->shell_count == 0) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    if (!increment_u64(&candidate.shell_attempted_steps) ||
        (step->step.accepted
             ? !increment_u64(&candidate.shell_accepted_steps)
             : !increment_u64(&candidate.shell_rejected_steps))) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  } else {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (step->configuration_changed &&
      !increment_u64(&candidate.configuration_changing_accepted_moves)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  if (isfinite(guide->log_floor) && guide->log_floor <= guide->log_guide) {
    if (!increment_u64(&candidate.positive_support_steps)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    if (zero_count != 0 &&
        !increment_u64(&candidate.floor_supported_zero_steps)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  } else if (!increment_u64(&candidate.support_violation_steps)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  if (!add_u64(&candidate.finite_proposal_components, finite_count) ||
      !add_u64(&candidate.exact_zero_proposal_components,
               exact_zero_count) ||
      !add_u64(&candidate.numeric_zero_proposal_components,
               numeric_zero_count) ||
      !add_u64(&candidate.terminal_amplitude_calls,
               step->step.proposal_krylov.statistics
                   .terminal_amplitude_calls)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }

  if (candidate.completed_steps == 1) {
    candidate.min_log_acceptance_ratio =
        acceptance->log_acceptance_ratio;
    candidate.max_log_acceptance_ratio =
        acceptance->log_acceptance_ratio;
    candidate.min_proposal_log_guide = guide->log_guide;
    candidate.max_proposal_log_guide = guide->log_guide;
  } else {
    candidate.min_log_acceptance_ratio =
        fmin(candidate.min_log_acceptance_ratio,
             acceptance->log_acceptance_ratio);
    candidate.max_log_acceptance_ratio =
        fmax(candidate.max_log_acceptance_ratio,
             acceptance->log_acceptance_ratio);
    candidate.min_proposal_log_guide =
        fmin(candidate.min_proposal_log_guide, guide->log_guide);
    candidate.max_proposal_log_guide =
        fmax(candidate.max_proposal_log_guide, guide->log_guide);
  }
  candidate.sum_log_acceptance_ratio +=
      acceptance->log_acceptance_ratio;
  if (!isfinite(candidate.sum_log_acceptance_ratio)) {
    invalidate_trace_statistics(MVMC_KRYLOV_STATUS_NONFINITE,
                                statistics);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }

  hash_u64(&candidate.trace_hash, candidate.completed_steps);
  hash_u64(&candidate.trace_hash, (uint64_t)(unsigned int)step->component);
  hash_u64(&candidate.trace_hash, (uint64_t)(unsigned int)step->self_proposal);
  hash_u64(&candidate.trace_hash,
           (uint64_t)(unsigned int)step->configuration_changed);
  hash_u64(&candidate.trace_hash, (uint64_t)step->selected_neighbor_index);
  hash_u64(&candidate.trace_hash, (uint64_t)step->component_draw_count);
  hash_u64(&candidate.trace_hash, (uint64_t)step->global_subset_draw_count);
  hash_u64(&candidate.trace_hash, (uint64_t)step->shell_draw_count);
  hash_u64(&candidate.trace_hash, (uint64_t)step->shell_max_distance);
  hash_u64(&candidate.trace_hash, (uint64_t)step->shell_distance);
  hash_u64(&candidate.trace_hash, (uint64_t)step->shell_up_distance);
  hash_u64(&candidate.trace_hash, (uint64_t)step->shell_down_distance);
  hash_u64(&candidate.trace_hash, (uint64_t)step->shell_count);
  hash_u64(&candidate.trace_hash, step->proposal_policy_hash);
  hash_u64(&candidate.trace_hash, step->proposal_model_hash);
  hash_u64(&candidate.trace_hash, double_bits(step->uniform));
  hash_u64(&candidate.trace_hash,
           double_bits(acceptance->log_acceptance_ratio));
  hash_u64(&candidate.trace_hash, (uint64_t)step->step.accepted);
  hash_u64(&candidate.trace_hash,
           hash_configuration_words(current_words, current_word_count));
  hash_u64(&candidate.trace_hash, finite_count);
  hash_u64(&candidate.trace_hash, exact_zero_count);
  hash_u64(&candidate.trace_hash, numeric_zero_count);
  if (step->rng_after.valid) {
    hash_u64(&candidate.trace_hash, step->rng_after.draws);
    hash_u64(&candidate.trace_hash, step->rng_after.state);
  }

  *statistics = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}
