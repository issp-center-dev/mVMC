/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "krylov_final_state_chain.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_P5_TESTING) ||                          \
    !defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                         \
    !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "krylov_final_state_chain.c is Testing-only"
#endif

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

static void hash_u64(uint64_t *hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    *hash ^= value & UINT64_C(0xff);
    *hash *= UINT64_C(1099511628211);
    value >>= 8;
  }
}

static uint64_t configuration_hash(const uint64_t *words,
                                   size_t word_count) {
  uint64_t hash = UINT64_C(1469598103934665603);
  size_t word;
  hash_u64(&hash, (uint64_t)word_count);
  for (word = 0; word < word_count; ++word) hash_u64(&hash, words[word]);
  return hash == 0 ? UINT64_C(1) : hash;
}

static uint64_t component_state_hash(
    const MVMCKrylovFinalStateSnapshot *current) {
  return mvmc_krylov_final_state_snapshot_integrity_hash(current);
}

static int rng_valid(const MVMCKrylovPositiveSamplerRng *rng) {
  return rng != NULL && rng->valid &&
         rng->state_version == MVMC_KRYLOV_POSITIVE_SAMPLER_STATE_VERSION &&
         rng->algorithm == MVMC_KRYLOV_POSITIVE_SAMPLER_RNG_SPLITMIX64;
}

static int limits_valid(const MVMCKrylovBoundedLimits *limits) {
  return limits != NULL && limits->amplitude_policy_hash != 0 &&
         limits->max_order >= 0 &&
         limits->max_order <= MVMC_KRYLOV_MAX_ORDER &&
         limits->cache_bytes != 0 && limits->max_row_transitions != 0 &&
         limits->max_workspace_bytes != 0 &&
         limits->max_node_expansions != 0 &&
         limits->max_terminal_amplitude_calls != 0 &&
         limits->max_total_row_transitions != 0;
}

static void invalidate_chain_step(
    MVMCKrylovStatus status, MVMCKrylovFinalStateChainStepResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
}

static void invalidate_manifest(
    MVMCKrylovStatus status, MVMCKrylovFinalStateChainManifest *manifest) {
  if (manifest == NULL) return;
  memset(manifest, 0, sizeof(*manifest));
  manifest->status = status;
}

static int manifest_equal(const MVMCKrylovFinalStateChainManifest *left,
                          const MVMCKrylovFinalStateChainManifest *right) {
  return left->valid == right->valid && left->status == right->status &&
         left->state_version == right->state_version &&
         left->rng_algorithm == right->rng_algorithm &&
         left->rng_state_version == right->rng_state_version &&
         left->rng_state == right->rng_state &&
         left->rng_stream == right->rng_stream &&
         left->rng_draws == right->rng_draws &&
         left->final_policy_hash == right->final_policy_hash &&
         left->coefficient_source == right->coefficient_source &&
         left->coefficient_provenance_hash ==
             right->coefficient_provenance_hash &&
         left->proposal_policy_hash == right->proposal_policy_hash &&
         left->proposal_model_hash == right->proposal_model_hash &&
         left->bounded_plan_hash == right->bounded_plan_hash &&
         left->amplitude_policy_hash == right->amplitude_policy_hash &&
         left->amplitude_generation_hash ==
             right->amplitude_generation_hash &&
         left->current_configuration_hash ==
             right->current_configuration_hash &&
         left->component_state_hash == right->component_state_hash &&
         left->accepted_generation == right->accepted_generation &&
         left->word_count == right->word_count &&
         left->final_order == right->final_order &&
         left->bounded_max_order == right->bounded_max_order;
}

MVMCKrylovStatus mvmc_krylov_final_state_observable_measure(
    MVMCKrylovFinalObservable observable,
    const MVMCKrylovFinalStateEvaluation *evaluation,
    double complex *value) {
  double complex local_energy = NAN + I * NAN;
  MVMCScaledComplexExportStatus export_status;
  if (value == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *value = NAN + I * NAN;
  if (evaluation == NULL || !evaluation->valid ||
      evaluation->status != MVMC_KRYLOV_STATUS_OK ||
      !evaluation->sampleable || !evaluation->local_energy_available) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (observable != MVMC_KRYLOV_FINAL_OBSERVABLE_ENERGY_COMPLEX &&
      observable !=
          MVMC_KRYLOV_FINAL_OBSERVABLE_LOCAL_ENERGY_ABS_SQUARED_DIAGNOSTIC) {
    return MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL;
  }
  export_status = mvmc_krylov_final_state_local_energy_export(
      evaluation, &local_energy);
  if (export_status != MVMC_SCALED_EXPORT_OK &&
      export_status != MVMC_SCALED_EXPORT_EXACT_ZERO) {
    return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  if (observable == MVMC_KRYLOV_FINAL_OBSERVABLE_ENERGY_COMPLEX) {
    *value = local_energy;
  } else {
    *value = creal(local_energy * conj(local_energy));
  }
  if (!isfinite(creal(*value)) || !isfinite(cimag(*value))) {
    *value = NAN + I * NAN;
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_final_state_chain_step_mixture_rng(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy,
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    uint64_t *current_words, size_t current_word_count,
    MVMCKrylovFinalStateSnapshot *current,
    MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words, size_t proposal_word_count,
    MVMCKrylovFinalStateChainStepResult *result) {
  MVMCKrylovPositiveSamplerProposalDrawResult proposal;
  MVMCKrylovFinalStateStepResult step;
  MVMCKrylovStatus status;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_chain_step(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (!rng_valid(rng)) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  status = mvmc_krylov_positive_sampler_draw_mixture_rng(
      proposal_model, proposal_policy, current_words, current_word_count,
      rng, proposal_words, proposal_word_count, &proposal);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_chain_step(status, result);
    return status;
  }
  status = mvmc_krylov_final_state_sampler_step(
      workspace, policy, current_words, current_word_count, current,
      proposal_words, proposal_word_count, proposal.log_proposal_ratio,
      proposal.uniform, amplitude, amplitude_context, &step);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    if (workspace != NULL && current != NULL && current->valid &&
        current->amplitude_generation_hash != 0 &&
        !mvmc_bounded_krylov_session_is_active(workspace)) {
      (void)mvmc_bounded_krylov_session_begin(
          workspace, amplitude, amplitude_context,
          current->amplitude_generation_hash);
    }
    invalidate_chain_step(status, result);
    return status;
  }
  result->valid = 1;
  result->status = MVMC_KRYLOV_STATUS_OK;
  result->configuration_changed = step.accepted && !proposal.self_proposal;
  result->proposal = proposal;
  result->step = step;
  result->rng_after = proposal.rng_after;
  *rng = proposal.rng_after;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_final_state_chain_manifest_create(
    const MVMCKrylovFinalStatePolicy *policy,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovBoundedLimits *limits, uint64_t bounded_plan_hash,
    const uint64_t *current_words, size_t current_word_count,
    const MVMCKrylovFinalStateSnapshot *current,
    const MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovFinalStateChainManifest *manifest) {
  MVMCKrylovFinalStateChainManifest candidate;
  MVMCKrylovFinalStateEvaluation reevaluated;
  MVMCKrylovFinalStateSnapshot canonical_snapshot;
  uint64_t proposal_model_hash = 0;
  MVMCKrylovStatus status;
  if (manifest == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_manifest(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, manifest);
  if (mvmc_krylov_final_state_policy_hash(policy) == 0 ||
      proposal_policy == NULL || !proposal_policy->valid ||
      proposal_policy->status != MVMC_KRYLOV_STATUS_OK ||
      !limits_valid(limits) || bounded_plan_hash == 0 ||
      current_words == NULL || current_word_count == 0 || current == NULL ||
      !current->valid || current->status != MVMC_KRYLOV_STATUS_OK ||
      current->word_count != current_word_count ||
      current->policy_hash != policy->policy_hash ||
      current->configuration_hash !=
          configuration_hash(current_words, current_word_count) ||
      current->integrity_hash == 0 ||
      current->integrity_hash != component_state_hash(current) ||
      current->amplitude_generation_hash == 0 ||
      !current->final_state.sampleable ||
      !current->final_state.local_energy_available || !rng_valid(rng) ||
      limits->max_order < policy->order + 1) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_final_state_evaluate(
      policy, &current->krylov, &reevaluated);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  canonical_snapshot = *current;
  canonical_snapshot.final_state = reevaluated;
  if (component_state_hash(current) !=
      component_state_hash(&canonical_snapshot)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_positive_sampler_proposal_model_policy_hash(
      proposal_policy, proposal_model, current_words, current_word_count,
      &proposal_model_hash);
  if (status != MVMC_KRYLOV_STATUS_OK || proposal_model_hash == 0) {
    invalidate_manifest(status, manifest);
    return status;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.state_version = MVMC_KRYLOV_FINAL_CHAIN_STATE_VERSION;
  candidate.rng_algorithm = rng->algorithm;
  candidate.rng_state_version = rng->state_version;
  candidate.rng_state = rng->state;
  candidate.rng_stream = rng->stream;
  candidate.rng_draws = rng->draws;
  candidate.final_policy_hash = policy->policy_hash;
  candidate.coefficient_source = (uint64_t)(unsigned int)policy->source;
  candidate.coefficient_provenance_hash = policy->provenance_hash;
  candidate.proposal_policy_hash = proposal_policy->policy_hash;
  candidate.proposal_model_hash = proposal_model_hash;
  candidate.bounded_plan_hash = bounded_plan_hash;
  candidate.amplitude_policy_hash = limits->amplitude_policy_hash;
  candidate.amplitude_generation_hash = current->amplitude_generation_hash;
  candidate.current_configuration_hash = current->configuration_hash;
  candidate.component_state_hash = component_state_hash(current);
  candidate.accepted_generation = current->accepted_generation;
  candidate.word_count = current_word_count;
  candidate.final_order = policy->order;
  candidate.bounded_max_order = limits->max_order;
  *manifest = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_final_state_chain_manifest_matches(
    const MVMCKrylovFinalStateChainManifest *expected,
    const MVMCKrylovFinalStatePolicy *policy,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovBoundedLimits *limits, uint64_t bounded_plan_hash,
    const uint64_t *current_words, size_t current_word_count,
    const MVMCKrylovFinalStateSnapshot *current,
    const MVMCKrylovPositiveSamplerRng *rng, int *matches) {
  MVMCKrylovFinalStateChainManifest candidate;
  MVMCKrylovStatus status;
  if (matches == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *matches = 0;
  if (expected == NULL || !expected->valid ||
      expected->status != MVMC_KRYLOV_STATUS_OK ||
      expected->state_version != MVMC_KRYLOV_FINAL_CHAIN_STATE_VERSION) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_final_state_chain_manifest_create(
      policy, proposal_policy, proposal_model, limits, bounded_plan_hash,
      current_words, current_word_count, current, rng, &candidate);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  *matches = manifest_equal(expected, &candidate);
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_final_state_chain_restart_restore(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStateChainManifest *expected,
    const MVMCKrylovFinalStatePolicy *policy,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovBoundedLimits *limits, uint64_t bounded_plan_hash,
    const uint64_t *current_words, size_t current_word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovFinalStateSnapshot *current,
    MVMCKrylovPositiveSamplerRng *rng) {
  MVMCKrylovFinalStateSnapshot restored;
  MVMCKrylovPositiveSamplerRng restored_rng;
  MVMCKrylovFinalStateChainManifest candidate;
  MVMCKrylovStatus status;
  if (current == NULL || rng == NULL || expected == NULL ||
      !expected->valid ||
      expected->status != MVMC_KRYLOV_STATUS_OK ||
      expected->state_version != MVMC_KRYLOV_FINAL_CHAIN_STATE_VERSION ||
      expected->amplitude_generation_hash == 0 ||
      expected->word_count != current_word_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&restored, 0, sizeof(restored));
  memset(&restored_rng, 0, sizeof(restored_rng));
  restored_rng.valid = 1;
  restored_rng.state_version = expected->rng_state_version;
  restored_rng.algorithm = expected->rng_algorithm;
  restored_rng.state = expected->rng_state;
  restored_rng.stream = expected->rng_stream;
  restored_rng.draws = expected->rng_draws;
  if (!rng_valid(&restored_rng)) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (!mvmc_bounded_krylov_session_is_active(workspace)) {
    status = mvmc_bounded_krylov_session_begin(
        workspace, amplitude, amplitude_context,
        expected->amplitude_generation_hash);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
  }
  status = mvmc_krylov_final_state_sampler_initialize(
      workspace, policy, current_words, current_word_count, amplitude,
      amplitude_context, &restored);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  restored.accepted_generation = expected->accepted_generation;
  status = mvmc_krylov_final_state_chain_manifest_create(
      policy, proposal_policy, proposal_model, limits, bounded_plan_hash,
      current_words, current_word_count, &restored, &restored_rng,
      &candidate);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (!manifest_equal(expected, &candidate)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *current = restored;
  *rng = restored_rng;
  return MVMC_KRYLOV_STATUS_OK;
}
