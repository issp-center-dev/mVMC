/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_KRYLOV_FINAL_STATE_CHAIN_H
#define MVMC_KRYLOV_FINAL_STATE_CHAIN_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_P5_TESTING)

#include "krylov_final_state.h"
#include "krylov_positive_sampler.h"

#define MVMC_KRYLOV_FINAL_CHAIN_STATE_VERSION UINT64_C(1)

typedef enum {
  MVMC_KRYLOV_FINAL_OBSERVABLE_ENERGY_COMPLEX = 1,
  MVMC_KRYLOV_FINAL_OBSERVABLE_LOCAL_ENERGY_ABS_SQUARED_DIAGNOSTIC = 2,
  MVMC_KRYLOV_FINAL_OBSERVABLE_ARBITRARY_UNSUPPORTED = 3
} MVMCKrylovFinalObservable;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int configuration_changed;
  MVMCKrylovPositiveSamplerProposalDrawResult proposal;
  MVMCKrylovFinalStateStepResult step;
  MVMCKrylovPositiveSamplerRng rng_after;
} MVMCKrylovFinalStateChainStepResult;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  uint64_t state_version;
  uint64_t rng_algorithm;
  uint64_t rng_state_version;
  uint64_t rng_state;
  uint64_t rng_stream;
  uint64_t rng_draws;
  uint64_t final_policy_hash;
  uint64_t coefficient_source;
  uint64_t coefficient_provenance_hash;
  uint64_t proposal_policy_hash;
  uint64_t proposal_model_hash;
  uint64_t bounded_plan_hash;
  uint64_t amplitude_policy_hash;
  uint64_t amplitude_generation_hash;
  uint64_t current_configuration_hash;
  uint64_t component_state_hash;
  uint64_t accepted_generation;
  size_t word_count;
  int final_order;
  int bounded_max_order;
} MVMCKrylovFinalStateChainManifest;

MVMCKrylovStatus mvmc_krylov_final_state_observable_measure(
    MVMCKrylovFinalObservable observable,
    const MVMCKrylovFinalStateEvaluation *evaluation,
    double complex *value);

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
    MVMCKrylovFinalStateChainStepResult *result);

MVMCKrylovStatus mvmc_krylov_final_state_chain_manifest_create(
    const MVMCKrylovFinalStatePolicy *policy,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovBoundedLimits *limits, uint64_t bounded_plan_hash,
    const uint64_t *current_words, size_t current_word_count,
    const MVMCKrylovFinalStateSnapshot *current,
    const MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovFinalStateChainManifest *manifest);

MVMCKrylovStatus mvmc_krylov_final_state_chain_manifest_matches(
    const MVMCKrylovFinalStateChainManifest *expected,
    const MVMCKrylovFinalStatePolicy *policy,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovBoundedLimits *limits, uint64_t bounded_plan_hash,
    const uint64_t *current_words, size_t current_word_count,
    const MVMCKrylovFinalStateSnapshot *current,
    const MVMCKrylovPositiveSamplerRng *rng, int *matches);

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
    MVMCKrylovPositiveSamplerRng *rng);

#endif /* MVMC_ENABLE_POWER_LANCZOS_P5_TESTING */

#endif /* MVMC_KRYLOV_FINAL_STATE_CHAIN_H */
