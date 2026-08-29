/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_KRYLOV_POSITIVE_SAMPLER_H
#define MVMC_KRYLOV_POSITIVE_SAMPLER_H

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) &&                         \
    defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "krylov_fock_proposal.h"
#include "krylov_positive_guide.h"

#define MVMC_KRYLOV_POSITIVE_SAMPLER_STATE_VERSION UINT64_C(3)
#define MVMC_KRYLOV_POSITIVE_SAMPLER_RNG_SPLITMIX64 UINT64_C(0x53504c49544d4958)
#define MVMC_KRYLOV_PROPOSAL_POLICY_VERSION UINT64_C(2)
#define MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_GLOBAL_MIXTURE                   \
  UINT64_C(0x4e474d4958545552)
#define MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_SHELL_MIXTURE                    \
  UINT64_C(0x4e534d4958545552)
#define MVMC_KRYLOV_PROPOSAL_SUBSET_CANONICAL_V1 UINT64_C(1)
#define MVMC_KRYLOV_PROPOSAL_DISTANCE_ROUND_HALF_UP_V1 UINT64_C(1)
#define MVMC_KRYLOV_PROPOSAL_SELF_REPORT_V1 UINT64_C(1)
#define MVMC_KRYLOV_PROPOSAL_MODEL_CONTRACT_VERSION UINT64_C(2)
#define MVMC_KRYLOV_SURROGATE_POLICY_VERSION UINT64_C(1)
#define MVMC_KRYLOV_SURROGATE_Q_ZEROTH_FLOOR_V1 UINT64_C(1)
#define MVMC_KRYLOV_SURROGATE_Q_PARTIAL_FLOOR_V1 UINT64_C(2)
#define MVMC_KRYLOV_SURROGATE_TRANSACTION_V1 UINT64_C(1)

typedef enum {
  MVMC_KRYLOV_PROPOSAL_COMPONENT_NEIGHBOR = 0,
  MVMC_KRYLOV_PROPOSAL_COMPONENT_GLOBAL = 1,
  MVMC_KRYLOV_PROPOSAL_COMPONENT_SHELL = 2
} MVMCKrylovPositiveSamplerProposalComponent;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  uint64_t version;
  uint64_t policy_hash;
  uint64_t kernel_id;
  size_t global_numerator;
  size_t global_denominator;
  size_t neighbor_numerator;
  size_t neighbor_denominator;
  size_t distance_numerator;
  size_t distance_denominator;
  uint64_t distance_rounding_rule;
  uint64_t subset_order;
  uint64_t self_report_rule;
  uint64_t rng_algorithm;
} MVMCKrylovPositiveSamplerProposalPolicy;

typedef struct {
  int valid;
  uint64_t state_version;
  uint64_t algorithm;
  uint64_t state;
  uint64_t stream;
  uint64_t draws;
} MVMCKrylovPositiveSamplerRng;

typedef MVMCKrylovStatus (*MVMCKrylovSurrogateLogWeightCallback)(
    const uint64_t *configuration_words, size_t word_count,
    void *context, double *log_weight);

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  uint64_t version;
  uint64_t policy_hash;
  uint64_t guide_shape_hash;
  uint64_t formula_version;
  uint64_t transaction_version;
  size_t step_count;
  int partial_order;
  size_t floor_multiplier;
  uint64_t rng_algorithm;
} MVMCKrylovPositiveSamplerSurrogatePolicy;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  double log_weight;
  double log_floor;
  double log_signal;
  int floor_only;
} MVMCKrylovPositiveSamplerSurrogateWeightResult;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  size_t step_count;
  size_t surrogate_evaluation_count;
  size_t inner_accepted;
  size_t inner_rejected;
  size_t inner_configuration_changes;
  int final_configuration_changed;
  double initial_log_weight;
  double final_log_weight;
  double log_proposal_ratio;
  double proposal_seconds;
  uint64_t surrogate_policy_hash;
  uint64_t proposal_model_hash;
  uint64_t final_configuration_hash;
  MVMCKrylovPositiveSamplerRng rng_after;
} MVMCKrylovPositiveSamplerSurrogateDrawResult;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  MVMCKrylovPositiveSamplerProposalComponent component;
  int self_proposal;
  size_t neighbor_count;
  size_t selected_neighbor_index;
  size_t component_draw_count;
  size_t global_subset_draw_count;
  size_t shell_draw_count;
  size_t shell_max_distance;
  size_t shell_distance;
  size_t shell_up_distance;
  size_t shell_down_distance;
  size_t shell_count;
  double uniform;
  double log_proposal_ratio;
  double proposal_seconds;
  double component_selection_seconds;
  double global_subset_seconds;
  double shell_generation_seconds;
  uint64_t proposal_policy_hash;
  uint64_t proposal_model_hash;
  MVMCKrylovPositiveSamplerRng rng_after;
} MVMCKrylovPositiveSamplerProposalDrawResult;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  size_t word_count;
  uint64_t policy_hash;
  uint64_t accepted_generation;
  MVMCKrylovBoundedResult krylov;
  MVMCKrylovPositiveGuideResult guide;
} MVMCKrylovPositiveSamplerSnapshot;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int accepted;
  uint64_t accepted_generation;
  MVMCKrylovBoundedResult proposal_krylov;
  MVMCKrylovPositiveGuideResult proposal_guide;
  MVMCKrylovPositiveGuideAcceptance acceptance;
} MVMCKrylovPositiveSamplerStepResult;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  MVMCKrylovPositiveSamplerProposalComponent component;
  int self_proposal;
  int configuration_changed;
  size_t neighbor_count;
  size_t selected_neighbor_index;
  size_t component_draw_count;
  size_t global_subset_draw_count;
  size_t shell_draw_count;
  size_t shell_max_distance;
  size_t shell_distance;
  size_t shell_up_distance;
  size_t shell_down_distance;
  size_t shell_count;
  double uniform;
  double log_proposal_ratio;
  double proposal_seconds;
  double component_selection_seconds;
  double global_subset_seconds;
  double shell_generation_seconds;
  double bounded_evaluation_seconds;
  double total_step_seconds;
  uint64_t proposal_policy_hash;
  uint64_t proposal_model_hash;
  MVMCKrylovPositiveSamplerRng rng_after;
  MVMCKrylovPositiveSamplerStepResult step;
} MVMCKrylovPositiveSamplerProposalStepResult;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  uint64_t state_version;
  uint64_t rng_algorithm;
  uint64_t rng_state;
  uint64_t rng_stream;
  uint64_t rng_draws;
  uint64_t policy_hash;
  uint64_t proposal_policy_hash;
  uint64_t proposal_model_hash;
  uint64_t surrogate_policy_hash;
  uint64_t guide_shape_hash;
  uint64_t bounded_plan_hash;
  uint64_t amplitude_policy_hash;
  uint64_t current_configuration_hash;
  uint64_t current_scale_hash;
  uint64_t accepted_generation;
  size_t word_count;
  size_t cache_bytes;
  int guide_order;
  int bounded_max_order;
} MVMCKrylovPositiveSamplerManifest;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  uint64_t attempted_steps;
  uint64_t completed_steps;
  uint64_t accepted_steps;
  uint64_t rejected_steps;
  uint64_t neighbor_attempted_steps;
  uint64_t neighbor_accepted_steps;
  uint64_t neighbor_rejected_steps;
  uint64_t global_attempted_steps;
  uint64_t global_accepted_steps;
  uint64_t global_rejected_steps;
  uint64_t global_self_proposals;
  uint64_t shell_attempted_steps;
  uint64_t shell_accepted_steps;
  uint64_t shell_rejected_steps;
  uint64_t configuration_changing_accepted_moves;
  uint64_t positive_support_steps;
  uint64_t support_violation_steps;
  uint64_t floor_supported_zero_steps;
  uint64_t finite_proposal_components;
  uint64_t exact_zero_proposal_components;
  uint64_t numeric_zero_proposal_components;
  uint64_t terminal_amplitude_calls;
  uint64_t trace_hash;
  double min_log_acceptance_ratio;
  double max_log_acceptance_ratio;
  double sum_log_acceptance_ratio;
  double min_proposal_log_guide;
  double max_proposal_log_guide;
} MVMCKrylovPositiveSamplerTraceStatistics;

MVMCKrylovStatus mvmc_krylov_positive_sampler_rng_seed(
    uint64_t seed, uint64_t stream,
    MVMCKrylovPositiveSamplerRng *rng);

MVMCKrylovStatus mvmc_krylov_positive_sampler_rng_draw_uint64(
    MVMCKrylovPositiveSamplerRng *rng, uint64_t *value);

MVMCKrylovStatus mvmc_krylov_positive_sampler_rng_draw_bounded(
    MVMCKrylovPositiveSamplerRng *rng, size_t bound, size_t *value);

MVMCKrylovStatus mvmc_krylov_positive_sampler_rng_draw_uniform(
    MVMCKrylovPositiveSamplerRng *rng, double *uniform);

MVMCKrylovStatus mvmc_krylov_positive_sampler_proposal_policy_create(
    size_t global_numerator, size_t global_denominator,
    MVMCKrylovPositiveSamplerProposalPolicy *policy);

MVMCKrylovStatus mvmc_krylov_positive_sampler_shell_policy_create(
    size_t neighbor_numerator, size_t neighbor_denominator,
    size_t distance_numerator, size_t distance_denominator,
    MVMCKrylovPositiveSamplerProposalPolicy *policy);

MVMCKrylovStatus mvmc_krylov_positive_sampler_surrogate_policy_create(
    size_t step_count, size_t floor_multiplier,
    const MVMCKrylovPositiveGuidePolicy *guide_policy,
    MVMCKrylovPositiveSamplerSurrogatePolicy *policy);

MVMCKrylovStatus mvmc_krylov_positive_sampler_surrogate_partial_policy_create(
    size_t step_count, int partial_order, size_t floor_multiplier,
    const MVMCKrylovPositiveGuidePolicy *guide_policy,
    MVMCKrylovPositiveSamplerSurrogatePolicy *policy);

MVMCKrylovStatus mvmc_krylov_positive_sampler_surrogate_weight_zeroth(
    const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy,
    const MVMCKrylovPositiveGuidePolicy *guide_policy,
    const MVMCScaledComplex *zeroth_value,
    MVMCKrylovPositiveSamplerSurrogateWeightResult *result);

MVMCKrylovStatus mvmc_krylov_positive_sampler_surrogate_weight_partial(
    const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy,
    const MVMCKrylovPositiveGuidePolicy *guide_policy,
    const MVMCScaledComplex *values, size_t value_count,
    MVMCKrylovPositiveSamplerSurrogateWeightResult *result);

MVMCKrylovStatus mvmc_krylov_positive_sampler_draw_surrogate_rng(
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy,
    const uint64_t *current_words, size_t current_word_count,
    double current_log_weight,
    const MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovSurrogateLogWeightCallback log_weight_callback,
    void *log_weight_context,
    uint64_t *proposal_words, size_t proposal_word_count,
    MVMCKrylovPositiveSamplerSurrogateDrawResult *result);

MVMCKrylovStatus mvmc_krylov_positive_sampler_proposal_model_hash(
    const MVMCKrylovFockModel *proposal_model,
    const uint64_t *configuration_words, size_t word_count,
    uint64_t *proposal_model_hash);

MVMCKrylovStatus mvmc_krylov_positive_sampler_proposal_model_policy_hash(
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const MVMCKrylovFockModel *proposal_model,
    const uint64_t *configuration_words, size_t word_count,
    uint64_t *proposal_model_hash);

MVMCKrylovStatus mvmc_krylov_positive_sampler_initialize(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *policy,
    const uint64_t *current_words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovPositiveSamplerSnapshot *snapshot);

MVMCKrylovStatus mvmc_krylov_positive_sampler_step(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *policy,
    uint64_t *current_words, size_t current_word_count,
    MVMCKrylovPositiveSamplerSnapshot *current,
    const uint64_t *proposal_words, size_t proposal_word_count,
    double log_proposal_ratio, double uniform,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovPositiveSamplerStepResult *result);

MVMCKrylovStatus mvmc_krylov_positive_sampler_step_selected_neighbor(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *policy,
    const MVMCKrylovFockModel *proposal_model,
    uint64_t *current_words, size_t current_word_count,
    MVMCKrylovPositiveSamplerSnapshot *current,
    size_t selected_neighbor_index, double uniform,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words, size_t proposal_word_count,
    MVMCKrylovPositiveSamplerProposalStepResult *result);

MVMCKrylovStatus mvmc_krylov_positive_sampler_step_rng(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *policy,
    const MVMCKrylovFockModel *proposal_model,
    uint64_t *current_words, size_t current_word_count,
    MVMCKrylovPositiveSamplerSnapshot *current,
    MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words, size_t proposal_word_count,
    MVMCKrylovPositiveSamplerProposalStepResult *result);

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
    MVMCKrylovPositiveSamplerProposalStepResult *result);

MVMCKrylovStatus mvmc_krylov_positive_sampler_draw_mixture_rng(
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const uint64_t *current_words, size_t current_word_count,
    const MVMCKrylovPositiveSamplerRng *rng,
    uint64_t *proposal_words, size_t proposal_word_count,
    MVMCKrylovPositiveSamplerProposalDrawResult *result);

MVMCKrylovStatus mvmc_krylov_positive_sampler_manifest_create(
    const MVMCKrylovPositiveGuidePolicy *policy,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovBoundedLimits *limits,
    uint64_t bounded_plan_hash,
    const uint64_t *current_words, size_t current_word_count,
    const MVMCKrylovPositiveSamplerSnapshot *current,
    const MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovPositiveSamplerManifest *manifest);

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
    int *matches);

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
    MVMCKrylovPositiveSamplerManifest *manifest);

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
    int *matches);

MVMCKrylovStatus mvmc_krylov_positive_sampler_trace_statistics_reset(
    MVMCKrylovPositiveSamplerTraceStatistics *statistics);

MVMCKrylovStatus mvmc_krylov_positive_sampler_trace_statistics_record_step(
    const MVMCKrylovPositiveSamplerProposalStepResult *step,
    const uint64_t *current_words, size_t current_word_count,
    MVMCKrylovPositiveSamplerTraceStatistics *statistics);

#endif

#endif /* MVMC_KRYLOV_POSITIVE_SAMPLER_H */
