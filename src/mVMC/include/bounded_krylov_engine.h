/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_BOUNDED_KRYLOV_ENGINE_H
#define MVMC_BOUNDED_KRYLOV_ENGINE_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "absolute_pfaffian.h"
#include "krylov_fock_reference.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define MVMC_BOUNDED_KRYLOV_CACHE_WAYS 4

typedef struct {
  uint64_t amplitude_policy_hash;
  size_t cache_bytes;
  size_t max_row_transitions;
  size_t max_workspace_bytes;
  uint64_t max_node_expansions;
  uint64_t max_terminal_amplitude_calls;
  uint64_t max_total_row_transitions;
  int max_order;
} MVMCKrylovBoundedLimits;

typedef struct {
  MVMCScaledComplex value;
  uint64_t well_pivoted_component_count;
  uint64_t near_pivot_component_count;
  uint64_t exact_zero_component_count;
  uint64_t numeric_zero_component_count;
  uint64_t local_factorization_count;
  uint64_t global_factorization_count;
} MVMCKrylovScaledAmplitudeResult;

typedef MVMCKrylovStatus (*MVMCKrylovScaledAmplitudeCallback)(
    const uint64_t *configuration_words, size_t word_count,
    void *context, MVMCKrylovScaledAmplitudeResult *result);

typedef struct {
  uint64_t root_evaluations;
  int requested_order;
  int completed_order;
  uint64_t node_expansions;
  uint64_t recursion_calls[MVMC_KRYLOV_MAX_ORDER + 1];
  uint64_t cache_hits[MVMC_KRYLOV_MAX_ORDER + 1];
  uint64_t cache_misses[MVMC_KRYLOV_MAX_ORDER + 1];
  uint64_t cache_insertions;
  uint64_t cache_evictions;
  uint64_t cache_epoch_full_clears;
  uint64_t raw_row_transitions;
  uint64_t merged_duplicate_transitions;
  uint64_t cancelled_zero_transitions;
  uint64_t terminal_amplitude_calls;
  uint64_t well_pivoted_component_count;
  uint64_t near_pivot_component_count;
  uint64_t exact_zero_component_count;
  uint64_t numeric_zero_component_count;
  uint64_t local_factorization_count;
  uint64_t global_factorization_count;
  uint64_t computed_exact_zero_values;
  uint64_t computed_numeric_zero_values;
  uint64_t trace_hash;
  size_t row_transition_peak;
  size_t cache_entries_peak;
  size_t plan_bytes;
  size_t model_bytes;
  size_t workspace_bytes;
  size_t frame_bytes;
  size_t scratch_bytes;
  size_t cache_requested_bytes;
  size_t cache_allocated_bytes;
  size_t cache_set_count;
  uint64_t engine_heap_allocations;
  double reset_seconds;
  double total_seconds;
  double evaluation_wall_seconds;
  double depth_wall_seconds[MVMC_KRYLOV_MAX_ORDER + 1];
  double amplitude_wall_seconds;
  double connectivity_wall_seconds;
  double cache_wall_seconds;
} MVMCKrylovBoundedStatistics;

typedef struct {
  MVMCKrylovStatus status;
  int depth;
  uint64_t configuration_hash;
} MVMCKrylovBoundedFailure;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int evaluated_order;
  MVMCScaledComplex value[MVMC_KRYLOV_MAX_ORDER + 1];
  MVMCKrylovBoundedStatistics statistics;
  MVMCKrylovBoundedFailure failure;
} MVMCKrylovBoundedResult;

typedef struct MVMCKrylovBoundedPlan MVMCKrylovBoundedPlan;
typedef struct MVMCKrylovBoundedWorkspace MVMCKrylovBoundedWorkspace;

MVMCKrylovStatus mvmc_bounded_krylov_plan_create(
    const MVMCKrylovFockModel *model,
    const MVMCKrylovBoundedLimits *limits,
    MVMCKrylovBoundedPlan **plan);
void mvmc_bounded_krylov_plan_destroy(MVMCKrylovBoundedPlan *plan);
size_t mvmc_bounded_krylov_plan_bytes(const MVMCKrylovBoundedPlan *plan);
size_t mvmc_bounded_krylov_plan_model_bytes(
    const MVMCKrylovBoundedPlan *plan);
size_t mvmc_bounded_krylov_plan_max_row_transitions(
    const MVMCKrylovBoundedPlan *plan);
size_t mvmc_bounded_krylov_cache_set_bytes(
    const MVMCKrylovBoundedPlan *plan);
uint64_t mvmc_bounded_krylov_plan_hash(
    const MVMCKrylovBoundedPlan *plan);

MVMCKrylovStatus mvmc_bounded_krylov_workspace_create(
    const MVMCKrylovBoundedPlan *plan,
    MVMCKrylovBoundedWorkspace **workspace);
void mvmc_bounded_krylov_workspace_destroy(
    MVMCKrylovBoundedWorkspace *workspace);
size_t mvmc_bounded_krylov_workspace_bytes(
    const MVMCKrylovBoundedWorkspace *workspace);

/*
 * Success commits v[0..max_order].  Failure commits an invalid diagnostic
 * snapshot whose value array is sanitized to exact zero; no completed prefix
 * is exposed as a mathematical result.
 */
MVMCKrylovStatus mvmc_bounded_krylov_evaluate(
    MVMCKrylovBoundedWorkspace *workspace,
    const uint64_t *root_configuration_words, size_t root_word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovBoundedResult *result);

/* Testing hooks for deterministic wrap and failure fixtures. */
MVMCKrylovStatus mvmc_bounded_krylov_testing_force_cache_counters(
    MVMCKrylovBoundedWorkspace *workspace,
    uint64_t epoch, uint64_t access_counter);

#endif /* MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE */

#endif /* MVMC_BOUNDED_KRYLOV_ENGINE_H */
