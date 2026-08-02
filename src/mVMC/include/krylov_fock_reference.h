/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_KRYLOV_FOCK_REFERENCE_H
#define MVMC_KRYLOV_FOCK_REFERENCE_H

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define MVMC_KRYLOV_MAX_ORDER 3

typedef enum {
  MVMC_KRYLOV_STATUS_OK = 0,
  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
  MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL,
  MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL,
  MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE,
  MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE,
  MVMC_KRYLOV_STATUS_NONFINITE,
  MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
  MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE,
  MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE
} MVMCKrylovStatus;

typedef enum {
  MVMC_KRYLOV_FERMION_CREATE = 0,
  MVMC_KRYLOV_FERMION_ANNIHILATE
} MVMCKrylovFermionOperatorKind;

typedef struct {
  MVMCKrylovFermionOperatorKind kind;
  size_t orbital;
} MVMCKrylovFermionOperator;

typedef struct {
  double complex coefficient;
  size_t operator_offset;
  size_t operator_count;
  int source_kind;
  size_t source_index;
} MVMCKrylovHamiltonianTerm;

typedef struct {
  size_t site_count;
  size_t up_electron_count;
  size_t down_electron_count;
  int pure_spin;
  int hermitian;
  const MVMCKrylovHamiltonianTerm *terms;
  size_t term_count;
  const MVMCKrylovFermionOperator *operators;
  size_t operator_count;
} MVMCKrylovFockModel;

typedef struct {
  size_t max_states;
  size_t max_transitions;
  size_t max_amplitude_evaluations;
  size_t max_bytes;
  int max_order;
} MVMCKrylovLimits;

typedef struct {
  double complex value;
  int total_zero;
  uint64_t regular_component_count;
  uint64_t near_pivot_component_count;
  uint64_t singular_component_count;
  uint64_t local_factorization_count;
  uint64_t global_factorization_count;
} MVMCKrylovAmplitudeResult;

typedef MVMCKrylovStatus (*MVMCKrylovAmplitudeCallback)(
    const uint64_t *configuration_words, size_t word_count,
    void *context, MVMCKrylovAmplitudeResult *result);

typedef struct {
  uint64_t root_evaluations;
  int requested_order;
  int completed_order;
  uint64_t recursion_calls[MVMC_KRYLOV_MAX_ORDER + 1];
  uint64_t unique_states[MVMC_KRYLOV_MAX_ORDER + 1];
  uint64_t memo_hits[MVMC_KRYLOV_MAX_ORDER + 1];
  uint64_t memo_misses[MVMC_KRYLOV_MAX_ORDER + 1];
  uint64_t raw_transitions;
  uint64_t merged_duplicate_transitions;
  uint64_t cancelled_zero_transitions;
  uint64_t terminal_amplitude_requests;
  uint64_t terminal_cache_hits;
  uint64_t regular_component_count;
  uint64_t near_pivot_component_count;
  uint64_t singular_component_count;
  uint64_t total_zero_count;
  uint64_t local_factorization_count;
  uint64_t global_factorization_count;
  size_t frontier_peak;
  size_t memo_entries_peak;
  size_t workspace_bytes;
  double depth_wall_seconds[MVMC_KRYLOV_MAX_ORDER + 1];
  double amplitude_wall_seconds;
  double connectivity_wall_seconds;
} MVMCKrylovStatistics;

typedef struct {
  MVMCKrylovStatus status;
  int depth;
  uint64_t configuration_hash;
} MVMCKrylovFailure;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int evaluated_order;
  double complex value[MVMC_KRYLOV_MAX_ORDER + 1];
  MVMCKrylovStatistics statistics;
  MVMCKrylovFailure failure;
} MVMCKrylovResult;

typedef struct MVMCKrylovWorkspace MVMCKrylovWorkspace;

size_t mvmc_krylov_fock_word_count(size_t site_count);

MVMCKrylovWorkspace *mvmc_krylov_workspace_create(
    size_t site_count, const MVMCKrylovLimits *limits,
    MVMCKrylovStatus *status);

void mvmc_krylov_workspace_destroy(MVMCKrylovWorkspace *workspace);

size_t mvmc_krylov_workspace_bytes(const MVMCKrylovWorkspace *workspace);

MVMCKrylovStatus mvmc_krylov_fock_validate(
    const MVMCKrylovFockModel *model,
    const uint64_t *configuration_words, size_t word_count);

MVMCKrylovStatus mvmc_krylov_evaluate(
    MVMCKrylovWorkspace *workspace, const MVMCKrylovFockModel *model,
    const uint64_t *root_configuration_words, size_t root_word_count,
    MVMCKrylovAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovResult *result);

const char *mvmc_krylov_status_string(MVMCKrylovStatus status);

#endif /* MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE */

#endif /* MVMC_KRYLOV_FOCK_REFERENCE_H */
