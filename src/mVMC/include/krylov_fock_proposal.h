/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_KRYLOV_FOCK_PROPOSAL_H
#define MVMC_KRYLOV_FOCK_PROPOSAL_H

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                         \
    defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "krylov_fock_reference.h"

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  size_t sector_count;
  size_t visited_count;
  int connected;
} MVMCKrylovFockProposalConnectivity;

typedef MVMCKrylovStatus (*MVMCKrylovFockProposalBoundedDraw)(
    void *context, size_t bound, size_t *value);

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  size_t word_count;
  size_t active_orbital_count;
  size_t draw_count;
} MVMCKrylovFockUniformProposalResult;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  size_t word_count;
  size_t max_distance;
  size_t distance;
  size_t up_distance;
  size_t down_distance;
  size_t shell_count;
  size_t draw_count;
} MVMCKrylovFockShellProposalResult;

MVMCKrylovStatus mvmc_krylov_fock_proposal_count_neighbors(
    const MVMCKrylovFockModel *model,
    const uint64_t *configuration_words, size_t word_count,
    size_t *neighbor_count);

MVMCKrylovStatus mvmc_krylov_fock_proposal_select_neighbor(
    const MVMCKrylovFockModel *model,
    const uint64_t *configuration_words, size_t word_count,
    size_t selected_index,
    uint64_t *proposal_words, size_t proposal_word_count,
    size_t *neighbor_count);

MVMCKrylovStatus mvmc_krylov_fock_proposal_log_ratio(
    const MVMCKrylovFockModel *model,
    const uint64_t *current_words, size_t current_word_count,
    const uint64_t *proposal_words, size_t proposal_word_count,
    double *log_proposal_ratio);

MVMCKrylovStatus mvmc_krylov_fock_proposal_draw_uniform_sector(
    const MVMCKrylovFockModel *model,
    MVMCKrylovFockProposalBoundedDraw draw_bounded, void *draw_context,
    uint64_t *proposal_words, size_t proposal_word_count,
    MVMCKrylovFockUniformProposalResult *result);

MVMCKrylovStatus mvmc_krylov_fock_proposal_max_shell_distance(
    const MVMCKrylovFockModel *model, size_t *max_distance);

MVMCKrylovStatus mvmc_krylov_fock_proposal_resolve_shell_distance(
    const MVMCKrylovFockModel *model, size_t fraction_numerator,
    size_t fraction_denominator, size_t *max_distance, size_t *distance);

MVMCKrylovStatus mvmc_krylov_fock_proposal_count_shell(
    const MVMCKrylovFockModel *model, size_t distance,
    size_t *shell_count);

MVMCKrylovStatus mvmc_krylov_fock_proposal_draw_shell(
    const MVMCKrylovFockModel *model,
    const uint64_t *configuration_words, size_t word_count,
    size_t distance, MVMCKrylovFockProposalBoundedDraw draw_bounded,
    void *draw_context, uint64_t *proposal_words,
    size_t proposal_word_count, MVMCKrylovFockShellProposalResult *result);

MVMCKrylovStatus mvmc_krylov_fock_proposal_check_connectivity(
    const MVMCKrylovFockModel *model, size_t max_states,
    MVMCKrylovFockProposalConnectivity *result);

#endif /* MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE ||
          MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE */

#endif /* MVMC_KRYLOV_FOCK_PROPOSAL_H */
