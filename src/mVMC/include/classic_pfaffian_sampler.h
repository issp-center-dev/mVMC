/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_CLASSIC_PFAFFIAN_SAMPLER_H
#define MVMC_CLASSIC_PFAFFIAN_SAMPLER_H

#include "classic_pfaffian_collective.h"
#if defined(MVMC_ENABLE_PFAFFIAN_STATE_AUDIT)
#include "classic_pfaffian_audit.h"
#endif

#include <stdint.h>

#if defined(MVMC_ENABLE_PFAFFIAN_SAMPLER_REFERENCE)

typedef enum {
  MVMC_CLASSIC_MOVE_HOPPING = 0,
  MVMC_CLASSIC_MOVE_SPIN_HOPPING,
  MVMC_CLASSIC_MOVE_EXCHANGE,
  MVMC_CLASSIC_MOVE_PAIR_HOPPING
} MVMCClassicMoveKind;

typedef enum {
  MVMC_CLASSIC_PROPOSAL_ERROR = 0,
  MVMC_CLASSIC_PROPOSAL_REJECTED,
  MVMC_CLASSIC_PROPOSAL_ACCEPTED
} MVMCClassicProposalDecision;

typedef double (*MVMCClassicUniformDraw)(void *context);

typedef struct {
  int valid;
  MVMCProjectedAmplitudeResult accepted_aggregate;
} MVMCClassicSamplerState;

typedef struct {
  int valid;
  MVMCPfaffianStatus status;
  MVMCClassicMoveKind move;
  MVMCClassicProposalDecision decision;
  int uniform_draw_count;
  double uniform;
  double log_acceptance_ratio;
  double complex current_total;
  double complex proposal_total;
  uint64_t accepted_generation;
#if defined(MVMC_ENABLE_PFAFFIAN_STATE_AUDIT)
  MVMCClassicPfaffianAuditReport audit;
#endif
} MVMCClassicProposalResult;

#if defined(MVMC_ENABLE_PFAFFIAN_STATE_AUDIT)

#define MVMC_CLASSIC_TRACE_INDEX_CAPACITY 4

typedef struct {
  /* Diagnostic identity only; it never affects state or acceptance. */
  uint64_t proposal_id;
  int index_count;
  int indices[MVMC_CLASSIC_TRACE_INDEX_CAPACITY];
} MVMCClassicProposalMetadata;

typedef struct {
  int valid;
  uint64_t proposal_id;
  MVMCClassicMoveKind move;
  int index_count;
  int indices[MVMC_CLASSIC_TRACE_INDEX_CAPACITY];
  MVMCClassicPfaffianAuditMode requested_audit_mode;
  MVMCClassicPfaffianAuditMode effective_audit_mode;
  int audit_executed;
  int audit_fallback;
  int uniform_draw_count;
  double uniform;
  double complex current_total;
  double complex proposal_total;
  size_t regular_count;
  size_t near_singular_count;
  size_t singular_count;
  MVMCClassicProposalDecision decision;
  uint64_t accepted_generation;
} MVMCClassicProposalTrace;

typedef struct {
  const MVMCClassicProposalMetadata *metadata;
  MVMCClassicPfaffianRealAuditWorkspace *workspace;
  MVMCClassicPfaffianRealObserver observer;
  void *observer_context;
  MVMCClassicProposalTrace *trace;
} MVMCClassicRealProposalAudit;

typedef struct {
  const MVMCClassicProposalMetadata *metadata;
  MVMCClassicPfaffianComplexAuditWorkspace *workspace;
  MVMCClassicPfaffianComplexObserver observer;
  void *observer_context;
  MVMCClassicProposalTrace *trace;
} MVMCClassicComplexProposalAudit;

#endif /* MVMC_ENABLE_PFAFFIAN_STATE_AUDIT */

/*
 * These entry points exist only in Testing builds.  The caller has already
 * applied a candidate configuration before propose().  On rejection it must
 * revert that configuration; on acceptance the absolute candidate and its
 * legacy mirror are published atomically.  A proposal that reaches this API
 * consumes exactly one uniform draw, including an exact-zero total.
 */
MVMCPfaffianStatus mvmc_classic_sampler_real_initialize(
    MVMCClassicPfaffianRealWorkspace *state_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const double *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance,
    double *legacy_storage, size_t legacy_element_count,
    MVMCClassicSamplerState *sampler_state);

MVMCPfaffianStatus mvmc_classic_sampler_complex_initialize(
    MVMCClassicPfaffianComplexWorkspace *state_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const double complex *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance,
    double complex *legacy_storage, size_t legacy_element_count,
    MVMCClassicSamplerState *sampler_state);

MVMCPfaffianStatus mvmc_classic_sampler_real_propose(
    MVMCClassicPfaffianRealWorkspace *state_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    MVMCClassicSamplerState *sampler_state, MVMCClassicMoveKind move,
    const double *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance,
    double log_projection_ratio, double log_rbm_ratio,
    MVMCClassicUniformDraw draw_uniform, void *draw_context,
    double *legacy_storage, size_t legacy_element_count,
    MVMCClassicProposalResult *result);

MVMCPfaffianStatus mvmc_classic_sampler_complex_propose(
    MVMCClassicPfaffianComplexWorkspace *state_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    MVMCClassicSamplerState *sampler_state, MVMCClassicMoveKind move,
    const double complex *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance,
    double log_projection_ratio, double log_rbm_ratio,
    MVMCClassicUniformDraw draw_uniform, void *draw_context,
    double complex *legacy_storage, size_t legacy_element_count,
    MVMCClassicProposalResult *result);

#if defined(MVMC_ENABLE_PFAFFIAN_STATE_AUDIT)

/*
 * Audit entry points are collective.  All ranks must supply the same audit
 * mode and valid metadata shape.  Metadata values are caller-owned diagnostic
 * labels and do not participate in the Metropolis decision.
 */

MVMCPfaffianStatus mvmc_classic_sampler_real_propose_with_audit(
    MVMCClassicPfaffianRealWorkspace *state_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    MVMCClassicSamplerState *sampler_state, MVMCClassicMoveKind move,
    const double *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance,
    double log_projection_ratio, double log_rbm_ratio,
    MVMCClassicUniformDraw draw_uniform, void *draw_context,
    const MVMCClassicRealProposalAudit *audit,
    double *legacy_storage, size_t legacy_element_count,
    MVMCClassicProposalResult *result);

MVMCPfaffianStatus mvmc_classic_sampler_complex_propose_with_audit(
    MVMCClassicPfaffianComplexWorkspace *state_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    MVMCClassicSamplerState *sampler_state, MVMCClassicMoveKind move,
    const double complex *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance,
    double log_projection_ratio, double log_rbm_ratio,
    MVMCClassicUniformDraw draw_uniform, void *draw_context,
    const MVMCClassicComplexProposalAudit *audit,
    double complex *legacy_storage, size_t legacy_element_count,
    MVMCClassicProposalResult *result);

/* Re-evaluate accepted configuration, compare, and always discard candidate. */
MVMCPfaffianStatus mvmc_classic_sampler_real_rebuild_audit(
    MVMCClassicPfaffianRealWorkspace *state_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    MVMCClassicSamplerState *sampler_state,
    const double *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCClassicPfaffianRealAuditWorkspace *audit_workspace,
    MVMCClassicPfaffianAuditReport *report);

MVMCPfaffianStatus mvmc_classic_sampler_complex_rebuild_audit(
    MVMCClassicPfaffianComplexWorkspace *state_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    MVMCClassicSamplerState *sampler_state,
    const double complex *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCClassicPfaffianComplexAuditWorkspace *audit_workspace,
    MVMCClassicPfaffianAuditReport *report);

#endif /* MVMC_ENABLE_PFAFFIAN_STATE_AUDIT */

#endif /* MVMC_ENABLE_PFAFFIAN_SAMPLER_REFERENCE */

#endif /* MVMC_CLASSIC_PFAFFIAN_SAMPLER_H */
