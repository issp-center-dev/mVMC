/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "classic_pfaffian_sampler.h"

#if !defined(MVMC_ENABLE_PFAFFIAN_SAMPLER_REFERENCE)
#error "classic_pfaffian_sampler.c is Testing-only"
#endif
#if !defined(MVMC_ENABLE_PFAFFIAN_STATE_AUDIT)
#error "classic_pfaffian_sampler.c requires the Testing-only audit module"
#endif

#include <math.h>
#include <string.h>

enum {
  MVMC_CLASSIC_BRANCH_INITIALIZE = 0,
  MVMC_CLASSIC_BRANCH_PROPOSE_BASE = 1,
  MVMC_CLASSIC_BRANCH_PROPOSE_STRIDE = 3,
  MVMC_CLASSIC_BRANCH_REBUILD_REAL = 13,
  MVMC_CLASSIC_BRANCH_REBUILD_COMPLEX = 16
};

static int valid_move(MVMCClassicMoveKind move) {
  return move >= MVMC_CLASSIC_MOVE_HOPPING &&
         move <= MVMC_CLASSIC_MOVE_PAIR_HOPPING;
}

static void reset_sampler_state(MVMCClassicSamplerState *state) {
  if (state != NULL) memset(state, 0, sizeof(*state));
}

static void reset_proposal_result(MVMCClassicProposalResult *result,
                                  MVMCClassicMoveKind move) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  result->move = move;
  result->decision = MVMC_CLASSIC_PROPOSAL_ERROR;
  result->uniform = NAN;
  result->log_acceptance_ratio = -INFINITY;
  result->audit.status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  result->audit.failure_rank = -1;
  result->audit.failure_local_qp = -1;
  result->audit.failure_global_qp = -1;
}

static int valid_metadata(const MVMCClassicProposalMetadata *metadata) {
  return metadata == NULL ||
         (metadata->index_count >= 0 &&
          metadata->index_count <= MVMC_CLASSIC_TRACE_INDEX_CAPACITY);
}

static void reset_trace(MVMCClassicProposalTrace *trace,
                        MVMCClassicMoveKind move,
                        const MVMCClassicProposalMetadata *metadata,
                        MVMCClassicPfaffianAuditMode requested_mode) {
  int index;

  if (trace == NULL) return;
  memset(trace, 0, sizeof(*trace));
  trace->move = move;
  trace->requested_audit_mode = requested_mode;
  trace->effective_audit_mode = MVMC_CLASSIC_AUDIT_OFF;
  trace->uniform = NAN;
  if (metadata == NULL || !valid_metadata(metadata)) return;
  trace->proposal_id = metadata->proposal_id;
  trace->index_count = metadata->index_count;
  for (index = 0; index < metadata->index_count; ++index) {
    trace->indices[index] = metadata->indices[index];
  }
}

static void finish_trace(MVMCClassicProposalTrace *trace,
                         const MVMCClassicProposalResult *result,
                         const MVMCProjectedAmplitudeResult *proposal) {
  if (trace == NULL || result == NULL || proposal == NULL) return;
  trace->effective_audit_mode = result->audit.effective_mode;
  trace->audit_executed = result->audit.executed;
  trace->audit_fallback = result->audit.fallback;
  trace->uniform_draw_count = result->uniform_draw_count;
  trace->uniform = result->uniform;
  trace->current_total = result->current_total;
  trace->proposal_total = result->proposal_total;
  trace->regular_count = proposal->regular_count;
  trace->near_singular_count = proposal->near_singular_count;
  trace->singular_count = proposal->singular_count;
  trace->decision = result->decision;
  trace->accepted_generation = result->accepted_generation;
  trace->valid = result->valid;
}

MVMCPfaffianStatus mvmc_classic_sampler_real_initialize(
    MVMCClassicPfaffianRealWorkspace *state_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const double *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance,
    double *legacy_storage, size_t legacy_element_count,
    MVMCClassicSamplerState *sampler_state) {
  MVMCClassicPfaffianCollectiveResult prepared;
  MVMCPfaffianStatus status;

  reset_sampler_state(sampler_state);
  status = mvmc_classic_pfaffian_collective_preflight(
      collective_workspace,
      sampler_state != NULL && state_workspace != NULL &&
          legacy_storage != NULL &&
          legacy_element_count >=
              mvmc_classic_pfaffian_real_legacy_element_count(
                  state_workspace),
      MVMC_CLASSIC_BRANCH_INITIALIZE);
  if (status != MVMC_PFAFFIAN_STATUS_OK) return status;
  status = mvmc_classic_pfaffian_real_prepare_collective(
      state_workspace, collective_workspace, slater_elm, ele_idx,
      global_weights, scaled_pivot_tolerance, residual_tolerance, &prepared);
  if (status != MVMC_PFAFFIAN_STATUS_OK) return status;
  if (cabs(prepared.aggregate.total) == 0.0) {
    mvmc_classic_pfaffian_real_discard_candidate(state_workspace);
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_classic_pfaffian_real_publish(
      state_workspace, legacy_storage, legacy_element_count);
  if (status != MVMC_PFAFFIAN_STATUS_OK) return status;
  sampler_state->accepted_aggregate = prepared.aggregate;
  sampler_state->valid = 1;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_classic_sampler_complex_initialize(
    MVMCClassicPfaffianComplexWorkspace *state_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const double complex *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance,
    double complex *legacy_storage, size_t legacy_element_count,
    MVMCClassicSamplerState *sampler_state) {
  MVMCClassicPfaffianCollectiveResult prepared;
  MVMCPfaffianStatus status;

  reset_sampler_state(sampler_state);
  status = mvmc_classic_pfaffian_collective_preflight(
      collective_workspace,
      sampler_state != NULL && state_workspace != NULL &&
          legacy_storage != NULL &&
          legacy_element_count >=
              mvmc_classic_pfaffian_complex_legacy_element_count(
                  state_workspace),
      MVMC_CLASSIC_BRANCH_INITIALIZE);
  if (status != MVMC_PFAFFIAN_STATUS_OK) return status;
  status = mvmc_classic_pfaffian_complex_prepare_collective(
      state_workspace, collective_workspace, slater_elm, ele_idx,
      global_weights, scaled_pivot_tolerance, residual_tolerance, &prepared);
  if (status != MVMC_PFAFFIAN_STATUS_OK) return status;
  if (cabs(prepared.aggregate.total) == 0.0) {
    mvmc_classic_pfaffian_complex_discard_candidate(state_workspace);
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_classic_pfaffian_complex_publish(
      state_workspace, legacy_storage, legacy_element_count);
  if (status != MVMC_PFAFFIAN_STATUS_OK) return status;
  sampler_state->accepted_aggregate = prepared.aggregate;
  sampler_state->valid = 1;
  return MVMC_PFAFFIAN_STATUS_OK;
}

#define DEFINE_PROPOSE_INTERNAL(                                               \
    NAME, TYPE, AUDIT_CONTROL, PREPARE, PUBLISH, DISCARD, ACCEPTED,            \
    CANDIDATE, CANDIDATE_INVERSE, LEGACY_COUNT, AUDIT_MODE,                    \
    AUDIT_COMPATIBLE, AUDIT_CANDIDATE)                                         \
static MVMCPfaffianStatus NAME(                                                \
    TYPE *state_workspace,                                                    \
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,             \
    MVMCClassicSamplerState *sampler_state, MVMCClassicMoveKind move,          \
    const ELEMENT_TYPE *slater_elm, const int *ele_idx,                        \
    const double complex *global_weights,                                     \
    double scaled_pivot_tolerance, double residual_tolerance,                  \
    double log_projection_ratio, double log_rbm_ratio,                         \
    MVMCClassicUniformDraw draw_uniform, void *draw_context,                   \
    const AUDIT_CONTROL *audit,                                                \
    ELEMENT_TYPE *legacy_storage, size_t legacy_element_count,                 \
    MVMCClassicProposalResult *result) {                                       \
  MVMCClassicPfaffianCollectiveResult prepared;                                \
  MVMCClassicPfaffianMetropolisResult metropolis;                              \
  const MVMCClassicPfaffianState *accepted;                                    \
  const MVMCClassicProposalMetadata *metadata =                                \
      audit == NULL ? NULL : audit->metadata;                                  \
  MVMCClassicProposalTrace *trace = audit == NULL ? NULL : audit->trace;       \
  MVMCClassicPfaffianAuditMode audit_mode =                                    \
      audit == NULL ? MVMC_CLASSIC_AUDIT_OFF : AUDIT_MODE(audit->workspace);   \
  MVMCPfaffianStatus status;                                                   \
  int branch_kind;                                                             \
  int audit_valid;                                                             \
  double uniform;                                                              \
  reset_proposal_result(result, move);                                         \
  reset_trace(trace, move, metadata, audit_mode);                              \
  audit_valid = audit == NULL ||                                               \
      (audit->workspace != NULL && valid_metadata(metadata) &&                 \
       AUDIT_COMPATIBLE(audit->workspace, state_workspace) &&                  \
       (audit_mode == MVMC_CLASSIC_AUDIT_OFF || audit->observer != NULL));      \
  branch_kind = valid_move(move)                                               \
      ? MVMC_CLASSIC_BRANCH_PROPOSE_BASE +                                     \
            (int)move * MVMC_CLASSIC_BRANCH_PROPOSE_STRIDE + (int)audit_mode   \
      : MVMC_CLASSIC_BRANCH_PROPOSE_BASE;                                      \
  status = mvmc_classic_pfaffian_collective_preflight(                         \
      collective_workspace,                                                   \
      state_workspace != NULL && sampler_state != NULL &&                      \
          sampler_state->valid && result != NULL && valid_move(move) &&        \
          draw_uniform != NULL && !isnan(log_projection_ratio) &&              \
          !isnan(log_rbm_ratio) && legacy_storage != NULL && audit_valid &&    \
          legacy_element_count >= LEGACY_COUNT(state_workspace),               \
      branch_kind);                                                            \
  if (status != MVMC_PFAFFIAN_STATUS_OK) return status;                        \
  /* Keep a local guard after collective preflight: it is both defensive and  \
   * makes the pointer contract explicit to static analysis. */                \
  if (state_workspace == NULL || sampler_state == NULL ||                      \
      !sampler_state->valid || result == NULL || !valid_move(move) ||          \
      draw_uniform == NULL || isnan(log_projection_ratio) ||                   \
      isnan(log_rbm_ratio) || legacy_storage == NULL || !audit_valid ||        \
      legacy_element_count < LEGACY_COUNT(state_workspace)) {                  \
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;                              \
  }                                                                            \
  result->current_total = sampler_state->accepted_aggregate.total;             \
  status = PREPARE(state_workspace, collective_workspace, slater_elm, ele_idx,\
                   global_weights, scaled_pivot_tolerance,                     \
                   residual_tolerance, &prepared);                             \
  if (status != MVMC_PFAFFIAN_STATUS_OK) {                                     \
    result->status = status;                                                   \
    return status;                                                             \
  }                                                                            \
  result->proposal_total = prepared.aggregate.total;                           \
  if (audit != NULL) {                                                         \
    status = AUDIT_CANDIDATE(                                                  \
        audit->workspace, collective_workspace, ACCEPTED(state_workspace),    \
        CANDIDATE(state_workspace), CANDIDATE_INVERSE(state_workspace),        \
        &prepared.aggregate, global_weights, audit->observer,                  \
        audit->observer_context, &result->audit);                              \
    if (status != MVMC_PFAFFIAN_STATUS_OK) {                                   \
      DISCARD(state_workspace);                                                \
      result->status = status;                                                 \
      return status;                                                           \
    }                                                                          \
  } else {                                                                     \
    result->audit.valid = 1;                                                   \
    result->audit.status = MVMC_PFAFFIAN_STATUS_OK;                            \
    result->audit.requested_mode = MVMC_CLASSIC_AUDIT_OFF;                     \
    result->audit.effective_mode = MVMC_CLASSIC_AUDIT_OFF;                     \
  }                                                                            \
  uniform = draw_uniform(draw_context);                                        \
  result->uniform = uniform;                                                   \
  result->uniform_draw_count = 1;                                              \
  status = mvmc_classic_pfaffian_collective_metropolis(                        \
      collective_workspace, result->current_total, result->proposal_total,    \
      log_projection_ratio + log_rbm_ratio, uniform, &metropolis);             \
  if (status != MVMC_PFAFFIAN_STATUS_OK) {                                     \
    DISCARD(state_workspace);                                                  \
    result->status = status;                                                   \
    return status;                                                             \
  }                                                                            \
  result->log_acceptance_ratio = metropolis.log_acceptance_ratio;              \
  if (metropolis.accepted) {                                                   \
    /* Preflight validated mirror capacity; prepare validated candidate. */    \
    status = PUBLISH(state_workspace, legacy_storage, legacy_element_count);   \
    if (status != MVMC_PFAFFIAN_STATUS_OK) {                                   \
      result->status = status;                                                 \
      return status;                                                           \
    }                                                                          \
    sampler_state->accepted_aggregate = prepared.aggregate;                    \
    result->decision = MVMC_CLASSIC_PROPOSAL_ACCEPTED;                         \
  } else {                                                                     \
    DISCARD(state_workspace);                                                  \
    result->decision = MVMC_CLASSIC_PROPOSAL_REJECTED;                         \
  }                                                                            \
  accepted = ACCEPTED(state_workspace);                                        \
  result->accepted_generation = accepted == NULL ? 0 : accepted->generation;   \
  result->valid = 1;                                                           \
  result->status = MVMC_PFAFFIAN_STATUS_OK;                                    \
  finish_trace(trace, result, &prepared.aggregate);                            \
  return MVMC_PFAFFIAN_STATUS_OK;                                              \
}

#define ELEMENT_TYPE double
DEFINE_PROPOSE_INTERNAL(
    propose_real_internal, MVMCClassicPfaffianRealWorkspace,
    MVMCClassicRealProposalAudit,
    mvmc_classic_pfaffian_real_prepare_collective,
    mvmc_classic_pfaffian_real_publish,
    mvmc_classic_pfaffian_real_discard_candidate,
    mvmc_classic_pfaffian_real_accepted,
    mvmc_classic_pfaffian_real_candidate,
    mvmc_classic_pfaffian_real_candidate_inverse,
    mvmc_classic_pfaffian_real_legacy_element_count,
    mvmc_classic_pfaffian_real_audit_mode,
    mvmc_classic_pfaffian_real_audit_compatible,
    mvmc_classic_pfaffian_real_audit_candidate)
#undef ELEMENT_TYPE

#define ELEMENT_TYPE double complex
DEFINE_PROPOSE_INTERNAL(
    propose_complex_internal, MVMCClassicPfaffianComplexWorkspace,
    MVMCClassicComplexProposalAudit,
    mvmc_classic_pfaffian_complex_prepare_collective,
    mvmc_classic_pfaffian_complex_publish,
    mvmc_classic_pfaffian_complex_discard_candidate,
    mvmc_classic_pfaffian_complex_accepted,
    mvmc_classic_pfaffian_complex_candidate,
    mvmc_classic_pfaffian_complex_candidate_inverse,
    mvmc_classic_pfaffian_complex_legacy_element_count,
    mvmc_classic_pfaffian_complex_audit_mode,
    mvmc_classic_pfaffian_complex_audit_compatible,
    mvmc_classic_pfaffian_complex_audit_candidate)
#undef ELEMENT_TYPE
#undef DEFINE_PROPOSE_INTERNAL

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
    MVMCClassicProposalResult *result) {
  return propose_real_internal(
      state_workspace, collective_workspace, sampler_state, move, slater_elm,
      ele_idx, global_weights, scaled_pivot_tolerance, residual_tolerance,
      log_projection_ratio, log_rbm_ratio, draw_uniform, draw_context, NULL,
      legacy_storage, legacy_element_count, result);
}

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
    MVMCClassicProposalResult *result) {
  return propose_real_internal(
      state_workspace, collective_workspace, sampler_state, move, slater_elm,
      ele_idx, global_weights, scaled_pivot_tolerance, residual_tolerance,
      log_projection_ratio, log_rbm_ratio, draw_uniform, draw_context, audit,
      legacy_storage, legacy_element_count, result);
}

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
    MVMCClassicProposalResult *result) {
  return propose_complex_internal(
      state_workspace, collective_workspace, sampler_state, move, slater_elm,
      ele_idx, global_weights, scaled_pivot_tolerance, residual_tolerance,
      log_projection_ratio, log_rbm_ratio, draw_uniform, draw_context, NULL,
      legacy_storage, legacy_element_count, result);
}

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
    MVMCClassicProposalResult *result) {
  return propose_complex_internal(
      state_workspace, collective_workspace, sampler_state, move, slater_elm,
      ele_idx, global_weights, scaled_pivot_tolerance, residual_tolerance,
      log_projection_ratio, log_rbm_ratio, draw_uniform, draw_context, audit,
      legacy_storage, legacy_element_count, result);
}

#define DEFINE_REBUILD_AUDIT(                                                  \
    NAME, TYPE, ELEMENT_TYPE, AUDIT_TYPE, PREPARE, ACCEPTED,                   \
    ACCEPTED_INVERSE, CANDIDATE, CANDIDATE_INVERSE, DISCARD,                   \
    AUDIT_MODE, AUDIT_COMPATIBLE, AUDIT_COMPARE, BRANCH_BASE)                  \
MVMCPfaffianStatus NAME(                                                       \
    TYPE *state_workspace,                                                     \
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,             \
    MVMCClassicSamplerState *sampler_state,                                    \
    const ELEMENT_TYPE *slater_elm, const int *ele_idx,                        \
    const double complex *global_weights,                                      \
    double scaled_pivot_tolerance, double residual_tolerance,                  \
    AUDIT_TYPE *audit_workspace,                                               \
    MVMCClassicPfaffianAuditReport *report) {                                  \
  MVMCClassicPfaffianCollectiveResult prepared;                                \
  MVMCClassicPfaffianAuditMode mode = AUDIT_MODE(audit_workspace);             \
  MVMCPfaffianStatus status;                                                   \
  status = mvmc_classic_pfaffian_collective_preflight(                         \
      collective_workspace,                                                   \
      state_workspace != NULL && sampler_state != NULL &&                      \
          sampler_state->valid && audit_workspace != NULL && report != NULL && \
          AUDIT_COMPATIBLE(audit_workspace, state_workspace),                  \
      BRANCH_BASE + (int)mode);                                                \
  if (status != MVMC_PFAFFIAN_STATUS_OK) return status;                        \
  if (state_workspace == NULL || sampler_state == NULL ||                      \
      !sampler_state->valid || audit_workspace == NULL || report == NULL ||    \
      !AUDIT_COMPATIBLE(audit_workspace, state_workspace)) {                   \
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;                              \
  }                                                                            \
  status = PREPARE(state_workspace, collective_workspace, slater_elm, ele_idx,\
                   global_weights, scaled_pivot_tolerance,                     \
                   residual_tolerance, &prepared);                             \
  if (status != MVMC_PFAFFIAN_STATUS_OK) return status;                        \
  status = AUDIT_COMPARE(                                                      \
      audit_workspace, collective_workspace, ACCEPTED(state_workspace),       \
      ACCEPTED_INVERSE(state_workspace),                                       \
      &sampler_state->accepted_aggregate, CANDIDATE(state_workspace),          \
      CANDIDATE_INVERSE(state_workspace), report);                             \
  DISCARD(state_workspace);                                                    \
  return status;                                                               \
}

DEFINE_REBUILD_AUDIT(
    mvmc_classic_sampler_real_rebuild_audit,
    MVMCClassicPfaffianRealWorkspace, double,
    MVMCClassicPfaffianRealAuditWorkspace,
    mvmc_classic_pfaffian_real_prepare_collective,
    mvmc_classic_pfaffian_real_accepted,
    mvmc_classic_pfaffian_real_accepted_inverse,
    mvmc_classic_pfaffian_real_candidate,
    mvmc_classic_pfaffian_real_candidate_inverse,
    mvmc_classic_pfaffian_real_discard_candidate,
    mvmc_classic_pfaffian_real_audit_mode,
    mvmc_classic_pfaffian_real_audit_compatible,
    mvmc_classic_pfaffian_real_audit_compare,
    MVMC_CLASSIC_BRANCH_REBUILD_REAL)

DEFINE_REBUILD_AUDIT(
    mvmc_classic_sampler_complex_rebuild_audit,
    MVMCClassicPfaffianComplexWorkspace, double complex,
    MVMCClassicPfaffianComplexAuditWorkspace,
    mvmc_classic_pfaffian_complex_prepare_collective,
    mvmc_classic_pfaffian_complex_accepted,
    mvmc_classic_pfaffian_complex_accepted_inverse,
    mvmc_classic_pfaffian_complex_candidate,
    mvmc_classic_pfaffian_complex_candidate_inverse,
    mvmc_classic_pfaffian_complex_discard_candidate,
    mvmc_classic_pfaffian_complex_audit_mode,
    mvmc_classic_pfaffian_complex_audit_compatible,
    mvmc_classic_pfaffian_complex_audit_compare,
    MVMC_CLASSIC_BRANCH_REBUILD_COMPLEX)

#undef DEFINE_REBUILD_AUDIT
