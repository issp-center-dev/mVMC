/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "classic_pfaffian_audit.h"

#if !defined(MVMC_ENABLE_PFAFFIAN_STATE_AUDIT)
#error "classic_pfaffian_audit.c is Testing-only"
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  MVMCClassicPfaffianAuditMode mode;
  int nsize;
  int qp_total;
  int qp_start;
  int qp_end;
  size_t local_qp_count;
  size_t matrix_count;
  size_t inverse_count;
  double absolute_tolerance;
  double relative_tolerance;
  MVMCAbsolutePfaffianResult *components;
} MVMCClassicPfaffianAuditLayout;

struct MVMCClassicPfaffianRealAuditWorkspace {
  MVMCClassicPfaffianAuditLayout layout;
  double *inverse;
};

struct MVMCClassicPfaffianComplexAuditWorkspace {
  MVMCClassicPfaffianAuditLayout layout;
  double complex *inverse;
};

static int valid_mode(MVMCClassicPfaffianAuditMode mode) {
  return mode >= MVMC_CLASSIC_AUDIT_OFF &&
         mode <= MVMC_CLASSIC_AUDIT_ALWAYS;
}

static int checked_multiply(size_t left, size_t right, size_t *product) {
  if (product == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *product = left * right;
  return 1;
}

static int initialize_layout(int nsize, int qp_total, int qp_start, int qp_end,
                             MVMCClassicPfaffianAuditMode mode,
                             double absolute_tolerance,
                             double relative_tolerance,
                             MVMCClassicPfaffianAuditLayout *layout) {
  if (layout == NULL || !valid_mode(mode) || nsize <= 0 ||
      qp_total <= 0 || qp_start < 0 || qp_end < qp_start ||
      qp_end > qp_total || !isfinite(absolute_tolerance) ||
      !isfinite(relative_tolerance) || absolute_tolerance < 0.0 ||
      relative_tolerance < 0.0 ||
      !checked_multiply((size_t)nsize, (size_t)nsize,
                        &layout->matrix_count) ||
      !checked_multiply((size_t)(qp_end - qp_start), layout->matrix_count,
                        &layout->inverse_count)) {
    return 0;
  }
  layout->mode = mode;
  layout->nsize = nsize;
  layout->qp_total = qp_total;
  layout->qp_start = qp_start;
  layout->qp_end = qp_end;
  layout->local_qp_count = (size_t)(qp_end - qp_start);
  layout->absolute_tolerance = absolute_tolerance;
  layout->relative_tolerance = relative_tolerance;
  return 1;
}

static void reset_report(MVMCClassicPfaffianAuditReport *report,
                         MVMCClassicPfaffianAuditMode requested) {
  if (report == NULL) return;
  memset(report, 0, sizeof(*report));
  report->status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  report->requested_mode = requested;
  report->effective_mode = MVMC_CLASSIC_AUDIT_OFF;
  report->failure_rank = -1;
  report->failure_local_qp = -1;
  report->failure_global_qp = -1;
}

static int close_value(double difference, double reference,
                       double observed,
                       const MVMCClassicPfaffianAuditLayout *layout) {
  const double scale = fmax(fabs(reference), fabs(observed));
  return isfinite(difference) &&
         difference <= layout->absolute_tolerance +
                           layout->relative_tolerance * scale;
}

static int close_complex(double complex reference, double complex observed,
                         const MVMCClassicPfaffianAuditLayout *layout) {
  return close_value(cabs(reference - observed), cabs(reference),
                     cabs(observed), layout);
}

static int valid_state(const MVMCClassicPfaffianState *state,
                       size_t local_qp_count) {
  return state != NULL && state->valid == 1 &&
         state->local_aggregate.valid == 1 &&
         (local_qp_count == 0 || state->components != NULL);
}

static MVMCPfaffianStatus project_observation(
    const MVMCClassicPfaffianAuditLayout *layout,
    const MVMCAbsolutePfaffianResult *components,
    const double complex *global_weights,
    MVMCProjectedAmplitudeResult *aggregate) {
  if (layout->local_qp_count == 0) {
    memset(aggregate, 0, sizeof(*aggregate));
    aggregate->valid = 1;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  return mvmc_projected_amplitude_slice(
      components, layout->local_qp_count, global_weights,
      layout->qp_total, layout->qp_start, layout->qp_end, aggregate);
}

static int locally_eligible(
    const MVMCClassicPfaffianAuditLayout *layout,
    const MVMCClassicPfaffianState *accepted_state) {
  size_t local_qp;

  if (!valid_state(accepted_state, layout->local_qp_count)) return 0;
  for (local_qp = 0; local_qp < layout->local_qp_count; ++local_qp) {
    const MVMCAbsolutePfaffianResult *component =
        &accepted_state->components[local_qp];
    if (component->state != MVMC_PFAFFIAN_REGULAR ||
        component->inverse_valid != 1 ||
        !isfinite(component->scaled_min_pivot) ||
        !isfinite(component->inverse_residual)) {
      return 0;
    }
  }
  return 1;
}

static int close_aggregate(
    const MVMCProjectedAmplitudeResult *reference,
    const MVMCProjectedAmplitudeResult *observed,
    const MVMCClassicPfaffianAuditLayout *layout) {
  return reference != NULL && observed != NULL &&
         reference->valid == 1 && observed->valid == 1 &&
         reference->regular_count == observed->regular_count &&
         reference->near_singular_count == observed->near_singular_count &&
         reference->singular_count == observed->singular_count &&
         close_complex(reference->total, observed->total, layout) &&
         close_value(fabs(reference->sum_abs - observed->sum_abs),
                     reference->sum_abs, observed->sum_abs, layout) &&
         close_value(fabs(reference->cancellation_ratio -
                              observed->cancellation_ratio),
                     reference->cancellation_ratio,
                     observed->cancellation_ratio, layout);
}

static int compare_real_components(
    const MVMCClassicPfaffianAuditLayout *layout,
    const MVMCClassicPfaffianState *reference_state,
    const double *reference_inverse,
    const MVMCClassicPfaffianState *observed_state,
    const double *observed_inverse,
    int *failure_local_qp, int *failure_global_qp) {
  size_t local_qp, index;

  if (!valid_state(reference_state, layout->local_qp_count) ||
      !valid_state(observed_state, layout->local_qp_count) ||
      (layout->inverse_count != 0 &&
       (reference_inverse == NULL || observed_inverse == NULL))) {
    return 0;
  }
  for (local_qp = 0; local_qp < layout->local_qp_count; ++local_qp) {
    const MVMCAbsolutePfaffianResult *reference =
        &reference_state->components[local_qp];
    const MVMCAbsolutePfaffianResult *observed =
        &observed_state->components[local_qp];
    int match = reference->state == observed->state &&
                reference->inverse_valid == observed->inverse_valid &&
                close_complex(reference->pfaffian, observed->pfaffian,
                              layout);
    for (index = 0; match && index < layout->matrix_count; ++index) {
      const size_t offset = local_qp * layout->matrix_count + index;
      if (reference->state == MVMC_PFAFFIAN_REGULAR) {
        match = close_value(fabs(reference_inverse[offset] -
                                     observed_inverse[offset]),
                            reference_inverse[offset],
                            observed_inverse[offset], layout);
      } else {
        match = reference_inverse[offset] == 0.0 &&
                observed_inverse[offset] == 0.0;
      }
    }
    if (!match) {
      *failure_local_qp = (int)local_qp;
      *failure_global_qp = layout->qp_start + (int)local_qp;
      return 0;
    }
  }
  return 1;
}

static int compare_complex_components(
    const MVMCClassicPfaffianAuditLayout *layout,
    const MVMCClassicPfaffianState *reference_state,
    const double complex *reference_inverse,
    const MVMCClassicPfaffianState *observed_state,
    const double complex *observed_inverse,
    int *failure_local_qp, int *failure_global_qp) {
  size_t local_qp, index;

  if (!valid_state(reference_state, layout->local_qp_count) ||
      !valid_state(observed_state, layout->local_qp_count) ||
      (layout->inverse_count != 0 &&
       (reference_inverse == NULL || observed_inverse == NULL))) {
    return 0;
  }
  for (local_qp = 0; local_qp < layout->local_qp_count; ++local_qp) {
    const MVMCAbsolutePfaffianResult *reference =
        &reference_state->components[local_qp];
    const MVMCAbsolutePfaffianResult *observed =
        &observed_state->components[local_qp];
    int match = reference->state == observed->state &&
                reference->inverse_valid == observed->inverse_valid &&
                close_complex(reference->pfaffian, observed->pfaffian,
                              layout);
    for (index = 0; match && index < layout->matrix_count; ++index) {
      const size_t offset = local_qp * layout->matrix_count + index;
      if (reference->state == MVMC_PFAFFIAN_REGULAR) {
        match = close_complex(reference_inverse[offset],
                              observed_inverse[offset], layout);
      } else {
        match = reference_inverse[offset] == 0.0 &&
                observed_inverse[offset] == 0.0;
      }
    }
    if (!match) {
      *failure_local_qp = (int)local_qp;
      *failure_global_qp = layout->qp_start + (int)local_qp;
      return 0;
    }
  }
  return 1;
}

static MVMCPfaffianStatus finish_collective_compare(
    const MVMCClassicPfaffianAuditLayout *layout,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    MVMCPfaffianStatus local_status,
    const MVMCClassicPfaffianState *observed_state,
    int failure_local_qp, int failure_global_qp,
    const MVMCProjectedAmplitudeResult *reference_global_aggregate,
    MVMCClassicPfaffianAuditReport *report) {
  MVMCClassicPfaffianCollectiveResult collected;
  MVMCPfaffianStatus status = mvmc_classic_pfaffian_collective_aggregate(
      collective_workspace, local_status,
      observed_state == NULL ? NULL : &observed_state->local_aggregate,
      layout->qp_total, layout->qp_start, layout->qp_end,
      failure_local_qp, failure_global_qp, &collected);

  report->executed = 1;
  report->effective_mode = report->requested_mode == MVMC_CLASSIC_AUDIT_OFF
                               ? MVMC_CLASSIC_AUDIT_ALWAYS
                               : report->requested_mode;
  report->valid = 1;
  report->status = status;
  if (status != MVMC_PFAFFIAN_STATUS_OK) {
    report->mismatch = 1;
    report->failure_rank = collected.failure_rank;
    report->failure_local_qp = collected.failure_local_qp;
    report->failure_global_qp = collected.failure_global_qp;
    return status;
  }
  report->observed_aggregate = collected.aggregate;
  if (!close_aggregate(reference_global_aggregate, &collected.aggregate,
                       layout)) {
    report->mismatch = 1;
    report->status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    return report->status;
  }
  report->status = MVMC_PFAFFIAN_STATUS_OK;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_real_audit_workspace_create(
    const MVMCClassicPfaffianRealWorkspace *state_workspace,
    MVMCClassicPfaffianAuditMode mode,
    double absolute_tolerance, double relative_tolerance,
    MVMCClassicPfaffianRealAuditWorkspace **workspace) {
  MVMCClassicPfaffianRealAuditWorkspace *created = NULL;
  int nsize, qp_total, qp_start, qp_end;

  if (workspace == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (!mvmc_classic_pfaffian_real_layout(
          state_workspace, &nsize, &qp_total, &qp_start, &qp_end)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  created = (MVMCClassicPfaffianRealAuditWorkspace *)calloc(
      1, sizeof(*created));
  if (created == NULL) return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  if (!initialize_layout(nsize, qp_total, qp_start, qp_end, mode,
                         absolute_tolerance, relative_tolerance,
                         &created->layout)) {
    free(created);
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  if (created->layout.local_qp_count != 0) {
    created->layout.components = (MVMCAbsolutePfaffianResult *)calloc(
        created->layout.local_qp_count, sizeof(*created->layout.components));
    created->inverse = (double *)calloc(created->layout.inverse_count,
                                        sizeof(*created->inverse));
  }
  if (created->layout.local_qp_count != 0 &&
      (created->layout.components == NULL || created->inverse == NULL)) {
    mvmc_classic_pfaffian_real_audit_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  }
  *workspace = created;
  return MVMC_PFAFFIAN_STATUS_OK;
}

void mvmc_classic_pfaffian_real_audit_workspace_destroy(
    MVMCClassicPfaffianRealAuditWorkspace *workspace) {
  if (workspace == NULL) return;
  free(workspace->inverse);
  free(workspace->layout.components);
  free(workspace);
}

MVMCPfaffianStatus mvmc_classic_pfaffian_complex_audit_workspace_create(
    const MVMCClassicPfaffianComplexWorkspace *state_workspace,
    MVMCClassicPfaffianAuditMode mode,
    double absolute_tolerance, double relative_tolerance,
    MVMCClassicPfaffianComplexAuditWorkspace **workspace) {
  MVMCClassicPfaffianComplexAuditWorkspace *created = NULL;
  int nsize, qp_total, qp_start, qp_end;

  if (workspace == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (!mvmc_classic_pfaffian_complex_layout(
          state_workspace, &nsize, &qp_total, &qp_start, &qp_end)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  created = (MVMCClassicPfaffianComplexAuditWorkspace *)calloc(
      1, sizeof(*created));
  if (created == NULL) return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  if (!initialize_layout(nsize, qp_total, qp_start, qp_end, mode,
                         absolute_tolerance, relative_tolerance,
                         &created->layout)) {
    free(created);
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  if (created->layout.local_qp_count != 0) {
    created->layout.components = (MVMCAbsolutePfaffianResult *)calloc(
        created->layout.local_qp_count, sizeof(*created->layout.components));
    created->inverse = (double complex *)calloc(
        created->layout.inverse_count, sizeof(*created->inverse));
  }
  if (created->layout.local_qp_count != 0 &&
      (created->layout.components == NULL || created->inverse == NULL)) {
    mvmc_classic_pfaffian_complex_audit_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  }
  *workspace = created;
  return MVMC_PFAFFIAN_STATUS_OK;
}

void mvmc_classic_pfaffian_complex_audit_workspace_destroy(
    MVMCClassicPfaffianComplexAuditWorkspace *workspace) {
  if (workspace == NULL) return;
  free(workspace->inverse);
  free(workspace->layout.components);
  free(workspace);
}

MVMCClassicPfaffianAuditMode mvmc_classic_pfaffian_real_audit_mode(
    const MVMCClassicPfaffianRealAuditWorkspace *workspace) {
  return workspace == NULL ? MVMC_CLASSIC_AUDIT_OFF : workspace->layout.mode;
}

MVMCClassicPfaffianAuditMode mvmc_classic_pfaffian_complex_audit_mode(
    const MVMCClassicPfaffianComplexAuditWorkspace *workspace) {
  return workspace == NULL ? MVMC_CLASSIC_AUDIT_OFF : workspace->layout.mode;
}

static int compatible_layout(const MVMCClassicPfaffianAuditLayout *audit,
                             int nsize, int qp_total,
                             int qp_start, int qp_end) {
  return audit != NULL && audit->nsize == nsize &&
         audit->qp_total == qp_total && audit->qp_start == qp_start &&
         audit->qp_end == qp_end;
}

int mvmc_classic_pfaffian_real_audit_compatible(
    const MVMCClassicPfaffianRealAuditWorkspace *audit_workspace,
    const MVMCClassicPfaffianRealWorkspace *state_workspace) {
  int nsize, qp_total, qp_start, qp_end;
  return audit_workspace != NULL &&
         mvmc_classic_pfaffian_real_layout(
             state_workspace, &nsize, &qp_total, &qp_start, &qp_end) &&
         compatible_layout(&audit_workspace->layout, nsize, qp_total,
                           qp_start, qp_end);
}

int mvmc_classic_pfaffian_complex_audit_compatible(
    const MVMCClassicPfaffianComplexAuditWorkspace *audit_workspace,
    const MVMCClassicPfaffianComplexWorkspace *state_workspace) {
  int nsize, qp_total, qp_start, qp_end;
  return audit_workspace != NULL &&
         mvmc_classic_pfaffian_complex_layout(
             state_workspace, &nsize, &qp_total, &qp_start, &qp_end) &&
         compatible_layout(&audit_workspace->layout, nsize, qp_total,
                           qp_start, qp_end);
}

MVMCPfaffianStatus mvmc_classic_pfaffian_real_audit_compare(
    MVMCClassicPfaffianRealAuditWorkspace *audit_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const MVMCClassicPfaffianState *reference_state,
    const double *reference_inverse,
    const MVMCProjectedAmplitudeResult *reference_global_aggregate,
    const MVMCClassicPfaffianState *observed_state,
    const double *observed_inverse,
    MVMCClassicPfaffianAuditReport *report) {
  MVMCPfaffianStatus local_status = MVMC_PFAFFIAN_STATUS_OK;
  int failure_local_qp = -1, failure_global_qp = -1;
  MVMCClassicPfaffianAuditMode mode =
      mvmc_classic_pfaffian_real_audit_mode(audit_workspace);

  reset_report(report, mode);
  if (audit_workspace == NULL || collective_workspace == NULL ||
      report == NULL || reference_global_aggregate == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  if (!compare_real_components(
          &audit_workspace->layout, reference_state, reference_inverse,
          observed_state, observed_inverse, &failure_local_qp,
          &failure_global_qp)) {
    local_status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  return finish_collective_compare(
      &audit_workspace->layout, collective_workspace, local_status,
      observed_state, failure_local_qp, failure_global_qp,
      reference_global_aggregate, report);
}

MVMCPfaffianStatus mvmc_classic_pfaffian_complex_audit_compare(
    MVMCClassicPfaffianComplexAuditWorkspace *audit_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const MVMCClassicPfaffianState *reference_state,
    const double complex *reference_inverse,
    const MVMCProjectedAmplitudeResult *reference_global_aggregate,
    const MVMCClassicPfaffianState *observed_state,
    const double complex *observed_inverse,
    MVMCClassicPfaffianAuditReport *report) {
  MVMCPfaffianStatus local_status = MVMC_PFAFFIAN_STATUS_OK;
  int failure_local_qp = -1, failure_global_qp = -1;
  MVMCClassicPfaffianAuditMode mode =
      mvmc_classic_pfaffian_complex_audit_mode(audit_workspace);

  reset_report(report, mode);
  if (audit_workspace == NULL || collective_workspace == NULL ||
      report == NULL || reference_global_aggregate == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  if (!compare_complex_components(
          &audit_workspace->layout, reference_state, reference_inverse,
          observed_state, observed_inverse, &failure_local_qp,
          &failure_global_qp)) {
    local_status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  return finish_collective_compare(
      &audit_workspace->layout, collective_workspace, local_status,
      observed_state, failure_local_qp, failure_global_qp,
      reference_global_aggregate, report);
}

static MVMCPfaffianStatus select_audit_mode(
    MVMCClassicPfaffianAuditLayout *layout,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const MVMCClassicPfaffianState *accepted_state,
    MVMCClassicPfaffianAuditReport *report) {
  int all_eligible = 1;
  MVMCPfaffianStatus status;

  if (layout->mode == MVMC_CLASSIC_AUDIT_OFF) {
    report->valid = 1;
    report->status = MVMC_PFAFFIAN_STATUS_OK;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  if (layout->mode == MVMC_CLASSIC_AUDIT_GUARDED) {
    status = mvmc_classic_pfaffian_collective_all_true(
        collective_workspace, locally_eligible(layout, accepted_state),
        &all_eligible);
    if (status != MVMC_PFAFFIAN_STATUS_OK) return status;
    if (!all_eligible) {
      report->valid = 1;
      report->status = MVMC_PFAFFIAN_STATUS_OK;
      report->fallback = 1;
      return MVMC_PFAFFIAN_STATUS_OK;
    }
  }
  report->effective_mode = layout->mode;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_real_audit_candidate(
    MVMCClassicPfaffianRealAuditWorkspace *audit_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const MVMCClassicPfaffianState *accepted_state,
    const MVMCClassicPfaffianState *absolute_candidate,
    const double *absolute_candidate_inverse,
    const MVMCProjectedAmplitudeResult *absolute_global_aggregate,
    const double complex *global_weights,
    MVMCClassicPfaffianRealObserver observer, void *observer_context,
    MVMCClassicPfaffianAuditReport *report) {
  MVMCClassicPfaffianRealObservation observation;
  MVMCClassicPfaffianState observed_state;
  MVMCPfaffianStatus local_status;
  MVMCPfaffianStatus status;

  reset_report(report, mvmc_classic_pfaffian_real_audit_mode(audit_workspace));
  if (audit_workspace == NULL || collective_workspace == NULL ||
      report == NULL || absolute_candidate == NULL ||
      absolute_global_aggregate == NULL || global_weights == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  status = select_audit_mode(&audit_workspace->layout, collective_workspace,
                             accepted_state, report);
  if (status != MVMC_PFAFFIAN_STATUS_OK || report->fallback ||
      audit_workspace->layout.mode == MVMC_CLASSIC_AUDIT_OFF) {
    return status;
  }
  if (observer == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  if (audit_workspace->layout.local_qp_count != 0) {
    memset(audit_workspace->layout.components, 0,
           audit_workspace->layout.local_qp_count *
               sizeof(*audit_workspace->layout.components));
  }
  if (audit_workspace->layout.inverse_count != 0) {
    memset(audit_workspace->inverse, 0,
           audit_workspace->layout.inverse_count *
               sizeof(*audit_workspace->inverse));
  }
  memset(&observation, 0, sizeof(observation));
  observation.component_count = audit_workspace->layout.local_qp_count;
  observation.matrix_element_count = audit_workspace->layout.matrix_count;
  observation.components = audit_workspace->layout.components;
  observation.inverse = audit_workspace->inverse;
  local_status = observer(observer_context, absolute_candidate,
                          absolute_candidate_inverse, &observation);
  memset(&observed_state, 0, sizeof(observed_state));
  observed_state.components = audit_workspace->layout.components;
  if (local_status == MVMC_PFAFFIAN_STATUS_OK &&
      observation.valid == 1 &&
      observation.components == audit_workspace->layout.components &&
      observation.inverse == audit_workspace->inverse &&
      observation.component_count == audit_workspace->layout.local_qp_count &&
      observation.matrix_element_count == audit_workspace->layout.matrix_count) {
    observed_state.valid = 1;
    observed_state.components = observation.components;
    local_status = project_observation(
        &audit_workspace->layout, observation.components, global_weights,
        &observed_state.local_aggregate);
  } else if (local_status == MVMC_PFAFFIAN_STATUS_OK) {
    local_status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  if (local_status != MVMC_PFAFFIAN_STATUS_OK) {
    observed_state.valid = 0;
  }
  {
    int failure_local_qp = -1, failure_global_qp = -1;
    if (local_status == MVMC_PFAFFIAN_STATUS_OK &&
        !compare_real_components(
            &audit_workspace->layout, absolute_candidate,
            absolute_candidate_inverse, &observed_state,
            audit_workspace->inverse, &failure_local_qp,
            &failure_global_qp)) {
      local_status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    return finish_collective_compare(
        &audit_workspace->layout, collective_workspace, local_status,
        &observed_state, failure_local_qp, failure_global_qp,
        absolute_global_aggregate, report);
  }
}

MVMCPfaffianStatus mvmc_classic_pfaffian_complex_audit_candidate(
    MVMCClassicPfaffianComplexAuditWorkspace *audit_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const MVMCClassicPfaffianState *accepted_state,
    const MVMCClassicPfaffianState *absolute_candidate,
    const double complex *absolute_candidate_inverse,
    const MVMCProjectedAmplitudeResult *absolute_global_aggregate,
    const double complex *global_weights,
    MVMCClassicPfaffianComplexObserver observer, void *observer_context,
    MVMCClassicPfaffianAuditReport *report) {
  MVMCClassicPfaffianComplexObservation observation;
  MVMCClassicPfaffianState observed_state;
  MVMCPfaffianStatus local_status;
  MVMCPfaffianStatus status;

  reset_report(report,
               mvmc_classic_pfaffian_complex_audit_mode(audit_workspace));
  if (audit_workspace == NULL || collective_workspace == NULL ||
      report == NULL || absolute_candidate == NULL ||
      absolute_global_aggregate == NULL || global_weights == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  status = select_audit_mode(&audit_workspace->layout, collective_workspace,
                             accepted_state, report);
  if (status != MVMC_PFAFFIAN_STATUS_OK || report->fallback ||
      audit_workspace->layout.mode == MVMC_CLASSIC_AUDIT_OFF) {
    return status;
  }
  if (observer == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  if (audit_workspace->layout.local_qp_count != 0) {
    memset(audit_workspace->layout.components, 0,
           audit_workspace->layout.local_qp_count *
               sizeof(*audit_workspace->layout.components));
  }
  if (audit_workspace->layout.inverse_count != 0) {
    memset(audit_workspace->inverse, 0,
           audit_workspace->layout.inverse_count *
               sizeof(*audit_workspace->inverse));
  }
  memset(&observation, 0, sizeof(observation));
  observation.component_count = audit_workspace->layout.local_qp_count;
  observation.matrix_element_count = audit_workspace->layout.matrix_count;
  observation.components = audit_workspace->layout.components;
  observation.inverse = audit_workspace->inverse;
  local_status = observer(observer_context, absolute_candidate,
                          absolute_candidate_inverse, &observation);
  memset(&observed_state, 0, sizeof(observed_state));
  observed_state.components = audit_workspace->layout.components;
  if (local_status == MVMC_PFAFFIAN_STATUS_OK &&
      observation.valid == 1 &&
      observation.components == audit_workspace->layout.components &&
      observation.inverse == audit_workspace->inverse &&
      observation.component_count == audit_workspace->layout.local_qp_count &&
      observation.matrix_element_count == audit_workspace->layout.matrix_count) {
    observed_state.valid = 1;
    observed_state.components = observation.components;
    local_status = project_observation(
        &audit_workspace->layout, observation.components, global_weights,
        &observed_state.local_aggregate);
  } else if (local_status == MVMC_PFAFFIAN_STATUS_OK) {
    local_status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  if (local_status != MVMC_PFAFFIAN_STATUS_OK) {
    observed_state.valid = 0;
  }
  {
    int failure_local_qp = -1, failure_global_qp = -1;
    if (local_status == MVMC_PFAFFIAN_STATUS_OK &&
        !compare_complex_components(
            &audit_workspace->layout, absolute_candidate,
            absolute_candidate_inverse, &observed_state,
            audit_workspace->inverse, &failure_local_qp,
            &failure_global_qp)) {
      local_status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    return finish_collective_compare(
        &audit_workspace->layout, collective_workspace, local_status,
        &observed_state, failure_local_qp, failure_global_qp,
        absolute_global_aggregate, report);
  }
}
