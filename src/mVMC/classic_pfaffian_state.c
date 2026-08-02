/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "classic_pfaffian_state.h"
#include "classic_pfaffian_matrix.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int nsite;
  int nsite2;
  int nsize;
  int qp_total;
  int qp_start;
  int qp_end;
  size_t local_qp_count;
  size_t matrix_count;
  size_t slater_matrix_count;
  size_t inverse_count;
  size_t legacy_element_count;
  int candidate_ready;
  int failure_local_qp;
  int failure_global_qp;
} MVMCClassicPfaffianLayout;

struct MVMCClassicPfaffianRealWorkspace {
  MVMCClassicPfaffianLayout layout;
  MVMCClassicPfaffianState accepted;
  MVMCClassicPfaffianState candidate;
  double *accepted_inverse;
  double *candidate_inverse;
  double *matrix_scratch;
  MVMCAbsolutePfaffianRealWorkspace *absolute_workspace;
};

struct MVMCClassicPfaffianComplexWorkspace {
  MVMCClassicPfaffianLayout layout;
  MVMCClassicPfaffianState accepted;
  MVMCClassicPfaffianState candidate;
  double complex *accepted_inverse;
  double complex *candidate_inverse;
  double complex *matrix_scratch;
  MVMCAbsolutePfaffianComplexWorkspace *absolute_workspace;
};

static int checked_multiply(size_t left, size_t right, size_t *product) {
  if (product == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *product = left * right;
  return 1;
}

static int initialize_layout(int nsite, int nsize, int qp_total,
                             int qp_start, int qp_end,
                             size_t element_size,
                             MVMCClassicPfaffianLayout *layout) {
  size_t nsite2;
  size_t total_matrix_elements;
  size_t legacy_per_qp;

  if (layout == NULL) return 0;
  memset(layout, 0, sizeof(*layout));
  layout->failure_local_qp = -1;
  layout->failure_global_qp = -1;
  if (nsite <= 0 || nsite > INT_MAX / 2 || nsize <= 0 ||
      (nsize % 2) != 0 || nsize / 2 > nsite || qp_total <= 0 ||
      qp_start < 0 || qp_end < qp_start || qp_end > qp_total) {
    return 0;
  }
  nsite2 = (size_t)nsite * 2;
  if (!checked_multiply((size_t)nsize, (size_t)nsize,
                        &layout->matrix_count) ||
      layout->matrix_count > SIZE_MAX / element_size ||
      !checked_multiply(nsite2, nsite2, &layout->slater_matrix_count) ||
      layout->slater_matrix_count > SIZE_MAX / element_size ||
      !checked_multiply((size_t)qp_total, layout->slater_matrix_count,
                        &total_matrix_elements) ||
      total_matrix_elements > SIZE_MAX / element_size ||
      !checked_multiply((size_t)(qp_end - qp_start),
                        layout->matrix_count, &layout->inverse_count) ||
      layout->inverse_count > SIZE_MAX / element_size ||
      layout->matrix_count == SIZE_MAX ||
      !checked_multiply((size_t)qp_total, layout->matrix_count + 1,
                        &layout->legacy_element_count) ||
      layout->legacy_element_count > SIZE_MAX / element_size) {
    return 0;
  }
  legacy_per_qp = layout->matrix_count + 1;
  if (legacy_per_qp == 0) return 0;
  layout->nsite = nsite;
  layout->nsite2 = (int)nsite2;
  layout->nsize = nsize;
  layout->qp_total = qp_total;
  layout->qp_start = qp_start;
  layout->qp_end = qp_end;
  layout->local_qp_count = (size_t)(qp_end - qp_start);
  return 1;
}

static void reset_state(MVMCClassicPfaffianState *state,
                        size_t component_count, uint64_t generation) {
  size_t index;
  MVMCAbsolutePfaffianResult *components;

  if (state == NULL) return;
  components = (MVMCAbsolutePfaffianResult *)state->components;
  memset(&state->local_aggregate, 0, sizeof(state->local_aggregate));
  state->valid = 0;
  state->generation = generation;
  for (index = 0; index < component_count; ++index) {
    memset(&components[index], 0, sizeof(components[index]));
    components[index].state = MVMC_PFAFFIAN_INVALID;
    components[index].rebuild_generation = generation;
  }
}

static MVMCAbsolutePfaffianResult *mutable_components(
    MVMCClassicPfaffianState *state) {
  return state == NULL ? NULL
                       : (MVMCAbsolutePfaffianResult *)state->components;
}

static int validate_ele_idx(const MVMCClassicPfaffianLayout *layout,
                            const int *ele_idx) {
  int index;

  if (layout == NULL || ele_idx == NULL) return 0;
  for (index = 0; index < layout->nsize; ++index) {
    if (ele_idx[index] < 0 || ele_idx[index] >= layout->nsite) return 0;
  }
  return 1;
}

static int validate_component(const MVMCAbsolutePfaffianResult *component) {
  if (component == NULL ||
      !isfinite(creal(component->pfaffian)) ||
      !isfinite(cimag(component->pfaffian))) {
    return 0;
  }
  switch (component->state) {
    case MVMC_PFAFFIAN_REGULAR:
      return component->inverse_valid == 1;
    case MVMC_PFAFFIAN_NEAR_SINGULAR:
      return component->inverse_valid == 0;
    case MVMC_PFAFFIAN_SINGULAR:
      return component->inverse_valid == 0 && component->pfaffian == 0.0;
    default:
      return 0;
  }
}

static int validate_aggregate(const MVMCProjectedAmplitudeResult *aggregate,
                              size_t component_count) {
  size_t counted;

  if (aggregate == NULL || !aggregate->valid ||
      !isfinite(creal(aggregate->total)) ||
      !isfinite(cimag(aggregate->total)) ||
      !isfinite(aggregate->sum_abs) || aggregate->sum_abs < 0.0 ||
      !isfinite(aggregate->cancellation_ratio) ||
      aggregate->cancellation_ratio < 0.0 ||
      aggregate->cancellation_ratio > 1.0 ||
      aggregate->regular_count > component_count ||
      aggregate->near_singular_count > component_count ||
      aggregate->singular_count > component_count) {
    return 0;
  }
  counted = aggregate->regular_count + aggregate->near_singular_count;
  if (counted < aggregate->regular_count ||
      aggregate->singular_count > SIZE_MAX - counted) {
    return 0;
  }
  counted += aggregate->singular_count;
  return counted == component_count;
}

static void set_failure(MVMCClassicPfaffianLayout *layout,
                        int local_qp, int global_qp) {
  layout->failure_local_qp = local_qp;
  layout->failure_global_qp = global_qp;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_real_workspace_create(
    int nsite, int nsize, int qp_total, int qp_start, int qp_end,
    MVMCClassicPfaffianRealWorkspace **workspace) {
  MVMCClassicPfaffianRealWorkspace *created;
  MVMCClassicPfaffianLayout layout;
  MVMCPfaffianStatus status;

  if (workspace == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (!initialize_layout(nsite, nsize, qp_total, qp_start, qp_end,
                         sizeof(double), &layout)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  created = (MVMCClassicPfaffianRealWorkspace *)calloc(1, sizeof(*created));
  if (created == NULL) return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  created->layout = layout;
  created->matrix_scratch =
      (double *)malloc(created->layout.matrix_count *
                       sizeof(*created->matrix_scratch));
  if (created->layout.local_qp_count != 0) {
    created->accepted.components = (MVMCAbsolutePfaffianResult *)calloc(
        created->layout.local_qp_count, sizeof(*created->accepted.components));
    created->candidate.components = (MVMCAbsolutePfaffianResult *)calloc(
        created->layout.local_qp_count, sizeof(*created->candidate.components));
    created->accepted_inverse = (double *)calloc(
        created->layout.inverse_count, sizeof(*created->accepted_inverse));
    created->candidate_inverse = (double *)calloc(
        created->layout.inverse_count, sizeof(*created->candidate_inverse));
  }
  if (created->matrix_scratch == NULL ||
      (created->layout.local_qp_count != 0 &&
       (created->accepted.components == NULL ||
        created->candidate.components == NULL ||
        created->accepted_inverse == NULL ||
        created->candidate_inverse == NULL))) {
    mvmc_classic_pfaffian_real_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  }
  status = mvmc_absolute_pfaffian_real_workspace_create(
      nsize, &created->absolute_workspace);
  if (status != MVMC_PFAFFIAN_STATUS_OK) {
    mvmc_classic_pfaffian_real_workspace_destroy(created);
    return status;
  }
  reset_state(&created->accepted, created->layout.local_qp_count, 0);
  reset_state(&created->candidate, created->layout.local_qp_count, 0);
  *workspace = created;
  return MVMC_PFAFFIAN_STATUS_OK;
}

void mvmc_classic_pfaffian_real_workspace_destroy(
    MVMCClassicPfaffianRealWorkspace *workspace) {
  if (workspace == NULL) return;
  mvmc_absolute_pfaffian_real_workspace_destroy(workspace->absolute_workspace);
  free(workspace->matrix_scratch);
  free(workspace->candidate_inverse);
  free(workspace->accepted_inverse);
  free((void *)workspace->candidate.components);
  free((void *)workspace->accepted.components);
  free(workspace);
}

MVMCPfaffianStatus mvmc_classic_pfaffian_complex_workspace_create(
    int nsite, int nsize, int qp_total, int qp_start, int qp_end,
    MVMCClassicPfaffianComplexWorkspace **workspace) {
  MVMCClassicPfaffianComplexWorkspace *created;
  MVMCClassicPfaffianLayout layout;
  MVMCPfaffianStatus status;

  if (workspace == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (!initialize_layout(nsite, nsize, qp_total, qp_start, qp_end,
                         sizeof(double complex), &layout)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  created =
      (MVMCClassicPfaffianComplexWorkspace *)calloc(1, sizeof(*created));
  if (created == NULL) return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  created->layout = layout;
  created->matrix_scratch = (double complex *)malloc(
      created->layout.matrix_count * sizeof(*created->matrix_scratch));
  if (created->layout.local_qp_count != 0) {
    created->accepted.components = (MVMCAbsolutePfaffianResult *)calloc(
        created->layout.local_qp_count, sizeof(*created->accepted.components));
    created->candidate.components = (MVMCAbsolutePfaffianResult *)calloc(
        created->layout.local_qp_count, sizeof(*created->candidate.components));
    created->accepted_inverse = (double complex *)calloc(
        created->layout.inverse_count, sizeof(*created->accepted_inverse));
    created->candidate_inverse = (double complex *)calloc(
        created->layout.inverse_count, sizeof(*created->candidate_inverse));
  }
  if (created->matrix_scratch == NULL ||
      (created->layout.local_qp_count != 0 &&
       (created->accepted.components == NULL ||
        created->candidate.components == NULL ||
        created->accepted_inverse == NULL ||
        created->candidate_inverse == NULL))) {
    mvmc_classic_pfaffian_complex_workspace_destroy(created);
    return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  }
  status = mvmc_absolute_pfaffian_complex_workspace_create(
      nsize, &created->absolute_workspace);
  if (status != MVMC_PFAFFIAN_STATUS_OK) {
    mvmc_classic_pfaffian_complex_workspace_destroy(created);
    return status;
  }
  reset_state(&created->accepted, created->layout.local_qp_count, 0);
  reset_state(&created->candidate, created->layout.local_qp_count, 0);
  *workspace = created;
  return MVMC_PFAFFIAN_STATUS_OK;
}

void mvmc_classic_pfaffian_complex_workspace_destroy(
    MVMCClassicPfaffianComplexWorkspace *workspace) {
  if (workspace == NULL) return;
  mvmc_absolute_pfaffian_complex_workspace_destroy(
      workspace->absolute_workspace);
  free(workspace->matrix_scratch);
  free(workspace->candidate_inverse);
  free(workspace->accepted_inverse);
  free((void *)workspace->candidate.components);
  free((void *)workspace->accepted.components);
  free(workspace);
}

MVMCPfaffianStatus mvmc_classic_pfaffian_real_prepare(
    MVMCClassicPfaffianRealWorkspace *workspace,
    const double *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance) {
  MVMCClassicPfaffianLayout *layout;
  uint64_t generation;
  size_t local_qp;
  int global_qp;

  if (workspace == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  layout = &workspace->layout;
  layout->candidate_ready = 0;
  set_failure(layout, -1, -1);
  if (slater_elm == NULL || global_weights == NULL ||
      !validate_ele_idx(layout, ele_idx) ||
      !isfinite(scaled_pivot_tolerance) ||
      !isfinite(residual_tolerance) ||
      workspace->accepted.generation == UINT64_MAX) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  for (global_qp = 0; global_qp < layout->qp_total; ++global_qp) {
    if (!isfinite(creal(global_weights[global_qp])) ||
        !isfinite(cimag(global_weights[global_qp])) ||
        cimag(global_weights[global_qp]) != 0.0) {
      set_failure(layout,
                  global_qp >= layout->qp_start && global_qp < layout->qp_end
                      ? global_qp - layout->qp_start
                      : -1,
                  global_qp);
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
  }
  generation = workspace->accepted.generation + 1;
  reset_state(&workspace->candidate, layout->local_qp_count, generation);
  if (layout->inverse_count != 0) {
    memset(workspace->candidate_inverse, 0,
           layout->inverse_count * sizeof(*workspace->candidate_inverse));
  }
  for (local_qp = 0; local_qp < layout->local_qp_count; ++local_qp) {
    MVMCPfaffianStatus status;
    global_qp = layout->qp_start + (int)local_qp;
    status = mvmc_classic_pfaffian_build_real_matrix(
        slater_elm + (size_t)global_qp * layout->slater_matrix_count,
        layout->nsite, layout->nsize, ele_idx, workspace->matrix_scratch);
    if (status != MVMC_PFAFFIAN_STATUS_OK) {
      set_failure(layout, (int)local_qp, global_qp);
      return status;
    }
    status = mvmc_absolute_pfaffian_real_with_workspace(
        workspace->absolute_workspace, workspace->matrix_scratch,
        layout->nsize, layout->nsize,
        workspace->candidate_inverse + local_qp * layout->matrix_count,
        layout->nsize, generation, scaled_pivot_tolerance,
        residual_tolerance,
        &mutable_components(&workspace->candidate)[local_qp]);
    if (status != MVMC_PFAFFIAN_STATUS_OK ||
        !validate_component(&workspace->candidate.components[local_qp]) ||
        cimag(workspace->candidate.components[local_qp].pfaffian) != 0.0) {
      set_failure(layout, (int)local_qp, global_qp);
      return status == MVMC_PFAFFIAN_STATUS_OK
                 ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT
                 : status;
    }
    {
      const double complex term =
          global_weights[global_qp] *
          workspace->candidate.components[local_qp].pfaffian;
      if (!isfinite(creal(term)) || !isfinite(cimag(term))) {
        set_failure(layout, (int)local_qp, global_qp);
        return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
      }
    }
  }
  if (layout->local_qp_count == 0) {
    memset(&workspace->candidate.local_aggregate, 0,
           sizeof(workspace->candidate.local_aggregate));
    workspace->candidate.local_aggregate.valid = 1;
  } else {
    MVMCPfaffianStatus status = mvmc_projected_amplitude_slice(
        workspace->candidate.components, layout->local_qp_count,
        global_weights, layout->qp_total, layout->qp_start, layout->qp_end,
        &workspace->candidate.local_aggregate);
    if (status != MVMC_PFAFFIAN_STATUS_OK) return status;
  }
  if (!validate_aggregate(&workspace->candidate.local_aggregate,
                          layout->local_qp_count) ||
      cimag(workspace->candidate.local_aggregate.total) != 0.0) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  workspace->candidate.valid = 1;
  layout->candidate_ready = 1;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_complex_prepare(
    MVMCClassicPfaffianComplexWorkspace *workspace,
    const double complex *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance) {
  MVMCClassicPfaffianLayout *layout;
  uint64_t generation;
  size_t local_qp;
  int global_qp;

  if (workspace == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  layout = &workspace->layout;
  layout->candidate_ready = 0;
  set_failure(layout, -1, -1);
  if (slater_elm == NULL || global_weights == NULL ||
      !validate_ele_idx(layout, ele_idx) ||
      !isfinite(scaled_pivot_tolerance) ||
      !isfinite(residual_tolerance) ||
      workspace->accepted.generation == UINT64_MAX) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  for (global_qp = 0; global_qp < layout->qp_total; ++global_qp) {
    if (!isfinite(creal(global_weights[global_qp])) ||
        !isfinite(cimag(global_weights[global_qp]))) {
      set_failure(layout,
                  global_qp >= layout->qp_start && global_qp < layout->qp_end
                      ? global_qp - layout->qp_start
                      : -1,
                  global_qp);
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
  }
  generation = workspace->accepted.generation + 1;
  reset_state(&workspace->candidate, layout->local_qp_count, generation);
  if (layout->inverse_count != 0) {
    memset(workspace->candidate_inverse, 0,
           layout->inverse_count * sizeof(*workspace->candidate_inverse));
  }
  for (local_qp = 0; local_qp < layout->local_qp_count; ++local_qp) {
    MVMCPfaffianStatus status;
    global_qp = layout->qp_start + (int)local_qp;
    status = mvmc_classic_pfaffian_build_complex_matrix(
        slater_elm + (size_t)global_qp * layout->slater_matrix_count,
        layout->nsite, layout->nsize, ele_idx, workspace->matrix_scratch);
    if (status != MVMC_PFAFFIAN_STATUS_OK) {
      set_failure(layout, (int)local_qp, global_qp);
      return status;
    }
    status = mvmc_absolute_pfaffian_complex_with_workspace(
        workspace->absolute_workspace, workspace->matrix_scratch,
        layout->nsize, layout->nsize,
        workspace->candidate_inverse + local_qp * layout->matrix_count,
        layout->nsize, generation, scaled_pivot_tolerance,
        residual_tolerance,
        &mutable_components(&workspace->candidate)[local_qp]);
    if (status != MVMC_PFAFFIAN_STATUS_OK ||
        !validate_component(&workspace->candidate.components[local_qp])) {
      set_failure(layout, (int)local_qp, global_qp);
      return status == MVMC_PFAFFIAN_STATUS_OK
                 ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT
                 : status;
    }
    {
      const double complex term =
          global_weights[global_qp] *
          workspace->candidate.components[local_qp].pfaffian;
      if (!isfinite(creal(term)) || !isfinite(cimag(term))) {
        set_failure(layout, (int)local_qp, global_qp);
        return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
      }
    }
  }
  if (layout->local_qp_count == 0) {
    memset(&workspace->candidate.local_aggregate, 0,
           sizeof(workspace->candidate.local_aggregate));
    workspace->candidate.local_aggregate.valid = 1;
  } else {
    MVMCPfaffianStatus status = mvmc_projected_amplitude_slice(
        workspace->candidate.components, layout->local_qp_count,
        global_weights, layout->qp_total, layout->qp_start, layout->qp_end,
        &workspace->candidate.local_aggregate);
    if (status != MVMC_PFAFFIAN_STATUS_OK) return status;
  }
  if (!validate_aggregate(&workspace->candidate.local_aggregate,
                          layout->local_qp_count)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  workspace->candidate.valid = 1;
  layout->candidate_ready = 1;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_real_publish(
    MVMCClassicPfaffianRealWorkspace *workspace,
    double *legacy_storage, size_t legacy_element_count) {
  MVMCClassicPfaffianState state_swap;
  double *inverse_swap;
  double *legacy_pfaffian;
  size_t local_qp, index;

  if (workspace == NULL || legacy_storage == NULL ||
      !workspace->layout.candidate_ready || !workspace->candidate.valid ||
      !validate_aggregate(&workspace->candidate.local_aggregate,
                          workspace->layout.local_qp_count) ||
      cimag(workspace->candidate.local_aggregate.total) != 0.0 ||
      legacy_element_count < workspace->layout.legacy_element_count) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  for (local_qp = 0; local_qp < workspace->layout.local_qp_count;
       ++local_qp) {
    const MVMCAbsolutePfaffianResult *component =
        &workspace->candidate.components[local_qp];
    if (!validate_component(component) || cimag(component->pfaffian) != 0.0 ||
        component->rebuild_generation != workspace->candidate.generation) {
      set_failure(&workspace->layout, (int)local_qp,
                  workspace->layout.qp_start + (int)local_qp);
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    if (component->state == MVMC_PFAFFIAN_REGULAR) {
      for (index = 0; index < workspace->layout.matrix_count; ++index) {
        if (!isfinite(workspace->candidate_inverse[
                local_qp * workspace->layout.matrix_count + index])) {
          set_failure(&workspace->layout, (int)local_qp,
                      workspace->layout.qp_start + (int)local_qp);
          return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
        }
      }
    }
  }
  legacy_pfaffian = legacy_storage +
                     (size_t)workspace->layout.qp_total *
                         workspace->layout.matrix_count;
  for (local_qp = 0; local_qp < workspace->layout.local_qp_count;
       ++local_qp) {
    const MVMCAbsolutePfaffianResult *component =
        &workspace->candidate.components[local_qp];
    double *legacy_inverse =
        legacy_storage + local_qp * workspace->layout.matrix_count;
    if (component->state == MVMC_PFAFFIAN_REGULAR) {
      for (index = 0; index < workspace->layout.matrix_count; ++index) {
        legacy_inverse[index] =
            -workspace->candidate_inverse[
                local_qp * workspace->layout.matrix_count + index];
      }
    } else {
      memset(legacy_inverse, 0,
             workspace->layout.matrix_count * sizeof(*legacy_inverse));
    }
    legacy_pfaffian[local_qp] = creal(component->pfaffian);
  }
  state_swap = workspace->accepted;
  workspace->accepted = workspace->candidate;
  workspace->candidate = state_swap;
  inverse_swap = workspace->accepted_inverse;
  workspace->accepted_inverse = workspace->candidate_inverse;
  workspace->candidate_inverse = inverse_swap;
  workspace->candidate.valid = 0;
  workspace->layout.candidate_ready = 0;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_complex_publish(
    MVMCClassicPfaffianComplexWorkspace *workspace,
    double complex *legacy_storage, size_t legacy_element_count) {
  MVMCClassicPfaffianState state_swap;
  double complex *inverse_swap;
  double complex *legacy_pfaffian;
  size_t local_qp, index;

  if (workspace == NULL || legacy_storage == NULL ||
      !workspace->layout.candidate_ready || !workspace->candidate.valid ||
      !validate_aggregate(&workspace->candidate.local_aggregate,
                          workspace->layout.local_qp_count) ||
      legacy_element_count < workspace->layout.legacy_element_count) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  for (local_qp = 0; local_qp < workspace->layout.local_qp_count;
       ++local_qp) {
    const MVMCAbsolutePfaffianResult *component =
        &workspace->candidate.components[local_qp];
    if (!validate_component(component) ||
        component->rebuild_generation != workspace->candidate.generation) {
      set_failure(&workspace->layout, (int)local_qp,
                  workspace->layout.qp_start + (int)local_qp);
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    if (component->state == MVMC_PFAFFIAN_REGULAR) {
      for (index = 0; index < workspace->layout.matrix_count; ++index) {
        const double complex value = workspace->candidate_inverse[
            local_qp * workspace->layout.matrix_count + index];
        if (!isfinite(creal(value)) || !isfinite(cimag(value))) {
          set_failure(&workspace->layout, (int)local_qp,
                      workspace->layout.qp_start + (int)local_qp);
          return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
        }
      }
    }
  }
  legacy_pfaffian = legacy_storage +
                     (size_t)workspace->layout.qp_total *
                         workspace->layout.matrix_count;
  for (local_qp = 0; local_qp < workspace->layout.local_qp_count;
       ++local_qp) {
    const MVMCAbsolutePfaffianResult *component =
        &workspace->candidate.components[local_qp];
    double complex *legacy_inverse =
        legacy_storage + local_qp * workspace->layout.matrix_count;
    if (component->state == MVMC_PFAFFIAN_REGULAR) {
      for (index = 0; index < workspace->layout.matrix_count; ++index) {
        legacy_inverse[index] =
            -workspace->candidate_inverse[
                local_qp * workspace->layout.matrix_count + index];
      }
    } else {
      memset(legacy_inverse, 0,
             workspace->layout.matrix_count * sizeof(*legacy_inverse));
    }
    legacy_pfaffian[local_qp] = component->pfaffian;
  }
  state_swap = workspace->accepted;
  workspace->accepted = workspace->candidate;
  workspace->candidate = state_swap;
  inverse_swap = workspace->accepted_inverse;
  workspace->accepted_inverse = workspace->candidate_inverse;
  workspace->candidate_inverse = inverse_swap;
  workspace->candidate.valid = 0;
  workspace->layout.candidate_ready = 0;
  return MVMC_PFAFFIAN_STATUS_OK;
}

void mvmc_classic_pfaffian_real_discard_candidate(
    MVMCClassicPfaffianRealWorkspace *workspace) {
  if (workspace == NULL) return;
  workspace->candidate.valid = 0;
  workspace->layout.candidate_ready = 0;
}

void mvmc_classic_pfaffian_complex_discard_candidate(
    MVMCClassicPfaffianComplexWorkspace *workspace) {
  if (workspace == NULL) return;
  workspace->candidate.valid = 0;
  workspace->layout.candidate_ready = 0;
}

const MVMCClassicPfaffianState *mvmc_classic_pfaffian_real_accepted(
    const MVMCClassicPfaffianRealWorkspace *workspace) {
  return workspace == NULL ? NULL : &workspace->accepted;
}

const MVMCClassicPfaffianState *mvmc_classic_pfaffian_real_candidate(
    const MVMCClassicPfaffianRealWorkspace *workspace) {
  return workspace == NULL || !workspace->layout.candidate_ready
             ? NULL
             : &workspace->candidate;
}

const double *mvmc_classic_pfaffian_real_accepted_inverse(
    const MVMCClassicPfaffianRealWorkspace *workspace) {
  return workspace == NULL || !workspace->accepted.valid
             ? NULL
             : workspace->accepted_inverse;
}

const double *mvmc_classic_pfaffian_real_candidate_inverse(
    const MVMCClassicPfaffianRealWorkspace *workspace) {
  return workspace == NULL || !workspace->layout.candidate_ready ||
                 !workspace->candidate.valid
             ? NULL
             : workspace->candidate_inverse;
}

const MVMCClassicPfaffianState *mvmc_classic_pfaffian_complex_accepted(
    const MVMCClassicPfaffianComplexWorkspace *workspace) {
  return workspace == NULL ? NULL : &workspace->accepted;
}

const MVMCClassicPfaffianState *mvmc_classic_pfaffian_complex_candidate(
    const MVMCClassicPfaffianComplexWorkspace *workspace) {
  return workspace == NULL || !workspace->layout.candidate_ready
             ? NULL
             : &workspace->candidate;
}

const double complex *mvmc_classic_pfaffian_complex_accepted_inverse(
    const MVMCClassicPfaffianComplexWorkspace *workspace) {
  return workspace == NULL || !workspace->accepted.valid
             ? NULL
             : workspace->accepted_inverse;
}

const double complex *mvmc_classic_pfaffian_complex_candidate_inverse(
    const MVMCClassicPfaffianComplexWorkspace *workspace) {
  return workspace == NULL || !workspace->layout.candidate_ready ||
                 !workspace->candidate.valid
             ? NULL
             : workspace->candidate_inverse;
}

int mvmc_classic_pfaffian_real_failure_local_qp(
    const MVMCClassicPfaffianRealWorkspace *workspace) {
  return workspace == NULL ? -1 : workspace->layout.failure_local_qp;
}

int mvmc_classic_pfaffian_real_failure_global_qp(
    const MVMCClassicPfaffianRealWorkspace *workspace) {
  return workspace == NULL ? -1 : workspace->layout.failure_global_qp;
}

int mvmc_classic_pfaffian_complex_failure_local_qp(
    const MVMCClassicPfaffianComplexWorkspace *workspace) {
  return workspace == NULL ? -1 : workspace->layout.failure_local_qp;
}

int mvmc_classic_pfaffian_complex_failure_global_qp(
    const MVMCClassicPfaffianComplexWorkspace *workspace) {
  return workspace == NULL ? -1 : workspace->layout.failure_global_qp;
}

size_t mvmc_classic_pfaffian_real_legacy_element_count(
    const MVMCClassicPfaffianRealWorkspace *workspace) {
  return workspace == NULL ? 0 : workspace->layout.legacy_element_count;
}

size_t mvmc_classic_pfaffian_complex_legacy_element_count(
    const MVMCClassicPfaffianComplexWorkspace *workspace) {
  return workspace == NULL ? 0 : workspace->layout.legacy_element_count;
}

static int copy_qp_range(const MVMCClassicPfaffianLayout *layout,
                         int *qp_total, int *qp_start, int *qp_end) {
  if (layout == NULL || qp_total == NULL || qp_start == NULL ||
      qp_end == NULL) {
    return 0;
  }
  *qp_total = layout->qp_total;
  *qp_start = layout->qp_start;
  *qp_end = layout->qp_end;
  return 1;
}

int mvmc_classic_pfaffian_real_qp_range(
    const MVMCClassicPfaffianRealWorkspace *workspace,
    int *qp_total, int *qp_start, int *qp_end) {
  return workspace != NULL &&
         copy_qp_range(&workspace->layout, qp_total, qp_start, qp_end);
}

int mvmc_classic_pfaffian_complex_qp_range(
    const MVMCClassicPfaffianComplexWorkspace *workspace,
    int *qp_total, int *qp_start, int *qp_end) {
  return workspace != NULL &&
         copy_qp_range(&workspace->layout, qp_total, qp_start, qp_end);
}

static int copy_layout(const MVMCClassicPfaffianLayout *layout,
                       int *nsize, int *qp_total,
                       int *qp_start, int *qp_end) {
  if (layout == NULL || nsize == NULL || qp_total == NULL ||
      qp_start == NULL || qp_end == NULL) {
    return 0;
  }
  *nsize = layout->nsize;
  return copy_qp_range(layout, qp_total, qp_start, qp_end);
}

int mvmc_classic_pfaffian_real_layout(
    const MVMCClassicPfaffianRealWorkspace *workspace,
    int *nsize, int *qp_total, int *qp_start, int *qp_end) {
  return workspace != NULL &&
         copy_layout(&workspace->layout, nsize, qp_total, qp_start, qp_end);
}

int mvmc_classic_pfaffian_complex_layout(
    const MVMCClassicPfaffianComplexWorkspace *workspace,
    int *nsize, int *qp_total, int *qp_start, int *qp_end) {
  return workspace != NULL &&
         copy_layout(&workspace->layout, nsize, qp_total, qp_start, qp_end);
}
