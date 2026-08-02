/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "classic_krylov_amplitude.h"

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)

#include "classic_pfaffian_matrix.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  size_t site_count;
  size_t up_electron_count;
  size_t down_electron_count;
  size_t orbital_count;
  size_t word_count;
  int pure_spin;
  int nsize;
  int qp_total;
  int qp_start;
  int qp_end;
  size_t local_qp_count;
  size_t slater_matrix_count;
  size_t matrix_count;
  double scaled_pivot_tolerance;
  MVMCClassicKrylovProjectionKind projection_kind;
  double projection_parameter;
  double complex *global_weights;
  int *ele_idx;
  int *ele_num;
  MVMCAbsolutePfaffianValueResult *components;
  MVMCAbsolutePfaffianValueResult *global_components;
  MVMCClassicPfaffianCollectiveWorkspace *collective_workspace;
} MVMCClassicKrylovAmplitudeCommon;

struct MVMCClassicKrylovRealAmplitudeWorkspace {
  MVMCClassicKrylovAmplitudeCommon common;
  double *slater_elm;
  double *matrix;
  MVMCAbsolutePfaffianRealValueWorkspace *value_workspace;
};

struct MVMCClassicKrylovComplexAmplitudeWorkspace {
  MVMCClassicKrylovAmplitudeCommon common;
  double complex *slater_elm;
  double complex *matrix;
  MVMCAbsolutePfaffianComplexValueWorkspace *value_workspace;
};

static int checked_multiply(size_t left, size_t right, size_t *product) {
  if (product == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *product = left * right;
  return 1;
}

static int checked_add(size_t left, size_t right, size_t *sum) {
  if (sum == NULL || right > SIZE_MAX - left) return 0;
  *sum = left + right;
  return 1;
}

static MVMCKrylovStatus map_create_status(MVMCPfaffianStatus status) {
  switch (status) {
    case MVMC_PFAFFIAN_STATUS_OK:
      return MVMC_KRYLOV_STATUS_OK;
    case MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE:
      return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    case MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT:
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    case MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE:
    case MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE:
      return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
}

static int valid_jastrow_binding(const int *const *binding,
                                  size_t site_count, int class_count) {
  size_t row, column;

  if (class_count == 0) return 1;
  if (binding == NULL) return 0;
  for (row = 0; row < site_count; ++row) {
    if (binding[row] == NULL) return 0;
    for (column = 0; column < site_count; ++column) {
      const int value = binding[row][column];
      if (row == column) {
        if (value != -1) return 0;
      } else if (value < 0 || value >= class_count ||
                 binding[column] == NULL ||
                 binding[column][row] != value) {
        return 0;
      }
    }
  }
  return 1;
}

static int valid_doublon_holon_binding(const int *const *binding,
                                       size_t site_count, int class_count,
                                       size_t neighbors_per_site) {
  int class_index;
  size_t value_index, value_count;

  if (class_count == 0) return 1;
  if (binding == NULL ||
      !checked_multiply(site_count, neighbors_per_site, &value_count)) {
    return 0;
  }
  for (class_index = 0; class_index < class_count; ++class_index) {
    if (binding[class_index] == NULL) return 0;
    for (value_index = 0; value_index < value_count; ++value_index) {
      const int value = binding[class_index][value_index];
      if (value < 0 || (size_t)value >= site_count) return 0;
    }
  }
  return 1;
}

static MVMCKrylovStatus initialize_common(
    const MVMCClassicKrylovAmplitudeLayout *layout,
    const double complex *global_weights,
    MVMCClassicKrylovAmplitudeCommon *common) {
  size_t nsize, orbital_count, word_count, matrix_count;
  size_t slater_dimension, slater_matrix_count;
  size_t expected_nproj;
  int index;
  int qp;

  if (layout == NULL || common == NULL || global_weights == NULL ||
      layout->site_count == 0 || layout->site_count > (size_t)INT_MAX / 2 ||
      layout->up_electron_count == 0 ||
      layout->up_electron_count != layout->down_electron_count ||
      layout->up_electron_count > layout->site_count ||
      layout->up_electron_count > (size_t)INT_MAX / 2 ||
      (layout->pure_spin != 0 && layout->pure_spin != 1) ||
      (layout->pure_spin &&
       layout->up_electron_count + layout->down_electron_count !=
           layout->site_count) ||
      (layout->qp_total != 1 && layout->qp_total != 4) ||
      layout->qp_start < 0 ||
      layout->qp_end < layout->qp_start ||
      layout->qp_end > layout->qp_total ||
      !isfinite(layout->scaled_pivot_tolerance) ||
      layout->nproj < 0 || layout->ngutzwiller_idx < 0 ||
      layout->njastrow_idx < 0 || layout->nspin_jastrow_idx < 0 ||
      layout->ndoublon_holon_2site_idx < 0 ||
      layout->ndoublon_holon_4site_idx < 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (layout->unsupported_features != 0) {
    return MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL;
  }
  expected_nproj = (size_t)layout->ngutzwiller_idx;
  if ((size_t)layout->njastrow_idx > SIZE_MAX - expected_nproj) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  expected_nproj += (size_t)layout->njastrow_idx;
  if ((size_t)layout->nspin_jastrow_idx > SIZE_MAX - expected_nproj) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  expected_nproj += (size_t)layout->nspin_jastrow_idx;
  if ((size_t)layout->ndoublon_holon_2site_idx >
          (SIZE_MAX - expected_nproj) / 6) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  expected_nproj += 6 * (size_t)layout->ndoublon_holon_2site_idx;
  if ((size_t)layout->ndoublon_holon_4site_idx >
          (SIZE_MAX - expected_nproj) / 10) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  expected_nproj += 10 * (size_t)layout->ndoublon_holon_4site_idx;
  if (expected_nproj > (size_t)INT_MAX ||
      (int)expected_nproj != layout->nproj) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (layout->ngutzwiller_idx != 0) {
    if (layout->gutzwiller_idx == NULL ||
        layout->projection_parameters == NULL) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    for (index = 0; index < (int)layout->site_count; ++index) {
      if (layout->gutzwiller_idx[index] < 0 ||
          layout->gutzwiller_idx[index] >= layout->ngutzwiller_idx) {
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      }
    }
  }
  for (index = 0; index < layout->nproj; ++index) {
    if (layout->projection_parameters == NULL ||
        !isfinite(creal(layout->projection_parameters[index])) ||
        !isfinite(cimag(layout->projection_parameters[index])) ||
        cimag(layout->projection_parameters[index]) != 0.0) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  if (!valid_jastrow_binding(layout->jastrow_idx, layout->site_count,
                              layout->njastrow_idx) ||
      !valid_jastrow_binding(layout->spin_jastrow_idx, layout->site_count,
                              layout->nspin_jastrow_idx) ||
      !valid_doublon_holon_binding(layout->doublon_holon_2site_idx,
                                    layout->site_count,
                                    layout->ndoublon_holon_2site_idx, 2) ||
      !valid_doublon_holon_binding(layout->doublon_holon_4site_idx,
                                    layout->site_count,
                                    layout->ndoublon_holon_4site_idx, 4)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (layout->njastrow_idx != 0 || layout->nspin_jastrow_idx != 0 ||
      layout->ndoublon_holon_2site_idx != 0 ||
      layout->ndoublon_holon_4site_idx != 0 ||
      layout->ngutzwiller_idx > 1) {
    return MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL;
  }
  if (layout->ngutzwiller_idx == 1) {
    for (index = 0; index < (int)layout->site_count; ++index) {
      if (layout->gutzwiller_idx[index] != 0) {
        return MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL;
      }
    }
  }
  nsize = layout->up_electron_count + layout->down_electron_count;
  if (nsize > (size_t)INT_MAX ||
      !checked_multiply(layout->site_count, 2, &orbital_count) ||
      orbital_count > SIZE_MAX - 63 ||
      !checked_multiply(nsize, nsize, &matrix_count) ||
      !checked_multiply(layout->site_count, 2, &slater_dimension) ||
      !checked_multiply(slater_dimension, slater_dimension,
                        &slater_matrix_count)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  word_count = (orbital_count + 63) / 64;
  for (qp = 0; qp < layout->qp_total; ++qp) {
    if (!isfinite(creal(global_weights[qp])) ||
        !isfinite(cimag(global_weights[qp]))) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  common->site_count = layout->site_count;
  common->up_electron_count = layout->up_electron_count;
  common->down_electron_count = layout->down_electron_count;
  common->orbital_count = orbital_count;
  common->word_count = word_count;
  common->pure_spin = layout->pure_spin;
  common->nsize = (int)nsize;
  common->qp_total = layout->qp_total;
  common->qp_start = layout->qp_start;
  common->qp_end = layout->qp_end;
  common->local_qp_count = (size_t)(layout->qp_end - layout->qp_start);
  common->slater_matrix_count = slater_matrix_count;
  common->matrix_count = matrix_count;
  common->scaled_pivot_tolerance = layout->scaled_pivot_tolerance;
  common->projection_kind = layout->ngutzwiller_idx == 0
                                ? MVMC_CLASSIC_KRYLOV_PROJECTION_NONE
                                : MVMC_CLASSIC_KRYLOV_PROJECTION_GUTZWILLER_ONLY;
  common->projection_parameter =
      layout->ngutzwiller_idx == 0
          ? 0.0 : creal(layout->projection_parameters[0]);
  return MVMC_KRYLOV_STATUS_OK;
}

static int allocate_common(MVMCClassicKrylovAmplitudeCommon *common) {
  if (common == NULL || common->nsize <= 0) return 0;
  common->ele_idx = (int *)calloc((size_t)common->nsize,
                                  sizeof(*common->ele_idx));
  common->ele_num = (int *)calloc(common->orbital_count,
                                  sizeof(*common->ele_num));
  common->global_weights = (double complex *)calloc(
      (size_t)common->qp_total, sizeof(*common->global_weights));
  if (common->local_qp_count != 0) {
    common->components = (MVMCAbsolutePfaffianValueResult *)calloc(
        common->local_qp_count, sizeof(*common->components));
  }
  common->global_components = (MVMCAbsolutePfaffianValueResult *)calloc(
      (size_t)common->qp_total, sizeof(*common->global_components));
  return common->ele_idx != NULL && common->ele_num != NULL &&
         common->global_weights != NULL &&
         common->global_components != NULL &&
         (common->local_qp_count == 0 || common->components != NULL);
}

static void destroy_common(MVMCClassicKrylovAmplitudeCommon *common) {
  if (common == NULL) return;
  mvmc_classic_pfaffian_collective_workspace_destroy(
      common->collective_workspace);
  free(common->global_components);
  free(common->components);
  free(common->global_weights);
  free(common->ele_num);
  free(common->ele_idx);
}

static MVMCPfaffianStatus audit_static_binding(
    MVMCClassicKrylovAmplitudeCommon *common, int scalar_kind,
    const void *slater, size_t slater_bytes,
    const double complex *global_weights) {
  uint64_t metadata[9];
  uint64_t tolerance_bits = 0;
  uint64_t projection_bits = 0;
  MVMCPfaffianStatus status;
  int all_equal = 0;

  memcpy(&tolerance_bits, &common->scaled_pivot_tolerance,
         sizeof(tolerance_bits));
  memcpy(&projection_bits, &common->projection_parameter,
         sizeof(projection_bits));
  metadata[0] = (uint64_t)common->site_count;
  metadata[1] = (uint64_t)common->up_electron_count;
  metadata[2] = (uint64_t)common->down_electron_count;
  metadata[3] = (uint64_t)common->pure_spin;
  metadata[4] = (uint64_t)common->qp_total;
  metadata[5] = (uint64_t)common->projection_kind;
  metadata[6] = (uint64_t)scalar_kind;
  metadata[7] = tolerance_bits;
  metadata[8] = projection_bits;
  status = mvmc_classic_pfaffian_collective_all_equal_u64(
      common->collective_workspace, metadata,
      sizeof(metadata) / sizeof(metadata[0]), &all_equal);
  if (status != MVMC_PFAFFIAN_STATUS_OK || !all_equal) {
    return status == MVMC_PFAFFIAN_STATUS_OK
               ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT : status;
  }
  status = mvmc_classic_pfaffian_collective_all_equal_bytes(
      common->collective_workspace, global_weights,
      (size_t)common->qp_total * sizeof(*global_weights), &all_equal);
  if (status != MVMC_PFAFFIAN_STATUS_OK || !all_equal) {
    return status == MVMC_PFAFFIAN_STATUS_OK
               ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT : status;
  }
  status = mvmc_classic_pfaffian_collective_all_equal_bytes(
      common->collective_workspace, slater, slater_bytes, &all_equal);
  if (status != MVMC_PFAFFIAN_STATUS_OK || !all_equal) {
    return status == MVMC_PFAFFIAN_STATUS_OK
               ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT : status;
  }
  return MVMC_PFAFFIAN_STATUS_OK;
}

static MVMCPfaffianStatus audit_qp_partition(
    MVMCClassicKrylovAmplitudeCommon *common) {
  MVMCProjectedAmplitudeResult local;
  MVMCClassicPfaffianCollectiveResult collective;

  memset(&local, 0, sizeof(local));
  local.valid = 1;
  local.regular_count = common->local_qp_count;
  return mvmc_classic_pfaffian_collective_aggregate(
      common->collective_workspace, MVMC_PFAFFIAN_STATUS_OK, &local,
      common->qp_total, common->qp_start, common->qp_end, -1, -1,
      &collective);
}

static int build_occupation(MVMCClassicKrylovAmplitudeCommon *common,
                            const uint64_t *words, size_t word_count) {
  size_t orbital;
  size_t up_count = 0;
  size_t down_count = 0;

  if (common == NULL || words == NULL || word_count != common->word_count) {
    return 0;
  }
  if ((common->orbital_count % 64) != 0 &&
      (words[word_count - 1] >> (common->orbital_count % 64)) != 0) {
    return 0;
  }
  memset(common->ele_num, 0,
         common->orbital_count * sizeof(*common->ele_num));
  for (orbital = 0; orbital < common->orbital_count; ++orbital) {
    const int occupied =
        (int)((words[orbital / 64] >> (orbital % 64)) & UINT64_C(1));
    common->ele_num[orbital] = occupied;
    if (!occupied) continue;
    if (orbital < common->site_count) {
      if (up_count >= common->up_electron_count) return 0;
      common->ele_idx[up_count++] = (int)orbital;
    } else {
      if (down_count >= common->down_electron_count) return 0;
      common->ele_idx[common->up_electron_count + down_count++] =
          (int)(orbital - common->site_count);
    }
  }
  if (up_count != common->up_electron_count ||
      down_count != common->down_electron_count) {
    return 0;
  }
  if (common->pure_spin) {
    size_t site;
    for (site = 0; site < common->site_count; ++site) {
      if (common->ele_num[site] +
              common->ele_num[common->site_count + site] != 1) {
        return 0;
      }
    }
  }
  return 1;
}

static MVMCKrylovStatus projection_log_value(
    MVMCClassicKrylovAmplitudeCommon *common, double *log_value) {
  size_t site;
  size_t double_occupancy = 0;

  *log_value = 0.0;
  if (common->projection_kind == MVMC_CLASSIC_KRYLOV_PROJECTION_NONE) {
    return MVMC_KRYLOV_STATUS_OK;
  }
  for (site = 0; site < common->site_count; ++site) {
    double_occupancy +=
        (size_t)(common->ele_num[site] *
                 common->ele_num[common->site_count + site]);
  }
  *log_value = common->projection_parameter * (double)double_occupancy;
  return isfinite(*log_value) ? MVMC_KRYLOV_STATUS_OK
                              : MVMC_KRYLOV_STATUS_NONFINITE;
}

static int neumaier_add(double value, double *sum, double *compensation) {
  const double updated = *sum + value;

  if (!isfinite(value) || !isfinite(updated)) return 0;
  if (fabs(*sum) >= fabs(value)) {
    *compensation += (*sum - updated) + value;
  } else {
    *compensation += (value - updated) + *sum;
  }
  if (!isfinite(*compensation)) return 0;
  *sum = updated;
  return 1;
}

static MVMCKrylovStatus deterministic_projected_amplitude(
    const MVMCClassicKrylovAmplitudeCommon *common,
    MVMCProjectedAmplitudeResult *result) {
  double real_sum = 0.0, real_compensation = 0.0;
  double imag_sum = 0.0, imag_compensation = 0.0;
  double abs_sum = 0.0, abs_compensation = 0.0;
  int qp;

  memset(result, 0, sizeof(*result));
  for (qp = 0; qp < common->qp_total; ++qp) {
    const MVMCAbsolutePfaffianValueResult *component =
        common->global_components + qp;
    double complex term;
    double magnitude;

    if (!isfinite(creal(component->pfaffian)) ||
        !isfinite(cimag(component->pfaffian))) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    switch (component->state) {
      case MVMC_PFAFFIAN_VALUE_WELL_PIVOTED:
        ++result->regular_count;
        break;
      case MVMC_PFAFFIAN_VALUE_NEAR_PIVOT:
        ++result->near_singular_count;
        break;
      case MVMC_PFAFFIAN_VALUE_SINGULAR:
        if (component->pfaffian != 0.0) {
          return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
        }
        ++result->singular_count;
        break;
      case MVMC_PFAFFIAN_VALUE_NONFINITE:
        return MVMC_KRYLOV_STATUS_NONFINITE;
      case MVMC_PFAFFIAN_VALUE_INVALID:
        return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      default:
        return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    term = common->global_weights[qp] * component->pfaffian;
    magnitude = cabs(term);
    if (!isfinite(creal(term)) || !isfinite(cimag(term)) ||
        !isfinite(magnitude) ||
        !neumaier_add(creal(term), &real_sum, &real_compensation) ||
        !neumaier_add(cimag(term), &imag_sum, &imag_compensation) ||
        !neumaier_add(magnitude, &abs_sum, &abs_compensation)) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  real_sum += real_compensation;
  imag_sum += imag_compensation;
  abs_sum += abs_compensation;
  if (!isfinite(real_sum) || !isfinite(imag_sum) || !isfinite(abs_sum) ||
      abs_sum < 0.0) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  result->total = real_sum + I * imag_sum;
  result->sum_abs = abs_sum;
  result->cancellation_ratio =
      abs_sum == 0.0 ? 0.0 : fmin(1.0, cabs(result->total) / abs_sum);
  if (!isfinite(result->cancellation_ratio)) {
    memset(result, 0, sizeof(*result));
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  result->valid = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus publish_amplitude(
    const MVMCClassicKrylovAmplitudeCommon *common,
    const MVMCProjectedAmplitudeResult *projected,
    double log_projection, uint64_t global_factorization_count,
    MVMCKrylovAmplitudeResult *result) {
  const double projection = exp(log_projection);
  double complex value;

  if (projection == 0.0 || !isfinite(projection)) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  value = projection * projected->total;
  if (!isfinite(creal(value)) || !isfinite(cimag(value))) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
#if SIZE_MAX > UINT64_MAX
  if (projected->regular_count > UINT64_MAX ||
      projected->near_singular_count > UINT64_MAX ||
      projected->singular_count > UINT64_MAX ||
      common->local_qp_count > UINT64_MAX) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
#endif
  result->value = value;
  result->total_zero = value == 0.0;
  result->regular_component_count =
      (uint64_t)projected->regular_count;
  result->near_pivot_component_count =
      (uint64_t)projected->near_singular_count;
  result->singular_component_count =
      (uint64_t)projected->singular_count;
  result->local_factorization_count = (uint64_t)common->local_qp_count;
  result->global_factorization_count = global_factorization_count;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_classic_krylov_real_amplitude_workspace_create(
    const MVMCClassicKrylovAmplitudeLayout *layout,
    const double *slater_elm, const double complex *global_weights,
    MVMCClassicKrylovRealAmplitudeWorkspace **workspace) {
  MVMCClassicKrylovRealAmplitudeWorkspace *created = NULL;
  MVMCClassicPfaffianCollectiveWorkspace *collective_workspace = NULL;
  MVMCKrylovStatus local_status = MVMC_KRYLOV_STATUS_OK;
  int global_status = (int)MVMC_KRYLOV_STATUS_OK;
  MVMCPfaffianStatus pf_status;
  int local_ready = 0;
  int all_ready = 0;
  size_t slater_count = 0;
  size_t slater_bytes = 0;
  int qp;

  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (layout == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  pf_status = mvmc_classic_pfaffian_collective_workspace_create(
      layout->communicator, &collective_workspace);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK) return map_create_status(pf_status);
  created = (MVMCClassicKrylovRealAmplitudeWorkspace *)calloc(
      1, sizeof(*created));
  pf_status = mvmc_classic_pfaffian_collective_all_true(
      collective_workspace, created != NULL, &all_ready);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK || !all_ready ||
      created == NULL) {
    mvmc_classic_pfaffian_collective_workspace_destroy(collective_workspace);
    free(created);
    return pf_status == MVMC_PFAFFIAN_STATUS_OK
               ? MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE
               : map_create_status(pf_status);
  }
  created->common.collective_workspace = collective_workspace;
  local_status = slater_elm == NULL
                     ? MVMC_KRYLOV_STATUS_INVALID_ARGUMENT
                     : initialize_common(layout, global_weights,
                                         &created->common);
  if (local_status == MVMC_KRYLOV_STATUS_OK) {
    for (qp = 0; qp < layout->qp_total; ++qp) {
      if (cimag(global_weights[qp]) != 0.0) {
        local_status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      }
    }
  }
  if (mvmc_classic_pfaffian_collective_max_int(
          collective_workspace, (int)local_status, &global_status) !=
      MVMC_PFAFFIAN_STATUS_OK) {
    mvmc_classic_krylov_real_amplitude_workspace_destroy(created);
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (global_status != (int)MVMC_KRYLOV_STATUS_OK) {
    mvmc_classic_krylov_real_amplitude_workspace_destroy(created);
    return (MVMCKrylovStatus)global_status;
  }
  local_ready =
      checked_multiply((size_t)created->common.qp_total,
                       created->common.slater_matrix_count, &slater_count) &&
      checked_multiply(slater_count, sizeof(*created->slater_elm),
                       &slater_bytes) &&
      allocate_common(&created->common);
  if (local_ready) {
    created->slater_elm = (double *)malloc(slater_bytes);
    created->matrix = (double *)malloc(
        created->common.matrix_count * sizeof(*created->matrix));
    local_ready = created->slater_elm != NULL && created->matrix != NULL;
  }
  if (local_ready) {
    pf_status = mvmc_absolute_pfaffian_real_value_workspace_create(
        created->common.nsize, &created->value_workspace);
    local_ready = pf_status == MVMC_PFAFFIAN_STATUS_OK;
  }
  pf_status = mvmc_classic_pfaffian_collective_all_true(
      created->common.collective_workspace, local_ready, &all_ready);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK || !all_ready) {
    mvmc_classic_krylov_real_amplitude_workspace_destroy(created);
    return pf_status != MVMC_PFAFFIAN_STATUS_OK
               ? map_create_status(pf_status)
               : MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
  pf_status = audit_static_binding(
      &created->common, 0, slater_elm, slater_bytes, global_weights);
  if (pf_status == MVMC_PFAFFIAN_STATUS_OK) {
    pf_status = audit_qp_partition(&created->common);
  }
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK) {
    mvmc_classic_krylov_real_amplitude_workspace_destroy(created);
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  memcpy(created->slater_elm, slater_elm, slater_bytes);
  memcpy(created->common.global_weights, global_weights,
         (size_t)created->common.qp_total *
             sizeof(*created->common.global_weights));
  *workspace = created;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_classic_krylov_real_amplitude_workspace_destroy(
    MVMCClassicKrylovRealAmplitudeWorkspace *workspace) {
  if (workspace == NULL) return;
  mvmc_absolute_pfaffian_real_value_workspace_destroy(
      workspace->value_workspace);
  free(workspace->matrix);
  free(workspace->slater_elm);
  destroy_common(&workspace->common);
  free(workspace);
}

size_t mvmc_classic_krylov_real_amplitude_workspace_bytes(
    const MVMCClassicKrylovRealAmplitudeWorkspace *workspace) {
  size_t bytes;
  size_t item;

  if (workspace == NULL) return 0;
  bytes = sizeof(*workspace);
#define ADD_REAL_BYTES(count, type)                                      \
  do {                                                                    \
    if (!checked_multiply((count), sizeof(type), &item) ||                \
        !checked_add(bytes, item, &bytes)) return 0;                       \
  } while (0)
  ADD_REAL_BYTES((size_t)workspace->common.nsize, int);
  ADD_REAL_BYTES(workspace->common.orbital_count, int);
  ADD_REAL_BYTES((size_t)workspace->common.qp_total, double complex);
  ADD_REAL_BYTES(workspace->common.local_qp_count,
                 MVMCAbsolutePfaffianValueResult);
  ADD_REAL_BYTES((size_t)workspace->common.qp_total,
                 MVMCAbsolutePfaffianValueResult);
  ADD_REAL_BYTES((size_t)workspace->common.qp_total *
                     workspace->common.slater_matrix_count,
                 double);
  ADD_REAL_BYTES(workspace->common.matrix_count, double);
#undef ADD_REAL_BYTES
  if (!checked_add(bytes,
                   mvmc_absolute_pfaffian_real_value_workspace_bytes(
                       workspace->value_workspace),
                   &bytes) ||
      !checked_add(bytes,
                   mvmc_classic_pfaffian_collective_workspace_bytes(
                       workspace->common.collective_workspace),
                   &bytes)) {
    return 0;
  }
  return bytes;
}

MVMCKrylovStatus mvmc_classic_krylov_complex_amplitude_workspace_create(
    const MVMCClassicKrylovAmplitudeLayout *layout,
    const double complex *slater_elm,
    const double complex *global_weights,
    MVMCClassicKrylovComplexAmplitudeWorkspace **workspace) {
  MVMCClassicKrylovComplexAmplitudeWorkspace *created = NULL;
  MVMCClassicPfaffianCollectiveWorkspace *collective_workspace = NULL;
  MVMCKrylovStatus local_status = MVMC_KRYLOV_STATUS_OK;
  int global_status = (int)MVMC_KRYLOV_STATUS_OK;
  MVMCPfaffianStatus pf_status;
  int local_ready = 0;
  int all_ready = 0;
  size_t slater_count = 0;
  size_t slater_bytes = 0;

  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (layout == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  pf_status = mvmc_classic_pfaffian_collective_workspace_create(
      layout->communicator, &collective_workspace);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK) return map_create_status(pf_status);
  created = (MVMCClassicKrylovComplexAmplitudeWorkspace *)calloc(
      1, sizeof(*created));
  pf_status = mvmc_classic_pfaffian_collective_all_true(
      collective_workspace, created != NULL, &all_ready);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK || !all_ready ||
      created == NULL) {
    mvmc_classic_pfaffian_collective_workspace_destroy(collective_workspace);
    free(created);
    return pf_status == MVMC_PFAFFIAN_STATUS_OK
               ? MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE
               : map_create_status(pf_status);
  }
  created->common.collective_workspace = collective_workspace;
  local_status = slater_elm == NULL
                     ? MVMC_KRYLOV_STATUS_INVALID_ARGUMENT
                     : initialize_common(layout, global_weights,
                                         &created->common);
  if (mvmc_classic_pfaffian_collective_max_int(
          collective_workspace, (int)local_status, &global_status) !=
      MVMC_PFAFFIAN_STATUS_OK) {
    mvmc_classic_krylov_complex_amplitude_workspace_destroy(created);
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (global_status != (int)MVMC_KRYLOV_STATUS_OK) {
    mvmc_classic_krylov_complex_amplitude_workspace_destroy(created);
    return (MVMCKrylovStatus)global_status;
  }
  local_ready =
      checked_multiply((size_t)created->common.qp_total,
                       created->common.slater_matrix_count, &slater_count) &&
      checked_multiply(slater_count, sizeof(*created->slater_elm),
                       &slater_bytes) &&
      allocate_common(&created->common);
  if (local_ready) {
    created->slater_elm =
        (double complex *)malloc(slater_bytes);
    created->matrix = (double complex *)malloc(
        created->common.matrix_count * sizeof(*created->matrix));
    local_ready = created->slater_elm != NULL && created->matrix != NULL;
  }
  if (local_ready) {
    pf_status = mvmc_absolute_pfaffian_complex_value_workspace_create(
        created->common.nsize, &created->value_workspace);
    local_ready = pf_status == MVMC_PFAFFIAN_STATUS_OK;
  }
  pf_status = mvmc_classic_pfaffian_collective_all_true(
      created->common.collective_workspace, local_ready, &all_ready);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK || !all_ready) {
    mvmc_classic_krylov_complex_amplitude_workspace_destroy(created);
    return pf_status != MVMC_PFAFFIAN_STATUS_OK
               ? map_create_status(pf_status)
               : MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
  pf_status = audit_static_binding(
      &created->common, 1, slater_elm, slater_bytes, global_weights);
  if (pf_status == MVMC_PFAFFIAN_STATUS_OK) {
    pf_status = audit_qp_partition(&created->common);
  }
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK) {
    mvmc_classic_krylov_complex_amplitude_workspace_destroy(created);
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  memcpy(created->slater_elm, slater_elm, slater_bytes);
  memcpy(created->common.global_weights, global_weights,
         (size_t)created->common.qp_total *
             sizeof(*created->common.global_weights));
  *workspace = created;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_classic_krylov_complex_amplitude_workspace_destroy(
    MVMCClassicKrylovComplexAmplitudeWorkspace *workspace) {
  if (workspace == NULL) return;
  mvmc_absolute_pfaffian_complex_value_workspace_destroy(
      workspace->value_workspace);
  free(workspace->matrix);
  free(workspace->slater_elm);
  destroy_common(&workspace->common);
  free(workspace);
}

size_t mvmc_classic_krylov_complex_amplitude_workspace_bytes(
    const MVMCClassicKrylovComplexAmplitudeWorkspace *workspace) {
  size_t bytes;
  size_t item;

  if (workspace == NULL) return 0;
  bytes = sizeof(*workspace);
#define ADD_COMPLEX_BYTES(count, type)                                   \
  do {                                                                    \
    if (!checked_multiply((count), sizeof(type), &item) ||                \
        !checked_add(bytes, item, &bytes)) return 0;                       \
  } while (0)
  ADD_COMPLEX_BYTES((size_t)workspace->common.nsize, int);
  ADD_COMPLEX_BYTES(workspace->common.orbital_count, int);
  ADD_COMPLEX_BYTES((size_t)workspace->common.qp_total, double complex);
  ADD_COMPLEX_BYTES(workspace->common.local_qp_count,
                    MVMCAbsolutePfaffianValueResult);
  ADD_COMPLEX_BYTES((size_t)workspace->common.qp_total,
                    MVMCAbsolutePfaffianValueResult);
  ADD_COMPLEX_BYTES((size_t)workspace->common.qp_total *
                        workspace->common.slater_matrix_count,
                    double complex);
  ADD_COMPLEX_BYTES(workspace->common.matrix_count, double complex);
#undef ADD_COMPLEX_BYTES
  if (!checked_add(bytes,
                   mvmc_absolute_pfaffian_complex_value_workspace_bytes(
                       workspace->value_workspace),
                   &bytes) ||
      !checked_add(bytes,
                   mvmc_classic_pfaffian_collective_workspace_bytes(
                       workspace->common.collective_workspace),
                   &bytes)) {
    return 0;
  }
  return bytes;
}

MVMCKrylovStatus mvmc_classic_krylov_real_amplitude(
    const uint64_t *configuration_words, size_t word_count,
    void *context, MVMCKrylovAmplitudeResult *result) {
  MVMCClassicKrylovRealAmplitudeWorkspace *workspace =
      (MVMCClassicKrylovRealAmplitudeWorkspace *)context;
  MVMCClassicKrylovAmplitudeCommon *common;
  MVMCProjectedAmplitudeResult projected;
  MVMCPfaffianStatus local_status = MVMC_PFAFFIAN_STATUS_OK;
  MVMCPfaffianStatus pf_status;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  double log_projection = 0.0;
  uint64_t local_factorizations = 0;
  uint64_t global_factorizations = 0;
  size_t local_qp;
  int all_equal = 0;
  int all_valid = 0;
  int global_status = 0;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(result, 0, sizeof(*result));
  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  common = &workspace->common;
  pf_status = mvmc_classic_pfaffian_collective_all_equal_u64(
      common->collective_workspace, configuration_words, word_count,
      &all_equal);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK || !all_equal) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  pf_status = mvmc_classic_pfaffian_collective_all_true(
      common->collective_workspace,
      build_occupation(common, configuration_words, word_count), &all_valid);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK || !all_valid) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = projection_log_value(common, &log_projection);
  pf_status = mvmc_classic_pfaffian_collective_max_int(
      common->collective_workspace, (int)status, &global_status);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (global_status != (int)MVMC_KRYLOV_STATUS_OK) {
    return (MVMCKrylovStatus)global_status;
  }
  for (local_qp = 0; local_qp < common->local_qp_count; ++local_qp) {
    const int global_qp = common->qp_start + (int)local_qp;
    local_status = mvmc_classic_pfaffian_build_real_matrix(
        workspace->slater_elm +
            (size_t)global_qp * common->slater_matrix_count,
        (int)common->site_count, common->nsize, common->ele_idx,
        workspace->matrix);
    if (local_status == MVMC_PFAFFIAN_STATUS_OK) {
      ++local_factorizations;
      local_status = mvmc_absolute_pfaffian_real_value_with_workspace(
          workspace->value_workspace, workspace->matrix, common->nsize,
          common->nsize, common->scaled_pivot_tolerance,
          &common->components[local_qp]);
    }
    if (local_status != MVMC_PFAFFIAN_STATUS_OK) {
      status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
      break;
    }
    if (cimag(common->components[local_qp].pfaffian) != 0.0 ||
        common->components[local_qp].state ==
            MVMC_PFAFFIAN_VALUE_INVALID) {
      status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
      break;
    }
    if (common->components[local_qp].state ==
        MVMC_PFAFFIAN_VALUE_NONFINITE) {
      status = MVMC_KRYLOV_STATUS_NONFINITE;
      break;
    }
    if (common->components[local_qp].state !=
            MVMC_PFAFFIAN_VALUE_WELL_PIVOTED &&
        common->components[local_qp].state !=
            MVMC_PFAFFIAN_VALUE_NEAR_PIVOT &&
        common->components[local_qp].state !=
            MVMC_PFAFFIAN_VALUE_SINGULAR) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      break;
    }
  }
  if (mvmc_classic_pfaffian_collective_max_int(
          common->collective_workspace, (int)status, &global_status) !=
          MVMC_PFAFFIAN_STATUS_OK ||
      mvmc_classic_pfaffian_collective_sum_u64(
          common->collective_workspace, local_factorizations,
          &global_factorizations) != MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (global_status != (int)MVMC_KRYLOV_STATUS_OK) {
    return (MVMCKrylovStatus)global_status;
  }
  if (global_factorizations != (uint64_t)common->qp_total) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  pf_status = mvmc_classic_pfaffian_collective_gather_value_components(
      common->collective_workspace, common->components,
      common->local_qp_count, common->qp_total, common->qp_start,
      common->qp_end, common->global_components);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  status = deterministic_projected_amplitude(common, &projected);
  if (mvmc_classic_pfaffian_collective_max_int(
          common->collective_workspace, (int)status, &global_status) !=
      MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (global_status != (int)MVMC_KRYLOV_STATUS_OK) {
    return (MVMCKrylovStatus)global_status;
  }
  return publish_amplitude(common, &projected, log_projection,
                           global_factorizations, result);
}

MVMCKrylovStatus mvmc_classic_krylov_complex_amplitude(
    const uint64_t *configuration_words, size_t word_count,
    void *context, MVMCKrylovAmplitudeResult *result) {
  MVMCClassicKrylovComplexAmplitudeWorkspace *workspace =
      (MVMCClassicKrylovComplexAmplitudeWorkspace *)context;
  MVMCClassicKrylovAmplitudeCommon *common;
  MVMCProjectedAmplitudeResult projected;
  MVMCPfaffianStatus local_status = MVMC_PFAFFIAN_STATUS_OK;
  MVMCPfaffianStatus pf_status;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  double log_projection = 0.0;
  uint64_t local_factorizations = 0;
  uint64_t global_factorizations = 0;
  size_t local_qp;
  int all_equal = 0;
  int all_valid = 0;
  int global_status = 0;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(result, 0, sizeof(*result));
  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  common = &workspace->common;
  pf_status = mvmc_classic_pfaffian_collective_all_equal_u64(
      common->collective_workspace, configuration_words, word_count,
      &all_equal);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK || !all_equal) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  pf_status = mvmc_classic_pfaffian_collective_all_true(
      common->collective_workspace,
      build_occupation(common, configuration_words, word_count), &all_valid);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK || !all_valid) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = projection_log_value(common, &log_projection);
  pf_status = mvmc_classic_pfaffian_collective_max_int(
      common->collective_workspace, (int)status, &global_status);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (global_status != (int)MVMC_KRYLOV_STATUS_OK) {
    return (MVMCKrylovStatus)global_status;
  }
  for (local_qp = 0; local_qp < common->local_qp_count; ++local_qp) {
    const int global_qp = common->qp_start + (int)local_qp;
    local_status = mvmc_classic_pfaffian_build_complex_matrix(
        workspace->slater_elm +
            (size_t)global_qp * common->slater_matrix_count,
        (int)common->site_count, common->nsize, common->ele_idx,
        workspace->matrix);
    if (local_status == MVMC_PFAFFIAN_STATUS_OK) {
      ++local_factorizations;
      local_status = mvmc_absolute_pfaffian_complex_value_with_workspace(
          workspace->value_workspace, workspace->matrix, common->nsize,
          common->nsize, common->scaled_pivot_tolerance,
          &common->components[local_qp]);
    }
    if (local_status != MVMC_PFAFFIAN_STATUS_OK) {
      status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
      break;
    }
    if (common->components[local_qp].state ==
        MVMC_PFAFFIAN_VALUE_NONFINITE) {
      status = MVMC_KRYLOV_STATUS_NONFINITE;
      break;
    }
    if (common->components[local_qp].state !=
            MVMC_PFAFFIAN_VALUE_WELL_PIVOTED &&
        common->components[local_qp].state !=
            MVMC_PFAFFIAN_VALUE_NEAR_PIVOT &&
        common->components[local_qp].state !=
            MVMC_PFAFFIAN_VALUE_SINGULAR) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      break;
    }
  }
  if (mvmc_classic_pfaffian_collective_max_int(
          common->collective_workspace, (int)status, &global_status) !=
          MVMC_PFAFFIAN_STATUS_OK ||
      mvmc_classic_pfaffian_collective_sum_u64(
          common->collective_workspace, local_factorizations,
          &global_factorizations) != MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (global_status != (int)MVMC_KRYLOV_STATUS_OK) {
    return (MVMCKrylovStatus)global_status;
  }
  if (global_factorizations != (uint64_t)common->qp_total) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  pf_status = mvmc_classic_pfaffian_collective_gather_value_components(
      common->collective_workspace, common->components,
      common->local_qp_count, common->qp_total, common->qp_start,
      common->qp_end, common->global_components);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  status = deterministic_projected_amplitude(common, &projected);
  if (mvmc_classic_pfaffian_collective_max_int(
          common->collective_workspace, (int)status, &global_status) !=
      MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (global_status != (int)MVMC_KRYLOV_STATUS_OK) {
    return (MVMCKrylovStatus)global_status;
  }
  return publish_amplitude(common, &projected, log_projection,
                           global_factorizations, result);
}

#endif /* MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE */
