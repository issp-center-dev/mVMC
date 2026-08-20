/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "bounded_krylov_collective.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "bounded_krylov_collective.c requires the power-Lanczos core"
#endif

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct MVMCKrylovBoundedCollectiveWorkspace {
  int rank;
  int size;
  size_t allocated_bytes;
  size_t max_qp_components;
  size_t max_accumulator_count;
  int *ranges;
  int *byte_counts;
  int *byte_displacements;
  uint64_t *rank_values;
  MVMCAbsolutePfaffianScaledValueResult *component_buffer;
  MVMCScaledComplex *projection_terms;
  MVMCScaledComplex *rank_accumulator_buffer;
  MVMCScaledComplex *rank_order_scratch;
  MVMCScaledComplex *merge_buffer;
#ifdef _mpi_use
  MPI_Comm communicator;
#endif
};

static int checked_add_size(size_t left, size_t right, size_t *result) {
  if (result == NULL || left > SIZE_MAX - right) return 0;
  *result = left + right;
  return 1;
}

static int checked_multiply_size(size_t left, size_t right, size_t *result) {
  if (result == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *result = left * right;
  return 1;
}

static int reserve_aligned(size_t *offset, size_t count, size_t item_size,
                           size_t alignment, size_t *start) {
  size_t aligned;
  size_t bytes;
  if (offset == NULL || start == NULL || alignment == 0 ||
      (alignment & (alignment - 1)) != 0 ||
      *offset > SIZE_MAX - (alignment - 1)) return 0;
  aligned = (*offset + alignment - 1) & ~(alignment - 1);
  if (!checked_multiply_size(count, item_size, &bytes) ||
      !checked_add_size(aligned, bytes, offset)) return 0;
  *start = aligned;
  return 1;
}

static int valid_status(MVMCKrylovStatus status) {
  return status >= MVMC_KRYLOV_STATUS_OK &&
         status <= MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

static uint64_t double_bits(double value) {
  uint64_t bits;
  memcpy(&bits, &value, sizeof(bits));
  return bits;
}

static void hash_u64(uint64_t *hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    *hash ^= value & UINT64_C(0xff);
    *hash *= UINT64_C(1099511628211);
    value >>= 8;
  }
}

static int finite_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static int exact_scaled_weight(double complex weight,
                               MVMCScaledComplex *value) {
  const double scale = fmax(fabs(creal(weight)), fabs(cimag(weight)));
  double complex scaled;
  double magnitude;
  if (value == NULL || !finite_complex(weight)) return 0;
  if (scale == 0.0) {
    return mvmc_scaled_complex_make_exact_zero(value) ==
           MVMC_PFAFFIAN_STATUS_OK;
  }
  scaled = creal(weight) / scale + I * (cimag(weight) / scale);
  magnitude = hypot(creal(scaled), cimag(scaled));
  if (!isfinite(magnitude) || magnitude == 0.0) return 0;
  return mvmc_scaled_complex_make_finite(
             scaled / magnitude, log(scale) + log(magnitude),
             -INFINITY, value) == MVMC_PFAFFIAN_STATUS_OK;
}

static void reset_result(MVMCKrylovBoundedCollectiveResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  result->failure_rank = -1;
}

static MVMCKrylovStatus synchronize_status(
    MVMCKrylovBoundedCollectiveWorkspace *workspace,
    MVMCKrylovStatus local_status, int result_available,
    MVMCKrylovBoundedCollectiveResult *result) {
  MVMCKrylovStatus effective = valid_status(local_status)
                                      ? local_status
                                      : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  int selected_rank = 0;
#ifndef _mpi_use
  (void)workspace;
#endif
  if (!result_available && effective == MVMC_KRYLOV_STATUS_OK) {
    effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#ifdef _mpi_use
  {
    struct {
      int value;
      int index;
    } local_selection, global_selection;
    local_selection.value = (int)effective;
    local_selection.index = workspace->rank;
    if (MPI_Allreduce(&local_selection, &global_selection, 1, MPI_2INT,
                      MPI_MAXLOC, workspace->communicator) != MPI_SUCCESS) {
      if (result_available) {
        result->status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
      }
      return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
    effective = (MVMCKrylovStatus)global_selection.value;
    selected_rank = global_selection.index;
  }
#endif
  if (result_available) {
    result->status = effective;
    result->failure_rank =
        effective == MVMC_KRYLOV_STATUS_OK ? -1 : selected_rank;
    result->valid = effective == MVMC_KRYLOV_STATUS_OK;
  }
  return effective;
}

static int valid_scaled_component(
    const MVMCAbsolutePfaffianScaledValueResult *component,
    MVMCKrylovStatus *status) {
  if (component == NULL || status == NULL) return 0;
  if (component->factor_state < MVMC_PFAFFIAN_VALUE_WELL_PIVOTED ||
      component->factor_state > MVMC_PFAFFIAN_VALUE_INVALID ||
      component->factor_state == MVMC_PFAFFIAN_VALUE_INVALID ||
      component->factor_info < 0 || component->matrix_scale < 0.0 ||
      component->scaled_min_pivot < 0.0) {
    *status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    return 0;
  }
  if (component->factor_state == MVMC_PFAFFIAN_VALUE_NONFINITE ||
      !isfinite(component->matrix_scale) ||
      !isfinite(component->scaled_min_pivot)) {
    *status = MVMC_KRYLOV_STATUS_NONFINITE;
    return 0;
  }
  if (!mvmc_scaled_complex_is_valid(&component->value) ||
      component->value.state < MVMC_SCALED_COMPLEX_FINITE_NONZERO ||
      component->value.state > MVMC_SCALED_COMPLEX_NONFINITE) {
    *status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    return 0;
  }
  if (component->value.state == MVMC_SCALED_COMPLEX_NONFINITE) {
    *status = MVMC_KRYLOV_STATUS_NONFINITE;
    return 0;
  }
  if ((component->factor_state == MVMC_PFAFFIAN_VALUE_WELL_PIVOTED &&
       component->value.state != MVMC_SCALED_COMPLEX_FINITE_NONZERO) ||
      (component->factor_state == MVMC_PFAFFIAN_VALUE_NEAR_PIVOT &&
       component->value.state != MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
       component->value.state != MVMC_SCALED_COMPLEX_NUMERIC_ZERO) ||
      (component->factor_state == MVMC_PFAFFIAN_VALUE_SINGULAR &&
       component->value.state != MVMC_SCALED_COMPLEX_EXACT_ZERO)) {
    *status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    return 0;
  }
  return 1;
}

static int first_invalid_partition(const int *ranges, int size,
                                   size_t maximum) {
  const int qp_total = ranges[0];
  int rank;
  if (qp_total < 0 || (size_t)qp_total > maximum) return 0;
  for (rank = 0; rank < size; ++rank) {
    const int *range = ranges + (size_t)rank * 3;
    if (range[0] != qp_total || range[1] < 0 || range[2] < range[1] ||
        range[2] > qp_total) return rank;
    if ((rank == 0 && range[1] != 0) ||
        (rank > 0 &&
         range[1] != ranges[(size_t)(rank - 1) * 3 + 2])) return rank;
  }
  return ranges[(size_t)(size - 1) * 3 + 2] == qp_total ? -1 : size - 1;
}

MVMCKrylovStatus mvmc_bounded_krylov_collective_workspace_create(
    MVMCKrylovBoundedCommunicator communicator,
    size_t max_qp_components, size_t max_accumulator_count,
    MVMCKrylovBoundedCollectiveWorkspace **workspace) {
  MVMCKrylovBoundedCollectiveWorkspace *created = NULL;
  int size = 1;
  int rank = 0;
  int local_valid = workspace != NULL && max_qp_components > 0 &&
                    max_accumulator_count > 0;
  size_t offset = sizeof(*created);
  size_t ranges_offset;
  size_t counts_offset;
  size_t displacements_offset;
  size_t rank_values_offset;
  size_t component_offset;
  size_t projection_term_offset;
  size_t rank_accumulator_offset;
  size_t rank_scratch_offset;
  size_t merge_offset;
  size_t range_count = 0;
  size_t rank_accumulator_count = 0;
#ifdef _mpi_use
  MPI_Comm duplicated = MPI_COMM_NULL;
  int initialized = 0;
  int all_valid = 0;
  uint64_t local_capacities[2];
  uint64_t minimum_capacities[2];
  uint64_t maximum_capacities[2];
  int allocation_ok;
  int all_allocation_ok;
  if (workspace != NULL) *workspace = NULL;
  if (communicator == MPI_COMM_NULL ||
      MPI_Initialized(&initialized) != MPI_SUCCESS || !initialized ||
      MPI_Comm_dup(communicator, &duplicated) != MPI_SUCCESS ||
      MPI_Comm_set_errhandler(duplicated, MPI_ERRORS_RETURN) != MPI_SUCCESS ||
      MPI_Comm_size(duplicated, &size) != MPI_SUCCESS || size <= 0 ||
      MPI_Comm_rank(duplicated, &rank) != MPI_SUCCESS) {
    if (duplicated != MPI_COMM_NULL) MPI_Comm_free(&duplicated);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#if SIZE_MAX > UINT64_MAX
  if (max_qp_components > UINT64_MAX ||
      max_accumulator_count > UINT64_MAX) local_valid = 0;
#endif
  local_capacities[0] = local_valid ? (uint64_t)max_qp_components : 0;
  local_capacities[1] = local_valid ? (uint64_t)max_accumulator_count : 0;
  if (MPI_Allreduce(&local_valid, &all_valid, 1, MPI_INT, MPI_MIN,
                    duplicated) != MPI_SUCCESS ||
      MPI_Allreduce(local_capacities, minimum_capacities, 2, MPI_UINT64_T,
                    MPI_MIN, duplicated) != MPI_SUCCESS ||
      MPI_Allreduce(local_capacities, maximum_capacities, 2, MPI_UINT64_T,
                    MPI_MAX, duplicated) != MPI_SUCCESS) {
    MPI_Comm_free(&duplicated);
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (!all_valid ||
      memcmp(minimum_capacities, maximum_capacities,
             sizeof(minimum_capacities)) != 0) {
    MPI_Comm_free(&duplicated);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#else
  (void)communicator;
  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (!local_valid) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
#endif
  if (max_qp_components > (size_t)INT_MAX /
                              sizeof(MVMCAbsolutePfaffianScaledValueResult) ||
      max_accumulator_count >
          (size_t)INT_MAX / sizeof(MVMCScaledComplex) ||
      !checked_multiply_size((size_t)size, 3, &range_count) ||
      !checked_multiply_size((size_t)size, max_accumulator_count,
                             &rank_accumulator_count) ||
      !reserve_aligned(&offset, range_count, sizeof(int),
                       _Alignof(int), &ranges_offset) ||
      !reserve_aligned(&offset, (size_t)size, sizeof(int),
                       _Alignof(int), &counts_offset) ||
      !reserve_aligned(&offset, (size_t)size, sizeof(int),
                       _Alignof(int), &displacements_offset) ||
      !reserve_aligned(&offset, (size_t)size, sizeof(uint64_t),
                       _Alignof(uint64_t), &rank_values_offset) ||
      !reserve_aligned(&offset, max_qp_components,
                       sizeof(MVMCAbsolutePfaffianScaledValueResult),
                       _Alignof(MVMCAbsolutePfaffianScaledValueResult),
                       &component_offset) ||
      !reserve_aligned(&offset, max_qp_components,
                       sizeof(MVMCScaledComplex),
                       _Alignof(MVMCScaledComplex),
                       &projection_term_offset) ||
      !reserve_aligned(&offset, rank_accumulator_count,
                       sizeof(MVMCScaledComplex),
                       _Alignof(MVMCScaledComplex),
                       &rank_accumulator_offset) ||
      !reserve_aligned(&offset, (size_t)size, sizeof(MVMCScaledComplex),
                       _Alignof(MVMCScaledComplex), &rank_scratch_offset) ||
      !reserve_aligned(&offset, max_accumulator_count,
                       sizeof(MVMCScaledComplex),
                       _Alignof(MVMCScaledComplex), &merge_offset)) {
#ifdef _mpi_use
    MPI_Comm_free(&duplicated);
#endif
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  created = (MVMCKrylovBoundedCollectiveWorkspace *)calloc(1, offset);
#ifdef _mpi_use
  allocation_ok = created != NULL;
  if (MPI_Allreduce(&allocation_ok, &all_allocation_ok, 1, MPI_INT, MPI_MIN,
                    duplicated) != MPI_SUCCESS) {
    free(created);
    MPI_Comm_free(&duplicated);
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (!all_allocation_ok) {
    free(created);
    MPI_Comm_free(&duplicated);
    return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
#else
  if (created == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
#endif
  created->rank = rank;
  created->size = size;
  created->allocated_bytes = offset;
  created->max_qp_components = max_qp_components;
  created->max_accumulator_count = max_accumulator_count;
  created->ranges = (int *)((unsigned char *)created + ranges_offset);
  created->byte_counts = (int *)((unsigned char *)created + counts_offset);
  created->byte_displacements =
      (int *)((unsigned char *)created + displacements_offset);
  created->rank_values =
      (uint64_t *)((unsigned char *)created + rank_values_offset);
  created->component_buffer =
      (MVMCAbsolutePfaffianScaledValueResult *)((unsigned char *)created +
                                                 component_offset);
  created->projection_terms =
      (MVMCScaledComplex *)((unsigned char *)created +
                            projection_term_offset);
  created->rank_accumulator_buffer =
      (MVMCScaledComplex *)((unsigned char *)created +
                            rank_accumulator_offset);
  created->rank_order_scratch =
      (MVMCScaledComplex *)((unsigned char *)created + rank_scratch_offset);
  created->merge_buffer =
      (MVMCScaledComplex *)((unsigned char *)created + merge_offset);
#ifdef _mpi_use
  created->communicator = duplicated;
#endif
  *workspace = created;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_bounded_krylov_collective_workspace_destroy(
    MVMCKrylovBoundedCollectiveWorkspace *workspace) {
  if (workspace == NULL) return;
#ifdef _mpi_use
  if (workspace->communicator != MPI_COMM_NULL) {
    MPI_Comm_free(&workspace->communicator);
  }
#endif
  free(workspace);
}

size_t mvmc_bounded_krylov_collective_workspace_bytes(
    const MVMCKrylovBoundedCollectiveWorkspace *workspace) {
  return workspace == NULL ? 0 : workspace->allocated_bytes;
}

MVMCKrylovStatus mvmc_bounded_krylov_collective_synchronize(
    MVMCKrylovBoundedCollectiveWorkspace *workspace,
    MVMCKrylovStatus local_status, uint64_t local_plan_hash,
    int local_allocation_ok, MVMCKrylovBoundedCollectiveResult *result) {
  MVMCKrylovBoundedCollectiveResult discarded;
  int mismatch_rank = -1;
  int rank;
  MVMCKrylovStatus effective = local_status;
  const int result_available = result != NULL;
  if (!result_available) result = &discarded;
  reset_result(result);
  if (workspace == NULL || workspace->rank_values == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (local_allocation_ok != 0 && local_allocation_ok != 1) {
    effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    local_allocation_ok = 1;
  }
#ifdef _mpi_use
  if (MPI_Allgather(&local_plan_hash, 1, MPI_UINT64_T,
                    workspace->rank_values, 1, MPI_UINT64_T,
                    workspace->communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#else
  workspace->rank_values[0] = local_plan_hash;
#endif
  for (rank = 1; rank < workspace->size; ++rank) {
    if (workspace->rank_values[rank] != workspace->rank_values[0]) {
      mismatch_rank = rank;
      break;
    }
  }
  if (!valid_status(effective)) effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (!local_allocation_ok &&
      effective < MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE) {
    effective = MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
  if (mismatch_rank == workspace->rank) {
    effective = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  effective = synchronize_status(workspace, effective, result_available,
                                 result);
  result->plan_hash = workspace->rank_values[0];
  return effective;
}

MVMCKrylovStatus mvmc_bounded_krylov_collective_max_u64(
    MVMCKrylovBoundedCollectiveWorkspace *workspace,
    uint64_t local_value, uint64_t *global_maximum) {
  uint64_t maximum = local_value;
  int local_valid;
  int all_valid;
  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  local_valid = global_maximum != NULL;
#ifdef _mpi_use
  if (MPI_Allreduce(&local_valid, &all_valid, 1, MPI_INT, MPI_MIN,
                    workspace->communicator) != MPI_SUCCESS ||
      MPI_Allreduce(&local_value, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                    workspace->communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#else
  all_valid = local_valid;
#endif
  if (!all_valid) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *global_maximum = maximum;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_bounded_krylov_collective_gather_scaled_components(
    MVMCKrylovBoundedCollectiveWorkspace *workspace,
    MVMCKrylovStatus local_status,
    const MVMCAbsolutePfaffianScaledValueResult *local_components,
    size_t local_component_count, int qp_total, int qp_start, int qp_end,
    MVMCAbsolutePfaffianScaledValueResult *global_components,
    MVMCKrylovBoundedCollectiveResult *result) {
  MVMCKrylovBoundedCollectiveResult discarded;
  const int result_available = result != NULL;
  MVMCKrylovStatus effective = valid_status(local_status)
                                      ? local_status
                                      : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  int local_range[3] = {qp_total, qp_start, qp_end};
  int partition_failure_rank;
  size_t index;
  if (!result_available) result = &discarded;
  reset_result(result);
  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
#ifdef _mpi_use
  if (MPI_Allgather(local_range, 3, MPI_INT, workspace->ranges, 3, MPI_INT,
                    workspace->communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#else
  memcpy(workspace->ranges, local_range, sizeof(local_range));
#endif
  partition_failure_rank = first_invalid_partition(
      workspace->ranges, workspace->size, workspace->max_qp_components);
  if (partition_failure_rank == workspace->rank) {
    effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (qp_start < 0 || qp_end < qp_start ||
      local_component_count != (size_t)(qp_end - qp_start) ||
      (local_component_count != 0 && local_components == NULL) ||
      (qp_total != 0 && global_components == NULL)) {
    effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (effective == MVMC_KRYLOV_STATUS_OK) {
    for (index = 0; index < local_component_count; ++index) {
      if (!valid_scaled_component(&local_components[index], &effective)) break;
    }
  }
  effective = synchronize_status(workspace, effective, result_available,
                                 result);
  if (effective != MVMC_KRYLOV_STATUS_OK) return effective;
  for (index = 0; index < (size_t)workspace->size; ++index) {
    const int start = workspace->ranges[index * 3 + 1];
    const int end = workspace->ranges[index * 3 + 2];
    workspace->byte_counts[index] =
        (end - start) * (int)sizeof(*local_components);
    workspace->byte_displacements[index] =
        start * (int)sizeof(*local_components);
  }
#ifdef _mpi_use
  if (MPI_Allgatherv(
          local_component_count == 0 ? workspace->component_buffer
                                     : local_components,
          (int)(local_component_count * sizeof(*local_components)), MPI_BYTE,
          workspace->component_buffer, workspace->byte_counts,
          workspace->byte_displacements, MPI_BYTE,
          workspace->communicator) != MPI_SUCCESS) {
    reset_result(result);
    result->status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    return result->status;
  }
#else
  if (local_component_count != 0) {
    memcpy(workspace->component_buffer, local_components,
           local_component_count * sizeof(*local_components));
  }
#endif
  if (global_components != workspace->component_buffer) {
    memcpy(global_components, workspace->component_buffer,
           (size_t)qp_total * sizeof(*global_components));
  }
  result->valid = 1;
  result->status = MVMC_KRYLOV_STATUS_OK;
  result->failure_rank = -1;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus
mvmc_bounded_krylov_collective_gather_projected_amplitude(
    MVMCKrylovBoundedCollectiveWorkspace *workspace,
    MVMCKrylovStatus local_status,
    const MVMCAbsolutePfaffianScaledValueResult *local_components,
    size_t local_component_count, const double complex *global_weights,
    int qp_total, int qp_start, int qp_end,
    MVMCKrylovScaledAmplitudeResult *amplitude,
    MVMCKrylovBoundedCollectiveResult *result) {
  MVMCKrylovScaledAmplitudeResult candidate;
  MVMCKrylovBoundedCollectiveResult preflight;
  MVMCKrylovStatus effective = valid_status(local_status)
                                      ? local_status
                                      : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  uint64_t weight_hash = UINT64_C(1469598103934665603);
  size_t index;
  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (qp_total <= 0 || (size_t)qp_total > workspace->max_qp_components ||
      global_weights == NULL || amplitude == NULL) {
    effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (qp_total > 0 && global_weights != NULL) {
    hash_u64(&weight_hash, (uint64_t)(unsigned int)qp_total);
    for (index = 0; index < (size_t)qp_total; ++index) {
      if (!finite_complex(global_weights[index])) {
        effective = MVMC_KRYLOV_STATUS_NONFINITE;
      }
      hash_u64(&weight_hash, double_bits(creal(global_weights[index])));
      hash_u64(&weight_hash, double_bits(cimag(global_weights[index])));
    }
  }
  effective = mvmc_bounded_krylov_collective_synchronize(
      workspace, effective, weight_hash, 1, &preflight);
  if (effective != MVMC_KRYLOV_STATUS_OK) {
    if (result != NULL) *result = preflight;
    return effective;
  }
  effective = mvmc_bounded_krylov_collective_gather_scaled_components(
      workspace, MVMC_KRYLOV_STATUS_OK, local_components,
      local_component_count, qp_total, qp_start, qp_end,
      workspace->component_buffer, result);
  if (effective != MVMC_KRYLOV_STATUS_OK) return effective;
  memset(&candidate, 0, sizeof(candidate));
  candidate.local_factorization_count = (uint64_t)local_component_count;
  candidate.global_factorization_count = (uint64_t)(unsigned int)qp_total;
  for (index = 0; index < (size_t)qp_total; ++index) {
    const MVMCAbsolutePfaffianScaledValueResult *component =
        &workspace->component_buffer[index];
    MVMCScaledComplex scaled_weight;
    switch (component->factor_state) {
      case MVMC_PFAFFIAN_VALUE_WELL_PIVOTED:
        if (component->value.state !=
            MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
          effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
        } else {
          ++candidate.well_pivoted_component_count;
        }
        break;
      case MVMC_PFAFFIAN_VALUE_NEAR_PIVOT:
        if (component->value.state !=
                MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
            component->value.state != MVMC_SCALED_COMPLEX_NUMERIC_ZERO) {
          effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
        } else {
          ++candidate.near_pivot_component_count;
        }
        break;
      case MVMC_PFAFFIAN_VALUE_SINGULAR:
        if (component->value.state != MVMC_SCALED_COMPLEX_EXACT_ZERO) {
          effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
        }
        break;
      case MVMC_PFAFFIAN_VALUE_NONFINITE:
      case MVMC_PFAFFIAN_VALUE_INVALID:
        effective = MVMC_KRYLOV_STATUS_NONFINITE;
        break;
    }
    if (component->value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO) {
      ++candidate.exact_zero_component_count;
    } else if (component->value.state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO) {
      ++candidate.numeric_zero_component_count;
    }
    if (effective != MVMC_KRYLOV_STATUS_OK ||
        !exact_scaled_weight(global_weights[index], &scaled_weight) ||
        mvmc_scaled_complex_multiply(
            &scaled_weight, &component->value,
            &workspace->projection_terms[index]) !=
            MVMC_PFAFFIAN_STATUS_OK ||
        workspace->projection_terms[index].state ==
            MVMC_SCALED_COMPLEX_NONFINITE) {
      if (effective == MVMC_KRYLOV_STATUS_OK) {
        effective = MVMC_KRYLOV_STATUS_NONFINITE;
      }
      break;
    }
  }
  if (effective == MVMC_KRYLOV_STATUS_OK &&
      (mvmc_scaled_complex_sum_ordered(
           workspace->projection_terms, (size_t)qp_total,
           &candidate.value) != MVMC_PFAFFIAN_STATUS_OK ||
       candidate.value.state == MVMC_SCALED_COMPLEX_NONFINITE)) {
    effective = MVMC_KRYLOV_STATUS_NONFINITE;
  }
  if (effective != MVMC_KRYLOV_STATUS_OK) {
    if (result != NULL) {
      reset_result(result);
      result->status = effective;
    }
    return effective;
  }
  *amplitude = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_bounded_krylov_collective_merge_scaled_accumulators(
    MVMCKrylovBoundedCollectiveWorkspace *workspace,
    MVMCKrylovStatus local_status,
    const MVMCScaledComplex *local_values, size_t value_count,
    MVMCScaledComplex *global_values,
    MVMCKrylovBoundedCollectiveResult *result) {
  MVMCKrylovBoundedCollectiveResult discarded;
  const int result_available = result != NULL;
  MVMCKrylovStatus effective = valid_status(local_status)
                                      ? local_status
                                      : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  uint64_t local_count;
  size_t index;
  int rank;
  if (!result_available) result = &discarded;
  reset_result(result);
  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
#if SIZE_MAX > UINT64_MAX
  if (value_count > UINT64_MAX) effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
#endif
  local_count = effective == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT
                    ? 0
                    : (uint64_t)value_count;
#ifdef _mpi_use
  if (MPI_Allgather(&local_count, 1, MPI_UINT64_T, workspace->rank_values, 1,
                    MPI_UINT64_T, workspace->communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#else
  workspace->rank_values[0] = local_count;
#endif
  if (value_count > workspace->max_accumulator_count ||
      (value_count != 0 && (local_values == NULL || global_values == NULL))) {
    effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (rank = 1; rank < workspace->size; ++rank) {
    if (workspace->rank_values[rank] != workspace->rank_values[0]) {
      effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      break;
    }
  }
  if (effective == MVMC_KRYLOV_STATUS_OK) {
    for (index = 0; index < value_count; ++index) {
      if (!mvmc_scaled_complex_is_valid(&local_values[index]) ||
          local_values[index].state == MVMC_SCALED_COMPLEX_NONFINITE) {
        effective = MVMC_KRYLOV_STATUS_NONFINITE;
        break;
      }
    }
  }
  effective = synchronize_status(workspace, effective, result_available,
                                 result);
  if (effective != MVMC_KRYLOV_STATUS_OK) return effective;
#ifdef _mpi_use
  if (MPI_Allgather(
          value_count == 0 ? workspace->rank_accumulator_buffer : local_values,
          (int)(value_count * sizeof(*local_values)), MPI_BYTE,
          workspace->rank_accumulator_buffer,
          (int)(value_count * sizeof(*local_values)), MPI_BYTE,
          workspace->communicator) != MPI_SUCCESS) {
    reset_result(result);
    result->status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    return result->status;
  }
#else
  if (value_count != 0) {
    memcpy(workspace->rank_accumulator_buffer, local_values,
           value_count * sizeof(*local_values));
  }
#endif
  for (index = 0; index < value_count; ++index) {
    for (rank = 0; rank < workspace->size; ++rank) {
      workspace->rank_order_scratch[rank] =
          workspace->rank_accumulator_buffer[
              (size_t)rank * value_count + index];
    }
    if (mvmc_scaled_complex_sum_ordered(
            workspace->rank_order_scratch, (size_t)workspace->size,
            &workspace->merge_buffer[index]) != MVMC_PFAFFIAN_STATUS_OK ||
        workspace->merge_buffer[index].state ==
            MVMC_SCALED_COMPLEX_NONFINITE) {
      reset_result(result);
      result->status = MVMC_KRYLOV_STATUS_NONFINITE;
      return result->status;
    }
  }
  if (value_count != 0) {
    memcpy(global_values, workspace->merge_buffer,
           value_count * sizeof(*global_values));
  }
  result->valid = 1;
  result->status = MVMC_KRYLOV_STATUS_OK;
  result->failure_rank = -1;
  return MVMC_KRYLOV_STATUS_OK;
}
