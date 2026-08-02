/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "classic_pfaffian_collective.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct MVMCClassicPfaffianCollectiveWorkspace {
  int rank;
  int size;
  int *ranges;
#ifdef _mpi_use
  MPI_Comm communicator;
#endif
};

static int valid_status(MVMCPfaffianStatus status) {
  return status >= MVMC_PFAFFIAN_STATUS_OK &&
         status <= MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE;
}

#ifdef _mpi_use
static int status_severity(MVMCPfaffianStatus status) {
  switch (status) {
    case MVMC_PFAFFIAN_STATUS_OK:
      return 0;
    case MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT:
      return 1;
    case MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE:
      return 2;
    case MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE:
      return 3;
    case MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE:
      return 4;
    default:
      return 1;
  }
}
#endif

static int valid_local_aggregate(
    const MVMCProjectedAmplitudeResult *aggregate, size_t expected_count) {
  size_t counted;

  if (aggregate == NULL || aggregate->valid != 1 ||
      !isfinite(creal(aggregate->total)) ||
      !isfinite(cimag(aggregate->total)) ||
      !isfinite(aggregate->sum_abs) || aggregate->sum_abs < 0.0 ||
      !isfinite(aggregate->cancellation_ratio) ||
      aggregate->cancellation_ratio < 0.0 ||
      aggregate->cancellation_ratio > 1.0 ||
      aggregate->regular_count > expected_count ||
      aggregate->near_singular_count > expected_count ||
      aggregate->singular_count > expected_count) {
    return 0;
  }
  counted = aggregate->regular_count + aggregate->near_singular_count;
  if (counted < aggregate->regular_count ||
      aggregate->singular_count > SIZE_MAX - counted) {
    return 0;
  }
  return counted + aggregate->singular_count == expected_count;
}

static void reset_collective_result(
    MVMCClassicPfaffianCollectiveResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  result->failure_rank = -1;
  result->failure_local_qp = -1;
  result->failure_global_qp = -1;
}

static int first_invalid_partition(const int *ranges, int size) {
  /*
   * Rank 0 supplies the comparison qp_total.  The returned blame rank is
   * therefore the first rank that disagrees with rank 0, not necessarily the
   * rank where an erroneous caller value originated.
   */
  const int qp_total = ranges[0];
  int rank;

  if (qp_total <= 0) return 0;
  for (rank = 0; rank < size; ++rank) {
    const int *range = ranges + (size_t)rank * 3;
    if (range[0] != qp_total || range[1] < 0 ||
        range[2] < range[1] || range[2] > qp_total) {
      return rank;
    }
    if ((rank == 0 && range[1] != 0) ||
        (rank > 0 && range[1] != ranges[(size_t)(rank - 1) * 3 + 2])) {
      return rank;
    }
  }
  return ranges[(size_t)(size - 1) * 3 + 2] == qp_total ? -1 : size - 1;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_collective_workspace_create(
    MVMCClassicPfaffianCommunicator communicator,
    MVMCClassicPfaffianCollectiveWorkspace **workspace) {
  MVMCClassicPfaffianCollectiveWorkspace *created = NULL;
  int size = 1;
#ifdef _mpi_use
  MPI_Comm duplicated = MPI_COMM_NULL;
  int initialized = 0;
  int allocation_ok, all_allocation_ok;

  if (workspace == NULL || communicator == MPI_COMM_NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  *workspace = NULL;
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  if (MPI_Initialized(&initialized) != MPI_SUCCESS || !initialized ||
      MPI_Comm_dup(communicator, &duplicated) != MPI_SUCCESS ||
      MPI_Comm_size(duplicated, &size) != MPI_SUCCESS || size <= 0) {
    if (duplicated != MPI_COMM_NULL) MPI_Comm_free(&duplicated);
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  created = (MVMCClassicPfaffianCollectiveWorkspace *)calloc(
      1, sizeof(*created));
  if (created != NULL &&
      (size_t)size <= SIZE_MAX / (3 * sizeof(*created->ranges))) {
    created->ranges = (int *)calloc((size_t)size * 3,
                                    sizeof(*created->ranges));
  }
  allocation_ok = created != NULL && created->ranges != NULL;
  if (MPI_Allreduce(&allocation_ok, &all_allocation_ok, 1, MPI_INT, MPI_MIN,
                    duplicated) != MPI_SUCCESS) {
    all_allocation_ok = 0;
  }
  if (!all_allocation_ok || created == NULL || created->ranges == NULL) {
    if (created != NULL) {
      free(created->ranges);
      free(created);
    }
    MPI_Comm_free(&duplicated);
    return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  }
  if (MPI_Comm_rank(duplicated, &created->rank) != MPI_SUCCESS) {
    free(created->ranges);
    free(created);
    MPI_Comm_free(&duplicated);
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  created->size = size;
  created->communicator = duplicated;
#else
  (void)communicator;
  if (workspace == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    return MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE;
  }
  created = (MVMCClassicPfaffianCollectiveWorkspace *)calloc(
      1, sizeof(*created));
  if (created != NULL) {
    created->ranges = (int *)calloc(3, sizeof(*created->ranges));
  }
  if (created == NULL || created->ranges == NULL) {
    if (created != NULL) free(created);
    return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
  }
  created->rank = 0;
  created->size = size;
#endif
  *workspace = created;
  return MVMC_PFAFFIAN_STATUS_OK;
}

void mvmc_classic_pfaffian_collective_workspace_destroy(
    MVMCClassicPfaffianCollectiveWorkspace *workspace) {
  if (workspace == NULL) return;
#ifdef _mpi_use
  if (workspace->communicator != MPI_COMM_NULL) {
    MPI_Comm_free(&workspace->communicator);
  }
#endif
  free(workspace->ranges);
  free(workspace);
}

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)
size_t mvmc_classic_pfaffian_collective_workspace_bytes(
    const MVMCClassicPfaffianCollectiveWorkspace *workspace) {
  if (workspace == NULL || workspace->size <= 0) return 0;
  return sizeof(*workspace) +
         (size_t)workspace->size * 3 * sizeof(*workspace->ranges);
}
#endif

MVMCPfaffianStatus mvmc_classic_pfaffian_collective_aggregate(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    MVMCPfaffianStatus local_status,
    const MVMCProjectedAmplitudeResult *local_aggregate,
    int qp_total, int qp_start, int qp_end,
    int failure_local_qp, int failure_global_qp,
    MVMCClassicPfaffianCollectiveResult *result) {
  MVMCClassicPfaffianCollectiveResult discarded_result;
  MVMCPfaffianStatus effective_status;
  const int result_available = result != NULL;
  int local_range[3] = {qp_total, qp_start, qp_end};
  int partition_failure_rank;
  int selected_rank;
  int failure_payload[3] = {
      (int)MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT, -1, -1};
  double complex local_total = 0.0;
  double complex global_total = 0.0;
  double local_sum_abs = 0.0;
  double global_sum_abs = 0.0;
  uint64_t local_counts[3] = {0, 0, 0};
  uint64_t global_counts[3] = {0, 0, 0};
  size_t expected_count = 0;
#ifdef _mpi_use
  struct {
    int value;
    int index;
  } local_selection, global_selection;
#endif

  if (!result_available) result = &discarded_result;
  reset_collective_result(result);
  if (workspace == NULL || workspace->ranges == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }

#ifdef _mpi_use
  if (MPI_Allgather(local_range, 3, MPI_INT, workspace->ranges, 3, MPI_INT,
                    workspace->communicator) != MPI_SUCCESS) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
#else
  memcpy(workspace->ranges, local_range, sizeof(local_range));
#endif
  partition_failure_rank =
      first_invalid_partition(workspace->ranges, workspace->size);

  effective_status = valid_status(local_status)
                         ? local_status
                         : MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  if (!result_available &&
      effective_status == MVMC_PFAFFIAN_STATUS_OK) {
    effective_status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    failure_local_qp = -1;
    failure_global_qp = -1;
  }
  if (partition_failure_rank == workspace->rank &&
      effective_status == MVMC_PFAFFIAN_STATUS_OK) {
    effective_status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    failure_local_qp = -1;
    failure_global_qp = -1;
  }
  if (qp_start >= 0 && qp_end >= qp_start) {
    expected_count = (size_t)(qp_end - qp_start);
  }
  if (effective_status == MVMC_PFAFFIAN_STATUS_OK &&
      local_aggregate == NULL) {
    effective_status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    failure_local_qp = -1;
    failure_global_qp = -1;
  }
  if (effective_status == MVMC_PFAFFIAN_STATUS_OK &&
      !valid_local_aggregate(local_aggregate, expected_count)) {
    effective_status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    failure_local_qp = -1;
    failure_global_qp = -1;
  }
#if SIZE_MAX > UINT64_MAX
  if (effective_status == MVMC_PFAFFIAN_STATUS_OK &&
      (local_aggregate->regular_count > UINT64_MAX ||
       local_aggregate->near_singular_count > UINT64_MAX ||
       local_aggregate->singular_count > UINT64_MAX)) {
    effective_status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    failure_local_qp = -1;
    failure_global_qp = -1;
  }
#endif

#ifdef _mpi_use
  local_selection.value = status_severity(effective_status);
  local_selection.index = workspace->rank;
  if (MPI_Allreduce(&local_selection, &global_selection, 1, MPI_2INT,
                    MPI_MAXLOC, workspace->communicator) != MPI_SUCCESS) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  selected_rank = global_selection.index;
#else
  selected_rank = 0;
#endif
  if (workspace->rank == selected_rank) {
    failure_payload[0] = (int)effective_status;
    failure_payload[1] = effective_status == MVMC_PFAFFIAN_STATUS_OK
                             ? -1
                             : failure_local_qp;
    failure_payload[2] = effective_status == MVMC_PFAFFIAN_STATUS_OK
                             ? -1
                             : failure_global_qp;
  }
#ifdef _mpi_use
  if (MPI_Bcast(failure_payload, 3, MPI_INT, selected_rank,
                workspace->communicator) != MPI_SUCCESS) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
#endif

  if (effective_status == MVMC_PFAFFIAN_STATUS_OK) {
    local_total = local_aggregate->total;
    local_sum_abs = local_aggregate->sum_abs;
    local_counts[0] = (uint64_t)local_aggregate->regular_count;
    local_counts[1] = (uint64_t)local_aggregate->near_singular_count;
    local_counts[2] = (uint64_t)local_aggregate->singular_count;
  }
#ifdef _mpi_use
  if (MPI_Allreduce(&local_total, &global_total, 1, MPI_C_DOUBLE_COMPLEX,
                    MPI_SUM, workspace->communicator) != MPI_SUCCESS ||
      MPI_Allreduce(&local_sum_abs, &global_sum_abs, 1, MPI_DOUBLE, MPI_SUM,
                    workspace->communicator) != MPI_SUCCESS ||
      MPI_Allreduce(local_counts, global_counts, 3, MPI_UINT64_T, MPI_SUM,
                    workspace->communicator) != MPI_SUCCESS) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
#else
  global_total = local_total;
  global_sum_abs = local_sum_abs;
  memcpy(global_counts, local_counts, sizeof(global_counts));
#endif

  result->status = (MVMCPfaffianStatus)failure_payload[0];
  if (result->status != MVMC_PFAFFIAN_STATUS_OK) {
    result->failure_rank = selected_rank;
    result->failure_local_qp = failure_payload[1];
    result->failure_global_qp = failure_payload[2];
    return result->status;
  }
  if (!isfinite(creal(global_total)) || !isfinite(cimag(global_total)) ||
      !isfinite(global_sum_abs) || global_sum_abs < 0.0 ||
      global_counts[0] + global_counts[1] < global_counts[0] ||
      global_counts[2] > UINT64_MAX -
                             (global_counts[0] + global_counts[1]) ||
      global_counts[0] + global_counts[1] + global_counts[2] !=
          (uint64_t)workspace->ranges[0]) {
    result->status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    return result->status;
  }
#if SIZE_MAX < UINT64_MAX
  if (global_counts[0] > (uint64_t)SIZE_MAX ||
      global_counts[1] > (uint64_t)SIZE_MAX ||
      global_counts[2] > (uint64_t)SIZE_MAX) {
    result->status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    return result->status;
  }
#endif

  result->aggregate.total = global_total;
  result->aggregate.sum_abs = global_sum_abs;
  result->aggregate.cancellation_ratio =
      global_sum_abs == 0.0
          ? 0.0
          : fmin(1.0, cabs(global_total) / global_sum_abs);
  result->aggregate.regular_count = (size_t)global_counts[0];
  result->aggregate.near_singular_count = (size_t)global_counts[1];
  result->aggregate.singular_count = (size_t)global_counts[2];
  result->aggregate.valid = 1;
  result->valid = 1;
  result->status = MVMC_PFAFFIAN_STATUS_OK;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_real_prepare_collective(
    MVMCClassicPfaffianRealWorkspace *state_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const double *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCClassicPfaffianCollectiveResult *result) {
  MVMCPfaffianStatus local_status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  MVMCPfaffianStatus global_status;
  const MVMCClassicPfaffianState *candidate = NULL;
  int qp_total = -1, qp_start = -1, qp_end = -1;
  int failure_local_qp = -1, failure_global_qp = -1;

  if (state_workspace != NULL &&
      mvmc_classic_pfaffian_real_qp_range(
          state_workspace, &qp_total, &qp_start, &qp_end)) {
    local_status = mvmc_classic_pfaffian_real_prepare(
        state_workspace, slater_elm, ele_idx, global_weights,
        scaled_pivot_tolerance, residual_tolerance);
    failure_local_qp =
        mvmc_classic_pfaffian_real_failure_local_qp(state_workspace);
    failure_global_qp =
        mvmc_classic_pfaffian_real_failure_global_qp(state_workspace);
    if (local_status == MVMC_PFAFFIAN_STATUS_OK) {
      candidate = mvmc_classic_pfaffian_real_candidate(state_workspace);
    }
  }
  global_status = mvmc_classic_pfaffian_collective_aggregate(
      collective_workspace, local_status,
      candidate == NULL ? NULL : &candidate->local_aggregate,
      qp_total, qp_start, qp_end, failure_local_qp, failure_global_qp,
      result);
  if (global_status != MVMC_PFAFFIAN_STATUS_OK && state_workspace != NULL) {
    mvmc_classic_pfaffian_real_discard_candidate(state_workspace);
  }
  return global_status;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_complex_prepare_collective(
    MVMCClassicPfaffianComplexWorkspace *state_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const double complex *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCClassicPfaffianCollectiveResult *result) {
  MVMCPfaffianStatus local_status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  MVMCPfaffianStatus global_status;
  const MVMCClassicPfaffianState *candidate = NULL;
  int qp_total = -1, qp_start = -1, qp_end = -1;
  int failure_local_qp = -1, failure_global_qp = -1;

  if (state_workspace != NULL &&
      mvmc_classic_pfaffian_complex_qp_range(
          state_workspace, &qp_total, &qp_start, &qp_end)) {
    local_status = mvmc_classic_pfaffian_complex_prepare(
        state_workspace, slater_elm, ele_idx, global_weights,
        scaled_pivot_tolerance, residual_tolerance);
    failure_local_qp =
        mvmc_classic_pfaffian_complex_failure_local_qp(state_workspace);
    failure_global_qp =
        mvmc_classic_pfaffian_complex_failure_global_qp(state_workspace);
    if (local_status == MVMC_PFAFFIAN_STATUS_OK) {
      candidate = mvmc_classic_pfaffian_complex_candidate(state_workspace);
    }
  }
  global_status = mvmc_classic_pfaffian_collective_aggregate(
      collective_workspace, local_status,
      candidate == NULL ? NULL : &candidate->local_aggregate,
      qp_total, qp_start, qp_end, failure_local_qp, failure_global_qp,
      result);
  if (global_status != MVMC_PFAFFIAN_STATUS_OK && state_workspace != NULL) {
    mvmc_classic_pfaffian_complex_discard_candidate(state_workspace);
  }
  return global_status;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_collective_metropolis(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    double complex current_total, double complex proposal_total,
    double log_non_pfaffian_ratio, double uniform,
    MVMCClassicPfaffianMetropolisResult *result) {
  MVMCClassicPfaffianMetropolisResult discarded_result;
  double values[6] = {creal(current_total), cimag(current_total),
                      creal(proposal_total), cimag(proposal_total),
                      log_non_pfaffian_ratio, uniform};
  double minimum[6];
  double maximum[6];
  int local_valid;
  int global_valid;
  int i;
  double log_acceptance;

  if (result == NULL) result = &discarded_result;
  memset(result, 0, sizeof(*result));
  result->status = MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  result->log_acceptance_ratio = -INFINITY;
  if (workspace == NULL || workspace->ranges == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }

  local_valid = isfinite(values[0]) && isfinite(values[1]) &&
                isfinite(values[2]) && isfinite(values[3]) &&
                !isnan(log_non_pfaffian_ratio) && isfinite(uniform) &&
                uniform >= 0.0 && uniform < 1.0 &&
                cabs(current_total) != 0.0 &&
                result != &discarded_result;
#ifdef _mpi_use
  if (MPI_Allreduce(&local_valid, &global_valid, 1, MPI_INT, MPI_MIN,
                    workspace->communicator) != MPI_SUCCESS ||
      MPI_Allreduce(values, minimum, 6, MPI_DOUBLE, MPI_MIN,
                    workspace->communicator) != MPI_SUCCESS ||
      MPI_Allreduce(values, maximum, 6, MPI_DOUBLE, MPI_MAX,
                    workspace->communicator) != MPI_SUCCESS) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
#else
  global_valid = local_valid;
  memcpy(minimum, values, sizeof(minimum));
  memcpy(maximum, values, sizeof(maximum));
#endif
  for (i = 0; i < 6; ++i) {
    if (minimum[i] != maximum[i]) global_valid = 0;
  }
  if (!global_valid) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;

  if (cabs(proposal_total) == 0.0) {
    result->valid = 1;
    result->status = MVMC_PFAFFIAN_STATUS_OK;
    return MVMC_PFAFFIAN_STATUS_OK;
  }
  log_acceptance = 2.0 *
      (log_non_pfaffian_ratio + log(cabs(proposal_total)) -
       log(cabs(current_total)));
  if (isnan(log_acceptance)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  result->accepted = log_acceptance != -INFINITY &&
      (uniform == 0.0 ||
       log(uniform) < fmin(0.0, log_acceptance));
  result->log_acceptance_ratio = log_acceptance;
  result->valid = 1;
  result->status = MVMC_PFAFFIAN_STATUS_OK;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_collective_preflight(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    int local_valid, int branch_kind) {
  int global_valid;
  int minimum_kind;
  int maximum_kind;

  if (workspace == NULL || workspace->ranges == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  local_valid = local_valid == 1;
#ifdef _mpi_use
  if (MPI_Allreduce(&local_valid, &global_valid, 1, MPI_INT, MPI_MIN,
                    workspace->communicator) != MPI_SUCCESS ||
      MPI_Allreduce(&branch_kind, &minimum_kind, 1, MPI_INT, MPI_MIN,
                    workspace->communicator) != MPI_SUCCESS ||
      MPI_Allreduce(&branch_kind, &maximum_kind, 1, MPI_INT, MPI_MAX,
                    workspace->communicator) != MPI_SUCCESS) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
#else
  global_valid = local_valid;
  minimum_kind = branch_kind;
  maximum_kind = branch_kind;
#endif
  return global_valid && minimum_kind == maximum_kind
             ? MVMC_PFAFFIAN_STATUS_OK
             : MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_collective_all_true(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    int local_true, int *all_true) {
  int local_valid;
  int global_valid;
  int global_true;

  if (workspace == NULL || workspace->ranges == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  local_valid = (local_true == 0 || local_true == 1) && all_true != NULL;
  if (!local_valid) local_true = 0;
#ifdef _mpi_use
  if (MPI_Allreduce(&local_valid, &global_valid, 1, MPI_INT, MPI_MIN,
                    workspace->communicator) != MPI_SUCCESS ||
      MPI_Allreduce(&local_true, &global_true, 1, MPI_INT, MPI_MIN,
                    workspace->communicator) != MPI_SUCCESS) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
#else
  global_valid = local_valid;
  global_true = local_true;
#endif
  if (!global_valid || all_true == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  *all_true = global_true;
  return MVMC_PFAFFIAN_STATUS_OK;
}

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)
MVMCPfaffianStatus mvmc_classic_pfaffian_collective_all_equal_bytes(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    const void *local_data, size_t byte_count, int *all_equal) {
  const unsigned char *bytes = (const unsigned char *)local_data;
  int local_valid;
  int global_valid;

  if (workspace == NULL || workspace->ranges == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  local_valid = all_equal != NULL && (byte_count == 0 || local_data != NULL);
#ifdef _mpi_use
  {
    uint64_t local_count = (uint64_t)byte_count;
    uint64_t minimum_count = 0;
    uint64_t maximum_count = 0;
    unsigned char *minimum = NULL;
    unsigned char *maximum = NULL;
    const size_t chunk_limit = 1024 * 1024;
    size_t buffer_size;
    size_t offset;
    int allocation_ok;
    int all_allocation_ok;
    int equal = 1;

#if SIZE_MAX > UINT64_MAX
    if (byte_count > UINT64_MAX) local_valid = 0;
#endif
    if (MPI_Allreduce(&local_valid, &global_valid, 1, MPI_INT, MPI_MIN,
                      workspace->communicator) != MPI_SUCCESS ||
        MPI_Allreduce(&local_count, &minimum_count, 1, MPI_UINT64_T, MPI_MIN,
                      workspace->communicator) != MPI_SUCCESS ||
        MPI_Allreduce(&local_count, &maximum_count, 1, MPI_UINT64_T, MPI_MAX,
                      workspace->communicator) != MPI_SUCCESS) {
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    if (!global_valid || minimum_count != maximum_count) {
      if (all_equal != NULL) *all_equal = 0;
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    buffer_size = byte_count < chunk_limit ? byte_count : chunk_limit;
    if (buffer_size != 0) {
      minimum = (unsigned char *)malloc(buffer_size);
      maximum = (unsigned char *)malloc(buffer_size);
    }
    allocation_ok = buffer_size == 0 ||
                    (minimum != NULL && maximum != NULL);
    if (MPI_Allreduce(&allocation_ok, &all_allocation_ok, 1, MPI_INT, MPI_MIN,
                      workspace->communicator) != MPI_SUCCESS) {
      free(maximum);
      free(minimum);
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    if (!all_allocation_ok ||
        (buffer_size != 0 && (minimum == NULL || maximum == NULL))) {
      free(maximum);
      free(minimum);
      return MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE;
    }
    for (offset = 0; offset < byte_count; offset += buffer_size) {
      const size_t remaining = byte_count - offset;
      const int count = (int)(remaining < buffer_size ? remaining
                                                       : buffer_size);
      if (MPI_Allreduce(bytes + offset, minimum, count, MPI_UNSIGNED_CHAR,
                        MPI_MIN, workspace->communicator) != MPI_SUCCESS ||
          MPI_Allreduce(bytes + offset, maximum, count, MPI_UNSIGNED_CHAR,
                        MPI_MAX, workspace->communicator) != MPI_SUCCESS) {
        free(maximum);
        free(minimum);
        return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
      }
      if (memcmp(minimum, maximum, (size_t)count) != 0) equal = 0;
    }
    free(maximum);
    free(minimum);
    if (all_equal == NULL) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    *all_equal = equal;
  }
#else
  global_valid = local_valid;
  if (!global_valid || all_equal == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  (void)bytes;
  *all_equal = 1;
#endif
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_collective_all_equal_u64(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    const uint64_t *local_values, size_t value_count, int *all_equal) {
  int local_valid;
  int global_valid;
  size_t index;
  int equal = 1;

  if (workspace == NULL || workspace->ranges == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  local_valid = all_equal != NULL &&
                (value_count == 0 || local_values != NULL);
#ifdef _mpi_use
  {
    uint64_t local_count = (uint64_t)value_count;
    uint64_t minimum_count = 0;
    uint64_t maximum_count = 0;
#if SIZE_MAX > UINT64_MAX
    if (value_count > UINT64_MAX) local_valid = 0;
#endif
    if (MPI_Allreduce(&local_valid, &global_valid, 1, MPI_INT, MPI_MIN,
                      workspace->communicator) != MPI_SUCCESS ||
        MPI_Allreduce(&local_count, &minimum_count, 1, MPI_UINT64_T, MPI_MIN,
                      workspace->communicator) != MPI_SUCCESS ||
        MPI_Allreduce(&local_count, &maximum_count, 1, MPI_UINT64_T, MPI_MAX,
                      workspace->communicator) != MPI_SUCCESS) {
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    if (!global_valid || minimum_count != maximum_count) {
      if (all_equal != NULL) *all_equal = 0;
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    for (index = 0; index < value_count; ++index) {
      uint64_t minimum = 0;
      uint64_t maximum = 0;
      if (MPI_Allreduce(local_values + index, &minimum, 1, MPI_UINT64_T,
                        MPI_MIN, workspace->communicator) != MPI_SUCCESS ||
          MPI_Allreduce(local_values + index, &maximum, 1, MPI_UINT64_T,
                        MPI_MAX, workspace->communicator) != MPI_SUCCESS) {
        return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
      }
      if (minimum != maximum) equal = 0;
    }
  }
#else
  global_valid = local_valid;
  (void)local_values;
  (void)value_count;
#endif
  if (!global_valid || all_equal == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  *all_equal = equal;
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_collective_sum_u64(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    uint64_t local_value, uint64_t *global_sum) {
  if (workspace == NULL || workspace->ranges == NULL || global_sum == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
#ifdef _mpi_use
  if (MPI_Allreduce(&local_value, global_sum, 1, MPI_UINT64_T, MPI_SUM,
                    workspace->communicator) != MPI_SUCCESS) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
#else
  *global_sum = local_value;
#endif
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus mvmc_classic_pfaffian_collective_max_int(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    int local_value, int *global_max) {
  if (workspace == NULL || workspace->ranges == NULL || global_max == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
#ifdef _mpi_use
  if (MPI_Allreduce(&local_value, global_max, 1, MPI_INT, MPI_MAX,
                    workspace->communicator) != MPI_SUCCESS) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
#else
  *global_max = local_value;
#endif
  return MVMC_PFAFFIAN_STATUS_OK;
}

MVMCPfaffianStatus
mvmc_classic_pfaffian_collective_gather_value_components(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    const MVMCAbsolutePfaffianValueResult *local_components,
    size_t local_component_count, int qp_total, int qp_start, int qp_end,
    MVMCAbsolutePfaffianValueResult *global_components) {
  int local_valid;
  int global_valid;
  int global_qp;

  if (workspace == NULL || workspace->ranges == NULL) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
  local_valid = qp_total > 0 && qp_start >= 0 && qp_end >= qp_start &&
                qp_end <= qp_total &&
                (size_t)(qp_end - qp_start) == local_component_count &&
                (local_component_count == 0 || local_components != NULL) &&
                global_components != NULL &&
                workspace->ranges[(size_t)workspace->rank * 3] == qp_total &&
                workspace->ranges[(size_t)workspace->rank * 3 + 1] ==
                    qp_start &&
                workspace->ranges[(size_t)workspace->rank * 3 + 2] == qp_end;
#ifdef _mpi_use
  if (MPI_Allreduce(&local_valid, &global_valid, 1, MPI_INT, MPI_MIN,
                    workspace->communicator) != MPI_SUCCESS) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }
#else
  global_valid = local_valid;
#endif
  if (!global_valid || global_components == NULL ||
      (local_component_count != 0 && local_components == NULL)) {
    return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  }

  for (global_qp = 0; global_qp < qp_total; ++global_qp) {
    int owner = -1;
    int rank;
    for (rank = 0; rank < workspace->size; ++rank) {
      const int *range = workspace->ranges + (size_t)rank * 3;
      if (global_qp >= range[1] && global_qp < range[2]) {
        owner = rank;
        break;
      }
    }
    if (owner < 0) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    if (workspace->rank == owner) {
      const size_t local_index = (size_t)(global_qp - qp_start);
      if (local_components != NULL) {
        global_components[global_qp] = local_components[local_index];
      }
    }
#ifdef _mpi_use
    if (MPI_Bcast(global_components + global_qp,
                  (int)sizeof(*global_components), MPI_BYTE, owner,
                  workspace->communicator) != MPI_SUCCESS) {
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
#endif
  }
  return MVMC_PFAFFIAN_STATUS_OK;
}
#endif /* MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE */
