#include "power_lanczos_block_collective.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "power_lanczos_block_collective.c requires the bounded Krylov core"
#endif

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct MVMCPowerLanczosBlockCollectiveWorkspace {
  int rank;
  int size;
  size_t max_block_count;
  size_t max_entries_per_block;
  size_t max_entry_count;
  size_t allocated_bytes;
  uint64_t *shape_buffer;
  uint64_t *rank_count_buffer;
  double complex *rank_sum_buffer;
  uint64_t *merge_counts;
  double complex *merge_sums;
  double complex *merge_compensation;
#ifdef _mpi_use
  MPI_Comm communicator;
#endif
};

static int checked_add(size_t left, size_t right, size_t *result) {
  if (result == NULL || left > SIZE_MAX - right) return 0;
  *result = left + right;
  return 1;
}

static int checked_multiply(size_t left, size_t right, size_t *result) {
  if (result == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *result = left * right;
  return 1;
}

static int reserve(size_t *offset, size_t count, size_t item_size,
                   size_t alignment, size_t *start) {
  size_t aligned;
  size_t bytes;
  if (offset == NULL || start == NULL || alignment == 0 ||
      (alignment & (alignment - 1)) != 0 ||
      *offset > SIZE_MAX - (alignment - 1)) {
    return 0;
  }
  aligned = (*offset + alignment - 1) & ~(alignment - 1);
  if (!checked_multiply(count, item_size, &bytes) ||
      !checked_add(aligned, bytes, offset)) {
    return 0;
  }
  *start = aligned;
  return 1;
}

static int valid_status(MVMCKrylovStatus status) {
  return status >= MVMC_KRYLOV_STATUS_OK &&
         status <= MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

static int finite_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static int ranges_overlap(const void *left, size_t left_size,
                          const void *right, size_t right_size) {
  const uintptr_t left_begin = (uintptr_t)left;
  const uintptr_t right_begin = (uintptr_t)right;
  if (left == NULL || right == NULL ||
      left_begin > UINTPTR_MAX - left_size ||
      right_begin > UINTPTR_MAX - right_size) {
    return 1;
  }
  return left_begin < right_begin + right_size &&
         right_begin < left_begin + left_size;
}

static void reset_result(MVMCPowerLanczosBlockCollectiveResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  result->failure_rank = -1;
  result->version = MVMC_POWER_LANCZOS_BLOCK_COLLECTIVE_VERSION;
}

static MVMCKrylovStatus synchronize_status(
    MVMCPowerLanczosBlockCollectiveWorkspace *workspace,
    MVMCKrylovStatus local_status, int result_available,
    MVMCPowerLanczosBlockCollectiveResult *result) {
  MVMCKrylovStatus effective = valid_status(local_status)
                                      ? local_status
                                      : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  int failure_rank = effective == MVMC_KRYLOV_STATUS_OK ? -1 : 0;
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
    } local, global;
    local.value = (int)effective;
    local.index = workspace->rank;
    if (MPI_Allreduce(&local, &global, 1, MPI_2INT, MPI_MAXLOC,
                      workspace->communicator) != MPI_SUCCESS) {
      if (result_available) {
        result->status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
      }
      return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
    effective = (MVMCKrylovStatus)global.value;
    failure_rank = effective == MVMC_KRYLOV_STATUS_OK ? -1 : global.index;
  }
#endif
  if (result_available) {
    result->status = effective;
    result->failure_rank = failure_rank;
    result->valid = effective == MVMC_KRYLOV_STATUS_OK;
  }
  return effective;
}

MVMCKrylovStatus mvmc_power_lanczos_block_collective_create(
    MVMCKrylovBoundedCommunicator communicator,
    size_t max_block_count, size_t max_entries_per_block,
    MVMCPowerLanczosBlockCollectiveWorkspace **workspace) {
  MVMCPowerLanczosBlockCollectiveWorkspace *created = NULL;
  int rank = 0;
  int size = 1;
  int local_valid = workspace != NULL && max_block_count > 0 &&
                    max_entries_per_block > 0;
  size_t entry_count = 0;
  size_t rank_block_count = 0;
  size_t rank_entry_count = 0;
  size_t shape_count = 0;
  size_t offset = sizeof(*created);
  size_t shape_offset;
  size_t rank_count_offset;
  size_t rank_sum_offset;
  size_t merge_count_offset;
  size_t merge_sum_offset;
  size_t compensation_offset;
#ifdef _mpi_use
  MPI_Comm duplicated = MPI_COMM_NULL;
  int initialized = 0;
  int all_valid = 0;
  int allocation_ok;
  int all_allocation_ok;
  uint64_t local_shape[2];
  uint64_t minimum_shape[2];
  uint64_t maximum_shape[2];
  if (workspace != NULL) *workspace = NULL;
  if (communicator == MPI_COMM_NULL ||
      MPI_Initialized(&initialized) != MPI_SUCCESS || !initialized ||
      MPI_Comm_dup(communicator, &duplicated) != MPI_SUCCESS ||
      MPI_Comm_set_errhandler(duplicated, MPI_ERRORS_RETURN) != MPI_SUCCESS ||
      MPI_Comm_rank(duplicated, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(duplicated, &size) != MPI_SUCCESS || size <= 0) {
    if (duplicated != MPI_COMM_NULL) MPI_Comm_free(&duplicated);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#if SIZE_MAX > UINT64_MAX
  if (max_block_count > UINT64_MAX ||
      max_entries_per_block > UINT64_MAX) {
    local_valid = 0;
  }
#endif
  local_shape[0] = local_valid ? (uint64_t)max_block_count : 0;
  local_shape[1] = local_valid ? (uint64_t)max_entries_per_block : 0;
  if (MPI_Allreduce(&local_valid, &all_valid, 1, MPI_INT, MPI_MIN,
                    duplicated) != MPI_SUCCESS ||
      MPI_Allreduce(local_shape, minimum_shape, 2, MPI_UINT64_T, MPI_MIN,
                    duplicated) != MPI_SUCCESS ||
      MPI_Allreduce(local_shape, maximum_shape, 2, MPI_UINT64_T, MPI_MAX,
                    duplicated) != MPI_SUCCESS) {
    MPI_Comm_free(&duplicated);
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (!all_valid ||
      memcmp(minimum_shape, maximum_shape, sizeof(minimum_shape)) != 0) {
    MPI_Comm_free(&duplicated);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#else
  (void)communicator;
  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (!local_valid) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
#endif
  if (!checked_multiply(max_block_count, max_entries_per_block,
                        &entry_count) ||
      !checked_multiply((size_t)size, max_block_count,
                        &rank_block_count) ||
      !checked_multiply((size_t)size, entry_count, &rank_entry_count) ||
      !checked_multiply((size_t)size, 2, &shape_count) ||
      max_block_count > (size_t)INT_MAX / sizeof(uint64_t) ||
      entry_count > (size_t)INT_MAX / sizeof(double complex) ||
      !reserve(&offset, shape_count, sizeof(uint64_t),
               _Alignof(uint64_t), &shape_offset) ||
      !reserve(&offset, rank_block_count, sizeof(uint64_t),
               _Alignof(uint64_t), &rank_count_offset) ||
      !reserve(&offset, rank_entry_count, sizeof(double complex),
               _Alignof(double complex), &rank_sum_offset) ||
      !reserve(&offset, max_block_count, sizeof(uint64_t),
               _Alignof(uint64_t), &merge_count_offset) ||
      !reserve(&offset, entry_count, sizeof(double complex),
               _Alignof(double complex), &merge_sum_offset) ||
      !reserve(&offset, entry_count, sizeof(double complex),
               _Alignof(double complex), &compensation_offset)) {
#ifdef _mpi_use
    MPI_Comm_free(&duplicated);
#endif
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  created = (MVMCPowerLanczosBlockCollectiveWorkspace *)calloc(1, offset);
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
  created->max_block_count = max_block_count;
  created->max_entries_per_block = max_entries_per_block;
  created->max_entry_count = entry_count;
  created->allocated_bytes = offset;
  created->shape_buffer =
      (uint64_t *)((unsigned char *)created + shape_offset);
  created->rank_count_buffer =
      (uint64_t *)((unsigned char *)created + rank_count_offset);
  created->rank_sum_buffer =
      (double complex *)((unsigned char *)created + rank_sum_offset);
  created->merge_counts =
      (uint64_t *)((unsigned char *)created + merge_count_offset);
  created->merge_sums =
      (double complex *)((unsigned char *)created + merge_sum_offset);
  created->merge_compensation =
      (double complex *)((unsigned char *)created + compensation_offset);
#ifdef _mpi_use
  created->communicator = duplicated;
#endif
  *workspace = created;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_power_lanczos_block_collective_destroy(
    MVMCPowerLanczosBlockCollectiveWorkspace *workspace) {
  if (workspace == NULL) return;
#ifdef _mpi_use
  if (workspace->communicator != MPI_COMM_NULL) {
    MPI_Comm_free(&workspace->communicator);
  }
#endif
  memset(workspace, 0, workspace->allocated_bytes);
  free(workspace);
}

size_t mvmc_power_lanczos_block_collective_allocated_bytes(
    const MVMCPowerLanczosBlockCollectiveWorkspace *workspace) {
  return workspace == NULL ? 0 : workspace->allocated_bytes;
}

static MVMCKrylovStatus compensated_add(
    double complex value, double complex *sum,
    double complex *compensation) {
  double real_sum = creal(*sum);
  double imag_sum = cimag(*sum);
  double real_compensation = creal(*compensation);
  double imag_compensation = cimag(*compensation);
  const double real_value = creal(value);
  const double imag_value = cimag(value);
  const double real_next = real_sum + real_value;
  const double imag_next = imag_sum + imag_value;
  if (!finite_complex(value) || !isfinite(real_next) ||
      !isfinite(imag_next)) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  if (fabs(real_sum) >= fabs(real_value)) {
    real_compensation += (real_sum - real_next) + real_value;
  } else {
    real_compensation += (real_value - real_next) + real_sum;
  }
  if (fabs(imag_sum) >= fabs(imag_value)) {
    imag_compensation += (imag_sum - imag_next) + imag_value;
  } else {
    imag_compensation += (imag_value - imag_next) + imag_sum;
  }
  if (!isfinite(real_compensation) || !isfinite(imag_compensation)) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  *sum = real_next + I * imag_next;
  *compensation = real_compensation + I * imag_compensation;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_block_collective_reduce(
    MVMCPowerLanczosBlockCollectiveWorkspace *workspace,
    MVMCKrylovStatus local_status, size_t block_count,
    size_t entries_per_block, const uint64_t *local_sample_counts,
    const double complex *local_sums,
    uint64_t *global_sample_counts, size_t global_count_capacity,
    double complex *global_sums, size_t global_sum_capacity,
    MVMCPowerLanczosBlockCollectiveResult *result) {
  MVMCKrylovStatus effective = valid_status(local_status)
                                      ? local_status
                                      : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  const int result_available = result != NULL;
  uint64_t local_shape[2] = {0, 0};
  size_t entry_count = 0;
  size_t count_bytes = 0;
  size_t sum_bytes = 0;
  size_t block;
  size_t entry;
  int mismatch_rank = -1;
  int rank;
  uint64_t local_block_length = 0;
  uint64_t global_block_length = 0;
  if (result_available) reset_result(result);
  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
#if SIZE_MAX > UINT64_MAX
  if (block_count > UINT64_MAX || entries_per_block > UINT64_MAX) {
    effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#endif
  local_shape[0] = (uint64_t)block_count;
  local_shape[1] = (uint64_t)entries_per_block;
#ifdef _mpi_use
  if (MPI_Allgather(local_shape, 2, MPI_UINT64_T, workspace->shape_buffer,
                    2, MPI_UINT64_T,
                    workspace->communicator) != MPI_SUCCESS) {
    if (result_available) {
      result->status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#else
  memcpy(workspace->shape_buffer, local_shape, sizeof(local_shape));
#endif
  if (block_count == 0 || entries_per_block == 0 ||
      block_count > workspace->max_block_count ||
      entries_per_block > workspace->max_entries_per_block ||
      !checked_multiply(block_count, entries_per_block, &entry_count) ||
      !checked_multiply(block_count, sizeof(uint64_t), &count_bytes) ||
      !checked_multiply(entry_count, sizeof(double complex), &sum_bytes) ||
      local_sample_counts == NULL || local_sums == NULL ||
      global_sample_counts == NULL || global_sums == NULL ||
      global_count_capacity != block_count ||
      global_sum_capacity != entry_count ||
      ranges_overlap(global_sample_counts, count_bytes,
                     global_sums, sum_bytes)) {
    effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (rank = 1; rank < workspace->size; ++rank) {
    const uint64_t *shape = workspace->shape_buffer + (size_t)rank * 2;
    if (shape[0] != workspace->shape_buffer[0] ||
        shape[1] != workspace->shape_buffer[1]) {
      mismatch_rank = rank;
      break;
    }
  }
  if (mismatch_rank == workspace->rank) {
    effective = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (effective == MVMC_KRYLOV_STATUS_OK) {
    local_block_length = local_sample_counts[0];
    if (local_block_length == 0 ||
        local_block_length >
            MVMC_POWER_LANCZOS_BLOCK_MAX_EXACT_SAMPLE_COUNT) {
      effective = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  }
  if (effective == MVMC_KRYLOV_STATUS_OK) {
    for (block = 0; block < block_count; ++block) {
      if (local_sample_counts[block] != local_block_length) {
        effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
        break;
      }
    }
  }
  if (effective == MVMC_KRYLOV_STATUS_OK) {
    for (entry = 0; entry < entry_count; ++entry) {
      if (!finite_complex(local_sums[entry])) {
        effective = MVMC_KRYLOV_STATUS_NONFINITE;
        break;
      }
    }
  }
  effective = synchronize_status(workspace, effective, result_available,
                                 result);
  if (effective != MVMC_KRYLOV_STATUS_OK) return effective;
#ifdef _mpi_use
  if (MPI_Allgather(local_sample_counts, (int)count_bytes, MPI_BYTE,
                    workspace->rank_count_buffer, (int)count_bytes,
                    MPI_BYTE, workspace->communicator) != MPI_SUCCESS ||
      MPI_Allgather(local_sums, (int)sum_bytes, MPI_BYTE,
                    workspace->rank_sum_buffer, (int)sum_bytes,
                    MPI_BYTE, workspace->communicator) != MPI_SUCCESS) {
    reset_result(result);
    result->status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    return result->status;
  }
#else
  memcpy(workspace->rank_count_buffer, local_sample_counts, count_bytes);
  memcpy(workspace->rank_sum_buffer, local_sums, sum_bytes);
#endif
  memset(workspace->merge_counts, 0, count_bytes);
  memset(workspace->merge_sums, 0, sum_bytes);
  memset(workspace->merge_compensation, 0, sum_bytes);
  local_block_length = workspace->rank_count_buffer[0];
  mismatch_rank = -1;
  for (rank = 1; rank < workspace->size && mismatch_rank < 0; ++rank) {
    const size_t rank_count_offset = (size_t)rank * block_count;
    for (block = 0; block < block_count; ++block) {
      if (workspace->rank_count_buffer[rank_count_offset + block] !=
          local_block_length) {
        mismatch_rank = rank;
        break;
      }
    }
  }
  if (mismatch_rank == workspace->rank) {
    effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (rank = 0; rank < workspace->size; ++rank) {
    const size_t rank_count_offset = (size_t)rank * block_count;
    const size_t rank_sum_offset = (size_t)rank * entry_count;
    for (block = 0; block < block_count; ++block) {
      const uint64_t count =
          workspace->rank_count_buffer[rank_count_offset + block];
      if (count > UINT64_MAX - workspace->merge_counts[block]) {
        effective = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
        break;
      }
      workspace->merge_counts[block] += count;
    }
    if (effective != MVMC_KRYLOV_STATUS_OK) break;
    for (entry = 0; entry < entry_count; ++entry) {
      effective = compensated_add(
          workspace->rank_sum_buffer[rank_sum_offset + entry],
          &workspace->merge_sums[entry],
          &workspace->merge_compensation[entry]);
      if (effective != MVMC_KRYLOV_STATUS_OK) break;
    }
    if (effective != MVMC_KRYLOV_STATUS_OK) break;
  }
  if (effective == MVMC_KRYLOV_STATUS_OK) {
    global_block_length = workspace->merge_counts[0];
    if (global_block_length >
            MVMC_POWER_LANCZOS_BLOCK_MAX_EXACT_SAMPLE_COUNT ||
        block_count > UINT64_MAX / global_block_length) {
      effective = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  }
  if (effective == MVMC_KRYLOV_STATUS_OK) {
    for (entry = 0; entry < entry_count; ++entry) {
      workspace->merge_sums[entry] +=
          workspace->merge_compensation[entry];
      if (!finite_complex(workspace->merge_sums[entry])) {
        effective = MVMC_KRYLOV_STATUS_NONFINITE;
        break;
      }
    }
  }
  effective = synchronize_status(workspace, effective, result_available,
                                 result);
  if (effective != MVMC_KRYLOV_STATUS_OK) return effective;
  memcpy(global_sample_counts, workspace->merge_counts, count_bytes);
  memcpy(global_sums, workspace->merge_sums, sum_bytes);
  result->valid = 1;
  result->status = MVMC_KRYLOV_STATUS_OK;
  result->failure_rank = -1;
  result->world_size = (size_t)workspace->size;
  result->block_count = block_count;
  result->entries_per_block = entries_per_block;
  result->local_block_length = local_block_length;
  result->global_block_length = global_block_length;
  result->global_sample_count = global_block_length * (uint64_t)block_count;
  return MVMC_KRYLOV_STATUS_OK;
}
