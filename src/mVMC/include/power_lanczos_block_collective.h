#ifndef MVMC_POWER_LANCZOS_BLOCK_COLLECTIVE_H
#define MVMC_POWER_LANCZOS_BLOCK_COLLECTIVE_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "bounded_krylov_collective.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_BLOCK_COLLECTIVE_VERSION UINT64_C(1)
#define MVMC_POWER_LANCZOS_BLOCK_MAX_EXACT_SAMPLE_COUNT \
  UINT64_C(9007199254740992)

typedef struct MVMCPowerLanczosBlockCollectiveWorkspace
    MVMCPowerLanczosBlockCollectiveWorkspace;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int failure_rank;
  uint64_t version;
  size_t world_size;
  size_t block_count;
  size_t entries_per_block;
  uint64_t local_block_length;
  uint64_t global_block_length;
  uint64_t global_sample_count;
} MVMCPowerLanczosBlockCollectiveResult;

/* Creation is collective and duplicates communicator. */
MVMCKrylovStatus mvmc_power_lanczos_block_collective_create(
    MVMCKrylovBoundedCommunicator communicator,
    size_t max_block_count, size_t max_entries_per_block,
    MVMCPowerLanczosBlockCollectiveWorkspace **workspace);

void mvmc_power_lanczos_block_collective_destroy(
    MVMCPowerLanczosBlockCollectiveWorkspace *workspace);

size_t mvmc_power_lanczos_block_collective_allocated_bytes(
    const MVMCPowerLanczosBlockCollectiveWorkspace *workspace);

/*
 * Every rank supplies the same block_count and entries_per_block.  Each
 * local block must have the same positive sample count on every rank.
 * Sums are gathered and merged block-major / entry-major / increasing-rank.
 * Output is identical on all ranks and remains untouched on failure.
 */
MVMCKrylovStatus mvmc_power_lanczos_block_collective_reduce(
    MVMCPowerLanczosBlockCollectiveWorkspace *workspace,
    MVMCKrylovStatus local_status, size_t block_count,
    size_t entries_per_block, const uint64_t *local_sample_counts,
    const double complex *local_sums,
    uint64_t *global_sample_counts, size_t global_count_capacity,
    double complex *global_sums, size_t global_sum_capacity,
    MVMCPowerLanczosBlockCollectiveResult *result);

#endif

#endif
