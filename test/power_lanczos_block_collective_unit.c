#include "power_lanczos_block_collective.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failures = 0;
static int world_rank = 0;
static int world_size = 1;

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "PowerLanczosBlockCollective_Unit FAIL rank %d: ",      \
              world_rank);                                                      \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, " (line %d)\n", __LINE__);                              \
      ++failures;                                                               \
    }                                                                           \
  } while (0)

static void fill_valid(uint64_t counts[2], double complex sums[6]) {
  size_t block;
  size_t entry;
  const double rank_factor = (double)(world_rank + 1);
  counts[0] = 4;
  counts[1] = 4;
  for (block = 0; block < 2; ++block) {
    for (entry = 0; entry < 3; ++entry) {
      const double base = (double)(block * 10 + entry + 1);
      sums[block * 3 + entry] =
          rank_factor * base - I * rank_factor * base * 0.5;
    }
  }
}

static void fill_sentinel(uint64_t counts[2], double complex sums[6]) {
  size_t index;
  counts[0] = UINT64_C(0x1111222233334444);
  counts[1] = UINT64_C(0x5555666677778888);
  for (index = 0; index < 6; ++index) {
    sums[index] = (1000.0 + (double)index) + I * (2000.0 + (double)index);
  }
}

static int outputs_unchanged(const uint64_t counts[2],
                             const double complex sums[6],
                             const uint64_t expected_counts[2],
                             const double complex expected_sums[6]) {
  return memcmp(counts, expected_counts, 2 * sizeof(*counts)) == 0 &&
         memcmp(sums, expected_sums, 6 * sizeof(*sums)) == 0;
}

static void test_workspace_creation(void) {
  MVMCPowerLanczosBlockCollectiveWorkspace *workspace = NULL;
  MVMCKrylovStatus status;
#ifdef _mpi_use
  status = mvmc_power_lanczos_block_collective_create(
      MPI_COMM_WORLD,
      world_rank == world_size - 1 && world_size > 1 ? 3U : 4U,
      3, &workspace);
  if (world_size > 1) {
    CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
          "rank-mismatched creation shape was accepted");
  } else {
    CHECK(status == MVMC_KRYLOV_STATUS_OK && workspace != NULL,
          "single-rank creation fixture failed");
    mvmc_power_lanczos_block_collective_destroy(workspace);
    workspace = NULL;
  }
  status = mvmc_power_lanczos_block_collective_create(
      MPI_COMM_WORLD, 0, 3, &workspace);
#else
  status = mvmc_power_lanczos_block_collective_create(0, 0, 3, &workspace);
#endif
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
        "zero block capacity was accepted");
}

static void test_valid_reduce(
    MVMCPowerLanczosBlockCollectiveWorkspace *workspace) {
  uint64_t local_counts[2];
  uint64_t global_counts[2];
  double complex local_sums[6];
  double complex global_sums[6];
  MVMCPowerLanczosBlockCollectiveResult result;
  MVMCKrylovStatus status;
  const double rank_sum =
      0.5 * (double)world_size * (double)(world_size + 1);
  size_t block;
  size_t entry;
  fill_valid(local_counts, local_sums);
  fill_sentinel(global_counts, global_sums);
  status = mvmc_power_lanczos_block_collective_reduce(
      workspace, MVMC_KRYLOV_STATUS_OK, 2, 3, local_counts, local_sums,
      global_counts, 2, global_sums, 6, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid &&
            result.status == MVMC_KRYLOV_STATUS_OK &&
            result.failure_rank == -1 && result.version == 1 &&
            result.world_size == (size_t)world_size &&
            result.block_count == 2 && result.entries_per_block == 3 &&
            result.local_block_length == 4 &&
            result.global_block_length == (uint64_t)(4 * world_size) &&
            result.global_sample_count == (uint64_t)(8 * world_size),
        "valid reduction metadata mismatch");
  for (block = 0; block < 2; ++block) {
    CHECK(global_counts[block] == (uint64_t)(4 * world_size),
          "global block count mismatch at %zu", block);
    for (entry = 0; entry < 3; ++entry) {
      const double base = (double)(block * 10 + entry + 1);
      const double complex expected =
          rank_sum * base - I * rank_sum * base * 0.5;
      CHECK(global_sums[block * 3 + entry] == expected,
            "rank-ordered sum mismatch at (%zu,%zu)", block, entry);
    }
  }
}

static void test_preflight_failures(
    MVMCPowerLanczosBlockCollectiveWorkspace *workspace) {
  const int target = world_size - 1;
  uint64_t local_counts[2];
  uint64_t global_counts[2];
  uint64_t sentinel_counts[2];
  double complex local_sums[6];
  double complex global_sums[6];
  double complex sentinel_sums[6];
  MVMCPowerLanczosBlockCollectiveResult result;
  MVMCKrylovStatus status;

  fill_valid(local_counts, local_sums);
  fill_sentinel(global_counts, global_sums);
  memcpy(sentinel_counts, global_counts, sizeof(sentinel_counts));
  memcpy(sentinel_sums, global_sums, sizeof(sentinel_sums));
  status = mvmc_power_lanczos_block_collective_reduce(
      workspace,
      world_rank == target ? MVMC_KRYLOV_STATUS_RESOURCE_LIMIT
                           : MVMC_KRYLOV_STATUS_OK,
      2, 3, local_counts, local_sums, global_counts, 2, global_sums, 6,
      &result);
  CHECK(status == MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
            result.failure_rank == target &&
            outputs_unchanged(global_counts, global_sums, sentinel_counts,
                              sentinel_sums),
        "rank-local status did not fail transactionally");

  fill_valid(local_counts, local_sums);
  if (world_rank == target) local_sums[2] = NAN + I;
  status = mvmc_power_lanczos_block_collective_reduce(
      workspace, MVMC_KRYLOV_STATUS_OK, 2, 3, local_counts, local_sums,
      global_counts, 2, global_sums, 6, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_NONFINITE &&
            result.failure_rank == target &&
            outputs_unchanged(global_counts, global_sums, sentinel_counts,
                              sentinel_sums),
        "rank-local nonfinite input did not fail transactionally");

  fill_valid(local_counts, local_sums);
  status = mvmc_power_lanczos_block_collective_reduce(
      workspace, MVMC_KRYLOV_STATUS_OK, 2, 3, local_counts,
      world_rank == target ? NULL : local_sums, global_counts, 2,
      global_sums, 6, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            result.failure_rank == target &&
            outputs_unchanged(global_counts, global_sums, sentinel_counts,
                              sentinel_sums),
        "rank-local null input did not fail transactionally");

  fill_valid(local_counts, local_sums);
  if (world_rank == target) {
    local_counts[0] =
        MVMC_POWER_LANCZOS_BLOCK_MAX_EXACT_SAMPLE_COUNT + UINT64_C(1);
    local_counts[1] = local_counts[0];
  }
  status = mvmc_power_lanczos_block_collective_reduce(
      workspace, MVMC_KRYLOV_STATUS_OK, 2, 3, local_counts, local_sums,
      global_counts, 2, global_sums, 6, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
            result.failure_rank == target &&
            outputs_unchanged(global_counts, global_sums, sentinel_counts,
                              sentinel_sums),
        "inexact local sample count was accepted");
}

static void test_shape_and_count_mismatch(
    MVMCPowerLanczosBlockCollectiveWorkspace *workspace) {
  const int target = world_size - 1;
  uint64_t local_counts[2];
  uint64_t global_counts[2];
  uint64_t sentinel_counts[2];
  double complex local_sums[6];
  double complex global_sums[6];
  double complex sentinel_sums[6];
  MVMCPowerLanczosBlockCollectiveResult result;
  MVMCKrylovStatus status;
  if (world_size <= 1) return;

  fill_valid(local_counts, local_sums);
  fill_sentinel(global_counts, global_sums);
  memcpy(sentinel_counts, global_counts, sizeof(sentinel_counts));
  memcpy(sentinel_sums, global_sums, sizeof(sentinel_sums));
  status = mvmc_power_lanczos_block_collective_reduce(
      workspace, MVMC_KRYLOV_STATUS_OK, 2,
      world_rank == target ? 2U : 3U, local_counts, local_sums,
      global_counts, 2, global_sums, world_rank == target ? 4U : 6U,
      &result);
  CHECK(status == MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE &&
            result.failure_rank == target &&
            outputs_unchanged(global_counts, global_sums, sentinel_counts,
                              sentinel_sums),
        "rank-mismatched shape was not attributed to target");

  fill_valid(local_counts, local_sums);
  if (world_rank == target) {
    local_counts[0] = 5;
    local_counts[1] = 5;
  }
  status = mvmc_power_lanczos_block_collective_reduce(
      workspace, MVMC_KRYLOV_STATUS_OK, 2, 3, local_counts, local_sums,
      global_counts, 2, global_sums, 6, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            result.failure_rank == target &&
            outputs_unchanged(global_counts, global_sums, sentinel_counts,
                              sentinel_sums),
        "rank-mismatched block length was not attributed to target");

  fill_valid(local_counts, local_sums);
  if (world_rank == target) local_counts[1] = 5;
  status = mvmc_power_lanczos_block_collective_reduce(
      workspace, MVMC_KRYLOV_STATUS_OK, 2, 3, local_counts, local_sums,
      global_counts, 2, global_sums, 6, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            result.failure_rank == target &&
            outputs_unchanged(global_counts, global_sums, sentinel_counts,
                              sentinel_sums),
        "rank-local unequal block lengths were accepted");
}

int main(int argc, char **argv) {
  MVMCPowerLanczosBlockCollectiveWorkspace *workspace = NULL;
  MVMCKrylovStatus status;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#else
  (void)argc;
  (void)argv;
#endif
  test_workspace_creation();
#ifdef _mpi_use
  status = mvmc_power_lanczos_block_collective_create(
      MPI_COMM_WORLD, 4, 3, &workspace);
#else
  status = mvmc_power_lanczos_block_collective_create(0, 4, 3, &workspace);
#endif
  CHECK(status == MVMC_KRYLOV_STATUS_OK && workspace != NULL &&
            mvmc_power_lanczos_block_collective_allocated_bytes(workspace) >
                0,
        "workspace creation failed");
  if (workspace != NULL) {
    test_valid_reduce(workspace);
    test_preflight_failures(workspace);
    test_shape_and_count_mismatch(workspace);
    test_valid_reduce(workspace);
  }
  mvmc_power_lanczos_block_collective_destroy(workspace);
#ifdef _mpi_use
  {
    int global_failures = 0;
    MPI_Allreduce(&failures, &global_failures, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    failures = global_failures;
  }
#endif
  if (failures != 0) {
    if (world_rank == 0) {
      fprintf(stderr, "PowerLanczosBlockCollective_Unit: %d failure(s)\n",
              failures);
    }
#ifdef _mpi_use
    MPI_Finalize();
#endif
    return EXIT_FAILURE;
  }
  if (world_rank == 0) {
    printf("PowerLanczosBlockCollective_Unit: PASS (%d rank%s)\n",
           world_size, world_size == 1 ? "" : "s");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}
