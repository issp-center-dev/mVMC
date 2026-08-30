#include "power_lanczos_chain_rng.h"

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
      fprintf(stderr, "PowerLanczosChainRng_Unit FAIL rank %d: ",             \
              world_rank);                                                      \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, " (line %d)\n", __LINE__);                              \
      ++failures;                                                               \
    }                                                                           \
  } while (0)

#ifdef _mpi_use
static void check_equal_u64(MPI_Comm communicator, uint64_t value,
                            const char *label) {
  uint64_t minimum = 0;
  uint64_t maximum = 0;
  CHECK(MPI_Allreduce(&value, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                      communicator) == MPI_SUCCESS &&
            MPI_Allreduce(&value, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                          communicator) == MPI_SUCCESS &&
            minimum == maximum,
        "%s differed within chain", label);
}
#endif

static void run_chain(MVMCKrylovBoundedCommunicator communicator,
                      size_t split_size) {
  MVMCPowerLanczosChainRngWorkspace *workspace = NULL;
  MVMCPowerLanczosRngDomain domain;
  MVMCPowerLanczosRngDomain bad_domain;
  MVMCKrylovPositiveSamplerRng rank_rng;
  MVMCKrylovPositiveSamplerRng proposal_rng;
  MVMCKrylovPositiveSamplerRng rank_snapshot;
  MVMCKrylovPositiveSamplerRng proposal_snapshot;
  MVMCPowerLanczosChainRngResult result;
  MVMCKrylovStatus status;
  uint64_t first_value = 0;
  int chain_rank = 0;
  int chain_size = 1;
  const size_t target_world_rank =
      ((size_t)world_rank / split_size) * split_size + split_size - 1 <
              (size_t)world_size
          ? ((size_t)world_rank / split_size) * split_size + split_size - 1
          : (size_t)world_size - 1;
#ifdef _mpi_use
  MPI_Comm_rank(communicator, &chain_rank);
  MPI_Comm_size(communicator, &chain_size);
#endif
  status = mvmc_power_lanczos_chain_rng_create(
      communicator, (size_t)world_rank, (size_t)world_size, split_size,
      &workspace);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && workspace != NULL &&
            mvmc_power_lanczos_chain_rng_allocated_bytes(workspace) > 0,
        "chain workspace creation failed");
  if (workspace == NULL) return;
  CHECK(mvmc_power_lanczos_rng_derive(
            UINT64_C(12345), 0, MVMC_POWER_LANCZOS_RNG_STAGE_COEFFICIENT,
            (size_t)world_rank, (size_t)world_size, split_size, &domain) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_rng_seed(
                domain.derived_seed, domain.stage_tag, &rank_rng) ==
                MVMC_KRYLOV_STATUS_OK,
        "rank domain initialization failed");
  memset(&proposal_rng, 0xa5, sizeof(proposal_rng));
  status = mvmc_power_lanczos_chain_rng_derive_proposal(
      workspace, &domain, &rank_rng, 0, &proposal_rng, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid &&
            result.status == MVMC_KRYLOV_STATUS_OK &&
            result.failure_world_rank == -1 && result.version == 1 &&
            result.stage == MVMC_POWER_LANCZOS_RNG_STAGE_COEFFICIENT &&
            result.chain_rank == (size_t)chain_rank &&
            result.chain_size == (size_t)chain_size &&
            result.proposal_ordinal == 0 &&
            result.domain_draws_before == 0 &&
            result.domain_draws_after == 1 && rank_rng.draws == 1 &&
            proposal_rng.draws == 0 &&
            proposal_rng.state_version ==
                MVMC_KRYLOV_POSITIVE_SAMPLER_STATE_VERSION,
        "first proposal derivation metadata mismatch");
  if (world_size == 1 && split_size == 1) {
    CHECK(domain.derived_seed == UINT64_C(0x8092f5366aaa14a3) &&
              result.proposal_seed == UINT64_C(0x96e5c74e335d1f06) &&
              result.proposal_stream == UINT64_C(0x71b286e9b8c31f5f) &&
              proposal_rng.state == UINT64_C(0x688350093aa6dfdc),
          "single-chain frozen vector mismatch domain=%016llx seed=%016llx "
          "stream=%016llx state=%016llx",
          (unsigned long long)domain.derived_seed,
          (unsigned long long)result.proposal_seed,
          (unsigned long long)result.proposal_stream,
          (unsigned long long)proposal_rng.state);
  }
#ifdef _mpi_use
  check_equal_u64(communicator, result.proposal_seed, "proposal seed");
  check_equal_u64(communicator, result.proposal_stream, "proposal stream");
  check_equal_u64(communicator, proposal_rng.state, "proposal RNG state");
#endif
  CHECK(mvmc_krylov_positive_sampler_rng_draw_uint64(
            &proposal_rng, &first_value) == MVMC_KRYLOV_STATUS_OK,
        "proposal RNG first draw failed");
#ifdef _mpi_use
  check_equal_u64(communicator, first_value, "proposal first value");
#endif

  status = mvmc_power_lanczos_chain_rng_derive_proposal(
      workspace, &domain, &rank_rng, 1, &proposal_rng, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && rank_rng.draws == 2 &&
            result.proposal_ordinal == 1 &&
            result.domain_draws_before == 1 &&
            result.domain_draws_after == 2,
        "second proposal derivation failed");

  rank_snapshot = rank_rng;
  proposal_snapshot = proposal_rng;
  status = mvmc_power_lanczos_chain_rng_derive_proposal(
      workspace, &domain, &rank_rng,
      (size_t)world_rank == target_world_rank ? 3 : 2,
      &proposal_rng, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            result.failure_world_rank == (int)target_world_rank &&
            memcmp(&rank_rng, &rank_snapshot, sizeof(rank_rng)) == 0 &&
            memcmp(&proposal_rng, &proposal_snapshot,
                   sizeof(proposal_rng)) == 0,
        "rank-local ordinal mismatch was not transactional");

  bad_domain = domain;
  if ((size_t)world_rank == target_world_rank) bad_domain.derived_seed ^= 1;
  status = mvmc_power_lanczos_chain_rng_derive_proposal(
      workspace, &bad_domain, &rank_rng, 2, &proposal_rng, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE &&
            result.failure_world_rank == (int)target_world_rank &&
            memcmp(&rank_rng, &rank_snapshot, sizeof(rank_rng)) == 0 &&
            memcmp(&proposal_rng, &proposal_snapshot,
                   sizeof(proposal_rng)) == 0,
        "rank-local domain mutation was not transactional");

  status = mvmc_power_lanczos_chain_rng_derive_proposal(
      workspace, &domain, &rank_rng, 2,
      (size_t)world_rank == target_world_rank ? NULL : &proposal_rng,
      &result);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            result.failure_world_rank == (int)target_world_rank &&
            memcmp(&rank_rng, &rank_snapshot, sizeof(rank_rng)) == 0,
        "rank-local null proposal output was not transactional");

  status = mvmc_power_lanczos_chain_rng_derive_proposal(
      workspace, &domain, &rank_rng, 2, &proposal_rng, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && rank_rng.draws == 3,
        "workspace did not recover after rejected calls");
  mvmc_power_lanczos_chain_rng_destroy(workspace);
}

static void test_invalid_topology(void) {
#ifdef _mpi_use
  MVMCPowerLanczosChainRngWorkspace *workspace = NULL;
  MVMCKrylovStatus status = mvmc_power_lanczos_chain_rng_create(
      MPI_COMM_WORLD, (size_t)world_rank, (size_t)world_size, 1,
      &workspace);
  if (world_size > 1) {
    CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
          "noncontiguous chain topology was accepted");
  } else {
    CHECK(status == MVMC_KRYLOV_STATUS_OK && workspace != NULL,
          "single-rank topology fixture failed");
    mvmc_power_lanczos_chain_rng_destroy(workspace);
  }
#endif
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Comm grouped = MPI_COMM_NULL;
  size_t grouped_split;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#else
  (void)argc;
  (void)argv;
#endif
  test_invalid_topology();
#ifdef _mpi_use
  run_chain(MPI_COMM_WORLD, (size_t)world_size);
  grouped_split = world_size >= 2 ? 2U : 1U;
  MPI_Comm_split(MPI_COMM_WORLD,
                 world_rank / (int)grouped_split, world_rank, &grouped);
  run_chain(grouped, grouped_split);
  MPI_Comm_free(&grouped);
#else
  run_chain(0, 1);
#endif
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
      fprintf(stderr, "PowerLanczosChainRng_Unit: %d failure(s)\n",
              failures);
    }
#ifdef _mpi_use
    MPI_Finalize();
#endif
    return EXIT_FAILURE;
  }
  if (world_rank == 0) {
    printf("PowerLanczosChainRng_Unit: PASS (%d rank%s)\n", world_size,
           world_size == 1 ? "" : "s");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}
