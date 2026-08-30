#include "power_lanczos_chain_rng.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "power_lanczos_chain_rng.c requires the bounded Krylov core"
#endif

#include <limits.h>
#include <stdlib.h>
#include <string.h>

struct MVMCPowerLanczosChainRngWorkspace {
  int chain_rank;
  int chain_size;
  size_t mpi_world_rank;
  size_t mpi_world_size;
  size_t split_size;
  size_t allocated_bytes;
  uint64_t *world_ranks;
  uint64_t *rank_words;
#ifdef _mpi_use
  MPI_Comm communicator;
#endif
};

static uint64_t splitmix64(uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

static int checked_multiply(size_t left, size_t right, size_t *result) {
  if (result == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *result = left * right;
  return 1;
}

static int valid_status(MVMCKrylovStatus status) {
  return status >= MVMC_KRYLOV_STATUS_OK &&
         status <= MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

static int valid_sampler_rng(const MVMCKrylovPositiveSamplerRng *rng) {
  return rng != NULL && rng->valid &&
         rng->state_version == MVMC_KRYLOV_POSITIVE_SAMPLER_STATE_VERSION &&
         rng->algorithm == MVMC_KRYLOV_POSITIVE_SAMPLER_RNG_SPLITMIX64;
}

static void reset_result(MVMCPowerLanczosChainRngResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  result->failure_world_rank = -1;
  result->version = MVMC_POWER_LANCZOS_CHAIN_RNG_VERSION;
}

static MVMCKrylovStatus synchronize_status(
    MVMCPowerLanczosChainRngWorkspace *workspace,
    MVMCKrylovStatus local_status, int result_available,
    MVMCPowerLanczosChainRngResult *result) {
  MVMCKrylovStatus effective = valid_status(local_status)
                                      ? local_status
                                      : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  int failure_world_rank =
      effective == MVMC_KRYLOV_STATUS_OK ? -1 : (int)workspace->mpi_world_rank;
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
    local.index = (int)workspace->mpi_world_rank;
    if (MPI_Allreduce(&local, &global, 1, MPI_2INT, MPI_MAXLOC,
                      workspace->communicator) != MPI_SUCCESS) {
      if (result_available) {
        result->status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
      }
      return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
    effective = (MVMCKrylovStatus)global.value;
    failure_world_rank =
        effective == MVMC_KRYLOV_STATUS_OK ? -1 : global.index;
  }
#endif
  if (result_available) {
    result->valid = effective == MVMC_KRYLOV_STATUS_OK;
    result->status = effective;
    result->failure_world_rank = failure_world_rank;
  }
  return effective;
}

MVMCKrylovStatus mvmc_power_lanczos_chain_rng_create(
    MVMCKrylovBoundedCommunicator chain_communicator,
    size_t mpi_world_rank, size_t mpi_world_size, size_t split_size,
    MVMCPowerLanczosChainRngWorkspace **workspace) {
  MVMCPowerLanczosChainRngWorkspace *created = NULL;
  int chain_rank = 0;
  int chain_size = 1;
  int local_valid = workspace != NULL && mpi_world_size > 0 &&
                    split_size > 0 && mpi_world_rank < mpi_world_size &&
                    mpi_world_rank <= (size_t)INT_MAX;
  size_t group_base = 0;
  size_t expected_chain_size = 0;
  size_t buffer_bytes = 0;
  size_t total_bytes = 0;
#ifdef _mpi_use
  MPI_Comm duplicated = MPI_COMM_NULL;
  int initialized = 0;
  int all_valid = 0;
  int allocation_ok;
  int all_allocation_ok;
#endif
  if (workspace != NULL) *workspace = NULL;
#ifdef _mpi_use
  if (chain_communicator == MPI_COMM_NULL ||
      MPI_Initialized(&initialized) != MPI_SUCCESS || !initialized ||
      MPI_Comm_dup(chain_communicator, &duplicated) != MPI_SUCCESS ||
      MPI_Comm_set_errhandler(duplicated, MPI_ERRORS_RETURN) != MPI_SUCCESS ||
      MPI_Comm_rank(duplicated, &chain_rank) != MPI_SUCCESS ||
      MPI_Comm_size(duplicated, &chain_size) != MPI_SUCCESS ||
      chain_size <= 0) {
    if (duplicated != MPI_COMM_NULL) MPI_Comm_free(&duplicated);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#if SIZE_MAX > UINT64_MAX
  if (mpi_world_rank > UINT64_MAX || mpi_world_size > UINT64_MAX ||
      split_size > UINT64_MAX) {
    local_valid = 0;
  }
#endif
  if (MPI_Allreduce(&local_valid, &all_valid, 1, MPI_INT, MPI_MIN,
                    duplicated) != MPI_SUCCESS) {
    MPI_Comm_free(&duplicated);
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (!all_valid) {
    MPI_Comm_free(&duplicated);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#else
  (void)chain_communicator;
  if (!local_valid || mpi_world_rank != 0 || mpi_world_size != 1 ||
      split_size != 1) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#endif
  group_base = (mpi_world_rank / split_size) * split_size;
  expected_chain_size = mpi_world_size - group_base;
  if (expected_chain_size > split_size) expected_chain_size = split_size;
  local_valid = local_valid &&
                (size_t)chain_rank == mpi_world_rank - group_base &&
                (size_t)chain_size == expected_chain_size &&
                checked_multiply((size_t)chain_size, sizeof(uint64_t),
                                 &buffer_bytes) &&
                buffer_bytes <= (SIZE_MAX - sizeof(*created)) / 2;
#ifdef _mpi_use
  if (MPI_Allreduce(&local_valid, &all_valid, 1, MPI_INT, MPI_MIN,
                    duplicated) != MPI_SUCCESS) {
    MPI_Comm_free(&duplicated);
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#endif
  if (!local_valid
#ifdef _mpi_use
      || !all_valid
#endif
  ) {
#ifdef _mpi_use
    MPI_Comm_free(&duplicated);
#endif
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  total_bytes = sizeof(*created) + 2 * buffer_bytes;
  created = (MVMCPowerLanczosChainRngWorkspace *)calloc(1, total_bytes);
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
  created->chain_rank = chain_rank;
  created->chain_size = chain_size;
  created->mpi_world_rank = mpi_world_rank;
  created->mpi_world_size = mpi_world_size;
  created->split_size = split_size;
  created->allocated_bytes = total_bytes;
  created->world_ranks =
      (uint64_t *)((unsigned char *)created + sizeof(*created));
  created->rank_words =
      (uint64_t *)((unsigned char *)created + sizeof(*created) + buffer_bytes);
#ifdef _mpi_use
  created->communicator = duplicated;
  {
    const uint64_t local_world_rank = (uint64_t)mpi_world_rank;
    if (MPI_Allgather(&local_world_rank, 1, MPI_UINT64_T,
                      created->world_ranks, 1, MPI_UINT64_T,
                      created->communicator) != MPI_SUCCESS) {
      mvmc_power_lanczos_chain_rng_destroy(created);
      return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
  }
#else
  created->world_ranks[0] = 0;
#endif
  {
    int rank;
    for (rank = 0; rank < chain_size; ++rank) {
      if (created->world_ranks[rank] != (uint64_t)(group_base + (size_t)rank)) {
        mvmc_power_lanczos_chain_rng_destroy(created);
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      }
    }
  }
  *workspace = created;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_power_lanczos_chain_rng_destroy(
    MVMCPowerLanczosChainRngWorkspace *workspace) {
  if (workspace == NULL) return;
#ifdef _mpi_use
  if (workspace->communicator != MPI_COMM_NULL) {
    MPI_Comm_free(&workspace->communicator);
  }
#endif
  memset(workspace, 0, workspace->allocated_bytes);
  free(workspace);
}

size_t mvmc_power_lanczos_chain_rng_allocated_bytes(
    const MVMCPowerLanczosChainRngWorkspace *workspace) {
  return workspace == NULL ? 0 : workspace->allocated_bytes;
}

MVMCKrylovStatus mvmc_power_lanczos_chain_rng_derive_proposal(
    MVMCPowerLanczosChainRngWorkspace *workspace,
    const MVMCPowerLanczosRngDomain *domain,
    MVMCKrylovPositiveSamplerRng *rank_domain_rng,
    uint64_t proposal_ordinal,
    MVMCKrylovPositiveSamplerRng *proposal_rng,
    MVMCPowerLanczosChainRngResult *result) {
  MVMCKrylovPositiveSamplerRng candidate_domain;
  MVMCKrylovPositiveSamplerRng candidate_proposal;
  MVMCKrylovStatus effective = MVMC_KRYLOV_STATUS_OK;
  const int result_available = result != NULL;
  int domain_matches = 0;
  uint64_t local_word = 0;
  uint64_t hash = MVMC_POWER_LANCZOS_CHAIN_RNG_INITIAL_DOMAIN;
  uint64_t proposal_seed;
  uint64_t proposal_stream;
  int rank;
  if (result_available) reset_result(result);
  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(&candidate_domain, 0, sizeof(candidate_domain));
  memset(&candidate_proposal, 0, sizeof(candidate_proposal));
  if (domain == NULL || !domain->valid ||
      domain->status != MVMC_KRYLOV_STATUS_OK ||
      !valid_sampler_rng(rank_domain_rng) || proposal_rng == NULL ||
      domain->stage == MVMC_POWER_LANCZOS_RNG_STAGE_BASE ||
      rank_domain_rng->stream != domain->stage_tag ||
      rank_domain_rng->draws != proposal_ordinal) {
    effective = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (effective == MVMC_KRYLOV_STATUS_OK) {
    effective = mvmc_power_lanczos_rng_domain_matches(
        domain, domain->resolved_base_seed, domain->run_index,
        domain->stage, workspace->mpi_world_rank,
        workspace->mpi_world_size, workspace->split_size,
        &domain_matches);
    if (effective == MVMC_KRYLOV_STATUS_OK && !domain_matches) {
      effective = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
  }
  effective = synchronize_status(workspace, effective, result_available,
                                 result);
  if (effective != MVMC_KRYLOV_STATUS_OK) return effective;
  candidate_domain = *rank_domain_rng;
  effective = mvmc_krylov_positive_sampler_rng_draw_uint64(
      &candidate_domain, &local_word);
  effective = synchronize_status(workspace, effective, result_available,
                                 result);
  if (effective != MVMC_KRYLOV_STATUS_OK) return effective;
#ifdef _mpi_use
  if (MPI_Allgather(&local_word, 1, MPI_UINT64_T, workspace->rank_words, 1,
                    MPI_UINT64_T, workspace->communicator) != MPI_SUCCESS) {
    reset_result(result);
    result->status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    return result->status;
  }
#else
  workspace->rank_words[0] = local_word;
#endif
  hash = splitmix64(hash ^ splitmix64(domain->stage_tag ^
                                      UINT64_C(0x5354414745544147)));
  hash = splitmix64(hash ^ splitmix64(proposal_ordinal ^
                                      UINT64_C(0x50524f504f524449)));
  for (rank = 0; rank < workspace->chain_size; ++rank) {
    hash = splitmix64(hash ^ splitmix64(
        workspace->world_ranks[rank] ^ UINT64_C(0x574f524c4452414e)));
    hash = splitmix64(hash ^ splitmix64(
        workspace->rank_words[rank] ^ UINT64_C(0x52414e4b574f5244)));
  }
  proposal_seed = splitmix64(hash);
  proposal_stream = splitmix64(
      domain->stage_tag ^ proposal_ordinal ^ UINT64_C(0x50524f505354524d));
  effective = mvmc_krylov_positive_sampler_rng_seed(
      proposal_seed, proposal_stream, &candidate_proposal);
  effective = synchronize_status(workspace, effective, result_available,
                                 result);
  if (effective != MVMC_KRYLOV_STATUS_OK) return effective;
  *rank_domain_rng = candidate_domain;
  *proposal_rng = candidate_proposal;
  result->valid = 1;
  result->status = MVMC_KRYLOV_STATUS_OK;
  result->failure_world_rank = -1;
  result->stage = domain->stage;
  result->chain_rank = (size_t)workspace->chain_rank;
  result->chain_size = (size_t)workspace->chain_size;
  result->proposal_ordinal = proposal_ordinal;
  result->domain_draws_before = proposal_ordinal;
  result->domain_draws_after = candidate_domain.draws;
  result->proposal_seed = proposal_seed;
  result->proposal_stream = proposal_stream;
  return MVMC_KRYLOV_STATUS_OK;
}
