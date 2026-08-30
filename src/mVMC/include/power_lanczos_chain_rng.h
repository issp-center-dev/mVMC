#ifndef MVMC_POWER_LANCZOS_CHAIN_RNG_H
#define MVMC_POWER_LANCZOS_CHAIN_RNG_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "bounded_krylov_collective.h"
#include "krylov_positive_sampler.h"
#include "power_lanczos_rng.h"

#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_CHAIN_RNG_VERSION UINT64_C(1)
#define MVMC_POWER_LANCZOS_CHAIN_RNG_INITIAL_DOMAIN \
  UINT64_C(0x5036434841494e31)

typedef struct MVMCPowerLanczosChainRngWorkspace
    MVMCPowerLanczosChainRngWorkspace;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int failure_world_rank;
  uint64_t version;
  MVMCPowerLanczosRngStage stage;
  size_t chain_rank;
  size_t chain_size;
  uint64_t proposal_ordinal;
  uint64_t domain_draws_before;
  uint64_t domain_draws_after;
  uint64_t proposal_seed;
  uint64_t proposal_stream;
} MVMCPowerLanczosChainRngResult;

/*
 * Creation is collective on chain_communicator.  The communicator must be
 * the existing contiguous comm_child1 group formed by
 * world_rank / split_size.  The final group may be shorter.
 */
MVMCKrylovStatus mvmc_power_lanczos_chain_rng_create(
    MVMCKrylovBoundedCommunicator chain_communicator,
    size_t mpi_world_rank, size_t mpi_world_size, size_t split_size,
    MVMCPowerLanczosChainRngWorkspace **workspace);

void mvmc_power_lanczos_chain_rng_destroy(
    MVMCPowerLanczosChainRngWorkspace *workspace);

size_t mvmc_power_lanczos_chain_rng_allocated_bytes(
    const MVMCPowerLanczosChainRngWorkspace *workspace);

/*
 * Consume exactly one word from every rank-domain RNG, gather those words
 * in chain-rank/world-rank order, and return an identical proposal-local RNG
 * on every chain rank.  Both RNG outputs remain untouched on failure.
 */
MVMCKrylovStatus mvmc_power_lanczos_chain_rng_derive_proposal(
    MVMCPowerLanczosChainRngWorkspace *workspace,
    const MVMCPowerLanczosRngDomain *domain,
    MVMCKrylovPositiveSamplerRng *rank_domain_rng,
    uint64_t proposal_ordinal,
    MVMCKrylovPositiveSamplerRng *proposal_rng,
    MVMCPowerLanczosChainRngResult *result);

#endif

#endif
