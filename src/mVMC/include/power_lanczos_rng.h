#ifndef MVMC_POWER_LANCZOS_RNG_H
#define MVMC_POWER_LANCZOS_RNG_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "krylov_fock_reference.h"

#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_RNG_DERIVATION_VERSION UINT64_C(2)
#define MVMC_POWER_LANCZOS_RNG_INITIAL_DOMAIN UINT64_C(0x50364c5a524e4732)

typedef enum {
  MVMC_POWER_LANCZOS_RNG_STAGE_BASE = 1,
  MVMC_POWER_LANCZOS_RNG_STAGE_COEFFICIENT = 2,
  MVMC_POWER_LANCZOS_RNG_STAGE_FINAL = 3
} MVMCPowerLanczosRngStage;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  uint64_t version;
  uint64_t resolved_base_seed;
  uint64_t run_index;
  MVMCPowerLanczosRngStage stage;
  uint64_t stage_tag;
  size_t mpi_world_rank;
  size_t mpi_world_size;
  size_t split_rank;
  size_t split_size;
  uint64_t derived_seed;
} MVMCPowerLanczosRngDomain;

MVMCKrylovStatus mvmc_power_lanczos_rng_derive(
    uint64_t resolved_base_seed, uint64_t run_index,
    MVMCPowerLanczosRngStage stage, size_t mpi_world_rank,
    size_t mpi_world_size, size_t split_size,
    MVMCPowerLanczosRngDomain *domain);

MVMCKrylovStatus mvmc_power_lanczos_rng_domain_matches(
    const MVMCPowerLanczosRngDomain *expected,
    uint64_t resolved_base_seed, uint64_t run_index,
    MVMCPowerLanczosRngStage stage, size_t mpi_world_rank,
    size_t mpi_world_size, size_t split_size, int *matches);

const char *mvmc_power_lanczos_rng_stage_string(
    MVMCPowerLanczosRngStage stage);

#endif

#endif
