#include "power_lanczos_rng.h"

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include <string.h>

enum {
  FIELD_COUNT = 5
};

static const uint64_t FieldTags[FIELD_COUNT] = {
    UINT64_C(0x4241534553454544), UINT64_C(0x52554e494e444558),
    UINT64_C(0x5354414745544147), UINT64_C(0x574f524c4452414e),
    UINT64_C(0x53504c495452414e)};

static uint64_t splitmix64(uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

static int stage_tag(MVMCPowerLanczosRngStage stage, uint64_t *tag) {
  if (tag == NULL) return 0;
  switch (stage) {
    case MVMC_POWER_LANCZOS_RNG_STAGE_BASE:
      *tag = UINT64_C(0x5036424153453032);
      return 1;
    case MVMC_POWER_LANCZOS_RNG_STAGE_COEFFICIENT:
      *tag = UINT64_C(0x5036434f45464632);
      return 1;
    case MVMC_POWER_LANCZOS_RNG_STAGE_FINAL:
      *tag = UINT64_C(0x503646494e414c32);
      return 1;
    default:
      return 0;
  }
}

MVMCKrylovStatus mvmc_power_lanczos_rng_derive(
    uint64_t resolved_base_seed, uint64_t run_index,
    MVMCPowerLanczosRngStage stage, size_t mpi_world_rank,
    size_t mpi_world_size, size_t split_size,
    MVMCPowerLanczosRngDomain *domain) {
  MVMCPowerLanczosRngDomain candidate;
  uint64_t fields[FIELD_COUNT];
  uint64_t tag = 0;
  uint64_t hash = MVMC_POWER_LANCZOS_RNG_INITIAL_DOMAIN;
  size_t split_rank;
  size_t field;
  if (domain == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(domain, 0, sizeof(*domain));
  domain->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (!stage_tag(stage, &tag) || mpi_world_size == 0 || split_size == 0 ||
      mpi_world_rank >= mpi_world_size ||
      mpi_world_rank > (size_t)UINT64_MAX ||
      mpi_world_size > (size_t)UINT64_MAX ||
      split_size > (size_t)UINT64_MAX) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  split_rank = mpi_world_rank % split_size;
  fields[0] = resolved_base_seed;
  fields[1] = run_index;
  fields[2] = tag;
  fields[3] = (uint64_t)mpi_world_rank;
  fields[4] = (uint64_t)split_rank;
  for (field = 0; field < FIELD_COUNT; ++field) {
    hash = splitmix64(hash ^ splitmix64(fields[field] ^ FieldTags[field]));
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.version = MVMC_POWER_LANCZOS_RNG_DERIVATION_VERSION;
  candidate.resolved_base_seed = resolved_base_seed;
  candidate.run_index = run_index;
  candidate.stage = stage;
  candidate.stage_tag = tag;
  candidate.mpi_world_rank = mpi_world_rank;
  candidate.mpi_world_size = mpi_world_size;
  candidate.split_rank = split_rank;
  candidate.split_size = split_size;
  candidate.derived_seed = splitmix64(hash);
  *domain = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_rng_domain_matches(
    const MVMCPowerLanczosRngDomain *expected,
    uint64_t resolved_base_seed, uint64_t run_index,
    MVMCPowerLanczosRngStage stage, size_t mpi_world_rank,
    size_t mpi_world_size, size_t split_size, int *matches) {
  MVMCPowerLanczosRngDomain actual;
  MVMCKrylovStatus status;
  if (matches == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *matches = 0;
  if (expected == NULL || !expected->valid ||
      expected->status != MVMC_KRYLOV_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_power_lanczos_rng_derive(
      resolved_base_seed, run_index, stage, mpi_world_rank,
      mpi_world_size, split_size, &actual);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  *matches = expected->valid == actual.valid &&
             expected->status == actual.status &&
             expected->version == actual.version &&
             expected->resolved_base_seed == actual.resolved_base_seed &&
             expected->run_index == actual.run_index &&
             expected->stage == actual.stage &&
             expected->stage_tag == actual.stage_tag &&
             expected->mpi_world_rank == actual.mpi_world_rank &&
             expected->mpi_world_size == actual.mpi_world_size &&
             expected->split_rank == actual.split_rank &&
             expected->split_size == actual.split_size &&
             expected->derived_seed == actual.derived_seed;
  return MVMC_KRYLOV_STATUS_OK;
}

const char *mvmc_power_lanczos_rng_stage_string(
    MVMCPowerLanczosRngStage stage) {
  switch (stage) {
    case MVMC_POWER_LANCZOS_RNG_STAGE_BASE:
      return "base";
    case MVMC_POWER_LANCZOS_RNG_STAGE_COEFFICIENT:
      return "coefficient";
    case MVMC_POWER_LANCZOS_RNG_STAGE_FINAL:
      return "final";
    default:
      return "invalid";
  }
}

#endif
