#ifndef MVMC_POWER_LANCZOS_STABILIZATION_OUTPUT_H
#define MVMC_POWER_LANCZOS_STABILIZATION_OUTPUT_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "power_lanczos_stabilization_statistics.h"

#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_STABILIZATION_OUTPUT_SCHEMA_VERSION UINT64_C(1)

typedef struct {
  const char *run_id;
  const char *source_commit;
  const char *input_sha256;
  const char *binary_sha256;
  const char *environment_id;
  const char *seed_id;
} MVMCPowerLanczosStabilizationOutputIdentity;

/*
 * Render one deterministic compact JSON record.  Sample-level primitive
 * values, bounds, support flags, and tails are intentionally never emitted.
 * The API is world-rank-0 only because gathered trace evidence is root-owned.
 */
MVMCKrylovStatus mvmc_power_lanczos_stabilization_output_render(
    const MVMCPowerLanczosProductionSession *session,
    const MVMCPowerLanczosStabilizationResult *statistics,
    const MVMCPowerLanczosStabilizationOutputIdentity *identity,
    char *output, size_t output_capacity, size_t *output_size);

#endif /* MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE */

#endif /* MVMC_POWER_LANCZOS_STABILIZATION_OUTPUT_H */
