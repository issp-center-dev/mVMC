#ifndef MVMC_POWER_LANCZOS_STABILIZED_H
#define MVMC_POWER_LANCZOS_STABILIZED_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "power_lanczos_classic_bridge.h"

#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_STABILIZED_VERSION UINT64_C(2)
#define MVMC_POWER_LANCZOS_STABILIZED_BLOCKS ((size_t)16)

static inline size_t mvmc_power_lanczos_stabilized_block_count(
    size_t sample_count) {
  size_t block_count = sample_count >= 32
                           ? MVMC_POWER_LANCZOS_STABILIZED_BLOCKS
                           : sample_count / 2;
  while (block_count >= 4 && sample_count % block_count != 0) {
    --block_count;
  }
  return block_count >= 4 ? block_count : 0;
}

typedef struct {
  const MVMCPowerLanczosClassicView *classic_view;
  MVMCClassicPfaffianCommunicator world_communicator;
  MVMCClassicPfaffianCommunicator chain_communicator;
  int power_step;
  uint64_t seed;
  size_t warm_up;
  size_t sample_count;
  size_t interval;
} MVMCPowerLanczosStabilizedInput;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  uint64_t version;
  int power_step;
  size_t block_count;
  uint64_t sampling_chains;
  uint64_t coefficient_samples;
  uint64_t final_samples;
  int retained_rank;
  double condition_estimate;
  double gevp_residual;
  double energy;
  double energy_standard_error;
  double variance;
  double variance_standard_error;
  double final_energy_imaginary;
  double variance_imaginary;
  double energy_tau_int;
  double effective_sample_count;
} MVMCPowerLanczosStabilizedResult;

MVMCKrylovStatus mvmc_power_lanczos_stabilized_run(
    const MVMCPowerLanczosStabilizedInput *input,
    MVMCPowerLanczosStabilizedResult *result);

#endif /* MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE */

#endif /* MVMC_POWER_LANCZOS_STABILIZED_H */
