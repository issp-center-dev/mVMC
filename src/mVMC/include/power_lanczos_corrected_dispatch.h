#ifndef MVMC_POWER_LANCZOS_CORRECTED_DISPATCH_H
#define MVMC_POWER_LANCZOS_CORRECTED_DISPATCH_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "power_lanczos_classic_bridge.h"
#include "power_lanczos_output_transaction.h"
#include "power_lanczos_stabilization_output.h"

#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_CORRECTED_DISPATCH_VERSION UINT64_C(1)

typedef enum {
  MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_PROVIDED_TERMINAL = 0,
  MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_UNIFORM_SECTOR = 1
} MVMCPowerLanczosCorrectedBootstrapMode;

typedef struct {
  size_t coefficient_warm_up;
  size_t coefficient_sample_count;
  size_t coefficient_interval;
  size_t final_warm_up;
  size_t final_sample_count;
  size_t final_interval;
} MVMCPowerLanczosCorrectedChainControls;

typedef struct {
  const MVMCPowerLanczosClassicView *classic_view;
  MVMCClassicPfaffianCommunicator world_communicator;
  MVMCClassicPfaffianCommunicator chain_communicator;
  int power_step;
  uint64_t resolved_base_seed;
  uint64_t run_index;
  size_t mpi_world_rank;
  size_t mpi_world_size;
  size_t split_size;
  uint64_t base_generation;
  MVMCPowerLanczosCorrectedBootstrapMode bootstrap_mode;
  const uint64_t *root_terminal_words;
  size_t root_terminal_word_count;
  uint64_t terminal_proposal_counter;
  MVMCPowerLanczosCorrectedChainControls controls;
  const MVMCPowerLanczosStabilizationOutputIdentity *root_identity;
  const char *root_output_directory;
  const char *root_output_basename;
} MVMCPowerLanczosCorrectedDispatchInput;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  uint64_t version;
  MVMCPowerLanczosCorrectedBootstrapMode bootstrap_mode;
  uint64_t bootstrap_proposal_counter;
  MVMCPowerLanczosStabilizationDecision decision;
  uint64_t failed_gates;
  double energy;
  double energy_standard_error;
  double variance;
  double variance_standard_error;
  size_t output_size;
  MVMCPowerLanczosOutputTransactionStatus output_status;
  char output_sha256[MVMC_POWER_LANCZOS_SHA256_HEX_CAPACITY];
} MVMCPowerLanczosCorrectedDispatchResult;

/*
 * Collective corrected energy/variance dispatch.  The root-only terminal,
 * identity, and output fields are ignored on non-root ranks.  No observable
 * plan is accepted and no legacy estimator is evaluated.
 */
MVMCKrylovStatus mvmc_power_lanczos_corrected_dispatch_run(
    const MVMCPowerLanczosCorrectedDispatchInput *input,
    MVMCPowerLanczosCorrectedDispatchResult *result);

#endif /* bounded engine */

#endif /* MVMC_POWER_LANCZOS_CORRECTED_DISPATCH_H */
