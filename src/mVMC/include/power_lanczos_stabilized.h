#ifndef MVMC_POWER_LANCZOS_STABILIZED_H
#define MVMC_POWER_LANCZOS_STABILIZED_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "power_lanczos_classic_bridge.h"
#include "power_lanczos_gc_bridge.h"

#include <complex.h>
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

/* Final-state particle numbers; populated only by the GC entry. The
 * blocks use the same global (chain-major) partition as energy. */
typedef struct {
  int valid;
  double mean;
  double standard_error;
  double mean_square;
  double mean_square_standard_error;
  double variance;
  uint64_t block_sample_count[MVMC_POWER_LANCZOS_STABILIZED_BLOCKS];
  double block_sum[MVMC_POWER_LANCZOS_STABILIZED_BLOCKS];
  double block_square_sum[MVMC_POWER_LANCZOS_STABILIZED_BLOCKS];
} MVMCPowerLanczosGCNumberResult;

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
  /* E_coeff = c^dagger K c / c^dagger S c of the coefficient chain. */
  int coefficient_energy_available;
  double complex coefficient_energy;
  MVMCPowerLanczosGCNumberResult gc_number;
} MVMCPowerLanczosStabilizedResult;

MVMCKrylovStatus mvmc_power_lanczos_stabilized_run(
    const MVMCPowerLanczosStabilizedInput *input,
    MVMCPowerLanczosStabilizedResult *result);

/*
 * Grand-canonical first-step entry (power step 1).  The coefficient chain
 * starts from a uniform configuration with initial_electron_count
 * electrons; the final chain tries initial-N and whole-domain uniform
 * configurations alternately (at most 1024 attempts) and redraws only
 * exact zero, numeric zero or unavailable local energy.  All proposals are
 * the symmetric GC kernels (15/16 local, 1/16 uniform) drawn by the chain
 * root and broadcast with the acceptance uniform.
 */
typedef struct {
  const MVMCPowerLanczosGCView *gc_view;
  MVMCClassicPfaffianCommunicator world_communicator;
  MVMCClassicPfaffianCommunicator chain_communicator;
  uint64_t seed;
  size_t warm_up;
  size_t sample_count;
  size_t interval;
  int initial_electron_count;
} MVMCPowerLanczosGCInput;

typedef enum {
  MVMC_POWER_LANCZOS_GC_STAGE_NONE = 0,
  MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_INITIAL,
  MVMC_POWER_LANCZOS_GC_STAGE_PILOT_WARM_UP,
  MVMC_POWER_LANCZOS_GC_STAGE_PILOT_SAMPLE,
  MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_WARM_UP,
  MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_SAMPLE,
  MVMC_POWER_LANCZOS_GC_STAGE_FINAL_BOOTSTRAP,
  MVMC_POWER_LANCZOS_GC_STAGE_FINAL_WARM_UP,
  MVMC_POWER_LANCZOS_GC_STAGE_FINAL_SAMPLE,
  MVMC_POWER_LANCZOS_GC_STAGE_COUNT
} MVMCPowerLanczosGCStage;

#define MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS 16

/*
 * World-root diagnostics of a GC run (valid on world rank 0 only).  The
 * per-chain arrays hold the first MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS
 * chains; totals cover every chain.  final_sector_samples, when provided by
 * the caller with final_sector_capacity >= 2L+1, receives the number of
 * final samples per electron count summed over the chains; it must be
 * provided on every rank or on none.
 */
typedef struct {
  int valid;
  uint64_t chain_count;
  int coefficient_available;
  double complex coefficient[2];
  int coefficient_energy_available;
  double complex coefficient_energy;
  uint64_t bootstrap_attempts[MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS];
  uint64_t bootstrap_exact_zero[MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS];
  uint64_t bootstrap_numeric_zero[MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS];
  uint64_t bootstrap_local_energy[MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS];
  int bootstrap_exhausted[MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS];
  uint64_t total_bootstrap_attempts;
  uint64_t total_bootstrap_exact_zero;
  uint64_t total_bootstrap_numeric_zero;
  uint64_t total_bootstrap_local_energy;
  int any_bootstrap_exhausted;
  /* First failing chain (lowest index), if any. */
  int failure_recorded;
  uint64_t failure_chain;
  MVMCPowerLanczosGCStage failure_stage;
  MVMCKrylovStatus failure_status;
  int failure_electron_count;
  int failure_amplitude_state;
  double failure_log_abs;
  double failure_log_error_bound;
  uint64_t *final_sector_samples;
  size_t final_sector_capacity;
} MVMCPowerLanczosGCDiagnostics;

MVMCKrylovStatus mvmc_power_lanczos_gc_run(
    const MVMCPowerLanczosGCInput *input,
    MVMCPowerLanczosStabilizedResult *result);

/* diagnostics may be NULL; when non-NULL it is required on every rank. */
MVMCKrylovStatus mvmc_power_lanczos_gc_run_diagnosed(
    const MVMCPowerLanczosGCInput *input,
    MVMCPowerLanczosStabilizedResult *result,
    MVMCPowerLanczosGCDiagnostics *diagnostics);

const char *mvmc_power_lanczos_gc_stage_name(MVMCPowerLanczosGCStage stage);

#if defined(MVMC_ENABLE_FAULT_INJECTION)
/* Testing only: proposal/step counts of the last GC run on this rank. */
typedef struct {
  uint64_t generated[MVMC_POWER_LANCZOS_GC_STAGE_COUNT];
  uint64_t steps[MVMC_POWER_LANCZOS_GC_STAGE_COUNT];
  uint64_t fixed_sector_calls;
  uint64_t bootstrap_retries;
} MVMCPowerLanczosGCTestingCounters;

typedef enum {
  MVMC_POWER_LANCZOS_GC_FAULT_NONE = 0,
  /* The first final bootstrap candidate is a numeric zero (retried). */
  MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_NUMERIC,
  /* Every final bootstrap candidate is a numeric zero (1024 exhausted). */
  MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_EXHAUSTED,
  /* The first final-sampling proposal is a numeric zero (run fails). */
  MVMC_POWER_LANCZOS_GC_FAULT_FINAL_NUMERIC
} MVMCPowerLanczosGCFault;

void mvmc_power_lanczos_gc_testing_counters(
    MVMCPowerLanczosGCTestingCounters *counters);
void mvmc_power_lanczos_gc_testing_set_fault(MVMCPowerLanczosGCFault fault);
#endif

#endif /* MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE */

#endif /* MVMC_POWER_LANCZOS_STABILIZED_H */
