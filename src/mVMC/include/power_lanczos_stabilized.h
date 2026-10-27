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
 * Grand-canonical first-step entry (power step 1).  Without a finite-S
 * projection the coefficient chain starts from a uniform configuration with
 * initial_electron_count electrons.  With active S > 0 the start sectors
 * come from mvmc_gc_krylov_initial_sectors (an explicit N is binding and a
 * forbidden one is INVALID_ARGUMENT before any draw; auto tries the allowed
 * N nearest first) and each N gets at most 100 structured draws, redrawn
 * only when the base Psi is an exact or numeric zero.  Each chain searches
 * independently with its bootstrap RNG only; the chains meet before the
 * pilot.  The final chain tries initial-N (the selected N, structured, for
 * S > 0) and whole-domain uniform configurations alternately (at most 1024
 * attempts) and redraws only exact zero, numeric zero or unavailable local
 * energy of Phi.  Every stage uses the same symmetric mixture drawn by the
 * chain root (mvmc_gc_krylov_draw_mixture: global 1/16; finite S two-pair
 * 8/16, spin swap 2/16, local 5/16; S = 0 triplet two-pair 8/16, local
 * 7/16; otherwise local 15/16) and broadcast with the acceptance uniform.
 * initial_electron_count and initial_electron_count_explicit (0/1) must
 * agree on every world rank.
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
  int initial_electron_count_explicit;
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
 * chains; totals cover every chain.  The initial_* fields describe the
 * finite-S coefficient start (selected N, -1 if none or S = 0/inactive;
 * attempts and exact/numeric-zero redraws summed over the tried N; whether
 * every allowed N was exhausted).  final_sector_samples, when provided by
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
  int initial_selected_n[MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS];
  uint64_t initial_attempts[MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS];
  uint64_t initial_exact_zero[MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS];
  uint64_t initial_numeric_zero[MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS];
  int initial_exhausted[MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS];
  uint64_t total_initial_attempts;
  uint64_t total_initial_exact_zero;
  uint64_t total_initial_numeric_zero;
  int any_initial_exhausted;
  /* First failing chain (lowest index), if any. */
  int failure_recorded;
  uint64_t failure_chain;
  MVMCPowerLanczosGCStage failure_stage;
  MVMCKrylovStatus failure_status;
  int failure_retry_reason; /* MVMCKrylovBootstrapReason; never implies success */
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
  /* RNG draws consumed when each stage starts on this rank: the finite-S
   * start draws from the bootstrap RNG only. */
  uint64_t bootstrap_rng_draws_after_start;
  uint64_t coefficient_rng_draws_at_pilot;
  uint64_t final_rng_draws_at_bootstrap;
  uint64_t initial_retries;
  /* Final bootstrap candidates on the chain root: initial-N side, whole
   * domain side, initial-N candidates outside the selected finite-S sector
   * (N or fewer than 2S singles), and the accepted start (first word). */
  uint64_t bootstrap_initial_draws;
  uint64_t bootstrap_uniform_draws;
  uint64_t bootstrap_initial_off_sector;
  uint64_t final_start_word;
} MVMCPowerLanczosGCTestingCounters;

typedef enum {
  MVMC_POWER_LANCZOS_GC_FAULT_NONE = 0,
  /* The first final bootstrap candidate is a numeric zero (retried). */
  MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_NUMERIC,
  /* Every final bootstrap candidate is a numeric zero (1024 exhausted). */
  MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_EXHAUSTED,
  /* The first final-sampling proposal is a numeric zero (run fails). */
  MVMC_POWER_LANCZOS_GC_FAULT_FINAL_NUMERIC,
  /* Rank 1 disagrees before the refinement collective (run fails). */
  MVMC_POWER_LANCZOS_GC_FAULT_REFINEMENT_DECISION,
  MVMC_POWER_LANCZOS_GC_FAULT_REFINEMENT_STATUS,
  /* The first 1023 final bootstrap candidates are numeric zeros. */
  MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_LAST,
  /* The first final bootstrap candidate has no local energy (retried). */
  MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_LOCAL_ENERGY,
  /* The first final bootstrap candidate is nonfinite (run fails). */
  MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_NONFINITE,
  /* Route the first 128 finite-S final sampling proposals to refinement,
   * preserving the proposed configuration and its acceptance uniform. */
  MVMC_POWER_LANCZOS_GC_FAULT_REFINE_FINAL
} MVMCPowerLanczosGCFault;

/* Finite-S coefficient start faults.  EXACT_ZERO/NUMERIC_ZERO classify the
 * first count attempts of chain chain_index (UINT64_MAX: every chain) as
 * that zero (redrawn); NONFINITE ends the start at attempt index count;
 * STATUS makes chain rank 1 disagree at the first attempt. */
typedef enum {
  MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_NONE = 0,
  MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_EXACT_ZERO,
  MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_NUMERIC_ZERO,
  MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_NONFINITE,
  MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_STATUS
} MVMCPowerLanczosGCInitialFault;

void mvmc_power_lanczos_gc_testing_counters(
    MVMCPowerLanczosGCTestingCounters *counters);
void mvmc_power_lanczos_gc_testing_set_fault(MVMCPowerLanczosGCFault fault);
void mvmc_power_lanczos_gc_testing_set_initial_fault(
    MVMCPowerLanczosGCInitialFault fault, uint64_t count,
    uint64_t chain_index);

/*
 * Testing final-step harness: fixed raw coefficients (SYNTHETIC policy, no
 * scaled basis, no GEVP/pilot/bootstrap), a fixed start configuration and a
 * script of proposals/uniforms run through the same private final-step
 * helper as production (ordinary step, refinement decision agreement,
 * refinement and transaction).  force replaces the ordinary step of that
 * entry by a synthesized NUMERIC_ZERO or unavailable-local-energy
 * diagnostic, which then enters the production refinement path.
 * Collective on the world; the chain root consumes each entry once and
 * broadcasts it.  The global fault may only be NONE or REFINEMENT_DECISION/
 * STATUS and must agree on every rank; the testing counters are reset.
 */
#include "krylov_final_state.h"

typedef enum {
  MVMC_GC_FINAL_FORCE_NONE = 0,
  MVMC_GC_FINAL_FORCE_NUMERIC_ZERO = 1,
  MVMC_GC_FINAL_FORCE_LOCAL_ENERGY = 2
} MVMCGCFinalForce;

typedef struct {
  const uint64_t *proposal_words; /* word_count elements */
  double uniform;                 /* finite, 0 <= u < 1 */
  MVMCGCFinalForce force;
} MVMCGCFinalScriptStep;

typedef struct {
  const MVMCPowerLanczosGCView *gc_view;
  MVMCClassicPfaffianCommunicator world_communicator, chain_communicator;
  double complex coefficient[2]; /* raw c0, c1; order 1 */
  const uint64_t *initial_words;
  size_t word_count, step_count;
  const MVMCGCFinalScriptStep *steps;
} MVMCGCFinalScript;

/* One step; before/after_words are caller buffers of word_count words that
 * alias neither the script nor another buffer.  Counters are cumulative;
 * precision/node_calls are of this step (0 without refinement).  N and N^2
 * are of the current configuration after the step. */
typedef struct MVMCGCFinalStepTrace {
  uint64_t *before_words;
  uint64_t *after_words;
  MVMCKrylovFinalStateSnapshot before;
  MVMCKrylovFinalStateSnapshot after;
  MVMCKrylovStatus status;
  int accepted;
  int used_refinement;
  int exact_zero_proposal;
  uint64_t refinement_calls;
  uint64_t refinement_accepted;
  uint64_t refinement_exact_zero;
  unsigned precision;
  uint64_t node_calls;
  uint64_t proposal_copies;
  uint64_t uniform_copies;
  uint64_t rng_draws;
  int electron_count;
  uint64_t electron_count_square;
} MVMCGCFinalStepTrace;

/* Invalid input (pointers, shapes, capacity, enums, fault mode, script
 * agreement) is INVALID_ARGUMENT on every rank before any step, with the
 * trace untouched and *written = 0.  A failing step is recorded (with its
 * unchanged snapshots) and ends the script; *written includes it. */
MVMCKrylovStatus mvmc_power_lanczos_gc_testing_final_steps(
    const MVMCGCFinalScript *script, MVMCGCFinalStepTrace *trace,
    size_t trace_capacity, size_t *written);
#endif

#endif /* MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE */

#endif /* MVMC_POWER_LANCZOS_STABILIZED_H */
