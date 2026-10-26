#include "power_lanczos_stabilized.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE) ||                     \
    !defined(MVMC_ENABLE_POWER_LANCZOS_P5_CORE)
#error "power_lanczos_stabilized.c requires the production Krylov core"
#endif

#include "gc_krylov_proposal.h"
#include "krylov_final_state_chain.h"
#include "krylov_gevp_solver.h"
#include "krylov_matrix_measurement.h"
#include "krylov_streaming_statistics.h"
#include "power_lanczos_sampling.h"

#include <complex.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(MVMC_ENABLE_FAULT_INJECTION)
#include <stdio.h>
#include <sys/time.h>
#endif

enum {
  MVMC_POWER_LANCZOS_PACKED_CAPACITY = 6,
  MVMC_POWER_LANCZOS_BOOTSTRAP_ATTEMPTS = 1024,
  MVMC_POWER_LANCZOS_INITIAL_ATTEMPTS = 100
};

typedef struct {
  uint64_t count;
  double complex overlap[MVMC_POWER_LANCZOS_PACKED_CAPACITY];
  double complex hamiltonian[MVMC_POWER_LANCZOS_PACKED_CAPACITY];
  double complex hamiltonian_adjoint[MVMC_POWER_LANCZOS_PACKED_CAPACITY];
  double complex hamiltonian_squared[MVMC_POWER_LANCZOS_PACKED_CAPACITY];
} MVMCPowerLanczosMatrixBlock;

static void invalidate_result(
    MVMCKrylovStatus status, MVMCPowerLanczosStabilizedResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->version = MVMC_POWER_LANCZOS_STABILIZED_VERSION;
  result->status = status;
  result->energy = NAN;
  result->energy_standard_error = NAN;
  result->variance = NAN;
  result->variance_standard_error = NAN;
  result->final_energy_imaginary = NAN;
  result->variance_imaginary = NAN;
  result->condition_estimate = NAN;
  result->gevp_residual = NAN;
  result->energy_tau_int = NAN;
  result->effective_sample_count = NAN;
}

static MVMCKrylovStatus synchronize_status(
    MVMCKrylovStatus local, MVMCClassicPfaffianCommunicator communicator) {
#ifdef _mpi_use
  int encoded = (int)local;
  int reduced = 0;
  if (MPI_Allreduce(&encoded, &reduced, 1, MPI_INT, MPI_MAX,
                    communicator) != MPI_SUCCESS ||
      reduced < (int)MVMC_KRYLOV_STATUS_OK ||
      reduced > (int)MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  return (MVMCKrylovStatus)reduced;
#else
  (void)communicator;
  return local;
#endif
}

static uint64_t mix_u64(uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

static uint64_t domain_seed(uint64_t seed, uint64_t domain) {
  const uint64_t derived = mix_u64(seed ^ domain);
  return derived == 0 ? UINT64_C(1) : derived;
}

static uint64_t chain_domain(uint64_t domain, uint64_t chain_index) {
  return domain ^
         (chain_index * UINT64_C(0x9e3779b97f4a7c15));
}

static uint64_t hash_u64(uint64_t hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    hash ^= (value >> (8 * byte)) & UINT64_C(0xff);
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

static MVMCKrylovStatus map_gevp_status(MVMCKrylovGEVPStatus status) {
  switch (status) {
    case MVMC_KRYLOV_GEVP_OK:
      return MVMC_KRYLOV_STATUS_OK;
    case MVMC_KRYLOV_GEVP_INVALID_ARGUMENT:
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    case MVMC_KRYLOV_GEVP_NONFINITE_INPUT:
      return MVMC_KRYLOV_STATUS_NONFINITE;
    case MVMC_KRYLOV_GEVP_LAPACK_FAILURE:
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    default:
      return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
}

static int checked_add_u64(uint64_t left, uint64_t right, uint64_t *sum) {
  if (sum == NULL || right > UINT64_MAX - left) return 0;
  *sum = left + right;
  return 1;
}

static MVMCKrylovStatus create_limits(
    const MVMCKrylovFockModel *model, int power_step,
    uint64_t generation_hash, MVMCKrylovBoundedLimits *limits) {
  const int maximum_order = power_step + 1;
  uint64_t frontier = 1;
  uint64_t nodes = 1;
  uint64_t transitions = 0;
  int depth;
  if (model == NULL || limits == NULL || model->term_count == 0 ||
      power_step < 1 || power_step > 2 || generation_hash == 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (depth = 0; depth < maximum_order; ++depth) {
    uint64_t next;
    if ((uint64_t)model->term_count > UINT64_MAX / frontier) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    next = frontier * (uint64_t)model->term_count;
    if (!checked_add_u64(nodes, next, &nodes) ||
        !checked_add_u64(transitions, next, &transitions)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    frontier = next;
  }
  if (nodes > UINT64_C(10000000) ||
      transitions > UINT64_C(100000000)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  memset(limits, 0, sizeof(*limits));
  limits->amplitude_policy_hash = generation_hash;
  limits->cache_bytes = (size_t)4 * 1024 * 1024;
  limits->max_row_transitions = model->term_count;
  limits->max_workspace_bytes = (size_t)512 * 1024 * 1024;
  limits->max_node_expansions = nodes;
  limits->max_terminal_amplitude_calls = frontier;
  limits->max_total_row_transitions = transitions;
  limits->max_order = maximum_order;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus draw_bounded(
    void *context, size_t bound, size_t *value) {
  return mvmc_krylov_positive_sampler_rng_draw_bounded(
      (MVMCKrylovPositiveSamplerRng *)context, bound, value);
}

static MVMCKrylovStatus draw_uniform_sector(
    const MVMCKrylovFockModel *model,
    MVMCClassicPfaffianCommunicator communicator,
    MVMCKrylovPositiveSamplerRng *rng, uint64_t *words,
    size_t word_count) {
  MVMCKrylovFockUniformProposalResult proposal;
  MVMCKrylovStatus status;
  memset(&proposal, 0, sizeof(proposal));
  status = mvmc_krylov_fock_proposal_draw_uniform_sector(
      model, draw_bounded, rng, words, word_count, &proposal);
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (!proposal.valid || proposal.status != MVMC_KRYLOV_STATUS_OK ||
       proposal.word_count != word_count)) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  status = synchronize_status(status, communicator);
#ifdef _mpi_use
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (word_count > (size_t)INT_MAX ||
       MPI_Bcast(words, (int)word_count, MPI_UINT64_T, 0,
                 communicator) != MPI_SUCCESS)) {
    status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#endif
  return synchronize_status(status, communicator);
}

static MVMCKrylovStatus coefficient_step(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *guide,
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    uint64_t *words, size_t word_count,
    MVMCKrylovPositiveSamplerSnapshot *snapshot,
    MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words) {
  MVMCKrylovPositiveSamplerProposalStepResult step;
  return mvmc_krylov_positive_sampler_step_mixture_rng(
      workspace, guide, model, proposal_policy, words, word_count, snapshot,
      rng, amplitude, amplitude_context, proposal_words, word_count, &step);
}

static MVMCKrylovStatus final_step(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy,
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    uint64_t *words, size_t word_count,
    MVMCKrylovFinalStateSnapshot *snapshot,
    MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words) {
  MVMCKrylovFinalStateChainStepResult step;
  return mvmc_krylov_final_state_chain_step_mixture_rng(
      workspace, policy, model, proposal_policy, words, word_count,
      snapshot, rng, amplitude, amplitude_context, proposal_words,
      word_count, &step);
}

/* One tried finite-S start sector of a chain. */
typedef struct {
  uint64_t attempts;
  uint64_t exact_zero;
  uint64_t numeric_zero;
  int selected;
  MVMCKrylovStatus status;
} MVMCPowerLanczosGCInitialRecord;

/*
 * Grand-canonical run state.  NULL selects the canonical fixed-sector
 * proposals; otherwise the chain root draws the symmetric GC proposal and
 * the acceptance uniform and broadcasts them before every candidate step.
 */
typedef struct {
  const MVMCKrylovFockModel *model;
  /* Reused transaction buffer; current is committed after chain agreement. */
  uint64_t *step_words;
  MVMCClassicPfaffianCommunicator chain_communicator;
  int topology_ready;
  int chain_rank;
  uint64_t chain_index;
  MVMCPowerLanczosGCStage stage;
  int initial_electron_count;
  int initial_electron_count_explicit;
  MVMCGCSpinSupport spin_support;
  /* Finite-S coefficient start: allowed sectors and per-N records. */
  int *initial_sectors;
  size_t initial_sector_count;
  size_t initial_tried;
  MVMCPowerLanczosGCInitialRecord *initial_records;
  int initial_selected_n;
  int initial_exhausted;
  uint64_t bootstrap_attempts;
  uint64_t bootstrap_exact_zero;
  uint64_t bootstrap_numeric_zero;
  uint64_t bootstrap_local_energy;
  uint64_t refinement_calls, refinement_accepted, refinement_exact_zero;
  unsigned refinement_max_precision;
  int bootstrap_exhausted;
  int failure_recorded;
  MVMCPowerLanczosGCStage failure_stage;
  MVMCKrylovStatus failure_status;
  MVMCKrylovBootstrapReason failure_retry_reason;
  int failure_electron_count;
  int failure_amplitude_state;
  double failure_log_abs;
  double failure_log_error_bound;
  int coefficient_available;
  double complex coefficient[2];
  int coefficient_energy_available;
  double complex coefficient_energy;
  uint64_t *sector_samples;
  size_t sector_capacity;
#if defined(MVMC_ENABLE_FAULT_INJECTION)
  int final_sample_seen;
  int b0_enabled, b0_io_failed, b0_kind, b0_different;
  FILE *b0_trace;
  double b0_started;
  double b0_seconds[MVMC_POWER_LANCZOS_GC_STAGE_COUNT];
  uint64_t b0_counts[MVMC_POWER_LANCZOS_GC_STAGE_COUNT][4][4];
#endif
} MVMCPowerLanczosGCRun;

/* Active S > 0 projection: structured start and finite-S proposals. */
static int gc_finite_spin(const MVMCPowerLanczosGCRun *gc) {
  return gc != NULL && gc->spin_support.active &&
         gc->spin_support.total_spin > 0;
}

#if defined(MVMC_ENABLE_FAULT_INJECTION)
static MVMCPowerLanczosGCTestingCounters gc_testing_counters;
static MVMCPowerLanczosGCFault gc_testing_fault =
    MVMC_POWER_LANCZOS_GC_FAULT_NONE;
static MVMCPowerLanczosGCInitialFault gc_testing_initial_fault =
    MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_NONE;
static uint64_t gc_testing_initial_count = 0;
static uint64_t gc_testing_initial_chain = 0;

void mvmc_power_lanczos_gc_testing_counters(
    MVMCPowerLanczosGCTestingCounters *counters) {
  if (counters != NULL) *counters = gc_testing_counters;
}

void mvmc_power_lanczos_gc_testing_set_fault(MVMCPowerLanczosGCFault fault) {
  gc_testing_fault = fault;
}

void mvmc_power_lanczos_gc_testing_set_initial_fault(
    MVMCPowerLanczosGCInitialFault fault, uint64_t count,
    uint64_t chain_index) {
  gc_testing_initial_fault = fault;
  gc_testing_initial_count = count;
  gc_testing_initial_chain = chain_index;
}
#define GC_COUNT(field) (++gc_testing_counters.field)
#else
#define GC_COUNT(field) ((void)0)
#endif

const char *mvmc_power_lanczos_gc_stage_name(MVMCPowerLanczosGCStage stage) {
  switch (stage) {
    case MVMC_POWER_LANCZOS_GC_STAGE_NONE:
      return "none";
    case MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_INITIAL:
      return "coefficient-initial";
    case MVMC_POWER_LANCZOS_GC_STAGE_PILOT_WARM_UP:
      return "pilot-warm-up";
    case MVMC_POWER_LANCZOS_GC_STAGE_PILOT_SAMPLE:
      return "pilot-sample";
    case MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_WARM_UP:
      return "coefficient-warm-up";
    case MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_SAMPLE:
      return "coefficient-sample";
    case MVMC_POWER_LANCZOS_GC_STAGE_FINAL_BOOTSTRAP:
      return "final-bootstrap";
    case MVMC_POWER_LANCZOS_GC_STAGE_FINAL_WARM_UP:
      return "final-warm-up";
    case MVMC_POWER_LANCZOS_GC_STAGE_FINAL_SAMPLE:
      return "final-sample";
    case MVMC_POWER_LANCZOS_GC_STAGE_COUNT:
      break;
  }
  return "invalid";
}

static int popcount_words(const uint64_t *words, size_t word_count) {
  int count = 0;
  size_t index;
  for (index = 0; index < word_count; ++index) {
    uint64_t value = words[index];
    while (value != 0) {
      count += (int)(value & UINT64_C(1));
      value >>= 1;
    }
  }
  return count;
}

/* B0 experiment only: no new RNG draws, amplitudes or collectives. */
#if defined(MVMC_ENABLE_FAULT_INJECTION)
static double gc_b0_now(void) {
  struct timeval now;
  if (gettimeofday(&now, NULL) != 0) return 0.0;
  return (double)now.tv_sec + (double)now.tv_usec * 1e-6;
}
static void gc_b0_stage(MVMCPowerLanczosGCRun *gc,
                        MVMCPowerLanczosGCStage stage) {
  if (gc->b0_enabled && gc->chain_rank == 0) {
    const double now = gc_b0_now();
    if (gc->b0_started != 0.0)
      gc->b0_seconds[gc->stage] += now - gc->b0_started;
    gc->b0_started = now;
  }
  gc->stage = stage;
}
static void gc_b0_sample(MVMCPowerLanczosGCRun *gc, size_t sample,
                         const uint64_t *words, size_t word_count) {
  if (gc == NULL || !gc->b0_enabled || gc->chain_rank != 0 || gc->b0_io_failed)
    return;
  if (gc->b0_trace == NULL) {
    char name[96];
    snprintf(name, sizeof(name), "gc-b0-lanczos.chain%llu",
             (unsigned long long)gc->chain_index);
    gc->b0_trace = fopen(name, "w");
    if (gc->b0_trace == NULL) { gc->b0_io_failed = 1; return; }
  }
  if (fprintf(gc->b0_trace, "S %s %zu %d\n",
              mvmc_power_lanczos_gc_stage_name(gc->stage), sample,
              popcount_words(words, word_count)) < 0)
    gc->b0_io_failed = 1;
}
static void gc_b0_step(MVMCPowerLanczosGCRun *gc, int accepted) {
  if (gc->b0_enabled && gc->chain_rank == 0) {
    uint64_t *count = gc->b0_counts[gc->stage][gc->b0_kind];
    ++count[0];
    count[1] += gc->b0_different;
    count[2] += accepted != 0;
    count[3] += accepted != 0 && gc->b0_different;
  }
}
static void gc_b0_finish(MVMCPowerLanczosGCRun *gc, MVMCKrylovStatus status) {
  int stage, kind;
  if (!gc->b0_enabled || gc->chain_rank != 0) return;
  gc_b0_stage(gc, MVMC_POWER_LANCZOS_GC_STAGE_NONE);
  if (gc->b0_trace != NULL) {
    for (stage = 1; stage < MVMC_POWER_LANCZOS_GC_STAGE_COUNT; ++stage) {
      if (fprintf(gc->b0_trace, "T %s %.9f\n",
                  mvmc_power_lanczos_gc_stage_name((MVMCPowerLanczosGCStage)stage),
                  gc->b0_seconds[stage]) < 0) gc->b0_io_failed = 1;
      /* kind 3 (swap) exists for finite S only; S = 0 keeps 3 rows. */
      for (kind = 0; kind < (gc_finite_spin(gc) ? 4 : 3); ++kind) {
        const uint64_t *c = gc->b0_counts[stage][kind];
        if (fprintf(gc->b0_trace, "C %s %d %llu %llu %llu %llu\n",
                    mvmc_power_lanczos_gc_stage_name((MVMCPowerLanczosGCStage)stage), kind,
                    (unsigned long long)c[0], (unsigned long long)c[1],
                    (unsigned long long)c[2], (unsigned long long)c[3]) < 0)
          gc->b0_io_failed = 1;
      }
    }
    if (!gc->b0_io_failed)
      fprintf(gc->b0_trace, "END %d\n", (int)status);
    if (fclose(gc->b0_trace) != 0) gc->b0_io_failed = 1;
    gc->b0_trace = NULL;
  }
  if (gc->b0_io_failed) fprintf(stderr, "GC B0 trace I/O failure\n");
}
#else
#define gc_b0_stage(gc, value) ((gc)->stage = (value))
#define gc_b0_sample(gc, sample, words, count) ((void)0)
#define gc_b0_step(gc, accepted) ((void)0)
#endif

/* The first failure of a chain: stage, candidate N and classification. */
static void gc_record_failure(MVMCPowerLanczosGCRun *gc,
                              MVMCKrylovStatus status,
                              const uint64_t *words, size_t word_count,
                              const MVMCKrylovFinalDiagnostic *diagnostic) {
  if (gc == NULL || gc->failure_recorded ||
      status == MVMC_KRYLOV_STATUS_OK) {
    return;
  }
  gc->failure_recorded = 1;
  gc->failure_stage = gc->stage;
  gc->failure_status = status;
  gc->failure_retry_reason = MVMC_KRYLOV_BOOTSTRAP_NO_RETRY;
  gc->failure_electron_count =
      words != NULL ? popcount_words(words, word_count) : -1;
  gc->failure_amplitude_state = -1;
  gc->failure_log_abs = NAN;
  gc->failure_log_error_bound = NAN;
  if (diagnostic != NULL && diagnostic->evaluation_available) {
    gc->failure_retry_reason = diagnostic->retry_reason;
    gc->failure_amplitude_state = (int)diagnostic->evaluation.amplitude.state;
    gc->failure_log_abs = diagnostic->evaluation.amplitude.log_abs;
    gc->failure_log_error_bound =
        diagnostic->evaluation.amplitude.log_abs_error_bound;
  }
}

/* Root status first, then the candidate words and the uniform. */
static MVMCKrylovStatus gc_broadcast(const MVMCPowerLanczosGCRun *gc,
                                     MVMCKrylovStatus root_status,
                                     uint64_t *words, size_t word_count,
                                     double *uniform) {
#ifdef _mpi_use
  int encoded = (int)root_status;
  if (MPI_Bcast(&encoded, 1, MPI_INT, 0, gc->chain_communicator) !=
      MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (encoded != (int)MVMC_KRYLOV_STATUS_OK) {
    return (MVMCKrylovStatus)encoded;
  }
  if (word_count > (size_t)INT_MAX ||
      MPI_Bcast(words, (int)word_count, MPI_UINT64_T, 0,
                gc->chain_communicator) != MPI_SUCCESS ||
      (uniform != NULL &&
       MPI_Bcast(uniform, 1, MPI_DOUBLE, 0, gc->chain_communicator) !=
           MPI_SUCCESS)) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  return MVMC_KRYLOV_STATUS_OK;
#else
  (void)gc;
  (void)words;
  (void)word_count;
  (void)uniform;
  return root_status;
#endif
}

/* Chain root: the fixed symmetric mixture of mvmc_gc_krylov_draw_mixture
 * (global 1/16; finite S: two-pair 8/16, swap 2/16, local 5/16; S = 0
 * triplet: two-pair 8/16, local 7/16; otherwise local 15/16), then the
 * acceptance uniform.  log ratio is 0. */
static MVMCKrylovStatus gc_propose(MVMCPowerLanczosGCRun *gc,
                                   MVMCKrylovPositiveSamplerRng *rng,
                                   const uint64_t *current,
                                   uint64_t *proposal, size_t word_count,
                                   double *uniform) {
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  *uniform = 0.0;
  if (gc->chain_rank == 0) {
    MVMCGCMoveKind kind = MVMC_GC_MOVE_LOCAL;
    status = mvmc_gc_krylov_draw_mixture(gc->model, &gc->spin_support,
                                         current, draw_bounded, rng, proposal,
                                         word_count, &kind);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_krylov_positive_sampler_rng_draw_uniform(rng, uniform);
    }
#if defined(MVMC_ENABLE_FAULT_INJECTION)
    if (gc->b0_enabled) {
      gc->b0_kind = (int)kind;
      gc->b0_different = status == MVMC_KRYLOV_STATUS_OK &&
          memcmp(current, proposal, word_count * sizeof(*current)) != 0;
    }
    if (gc->stage < MVMC_POWER_LANCZOS_GC_STAGE_COUNT) {
      ++gc_testing_counters.generated[gc->stage];
    }
#endif
  }
  return gc_broadcast(gc, status, proposal, word_count, uniform);
}

/* Chain root: a start configuration (initial N or the whole domain).  With
 * S > 0 the initial-N side is the structured draw at the N selected by the
 * coefficient start of this chain. */
static MVMCKrylovStatus gc_draw_start(MVMCPowerLanczosGCRun *gc,
                                      MVMCKrylovPositiveSamplerRng *rng,
                                      uint64_t *words, size_t word_count,
                                      int initial_n) {
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  if (gc->chain_rank == 0) {
    if (!initial_n) {
      status = mvmc_gc_krylov_draw_uniform(gc->model, draw_bounded, rng,
                                           words, word_count);
    } else if (gc_finite_spin(gc)) {
      status = mvmc_gc_krylov_draw_spin_initial(
          gc->model, &gc->spin_support, gc->initial_selected_n, draw_bounded,
          rng, words, word_count);
    } else {
      status = mvmc_gc_krylov_draw_initial(gc->model,
                                           gc->initial_electron_count,
                                           draw_bounded, rng, words,
                                           word_count);
    }
#if defined(MVMC_ENABLE_FAULT_INJECTION)
    if (gc->stage < MVMC_POWER_LANCZOS_GC_STAGE_COUNT) {
      ++gc_testing_counters.generated[gc->stage];
    }
#endif
  }
  return gc_broadcast(gc, status, words, word_count, NULL);
}

static MVMCKrylovStatus coefficient_step_any(
    MVMCPowerLanczosGCRun *gc, MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *guide,
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    uint64_t *words, size_t word_count,
    MVMCKrylovPositiveSamplerSnapshot *snapshot,
    MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words) {
  MVMCKrylovPositiveSamplerStepResult step;
  double uniform = 0.0;
  MVMCKrylovStatus status;
  if (gc == NULL) {
    GC_COUNT(fixed_sector_calls);
    return coefficient_step(workspace, guide, model, proposal_policy, words,
                            word_count, snapshot, rng, amplitude,
                            amplitude_context, proposal_words);
  }
  status = gc_propose(gc, rng, words, proposal_words, word_count, &uniform);
  if (status == MVMC_KRYLOV_STATUS_OK) {
#if defined(MVMC_ENABLE_FAULT_INJECTION)
    ++gc_testing_counters.steps[gc->stage];
#endif
    status = mvmc_krylov_positive_sampler_step(
        workspace, guide, words, word_count, snapshot, proposal_words,
        word_count, 0.0, uniform, amplitude, amplitude_context, &step);
    if (status == MVMC_KRYLOV_STATUS_OK) gc_b0_step(gc, step.accepted);
  }
  gc_record_failure(gc, status, proposal_words, word_count, NULL);
  return status;
}

/* Agree before entering a conditional collective.  Return OK only when
 * every rank has the same status and decision, even for an error status. */
static MVMCKrylovStatus gc_agree_value(const MVMCPowerLanczosGCRun *gc,
                                 MVMCKrylovStatus status,
                                 int value) {
#ifdef _mpi_use
  int local[2];
  int minimum[2];
  int maximum[2];
  local[0] = (int)status;
  local[1] = value;
  if (MPI_Allreduce(local, minimum, 2, MPI_INT, MPI_MIN,
                    gc->chain_communicator) != MPI_SUCCESS ||
      MPI_Allreduce(local, maximum, 2, MPI_INT, MPI_MAX,
                    gc->chain_communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (minimum[0] != maximum[0] || minimum[1] != maximum[1]) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
#else
  (void)gc;
  (void)value;
  (void)status;
#endif
  return MVMC_KRYLOV_STATUS_OK;
}

/* Decision of one final proposal step (harness trace). */
typedef struct {
  int accepted;
  int used_refinement;
  int exact_zero_proposal;
  unsigned precision;
  uint64_t node_calls;
} MVMCPowerLanczosGCStepOutcome;

#if defined(MVMC_ENABLE_FAULT_INJECTION)
/* Testing hook: a synthesized ordinary-step failure with a valid
 * evaluation (NUMERIC_ZERO amplitude, or a finite amplitude without a local
 * energy).  words/snapshot are not touched. */
static MVMCKrylovStatus gc_forced_diagnostic(
    const MVMCKrylovFinalStatePolicy *policy, MVMCGCFinalForce force,
    MVMCKrylovFinalDiagnostic *diagnostic) {
  MVMCKrylovFinalStateEvaluation *evaluation = &diagnostic->evaluation;
  evaluation->status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  evaluation->order = policy->order;
  evaluation->policy_hash = policy->policy_hash;
  if (mvmc_scaled_complex_make_exact_zero(&evaluation->hamiltonian_amplitude) !=
          MVMC_PFAFFIAN_STATUS_OK ||
      mvmc_scaled_complex_make_exact_zero(&evaluation->local_energy) !=
          MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (force == MVMC_GC_FINAL_FORCE_NUMERIC_ZERO) {
    if (mvmc_scaled_complex_make_numeric_zero(-30.0, -INFINITY, -INFINITY,
                                              0.0, &evaluation->amplitude) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    diagnostic->retry_reason = MVMC_KRYLOV_BOOTSTRAP_NUMERIC_ZERO;
  } else {
    if (mvmc_scaled_complex_make_finite(1.0 + 0.0 * I, 0.0, -INFINITY,
                                        &evaluation->amplitude) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    evaluation->sampleable = 1;
    evaluation->local_energy_available = 0;
    diagnostic->retry_reason = MVMC_KRYLOV_BOOTSTRAP_LOCAL_ENERGY;
  }
  diagnostic->evaluation_available = 1;
  diagnostic->status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
}
#endif

/*
 * One final-stage step for a broadcast proposal and uniform, shared by
 * production (force always 0) and the Testing final-step harness: the
 * ordinary diagnosed step (replaced by a synthesized diagnostic when the
 * harness forces one), the chain agreement on status and refinement entry,
 * the exact-input refinement with its transactional snapshot replacement,
 * and the B0 step count.  No world collective.
 */
static MVMCKrylovStatus gc_final_proposal_step(
    MVMCPowerLanczosGCRun *gc, MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy,
    const MVMCKrylovFockModel *model, uint64_t *words, size_t word_count,
    MVMCKrylovFinalStateSnapshot *snapshot,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    const uint64_t *proposal_words, double uniform, int force,
    MVMCKrylovFinalDiagnostic *diagnostic,
    MVMCPowerLanczosGCStepOutcome *outcome) {
  MVMCKrylovFinalStateStepResult step;
  MVMCKrylovFinalStateSnapshot pending = *snapshot;
  MVMCKrylovStatus status;
  memcpy(gc->step_words, words, word_count * sizeof(*words));
  memset(&step, 0, sizeof(step));
  memset(diagnostic, 0, sizeof(*diagnostic));
  memset(outcome, 0, sizeof(*outcome));
#if defined(MVMC_ENABLE_FAULT_INJECTION)
  ++gc_testing_counters.steps[gc->stage];
  if (force != MVMC_GC_FINAL_FORCE_NONE) {
    /* Synthesized on every rank instead of the ordinary step. */
    status = gc_forced_diagnostic(policy, (MVMCGCFinalForce)force,
                                  diagnostic);
  } else
#else
  (void)force;
#endif
  {
    status = mvmc_krylov_final_state_sampler_step_diagnosed(
        workspace, policy, gc->step_words, word_count, &pending, proposal_words,
        word_count, 0.0, uniform, amplitude, amplitude_context, &step,
        diagnostic);
  }
  outcome->exact_zero_proposal =
      status == MVMC_KRYLOV_STATUS_OK && step.acceptance.exact_zero_proposal;
#if !defined(MVMC_HAVE_MPFR_REFINEMENT)
  (void)model;
#if !defined(MVMC_ENABLE_FAULT_INJECTION)
  (void)gc;
#endif
#else
  {
    /* The failed ordinary step is transactional. Reuse its proposal and
     * uniform; do not redraw or reinterpret NUMERIC_ZERO as a node.
     * PFAPACK singularity alone is not an exact-zero proof either. */
    int needs_refinement =
        (status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
         diagnostic->evaluation_available &&
         (diagnostic->retry_reason == MVMC_KRYLOV_BOOTSTRAP_NUMERIC_ZERO ||
          diagnostic->retry_reason == MVMC_KRYLOV_BOOTSTRAP_LOCAL_ENERGY)) ||
        (status == MVMC_KRYLOV_STATUS_OK &&
         step.acceptance.exact_zero_proposal);
    MVMCKrylovStatus agreement;
#if defined(MVMC_ENABLE_FAULT_INJECTION)
    if (gc->chain_rank == 1) {
      if (gc_testing_fault == MVMC_POWER_LANCZOS_GC_FAULT_REFINEMENT_DECISION)
        needs_refinement = !needs_refinement;
      if (gc_testing_fault == MVMC_POWER_LANCZOS_GC_FAULT_REFINEMENT_STATUS)
        status = MVMC_KRYLOV_STATUS_NONFINITE;
    }
#endif
    agreement = gc_agree_value(gc, status, needs_refinement);
    if (agreement != MVMC_KRYLOV_STATUS_OK) {
      status = agreement;
      needs_refinement = 0;
    }
    if (needs_refinement) {
      MVMCKrylovRefinementResult refined;
      MVMCKrylovFinalStateSnapshot updated;
      outcome->used_refinement = 1;
      status = mvmc_gc_krylov_refine_final(amplitude_context, model, policy,
          words, proposal_words, word_count, uniform, &refined);
      ++gc->refinement_calls;
      if (status != MVMC_KRYLOV_STATUS_OK && gc->chain_rank == 0)
        fprintf(stderr, "GC Lanczos refinement failed: %s (chain=%llu).\n",
                mvmc_krylov_status_string(status),
                (unsigned long long)gc->chain_index);
      if (status == MVMC_KRYLOV_STATUS_OK) {
        gc->refinement_accepted += refined.accepted != 0;
        gc->refinement_exact_zero += refined.exact_zero_proposal != 0;
        if (refined.precision > gc->refinement_max_precision)
          gc->refinement_max_precision = refined.precision;
        outcome->exact_zero_proposal = refined.exact_zero_proposal != 0;
        outcome->precision = refined.precision;
        outcome->node_calls = refined.node_calls;
      }
      if (status == MVMC_KRYLOV_STATUS_OK && refined.accepted) {
        status = snapshot->accepted_generation == UINT64_MAX
            ? MVMC_KRYLOV_STATUS_RESOURCE_LIMIT
            : mvmc_krylov_final_state_snapshot_from_evaluation(policy,
                  proposal_words, word_count, &refined.krylov,
                  &refined.evaluation, &updated);
        status = synchronize_status(status, gc->chain_communicator);
        if (status == MVMC_KRYLOV_STATUS_OK) {
          updated.accepted_generation = snapshot->accepted_generation + 1;
          memcpy(gc->step_words, proposal_words, word_count * sizeof(*words));
          pending = updated;
        }
      }
      if (status == MVMC_KRYLOV_STATUS_OK) step.accepted = refined.accepted;
    }
  }
#endif
  if (status == MVMC_KRYLOV_STATUS_OK) {
    memcpy(words, gc->step_words, word_count * sizeof(*words));
    *snapshot = pending;
    gc_b0_step(gc, step.accepted);
    outcome->accepted = step.accepted;
  }
  return status;
}

static MVMCKrylovStatus final_step_any(
    MVMCPowerLanczosGCRun *gc, MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy,
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    uint64_t *words, size_t word_count,
    MVMCKrylovFinalStateSnapshot *snapshot,
    MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words) {
  MVMCKrylovFinalDiagnostic diagnostic;
  MVMCPowerLanczosGCStepOutcome outcome;
  double uniform = 0.0;
  MVMCKrylovStatus status;
  if (gc == NULL) {
    GC_COUNT(fixed_sector_calls);
    return final_step(workspace, policy, model, proposal_policy, words,
                      word_count, snapshot, rng, amplitude,
                      amplitude_context, proposal_words);
  }
  memset(&diagnostic, 0, sizeof(diagnostic));
  status = gc_propose(gc, rng, words, proposal_words, word_count, &uniform);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = gc_final_proposal_step(gc, workspace, policy, model, words,
                                    word_count, snapshot, amplitude,
                                    amplitude_context, proposal_words,
                                    uniform, 0, &diagnostic, &outcome);
#if defined(MVMC_ENABLE_FAULT_INJECTION)
    /* A genuine numeric zero with a finite error bound, never EXACT_ZERO.
     * Kept outside the shared helper (the harness never reaches it). */
    if (status == MVMC_KRYLOV_STATUS_OK &&
        gc_testing_fault == MVMC_POWER_LANCZOS_GC_FAULT_FINAL_NUMERIC &&
        gc->stage == MVMC_POWER_LANCZOS_GC_STAGE_FINAL_SAMPLE &&
        !gc->final_sample_seen &&
        mvmc_scaled_complex_make_numeric_zero(
            -30.0, -INFINITY, -INFINITY, 0.0,
            &diagnostic.evaluation.amplitude) == MVMC_PFAFFIAN_STATUS_OK) {
      diagnostic.evaluation_available = 1;
      diagnostic.retry_reason = MVMC_KRYLOV_BOOTSTRAP_NUMERIC_ZERO;
      status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    }
    if (gc->stage == MVMC_POWER_LANCZOS_GC_STAGE_FINAL_SAMPLE) {
      gc->final_sample_seen = 1;
    }
#endif
  }
  gc_record_failure(gc, status, proposal_words, word_count, &diagnostic);
  return status;
}

/* Bootstrap reasons and statuses must also agree inside the chain. */
static MVMCKrylovStatus gc_agree(const MVMCPowerLanczosGCRun *gc,
                                MVMCKrylovStatus status,
                                MVMCKrylovBootstrapReason reason) {
  MVMCKrylovStatus agreement = gc_agree_value(gc, status, (int)reason);
  return agreement == MVMC_KRYLOV_STATUS_OK ? status : agreement;
}

/* World agreement of the start request; called on every rank. */
static MVMCKrylovStatus gc_agree_inputs(
    const MVMCPowerLanczosGCRun *gc,
    MVMCClassicPfaffianCommunicator world_communicator) {
  MVMCKrylovStatus status =
      gc->initial_electron_count_explicit == 0 ||
              gc->initial_electron_count_explicit == 1
          ? MVMC_KRYLOV_STATUS_OK
          : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
#ifdef _mpi_use
  {
    int local[2];
    int minimum[2] = {0, 0};
    int maximum[2] = {0, 0};
    local[0] = gc->initial_electron_count;
    local[1] = gc->initial_electron_count_explicit;
    if (MPI_Allreduce(local, minimum, 2, MPI_INT, MPI_MIN,
                      world_communicator) != MPI_SUCCESS ||
        MPI_Allreduce(local, maximum, 2, MPI_INT, MPI_MAX,
                      world_communicator) != MPI_SUCCESS) {
      return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
    if (minimum[0] != maximum[0] || minimum[1] != maximum[1]) {
      status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
#endif
  return synchronize_status(status, world_communicator);
}

/* Local: the allowed start sectors and their records (L + 1 slots).  The
 * caller agrees the status on the world before any draw. */
static MVMCKrylovStatus gc_prepare_initial(MVMCPowerLanczosGCRun *gc) {
  const size_t slots = gc->model->site_count + 1;
  MVMCKrylovStatus status;
  if (slots == 0 || slots > SIZE_MAX / MVMC_POWER_LANCZOS_INITIAL_ATTEMPTS ||
      slots > SIZE_MAX / sizeof(*gc->initial_records)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  gc->initial_sectors = (int *)calloc(slots, sizeof(*gc->initial_sectors));
  gc->initial_records = (MVMCPowerLanczosGCInitialRecord *)calloc(
      slots, sizeof(*gc->initial_records));
  if (gc->initial_sectors == NULL || gc->initial_records == NULL) {
    return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
  status = mvmc_gc_krylov_initial_sectors(
      gc->model, &gc->spin_support, gc->initial_electron_count,
      gc->initial_electron_count_explicit, gc->initial_sectors, slots,
      &gc->initial_sector_count);
  return status;
}

#if defined(MVMC_ENABLE_FAULT_INJECTION)
static void gc_initial_fault(const MVMCPowerLanczosGCRun *gc,
                             uint64_t attempt, MVMCKrylovStatus *status,
                             MVMCKrylovBootstrapReason *reason) {
  const int chain = gc_testing_initial_chain == UINT64_MAX ||
                    gc_testing_initial_chain == gc->chain_index;
  const int recoverable = *status == MVMC_KRYLOV_STATUS_OK ||
                          (*status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
                           *reason != MVMC_KRYLOV_BOOTSTRAP_NO_RETRY);
  switch (gc_testing_initial_fault) {
    case MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_EXACT_ZERO:
    case MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_NUMERIC_ZERO:
      if (chain && recoverable && attempt < gc_testing_initial_count) {
        *status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
        *reason = gc_testing_initial_fault ==
                          MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_EXACT_ZERO
                      ? MVMC_KRYLOV_BOOTSTRAP_EXACT_ZERO
                      : MVMC_KRYLOV_BOOTSTRAP_NUMERIC_ZERO;
      }
      break;
    case MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_NONFINITE:
      if (chain && attempt == gc_testing_initial_count) {
        *status = MVMC_KRYLOV_STATUS_NONFINITE;
        *reason = MVMC_KRYLOV_BOOTSTRAP_NO_RETRY;
      }
      break;
    case MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_STATUS:
      if (gc->chain_rank == 1 && attempt == 0) {
        *status = MVMC_KRYLOV_STATUS_NONFINITE;
        *reason = MVMC_KRYLOV_BOOTSTRAP_NO_RETRY;
      }
      break;
    case MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_NONE:
      break;
  }
}
#endif

/*
 * Finite-S coefficient start (chain collectives only).  The allowed sectors
 * are tried in order with at most 100 structured draws each from the
 * bootstrap RNG; the root draws into the candidate buffer and broadcasts
 * it.  A finite nonzero base Psi enters the positive-sampler initialize and
 * only a success becomes the current configuration.  Exact and numeric
 * zeros of Psi are redrawn (a guide floor is not a nonzero proof); every
 * other status ends the start after the chain agreed on status and reason.
 * Chains may select different N and use different retry counts; the world
 * meets at the synchronization before the pilot.
 */
static MVMCKrylovStatus gc_initialize_coefficient(
    MVMCPowerLanczosGCRun *gc, MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *guide, uint64_t *words,
    uint64_t *candidate, size_t word_count,
    MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovPositiveSamplerSnapshot *snapshot) {
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
#if defined(MVMC_ENABLE_FAULT_INJECTION)
  uint64_t attempt_index = 0;
#endif
  size_t sector;
  int selected = 0, stopped = 0;
  gc->initial_tried = 0;
  for (sector = 0; !selected && !stopped && sector < gc->initial_sector_count;
       ++sector) {
    MVMCPowerLanczosGCInitialRecord *record = gc->initial_records + sector;
    const int n = gc->initial_sectors[sector];
    size_t attempt;
    gc->initial_tried = sector + 1;
    for (attempt = 0; attempt < MVMC_POWER_LANCZOS_INITIAL_ATTEMPTS;
         ++attempt) {
      MVMCKrylovScaledAmplitudeResult psi;
      MVMCKrylovBootstrapReason reason = MVMC_KRYLOV_BOOTSTRAP_NO_RETRY;
      status = MVMC_KRYLOV_STATUS_OK;
      if (gc->chain_rank == 0) {
        status = mvmc_gc_krylov_draw_spin_initial(
            gc->model, &gc->spin_support, n, draw_bounded, rng, candidate,
            word_count);
#if defined(MVMC_ENABLE_FAULT_INJECTION)
        ++gc_testing_counters.generated[gc->stage];
#endif
      }
      status = gc_broadcast(gc, status, candidate, word_count, NULL);
      if (status != MVMC_KRYLOV_STATUS_OK) break;
      ++record->attempts;
      memset(&psi, 0, sizeof(psi));
      status = amplitude(candidate, word_count, amplitude_context, &psi);
      if (status == MVMC_KRYLOV_STATUS_OK) {
        switch (psi.value.state) {
          case MVMC_SCALED_COMPLEX_FINITE_NONZERO:
            status = mvmc_krylov_positive_sampler_initialize(
                workspace, guide, candidate, word_count, amplitude,
                amplitude_context, snapshot);
            break;
          case MVMC_SCALED_COMPLEX_EXACT_ZERO:
            reason = MVMC_KRYLOV_BOOTSTRAP_EXACT_ZERO;
            status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
            break;
          case MVMC_SCALED_COMPLEX_NUMERIC_ZERO:
            reason = MVMC_KRYLOV_BOOTSTRAP_NUMERIC_ZERO;
            status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
            break;
          default:
            status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
            break;
        }
      }
#if defined(MVMC_ENABLE_FAULT_INJECTION)
      gc_initial_fault(gc, attempt_index++, &status, &reason);
#endif
      status = gc_agree(gc, status, reason);
      if (status == MVMC_KRYLOV_STATUS_OK) {
        memcpy(words, candidate, word_count * sizeof(*words));
        record->selected = 1;
        gc->initial_selected_n = n;
        selected = 1;
        break;
      }
      if (status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
          reason != MVMC_KRYLOV_BOOTSTRAP_NO_RETRY) {
        if (reason == MVMC_KRYLOV_BOOTSTRAP_EXACT_ZERO) {
          ++record->exact_zero;
        } else {
          ++record->numeric_zero;
        }
#if defined(MVMC_ENABLE_FAULT_INJECTION)
        ++gc_testing_counters.initial_retries;
#endif
        continue;
      }
      break;
    }
    record->status = status;
    stopped = !selected && attempt < MVMC_POWER_LANCZOS_INITIAL_ATTEMPTS;
  }
  if (!selected && !stopped) {
    gc->initial_exhausted = 1;
    status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  if (gc->chain_rank == 0) {
    for (sector = 0; sector < gc->initial_tried; ++sector) {
      const MVMCPowerLanczosGCInitialRecord *record =
          gc->initial_records + sector;
      fprintf(stdout,
              "GC Lanczos initialization: chain=%llu N=%d attempts=%llu "
              "exact_zero=%llu numeric_zero=%llu selected=%d status=%s\n",
              (unsigned long long)gc->chain_index, gc->initial_sectors[sector],
              (unsigned long long)record->attempts,
              (unsigned long long)record->exact_zero,
              (unsigned long long)record->numeric_zero, record->selected,
              mvmc_krylov_status_string(record->status));
    }
  }
  if (!selected) gc_record_failure(gc, status, candidate, word_count, NULL);
  return selected ? MVMC_KRYLOV_STATUS_OK : status;
}

/*
 * Final-chain start: attempt a uses an initial-N configuration for even a
 * and a whole-domain configuration for odd a.  Only EXACT_ZERO,
 * NUMERIC_ZERO and an unavailable local energy are redrawn; everything else
 * ends the run.  The transition kernel and the target are not changed.
 */
static MVMCKrylovStatus gc_final_bootstrap(
    MVMCPowerLanczosGCRun *gc, MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy, uint64_t *words,
    size_t word_count, MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovFinalStateSnapshot *snapshot) {
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  MVMCKrylovFinalDiagnostic diagnostic;
  size_t attempt;
  gc_b0_stage(gc, MVMC_POWER_LANCZOS_GC_STAGE_FINAL_BOOTSTRAP);
  for (attempt = 0; attempt < MVMC_POWER_LANCZOS_BOOTSTRAP_ATTEMPTS;
       ++attempt) {
    MVMCKrylovBootstrapReason reason = MVMC_KRYLOV_BOOTSTRAP_NO_RETRY;
    /* A failed draw must not report the previous candidate. */
    memset(&diagnostic, 0, sizeof(diagnostic));
    status = gc_draw_start(gc, rng, words, word_count, attempt % 2 == 0);
    if (status != MVMC_KRYLOV_STATUS_OK) break;
    ++gc->bootstrap_attempts;
#if defined(MVMC_ENABLE_FAULT_INJECTION)
    ++gc_testing_counters.steps[gc->stage];
    if (gc->chain_rank == 0 && attempt % 2 == 0) {
      ++gc_testing_counters.bootstrap_initial_draws;
      if (gc_finite_spin(gc)) {
        MVMCGCSpinSupportResult support;
        int occupation[2 * 64];
        const size_t L = gc->model->site_count;
        size_t orbital;
        int off = 2 * L > sizeof(occupation) / sizeof(occupation[0]);
        for (orbital = 0; !off && orbital < 2 * L; ++orbital) {
          occupation[orbital] =
              (int)((words[orbital / 64] >> (orbital % 64)) & UINT64_C(1));
        }
        off = off ||
              mvmc_gc_spin_support_evaluate(
                  &gc->spin_support, (int)L, gc->model->domain, occupation,
                  2 * L, &support) != MVMC_KRYLOV_STATUS_OK ||
              support.electron_count != gc->initial_selected_n ||
              support.single_count < 2 * gc->spin_support.total_spin;
        gc_testing_counters.bootstrap_initial_off_sector += off != 0;
      }
    } else if (gc->chain_rank == 0) {
      ++gc_testing_counters.bootstrap_uniform_draws;
    }
#endif
    status = mvmc_krylov_final_state_sampler_initialize_diagnosed(
        workspace, policy, words, word_count, amplitude, amplitude_context,
        snapshot, &diagnostic);
#if defined(MVMC_ENABLE_FAULT_INJECTION)
    if (attempt == 0 &&
        (gc_testing_fault ==
             MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_LOCAL_ENERGY ||
         gc_testing_fault ==
             MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_NONFINITE) &&
        (status == MVMC_KRYLOV_STATUS_OK ||
         status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE)) {
      memset(snapshot, 0, sizeof(*snapshot));
      if (gc_testing_fault ==
          MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_LOCAL_ENERGY) {
        diagnostic.evaluation_available = 1;
        diagnostic.retry_reason = MVMC_KRYLOV_BOOTSTRAP_LOCAL_ENERGY;
        status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
      } else {
        diagnostic.retry_reason = MVMC_KRYLOV_BOOTSTRAP_NO_RETRY;
        status = MVMC_KRYLOV_STATUS_NONFINITE;
      }
    }
    if ((gc_testing_fault ==
             MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_NUMERIC &&
         attempt == 0) ||
        (gc_testing_fault == MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_LAST &&
         attempt + 1 < MVMC_POWER_LANCZOS_BOOTSTRAP_ATTEMPTS) ||
        gc_testing_fault ==
            MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_EXHAUSTED) {
      if (status == MVMC_KRYLOV_STATUS_OK ||
          status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE) {
        memset(snapshot, 0, sizeof(*snapshot));
        mvmc_scaled_complex_make_numeric_zero(
            -30.0, -INFINITY, -INFINITY, 0.0,
            &diagnostic.evaluation.amplitude);
        diagnostic.evaluation_available = 1;
        diagnostic.retry_reason = MVMC_KRYLOV_BOOTSTRAP_NUMERIC_ZERO;
        status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
      }
    }
#endif
    if (status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE) {
      reason = diagnostic.retry_reason;
    }
    status = gc_agree(gc, status, reason);
    if (status == MVMC_KRYLOV_STATUS_OK) {
#if defined(MVMC_ENABLE_FAULT_INJECTION)
      gc_testing_counters.final_start_word = words[0];
#endif
      break;
    }
    if (status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
        reason != MVMC_KRYLOV_BOOTSTRAP_NO_RETRY) {
      if (reason == MVMC_KRYLOV_BOOTSTRAP_EXACT_ZERO) {
        ++gc->bootstrap_exact_zero;
      } else if (reason == MVMC_KRYLOV_BOOTSTRAP_NUMERIC_ZERO) {
        ++gc->bootstrap_numeric_zero;
      } else {
        ++gc->bootstrap_local_energy;
      }
      GC_COUNT(bootstrap_retries);
      continue;
    }
    break;
  }
  if (status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
      attempt == MVMC_POWER_LANCZOS_BOOTSTRAP_ATTEMPTS) {
    gc->bootstrap_exhausted = 1;
  }
  gc_record_failure(gc, status, words, word_count, &diagnostic);
  return status;
}

static void log_sum_add(double value, double *maximum, double *scaled_sum) {
  if (value > *maximum) {
    *scaled_sum =
        *maximum == -INFINITY ? 1.0
                             : *scaled_sum * exp(*maximum - value) + 1.0;
    *maximum = value;
  } else {
    *scaled_sum += exp(value - *maximum);
  }
}

static uint64_t guide_policy_hash(
    int order, double eta, const double *log_basis_scale) {
  uint64_t hash = UINT64_C(1469598103934665603);
  uint64_t bits = 0;
  int index;
  hash = hash_u64(hash, UINT64_C(0x504c534c494d4731));
  hash = hash_u64(hash, (uint64_t)(unsigned)order);
  memcpy(&bits, &eta, sizeof(bits));
  hash = hash_u64(hash, bits);
  for (index = 0; index <= order; ++index) {
    memcpy(&bits, log_basis_scale + index, sizeof(bits));
    hash = hash_u64(hash, bits);
  }
  return hash == 0 ? UINT64_C(1) : hash;
}

static MVMCKrylovStatus scale_pilot(
    MVMCPowerLanczosGCRun *gc, MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    size_t pilot_warm_up, size_t pilot_sample_count,
    uint64_t *words, size_t word_count,
    MVMCKrylovPositiveSamplerSnapshot *snapshot,
    MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words, MVMCKrylovPositiveGuidePolicy *guide,
    MVMCKrylovMatrixMeasurementPolicy *matrix,
    const MVMCPowerLanczosSamplingTopology *sampling) {
  double local_maximum[MVMC_KRYLOV_MAX_ORDER + 1];
  double global_maximum[MVMC_KRYLOV_MAX_ORDER + 1];
  double local_scaled_sum[MVMC_KRYLOV_MAX_ORDER + 1];
  double global_scaled_sum[MVMC_KRYLOV_MAX_ORDER + 1];
  size_t step;
  size_t sample;
  int index;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  for (index = 0; index <= guide->order; ++index) {
    local_maximum[index] = -INFINITY;
    local_scaled_sum[index] = 0.0;
  }
  if (gc != NULL) gc_b0_stage(gc, MVMC_POWER_LANCZOS_GC_STAGE_PILOT_WARM_UP);
  for (step = 0; status == MVMC_KRYLOV_STATUS_OK &&
                 step < pilot_warm_up; ++step) {
    status = coefficient_step_any(
        gc, workspace, guide, model, proposal_policy, words, word_count,
        snapshot, rng, amplitude, amplitude_context, proposal_words);
  }
  if (gc != NULL) gc_b0_stage(gc, MVMC_POWER_LANCZOS_GC_STAGE_PILOT_SAMPLE);
  for (sample = 0; status == MVMC_KRYLOV_STATUS_OK &&
                   sample < pilot_sample_count; ++sample) {
    status = coefficient_step_any(
        gc, workspace, guide, model, proposal_policy, words, word_count,
        snapshot, rng, amplitude, amplitude_context, proposal_words);
    if (status == MVMC_KRYLOV_STATUS_OK) gc_b0_sample(gc, sample, words, word_count);
    for (index = 0; status == MVMC_KRYLOV_STATUS_OK &&
                    index <= guide->order; ++index) {
      const MVMCScaledComplex *value = snapshot->krylov.value + index;
      if (!mvmc_scaled_complex_is_valid(value) ||
          value->state == MVMC_SCALED_COMPLEX_NONFINITE) {
        status = MVMC_KRYLOV_STATUS_NONFINITE;
      } else if (value->state == MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
        const double log_squared = 2.0 * value->log_abs;
        if (!isfinite(log_squared)) {
          status = MVMC_KRYLOV_STATUS_NONFINITE;
        } else {
          log_sum_add(log_squared, local_maximum + index,
                      local_scaled_sum + index);
        }
      }
    }
  }
  status = mvmc_power_lanczos_sampling_synchronize(sampling, status);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    memcpy(global_maximum, local_maximum,
           sizeof(global_maximum));
    status = mvmc_power_lanczos_sampling_max_double(
        sampling, global_maximum, (size_t)guide->order + 1);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    for (index = 0; index <= guide->order; ++index) {
      global_scaled_sum[index] =
          local_maximum[index] == -INFINITY
              ? 0.0
              : local_scaled_sum[index] *
                    exp(local_maximum[index] - global_maximum[index]);
    }
    status = mvmc_power_lanczos_sampling_sum_double(
        sampling, global_scaled_sum, (size_t)guide->order + 1);
  }
  for (index = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  index <= guide->order; ++index) {
    double scale = 0.0;
    if (global_maximum[index] != -INFINITY) {
      const double pilot_total =
          (double)pilot_sample_count * (double)sampling->chain_count;
      const double log_mean =
          global_maximum[index] + log(global_scaled_sum[index]) -
          log(pilot_total);
      scale = -0.5 * log_mean;
      if (!(pilot_total > 0.0) || !isfinite(pilot_total) ||
          !(global_scaled_sum[index] > 0.0) ||
          !isfinite(log_mean) || !isfinite(scale)) {
        status = MVMC_KRYLOV_STATUS_NONFINITE;
      }
    }
    guide->lambda[index] = 1.0;
    guide->log_basis_scale[index] = scale;
    matrix->guide_lambda[index] = 1.0;
    matrix->target_weight[index] = 1.0;
    matrix->log_basis_scale[index] = scale;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    const double eta = 0x1p-40 * (double)(guide->order + 1);
    guide->eta = eta;
    matrix->eta = eta;
    guide->policy_hash =
        guide_policy_hash(guide->order, eta, guide->log_basis_scale);
    status = mvmc_krylov_positive_sampler_initialize(
        workspace, guide, words, word_count, amplitude, amplitude_context,
        snapshot);
  }
  return status;
}

static MVMCKrylovStatus initialize_matrix_blocks(
    size_t block_count, int dimension, size_t upper_count,
    MVMCKrylovMatrixMeasurementAccumulator *accumulators,
    MVMCKrylovStreamingComplexSum overlap[][6],
    MVMCKrylovStreamingComplexSum hamiltonian[][6],
    MVMCKrylovStreamingComplexSum hamiltonian_adjoint[][6],
    MVMCKrylovStreamingComplexSum hamiltonian_squared[][6]) {
  size_t block;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  for (block = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  block < block_count; ++block) {
    status =
        mvmc_krylov_matrix_measurement_accumulator_init_with_diagnostics(
            dimension, overlap[block], hamiltonian[block],
            hamiltonian_adjoint[block], hamiltonian_squared[block],
            upper_count, accumulators + block);
  }
  return status;
}

static MVMCKrylovStatus extract_matrix_block(
    const MVMCKrylovMatrixMeasurementAccumulator *accumulator,
    MVMCPowerLanczosMatrixBlock *block) {
  size_t row;
  size_t column;
  size_t entry = 0;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  if (block == NULL || accumulator == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(block, 0, sizeof(*block));
  block->count = accumulator->sample_count;
  for (row = 0; status == MVMC_KRYLOV_STATUS_OK &&
                row < accumulator->dimension; ++row) {
    for (column = row; status == MVMC_KRYLOV_STATUS_OK &&
                         column < accumulator->dimension;
         ++column, ++entry) {
      MVMCKrylovJackknifeBlock extracted;
      status = mvmc_krylov_matrix_measurement_extract_block(
          accumulator, MVMC_KRYLOV_MATRIX_OVERLAP, row, column,
          &extracted);
      if (status == MVMC_KRYLOV_STATUS_OK) {
        block->overlap[entry] = extracted.numerator;
        status = mvmc_krylov_matrix_measurement_extract_block(
            accumulator, MVMC_KRYLOV_MATRIX_HAMILTONIAN, row, column,
            &extracted);
      }
      if (status == MVMC_KRYLOV_STATUS_OK) {
        block->hamiltonian[entry] = extracted.numerator;
        status = mvmc_krylov_matrix_measurement_extract_block(
            accumulator, MVMC_KRYLOV_MATRIX_HAMILTONIAN_ADJOINT,
            row, column, &extracted);
      }
      if (status == MVMC_KRYLOV_STATUS_OK) {
        block->hamiltonian_adjoint[entry] = extracted.numerator;
        status = mvmc_krylov_matrix_measurement_extract_block(
            accumulator, MVMC_KRYLOV_MATRIX_HAMILTONIAN_SQUARED,
            row, column, &extracted);
      }
      if (status == MVMC_KRYLOV_STATUS_OK) {
        block->hamiltonian_squared[entry] = extracted.numerator;
      }
    }
  }
  return status;
}

static void matrix_block_add(
    MVMCPowerLanczosMatrixBlock *total,
    const MVMCPowerLanczosMatrixBlock *block, size_t upper_count) {
  size_t entry;
  total->count += block->count;
  for (entry = 0; entry < upper_count; ++entry) {
    total->overlap[entry] += block->overlap[entry];
    total->hamiltonian[entry] += block->hamiltonian[entry];
    total->hamiltonian_adjoint[entry] +=
        block->hamiltonian_adjoint[entry];
    total->hamiltonian_squared[entry] +=
        block->hamiltonian_squared[entry];
  }
}

static MVMCKrylovStatus reduce_matrix_blocks(
    const MVMCPowerLanczosSamplingTopology *sampling,
    MVMCPowerLanczosMatrixBlock *blocks, size_t block_count,
    size_t upper_count) {
  double complex values[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS * 4 *
      MVMC_POWER_LANCZOS_PACKED_CAPACITY];
  uint64_t counts[MVMC_POWER_LANCZOS_STABILIZED_BLOCKS];
  size_t cursor = 0;
  size_t block;
  size_t entry;
  MVMCKrylovStatus status;
  if (sampling == NULL || blocks == NULL || block_count == 0 ||
      block_count > MVMC_POWER_LANCZOS_STABILIZED_BLOCKS ||
      upper_count > MVMC_POWER_LANCZOS_PACKED_CAPACITY) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(values, 0, sizeof(values));
  memset(counts, 0, sizeof(counts));
  for (block = 0; block < block_count; ++block) {
    counts[block] = blocks[block].count;
    for (entry = 0; entry < upper_count; ++entry) {
      values[cursor++] = blocks[block].overlap[entry];
      values[cursor++] = blocks[block].hamiltonian[entry];
      values[cursor++] = blocks[block].hamiltonian_adjoint[entry];
      values[cursor++] = blocks[block].hamiltonian_squared[entry];
    }
  }
  status = mvmc_power_lanczos_sampling_sum_complex(
      sampling, values, cursor);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_sampling_sum_u64(
        sampling, counts, block_count);
  }
  cursor = 0;
  for (block = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  block < block_count; ++block) {
    blocks[block].count = counts[block];
    for (entry = 0; entry < upper_count; ++entry) {
      blocks[block].overlap[entry] = values[cursor++];
      blocks[block].hamiltonian[entry] = values[cursor++];
      blocks[block].hamiltonian_adjoint[entry] = values[cursor++];
      blocks[block].hamiltonian_squared[entry] = values[cursor++];
    }
  }
  return status;
}

static MVMCPowerLanczosMatrixBlock matrix_block_without(
    const MVMCPowerLanczosMatrixBlock *total,
    const MVMCPowerLanczosMatrixBlock *removed, size_t upper_count) {
  MVMCPowerLanczosMatrixBlock result = *total;
  size_t entry;
  result.count -= removed->count;
  for (entry = 0; entry < upper_count; ++entry) {
    result.overlap[entry] -= removed->overlap[entry];
    result.hamiltonian[entry] -= removed->hamiltonian[entry];
    result.hamiltonian_adjoint[entry] -=
        removed->hamiltonian_adjoint[entry];
    result.hamiltonian_squared[entry] -=
        removed->hamiltonian_squared[entry];
  }
  return result;
}

static double complex packed_quadratic(
    const double complex *packed, int dimension,
    const double complex *coefficient) {
  double complex value = 0.0;
  size_t entry = 0;
  int row;
  int column;
  for (row = 0; row < dimension; ++row) {
    for (column = row; column < dimension; ++column, ++entry) {
      value += conj(coefficient[row]) * packed[entry] *
               coefficient[column];
      if (row != column) {
        value += conj(coefficient[column]) * conj(packed[entry]) *
                 coefficient[row];
      }
    }
  }
  return value;
}

static MVMCKrylovStatus solve_gevp(
    const MVMCKrylovGEVPPolicy *policy,
    const MVMCPowerLanczosMatrixBlock *matrix, int dimension,
    size_t upper_count, MVMCKrylovGEVPResult *result) {
  MVMCKrylovGEVPStatus status =
      mvmc_krylov_gevp_solve_complex_packed(
          policy, dimension, matrix->overlap, matrix->hamiltonian,
          matrix->hamiltonian_adjoint, matrix->hamiltonian_squared,
          upper_count, result);
  if (status != MVMC_KRYLOV_GEVP_OK) {
    fprintf(stderr, "Lanczos GEVP failed: %s (count=%llu, residual=%.17g, "
            "variance=%.17g).\n", mvmc_krylov_gevp_status_string(status),
            (unsigned long long)matrix->count,
            result->gevp_relative_residual, result->variance);

  }
  return map_gevp_status(status);
}

static MVMCKrylovStatus physical_second_moment(
    const MVMCPowerLanczosMatrixBlock *matrix,
    const MVMCKrylovGEVPResult *gevp, int dimension,
    double complex *second_moment) {
  const double complex normalization =
      packed_quadratic(matrix->overlap, dimension, gevp->coefficient);
  const double complex numerator =
      packed_quadratic(
          matrix->hamiltonian_squared, dimension, gevp->coefficient);
  if (second_moment == NULL ||
      !isfinite(creal(normalization)) ||
      !isfinite(cimag(normalization)) ||
      !(creal(normalization) > 0.0) ||
      fabs(cimag(normalization)) >
          1.0e-10 * fmax(1.0, fabs(creal(normalization))) ||
      !isfinite(creal(numerator)) || !isfinite(cimag(numerator))) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  *second_moment = numerator / normalization;
  return isfinite(creal(*second_moment)) &&
                 isfinite(cimag(*second_moment))
             ? MVMC_KRYLOV_STATUS_OK
             : MVMC_KRYLOV_STATUS_NONFINITE;
}

static double jackknife_standard_error(
    const double *leave_one, size_t block_count) {
  double mean = 0.0;
  double sum = 0.0;
  size_t block;
  for (block = 0; block < block_count; ++block) mean += leave_one[block];
  mean /= (double)block_count;
  for (block = 0; block < block_count; ++block) {
    const double delta = leave_one[block] - mean;
    sum += delta * delta;
  }
  return sqrt(((double)block_count - 1.0) /
              (double)block_count * sum);
}

static double block_mean_standard_error(
    const double complex *block_sum, const uint64_t *block_sample_count,
    size_t block_count, double overall_mean) {
  double squared_deviation = 0.0;
  size_t block;
  for (block = 0; block < block_count; ++block) {
    const double block_mean =
        creal(block_sum[block]) / (double)block_sample_count[block];
    const double delta = block_mean - overall_mean;
    squared_deviation += delta * delta;
  }
  return sqrt(squared_deviation /
              ((double)block_count * ((double)block_count - 1.0)));
}

/*
 * Shared stabilized core.  gc == NULL is the canonical path with its
 * historical bridge, proposals, RNG order and failure statuses; otherwise
 * the GC bridge, initial configurations, proposals and bootstrap are used
 * and every statistical step is the same code.
 */
static MVMCKrylovStatus stabilized_core(
    const MVMCPowerLanczosStabilizedInput *input,
    const MVMCPowerLanczosGCView *gc_view, MVMCPowerLanczosGCRun *gc,
    MVMCPowerLanczosStabilizedResult *result) {
  MVMCPowerLanczosClassicBridge *bridge = NULL;
  MVMCPowerLanczosGCBridge *gc_bridge = NULL;
  const MVMCKrylovFockModel *model = NULL;
  MVMCKrylovScaledAmplitudeCallback amplitude = NULL;
  void *amplitude_context = NULL;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  MVMCKrylovPositiveGuidePolicy guide;
  MVMCKrylovMatrixMeasurementPolicy matrix_policy;
  MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
  MVMCKrylovPositiveSamplerRng bootstrap_rng;
  MVMCKrylovPositiveSamplerRng coefficient_rng;
  MVMCKrylovPositiveSamplerRng final_rng;
  MVMCKrylovPositiveSamplerSnapshot coefficient_snapshot;
  MVMCKrylovFinalStateSnapshot final_snapshot;
  MVMCKrylovFinalStatePolicy final_policy;
  MVMCKrylovGEVPPolicy gevp_policy;
  MVMCKrylovGEVPResult gevp;
  MVMCPowerLanczosSamplingTopology sampling;
  MVMCKrylovMatrixMeasurementAccumulator accumulators[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS];
  MVMCKrylovStreamingComplexSum overlap_entries[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS][6];
  MVMCKrylovStreamingComplexSum hamiltonian_entries[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS][6];
  MVMCKrylovStreamingComplexSum hamiltonian_adjoint_entries[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS][6];
  MVMCKrylovStreamingComplexSum hamiltonian_squared_entries[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS][6];
  MVMCPowerLanczosMatrixBlock matrix_blocks[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS];
  MVMCPowerLanczosMatrixBlock matrix_total;
  MVMCPowerLanczosGCNumberResult number = {0};
  double complex final_number_block_sum[MVMC_POWER_LANCZOS_STABILIZED_BLOCKS] = {0};
  double complex final_number_square_block_sum[MVMC_POWER_LANCZOS_STABILIZED_BLOCKS] = {0};
  double complex final_energy_block_sum[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS];
  uint64_t final_energy_block_count[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS];
  double leave_one_second_moment[
      MVMC_POWER_LANCZOS_STABILIZED_BLOCKS];
  double *energy_series = NULL;
  uint64_t *words = NULL;
  uint64_t *proposal_words = NULL;
  size_t dimension = 0;
  size_t upper_count = 0;
  size_t word_count = 0;
  size_t block_count = 0;
  size_t pilot_warm_up = 0;
  size_t pilot_sample_count = 0;
  uint64_t generation_hash = 0;
  uint64_t provenance_hash = 0;
  double complex energy_sum = 0.0;
  double complex energy = NAN + I * NAN;
  double complex second_moment = NAN + I * NAN;
  double complex variance = NAN + I * NAN;
  double energy_standard_error = NAN;
  double second_moment_standard_error = NAN;
  double effective_sample_count = NAN;
  double complex coefficient_energy = NAN + I * NAN;
  int coefficient_energy_available = 0;
  MVMCKrylovTauIntResult tau;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  size_t step;
  size_t sample;
  size_t block;
  int index;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_result(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  memset(&bootstrap_rng, 0, sizeof(bootstrap_rng));
  memset(&coefficient_rng, 0, sizeof(coefficient_rng));
  memset(&final_rng, 0, sizeof(final_rng));
  if (input == NULL ||
      (gc == NULL ? input->classic_view == NULL
                  : (gc_view == NULL || input->power_step != 1)) ||
      (input->power_step != 1 && input->power_step != 2) ||
      input->seed == 0 || input->sample_count < 8 || input->interval == 0 ||
      input->sample_count > SIZE_MAX / sizeof(*energy_series)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  dimension = (size_t)input->power_step + 1;
  status = mvmc_power_lanczos_sampling_topology_create(
      input->world_communicator, input->chain_communicator,
      input->sample_count, &sampling);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_result(status, result);
    return status;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_streaming_upper_count(dimension, &upper_count);
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      upper_count > MVMC_POWER_LANCZOS_PACKED_CAPACITY) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  block_count = mvmc_power_lanczos_stabilized_block_count(
      (size_t)sampling.total_samples);
  if (block_count == 0) status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;

  if (status == MVMC_KRYLOV_STATUS_OK && gc == NULL) {
    status = mvmc_power_lanczos_classic_bridge_create(
        input->classic_view, input->world_communicator,
        input->chain_communicator, &bridge);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      model = mvmc_power_lanczos_classic_bridge_model(bridge);
      amplitude = mvmc_power_lanczos_classic_bridge_amplitude(bridge);
      amplitude_context =
          mvmc_power_lanczos_classic_bridge_amplitude_context(bridge);
      generation_hash =
          mvmc_power_lanczos_classic_bridge_generation_hash(bridge);
    }
  } else if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_gc_bridge_create(
        gc_view, input->world_communicator, input->chain_communicator,
        &gc_bridge);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      model = mvmc_power_lanczos_gc_bridge_model(gc_bridge);
      amplitude = mvmc_power_lanczos_gc_bridge_amplitude(gc_bridge);
      amplitude_context =
          mvmc_power_lanczos_gc_bridge_amplitude_context(gc_bridge);
      generation_hash =
          mvmc_power_lanczos_gc_bridge_generation_hash(gc_bridge);
      gc->model = model;
      gc->spin_support = gc_view->spin_support;
      gc->chain_communicator = input->chain_communicator;
      gc->chain_rank = sampling.chain_rank;
      gc->chain_index = sampling.chain_index;
      gc->topology_ready = 1;
      /* World-agreed bridge status: every rank enters this collective. */
      status = gc_agree_inputs(gc, input->world_communicator);
      /* A forbidden explicit N is agreed with the synchronization below,
       * before any start draw. */
      if (status == MVMC_KRYLOV_STATUS_OK && gc_finite_spin(gc)) {
        status = gc_prepare_initial(gc);
      }
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (model == NULL || amplitude == NULL || amplitude_context == NULL ||
       generation_hash == 0)) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    word_count = mvmc_krylov_fock_word_count(model->site_count);
    if (word_count == 0 ||
        word_count > SIZE_MAX / sizeof(*words)) {
      status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    words = (uint64_t *)calloc(word_count, sizeof(*words));
    proposal_words =
        (uint64_t *)calloc(word_count, sizeof(*proposal_words));
    if (gc != NULL) gc->step_words = (uint64_t *)calloc(word_count, sizeof(*words));
    energy_series =
        (double *)calloc(input->sample_count, sizeof(*energy_series));
    if (words == NULL || proposal_words == NULL || energy_series == NULL ||
        (gc != NULL && gc->step_words == NULL)) {
      status = MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = create_limits(
        model, input->power_step, generation_hash, &limits);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_plan_create(model, &limits, &plan);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_workspace_create(plan, &workspace);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_proposal_policy_create(
        1, 16, &proposal_policy);
  }
  status = mvmc_power_lanczos_sampling_synchronize(&sampling, status);
  memset(&guide, 0, sizeof(guide));
  memset(&matrix_policy, 0, sizeof(matrix_policy));
  guide.order = (int)dimension;
  guide.eta = 0x1p-40;
  matrix_policy.order = (int)dimension;
  matrix_policy.eta = 0x1p-40;
  for (index = 0; index <= (int)dimension; ++index) {
    const double initial_weight = index == 0 ? 0x1p-40 : 1.0;
    guide.lambda[index] = initial_weight;
    matrix_policy.guide_lambda[index] = initial_weight;
    matrix_policy.target_weight[index] = 1.0;
  }
  guide.policy_hash =
      guide_policy_hash(guide.order, guide.eta, guide.log_basis_scale);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_rng_seed(
        domain_seed(input->seed,
                    chain_domain(UINT64_C(0x504c424f4f545354),
                                 sampling.chain_index)),
        UINT64_C(0x504c424f4f545354), &bootstrap_rng);
  }
  if (status == MVMC_KRYLOV_STATUS_OK && gc == NULL) {
    GC_COUNT(fixed_sector_calls);
    status = draw_uniform_sector(
        model, input->chain_communicator, &bootstrap_rng, words,
        word_count);
  } else if (status == MVMC_KRYLOV_STATUS_OK) {
    gc_b0_stage(gc, MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_INITIAL);
    if (!gc_finite_spin(gc)) {
      status = gc_draw_start(gc, &bootstrap_rng, words, word_count, 1);
      gc_record_failure(gc, status, words, word_count, NULL);
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_session_begin(
        workspace, amplitude, amplitude_context, generation_hash);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_rng_seed(
        domain_seed(input->seed,
                    chain_domain(UINT64_C(0x504c434f45464631),
                                 sampling.chain_index)),
        UINT64_C(0x504c434f45464631), &coefficient_rng);
  }
  if (status == MVMC_KRYLOV_STATUS_OK && gc_finite_spin(gc)) {
    /* Draws from the bootstrap RNG only; coefficient_rng is untouched. */
    status = gc_initialize_coefficient(
        gc, workspace, &guide, words, proposal_words, word_count,
        &bootstrap_rng, amplitude, amplitude_context, &coefficient_snapshot);
  } else if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_initialize(
        workspace, &guide, words, word_count, amplitude,
        amplitude_context, &coefficient_snapshot);
  }
#if defined(MVMC_ENABLE_FAULT_INJECTION)
  gc_testing_counters.bootstrap_rng_draws_after_start = bootstrap_rng.draws;
  gc_testing_counters.coefficient_rng_draws_at_pilot = coefficient_rng.draws;
#endif
  gc_record_failure(gc, status, words, word_count, NULL);
  /* The start is chain-local; the pilot below has world collectives, so
   * every chain must enter it or none. */
  status = mvmc_power_lanczos_sampling_synchronize(&sampling, status);
  pilot_warm_up = input->warm_up < 256 ? input->warm_up : 256;
  pilot_sample_count =
      input->sample_count < 256 ? input->sample_count : 256;
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = scale_pilot(
        gc, workspace, model, &proposal_policy, pilot_warm_up,
        pilot_sample_count, words, word_count, &coefficient_snapshot,
        &coefficient_rng, amplitude, amplitude_context, proposal_words,
        &guide, &matrix_policy, &sampling);
  }
  if (gc != NULL) gc_b0_stage(gc, MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_WARM_UP);
  for (step = 0; status == MVMC_KRYLOV_STATUS_OK &&
                 step < input->warm_up; ++step) {
    status = coefficient_step_any(
        gc, workspace, &guide, model, &proposal_policy, words, word_count,
        &coefficient_snapshot, &coefficient_rng, amplitude,
        amplitude_context, proposal_words);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = initialize_matrix_blocks(
        block_count, (int)dimension, upper_count, accumulators,
        overlap_entries, hamiltonian_entries,
        hamiltonian_adjoint_entries, hamiltonian_squared_entries);
  }
  status = mvmc_power_lanczos_sampling_synchronize(&sampling, status);
  if (gc != NULL) gc_b0_stage(gc, MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_SAMPLE);
  for (sample = 0; status == MVMC_KRYLOV_STATUS_OK &&
                   sample < input->sample_count; ++sample) {
    MVMCKrylovMatrixMeasurementSampleDiagnostics diagnostics;
    block = mvmc_power_lanczos_sampling_block_index(
        &sampling, sample, block_count);
    if (block >= block_count) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    for (step = 0; status == MVMC_KRYLOV_STATUS_OK &&
                   step < input->interval; ++step) {
      status = coefficient_step_any(
          gc, workspace, &guide, model, &proposal_policy, words, word_count,
          &coefficient_snapshot, &coefficient_rng, amplitude,
          amplitude_context, proposal_words);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      gc_b0_sample(gc, sample, words, word_count);
      status = mvmc_krylov_matrix_measurement_accumulator_add_sample(
          accumulators + block, &matrix_policy,
          coefficient_snapshot.krylov.value,
          (size_t)coefficient_snapshot.krylov.evaluated_order + 1,
          &diagnostics);
    }
  }
  if (workspace != NULL &&
      mvmc_bounded_krylov_session_is_active(workspace)) {
    const MVMCKrylovStatus end_status =
        mvmc_bounded_krylov_session_end(workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) status = end_status;
  }
  status = mvmc_power_lanczos_sampling_synchronize(&sampling, status);

  memset(&matrix_total, 0, sizeof(matrix_total));
  for (block = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  block < block_count; ++block) {
    status = extract_matrix_block(
        accumulators + block, matrix_blocks + block);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = reduce_matrix_blocks(
        &sampling, matrix_blocks, block_count, upper_count);
  }
  for (block = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  block < block_count; ++block) {
    if (matrix_blocks[block].count == 0) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    } else {
      matrix_block_add(
          &matrix_total, matrix_blocks + block, upper_count);
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      matrix_total.count != sampling.total_samples) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      mvmc_krylov_gevp_default_policy(0x1p-40, &gevp_policy) !=
          MVMC_KRYLOV_GEVP_OK) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = solve_gevp(
        &gevp_policy, &matrix_total, (int)dimension, upper_count, &gevp);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    /* E_coeff with the GEVP's symmetrized K and the same scaled basis. */
    double complex symmetric[MVMC_POWER_LANCZOS_PACKED_CAPACITY];
    double complex numerator, denominator;
    size_t entry;
    for (entry = 0; entry < upper_count; ++entry) {
      symmetric[entry] = 0.5 * (matrix_total.hamiltonian[entry] +
                                conj(matrix_total.hamiltonian_adjoint[entry]));
    }
    numerator = packed_quadratic(symmetric, (int)dimension, gevp.coefficient);
    denominator =
        packed_quadratic(matrix_total.overlap, (int)dimension,
                         gevp.coefficient);
    if (isfinite(creal(numerator)) && isfinite(cimag(numerator)) &&
        isfinite(creal(denominator)) && creal(denominator) > 0.0) {
      coefficient_energy = numerator / denominator;
      coefficient_energy_available =
          isfinite(creal(coefficient_energy)) &&
          isfinite(cimag(coefficient_energy));
    }
  }
  provenance_hash = hash_u64(generation_hash, input->seed);
  if (provenance_hash == 0) provenance_hash = UINT64_C(1);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_final_state_policy_create_scaled_basis(
        input->power_step,
        MVMC_KRYLOV_FINAL_COEFFICIENT_PRODUCTION_NOISY,
        provenance_hash, gevp.coefficient, matrix_policy.log_basis_scale,
        dimension, &final_policy);
  }
  if (status == MVMC_KRYLOV_STATUS_OK && gc != NULL) {
    gc->coefficient_available = 1;
    gc->coefficient[0] = final_policy.coefficient[0];
    gc->coefficient[1] = final_policy.coefficient[1];
    gc->coefficient_energy_available = coefficient_energy_available;
    gc->coefficient_energy = coefficient_energy;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_positive_sampler_rng_seed(
        domain_seed(input->seed,
                    chain_domain(UINT64_C(0x504c46494e414c31),
                                 sampling.chain_index)),
        UINT64_C(0x504c46494e414c31), &final_rng);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_session_begin(
        workspace, amplitude, amplitude_context, generation_hash);
  }
#if defined(MVMC_ENABLE_FAULT_INJECTION)
  if (status == MVMC_KRYLOV_STATUS_OK) {
    gc_testing_counters.final_rng_draws_at_bootstrap = final_rng.draws;
  }
#endif
  if (status == MVMC_KRYLOV_STATUS_OK && gc != NULL) {
    status = gc_final_bootstrap(gc, workspace, &final_policy, words,
                                word_count, &final_rng, amplitude,
                                amplitude_context, &final_snapshot);
  } else if (status == MVMC_KRYLOV_STATUS_OK) {
    size_t attempt;
    status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    for (attempt = 0;
         status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
         attempt < MVMC_POWER_LANCZOS_BOOTSTRAP_ATTEMPTS;
         ++attempt) {
      GC_COUNT(fixed_sector_calls);
      status = draw_uniform_sector(
          model, input->chain_communicator, &final_rng, words,
          word_count);
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = mvmc_krylov_final_state_sampler_initialize(
            workspace, &final_policy, words, word_count, amplitude,
            amplitude_context, &final_snapshot);
      }
    }
  }
  if (gc != NULL) gc_b0_stage(gc, MVMC_POWER_LANCZOS_GC_STAGE_FINAL_WARM_UP);
  for (step = 0; status == MVMC_KRYLOV_STATUS_OK &&
                 step < input->warm_up; ++step) {
    status = final_step_any(
        gc, workspace, &final_policy, model, &proposal_policy, words,
        word_count, &final_snapshot, &final_rng, amplitude,
        amplitude_context, proposal_words);
  }
  if (gc != NULL) gc_b0_stage(gc, MVMC_POWER_LANCZOS_GC_STAGE_FINAL_SAMPLE);
  memset(final_energy_block_sum, 0, sizeof(final_energy_block_sum));
  memset(final_energy_block_count, 0, sizeof(final_energy_block_count));
  for (sample = 0; status == MVMC_KRYLOV_STATUS_OK &&
                   sample < input->sample_count; ++sample) {
    double complex local_energy = NAN + I * NAN;
    MVMCScaledComplexExportStatus export_status;
    block = mvmc_power_lanczos_sampling_block_index(
        &sampling, sample, block_count);
    if (block >= block_count) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    for (step = 0; status == MVMC_KRYLOV_STATUS_OK &&
                   step < input->interval; ++step) {
      status = final_step_any(
          gc, workspace, &final_policy, model, &proposal_policy, words,
          word_count, &final_snapshot, &final_rng, amplitude,
          amplitude_context, proposal_words);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      export_status = mvmc_krylov_final_state_local_energy_export(
          &final_snapshot.final_state, &local_energy);
      if ((export_status != MVMC_SCALED_EXPORT_OK &&
           export_status != MVMC_SCALED_EXPORT_EXACT_ZERO) ||
          !isfinite(creal(local_energy)) ||
          !isfinite(cimag(local_energy))) {
        status = MVMC_KRYLOV_STATUS_NONFINITE;
      }
    }
    if (status == MVMC_KRYLOV_STATUS_OK && gc != NULL) {
      gc_b0_sample(gc, sample, words, word_count);
      const int count = popcount_words(words, word_count);
      const double n = (double)count;
      if (count < 0 || (size_t)count > 2 * model->site_count) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      } else {
        final_number_block_sum[block] += n;
        final_number_square_block_sum[block] += n * n;
        if (gc->sector_samples != NULL && (size_t)count < gc->sector_capacity)
          ++gc->sector_samples[count];
      }
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      energy_series[sample] = creal(local_energy);
      energy_sum += local_energy;
      final_energy_block_sum[block] += local_energy;
      ++final_energy_block_count[block];
    }
  }
  if (workspace != NULL &&
      mvmc_bounded_krylov_session_is_active(workspace)) {
    const MVMCKrylovStatus end_status =
        mvmc_bounded_krylov_session_end(workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) status = end_status;
  }
  status = mvmc_power_lanczos_sampling_synchronize(&sampling, status);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_sampling_sum_complex(
        &sampling, &energy_sum, 1);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_sampling_sum_complex(
        &sampling, final_energy_block_sum, block_count);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_sampling_sum_u64(
        &sampling, final_energy_block_count, block_count);
  }
  if (status == MVMC_KRYLOV_STATUS_OK && gc != NULL) {
    status = mvmc_power_lanczos_sampling_sum_complex(
        &sampling, final_number_block_sum, block_count);
  }
  if (status == MVMC_KRYLOV_STATUS_OK && gc != NULL) {
    status = mvmc_power_lanczos_sampling_sum_complex(
        &sampling, final_number_square_block_sum, block_count);
  }
  if (status == MVMC_KRYLOV_STATUS_OK && gc != NULL &&
      gc->sector_samples != NULL) {
    status = mvmc_power_lanczos_sampling_sum_u64(
        &sampling, gc->sector_samples, gc->sector_capacity);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    energy = energy_sum / (double)sampling.total_samples;
    status = physical_second_moment(
        &matrix_total, &gevp, (int)dimension, &second_moment);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    variance = second_moment - energy * energy;
    if (!isfinite(creal(variance)) || !isfinite(cimag(variance))) {
      status = MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  for (block = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  block < block_count; ++block) {
    const MVMCPowerLanczosMatrixBlock leave_matrix =
        matrix_block_without(&matrix_total, matrix_blocks + block,
                             upper_count);
    MVMCKrylovGEVPResult leave_gevp;
    double complex leave_second = NAN + I * NAN;
    status = solve_gevp(
        &gevp_policy, &leave_matrix, (int)dimension, upper_count,
        &leave_gevp);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = physical_second_moment(
          &leave_matrix, &leave_gevp, (int)dimension, &leave_second);
    }
    if (status == MVMC_KRYLOV_STATUS_OK &&
        !isfinite(creal(leave_second))) {
      status = MVMC_KRYLOV_STATUS_NONFINITE;
    }
    leave_one_second_moment[block] = creal(leave_second);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    energy_standard_error = block_mean_standard_error(
        final_energy_block_sum, final_energy_block_count, block_count,
        creal(energy));
    second_moment_standard_error = jackknife_standard_error(
        leave_one_second_moment, block_count);
    if (!isfinite(energy_standard_error) ||
        !isfinite(second_moment_standard_error)) {
      status = MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_tau_int_geyer_initial_positive(
        energy_series, input->sample_count, &tau);
  }
  status = mvmc_power_lanczos_sampling_synchronize(&sampling, status);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    /* The Geyer helper returns g = 1 + 2 sum(rho), not g/2. */
    effective_sample_count =
        (double)input->sample_count / tau.tau_int;
    status = mvmc_power_lanczos_sampling_sum_double(
        &sampling, &effective_sample_count, 1);
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (!(effective_sample_count > 0.0) ||
       !isfinite(effective_sample_count))) {
    status = MVMC_KRYLOV_STATUS_NONFINITE;
  }
  if (status == MVMC_KRYLOV_STATUS_OK && gc != NULL) {
    double sum = 0.0, square_sum = 0.0;
    double roundoff;
    for (block = 0; block < block_count; ++block) {
      number.block_sample_count[block] = final_energy_block_count[block];
      number.block_sum[block] = creal(final_number_block_sum[block]);
      number.block_square_sum[block] = creal(final_number_square_block_sum[block]);
      if (final_energy_block_count[block] != sampling.total_samples / block_count)
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      sum += number.block_sum[block];
      square_sum += number.block_square_sum[block];
    }
    number.mean = sum / (double)sampling.total_samples;
    number.mean_square = square_sum / (double)sampling.total_samples;
    number.variance = number.mean_square - number.mean * number.mean;
    roundoff = 64.0 * DBL_EPSILON *
        fmax(1.0, fmax(number.mean_square, number.mean * number.mean));
    if (number.variance < 0.0 && number.variance >= -roundoff)
      number.variance = 0.0;
    number.standard_error = block_mean_standard_error(
        final_number_block_sum, final_energy_block_count, block_count, number.mean);
    number.mean_square_standard_error = block_mean_standard_error(
        final_number_square_block_sum, final_energy_block_count, block_count,
        number.mean_square);
    if (!isfinite(number.mean) || !isfinite(number.mean_square) ||
        !isfinite(number.variance) || !isfinite(number.standard_error) ||
        !isfinite(number.mean_square_standard_error))
      status = MVMC_KRYLOV_STATUS_NONFINITE;
    else if (number.variance < 0.0)
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    number.valid = status == MVMC_KRYLOV_STATUS_OK;
  }
  if (gc != NULL)
    status = mvmc_power_lanczos_sampling_synchronize(&sampling, status);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    result->valid = 1;
    result->status = MVMC_KRYLOV_STATUS_OK;
    result->power_step = input->power_step;
    result->block_count = block_count;
    result->sampling_chains = sampling.chain_count;
    result->coefficient_samples = sampling.total_samples;
    result->final_samples = sampling.total_samples;
    result->retained_rank = gevp.retained_rank;
    result->condition_estimate = gevp.condition_estimate;
    result->gevp_residual = gevp.gevp_relative_residual;
    result->energy = creal(energy);
    result->energy_standard_error = energy_standard_error;
    result->variance = creal(variance);
    result->variance_standard_error = hypot(
        second_moment_standard_error,
        2.0 * cabs(energy) * energy_standard_error);
    result->final_energy_imaginary = cimag(energy);
    result->variance_imaginary = cimag(variance);
    result->energy_tau_int =
        (double)sampling.total_samples /
        (2.0 * effective_sample_count);
    result->effective_sample_count = effective_sample_count;
    result->coefficient_energy_available = coefficient_energy_available;
    result->coefficient_energy = coefficient_energy;
    result->gc_number = number;
  } else {
    invalidate_result(status, result);
  }

  free(energy_series);
  free(proposal_words);
  if (gc != NULL) {
    free(gc->step_words);
    gc->step_words = NULL;
  }
  free(words);
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
  mvmc_power_lanczos_classic_bridge_destroy(bridge);
  mvmc_power_lanczos_gc_bridge_destroy(gc_bridge);
  return status;
}

MVMCKrylovStatus mvmc_power_lanczos_stabilized_run(
    const MVMCPowerLanczosStabilizedInput *input,
    MVMCPowerLanczosStabilizedResult *result) {
  return stabilized_core(input, NULL, NULL, result);
}

enum {
  GC_CHAIN_FIELDS = 10,
  GC_FAILURE_FIELDS = 8,
  GC_REPORT_VALUES = MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS * GC_CHAIN_FIELDS +
                     GC_CHAIN_FIELDS + 1 + GC_FAILURE_FIELDS
};

/*
 * World aggregation without heap memory: every chain root writes its own
 * slot (others contribute zero) and the first failing chain (lowest index)
 * contributes the failure fields; two reductions in a fixed order.
 */
static void gc_collect_diagnostics(
    const MVMCPowerLanczosGCRun *gc,
    MVMCClassicPfaffianCommunicator world_communicator,
    MVMCPowerLanczosGCDiagnostics *diagnostics) {
  double values[GC_REPORT_VALUES];
  double reduced[GC_REPORT_VALUES];
  double *totals = values + MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS *
                                GC_CHAIN_FIELDS;
  double *failure = totals + GC_CHAIN_FIELDS + 1;
  uint64_t *sectors = diagnostics->final_sector_samples;
  const size_t capacity = diagnostics->final_sector_capacity;
  const int contributes = gc->topology_ready && gc->chain_rank == 0;
  double local_failure_chain =
      contributes && gc->failure_recorded ? (double)gc->chain_index
                                          : INFINITY;
  double first_failure = local_failure_chain;
  int rank = 0;
  size_t chain;
  memset(values, 0, sizeof(values));
#ifdef _mpi_use
  MPI_Comm_rank(world_communicator, &rank);
  if (MPI_Allreduce(&local_failure_chain, &first_failure, 1, MPI_DOUBLE,
                    MPI_MIN, world_communicator) != MPI_SUCCESS) {
    first_failure = INFINITY;
  }
#else
  (void)world_communicator;
#endif
  if (contributes) {
    double record[GC_CHAIN_FIELDS] = {
        (double)gc->bootstrap_attempts, (double)gc->bootstrap_exact_zero,
        (double)gc->bootstrap_numeric_zero,
        (double)gc->bootstrap_local_energy, (double)gc->bootstrap_exhausted,
        /* selected N + 1, so that an absent chain decodes as -1 */
        (double)gc->initial_selected_n + 1.0, 0.0, 0.0, 0.0,
        (double)gc->initial_exhausted};
    size_t field;
    for (field = 0; field < gc->initial_tried; ++field) {
      record[6] += (double)gc->initial_records[field].attempts;
      record[7] += (double)gc->initial_records[field].exact_zero;
      record[8] += (double)gc->initial_records[field].numeric_zero;
    }
    if (gc->chain_index < MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS) {
      memcpy(values + gc->chain_index * GC_CHAIN_FIELDS, record,
             sizeof(record));
    }
    for (field = 0; field < GC_CHAIN_FIELDS; ++field) {
      totals[field] = record[field];
    }
    totals[GC_CHAIN_FIELDS] = (double)(gc->chain_index + 1);
    if (gc->failure_recorded &&
        (double)gc->chain_index == first_failure) {
      failure[0] = 1.0;
      failure[1] = (double)gc->failure_stage;
      failure[2] = (double)gc->failure_status;
      failure[3] = (double)gc->failure_electron_count;
      failure[4] = (double)gc->failure_amplitude_state;
      failure[5] = gc->failure_log_abs;
      failure[6] = gc->failure_log_error_bound;
      failure[7] = (double)gc->failure_retry_reason;
    }
  }
  memcpy(reduced, values, sizeof(values));
#ifdef _mpi_use
  /* Chain count is a maximum; everything else is a disjoint sum. */
  {
    double count = totals[GC_CHAIN_FIELDS];
    double max_count = count;
    totals[GC_CHAIN_FIELDS] = 0.0;
    if (MPI_Reduce(&count, &max_count, 1, MPI_DOUBLE, MPI_MAX, 0,
                   world_communicator) != MPI_SUCCESS ||
        MPI_Reduce(values, reduced, GC_REPORT_VALUES, MPI_DOUBLE, MPI_SUM, 0,
                   world_communicator) != MPI_SUCCESS) {
      memset(reduced, 0, sizeof(reduced));
      max_count = 0.0;
    }
    reduced[MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS * GC_CHAIN_FIELDS +
            GC_CHAIN_FIELDS] = max_count;
  }
#endif
  memset(diagnostics, 0, sizeof(*diagnostics));
  diagnostics->final_sector_samples = sectors;
  diagnostics->final_sector_capacity = capacity;
  if (rank != 0) return;
  totals = reduced + MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS * GC_CHAIN_FIELDS;
  failure = totals + GC_CHAIN_FIELDS + 1;
  diagnostics->chain_count = (uint64_t)totals[GC_CHAIN_FIELDS];
  for (chain = 0; chain < MVMC_POWER_LANCZOS_GC_REPORTED_CHAINS; ++chain) {
    const double *record = reduced + chain * GC_CHAIN_FIELDS;
    diagnostics->bootstrap_attempts[chain] = (uint64_t)record[0];
    diagnostics->bootstrap_exact_zero[chain] = (uint64_t)record[1];
    diagnostics->bootstrap_numeric_zero[chain] = (uint64_t)record[2];
    diagnostics->bootstrap_local_energy[chain] = (uint64_t)record[3];
    diagnostics->bootstrap_exhausted[chain] = record[4] != 0.0;
    diagnostics->initial_selected_n[chain] = (int)record[5] - 1;
    diagnostics->initial_attempts[chain] = (uint64_t)record[6];
    diagnostics->initial_exact_zero[chain] = (uint64_t)record[7];
    diagnostics->initial_numeric_zero[chain] = (uint64_t)record[8];
    diagnostics->initial_exhausted[chain] = record[9] != 0.0;
  }
  diagnostics->total_bootstrap_attempts = (uint64_t)totals[0];
  diagnostics->total_bootstrap_exact_zero = (uint64_t)totals[1];
  diagnostics->total_bootstrap_numeric_zero = (uint64_t)totals[2];
  diagnostics->total_bootstrap_local_energy = (uint64_t)totals[3];
  diagnostics->any_bootstrap_exhausted = totals[4] != 0.0;
  diagnostics->total_initial_attempts = (uint64_t)totals[6];
  diagnostics->total_initial_exact_zero = (uint64_t)totals[7];
  diagnostics->total_initial_numeric_zero = (uint64_t)totals[8];
  diagnostics->any_initial_exhausted = totals[9] != 0.0;
  if (failure[0] != 0.0 && isfinite(first_failure)) {
    diagnostics->failure_recorded = 1;
    diagnostics->failure_chain = (uint64_t)first_failure;
    diagnostics->failure_stage = (MVMCPowerLanczosGCStage)failure[1];
    diagnostics->failure_status = (MVMCKrylovStatus)failure[2];
    diagnostics->failure_electron_count = (int)failure[3];
    diagnostics->failure_amplitude_state = (int)failure[4];
    diagnostics->failure_log_abs = failure[5];
    diagnostics->failure_log_error_bound = failure[6];
    diagnostics->failure_retry_reason = (int)failure[7];
  }
  diagnostics->coefficient_available = gc->coefficient_available;
  diagnostics->coefficient[0] = gc->coefficient[0];
  diagnostics->coefficient[1] = gc->coefficient[1];
  diagnostics->coefficient_energy_available =
      gc->coefficient_energy_available;
  diagnostics->coefficient_energy = gc->coefficient_energy;
  diagnostics->valid = 1;
}

MVMCKrylovStatus mvmc_power_lanczos_gc_run_diagnosed(
    const MVMCPowerLanczosGCInput *input,
    MVMCPowerLanczosStabilizedResult *result,
    MVMCPowerLanczosGCDiagnostics *diagnostics) {
  MVMCPowerLanczosStabilizedInput core;
  MVMCPowerLanczosGCRun gc;
  MVMCKrylovStatus status;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (input == NULL) {
    invalidate_result(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&gc, 0, sizeof(gc));
  gc.initial_electron_count = input->initial_electron_count;
  gc.initial_electron_count_explicit = input->initial_electron_count_explicit;
  gc.initial_selected_n = -1;
  gc.failure_log_abs = NAN;
  gc.failure_log_error_bound = NAN;
  if (diagnostics != NULL && diagnostics->final_sector_samples != NULL &&
      diagnostics->final_sector_capacity != 0) {
    gc.sector_samples = diagnostics->final_sector_samples;
    gc.sector_capacity = diagnostics->final_sector_capacity;
    memset(gc.sector_samples, 0,
           gc.sector_capacity * sizeof(*gc.sector_samples));
  }
#if defined(MVMC_ENABLE_FAULT_INJECTION)
  memset(&gc_testing_counters, 0, sizeof(gc_testing_counters));
  {
    const char *value = getenv("MVMC_GC_B0_LANCZOS_TRACE");
    gc.b0_enabled = value != NULL && strcmp(value, "1") == 0;
  }
#endif
  memset(&core, 0, sizeof(core));
  core.world_communicator = input->world_communicator;
  core.chain_communicator = input->chain_communicator;
  core.power_step = 1;
  core.seed = input->seed;
  core.warm_up = input->warm_up;
  core.sample_count = input->sample_count;
  core.interval = input->interval;
  status = stabilized_core(&core, input->gc_view, &gc, result);
#if defined(MVMC_HAVE_MPFR_REFINEMENT)
  if (gc.topology_ready && gc.chain_rank == 0 && gc.refinement_calls != 0)
    fprintf(stdout, "GC Lanczos refinement: chain=%llu calls=%llu "
            "accepted=%llu exact_zero=%llu max_precision=%u\n",
            (unsigned long long)gc.chain_index,
            (unsigned long long)gc.refinement_calls,
            (unsigned long long)gc.refinement_accepted,
            (unsigned long long)gc.refinement_exact_zero,
            gc.refinement_max_precision);
#endif
#if defined(MVMC_ENABLE_FAULT_INJECTION)
  gc_b0_finish(&gc, status);
#endif
  if (diagnostics != NULL) {
    gc_collect_diagnostics(&gc, input->world_communicator, diagnostics);
  }
  free(gc.initial_records);
  free(gc.initial_sectors);
  return status;
}

MVMCKrylovStatus mvmc_power_lanczos_gc_run(
    const MVMCPowerLanczosGCInput *input,
    MVMCPowerLanczosStabilizedResult *result) {
  return mvmc_power_lanczos_gc_run_diagnosed(input, result, NULL);
}

#if defined(MVMC_ENABLE_FAULT_INJECTION)
static int words_overlap(const uint64_t *a, const uint64_t *b, size_t count) {
  const uintptr_t x = (uintptr_t)a, y = (uintptr_t)b;
  const uintptr_t bytes = (uintptr_t)(count * sizeof(*a));
  return a != NULL && b != NULL && x < y + bytes && y < x + bytes;
}

/* Local shape/pointer/enum/alias checks of the script and trace. */
static MVMCKrylovStatus final_script_check(const MVMCGCFinalScript *script,
                                           const MVMCGCFinalStepTrace *trace,
                                           size_t capacity) {
  const MVMCPowerLanczosGCView *view = script->gc_view;
  size_t k, j, word_count;
  int c;
  if (view == NULL || view->data.site_count <= 0 ||
      view->data.site_count > INT_MAX / 2 || script->initial_words == NULL ||
      script->step_count == 0 || script->step_count > capacity ||
      script->steps == NULL || trace == NULL ||
      script->step_count > SIZE_MAX / 4) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  word_count = mvmc_krylov_fock_word_count((size_t)view->data.site_count);
  if (word_count == 0 || script->word_count != word_count ||
      word_count > SIZE_MAX / sizeof(uint64_t) / 4) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (c = 0; c < 2; ++c) {
    if (!isfinite(creal(script->coefficient[c])) ||
        !isfinite(cimag(script->coefficient[c]))) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  for (k = 0; k < script->step_count; ++k) {
    const MVMCGCFinalScriptStep *step = script->steps + k;
    if (step->proposal_words == NULL || !isfinite(step->uniform) ||
        step->uniform < 0.0 || step->uniform >= 1.0 ||
        (step->force != MVMC_GC_FINAL_FORCE_NONE &&
         step->force != MVMC_GC_FINAL_FORCE_NUMERIC_ZERO &&
         step->force != MVMC_GC_FINAL_FORCE_LOCAL_ENERGY) ||
        trace[k].before_words == NULL || trace[k].after_words == NULL ||
        words_overlap(trace[k].before_words, trace[k].after_words,
                      word_count) ||
        words_overlap(trace[k].before_words, script->initial_words,
                      word_count) ||
        words_overlap(trace[k].after_words, script->initial_words,
                      word_count)) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    for (j = 0; j < script->step_count; ++j) {
      if (words_overlap(trace[k].before_words, script->steps[j].proposal_words,
                        word_count) ||
          words_overlap(trace[k].after_words, script->steps[j].proposal_words,
                        word_count) ||
          (j != k &&
           (words_overlap(trace[k].before_words, trace[j].before_words,
                          word_count) ||
            words_overlap(trace[k].before_words, trace[j].after_words,
                          word_count) ||
            words_overlap(trace[k].after_words, trace[j].after_words,
                          word_count)))) {
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      }
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

/* Coefficients, shape, configurations, uniform bits and forces; never a
 * pointer or padding. */
static uint64_t final_script_hash(const MVMCGCFinalScript *script) {
  uint64_t hash = UINT64_C(1469598103934665603), bits;
  double parts[4];
  size_t k, w;
  parts[0] = creal(script->coefficient[0]);
  parts[1] = cimag(script->coefficient[0]);
  parts[2] = creal(script->coefficient[1]);
  parts[3] = cimag(script->coefficient[1]);
  hash = hash_u64(hash, UINT64_C(0x4743464e414c5331)); /* GCFNALS1 */
  for (k = 0; k < 4; ++k) {
    memcpy(&bits, parts + k, sizeof(bits));
    hash = hash_u64(hash, bits);
  }
  hash = hash_u64(hash, (uint64_t)script->word_count);
  hash = hash_u64(hash, (uint64_t)script->step_count);
  for (w = 0; w < script->word_count; ++w) {
    hash = hash_u64(hash, script->initial_words[w]);
  }
  for (k = 0; k < script->step_count; ++k) {
    for (w = 0; w < script->word_count; ++w) {
      hash = hash_u64(hash, script->steps[k].proposal_words[w]);
    }
    memcpy(&bits, &script->steps[k].uniform, sizeof(bits));
    hash = hash_u64(hash, bits);
    hash = hash_u64(hash, (uint64_t)(unsigned)script->steps[k].force);
  }
  return hash == 0 ? UINT64_C(1) : hash;
}

/* Every rank holds the same value (u64 pair through MIN/MAX). */
static MVMCKrylovStatus world_all_equal_u64(
    uint64_t value, MVMCClassicPfaffianCommunicator world,
    int *all_equal) {
#ifdef _mpi_use
  uint64_t minimum = 0, maximum = 0;
  if (MPI_Allreduce(&value, &minimum, 1, MPI_UINT64_T, MPI_MIN, world) !=
          MPI_SUCCESS ||
      MPI_Allreduce(&value, &maximum, 1, MPI_UINT64_T, MPI_MAX, world) !=
          MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  *all_equal = minimum == maximum;
#else
  (void)value;
  (void)world;
  *all_equal = 1;
#endif
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_gc_testing_final_steps(
    const MVMCGCFinalScript *script, MVMCGCFinalStepTrace *trace,
    size_t trace_capacity, size_t *written) {
  MVMCPowerLanczosGCRun gc;
  MVMCPowerLanczosSamplingTopology sampling;
  MVMCPowerLanczosGCBridge *bridge = NULL;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  MVMCKrylovFinalStatePolicy policy;
  MVMCKrylovFinalStateSnapshot snapshot;
  const MVMCKrylovFockModel *model = NULL;
  MVMCKrylovScaledAmplitudeCallback amplitude = NULL;
  void *amplitude_context = NULL;
  uint64_t *words = NULL, *proposal = NULL;
  uint64_t generation_hash = 0, script_hash = 0;
  MVMCKrylovStatus status;
  size_t k, count = 0;
  int equal = 0;
  /* RR-1: counters start at zero; only refinement-entry faults apply.
   * Without a script there is no communicator to agree on. */
  memset(&gc_testing_counters, 0, sizeof(gc_testing_counters));
  if (script == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  status = written == NULL ? MVMC_KRYLOV_STATUS_INVALID_ARGUMENT
                           : final_script_check(script, trace,
                                                trace_capacity);
  if (gc_testing_fault != MVMC_POWER_LANCZOS_GC_FAULT_NONE &&
      gc_testing_fault != MVMC_POWER_LANCZOS_GC_FAULT_REFINEMENT_DECISION &&
      gc_testing_fault != MVMC_POWER_LANCZOS_GC_FAULT_REFINEMENT_STATUS) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (written != NULL) *written = 0;
  status = synchronize_status(status, script->world_communicator);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = world_all_equal_u64((uint64_t)(unsigned)gc_testing_fault,
                                 script->world_communicator, &equal);
    if (status == MVMC_KRYLOV_STATUS_OK && !equal) {
      status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    script_hash = final_script_hash(script);
    status = world_all_equal_u64(script_hash, script->world_communicator,
                                 &equal);
    if (status == MVMC_KRYLOV_STATUS_OK && !equal) {
      status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  if (status != MVMC_KRYLOV_STATUS_OK) return status;

  memset(&gc, 0, sizeof(gc));
  gc.initial_selected_n = -1;
  gc.failure_log_abs = NAN;
  gc.failure_log_error_bound = NAN;
  status = mvmc_power_lanczos_sampling_topology_create(
      script->world_communicator, script->chain_communicator, 1, &sampling);
  status = synchronize_status(status, script->world_communicator);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_gc_bridge_create(
        script->gc_view, script->world_communicator,
        script->chain_communicator, &bridge);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    model = mvmc_power_lanczos_gc_bridge_model(bridge);
    amplitude = mvmc_power_lanczos_gc_bridge_amplitude(bridge);
    amplitude_context = mvmc_power_lanczos_gc_bridge_amplitude_context(bridge);
    generation_hash = mvmc_power_lanczos_gc_bridge_generation_hash(bridge);
    gc.model = model;
    gc.spin_support = script->gc_view->spin_support;
    gc.chain_communicator = script->chain_communicator;
    gc.chain_rank = sampling.chain_rank;
    gc.chain_index = sampling.chain_index;
    gc.topology_ready = 1;
    gc.stage = MVMC_POWER_LANCZOS_GC_STAGE_FINAL_SAMPLE;
    words = (uint64_t *)calloc(script->word_count, sizeof(*words));
    proposal = (uint64_t *)calloc(script->word_count, sizeof(*proposal));
    gc.step_words = (uint64_t *)calloc(script->word_count, sizeof(*words));
    status = words == NULL || proposal == NULL || gc.step_words == NULL
                 ? MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE
                 : create_limits(model, 1, generation_hash, &limits);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_bounded_krylov_plan_create(model, &limits, &plan);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_bounded_krylov_workspace_create(plan, &workspace);
    }
    /* Raw coefficients as given: no scaled basis. */
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_krylov_final_state_policy_create(
          1, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC, script_hash,
          script->coefficient, 2, &policy);
    }
    status = synchronize_status(status, script->world_communicator);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_session_begin(
        workspace, amplitude, amplitude_context, generation_hash);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    memcpy(words, script->initial_words,
           script->word_count * sizeof(*words));
    status = mvmc_krylov_final_state_sampler_initialize(
        workspace, &policy, words, script->word_count, amplitude,
        amplitude_context, &snapshot);
  }
  status = synchronize_status(status, script->world_communicator);
  /* Chain-local loop: no world collective until it ends. */
  for (k = 0; status == MVMC_KRYLOV_STATUS_OK && k < script->step_count;
       ++k) {
    MVMCGCFinalStepTrace *entry = trace + k;
    MVMCKrylovFinalDiagnostic diagnostic;
    MVMCPowerLanczosGCStepOutcome outcome;
    double uniform = 0.0;
    int force = MVMC_GC_FINAL_FORCE_NONE;
    if (gc.chain_rank == 0) {
      memcpy(proposal, script->steps[k].proposal_words,
             script->word_count * sizeof(*proposal));
      uniform = script->steps[k].uniform;
      force = (int)script->steps[k].force;
    }
    status = gc_broadcast(&gc, MVMC_KRYLOV_STATUS_OK, proposal,
                          script->word_count, &uniform);
#ifdef _mpi_use
    if (status == MVMC_KRYLOV_STATUS_OK &&
        MPI_Bcast(&force, 1, MPI_INT, 0, gc.chain_communicator) !=
            MPI_SUCCESS) {
      status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
#endif
    /* The step index and the force are agreed inside the chain. */
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = gc_agree_value(&gc, MVMC_KRYLOV_STATUS_OK, (int)k);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = gc_agree_value(&gc, MVMC_KRYLOV_STATUS_OK, force);
    }
    if (status != MVMC_KRYLOV_STATUS_OK) break;
    ++gc_testing_counters.generated[gc.stage];
    memcpy(entry->before_words, words, script->word_count * sizeof(*words));
    entry->before = snapshot;
    status = gc_final_proposal_step(&gc, workspace, &policy, model, words,
                                    script->word_count, &snapshot, amplitude,
                                    amplitude_context, proposal, uniform,
                                    force, &diagnostic, &outcome);
    memcpy(entry->after_words, words, script->word_count * sizeof(*words));
    entry->after = snapshot;
    entry->status = status;
    entry->accepted = outcome.accepted;
    entry->used_refinement = outcome.used_refinement;
    entry->exact_zero_proposal = outcome.exact_zero_proposal;
    entry->refinement_calls = gc.refinement_calls;
    entry->refinement_accepted = gc.refinement_accepted;
    entry->refinement_exact_zero = gc.refinement_exact_zero;
    entry->precision = outcome.precision;
    entry->node_calls = outcome.node_calls;
    entry->proposal_copies = (uint64_t)k + 1;
    entry->uniform_copies = (uint64_t)k + 1;
    entry->rng_draws = 0;
    entry->electron_count = popcount_words(words, script->word_count);
    entry->electron_count_square =
        (uint64_t)entry->electron_count * (uint64_t)entry->electron_count;
    count = k + 1;
  }
  if (workspace != NULL && mvmc_bounded_krylov_session_is_active(workspace)) {
    const MVMCKrylovStatus end = mvmc_bounded_krylov_session_end(workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) status = end;
  }
  status = synchronize_status(status, script->world_communicator);
  *written = count;
  free(proposal);
  free(gc.step_words);
  free(words);
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
  mvmc_power_lanczos_gc_bridge_destroy(bridge);
  return status;
}
#endif
