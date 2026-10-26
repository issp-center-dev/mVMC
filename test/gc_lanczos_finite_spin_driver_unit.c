/*
 * GC first-step Lanczos driver with a finite-S (S > 0) spin support, on the
 * frozen Task 5 fixtures (E1: L4/S1/G8 general, R1_G2: L2/S1/G2 triplet):
 * auto/explicit start sectors (a forbidden explicit N is agreed before any
 * draw), the 100-per-N structured start with exact/numeric-zero redraws
 * (100th success, next N, exhaustion, non-retry stop, rank disagreement),
 * bootstrap-RNG-only consumption, chain-local retry counts merged before
 * the pilot, diagnostics, and the final bootstrap (structured selected-N /
 * whole-domain alternation, 1024th success, local-energy retry, stop on
 * other errors, Psi = 0 starts with Phi != 0).  Topologies: serial, np2
 * (one chain of 2 or two of 1), np3 (2+1 with an empty QP rank), np4 (2x2).
 */
#include "gc_lanczos_task5_case.h"
#include "krylov_final_state.h"
#include "power_lanczos_stabilized.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failure_count = 0;

#define CHECK(condition, message)                                    \
  do {                                                               \
    if (!(condition)) {                                              \
      fprintf(stderr, "FAIL: %s (line %d)\n", (message), __LINE__); \
      ++failure_count;                                               \
    }                                                                \
  } while (0)

static int world_rank(void) {
  int rank = 0;
#ifdef _mpi_use
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  return rank;
}

static int world_size(void) {
  int size = 1;
#ifdef _mpi_use
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  return size;
}

static MVMCClassicPfaffianCommunicator world(void) {
#ifdef _mpi_use
  return MPI_COMM_WORLD;
#else
  return 0;
#endif
}

typedef struct {
  MVMCClassicPfaffianCommunicator chain;
  int chain_rank;
  int chain_size;
  int owns;
} Chain;

/* 0: one chain of every rank; 1: one chain per rank; 2: chains of two. */
static Chain make_chain(int mode) {
  Chain chain;
  memset(&chain, 0, sizeof(chain));
  chain.chain_size = 1;
#ifdef _mpi_use
  if (mode == 0) {
    chain.chain = MPI_COMM_WORLD;
  } else if (mode == 1) {
    chain.chain = MPI_COMM_SELF;
  } else {
    MPI_Comm_split(MPI_COMM_WORLD, world_rank() / 2, world_rank(),
                   &chain.chain);
    chain.owns = 1;
  }
  MPI_Comm_rank(chain.chain, &chain.chain_rank);
  MPI_Comm_size(chain.chain, &chain.chain_size);
#else
  (void)mode;
  chain.chain = 0;
#endif
  return chain;
}

static void free_chain(Chain *chain) {
#ifdef _mpi_use
  if (chain->owns) MPI_Comm_free(&chain->chain);
#else
  (void)chain;
#endif
}

static uint64_t chain_count_for(int mode) {
  if (mode == 0) return 1;
  if (mode == 1) return (uint64_t)world_size();
  return (uint64_t)(world_size() + 1) / 2;
}

typedef struct {
  MVMCKrylovStatus status;
  MVMCPowerLanczosStabilizedResult result;
  MVMCPowerLanczosGCDiagnostics diagnostics;
  MVMCPowerLanczosGCTestingCounters counters;
  uint64_t sectors[9];
} Run;

typedef struct {
  const Task5Case *source;
  int symmetry; /* -1: the fixture's own */
  int initial_n;
  int explicit_n;
  uint64_t seed;
  size_t warm_up, samples, interval;
  MVMCPowerLanczosGCFault fault;
  MVMCPowerLanczosGCInitialFault initial_fault;
  uint64_t initial_count, initial_chain;
  int rank1_initial_n; /* -1: same as everyone */
} RunSpec;

static RunSpec spec_for(const Task5Case *source, int initial_n,
                        int explicit_n) {
  RunSpec spec;
  memset(&spec, 0, sizeof(spec));
  spec.source = source;
  spec.symmetry = -1;
  spec.initial_n = initial_n;
  spec.explicit_n = explicit_n;
  spec.seed = 17;
  spec.warm_up = 4;
  spec.samples = 16;
  spec.interval = 1;
  spec.initial_chain = UINT64_MAX;
  spec.rank1_initial_n = -1;
  return spec;
}

static void run_gc(const Chain *chain, const RunSpec *spec, Run *run) {
  Task5Built built;
  MVMCPowerLanczosGCView view;
  MVMCPowerLanczosGCInput input;
  memset(run, 0, sizeof(*run));
  if (!task5_build(spec->source, &built)) {
    run->status = MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    return;
  }
  if (spec->symmetry >= 0) {
    built.support.pair_symmetry = (MVMCGCPairSymmetry)spec->symmetry;
  }
  task5_view(&built, world_rank(), chain->chain_rank, chain->chain_size,
             &view);
  memset(&input, 0, sizeof(input));
  input.gc_view = &view;
  input.world_communicator = world();
  input.chain_communicator = chain->chain;
  input.seed = spec->seed;
  input.warm_up = spec->warm_up;
  input.sample_count = spec->samples;
  input.interval = spec->interval;
  input.initial_electron_count =
      spec->rank1_initial_n >= 0 && world_rank() == 1 ? spec->rank1_initial_n
                                                      : spec->initial_n;
  input.initial_electron_count_explicit = spec->explicit_n;
  run->diagnostics.final_sector_samples = run->sectors;
  run->diagnostics.final_sector_capacity =
      2 * (size_t)spec->source->site_count + 1;
  mvmc_power_lanczos_gc_testing_set_fault(spec->fault);
  mvmc_power_lanczos_gc_testing_set_initial_fault(
      spec->initial_fault, spec->initial_count, spec->initial_chain);
  run->status = mvmc_power_lanczos_gc_run_diagnosed(&input, &run->result,
                                                    &run->diagnostics);
  mvmc_power_lanczos_gc_testing_set_fault(MVMC_POWER_LANCZOS_GC_FAULT_NONE);
  mvmc_power_lanczos_gc_testing_set_initial_fault(
      MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_NONE, 0, 0);
  mvmc_power_lanczos_gc_testing_counters(&run->counters);
  task5_free(&built);
}

static int all_ranks_equal_int(int value) {
#ifdef _mpi_use
  int minimum = value, maximum = value;
  MPI_Allreduce(&value, &minimum, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&value, &maximum, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return minimum == maximum;
#else
  (void)value;
  return 1;
#endif
}

/* A full run completed.  Without the refinement backend a NUMERIC_ZERO or
 * unavailable local energy of a final proposal ends the run (fail closed,
 * as before Task 5); the start and bootstrap under test happened earlier. */
static int completed(const Run *run) {
#if defined(MVMC_HAVE_MPFR_REFINEMENT)
  return run->status == MVMC_KRYLOV_STATUS_OK && run->result.valid;
#else
  if (run->status == MVMC_KRYLOV_STATUS_OK) return run->result.valid;
  return run->status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
         (world_rank() != 0 ||
          run->diagnostics.failure_stage ==
              MVMC_POWER_LANCZOS_GC_STAGE_FINAL_WARM_UP ||
          run->diagnostics.failure_stage ==
              MVMC_POWER_LANCZOS_GC_STAGE_FINAL_SAMPLE);
#endif
}

static void report(const char *label, const Run *run) {
  if (world_rank() != 0) return;
  fprintf(stderr,
          "%s: status=%s chains=%llu selected0=%d attempts0=%llu "
          "exact0=%llu numeric0=%llu exhausted=%d failure=%d/%s/%s\n",
          label, mvmc_krylov_status_string(run->status),
          (unsigned long long)run->diagnostics.chain_count,
          run->diagnostics.initial_selected_n[0],
          (unsigned long long)run->diagnostics.initial_attempts[0],
          (unsigned long long)run->diagnostics.initial_exact_zero[0],
          (unsigned long long)run->diagnostics.initial_numeric_zero[0],
          run->diagnostics.any_initial_exhausted,
          run->diagnostics.failure_recorded,
          mvmc_power_lanczos_gc_stage_name(run->diagnostics.failure_stage),
          mvmc_krylov_status_string(run->diagnostics.failure_status));
}

/* Sector selection: auto order, explicit keep/reject before any draw. */
static void test_sectors(int mode) {
  Chain chain = make_chain(mode);
  const uint64_t chains = chain_count_for(mode);
  const int root = chain.chain_rank == 0;
  Run run;
  RunSpec spec = spec_for(&task5_e1, 4, 0);
  uint64_t chain_index;
  run_gc(&chain, &spec, &run);
  if (!completed(&run)) report("E1 auto N4", &run);
  CHECK(completed(&run) &&
            (run.status != MVMC_KRYLOV_STATUS_OK ||
             run.result.sampling_chains == chains),
        "E1 finite-S run with auto N4 completes");
  CHECK(run.counters.coefficient_rng_draws_at_pilot == 0 &&
            run.counters.final_rng_draws_at_bootstrap == 0,
        "the start uses neither the coefficient nor the final RNG");
  CHECK(run.counters.generated[MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_INITIAL] ==
            (root ? run.counters.initial_retries + 1 : 0u),
        "root draws one candidate per start attempt");
  CHECK(run.counters.bootstrap_rng_draws_after_start > 0 || !root,
        "the start draws from the bootstrap RNG on the chain root");
  CHECK(run.counters.bootstrap_initial_off_sector == 0 &&
            (!root || (run.counters.bootstrap_initial_draws >= 1 &&
                       run.counters.bootstrap_uniform_draws + 1 >=
                           run.counters.bootstrap_initial_draws)),
        "final bootstrap initial side stays in the selected sector");
  if (world_rank() == 0) {
    for (chain_index = 0; chain_index < chains && chain_index < 16;
         ++chain_index) {
      CHECK(run.diagnostics.initial_selected_n[chain_index] == 4 &&
                run.diagnostics.initial_attempts[chain_index] >= 1 &&
                !run.diagnostics.initial_exhausted[chain_index],
            "every chain selects the nearest allowed N4");
    }
    CHECK(run.diagnostics.initial_selected_n[chains < 16 ? chains : 15] ==
              (chains < 16 ? -1 : 4),
          "an absent chain slot decodes as unselected");
  }
  /* A triplet descriptor forbids N4: auto takes [2, 6], explicit N4 fails
   * before any start draw on every rank, explicit N2 is kept. */
  spec.symmetry = MVMC_GC_PAIR_TRIPLET;
  run_gc(&chain, &spec, &run);
  CHECK(completed(&run), "triplet auto N4 run");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.initial_selected_n[0] == 2,
          "triplet auto N4 selects N2 first");
  }
  spec.explicit_n = 1;
  run_gc(&chain, &spec, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !run.result.valid &&
            all_ranks_equal_int((int)run.status) &&
            run.counters.generated[MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_INITIAL] ==
                0 &&
            run.counters.bootstrap_rng_draws_after_start == 0 &&
            run.counters.steps[MVMC_POWER_LANCZOS_GC_STAGE_PILOT_WARM_UP] == 0,
        "explicit forbidden N4 is rejected on every rank before any draw");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.initial_selected_n[0] == -1 &&
              run.diagnostics.total_initial_attempts == 0,
          "no start attempt for a rejected explicit N");
  }
  spec.initial_n = 2;
  run_gc(&chain, &spec, &run);
  CHECK(completed(&run), "explicit N2 run");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.initial_selected_n[0] == 2,
          "explicit allowed N2 is kept");
  }
  /* Range violation (2S <= N <= 2L-2S) with the general descriptor. */
  spec = spec_for(&task5_e1, 0, 1);
  run_gc(&chain, &spec, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            all_ranks_equal_int((int)run.status),
        "explicit N0 outside 2S <= N is rejected");
  spec = spec_for(&task5_e1, 3, 0);
  run_gc(&chain, &spec, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "odd requested N is rejected, never rounded");
  free_chain(&chain);
}

/* 100 attempts per N, next N, exhaustion, non-retry stop. */
static void test_retries(void) {
  Chain chain = make_chain(0);
  const int root = chain.chain_rank == 0;
  Run base, run;
  RunSpec spec = spec_for(&task5_e1, 4, 0);
  run_gc(&chain, &spec, &base);
  CHECK(completed(&base), "E1 base run");

  spec.initial_fault = MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_EXACT_ZERO;
  spec.initial_count = 99;
  run_gc(&chain, &spec, &run);
  if (!completed(&run)) report("99 exact zeros", &run);
  CHECK(completed(&run), "the 100th attempt succeeds");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.initial_selected_n[0] == 4 &&
              run.diagnostics.initial_attempts[0] == 100 &&
              run.diagnostics.initial_exact_zero[0] == 99 &&
              run.diagnostics.initial_numeric_zero[0] == 0,
          "99 exact-zero redraws then N4 at the 100th attempt");
  }
  CHECK(run.counters.coefficient_rng_draws_at_pilot == 0 &&
            run.counters.final_rng_draws_at_bootstrap == 0 &&
            (!root || run.counters.bootstrap_rng_draws_after_start >
                          base.counters.bootstrap_rng_draws_after_start),
        "retries consume the bootstrap RNG only");

  spec.initial_fault = MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_NUMERIC_ZERO;
  spec.initial_count = 100;
  run_gc(&chain, &spec, &run);
  CHECK(completed(&run), "100 numeric zeros, next N");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.initial_selected_n[0] == 2 &&
              run.diagnostics.initial_attempts[0] == 101 &&
              run.diagnostics.initial_numeric_zero[0] == 100 &&
              !run.diagnostics.initial_exhausted[0],
          "N4 exhausted after 100 attempts, N2 selected next");
  }

  spec.initial_fault = MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_EXACT_ZERO;
  spec.initial_count = 300;
  run_gc(&chain, &spec, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
            !run.result.valid && all_ranks_equal_int((int)run.status) &&
            run.counters.steps[MVMC_POWER_LANCZOS_GC_STAGE_PILOT_WARM_UP] == 0,
        "every allowed N exhausted: AMPLITUDE_FAILURE before the pilot");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.initial_exhausted[0] &&
              run.diagnostics.any_initial_exhausted &&
              run.diagnostics.initial_attempts[0] == 300 &&
              run.diagnostics.initial_exact_zero[0] == 300 &&
              run.diagnostics.initial_selected_n[0] == -1 &&
              run.diagnostics.failure_recorded &&
              run.diagnostics.failure_stage ==
                  MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_INITIAL,
          "exhaustion diagnostics over [4, 2, 6]");
  }

  spec.initial_fault = MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_NONFINITE;
  spec.initial_count = 0;
  run_gc(&chain, &spec, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_NONFINITE && !run.result.valid &&
            all_ranks_equal_int((int)run.status),
        "a nonfinite start candidate is not retried");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.initial_attempts[0] == 1 &&
              !run.diagnostics.initial_exhausted[0] &&
              run.diagnostics.initial_selected_n[0] == -1,
          "the stop happens at that attempt");
  }
  if (chain.chain_size > 1) {
    spec.initial_fault = MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_STATUS;
    run_gc(&chain, &spec, &run);
    CHECK(run.status == MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE &&
              all_ranks_equal_int((int)run.status),
          "chain rank 1 disagreement fails on every rank without a hang");
  }
  /* World input agreement. */
  if (world_size() > 1) {
    spec = spec_for(&task5_e1, 4, 0);
    spec.rank1_initial_n = 2;
    run_gc(&chain, &spec, &run);
    CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              all_ranks_equal_int((int)run.status),
          "rank-mismatched initial N rejected on every rank");
  }
  spec = spec_for(&task5_e1, 4, 2);
  run_gc(&chain, &spec, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "explicit flag must be 0 or 1");
  free_chain(&chain);
}

/* Different chains retry independently and meet before the pilot. */
static void test_chain_local(void) {
  Chain chain;
  Run run;
  RunSpec spec = spec_for(&task5_e1, 4, 0);
  if (world_size() < 2) return;
  chain = make_chain(world_size() >= 4 ? 2 : 1);
  spec.initial_fault = MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_EXACT_ZERO;
  spec.initial_count = 5;
  spec.initial_chain = 1;
  run_gc(&chain, &spec, &run);
  CHECK(completed(&run), "chain-local retries merge at the pilot");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.initial_attempts[0] == 1 &&
              run.diagnostics.initial_attempts[1] == 6 &&
              run.diagnostics.initial_exact_zero[1] == 5 &&
              run.diagnostics.total_initial_attempts >= 7,
          "chain 1 retried five times, chain 0 not");
  }
  /* A non-retry stop on chain 1 only still fails every rank. */
  spec.initial_fault = MVMC_POWER_LANCZOS_GC_INITIAL_FAULT_NONFINITE;
  spec.initial_count = 0;
  run_gc(&chain, &spec, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_NONFINITE &&
            all_ranks_equal_int((int)run.status),
        "one chain's start failure is agreed by the world");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.failure_recorded && run.diagnostics.failure_chain == 1,
          "the failing chain is reported");
  }
  free_chain(&chain);
}

static void test_final_bootstrap(void) {
  Chain chain = make_chain(0);
  const int root = chain.chain_rank == 0;
  Run run;
  RunSpec spec = spec_for(&task5_e1, 4, 0);
  uint64_t seed;
  int psi_zero_starts = 0;
  spec.fault = MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_EXHAUSTED;
  run_gc(&chain, &spec, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
            (!root || (run.counters.bootstrap_initial_draws == 512 &&
                       run.counters.bootstrap_uniform_draws == 512 &&
                       run.counters.bootstrap_initial_off_sector == 0)),
        "1024 exhausted: alternating structured N4 and whole-domain draws");
  spec.fault = MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_LAST;
  run_gc(&chain, &spec, &run);
  CHECK(completed(&run), "the 1024th bootstrap candidate succeeds");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.total_bootstrap_attempts == 1024 &&
              run.diagnostics.total_bootstrap_numeric_zero == 1023 &&
              !run.diagnostics.any_bootstrap_exhausted,
          "1023 numeric-zero redraws then success");
  }
  spec.fault = MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_LOCAL_ENERGY;
  run_gc(&chain, &spec, &run);
  CHECK(completed(&run), "local-energy retry");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.total_bootstrap_local_energy == 1 &&
              run.diagnostics.total_bootstrap_attempts >= 2,
          "an unavailable local energy is redrawn");
  }
  spec.fault = MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_NONFINITE;
  run_gc(&chain, &spec, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_NONFINITE && !run.result.valid &&
            all_ranks_equal_int((int)run.status),
        "another bootstrap error stops the run");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.total_bootstrap_attempts == 1 &&
              run.diagnostics.failure_stage ==
                  MVMC_POWER_LANCZOS_GC_STAGE_FINAL_BOOTSTRAP,
          "the stop is at the first bootstrap candidate");
  }
  /* E1 (H = N): Phi = 0 exactly where Psi has fewer than 2S singles.  The
   * first (structured) candidate is made a numeric zero so that the next
   * ones come from the whole domain; Phi = 0 candidates are redrawn and
   * the accepted start always has Phi != 0. */
  for (seed = 1; seed <= 24; ++seed) {
    spec = spec_for(&task5_e1, 4, 0);
    spec.seed = seed;
    spec.fault = MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_NUMERIC;
    run_gc(&chain, &spec, &run);
    CHECK(completed(&run), "E1 finite-S run");
    if (!completed(&run)) {
      report("E1 bootstrap", &run);
      continue;
    }
    {
      const uint64_t w = run.counters.final_start_word;
      const int singles = __builtin_popcountll((w ^ (w >> 4)) & 0xfu);
      CHECK(singles >= 2, "the final start is inside the S=1 support");
    }
    if (world_rank() == 0) {
      psi_zero_starts += run.diagnostics.total_bootstrap_exact_zero > 0;
    }
  }
  if (world_rank() == 0) {
    CHECK(psi_zero_starts > 0, "Phi = 0 whole-domain candidates are redrawn");
  }
  /* The bootstrap judges Phi, not Psi: R1_G2 Psi(5) = 0 but Phi(5) != 0. */
  {
    Task5Built built;
    MVMCGCKrylovAmplitudeLayout layout;
    MVMCGCKrylovAmplitudeWorkspace *workspace = NULL;
    MVMCKrylovBoundedLimits limits;
    MVMCKrylovBoundedPlan *plan = NULL;
    MVMCKrylovBoundedWorkspace *engine = NULL;
    MVMCKrylovFinalStatePolicy policy;
    MVMCKrylovFinalStateSnapshot snapshot;
    MVMCKrylovFinalDiagnostic diagnostic;
    uint64_t words;
    CHECK(task5_build(&task5_r1_g2, &built), "R1_G2 build");
    task5_layout(&built, world(), world_rank(), world_size(), &layout);
    CHECK(mvmc_gc_krylov_amplitude_workspace_create(
              &layout, built.slater, built.weight, &workspace) ==
              MVMC_KRYLOV_STATUS_OK,
          "R1_G2 workspace");
    memset(&limits, 0, sizeof(limits));
    limits.amplitude_policy_hash =
        mvmc_gc_krylov_amplitude_generation_hash(workspace);
    limits.cache_bytes = (size_t)1 << 16;
    limits.max_row_transitions = built.model.term_count;
    limits.max_workspace_bytes = (size_t)64 * 1024 * 1024;
    limits.max_node_expansions = 1000;
    limits.max_terminal_amplitude_calls = 1000;
    limits.max_total_row_transitions = 10000;
    limits.max_order = 2;
    CHECK(workspace != NULL &&
              mvmc_bounded_krylov_plan_create(&built.model, &limits, &plan) ==
                  MVMC_KRYLOV_STATUS_OK &&
              mvmc_bounded_krylov_workspace_create(plan, &engine) ==
                  MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_final_state_policy_create(
                  1, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC, 47,
                  built.coefficient, 2, &policy) == MVMC_KRYLOV_STATUS_OK,
          "R1_G2 engine and raw policy");
    if (engine != NULL) {
      words = 5;
      memset(&diagnostic, 0, sizeof(diagnostic));
      CHECK(mvmc_krylov_final_state_sampler_initialize_diagnosed(
                engine, &policy, &words, 1, mvmc_gc_krylov_scaled_amplitude,
                workspace, &snapshot, &diagnostic) == MVMC_KRYLOV_STATUS_OK &&
                snapshot.krylov.value[0].state ==
                    MVMC_SCALED_COMPLEX_EXACT_ZERO &&
                snapshot.final_state.amplitude.state ==
                    MVMC_SCALED_COMPLEX_FINITE_NONZERO,
            "a Psi = 0, Phi != 0 start is accepted");
      words = 0;
      memset(&diagnostic, 0, sizeof(diagnostic));
      CHECK(mvmc_krylov_final_state_sampler_initialize_diagnosed(
                engine, &policy, &words, 1, mvmc_gc_krylov_scaled_amplitude,
                workspace, &snapshot, &diagnostic) ==
                    MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
                diagnostic.retry_reason == MVMC_KRYLOV_BOOTSTRAP_EXACT_ZERO,
            "a Phi = 0 start is an exact-zero retry");
    }
    mvmc_bounded_krylov_workspace_destroy(engine);
    mvmc_bounded_krylov_plan_destroy(plan);
    mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
    task5_free(&built);
  }
  free_chain(&chain);
}

int main(int argc, char **argv) {
  int total = 0;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  test_sectors(0);
  if (world_size() > 1) test_sectors(1);
  if (world_size() > 3) test_sectors(2);
  test_retries();
  test_chain_local();
  test_final_bootstrap();
#ifdef _mpi_use
  MPI_Allreduce(&failure_count, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  total = failure_count;
#endif
  if (world_rank() == 0) {
    printf(total == 0 ? "gc_lanczos_finite_spin_driver_unit: PASS\n"
                      : "gc_lanczos_finite_spin_driver_unit: FAIL\n");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return total == 0 ? 0 : 1;
}
