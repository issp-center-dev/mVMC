/*
 * Testing final-step harness of the GC Lanczos driver on the frozen Task 5
 * fixtures (raw policy c = (1, 1/2), no GEVP/pilot/bootstrap): every
 * scripted variant of the specification table (R1, R1-LE, R1 boundary and
 * rejection, R1C, E1 mixed/all-refined/rejection, S2_G2 and its ordinary/
 * refined mask mutations, R1_G2 with an empty QP rank), the snapshot rule
 * (exact state/decision/configuration/generation/N/N^2 and counters; phase
 * and log_abs within 16 eps refined / 64 eps ordinary; E_local within
 * 1e-14 / 1e-13), comparator mutations, RR-1 fault limits and counter
 * reset, invalid scripts before any step, chain-rank-1 refinement faults
 * (one chain, 2x2 chains) and refinement resource limits.
 */
#include "gc_lanczos_task5_case.h"
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

#define MAX_STEPS 4

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
  int chain_rank, chain_size, owns;
} Chain;

/* 0: one chain of every rank; 2: chains of two ranks. */
static Chain make_chain(int mode) {
  Chain chain;
  memset(&chain, 0, sizeof(chain));
  chain.chain_size = 1;
#ifdef _mpi_use
  if (mode == 0) {
    chain.chain = MPI_COMM_WORLD;
  } else {
    MPI_Comm_split(MPI_COMM_WORLD, world_rank() / 2, world_rank(),
                   &chain.chain);
    chain.owns = 1;
  }
  MPI_Comm_rank(chain.chain, &chain.chain_rank);
  MPI_Comm_size(chain.chain, &chain.chain_size);
#else
  (void)mode;
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

static int all_ranks_equal_u64(uint64_t value) {
#ifdef _mpi_use
  uint64_t minimum = value, maximum = value;
  MPI_Allreduce(&value, &minimum, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&value, &maximum, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  return minimum == maximum;
#else
  (void)value;
  return 1;
#endif
}

typedef struct {
  unsigned to;
  double uniform;
  MVMCGCFinalForce force;
} Step;

typedef struct {
  MVMCKrylovStatus status;
  size_t written;
  MVMCGCFinalStepTrace trace[MAX_STEPS];
  uint64_t before[MAX_STEPS], after[MAX_STEPS];
  MVMCPowerLanczosGCTestingCounters counters;
} Run;

static void run_script(const Task5Built *built, const Chain *chain,
                       double complex c1, unsigned start, const Step *steps,
                       size_t count, size_t capacity, Run *run) {
  MVMCPowerLanczosGCView view;
  MVMCGCFinalScript script;
  MVMCGCFinalScriptStep entries[MAX_STEPS];
  uint64_t proposals[MAX_STEPS];
  uint64_t initial = start;
  size_t k;
  memset(run, 0x5a, sizeof(*run));
  for (k = 0; k < MAX_STEPS; ++k) {
    run->trace[k].before_words = run->before + k;
    run->trace[k].after_words = run->after + k;
  }
  task5_view(built, world_rank(), chain->chain_rank, chain->chain_size,
             &view);
  memset(&script, 0, sizeof(script));
  script.gc_view = &view;
  script.world_communicator = world();
  script.chain_communicator = chain->chain;
  script.coefficient[0] = built->coefficient[0];
  script.coefficient[1] = c1;
  script.initial_words = &initial;
  script.word_count = 1;
  script.step_count = count;
  script.steps = entries;
  for (k = 0; k < count; ++k) {
    proposals[k] = steps[k].to;
    entries[k].proposal_words = proposals + k;
    entries[k].uniform = steps[k].uniform;
    entries[k].force = steps[k].force;
  }
  run->status = mvmc_power_lanczos_gc_testing_final_steps(
      &script, run->trace, capacity, &run->written);
  mvmc_power_lanczos_gc_testing_counters(&run->counters);
  /* One chain of every rank: the same trace on every rank. */
  if (chain->chain_size == world_size()) {
    int same = all_ranks_equal_int((int)run->status) &&
               all_ranks_equal_u64((uint64_t)run->written);
    for (k = 0; same && k < run->written && k < MAX_STEPS; ++k) {
      same = all_ranks_equal_u64(run->trace[k].after.integrity_hash) &&
             all_ranks_equal_u64(run->trace[k].after_words[0]) &&
             all_ranks_equal_int(run->trace[k].accepted) &&
             all_ranks_equal_u64(run->trace[k].refinement_calls);
    }
    CHECK(same, "status, decisions and snapshots agree on every rank");
  }
}

/* ---- snapshot comparison rule (R-5) -------------------------------------- */

typedef struct {
  unsigned words;
  uint64_t generation;
  int electron_count;
  uint64_t electron_count_square;
  int accepted;
  uint64_t calls, accepted_calls, exact_zero_calls;
  int refined;           /* 16 eps / 1e-14 when 1, else 64 eps / 1e-13 */
  const Task5Scaled *table;
  size_t table_count;
  double complex energy;
} Expect;

static int value_matches(const MVMCScaledComplex *value,
                         const Task5Scaled *table, size_t count,
                         unsigned mask, int index, double k) {
  const Task5Scaled *ref = task5_reference(table, count, mask, index);
  return ref != NULL && task5_scaled_matches(value, ref, k);
}

/* Returns 1 when the after-step state satisfies every exact and tolerance
 * condition; the caller compares the decision separately. */
static int step_matches(const MVMCGCFinalStepTrace *entry,
                        const Expect *expect) {
  const MVMCKrylovFinalStateSnapshot *s = &entry->after;
  const double k = expect->refined ? 16.0 : 64.0;
  const double tolerance = expect->refined ? 1e-14 : 1e-13;
  double complex energy = NAN;
  int ok = entry->status == MVMC_KRYLOV_STATUS_OK &&
           entry->after_words[0] == expect->words &&
           s->accepted_generation == expect->generation &&
           entry->electron_count == expect->electron_count &&
           entry->electron_count_square == expect->electron_count_square &&
           entry->accepted == expect->accepted &&
           entry->refinement_calls == expect->calls &&
           entry->refinement_accepted == expect->accepted_calls &&
           entry->refinement_exact_zero == expect->exact_zero_calls &&
           s->integrity_hash ==
               mvmc_krylov_final_state_snapshot_integrity_hash(s);
  if (expect->table != NULL) {
    int v;
    for (v = 0; v < 3; ++v) {
      ok = ok && value_matches(s->krylov.value + v, expect->table,
                               expect->table_count, expect->words, v, k);
    }
    ok = ok && value_matches(&s->final_state.amplitude, expect->table,
                             expect->table_count, expect->words, 3, k) &&
         value_matches(&s->final_state.hamiltonian_amplitude, expect->table,
                       expect->table_count, expect->words, 4, k);
  }
  ok = ok && mvmc_krylov_final_state_local_energy_export(&s->final_state,
                                                         &energy) ==
                 MVMC_SCALED_EXPORT_OK &&
       cabs(energy - expect->energy) <= tolerance * (1.0 + cabs(expect->energy));
  return ok;
}

/* Rejected or failed: the current state is unchanged field by field. */
static int unchanged(const MVMCGCFinalStepTrace *entry) {
  const MVMCKrylovFinalStateSnapshot *a = &entry->before, *b = &entry->after;
  return entry->before_words[0] == entry->after_words[0] &&
         a->accepted_generation == b->accepted_generation &&
         a->integrity_hash == b->integrity_hash &&
         a->configuration_hash == b->configuration_hash &&
         memcmp(&a->final_state.amplitude, &b->final_state.amplitude,
                sizeof(a->final_state.amplitude)) == 0 &&
         memcmp(&a->final_state.local_energy, &b->final_state.local_energy,
                sizeof(a->final_state.local_energy)) == 0 &&
         memcmp(a->krylov.value, b->krylov.value, sizeof(a->krylov.value)) ==
             0;
}

static Expect expect_of(unsigned words, uint64_t generation, int accepted,
                        uint64_t calls, uint64_t accepted_calls,
                        uint64_t exact_zero_calls, int refined,
                        const Task5Scaled *table, size_t count,
                        double complex energy) {
  Expect e;
  int n = 0;
  unsigned w = words;
  while (w) {
    n += (int)(w & 1u);
    w >>= 1;
  }
  e.words = words;
  e.generation = generation;
  e.electron_count = n;
  e.electron_count_square = (uint64_t)n * (uint64_t)n;
  e.accepted = accepted;
  e.calls = calls;
  e.accepted_calls = accepted_calls;
  e.exact_zero_calls = exact_zero_calls;
  e.refined = refined;
  e.table = table;
  e.table_count = count;
  e.energy = energy;
  return e;
}

static void report(const char *label, const Run *run) {
  size_t k;
  if (world_rank() != 0) return;
  fprintf(stderr, "%s: status=%s written=%zu\n", label,
          mvmc_krylov_status_string(run->status), run->written);
  for (k = 0; k < run->written && k < MAX_STEPS; ++k) {
    const MVMCGCFinalStepTrace *e = run->trace + k;
    fprintf(stderr,
            "  step %zu: %s acc=%d ref=%d ez=%d calls=%llu/%llu/%llu "
            "%llu->%llu gen=%llu N=%d prec=%u nodes=%llu\n",
            k + 1, mvmc_krylov_status_string(e->status), e->accepted,
            e->used_refinement, e->exact_zero_proposal,
            (unsigned long long)e->refinement_calls,
            (unsigned long long)e->refinement_accepted,
            (unsigned long long)e->refinement_exact_zero,
            (unsigned long long)e->before_words[0],
            (unsigned long long)e->after_words[0],
            (unsigned long long)e->after.accepted_generation,
            e->electron_count, e->precision,
            (unsigned long long)e->node_calls);
  }
}

static void test_r1_variants(const Task5Case *source, const Chain *chain,
                             const char *label) {
  Task5Built built;
  Run run;
  const Step r1[3] = {{5, 0.125, MVMC_GC_FINAL_FORCE_NUMERIC_ZERO},
                      {6, 0.75, MVMC_GC_FINAL_FORCE_NONE},
                      {0, 0.0, MVMC_GC_FINAL_FORCE_NONE}};
  Step le[3];
  const Task5Scaled *table = source->scaled;
  const size_t count = source->scaled_count;
  int variant;
  CHECK(task5_build(source, &built), "build");
  memcpy(le, r1, sizeof(le));
  le[0].force = MVMC_GC_FINAL_FORCE_LOCAL_ENERGY;
  for (variant = 0; variant < 2; ++variant) {
    run_script(&built, chain, built.coefficient[1], 6, variant ? le : r1, 3,
               MAX_STEPS, &run);
    if (run.status != MVMC_KRYLOV_STATUS_OK) report(label, &run);
    CHECK(run.status == MVMC_KRYLOV_STATUS_OK && run.written == 3,
          "R1 script completes");
    if (run.status != MVMC_KRYLOV_STATUS_OK) continue;
    {
      const Expect e1 = expect_of(5, 1, 1, 1, 1, 0, 1, table, count, 2.0);
      const Expect e2 = expect_of(6, 2, 1, 1, 1, 0, 0, table, count, 0.5);
      CHECK(step_matches(run.trace + 0, &e1) && run.trace[0].used_refinement &&
                !run.trace[0].exact_zero_proposal &&
                run.trace[0].precision >= 128 && run.trace[0].node_calls > 0,
            "R1 step 1: forced entry, refined accept u=1/8 -> 5");
      CHECK(step_matches(run.trace + 1, &e2) && !run.trace[1].used_refinement,
            "R1 step 2: ordinary accept from the refined snapshot -> 6");
      CHECK(run.trace[2].status == MVMC_KRYLOV_STATUS_OK &&
                !run.trace[2].accepted && run.trace[2].used_refinement &&
                run.trace[2].exact_zero_proposal && unchanged(run.trace + 2) &&
                run.trace[2].after_words[0] == 6 &&
                run.trace[2].after.accepted_generation == 2 &&
                run.trace[2].refinement_calls == 2 &&
                run.trace[2].refinement_accepted == 1 &&
                run.trace[2].refinement_exact_zero == 1,
            "R1 step 3: Phi(0)=0 natural exact-zero entry, rejected (2/1/1)");
    }
    CHECK(run.counters.steps[MVMC_POWER_LANCZOS_GC_STAGE_FINAL_SAMPLE] == 3 &&
              run.counters.generated[MVMC_POWER_LANCZOS_GC_STAGE_FINAL_SAMPLE] ==
                  3 &&
              run.trace[2].proposal_copies == 3 &&
              run.trace[2].uniform_copies == 3 && run.trace[2].rng_draws == 0,
          "harness counters: one step per entry, no RNG");
  }
  /* Boundary u = 1/4 (strict) and u = 1/2: rejected, current kept. */
  {
    const double uniforms[2] = {0.25, 0.5};
    int u;
    for (u = 0; u < 2; ++u) {
      const Step step = {5, uniforms[u], MVMC_GC_FINAL_FORCE_NUMERIC_ZERO};
      run_script(&built, chain, built.coefficient[1], 6, &step, 1, MAX_STEPS,
                 &run);
      CHECK(run.status == MVMC_KRYLOV_STATUS_OK && run.written == 1 &&
                !run.trace[0].accepted && run.trace[0].used_refinement &&
                !run.trace[0].exact_zero_proposal && unchanged(run.trace) &&
                run.trace[0].after_words[0] == 6 &&
                run.trace[0].after.accepted_generation == 0 &&
                run.trace[0].refinement_calls == 1 &&
                run.trace[0].refinement_accepted == 0 &&
                run.trace[0].refinement_exact_zero == 0,
            "R1 boundary/rejection: ratio^2 = 1/4 is not > u (1/0/0)");
    }
  }
  /* R1C: c1 = i/2, complex local energy. */
  if (source == &task5_r1) {
    const Step step = {5, 0.125, MVMC_GC_FINAL_FORCE_NUMERIC_ZERO};
    const Expect e = expect_of(
        5, 1, 1, 1, 1, 0, 1, task5_r1c_scaled,
        sizeof(task5_r1c_scaled) / sizeof(task5_r1c_scaled[0]), -2.0 * I);
    run_script(&built, chain, task5_complex(task5_r1c_c1), 6, &step, 1,
               MAX_STEPS, &run);
    CHECK(run.status == MVMC_KRYLOV_STATUS_OK && step_matches(run.trace, &e),
          "R1C: refined accept with E_local(5) = -2i (1/1/0)");
  }
  task5_free(&built);
}

static void test_e1_variants(const Chain *chain) {
  Task5Built built;
  Run run;
  const Step mixed[3] = {{53, 0.125, MVMC_GC_FINAL_FORCE_NUMERIC_ZERO},
                         {18, 0.0, MVMC_GC_FINAL_FORCE_NONE},
                         {126, 0.0, MVMC_GC_FINAL_FORCE_NONE}};
  Step all[3];
  const unsigned path[3] = {53, 18, 126};
  int variant;
  size_t k;
  CHECK(task5_build(&task5_e1, &built), "E1 build");
  memcpy(all, mixed, sizeof(all));
  all[1].force = all[2].force = MVMC_GC_FINAL_FORCE_NUMERIC_ZERO;
  for (variant = 0; variant < 2; ++variant) {
    run_script(&built, chain, built.coefficient[1], 60,
               variant ? all : mixed, 3, MAX_STEPS, &run);
    if (run.status != MVMC_KRYLOV_STATUS_OK) report("E1", &run);
    CHECK(run.status == MVMC_KRYLOV_STATUS_OK && run.written == 3,
          "E1 script completes");
    for (k = 0; run.status == MVMC_KRYLOV_STATUS_OK && k < 3; ++k) {
      const int refined = variant || k == 0;
      const uint64_t calls = variant ? k + 1 : 1;
      const Expect e = expect_of(path[k], k + 1, 1, calls, calls, 0, refined,
                                 task5_e1.scaled, task5_e1.scaled_count,
                                 (double)(k == 0 ? 4 : k == 1 ? 2 : 6));
      CHECK(step_matches(run.trace + k, &e) &&
                run.trace[k].used_refinement == refined,
            variant ? "E1 all refined: N 4->2->6, N^2, energy (3/3/0)"
                    : "E1 mixed: refined then ordinary, N 4->2->6 (1/1/0)");
      if (refined && world_rank() == 0) {
        printf("E1 %s step %zu: precision=%u node_calls=%llu\n",
               variant ? "all-refined" : "mixed", k + 1,
               run.trace[k].precision,
               (unsigned long long)run.trace[k].node_calls);
      }
    }
  }
  {
    const Step reject = {53, 0.75, MVMC_GC_FINAL_FORCE_NUMERIC_ZERO};
    run_script(&built, chain, built.coefficient[1], 60, &reject, 1, MAX_STEPS,
               &run);
    CHECK(run.status == MVMC_KRYLOV_STATUS_OK && !run.trace[0].accepted &&
              run.trace[0].used_refinement && unchanged(run.trace) &&
              run.trace[0].after_words[0] == 60,
          "E1 u=3/4 rejected by the exponential interval");
  }
  task5_free(&built);
}

static void test_s2_g2(const Chain *chain) {
  Task5Built built;
  Run run;
  const Step step = {53, 0.5, MVMC_GC_FINAL_FORCE_NONE};
  const MVMCGCKrylovAmplitudeFault mutations[2] = {
      MVMC_GC_KRYLOV_AMPLITUDE_FAULT_ORDINARY_UNMASKED,
      MVMC_GC_KRYLOV_AMPLITUDE_FAULT_REFINED_UNMASKED};
  int m;
  CHECK(task5_build(&task5_s2_g2, &built), "S2_G2 build");
  run_script(&built, chain, built.coefficient[1], 60, &step, 1, MAX_STEPS,
             &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_OK && !run.trace[0].accepted &&
            run.trace[0].used_refinement && run.trace[0].exact_zero_proposal &&
            unchanged(run.trace) && run.trace[0].after_words[0] == 60 &&
            run.trace[0].refinement_calls == 1 &&
            run.trace[0].refinement_accepted == 0 &&
            run.trace[0].refinement_exact_zero == 1,
        "S2_G2: ordinary exact zero, refined zero, current 60 (1/0/1)");
  for (m = 0; m < 2; ++m) {
    mvmc_gc_krylov_amplitude_testing_set_fault(mutations[m], world_rank());
    run_script(&built, chain, built.coefficient[1], 60, &step, 1, MAX_STEPS,
               &run);
    mvmc_gc_krylov_amplitude_testing_set_fault(
        MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONE, -1);
    CHECK(run.status == MVMC_KRYLOV_STATUS_OK && run.trace[0].accepted &&
              run.trace[0].after_words[0] == 53 &&
              run.trace[0].used_refinement == m,
          m ? "refined-only mask removal flips the decision to accept"
            : "ordinary-only mask removal flips the decision to accept");
  }
  task5_free(&built);
}

/* The comparator rejects each mutation of a correct trace entry. */
static void test_comparator_mutations(const Chain *chain) {
  Task5Built built;
  Run run;
  const Step step = {5, 0.125, MVMC_GC_FINAL_FORCE_NUMERIC_ZERO};
  Expect e;
  MVMCGCFinalStepTrace mutated;
  CHECK(task5_build(&task5_r1, &built), "R1 build");
  e = expect_of(5, 1, 1, 1, 1, 0, 1, task5_r1.scaled, task5_r1.scaled_count,
                2.0);
  run_script(&built, chain, built.coefficient[1], 6, &step, 1, MAX_STEPS,
             &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_OK && step_matches(run.trace, &e),
        "reference entry");
  /* Each mutant gets a consistent integrity hash, so only the rule under
   * test can reject it. */
#define REHASH(t) \
  ((t).after.integrity_hash = \
       mvmc_krylov_final_state_snapshot_integrity_hash(&(t).after))
  mutated = run.trace[0];
  REHASH(mutated);
  CHECK(step_matches(&mutated, &e), "rehashed copy still matches");
  mutated = run.trace[0];
  mutated.after.final_state.amplitude.phase =
      -mutated.after.final_state.amplitude.phase;
  REHASH(mutated);
  CHECK(!step_matches(&mutated, &e), "phase sign flip rejected");
  mutated = run.trace[0];
  mutated.after.krylov.value[1].log_abs += 1e-10;
  REHASH(mutated);
  CHECK(!step_matches(&mutated, &e), "log_abs shift of 1e-10 rejected");
  mutated = run.trace[0];
  mutated.after.accepted_generation = 0;
  REHASH(mutated);
  CHECK(!step_matches(&mutated, &e), "generation not updated rejected");
  mutated = run.trace[0];
  mutated.electron_count_square = 2;
  CHECK(!step_matches(&mutated, &e), "N^2 miscount rejected");
  mutated = run.trace[0];
  mutated.after.krylov.value[0].state = MVMC_SCALED_COMPLEX_NUMERIC_ZERO;
  REHASH(mutated);
  CHECK(!step_matches(&mutated, &e),
        "NUMERIC_ZERO in place of EXACT_ZERO rejected");
  mutated = run.trace[0];
  mutated.after.final_state.local_energy.log_abs += 1e-12;
  REHASH(mutated);
  CHECK(!step_matches(&mutated, &e), "E_local beyond 1e-14 rejected");
#undef REHASH
  task5_free(&built);
}

/* RR-1 and invalid scripts: agreed before any step, trace untouched. */
static void test_invalid(const Chain *chain) {
  Task5Built built;
  Run run, sentinel;
  const Step step = {5, 0.125, MVMC_GC_FINAL_FORCE_NUMERIC_ZERO};
  const MVMCPowerLanczosGCFault forbidden[] = {
      MVMC_POWER_LANCZOS_GC_FAULT_FINAL_NUMERIC,
      MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_NUMERIC,
      MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_EXHAUSTED};
  size_t k;
  CHECK(task5_build(&task5_r1, &built), "R1 build");
  memset(&sentinel, 0x5a, sizeof(sentinel));
  for (k = 0; k < sizeof(forbidden) / sizeof(forbidden[0]); ++k) {
    mvmc_power_lanczos_gc_testing_set_fault(forbidden[k]);
    run_script(&built, chain, built.coefficient[1], 6, &step, 1, MAX_STEPS,
               &run);
    mvmc_power_lanczos_gc_testing_set_fault(MVMC_POWER_LANCZOS_GC_FAULT_NONE);
    CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              run.written == 0 &&
              memcmp(&run.trace[0].before, &sentinel.trace[0].before,
                     sizeof(run.trace[0].before)) == 0 &&
              run.counters.steps[MVMC_POWER_LANCZOS_GC_STAGE_FINAL_SAMPLE] ==
                  0 &&
              all_ranks_equal_int((int)run.status),
          "a forbidden global fault is rejected before any step");
  }
  if (world_size() > 1) {
    mvmc_power_lanczos_gc_testing_set_fault(
        world_rank() == 1 ? MVMC_POWER_LANCZOS_GC_FAULT_REFINEMENT_DECISION
                          : MVMC_POWER_LANCZOS_GC_FAULT_NONE);
    run_script(&built, chain, built.coefficient[1], 6, &step, 1, MAX_STEPS,
               &run);
    mvmc_power_lanczos_gc_testing_set_fault(MVMC_POWER_LANCZOS_GC_FAULT_NONE);
    CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              run.written == 0 && all_ranks_equal_int((int)run.status),
          "a rank-mismatched fault mode is rejected");
  }
  run_script(&built, chain, built.coefficient[1], 6, &step, 1, 0, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && run.written == 0,
        "capacity below the step count rejected");
  {
    const double bad_uniforms[3] = {1.0, -0.125, NAN};
    for (k = 0; k < 3; ++k) {
      Step bad = step;
      bad.uniform = bad_uniforms[k];
      run_script(&built, chain, built.coefficient[1], 6, &bad, 1, MAX_STEPS,
                 &run);
      CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                run.written == 0,
            "uniform outside [0, 1) or NaN rejected");
    }
  }
  {
    Step bad = step;
    bad.force = (MVMCGCFinalForce)7;
    run_script(&built, chain, built.coefficient[1], 6, &bad, 1, MAX_STEPS,
               &run);
    CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "unknown force rejected");
  }
  if (world_size() > 1) {
    Step bad = step;
    if (world_rank() == world_size() - 1) bad.uniform = 0.25;
    run_script(&built, chain, built.coefficient[1], 6, &bad, 1, MAX_STEPS,
               &run);
    CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              run.written == 0 && all_ranks_equal_int((int)run.status),
          "a script differing on one rank is rejected by the world hash");
    bad = step;
    if (world_rank() == 1) bad.force = MVMC_GC_FINAL_FORCE_NONE;
    run_script(&built, chain, built.coefficient[1], 6, &bad, 1, MAX_STEPS,
               &run);
    CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              all_ranks_equal_int((int)run.status),
          "a force differing on one rank is rejected");
  }
  {
    MVMCPowerLanczosGCView view;
    MVMCGCFinalScript script;
    MVMCGCFinalScriptStep entry;
    MVMCGCFinalStepTrace trace;
    uint64_t initial = 6, proposal = 5, shared = 0;
    size_t written = 77;
    task5_view(&built, world_rank(), chain->chain_rank, chain->chain_size,
               &view);
    memset(&script, 0, sizeof(script));
    script.gc_view = &view;
    script.world_communicator = world();
    script.chain_communicator = chain->chain;
    script.coefficient[0] = 1.0;
    script.coefficient[1] = 0.5;
    script.initial_words = &initial;
    script.word_count = 1;
    script.step_count = 1;
    script.steps = &entry;
    entry.proposal_words = &proposal;
    entry.uniform = 0.125;
    entry.force = MVMC_GC_FINAL_FORCE_NONE;
    trace.before_words = &shared;
    trace.after_words = &shared;
    CHECK(mvmc_power_lanczos_gc_testing_final_steps(&script, &trace, 1,
                                                    &written) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              written == 0,
          "aliased trace buffers rejected");
    trace.after_words = &initial;
    CHECK(mvmc_power_lanczos_gc_testing_final_steps(&script, &trace, 1,
                                                    &written) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "a trace buffer aliasing the script rejected");
    trace.after_words = &proposal;
    script.word_count = 2;
    CHECK(mvmc_power_lanczos_gc_testing_final_steps(&script, &trace, 1,
                                                    &written) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "wrong word count rejected");
  }
  task5_free(&built);
}

/* Allowed RR-1 faults on chain rank 1: disagreement, no hang. */
static void test_refinement_faults(int mode) {
  Chain chain = make_chain(mode);
  Task5Built built;
  Run run;
  const Step step = {5, 0.125, MVMC_GC_FINAL_FORCE_NUMERIC_ZERO};
  const MVMCPowerLanczosGCFault faults[2] = {
      MVMC_POWER_LANCZOS_GC_FAULT_REFINEMENT_DECISION,
      MVMC_POWER_LANCZOS_GC_FAULT_REFINEMENT_STATUS};
  int f;
  if (chain.chain_size < 2) {
    free_chain(&chain);
    return;
  }
  CHECK(task5_build(&task5_r1, &built), "R1 build");
  for (f = 0; f < 2; ++f) {
    mvmc_power_lanczos_gc_testing_set_fault(faults[f]);
    run_script(&built, &chain, built.coefficient[1], 6, &step, 1, MAX_STEPS,
               &run);
    mvmc_power_lanczos_gc_testing_set_fault(MVMC_POWER_LANCZOS_GC_FAULT_NONE);
    CHECK(run.status == MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE &&
              run.written == 1 && unchanged(run.trace) &&
              run.trace[0].after_words[0] == 6 &&
              all_ranks_equal_int((int)run.status),
          "chain rank 1 refinement fault fails every rank, current kept");
  }
  task5_free(&built);
  free_chain(&chain);
}

/* Insufficient refinement limits: RESOURCE_LIMIT, current kept, the same
 * script succeeds again after the defaults are restored. */
static void test_resource_limits(const Chain *chain) {
  Task5Built built;
  Run run;
  const Step step = {53, 0.125, MVMC_GC_FINAL_FORCE_NUMERIC_ZERO};
  const MVMCKrylovRefinementLimits defaults = {16384, 8192, 8192, 65536,
                                               UINT64_C(2000000), 4096};
  MVMCKrylovRefinementLimits tight[5];
  const char *names[5] = {"matrix", "cache", "terms", "bits", "node calls"};
  int k;
  CHECK(task5_build(&task5_e1, &built), "E1 build");
  for (k = 0; k < 5; ++k) tight[k] = defaults;
  tight[0].max_matrix_entries = 1;
  tight[1].max_cache_entries = 1;
  tight[2].max_expression_terms = 1;
  tight[3].max_rational_bits = 16;
  tight[4].max_node_calls = 1;
  for (k = 0; k < 5; ++k) {
    mvmc_gc_krylov_amplitude_testing_set_refinement_limits(tight + k);
    run_script(&built, chain, built.coefficient[1], 60, &step, 1, MAX_STEPS,
               &run);
    mvmc_gc_krylov_amplitude_testing_set_refinement_limits(NULL);
    if (run.status != MVMC_KRYLOV_STATUS_RESOURCE_LIMIT && world_rank() == 0)
      fprintf(stderr, "limit %s: %s\n", names[k],
              mvmc_krylov_status_string(run.status));
    CHECK(run.status == MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
              run.written == 1 &&
              run.trace[0].status == MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
              unchanged(run.trace) && run.trace[0].after_words[0] == 60 &&
              run.trace[0].after.accepted_generation == 0,
          "an insufficient refinement limit is RESOURCE_LIMIT, current kept");
  }
  run_script(&built, chain, built.coefficient[1], 60, &step, 1, MAX_STEPS,
             &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_OK && run.trace[0].accepted,
        "the default limits resolve the same step");
  task5_free(&built);
}

/* The refinement API directly (no harness): decisions and the refined
 * Krylov/Phi/E_local of an accepted proposal against the 16-eps rule. */
static void test_direct_api(void) {
  static const struct {
    const Task5Case *source;
    int complex_c1;
    unsigned current, proposal;
    double uniform;
    int accepted;
    double complex energy;
  } cases[] = {{&task5_r1, 0, 6, 5, 0.125, 1, 2.0},
               {&task5_r1, 0, 6, 5, 0.25, 0, 0.0},
               {&task5_r1, 0, 6, 5, 0.5, 0, 0.0},
               {&task5_r1, 1, 6, 5, 0.125, 1, -2.0 * I},
               {&task5_e1, 0, 60, 53, 0.125, 1, 4.0},
               {&task5_e1, 0, 60, 53, 0.75, 0, 0.0},
               {&task5_r1_g2, 0, 6, 5, 0.125, 1, 2.0}};
  size_t c;
  for (c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
    Task5Built built;
    MVMCGCKrylovAmplitudeLayout layout;
    MVMCGCKrylovAmplitudeWorkspace *workspace = NULL;
    MVMCKrylovFinalStatePolicy policy;
    MVMCKrylovRefinementResult refined;
    double complex coefficient[2];
    const uint64_t current = cases[c].current, proposal = cases[c].proposal;
    const Task5Scaled *table;
    size_t count;
    CHECK(task5_build(cases[c].source, &built), "build");
    coefficient[0] = built.coefficient[0];
    coefficient[1] = cases[c].complex_c1 ? task5_complex(task5_r1c_c1)
                                         : built.coefficient[1];
    table = cases[c].complex_c1 ? task5_r1c_scaled : cases[c].source->scaled;
    count = cases[c].complex_c1
                ? sizeof(task5_r1c_scaled) / sizeof(task5_r1c_scaled[0])
                : cases[c].source->scaled_count;
    task5_layout(&built, world(), world_rank(), world_size(), &layout);
    CHECK(mvmc_gc_krylov_amplitude_workspace_create(
              &layout, built.slater, built.weight, &workspace) ==
                  MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_final_state_policy_create(
                  1, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC, 53, coefficient,
                  2, &policy) == MVMC_KRYLOV_STATUS_OK,
          "direct API setup");
    if (workspace != NULL) {
      memset(&refined, 0, sizeof(refined));
      CHECK(mvmc_gc_krylov_refine_final(workspace, &built.model, &policy,
                                        &current, &proposal, 1,
                                        cases[c].uniform, &refined) ==
                    MVMC_KRYLOV_STATUS_OK &&
                refined.valid && refined.accepted == cases[c].accepted &&
                !refined.exact_zero_proposal,
            "direct refinement decision (u=1/4 strict)");
      if (refined.accepted && table != NULL && count != 0) {
        double complex energy = NAN;
        int v, ok = 1;
        for (v = 0; v < 3; ++v)
          ok = ok && value_matches(refined.krylov.value + v, table, count,
                                   proposal, v, 16.0);
        ok = ok && value_matches(&refined.evaluation.amplitude, table, count,
                                 proposal, 3, 16.0) &&
             value_matches(&refined.evaluation.hamiltonian_amplitude, table,
                           count, proposal, 4, 16.0);
        CHECK(ok, "refined Psi/H Psi/H^2 Psi/Phi/H Phi within 16 eps");
        CHECK(mvmc_krylov_final_state_local_energy_export(
                  &refined.evaluation, &energy) == MVMC_SCALED_EXPORT_OK &&
                  cabs(energy - cases[c].energy) <=
                      1e-14 * (1.0 + cabs(cases[c].energy)),
              "refined E_local within 1e-14");
      }
      mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
    }
    task5_free(&built);
  }
}

int main(int argc, char **argv) {
  int total = 0;
  Chain chain;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  test_direct_api();
  chain = make_chain(0);
  test_r1_variants(&task5_r1, &chain, "R1");
  /* R1_G2: 2 QPs, an empty QP rank at np >= 3, same decisions. */
  {
    Task5Built built;
    Run run;
    const Step r1[3] = {{5, 0.125, MVMC_GC_FINAL_FORCE_NUMERIC_ZERO},
                        {6, 0.75, MVMC_GC_FINAL_FORCE_NONE},
                        {0, 0.0, MVMC_GC_FINAL_FORCE_NONE}};
    CHECK(task5_build(&task5_r1_g2, &built), "R1_G2 build");
    run_script(&built, &chain, built.coefficient[1], 6, r1, 3, MAX_STEPS,
               &run);
    if (run.status != MVMC_KRYLOV_STATUS_OK) report("R1_G2", &run);
    CHECK(run.status == MVMC_KRYLOV_STATUS_OK && run.written == 3 &&
              run.trace[0].accepted && run.trace[0].after_words[0] == 5 &&
              run.trace[1].accepted && run.trace[1].after_words[0] == 6 &&
              !run.trace[2].accepted && run.trace[2].refinement_calls == 2 &&
              run.trace[2].refinement_accepted == 1 &&
              run.trace[2].refinement_exact_zero == 1,
          "R1_G2 (empty QP rank): ratio^2 = 1/4 accepted at u=1/8 (2/1/1)");
    task5_free(&built);
  }
  test_e1_variants(&chain);
  test_s2_g2(&chain);
  test_comparator_mutations(&chain);
  test_invalid(&chain);
  test_resource_limits(&chain);
  free_chain(&chain);
  test_refinement_faults(0);
  if (world_size() >= 4) test_refinement_faults(2);
#ifdef _mpi_use
  MPI_Allreduce(&failure_count, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  total = failure_count;
#endif
  if (world_rank() == 0) {
    printf(total == 0 ? "gc_lanczos_final_steps_unit: PASS\n"
                      : "gc_lanczos_final_steps_unit: FAIL\n");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return total == 0 ? 0 : 1;
}
