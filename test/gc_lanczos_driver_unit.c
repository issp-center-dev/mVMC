/*
 * GC first-step Lanczos driver contract: every proposal stage (pilot,
 * coefficient warm-up/sampling, final bootstrap/warm-up/sampling) uses the
 * root-drawn symmetric GC kernel and never the fixed-sector helpers;
 * generation counts on the chain root and step counts on every rank;
 * reason-aware final bootstrap and its faults; zero H rejection; chain
 * topologies (serial, np2 with chain size 1 and 2, np3 with chains 2+1).
 */
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

#define NSITE 2
#define ORB 4

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

typedef struct {
  MVMCClassicPfaffianCommunicator chain;
  int chain_rank;
  int chain_size;
  int owns;
} Chain;

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
    MPI_Comm_split(MPI_COMM_WORLD, world_rank() < 2 ? 0 : 1, world_rank(),
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

static MVMCClassicPfaffianCommunicator world(void) {
#ifdef _mpi_use
  return MPI_COMM_WORLD;
#else
  return 0;
#endif
}

typedef struct {
  int transfer_rows[6][4];
  int *transfer_index[6];
  double complex transfer_value[6];
  int intra_index[2];
  double intra_value[2];
  double complex translation[2];
  double complex weights[2];
  double complex slater[2 * ORB * ORB];
  MVMCPowerLanczosGCView view;
} Fixture;

/* F2 anti-parallel pairing, PBC k=0 momentum (2 QPs), t=1, mu=1/4, U=3/2. */
static void make_fixture(Fixture *f, const Chain *chain, int zero_h) {
  static const double complex F[NSITE][NSITE] = {
      {1.0, (1.0 + 2.0 * I) / 8.0}, {(3.0 - 1.0 * I) / 16.0, 0.75}};
  const int rows[6][4] = {{0, 0, 1, 0}, {1, 0, 0, 0}, {0, 1, 1, 1},
                          {1, 1, 0, 1}, {0, 0, 0, 0}, {1, 1, 1, 1}};
  int i, j, r;
  memset(f, 0, sizeof(*f));
  for (r = 0; r < 6; ++r) {
    memcpy(f->transfer_rows[r], rows[r], sizeof(rows[r]));
    f->transfer_index[r] = f->transfer_rows[r];
    f->transfer_value[r] = r < 4 ? 1.0 : 0.25;
  }
  f->intra_index[0] = 0;
  f->intra_index[1] = 1;
  f->intra_value[0] = f->intra_value[1] = 1.5;
  for (r = 0; r < 2; ++r) {
    double complex *s = f->slater + r * ORB * ORB;
    for (i = 0; i < NSITE; ++i) {
      for (j = 0; j < NSITE; ++j) {
        /* Translation by r: F_r[i][j] = F[i+r][j+r] (PBC, no sign). */
        const double complex value = F[(i + r) % NSITE][(j + r) % NSITE];
        s[i * ORB + (j + NSITE)] = value;
        s[(j + NSITE) * ORB + i] = -value;
      }
    }
    f->translation[r] = 0.5;
    f->weights[r] = 0.5;
  }
  f->view.data.site_count = NSITE;
  f->view.data.arithmetic = MVMC_POWER_LANCZOS_CLASSIC_COMPLEX;
  if (!zero_h) {
    f->view.data.transfer_count = 6;
    f->view.data.coulomb_intra_count = 2;
  }
  if (world_rank() == 0) {
    f->view.data.transfer_indices = f->transfer_index;
    f->view.data.transfer_parameters = f->transfer_value;
    f->view.data.coulomb_intra_indices = f->intra_index;
    f->view.data.coulomb_intra_parameters = f->intra_value;
  }
  f->view.data.qp_total = 2;
  f->view.data.qp_start = 2 * chain->chain_rank / chain->chain_size;
  f->view.data.qp_end = 2 * (chain->chain_rank + 1) / chain->chain_size;
  f->view.data.qp_weights = f->weights;
  f->view.data.slater_complex = f->slater;
  f->view.domain = MVMC_KRYLOV_DOMAIN_GC_BALANCED;
  f->view.translation_count = 2;
  f->view.translation_weights = f->translation;
}

typedef struct {
  MVMCKrylovStatus status;
  MVMCPowerLanczosStabilizedResult result;
  MVMCPowerLanczosGCDiagnostics diagnostics;
  MVMCPowerLanczosGCTestingCounters counters;
  uint64_t sectors[ORB + 1];
} Run;

static void run_gc(const Chain *chain, int zero_h, int initial_n,
                   size_t warm_up, size_t samples, size_t interval,
                   MVMCPowerLanczosGCFault fault, Run *run) {
  Fixture fixture;
  MVMCPowerLanczosGCInput input;
  make_fixture(&fixture, chain, zero_h);
  memset(run, 0, sizeof(*run));
  memset(&input, 0, sizeof(input));
  input.gc_view = &fixture.view;
  input.world_communicator = world();
  input.chain_communicator = chain->chain;
  input.seed = UINT64_C(17);
  input.warm_up = warm_up;
  input.sample_count = samples;
  input.interval = interval;
  input.initial_electron_count = initial_n;
  run->diagnostics.final_sector_samples = run->sectors;
  run->diagnostics.final_sector_capacity = ORB + 1;
  mvmc_power_lanczos_gc_testing_set_fault(fault);
  run->status = mvmc_power_lanczos_gc_run_diagnosed(&input, &run->result,
                                                    &run->diagnostics);
  mvmc_power_lanczos_gc_testing_set_fault(MVMC_POWER_LANCZOS_GC_FAULT_NONE);
  mvmc_power_lanczos_gc_testing_counters(&run->counters);
}

static uint64_t chain_count_for(int mode) {
  const int size = world_size();
  if (mode == 0) return 1;
  if (mode == 1) return (uint64_t)size;
  return size >= 3 ? 2 : 1;
}

static int all_ranks_equal_double(double value) {
#ifdef _mpi_use
  double minimum = value, maximum = value;
  MPI_Allreduce(&value, &minimum, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&value, &maximum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return minimum == maximum;
#else
  (void)value;
  return 1;
#endif
}

static void test_stage_counts(int mode) {
  const size_t warm_up = 8, samples = 256, interval = 2;
  Chain chain = make_chain(mode);
  const uint64_t chains = chain_count_for(mode);
  Run run, again;
  const MVMCPowerLanczosGCTestingCounters *c = &run.counters;
  const int root = chain.chain_rank == 0;
  run_gc(&chain, 0, 2, warm_up, samples, interval,
         MVMC_POWER_LANCZOS_GC_FAULT_NONE, &run);
  if (run.status != MVMC_KRYLOV_STATUS_OK && world_rank() == 0) {
    fprintf(stderr, "GC run status %s, failure %d stage %s status %s N=%d\n",
            mvmc_krylov_status_string(run.status),
            run.diagnostics.failure_recorded,
            mvmc_power_lanczos_gc_stage_name(run.diagnostics.failure_stage),
            mvmc_krylov_status_string(run.diagnostics.failure_status),
            run.diagnostics.failure_electron_count);
    {
      int s;
      for (s = 0; s < MVMC_POWER_LANCZOS_GC_STAGE_COUNT; ++s) {
        fprintf(stderr, "  stage %-20s generated %llu steps %llu\n",
                mvmc_power_lanczos_gc_stage_name((MVMCPowerLanczosGCStage)s),
                (unsigned long long)run.counters.generated[s],
                (unsigned long long)run.counters.steps[s]);
      }
    }
  }
  CHECK(run.status == MVMC_KRYLOV_STATUS_OK && run.result.valid &&
            isfinite(run.result.energy) && isfinite(run.result.variance) &&
            run.result.power_step == 1 &&
            run.result.sampling_chains == chains &&
            run.result.coefficient_samples == chains * samples &&
            run.result.coefficient_energy_available,
        "GC run succeeds with the expected topology");
  CHECK(c->fixed_sector_calls == 0, "GC never calls fixed-sector proposals");
  CHECK(c->generated[MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_INITIAL] ==
                (root ? 1u : 0u) &&
            c->generated[MVMC_POWER_LANCZOS_GC_STAGE_PILOT_WARM_UP] ==
                (root ? warm_up : 0u) &&
            c->generated[MVMC_POWER_LANCZOS_GC_STAGE_PILOT_SAMPLE] ==
                (root ? samples : 0u) &&
            c->generated[MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_WARM_UP] ==
                (root ? warm_up : 0u) &&
            c->generated[MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_SAMPLE] ==
                (root ? samples * interval : 0u) &&
            c->generated[MVMC_POWER_LANCZOS_GC_STAGE_FINAL_WARM_UP] ==
                (root ? warm_up : 0u) &&
            c->generated[MVMC_POWER_LANCZOS_GC_STAGE_FINAL_SAMPLE] ==
                (root ? samples * interval : 0u),
        "candidates are generated on the chain root only, per stage");
  CHECK(c->steps[MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_INITIAL] == 0 &&
            c->steps[MVMC_POWER_LANCZOS_GC_STAGE_PILOT_WARM_UP] == warm_up &&
            c->steps[MVMC_POWER_LANCZOS_GC_STAGE_PILOT_SAMPLE] == samples &&
            c->steps[MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_WARM_UP] ==
                warm_up &&
            c->steps[MVMC_POWER_LANCZOS_GC_STAGE_COEFFICIENT_SAMPLE] ==
                samples * interval &&
            c->steps[MVMC_POWER_LANCZOS_GC_STAGE_FINAL_BOOTSTRAP] >= 1 &&
            c->steps[MVMC_POWER_LANCZOS_GC_STAGE_FINAL_WARM_UP] == warm_up &&
            c->steps[MVMC_POWER_LANCZOS_GC_STAGE_FINAL_SAMPLE] ==
                samples * interval,
        "every rank enters each candidate step");
  CHECK(c->generated[MVMC_POWER_LANCZOS_GC_STAGE_FINAL_BOOTSTRAP] ==
            (root ? c->steps[MVMC_POWER_LANCZOS_GC_STAGE_FINAL_BOOTSTRAP]
                  : 0u),
        "bootstrap candidates are drawn on the root");
  if (world_rank() == 0) {
    uint64_t total = 0;
    size_t n;
    for (n = 0; n <= ORB; ++n) total += run.sectors[n];
    CHECK(run.diagnostics.valid && run.diagnostics.chain_count == chains &&
              run.diagnostics.coefficient_available &&
              run.diagnostics.coefficient_energy_available &&
              run.diagnostics.total_bootstrap_attempts >= chains &&
              !run.diagnostics.failure_recorded,
          "world-root diagnostics");
    CHECK(total == chains * samples, "final N histogram counts every sample");
    {
      const MVMCPowerLanczosGCNumberResult *number = &run.result.gc_number;
      double sum = 0.0, square_sum = 0.0;
      double block_sum = 0.0, block_square_sum = 0.0;
      double se_sum = 0.0, square_se_sum = 0.0;
      size_t b;
      for (n = 0; n <= ORB; ++n) {
        sum += (double)n * (double)run.sectors[n];
        square_sum += (double)n * (double)n * (double)run.sectors[n];
      }
      CHECK(number->valid && number->mean == sum / (double)total &&
                number->mean_square == square_sum / (double)total &&
                fabs(number->variance - (square_sum / (double)total -
                     number->mean * number->mean)) < 1e-12,
            "final number moments equal the independently accumulated histogram");
      for (b = 0; b < run.result.block_count; ++b) {
        double delta, square_delta;
        CHECK(number->block_sample_count[b] == total / run.result.block_count,
              "number blocks count each sample once");
        block_sum += number->block_sum[b];
        block_square_sum += number->block_square_sum[b];
        delta = number->block_sum[b] / (double)number->block_sample_count[b] -
                number->mean;
        square_delta = number->block_square_sum[b] /
                (double)number->block_sample_count[b] - number->mean_square;
        se_sum += delta * delta;
        square_se_sum += square_delta * square_delta;
      }
      CHECK(block_sum == sum && block_square_sum == square_sum,
            "number block totals equal histogram moments");
      CHECK(fabs(number->standard_error - sqrt(se_sum /
                (run.result.block_count * (run.result.block_count - 1)))) < 1e-12 &&
            fabs(number->mean_square_standard_error - sqrt(square_se_sum /
                (run.result.block_count * (run.result.block_count - 1)))) < 1e-12,
            "number SE agrees with the independent block calculation");
    }
    CHECK(cabs(run.diagnostics.coefficient_energy -
               run.result.coefficient_energy) == 0.0,
          "diagnostic E_coeff equals the result field");
  }
  CHECK(all_ranks_equal_double(run.result.energy) &&
            all_ranks_equal_double(run.result.variance) &&
            all_ranks_equal_double(run.result.gc_number.mean) &&
            all_ranks_equal_double(run.result.gc_number.mean_square_standard_error) &&
            all_ranks_equal_double(creal(run.result.coefficient_energy)),
        "results are identical on every rank");
  run_gc(&chain, 0, 2, warm_up, samples, interval,
         MVMC_POWER_LANCZOS_GC_FAULT_NONE, &again);
  CHECK(again.status == MVMC_KRYLOV_STATUS_OK &&
            memcmp(&again.result.energy, &run.result.energy,
                   sizeof(double)) == 0 &&
            memcmp(&again.result.variance, &run.result.variance,
                   sizeof(double)) == 0 &&
            memcmp(again.result.gc_number.block_sum, run.result.gc_number.block_sum,
                   sizeof(run.result.gc_number.block_sum)) == 0 &&
            memcmp(again.result.gc_number.block_square_sum, run.result.gc_number.block_square_sum,
                   sizeof(run.result.gc_number.block_square_sum)) == 0,
        "the same seed reproduces the run bit for bit");
  free_chain(&chain);
}

static void test_faults_and_rejections(void) {
  Chain chain = make_chain(0);
  Run run;
  run_gc(&chain, 0, 2, 4, 256, 1,
         MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_NUMERIC, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_OK && run.result.valid,
        "a numeric-zero start candidate is redrawn");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.total_bootstrap_numeric_zero == 1 &&
              run.diagnostics.total_bootstrap_attempts >= 2 &&
              !run.diagnostics.any_bootstrap_exhausted,
          "bootstrap numeric-zero breakdown");
  }
  CHECK(run.counters.bootstrap_retries == 1, "one bootstrap retry");

  run_gc(&chain, 0, 2, 4, 256, 1,
         MVMC_POWER_LANCZOS_GC_FAULT_BOOTSTRAP_EXHAUSTED, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
            !run.result.valid && isnan(run.result.energy),
        "1024 unusable starts end the run without a result");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.any_bootstrap_exhausted &&
              run.diagnostics.total_bootstrap_attempts == 1024 &&
              run.diagnostics.failure_recorded &&
              run.diagnostics.failure_stage ==
                  MVMC_POWER_LANCZOS_GC_STAGE_FINAL_BOOTSTRAP,
          "exhausted bootstrap diagnostics");
  }

  run_gc(&chain, 0, 2, 4, 256, 1, MVMC_POWER_LANCZOS_GC_FAULT_FINAL_NUMERIC,
         &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
            !run.result.valid && !run.result.gc_number.valid,
        "a numeric-zero final proposal invalidates both energy and number");
  if (world_rank() == 0) {
    CHECK(run.diagnostics.failure_recorded &&
              run.diagnostics.failure_stage ==
                  MVMC_POWER_LANCZOS_GC_STAGE_FINAL_SAMPLE &&
              run.diagnostics.failure_amplitude_state ==
                  (int)MVMC_SCALED_COMPLEX_NUMERIC_ZERO &&
              isfinite(run.diagnostics.failure_log_error_bound) &&
              run.diagnostics.failure_electron_count >= 0,
          "final numeric-zero diagnostics");
  }

  {
    /* Rank-1 amplitude faults end every rank's run (no hang, no result). */
    const int target = world_size() > 1 ? 1 : 0;
    mvmc_gc_krylov_amplitude_testing_set_fault(
        MVMC_GC_KRYLOV_AMPLITUDE_FAULT_ALLOCATION, target);
    run_gc(&chain, 0, 2, 4, 256, 1, MVMC_POWER_LANCZOS_GC_FAULT_NONE, &run);
    CHECK(run.status == MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE &&
              !run.result.valid,
          "rank-1 amplitude allocation fault");
    mvmc_gc_krylov_amplitude_testing_set_fault(
        MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONFINITE, target);
    run_gc(&chain, 0, 2, 4, 256, 1, MVMC_POWER_LANCZOS_GC_FAULT_NONE, &run);
    CHECK(run.status == MVMC_KRYLOV_STATUS_NONFINITE && !run.result.valid,
          "rank-1 nonfinite amplitude fault");
    mvmc_gc_krylov_amplitude_testing_set_fault(
        MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONE, -1);
  }

#if defined(MVMC_HAVE_MPFR_REFINEMENT)
  if (world_size() > 1) {
    const MVMCPowerLanczosGCFault faults[] = {
        MVMC_POWER_LANCZOS_GC_FAULT_REFINEMENT_DECISION,
        MVMC_POWER_LANCZOS_GC_FAULT_REFINEMENT_STATUS};
    size_t i;
    for (i = 0; i < sizeof(faults) / sizeof(faults[0]); ++i) {
      run_gc(&chain, 0, 2, 4, 256, 1, faults[i], &run);
      CHECK(run.status == MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE &&
                !run.result.valid && !run.result.gc_number.valid,
            "rank-1 refinement entry disagreement fails without a hang");
    }
  }
#endif

  /* Decision A: zero H is rejected as in the canonical path. */
  run_gc(&chain, 1, 2, 4, 16, 1, MVMC_POWER_LANCZOS_GC_FAULT_NONE, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !run.result.valid,
        "zero H is rejected");

  run_gc(&chain, 0, 3, 4, 16, 1, MVMC_POWER_LANCZOS_GC_FAULT_NONE, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "odd initial N rejected");
  run_gc(&chain, 0, 6, 4, 16, 1, MVMC_POWER_LANCZOS_GC_FAULT_NONE, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "initial N above 2L rejected");
  run_gc(&chain, 0, 2, 4, 7, 1, MVMC_POWER_LANCZOS_GC_FAULT_NONE, &run);
  CHECK(run.status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "fewer than 8 samples rejected");
  {
    MVMCPowerLanczosGCInput input;
    MVMCPowerLanczosStabilizedResult result;
    memset(&input, 0, sizeof(input));
    input.world_communicator = world();
    input.chain_communicator = chain.chain;
    input.seed = 1;
    input.sample_count = 16;
    input.interval = 1;
    CHECK(mvmc_power_lanczos_gc_run(&input, &result) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              !result.valid,
          "missing GC view rejected");
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
  test_stage_counts(0);
  if (world_size() > 1) test_stage_counts(1);
  if (world_size() > 2) test_stage_counts(2);
  test_faults_and_rejections();
#ifdef _mpi_use
  MPI_Allreduce(&failure_count, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  total = failure_count;
#endif
  if (world_rank() == 0) {
    printf(total == 0 ? "gc_lanczos_driver_unit: PASS\n"
                      : "gc_lanczos_driver_unit: FAIL\n");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return total == 0 ? 0 : 1;
}
