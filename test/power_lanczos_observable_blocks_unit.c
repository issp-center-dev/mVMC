#include "power_lanczos_observable_blocks.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

static void Check(int condition, const char *message) {
  if (!condition) {
    fprintf(stderr, "FAIL: %s\n", message);
    ++failures;
  }
}

static int Close(double complex left, double complex right) {
  return cabs(left - right) <= 1.0e-14 * (1.0 + cabs(right));
}

int main(void) {
  MVMCPowerLanczosObservableBlockAccumulator *rank0 = NULL;
  MVMCPowerLanczosObservableBlockAccumulator *rank1 = NULL;
  MVMCPowerLanczosObservableBlockAccumulator *reduced = NULL;
  MVMCPowerLanczosObservableBlockAccumulator *partial = NULL;
  MVMCPowerLanczosObservableBlockSummary summary;
  const MVMCPowerLanczosObservableBlockAccumulator *inputs[2];
  double complex sample[4];
  double complex exported[4];
  uint64_t counts[1];
  size_t bytes = 0;
  MVMCKrylovStatus status;
  status = mvmc_power_lanczos_observable_block_payload_bytes(
      MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT, 7, 11, &bytes);
  Check(status == MVMC_KRYLOV_STATUS_OK && bytes == 64 * 7 * 11,
        "coefficient payload formula");
  status = mvmc_power_lanczos_observable_block_payload_bytes(
      MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_FINAL, 7, 11, &bytes);
  Check(status == MVMC_KRYLOV_STATUS_OK && bytes == 16 * 7 * 11,
        "final payload formula");

  status = mvmc_power_lanczos_observable_block_create(
      MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT, 1, 1, 3,
      &rank0);
  Check(status == MVMC_KRYLOV_STATUS_OK, "rank0 create");
  status = mvmc_power_lanczos_observable_block_create(
      MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT, 1, 1, 3,
      &rank1);
  Check(status == MVMC_KRYLOV_STATUS_OK, "rank1 create");
  status = mvmc_power_lanczos_observable_block_create(
      MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT, 1, 1, 3,
      &reduced);
  Check(status == MVMC_KRYLOV_STATUS_OK, "reduced create");

  sample[0] = 1.0e16 + 1.0e16 * I;
  sample[1] = 2.0;
  sample[2] = 3.0;
  sample[3] = 4.0;
  Check(mvmc_power_lanczos_observable_block_add_sample(rank0, sample, 4) ==
            MVMC_KRYLOV_STATUS_OK,
        "rank0 sample 1");
  sample[0] = 1.0 + I;
  Check(mvmc_power_lanczos_observable_block_add_sample(rank0, sample, 4) ==
            MVMC_KRYLOV_STATUS_OK,
        "rank0 sample 2");
  sample[0] = -1.0e16 - 1.0e16 * I;
  Check(mvmc_power_lanczos_observable_block_add_sample(rank0, sample, 4) ==
            MVMC_KRYLOV_STATUS_OK,
        "rank0 sample 3");

  sample[0] = 2.0 + 3.0 * I;
  sample[1] = -1.0;
  sample[2] = 0.5;
  sample[3] = -4.0;
  Check(mvmc_power_lanczos_observable_block_add_sample(rank1, sample, 4) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_observable_block_add_sample(rank1, sample, 4) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_observable_block_add_sample(rank1, sample, 4) ==
                MVMC_KRYLOV_STATUS_OK,
        "rank1 block samples");

  inputs[0] = rank0;
  inputs[1] = rank1;
  status = mvmc_power_lanczos_observable_block_reduce_rank_ordered(
      inputs, 2, reduced);
  Check(status == MVMC_KRYLOV_STATUS_OK, "rank-ordered reduction");
  status = mvmc_power_lanczos_observable_block_export(
      reduced, exported, 4, counts, 1);
  Check(status == MVMC_KRYLOV_STATUS_OK && counts[0] == 6,
        "reduced block count");
  Check(Close(exported[0], 7.0 + 10.0 * I),
        "compensated and rank-ordered entry 0");
  Check(Close(exported[1], 3.0) && Close(exported[2], 10.5) &&
            Close(exported[3], 0.0),
        "reduced remaining entries");
  status = mvmc_power_lanczos_observable_block_summary(reduced, &summary);
  Check(status == MVMC_KRYLOV_STATUS_OK && summary.valid &&
            summary.completed_block_count == 1 &&
            summary.completed_sample_count == 6 &&
            summary.payload_bytes == 64 && summary.allocated_bytes > 64,
        "reduced summary");

  status = mvmc_power_lanczos_observable_block_create(
      MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_FINAL, 1, 2, 3, &partial);
  Check(status == MVMC_KRYLOV_STATUS_OK, "partial create");
  sample[0] = 5.0;
  Check(mvmc_power_lanczos_observable_block_add_sample(partial, sample, 1) ==
            MVMC_KRYLOV_STATUS_OK,
        "partial sample");
  sample[0] = NAN;
  Check(mvmc_power_lanczos_observable_block_add_sample(partial, sample, 1) ==
            MVMC_KRYLOV_STATUS_NONFINITE,
        "nonfinite sample accepted");
  Check(mvmc_power_lanczos_observable_block_discard_partial(partial) ==
            MVMC_KRYLOV_STATUS_OK,
        "partial discard");
  Check(mvmc_power_lanczos_observable_block_summary(partial, &summary) ==
                MVMC_KRYLOV_STATUS_OK &&
            summary.current_block_sample_count == 0 &&
            summary.discarded_partial_sample_count == 1 &&
            summary.completed_block_count == 0,
        "partial discard summary");

  mvmc_power_lanczos_observable_block_destroy(partial);
  mvmc_power_lanczos_observable_block_destroy(reduced);
  mvmc_power_lanczos_observable_block_destroy(rank1);
  mvmc_power_lanczos_observable_block_destroy(rank0);
  if (failures != 0) {
    fprintf(stderr, "%d observable block assertion(s) failed\n", failures);
    return EXIT_FAILURE;
  }
  puts("power_lanczos_observable_blocks_unit: PASS");
  return EXIT_SUCCESS;
}
