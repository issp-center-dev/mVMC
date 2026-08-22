#ifndef MVMC_POWER_LANCZOS_OBSERVABLE_BLOCKS_H
#define MVMC_POWER_LANCZOS_OBSERVABLE_BLOCKS_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "krylov_fock_reference.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

typedef enum {
  MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT = 0,
  MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_FINAL = 1
} MVMCPowerLanczosObservableBlockStage;

typedef struct MVMCPowerLanczosObservableBlockAccumulator
    MVMCPowerLanczosObservableBlockAccumulator;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  MVMCPowerLanczosObservableBlockStage stage;
  size_t request_count;
  size_t entries_per_request;
  size_t block_capacity;
  uint64_t block_length;
  size_t completed_block_count;
  uint64_t completed_sample_count;
  uint64_t current_block_sample_count;
  uint64_t discarded_partial_sample_count;
  size_t payload_bytes;
  size_t allocated_bytes;
} MVMCPowerLanczosObservableBlockSummary;

MVMCKrylovStatus mvmc_power_lanczos_observable_block_create(
    MVMCPowerLanczosObservableBlockStage stage, size_t request_count,
    size_t block_capacity, uint64_t block_length,
    MVMCPowerLanczosObservableBlockAccumulator **accumulator);

void mvmc_power_lanczos_observable_block_destroy(
    MVMCPowerLanczosObservableBlockAccumulator *accumulator);

MVMCKrylovStatus mvmc_power_lanczos_observable_block_add_sample(
    MVMCPowerLanczosObservableBlockAccumulator *accumulator,
    const double complex *entries, size_t entry_count);

/* Discard the current incomplete block and retain an auditable count. */
MVMCKrylovStatus mvmc_power_lanczos_observable_block_discard_partial(
    MVMCPowerLanczosObservableBlockAccumulator *accumulator);

MVMCKrylovStatus mvmc_power_lanczos_observable_block_summary(
    const MVMCPowerLanczosObservableBlockAccumulator *accumulator,
    MVMCPowerLanczosObservableBlockSummary *summary);

/* Export completed blocks in block-major, request-major entry order. */
MVMCKrylovStatus mvmc_power_lanczos_observable_block_export(
    const MVMCPowerLanczosObservableBlockAccumulator *accumulator,
    double complex *block_entries, size_t block_entry_capacity,
    uint64_t *block_sample_counts, size_t block_count_capacity);

/*
 * Deterministically merge completed, aligned rank-local blocks in the given
 * rank order.  inputs must have no partial block; output must be empty and
 * shape-compatible.  Failure does not modify output.
 */
MVMCKrylovStatus mvmc_power_lanczos_observable_block_reduce_rank_ordered(
    const MVMCPowerLanczosObservableBlockAccumulator *const *inputs,
    size_t input_count,
    MVMCPowerLanczosObservableBlockAccumulator *output);

MVMCKrylovStatus mvmc_power_lanczos_observable_block_payload_bytes(
    MVMCPowerLanczosObservableBlockStage stage, size_t request_count,
    size_t block_count, size_t *payload_bytes);

#endif

#endif
