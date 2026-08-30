#ifndef MVMC_POWER_LANCZOS_SNAPSHOT_H
#define MVMC_POWER_LANCZOS_SNAPSHOT_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "krylov_fock_reference.h"
#include "power_lanczos_sha256.h"

#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_SNAPSHOT_VERSION UINT64_C(1)

typedef enum {
  MVMC_POWER_LANCZOS_SNAPSHOT_READY = 0,
  MVMC_POWER_LANCZOS_SNAPSHOT_COEFFICIENT_ACTIVE,
  MVMC_POWER_LANCZOS_SNAPSHOT_COEFFICIENT_COMPLETE,
  MVMC_POWER_LANCZOS_SNAPSHOT_FINAL_ACTIVE,
  MVMC_POWER_LANCZOS_SNAPSHOT_FINAL_COMPLETE,
  MVMC_POWER_LANCZOS_SNAPSHOT_FAILED
} MVMCPowerLanczosSnapshotState;

#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
typedef enum {
  MVMC_POWER_LANCZOS_SNAPSHOT_TEST_BASE = 0,
  MVMC_POWER_LANCZOS_SNAPSHOT_TEST_COEFFICIENT,
  MVMC_POWER_LANCZOS_SNAPSHOT_TEST_FINAL
} MVMCPowerLanczosSnapshotTestingTarget;
#endif

typedef struct MVMCPowerLanczosSnapshot MVMCPowerLanczosSnapshot;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  MVMCPowerLanczosSnapshotState state;
  uint64_t version;
  uint64_t generation;
  size_t word_count;
  size_t configuration_bit_count;
  size_t allocated_bytes;
  unsigned char base_sha256[MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES];
  unsigned char coefficient_sha256[MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES];
  unsigned char final_sha256[MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES];
} MVMCPowerLanczosSnapshotSummary;

/*
 * The snapshot SHA-256 covers the NUL-terminated domain
 * "MVMC-P6-SNAPSHOT-V1", followed by version, word_count, bit_count,
 * generation and every canonical word as unsigned 64-bit big-endian values.
 * Unused bits in the last word must be zero.
 */
MVMCKrylovStatus mvmc_power_lanczos_snapshot_digest(
    const uint64_t *words, size_t word_count, size_t bit_count,
    uint64_t generation,
    unsigned char digest[MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES]);

MVMCKrylovStatus mvmc_power_lanczos_snapshot_create(
    const uint64_t *canonical_words, size_t word_count,
    size_t configuration_bit_count, uint64_t generation,
    MVMCPowerLanczosSnapshot **snapshot);

void mvmc_power_lanczos_snapshot_destroy(
    MVMCPowerLanczosSnapshot *snapshot);

MVMCKrylovStatus mvmc_power_lanczos_snapshot_verify(
    MVMCPowerLanczosSnapshot *snapshot);

MVMCKrylovStatus mvmc_power_lanczos_snapshot_coefficient_begin(
    MVMCPowerLanczosSnapshot *snapshot, uint64_t **mutable_words,
    size_t *word_count);

MVMCKrylovStatus mvmc_power_lanczos_snapshot_coefficient_complete(
    MVMCPowerLanczosSnapshot *snapshot);

MVMCKrylovStatus mvmc_power_lanczos_snapshot_final_begin(
    MVMCPowerLanczosSnapshot *snapshot, uint64_t **mutable_words,
    size_t *word_count);

MVMCKrylovStatus mvmc_power_lanczos_snapshot_final_complete(
    MVMCPowerLanczosSnapshot *snapshot);

MVMCKrylovStatus mvmc_power_lanczos_snapshot_summary(
    const MVMCPowerLanczosSnapshot *snapshot,
    MVMCPowerLanczosSnapshotSummary *summary);

#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
MVMCKrylovStatus mvmc_power_lanczos_snapshot_testing_corrupt(
    MVMCPowerLanczosSnapshot *snapshot,
    MVMCPowerLanczosSnapshotTestingTarget target,
    size_t word_index, uint64_t xor_mask);
#endif

#endif

#endif
