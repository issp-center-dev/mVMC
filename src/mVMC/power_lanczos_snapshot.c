#include "power_lanczos_snapshot.h"

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct MVMCPowerLanczosSnapshot {
  int valid;
  MVMCKrylovStatus status;
  MVMCPowerLanczosSnapshotState state;
  uint64_t generation;
  size_t word_count;
  size_t configuration_bit_count;
  size_t allocated_bytes;
  uint64_t *base_words;
  uint64_t *coefficient_words;
  uint64_t *final_words;
  unsigned char base_sha256[32];
  unsigned char coefficient_sha256[32];
  unsigned char final_sha256[32];
};

static const unsigned char SnapshotDomain[] = "MVMC-P6-SNAPSHOT-V1";

static void store_u64_be(uint64_t value, unsigned char output[8]) {
  int index;
  for (index = 7; index >= 0; --index) {
    output[index] = (unsigned char)value;
    value >>= 8;
  }
}

static int range_overlaps(const void *left, size_t left_size,
                          const void *right, size_t right_size) {
  const uintptr_t left_begin = (uintptr_t)left;
  const uintptr_t right_begin = (uintptr_t)right;
  uintptr_t left_end;
  uintptr_t right_end;
  if (left == NULL || right == NULL ||
      left_begin > UINTPTR_MAX - left_size ||
      right_begin > UINTPTR_MAX - right_size) {
    return 1;
  }
  left_end = left_begin + left_size;
  right_end = right_begin + right_size;
  return left_begin < right_end && right_begin < left_end;
}

static int valid_shape(const uint64_t *words, size_t word_count,
                       size_t bit_count, uint64_t generation) {
  size_t required_words;
  size_t remainder;
  if (words == NULL || word_count == 0 || bit_count == 0 ||
      generation == 0 || bit_count > (size_t)UINT64_MAX ||
      word_count > (size_t)UINT64_MAX) {
    return 0;
  }
  required_words = (bit_count - 1) / 64 + 1;
  if (required_words != word_count) return 0;
  remainder = bit_count % 64;
  if (remainder != 0) {
    const uint64_t mask = (UINT64_C(1) << remainder) - UINT64_C(1);
    if ((words[word_count - 1] & ~mask) != 0) return 0;
  }
  return 1;
}

MVMCKrylovStatus mvmc_power_lanczos_snapshot_digest(
    const uint64_t *words, size_t word_count, size_t bit_count,
    uint64_t generation, unsigned char digest[32]) {
  size_t word_bytes;
  size_t header_bytes = sizeof(SnapshotDomain) + 4 * sizeof(uint64_t);
  size_t total_bytes;
  unsigned char *record;
  unsigned char *cursor;
  size_t word;
  if (!valid_shape(words, word_count, bit_count, generation) ||
      digest == NULL || word_count > SIZE_MAX / sizeof(uint64_t)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  word_bytes = word_count * sizeof(uint64_t);
  if (word_bytes > SIZE_MAX - header_bytes) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  total_bytes = header_bytes + word_bytes;
  record = (unsigned char *)malloc(total_bytes);
  if (record == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  cursor = record;
  memcpy(cursor, SnapshotDomain, sizeof(SnapshotDomain));
  cursor += sizeof(SnapshotDomain);
  store_u64_be(MVMC_POWER_LANCZOS_SNAPSHOT_VERSION, cursor);
  cursor += sizeof(uint64_t);
  store_u64_be((uint64_t)word_count, cursor);
  cursor += sizeof(uint64_t);
  store_u64_be((uint64_t)bit_count, cursor);
  cursor += sizeof(uint64_t);
  store_u64_be(generation, cursor);
  cursor += sizeof(uint64_t);
  for (word = 0; word < word_count; ++word) {
    store_u64_be(words[word], cursor);
    cursor += sizeof(uint64_t);
  }
  if (!mvmc_power_lanczos_sha256_bytes(record, total_bytes, digest)) {
    memset(record, 0, total_bytes);
    free(record);
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  memset(record, 0, total_bytes);
  free(record);
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus verify_digest(
    const MVMCPowerLanczosSnapshot *snapshot, const uint64_t *words,
    const unsigned char expected[32]) {
  unsigned char actual[32];
  MVMCKrylovStatus status = mvmc_power_lanczos_snapshot_digest(
      words, snapshot->word_count, snapshot->configuration_bit_count,
      snapshot->generation, actual);
  if (status == MVMC_KRYLOV_STATUS_OK &&
      memcmp(actual, expected, sizeof(actual)) != 0) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  memset(actual, 0, sizeof(actual));
  return status;
}

static MVMCKrylovStatus fail_snapshot(
    MVMCPowerLanczosSnapshot *snapshot, MVMCKrylovStatus status) {
  if (snapshot != NULL) {
    snapshot->valid = 0;
    snapshot->status = status;
    snapshot->state = MVMC_POWER_LANCZOS_SNAPSHOT_FAILED;
  }
  return status;
}

MVMCKrylovStatus mvmc_power_lanczos_snapshot_create(
    const uint64_t *canonical_words, size_t word_count,
    size_t configuration_bit_count, uint64_t generation,
    MVMCPowerLanczosSnapshot **snapshot) {
  MVMCPowerLanczosSnapshot *candidate;
  MVMCKrylovStatus status;
  size_t bytes;
  if (snapshot == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *snapshot = NULL;
  if (!valid_shape(canonical_words, word_count, configuration_bit_count,
                   generation) ||
      word_count > SIZE_MAX / sizeof(uint64_t)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  bytes = word_count * sizeof(uint64_t);
  if (bytes > (SIZE_MAX - sizeof(*candidate)) / 3) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  candidate = (MVMCPowerLanczosSnapshot *)calloc(1, sizeof(*candidate));
  if (candidate == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  candidate->base_words = (uint64_t *)malloc(bytes);
  candidate->coefficient_words = (uint64_t *)malloc(bytes);
  candidate->final_words = (uint64_t *)malloc(bytes);
  if (candidate->base_words == NULL || candidate->coefficient_words == NULL ||
      candidate->final_words == NULL) {
    mvmc_power_lanczos_snapshot_destroy(candidate);
    return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
  if (range_overlaps(candidate->base_words, bytes,
                     candidate->coefficient_words, bytes) ||
      range_overlaps(candidate->base_words, bytes,
                     candidate->final_words, bytes) ||
      range_overlaps(candidate->coefficient_words, bytes,
                     candidate->final_words, bytes)) {
    mvmc_power_lanczos_snapshot_destroy(candidate);
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  memcpy(candidate->base_words, canonical_words, bytes);
  memcpy(candidate->coefficient_words, canonical_words, bytes);
  memcpy(candidate->final_words, canonical_words, bytes);
  candidate->generation = generation;
  candidate->word_count = word_count;
  candidate->configuration_bit_count = configuration_bit_count;
  candidate->allocated_bytes = sizeof(*candidate) + 3 * bytes;
  status = mvmc_power_lanczos_snapshot_digest(
      candidate->base_words, word_count, configuration_bit_count,
      generation, candidate->base_sha256);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    mvmc_power_lanczos_snapshot_destroy(candidate);
    return status;
  }
  memcpy(candidate->coefficient_sha256, candidate->base_sha256,
         sizeof(candidate->base_sha256));
  memcpy(candidate->final_sha256, candidate->base_sha256,
         sizeof(candidate->base_sha256));
  candidate->valid = 1;
  candidate->status = MVMC_KRYLOV_STATUS_OK;
  candidate->state = MVMC_POWER_LANCZOS_SNAPSHOT_READY;
  *snapshot = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_power_lanczos_snapshot_destroy(
    MVMCPowerLanczosSnapshot *snapshot) {
  size_t bytes = 0;
  if (snapshot == NULL) return;
  if (snapshot->word_count <= SIZE_MAX / sizeof(uint64_t)) {
    bytes = snapshot->word_count * sizeof(uint64_t);
  }
  if (snapshot->base_words != NULL) {
    memset(snapshot->base_words, 0, bytes);
    free(snapshot->base_words);
  }
  if (snapshot->coefficient_words != NULL) {
    memset(snapshot->coefficient_words, 0, bytes);
    free(snapshot->coefficient_words);
  }
  if (snapshot->final_words != NULL) {
    memset(snapshot->final_words, 0, bytes);
    free(snapshot->final_words);
  }
  memset(snapshot, 0, sizeof(*snapshot));
  free(snapshot);
}

MVMCKrylovStatus mvmc_power_lanczos_snapshot_verify(
    MVMCPowerLanczosSnapshot *snapshot) {
  MVMCKrylovStatus status;
  if (snapshot == NULL || !snapshot->valid ||
      snapshot->status != MVMC_KRYLOV_STATUS_OK ||
      snapshot->state == MVMC_POWER_LANCZOS_SNAPSHOT_FAILED) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = verify_digest(snapshot, snapshot->base_words,
                         snapshot->base_sha256);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return fail_snapshot(snapshot, status);
  }
  if (snapshot->state != MVMC_POWER_LANCZOS_SNAPSHOT_FINAL_ACTIVE) {
    status = verify_digest(snapshot, snapshot->final_words,
                           snapshot->final_sha256);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      return fail_snapshot(snapshot, status);
    }
  }
  if (snapshot->state == MVMC_POWER_LANCZOS_SNAPSHOT_READY ||
      snapshot->state == MVMC_POWER_LANCZOS_SNAPSHOT_COEFFICIENT_COMPLETE ||
      snapshot->state == MVMC_POWER_LANCZOS_SNAPSHOT_FINAL_ACTIVE ||
      snapshot->state == MVMC_POWER_LANCZOS_SNAPSHOT_FINAL_COMPLETE) {
    status = verify_digest(snapshot, snapshot->coefficient_words,
                           snapshot->coefficient_sha256);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      return fail_snapshot(snapshot, status);
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_snapshot_coefficient_begin(
    MVMCPowerLanczosSnapshot *snapshot, uint64_t **mutable_words,
    size_t *word_count) {
  MVMCKrylovStatus status;
  if (mutable_words == NULL || word_count == NULL || snapshot == NULL ||
      snapshot->state != MVMC_POWER_LANCZOS_SNAPSHOT_READY) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *mutable_words = NULL;
  *word_count = 0;
  status = mvmc_power_lanczos_snapshot_verify(snapshot);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  snapshot->state = MVMC_POWER_LANCZOS_SNAPSHOT_COEFFICIENT_ACTIVE;
  *mutable_words = snapshot->coefficient_words;
  *word_count = snapshot->word_count;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_snapshot_coefficient_complete(
    MVMCPowerLanczosSnapshot *snapshot) {
  MVMCKrylovStatus status;
  if (snapshot == NULL || !snapshot->valid ||
      snapshot->state != MVMC_POWER_LANCZOS_SNAPSHOT_COEFFICIENT_ACTIVE) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = verify_digest(snapshot, snapshot->base_words,
                         snapshot->base_sha256);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = verify_digest(snapshot, snapshot->final_words,
                           snapshot->final_sha256);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_snapshot_digest(
        snapshot->coefficient_words, snapshot->word_count,
        snapshot->configuration_bit_count, snapshot->generation,
        snapshot->coefficient_sha256);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return fail_snapshot(snapshot, status);
  }
  snapshot->state = MVMC_POWER_LANCZOS_SNAPSHOT_COEFFICIENT_COMPLETE;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_snapshot_final_begin(
    MVMCPowerLanczosSnapshot *snapshot, uint64_t **mutable_words,
    size_t *word_count) {
  MVMCKrylovStatus status;
  if (mutable_words == NULL || word_count == NULL || snapshot == NULL ||
      snapshot->state != MVMC_POWER_LANCZOS_SNAPSHOT_COEFFICIENT_COMPLETE) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *mutable_words = NULL;
  *word_count = 0;
  status = mvmc_power_lanczos_snapshot_verify(snapshot);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  snapshot->state = MVMC_POWER_LANCZOS_SNAPSHOT_FINAL_ACTIVE;
  *mutable_words = snapshot->final_words;
  *word_count = snapshot->word_count;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_snapshot_final_complete(
    MVMCPowerLanczosSnapshot *snapshot) {
  MVMCKrylovStatus status;
  if (snapshot == NULL || !snapshot->valid ||
      snapshot->state != MVMC_POWER_LANCZOS_SNAPSHOT_FINAL_ACTIVE) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = verify_digest(snapshot, snapshot->base_words,
                         snapshot->base_sha256);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = verify_digest(snapshot, snapshot->coefficient_words,
                           snapshot->coefficient_sha256);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_snapshot_digest(
        snapshot->final_words, snapshot->word_count,
        snapshot->configuration_bit_count, snapshot->generation,
        snapshot->final_sha256);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return fail_snapshot(snapshot, status);
  }
  snapshot->state = MVMC_POWER_LANCZOS_SNAPSHOT_FINAL_COMPLETE;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_snapshot_summary(
    const MVMCPowerLanczosSnapshot *snapshot,
    MVMCPowerLanczosSnapshotSummary *summary) {
  MVMCPowerLanczosSnapshotSummary candidate;
  if (summary == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(summary, 0, sizeof(*summary));
  summary->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (snapshot == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = snapshot->valid;
  candidate.status = snapshot->status;
  candidate.state = snapshot->state;
  candidate.version = MVMC_POWER_LANCZOS_SNAPSHOT_VERSION;
  candidate.generation = snapshot->generation;
  candidate.word_count = snapshot->word_count;
  candidate.configuration_bit_count = snapshot->configuration_bit_count;
  candidate.allocated_bytes = snapshot->allocated_bytes;
  memcpy(candidate.base_sha256, snapshot->base_sha256,
         sizeof(candidate.base_sha256));
  memcpy(candidate.coefficient_sha256, snapshot->coefficient_sha256,
         sizeof(candidate.coefficient_sha256));
  memcpy(candidate.final_sha256, snapshot->final_sha256,
         sizeof(candidate.final_sha256));
  *summary = candidate;
  return candidate.status;
}

#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
MVMCKrylovStatus mvmc_power_lanczos_snapshot_testing_corrupt(
    MVMCPowerLanczosSnapshot *snapshot,
    MVMCPowerLanczosSnapshotTestingTarget target,
    size_t word_index, uint64_t xor_mask) {
  uint64_t *words = NULL;
  if (snapshot == NULL || !snapshot->valid ||
      word_index >= snapshot->word_count || xor_mask == 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  switch (target) {
    case MVMC_POWER_LANCZOS_SNAPSHOT_TEST_BASE:
      words = snapshot->base_words;
      break;
    case MVMC_POWER_LANCZOS_SNAPSHOT_TEST_COEFFICIENT:
      words = snapshot->coefficient_words;
      break;
    case MVMC_POWER_LANCZOS_SNAPSHOT_TEST_FINAL:
      words = snapshot->final_words;
      break;
    default:
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  words[word_index] ^= xor_mask;
  return MVMC_KRYLOV_STATUS_OK;
}
#endif

#endif
