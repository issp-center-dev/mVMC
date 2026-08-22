#include "power_lanczos_observable_blocks.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "power_lanczos_observable_blocks.c requires the bounded engine"
#endif

#include "power_lanczos_observable_census.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  double complex sum;
  double complex compensation;
} CompensatedComplex;

struct MVMCPowerLanczosObservableBlockAccumulator {
  MVMCPowerLanczosObservableBlockStage stage;
  size_t request_count;
  size_t entries_per_request;
  size_t entry_count;
  size_t block_capacity;
  uint64_t block_length;
  size_t completed_block_count;
  uint64_t completed_sample_count;
  uint64_t current_block_sample_count;
  uint64_t discarded_partial_sample_count;
  size_t payload_bytes;
  size_t allocated_bytes;
  CompensatedComplex *completed;
  uint64_t *completed_counts;
  CompensatedComplex *current;
  CompensatedComplex *transaction;
};

_Static_assert(sizeof(double complex) == 16,
               "observable block payload contract requires complex128");

static int CheckedAddSize(size_t *total, size_t addition) {
  if (addition > SIZE_MAX - *total) return 0;
  *total += addition;
  return 1;
}

static int CheckedMultiplySize(size_t left, size_t right, size_t *product) {
  if (left != 0 && right > SIZE_MAX / left) return 0;
  *product = left * right;
  return 1;
}

static int CheckedAddU64(uint64_t *total, uint64_t addition) {
  if (addition > UINT64_MAX - *total) return 0;
  *total += addition;
  return 1;
}

static int FiniteComplex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static int AddReal(double value, double *sum, double *compensation) {
  double updated;
  if (!isfinite(value) || !isfinite(*sum) || !isfinite(*compensation)) {
    return 0;
  }
  updated = *sum + value;
  if (!isfinite(updated)) return 0;
  if (fabs(*sum) >= fabs(value)) {
    *compensation += (*sum - updated) + value;
  } else {
    *compensation += (value - updated) + *sum;
  }
  if (!isfinite(*compensation)) return 0;
  *sum = updated;
  return 1;
}

static int AddComplex(double complex value, CompensatedComplex *target) {
  double real_sum = creal(target->sum);
  double imaginary_sum = cimag(target->sum);
  double real_compensation = creal(target->compensation);
  double imaginary_compensation = cimag(target->compensation);
  if (!FiniteComplex(value) ||
      !AddReal(creal(value), &real_sum, &real_compensation) ||
      !AddReal(cimag(value), &imaginary_sum, &imaginary_compensation)) {
    return 0;
  }
  target->sum = real_sum + I * imaginary_sum;
  target->compensation =
      real_compensation + I * imaginary_compensation;
  return 1;
}

static double complex SumValue(const CompensatedComplex *value) {
  return value->sum + value->compensation;
}

static int ValidStage(MVMCPowerLanczosObservableBlockStage stage) {
  return stage == MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT ||
         stage == MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_FINAL;
}

static size_t EntriesPerRequest(
    MVMCPowerLanczosObservableBlockStage stage) {
  return stage == MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT ? 4 : 1;
}

MVMCKrylovStatus mvmc_power_lanczos_observable_block_payload_bytes(
    MVMCPowerLanczosObservableBlockStage stage, size_t request_count,
    size_t block_count, size_t *payload_bytes) {
  size_t entries;
  if (payload_bytes == NULL || !ValidStage(stage) || request_count == 0 ||
      request_count > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS ||
      block_count == 0 ||
      !CheckedMultiplySize(request_count, EntriesPerRequest(stage),
                           &entries) ||
      !CheckedMultiplySize(entries, block_count, &entries) ||
      !CheckedMultiplySize(entries, sizeof(double complex), payload_bytes)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_observable_block_create(
    MVMCPowerLanczosObservableBlockStage stage, size_t request_count,
    size_t block_capacity, uint64_t block_length,
    MVMCPowerLanczosObservableBlockAccumulator **accumulator) {
  MVMCPowerLanczosObservableBlockAccumulator *candidate;
  size_t completed_entry_count;
  size_t bytes;
  if (accumulator == NULL || *accumulator != NULL || !ValidStage(stage) ||
      request_count == 0 ||
      request_count > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS ||
      block_capacity == 0 || block_length == 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  candidate = (MVMCPowerLanczosObservableBlockAccumulator *)calloc(
      1, sizeof(*candidate));
  if (candidate == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  candidate->stage = stage;
  candidate->request_count = request_count;
  candidate->entries_per_request = EntriesPerRequest(stage);
  candidate->block_capacity = block_capacity;
  candidate->block_length = block_length;
  candidate->allocated_bytes = sizeof(*candidate);
  if (!CheckedMultiplySize(request_count, candidate->entries_per_request,
                           &candidate->entry_count) ||
      !CheckedMultiplySize(block_capacity, candidate->entry_count,
                           &completed_entry_count)) {
    free(candidate);
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  if (mvmc_power_lanczos_observable_block_payload_bytes(
          stage, request_count, block_capacity,
          &candidate->payload_bytes) != MVMC_KRYLOV_STATUS_OK) {
    free(candidate);
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  if (!CheckedMultiplySize(completed_entry_count,
                           sizeof(candidate->completed[0]), &bytes) ||
      !CheckedAddSize(&candidate->allocated_bytes, bytes)) {
    free(candidate);
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  candidate->completed = (CompensatedComplex *)calloc(
      completed_entry_count, sizeof(candidate->completed[0]));
  if (!CheckedMultiplySize(block_capacity,
                           sizeof(candidate->completed_counts[0]), &bytes) ||
      !CheckedAddSize(&candidate->allocated_bytes, bytes)) {
    mvmc_power_lanczos_observable_block_destroy(candidate);
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  candidate->completed_counts =
      (uint64_t *)calloc(block_capacity, sizeof(candidate->completed_counts[0]));
  if (!CheckedMultiplySize(candidate->entry_count,
                           sizeof(candidate->current[0]), &bytes) ||
      !CheckedAddSize(&candidate->allocated_bytes, bytes)) {
    mvmc_power_lanczos_observable_block_destroy(candidate);
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  candidate->current = (CompensatedComplex *)calloc(
      candidate->entry_count, sizeof(candidate->current[0]));
  if (!CheckedAddSize(&candidate->allocated_bytes, bytes)) {
    mvmc_power_lanczos_observable_block_destroy(candidate);
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  candidate->transaction = (CompensatedComplex *)calloc(
      candidate->entry_count, sizeof(candidate->transaction[0]));
  if (candidate->completed == NULL || candidate->completed_counts == NULL ||
      candidate->current == NULL || candidate->transaction == NULL) {
    mvmc_power_lanczos_observable_block_destroy(candidate);
    return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
  *accumulator = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_power_lanczos_observable_block_destroy(
    MVMCPowerLanczosObservableBlockAccumulator *accumulator) {
  if (accumulator == NULL) return;
  free(accumulator->transaction);
  free(accumulator->current);
  free(accumulator->completed_counts);
  free(accumulator->completed);
  memset(accumulator, 0, sizeof(*accumulator));
  free(accumulator);
}

MVMCKrylovStatus mvmc_power_lanczos_observable_block_add_sample(
    MVMCPowerLanczosObservableBlockAccumulator *accumulator,
    const double complex *entries, size_t entry_count) {
  size_t entry;
  uint64_t next_completed_sample_count;
  if (accumulator == NULL || entries == NULL ||
      entry_count != accumulator->entry_count ||
      accumulator->completed_block_count >= accumulator->block_capacity ||
      accumulator->current_block_sample_count >= accumulator->block_length) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memcpy(accumulator->transaction, accumulator->current,
         entry_count * sizeof(accumulator->transaction[0]));
  for (entry = 0; entry < entry_count; ++entry) {
    if (!AddComplex(entries[entry], &accumulator->transaction[entry])) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  next_completed_sample_count = accumulator->completed_sample_count;
  if (accumulator->current_block_sample_count + 1 ==
          accumulator->block_length &&
      !CheckedAddU64(&next_completed_sample_count,
                     accumulator->block_length)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  memcpy(accumulator->current, accumulator->transaction,
         entry_count * sizeof(accumulator->current[0]));
  ++accumulator->current_block_sample_count;
  if (accumulator->current_block_sample_count == accumulator->block_length) {
    const size_t block = accumulator->completed_block_count;
    memcpy(accumulator->completed + block * entry_count,
           accumulator->current,
           entry_count * sizeof(accumulator->completed[0]));
    accumulator->completed_counts[block] = accumulator->block_length;
    ++accumulator->completed_block_count;
    accumulator->completed_sample_count = next_completed_sample_count;
    accumulator->current_block_sample_count = 0;
    memset(accumulator->current, 0,
           entry_count * sizeof(accumulator->current[0]));
  }
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_observable_block_discard_partial(
    MVMCPowerLanczosObservableBlockAccumulator *accumulator) {
  uint64_t discarded;
  if (accumulator == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  discarded = accumulator->discarded_partial_sample_count;
  if (!CheckedAddU64(&discarded,
                     accumulator->current_block_sample_count)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  accumulator->discarded_partial_sample_count = discarded;
  accumulator->current_block_sample_count = 0;
  memset(accumulator->current, 0,
         accumulator->entry_count * sizeof(accumulator->current[0]));
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_observable_block_summary(
    const MVMCPowerLanczosObservableBlockAccumulator *accumulator,
    MVMCPowerLanczosObservableBlockSummary *summary) {
  MVMCPowerLanczosObservableBlockSummary candidate;
  if (accumulator == NULL || summary == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.stage = accumulator->stage;
  candidate.request_count = accumulator->request_count;
  candidate.entries_per_request = accumulator->entries_per_request;
  candidate.block_capacity = accumulator->block_capacity;
  candidate.block_length = accumulator->block_length;
  candidate.completed_block_count = accumulator->completed_block_count;
  candidate.completed_sample_count = accumulator->completed_sample_count;
  candidate.current_block_sample_count =
      accumulator->current_block_sample_count;
  candidate.discarded_partial_sample_count =
      accumulator->discarded_partial_sample_count;
  candidate.payload_bytes = accumulator->payload_bytes;
  candidate.allocated_bytes = accumulator->allocated_bytes;
  *summary = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_observable_block_export(
    const MVMCPowerLanczosObservableBlockAccumulator *accumulator,
    double complex *block_entries, size_t block_entry_capacity,
    uint64_t *block_sample_counts, size_t block_count_capacity) {
  size_t required_entries;
  size_t entry;
  if (accumulator == NULL || block_entries == NULL ||
      block_sample_counts == NULL ||
      !CheckedMultiplySize(accumulator->completed_block_count,
                           accumulator->entry_count, &required_entries) ||
      block_entry_capacity < required_entries ||
      block_count_capacity < accumulator->completed_block_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (entry = 0; entry < required_entries; ++entry) {
    if (!FiniteComplex(SumValue(&accumulator->completed[entry]))) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  for (entry = 0; entry < required_entries; ++entry) {
    block_entries[entry] = SumValue(&accumulator->completed[entry]);
  }
  memcpy(block_sample_counts, accumulator->completed_counts,
         accumulator->completed_block_count * sizeof(block_sample_counts[0]));
  return MVMC_KRYLOV_STATUS_OK;
}

static int SameShape(
    const MVMCPowerLanczosObservableBlockAccumulator *left,
    const MVMCPowerLanczosObservableBlockAccumulator *right) {
  return left != NULL && right != NULL && left->stage == right->stage &&
         left->request_count == right->request_count &&
         left->entries_per_request == right->entries_per_request &&
         left->entry_count == right->entry_count &&
         left->block_capacity == right->block_capacity &&
         left->block_length == right->block_length;
}

static MVMCKrylovStatus SimulateReduction(
    const MVMCPowerLanczosObservableBlockAccumulator *const *inputs,
    size_t input_count,
    MVMCPowerLanczosObservableBlockAccumulator *output,
    int commit) {
  const size_t block_count = inputs[0]->completed_block_count;
  size_t block;
  uint64_t completed_samples = 0;
  uint64_t discarded_samples = 0;
  for (block = 0; block < block_count; ++block) {
    uint64_t block_samples = 0;
    size_t rank;
    size_t entry;
    memset(output->transaction, 0,
           output->entry_count * sizeof(output->transaction[0]));
    for (rank = 0; rank < input_count; ++rank) {
      const CompensatedComplex *source =
          inputs[rank]->completed + block * output->entry_count;
      if (!CheckedAddU64(&block_samples,
                         inputs[rank]->completed_counts[block])) {
        return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      }
      for (entry = 0; entry < output->entry_count; ++entry) {
        if (!AddComplex(SumValue(&source[entry]),
                        &output->transaction[entry])) {
          return MVMC_KRYLOV_STATUS_NONFINITE;
        }
      }
    }
    if (!CheckedAddU64(&completed_samples, block_samples)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    if (commit) {
      memcpy(output->completed + block * output->entry_count,
             output->transaction,
             output->entry_count * sizeof(output->completed[0]));
      output->completed_counts[block] = block_samples;
    }
  }
  {
    size_t rank;
    for (rank = 0; rank < input_count; ++rank) {
      if (!CheckedAddU64(&discarded_samples,
                         inputs[rank]->discarded_partial_sample_count)) {
        return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      }
    }
  }
  if (commit) {
    output->completed_block_count = block_count;
    output->completed_sample_count = completed_samples;
    output->discarded_partial_sample_count = discarded_samples;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_observable_block_reduce_rank_ordered(
    const MVMCPowerLanczosObservableBlockAccumulator *const *inputs,
    size_t input_count,
    MVMCPowerLanczosObservableBlockAccumulator *output) {
  size_t rank;
  MVMCKrylovStatus status;
  if (inputs == NULL || input_count == 0 || output == NULL ||
      output->completed_block_count != 0 ||
      output->current_block_sample_count != 0 ||
      output->completed_sample_count != 0 ||
      output->discarded_partial_sample_count != 0 || inputs[0] == NULL ||
      !SameShape(inputs[0], output) ||
      inputs[0]->completed_block_count > output->block_capacity ||
      inputs[0]->current_block_sample_count != 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (rank = 0; rank < input_count; ++rank) {
    if (!SameShape(inputs[rank], output) ||
        inputs[rank]->completed_block_count !=
            inputs[0]->completed_block_count ||
        inputs[rank]->current_block_sample_count != 0) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  status = SimulateReduction(inputs, input_count, output, 0);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  return SimulateReduction(inputs, input_count, output, 1);
}
