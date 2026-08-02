/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "krylov_fock_reference.h"

#if !defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)
#error "krylov_fock_reference.c is Testing-only"
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef enum {
  MVMC_MEMO_EMPTY = 0,
  MVMC_MEMO_COMPUTING,
  MVMC_MEMO_COMPLETE
} MVMCMemoState;

typedef struct {
  MVMCMemoState state;
  int depth;
  uint64_t hash;
  double complex value;
} MVMCMemoEntry;

typedef struct {
  size_t key_index;
  double complex row_weight;
  size_t term_index;
} MVMCNeighbor;

struct MVMCKrylovWorkspace {
  size_t site_count;
  size_t orbital_count;
  size_t word_count;
  MVMCKrylovLimits limits;
  size_t allocated_bytes;
  size_t memo_slot_count;
  MVMCMemoEntry *memo_entries;
  uint64_t *memo_words;
  MVMCNeighbor *neighbors;
  uint64_t *neighbor_words;
  uint64_t *scratch_words;
  size_t memo_entry_count;
  size_t transition_used;
};

typedef struct {
  double real_sum;
  double real_compensation;
  double imag_sum;
  double imag_compensation;
} MVMCComplexAccumulator;

typedef struct {
  MVMCKrylovWorkspace *workspace;
  const MVMCKrylovFockModel *model;
  MVMCKrylovAmplitudeCallback amplitude;
  void *amplitude_context;
  MVMCKrylovStatistics *statistics;
  MVMCKrylovFailure *failure;
} MVMCEvaluation;

static int checked_add_size(size_t lhs, size_t rhs, size_t *result) {
  if (result == NULL || lhs > SIZE_MAX - rhs) return 0;
  *result = lhs + rhs;
  return 1;
}

static int checked_multiply_size(size_t lhs, size_t rhs, size_t *result) {
  if (result == NULL || (lhs != 0 && rhs > SIZE_MAX / lhs)) return 0;
  *result = lhs * rhs;
  return 1;
}

static int finite_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static double wall_seconds(void) {
  struct timespec now;
  if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return -1.0;
  return (double)now.tv_sec + 1.0e-9 * (double)now.tv_nsec;
}

static double elapsed_seconds(double start) {
  const double end = wall_seconds();
  return start >= 0.0 && end >= start ? end - start : 0.0;
}

static int exact_complex_zero(double complex value) {
  return creal(value) == 0.0 && cimag(value) == 0.0;
}

static int increment_u64(uint64_t *value) {
  if (value == NULL || *value == UINT64_MAX) return 0;
  ++(*value);
  return 1;
}

static int add_u64(uint64_t *value, uint64_t increment) {
  if (value == NULL || increment > UINT64_MAX - *value) return 0;
  *value += increment;
  return 1;
}

static uint64_t popcount_word(uint64_t value) {
  uint64_t count = 0;
  while (value != 0) {
    value &= value - 1;
    ++count;
  }
  return count;
}

static int configuration_bit(const uint64_t *words, size_t orbital) {
  return (int)((words[orbital / 64] >> (orbital % 64)) & UINT64_C(1));
}

static void configuration_set_bit(uint64_t *words, size_t orbital, int value) {
  const uint64_t mask = UINT64_C(1) << (orbital % 64);
  if (value) {
    words[orbital / 64] |= mask;
  } else {
    words[orbital / 64] &= ~mask;
  }
}

static int configuration_compare(const uint64_t *lhs, const uint64_t *rhs,
                                 size_t word_count) {
  size_t word;
  for (word = 0; word < word_count; ++word) {
    if (lhs[word] < rhs[word]) return -1;
    if (lhs[word] > rhs[word]) return 1;
  }
  return 0;
}

static uint64_t configuration_hash(const uint64_t *words, size_t word_count,
                                   int depth) {
  uint64_t hash = UINT64_C(1469598103934665603);
  size_t word;
  hash ^= (uint64_t)(unsigned int)depth;
  hash *= UINT64_C(1099511628211);
  for (word = 0; word < word_count; ++word) {
    uint64_t value = words[word];
    int byte;
    for (byte = 0; byte < 8; ++byte) {
      hash ^= value & UINT64_C(0xff);
      hash *= UINT64_C(1099511628211);
      value >>= 8;
    }
  }
  return hash;
}

static void accumulator_add_component(double value, double *sum,
                                      double *compensation) {
  const double updated = *sum + value;
  if (fabs(*sum) >= fabs(value)) {
    *compensation += (*sum - updated) + value;
  } else {
    *compensation += (value - updated) + *sum;
  }
  *sum = updated;
}

static int accumulator_add(MVMCComplexAccumulator *accumulator,
                           double complex value) {
  if (accumulator == NULL || !finite_complex(value)) return 0;
  accumulator_add_component(creal(value), &accumulator->real_sum,
                            &accumulator->real_compensation);
  accumulator_add_component(cimag(value), &accumulator->imag_sum,
                            &accumulator->imag_compensation);
  return isfinite(accumulator->real_sum) &&
         isfinite(accumulator->real_compensation) &&
         isfinite(accumulator->imag_sum) &&
         isfinite(accumulator->imag_compensation);
}

static double complex accumulator_value(
    const MVMCComplexAccumulator *accumulator) {
  return (accumulator->real_sum + accumulator->real_compensation) +
         I * (accumulator->imag_sum + accumulator->imag_compensation);
}

size_t mvmc_krylov_fock_word_count(size_t site_count) {
  size_t orbital_count;
  if (site_count == 0 || site_count > SIZE_MAX / 2) return 0;
  orbital_count = 2 * site_count;
  return orbital_count / 64 + (orbital_count % 64 != 0 ? 1 : 0);
}

static MVMCKrylovStatus validate_configuration(
    const MVMCKrylovFockModel *model, const uint64_t *words,
    size_t word_count) {
  size_t site;
  size_t up_count = 0;
  size_t down_count = 0;
  size_t required_words;
  size_t orbital_count;
  size_t used_bits;

  if (model == NULL || words == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  required_words = mvmc_krylov_fock_word_count(model->site_count);
  if (required_words == 0 || word_count != required_words) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  orbital_count = 2 * model->site_count;
  used_bits = orbital_count % 64;
  if (used_bits != 0) {
    const uint64_t used_mask = (UINT64_C(1) << used_bits) - UINT64_C(1);
    if ((words[word_count - 1] & ~used_mask) != 0) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  for (site = 0; site < model->site_count; ++site) {
    const int up = configuration_bit(words, site);
    const int down = configuration_bit(words, site + model->site_count);
    up_count += (size_t)up;
    down_count += (size_t)down;
    if (model->pure_spin && up + down != 1) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  if (up_count != model->up_electron_count ||
      down_count != model->down_electron_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus validate_model(const MVMCKrylovFockModel *model) {
  size_t term_index;
  if (model == NULL || mvmc_krylov_fock_word_count(model->site_count) == 0 ||
      model->up_electron_count > model->site_count ||
      model->down_electron_count > model->site_count ||
      (model->pure_spin != 0 && model->pure_spin != 1) ||
      (model->pure_spin &&
       (model->up_electron_count + model->down_electron_count !=
        model->site_count)) ||
      (model->term_count != 0 && model->terms == NULL) ||
      (model->operator_count != 0 && model->operators == NULL)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (!model->hermitian) return MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL;
  for (term_index = 0; term_index < model->term_count; ++term_index) {
    const MVMCKrylovHamiltonianTerm *term = &model->terms[term_index];
    size_t operator_index;
    if (!finite_complex(term->coefficient)) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    if ((term->operator_count != 2 && term->operator_count != 4) ||
        term->operator_offset > model->operator_count ||
        term->operator_count > model->operator_count - term->operator_offset) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    for (operator_index = 0; operator_index < term->operator_count;
         ++operator_index) {
      const MVMCKrylovFermionOperator *operator_item =
          &model->operators[term->operator_offset + operator_index];
      if ((operator_item->kind != MVMC_KRYLOV_FERMION_CREATE &&
           operator_item->kind != MVMC_KRYLOV_FERMION_ANNIHILATE) ||
          operator_item->orbital >= 2 * model->site_count) {
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      }
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_fock_validate(
    const MVMCKrylovFockModel *model, const uint64_t *configuration_words,
    size_t word_count) {
  MVMCKrylovStatus status = validate_model(model);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  return validate_configuration(model, configuration_words, word_count);
}

static int workspace_allocation_bytes(size_t word_count,
                                      const MVMCKrylovLimits *limits,
                                      size_t *memo_slot_count,
                                      size_t *allocated_bytes) {
  size_t slots;
  size_t bytes = sizeof(MVMCKrylovWorkspace);
  size_t part;

  if (word_count == 0 || limits == NULL || memo_slot_count == NULL ||
      allocated_bytes == NULL || limits->max_states == 0 ||
      limits->max_transitions == 0 || limits->max_amplitude_evaluations == 0 ||
      limits->max_bytes == 0 || limits->max_order < 0 ||
      limits->max_order > MVMC_KRYLOV_MAX_ORDER ||
      limits->max_states > (SIZE_MAX - 1) / 2) {
    return 0;
  }
  slots = 2 * limits->max_states + 1;
  if (!checked_multiply_size(slots, sizeof(MVMCMemoEntry), &part) ||
      !checked_add_size(bytes, part, &bytes) ||
      !checked_multiply_size(slots, word_count, &part) ||
      !checked_multiply_size(part, sizeof(uint64_t), &part) ||
      !checked_add_size(bytes, part, &bytes) ||
      !checked_multiply_size(limits->max_transitions,
                             sizeof(MVMCNeighbor), &part) ||
      !checked_add_size(bytes, part, &bytes) ||
      !checked_multiply_size(limits->max_transitions, word_count, &part) ||
      !checked_multiply_size(part, sizeof(uint64_t), &part) ||
      !checked_add_size(bytes, part, &bytes) ||
      !checked_multiply_size(word_count, sizeof(uint64_t), &part) ||
      !checked_add_size(bytes, part, &bytes)) {
    return 0;
  }
  *memo_slot_count = slots;
  *allocated_bytes = bytes;
  return 1;
}

MVMCKrylovWorkspace *mvmc_krylov_workspace_create(
    size_t site_count, const MVMCKrylovLimits *limits,
    MVMCKrylovStatus *status) {
  MVMCKrylovWorkspace *workspace = NULL;
  size_t memo_slot_count = 0;
  size_t allocated_bytes = 0;
  size_t word_count = mvmc_krylov_fock_word_count(site_count);

  if (status != NULL) *status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (!workspace_allocation_bytes(word_count, limits, &memo_slot_count,
                                  &allocated_bytes)) {
    return NULL;
  }
  if (allocated_bytes > limits->max_bytes) {
    if (status != NULL) *status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    return NULL;
  }
  workspace = (MVMCKrylovWorkspace *)calloc(1, sizeof(*workspace));
  if (workspace == NULL) {
    if (status != NULL) *status = MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    return NULL;
  }
  workspace->memo_entries = (MVMCMemoEntry *)calloc(
      memo_slot_count, sizeof(*workspace->memo_entries));
  workspace->memo_words = (uint64_t *)calloc(
      memo_slot_count * word_count, sizeof(*workspace->memo_words));
  workspace->neighbors = (MVMCNeighbor *)calloc(
      limits->max_transitions, sizeof(*workspace->neighbors));
  workspace->neighbor_words = (uint64_t *)calloc(
      limits->max_transitions * word_count,
      sizeof(*workspace->neighbor_words));
  workspace->scratch_words =
      (uint64_t *)calloc(word_count, sizeof(*workspace->scratch_words));
  if (workspace->memo_entries == NULL || workspace->memo_words == NULL ||
      workspace->neighbors == NULL || workspace->neighbor_words == NULL ||
      workspace->scratch_words == NULL) {
    mvmc_krylov_workspace_destroy(workspace);
    if (status != NULL) *status = MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    return NULL;
  }
  workspace->site_count = site_count;
  workspace->orbital_count = 2 * site_count;
  workspace->word_count = word_count;
  workspace->limits = *limits;
  workspace->allocated_bytes = allocated_bytes;
  workspace->memo_slot_count = memo_slot_count;
  if (status != NULL) *status = MVMC_KRYLOV_STATUS_OK;
  return workspace;
}

void mvmc_krylov_workspace_destroy(MVMCKrylovWorkspace *workspace) {
  if (workspace == NULL) return;
  free(workspace->scratch_words);
  free(workspace->neighbor_words);
  free(workspace->neighbors);
  free(workspace->memo_words);
  free(workspace->memo_entries);
  free(workspace);
}

size_t mvmc_krylov_workspace_bytes(const MVMCKrylovWorkspace *workspace) {
  return workspace == NULL ? 0 : workspace->allocated_bytes;
}

static void reset_result(MVMCKrylovResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  result->evaluated_order = -1;
  result->statistics.requested_order = -1;
  result->statistics.completed_order = -1;
  result->failure.status = MVMC_KRYLOV_STATUS_OK;
  result->failure.depth = -1;
}

static MVMCKrylovStatus record_failure(MVMCEvaluation *evaluation,
                                      MVMCKrylovStatus status, int depth,
                                      const uint64_t *words) {
  if (evaluation->failure->status == MVMC_KRYLOV_STATUS_OK) {
    evaluation->failure->status = status;
    evaluation->failure->depth = depth;
    evaluation->failure->configuration_hash =
        words == NULL ? 0 : configuration_hash(
                                  words, evaluation->workspace->word_count,
                                  depth < 0 ? 0 : depth);
  }
  return status;
}

static uint64_t count_occupied_before(const uint64_t *words, size_t orbital) {
  const size_t full_words = orbital / 64;
  const size_t partial_bits = orbital % 64;
  uint64_t count = 0;
  size_t word;
  for (word = 0; word < full_words; ++word) {
    count += popcount_word(words[word]);
  }
  if (partial_bits != 0) {
    const uint64_t mask = (UINT64_C(1) << partial_bits) - UINT64_C(1);
    count += popcount_word(words[full_words] & mask);
  }
  return count;
}

static MVMCKrylovStatus apply_term(
    const MVMCKrylovFockModel *model,
    const MVMCKrylovHamiltonianTerm *term, const uint64_t *input,
    size_t word_count, uint64_t *output, int *applied, int *fermion_sign) {
  size_t reverse_index;
  int sign = 1;
  memcpy(output, input, word_count * sizeof(*output));
  *applied = 0;
  *fermion_sign = 1;
  for (reverse_index = term->operator_count; reverse_index > 0;
       --reverse_index) {
    const MVMCKrylovFermionOperator *operator_item =
        &model->operators[term->operator_offset + reverse_index - 1];
    const int occupied = configuration_bit(output, operator_item->orbital);
    if ((operator_item->kind == MVMC_KRYLOV_FERMION_CREATE && occupied) ||
        (operator_item->kind == MVMC_KRYLOV_FERMION_ANNIHILATE &&
         !occupied)) {
      return MVMC_KRYLOV_STATUS_OK;
    }
    if ((count_occupied_before(output, operator_item->orbital) & UINT64_C(1)) !=
        0) {
      sign = -sign;
    }
    configuration_set_bit(
        output, operator_item->orbital,
        operator_item->kind == MVMC_KRYLOV_FERMION_CREATE);
  }
  if (validate_configuration(model, output, word_count) !=
      MVMC_KRYLOV_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *applied = 1;
  *fermion_sign = sign;
  return MVMC_KRYLOV_STATUS_OK;
}

static const uint64_t *neighbor_key(const MVMCKrylovWorkspace *workspace,
                                    const MVMCNeighbor *neighbor) {
  return workspace->neighbor_words + neighbor->key_index * workspace->word_count;
}

static uint64_t double_bits(double value) {
  uint64_t bits;
  memcpy(&bits, &value, sizeof(bits));
  return bits;
}

static int term_compare(const MVMCKrylovFockModel *model, size_t lhs_index,
                        size_t rhs_index) {
  const MVMCKrylovHamiltonianTerm *lhs = &model->terms[lhs_index];
  const MVMCKrylovHamiltonianTerm *rhs = &model->terms[rhs_index];
  size_t operator_index;
  if (lhs->operator_count < rhs->operator_count) return -1;
  if (lhs->operator_count > rhs->operator_count) return 1;
  for (operator_index = 0; operator_index < lhs->operator_count;
       ++operator_index) {
    const MVMCKrylovFermionOperator *lhs_operator =
        &model->operators[lhs->operator_offset + operator_index];
    const MVMCKrylovFermionOperator *rhs_operator =
        &model->operators[rhs->operator_offset + operator_index];
    if (lhs_operator->kind < rhs_operator->kind) return -1;
    if (lhs_operator->kind > rhs_operator->kind) return 1;
    if (lhs_operator->orbital < rhs_operator->orbital) return -1;
    if (lhs_operator->orbital > rhs_operator->orbital) return 1;
  }
  if (double_bits(creal(lhs->coefficient)) <
      double_bits(creal(rhs->coefficient))) {
    return -1;
  }
  if (double_bits(creal(lhs->coefficient)) >
      double_bits(creal(rhs->coefficient))) {
    return 1;
  }
  if (double_bits(cimag(lhs->coefficient)) <
      double_bits(cimag(rhs->coefficient))) {
    return -1;
  }
  if (double_bits(cimag(lhs->coefficient)) >
      double_bits(cimag(rhs->coefficient))) {
    return 1;
  }
  if (lhs->source_kind < rhs->source_kind) return -1;
  if (lhs->source_kind > rhs->source_kind) return 1;
  if (lhs->source_index < rhs->source_index) return -1;
  if (lhs->source_index > rhs->source_index) return 1;
  return 0;
}

static int neighbor_compare(const MVMCEvaluation *evaluation,
                            const MVMCNeighbor *lhs,
                            const MVMCNeighbor *rhs) {
  const int key_order = configuration_compare(
      neighbor_key(evaluation->workspace, lhs),
      neighbor_key(evaluation->workspace, rhs),
      evaluation->workspace->word_count);
  if (key_order != 0) return key_order;
  return term_compare(evaluation->model, lhs->term_index, rhs->term_index);
}

static void sort_neighbors(MVMCEvaluation *evaluation, size_t base,
                           size_t count) {
  size_t index;
  for (index = 1; index < count; ++index) {
    MVMCNeighbor value = evaluation->workspace->neighbors[base + index];
    size_t insertion = index;
    while (insertion > 0 &&
           neighbor_compare(evaluation, &value,
                            &evaluation->workspace->neighbors[
                                base + insertion - 1]) < 0) {
      evaluation->workspace->neighbors[base + insertion] =
          evaluation->workspace->neighbors[base + insertion - 1];
      --insertion;
    }
    evaluation->workspace->neighbors[base + insertion] = value;
  }
}

static MVMCKrylovStatus merge_neighbors(MVMCEvaluation *evaluation,
                                       size_t base, size_t raw_count,
                                       size_t *merged_count) {
  size_t read = 0;
  size_t write = 0;
  while (read < raw_count) {
    MVMCComplexAccumulator accumulator;
    MVMCNeighbor merged = evaluation->workspace->neighbors[base + read];
    size_t next = read + 1;
    memset(&accumulator, 0, sizeof(accumulator));
    if (!accumulator_add(&accumulator, merged.row_weight)) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    while (next < raw_count &&
           configuration_compare(
               neighbor_key(evaluation->workspace, &merged),
               neighbor_key(evaluation->workspace,
                            &evaluation->workspace->neighbors[base + next]),
               evaluation->workspace->word_count) == 0) {
      if (!increment_u64(
              &evaluation->statistics->merged_duplicate_transitions)) {
        return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      }
      if (!accumulator_add(
              &accumulator,
              evaluation->workspace->neighbors[base + next].row_weight)) {
        return MVMC_KRYLOV_STATUS_NONFINITE;
      }
      ++next;
    }
    merged.row_weight = accumulator_value(&accumulator);
    if (!finite_complex(merged.row_weight)) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    if (exact_complex_zero(merged.row_weight)) {
      if (!increment_u64(
              &evaluation->statistics->cancelled_zero_transitions)) {
        return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      }
    } else {
      evaluation->workspace->neighbors[base + write] = merged;
      ++write;
    }
    read = next;
  }
  *merged_count = write;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus memo_entry(MVMCEvaluation *evaluation, int depth,
                                  const uint64_t *words,
                                  MVMCMemoEntry **entry,
                                  const uint64_t **canonical_words,
                                  int *inserted) {
  MVMCKrylovWorkspace *workspace = evaluation->workspace;
  const uint64_t hash = configuration_hash(words, workspace->word_count, depth);
  const size_t start = (size_t)(hash % workspace->memo_slot_count);
  size_t probe;
  for (probe = 0; probe < workspace->memo_slot_count; ++probe) {
    const size_t slot = (start + probe) % workspace->memo_slot_count;
    MVMCMemoEntry *candidate = &workspace->memo_entries[slot];
    uint64_t *candidate_words =
        workspace->memo_words + slot * workspace->word_count;
    if (candidate->state == MVMC_MEMO_EMPTY) {
      if (workspace->memo_entry_count >= workspace->limits.max_states) {
        return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      }
      memcpy(candidate_words, words,
             workspace->word_count * sizeof(*candidate_words));
      candidate->state = MVMC_MEMO_COMPUTING;
      candidate->depth = depth;
      candidate->hash = hash;
      candidate->value = 0.0;
      ++workspace->memo_entry_count;
      if (workspace->memo_entry_count >
          evaluation->statistics->memo_entries_peak) {
        evaluation->statistics->memo_entries_peak =
            workspace->memo_entry_count;
      }
      *entry = candidate;
      *canonical_words = candidate_words;
      *inserted = 1;
      return MVMC_KRYLOV_STATUS_OK;
    }
    if (candidate->hash == hash && candidate->depth == depth &&
        configuration_compare(candidate_words, words, workspace->word_count) ==
            0) {
      *entry = candidate;
      *canonical_words = candidate_words;
      *inserted = 0;
      return MVMC_KRYLOV_STATUS_OK;
    }
  }
  return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
}

static MVMCKrylovStatus accumulate_amplitude_statistics(
    MVMCEvaluation *evaluation, const MVMCKrylovAmplitudeResult *amplitude) {
  if (!add_u64(&evaluation->statistics->regular_component_count,
               amplitude->regular_component_count) ||
      !add_u64(&evaluation->statistics->near_pivot_component_count,
               amplitude->near_pivot_component_count) ||
      !add_u64(&evaluation->statistics->singular_component_count,
               amplitude->singular_component_count) ||
      !add_u64(&evaluation->statistics->local_factorization_count,
               amplitude->local_factorization_count) ||
      !add_u64(&evaluation->statistics->global_factorization_count,
               amplitude->global_factorization_count) ||
      (amplitude->total_zero &&
       !increment_u64(&evaluation->statistics->total_zero_count))) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus evaluate_node(MVMCEvaluation *evaluation, int depth,
                                     size_t stack_level,
                                     const uint64_t *configuration,
                                     double complex *value) {
  MVMCKrylovWorkspace *workspace = evaluation->workspace;
  MVMCMemoEntry *entry = NULL;
  const uint64_t *canonical = NULL;
  int inserted = 0;
  MVMCKrylovStatus status;

  if (depth < 0 || depth > MVMC_KRYLOV_MAX_ORDER ||
      !increment_u64(&evaluation->statistics->recursion_calls[depth])) {
    return record_failure(evaluation, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                          depth, configuration);
  }
  if (stack_level + 1 > evaluation->statistics->frontier_peak) {
    evaluation->statistics->frontier_peak = stack_level + 1;
  }
  status = memo_entry(evaluation, depth, configuration, &entry, &canonical,
                      &inserted);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return record_failure(evaluation, status, depth, configuration);
  }
  if (!inserted) {
    if (entry->state == MVMC_MEMO_COMPUTING) {
      return record_failure(evaluation,
                            MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE,
                            depth, canonical);
    }
    if (entry->state != MVMC_MEMO_COMPLETE ||
        !increment_u64(&evaluation->statistics->memo_hits[depth]) ||
        (depth == 0 &&
         !increment_u64(&evaluation->statistics->terminal_cache_hits))) {
      return record_failure(evaluation, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                            depth, canonical);
    }
    *value = entry->value;
    return MVMC_KRYLOV_STATUS_OK;
  }
  if (!increment_u64(&evaluation->statistics->memo_misses[depth]) ||
      !increment_u64(&evaluation->statistics->unique_states[depth])) {
    return record_failure(evaluation, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                          depth, canonical);
  }

  if (depth == 0) {
    MVMCKrylovAmplitudeResult amplitude;
    double amplitude_start;
    if (evaluation->statistics->terminal_amplitude_requests >=
        workspace->limits.max_amplitude_evaluations) {
      return record_failure(evaluation, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                            depth, canonical);
    }
    if (!increment_u64(
            &evaluation->statistics->terminal_amplitude_requests)) {
      return record_failure(evaluation, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                            depth, canonical);
    }
    memset(&amplitude, 0, sizeof(amplitude));
    amplitude_start = wall_seconds();
    status = evaluation->amplitude(canonical, workspace->word_count,
                                   evaluation->amplitude_context, &amplitude);
    evaluation->statistics->amplitude_wall_seconds +=
        elapsed_seconds(amplitude_start);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      if (status >= MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
          status <= MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE) {
        return record_failure(evaluation, status, depth, canonical);
      }
      return record_failure(evaluation, MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE,
                            depth, canonical);
    }
    if (!finite_complex(amplitude.value) ||
        (amplitude.total_zero != 0 && amplitude.total_zero != 1) ||
        (amplitude.total_zero != exact_complex_zero(amplitude.value))) {
      return record_failure(evaluation, MVMC_KRYLOV_STATUS_NONFINITE, depth,
                            canonical);
    }
    status = accumulate_amplitude_statistics(evaluation, &amplitude);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      return record_failure(evaluation, status, depth, canonical);
    }
    entry->value = amplitude.value;
  } else {
    const size_t neighbor_base = workspace->transition_used;
    size_t raw_count;
    size_t merged_count = 0;
    size_t term_index;
    size_t neighbor_index;
    MVMCComplexAccumulator accumulator;
    double connectivity_start;
    memset(&accumulator, 0, sizeof(accumulator));

    connectivity_start = wall_seconds();
    for (term_index = 0; term_index < evaluation->model->term_count;
         ++term_index) {
      const MVMCKrylovHamiltonianTerm *term =
          &evaluation->model->terms[term_index];
      int applied;
      int fermion_sign;
      double complex row_weight;
      if (exact_complex_zero(term->coefficient)) continue;
      status = apply_term(evaluation->model, term, canonical,
                          workspace->word_count, workspace->scratch_words,
                          &applied, &fermion_sign);
      if (status != MVMC_KRYLOV_STATUS_OK) {
        return record_failure(evaluation, status, depth, canonical);
      }
      if (!applied) continue;
      if (workspace->transition_used >= workspace->limits.max_transitions) {
        return record_failure(evaluation, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                              depth, canonical);
      }
      row_weight = conj(term->coefficient * (double)fermion_sign);
      if (!finite_complex(row_weight)) {
        return record_failure(evaluation, MVMC_KRYLOV_STATUS_NONFINITE, depth,
                              canonical);
      }
      memcpy(workspace->neighbor_words +
                 workspace->transition_used * workspace->word_count,
             workspace->scratch_words,
             workspace->word_count * sizeof(*workspace->scratch_words));
      workspace->neighbors[workspace->transition_used].key_index =
          workspace->transition_used;
      workspace->neighbors[workspace->transition_used].row_weight = row_weight;
      workspace->neighbors[workspace->transition_used].term_index = term_index;
      ++workspace->transition_used;
      if (!increment_u64(&evaluation->statistics->raw_transitions)) {
        return record_failure(evaluation, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                              depth, canonical);
      }
    }
    raw_count = workspace->transition_used - neighbor_base;
    sort_neighbors(evaluation, neighbor_base, raw_count);
    status = merge_neighbors(evaluation, neighbor_base, raw_count,
                             &merged_count);
    evaluation->statistics->connectivity_wall_seconds +=
        elapsed_seconds(connectivity_start);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      return record_failure(evaluation, status, depth, canonical);
    }
    for (neighbor_index = 0; neighbor_index < merged_count;
         ++neighbor_index) {
      const MVMCNeighbor neighbor =
          workspace->neighbors[neighbor_base + neighbor_index];
      const uint64_t *neighbor_configuration =
          neighbor_key(workspace, &neighbor);
      double complex child_value;
      double complex contribution;
      status = evaluate_node(evaluation, depth - 1, stack_level + 1,
                             neighbor_configuration, &child_value);
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      contribution = neighbor.row_weight * child_value;
      if (!finite_complex(contribution) ||
          !accumulator_add(&accumulator, contribution)) {
        return record_failure(evaluation, MVMC_KRYLOV_STATUS_NONFINITE, depth,
                              canonical);
      }
    }
    entry->value = accumulator_value(&accumulator);
    if (!finite_complex(entry->value)) {
      return record_failure(evaluation, MVMC_KRYLOV_STATUS_NONFINITE, depth,
                            canonical);
    }
  }
  entry->state = MVMC_MEMO_COMPLETE;
  *value = entry->value;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_evaluate(
    MVMCKrylovWorkspace *workspace, const MVMCKrylovFockModel *model,
    const uint64_t *root_configuration_words, size_t root_word_count,
    MVMCKrylovAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovResult *result) {
  MVMCEvaluation evaluation;
  MVMCKrylovStatus status;
  int order;

  reset_result(result);
  if (workspace == NULL || model == NULL || root_configuration_words == NULL ||
      amplitude == NULL || result == NULL ||
      model->site_count != workspace->site_count ||
      root_word_count != workspace->word_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_fock_validate(model, root_configuration_words,
                                     root_word_count);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    result->status = status;
    result->failure.status = status;
    return status;
  }
  memset(workspace->memo_entries, 0,
         workspace->memo_slot_count * sizeof(*workspace->memo_entries));
  workspace->memo_entry_count = 0;
  workspace->transition_used = 0;
  result->statistics.root_evaluations = 1;
  result->statistics.requested_order = workspace->limits.max_order;
  result->statistics.completed_order = -1;
  result->statistics.workspace_bytes = workspace->allocated_bytes;

  evaluation.workspace = workspace;
  evaluation.model = model;
  evaluation.amplitude = amplitude;
  evaluation.amplitude_context = amplitude_context;
  evaluation.statistics = &result->statistics;
  evaluation.failure = &result->failure;

  for (order = 0; order <= workspace->limits.max_order; ++order) {
    const double depth_start = wall_seconds();
    status = evaluate_node(&evaluation, order, 0, root_configuration_words,
                           &result->value[order]);
    result->statistics.depth_wall_seconds[order] =
        elapsed_seconds(depth_start);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      memset(result->value, 0, sizeof(result->value));
      result->valid = 0;
      result->status = status;
      result->evaluated_order = -1;
      return status;
    }
    result->statistics.completed_order = order;
  }
  result->valid = 1;
  result->status = MVMC_KRYLOV_STATUS_OK;
  result->evaluated_order = workspace->limits.max_order;
  return MVMC_KRYLOV_STATUS_OK;
}

const char *mvmc_krylov_status_string(MVMCKrylovStatus status) {
  switch (status) {
    case MVMC_KRYLOV_STATUS_OK:
      return "ok";
    case MVMC_KRYLOV_STATUS_INVALID_ARGUMENT:
      return "invalid argument";
    case MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL:
      return "unsupported model";
    case MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL:
      return "non-Hermitian model";
    case MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE:
      return "allocation failure";
    case MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE:
      return "amplitude failure";
    case MVMC_KRYLOV_STATUS_NONFINITE:
      return "nonfinite value";
    case MVMC_KRYLOV_STATUS_RESOURCE_LIMIT:
      return "resource limit";
    case MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE:
      return "collective failure";
    case MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE:
      return "internal invariant failure";
  }
  return "unknown status";
}
