/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "bounded_krylov_engine.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "bounded_krylov_engine.c requires the power-Lanczos core"
#endif

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
  size_t key_index;
  size_t canonical_term_index;
  double complex row_weight;
} MVMCBoundedNeighbor;

typedef struct {
  uint64_t epoch;
  uint64_t access;
  uint64_t hash;
  int depth;
  MVMCScaledComplex value;
} MVMCBoundedCacheEntry;

struct MVMCKrylovBoundedPlan {
  MVMCKrylovFockModel model;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovHamiltonianTerm *terms;
  MVMCKrylovFermionOperator *operators;
  size_t *term_order;
  size_t word_count;
  size_t orbital_count;
  size_t max_row_transitions;
  size_t plan_bytes;
  size_t model_bytes;
  size_t cache_set_bytes;
  uint64_t plan_hash;
};

struct MVMCKrylovBoundedWorkspace {
  const MVMCKrylovBoundedPlan *plan;
  size_t allocated_bytes;
  size_t frame_bytes;
  size_t scratch_bytes;
  size_t cache_allocated_bytes;
  size_t cache_set_count;
  size_t cache_live_count;
  uint64_t cache_epoch;
  uint64_t cache_access;
  MVMCKrylovScaledAmplitudeCallback session_amplitude;
  void *session_amplitude_context;
  uint64_t session_amplitude_generation_hash;
  uint64_t session_root_evaluations;
  uint64_t session_epoch_full_clears_pending;
  double session_reset_seconds_pending;
  int session_active;
  int session_reset_pending;
  MVMCKrylovScaledAmplitudeCallback last_session_amplitude;
  void *last_session_amplitude_context;
  uint64_t last_session_amplitude_generation_hash;
  int last_session_binding_valid;
  MVMCBoundedNeighbor *neighbors;
  uint64_t *neighbor_words;
  MVMCScaledComplex *contributions;
  uint64_t *scratch_words;
  MVMCBoundedCacheEntry *cache_entries;
  uint64_t *cache_words;
};

typedef struct {
  MVMCKrylovBoundedWorkspace *workspace;
  MVMCKrylovScaledAmplitudeCallback amplitude;
  void *amplitude_context;
  MVMCKrylovBoundedStatistics *statistics;
  MVMCKrylovBoundedFailure *failure;
} MVMCBoundedEvaluation;

typedef struct {
  double real_sum;
  double real_compensation;
  double imag_sum;
  double imag_compensation;
} MVMCRawComplexAccumulator;

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

static int reserve_aligned(size_t *offset, size_t count, size_t item_size,
                           size_t alignment, size_t *start) {
  size_t aligned;
  size_t bytes;
  if (offset == NULL || start == NULL || alignment == 0 ||
      (alignment & (alignment - 1)) != 0 ||
      *offset > SIZE_MAX - (alignment - 1)) {
    return 0;
  }
  aligned = (*offset + alignment - 1) & ~(alignment - 1);
  if (!checked_multiply_size(count, item_size, &bytes) ||
      !checked_add_size(aligned, bytes, offset)) {
    return 0;
  }
  *start = aligned;
  return 1;
}

static int finite_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static int exact_complex_zero(double complex value) {
  return creal(value) == 0.0 && cimag(value) == 0.0;
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

static void configuration_set_bit(uint64_t *words, size_t orbital,
                                  int value) {
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

static void hash_byte(uint64_t *hash, unsigned char value) {
  *hash ^= (uint64_t)value;
  *hash *= UINT64_C(1099511628211);
}

static void hash_u64(uint64_t *hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    hash_byte(hash, (unsigned char)(value & UINT64_C(0xff)));
    value >>= 8;
  }
}

static uint64_t double_bits(double value) {
  uint64_t bits;
  memcpy(&bits, &value, sizeof(bits));
  return bits;
}

static int term_compare_model(const MVMCKrylovFockModel *model,
                              size_t lhs_index, size_t rhs_index) {
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

static MVMCKrylovStatus validate_model_only(
    const MVMCKrylovFockModel *model) {
  size_t term_index;
  size_t orbital_count;
  if (model == NULL || model->site_count == 0 ||
      model->site_count > SIZE_MAX / 2 ||
      model->up_electron_count > model->site_count ||
      model->down_electron_count > model->site_count ||
      (model->pure_spin != 0 && model->pure_spin != 1) ||
      (model->pure_spin &&
       model->up_electron_count + model->down_electron_count !=
           model->site_count) ||
      (model->term_count != 0 && model->terms == NULL) ||
      (model->operator_count != 0 && model->operators == NULL)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (!model->hermitian) return MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL;
  orbital_count = 2 * model->site_count;
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
          operator_item->orbital >= orbital_count) {
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      }
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static uint64_t plan_hash_value(const MVMCKrylovBoundedPlan *plan) {
  uint64_t hash = UINT64_C(1469598103934665603);
  size_t canonical_index;
  hash_u64(&hash, (uint64_t)plan->model.site_count);
  hash_u64(&hash, (uint64_t)plan->model.up_electron_count);
  hash_u64(&hash, (uint64_t)plan->model.down_electron_count);
  hash_u64(&hash, (uint64_t)(unsigned int)plan->model.pure_spin);
  hash_u64(&hash, (uint64_t)(unsigned int)plan->model.hermitian);
  hash_u64(&hash, (uint64_t)plan->model.term_count);
  hash_u64(&hash, plan->limits.amplitude_policy_hash);
  hash_u64(&hash, (uint64_t)plan->limits.cache_bytes);
  hash_u64(&hash, (uint64_t)plan->limits.max_row_transitions);
  hash_u64(&hash, (uint64_t)plan->limits.max_workspace_bytes);
  hash_u64(&hash, plan->limits.max_node_expansions);
  hash_u64(&hash, plan->limits.max_terminal_amplitude_calls);
  hash_u64(&hash, plan->limits.max_total_row_transitions);
  hash_u64(&hash, (uint64_t)(unsigned int)plan->limits.max_order);
  for (canonical_index = 0; canonical_index < plan->model.term_count;
       ++canonical_index) {
    const MVMCKrylovHamiltonianTerm *term =
        &plan->terms[plan->term_order[canonical_index]];
    size_t operator_index;
    hash_u64(&hash, double_bits(creal(term->coefficient)));
    hash_u64(&hash, double_bits(cimag(term->coefficient)));
    hash_u64(&hash, (uint64_t)term->operator_count);
    hash_u64(&hash, (uint64_t)(uint32_t)term->source_kind);
    hash_u64(&hash, (uint64_t)term->source_index);
    for (operator_index = 0; operator_index < term->operator_count;
         ++operator_index) {
      const MVMCKrylovFermionOperator *operator_item =
          &plan->operators[term->operator_offset + operator_index];
      hash_u64(&hash, (uint64_t)(unsigned int)operator_item->kind);
      hash_u64(&hash, (uint64_t)operator_item->orbital);
    }
  }
  return hash;
}

MVMCKrylovStatus mvmc_bounded_krylov_plan_create(
    const MVMCKrylovFockModel *model,
    const MVMCKrylovBoundedLimits *limits,
    MVMCKrylovBoundedPlan **plan) {
  MVMCKrylovBoundedPlan *created = NULL;
  MVMCKrylovStatus status;
  size_t term_bytes;
  size_t operator_bytes;
  size_t order_bytes;
  size_t plan_bytes = sizeof(*created);
  size_t max_row = 0;
  size_t cache_entry_bytes;
  size_t cache_word_bytes;
  size_t index;
  if (plan == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *plan = NULL;
  status = validate_model_only(model);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (limits == NULL || limits->max_order < 0 ||
      limits->max_order > MVMC_KRYLOV_MAX_ORDER ||
      limits->max_workspace_bytes == 0 ||
      limits->max_node_expansions == 0 ||
      limits->max_terminal_amplitude_calls == 0 ||
      limits->max_total_row_transitions == 0 ||
      !checked_multiply_size(model->term_count,
                             sizeof(MVMCKrylovHamiltonianTerm), &term_bytes) ||
      !checked_multiply_size(model->operator_count,
                             sizeof(MVMCKrylovFermionOperator),
                             &operator_bytes) ||
      !checked_multiply_size(model->term_count, sizeof(size_t), &order_bytes) ||
      !checked_add_size(plan_bytes, term_bytes, &plan_bytes) ||
      !checked_add_size(plan_bytes, operator_bytes, &plan_bytes) ||
      !checked_add_size(plan_bytes, order_bytes, &plan_bytes)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (index = 0; index < model->term_count; ++index) {
    if (!exact_complex_zero(model->terms[index].coefficient)) ++max_row;
  }
  if (max_row > limits->max_row_transitions) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  created = (MVMCKrylovBoundedPlan *)calloc(1, sizeof(*created));
  if (created == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  if (model->term_count != 0) {
    created->terms = (MVMCKrylovHamiltonianTerm *)malloc(term_bytes);
    created->term_order = (size_t *)malloc(order_bytes);
  }
  if (model->operator_count != 0) {
    created->operators = (MVMCKrylovFermionOperator *)malloc(operator_bytes);
  }
  if ((model->term_count != 0 &&
       (created->terms == NULL || created->term_order == NULL)) ||
      (model->operator_count != 0 && created->operators == NULL)) {
    mvmc_bounded_krylov_plan_destroy(created);
    return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
  if (term_bytes != 0) memcpy(created->terms, model->terms, term_bytes);
  if (operator_bytes != 0) {
    memcpy(created->operators, model->operators, operator_bytes);
  }
  created->model = *model;
  created->model.terms = created->terms;
  created->model.operators = created->operators;
  created->limits = *limits;
  created->word_count = mvmc_krylov_fock_word_count(model->site_count);
  created->orbital_count = 2 * model->site_count;
  created->max_row_transitions = max_row;
  created->plan_bytes = plan_bytes;
  created->model_bytes = term_bytes + operator_bytes;
  if (!checked_multiply_size(created->word_count, sizeof(uint64_t),
                             &cache_word_bytes) ||
      !checked_add_size(sizeof(MVMCBoundedCacheEntry), cache_word_bytes,
                        &cache_entry_bytes) ||
      !checked_multiply_size(MVMC_BOUNDED_KRYLOV_CACHE_WAYS,
                             cache_entry_bytes,
                             &created->cache_set_bytes)) {
    mvmc_bounded_krylov_plan_destroy(created);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (index = 0; index < model->term_count; ++index) {
    size_t insertion = index;
    const size_t value = index;
    while (insertion > 0 &&
           term_compare_model(&created->model, value,
                              created->term_order[insertion - 1]) < 0) {
      created->term_order[insertion] =
          created->term_order[insertion - 1];
      --insertion;
    }
    created->term_order[insertion] = value;
  }
  created->plan_hash = plan_hash_value(created);
  *plan = created;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_bounded_krylov_plan_destroy(MVMCKrylovBoundedPlan *plan) {
  if (plan == NULL) return;
  free(plan->term_order);
  free(plan->operators);
  free(plan->terms);
  free(plan);
}

size_t mvmc_bounded_krylov_plan_bytes(const MVMCKrylovBoundedPlan *plan) {
  return plan == NULL ? 0 : plan->plan_bytes;
}

size_t mvmc_bounded_krylov_plan_model_bytes(
    const MVMCKrylovBoundedPlan *plan) {
  return plan == NULL ? 0 : plan->model_bytes;
}

size_t mvmc_bounded_krylov_plan_max_row_transitions(
    const MVMCKrylovBoundedPlan *plan) {
  return plan == NULL ? 0 : plan->max_row_transitions;
}

size_t mvmc_bounded_krylov_cache_set_bytes(
    const MVMCKrylovBoundedPlan *plan) {
  return plan == NULL ? 0 : plan->cache_set_bytes;
}

uint64_t mvmc_bounded_krylov_plan_hash(
    const MVMCKrylovBoundedPlan *plan) {
  return plan == NULL ? 0 : plan->plan_hash;
}

MVMCKrylovStatus mvmc_bounded_krylov_workspace_create(
    const MVMCKrylovBoundedPlan *plan,
    MVMCKrylovBoundedWorkspace **workspace) {
  MVMCKrylovBoundedWorkspace *created;
  size_t offset = sizeof(MVMCKrylovBoundedWorkspace);
  size_t neighbor_offset;
  size_t neighbor_word_offset;
  size_t contribution_offset;
  size_t scratch_offset;
  size_t cache_entry_offset;
  size_t cache_word_offset;
  size_t frame_slots;
  size_t cache_slots;
  size_t frame_words;
  size_t cache_words;
  size_t frame_bytes;
  size_t scratch_bytes;
  size_t cache_entry_bytes;
  size_t cache_word_bytes;
  size_t cache_allocated_bytes;
  size_t set_count;
  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (plan == NULL || plan->cache_set_bytes == 0 ||
      !checked_multiply_size((size_t)plan->limits.max_order + 1,
                             plan->max_row_transitions, &frame_slots) ||
      !checked_multiply_size(frame_slots, plan->word_count, &frame_words)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  set_count = plan->limits.cache_bytes / plan->cache_set_bytes;
  if (!checked_multiply_size(set_count,
                             MVMC_BOUNDED_KRYLOV_CACHE_WAYS,
                             &cache_slots) ||
      !checked_multiply_size(cache_slots, plan->word_count, &cache_words) ||
      !checked_multiply_size(frame_slots, sizeof(MVMCBoundedNeighbor),
                             &frame_bytes) ||
      !checked_multiply_size(frame_words, sizeof(uint64_t), &scratch_bytes) ||
      !checked_add_size(frame_bytes, scratch_bytes, &frame_bytes) ||
      !checked_multiply_size(frame_slots, sizeof(MVMCScaledComplex),
                             &scratch_bytes) ||
      !checked_add_size(frame_bytes, scratch_bytes, &frame_bytes) ||
      !checked_multiply_size(plan->word_count, sizeof(uint64_t),
                             &scratch_bytes) ||
      !checked_multiply_size(cache_slots, sizeof(MVMCBoundedCacheEntry),
                             &cache_entry_bytes) ||
      !checked_multiply_size(cache_words, sizeof(uint64_t),
                             &cache_word_bytes) ||
      !checked_add_size(cache_entry_bytes, cache_word_bytes,
                        &cache_allocated_bytes) ||
      !reserve_aligned(&offset, frame_slots, sizeof(MVMCBoundedNeighbor),
                       _Alignof(MVMCBoundedNeighbor), &neighbor_offset) ||
      !reserve_aligned(&offset, frame_words, sizeof(uint64_t),
                       _Alignof(uint64_t), &neighbor_word_offset) ||
      !reserve_aligned(&offset, frame_slots, sizeof(MVMCScaledComplex),
                       _Alignof(MVMCScaledComplex), &contribution_offset) ||
      !reserve_aligned(&offset, plan->word_count, sizeof(uint64_t),
                       _Alignof(uint64_t), &scratch_offset) ||
      !reserve_aligned(&offset, cache_slots, sizeof(MVMCBoundedCacheEntry),
                       _Alignof(MVMCBoundedCacheEntry), &cache_entry_offset) ||
      !reserve_aligned(&offset, cache_words, sizeof(uint64_t),
                       _Alignof(uint64_t), &cache_word_offset)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (offset > plan->limits.max_workspace_bytes) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  created = (MVMCKrylovBoundedWorkspace *)calloc(1, offset);
  if (created == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  created->plan = plan;
  created->allocated_bytes = offset;
  created->frame_bytes = frame_bytes;
  created->scratch_bytes = scratch_bytes;
  created->cache_allocated_bytes = cache_allocated_bytes;
  created->cache_set_count = set_count;
  created->neighbors = (MVMCBoundedNeighbor *)((unsigned char *)created +
                                                neighbor_offset);
  created->neighbor_words = (uint64_t *)((unsigned char *)created +
                                         neighbor_word_offset);
  created->contributions = (MVMCScaledComplex *)((unsigned char *)created +
                                                  contribution_offset);
  created->scratch_words = (uint64_t *)((unsigned char *)created +
                                        scratch_offset);
  created->cache_entries = (MVMCBoundedCacheEntry *)((unsigned char *)created +
                                                      cache_entry_offset);
  created->cache_words = (uint64_t *)((unsigned char *)created +
                                      cache_word_offset);
  *workspace = created;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_bounded_krylov_workspace_destroy(
    MVMCKrylovBoundedWorkspace *workspace) {
  free(workspace);
}

size_t mvmc_bounded_krylov_workspace_bytes(
    const MVMCKrylovBoundedWorkspace *workspace) {
  return workspace == NULL ? 0 : workspace->allocated_bytes;
}

#if defined(MVMC_ENABLE_POWER_LANCZOS_TESTING_HOOKS)
MVMCKrylovStatus mvmc_bounded_krylov_testing_force_cache_counters(
    MVMCKrylovBoundedWorkspace *workspace,
    uint64_t epoch, uint64_t access_counter) {
  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  workspace->cache_epoch = epoch;
  workspace->cache_access = access_counter;
  return MVMC_KRYLOV_STATUS_OK;
}
#endif

static uint64_t count_occupied_before(const uint64_t *words,
                                      size_t orbital) {
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

static MVMCKrylovStatus validate_configuration_only(
    const MVMCKrylovFockModel *model, const uint64_t *words,
    size_t word_count) {
  size_t site;
  size_t up_count = 0;
  size_t down_count = 0;
  const size_t orbital_count = 2 * model->site_count;
  const size_t used_bits = orbital_count % 64;
  if (words == NULL || word_count != mvmc_krylov_fock_word_count(
                                        model->site_count)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (used_bits != 0) {
    const uint64_t used_mask =
        (UINT64_C(1) << used_bits) - UINT64_C(1);
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
    if ((count_occupied_before(output, operator_item->orbital) &
         UINT64_C(1)) != 0) {
      sign = -sign;
    }
    configuration_set_bit(
        output, operator_item->orbital,
        operator_item->kind == MVMC_KRYLOV_FERMION_CREATE);
  }
  if (validate_configuration_only(model, output, word_count) !=
      MVMC_KRYLOV_STATUS_OK) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *applied = 1;
  *fermion_sign = sign;
  return MVMC_KRYLOV_STATUS_OK;
}

static const uint64_t *neighbor_key(
    const MVMCKrylovBoundedWorkspace *workspace,
    const MVMCBoundedNeighbor *neighbor) {
  return workspace->neighbor_words +
         neighbor->key_index * workspace->plan->word_count;
}

static int neighbor_compare(const MVMCBoundedEvaluation *evaluation,
                            const MVMCBoundedNeighbor *lhs,
                            const MVMCBoundedNeighbor *rhs) {
  const int key_order = configuration_compare(
      neighbor_key(evaluation->workspace, lhs),
      neighbor_key(evaluation->workspace, rhs),
      evaluation->workspace->plan->word_count);
  if (key_order != 0) return key_order;
  if (lhs->canonical_term_index < rhs->canonical_term_index) return -1;
  if (lhs->canonical_term_index > rhs->canonical_term_index) return 1;
  return 0;
}

static void sort_neighbors(MVMCBoundedEvaluation *evaluation, size_t base,
                           size_t count) {
  size_t index;
  for (index = 1; index < count; ++index) {
    MVMCBoundedNeighbor value =
        evaluation->workspace->neighbors[base + index];
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

static void raw_accumulator_add_component(double value, double *sum,
                                          double *compensation) {
  const double updated = *sum + value;
  if (fabs(*sum) >= fabs(value)) {
    *compensation += (*sum - updated) + value;
  } else {
    *compensation += (value - updated) + *sum;
  }
  *sum = updated;
}

static int raw_accumulator_add(MVMCRawComplexAccumulator *accumulator,
                               double complex value) {
  if (accumulator == NULL || !finite_complex(value)) return 0;
  raw_accumulator_add_component(creal(value), &accumulator->real_sum,
                                &accumulator->real_compensation);
  raw_accumulator_add_component(cimag(value), &accumulator->imag_sum,
                                &accumulator->imag_compensation);
  return isfinite(accumulator->real_sum) &&
         isfinite(accumulator->real_compensation) &&
         isfinite(accumulator->imag_sum) &&
         isfinite(accumulator->imag_compensation);
}

static double complex raw_accumulator_value(
    const MVMCRawComplexAccumulator *accumulator) {
  return (accumulator->real_sum + accumulator->real_compensation) +
         I * (accumulator->imag_sum + accumulator->imag_compensation);
}

static MVMCKrylovStatus merge_neighbors(MVMCBoundedEvaluation *evaluation,
                                       size_t base, size_t raw_count,
                                       size_t *merged_count) {
  size_t read = 0;
  size_t write = 0;
  while (read < raw_count) {
    MVMCRawComplexAccumulator accumulator;
    MVMCBoundedNeighbor merged =
        evaluation->workspace->neighbors[base + read];
    size_t next = read + 1;
    memset(&accumulator, 0, sizeof(accumulator));
    if (!raw_accumulator_add(&accumulator, merged.row_weight)) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    while (next < raw_count &&
           configuration_compare(
               neighbor_key(evaluation->workspace, &merged),
               neighbor_key(evaluation->workspace,
                            &evaluation->workspace->neighbors[base + next]),
               evaluation->workspace->plan->word_count) == 0) {
      if (!increment_u64(
              &evaluation->statistics->merged_duplicate_transitions)) {
        return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      }
      if (!raw_accumulator_add(
              &accumulator,
              evaluation->workspace->neighbors[base + next].row_weight)) {
        return MVMC_KRYLOV_STATUS_NONFINITE;
      }
      ++next;
    }
    merged.row_weight = raw_accumulator_value(&accumulator);
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

static int decompose_exact_coefficient(double complex coefficient,
                                       MVMCScaledComplex *value) {
  const double scale = fmax(fabs(creal(coefficient)),
                            fabs(cimag(coefficient)));
  double complex scaled;
  double scaled_abs;
  if (value == NULL || !finite_complex(coefficient)) return 0;
  if (scale == 0.0) {
    return mvmc_scaled_complex_make_exact_zero(value) ==
           MVMC_PFAFFIAN_STATUS_OK;
  }
  scaled = (creal(coefficient) / scale) +
           I * (cimag(coefficient) / scale);
  scaled_abs = hypot(creal(scaled), cimag(scaled));
  if (!isfinite(scaled_abs) || scaled_abs == 0.0) return 0;
  return mvmc_scaled_complex_make_finite(
             scaled / scaled_abs, log(scale) + log(scaled_abs),
             -INFINITY, value) == MVMC_PFAFFIAN_STATUS_OK;
}

static void trace_event(MVMCBoundedEvaluation *evaluation, uint64_t tag,
                        uint64_t first, uint64_t second) {
  hash_u64(&evaluation->statistics->trace_hash, tag);
  hash_u64(&evaluation->statistics->trace_hash, first);
  hash_u64(&evaluation->statistics->trace_hash, second);
}

static MVMCKrylovStatus record_failure(
    MVMCBoundedEvaluation *evaluation, MVMCKrylovStatus status,
    int depth, const uint64_t *words) {
  if (evaluation->failure->status == MVMC_KRYLOV_STATUS_OK) {
    evaluation->failure->status = status;
    evaluation->failure->depth = depth;
    evaluation->failure->configuration_hash =
        words == NULL
            ? 0
            : configuration_hash(words, evaluation->workspace->plan->word_count,
                                 depth < 0 ? 0 : depth);
  }
  return status;
}

static uint64_t *cache_key(MVMCKrylovBoundedWorkspace *workspace,
                           size_t slot) {
  return workspace->cache_words + slot * workspace->plan->word_count;
}

static MVMCKrylovStatus next_cache_access(
    MVMCKrylovBoundedWorkspace *workspace, uint64_t *access) {
  if (workspace->cache_access == UINT64_MAX) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  ++workspace->cache_access;
  *access = workspace->cache_access;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus cache_lookup(MVMCBoundedEvaluation *evaluation,
                                    int depth, const uint64_t *words,
                                    MVMCScaledComplex *value, int *hit) {
  MVMCKrylovBoundedWorkspace *workspace = evaluation->workspace;
  const uint64_t hash = configuration_hash(
      words, workspace->plan->word_count, depth);
  const double start = wall_seconds();
  size_t way;
  *hit = 0;
  if (workspace->cache_set_count == 0) {
    evaluation->statistics->cache_wall_seconds += elapsed_seconds(start);
    return MVMC_KRYLOV_STATUS_OK;
  }
  {
    const size_t set = (size_t)(hash % workspace->cache_set_count);
    const size_t base = set * MVMC_BOUNDED_KRYLOV_CACHE_WAYS;
    for (way = 0; way < MVMC_BOUNDED_KRYLOV_CACHE_WAYS; ++way) {
      const size_t slot = base + way;
      MVMCBoundedCacheEntry *entry = &workspace->cache_entries[slot];
      if (entry->epoch == workspace->cache_epoch &&
          entry->hash == hash && entry->depth == depth &&
          configuration_compare(cache_key(workspace, slot), words,
                                workspace->plan->word_count) == 0) {
        MVMCKrylovStatus status =
            next_cache_access(workspace, &entry->access);
        if (status != MVMC_KRYLOV_STATUS_OK) {
          evaluation->statistics->cache_wall_seconds +=
              elapsed_seconds(start);
          return status;
        }
        *value = entry->value;
        *hit = 1;
        evaluation->statistics->cache_wall_seconds += elapsed_seconds(start);
        return MVMC_KRYLOV_STATUS_OK;
      }
    }
  }
  evaluation->statistics->cache_wall_seconds += elapsed_seconds(start);
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus cache_insert(MVMCBoundedEvaluation *evaluation,
                                    int depth, const uint64_t *words,
                                    const MVMCScaledComplex *value) {
  MVMCKrylovBoundedWorkspace *workspace = evaluation->workspace;
  const uint64_t hash = configuration_hash(
      words, workspace->plan->word_count, depth);
  const double start = wall_seconds();
  size_t victim = 0;
  size_t way;
  uint64_t access;
  int found_empty = 0;
  if (workspace->cache_set_count == 0) {
    evaluation->statistics->cache_wall_seconds += elapsed_seconds(start);
    return MVMC_KRYLOV_STATUS_OK;
  }
  {
    const size_t set = (size_t)(hash % workspace->cache_set_count);
    const size_t base = set * MVMC_BOUNDED_KRYLOV_CACHE_WAYS;
    victim = base;
    for (way = 0; way < MVMC_BOUNDED_KRYLOV_CACHE_WAYS; ++way) {
      const size_t slot = base + way;
      MVMCBoundedCacheEntry *entry = &workspace->cache_entries[slot];
      if (entry->epoch != workspace->cache_epoch) {
        victim = slot;
        found_empty = 1;
        break;
      }
      if (entry->access < workspace->cache_entries[victim].access) {
        victim = slot;
      }
    }
  }
  if (next_cache_access(workspace, &access) != MVMC_KRYLOV_STATUS_OK) {
    evaluation->statistics->cache_wall_seconds += elapsed_seconds(start);
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  if (!increment_u64(&evaluation->statistics->cache_insertions) ||
      (!found_empty &&
       !increment_u64(&evaluation->statistics->cache_evictions))) {
    evaluation->statistics->cache_wall_seconds += elapsed_seconds(start);
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  if (found_empty) {
    ++workspace->cache_live_count;
    if (workspace->cache_live_count >
        evaluation->statistics->cache_entries_peak) {
      evaluation->statistics->cache_entries_peak =
          workspace->cache_live_count;
    }
  }
  workspace->cache_entries[victim].epoch = workspace->cache_epoch;
  workspace->cache_entries[victim].access = access;
  workspace->cache_entries[victim].hash = hash;
  workspace->cache_entries[victim].depth = depth;
  workspace->cache_entries[victim].value = *value;
  memcpy(cache_key(workspace, victim), words,
         workspace->plan->word_count * sizeof(*words));
  evaluation->statistics->cache_wall_seconds += elapsed_seconds(start);
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus accumulate_amplitude_statistics(
    MVMCBoundedEvaluation *evaluation,
    const MVMCKrylovScaledAmplitudeResult *amplitude) {
  if (!add_u64(&evaluation->statistics->well_pivoted_component_count,
               amplitude->well_pivoted_component_count) ||
      !add_u64(&evaluation->statistics->near_pivot_component_count,
               amplitude->near_pivot_component_count) ||
      !add_u64(&evaluation->statistics->exact_zero_component_count,
               amplitude->exact_zero_component_count) ||
      !add_u64(&evaluation->statistics->numeric_zero_component_count,
               amplitude->numeric_zero_component_count) ||
      !add_u64(&evaluation->statistics->local_factorization_count,
               amplitude->local_factorization_count) ||
      !add_u64(&evaluation->statistics->global_factorization_count,
               amplitude->global_factorization_count)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus evaluate_node(MVMCBoundedEvaluation *evaluation,
                                     int depth,
                                     const uint64_t *configuration,
                                     MVMCScaledComplex *value) {
  MVMCKrylovBoundedWorkspace *workspace = evaluation->workspace;
  const MVMCKrylovBoundedPlan *plan = workspace->plan;
  MVMCKrylovStatus status;
  int cache_hit = 0;
  if (depth < 0 || depth > plan->limits.max_order ||
      !increment_u64(&evaluation->statistics->recursion_calls[depth])) {
    return record_failure(evaluation, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                          depth, configuration);
  }
  status = cache_lookup(evaluation, depth, configuration, value, &cache_hit);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return record_failure(evaluation, status, depth, configuration);
  }
  if (cache_hit) {
    if (!increment_u64(&evaluation->statistics->cache_hits[depth])) {
      return record_failure(evaluation, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                            depth, configuration);
    }
    trace_event(evaluation, UINT64_C(0x484954),
                (uint64_t)(unsigned int)depth,
                configuration_hash(configuration, plan->word_count, depth));
    return MVMC_KRYLOV_STATUS_OK;
  }
  if (!increment_u64(&evaluation->statistics->cache_misses[depth]) ||
      evaluation->statistics->node_expansions >=
          plan->limits.max_node_expansions ||
      !increment_u64(&evaluation->statistics->node_expansions)) {
    return record_failure(evaluation, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                          depth, configuration);
  }
  trace_event(evaluation, UINT64_C(0x4d495353),
              (uint64_t)(unsigned int)depth,
              configuration_hash(configuration, plan->word_count, depth));

  if (depth == 0) {
    MVMCKrylovScaledAmplitudeResult amplitude;
    const double start = wall_seconds();
    if (evaluation->statistics->terminal_amplitude_calls >=
            plan->limits.max_terminal_amplitude_calls ||
        !increment_u64(
            &evaluation->statistics->terminal_amplitude_calls)) {
      return record_failure(evaluation, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                            depth, configuration);
    }
    memset(&amplitude, 0, sizeof(amplitude));
    status = evaluation->amplitude(configuration, plan->word_count,
                                   evaluation->amplitude_context, &amplitude);
    evaluation->statistics->amplitude_wall_seconds += elapsed_seconds(start);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      if (status < MVMC_KRYLOV_STATUS_INVALID_ARGUMENT ||
          status > MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE) {
        status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
      }
      return record_failure(evaluation, status, depth, configuration);
    }
    if (!mvmc_scaled_complex_is_valid(&amplitude.value) ||
        amplitude.value.state == MVMC_SCALED_COMPLEX_NONFINITE) {
      return record_failure(evaluation, MVMC_KRYLOV_STATUS_NONFINITE,
                            depth, configuration);
    }
    status = accumulate_amplitude_statistics(evaluation, &amplitude);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      return record_failure(evaluation, status, depth, configuration);
    }
    trace_event(evaluation, UINT64_C(0x414d50),
                configuration_hash(configuration, plan->word_count, 0),
                (uint64_t)amplitude.value.state);
    *value = amplitude.value;
  } else {
    const size_t base = (size_t)depth * plan->max_row_transitions;
    size_t raw_count = 0;
    size_t merged_count = 0;
    size_t canonical_index;
    size_t neighbor_index;
    const double connectivity_start = wall_seconds();
    for (canonical_index = 0; canonical_index < plan->model.term_count;
         ++canonical_index) {
      const MVMCKrylovHamiltonianTerm *term =
          &plan->terms[plan->term_order[canonical_index]];
      MVMCBoundedNeighbor *neighbor;
      int applied;
      int fermion_sign;
      if (exact_complex_zero(term->coefficient)) continue;
      status = apply_term(&plan->model, term, configuration,
                          plan->word_count, workspace->scratch_words,
                          &applied, &fermion_sign);
      if (status != MVMC_KRYLOV_STATUS_OK) {
        return record_failure(evaluation, status, depth, configuration);
      }
      if (!applied) continue;
      if (raw_count >= plan->max_row_transitions ||
          evaluation->statistics->raw_row_transitions >=
              plan->limits.max_total_row_transitions ||
          !increment_u64(
              &evaluation->statistics->raw_row_transitions)) {
        return record_failure(evaluation, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                              depth, configuration);
      }
      neighbor = &workspace->neighbors[base + raw_count];
      neighbor->key_index = base + raw_count;
      neighbor->canonical_term_index = canonical_index;
      neighbor->row_weight = conj(term->coefficient * (double)fermion_sign);
      if (!finite_complex(neighbor->row_weight)) {
        return record_failure(evaluation, MVMC_KRYLOV_STATUS_NONFINITE,
                              depth, configuration);
      }
      memcpy(workspace->neighbor_words +
                 (base + raw_count) * plan->word_count,
             workspace->scratch_words,
             plan->word_count * sizeof(*workspace->scratch_words));
      ++raw_count;
    }
    if (raw_count > evaluation->statistics->row_transition_peak) {
      evaluation->statistics->row_transition_peak = raw_count;
    }
    sort_neighbors(evaluation, base, raw_count);
    status = merge_neighbors(evaluation, base, raw_count, &merged_count);
    evaluation->statistics->connectivity_wall_seconds +=
        elapsed_seconds(connectivity_start);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      return record_failure(evaluation, status, depth, configuration);
    }
    trace_event(evaluation, UINT64_C(0x524f57),
                configuration_hash(configuration, plan->word_count, depth),
                (uint64_t)merged_count);
    for (neighbor_index = 0; neighbor_index < merged_count;
         ++neighbor_index) {
      const MVMCBoundedNeighbor neighbor =
          workspace->neighbors[base + neighbor_index];
      const uint64_t *neighbor_configuration =
          neighbor_key(workspace, &neighbor);
      MVMCScaledComplex child;
      MVMCScaledComplex row_weight;
      status = evaluate_node(evaluation, depth - 1,
                             neighbor_configuration, &child);
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      if (!decompose_exact_coefficient(neighbor.row_weight, &row_weight) ||
          mvmc_scaled_complex_multiply(
              &row_weight, &child,
              &workspace->contributions[base + neighbor_index]) !=
              MVMC_PFAFFIAN_STATUS_OK ||
          workspace->contributions[base + neighbor_index].state ==
              MVMC_SCALED_COMPLEX_NONFINITE) {
        return record_failure(evaluation, MVMC_KRYLOV_STATUS_NONFINITE,
                              depth, configuration);
      }
    }
    if (merged_count == 0) {
      if (mvmc_scaled_complex_make_exact_zero(value) !=
          MVMC_PFAFFIAN_STATUS_OK) {
        return record_failure(evaluation,
                              MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE,
                              depth, configuration);
      }
    } else if (mvmc_scaled_complex_sum_ordered(
                   workspace->contributions + base, merged_count, value) !=
                   MVMC_PFAFFIAN_STATUS_OK ||
               value->state == MVMC_SCALED_COMPLEX_NONFINITE) {
      return record_failure(evaluation, MVMC_KRYLOV_STATUS_NONFINITE,
                            depth, configuration);
    }
  }
  if ((value->state == MVMC_SCALED_COMPLEX_EXACT_ZERO &&
       !increment_u64(
           &evaluation->statistics->computed_exact_zero_values)) ||
      (value->state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO &&
       !increment_u64(
           &evaluation->statistics->computed_numeric_zero_values))) {
    return record_failure(evaluation, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                          depth, configuration);
  }
  status = cache_insert(evaluation, depth, configuration, value);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return record_failure(evaluation, status, depth, configuration);
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static void initialize_result(MVMCKrylovBoundedResult *result) {
  int order;
  memset(result, 0, sizeof(*result));
  result->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  result->evaluated_order = -1;
  result->statistics.requested_order = -1;
  result->statistics.completed_order = -1;
  result->failure.status = MVMC_KRYLOV_STATUS_OK;
  result->failure.depth = -1;
  for (order = 0; order <= MVMC_KRYLOV_MAX_ORDER; ++order) {
    (void)mvmc_scaled_complex_make_exact_zero(&result->value[order]);
  }
}

static void invalidate_partial_values(MVMCKrylovBoundedResult *result,
                                      MVMCKrylovStatus status) {
  int order;
  for (order = 0; order <= MVMC_KRYLOV_MAX_ORDER; ++order) {
    (void)mvmc_scaled_complex_make_exact_zero(&result->value[order]);
  }
  result->valid = 0;
  result->status = status;
  result->evaluated_order = -1;
  result->statistics.completed_order = -1;
}

static void close_session(MVMCKrylovBoundedWorkspace *workspace) {
  if (workspace == NULL) return;
  workspace->session_active = 0;
  workspace->session_amplitude = NULL;
  workspace->session_amplitude_context = NULL;
  workspace->session_amplitude_generation_hash = 0;
  workspace->session_root_evaluations = 0;
  workspace->session_epoch_full_clears_pending = 0;
  workspace->session_reset_seconds_pending = 0.0;
  workspace->session_reset_pending = 0;
}

static MVMCKrylovStatus reset_cache(MVMCBoundedEvaluation *evaluation) {
  MVMCKrylovBoundedWorkspace *workspace = evaluation->workspace;
  const double start = wall_seconds();
  if (workspace->cache_epoch == UINT64_MAX) {
    size_t cache_slots;
    if (!checked_multiply_size(workspace->cache_set_count,
                               MVMC_BOUNDED_KRYLOV_CACHE_WAYS,
                               &cache_slots)) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    memset(workspace->cache_entries, 0,
           cache_slots * sizeof(*workspace->cache_entries));
    workspace->cache_epoch = 1;
    if (!increment_u64(
            &evaluation->statistics->cache_epoch_full_clears)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  } else {
    ++workspace->cache_epoch;
    if (workspace->cache_epoch == 0) workspace->cache_epoch = 1;
  }
  workspace->cache_live_count = 0;
  evaluation->statistics->reset_seconds += elapsed_seconds(start);
  return MVMC_KRYLOV_STATUS_OK;
}

static void initialize_evaluation_statistics(
    MVMCKrylovBoundedWorkspace *workspace,
    MVMCKrylovBoundedStatistics *statistics) {
  statistics->root_evaluations = 1;
  statistics->requested_order = workspace->plan->limits.max_order;
  statistics->plan_bytes = workspace->plan->plan_bytes;
  statistics->model_bytes = workspace->plan->model_bytes;
  statistics->workspace_bytes = workspace->allocated_bytes;
  statistics->frame_bytes = workspace->frame_bytes;
  statistics->scratch_bytes = workspace->scratch_bytes;
  statistics->cache_requested_bytes = workspace->plan->limits.cache_bytes;
  statistics->cache_allocated_bytes = workspace->cache_allocated_bytes;
  statistics->cache_set_count = workspace->cache_set_count;
}

static MVMCKrylovStatus evaluate_bound_amplitude(
    MVMCKrylovBoundedWorkspace *workspace,
    const uint64_t *root_configuration_words, size_t root_word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    int reset_before_evaluation, uint64_t amplitude_generation_hash,
    uint64_t session_root_evaluation, MVMCKrylovBoundedResult *result) {
  MVMCKrylovBoundedResult candidate;
  MVMCBoundedEvaluation evaluation;
  MVMCKrylovStatus status;
  double total_start;
  double evaluation_start;
  int order;
  size_t word;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  initialize_result(&candidate);
  if (workspace == NULL || workspace->plan == NULL ||
      root_configuration_words == NULL || amplitude == NULL ||
      root_word_count != workspace->plan->word_count) {
    *result = candidate;
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = validate_configuration_only(
      &workspace->plan->model, root_configuration_words, root_word_count);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    candidate.status = status;
    candidate.failure.status = status;
    *result = candidate;
    return status;
  }
  initialize_evaluation_statistics(workspace, &candidate.statistics);
  candidate.statistics.persistent_session_active =
      reset_before_evaluation ? 0 : 1;
  candidate.statistics.amplitude_generation_hash =
      amplitude_generation_hash;
  candidate.statistics.session_root_evaluation = session_root_evaluation;
  candidate.statistics.trace_hash = UINT64_C(1469598103934665603);
  hash_u64(&candidate.statistics.trace_hash, workspace->plan->plan_hash);
  for (word = 0; word < root_word_count; ++word) {
    hash_u64(&candidate.statistics.trace_hash,
             root_configuration_words[word]);
  }
  evaluation.workspace = workspace;
  evaluation.amplitude = amplitude;
  evaluation.amplitude_context = amplitude_context;
  evaluation.statistics = &candidate.statistics;
  evaluation.failure = &candidate.failure;
  total_start = wall_seconds();
  if (reset_before_evaluation) {
    status = reset_cache(&evaluation);
    candidate.statistics.cache_reset_performed = 1;
    if (status != MVMC_KRYLOV_STATUS_OK) {
      candidate.failure.status = status;
      candidate.failure.depth = -1;
      invalidate_partial_values(&candidate, status);
      candidate.statistics.total_seconds = elapsed_seconds(total_start);
      *result = candidate;
      return status;
    }
  } else if (workspace->session_reset_pending) {
    candidate.statistics.cache_reset_performed = 1;
    candidate.statistics.cache_epoch_full_clears =
        workspace->session_epoch_full_clears_pending;
    candidate.statistics.reset_seconds =
        workspace->session_reset_seconds_pending;
    workspace->session_epoch_full_clears_pending = 0;
    workspace->session_reset_seconds_pending = 0.0;
    workspace->session_reset_pending = 0;
  }
  candidate.statistics.cache_entries_resident_start =
      workspace->cache_live_count;
  candidate.statistics.cache_entries_peak = workspace->cache_live_count;
  evaluation_start = wall_seconds();
  for (order = 0; order <= workspace->plan->limits.max_order; ++order) {
    const double depth_start = wall_seconds();
    status = evaluate_node(&evaluation, order, root_configuration_words,
                           &candidate.value[order]);
    candidate.statistics.depth_wall_seconds[order] =
        elapsed_seconds(depth_start);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      candidate.statistics.evaluation_wall_seconds =
          elapsed_seconds(evaluation_start);
      invalidate_partial_values(&candidate, status);
      candidate.statistics.cache_entries_resident_end =
          workspace->cache_live_count;
      candidate.statistics.total_seconds = elapsed_seconds(total_start);
      *result = candidate;
      return status;
    }
    candidate.statistics.completed_order = order;
  }
  candidate.statistics.evaluation_wall_seconds =
      elapsed_seconds(evaluation_start);
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.evaluated_order = workspace->plan->limits.max_order;
  candidate.statistics.cache_entries_resident_end =
      workspace->cache_live_count;
  candidate.statistics.total_seconds = elapsed_seconds(total_start);
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_bounded_krylov_evaluate(
    MVMCKrylovBoundedWorkspace *workspace,
    const uint64_t *root_configuration_words, size_t root_word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovBoundedResult *result) {
  if (workspace != NULL && workspace->session_active) close_session(workspace);
  return evaluate_bound_amplitude(
      workspace, root_configuration_words, root_word_count, amplitude,
      amplitude_context, 1, 0, 0, result);
}

MVMCKrylovStatus mvmc_bounded_krylov_session_begin(
    MVMCKrylovBoundedWorkspace *workspace,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t amplitude_generation_hash) {
  MVMCKrylovBoundedStatistics statistics;
  MVMCKrylovBoundedFailure failure;
  MVMCBoundedEvaluation evaluation;
  MVMCKrylovStatus status;
  if (workspace == NULL || workspace->plan == NULL || amplitude == NULL ||
      amplitude_generation_hash == 0 || workspace->session_active ||
      (workspace->last_session_binding_valid &&
       amplitude_generation_hash ==
           workspace->last_session_amplitude_generation_hash &&
       (amplitude != workspace->last_session_amplitude ||
        amplitude_context != workspace->last_session_amplitude_context))) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&statistics, 0, sizeof(statistics));
  memset(&failure, 0, sizeof(failure));
  evaluation.workspace = workspace;
  evaluation.amplitude = amplitude;
  evaluation.amplitude_context = amplitude_context;
  evaluation.statistics = &statistics;
  evaluation.failure = &failure;
  status = reset_cache(&evaluation);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    close_session(workspace);
    return status;
  }
  workspace->session_amplitude = amplitude;
  workspace->session_amplitude_context = amplitude_context;
  workspace->session_amplitude_generation_hash = amplitude_generation_hash;
  workspace->session_root_evaluations = 0;
  workspace->session_epoch_full_clears_pending =
      statistics.cache_epoch_full_clears;
  workspace->session_reset_seconds_pending = statistics.reset_seconds;
  workspace->session_active = 1;
  workspace->session_reset_pending = 1;
  workspace->last_session_amplitude = amplitude;
  workspace->last_session_amplitude_context = amplitude_context;
  workspace->last_session_amplitude_generation_hash =
      amplitude_generation_hash;
  workspace->last_session_binding_valid = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

int mvmc_bounded_krylov_session_is_active(
    const MVMCKrylovBoundedWorkspace *workspace) {
  return workspace != NULL && workspace->session_active;
}

static MVMCKrylovStatus session_evaluate_checked(
    MVMCKrylovBoundedWorkspace *workspace,
    const uint64_t *root_configuration_words, size_t root_word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovBoundedResult *result) {
  MVMCKrylovStatus status;
  const int counter_overflow =
      workspace != NULL && workspace->session_active &&
      workspace->session_root_evaluations == UINT64_MAX;
  if (workspace == NULL || !workspace->session_active ||
      amplitude != workspace->session_amplitude ||
      amplitude_context != workspace->session_amplitude_context ||
      counter_overflow) {
    if (workspace != NULL && workspace->session_active) close_session(workspace);
    status = counter_overflow ? MVMC_KRYLOV_STATUS_RESOURCE_LIMIT
                              : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    if (result != NULL) {
      initialize_result(result);
      result->status = status;
      result->failure.status = status;
    }
    return status;
  }
  ++workspace->session_root_evaluations;
  status = evaluate_bound_amplitude(
      workspace, root_configuration_words, root_word_count,
      workspace->session_amplitude, workspace->session_amplitude_context, 0,
      workspace->session_amplitude_generation_hash,
      workspace->session_root_evaluations, result);
  if (status != MVMC_KRYLOV_STATUS_OK) close_session(workspace);
  return status;
}

MVMCKrylovStatus mvmc_bounded_krylov_session_evaluate(
    MVMCKrylovBoundedWorkspace *workspace,
    const uint64_t *root_configuration_words, size_t root_word_count,
    MVMCKrylovBoundedResult *result) {
  if (workspace == NULL || !workspace->session_active) {
    if (result != NULL) initialize_result(result);
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return session_evaluate_checked(
      workspace, root_configuration_words, root_word_count,
      workspace->session_amplitude, workspace->session_amplitude_context,
      result);
}

MVMCKrylovStatus mvmc_bounded_krylov_session_evaluate_bound(
    MVMCKrylovBoundedWorkspace *workspace,
    const uint64_t *root_configuration_words, size_t root_word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovBoundedResult *result) {
  return session_evaluate_checked(
      workspace, root_configuration_words, root_word_count, amplitude,
      amplitude_context, result);
}

MVMCKrylovStatus mvmc_bounded_krylov_session_end(
    MVMCKrylovBoundedWorkspace *workspace) {
  if (workspace == NULL || !workspace->session_active) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  close_session(workspace);
  return MVMC_KRYLOV_STATUS_OK;
}
