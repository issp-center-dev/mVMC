/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "krylov_fock_proposal.h"

#if !defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)
#error "krylov_fock_proposal.c is Testing-only"
#endif

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

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

static int binomial_limited(size_t n, size_t k, size_t limit,
                            size_t *value, int *exceeded);

static int finite_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
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

static void invalidate_uniform_proposal(
    MVMCKrylovStatus status, uint64_t *proposal_words,
    size_t proposal_word_count, MVMCKrylovFockUniformProposalResult *result) {
  if (proposal_words != NULL &&
      proposal_word_count <= SIZE_MAX / sizeof(proposal_words[0])) {
    memset(proposal_words, 0,
           proposal_word_count * sizeof(proposal_words[0]));
  }
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
}

static void invalidate_shell_proposal(
    MVMCKrylovStatus status, uint64_t *proposal_words,
    size_t proposal_word_count, MVMCKrylovFockShellProposalResult *result) {
  if (proposal_words != NULL &&
      proposal_word_count <= SIZE_MAX / sizeof(proposal_words[0])) {
    memset(proposal_words, 0,
           proposal_word_count * sizeof(proposal_words[0]));
  }
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
}

static MVMCKrylovStatus validate_model_only(
    const MVMCKrylovFockModel *model) {
  size_t term_index;
  const size_t word_count =
      model != NULL ? mvmc_krylov_fock_word_count(model->site_count) : 0;
  if (model == NULL || word_count == 0 ||
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

static MVMCKrylovStatus validate_configuration(
    const MVMCKrylovFockModel *model, const uint64_t *words,
    size_t word_count) {
  MVMCKrylovStatus status = validate_model_only(model);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  return mvmc_krylov_fock_validate(model, words, word_count);
}

static MVMCKrylovStatus compute_neighbor_count(
    const MVMCKrylovFockModel *model, size_t *neighbor_count) {
  size_t up_choices;
  size_t down_choices;
  size_t total;
  if (neighbor_count == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *neighbor_count = 0;
  if (model->pure_spin) {
    if (!checked_multiply_size(model->up_electron_count,
                               model->down_electron_count, &total)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    *neighbor_count = total;
    return MVMC_KRYLOV_STATUS_OK;
  }
  if (!checked_multiply_size(
          model->up_electron_count,
          model->site_count - model->up_electron_count, &up_choices) ||
      !checked_multiply_size(
          model->down_electron_count,
          model->site_count - model->down_electron_count, &down_choices) ||
      !checked_add_size(up_choices, down_choices, &total)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  *neighbor_count = total;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_fock_proposal_count_neighbors(
    const MVMCKrylovFockModel *model,
    const uint64_t *configuration_words, size_t word_count,
    size_t *neighbor_count) {
  MVMCKrylovStatus status;
  if (neighbor_count == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *neighbor_count = 0;
  status = validate_configuration(model, configuration_words, word_count);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  return compute_neighbor_count(model, neighbor_count);
}

static void copy_words(const uint64_t *source, uint64_t *destination,
                       size_t word_count) {
  memmove(destination, source, word_count * sizeof(source[0]));
}

static void apply_electronic_hop(const MVMCKrylovFockModel *model,
                                 const uint64_t *configuration_words,
                                 size_t word_count, int spin,
                                 size_t source_site, size_t target_site,
                                 uint64_t *proposal_words) {
  const size_t offset = spin != 0 ? model->site_count : 0;
  copy_words(configuration_words, proposal_words, word_count);
  configuration_set_bit(proposal_words, offset + source_site, 0);
  configuration_set_bit(proposal_words, offset + target_site, 1);
}

static void apply_pure_spin_exchange(const MVMCKrylovFockModel *model,
                                     const uint64_t *configuration_words,
                                     size_t word_count, size_t up_site,
                                     size_t down_site,
                                     uint64_t *proposal_words) {
  copy_words(configuration_words, proposal_words, word_count);
  configuration_set_bit(proposal_words, up_site, 0);
  configuration_set_bit(proposal_words, model->site_count + up_site, 1);
  configuration_set_bit(proposal_words, down_site, 1);
  configuration_set_bit(proposal_words, model->site_count + down_site, 0);
}

MVMCKrylovStatus mvmc_krylov_fock_proposal_select_neighbor(
    const MVMCKrylovFockModel *model,
    const uint64_t *configuration_words, size_t word_count,
    size_t selected_index,
    uint64_t *proposal_words, size_t proposal_word_count,
    size_t *neighbor_count) {
  MVMCKrylovStatus status;
  size_t count;
  size_t ordinal = 0;
  size_t first_site;

  if (proposal_words == NULL || proposal_word_count != word_count) {
    if (neighbor_count != NULL) *neighbor_count = 0;
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_fock_proposal_count_neighbors(
      model, configuration_words, word_count, &count);
  if (neighbor_count != NULL) *neighbor_count = count;
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (selected_index >= count) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;

  if (model->pure_spin) {
    for (first_site = 0; first_site < model->site_count; ++first_site) {
      size_t second_site;
      if (!configuration_bit(configuration_words, first_site)) continue;
      for (second_site = 0; second_site < model->site_count; ++second_site) {
        if (!configuration_bit(configuration_words,
                               model->site_count + second_site)) {
          continue;
        }
        if (ordinal == selected_index) {
          apply_pure_spin_exchange(model, configuration_words, word_count,
                                   first_site, second_site, proposal_words);
          return validate_configuration(model, proposal_words,
                                        proposal_word_count);
        }
        ++ordinal;
      }
    }
  } else {
    int spin;
    for (spin = 0; spin < 2; ++spin) {
      const size_t offset = spin != 0 ? model->site_count : 0;
      for (first_site = 0; first_site < model->site_count; ++first_site) {
        size_t second_site;
        if (!configuration_bit(configuration_words, offset + first_site)) {
          continue;
        }
        for (second_site = 0; second_site < model->site_count;
             ++second_site) {
          if (configuration_bit(configuration_words,
                                offset + second_site)) {
            continue;
          }
          if (ordinal == selected_index) {
            apply_electronic_hop(model, configuration_words, word_count,
                                 spin, first_site, second_site,
                                 proposal_words);
            return validate_configuration(model, proposal_words,
                                          proposal_word_count);
          }
          ++ordinal;
        }
      }
    }
  }
  return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

static int configurations_are_neighbors(
    const MVMCKrylovFockModel *model,
    const uint64_t *current_words, const uint64_t *proposal_words) {
  size_t site;
  size_t up_changes = 0;
  size_t down_changes = 0;
  size_t pure_changed_sites = 0;
  for (site = 0; site < model->site_count; ++site) {
    const int current_up = configuration_bit(current_words, site);
    const int proposal_up = configuration_bit(proposal_words, site);
    const int current_down =
        configuration_bit(current_words, model->site_count + site);
    const int proposal_down =
        configuration_bit(proposal_words, model->site_count + site);
    if (current_up != proposal_up) ++up_changes;
    if (current_down != proposal_down) ++down_changes;
    if (model->pure_spin &&
        (current_up != proposal_up || current_down != proposal_down)) {
      if (current_up == proposal_up || current_down == proposal_down) {
        return 0;
      }
      ++pure_changed_sites;
    }
  }
  if (model->pure_spin) {
    return up_changes == 2 && down_changes == 2 && pure_changed_sites == 2;
  }
  return (up_changes == 2 && down_changes == 0) ||
         (up_changes == 0 && down_changes == 2);
}

MVMCKrylovStatus mvmc_krylov_fock_proposal_log_ratio(
    const MVMCKrylovFockModel *model,
    const uint64_t *current_words, size_t current_word_count,
    const uint64_t *proposal_words, size_t proposal_word_count,
    double *log_proposal_ratio) {
  MVMCKrylovStatus status;
  size_t current_count = 0;
  size_t proposal_count = 0;
  if (log_proposal_ratio == NULL || current_word_count != proposal_word_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *log_proposal_ratio = 0.0;
  status = mvmc_krylov_fock_proposal_count_neighbors(
      model, current_words, current_word_count, &current_count);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  status = mvmc_krylov_fock_proposal_count_neighbors(
      model, proposal_words, proposal_word_count, &proposal_count);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (current_count == 0 || proposal_count == 0 ||
      !configurations_are_neighbors(model, current_words, proposal_words)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *log_proposal_ratio =
      log((double)current_count) - log((double)proposal_count);
  return isfinite(*log_proposal_ratio) ? MVMC_KRYLOV_STATUS_OK
                                      : MVMC_KRYLOV_STATUS_NONFINITE;
}

static MVMCKrylovStatus draw_uniform_subset(
    uint64_t *proposal_words, size_t orbital_offset, size_t site_count,
    size_t occupied_count, MVMCKrylovFockProposalBoundedDraw draw_bounded,
    void *draw_context, size_t *draw_count) {
  size_t site;
  size_t needed = occupied_count;
  if (proposal_words == NULL || draw_bounded == NULL || draw_count == NULL ||
      occupied_count > site_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (site = 0; site < site_count; ++site) {
    const size_t remaining = site_count - site;
    size_t selected = 0;
    MVMCKrylovStatus status;
    if (needed == 0) break;
    if (needed == remaining) {
      for (; site < site_count; ++site) {
        configuration_set_bit(proposal_words, orbital_offset + site, 1);
      }
      needed = 0;
      break;
    }
    status = draw_bounded(draw_context, remaining, &selected);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    if (selected >= remaining) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    if (*draw_count == SIZE_MAX) return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    ++(*draw_count);
    if (selected < needed) {
      configuration_set_bit(proposal_words, orbital_offset + site, 1);
      --needed;
    }
  }
  return needed == 0 ? MVMC_KRYLOV_STATUS_OK
                     : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

MVMCKrylovStatus mvmc_krylov_fock_proposal_draw_uniform_sector(
    const MVMCKrylovFockModel *model,
    MVMCKrylovFockProposalBoundedDraw draw_bounded, void *draw_context,
    uint64_t *proposal_words, size_t proposal_word_count,
    MVMCKrylovFockUniformProposalResult *result) {
  MVMCKrylovStatus status;
  size_t expected_word_count;
  size_t draw_count = 0;
  size_t site;
  MVMCKrylovFockUniformProposalResult candidate;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_uniform_proposal(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                              proposal_words, proposal_word_count, result);
  status = validate_model_only(model);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_uniform_proposal(status, proposal_words, proposal_word_count,
                                result);
    return status;
  }
  expected_word_count = mvmc_krylov_fock_word_count(model->site_count);
  if (draw_bounded == NULL || proposal_words == NULL ||
      proposal_word_count != expected_word_count ||
      proposal_word_count > SIZE_MAX / sizeof(proposal_words[0])) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(proposal_words, 0,
         proposal_word_count * sizeof(proposal_words[0]));

  status = draw_uniform_subset(
      proposal_words, 0, model->site_count, model->up_electron_count,
      draw_bounded, draw_context, &draw_count);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_uniform_proposal(status, proposal_words, proposal_word_count,
                                result);
    return status;
  }
  if (model->pure_spin) {
    for (site = 0; site < model->site_count; ++site) {
      configuration_set_bit(
          proposal_words, model->site_count + site,
          !configuration_bit(proposal_words, site));
    }
  } else {
    status = draw_uniform_subset(
        proposal_words, model->site_count, model->site_count,
        model->down_electron_count, draw_bounded, draw_context, &draw_count);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      invalidate_uniform_proposal(status, proposal_words,
                                  proposal_word_count, result);
      return status;
    }
  }
  status = validate_configuration(model, proposal_words, proposal_word_count);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_uniform_proposal(status, proposal_words, proposal_word_count,
                                result);
    return status;
  }

  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.word_count = proposal_word_count;
  candidate.active_orbital_count = 2 * model->site_count;
  candidate.draw_count = draw_count;
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

static size_t min_size(size_t lhs, size_t rhs) {
  return lhs < rhs ? lhs : rhs;
}

MVMCKrylovStatus mvmc_krylov_fock_proposal_max_shell_distance(
    const MVMCKrylovFockModel *model, size_t *max_distance) {
  MVMCKrylovStatus status;
  size_t up_max;
  size_t down_max;
  if (max_distance == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *max_distance = 0;
  status = validate_model_only(model);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (model->pure_spin) {
    *max_distance = min_size(model->up_electron_count,
                             model->down_electron_count);
    return MVMC_KRYLOV_STATUS_OK;
  }
  up_max = min_size(model->up_electron_count,
                    model->site_count - model->up_electron_count);
  down_max = min_size(model->down_electron_count,
                      model->site_count - model->down_electron_count);
  if (!checked_add_size(up_max, down_max, max_distance)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_fock_proposal_resolve_shell_distance(
    const MVMCKrylovFockModel *model, size_t fraction_numerator,
    size_t fraction_denominator, size_t *max_distance, size_t *distance) {
  MVMCKrylovStatus status;
  size_t maximum = 0;
  size_t product;
  size_t resolved;
  size_t remainder;
  const size_t rounding_threshold =
      fraction_denominator / 2 + fraction_denominator % 2;
  if (max_distance == NULL || distance == NULL || max_distance == distance ||
      fraction_numerator == 0 ||
      fraction_denominator == 0 ||
      fraction_numerator > fraction_denominator) {
    if (max_distance != NULL) *max_distance = 0;
    if (distance != NULL) *distance = 0;
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *max_distance = 0;
  *distance = 0;
  status = mvmc_krylov_fock_proposal_max_shell_distance(model, &maximum);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (maximum == 0) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (!checked_multiply_size(maximum, fraction_numerator, &product)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  resolved = product / fraction_denominator;
  remainder = product % fraction_denominator;
  if (remainder >= rounding_threshold) {
    if (!checked_add_size(resolved, 1, &resolved)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  }
  if (resolved == 0) resolved = 1;
  if (resolved > maximum) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  *max_distance = maximum;
  *distance = resolved;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus checked_binomial_size(size_t n, size_t k,
                                              size_t *value) {
  int exceeded = 0;
  if (value == NULL || k > n) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *value = 0;
  if (!binomial_limited(n, k, SIZE_MAX, value, &exceeded)) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  return exceeded ? MVMC_KRYLOV_STATUS_RESOURCE_LIMIT
                  : MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus electronic_shell_split_count(
    const MVMCKrylovFockModel *model, size_t distance,
    size_t up_distance, size_t *split_count) {
  MVMCKrylovStatus status;
  size_t down_distance;
  size_t factors[4];
  size_t product = 1;
  size_t factor_index;
  if (model == NULL || split_count == NULL || up_distance > distance) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *split_count = 0;
  down_distance = distance - up_distance;
  if (up_distance > model->up_electron_count ||
      up_distance > model->site_count - model->up_electron_count ||
      down_distance > model->down_electron_count ||
      down_distance > model->site_count - model->down_electron_count) {
    return MVMC_KRYLOV_STATUS_OK;
  }
  status = checked_binomial_size(model->up_electron_count, up_distance,
                                 &factors[0]);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = checked_binomial_size(
        model->site_count - model->up_electron_count, up_distance,
        &factors[1]);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = checked_binomial_size(model->down_electron_count,
                                   down_distance, &factors[2]);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = checked_binomial_size(
        model->site_count - model->down_electron_count, down_distance,
        &factors[3]);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  for (factor_index = 0; factor_index < 4; ++factor_index) {
    if (!checked_multiply_size(product, factors[factor_index], &product)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  }
  *split_count = product;
  return MVMC_KRYLOV_STATUS_OK;
}

static void electronic_shell_split_interval(
    const MVMCKrylovFockModel *model, size_t distance,
    size_t *minimum_up_distance, size_t *maximum_up_distance) {
  const size_t up_max =
      min_size(model->up_electron_count,
               model->site_count - model->up_electron_count);
  const size_t down_max =
      min_size(model->down_electron_count,
               model->site_count - model->down_electron_count);
  *minimum_up_distance = distance > down_max ? distance - down_max : 0;
  *maximum_up_distance = min_size(distance, up_max);
}

static MVMCKrylovStatus compute_shell_count(
    const MVMCKrylovFockModel *model, size_t distance,
    size_t *shell_count) {
  MVMCKrylovStatus status;
  size_t maximum = 0;
  if (shell_count == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *shell_count = 0;
  status = mvmc_krylov_fock_proposal_max_shell_distance(model, &maximum);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (distance == 0 || distance > maximum) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (model->pure_spin) {
    size_t up_count;
    size_t down_count;
    status = checked_binomial_size(model->up_electron_count, distance,
                                   &up_count);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = checked_binomial_size(model->down_electron_count, distance,
                                     &down_count);
    }
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    if (!checked_multiply_size(up_count, down_count, shell_count)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  } else {
    size_t up_distance;
    size_t minimum_up_distance;
    size_t maximum_up_distance;
    size_t total = 0;
    electronic_shell_split_interval(model, distance, &minimum_up_distance,
                                    &maximum_up_distance);
    for (up_distance = minimum_up_distance;
         up_distance <= maximum_up_distance; ++up_distance) {
      size_t split_count;
      status = electronic_shell_split_count(model, distance, up_distance,
                                            &split_count);
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      if (!checked_add_size(total, split_count, &total)) {
        return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      }
      if (up_distance == maximum_up_distance) break;
      if (up_distance == SIZE_MAX) {
        return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      }
    }
    *shell_count = total;
  }
  return *shell_count != 0 ? MVMC_KRYLOV_STATUS_OK
                           : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

MVMCKrylovStatus mvmc_krylov_fock_proposal_count_shell(
    const MVMCKrylovFockModel *model, size_t distance,
    size_t *shell_count) {
  return compute_shell_count(model, distance, shell_count);
}

static int word_buffers_overlap(const uint64_t *first, size_t first_count,
                                const uint64_t *second,
                                size_t second_count) {
  uintptr_t first_start;
  uintptr_t second_start;
  size_t first_bytes;
  size_t second_bytes;
  if (first == NULL || second == NULL ||
      !checked_multiply_size(first_count, sizeof(first[0]), &first_bytes) ||
      !checked_multiply_size(second_count, sizeof(second[0]), &second_bytes)) {
    return 1;
  }
  first_start = (uintptr_t)first;
  second_start = (uintptr_t)second;
  if (first_start > UINTPTR_MAX - first_bytes ||
      second_start > UINTPTR_MAX - second_bytes) {
    return 1;
  }
  return first_start < second_start + second_bytes &&
         second_start < first_start + first_bytes;
}

static MVMCKrylovStatus draw_eligible_subset(
    const uint64_t *configuration_words, uint64_t *proposal_words,
    size_t orbital_offset, size_t site_count, int eligible_value,
    size_t eligible_count, size_t selected_count, int proposal_value,
    MVMCKrylovFockProposalBoundedDraw draw_bounded, void *draw_context,
    size_t *draw_count) {
  size_t site;
  size_t remaining = eligible_count;
  size_t needed = selected_count;
  if (configuration_words == NULL || proposal_words == NULL ||
      draw_bounded == NULL || draw_count == NULL ||
      selected_count > eligible_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (site = 0; site < site_count; ++site) {
    size_t selected;
    MVMCKrylovStatus status;
    if (configuration_bit(configuration_words, orbital_offset + site) !=
        eligible_value) {
      continue;
    }
    if (remaining == 0) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    if (needed == 0) break;
    if (needed == remaining) {
      configuration_set_bit(proposal_words, orbital_offset + site,
                            proposal_value);
      --needed;
      --remaining;
      continue;
    }
    status = draw_bounded(draw_context, remaining, &selected);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    if (selected >= remaining) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    if (*draw_count == SIZE_MAX) return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    ++(*draw_count);
    if (selected < needed) {
      configuration_set_bit(proposal_words, orbital_offset + site,
                            proposal_value);
      --needed;
    }
    --remaining;
  }
  return needed == 0 ? MVMC_KRYLOV_STATUS_OK
                     : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

static MVMCKrylovStatus draw_electronic_split(
    const MVMCKrylovFockModel *model, size_t distance, size_t shell_count,
    MVMCKrylovFockProposalBoundedDraw draw_bounded, void *draw_context,
    size_t *up_distance, size_t *draw_count) {
  MVMCKrylovStatus status;
  size_t selected;
  size_t cumulative = 0;
  size_t split;
  size_t minimum_split;
  size_t maximum_split;
  if (up_distance == NULL || draw_count == NULL || shell_count == 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = draw_bounded(draw_context, shell_count, &selected);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (selected >= shell_count) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (*draw_count == SIZE_MAX) return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  ++(*draw_count);
  electronic_shell_split_interval(model, distance, &minimum_split,
                                  &maximum_split);
  for (split = minimum_split; split <= maximum_split; ++split) {
    size_t split_count;
    status = electronic_shell_split_count(model, distance, split,
                                          &split_count);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    if (!checked_add_size(cumulative, split_count, &cumulative)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    if (selected < cumulative) {
      *up_distance = split;
      return MVMC_KRYLOV_STATUS_OK;
    }
    if (split == maximum_split) break;
    if (split == SIZE_MAX) return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

static void configuration_replacement_distances(
    const MVMCKrylovFockModel *model, const uint64_t *current_words,
    const uint64_t *proposal_words, size_t *up_distance,
    size_t *down_distance) {
  size_t site;
  *up_distance = 0;
  *down_distance = 0;
  for (site = 0; site < model->site_count; ++site) {
    if (configuration_bit(current_words, site) &&
        !configuration_bit(proposal_words, site)) {
      ++(*up_distance);
    }
    if (configuration_bit(current_words, model->site_count + site) &&
        !configuration_bit(proposal_words, model->site_count + site)) {
      ++(*down_distance);
    }
  }
}

MVMCKrylovStatus mvmc_krylov_fock_proposal_draw_shell(
    const MVMCKrylovFockModel *model,
    const uint64_t *configuration_words, size_t word_count,
    size_t distance, MVMCKrylovFockProposalBoundedDraw draw_bounded,
    void *draw_context, uint64_t *proposal_words,
    size_t proposal_word_count, MVMCKrylovFockShellProposalResult *result) {
  MVMCKrylovStatus status;
  size_t expected_word_count;
  size_t maximum = 0;
  size_t shell_count = 0;
  size_t up_distance = 0;
  size_t down_distance = 0;
  size_t draw_count = 0;
  size_t site;
  MVMCKrylovFockShellProposalResult candidate;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(result, 0, sizeof(*result));
  result->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  status = validate_model_only(model);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    result->status = status;
    return status;
  }
  expected_word_count = mvmc_krylov_fock_word_count(model->site_count);
  if (configuration_words == NULL || proposal_words == NULL ||
      draw_bounded == NULL || word_count != expected_word_count ||
      proposal_word_count != expected_word_count) {
    if (proposal_words != NULL &&
        proposal_word_count == expected_word_count) {
      invalidate_shell_proposal(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                                proposal_words, proposal_word_count, result);
    }
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (word_buffers_overlap(configuration_words, word_count, proposal_words,
                           proposal_word_count)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  invalidate_shell_proposal(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                            proposal_words, proposal_word_count, result);
  status = validate_configuration(model, configuration_words, word_count);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_shell_proposal(status, proposal_words, proposal_word_count,
                              result);
    return status;
  }
  status = mvmc_krylov_fock_proposal_max_shell_distance(model, &maximum);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = compute_shell_count(model, distance, &shell_count);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_shell_proposal(status, proposal_words, proposal_word_count,
                              result);
    return status;
  }

  copy_words(configuration_words, proposal_words, word_count);
  if (model->pure_spin) {
    up_distance = distance;
    down_distance = distance;
    status = draw_eligible_subset(
        configuration_words, proposal_words, 0, model->site_count, 1,
        model->up_electron_count, distance, 0, draw_bounded, draw_context,
        &draw_count);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = draw_eligible_subset(
          configuration_words, proposal_words, 0, model->site_count, 0,
          model->down_electron_count, distance, 1, draw_bounded, draw_context,
          &draw_count);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      for (site = 0; site < model->site_count; ++site) {
        configuration_set_bit(
            proposal_words, model->site_count + site,
            !configuration_bit(proposal_words, site));
      }
    }
  } else {
    status = draw_electronic_split(
        model, distance, shell_count, draw_bounded, draw_context,
        &up_distance, &draw_count);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      down_distance = distance - up_distance;
      status = draw_eligible_subset(
          configuration_words, proposal_words, 0, model->site_count, 1,
          model->up_electron_count, up_distance, 0, draw_bounded,
          draw_context, &draw_count);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = draw_eligible_subset(
          configuration_words, proposal_words, 0, model->site_count, 0,
          model->site_count - model->up_electron_count, up_distance, 1,
          draw_bounded, draw_context, &draw_count);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = draw_eligible_subset(
          configuration_words, proposal_words, model->site_count,
          model->site_count, 1, model->down_electron_count, down_distance, 0,
          draw_bounded, draw_context, &draw_count);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = draw_eligible_subset(
          configuration_words, proposal_words, model->site_count,
          model->site_count, 0,
          model->site_count - model->down_electron_count, down_distance, 1,
          draw_bounded, draw_context, &draw_count);
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = validate_configuration(model, proposal_words,
                                    proposal_word_count);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    size_t actual_up_distance;
    size_t actual_down_distance;
    configuration_replacement_distances(
        model, configuration_words, proposal_words, &actual_up_distance,
        &actual_down_distance);
    if (actual_up_distance != up_distance ||
        actual_down_distance != down_distance) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_shell_proposal(status, proposal_words, proposal_word_count,
                              result);
    return status;
  }

  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.word_count = word_count;
  candidate.max_distance = maximum;
  candidate.distance = distance;
  candidate.up_distance = up_distance;
  candidate.down_distance = down_distance;
  candidate.shell_count = shell_count;
  candidate.draw_count = draw_count;
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

static size_t gcd_size(size_t lhs, size_t rhs) {
  while (rhs != 0) {
    const size_t next = lhs % rhs;
    lhs = rhs;
    rhs = next;
  }
  return lhs;
}

static int binomial_limited(size_t n, size_t k, size_t limit,
                            size_t *value, int *exceeded) {
  size_t result = 1;
  size_t index;
  if (value == NULL || exceeded == NULL || k > n) return 0;
  if (k > n - k) k = n - k;
  *exceeded = 0;
  for (index = 1; index <= k; ++index) {
    size_t numerator = n - k + index;
    size_t denominator = index;
    size_t divisor = gcd_size(numerator, denominator);
    numerator /= divisor;
    denominator /= divisor;
    divisor = gcd_size(result, denominator);
    result /= divisor;
    denominator /= divisor;
    if (denominator != 1 ||
        (numerator != 0 && result > limit / numerator)) {
      *exceeded = 1;
      *value = result;
      return 1;
    }
    result *= numerator;
  }
  *value = result;
  if (result > limit) *exceeded = 1;
  return 1;
}

static MVMCKrylovStatus exact_sector_count(
    const MVMCKrylovFockModel *model, size_t max_states,
    size_t *sector_count) {
  size_t up_sector_count;
  int exceeded;
  MVMCKrylovStatus status = validate_model_only(model);
  if (sector_count == NULL || max_states == 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *sector_count = 0;
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (model->site_count > 32) return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  if (!binomial_limited(model->site_count, model->up_electron_count,
                        max_states, &up_sector_count, &exceeded)) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (exceeded) return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  if (model->pure_spin) {
    *sector_count = up_sector_count;
    return MVMC_KRYLOV_STATUS_OK;
  }
  {
    size_t down_sector_count;
    if (!binomial_limited(model->site_count, model->down_electron_count,
                          max_states, &down_sector_count, &exceeded)) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    if (exceeded ||
        (down_sector_count != 0 &&
         up_sector_count > max_states / down_sector_count)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    *sector_count = up_sector_count * down_sector_count;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static uint64_t lowest_combination(size_t bit_count) {
  return bit_count == 0 ? UINT64_C(0)
                        : ((UINT64_C(1) << bit_count) - UINT64_C(1));
}

static int next_combination(uint64_t *mask, size_t site_count) {
  const uint64_t limit = UINT64_C(1) << site_count;
  const uint64_t value = *mask;
  const uint64_t smallest = value & (~value + UINT64_C(1));
  const uint64_t ripple = value + smallest;
  uint64_t next;
  if (value == 0 || ripple >= limit) return 0;
  next = ripple | (((value ^ ripple) >> 2) / smallest);
  if (next >= limit) return 0;
  *mask = next;
  return 1;
}

static MVMCKrylovStatus fill_sector_states(
    const MVMCKrylovFockModel *model, uint64_t *states,
    size_t state_capacity, size_t *state_count) {
  const uint64_t all_sites = model->site_count == 64
                                 ? UINT64_MAX
                                 : ((UINT64_C(1) << model->site_count) -
                                    UINT64_C(1));
  uint64_t up_mask = lowest_combination(model->up_electron_count);
  size_t index = 0;

  if (states == NULL || state_count == NULL || model->site_count > 32) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (model->up_electron_count == 0) up_mask = UINT64_C(0);
  for (;;) {
    if (model->pure_spin) {
      const uint64_t down_mask = all_sites ^ up_mask;
      if (index >= state_capacity) return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      states[index++] = up_mask | (down_mask << model->site_count);
    } else {
      uint64_t down_mask = lowest_combination(model->down_electron_count);
      if (model->down_electron_count == 0) down_mask = UINT64_C(0);
      for (;;) {
        if (index >= state_capacity) return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
        states[index++] = up_mask | (down_mask << model->site_count);
        if (!next_combination(&down_mask, model->site_count)) break;
      }
    }
    if (!next_combination(&up_mask, model->site_count)) break;
  }
  *state_count = index;
  return MVMC_KRYLOV_STATUS_OK;
}

static size_t find_state_index(const uint64_t *states, size_t state_count,
                               uint64_t state, int *found) {
  size_t index;
  for (index = 0; index < state_count; ++index) {
    if (states[index] == state) {
      *found = 1;
      return index;
    }
  }
  *found = 0;
  return 0;
}

static void invalidate_connectivity(
    MVMCKrylovStatus status, MVMCKrylovFockProposalConnectivity *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
}

MVMCKrylovStatus mvmc_krylov_fock_proposal_check_connectivity(
    const MVMCKrylovFockModel *model, size_t max_states,
    MVMCKrylovFockProposalConnectivity *result) {
  MVMCKrylovStatus status;
  size_t sector_count = 0;
  uint64_t *states = NULL;
  unsigned char *visited = NULL;
  size_t *queue = NULL;
  size_t generated_count = 0;
  size_t visited_count = 0;
  size_t queue_head = 0;
  size_t queue_tail = 0;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_connectivity(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  status = exact_sector_count(model, max_states, &sector_count);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_connectivity(status, result);
    return status;
  }
  states = (uint64_t *)calloc(sector_count, sizeof(*states));
  visited = (unsigned char *)calloc(sector_count, sizeof(*visited));
  queue = (size_t *)calloc(sector_count, sizeof(*queue));
  if (states == NULL || visited == NULL || queue == NULL) {
    free(queue);
    free(visited);
    free(states);
    invalidate_connectivity(MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE, result);
    return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
  status = fill_sector_states(model, states, sector_count, &generated_count);
  if (status == MVMC_KRYLOV_STATUS_OK && generated_count != sector_count) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    free(queue);
    free(visited);
    free(states);
    invalidate_connectivity(status, result);
    return status;
  }

  visited[0] = 1;
  queue[queue_tail++] = 0;
  visited_count = 1;
  while (queue_head < queue_tail) {
    const size_t state_index = queue[queue_head++];
    const uint64_t current = states[state_index];
    size_t neighbor_count = 0;
    size_t neighbor_index;
    status = mvmc_krylov_fock_proposal_count_neighbors(
        model, &current, 1, &neighbor_count);
    if (status != MVMC_KRYLOV_STATUS_OK) break;
    for (neighbor_index = 0; neighbor_index < neighbor_count;
         ++neighbor_index) {
      uint64_t neighbor = 0;
      size_t returned_count = 0;
      int found = 0;
      size_t next_index;
      status = mvmc_krylov_fock_proposal_select_neighbor(
          model, &current, 1, neighbor_index, &neighbor, 1,
          &returned_count);
      if (status != MVMC_KRYLOV_STATUS_OK ||
          returned_count != neighbor_count) {
        status = status == MVMC_KRYLOV_STATUS_OK
                     ? MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE
                     : status;
        break;
      }
      next_index =
          find_state_index(states, sector_count, neighbor, &found);
      if (!found) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
        break;
      }
      if (!visited[next_index]) {
        visited[next_index] = 1;
        queue[queue_tail++] = next_index;
        ++visited_count;
      }
    }
    if (status != MVMC_KRYLOV_STATUS_OK) break;
  }

  free(queue);
  free(visited);
  free(states);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_connectivity(status, result);
    return status;
  }
  result->valid = 1;
  result->status = MVMC_KRYLOV_STATUS_OK;
  result->sector_count = sector_count;
  result->visited_count = visited_count;
  result->connected = visited_count == sector_count;
  return MVMC_KRYLOV_STATUS_OK;
}
