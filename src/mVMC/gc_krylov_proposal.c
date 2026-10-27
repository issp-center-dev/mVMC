/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "gc_krylov_proposal.h"

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                         \
    defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include <limits.h>
#include <stdlib.h>
#include <string.h>

static int gc_domain(const MVMCKrylovFockModel *model) {
  return model->domain == MVMC_KRYLOV_DOMAIN_GC_BALANCED ||
         model->domain == MVMC_KRYLOV_DOMAIN_GC_EVEN;
}

static MVMCKrylovStatus validate_target(const MVMCKrylovFockModel *model,
                                        MVMCKrylovFockProposalBoundedDraw draw,
                                        const uint64_t *out,
                                        size_t word_count) {
  if (model == NULL || draw == NULL || out == NULL || !gc_domain(model) ||
      model->up_electron_count != 0 || model->down_electron_count != 0 ||
      model->pure_spin != 0 ||
      mvmc_krylov_fock_word_count(model->site_count) == 0 ||
      word_count != mvmc_krylov_fock_word_count(model->site_count)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

/* One bounded integer; a value outside [0, bound) is an invariant failure. */
static MVMCKrylovStatus bounded(MVMCKrylovFockProposalBoundedDraw draw,
                                void *context, size_t bound, size_t *value) {
  MVMCKrylovStatus status;
  *value = bound;
  status = draw(context, bound, value);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  return *value < bound ? MVMC_KRYLOV_STATUS_OK
                        : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

static int get_bit(const uint64_t *words, size_t orbital) {
  return (int)((words[orbital / 64] >> (orbital % 64)) & UINT64_C(1));
}

static void set_bit(uint64_t *words, size_t orbital) {
  words[orbital / 64] |= UINT64_C(1) << (orbital % 64);
}

static void flip_bit(uint64_t *words, size_t orbital) {
  words[orbital / 64] ^= UINT64_C(1) << (orbital % 64);
}

MVMCKrylovStatus mvmc_gc_krylov_draw_uniform(
    const MVMCKrylovFockModel *model,
    MVMCKrylovFockProposalBoundedDraw draw, void *context,
    uint64_t *out, size_t word_count) {
  MVMCKrylovStatus status = validate_target(model, draw, out, word_count);
  size_t nsite, orbital;
  int parity = 0;
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  nsite = model->site_count;
  memset(out, 0, word_count * sizeof(*out));
  if (model->domain == MVMC_KRYLOV_DOMAIN_GC_EVEN) {
    /* 2L-1 free bits; the last orbital restores even parity. */
    for (orbital = 0; orbital + 1 < 2 * nsite; ++orbital) {
      size_t value;
      status = bounded(draw, context, 2, &value);
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      if (value) {
        set_bit(out, orbital);
        parity ^= 1;
      }
    }
    if (parity) set_bit(out, 2 * nsite - 1);
    return MVMC_KRYLOV_STATUS_OK;
  }
  /* Labels U_0..U_{L-1}, Dbar_0..Dbar_{L-1}; select L of them. */
  {
    size_t label;
    size_t remaining = nsite;
    for (label = 0; label < 2 * nsite; ++label) {
      size_t value;
      int taken;
      status = bounded(draw, context, 2 * nsite - label, &value);
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      taken = value < remaining;
      remaining -= (size_t)taken;
      if (label < nsite) {
        if (taken) set_bit(out, label);
      } else if (!taken) {
        set_bit(out, label);
      }
    }
    return remaining == 0 ? MVMC_KRYLOV_STATUS_OK
                          : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
}

MVMCKrylovStatus mvmc_gc_krylov_draw_initial(
    const MVMCKrylovFockModel *model, int initial_n,
    MVMCKrylovFockProposalBoundedDraw draw, void *context,
    uint64_t *out, size_t word_count) {
  MVMCKrylovStatus status = validate_target(model, draw, out, word_count);
  size_t nsite, item;
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  nsite = model->site_count;
  if (initial_n < 0 || (initial_n % 2) != 0 ||
      (size_t)initial_n > 2 * nsite) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(out, 0, word_count * sizeof(*out));
  /* Sequential subset selection writes one bit per item, no scratch. */
  if (model->domain == MVMC_KRYLOV_DOMAIN_GC_EVEN) {
    size_t remaining = (size_t)initial_n;
    for (item = 0; item < 2 * nsite; ++item) {
      size_t value;
      status = bounded(draw, context, 2 * nsite - item, &value);
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      if (value < remaining) {
        set_bit(out, item);
        --remaining;
      }
    }
    return remaining == 0 ? MVMC_KRYLOV_STATUS_OK
                          : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  {
    int spin;
    for (spin = 0; spin < 2; ++spin) {
      size_t remaining = (size_t)initial_n / 2;
      for (item = 0; item < nsite; ++item) {
        size_t value;
        status = bounded(draw, context, nsite - item, &value);
        if (status != MVMC_KRYLOV_STATUS_OK) return status;
        if (value < remaining) {
          set_bit(out, item + (size_t)spin * nsite);
          --remaining;
        }
      }
      if (remaining != 0) return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

/* Second element of an ordered pair of distinct items: skip the first. */
static MVMCKrylovStatus distinct_pair(
    MVMCKrylovFockProposalBoundedDraw draw, void *context, size_t count,
    size_t *first, size_t *second) {
  size_t value;
  MVMCKrylovStatus status = bounded(draw, context, count, first);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  status = bounded(draw, context, count - 1, &value);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  *second = value < *first ? value : value + 1;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_gc_krylov_draw_local(
    const MVMCKrylovFockModel *model, const uint64_t *current,
    MVMCKrylovFockProposalBoundedDraw draw, void *context,
    uint64_t *out, size_t word_count) {
  MVMCKrylovStatus status = validate_target(model, draw, out, word_count);
  size_t nsite;
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (current == NULL || current == out) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  {
    const uintptr_t a = (uintptr_t)current;
    const uintptr_t b = (uintptr_t)out;
    const uintptr_t bytes = (uintptr_t)(word_count * sizeof(*out));
    if (a < b + bytes && b < a + bytes) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  status = mvmc_krylov_fock_validate(model, current, word_count);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  nsite = model->site_count;
  memcpy(out, current, word_count * sizeof(*out));
  if (model->domain == MVMC_KRYLOV_DOMAIN_GC_EVEN) {
    size_t first, second;
    status = distinct_pair(draw, context, 2 * nsite, &first, &second);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    flip_bit(out, first);
    flip_bit(out, second);
    return MVMC_KRYLOV_STATUS_OK;
  }
  {
    size_t kind;
    status = bounded(draw, context, 2, &kind);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    if (kind == 0) {
      /* Pair flip (up_i, down_j) when both are empty or both occupied. */
      size_t up_site, down_site;
      status = bounded(draw, context, nsite, &up_site);
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = bounded(draw, context, nsite, &down_site);
      }
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      if (get_bit(current, up_site) == get_bit(current, down_site + nsite)) {
        flip_bit(out, up_site);
        flip_bit(out, down_site + nsite);
      }
      return MVMC_KRYLOV_STATUS_OK;
    }
    /* Same-spin exchange of two distinct sites with different occupation. */
    {
      size_t spin, first, second;
      status = bounded(draw, context, 2, &spin);
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      if (nsite < 2) return MVMC_KRYLOV_STATUS_OK;
      status = distinct_pair(draw, context, nsite, &first, &second);
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      first += spin * nsite;
      second += spin * nsite;
      if (get_bit(current, first) != get_bit(current, second)) {
        flip_bit(out, first);
        flip_bit(out, second);
      }
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_gc_krylov_draw_two_pairs(
    const MVMCKrylovFockModel *model, const uint64_t *current,
    MVMCKrylovFockProposalBoundedDraw draw, void *context,
    uint64_t *out, size_t word_count) {
  MVMCKrylovStatus status = validate_target(model, draw, out, word_count);
  size_t up0, up1, down0, down1, nsite;
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (model->domain != MVMC_KRYLOV_DOMAIN_GC_BALANCED || current == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  {
    const uintptr_t a = (uintptr_t)current, b = (uintptr_t)out;
    const uintptr_t bytes = (uintptr_t)(word_count * sizeof(*out));
    if (a < b + bytes && b < a + bytes)
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_fock_validate(model, current, word_count);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  memcpy(out, current, word_count * sizeof(*out));
  nsite = model->site_count;
  if (nsite < 2) return MVMC_KRYLOV_STATUS_OK;
  status = distinct_pair(draw, context, nsite, &up0, &up1);
  if (status == MVMC_KRYLOV_STATUS_OK)
    status = distinct_pair(draw, context, nsite, &down0, &down1);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  down0 += nsite;
  down1 += nsite;
  if (get_bit(current, up0) == get_bit(current, up1) &&
      get_bit(current, up0) == get_bit(current, down0) &&
      get_bit(current, up0) == get_bit(current, down1)) {
    flip_bit(out, up0);
    flip_bit(out, up1);
    flip_bit(out, down0);
    flip_bit(out, down1);
  }
  return MVMC_KRYLOV_STATUS_OK;
}

/* The index-th site (ascending) whose up/down occupation is (up, down). */
static size_t nth_site(const uint64_t *words, size_t nsite, int up, int down,
                       size_t index) {
  size_t site;
  for (site = 0; site < nsite; ++site) {
    if (get_bit(words, site) == up && get_bit(words, site + nsite) == down) {
      if (index == 0) return site;
      --index;
    }
  }
  return nsite;
}

MVMCKrylovStatus mvmc_gc_krylov_draw_spin_swap(
    const MVMCKrylovFockModel *model, const uint64_t *current,
    MVMCKrylovFockProposalBoundedDraw draw, void *context,
    uint64_t *out, size_t word_count) {
  MVMCKrylovStatus status = validate_target(model, draw, out, word_count);
  size_t nsite, site, up_only = 0, down_only = 0, a, b, up_site, down_site;
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (model->domain != MVMC_KRYLOV_DOMAIN_GC_BALANCED || current == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  {
    const uintptr_t x = (uintptr_t)current, y = (uintptr_t)out;
    const uintptr_t bytes = (uintptr_t)(word_count * sizeof(*out));
    if (x < y + bytes && y < x + bytes)
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_krylov_fock_validate(model, current, word_count);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  memcpy(out, current, word_count * sizeof(*out));
  nsite = model->site_count;
  for (site = 0; site < nsite; ++site) {
    const int up = get_bit(current, site), down = get_bit(current, site + nsite);
    up_only += (size_t)(up && !down);
    down_only += (size_t)(down && !up);
  }
  if (up_only == 0 || down_only == 0) return MVMC_KRYLOV_STATUS_OK;
  status = bounded(draw, context, up_only, &a);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = bounded(draw, context, down_only, &b);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  up_site = nth_site(current, nsite, 1, 0, a);
  down_site = nth_site(current, nsite, 0, 1, b);
  if (up_site >= nsite || down_site >= nsite) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  flip_bit(out, up_site);
  flip_bit(out, up_site + nsite);
  flip_bit(out, down_site);
  flip_bit(out, down_site + nsite);
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_gc_krylov_draw_mixture(
    const MVMCKrylovFockModel *model, const MVMCGCSpinSupport *support,
    const uint64_t *current, MVMCKrylovFockProposalBoundedDraw draw,
    void *context, uint64_t *out, size_t word_count, MVMCGCMoveKind *kind) {
  MVMCKrylovStatus status = validate_target(model, draw, out, word_count);
  MVMCGCMoveKind move = MVMC_GC_MOVE_LOCAL;
  int finite = 0, s0_triplet = 0;
  size_t value;
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (support != NULL) {
    if (mvmc_gc_spin_support_validate(
            support, model->site_count > (size_t)(INT_MAX / 2)
                         ? 0
                         : (int)model->site_count,
            model->domain) != MVMC_KRYLOV_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    finite = support->active && support->total_spin > 0;
    s0_triplet = support->active && support->total_spin == 0 &&
                 support->pair_symmetry == MVMC_GC_PAIR_TRIPLET;
  }
  status = bounded(draw, context, 16, &value);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (value == 0) {
    move = MVMC_GC_MOVE_GLOBAL;
  } else if ((finite || s0_triplet) && value <= 8) {
    move = MVMC_GC_MOVE_TWO_PAIR;
  } else if (finite && value <= 10) {
    move = MVMC_GC_MOVE_SWAP;
  }
  if (kind != NULL) *kind = move;
  switch (move) {
    case MVMC_GC_MOVE_GLOBAL:
      return mvmc_gc_krylov_draw_uniform(model, draw, context, out,
                                         word_count);
    case MVMC_GC_MOVE_TWO_PAIR:
      return mvmc_gc_krylov_draw_two_pairs(model, current, draw, context, out,
                                           word_count);
    case MVMC_GC_MOVE_SWAP:
      return mvmc_gc_krylov_draw_spin_swap(model, current, draw, context, out,
                                           word_count);
    case MVMC_GC_MOVE_LOCAL:
      break;
  }
  return mvmc_gc_krylov_draw_local(model, current, draw, context, out,
                                   word_count);
}

/* Validated finite-S support of a balanced model; SINGLET/ZERO vanish. */
static MVMCKrylovStatus finite_spin_support(
    const MVMCKrylovFockModel *model, const MVMCGCSpinSupport *support) {
  if (model == NULL || support == NULL ||
      model->domain != MVMC_KRYLOV_DOMAIN_GC_BALANCED ||
      model->up_electron_count != 0 || model->down_electron_count != 0 ||
      model->pure_spin != 0 || model->site_count == 0 ||
      model->site_count > (size_t)(INT_MAX / 2) ||
      mvmc_gc_spin_support_validate(support, (int)model->site_count,
                                    model->domain) != MVMC_KRYLOV_STATUS_OK ||
      !support->active || support->total_spin <= 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return mvmc_gc_spin_support_all_zero(support)
             ? MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE
             : MVMC_KRYLOV_STATUS_OK;
}

/* Whether an even n in [0, 2L] is an allowed finite-S start sector. */
static int spin_sector_allowed(int site_count,
                               const MVMCGCSpinSupport *support, int n) {
  return GCSpinSectorPossible(
      site_count, support->total_spin,
      support->pair_symmetry == MVMC_GC_PAIR_TRIPLET, n);
}

MVMCKrylovStatus mvmc_gc_krylov_initial_sectors(
    const MVMCKrylovFockModel *model, const MVMCGCSpinSupport *support,
    int requested_n, int explicit_n,
    int *sectors, size_t capacity, size_t *count) {
  MVMCKrylovStatus status = finite_spin_support(model, support);
  int L;
  int64_t distance, first_distance, last_distance, lower, upper;
  size_t written = 0;
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  L = (int)model->site_count;
  if (sectors == NULL || count == NULL || capacity < (size_t)L + 1 ||
      (explicit_n != 0 && explicit_n != 1) || requested_n < 0 ||
      requested_n % 2 != 0 || requested_n > 2 * L) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (explicit_n) {
    if (!spin_sector_allowed(L, support, requested_n)) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    sectors[0] = requested_n;
    *count = 1;
    return MVMC_KRYLOV_STATUS_OK;
  }
  /* Restrict the search to the support interval. Use wide arithmetic for
   * both candidate sums and the loop's final increment at L=INT_MAX/2. */
  lower = 2 * (int64_t)support->total_spin;
  upper = 2 * ((int64_t)L - support->total_spin);
  first_distance = requested_n < lower ? lower - requested_n
                    : requested_n > upper ? requested_n - upper : 0;
  last_distance = llabs((int64_t)requested_n - lower);
  if (llabs((int64_t)requested_n - upper) > last_distance)
    last_distance = llabs((int64_t)requested_n - upper);
  /* Nearest first; the smaller N first at equal distance. */
  for (distance = first_distance; distance <= last_distance; distance += 2) {
    const int64_t below = (int64_t)requested_n - distance;
    const int64_t above = (int64_t)requested_n + distance;
    if (below >= lower && below <= upper &&
        spin_sector_allowed(L, support, (int)below)) {
      sectors[written++] = (int)below;
    }
    if (distance != 0 && above >= lower && above <= upper &&
        spin_sector_allowed(L, support, (int)above)) {
      sectors[written++] = (int)above;
    }
  }
  if (written == 0 || written > (size_t)L + 1) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  *count = written;
  return MVMC_KRYLOV_STATUS_OK;
}

/* k uniform sites among the sites unused by out, ascending, by sequential
 * selection; no draw when k is 0 or equals the remaining count. */
static MVMCKrylovStatus select_unused_sites(
    MVMCKrylovFockProposalBoundedDraw draw, void *context, size_t nsite,
    size_t k, int up, int down, uint64_t *out) {
  size_t site, remaining = 0;
  for (site = 0; site < nsite; ++site) {
    remaining += !get_bit(out, site) && !get_bit(out, site + nsite);
  }
  for (site = 0; site < nsite && k != 0; ++site) {
    int selected;
    if (get_bit(out, site) || get_bit(out, site + nsite)) continue;
    if (k == remaining) {
      selected = 1;
    } else {
      size_t value;
      const MVMCKrylovStatus status =
          bounded(draw, context, remaining, &value);
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      selected = value < k;
    }
    if (selected) {
      if (up) set_bit(out, site);
      if (down) set_bit(out, site + nsite);
      --k;
    }
    --remaining;
  }
  return k == 0 ? MVMC_KRYLOV_STATUS_OK
                : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

MVMCKrylovStatus mvmc_gc_krylov_draw_spin_initial(
    const MVMCKrylovFockModel *model, const MVMCGCSpinSupport *support,
    int initial_n, MVMCKrylovFockProposalBoundedDraw draw, void *context,
    uint64_t *out, size_t word_count) {
  MVMCKrylovStatus status = validate_target(model, draw, out, word_count);
  int L, m, low, high, doublons;
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = finite_spin_support(model, support);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  L = (int)model->site_count;
  if (initial_n < 0 || initial_n % 2 != 0 || initial_n > 2 * L ||
      !spin_sector_allowed(L, support, initial_n) ||
      GCSpinDoublonRange(L, support->total_spin, initial_n, &low, &high) !=
          1) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  m = initial_n / 2;
  doublons = low;
  if (high > low) {
    size_t value;
    status = bounded(draw, context, (size_t)(high - low) + 1, &value);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    doublons = low + (int)value;
  }
  memset(out, 0, word_count * sizeof(*out));
  status = select_unused_sites(draw, context, (size_t)L, (size_t)doublons, 1,
                               1, out);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = select_unused_sites(draw, context, (size_t)L,
                                 (size_t)(m - doublons), 1, 0, out);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = select_unused_sites(draw, context, (size_t)L,
                                 (size_t)(m - doublons), 0, 1, out);
  }
  return status;
}

#endif /* reference or bounded engine */
