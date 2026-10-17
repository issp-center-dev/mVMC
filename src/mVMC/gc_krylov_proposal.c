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

#endif /* reference or bounded engine */
