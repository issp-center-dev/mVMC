/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_GC_KRYLOV_SPIN_SUPPORT_H
#define MVMC_GC_KRYLOV_SPIN_SUPPORT_H

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                         \
    defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "gc_spin_support.h"
#include "krylov_fock_reference.h"

#include <limits.h>
#include <stddef.h>

/*
 * Immutable spin-projection support of one GC Lanczos generation.  The
 * pure integer rules of gc_spin_support.h are reused; nothing here reads
 * production globals, allocates, draws random numbers or communicates.
 *
 * inactive: the only representation is {0, 1, 0, GENERAL} (no mask).
 * active:   GC_BALANCED, G >= 2, 0 <= S <= L/2.  S = 0 uses GENERAL or
 *           TRIPLET only (ZERO folds to TRIPLET and SINGLET to GENERAL when
 *           the descriptor is built); S > 0 uses all four classes of the
 *           exact real-F tests.  validate never folds a value silently.
 */
typedef enum {
  MVMC_GC_PAIR_GENERAL = 0,
  MVMC_GC_PAIR_SINGLET = 1,
  MVMC_GC_PAIR_TRIPLET = 2,
  MVMC_GC_PAIR_ZERO = 3
} MVMCGCPairSymmetry;

typedef struct {
  int active;                 /* 0 or 1 */
  int gauss_count;            /* G */
  int total_spin;             /* S, not 2S */
  MVMCGCPairSymmetry pair_symmetry;
} MVMCGCSpinSupport;

typedef struct {
  int electron_count;
  int single_count;
  unsigned zero_reasons;      /* CONFIGURATION=1, PARITY=2, SINGLET=4 */
} MVMCGCSpinSupportResult;

#define MVMC_GC_SPIN_SUPPORT_SCHEMA 1
#define MVMC_GC_SPIN_ZERO_CONFIGURATION 1u
#define MVMC_GC_SPIN_ZERO_PARITY 2u
#define MVMC_GC_SPIN_ZERO_SINGLET 4u

static inline MVMCGCSpinSupport mvmc_gc_spin_support_inactive(void) {
  MVMCGCSpinSupport support;
  support.active = 0;
  support.gauss_count = 1;
  support.total_spin = 0;
  support.pair_symmetry = MVMC_GC_PAIR_GENERAL;
  return support;
}

/*
 * Canonical descriptor from the exact classification of the unrotated
 * pair function of this generation (GCSpinPureTriplet/Singlet, 1 or 0).
 * Inactive projection gives the canonical inactive value; an inactive G
 * other than 1 is kept so that validate rejects it.  At S = 0 singlet is
 * not read.  On INVALID_ARGUMENT *support is unchanged.
 */
static inline MVMCKrylovStatus mvmc_gc_spin_support_classify(
    int active, int gauss_count, int total_spin, int triplet, int singlet,
    MVMCGCSpinSupport *support) {
  MVMCGCSpinSupport made = mvmc_gc_spin_support_inactive();
  if (support == NULL || (active != 0 && active != 1)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (!active) {
    made.gauss_count = gauss_count;
    *support = made;
    return MVMC_KRYLOV_STATUS_OK;
  }
  if ((triplet != 0 && triplet != 1) ||
      (total_spin > 0 && singlet != 0 && singlet != 1)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  made.active = 1;
  made.gauss_count = gauss_count;
  made.total_spin = total_spin;
  if (total_spin == 0) {
    made.pair_symmetry = triplet ? MVMC_GC_PAIR_TRIPLET : MVMC_GC_PAIR_GENERAL;
  } else if (triplet && singlet) {
    made.pair_symmetry = MVMC_GC_PAIR_ZERO;
  } else if (triplet) {
    made.pair_symmetry = MVMC_GC_PAIR_TRIPLET;
  } else if (singlet) {
    made.pair_symmetry = MVMC_GC_PAIR_SINGLET;
  } else {
    made.pair_symmetry = MVMC_GC_PAIR_GENERAL;
  }
  *support = made;
  return MVMC_KRYLOV_STATUS_OK;
}

static inline MVMCKrylovStatus mvmc_gc_spin_support_validate(
    const MVMCGCSpinSupport *support, int site_count,
    MVMCKrylovConfigurationDomain domain) {
  if (support == NULL || site_count < 1 || site_count > INT_MAX / 2 ||
      (domain != MVMC_KRYLOV_DOMAIN_GC_BALANCED &&
       domain != MVMC_KRYLOV_DOMAIN_GC_EVEN)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (support->active == 0) {
    return support->gauss_count == 1 && support->total_spin == 0 &&
                   support->pair_symmetry == MVMC_GC_PAIR_GENERAL
               ? MVMC_KRYLOV_STATUS_OK
               : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (support->active != 1 || domain != MVMC_KRYLOV_DOMAIN_GC_BALANCED ||
      support->gauss_count < 2 || support->total_spin < 0 ||
      support->total_spin > site_count / 2) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  switch ((int)support->pair_symmetry) {
    case MVMC_GC_PAIR_GENERAL:
    case MVMC_GC_PAIR_TRIPLET:
      return MVMC_KRYLOV_STATUS_OK;
    case MVMC_GC_PAIR_SINGLET:
    case MVMC_GC_PAIR_ZERO:
      return support->total_spin > 0 ? MVMC_KRYLOV_STATUS_OK
                                     : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    default:
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
}

/* Whether this descriptor annihilates every configuration (S > 0 singlet). */
static inline int mvmc_gc_spin_support_all_zero(
    const MVMCGCSpinSupport *support) {
  return support != NULL && support->active == 1 && support->total_spin > 0 &&
         (support->pair_symmetry == MVMC_GC_PAIR_SINGLET ||
          support->pair_symmetry == MVMC_GC_PAIR_ZERO);
}

/*
 * Structural zero of one configuration: occupation is the eleNum layout
 * (up 0..L-1, down L..2L-1) with occupation_count = 2L entries of 0/1,
 * balanced in GC_BALANCED and of even N in GC_EVEN.  The configuration is
 * checked even when inactive (zero_reasons is then 0).  On failure *result
 * is unchanged.
 */
static inline MVMCKrylovStatus mvmc_gc_spin_support_evaluate(
    const MVMCGCSpinSupport *support, int site_count,
    MVMCKrylovConfigurationDomain domain,
    const int *occupation, size_t occupation_count,
    MVMCGCSpinSupportResult *result) {
  MVMCKrylovStatus status =
      mvmc_gc_spin_support_validate(support, site_count, domain);
  int site, up = 0, down = 0, singles, electrons;
  unsigned reasons = 0;
  int triplet;
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (occupation == NULL || result == NULL ||
      occupation_count != 2 * (size_t)site_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (site = 0; site < 2 * site_count; ++site) {
    if (occupation[site] != 0 && occupation[site] != 1) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    if (site < site_count) {
      up += occupation[site];
    } else {
      down += occupation[site];
    }
  }
  electrons = up + down;
  if ((domain == MVMC_KRYLOV_DOMAIN_GC_BALANCED && up != down) ||
      electrons % 2 != 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  singles = GCSpinCountSingles(occupation, site_count);
  if (singles < 0) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (support->active) {
    triplet = support->pair_symmetry == MVMC_GC_PAIR_TRIPLET ||
              support->pair_symmetry == MVMC_GC_PAIR_ZERO;
    if (GCSpinConfigurationZero(support->total_spin, singles)) {
      reasons |= MVMC_GC_SPIN_ZERO_CONFIGURATION;
    }
    if (GCSpinForbiddenSector(1, 0, support->gauss_count,
                              support->total_spin, triplet, electrons)) {
      reasons |= MVMC_GC_SPIN_ZERO_PARITY;
    }
    if (mvmc_gc_spin_support_all_zero(support)) {
      reasons |= MVMC_GC_SPIN_ZERO_SINGLET;
    }
  }
  result->electron_count = electrons;
  result->single_count = singles;
  result->zero_reasons = reasons;
  return MVMC_KRYLOV_STATUS_OK;
}

#endif /* reference or bounded engine */

#endif /* MVMC_GC_KRYLOV_SPIN_SUPPORT_H */
