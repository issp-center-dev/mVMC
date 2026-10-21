#ifndef MVMC_GC_SPIN_SUPPORT_H
#define MVMC_GC_SPIN_SUPPORT_H

#include <complex.h>
#include <math.h>
#include <stddef.h>

/* Test the actual, unrotated up-down pair function, including orbital signs.
 * Exact equality is essential: a nearby general pair function has nonzero
 * weight in both pair-number parities. Return -1 for invalid/nonfinite input.
 * Call again whenever the variational parameters change, before sampling. */
static inline int GCSpinPureTriplet(int sites, int parameters,
                                    int *const *index, int *const *sign,
                                    const double complex *value) {
  int i, j, triplet = 1;
  if (sites <= 0 || parameters <= 0 || index == NULL || sign == NULL ||
      value == NULL) return -1;
  for (i = 0; i < sites; ++i) {
    if (index[i] == NULL || sign[i] == NULL) return -1;
    for (j = 0; j < sites; ++j) {
      const int k = index[i][j];
      if (k < 0 || k >= parameters ||
          (sign[i][j] != 1 && sign[i][j] != -1) ||
          !isfinite(creal(value[k])) || !isfinite(cimag(value[k]))) return -1;
    }
  }
  for (i = 0; i < sites; ++i) {
    for (j = i; j < sites; ++j) {
      const double complex a = value[index[i][j]] * (double)sign[i][j];
      const double complex b = value[index[j][i]] * (double)sign[j][i];
      if (a != -b) triplet = 0;
    }
  }
  return triplet;
}

/* The same exact test for F^T=F. A pure singlet has no S>0 component; the
 * zero table satisfies both tests. Return -1 for invalid/nonfinite input. */
static inline int GCSpinPureSinglet(int sites, int parameters,
                                    int *const *index, int *const *sign,
                                    const double complex *value) {
  int i, j, singlet = 1;
  /* Shares the index/sign/finiteness validation of the triplet test. */
  if (GCSpinPureTriplet(sites, parameters, index, sign, value) < 0)
    return -1;
  for (i = 0; i < sites; ++i) {
    for (j = i + 1; j < sites; ++j) {
      const double complex a = value[index[i][j]] * (double)sign[i][j];
      const double complex b = value[index[j][i]] * (double)sign[j][i];
      if (a != b) singlet = 0;
    }
  }
  return singlet;
}

/* pi spin rotation multiplies an m-triplet-pair state by (-1)^m and a
 * total-spin-S, M=0 state by (-1)^S, so m-S must be even. At S=0 this is
 * the N%4==2 rule. G=1 means no spin projection in mVMC. */
static inline int GCSpinForbiddenSector(int grand_canonical, int general,
                                        int gauss, int spin, int triplet,
                                        int electrons) {
  return grand_canonical && !general && gauss > 1 && spin >= 0 &&
         triplet == 1 && electrons >= 0 && electrons % 2 == 0 &&
         (electrons / 2 - spin) % 2 != 0;
}

/* Number of singly occupied sites, or -1 for a non-0/1 occupation. */
static inline int GCSpinCountSingles(const int *eleNum, int nsite) {
  int site, singles = 0;
  if (eleNum == NULL || nsite <= 0) return -1;
  for (site = 0; site < nsite; ++site) {
    const int up = eleNum[site], down = eleNum[site + nsite];
    if ((up != 0 && up != 1) || (down != 0 && down != 1)) return -1;
    singles += up != down;
  }
  return singles;
}

/* Doubly occupied and empty sites are spin singlets, so a configuration
 * with n singles carries total spin at most n/2 (parameter independent). */
static inline int GCSpinConfigurationZero(int spin, int singles) {
  return spin > 0 && singles < 2 * spin;
}

/* Doublon counts D of an N=2m, M=0 configuration with at least 2S singles:
 * max(0, 2m-L) <= D <= m-S. Return 1 with the range, 0 when empty, -1 for
 * invalid input. A nonempty range is exactly 2S <= N <= 2L-2S. */
static inline int GCSpinDoublonRange(int nsite, int spin, int electrons,
                                     int *low, int *high) {
  int m;
  if (nsite <= 0 || spin < 0 || electrons < 0 || electrons % 2 != 0 ||
      electrons > 2 * nsite || low == NULL || high == NULL)
    return -1;
  m = electrons / 2;
  *low = 2 * m - nsite > 0 ? 2 * m - nsite : 0;
  *high = m - spin;
  return *low <= *high;
}

/* Whether an even N sector can carry total spin S at all, including the
 * pair-parity rule of an exactly antisymmetric F. A general F may still
 * cancel there; that is found by evaluation, not by this test. */
static inline int GCSpinSectorPossible(int nsite, int spin, int triplet,
                                       int electrons) {
  int low, high;
  if (GCSpinDoublonRange(nsite, spin, electrons, &low, &high) != 1)
    return 0;
  return !(triplet == 1 && (electrons / 2 - spin) % 2 != 0);
}

#endif
