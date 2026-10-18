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

/* pi spin rotation multiplies an m-triplet-pair state by (-1)^m.
 * At S=0, odd m is forbidden. G=1 means no spin projection in mVMC.
 * Keep this fix restricted to S=0; it does not enable finite-S input. */
static inline int GCSpinForbiddenSector(int grand_canonical, int general,
                                        int gauss, int spin, int triplet,
                                        int electrons) {
  return grand_canonical && !general && gauss > 1 && spin == 0 &&
         triplet == 1 && electrons >= 0 && electrons % 4 == 2;
}

#endif
