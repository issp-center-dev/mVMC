#include <stddef.h>

#include "nbody_operator.h"

static void SetReduction(NBodyReduction *reduction, NBodyReducedKind kind,
                         int order, int sign) {
  if (reduction == NULL) return;
  reduction->kind = kind;
  reduction->order = order;
  reduction->sign = sign;
}

static int ReductionAliasesArray(const NBodyReduction *reduction,
                                 const int *array) {
  return (const void *)reduction == (const void *)array;
}

static void RemoveFactor(int index, int order, int *rsiWork, int *rsjWork) {
  int k;
  for (k = index; k < order - 1; k++) {
    rsiWork[k] = rsiWork[k + 1];
    rsjWork[k] = rsjWork[k + 1];
  }
}

int ReduceNBodyTerm(int n, const int *rsi, const int *rsj,
                    const int *occupation, int nOrbitals,
                    int *rsiWork, int *rsjWork, int workCapacity,
                    NBodyReduction *reduction) {
  int k;
  int order;
  int sign = 1;

  if (reduction == NULL) return -1;
  if (ReductionAliasesArray(reduction, rsi)
      || ReductionAliasesArray(reduction, rsj)
      || ReductionAliasesArray(reduction, occupation)
      || ReductionAliasesArray(reduction, rsiWork)
      || ReductionAliasesArray(reduction, rsjWork)) {
    return -1;
  }
  SetReduction(reduction, NBODY_REDUCED_INVALID, 0, 1);

  if (n <= 0 || nOrbitals <= 0 || workCapacity < n ||
      rsi == NULL || rsj == NULL || occupation == NULL ||
      rsiWork == NULL || rsjWork == NULL ||
      rsiWork == rsjWork ||
      rsiWork == rsi || rsiWork == rsj || rsiWork == occupation ||
      rsjWork == rsi || rsjWork == rsj || rsjWork == occupation) {
    return -1;
  }

  for (k = 0; k < nOrbitals; k++) {
    if (occupation[k] != 0 && occupation[k] != 1) return -1;
  }
  for (k = 0; k < n; k++) {
    if (rsi[k] < 0 || rsi[k] >= nOrbitals ||
        rsj[k] < 0 || rsj[k] >= nOrbitals) {
      return -1;
    }
  }

  for (k = 0; k < n; k++) {
    rsiWork[k] = rsi[k];
    rsjWork[k] = rsj[k];
  }
  order = n;

  while (order > 0) {
    int changed = 0;

    for (k = order - 1; k >= 0; k--) {
      int l;
      int orbital = rsjWork[k];

      for (l = k + 1; l < order; l++) {
        if (orbital == rsiWork[l]) {
          rsjWork[k] = rsjWork[l];
          RemoveFactor(l, order, rsiWork, rsjWork);
          order--;
          changed = 1;
          break;
        }
        if (orbital == rsjWork[l]) {
          SetReduction(reduction, NBODY_REDUCED_ZERO, 0, sign);
          return 0;
        }
      }
      if (changed) break;
      if (occupation[orbital] == 0) {
        SetReduction(reduction, NBODY_REDUCED_ZERO, 0, sign);
        return 0;
      }

      orbital = rsiWork[k];
      if (orbital == rsjWork[k]) {
        RemoveFactor(k, order, rsiWork, rsjWork);
        order--;
        changed = 1;
        break;
      }

      for (l = k + 1; l < order; l++) {
        if (orbital == rsiWork[l]) {
          SetReduction(reduction, NBODY_REDUCED_ZERO, 0, sign);
          return 0;
        }
        if (orbital == rsjWork[l]) {
          rsiWork[k] = rsiWork[l];
          RemoveFactor(l, order, rsiWork, rsjWork);
          order--;
          sign = -sign;
          changed = 1;
          break;
        }
      }
      if (changed) break;
      if (occupation[orbital] == 1) {
        SetReduction(reduction, NBODY_REDUCED_ZERO, 0, sign);
        return 0;
      }
    }

    if (!changed) break;
  }

  if (order == 0) {
    SetReduction(reduction, NBODY_REDUCED_SCALAR, 0, sign);
  } else {
    SetReduction(reduction, NBODY_REDUCED_HOPS, order, sign);
  }
  return 0;
}
