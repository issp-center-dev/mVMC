#include <complex.h>
#include <limits.h>
#include <stdio.h>

#include "include/gc_component.h"
#include "include/global.h"
#include "include/matrix_gc.h"

#ifndef _SRC_MATRIX_GC
#define _SRC_MATRIX_GC

static const char *GCMatrixStageName(const int status) {
  switch (status) {
    case GC_MALL_PFAFFIAN:
      return "Pfaffian";
    case GC_MALL_GETRF:
      return "GETRF";
    case GC_MALL_GETRI:
      return "GETRI";
    case GC_MALL_NONFINITE:
      return "nonfinite";
    case GC_MALL_FP_MODE:
      return "strict-FP build required";
    case GC_MALL_ALLOCATION:
      return "allocation";
    default:
      return "unknown";
  }
}

int CalculateMAllGC_fcmp(const int ncur, const int *eleIdx,
                         const int qpStart, const int qpEnd) {
  int failureStatus = GC_MALL_OK;
  int failureQp = INT_MAX;
  int qpidx;
  int i;

  if (ncur < 0 || ncur > NsizeMax || (ncur % 2) != 0 ||
      NsizeMax <= 0 || Nsite2 <= 0 || NsizeMax > Nsite2 ||
      qpStart < 0 || qpEnd < qpStart || qpEnd > NQPFull ||
      SlaterElm == NULL || InvM == NULL || PfM == NULL ||
      GCQpState == NULL || (ncur > 0 && eleIdx == NULL)) {
    return GC_MALL_INVALID_ARGUMENT;
  }
  for (i = 0; i < ncur; i++) {
    if (eleIdx[i] < 0 || eleIdx[i] >= Nsite2) {
      return GC_MALL_INVALID_ARGUMENT;
    }
  }
  if (qpStart == qpEnd) return GC_MALL_OK;
  if (ncur == 0) {
    for (qpidx = 0; qpidx < qpEnd - qpStart; qpidx++) {
      PfM[qpidx] = 1.0;
      GCQpState[qpidx] = GC_QP_INVERSE;
    }
    return GC_MALL_OK;
  }

#pragma omp parallel default(shared)
  {
    GCComponentCache *cache = GCComponentThreadCache();

#pragma omp for
    for (qpidx = 0; qpidx < qpEnd - qpStart; qpidx++) {
      const int myStatus = GCComponentRebuild(
          qpStart + qpidx, qpidx, ncur, eleIdx, cache);
      if (myStatus != GC_MALL_OK) {
#pragma omp critical
        {
          if (qpStart + qpidx < failureQp) {
            failureStatus = myStatus;
            failureQp = qpStart + qpidx;
          }
        }
      }
    }
  }
  if (failureStatus != GC_MALL_OK) {
    fprintf(stderr, "Error: GC matrix rebuild failed at qp=%d stage=%s.\n",
            failureQp, GCMatrixStageName(failureStatus));
  }
  return failureStatus;
}

#endif
