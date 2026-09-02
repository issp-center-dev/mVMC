#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "include/blas_externs.h"
#include "include/gc_size.h"
#include "include/global.h"
#include "include/matrix_gc.h"
#include "include/workspace.h"

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
    default:
      return "unknown";
  }
}

static int CalculateMAllGCChild(
    const int ncur, const int *eleIdx, const int qpStart, const int qpidx,
    double complex *rawPf, double complex *rawInv, int *iwork,
    double complex *work, const int lwork, double *rwork) {
  const int globalQpidx = qpStart + qpidx;
  const size_t slaterStride =
      (size_t)Nsite2 * (size_t)Nsite2;
  const double complex *slater =
      SlaterElm + (size_t)globalQpidx * slaterStride;
  char uplo = 'U';
  char method = 'P';
  int n = ncur;
  int lda = ncur;
  int info = 0;
  int row;
  double complex pfaffian;

  for (row = 0; row < ncur; row++) {
    const int rsRow = eleIdx[row];
    int column;
    for (column = 0; column < ncur; column++) {
      const int rsColumn = eleIdx[column];
      rawPf[(size_t)column * (size_t)ncur + (size_t)row] =
          -slater[(size_t)rsRow * (size_t)Nsite2 + (size_t)rsColumn];
    }
  }
  memcpy(rawInv, rawPf,
         (size_t)ncur * (size_t)ncur * sizeof(*rawInv));

  M_ZSKPFA(&uplo, &method, &n, rawPf, &lda, &pfaffian, iwork, work,
           &lwork, rwork, &info);
  if (info != 0) return GC_MALL_PFAFFIAN;
  if (!(isfinite(creal(pfaffian)) && isfinite(cimag(pfaffian)))) {
    return GC_MALL_NONFINITE;
  }
  if (cabs(pfaffian) == 0.0) return GC_MALL_PFAFFIAN;

  M_ZGETRF(&n, &n, rawInv, &lda, iwork, &info);
  if (info != 0) return GC_MALL_GETRF;
  M_ZGETRI(&n, rawInv, &lda, iwork, work, &lwork, &info);
  if (info != 0) return GC_MALL_GETRI;

  for (row = 0; row < ncur; row++) {
    int column;
    for (column = 0; column < ncur; column++) {
      const size_t index =
          (size_t)row * (size_t)ncur + (size_t)column;
      rawInv[index] = -rawInv[index];
      if (!(isfinite(creal(rawInv[index])) &&
            isfinite(cimag(rawInv[index])))) {
        return GC_MALL_NONFINITE;
      }
    }
  }

  for (row = 0; row < ncur; row++) {
    memcpy(InvM + (size_t)qpidx * (size_t)NsizeMax * (size_t)NsizeMax +
               (size_t)row * (size_t)NsizeMax,
           rawInv + (size_t)row * (size_t)ncur,
           (size_t)ncur * sizeof(*rawInv));
  }
  /* GC amplitude convention: <x|phi_GC> = Pf(SlaterElm_x) so that
   * |phi_GC> = exp[sum_{IJ} f_IJ c^+_I c^+_J]|0>.  The Pfaffian above is
   * Pf(X) with X = -SlaterElm_x, and Pf(-A) = (-1)^{n/2} Pf(A), so the
   * sector-dependent sign is restored here.  InvM keeps the X convention. */
  PfM[qpidx] = ((ncur / 2) % 2 == 0) ? pfaffian : -pfaffian;
  return GC_MALL_OK;
}

int CalculateMAllGC_fcmp(const int ncur, const int *eleIdx,
                         const int qpStart, const int qpEnd) {
  size_t matrixCount;
  size_t complexCount;
  int matrixWorkspaceCount;
  int complexWorkspaceCount;
  int failureStatus = GC_MALL_OK;
  int failureQp = INT_MAX;
  int qpidx;
  int i;

  if (ncur < 0 || ncur > NsizeMax || (ncur % 2) != 0 ||
      NsizeMax <= 0 || Nsite2 <= 0 || NsizeMax > Nsite2 ||
      qpStart < 0 || qpEnd < qpStart || qpEnd > NQPFull ||
      SlaterElm == NULL || InvM == NULL || PfM == NULL ||
      LapackLWork < NsizeMax || (ncur > 0 && eleIdx == NULL)) {
    return GC_MALL_INVALID_ARGUMENT;
  }
  for (i = 0; i < ncur; i++) {
    if (eleIdx[i] < 0 || eleIdx[i] >= Nsite2) {
      return GC_MALL_INVALID_ARGUMENT;
    }
  }
  if (qpStart == qpEnd) return GC_MALL_OK;
  if (ncur == 0) {
    for (qpidx = 0; qpidx < qpEnd - qpStart; qpidx++) PfM[qpidx] = 1.0;
    return GC_MALL_OK;
  }

  matrixCount = GCCheckedMulSize((size_t)NsizeMax, (size_t)NsizeMax);
  complexCount = GCCheckedMulSize(2, matrixCount);
  complexCount = GCCheckedAddSize(complexCount, (size_t)LapackLWork);
  matrixWorkspaceCount =
      GCCheckedSizeToInt(matrixCount, "GC matrix workspace");
  complexWorkspaceCount =
      GCCheckedSizeToInt(complexCount, "GC complex workspace");
  RequestWorkSpaceThreadInt(NsizeMax);
  RequestWorkSpaceThreadComplex(complexWorkspaceCount);
  RequestWorkSpaceThreadDouble(LapackLWork);

#pragma omp parallel default(shared)
  {
    int *myIWork = GetWorkSpaceThreadInt(NsizeMax);
    double complex *myRawPf =
        GetWorkSpaceThreadComplex(matrixWorkspaceCount);
    double complex *myRawInv =
        GetWorkSpaceThreadComplex(matrixWorkspaceCount);
    double complex *myWork = GetWorkSpaceThreadComplex(LapackLWork);
    double *myRWork = GetWorkSpaceThreadDouble(LapackLWork);

#pragma omp for
    for (qpidx = 0; qpidx < qpEnd - qpStart; qpidx++) {
      const int myStatus = CalculateMAllGCChild(
          ncur, eleIdx, qpStart, qpidx, myRawPf, myRawInv, myIWork,
          myWork, LapackLWork, myRWork);
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

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  ReleaseWorkSpaceThreadDouble();
  if (failureStatus != GC_MALL_OK) {
    fprintf(stderr, "Error: GC matrix rebuild failed at qp=%d stage=%s.\n",
            failureQp, GCMatrixStageName(failureStatus));
  }
  return failureStatus;
}

#endif
