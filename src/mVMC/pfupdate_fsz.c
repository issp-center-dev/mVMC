/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

his program is developed based on the mVMC-mini program
(https://github.com/fiber-miniapp/mVMC-mini)
which follows "The BSD 3-Clause License".

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details. 

You should have received a copy of the GNU General Public License 
along with this program. If not, see http://www.gnu.org/licenses/. 
*/
/*-------------------------------------------------------------
 * Variational Monte Carlo
 * fast Pfaffian update
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

#include "pfupdate_fsz.h"
#include <stdint.h>
#include <time.h>
#ifndef _PFUDATE_FSZ_SRC
#define _PFUDATE_FSZ_SRC

static int BF_FSZ_MultiplySize(const size_t a, const size_t b, size_t *result) {
  if(result == NULL || (a != 0 && b > SIZE_MAX/a)) return 1;
  *result = a*b;
  return 0;
}

static int BF_FSZ_MemoryRangesOverlap(
    const void *rangeA, const size_t countA, const size_t elementSizeA,
    const void *rangeB, const size_t countB, const size_t elementSizeB) {
  size_t bytesA, bytesB;
  uintptr_t beginA, beginB;
  if(countA == 0 || countB == 0) return 0;
  if(rangeA == NULL || rangeB == NULL || elementSizeA == 0
      || elementSizeB == 0 || countA > SIZE_MAX/elementSizeA
      || countB > SIZE_MAX/elementSizeB) return 1;
  bytesA = countA*elementSizeA;
  bytesB = countB*elementSizeB;
  beginA = (uintptr_t)rangeA;
  beginB = (uintptr_t)rangeB;
  if(beginA > UINTPTR_MAX-bytesA || beginB > UINTPTR_MAX-bytesB) return 1;
  return beginA < beginB+bytesB && beginB < beginA+bytesA;
}

static int BF_FSZ_ComplexRangesOverlap(
    const double complex *rangeA, const size_t countA,
    const double complex *rangeB, const size_t countB) {
  return BF_FSZ_MemoryRangesOverlap(
      rangeA,countA,sizeof(double complex),
      rangeB,countB,sizeof(double complex));
}

static int BF_FSZ_PfUpdateTailWorkSizeForK(
    const int k, size_t *complexCount,
    size_t *intCount, size_t *doubleCount) {
  size_t kk, n, smallCount, total;
  if(k < 1 || k > Nsize || complexCount == NULL || intCount == NULL
      || doubleCount == NULL) return 1;
  kk = (size_t)k;
  n = (size_t)Nsize;
  if(BF_FSZ_MultiplySize(kk, kk, &smallCount) != 0
      || smallCount > SIZE_MAX/4) return 1;
  smallCount *= 4; /* (2*k)^2 */
  if(n > SIZE_MAX-smallCount || n+smallCount > SIZE_MAX-smallCount) return 1;
  total = n+smallCount+smallCount;
  if(kk > SIZE_MAX/2) return 1;
  *complexCount = total;
  *intCount = 2*kk;
  *doubleCount = smallCount;
  return 0;
}

static int BF_FSZ_PfUpdateWorkSizeForK(const int k, size_t *complexCount,
                                       size_t *intCount, size_t *doubleCount) {
  size_t vecCount, tailCount;
  if(BF_FSZ_PfUpdateTailWorkSizeForK(
      k, &tailCount, intCount, doubleCount) != 0
      || BF_FSZ_MultiplySize((size_t)k, (size_t)Nsize, &vecCount) != 0
      || vecCount > SIZE_MAX-tailCount) return 1;
  *complexCount = vecCount+tailCount;
  return 0;
}

int GetCalculateNewPfMBF_fsz_rows_work_size(size_t *complexCount,
    size_t *intCount, size_t *doubleCount) {
  int kMax;
  if(Nsize < 2 || complexCount == NULL || intCount == NULL
      || doubleCount == NULL) return 1;
  kMax = Nsize - 1; /* nAffected==Nsize always selects the full path. */
  if(kMax >= BF_FSZ_PF_UPDATE_KFULL_DEFAULT) {
    kMax = BF_FSZ_PF_UPDATE_KFULL_DEFAULT - 1;
  }
  return BF_FSZ_PfUpdateWorkSizeForK(kMax, complexCount, intCount, doubleCount);
}

int GetCalculateNewPfMBF_fsz_row_values_work_size(size_t *complexCount,
    size_t *intCount, size_t *doubleCount) {
  int kMax;
  if(Nsize < 2 || complexCount == NULL || intCount == NULL
      || doubleCount == NULL) return 1;
  kMax = Nsize-1;
  if(kMax >= BF_FSZ_PF_UPDATE_KFULL_DEFAULT) {
    kMax = BF_FSZ_PF_UPDATE_KFULL_DEFAULT-1;
  }
  return BF_FSZ_PfUpdateTailWorkSizeForK(
      kMax, complexCount, intCount, doubleCount);
}

int BF_FSZ_ShouldUseFullPfaffian(const int nAffected) {
  return nAffected >= BFFSZPfUpdateKFull || nAffected == Nsize;
}

static int CalculateNewPfMBF_fsz_row_values_child(
    const int nAffected, const int *affected,
    double complex *pfMNew, const double complex *oldPfM,
    const double complex *oldInvM, const size_t invMQpStride,
    const int qpStart, const int qpidx,
    const double complex *candidateRows, const size_t rowAffectedStride,
    int *failureDetail, double complex *w,
    double complex *smallMatrix, int *iwork, double complex *pfWork,
    double *rwork) {
  const int nsize = Nsize;
  const int n2 = 2*nAffected;
  const int globalQpidx = qpStart + qpidx;
  const double complex *invM = oldInvM + (size_t)qpidx*invMQpStride;
  char uplo = 'U', method = 'P';
  int k,l,i,j,nn=n2,lda=n2,info=0,lwork=n2*n2;
  double complex pfaff, result;
  double sign;

  /* w = InvM * vec_k. Keeping one w row at a time makes this scalar
     implementation O(k*Nsize^2 + k^2*Nsize) without a k*Nsize cache. */
  for(k=0;k<nAffected;k++) {
    const double complex *vec_k
        = candidateRows+(size_t)k*rowAffectedStride;
    double complex *matY_k = smallMatrix + (size_t)n2*(size_t)k + nAffected;
    for(i=0;i<nsize;i++) {
      double complex value = 0.0 + 0.0*I;
      const double complex *invM_i = invM + (size_t)i*(size_t)nsize;
      for(j=0;j<nsize;j++) value += invM_i[j]*vec_k[j];
      w[i] = value;
    }

    for(l=0;l<k;l++) {
      const double complex *vec_l
          = candidateRows+(size_t)l*rowAffectedStride;
      double complex value = 0.0 + 0.0*I;
      for(i=0;i<nsize;i++) value += w[i]*vec_l[i];
      smallMatrix[(size_t)n2*(size_t)l+k] = value + vec_l[affected[k]];
    }
    for(l=0;l<nAffected;l++) matY_k[l] = w[affected[l]];
  }

  for(k=0;k<nAffected;k++) {
    const double complex *invM_k = invM + (size_t)affected[k]*(size_t)nsize;
    double complex *matZ_k = smallMatrix
        + (size_t)n2*(size_t)(k+nAffected) + nAffected;
    for(l=k+1;l<nAffected;l++) matZ_k[l] = invM_k[affected[l]];
  }

  for(k=0;k<n2;k++) {
    for(l=0;l<k;l++) {
      smallMatrix[(size_t)n2*(size_t)k+l]
          = -smallMatrix[(size_t)n2*(size_t)l+k];
    }
    smallMatrix[(size_t)n2*(size_t)k+k] = 0.0 + 0.0*I;
  }

  M_ZSKPFA(&uplo, &method, &nn, smallMatrix, &lda, &pfaff,
      iwork, pfWork, &lwork, rwork, &info);
  if(info != 0) {
    if(failureDetail != NULL) *failureDetail = info;
    return BF_FSZ_PF_UPDATE_LAPACK_FAILURE;
  }
  if(!(isfinite(creal(pfaff)) && isfinite(cimag(pfaff)))) {
    if(failureDetail != NULL) *failureDetail = globalQpidx;
    return BF_FSZ_PF_UPDATE_NONFINITE;
  }
  if(pfaff == 0.0 + 0.0*I || oldPfM[qpidx] == 0.0 + 0.0*I) {
    if(failureDetail != NULL) *failureDetail = globalQpidx;
    return BF_FSZ_PF_UPDATE_EXACT_ZERO;
  }
  if(!(isfinite(creal(oldPfM[qpidx])) && isfinite(cimag(oldPfM[qpidx])))) {
    if(failureDetail != NULL) *failureDetail = globalQpidx;
    return BF_FSZ_PF_UPDATE_NONFINITE;
  }

  sign = (((nAffected*(nAffected-1)/2) & 1) == 0) ? 1.0 : -1.0;
  result = sign*pfaff*oldPfM[qpidx];
  if(!(isfinite(creal(result)) && isfinite(cimag(result)))) {
    if(failureDetail != NULL) *failureDetail = globalQpidx;
    return BF_FSZ_PF_UPDATE_NONFINITE;
  }
  pfMNew[qpidx] = result;
  return BF_FSZ_PF_UPDATE_OK;
}

int CalculateNewPfMBF_fsz_rows_workspace(
    const int nAffected, const int *affected,
    double complex *pfMNew, const double complex *oldPfM,
    const double complex *oldInvM, const size_t invMQpStride,
    const int *eleIdx, const int *eleSpn,
    const int qpStart, const int qpEnd,
    const double complex *candidateSlater, int *failureDetail,
    double complex *complexWork, const size_t complexWorkCount,
    int *iwork, const size_t intWorkCount,
    double *rwork, const size_t rworkCount) {
  size_t requiredComplex, requiredInt, requiredDouble;
  size_t vecCount, smallCount;
  size_t qpNum;
  double complex *vec, *w, *smallMatrix, *pfWork;
  int i,j,qpidx;

  if(failureDetail != NULL) *failureDetail = 0;
  if(nAffected < 1 || nAffected > Nsize
      || qpStart < 0 || qpEnd < qpStart || qpEnd > NQPFull
      || BF_FSZ_ShouldUseFullPfaffian(nAffected)) {
    return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
  }
  if(affected == NULL || pfMNew == NULL || oldPfM == NULL || oldInvM == NULL
      || eleIdx == NULL || eleSpn == NULL || candidateSlater == NULL
      || complexWork == NULL || iwork == NULL || rwork == NULL
      || invMQpStride < (size_t)Nsize*(size_t)Nsize
      || invMQpStride > (size_t)PTRDIFF_MAX
      || pfMNew == oldPfM || candidateSlater == SlaterElmBF) {
    return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
  }
  qpNum = (size_t)(qpEnd-qpStart);
  if(qpNum > 1
      && (qpNum-1) > (SIZE_MAX-(size_t)Nsize*(size_t)Nsize)/invMQpStride) {
    return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
  }
  if(BF_FSZ_PfUpdateWorkSizeForK(nAffected, &requiredComplex,
      &requiredInt, &requiredDouble) != 0
      || complexWorkCount < requiredComplex || intWorkCount < requiredInt
      || rworkCount < requiredDouble) {
    return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
  }
  for(i=0;i<nAffected;i++) {
    if(affected[i] < 0 || affected[i] >= Nsize) {
      return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
    }
    for(j=0;j<i;j++) if(affected[j] == affected[i]) {
      return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
    }
  }
  for(i=0;i<Nsize;i++) {
    if(eleIdx[i] < 0 || eleIdx[i] >= Nsite
        || (eleSpn[i] != 0 && eleSpn[i] != 1)) {
      return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
    }
  }

  vecCount = (size_t)nAffected*(size_t)Nsize;
  smallCount = 4*(size_t)nAffected*(size_t)nAffected;
  vec = complexWork;
  w = vec + vecCount;
  smallMatrix = w + Nsize;
  pfWork = smallMatrix + smallCount;

  for(qpidx=0;qpidx<qpEnd-qpStart;qpidx++) {
    const int globalQpidx = qpStart+qpidx;
    const double complex *sltE = candidateSlater
        +(size_t)globalQpidx*(size_t)Nsite2*(size_t)Nsite2;
    int k;
    for(k=0;k<nAffected;k++) {
      const int rsk = eleIdx[affected[k]]+eleSpn[affected[k]]*Nsite;
      const double complex *sltE_k = sltE+(size_t)rsk*(size_t)Nsite2;
      double complex *vec_k = vec+(size_t)k*(size_t)Nsize;
      for(i=0;i<Nsize;i++) {
        const int rsi = eleIdx[i]+eleSpn[i]*Nsite;
        vec_k[i] = sltE_k[rsi];
      }
    }
    {
      const int status = CalculateNewPfMBF_fsz_row_values_child(
        nAffected, affected, pfMNew, oldPfM, oldInvM, invMQpStride,
        qpStart, qpidx, vec, (size_t)Nsize, failureDetail,
        w, smallMatrix, iwork, pfWork, rwork);
      if(status != BF_FSZ_PF_UPDATE_OK) return status;
    }
  }
  return BF_FSZ_PF_UPDATE_OK;
}

int CalculateNewPfMBF_fsz_row_values_workspace(
    const int nAffected, const int *affected,
    double complex *pfMNew, const double complex *oldPfM,
    const double complex *oldInvM, const size_t invMQpStride,
    const int qpStart, const int qpEnd,
    const double complex *candidateRows, const size_t rowQpStride,
    const size_t rowAffectedStride, int *failureDetail,
    double complex *complexWork, const size_t complexWorkCount,
    int *iwork, const size_t intWorkCount,
    double *rwork, const size_t rworkCount) {
  const size_t qpNum = (size_t)(qpEnd-qpStart);
  size_t requiredComplex, requiredInt, requiredDouble;
  size_t smallCount, rowDataCount = 0, rowSpan = 0;
  double complex *w, *smallMatrix, *pfWork;
  int i, j, qpidx;

  if(failureDetail != NULL) *failureDetail = 0;
  if(nAffected < 1 || nAffected > Nsize
      || qpStart < 0 || qpEnd < qpStart || qpEnd > NQPFull
      || BF_FSZ_ShouldUseFullPfaffian(nAffected)
      || affected == NULL || pfMNew == NULL || oldPfM == NULL
      || oldInvM == NULL || complexWork == NULL || iwork == NULL
      || rwork == NULL || invMQpStride < (size_t)Nsize*(size_t)Nsize
      || invMQpStride > (size_t)PTRDIFF_MAX
      || pfMNew == oldPfM
      || (qpNum > 1
          && (qpNum-1)
              > (SIZE_MAX-(size_t)Nsize*(size_t)Nsize)/invMQpStride)) {
    return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
  }
  if(BF_FSZ_PfUpdateTailWorkSizeForK(nAffected, &requiredComplex,
      &requiredInt, &requiredDouble) != 0
      || complexWorkCount < requiredComplex || intWorkCount < requiredInt
      || rworkCount < requiredDouble) {
    return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
  }
  for(i=0;i<nAffected;i++) {
    if(affected[i] < 0 || affected[i] >= Nsize) {
      return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
    }
    for(j=0;j<i;j++) if(affected[j] == affected[i]) {
      return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
    }
  }
  if(qpNum > 0) {
    if(candidateRows == NULL || rowAffectedStride < (size_t)Nsize
        || (size_t)(nAffected-1)
            > (SIZE_MAX-(size_t)Nsize)/rowAffectedStride) {
      return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
    }
    rowDataCount = (size_t)(nAffected-1)*rowAffectedStride+(size_t)Nsize;
    if(rowQpStride < rowDataCount
        || (qpNum-1) > (SIZE_MAX-rowDataCount)/rowQpStride) {
      return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
    }
    rowSpan = (qpNum-1)*rowQpStride+rowDataCount;
    if(rowSpan > (size_t)PTRDIFF_MAX
        || BF_FSZ_ComplexRangesOverlap(
            candidateRows, rowSpan, complexWork, complexWorkCount)
        || BF_FSZ_ComplexRangesOverlap(
            candidateRows, rowSpan, pfMNew, qpNum)
        || BF_FSZ_MemoryRangesOverlap(
            candidateRows,rowSpan,sizeof(double complex),
            iwork,intWorkCount,sizeof(int))
        || BF_FSZ_MemoryRangesOverlap(
            candidateRows,rowSpan,sizeof(double complex),
            rwork,rworkCount,sizeof(double))) {
      return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
    }
  }

  smallCount = 4*(size_t)nAffected*(size_t)nAffected;
  w = complexWork;
  smallMatrix = w+Nsize;
  pfWork = smallMatrix+smallCount;
  for(qpidx=0;qpidx<qpEnd-qpStart;qpidx++) {
    const int status = CalculateNewPfMBF_fsz_row_values_child(
        nAffected, affected, pfMNew, oldPfM, oldInvM, invMQpStride,
        qpStart, qpidx, candidateRows+(size_t)qpidx*rowQpStride,
        rowAffectedStride, failureDetail, w, smallMatrix, iwork, pfWork, rwork);
    if(status != BF_FSZ_PF_UPDATE_OK) return status;
  }
  return BF_FSZ_PF_UPDATE_OK;
}

int CalculateNewPfMBF_fsz_rows(
    const int nAffected, const int *affected,
    double complex *pfMNew, const double complex *oldPfM,
    const double complex *oldInvM, const size_t invMQpStride,
    const int *eleIdx, const int *eleSpn,
    const int qpStart, const int qpEnd,
    const double complex *candidateSlater, int *failureDetail) {
  size_t complexCount, intCount, doubleCount;
  double complex *complexWork = NULL;
  int *iwork = NULL;
  double *rwork = NULL;
  int status;

  if(BF_FSZ_PfUpdateWorkSizeForK(nAffected, &complexCount, &intCount,
      &doubleCount) != 0 || complexCount > SIZE_MAX/sizeof(double complex)
      || intCount > SIZE_MAX/sizeof(int)
      || doubleCount > SIZE_MAX/sizeof(double)) {
    return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
  }
  complexWork = (double complex *)malloc(sizeof(double complex)*complexCount);
  iwork = (int *)malloc(sizeof(int)*intCount);
  rwork = (double *)malloc(sizeof(double)*doubleCount);
  if(complexWork == NULL || iwork == NULL || rwork == NULL) {
    free(complexWork);
    free(iwork);
    free(rwork);
    return BF_FSZ_PF_UPDATE_INVALID_ARGUMENT;
  }
  status = CalculateNewPfMBF_fsz_rows_workspace(
      nAffected, affected, pfMNew, oldPfM, oldInvM, invMQpStride,
      eleIdx, eleSpn, qpStart, qpEnd, candidateSlater, failureDetail,
      complexWork, complexCount, iwork, intCount, rwork, doubleCount);
  free(complexWork);
  free(iwork);
  free(rwork);
  return status;
}

static BF_FSZ_InvUpdateResult BF_FSZ_InvResultValue(const int status,
    const int stage, const int qpidx, const int lapackInfo,
    const double antisymmetryResidual, const double affectedResidual,
    const double checkSeconds) {
  BF_FSZ_InvUpdateResult result;
  result.status = status;
  result.stage = stage;
  result.qpidx = qpidx;
  result.lapackInfo = lapackInfo;
  result.antisymmetryResidual = antisymmetryResidual;
  result.affectedResidual = affectedResidual;
  result.checkSeconds = checkSeconds;
  return result;
}

static double BF_FSZ_TimeDifference(const struct timespec *start,
                                    const struct timespec *end) {
  return (double)(end->tv_sec-start->tv_sec)
      + 1.0e-9*(double)(end->tv_nsec-start->tv_nsec);
}

static int BF_FSZ_InvUpdateTailWorkSizeForK(
    const int k, size_t *complexCount, size_t *intCount) {
  size_t kk,n,kn,small,n2,total;
  if(k < 1 || k > Nsize || complexCount == NULL || intCount == NULL) return 1;
  kk = (size_t)k;
  n = (size_t)Nsize;
  if(BF_FSZ_MultiplySize(kk,n,&kn) != 0
      || BF_FSZ_MultiplySize(kk,kk,&small) != 0
      || small > SIZE_MAX/4 || kk > SIZE_MAX/2) return 1;
  small *= 4;
  n2 = 2*kk;
  /* w + matUV + correctionTmp + small matrix + inverse + GETRI work */
  if(kn > SIZE_MAX/5 || small > SIZE_MAX/2
      || 5*kn > SIZE_MAX-2*small
      || 5*kn+2*small > SIZE_MAX-n2) return 1;
  total = 5*kn+2*small+n2;
  *complexCount = total;
  *intCount = n2;
  return 0;
}

static int BF_FSZ_InvUpdateWorkSizeForK(const int k, size_t *complexCount,
                                        size_t *intCount) {
  size_t vecCount, tailCount;
  if(BF_FSZ_InvUpdateTailWorkSizeForK(k, &tailCount, intCount) != 0
      || BF_FSZ_MultiplySize((size_t)k, (size_t)Nsize, &vecCount) != 0
      || vecCount > SIZE_MAX-tailCount) return 1;
  *complexCount = vecCount+tailCount;
  return 0;
}

int GetPrepareInvMBF_fsz_rows_work_size(size_t *complexCount,
    size_t *intCount) {
  int kMax;
  if(Nsize < 2 || complexCount == NULL || intCount == NULL) return 1;
  kMax = Nsize-1;
  if(kMax >= BF_FSZ_PF_UPDATE_KFULL_DEFAULT) {
    kMax = BF_FSZ_PF_UPDATE_KFULL_DEFAULT-1;
  }
  return BF_FSZ_InvUpdateWorkSizeForK(kMax,complexCount,intCount);
}

int GetPrepareInvMBF_fsz_row_values_work_size(size_t *complexCount,
    size_t *intCount) {
  int kMax;
  if(Nsize < 2 || complexCount == NULL || intCount == NULL) return 1;
  kMax = Nsize-1;
  if(kMax >= BF_FSZ_PF_UPDATE_KFULL_DEFAULT) {
    kMax = BF_FSZ_PF_UPDATE_KFULL_DEFAULT-1;
  }
  return BF_FSZ_InvUpdateTailWorkSizeForK(kMax,complexCount,intCount);
}

static BF_FSZ_InvUpdateResult PrepareInvMBF_fsz_row_values_child(
    const int nAffected, const int *affected,
    const double complex *oldInvM, const size_t oldInvMQpStride,
    double complex *invMNew, const size_t newInvMQpStride,
    const int qpStart, const int qpidx,
    const double complex *candidateRows, const size_t rowAffectedStride,
    double complex *w, double complex *smallMatrix,
    double complex *smallInverse, double complex *matUV,
    double complex *correctionTmp, double complex *lapackWork, int *iwork) {
  const int nsize = Nsize;
  const int n2 = 2*nAffected;
  const int globalQpidx = qpStart+qpidx;
  const double complex *oldInv = oldInvM+(size_t)qpidx*oldInvMQpStride;
  double complex *newInv = invMNew+(size_t)qpidx*newInvMQpStride;
  int i,j,k,l,m=n2,n=n2,lda=n2,info=0,lwork=n2;
  double maxAbs=0.0,maxSkew=0.0,maxResidual=0.0;
  struct timespec checkStart,checkEnd;

  for(k=0;k<nAffected;k++) {
    const double complex *vec_k
        = candidateRows+(size_t)k*rowAffectedStride;
    double complex *w_k = w+(size_t)k*(size_t)nsize;
    for(i=0;i<nsize;i++) {
      double complex value = 0.0+0.0*I;
      const double complex *oldInv_i = oldInv+(size_t)i*(size_t)nsize;
      for(j=0;j<nsize;j++) value += oldInv_i[j]*vec_k[j];
      w_k[i] = value;
    }
  }

  for(k=0;k<nAffected;k++) {
    const double complex *vec_k
        = candidateRows+(size_t)k*rowAffectedStride;
    const double complex *w_k = w+(size_t)k*(size_t)nsize;
    double complex *matY_k = smallMatrix+(size_t)n2*(size_t)k+nAffected;
    for(l=k+1;l<nAffected;l++) {
      const double complex *w_l = w+(size_t)l*(size_t)nsize;
      double complex value = 0.0+0.0*I;
      for(i=0;i<nsize;i++) value += w_l[i]*vec_k[i];
      smallMatrix[(size_t)n2*(size_t)k+l]
          = value+vec_k[affected[l]];
    }
    for(l=0;l<nAffected;l++) matY_k[l] = w_k[affected[l]];
  }
  for(k=0;k<nAffected;k++) {
    const double complex *oldInv_k
        = oldInv+(size_t)affected[k]*(size_t)nsize;
    double complex *matZ_k = smallMatrix
        +(size_t)n2*(size_t)(k+nAffected)+nAffected;
    for(l=k+1;l<nAffected;l++) matZ_k[l] = oldInv_k[affected[l]];
  }
  for(k=0;k<n2;k++) {
    for(l=0;l<k;l++) {
      smallMatrix[(size_t)n2*(size_t)k+l]
          = -smallMatrix[(size_t)n2*(size_t)l+k];
    }
    smallMatrix[(size_t)n2*(size_t)k+k] = 0.0+0.0*I;
  }
  for(k=0;k<n2;k++) for(l=0;l<n2;l++) {
    smallInverse[(size_t)n2*(size_t)k+l]
        = smallMatrix[(size_t)n2*(size_t)l+k];
  }

  for(i=0;i<nsize;i++) {
    double complex *matUV_i = matUV+(size_t)i*(size_t)n2;
    const double complex *oldInv_i = oldInv+(size_t)i*(size_t)nsize;
    for(k=0;k<nAffected;k++) {
      matUV_i[k] = -w[(size_t)k*(size_t)nsize+i]
          - ((affected[k] == i) ? 1.0 : 0.0);
      matUV_i[k+nAffected] = oldInv_i[affected[k]];
    }
  }

  M_ZGETRF(&m,&n,smallInverse,&lda,iwork,&info);
  if(info != 0) {
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_LAPACK_FAILURE,
        BF_FSZ_INV_STAGE_GETRF,globalQpidx,info,0.0,0.0,0.0);
  }
  M_ZGETRI(&n,smallInverse,&lda,iwork,lapackWork,&lwork,&info);
  if(info != 0) {
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_LAPACK_FAILURE,
        BF_FSZ_INV_STAGE_GETRI,globalQpidx,info,0.0,0.0,0.0);
  }

  for(j=0;j<nsize;j++) {
    const double complex *matUV_j = matUV+(size_t)j*(size_t)n2;
    double complex *tmp_j = correctionTmp+(size_t)j*(size_t)n2;
    for(k=0;k<n2;k++) {
      double complex value = 0.0+0.0*I;
      const double complex *smallInverse_k
          = smallInverse+(size_t)k*(size_t)n2;
      for(l=0;l<n2;l++) value += smallInverse_k[l]*matUV_j[l];
      tmp_j[k] = value;
    }
  }
  for(i=0;i<nsize;i++) {
    const double complex *matUV_i = matUV+(size_t)i*(size_t)n2;
    double complex *newInv_i = newInv+(size_t)i*(size_t)nsize;
    const double complex *oldInv_i = oldInv+(size_t)i*(size_t)nsize;
    for(j=0;j<nsize;j++) {
      const double complex *tmp_j = correctionTmp+(size_t)j*(size_t)n2;
      double complex value = 0.0+0.0*I;
      for(k=0;k<n2;k++) value += matUV_i[k]*tmp_j[k];
      newInv_i[j] = oldInv_i[j]-value;
    }
  }

  clock_gettime(CLOCK_MONOTONIC,&checkStart);
  for(i=0;i<nsize;i++) for(j=0;j<nsize;j++) {
    const double valueAbs = cabs(newInv[(size_t)i*(size_t)nsize+j]);
    const double skewAbs = cabs(newInv[(size_t)i*(size_t)nsize+j]
        +newInv[(size_t)j*(size_t)nsize+i]);
    if(!isfinite(valueAbs)) {
      clock_gettime(CLOCK_MONOTONIC,&checkEnd);
      return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_NONFINITE,
          BF_FSZ_INV_STAGE_FINITE,globalQpidx,0,maxSkew,maxResidual,
          BF_FSZ_TimeDifference(&checkStart,&checkEnd));
    }
    if(valueAbs > maxAbs) maxAbs = valueAbs;
    if(skewAbs > maxSkew) maxSkew = skewAbs;
  }
  if(maxSkew/fmax(1.0,maxAbs) > 1.0e-9) {
    clock_gettime(CLOCK_MONOTONIC,&checkEnd);
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_ANTISYMMETRY,
        BF_FSZ_INV_STAGE_ANTISYMMETRY,globalQpidx,0,
        maxSkew/fmax(1.0,maxAbs),0.0,
        BF_FSZ_TimeDifference(&checkStart,&checkEnd));
  }
  for(i=0;i<nsize;i++) {
    newInv[(size_t)i*(size_t)nsize+i] = 0.0+0.0*I;
    for(j=i+1;j<nsize;j++) {
      const double complex value = 0.5
          *(newInv[(size_t)i*(size_t)nsize+j]
            -newInv[(size_t)j*(size_t)nsize+i]);
      newInv[(size_t)i*(size_t)nsize+j] = value;
      newInv[(size_t)j*(size_t)nsize+i] = -value;
    }
  }
  for(k=0;k<nAffected;k++) {
    const int row = affected[k];
    const double complex *candidateRow
        = candidateRows+(size_t)k*rowAffectedStride;
    for(j=0;j<nsize;j++) {
      double complex value = 0.0+0.0*I;
      for(i=0;i<nsize;i++) {
        value += candidateRow[i]*newInv[(size_t)i*(size_t)nsize+j];
      }
      value -= (row == j) ? 1.0 : 0.0;
      {
        const double residualAbs = cabs(value);
        if(!isfinite(residualAbs)) {
          maxResidual = INFINITY;
        } else if(residualAbs > maxResidual) {
          maxResidual = residualAbs;
        }
      }
    }
  }
  clock_gettime(CLOCK_MONOTONIC,&checkEnd);
  if(!(isfinite(maxResidual) && maxResidual <= 1.0e-9)) {
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_RESIDUAL,
        BF_FSZ_INV_STAGE_RESIDUAL,globalQpidx,0,
        maxSkew/fmax(1.0,maxAbs),maxResidual,
        BF_FSZ_TimeDifference(&checkStart,&checkEnd));
  }
  return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_OK,BF_FSZ_INV_STAGE_NONE,
      globalQpidx,0,maxSkew/fmax(1.0,maxAbs),maxResidual,
      BF_FSZ_TimeDifference(&checkStart,&checkEnd));
}

BF_FSZ_InvUpdateResult PrepareInvMBF_fsz_rows_workspace(
    const int nAffected, const int *affected,
    const double complex *oldInvM, const size_t oldInvMQpStride,
    double complex *invMNew, const size_t newInvMQpStride,
    const int *eleIdx, const int *eleSpn,
    const int qpStart, const int qpEnd,
    const double complex *candidateSlater,
    double complex *complexWork, const size_t complexWorkCount,
    int *iwork, const size_t intWorkCount) {
  const size_t nsizeSquared = (size_t)Nsize*(size_t)Nsize;
  const size_t qpNum = (size_t)(qpEnd-qpStart);
  size_t requiredComplex,requiredInt,kn,small,n2;
  double complex *vec,*w,*smallMatrix,*smallInverse,*matUV;
  double complex *correctionTmp,*lapackWork;
  double totalCheckSeconds=0.0,maxAntisymmetry=0.0,maxResidual=0.0;
  int i,j,qpidx;

  if(nAffected < 1 || nAffected > Nsize
      || qpStart < 0 || qpEnd < qpStart || qpEnd > NQPFull
      || BF_FSZ_ShouldUseFullPfaffian(nAffected)) {
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_INVALID_ARGUMENT,
        BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
  }
  if(qpNum == 0) {
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_OK,
        BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
  }
  if(affected == NULL || oldInvM == NULL || invMNew == NULL
      || eleIdx == NULL || eleSpn == NULL || candidateSlater == NULL
      || complexWork == NULL || iwork == NULL || invMNew == oldInvM
      || candidateSlater == SlaterElmBF
      || oldInvMQpStride < nsizeSquared || newInvMQpStride < nsizeSquared
      || oldInvMQpStride > (size_t)PTRDIFF_MAX
      || newInvMQpStride > (size_t)PTRDIFF_MAX
      || (qpNum > 1
          && ((qpNum-1) > (SIZE_MAX-nsizeSquared)/oldInvMQpStride
            || (qpNum-1) > (SIZE_MAX-nsizeSquared)/newInvMQpStride))) {
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_INVALID_ARGUMENT,
        BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
  }
  if(BF_FSZ_InvUpdateWorkSizeForK(nAffected,&requiredComplex,&requiredInt) != 0
      || complexWorkCount < requiredComplex || intWorkCount < requiredInt) {
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_INVALID_ARGUMENT,
        BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
  }
  for(i=0;i<nAffected;i++) {
    if(affected[i] < 0 || affected[i] >= Nsize) {
      return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_INVALID_ARGUMENT,
          BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
    }
    for(j=0;j<i;j++) if(affected[i] == affected[j]) {
      return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_INVALID_ARGUMENT,
          BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
    }
  }
  for(i=0;i<Nsize;i++) {
    if(eleIdx[i] < 0 || eleIdx[i] >= Nsite
        || (eleSpn[i] != 0 && eleSpn[i] != 1)) {
      return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_INVALID_ARGUMENT,
          BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
    }
  }

  kn = (size_t)nAffected*(size_t)Nsize;
  small = 4*(size_t)nAffected*(size_t)nAffected;
  n2 = 2*(size_t)nAffected;
  vec = complexWork;
  w = vec+kn;
  smallMatrix = w+kn;
  smallInverse = smallMatrix+small;
  matUV = smallInverse+small;
  correctionTmp = matUV+2*kn;
  lapackWork = correctionTmp+2*kn;

  for(qpidx=0;qpidx<qpEnd-qpStart;qpidx++) {
    const int globalQpidx = qpStart+qpidx;
    const double complex *sltE = candidateSlater
        +(size_t)globalQpidx*(size_t)Nsite2*(size_t)Nsite2;
    int k;
    for(k=0;k<nAffected;k++) {
      const int rsk = eleIdx[affected[k]]+eleSpn[affected[k]]*Nsite;
      const double complex *sltE_k = sltE+(size_t)rsk*(size_t)Nsite2;
      double complex *vec_k = vec+(size_t)k*(size_t)Nsize;
      for(i=0;i<Nsize;i++) {
        const int rsi = eleIdx[i]+eleSpn[i]*Nsite;
        vec_k[i] = sltE_k[rsi];
      }
    }
    BF_FSZ_InvUpdateResult result = PrepareInvMBF_fsz_row_values_child(
        nAffected,affected,oldInvM,oldInvMQpStride,invMNew,newInvMQpStride,
        qpStart,qpidx,vec,(size_t)Nsize,w,smallMatrix,
        smallInverse,matUV,correctionTmp,lapackWork,iwork);
    totalCheckSeconds += result.checkSeconds;
    if(result.antisymmetryResidual > maxAntisymmetry) {
      maxAntisymmetry = result.antisymmetryResidual;
    }
    if(result.affectedResidual > maxResidual) {
      maxResidual = result.affectedResidual;
    }
    if(result.status != BF_FSZ_INV_UPDATE_OK) {
      result.checkSeconds = totalCheckSeconds;
      result.antisymmetryResidual = maxAntisymmetry;
      result.affectedResidual = maxResidual;
      return result;
    }
  }
  return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_OK,BF_FSZ_INV_STAGE_NONE,
      qpStart,0,maxAntisymmetry,maxResidual,totalCheckSeconds);
}

BF_FSZ_InvUpdateResult PrepareInvMBF_fsz_row_values_workspace(
    const int nAffected, const int *affected,
    const double complex *oldInvM, const size_t oldInvMQpStride,
    double complex *invMNew, const size_t newInvMQpStride,
    const int qpStart, const int qpEnd,
    const double complex *candidateRows, const size_t rowQpStride,
    const size_t rowAffectedStride,
    double complex *complexWork, const size_t complexWorkCount,
    int *iwork, const size_t intWorkCount) {
  const size_t nsizeSquared = (size_t)Nsize*(size_t)Nsize;
  const size_t qpNum = (size_t)(qpEnd-qpStart);
  size_t requiredComplex,requiredInt,kn,small,n2;
  size_t rowDataCount,rowSpan,newInvSpan;
  double complex *w,*smallMatrix,*smallInverse,*matUV;
  double complex *correctionTmp,*lapackWork;
  double totalCheckSeconds=0.0,maxAntisymmetry=0.0,maxResidual=0.0;
  int i,j,qpidx;

  if(nAffected < 1 || nAffected > Nsize
      || qpStart < 0 || qpEnd < qpStart || qpEnd > NQPFull
      || BF_FSZ_ShouldUseFullPfaffian(nAffected)) {
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_INVALID_ARGUMENT,
        BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
  }
  if(qpNum == 0) {
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_OK,
        BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
  }
  if(affected == NULL || oldInvM == NULL || invMNew == NULL
      || candidateRows == NULL || complexWork == NULL || iwork == NULL
      || invMNew == oldInvM
      || oldInvMQpStride < nsizeSquared || newInvMQpStride < nsizeSquared
      || oldInvMQpStride > (size_t)PTRDIFF_MAX
      || newInvMQpStride > (size_t)PTRDIFF_MAX
      || (qpNum > 1
          && ((qpNum-1) > (SIZE_MAX-nsizeSquared)/oldInvMQpStride
            || (qpNum-1) > (SIZE_MAX-nsizeSquared)/newInvMQpStride))) {
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_INVALID_ARGUMENT,
        BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
  }
  if(BF_FSZ_InvUpdateTailWorkSizeForK(
      nAffected,&requiredComplex,&requiredInt) != 0
      || complexWorkCount < requiredComplex || intWorkCount < requiredInt) {
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_INVALID_ARGUMENT,
        BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
  }
  for(i=0;i<nAffected;i++) {
    if(affected[i] < 0 || affected[i] >= Nsize) {
      return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_INVALID_ARGUMENT,
          BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
    }
    for(j=0;j<i;j++) if(affected[i] == affected[j]) {
      return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_INVALID_ARGUMENT,
          BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
    }
  }
  if(rowAffectedStride < (size_t)Nsize
      || (size_t)(nAffected-1)
          > (SIZE_MAX-(size_t)Nsize)/rowAffectedStride) {
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_INVALID_ARGUMENT,
        BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
  }
  rowDataCount = (size_t)(nAffected-1)*rowAffectedStride+(size_t)Nsize;
  if(rowQpStride < rowDataCount
      || (qpNum-1) > (SIZE_MAX-rowDataCount)/rowQpStride) {
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_INVALID_ARGUMENT,
        BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
  }
  rowSpan = (qpNum-1)*rowQpStride+rowDataCount;
  newInvSpan = (qpNum-1)*newInvMQpStride+nsizeSquared;
  if(rowSpan > (size_t)PTRDIFF_MAX || newInvSpan > (size_t)PTRDIFF_MAX
      || BF_FSZ_ComplexRangesOverlap(
          candidateRows,rowSpan,complexWork,complexWorkCount)
      || BF_FSZ_ComplexRangesOverlap(
          candidateRows,rowSpan,invMNew,newInvSpan)
      || BF_FSZ_MemoryRangesOverlap(
          candidateRows,rowSpan,sizeof(double complex),
          iwork,intWorkCount,sizeof(int))) {
    return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_INVALID_ARGUMENT,
        BF_FSZ_INV_STAGE_NONE,qpStart,0,0.0,0.0,0.0);
  }

  kn = (size_t)nAffected*(size_t)Nsize;
  small = 4*(size_t)nAffected*(size_t)nAffected;
  n2 = 2*(size_t)nAffected;
  w = complexWork;
  smallMatrix = w+kn;
  smallInverse = smallMatrix+small;
  matUV = smallInverse+small;
  correctionTmp = matUV+2*kn;
  lapackWork = correctionTmp+2*kn;

  for(qpidx=0;qpidx<qpEnd-qpStart;qpidx++) {
    BF_FSZ_InvUpdateResult result
        = PrepareInvMBF_fsz_row_values_child(
            nAffected,affected,oldInvM,oldInvMQpStride,
            invMNew,newInvMQpStride,qpStart,qpidx,
            candidateRows+(size_t)qpidx*rowQpStride,rowAffectedStride,
            w,smallMatrix,smallInverse,matUV,correctionTmp,lapackWork,iwork);
    totalCheckSeconds += result.checkSeconds;
    if(result.antisymmetryResidual > maxAntisymmetry) {
      maxAntisymmetry = result.antisymmetryResidual;
    }
    if(result.affectedResidual > maxResidual) {
      maxResidual = result.affectedResidual;
    }
    if(result.status != BF_FSZ_INV_UPDATE_OK) {
      result.checkSeconds = totalCheckSeconds;
      result.antisymmetryResidual = maxAntisymmetry;
      result.affectedResidual = maxResidual;
      return result;
    }
  }
  return BF_FSZ_InvResultValue(BF_FSZ_INV_UPDATE_OK,BF_FSZ_INV_STAGE_NONE,
      qpStart,0,maxAntisymmetry,maxResidual,totalCheckSeconds);
}


/* Calculate new pfaffian. The ma-th electron with spin s hops. */
void CalculateNewPfM_fsz(const int ma, const int s, double complex *pfMNew, const int *eleIdx,const int *eleSpn,
                     const int qpStart, const int qpEnd) {
  #pragma procedure serial
  const int qpNum = qpEnd-qpStart;
  //const int msa = ma+s*Ne;//fsz
  const int msa = ma;
  const int rsa = eleIdx[msa] + s*Nsite;

  int qpidx;
  int msj,rsj;
  const double complex *sltE_a; /* update elements of msa-th row */
  const double complex *invM_a;
  double complex ratio;

  /* optimization for Kei */
  const int nsize = Nsize;

  #pragma loop noalias
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    sltE_a = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2 + rsa*Nsite2;
    invM_a = InvM + qpidx*Nsize*Nsize + msa*Nsize;

    ratio = 0.0;
    #pragma omp simd private(rsj) reduction(+:ratio)
    for(msj=0;msj<nsize;msj++) { //fsz
      rsj = eleIdx[msj]+eleSpn[msj]*Nsite;//fsz
      ratio += invM_a[msj] * sltE_a[rsj];
    }
//    for(msj=ne;msj<nsize;msj++) {
//      rsj = eleIdx[msj] + Nsite;
//      ratio += invM_a[msj] * sltE_a[rsj];
 //   }

    pfMNew[qpidx] = -ratio*PfM[qpidx];
  }

  return;
}

/* thread parallel version of CalculateNewPfM_fsz */
void CalculateNewPfM2_fsz(const int ma, const int s, double complex *pfMNew, const int *eleIdx,const int *eleSpn,
                     const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  //const int msa = ma+s*Ne;
  const int msa = ma;
  const int rsa = eleIdx[msa] + s*Nsite;

  int qpidx;
  int msj,rsj;
  const double complex *sltE_a; /* update elements of msa-th row */
  const double complex *invM_a;
  double complex ratio;

  /* optimization for Kei */
  const int nsize = Nsize;

  #pragma omp parallel for default(shared)        \
    private(qpidx,msj,sltE_a,invM_a,ratio,rsj)
  #pragma loop noalias
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    sltE_a = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2 + rsa*Nsite2;
    invM_a = InvM + qpidx*Nsize*Nsize + msa*Nsize;

    ratio = 0.0;
    #pragma omp simd private(rsj) reduction(+:ratio)
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj]+eleSpn[msj]*Nsite;//fsz
      ratio += invM_a[msj] * sltE_a[rsj];
    }
/*    for(msj=ne;msj<nsize;msj++) {
      rsj = eleIdx[msj] + Nsite;
      ratio += invM_a[msj] * sltE_a[rsj];
    }*/

    pfMNew[qpidx] = -ratio*PfM[qpidx];
  }

  return;
}

/* Update PfM and InvM. The ma-th electron with spin s hops to site ra=eleIdx[msi] */
void UpdateMAll_fsz(const int ma, const int s, const int *eleIdx,const int *eleSpn,
                const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  int qpidx;
  double complex *vec1,*vec2;

  RequestWorkSpaceThreadComplex(2*Nsize);

  #pragma omp parallel default(shared) private(vec1,vec2)
  {
    vec1 = GetWorkSpaceThreadComplex(Nsize);
    vec2 = GetWorkSpaceThreadComplex(Nsize);
   
    #pragma omp for private(qpidx)
    #pragma loop nounroll
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      updateMAll_child_fsz(ma, s, eleIdx,eleSpn, qpStart, qpEnd, qpidx, vec1, vec2);
    }
  }

  ReleaseWorkSpaceThreadComplex();
  return;
}

void updateMAll_child_fsz(const int ma, const int s, const int *eleIdx,const int *eleSpn,
                      const int qpStart, const int qpEnd, const int qpidx,
                      double complex *vec1, double complex *vec2) {
  #pragma procedure serial
  /* const int qpNum = qpEnd-qpStart; */
  //const int msa = ma+s*Ne;
  const int msa = ma;
  const int rsa = eleIdx[msa] + s*Nsite;
  const int nsize = Nsize; /* optimization for Kei */

  int msi,msj,rsj;

  const double complex *sltE_a; /* update elements of msa-th row */
  double complex sltE_aj;
  double complex *invM;
  double complex *invM_i,*invM_j,*invM_a;

  double complex vec1_i,vec2_i;
  double complex invVec1_a;
  double complex tmp;

  sltE_a = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2 + rsa*Nsite2;

  invM = InvM + qpidx*Nsize*Nsize;
  invM_a = invM + msa*Nsize;

  #pragma omp simd
  for(msi=0;msi<nsize;msi++) vec1[msi] = 0.0+0.0*I; //TBC

  /* Calculate vec1[i] = sum_j invM[i][j] sltE[a][j] */
  /* Note that invM[i][j] = -invM[j][i] */
  #pragma loop noalias
  for(msj=0;msj<nsize;msj++) {
    //rsj = eleIdx[msj] + (msj/Ne)*Nsite;
    rsj = eleIdx[msj] + eleSpn[msj]*Nsite; //fsz
    sltE_aj = sltE_a[rsj];
    invM_j = invM + msj*Nsize;

    #pragma omp simd
    for(msi=0;msi<nsize;msi++) {
      vec1[msi] += -invM_j[msi] * sltE_aj;
    }
  }

  /* Update Pfaffian */
  /* Calculate -1.0/bufV_a to reduce devision */
  tmp = vec1[msa];
  PfM[qpidx] *= -tmp;
  invVec1_a = -1.0/tmp;

  /* Calculate vec2[i] = -InvM[a][i]/vec1[a] */
  #pragma loop noalias
  #pragma omp simd
  for(msi=0;msi<nsize;msi++) {
    vec2[msi] = invM_a[msi] * invVec1_a;
  }

  /* Update InvM */
  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    invM_i = invM + msi*Nsize;
    vec1_i = vec1[msi];
    vec2_i = vec2[msi];

    #pragma omp simd
    for(msj=0;msj<nsize;msj++) {
      invM_i[msj] += vec1_i * vec2[msj] - vec1[msj] * vec2_i;
    }

    invM_i[msa] -= vec2_i;
  }

  #pragma loop noalias
  #pragma omp simd
  for(msj=0;msj<nsize;msj++) {
    invM_a[msj] += vec2[msj];
  }
  /* end of update invM */

  return;
}

#endif
