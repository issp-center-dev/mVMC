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
#ifndef _PFUDATE_FSZ_SRC
#define _PFUDATE_FSZ_SRC

static int BF_FSZ_MultiplySize(const size_t a, const size_t b, size_t *result) {
  if(result == NULL || (a != 0 && b > SIZE_MAX/a)) return 1;
  *result = a*b;
  return 0;
}

static int BF_FSZ_PfUpdateWorkSizeForK(const int k, size_t *complexCount,
                                       size_t *intCount, size_t *doubleCount) {
  size_t kk, n, vecCount, smallCount, total;
  if(k < 1 || k > Nsize || complexCount == NULL || intCount == NULL
      || doubleCount == NULL) return 1;
  kk = (size_t)k;
  n = (size_t)Nsize;
  if(BF_FSZ_MultiplySize(kk, n, &vecCount) != 0
      || BF_FSZ_MultiplySize(kk, kk, &smallCount) != 0
      || smallCount > SIZE_MAX/4) return 1;
  smallCount *= 4; /* (2*k)^2 */
  if(vecCount > SIZE_MAX-n || vecCount+n > SIZE_MAX-smallCount
      || vecCount+n+smallCount > SIZE_MAX-smallCount) return 1;
  total = vecCount + n + smallCount + smallCount;
  if(kk > SIZE_MAX/2) return 1;
  *complexCount = total;
  *intCount = 2*kk;
  *doubleCount = smallCount;
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

int BF_FSZ_ShouldUseFullPfaffian(const int nAffected) {
  return nAffected >= BFFSZPfUpdateKFull || nAffected == Nsize;
}

static int CalculateNewPfMBF_fsz_rows_child(
    const int nAffected, const int *affected,
    double complex *pfMNew, const double complex *oldPfM,
    const double complex *oldInvM, const size_t invMQpStride,
    const int *eleIdx, const int *eleSpn, const int qpStart,
    const int qpidx, const double complex *candidateSlater,
    int *failureDetail, double complex *vec, double complex *w,
    double complex *smallMatrix, int *iwork, double complex *pfWork,
    double *rwork) {
  const int nsize = Nsize;
  const int n2 = 2*nAffected;
  const int globalQpidx = qpStart + qpidx;
  const double complex *sltE = candidateSlater
      + (size_t)globalQpidx*(size_t)Nsite2*(size_t)Nsite2;
  const double complex *invM = oldInvM + (size_t)qpidx*invMQpStride;
  char uplo = 'U', method = 'P';
  int k,l,i,j,nn=n2,lda=n2,info=0,lwork=n2*n2;
  double complex pfaff, result;
  double sign;

  for(k=0;k<nAffected;k++) {
    const int rsk = eleIdx[affected[k]] + eleSpn[affected[k]]*Nsite;
    const double complex *sltE_k = sltE + (size_t)rsk*(size_t)Nsite2;
    double complex *vec_k = vec + (size_t)k*(size_t)nsize;
    for(i=0;i<nsize;i++) {
      const int rsi = eleIdx[i] + eleSpn[i]*Nsite;
      vec_k[i] = sltE_k[rsi];
    }
  }

  /* w = InvM * vec_k. Keeping one w row at a time makes this scalar
     implementation O(k*Nsize^2 + k^2*Nsize) without a k*Nsize cache. */
  for(k=0;k<nAffected;k++) {
    const double complex *vec_k = vec + (size_t)k*(size_t)nsize;
    double complex *matY_k = smallMatrix + (size_t)n2*(size_t)k + nAffected;
    for(i=0;i<nsize;i++) {
      double complex value = 0.0 + 0.0*I;
      const double complex *invM_i = invM + (size_t)i*(size_t)nsize;
      for(j=0;j<nsize;j++) value += invM_i[j]*vec_k[j];
      w[i] = value;
    }

    for(l=0;l<k;l++) {
      const double complex *vec_l = vec + (size_t)l*(size_t)nsize;
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
  if(qpNum > 1 && invMQpStride > SIZE_MAX/(qpNum-1)) {
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
    const int status = CalculateNewPfMBF_fsz_rows_child(
        nAffected, affected, pfMNew, oldPfM, oldInvM, invMQpStride,
        eleIdx, eleSpn, qpStart, qpidx, candidateSlater, failureDetail,
        vec, w, smallMatrix, iwork, pfWork, rwork);
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
