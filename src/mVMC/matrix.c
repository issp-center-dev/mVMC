/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is developed based on the mVMC-mini program
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
 * matrix Package (LAPACK and Pfapack)
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "./include/matrix.h"
#ifndef _SRC_MATRIX
#define _SRC_MATRIX
#include <complex.h>
#include <stdint.h>
#include "./include/global.h"
#include "workspace.c"

#define D_PfLimit 1.0e-100

/* Full skew factorization can report a zero separator even for an
 * invertible matrix. For UPLO=U, only T(0,1), T(2,3), ... are divisors
 * in utu2inv. Inspect all of them: INFO records only the first zero and
 * its parity cannot rule out a later zero paired pivot. */
static int fullSkewPairedPivotsReal(const double *a, int n, int lda) {
  if(n <= 0 || n % 2 != 0) return 0;
  for(int i=0; i<n; i+=2) {
    const double pivot = a[i+(size_t)(i+1)*lda];
    if(!isfinite(pivot) || pivot == 0.0) return 0;
  }
  return 1;
}

static int fullSkewPairedPivotsComplex(const double complex *a, int n, int lda) {
  if(n <= 0 || n % 2 != 0) return 0;
  for(int i=0; i<n; i+=2) {
    const double complex pivot = a[i+(size_t)(i+1)*lda];
    if(!isfinite(creal(pivot)) || !isfinite(cimag(pivot)) ||
       (creal(pivot) == 0.0 && cimag(pivot) == 0.0)) return 0;
  }
  return 1;
}

/* Preserve the INFO=0 arithmetic and publication order. On recovery,
 * publish PfM only after the inverse is finite as well. Failure leaves
 * the matrix workspace overwritten, as in the original implementation. */
static int fullSkewPfaffianInverseReal(int n, double *a, int lda,
    int *ipiv, double *work, double *buf, double *pfOut, int failure) {
  int info=0, nsq=n*n, one=1;
  double pfaff, minus_one=-1.0;
  M_DSKTRF("U", "N", &n, a, &lda, ipiv, buf, &nsq, &info);
  if(info < 0) return info;
  if(info > 0 && !fullSkewPairedPivotsReal(a,n,lda)) return info;

  utu2pfa_d(n,a,lda,ipiv,&pfaff);
  if(!isfinite(pfaff)) return failure;
  if(info == 0) *pfOut = pfaff;
  utu2inv_d(n,a,lda,ipiv,work,buf,lda);
  /* mVMC stores the inverse as row-major: transpose the skew inverse. */
  M_DSCAL(&nsq,&minus_one,a,&one);
  if(info > 0) {
    for(int i=0; i<nsq; ++i) if(!isfinite(a[i])) return failure;
    *pfOut = pfaff;
  }
  return 0;
}

static int fullSkewPfaffianInverseComplex(int n, double complex *a, int lda,
    int *ipiv, double complex *work, double complex *buf,
    double complex *pfOut, int failure) {
  int info=0, nsq=n*n, one=1;
  double complex pfaff, minus_one=-1.0;
  M_ZSKTRF("U", "N", &n, a, &lda, ipiv, buf, &nsq, &info);
  if(info < 0) return info;
  if(info > 0 && !fullSkewPairedPivotsComplex(a,n,lda)) return info;

  utu2pfa_z(n,a,lda,ipiv,&pfaff);
  if(info == 0) {
    /* Keep the historical normal-path check unchanged. */
    if(!isfinite(creal(pfaff) + cimag(pfaff))) return failure;
    *pfOut = pfaff;
  } else if(!isfinite(creal(pfaff)) || !isfinite(cimag(pfaff))) {
    return failure;
  }
  utu2inv_z(n,a,lda,ipiv,work,buf,lda);
  M_ZSCAL(&nsq,&minus_one,a,&one);
  if(info > 0) {
    for(int i=0; i<nsq; ++i) {
      if(!isfinite(creal(a[i])) || !isfinite(cimag(a[i]))) return failure;
    }
    *pfOut = pfaff;
  }
  return 0;
}

int getLWork() {
  char uplo='U', mthd='P';
  int n,lda,lwork,info=0;
  double pfaff;
  int iwork;
  double a;
  double optSize1,optSize2;

  /* ask the optimal size of work */
  n=lda=Nsize;
  lwork=-1;
  M_DGETRI(&n, &a, &lda, &iwork, &optSize1, &lwork, &info);
  lwork=-1;
  M_DSKPFA(&uplo, &mthd, &n, &a, &lda, &pfaff, &iwork, &optSize2, &lwork, &info);

  lwork = (optSize1>optSize2) ? (int)optSize1 : (int)optSize2;
  return lwork;
}

int getLWork_fcmp_dim(const int nDim) {
  char uplo='U', mthd='P';
  int n,lda,lwork,info=0;
  double rwork;
  double complex pfaff;
  int iwork;
  double complex a;
  double complex optSize1,optSize2;

  /* ask the optimal size of work */
  n=lda=nDim;
  lwork=-1;
  M_ZGETRI(&n, &a, &lda, &iwork, &optSize1, &lwork, &info);
  lwork=-1;
  M_ZSKPFA(&uplo, &mthd, &n, &a, &lda, &pfaff, &iwork, &optSize2, &lwork, &rwork, &info);

  lwork = (creal(optSize1)>creal(optSize2)) ? (int)creal(optSize1) : (int)creal(optSize2);
  return lwork;
}

int getLWork_fcmp() {
  return getLWork_fcmp_dim(Nsize);
}

//==============s fsz =============//
/* Calculate PfM and InvM from qpidx=qpStart to qpEnd */
int CalculateMAll_fsz(const int *eleIdx,const int *eleSpn, const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  int qpidx;

  int info = 0;

  double complex *myBufM;
  double complex *myWork;
  int            *myIWork;
  int             myInfo;
  double         *myRWork;

  RequestWorkSpaceThreadInt(Nsize);         //int

  RequestWorkSpaceThreadComplex(Nsize*Nsize+LapackLWork);

  RequestWorkSpaceThreadDouble(LapackLWork); // TBC for rwork

#pragma omp parallel default(shared)              \
  private(myIWork,myWork,myInfo,myBufM)
  {
    myIWork = GetWorkSpaceThreadInt(Nsize); // int

    myBufM  = GetWorkSpaceThreadComplex(Nsize*Nsize); //comp
    myWork  = GetWorkSpaceThreadComplex(LapackLWork); // comp

    myRWork = GetWorkSpaceThreadDouble(LapackLWork); //TBC for rwork

#pragma omp for private(qpidx)
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      if(info!=0) continue;

      myInfo = calculateMAll_child_fsz(eleIdx,eleSpn, qpStart, qpEnd, qpidx,
          myBufM, myIWork, myWork, LapackLWork,myRWork);
      if(myInfo!=0) {
#pragma omp critical
        info=myInfo;
      }
    }
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  ReleaseWorkSpaceThreadDouble();
  return info;
}

int calculateMAll_child_fsz(const int *eleIdx,const int *eleSpn, const int qpStart, const int qpEnd, const int qpidx,
    double complex *bufM, int *iwork, double complex *work, int lwork,double *rwork) {
#pragma procedure serial
  /* const int qpNum = qpEnd-qpStart; */
  int msi,msj;
  int rsi,rsj;

  int n,lda;
  //int nspn = 2*Ne+2*Nsite+2*Nsite+NProj; this is useful?

  /* optimization for Kei */
  const int nsize = Nsize;

  const double complex *sltE = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2;
  const double complex *sltE_i;

  double complex *invM = InvM + qpidx*Nsize*Nsize;
  double complex *invM_i, *invM_i2;

  n=lda=Nsize;

  /* store invM */
  /* Note that invM is column-major and skew-symmetric. */
  /* invM[msj][msi] = -sltE[rsi][rsj] */
#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    rsi = eleIdx[msi] + eleSpn[msi]*Nsite;//fsz
    invM_i = invM + msi*Nsize;
    sltE_i = sltE + rsi*Nsite2;
#pragma loop norecurrence
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj] + eleSpn[msj]*Nsite;//fsz
      invM_i[msj] = -sltE_i[rsj];
    }
  }

  return fullSkewPfaffianInverseComplex(n, invM, lda, iwork, work,
      bufM, &PfM[qpidx], qpidx+1);
}

/* Calculate PfM and InvM from qpidx=qpStart to qpEnd */
int CalculateMAll_fsz_real(const int *eleIdx,const int *eleSpn, const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  int qpidx;

  int info = 0;

  double *myBufM;
  double *myWork;
  int            *myIWork;
  int             myInfo;
  double         *myRWork;

  RequestWorkSpaceThreadInt(Nsize);         //int
  RequestWorkSpaceThreadDouble(Nsize*Nsize+LapackLWork);

#pragma omp parallel default(shared)              \
  private(myIWork,myWork,myInfo,myBufM)
  {
    myIWork = GetWorkSpaceThreadInt(Nsize); // int

    myBufM  = GetWorkSpaceThreadDouble(Nsize*Nsize); //comp
    myWork  = GetWorkSpaceThreadDouble(LapackLWork); // comp

#pragma omp for private(qpidx)
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      if(info!=0) continue;

      myInfo = calculateMAll_child_fsz_real(eleIdx,eleSpn, qpStart, qpEnd, qpidx,
               myBufM, myIWork, myWork, LapackLWork);
      if(myInfo!=0) {
#pragma omp critical
        info=myInfo;
      }
    }
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadDouble();
  return info;
}

int calculateMAll_child_fsz_real(const int *eleIdx,const int *eleSpn, const int qpStart, const int qpEnd, const int qpidx,
    double *bufM, int *iwork, double *work, int lwork) {
#pragma procedure serial
  /* const int qpNum = qpEnd-qpStart; */
  int msi,msj;
  int rsi,rsj;

  int n,lda;
  //int nspn = 2*Ne+2*Nsite+2*Nsite+NProj; this is useful?

  /* optimization for Kei */
  const int nsize = Nsize;

  const double *sltE = SlaterElm_real + (qpidx+qpStart)*Nsite2*Nsite2;
  const double *sltE_i;

  double *invM = InvM_real + qpidx*Nsize*Nsize;
  double *invM_i, *invM_i2;

  n=lda=Nsize;

  /* store invM */
  /* Note that invM is column-major and skew-symmetric. */
  /* invM[msj][msi] = -sltE[rsi][rsj] */
#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    rsi = eleIdx[msi] + eleSpn[msi]*Nsite;//fsz
    invM_i = invM + msi*Nsize;
    sltE_i = sltE + rsi*Nsite2;
#pragma loop norecurrence
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj] + eleSpn[msj]*Nsite;//fsz
      invM_i[msj] = -sltE_i[rsj];
    }
  }

  /* Calculate Pf M */
  return fullSkewPfaffianInverseReal(n, invM, lda, iwork, work,
      bufM, &PfM_real[qpidx], qpidx+1);
}




//==============s fcmp =============//
/* Calculate PfM and InvM from qpidx=qpStart to qpEnd */
int CalculateMAll_fcmp(const int *eleIdx, const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  int qpidx;

  int info = 0;

  double complex *myBufM;
  double complex *myWork;
  int            *myIWork;
  int             myInfo;
  double         *myRWork;

  RequestWorkSpaceThreadInt(Nsize);         //int

  RequestWorkSpaceThreadComplex(Nsize*Nsize+LapackLWork);

  RequestWorkSpaceThreadDouble(LapackLWork); // TBC for rwork

#pragma omp parallel default(shared)              \
  private(myIWork,myWork,myInfo,myBufM, myRWork) //TODO: Check to add myRWork is correct or not.
  {
    myIWork = GetWorkSpaceThreadInt(Nsize); // int

    myBufM  = GetWorkSpaceThreadComplex(Nsize*Nsize); //comp
    myWork  = GetWorkSpaceThreadComplex(LapackLWork); // comp

    myRWork = GetWorkSpaceThreadDouble(LapackLWork); //TBC for rwork

#pragma omp for private(qpidx)
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      if(info!=0) continue;

      myInfo = calculateMAll_child_fcmp(eleIdx, qpStart, qpEnd, qpidx,
          myBufM, myIWork, myWork, LapackLWork,myRWork);
      if(myInfo!=0) {
#pragma omp critical
        info=myInfo;
      }
    }
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  ReleaseWorkSpaceThreadDouble();
  return info;
}

int calculateMAll_child_fcmp(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
    double complex *bufM, int *iwork, double complex *work, int lwork,double *rwork) {
#pragma procedure serial
  /* const int qpNum = qpEnd-qpStart; */
  int msi,msj;
  int rsi,rsj;

  int n,lda;

  /* optimization for Kei */
  const int nsize = Nsize;

  const double complex *sltE = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2;
  const double complex *sltE_i;

  double complex *invM = InvM + qpidx*Nsize*Nsize;
  double complex *invM_i, *invM_i2;

  n=lda=Nsize;

  /* store invM */
  /* Note that invM is column-major and skew-symmetric. */
  /* invM[msj][msi] = -sltE[rsi][rsj] */
#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    rsi = eleIdx[msi] + (msi/Ne)*Nsite;
    invM_i = invM + msi*Nsize;
    sltE_i = sltE + rsi*Nsite2;
#pragma loop norecurrence
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj] + (msj/Ne)*Nsite;
      invM_i[msj] = -sltE_i[rsj];
    }
  }

  /* Calculate Pf M */
  return fullSkewPfaffianInverseComplex(n, invM, lda, iwork, work,
      bufM, &PfM[qpidx], qpidx+1);
}

int CalculateMAll_BF_fsz(const int *eleIdx, const int *eleSpn, const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  int qpidx;

  int info = 0;

  double complex *myBufM;
  double complex *myWork;
  int *myIWork;
  int myInfo;
  double *myRWork;

  RequestWorkSpaceThreadInt(Nsize);
  RequestWorkSpaceThreadComplex(Nsize*Nsize+LapackLWork);
  RequestWorkSpaceThreadDouble(LapackLWork);

#pragma omp parallel default(shared)              \
  private(myIWork,myWork,myRWork,myInfo,myBufM)
  {
    myIWork = GetWorkSpaceThreadInt(Nsize);
    myBufM  = GetWorkSpaceThreadComplex(Nsize*Nsize);
    myWork  = GetWorkSpaceThreadComplex(LapackLWork);
    myRWork = GetWorkSpaceThreadDouble(LapackLWork);

#pragma omp for private(qpidx)
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      if(info!=0) continue;

      myInfo = calculateMAll_BF_fsz_child(eleIdx, eleSpn, qpStart, qpEnd, qpidx,
          myBufM, myIWork, myWork, LapackLWork, myRWork, PfM, InvM);
      if(myInfo!=0) {
#pragma omp critical
        info=myInfo;
      }
    }
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  ReleaseWorkSpaceThreadDouble();
  return info;
}

int CalculatePfM_BF_fsz(const int *eleIdx, const int *eleSpn, const int qpStart, const int qpEnd,
    double complex *pfMOut, int *failureDetail) {
  return CalculatePfM_BF_fsz_from(SlaterElmBF, eleIdx, eleSpn, qpStart, qpEnd,
      pfMOut, failureDetail);
}

int CalculatePfM_BF_fsz_from(const double complex *sltElmBF, const int *eleIdx,
    const int *eleSpn, const int qpStart, const int qpEnd,
    double complex *pfMOut, int *failureDetail) {
  const int qpNum = qpEnd-qpStart;
  int qpidx;

  int status = BF_FSZ_PF_OK;
  int detail = 0;

  double complex *myBufM;
  double complex *myWork;
  int *myIWork;
  int myStatus;
  int myDetail;
  double *myRWork;

  RequestWorkSpaceThreadInt(Nsize);
  RequestWorkSpaceThreadComplex(Nsize*Nsize+LapackLWork);
  RequestWorkSpaceThreadDouble(LapackLWork);

#pragma omp parallel default(shared)              \
  private(myIWork,myWork,myRWork,myStatus,myDetail,myBufM)
  {
    myIWork = GetWorkSpaceThreadInt(Nsize);
    myBufM  = GetWorkSpaceThreadComplex(Nsize*Nsize);
    myWork  = GetWorkSpaceThreadComplex(LapackLWork);
    myRWork = GetWorkSpaceThreadDouble(LapackLWork);

#pragma omp for private(qpidx)
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      myDetail = 0;
      myStatus = calculatePfM_BF_fsz_child_from(sltElmBF, eleIdx, eleSpn,
          qpStart, qpidx, myBufM, myIWork, myWork, LapackLWork, myRWork,
          pfMOut, &myDetail);
      if(myStatus!=BF_FSZ_PF_OK) {
#pragma omp critical
        {
          if(myStatus > status) {
            status = myStatus;
            detail = myDetail;
          }
        }
      }
    }
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  ReleaseWorkSpaceThreadDouble();
  if(failureDetail != NULL) *failureDetail = detail;
  return status;
}

int CalculatePfM_BF_fsz_from_workspace(const double complex *sltElmBF,
    const int *eleIdx, const int *eleSpn, const int qpStart, const int qpEnd,
    double complex *pfMOut, int *failureDetail, double complex *bufM,
    int *iwork, double complex *work, int lwork, double *rwork) {
  if(qpStart == qpEnd && qpStart >= 0 && qpEnd <= NQPFull) {
    if(failureDetail != NULL) *failureDetail = 0;
    return BF_FSZ_PF_OK;
  }
  return CalculatePfM_BF_from_workspace(sltElmBF, eleIdx, eleSpn, qpStart,
      qpEnd, pfMOut, failureDetail, bufM, iwork, work, lwork, rwork);
}

int CalculatePfM_BF_from_workspace(const double complex *sltElmBF,
    const int *eleIdx, const int *eleSpn, int qpStart, int qpEnd,
    double complex *pfMOut, int *failureDetail, double complex *bufM,
    int *iwork, double complex *work, int lwork, double *rwork) {
  int qpNum;
  int qpidx;
  int status = BF_PF_OK;
  int detail = 0;

  if(failureDetail != NULL) *failureDetail = 0;
  if(sltElmBF == NULL || eleIdx == NULL || eleSpn == NULL || pfMOut == NULL
     || bufM == NULL || iwork == NULL || work == NULL || rwork == NULL
     || qpStart < 0 || qpStart >= qpEnd || qpEnd > NQPFull
     || Nsize <= 0 || Nsite2 <= 0 || lwork < LapackLWork) {
    return BF_PF_INVALID_ARGUMENT;
  }

  qpNum = qpEnd-qpStart;
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    int myDetail = 0;
    int myStatus = calculatePfM_BF_fsz_child_from(sltElmBF, eleIdx, eleSpn,
        qpStart, qpidx, bufM, iwork, work, lwork, rwork, pfMOut, &myDetail);
    if(myStatus > status) {
      status = myStatus;
      detail = myDetail;
    }
  }

  if(failureDetail != NULL) *failureDetail = detail;
  return status;
}

static BF_FSZ_MAllResult BF_FSZ_MAllResultValue(const int status,
    const int stage, const int qpidx, const int lapackInfo) {
  BF_FSZ_MAllResult result;
  result.status = status;
  result.stage = stage;
  result.qpidx = qpidx;
  result.lapackInfo = lapackInfo;
  return result;
}

static BF_FSZ_MAllResult calculateMAll_BF_fsz_child_from(
    const double complex *sltElmBF, const int *eleIdx, const int *eleSpn,
    const int qpStart, const int qpidx, double complex *pfMOut,
    double complex *invMOut, const size_t invMQpStride,
    double complex *bufM, int *iwork, double complex *work,
    int lwork, double *rwork) {
  const int nsize = Nsize;
  const int globalQpidx = qpStart+qpidx;
  const double complex *sltE = sltElmBF
      + (size_t)globalQpidx*(size_t)Nsite2*(size_t)Nsite2;
  double complex *invM = invMOut+(size_t)qpidx*invMQpStride;
  char uplo='U',mthd='P';
  int m=nsize,n=nsize,lda=nsize,info=0;
  int i,j;
  double complex pfaff;

  for(i=0;i<nsize;i++) {
    const int rsi = eleIdx[i]+eleSpn[i]*Nsite;
    const double complex *sltE_i = sltE+(size_t)rsi*(size_t)Nsite2;
    double complex *bufM_i = bufM+(size_t)i*(size_t)nsize;
    for(j=0;j<nsize;j++) {
      const int rsj = eleIdx[j]+eleSpn[j]*Nsite;
      bufM_i[j] = -sltE_i[rsj];
    }
  }
  memcpy(invM,bufM,sizeof(double complex)*(size_t)nsize*(size_t)nsize);

  M_ZSKPFA(&uplo,&mthd,&n,bufM,&lda,&pfaff,iwork,work,&lwork,rwork,&info);
  if(info != 0) {
    return BF_FSZ_MAllResultValue(BF_FSZ_MALL_LAPACK_FAILURE,
        BF_FSZ_MALL_STAGE_PFAFFIAN,globalQpidx,info);
  }
  if(!(isfinite(creal(pfaff)) && isfinite(cimag(pfaff)))) {
    return BF_FSZ_MAllResultValue(BF_FSZ_MALL_NONFINITE,
        BF_FSZ_MALL_STAGE_PFAFFIAN,globalQpidx,0);
  }

  M_ZGETRF(&m,&n,invM,&lda,iwork,&info);
  if(info != 0) {
    return BF_FSZ_MAllResultValue(BF_FSZ_MALL_LAPACK_FAILURE,
        BF_FSZ_MALL_STAGE_GETRF,globalQpidx,info);
  }
  M_ZGETRI(&n,invM,&lda,iwork,work,&lwork,&info);
  if(info != 0) {
    return BF_FSZ_MAllResultValue(BF_FSZ_MALL_LAPACK_FAILURE,
        BF_FSZ_MALL_STAGE_GETRI,globalQpidx,info);
  }
  for(i=0;i<nsize*nsize;i++) {
    invM[i] = -invM[i];
    if(!(isfinite(creal(invM[i])) && isfinite(cimag(invM[i])))) {
      return BF_FSZ_MAllResultValue(BF_FSZ_MALL_NONFINITE,
          BF_FSZ_MALL_STAGE_INVERSE,globalQpidx,0);
    }
  }
  pfMOut[qpidx] = pfaff;
  return BF_FSZ_MAllResultValue(BF_FSZ_MALL_OK,BF_FSZ_MALL_STAGE_NONE,
      globalQpidx,0);
}

BF_FSZ_MAllResult CalculateMAll_BF_fsz_from_workspace(
    const double complex *sltElmBF, const int *eleIdx, const int *eleSpn,
    const int qpStart, const int qpEnd, double complex *pfMOut,
    double complex *invMOut, const size_t invMQpStride,
    double complex *bufM, int *iwork, double complex *work,
    int lwork, double *rwork) {
  const size_t nsizeSquared = (size_t)Nsize*(size_t)Nsize;
  const int qpNum = qpEnd-qpStart;
  int qpidx,i;

  if(qpStart < 0 || qpEnd < qpStart || qpEnd > NQPFull) {
    return BF_FSZ_MAllResultValue(BF_FSZ_MALL_INVALID_ARGUMENT,
        BF_FSZ_MALL_STAGE_NONE,qpStart,0);
  }
  if(qpNum == 0) {
    return BF_FSZ_MAllResultValue(BF_FSZ_MALL_OK,BF_FSZ_MALL_STAGE_NONE,
        qpStart,0);
  }
  if(sltElmBF == NULL || eleIdx == NULL || eleSpn == NULL || pfMOut == NULL
      || invMOut == NULL || bufM == NULL || iwork == NULL || work == NULL
      || rwork == NULL || lwork < Nsize || invMQpStride < nsizeSquared
      || invMQpStride > (size_t)PTRDIFF_MAX
      || (qpNum > 1 && invMQpStride > SIZE_MAX/(size_t)(qpNum-1))) {
    return BF_FSZ_MAllResultValue(BF_FSZ_MALL_INVALID_ARGUMENT,
        BF_FSZ_MALL_STAGE_NONE,qpStart,0);
  }
  for(i=0;i<Nsize;i++) {
    if(eleIdx[i] < 0 || eleIdx[i] >= Nsite
        || (eleSpn[i] != 0 && eleSpn[i] != 1)) {
      return BF_FSZ_MAllResultValue(BF_FSZ_MALL_INVALID_ARGUMENT,
          BF_FSZ_MALL_STAGE_NONE,qpStart,0);
    }
  }
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    const BF_FSZ_MAllResult result = calculateMAll_BF_fsz_child_from(
        sltElmBF,eleIdx,eleSpn,qpStart,qpidx,pfMOut,invMOut,invMQpStride,
        bufM,iwork,work,lwork,rwork);
    if(result.status != BF_FSZ_MALL_OK) return result;
  }
  return BF_FSZ_MAllResultValue(BF_FSZ_MALL_OK,BF_FSZ_MALL_STAGE_NONE,
      qpStart,0);
}

int calculatePfM_BF_fsz_child_from(
       const double complex *sltElmBF,
       const int *eleIdx,
       const int *eleSpn,
       const int qpStart,
       const int qpidx,
       double complex *bufM,
       int *iwork,
       double complex *work,
       int lwork,
       double *rwork,
       double complex *pfMOut,
       int *failureDetail
       )
{
#pragma procedure serial
  int msi,msj;
  int rsi,rsj;

  char uplo='U', mthd='P';
  int n,lda,info=0;
  double complex pfaff;

  const int nsize = Nsize;

  const double complex *sltE = sltElmBF + (qpidx+qpStart)*Nsite2*Nsite2;
  const double complex *sltE_i;

  double complex *bufM_i;

  n=lda=Nsize;

#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    rsi = eleIdx[msi] + eleSpn[msi]*Nsite;
    bufM_i = bufM + msi*Nsize;
    sltE_i = sltE + rsi*Nsite2;
#pragma loop norecurrence
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj] + eleSpn[msj]*Nsite;
      bufM_i[msj] = -sltE_i[rsj];
    }
  }

  M_ZSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff, iwork, work, &lwork, rwork, &info);
  if(info!=0) {
    if(failureDetail != NULL) *failureDetail = info;
    return info < 0 ? BF_FSZ_PF_INVALID_ARGUMENT : BF_FSZ_PF_LAPACK_FAILURE;
  }
  if(!(isfinite(creal(pfaff)) && isfinite(cimag(pfaff)))) {
    if(failureDetail != NULL) *failureDetail = qpidx+qpStart;
    return BF_FSZ_PF_NONFINITE;
  }
  pfMOut[qpidx] = pfaff;
  return BF_FSZ_PF_OK;
}

int calculateMAll_BF_fsz_child(
       const int *eleIdx,
       const int *eleSpn,
       const int qpStart,
       const int qpEnd,
       const int qpidx,
       double complex *bufM,
       int *iwork,
       double complex *work,
       int lwork,
       double *rwork,
       double complex *pfM,
       double complex *invMAll
       )
{
#pragma procedure serial
  int msi,msj;
  int rsi,rsj;

  char uplo='U', mthd='P';
  int m,n,lda,nsq,one,info=0;
  double complex pfaff,minus_one;

  const int nsize = Nsize;

  const double complex *sltE = SlaterElmBF + (qpidx+qpStart)*Nsite2*Nsite2;
  const double complex *sltE_i;

  double complex *invM = invMAll + qpidx*Nsize*Nsize;
  double complex *bufM_i;

  (void)qpEnd;

  m=n=lda=Nsize;
  nsq=n*n;
  one=1;
  minus_one=-1.0;

#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    rsi = eleIdx[msi] + eleSpn[msi]*Nsite;
    bufM_i = bufM + msi*Nsize;
    sltE_i = sltE + rsi*Nsite2;
#pragma loop norecurrence
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj] + eleSpn[msj]*Nsite;
      bufM_i[msj] = -sltE_i[rsj];
    }
  }

  for(msi=0;msi<nsize*nsize;msi++) invM[msi] = bufM[msi];

  M_ZSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff, iwork, work, &lwork, rwork, &info);

  if(info!=0) return info;
  if(!(isfinite(creal(pfaff)) && isfinite(cimag(pfaff)))) return qpidx+1;
  pfM[qpidx] = pfaff;

  M_ZGETRF(&m, &n, invM, &lda, iwork, &info);
  if(info!=0) return info;
  M_ZGETRI(&n, invM, &lda, iwork, work, &lwork, &info);
  if(info!=0) return info;

  M_ZSCAL(&nsq, &minus_one, invM, &one);

  return info;
}

//==============e fcmp =============//

//==============s real =============//
/* Calculate PfM and InvM from qpidx=qpStart to qpEnd */
int CalculateMAll_real(const int *eleIdx, const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  int qpidx;

  int info = 0;

  double *myBufM;
  double *myWork;
  int *myIWork;
  int myInfo;

  RequestWorkSpaceThreadInt(Nsize);
  RequestWorkSpaceThreadDouble(Nsize*Nsize+LapackLWork);

#pragma omp parallel default(shared)              \
  private(myIWork,myWork,myInfo,myBufM)
  {
    myIWork = GetWorkSpaceThreadInt(Nsize);
    myBufM = GetWorkSpaceThreadDouble(Nsize*Nsize);
    myWork = GetWorkSpaceThreadDouble(LapackLWork);
#pragma omp for private(qpidx)
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      if(info!=0) continue;

      //      myInfo = calculateMAll_child_real(eleIdx, qpStart, qpEnd, qpidx,
      //                                  myBufM, myIWork, myWork, LapackLWork);
      myInfo = calculateMAll_child_real(eleIdx, qpStart, qpEnd, qpidx,
                                        myBufM, myIWork, myWork, LapackLWork, PfM_real, InvM_real);

      if(myInfo!=0) {
#pragma omp critical
        info=myInfo;
      }
    }
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadDouble();
  return info;
}

int calculateMAll_child_real(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
    double *bufM, int *iwork, double *work, int lwork, double* PfM_real, double *InvM_real) {
#pragma procedure serial
  /* const int qpNum = qpEnd-qpStart; */
  int msi,msj;
  int rsi,rsj;

  int n,lda;

  /* optimization for Kei */
  const int nsize = Nsize;

  const double *sltE = SlaterElm_real + (qpidx+qpStart)*Nsite2*Nsite2;
  const double *sltE_i;

  double *invM = InvM_real + qpidx*Nsize*Nsize;
  double *invM_i, *invM_i2;

  n=lda=Nsize;

  /* store invM */
  /* Note that invM is column-major and skew-symmetric. */
  /* invM[msj][msi] = -sltE[rsi][rsj] */
#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    rsi = eleIdx[msi] + (msi/Ne)*Nsite;
    invM_i = invM + msi*Nsize;
    sltE_i = sltE + rsi*Nsite2;
#pragma loop norecurrence
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj] + (msj/Ne)*Nsite;
      invM_i[msj] = -sltE_i[rsj];

    }
  }

  /* calculate Pf M */
  return fullSkewPfaffianInverseReal(n, invM, lda, iwork, work,
      bufM, &PfM_real[qpidx], qpidx+1);
}


#endif
