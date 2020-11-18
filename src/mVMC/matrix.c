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
#include "./include/global.h"
#include "workspace.c"

//int calculateMAll_child_real(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
//                        double *bufM, int *iwork, double *work, int lwork);
#define D_PfLimit 1.0e-100

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

int getLWork_fcmp() {
  char uplo='U', mthd='P';
  int n,lda,lwork,info=0;
  double rwork;
  double complex pfaff;
  int iwork;
  double complex a;
  double complex optSize1,optSize2;

  /* ask the optimal size of work */
  n=lda=Nsize;
  lwork=-1;
  M_ZGETRI(&n, &a, &lda, &iwork, &optSize1, &lwork, &info);
  lwork=-1;
#ifdef _pfaffine
  M_ZSKPFA(&uplo, &mthd, &n, &a, &lda, &pfaff, &iwork, &optSize2, &lwork/*, &rwork*/, &info);
#else
  M_ZSKPFA(&uplo, &mthd, &n, &a, &lda, &pfaff, &iwork, &optSize2, &lwork, &rwork, &info);
#endif

  lwork = (creal(optSize1)>creal(optSize2)) ? (int)creal(optSize1) : (int)creal(optSize2);
  return lwork;
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

  char uplo='U', mthd='P';
  int m,n,nsq,one,lda,info=0;
  //int nspn = 2*Ne+2*Nsite+2*Nsite+NProj; this is useful?
  double complex pfaff,minus_one;

  /* optimization for Kei */
  const int nsize = Nsize;

  const double complex *sltE = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2;
  const double complex *sltE_i;

  double complex *invM = InvM + qpidx*Nsize*Nsize;
  double complex *invM_i;

  double complex *bufM_i, *bufM_i2;

  m=n=lda=Nsize;
  nsq=n*n;
  one=1;
  minus_one=-1.0;

  /* store bufM */
  /* Note that bufM is column-major and skew-symmetric. */
  /* bufM[msj][msi] = -sltE[rsi][rsj] */
#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    rsi = eleIdx[msi] + eleSpn[msi]*Nsite;//fsz
    bufM_i = bufM + msi*Nsize;
    sltE_i = sltE + rsi*Nsite2;
#pragma loop norecurrence
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj] + eleSpn[msj]*Nsite;//fsz
      bufM_i[msj] = -sltE_i[rsj];
    }
  }

#ifdef _pfaffine
  info=0; /* Fused Pfaffian/inverse computation. */
  M_ZSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff, iwork, work/*, rwork*/, &lwork, &info);
#else
  /* Pfaffian/inverse computed separately. */
  /* Copy bufM to invM before using bufM to compute Pfaffian. */
  for(msi=0;msi<nsize*nsize;msi++)
    invM[msi] = bufM[msi];

  /* Calculate Pf M */
  M_ZSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff, iwork, work, &lwork, rwork, &info);
#endif

  if(info!=0) return info;
  if(!isfinite(creal(pfaff) + cimag(pfaff))) return qpidx+1;
  PfM[qpidx] = pfaff;

#ifdef _pfaffine
  /* For fused Pfaffian/inverse, inv(M) is already stored in bufM.
   * Now to transpose (.* -1) it to invM. */
  for(msi=0;msi<nsize*nsize;msi++)
    invM[msi] = -bufM[msi];
#else
  /* Calculate inverse. */
  M_ZGETRF(&m, &n, invM, &lda, iwork, &info); /* ipiv = iwork */
  M_ZGETRI(&n, invM, &lda, iwork, work, &lwork, &info);

  /* InvM -> InvM(T) -> -InvM */
  M_ZSCAL(&nsq, &minus_one, invM, &one);
#endif

  return info;
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

  char uplo='U', mthd='P';
  int m,n,lda,one,nsq,info=0;
  //int nspn = 2*Ne+2*Nsite+2*Nsite+NProj; this is useful?
  double pfaff,minus_one;

  /* optimization for Kei */
  const int nsize = Nsize;

  const double *sltE = SlaterElm_real + (qpidx+qpStart)*Nsite2*Nsite2;
  const double *sltE_i;

  double *invM = InvM_real + qpidx*Nsize*Nsize;
  double *invM_i;

  double *bufM_i, *bufM_i2;

  m=n=lda=Nsize;
  nsq=n*n;
  one=1;
  minus_one=-1.0;

  /* store bufM */
  /* Note that bufM is column-major and skew-symmetric. */
  /* bufM[msj][msi] = -sltE[rsi][rsj] */
#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    rsi = eleIdx[msi] + eleSpn[msi]*Nsite;//fsz
    bufM_i = bufM + msi*Nsize;
    sltE_i = sltE + rsi*Nsite2;
#pragma loop norecurrence
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj] + eleSpn[msj]*Nsite;//fsz
      bufM_i[msj] = -sltE_i[rsj];
    }
  }

#ifdef _pfaffine
  info=0; /* Fused Pfaffian/inverse computation. */
#else
  /* Pfaffian/inverse computed separately. */
  /* Copy bufM to invM before using bufM to compute Pfaffian. */
  for(msi=0;msi<nsize*nsize;msi++)
    invM[msi] = bufM[msi];
#endif

  /* Calculate Pf M */
  M_DSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff, iwork, work, &lwork, &info);

  if(info!=0) return info;
  if(!isfinite(pfaff)) return qpidx+1;
  PfM_real[qpidx] = pfaff;

#ifdef _pfaffine
  /* inv(M) already stored in bufM.
   * Transpose (.* -1) to invM. */
  for(msi=0;msi<nsize*nsize;msi++)
    invM[msi] = -bufM[msi];
#else
  /* Calculate inverse. */
  M_DGETRF(&m, &n, invM, &lda, iwork, &info); /* ipiv = iwork */
  if(info!=0) return info;
  M_DGETRI(&n, invM, &lda, iwork, work, &lwork, &info);
  if(info!=0) return info;

  /* InvM -> InvM(T) -> -InvM */
  M_DSCAL(&nsq, &minus_one, invM, &one);
#endif

  return info;
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

  char uplo='U', mthd='P';
  int m,n,nsq,lda,one,info=0;
  double complex pfaff,minus_one;

  /* optimization for Kei */
  const int nsize = Nsize;

  const double complex *sltE = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2;
  const double complex *sltE_i;

  double complex *invM = InvM + qpidx*Nsize*Nsize;
  double complex *invM_i;

  double complex *bufM_i, *bufM_i2;

  m=n=lda=Nsize;
  nsq=n*n;
  one=1;
  minus_one=-1.0;

  /* store bufM */
  /* Note that bufM is column-major and skew-symmetric. */
  /* bufM[msj][msi] = -sltE[rsi][rsj] */
#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    rsi = eleIdx[msi] + (msi/Ne)*Nsite;
    bufM_i = bufM + msi*Nsize;
    sltE_i = sltE + rsi*Nsite2;
#pragma loop norecurrence
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj] + (msj/Ne)*Nsite;
      bufM_i[msj] = -sltE_i[rsj];
    }
  }

#ifdef _pfaffine
  info=0; /* Fused Pfaffian/inverse computation. */
  M_ZSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff, iwork, work/*, rwork*/, &lwork, &info);
#else
  /* Pfaffian/inverse computed separately. */
  /* Copy bufM to invM before using bufM to compute Pfaffian. */
  for(msi=0;msi<nsize*nsize;msi++)
    invM[msi] = bufM[msi];

  /* Calculate Pf M */
  M_ZSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff, iwork, work, &lwork, rwork, &info);
#endif

  if(info!=0) return info;
  if(!isfinite(creal(pfaff) + cimag(pfaff))) return qpidx+1;
  PfM[qpidx] = pfaff;

#ifdef _pfaffine
  /* inv(M) already stored in bufM.
   * Transpose (.* -1) to invM. */
  for(msi=0;msi<nsize*nsize;msi++)
    invM[msi] = -bufM[msi];
#else
  /* Calculate inverse. */
  M_ZGETRF(&m, &n, invM, &lda, iwork, &info); /* ipiv = iwork */
  M_ZGETRI(&n, invM, &lda, iwork, work, &lwork, &info);

  /* mVMC's handling InvM as row-major,
   * i.e. InvM needs a transpose, InvM -> -InvM according antisymmetric properties. */
  M_ZSCAL(&nsq, &minus_one, invM, &one);
#endif

  return info;
}

int CalculateMAll_BF_fcmp(const int *eleIdx, const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  int qpidx;

  int info = 0;

  double complex*myBufM;
  double complex*myWork;
  int *myIWork;
  int myInfo;
  double *myRWork;

  RequestWorkSpaceThreadInt(Nsize);
  RequestWorkSpaceThreadComplex(Nsize*Nsize+LapackLWork);

  RequestWorkSpaceThreadDouble(LapackLWork); // TBC for rwork

#pragma omp parallel default(shared)              \
  private(myIWork,myWork,myRWork, myInfo,myBufM) //TODO: Check to add myRWork is correct or not.
  {
    myIWork = GetWorkSpaceThreadInt(Nsize);
    myBufM = GetWorkSpaceThreadComplex(Nsize*Nsize);
    myWork = GetWorkSpaceThreadComplex(LapackLWork);

    myRWork = GetWorkSpaceThreadDouble(LapackLWork); //TBC for rwork

#pragma omp for private(qpidx)
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      if(info!=0) continue;

      myInfo = calculateMAll_BF_fcmp_child(eleIdx, qpStart, qpEnd, qpidx,
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

int calculateMAll_BF_fcmp_child(
       const int *eleIdx,
       const int qpStart,
       const int qpEnd,
       const int qpidx,
       double complex*bufM,
       int *iwork,
       double complex*work,
       int lwork,
       double* rwork,
       double complex* PfM,
       double complex* InvM
       )
{
#pragma procedure serial
  /* const int qpNum = qpEnd-qpStart; */
  int msi,msj;
  int rsi,rsj;

  char uplo='U', mthd='P';
  int m,n,lda,nsq,one,info=0;
  double complex pfaff,minus_one;

  /* optimization for Kei */
  const int nsize = Nsize;

  const double complex *sltE = SlaterElmBF + (qpidx+qpStart)*Nsite2*Nsite2;
  const double complex *sltE_i;

  double complex*invM = InvM + qpidx*Nsize*Nsize;
  double complex*invM_i;

  double complex*bufM_i, *bufM_i2;

  m=n=lda=Nsize;
  nsq=n*n;
  one=1;
  minus_one=-1.0;

  /* store bufM */
  /* Note that bufM is column-major and skew-symmetric. */
  /* bufM[msj][msi] = -sltE[rsi][rsj] */
#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    rsi = eleIdx[msi] + (msi/Ne)*Nsite;
    bufM_i = bufM + msi*Nsize;
    sltE_i = sltE + rsi*Nsite2;
#pragma loop norecurrence
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj] + (msj/Ne)*Nsite;
      bufM_i[msj] = -sltE_i[rsj];
    }
  }

#ifdef _pfaffine
  info=0; /* Fused Pfaffian/inverse computation. */
  M_ZSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff, iwork, work, &lwork/*, rwork*/, &info);
#else
  /* Pfaffian/inverse computed separately. */
  /* Copy bufM to invM before using bufM to compute Pfaffian. */
  for(msi=0;msi<nsize*nsize;msi++)
    invM[msi] = bufM[msi];

  /* Calculate Pf M */
  M_ZSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff, iwork, work, &lwork, rwork, &info);
#endif

  if(info!=0) return info;
  if(!(isfinite(creal(pfaff)) && isfinite(cimag(pfaff)))) return qpidx+1;
  PfM[qpidx] = pfaff;

#ifdef _pfaffine
  /* inv(M) already stored in bufM.
   * Transpose (.* -1) to invM. */
  for(msi=0;msi<nsize*nsize;msi++)
    invM[msi] = -bufM[msi];
#else
  /* Calculate inverse. */
  M_ZGETRF(&m, &n, invM, &lda, iwork, &info); /* ipiv = iwork */
  if(info!=0) return info;
  M_ZGETRI(&n, invM, &lda, iwork, work, &lwork, &info);
  if(info!=0) return info;

  /* InvM -> InvM(T) -> -InvM */
  M_ZSCAL(&nsq, &minus_one, invM, &one);
#endif

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

//int calculateMAll_child_real(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
//                        double *bufM, int *iwork, double *work, int lwork) {
int calculateMAll_child_real(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
    double *bufM, int *iwork, double *work, int lwork, double* PfM_real, double *InvM_real) {
#pragma procedure serial
  /* const int qpNum = qpEnd-qpStart; */
  int msi,msj;
  int rsi,rsj;

  char uplo='U', mthd='P';
  int m,n,nsq,one,lda,info=0;
  double pfaff,minus_one;

  /* optimization for Kei */
  const int nsize = Nsize;

  const double *sltE = SlaterElm_real + (qpidx+qpStart)*Nsite2*Nsite2;
  const double *sltE_i;

  double *invM = InvM_real + qpidx*Nsize*Nsize;
  double *invM_i;

  double *bufM_i, *bufM_i2;

  m=n=lda=Nsize;
  nsq=n*n;
  one=1;
  minus_one=-1.0;

  /* store bufM */
  /* Note that bufM is column-major and skew-symmetric. */
  /* bufM[msj][msi] = -sltE[rsi][rsj] */
#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    rsi = eleIdx[msi] + (msi/Ne)*Nsite;
    bufM_i = bufM + msi*Nsize;
    sltE_i = sltE + rsi*Nsite2;
#pragma loop norecurrence
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj] + (msj/Ne)*Nsite;
      bufM_i[msj] = -sltE_i[rsj];

    }
  }

#ifdef _pfaffine
  info=0; /* Fused Pfaffian/inverse computation. */
#else
  /* Pfaffian/inverse computed separately. */
  /* Copy bufM to invM before using bufM to compute Pfaffian. */
  for(msi=0;msi<nsize*nsize;msi++)
    invM[msi] = bufM[msi];
#endif

  /* calculate Pf M */
  M_DSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff, iwork, work, &lwork, &info);

  if(info!=0) return info;
  if(!isfinite(pfaff)) return qpidx+1;
  PfM_real[qpidx] = pfaff;

#ifdef _pfaffine
  /* inv(M) already stored in bufM.
   * Transpose (.* -1) to invM. */
  for(msi=0;msi<nsize*nsize;msi++)
    invM[msi] = -bufM[msi];
#else
  /* Compute inverse */
  M_DGETRF(&m, &n, invM, &lda, iwork, &info); /* ipiv = iwork */
  M_DGETRI(&n, invM, &lda, iwork, work, &lwork, &info);

  // InvM -> InvM' = -InvM
  M_DSCAL(&nsq, &minus_one, invM, &one);
#endif

  return info;
}

int CalculateMAll_BF_real(const int *eleIdx, const int qpStart, const int qpEnd){
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

      myInfo = calculateMAll_BF_real_child(eleIdx, qpStart, qpEnd, qpidx,
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

int calculateMAll_BF_real_child(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
    double *bufM, int *iwork, double *work, int lwork, double* PfM_real, double* InvM_real) {
#pragma procedure serial
  /* const int qpNum = qpEnd-qpStart; */
  int msi,msj;
  int rsi,rsj;

  char uplo='U', mthd='P';
  int m,n,nsq,one,lda,info=0;
  double pfaff,minus_one;

  /* optimization for Kei */
  const int nsize = Nsize;

  const double *sltE = SlaterElmBF_real + (qpidx+qpStart)*Nsite2*Nsite2;
  const double *sltE_i;

  double *invM = InvM_real + qpidx*Nsize*Nsize;
  double *invM_i;

  double *bufM_i, *bufM_i2;

  m=n=lda=Nsize;
  nsq=n*n;
  one=1;
  minus_one=-1.0;

  /* store bufM */
  /* Note that bufM is column-major and skew-symmetric. */
  /* bufM[msj][msi] = -sltE[rsi][rsj] */
#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    rsi = eleIdx[msi] + (msi/Ne)*Nsite;
    bufM_i = bufM + msi*Nsize;
    sltE_i = sltE + rsi*Nsite2;
#pragma loop norecurrence
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj] + (msj/Ne)*Nsite;
      bufM_i[msj] = -sltE_i[rsj];
    }
  }

#ifdef _pfaffine
  info=0; /* Fused Pfaffian/inverse computation. */
#else
  /* Pfaffian/inverse computed separately. */
  /* Copy bufM to invM before using bufM to compute Pfaffian. */
  for(msi=0;msi<nsize*nsize;msi++)
    invM[msi] = bufM[msi];
#endif

  /* Calculate Pf M */
  M_DSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff, iwork, work, &lwork, &info);

  if(info!=0) return info;
  if(!isfinite(pfaff)) return qpidx+1;
  PfM_real[qpidx] = pfaff;

#ifdef _pfaffine
  /* inv(M) already stored in bufM.
   * Transpose (.* -1) to invM. */
  for(msi=0;msi<nsize*nsize;msi++)
    invM[msi] = -bufM[msi];
#else
  /* Compute inverse. */
  M_DGETRF(&m, &n, invM, &lda, iwork, &info); /* ipiv = iwork */
  if(info!=0) return info;
  M_DGETRI(&n, invM, &lda, iwork, work, &lwork, &info);
  if(info!=0) return info;

  // InvM -> InvM' = -InvM
  M_DSCAL(&nsq, &minus_one, invM, &one);
#endif

  return info;
}

//==============e real =============//

#endif
