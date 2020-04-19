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
  M_ZSKPFA(&uplo, &mthd, &n, &a, &lda, &pfaff, &iwork, &optSize2, &lwork/*, &rwork*/, &info);

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
      //printf("DEBUG: msi=%d msj=%d: rsi=%d rsj=%d :bufM=%lf %lf \n",msi,msj,rsi,rsj,creal(bufM_i[msj]),cimag(bufM_i[msj]));
    }
  }

  /* copy bufM to invM */
  /* For Pfaffian calculation, invM is used as second buffer */
#pragma loop noalias
  for(msi=0;msi<nsize*nsize;msi++) {
    invM[msi] = bufM[msi];
  }
  /* calculate Pf M */
  //printf("DEBUG: n=%d \n",n);
  //for(msi=0;msi<nsize;msi++){
  //  for(msj=0;msj<nsize;msj++){
  //    printf("DEBUG: msi=%d msj=%d bufM %lf %lf \n",msi,msj,creal(bufM[msi+msj*n]),cimag(bufM[msi+msj*n]));
  //  }
  //}
  M_ZSKPFA(&uplo, &mthd, &n, invM, &lda, &pfaff, iwork, work/*, rwork*/, &lwork, &info);
  //printf("DEBUG: pfaff=%lf %lf\n",creal(pfaff),cimag(pfaff));
  if(info!=0) return info;
  if(!isfinite(creal(pfaff) + cimag(pfaff))) return qpidx+1;
  PfM[qpidx] = pfaff;

  /* DInv */
  // M_ZGETRF(&m, &n, bufM, &lda, iwork, &info); /* ipiv = iwork */
  // M_ZGETRI(&n, bufM, &lda, iwork, work, &lwork, &info);

  // InvM -> InvM(T) -> -InvM
  M_ZSCAL(&nsq, &minus_one, InvM, &one);

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
  int m,n,lda,info=0;
  //int nspn = 2*Ne+2*Nsite+2*Nsite+NProj; this is useful?
  double pfaff;

  /* optimization for Kei */
  const int nsize = Nsize;

  const double *sltE = SlaterElm_real + (qpidx+qpStart)*Nsite2*Nsite2;
  const double *sltE_i;

  double *invM = InvM_real + qpidx*Nsize*Nsize;
  double *invM_i;

  double *bufM_i, *bufM_i2;

  m=n=lda=Nsize;

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

  /* copy bufM to invM */
  /* For Pfaffian calculation, invM is used as second buffer */
#pragma loop noalias
  for(msi=0;msi<nsize*nsize;msi++) {
    invM[msi] = bufM[msi];
  }
  /* calculate Pf M */
  M_DSKPFA(&uplo, &mthd, &n, invM, &lda, &pfaff, iwork, work, &lwork, &info);
  if(info!=0) return info;
  if(!isfinite(pfaff)) return qpidx+1;
  PfM_real[qpidx] = pfaff;

  /* DInv */
  M_DGETRF(&m, &n, bufM, &lda, iwork, &info); /* ipiv = iwork */
  if(info!=0) return info;

  M_DGETRI(&n, bufM, &lda, iwork, work, &lwork, &info);
  if(info!=0) return info;

  /* store InvM */
  /* BufM is column-major, InvM is row-major */
#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    invM_i = invM + msi*Nsize;
    bufM_i = bufM + msi*Nsize;
    bufM_i2 = bufM + msi;
    for(msj=0;msj<nsize;msj++) {
      invM_i[msj] = 0.5*(bufM_i2[msj*nsize] - bufM_i[msj]);
      /* invM[i][j] = 0.5*(bufM[i][j]-bufM[j][i]) */
    }
  }

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
      //      printf("DEBUG: msi=%d msj=%d bufM=%lf %lf \n",msi,msj,creal(bufM_i[msj]),cimag(bufM_i[msj]));
    }
  }

  /* copy bufM to invM */
  /* For Pfaffian calculation, invM is used as second buffer */
#pragma loop noalias
  for(msi=0;msi<nsize*nsize;msi++) {
    invM[msi] = bufM[msi];
  }
  // [TODO] Now bufM is not needed. Remove.
  /* calculate Pf M */
  //printf("DEBUG: n=%d \n",n);
  //for(msi=0;msi<nsize;msi++){
  //  for(msj=0;msj<nsize;msj++){
  //    printf("DEBUG: msi=%d msj=%d bufM %lf %lf \n",msi,msj,creal(bufM[msi+msj*n]),cimag(bufM[msi+msj*n]));
  //  }
  //}
  M_ZSKPFA(&uplo, &mthd, &n, invM, &lda, &pfaff, iwork, work/*, rwork*/, &lwork, &info);
  // printf("DEBUG: pfaff=%lf %lf\n",creal(pfaff),cimag(pfaff));
  if(info!=0) return info;
  if(!isfinite(creal(pfaff) + cimag(pfaff))) return qpidx+1;
  PfM[qpidx] = pfaff;

  /* DInv */
  // M_ZGETRF(&m, &n, bufM, &lda, iwork, &info); /* ipiv = iwork */
  // M_ZGETRI(&n, bufM, &lda, iwork, work, &lwork, &info);

  /* mVMC's handling InvM as row-major,
   * i.e. InvM needs a transpose, InvM -> -InvM according antisymmetric properties. */
  M_ZSCAL(&nsq, &minus_one, InvM, &one);

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
  int m,n,lda,info=0;
  double complex pfaff;

  /* optimization for Kei */
  const int nsize = Nsize;

  const double complex *sltE = SlaterElmBF + (qpidx+qpStart)*Nsite2*Nsite2;
  const double complex *sltE_i;

  double complex*invM = InvM + qpidx*Nsize*Nsize;
  double complex*invM_i;

  double complex*bufM_i, *bufM_i2;

  m=n=lda=Nsize;

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

  /* copy bufM to invM */
  /* For Pfaffian calculation, invM is used as second buffer */
#pragma loop noalias
  for(msi=0;msi<nsize*nsize;msi++) {
    invM[msi] = bufM[msi];
  }
  /* calculate Pf M */
  // [R-Xu] Coverage-0 code. Skipping inverse only.
  // NOTE: One may want to modify this according to calculateMAll_fcmp_child if he needs BF.
  info = 1; // Skip inverse.
  M_ZSKPFA(&uplo, &mthd, &n, invM, &lda, &pfaff, iwork, work, &lwork/*, rwork*/, &info);
  if(info!=0) return info;
  if(!(isfinite(creal(pfaff)) && isfinite(cimag(pfaff)))) return qpidx+1;
  PfM[qpidx] = pfaff;

  /* DInv */
  M_ZGETRF(&m, &n, bufM, &lda, iwork, &info); /* ipiv = iwork */
  if(info!=0) return info;

  M_ZGETRI(&n, bufM, &lda, iwork, work, &lwork, &info);
  if(info!=0) return info;

  /* store InvM */
  /* BufM is column-major, InvM is row-major */
#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    invM_i = invM + msi*Nsize;
    bufM_i = bufM + msi*Nsize;
    bufM_i2 = bufM + msi;
    for(msj=0;msj<nsize;msj++) {
      invM_i[msj] = 0.5*(bufM_i2[msj*nsize] - bufM_i[msj]);
    }
  }

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
    //printf("Debug: bufM=%ld, sltE=%ld\n", bufM, sltE);
    //    printf("Debug: bufM_i=%ld, sltE_i=%ld, rsi=%ld\n", bufM_i, sltE_i, rsi);
#pragma loop norecurrence
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj] + (msj/Ne)*Nsite;
      bufM_i[msj] = -sltE_i[rsj];

    }
  }
  /*
     printf("DEBUG: n=%d \n",n);
     for(msi=0;msi<nsize;msi++){
     for(msj=0;msj<nsize;msj++){
     printf("DEBUG: msi=%d msj=%d bufM %lf \n",msi,msj,bufM[msi+msj*n]);
     }
     }
     */

  /* copy bufM to invM */
  /* For Pfaffian calculation, invM is used as second buffer */
#pragma loop noalias
  for(msi=0;msi<nsize*nsize;msi++) {
    invM[msi] = bufM[msi];
  }

  /* calculate Pf M */
  M_DSKPFA(&uplo, &mthd, &n, invM, &lda, &pfaff, iwork, work, &lwork, &info);
  if(info!=0) return info;
  if(!isfinite(pfaff)) return qpidx+1;
  PfM_real[qpidx] = pfaff;

  //  printf("Debug: M_DGETRF\n");
  /* DInv */
  // M_DGETRF(&m, &n, bufM, &lda, iwork, &info); /* ipiv = iwork */
  // M_DGETRI(&n, bufM, &lda, iwork, work, &lwork, &info);

  // InvM -> InvM' = -InvM
  M_DSCAL(&nsq, &minus_one, invM, &one);

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
  int m,n,lda,info=0;
  double pfaff;

  /* optimization for Kei */
  const int nsize = Nsize;

  const double *sltE = SlaterElmBF_real + (qpidx+qpStart)*Nsite2*Nsite2;
  const double *sltE_i;

  double *invM = InvM_real + qpidx*Nsize*Nsize;
  double *invM_i;

  double *bufM_i, *bufM_i2;

  m=n=lda=Nsize;

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

  /* copy bufM to invM */
  /* For Pfaffian calculation, invM is used as second buffer */
#pragma loop noalias
  for(msi=0;msi<nsize*nsize;msi++) {
    invM[msi] = bufM[msi];
  }
  /* calculate Pf M */
  // [R-Xu]
  // NOTE: This piece of code is coverage 0, i.e. cannot be tested
  //       hence I'm keeping GETRF and GETRI to safety.
  // TODO: If anyone'd like to utilize this, try to use Pfaffine's inv-mode and remove GETRF&GETRI.
  info = 1; // Skip inverse.
  M_DSKPFA(&uplo, &mthd, &n, invM, &lda, &pfaff, iwork, work, &lwork, &info);
  if(info!=0) return info;
  if(!isfinite(pfaff)) return qpidx+1;
  PfM_real[qpidx] = pfaff;

  /* DInv */
  M_DGETRF(&m, &n, bufM, &lda, iwork, &info); /* ipiv = iwork */
  if(info!=0) return info;

  M_DGETRI(&n, bufM, &lda, iwork, work, &lwork, &info);
  if(info!=0) return info;

  /* store InvM */
  /* BufM is column-major, InvM is row-major */
#pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    invM_i = invM + msi*Nsize;
    bufM_i = bufM + msi*Nsize;
    bufM_i2 = bufM + msi;
    for(msj=0;msj<nsize;msj++) {
      invM_i[msj] = 0.5*(bufM_i2[msj*nsize] - bufM_i[msj]);
      /* invM[i][j] = 0.5*(bufM[i][j]-bufM[j][i]) */
    }
  }

  return info;
}

//==============e real =============//

#endif
