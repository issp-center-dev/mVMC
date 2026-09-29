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
#include <complex.h>
#include <stdint.h>
#include <limits.h>
#include <string.h>
#include "./include/global.h"
#include "./include/pfupdate.h"
#include "./include/pfupdate_real.h"


enum { BF_FAIL_PROPOSAL=1, BF_FAIL_ACCEPT_LU, BF_FAIL_ACCEPT_INVERSE,
       BF_FAIL_ACCEPT_PF, BF_FAIL_FULL_PF, BF_FAIL_FULL_INVERSE, BF_FAIL_PERIODIC };
#ifdef MVMC_ENABLE_FAULT_INJECTION
static unsigned BFFailureMask, BFFailureNegativeMask, BFFailureNanMask, BFFailureConsumed;
static int BFFailureForceAccept;
static int BFInjectNumericalInfo(int stage, int globalQp, int info) {
  const unsigned bit=1u<<stage;
  if(globalQp != 0 || !(BFFailureMask & bit) || (BFFailureConsumed & bit)
      || (BFFailureNanMask & bit)) return info;
  BFFailureConsumed |= bit;
  return BFFailureNegativeMask & bit ? -4 : 1;
}
static int BFInjectNumericalNan(int stage, int globalQp) {
  const unsigned bit=1u<<stage;
  if(globalQp != 0 || !(BFFailureNanMask & bit) || (BFFailureConsumed & bit)) return 0;
  BFFailureConsumed |= bit;
  return 1;
}
#else
#define BFInjectNumericalInfo(stage,globalQp,info) (info)
#define BFInjectNumericalNan(stage,globalQp) 0
#endif

/* A normalized status is safe to aggregate with MPI_MAX; raw LAPACK info is
 * signed and must never be used as a collective severity. */
static int BFStatusFromLapack(int info) {
  return info < 0 ? BF_PF_INVALID_ARGUMENT : (info > 0 ? BF_PF_LAPACK_FAILURE : BF_PF_OK);
}
static void BFRequirePfSuccess(int status, const char *stage) {
  if(status != BF_PF_OK) {
    fprintf(stderr,"Error: BackFlow %s failed (status=%d).\n",stage,status);
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
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
  reduction(max:info) private(myIWork,myWork,myRWork, myInfo,myBufM) //TODO: Check to add myRWork is correct or not.
  {
    myIWork = GetWorkSpaceThreadInt(Nsize);
    myBufM = GetWorkSpaceThreadComplex(Nsize*Nsize);
    myWork = GetWorkSpaceThreadComplex(LapackLWork);

    myRWork = GetWorkSpaceThreadDouble(LapackLWork); //TBC for rwork

#pragma omp for private(qpidx)
    for(qpidx=0;qpidx<qpNum;qpidx++) {

      myInfo = calculateMAll_BF_fcmp_child(eleIdx, qpStart, qpEnd, qpidx,
          myBufM, myIWork, myWork, LapackLWork, myRWork, PfM, InvM);
      if(myInfo!=0) {
        if(myInfo > info) info=myInfo;
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

  /* Pfaffian/inverse computed separately. */
  /* Copy bufM to invM before using bufM to compute Pfaffian. */
  for(msi=0;msi<nsize*nsize;msi++)
    invM[msi] = bufM[msi];

  /* Calculate Pf M */
  M_ZSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff, iwork, work, &lwork, rwork, &info);

  if(info!=0) return BFStatusFromLapack(info);
  if(!(isfinite(creal(pfaff)) && isfinite(cimag(pfaff)))) return BF_PF_NONFINITE;
  PfM[qpidx] = pfaff;

  /* Calculate inverse. */
  M_ZGETRF(&m, &n, invM, &lda, iwork, &info); /* ipiv = iwork */
  if(info!=0) return BFStatusFromLapack(info);
  M_ZGETRI(&n, invM, &lda, iwork, work, &lwork, &info);
  if(info!=0) return BFStatusFromLapack(info);

  /* InvM -> InvM(T) -> -InvM */
  M_ZSCAL(&nsq, &minus_one, invM, &one);

  for(msi=0;msi<nsize*nsize;msi++)
    if(!(isfinite(creal(invM[msi])) && isfinite(cimag(invM[msi])))) return BF_PF_NONFINITE;
  return BF_PF_OK;
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

  /* Keep real BackFlow MultiQP rebuilds serialized; Linux CI shows the
     concurrent QP rebuild path can corrupt the identity-energy run. */
#pragma omp parallel if(qpNum == 1) default(shared)              \
  reduction(max:info) private(myIWork,myWork,myInfo,myBufM)
  {
    myIWork = GetWorkSpaceThreadInt(Nsize);
    myBufM = GetWorkSpaceThreadDouble(Nsize*Nsize);
    myWork = GetWorkSpaceThreadDouble(LapackLWork);

#pragma omp for private(qpidx)
    for(qpidx=0;qpidx<qpNum;qpidx++) {

      myInfo = calculateMAll_BF_real_child(eleIdx, qpStart, qpEnd, qpidx,
          myBufM, myIWork, myWork, LapackLWork, PfM_real, InvM_real);
      if(myInfo!=0) {
        if(myInfo > info) info=myInfo;
      }
    }
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadDouble();
  return info;
}

int calculateMAll_BF_real_child(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
    double *bufM, int *iwork, double *work, int lwork, double* PfM_real, double* InvM_real) {
  const size_t invMQpStride = (size_t)Nsize*(size_t)Nsize;

#pragma procedure serial
  return CalculateMAll_BF_real_from_workspace(
      SlaterElmBF_real, eleIdx, qpStart+qpidx, qpStart+qpidx+1,
      PfM_real+qpidx, InvM_real+(size_t)qpidx*invMQpStride,
      invMQpStride, bufM, iwork, work, lwork);
}

int CalculateMAll_BF_real_from_workspace(const double *sltElmBF,
    const int *eleIdx, const int qpStart, const int qpEnd,
    double *pfMOut, double *invMOut, const size_t invMQpStride,
    double *bufM, int *iwork, double *work, int lwork) {
#pragma procedure serial
  const size_t nsizeSquared = (size_t)Nsize*(size_t)Nsize;
  int qpidx;

  if(qpStart < 0 || qpEnd < qpStart || qpEnd > NQPFull) return BF_PF_INVALID_ARGUMENT;
  if(qpStart == qpEnd) return 0;
  if(sltElmBF == NULL || eleIdx == NULL || pfMOut == NULL
     || invMOut == NULL || bufM == NULL || iwork == NULL || work == NULL
     || Nsize <= 0 || Nsite <= 0 || Nsite2 <= 0 || Ne <= 0
     || lwork < LapackLWork || invMQpStride < nsizeSquared
     || invMQpStride > (size_t)PTRDIFF_MAX
     || (qpEnd-qpStart > 1
         && invMQpStride > SIZE_MAX/(size_t)(qpEnd-qpStart-1))) {
    return BF_PF_INVALID_ARGUMENT;
  }
  for(qpidx=0;qpidx<Nsize;qpidx++) {
    if(eleIdx[qpidx] < 0 || eleIdx[qpidx] >= Nsite) return BF_PF_INVALID_ARGUMENT;
  }

  for(qpidx=0;qpidx<qpEnd-qpStart;qpidx++) {
    int msi,msj;
    int rsi,rsj;

    char uplo='U', mthd='P';
    int m,n,nsq,one,lda,info=0;
    double pfaff,minus_one;

    const int nsize = Nsize;
    const int globalQpidx = qpStart+qpidx;

    const double *sltE = sltElmBF
        + (size_t)globalQpidx*(size_t)Nsite2*(size_t)Nsite2;
    const double *sltE_i;

    double *invM = invMOut+(size_t)qpidx*invMQpStride;

    double *bufM_i;

    m=n=lda=Nsize;
    nsq=n*n;
    one=1;
    minus_one=-1.0;

    /* bufM is column-major and skew-symmetric. */
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

    for(msi=0;msi<nsize*nsize;msi++) invM[msi] = bufM[msi];

    M_DSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff,
             iwork, work, &lwork, &info);

    if(info!=0) return BFStatusFromLapack(info);
    if(!isfinite(pfaff)) return BF_PF_NONFINITE;
    pfMOut[qpidx] = pfaff;

    M_DGETRF(&m, &n, invM, &lda, iwork, &info);
    if(info!=0) return BFStatusFromLapack(info);
    M_DGETRI(&n, invM, &lda, iwork, work, &lwork, &info);
    if(info!=0) return BFStatusFromLapack(info);

    M_DSCAL(&nsq, &minus_one, invM, &one);
    for(msi=0;msi<nsize*nsize;msi++) {
      if(!isfinite(invM[msi])) return BF_PF_NONFINITE;
    }
  }

  return 0;
}

//==============e real =============//
double complex updateMAll_BF_fcmp_child(
        const int qpidx,
        const int globalQpidx,
        const int n, const int *msa,
        const int *eleIdx, double complex *candidateInv, int *status);


static int CalculateNewPfMBFWithStrideChecked(const int *icount, const int *msaTmp, const int msaStride,
                       double complex* pfMNew, const int *eleIdx,
                       const int qpStart, const int qpEnd, const double complex* bufM) {
  if(qpStart < 0 || qpEnd < qpStart || qpEnd > NQPFull || !icount || !msaTmp || !pfMNew || !eleIdx)
    return BF_PF_INVALID_ARGUMENT;
  if(msaStride < Nsize || !bufM) return BF_PF_INVALID_ARGUMENT;

  //#pragma procedure serial
  int i, childStatus;
  const int qpNum = qpEnd-qpStart;
  int qpidx;
  int globalQpidx;
  int *msa;

  for(qpidx=0;qpidx<qpNum;qpidx++) {
    globalQpidx = qpidx + qpStart;
    if(icount[globalQpidx] < 0 || icount[globalQpidx] > Nsize) return BF_PF_INVALID_ARGUMENT;
    //Store msa//
    msa=(int *)malloc(sizeof(int)*icount[globalQpidx]);
    //printf("Total=%d\n",icount[qpidx]);
    if(!msa && icount[globalQpidx] != 0) return BF_PF_INVALID_ARGUMENT;
    for(i=0;i<icount[globalQpidx];i++){
      msa[i] = msaTmp[i+globalQpidx*msaStride];
      //printf("hop[%d]=%d\n",i,msa[i]);
    }

    /* calculateNewPfM */
    pfMNew[qpidx] = calculateNewPfMBFN4_child(qpidx,globalQpidx,icount[globalQpidx],msa,eleIdx,bufM, &childStatus);

    free(msa);
    if(childStatus != BF_PF_OK) {
      fprintf(stderr,"BackFlow kernel failure: qp=%d status=%d\n",globalQpidx,childStatus);
      return childStatus;
    }
  }

  return BF_PF_OK;
}


void CalculateNewPfMBFWithStride(const int *icount, const int *msaTmp, const int msaStride,
    double complex *pfMNew, const int *eleIdx, int qpStart, int qpEnd, const double complex *bufM) {
  BFRequirePfSuccess(CalculateNewPfMBFWithStrideChecked(icount,msaTmp,msaStride,pfMNew,eleIdx,qpStart,qpEnd,bufM),"Green proposal");
}
int CalculateNewPfMBFChecked(const int *icount, const int *msaTmp,
    double complex *pfMNew, const int *eleIdx, int qpStart, int qpEnd, const double complex *bufM) {
  return CalculateNewPfMBFWithStrideChecked(icount,msaTmp,Nsize,pfMNew,eleIdx,qpStart,qpEnd,bufM);
}

void CalculateNewPfMBF(const int *icount, const int *msaTmp,
                       double complex* pfMNew, const int *eleIdx,
                       const int qpStart, const int qpEnd, const double complex* bufM) {
  CalculateNewPfMBFWithStride(icount, msaTmp, Nsize, pfMNew, eleIdx, qpStart, qpEnd, bufM);
}

double complex calculateNewPfMBFN4_child(const int qpidx, const int globalQpidx, const int n, const int *msa,
                                         const int *eleIdx, const double complex* bufM, int *status)
{
  if(n < 0 || n > Nsize || 4LL*n*n > INT_MAX) { *status=BF_PF_INVALID_ARGUMENT; return 0.0; }
  for(int row=0;row<n;row++) if(msa[row] < 0 || msa[row] >= Nsize) {
    *status=BF_PF_INVALID_ARGUMENT;return 0.0;
  }

  *status = BF_PF_OK;
  if(n == 0) return PfM[qpidx];
  const int nsize = Nsize;
  const int n2 = 2*n;
  const double complex *sltE;
  const double complex *sltE_k;
  double complex *invM;
  double complex *invM_i, *invM_k, *invM_l;

  double complex *vec; /* vec[n][nsize] */
  //double complex vec[n*nsize]; /* vec[n][nsize] */
  double complex *vec_k, *vec_l;
  //double complex mat[n2*n2]; /* mat[n2][n2] */
  double complex *mat; /* mat[n2][n2] */
  double complex *mat_k;
  //double complex matUV[n2*nsize]; /* mat[n2][nsize] */
  double sgn;

  int rsi,rsk,msi,msj,k,l;
  double complex val,tmp;

  /* for ZSKPFA */
  char uplo='U', mthd='P';
  int nn,lda,info=0;
  double complex pfaff;
  int iwork[n2];
  //double complex work[n2*n2]; /* [n2][n2] */
  double complex *work; /* [n2][n2] */
  double rwork[n2*n2]; /* [n2][n2] */
  int lwork = n2*n2;
  nn=lda=n2;

  sltE = bufM + globalQpidx*Nsite2*Nsite2;
  invM = InvM + qpidx*Nsize*Nsize;

  vec = (double complex*)malloc(sizeof(double complex)*n*nsize);
  //invMat = (double *)malloc(sizeof(double)*n2*n2);
  mat = (double complex*)malloc(sizeof(double complex)*n2*n2);
  //matUV = (double *)malloc(sizeof(double)*n2*nsize);
  work = (double complex*)malloc(sizeof(double complex)*n2*n2);
  if(!vec || !mat || !work) {
    free(vec);free(mat);free(work);
    *status=BF_PF_INVALID_ARGUMENT;return 0.0;
  }

  //#pragma loop noalias
  for(k=0;k<n;k++) {
    rsk = eleIdx[msa[k]] + (msa[k]/Ne)*Nsite;
    sltE_k = sltE + rsk*Nsite2;
    vec_k = vec + k*nsize;
    //#pragma loop norecurrence
    for(msi=0;msi<nsize;msi++) {
      rsi = eleIdx[msi] + (msi/Ne)*Nsite;
      vec_k[msi] = sltE_k[rsi];
    }
  }

  /* X_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k;
    vec_k = vec + k*nsize;
    for(l=k+1;l<n;l++) {
      vec_l = vec + l*nsize;
      val = 0.0;
      for(msi=0;msi<nsize;msi++) {
        invM_i = invM + msi*nsize;
        tmp = 0.0;
        for(msj=0;msj<nsize;msj++) {
          tmp += invM_i[msj] * vec_l[msj];
        }
        val += tmp * vec_k[msi];
      }
      mat_k[l] = val + vec_k[msa[l]];
    }
  }

  /* Y_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k + n;
    vec_k = vec + k*nsize;
    for(l=0;l<n;l++) {
      invM_l = invM + msa[l]*nsize;
      val = 0.0;
      for(msi=0;msi<nsize;msi++) {
        val += vec_k[msi] * invM_l[msi];
      }
      mat_k[l] = val;
    }
  }

  /* Z_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*(k+n) + n;
    invM_k = invM + msa[k]*nsize;
    for(l=k+1;l<n;l++) {
      mat_k[l] = invM_k[msa[l]];
    }
  }

  //#pragma loop noalias
  for(k=0;k<n2;k++) {
    //#pragma loop norecurrence
    for(l=0;l<k;l++) {
      mat[n2*k + l] = -mat[n2*l + k]; /* transpose */
    }
    mat[n2*k + k] = 0.0; /* diagonal elements */
  }

  /* Calculate Pf M */
  //TODO: CHECK rwork is real <-- [R-Xu] Not needed anymore?
  M_ZSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, rwork, &info);
  info = BFInjectNumericalInfo(BF_FAIL_PROPOSAL,globalQpidx,info);
  if(BFInjectNumericalNan(BF_FAIL_PROPOSAL,globalQpidx)) pfaff=NAN;
  *status = BFStatusFromLapack(info);

  //M_DSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, &info);

  sgn = ( (n*(n-1)/2)%2==0 ) ? 1.0 : -1.0;

  //free(matUV);
  free(mat);
  free(vec);
  //free(invMat);
  free(work);
  if(*status != BF_PF_OK) return 0.0;
  const double complex value = sgn * pfaff * PfM[qpidx];
  if(!(isfinite(creal(value)) && isfinite(cimag(value)))) {
    *status = BF_PF_NONFINITE;
    return 0.0;
  }
  return value;
}

int UpdateMAll_BF_fcmp(const int *icount, const int *msaTmp,
                        double complex* pfMNew, const int *eleIdx,
                        const int qpStart, const int qpEnd, double complex *candidateInv)
{
  if(qpStart < 0 || qpEnd < qpStart || qpEnd > NQPFull || !icount || !msaTmp || !pfMNew || !eleIdx)
    return BF_PF_INVALID_ARGUMENT;
  if(!candidateInv) return BF_PF_INVALID_ARGUMENT;

#pragma procedure serial
  const int qpNum = qpEnd-qpStart;
  int qpidx;
  int globalQpidx;
  //double complex *sltE;
  //double complex *sltE_i;
  int *msa;
  int i, childStatus;
  //int *hop;
  //double complex diff;

  memcpy(candidateInv, InvM, (size_t)qpNum*Nsize*Nsize*sizeof(double complex));
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    globalQpidx = qpidx + qpStart;
    if(icount[globalQpidx] < 0 || icount[globalQpidx] > Nsize) return BF_PF_INVALID_ARGUMENT;
    //Store msa//
    msa=(int *)malloc(sizeof(int)*icount[globalQpidx]);
    if(!msa && icount[globalQpidx] != 0) return BF_PF_INVALID_ARGUMENT;
    for(i=0;i<icount[globalQpidx];i++){
      msa[i] = msaTmp[i+globalQpidx*Nsize];
    }

    /* calculateNewPfM */
    pfMNew[qpidx] = updateMAll_BF_fcmp_child(qpidx,globalQpidx,icount[globalQpidx],msa,eleIdx, candidateInv, &childStatus);

    free(msa);
    if(childStatus != BF_PF_OK) {
      fprintf(stderr,"BackFlow kernel failure: qp=%d status=%d\n",globalQpidx,childStatus);
      return childStatus;
    }
  }

  return BF_PF_OK;
}




/* msa[k]-th electron hops from rsa[k] to eleIdx[msa[k]] */
/* buffer size = n*Nsize */
//void updateMAllBF3_child(const int qpidx, const int n, const int *msa,
double complex updateMAll_BF_fcmp_child(
        const int qpidx,
        const int globalQpidx,
        const int n, const int *msa,
        const int *eleIdx, double complex *candidateInv, int *status)
{
  if(n < 0 || n > Nsize || 4LL*n*n > INT_MAX) { *status=BF_PF_INVALID_ARGUMENT; return 0.0; }
  for(int row=0;row<n;row++) if(msa[row] < 0 || msa[row] >= Nsize) {
    *status=BF_PF_INVALID_ARGUMENT;return 0.0;
  }

  *status = BF_PF_OK;
  if(n == 0) return PfM[qpidx];
  const int nsize = Nsize;
  const int n2 = 2*n;
  const double complex*sltE;
  const double complex*sltE_k;
  double complex*invM;
  double complex *invM_i, *invM_k, *invM_l;

  double complex*vec; /* vec[n][nsize] */
  double complex*vec_k, *vec_l;
  //double complex mat[n2*n2]; /* mat[n2][n2] */
  double complex*mat; /* mat[n2][n2] */
  double complex*mat_k;
  //double complex invMat[n2*n2]; /* mat[n2][n2] */
  double complex*invMat; /* mat[n2][n2] */
  //double complex matUV[n2*nsize]; /* mat[n2][nsize] */
  double complex*matUV; /* mat[n2][nsize] */
  double complex*invMat_k, *matUV_i, *matUV_j;
  double complex sgn;

  int rsi,rsk,msi,msj,k,l;
  double complex val,tmp;

  /* for DSKPFA */
  char uplo='U', mthd='P';
  int m,nn,lda,info=0;
  double complex pfaff;
  int iwork[n2];
  int iwork2[n2];
  //double complex work[n2*n2]; /* [n2][n2] */
  double rwork[n2*n2]; /* [n2][n2] */
  double complex*work; /* [n2][n2] */
  double complex*work2; /* [n2][n2] */
  int lwork = n2*n2;
  int lwork2 = n2;
  m=nn=lda=n2;

  sltE = SlaterElmBF + globalQpidx*Nsite2*Nsite2;
  invM = candidateInv + qpidx*Nsize*Nsize;

  //vec = bufferc; /* n*nsize */
  vec = (double complex*)malloc(sizeof(double complex)*n*nsize);
  invMat = (double complex*)malloc(sizeof(double complex)*n2*n2);
  mat = (double complex*)malloc(sizeof(double complex)*n2*n2);
  matUV = (double complex*)malloc(sizeof(double complex)*n2*nsize);
  work = (double complex*)malloc(sizeof(double complex)*n2*n2);
  work2 = (double complex*)malloc(sizeof(double complex)*n2);
  if(!vec || !invMat || !mat || !matUV || !work || !work2) {
    free(vec);free(invMat);free(mat);free(matUV);free(work);free(work2);
    *status=BF_PF_INVALID_ARGUMENT;return 0.0;
  }

#pragma loop noalias
  for(k=0;k<n;k++) {
    rsk = eleIdx[msa[k]] + (msa[k]/Ne)*Nsite;
    sltE_k = sltE + rsk*Nsite2;
    vec_k = vec + k*nsize;
#pragma loop norecurrence
    for(msi=0;msi<nsize;msi++) {
      rsi = eleIdx[msi] + (msi/Ne)*Nsite;
      vec_k[msi] = sltE_k[rsi];
    }
  }

  /* X_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k;
    vec_k = vec + k*nsize;
    for(l=k+1;l<n;l++) {
      vec_l = vec + l*nsize;
      val = 0.0;
      for(msi=0;msi<nsize;msi++) {
        invM_i = invM + msi*nsize;
        tmp = 0.0;
        for(msj=0;msj<nsize;msj++) {
          tmp += invM_i[msj] * vec_l[msj];
        }
        val += tmp * vec_k[msi];
      }
      mat_k[l] = val + vec_k[msa[l]];
    }
  }

  /* Y_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k + n;
    vec_k = vec + k*nsize;
    for(l=0;l<n;l++) {
      invM_l = invM + msa[l]*nsize;
      val = 0.0;
      for(msi=0;msi<nsize;msi++) {
        val += vec_k[msi] * invM_l[msi];
      }
      mat_k[l] = val;
    }
  }

  /* Z_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*(k+n) + n;
    invM_k = invM + msa[k]*nsize;
    for(l=k+1;l<n;l++) {
      mat_k[l] = invM_k[msa[l]];
    }
  }

#pragma loop noalias
  for(k=0;k<n2;k++) {
#pragma loop norecurrence
    for(l=0;l<k;l++) {
      mat[n2*k + l] = -mat[n2*l + k]; /* transpose */
    }
    mat[n2*k + k] = 0.0; /* diagonal elements */
  }

  //#pragma loop noalias
  for(k=0;k<n2;k++) {
    //#pragma loop norecurrence
    for(l=0;l<n2;l++) {
      invMat[n2*k + l] = mat[n2*l + k];
    }
  }

  /* matUV_ik */
  for(msi=0;msi<nsize;msi++) {
    matUV_i = matUV + n2*msi;
    invM_i = invM + msi*nsize;
    for(k=0;k<n;k++) {
      vec_k = vec + k*nsize;
      val = 0.0;
      for(msj=0;msj<nsize;msj++) {
        val -= invM_i[msj] * vec_k[msj];
      }
      if(msa[k]==msi){val -= 1.0;}
      matUV_i[k] = val;
      matUV_i[k+n] = invM_i[msa[k]];
    }
  }
  /* DInv */
  m=nn=lda=n2;
  M_ZGETRF(&m, &nn, invMat, &lda, iwork2, &info);
  info = BFInjectNumericalInfo(BF_FAIL_ACCEPT_LU,globalQpidx,info); /* ipiv = iwork */
  //M_DGETRF(&m, &nn, invMat, &lda, iwork2, &info); /* ipiv = iwork */
  if(info!=0) { *status = BFStatusFromLapack(info); goto cleanup; }

  M_ZGETRI(&nn, invMat, &lda, iwork2, work2, &lwork2, &info);
  info = BFInjectNumericalInfo(BF_FAIL_ACCEPT_INVERSE,globalQpidx,info);
  //M_DGETRI(&nn, invMat, &lda, iwork2, work2, &lwork2, &info);
  if(info!=0) { *status = BFStatusFromLapack(info); goto cleanup; }

  for(msi=0;msi<nsize;msi++) {
    matUV_i = matUV + n2*msi;
    invM_i = invM + msi*nsize;
    for(msj=0;msj<nsize;msj++) {
      matUV_j = matUV + n2*msj;
      val = 0.0;
      for(k=0;k<n2;k++) {
        invMat_k = invMat + k*n2;
        tmp=0.0;
        for(l=0;l<n2;l++) {
          tmp += invMat_k[l]*matUV_j[l];
        }
        val += matUV_i[k]*tmp;
      }
      invM_i[msj] -= val;
    }
  }

  for(msi=0;msi<Ne;msi++) {
    //invM_i = invM + msi*nsize;
    for(msj=Ne;msj<Nsize;msj++) {
      invM[msj*nsize+msi] = -invM[msi*nsize+msj];
    }
  }
  for(msi=Ne;msi<Nsize;msi++) {
    for(msj=Ne;msj<Nsize;msj++) {
      invM[msi*nsize+msj] = 0.0;
    }
  }
  for(msi=0;msi<Ne;msi++) {
    for(msj=0;msj<Ne;msj++) {
      invM[msi*nsize+msj] = 0.0;
    }
  }

  /* Calculate Pf M */
  M_ZSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, rwork, &info);
  info = BFInjectNumericalInfo(BF_FAIL_ACCEPT_PF,globalQpidx,info);
  *status = BFStatusFromLapack(info);
  if(*status != BF_PF_OK) goto cleanup;
  //M_DSKPFA(&uplo, &mthd, &nn, invM, &lda, &pfaff, iwork, work, &lwork, &info);
  //M_DSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, &info);

  sgn = ( (n*(n-1)/2)%2==0 ) ? 1.0 : -1.0;

cleanup:
  free(matUV);
  free(mat);
  free(vec);
  free(invMat);
  free(work);
  free(work2);

  if(*status != BF_PF_OK) return 0.0;
  const double complex value = sgn * pfaff * PfM[qpidx];
  if(!(isfinite(creal(value)) && isfinite(cimag(value)))) *status = BF_PF_NONFINITE;
  for(msi=0;msi<nsize*nsize;msi++)
    if(!(isfinite(creal(invM[msi])) && isfinite(cimag(invM[msi])))) *status = BF_PF_NONFINITE;
  return *status == BF_PF_OK ? value : 0.0;
}

double calculateNewPfMBFN4_real_child(const int qpidx, const int globalQpidx, const int n, const int *msa,
                                      const int *eleIdx, const double *bufM, int *status);
double calculateNewPfMBFN4_real_child_vec(const int qpidx, const int globalQpidx, const int n, const int *msa, const double *vec, int *status);
static double calculateNewPfMBFN4_real_child_vec_w(const int qpidx, const int globalQpidx, const int n, const int *msa,
                                                   const double *vec, const double *w, int *status);

double updateMAll_BF_real_child(const int qpidx, const int globalQpidx, const int n, const int *msa,
                                const int *eleIdx, double *candidateInv, int *status);

static void DgemmRowMajorNN(const int m, const int n, const int k,
                            const double *a, const double *b,
                            const double alpha, const double beta, double *c) {
  char trans = 'N';
  int fm = n, fn = m, fk = k;
  int lda = n, ldb = k, ldc = n;
  M_DGEMM(&trans, &trans, &fm, &fn, &fk, &alpha, b, &lda, a, &ldb, &beta, c, &ldc);
}

static void DgemmRowMajorNT(const int m, const int n, const int k,
                            const double *a, const double *b,
                            const double alpha, const double beta, double *c) {
  char transa = 'T', transb = 'N';
  int fm = n, fn = m, fk = k;
  int lda = k, ldb = k, ldc = n;
  M_DGEMM(&transa, &transb, &fm, &fn, &fk, &alpha, b, &lda, a, &ldb, &beta, c, &ldc);
}

static int CalculateNewPfMBFWithStride_realChecked(const int *icount, const int *msaTmp, const int msaStride,
                       double *pfMNew, const int *eleIdx,
                       const int qpStart, const int qpEnd, const double *bufM) {
  if(qpStart < 0 || qpEnd < qpStart || qpEnd > NQPFull || !icount || !msaTmp || !pfMNew || !eleIdx)
    return BF_PF_INVALID_ARGUMENT;
  if(msaStride < Nsize || !bufM) return BF_PF_INVALID_ARGUMENT;

  //#pragma procedure serial
  int i, childStatus;
  const int qpNum = qpEnd-qpStart;
  int qpidx;
  int globalQpidx;
  //double *sltE;
  //double *sltE_i;
  int *msa;
  //int msi,msj,rsi,rsj,i;
  //int *hop;
  //double complex diff;

  for(qpidx=0;qpidx<qpNum;qpidx++) {
    globalQpidx = qpidx + qpStart;
    if(icount[globalQpidx] < 0 || icount[globalQpidx] > Nsize) return BF_PF_INVALID_ARGUMENT;
    //Store msa//
    msa=(int *)malloc(sizeof(int)*icount[globalQpidx]);
    //printf("Total=%d\n",icount[qpidx]);
    if(!msa && icount[globalQpidx] != 0) return BF_PF_INVALID_ARGUMENT;
    for(i=0;i<icount[globalQpidx];i++){
      msa[i] = msaTmp[i+globalQpidx*msaStride];
      //printf("hop[%d]=%d\n",i,msa[i]);
    }

    /* calculateNewPfM */
    pfMNew[qpidx] = calculateNewPfMBFN4_real_child(qpidx,globalQpidx,icount[globalQpidx],msa,eleIdx,bufM, &childStatus);

    free(msa);
    if(childStatus != BF_PF_OK) {
      fprintf(stderr,"BackFlow kernel failure: qp=%d status=%d\n",globalQpidx,childStatus);
      return childStatus;
    }
  }
  //icount = UpdateSlaterElmBFTmp3(ma, rb, ra, s, eleCfg, eleNum, msaTmp);

  return BF_PF_OK;
}


void CalculateNewPfMBFWithStride_real(const int *icount, const int *msaTmp, const int msaStride,
    double *pfMNew, const int *eleIdx, int qpStart, int qpEnd, const double *bufM) {
  BFRequirePfSuccess(CalculateNewPfMBFWithStride_realChecked(icount,msaTmp,msaStride,pfMNew,eleIdx,qpStart,qpEnd,bufM),"Green proposal");
}
int CalculateNewPfMBF_realChecked(const int *icount, const int *msaTmp,
    double *pfMNew, const int *eleIdx, int qpStart, int qpEnd, const double *bufM) {
  return CalculateNewPfMBFWithStride_realChecked(icount,msaTmp,Nsize,pfMNew,eleIdx,qpStart,qpEnd,bufM);
}

void CalculateNewPfMBF_real(const int *icount, const int *msaTmp,
                       double *pfMNew, const int *eleIdx,
                       const int qpStart, const int qpEnd, const double *bufM) {
  CalculateNewPfMBFWithStride_real(icount, msaTmp, Nsize, pfMNew, eleIdx, qpStart, qpEnd, bufM);
}

void CalculateNewPfMBFVecWithStride_real(const int *icount, const int *msaTmp, const int msaStride,
                       double *pfMNew, const int qpStart, const int qpEnd,
                       const double *vecM, const int vecStride) {
  int childStatus = BF_PF_OK;
  int i;
  const int qpNum = qpEnd-qpStart;
  int qpidx;
  int globalQpidx;
  int *msa;
  const double *vec;

  for(qpidx=0;qpidx<qpNum;qpidx++) {
    globalQpidx = qpidx + qpStart;
    msa=(int *)malloc(sizeof(int)*icount[globalQpidx]);
    for(i=0;i<icount[globalQpidx];i++){
      msa[i] = msaTmp[i+globalQpidx*msaStride];
    }

    vec = vecM + globalQpidx*vecStride*Nsize;
    pfMNew[qpidx] = calculateNewPfMBFN4_real_child_vec(qpidx,globalQpidx,icount[globalQpidx],msa,vec, &childStatus);
    BFRequirePfSuccess(childStatus,"Green vector proposal");

    free(msa);
  }

  return;
}

void CalculateNewPfMBFVec_real(const int *icount, const int *msaTmp,
                       double *pfMNew, const int qpStart, const int qpEnd,
                       const double *vecM) {
  CalculateNewPfMBFVecWithStride_real(icount, msaTmp, Nsize, pfMNew, qpStart, qpEnd, vecM, Nsize);
}

void CalculateNewPfMBFVecBatched_real(const int batchSize, const int *icount, const int *msaTmp,
                       double *pfMNew, const int qpStart, const int qpEnd,
                       const double *vecM, double *vecStack, double *wStack) {
  int childStatus = BF_PF_OK;
  const int nsize = Nsize;
  const int qpNum = qpEnd-qpStart;
  int qpidx, globalQpidx, batchIdx, row, totalRows, offset;
  int rowOffset[batchSize];
  int n;
  const int *msa;
  const double *termVec;
  const double *invM;

  for(qpidx=0;qpidx<qpNum;qpidx++) {
    globalQpidx = qpidx + qpStart;
    totalRows = 0;
    for(batchIdx=0;batchIdx<batchSize;batchIdx++) {
      n = icount[batchIdx*NQPFull + globalQpidx];
      rowOffset[batchIdx] = totalRows;
      termVec = vecM + ((size_t)batchIdx*NQPFull + globalQpidx)*nsize*nsize;
      for(row=0;row<n;row++) {
        memcpy(vecStack + (totalRows+row)*nsize, termVec + row*nsize, sizeof(double)*nsize);
      }
      totalRows += n;
    }

    if(totalRows > 0) {
      invM = InvM_real + qpidx*nsize*nsize;
      DgemmRowMajorNT(totalRows, nsize, nsize, vecStack, invM, 1.0, 0.0, wStack);
    }

    for(batchIdx=0;batchIdx<batchSize;batchIdx++) {
      n = icount[batchIdx*NQPFull + globalQpidx];
      offset = rowOffset[batchIdx];
      msa = msaTmp + ((size_t)batchIdx*NQPFull + globalQpidx)*nsize;
      pfMNew[batchIdx*NQPFull + globalQpidx] =
        calculateNewPfMBFN4_real_child_vec_w(qpidx,globalQpidx, n, msa,
                                             vecStack + offset*nsize,
                                             wStack + offset*nsize, &childStatus);
    BFRequirePfSuccess(childStatus,"Green vector proposal");
    }
  }

  return;
}

/* msa[k]-th electron hops from rsa[k] to eleIdx[msa[k]] */
/* buffer size = n*Nsize */
//double complex calculateNewPfMBFN_child(const int qpidx, const int n, const int *msa,
//                              const int *eleIdx, double *rwork, double complex *bufferc) {
double calculateNewPfMBFN4_real_child(const int qpidx, const int globalQpidx, const int n, const int *msa,
                                      const int *eleIdx, const double *bufM, int *status) {
  if(n < 0 || n > Nsize || 4LL*n*n > INT_MAX) { *status=BF_PF_INVALID_ARGUMENT; return 0.0; }
  for(int row=0;row<n;row++) if(msa[row] < 0 || msa[row] >= Nsize) {
    *status=BF_PF_INVALID_ARGUMENT;return 0.0;
  }

  const int nsize = Nsize;
  const double *sltE;
  const double *sltE_k;

  double *vec; /* vec[n][nsize] */
  double *vec_k;

  int rsi,rsk,msi,k;
  double pfnew;

  sltE = bufM + globalQpidx*Nsite2*Nsite2;

  vec = (double *)malloc(sizeof(double)*n*nsize);
  if(!vec) {
    free(vec);
    *status=BF_PF_INVALID_ARGUMENT;return 0.0;
  }

  //#pragma loop noalias
  for(k=0;k<n;k++) {
    rsk = eleIdx[msa[k]] + (msa[k]/Ne)*Nsite;
    sltE_k = sltE + rsk*Nsite2;
    vec_k = vec + k*nsize;
    //#pragma loop norecurrence
    for(msi=0;msi<nsize;msi++) {
      rsi = eleIdx[msi] + (msi/Ne)*Nsite;
      vec_k[msi] = sltE_k[rsi];
    }
  }

  pfnew = calculateNewPfMBFN4_real_child_vec(qpidx,globalQpidx,n,msa,vec, status);

  free(vec);
  return pfnew;
}

double calculateNewPfMBFN4_real_child_vec(const int qpidx, const int globalQpidx, const int n, const int *msa, const double *vec, int *status) {
  if(n < 0 || n > Nsize || 4LL*n*n > INT_MAX) { *status=BF_PF_INVALID_ARGUMENT; return 0.0; }
  for(int row=0;row<n;row++) if(msa[row] < 0 || msa[row] >= Nsize) {
    *status=BF_PF_INVALID_ARGUMENT;return 0.0;
  }

  const int nsize = Nsize;
  double *invM;
  double *w; /* w[n][nsize] = vec * invM^T */
  double pfnew;

  invM = InvM_real + qpidx*Nsize*Nsize;

  w = (double *)malloc(sizeof(double)*n*nsize);
  if(!w) {
    free(w);
    *status=BF_PF_INVALID_ARGUMENT;return 0.0;
  }

  DgemmRowMajorNT(n, nsize, nsize, vec, invM, 1.0, 0.0, w);

  pfnew = calculateNewPfMBFN4_real_child_vec_w(qpidx,globalQpidx,n,msa,vec,w, status);

  free(w);
  return pfnew;
}

static double calculateNewPfMBFN4_real_child_vec_w(const int qpidx, const int globalQpidx, const int n, const int *msa,
                                                   const double *vec, const double *w, int *status) {
  if(n < 0 || n > Nsize || 4LL*n*n > INT_MAX) { *status=BF_PF_INVALID_ARGUMENT; return 0.0; }
  for(int row=0;row<n;row++) if(msa[row] < 0 || msa[row] >= Nsize) {
    *status=BF_PF_INVALID_ARGUMENT;return 0.0;
  }

  *status = BF_PF_OK;
  if(n == 0) return PfM_real[qpidx];
  const int nsize = Nsize;
  const int n2 = 2*n;
  double *invM;
  double *invM_k;

  const double *vec_k;
  double *mat; /* mat[n2][n2] */
  double *mat_k;
  double sgn;

  int msi,k,l;
  double val;

  /* for ZSKPFA */
  char uplo='U', mthd='P';
  int nn,lda,info=0;
  double pfaff;
  int iwork[n2];
  double *work; /* [n2][n2] */
  int lwork = n2*n2;
  nn=lda=n2;

  invM = InvM_real + qpidx*Nsize*Nsize;

  mat = (double *)malloc(sizeof(double)*n2*n2);
  work = (double *)malloc(sizeof(double)*n2*n2);
  if(!mat || !work) {
    free(mat);free(work);
    *status=BF_PF_INVALID_ARGUMENT;return 0.0;
  }

  /* X_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k;
    vec_k = vec + k*nsize;
    for(l=k+1;l<n;l++) {
      val = 0.0;
      for(msi=0;msi<nsize;msi++) {
        val += w[l*nsize + msi] * vec_k[msi];
      }
      mat_k[l] = val + vec_k[msa[l]];
    }
  }

  /* Y_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k + n;
    for(l=0;l<n;l++) {
      mat_k[l] = w[k*nsize + msa[l]];
    }
  }

  /* Z_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*(k+n) + n;
    invM_k = invM + msa[k]*nsize;
    for(l=k+1;l<n;l++) {
      mat_k[l] = invM_k[msa[l]];
    }
  }

  //#pragma loop noalias
  for(k=0;k<n2;k++) {
    //#pragma loop norecurrence
    for(l=0;l<k;l++) {
      mat[n2*k + l] = -mat[n2*l + k]; /* transpose */
    }
    mat[n2*k + k] = 0.0; /* diagonal elements */
  }

  /* calculate Pf M */
  M_DSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, &info);
  info = BFInjectNumericalInfo(BF_FAIL_PROPOSAL,globalQpidx,info);
  if(BFInjectNumericalNan(BF_FAIL_PROPOSAL,globalQpidx)) pfaff=NAN;
  *status = BFStatusFromLapack(info);


  sgn = ( (n*(n-1)/2)%2==0 ) ? 1.0 : -1.0;

  //free(matUV);
  free(mat);
  //free(invMat);
  free(work);
  if(*status != BF_PF_OK) return 0.0;
  const double value = sgn * pfaff * PfM_real[qpidx];
  if(!(isfinite(creal(value)) && isfinite(cimag(value)))) {
    *status = BF_PF_NONFINITE;
    return 0.0;
  }
  return value;
}

/* Calculate new pfaffian with Backflow effects.
   The ma-th electron with spin s hops from ra to rb */
int UpdateMAll_BF_real(const int *icount, const int *msaTmp,
                  double *pfMNew, const int *eleIdx,
                  const int qpStart, const int qpEnd, double *candidateInv) {
  if(qpStart < 0 || qpEnd < qpStart || qpEnd > NQPFull || !icount || !msaTmp || !pfMNew || !eleIdx)
    return BF_PF_INVALID_ARGUMENT;
  if(!candidateInv) return BF_PF_INVALID_ARGUMENT;

#pragma procedure serial
  const int qpNum = qpEnd-qpStart;
  int qpidx;
  int globalQpidx;
  //double complex *sltE;
  //double complex *sltE_i;
  int *msa;
  int i, childStatus;
  //int *hop;
  //double complex diff;

  memcpy(candidateInv, InvM_real, (size_t)qpNum*Nsize*Nsize*sizeof(double));
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    globalQpidx = qpidx + qpStart;
    if(icount[globalQpidx] < 0 || icount[globalQpidx] > Nsize) return BF_PF_INVALID_ARGUMENT;
    //Store msa//
    msa=(int *)malloc(sizeof(int)*icount[globalQpidx]);
    if(!msa && icount[globalQpidx] != 0) return BF_PF_INVALID_ARGUMENT;
    for(i=0;i<icount[globalQpidx];i++){
      msa[i] = msaTmp[i+globalQpidx*Nsize];
    }

    /* calculateNewPfM */
    pfMNew[qpidx] = updateMAll_BF_real_child(qpidx,globalQpidx,icount[globalQpidx],msa,eleIdx, candidateInv, &childStatus);

    free(msa);
    if(childStatus != BF_PF_OK) {
      fprintf(stderr,"BackFlow kernel failure: qp=%d status=%d\n",globalQpidx,childStatus);
      return childStatus;
    }
  }

  return BF_PF_OK;
}




/* msa[k]-th electron hops from rsa[k] to eleIdx[msa[k]] */
/* buffer size = n*Nsize */
//void updateMAllBF3_child(const int qpidx, const int n, const int *msa,
double updateMAll_BF_real_child(const int qpidx, const int globalQpidx, const int n, const int *msa,
                          const int *eleIdx, double *candidateInv, int *status) {
  if(n < 0 || n > Nsize || 4LL*n*n > INT_MAX) { *status=BF_PF_INVALID_ARGUMENT; return 0.0; }
  for(int row=0;row<n;row++) if(msa[row] < 0 || msa[row] >= Nsize) {
    *status=BF_PF_INVALID_ARGUMENT;return 0.0;
  }

  *status = BF_PF_OK;
  if(n == 0) return PfM_real[qpidx];
  const int nsize = Nsize;
  const int n2 = 2*n;
  const double *sltE;
  const double *sltE_k;
  double *invM;
  double *invM_i, *invM_k;

  double *vec; /* vec[n][nsize] */
  double *vec_k;
  double *w; /* w[n][nsize] = vec * invM^T */
  //double complex mat[n2*n2]; /* mat[n2][n2] */
  double *mat; /* mat[n2][n2] */
  double *mat_k;
  //double complex invMat[n2*n2]; /* mat[n2][n2] */
  double *invMat; /* mat[n2][n2] */
  //double complex matUV[n2*nsize]; /* mat[n2][nsize] */
  double *matUV; /* mat[n2][nsize] */
  double *woodburyTmp; /* [nsize][n2] */
  double *matUV_i;
  double sgn;

  int rsi,rsk,msi,msj,k,l;
  double val;

  /* for DSKPFA */
  char uplo='U', mthd='P';
  int m,nn,lda,info=0;
  double pfaff;
  int iwork[n2];
  int iwork2[n2];
  //double complex work[n2*n2]; /* [n2][n2] */
  double *work; /* [n2][n2] */
  double *work2; /* [n2][n2] */
  int lwork = n2*n2;
  int lwork2 = n2;
  m=nn=lda=n2;

  sltE = SlaterElmBF_real + globalQpidx*Nsite2*Nsite2;
  invM = candidateInv + qpidx*Nsize*Nsize;

  //vec = bufferc; /* n*nsize */
  vec = (double *)malloc(sizeof(double)*n*nsize);
  w = (double *)malloc(sizeof(double)*n*nsize);
  invMat = (double *)malloc(sizeof(double)*n2*n2);
  mat = (double *)malloc(sizeof(double)*n2*n2);
  matUV = (double *)malloc(sizeof(double)*n2*nsize);
  woodburyTmp = (double *)malloc(sizeof(double)*nsize*n2);
  work = (double *)malloc(sizeof(double)*n2*n2);
  work2 = (double *)malloc(sizeof(double)*n2);
  if(!vec || !w || !invMat || !mat || !matUV || !woodburyTmp || !work || !work2) {
    free(vec);free(w);free(invMat);free(mat);free(matUV);free(woodburyTmp);free(work);free(work2);
    *status=BF_PF_INVALID_ARGUMENT;return 0.0;
  }

#pragma loop noalias
  for(k=0;k<n;k++) {
    rsk = eleIdx[msa[k]] + (msa[k]/Ne)*Nsite;
    sltE_k = sltE + rsk*Nsite2;
    vec_k = vec + k*nsize;
#pragma loop norecurrence
    for(msi=0;msi<nsize;msi++) {
      rsi = eleIdx[msi] + (msi/Ne)*Nsite;
      vec_k[msi] = sltE_k[rsi];
    }
  }

  DgemmRowMajorNT(n, nsize, nsize, vec, invM, 1.0, 0.0, w);

  /* X_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k;
    vec_k = vec + k*nsize;
    for(l=k+1;l<n;l++) {
      val = 0.0;
      for(msi=0;msi<nsize;msi++) {
        val += w[l*nsize + msi] * vec_k[msi];
      }
      mat_k[l] = val + vec_k[msa[l]];
    }
  }

  /* Y_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k + n;
    for(l=0;l<n;l++) {
      mat_k[l] = w[k*nsize + msa[l]];
    }
  }

  /* Z_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*(k+n) + n;
    invM_k = invM + msa[k]*nsize;
    for(l=k+1;l<n;l++) {
      mat_k[l] = invM_k[msa[l]];
    }
  }

#pragma loop noalias
  for(k=0;k<n2;k++) {
#pragma loop norecurrence
    for(l=0;l<k;l++) {
      mat[n2*k + l] = -mat[n2*l + k]; /* transpose */
    }
    mat[n2*k + k] = 0.0; /* diagonal elements */
  }

  //#pragma loop noalias
  for(k=0;k<n2;k++) {
    //#pragma loop norecurrence
    for(l=0;l<n2;l++) {
      invMat[n2*k + l] = mat[n2*l + k];
    }
  }

  /* matUV_ik */
  for(msi=0;msi<nsize;msi++) {
    matUV_i = matUV + n2*msi;
    invM_i = invM + msi*nsize;
    for(k=0;k<n;k++) {
      val = -w[k*nsize + msi];
      if(msa[k]==msi){val -= 1.0;}
      matUV_i[k] = val;
      matUV_i[k+n] = invM_i[msa[k]];
    }
  }
  /* DInv */
  // NOTE: This routines is not called in current version (i.e. lcoal coverage 0),
  //       hence cannot be tested. I'm keeping GETRF and GETRI here.
  // TODO: If anyone's interested in utilizing this routine, she or he might want to
  //       utilize DSKTRF+LTL2INV
  m=nn=lda=n2;
  //M_ZGETRF(&m, &nn, invMat, &lda, iwork2, &info); /* ipiv = iwork */
  M_DGETRF(&m, &nn, invMat, &lda, iwork2, &info);
  info = BFInjectNumericalInfo(BF_FAIL_ACCEPT_LU,globalQpidx,info); /* ipiv = iwork */
  if(info!=0) { *status = BFStatusFromLapack(info); goto cleanup; }

  //M_ZGETRI(&nn, invMat, &lda, iwork2, work2, &lwork2, &info);
  M_DGETRI(&nn, invMat, &lda, iwork2, work2, &lwork2, &info);
  info = BFInjectNumericalInfo(BF_FAIL_ACCEPT_INVERSE,globalQpidx,info);
  if(info!=0) { *status = BFStatusFromLapack(info); goto cleanup; }

  DgemmRowMajorNN(nsize, n2, n2, matUV, invMat, 1.0, 0.0, woodburyTmp);
  DgemmRowMajorNT(nsize, nsize, n2, woodburyTmp, matUV, -1.0, 1.0, invM);

  for(msi=0;msi<Ne;msi++) {
    //invM_i = invM + msi*nsize;
    for(msj=Ne;msj<Nsize;msj++) {
      invM[msj*nsize+msi] = -invM[msi*nsize+msj];
    }
  }
  for(msi=Ne;msi<Nsize;msi++) {
    for(msj=Ne;msj<Nsize;msj++) {
      invM[msi*nsize+msj] = 0.0;
    }
  }
  for(msi=0;msi<Ne;msi++) {
    for(msj=0;msj<Ne;msj++) {
      invM[msi*nsize+msj] = 0.0;
    }
  }

  /* calculate Pf M */
  //M_ZSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, rwork, &info);
  //M_DSKPFA(&uplo, &mthd, &nn, invM, &lda, &pfaff, iwork, work, &lwork, &info);
  M_DSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, &info);
  info = BFInjectNumericalInfo(BF_FAIL_ACCEPT_PF,globalQpidx,info);
  *status = BFStatusFromLapack(info);
  if(*status != BF_PF_OK) goto cleanup;

  sgn = ( (n*(n-1)/2)%2==0 ) ? 1.0 : -1.0;

cleanup:
  free(matUV);
  free(mat);
  free(w);
  free(vec);
  free(invMat);
  free(woodburyTmp);
  free(work);
  free(work2);

  if(*status != BF_PF_OK) return 0.0;
  const double value = sgn * pfaff * PfM_real[qpidx];
  if(!(isfinite(creal(value)) && isfinite(cimag(value)))) *status = BF_PF_NONFINITE;
  for(msi=0;msi<nsize*nsize;msi++)
    if(!(isfinite(creal(invM[msi])) && isfinite(cimag(invM[msi])))) *status = BF_PF_NONFINITE;
  return *status == BF_PF_OK ? value : 0.0;
}
