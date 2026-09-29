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

#include "backflow_stable.c"

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
  if(qpidx<0 || qpidx>=qpEnd-qpStart) return BF_PF_INVALID_ARGUMENT;
  return CalculateMAll_BF_fcmp_from_workspace(SlaterElmBF,eleIdx,qpStart+qpidx,qpStart+qpidx+1,
      PfM+qpidx,InvM+(size_t)qpidx*Nsize*Nsize,(size_t)Nsize*Nsize,
      bufM,iwork,work,lwork,rwork);
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
  if(qpStart<0 || qpEnd<qpStart || qpEnd>NQPFull) return BF_PF_INVALID_ARGUMENT;
  if(qpStart==qpEnd) return BF_PF_OK;
  if(Nsize<=0 || !pfMOut || !invMOut || !bufM || !iwork || !work
     || lwork<Nsize || invMQpStride<(size_t)Nsize*Nsize
     || invMQpStride>(size_t)PTRDIFF_MAX/sizeof(*invMOut)/(size_t)(qpEnd-qpStart))
    return BF_PF_INVALID_ARGUMENT;
  for(int q=qpStart;q<qpEnd;q++) {
    int status=BFStableBuild_real(sltElmBF,eleIdx,q,0,NULL,NULL,bufM);
    if(status!=BF_PF_OK) return status;
    status=BFStableFactor_real(bufM,pfMOut+q-qpStart,
        invMOut+(size_t)(q-qpStart)*invMQpStride,iwork,work,lwork,q,0);
    if(status!=BF_PF_OK) return status;
  }
  return BF_PF_OK;
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
  if(n<0 || n>Nsize || (n && !msa)) { *status=BF_PF_INVALID_ARGUMENT; return 0.0; }
  for(int k=0;k<n;k++) if(msa[k]<0 || msa[k]>=Nsize) { *status=BF_PF_INVALID_ARGUMENT; return 0.0; }
  *status=BF_PF_OK;
  if(n==0) return PfM[qpidx];
  double complex pf=0.0;
  *status=BFStableEvaluate_fcmp(bufM,eleIdx,globalQpidx,0,NULL,NULL,&pf,NULL,BF_FAIL_PROPOSAL);
  return *status==BF_PF_OK ? pf : 0.0;
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
  if(n<0 || n>Nsize || (n && !msa)) { *status=BF_PF_INVALID_ARGUMENT; return 0.0; }
  for(int k=0;k<n;k++) if(msa[k]<0 || msa[k]>=Nsize) { *status=BF_PF_INVALID_ARGUMENT; return 0.0; }
  *status=BF_PF_OK;
  if(n==0) return PfM[qpidx];
  double complex pf=0.0;
  *status=BFStableEvaluate_fcmp(SlaterElmBF,eleIdx,globalQpidx,0,NULL,NULL,&pf,
      candidateInv+(size_t)qpidx*Nsize*Nsize,BF_FAIL_ACCEPT_PF);
  return *status==BF_PF_OK ? pf : 0.0;
}

double calculateNewPfMBFN4_real_child(const int qpidx, const int globalQpidx, const int n, const int *msa,
                                      const int *eleIdx, const double *bufM, int *status);
static double calculateNewPfMBFN4_real_child_vec(const int qpidx, const int globalQpidx,
    const int n, const int *msa, const double *vec, const int *eleIdx, int *status) {
  *status=BF_PF_OK;
  if(n==0) return PfM_real[qpidx];
  double pf=0.0;
  *status=BFStableEvaluate_real(SlaterElmBF_real,eleIdx,globalQpidx,n,msa,vec,&pf,NULL,BF_FAIL_PROPOSAL);
  return *status==BF_PF_OK ? pf : 0.0;
}

double updateMAll_BF_real_child(const int qpidx, const int globalQpidx, const int n, const int *msa,
                                const int *eleIdx, double *candidateInv, int *status);

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
                       const double *vecM, const int vecStride, const int *eleIdx) {
  if(qpStart<0 || qpEnd<qpStart || qpEnd>NQPFull || msaStride<Nsize || vecStride<Nsize
     || !icount || !msaTmp || !pfMNew || !vecM || !eleIdx) {
    BFRequirePfSuccess(BF_PF_INVALID_ARGUMENT,"Green vector input"); return;
  }
  for(int q=qpStart;q<qpEnd;q++) {
    int status=BF_PF_OK;
    pfMNew[q-qpStart]=calculateNewPfMBFN4_real_child_vec(q-qpStart,q,icount[q],
        msaTmp+(size_t)q*msaStride,vecM+(size_t)q*vecStride*Nsize,eleIdx,&status);
    BFRequirePfSuccess(status,"Green vector proposal");
  }
}

void CalculateNewPfMBFVec_real(const int *icount, const int *msaTmp,
                       double *pfMNew, const int qpStart, const int qpEnd,
                       const double *vecM, const int *eleIdx) {
  CalculateNewPfMBFVecWithStride_real(icount, msaTmp, Nsize, pfMNew, qpStart, qpEnd, vecM, Nsize, eleIdx);
}

void CalculateNewPfMBFVecBatched_real(const int batchSize, const int *icount, const int *msaTmp,
                       double *pfMNew, const int qpStart, const int qpEnd,
                       const double *vecM, const int *eleIdx) {
  if(batchSize<0 || qpStart<0 || qpEnd<qpStart || qpEnd>NQPFull
     || !icount || !msaTmp || !pfMNew || !vecM || !eleIdx) {
    BFRequirePfSuccess(BF_PF_INVALID_ARGUMENT,"Green batch input"); return;
  }
  for(int b=0;b<batchSize;b++) for(int q=qpStart;q<qpEnd;q++) {
    const size_t offset=(size_t)b*NQPFull+q;
    int status=BF_PF_OK;
    pfMNew[offset]=calculateNewPfMBFN4_real_child_vec(q-qpStart,q,icount[offset],
        msaTmp+offset*Nsize,vecM+offset*Nsize*Nsize,eleIdx,&status);
    BFRequirePfSuccess(status,"Green batch proposal");
  }
}

/* msa[k]-th electron hops from rsa[k] to eleIdx[msa[k]] */
/* buffer size = n*Nsize */
//double complex calculateNewPfMBFN_child(const int qpidx, const int n, const int *msa,
//                              const int *eleIdx, double *rwork, double complex *bufferc) {
double calculateNewPfMBFN4_real_child(const int qpidx, const int globalQpidx, const int n, const int *msa,
                                      const int *eleIdx, const double *bufM, int *status) {
  if(n<0 || n>Nsize || (n && !msa)) { *status=BF_PF_INVALID_ARGUMENT; return 0.0; }
  for(int k=0;k<n;k++) if(msa[k]<0 || msa[k]>=Nsize) { *status=BF_PF_INVALID_ARGUMENT; return 0.0; }
  *status=BF_PF_OK;
  if(n==0) return PfM_real[qpidx];
  double pf=0.0;
  *status=BFStableEvaluate_real(bufM,eleIdx,globalQpidx,0,NULL,NULL,&pf,NULL,BF_FAIL_PROPOSAL);
  return *status==BF_PF_OK ? pf : 0.0;
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
  if(n<0 || n>Nsize || (n && !msa)) { *status=BF_PF_INVALID_ARGUMENT; return 0.0; }
  for(int k=0;k<n;k++) if(msa[k]<0 || msa[k]>=Nsize) { *status=BF_PF_INVALID_ARGUMENT; return 0.0; }
  *status=BF_PF_OK;
  if(n==0) return PfM_real[qpidx];
  double pf=0.0;
  *status=BFStableEvaluate_real(SlaterElmBF_real,eleIdx,globalQpidx,0,NULL,NULL,&pf,
      candidateInv+(size_t)qpidx*Nsize*Nsize,BF_FAIL_ACCEPT_PF);
  return *status==BF_PF_OK ? pf : 0.0;
}

int CalculateMAll_BF_fcmp_from_workspace(const double complex *sltElmBF,
    const int *eleIdx, int qpStart, int qpEnd,
    double complex *pfMOut, double complex *invMOut, size_t invMQpStride,
    double complex *bufM, int *iwork, double complex *work, int lwork, double *rwork) {
  if(qpStart<0 || qpEnd<qpStart || qpEnd>NQPFull) return BF_PF_INVALID_ARGUMENT;
  if(qpStart==qpEnd) return BF_PF_OK;
  if(Nsize<=0 || !pfMOut || !invMOut || !bufM || !iwork || !work || !rwork
     || lwork<Nsize || invMQpStride<(size_t)Nsize*Nsize
     || invMQpStride>(size_t)PTRDIFF_MAX/sizeof(*invMOut)/(size_t)(qpEnd-qpStart))
    return BF_PF_INVALID_ARGUMENT;
  for(int q=qpStart;q<qpEnd;q++) {
    int status=BFStableBuild_fcmp(sltElmBF,eleIdx,q,0,NULL,NULL,bufM);
    if(status!=BF_PF_OK) return status;
    status=BFStableFactor_fcmp(bufM,pfMOut+q-qpStart,
        invMOut+(size_t)(q-qpStart)*invMQpStride,iwork,work,lwork, rwork,q,0);
    if(status!=BF_PF_OK) return status;
  }
  return BF_PF_OK;
}
