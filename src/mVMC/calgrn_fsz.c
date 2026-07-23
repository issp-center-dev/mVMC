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
 * Cauculate Green Functions
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "calgrn_fsz.h"
#include "backflow_nbody.h"
#ifndef _CALGRN_FSZ_SRC
#define _CALGRN_FSZ_SRC

double complex GreenFunc1BF_fsz(const int ri, const int rj, const int s, const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork);

double complex GreenFunc1BF_fsz_workspace(const int ri, const int rj, const int s,
                  const double complex ip, int *eleIdx, int *eleCfg, int *eleNum,
                  const int *eleProjCnt, int *eleSpn, int *projCntNew,
                  const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, double complex *pfBufM,
                  int *affected, int *hopIntWork, int hopIntWorkSize, int *pfIWork,
                  double complex *pfWork, double *pfRWork);

double complex GreenFunc1BF_fsz2_workspace(const int ri, const int rj,
                  const int s, const int t, const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum,
                  const int *eleProjCnt, int *eleSpn, int *projCntNew,
                  const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, double complex *pfBufM,
                  int *affected, int *hopIntWork, int hopIntWorkSize, int *pfIWork,
                  double complex *pfWork, double *pfRWork);

double complex GreenFunc2BF_fsz(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork);

double complex GreenFunc2BF_fsz2(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const int u, const int v,
                  const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork);

double complex GreenFunc2BF_fsz2WithProfile(
                  const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const int u, const int v,
                  const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork,
                  BFFSZC2DetailContext *detailProfile);
int GetGreenFuncBF_fsz_buffer_work_size(
    size_t *bufferComplexCount, size_t *pfUpdateIntCount,
    size_t *pfUpdateDoubleCount);

typedef struct {
  BFNBodyStatus status;
  int term;
  BFNBodyStage stage;
  int detail;
} BFFSZGreenNBodyFailure;

static void RecordBFFSZGreenNBodyFailure(
    BFFSZGreenNBodyFailure *failure, int term,
    const BFNBodyResult *result) {
  if(failure == NULL || result == NULL
     || result->status == BF_NBODY_OK
     || result->status == BF_NBODY_PHYSICAL_ZERO) {
    return;
  }
#pragma omp critical(BFFSZGreenNBodyFailureRecord)
  {
    if(failure->status == BF_NBODY_OK) {
      failure->status = result->status;
      failure->term = term;
      failure->stage = result->stage;
      failure->detail = result->detail;
    }
  }
}

static void AbortBFFSZGreenNBodyFailure(
    const char *context, const BFFSZGreenNBodyFailure *failure) {
  fprintf(stderr,
          "Error: %s failed: status=%d term=%d stage=%d detail=%d.\n",
          context, (int)failure->status, failure->term,
          (int)failure->stage, failure->detail);
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}


void CalculateGreenFunc_fsz(const double w, const double complex ip, int *eleIdx, int *eleCfg,
                         int *eleNum, int *eleSpn,int *eleProjCnt) {

  int idx,idx0,idx1;
  int k, offset, nbody;
  int ri,rj,s,rk,rl,t,u,v;
  double complex tmp;
  int *myEleIdx, *myEleNum, *myProjCntNew,*myEleSpn, *myRsi, *myRsj;
  double complex *myBuffer;
  double *myRWork;
  double two_pi = 2.0*acos(-1.0);
  int site_idx, spin_idx, rsi_idx, rsi;
  int x, y, z;
  double Wx, Wy, Wz;
  double weight;
  int maxNBodyG = (NBodyGMaxN > 2) ? NBodyGMaxN : 2;
  int needNBodyRWork = (NNBodyG > 0 && NBodyGMaxN > 2);

  RequestWorkSpaceThreadInt(Nsize+Nsize+Nsite2+NProj+2*maxNBodyG);
  RequestWorkSpaceThreadComplex(NQPFull+maxNBodyG*Nsize);
  if (needNBodyRWork) RequestWorkSpaceThreadDouble(LapackLWork);
  /* GreenFunc1: NQPFull, GreenFunc2/NBodyG: NQPFull+maxNBodyG*Nsize */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myEleSpn,myProjCntNew,myRsi,myRsj,myBuffer,myRWork,idx)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleSpn = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myRsi = GetWorkSpaceThreadInt(maxNBodyG);
    myRsj = GetWorkSpaceThreadInt(maxNBodyG);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull+maxNBodyG*Nsize);
    myRWork = needNBodyRWork ? GetWorkSpaceThreadDouble(LapackLWork) : NULL;

    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleSpn[idx] = eleSpn[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];

    #pragma omp master
    {StartTimer(50);}

    #pragma omp for private(idx,ri,rj,s,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NCisAjs;idx++) {
      ri = CisAjsIdx[idx][0];
      s  = CisAjsIdx[idx][1];
      rj = CisAjsIdx[idx][2];
      t  = CisAjsIdx[idx][3];
      if(s==t){
        tmp = GreenFunc1_fsz(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,
                       myProjCntNew,myBuffer);
      }else{
        tmp = GreenFunc1_fsz2(ri,rj,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,
                       myProjCntNew,myBuffer);
      }
      LocalCisAjs[idx] = tmp;
    }

    #pragma omp master
    {StopTimer(50);StartTimer(51);}

    #pragma omp for private(idx,ri,rj,s,rk,rl,t,tmp) schedule(dynamic)
    for(idx=0;idx<NCisAjsCktAltDC;idx++) {
      /*
      ri = CisAjsCktAltDCIdx[idx][0];
      rj = CisAjsCktAltDCIdx[idx][1];
      s  = CisAjsCktAltDCIdx[idx][2];
      rk = CisAjsCktAltDCIdx[idx][3];
      rl = CisAjsCktAltDCIdx[idx][4];
      t  = CisAjsCktAltDCIdx[idx][5];
      */
      ri = CisAjsCktAltDCIdx[idx][0];
      s  = CisAjsCktAltDCIdx[idx][1];
      rj = CisAjsCktAltDCIdx[idx][2];
      t  = CisAjsCktAltDCIdx[idx][3];
//
      rk = CisAjsCktAltDCIdx[idx][4];
      u  = CisAjsCktAltDCIdx[idx][5];
      rl = CisAjsCktAltDCIdx[idx][6];
      v  = CisAjsCktAltDCIdx[idx][7];

      if(s==t && u==v){
        tmp = GreenFunc2_fsz(ri,rj,rk,rl,s,u,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,
                       myProjCntNew,myBuffer);
      }else{
        tmp = GreenFunc2_fsz2(ri,rj,rk,rl,s,t,u,v,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,
                       myProjCntNew,myBuffer);
      }
      PhysCisAjsCktAltDC[idx] += w*tmp;
    }

    #pragma omp master
    {StopTimer(51);StartTimer(52);}

    #pragma omp for private(idx) nowait
    for(idx=0;idx<NCisAjs;idx++) {
      PhysCisAjs[idx] += w*LocalCisAjs[idx];
    }

    #pragma omp master
    {StopTimer(52);StartTimer(53);}

    #pragma omp for private(idx,idx0,idx1) nowait
    for(idx=0;idx<NCisAjsCktAlt;idx++) {
      idx0 = CisAjsCktAltIdx[idx][0];
      idx1 = CisAjsCktAltIdx[idx][1];
      PhysCisAjsCktAlt[idx] += w*LocalCisAjs[idx0]*conj(LocalCisAjs[idx1]);// TBC conj ok?
    }

    #pragma omp master
    {StopTimer(53);StartTimer(54);}

    for(idx=0;idx<NTwist;idx++) {
      #pragma omp single
      weight = 0.0;
      #pragma omp for private(site_idx, ri, s, rsi, x,y,z, Wx, Wy,Wz,tmp) reduction(+:weight)
      for(site_idx=0;site_idx<2*Nsite;site_idx++) {
        ri  = TwistIdx[idx][2*site_idx];
        s   = TwistIdx[idx][2*site_idx+1];
        rsi = ri + s*Nsite;

        x = LatticeIdx[ri][0];
        y = LatticeIdx[ri][1];
        z = LatticeIdx[ri][2];

        Wx = ParaTwist[idx][3*site_idx + 0];
        Wy = ParaTwist[idx][3*site_idx + 1];
        Wz = ParaTwist[idx][3*site_idx + 2];
        tmp = x*Wx + y*Wy + z*Wz;

        weight += tmp*(double)myEleNum[rsi];
      }
      #pragma omp single
      PhysTwist[idx] += w*cexp(I*two_pi*weight);
    }

    #pragma omp for private(idx,k,offset,nbody,tmp) schedule(dynamic)
    for(idx=0;idx<NNBodyG;idx++) {
      nbody = NBodyGN[idx];
      offset = NBodyGOffset[idx];
      for(k=0;k<nbody;k++) {
        myRsi[k] = NBodyGIdx[offset+k][0] + NBodyGIdx[offset+k][1]*Nsite;
        myRsj[k] = NBodyGIdx[offset+k][2] + NBodyGIdx[offset+k][3]*Nsite;
      }
      tmp = GreenFuncN_fsz(nbody,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,
                           eleProjCnt,myEleSpn,myBuffer,myProjCntNew,myRWork);
      PhysNBodyG[idx] += w*tmp;
    }

    #pragma omp master
    {StopTimer(54);}
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  if (needNBodyRWork) ReleaseWorkSpaceThreadDouble();
  return;
}

void CalculateGreenFuncBF_fsz(const double w, const double complex ip, int *eleIdx, int *eleCfg,
                         int *eleNum, int *eleSpn, int *eleProjCnt,
                         const int *eleProjBFCnt) {
  int idx,idx0,idx1;
  int ri,rj,s,t,rk,rl,u,v;
  double complex tmp;
  int *myEleIdx, *myEleCfg, *myEleNum, *myEleSpn;
  int *myProjCntNew, *myProjBFCntNew, *myAffected, *myHopIntWork;
  int *myPfIWork;
  double complex *myBuffer;
  double complex *myPfBufM, *myPfWork;
  double *myPfRWork;
  size_t pfUpdateIntCount,pfUpdateDoubleCount;
  size_t bfGreen1BufferSizeValue,bfGreen1ComplexSizeValue;
  int bfGreen1BufferSize,bfGreen1ComplexSize,bfPfIntSize,bfPfDoubleSize;
  long long bfSerialBaseIntSizeLL, bfGreen1IntSizeLL;
  int bfSerialBaseIntSize;
  int bfHopIntSize, bfSerialIntSize, bfGreen1IntSize;
  int reuseCensusCapacity=0, reuseCensusIntSize=0;
  int *reuseCensusKeys = NULL;
  int *thEleIdx, *thEleCfg, *thEleNum, *thEleSpn;
  int *thProjCntNew, *thProjBFCntNew, *thAffected, *thHopIntWork, *thPfIWork;
  double complex *thBuffer, *thPfBufM, *thPfWork;
  double *thPfRWork;
  BFFSZC2DetailContext c2DetailContext;
  BFFSZC2DetailContext *c2DetailProfile = NULL;
  BFFSZC2ReuseCensus reuseCensus;
  const int maxNBody = NBodyGMaxN > 2 ? NBodyGMaxN : 2;
  BFNBodyScratchSizes nbodyScratchSizes;
  BFFSZGreenNBodyFailure nbodyFailure = {
      BF_NBODY_OK, -1, BF_NBODY_STAGE_NONE, BF_NBODY_DETAIL_NONE
  };
  int workspaceBindingFailed = 0;
  int *thIntBase;
  double complex *thComplexBase;
  double *thDoubleBase;
  BFNBodyScratch thNBodyScratch;
  BFNBodyResult nbodyResult;
  int bindStatus;
  int k, offset, nbody;

  if(BFFSZC2DetailProfileEnabled) {
    InitBFFSZC2DetailContext(&c2DetailContext,
        BFFSZ_C2_DETAIL_SOURCE_MEASUREMENT);
    c2DetailProfile = &c2DetailContext;
  }

  if(NTwist > 0) {
    fprintf(stderr, "Error: CalculateGreenFuncBF_fsz does not support Twist yet.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(NNBodyG > 0
     && GetBFNBodyScratchSizes(maxNBody, 1, &nbodyScratchSizes)
          != BF_NBODY_OK) {
    nbodyFailure.status = BF_NBODY_WORKSPACE_ERROR;
    AbortBFFSZGreenNBodyFailure(
        "BackFlow FSZ NBodyG workspace sizing", &nbodyFailure);
    return;
  }
  if(GetSlaterElmBF_fsz_hop_int_work_size(&bfHopIntSize) != 0) {
    fprintf(stderr, "Error: invalid BF-FSZ Green hop integer workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(GetGreenFuncBF_fsz_buffer_work_size(
      &bfGreen1BufferSizeValue,&pfUpdateIntCount,&pfUpdateDoubleCount) != 0) {
    fprintf(stderr,"Error: invalid BF-FSZ Green Pfaffian-update workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  if((size_t)Nsize*(size_t)Nsize > SIZE_MAX-(size_t)LapackLWork
      || bfGreen1BufferSizeValue > SIZE_MAX-(size_t)Nsize*(size_t)Nsize-(size_t)LapackLWork) {
    fprintf(stderr,"Error: BF-FSZ Green complex workspace size overflow.\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  bfGreen1ComplexSizeValue = bfGreen1BufferSizeValue
      + (size_t)Nsize*(size_t)Nsize + (size_t)LapackLWork;
  if(bfGreen1BufferSizeValue > INT_MAX || bfGreen1ComplexSizeValue > INT_MAX) {
    fprintf(stderr,"Error: BF-FSZ Green complex workspace is too large.\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  bfGreen1BufferSize = (int)bfGreen1BufferSizeValue;
  bfGreen1ComplexSize = (int)bfGreen1ComplexSizeValue;
  bfPfIntSize = (Nsize > (int)pfUpdateIntCount) ? Nsize : (int)pfUpdateIntCount;
  bfPfDoubleSize = (LapackLWork > (int)pfUpdateDoubleCount)
      ? LapackLWork : (int)pfUpdateDoubleCount;
  bfSerialBaseIntSizeLL = 2LL*Nsize + 2LL*Nsite2 + (long long)NProj
      + 16LL*Nsite*Nrange;
  if(bfSerialBaseIntSizeLL < 0 || bfSerialBaseIntSizeLL > INT_MAX) {
    fprintf(stderr, "Error: invalid BF-FSZ Green base integer workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  bfSerialBaseIntSize = (int)bfSerialBaseIntSizeLL;
  if(Nsize > INT_MAX - bfSerialBaseIntSize) {
    fprintf(stderr, "Error: invalid BF-FSZ Green LAPACK workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  bfGreen1IntSizeLL = (long long)bfSerialBaseIntSize + bfPfIntSize
      + Nsize + bfHopIntSize;
  if(bfGreen1IntSizeLL < 0 || bfGreen1IntSizeLL > INT_MAX) {
    fprintf(stderr, "Error: invalid BF-FSZ Green integer workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  bfGreen1IntSize = (int)bfGreen1IntSizeLL;
  if(BFFSZC2ReuseCensusEnabled
      && GetBFFSZC2ReuseCensusWorkSize(
          (long long)NCisAjsCktAltDC,&reuseCensusCapacity,
          &reuseCensusIntSize) != 0) {
    fprintf(stderr,"Error: invalid BF-FSZ Green C2 reuse census workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  if(reuseCensusIntSize > INT_MAX-bfGreen1IntSize) {
    fprintf(stderr,"Error: BF-FSZ Green C2 reuse census workspace is too large.\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  bfSerialIntSize = bfGreen1IntSize+reuseCensusIntSize;

  RequestWorkSpaceInt(bfSerialIntSize);
  RequestWorkSpaceComplex(bfGreen1ComplexSize);
  RequestWorkSpaceDouble(bfPfDoubleSize);
  if(NNBodyG > 0) {
    RequestWorkSpaceThreadInt((int)nbodyScratchSizes.intCount);
    RequestWorkSpaceThreadComplex((int)nbodyScratchSizes.complexCount);
    RequestWorkSpaceThreadDouble((int)nbodyScratchSizes.doubleCount);
  } else {
    RequestWorkSpaceThreadInt(bfGreen1IntSize);
    RequestWorkSpaceThreadComplex(bfGreen1ComplexSize);
    RequestWorkSpaceThreadDouble(bfPfDoubleSize);
  }

  myEleIdx = GetWorkSpaceInt(Nsize);
  myEleCfg = GetWorkSpaceInt(Nsite2);
  myEleNum = GetWorkSpaceInt(Nsite2);
  myEleSpn = GetWorkSpaceInt(Nsize);
  myProjCntNew = GetWorkSpaceInt(NProj);
  myProjBFCntNew = GetWorkSpaceInt(16*Nsite*Nrange);
  myAffected = GetWorkSpaceInt(Nsize);
  myHopIntWork = GetWorkSpaceInt(bfHopIntSize);
  myPfIWork = GetWorkSpaceInt(bfPfIntSize);
  if(reuseCensusIntSize > 0) {
    reuseCensusKeys = GetWorkSpaceInt(reuseCensusIntSize);
  }
  myBuffer = GetWorkSpaceComplex(bfGreen1BufferSize);
  myPfBufM = GetWorkSpaceComplex(Nsize*Nsize);
  myPfWork = GetWorkSpaceComplex(LapackLWork);
  myPfRWork = GetWorkSpaceDouble(bfPfDoubleSize);

  if(BFFSZC2ReuseCensusEnabled) {
    if(InitBFFSZC2ReuseCensus(
        &reuseCensus,reuseCensusKeys,reuseCensusCapacity) != 0) {
      fprintf(stderr,"Error: BF-FSZ Green C2 reuse census initialization failed.\n");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    c2DetailContext.reuseCensus = &reuseCensus;
  }

  for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
  for(idx=0;idx<Nsite2;idx++) myEleCfg[idx] = eleCfg[idx];
  for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
  for(idx=0;idx<Nsize;idx++) myEleSpn[idx] = eleSpn[idx];

  #pragma omp parallel default(shared)                                      \
    private(thEleIdx,thEleCfg,thEleNum,thEleSpn,thProjCntNew,thProjBFCntNew, \
            thAffected,thHopIntWork,thPfIWork,thBuffer,thPfBufM,thPfWork,thPfRWork, \
            thIntBase,thComplexBase,thDoubleBase,thNBodyScratch,nbodyResult, \
            bindStatus,idx,ri,rj,s,t,tmp,k,offset,nbody)
  {
    if(NNBodyG > 0) {
      thIntBase =
          GetWorkSpaceThreadInt((int)nbodyScratchSizes.intCount);
      thComplexBase =
          GetWorkSpaceThreadComplex((int)nbodyScratchSizes.complexCount);
      thDoubleBase =
          GetWorkSpaceThreadDouble((int)nbodyScratchSizes.doubleCount);
      bindStatus = BindBFNBodyScratch(
          &nbodyScratchSizes,
          thIntBase, nbodyScratchSizes.intCount,
          thComplexBase, nbodyScratchSizes.complexCount,
          thDoubleBase, nbodyScratchSizes.doubleCount,
          &thNBodyScratch);
      if(bindStatus != BF_NBODY_OK) {
        const BFNBodyResult bindResult = {
            BF_NBODY_WORKSPACE_ERROR, BF_NBODY_STAGE_NONE,
            bindStatus, 0, 0.0+0.0*I
        };
        RecordBFFSZGreenNBodyFailure(&nbodyFailure, -1, &bindResult);
#pragma omp atomic write
        workspaceBindingFailed = 1;
      } else {
        thEleIdx = thNBodyScratch.eleIdx;
        thEleCfg = thNBodyScratch.eleCfg;
        thEleNum = thNBodyScratch.eleNum;
        thEleSpn = thNBodyScratch.eleSpn;
        thProjCntNew = thNBodyScratch.projCnt;
        thProjBFCntNew = thNBodyScratch.projBFCnt;
        thPfIWork = thNBodyScratch.pfIWork;
        thAffected = thNBodyScratch.affected;
        thHopIntWork = thNBodyScratch.hopIntWork;
        thBuffer = thNBodyScratch.greenBuffer;
        thPfBufM = thNBodyScratch.pfBufM;
        thPfWork = thNBodyScratch.pfWork;
        thPfRWork = thNBodyScratch.pfRWork;
      }
    } else {
      thEleIdx = GetWorkSpaceThreadInt(Nsize);
      thEleCfg = GetWorkSpaceThreadInt(Nsite2);
      thEleNum = GetWorkSpaceThreadInt(Nsite2);
      thEleSpn = GetWorkSpaceThreadInt(Nsize);
      thProjCntNew = GetWorkSpaceThreadInt(NProj);
      thProjBFCntNew = GetWorkSpaceThreadInt(16*Nsite*Nrange);
      thPfIWork = GetWorkSpaceThreadInt(bfPfIntSize);
      thAffected = GetWorkSpaceThreadInt(Nsize);
      thHopIntWork = GetWorkSpaceThreadInt(bfHopIntSize);
      thBuffer = GetWorkSpaceThreadComplex(bfGreen1BufferSize);
      thPfBufM = GetWorkSpaceThreadComplex(Nsize*Nsize);
      thPfWork = GetWorkSpaceThreadComplex(LapackLWork);
      thPfRWork = GetWorkSpaceThreadDouble(bfPfDoubleSize);
    }

#pragma omp barrier
    if(!workspaceBindingFailed) {
#pragma loop noalias
      for(idx=0;idx<Nsize;idx++) thEleIdx[idx] = eleIdx[idx];
#pragma loop noalias
      for(idx=0;idx<Nsite2;idx++) thEleCfg[idx] = eleCfg[idx];
#pragma loop noalias
      for(idx=0;idx<Nsite2;idx++) thEleNum[idx] = eleNum[idx];
#pragma loop noalias
      for(idx=0;idx<Nsize;idx++) thEleSpn[idx] = eleSpn[idx];

#pragma omp barrier
#pragma omp master
      { StartTimer(50); }
#pragma omp barrier

#pragma omp for schedule(dynamic)
      for(idx=0;idx<NCisAjs;idx++) {
        ri = CisAjsIdx[idx][0];
        s  = CisAjsIdx[idx][1];
        rj = CisAjsIdx[idx][2];
        t  = CisAjsIdx[idx][3];
        tmp = GreenFunc1BF_fsz2_workspace(
            ri,rj,s,t,ip,
            thEleIdx,thEleCfg,thEleNum,eleProjCnt,thEleSpn,
            thProjCntNew,eleProjBFCnt,thProjBFCntNew,
            thBuffer,thPfBufM,thAffected,thHopIntWork,
            bfHopIntSize,thPfIWork,thPfWork,thPfRWork);
        LocalCisAjs[idx] = tmp;
      }

#pragma omp barrier
#pragma omp master
      { StopTimer(50); }

      if(NNBodyG > 0) {
#pragma omp barrier
#pragma omp master
        { StartTimer(54); }
#pragma omp barrier
#pragma omp for private(idx,k,offset,nbody,nbodyResult) schedule(dynamic)
        for(idx=0;idx<NNBodyG;idx++) {
          nbody = NBodyGN[idx];
          offset = NBodyGOffset[idx];
          for(k=0;k<nbody;k++) {
            thNBodyScratch.inputRsi[k] =
                NBodyGIdx[offset+k][0] + NBodyGIdx[offset+k][1]*Nsite;
            thNBodyScratch.inputRsj[k] =
                NBodyGIdx[offset+k][2] + NBodyGIdx[offset+k][3]*Nsite;
          }
          nbodyResult = GreenFuncNBF_fsz(
              nbody,
              thNBodyScratch.inputRsi, thNBodyScratch.inputRsj, ip,
              eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn,
              eleProjBFCnt, &thNBodyScratch);
          if(nbodyResult.status == BF_NBODY_OK
             || nbodyResult.status == BF_NBODY_PHYSICAL_ZERO) {
            PhysNBodyG[idx] += w*nbodyResult.value;
          } else {
            RecordBFFSZGreenNBodyFailure(
                &nbodyFailure, idx, &nbodyResult);
          }
        }
#pragma omp barrier
#pragma omp master
        { StopTimer(54); }
      }
    }
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  ReleaseWorkSpaceThreadDouble();
  if(nbodyFailure.status != BF_NBODY_OK) {
    ReleaseWorkSpaceInt();
    ReleaseWorkSpaceComplex();
    ReleaseWorkSpaceDouble();
    AbortBFFSZGreenNBodyFailure(
        "BackFlow FSZ NBodyG evaluation", &nbodyFailure);
    return;
  }

  StartTimer(51);
  for(idx=0;idx<NCisAjsCktAltDC;idx++) {
    ri = CisAjsCktAltDCIdx[idx][0];
    s  = CisAjsCktAltDCIdx[idx][1];
    rj = CisAjsCktAltDCIdx[idx][2];
    t  = CisAjsCktAltDCIdx[idx][3];
    rk = CisAjsCktAltDCIdx[idx][4];
    u  = CisAjsCktAltDCIdx[idx][5];
    rl = CisAjsCktAltDCIdx[idx][6];
    v  = CisAjsCktAltDCIdx[idx][7];

    tmp = GreenFunc2BF_fsz2WithProfile(ri,rj,rk,rl,s,t,u,v,ip,
                           myEleIdx,myEleCfg,myEleNum,eleProjCnt,myEleSpn,
                           myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer,
                           myAffected,myHopIntWork,bfHopIntSize,
                           myPfIWork,myPfRWork,myPfBufM,myPfWork,c2DetailProfile);
    PhysCisAjsCktAltDC[idx] += w*tmp;
  }
  StopTimer(51);
  if(c2DetailProfile != NULL) MergeBFFSZC2DetailContext(c2DetailProfile);
  if(BFFSZC2ReuseCensusEnabled) {
    MergeBFFSZC2ReuseCensus(BFFSZ_C2_REUSE_SCOPE_MEASUREMENT,&reuseCensus);
  }

  StartTimer(52);
  for(idx=0;idx<NCisAjs;idx++) {
    PhysCisAjs[idx] += w*LocalCisAjs[idx];
  }
  StopTimer(52);

  StartTimer(53);
  for(idx=0;idx<NCisAjsCktAlt;idx++) {
    idx0 = CisAjsCktAltIdx[idx][0];
    idx1 = CisAjsCktAltIdx[idx][1];
    PhysCisAjsCktAlt[idx] += w*LocalCisAjs[idx0]*conj(LocalCisAjs[idx1]);
  }
  StopTimer(53);

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceComplex();
  ReleaseWorkSpaceDouble();
  return;
}
#endif
