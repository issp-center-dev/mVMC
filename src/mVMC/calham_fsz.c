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
 * calculate Hamiltonian
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "calham_fsz.h"
#include "backflow_nbody.h"

#pragma once

double complex GreenFuncN_fsz(const int n, int *rsi, int *rsj,
                  const double complex ip, int *eleIdx, const int *eleCfg,
                  int *eleNum, const int *eleProjCnt, int *eleSpn,
                  double complex *buffer, int *bufferInt, double *rwork);

double complex GreenFunc1BF_fsz(const int ri, const int rj, const int s, const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork);

double complex GreenFunc1BF_fsz2(const int ri, const int rj, const int s, const int t,
                  const double complex ip,
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
} BFFSZHamiltonianNBodyFailure;

static void RecordBFFSZHamiltonianNBodyFailure(
    BFFSZHamiltonianNBodyFailure *failure, int term,
    const BFNBodyResult *result) {
  if(failure == NULL || result == NULL
     || result->status == BF_NBODY_OK
     || result->status == BF_NBODY_PHYSICAL_ZERO
     || failure->status != BF_NBODY_OK) {
    return;
  }
  failure->status = result->status;
  failure->term = term;
  failure->stage = result->stage;
  failure->detail = result->detail;
}

static void AbortBFFSZHamiltonianNBodyFailure(
    const char *context, const BFFSZHamiltonianNBodyFailure *failure) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) {
    fprintf(stderr,
            "Error: %s failed: status=%d term=%d stage=%s(%d) "
            "detail=%s(%d).\n",
            context, (int)failure->status, failure->term,
            BFNBodyStageName(failure->stage), (int)failure->stage,
            BFNBodyDetailName(failure->detail), failure->detail);
  }
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

double CalculateSz_fsz(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt,int *eleSpn) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  double  Sz=0.0;
  int ri;
    
  Sz = 0.0;
//#pragma omp for default(none) firstprivate(Nsite,n0,n1) private(idx,ri) reduction(+:Sz)
  for(ri=0;ri<Nsite;ri++) {
    Sz += n0[ri]-n1[ri];
  }
  return Sz;
}



double complex CalculateHamiltonian_fsz(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt,int *eleSpn) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  double complex e=0.0, tmp;
  int idx;
  int ri,rj,s,rk,rl,t,u,v;
  int k, offset, nbody;
  int *myEleIdx, *myEleNum, *myProjCntNew,*myEleSpn;
  int *myNBodyRsi, *myNBodyRsj;
  double complex *myBuffer;
  double *myRWork;
  double complex myEnergy;
  int maxNBodyInterAll = (NBodyInterAllMaxN > 2) ? NBodyInterAllMaxN : 2;
  int needNBodyInterAllRWork = (NNBodyInterAll > 0 && NBodyInterAllMaxN > 2);

  RequestWorkSpaceThreadInt(Nsize+Nsize+Nsite2+NProj+2*maxNBodyInterAll);
  RequestWorkSpaceThreadComplex(NQPFull+maxNBodyInterAll*Nsize);
  if (needNBodyInterAllRWork) RequestWorkSpaceThreadDouble(LapackLWork);
  /* GreenFunc1: NQPFull, GreenFunc2/NBodyInterAll: NQPFull+maxNBodyInterAll*Nsize */

  /*
#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myProjCntNew,myBuffer,myEnergy,idx)  \
  reduction(+:e)
  */
#pragma omp parallel default(none)                                      \
  private(myEleIdx,myEleSpn,myEleNum,myProjCntNew,myNBodyRsi,myNBodyRsj,myBuffer,myRWork,myEnergy, idx, ri, rj, rk, rl, s, t,u,v) \
  firstprivate(ip, Nsize, Nsite2, NProj, NQPFull, NCoulombIntra, CoulombIntra, ParaCoulombIntra, \
               NCoulombInter, CoulombInter, ParaCoulombInter, NHundCoupling, HundCoupling, ParaHundCoupling, \
               NTransfer, Transfer, ParaTransfer, NPairHopping, PairHopping, ParaPairHopping, \
               NExchangeCoupling, ExchangeCoupling, ParaExchangeCoupling, NInterAll, InterAll, ParaInterAll, \
               Nsite, maxNBodyInterAll, needNBodyInterAllRWork, LapackLWork, \
               NNBodyInterAll, NBodyInterAllN, NBodyInterAllOffset, NBodyInterAllIdx, ParaNBodyInterAll, n0, n1) \
  shared(eleCfg, eleProjCnt, eleIdx, eleNum,eleSpn) reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleSpn = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myNBodyRsi = GetWorkSpaceThreadInt(maxNBodyInterAll);
    myNBodyRsj = GetWorkSpaceThreadInt(maxNBodyInterAll);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull+maxNBodyInterAll*Nsize);
    myRWork = needNBodyInterAllRWork ? GetWorkSpaceThreadDouble(LapackLWork) : NULL;

    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleSpn[idx] = eleSpn[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
    #pragma omp barrier
    
    myEnergy = 0.0;

    #pragma omp master
    {StartTimer(70);}

    /* CoulombIntra */
    #pragma omp for private(idx,ri) nowait
    for(idx=0;idx<NCoulombIntra;idx++) {
      ri = CoulombIntra[idx];
      myEnergy += ParaCoulombIntra[idx] * n0[ri] * n1[ri];
    }
  
    /* CoulombInter */
    #pragma omp for private(idx,ri,rj) nowait
    for(idx=0;idx<NCoulombInter;idx++) {
      ri = CoulombInter[idx][0];
      rj = CoulombInter[idx][1];
      myEnergy += ParaCoulombInter[idx] * (n0[ri]+n1[ri]) * (n0[rj]+n1[rj]);
    }

    /* HundCoupling */
    #pragma omp for private(idx,ri,rj) nowait
    for(idx=0;idx<NHundCoupling;idx++) {
      ri = HundCoupling[idx][0];
      rj = HundCoupling[idx][1];
      myEnergy -= ParaHundCoupling[idx] * (n0[ri]*n0[rj] + n1[ri]*n1[rj]);
      /* Caution: negative sign */
    }

    #pragma omp master
    {StopTimer(70);StartTimer(71);}

    /* Transfer */
    #pragma omp for private(idx,ri,rj,s,t) schedule(dynamic) nowait
    for(idx=0;idx<NTransfer;idx++) {
      ri = Transfer[idx][0];
      s  = Transfer[idx][1];
      rj = Transfer[idx][2];
      t  = Transfer[idx][3];
      if(s==t){
        myEnergy -= ParaTransfer[idx]
         * GreenFunc1_fsz(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      }else{ 
        myEnergy -= ParaTransfer[idx]
        * GreenFunc1_fsz2(ri,rj,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      }
      /* Caution: negative sign */
    }

    #pragma omp master
    {StopTimer(71);StartTimer(72);}

    /* Pair Hopping */
    #pragma omp for private(idx,ri,rj) schedule(dynamic) nowait
    for(idx=0;idx<NPairHopping;idx++) {
      ri = PairHopping[idx][0];
      rj = PairHopping[idx][1];
    
      myEnergy += ParaPairHopping[idx]
        * GreenFunc2_fsz(ri,rj,ri,rj,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
    }

    /* Exchange Coupling */
    #pragma omp for private(idx,ri,rj,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NExchangeCoupling;idx++) {
      ri = ExchangeCoupling[idx][0];
      rj = ExchangeCoupling[idx][1];
    
      tmp =  GreenFunc2_fsz(ri,rj,rj,ri,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      //printf("idx=%d ri=%d rj=%d:  tmp=%lf \n",idx,ri,rj,creal(tmp));
      tmp += GreenFunc2_fsz(ri,rj,rj,ri,1,0,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      //tmp =  GreenFunc2_fsz2(ri,rj,rj,ri,0,0,1,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      //tmp += GreenFunc2_fsz2(ri,rj,rj,ri,1,1,0,0,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      myEnergy += ParaExchangeCoupling[idx] * tmp;
    }

    /* Inter All */
    #pragma omp for private(idx,ri,rj,s,rk,rl,t) schedule(dynamic) nowait
    for(idx=0;idx<NInterAll;idx++) {
      ri = InterAll[idx][0];
      s  = InterAll[idx][1];
      rj = InterAll[idx][2];
      t  = InterAll[idx][3];
      rk = InterAll[idx][4];
      u  = InterAll[idx][5];
      rl = InterAll[idx][6];
      v  = InterAll[idx][7];
      
      if(s==t && u==v){
        myEnergy += ParaInterAll[idx]
          * GreenFunc2_fsz(ri,rj,rk,rl,s,u,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      }else{
        myEnergy += ParaInterAll[idx]
          * GreenFunc2_fsz2(ri,rj,rk,rl,s,t,u,v,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      } 
    }

    /* NBodyInterAll */
    #pragma omp for private(idx,k,offset,nbody,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NNBodyInterAll;idx++) {
      nbody = NBodyInterAllN[idx];
      offset = NBodyInterAllOffset[idx];
      for(k=0;k<nbody;k++) {
        myNBodyRsi[k] = NBodyInterAllIdx[offset+k][0]
                      + NBodyInterAllIdx[offset+k][1]*Nsite;
        myNBodyRsj[k] = NBodyInterAllIdx[offset+k][2]
                      + NBodyInterAllIdx[offset+k][3]*Nsite;
      }
      tmp = GreenFuncN_fsz(nbody,myNBodyRsi,myNBodyRsj,ip,myEleIdx,eleCfg,
                           myEleNum,eleProjCnt,myEleSpn,myBuffer,
                           myProjCntNew,myRWork);
      myEnergy += ParaNBodyInterAll[idx] * tmp;
    }

    #pragma omp master
    {StopTimer(72);}

    e += myEnergy;
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  if (needNBodyInterAllRWork) ReleaseWorkSpaceThreadDouble();
  return e;
}

double complex CalculateHamiltonianBF_fsz(const double complex ip, int *eleIdx, int *eleCfg,
                             int *eleNum, const int *eleProjCnt, int *eleSpn,
                             const int *eleProjBFCnt) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  double complex e=0.0,tmp;
  int idx,k,offset,nbody;
  int ri,rj,s,t,rk,rl,u,v;
  int *myEleIdx, *myEleCfg, *myEleNum, *myEleSpn;
  int *myProjCntNew, *myProjBFCntNew, *myAffected, *myHopIntWork;
  int *myPfIWork;
  double complex *myBuffer, *myPfBufM, *myPfWork;
  double *myPfRWork;
  size_t pfUpdateIntCount,pfUpdateDoubleCount;
  size_t bfBufferSize,bfComplexSize,nsizeSquared;
  long long bfBaseIntSizeLL,bfIntSizeLL;
  int bfBaseIntSize, bfHopIntSize, bfIntSize, bfPfIntSize, bfPfDoubleSize;
  int reuseCensusCapacity=0, reuseCensusIntSize=0;
  int *reuseCensusKeys = NULL;
  BFFSZC2DetailContext pairHopDetailContext;
  BFFSZC2DetailContext exchangeDetailContext;
  BFFSZC2DetailContext interAllDetailContext;
  BFFSZC2DetailContext *pairHopDetailProfile = NULL;
  BFFSZC2DetailContext *exchangeDetailProfile = NULL;
  BFFSZC2DetailContext *interAllDetailProfile = NULL;
  BFFSZC2DetailTermContext termDetailContext;
  BFFSZC2DetailTermContext *termDetailProfile = NULL;
  BFFSZC2ReuseCensus reuseCensus;
  double termStart = 0.0;
  const int maxNBody =
      NBodyInterAllMaxN > 2 ? NBodyInterAllMaxN : 2;
  BFNBodyScratchSizes nbodyScratchSizes;
  BFNBodyScratch nbodyScratch;
  BFNBodyResult nbodyResult;
  BFFSZHamiltonianNBodyFailure nbodyFailure = {
      BF_NBODY_OK, -1, BF_NBODY_STAGE_NONE, BF_NBODY_DETAIL_NONE
  };
  int *nbodyIntBase;
  double complex *nbodyComplexBase;
  double *nbodyDoubleBase;
  int bindStatus;
  int nbodyIntSize;

  if(BFFSZC2DetailProfileEnabled) {
    InitBFFSZC2DetailContext(&pairHopDetailContext,
        BFFSZ_C2_DETAIL_SOURCE_PAIR_HOP);
    InitBFFSZC2DetailContext(&exchangeDetailContext,
        BFFSZ_C2_DETAIL_SOURCE_EXCHANGE);
    InitBFFSZC2DetailContext(&interAllDetailContext,
        BFFSZ_C2_DETAIL_SOURCE_INTER_ALL);
    InitBFFSZC2DetailTermContext(&termDetailContext);
    pairHopDetailProfile = &pairHopDetailContext;
    exchangeDetailProfile = &exchangeDetailContext;
    interAllDetailProfile = &interAllDetailContext;
    termDetailProfile = &termDetailContext;
  }

  if(NNBodyInterAll > 0
     && GetBFNBodyScratchSizes(maxNBody, 1, &nbodyScratchSizes)
          != BF_NBODY_OK) {
    nbodyFailure.status = BF_NBODY_WORKSPACE_ERROR;
    AbortBFFSZHamiltonianNBodyFailure(
        "BackFlow FSZ NBodyInterAll workspace sizing", &nbodyFailure);
    return 0.0+0.0*I;
  }
  if(GetSlaterElmBF_fsz_hop_int_work_size(&bfHopIntSize) != 0) {
    fprintf(stderr, "Error: invalid BF-FSZ Hamiltonian hop integer workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(GetGreenFuncBF_fsz_buffer_work_size(
      &bfBufferSize,&pfUpdateIntCount,&pfUpdateDoubleCount) != 0) {
    fprintf(stderr,"Error: invalid BF-FSZ Hamiltonian Pfaffian workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  nsizeSquared = (size_t)Nsize*(size_t)Nsize;
  if(nsizeSquared > SIZE_MAX-(size_t)LapackLWork
      || bfBufferSize > SIZE_MAX-nsizeSquared-(size_t)LapackLWork) {
    fprintf(stderr,"Error: BF-FSZ Hamiltonian complex workspace size overflow.\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  bfComplexSize = bfBufferSize+nsizeSquared+(size_t)LapackLWork;
  if(bfComplexSize > INT_MAX) {
    fprintf(stderr,"Error: BF-FSZ Hamiltonian complex workspace is too large.\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  bfPfIntSize = (Nsize > (int)pfUpdateIntCount)
      ? Nsize : (int)pfUpdateIntCount;
  bfPfDoubleSize = (LapackLWork > (int)pfUpdateDoubleCount)
      ? LapackLWork : (int)pfUpdateDoubleCount;
  bfBaseIntSizeLL = 2LL*Nsize + 2LL*Nsite2 + (long long)NProj
      + 16LL*Nsite*Nrange;
  if(bfBaseIntSizeLL < 0 || bfBaseIntSizeLL > INT_MAX) {
    fprintf(stderr, "Error: invalid BF-FSZ Hamiltonian base integer workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  bfBaseIntSize = (int)bfBaseIntSizeLL;
  if(BFFSZC2ReuseCensusEnabled
      && (NPairHopping < 0 || NExchangeCoupling < 0 || NInterAll < 0
          || GetBFFSZC2ReuseCensusWorkSize(
              (long long)NPairHopping+2LL*NExchangeCoupling+NInterAll,
              &reuseCensusCapacity,&reuseCensusIntSize) != 0)) {
    fprintf(stderr,"Error: invalid BF-FSZ Hamiltonian C2 reuse census workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  bfIntSizeLL = (long long)bfBaseIntSize+Nsize+bfHopIntSize+bfPfIntSize
      + reuseCensusIntSize;
  if(bfIntSizeLL < 0 || bfIntSizeLL > INT_MAX) {
    fprintf(stderr, "Error: invalid BF-FSZ Hamiltonian integer workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  bfIntSize = (int)bfIntSizeLL;

  if(NNBodyInterAll > 0) {
    if(nbodyScratchSizes.intCount
         > (size_t)INT_MAX-(size_t)reuseCensusIntSize) {
      nbodyFailure.status = BF_NBODY_WORKSPACE_ERROR;
      AbortBFFSZHamiltonianNBodyFailure(
          "BackFlow FSZ NBodyInterAll workspace sizing", &nbodyFailure);
      return 0.0+0.0*I;
    }
    nbodyIntSize =
        (int)nbodyScratchSizes.intCount+reuseCensusIntSize;
    RequestWorkSpaceInt(nbodyIntSize);
    RequestWorkSpaceComplex((int)nbodyScratchSizes.complexCount);
    RequestWorkSpaceDouble((int)nbodyScratchSizes.doubleCount);
    nbodyIntBase =
        GetWorkSpaceInt((int)nbodyScratchSizes.intCount);
    if(reuseCensusIntSize > 0) {
      reuseCensusKeys = GetWorkSpaceInt(reuseCensusIntSize);
    }
    nbodyComplexBase =
        GetWorkSpaceComplex((int)nbodyScratchSizes.complexCount);
    nbodyDoubleBase =
        GetWorkSpaceDouble((int)nbodyScratchSizes.doubleCount);
    bindStatus = BindBFNBodyScratch(
        &nbodyScratchSizes,
        nbodyIntBase, nbodyScratchSizes.intCount,
        nbodyComplexBase, nbodyScratchSizes.complexCount,
        nbodyDoubleBase, nbodyScratchSizes.doubleCount,
        &nbodyScratch);
    if(bindStatus != BF_NBODY_OK) {
      const BFNBodyResult bindResult = {
          BF_NBODY_WORKSPACE_ERROR, BF_NBODY_STAGE_NONE,
          bindStatus, 0, 0.0+0.0*I
      };
      RecordBFFSZHamiltonianNBodyFailure(
          &nbodyFailure, -1, &bindResult);
      ReleaseWorkSpaceInt();
      ReleaseWorkSpaceComplex();
      ReleaseWorkSpaceDouble();
      AbortBFFSZHamiltonianNBodyFailure(
          "BackFlow FSZ NBodyInterAll workspace binding", &nbodyFailure);
      return 0.0+0.0*I;
    }
    myEleIdx = nbodyScratch.eleIdx;
    myEleCfg = nbodyScratch.eleCfg;
    myEleNum = nbodyScratch.eleNum;
    myEleSpn = nbodyScratch.eleSpn;
    myProjCntNew = nbodyScratch.projCnt;
    myProjBFCntNew = nbodyScratch.projBFCnt;
    myAffected = nbodyScratch.affected;
    myHopIntWork = nbodyScratch.hopIntWork;
    myPfIWork = nbodyScratch.pfIWork;
    myBuffer = nbodyScratch.greenBuffer;
    myPfBufM = nbodyScratch.pfBufM;
    myPfWork = nbodyScratch.pfWork;
    myPfRWork = nbodyScratch.pfRWork;
  } else {
    RequestWorkSpaceInt(bfIntSize);
    RequestWorkSpaceComplex((int)bfComplexSize);
    RequestWorkSpaceDouble(bfPfDoubleSize);

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
    myBuffer = GetWorkSpaceComplex((int)bfBufferSize);
    myPfBufM = GetWorkSpaceComplex(Nsize*Nsize);
    myPfWork = GetWorkSpaceComplex(LapackLWork);
    myPfRWork = GetWorkSpaceDouble(bfPfDoubleSize);
  }

  if(BFFSZC2ReuseCensusEnabled) {
    if(InitBFFSZC2ReuseCensus(
        &reuseCensus,reuseCensusKeys,reuseCensusCapacity) != 0) {
      fprintf(stderr,"Error: BF-FSZ Hamiltonian C2 reuse census initialization failed.\n");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    pairHopDetailContext.reuseCensus = &reuseCensus;
    exchangeDetailContext.reuseCensus = &reuseCensus;
    interAllDetailContext.reuseCensus = &reuseCensus;
  }

  for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
  for(idx=0;idx<Nsite2;idx++) myEleCfg[idx] = eleCfg[idx];
  for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
  for(idx=0;idx<Nsize;idx++) myEleSpn[idx] = eleSpn[idx];

  if(termDetailProfile != NULL) termStart = BFFSZC2DetailMonotonicSeconds();
  for(idx=0;idx<NCoulombIntra;idx++) {
    ri = CoulombIntra[idx];
    e += ParaCoulombIntra[idx] * n0[ri] * n1[ri];
  }

  for(idx=0;idx<NCoulombInter;idx++) {
    ri = CoulombInter[idx][0];
    rj = CoulombInter[idx][1];
    e += ParaCoulombInter[idx] * (n0[ri]+n1[ri]) * (n0[rj]+n1[rj]);
  }

  for(idx=0;idx<NHundCoupling;idx++) {
    ri = HundCoupling[idx][0];
    rj = HundCoupling[idx][1];
    e -= ParaHundCoupling[idx] * (n0[ri]*n0[rj] + n1[ri]*n1[rj]);
  }
  if(termDetailProfile != NULL) {
    termDetailProfile->seconds[BFFSZ_C2_DETAIL_TERM_NUMBER]
        += BFFSZC2DetailMonotonicSeconds()-termStart;
    termStart = BFFSZC2DetailMonotonicSeconds();
  }

  for(idx=0;idx<NTransfer;idx++) {
    ri = Transfer[idx][0];
    s  = Transfer[idx][1];
    rj = Transfer[idx][2];
    t  = Transfer[idx][3];
    e -= ParaTransfer[idx]
      * GreenFunc1BF_fsz2(ri,rj,s,t,ip,myEleIdx,myEleCfg,myEleNum,eleProjCnt,myEleSpn,
                         myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer,
                         myAffected,myHopIntWork,bfHopIntSize,
                         myPfIWork,myPfRWork,myPfBufM,myPfWork);
  }
  if(termDetailProfile != NULL) {
    termDetailProfile->seconds[BFFSZ_C2_DETAIL_TERM_TRANSFER]
        += BFFSZC2DetailMonotonicSeconds()-termStart;
    termStart = BFFSZC2DetailMonotonicSeconds();
  }

  for(idx=0;idx<NPairHopping;idx++) {
    ri = PairHopping[idx][0];
    rj = PairHopping[idx][1];
    e += ParaPairHopping[idx]
      * GreenFunc2BF_fsz2WithProfile(ri,rj,ri,rj,0,0,1,1,ip,
                          myEleIdx,myEleCfg,myEleNum,eleProjCnt,myEleSpn,
                          myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer,
                          myAffected,myHopIntWork,bfHopIntSize,
                          myPfIWork,myPfRWork,myPfBufM,myPfWork,
                          pairHopDetailProfile);
  }
  if(termDetailProfile != NULL) {
    termDetailProfile->seconds[BFFSZ_C2_DETAIL_TERM_PAIR_HOP]
        += BFFSZC2DetailMonotonicSeconds()-termStart;
    termStart = BFFSZC2DetailMonotonicSeconds();
  }

  for(idx=0;idx<NExchangeCoupling;idx++) {
    ri = ExchangeCoupling[idx][0];
    rj = ExchangeCoupling[idx][1];
    tmp = GreenFunc2BF_fsz2WithProfile(ri,rj,rj,ri,0,0,1,1,ip,
                            myEleIdx,myEleCfg,myEleNum,eleProjCnt,myEleSpn,
                            myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer,
                            myAffected,myHopIntWork,bfHopIntSize,
                            myPfIWork,myPfRWork,myPfBufM,myPfWork,
                            exchangeDetailProfile);
    tmp += GreenFunc2BF_fsz2WithProfile(ri,rj,rj,ri,1,1,0,0,ip,
                             myEleIdx,myEleCfg,myEleNum,eleProjCnt,myEleSpn,
                             myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer,
                             myAffected,myHopIntWork,bfHopIntSize,
                             myPfIWork,myPfRWork,myPfBufM,myPfWork,
                             exchangeDetailProfile);
    e += ParaExchangeCoupling[idx]*tmp;
  }
  if(termDetailProfile != NULL) {
    termDetailProfile->seconds[BFFSZ_C2_DETAIL_TERM_EXCHANGE]
        += BFFSZC2DetailMonotonicSeconds()-termStart;
    termStart = BFFSZC2DetailMonotonicSeconds();
  }

  for(idx=0;idx<NInterAll;idx++) {
    ri = InterAll[idx][0];
    s  = InterAll[idx][1];
    rj = InterAll[idx][2];
    t  = InterAll[idx][3];
    rk = InterAll[idx][4];
    u  = InterAll[idx][5];
    rl = InterAll[idx][6];
    v  = InterAll[idx][7];
    e += ParaInterAll[idx]
      * GreenFunc2BF_fsz2WithProfile(ri,rj,rk,rl,s,t,u,v,ip,
                          myEleIdx,myEleCfg,myEleNum,eleProjCnt,myEleSpn,
                          myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer,
                          myAffected,myHopIntWork,bfHopIntSize,
                          myPfIWork,myPfRWork,myPfBufM,myPfWork,
                          interAllDetailProfile);
  }
  if(termDetailProfile != NULL) {
    termDetailProfile->seconds[BFFSZ_C2_DETAIL_TERM_INTER_ALL]
        += BFFSZC2DetailMonotonicSeconds()-termStart;
  }
  for(idx=0;idx<NNBodyInterAll;idx++) {
    nbodyScratch.termIndex = idx;
    nbody = NBodyInterAllN[idx];
    offset = NBodyInterAllOffset[idx];
    for(k=0;k<nbody;k++) {
      nbodyScratch.inputRsi[k] =
          NBodyInterAllIdx[offset+k][0]
          + NBodyInterAllIdx[offset+k][1]*Nsite;
      nbodyScratch.inputRsj[k] =
          NBodyInterAllIdx[offset+k][2]
          + NBodyInterAllIdx[offset+k][3]*Nsite;
    }
    nbodyResult = GreenFuncNBF_fsz(
        nbody, nbodyScratch.inputRsi, nbodyScratch.inputRsj, ip,
        eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn,
        eleProjBFCnt, &nbodyScratch);
    if(nbodyResult.status == BF_NBODY_OK
       || nbodyResult.status == BF_NBODY_PHYSICAL_ZERO) {
      e += ParaNBodyInterAll[idx]*nbodyResult.value;
    } else {
      RecordBFFSZHamiltonianNBodyFailure(
          &nbodyFailure, idx, &nbodyResult);
      break;
    }
  }
  if(termDetailProfile != NULL) {
    MergeBFFSZC2DetailContext(pairHopDetailProfile);
    MergeBFFSZC2DetailContext(exchangeDetailProfile);
    MergeBFFSZC2DetailContext(interAllDetailProfile);
    MergeBFFSZC2DetailTermContext(termDetailProfile);
  }
  if(BFFSZC2ReuseCensusEnabled) {
    MergeBFFSZC2ReuseCensus(BFFSZ_C2_REUSE_SCOPE_HAMILTONIAN,&reuseCensus);
  }

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceComplex();
  ReleaseWorkSpaceDouble();
  if(nbodyFailure.status != BF_NBODY_OK) {
    AbortBFFSZHamiltonianNBodyFailure(
        "BackFlow FSZ NBodyInterAll evaluation", &nbodyFailure);
  }
  return e;
}

/* Calculate the CoulombIntra, CoulombInter, Hund terms, */
/* which can be calculated by number operators. */
/* This function will be used in the Lanczos mode */
double complex CalculateHamiltonian0_fsz(const int *eleNum) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  double complex e=0.0;
  int idx;
  int ri,rj;
  double complex myEnergy;

#pragma omp parallel default(shared)\
  private(myEnergy) reduction(+:e)
  {
    myEnergy = 0.0;

    /* CoulombIntra */
    #pragma omp for private(idx,ri)
    for(idx=0;idx<NCoulombIntra;idx++) {
      ri = CoulombIntra[idx];
      myEnergy += ParaCoulombIntra[idx] * n0[ri] * n1[ri];
    }
  
    /* CoulombInter */
    #pragma omp for private(idx,ri,rj)
    for(idx=0;idx<NCoulombInter;idx++) {
      ri = CoulombInter[idx][0];
      rj = CoulombInter[idx][1];
      myEnergy += ParaCoulombInter[idx] * (n0[ri]+n1[ri]) * (n0[rj]+n1[rj]);
    }

    /* HundCoupling */
    #pragma omp for private(idx,ri,rj)
    for(idx=0;idx<NHundCoupling;idx++) {
      ri = HundCoupling[idx][0];
      rj = HundCoupling[idx][1];
      myEnergy -= ParaHundCoupling[idx] * (n0[ri]*n0[rj] + n1[ri]*n1[rj]);
      /* Caution: negative sign */
    }

    e += myEnergy;
  }

  return e;
}

/* Calculate the transfer terms, */
/* which can be calculated by 1-body Green function. */
/* This function will be used in the Lanczos mode */
double complex CalculateHamiltonian1_fsz(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt,int *eleSpn) {
  double complex e=0.0;
  int idx;
  int ri,rj,s;
  int *myEleIdx, *myEleNum, *myProjCntNew,*myEleSpn;
  double complex *myBuffer;
  double complex myEnergy;

  RequestWorkSpaceThreadInt(Nsize+Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadComplex(NQPFull);
  /* GreenFunc1: NQPFull */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleSpn,myEleNum,myProjCntNew,myBuffer,myEnergy,idx)  \
  reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleSpn = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull);

    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleSpn[idx] = eleSpn[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];

    myEnergy = 0.0;

    /* Transfer */
    #pragma omp for private(idx,ri,rj,s) schedule(dynamic) nowait
    for(idx=0;idx<NTransfer;idx++) {
      ri = Transfer[idx][0];
      s  = Transfer[idx][1];
      rj = Transfer[idx][2];
      
      myEnergy -= ParaTransfer[idx]
        * GreenFunc1_fsz(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      /* Caution: negative sign */
    }

    e += myEnergy;
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  return e;
}

/* Calculate the exchange coupling, pair hopping, interAll terms, */
/* which can be calculated by 2-body Green function. */
/* This function will be used in the Lanczos mode */
double complex CalculateHamiltonian2_fsz(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt,int *eleSpn) {
  double e=0.0, tmp;
  int idx;
  int ri,rj,s,rk,rl,t;
  int *myEleIdx, *myEleNum, *myProjCntNew,*myEleSpn;
  double complex *myBuffer;
  double complex myEnergy;

  RequestWorkSpaceThreadInt(Nsize+Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadComplex(NQPFull+2*Nsize);
  /* GreenFunc2: NQPFull+2*Nsize */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleSpn,myEleNum,myProjCntNew,myBuffer,myEnergy,idx)  \
  reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleSpn = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull+2*Nsize);

    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleSpn[idx] = eleSpn[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];

    myEnergy = 0.0;

    /* Pair Hopping */
    #pragma omp for private(idx,ri,rj) schedule(dynamic) nowait
    for(idx=0;idx<NPairHopping;idx++) {
      ri = PairHopping[idx][0];
      rj = PairHopping[idx][1];
    
      myEnergy += ParaPairHopping[idx]
        * GreenFunc2_fsz(ri,rj,ri,rj,0,1,ip,myEleIdx,eleCfg,myEleNum,
                     eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
    }

    /* Exchange Coupling */
    #pragma omp for private(idx,ri,rj,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NExchangeCoupling;idx++) {
      ri = ExchangeCoupling[idx][0];
      rj = ExchangeCoupling[idx][1];
    
      tmp =  GreenFunc2_fsz(ri,rj,rj,ri,0,1,ip,myEleIdx,eleCfg,myEleNum,
                        eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      tmp += GreenFunc2_fsz(ri,rj,rj,ri,1,0,ip,myEleIdx,eleCfg,myEleNum,
                        eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      myEnergy += ParaExchangeCoupling[idx] * tmp;
    }

    /* Inter All */
    #pragma omp for private(idx,ri,rj,s,rk,rl,t) schedule(dynamic) nowait
    for(idx=0;idx<NInterAll;idx++) {
      ri = InterAll[idx][0];
      rj = InterAll[idx][2];
      s  = InterAll[idx][3];
      rk = InterAll[idx][4];
      rl = InterAll[idx][6];
      t  = InterAll[idx][7];
      
      myEnergy += ParaInterAll[idx]
        * GreenFunc2_fsz(ri,rj,rk,rl,s,t,ip,myEleIdx,eleCfg,myEleNum,
                     eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
    }

    e += myEnergy;
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  return e;
}
