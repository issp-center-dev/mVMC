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
 * make sample
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "global.h"
#include "vmcmake.h"
#include "slater.h"
#ifndef _SRC_VMCMAKE
#define _SRC_VMCMAKE
#include "projection.h"
#include "pfupdate_two_fcmp.h"
#include "matrix.c"
#include "splitloop.h"
#include "qp.h"

void VMCMakeSample(MPI_Comm comm) {
  int outStep,nOutStep;
  int inStep,nInStep;
  UpdateType updateType;
  int mi,mj,ri,rj,s,t,i;
  int nAccept=0;
  int sample;

  double complex logIpOld,logIpNew; /* logarithm of inner product <phi|L|x> */ // is this ok ? TBC
  int projCntNew[NProj];
  double complex pfMNew[NQPFull];
  double x,w; // TBC x will be complex number

  int qpStart,qpEnd;
  int rejectFlag;
  int rank,size;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  SplitLoop(&qpStart,&qpEnd,NQPFull,rank,size);

  StartTimer(30);
  if(BurnFlag==0) {
    makeInitialSample(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,
                      qpStart,qpEnd,comm);
  } else {
    copyFromBurnSample(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt);
  }
  
#ifdef _pf_block_update
  // TODO: Compute from qpStart to qpEnd to support loop splitting.
  void *pfOrbital[NQPFull];
  void *pfUpdator[NQPFull];
  void *pfMat[NQPFull];
  void *pfMap[NQPFull];
  // Read block size from input.
  const char *optBlockSize = getenv("VMC_BLOCK_UPDATE_SIZE");
  if (optBlockSize)
    NBlockUpdateSize = atoi(optBlockSize);
  // Fall back to default if input is invalid.
  if (NBlockUpdateSize < 1 || NBlockUpdateSize > 100)
    if (NExUpdatePath == 0)
      NBlockUpdateSize = 4;
    else
      NBlockUpdateSize = 20;

  // Set one universal EleSpn.
  for (mi=0; mi<Ne;  mi++) EleSpn[mi] = 0;
  for (mi=Ne;mi<Ne*2;mi++) EleSpn[mi] = 1;
  // Initialize.
  updated_tdi_v_init_z(NQPFull, Nsite, Nsite2, Nsize,
                       SlaterElm, Nsite2*Nsite2,
                       InvM, Nsize*Nsize,
                       TmpEleIdx, EleSpn,
                       NBlockUpdateSize,
                       pfUpdator, pfOrbital, pfMat, pfMap);
  updated_tdi_v_get_pfa_z(NQPFull, PfM, pfUpdator);
#else
  CalculateMAll_fcmp(TmpEleIdx,qpStart,qpEnd);
#endif
  // printf("DEBUG: maker1: PfM=%lf\n",creal(PfM[0]));
  logIpOld = CalculateLogIP_fcmp(PfM,qpStart,qpEnd,comm);

  if( !isfinite(creal(logIpOld) + cimag(logIpOld)) ) {
    if(rank==0) fprintf(stderr,"waring: VMCMakeSample remakeSample logIpOld=%e\n",creal(logIpOld)); //TBC
    makeInitialSample(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,
                      qpStart,qpEnd,comm);
#ifdef _pf_block_update
    // Clear and reinitialize.
    updated_tdi_v_free_z(NQPFull, pfUpdator, pfOrbital, pfMat, pfMap);
    updated_tdi_v_init_z(NQPFull, Nsite, Nsite2, Nsize,
                         SlaterElm, Nsite2*Nsite2,
                         InvM, Nsize*Nsize,
                         TmpEleIdx, EleSpn,
                         NBlockUpdateSize,
                         pfUpdator, pfOrbital, pfMat, pfMap);
    updated_tdi_v_get_pfa_z(NQPFull, PfM, pfUpdator);
#else
    CalculateMAll_fcmp(TmpEleIdx,qpStart,qpEnd);
#endif
    //printf("DEBUG: maker2: PfM=%lf\n",creal(PfM[0]));
    logIpOld = CalculateLogIP_fcmp(PfM,qpStart,qpEnd,comm);
    BurnFlag = 0;
  }
  StopTimer(30);

  nOutStep = (BurnFlag==0) ? NVMCWarmUp+NVMCSample : NVMCSample+1;
  nInStep = NVMCInterval * Nsite;

  for(i=0;i<Counter_max;i++) Counter[i]=0;  /* reset counter */

  for(outStep=0;outStep<nOutStep;outStep++) {
    for(inStep=0;inStep<nInStep;inStep++) {

      updateType = getUpdateType(NExUpdatePath);

      if(updateType==HOPPING) { /* hopping */
        Counter[0]++;

        StartTimer(31);
        makeCandidate_hopping(&mi, &ri, &rj, &s, &rejectFlag,
                              TmpEleIdx, TmpEleCfg);
        StopTimer(31);

        if(rejectFlag) continue;

        StartTimer(32);
        StartTimer(60);
        /* The mi-th electron with spin s hops to site rj */
        updateEleConfig(mi,ri,rj,s,TmpEleIdx,TmpEleCfg,TmpEleNum);
        UpdateProjCnt(ri,rj,s,projCntNew,TmpEleProjCnt,TmpEleNum);
        StopTimer(60);

        StartTimer(61);
#ifdef _pf_block_update
        updated_tdi_v_push_z(NQPFull, rj+s*Nsite, mi+s*Ne, 1, pfUpdator);
        updated_tdi_v_get_pfa_z(NQPFull, pfMNew, pfUpdator);
#else
        //CalculateNewPfM2(mi,s,pfMNew,TmpEleIdx,qpStart,qpEnd);
        CalculateNewPfM2(mi,s,pfMNew,TmpEleIdx,qpStart,qpEnd);
#endif
        //printf("DEBUG: out %d in %d pfMNew=%lf \n",outStep,inStep,creal(pfMNew[0]));
        StopTimer(61);

        StartTimer(62);
        /* calculate inner product <phi|L|x> */
        //logIpNew = CalculateLogIP_fcmp(pfMNew,qpStart,qpEnd,comm);
        logIpNew = CalculateLogIP_fcmp(pfMNew,qpStart,qpEnd,comm);
        StopTimer(62);

        /* Metroplis */
        x = LogProjRatio(projCntNew,TmpEleProjCnt);
        w = exp(2.0*(x+creal(logIpNew-logIpOld)));
        if( !isfinite(w) ) w = -1.0; /* should be rejected */

        if(w > genrand_real2()) { /* accept */
          StartTimer(63);
#ifdef _pf_block_update
          // Inv already updated. Only need to get PfM again.
          updated_tdi_v_get_pfa_z(NQPFull, PfM, pfUpdator);
#else
          // UpdateMAll will change SlaterElm, InvM (including PfM)
          UpdateMAll(mi,s,TmpEleIdx,qpStart,qpEnd);
#endif
          StopTimer(63);

          for(i=0;i<NProj;i++) TmpEleProjCnt[i] = projCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          Counter[1]++;
        } else { /* reject */
#ifdef _pf_block_update
          updated_tdi_v_pop_z(NQPFull, 0, pfUpdator);
#endif
          revertEleConfig(mi,ri,rj,s,TmpEleIdx,TmpEleCfg,TmpEleNum);
        }
        StopTimer(32);

      } else if(updateType==EXCHANGE) { /* exchange */
        Counter[2]++;

        StartTimer(31);
        makeCandidate_exchange(&mi, &ri, &rj, &s, &rejectFlag,
                               TmpEleIdx, TmpEleCfg, TmpEleNum);
        StopTimer(31);

        if(rejectFlag) continue;

        StartTimer(33);
        StartTimer(65);

        /* The mi-th electron with spin s exchanges with the electron on site rj with spin 1-s */
        t = 1-s;
        mj = TmpEleCfg[rj+t*Nsite];

        /* The mi-th electron with spin s hops to rj */
        updateEleConfig(mi,ri,rj,s,TmpEleIdx,TmpEleCfg,TmpEleNum);
        UpdateProjCnt(ri,rj,s,projCntNew,TmpEleProjCnt,TmpEleNum);
        /* The mj-th electron with spin t hops to ri */
        updateEleConfig(mj,rj,ri,t,TmpEleIdx,TmpEleCfg,TmpEleNum);
        UpdateProjCnt(rj,ri,t,projCntNew,projCntNew,TmpEleNum);

        StopTimer(65);
        StartTimer(66);

#ifdef _pf_block_update
        updated_tdi_v_push_pair_z(NQPFull,
                                  rj+s*Nsite, mi+s*Ne,
                                  ri+t*Nsite, mj+t*Ne,
                                  1, pfUpdator);
        updated_tdi_v_get_pfa_z(NQPFull, pfMNew, pfUpdator);
#else
        CalculateNewPfMTwo2_fcmp(mi, s, mj, t, pfMNew, TmpEleIdx, qpStart, qpEnd);
#endif
        StopTimer(66);
        StartTimer(67);

        /* calculate inner product <phi|L|x> */
        logIpNew = CalculateLogIP_fcmp(pfMNew,qpStart,qpEnd,comm);

        StopTimer(67);

        /* Metroplis */
        x = LogProjRatio(projCntNew,TmpEleProjCnt);
        w = exp(2.0*(x+creal(logIpNew-logIpOld))); //TBC
        if( !isfinite(w) ) w = -1.0; /* should be rejected */

        if(w > genrand_real2()) { /* accept */
          StartTimer(68);
#ifdef _pf_block_update
          // Inv already updated. Only need to get PfM again.
          updated_tdi_v_get_pfa_z(NQPFull, PfM, pfUpdator);
#else
          UpdateMAllTwo_fcmp(mi, s, mj, t, ri, rj, TmpEleIdx,qpStart,qpEnd);
#endif
          StopTimer(68);

          for(i=0;i<NProj;i++) TmpEleProjCnt[i] = projCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          Counter[3]++;
        } else { /* reject */
#ifdef _pf_block_update
          updated_tdi_v_pop_z(NQPFull, 0, pfUpdator);
          updated_tdi_v_pop_z(NQPFull, 0, pfUpdator);
#endif
          revertEleConfig(mj,rj,ri,t,TmpEleIdx,TmpEleCfg,TmpEleNum);
          revertEleConfig(mi,ri,rj,s,TmpEleIdx,TmpEleCfg,TmpEleNum);
        }
        StopTimer(33);
      }

      if(nAccept>Nsite) {
        // Recalculate PfM and InvM.
        StartTimer(34);
#ifdef _pf_block_update
        // Clear and reinitialize.
        updated_tdi_v_free_z(NQPFull, pfUpdator, pfOrbital, pfMat, pfMap);
        updated_tdi_v_init_z(NQPFull, Nsite, Nsite2, Nsize,
                             SlaterElm, Nsite2*Nsite2,
                             InvM, Nsize*Nsize,
                             TmpEleIdx, EleSpn,
                             NBlockUpdateSize,
                             pfUpdator, pfOrbital, pfMat, pfMap);
        updated_tdi_v_get_pfa_z(NQPFull, PfM, pfUpdator);
#else
        CalculateMAll_fcmp(TmpEleIdx,qpStart,qpEnd);
#endif
        //printf("DEBUG: maker3: PfM=%lf\n",creal(PfM[0]));
        logIpOld = CalculateLogIP_fcmp(PfM,qpStart,qpEnd,comm);
        StopTimer(34);
        nAccept=0;
      }
    } /* end of instep */

    StartTimer(35);
    /* save Electron Configuration */
    if(outStep >= nOutStep-NVMCSample) {
      sample = outStep-(nOutStep-NVMCSample);
      saveEleConfig(sample,logIpOld,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt);
    }
    StopTimer(35);

  } /* end of outstep */

  copyToBurnSample(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt);
  BurnFlag=1;

#ifdef _pf_block_update
  // Free-up updator space.
  updated_tdi_v_free_z(NQPFull, pfUpdator, pfOrbital, pfMat, pfMap);
#endif

  return;
}

int makeInitialSample(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,
                      const int qpStart, const int qpEnd, MPI_Comm comm) {
  const int nsize = Nsize;
  const int nsite2 = Nsite2;
  int flag=1,flagRdc,loop=0;
  int ri,mi,si,msi,rsi;
  int rank,size;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  
  do {
    /* initialize */
    #pragma omp parallel for default(shared) private(msi)
    for(msi=0;msi<nsize;msi++) eleIdx[msi] = -1;
    #pragma omp parallel for default(shared) private(rsi)
    for(rsi=0;rsi<nsite2;rsi++) eleCfg[rsi] = -1;
    
    /* local spin */
    for(ri=0;ri<Nsite;ri++) {
      if(LocSpn[ri]==1) {
        do {
          mi = gen_rand32()%Ne;
          si = (genrand_real2()<0.5) ? 0 : 1;
        } while(eleIdx[mi+si*Ne]!=-1);
        eleCfg[ri+si*Nsite] = mi;
        eleIdx[mi+si*Ne] = ri;
      }
    }
    
    /* itinerant electron */
    for(si=0;si<2;si++) {
      for(mi=0;mi<Ne;mi++) {
        if(eleIdx[mi+si*Ne]== -1) {
          do {
            ri = gen_rand32()%Nsite;
          } while (eleCfg[ri+si*Nsite]!= -1 || LocSpn[ri]==1);
          eleCfg[ri+si*Nsite] = mi;
          eleIdx[mi+si*Ne] = ri;
        }
      }
    }
    
    /* EleNum */
    #pragma omp parallel for default(shared) private(rsi)
    #pragma loop noalias
    for(rsi=0;rsi<nsite2;rsi++) {
      eleNum[rsi] = (eleCfg[rsi] < 0) ? 0 : 1;
    }
    
    MakeProjCnt(eleProjCnt,eleNum);

    flag = CalculateMAll_fcmp(eleIdx,qpStart,qpEnd);
    //printf("DEBUG: maker4: PfM=%lf\n",creal(PfM[0]));
    if(size>1) {
      MPI_Allreduce(&flag,&flagRdc,1,MPI_INT,MPI_MAX,comm);
      flag = flagRdc;
    }

    loop++;
    if(loop>100) {
      if(rank==0) fprintf(stderr, "error: makeInitialSample: Too many loops\n");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
  } while (flag>0);
 
  return 0;
}

void copyFromBurnSample(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt) {
  int i,n;
  const int *burnEleIdx = BurnEleIdx;
  n = Nsize + 2*Nsite + 2*Nsite + NProj;
  #pragma loop noalias
  for(i=0;i<n;i++) eleIdx[i] = burnEleIdx[i];
  return;
}

void copyToBurnSample(const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt) {
  int i,n;
  int *burnEleIdx = BurnEleIdx;
  n = Nsize + 2*Nsite + 2*Nsite + NProj;
  #pragma loop noalias
  for(i=0;i<n;i++) burnEleIdx[i] = eleIdx[i];
  return;
}

void saveEleConfig(const int sample, const double complex logIp,
                   const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt) {
  int i,offset;
  double x;
  const int nsize=Nsize;
  const int nsite2 = Nsite2;
  const int nProj = NProj;

  offset = sample*nsize;
  #pragma loop noalias
  for(i=0;i<nsize;i++) EleIdx[offset+i] = eleIdx[i];
  offset = sample*nsite2;
  #pragma loop noalias
  for(i=0;i<nsite2;i++) EleCfg[offset+i] = eleCfg[i];
  #pragma loop noalias
  for(i=0;i<nsite2;i++) EleNum[offset+i] = eleNum[i];
  offset = sample*nProj;
  #pragma loop noalias
  for(i=0;i<nProj;i++) EleProjCnt[offset+i] = eleProjCnt[i];
  
  x = LogProjVal(eleProjCnt);
  logSqPfFullSlater[sample] = 2.0*(x+creal(logIp));//TBC
  
  return;
}

void sortEleConfig(int *eleIdx, int *eleCfg, const int *eleNum) {
/*   int ri,mi=0; */
/*   for(ri=0;ri<Nsite;ri++) { */
/*     if(eleNum[ri]>0) { */
/*       eleCfg[ri]=mi; */
/*       eleIdx[mi]=ri; */
/*       mi++; */
/*     } else { */
/*       eleCfg[ri]=-1; */
/*     } */
/*   } */
/*   mi=0; */
/*   for(ri=0;ri<Nsite;ri++) { */
/*     if(eleNum[ri+Nsite]>0) { */
/*       eleCfg[ri+Nsite]=mi; */
/*       eleIdx[mi+Ne]=ri; */
/*       mi++; */
/*     } else { */
/*       eleCfg[ri+Nsite]=-1; */
/*     } */
/*   } */

  return;
}

void ReduceCounter(MPI_Comm comm) {
  #ifdef _mpi_use
  int n=Counter_max;
  int recv[n];
  int i;
  int rank,size;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  MPI_Allreduce(Counter,recv,n,MPI_INT,MPI_SUM,comm);
  if(rank==0) {
    for(i=0;i<n;i++) Counter[i] = recv[i];
  }
  #endif
  return;
}


/* The mi-th electron with spin s hops to site rj */
void makeCandidate_hopping(int *mi_, int *ri_, int *rj_, int *s_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg) {
  const int icnt_max = Nsite*Nsite;
  int icnt;
  int mi, ri, rj, s, flag;

  flag = 0; // FALSE
  do {
    mi = gen_rand32()%Ne;
    s = (genrand_real2()<0.5) ? 0 : 1;
    ri = eleIdx[mi+s*Ne];
  } while (LocSpn[ri] == 1);

  icnt = 0;
  do {
    rj = gen_rand32()%Nsite;
    if(icnt> icnt_max){
      flag = 1; // TRUE
      break;
    }
    icnt+=1;
  } while (eleCfg[rj+s*Nsite] != -1 || LocSpn[rj]==1);

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_ = s;
  *rejectFlag_ = flag;

  return;
}

/* The mi-th electron with spin s exchanges with the electron on site rj with spin 1-s */
void makeCandidate_exchange(int *mi_, int *ri_, int *rj_, int *s_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg, const int *eleNum) {
  int mi, mj, ri, rj, s, t, flag;

  flag = 1; // TRUE
  for(ri=0;ri<Nsite;ri++){
    if((eleNum[ri]+eleNum[ri+Nsite]) == 1){
      flag = 0; // FALSE
      break;
    }
  }
  if(flag) {
    *rejectFlag_ = flag;
    return;
  }

  do {
    mi = gen_rand32()%Ne;
    s = (genrand_real2()<0.5) ? 0 : 1;
    ri = eleIdx[mi+s*Ne];
  } while (eleCfg[ri+(1-s)*Nsite] != -1);
  do {
    mj = gen_rand32()%Ne;
    t = 1-s;
    rj = eleIdx[mj+t*Ne];
  } while (eleCfg[rj+(1-t)*Nsite] != -1);

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_ = s;
  *rejectFlag_ = flag;

  return;
}

/* The mi-th electron with spin s hops to site rj */
void updateEleConfig(int mi, int ri, int rj, int s,
                     int *eleIdx, int *eleCfg, int *eleNum) {
  eleIdx[mi+s*Ne] = rj;
  eleCfg[ri+s*Nsite] = -1;
  eleCfg[rj+s*Nsite] = mi;
  eleNum[ri+s*Nsite] = 0;
  eleNum[rj+s*Nsite] = 1;
  return;
}

void revertEleConfig(int mi, int ri, int rj, int s,
                     int *eleIdx, int *eleCfg, int *eleNum) {
  eleIdx[mi+s*Ne] = ri;
  eleCfg[ri+s*Nsite] = mi;
  eleCfg[rj+s*Nsite] = -1;
  eleNum[ri+s*Nsite] = 1;
  eleNum[rj+s*Nsite] = 0;
  return;
}


UpdateType getUpdateType(int path) {
  if(path==0) {
    return HOPPING;
  } else if (path==1) {
    return (genrand_real2()<0.5) ? EXCHANGE : HOPPING; /* exchange or hopping */
  } else if (path==2) {
    if(iFlgOrbitalGeneral==0){
      return EXCHANGE;
    }else{
      if(TwoSz==-1){ //Sz is not conserved
        return (genrand_real2()<0.5) ? EXCHANGE : LOCALSPINFLIP; /* exchange or localspinflip */ //fsz
      }else{
        return EXCHANGE ; /* exchange */
      } 
    }
  }else if(path==3){ //for KondoGC
    if(genrand_real2()<0.5){ // for conduction electrons
      return HOPPING; /* hopping */
    }else{/* Exchange for conductions and local spins, localspinflip for local spins　*/
      return (genrand_real2()<0.5) ? EXCHANGE : LOCALSPINFLIP; /* exchange or localspinflip */
    }
  }
  return NONE;
}

void VMC_BF_MakeSample(MPI_Comm comm)
{
  int outStep, nOutStep;
  int inStep, nInStep;
  UpdateType updateType;
  int mi, mj, ri, rj, s, t, i;
  int nAccept = 0;
  int sample;

  double complex logIpOld, logIpNew; /* logarithm of inner product <phi|L|x> */ // is this ok ? TBC
  int projCntNew[NProj];
  int projBFCntNew[16 * Nsite * Nrange]; // For BackFlow
  int msaTmp[NQPFull * Nsite], icount[NQPFull]; // For BackFlow
  double complex pfMNew[NQPFull];
  double x, w; // TBC x will be complex number

  int qpStart, qpEnd;
  int rejectFlag;
  int rank, size;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  SplitLoop(&qpStart, &qpEnd, NQPFull, rank, size);

  StartTimer(30);
  if (BurnFlag == 0) {
    makeInitialSampleBF(TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt, TmpEleProjBFCnt,
                        qpStart, qpEnd, comm);
  } else {
    copyFromBurnSampleBF(TmpEleIdx);
    MakeSlaterElmBF_fcmp(TmpEleNum, TmpEleProjBFCnt);
  }

  CalculateMAll_BF_fcmp(TmpEleIdx, qpStart, qpEnd);
  logIpOld = CalculateLogIP_fcmp(PfM, qpStart, qpEnd, comm);
  if (! (isfinite(creal(logIpOld)) && isfinite(cimag(logIpOld)))) {
    if (rank == 0) fprintf(stderr, "waring: VMCMakeSample remakeSample logIpOld=%e\n", creal(logIpOld)); //TBC
    makeInitialSampleBF(TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt, TmpEleProjBFCnt,
                        qpStart, qpEnd, comm);
    CalculateMAll_BF_fcmp(TmpEleIdx, qpStart, qpEnd);
    logIpOld = CalculateLogIP_fcmp(PfM, qpStart, qpEnd, comm);
    BurnFlag = 0;
  }
  StopTimer(30);

  nOutStep = (BurnFlag == 0) ? NVMCWarmUp + NVMCSample : NVMCSample + 1;
  nInStep = NVMCInterval * Nsite;

  for (i = 0; i < Counter_max; i++) Counter[i] = 0;  /* reset counter */

  for (outStep = 0; outStep < nOutStep; outStep++) {
    for (inStep = 0; inStep < nInStep; inStep++) {

      updateType = getUpdateType(NExUpdatePath);

      if (updateType == HOPPING) { /* hopping */
        Counter[0]++;

        StartTimer(31);
        makeCandidate_hopping(&mi, &ri, &rj, &s, &rejectFlag,
                              TmpEleIdx, TmpEleCfg);
        StopTimer(31);

        if (rejectFlag) continue;

        StartTimer(32);
        StartTimer(60);
        /* The mi-th electron with spin s hops to site rj */
        updateEleConfig(mi, ri, rj, s, TmpEleIdx, TmpEleCfg, TmpEleNum);
        UpdateProjCnt(ri, rj, s, projCntNew, TmpEleProjCnt, TmpEleNum);
        MakeProjBFCnt(projBFCntNew, TmpEleNum);
        StopTimer(60);
        UpdateSlaterElmBF_fcmp(mi, ri, rj, s, TmpEleCfg, TmpEleNum, projBFCntNew, msaTmp, icount,
                               SlaterElmBF);
        StartTimer(61);
        //CalculateNewPfM2(mi,s,pfMNew,TmpEleIdx,qpStart,qpEnd);
        //CalculateNewPfM2_real(mi,s,pfMNew_real,TmpEleIdx,qpStart,qpEnd);
        CalculateNewPfMBF(icount, msaTmp, pfMNew, TmpEleIdx, qpStart, qpEnd, SlaterElmBF);

        //printf("DEBUG: out %d in %d pfMNew=%lf \n",outStep,inStep,creal(pfMNew[0]));
        StopTimer(61);

        StartTimer(62);
        /* calculate inner product <phi|L|x> */
        //logIpNew = CalculateLogIP_fcmp(pfMNew,qpStart,qpEnd,comm);
        logIpNew = CalculateLogIP_fcmp(pfMNew, qpStart, qpEnd, comm);
        StopTimer(62);

        /* Metroplis */
        x = LogProjRatio(projCntNew, TmpEleProjCnt);
        w = exp(2.0 * (x + (logIpNew - logIpOld)));
        if (!isfinite(w)) w = -1.0; /* should be rejected */

        if (w > genrand_real2()) { /* accept */
          // UpdateMAll will change SlaterElm, InvM (including PfM)
          StartTimer(63);
          UpdateMAll_BF_fcmp(icount, msaTmp, PfM, TmpEleIdx, qpStart, qpEnd);
          //          UpdateMAll_real(mi,s,TmpEleIdx,qpStart,qpEnd);
          //            UpdateMAll(mi,s,TmpEleIdx,qpStart,qpEnd);
          StopTimer(63);

          for (i = 0; i < NProj; i++) TmpEleProjCnt[i] = projCntNew[i];
          for (i = 0; i < 16 * Nsite * Nrange; i++) TmpEleProjBFCnt[i] = projBFCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          Counter[1]++;
        } else { /* reject */
          revertEleConfig(mi, ri, rj, s, TmpEleIdx, TmpEleCfg, TmpEleNum);
          //TODO: Add Timer
          UpdateSlaterElmBF_fcmp(mi, rj, ri, s, TmpEleCfg, TmpEleNum, TmpEleProjBFCnt, msaTmp, icount,
                                 SlaterElmBF);
        }
        StopTimer(32);

      } else if (updateType == EXCHANGE) { /* exchange */
        Counter[2]++;

        StartTimer(31);
        makeCandidate_exchange(&mi, &ri, &rj, &s, &rejectFlag,
                               TmpEleIdx, TmpEleCfg, TmpEleNum);
        StopTimer(31);

        if (rejectFlag) continue;

        StartTimer(33);
        StartTimer(65);

        /* The mi-th electron with spin s exchanges with the electron on site rj with spin 1-s */
        t = 1 - s;
        mj = TmpEleCfg[rj + t * Nsite];

        /* The mi-th electron with spin s hops to rj */
        updateEleConfig(mi, ri, rj, s, TmpEleIdx, TmpEleCfg, TmpEleNum);
        UpdateProjCnt(ri, rj, s, projCntNew, TmpEleProjCnt, TmpEleNum);
        /* The mj-th electron with spin t hops to ri */
        updateEleConfig(mj, rj, ri, t, TmpEleIdx, TmpEleCfg, TmpEleNum);
        UpdateProjCnt(rj, ri, t, projCntNew, projCntNew, TmpEleNum);

        StopTimer(65);
        StartTimer(66);

        CalculateNewPfMTwo2_fcmp(mi, s, mj, t, pfMNew, TmpEleIdx, qpStart, qpEnd);
        StopTimer(66);
        StartTimer(67);

        /* calculate inner product <phi|L|x> */
        logIpNew = CalculateLogIP_fcmp(pfMNew, qpStart, qpEnd, comm);

        StopTimer(67);

        /* Metroplis */
        x = LogProjRatio(projCntNew, TmpEleProjCnt);
        w = exp(2.0 * (x + (logIpNew - logIpOld))); //TBC
        if (!isfinite(w)) w = -1.0; /* should be rejected */

        if (w > genrand_real2()) { /* accept */
          StartTimer(68);
          UpdateMAllTwo_fcmp(mi, s, mj, t, ri, rj, TmpEleIdx, qpStart, qpEnd);
          StopTimer(68);

          for (i = 0; i < NProj; i++) TmpEleProjCnt[i] = projCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          Counter[3]++;
        } else { /* reject */
          revertEleConfig(mj, rj, ri, t, TmpEleIdx, TmpEleCfg, TmpEleNum);
          revertEleConfig(mi, ri, rj, s, TmpEleIdx, TmpEleCfg, TmpEleNum);
        }
        StopTimer(33);
      }

      if (nAccept > Nsite) {
        StartTimer(34);
        /* recal PfM and InvM */
        //CalculateMAll_real(TmpEleIdx,qpStart,qpEnd);
        //printf("DEBUG: maker3: PfM=%lf\n",creal(PfM[0]));
        CalculateMAll_BF_fcmp(TmpEleIdx, qpStart, qpEnd);
        logIpOld = CalculateLogIP_fcmp(PfM, qpStart, qpEnd, comm);
        StopTimer(34);
        nAccept = 0;
      }
    } /* end of instep */

    StartTimer(35);
    /* save Electron Configuration */
    if (outStep >= nOutStep - NVMCSample) {
      sample = outStep - (nOutStep - NVMCSample);
      saveEleConfigBF(sample, logIpOld, TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt, TmpEleProjBFCnt);
    }
    StopTimer(35);
  } /* end of outstep */

  copyToBurnSampleBF(TmpEleIdx);
  BurnFlag = 1;
  return;
}

int makeInitialSampleBF(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt, int *eleProjBFCnt,
                        const int qpStart, const int qpEnd, MPI_Comm comm)
{
  const int nsize = Nsize;
  const int nsite2 = Nsite2;
  int flag = 1, flagRdc, loop = 0;
  int ri, mi, si, msi, rsi;
  int rank, size;

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  do {
    /* initialize */
#pragma omp parallel for default(shared) private(msi)
    for (msi = 0; msi < nsize; msi++) eleIdx[msi] = -1;
#pragma omp parallel for default(shared) private(rsi)
    for (rsi = 0; rsi < nsite2; rsi++) eleCfg[rsi] = -1;

    /* local spin */
    for (ri = 0; ri < Nsite; ri++) {
      if (LocSpn[ri] == 1) {
        do {
          mi = gen_rand32() % Ne;
          si = (genrand_real2() < 0.5) ? 0 : 1;
        } while (eleIdx[mi + si * Ne] != -1);
        eleCfg[ri + si * Nsite] = mi;
        eleIdx[mi + si * Ne] = ri;
      }
    }

    /* itinerant electron */
    for (si = 0; si < 2; si++) {
      for (mi = 0; mi < Ne; mi++) {
        if (eleIdx[mi + si * Ne] == -1) {
          do {
            ri = gen_rand32() % Nsite;
          } while (eleCfg[ri + si * Nsite] != -1 || LocSpn[ri] == 1);
          eleCfg[ri + si * Nsite] = mi;
          eleIdx[mi + si * Ne] = ri;
        }
      }
    }

    /* EleNum */
#pragma omp parallel for default(shared) private(rsi)
#pragma loop noalias
    for (rsi = 0; rsi < nsite2; rsi++) {
      eleNum[rsi] = (eleCfg[rsi] < 0) ? 0 : 1;
    }

    MakeProjCnt(eleProjCnt, eleNum);
    MakeProjBFCnt(eleProjBFCnt, eleNum);

    MakeSlaterElmBF_fcmp(eleNum, eleProjBFCnt);

    flag = CalculateMAll_BF_fcmp(eleIdx, qpStart, qpEnd);
    //printf("DEBUG: maker4: PfM=%lf\n",creal(PfM[0]));
    if (size > 1) {
      MPI_Allreduce(&flag, &flagRdc, 1, MPI_INT, MPI_MAX, comm);
      flag = flagRdc;
    }

    loop++;
    if (loop > 100) {
      if (rank == 0) fprintf(stderr, "error: makeInitialSample: Too many loops\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  } while (flag > 0);
  return 0;
}

void copyFromBurnSampleBF(int *eleIdx) {
  int i, n;
  const int *burnEleIdx = BurnEleIdx;
  n = Nsize + 2 * Nsite + 2 * Nsite + NProj + 16 * Nsite * Nrange;
#pragma loop noalias
  for (i = 0; i < n; i++) eleIdx[i] = burnEleIdx[i];
  return;
}

void copyToBurnSampleBF(const int *eleIdx) {
  int i, n;
  int *burnEleIdx = BurnEleIdx;
  n = Nsize + 2 * Nsite + 2 * Nsite + NProj + 16 * Nsite * Nrange;
#pragma loop noalias
  for (i = 0; i < n; i++) burnEleIdx[i] = eleIdx[i];
  return;
}

void saveEleConfigBF(const int sample, const double logIp,
                     const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt,
                     const int *eleProjBFCnt) {
  int i, offset;
  double x;
  const int nsize = Nsize;
  const int nsite2 = Nsite2;
  const int nProj = NProj;
  //const int nQPFull = NQPFull;

  offset = sample * nsize;
#pragma loop noalias
  for (i = 0; i < nsize; i++) EleIdx[offset + i] = eleIdx[i];
  offset = sample * nsite2;
#pragma loop noalias
  for (i = 0; i < nsite2; i++) EleCfg[offset + i] = eleCfg[i];
#pragma loop noalias
  for (i = 0; i < nsite2; i++) EleNum[offset + i] = eleNum[i];
  offset = sample * nProj;
#pragma loop noalias
  for (i = 0; i < nProj; i++) EleProjCnt[offset + i] = eleProjCnt[i];
  offset = sample * 16 * Nsite * Nrange;
#pragma loop noalias
  for (i = 0; i < 16 * Nsite * Nrange; i++) EleProjBFCnt[offset + i] = eleProjBFCnt[i];

  x = LogProjVal(eleProjCnt);
  logSqPfFullSlater[sample] = 2.0 * (x + logIp);

  /*
  if (NStoreM != 0) {
    offset = sample * nQPFull * nsize * nsize;
#pragma loop noalias
    for (i = 0; i < nQPFull * nsize * nsize; i++) InvM_Store[offset + i] = InvM[i];
    offset = sample * nQPFull;
#pragma loop noalias
    for (i = 0; i < nQPFull; i++) PfM_Store[offset + i] = PfM[i];
    offset = sample * nsite2 * nsite2 * NQPFull;
#pragma loop noalias
    for (i = 0; i < nsite2 * nsite2 * nQPFull; i++) {
      SmpSltElmBF_real[offset + i] = SlaterElmBF_real[i];
    }
    offset = sample * NQPFull * nsite * nsite;
#pragma loop noalias
    for (i = 0; i < nsite; i++) {
      for (j = 0; j < nsite; j++) {
        SmpEta[offset + i * nsite + j] = eta[i][j];
        SmpEtaFlag[offset + i * nsite + j] = etaFlag[i][j];
      }
    }
  }
  */
  return;
}
#endif
