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
 * make sample
 *-------------------------------------------------------------
 *-------------------------------------------------------------*/
#include "vmcmake_fsz_real.h"
#ifndef _SRC_VMCMAKE_FSZ_REAL
#define _SRC_VMCMAKE_FSZ_REAL

#include "global.h"
#include "slater.h"
#include "matrix.h"
#include "pfupdate_fsz_real.h"
#include "qp_real.h"
#include "splitloop.h"
#include "vmcmake_fsz.h"

void VMCMakeSample_fsz_real(MPI_Comm comm) {
  int outStep,nOutStep;
  int inStep,nInStep;
  UpdateType updateType;
  int mi,mj,ri,rj,s,t,i;
  int nAccept=0;
  int sample;

  double logIpOld,logIpNew; /* logarithm of inner product <phi|L|x> */ // is this ok ? TBC
  int projCntNew[NProj];
  double pfMNew[NQPFull];
  double x,w; // TBC x will be complex number

  int qpStart,qpEnd;
  int rejectFlag;
  int rank,size;
  int flag_hop;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  SplitLoop(&qpStart,&qpEnd,NQPFull,rank,size);

  StartTimer(30);
  if(BurnFlag==0) {
#ifdef _DEBUG_DETAIL
    printf("DEBUG: make1: \n");
#endif
    makeInitialSample_fsz_real(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,TmpEleSpn,
                      qpStart,qpEnd,comm);
//DEBUG
    //int total_num;
    //CheckEleConfig_fsz(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn,comm);
    //total_num= CheckEleNum_fsz(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn,comm);
    //printf("%d \n",total_num);
//DEBUG
  } else {
    copyFromBurnSample_fsz(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,TmpEleSpn) ;//fsz
  }
  
  CalculateMAll_fsz_real(TmpEleIdx,TmpEleSpn,qpStart,qpEnd);
#ifdef _DEBUG_DETAIL
  printf("DEBUG: maker1: PfM=%lf\n", PfM_real[0]);
#endif
  logIpOld = CalculateLogIP_real(PfM_real,qpStart,qpEnd,comm);
  if( !isfinite(logIpOld) ) {
    if(rank==0) fprintf(stderr,"waring: VMCMakeSample remakeSample logIpOld=%e\n",logIpOld); //TBC
    makeInitialSample_fsz_real(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,TmpEleSpn,
                      qpStart,qpEnd,comm);
    CalculateMAll_fsz_real(TmpEleIdx,TmpEleSpn,qpStart,qpEnd);
#ifdef _DEBUG_DETAIL
    printf("DEBUG: maker2: PfM=%lf\n",creal(PfM_real[0]));
#endif
    logIpOld = CalculateLogIP_real(PfM_real,qpStart,qpEnd,comm);
    BurnFlag = 0;
  }
  StopTimer(30);

  nOutStep = (BurnFlag==0) ? NVMCWarmUp+NVMCSample : NVMCSample+1;
  nInStep = NVMCInterval * Nsite;

  for(i=0;i<Counter_max;i++) Counter[i]=0;  /* reset counter */
  // Counter[0] ->  Hopping  all,      Counter[1] ->  Hopping accept
  // Counter[2] ->  exchange all,      Counter[3] ->  exchange accept
  // Counter[4] ->  localspinflip all, Counter[5] ->  localspin flip accept

  for(outStep=0;outStep<nOutStep;outStep++) {
    for(inStep=0;inStep<nInStep;inStep++) {
#ifdef _DEBUG_DETAIL
      fprintf(stdout, "instep=%d/%d, outstep=%d/%d\n", inStep,nInStep, outStep, nOutStep);
#endif
//DEBUG
      //CheckEleConfig_fsz(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn,comm);
      //total_num= CheckEleNum_fsz(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn,comm);
      //printf("%d \n",total_num);
//DEBUG
      updateType = getUpdateType(NExUpdatePath);

      if(updateType==HOPPING) { /* hopping */
        
        StartTimer(31);
        flag_hop = 0;
        if(TwoSz==-1){//total spin is not conserved
          if(genrand_real2()<0.5){ // this ratio can be changed
            flag_hop = 1;
            Counter[0]++;
            makeCandidate_hopping_fsz(&mi, &ri, &rj, &s,&t, &rejectFlag,
                              TmpEleIdx, TmpEleCfg,TmpEleNum,TmpEleSpn);
          }else{
            Counter[4]++;
            makeCandidate_LocalSpinFlip_conduction(&mi, &ri, &rj, &s,&t, &rejectFlag,
                              TmpEleIdx, TmpEleCfg,TmpEleNum,TmpEleSpn);
          } 
        }else{ //csz : t=s
          flag_hop = 1;
          Counter[0]++;
          makeCandidate_hopping_csz(&mi, &ri, &rj, &s,&t, &rejectFlag,
                              TmpEleIdx, TmpEleCfg,TmpEleNum,TmpEleSpn);
        } 
        StopTimer(31);

        if(rejectFlag) continue; 

        StartTimer(32);
        StartTimer(60);
        /* The mi-th electron with spin s hops to site rj with t */
        updateEleConfig_fsz(mi,ri,rj,s,t,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn);
        if(s==t){
          UpdateProjCnt(ri,rj,s,projCntNew,TmpEleProjCnt,TmpEleNum);
        }else{
          UpdateProjCnt_fsz(ri,rj,s,t,projCntNew,TmpEleProjCnt,TmpEleNum);
        }   
        StopTimer(60);
        StartTimer(61);
        CalculateNewPfM2_fsz_real(mi,t,pfMNew,TmpEleIdx,TmpEleSpn,qpStart,qpEnd); // fsz: s->t 
        StopTimer(61);

        StartTimer(62);
        /* calculate inner product <phi|L|x> */
        //logIpNew = CalculateLogIP_fcmp(pfMNew,qpStart,qpEnd,comm);
        logIpNew = CalculateLogIP_real(pfMNew,qpStart,qpEnd,comm);
        StopTimer(62);

        /* Metroplis */
        x = LogProjRatio(projCntNew,TmpEleProjCnt);
        w = exp(2.0*(x+(logIpNew-logIpOld)));
        if( !isfinite(w) ) w = -1.0; /* should be rejected */

        if(w > genrand_real2()) { /* accept */
          // UpdateMAll will change SlaterElm, InvM (including PfM)
          StartTimer(63);
          UpdateMAll_fsz_real(mi,t,TmpEleIdx,TmpEleSpn,qpStart,qpEnd); // fsz : s->t
          StopTimer(63);

          for(i=0;i<NProj;i++) TmpEleProjCnt[i] = projCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          if(flag_hop==1){// hopping
            Counter[1]++;
          }else{// local spin flip
            Counter[5]++;
          }
        } else { /* reject */ //(ri,s) <- (rj,t)
          revertEleConfig_fsz(mi,ri,rj,s,t,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn);
        }
        StopTimer(32);
      } else if(updateType==EXCHANGE) { /* exchange */
        Counter[2]++;

        StartTimer(31);
        makeCandidate_exchange_fsz(&mi, &ri, &rj, &s, &rejectFlag,
                               TmpEleIdx, TmpEleCfg, TmpEleNum,TmpEleSpn);
        StopTimer(31);
        if(rejectFlag) continue;

        StartTimer(33); //LSF
        StartTimer(65);

        /* The mi-th electron with spin s exchanges with the electron on site rj with spin 1-s */
        t = 1-s;
        mj = TmpEleCfg[rj+t*Nsite];

        /* The mi-th electron with spin s hops to rj */
        updateEleConfig_fsz(mi,ri,rj,s,s,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn);
        UpdateProjCnt(ri,rj,s,projCntNew,TmpEleProjCnt,TmpEleNum);
        /* The mj-th electron with spin t hops to ri */
        updateEleConfig_fsz(mj,rj,ri,t,t,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn);
        UpdateProjCnt(rj,ri,t,projCntNew,projCntNew,TmpEleNum);

        StopTimer(65);
        StartTimer(66);

        CalculateNewPfMTwo2_fsz_real(mi, s, mj, t, pfMNew, TmpEleIdx,TmpEleSpn, qpStart, qpEnd);
        StopTimer(66);
        StartTimer(67);

        /* calculate inner product <phi|L|x> */
        logIpNew = CalculateLogIP_real(pfMNew,qpStart,qpEnd,comm);

        StopTimer(67);

        /* Metroplis */
        x = LogProjRatio(projCntNew,TmpEleProjCnt);
        w = exp(2.0*(x+(logIpNew-logIpOld))); //TBC
        if( !isfinite(w) ) w = -1.0; /* should be rejected */

        if(w > genrand_real2()) { /* accept */
          StartTimer(68);
          UpdateMAllTwo_fsz_real(mi, s, mj, t, ri, rj, TmpEleIdx,TmpEleSpn,qpStart,qpEnd);
          StopTimer(68);

          for(i=0;i<NProj;i++) TmpEleProjCnt[i] = projCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          Counter[3]++;
        } else { /* reject */
          revertEleConfig_fsz(mj,rj,ri,t,t,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn);
          revertEleConfig_fsz(mi,ri,rj,s,s,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn);
        }
        StopTimer(33);
      }else if (updateType==LOCALSPINFLIP){
        Counter[4]++;

        StartTimer(31);
        makeCandidate_LocalSpinFlip_localspin(&mi, &ri, &rj, &s,&t, &rejectFlag,
                              TmpEleIdx, TmpEleCfg,TmpEleNum,TmpEleSpn);
        StopTimer(31);

        if(rejectFlag) continue; 

        StartTimer(36);
        StartTimer(600);
        /* The mi-th electron with spin s hops to site rj with t */
        // note we assume t=1-s,rj = ri
        //
        updateEleConfig_fsz(mi,ri,rj,s,t,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn);
        UpdateProjCnt_fsz(ri,rj,s,t,projCntNew,TmpEleProjCnt,TmpEleNum);
        StopTimer(600);
        StartTimer(601);
        CalculateNewPfM2_fsz_real(mi,t,pfMNew,TmpEleIdx,TmpEleSpn,qpStart,qpEnd); // fsz: s->t 
        StopTimer(610);

        StartTimer(602);
        /* calculate inner product <phi|L|x> */
        //logIpNew = CalculateLogIP_fcmp(pfMNew,qpStart,qpEnd,comm);
        logIpNew = CalculateLogIP_real(pfMNew,qpStart,qpEnd,comm);
        StopTimer(602);

        /* Metroplis */
        x = LogProjRatio(projCntNew,TmpEleProjCnt);
        w = exp(2.0*(x+(logIpNew-logIpOld)));
        if( !isfinite(w) ) w = -1.0; /* should be rejected */

        if(w > genrand_real2()) { /* accept */
          // UpdateMAll will change SlaterElm, InvM (including PfM)
          StartTimer(603);
          UpdateMAll_fsz_real(mi,t,TmpEleIdx,TmpEleSpn,qpStart,qpEnd); // fsz : s->t
          StopTimer(603);

          for(i=0;i<NProj;i++) TmpEleProjCnt[i] = projCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          Counter[5]++;
        } else { /* reject */ //(ri,s) <- (rj,t)
          revertEleConfig_fsz(mi,ri,rj,s,t,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn);
        }
        StopTimer(36);
      }

      if(nAccept>Nsite) {
        StartTimer(34);
        /* recal PfM and InvM */
        CalculateMAll_fsz_real(TmpEleIdx,TmpEleSpn,qpStart,qpEnd);
        //printf("DEBUG: maker3: PfM=%lf\n",creal(PfM[0]));
        logIpOld = CalculateLogIP_real(PfM_real,qpStart,qpEnd,comm);
        StopTimer(34);
        nAccept=0;
      }
    } /* end of instep */

    StartTimer(35);
    /* save Electron Configuration */
    if(outStep >= nOutStep-NVMCSample) {
      sample = outStep-(nOutStep-NVMCSample);
      #ifdef _DEBUG_DETAIL
      fprintf(stdout, "Debug: save Electron Configuration.\n");
      #endif
      saveEleConfig_fsz(sample,logIpOld,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,TmpEleSpn);
    }
    StopTimer(35);
  } /* end of outstep */
#ifdef _DEBUG_DETAIL
  fprintf(stdout, "Debug: finish step\n");
#endif

#ifdef _DEBUG_DETAIL
  fprintf(stdout, "Debug: copyToBurnSample_fsz\n");
#endif
  copyToBurnSample_fsz(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,TmpEleSpn);
#ifdef _DEBUG_DETAIL
  fprintf(stdout, "Debug: Finish copyToBurnSample_fsz\n");
#endif

  BurnFlag=1;
  return;
}

int makeInitialSample_fsz_real(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,int *eleSpn,
                      const int qpStart, const int qpEnd, MPI_Comm comm) {
  const int nsize = Nsize;
  const int nsite2 = Nsite2;
  int flag=1,flagRdc,loop=0;
  int ri,si,msi,rsi;
  int X_mi,tmp_TwoSz;
  int rank,size;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  do {
    /* initialize */
    #pragma omp parallel for default(shared) private(msi)
    for(msi=0;msi<nsize;msi++) eleIdx[msi] = -1;
    #pragma omp parallel for default(shared) private(msi)
    for(msi=0;msi<nsize;msi++) eleSpn[msi] = -1;
    #pragma omp parallel for default(shared) private(rsi)
    for(rsi=0;rsi<nsite2;rsi++) eleCfg[rsi] = -1;
    
    if(TwoSz==-1){
      tmp_TwoSz = 0;  //note: Sz is not conserved quantity but initially we take Sz=0 
    }else{
      tmp_TwoSz = TwoSz/2; // if TwoSz is not even, mVMC does not work 
    }
    //note:  2Sz=TwoSz X_mi=0-2*Ne=Nsize
    for(X_mi=0;X_mi<Nsize;X_mi++) {
      if(X_mi<Ne+tmp_TwoSz){
        eleSpn[X_mi]   = 0;
      }else{
        eleSpn[X_mi]   = 1;
      }
    }  
    /* local spin */
    for(ri=0;ri<Nsite;ri++) {
      if(LocSpn[ri]==1) {
        do {
          X_mi = gen_rand32()%Nsize;
          si = eleSpn[X_mi];
          //si = (genrand_real2()<0.5) ? 0 : 1;
        } while(eleIdx[X_mi]!=-1); // seeking empty site
        eleCfg[ri+si*Nsite] = X_mi;//;+si*Ne;
        eleIdx[X_mi]        = ri;
        //eleSpn[mi]    = si;
      }
    }
    /* itinerant electron */
    for(X_mi=0;X_mi<Nsize;X_mi++) { 
      si = eleSpn[X_mi];
      if(eleIdx[X_mi]== -1) {
        do {
          ri = gen_rand32()%Nsite;
        } while (eleCfg[ri+si*Nsite]!= -1 || LocSpn[ri]==1); // seeking empty and itinerant site
        eleCfg[ri+si*Nsite]     = X_mi; // buggged 4/26
        eleIdx[X_mi]            = ri;
        //eleSpn[mi+si*Ne]        = si;
      }
    }
    /* EleNum */
    #pragma omp parallel for default(shared) private(rsi)
    #pragma loop noalias
    for(rsi=0;rsi<nsite2;rsi++) {
      eleNum[rsi] = (eleCfg[rsi] < 0) ? 0 : 1;
    }
    
    MakeProjCnt(eleProjCnt,eleNum); // this function does not change even for fsz

    flag = CalculateMAll_fsz_real(eleIdx,eleSpn,qpStart,qpEnd);
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

#endif
