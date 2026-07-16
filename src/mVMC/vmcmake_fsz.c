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
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "vmcmake_fsz.h"
#ifndef _SRC_VMCMAKE_FSZ
#define _SRC_VMCMAKE_FSZ

#include "global.h"
#include "slater.h"
#include "matrix.h"
#include "pfupdate_fsz.h"
#include "qp.h"
#include "splitloop.h"

//typedef enum {HOPPING, EXCHANGE, NONE} UpdateType;
//UpdateType getUpdateType(int path);

void VMCMakeSample_fsz(MPI_Comm comm) {
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
  int flag_hop;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  SplitLoop(&qpStart,&qpEnd,NQPFull,rank,size);

  StartTimer(30);
  if(BurnFlag==0) {
#ifdef _DEBUG_DETAIL
    printf("DEBUG: make1: \n");
#endif
    makeInitialSample_fsz(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,TmpEleSpn,
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


#ifdef _pf_block_update
  // TODO: Compute from qpStart to qpEnd to support loop splitting.
  void *pfOrbital[NQPFull];
  void *pfUpdator[NQPFull];
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

  // Initialize with free spin configuration.
  updated_tdi_v_init_z(NQPFull, Nsite, Nsite2, Nsize,
                       SlaterElm, Nsite2*Nsite2,
                       InvM, Nsize*Nsize,
                       TmpEleIdx, TmpEleSpn,
                       NBlockUpdateSize,
                       pfUpdator, pfOrbital);
  updated_tdi_v_get_pfa_z(NQPFull, PfM, pfUpdator);
#else
  CalculateMAll_fsz(TmpEleIdx,TmpEleSpn,qpStart,qpEnd);
#endif
#ifdef _DEBUG_DETAIL
  printf("DEBUG: maker1: PfM=%lf\n",creal(PfM[0]));
#endif
  logIpOld = CalculateLogIP_fcmp(PfM,qpStart,qpEnd,comm);

  if( !isfinite(creal(logIpOld) + cimag(logIpOld)) ) {
    if(rank==0) fprintf(stderr,"waring: VMCMakeSample remakeSample logIpOld=%e\n",creal(logIpOld)); //TBC
    makeInitialSample_fsz(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,TmpEleSpn,
                      qpStart,qpEnd,comm);

#ifdef _pf_block_update
    // Clear and reinitialize.
    updated_tdi_v_free_z(NQPFull, pfUpdator, pfOrbital);
    updated_tdi_v_init_z(NQPFull, Nsite, Nsite2, Nsize,
                         SlaterElm, Nsite2*Nsite2,
                         InvM, Nsize*Nsize,
                         TmpEleIdx, TmpEleSpn,
                         NBlockUpdateSize,
                         pfUpdator, pfOrbital);
    updated_tdi_v_get_pfa_z(NQPFull, PfM, pfUpdator);
#else
    CalculateMAll_fsz(TmpEleIdx,TmpEleSpn,qpStart,qpEnd);
#endif
#ifdef _DEBUG_DETAIL
    printf("DEBUG: maker2: PfM=%lf\n",creal(PfM[0]));
#endif
    logIpOld = CalculateLogIP_fcmp(PfM,qpStart,qpEnd,comm);
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

      if (updateType == HOPPING || updateType == SPINHOPPING) { /* hopping or localized spin hopping */

        StartTimer(31);
        flag_hop = 0;
	switch(updateType){
	  case HOPPING:
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
	    break;
	  case SPINHOPPING:
            if(TwoSz==-1){//total spin is not conserved
              if(genrand_real2()<0.5){ // this ratio can be changed
                flag_hop = 1;
                Counter[0]++;
                makeCandidate_spin_hopping_fsz(&mi, &ri, &rj, &s,&t, &rejectFlag,
                              TmpEleIdx, TmpEleCfg,TmpEleNum,TmpEleSpn);
              }else{
                Counter[4]++;
                makeCandidate_LocalSpinFlip_conduction(&mi, &ri, &rj, &s,&t, &rejectFlag,
                              TmpEleIdx, TmpEleCfg,TmpEleNum,TmpEleSpn);
              }
            }else{ //csz : t=s
              flag_hop = 1;
              Counter[0]++;
              makeCandidate_spin_hopping_csz(&mi, &ri, &rj, &s,&t, &rejectFlag,
                              TmpEleIdx, TmpEleCfg,TmpEleNum,TmpEleSpn);
            }
	    break;
	  default:
	    if (rank == 0) fprintf(stderr, "error: strange updateType \n");
	    exit(1);
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
#ifdef _pf_block_update
        updated_tdi_v_push_z(NQPFull, rj+t*Nsite, mi, 1, pfUpdator);
        updated_tdi_v_get_pfa_z(NQPFull, pfMNew, pfUpdator);
#else
        CalculateNewPfM2_fsz(mi,t,pfMNew,TmpEleIdx,TmpEleSpn,qpStart,qpEnd); // fsz: s->t 
#endif
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
          UpdateMAll_fsz(mi,t,TmpEleIdx,TmpEleSpn,qpStart,qpEnd); // fsz : s->t
#endif
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
#ifdef _pf_block_update
          updated_tdi_v_pop_z(NQPFull, 0, pfUpdator);
#endif
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

#ifdef _pf_block_update
        updated_tdi_v_push_pair_z(NQPFull,
                                  rj+s*Nsite, mi,
                                  ri+t*Nsite, mj,
                                  1, pfUpdator);
        updated_tdi_v_get_pfa_z(NQPFull, pfMNew, pfUpdator);
#else
        CalculateNewPfMTwo2_fsz(mi, s, mj, t, pfMNew, TmpEleIdx,TmpEleSpn, qpStart, qpEnd);
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
          UpdateMAllTwo_fsz(mi, s, mj, t, ri, rj, TmpEleIdx,TmpEleSpn,qpStart,qpEnd);
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
#ifdef _pf_block_update
        updated_tdi_v_push_z(NQPFull, rj+t*Nsite, mi, 1, pfUpdator);
        updated_tdi_v_get_pfa_z(NQPFull, pfMNew, pfUpdator);
#else
        CalculateNewPfM2_fsz(mi,t,pfMNew,TmpEleIdx,TmpEleSpn,qpStart,qpEnd); // fsz: s->t 
#endif
        StopTimer(610);

        StartTimer(602);
        /* calculate inner product <phi|L|x> */
        //logIpNew = CalculateLogIP_fcmp(pfMNew,qpStart,qpEnd,comm);
        logIpNew = CalculateLogIP_fcmp(pfMNew,qpStart,qpEnd,comm);
        StopTimer(602);

        /* Metroplis */
        x = LogProjRatio(projCntNew,TmpEleProjCnt);
        w = exp(2.0*(x+creal(logIpNew-logIpOld)));
        if( !isfinite(w) ) w = -1.0; /* should be rejected */
        

        //printf("\n");
        //printf("MDEBUG: mi=%d: ri=%d rj=%d s=%d t=%d\n",mi,ri,rj,s,t);
        //printf("%lf: %d %d, %d %d, %d %d, %d %d \n",w,TmpEleNum[0+0*Nsite],TmpEleNum[0+1*Nsite],TmpEleNum[1+0*Nsite],TmpEleNum[1+1*Nsite],TmpEleNum[2+0*Nsite],TmpEleNum[2+1*Nsite],TmpEleNum[3+0*Nsite],TmpEleNum[3+1*Nsite]);
        //printf("\n");
        if(w > genrand_real2()) { /* accept */
          StartTimer(603);
#ifdef _pf_block_update
          // Inv already updated. Only need to get PfM again.
          updated_tdi_v_get_pfa_z(NQPFull, PfM, pfUpdator);
#else
          // UpdateMAll will change SlaterElm, InvM (including PfM)
          UpdateMAll_fsz(mi,t,TmpEleIdx,TmpEleSpn,qpStart,qpEnd); // fsz : s->t
#endif
          StopTimer(603);

          for(i=0;i<NProj;i++) TmpEleProjCnt[i] = projCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          Counter[5]++;
        } else { /* reject */ //(ri,s) <- (rj,t)
#ifdef _pf_block_update
          updated_tdi_v_pop_z(NQPFull, 0, pfUpdator);
#endif
          revertEleConfig_fsz(mi,ri,rj,s,t,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn);
        }
        StopTimer(36);
      }

      if(nAccept>Nsite) {
        // Recalculate PfM and InvM.
        StartTimer(34);
#ifdef _pf_block_update
        // Clear and reinitialize.
        updated_tdi_v_free_z(NQPFull, pfUpdator, pfOrbital);
        updated_tdi_v_init_z(NQPFull, Nsite, Nsite2, Nsize,
                             SlaterElm, Nsite2*Nsite2,
                             InvM, Nsize*Nsize,
                             TmpEleIdx, TmpEleSpn,
                             NBlockUpdateSize,
                             pfUpdator, pfOrbital);
        updated_tdi_v_get_pfa_z(NQPFull, PfM, pfUpdator);
#else
        CalculateMAll_fsz(TmpEleIdx,TmpEleSpn,qpStart,qpEnd);
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

#ifdef _pf_block_update
  // Free-up updator space.
  updated_tdi_v_free_z(NQPFull, pfUpdator, pfOrbital);
#endif

  return;
}

static int reduceBFInfo_fsz(const int info, MPI_Comm comm) {
  int size, infoSend, infoRdc;
  MPI_Comm_size(comm,&size);
  if(size <= 1) return info;
  infoSend = info;
  MPI_Allreduce(&infoSend,&infoRdc,1,MPI_INT,MPI_MAX,comm);
  return infoRdc;
}

static int BF_FSZ_PermuteParticleLabelsForCheck(int *eleIdx, int *eleCfg,
                                                int *eleSpn) {
  int a,b,siteA,siteB,spinA,spinB;
  if(!BFFSZPermuteParticleLabelsCheckEnabled) return 0;
  if(eleIdx == NULL || eleCfg == NULL || eleSpn == NULL || Nsize < 2) return 1;
  a = 0;
  for(b=1;b<Nsize;b++) if(eleSpn[b] != eleSpn[a]) break;
  if(b >= Nsize) return 1;
  siteA = eleIdx[a];
  siteB = eleIdx[b];
  spinA = eleSpn[a];
  spinB = eleSpn[b];
  if(siteA < 0 || siteA >= Nsite || siteB < 0 || siteB >= Nsite
      || (spinA != 0 && spinA != 1) || (spinB != 0 && spinB != 1)) return 1;
  eleIdx[a] = siteB;
  eleSpn[a] = spinB;
  eleIdx[b] = siteA;
  eleSpn[b] = spinA;
  eleCfg[siteB+spinB*Nsite] = a;
  eleCfg[siteA+spinA*Nsite] = b;
  return 0;
}

static int BF_FSZ_PfUpdateMatchesFull(const double complex *updated,
                                     const double complex *full,
                                     const int qpNum, int *badQp,
                                     double *badError) {
  int qpidx;
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    const double scale = fmax(1.0,cabs(full[qpidx]));
    const double error = cabs(updated[qpidx]-full[qpidx])/scale;
    if(!(isfinite(error) && error <= 1.0e-10)) {
      if(badQp != NULL) *badQp = qpidx;
      if(badError != NULL) *badError = error;
      return 0;
    }
  }
  return 1;
}

static int BF_FSZ_RunPfUpdateArgumentChecks(
    const int nAffected, const int *affected,
    double complex *pfMNew, double complex *oldPfMScratch,
    const double complex *oldInvM, const size_t invMQpStride,
    const int *eleIdx, const int *eleSpn,
    const int qpStart, const int qpEnd,
    const double complex *candidateSlater,
    double complex *complexWork, const size_t complexWorkCount,
    int *iwork, const size_t intWorkCount,
    double *rwork, const size_t rworkCount) {
  static int checked = 0;
  int badAffected = -1;
  int status;
  const size_t qpNum = (size_t)(qpEnd-qpStart);

  if(!BFFSZPfUpdateArgumentCheckEnabled || checked || qpNum == 0) return 0;
#define BF_FSZ_EXPECT_INVALID(call) do { \
    status = (call); \
    if(status != BF_FSZ_PF_UPDATE_INVALID_ARGUMENT) return 1; \
  } while(0)
  BF_FSZ_EXPECT_INVALID(CalculateNewPfMBF_fsz_rows_workspace(
      0,affected,pfMNew,PfM,oldInvM,invMQpStride,eleIdx,eleSpn,
      qpStart,qpEnd,candidateSlater,NULL,complexWork,complexWorkCount,
      iwork,intWorkCount,rwork,rworkCount));
  BF_FSZ_EXPECT_INVALID(CalculateNewPfMBF_fsz_rows_workspace(
      -1,affected,pfMNew,PfM,oldInvM,invMQpStride,eleIdx,eleSpn,
      qpStart,qpEnd,candidateSlater,NULL,complexWork,complexWorkCount,
      iwork,intWorkCount,rwork,rworkCount));
  BF_FSZ_EXPECT_INVALID(CalculateNewPfMBF_fsz_rows_workspace(
      Nsize+1,affected,pfMNew,PfM,oldInvM,invMQpStride,eleIdx,eleSpn,
      qpStart,qpEnd,candidateSlater,NULL,complexWork,complexWorkCount,
      iwork,intWorkCount,rwork,rworkCount));
  BF_FSZ_EXPECT_INVALID(CalculateNewPfMBF_fsz_rows_workspace(
      1,&badAffected,pfMNew,PfM,oldInvM,invMQpStride,eleIdx,eleSpn,
      qpStart,qpEnd,candidateSlater,NULL,complexWork,complexWorkCount,
      iwork,intWorkCount,rwork,rworkCount));
  BF_FSZ_EXPECT_INVALID(CalculateNewPfMBF_fsz_rows_workspace(
      nAffected,affected,pfMNew,PfM,oldInvM,SIZE_MAX,eleIdx,eleSpn,
      qpStart,qpEnd,candidateSlater,NULL,complexWork,complexWorkCount,
      iwork,intWorkCount,rwork,rworkCount));
  BF_FSZ_EXPECT_INVALID(CalculateNewPfMBF_fsz_rows_workspace(
      nAffected,affected,pfMNew,PfM,oldInvM,invMQpStride,eleIdx,eleSpn,
      qpStart,qpEnd,candidateSlater,NULL,complexWork,complexWorkCount-1,
      iwork,intWorkCount,rwork,rworkCount));
  BF_FSZ_EXPECT_INVALID(CalculateNewPfMBF_fsz_rows_workspace(
      nAffected,affected,(double complex *)PfM,PfM,oldInvM,invMQpStride,
      eleIdx,eleSpn,qpStart,qpEnd,candidateSlater,NULL,complexWork,
      complexWorkCount,iwork,intWorkCount,rwork,rworkCount));
  BF_FSZ_EXPECT_INVALID(CalculateNewPfMBF_fsz_rows_workspace(
      nAffected,affected,pfMNew,PfM,oldInvM,invMQpStride,eleIdx,eleSpn,
      qpStart,qpEnd,SlaterElmBF,NULL,complexWork,complexWorkCount,
      iwork,intWorkCount,rwork,rworkCount));
#undef BF_FSZ_EXPECT_INVALID

  memcpy(oldPfMScratch,PfM,sizeof(double complex)*qpNum);
  oldPfMScratch[0] = 0.0 + 0.0*I;
  status = CalculateNewPfMBF_fsz_rows_workspace(
      nAffected,affected,pfMNew,oldPfMScratch,oldInvM,invMQpStride,
      eleIdx,eleSpn,qpStart,qpEnd,candidateSlater,NULL,complexWork,
      complexWorkCount,iwork,intWorkCount,rwork,rworkCount);
  if(status != BF_FSZ_PF_UPDATE_EXACT_ZERO) return 1;
  checked = 1;
  return 0;
}

static int BF_FSZ_InvUpdateMatchesFull(const double complex *updated,
                                       const double complex *full,
                                       const size_t count, size_t *badIndex,
                                       double *badError) {
  size_t idx;
  for(idx=0;idx<count;idx++) {
    const double scale = fmax(1.0,cabs(full[idx]));
    const double error = cabs(updated[idx]-full[idx])/scale;
    if(!(isfinite(error) && error <= 1.0e-9)) {
      if(badIndex != NULL) *badIndex = idx;
      if(badError != NULL) *badError = error;
      return 0;
    }
  }
  return 1;
}

static int BF_FSZ_RunInvUpdateArgumentChecks(
    const int nAffected, const int *affected,
    const double complex *oldInvM, const size_t invMQpStride,
    double complex *pfMNew, double complex *invMNew,
    const int *eleIdx, const int *eleSpn,
    const int qpStart, const int qpEnd,
    const double complex *candidateSlater,
    double complex *complexWork, const size_t complexWorkCount,
    int *iwork, const size_t intWorkCount,
    double complex *fullBufM, int *fullIWork,
    double complex *fullWork, const int fullLWork, double *fullRWork) {
  static int checked = 0;
  int badAffected = -1;
  BF_FSZ_InvUpdateResult result;
  BF_FSZ_MAllResult fullResult;
  if(!BFFSZInvUpdateArgumentCheckEnabled || checked || qpEnd == qpStart) return 0;
#define BF_FSZ_EXPECT_INV_INVALID(call) do { \
    result = (call); \
    if(result.status != BF_FSZ_INV_UPDATE_INVALID_ARGUMENT) return 1; \
  } while(0)
  BF_FSZ_EXPECT_INV_INVALID(PrepareInvMBF_fsz_rows_workspace(
      0,affected,oldInvM,invMQpStride,invMNew,invMQpStride,eleIdx,eleSpn,
      qpStart,qpEnd,candidateSlater,complexWork,complexWorkCount,
      iwork,intWorkCount));
  BF_FSZ_EXPECT_INV_INVALID(PrepareInvMBF_fsz_rows_workspace(
      -1,affected,oldInvM,invMQpStride,invMNew,invMQpStride,eleIdx,eleSpn,
      qpStart,qpEnd,candidateSlater,complexWork,complexWorkCount,
      iwork,intWorkCount));
  BF_FSZ_EXPECT_INV_INVALID(PrepareInvMBF_fsz_rows_workspace(
      Nsize+1,affected,oldInvM,invMQpStride,invMNew,invMQpStride,eleIdx,eleSpn,
      qpStart,qpEnd,candidateSlater,complexWork,complexWorkCount,
      iwork,intWorkCount));
  BF_FSZ_EXPECT_INV_INVALID(PrepareInvMBF_fsz_rows_workspace(
      1,&badAffected,oldInvM,invMQpStride,invMNew,invMQpStride,eleIdx,eleSpn,
      qpStart,qpEnd,candidateSlater,complexWork,complexWorkCount,
      iwork,intWorkCount));
  BF_FSZ_EXPECT_INV_INVALID(PrepareInvMBF_fsz_rows_workspace(
      nAffected,affected,oldInvM,SIZE_MAX,invMNew,invMQpStride,eleIdx,eleSpn,
      qpStart,qpEnd,candidateSlater,complexWork,complexWorkCount,
      iwork,intWorkCount));
  BF_FSZ_EXPECT_INV_INVALID(PrepareInvMBF_fsz_rows_workspace(
      nAffected,affected,oldInvM,invMQpStride,(double complex *)oldInvM,
      invMQpStride,eleIdx,eleSpn,qpStart,qpEnd,candidateSlater,
      complexWork,complexWorkCount,iwork,intWorkCount));
  BF_FSZ_EXPECT_INV_INVALID(PrepareInvMBF_fsz_rows_workspace(
      nAffected,affected,oldInvM,invMQpStride,invMNew,invMQpStride,
      eleIdx,eleSpn,qpStart,qpEnd,SlaterElmBF,complexWork,
      complexWorkCount,iwork,intWorkCount));
  BF_FSZ_EXPECT_INV_INVALID(PrepareInvMBF_fsz_rows_workspace(
      nAffected,affected,oldInvM,invMQpStride,invMNew,invMQpStride,
      eleIdx,eleSpn,qpStart,qpEnd,candidateSlater,complexWork,
      0,iwork,intWorkCount));
#undef BF_FSZ_EXPECT_INV_INVALID
  fullResult = CalculateMAll_BF_fsz_from_workspace(
      candidateSlater,eleIdx,eleSpn,qpStart,qpEnd,pfMNew,invMNew,SIZE_MAX,
      fullBufM,fullIWork,fullWork,fullLWork,fullRWork);
  if(fullResult.status != BF_FSZ_MALL_INVALID_ARGUMENT) return 1;
  fullResult = CalculateMAll_BF_fsz_from_workspace(
      candidateSlater,eleIdx,eleSpn,qpStart,qpEnd,pfMNew,invMNew,
      invMQpStride,fullBufM,fullIWork,fullWork,0,fullRWork);
  if(fullResult.status != BF_FSZ_MALL_INVALID_ARGUMENT) return 1;
  checked = 1;
  return 0;
}

static BF_FSZ_InvUpdateResult BF_FSZ_InjectedInvResult(
    BF_FSZ_InvUpdateResult result, const int rank) {
  if(result.status != BF_FSZ_INV_UPDATE_OK) return result;
  if(BFFSZInvUpdateInjectedRank >= 0 && rank != BFFSZInvUpdateInjectedRank) {
    return result;
  }
  if(BFFSZInvUpdateForceFallback) {
    result.status = BF_FSZ_INV_UPDATE_NONFINITE;
    result.stage = BF_FSZ_INV_STAGE_FINITE;
    return result;
  }
  if(BFFSZInvUpdateInjectedStage == BF_FSZ_INV_STAGE_GETRF
      || BFFSZInvUpdateInjectedStage == BF_FSZ_INV_STAGE_GETRI) {
    result.status = BF_FSZ_INV_UPDATE_LAPACK_FAILURE;
    result.stage = BFFSZInvUpdateInjectedStage;
    result.lapackInfo = 1;
  } else if(BFFSZInvUpdateInjectedStage == BF_FSZ_INV_STAGE_FINITE) {
    result.status = BF_FSZ_INV_UPDATE_NONFINITE;
    result.stage = BF_FSZ_INV_STAGE_FINITE;
  } else if(BFFSZInvUpdateInjectedStage == BF_FSZ_INV_STAGE_ANTISYMMETRY) {
    result.status = BF_FSZ_INV_UPDATE_ANTISYMMETRY;
    result.stage = BF_FSZ_INV_STAGE_ANTISYMMETRY;
  } else if(BFFSZInvUpdateInjectedStage == BF_FSZ_INV_STAGE_RESIDUAL) {
    result.status = BF_FSZ_INV_UPDATE_RESIDUAL;
    result.stage = BF_FSZ_INV_STAGE_RESIDUAL;
  }
  return result;
}

void VMC_BF_MakeSample_fsz(MPI_Comm comm) {
  const size_t bfSlaterSize = (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
  int outStep,nOutStep;
  int inStep,nInStep;
  int mi,ri,rj,s,t,i;
  int nAccept=0;
  int sample;

  double complex logIpOld,logIpNew;
  double complex *candidateSlater;
  double complex *pfMNew;
  double complex *pfMCheck;
  double complex *pfUpdateWork;
  double complex *invMNew, *invMCheck, *invUpdateWork;
  double complex *pfFullBufM, *pfFullWork;
  double complex *oldPfMCheck, *oldInvMCheck;
  int *pfUpdateIWork, *invUpdateIWork;
  int *pfFullIWork;
  double *pfUpdateRWork, *pfFullRWork;
  size_t pfUpdateComplexCount, pfUpdateIntCount, pfUpdateDoubleCount;
  size_t invUpdateComplexCount, invUpdateIntCount;
  size_t qpNum, nsizeSquared, oldInvCount;
  int *hopIntStorage, *affected, *hopIntWork;
  int hopIntWorkSize, hopIntStorageSize;
  int nChanged,nAffected;
  int projCntNew[NProj];
  int projBFCntNew[16*Nsite*Nrange];
  double x,w;

  int qpStart,qpEnd;
  int rejectFlag;
  int info,infoGlobal;
  int rank,size;

  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  SplitLoop(&qpStart,&qpEnd,NQPFull,rank,size);
  qpNum = (size_t)(qpEnd-qpStart);
  nsizeSquared = (size_t)Nsize*(size_t)Nsize;

  if(GetSlaterElmBF_fsz_hop_int_work_size(&hopIntWorkSize) != 0
      || Nsize > INT_MAX - hopIntWorkSize) {
    fprintf(stderr, "error: invalid BF-FSZ sampling integer workspace size\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(GetCalculateNewPfMBF_fsz_rows_work_size(&pfUpdateComplexCount,
      &pfUpdateIntCount, &pfUpdateDoubleCount) != 0
      || GetPrepareInvMBF_fsz_rows_work_size(&invUpdateComplexCount,
          &invUpdateIntCount) != 0
      || pfUpdateComplexCount > SIZE_MAX/sizeof(double complex)
      || pfUpdateIntCount > SIZE_MAX/sizeof(int)
      || pfUpdateDoubleCount > SIZE_MAX/sizeof(double)
      || invUpdateComplexCount > SIZE_MAX/sizeof(double complex)
      || invUpdateIntCount > SIZE_MAX/sizeof(int)
      || nsizeSquared > SIZE_MAX/sizeof(double complex)
      || (size_t)LapackLWork > SIZE_MAX/sizeof(double complex)
      || (size_t)LapackLWork > SIZE_MAX/sizeof(double)
      || (qpNum != 0 && nsizeSquared > SIZE_MAX/qpNum)) {
    fprintf(stderr, "error: invalid BF-FSZ sampling Pfaffian-update workspace size\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  hopIntStorageSize = Nsize + hopIntWorkSize;
  if((size_t)hopIntStorageSize > (size_t)-1/sizeof(int)) {
    fprintf(stderr, "error: BF-FSZ sampling integer workspace byte size overflow\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  candidateSlater = (double complex *)malloc(sizeof(double complex)*bfSlaterSize);
  pfMNew = (double complex *)malloc(sizeof(double complex)*(size_t)NQPFull);
  pfMCheck = (double complex *)malloc(sizeof(double complex)*(size_t)NQPFull);
  hopIntStorage = (int *)malloc(sizeof(int)*(size_t)hopIntStorageSize);
  pfUpdateWork = (double complex *)malloc(sizeof(double complex)*pfUpdateComplexCount);
  pfUpdateIWork = (int *)malloc(sizeof(int)*pfUpdateIntCount);
  pfUpdateRWork = (double *)malloc(sizeof(double)*pfUpdateDoubleCount);
  invUpdateWork = (double complex *)malloc(
      sizeof(double complex)*invUpdateComplexCount);
  invUpdateIWork = (int *)malloc(sizeof(int)*invUpdateIntCount);
  pfFullBufM = (double complex *)malloc(sizeof(double complex)*nsizeSquared);
  pfFullIWork = (int *)malloc(sizeof(int)*(size_t)Nsize);
  pfFullWork = (double complex *)malloc(sizeof(double complex)*(size_t)LapackLWork);
  pfFullRWork = (double *)malloc(sizeof(double)*(size_t)LapackLWork);
  oldInvCount = qpNum*nsizeSquared;
  invMNew = (double complex *)malloc(sizeof(double complex)
      *((oldInvCount > 0) ? oldInvCount : 1));
  invMCheck = NULL;
  if(qpNum > 0 && BFFSZInvUpdateCheckEnabled) {
    invMCheck = (double complex *)malloc(sizeof(double complex)*oldInvCount);
  }
  oldPfMCheck = NULL;
  oldInvMCheck = NULL;
  if(qpNum > 0 && (BFFSZPfUpdateExplicitStateCheckEnabled
      || BFFSZPfUpdateArgumentCheckEnabled
      || BFFSZInvUpdateExplicitStateCheckEnabled)) {
    oldPfMCheck = (double complex *)malloc(sizeof(double complex)*qpNum);
    oldInvMCheck = (double complex *)malloc(sizeof(double complex)*oldInvCount);
  }
  if(candidateSlater == NULL || pfMNew == NULL || pfMCheck == NULL
      || hopIntStorage == NULL
      || pfUpdateWork == NULL || pfUpdateIWork == NULL || pfUpdateRWork == NULL
      || invUpdateWork == NULL || invUpdateIWork == NULL
      || invMNew == NULL
      || (qpNum > 0 && BFFSZInvUpdateCheckEnabled && invMCheck == NULL)
      || pfFullBufM == NULL || pfFullIWork == NULL || pfFullWork == NULL
      || pfFullRWork == NULL
      || (qpNum > 0 && (BFFSZPfUpdateExplicitStateCheckEnabled
          || BFFSZPfUpdateArgumentCheckEnabled
          || BFFSZInvUpdateExplicitStateCheckEnabled)
          && (oldPfMCheck == NULL || oldInvMCheck == NULL))) {
    fprintf(stderr, "error: failed to allocate BF-FSZ sampling candidate workspace\n");
    free(candidateSlater);
    free(pfMNew);
    free(pfMCheck);
    free(hopIntStorage);
    free(pfUpdateWork);
    free(pfUpdateIWork);
    free(pfUpdateRWork);
    free(invUpdateWork);
    free(invUpdateIWork);
    free(invMNew);
    free(invMCheck);
    free(pfFullBufM);
    free(pfFullIWork);
    free(pfFullWork);
    free(pfFullRWork);
    free(oldPfMCheck);
    free(oldInvMCheck);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  affected = hopIntStorage;
  hopIntWork = affected + Nsize;

  StartTimer(30);
  if(BurnFlag==0) {
    makeInitialSampleBF_fsz(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,TmpEleProjBFCnt,TmpEleSpn,
                            qpStart,qpEnd,comm);
  } else {
    copyFromBurnSampleBF_fsz(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,TmpEleProjBFCnt,TmpEleSpn);
    MakeSlaterElmBF_fsz(TmpEleNum,TmpEleProjBFCnt);
  }

  info = CalculateMAll_BF_fsz(TmpEleIdx,TmpEleSpn,qpStart,qpEnd);
  infoGlobal = reduceBFInfo_fsz(info,comm);
  logIpOld = (infoGlobal == 0) ? CalculateLogIP_fcmp(PfM,qpStart,qpEnd,comm) : 0.0 + 0.0*I;
  if(infoGlobal != 0 || !isfinite(creal(logIpOld) + cimag(logIpOld))) {
    if(rank==0) fprintf(stderr,"warning: VMC_BF_MakeSample_fsz remakeSample logIpOld=%e\n",creal(logIpOld));
    makeInitialSampleBF_fsz(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,TmpEleProjBFCnt,TmpEleSpn,
                            qpStart,qpEnd,comm);
    info = CalculateMAll_BF_fsz(TmpEleIdx,TmpEleSpn,qpStart,qpEnd);
    infoGlobal = reduceBFInfo_fsz(info,comm);
    if(infoGlobal != 0) {
      if(rank==0) fprintf(stderr, "error: VMC_BF_MakeSample_fsz failed to rebuild initial sample\n");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    logIpOld = CalculateLogIP_fcmp(PfM,qpStart,qpEnd,comm);
    BurnFlag = 0;
  }
  StopTimer(30);

  nOutStep = (BurnFlag==0) ? NVMCWarmUp+NVMCSample : NVMCSample+1;
  nInStep = NVMCInterval * Nsite;

  for(i=0;i<Counter_max;i++) Counter[i]=0;

  for(outStep=0;outStep<nOutStep;outStep++) {
    for(inStep=0;inStep<nInStep;inStep++) {
      StartTimer(31);
      Counter[0]++;
      makeCandidate_hopping_csz(&mi, &ri, &rj, &s, &t, &rejectFlag,
                                TmpEleIdx, TmpEleCfg,TmpEleNum,TmpEleSpn);
      StopTimer(31);

      if(rejectFlag) continue;

      StartTimer(32);
      StartTimer(60);
      updateEleConfig_fsz(mi,ri,rj,s,t,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn);
      UpdateProjCnt(ri,rj,s,projCntNew,TmpEleProjCnt,TmpEleNum);
      MakeProjBFCnt(projBFCntNew,TmpEleNum);
      StopTimer(60);

      StartTimer(64);
      info = MakeSlaterElmBF_fsz_hop_to_with_rows_workspace(
          candidateSlater, SlaterElmBF, mi, TmpEleIdx, TmpEleSpn,
          TmpEleProjBFCnt, projBFCntNew, affected, &nChanged, &nAffected,
          hopIntWork, hopIntWorkSize);
      if(info != 0) {
        if(rank==0) fprintf(stderr, "error: BF-FSZ sampling affected-row collection failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
      RecordBFFSZProfile(BFFSZ_PROFILE_SAMPLE,nChanged,nAffected);
      StopTimer(64);

      StartTimer(61);
      if(BF_FSZ_ShouldUseFullPfaffian(nAffected)) {
        RecordBFFSZPfPath(BFFSZ_PROFILE_SAMPLE,BFFSZ_PF_PATH_DIRECT_FULL);
        info = CalculatePfM_BF_fsz_from_workspace(candidateSlater,TmpEleIdx,
            TmpEleSpn,qpStart,qpEnd,pfMNew,NULL,pfFullBufM,pfFullIWork,
            pfFullWork,LapackLWork,pfFullRWork);
      } else {
        int updateDetail = 0;
        int updateStatus;
        const double complex *oldPfMInput = PfM;
        const double complex *oldInvMInput = InvM;
        double complex savedPfM = 0.0 + 0.0*I;
        double complex savedInvM = 0.0 + 0.0*I;
        if(BF_FSZ_RunPfUpdateArgumentChecks(nAffected,affected,pfMNew,
            oldPfMCheck,InvM,nsizeSquared,TmpEleIdx,TmpEleSpn,qpStart,qpEnd,
            candidateSlater,pfUpdateWork,pfUpdateComplexCount,pfUpdateIWork,
            pfUpdateIntCount,pfUpdateRWork,pfUpdateDoubleCount) != 0) {
          if(rank==0) fprintf(stderr,
              "error: BF-FSZ Pfaffian-update argument/status check failed\n");
          MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
        }
        if(BFFSZPfUpdateExplicitStateCheckEnabled && qpNum > 0) {
          memcpy(oldPfMCheck,PfM,sizeof(double complex)*qpNum);
          memcpy(oldInvMCheck,InvM,sizeof(double complex)*oldInvCount);
          oldPfMInput = oldPfMCheck;
          oldInvMInput = oldInvMCheck;
          savedPfM = PfM[0];
          savedInvM = InvM[0];
          PfM[0] = NAN + NAN*I;
          InvM[0] = NAN + NAN*I;
        }
        updateStatus = CalculateNewPfMBF_fsz_rows_workspace(
            nAffected, affected, pfMNew, oldPfMInput, oldInvMInput, nsizeSquared,
            TmpEleIdx, TmpEleSpn, qpStart, qpEnd, candidateSlater, &updateDetail,
            pfUpdateWork, pfUpdateComplexCount, pfUpdateIWork, pfUpdateIntCount,
            pfUpdateRWork, pfUpdateDoubleCount);
        if(BFFSZPfUpdateExplicitStateCheckEnabled && qpNum > 0) {
          PfM[0] = savedPfM;
          InvM[0] = savedInvM;
        }
        if(updateStatus == BF_FSZ_PF_UPDATE_OK) {
          if(BFFSZPfUpdateInjectedStatus != BF_FSZ_PF_UPDATE_OK) {
            updateStatus = BFFSZPfUpdateInjectedStatus;
          } else if(BFFSZPfUpdateForceFallback) {
            updateStatus = BF_FSZ_PF_UPDATE_NONFINITE;
          }
        }
        if(updateStatus == BF_FSZ_PF_UPDATE_INVALID_ARGUMENT) {
          if(rank==0) fprintf(stderr,
              "error: invalid BF-FSZ sampling Pfaffian-update arguments\n");
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        if(updateStatus == BF_FSZ_PF_UPDATE_OK) {
          RecordBFFSZPfPath(BFFSZ_PROFILE_SAMPLE,BFFSZ_PF_PATH_OPTIMIZED);
          info = BF_FSZ_PF_OK;
          if(BFFSZPfUpdateCheckEnabled) {
            int badQp = -1;
            double badError = 0.0;
            double ipError;
            double complex ipUpdated,ipFull;
            info = CalculatePfM_BF_fsz_from_workspace(candidateSlater,
                TmpEleIdx,TmpEleSpn,qpStart,qpEnd,pfMCheck,NULL,pfFullBufM,
                pfFullIWork,pfFullWork,LapackLWork,pfFullRWork);
            if(info != BF_FSZ_PF_OK || !BF_FSZ_PfUpdateMatchesFull(
                pfMNew,pfMCheck,qpEnd-qpStart,&badQp,&badError)) {
              if(rank==0) fprintf(stderr,
                  "error: BF-FSZ sampling Pfaffian-update/full mismatch at qp=%d relative_error=%.17e\n",
                  badQp,badError);
              MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            ipUpdated = CalculateIP_fcmp(pfMNew,qpStart,qpEnd,comm);
            ipFull = CalculateIP_fcmp(pfMCheck,qpStart,qpEnd,comm);
            ipError = cabs(ipUpdated-ipFull)/fmax(1.0,cabs(ipFull));
            if(!(isfinite(ipError) && ipError <= 1.0e-10)) {
              if(rank==0) fprintf(stderr,
                  "error: BF-FSZ sampling Pfaffian-update/full IP mismatch relative_error=%.17e\n",
                  ipError);
              MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
            }
          }
        } else {
          RecordBFFSZPfPath(BFFSZ_PROFILE_SAMPLE,BFFSZ_PF_PATH_FALLBACK);
          info = CalculatePfM_BF_fsz_from_workspace(candidateSlater,TmpEleIdx,
              TmpEleSpn,qpStart,qpEnd,pfMNew,NULL,pfFullBufM,pfFullIWork,
              pfFullWork,LapackLWork,pfFullRWork);
        }
      }
      StopTimer(61);

      StartTimer(62);
      infoGlobal = reduceBFInfo_fsz(info,comm);
      logIpNew = (infoGlobal == 0) ? CalculateLogIP_fcmp(pfMNew,qpStart,qpEnd,comm) : 0.0 + 0.0*I;
      StopTimer(62);

      x = LogProjRatio(projCntNew,TmpEleProjCnt);
      w = (infoGlobal == 0) ? exp(2.0*(x+creal(logIpNew-logIpOld))) : -1.0;
      if(!isfinite(w)) w = -1.0;

      if(w > genrand_real2()) {
        BF_FSZ_MAllResult fullResult;
        BF_FSZ_InvUpdateResult invResult;
        int invStatusGlobal;
        StartTimer(63);
        if(BF_FSZ_ShouldUseFullPfaffian(nAffected)) {
          fullResult = CalculateMAll_BF_fsz_from_workspace(candidateSlater,
              TmpEleIdx,TmpEleSpn,qpStart,qpEnd,pfMNew,invMNew,nsizeSquared,
              pfFullBufM,pfFullIWork,pfFullWork,LapackLWork,pfFullRWork);
          invStatusGlobal = reduceBFInfo_fsz(fullResult.status,comm);
          if(invStatusGlobal != BF_FSZ_MALL_OK) {
            if(rank==0) fprintf(stderr,
                "error: BF-FSZ accepted direct full inverse prepare failed\n");
            MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
          }
          RecordBFFSZInvPath(BFFSZ_PF_PATH_DIRECT_FULL);
        } else {
          const double complex *oldInvInput = InvM;
          double complex savedInvM = 0.0+0.0*I;
          if(BF_FSZ_RunInvUpdateArgumentChecks(nAffected,affected,InvM,
              nsizeSquared,pfMNew,invMNew,TmpEleIdx,TmpEleSpn,qpStart,qpEnd,
              candidateSlater,invUpdateWork,invUpdateComplexCount,
              invUpdateIWork,invUpdateIntCount,pfFullBufM,pfFullIWork,
              pfFullWork,LapackLWork,pfFullRWork) != 0) {
            if(rank==0) fprintf(stderr,
                "error: BF-FSZ inverse-update argument check failed\n");
            MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
          }
          if(BFFSZInvUpdateExplicitStateCheckEnabled && qpNum > 0) {
            memcpy(oldInvMCheck,InvM,sizeof(double complex)*oldInvCount);
            oldInvInput = oldInvMCheck;
            savedInvM = InvM[0];
            InvM[0] = NAN+NAN*I;
          }
          invResult = PrepareInvMBF_fsz_rows_workspace(
              nAffected,affected,oldInvInput,nsizeSquared,invMNew,
              nsizeSquared,TmpEleIdx,TmpEleSpn,qpStart,qpEnd,
              candidateSlater,invUpdateWork,invUpdateComplexCount,
              invUpdateIWork,invUpdateIntCount);
          if(BFFSZInvUpdateExplicitStateCheckEnabled && qpNum > 0) {
            InvM[0] = savedInvM;
          }
          RecordBFFSZInvChecks(invResult.checkSeconds,
              invResult.antisymmetryResidual,invResult.affectedResidual);
          invResult = BF_FSZ_InjectedInvResult(invResult,rank);
          if(invResult.status == BF_FSZ_INV_UPDATE_INVALID_ARGUMENT) {
            if(rank==0) fprintf(stderr,
                "error: invalid BF-FSZ accepted inverse-update arguments\n");
            MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
          }
          invStatusGlobal = reduceBFInfo_fsz(invResult.status,comm);
          if(invStatusGlobal == BF_FSZ_INV_UPDATE_OK) {
            if(BFFSZInvUpdateCheckEnabled) {
              size_t badIndex = 0;
              double badError = 0.0;
              fullResult = CalculateMAll_BF_fsz_from_workspace(
                  candidateSlater,TmpEleIdx,TmpEleSpn,qpStart,qpEnd,
                  pfMCheck,invMCheck,nsizeSquared,pfFullBufM,pfFullIWork,
                  pfFullWork,LapackLWork,pfFullRWork);
              infoGlobal = reduceBFInfo_fsz(fullResult.status,comm);
              if(infoGlobal != BF_FSZ_MALL_OK
                  || !BF_FSZ_PfUpdateMatchesFull(pfMNew,pfMCheck,
                      qpEnd-qpStart,NULL,NULL)
                  || !BF_FSZ_InvUpdateMatchesFull(invMNew,invMCheck,
                      oldInvCount,&badIndex,&badError)) {
                if(rank==0) fprintf(stderr,
                    "error: BF-FSZ accepted inverse-update/full mismatch at index=%zu relative_error=%.17e\n",
                    badIndex,badError);
                MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
              }
            }
            RecordBFFSZInvPath(BFFSZ_PF_PATH_OPTIMIZED);
          } else {
            fullResult = CalculateMAll_BF_fsz_from_workspace(
                candidateSlater,TmpEleIdx,TmpEleSpn,qpStart,qpEnd,
                pfMNew,invMNew,nsizeSquared,pfFullBufM,pfFullIWork,
                pfFullWork,LapackLWork,pfFullRWork);
            infoGlobal = reduceBFInfo_fsz(fullResult.status,comm);
            if(infoGlobal != BF_FSZ_MALL_OK) {
              if(rank==0) fprintf(stderr,
                  "error: BF-FSZ accepted full inverse fallback failed\n");
              MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
            }
            RecordBFFSZInvPath(BFFSZ_PF_PATH_FALLBACK);
          }
        }
        memcpy(SlaterElmBF,candidateSlater,
            sizeof(double complex)*bfSlaterSize);
        if(qpNum > 0) {
          memcpy(PfM,pfMNew,sizeof(double complex)*qpNum);
          memcpy(InvM,invMNew,sizeof(double complex)*oldInvCount);
        }
        StopTimer(63);
        for(i=0;i<NProj;i++) TmpEleProjCnt[i] = projCntNew[i];
        for(i=0;i<16*Nsite*Nrange;i++) TmpEleProjBFCnt[i] = projBFCntNew[i];
        logIpOld = CalculateLogIP_fcmp(PfM,qpStart,qpEnd,comm);
        nAccept++;
        Counter[1]++;
      } else {
        revertEleConfig_fsz(mi,ri,rj,s,t,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleSpn);
      }
      StopTimer(32);

      if(nAccept>Nsite) {
        StartTimer(34);
        MakeSlaterElmBF_fsz(TmpEleNum,TmpEleProjBFCnt);
        info = CalculateMAll_BF_fsz(TmpEleIdx,TmpEleSpn,qpStart,qpEnd);
        infoGlobal = reduceBFInfo_fsz(info,comm);
        if(infoGlobal != 0) {
          if(rank==0) fprintf(stderr, "error: VMC_BF_MakeSample_fsz failed to refresh accepted sample\n");
          MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
        }
        logIpOld = CalculateLogIP_fcmp(PfM,qpStart,qpEnd,comm);
        StopTimer(34);
        nAccept=0;
      }
    }

    StartTimer(35);
    if(outStep >= nOutStep-NVMCSample) {
      sample = outStep-(nOutStep-NVMCSample);
      saveEleConfigBF_fsz(sample,logIpOld,TmpEleIdx,TmpEleCfg,TmpEleNum,
                          TmpEleProjCnt,TmpEleProjBFCnt,TmpEleSpn);
    }
    StopTimer(35);
  }

  copyToBurnSampleBF_fsz(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,TmpEleProjBFCnt,TmpEleSpn);
  BurnFlag=1;

  free(candidateSlater);
  free(pfMNew);
  free(pfMCheck);
  free(hopIntStorage);
  free(pfUpdateWork);
  free(pfUpdateIWork);
  free(pfUpdateRWork);
  free(invUpdateWork);
  free(invUpdateIWork);
  free(invMNew);
  free(invMCheck);
  free(pfFullBufM);
  free(pfFullIWork);
  free(pfFullWork);
  free(pfFullRWork);
  free(oldPfMCheck);
  free(oldInvMCheck);

  return;
}

int makeInitialSampleBF_fsz(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,
                            int *eleProjBFCnt, int *eleSpn,
                            const int qpStart, const int qpEnd, MPI_Comm comm) {
  int flag=1,flagRdc,loop=0;
  int rank,size;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  do {
    makeInitialSample_fsz(eleIdx,eleCfg,eleNum,eleProjCnt,eleSpn,qpStart,qpEnd,comm);
    if(BF_FSZ_PermuteParticleLabelsForCheck(eleIdx,eleCfg,eleSpn) != 0) {
      if(rank==0) fprintf(stderr,"error: failed to permute BF-FSZ particle labels for check\n");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    MakeProjBFCnt(eleProjBFCnt,eleNum);
    MakeSlaterElmBF_fsz(eleNum,eleProjBFCnt);

    flag = CalculateMAll_BF_fsz(eleIdx,eleSpn,qpStart,qpEnd);
    if(size>1) {
      MPI_Allreduce(&flag,&flagRdc,1,MPI_INT,MPI_MAX,comm);
      flag = flagRdc;
    }

    loop++;
    if(loop>100) {
      if(rank==0) fprintf(stderr, "error: makeInitialSampleBF_fsz: Too many loops\n");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
  } while (flag>0);

  return 0;
}

int makeInitialSample_fsz(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,int *eleSpn,
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
    if (NExUpdatePath == 4 || NExUpdatePath == 5){ //for t-J
      for(X_mi=0;X_mi<Nsize;X_mi++) {
        si = eleSpn[X_mi];
        if(eleIdx[X_mi]== -1) {
          do {
            ri = gen_rand32()%Nsite;
          } while (eleCfg[ri+si*Nsite]!= -1 || eleCfg[ri+(1-si)*Nsite]!= -1 || LocSpn[ri]==1); // seeking empty and itinerant site
          eleCfg[ri+si*Nsite]     = X_mi;
          eleIdx[X_mi]            = ri;
        }
      }
    }else{
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
    }
    /* EleNum */
    #pragma omp parallel for default(shared) private(rsi)
    #pragma loop noalias
    for(rsi=0;rsi<nsite2;rsi++) {
      eleNum[rsi] = (eleCfg[rsi] < 0) ? 0 : 1;
    }
    
    MakeProjCnt(eleProjCnt,eleNum); // this function does not change even for fsz

    flag = CalculateMAll_fsz(eleIdx,eleSpn,qpStart,qpEnd);
    //printf("DEBUG: make4: PfM=%lf\n",creal(PfM[0]));
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

void copyFromBurnSample_fsz(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,int *eleSpn) {
  int i,n;
  const int *burnEleIdx = BurnEleIdx;// BurnEleIdx is global
//  n = Nsize + 2*Nsite + 2*Nsite + NProj+Nsite;//fsz
  n = Nsize + 2*Nsite + 2*Nsite + NProj+Nsize;//fsz
  #pragma loop noalias
  for(i=0;i<n;i++) eleIdx[i] = burnEleIdx[i]; 
  return;
}

void copyToBurnSample_fsz(const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt,const int *eleSpn) {
  int i,n;
  int *burnEleIdx = BurnEleIdx;
  //n = Nsize + 2*Nsite + 2*Nsite + NProj+Nsite;//fsz
  n = Nsize + 2*Nsite + 2*Nsite + NProj+Nsize;//fsz
  #pragma loop noalias
  for(i=0;i<n;i++) burnEleIdx[i] = eleIdx[i];
  return;
}

void saveEleConfig_fsz(const int sample, const double complex logIp,
                   const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt,const int *eleSpn) {
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
  offset = sample*nsize;
  #pragma loop noalias
  for(i=0;i<nsize;i++) EleSpn[offset+i] = eleSpn[i];
  
  x = LogProjVal(eleProjCnt);
  logSqPfFullSlater[sample] = 2.0*(x+creal(logIp));//TBC
  
  return;
}

void copyFromBurnSampleBF_fsz(int *eleIdx, int *eleCfg, int *eleNum,
                              int *eleProjCnt, int *eleProjBFCnt,
                              int *eleSpn) {
  int i;
  for(i=0;i<Nsize;i++) eleIdx[i] = BurnEleIdx[i];
  for(i=0;i<Nsite2;i++) eleCfg[i] = BurnEleCfg[i];
  for(i=0;i<Nsite2;i++) eleNum[i] = BurnEleNum[i];
  for(i=0;i<NProj;i++) eleProjCnt[i] = BurnEleProjCnt[i];
  for(i=0;i<16*Nsite*Nrange;i++) eleProjBFCnt[i] = BurnEleProjCnt[NProj+i];
  for(i=0;i<Nsize;i++) eleSpn[i] = BurnEleSpn[i];
  return;
}

void copyToBurnSampleBF_fsz(const int *eleIdx, const int *eleCfg,
                            const int *eleNum, const int *eleProjCnt,
                            const int *eleProjBFCnt, const int *eleSpn) {
  int i;
  for(i=0;i<Nsize;i++) BurnEleIdx[i] = eleIdx[i];
  for(i=0;i<Nsite2;i++) BurnEleCfg[i] = eleCfg[i];
  for(i=0;i<Nsite2;i++) BurnEleNum[i] = eleNum[i];
  for(i=0;i<NProj;i++) BurnEleProjCnt[i] = eleProjCnt[i];
  for(i=0;i<16*Nsite*Nrange;i++) BurnEleProjCnt[NProj+i] = eleProjBFCnt[i];
  for(i=0;i<Nsize;i++) BurnEleSpn[i] = eleSpn[i];
  return;
}

void saveEleConfigBF_fsz(const int sample, const double complex logIp,
                         const int *eleIdx, const int *eleCfg,
                         const int *eleNum, const int *eleProjCnt,
                         const int *eleProjBFCnt, const int *eleSpn) {
  int i,offset;
  double x;

  offset = sample*Nsize;
  for(i=0;i<Nsize;i++) EleIdx[offset+i] = eleIdx[i];
  for(i=0;i<Nsize;i++) EleSpn[offset+i] = eleSpn[i];

  offset = sample*Nsite2;
  for(i=0;i<Nsite2;i++) EleCfg[offset+i] = eleCfg[i];
  for(i=0;i<Nsite2;i++) EleNum[offset+i] = eleNum[i];

  offset = sample*NProj;
  for(i=0;i<NProj;i++) EleProjCnt[offset+i] = eleProjCnt[i];

  offset = sample*16*Nsite*Nrange;
  for(i=0;i<16*Nsite*Nrange;i++) EleProjBFCnt[offset+i] = eleProjBFCnt[i];

  x = LogProjVal(eleProjCnt);
  logSqPfFullSlater[sample] = 2.0*(x+creal(logIp));

  return;
}

//void sortEleConfig(int *eleIdx, int *eleCfg, const int *eleNum) {
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

//  return;
//}

// mi (ri,s) -> mi (rj,t)
void makeCandidate_hopping_fsz(int *mi_, int *ri_, int *rj_, int *s_,int *t_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg,const int *eleNum,const int *eleSpn) {
  const int icnt_max = Nsite*Nsite;
  int icnt;
  int mi, ri, rj, s, flag;
  int t; //fsz

  flag = 0; // FALSE
  do {
    mi = gen_rand32()%Nsize;
    s  = eleSpn[mi] ; //fsz 
    //t  = (genrand_real2()<0.5) ? s : 1-s; //fsz
    ri = eleIdx[mi];  //fsz
  } while (LocSpn[ri] == 1);

  icnt = 0;
  do {
    rj = gen_rand32()%Nsite;
    t  = (genrand_real2()<0.5) ? 0 : 1; //fsz
    if(icnt> icnt_max){
      flag = 1; // TRUE
      break;
    }
    icnt+=1;
  } while (eleCfg[rj+t*Nsite] != -1 || LocSpn[rj]==1);

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_  = s;
  *t_  = t;
  *rejectFlag_ = flag;

  return;
}
// mi (ri,s) -> mi (rj,s)
void makeCandidate_hopping_csz(int *mi_, int *ri_, int *rj_, int *s_,int *t_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg,const int *eleNum,const int *eleSpn) {
  const int icnt_max = Nsite*Nsite;
  int icnt;
  int mi, ri, rj, s, flag;
  int t; //fsz

  flag = 0; // FALSE
  do {
    mi = gen_rand32()%Nsize;
    s  = eleSpn[mi] ; //fsz 
    //t  = (genrand_real2()<0.5) ? s : 1-s; //fsz
    t  = s;//csz
    ri = eleIdx[mi];  //fsz
  } while (LocSpn[ri] == 1);

  icnt = 0;
  do {
    rj = gen_rand32()%Nsite;
    if(icnt> icnt_max){
      flag = 1; // TRUE
      break;
    }
    icnt+=1;
  } while (eleCfg[rj+t*Nsite] != -1 || LocSpn[rj]==1);

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_  = s;
  *t_  = t;
  *rejectFlag_ = flag;

  return;
}

// mi (ri,s) -> mi (rj,t)
void makeCandidate_spin_hopping_fsz(int *mi_, int *ri_, int *rj_, int *s_,int *t_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg,const int *eleNum,const int *eleSpn) {
  const int icnt_max = Nsite*Nsite;
  int icnt;
  int mi, ri, rj, s, flag;
  int t; //fsz

  flag = 0; // FALSE
  //do {
    mi = gen_rand32()%Nsize;
    s  = eleSpn[mi] ; //fsz
    //t  = (genrand_real2()<0.5) ? s : 1-s; //fsz
    ri = eleIdx[mi];  //fsz
  //} while (LocSpn[ri] == 1);

  icnt = 0;
  do {
    rj = gen_rand32()%Nsite;
    t  = (genrand_real2()<0.5) ? 0 : 1; //fsz
    if(icnt> icnt_max){
      flag = 1; // TRUE
      break;
    }
    icnt+=1;
  } while (eleCfg[rj+t*Nsite] != -1 || eleCfg[rj+(1-t)*Nsite] != -1 || LocSpn[rj]==1);

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_  = s;
  *t_  = t;
  *rejectFlag_ = flag;

  return;
}
// mi (ri,s) -> mi (rj,s)
void makeCandidate_spin_hopping_csz(int *mi_, int *ri_, int *rj_, int *s_,int *t_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg,const int *eleNum,const int *eleSpn) {
  const int icnt_max = Nsite*Nsite;
  int icnt;
  int mi, ri, rj, s, flag;
  int t; //fsz

  flag = 0; // FALSE
  //do {
    mi = gen_rand32()%Nsize;
    s  = eleSpn[mi] ; //fsz
    //t  = (genrand_real2()<0.5) ? s : 1-s; //fsz
    t  = s;//csz
    ri = eleIdx[mi];  //fsz
  //} while (LocSpn[ri] == 1);

  icnt = 0;
  do {
    rj = gen_rand32()%Nsite;
    if(icnt> icnt_max){
      flag = 1; // TRUE
      break;
    }
    icnt+=1;
  } while (eleCfg[rj+t*Nsite] != -1 || eleCfg[rj+(1-s)*Nsite] != -1 || LocSpn[rj]==1);

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_  = s;
  *t_  = t;
  *rejectFlag_ = flag;

  return;
}


/* The mi-th electron with spin s exchanges with the electron on site rj with spin 1-s */
void makeCandidate_exchange_fsz(int *mi_, int *ri_, int *rj_, int *s_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg, const int *eleNum,const int *eleSpn) {
  int mi, mj, ri, rj, s, t, flag,spn_0,spn_1;

// DEBUG!!!!!!!!!!!!!!!!!!!!!
/*
  for(mi=0;mi<Nsize;mi++){
    printf("XDEBUG: mi=%d spn=%d idx=%d\n",mi,eleSpn[mi],eleIdx[mi]);
  }
  for(ri=0;ri<Nsite;ri++){
    printf("XDEBUG: ri=%d up=%d down=%d\n",ri,eleNum[ri],eleNum[ri+Nsite]);
  }
*/
// DEBUG!!!!!!!!!!!!!!!!!!!!!

  flag = 1; // TRUE
  spn_0 = 0;//
  spn_1 = 0;//
  for(ri=0;ri<Nsite;ri++){
    if((eleNum[ri]+eleNum[ri+Nsite]) == 1  ){// up or down exists
      if(spn_0==0){
        spn_0  = 2*eleNum[ri]-1;// 0 (up)-> 1, 1(down)-> -1
      }else{
        spn_1 =  2*eleNum[ri]-1;// 0 (up)-> 1, 1(down)-> -1
      }
      //printf("ri =%d %d %d : spn0 %d spn1 %d\n",ri,eleNum[ri],eleNum[ri+Nsite],spn_0,spn_1);
      if(spn_0*spn_1<0){
        flag = 0; // FALSE
        break;
      }
    }
  }
  //printf("flag= %d spn_0=%d spn_1=%d \n",flag,spn_0,spn_1);
  if(flag) {
    *rejectFlag_ = flag;
    return;
  }

  do {
    mi = gen_rand32()%Nsize;//fsz
    s  = eleSpn[mi];// fsz //s = (genrand_real2()<0.5) ? 0 : 1;
    ri = eleIdx[mi]; //fsz
  } while (eleCfg[ri+(1-s)*Nsite] != -1);
  t = 1-s;
  do {
    mj = gen_rand32()%Nsize; //fsz
    rj = eleIdx[mj]; //fsz
  } while (eleCfg[rj+(1-t)*Nsite] != -1 || eleSpn[mj]!=t); // is it OK ?

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_ = s;
  *rejectFlag_ = flag;
  return;
}
//
// mi (ri,s) -> mi (ri,1-s) // local spin flip for conduction
void makeCandidate_LocalSpinFlip_localspin(int *mi_, int *ri_, int *rj_, int *s_,int *t_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg,const int *eleNum,const int *eleSpn) {
  int mi, ri, rj, s, flag;
  int t; //fsz

  flag = 0; // FALSE
  do {
    mi = gen_rand32()%Nsize;
    s  = eleSpn[mi] ; //fsz 
    t  = 1-s;
    //t  = (genrand_real2()<0.5) ? s : 1-s; //fsz
    ri = eleIdx[mi];  //fsz
    rj = ri;  //fsz // note ! we assume local spin
  } while (LocSpn[ri] == 0);

  //if(t==s) flag=1;

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_  = s;
  *t_  = t;
  *rejectFlag_ = flag;

  return;
}
//
// mi (ri,s) -> mi (ri,1-s) // local spin flip for conduction electrons
void makeCandidate_LocalSpinFlip_conduction(int *mi_, int *ri_, int *rj_, int *s_,int *t_, int *rejectFlag_,
                           const int *eleIdx, const int *eleCfg,const int *eleNum,const int *eleSpn) {
  const int icnt_max = Nsite*Nsite;
  int icnt=0;
  int mi, ri, rj, s, flag;
  int t; //fsz

  flag = 0; // FALSE
  do {
    mi = gen_rand32()%Nsize;
    s  = eleSpn[mi] ; //fsz 
    t  = 1-s;
    //t  = (genrand_real2()<0.5) ? s : 1-s; //fsz
    ri = eleIdx[mi];  //fsz
    rj = ri;  //fsz // note ! we assume local spin
    if(icnt> icnt_max){ // all doublons can not be accepted
      flag = 1; // TRUE
      break;
    }
    icnt+=1;
  } while (LocSpn[ri] == 1 || eleCfg[ri+t*Nsite] != -1);

  //if(t==s) flag=1;

  *mi_ = mi;
  *ri_ = ri;
  *rj_ = rj;
  *s_  = s;
  *t_  = t;
  *rejectFlag_ = flag;

  return;
}





/* The mi-th electron with spin s hops to site rj and t */
void updateEleConfig_fsz(int mi, int org_r, int dst_r, int org_spn,int dst_spn,
                     int *eleIdx, int *eleCfg, int *eleNum, int *eleSpn) {
  eleIdx[mi]         = dst_r; 
  eleSpn[mi]         = dst_spn;  //fsz 
//
  eleCfg[org_r+org_spn*Nsite] = -1;
  eleCfg[dst_r+dst_spn*Nsite] = mi;
//
  eleNum[org_r+org_spn*Nsite] = 0;
  eleNum[dst_r+dst_spn*Nsite] = 1;
  return;
}

void revertEleConfig_fsz(int mi, int org_r, int dst_r, int org_spn,int dst_spn,
                     int *eleIdx, int *eleCfg, int *eleNum,int *eleSpn) {
  eleIdx[mi]         = org_r; 
  eleSpn[mi]         = org_spn; //fsz 
//
  eleCfg[org_r+org_spn*Nsite] = mi;
  eleCfg[dst_r+dst_spn*Nsite] = -1;
//
  eleNum[org_r+org_spn*Nsite] = 1;
  eleNum[dst_r+dst_spn*Nsite] = 0;
  return;
}

int CheckEleNum_fsz(int *eleIdx, int *eleCfg, int *eleNum,int *eleSpn,MPI_Comm comm){
  int ri,si;
  int total_num;

  total_num=0;
  for(ri=0;ri<Nsite;ri++){
    for(si=0;si<2;si++){
      total_num+=eleNum[ri+si*Nsite];
    }
  }
  return total_num;
}

void CheckEleConfig_fsz(int *eleIdx, int *eleCfg, int *eleNum,int *eleSpn,MPI_Comm comm){
  int mi,ri,si;
  int check_ri,check_si;
  int rank;
  MPI_Comm_rank(comm,&rank);

  for(ri=0;ri<Nsite;ri++){
    for(si=0;si<2;si++){
      mi = eleCfg[ri+si*Nsite];
      if(mi>=0){
        check_ri = eleIdx[mi];
        check_si = eleSpn[mi];
        if(ri!=check_ri || si!=check_si){
          if(rank==0) fprintf(stderr, "error: vmcmakesample: fatal error in making sample: mi %d :ri %d %d: si %d %d\n",mi,ri,check_ri,si,check_si);
          MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
        }
      }
    }
  }
}

#endif
