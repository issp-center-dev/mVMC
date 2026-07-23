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
 * BackFlow sampling and physical quantities
 *-------------------------------------------------------------*/
#pragma once

#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "backflow_nbody.h"
#include "nbody_operator.h"
#include "vmcmake.h"
#include "vmcmake_real.h"
#include "vmccal.h"
#include "slater.h"
#include "projection.h"
#include "pfupdate_two_fcmp.h"
#include "pfupdate_real.h"
#include "matrix.h"
#include "splitloop.h"
#include "qp.h"
#include "qp_real.h"
#include "rbm.h"
#include "calham.h"
#include "calham_real.h"
#include "lslocgrn_bf.h"

#ifdef _pf_block_update
#include "../pfupdates/pf_interface.h"
#endif

static int lsbfLocalVectorFinite(void) {
  int idx;
  if(AllComplexFlag == 0) {
    for(idx=0; idx<NLSHam*NLSHam; idx++) {
      if(!isfinite(LSLQ_real[idx])) return 0;
    }
  } else {
    for(idx=0; idx<NLSHam*NLSHam; idx++) {
      if(!isfinite(creal(LSLQ[idx])) || !isfinite(cimag(LSLQ[idx]))) return 0;
    }
  }
  return 1;
}

static long lsbfNonfiniteInjectionSample(void) {
  const char *value = getenv("MVMC_LANCZOS_TEST_NONFINITE_SAMPLE");
  char *end = NULL;
  long sample;
  if(value == NULL || value[0] == '\0') return -1;
  errno = 0;
  sample = strtol(value, &end, 10);
  if(errno != 0 || end == value || *end != '\0' || sample < 0) return -2;
  return sample;
}

static void lsbfFinalizeAccounting(MPI_Comm comm, int rank,
                                   long long checkedLocal,
                                   long long rejectedLocal) {
  long long local[2] = {checkedLocal, rejectedLocal};
  long long global[2] = {checkedLocal, rejectedLocal};
#ifdef _mpi_use
  MPI_Allreduce(local, global, 2, MPI_LONG_LONG, MPI_SUM, comm);
#endif
  if(global[0] <= 0) {
    if(rank == 0) fprintf(stderr, "Error: BackFlow Lanczos checked zero samples.\n");
    MPI_Abort(comm, EXIT_FAILURE);
  }
  if(global[1] > 0 && global[1]*100 <= global[0]) {
    if(rank == 0) {
      fprintf(stderr,
              "warning: BackFlow Lanczos rejected %lld/%lld samples due to "
              "non-finite local vectors or numerical candidate rebuild failures.\n",
              global[1], global[0]);
    }
  } else if(global[1]*100 > global[0]) {
    if(rank == 0) {
      fprintf(stderr,
              "Error: BackFlow Lanczos rejected %lld/%lld samples (>1%%); biased output will not be written.\n",
              global[1], global[0]);
    }
    MPI_Abort(comm, EXIT_FAILURE);
  }
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
  int msaTmp[NQPFull * Nsize], icount[NQPFull]; // For BackFlow
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
        StartTimer(64);
        UpdateSlaterElmBF_fcmp(mi, ri, rj, s, TmpEleCfg, TmpEleNum, projBFCntNew, msaTmp, icount,
                               SlaterElmBF);
        StopTimer(64);
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
          StartTimer(64);
          UpdateSlaterElmBF_fcmp(mi, rj, ri, s, TmpEleCfg, TmpEleNum, TmpEleProjBFCnt, msaTmp, icount,
                                 SlaterElmBF);
          StopTimer(64);
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

#ifdef MVMC_DEBUG_BF_REAL_UPDATE
static void CheckBFRealUpdate(const char *label, const int *msa, const int *hopNum,
                              const int *msaRef, const int *hopNumRef, const int rank) {
  int qpidx, hop, idx;

  for (qpidx = 0; qpidx < NQPFull; qpidx++) {
    if (hopNum[qpidx] != hopNumRef[qpidx]) {
      fprintf(stderr,
              "error: BF real update mismatch (%s): rank=%d qp=%d hopNum=%d ref=%d\n",
              label, rank, qpidx, hopNum[qpidx], hopNumRef[qpidx]);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    for (hop = 0; hop < hopNum[qpidx]; hop++) {
      idx = qpidx*Nsize + hop;
      if (msa[idx] != msaRef[idx]) {
        fprintf(stderr,
                "error: BF real update mismatch (%s): rank=%d qp=%d hop=%d msa=%d ref=%d\n",
                label, rank, qpidx, hop, msa[idx], msaRef[idx]);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
    }
  }

  for (idx = 0; idx < NQPFull*Nsite2*Nsite2; idx++) {
    const double ref = creal(SlaterElmBF[idx]);
    if (SlaterElmBF_real[idx] != ref) {
      fprintf(stderr,
              "error: BF real update mismatch (%s): rank=%d idx=%d slt=%.17e ref=%.17e\n",
              label, rank, idx, SlaterElmBF_real[idx], ref);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }
}

static void UpdateSlaterElmBF_real_checked(const char *label,
                                           const int ma, const int ra, const int rb, const int u,
                                           const int *eleCfg, const int *eleNum,
                                           const int *eleProjBFCntOld, const int *eleProjBFCnt,
                                           int *msa, int *hopNum,
                                           int *msaRef, int *hopNumRef, const int rank) {
  UpdateSlaterElmBF_fcmp(ma, ra, rb, u, eleCfg, eleNum, eleProjBFCnt, msaRef, hopNumRef,
                         SlaterElmBF);
  UpdateSlaterElmBF_real(ma, ra, rb, u, eleCfg, eleNum, eleProjBFCntOld, eleProjBFCnt, msa, hopNum,
                         SlaterElmBF_real);
  CheckBFRealUpdate(label, msa, hopNum, msaRef, hopNumRef, rank);
}
#endif

void VMC_BF_MakeSample_real(MPI_Comm comm) {
  int outStep, nOutStep;
  int inStep, nInStep;
  UpdateType updateType;
  int mi, mj, ri, rj, s, t, i;
  int nAccept = 0;
  int sample;
  int tmp_i;

  double logIpOld, logIpNew; /* logarithm of inner product <phi|L|x> */ // is this ok ? TBC
  int projCntNew[NProj];
  int projBFCntNew[16 * Nsite * Nrange]; // For BackFlow
  int msaTmp[NQPFull * Nsize], icount[NQPFull]; // For BackFlow
#ifdef MVMC_DEBUG_BF_REAL_UPDATE
  int msaRef[NQPFull * Nsize], icountRef[NQPFull];
#endif
  double pfMNew_real[NQPFull];
  double x, w; // TBC x will be complex number

  int qpStart, qpEnd;
  int rejectFlag;
  int rank, size;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  SplitLoop(&qpStart, &qpEnd, NQPFull, rank, size);

  StartTimer(30);
  if (BurnFlag == 0) {
    makeInitialSampleBF_real(TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt, TmpEleProjBFCnt,
                             qpStart, qpEnd, comm);
    //makeInitialSample(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,
    //                  qpStart,qpEnd,comm);
  } else {
    //copyFromBurnSample(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt);
    copyFromBurnSampleBF(TmpEleIdx);
    MakeSlaterElmBF_fcmp(TmpEleNum, TmpEleProjBFCnt);
#pragma omp parallel for default(shared) private(tmp_i)
    for(tmp_i=0;tmp_i<NQPFull*(2*Nsite)*(2*Nsite);tmp_i++) SlaterElmBF_real[tmp_i]= creal(SlaterElmBF[tmp_i]);
  }

  CalculateMAll_BF_real(TmpEleIdx, qpStart, qpEnd);
  // printf("DEBUG: maker1: PfM=%lf\n",creal(PfM[0]));
  logIpOld = CalculateLogIP_real(PfM_real, qpStart, qpEnd, comm);
  if (!isfinite(logIpOld)) {
    if (rank == 0) fprintf(stderr, "waring: VMCMakeSample remakeSample logIpOld=%e\n", creal(logIpOld)); //TBC
    //    makeInitialSample(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt,
    //                    qpStart,qpEnd,comm);
    makeInitialSampleBF_real(TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt, TmpEleProjBFCnt,
                             qpStart, qpEnd, comm);

    CalculateMAll_BF_real(TmpEleIdx, qpStart, qpEnd);
    //printf("DEBUG: maker2: PfM=%lf\n",creal(PfM[0]));
    logIpOld = CalculateLogIP_real(PfM_real, qpStart, qpEnd, comm);
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
        AddBFProfileCounter(BFPROF_HOP_TRY, 1);

        StartTimer(31);
        makeCandidate_hopping(&mi, &ri, &rj, &s, &rejectFlag,
                              TmpEleIdx, TmpEleCfg);
        StopTimer(31);

        if (rejectFlag) {
          AddBFProfileCounter(BFPROF_HOP_CANDIDATE_REJECT, 1);
          continue;
        }
        AddBFProfileCounter(BFPROF_HOP_VALID, 1);

        StartTimer(32);
        StartTimer(60);
        /* The mi-th electron with spin s hops to site rj */
        updateEleConfig(mi, ri, rj, s, TmpEleIdx, TmpEleCfg, TmpEleNum);
        UpdateProjCnt(ri, rj, s, projCntNew, TmpEleProjCnt, TmpEleNum);
        MakeProjBFCnt(projBFCntNew, TmpEleNum);
        StopTimer(60);
        StartTimer(64);
#ifdef MVMC_DEBUG_BF_REAL_UPDATE
        UpdateSlaterElmBF_real_checked("hopping proposal", mi, ri, rj, s, TmpEleCfg, TmpEleNum,
                                       TmpEleProjBFCnt, projBFCntNew, msaTmp, icount, msaRef, icountRef, rank);
#else
        /* Real BF sampling updates SlaterElmBF_real directly; SlaterElmBF is not kept current on this hot path. */
        UpdateSlaterElmBF_real(mi, ri, rj, s, TmpEleCfg, TmpEleNum, TmpEleProjBFCnt, projBFCntNew, msaTmp, icount,
                               SlaterElmBF_real);
#endif
        StopTimer(64);

        StartTimer(61);
        //CalculateNewPfM2(mi,s,pfMNew,TmpEleIdx,qpStart,qpEnd);
        //CalculateNewPfM2_real(mi,s,pfMNew_real,TmpEleIdx,qpStart,qpEnd);
        CalculateNewPfMBF_real(icount, msaTmp, pfMNew_real, TmpEleIdx, qpStart, qpEnd, SlaterElmBF_real);

        //printf("DEBUG: out %d in %d pfMNew=%lf \n",outStep,inStep,creal(pfMNew[0]));
        StopTimer(61);

        StartTimer(62);
        /* calculate inner product <phi|L|x> */
        //logIpNew = CalculateLogIP_fcmp(pfMNew,qpStart,qpEnd,comm);
        logIpNew = CalculateLogIP_real(pfMNew_real, qpStart, qpEnd, comm);
        StopTimer(62);

        /* Metroplis */
        x = LogProjRatio(projCntNew, TmpEleProjCnt);
        w = exp(2.0 * (x + (logIpNew - logIpOld)));
        if (!isfinite(w)) w = -1.0; /* should be rejected */

        if (w > genrand_real2()) { /* accept */
          // UpdateMAll will change SlaterElm, InvM (including PfM)
          StartTimer(63);
          UpdateMAll_BF_real(icount, msaTmp, PfM_real, TmpEleIdx, qpStart, qpEnd);
          //          UpdateMAll_real(mi,s,TmpEleIdx,qpStart,qpEnd);
          //            UpdateMAll(mi,s,TmpEleIdx,qpStart,qpEnd);
          StopTimer(63);

          for (i = 0; i < NProj; i++) TmpEleProjCnt[i] = projCntNew[i];
          for (i = 0; i < 16 * Nsite * Nrange; i++) TmpEleProjBFCnt[i] = projBFCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          Counter[1]++;
          AddBFProfileCounter(BFPROF_HOP_ACCEPT, 1);
        } else { /* reject */
          AddBFProfileCounter(BFPROF_HOP_METROPOLIS_REJECT, 1);
          revertEleConfig(mi, ri, rj, s, TmpEleIdx, TmpEleCfg, TmpEleNum);
          StartTimer(64);
#ifdef MVMC_DEBUG_BF_REAL_UPDATE
          UpdateSlaterElmBF_real_checked("hopping reject", mi, rj, ri, s, TmpEleCfg, TmpEleNum,
                                         projBFCntNew, TmpEleProjBFCnt, msaTmp, icount, msaRef, icountRef, rank);
#else
          UpdateSlaterElmBF_real(mi, rj, ri, s, TmpEleCfg, TmpEleNum, projBFCntNew, TmpEleProjBFCnt, msaTmp, icount,
                                 SlaterElmBF_real);
#endif
          StopTimer(64);
        }
        StopTimer(32);

      } else if (updateType == EXCHANGE) { /* exchange */
        Counter[2]++;
        AddBFProfileCounter(BFPROF_EXCHANGE_TRY, 1);

        StartTimer(31);
        makeCandidate_exchange(&mi, &ri, &rj, &s, &rejectFlag,
                               TmpEleIdx, TmpEleCfg, TmpEleNum);
        StopTimer(31);

        if (rejectFlag) {
          AddBFProfileCounter(BFPROF_EXCHANGE_CANDIDATE_REJECT, 1);
          continue;
        }
        AddBFProfileCounter(BFPROF_EXCHANGE_VALID, 1);

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

        CalculateNewPfMTwo2_real(mi, s, mj, t, pfMNew_real, TmpEleIdx, qpStart, qpEnd);
        StopTimer(66);
        StartTimer(67);

        /* calculate inner product <phi|L|x> */
        logIpNew = CalculateLogIP_real(pfMNew_real, qpStart, qpEnd, comm);

        StopTimer(67);

        /* Metroplis */
        x = LogProjRatio(projCntNew, TmpEleProjCnt);
        w = exp(2.0 * (x + (logIpNew - logIpOld))); //TBC
        if (!isfinite(w)) w = -1.0; /* should be rejected */

        if (w > genrand_real2()) { /* accept */
          StartTimer(68);
          UpdateMAllTwo_real(mi, s, mj, t, ri, rj, TmpEleIdx, qpStart, qpEnd);
          StopTimer(68);

          for (i = 0; i < NProj; i++) TmpEleProjCnt[i] = projCntNew[i];
          logIpOld = logIpNew;
          nAccept++;
          Counter[3]++;
          AddBFProfileCounter(BFPROF_EXCHANGE_ACCEPT, 1);
        } else { /* reject */
          AddBFProfileCounter(BFPROF_EXCHANGE_METROPOLIS_REJECT, 1);
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
        CalculateMAll_BF_real(TmpEleIdx, qpStart, qpEnd);
        logIpOld = CalculateLogIP_real(PfM_real, qpStart, qpEnd, comm);
        StopTimer(34);
        nAccept = 0;
      }
    } /* end of instep */

    StartTimer(35);
    /* save Electron Configuration */
    if (outStep >= nOutStep - NVMCSample) {
      sample = outStep - (nOutStep - NVMCSample);
      saveEleConfigBF(sample, logIpOld, TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt, TmpEleProjBFCnt);
      //      saveEleConfig(sample,logIpOld,TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt);
    }
    StopTimer(35);

  } /* end of outstep */

  //copyToBurnSample(TmpEleIdx,TmpEleCfg,TmpEleNum,TmpEleProjCnt);
  copyToBurnSampleBF(TmpEleIdx);
  BurnFlag = 1;
  return;
}

int makeInitialSampleBF_real(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt, int *eleProjBFCnt,
                             const int qpStart, const int qpEnd, MPI_Comm comm) {
  const int nsize = Nsize;
  const int nsite2 = Nsite2;
  int flag = 1, flagRdc, loop = 0;
  int ri, mi, si, msi, rsi;
  int rank, size;
  int tmp_i;

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
#pragma omp parallel for default(shared) private(tmp_i)
    for(tmp_i=0;tmp_i<NQPFull*(2*Nsite)*(2*Nsite);tmp_i++) SlaterElmBF_real[tmp_i]= creal(SlaterElmBF[tmp_i]);

    flag = CalculateMAll_BF_real(eleIdx, qpStart, qpEnd);
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


static void dumpBFIdentityCheck(const char *path, const int *eleIdx,
                                const int qpStart, const int qpEnd) {
  FILE *fp;
  double complex *pfNoBF;
  double complex *pfBF;
  double maxSlaterDiff = 0.0;
  double maxPfDiff = 0.0;
  int nanCount = 0;
  int infoNoBF;
  int infoBF;
  int i;
  int totalSlater = NQPFull * Nsite2 * Nsite2;

  if (path == NULL || path[0] == '\0') return;

  pfNoBF = (double complex *)malloc(sizeof(double complex) * (size_t)NQPFull);
  pfBF = (double complex *)malloc(sizeof(double complex) * (size_t)NQPFull);
  if (pfNoBF == NULL || pfBF == NULL) {
    fprintf(stderr, "Error: memory allocation failed for BackFlow identity dump.\n");
    free(pfNoBF);
    free(pfBF);
    return;
  }

  for (i = 0; i < totalSlater; i++) {
    double diff = cabs(SlaterElm[i] - SlaterElmBF[i]);
    if (!isfinite(diff)) {
      nanCount++;
    } else if (diff > maxSlaterDiff) {
      maxSlaterDiff = diff;
    }
  }

  infoNoBF = CalculateMAll_fcmp(eleIdx, qpStart, qpEnd);
  for (i = 0; i < NQPFull; i++) pfNoBF[i] = PfM[i];

  infoBF = CalculateMAll_BF_fcmp(eleIdx, qpStart, qpEnd);
  for (i = 0; i < NQPFull; i++) pfBF[i] = PfM[i];

  for (i = 0; i < NQPFull; i++) {
    double diff = cabs(pfNoBF[i] - pfBF[i]);
    if (!isfinite(diff)) {
      nanCount++;
    } else if (diff > maxPfDiff) {
      maxPfDiff = diff;
    }
  }

  fp = fopen(path, "w");
  if (fp == NULL) {
    fprintf(stderr, "Error: failed to open BackFlow identity dump file: %s\n", path);
  } else {
    fprintf(fp, "info_no_bf %d\n", infoNoBF);
    fprintf(fp, "info_bf %d\n", infoBF);
    fprintf(fp, "nqp_full %d\n", NQPFull);
    fprintf(fp, "nsite2 %d\n", Nsite2);
    fprintf(fp, "max_abs_slater_diff %.17e\n", maxSlaterDiff);
    fprintf(fp, "max_abs_pf_diff %.17e\n", maxPfDiff);
    fprintf(fp, "nan_count %d\n", nanCount);
    fclose(fp);
  }

  free(pfNoBF);
  free(pfBF);
}

static void dumpBFSlaterDiffCheck(const char *path, int *eleIdx, int *eleNum,
                                  int *eleCfg, int *eleProjCnt,
                                  const int *eleProjBFCnt,
                                  const int qpStart, const int qpEnd,
                                  const int sample, const int append) {
  FILE *fp;
  double complex *diffNoBF;
  double complex *diffBF;
  double complex *diffBFMainOrder;
  double complex ipNoBF;
  double complex ipBF;
  double complex eNoBF = 0.0 + 0.0*I;
  double complex eBF = 0.0 + 0.0*I;
  double maxDiff = 0.0;
  double maxMainOrderDiff = 0.0;
  int nanCount = 0;
  int maxIdx = -1;
  int maxMainOrderIdx = -1;
  int infoNoBF;
  int infoBF;
  int i;

  if (path == NULL || path[0] == '\0') return;

  diffNoBF = (double complex *)calloc((size_t)2 * (size_t)NSlater, sizeof(double complex));
  diffBF = (double complex *)calloc((size_t)2 * (size_t)NSlater, sizeof(double complex));
  diffBFMainOrder = (double complex *)calloc((size_t)2 * (size_t)(NProjBF + NSlater), sizeof(double complex));
  if (diffNoBF == NULL || diffBF == NULL || diffBFMainOrder == NULL) {
    fprintf(stderr, "Error: memory allocation failed for BackFlow diff dump.\n");
    free(diffNoBF);
    free(diffBF);
    free(diffBFMainOrder);
    return;
  }

  infoNoBF = CalculateMAll_fcmp(eleIdx, qpStart, qpEnd);
  ipNoBF = CalculateIP_fcmp(PfM, qpStart, qpEnd, MPI_COMM_SELF);
  if (infoNoBF == 0) {
    SlaterElmDiff_fcmp(diffNoBF, ipNoBF, eleIdx);
    eNoBF = CalculateHamiltonian(ipNoBF, eleIdx, eleCfg, eleNum, eleProjCnt, NULL);
  }

  infoBF = CalculateMAll_BF_fcmp(eleIdx, qpStart, qpEnd);
  ipBF = CalculateIP_fcmp(PfM, qpStart, qpEnd, MPI_COMM_SELF);
  if (infoBF == 0) {
    SlaterElmBFDiff_fcmp(diffBF, ipBF, eleIdx, eleNum, eleCfg, eleProjCnt, eleProjBFCnt);
    BackFlowDiff_fcmp(diffBFMainOrder, ipBF, eleIdx, eleNum, eleProjCnt, eleProjBFCnt);
    SlaterElmBFDiff_fcmp(diffBFMainOrder + 2 * NProjBF, ipBF, eleIdx, eleNum, eleCfg, eleProjCnt, eleProjBFCnt);
    eBF = CalculateHamiltonianBF_fcmp(ipBF, eleIdx, eleCfg, eleNum, eleProjCnt, eleProjBFCnt);
  }

  for (i = 0; i < 2 * NSlater; i++) {
    double diff = cabs(diffBF[i] - diffNoBF[i]);
    if (!isfinite(diff)) {
      nanCount++;
    } else if (diff > maxDiff) {
      maxDiff = diff;
      maxIdx = i;
    }
    diff = cabs(diffBFMainOrder[2 * NProjBF + i] - diffNoBF[i]);
    if (!isfinite(diff)) {
      nanCount++;
    } else if (diff > maxMainOrderDiff) {
      maxMainOrderDiff = diff;
      maxMainOrderIdx = i;
    }
  }

  fp = fopen(path, append ? "a" : "w");
  if (fp == NULL) {
    fprintf(stderr, "Error: failed to open BackFlow diff dump file: %s\n", path);
  } else {
    fprintf(fp, "sample %d\n", sample);
    fprintf(fp, "info_no_bf %d\n", infoNoBF);
    fprintf(fp, "info_bf %d\n", infoBF);
    fprintf(fp, "nslater %d\n", NSlater);
    fprintf(fp, "max_abs_slater_diff_o %.17e\n", maxDiff);
    fprintf(fp, "max_idx %d\n", maxIdx);
    fprintf(fp, "max_abs_main_order_slater_diff_o %.17e\n", maxMainOrderDiff);
    fprintf(fp, "max_main_order_idx %d\n", maxMainOrderIdx);
    fprintf(fp, "nan_count %d\n", nanCount);
    fprintf(fp, "e_no_bf %.17e %.17e\n", creal(eNoBF), cimag(eNoBF));
    fprintf(fp, "e_bf %.17e %.17e\n", creal(eBF), cimag(eBF));
    fprintf(fp, "abs_energy_diff %.17e\n", cabs(eBF - eNoBF));
    if (maxIdx >= 0) {
      fprintf(fp, "no_bf_at_max %.17e %.17e\n",
              creal(diffNoBF[maxIdx]), cimag(diffNoBF[maxIdx]));
      fprintf(fp, "bf_at_max %.17e %.17e\n",
              creal(diffBF[maxIdx]), cimag(diffBF[maxIdx]));
    }
    fprintf(fp, "\n");
    fclose(fp);
  }

  free(diffNoBF);
  free(diffBF);
  free(diffBFMainOrder);
}

static int calculateBFIPForFD(int *eleIdx, int *eleNum, const int *eleProjBFCnt,
                              const int qpStart, const int qpEnd,
                              double complex *ip) {
  int info;

  MakeSlaterElmBF_fcmp(eleNum, eleProjBFCnt);
  info = CalculateMAll_BF_fcmp(eleIdx, qpStart, qpEnd);
  *ip = CalculateIP_fcmp(PfM, qpStart, qpEnd, MPI_COMM_SELF);
  return info;
}

static void dumpBFProjBFFiniteDiffCheck(const char *path, int *eleIdx,
                                        int *eleNum, int *eleProjCnt,
                                        const int *eleProjBFCnt,
                                        const int qpStart, const int qpEnd,
                                        const int sample, const int append) {
  FILE *fp;
  double complex *analytic;
  double complex *projBFStore;
  double complex ip0 = 0.0 + 0.0*I;
  double complex ipPlus = 0.0 + 0.0*I;
  double complex ipMinus = 0.0 + 0.0*I;
  double complex fd;
  double complex analyticAtMax = 0.0 + 0.0*I;
  double complex fdAtMax = 0.0 + 0.0*I;
  const double h = 1.0e-6;
  double maxRealDiff = 0.0;
  double maxImagDiff = 0.0;
  double maxDiff = 0.0;
  double maxFDAbs = 0.0;
  int maxRealIdx = -1;
  int maxImagIdx = -1;
  int maxIdx = -1;
  int maxIsImag = 0;
  int nonzeroFDCount = 0;
  int nanCount = 0;
  int fdFailCount = 0;
  int infoBase;
  int infoPlus;
  int infoMinus;
  int idx;

  if (path == NULL || path[0] == '\0') return;
  if (NProjBF <= 0 || ProjBF == NULL) return;

  analytic = (double complex *)calloc((size_t)2 * (size_t)NProjBF, sizeof(double complex));
  projBFStore = (double complex *)calloc((size_t)NProjBF, sizeof(double complex));
  if (analytic == NULL || projBFStore == NULL) {
    fprintf(stderr, "Error: memory allocation failed for BackFlow finite-difference dump.\n");
    free(analytic);
    free(projBFStore);
    return;
  }

  for (idx = 0; idx < NProjBF; idx++) projBFStore[idx] = ProjBF[idx];

  infoBase = calculateBFIPForFD(eleIdx, eleNum, eleProjBFCnt, qpStart, qpEnd, &ip0);
  if (infoBase == 0 && isfinite(creal(ip0)) && isfinite(cimag(ip0)) && cabs(ip0) > 0.0) {
    BackFlowDiff_fcmp(analytic, ip0, eleIdx, eleNum, eleProjCnt, eleProjBFCnt);
    for (idx = 0; idx < NProjBF; idx++) {
      double diff;

      ProjBF[idx] = projBFStore[idx] + h;
      infoPlus = calculateBFIPForFD(eleIdx, eleNum, eleProjBFCnt, qpStart, qpEnd, &ipPlus);
      ProjBF[idx] = projBFStore[idx] - h;
      infoMinus = calculateBFIPForFD(eleIdx, eleNum, eleProjBFCnt, qpStart, qpEnd, &ipMinus);
      ProjBF[idx] = projBFStore[idx];
      if (infoPlus != 0 || infoMinus != 0 ||
          !isfinite(creal(ipPlus)) || !isfinite(cimag(ipPlus)) ||
          !isfinite(creal(ipMinus)) || !isfinite(cimag(ipMinus))) {
        fdFailCount++;
      } else {
        fd = ((ipPlus - ipMinus) / (2.0 * h)) / ip0;
        diff = cabs(fd - analytic[2 * idx]);
        if (!isfinite(diff) || !isfinite(cabs(fd))) {
          nanCount++;
        } else {
          if (cabs(fd) > 1.0e-12) nonzeroFDCount++;
          if (cabs(fd) > maxFDAbs) maxFDAbs = cabs(fd);
          if (diff > maxRealDiff) {
            maxRealDiff = diff;
            maxRealIdx = idx;
          }
          if (diff > maxDiff) {
            maxDiff = diff;
            maxIdx = idx;
            maxIsImag = 0;
            analyticAtMax = analytic[2 * idx];
            fdAtMax = fd;
          }
        }
      }

      if (idx == 0) continue; /* ProjBF[0] imaginary part is fixed by the input contract. */
      ProjBF[idx] = projBFStore[idx] + I * h;
      infoPlus = calculateBFIPForFD(eleIdx, eleNum, eleProjBFCnt, qpStart, qpEnd, &ipPlus);
      ProjBF[idx] = projBFStore[idx] - I * h;
      infoMinus = calculateBFIPForFD(eleIdx, eleNum, eleProjBFCnt, qpStart, qpEnd, &ipMinus);
      ProjBF[idx] = projBFStore[idx];
      if (infoPlus != 0 || infoMinus != 0 ||
          !isfinite(creal(ipPlus)) || !isfinite(cimag(ipPlus)) ||
          !isfinite(creal(ipMinus)) || !isfinite(cimag(ipMinus))) {
        fdFailCount++;
      } else {
        fd = ((ipPlus - ipMinus) / (2.0 * h)) / ip0;
        diff = cabs(fd - analytic[2 * idx + 1]);
        if (!isfinite(diff) || !isfinite(cabs(fd))) {
          nanCount++;
        } else {
          if (cabs(fd) > 1.0e-12) nonzeroFDCount++;
          if (cabs(fd) > maxFDAbs) maxFDAbs = cabs(fd);
          if (diff > maxImagDiff) {
            maxImagDiff = diff;
            maxImagIdx = idx;
          }
          if (diff > maxDiff) {
            maxDiff = diff;
            maxIdx = idx;
            maxIsImag = 1;
            analyticAtMax = analytic[2 * idx + 1];
            fdAtMax = fd;
          }
        }
      }
    }
  } else {
    fdFailCount++;
  }

  for (idx = 0; idx < NProjBF; idx++) ProjBF[idx] = projBFStore[idx];
  MakeSlaterElmBF_fcmp(eleNum, eleProjBFCnt);

  fp = fopen(path, append ? "a" : "w");
  if (fp == NULL) {
    fprintf(stderr, "Error: failed to open BackFlow finite-difference dump file: %s\n", path);
  } else {
    fprintf(fp, "sample %d\n", sample);
    fprintf(fp, "info_base %d\n", infoBase);
    fprintf(fp, "nprojbf %d\n", NProjBF);
    fprintf(fp, "step %.17e\n", h);
    fprintf(fp, "fd_fail_count %d\n", fdFailCount);
    fprintf(fp, "nan_count %d\n", nanCount);
    fprintf(fp, "nonzero_fd_count %d\n", nonzeroFDCount);
    fprintf(fp, "max_abs_fd_value %.17e\n", maxFDAbs);
    fprintf(fp, "max_abs_projbf_fd_real %.17e\n", maxRealDiff);
    fprintf(fp, "max_real_idx %d\n", maxRealIdx);
    fprintf(fp, "max_abs_projbf_fd_imag %.17e\n", maxImagDiff);
    fprintf(fp, "max_imag_idx %d\n", maxImagIdx);
    fprintf(fp, "max_abs_projbf_fd_diff %.17e\n", maxDiff);
    fprintf(fp, "max_idx %d\n", maxIdx);
    fprintf(fp, "max_is_imag %d\n", maxIsImag);
    fprintf(fp, "analytic_at_max %.17e %.17e\n", creal(analyticAtMax), cimag(analyticAtMax));
    fprintf(fp, "fd_at_max %.17e %.17e\n", creal(fdAtMax), cimag(fdAtMax));
    fprintf(fp, "\n");
    fclose(fp);
  }

  free(analytic);
  free(projBFStore);
}

static void copySlaterElmBFToReal(void) {
  int idx;
  const int n = NQPFull * Nsite2 * Nsite2;
  for (idx = 0; idx < n; idx++) SlaterElmBF_real[idx] = creal(SlaterElmBF[idx]);
}

static void copyIntArray(int *dst, const int *src, const int n) {
  int idx;
  for (idx = 0; idx < n; idx++) dst[idx] = src[idx];
}

static int isGreen2BFGeneralCase(const int ri, const int rj, const int rk,
                                 const int rl, const int s, const int t) {
  if (s == t) {
    return !(rk == rl || rj == rl || ri == rl ||
             rj == rk || ri == rk || ri == rj);
  }
  return !(rk == rl || ri == rj);
}

static int setupGreen2BFGeneralHop(const int ri, const int rj, const int rk,
                                   const int rl, const int s, const int t,
                                   int *eleIdx, int *eleCfg, int *eleNum,
                                   const int *eleProjCnt, int *projCntNew,
                                   int *projBFCntNew, double *projRatio) {
  int mj, msj, ml, mtl;
  int rsi, rsj, rtk, rtl;

  if (!isGreen2BFGeneralCase(ri, rj, rk, rl, s, t)) return 0;

  rsi = ri + s * Nsite;
  rsj = rj + s * Nsite;
  rtk = rk + t * Nsite;
  rtl = rl + t * Nsite;

  if (eleNum[rsi] == 1 || eleNum[rsj] == 0 ||
      eleNum[rtk] == 1 || eleNum[rtl] == 0) {
    return 0;
  }

  mj = eleCfg[rsj];
  ml = eleCfg[rtl];
  if (mj < 0 || ml < 0) return 0;
  msj = mj + s * Ne;
  mtl = ml + t * Ne;

  eleCfg[rtl] = -1;
  eleCfg[rtk] = ml;
  eleIdx[mtl] = rk;
  eleNum[rtl] = 0;
  eleNum[rtk] = 1;
  UpdateProjCnt(rl, rk, t, projCntNew, eleProjCnt, eleNum);

  eleCfg[rsj] = -1;
  eleCfg[rsi] = mj;
  eleIdx[msj] = ri;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt(rj, ri, s, projCntNew, projCntNew, eleNum);

  if (!IsSectorStateAllowed(eleNum)) return 0;

  *projRatio = ProjRatio(projCntNew, eleProjCnt);
  MakeProjBFCnt(projBFCntNew, eleNum);
  return 1;
}

static void dumpBFGreen2BruteForceCheck(const char *path, int *eleIdx,
                                        int *eleCfg, int *eleNum,
                                        const int *eleProjCnt,
                                        const int *eleProjBFCnt,
                                        const int qpStart, const int qpEnd,
                                        const int sample,
                                        const int append) {
  FILE *fp;
  int *idxCopy;
  int *cfgCopy;
  int *numCopy;
  int *projCntNew;
  int *projBFCntNew;
  double complex *buffer;
  double complex *sltBFTmp;
  double *bufferReal;
  double *sltBFTmpReal;
  double maxDiff = 0.0;
  double maxBrute = 0.0;
  double complex fastAtMax = 0.0 + 0.0 * I;
  double complex bruteAtMax = 0.0 + 0.0 * I;
  int maxIdx = -1;
  int compared = 0;
  int nonzeroBrute = 0;
  int skipped = 0;
  int infoFail = 0;
  int nanCount = 0;
  int idx;
  const int projSize = (NProj > 0) ? NProj : 1;
  const int bfCntSize = 16 * Nsite * Nrange;
  const int slaterSize = NQPFull * Nsite2 * Nsite2;

  if (path == NULL || path[0] == '\0') return;
  if (NCisAjsCktAltDC <= 0) return;

  idxCopy = (int *)malloc(sizeof(int) * (size_t)Nsize);
  cfgCopy = (int *)malloc(sizeof(int) * (size_t)Nsite2);
  numCopy = (int *)malloc(sizeof(int) * (size_t)Nsite2);
  projCntNew = (int *)malloc(sizeof(int) * (size_t)projSize);
  projBFCntNew = (int *)malloc(sizeof(int) * (size_t)bfCntSize);
  buffer = (double complex *)malloc(sizeof(double complex) * (size_t)(NQPFull + 2 * Nsize));
  sltBFTmp = (double complex *)malloc(sizeof(double complex) * (size_t)slaterSize);
  bufferReal = (double *)malloc(sizeof(double) * (size_t)(NQPFull + 2 * Nsize));
  sltBFTmpReal = (double *)malloc(sizeof(double) * (size_t)slaterSize);
  if (idxCopy == NULL || cfgCopy == NULL || numCopy == NULL ||
      projCntNew == NULL || projBFCntNew == NULL || buffer == NULL ||
      sltBFTmp == NULL || bufferReal == NULL || sltBFTmpReal == NULL) {
    fprintf(stderr, "Error: memory allocation failed for BackFlow GreenFunc2 dump.\n");
    free(idxCopy);
    free(cfgCopy);
    free(numCopy);
    free(projCntNew);
    free(projBFCntNew);
    free(buffer);
    free(sltBFTmp);
    free(bufferReal);
    free(sltBFTmpReal);
    return;
  }

  for (idx = 0; idx < NCisAjsCktAltDC; idx++) {
    const int ri = CisAjsCktAltDCIdx[idx][0];
    const int rj = CisAjsCktAltDCIdx[idx][2];
    const int s  = CisAjsCktAltDCIdx[idx][1];
    const int rk = CisAjsCktAltDCIdx[idx][4];
    const int rl = CisAjsCktAltDCIdx[idx][6];
    const int t  = CisAjsCktAltDCIdx[idx][5];
    double projRatio;
    double diff;
    double bruteAbs;
    int info;
    double complex ip;
    double complex fast;
    double complex brute;

    if (!isGreen2BFGeneralCase(ri, rj, rk, rl, s, t)) {
      skipped++;
      continue;
    }

    copyIntArray(idxCopy, eleIdx, Nsize);
    copyIntArray(cfgCopy, eleCfg, Nsite2);
    copyIntArray(numCopy, eleNum, Nsite2);
    MakeSlaterElmBF_fcmp(eleNum, eleProjBFCnt);
    if (AllComplexFlag == 0) {
      double ipReal;
      copySlaterElmBFToReal();
      info = CalculateMAll_BF_real(idxCopy, qpStart, qpEnd);
      ipReal = CalculateIP_real(PfM_real, qpStart, qpEnd, MPI_COMM_SELF);
      ip = ipReal + 0.0 * I;
    } else {
      info = CalculateMAll_BF_fcmp(idxCopy, qpStart, qpEnd);
      ip = CalculateIP_fcmp(PfM, qpStart, qpEnd, MPI_COMM_SELF);
    }
    if (info != 0 || !isfinite(creal(ip)) || !isfinite(cimag(ip)) || cabs(ip) <= 0.0) {
      infoFail++;
      continue;
    }
    if (AllComplexFlag == 0) {
      StoreSlaterElmBF_real(sltBFTmpReal);
      fast = GreenFunc2BF_real(ri, rj, rk, rl, s, t, creal(ip), sltBFTmpReal,
                               idxCopy, cfgCopy, numCopy, eleProjCnt,
                               projCntNew, eleProjBFCnt, projBFCntNew,
                               bufferReal) + 0.0 * I;
    } else {
      StoreSlaterElmBF_fcmp(sltBFTmp);
      fast = GreenFunc2BF(ri, rj, rk, rl, s, t, ip, sltBFTmp,
                          idxCopy, cfgCopy, numCopy, eleProjCnt,
                          projCntNew, eleProjBFCnt, projBFCntNew, buffer);
    }

    copyIntArray(idxCopy, eleIdx, Nsize);
    copyIntArray(cfgCopy, eleCfg, Nsite2);
    copyIntArray(numCopy, eleNum, Nsite2);
    if (!setupGreen2BFGeneralHop(ri, rj, rk, rl, s, t, idxCopy, cfgCopy,
                                 numCopy, eleProjCnt, projCntNew,
                                 projBFCntNew, &projRatio)) {
      skipped++;
      continue;
    }

    MakeSlaterElmBF_fcmp(numCopy, projBFCntNew);
    if (AllComplexFlag == 0) {
      double ipNewReal;
      copySlaterElmBFToReal();
      info = CalculateMAll_BF_real(idxCopy, qpStart, qpEnd);
      ipNewReal = CalculateIP_real(PfM_real, qpStart, qpEnd, MPI_COMM_SELF);
      brute = (projRatio * ipNewReal / creal(ip)) + 0.0 * I;
    } else {
      double complex ipNew;
      info = CalculateMAll_BF_fcmp(idxCopy, qpStart, qpEnd);
      ipNew = CalculateIP_fcmp(PfM, qpStart, qpEnd, MPI_COMM_SELF);
      brute = conj((projRatio * ipNew) / ip);
    }
    if (info != 0) {
      infoFail++;
      continue;
    }

    diff = cabs(fast - brute);
    bruteAbs = cabs(brute);
    if (!isfinite(diff) || !isfinite(bruteAbs)) {
      nanCount++;
      continue;
    }
    compared++;
    if (bruteAbs > 1.0e-12) nonzeroBrute++;
    if (bruteAbs > maxBrute) maxBrute = bruteAbs;
    if (diff > maxDiff) {
      maxDiff = diff;
      maxIdx = idx;
      fastAtMax = fast;
      bruteAtMax = brute;
    }
  }

  MakeSlaterElmBF_fcmp(eleNum, eleProjBFCnt);
  if (AllComplexFlag == 0) copySlaterElmBFToReal();

  fp = fopen(path, append ? "a" : "w");
  if (fp == NULL) {
    fprintf(stderr, "Error: failed to open BackFlow GreenFunc2 dump file: %s\n", path);
  } else {
    fprintf(fp, "sample %d\n", sample);
    fprintf(fp, "all_complex_flag %d\n", AllComplexFlag);
    fprintf(fp, "n_twobodyg %d\n", NCisAjsCktAltDC);
    fprintf(fp, "compared_count %d\n", compared);
    fprintf(fp, "nonzero_bruteforce_count %d\n", nonzeroBrute);
    fprintf(fp, "skipped_count %d\n", skipped);
    fprintf(fp, "info_fail_count %d\n", infoFail);
    fprintf(fp, "nan_count %d\n", nanCount);
    fprintf(fp, "max_abs_bruteforce %.17e\n", maxBrute);
    fprintf(fp, "max_abs_green2_diff %.17e\n", maxDiff);
    fprintf(fp, "max_idx %d\n", maxIdx);
    fprintf(fp, "fast_at_max %.17e %.17e\n", creal(fastAtMax), cimag(fastAtMax));
    fprintf(fp, "bruteforce_at_max %.17e %.17e\n", creal(bruteAtMax), cimag(bruteAtMax));
    fprintf(fp, "\n");
    fclose(fp);
  }

  free(idxCopy);
  free(cfgCopy);
  free(numCopy);
  free(projCntNew);
  free(projBFCntNew);
  free(buffer);
  free(sltBFTmp);
  free(bufferReal);
  free(sltBFTmpReal);
}

static int dumpBFNBodyComponentCheck(const char *path,
                                     const int *eleIdx,
                                     const int *eleNum,
                                     const int *eleProjBFCnt,
                                     int parentRank,
                                     int childRank) {
  const size_t slaterCount = (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
  const size_t invPfCount =
      (size_t)NQPFull*((size_t)Nsize*(size_t)Nsize + 1u);
  const size_t etaCount = (size_t)Nsite*(size_t)Nsite;
  double complex *candidate = NULL;
  double complex *zeroSlater = NULL;
  double complex *slaterBefore = NULL;
  double *slaterRealBefore = NULL;
  double complex *invPfBefore = NULL;
  double complex *invPfAfterReference = NULL;
  double complex *etaBefore = NULL;
  double complex *pfCandidate = NULL;
  double complex *pfReference = NULL;
  double complex *pfZero = NULL;
  double complex *pfBufM = NULL;
  double complex *pfWork = NULL;
  int *etaFlagBefore = NULL;
  int *eleSpn = NULL;
  int *pfIWork = NULL;
  double *pfRWork = NULL;
  double slaterMaxDiff = 0.0;
  double pfMaxDiff = 0.0;
  int slaterGlobalStateChanged;
  int etaStateChanged = 0;
  int etaFlagStateChanged = 0;
  int pfReferenceStatus;
  int pfStatus;
  int pfFailureDetail = -1;
  int pfGlobalStateChanged;
  int pfZeroStatus;
  int pfZeroCandidateNonzero = 0;
  int pfInvalidNullStatus;
  int pfInvalidRangeStatus;
  int pfInvalidWorkspaceStatus;
  int pfInvalidFailureDetail = 0;
  int invalidDetail;
  size_t idx;
  int ri;
  FILE *fp = NULL;
  int status = -1;

  candidate = (double complex *)malloc(slaterCount*sizeof(double complex));
  zeroSlater = (double complex *)calloc(slaterCount, sizeof(double complex));
  slaterBefore = (double complex *)malloc(slaterCount*sizeof(double complex));
  slaterRealBefore = (double *)malloc(slaterCount*sizeof(double));
  invPfBefore = (double complex *)malloc(invPfCount*sizeof(double complex));
  invPfAfterReference =
      (double complex *)malloc(invPfCount*sizeof(double complex));
  etaBefore = (double complex *)malloc(etaCount*sizeof(double complex));
  etaFlagBefore = (int *)malloc(etaCount*sizeof(int));
  pfCandidate =
      (double complex *)malloc((size_t)NQPFull*sizeof(double complex));
  pfReference =
      (double complex *)malloc((size_t)NQPFull*sizeof(double complex));
  pfZero = (double complex *)malloc((size_t)NQPFull*sizeof(double complex));
  pfBufM = (double complex *)malloc(
      (size_t)Nsize*(size_t)Nsize*sizeof(double complex));
  pfWork =
      (double complex *)malloc((size_t)LapackLWork*sizeof(double complex));
  eleSpn = (int *)malloc((size_t)Nsize*sizeof(int));
  pfIWork = (int *)malloc((size_t)Nsize*sizeof(int));
  pfRWork = (double *)malloc((size_t)LapackLWork*sizeof(double));
  if(candidate == NULL || zeroSlater == NULL || slaterBefore == NULL
     || slaterRealBefore == NULL || invPfBefore == NULL
     || invPfAfterReference == NULL || etaBefore == NULL
     || etaFlagBefore == NULL || pfCandidate == NULL || pfReference == NULL
     || pfZero == NULL || pfBufM == NULL || pfWork == NULL || eleSpn == NULL
     || pfIWork == NULL || pfRWork == NULL) {
    fprintf(stderr,
            "Error: failed to allocate BackFlow N-body component diagnostic buffers.\n");
    goto cleanup;
  }

  memcpy(slaterBefore, SlaterElmBF, slaterCount*sizeof(double complex));
  memcpy(slaterRealBefore, SlaterElmBF_real, slaterCount*sizeof(double));
  memcpy(invPfBefore, InvM, invPfCount*sizeof(double complex));
  for(ri=0;ri<Nsite;ri++) {
    memcpy(etaBefore + (size_t)ri*(size_t)Nsite, eta[ri],
           (size_t)Nsite*sizeof(double complex));
    memcpy(etaFlagBefore + (size_t)ri*(size_t)Nsite, etaFlag[ri],
           (size_t)Nsite*sizeof(int));
  }

  MakeSlaterElmBF_fcmp_to_serial(candidate, eleNum, eleProjBFCnt);
  for(idx=0;idx<slaterCount;idx++) {
    const double diff = cabs(candidate[idx] - slaterBefore[idx]);
    if(diff > slaterMaxDiff) slaterMaxDiff = diff;
  }

  slaterGlobalStateChanged =
      memcmp(slaterBefore, SlaterElmBF,
             slaterCount*sizeof(double complex)) != 0
      || memcmp(slaterRealBefore, SlaterElmBF_real,
                slaterCount*sizeof(double)) != 0
      || memcmp(invPfBefore, InvM,
                invPfCount*sizeof(double complex)) != 0;
  for(ri=0;ri<Nsite;ri++) {
    if(memcmp(etaBefore + (size_t)ri*(size_t)Nsite, eta[ri],
              (size_t)Nsite*sizeof(double complex)) != 0) {
      etaStateChanged = 1;
    }
    if(memcmp(etaFlagBefore + (size_t)ri*(size_t)Nsite, etaFlag[ri],
              (size_t)Nsite*sizeof(int)) != 0) {
      etaFlagStateChanged = 1;
    }
  }

  for(idx=0;idx<(size_t)Nsize;idx++) {
    eleSpn[idx] = idx < (size_t)Ne ? 0 : 1;
  }
  for(idx=0;idx<(size_t)NQPFull;idx++) {
    pfCandidate[idx] = 1.0 + 0.0*I;
    pfReference[idx] = 1.0 + 0.0*I;
    pfZero[idx] = 1.0 + 0.0*I;
  }
  pfReferenceStatus = CalculateMAll_BF_fcmp(eleIdx, 0, NQPFull);
  memcpy(pfReference, PfM, (size_t)NQPFull*sizeof(double complex));
  memcpy(invPfAfterReference, InvM,
         invPfCount*sizeof(double complex));

  pfStatus = CalculatePfM_BF_from_workspace(
      candidate, eleIdx, eleSpn, 0, NQPFull, pfCandidate, &pfFailureDetail,
      pfBufM, pfIWork, pfWork, LapackLWork, pfRWork);
  for(idx=0;idx<(size_t)NQPFull;idx++) {
    const double diff = cabs(pfCandidate[idx] - pfReference[idx]);
    if(diff > pfMaxDiff) pfMaxDiff = diff;
  }

  pfZeroStatus = CalculatePfM_BF_from_workspace(
      zeroSlater, eleIdx, eleSpn, 0, NQPFull, pfZero, NULL,
      pfBufM, pfIWork, pfWork, LapackLWork, pfRWork);
  for(idx=0;idx<(size_t)NQPFull;idx++) {
    if(creal(pfZero[idx]) != 0.0 || cimag(pfZero[idx]) != 0.0) {
      pfZeroCandidateNonzero = 1;
    }
  }

  invalidDetail = -1;
  pfInvalidNullStatus = CalculatePfM_BF_from_workspace(
      NULL, eleIdx, eleSpn, 0, NQPFull, pfCandidate,
      &invalidDetail, pfBufM, pfIWork, pfWork, LapackLWork,
      pfRWork);
  if(invalidDetail != 0) pfInvalidFailureDetail = 1;
  invalidDetail = -1;
  pfInvalidRangeStatus = CalculatePfM_BF_from_workspace(
      candidate, eleIdx, eleSpn, 0, 0, pfCandidate,
      &invalidDetail, pfBufM, pfIWork, pfWork, LapackLWork,
      pfRWork);
  if(invalidDetail != 0) pfInvalidFailureDetail = 1;
  invalidDetail = -1;
  pfInvalidWorkspaceStatus = CalculatePfM_BF_from_workspace(
      candidate, eleIdx, eleSpn, 0, NQPFull, pfCandidate,
      &invalidDetail, pfBufM, pfIWork, pfWork, LapackLWork-1,
      pfRWork);
  if(invalidDetail != 0) pfInvalidFailureDetail = 1;

  pfGlobalStateChanged =
      memcmp(slaterBefore, SlaterElmBF,
             slaterCount*sizeof(double complex)) != 0
      || memcmp(slaterRealBefore, SlaterElmBF_real,
                slaterCount*sizeof(double)) != 0
      || memcmp(invPfAfterReference, InvM,
                invPfCount*sizeof(double complex)) != 0;

  /* Keep this test-only diagnostic non-invasive if a future regression occurs. */
  memcpy(SlaterElmBF, slaterBefore, slaterCount*sizeof(double complex));
  memcpy(SlaterElmBF_real, slaterRealBefore, slaterCount*sizeof(double));
  memcpy(InvM, invPfBefore, invPfCount*sizeof(double complex));
  for(ri=0;ri<Nsite;ri++) {
    memcpy(eta[ri], etaBefore + (size_t)ri*(size_t)Nsite,
           (size_t)Nsite*sizeof(double complex));
    memcpy(etaFlag[ri], etaFlagBefore + (size_t)ri*(size_t)Nsite,
           (size_t)Nsite*sizeof(int));
  }

  fp = fopen(path, "w");
  if(fp == NULL) {
    fprintf(stderr,
            "Error: failed to open BackFlow N-body component dump file: %s\n",
            path);
    goto cleanup;
  }
  fprintf(fp, "slater_max_abs_diff %.17e\n", slaterMaxDiff);
  fprintf(fp, "slater_global_state_changed %d\n",
          slaterGlobalStateChanged);
  fprintf(fp, "eta_state_changed %d\n", etaStateChanged);
  fprintf(fp, "eta_flag_state_changed %d\n", etaFlagStateChanged);
  fprintf(fp, "pf_max_abs_diff %.17e\n", pfMaxDiff);
  fprintf(fp, "pf_reference_status %d\n", pfReferenceStatus);
  fprintf(fp, "pf_status %d\n", pfStatus);
  fprintf(fp, "pf_failure_detail %d\n", pfFailureDetail);
  fprintf(fp, "pf_global_state_changed %d\n", pfGlobalStateChanged);
  fprintf(fp, "pf_zero_status %d\n", pfZeroStatus);
  fprintf(fp, "pf_zero_candidate_nonzero %d\n",
          pfZeroCandidateNonzero);
  fprintf(fp, "pf_invalid_null_status %d\n", pfInvalidNullStatus);
  fprintf(fp, "pf_invalid_range_status %d\n", pfInvalidRangeStatus);
  fprintf(fp, "pf_invalid_workspace_status %d\n",
          pfInvalidWorkspaceStatus);
  fprintf(fp, "pf_invalid_failure_detail %d\n",
          pfInvalidFailureDetail);
  fprintf(fp, "writer_parent_rank %d\n", parentRank);
  fprintf(fp, "writer_child_rank %d\n", childRank);
  status = 0;

cleanup:
  if(fp != NULL) fclose(fp);
  free(candidate);
  free(zeroSlater);
  free(slaterBefore);
  free(slaterRealBefore);
  free(invPfBefore);
  free(invPfAfterReference);
  free(etaBefore);
  free(etaFlagBefore);
  free(pfCandidate);
  free(pfReference);
  free(pfZero);
  free(pfBufM);
  free(pfWork);
  free(eleSpn);
  free(pfIWork);
  free(pfRWork);
  return status;
}

typedef struct {
  BFNBodyScratchSizes sizes;
  BFNBodyScratch scratch;
  int *intBase;
  double complex *complexBase;
  double *doubleBase;
} BFDiagScratch;

static int initBFDiagScratch(BFDiagScratch *owned) {
  if(owned == NULL) return -1;
  memset(owned, 0, sizeof(*owned));
  if(GetBFNBodyScratchSizes(4, 0, &owned->sizes) != BF_NBODY_OK) {
    return -1;
  }
  owned->intBase =
      (int *)malloc(owned->sizes.intCount*sizeof(int));
  owned->complexBase = (double complex *)malloc(
      owned->sizes.complexCount*sizeof(double complex));
  owned->doubleBase =
      (double *)malloc(owned->sizes.doubleCount*sizeof(double));
  if(owned->intBase == NULL || owned->complexBase == NULL
     || owned->doubleBase == NULL) {
    return -1;
  }
  memset(owned->intBase, 0x5a,
         owned->sizes.intCount*sizeof(int));
  memset(owned->complexBase, 0xa5,
         owned->sizes.complexCount*sizeof(double complex));
  memset(owned->doubleBase, 0x3c,
         owned->sizes.doubleCount*sizeof(double));
  return BindBFNBodyScratch(
      &owned->sizes, owned->intBase, owned->sizes.intCount,
      owned->complexBase, owned->sizes.complexCount,
      owned->doubleBase, owned->sizes.doubleCount, &owned->scratch)
      == BF_NBODY_OK ? 0 : -1;
}

static void freeBFDiagScratch(BFDiagScratch *owned) {
  if(owned == NULL) return;
  free(owned->intBase);
  free(owned->complexBase);
  free(owned->doubleBase);
  memset(owned, 0, sizeof(*owned));
}

static int BFDiagResultMatches(BFNBodyResult result,
                               BFNBodyStatus status,
                               BFNBodyStage stage, int detail,
                               int reducedOrder,
                               double complex value, double tolerance) {
  return result.status == status && result.stage == stage
      && result.detail == detail && result.reducedOrder == reducedOrder
      && isfinite(creal(result.value)) && isfinite(cimag(result.value))
      && cabs(result.value-value) <= tolerance;
}

static int BFDiagDirectValue(
    int n, const int *rsi, const int *rsj, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const int *eleProjBFCnt,
    BFNBodyScratch *scratch, NBodyReduction *reductionOut,
    double complex *valueOut) {
  NBodyReduction reduction;
  double complex value;
  const size_t slaterCount =
      (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;

  if(scratch == NULL || reductionOut == NULL || valueOut == NULL
     || ReduceNBodyTerm(n, rsi, rsj, eleNum, Nsite2,
                        scratch->rsi, scratch->rsj,
                        scratch->maxOrder, &reduction) != 0) {
    return -1;
  }
  *reductionOut = reduction;
  if(reduction.kind == NBODY_REDUCED_ZERO) {
    *valueOut = 0.0+0.0*I;
    return 0;
  }
  if(reduction.kind == NBODY_REDUCED_SCALAR) {
    *valueOut = (double)reduction.sign+0.0*I;
    return 0;
  }
  if(reduction.kind != NBODY_REDUCED_HOPS
     || reduction.order < 1 || reduction.order > 2) {
    return -1;
  }

  memcpy(scratch->eleIdx, eleIdx, (size_t)Nsize*sizeof(int));
  memcpy(scratch->eleCfg, eleCfg, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->eleNum, eleNum, (size_t)Nsite2*sizeof(int));
  if(NProj > 0) {
    memcpy(scratch->projCnt, eleProjCnt, (size_t)NProj*sizeof(int));
  }
  memcpy(scratch->projBFCnt, eleProjBFCnt,
         (size_t)16*(size_t)Nsite*(size_t)Nrange*sizeof(int));
  memcpy(scratch->slater, SlaterElmBF,
         slaterCount*sizeof(double complex));

  if(reduction.order == 1) {
    value = GreenFunc1BF(
        scratch->rsi[0]%Nsite, scratch->rsj[0]%Nsite,
        scratch->rsj[0]/Nsite, ip, scratch->slater,
        scratch->eleIdx, scratch->eleCfg, scratch->eleNum,
        eleProjCnt, scratch->projCnt, eleProjBFCnt,
        scratch->projBFCnt, scratch->greenBuffer);
  } else {
    value = GreenFunc2BF(
        scratch->rsi[0]%Nsite, scratch->rsj[0]%Nsite,
        scratch->rsi[1]%Nsite, scratch->rsj[1]%Nsite,
        scratch->rsj[0]/Nsite, scratch->rsj[1]/Nsite, ip,
        scratch->slater, scratch->eleIdx, scratch->eleCfg,
        scratch->eleNum, eleProjCnt, scratch->projCnt,
        eleProjBFCnt, scratch->projBFCnt, scratch->greenBuffer);
  }
  if(memcmp(scratch->eleIdx, eleIdx, (size_t)Nsize*sizeof(int)) != 0
     || memcmp(scratch->eleCfg, eleCfg,
               (size_t)Nsite2*sizeof(int)) != 0
     || memcmp(scratch->eleNum, eleNum,
               (size_t)Nsite2*sizeof(int)) != 0
     || memcmp(scratch->slater, SlaterElmBF,
               slaterCount*sizeof(double complex)) != 0) {
    return -1;
  }
  *valueOut = (double)reduction.sign*value;
  return isfinite(creal(*valueOut)) && isfinite(cimag(*valueOut)) ? 0 : -1;
}

static void BFDiagRestoreGlobals(
    const double complex *slater, const double *slaterReal,
    const double complex *invM, const double complex *pfM,
    const double complex *etaValues, const int *etaFlags) {
  const size_t slaterCount =
      (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
  const size_t invCount =
      (size_t)NQPFull*(size_t)Nsize*(size_t)Nsize;
  int ri;

  memcpy(SlaterElmBF, slater, slaterCount*sizeof(double complex));
  memcpy(SlaterElmBF_real, slaterReal, slaterCount*sizeof(double));
  memcpy(InvM, invM, invCount*sizeof(double complex));
  memcpy(PfM, pfM, (size_t)NQPFull*sizeof(double complex));
  for(ri=0;ri<Nsite;ri++) {
    memcpy(eta[ri], etaValues+(size_t)ri*(size_t)Nsite,
           (size_t)Nsite*sizeof(double complex));
    memcpy(etaFlag[ri], etaFlags+(size_t)ri*(size_t)Nsite,
           (size_t)Nsite*sizeof(int));
  }
}

static int BFDiagGlobalsDiffer(
    const double complex *slater, const double *slaterReal,
    const double complex *invM, const double complex *pfM,
    const double complex *etaValues, const int *etaFlags) {
  const size_t slaterCount =
      (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
  const size_t invCount =
      (size_t)NQPFull*(size_t)Nsize*(size_t)Nsize;
  int ri;

  if(memcmp(SlaterElmBF, slater,
            slaterCount*sizeof(double complex)) != 0
     || memcmp(SlaterElmBF_real, slaterReal,
               slaterCount*sizeof(double)) != 0
     || memcmp(InvM, invM, invCount*sizeof(double complex)) != 0
     || memcmp(PfM, pfM,
               (size_t)NQPFull*sizeof(double complex)) != 0) {
    return 1;
  }
  for(ri=0;ri<Nsite;ri++) {
    if(memcmp(eta[ri], etaValues+(size_t)ri*(size_t)Nsite,
              (size_t)Nsite*sizeof(double complex)) != 0
       || memcmp(etaFlag[ri], etaFlags+(size_t)ri*(size_t)Nsite,
                 (size_t)Nsite*sizeof(int)) != 0) {
      return 1;
    }
  }
  return 0;
}

static int BFDiagFullRebuildValue(
    int n, const int *rsi, const int *rsj, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, BFNBodyScratch *scratch,
    const double complex *slaterBefore, const double *slaterRealBefore,
    const double complex *invBefore, const double complex *pfBefore,
    const double complex *etaBefore, const int *etaFlagBefore,
    double complex *valueOut) {
  double complex ipNew;
  double projRatio;
  int k;
  int info = -1;

  memcpy(scratch->eleIdx, eleIdx, (size_t)Nsize*sizeof(int));
  memcpy(scratch->eleCfg, eleCfg, (size_t)Nsite2*sizeof(int));
  memcpy(scratch->eleNum, eleNum, (size_t)Nsite2*sizeof(int));
  for(k=0;k<n;k++) {
    const int source = rsj[k];
    const int spin = source/Nsite;
    const int localLabel = eleCfg[source];
    const int particle = localLabel+spin*Ne;
    if(rsi[k]/Nsite != spin || localLabel < 0 || localLabel >= Ne
       || particle < 0 || particle >= Nsize
       || eleNum[source] != 1 || eleNum[rsi[k]] != 0
       || eleIdx[particle] != source%Nsite) {
      goto restore;
    }
    scratch->moved[k] = particle;
  }
  for(k=0;k<n;k++) {
    scratch->eleCfg[rsj[k]] = -1;
    scratch->eleNum[rsj[k]] = 0;
  }
  for(k=0;k<n;k++) {
    const int particle = scratch->moved[k];
    scratch->eleCfg[rsi[k]] = particle%Ne;
    scratch->eleNum[rsi[k]] = 1;
    scratch->eleIdx[particle] = rsi[k]%Nsite;
  }
  MakeProjCnt(scratch->projCnt, scratch->eleNum);
  MakeProjBFCnt(scratch->projBFCnt, scratch->eleNum);
  if(!IsSectorStateAllowed(scratch->eleNum)) {
    *valueOut = 0.0+0.0*I;
    info = 0;
    goto restore;
  }
  MakeSlaterElmBF_fcmp(scratch->eleNum, scratch->projBFCnt);
  if(CalculateMAll_BF_fcmp(scratch->eleIdx, 0, NQPFull) != 0) {
    goto restore;
  }
  ipNew = CalculateIP_fcmp(PfM, 0, NQPFull, MPI_COMM_SELF);
  projRatio = ProjRatio(scratch->projCnt, eleProjCnt);
  *valueOut = conj(projRatio*ipNew/ip);
  if(isfinite(creal(*valueOut)) && isfinite(cimag(*valueOut))) info = 0;

restore:
  BFDiagRestoreGlobals(
      slaterBefore, slaterRealBefore, invBefore, pfBefore,
      etaBefore, etaFlagBefore);
  return info;
}

static int dumpBFNBodyDispatchCheck(const char *path, const int *eleIdx,
                                    const int *eleCfg, const int *eleNum,
                                    const int *eleProjCnt,
                                    const int *eleProjBFCnt) {
  const double tolerance = 1.0e-11;
  const size_t slaterCount =
      (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
  const size_t invCount =
      (size_t)NQPFull*(size_t)Nsize*(size_t)Nsize;
  const size_t etaCount = (size_t)Nsite*(size_t)Nsite;
  const size_t projBFCount =
      (size_t)16*(size_t)Nsite*(size_t)Nrange;
  BFDiagScratch evaluator;
  BFDiagScratch reference;
  double complex *slaterBefore = NULL;
  double *slaterRealBefore = NULL;
  double complex *invBefore = NULL;
  double complex *pfBefore = NULL;
  double complex *etaBefore = NULL;
  int *etaFlagBefore = NULL;
  int *idxBefore = NULL;
  int *cfgBefore = NULL;
  int *numBefore = NULL;
  int *projBefore = NULL;
  int *projBFBefore = NULL;
  int *badCfg = NULL;
  int *occupied = NULL;
  int *empty = NULL;
  int source[4], target[4];
  int pairStart[2], pairCount[2];
  int scalarCount = 0;
  int zeroCount = 0;
  int dispatch1Count = 0;
  int dispatch2Count = 0;
  int fullRebuildCount = 0;
  int contractFailures = 0;
  int setupFailures = 0;
  int mixedOrderFailures = 0;
  int callerStateChanged = 0;
  int globalStateChanged = 0;
  double maxDirectDiff = 0.0;
  double maxFullRebuildDiff = 0.0;
  double maxN4Imag = 0.0;
  double complex ip = 0.0+0.0*I;
  FILE *fp = NULL;
  int spin;
  int k;
  int status = -1;

  memset(&evaluator, 0, sizeof(evaluator));
  memset(&reference, 0, sizeof(reference));
  if(path == NULL || initBFDiagScratch(&evaluator) != 0
     || initBFDiagScratch(&reference) != 0) {
    goto cleanup;
  }
  slaterBefore =
      (double complex *)malloc(slaterCount*sizeof(double complex));
  slaterRealBefore = (double *)malloc(slaterCount*sizeof(double));
  invBefore = (double complex *)malloc(invCount*sizeof(double complex));
  pfBefore =
      (double complex *)malloc((size_t)NQPFull*sizeof(double complex));
  etaBefore = (double complex *)malloc(etaCount*sizeof(double complex));
  etaFlagBefore = (int *)malloc(etaCount*sizeof(int));
  idxBefore = (int *)malloc((size_t)Nsize*sizeof(int));
  cfgBefore = (int *)malloc((size_t)Nsite2*sizeof(int));
  numBefore = (int *)malloc((size_t)Nsite2*sizeof(int));
  if(NProj > 0) {
    projBefore = (int *)malloc((size_t)NProj*sizeof(int));
  }
  projBFBefore = (int *)malloc(projBFCount*sizeof(int));
  badCfg = (int *)malloc((size_t)Nsite2*sizeof(int));
  occupied = (int *)malloc((size_t)Nsite2*sizeof(int));
  empty = (int *)malloc((size_t)Nsite2*sizeof(int));
  if(slaterBefore == NULL || slaterRealBefore == NULL
     || invBefore == NULL || pfBefore == NULL || etaBefore == NULL
     || etaFlagBefore == NULL || idxBefore == NULL || cfgBefore == NULL
     || numBefore == NULL || (NProj > 0 && projBefore == NULL)
     || projBFBefore == NULL || badCfg == NULL
     || occupied == NULL || empty == NULL) {
    goto cleanup;
  }

  if(CalculateMAll_BF_fcmp(eleIdx, 0, NQPFull) != 0) goto cleanup;
  ip = CalculateIP_fcmp(PfM, 0, NQPFull, MPI_COMM_SELF);
  if(!isfinite(creal(ip)) || !isfinite(cimag(ip)) || cabs(ip) == 0.0) {
    goto cleanup;
  }
  memcpy(slaterBefore, SlaterElmBF,
         slaterCount*sizeof(double complex));
  memcpy(slaterRealBefore, SlaterElmBF_real,
         slaterCount*sizeof(double));
  memcpy(invBefore, InvM, invCount*sizeof(double complex));
  memcpy(pfBefore, PfM, (size_t)NQPFull*sizeof(double complex));
  memcpy(idxBefore, eleIdx, (size_t)Nsize*sizeof(int));
  memcpy(cfgBefore, eleCfg, (size_t)Nsite2*sizeof(int));
  memcpy(numBefore, eleNum, (size_t)Nsite2*sizeof(int));
  if(NProj > 0) {
    memcpy(projBefore, eleProjCnt, (size_t)NProj*sizeof(int));
  }
  memcpy(projBFBefore, eleProjBFCnt, projBFCount*sizeof(int));
  for(k=0;k<Nsite;k++) {
    memcpy(etaBefore+(size_t)k*(size_t)Nsite, eta[k],
           (size_t)Nsite*sizeof(double complex));
    memcpy(etaFlagBefore+(size_t)k*(size_t)Nsite, etaFlag[k],
           (size_t)Nsite*sizeof(int));
  }

  k = 0;
  for(spin=0;spin<2;spin++) {
    int occupiedCount = 0;
    int emptyCount = 0;
    int orbital;
    pairStart[spin] = k;
    for(orbital=spin*Nsite;orbital<(spin+1)*Nsite;orbital++) {
      if(eleNum[orbital] == 1) occupied[occupiedCount++] = orbital;
      else empty[emptyCount++] = orbital;
    }
    pairCount[spin] =
        occupiedCount < emptyCount ? occupiedCount : emptyCount;
    for(orbital=0;orbital<pairCount[spin] && k<4;orbital++) {
      source[k] = occupied[orbital];
      target[k] = empty[orbital];
      k++;
    }
  }
  if(k != 4 || pairCount[0] != 2 || pairCount[1] != 2) {
    setupFailures++;
    goto write_dump;
  }

  {
    const int rsi[1] = {source[0]};
    const int rsj[1] = {source[0]};
    BFNBodyResult result = GreenFuncNBF(
        1, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
        eleProjCnt, eleProjBFCnt, &evaluator.scratch);
    scalarCount++;
    if(!BFDiagResultMatches(
           result, BF_NBODY_OK, BF_NBODY_STAGE_REDUCE,
           BF_NBODY_DETAIL_NONE, 0, 1.0+0.0*I, 0.0)) {
      contractFailures++;
    }
  }
  {
    const int rsi[1] = {target[0]};
    const int rsj[1] = {target[0]};
    BFNBodyResult result = GreenFuncNBF(
        1, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
        eleProjCnt, eleProjBFCnt, &evaluator.scratch);
    zeroCount++;
    if(!BFDiagResultMatches(
           result, BF_NBODY_PHYSICAL_ZERO, BF_NBODY_STAGE_REDUCE,
           BF_NBODY_DETAIL_NONE, 0, 0.0+0.0*I, 0.0)) {
      contractFailures++;
    }
  }
  {
    const int rsi[1] = {target[0]};
    const int rsj[1] = {source[0]};
    NBodyReduction reduction;
    double complex direct;
    BFNBodyResult result;
    if(BFDiagDirectValue(
           1, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
           eleProjCnt, eleProjBFCnt, &reference.scratch,
           &reduction, &direct) != 0) {
      setupFailures++;
    } else {
      const BFNBodyStatus expectedStatus =
          cabs(direct) == 0.0 ? BF_NBODY_PHYSICAL_ZERO : BF_NBODY_OK;
      const BFNBodyStage expectedStage =
          cabs(direct) == 0.0 ? BF_NBODY_STAGE_RATIO : BF_NBODY_STAGE_NONE;
      result = GreenFuncNBF(
          1, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
          eleProjCnt, eleProjBFCnt, &evaluator.scratch);
      dispatch1Count++;
      if(!BFDiagResultMatches(
             result, expectedStatus, expectedStage, BF_NBODY_DETAIL_NONE,
             1, direct, tolerance)) {
        contractFailures++;
      }
      if(cabs(result.value-direct) > maxDirectDiff) {
        maxDirectDiff = cabs(result.value-direct);
      }
    }
  }
  {
    const int rsi[2] = {target[0], target[1]};
    const int rsj[2] = {source[0], source[1]};
    NBodyReduction reduction;
    double complex direct;
    BFNBodyResult result;
    if(BFDiagDirectValue(
           2, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
           eleProjCnt, eleProjBFCnt, &reference.scratch,
           &reduction, &direct) != 0) {
      setupFailures++;
    } else {
      const BFNBodyStatus expectedStatus =
          cabs(direct) == 0.0 ? BF_NBODY_PHYSICAL_ZERO : BF_NBODY_OK;
      const BFNBodyStage expectedStage =
          cabs(direct) == 0.0 ? BF_NBODY_STAGE_RATIO : BF_NBODY_STAGE_NONE;
      result = GreenFuncNBF(
          2, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
          eleProjCnt, eleProjBFCnt, &evaluator.scratch);
      dispatch2Count++;
      if(!BFDiagResultMatches(
             result, expectedStatus, expectedStage, BF_NBODY_DETAIL_NONE,
             2, direct, tolerance)) {
        contractFailures++;
      }
      if(cabs(result.value-direct) > maxDirectDiff) {
        maxDirectDiff = cabs(result.value-direct);
      }
    }
  }
  {
    const int rsi[3] = {target[0], target[2], target[1]};
    const int rsj[3] = {target[2], source[0], source[1]};
    NBodyReduction reduction;
    double complex direct;
    BFNBodyResult result;
    if(BFDiagDirectValue(
           3, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
           eleProjCnt, eleProjBFCnt, &reference.scratch,
           &reduction, &direct) != 0 || reduction.order != 2
       || reduction.sign != 1) {
      setupFailures++;
    } else {
      const BFNBodyStatus expectedStatus =
          cabs(direct) == 0.0 ? BF_NBODY_PHYSICAL_ZERO : BF_NBODY_OK;
      const BFNBodyStage expectedStage =
          cabs(direct) == 0.0 ? BF_NBODY_STAGE_RATIO : BF_NBODY_STAGE_NONE;
      result = GreenFuncNBF(
          3, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
          eleProjCnt, eleProjBFCnt, &evaluator.scratch);
      dispatch2Count++;
      if(!BFDiagResultMatches(
             result, expectedStatus, expectedStage, BF_NBODY_DETAIL_NONE,
             2, direct, tolerance)) {
        contractFailures++;
      }
      if(cabs(result.value-direct) > maxDirectDiff) {
        maxDirectDiff = cabs(result.value-direct);
      }
    }
  }
  {
    const int rsi[3] = {source[2], target[0], source[1]};
    const int rsj[3] = {source[0], source[2], source[1]};
    NBodyReduction reduction;
    double complex direct;
    BFNBodyResult result;
    if(BFDiagDirectValue(
           3, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
           eleProjCnt, eleProjBFCnt, &reference.scratch,
           &reduction, &direct) != 0 || reduction.order != 1
       || reduction.sign != -1) {
      setupFailures++;
    } else {
      const BFNBodyStatus expectedStatus =
          cabs(direct) == 0.0 ? BF_NBODY_PHYSICAL_ZERO : BF_NBODY_OK;
      const BFNBodyStage expectedStage =
          cabs(direct) == 0.0 ? BF_NBODY_STAGE_RATIO : BF_NBODY_STAGE_NONE;
      result = GreenFuncNBF(
          3, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
          eleProjCnt, eleProjBFCnt, &evaluator.scratch);
      dispatch1Count++;
      if(!BFDiagResultMatches(
             result, expectedStatus, expectedStage, BF_NBODY_DETAIL_NONE,
             1, direct, tolerance)) {
        contractFailures++;
      }
      if(cabs(result.value-direct) > maxDirectDiff) {
        maxDirectDiff = cabs(result.value-direct);
      }
    }
  }
  {
    const int rsi[3] = {target[0], target[1], target[2]};
    const int rsj[3] = {source[0], source[1], source[2]};
    double complex full;
    BFNBodyResult result;
    globalStateChanged |= BFDiagGlobalsDiffer(
        slaterBefore, slaterRealBefore, invBefore, pfBefore,
        etaBefore, etaFlagBefore);
    if(BFDiagFullRebuildValue(
           3, rsi, rsj, ip, eleIdx, eleCfg, eleNum, eleProjCnt,
           &reference.scratch, slaterBefore, slaterRealBefore,
           invBefore, pfBefore, etaBefore, etaFlagBefore, &full) != 0) {
      setupFailures++;
    } else {
      result = GreenFuncNBF(
          3, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
          eleProjCnt, eleProjBFCnt, &evaluator.scratch);
      fullRebuildCount++;
      if(!BFDiagResultMatches(
             result, cabs(full) == 0.0 ? BF_NBODY_PHYSICAL_ZERO
                                       : BF_NBODY_OK,
             cabs(full) == 0.0 ? BF_NBODY_STAGE_RATIO
                               : BF_NBODY_STAGE_NONE,
             BF_NBODY_DETAIL_NONE, 3, full, tolerance)) {
        contractFailures++;
      }
      if(cabs(result.value-full) > maxFullRebuildDiff) {
        maxFullRebuildDiff = cabs(result.value-full);
      }
    }
  }
  {
    int permutation;
    for(permutation=0;permutation<4;permutation++) {
      int rsi[4] = {target[0], target[1], target[2], target[3]};
      const int rsj[4] = {source[0], source[1], source[2], source[3]};
      double complex full;
      BFNBodyResult result;
      double diff;
      if(permutation & 1) {
        const int tmp = rsi[pairStart[0]];
        rsi[pairStart[0]] = rsi[pairStart[0]+1];
        rsi[pairStart[0]+1] = tmp;
      }
      if(permutation & 2) {
        const int tmp = rsi[pairStart[1]];
        rsi[pairStart[1]] = rsi[pairStart[1]+1];
        rsi[pairStart[1]+1] = tmp;
      }
      globalStateChanged |= BFDiagGlobalsDiffer(
          slaterBefore, slaterRealBefore, invBefore, pfBefore,
          etaBefore, etaFlagBefore);
      if(BFDiagFullRebuildValue(
             4, rsi, rsj, ip, eleIdx, eleCfg, eleNum, eleProjCnt,
             &reference.scratch, slaterBefore, slaterRealBefore,
             invBefore, pfBefore, etaBefore, etaFlagBefore, &full) != 0) {
        setupFailures++;
        continue;
      }
      result = GreenFuncNBF(
          4, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
          eleProjCnt, eleProjBFCnt, &evaluator.scratch);
      fullRebuildCount++;
      if(!BFDiagResultMatches(
             result, cabs(full) == 0.0 ? BF_NBODY_PHYSICAL_ZERO
                                       : BF_NBODY_OK,
             cabs(full) == 0.0 ? BF_NBODY_STAGE_RATIO
                               : BF_NBODY_STAGE_NONE,
             BF_NBODY_DETAIL_NONE, 4, full, tolerance)) {
        contractFailures++;
      }
      diff = cabs(result.value-full);
      if(diff > maxFullRebuildDiff) maxFullRebuildDiff = diff;
      if(fabs(cimag(result.value)) > maxN4Imag) {
        maxN4Imag = fabs(cimag(result.value));
      }
    }
  }
  {
    const int rsi[1] = {target[0]};
    const int rsj[1] = {source[0]};
    NBodyReduction reduction;
    double complex direct;
    BFNBodyResult result;
    if(BFDiagDirectValue(
           1, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
           eleProjCnt, eleProjBFCnt, &reference.scratch,
           &reduction, &direct) != 0) {
      mixedOrderFailures++;
    } else {
      result = GreenFuncNBF(
          1, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
          eleProjCnt, eleProjBFCnt, &evaluator.scratch);
      dispatch1Count++;
      if(cabs(result.value-direct) > tolerance
         || result.reducedOrder != 1
         || (cabs(direct) == 0.0
             ? result.status != BF_NBODY_PHYSICAL_ZERO
             : result.status != BF_NBODY_OK)) {
        mixedOrderFailures++;
      }
    }
  }
  {
    const int rsi[1] = {target[0]};
    const int rsj[1] = {source[0]};
    const int crossTarget = target[pairStart[1]];
    BFNBodyResult result;
    result = GreenFuncNBF(
        0, rsi, rsj, ip, eleIdx, eleCfg, eleNum,
        eleProjCnt, eleProjBFCnt, &evaluator.scratch);
    if(!BFDiagResultMatches(
           result, BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_REDUCE,
           BF_NBODY_DETAIL_REDUCER, 0, 0.0+0.0*I, 0.0)) {
      contractFailures++;
    }
    result = GreenFuncNBF(
        1, rsi, rsj, 0.0+0.0*I, eleIdx, eleCfg, eleNum,
        eleProjCnt, eleProjBFCnt, &evaluator.scratch);
    if(!BFDiagResultMatches(
           result, BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_RATIO,
           BF_NBODY_DETAIL_NONE, 1, 0.0+0.0*I, 0.0)) {
      contractFailures++;
    }
    result = GreenFuncNBF(
        1, rsi, rsj, NAN+0.0*I, eleIdx, eleCfg, eleNum,
        eleProjCnt, eleProjBFCnt, &evaluator.scratch);
    if(!BFDiagResultMatches(
           result, BF_NBODY_NONFINITE, BF_NBODY_STAGE_RATIO,
           BF_NBODY_DETAIL_NONE, 1, 0.0+0.0*I, 0.0)) {
      contractFailures++;
    }
    result = GreenFuncNBF(
        1, &crossTarget, rsj, ip, eleIdx, eleCfg, eleNum,
        eleProjCnt, eleProjBFCnt, &evaluator.scratch);
    if(!BFDiagResultMatches(
           result, BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_DISPATCH,
           BF_NBODY_DETAIL_SPIN_CHANGE, 1, 0.0+0.0*I, 0.0)) {
      contractFailures++;
    }
    memcpy(badCfg, eleCfg, (size_t)Nsite2*sizeof(int));
    badCfg[source[0]] = -1;
    result = GreenFuncNBF(
        1, rsi, rsj, ip, eleIdx, badCfg, eleNum,
        eleProjCnt, eleProjBFCnt, &evaluator.scratch);
    if(!BFDiagResultMatches(
           result, BF_NBODY_INVALID_ARGUMENT, BF_NBODY_STAGE_CANDIDATE,
           BF_NBODY_DETAIL_BAD_ELECTRON_LABEL, 1,
           0.0+0.0*I, 0.0)) {
      contractFailures++;
    }
  }

  callerStateChanged =
      memcmp(eleIdx, idxBefore, (size_t)Nsize*sizeof(int)) != 0
      || memcmp(eleCfg, cfgBefore, (size_t)Nsite2*sizeof(int)) != 0
      || memcmp(eleNum, numBefore, (size_t)Nsite2*sizeof(int)) != 0
      || (NProj > 0
          && memcmp(eleProjCnt, projBefore,
                    (size_t)NProj*sizeof(int)) != 0)
      || memcmp(eleProjBFCnt, projBFBefore,
                projBFCount*sizeof(int)) != 0;
  globalStateChanged |= BFDiagGlobalsDiffer(
      slaterBefore, slaterRealBefore, invBefore, pfBefore,
      etaBefore, etaFlagBefore);

write_dump:
  fp = fopen(path, "w");
  if(fp == NULL) goto cleanup;
  fprintf(fp, "implemented 1\n");
  fprintf(fp, "scalar %d\n", scalarCount);
  fprintf(fp, "zero %d\n", zeroCount);
  fprintf(fp, "dispatch1 %d\n", dispatch1Count);
  fprintf(fp, "dispatch2 %d\n", dispatch2Count);
  fprintf(fp, "full_rebuild %d\n", fullRebuildCount);
  fprintf(fp, "contract_failures %d\n", contractFailures);
  fprintf(fp, "setup_failures %d\n", setupFailures);
  fprintf(fp, "mixed_order_failures %d\n", mixedOrderFailures);
  fprintf(fp, "normal_projection_count %d\n", NProj);
  fprintf(fp, "caller_state_changed %d\n", callerStateChanged);
  fprintf(fp, "global_state_changed %d\n", globalStateChanged);
  fprintf(fp, "max_direct_diff %.17e\n", maxDirectDiff);
  fprintf(fp, "max_full_rebuild_diff %.17e\n", maxFullRebuildDiff);
  fprintf(fp, "max_n4_imag %.17e\n", maxN4Imag);
  status = 0;

cleanup:
  if(fp != NULL) fclose(fp);
  freeBFDiagScratch(&evaluator);
  freeBFDiagScratch(&reference);
  free(slaterBefore);
  free(slaterRealBefore);
  free(invBefore);
  free(pfBefore);
  free(etaBefore);
  free(etaFlagBefore);
  free(idxBefore);
  free(cfgBefore);
  free(numBefore);
  free(projBefore);
  free(projBFBefore);
  free(badCfg);
  free(occupied);
  free(empty);
  return status;
}

void VMC_BF_MainCal(MPI_Comm comm_parent, MPI_Comm comm) {
  int *eleIdx, *eleCfg, *eleNum, *eleProjCnt, *eleProjBFCnt;
  double complex e, ip; //db is double?
  double w, db;
  double sqrtw;
  double complex we;
  int int_i, sampleSize, tmp_i;
  int bfIdentityDumped = 0;
  const char *bfIdentityDumpPath = getenv("MVMC_BF_IDENTITY_DUMP");
  const char *bfDiffDumpPath = getenv("MVMC_BF_DIFF_DUMP");
  const char *bfFDDumpPath = getenv("MVMC_BF_FD_DUMP");
  const char *bfGreen2DumpPath = getenv("MVMC_BF_GREEN2_DUMP");
  const char *bfNBodyComponentDumpPath =
      getenv("MVMC_BF_NBODY_COMPONENT_DUMP");
  int bfNBodyComponentDumped = 0;
  const char *bfNBodyDispatchDumpPath =
      getenv("MVMC_BF_NBODY_DISPATCH_DUMP");
  int bfNBodyDispatchDumped = 0;
  const int qpStart = 0;
  const int qpEnd = NQPFull;
  int sample, sampleStart, sampleEnd;
  int i, info;
  double complex *InvM_Moto, *PfM_Moto;
  LSLanczosBFScratch lanczosScratch;
  int lanczosScratchReady = 0;
  int lanczosWarnings = 0;
  long lanczosInjectSample = -1;
  long lanczosInjectParentRank = -1;
  long long lanczosCheckedSamples = 0;
  long long lanczosRejectedSamples = 0;
  FILE *lanczosOracleDump = NULL;
//  double *InvM_real_Moto, *PfM_real_Moto;

  /* optimazation for Kei */
  const int nProj = NProj;
  double complex *srOptO = SROptO;
  double *srOptO_real = SROptO_real;

//  double rtmp;

  int rank, size, parentRank, parentSize;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm_parent, &parentSize);
  MPI_Comm_rank(comm_parent, &parentRank);

  SplitLoop(&sampleStart, &sampleEnd, NVMCSample, rank, size);
  //SplitLoop(&sampleStart,&sampleEnd,NExactSample,rank,size);

  /* initialization */
  StartTimer(24);
  clearPhysQuantity();
  StopTimer(24);

  memset(&lanczosScratch, 0, sizeof(lanczosScratch));
  if(NVMCCalMode == 1 && NLanczosMode == 1) {
    const char *dumpValue = getenv("MVMC_LANCZOS_ORACLE_DUMP");
    lanczosInjectSample = lsbfNonfiniteInjectionSample();
    lanczosInjectParentRank = LSLanczosTestNonfiniteParentRank();
    if(lanczosInjectSample == -2) {
      if(rank == 0) {
        fprintf(stderr,
                "Error: MVMC_LANCZOS_TEST_NONFINITE_SAMPLE must be a non-negative integer.\n");
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
    if(lanczosInjectParentRank == -2 ||
       lanczosInjectParentRank >= parentSize) {
      if(parentRank == 0) {
        fprintf(stderr,
                "Error: MVMC_LANCZOS_TEST_NONFINITE_PARENT_RANK must identify a parent communicator rank.\n");
      }
      MPI_Abort(comm_parent, EXIT_FAILURE);
    }
    if(lanczosInjectParentRank >= 0 &&
       lanczosInjectParentRank != parentRank) {
      lanczosInjectSample = -1;
    }
    if(LSLanczosBFScratchInit(&lanczosScratch, AllComplexFlag == 0) != 0) {
      if(rank == 0) {
        fprintf(stderr, "Error: failed to allocate BackFlow Lanczos scratch.\n");
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
    lanczosScratchReady = 1;
    if(dumpValue != NULL && dumpValue[0] != '\0' && strcmp(dumpValue, "0") != 0) {
      char dumpPath[1024];
      int pathLength;
      const char *basePath = strcmp(dumpValue, "1") == 0
                           ? "lanczos_oracle_bf.dat" : dumpValue;
      if(parentSize > 1) {
        pathLength = snprintf(dumpPath, sizeof(dumpPath), "%s.rank%04d",
                              basePath, parentRank);
      } else {
        pathLength = snprintf(dumpPath, sizeof(dumpPath), "%s", basePath);
      }
      if(pathLength < 0 || (size_t)pathLength >= sizeof(dumpPath)) {
        fprintf(stderr,
                "Error: BackFlow Lanczos oracle dump path is too long on rank %d.\n",
                parentRank);
        MPI_Abort(comm_parent, EXIT_FAILURE);
      }
      lanczosOracleDump = fopen(dumpPath, "w");
      if(lanczosOracleDump == NULL) {
        fprintf(stderr,
                "Error: failed to open BackFlow Lanczos oracle dump '%s' on rank %d.\n",
                dumpPath, parentRank);
        MPI_Abort(comm_parent, EXIT_FAILURE);
      }
    }
  }

  if(NVMCCalMode==0 && NStoreO!=0 && NSRCG==0) {
    clearStoredOSampleRange(sampleStart, sampleEnd);
  }

  InvM_Moto = InvM;
  PfM_Moto = PfM;

//  InvM_real_Moto = InvM_real;
//  PfM_real_Moto = PfM_real;

  for (sample = sampleStart; sample < sampleEnd; sample++) {
    eleIdx = EleIdx + sample * Nsize;
    eleCfg = EleCfg + sample * Nsite2;
    eleNum = EleNum + sample * Nsite2;
    eleProjCnt = EleProjCnt + sample * NProj;
    eleProjBFCnt = EleProjBFCnt + sample * 16 * Nsite * Nrange;

    StartTimer(45);
    MakeSlaterElmBF_fcmp(eleNum, eleProjBFCnt);
    StopTimer(45);

    if(parentRank == 0 && !bfNBodyComponentDumped
       && bfNBodyComponentDumpPath != NULL
       && bfNBodyComponentDumpPath[0] != '\0') {
      if(dumpBFNBodyComponentCheck(bfNBodyComponentDumpPath, eleIdx,
                                   eleNum, eleProjBFCnt,
                                   parentRank, rank) != 0) {
        MPI_Abort(comm_parent, EXIT_FAILURE);
      }
      bfNBodyComponentDumped = 1;
    }
    if(parentRank == 0 && !bfNBodyDispatchDumped
       && bfNBodyDispatchDumpPath != NULL
       && bfNBodyDispatchDumpPath[0] != '\0') {
      if(dumpBFNBodyDispatchCheck(bfNBodyDispatchDumpPath, eleIdx, eleCfg,
                                  eleNum, eleProjCnt, eleProjBFCnt) != 0) {
        MPI_Abort(comm_parent, EXIT_FAILURE);
      }
      bfNBodyDispatchDumped = 1;
    }
    if (rank == 0 && !bfIdentityDumped &&
        bfIdentityDumpPath != NULL && bfIdentityDumpPath[0] != '\0') {
      dumpBFIdentityCheck(bfIdentityDumpPath, eleIdx, qpStart, qpEnd);
      bfIdentityDumped = 1;
    }
    if (rank == 0 && bfDiffDumpPath != NULL && bfDiffDumpPath[0] != '\0') {
      dumpBFSlaterDiffCheck(bfDiffDumpPath, eleIdx, eleNum, eleCfg, eleProjCnt,
                            eleProjBFCnt, qpStart, qpEnd,
                            sample, sample != sampleStart);
    }
    if (rank == 0 && bfFDDumpPath != NULL && bfFDDumpPath[0] != '\0') {
      dumpBFProjBFFiniteDiffCheck(bfFDDumpPath, eleIdx, eleNum, eleProjCnt,
                                  eleProjBFCnt, qpStart, qpEnd,
                                  sample, sample != sampleStart);
    }
    if (rank == 0 && bfGreen2DumpPath != NULL && bfGreen2DumpPath[0] != '\0') {
      dumpBFGreen2BruteForceCheck(bfGreen2DumpPath, eleIdx, eleCfg, eleNum,
                                  eleProjCnt, eleProjBFCnt, qpStart, qpEnd,
                                  sample, sample != sampleStart);
    }
    StartTimer(40);
    if (AllComplexFlag == 0) {
#pragma omp parallel for default(shared) private(tmp_i)
      for (tmp_i = 0; tmp_i < NQPFull * (2 * Nsite) * (2 * Nsite); tmp_i++)
        SlaterElmBF_real[tmp_i] = creal(SlaterElmBF[tmp_i]);

      info = CalculateMAll_BF_real(eleIdx, qpStart, qpEnd);  // InvM_real,PfM_real will change
#pragma omp parallel for default(shared) private(tmp_i)
      for (tmp_i = 0; tmp_i < NQPFull * (Nsize * Nsize + 1); tmp_i++)
        InvM[tmp_i] = InvM_real[tmp_i]; // InvM will be used in  SlaterElmDiff_fcmp
    } else {//complex
      info = CalculateMAll_BF_fcmp(eleIdx,qpStart,qpEnd);
    }
    StopTimer(40);

    if (info != 0) {
      fprintf(stderr, "waring: VMCMainCal rank:%d sample:%d info:%d (CalculateMAll)\n", rank, sample, info);
      continue;
    }

    if (AllComplexFlag == 0) {
      ip = CalculateIP_real(PfM_real, qpStart, qpEnd, MPI_COMM_SELF);
    } else {
      ip = CalculateIP_fcmp(PfM, qpStart, qpEnd, MPI_COMM_SELF);
    }

    LogProjVal(eleProjCnt);
    /* calculate reweight */
    // TODO: check ~ Is it OK to fix w ?
    //w = exp(2.0*(log(fabs(ip))+x) - logSqPfFullSlater[sample]);
    w = 1.0;
    /*
       if(log(fabs(1.0-w)) > -5) {
       printf("w=%.3e\n",w);
       }
       */
    if (!isfinite(w)) {
      fprintf(stderr, "waring: VMCMainCal rank:%d sample:%d w=%e\n", rank, sample, w);
      continue;
    }

    StartTimer(41);
    if (AllComplexFlag == 0) {
      e = CalculateHamiltonianBF_real(creal(ip), eleIdx, eleCfg, eleNum, eleProjCnt, eleProjBFCnt);
    } else {/* calculate energy */
      e = CalculateHamiltonianBF_fcmp(ip, eleIdx, eleCfg, eleNum, eleProjCnt, eleProjBFCnt);
    }

    /* calculate double occupation D */
    db = CalculateDoubleOccupation(eleIdx, eleCfg, eleNum, eleProjCnt);
    StopTimer(41);
    if (! (isfinite(creal(e)) && isfinite(cimag(e)))) {
      fprintf(stderr, "waring: VMCMainCal rank:%d sample:%d e=%e\n", rank, sample, creal(e));
      continue;
    }

    if(NVMCCalMode == 1 && NLanczosMode == 1) {
      int lanczosInfo;
      lanczosCheckedSamples++;
      StartTimer(43);
      if(AllComplexFlag == 0) {
        lanczosInfo = LSLocalQBF_real(
            creal(e), creal(ip), eleIdx, eleCfg, eleNum,
            eleProjCnt, eleProjBFCnt, &lanczosScratch, LSLQ_real);
      } else {
        lanczosInfo = LSLocalQBF(
            e, ip, eleIdx, eleCfg, eleNum,
            eleProjCnt, eleProjBFCnt, &lanczosScratch, LSLQ);
      }
      StopTimer(43);
      if(lanczosInfo == LSLANCZOS_NUMERIC_REJECT) {
        lanczosRejectedSamples++;
        if(lanczosWarnings < 3) {
          fprintf(stderr,
                  "warning: BackFlow Lanczos rejected numerical candidate "
                  "rebuild rank:%d sample:%d.\n", rank, sample);
          lanczosWarnings++;
        }
        continue;
      }
      if(lanczosInfo != LSLANCZOS_OK) {
        fprintf(stderr,
                "Error: BackFlow Lanczos state or contract failure rank:%d "
                "sample:%d status:%d.\n", rank, sample, lanczosInfo);
        MPI_Abort(comm, EXIT_FAILURE);
      }
      if(sample == lanczosInjectSample) {
        if(AllComplexFlag == 0) LSLQ_real[3] = NAN;
        else LSLQ[3] = NAN + 0.0*I;
      }
      if(!lsbfLocalVectorFinite()) {
        lanczosRejectedSamples++;
        if(lanczosWarnings < 3) {
          fprintf(stderr,
                  "warning: BackFlow Lanczos rejected non-finite local vector rank:%d sample:%d.\n",
                  rank, sample);
          lanczosWarnings++;
        }
        continue;
      }
      if(lanczosOracleDump != NULL) {
        fprintf(lanczosOracleDump, "sample %d occ", sample);
        for(i=0; i<Nsite2; i++) fprintf(lanczosOracleDump, " %d", eleNum[i]);
        fprintf(lanczosOracleDump, " lslq");
        if(AllComplexFlag == 0) {
          for(i=0; i<NLSHam*NLSHam; i++) {
            fprintf(lanczosOracleDump, " %.17e %.17e", LSLQ_real[i], 0.0);
          }
        } else {
          for(i=0; i<NLSHam*NLSHam; i++) {
            fprintf(lanczosOracleDump, " %.17e %.17e",
                    creal(LSLQ[i]), cimag(LSLQ[i]));
          }
        }
        fprintf(lanczosOracleDump, "\n");
      }
    }

    Wc += w;
    Etot += w * e;
    Etot2 += w * conj(e) * e;
    Dbtot += w * db;
    Dbtot2 += w * db * db;

    if (NVMCCalMode == 0) {
      /* Calculate O for correlation fauctors */
      /*SROptO[0] = 1.0;
#pragma loop noalias
for(i=0;i<nProj;i++) srOptO[i+1] = (double)(eleProjCnt[i]);
*/
      srOptO[0] = 1.0 + 0.0 * I;//   real
      srOptO[1] = 0.0 + 0.0 * I;//   real
#pragma loop noalias
      for (i = 0; i < nProj; i++) {
        srOptO[(i + 1) * 2] = (double) (eleProjCnt[i]); // even real
        srOptO[(i + 1) * 2 + 1] = 0.0 + 0.0 * I;               // odd  comp
      }

      //StartTimer(44);
      /* BackflowDiff */
      //BackFlowDiff(SROptO+NProj+1,ip,eleIdx,eleNum,eleProjCnt,eleProjBFCnt);
      BackFlowDiff_fcmp(SROptO + 2 * NProj + 2, ip, eleIdx, eleNum, eleProjCnt, eleProjBFCnt);
      //StopTimer(44);

      StartTimer(42);
      /* SlaterElmDiff */
      //SlaterElmBFDiff_fcmp(SROptO+NProj+NProjBF+1,ip,eleIdx,eleNum,eleCfg,eleProjCnt,eleProjBFCnt);
      SlaterElmBFDiff_fcmp(SROptO + 2 * NProj + 2 * NProjBF + 2, ip, eleIdx, eleNum, eleCfg, eleProjCnt,
          eleProjBFCnt);
      StopTimer(42);

      if (FlagOptTrans > 0) {
        calculateOptTransDiff(SROptO + 2 * NProj + 2 * NProjBF + 2 * NSlater + 2, ip);
      }

      //[s] this part will be used for real varaibles
      if (AllComplexFlag == 0) {
#pragma loop noalias
        for (i = 0; i < SROptSize; i++) {
          srOptO_real[i] = creal(srOptO[2 * i]);
        }
      }
      //[e]

      StartTimer(43);
      /* Calculate OO and HO */
      if (NSRCG==0 && NStoreO == 0) {
        if (AllComplexFlag == 0) {
          calculateOO_real(SROptOO_real, SROptHO_real, SROptO_real, w, creal(e), SROptSize);
        } else {
          calculateOO(SROptOO, SROptHO, SROptO, w, e, SROptSize);
        }
      } else {
        we = w * e;
        sqrtw = sqrt(w);
        if (AllComplexFlag == 0) {
#pragma omp parallel for default(shared) private(int_i)
          for (int_i = 0; int_i < SROptSize; int_i++) {
            // SROptO_Store for fortran
            SROptO_Store_real[int_i + sample * SROptSize] = sqrtw * SROptO_real[int_i];
            SROptHO_real[int_i] += creal(we) * SROptO_real[int_i];
          }
        } else {
#pragma omp parallel for default(shared) private(int_i)
          for (int_i = 0; int_i < SROptSize * 2; int_i++) {
            // SROptO_Store for fortran
            SROptO_Store[int_i + sample * (2 * SROptSize)] = sqrtw * SROptO[int_i];
            SROptHO[int_i] += we * SROptO[int_i];
          }
        }
      }
      StopTimer(43);

    } else if (NVMCCalMode == 1) {
      StartTimer(42);
      /* Calculate Green Function */
      CalculateGreenFuncBF(w, ip, eleIdx, eleCfg, eleNum, eleProjCnt, eleProjBFCnt);
      StopTimer(42);
      if (NLanczosMode == 1) {
        StartTimer(43);
        if (AllComplexFlag == 0) {
          calculateQQQQ_real(QQQQ_real, LSLQ_real, w, NLSHam);
        } else {
          calculateQQQQ(QQQQ, LSLQ, w, NLSHam);
        }
        StopTimer(43);
      }
    }
  } /* end of for(sample) */

  if(NVMCCalMode == 1 && NLanczosMode == 1) {
    lsbfFinalizeAccounting(comm_parent, parentRank, lanczosCheckedSamples,
                           lanczosRejectedSamples);
  }
  if(lanczosOracleDump != NULL) fclose(lanczosOracleDump);
  if(lanczosScratchReady) LSLanczosBFScratchFree(&lanczosScratch);

  // calculate OO and HO at NVMCCalMode==0
  if(NVMCCalMode==0){
    if(NStoreO!=0 || NSRCG!=0){
      sampleSize=sampleEnd-sampleStart;
      if(NSRCG!=0 || sampleSize>0){
        if(AllComplexFlag==0){
          double *srOptO_Store_ptr = SROptO_Store_real;
          if(NSRCG==0) {
            srOptO_Store_ptr += (size_t)sampleStart * (size_t)SROptSize;
          }
          StartTimer(45);
          calculateOO_Store_real(SROptOO_real,SROptHO_real,srOptO_Store_ptr,creal(w),creal(e),SROptSize,sampleSize);
          StopTimer(45);
        }else{
          double complex *srOptO_Store_ptr = SROptO_Store;
          if(NSRCG==0) {
            srOptO_Store_ptr += (size_t)sampleStart * (2 * (size_t)SROptSize);
          }
          StartTimer(45);
          calculateOO_Store(SROptOO,SROptHO,srOptO_Store_ptr,w,e,2*SROptSize,sampleSize);
          StopTimer(45);
        }
      }
    }
  }

  InvM = InvM_Moto;
  PfM = PfM_Moto;

  return;
}
