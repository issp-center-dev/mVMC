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
 * timer program
 * "-lrt" option is needed for clock_gettime().
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include <time.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include "setmemory.h"
#ifndef _SRC_TIME
#define _SRC_TIME

static void AbortInvalidBFNBodyEnvironment(const char *name) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) {
    fprintf(stderr, "Error: invalid %s environment value.\n", name);
  }
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

static int ParseBFNBodyBooleanEnvironment(const char *name) {
  const char *value = getenv(name);
  if(value == NULL) return 0;
  if(strcmp(value, "0") == 0) return 0;
  if(strcmp(value, "1") == 0) return 1;
  AbortInvalidBFNBodyEnvironment(name);
  return 0;
}

static int ParseBFNBodyInjectionTerm(const char *value) {
  char *end = NULL;
  long parsed;
  if(value == NULL || value[0] < '0' || value[0] > '9') {
    AbortInvalidBFNBodyEnvironment("MVMC_BF_NBODY_INJECT_TERM");
    return -1;
  }
  errno = 0;
  parsed = strtol(value, &end, 10);
  if(errno != 0 || end == value || *end != '\0'
     || parsed < 0 || parsed > INT_MAX) {
    AbortInvalidBFNBodyEnvironment("MVMC_BF_NBODY_INJECT_TERM");
    return -1;
  }
  return (int)parsed;
}

void OutputTime(int step) {
  time_t tx;
  double pHop,pEx,pLSF;

  tx = time(NULL);
  if(step==0) {
    fprintf(FileTime, "%05d  acc_hop acc_ex  acc_lsf n_hop    n_ex      n_lsf   : %s", step, ctime(&tx));
  } else {
    pHop = (Counter[0] == 0) ? 0.0 : (double)Counter[1] / (double)Counter[0];
    pEx  = (Counter[2] == 0) ? 0.0 : (double)Counter[3] / (double)Counter[2];
    pLSF = (Counter[4] == 0) ? 0.0 : (double)Counter[5] / (double)Counter[4];
    fprintf(FileTime, "%05d  %.5lf %.5lf %.5lf %-8d %-8d  %-8d: %s", step, pHop,pEx,pLSF,
            Counter[0], Counter[2],Counter[4], ctime(&tx));
  }
}

double BFFSZC2DetailMonotonicSeconds(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC,&ts);
  return ts.tv_sec + ts.tv_nsec*1.0e-9;
}

static unsigned long long BFFSZC2ReuseCensusHash(
    int XI, int XJ, int XK, int XL) {
  const int values[4] = {XI,XJ,XK,XL};
  unsigned long long hash = 1469598103934665603ULL;
  int i;
  for(i=0;i<4;i++) {
    hash ^= (unsigned int)values[i];
    hash *= 1099511628211ULL;
  }
  return hash;
}

int GetBFFSZC2ReuseCensusWorkSize(
    long long maxCalls, int *capacity, int *intWorkSize) {
  long long required, value;
  if(capacity == NULL || intWorkSize == NULL || maxCalls < 0) return -1;
  *capacity = 0;
  *intWorkSize = 0;
  if(maxCalls == 0) return 0;
  if(maxCalls > INT_MAX/8) return -1;
  required = 2*maxCalls;
  value = 1;
  while(value < required) value *= 2;
  if(value > INT_MAX/4) return -1;
  *capacity = (int)value;
  *intWorkSize = 4*(int)value;
  return 0;
}

int InitBFFSZC2ReuseCensus(
    BFFSZC2ReuseCensus *census, int *keys, int capacity) {
  int i;
  if(census == NULL || capacity < 0
      || (capacity > 0 && (keys == NULL || (capacity & (capacity-1)) != 0))) {
    return -1;
  }
  census->keys = keys;
  census->capacity = capacity;
  census->trueCalls = 0;
  census->uniqueExactOrderedMoves = 0;
  census->duplicateTrueCalls = 0;
  census->overflowCalls = 0;
  for(i=0;i<capacity;i++) keys[4*i] = -1;
  return 0;
}

int RecordBFFSZC2ReuseCensus(
    BFFSZC2ReuseCensus *census, int XI, int XJ, int XK, int XL) {
  int probe;
  unsigned long long slot;
  int *key;
  if(census == NULL || census->capacity <= 0 || census->keys == NULL
      || (census->capacity & (census->capacity-1)) != 0
      || XI < 0 || XJ < 0 || XK < 0 || XL < 0) {
    return -1;
  }
  census->trueCalls++;
  slot = BFFSZC2ReuseCensusHash(XI,XJ,XK,XL)
      & (unsigned long long)(census->capacity-1);
  for(probe=0;probe<census->capacity;probe++) {
    key = census->keys+4*(int)slot;
    if(key[0] < 0) {
      key[0] = XI;
      key[1] = XJ;
      key[2] = XK;
      key[3] = XL;
      census->uniqueExactOrderedMoves++;
      return 0;
    }
    if(key[0] == XI && key[1] == XJ && key[2] == XK && key[3] == XL) {
      census->duplicateTrueCalls++;
      return 0;
    }
    slot = (slot+1) & (unsigned long long)(census->capacity-1);
  }
  census->overflowCalls++;
  return -1;
}

void MergeBFFSZC2ReuseCensus(
    int scope, const BFFSZC2ReuseCensus *census) {
  if(!BFFSZC2ReuseCensusEnabled || census == NULL
      || scope < 0 || scope >= NBFFSZC2ReuseScope) return;
  BFFSZC2ReuseCensusInvocations[scope]++;
  BFFSZC2ReuseCensusTrueCalls[scope] += census->trueCalls;
  BFFSZC2ReuseCensusUniqueExactOrderedMoves[scope]
      += census->uniqueExactOrderedMoves;
  BFFSZC2ReuseCensusDuplicateTrueCalls[scope] += census->duplicateTrueCalls;
  BFFSZC2ReuseCensusOverflowCalls[scope] += census->overflowCalls;
}

static int CheckBFFSZC2ReuseCensus(void) {
  BFFSZC2ReuseCensus census;
  int keys[32];
  int first[8];
  int a,slot,firstValue=-1,secondValue=-1;
  for(a=0;a<8;a++) first[a] = -1;
  for(a=0;a<9;a++) {
    slot = (int)(BFFSZC2ReuseCensusHash(a,a+1,a+2,a+3) & 7ULL);
    if(first[slot] >= 0) {
      firstValue = first[slot];
      secondValue = a;
      break;
    }
    first[slot] = a;
  }
  if(firstValue < 0 || InitBFFSZC2ReuseCensus(&census,keys,8) != 0
      || RecordBFFSZC2ReuseCensus(&census,firstValue,firstValue+1,
          firstValue+2,firstValue+3) != 0
      || RecordBFFSZC2ReuseCensus(&census,secondValue,secondValue+1,
          secondValue+2,secondValue+3) != 0
      || RecordBFFSZC2ReuseCensus(&census,firstValue,firstValue+1,
          firstValue+2,firstValue+3) != 0
      || census.trueCalls != 3 || census.uniqueExactOrderedMoves != 2
      || census.duplicateTrueCalls != 1 || census.overflowCalls != 0) {
    return -1;
  }
  if(InitBFFSZC2ReuseCensus(&census,keys,8) != 0) return -1;
  for(a=0;a<4;a++) {
    if(RecordBFFSZC2ReuseCensus(&census,a,10+a,20+a,30+a) != 0) return -1;
  }
  if(census.trueCalls != 4 || census.uniqueExactOrderedMoves != 4
      || census.duplicateTrueCalls != 0 || census.overflowCalls != 0) {
    return -1;
  }
  return 0;
}

void InitBFFSZC2DetailContext(BFFSZC2DetailContext *context, int source) {
  int i;
  if(context == NULL) return;
  context->source = source;
  context->evaluatedCalls = 0;
  context->changedSum = 0;
  context->changedMax = 0;
  context->affectedSum = 0;
  context->affectedMax = 0;
  context->affectedAtOrAboveKFull = 0;
  context->orderedDescriptorTotal = 0;
  context->reuseCensus = NULL;
  for(i=0;i<NBFFSZC2DetailClass;i++) context->classCall[i] = 0;
  for(i=0;i<NBFFSZC2DetailOutcome;i++) context->outcome[i] = 0;
  for(i=0;i<NBFFSZC2DetailPath;i++) context->path[i] = 0;
  for(i=0;i<NBFFSZProfileHist;i++) {
    context->changedHist[i] = 0;
    context->affectedHist[i] = 0;
  }
  for(i=0;i<NBFFSZC2DetailComponent;i++) context->componentSeconds[i] = 0.0;
}

void MergeBFFSZC2DetailContext(const BFFSZC2DetailContext *context) {
  int i;
  const int source = (context == NULL) ? -1 : context->source;
  if(!BFFSZC2DetailProfileEnabled || source < 0
      || source >= NBFFSZC2DetailSource) return;
  for(i=0;i<NBFFSZC2DetailClass;i++) {
    BFFSZC2DetailClassCall[source][i] += context->classCall[i];
  }
  for(i=0;i<NBFFSZC2DetailOutcome;i++) {
    BFFSZC2DetailOutcome[source][i] += context->outcome[i];
  }
  for(i=0;i<NBFFSZC2DetailPath;i++) {
    BFFSZC2DetailPath[source][i] += context->path[i];
  }
  BFFSZC2DetailEvaluatedCalls[source] += context->evaluatedCalls;
  BFFSZC2DetailChangedSum[source] += context->changedSum;
  BFFSZC2DetailAffectedSum[source] += context->affectedSum;
  if(context->changedMax > BFFSZC2DetailChangedMax[source]) {
    BFFSZC2DetailChangedMax[source] = context->changedMax;
  }
  if(context->affectedMax > BFFSZC2DetailAffectedMax[source]) {
    BFFSZC2DetailAffectedMax[source] = context->affectedMax;
  }
  for(i=0;i<NBFFSZProfileHist;i++) {
    BFFSZC2DetailChangedHist[source][i] += context->changedHist[i];
    BFFSZC2DetailAffectedHist[source][i] += context->affectedHist[i];
  }
  BFFSZC2DetailAffectedAtOrAboveKFull[source]
      += context->affectedAtOrAboveKFull;
  BFFSZC2DetailOrderedDescriptorTotal[source]
      += context->orderedDescriptorTotal;
  for(i=0;i<NBFFSZC2DetailComponent;i++) {
    BFFSZC2DetailComponentSeconds[source][i] += context->componentSeconds[i];
  }
}

void InitBFFSZC2DetailTermContext(BFFSZC2DetailTermContext *context) {
  int i;
  if(context == NULL) return;
  for(i=0;i<NBFFSZC2DetailTerm;i++) context->seconds[i] = 0.0;
}

void MergeBFFSZC2DetailTermContext(const BFFSZC2DetailTermContext *context) {
  int i;
  if(!BFFSZC2DetailProfileEnabled || context == NULL) return;
  for(i=0;i<NBFFSZC2DetailTerm;i++) {
    BFFSZC2DetailTermSeconds[i] += context->seconds[i];
  }
}

void InitTimer() {
  int i, j;
  const char *bfProfileEnv, *greenCheckEnv, *affectedCheckEnv;
  const char *pfUpdateCheckEnv, *pfUpdateFallbackEnv, *pfUpdateInjectEnv;
  const char *pfUpdateInjectRankEnv;
  const char *pfUpdateExplicitStateEnv, *pfUpdateArgumentEnv, *pfUpdateKFullEnv;
  const char *permuteParticleLabelsEnv;
  const char *invUpdateCheckEnv, *invUpdateFallbackEnv, *invUpdateInjectEnv;
  const char *invUpdateInjectRankEnv, *invUpdateExplicitStateEnv;
  const char *invUpdateArgumentEnv, *invGemmCheckEnv, *invDetailProfileEnv;
  const char *matrixFreeCheckEnv, *matrixFreeArgumentEnv, *multiMoveRowCheckEnv;
  const char *c2StateCheckEnv, *c2BufferCheckEnv, *c2ForceSectorZeroEnv;
  const char *samplingRejectCheckEnv, *c2DetailProfileEnv, *c2ReuseCensusEnv;
  const char *nbodyInjectStageEnv, *nbodyInjectTermEnv;
  for(i=0;i<NTimer;i++) Timer[i]=0.0;
  for(i=0;i<NTimer;i++) TimerStart[i]=0.0;
  for(i=0;i<NBFProfileCounter;i++) BFProfileCounter[i]=0;
  for(i=0;i<NBFFSZProfileSource;i++) {
    BFFSZProfileCall[i] = 0;
    BFFSZProfileChangedSum[i] = 0;
    BFFSZProfileChangedMax[i] = 0;
    BFFSZProfileAffectedSum[i] = 0;
    BFFSZProfileAffectedMax[i] = 0;
    for(j=0;j<NBFFSZProfileHist;j++) {
      BFFSZProfileChangedHist[i][j] = 0;
      BFFSZProfileAffectedHist[i][j] = 0;
    }
    for(j=0;j<NBFFSZProfileRatioHist;j++) {
      BFFSZProfileAffectedRatioHist[i][j] = 0;
    }
    for(j=0;j<NBFFSZPfPath;j++) BFFSZProfilePfPath[i][j] = 0;
  }
  for(j=0;j<NBFFSZPfPath;j++) BFFSZProfileInvPath[j] = 0;
  for(j=0;j<NBFFSZGreenMaterialize;j++) {
    BFFSZProfileGreenMaterialize[j] = 0;
  }
  for(j=0;j<NBFFSZSampleMaterialize;j++) {
    BFFSZProfileSampleMaterialize[j] = 0;
  }
  for(j=0;j<NBFFSZPfPath;j++) BFFSZProfileSampleCommit[j] = 0;
  BFFSZProfileInvCheckSeconds = 0.0;
  BFFSZProfileInvAntisymmetryMax = 0.0;
  BFFSZProfileInvResidualMax = 0.0;
  BFFSZProfileInvGemmCheckCount = 0;
  BFFSZProfileInvGemmDifferenceMax = 0.0;
  for(j=0;j<NBFFSZInvDetail;j++) {
    BFFSZProfileInvDetailSeconds[j] = 0.0;
    BFFSZProfileInvDetailMaxSeconds[j] = 0.0;
  }
  for(i=0;i<NBFFSZC2DetailSource;i++) {
    BFFSZC2DetailEvaluatedCalls[i] = 0;
    BFFSZC2DetailChangedSum[i] = 0;
    BFFSZC2DetailChangedMax[i] = 0;
    BFFSZC2DetailAffectedSum[i] = 0;
    BFFSZC2DetailAffectedMax[i] = 0;
    BFFSZC2DetailAffectedAtOrAboveKFull[i] = 0;
    BFFSZC2DetailOrderedDescriptorTotal[i] = 0;
    for(j=0;j<NBFFSZC2DetailClass;j++) BFFSZC2DetailClassCall[i][j] = 0;
    for(j=0;j<NBFFSZC2DetailOutcome;j++) BFFSZC2DetailOutcome[i][j] = 0;
    for(j=0;j<NBFFSZC2DetailPath;j++) BFFSZC2DetailPath[i][j] = 0;
    for(j=0;j<NBFFSZProfileHist;j++) {
      BFFSZC2DetailChangedHist[i][j] = 0;
      BFFSZC2DetailAffectedHist[i][j] = 0;
    }
    for(j=0;j<NBFFSZC2DetailComponent;j++) {
      BFFSZC2DetailComponentSeconds[i][j] = 0.0;
    }
  }
  for(i=0;i<NBFFSZC2DetailTerm;i++) BFFSZC2DetailTermSeconds[i] = 0.0;
  BFFSZC2DetailSzSecondsRank0 = 0.0;
  for(i=0;i<NBFFSZC2ReuseScope;i++) {
    BFFSZC2ReuseCensusInvocations[i] = 0;
    BFFSZC2ReuseCensusTrueCalls[i] = 0;
    BFFSZC2ReuseCensusUniqueExactOrderedMoves[i] = 0;
    BFFSZC2ReuseCensusDuplicateTrueCalls[i] = 0;
    BFFSZC2ReuseCensusOverflowCalls[i] = 0;
  }
  BFNBodyStateCheckEnabled =
      ParseBFNBodyBooleanEnvironment("MVMC_BF_NBODY_STATE_CHECK");
  BFNBodyInjectStage = BF_NBODY_INJECT_NONE;
  BFNBodyInjectTerm = -1;
  BFNBodyStateCheckCount = 0;
  BFNBodyStateCheckFailureCount = 0;
  nbodyInjectStageEnv = getenv("MVMC_BF_NBODY_INJECT_STAGE");
  nbodyInjectTermEnv = getenv("MVMC_BF_NBODY_INJECT_TERM");
  if((nbodyInjectStageEnv == NULL) != (nbodyInjectTermEnv == NULL)) {
    AbortInvalidBFNBodyEnvironment(
        nbodyInjectStageEnv == NULL
        ? "MVMC_BF_NBODY_INJECT_STAGE"
        : "MVMC_BF_NBODY_INJECT_TERM");
  }
  if(nbodyInjectStageEnv != NULL) {
    if(strcmp(nbodyInjectStageEnv, "workspace") == 0) {
      BFNBodyInjectStage = BF_NBODY_INJECT_WORKSPACE;
    } else if(strcmp(nbodyInjectStageEnv, "candidate") == 0) {
      BFNBodyInjectStage = BF_NBODY_INJECT_CANDIDATE;
    } else if(strcmp(nbodyInjectStageEnv, "pfaffian") == 0) {
      BFNBodyInjectStage = BF_NBODY_INJECT_PFAFFIAN;
    } else {
      AbortInvalidBFNBodyEnvironment("MVMC_BF_NBODY_INJECT_STAGE");
    }
    BFNBodyInjectTerm = ParseBFNBodyInjectionTerm(nbodyInjectTermEnv);
  }
  bfProfileEnv = getenv("MVMC_BF_PROFILE");
  BFProfileEnabled = (bfProfileEnv != NULL && atoi(bfProfileEnv) != 0);
  greenCheckEnv = getenv("MVMC_BF_FSZ_GREEN_REBUILD_CHECK");
  BFFSZGreenRebuildCheckEnabled = (greenCheckEnv != NULL && atoi(greenCheckEnv) != 0);
  affectedCheckEnv = getenv("MVMC_BF_FSZ_AFFECTED_CHECK");
  BFFSZAffectedCheckEnabled = (affectedCheckEnv != NULL && atoi(affectedCheckEnv) != 0);
  pfUpdateCheckEnv = getenv("MVMC_BF_FSZ_PF_UPDATE_CHECK");
  BFFSZPfUpdateCheckEnabled = (pfUpdateCheckEnv != NULL && atoi(pfUpdateCheckEnv) != 0);
  pfUpdateFallbackEnv = getenv("MVMC_BF_FSZ_PF_UPDATE_FORCE_FALLBACK");
  BFFSZPfUpdateForceFallback = (pfUpdateFallbackEnv != NULL
      && atoi(pfUpdateFallbackEnv) != 0);
  BFFSZPfUpdateInjectedStatus = BF_FSZ_PF_UPDATE_OK;
  pfUpdateInjectEnv = getenv("MVMC_BF_FSZ_PF_UPDATE_INJECT_STATUS");
  if(pfUpdateInjectEnv != NULL) {
    const int value = atoi(pfUpdateInjectEnv);
    if(value >= BF_FSZ_PF_UPDATE_EXACT_ZERO
        && value <= BF_FSZ_PF_UPDATE_NONFINITE) {
      BFFSZPfUpdateInjectedStatus = value;
    }
  }
  BFFSZPfUpdateInjectedRank = -1;
  pfUpdateInjectRankEnv = getenv("MVMC_BF_FSZ_PF_UPDATE_INJECT_RANK");
  if(pfUpdateInjectRankEnv != NULL) {
    int rank = 0;
    BFFSZPfUpdateInjectedRank = atoi(pfUpdateInjectRankEnv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(BFFSZPfUpdateInjectedRank >= 0 && rank != BFFSZPfUpdateInjectedRank) {
      BFFSZPfUpdateInjectedStatus = BF_FSZ_PF_UPDATE_OK;
    }
  }
  pfUpdateExplicitStateEnv = getenv("MVMC_BF_FSZ_PF_UPDATE_EXPLICIT_STATE_CHECK");
  BFFSZPfUpdateExplicitStateCheckEnabled = (pfUpdateExplicitStateEnv != NULL
      && atoi(pfUpdateExplicitStateEnv) != 0);
  pfUpdateArgumentEnv = getenv("MVMC_BF_FSZ_PF_UPDATE_ARGUMENT_CHECK");
  BFFSZPfUpdateArgumentCheckEnabled = (pfUpdateArgumentEnv != NULL
      && atoi(pfUpdateArgumentEnv) != 0);
  BFFSZPfUpdateKFull = BF_FSZ_PF_UPDATE_KFULL_DEFAULT;
  pfUpdateKFullEnv = getenv("MVMC_BF_FSZ_PF_UPDATE_KFULL");
  if(pfUpdateKFullEnv != NULL) {
    const int value = atoi(pfUpdateKFullEnv);
    if(value >= 1 && value <= BF_FSZ_PF_UPDATE_KFULL_DEFAULT) {
      BFFSZPfUpdateKFull = value;
    }
  }
  permuteParticleLabelsEnv = getenv("MVMC_BF_FSZ_PERMUTE_PARTICLE_LABELS");
  BFFSZPermuteParticleLabelsCheckEnabled = (permuteParticleLabelsEnv != NULL
      && atoi(permuteParticleLabelsEnv) != 0);
  invUpdateCheckEnv = getenv("MVMC_BF_FSZ_INV_UPDATE_CHECK");
  BFFSZInvUpdateCheckEnabled = (invUpdateCheckEnv != NULL
      && atoi(invUpdateCheckEnv) != 0);
  invUpdateFallbackEnv = getenv("MVMC_BF_FSZ_INV_UPDATE_FORCE_FALLBACK");
  BFFSZInvUpdateForceFallback = (invUpdateFallbackEnv != NULL
      && atoi(invUpdateFallbackEnv) != 0);
  BFFSZInvUpdateInjectedStage = BF_FSZ_INV_STAGE_NONE;
  invUpdateInjectEnv = getenv("MVMC_BF_FSZ_INV_UPDATE_INJECT_STAGE");
  if(invUpdateInjectEnv != NULL) {
    const int value = atoi(invUpdateInjectEnv);
    if(value >= BF_FSZ_INV_STAGE_GETRF && value <= BF_FSZ_INV_STAGE_RESIDUAL) {
      BFFSZInvUpdateInjectedStage = value;
    }
  }
  BFFSZInvUpdateInjectedRank = -1;
  invUpdateInjectRankEnv = getenv("MVMC_BF_FSZ_INV_UPDATE_INJECT_RANK");
  if(invUpdateInjectRankEnv != NULL) {
    BFFSZInvUpdateInjectedRank = atoi(invUpdateInjectRankEnv);
  }
  invUpdateExplicitStateEnv = getenv("MVMC_BF_FSZ_INV_UPDATE_EXPLICIT_STATE_CHECK");
  BFFSZInvUpdateExplicitStateCheckEnabled = (invUpdateExplicitStateEnv != NULL
      && atoi(invUpdateExplicitStateEnv) != 0);
  invUpdateArgumentEnv = getenv("MVMC_BF_FSZ_INV_UPDATE_ARGUMENT_CHECK");
  BFFSZInvUpdateArgumentCheckEnabled = (invUpdateArgumentEnv != NULL
      && atoi(invUpdateArgumentEnv) != 0);
  invGemmCheckEnv = getenv("MVMC_BF_FSZ_INV_GEMM_CHECK");
  BFFSZInvGemmCheckEnabled = (invGemmCheckEnv != NULL
      && atoi(invGemmCheckEnv) != 0);
  invDetailProfileEnv = getenv("MVMC_BF_FSZ_INV_DETAIL_PROFILE");
  BFFSZInvDetailProfileEnabled = (invDetailProfileEnv != NULL
      && atoi(invDetailProfileEnv) != 0);
  c2DetailProfileEnv = getenv("MVMC_BF_FSZ_C2_DETAIL_PROFILE");
  BFFSZC2DetailProfileEnabled = (c2DetailProfileEnv != NULL
      && atoi(c2DetailProfileEnv) != 0);
  c2ReuseCensusEnv = getenv("MVMC_BF_FSZ_C2_REUSE_CENSUS");
  BFFSZC2ReuseCensusEnabled = BFFSZC2DetailProfileEnabled
      && c2ReuseCensusEnv != NULL && atoi(c2ReuseCensusEnv) != 0;
  if(BFFSZC2ReuseCensusEnabled && CheckBFFSZC2ReuseCensus() != 0) {
    fprintf(stderr,"error: BF-FSZ C2 reuse census self-check failed\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  matrixFreeCheckEnv = getenv("MVMC_BF_FSZ_MATRIX_FREE_CHECK");
  BFFSZMatrixFreeCheckEnabled = (matrixFreeCheckEnv != NULL
      && atoi(matrixFreeCheckEnv) != 0);
  matrixFreeArgumentEnv = getenv("MVMC_BF_FSZ_MATRIX_FREE_ARGUMENT_CHECK");
  BFFSZMatrixFreeArgumentCheckEnabled = (matrixFreeArgumentEnv != NULL
      && atoi(matrixFreeArgumentEnv) != 0);
  if(BFFSZMatrixFreeArgumentCheckEnabled) {
    BFFSZMatrixFreeCheckEnabled = 1;
  }
  multiMoveRowCheckEnv = getenv("MVMC_BF_FSZ_MULTI_MOVE_ROW_CHECK");
  BFFSZMultiMoveRowCheckEnabled = (multiMoveRowCheckEnv != NULL
      && atoi(multiMoveRowCheckEnv) != 0);
  c2StateCheckEnv = getenv("MVMC_BF_FSZ_C2_STATE_CHECK");
  BFFSZC2StateCheckEnabled = (c2StateCheckEnv != NULL
      && atoi(c2StateCheckEnv) != 0);
  c2BufferCheckEnv = getenv("MVMC_BF_FSZ_C2_BUFFER_CHECK");
  BFFSZC2BufferCheckEnabled = (c2BufferCheckEnv != NULL
      && atoi(c2BufferCheckEnv) != 0);
  c2ForceSectorZeroEnv = getenv("MVMC_BF_FSZ_C2_FORCE_SECTOR_ZERO");
  BFFSZC2ForceSectorZero = (c2ForceSectorZeroEnv != NULL
      && atoi(c2ForceSectorZeroEnv) != 0);
  samplingRejectCheckEnv = getenv("MVMC_BF_FSZ_SAMPLING_REJECT_CHECK");
  BFFSZSamplingRejectCheckEnabled = (samplingRejectCheckEnv != NULL
      && atoi(samplingRejectCheckEnv) != 0);
  return;
}

void StartTimer(int n) {
#ifdef _mpi_use
  TimerStart[n]=MPI_Wtime();
#else
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME,&ts);
  TimerStart[n]=ts.tv_sec + ts.tv_nsec*1.0e-9;
#endif
  return;
}

void StopTimer(int n) {
#ifdef _mpi_use
  Timer[n] += MPI_Wtime() - TimerStart[n];
#else
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME,&ts);
  Timer[n] += ts.tv_sec + ts.tv_nsec*1.0e-9 - TimerStart[n];
#endif
  return;
}

static void OutputBFProfileCounters(FILE *fp) {
  int i, hasData = 0;
  if(!BFProfileEnabled && !BFFSZInvGemmCheckEnabled
      && !BFFSZInvDetailProfileEnabled
      && !BFFSZC2DetailProfileEnabled) return;
  if(BFProfileEnabled) {
    for(i=0;i<NBFProfileCounter;i++) {
      if(BFProfileCounter[i] != 0) {
        hasData = 1;
        break;
      }
    }
    if(hasData || BFFSZProfileCall[BFFSZ_PROFILE_SAMPLE] != 0
        || BFFSZProfileCall[BFFSZ_PROFILE_GREEN] != 0) {
      fprintf(fp,"  BF profile counters (MVMC_BF_PROFILE=1)\n");
  fprintf(fp,"    BF sample row requests          [910] %12lld\n",BFProfileCounter[BFPROF_SAMPLE_ROW_REQUEST]);
  fprintf(fp,"    BF sample row recompute         [911] %12lld\n",BFProfileCounter[BFPROF_SAMPLE_ROW_RECOMPUTE]);
  fprintf(fp,"    BF sample row reuse             [912] %12lld\n",BFProfileCounter[BFPROF_SAMPLE_ROW_REUSE]);
  fprintf(fp,"    BF sample pair requests         [913] %12lld\n",BFProfileCounter[BFPROF_SAMPLE_PAIR_REQUEST]);
  fprintf(fp,"    BF sample pair recompute        [914] %12lld\n",BFProfileCounter[BFPROF_SAMPLE_PAIR_RECOMPUTE]);
  fprintf(fp,"    BF sample term dense candidates [950] %12lld\n",BFProfileCounter[BFPROF_SAMPLE_TERM_DENSE_CANDIDATE]);
  fprintf(fp,"    BF sample term geometry valid   [951] %12lld\n",BFProfileCounter[BFPROF_SAMPLE_TERM_GEOMETRY_VALID]);
  fprintf(fp,"    BF sample term sparse pairs     [952] %12lld\n",BFProfileCounter[BFPROF_SAMPLE_TERM_SPARSE_PAIR]);
  fprintf(fp,"    BF sample term actual adds      [953] %12lld\n",BFProfileCounter[BFPROF_SAMPLE_TERM_ACTUAL_ADD]);
  fprintf(fp,"    BF green row requests           [920] %12lld\n",BFProfileCounter[BFPROF_GREEN_ROW_REQUEST]);
  fprintf(fp,"    BF green row recompute          [921] %12lld\n",BFProfileCounter[BFPROF_GREEN_ROW_RECOMPUTE]);
  fprintf(fp,"    BF green row reuse              [922] %12lld\n",BFProfileCounter[BFPROF_GREEN_ROW_REUSE]);
  fprintf(fp,"    BF green pair requests          [923] %12lld\n",BFProfileCounter[BFPROF_GREEN_PAIR_REQUEST]);
  fprintf(fp,"    BF green pair recompute         [924] %12lld\n",BFProfileCounter[BFPROF_GREEN_PAIR_RECOMPUTE]);
  fprintf(fp,"    BF green term dense candidates  [960] %12lld\n",BFProfileCounter[BFPROF_GREEN_TERM_DENSE_CANDIDATE]);
  fprintf(fp,"    BF green term geometry valid    [961] %12lld\n",BFProfileCounter[BFPROF_GREEN_TERM_GEOMETRY_VALID]);
  fprintf(fp,"    BF green term sparse pairs      [962] %12lld\n",BFProfileCounter[BFPROF_GREEN_TERM_SPARSE_PAIR]);
  fprintf(fp,"    BF green term actual adds       [963] %12lld\n",BFProfileCounter[BFPROF_GREEN_TERM_ACTUAL_ADD]);
  fprintf(fp,"    BF sample delta cnt total       [970] %12lld\n",BFProfileCounter[BFPROF_SAMPLE_DELTA_CNT_TOTAL]);
  fprintf(fp,"    BF sample delta cnt changed     [971] %12lld\n",BFProfileCounter[BFPROF_SAMPLE_DELTA_CNT_CHANGED]);
  fprintf(fp,"    BF sample delta pair new        [972] %12lld\n",BFProfileCounter[BFPROF_SAMPLE_DELTA_PAIR_NEW]);
  fprintf(fp,"    BF sample delta pair old        [973] %12lld\n",BFProfileCounter[BFPROF_SAMPLE_DELTA_PAIR_OLD]);
  fprintf(fp,"    BF sample delta pair total      [974] %12lld\n",BFProfileCounter[BFPROF_SAMPLE_DELTA_PAIR_TOTAL]);
  fprintf(fp,"    BF green delta cnt total        [980] %12lld\n",BFProfileCounter[BFPROF_GREEN_DELTA_CNT_TOTAL]);
  fprintf(fp,"    BF green delta cnt changed      [981] %12lld\n",BFProfileCounter[BFPROF_GREEN_DELTA_CNT_CHANGED]);
  fprintf(fp,"    BF green delta pair new         [982] %12lld\n",BFProfileCounter[BFPROF_GREEN_DELTA_PAIR_NEW]);
  fprintf(fp,"    BF green delta pair old         [983] %12lld\n",BFProfileCounter[BFPROF_GREEN_DELTA_PAIR_OLD]);
  fprintf(fp,"    BF green delta pair total       [984] %12lld\n",BFProfileCounter[BFPROF_GREEN_DELTA_PAIR_TOTAL]);
  fprintf(fp,"    BF cnt snapshots                [930] %12lld\n",BFProfileCounter[BFPROF_BFCNT_SNAPSHOT]);
  fprintf(fp,"    BF cnt total entries            [931] %12lld\n",BFProfileCounter[BFPROF_BFCNT_TOTAL_ENTRY]);
  fprintf(fp,"    BF cnt group0 nnz               [932] %12lld\n",BFProfileCounter[BFPROF_BFCNT_GROUP0_NNZ]);
  fprintf(fp,"    BF cnt group1 nnz               [933] %12lld\n",BFProfileCounter[BFPROF_BFCNT_GROUP1_NNZ]);
  fprintf(fp,"    BF cnt group2 nnz               [934] %12lld\n",BFProfileCounter[BFPROF_BFCNT_GROUP2_NNZ]);
  fprintf(fp,"    BF cnt group3 nnz               [935] %12lld\n",BFProfileCounter[BFPROF_BFCNT_GROUP3_NNZ]);
  fprintf(fp,"    BF cnt state0 nnz               [936] %12lld\n",BFProfileCounter[BFPROF_BFCNT_STATE0_NNZ]);
  fprintf(fp,"    BF cnt state1 nnz               [937] %12lld\n",BFProfileCounter[BFPROF_BFCNT_STATE1_NNZ]);
  fprintf(fp,"    BF cnt state2 nnz               [938] %12lld\n",BFProfileCounter[BFPROF_BFCNT_STATE2_NNZ]);
  fprintf(fp,"    BF cnt state3 nnz               [939] %12lld\n",BFProfileCounter[BFPROF_BFCNT_STATE3_NNZ]);
  fprintf(fp,"    BF hop try                      [940] %12lld\n",BFProfileCounter[BFPROF_HOP_TRY]);
  fprintf(fp,"    BF hop candidate reject         [941] %12lld\n",BFProfileCounter[BFPROF_HOP_CANDIDATE_REJECT]);
  fprintf(fp,"    BF hop valid                    [942] %12lld\n",BFProfileCounter[BFPROF_HOP_VALID]);
  fprintf(fp,"    BF hop accept                   [943] %12lld\n",BFProfileCounter[BFPROF_HOP_ACCEPT]);
  fprintf(fp,"    BF hop metropolis reject        [944] %12lld\n",BFProfileCounter[BFPROF_HOP_METROPOLIS_REJECT]);
  fprintf(fp,"    BF exchange try                 [945] %12lld\n",BFProfileCounter[BFPROF_EXCHANGE_TRY]);
  fprintf(fp,"    BF exchange candidate reject    [946] %12lld\n",BFProfileCounter[BFPROF_EXCHANGE_CANDIDATE_REJECT]);
  fprintf(fp,"    BF exchange valid               [947] %12lld\n",BFProfileCounter[BFPROF_EXCHANGE_VALID]);
  fprintf(fp,"    BF exchange accept              [948] %12lld\n",BFProfileCounter[BFPROF_EXCHANGE_ACCEPT]);
  fprintf(fp,"    BF exchange metropolis reject   [949] %12lld\n",BFProfileCounter[BFPROF_EXCHANGE_METROPOLIS_REJECT]);

  for(i=0;i<NBFFSZProfileSource;i++) {
    const char *label = (i == BFFSZ_PROFILE_SAMPLE) ? "sample" : "green";
    int k;
    double changedMean, affectedMean;
    if(BFFSZProfileCall[i] == 0) continue;
    changedMean = (double)BFFSZProfileChangedSum[i] / (double)BFFSZProfileCall[i];
    affectedMean = (double)BFFSZProfileAffectedSum[i] / (double)BFFSZProfileCall[i];
    fprintf(fp,"    BF-FSZ %s changed stats         calls=%lld mean=%.6f max=%lld\n",
            label, BFFSZProfileCall[i], changedMean, BFFSZProfileChangedMax[i]);
    fprintf(fp,"    BF-FSZ %s affected stats        calls=%lld mean=%.6f max=%lld\n",
            label, BFFSZProfileCall[i], affectedMean, BFFSZProfileAffectedMax[i]);
    fprintf(fp,"    BF-FSZ %s changed hist          ",label);
    for(k=0;k<=16;k++) fprintf(fp," %d:%lld",k,BFFSZProfileChangedHist[i][k]);
    fprintf(fp," 17-32:%lld 33-64:%lld 65-128:%lld 129-256:%lld 257+:%lld\n",
            BFFSZProfileChangedHist[i][17],BFFSZProfileChangedHist[i][18],
            BFFSZProfileChangedHist[i][19],BFFSZProfileChangedHist[i][20],
            BFFSZProfileChangedHist[i][21]);
    fprintf(fp,"    BF-FSZ %s affected hist         ",label);
    for(k=0;k<=16;k++) fprintf(fp," %d:%lld",k,BFFSZProfileAffectedHist[i][k]);
    fprintf(fp," 17-32:%lld 33-64:%lld 65-128:%lld 129-256:%lld 257+:%lld\n",
            BFFSZProfileAffectedHist[i][17],BFFSZProfileAffectedHist[i][18],
            BFFSZProfileAffectedHist[i][19],BFFSZProfileAffectedHist[i][20],
            BFFSZProfileAffectedHist[i][21]);
    fprintf(fp,"    BF-FSZ %s affected ratio hist   0:%lld <=1/16:%lld <=1/8:%lld <=1/4:%lld <=1/2:%lld <=3/4:%lld <1:%lld =1:%lld\n",
            label,BFFSZProfileAffectedRatioHist[i][0],BFFSZProfileAffectedRatioHist[i][1],
            BFFSZProfileAffectedRatioHist[i][2],BFFSZProfileAffectedRatioHist[i][3],
            BFFSZProfileAffectedRatioHist[i][4],BFFSZProfileAffectedRatioHist[i][5],
            BFFSZProfileAffectedRatioHist[i][6],BFFSZProfileAffectedRatioHist[i][7]);
    fprintf(fp,"    BF-FSZ %s Pfaffian paths        optimized:%lld direct-full:%lld fallback:%lld kFull:%d\n",
            label,BFFSZProfilePfPath[i][BFFSZ_PF_PATH_OPTIMIZED],
            BFFSZProfilePfPath[i][BFFSZ_PF_PATH_DIRECT_FULL],
            BFFSZProfilePfPath[i][BFFSZ_PF_PATH_FALLBACK],BFFSZPfUpdateKFull);
  }
  fprintf(fp,"    BF-FSZ accepted inverse paths   optimized:%lld direct-full:%lld fallback:%lld kFull:%d\n",
          BFFSZProfileInvPath[BFFSZ_PF_PATH_OPTIMIZED],
          BFFSZProfileInvPath[BFFSZ_PF_PATH_DIRECT_FULL],
          BFFSZProfileInvPath[BFFSZ_PF_PATH_FALLBACK],BFFSZPfUpdateKFull);
  fprintf(fp,"    BF-FSZ green full materialize   direct-full:%lld fallback:%lld oracle:%lld\n",
          BFFSZProfileGreenMaterialize[BFFSZ_GREEN_MATERIALIZE_DIRECT_FULL],
          BFFSZProfileGreenMaterialize[BFFSZ_GREEN_MATERIALIZE_FALLBACK],
          BFFSZProfileGreenMaterialize[BFFSZ_GREEN_MATERIALIZE_ORACLE]);
  fprintf(fp,"    BF-FSZ sample full materialize  direct-full:%lld pf-fallback:%lld inv-fallback:%lld oracle:%lld\n",
          BFFSZProfileSampleMaterialize[BFFSZ_SAMPLE_MATERIALIZE_DIRECT_FULL],
          BFFSZProfileSampleMaterialize[BFFSZ_SAMPLE_MATERIALIZE_PF_FALLBACK],
          BFFSZProfileSampleMaterialize[BFFSZ_SAMPLE_MATERIALIZE_INV_FALLBACK],
          BFFSZProfileSampleMaterialize[BFFSZ_SAMPLE_MATERIALIZE_ORACLE]);
  fprintf(fp,"    BF-FSZ sample Slater commits    optimized:%lld direct-full:%lld fallback:%lld\n",
          BFFSZProfileSampleCommit[BFFSZ_PF_PATH_OPTIMIZED],
          BFFSZProfileSampleCommit[BFFSZ_PF_PATH_DIRECT_FULL],
          BFFSZProfileSampleCommit[BFFSZ_PF_PATH_FALLBACK]);
      fprintf(fp,"    BF-FSZ inverse checks           seconds_max_rank:%.9f antisymmetry_max:%.17e affected_residual_max:%.17e\n",
              BFFSZProfileInvCheckSeconds,BFFSZProfileInvAntisymmetryMax,
              BFFSZProfileInvResidualMax);
    }
  }
  if(BFFSZInvGemmCheckEnabled) {
    fprintf(fp,"  BF-FSZ inverse GEMM check (MVMC_BF_FSZ_INV_GEMM_CHECK=1)\n");
    fprintf(fp,"    BF-FSZ inverse GEMM checks      calls=%lld scaled_max=%.17e\n",
            BFFSZProfileInvGemmCheckCount,
            BFFSZProfileInvGemmDifferenceMax);
  }
  if(BFFSZInvDetailProfileEnabled) {
    const char *labels[NBFFSZInvDetail] = {
      "w",
      "small-transpose",
      "u",
      "lapack",
      "correction",
      "scan-antisymmetrize",
      "affected-residual",
      "mpi-agreement",
      "commit-copy"
    };
    double classifiedSeconds = 0.0;
    double unclassifiedSeconds;
    fprintf(fp,"  BF-FSZ inverse detail profile (MVMC_BF_FSZ_INV_DETAIL_PROFILE=1)\n");
    for(i=0;i<NBFFSZInvDetail;i++) {
      const double share = (Timer[63] > 0.0)
          ? BFFSZProfileInvDetailSeconds[i]/Timer[63] : 0.0;
      classifiedSeconds += BFFSZProfileInvDetailSeconds[i];
      fprintf(fp,"    BF-FSZ inverse detail component=%s seconds_rank0=%.9f seconds_max_rank=%.9f share_rank0=%.9f\n",
              labels[i],BFFSZProfileInvDetailSeconds[i],
              BFFSZProfileInvDetailMaxSeconds[i],share);
    }
    unclassifiedSeconds = Timer[63]-classifiedSeconds;
    fprintf(fp,"    BF-FSZ inverse detail total timer63_seconds_rank0=%.9f classified_seconds_rank0=%.9f unclassified_seconds_rank0=%.9f\n",
            Timer[63],classifiedSeconds,unclassifiedSeconds);
  }
  if(BFFSZC2DetailProfileEnabled) {
    const char *sourceLabels[NBFFSZC2DetailSource] = {
      "measurement", "pair_hop", "exchange", "inter_all"
    };
    const char *outcomeLabels[NBFFSZC2DetailOutcome] = {
      "occupancy_zero", "sector_zero", "evaluated"
    };
    const char *pathLabels[NBFFSZC2DetailPath] = {
      "legacy_full", "optimized_row", "direct_full", "fallback", "debug_oracle"
    };
    const char *componentLabels[NBFFSZC2DetailComponent] = {
      "dispatch", "state_projection", "bf_count", "candidate_build",
      "pfaffian", "restore", "affected_collect"
    };
    const char *termLabels[NBFFSZC2DetailTerm] = {
      "number", "transfer", "pair_hop", "exchange", "inter_all"
    };
    const char *reuseScopeLabels[NBFFSZC2ReuseScope] = {
      "measurement", "hamiltonian"
    };
    double termSeconds = 0.0;
    long long gateEvaluated = 0;
    long long gateFailures = 0;
    int j;
    fprintf(fp,"BF_FSZ_C2_DETAIL_PROFILE_BEGIN version=1 rank=0 env=MVMC_BF_FSZ_C2_DETAIL_PROFILE\n");
    for(i=0;i<NBFFSZC2DetailSource;i++) {
      for(j=0;j<NBFFSZC2DetailClass;j++) {
        fprintf(fp,"source=%s class=%d calls=%lld\n",sourceLabels[i],j+1,
                BFFSZC2DetailClassCall[i][j]);
      }
      for(j=0;j<NBFFSZC2DetailOutcome;j++) {
        fprintf(fp,"source=%s outcome=%s calls=%lld\n",sourceLabels[i],
                outcomeLabels[j],BFFSZC2DetailOutcome[i][j]);
      }
      for(j=0;j<NBFFSZC2DetailPath;j++) {
        fprintf(fp,"source=%s path=%s calls=%lld\n",sourceLabels[i],
                pathLabels[j],BFFSZC2DetailPath[i][j]);
      }
      fprintf(fp,"source=%s changed calls=%lld sum=%lld max=%lld\n",sourceLabels[i],
              BFFSZC2DetailEvaluatedCalls[i],BFFSZC2DetailChangedSum[i],
              BFFSZC2DetailChangedMax[i]);
      fprintf(fp,"source=%s affected calls=%lld sum=%lld max=%lld at_or_above_k_full=%lld\n",
              sourceLabels[i],BFFSZC2DetailEvaluatedCalls[i],
              BFFSZC2DetailAffectedSum[i],BFFSZC2DetailAffectedMax[i],
              BFFSZC2DetailAffectedAtOrAboveKFull[i]);
      fprintf(fp,"source=%s ordered_descriptors=%lld\n",sourceLabels[i],
              BFFSZC2DetailOrderedDescriptorTotal[i]);
      for(j=0;j<NBFFSZProfileHist;j++) {
        fprintf(fp,"source=%s changed_hist bin=%d calls=%lld\n",sourceLabels[i],j,
                BFFSZC2DetailChangedHist[i][j]);
        fprintf(fp,"source=%s affected_hist bin=%d calls=%lld\n",sourceLabels[i],j,
                BFFSZC2DetailAffectedHist[i][j]);
      }
      for(j=0;j<NBFFSZC2DetailComponent;j++) {
        fprintf(fp,"source=%s component=%s seconds_rank0=%.9f\n",sourceLabels[i],
                componentLabels[j],BFFSZC2DetailComponentSeconds[i][j]);
      }
      gateEvaluated += BFFSZC2DetailEvaluatedCalls[i];
      gateFailures += BFFSZC2DetailAffectedAtOrAboveKFull[i];
    }
    for(i=0;i<NBFFSZC2DetailTerm;i++) {
      fprintf(fp,"hamiltonian_term=%s seconds_rank0=%.9f\n",termLabels[i],
              BFFSZC2DetailTermSeconds[i]);
      termSeconds += BFFSZC2DetailTermSeconds[i];
    }
    fprintf(fp,"sz_seconds_rank0=%.9f\n",BFFSZC2DetailSzSecondsRank0);
    fprintf(fp,"timer41_seconds_rank0=%.9f hamiltonian_term_seconds_rank0=%.9f sz_seconds_rank0=%.9f wrapper_overhead_seconds_rank0=%.9f\n",
            Timer[41],termSeconds,BFFSZC2DetailSzSecondsRank0,
            Timer[41]-termSeconds-BFFSZC2DetailSzSecondsRank0);
    fprintf(fp,"affected_gate k_full=%d evaluated_calls=%lld at_or_above_k_full=%lld pass=%d\n",
            BFFSZPfUpdateKFull,gateEvaluated,gateFailures,(gateFailures == 0));
    if(BFFSZC2ReuseCensusEnabled) {
      for(i=0;i<NBFFSZC2ReuseScope;i++) {
        const double reuseRatio = (BFFSZC2ReuseCensusUniqueExactOrderedMoves[i] > 0)
            ? (double)BFFSZC2ReuseCensusTrueCalls[i]
                /(double)BFFSZC2ReuseCensusUniqueExactOrderedMoves[i] : 0.0;
        fprintf(fp,"reuse_census scope=%s invocations=%lld true_calls=%lld unique_exact_ordered_moves=%lld duplicate_true_calls=%lld overflow_calls=%lld reuse_ratio=%.9f\n",
                reuseScopeLabels[i],BFFSZC2ReuseCensusInvocations[i],
                BFFSZC2ReuseCensusTrueCalls[i],
                BFFSZC2ReuseCensusUniqueExactOrderedMoves[i],
                BFFSZC2ReuseCensusDuplicateTrueCalls[i],
                BFFSZC2ReuseCensusOverflowCalls[i],reuseRatio);
      }
    }
    fprintf(fp,"BF_FSZ_C2_DETAIL_PROFILE_END\n");
  }
}

void OutputTimerParaOpt() {
  char fileName[D_FileNameMax];
  FILE *fp;
  sprintf(fileName, "%s_CalcTimer.dat", CDataFileHead); 
  fp = fopen(fileName, "w");

  fprintf(fp,"All                         [0] %12.5lf\n",Timer[0]);
  fprintf(fp,"Initialization              [1] %12.5lf\n",Timer[1]);
  fprintf(fp,"  read options             [10] %12.5lf\n",Timer[10]);
  fprintf(fp,"  ReadDefFile              [11] %12.5lf\n",Timer[11]);
  fprintf(fp,"  SetMemory                [12] %12.5lf\n",Timer[12]);
  fprintf(fp,"  InitParameter            [13] %12.5lf\n",Timer[13]);
  fprintf(fp,"VMCParaOpt                  [2] %12.5lf\n",Timer[2]);
  fprintf(fp,"  VMCMakeSample             [3] %12.5lf\n",Timer[3]);
  fprintf(fp,"    makeInitialSample      [30] %12.5lf\n",Timer[30]);
  fprintf(fp,"    make candidate         [31] %12.5lf\n",Timer[31]);
  fprintf(fp,"    hopping update         [32] %12.5lf\n",Timer[32]);
  fprintf(fp,"      UpdateProjCnt        [60] %12.5lf\n",Timer[60]);
  fprintf(fp,"      CalculateNewPfM2     [61] %12.5lf\n",Timer[61]);
  fprintf(fp,"      CalculateLogIP       [62] %12.5lf\n",Timer[62]);
  fprintf(fp,"      UpdateMAll           [63] %12.5lf\n",Timer[63]);
  fprintf(fp,"      UpdateSlaterElmBF    [64] %12.5lf\n",Timer[64]);
  fprintf(fp,"      CommitSlaterElmBF    [94] %12.5lf\n",Timer[94]);
  fprintf(fp,"    exchange update        [33] %12.5lf\n",Timer[33]);
  fprintf(fp,"      UpdateProjCnt        [65] %12.5lf\n",Timer[65]);
  fprintf(fp,"      CalculateNewPfMTwo2  [66] %12.5lf\n",Timer[66]);
  fprintf(fp,"      CalculateLogIP       [67] %12.5lf\n",Timer[67]);
  fprintf(fp,"      UpdateMAllTwo        [68] %12.5lf\n",Timer[68]);
  fprintf(fp,"    lspinflip update       [36] %12.5lf\n",Timer[36]);
  fprintf(fp,"      UpdateProjCnt       [600] %12.5lf\n",Timer[600]);
  fprintf(fp,"      CalculateNewPfMTwo2 [601] %12.5lf\n",Timer[601]);
  fprintf(fp,"      CalculateLogIP      [602] %12.5lf\n",Timer[602]);
  fprintf(fp,"      UpdateMAllTwo       [603] %12.5lf\n",Timer[603]);
  fprintf(fp,"    recal PfM and InvM     [34] %12.5lf\n",Timer[34]);
  fprintf(fp,"    save electron config   [35] %12.5lf\n",Timer[35]);
  fprintf(fp,"  VMCMainCal                [4] %12.5lf\n",Timer[4]);
  fprintf(fp,"    CalculateMAll          [40] %12.5lf\n",Timer[40]);
  fprintf(fp,"    LocEnergyCal           [41] %12.5lf\n",Timer[41]);
  fprintf(fp,"      CalHamiltonian0      [70] %12.5lf\n",Timer[70]);
  fprintf(fp,"      CalHamiltonian1      [71] %12.5lf\n",Timer[71]);
  fprintf(fp,"      CalHamiltonian2      [72] %12.5lf\n",Timer[72]);
  fprintf(fp,"    ReturnSlaterElmDiff    [42] %12.5lf\n",Timer[42]);
  fprintf(fp,"    calculate OO and HO    [43] %12.5lf\n",Timer[43]);
  fprintf(fp,"    multiply store OO      [45] %12.5lf\n",Timer[45]);
  fprintf(fp,"  BF detail timers\n");
  fprintf(fp,"    BF Green MakeProjBFCnt [81] %12.5lf\n",Timer[81]);
  fprintf(fp,"    BF Green UpdateSlt     [82] %12.5lf\n",Timer[82]);
  fprintf(fp,"    BF Green CalcNewPfM    [83] %12.5lf\n",Timer[83]);
  fprintf(fp,"    BF GreenFunc1 real     [84] %12.5lf\n",Timer[84]);
  fprintf(fp,"    BF GreenFunc2 real     [85] %12.5lf\n",Timer[85]);
  fprintf(fp,"    BF Green setup/proj    [86] %12.5lf\n",Timer[86]);
  fprintf(fp,"    BF Green CalculateIP   [87] %12.5lf\n",Timer[87]);
  fprintf(fp,"    BF Green restore cfg   [88] %12.5lf\n",Timer[88]);
  fprintf(fp,"    BF Green StoreSlater   [89] %12.5lf\n",Timer[89]);
  fprintf(fp,"    BF Green merge lists   [90] %12.5lf\n",Timer[90]);
  fprintf(fp,"    BF Update collect rows [91] %12.5lf\n",Timer[91]);
  fprintf(fp,"    BF Update recompute    [92] %12.5lf\n",Timer[92]);
  fprintf(fp,"    BF Update copy rows    [93] %12.5lf\n",Timer[93]);
  OutputBFProfileCounters(fp);
  fprintf(fp,"  StochasticOpt             [5] %12.5lf\n",Timer[5]);
  fprintf(fp,"    preprocess             [50] %12.5lf\n",Timer[50]);
  fprintf(fp,"    stcOptMain             [51] %12.5lf\n",Timer[51]);
  fprintf(fp,"      initBLACS            [55] %12.5lf\n",Timer[55]);
  fprintf(fp,"      calculate S and g    [56] %12.5lf\n",Timer[56]);
  fprintf(fp,"      DPOSV                [57] %12.5lf\n",Timer[57]);
  fprintf(fp,"      gatherParaChange     [58] %12.5lf\n",Timer[58]);
  fprintf(fp,"    postprocess            [52] %12.5lf\n",Timer[52]);
  fprintf(fp,"  UpdateSlaterElm          [20] %12.5lf\n",Timer[20]);
  fprintf(fp,"  WeightAverage            [21] %12.5lf\n",Timer[21]);
  fprintf(fp,"  outputData               [22] %12.5lf\n",Timer[22]);
  fprintf(fp,"  SyncModifiedParameter    [23] %12.5lf\n",Timer[23]);
  fprintf(fp,"  cal                      [24] %12.5lf\n",Timer[24]);
  fprintf(fp,"  SR                       [25] %12.5lf\n",Timer[25]);
  fprintf(fp,"  MAll                     [69] %12.5lf\n",Timer[69]);

  fclose(fp);
}

void OutputTimerPhysCal() {
  char fileName[D_FileNameMax];
  FILE *fp;
  sprintf(fileName, "%s_CalcTimer.dat", CDataFileHead); 
  fp = fopen(fileName, "w");

  fprintf(fp,"All                         [0] %12.5lf\n",Timer[0]);
  fprintf(fp,"Initialization              [1] %12.5lf\n",Timer[1]);
  fprintf(fp,"  read options             [10] %12.5lf\n",Timer[10]);
  fprintf(fp,"  ReadDefFile              [11] %12.5lf\n",Timer[11]);
  fprintf(fp,"  SetMemory                [12] %12.5lf\n",Timer[12]);
  fprintf(fp,"  InitParameter            [13] %12.5lf\n",Timer[13]);
  fprintf(fp,"VMCPhysCal                  [2] %12.5lf\n",Timer[2]);
  fprintf(fp,"  VMCMakeSample             [3] %12.5lf\n",Timer[3]);
  fprintf(fp,"    makeInitialSample      [30] %12.5lf\n",Timer[30]);
  fprintf(fp,"    make candidate         [31] %12.5lf\n",Timer[31]);
  fprintf(fp,"    hopping update         [32] %12.5lf\n",Timer[32]);
  fprintf(fp,"      UpdateProjCnt        [60] %12.5lf\n",Timer[60]);
  fprintf(fp,"      CalculateNewPfM2     [61] %12.5lf\n",Timer[61]);
  fprintf(fp,"      CalculateLogIP       [62] %12.5lf\n",Timer[62]);
  fprintf(fp,"      UpdateMAll           [63] %12.5lf\n",Timer[63]);
  fprintf(fp,"      UpdateSlaterElmBF    [64] %12.5lf\n",Timer[64]);
  fprintf(fp,"      CommitSlaterElmBF    [94] %12.5lf\n",Timer[94]);
  fprintf(fp,"    exchange update        [33] %12.5lf\n",Timer[33]);
  fprintf(fp,"      UpdateProjCnt        [65] %12.5lf\n",Timer[65]);
  fprintf(fp,"      CalculateNewPfMTwo2  [66] %12.5lf\n",Timer[66]);
  fprintf(fp,"      CalculateLogIP       [67] %12.5lf\n",Timer[67]);
  fprintf(fp,"      UpdateMAllTwo        [68] %12.5lf\n",Timer[68]);
  fprintf(fp,"    lspinflip update       [36] %12.5lf\n",Timer[36]);
  fprintf(fp,"      UpdateProjCnt       [600] %12.5lf\n",Timer[600]);
  fprintf(fp,"      CalculateNewPfMTwo2 [601] %12.5lf\n",Timer[601]);
  fprintf(fp,"      CalculateLogIP      [602] %12.5lf\n",Timer[602]);
  fprintf(fp,"      UpdateMAllTwo       [603] %12.5lf\n",Timer[603]);
  fprintf(fp,"    recal PfM and InvM     [34] %12.5lf\n",Timer[34]);
  fprintf(fp,"    save electron config   [35] %12.5lf\n",Timer[35]);
  fprintf(fp,"  VMCMainCal                [4] %12.5lf\n",Timer[4]);
  fprintf(fp,"    CalculateMAll          [40] %12.5lf\n",Timer[40]);
  fprintf(fp,"    LocEnergyCal           [41] %12.5lf\n",Timer[41]);
  fprintf(fp,"      CalHamiltonian0      [70] %12.5lf\n",Timer[70]);
  fprintf(fp,"      CalHamiltonian1      [71] %12.5lf\n",Timer[71]);
  fprintf(fp,"      CalHamiltonian2      [72] %12.5lf\n",Timer[72]);
  fprintf(fp,"    CalculateGreenFunc     [42] %12.5lf\n",Timer[42]);
  fprintf(fp,"      GreenFunc1           [50] %12.5lf\n",Timer[50]);
  fprintf(fp,"      GreenFunc2           [51] %12.5lf\n",Timer[51]);
  fprintf(fp,"      addPhysCA            [52] %12.5lf\n",Timer[52]);
  fprintf(fp,"      addPhysCACA          [53] %12.5lf\n",Timer[53]);
  fprintf(fp,"      addPhysTwist         [54] %12.5lf\n",Timer[54]);
  fprintf(fp,"  BF detail timers\n");
  fprintf(fp,"    BF Green MakeProjBFCnt [81] %12.5lf\n",Timer[81]);
  fprintf(fp,"    BF Green UpdateSlt     [82] %12.5lf\n",Timer[82]);
  fprintf(fp,"    BF Green CalcNewPfM    [83] %12.5lf\n",Timer[83]);
  fprintf(fp,"    BF GreenFunc1 real     [84] %12.5lf\n",Timer[84]);
  fprintf(fp,"    BF GreenFunc2 real     [85] %12.5lf\n",Timer[85]);
  fprintf(fp,"    BF Green setup/proj    [86] %12.5lf\n",Timer[86]);
  fprintf(fp,"    BF Green CalculateIP   [87] %12.5lf\n",Timer[87]);
  fprintf(fp,"    BF Green restore cfg   [88] %12.5lf\n",Timer[88]);
  fprintf(fp,"    BF Green StoreSlater   [89] %12.5lf\n",Timer[89]);
  fprintf(fp,"    BF Green merge lists   [90] %12.5lf\n",Timer[90]);
  fprintf(fp,"    BF Update collect rows [91] %12.5lf\n",Timer[91]);
  fprintf(fp,"    BF Update recompute    [92] %12.5lf\n",Timer[92]);
  fprintf(fp,"    BF Update copy rows    [93] %12.5lf\n",Timer[93]);
  OutputBFProfileCounters(fp);
  fprintf(fp,"    Lanczos1               [43] %12.5lf\n",Timer[43]);
  fprintf(fp,"    Lanczos2               [44] %12.5lf\n",Timer[44]);
  fprintf(fp,"  UpdateSlaterElm          [20] %12.5lf\n",Timer[20]);
  fprintf(fp,"  WeightAverage            [21] %12.5lf\n",Timer[21]);
  fprintf(fp,"  outputData               [22] %12.5lf\n",Timer[22]);

  fclose(fp);
}

#endif
