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
 * local Green Functions
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#pragma once
#include "backflow.h"
#include "backflow_nbody.h"
#include "locgrn.h"
#include "projection.h"
#include "rbm.h"
#include "pfupdate.h"
#include "pfupdate_two_fcmp.h"
#include "qp.h"
#include "sector_projection.h"

#include "locgrn_real.h"
#include "calham.h"
#include "calham_real.h"
#include "calgrn.h"
#include "workspace.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int BFCanonicalGreenInjectFailure(const char *target) {
  const char *value = getenv("MVMC_BF_CANONICAL_GREEN_INJECT_FAILURE");
  return value != NULL && target != NULL && strcmp(value, target) == 0;
}

static int AppendBFHopIndex(int *dst, int *count, const int stride, const int value) {
  int i;
  for (i = 0; i < *count; i++) {
    if (dst[i] == value) return 0;
  }
  if (*count >= stride) return 1;
  dst[*count] = value;
  (*count)++;
  return 0;
}

static int MergeBFHopLists(const int *left, const int *leftCount,
                           const int *right, const int *rightCount,
                           int *merged, int *mergedCount) {
  int qpidx, i;
  for (qpidx = 0; qpidx < NQPFull; qpidx++) {
    int *dst = merged + qpidx * Nsize;
    mergedCount[qpidx] = 0;
    for (i = 0; i < leftCount[qpidx]; i++) {
      if (AppendBFHopIndex(dst, &mergedCount[qpidx], Nsize, left[qpidx * Nsize + i]) != 0) return 1;
    }
    for (i = 0; i < rightCount[qpidx]; i++) {
      if (AppendBFHopIndex(dst, &mergedCount[qpidx], Nsize, right[qpidx * Nsize + i]) != 0) return 1;
    }
  }
  return 0;
}

/* Calculate 1-body Green function <CisAjs> */
/* buffer size = NQPFull */
double complex GreenFunc1BF(const int ri, const int rj, const int s, const double complex ip, double complex *bufM,
                    int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,
                    int *projCntNew, const int *eleProjBFCnt,int *projBFCntNew,
                    double complex *buffer, int *greenStatus) {
  double complex z;
  int mj,msj,rsi,rsj;
  //double complex bufM[NQPFull*Nsize*Nsize];
  double complex* pfMNew = buffer; /* NQPFull */
  int msaTmp[NQPFull*Nsize],icount[NQPFull];

  if(greenStatus == NULL) return NAN + NAN*I;
  *greenStatus = BF_PF_OK;
  if(ri==rj) return eleNum[ri+s*Nsite];
  if(eleNum[ri+s*Nsite]==1 || eleNum[rj+s*Nsite]==0) return 0.0;

  /* store SlaterElmBF before hopping */
  mj = eleCfg[rj+s*Nsite];
  msj = mj + s*Ne;
  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;

  /* hopping */
  eleCfg[rsj] = -1;
  eleCfg[rsi] = mj;
  eleIdx[msj] = ri;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt(rj, ri, s, projCntNew, eleProjCnt, eleNum);
  z = ProjRatio(projCntNew,eleProjCnt);

  /* calculate Pfaffian */
  //printf("1");
  StartTimer(81);
  MakeProjBFCnt(projBFCntNew, eleNum);
  StopTimer(81);
  if(BFUseCanonicalNonFszPath()) {
    const int rebuildStatus = BFCanonicalGreenInjectFailure("complex1")
        ? BF_PF_INVALID_ARGUMENT
        : CalculateBFCanonicalPf_fcmp(
              eleIdx, eleNum, projBFCntNew, 0, NQPFull, pfMNew);
    if(rebuildStatus == 0) {
      z *= CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);
    } else {
      *greenStatus = rebuildStatus;
    }

    eleCfg[rsj] = mj;
    eleCfg[rsi] = -1;
    eleIdx[msj] = rj;
    eleNum[rsj] = 1;
    eleNum[rsi] = 0;
    return rebuildStatus == 0 ? conj(z/ip) : 0.0+0.0*I;
  }
  StartTimer(82);
  UpdateSlaterElmBFGrn(mj, rj, ri, s, eleCfg, eleNum, projBFCntNew, msaTmp, icount, bufM);
  StopTimer(82);
  StartTimer(83);
  CalculateNewPfMBF(icount, msaTmp, pfMNew, eleIdx, 0, NQPFull, bufM);
  StopTimer(83);
  z *= CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleCfg[rsj] = mj;
  eleCfg[rsi] = -1;
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  StoreSlaterElmBF_fcmp(bufM);

  return conj(z/ip);
}

/* Calculate 2-body Green function with BackFlow. */
double complex GreenFunc2BF(const int ri, const int rj, const int rk, const int rl,
                    const int s, const int t, const double complex ip, double complex *bufM,
                    int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,
                    int *projCntNew, const int *eleProjBFCnt,int *projBFCntNew,
                    double complex *buffer, int *greenStatus) {
  double complex z;
  int mj,msj,ml,mtl;
  int rsi,rsj,rtk,rtl;
  double complex *pfMNew = buffer; /* [NQPFull] */
  int msaTmp0[NQPFull*Nsize], msaTmp1[NQPFull*Nsize], msaTmp[NQPFull*Nsize];
  int icount0[NQPFull], icount1[NQPFull], icount[NQPFull];

  if(greenStatus == NULL) return NAN + NAN*I;
  *greenStatus = BF_PF_OK;
  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;
  rtk = rk + t*Nsite;
  rtl = rl + t*Nsite;

  if(s==t) {
    if(rk==rl) {
      if(eleNum[rtk]==0) return 0.0;
      else return GreenFunc1BF(ri,rj,s,ip,bufM,eleIdx,eleCfg,eleNum,
                               eleProjCnt,projCntNew,eleProjBFCnt,projBFCntNew,
                               buffer,greenStatus);
    }else if(rj==rl) {
      return 0.0;
    }else if(ri==rl) {
      if(eleNum[rsi]==0) return 0.0;
      else if(rj==rk) return 1.0-eleNum[rsj];
      else return -GreenFunc1BF(rk,rj,s,ip,bufM,eleIdx,eleCfg,eleNum,
                                eleProjCnt,projCntNew,eleProjBFCnt,projBFCntNew,
                                buffer,greenStatus);
    }else if(rj==rk) {
      if(eleNum[rsj]==1) return 0.0;
      else if(ri==rl) return eleNum[rsi];
      else return GreenFunc1BF(ri,rl,s,ip,bufM,eleIdx,eleCfg,eleNum,
                               eleProjCnt,projCntNew,eleProjBFCnt,projBFCntNew,
                               buffer,greenStatus);
    }else if(ri==rk) {
      return 0.0;
    }else if(ri==rj) {
      if(eleNum[rsi]==0) return 0.0;
      else return GreenFunc1BF(rk,rl,s,ip,bufM,eleIdx,eleCfg,eleNum,
                               eleProjCnt,projCntNew,eleProjBFCnt,projBFCntNew,
                               buffer,greenStatus);
    }
  }else{
    if(rk==rl) {
      if(eleNum[rtk]==0) return 0.0;
      else if(ri==rj) return eleNum[rsi];
      else return GreenFunc1BF(ri,rj,s,ip,bufM,eleIdx,eleCfg,eleNum,
                               eleProjCnt,projCntNew,eleProjBFCnt,projBFCntNew,
                               buffer,greenStatus);
    }else if(ri==rj) {
      if(eleNum[rsi]==0) return 0.0;
      else return GreenFunc1BF(rk,rl,t,ip,bufM,eleIdx,eleCfg,eleNum,
                               eleProjCnt,projCntNew,eleProjBFCnt,projBFCntNew,
                               buffer,greenStatus);
    }
  }

  if(eleNum[rsi]==1 || eleNum[rsj]==0 || eleNum[rtk]==1 || eleNum[rtl]==0) return 0.0;

  mj = eleCfg[rsj];
  ml = eleCfg[rtl];
  msj = mj + s*Ne;
  mtl = ml + t*Ne;

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

  if(!IsSectorStateAllowed(eleNum)) {
    eleCfg[rtl] = ml;
    eleCfg[rtk] = -1;
    eleIdx[mtl] = rl;
    eleNum[rtl] = 1;
    eleNum[rtk] = 0;
    eleCfg[rsj] = mj;
    eleCfg[rsi] = -1;
    eleIdx[msj] = rj;
    eleNum[rsj] = 1;
    eleNum[rsi] = 0;
    return 0.0;
  }

  z = ProjRatio(projCntNew,eleProjCnt);

  StartTimer(81);
  MakeProjBFCnt(projBFCntNew, eleNum);
  StopTimer(81);
  if(BFUseCanonicalNonFszPath()) {
    const int rebuildStatus = BFCanonicalGreenInjectFailure("complex2")
        ? BF_PF_INVALID_ARGUMENT
        : CalculateBFCanonicalPf_fcmp(
              eleIdx, eleNum, projBFCntNew, 0, NQPFull, pfMNew);
    if(rebuildStatus == 0) {
      z *= CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);
    } else {
      *greenStatus = rebuildStatus;
    }

    eleCfg[rtl] = ml;
    eleCfg[rtk] = -1;
    eleIdx[mtl] = rl;
    eleNum[rtl] = 1;
    eleNum[rtk] = 0;
    eleCfg[rsj] = mj;
    eleCfg[rsi] = -1;
    eleIdx[msj] = rj;
    eleNum[rsj] = 1;
    eleNum[rsi] = 0;
    return rebuildStatus == 0 ? conj(z/ip) : 0.0+0.0*I;
  }
  StartTimer(82);
  UpdateSlaterElmBFGrn(ml, rl, rk, t, eleCfg, eleNum, projBFCntNew, msaTmp0, icount0, bufM);
  UpdateSlaterElmBFGrn(mj, rj, ri, s, eleCfg, eleNum, projBFCntNew, msaTmp1, icount1, bufM);
  StopTimer(82);

  if (MergeBFHopLists(msaTmp0, icount0, msaTmp1, icount1, msaTmp, icount) != 0) {
    eleCfg[rtl] = ml;
    eleCfg[rtk] = -1;
    eleIdx[mtl] = rl;
    eleNum[rtl] = 1;
    eleNum[rtk] = 0;
    eleCfg[rsj] = mj;
    eleCfg[rsi] = -1;
    eleIdx[msj] = rj;
    eleNum[rsj] = 1;
    eleNum[rsi] = 0;
    StoreSlaterElmBF_fcmp(bufM);
    *greenStatus = BF_PF_INVALID_ARGUMENT;
    return 0.0;
  }

  StartTimer(83);
  CalculateNewPfMBFWithStride(icount, msaTmp, Nsize, pfMNew, eleIdx, 0, NQPFull, bufM);
  StopTimer(83);
  z *= CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  eleCfg[rtl] = ml;
  eleCfg[rtk] = -1;
  eleIdx[mtl] = rl;
  eleNum[rtl] = 1;
  eleNum[rtk] = 0;
  eleCfg[rsj] = mj;
  eleCfg[rsi] = -1;
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  StoreSlaterElmBF_fcmp(bufM);

  return conj(z/ip);
}
static int FindBFHopIndex_real(const int *list, const int count, const int value) {
  int i;
  for (i = 0; i < count; i++) {
    if (list[i] == value) return i;
  }
  return -1;
}

static void CopyBFVecRow_real(double *dstVec, const int dstRow,
                              const double *srcVec, const int srcRow) {
  int i;
  double *dst = dstVec + dstRow*Nsize;
  const double *src = srcVec + srcRow*Nsize;
  if (dst == src) return;
  for (i = 0; i < Nsize; i++) {
    dst[i] = src[i];
  }
}

static int MergeBFHopListsVec_real(const int *left, const int *leftCount, const double *leftVec,
                                   const int *right, const int *rightCount, const double *rightVec,
                                   int *merged, int *mergedCount, double *mergedVec) {
  int qpidx, i, outIdx;
  for (qpidx = 0; qpidx < NQPFull; qpidx++) {
    int *dst = merged + qpidx * Nsize;
    const int *leftQp = left + qpidx * Nsize;
    const int *rightQp = right + qpidx * Nsize;
    const double *leftVecQp = leftVec + qpidx * Nsize * Nsize;
    const double *rightVecQp = rightVec + qpidx * Nsize * Nsize;
    double *mergedVecQp = mergedVec + qpidx * Nsize * Nsize;
    mergedCount[qpidx] = 0;

    for (i = 0; i < leftCount[qpidx]; i++) {
      if (FindBFHopIndex_real(dst, mergedCount[qpidx], leftQp[i]) >= 0) continue;
      if (mergedCount[qpidx] >= Nsize) return 1;
      outIdx = mergedCount[qpidx];
      dst[outIdx] = leftQp[i];
      CopyBFVecRow_real(mergedVecQp, outIdx, leftVecQp, i);
      mergedCount[qpidx]++;
    }
    for (i = 0; i < rightCount[qpidx]; i++) {
      if (FindBFHopIndex_real(dst, mergedCount[qpidx], rightQp[i]) >= 0) continue;
      if (mergedCount[qpidx] >= Nsize) return 1;
      outIdx = mergedCount[qpidx];
      dst[outIdx] = rightQp[i];
      CopyBFVecRow_real(mergedVecQp, outIdx, rightVecQp, i);
      mergedCount[qpidx]++;
    }
  }
  return 0;
}


/* Calculate 1-body Green function <CisAjs> */
/* buffer size = NQPFull */
int GreenFunc1BF_real_prepare(const int ri, const int rj, const int s, double *greenValue,
                    double *projRatio, double *vecM, int *msaTmp, int *icount,
                    int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,
                    int *projCntNew, const int *eleProjBFCnt,int *projBFCntNew) {
  int mj,msj,rsi,rsj;

  StartTimer(84);
  if(ri==rj) {
    *greenValue = eleNum[ri+s*Nsite];
    StopTimer(84);
    return 0;
  }
  if(eleNum[ri+s*Nsite]==1 || eleNum[rj+s*Nsite]==0) {
    *greenValue = 0.0;
    StopTimer(84);
    return 0;
  }

  mj = eleCfg[rj+s*Nsite];
  msj = mj + s*Ne;
  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;

  StartTimer(86);
  eleCfg[rsj] = -1;
  eleCfg[rsi] = mj;
  eleIdx[msj] = ri;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt(rj, ri, s, projCntNew, eleProjCnt, eleNum);
  *projRatio = ProjRatio(projCntNew,eleProjCnt);
  StopTimer(86);

  StartTimer(81);
  MakeProjBFCnt(projBFCntNew, eleNum);
  StopTimer(81);
  StartTimer(82);
  UpdateSlaterElmBFGrnVec_real(mj, rj, ri, s, eleIdx, eleCfg, eleNum,
                               eleProjBFCnt, projBFCntNew, msaTmp, icount, vecM);
  StopTimer(82);

  StartTimer(88);
  eleCfg[rsj] = mj;
  eleCfg[rsi] = -1;
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;
  StopTimer(88);

  StopTimer(84);
  return 1;
}

void GreenFunc1BF_real_finish_batch(const int batchSize, const double ip,
                    const double *projRatio, const int *icount, const int *msaTmp,
                    const double *vecM, double *greenValue, double *pfMNew,
                    double *vecStack, double *wStack) {
  int batchIdx;

  StartTimer(84);
  StartTimer(83);
  CalculateNewPfMBFVecBatched_real(batchSize, icount, msaTmp, pfMNew, 0, NQPFull,
                                   vecM, vecStack, wStack);
  StopTimer(83);
  StartTimer(87);
  for(batchIdx=0;batchIdx<batchSize;batchIdx++) {
    greenValue[batchIdx] = projRatio[batchIdx]
      * CalculateIP_real(pfMNew + batchIdx*NQPFull, 0, NQPFull, MPI_COMM_SELF) / ip;
  }
  StopTimer(87);
  StopTimer(84);

  return;
}

double GreenFunc1BF_real(const int ri, const int rj, const int s, const double ip, double *bufM,
                    int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,
                    int *projCntNew, const int *eleProjBFCnt,int *projBFCntNew,
                    double *buffer, int *greenStatus) {
  double z;
  int mj,msj,rsi,rsj;
  //double complex bufM[NQPFull*Nsize*Nsize];
  double *pfMNew_real = buffer; /* NQPFull */
  int msaTmp[NQPFull*Nsize],icount[NQPFull];

#define RETURN_GREENFUNC1BF_REAL(value) do { z = (value); StopTimer(84); return z; } while (0)

  if(greenStatus == NULL) return NAN;
  *greenStatus = BF_PF_OK;
  StartTimer(84);
  if(ri==rj) RETURN_GREENFUNC1BF_REAL(eleNum[ri+s*Nsite]);
  if(eleNum[ri+s*Nsite]==1 || eleNum[rj+s*Nsite]==0) RETURN_GREENFUNC1BF_REAL(0.0);

  /* store SlaterElmBF before hopping */
  mj = eleCfg[rj+s*Nsite];
  msj = mj + s*Ne;
  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;

  StartTimer(86);
  /* hopping */
  eleCfg[rsj] = -1;
  eleCfg[rsi] = mj;
  eleIdx[msj] = ri;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt(rj, ri, s, projCntNew, eleProjCnt, eleNum);
  z = ProjRatio(projCntNew,eleProjCnt);
  StopTimer(86);

  /* calculate Pfaffian */
  //printf("1");
  StartTimer(81);
  MakeProjBFCnt(projBFCntNew, eleNum);
  StopTimer(81);
  if(BFUseCanonicalNonFszPath()) {
    const int rebuildStatus = BFCanonicalGreenInjectFailure("real1")
        ? BF_PF_INVALID_ARGUMENT
        : CalculateBFCanonicalPf_real(
              eleIdx, eleNum, projBFCntNew, 0, NQPFull, pfMNew_real);
    if(rebuildStatus == 0) {
      z *= CalculateIP_real(pfMNew_real, 0, NQPFull, MPI_COMM_SELF);
    } else {
      *greenStatus = rebuildStatus;
    }

    eleCfg[rsj] = mj;
    eleCfg[rsi] = -1;
    eleIdx[msj] = rj;
    eleNum[rsj] = 1;
    eleNum[rsi] = 0;
    RETURN_GREENFUNC1BF_REAL(rebuildStatus == 0 ? z/ip : 0.0);
  }
  StartTimer(82);
  UpdateSlaterElmBFGrnVec_real(mj, rj, ri, s, eleIdx, eleCfg, eleNum,
                               eleProjBFCnt, projBFCntNew, msaTmp, icount, bufM);
  StopTimer(82);
  StartTimer(83);
  CalculateNewPfMBFVec_real(icount, msaTmp, pfMNew_real, 0, NQPFull, bufM);
  StopTimer(83);
  StartTimer(87);
  z *= CalculateIP_real(pfMNew_real, 0, NQPFull, MPI_COMM_SELF);
  StopTimer(87);
  //printf("1");
  //MakeSlaterElmBF(eleNum,projBFCntNew);
  //UpdateSlaterElmBF2(eleNum,projBFCntNew);
  //printf("1");
  //CalculateMAllBF_NoThread(eleIdx,0,NQPFull,bufM);
  //CalculateMAllBF_NoThread(eleIdx,0,NQPFull,SlaterElmBF);
  //printf("1\n");
  //z *= CalculateIP(PfM,0,NQPFull,MPI_COMM_SELF);

  StartTimer(88);
  /* revert hopping */
  eleCfg[rsj] = mj;
  eleCfg[rsi] = -1;
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;
  StopTimer(88);

  //UpdateSlaterElmBFTmp5(mj, ri, rj, s, eleCfg, eleNum, eleProjBFCnt, msaTmp, icount, bufM);
  //MakeProjBFCnt(projBFCntNew, eleNum);
  //MakeSlaterElmBF(eleNum,projBFCntNew);
  //CalculateMAllBF_NoThread(eleIdx,0,NQPFull,bufM);
  //CalculateMAllBF_NoThread(eleIdx,0,NQPFull,SlaterElmBF);

  RETURN_GREENFUNC1BF_REAL(z/ip);
#undef RETURN_GREENFUNC1BF_REAL
}

/* Calculate 2-body Green function with BackFlow. */
/* Two-body Green function with BackFlow, real variational parameters.
 * vecTmp0 / vecTmp1: caller-owned packed row-vector workspaces of
 * NQPFull*Nsize*Nsize doubles each (no allocation inside).
 * Returns 0 and stores the Green function in *value, or a nonzero status when
 * an internal invariant fails (hop-list merge overflow).  A physical zero is
 * reported as status 0 with *value == 0. */
int GreenFunc2BF_real_ws(const int ri, const int rj, const int rk, const int rl,
                    const int s, const int t, const double ip,
                    double *vecTmp0, double *vecTmp1,
                    int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,
                    int *projCntNew, const int *eleProjBFCnt,int *projBFCntNew, double *buffer,
                    double *value) {
  double z;
  double nestedValue;
  int mj,msj,ml,mtl;
  int rsi,rsj,rtk,rtl;
  int nestedStatus = BF_PF_OK;
  double *pfMNew_real = buffer; /* [NQPFull] */
  int msaTmp0[NQPFull*Nsize], msaTmp1[NQPFull*Nsize], msaTmp[NQPFull*Nsize];
  int icount0[NQPFull], icount1[NQPFull], icount[NQPFull];

#define RETURN_GREENFUNC2BF_REAL(v) do { *value = (v); StopTimer(85); return 0; } while (0)
#define RETURN_GREENFUNC1BF_REAL_STATUS(v) do { \
    nestedValue = (v); \
    if(nestedStatus != BF_PF_OK) { \
      *value = 0.0; \
      StopTimer(85); \
      return nestedStatus; \
    } \
    RETURN_GREENFUNC2BF_REAL(nestedValue); \
  } while (0)

  if(value == NULL || vecTmp0 == NULL || vecTmp1 == NULL || buffer == NULL) {
    if(value != NULL) *value = 0.0;
    return BF_PF_INVALID_ARGUMENT;
  }
  StartTimer(85);

  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;
  rtk = rk + t*Nsite;
  rtl = rl + t*Nsite;

  if(s==t) {
    if(rk==rl) {
      if(eleNum[rtk]==0) RETURN_GREENFUNC2BF_REAL(0.0);
      else RETURN_GREENFUNC1BF_REAL_STATUS(GreenFunc1BF_real(ri,rj,s,ip,vecTmp1,eleIdx,eleCfg,eleNum,
                                    eleProjCnt,projCntNew,eleProjBFCnt,projBFCntNew,
                                    buffer,&nestedStatus));
    }else if(rj==rl) {
      RETURN_GREENFUNC2BF_REAL(0.0);
    }else if(ri==rl) {
      if(eleNum[rsi]==0) RETURN_GREENFUNC2BF_REAL(0.0);
      else if(rj==rk) RETURN_GREENFUNC2BF_REAL(1.0-eleNum[rsj]);
      else RETURN_GREENFUNC1BF_REAL_STATUS(-GreenFunc1BF_real(rk,rj,s,ip,vecTmp1,eleIdx,eleCfg,eleNum,
                                     eleProjCnt,projCntNew,eleProjBFCnt,projBFCntNew,
                                     buffer,&nestedStatus));
    }else if(rj==rk) {
      if(eleNum[rsj]==1) RETURN_GREENFUNC2BF_REAL(0.0);
      else if(ri==rl) RETURN_GREENFUNC2BF_REAL(eleNum[rsi]);
      else RETURN_GREENFUNC1BF_REAL_STATUS(GreenFunc1BF_real(ri,rl,s,ip,vecTmp1,eleIdx,eleCfg,eleNum,
                                    eleProjCnt,projCntNew,eleProjBFCnt,projBFCntNew,
                                    buffer,&nestedStatus));
    }else if(ri==rk) {
      RETURN_GREENFUNC2BF_REAL(0.0);
    }else if(ri==rj) {
      if(eleNum[rsi]==0) RETURN_GREENFUNC2BF_REAL(0.0);
      else RETURN_GREENFUNC1BF_REAL_STATUS(GreenFunc1BF_real(rk,rl,s,ip,vecTmp1,eleIdx,eleCfg,eleNum,
                                    eleProjCnt,projCntNew,eleProjBFCnt,projBFCntNew,
                                    buffer,&nestedStatus));
    }
  }else{
    if(rk==rl) {
      if(eleNum[rtk]==0) RETURN_GREENFUNC2BF_REAL(0.0);
      else if(ri==rj) RETURN_GREENFUNC2BF_REAL(eleNum[rsi]);
      else RETURN_GREENFUNC1BF_REAL_STATUS(GreenFunc1BF_real(ri,rj,s,ip,vecTmp1,eleIdx,eleCfg,eleNum,
                                    eleProjCnt,projCntNew,eleProjBFCnt,projBFCntNew,
                                    buffer,&nestedStatus));
    }else if(ri==rj) {
      if(eleNum[rsi]==0) RETURN_GREENFUNC2BF_REAL(0.0);
      else RETURN_GREENFUNC1BF_REAL_STATUS(GreenFunc1BF_real(rk,rl,t,ip,vecTmp1,eleIdx,eleCfg,eleNum,
                                    eleProjCnt,projCntNew,eleProjBFCnt,projBFCntNew,
                                    buffer,&nestedStatus));
    }
  }

  if(eleNum[rsi]==1 || eleNum[rsj]==0 || eleNum[rtk]==1 || eleNum[rtl]==0) RETURN_GREENFUNC2BF_REAL(0.0);

  mj = eleCfg[rsj];
  ml = eleCfg[rtl];
  msj = mj + s*Ne;
  mtl = ml + t*Ne;

  if(NQPFull <= 0 || Nsize <= 0) {
    *value = 0.0;
    StopTimer(85);
    return BF_PF_INVALID_ARGUMENT;
  }
  StartTimer(86);
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

  if(!IsSectorStateAllowed(eleNum)) {
    eleCfg[rtl] = ml;
    eleCfg[rtk] = -1;
    eleIdx[mtl] = rl;
    eleNum[rtl] = 1;
    eleNum[rtk] = 0;
    eleCfg[rsj] = mj;
    eleCfg[rsi] = -1;
    eleIdx[msj] = rj;
    eleNum[rsj] = 1;
    eleNum[rsi] = 0;
    StopTimer(86);
    RETURN_GREENFUNC2BF_REAL(0.0);
  }

  z = ProjRatio(projCntNew,eleProjCnt);
  StopTimer(86);

  StartTimer(81);
  MakeProjBFCnt(projBFCntNew, eleNum);
  StopTimer(81);
  if(BFUseCanonicalNonFszPath()) {
    const int rebuildStatus = BFCanonicalGreenInjectFailure("real2")
        ? BF_PF_INVALID_ARGUMENT
        : CalculateBFCanonicalPf_real(
              eleIdx, eleNum, projBFCntNew, 0, NQPFull, pfMNew_real);
    if(rebuildStatus == BF_PF_OK) {
      z *= CalculateIP_real(pfMNew_real, 0, NQPFull, MPI_COMM_SELF);
    }

    eleCfg[rtl] = ml;
    eleCfg[rtk] = -1;
    eleIdx[mtl] = rl;
    eleNum[rtl] = 1;
    eleNum[rtk] = 0;
    eleCfg[rsj] = mj;
    eleCfg[rsi] = -1;
    eleIdx[msj] = rj;
    eleNum[rsj] = 1;
    eleNum[rsi] = 0;
    if(rebuildStatus != BF_PF_OK) {
      *value = 0.0;
      StopTimer(85);
      return rebuildStatus;
    }
    RETURN_GREENFUNC2BF_REAL(z/ip);
  }
  StartTimer(82);
  UpdateSlaterElmBFGrnVec_real(ml, rl, rk, t, eleIdx, eleCfg, eleNum,
                               eleProjBFCnt, projBFCntNew, msaTmp0, icount0, vecTmp0);
  UpdateSlaterElmBFGrnVec_real(mj, rj, ri, s, eleIdx, eleCfg, eleNum,
                               eleProjBFCnt, projBFCntNew, msaTmp1, icount1, vecTmp1);
  StopTimer(82);

  StartTimer(90);
  if (MergeBFHopListsVec_real(msaTmp0, icount0, vecTmp0, msaTmp1, icount1, vecTmp1, msaTmp, icount, vecTmp0) != 0) {
    StopTimer(90);
    StartTimer(88);
    eleCfg[rtl] = ml;
    eleCfg[rtk] = -1;
    eleIdx[mtl] = rl;
    eleNum[rtl] = 1;
    eleNum[rtk] = 0;
    eleCfg[rsj] = mj;
    eleCfg[rsi] = -1;
    eleIdx[msj] = rj;
    eleNum[rsj] = 1;
    eleNum[rsi] = 0;
    StopTimer(88);
    *value = 0.0;
    StopTimer(85);
    return BF_PF_INVALID_ARGUMENT; /* invariant violation, not a physical zero */
  }
  StopTimer(90);

  StartTimer(83);
  CalculateNewPfMBFVecWithStride_real(icount, msaTmp, Nsize, pfMNew_real, 0, NQPFull, vecTmp0, Nsize);
  StopTimer(83);
  StartTimer(87);
  z *= CalculateIP_real(pfMNew_real, 0, NQPFull, MPI_COMM_SELF);
  StopTimer(87);

  StartTimer(88);
  eleCfg[rtl] = ml;
  eleCfg[rtk] = -1;
  eleIdx[mtl] = rl;
  eleNum[rtl] = 1;
  eleNum[rtk] = 0;
  eleCfg[rsj] = mj;
  eleCfg[rsi] = -1;
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;
  StopTimer(88);

  RETURN_GREENFUNC2BF_REAL(z/ip);
#undef RETURN_GREENFUNC1BF_REAL_STATUS
#undef RETURN_GREENFUNC2BF_REAL
}

/* Convenience wrapper with the historical signature: bufM is one packed
 * row-vector workspace (NQPFull*Nsize*Nsize doubles); the second one is
 * allocated here.  Not for hot loops -- CalculateHamiltonianBF_real uses
 * GreenFunc2BF_real_ws with thread-owned workspaces. */
double GreenFunc2BF_real(const int ri, const int rj, const int rk, const int rl,
                    const int s, const int t, const double ip, double *bufM,
                    int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,
                    int *projCntNew, const int *eleProjBFCnt,int *projBFCntNew,
                    double *buffer, int *greenStatus) {
  double value = 0.0;
  double *vecTmp0;
  size_t vecTmpSize;
  int status;

  if(greenStatus == NULL) return NAN;
  *greenStatus = BF_PF_OK;
  if(NQPFull <= 0 || Nsize <= 0
     || (size_t)NQPFull > SIZE_MAX / (size_t)Nsize
     || (size_t)NQPFull * (size_t)Nsize > SIZE_MAX / (size_t)Nsize) {
    *greenStatus = BF_PF_INVALID_ARGUMENT;
    return 0.0;
  }
  vecTmpSize = (size_t)NQPFull * (size_t)Nsize * (size_t)Nsize;
  if(vecTmpSize > SIZE_MAX / sizeof(double)) {
    *greenStatus = BF_PF_INVALID_ARGUMENT;
    return 0.0;
  }
  vecTmp0 = (double *)malloc(sizeof(double)*vecTmpSize);
  if(vecTmp0 == NULL) {
    *greenStatus = BF_PF_INVALID_ARGUMENT;
    return 0.0;
  }
  status = GreenFunc2BF_real_ws(ri, rj, rk, rl, s, t, ip, vecTmp0, bufM,
                                eleIdx, eleCfg, eleNum, eleProjCnt, projCntNew,
                                eleProjBFCnt, projBFCntNew, buffer, &value);
  free(vecTmp0);
  *greenStatus = status;
  return status == BF_PF_OK ? value : 0.0;
}

#ifndef BF_TRANSFER_BATCH_SIZE
#define BF_TRANSFER_BATCH_SIZE 8
#endif

typedef struct {
  BFNBodyStatus status;
  int term;
  BFNBodyStage stage;
  int detail;
} BFNBodyLoopFailure;

typedef struct {
  int status;
  int term;
} BFGreenLoopFailure;

static void RecordBFGreenLoopFailure(BFGreenLoopFailure *failure,
                                     int term, int status) {
  if(failure == NULL || status == BF_PF_OK) return;
#pragma omp critical(BFGreenLoopFailureRecord)
  {
    if(failure->status == BF_PF_OK) {
      failure->status = status;
      failure->term = term;
    }
  }
}

static void AbortBFGreenLoopFailure(const char *context,
                                    const BFGreenLoopFailure *failure) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) {
    fprintf(stderr,
            "Error: BackFlow canonical Green evaluation failed: "
            "context=%s status=%d term=%d.\n",
            context, failure->status, failure->term);
  }
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

static void RecordBFNBodyLoopFailure(BFNBodyLoopFailure *failure, int term,
                                     const BFNBodyResult *result) {
  if(failure == NULL || result == NULL
     || result->status == BF_NBODY_OK
     || result->status == BF_NBODY_PHYSICAL_ZERO) {
    return;
  }
#pragma omp critical(BFNBodyLoopFailureRecord)
  {
    if(failure->status == BF_NBODY_OK) {
      failure->status = result->status;
      failure->term = term;
      failure->stage = result->stage;
      failure->detail = result->detail;
    }
  }
}

static void AbortBFNBodyLoopFailure(const char *context,
                                    const BFNBodyLoopFailure *failure) {
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

double complex CalculateHamiltonianBF_fcmp(const double complex ip, int *eleIdx, const int *eleCfg,
                              int *eleNum, const int *eleProjCnt, const int *eleProjBFCnt)
{
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  const int maxNBody =
      NBodyInterAllMaxN > 2 ? NBodyInterAllMaxN : 2;
  BFNBodyScratchSizes scratchSizes;
  BFNBodyLoopFailure nbodyFailure = {
      BF_NBODY_OK, -1, BF_NBODY_STAGE_NONE, BF_NBODY_DETAIL_NONE
  };
  BFGreenLoopFailure greenFailure = {BF_PF_OK, -1};
  BFGreenLoopFailure twoBodyFailure = {BF_PF_OK, -1};
  double complex e=0.0, tmp;
  double complex greenValue;
  int idx, k, offset, nbody;
  int ri,rj,s,rk,rl,t;
  int workspaceBindingFailed = 0;
  int *myEleIdx, *myEleNum, *myEleCfg, *myProjCntNew, *myProjBFCntNew;
  int *myIntBase;
  double complex *mySltBFTmp;
  double complex *myBuffer;
  double complex *myComplexBase;
  double *myDoubleBase;
  double complex myEnergy;
  BFNBodyScratch scratch;
  BFNBodyResult result;
  int bindStatus;
  int greenStatus;

  if(GetBFNBodyScratchSizes(maxNBody, 0, &scratchSizes) != BF_NBODY_OK) {
    nbodyFailure.status = BF_NBODY_WORKSPACE_ERROR;
    AbortBFNBodyLoopFailure("BackFlow NBodyInterAll workspace sizing",
                            &nbodyFailure);
    return 0.0+0.0*I;
  }
  RequestWorkSpaceThreadInt((int)scratchSizes.intCount);
  RequestWorkSpaceThreadComplex((int)scratchSizes.complexCount);
  RequestWorkSpaceThreadDouble((int)scratchSizes.doubleCount);

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myEleCfg,myProjCntNew,myProjBFCntNew,myIntBase,\
          myBuffer,mySltBFTmp,myComplexBase,myDoubleBase,myEnergy,scratch,\
          result,bindStatus,greenStatus,greenValue,idx)\
  reduction(+:e)
  {
    myEnergy = 0.0;
    myIntBase = GetWorkSpaceThreadInt((int)scratchSizes.intCount);
    myComplexBase =
        GetWorkSpaceThreadComplex((int)scratchSizes.complexCount);
    myDoubleBase = GetWorkSpaceThreadDouble((int)scratchSizes.doubleCount);
    bindStatus = BindBFNBodyScratch(
        &scratchSizes, myIntBase, scratchSizes.intCount,
        myComplexBase, scratchSizes.complexCount,
        myDoubleBase, scratchSizes.doubleCount, &scratch);
    if(bindStatus != BF_NBODY_OK) {
      BFNBodyResult bindResult = {
          BF_NBODY_WORKSPACE_ERROR, BF_NBODY_STAGE_NONE,
          bindStatus, 0, 0.0+0.0*I
      };
      RecordBFNBodyLoopFailure(&nbodyFailure, -1, &bindResult);
#pragma omp atomic write
      workspaceBindingFailed = 1;
    }
#pragma omp barrier

    if(!workspaceBindingFailed) {
      myEleIdx = scratch.eleIdx;
      myEleNum = scratch.eleNum;
      myEleCfg = scratch.eleCfg;
      myProjCntNew = scratch.projCnt;
      myProjBFCntNew = scratch.projBFCnt;
      myBuffer = scratch.greenBuffer;
      mySltBFTmp = scratch.slater;

#pragma loop noalias
      for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
#pragma loop noalias
      for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
#pragma loop noalias
      for(idx=0;idx<Nsite2;idx++) myEleCfg[idx] = eleCfg[idx];

      StoreSlaterElmBF_fcmp(mySltBFTmp);
#pragma omp barrier

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
        myEnergy +=
            ParaCoulombInter[idx] * (n0[ri]+n1[ri]) * (n0[rj]+n1[rj]);
      }

      /* HundCoupling */
#pragma omp for private(idx,ri,rj) nowait
      for(idx=0;idx<NHundCoupling;idx++) {
        ri = HundCoupling[idx][0];
        rj = HundCoupling[idx][1];
        myEnergy -=
            ParaHundCoupling[idx] * (n0[ri]*n0[rj] + n1[ri]*n1[rj]);
        /* Caution: negative sign */
      }

#pragma omp master
      {StopTimer(70);StartTimer(71);}

      /* Transfer */
#pragma omp for private(idx,ri,rj,s) schedule(dynamic) nowait
      for(idx=0;idx<NTransfer;idx++) {
        ri = Transfer[idx][0];
        rj = Transfer[idx][2];
        s  = Transfer[idx][3];

        greenValue = GreenFunc1BF(
            ri,rj,s,ip,mySltBFTmp,myEleIdx,myEleCfg,myEleNum,eleProjCnt,
            myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer,
            &greenStatus);
        if(greenStatus == BF_PF_OK) {
          myEnergy -= ParaTransfer[idx]*greenValue;
        } else {
          RecordBFGreenLoopFailure(&greenFailure, idx, greenStatus);
        }
        /* Caution: negative sign */
      }

#pragma omp master
      {StopTimer(71);StartTimer(72);}

      /* Pair Hopping (BackFlow two-body Green function) */
#pragma omp for private(idx,ri,rj) schedule(dynamic) nowait
      for(idx=0;idx<NPairHopping;idx++) {
        ri = PairHopping[idx][0];
        rj = PairHopping[idx][1];

        greenValue = GreenFunc2BF(ri,rj,ri,rj,0,1,ip,mySltBFTmp,
                                 myEleIdx,myEleCfg,myEleNum,eleProjCnt,
                                 myProjCntNew,eleProjBFCnt,myProjBFCntNew,
                                 myBuffer,&greenStatus);
        if(greenStatus == BF_PF_OK) {
          myEnergy += ParaPairHopping[idx] * greenValue;
        } else {
          RecordBFGreenLoopFailure(&twoBodyFailure, idx, greenStatus);
        }
      }

      /* Exchange Coupling (BackFlow two-body Green function) */
#pragma omp for private(idx,ri,rj,tmp) schedule(dynamic) nowait
      for(idx=0;idx<NExchangeCoupling;idx++) {
        ri = ExchangeCoupling[idx][0];
        rj = ExchangeCoupling[idx][1];

        tmp = GreenFunc2BF(ri,rj,rj,ri,0,1,ip,mySltBFTmp,
                           myEleIdx,myEleCfg,myEleNum,eleProjCnt,
                           myProjCntNew,eleProjBFCnt,myProjBFCntNew,
                           myBuffer,&greenStatus);
        if(greenStatus == BF_PF_OK) {
          tmp += GreenFunc2BF(ri,rj,rj,ri,1,0,ip,mySltBFTmp,
                              myEleIdx,myEleCfg,myEleNum,eleProjCnt,
                              myProjCntNew,eleProjBFCnt,myProjBFCntNew,
                              myBuffer,&greenStatus);
        }
        if(greenStatus == BF_PF_OK) {
          myEnergy += ParaExchangeCoupling[idx] * tmp;
        } else {
          RecordBFGreenLoopFailure(&twoBodyFailure,
                                   NPairHopping+idx, greenStatus);
        }
      }

      /* Inter All (BackFlow two-body Green function) */
#pragma omp for private(idx,ri,rj,s,rk,rl,t) schedule(dynamic) nowait
      for(idx=0;idx<NInterAll;idx++) {
        ri = InterAll[idx][0];
        rj = InterAll[idx][2];
        s  = InterAll[idx][3];
        rk = InterAll[idx][4];
        rl = InterAll[idx][6];
        t  = InterAll[idx][7];
        greenValue = GreenFunc2BF(ri,rj,rk,rl,s,t,ip,mySltBFTmp,
                                 myEleIdx,myEleCfg,myEleNum,eleProjCnt,
                                 myProjCntNew,eleProjBFCnt,myProjBFCntNew,
                                 myBuffer,&greenStatus);
        if(greenStatus == BF_PF_OK) {
          myEnergy += ParaInterAll[idx] * greenValue;
        } else {
          RecordBFGreenLoopFailure(&twoBodyFailure,
                                   NPairHopping+NExchangeCoupling+idx,
                                   greenStatus);
        }
      }

      /*
       * scratch.slater aliases mySltBFTmp. N>=2 evaluation replaces that
       * base snapshot, so this loop must remain after all legacy consumers.
       * GreenFuncNBF reconstructs the state needed by every term.
       */
#pragma omp for private(idx,k,offset,nbody,result) schedule(dynamic)
      for(idx=0;idx<NNBodyInterAll;idx++) {
        scratch.termIndex = idx;
        nbody = NBodyInterAllN[idx];
        offset = NBodyInterAllOffset[idx];
        for(k=0;k<nbody;k++) {
          scratch.inputRsi[k] =
              NBodyInterAllIdx[offset+k][0]
              + NBodyInterAllIdx[offset+k][1]*Nsite;
          scratch.inputRsj[k] =
              NBodyInterAllIdx[offset+k][2]
              + NBodyInterAllIdx[offset+k][3]*Nsite;
        }
        result = GreenFuncNBF(
            nbody, scratch.inputRsi, scratch.inputRsj, ip,
            eleIdx, eleCfg, eleNum, eleProjCnt, eleProjBFCnt,
            &scratch);
        if(result.status == BF_NBODY_OK
           || result.status == BF_NBODY_PHYSICAL_ZERO) {
          myEnergy += ParaNBodyInterAll[idx] * result.value;
        } else {
          RecordBFNBodyLoopFailure(&nbodyFailure, idx, &result);
        }
      }

#pragma omp master
      {StopTimer(72);}
    }

    e += myEnergy;
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  ReleaseWorkSpaceThreadDouble();
  if(greenFailure.status != BF_PF_OK) {
    AbortBFGreenLoopFailure("Hamiltonian complex Transfer", &greenFailure);
  }
  if(twoBodyFailure.status != BF_PF_OK) {
    AbortBFGreenLoopFailure("Hamiltonian complex two-body", &twoBodyFailure);
  }
  if(nbodyFailure.status != BF_NBODY_OK) {
    AbortBFNBodyLoopFailure("BackFlow NBodyInterAll evaluation",
                            &nbodyFailure);
  }
  return e;
}
double CalculateHamiltonianBF_real(const double ip, int *eleIdx, const int *eleCfg,
                                   int *eleNum, const int *eleProjCnt, const int *eleProjBFCnt) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  BFGreenLoopFailure greenFailure = {BF_PF_OK, -1};
  double e=0.0, tmp;
  double greenValue;
  int idx;
  int ri,rj,s,rk,rl,t;
  int batchNo,batchStart,batchEnd,batchIdx,batchCount,nTransferBatch;
  int *myEleIdx, *myEleNum, *myEleCfg, *myProjCntNew, *myProjBFCntNew;
  int *myBatchIdx, *myBatchMsa, *myBatchIcount;
  //double sltTmp[NThread*NQPFull*Nsite2*Nsite2];
  double *myBatchVec, *myBatchPfM, *myBatchProj, *myBatchGreen, *myBatchVecStack, *myBatchWStack;
  double *myBuffer;
  double *myVecTmp0, *myVecTmp1; /* per-thread packed row-vector workspaces for GreenFunc2BF_real_ws */
  const int needTwoBody = (NPairHopping > 0 || NExchangeCoupling > 0 || NInterAll > 0);
  long long vecTmpSizeLL = needTwoBody ? (long long)NQPFull*(long long)Nsize*(long long)Nsize : 0;
  int vecTmpSize;
  BFGreenLoopFailure twoBodyFailure = {BF_PF_OK, -1};
  int greenStatus;
  double green;
  double myEnergy;
  if(vecTmpSizeLL < 0 || vecTmpSizeLL > INT_MAX/2) {
    fprintf(stderr, "error: CalculateHamiltonianBF_real workspace size overflow (NQPFull=%d, Nsize=%d)\n", NQPFull, Nsize);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  vecTmpSize = (int)vecTmpSizeLL;

  RequestWorkSpaceThreadInt(Nsize+2*Nsite2+NProj+16*Nsite*Nrange
                            + BF_TRANSFER_BATCH_SIZE
                            + BF_TRANSFER_BATCH_SIZE*NQPFull
                            + BF_TRANSFER_BATCH_SIZE*NQPFull*Nsize);
  RequestWorkSpaceThreadDouble(NQPFull+2*Nsize + 2*vecTmpSize
                               + BF_TRANSFER_BATCH_SIZE*NQPFull*Nsize*Nsize
                               + BF_TRANSFER_BATCH_SIZE*NQPFull
                               + 2*BF_TRANSFER_BATCH_SIZE
                               + 2*BF_TRANSFER_BATCH_SIZE*Nsize*Nsize);
  nTransferBatch = (NTransfer + BF_TRANSFER_BATCH_SIZE - 1) / BF_TRANSFER_BATCH_SIZE;
  /* GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize, BF transfer batch scratch */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myEleCfg,myProjCntNew,myProjBFCntNew,myBuffer,myVecTmp0,myVecTmp1,greenStatus,green,\
          myBatchIdx,myBatchMsa,myBatchIcount,myBatchVec,myBatchPfM,myBatchProj,myBatchGreen,\
          myBatchVecStack,myBatchWStack,myEnergy,idx,ri,rj,s,rk,rl,t,tmp,\
          greenValue,batchNo,batchStart,batchEnd,batchIdx,batchCount)\
  reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myEleCfg = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew   = GetWorkSpaceThreadInt(NProj);
    myProjBFCntNew = GetWorkSpaceThreadInt(16*Nsite*Nrange);
    myBatchIdx     = GetWorkSpaceThreadInt(BF_TRANSFER_BATCH_SIZE);
    myBatchIcount  = GetWorkSpaceThreadInt(BF_TRANSFER_BATCH_SIZE*NQPFull);
    myBatchMsa     = GetWorkSpaceThreadInt(BF_TRANSFER_BATCH_SIZE*NQPFull*Nsize);
    myBuffer   = GetWorkSpaceThreadDouble(NQPFull+2*Nsize);
    myVecTmp0 = needTwoBody ? GetWorkSpaceThreadDouble(vecTmpSize) : NULL;
    myVecTmp1 = needTwoBody ? GetWorkSpaceThreadDouble(vecTmpSize) : NULL;
    myBatchVec = GetWorkSpaceThreadDouble(BF_TRANSFER_BATCH_SIZE*NQPFull*Nsize*Nsize);
    myBatchPfM = GetWorkSpaceThreadDouble(BF_TRANSFER_BATCH_SIZE*NQPFull);
    myBatchProj = GetWorkSpaceThreadDouble(BF_TRANSFER_BATCH_SIZE);
    myBatchGreen = GetWorkSpaceThreadDouble(BF_TRANSFER_BATCH_SIZE);
    myBatchVecStack = GetWorkSpaceThreadDouble(BF_TRANSFER_BATCH_SIZE*Nsize*Nsize);
    myBatchWStack = GetWorkSpaceThreadDouble(BF_TRANSFER_BATCH_SIZE*Nsize*Nsize);

#pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
#pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
#pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleCfg[idx] = eleCfg[idx];

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
#pragma omp for private(idx,ri,rj,s) schedule(dynamic) nowait
    for(batchNo=0;batchNo<nTransferBatch;batchNo++) {
      batchStart = batchNo * BF_TRANSFER_BATCH_SIZE;
      batchEnd = batchStart + BF_TRANSFER_BATCH_SIZE;
      if(batchEnd > NTransfer) batchEnd = NTransfer;
      batchCount = 0;

      for(idx=batchStart;idx<batchEnd;idx++) {
        ri = Transfer[idx][0];
        rj = Transfer[idx][2];
        s  = Transfer[idx][3];

        if(BFUseCanonicalNonFszPath()) {
          greenValue = GreenFunc1BF_real(
              ri, rj, s, ip, SlaterElmBF_real,
              myEleIdx, myEleCfg, myEleNum, eleProjCnt,
              myProjCntNew, eleProjBFCnt, myProjBFCntNew, myBuffer,
              &greenStatus);
          if(greenStatus == BF_PF_OK) {
            myEnergy -= creal(ParaTransfer[idx])*greenValue;
          } else {
            RecordBFGreenLoopFailure(&greenFailure, idx, greenStatus);
          }
          continue;
        }

        if(GreenFunc1BF_real_prepare(ri,rj,s,myBatchGreen+batchCount,myBatchProj+batchCount,
                                     myBatchVec + batchCount*NQPFull*Nsize*Nsize,
                                     myBatchMsa + batchCount*NQPFull*Nsize,
                                     myBatchIcount + batchCount*NQPFull,
                                     myEleIdx,myEleCfg,myEleNum,eleProjCnt,myProjCntNew,
                                     eleProjBFCnt,myProjBFCntNew)) {
          myBatchIdx[batchCount] = idx;
          batchCount++;
        } else {
          myEnergy -= creal(ParaTransfer[idx]) * myBatchGreen[batchCount];
        }
        /* Caution: negative sign */
      }

      if(batchCount > 0) {
        GreenFunc1BF_real_finish_batch(batchCount, ip, myBatchProj, myBatchIcount, myBatchMsa,
                                       myBatchVec, myBatchGreen, myBatchPfM,
                                       myBatchVecStack, myBatchWStack);
        for(batchIdx=0;batchIdx<batchCount;batchIdx++) {
          idx = myBatchIdx[batchIdx];
          myEnergy -= creal(ParaTransfer[idx]) * myBatchGreen[batchIdx];
        }
        /* Caution: negative sign */
      }
    }

#pragma omp master
    {StopTimer(71);StartTimer(72);}

    /* Pair Hopping (BackFlow two-body Green function) */
#pragma omp for private(idx,ri,rj) schedule(dynamic) nowait
    for(idx=0;idx<NPairHopping;idx++) {
      ri = PairHopping[idx][0];
      rj = PairHopping[idx][1];

      greenStatus = GreenFunc2BF_real_ws(ri,rj,ri,rj,0,1,ip,myVecTmp0,myVecTmp1,myEleIdx,myEleCfg,myEleNum,
                                         eleProjCnt,myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer,&green);
      if(greenStatus != BF_PF_OK) {
        RecordBFGreenLoopFailure(&twoBodyFailure, idx, greenStatus);
        continue;
      }
      myEnergy += ParaPairHopping[idx] * green;
    }

    /* Exchange Coupling (BackFlow two-body Green function) */
#pragma omp for private(idx,ri,rj,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NExchangeCoupling;idx++) {
      ri = ExchangeCoupling[idx][0];
      rj = ExchangeCoupling[idx][1];

      greenStatus = GreenFunc2BF_real_ws(ri,rj,rj,ri,0,1,ip,myVecTmp0,myVecTmp1,myEleIdx,myEleCfg,myEleNum,
                                         eleProjCnt,myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer,&green);
      tmp = green;
      if(greenStatus == BF_PF_OK) {
        greenStatus = GreenFunc2BF_real_ws(ri,rj,rj,ri,1,0,ip,myVecTmp0,myVecTmp1,myEleIdx,myEleCfg,myEleNum,
                                            eleProjCnt,myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer,&green);
      }
      if(greenStatus != BF_PF_OK) {
        RecordBFGreenLoopFailure(&twoBodyFailure, NPairHopping+idx, greenStatus);
        continue;
      }
      tmp += green;
      myEnergy += ParaExchangeCoupling[idx] * tmp;
    }

    /* Inter All (BackFlow two-body Green function) */
#pragma omp for private(idx,ri,rj,s,rk,rl,t) schedule(dynamic) nowait
    for(idx=0;idx<NInterAll;idx++) {
      ri = InterAll[idx][0];
      rj = InterAll[idx][2];
      s  = InterAll[idx][3];
      rk = InterAll[idx][4];
      rl = InterAll[idx][6];
      t  = InterAll[idx][7];
      greenStatus = GreenFunc2BF_real_ws(ri,rj,rk,rl,s,t,ip,myVecTmp0,myVecTmp1,myEleIdx,myEleCfg,myEleNum,
                                         eleProjCnt,myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer,&green);
      if(greenStatus != BF_PF_OK) {
        RecordBFGreenLoopFailure(&twoBodyFailure,
                                 NPairHopping+NExchangeCoupling+idx,
                                 greenStatus);
        continue;
      }
      myEnergy += creal(ParaInterAll[idx]) * green;
    }

#pragma omp master
    {StopTimer(72);}

    e += myEnergy;
  }
  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadDouble();
  if(greenFailure.status != BF_PF_OK) {
    AbortBFGreenLoopFailure("Hamiltonian real Transfer", &greenFailure);
  }
  if(twoBodyFailure.status != BF_PF_OK) {
    AbortBFGreenLoopFailure("Hamiltonian real two-body", &twoBodyFailure);
  }
  return e;
}

/* Calculate the CoulombIntra, CoulombInter, Hund terms, */
/* which can be calculated by number operators. */
///
/// \param eleNum [in]
/// \return myEnergy
/// \version 1.0
void CalculateGreenFuncBF(const double w, const double complex ip, int *eleIdx, int *eleCfg,
                          int *eleNum, int *eleProjCnt, const int *eleProjBFCnt) {

  const int maxNBody = NBodyGMaxN > 2 ? NBodyGMaxN : 2;
  BFNBodyScratchSizes scratchSizes;
  BFNBodyLoopFailure nbodyFailure = {
      BF_NBODY_OK, -1, BF_NBODY_STAGE_NONE, BF_NBODY_DETAIL_NONE
  };
  BFGreenLoopFailure greenFailure = {BF_PF_OK, -1};
  int idx,idx0,idx1,k,offset,nbody;
  int ri,rj,s,rk,rl,t;
  int workspaceBindingFailed = 0;
  double complex tmp;
  int *myEleIdx, *myEleNum, *myEleCfg, *myProjCntNew, *myProjBFCntNew;
  int *myIntBase;
  double complex* mySltBFTmp;
  double complex* myBuffer;
  double complex *myComplexBase;
  double *myDoubleBase;
  BFNBodyScratch scratch;
  BFNBodyResult result;
  int bindStatus;
  int greenStatus;

  if(GetBFNBodyScratchSizes(maxNBody, 0, &scratchSizes) != BF_NBODY_OK) {
    nbodyFailure.status = BF_NBODY_WORKSPACE_ERROR;
    AbortBFNBodyLoopFailure("BackFlow NBodyG workspace sizing",
                            &nbodyFailure);
    return;
  }
  RequestWorkSpaceThreadInt((int)scratchSizes.intCount);
  RequestWorkSpaceThreadComplex((int)scratchSizes.complexCount);
  RequestWorkSpaceThreadDouble((int)scratchSizes.doubleCount);

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myEleCfg,myProjCntNew,myProjBFCntNew,myIntBase,\
          myBuffer,mySltBFTmp,myComplexBase,myDoubleBase,scratch,result,\
          bindStatus,greenStatus,idx)
  {
    myIntBase = GetWorkSpaceThreadInt((int)scratchSizes.intCount);
    myComplexBase =
        GetWorkSpaceThreadComplex((int)scratchSizes.complexCount);
    myDoubleBase = GetWorkSpaceThreadDouble((int)scratchSizes.doubleCount);
    bindStatus = BindBFNBodyScratch(
        &scratchSizes, myIntBase, scratchSizes.intCount,
        myComplexBase, scratchSizes.complexCount,
        myDoubleBase, scratchSizes.doubleCount, &scratch);
    if(bindStatus != BF_NBODY_OK) {
      BFNBodyResult bindResult = {
          BF_NBODY_WORKSPACE_ERROR, BF_NBODY_STAGE_NONE,
          bindStatus, 0, 0.0+0.0*I
      };
      RecordBFNBodyLoopFailure(&nbodyFailure, -1, &bindResult);
#pragma omp atomic write
      workspaceBindingFailed = 1;
    }
#pragma omp barrier

    if(!workspaceBindingFailed) {
      myEleIdx = scratch.eleIdx;
      myEleNum = scratch.eleNum;
      myEleCfg = scratch.eleCfg;
      myProjCntNew = scratch.projCnt;
      myProjBFCntNew = scratch.projBFCnt;
      myBuffer = scratch.greenBuffer;
      mySltBFTmp = scratch.slater;

#pragma loop noalias
      for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
#pragma loop noalias
      for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
#pragma loop noalias
      for(idx=0;idx<Nsite2;idx++) myEleCfg[idx] = eleCfg[idx];

      StoreSlaterElmBF_fcmp(mySltBFTmp);
#pragma omp master
      {StartTimer(50);}

#pragma omp for private(idx,ri,rj,s,tmp) schedule(dynamic) nowait
      for(idx=0;idx<NCisAjs;idx++) {
        ri = CisAjsIdx[idx][0];
        rj = CisAjsIdx[idx][2];
        s  = CisAjsIdx[idx][3];
        tmp = GreenFunc1BF(ri,rj,s,ip,mySltBFTmp,myEleIdx,myEleCfg,
                           myEleNum,eleProjCnt,myProjCntNew,
                           eleProjBFCnt,myProjBFCntNew,myBuffer,
                           &greenStatus);
        if(greenStatus == BF_PF_OK) {
          LocalCisAjs[idx] = tmp;
        } else {
          LocalCisAjs[idx] = 0.0+0.0*I;
          RecordBFGreenLoopFailure(&greenFailure, idx, greenStatus);
        }
      }

#pragma omp master
      {StopTimer(50);StartTimer(51);}

#pragma omp for private(idx,ri,rj,s,rk,rl,t,tmp) schedule(dynamic)
      for(idx=0;idx<NCisAjsCktAltDC;idx++) {
        ri = CisAjsCktAltDCIdx[idx][0];
        rj = CisAjsCktAltDCIdx[idx][2];
        s  = CisAjsCktAltDCIdx[idx][1];
        rk = CisAjsCktAltDCIdx[idx][4];
        rl = CisAjsCktAltDCIdx[idx][6];
        t  = CisAjsCktAltDCIdx[idx][5];
        tmp = GreenFunc2BF(ri,rj,rk,rl,s,t,ip,mySltBFTmp,
                           myEleIdx,myEleCfg,myEleNum,eleProjCnt,
                           myProjCntNew,eleProjBFCnt,myProjBFCntNew,
                           myBuffer,&greenStatus);
        if(greenStatus == BF_PF_OK) {
          LocalCisAjsCktAltDC[idx] = tmp;
        } else {
          LocalCisAjsCktAltDC[idx] = 0.0+0.0*I;
          RecordBFGreenLoopFailure(&greenFailure, idx, greenStatus);
        }
      }

#pragma omp master
      {StopTimer(51);StartTimer(52);}

#pragma omp for private(idx) nowait
      for(idx=0;idx<NCisAjs;idx++) {
        PhysCisAjs[idx] += w*LocalCisAjs[idx];
      }

#pragma omp for private(idx) nowait
      for (idx=0;idx<NCisAjsCktAltDC;idx++) {
        PhysCisAjsCktAltDC[idx] += w*LocalCisAjsCktAltDC[idx];
      }

#pragma omp master
      {StopTimer(52);StartTimer(53);}

#pragma omp for private(idx,idx0,idx1) nowait
      for(idx=0;idx<NCisAjsCktAlt;idx++) {
        idx0 = CisAjsCktAltIdx[idx][0];
        idx1 = CisAjsCktAltIdx[idx][1];
        PhysCisAjsCktAlt[idx] +=
            w*LocalCisAjs[idx0]*conj(LocalCisAjs[idx1]);
      }

#pragma omp master
      {StopTimer(53);StartTimer(54);}

      /*
       * scratch.slater aliases mySltBFTmp. N>=2 evaluation replaces that
       * base snapshot, so this loop must remain after every legacy and
       * derived-observable consumer. GreenFuncNBF reconstructs the state
       * needed by every term.
       */
#pragma omp for private(idx,k,offset,nbody,result) schedule(dynamic)
      for(idx=0;idx<NNBodyG;idx++) {
        scratch.termIndex = idx;
        nbody = NBodyGN[idx];
        offset = NBodyGOffset[idx];
        for(k=0;k<nbody;k++) {
          scratch.inputRsi[k] =
              NBodyGIdx[offset+k][0] + NBodyGIdx[offset+k][1]*Nsite;
          scratch.inputRsj[k] =
              NBodyGIdx[offset+k][2] + NBodyGIdx[offset+k][3]*Nsite;
        }
        result = GreenFuncNBF(
            nbody, scratch.inputRsi, scratch.inputRsj, ip,
            eleIdx, eleCfg, eleNum, eleProjCnt, eleProjBFCnt,
            &scratch);
        if(result.status == BF_NBODY_OK
           || result.status == BF_NBODY_PHYSICAL_ZERO) {
          PhysNBodyG[idx] += w*result.value;
        } else {
          RecordBFNBodyLoopFailure(&nbodyFailure, idx, &result);
        }
      }

#pragma omp master
      {StopTimer(54);}
    }
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  ReleaseWorkSpaceThreadDouble();
  if(greenFailure.status != BF_PF_OK) {
    AbortBFGreenLoopFailure("observable Green function", &greenFailure);
  }
  if(nbodyFailure.status != BF_NBODY_OK) {
    AbortBFNBodyLoopFailure("BackFlow NBodyG evaluation", &nbodyFailure);
  }
  return;
}
