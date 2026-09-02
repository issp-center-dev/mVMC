/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
*/

#include <complex.h>

#include "include/gc_size.h"
#include "include/global.h"
#include "include/locgrn_gc.h"
#include "include/workspace.h"

#ifndef _SRC_CALHAM_GC
#define _SRC_CALHAM_GC

static void GCHamiltonianWorkspaceCounts(
    const int maxN, int *intCount, int *complexCount, int *doubleCount,
    int *smallCount, int *twiceMaxCount, int *vecCount) {
  const size_t maxNSize = (size_t)maxN;
  const size_t twiceMax = GCCheckedMulSize(2, (size_t)maxN);
  const size_t small = GCCheckedMulSize(twiceMax, twiceMax);
  const size_t vec =
      GCCheckedMulSize(maxNSize, (size_t)NsizeMax);
  size_t ints = (size_t)NsizeMax;
  size_t complexes = (size_t)NQPFull;
  ints = GCCheckedAddSize(ints,
                          GCCheckedMulSize(2, (size_t)Nsite2));
  ints = GCCheckedAddSize(ints, (size_t)NProj);
  ints = GCCheckedAddSize(ints, GCCheckedMulSize(5, (size_t)maxN));
  complexes = GCCheckedAddSize(
      complexes,
      GCCheckedMulSize(GCCheckedAddSize(maxNSize, 1),
                       (size_t)NsizeMax));
  complexes = GCCheckedAddSize(complexes, GCCheckedMulSize(2, small));
  *intCount = GCCheckedSizeToInt(ints, "GC Hamiltonian integer workspace");
  *complexCount =
      GCCheckedSizeToInt(complexes, "GC Hamiltonian complex workspace");
  *doubleCount =
      GCCheckedSizeToInt(small, "GC Hamiltonian real workspace");
  *smallCount = GCCheckedSizeToInt(small, "GC Hamiltonian small matrix");
  *twiceMaxCount =
      GCCheckedSizeToInt(twiceMax, "GC Hamiltonian Pfaffian integer work");
  *vecCount = GCCheckedSizeToInt(vec, "GC Hamiltonian vector workspace");
}

double CalculateSzGC(const int *eleNum) {
  double sz = 0.0;
  int site;
  for (site = 0; site < Nsite; site++) {
    sz += (double)(eleNum[site] - eleNum[site + Nsite]);
  }
  return 0.5 * sz;
}

double complex CalculateHamiltonianGC(const double complex ip,
                                      const int ncur, int *eleIdx,
                                      int *eleCfg, int *eleNum,
                                      const int *eleProjCnt) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  const int maxN = NBodyInterAllMaxN > 2 ? NBodyInterAllMaxN : 2;
  int intCount;
  int complexCount;
  int doubleCount;
  int smallCount;
  int twiceMaxCount;
  int vecCount;
  double complex energy = 0.0;

  GCHamiltonianWorkspaceCounts(
      maxN, &intCount, &complexCount, &doubleCount, &smallCount,
      &twiceMaxCount, &vecCount);
  RequestWorkSpaceThreadInt(intCount);
  RequestWorkSpaceThreadComplex(complexCount);
  RequestWorkSpaceThreadDouble(doubleCount);

#pragma omp parallel reduction(+ : energy)
  {
    int *myEleIdx = GetWorkSpaceThreadInt(NsizeMax);
    int *myEleCfg = GetWorkSpaceThreadInt(Nsite2);
    int *myEleNum = GetWorkSpaceThreadInt(Nsite2);
    GCGreenScratch scratch;
    double complex myEnergy = 0.0;
    int idx;
    scratch.maxN = maxN;
    scratch.pfLWork = smallCount;
    scratch.projCntNew = GetWorkSpaceThreadInt(NProj);
    scratch.rsi = GetWorkSpaceThreadInt(maxN);
    scratch.rsj = GetWorkSpaceThreadInt(maxN);
    scratch.msa = GetWorkSpaceThreadInt(maxN);
    scratch.pfIWork = GetWorkSpaceThreadInt(twiceMaxCount);
    scratch.pfMNew = GetWorkSpaceThreadComplex(NQPFull);
    scratch.vec = GetWorkSpaceThreadComplex(vecCount);
    scratch.w = GetWorkSpaceThreadComplex(NsizeMax);
    scratch.smallMat = GetWorkSpaceThreadComplex(smallCount);
    scratch.pfWork = GetWorkSpaceThreadComplex(smallCount);
    scratch.pfRWork = GetWorkSpaceThreadDouble(smallCount);
    for (idx = 0; idx < NsizeMax; idx++) myEleIdx[idx] = eleIdx[idx];
    for (idx = 0; idx < Nsite2; idx++) {
      myEleCfg[idx] = eleCfg[idx];
      myEleNum[idx] = eleNum[idx];
    }

#pragma omp master
    { StartTimer(70); }
#pragma omp for
    for (idx = 0; idx < NCoulombIntra; idx++) {
      const int ri = CoulombIntra[idx];
      myEnergy += ParaCoulombIntra[idx] * n0[ri] * n1[ri];
    }
#pragma omp for
    for (idx = 0; idx < NCoulombInter; idx++) {
      const int ri = CoulombInter[idx][0];
      const int rj = CoulombInter[idx][1];
      myEnergy += ParaCoulombInter[idx] * (n0[ri] + n1[ri]) *
                  (n0[rj] + n1[rj]);
    }
#pragma omp for
    for (idx = 0; idx < NHundCoupling; idx++) {
      const int ri = HundCoupling[idx][0];
      const int rj = HundCoupling[idx][1];
      myEnergy -= ParaHundCoupling[idx] *
                  (n0[ri] * n0[rj] + n1[ri] * n1[rj]);
    }
#pragma omp master
    {
      StopTimer(70);
      StartTimer(71);
    }
#pragma omp for schedule(dynamic)
    for (idx = 0; idx < NTransfer; idx++) {
      const int rsi = Transfer[idx][0] + Transfer[idx][1] * Nsite;
      const int rsj = Transfer[idx][2] + Transfer[idx][3] * Nsite;
      myEnergy -= ParaTransfer[idx] *
                  GreenFunc1GC(rsi, rsj, ip, ncur, myEleIdx, myEleCfg,
                               myEleNum, eleProjCnt, &scratch);
    }
#pragma omp master
    {
      StopTimer(71);
      StartTimer(72);
    }
#pragma omp for schedule(dynamic)
    for (idx = 0; idx < NPairHopping; idx++) {
      const int ri = PairHopping[idx][0];
      const int rj = PairHopping[idx][1];
      myEnergy += ParaPairHopping[idx] * GreenFunc2GC(
          ri, rj, ri + Nsite, rj + Nsite, ip, ncur, myEleIdx, myEleCfg,
          myEleNum, eleProjCnt, &scratch);
    }
#pragma omp for schedule(dynamic)
    for (idx = 0; idx < NExchangeCoupling; idx++) {
      const int ri = ExchangeCoupling[idx][0];
      const int rj = ExchangeCoupling[idx][1];
      double complex value = GreenFunc2GC(
          ri, rj, rj + Nsite, ri + Nsite, ip, ncur, myEleIdx, myEleCfg,
          myEleNum, eleProjCnt, &scratch);
      value += GreenFunc2GC(ri + Nsite, rj + Nsite, rj, ri, ip, ncur,
                            myEleIdx, myEleCfg, myEleNum, eleProjCnt,
                            &scratch);
      myEnergy += ParaExchangeCoupling[idx] * value;
    }
#pragma omp for schedule(dynamic)
    for (idx = 0; idx < NInterAll; idx++) {
      const int rsi = InterAll[idx][0] + InterAll[idx][1] * Nsite;
      const int rsj = InterAll[idx][2] + InterAll[idx][3] * Nsite;
      const int rsk = InterAll[idx][4] + InterAll[idx][5] * Nsite;
      const int rsl = InterAll[idx][6] + InterAll[idx][7] * Nsite;
      myEnergy += ParaInterAll[idx] * GreenFunc2GC(
          rsi, rsj, rsk, rsl, ip, ncur, myEleIdx, myEleCfg, myEleNum,
          eleProjCnt, &scratch);
    }
#pragma omp for schedule(dynamic)
    for (idx = 0; idx < NNBodyInterAll; idx++) {
      const int nbody = NBodyInterAllN[idx];
      const int offset = NBodyInterAllOffset[idx];
      int k;
      for (k = 0; k < nbody; k++) {
        scratch.rsi[k] = NBodyInterAllIdx[offset + k][0] +
                         NBodyInterAllIdx[offset + k][1] * Nsite;
        scratch.rsj[k] = NBodyInterAllIdx[offset + k][2] +
                         NBodyInterAllIdx[offset + k][3] * Nsite;
      }
      myEnergy += ParaNBodyInterAll[idx] * GreenFuncNGC(
          nbody, scratch.rsi, scratch.rsj, ip, ncur, myEleIdx, myEleCfg,
          myEleNum, eleProjCnt, &scratch);
    }
#pragma omp for schedule(dynamic)
    for (idx = 0; idx < NAnomalousTerm; idx++) {
      const int type = AnomalousTerm[idx][0];
      const int rs1 =
          AnomalousTerm[idx][1] + AnomalousTerm[idx][2] * Nsite;
      const int rs2 =
          AnomalousTerm[idx][3] + AnomalousTerm[idx][4] * Nsite;
      const double complex green =
          type == 1
              ? GreenFuncPairAddGC(rs1, rs2, ip, ncur, myEleIdx,
                                   myEleCfg, myEleNum, eleProjCnt, &scratch)
              : GreenFuncPairRemoveGC(rs1, rs2, ip, ncur, myEleIdx,
                                      myEleCfg, myEleNum, eleProjCnt,
                                      &scratch);
      myEnergy += ParaAnomalousTerm[idx] * green;
    }
#pragma omp master
    { StopTimer(72); }
    energy += myEnergy;
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  ReleaseWorkSpaceThreadDouble();
  return energy;
}

#endif
