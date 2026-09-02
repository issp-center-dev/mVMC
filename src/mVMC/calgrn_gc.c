/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
*/

#include <complex.h>
#include <math.h>

#include "include/gc_size.h"
#include "include/global.h"
#include "include/locgrn_gc.h"
#include "include/workspace.h"

#ifndef _SRC_CALGRN_GC
#define _SRC_CALGRN_GC

static void GCGreenWorkspaceCounts(const int maxN, int *intCount,
                                   int *complexCount, int *doubleCount,
                                   int *smallCount, int *twiceMaxCount,
                                   int *vecCount) {
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
  *intCount = GCCheckedSizeToInt(ints, "GC measurement integer workspace");
  *complexCount =
      GCCheckedSizeToInt(complexes, "GC measurement complex workspace");
  *doubleCount =
      GCCheckedSizeToInt(small, "GC measurement real workspace");
  *smallCount = GCCheckedSizeToInt(small, "GC measurement small matrix");
  *twiceMaxCount =
      GCCheckedSizeToInt(twiceMax, "GC measurement Pfaffian integer work");
  *vecCount = GCCheckedSizeToInt(vec, "GC measurement vector workspace");
}

void CalculateGreenFuncGC(const double w, const double complex ip,
                          const int ncur, int *eleIdx, int *eleCfg,
                          int *eleNum, int *eleProjCnt) {
  const int maxN = NBodyGMaxN > 2 ? NBodyGMaxN : 2;
  const double twoPi = 2.0 * acos(-1.0);
  int intCount;
  int complexCount;
  int doubleCount;
  int smallCount;
  int twiceMaxCount;
  int vecCount;

  GCGreenWorkspaceCounts(maxN, &intCount, &complexCount, &doubleCount,
                         &smallCount, &twiceMaxCount, &vecCount);
  RequestWorkSpaceThreadInt(intCount);
  RequestWorkSpaceThreadComplex(complexCount);
  RequestWorkSpaceThreadDouble(doubleCount);

#pragma omp parallel
  {
    int *myEleIdx = GetWorkSpaceThreadInt(NsizeMax);
    int *myEleCfg = GetWorkSpaceThreadInt(Nsite2);
    int *myEleNum = GetWorkSpaceThreadInt(Nsite2);
    GCGreenScratch scratch;
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
    { StartTimer(50); }
#pragma omp for schedule(dynamic)
    for (idx = 0; idx < NCisAjs; idx++) {
      const int rsi = CisAjsIdx[idx][0] + CisAjsIdx[idx][1] * Nsite;
      const int rsj = CisAjsIdx[idx][2] + CisAjsIdx[idx][3] * Nsite;
      LocalCisAjs[idx] = GreenFunc1GC(
          rsi, rsj, ip, ncur, myEleIdx, myEleCfg, myEleNum, eleProjCnt,
          &scratch);
    }
#pragma omp master
    {
      StopTimer(50);
      StartTimer(51);
    }
#pragma omp for schedule(dynamic)
    for (idx = 0; idx < NCisAjsCktAltDC; idx++) {
      const int rsi = CisAjsCktAltDCIdx[idx][0] +
                      CisAjsCktAltDCIdx[idx][1] * Nsite;
      const int rsj = CisAjsCktAltDCIdx[idx][2] +
                      CisAjsCktAltDCIdx[idx][3] * Nsite;
      const int rsk = CisAjsCktAltDCIdx[idx][4] +
                      CisAjsCktAltDCIdx[idx][5] * Nsite;
      const int rsl = CisAjsCktAltDCIdx[idx][6] +
                      CisAjsCktAltDCIdx[idx][7] * Nsite;
      PhysCisAjsCktAltDC[idx] +=
          w * GreenFunc2GC(rsi, rsj, rsk, rsl, ip, ncur, myEleIdx,
                           myEleCfg, myEleNum, eleProjCnt, &scratch);
    }
#pragma omp master
    {
      StopTimer(51);
      StartTimer(52);
    }
#pragma omp for
    for (idx = 0; idx < NCisAjs; idx++) {
      PhysCisAjs[idx] += w * LocalCisAjs[idx];
    }
#pragma omp master
    {
      StopTimer(52);
      StartTimer(53);
    }
#pragma omp for
    for (idx = 0; idx < NCisAjsCktAlt; idx++) {
      const int idx0 = CisAjsCktAltIdx[idx][0];
      const int idx1 = CisAjsCktAltIdx[idx][1];
      PhysCisAjsCktAlt[idx] +=
          w * LocalCisAjs[idx0] * conj(LocalCisAjs[idx1]);
    }
#pragma omp master
    {
      StopTimer(53);
      StartTimer(54);
    }
#pragma omp for
    for (idx = 0; idx < NTwist; idx++) {
      double twistWeight = 0.0;
      int fused;
      for (fused = 0; fused < Nsite2; fused++) {
        const int site = TwistIdx[idx][2 * fused];
        const int spin = TwistIdx[idx][2 * fused + 1];
        const int rsi = site + spin * Nsite;
        const int x = LatticeIdx[site][0];
        const int y = LatticeIdx[site][1];
        const int z = LatticeIdx[site][2];
        twistWeight +=
            (x * ParaTwist[idx][3 * fused] +
             y * ParaTwist[idx][3 * fused + 1] +
             z * ParaTwist[idx][3 * fused + 2]) *
            (double)myEleNum[rsi];
      }
      PhysTwist[idx] += w * cexp(I * twoPi * twistWeight);
    }
#pragma omp for schedule(dynamic)
    for (idx = 0; idx < NNBodyG; idx++) {
      const int nbody = NBodyGN[idx];
      const int offset = NBodyGOffset[idx];
      int k;
      for (k = 0; k < nbody; k++) {
        scratch.rsi[k] = NBodyGIdx[offset + k][0] +
                         NBodyGIdx[offset + k][1] * Nsite;
        scratch.rsj[k] = NBodyGIdx[offset + k][2] +
                         NBodyGIdx[offset + k][3] * Nsite;
      }
      PhysNBodyG[idx] += w * GreenFuncNGC(
          nbody, scratch.rsi, scratch.rsj, ip, ncur, myEleIdx, myEleCfg,
          myEleNum, eleProjCnt, &scratch);
    }
#pragma omp master
    { StopTimer(54); }
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  ReleaseWorkSpaceThreadDouble();
}

#endif
