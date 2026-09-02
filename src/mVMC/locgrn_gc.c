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
#include "include/pfupdate_gc.h"
#include "include/projection.h"
#include "include/qp.h"

#ifndef _SRC_LOCGRN_GC
#define _SRC_LOCGRN_GC

static double complex GCProjectedCandidate(const double complex *pfMNew) {
  return CalculateIP_fcmp((double complex *)pfMNew, 0, NQPFull,
                          MPI_COMM_SELF);
}

double complex GreenFunc1GC(const int rsi, const int rsj,
                            const double complex ip, const int ncur,
                            int *eleIdx, int *eleCfg, int *eleNum,
                            const int *eleProjCnt,
                            GCGreenScratch *scratch) {
  int particle;
  double projectionRatio;
  double complex candidateOverlap;

  if (rsi == rsj) return (double)eleNum[rsi];
  if (eleNum[rsi] != 0 || eleNum[rsj] == 0) return 0.0;
  particle = eleCfg[rsj];

  eleIdx[particle] = rsi;
  eleCfg[rsj] = -1;
  eleCfg[rsi] = particle;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt_fsz(rsj % Nsite, rsi % Nsite, rsj / Nsite, rsi / Nsite,
                    scratch->projCntNew, eleProjCnt, eleNum);
  projectionRatio = ProjRatio(scratch->projCntNew, eleProjCnt);
  CalculateNewPfMHopGC(particle, scratch->pfMNew, eleIdx, ncur, 0,
                       NQPFull);
  candidateOverlap = GCProjectedCandidate(scratch->pfMNew);

  eleIdx[particle] = rsj;
  eleCfg[rsi] = -1;
  eleCfg[rsj] = particle;
  eleNum[rsi] = 0;
  eleNum[rsj] = 1;
  return conj(projectionRatio * candidateOverlap / ip);
}

double complex GreenFunc2GC(const int rsi, const int rsj, const int rsk,
                            const int rsl, const double complex ip,
                            const int ncur, int *eleIdx, int *eleCfg,
                            int *eleNum, const int *eleProjCnt,
                            GCGreenScratch *scratch) {
  int particleJ;
  int particleL;
  double projectionRatio;
  double complex candidateOverlap;

  if (rsi == rsj) {
    if (rsj == rsk) {
      if (rsk == rsl) return (double)eleNum[rsi];
      if (eleNum[rsi] != 0) return 0.0;
      return GreenFunc1GC(rsk, rsl, ip, ncur, eleIdx, eleCfg, eleNum,
                          eleProjCnt, scratch);
    }
    if (rsj == rsl) return 0.0;
    if (rsk == rsl) return (double)(eleNum[rsi] * eleNum[rsk]);
    if (eleNum[rsi] == 0) return 0.0;
    return GreenFunc1GC(rsk, rsl, ip, ncur, eleIdx, eleCfg, eleNum,
                        eleProjCnt, scratch);
  }
  if (rsi == rsk) return 0.0;
  if (rsi == rsl) {
    if (rsj == rsk) return (double)(eleNum[rsi] * (1 - eleNum[rsj]));
    if (eleNum[rsi] == 0) return 0.0;
    return -GreenFunc1GC(rsk, rsj, ip, ncur, eleIdx, eleCfg, eleNum,
                         eleProjCnt, scratch);
  }
  if (rsj == rsk) {
    if (rsk == rsl) {
      if (eleNum[rsj] == 0) return 0.0;
      return GreenFunc1GC(rsi, rsj, ip, ncur, eleIdx, eleCfg, eleNum,
                          eleProjCnt, scratch);
    }
    if (eleNum[rsj] != 0) return 0.0;
    return GreenFunc1GC(rsi, rsl, ip, ncur, eleIdx, eleCfg, eleNum,
                        eleProjCnt, scratch);
  }
  if (rsj == rsl) return 0.0;
  if (rsk == rsl) {
    if (eleNum[rsk] == 0) return 0.0;
    return GreenFunc1GC(rsi, rsj, ip, ncur, eleIdx, eleCfg, eleNum,
                        eleProjCnt, scratch);
  }
  if (eleNum[rsi] != 0 || eleNum[rsj] == 0 || eleNum[rsk] != 0 ||
      eleNum[rsl] == 0) {
    return 0.0;
  }

  particleJ = eleCfg[rsj];
  particleL = eleCfg[rsl];
  eleIdx[particleL] = rsk;
  eleCfg[rsl] = -1;
  eleCfg[rsk] = particleL;
  eleNum[rsl] = 0;
  eleNum[rsk] = 1;
  UpdateProjCnt_fsz(rsl % Nsite, rsk % Nsite, rsl / Nsite, rsk / Nsite,
                    scratch->projCntNew, eleProjCnt, eleNum);
  eleIdx[particleJ] = rsi;
  eleCfg[rsj] = -1;
  eleCfg[rsi] = particleJ;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt_fsz(rsj % Nsite, rsi % Nsite, rsj / Nsite, rsi / Nsite,
                    scratch->projCntNew, scratch->projCntNew, eleNum);
  projectionRatio = ProjRatio(scratch->projCntNew, eleProjCnt);
  CalculateNewPfMTwoHopGCWorkspace(
      particleL, particleJ, scratch->pfMNew, eleIdx, ncur, 0, NQPFull,
      scratch->vec, scratch->vec + NsizeMax);
  candidateOverlap = GCProjectedCandidate(scratch->pfMNew);

  eleIdx[particleJ] = rsj;
  eleCfg[rsi] = -1;
  eleCfg[rsj] = particleJ;
  eleNum[rsi] = 0;
  eleNum[rsj] = 1;
  eleIdx[particleL] = rsl;
  eleCfg[rsk] = -1;
  eleCfg[rsl] = particleL;
  eleNum[rsk] = 0;
  eleNum[rsl] = 1;
  return conj(projectionRatio * candidateOverlap / ip);
}

double complex GreenFuncNGC(const int n, int *rsi, int *rsj,
                            const double complex ip, const int ncur,
                            int *eleIdx, int *eleCfg, int *eleNum,
                            const int *eleProjCnt,
                            GCGreenScratch *scratch) {
  int k;
  int l;
  int m;
  double projectionRatio;
  double complex candidateOverlap;

  if (n <= 0 || n > scratch->maxN) return 0.0;
  if (n == 1) {
    return GreenFunc1GC(rsi[0], rsj[0], ip, ncur, eleIdx, eleCfg, eleNum,
                        eleProjCnt, scratch);
  }
  if (n == 2) {
    return GreenFunc2GC(rsi[0], rsj[0], rsi[1], rsj[1], ip, ncur,
                        eleIdx, eleCfg, eleNum, eleProjCnt, scratch);
  }

  for (k = n - 1; k >= 0; k--) {
    const int annihilator = rsj[k];
    for (l = k + 1; l < n; l++) {
      if (annihilator == rsi[l]) {
        rsj[k] = rsj[l];
        for (m = l; m < n - 1; m++) {
          rsi[m] = rsi[m + 1];
          rsj[m] = rsj[m + 1];
        }
        return GreenFuncNGC(n - 1, rsi, rsj, ip, ncur, eleIdx, eleCfg,
                            eleNum, eleProjCnt, scratch);
      }
      if (annihilator == rsj[l]) return 0.0;
    }
    if (eleNum[annihilator] == 0) return 0.0;

    if (rsi[k] == rsj[k]) {
      for (m = k; m < n - 1; m++) {
        rsi[m] = rsi[m + 1];
        rsj[m] = rsj[m + 1];
      }
      return GreenFuncNGC(n - 1, rsi, rsj, ip, ncur, eleIdx, eleCfg,
                          eleNum, eleProjCnt, scratch);
    }
    for (l = k + 1; l < n; l++) {
      if (rsi[k] == rsi[l]) return 0.0;
      if (rsi[k] == rsj[l]) {
        rsi[k] = rsi[l];
        for (m = l; m < n - 1; m++) {
          rsi[m] = rsi[m + 1];
          rsj[m] = rsj[m + 1];
        }
        return -GreenFuncNGC(n - 1, rsi, rsj, ip, ncur, eleIdx,
                             eleCfg, eleNum, eleProjCnt, scratch);
      }
    }
    if (eleNum[rsi[k]] != 0) return 0.0;
  }

  for (k = 0; k < n; k++) {
    const int source = rsj[k];
    const int destination = rsi[k];
    const int particle = eleCfg[source];
    scratch->msa[k] = particle;
    eleIdx[particle] = destination;
    eleCfg[source] = -1;
    eleCfg[destination] = particle;
    eleNum[source] = 0;
    eleNum[destination] = 1;
    UpdateProjCnt_fsz(source % Nsite, destination % Nsite, source / Nsite,
                      destination / Nsite, scratch->projCntNew,
                      k == 0 ? eleProjCnt : scratch->projCntNew, eleNum);
  }
  projectionRatio = ProjRatio(scratch->projCntNew, eleProjCnt);
  for (k = 0; k < NQPFull; k++) {
    scratch->pfMNew[k] = CalculateNewPfMNGC(
        k, n, scratch->msa, rsi, eleIdx, ncur, scratch->vec, scratch->w,
        scratch->smallMat, scratch->pfIWork, scratch->pfWork,
        scratch->pfLWork, scratch->pfRWork);
  }
  candidateOverlap = GCProjectedCandidate(scratch->pfMNew);

  for (k = n - 1; k >= 0; k--) {
    const int source = rsj[k];
    const int destination = rsi[k];
    const int particle = scratch->msa[k];
    eleIdx[particle] = source;
    eleCfg[destination] = -1;
    eleCfg[source] = particle;
    eleNum[destination] = 0;
    eleNum[source] = 1;
  }
  return conj(projectionRatio * candidateOverlap / ip);
}

double complex GreenFuncPairAddGC(
    const int rsa, const int rsb, const double complex ip,
    const int ncur, const int *eleIdx, const int *eleCfg, int *eleNum,
    const int *eleProjCnt, GCGreenScratch *scratch) {
  double projectionRatio;
  double complex candidateOverlap;
  (void)eleCfg;
  if (rsa == rsb) return 0.0;
  if (eleNum[rsa] != 0 || eleNum[rsb] != 0) return 0.0;
  if (ncur + 2 > NsizeMax) return 0.0;

  eleNum[rsa] = 1;
  eleNum[rsb] = 1;
  MakeProjCnt(scratch->projCntNew, eleNum);
  projectionRatio = ProjRatio(scratch->projCntNew, eleProjCnt);
  eleNum[rsa] = 0;
  eleNum[rsb] = 0;

  CalculateNewPfMAddGCWorkspace(rsa, rsb, scratch->pfMNew, eleIdx, ncur,
                                0, NQPFull, scratch->vec,
                                scratch->vec + NsizeMax);
  candidateOverlap = GCProjectedCandidate(scratch->pfMNew);
  return conj(projectionRatio * candidateOverlap / ip);
}

double complex GreenFuncPairRemoveGC(
    const int rsi, const int rsj, const double complex ip,
    const int ncur, const int *eleIdx, const int *eleCfg, int *eleNum,
    const int *eleProjCnt, GCGreenScratch *scratch) {
  const size_t invStride =
      GCCheckedMulSize((size_t)NsizeMax, (size_t)NsizeMax);
  double projectionRatio;
  double complex candidateOverlap;
  double sign;
  int pi;
  int pj;
  int p0;
  int p1;
  int qpidx;
  (void)eleIdx;
  if (rsi == rsj) return 0.0;
  if (eleNum[rsi] == 0 || eleNum[rsj] == 0) return 0.0;
  if (ncur < 2) return 0.0;
  pi = eleCfg[rsi];
  pj = eleCfg[rsj];
  sign = (pi < pj) ? 1.0 : -1.0;
  p0 = (pi < pj) ? pi : pj;
  p1 = (pi < pj) ? pj : pi;

  eleNum[rsi] = 0;
  eleNum[rsj] = 0;
  MakeProjCnt(scratch->projCntNew, eleNum);
  projectionRatio = ProjRatio(scratch->projCntNew, eleProjCnt);
  eleNum[rsi] = 1;
  eleNum[rsj] = 1;

  for (qpidx = 0; qpidx < NQPFull; qpidx++) {
    const double complex *invM = InvM + (size_t)qpidx * invStride;
    scratch->pfMNew[qpidx] =
        sign * invM[(size_t)p0 * (size_t)NsizeMax + (size_t)p1] *
        PfM[qpidx];
  }
  candidateOverlap = GCProjectedCandidate(scratch->pfMNew);
  return conj(projectionRatio * candidateOverlap / ip);
}

#endif
