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
#include <stddef.h>
#include <stdio.h>

#include "include/global.h"
#include "include/locgrn_gc.h"
#include "include/matrix_gc.h"
#include "include/projection.h"
#include "include/qp.h"
#include "include/slater.h"
#include "include/vmccal.h"

#ifndef _SRC_VMCCAL_GC
#define _SRC_VMCCAL_GC

static void clearStoredOSampleRangeGC(const int sampleStart,
                                      const int sampleEnd) {
  const int sampleSize = sampleEnd - sampleStart;
  const size_t stride = 2 * (size_t)SROptSize;
  const size_t offset = (size_t)sampleStart * stride;
  const size_t count = (size_t)sampleSize * stride;
  size_t i;
  if (sampleSize <= 0) return;
  for (i = 0; i < count; i++) SROptO_Store[offset + i] = 0.0;
}

void VMCMainCalGC(MPI_Comm comm_parent, MPI_Comm comm) {
  int sampleStart;
  int sampleEnd;
  int sample;
  int rank;
  int size;
  double w = 0.0;
  double complex e = 0.0;
  (void)comm_parent;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  SplitLoop(&sampleStart, &sampleEnd, NVMCSample, rank, size);

  StartTimer(24);
  clearPhysQuantity();
  StopTimer(24);
  if (NVMCCalMode == 0 && NStoreO != 0) {
    clearStoredOSampleRangeGC(sampleStart, sampleEnd);
  }

  for (sample = sampleStart; sample < sampleEnd; sample++) {
    const int ncur = EleNumSample[sample];
    int *eleIdx = EleIdx + (size_t)sample * (size_t)NsizeMax;
    int *eleCfg = EleCfg + (size_t)sample * (size_t)Nsite2;
    int *eleNum = EleNum + (size_t)sample * (size_t)Nsite2;
    int *eleProjCnt = EleProjCnt + (size_t)sample * (size_t)NProj;
    double complex ip;
    double x;
    double sz;
    int info;

    StartTimer(40);
    info = CalculateMAllGC_fcmp(ncur, eleIdx, 0, NQPFull);
    StopTimer(40);
    if (info != GC_MALL_OK) {
      fprintf(stderr,
              "warning: VMCMainCalGC rank:%d sample:%d status:%d "
              "(CalculateMAllGC_fcmp)\n",
              rank, sample, info);
      continue;
    }
    ip = CalculateIP_fcmp(PfM, 0, NQPFull, MPI_COMM_SELF);
    x = LogProjVal(eleProjCnt);
    if (reweight == 1) {
      w = exp(2.0 * (log(cabs(ip)) + x) -
              logSqPfFullSlater[sample]);
    } else {
      w = 1.0;
    }
    if (!isfinite(w)) {
      fprintf(stderr, "warning: VMCMainCalGC rank:%d sample:%d w=%e\n",
              rank, sample, w);
      continue;
    }

    StartTimer(41);
    e = CalculateHamiltonianGC(ip, ncur, eleIdx, eleCfg, eleNum,
                               eleProjCnt);
    sz = CalculateSzGC(eleNum);
    StopTimer(41);
    if (!isfinite(creal(e)) || !isfinite(cimag(e))) {
      fprintf(stderr,
              "warning: VMCMainCalGC rank:%d sample:%d e=(%e,%e)\n",
              rank, sample, creal(e), cimag(e));
      continue;
    }

    Wc += w;
    Etot += w * e;
    Etot2 += w * conj(e) * e;
    Sztot += w * sz;
    Sztot2 += w * sz * sz;
    Ntot += w * (double)ncur;
    Ntot2 += w * (double)ncur * (double)ncur;

    if (NVMCCalMode == 0) {
      const int nProj = NProj;
      int i;
      SROptO[0] = 1.0;
      SROptO[1] = 0.0;
      for (i = 0; i < nProj; i++) {
        SROptO[2 * (i + 1)] = (double)eleProjCnt[i];
        SROptO[2 * (i + 1) + 1] = 0.0;
      }
      StartTimer(42);
      SlaterElmDiffGC_fcmp(SROptO + 2 * NProj + 2, ip, eleIdx, ncur);
      StopTimer(42);
      StartTimer(43);
      if (NStoreO == 0) {
        calculateOO(SROptOO, SROptHO, SROptO, w, e, SROptSize);
      } else {
        const double sqrtw = sqrt(w);
        const double complex we = w * e;
#pragma omp parallel for
        for (i = 0; i < 2 * SROptSize; i++) {
          SROptO_Store[i + (size_t)sample * 2 * (size_t)SROptSize] =
              sqrtw * SROptO[i];
          SROptHO[i] += we * SROptO[i];
        }
      }
      StopTimer(43);
    } else {
      StartTimer(42);
      CalculateGreenFuncGC(w, ip, ncur, eleIdx, eleCfg, eleNum,
                           eleProjCnt);
      StopTimer(42);
    }
  }

  if (NVMCCalMode == 0 && NStoreO != 0) {
    const int sampleSize = sampleEnd - sampleStart;
    if (sampleSize > 0) {
      double complex *store =
          SROptO_Store + (size_t)sampleStart * 2 * (size_t)SROptSize;
      StartTimer(45);
      calculateOO_Store(SROptOO, SROptHO, store, w, e,
                        2 * SROptSize, sampleSize);
      StopTimer(45);
    }
  }
}

#endif
