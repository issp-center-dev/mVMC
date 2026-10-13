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
#include <stdlib.h>
#include <string.h>

#include "include/gc_antiparallel.h"
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

static void GCMeasurementAbort(const char *what, const int rank,
                               const int sample, MPI_Comm comm) {
  fprintf(stderr, "Error: GC measurement rank:%d sample:%d %s.\n", rank,
          sample, what);
  MPI_Abort(comm, EXIT_FAILURE);
}

#ifdef MVMC_ENABLE_FAULT_INJECTION
static int GCMeasureFaultActive(const char *kind, const int sample,
                                const int ncur) {
  const char *spec = getenv("MVMC_GC_FAULT_MEASURE");
  char name[32];
  int target = -1;
  int targetRank = -1;
  int worldRank = 0;
  if (spec == NULL || *spec == '\0') return 0;
  if (sscanf(spec, "%31[a-z_]:%d:rank=%d", name, &target, &targetRank) < 2)
    return 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  if (strcmp(name, kind) != 0 || target != sample ||
      (targetRank >= 0 && targetRank != worldRank))
    return 0;
  fprintf(stderr,
          "GC fault injected: kind=%s world_rank=%d sample=%d ncur=%d\n",
          kind, worldRank, sample, ncur);
  fflush(stderr);
  return 1;
}
#endif

void VMCMainCalGC(MPI_Comm comm_parent, MPI_Comm comm) {
  int sampleStart;
  int sampleEnd;
  int sample;
  int rank;
  int size;
  double w = 0.0;
  double complex e = 0.0;
  const int anti = GCAntiEnabled(FlagGrandCanonical, iFlgOrbitalGeneral);
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
#ifdef MVMC_ENABLE_FAULT_INJECTION
    const int injectRebuild =
        GCMeasureFaultActive("rebuild_nonfinite", sample, ncur);
#endif

    if (anti && !GCAntiValidateConfig(eleIdx, eleCfg, eleNum, Nsite, ncur)) {
      GCMeasurementAbort("is not a consistent Sz=0 configuration", rank,
                         sample, comm);
    }
#ifdef MVMC_ENABLE_FAULT_INJECTION
    if (injectRebuild && ncur >= 2) {
      SlaterElm[(size_t)eleIdx[0] * (size_t)Nsite2 +
                (size_t)eleIdx[1]] = NAN;
      SlaterElm[(size_t)eleIdx[1] * (size_t)Nsite2 +
                (size_t)eleIdx[0]] = NAN;
    }
#endif
    StartTimer(40);
    info = CalculateMAllGC_fcmp(ncur, eleIdx, 0, NQPFull);
    StopTimer(40);
    if (info != GC_MALL_OK) {
      GCMeasurementAbort("failed the matrix rebuild", rank, sample, comm);
    }
#ifdef MVMC_ENABLE_FAULT_INJECTION
    if (injectRebuild && ncur < 2) PfM[0] = NAN;
    if (GCMeasureFaultActive("zero_ip", sample, ncur)) {
      int qp;
      for (qp = 0; qp < NQPFull; qp++) PfM[qp] = 0.0;
    }
#endif
    ip = CalculateIP_fcmp(PfM, 0, NQPFull, MPI_COMM_SELF);
    if (!isfinite(creal(ip)) || !isfinite(cimag(ip)) || cabs(ip) == 0.0) {
      GCMeasurementAbort("has a zero or nonfinite amplitude", rank, sample,
                         comm);
    }
    x = LogProjVal(eleProjCnt);
    if (reweight == 1) {
      w = exp(2.0 * (log(cabs(ip)) + x) -
              logSqPfFullSlater[sample]);
    } else {
      w = 1.0;
    }
    if (!isfinite(w)) {
      GCMeasurementAbort("has a nonfinite reweight factor", rank, sample,
                         comm);
    }

    StartTimer(41);
    e = CalculateHamiltonianGC(ip, ncur, eleIdx, eleCfg, eleNum,
                               eleProjCnt);
    sz = CalculateSzGC(eleNum);
    StopTimer(41);
    if (!isfinite(creal(e)) || !isfinite(cimag(e))) {
      GCMeasurementAbort("has a nonfinite local energy", rank, sample,
                         comm);
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
