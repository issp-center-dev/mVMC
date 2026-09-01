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
#include "include/slater.h"
#include "include/workspace.h"

#ifndef _SRC_SLATER_GC
#define _SRC_SLATER_GC

void SlaterElmDiffGC_fcmp(double complex *srOptO, const double complex ip,
                          const int *eleIdx, const int ncur) {
  const size_t bufferCount =
      GCCheckedMulSize((size_t)NQPFull, (size_t)NSlater);
  const int bufferWorkCount =
      GCCheckedSizeToInt(bufferCount, "GC Slater derivative workspace");
  const size_t inverseStride =
      GCCheckedMulSize((size_t)NsizeMax, (size_t)NsizeMax);
  const double complex invIP = 1.0 / ip;
  double complex *buffer;
  int orbital;
  int qpidx;

  RequestWorkSpaceComplex(bufferWorkCount);
  buffer = GetWorkSpaceComplex(bufferWorkCount);
  for (orbital = 0; orbital < bufferWorkCount; orbital++) buffer[orbital] = 0.0;

  for (qpidx = 0; qpidx < NQPFull; qpidx++) {
    const int optidx = qpidx / NQPFix;
    const int mpidx = (qpidx % NQPFix) / NSPGaussLeg;
    const int *xqpOpt = QPOptTrans[optidx];
    const int *xqpOptSgn = QPOptTransSgn[optidx];
    const int *xqp = QPTrans[mpidx];
    const int *xqpSgn = QPTransSgn[mpidx];
    const double complex pfaffian = PfM[qpidx];
    const double complex *inverse = InvM + (size_t)qpidx * inverseStride;
    double complex *qpBuffer = buffer + (size_t)qpidx * (size_t)NSlater;
    int msi;

    for (msi = 0; msi < ncur; msi++) {
      const int rsi = eleIdx[msi];
      const int ri = rsi % Nsite;
      const int si = rsi / Nsite;
      const int ori = xqpOpt[ri];
      const int tri = xqp[ori] + si * Nsite;
      const int sgni = xqpSgn[ori] * xqpOptSgn[ri];
      int msj;
      for (msj = 0; msj < ncur; msj++) {
        const int rsj = eleIdx[msj];
        const int rj = rsj % Nsite;
        const int sj = rsj / Nsite;
        const int orj = xqpOpt[rj];
        const int trj = xqp[orj] + sj * Nsite;
        const int sgnj = xqpSgn[orj] * xqpOptSgn[rj];
        const int orbitalIndex = OrbitalIdx[tri][trj];
        /* matrix_gc stores (-SlaterElm)^-1, while the legacy fixed-Sz
         * derivative receives the opposite inverse convention. */
        qpBuffer[orbitalIndex] +=
            inverse[(size_t)msi * (size_t)NsizeMax + (size_t)msj] *
            pfaffian * (double)(sgni * sgnj * OrbitalSgn[tri][trj]);
      }
    }
  }

  for (orbital = 0; orbital < NSlater; orbital++) {
    srOptO[2 * orbital] = 0.0;
    srOptO[2 * orbital + 1] = 0.0;
  }
  for (qpidx = 0; qpidx < NQPFull; qpidx++) {
    const double complex weight = QPFullWeight[qpidx];
    const double complex *qpBuffer =
        buffer + (size_t)qpidx * (size_t)NSlater;
    for (orbital = 0; orbital < NSlater; orbital++) {
      srOptO[2 * orbital] += weight * qpBuffer[orbital];
      srOptO[2 * orbital + 1] += weight * qpBuffer[orbital] * I;
    }
  }
  for (orbital = 0; orbital < NSlater; orbital++) {
    srOptO[2 * orbital] *= invIP;
    srOptO[2 * orbital + 1] *= invIP;
  }
  ReleaseWorkSpaceComplex();
}

#endif
