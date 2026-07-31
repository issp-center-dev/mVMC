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
#ifndef _LSLOCGRN_HEISENBERG
#define _LSLOCGRN_HEISENBERG

#include <complex.h>
#include <stddef.h>

typedef struct {
  int initialized;
  int complexMode;
  int nsize;
  int nsite2;
  int nproj;
  int nqpFull;
  int lapackLWork;
  int *pathEleNumA;
  int *pathEleNumBA;
  int *workEleIdx;
  int *workEleNum;
  int *projCntNew;
  double *realBuffer;
  double complex *complexBuffer;
  double *complexRWork;
  long long greenCalls[3];
#ifdef MVMC_ENABLE_FAULT_INJECTION
  int auditEnabled;
  long long auditedSamples;
#endif
} Lanczos2HeisenbergScratch;

/* The caller must zero-initialize scratch before its first initialization. */
int Lanczos2HeisenbergScratchInit(
    Lanczos2HeisenbergScratch *scratch, int complexMode);

void Lanczos2HeisenbergScratchFree(
    Lanczos2HeisenbergScratch *scratch);

int CalculateLS2HeisenbergLocalPowerReal(
    const double h1, const double ip, const int *eleIdx, const int *eleCfg,
    const int *eleNum, const int *eleProjCnt,
    Lanczos2HeisenbergScratch *scratch, double *localPower);

int CalculateLS2HeisenbergLocalPowerComplex(
    const double complex h1, const double complex ip, const int *eleIdx,
    const int *eleCfg, const int *eleNum, const int *eleProjCnt,
    const double complex *rbmCnt, Lanczos2HeisenbergScratch *scratch,
    double complex *localPower);

#endif
