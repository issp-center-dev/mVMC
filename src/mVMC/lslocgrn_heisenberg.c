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
#include "lslocgrn_heisenberg.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"

static void *Lanczos2HeisenbergCheckedCalloc(size_t count, size_t size) {
  if (count == 0) count = 1;
  if (size == 0 || count > SIZE_MAX / size) return NULL;
  return calloc(count, size);
}

void Lanczos2HeisenbergScratchFree(
    Lanczos2HeisenbergScratch *scratch) {
  if (scratch == NULL) return;
  free(scratch->pathEleNumA);
  free(scratch->pathEleNumBA);
  free(scratch->workEleIdx);
  free(scratch->workEleNum);
  free(scratch->projCntNew);
  free(scratch->realBuffer);
  free(scratch->complexBuffer);
  free(scratch->complexRWork);
  memset(scratch, 0, sizeof(*scratch));
}

int Lanczos2HeisenbergScratchInit(
    Lanczos2HeisenbergScratch *scratch, int complexMode) {
  size_t bufferCount;

  if (scratch == NULL || (complexMode != 0 && complexMode != 1) ||
      Nsize <= 0 || Nsite2 <= 0 || NProj < 0 || NQPFull != 1 ||
      LapackLWork <= 0 ||
      (size_t)Nsize > (SIZE_MAX - (size_t)NQPFull) / 6) {
    return -1;
  }
  if (scratch->initialized != 0) return -1;
  memset(scratch, 0, sizeof(*scratch));
  scratch->complexMode = complexMode;
  scratch->nsize = Nsize;
  scratch->nsite2 = Nsite2;
  scratch->nproj = NProj;
  scratch->nqpFull = NQPFull;
  scratch->lapackLWork = LapackLWork;
  bufferCount = (size_t)NQPFull + 6 * (size_t)Nsize;

  scratch->pathEleNumA = Lanczos2HeisenbergCheckedCalloc(
      (size_t)Nsite2, sizeof(*scratch->pathEleNumA));
  scratch->pathEleNumBA = Lanczos2HeisenbergCheckedCalloc(
      (size_t)Nsite2, sizeof(*scratch->pathEleNumBA));
  scratch->workEleIdx = Lanczos2HeisenbergCheckedCalloc(
      (size_t)Nsize, sizeof(*scratch->workEleIdx));
  scratch->workEleNum = Lanczos2HeisenbergCheckedCalloc(
      (size_t)Nsite2, sizeof(*scratch->workEleNum));
  scratch->projCntNew = Lanczos2HeisenbergCheckedCalloc(
      (size_t)(NProj > 0 ? NProj : 1), sizeof(*scratch->projCntNew));
  if (complexMode == 0) {
    scratch->realBuffer = Lanczos2HeisenbergCheckedCalloc(
        bufferCount, sizeof(*scratch->realBuffer));
  } else {
    scratch->complexBuffer = Lanczos2HeisenbergCheckedCalloc(
        bufferCount, sizeof(*scratch->complexBuffer));
    scratch->complexRWork = Lanczos2HeisenbergCheckedCalloc(
        (size_t)LapackLWork, sizeof(*scratch->complexRWork));
  }

  if (scratch->pathEleNumA == NULL || scratch->pathEleNumBA == NULL ||
      scratch->workEleIdx == NULL || scratch->workEleNum == NULL ||
      scratch->projCntNew == NULL ||
      (complexMode == 0 && scratch->realBuffer == NULL) ||
      (complexMode != 0 &&
       (scratch->complexBuffer == NULL || scratch->complexRWork == NULL))) {
    Lanczos2HeisenbergScratchFree(scratch);
    return -1;
  }
#ifdef MVMC_ENABLE_FAULT_INJECTION
  {
    const char *auditValue = getenv("MVMC_LANCZOS2_HEISENBERG_AUDIT");
    scratch->auditEnabled =
        auditValue != NULL && auditValue[0] != '\0' &&
        strcmp(auditValue, "0") != 0;
  }
#endif
  scratch->initialized = 1;
  return 0;
}

static int Lanczos2HeisenbergScratchMatches(
    const Lanczos2HeisenbergScratch *scratch, int complexMode) {
  return scratch != NULL && scratch->initialized != 0 &&
         scratch->complexMode == complexMode &&
         scratch->nsize == Nsize && scratch->nsite2 == Nsite2 &&
         scratch->nproj == NProj && scratch->nqpFull == NQPFull &&
         scratch->lapackLWork == LapackLWork;
}

/*
 * Build one oriented Exchange operator and its resulting occupation.
 * orientation=0 is c^+_{i up} c_{j up} c^+_{j down} c_{i down};
 * orientation=1 is its spin-reversed partner.
 */
static int Lanczos2BuildExchangeOperator(
    int bond, int orientation, const int *sourceEleNum, int *targetEleNum,
    int *rsi, int *rsj) {
  int ri;
  int rj;
  int riUp;
  int rjUp;
  int riDown;
  int rjDown;
  int source0;
  int target0;
  int source1;
  int target1;

  if (bond < 0 || bond >= NExchangeCoupling ||
      (orientation != 0 && orientation != 1) ||
      sourceEleNum == NULL || rsi == NULL || rsj == NULL) {
    return 0;
  }
  ri = ExchangeCoupling[bond][0];
  rj = ExchangeCoupling[bond][1];
  if (ri < 0 || ri >= Nsite || rj < 0 || rj >= Nsite) return 0;
  riUp = ri;
  rjUp = rj;
  riDown = ri + Nsite;
  rjDown = rj + Nsite;
  if (orientation == 0) {
    source0 = rjUp;
    target0 = riUp;
    source1 = riDown;
    target1 = rjDown;
  } else {
    source0 = rjDown;
    target0 = riDown;
    source1 = riUp;
    target1 = rjUp;
  }
  if (sourceEleNum[source0] == 0 || sourceEleNum[target0] != 0 ||
      sourceEleNum[source1] == 0 || sourceEleNum[target1] != 0) {
    return 0;
  }

  if (targetEleNum != NULL) {
    memcpy(targetEleNum, sourceEleNum, (size_t)Nsite2 * sizeof(int));
    targetEleNum[source0] = 0;
    targetEleNum[target0] = 1;
    targetEleNum[source1] = 0;
    targetEleNum[target1] = 1;
  }
  rsi[0] = target0;
  rsj[0] = source0;
  rsi[1] = target1;
  rsj[1] = source1;
  return 1;
}

static void Lanczos2ComposeOperators(
    int depth, const int *leftRsi, const int *leftRsj,
    const int *middleRsi, const int *middleRsj,
    const int *rightRsi, const int *rightRsj,
    int *rsi, int *rsj) {
  int offset = 0;
  int index;
  const int *sourcesRsi[3] = {leftRsi, middleRsi, rightRsi};
  const int *sourcesRsj[3] = {leftRsj, middleRsj, rightRsj};
  for (index = 0; index < depth; index++) {
    rsi[offset] = sourcesRsi[index][0];
    rsj[offset++] = sourcesRsj[index][0];
    rsi[offset] = sourcesRsi[index][1];
    rsj[offset++] = sourcesRsj[index][1];
  }
}

static double Lanczos2HeisenbergGreenReal(
    int n, const int *rsi, const int *rsj, double ip, const int *eleIdx,
    const int *eleCfg, const int *eleNum, const int *eleProjCnt,
    Lanczos2HeisenbergScratch *scratch) {
  int callRsi[6];
  int callRsj[6];
  memcpy(callRsi, rsi, (size_t)n * sizeof(int));
  memcpy(callRsj, rsj, (size_t)n * sizeof(int));
  memcpy(scratch->workEleIdx, eleIdx, (size_t)Nsize * sizeof(int));
  memcpy(scratch->workEleNum, eleNum, (size_t)Nsite2 * sizeof(int));
  scratch->greenCalls[n / 2 - 1]++;
  return GreenFuncN_real(
      n, callRsi, callRsj, ip, scratch->workEleIdx, eleCfg,
      scratch->workEleNum, eleProjCnt, scratch->realBuffer,
      scratch->projCntNew);
}

static double complex Lanczos2HeisenbergGreenComplex(
    int n, const int *rsi, const int *rsj, double complex ip,
    const int *eleIdx, const int *eleCfg, const int *eleNum,
    const int *eleProjCnt, const double complex *rbmCnt,
    Lanczos2HeisenbergScratch *scratch) {
  int callRsi[6];
  int callRsj[6];
  memcpy(callRsi, rsi, (size_t)n * sizeof(int));
  memcpy(callRsj, rsj, (size_t)n * sizeof(int));
  memcpy(scratch->workEleIdx, eleIdx, (size_t)Nsize * sizeof(int));
  memcpy(scratch->workEleNum, eleNum, (size_t)Nsite2 * sizeof(int));
  scratch->greenCalls[n / 2 - 1]++;
  return GreenFuncN(
      n, callRsi, callRsj, ip, scratch->workEleIdx, eleCfg,
      scratch->workEleNum, eleProjCnt, rbmCnt, NULL,
      scratch->complexBuffer, scratch->projCntNew, scratch->complexRWork);
}

#ifdef MVMC_ENABLE_FAULT_INJECTION
static int Lanczos2ScaledCloseReal(double left, double right) {
  const double scale = fmax(fmax(fabs(left), fabs(right)), 1.0);
  return fabs(left - right) <= 2.0e-10 + 2.0e-9 * scale;
}

static int Lanczos2ScaledCloseComplex(
    double complex left, double complex right) {
  const double scale = fmax(fmax(cabs(left), cabs(right)), 1.0);
  return cabs(left - right) <= 2.0e-10 + 2.0e-9 * scale;
}
#endif

int CalculateLS2HeisenbergLocalPowerReal(
    const double h1, const double ip, const int *eleIdx, const int *eleCfg,
    const int *eleNum, const int *eleProjCnt,
    Lanczos2HeisenbergScratch *scratch, double *localPower) {
  double legacyLocalPower[4];
  double directF1;
  double directF2;
  double directF3;
  double v0;
  int aBond;
  int aOrientation;

  if (!Lanczos2HeisenbergScratchMatches(scratch, 0) ||
      eleIdx == NULL || eleCfg == NULL || eleNum == NULL ||
      localPower == NULL || !isfinite(ip) || fabs(ip) == 0.0) {
    return -1;
  }
  v0 = CalculateHamiltonian0_real(eleNum);
  LSLocalQ_real(h1, ip, (int *)eleIdx, (int *)eleCfg, (int *)eleNum,
                (int *)eleProjCnt, legacyLocalPower);
  localPower[0] = 1.0;
  localPower[1] = h1;
  localPower[2] = legacyLocalPower[3];
  directF1 = v0;
  directF2 = v0 * v0;
  directF3 = v0 * v0 * v0;

  for (aBond = 0; aBond < NExchangeCoupling; aBond++) {
    for (aOrientation = 0; aOrientation < 2; aOrientation++) {
      int aRsi[2];
      int aRsj[2];
      int composedRsi[6];
      int composedRsj[6];
      const double ja = ParaExchangeCoupling[aBond];
      double ga;
      double va;
      int bBond;
      int bOrientation;

      if (!Lanczos2BuildExchangeOperator(
              aBond, aOrientation, eleNum, scratch->pathEleNumA,
              aRsi, aRsj)) {
        continue;
      }
      va = CalculateHamiltonian0_real(scratch->pathEleNumA);
      Lanczos2ComposeOperators(
          1, aRsi, aRsj, NULL, NULL, NULL, NULL,
          composedRsi, composedRsj);
      ga = Lanczos2HeisenbergGreenReal(
          2, composedRsi, composedRsj, ip, eleIdx, eleCfg, eleNum,
          eleProjCnt, scratch);
      directF1 += ja * ga;
      directF2 += ja * (v0 + va) * ga;
      directF3 += ja * (v0 * v0 + v0 * va + va * va) * ga;

      for (bBond = 0; bBond < NExchangeCoupling; bBond++) {
        for (bOrientation = 0; bOrientation < 2; bOrientation++) {
          int bRsi[2];
          int bRsj[2];
          const double jb = ParaExchangeCoupling[bBond];
          double gba;
          double vba;
          int cBond;
          int cOrientation;

          if (!Lanczos2BuildExchangeOperator(
                  bBond, bOrientation, scratch->pathEleNumA,
                  scratch->pathEleNumBA, bRsi, bRsj)) {
            continue;
          }
          vba = CalculateHamiltonian0_real(scratch->pathEleNumBA);
          Lanczos2ComposeOperators(
              2, bRsi, bRsj, aRsi, aRsj, NULL, NULL,
              composedRsi, composedRsj);
          gba = Lanczos2HeisenbergGreenReal(
              4, composedRsi, composedRsj, ip, eleIdx, eleCfg, eleNum,
              eleProjCnt, scratch);
          directF2 += jb * ja * gba;
          directF3 += jb * ja * (v0 + va + vba) * gba;

          for (cBond = 0; cBond < NExchangeCoupling; cBond++) {
            for (cOrientation = 0; cOrientation < 2; cOrientation++) {
              int cRsi[2];
              int cRsj[2];
              const double jc = ParaExchangeCoupling[cBond];
              double gcba;

              if (!Lanczos2BuildExchangeOperator(
                      cBond, cOrientation, scratch->pathEleNumBA,
                      NULL, cRsi, cRsj)) {
                continue;
              }
              Lanczos2ComposeOperators(
                  3, cRsi, cRsj, bRsi, bRsj, aRsi, aRsj,
                  composedRsi, composedRsj);
              gcba = Lanczos2HeisenbergGreenReal(
                  6, composedRsi, composedRsj, ip, eleIdx, eleCfg,
                  eleNum, eleProjCnt, scratch);
              directF3 += jc * jb * ja * gcba;
            }
          }
        }
      }
    }
  }
  localPower[3] = directF3;
#ifdef MVMC_ENABLE_FAULT_INJECTION
  if (scratch->auditEnabled) {
    if (!Lanczos2ScaledCloseReal(h1, directF1) ||
        !Lanczos2ScaledCloseReal(localPower[2], directF2)) {
      return -1;
    }
    scratch->auditedSamples++;
  }
#endif
  return 0;
}

int CalculateLS2HeisenbergLocalPowerComplex(
    const double complex h1, const double complex ip, const int *eleIdx,
    const int *eleCfg, const int *eleNum, const int *eleProjCnt,
    const double complex *rbmCnt, Lanczos2HeisenbergScratch *scratch,
    double complex *localPower) {
  double complex legacyLocalPower[4];
  double complex directF1;
  double complex directF2;
  double complex directF3;
  double complex v0;
  int aBond;
  int aOrientation;

  if (!Lanczos2HeisenbergScratchMatches(scratch, 1) ||
      eleIdx == NULL || eleCfg == NULL || eleNum == NULL ||
      localPower == NULL ||
      !isfinite(creal(ip)) || !isfinite(cimag(ip)) || cabs(ip) == 0.0) {
    return -1;
  }
  v0 = CalculateHamiltonian0(eleNum);
  LSLocalQ(h1, ip, (int *)eleIdx, (int *)eleCfg, (int *)eleNum,
           (int *)eleProjCnt, (double complex *)rbmCnt, legacyLocalPower);
  localPower[0] = 1.0;
  localPower[1] = h1;
  localPower[2] = legacyLocalPower[3];
  directF1 = v0;
  directF2 = v0 * v0;
  directF3 = v0 * v0 * v0;

  for (aBond = 0; aBond < NExchangeCoupling; aBond++) {
    for (aOrientation = 0; aOrientation < 2; aOrientation++) {
      int aRsi[2];
      int aRsj[2];
      int composedRsi[6];
      int composedRsj[6];
      const double ja = ParaExchangeCoupling[aBond];
      double complex ga;
      double complex va;
      int bBond;
      int bOrientation;

      if (!Lanczos2BuildExchangeOperator(
              aBond, aOrientation, eleNum, scratch->pathEleNumA,
              aRsi, aRsj)) {
        continue;
      }
      va = CalculateHamiltonian0(scratch->pathEleNumA);
      Lanczos2ComposeOperators(
          1, aRsi, aRsj, NULL, NULL, NULL, NULL,
          composedRsi, composedRsj);
      ga = Lanczos2HeisenbergGreenComplex(
          2, composedRsi, composedRsj, ip, eleIdx, eleCfg, eleNum,
          eleProjCnt, rbmCnt, scratch);
      directF1 += ja * ga;
      directF2 += ja * (v0 + va) * ga;
      directF3 += ja * (v0 * v0 + v0 * va + va * va) * ga;

      for (bBond = 0; bBond < NExchangeCoupling; bBond++) {
        for (bOrientation = 0; bOrientation < 2; bOrientation++) {
          int bRsi[2];
          int bRsj[2];
          const double jb = ParaExchangeCoupling[bBond];
          double complex gba;
          double complex vba;
          int cBond;
          int cOrientation;

          if (!Lanczos2BuildExchangeOperator(
                  bBond, bOrientation, scratch->pathEleNumA,
                  scratch->pathEleNumBA, bRsi, bRsj)) {
            continue;
          }
          vba = CalculateHamiltonian0(scratch->pathEleNumBA);
          Lanczos2ComposeOperators(
              2, bRsi, bRsj, aRsi, aRsj, NULL, NULL,
              composedRsi, composedRsj);
          gba = Lanczos2HeisenbergGreenComplex(
              4, composedRsi, composedRsj, ip, eleIdx, eleCfg, eleNum,
              eleProjCnt, rbmCnt, scratch);
          directF2 += jb * ja * gba;
          directF3 += jb * ja * (v0 + va + vba) * gba;

          for (cBond = 0; cBond < NExchangeCoupling; cBond++) {
            for (cOrientation = 0; cOrientation < 2; cOrientation++) {
              int cRsi[2];
              int cRsj[2];
              const double jc = ParaExchangeCoupling[cBond];
              double complex gcba;

              if (!Lanczos2BuildExchangeOperator(
                      cBond, cOrientation, scratch->pathEleNumBA,
                      NULL, cRsi, cRsj)) {
                continue;
              }
              Lanczos2ComposeOperators(
                  3, cRsi, cRsj, bRsi, bRsj, aRsi, aRsj,
                  composedRsi, composedRsj);
              gcba = Lanczos2HeisenbergGreenComplex(
                  6, composedRsi, composedRsj, ip, eleIdx, eleCfg,
                  eleNum, eleProjCnt, rbmCnt, scratch);
              directF3 += jc * jb * ja * gcba;
            }
          }
        }
      }
    }
  }
  localPower[3] = directF3;
#ifdef MVMC_ENABLE_FAULT_INJECTION
  if (scratch->auditEnabled) {
    if (!Lanczos2ScaledCloseComplex(h1, directF1) ||
        !Lanczos2ScaledCloseComplex(localPower[2], directF2)) {
      return -1;
    }
    scratch->auditedSamples++;
  }
#endif
  return 0;
}
