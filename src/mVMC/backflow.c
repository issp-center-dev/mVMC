/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/
/*-------------------------------------------------------------
 * Backflow plumbing and validation
 *-------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "./include/backflow.h"
#include "./include/global.h"

static int BFCheckedMulLL(long long a, long long b, const char *name, long long *out) {
  if (a < 0 || b < 0) {
    fprintf(stderr, "Error: negative size while computing %s.\n", name);
    return 1;
  }
  if (a != 0 && b > LLONG_MAX / a) {
    fprintf(stderr, "Error: integer overflow while computing %s.\n", name);
    return 1;
  }
  *out = a * b;
  return 0;
}

static int BFCheckedAddLL(long long a, long long b, const char *name, long long *out) {
  if (a < 0 || b < 0 || b > LLONG_MAX - a) {
    fprintf(stderr, "Error: integer overflow while computing %s.\n", name);
    return 1;
  }
  *out = a + b;
  return 0;
}

static int BFIntFromLL(long long value, const char *name, int *out) {
  if (value < 0 || value > INT_MAX) {
    fprintf(stderr, "Error: %s is out of int range (%lld).\n", name, value);
    return 1;
  }
  *out = (int)value;
  return 0;
}

static void *BFMallocArray(size_t count, size_t elemSize, const char *name) {
  void *ptr;
  if (count == 0) return NULL;
  if (elemSize != 0 && count > ((size_t)-1) / elemSize) {
    fprintf(stderr, "Error: allocation size overflow for %s.\n", name);
    exit(EXIT_FAILURE);
  }
  ptr = malloc(count * elemSize);
  if (ptr == NULL) {
    fprintf(stderr, "Error: memory allocation failed for %s.\n", name);
    exit(EXIT_FAILURE);
  }
  return ptr;
}

int BFComputeSizes(int nBackFlowIdx, int nrange, int nzBF,
                   int *nrangeIdx, int *nBFIdxTotal, int *nProjBF) {
  long long rangeIdx;
  long long bfIdxTotal;
  long long projBF;

  if (nBackFlowIdx < 0) {
    fprintf(stderr, "Error: NBackFlowIdx must be non-negative (got %d).\n", nBackFlowIdx);
    return 1;
  }
  if (nBackFlowIdx == 0) {
    *nrangeIdx = 0;
    *nBFIdxTotal = 0;
    *nProjBF = 0;
    return 0;
  }
  if (nrange <= 0) {
    fprintf(stderr, "Error: Nrange must be positive when BackFlow is enabled (got %d).\n", nrange);
    return 1;
  }
  if (nzBF <= 0) {
    fprintf(stderr, "Error: NzBF must be positive when BackFlow is enabled (got %d).\n", nzBF);
    return 1;
  }
  if ((nrange - 1) % nzBF != 0) {
    fprintf(stderr, "Error: Nrange-1 must be divisible by NzBF for BackFlow (Nrange=%d, NzBF=%d).\n",
            nrange, nzBF);
    return 1;
  }

  rangeIdx = 3LL * (long long)(nrange - 1) / (long long)nzBF + 1LL;
  bfIdxTotal = rangeIdx * (rangeIdx + 1LL) / 2LL;
  if (BFCheckedMulLL(bfIdxTotal, nBackFlowIdx, "NProjBF", &projBF) != 0) return 1;
  if (BFIntFromLL(rangeIdx, "NrangeIdx", nrangeIdx) != 0) return 1;
  if (BFIntFromLL(bfIdxTotal, "NBFIdxTotal", nBFIdxTotal) != 0) return 1;
  if (BFIntFromLL(projBF, "NProjBF", nProjBF) != 0) return 1;
  return 0;
}

int BFValidateSettings(int hasBF, int hasBFRange, int backflowSupported) {
  int dummyRangeIdx, dummyBFIdxTotal, dummyProjBF;
  const int hasIntent = (hasBF != 0 || hasBFRange != 0 || NBackFlowIdx != 0 ||
                         Nrange != 0 || NzBF != 0);

  if (NBackFlowIdx < 0) {
    fprintf(stderr, "Error: NBackFlowIdx must be non-negative (got %d).\n", NBackFlowIdx);
    return 1;
  }
  if (Nrange < 0) {
    fprintf(stderr, "Error: Nrange must be non-negative (got %d).\n", Nrange);
    return 1;
  }
  if (NzBF < 0) {
    fprintf(stderr, "Error: NzBF must be non-negative (got %d).\n", NzBF);
    return 1;
  }
  if (!hasIntent) return 0;

  if (!backflowSupported) {
    fprintf(stderr, "Error: Back Flow is not supported in this build.\n");
    return 1;
  }

  if (hasBF != hasBFRange) {
    fprintf(stderr, "Error: BackFlow requires both BF and BFRange definition files.\n");
    return 1;
  }
  if (NBackFlowIdx <= 0) {
    fprintf(stderr, "Error: BackFlow requires NBackFlowIdx > 0.\n");
    return 1;
  }
  if (NBackFlowIdx != 1) {
    fprintf(stderr, "Error: BackFlow MVP supports only NBackFlowIdx==1 (got %d).\n", NBackFlowIdx);
    return 1;
  }
  if (BFComputeSizes(NBackFlowIdx, Nrange, NzBF,
                     &dummyRangeIdx, &dummyBFIdxTotal, &dummyProjBF) != 0) {
    return 1;
  }
  if (NExUpdatePath != 0) {
    fprintf(stderr, "Error: BackFlow MVP supports only NExUpdatePath==0 (got %d).\n", NExUpdatePath);
    return 1;
  }
  if (iFlgOrbitalGeneral != 0) {
    fprintf(stderr, "Error: BackFlow MVP does not support iFlgOrbitalGeneral=1 (FSZ).\n");
    return 1;
  }
  if (NSPGaussLeg != 1) {
    fprintf(stderr, "Error: BackFlow MVP supports only NSPGaussLeg==1 (got %d).\n", NSPGaussLeg);
    return 1;
  }
  if (FlagRBM != 0) {
    fprintf(stderr, "Error: BackFlow MVP does not support RBM.\n");
    return 1;
  }
  if (NTwist > 0) {
    fprintf(stderr, "Error: BackFlow MVP does not support Twist (NTwist=%d).\n", NTwist);
    return 1;
  }
  if (reweight == 1) {
    fprintf(stderr, "Error: BackFlow MVP does not support reweight.\n");
    return 1;
  }
  if (APFlag != 0) {
    fprintf(stderr, "Error: BackFlow MVP does not support APFlag=1.\n");
    return 1;
  }
  if (NQPOptTrans > 1) {
    fprintf(stderr, "Error: BackFlow MVP supports only NQPOptTrans==1 (got %d).\n", NQPOptTrans);
    return 1;
  }
  if (NLanczosMode > 0) {
    fprintf(stderr, "Error: BackFlow MVP does not support NLanczosMode > 0.\n");
    return 1;
  }
  if (NPairHopping > 0 || NExchangeCoupling > 0 || NInterAll > 0 || NCisAjsCktAltDC > 0) {
    fprintf(stderr, "Error: BackFlow MVP does not support two-body terms or GreenFunc2 definitions.\n");
    return 1;
  }
  return 0;
}

int BFDefIntCount(void) {
  long long posCount, rangeCount, nsite2, backflowCount, total;
  int count;
  if (NBackFlowIdx <= 0) return 0;
  if (BFCheckedMulLL(Nsite, Nrange, "PosBF", &posCount) != 0) exit(EXIT_FAILURE);
  if (BFCheckedMulLL(Nsite, Nsite, "RangeIdx", &rangeCount) != 0) exit(EXIT_FAILURE);
  if (BFCheckedMulLL(Nsite, Nsite, "BackFlowIdx rows", &nsite2) != 0) exit(EXIT_FAILURE);
  if (BFCheckedMulLL(nsite2, nsite2, "BackFlowIdx", &backflowCount) != 0) exit(EXIT_FAILURE);
  if (BFCheckedAddLL(posCount, rangeCount, "BackFlow definition table size", &total) != 0) exit(EXIT_FAILURE);
  if (BFCheckedAddLL(total, backflowCount, "BackFlow definition table size", &total) != 0) exit(EXIT_FAILURE);
  if (BFIntFromLL(total, "BackFlow definition table size", &count) != 0) exit(EXIT_FAILURE);
  return count;
}

int BFWorkIntCount(void) {
  long long count;
  int intCount;
  if (NBackFlowIdx <= 0) return 0;
  if (BFCheckedMulLL(16LL * (long long)Nsite, Nrange, "BackFlow work count", &count) != 0) {
    exit(EXIT_FAILURE);
  }
  if (BFIntFromLL(count, "BackFlow work count", &intCount) != 0) exit(EXIT_FAILURE);
  return intCount;
}

void BFBindDefTables(int **pInt) {
  int i;
  int *cursor = *pInt;
  if (NBackFlowIdx <= 0) return;

  PosBF = (int **)BFMallocArray((size_t)Nsite, sizeof(int *), "PosBF");
  for (i = 0; i < Nsite; i++) {
    PosBF[i] = cursor;
    cursor += Nrange;
  }
  RangeIdx = (int **)BFMallocArray((size_t)Nsite, sizeof(int *), "RangeIdx");
  for (i = 0; i < Nsite; i++) {
    RangeIdx[i] = cursor;
    cursor += Nsite;
  }
  BackFlowIdx = (int **)BFMallocArray((size_t)Nsite * (size_t)Nsite, sizeof(int *), "BackFlowIdx");
  for (i = 0; i < Nsite * Nsite; i++) {
    BackFlowIdx[i] = cursor;
    cursor += Nsite * Nsite;
  }
  *pInt = cursor;
}

void BFFreeDefTables(void) {
  free(PosBF);
  free(RangeIdx);
  free(BackFlowIdx);
  PosBF = NULL;
  RangeIdx = NULL;
  BackFlowIdx = NULL;
}

void BFSetupIndex(void) {
  int a, b, cnt;
  if (NBackFlowIdx <= 0 || BFSubIdx == NULL) return;
  cnt = 0;
  for (a = 0; a < NrangeIdx; a++) {
    for (b = a; b < NrangeIdx; b++) {
      BFSubIdx[a][b] = cnt;
      BFSubIdx[b][a] = cnt;
      cnt++;
    }
  }
}

void BFAllocRuntime(void) {
  int i;
  size_t slaterCount;
  if (NBackFlowIdx <= 0) {
    EleProjBFCnt = NULL;
    SmpSltElmBF_real = NULL;
    SmpEta = NULL;
    SmpEtaFlag = NULL;
    SlaterElmBF = NULL;
    SlaterElmBF_real = NULL;
    eta = NULL;
    etaFlag = NULL;
    BFSubIdx = NULL;
    return;
  }

  EleProjBFCnt = (int *)BFMallocArray((size_t)NVMCSample * (size_t)BFWorkIntCount(),
                                      sizeof(int), "EleProjBFCnt");
  slaterCount = (size_t)NQPFull * (size_t)Nsite2 * (size_t)Nsite2;
  SlaterElmBF = (double complex *)BFMallocArray(slaterCount, sizeof(double complex), "SlaterElmBF");
  SlaterElmBF_real = (double *)BFMallocArray(slaterCount, sizeof(double), "SlaterElmBF_real");
  SmpSltElmBF_real = (double *)BFMallocArray((size_t)NVMCSample * slaterCount,
                                             sizeof(double), "SmpSltElmBF_real");
  SmpEta = (double *)BFMallocArray((size_t)NVMCSample * (size_t)NQPFull *
                                   (size_t)Nsite * (size_t)Nsite,
                                   sizeof(double), "SmpEta");
  SmpEtaFlag = (int *)BFMallocArray((size_t)NVMCSample * (size_t)NQPFull *
                                    (size_t)Nsite * (size_t)Nsite,
                                    sizeof(int), "SmpEtaFlag");

  eta = (double complex **)BFMallocArray((size_t)Nsite, sizeof(double complex *), "eta");
  for (i = 0; i < Nsite; i++) {
    eta[i] = (double complex *)BFMallocArray((size_t)Nsite, sizeof(double complex), "eta row");
  }
  etaFlag = (int **)BFMallocArray((size_t)Nsite, sizeof(int *), "etaFlag");
  for (i = 0; i < Nsite; i++) {
    etaFlag[i] = (int *)BFMallocArray((size_t)Nsite, sizeof(int), "etaFlag row");
  }
  BFSubIdx = (int **)BFMallocArray((size_t)NrangeIdx, sizeof(int *), "BFSubIdx");
  for (i = 0; i < NrangeIdx; i++) {
    BFSubIdx[i] = (int *)BFMallocArray((size_t)NrangeIdx, sizeof(int), "BFSubIdx row");
  }
  BFSetupIndex();
}

void BFFreeRuntime(void) {
  int i;
  free(EleProjBFCnt);
  free(SmpSltElmBF_real);
  free(SmpEta);
  free(SmpEtaFlag);
  free(SlaterElmBF);
  free(SlaterElmBF_real);
  if (eta != NULL) {
    for (i = 0; i < Nsite; i++) free(eta[i]);
  }
  free(eta);
  if (etaFlag != NULL) {
    for (i = 0; i < Nsite; i++) free(etaFlag[i]);
  }
  free(etaFlag);
  if (BFSubIdx != NULL) {
    for (i = 0; i < NrangeIdx; i++) free(BFSubIdx[i]);
  }
  free(BFSubIdx);

  EleProjBFCnt = NULL;
  SmpSltElmBF_real = NULL;
  SmpEta = NULL;
  SmpEtaFlag = NULL;
  SlaterElmBF = NULL;
  SlaterElmBF_real = NULL;
  eta = NULL;
  etaFlag = NULL;
  BFSubIdx = NULL;
}
