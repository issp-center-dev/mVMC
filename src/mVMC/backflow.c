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

static int BFCheckSiteIndex(int site) {
  return (site < 0 || site >= Nsite);
}

static int BFSkipBodyHeader(FILE *fp, const char *defname) {
  char ctmp[256];
  int line;
  for (line = 0; line < 5; line++) {
    if (fgets(ctmp, sizeof(ctmp), fp) == NULL) {
      fprintf(stderr, "Error in %s: BackFlow body header is incomplete.\n", defname);
      return 1;
    }
  }
  return 0;
}

static int BFParameterOffset(void) {
  return NProj + FlagRBM * NRBM;
}

static int BFRangeShellMax(void) {
  if (NzBF <= 0 || Nrange <= 0 || (Nrange - 1) % NzBF != 0) return -1;
  return (Nrange - 1) / NzBF;
}

static int BFRangeExpectedRows(long long *expectedRows) {
  return BFCheckedMulLL(Nsite, Nrange, "BFRange row count", expectedRows);
}

static int BFDefinitionExpectedRows(long long *expectedRows) {
  long long nsite2;
  long long nrange2;
  if (BFCheckedMulLL(Nsite, Nsite, "BF row count", &nsite2) != 0) return 1;
  if (BFCheckedMulLL(Nrange, Nrange, "BF row count", &nrange2) != 0) return 1;
  return BFCheckedMulLL(nsite2, nrange2, "BF row count", expectedRows);
}

static int BFContainsRangeSite(int center, int site) {
  int r;
  if (BFCheckSiteIndex(center) != 0 || BFCheckSiteIndex(site) != 0 || PosBF == NULL) return 0;
  for (r = 0; r < Nrange; r++) {
    if (PosBF[center][r] == site) return 1;
  }
  return 0;
}

static int BFValidOptFlagValue(int flag) {
  return (flag == 0 || flag == 1);
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
  if (NPairHopping > 0 || NExchangeCoupling > 0 || NInterAll > 0) {
    fprintf(stderr, "Error: BackFlow MVP does not support two-body Hamiltonian terms.\n");
    return 1;
  }
  return 0;
}

int BFReadRange(FILE *fp, const char *defname) {
  int i;
  int j;
  int n;
  int info = 0;
  int shellMax;
  int *seenCount = NULL;
  int *seenPair = NULL;
  long long row;
  long long expectedRows;

  if (Nrange <= 0) return 0;
  if (BFSkipBodyHeader(fp, defname) != 0) return 1;
  if (PosBF == NULL || RangeIdx == NULL) {
    fprintf(stderr, "Error in %s: BackFlow range tables are not allocated.\n", defname);
    return 1;
  }
  if (BFRangeExpectedRows(&expectedRows) != 0) return 1;
  shellMax = BFRangeShellMax();
  if (shellMax < 0) {
    fprintf(stderr, "Error in %s: invalid BFRange size parameters (Nrange=%d, NzBF=%d).\n",
            defname, Nrange, NzBF);
    return 1;
  }

  seenCount = (int *)calloc((size_t)Nsite, sizeof(int));
  seenPair = (int *)calloc((size_t)Nsite * (size_t)Nsite, sizeof(int));
  if (seenCount == NULL || seenPair == NULL) {
    fprintf(stderr, "Error in %s: memory allocation failed during BFRange validation.\n", defname);
    info = 1;
    goto cleanup;
  }

  for (i = 0; i < Nsite; i++) {
    for (j = 0; j < Nrange; j++) PosBF[i][j] = -1;
  }
  for (i = 0; i < Nsite; i++) {
    for (j = 0; j < Nsite; j++) RangeIdx[i][j] = -1;
  }

  for (row = 0; row < expectedRows; row++) {
    if (fscanf(fp, "%d %d %d", &i, &j, &n) != 3) {
      fprintf(stderr, "Error in %s: failed to read BFRange row %lld of %lld.\n",
              defname, row + 1, expectedRows);
      info = 1;
      goto cleanup;
    }
    if (BFCheckSiteIndex(i) != 0 || BFCheckSiteIndex(j) != 0) {
      fprintf(stderr, "Error in %s: BFRange site index out of range at row %lld (i=%d, j=%d).\n",
              defname, row + 1, i, j);
      info = 1;
      goto cleanup;
    }
    if (n < 0 || n > shellMax) {
      fprintf(stderr, "Error in %s: BFRange shell index %d out of range [0,%d] at row %lld.\n",
              defname, n, shellMax, row + 1);
      info = 1;
      goto cleanup;
    }
    if ((i == j && n != 0) || (i != j && n == 0)) {
      fprintf(stderr, "Error in %s: BFRange self shell must be exactly i==j and n==0 at row %lld.\n",
              defname, row + 1);
      info = 1;
      goto cleanup;
    }
    if (seenPair[i * Nsite + j] != 0) {
      fprintf(stderr, "Error in %s: duplicated BFRange pair (i=%d, j=%d) at row %lld.\n",
              defname, i, j, row + 1);
      info = 1;
      goto cleanup;
    }
    if (seenCount[i] >= Nrange) {
      fprintf(stderr, "Error in %s: too many BFRange entries for site %d.\n", defname, i);
      info = 1;
      goto cleanup;
    }

    PosBF[i][seenCount[i]] = j;
    RangeIdx[i][j] = n;
    seenCount[i]++;
    seenPair[i * Nsite + j] = 1;
  }

  for (i = 0; i < Nsite; i++) {
    if (seenCount[i] != Nrange) {
      fprintf(stderr, "Error in %s: BFRange site %d has %d rows, expected %d.\n",
              defname, i, seenCount[i], Nrange);
      info = 1;
      goto cleanup;
    }
    if (seenPair[i * Nsite + i] == 0 || RangeIdx[i][i] != 0) {
      fprintf(stderr, "Error in %s: BFRange site %d is missing the self row.\n", defname, i);
      info = 1;
      goto cleanup;
    }
  }

cleanup:
  free(seenCount);
  free(seenPair);
  return info;
}

int BFReadDefinition(FILE *fp, int *optFlag, int *countIdx, const char *defname) {
  int i;
  int j;
  int x0;
  int x1;
  int n;
  int idx;
  int flag;
  int offset;
  int rowIndex;
  int colIndex;
  int info = 0;
  int *seenOpt = NULL;
  long long row;
  long long expectedRows;

  if (NBackFlowIdx <= 0) return 0;
  if (BFSkipBodyHeader(fp, defname) != 0) return 1;
  if (BackFlowIdx == NULL) {
    fprintf(stderr, "Error in %s: BackFlow index table is not allocated.\n", defname);
    return 1;
  }
  if (PosBF == NULL || RangeIdx == NULL) {
    fprintf(stderr, "Error in %s: BFRange must be loaded before BF.\n", defname);
    return 1;
  }
  if (BFDefinitionExpectedRows(&expectedRows) != 0) return 1;

  for (i = 0; i < Nsite * Nsite; i++) {
    for (j = 0; j < Nsite * Nsite; j++) BackFlowIdx[i][j] = -1;
  }

  for (row = 0; row < expectedRows; row++) {
    if (fscanf(fp, "%d %d %d %d %d", &i, &j, &x0, &x1, &n) != 5) {
      fprintf(stderr, "Error in %s: failed to read BF row %lld of %lld.\n",
              defname, row + 1, expectedRows);
      info = 1;
      goto cleanup;
    }
    if (BFCheckSiteIndex(i) != 0 || BFCheckSiteIndex(j) != 0 ||
        BFCheckSiteIndex(x0) != 0 || BFCheckSiteIndex(x1) != 0) {
      fprintf(stderr,
              "Error in %s: BF site index out of range at row %lld (i=%d, j=%d, x0=%d, x1=%d).\n",
              defname, row + 1, i, j, x0, x1);
      info = 1;
      goto cleanup;
    }
    if (n < 0 || n >= NBackFlowIdx) {
      fprintf(stderr, "Error in %s: BF group index %d out of range [0,%d) at row %lld.\n",
              defname, n, NBackFlowIdx, row + 1);
      info = 1;
      goto cleanup;
    }
    if (NBackFlowIdx == 1 && n != 0) {
      fprintf(stderr, "Error in %s: BackFlow MVP supports only group n==0 at row %lld.\n",
              defname, row + 1);
      info = 1;
      goto cleanup;
    }
    if (BFContainsRangeSite(i, x0) == 0 || BFContainsRangeSite(j, x1) == 0) {
      fprintf(stderr,
              "Error in %s: BF row %lld uses x0/x1 outside the corresponding BFRange set.\n",
              defname, row + 1);
      info = 1;
      goto cleanup;
    }

    rowIndex = i * Nsite + j;
    colIndex = x0 * Nsite + x1;
    if (BackFlowIdx[rowIndex][colIndex] != -1) {
      fprintf(stderr,
              "Error in %s: duplicated BF entry for (i=%d, j=%d, x0=%d, x1=%d) at row %lld.\n",
              defname, i, j, x0, x1, row + 1);
      info = 1;
      goto cleanup;
    }
    BackFlowIdx[rowIndex][colIndex] = n;
  }

  offset = BFParameterOffset();
  seenOpt = (int *)calloc((size_t)NProjBF, sizeof(int));
  if (seenOpt == NULL) {
    fprintf(stderr, "Error in %s: memory allocation failed during BF OptFlag validation.\n", defname);
    info = 1;
    goto cleanup;
  }

  for (row = 0; row < NProjBF; row++) {
    if (fscanf(fp, "%d %d", &idx, &flag) != 2) {
      fprintf(stderr, "Error in %s: failed to read BF OptFlag row %lld of %d.\n",
              defname, row + 1, NProjBF);
      info = 1;
      goto cleanup;
    }
    if (idx < 0 || idx >= NProjBF) {
      fprintf(stderr, "Error in %s: BF OptFlag index %d out of range [0,%d) at row %lld.\n",
              defname, idx, NProjBF, row + 1);
      info = 1;
      goto cleanup;
    }
    if (!BFValidOptFlagValue(flag)) {
      fprintf(stderr, "Error in %s: BF OptFlag value must be 0 or 1 at row %lld (got %d).\n",
              defname, row + 1, flag);
      info = 1;
      goto cleanup;
    }
    if (seenOpt[idx] != 0) {
      fprintf(stderr, "Error in %s: duplicated BF OptFlag index %d at row %lld.\n",
              defname, idx, row + 1);
      info = 1;
      goto cleanup;
    }

    seenOpt[idx] = 1;
    optFlag[2 * (offset + idx)] = flag;
    optFlag[2 * (offset + idx) + 1] = (idx == 0) ? 0 : ((AllComplexFlag > 0) ? flag : 0);
    (*countIdx)++;
  }

  for (idx = 0; idx < NProjBF; idx++) {
    if (seenOpt[idx] == 0) {
      fprintf(stderr, "Error in %s: missing BF OptFlag index %d.\n", defname, idx);
      info = 1;
      goto cleanup;
    }
  }

cleanup:
  free(seenOpt);
  return info;
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

void BFInitParameters(void) {
  int i;
  if (NProjBF <= 0 || ProjBF == NULL) return;
  ProjBF[0] = 1.0;
  for (i = 1; i < NProjBF; i++) ProjBF[i] = 0.0;
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

void BFRefreshRealLookupTables(void) {
  int a, b, ri, rj, sgn;
  if (NBackFlowIdx <= 0 || BFRealProj == NULL || BFRealSlater == NULL ||
      BFRealSlaterSign == NULL ||
      BFSubIdx == NULL || ProjBF == NULL || Slater == NULL) {
    return;
  }

  BFRealEta = creal(ProjBF[0]);

  for (a = 0; a < NrangeIdx; a++) {
    for (b = 0; b < NrangeIdx; b++) {
      BFRealProj[a * NrangeIdx + b] = -creal(ProjBF[BFSubIdx[a][b]]);
    }
  }

  for (ri = 0; ri < Nsite; ri++) {
    for (rj = 0; rj < Nsite; rj++) {
      sgn = OrbitalSgn[ri][rj];
      if (sgn != 1 && sgn != -1) {
        fprintf(stderr, "Error: OrbitalSgn[%d][%d]=%d is invalid for real BackFlow table.\n",
                ri, rj, sgn);
        exit(EXIT_FAILURE);
      }
      BFRealSlater[ri * Nsite + rj] = creal(Slater[OrbitalIdx[ri][rj]]);
      BFRealSlaterSign[ri * Nsite + rj] = (double)sgn;
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
    BFRealProj = NULL;
    BFRealSlater = NULL;
    BFRealSlaterSign = NULL;
    BFRealEta = 1.0;
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
  BFRealProj = (double *)BFMallocArray((size_t)NrangeIdx * (size_t)NrangeIdx,
                                       sizeof(double), "BFRealProj");
  BFRealSlater = (double *)BFMallocArray((size_t)Nsite * (size_t)Nsite,
                                         sizeof(double), "BFRealSlater");
  BFRealSlaterSign = (double *)BFMallocArray((size_t)Nsite * (size_t)Nsite,
                                             sizeof(double), "BFRealSlaterSign");
  BFRealEta = 1.0;

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
  free(BFRealProj);
  free(BFRealSlater);
  free(BFRealSlaterSign);
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
  BFRealProj = NULL;
  BFRealSlater = NULL;
  BFRealSlaterSign = NULL;
  BFRealEta = 1.0;
}
