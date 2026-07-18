/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

his program is developed based on the mVMC-mini program
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
/*-------------------------------------------------------------
 * Variational Monte Carlo
 * local Green Functions
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

#include "sector_projection.c"
#include <stdint.h>

int CollectBFAffectedParticlesTwoMove_fsz(
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    int movedParticle0, int movedParticle1,
    const int *eleIdx, const int *eleSpn,
    int *affected, int *intWork, int intWorkSize,
    int *nChangedOut, int *nAffectedOut);

int ClassifyGreenFunc2FSZ(const int XI, const int XJ,
                          const int XK, const int XL) {
  if(XI == XJ) {
    if(XJ == XK) return (XK == XL) ? 1 : 2;
    if(XJ == XL) return 3;
    if(XK == XL) return 4;
    return 5;
  }
  if(XI == XK) {
    if(XJ == XL) return 6;
    if(XK == XL) return 7;
    return 8;
  }
  if(XI == XL) return (XJ == XK) ? 9 : 10;
  if(XJ == XK) return (XK == XL) ? 11 : 12;
  if(XJ == XL) return 13;
  if(XK == XL) return 14;
  return 15;
}

static void RecordBFFSZC2DispatchSeconds(BFFSZC2DetailContext *context,
                                         const double start) {
  if(context == NULL) return;
  context->componentSeconds[BFFSZ_C2_DETAIL_COMPONENT_DISPATCH]
      += BFFSZC2DetailMonotonicSeconds()-start;
}

static void RecordBFFSZC2Path(BFFSZC2DetailContext *context, const int path) {
  if(context == NULL || path < 0 || path >= NBFFSZC2DetailPath) return;
  context->path[path]++;
}

static int *BF_FSZ_CaptureGreen2State(const int *eleIdx, const int *eleSpn,
                                      const int *eleCfg, const int *eleNum) {
  size_t count;
  int *snapshot;

  if(eleIdx == NULL || eleSpn == NULL || eleCfg == NULL || eleNum == NULL
      || Nsize <= 0 || Nsite2 <= 0
      || (size_t)Nsize > SIZE_MAX/2
      || 2*(size_t)Nsize > SIZE_MAX-2*(size_t)Nsite2) return NULL;
  count = 2*(size_t)Nsize+2*(size_t)Nsite2;
  if(count > SIZE_MAX/sizeof(int)) return NULL;
  snapshot = (int *)malloc(sizeof(int)*count);
  if(snapshot == NULL) return NULL;
  memcpy(snapshot,eleIdx,sizeof(int)*(size_t)Nsize);
  memcpy(snapshot+Nsize,eleSpn,sizeof(int)*(size_t)Nsize);
  memcpy(snapshot+2*Nsize,eleCfg,sizeof(int)*(size_t)Nsite2);
  memcpy(snapshot+2*Nsize+Nsite2,eleNum,sizeof(int)*(size_t)Nsite2);
  return snapshot;
}

static int BF_FSZ_Green2StateMatches(const int *snapshot,
                                     const int *eleIdx, const int *eleSpn,
                                     const int *eleCfg, const int *eleNum) {
  if(snapshot == NULL || eleIdx == NULL || eleSpn == NULL
      || eleCfg == NULL || eleNum == NULL) return 0;
  return memcmp(snapshot,eleIdx,sizeof(int)*(size_t)Nsize) == 0
      && memcmp(snapshot+Nsize,eleSpn,sizeof(int)*(size_t)Nsize) == 0
      && memcmp(snapshot+2*Nsize,eleCfg,sizeof(int)*(size_t)Nsite2) == 0
      && memcmp(snapshot+2*Nsize+Nsite2,eleNum,
          sizeof(int)*(size_t)Nsite2) == 0;
}

static int BF_FSZ_IntWorkspacesOverlap(const int *workA, const int workACount,
                                       const int *workB, const int workBCount) {
  size_t workABytes, workBBytes;
  uintptr_t workABegin, workBBegin;

  if(workA == NULL || workB == NULL || workACount < 0 || workBCount < 0
      || (size_t)workACount > SIZE_MAX/sizeof(int)
      || (size_t)workBCount > SIZE_MAX/sizeof(int)) {
    return 1;
  }
  workABytes = sizeof(int)*(size_t)workACount;
  workBBytes = sizeof(int)*(size_t)workBCount;
  workABegin = (uintptr_t)(const void *)workA;
  workBBegin = (uintptr_t)(const void *)workB;
  if(workABegin > UINTPTR_MAX - workABytes
      || workBBegin > UINTPTR_MAX - workBBytes) return 1;
  return workABegin < workBBegin + workBBytes
      && workBBegin < workABegin + workABytes;
}

double complex GreenFunc1_fsz(const int ri, const int rj, const int s, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double complex *buffer);

double complex GreenFunc1_fsz2(const int ri, const int rj, const int s,const int t, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double complex *buffer);

double complex GreenFunc2_fsz(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double complex  ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double complex *buffer);

double complex GreenFunc1BF_fsz(const int ri, const int rj, const int s, const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork);

double complex GreenFunc1BF_fsz2(const int ri, const int rj, const int s, const int t,
                  const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork);

double complex GreenFunc1BF_fsz_workspace(const int ri, const int rj, const int s,
                  const double complex ip, int *eleIdx, int *eleCfg, int *eleNum,
                  const int *eleProjCnt, int *eleSpn, int *projCntNew,
                  const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, double complex *pfBufM,
                  int *affected, int *hopIntWork, int hopIntWorkSize, int *pfIWork,
                  double complex *pfWork, double *pfRWork);

double complex GreenFunc1BF_fsz2_workspace(const int ri, const int rj,
                  const int s, const int t, const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum,
                  const int *eleProjCnt, int *eleSpn, int *projCntNew,
                  const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, double complex *pfBufM,
                  int *affected, int *hopIntWork, int hopIntWorkSize, int *pfIWork,
                  double complex *pfWork, double *pfRWork);

double complex GreenFunc2BF_fsz(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork);

double complex GreenFunc2BF_fsz2(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const int u, const int v,
                  const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork);

double complex GreenFunc2BF_fsz2WithProfile(
                  const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const int u, const int v,
                  const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork,
                  BFFSZC2DetailContext *detailProfile);

double complex GreenFuncN_fsz(const int n, int *rsi, int *rsj, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  int *eleSpn, double complex *buffer, int *bufferInt, double *rwork);

double complex calculateNewPfMN_child_fsz(const int qpidx, const int n, const int *msa,
                              const int *eleIdx, const int *eleSpn,
                              double complex *buffer, double *rwork);

static double complex BF_FSZ_NaNValue(void) {
  const double nanValue = NAN;
  return nanValue + I*nanValue;
}

static double complex CalculateBF_FSZ_CandidateIP(const char *caller,
    const double complex *sltElmBF, const int *eleIdx, const int *eleSpn,
    double complex *pfMNew, double complex *pfBufM, int *pfIWork,
    double complex *pfWork, double *pfRWork) {
  int status;
  int detail = 0;
  double complex ipNew;

  status = CalculatePfM_BF_fsz_from_workspace(sltElmBF,eleIdx,eleSpn,
      0,NQPFull,pfMNew,&detail,pfBufM,pfIWork,pfWork,LapackLWork,pfRWork);
  if(status == BF_FSZ_PF_LAPACK_FAILURE) {
    fprintf(stderr, "warning: %s: BF-FSZ candidate Pfaffian failed in M_ZSKPFA (info=%d); propagating NaN.\n",
        caller, detail);
    return BF_FSZ_NaNValue();
  }
  if(status == BF_FSZ_PF_NONFINITE) {
    fprintf(stderr, "warning: %s: BF-FSZ candidate Pfaffian is non-finite (qpidx=%d); propagating NaN.\n",
        caller, detail);
    return BF_FSZ_NaNValue();
  }
  if(status != BF_FSZ_PF_OK) {
    fprintf(stderr, "warning: %s: BF-FSZ candidate Pfaffian returned unknown status %d; propagating NaN.\n",
        caller, status);
    return BF_FSZ_NaNValue();
  }

  ipNew = CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);
  if(!(isfinite(creal(ipNew)) && isfinite(cimag(ipNew)))) {
    fprintf(stderr, "warning: %s: BF-FSZ candidate projected overlap is non-finite; propagating NaN.\n",
        caller);
    return BF_FSZ_NaNValue();
  }
  return ipNew;
}

static int BF_FSZ_GreenPfUpdateMatchesFull(const double complex *updated,
                                           const double complex *full,
                                           const int qpNum, int *badQp,
                                           double *badError) {
  int qpidx;
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    const double scale = fmax(1.0,cabs(full[qpidx]));
    const double error = cabs(updated[qpidx]-full[qpidx])/scale;
    if(!(isfinite(error) && error <= 1.0e-10)) {
      if(badQp != NULL) *badQp = qpidx;
      if(badError != NULL) *badError = error;
      return 0;
    }
  }
  return 1;
}

static int CalculateBF_FSZ_FullCandidatePfM(const double complex *sltElmBF,
    const int *eleIdx, const int *eleSpn, double complex *pfMOut,
    int *failureDetail, double complex *pfBufM, int *pfIWork,
    double complex *pfWork, double *pfRWork) {
  if(pfBufM != NULL && pfIWork != NULL && pfWork != NULL && pfRWork != NULL) {
    return CalculatePfM_BF_fsz_from_workspace(sltElmBF, eleIdx, eleSpn, 0,
        NQPFull, pfMOut, failureDetail, pfBufM, pfIWork, pfWork,
        LapackLWork, pfRWork);
  }
  return CalculatePfM_BF_fsz_from(sltElmBF, eleIdx, eleSpn, 0, NQPFull,
      pfMOut, failureDetail);
}

static int BF_FSZ_GreenRowsMatchFull(
    const double complex *candidateRows, const size_t rowQpStride,
    const size_t rowAffectedStride, const int nAffected,
    const int *affected, const int *eleIdx, const int *eleSpn,
    const double complex *candidateSlater, size_t *badIndex) {
  int qpidx,k,i;
  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    const double complex *sltE = candidateSlater
        +(size_t)qpidx*(size_t)Nsite2*(size_t)Nsite2;
    const double complex *rowsQp
        = candidateRows+(size_t)qpidx*rowQpStride;
    for(k=0;k<nAffected;k++) {
      const int rsk = eleIdx[affected[k]]+eleSpn[affected[k]]*Nsite;
      const double complex *row
          = rowsQp+(size_t)k*rowAffectedStride;
      for(i=0;i<Nsize;i++) {
        const int rsi = eleIdx[i]+eleSpn[i]*Nsite;
        if(row[i] != sltE[(size_t)rsk*(size_t)Nsite2+(size_t)rsi]) {
          if(badIndex != NULL) {
            *badIndex = (size_t)qpidx*rowQpStride
                +(size_t)k*rowAffectedStride+(size_t)i;
          }
          return 0;
        }
      }
    }
  }
  return 1;
}

static int BF_FSZ_GreenAffectedCoversFull(
    const double complex *candidateSlater,
    const int *eleIdx, const int *eleSpn,
    const int nAffected, const int *affected,
    int *badQp, int *badM, int *badN) {
  int qpidx,m,n,k;
  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    const size_t qpOffset
        = (size_t)qpidx*(size_t)Nsite2*(size_t)Nsite2;
    const double complex *baseSltE = SlaterElmBF+qpOffset;
    const double complex *candidateSltE = candidateSlater+qpOffset;
    for(m=0;m<Nsize;m++) {
      const int rsm = eleIdx[m]+eleSpn[m]*Nsite;
      int affectedM = 0;
      for(k=0;k<nAffected;k++) if(affected[k] == m) affectedM = 1;
      for(n=m+1;n<Nsize;n++) {
        const int rsn = eleIdx[n]+eleSpn[n]*Nsite;
        int affectedN = 0;
        for(k=0;k<nAffected;k++) if(affected[k] == n) affectedN = 1;
        if(baseSltE[(size_t)rsm*(size_t)Nsite2+(size_t)rsn]
            != candidateSltE[(size_t)rsm*(size_t)Nsite2+(size_t)rsn]
            && !affectedM && !affectedN) {
          if(badQp != NULL) *badQp = qpidx;
          if(badM != NULL) *badM = m;
          if(badN != NULL) *badN = n;
          return 0;
        }
      }
    }
  }
  return 1;
}

static int BF_FSZ_MultiMoveRowsMatchFullSlice(
    const double complex *candidateRows, const size_t rowQpStride,
    const size_t rowAffectedStride, const int nAffected,
    const int *affected, const int *eleIdx, const int *eleSpn,
    const int qpStart, const int qpEnd,
    const double complex *candidateSlater) {
  int qpidx,k,i;
  for(qpidx=qpStart;qpidx<qpEnd;qpidx++) {
    const double complex *sltE = candidateSlater
        +(size_t)qpidx*(size_t)Nsite2*(size_t)Nsite2;
    const double complex *rowsQp = candidateRows
        +(size_t)(qpidx-qpStart)*rowQpStride;
    for(k=0;k<nAffected;k++) {
      const int rsk = eleIdx[affected[k]]+eleSpn[affected[k]]*Nsite;
      const double complex *row = rowsQp+(size_t)k*rowAffectedStride;
      for(i=0;i<Nsize;i++) {
        const int rsi = eleIdx[i]+eleSpn[i]*Nsite;
        if(row[i] != sltE[(size_t)rsk*(size_t)Nsite2+(size_t)rsi]) {
          return 0;
        }
      }
    }
  }
  return 1;
}

static int BF_FSZ_MultiMoveRowsEqual(
    const double complex *rowsA, const double complex *rowsB,
    const size_t rowQpStride, const size_t rowAffectedStride,
    const int nAffected, const int qpNum) {
  int qpidx,k;
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    for(k=0;k<nAffected;k++) {
      const double complex *rowA = rowsA+(size_t)qpidx*rowQpStride
          +(size_t)k*rowAffectedStride;
      const double complex *rowB = rowsB+(size_t)qpidx*rowQpStride
          +(size_t)k*rowAffectedStride;
      if(memcmp(rowA,rowB,sizeof(double complex)*(size_t)Nsize) != 0) {
        return 0;
      }
    }
  }
  return 1;
}

static int BF_FSZ_AffectedContainsParticle(const int *affected,
                                           const int nAffected,
                                           const int particle) {
  int k;
  for(k=0;k<nAffected;k++) if(affected[k] == particle) return 1;
  return 0;
}

static int BF_FSZ_RunMultiMoveArgumentChecks(
    double complex *candidateRows, const size_t rowQpStride,
    const size_t rowAffectedStride, const double complex *baseSlater,
    const int *movedParticles, const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    const int qpStart, const int qpEnd,
    int *affected, int *intWork, const int intWorkSize) {
  const double complex sentinel = 1234.5-678.25*I;
  const double complex savedBase = baseSlater[0];
  int badNegative[2] = {-1,movedParticles[1]};
  int badLarge[2] = {movedParticles[0],Nsize};
  int duplicate[2] = {movedParticles[0],movedParticles[0]};
  int status,nChanged,nAffected;

#define BF_FSZ_EXPECT_MULTI_INVALID(call) do { \
    candidateRows[0] = sentinel; \
    status = (call); \
    if(status != BF_FSZ_ROW_BUILD_INVALID_ARGUMENT \
        || candidateRows[0] != sentinel) return 1; \
  } while(0)
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      0,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      -1,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      Nsize+1,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      2,NULL,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      2,badNegative,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      2,badLarge,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      2,duplicate,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,-1,baseSlater,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,(size_t)Nsize-1,Nsize,baseSlater,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,(size_t)Nsize-1,rowAffectedStride,Nsize,baseSlater,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,SIZE_MAX,Nsize,baseSlater,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      -1,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpEnd,qpStart,affected,&nChanged,&nAffected,intWork,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,NQPFull+1,affected,&nChanged,&nAffected,intWork,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,NULL,intWorkSize));
  BF_FSZ_EXPECT_MULTI_INVALID(MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize-1));
#undef BF_FSZ_EXPECT_MULTI_INVALID

  status = MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      NULL,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize);
  if(status != BF_FSZ_ROW_BUILD_INVALID_ARGUMENT) return 1;

  candidateRows[0] = sentinel;
  status = MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,1,baseSlater,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize);
  if(status != BF_FSZ_ROW_BUILD_NEEDS_FULL || candidateRows[0] != sentinel
      || nAffected < 2) return 1;

  status = MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      NULL,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpStart,affected,&nChanged,&nAffected,intWork,intWorkSize);
  if(status != BF_FSZ_ROW_BUILD_OK) return 1;

  status = MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      (double complex *)baseSlater,rowQpStride,rowAffectedStride,Nsize,
      baseSlater,2,movedParticles,eleIdx,eleSpn,
      oldEleProjBFCnt,newEleProjBFCnt,qpStart,qpEnd,
      affected,&nChanged,&nAffected,intWork,intWorkSize);
  if(status != BF_FSZ_ROW_BUILD_INVALID_ARGUMENT || baseSlater[0] != savedBase) {
    return 1;
  }
  status = MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      (double complex *)intWork,rowQpStride,rowAffectedStride,Nsize,
      baseSlater,2,movedParticles,eleIdx,eleSpn,
      oldEleProjBFCnt,newEleProjBFCnt,qpStart,qpEnd,
      affected,&nChanged,&nAffected,intWork,intWorkSize);
  if(status != BF_FSZ_ROW_BUILD_INVALID_ARGUMENT) return 1;
  status = MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      (double complex *)affected,rowQpStride,rowAffectedStride,Nsize,
      baseSlater,2,movedParticles,eleIdx,eleSpn,
      oldEleProjBFCnt,newEleProjBFCnt,qpStart,qpEnd,
      affected,&nChanged,&nAffected,intWork,intWorkSize);
  if(status != BF_FSZ_ROW_BUILD_INVALID_ARGUMENT) return 1;
  status = MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      (double complex *)movedParticles,rowQpStride,rowAffectedStride,Nsize,
      baseSlater,2,movedParticles,eleIdx,eleSpn,
      oldEleProjBFCnt,newEleProjBFCnt,qpStart,qpEnd,
      affected,&nChanged,&nAffected,intWork,intWorkSize);
  if(status != BF_FSZ_ROW_BUILD_INVALID_ARGUMENT) return 1;

  affected[0] = movedParticles[0];
  affected[1] = movedParticles[1];
  status = MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      2,affected,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize);
  if(status != BF_FSZ_ROW_BUILD_INVALID_ARGUMENT) return 1;

  intWork[0] = movedParticles[0];
  intWork[1] = movedParticles[1];
  status = MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,Nsize,baseSlater,
      2,intWork,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,intWork,intWorkSize);
  if(status != BF_FSZ_ROW_BUILD_INVALID_ARGUMENT) return 1;
  return 0;
}

static int BF_FSZ_RunMultiMoveRowCheck(
    const int movedParticle0, const int movedParticle1,
    const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    const double complex *candidateSlater,
    int *hopIntWork, const int hopIntWorkSize) {
  const int qpStart = (NQPFull > 1) ? 1 : 0;
  const int qpEnd = NQPFull;
  const int qpNum = qpEnd-qpStart;
  const size_t rowAffectedStride = (size_t)Nsize;
  const size_t rowQpStride = (size_t)Nsize*rowAffectedStride;
  size_t rowCount;
  double complex *rows = NULL, *serialRows = NULL, *reverseRows = NULL;
  int *affected = NULL, *serialAffected = NULL, *reverseAffected = NULL;
  int movedParticles[2] = {movedParticle0,movedParticle1};
  int reverseMovedParticles[2] = {movedParticle1,movedParticle0};
  int status,nChanged,nAffected,serialStatus,serialChanged,serialNAffected;
  int reverseStatus,reverseChanged,reverseNAffected,k;
  int result = 1;

  if(movedParticle0 == movedParticle1
      || GetSlaterElmBF_fsz_hop_row_work_size(&rowCount,qpNum,Nsize)
          != BF_FSZ_ROW_BUILD_OK
      || rowCount > SIZE_MAX/sizeof(double complex)) return 1;
  rows = (double complex *)malloc(sizeof(double complex)*rowCount);
  serialRows = (double complex *)malloc(sizeof(double complex)*rowCount);
  reverseRows = (double complex *)malloc(sizeof(double complex)*rowCount);
  affected = (int *)malloc(sizeof(int)*(size_t)Nsize);
  serialAffected = (int *)malloc(sizeof(int)*(size_t)Nsize);
  reverseAffected = (int *)malloc(sizeof(int)*(size_t)Nsize);
  if(rows == NULL || serialRows == NULL || reverseRows == NULL
      || affected == NULL || serialAffected == NULL || reverseAffected == NULL) {
    goto cleanup;
  }

  status = MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      rows,rowQpStride,rowAffectedStride,Nsize,SlaterElmBF,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,
      hopIntWork,hopIntWorkSize);
  serialStatus = MakeSlaterElmBF_fsz_multi_move_rows_workspace_serial(
      serialRows,rowQpStride,rowAffectedStride,Nsize,SlaterElmBF,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,serialAffected,&serialChanged,&serialNAffected,
      hopIntWork,hopIntWorkSize);
  reverseStatus = MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      reverseRows,rowQpStride,rowAffectedStride,Nsize,SlaterElmBF,
      2,reverseMovedParticles,eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,
      qpStart,qpEnd,reverseAffected,&reverseChanged,&reverseNAffected,
      hopIntWork,hopIntWorkSize);
  if(status != BF_FSZ_ROW_BUILD_OK || serialStatus != status
      || reverseStatus != status || serialChanged != nChanged
      || reverseChanged != nChanged || serialNAffected != nAffected
      || reverseNAffected != nAffected || nAffected < 2
      || !BF_FSZ_AffectedContainsParticle(affected,nAffected,movedParticle0)
      || !BF_FSZ_AffectedContainsParticle(affected,nAffected,movedParticle1)
      || !BF_FSZ_MultiMoveRowsMatchFullSlice(
          rows,rowQpStride,rowAffectedStride,nAffected,affected,
          eleIdx,eleSpn,qpStart,qpEnd,candidateSlater)
      || !BF_FSZ_MultiMoveRowsEqual(
          rows,serialRows,rowQpStride,rowAffectedStride,nAffected,qpNum)
      || !BF_FSZ_MultiMoveRowsEqual(
          rows,reverseRows,rowQpStride,rowAffectedStride,nAffected,qpNum)) {
    goto cleanup;
  }
  for(k=0;k<nAffected;k++) {
    if((k > 0 && affected[k-1] >= affected[k])
        || serialAffected[k] != affected[k]
        || reverseAffected[k] != affected[k]) goto cleanup;
  }

  /* With unchanged BF counts, only explicit moved-particle forcing remains.
   * Omitting either moved particle must therefore fail the coverage oracle. */
  status = MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      rows,rowQpStride,rowAffectedStride,Nsize,SlaterElmBF,
      2,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,oldEleProjBFCnt,
      qpStart,qpEnd,affected,&nChanged,&nAffected,
      hopIntWork,hopIntWorkSize);
  serialStatus = MakeSlaterElmBF_fsz_hop_rows_workspace(
      serialRows,rowQpStride,rowAffectedStride,Nsize,SlaterElmBF,
      movedParticle0,eleIdx,eleSpn,oldEleProjBFCnt,oldEleProjBFCnt,
      qpStart,qpEnd,serialAffected,&serialChanged,&serialNAffected,
      hopIntWork,hopIntWorkSize);
  reverseStatus = MakeSlaterElmBF_fsz_multi_move_rows_workspace(
      reverseRows,rowQpStride,rowAffectedStride,Nsize,SlaterElmBF,
      1,movedParticles,eleIdx,eleSpn,oldEleProjBFCnt,oldEleProjBFCnt,
      qpStart,qpEnd,reverseAffected,&reverseChanged,&reverseNAffected,
      hopIntWork,hopIntWorkSize);
  if(status != BF_FSZ_ROW_BUILD_OK || nChanged != 0 || nAffected != 2
      || !BF_FSZ_MultiMoveRowsMatchFullSlice(
          rows,rowQpStride,rowAffectedStride,nAffected,affected,
          eleIdx,eleSpn,qpStart,qpEnd,SlaterElmBF)
      || serialStatus != BF_FSZ_ROW_BUILD_OK || reverseStatus != serialStatus
      || serialChanged != 0 || reverseChanged != serialChanged
      || serialNAffected != 1 || reverseNAffected != serialNAffected
      || serialAffected[0] != movedParticle0
      || reverseAffected[0] != serialAffected[0]
      || BF_FSZ_AffectedContainsParticle(
          serialAffected,serialNAffected,movedParticle1)
      || !BF_FSZ_MultiMoveRowsEqual(
          serialRows,reverseRows,rowQpStride,rowAffectedStride,
          serialNAffected,qpNum)) {
    goto cleanup;
  }

  if(BF_FSZ_RunMultiMoveArgumentChecks(
      rows,rowQpStride,rowAffectedStride,SlaterElmBF,movedParticles,
      eleIdx,eleSpn,oldEleProjBFCnt,newEleProjBFCnt,qpStart,qpEnd,
      affected,hopIntWork,hopIntWorkSize) != 0) goto cleanup;
  result = 0;

cleanup:
  free(rows);
  free(serialRows);
  free(reverseRows);
  free(affected);
  free(serialAffected);
  free(reverseAffected);
  return result;
}

static void BF_FSZ_CheckMultiMoveRowsOnce(
    const int movedParticle0, const int movedParticle1,
    const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    const double complex *candidateSlater,
    int *hopIntWork, const int hopIntWorkSize) {
  static int checked = 0;
  if(!BFFSZMultiMoveRowCheckEnabled) return;
  #pragma omp critical(BF_FSZ_multi_move_row_check)
  {
    if(!checked) {
      if(BF_FSZ_RunMultiMoveRowCheck(
          movedParticle0,movedParticle1,eleIdx,eleSpn,
          oldEleProjBFCnt,newEleProjBFCnt,candidateSlater,
          hopIntWork,hopIntWorkSize) != 0) {
        fprintf(stderr,"error: BF-FSZ multi-move candidate-row check failed\n");
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
      }
      checked = 1;
    }
  }
}

int GetGreenFuncBF_fsz_buffer_work_size(
    size_t *bufferComplexCount, size_t *pfUpdateIntCount,
    size_t *pfUpdateDoubleCount) {
  const int rowCapacity = (BFFSZPfUpdateKFull-1 < Nsize-1)
      ? BFFSZPfUpdateKFull-1 : Nsize-1;
  size_t slaterSquare,bfSlaterSize,rowCount,pfUpdateComplexCount;
  size_t bufferCount;

  if(bufferComplexCount == NULL || pfUpdateIntCount == NULL
      || pfUpdateDoubleCount == NULL || NQPFull <= 0 || Nsite2 <= 0
      || Nsize <= 0 || rowCapacity < 0
      || (size_t)Nsite2 > SIZE_MAX/(size_t)Nsite2) return -1;
  slaterSquare = (size_t)Nsite2*(size_t)Nsite2;
  if((size_t)NQPFull > SIZE_MAX/slaterSquare) return -1;
  bfSlaterSize = (size_t)NQPFull*slaterSquare;
  if(GetCalculateNewPfMBF_fsz_row_values_work_size(
        &pfUpdateComplexCount,pfUpdateIntCount,pfUpdateDoubleCount) != 0
      || GetSlaterElmBF_fsz_hop_row_work_size(
        &rowCount,NQPFull,rowCapacity) != BF_FSZ_ROW_BUILD_OK
      || *pfUpdateIntCount > INT_MAX || *pfUpdateDoubleCount > INT_MAX
      || (size_t)NQPFull > SIZE_MAX/2
      || bfSlaterSize > SIZE_MAX-2*(size_t)NQPFull
      || bfSlaterSize+2*(size_t)NQPFull > SIZE_MAX-rowCount
      || bfSlaterSize+2*(size_t)NQPFull+rowCount
          > SIZE_MAX-pfUpdateComplexCount) return -1;
  bufferCount = 2*(size_t)NQPFull+bfSlaterSize+rowCount
      +pfUpdateComplexCount;
  if(BFFSZGreenRebuildCheckEnabled) {
    if(bufferCount > SIZE_MAX-bfSlaterSize) return -1;
    bufferCount += bfSlaterSize;
  }
  *bufferComplexCount = bufferCount;
  return 0;
}

static void BF_FSZ_MaterializeGreenCandidate(
    const char *caller, double complex *candidateSlater,
    double complex *debugFull, const int *eleNum,
    const int *eleIdx, const int *eleSpn,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    const double complex *candidateRows, const size_t rowQpStride,
    const size_t rowAffectedStride, const int rowsBuilt,
    const int nAffected, const int *affected,
    int *hopIntWork, const int hopIntWorkSize,
    const int useOMP, const int detailTimers, const int reason,
    double *materializeSeconds) {
  const size_t slaterCount
      = (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
  int status;
  double materializeStart = 0.0;

  if(materializeSeconds != NULL) {
    materializeStart = BFFSZC2DetailMonotonicSeconds();
  }

  if(detailTimers) {
    StopTimer(83);
    StartTimer(82);
  }
  memcpy(candidateSlater,SlaterElmBF,
      sizeof(double complex)*slaterCount);
  status = useOMP
      ? CommitSlaterElmBF_fsz_hop_workspace(
          candidateSlater,oldEleProjBFCnt,newEleProjBFCnt,
          hopIntWork,hopIntWorkSize)
      : CommitSlaterElmBF_fsz_hop_workspace_serial(
          candidateSlater,oldEleProjBFCnt,newEleProjBFCnt,
          hopIntWork,hopIntWorkSize);
  if(detailTimers) {
    StopTimer(82);
    StartTimer(83);
  }
  if(status != BF_FSZ_ROW_BUILD_OK) {
    fprintf(stderr,"error: %s: BF-FSZ full candidate materialization failed\n",
        caller);
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  RecordBFFSZGreenMaterialize(reason);

  if(rowsBuilt && BFFSZMatrixFreeCheckEnabled) {
    size_t badIndex=0;
    if(!BF_FSZ_GreenRowsMatchFull(
        candidateRows,rowQpStride,rowAffectedStride,nAffected,affected,
        eleIdx,eleSpn,candidateSlater,&badIndex)) {
      fprintf(stderr,
          "error: %s: BF-FSZ Green matrix-free row/full mismatch at index %zu\n",
          caller,badIndex);
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
  }
  if(BFFSZAffectedCheckEnabled) {
    int badQp=-1,badM=-1,badN=-1;
    if(!BF_FSZ_GreenAffectedCoversFull(
        candidateSlater,eleIdx,eleSpn,nAffected,affected,
        &badQp,&badM,&badN)) {
      fprintf(stderr,
          "error: %s: BF-FSZ Green affected coverage mismatch "
          "at qp=%d particles=(%d,%d)\n",caller,badQp,badM,badN);
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
  }
  if(BFFSZGreenRebuildCheckEnabled) {
    size_t idx;
    if(debugFull == NULL) {
      fprintf(stderr,"error: %s: missing BF-FSZ Green rebuild workspace\n",
          caller);
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    MakeSlaterElmBF_fsz_to_serial(debugFull,eleNum,newEleProjBFCnt);
    for(idx=0;idx<slaterCount;idx++) {
      if(candidateSlater[idx] != debugFull[idx]) {
        fprintf(stderr,
            "error: %s: BF-FSZ Green incremental/full rebuild mismatch "
            "at index %zu\n",caller,idx);
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
      }
    }
  }
  if(materializeSeconds != NULL) {
    *materializeSeconds += BFFSZC2DetailMonotonicSeconds()-materializeStart;
  }
}

static double complex CalculateBF_FSZ_UpdatedCandidateIP(const char *caller,
    const double complex *candidateRows, const size_t rowQpStride,
    const size_t rowAffectedStride, const int rowsBuilt,
    double complex *candidateSlater, double complex *debugFull,
    const int *eleIdx, const int *eleSpn, const int *eleNum,
    const int *oldEleProjBFCnt, const int *newEleProjBFCnt,
    const int nAffected, const int *affected,
    double complex *pfMNew, double complex *pfMFull,
    double complex *pfUpdateWork,
    const size_t pfUpdateComplexCount, int *pfUpdateIWork,
    const size_t pfUpdateIntCount, double *pfUpdateRWork,
    const size_t pfUpdateDoubleCount, double complex *pfBufM,
    int *pfIWork, double complex *pfWork, double *pfRWork,
    int *hopIntWork, const int hopIntWorkSize,
    const int useOMP, const int detailTimers, const int forceOracle,
    BFFSZC2DetailContext *detailProfile, double *materializeSeconds) {
  const unsigned char candidateSentinel = 0xa5;
  size_t candidateCount = 0;
  int status, detail = 0;
  double complex ipNew;

  if(BFFSZC2BufferCheckEnabled && detailProfile != NULL) {
    if(NQPFull <= 0 || Nsite2 <= 0
        || (size_t)Nsite2 > SIZE_MAX/(size_t)Nsite2
        || (size_t)NQPFull
            > SIZE_MAX/((size_t)Nsite2*(size_t)Nsite2)) {
      fprintf(stderr,"error: %s: BF-FSZ candidate sentinel size overflow\n",
          caller);
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    candidateCount = (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
    if(candidateCount > SIZE_MAX/sizeof(double complex)) {
      fprintf(stderr,"error: %s: BF-FSZ candidate sentinel byte overflow\n",
          caller);
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    memset(candidateSlater,candidateSentinel,
        sizeof(double complex)*candidateCount);
  }

  if(!rowsBuilt) {
    RecordBFFSZPfPath(BFFSZ_PROFILE_GREEN,BFFSZ_PF_PATH_DIRECT_FULL);
    RecordBFFSZC2Path(detailProfile,BFFSZ_C2_DETAIL_PATH_DIRECT_FULL);
    BF_FSZ_MaterializeGreenCandidate(
        caller,candidateSlater,debugFull,eleNum,eleIdx,eleSpn,
        oldEleProjBFCnt,newEleProjBFCnt,candidateRows,rowQpStride,
        rowAffectedStride,rowsBuilt,nAffected,affected,
        hopIntWork,hopIntWorkSize,useOMP,detailTimers,
        BFFSZ_GREEN_MATERIALIZE_DIRECT_FULL,materializeSeconds);
    status = CalculateBF_FSZ_FullCandidatePfM(candidateSlater,eleIdx,eleSpn,
        pfMNew,&detail,pfBufM,pfIWork,pfWork,pfRWork);
  } else {
    int updateStatus = CalculateNewPfMBF_fsz_row_values_workspace(
        nAffected,affected,pfMNew,PfM,InvM,(size_t)Nsize*(size_t)Nsize,
        0,NQPFull,candidateRows,rowQpStride,rowAffectedStride,
        &detail,pfUpdateWork,
        pfUpdateComplexCount,pfUpdateIWork,pfUpdateIntCount,pfUpdateRWork,
        pfUpdateDoubleCount);
    if(updateStatus == BF_FSZ_PF_UPDATE_OK) {
      if(BFFSZPfUpdateInjectedStatus != BF_FSZ_PF_UPDATE_OK) {
        updateStatus = BFFSZPfUpdateInjectedStatus;
      } else if(BFFSZPfUpdateForceFallback) {
        updateStatus = BF_FSZ_PF_UPDATE_NONFINITE;
      }
    }
    if(updateStatus == BF_FSZ_PF_UPDATE_INVALID_ARGUMENT) {
      fprintf(stderr,"error: %s: invalid BF-FSZ Pfaffian-update arguments\n",caller);
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    if(updateStatus == BF_FSZ_PF_UPDATE_OK) {
      if(BFFSZC2BufferCheckEnabled && detailProfile != NULL) {
        const unsigned char *candidateBytes
            = (const unsigned char *)(const void *)candidateSlater;
        const size_t candidateByteCount
            = sizeof(double complex)*candidateCount;
        size_t candidateByte;
        for(candidateByte=0;candidateByte<candidateByteCount;candidateByte++) {
          if(candidateBytes[candidateByte] != candidateSentinel) {
            fprintf(stderr,
                "error: %s: BF-FSZ optimized row path wrote full candidate "
                "at byte %zu\n",caller,candidateByte);
            MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
          }
        }
      }
      RecordBFFSZPfPath(BFFSZ_PROFILE_GREEN,BFFSZ_PF_PATH_OPTIMIZED);
      RecordBFFSZC2Path(detailProfile,BFFSZ_C2_DETAIL_PATH_OPTIMIZED_ROW);
      status = BF_FSZ_PF_OK;
      if(BFFSZPfUpdateCheckEnabled || BFFSZGreenRebuildCheckEnabled
          || BFFSZAffectedCheckEnabled || BFFSZMatrixFreeCheckEnabled
          || forceOracle) {
        int badQp = -1;
        double badError = 0.0;
        double ipError;
        double complex ipUpdated,ipFull;
        BF_FSZ_MaterializeGreenCandidate(
            caller,candidateSlater,debugFull,eleNum,eleIdx,eleSpn,
            oldEleProjBFCnt,newEleProjBFCnt,candidateRows,rowQpStride,
            rowAffectedStride,rowsBuilt,nAffected,affected,
            hopIntWork,hopIntWorkSize,useOMP,detailTimers,
            BFFSZ_GREEN_MATERIALIZE_ORACLE,materializeSeconds);
        RecordBFFSZC2Path(detailProfile,BFFSZ_C2_DETAIL_PATH_DEBUG_ORACLE);
        if(BFFSZPfUpdateCheckEnabled || forceOracle) {
          status = CalculateBF_FSZ_FullCandidatePfM(
              candidateSlater,eleIdx,eleSpn,
              pfMFull,&detail,pfBufM,pfIWork,pfWork,pfRWork);
          if(status != BF_FSZ_PF_OK || !BF_FSZ_GreenPfUpdateMatchesFull(
              pfMNew,pfMFull,NQPFull,&badQp,&badError)) {
            fprintf(stderr,
                "error: %s: BF-FSZ Pfaffian-update/full mismatch "
                "at qp=%d relative_error=%.17e\n",
                caller,badQp,badError);
            MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
          }
          ipUpdated = CalculateIP_fcmp(pfMNew,0,NQPFull,MPI_COMM_SELF);
          ipFull = CalculateIP_fcmp(pfMFull,0,NQPFull,MPI_COMM_SELF);
          ipError = cabs(ipUpdated-ipFull)/fmax(1.0,cabs(ipFull));
          if(!(isfinite(ipError) && ipError <= 1.0e-10)) {
            fprintf(stderr,
                "error: %s: BF-FSZ Pfaffian-update/full IP mismatch "
                "relative_error=%.17e\n",caller,ipError);
            MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
          }
        }
      }
    } else {
      RecordBFFSZPfPath(BFFSZ_PROFILE_GREEN,BFFSZ_PF_PATH_FALLBACK);
      RecordBFFSZC2Path(detailProfile,BFFSZ_C2_DETAIL_PATH_FALLBACK);
      BF_FSZ_MaterializeGreenCandidate(
          caller,candidateSlater,debugFull,eleNum,eleIdx,eleSpn,
          oldEleProjBFCnt,newEleProjBFCnt,candidateRows,rowQpStride,
          rowAffectedStride,rowsBuilt,nAffected,affected,
          hopIntWork,hopIntWorkSize,useOMP,detailTimers,
          BFFSZ_GREEN_MATERIALIZE_FALLBACK,materializeSeconds);
      status = CalculateBF_FSZ_FullCandidatePfM(candidateSlater,eleIdx,eleSpn,
          pfMNew,&detail,pfBufM,pfIWork,pfWork,pfRWork);
    }
  }

  if(status == BF_FSZ_PF_LAPACK_FAILURE) {
    fprintf(stderr,"warning: %s: BF-FSZ candidate Pfaffian failed in M_ZSKPFA (info=%d); propagating NaN.\n",
        caller,detail);
    return BF_FSZ_NaNValue();
  }
  if(status == BF_FSZ_PF_NONFINITE) {
    fprintf(stderr,"warning: %s: BF-FSZ candidate Pfaffian is non-finite (qpidx=%d); propagating NaN.\n",
        caller,detail);
    return BF_FSZ_NaNValue();
  }
  if(status != BF_FSZ_PF_OK) {
    fprintf(stderr,"warning: %s: BF-FSZ candidate Pfaffian returned unknown status %d; propagating NaN.\n",
        caller,status);
    return BF_FSZ_NaNValue();
  }
  ipNew = CalculateIP_fcmp(pfMNew,0,NQPFull,MPI_COMM_SELF);
  if(!(isfinite(creal(ipNew)) && isfinite(cimag(ipNew)))) {
    fprintf(stderr,"warning: %s: BF-FSZ candidate projected overlap is non-finite; propagating NaN.\n",
        caller);
    return BF_FSZ_NaNValue();
  }
  return ipNew;
}
/*
double complex GreenFuncN(const int n, int *rsi, int *rsj, const double complex  ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  double complex *buffer, int *bufferInt);
double complex calculateNewPfMN_child(const int qpidx, const int n, const int *msa, const int *rsa,
                              const int *eleIdx, double complex *buffer);
*/
/* Calculate 1-body Green function <CisAjs> */
/* buffer size = NQPFull */
double complex GreenFunc1_fsz(const int ri, const int rj, const int s, const double complex  ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, complex double *buffer) {
  double complex z;
  int mj,msj,rsi,rsj;
  double complex *pfMNew = buffer; /* NQPFull */

  if(ri==rj) return eleNum[ri+s*Nsite];
  if(eleNum[ri+s*Nsite]==1 || eleNum[rj+s*Nsite]==0) return 0.0;
  if(NExUpdatePath==4 || NExUpdatePath==5){ //For t-J
    if(eleNum[ri+(1-s)*Nsite]==1) return 0.0;
  }

  mj  = eleCfg[rj+s*Nsite];
  msj = mj;// + s*Ne;
  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;

  /* hopping */
  eleIdx[msj] = ri;
  eleSpn[msj] = s;//fsz
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt(rj, ri, s, projCntNew, eleProjCnt, eleNum);
  z = ProjRatio(projCntNew,eleProjCnt);

  /* calculate Pfaffian */
  CalculateNewPfM_fsz(mj, s, pfMNew, eleIdx,eleSpn, 0, NQPFull);//fsz
  z *= CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[msj] = rj;
  eleSpn[msj] = s; //fsz
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  return conj(z/ip);//TBC
}

/* Calculate 1-body Green function <CisAjt> */ // s!=t
/* buffer size = NQPFull */
double complex GreenFunc1_fsz2(const int ri, const int rj, const int s,const int t, const double complex  ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, complex double *buffer) {
  double complex z;
  int mj,rsi,rtj;//msj
  double complex *pfMNew = buffer; /* NQPFull */

  //if(ri==rj) return eleNum[ri+s*Nsite]; //fsz
  if(eleNum[ri+s*Nsite]==1 || eleNum[rj+t*Nsite]==0) return 0.0;
  if(NExUpdatePath==4 || NExUpdatePath==5){ //For t-J
    if(ri!=rj && eleNum[ri+(1-s)*Nsite]==1) return 0.0;
  }

  mj  = eleCfg[rj+t*Nsite];
  //mtj = mj;// + s*Ne;
  rsi = ri + s*Nsite;
  rtj = rj + t*Nsite;

  /* hopping */
  // (j,t) -> (i,s)
  eleIdx[mj] = ri;
  eleSpn[mj] = s;//fsz
  eleNum[rtj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt_fsz(rj, ri, t,s, projCntNew, eleProjCnt, eleNum);
  z = ProjRatio(projCntNew,eleProjCnt);

  /* calculate Pfaffian */
  CalculateNewPfM_fsz(mj, s, pfMNew, eleIdx,eleSpn, 0, NQPFull);//fsz: note EleSpn[mj]=s
  z *= CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[mj] = rj;
  eleSpn[mj] = t; //fsz
  eleNum[rtj] = 1;
  eleNum[rsi] = 0;

  return conj(z/ip);//TBC
}


/* Calculate 2-body Green function <psi|CisAjsCktAlt|x>/<psi|x> */
/* buffer size = NQPFull+2*Nsize */
double complex GreenFunc2_fsz(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double complex *buffer) {
  double complex z;
  int mj,msj,ml,mtl;
  int rsi,rsj,rtk,rtl;
  double complex *pfMNew = buffer; /* [NQPFull] */
  double complex *bufV   = buffer+NQPFull; /* 2*Nsize */

  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;
  rtk = rk + t*Nsite;
  rtl = rl + t*Nsite;

  if(s==t) {
    if(rk==rl) { /* CisAjsNks */
      if(eleNum[rtk]==0) return 0.0;
      else return GreenFunc1_fsz(ri,rj,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CisAjs */
    }else if(rj==rl) {
      return 0.0; /* CisAjsCksAjs (j!=k) */
    }else if(ri==rl) { /* AjsCksNis */
      if(eleNum[rsi]==0) return 0.0;
      else if(rj==rk) return 1.0-eleNum[rsj];
      else return -GreenFunc1_fsz(rk,rj,s,ip,eleIdx,eleCfg,eleNum,
                              eleProjCnt,eleSpn,projCntNew,buffer); /* -CksAjs */
    }else if(rj==rk) { /* CisAls(1-Njs) */
      if(eleNum[rsj]==1) return 0.0;
      else if(ri==rl) return eleNum[rsi];
      else return GreenFunc1_fsz(ri,rl,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CisAls */
    }else if(ri==rk) {
      return 0.0; /* CisAjsCisAls (i!=j) */
    }else if(ri==rj) { /* NisCksAls (i!=k,l) */
      if(eleNum[rsi]==0) return 0.0;
      else return GreenFunc1_fsz(rk,rl,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CksAls */
    }
  }else{
    if(rk==rl) { /* CisAjsNkt */
      if(eleNum[rtk]==0) return 0.0;
      else if(ri==rj) return eleNum[rsi];
      else return GreenFunc1_fsz(ri,rj,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CisAjs */
    }else if(ri==rj) { /* NisCktAlt */
      if(eleNum[rsi]==0) return 0.0;
      else return GreenFunc1_fsz(rk,rl,t,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CktAlt */
    }
  }

  if(eleNum[rsi]==1 || eleNum[rsj]==0 || eleNum[rtk]==1 || eleNum[rtl]==0) return 0.0;

  mj = eleCfg[rj+s*Nsite];
  ml = eleCfg[rl+t*Nsite];
  msj = mj;// + s*Ne;
  mtl = ml;// + t*Ne;

  /* hopping */
  eleIdx[mtl] = rk;
  eleSpn[mtl] = t;
  eleNum[rtl] = 0;
  eleNum[rtk] = 1;
  UpdateProjCnt(rl, rk, t, projCntNew, eleProjCnt, eleNum);
  eleIdx[msj] = ri;
  eleSpn[msj] = s;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt(rj, ri, s, projCntNew, projCntNew, eleNum);

  if(NExUpdatePath==4 || NExUpdatePath==5){ //For t-J
    int doublon_i = eleNum[ri + s*Nsite] *eleNum[ri + (1-s)*Nsite];
    int doublon_k = eleNum[rk + t*Nsite] *eleNum[rk + (1-t)*Nsite];
    if(doublon_i == 1 || doublon_k ==1 ){
      /* revert hopping */
      eleIdx[mtl] = rl;
      eleSpn[mtl] = t;
      eleNum[rtl] = 1;
      eleNum[rtk] = 0;
      eleIdx[msj] = rj;
      eleSpn[msj] = s;
      eleNum[rsj] = 1;
      eleNum[rsi] = 0;
      return 0.0;
    }
  }

  z = ProjRatio(projCntNew,eleProjCnt);

  /* calculate Pfaffian */
  CalculateNewPfMTwo_fsz(ml, t, mj, s, pfMNew, eleIdx,eleSpn, 0, NQPFull, bufV);
  z *= CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[mtl] = rl;
  eleSpn[mtl] = t;
  eleNum[rtl] = 1;
  eleNum[rtk] = 0;
  eleIdx[msj] = rj;
  eleSpn[msj] = s;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  return conj(z/ip);//TBC
}

/* Calculate 2-body Green function <psi|CisAjtCkuAlv|x>/<psi|x> */
/* buffer size = NQPFull+2*Nsize */
double complex GreenFunc2_fsz2(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t,const int u,const int v, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double complex *buffer) {
  double complex z;
  //int mj,msj,ml,mtl;
  int mj,ml; // mj: (rj,t) -> (ri,s) ml:(rl,v) -> (rk,u)
  //int rsi,rtj,ruk,rvl;
  int XI,XJ,XK,XL; //fsz
  double complex *pfMNew = buffer; /* [NQPFull] */
  double complex *bufV   = buffer+NQPFull; /* 2*Nsize */

  XI = ri + s*Nsite;
  XJ = rj + t*Nsite;
  XK = rk + u*Nsite;
  XL = rl + v*Nsite;


  if(XI == XJ) {
    if(XJ == XK) {
      if(XK == XL) {
        // I=J=K=L #1
        return eleNum[XI];
      }else{
        // I=J=K !=L #2
        if(eleNum[XI]==1) return 0.0; //
        else return GreenFunc1_fsz2(rk,rl,u,v,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CkuAlv */
      }
    }else if(XJ == XL){
      // I=J=L !=K #3
      return 0.0;
    }else if(XK == XL){
      // I=J ! K=L #4
      return eleNum[XI]*eleNum[XK];
    }else {
      // I=J   K!=L #5
      if(eleNum[XI]==0) return 0.0; //
      else return GreenFunc1_fsz2(rk,rl,u,v,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer); /* CkuAlv */
    }
  }else if(XI == XK){
    if(XJ == XL){
      // I=K != J=L  #6
      return 0.0;
    }else if(XK == XL){
      // I=K=L !=J #7
      return 0.0;
    }else{
      // I=K != L!=J #8
      return 0.0;
    }
  }else if(XI == XL){
    if(XJ == XK){
      // I=L != J==K #9
      return eleNum[XI]*(1-eleNum[XJ]);
    }else{
      // I=L != J!=K #10
      if(eleNum[XI]==0) return 0.0; //
      else return -1.0*GreenFunc1_fsz2(rk,rj,u,t,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer); /* CkuAjt */
    }
  }else if(XJ == XK){
   if(XK == XL){
     // I != J=K=L  #11
     if(eleNum[XJ]==0) return 0.0; //
     else return GreenFunc1_fsz2(ri,rj,s,t,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer); /* CisAjt */
   }else{
     // I != J=K !=L  #12
     if(eleNum[XJ]==1) return 0.0; //
     else return GreenFunc1_fsz2(ri,rl,s,v,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer); /* CisAlv */
   }
  }else if(XJ == XL){
    // I != J=L !=K  #13
    return 0.0; //
  }else if(XK == XL){
    // I != J != K =L  #14
    if(eleNum[XK]==0) return 0.0; //
    else return GreenFunc1_fsz2(ri,rj,s,t,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer); /* CisAjt */
  }

  // from here, no pair exists
  if(eleNum[XI]==1 || eleNum[XJ]==0 || eleNum[XK]==1 || eleNum[XL]==0) return 0.0;


  mj = eleCfg[XJ]; // electron exists
  ml = eleCfg[XL]; // electron exists
  //msj = mj;// + s*Ne;
  //mtl = ml;// + t*Ne;

  /* hopping */
  eleIdx[ml] = rk; // ml : (rk,u)
  eleSpn[ml] = u;
  eleNum[XL] = 0;
  eleNum[XK] = 1;
  UpdateProjCnt_fsz(rl, rk, v,u, projCntNew, eleProjCnt, eleNum); // (rl,v) -> (rk,u) v!=u
  eleIdx[mj]  = ri; // mj : (ri,s)
  eleSpn[mj]  = s;
  eleNum[XJ] = 0;
  eleNum[XI] = 1;
  UpdateProjCnt_fsz(rj,ri, t,s, projCntNew, projCntNew, eleNum); // (rj,t) -> (ri,s) t!=s

  if(NExUpdatePath==4 || NExUpdatePath==5){
    int doublon_i = eleNum[ri + s*Nsite] *eleNum[ri + (1-s)*Nsite];
    int doublon_k = eleNum[rk + u*Nsite] *eleNum[rk + (1-u)*Nsite];
    if(doublon_i == 1 || doublon_k ==1 ){
      /* revert hopping */
      eleIdx[ml]  = rl;  // ml : (rl,v)
      eleSpn[ml]  = v;
      eleNum[XL] = 1;
      eleNum[XK] = 0;
      eleIdx[mj]  = rj;  // mj : (rj,t)
      eleSpn[mj]  = t;
      eleNum[XJ] = 1;
      eleNum[XI] = 0;
      return 0.0;
    }
  }

  z = ProjRatio(projCntNew,eleProjCnt);

  /* calculate Pfaffian */
  CalculateNewPfMTwo_fsz(ml, u, mj, s, pfMNew, eleIdx,eleSpn, 0, NQPFull, bufV); // ml -> rk,u; mj -> ri,s
  z *= CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[ml]  = rl;  // ml : (rl,v)
  eleSpn[ml]  = v;
  eleNum[XL] = 1;
  eleNum[XK] = 0;
  eleIdx[mj]  = rj;  // mj : (rj,t)
  eleSpn[mj]  = t;
  eleNum[XJ] = 1;
  eleNum[XI] = 0;

  return conj(z/ip);//TBC
}

double complex GreenFunc1BF_fsz(const int ri, const int rj, const int s, const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork) {
  return GreenFunc1BF_fsz2(ri,rj,s,s,ip,eleIdx,eleCfg,eleNum,eleProjCnt,
      eleSpn,projCntNew,eleProjBFCnt,projBFCntNew,buffer,affected,
      hopIntWork,hopIntWorkSize,pfIWork,pfRWork,pfBufM,pfWork);
}

double complex GreenFunc1BF_fsz2(const int ri, const int rj, const int s, const int t,
                  const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork) {
  double complex z;
  int mj,rsi,rtj,nChanged,nAffected,rowStatus,rowsBuilt;
  size_t pfUpdateComplexCount,pfUpdateIntCount,pfUpdateDoubleCount;
  const size_t bfSlaterSize = (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
  const int rowCapacity = (BFFSZPfUpdateKFull-1 < Nsize-1)
      ? BFFSZPfUpdateKFull-1 : Nsize-1;
  const size_t rowAffectedStride = (size_t)Nsize;
  const size_t rowQpStride = (size_t)rowCapacity*rowAffectedStride;
  size_t rowCount;
  /* buffer: new/full Pfaffians, full fallback Slater, candidate rows,
     row-update tail workspace, optional full-rebuild oracle. */
  double complex *pfMNew = buffer;
  double complex *pfMFull = pfMNew + NQPFull;
  double complex *sltElmBFNew = pfMFull + NQPFull;
  double complex *candidateRows = sltElmBFNew + bfSlaterSize;
  double complex *pfUpdateWork;
  double complex *debugFull;

  rsi = ri+s*Nsite;
  rtj = rj+t*Nsite;
  if(rsi==rtj) return eleNum[rsi];
  if(eleNum[rsi]==1 || eleNum[rtj]==0) return 0.0;
  if(NExUpdatePath==4 || NExUpdatePath==5){
    if(ri!=rj && eleNum[ri+(1-s)*Nsite]==1) return 0.0;
  }

  mj  = eleCfg[rtj];

  eleCfg[rtj] = -1;
  eleCfg[rsi] = mj;
  eleIdx[mj] = ri;
  eleSpn[mj] = s;
  eleNum[rtj] = 0;
  eleNum[rsi] = 1;
  if(s == t) {
    UpdateProjCnt(rj, ri, s, projCntNew, eleProjCnt, eleNum);
  } else {
    UpdateProjCnt_fsz(rj,ri,t,s,projCntNew,eleProjCnt,eleNum);
  }
  z = ProjRatio(projCntNew,eleProjCnt);

  StartTimer(81);
  MakeProjBFCnt(projBFCntNew, eleNum);
  StopTimer(81);
  StartTimer(82);
  if(GetSlaterElmBF_fsz_hop_row_work_size(
      &rowCount,NQPFull,rowCapacity) != BF_FSZ_ROW_BUILD_OK) {
    fprintf(stderr, "error: GreenFunc1BF_fsz invalid row workspace size\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  pfUpdateWork = candidateRows+rowCount;
  rowStatus = MakeSlaterElmBF_fsz_hop_rows_workspace(
      candidateRows,rowQpStride,rowAffectedStride,rowCapacity,
      SlaterElmBF,mj,eleIdx,eleSpn,
      eleProjBFCnt,projBFCntNew,0,NQPFull,
      affected,&nChanged,&nAffected,
      hopIntWork, hopIntWorkSize);
  if(rowStatus == BF_FSZ_ROW_BUILD_INVALID_ARGUMENT) {
    fprintf(stderr, "error: GreenFunc1BF_fsz candidate-row construction failed\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  rowsBuilt = (rowStatus == BF_FSZ_ROW_BUILD_OK);
  RecordBFFSZProfile(BFFSZ_PROFILE_GREEN,nChanged,nAffected);
  StopTimer(82);
  StartTimer(83);
  if(GetCalculateNewPfMBF_fsz_row_values_work_size(&pfUpdateComplexCount,
      &pfUpdateIntCount,&pfUpdateDoubleCount) != 0
      || pfIWork == NULL || pfRWork == NULL || pfBufM == NULL
      || pfWork == NULL) {
    fprintf(stderr,"error: invalid GreenFunc1BF_fsz Pfaffian workspace\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  debugFull = BFFSZGreenRebuildCheckEnabled
      ? pfUpdateWork+pfUpdateComplexCount : NULL;
  z *= CalculateBF_FSZ_UpdatedCandidateIP(
      "GreenFunc1BF_fsz2",candidateRows,rowQpStride,rowAffectedStride,
      rowsBuilt,sltElmBFNew,debugFull,eleIdx,eleSpn,eleNum,
      eleProjBFCnt,projBFCntNew,nAffected,affected,pfMNew,pfMFull,pfUpdateWork,
      pfUpdateComplexCount,pfIWork,pfUpdateIntCount,pfRWork,
      pfUpdateDoubleCount,pfBufM,pfIWork,pfWork,pfRWork,
      hopIntWork,hopIntWorkSize,1,1,0,NULL,NULL);
  StopTimer(83);

  eleCfg[rtj] = mj;
  eleCfg[rsi] = -1;
  eleIdx[mj] = rj;
  eleSpn[mj] = t;
  eleNum[rtj] = 1;
  eleNum[rsi] = 0;

  return conj(z/ip);
}

double complex GreenFunc1BF_fsz_workspace(const int ri, const int rj, const int s,
                  const double complex ip, int *eleIdx, int *eleCfg, int *eleNum,
                  const int *eleProjCnt, int *eleSpn, int *projCntNew,
                  const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, double complex *pfBufM,
                  int *affected, int *hopIntWork, int hopIntWorkSize, int *pfIWork,
                  double complex *pfWork, double *pfRWork) {
  return GreenFunc1BF_fsz2_workspace(ri,rj,s,s,ip,eleIdx,eleCfg,eleNum,
      eleProjCnt,eleSpn,projCntNew,eleProjBFCnt,projBFCntNew,buffer,pfBufM,
      affected,hopIntWork,hopIntWorkSize,pfIWork,pfWork,pfRWork);
}

double complex GreenFunc1BF_fsz2_workspace(const int ri, const int rj,
                  const int s, const int t, const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum,
                  const int *eleProjCnt, int *eleSpn, int *projCntNew,
                  const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, double complex *pfBufM,
                  int *affected, int *hopIntWork, int hopIntWorkSize, int *pfIWork,
                  double complex *pfWork, double *pfRWork) {
  double complex z;
  int mj,rsi,rtj,nChanged,nAffected,rowStatus,rowsBuilt;
  size_t pfUpdateComplexCount,pfUpdateIntCount,pfUpdateDoubleCount;
  const size_t bfSlaterSize = (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
  const int rowCapacity = (BFFSZPfUpdateKFull-1 < Nsize-1)
      ? BFFSZPfUpdateKFull-1 : Nsize-1;
  const size_t rowAffectedStride = (size_t)Nsize;
  const size_t rowQpStride = (size_t)rowCapacity*rowAffectedStride;
  size_t rowCount;
  /* buffer: new/full Pfaffians, full fallback Slater, candidate rows,
     row-update tail workspace, optional full-rebuild oracle. */
  double complex *pfMNew = buffer;
  double complex *pfMFull = pfMNew + NQPFull;
  double complex *sltElmBFNew = pfMFull + NQPFull;
  double complex *candidateRows = sltElmBFNew + bfSlaterSize;
  double complex *pfUpdateWork;
  double complex *debugFull;

  rsi = ri+s*Nsite;
  rtj = rj+t*Nsite;
  if(rsi==rtj) return eleNum[rsi];
  if(eleNum[rsi]==1 || eleNum[rtj]==0) return 0.0;
  if(NExUpdatePath==4 || NExUpdatePath==5){
    if(ri!=rj && eleNum[ri+(1-s)*Nsite]==1) return 0.0;
  }

  mj  = eleCfg[rtj];

  eleCfg[rtj] = -1;
  eleCfg[rsi] = mj;
  eleIdx[mj] = ri;
  eleSpn[mj] = s;
  eleNum[rtj] = 0;
  eleNum[rsi] = 1;
  if(s == t) {
    UpdateProjCnt(rj,ri,s,projCntNew,eleProjCnt,eleNum);
  } else {
    UpdateProjCnt_fsz(rj,ri,t,s,projCntNew,eleProjCnt,eleNum);
  }
  z = ProjRatio(projCntNew,eleProjCnt);

  MakeProjBFCnt(projBFCntNew, eleNum);
  if(hopIntWork == NULL
      || BF_FSZ_IntWorkspacesOverlap(affected, Nsize, pfIWork, Nsize)) {
    fprintf(stderr, "error: invalid BF-FSZ Green workspace alias\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(GetSlaterElmBF_fsz_hop_row_work_size(
      &rowCount,NQPFull,rowCapacity) != BF_FSZ_ROW_BUILD_OK) {
    fprintf(stderr, "error: GreenFunc1BF_fsz_workspace invalid row workspace size\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  pfUpdateWork = candidateRows+rowCount;
  rowStatus = MakeSlaterElmBF_fsz_hop_rows_workspace_serial(
      candidateRows,rowQpStride,rowAffectedStride,rowCapacity,
      SlaterElmBF,mj,eleIdx,eleSpn,
      eleProjBFCnt,projBFCntNew,0,NQPFull,
      affected,&nChanged,&nAffected,
      hopIntWork, hopIntWorkSize);
  if(rowStatus == BF_FSZ_ROW_BUILD_INVALID_ARGUMENT) {
    fprintf(stderr, "error: GreenFunc1BF_fsz_workspace candidate-row construction failed\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  rowsBuilt = (rowStatus == BF_FSZ_ROW_BUILD_OK);
  RecordBFFSZProfile(BFFSZ_PROFILE_GREEN,nChanged,nAffected);
  if(GetCalculateNewPfMBF_fsz_row_values_work_size(&pfUpdateComplexCount,
      &pfUpdateIntCount,&pfUpdateDoubleCount) != 0) {
    fprintf(stderr,"error: invalid GreenFunc1BF_fsz_workspace Pfaffian workspace\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  debugFull = BFFSZGreenRebuildCheckEnabled
      ? pfUpdateWork+pfUpdateComplexCount : NULL;
  z *= CalculateBF_FSZ_UpdatedCandidateIP(
      "GreenFunc1BF_fsz2",candidateRows,rowQpStride,rowAffectedStride,
      rowsBuilt,sltElmBFNew,debugFull,eleIdx,eleSpn,eleNum,
      eleProjBFCnt,projBFCntNew,nAffected,affected,pfMNew,pfMFull,pfUpdateWork,
      pfUpdateComplexCount,pfIWork,pfUpdateIntCount,pfRWork,
      pfUpdateDoubleCount,pfBufM,pfIWork,pfWork,pfRWork,
      hopIntWork,hopIntWorkSize,0,0,0,NULL,NULL);

  eleCfg[rtj] = mj;
  eleCfg[rsi] = -1;
  eleIdx[mj] = rj;
  eleSpn[mj] = t;
  eleNum[rtj] = 1;
  eleNum[rsi] = 0;

  return conj(z/ip);
}

double complex GreenFunc2BF_fsz(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork) {
  return GreenFunc2BF_fsz2(ri,rj,rk,rl,s,s,t,t,ip,
      eleIdx,eleCfg,eleNum,eleProjCnt,eleSpn,projCntNew,
      eleProjBFCnt,projBFCntNew,buffer,affected,hopIntWork,
      hopIntWorkSize,pfIWork,pfRWork,pfBufM,pfWork);
}

double complex GreenFunc2BF_fsz2WithProfile(
                  const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const int u, const int v,
                  const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork,
                  BFFSZC2DetailContext *detailProfile) {
  double complex z,result,ipNew;
  int mj,ml,nChanged,nAffected,rowStatus,rowsBuilt;
  int profileNChanged=0,profileNAffected=0;
  int movedParticles[2];
  int c2Class = 0;
  int stateAllowed;
  int XI,XJ,XK,XL;
  int *stateSnapshot = NULL;
  double dispatchStart = 0.0;
  double componentStart = 0.0;
  double materializeSeconds = 0.0;
  double pfaffianSeconds = 0.0;
  size_t pfUpdateComplexCount,pfUpdateIntCount,pfUpdateDoubleCount,rowCount;
  const size_t bfSlaterSize
      = (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
  const int rowCapacity = (BFFSZPfUpdateKFull-1 < Nsize-1)
      ? BFFSZPfUpdateKFull-1 : Nsize-1;
  const size_t rowAffectedStride = (size_t)Nsize;
  const size_t rowQpStride = (size_t)rowCapacity*rowAffectedStride;
  double complex *pfMNew = buffer;
  double complex *pfMFull = pfMNew+NQPFull;
  double complex *candidateSlater = buffer+2*(size_t)NQPFull;
  double complex *candidateRows = candidateSlater+bfSlaterSize;
  double complex *pfUpdateWork;
  double complex *debugFull;

  XI = ri+s*Nsite;
  XJ = rj+t*Nsite;
  XK = rk+u*Nsite;
  XL = rl+v*Nsite;

  if(detailProfile != NULL) {
    dispatchStart = BFFSZC2DetailMonotonicSeconds();
    c2Class = ClassifyGreenFunc2FSZ(XI,XJ,XK,XL);
    detailProfile->classCall[c2Class-1]++;
  }

  if(XI == XJ) {
    if(XJ == XK) {
      if(XK == XL) {
        RecordBFFSZC2DispatchSeconds(detailProfile,dispatchStart);
        return eleNum[XI];
      } else {
        RecordBFFSZC2DispatchSeconds(detailProfile,dispatchStart);
        if(eleNum[XI] == 1) return 0.0;
        return GreenFunc1BF_fsz2(rk,rl,u,v,ip,eleIdx,eleCfg,eleNum,
            eleProjCnt,eleSpn,projCntNew,eleProjBFCnt,projBFCntNew,
            buffer,affected,hopIntWork,hopIntWorkSize,
            pfIWork,pfRWork,pfBufM,pfWork);
      }
    } else if(XJ == XL) {
      RecordBFFSZC2DispatchSeconds(detailProfile,dispatchStart);
      return 0.0;
    } else if(XK == XL) {
      RecordBFFSZC2DispatchSeconds(detailProfile,dispatchStart);
      return eleNum[XI]*eleNum[XK];
    } else {
      RecordBFFSZC2DispatchSeconds(detailProfile,dispatchStart);
      if(eleNum[XI] == 0) return 0.0;
      return GreenFunc1BF_fsz2(rk,rl,u,v,ip,eleIdx,eleCfg,eleNum,
          eleProjCnt,eleSpn,projCntNew,eleProjBFCnt,projBFCntNew,
          buffer,affected,hopIntWork,hopIntWorkSize,
          pfIWork,pfRWork,pfBufM,pfWork);
    }
  } else if(XI == XK) {
    RecordBFFSZC2DispatchSeconds(detailProfile,dispatchStart);
    return 0.0;
  } else if(XI == XL) {
    if(XJ == XK) {
      RecordBFFSZC2DispatchSeconds(detailProfile,dispatchStart);
      return eleNum[XI]*(1-eleNum[XJ]);
    }
    RecordBFFSZC2DispatchSeconds(detailProfile,dispatchStart);
    if(eleNum[XI] == 0) return 0.0;
    return -GreenFunc1BF_fsz2(rk,rj,u,t,ip,eleIdx,eleCfg,eleNum,
        eleProjCnt,eleSpn,projCntNew,eleProjBFCnt,projBFCntNew,
        buffer,affected,hopIntWork,hopIntWorkSize,
        pfIWork,pfRWork,pfBufM,pfWork);
  } else if(XJ == XK) {
    if(XK == XL) {
      RecordBFFSZC2DispatchSeconds(detailProfile,dispatchStart);
      if(eleNum[XJ] == 0) return 0.0;
      return GreenFunc1BF_fsz2(ri,rj,s,t,ip,eleIdx,eleCfg,eleNum,
          eleProjCnt,eleSpn,projCntNew,eleProjBFCnt,projBFCntNew,
          buffer,affected,hopIntWork,hopIntWorkSize,
          pfIWork,pfRWork,pfBufM,pfWork);
    }
    RecordBFFSZC2DispatchSeconds(detailProfile,dispatchStart);
    if(eleNum[XJ] == 1) return 0.0;
    return GreenFunc1BF_fsz2(ri,rl,s,v,ip,eleIdx,eleCfg,eleNum,
        eleProjCnt,eleSpn,projCntNew,eleProjBFCnt,projBFCntNew,
        buffer,affected,hopIntWork,hopIntWorkSize,
        pfIWork,pfRWork,pfBufM,pfWork);
  } else if(XJ == XL) {
    RecordBFFSZC2DispatchSeconds(detailProfile,dispatchStart);
    return 0.0;
  } else if(XK == XL) {
    RecordBFFSZC2DispatchSeconds(detailProfile,dispatchStart);
    if(eleNum[XK] == 0) return 0.0;
    return GreenFunc1BF_fsz2(ri,rj,s,t,ip,eleIdx,eleCfg,eleNum,
        eleProjCnt,eleSpn,projCntNew,eleProjBFCnt,projBFCntNew,
        buffer,affected,hopIntWork,hopIntWorkSize,
        pfIWork,pfRWork,pfBufM,pfWork);
  }

  if(BFFSZC2StateCheckEnabled) {
    stateSnapshot = BF_FSZ_CaptureGreen2State(eleIdx,eleSpn,eleCfg,eleNum);
    if(stateSnapshot == NULL) {
      fprintf(stderr,"error: GreenFunc2BF_fsz2 state snapshot allocation failed\n");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
  }

  if(eleNum[XI] == 1 || eleNum[XJ] == 0
      || eleNum[XK] == 1 || eleNum[XL] == 0) {
    if(detailProfile != NULL) {
      RecordBFFSZC2DispatchSeconds(detailProfile,dispatchStart);
      detailProfile->orderedDescriptorTotal++;
      detailProfile->outcome[BFFSZ_C2_DETAIL_OUTCOME_OCCUPANCY_ZERO]++;
    }
    if(stateSnapshot != NULL && !BF_FSZ_Green2StateMatches(
        stateSnapshot,eleIdx,eleSpn,eleCfg,eleNum)) {
      fprintf(stderr,"error: GreenFunc2BF_fsz2 occupancy-zero changed local state\n");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    free(stateSnapshot);
    return 0.0;
  }
  if(detailProfile != NULL) {
    RecordBFFSZC2DispatchSeconds(detailProfile,dispatchStart);
    detailProfile->orderedDescriptorTotal++;
  }

  mj = eleCfg[XJ];
  ml = eleCfg[XL];

  if(detailProfile != NULL) componentStart = BFFSZC2DetailMonotonicSeconds();

  eleCfg[XL] = -1;
  eleCfg[XK] = ml;
  eleIdx[ml] = rk;
  eleSpn[ml] = u;
  eleNum[XL] = 0;
  eleNum[XK] = 1;
  UpdateProjCnt_fsz(rl,rk,v,u,projCntNew,eleProjCnt,eleNum);

  eleCfg[XJ] = -1;
  eleCfg[XI] = mj;
  eleIdx[mj] = ri;
  eleSpn[mj] = s;
  eleNum[XJ] = 0;
  eleNum[XI] = 1;
  UpdateProjCnt_fsz(rj,ri,t,s,projCntNew,projCntNew,eleNum);

  result = 0.0;
  stateAllowed = IsSectorStateAllowed(eleNum) && !BFFSZC2ForceSectorZero;
  if(stateAllowed) {
    z = ProjRatio(projCntNew,eleProjCnt);
  }
  if(detailProfile != NULL) {
    detailProfile->componentSeconds[BFFSZ_C2_DETAIL_COMPONENT_STATE_PROJECTION]
        += BFFSZC2DetailMonotonicSeconds()-componentStart;
  }
  if(stateAllowed) {
    if(detailProfile != NULL) componentStart = BFFSZC2DetailMonotonicSeconds();
    MakeProjBFCnt(projBFCntNew,eleNum);
    if(detailProfile != NULL) {
      detailProfile->componentSeconds[BFFSZ_C2_DETAIL_COMPONENT_BF_COUNT]
          += BFFSZC2DetailMonotonicSeconds()-componentStart;
      componentStart = BFFSZC2DetailMonotonicSeconds();
      if(CollectBFAffectedParticlesTwoMove_fsz(
          eleProjBFCnt,projBFCntNew,ml,mj,eleIdx,eleSpn,
          affected,hopIntWork,hopIntWorkSize,
          &profileNChanged,&profileNAffected) != 0) {
        fprintf(stderr,"error: BF-FSZ C2 detail affected collector failed\n");
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
      }
      detailProfile->componentSeconds[BFFSZ_C2_DETAIL_COMPONENT_AFFECTED_COLLECT]
          += BFFSZC2DetailMonotonicSeconds()-componentStart;
      componentStart = BFFSZC2DetailMonotonicSeconds();
    }
    if(GetSlaterElmBF_fsz_hop_row_work_size(
          &rowCount,NQPFull,rowCapacity) != BF_FSZ_ROW_BUILD_OK
        || GetCalculateNewPfMBF_fsz_row_values_work_size(
          &pfUpdateComplexCount,&pfUpdateIntCount,&pfUpdateDoubleCount) != 0
        || pfIWork == NULL || pfRWork == NULL || pfBufM == NULL
        || pfWork == NULL) {
      fprintf(stderr,"error: invalid GreenFunc2BF_fsz Pfaffian workspace\n");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    movedParticles[0] = ml;
    movedParticles[1] = mj;
    pfUpdateWork = candidateRows+rowCount;
    debugFull = BFFSZGreenRebuildCheckEnabled
        ? pfUpdateWork+pfUpdateComplexCount : NULL;
    rowStatus = MakeSlaterElmBF_fsz_multi_move_rows_workspace_serial(
        candidateRows,rowQpStride,rowAffectedStride,rowCapacity,
        SlaterElmBF,2,movedParticles,eleIdx,eleSpn,
        eleProjBFCnt,projBFCntNew,0,NQPFull,
        affected,&nChanged,&nAffected,hopIntWork,hopIntWorkSize);
    if(rowStatus == BF_FSZ_ROW_BUILD_INVALID_ARGUMENT) {
      fprintf(stderr,"error: GreenFunc2BF_fsz candidate-row construction failed\n");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    rowsBuilt = (rowStatus == BF_FSZ_ROW_BUILD_OK);
    if(detailProfile != NULL) {
      if(nChanged != profileNChanged || nAffected != profileNAffected) {
        fprintf(stderr,"error: BF-FSZ C2 detail/row affected mismatch\n");
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
      }
      detailProfile->evaluatedCalls++;
      detailProfile->changedSum += nChanged;
      detailProfile->affectedSum += nAffected;
      if(nChanged > detailProfile->changedMax) detailProfile->changedMax = nChanged;
      if(nAffected > detailProfile->affectedMax) detailProfile->affectedMax = nAffected;
      detailProfile->changedHist[BFFSZProfileHistBin(nChanged)]++;
      detailProfile->affectedHist[BFFSZProfileHistBin(nAffected)]++;
      if(nAffected >= BFFSZPfUpdateKFull) detailProfile->affectedAtOrAboveKFull++;
      detailProfile->outcome[BFFSZ_C2_DETAIL_OUTCOME_EVALUATED]++;
      detailProfile->componentSeconds[BFFSZ_C2_DETAIL_COMPONENT_CANDIDATE_BUILD]
          += BFFSZC2DetailMonotonicSeconds()-componentStart;
      if(detailProfile->reuseCensus != NULL
          && RecordBFFSZC2ReuseCensus(
              detailProfile->reuseCensus,XI,XJ,XK,XL) != 0) {
        fprintf(stderr,"error: BF-FSZ C2 reuse census insertion failed\n");
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
      }
    }
    RecordBFFSZProfile(BFFSZ_PROFILE_GREEN,nChanged,nAffected);
    if(detailProfile != NULL) componentStart = BFFSZC2DetailMonotonicSeconds();
    ipNew = CalculateBF_FSZ_UpdatedCandidateIP(
        "GreenFunc2BF_fsz2",candidateRows,rowQpStride,rowAffectedStride,
        rowsBuilt,candidateSlater,debugFull,eleIdx,eleSpn,eleNum,
        eleProjBFCnt,projBFCntNew,nAffected,affected,pfMNew,pfMFull,
        pfUpdateWork,pfUpdateComplexCount,pfIWork,pfUpdateIntCount,
        pfRWork,pfUpdateDoubleCount,pfBufM,pfIWork,pfWork,pfRWork,
        hopIntWork,hopIntWorkSize,0,0,BFFSZMultiMoveRowCheckEnabled,
        detailProfile,(detailProfile != NULL) ? &materializeSeconds : NULL);
    if(BFFSZMultiMoveRowCheckEnabled) {
      BF_FSZ_CheckMultiMoveRowsOnce(
          ml,mj,eleIdx,eleSpn,eleProjBFCnt,projBFCntNew,
          candidateSlater,hopIntWork,hopIntWorkSize);
    }
    if(detailProfile != NULL) {
      pfaffianSeconds = BFFSZC2DetailMonotonicSeconds()-componentStart;
      if(materializeSeconds > pfaffianSeconds) {
        materializeSeconds = pfaffianSeconds;
      }
      detailProfile->componentSeconds[BFFSZ_C2_DETAIL_COMPONENT_CANDIDATE_BUILD]
          += materializeSeconds;
      detailProfile->componentSeconds[BFFSZ_C2_DETAIL_COMPONENT_PFAFFIAN]
          += pfaffianSeconds-materializeSeconds;
    }
    z *= ipNew;
    result = conj(z/ip);
  } else if(detailProfile != NULL) {
    detailProfile->outcome[BFFSZ_C2_DETAIL_OUTCOME_SECTOR_ZERO]++;
  }

  if(detailProfile != NULL) componentStart = BFFSZC2DetailMonotonicSeconds();
  eleCfg[XL] = ml;
  eleCfg[XK] = -1;
  eleIdx[ml] = rl;
  eleSpn[ml] = v;
  eleNum[XL] = 1;
  eleNum[XK] = 0;
  eleCfg[XJ] = mj;
  eleCfg[XI] = -1;
  eleIdx[mj] = rj;
  eleSpn[mj] = t;
  eleNum[XJ] = 1;
  eleNum[XI] = 0;
  if(detailProfile != NULL) {
    detailProfile->componentSeconds[BFFSZ_C2_DETAIL_COMPONENT_RESTORE]
        += BFFSZC2DetailMonotonicSeconds()-componentStart;
  }

  if(stateSnapshot != NULL && !BF_FSZ_Green2StateMatches(
      stateSnapshot,eleIdx,eleSpn,eleCfg,eleNum)) {
    fprintf(stderr,"error: GreenFunc2BF_fsz2 failed to restore local state\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  free(stateSnapshot);

  return result;
}

double complex GreenFunc2BF_fsz2(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const int u, const int v,
                  const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork) {
  return GreenFunc2BF_fsz2WithProfile(ri,rj,rk,rl,s,t,u,v,ip,
      eleIdx,eleCfg,eleNum,eleProjCnt,eleSpn,projCntNew,
      eleProjBFCnt,projBFCntNew,buffer,affected,hopIntWork,hopIntWorkSize,
      pfIWork,pfRWork,pfBufM,pfWork,NULL);
}

/* Calculate n-body Green function in orbital-general (fsz) mode. */
double complex GreenFuncN_fsz(const int n, int *rsi, int *rsj, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  int *eleSpn, double complex *buffer, int *bufferInt, double *rwork) {
  int ri,rj,rk,rl,si,sj,sk,sl;
  int k,l,m,rsk;
  double complex z,x;
  int qpidx;

  int *projCntNew = bufferInt; /* [NProj] */

  int msj[n];
  int oldSite[n], oldSpin[n], newSite[n], newSpin[n];
  double complex *pfMNew = buffer; /* [NQPFull] */
  double complex *bufV = buffer+NQPFull; /* [n*Nsize] */

  if(n<=0) return 0;
  else if(n==1) {
    ri = rsi[0]%Nsite;
    si = rsi[0]/Nsite;
    rj = rsj[0]%Nsite;
    sj = rsj[0]/Nsite;
    if(si==sj) {
      return GreenFunc1_fsz(ri,rj,si,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer);
    }
    return GreenFunc1_fsz2(ri,rj,si,sj,ip,eleIdx,eleCfg,eleNum,
                           eleProjCnt,eleSpn,projCntNew,buffer);
  } else if(n==2) {
    ri = rsi[0]%Nsite;
    si = rsi[0]/Nsite;
    rj = rsj[0]%Nsite;
    sj = rsj[0]/Nsite;
    rk = rsi[1]%Nsite;
    sk = rsi[1]/Nsite;
    rl = rsj[1]%Nsite;
    sl = rsj[1]/Nsite;
    if(si==sj && sk==sl) {
      return GreenFunc2_fsz(ri,rj,rk,rl,si,sk,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer);
    }
    return GreenFunc2_fsz2(ri,rj,rk,rl,si,sj,sk,sl,ip,eleIdx,eleCfg,eleNum,
                           eleProjCnt,eleSpn,projCntNew,buffer);
  }

  /* reduction */
  for(k=n-1;k>=0;k--) {
    /* ** check for an annihilation operator at rsj[k] ** */
    rsk = rsj[k];
    for(l=k+1;l<n;l++) {
      /* rsj[k] == rsi[l] */
      if(rsk==rsi[l]) {
        rsj[k] = rsj[l];
        for(m=l;m<n-1;m++) { /* shift */
          rsi[m] = rsi[m+1];
          rsj[m] = rsj[m+1];
        }
        return GreenFuncN_fsz(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,
                              eleProjCnt,eleSpn,buffer,bufferInt,rwork);
      }
      /* rsj[k] == rsj[l] */
      if(rsk==rsj[l]) return 0;
    }
    /* check electron number */
    if(eleNum[rsk]==0) return 0;

    /* ** check for a creation operator at rsi[k] ** */
    rsk = rsi[k];
    /* rsi[k] == rsj[k] */
    if(rsk==rsj[k]) {
      for(m=k;m<n-1;m++) { /* shift */
        rsi[m] = rsi[m+1];
        rsj[m] = rsj[m+1];
      }
      return GreenFuncN_fsz(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,buffer,bufferInt,rwork);
    }
    for(l=k+1;l<n;l++) {
      /* rsi[k] == rsi[l] */
      if(rsk==rsi[l]) return 0;
      /* rsi[k] == rsj[l] (k<l) */
      if(rsk==rsj[l]) {
        rsi[k] = rsi[l];
        for(m=l;m<n-1;m++) { /* shift */
          rsi[m] = rsi[m+1];
          rsj[m] = rsj[m+1];
        }
        return (-1.0)*GreenFuncN_fsz(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,
                                     eleProjCnt,eleSpn,buffer,bufferInt,rwork);
      }
    }
    /* check electron number */
    if(eleNum[rsk]==1) return 0;
  }

  /* hopping */
  #pragma loop noalias
  for(k=0;k<n;k++) {
    oldSite[k] = rsj[k]%Nsite;
    oldSpin[k] = rsj[k]/Nsite;
    newSite[k] = rsi[k]%Nsite;
    newSpin[k] = rsi[k]/Nsite;
    msj[k] = eleCfg[rsj[k]];

    eleIdx[msj[k]] = newSite[k];
    eleSpn[msj[k]] = newSpin[k];
    eleNum[rsj[k]] = 0;
    eleNum[rsi[k]] = 1;
    if(k==0) {
      UpdateProjCnt_fsz(oldSite[k], newSite[k], oldSpin[k], newSpin[k],
                        projCntNew, eleProjCnt, eleNum);
    } else {
      UpdateProjCnt_fsz(oldSite[k], newSite[k], oldSpin[k], newSpin[k],
                        projCntNew, projCntNew, eleNum);
    }
  }

  if(!IsSectorStateAllowed(eleNum)) {
    #pragma loop noalias
    for(k=0;k<n;k++) {
      eleIdx[msj[k]] = oldSite[k];
      eleSpn[msj[k]] = oldSpin[k];
      eleNum[rsj[k]] = 1;
      eleNum[rsi[k]] = 0;
    }
    return 0.0;
  }

  z = ProjRatio(projCntNew,eleProjCnt);

  /* calculateNewPfM */
  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    pfMNew[qpidx] = calculateNewPfMN_child_fsz(qpidx,n,msj,eleIdx,eleSpn,bufV,rwork);
  }
  x = CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  #pragma loop noalias
  for(k=0;k<n;k++) {
    eleIdx[msj[k]] = oldSite[k];
    eleSpn[msj[k]] = oldSpin[k];
    eleNum[rsj[k]] = 1;
    eleNum[rsi[k]] = 0;
  }

  return conj(z*x/ip);
}

/* msa[k]-th electron has already hopped to (eleIdx[msa[k]], eleSpn[msa[k]]). */
double complex calculateNewPfMN_child_fsz(const int qpidx, const int n, const int *msa,
                              const int *eleIdx, const int *eleSpn,
                              double complex *buffer, double *rwork) {
  const int nsize = Nsize;
  const int n2 = 2*n;
  const double complex *sltE;
  const double complex *sltE_k;
  const double complex *invM;
  const double complex *invM_i, *invM_k, *invM_l;

  double complex *vec; /* vec[n][nsize] */
  double complex *vec_k, *vec_l;
  double complex mat[n2*n2]; /* mat[n2][n2] */
  double complex *mat_k;
  double sgn;

  int rsi,rsk,msi,msj,k,l;
  double complex val,tmp;

  char uplo='U', mthd='P';
  int lda,info=0;
  int nn;
  double complex pfaff;
  int iwork[n2];
  double complex work[n2*n2]; /* [n2][n2] */
  int lwork = n2*n2;

  nn=lda=n2;

  sltE = SlaterElm + qpidx*Nsite2*Nsite2;
  invM = InvM + qpidx*Nsize*Nsize;

  vec = buffer; /* n*nsize */

  #pragma loop noalias
  for(k=0;k<n;k++) {
    rsk = eleIdx[msa[k]] + eleSpn[msa[k]]*Nsite;
    sltE_k = sltE + rsk*Nsite2;
    vec_k = vec + k*nsize;
    #pragma loop norecurrence
    for(msi=0;msi<nsize;msi++) {
      rsi = eleIdx[msi] + eleSpn[msi]*Nsite;
      vec_k[msi] = sltE_k[rsi];
    }
  }

  /* X_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k;
    vec_k = vec + k*nsize;
    for(l=k+1;l<n;l++) {
      vec_l = vec + l*nsize;
      val = 0.0;
      for(msi=0;msi<nsize;msi++) {
        invM_i = invM + msi*nsize;
        tmp = 0.0;
        for(msj=0;msj<nsize;msj++) {
          tmp += invM_i[msj] * vec_l[msj];
        }
        val += tmp * vec_k[msi];
      }
      mat_k[l] = val + vec_k[msa[l]];
    }
  }

  /* Y_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k + n;
    vec_k = vec + k*nsize;
    for(l=0;l<n;l++) {
      invM_l = invM + msa[l]*nsize;
      val = 0.0;
      for(msi=0;msi<nsize;msi++) {
        val += vec_k[msi] * invM_l[msi];
      }
      mat_k[l] = val;
    }
  }

  /* Z_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*(k+n) + n;
    invM_k = invM + msa[k]*nsize;
    for(l=k+1;l<n;l++) {
      mat_k[l] = invM_k[msa[l]];
    }
  }

  #pragma loop noalias
  for(k=0;k<n2;k++) {
    #pragma loop norecurrence
    for(l=0;l<k;l++) {
      mat[n2*k + l] = -mat[n2*l + k]; /* transpose */
    }
    mat[n2*k + k] = 0.0; /* diagonal elements */
  }

  M_ZSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, rwork, &info);
  sgn = ( (n*(n-1)/2)%2==0 ) ? 1.0 : -1.0;

  return sgn * pfaff * PfM[qpidx];
}

/* Calculate n-body Green function */
/* <phi| c1 a1 c2 a2 ... cn an |x> */
/* c_k = c_rsi[k], a_k = a_rsj[k] */
/* rsi is the array of indices of creation operators */
/* rsj is the array of indices of annihilation operators */
/* buffer size    = NQPFull + n*Nsize */
/* bufferInt size = NProj */

//double complex GreenFuncN(const int n, int *rsi, int *rsj, const double complex ip,
//                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
//                  double complex *buffer, int *bufferInt){
//  int ri,rj,rk,rl,si,sj,sk,mj;
//  int k,l,m,rsk;
//  double complex z,x;
//  int qpidx;
//
//  int *projCntNew = bufferInt; /* [NProj] */
//
//  int msj[n];
//  double complex *pfMNew = buffer; /* [NQPFull] */
//  double complex *bufV = buffer+NQPFull; /* [n*Nsize] */
//
//  for(k=0;k<n;k++) {
//    si = rsi[k]/Nsite;
//    sj = rsj[k]/Nsite;
//    if(si!=sj) return 0;
//  }
//
//  if(n<=0) return 0;
//  else if(n==1) {
//    ri = rsi[0]%Nsite;
//    rj = rsj[0]%Nsite;
//    si = rsi[0]/Nsite;
//    return GreenFunc1(ri,rj,si,ip,eleIdx,eleCfg,eleNum,eleProjCnt,projCntNew,buffer);
//  } else if(n==2) {
//    ri = rsi[0]%Nsite;
//    rj = rsj[0]%Nsite;
//    si = rsi[0]/Nsite;
//    rk = rsi[1]%Nsite;
//    rl = rsj[1]%Nsite;
//    sk = rsi[1]/Nsite;
//    return GreenFunc2(ri,rj,rk,rl,si,sk,ip,eleIdx,eleCfg,eleNum,eleProjCnt,projCntNew,buffer);
//  }
//
//  /* reduction */
//  for(k=n-1;k>=0;k--) {
//    /* ** check for an annihilation operator at rsj[k] ** */
//    rsk = rsj[k];
//    for(l=k+1;l<n;l++) {
//      /* rsj[k] == rsi[l] */
//      if(rsk==rsi[l]) {
//        rsj[k] = rsj[l];
//        for(m=l;m<n-1;m++) { /* shift */
//          rsi[m] = rsi[m+1];
//          rsj[m] = rsj[m+1];
//        }
//        return GreenFuncN(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,eleProjCnt,buffer,bufferInt);
//      }
//      /* rsj[k] == rsj[l] */
//      if(rsk==rsj[l]) return 0;
//    }
//    /* check electron number */
//    if(eleNum[rsk]==0) return 0;
//
//    /* ** check for a creation operator at rsi[k] ** */
//    rsk = rsi[k];
//    /* rsi[k] == rsj[k] */
//    if(rsk==rsj[k]) {
//      for(m=k;m<n-1;m++) { /* shift */
//        rsi[m] = rsi[m+1];
//        rsj[m] = rsj[m+1];
//      }
//      return GreenFuncN(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,eleProjCnt,buffer,bufferInt);
//    }
//    for(l=k+1;l<n;l++) {
//      /* rsi[k] == rsi[l] */
//      if(rsk==rsi[l]) return 0;
//      /* rsi[k] == rsj[l] (k<l) */
//      if(rsk==rsj[l]) {
//        rsi[k] = rsi[l];
//        for(m=l;m<n-1;m++) { /* shift */
//          rsi[m] = rsi[m+1];
//          rsj[m] = rsj[m+1];
//        }
//        return (-1.0)*GreenFuncN(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,eleProjCnt,buffer,bufferInt);
//      }
//    }
//    /* check electron number */
//    if(eleNum[rsk]==1) return 0;
//  }
//
//  /* hopping */
//  #pragma loop noalias
//  for(k=0;k<n;k++) {
//    ri = rsi[k]%Nsite;
//    rj = rsj[k]%Nsite;
//    sj = rsj[k]/Nsite;
//    mj = eleCfg[rsj[k]];
//    msj[k] = mj + sj*Ne;
//
//    eleIdx[msj[k]] = ri;
//    eleNum[rsj[k]] = 0;
//    eleNum[rsi[k]] = 1;
//  }
//
//  MakeProjCnt(projCntNew,eleNum);
//  x = ProjRatio(projCntNew,eleProjCnt);
//
//  /* calculateNewPfM */
//  for(qpidx=0;qpidx<NQPFull;qpidx++) {
//    pfMNew[qpidx] = calculateNewPfMN_child(qpidx,n,msj,rsj,eleIdx,bufV);
//  }
//  z = CalculateIP(pfMNew, 0, NQPFull, MPI_COMM_SELF);
//
//  /* revert hoppint */
//  #pragma loop noalias
//  for(k=0;k<n;k++) {
//    rj = rsj[k]%Nsite;
//    eleIdx[msj[k]] = rj;
//    eleNum[rsj[k]] = 1;
//    eleNum[rsi[k]] = 0;
//  }
//
//  return x*z/ip;
//}
//
///* msa[k]-th electron hops from rsa[k] to eleIdx[msa[k]] */
///* buffer size = n*Nsize */
//double complex calculateNewPfMN_child(const int qpidx, const int n, const int *msa, const int *rsa,
//                              const int *eleIdx, double complex *buffer) {
//  const int nsize = Nsize;
//  const int n2 = 2*n;
//  const double complex *sltE;
//  const double complex *sltE_k;
//  const double complex *invM;
//  const double complex *invM_i, *invM_k, *invM_l;
//
//  double complex *vec; /* vec[n][nsize] */
//  double complex *vec_k, *vec_l;
//  double complex mat[n2*n2]; /* mat[n2][n2] */
//  double complex *mat_k;
//  double sgn;
//
//  int rsi,rsk,msi,msj,k,l;
//  double complex val,tmp;
//
//  /* for DSKPFA */
//  char uplo='U', mthd='P';
//  int nn,lda,info=0;
//  double complex pfaff;
//  int iwork[n2];
//  double complex work[n2*n2]; /* [n2][n2] */
//  int lwork = n2*n2;
//  nn=lda=n2;
//
//  sltE = SlaterElm + qpidx*Nsite2*Nsite2;
//  invM = InvM + qpidx*Nsize*Nsize;
//
//  vec = buffer; /* n*nsize */
//
//  #pragma loop noalias
//  for(k=0;k<n;k++) {
//    rsk = eleIdx[msa[k]] + (msa[k]/Ne)*Nsite;
//    sltE_k = sltE + rsk*Nsite2;
//    vec_k = vec + k*nsize;
//    #pragma loop norecurrence
//    for(msi=0;msi<nsize;msi++) {
//      rsi = eleIdx[msi] + (msi/Ne)*Nsite;
//      vec_k[msi] = sltE_k[rsi];
//    }
//  }
//
//  /* X_kl */
//  for(k=0;k<n;k++) {
//    mat_k = mat + n2*k;
//    vec_k = vec + k*nsize;
//    for(l=k+1;l<n;l++) {
//      vec_l = vec + l*nsize;
//      val = 0.0;
//      for(msi=0;msi<nsize;msi++) {
//        invM_i = invM + msi*nsize;
//        tmp = 0.0;
//        for(msj=0;msj<nsize;msj++) {
//          tmp += invM_i[msj] * vec_l[msj];
//        }
//        val += tmp * vec_k[msi];
//      }
//      mat_k[l] = val + vec_k[msa[l]];
//    }
//  }
//
//  /* Y_kl */
//  for(k=0;k<n;k++) {
//    mat_k = mat + n2*k + n;
//    vec_k = vec + k*nsize;
//    for(l=0;l<n;l++) {
//      invM_l = invM + msa[l]*nsize;
//      val = 0.0;
//      for(msi=0;msi<nsize;msi++) {
//        val += vec_k[msi] * invM_l[msi];
//      }
//      mat_k[l] = val;
//    }
//  }
//
//  /* Z_kl */
//  for(k=0;k<n;k++) {
//    mat_k = mat + n2*(k+n) + n;
//    invM_k = invM + msa[k]*nsize;
//    for(l=k+1;l<n;l++) {
//      mat_k[l] = invM_k[msa[l]];
//    }
//  }
//
//  #pragma loop noalias
//  for(k=0;k<n2;k++) {
//    #pragma loop norecurrence
//    for(l=0;l<k;l++) {
//      mat[n2*k + l] = -mat[n2*l + k]; /* transpose */
//    }
//    mat[n2*k + k] = 0.0; /* diagonal elements */
//  }
//
//  /* calculate Pf M */
//  // Coverage-0 and commented out,
//  // but should be useful sometimes.
//  //M_DSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, &info);
//  info = 1; // Skipping inverse.
//  M_ZSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff, iwork, work, &lwork/*, rwork*/, &info);
//
//  sgn = ( (n*(n-1)/2)%2==0 ) ? 1.0 : -1.0;
//
//  return sgn * pfaff * PfM[qpidx];
//}
