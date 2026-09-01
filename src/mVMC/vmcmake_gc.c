#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../sfmt/SFMT.h"
#include "include/gc_config.h"
#include "include/gc_size.h"
#include "include/global.h"
#include "include/matrix_gc.h"
#include "include/pfupdate_gc.h"
#include "include/projection.h"
#include "include/qp.h"
#include "include/splitloop.h"
#include "include/vmcmake_gc.h"
#include "include/workspace.h"

#ifndef _SRC_VMCMAKE_GC
#define _SRC_VMCMAKE_GC

void UpdateProjCnt_fsz(int ri, int rj, int s, int t, int *projCntNew,
                       const int *projCntOld, const int *eleNum);

static double GCMoveClassProbability(const enum GCMoveClass moveClass,
                                     const int ncur,
                                     const int nsite2) {
  const int hopPossible = ncur > 0 && ncur < nsite2;
  const int addPossible = nsite2 - ncur >= 2;
  const int removePossible = ncur >= 2;
  const double normalization = 0.5 * (double)hopPossible +
                               0.25 * (double)addPossible +
                               0.25 * (double)removePossible;
  if (normalization == 0.0) return 0.0;
  if (moveClass == GC_MOVE_HOP) {
    return hopPossible ? 0.5 / normalization : 0.0;
  }
  if (moveClass == GC_MOVE_ADD) {
    return addPossible ? 0.25 / normalization : 0.0;
  }
  return removePossible ? 0.25 / normalization : 0.0;
}

enum GCMoveClass GCSelectMoveClass(const double classDraw, const int ncur,
                                   const int nsite2) {
  const double pHop =
      GCMoveClassProbability(GC_MOVE_HOP, ncur, nsite2);
  const double pAdd =
      GCMoveClassProbability(GC_MOVE_ADD, ncur, nsite2);
  if (classDraw < pHop) return GC_MOVE_HOP;
  if (classDraw < pHop + pAdd) return GC_MOVE_ADD;
  return GC_MOVE_REMOVE;
}

static void GCMutateHop(const int ma, const int rsa, int *eleIdx,
                        int *eleCfg, int *eleNum, int *oldRs) {
  *oldRs = eleIdx[ma];
  eleIdx[ma] = rsa;
  eleCfg[*oldRs] = -1;
  eleCfg[rsa] = ma;
  eleNum[*oldRs] = 0;
  eleNum[rsa] = 1;
}

static void GCRestoreHop(const int ma, const int rsa, const int oldRs,
                         int *eleIdx, int *eleCfg, int *eleNum) {
  eleIdx[ma] = oldRs;
  eleCfg[rsa] = -1;
  eleCfg[oldRs] = ma;
  eleNum[rsa] = 0;
  eleNum[oldRs] = 1;
}

int GCAttemptMove(const enum GCMoveClass moveClass, const int arg0,
                  const int arg1, const double acceptDraw, int *eleIdx,
                  int *eleCfg, int *eleNum, int *eleProjCnt,
                  double complex *logIpOld, double complex *pfMNew,
                  int *projCntNew, const int qpStart, const int qpEnd,
                  MPI_Comm comm) {
  const int ncurOld = Ncur;
  int ncurProposal = ncurOld;
  int oldRs = -1;
  int oldTail0 = 0;
  int oldTail1 = 0;
  double proposalRatio = 1.0;
  double complex logIpNew;
  double x;
  double weight;
  int accepted;

  if (eleIdx == NULL || eleCfg == NULL || eleNum == NULL ||
      (NProj > 0 && (eleProjCnt == NULL || projCntNew == NULL)) ||
      logIpOld == NULL || pfMNew == NULL || acceptDraw < 0.0 ||
      acceptDraw >= 1.0) {
    return 0;
  }

  if (moveClass == GC_MOVE_HOP) {
    int oldSite;
    int newSite;
    int oldSpin;
    int newSpin;
    if (arg0 < 0 || arg0 >= ncurOld || arg1 < 0 || arg1 >= Nsite2 ||
        eleNum[arg1] != 0) {
      return 0;
    }
    GCMutateHop(arg0, arg1, eleIdx, eleCfg, eleNum, &oldRs);
    oldSite = oldRs % Nsite;
    newSite = arg1 % Nsite;
    oldSpin = oldRs / Nsite;
    newSpin = arg1 / Nsite;
    if (oldSpin == newSpin) {
      UpdateProjCnt(oldSite, newSite, oldSpin, projCntNew, eleProjCnt,
                    eleNum);
    } else {
      UpdateProjCnt_fsz(oldSite, newSite, oldSpin, newSpin, projCntNew,
                        eleProjCnt, eleNum);
    }
    CalculateNewPfMHopGC(arg0, pfMNew, eleIdx, ncurOld, qpStart, qpEnd);
  } else if (moveClass == GC_MOVE_ADD) {
    const double pAddX =
        GCMoveClassProbability(GC_MOVE_ADD, ncurOld, Nsite2);
    const double pRemoveY =
        GCMoveClassProbability(GC_MOVE_REMOVE, ncurOld + 2, Nsite2);
    if (arg0 < 0 || arg0 >= Nsite2 || arg1 < 0 || arg1 >= Nsite2 ||
        arg0 == arg1 || eleNum[arg0] != 0 || eleNum[arg1] != 0 ||
        ncurOld + 2 > NsizeMax) {
      return 0;
    }
    oldTail0 = eleIdx[ncurOld];
    oldTail1 = eleIdx[ncurOld + 1];
    GCAddPair(arg0, arg1, eleIdx, eleCfg, eleNum, &ncurProposal);
    MakeProjCnt(projCntNew, eleNum);
    CalculateNewPfMAddGC(arg0, arg1, pfMNew, eleIdx, ncurOld, qpStart,
                         qpEnd);
    proposalRatio =
        GCProposalRatioAdd(ncurOld, Nsite2, pAddX, pRemoveY);
  } else if (moveClass == GC_MOVE_REMOVE) {
    const double pRemoveX =
        GCMoveClassProbability(GC_MOVE_REMOVE, ncurOld, Nsite2);
    const double pAddY =
        GCMoveClassProbability(GC_MOVE_ADD, ncurOld - 2, Nsite2);
    if (arg0 < 0 || arg0 >= ncurOld || arg1 <= arg0 ||
        arg1 >= ncurOld || ncurOld < 2) {
      return 0;
    }
    (void)GCRemovePair(arg0, arg1, eleIdx, eleCfg, eleNum,
                       &ncurProposal);
    MakeProjCnt(projCntNew, eleNum);
    CalculateNewPfMRemoveGC(arg0, arg1, pfMNew, eleIdx, ncurOld,
                            qpStart, qpEnd);
    proposalRatio =
        GCProposalRatioRemove(ncurOld, Nsite2, pRemoveX, pAddY);
  } else {
    return 0;
  }

  logIpNew = CalculateLogIP_fcmp(pfMNew, qpStart, qpEnd, comm);
  x = LogProjRatio(projCntNew, eleProjCnt);
  weight = exp(2.0 * (x + creal(logIpNew - *logIpOld))) * proposalRatio;
  accepted = isfinite(weight) && weight > acceptDraw;
  if (accepted) {
    if (moveClass == GC_MOVE_HOP) {
      UpdateMAllHopGC(arg0, eleIdx, ncurOld, qpStart, qpEnd);
    } else if (moveClass == GC_MOVE_ADD) {
      UpdateMAllAddGC(arg0, arg1, eleIdx, ncurOld, qpStart, qpEnd);
    } else {
      UpdateMAllRemoveGC(arg0, arg1, eleIdx, ncurOld, qpStart, qpEnd);
    }
    if (NProj > 0) {
      memcpy(eleProjCnt, projCntNew, (size_t)NProj * sizeof(*eleProjCnt));
    }
    Ncur = ncurProposal;
    *logIpOld = logIpNew;
    return 1;
  }

  if (moveClass == GC_MOVE_HOP) {
    GCRestoreHop(arg0, arg1, oldRs, eleIdx, eleCfg, eleNum);
  } else if (moveClass == GC_MOVE_ADD) {
    eleCfg[arg0] = -1;
    eleCfg[arg1] = -1;
    eleNum[arg0] = 0;
    eleNum[arg1] = 0;
    eleIdx[ncurOld] = oldTail0;
    eleIdx[ncurOld + 1] = oldTail1;
  } else {
    GCRestoreRemovedPair(arg0, arg1, eleIdx, eleCfg, eleNum,
                         &ncurProposal);
  }
  return 0;
}

static int GCSelectTwoDistinct(const int count, int *first, int *second) {
  int a;
  int b;
  if (count < 2) return 0;
  a = (int)(gen_rand32() % (uint32_t)count);
  b = (int)(gen_rand32() % (uint32_t)(count - 1));
  if (b >= a) b++;
  if (a > b) {
    const int temporary = a;
    a = b;
    b = temporary;
  }
  *first = a;
  *second = b;
  return 1;
}

static int GCDebugRebuildInterval(void) {
  static int initialized = 0;
  static int interval = 0;
  if (!initialized) {
    const char *value = getenv("MVMC_GC_DEBUG_REBUILD_INTERVAL");
    if (value != NULL && *value != '\0') {
      char *end = NULL;
      long parsed;
      errno = 0;
      parsed = strtol(value, &end, 10);
      if (errno == 0 && end != value && *end == '\0' && parsed > 0 &&
          parsed <= 1000000000L) {
        interval = (int)parsed;
      }
    }
    initialized = 1;
  }
  return interval;
}

static FILE *GCOpenStateDump(void) {
  const char *base = getenv("MVMC_GC_STATE_DUMP");
  char path[D_FileNameMax + 64];
  int worldRank = 0;
  int written;
  if (base == NULL || *base == '\0') return NULL;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  written = snprintf(path, sizeof(path), "%s.rank%d", base, worldRank);
  if (written < 0 || (size_t)written >= sizeof(path)) {
    fprintf(stderr, "Error: MVMC_GC_STATE_DUMP path is too long.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  return fopen(path, "a");
}

static void GCWriteStateDump(FILE *fp, const char *label, const int sample,
                             const int *eleIdx, const int *eleCfg,
                             const int *eleNum, const int *eleProjCnt) {
  int i;
  if (fp == NULL) return;
  fprintf(fp, "%s %d %d", label, sample, Ncur);
  for (i = 0; i < NsizeMax; i++) fprintf(fp, " %d", eleIdx[i]);
  fprintf(fp, " |");
  for (i = 0; i < Nsite2; i++) fprintf(fp, " %d", eleCfg[i]);
  fprintf(fp, " |");
  for (i = 0; i < Nsite2; i++) fprintf(fp, " %d", eleNum[i]);
  fprintf(fp, " |");
  for (i = 0; i < NProj; i++) fprintf(fp, " %d", eleProjCnt[i]);
  fprintf(fp, "\n");
}

static void GCDebugCheckFastState(const int *eleIdx,
                                  double complex *logIpOld,
                                  const int qpStart, const int qpEnd,
                                  MPI_Comm comm) {
  const size_t invStride =
      GCCheckedMulSize((size_t)NsizeMax, (size_t)NsizeMax);
  const size_t qpCount = (size_t)(qpEnd - qpStart);
  const size_t invCount = GCCheckedMulSize(qpCount, invStride);
  const size_t totalCount = GCCheckedAddSize(invCount, qpCount);
  const int workspaceCount =
      GCCheckedSizeToInt(totalCount, "GC debug rebuild workspace");
  size_t index;
  RequestWorkSpaceComplex(workspaceCount);
  {
    double complex *fastInv = GetWorkSpaceComplex((int)invCount);
    double complex *fastPf = GetWorkSpaceComplex((int)qpCount);
    memcpy(fastInv, InvM, invCount * sizeof(*fastInv));
    memcpy(fastPf, PfM, qpCount * sizeof(*fastPf));
    if (CalculateMAllGC_fcmp(Ncur, eleIdx, qpStart, qpEnd) != GC_MALL_OK) {
      fprintf(stderr, "Error: GC debug full rebuild failed.\n");
      MPI_Abort(comm, EXIT_FAILURE);
    }
    for (index = 0; index < qpCount; index++) {
      const double error = cabs(fastPf[index] - PfM[index]) /
                           fmax(1.0, cabs(PfM[index]));
      if (!(isfinite(error) && error <= 1.0e-8)) {
        fprintf(stderr,
                "Error: GC debug Pfaffian mismatch qp=%zu error=%.17e.\n",
                index + (size_t)qpStart, error);
        MPI_Abort(comm, EXIT_FAILURE);
      }
    }
    for (index = 0; index < qpCount; index++) {
      int row;
      for (row = 0; row < Ncur; row++) {
        int column;
        for (column = 0; column < Ncur; column++) {
          const size_t offset = index * invStride +
                                (size_t)row * (size_t)NsizeMax +
                                (size_t)column;
          const double error = cabs(fastInv[offset] - InvM[offset]) /
                               fmax(1.0, cabs(InvM[offset]));
          if (!(isfinite(error) && error <= 1.0e-8)) {
            fprintf(stderr,
                    "Error: GC debug inverse mismatch qp=%zu row=%d "
                    "column=%d error=%.17e.\n",
                    index + (size_t)qpStart, row, column, error);
            MPI_Abort(comm, EXIT_FAILURE);
          }
        }
      }
    }
  }
  ReleaseWorkSpaceComplex();
  *logIpOld = CalculateLogIP_fcmp(PfM, qpStart, qpEnd, comm);
}

int GCMakeOneStep(int *eleIdx, int *eleCfg, int *eleNum,
                  int *eleProjCnt, double complex *logIpOld,
                  double complex *pfMNew, int *projCntNew,
                  const int qpStart, const int qpEnd, MPI_Comm comm) {
  static long long debugAccepted = 0;
  const enum GCMoveClass moveClass =
      GCSelectMoveClass(genrand_real2(), Ncur, Nsite2);
  int arg0;
  int arg1;
  int accepted;
  if (moveClass == GC_MOVE_HOP) {
    const int emptyCount = Nsite2 - Ncur;
    arg0 = (int)(gen_rand32() % (uint32_t)Ncur);
    arg1 = GCFindKthEmpty(eleNum, Nsite2,
                          (int)(gen_rand32() % (uint32_t)emptyCount));
  } else if (moveClass == GC_MOVE_ADD) {
    int empty0;
    int empty1;
    const int emptyCount = Nsite2 - Ncur;
    if (!GCSelectTwoDistinct(emptyCount, &empty0, &empty1)) return 0;
    arg0 = GCFindKthEmpty(eleNum, Nsite2, empty0);
    arg1 = GCFindKthEmpty(eleNum, Nsite2, empty1);
  } else {
    if (!GCSelectTwoDistinct(Ncur, &arg0, &arg1)) return 0;
  }
  accepted = GCAttemptMove(moveClass, arg0, arg1, genrand_real2(), eleIdx,
                           eleCfg, eleNum, eleProjCnt, logIpOld, pfMNew,
                           projCntNew, qpStart, qpEnd, comm);
  if (accepted && GCDebugRebuildInterval() > 0) {
    debugAccepted++;
    if (debugAccepted % GCDebugRebuildInterval() == 0) {
      GCDebugCheckFastState(eleIdx, logIpOld, qpStart, qpEnd, comm);
    }
  }
  return accepted;
}

static int GCCollectiveRebuildStatus(const int localStatus, MPI_Comm comm) {
  int globalStatus = GC_MALL_OK;
  MPI_Allreduce(&localStatus, &globalStatus, 1, MPI_INT, MPI_MAX, comm);
  return globalStatus;
}

int makeInitialSampleGC(int *eleIdx, int *eleCfg, int *eleNum,
                        int *eleProjCnt, const int qpStart,
                        const int qpEnd, MPI_Comm comm) {
  int attempt;
  int rank;
  MPI_Comm_rank(comm, &rank);
  for (attempt = 0; attempt < 100; attempt++) {
    int position;
    int status;
    int globalStatus;
    for (position = 0; position < NsizeMax; position++) eleIdx[position] = -1;
    for (position = 0; position < Nsite2; position++) {
      eleCfg[position] = -1;
      eleNum[position] = 0;
    }
    for (position = 0; position < Ncur; position++) {
      const int emptyCount = Nsite2 - position;
      const int kth = (int)(gen_rand32() % (uint32_t)emptyCount);
      const int rs = GCFindKthEmpty(eleNum, Nsite2, kth);
      eleIdx[position] = rs;
      eleCfg[rs] = position;
      eleNum[rs] = 1;
    }
    MakeProjCnt(eleProjCnt, eleNum);
    status = CalculateMAllGC_fcmp(Ncur, eleIdx, qpStart, qpEnd);
    globalStatus = GCCollectiveRebuildStatus(status, comm);
    if (globalStatus == GC_MALL_OK) return 0;
  }
  if (rank == 0) {
    fprintf(stderr, "Error: makeInitialSampleGC exceeded 100 attempts.\n");
  }
  MPI_Abort(comm, EXIT_FAILURE);
  return 1;
}

void copyFromBurnSampleGC(int *eleIdx, int *eleCfg, int *eleNum,
                          int *eleProjCnt) {
  memcpy(eleIdx, BurnEleIdx, (size_t)NsizeMax * sizeof(*eleIdx));
  memcpy(eleCfg, BurnEleCfg, (size_t)Nsite2 * sizeof(*eleCfg));
  memcpy(eleNum, BurnEleNum, (size_t)Nsite2 * sizeof(*eleNum));
  if (NProj > 0) {
    memcpy(eleProjCnt, BurnEleProjCnt, (size_t)NProj * sizeof(*eleProjCnt));
  }
  Ncur = BurnNcur;
}

void copyToBurnSampleGC(const int *eleIdx, const int *eleCfg,
                        const int *eleNum, const int *eleProjCnt) {
  memcpy(BurnEleIdx, eleIdx, (size_t)NsizeMax * sizeof(*eleIdx));
  memcpy(BurnEleCfg, eleCfg, (size_t)Nsite2 * sizeof(*eleCfg));
  memcpy(BurnEleNum, eleNum, (size_t)Nsite2 * sizeof(*eleNum));
  if (NProj > 0) {
    memcpy(BurnEleProjCnt, eleProjCnt,
           (size_t)NProj * sizeof(*eleProjCnt));
  }
  BurnNcur = Ncur;
}

void saveEleConfigGC(const int sample, const double complex logIp,
                     const int *eleIdx, const int *eleCfg,
                     const int *eleNum, const int *eleProjCnt,
                     const int ncur) {
  memcpy(EleIdx + (size_t)sample * (size_t)NsizeMax, eleIdx,
         (size_t)NsizeMax * sizeof(*EleIdx));
  memcpy(EleCfg + (size_t)sample * (size_t)Nsite2, eleCfg,
         (size_t)Nsite2 * sizeof(*EleCfg));
  memcpy(EleNum + (size_t)sample * (size_t)Nsite2, eleNum,
         (size_t)Nsite2 * sizeof(*EleNum));
  if (NProj > 0) {
    memcpy(EleProjCnt + (size_t)sample * (size_t)NProj, eleProjCnt,
           (size_t)NProj * sizeof(*EleProjCnt));
  }
  EleNumSample[sample] = ncur;
  logSqPfFullSlater[sample] =
      2.0 * (LogProjVal(eleProjCnt) + creal(logIp));
}

void VMCMakeSampleGC(MPI_Comm comm) {
  const size_t qpCount = GCNonnegativeSize(NQPFull, "GC sampler QP count");
  double complex *pfMNew = (double complex *)GCCheckedMalloc(
      GCCheckedMulSize(sizeof(*pfMNew), qpCount), "GC sampler Pf scratch");
  double complex logIpOld;
  int *projCntNew;
  int qpStart;
  int qpEnd;
  int rank;
  int size;
  int nAccept = 0;
  int nOutStep;
  int nInStep;
  int outStep;
  const char *stateDumpPath = getenv("MVMC_GC_STATE_DUMP");
  FILE *stateDump = GCOpenStateDump();
  if (stateDumpPath != NULL && *stateDumpPath != '\0' && stateDump == NULL) {
    fprintf(stderr, "Error: failed to open MVMC_GC_STATE_DUMP.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  SplitLoop(&qpStart, &qpEnd, NQPFull, rank, size);
  RequestWorkSpaceInt(NProj);
  projCntNew = GetWorkSpaceInt(NProj);

  if (BurnFlag == 0) {
    (void)makeInitialSampleGC(TmpEleIdx, TmpEleCfg, TmpEleNum,
                              TmpEleProjCnt, qpStart, qpEnd, comm);
  } else {
    int rebuildStatus;
    copyFromBurnSampleGC(TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt);
    GCWriteStateDump(stateDump, "RESTORE", -1, TmpEleIdx, TmpEleCfg,
                     TmpEleNum, TmpEleProjCnt);
    rebuildStatus = CalculateMAllGC_fcmp(Ncur, TmpEleIdx, qpStart, qpEnd);
    rebuildStatus = GCCollectiveRebuildStatus(rebuildStatus, comm);
    if (rebuildStatus != GC_MALL_OK) {
      (void)makeInitialSampleGC(TmpEleIdx, TmpEleCfg, TmpEleNum,
                                TmpEleProjCnt, qpStart, qpEnd, comm);
      BurnFlag = 0;
    }
  }
  logIpOld = CalculateLogIP_fcmp(PfM, qpStart, qpEnd, comm);
  if (!isfinite(creal(logIpOld)) || !isfinite(cimag(logIpOld))) {
    (void)makeInitialSampleGC(TmpEleIdx, TmpEleCfg, TmpEleNum,
                              TmpEleProjCnt, qpStart, qpEnd, comm);
    logIpOld = CalculateLogIP_fcmp(PfM, qpStart, qpEnd, comm);
    BurnFlag = 0;
  }

  nOutStep = (BurnFlag == 0) ? NVMCWarmUp + NVMCSample : NVMCSample + 1;
  nInStep = NVMCInterval * Nsite;
  for (outStep = 0; outStep < nOutStep; outStep++) {
    int inStep;
    for (inStep = 0; inStep < nInStep; inStep++) {
      nAccept += GCMakeOneStep(TmpEleIdx, TmpEleCfg, TmpEleNum,
                               TmpEleProjCnt, &logIpOld, pfMNew,
                               projCntNew, qpStart, qpEnd, comm);
      if (nAccept > Nsite) {
        if (CalculateMAllGC_fcmp(Ncur, TmpEleIdx, qpStart, qpEnd) !=
            GC_MALL_OK) {
          fprintf(stderr, "Error: periodic GC sampler rebuild failed.\n");
          MPI_Abort(comm, EXIT_FAILURE);
        }
        logIpOld = CalculateLogIP_fcmp(PfM, qpStart, qpEnd, comm);
        nAccept = 0;
      }
    }
    if (outStep >= nOutStep - NVMCSample) {
      const int sample = outStep - (nOutStep - NVMCSample);
      saveEleConfigGC(sample, logIpOld, TmpEleIdx, TmpEleCfg, TmpEleNum,
                      TmpEleProjCnt, Ncur);
      GCWriteStateDump(stateDump, "SAMPLE", sample, TmpEleIdx, TmpEleCfg,
                       TmpEleNum, TmpEleProjCnt);
    }
  }
  copyToBurnSampleGC(TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt);
  GCWriteStateDump(stateDump, "STORE", -1, TmpEleIdx, TmpEleCfg,
                   TmpEleNum, TmpEleProjCnt);
  BurnFlag = 1;
  if (stateDump != NULL) fclose(stateDump);
  ReleaseWorkSpaceInt();
  free(pfMNew);
}

#endif
