#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../sfmt/SFMT.h"
#include "include/gc_antiparallel.h"
#include "include/gc_config.h"
#include "include/gc_size.h"
#include "include/gc_spin_support.h"
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

static int GCAntiMode(void) {
  return GCAntiEnabled(FlagGrandCanonical, iFlgOrbitalGeneral);
}

static uint32_t GCAntiNext32(void) { return gen_rand32(); }

static int GCActiveSpinProjection(void) {
  return GCAntiMode() && NSPGaussLeg > 1 && NSPStot == 0;
}

static int GCForbiddenSpinSector(const int n) {
  return GCSpinForbiddenSector(FlagGrandCanonical, iFlgOrbitalGeneral,
                               NSPGaussLeg, NSPStot, GCSpinTriplet, n);
}

static void GCRequireSpinSupport(MPI_Comm comm) {
  int low = GCSpinTriplet, high = GCSpinTriplet;
  int invalid, singlet = 0;
  if (!GCActiveSpinProjection() && !GCFiniteSpinActive()) return;
  MPI_Allreduce(&GCSpinTriplet, &low, 1, MPI_INT, MPI_MIN, comm);
  MPI_Allreduce(&GCSpinTriplet, &high, 1, MPI_INT, MPI_MAX, comm);
  invalid = low < 0 || high > 1 || low != high;
  if (GCFiniteSpinActive()) {
    int singletLow = GCSpinSinglet, singletHigh = GCSpinSinglet;
    MPI_Allreduce(&GCSpinSinglet, &singletLow, 1, MPI_INT, MPI_MIN, comm);
    MPI_Allreduce(&GCSpinSinglet, &singletHigh, 1, MPI_INT, MPI_MAX, comm);
    invalid = invalid || singletLow < 0 || singletHigh > 1 ||
              singletLow != singletHigh;
    singlet = !invalid && singletLow == 1;
  }
  if (invalid || singlet) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
      if (invalid)
        fprintf(stderr,
                "Error: invalid or rank-inconsistent GC spin support.\n");
      else
        fprintf(stderr,
                "Error: the pair function is exactly symmetric (pure "
                "singlet) and has no total spin S=%d component.\n",
                NSPStot);
    }
    MPI_Barrier(comm);
    MPI_Abort(comm, EXIT_FAILURE);
  }
}

static void GCRecordCompositeZero(MPI_Comm comm) {
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0)
    GCComponentCountersOf(GCComponentThreadCache())->compositeZero++;
}

static void GCRequireRebuiltAmplitude(const double complex logIp,
                                      const char *label, MPI_Comm comm) {
  int rank = 0;
  if (isfinite(creal(logIp)) && isfinite(cimag(logIp))) return;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    if (creal(logIp) == -INFINITY && isfinite(cimag(logIp))) {
      fprintf(stderr,
              "Error: GC sampler state has zero projected amplitude "
              "after %s rebuild.\n",
              label);
    } else {
      fprintf(stderr,
              "Error: GC sampler state has a nonfinite projected "
              "amplitude after %s rebuild.\n",
              label);
    }
  }
  MPI_Abort(comm, EXIT_FAILURE);
}

/* Uniform index in [0, bound); -1 is checked before any int conversion of
 * the UINT32_MAX failure value. */
static int GCAntiDrawIndex(const int bound) {
  uint32_t draw;
  if (bound <= 0) return -1;
  draw = GCAntiDrawBelow((uint32_t)bound, GCAntiNext32);
  if (draw == UINT32_MAX || draw >= (uint32_t)bound) return -1;
  return (int)draw; /* draw < positive int bound: representable */
}

static void GCAntiAbortState(const char *label, MPI_Comm comm) {
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    fprintf(stderr,
            "Error: GC anti-parallel %s state is not a consistent Sz=0 "
            "configuration (Ncur=%d).\n",
            label, Ncur);
  }
  MPI_Abort(comm, EXIT_FAILURE);
}

/* All ranks of a chain hold the same configuration; the verdict is shared
 * so that every rank reaches the same abort before any collective work. */
static void GCAntiRequireConfig(const int *eleIdx, const int *eleCfg,
                                const int *eleNum, const char *label,
                                MPI_Comm comm) {
  const int invalid =
      GCAntiValidateConfig(eleIdx, eleCfg, eleNum, Nsite, Ncur) ? 0 : 1;
  int globalInvalid = invalid;
  MPI_Allreduce(&invalid, &globalInvalid, 1, MPI_INT, MPI_MAX, comm);
  if (globalInvalid != 0) GCAntiAbortState(label, comm);
}

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
  int accepted;
  const int anti = GCAntiMode();

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
    /* Anti-parallel pairs vanish outside Sz=0: a hop keeps its spin. */
    if (anti && eleIdx[arg0] / Nsite != arg1 / Nsite) return 0;
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
    if (anti && arg0 / Nsite == arg1 / Nsite) return 0;
    oldTail0 = eleIdx[ncurOld];
    oldTail1 = eleIdx[ncurOld + 1];
    GCAddPair(arg0, arg1, eleIdx, eleCfg, eleNum, &ncurProposal);
    MakeProjCnt(projCntNew, eleNum);
    CalculateNewPfMAddGC(arg0, arg1, pfMNew, eleIdx, ncurOld, qpStart,
                         qpEnd);
    proposalRatio =
        anti ? GCAntiRatioAdd(ncurOld / 2, Nsite, pAddX, pRemoveY)
             : GCProposalRatioAdd(ncurOld, Nsite2, pAddX, pRemoveY);
  } else if (moveClass == GC_MOVE_REMOVE) {
    const double pRemoveX =
        GCMoveClassProbability(GC_MOVE_REMOVE, ncurOld, Nsite2);
    const double pAddY =
        GCMoveClassProbability(GC_MOVE_ADD, ncurOld - 2, Nsite2);
    if (arg0 < 0 || arg0 >= ncurOld || arg1 <= arg0 ||
        arg1 >= ncurOld || ncurOld < 2) {
      return 0;
    }
    if (anti && eleIdx[arg0] / Nsite == eleIdx[arg1] / Nsite) return 0;
    (void)GCRemovePair(arg0, arg1, eleIdx, eleCfg, eleNum,
                       &ncurProposal);
    MakeProjCnt(projCntNew, eleNum);
    CalculateNewPfMRemoveGC(arg0, arg1, pfMNew, eleIdx, ncurOld,
                            qpStart, qpEnd);
    proposalRatio =
        anti ? GCAntiRatioRemove(ncurOld / 2, Nsite, pRemoveX, pAddY)
             : GCProposalRatioRemove(ncurOld, Nsite2, pRemoveX, pAddY);
  } else {
    return 0;
  }

  /* eleNum already holds the candidate occupation. */
  logIpNew = CalculateLogIPGC_fcmp(pfMNew, ncurProposal, eleNum, qpStart,
                                   qpEnd, comm);
  x = LogProjRatio(projCntNew, eleProjCnt);
  if (GCAntiAcceptLog(*logIpOld, logIpNew, x, proposalRatio, acceptDraw,
                      &accepted) != 0) {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
      fprintf(stderr,
              "Error: GC sampler met a nonfinite or invalid acceptance "
              "input (%s mode, move=%d oldLog=(%.17g,%.17g) "
              "newLog=(%.17g,%.17g) logProjRatio=%.17g ratio=%.17g).\n",
              anti ? "anti-parallel" : "general", (int)moveClass,
              creal(*logIpOld), cimag(*logIpOld), creal(logIpNew),
              cimag(logIpNew), x, proposalRatio);
    }
    MPI_Abort(comm, EXIT_FAILURE);
  }
  if (creal(logIpNew) == -INFINITY) GCRecordCompositeZero(comm);
  if (accepted) {
    const int kind = moveClass == GC_MOVE_HOP
                         ? GC_COMMIT_HOP
                         : moveClass == GC_MOVE_ADD ? GC_COMMIT_ADD
                                                    : GC_COMMIT_REMOVE;
    const int status = GCCommitAccepted(kind, arg0, arg1, pfMNew, eleIdx,
                                        ncurOld, ncurProposal, qpStart,
                                        qpEnd);
    if (status != GC_MALL_OK) {
      int rank = 0;
      MPI_Comm_rank(comm, &rank);
      if (rank == 0)
        fprintf(stderr,
                "Error: GC component commit failed (move=%d status=%d).\n",
                (int)moveClass, status);
      MPI_Abort(comm, EXIT_FAILURE);
    }
    if (NProj > 0)
      memcpy(eleProjCnt, projCntNew, (size_t)NProj * sizeof(*eleProjCnt));
    Ncur = ncurProposal;
    *logIpOld = CalculateLogIPGC_fcmp(PfM, Ncur, eleNum, qpStart, qpEnd,
                                      comm);
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
                             const int *eleNum, const int *eleProjCnt,
                             const double complex *logIp) {
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
#ifdef MVMC_ENABLE_FAULT_INJECTION
  /* Optional lossless sampler value, without another projection/reduction. */
  if (logIp != NULL) {
    const char *dumpAmplitude = getenv("MVMC_GC_STATE_AMPLITUDE");
    if (dumpAmplitude != NULL && strcmp(dumpAmplitude, "1") == 0)
      fprintf(fp, " | LOGIP %a %a", creal(*logIp), cimag(*logIp));
  }
#else
  (void)logIp;
#endif
  fprintf(fp, "\n");
}

static void GCDebugCheckFastState(const int *eleIdx, const int *eleNum,
                                  double complex *logIpOld,
                                  const int qpStart, const int qpEnd,
                                  MPI_Comm comm) {
  const size_t invStride =
      GCCheckedMulSize((size_t)NsizeMax, (size_t)NsizeMax);
  const size_t qpCount = (size_t)(qpEnd - qpStart);
  const size_t invCount = GCCheckedMulSize(qpCount, invStride);
  double complex *fastInv = (double complex *)GCCheckedMalloc(
      GCCheckedMulSize(invCount, sizeof(*fastInv)),
      "GC debug inverse snapshot");
  double complex *fastPf = (double complex *)GCCheckedMalloc(
      GCCheckedMulSize(qpCount, sizeof(*fastPf)), "GC debug Pf snapshot");
  int *fastState = (int *)GCCheckedMalloc(
      GCCheckedMulSize(qpCount, sizeof(*fastState)),
      "GC debug state snapshot");
  size_t index;
  {
    memcpy(fastInv, InvM, invCount * sizeof(*fastInv));
    memcpy(fastPf, PfM, qpCount * sizeof(*fastPf));
    memcpy(fastState, GCQpState, qpCount * sizeof(*fastState));
    if (CalculateMAllGC_fcmp(Ncur, eleIdx, qpStart, qpEnd) != GC_MALL_OK) {
      fprintf(stderr, "Error: GC debug full rebuild failed.\n");
      MPI_Abort(comm, EXIT_FAILURE);
    }
    for (index = 0; index < qpCount; index++) {
      const int fast = fastState[index];
      const int rebuilt = GCQpState[index];
      const double error = cabs(fastPf[index] - PfM[index]) /
                           fmax(1.0, cabs(PfM[index]));
      if (fast == GC_QP_PFONLY && rebuilt == GC_QP_INVERSE) {
        fprintf(stderr,
                "Error: GC debug state mismatch qp=%zu fast=PFONLY "
                "rebuild=INVERSE.\n",
                index + (size_t)qpStart);
        MPI_Abort(comm, EXIT_FAILURE);
      }
      if (fast == GC_QP_INVERSE && rebuilt == GC_QP_PFONLY)
        GCComponentCountersOf(GCComponentThreadCache())
            ->debugNearSingular++;
      if (!(isfinite(error) && error <= 1.0e-8)) {
        fprintf(stderr,
                "Error: GC debug Pfaffian mismatch qp=%zu error=%.17e.\n",
                index + (size_t)qpStart, error);
        MPI_Abort(comm, EXIT_FAILURE);
      }
    }
    for (index = 0; index < qpCount; index++) {
      int row;
      if (fastState[index] != GC_QP_INVERSE ||
          GCQpState[index] != GC_QP_INVERSE)
        continue;
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
  free(fastState);
  free(fastPf);
  free(fastInv);
  *logIpOld = CalculateLogIPGC_fcmp(PfM, Ncur, eleNum, qpStart, qpEnd, comm);
  GCRequireRebuiltAmplitude(*logIpOld, "debug", comm);
}

/* Full endpoint evaluation: never accept two successive N+/-2 moves through
 * a zero intermediate sector. Candidate values use private configuration
 * arrays and component scratch; rejected proposals leave the sampler's
 * arrays, Pfaffians, inverses, and QP states untouched. */
static int GCAttemptBridge(const enum GCAntiBridgeMove move,
                           int *eleIdx, int *eleCfg, int *eleNum,
                           int *eleProjCnt, double complex *logIpOld,
                           double complex *pfMNew, int *projCntNew,
                           const int qpStart, const int qpEnd, MPI_Comm comm) {
  int *candidate = (int *)GCCheckedMalloc(
      GCCheckedMulSize(GCNonnegativeSize(Nsite2, "GC spin bridge orbitals"),
                        2 * sizeof(*candidate)),
      "GC spin bridge configuration");
  int *indices = candidate + Nsite2;
  int count = 0, rs, qp, accepted = 0;
  int localStatus = GC_MALL_OK, status = GC_MALL_OK;
  double ratio = 1.0;
  double complex logIpNew;
  const int drawn = GCAntiDrawBridge(move, eleNum, Nsite, Ncur / 2,
                                     candidate, GCAntiNext32, &ratio);
  if (drawn == 0) {
    free(candidate);
    return 0;
  }
  if (drawn < 0) {
    fprintf(stderr, "Error: invalid GC spin bridge proposal.\n");
    MPI_Abort(comm, EXIT_FAILURE);
  }
  for (rs = 0; rs < Nsite2; ++rs)
    if (candidate[rs]) indices[count++] = rs;
  MakeProjCnt(projCntNew, candidate);
#pragma omp parallel default(shared)
  {
    GCComponentCache *cache = GCComponentThreadCache();
#pragma omp for
    for (qp = qpStart; qp < qpEnd; ++qp) {
      const int result = GCComponentValue(qp, count, indices, cache,
                                          pfMNew + qp - qpStart);
      if (result != GC_MALL_OK) {
#pragma omp critical
        { if (result > localStatus) localStatus = result; }
      }
    }
  }
  MPI_Allreduce(&localStatus, &status, 1, MPI_INT, MPI_MAX, comm);
  if (status != GC_MALL_OK) {
    fprintf(stderr, "Error: GC spin bridge evaluation failed (status=%d).\n",
            status);
    MPI_Abort(comm, EXIT_FAILURE);
  }
  logIpNew = CalculateLogIPGC_fcmp(pfMNew, count, candidate, qpStart, qpEnd,
                                   comm);
  if (GCAntiAcceptLog(*logIpOld, logIpNew,
                      LogProjRatio(projCntNew, eleProjCnt), ratio,
                      genrand_real2(), &accepted) != 0) {
    fprintf(stderr, "Error: invalid GC spin bridge acceptance input.\n");
    MPI_Abort(comm, EXIT_FAILURE);
  }
  if (creal(logIpNew) == -INFINITY) GCRecordCompositeZero(comm);
  if (accepted) {
    for (rs = 0; rs < Nsite2; ++rs) {
      eleNum[rs] = candidate[rs];
      eleCfg[rs] = -1;
    }
    for (rs = 0; rs < NsizeMax; ++rs) eleIdx[rs] = -1;
    for (rs = 0; rs < count; ++rs) {
      eleIdx[rs] = indices[rs];
      eleCfg[indices[rs]] = rs;
    }
    if (NProj > 0)
      memcpy(eleProjCnt, projCntNew, (size_t)NProj * sizeof(*eleProjCnt));
    Ncur = count;
    localStatus = CalculateMAllGC_fcmp(Ncur, eleIdx, qpStart, qpEnd);
    MPI_Allreduce(&localStatus, &status, 1, MPI_INT, MPI_MAX, comm);
    if (status != GC_MALL_OK) {
      fprintf(stderr, "Error: GC spin bridge commit failed (status=%d).\n",
              status);
      MPI_Abort(comm, EXIT_FAILURE);
    }
    *logIpOld = CalculateLogIPGC_fcmp(PfM, Ncur, eleNum, qpStart, qpEnd,
                                      comm);
    GCRequireRebuiltAmplitude(*logIpOld, "spin bridge", comm);
  }
  free(candidate);
  return accepted;
}

int GCMakeOneStep(int *eleIdx, int *eleCfg, int *eleNum,
                  int *eleProjCnt, double complex *logIpOld,
                  double complex *pfMNew, int *projCntNew,
                  const int qpStart, const int qpEnd, MPI_Comm comm) {
  static long long debugAccepted = 0;
  if (GCActiveSpinProjection()) {
    const int bridge = GCAntiSelectBridge(genrand_real2());
    if (bridge >= 0)
      return GCAttemptBridge((enum GCAntiBridgeMove)bridge, eleIdx, eleCfg,
                             eleNum, eleProjCnt, logIpOld, pfMNew,
                             projCntNew, qpStart, qpEnd, comm);
  } else if (GCFiniteSpinActive()) {
    /* Each kernel keeps its own Metropolis-Hastings ratio, so the fixed
     * mixture preserves the common stationary distribution. */
    const int bridge = GCAntiSelectFiniteSpinBridge(genrand_real2());
    if (bridge >= 0)
      return GCAttemptBridge((enum GCAntiBridgeMove)bridge, eleIdx, eleCfg,
                             eleNum, eleProjCnt, logIpOld, pfMNew,
                             projCntNew, qpStart, qpEnd, comm);
  }
  const enum GCMoveClass moveClass =
      GCSelectMoveClass(genrand_real2(), Ncur, Nsite2);
  int arg0;
  int arg1;
  int accepted;
  if (GCAntiMode()) {
    /* Spin-resolved candidates; m pairs occupy m up and m down sites. */
    const int m = Ncur / 2;
    if (moveClass == GC_MOVE_HOP) {
      const int s = GCAntiDrawIndex(2);
      const int oldK = GCAntiDrawIndex(m);
      const int newK = GCAntiDrawIndex(Nsite - m);
      int oldRs;
      int site;
      if (s < 0 || oldK < 0 || newK < 0) return 0;
      oldRs = GCAntiFindOccupied(eleNum, Nsite, s, oldK);
      site = GCFindKthEmpty(eleNum + s * Nsite, Nsite, newK);
      if (oldRs < 0 || oldRs >= Nsite2 || site < 0 || site >= Nsite) {
        return 0;
      }
      arg0 = eleCfg[oldRs];
      if (arg0 < 0 || arg0 >= Ncur) return 0;
      arg1 = s * Nsite + site;
    } else if (moveClass == GC_MOVE_ADD) {
      const int upK = GCAntiDrawIndex(Nsite - m);
      const int downK = GCAntiDrawIndex(Nsite - m);
      int up;
      int down;
      if (upK < 0 || downK < 0) return 0;
      up = GCFindKthEmpty(eleNum, Nsite, upK);
      down = GCFindKthEmpty(eleNum + Nsite, Nsite, downK);
      if (up < 0 || up >= Nsite || down < 0 || down >= Nsite) return 0;
      arg0 = up;
      arg1 = Nsite + down;
    } else if (moveClass == GC_MOVE_REMOVE) {
      const int upK = GCAntiDrawIndex(m);
      const int downK = GCAntiDrawIndex(m);
      int up;
      int down;
      if (upK < 0 || downK < 0) return 0;
      up = GCAntiFindOccupied(eleNum, Nsite, 0, upK);
      down = GCAntiFindOccupied(eleNum, Nsite, 1, downK);
      if (up < 0 || up >= Nsite2 || down < 0 || down >= Nsite2) return 0;
      arg0 = eleCfg[up];
      arg1 = eleCfg[down];
      if (arg0 < 0 || arg0 >= Ncur || arg1 < 0 || arg1 >= Ncur) return 0;
      if (arg0 > arg1) {
        const int temporary = arg0;
        arg0 = arg1;
        arg1 = temporary;
      }
    } else {
      return 0;
    }
  } else if (moveClass == GC_MOVE_HOP) {
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
      GCDebugCheckFastState(eleIdx, eleNum, logIpOld, qpStart, qpEnd,
                            comm);
    }
  }
  return accepted;
}

static int GCCollectiveRebuildStatus(const int localStatus, MPI_Comm comm) {
  int globalStatus = GC_MALL_OK;
  MPI_Allreduce(&localStatus, &globalStatus, 1, MPI_INT, MPI_MAX, comm);
  return globalStatus;
}

static int GCFiniteSpinPermitted(const int electrons) {
  return GCSpinSectorPossible(Nsite, NSPStot, GCSpinTriplet, electrons);
}

/* Append the permitted even N not yet listed, by distance from center (the
 * smaller N first on a tie). Returns the new count. */
static int GCFiniteSpinAppendByDistance(const int center, int *order,
                                        int count) {
  int distance;
  for (distance = 0; distance <= Nsite2; ++distance) {
    int side;
    for (side = 0; side < (distance == 0 ? 1 : 2); ++side) {
      const int n = side == 0 ? center - distance : center + distance;
      int k, listed = 0;
      if (n < 0 || n > Nsite2 || n % 2 != 0 || !GCFiniteSpinPermitted(n))
        continue;
      for (k = 0; k < count; ++k) listed = listed || order[k] == n;
      if (!listed) order[count++] = n;
    }
  }
  return count;
}

/* Finite S>0 start: structured draws with at least 2S singles, at most 100
 * per N and 100*(L+1) in total, trying N in the given order. The first
 * nonzero amplitude is the start; Ncur is set to its N. */
static int GCFiniteSpinInitialize(const int *order, const int count,
                                  int *eleIdx, int *eleCfg, int *eleNum,
                                  int *eleProjCnt, const int qpStart,
                                  const int qpEnd, MPI_Comm comm) {
  const int totalLimit = 100 * (Nsite + 1);
  int rank = 0, k, total = 0;
  MPI_Comm_rank(comm, &rank);
  for (k = 0; k < count && total < totalLimit; ++k) {
    int attempt;
    for (attempt = 0; attempt < 100 && total < totalLimit;
         ++attempt, ++total) {
      int rs, position = 0, status;
      double complex logIp;
      if (GCAntiDrawSpinConfiguration(Nsite, order[k] / 2, NSPStot,
                                      GCAntiNext32, eleNum) != 1) {
        if (rank == 0)
          fprintf(stderr, "Error: invalid GC finite-spin initial draw "
                  "(S=%d, N=%d).\n", NSPStot, order[k]);
        MPI_Abort(comm, EXIT_FAILURE);
      }
      for (rs = 0; rs < NsizeMax; ++rs) eleIdx[rs] = -1;
      for (rs = 0; rs < Nsite2; ++rs) {
        eleCfg[rs] = -1;
        if (eleNum[rs] != 0) {
          eleIdx[position] = rs;
          eleCfg[rs] = position++;
        }
      }
      Ncur = order[k];
      GCAntiRequireConfig(eleIdx, eleCfg, eleNum, "initial", comm);
      MakeProjCnt(eleProjCnt, eleNum);
      status = GCCollectiveRebuildStatus(
          CalculateMAllGC_fcmp(Ncur, eleIdx, qpStart, qpEnd), comm);
      if (status != GC_MALL_OK) {
        if (rank == 0)
          fprintf(stderr, "Error: GC initial rebuild failed (status=%d).\n",
                  status);
        MPI_Abort(comm, EXIT_FAILURE);
      }
      logIp = CalculateLogIPGC_fcmp(PfM, Ncur, eleNum, qpStart, qpEnd, comm);
      if (isfinite(creal(logIp)) && isfinite(cimag(logIp))) return order[k];
      if (!(creal(logIp) == -INFINITY && isfinite(cimag(logIp)))) {
        if (rank == 0)
          fprintf(stderr, "Error: GC initial configuration has a nonfinite "
                  "projected amplitude.\n");
        MPI_Abort(comm, EXIT_FAILURE);
      }
      if (rank == 0)
        GCComponentCountersOf(GCComponentThreadCache())->compositeZero++;
    }
  }
  if (rank == 0)
    fprintf(stderr,
            "Error: GC total-spin S=%d initial configuration has zero "
            "projected amplitude after %d attempt(s) over %d particle "
            "number(s). Check the pair function and NGCInitNelec.\n",
            NSPStot, total, k);
  MPI_Abort(comm, EXIT_FAILURE);
  return -1;
}

/* restartN < 0: first start, where an explicit NGCInitNelec (after the
 * reader's even rounding, i.e. Ncur) is binding and an unspecified one picks
 * the nearest permitted N. restartN >= 0: an SR generation made the saved
 * state zero; try the saved N, then NGCInitNelec, then the rest. */
static void GCFiniteSpinStart(const int restartN, int *eleIdx, int *eleCfg,
                              int *eleNum, int *eleProjCnt, const int qpStart,
                              const int qpEnd, MPI_Comm comm) {
  const size_t slots = GCCheckedAddSize(
      GCNonnegativeSize(Nsite2, "GC finite-spin sectors"), 1);
  int *order = (int *)GCCheckedMalloc(
      GCCheckedMulSize(slots, sizeof(*order)), "GC finite-spin sectors");
  const int requested = Ncur;
  int rank = 0, count = 0, chosen;
  MPI_Comm_rank(comm, &rank);
  if (restartN < 0 && NGCInitNelec >= 0) {
    if (!GCFiniteSpinPermitted(requested)) {
      if (rank == 0)
        fprintf(stderr,
                "Error: NGCInitNelec=%d cannot carry total spin S=%d on "
                "%d sites (%s).\n",
                requested, NSPStot, Nsite,
                2 * NSPStot > requested || requested > Nsite2 - 2 * NSPStot
                    ? "2S <= N <= 2*Nsite-2S is required"
                    : "an exactly antisymmetric pair function requires "
                      "N/2-S even");
      MPI_Barrier(comm);
      MPI_Abort(comm, EXIT_FAILURE);
    }
    order[count++] = requested;
  } else if (restartN < 0) {
    count = GCFiniteSpinAppendByDistance(requested, order, 0);
  } else {
    const int initial = NGCInitNelec - NGCInitNelec % 2;
    if (GCFiniteSpinPermitted(restartN)) order[count++] = restartN;
    if (NGCInitNelec >= 0 && initial != restartN &&
        GCFiniteSpinPermitted(initial))
      order[count++] = initial;
    count = GCFiniteSpinAppendByDistance(restartN, order, count);
  }
  if (count == 0) {
    if (rank == 0)
      fprintf(stderr,
              "Error: no particle number can carry total spin S=%d on %d "
              "sites with this pair function.\n",
              NSPStot, Nsite);
    MPI_Barrier(comm);
    MPI_Abort(comm, EXIT_FAILURE);
  }
  chosen = GCFiniteSpinInitialize(order, count, eleIdx, eleCfg, eleNum,
                                  eleProjCnt, qpStart, qpEnd, comm);
  if (rank == 0 && chosen != requested)
    fprintf(stdout, "  Note: GC S=%d %s particle number is %d (was %d).\n",
            NSPStot, restartN < 0 ? "initial" : "restart", chosen,
            requested);
  free(order);
}

int makeInitialSampleGC(int *eleIdx, int *eleCfg, int *eleNum,
                        int *eleProjCnt, const int qpStart,
                        const int qpEnd, MPI_Comm comm) {
  int attempt;
  int rank;
  double complex logIp;
  MPI_Comm_rank(comm, &rank);
  if (GCFiniteSpinActive()) {
    GCFiniteSpinStart(-1, eleIdx, eleCfg, eleNum, eleProjCnt, qpStart, qpEnd,
                      comm);
    return 0;
  }
  if (GCForbiddenSpinSector(Ncur)) {
    if (rank == 0)
      fprintf(stderr, "Error: NGCInitNelec=%d is forbidden by S=0 spin "
              "projection of an exactly antisymmetric pair function "
              "(N must be a multiple of 4).\n", Ncur);
    MPI_Barrier(comm);
    MPI_Abort(comm, EXIT_FAILURE);
  }
  for (attempt = 0; attempt < 100; attempt++) {
    int position;
    int status;
    int globalStatus;
    for (position = 0; position < NsizeMax; position++) eleIdx[position] = -1;
    for (position = 0; position < Nsite2; position++) {
      eleCfg[position] = -1;
      eleNum[position] = 0;
    }
    if (GCAntiMode()) {
      /* m up and m down electrons, each spin drawn without repetition. */
      int s;
      for (s = 0; s < 2; s++) {
        int k;
        for (k = 0; k < Ncur / 2; k++) {
          const int q = GCAntiDrawIndex(Nsite - k);
          const int site =
              q < 0 ? -1 : GCFindKthEmpty(eleNum + s * Nsite, Nsite, q);
          int rs;
          int pos;
          if (site < 0 || site >= Nsite) {
            if (rank == 0) {
              fprintf(stderr,
                      "Error: GC anti-parallel initialization could not "
                      "place a spin-%d electron (Ncur=%d, Nsite=%d).\n",
                      s, Ncur, Nsite);
            }
            MPI_Abort(comm, EXIT_FAILURE);
          }
          rs = s * Nsite + site;
          pos = s * (Ncur / 2) + k;
          eleIdx[pos] = rs;
          eleCfg[rs] = pos;
          eleNum[rs] = 1;
        }
      }
      GCAntiRequireConfig(eleIdx, eleCfg, eleNum, "initial", comm);
    } else {
      for (position = 0; position < Ncur; position++) {
        const int emptyCount = Nsite2 - position;
        const int kth = (int)(gen_rand32() % (uint32_t)emptyCount);
        const int rs = GCFindKthEmpty(eleNum, Nsite2, kth);
        eleIdx[position] = rs;
        eleCfg[rs] = position;
        eleNum[rs] = 1;
      }
    }
    MakeProjCnt(eleProjCnt, eleNum);
    status = CalculateMAllGC_fcmp(Ncur, eleIdx, qpStart, qpEnd);
    globalStatus = GCCollectiveRebuildStatus(status, comm);
    if (globalStatus != GC_MALL_OK) {
      if (rank == 0) {
        fprintf(stderr,
                "Error: GC initial rebuild failed (status=%d).\n",
                globalStatus);
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
    logIp = CalculateLogIPGC_fcmp(PfM, Ncur, eleNum, qpStart, qpEnd, comm);
    if (isfinite(creal(logIp)) && isfinite(cimag(logIp))) return 0;
    if (!(creal(logIp) == -INFINITY && isfinite(cimag(logIp)))) {
      if (rank == 0) {
        fprintf(stderr,
                "Error: GC initial configuration has a nonfinite projected "
                "amplitude.\n");
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
    if (rank == 0)
      GCComponentCountersOf(GCComponentThreadCache())->compositeZero++;
    if (Ncur == 0 || Ncur == Nsite2) {
      attempt = 1;
      break;
    }
  }
  if (rank == 0) {
    fprintf(stderr,
            "Error: GC initial configuration has zero projected amplitude "
            "after %d attempt(s) (%s, Ncur=%d). Change NGCInitNelec or "
            "check the TransSym weights.\n",
            attempt, GCAntiMode() ? "anti-parallel mode" : "general mode",
            Ncur);
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

/* Finite S>0: the pair-function generation of the stored burn sample. Only
 * the Slater parameters can create zeros (Proj factors are exponentials). */
static uint64_t GCBurnSlaterHash = 0;

static uint64_t GCSlaterHash(void) {
  const unsigned char *byte = (const unsigned char *)Slater;
  const size_t count = (size_t)NSlater * sizeof(*Slater);
  uint64_t hash = UINT64_C(1469598103934665603);
  size_t k;
  for (k = 0; Slater != NULL && k < count; ++k) {
    hash ^= byte[k];
    hash *= UINT64_C(1099511628211);
  }
  return hash;
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
  if (GCFiniteSpinActive()) GCBurnSlaterHash = GCSlaterHash();
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
  GCRequireSpinSupport(comm);
  RequestWorkSpaceInt(NProj);
  projCntNew = GetWorkSpaceInt(NProj);

  if (BurnFlag == 0) {
    (void)makeInitialSampleGC(TmpEleIdx, TmpEleCfg, TmpEleNum,
                              TmpEleProjCnt, qpStart, qpEnd, comm);
  } else {
    int rebuildStatus;
    copyFromBurnSampleGC(TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt);
    GCWriteStateDump(stateDump, "RESTORE", -1, TmpEleIdx, TmpEleCfg,
                     TmpEleNum, TmpEleProjCnt, NULL);
    if (GCAntiMode()) {
      GCAntiRequireConfig(TmpEleIdx, TmpEleCfg, TmpEleNum, "restored", comm);
    }
    rebuildStatus = CalculateMAllGC_fcmp(Ncur, TmpEleIdx, qpStart, qpEnd);
    rebuildStatus = GCCollectiveRebuildStatus(rebuildStatus, comm);
    if (rebuildStatus != GC_MALL_OK) {
      if (rank == 0)
        fprintf(stderr,
                "Error: GC burn-sample rebuild failed (status=%d).\n",
                rebuildStatus);
      MPI_Abort(comm, EXIT_FAILURE);
    }
  }
  logIpOld = CalculateLogIPGC_fcmp(PfM, Ncur, TmpEleNum, qpStart, qpEnd,
                                   comm);
  if (!isfinite(creal(logIpOld)) || !isfinite(cimag(logIpOld))) {
    if (!(creal(logIpOld) == -INFINITY && isfinite(cimag(logIpOld)))) {
      if (rank == 0)
        fprintf(stderr,
                "Error: GC restored configuration has a nonfinite "
                "amplitude.\n");
      MPI_Abort(comm, EXIT_FAILURE);
    }
    /* An SR update may enter the pure-triplet manifold and invalidate the
     * saved N. Restart in a permitted sector and perform the full warm-up.
     * Initial user input is checked separately by makeInitialSampleGC. */
    if (GCFiniteSpinActive()) {
      /* Same pair-function generation: the stored state was nonzero, so a
       * zero now is a state or numerical anomaly, not a support change. */
      if (GCSlaterHash() == GCBurnSlaterHash) {
        if (rank == 0)
          fprintf(stderr,
                  "Error: GC total-spin S=%d restored state is zero without "
                  "a pair-function change.\n",
                  NSPStot);
        MPI_Abort(comm, EXIT_FAILURE);
      }
      GCFiniteSpinStart(Ncur, TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt,
                        qpStart, qpEnd, comm);
    } else {
      if (GCForbiddenSpinSector(Ncur)) {
        Ncur = NGCInitNelec;
        if (GCForbiddenSpinSector(Ncur)) Ncur = 0;
        if (Ncur == 0 && GCProjectionSectorZero(0, NULL)) Ncur = 4;
        if (Ncur > NsizeMax) {
          if (rank == 0)
            fprintf(stderr,
                    "Error: no GC spin sector available for restart.\n");
          MPI_Abort(comm, EXIT_FAILURE);
        }
      }
      (void)makeInitialSampleGC(TmpEleIdx, TmpEleCfg, TmpEleNum,
                                TmpEleProjCnt, qpStart, qpEnd, comm);
    }
    logIpOld = CalculateLogIPGC_fcmp(PfM, Ncur, TmpEleNum, qpStart, qpEnd,
                                     comm);
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
        logIpOld = CalculateLogIPGC_fcmp(PfM, Ncur, TmpEleNum, qpStart, qpEnd,
                                         comm);
        GCRequireRebuiltAmplitude(logIpOld, "periodic", comm);
        nAccept = 0;
      }
    }
    if (outStep >= nOutStep - NVMCSample) {
      const int sample = outStep - (nOutStep - NVMCSample);
      if (GCAntiMode() && !GCAntiValidateConfig(TmpEleIdx, TmpEleCfg,
                                                 TmpEleNum, Nsite, Ncur)) {
        GCAntiAbortState("sampled", comm);
      }
      saveEleConfigGC(sample, logIpOld, TmpEleIdx, TmpEleCfg, TmpEleNum,
                      TmpEleProjCnt, Ncur);
      GCWriteStateDump(stateDump, "SAMPLE", sample, TmpEleIdx, TmpEleCfg,
                       TmpEleNum, TmpEleProjCnt, &logIpOld);
    }
  }
  copyToBurnSampleGC(TmpEleIdx, TmpEleCfg, TmpEleNum, TmpEleProjCnt);
  GCWriteStateDump(stateDump, "STORE", -1, TmpEleIdx, TmpEleCfg,
                   TmpEleNum, TmpEleProjCnt, NULL);
  BurnFlag = 1;
  if (stateDump != NULL) fclose(stateDump);
  ReleaseWorkSpaceInt();
  free(pfMNew);
}

#endif
