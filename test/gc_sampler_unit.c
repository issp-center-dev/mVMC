#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

extern int omp_get_thread_num(void);
extern int omp_get_max_threads(void);

#include "global.h"
#include "blas_externs.h"
#include "SFMT.h"

#include "../src/mVMC/gc_size.c"
#include "../src/mVMC/workspace.c"
#include "../src/mVMC/projection.c"
#include "../src/mVMC/gc_config.c"
#include "../src/mVMC/matrix_gc.c"
#include "../src/mVMC/pfupdate_gc.c"
#include "../src/mVMC/splitloop.c"

double complex CalculateLogIP_fcmp(double complex *const pfM,
                                   const int qpStart, const int qpEnd,
                                   MPI_Comm comm) {
  double complex ip = 0.0;
  double complex reduced;
  int qpidx;
  int size;
  MPI_Comm_size(comm, &size);
  for (qpidx = 0; qpidx < qpEnd - qpStart; qpidx++) {
    ip += QPFullWeight[qpStart + qpidx] * pfM[qpidx];
  }
  if (size > 1) {
    MPI_Allreduce(&ip, &reduced, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);
    ip = reduced;
  }
  return clog(ip);
}

#include "../src/mVMC/gc_antiparallel.c"
#include "../src/mVMC/vmcmake_gc.c"

#define ORBITALS 4
#define STATE_COUNT 8
#define CHAIN_COUNT 32
#define CHAIN_WARMUP 10000
#define CHAIN_SAMPLES 100000

static const unsigned int evenStates[STATE_COUNT] = {
    0U, 3U, 5U, 6U, 9U, 10U, 12U, 15U};
static int failures = 0;

#define CHECK(condition, ...)                                                  \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "GCSampler_Unit FAIL: ");                              \
      fprintf(stderr, __VA_ARGS__);                                            \
      fprintf(stderr, "\n");                                                 \
      failures++;                                                              \
    }                                                                          \
  } while (0)

typedef struct {
  int ncur;
  int eleIdx[ORBITALS];
  int eleCfg[ORBITALS];
  int eleNum[ORBITALS];
  int eleProjCnt[1];
  double complex inv[ORBITALS * ORBITALS];
  double complex pf;
  double complex logIp;
} SamplerState;

static int popcount(unsigned int value) {
  int count = 0;
  while (value != 0U) {
    count += (int)(value & 1U);
    value >>= 1;
  }
  return count;
}

static int state_index(const unsigned int mask) {
  int i;
  for (i = 0; i < STATE_COUNT; i++) {
    if (evenStates[i] == mask) return i;
  }
  return -1;
}

static unsigned int current_mask(const int *eleNum) {
  unsigned int mask = 0U;
  int rs;
  for (rs = 0; rs < ORBITALS; rs++) {
    if (eleNum[rs] != 0) mask |= 1U << rs;
  }
  return mask;
}

static int close_complex(const double complex actual,
                         const double complex expected) {
  return cabs(actual - expected) <= 2.0e-9 * (1.0 + cabs(expected));
}

static void test_collective_rebuild_status(void) {
  int rank = 0;
  int size = 1;
  int localStatus;
  int globalStatus;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  localStatus = (size > 1 && rank == 0) ? GC_MALL_GETRF : GC_MALL_OK;
  globalStatus = GCCollectiveRebuildStatus(localStatus, MPI_COMM_WORLD);
  CHECK(globalStatus == (size > 1 ? GC_MALL_GETRF : GC_MALL_OK),
        "collective rebuild status mismatch: rank=%d size=%d got=%d",
        rank, size, globalStatus);
  globalStatus = GCCollectiveRebuildStatus(GC_MALL_OK, MPI_COMM_WORLD);
  CHECK(globalStatus == GC_MALL_OK,
        "collective rebuild success changed: rank=%d got=%d", rank,
        globalStatus);
}

static void fill_slater(void) {
  int row;
  memset(SlaterElm, 0, ORBITALS * ORBITALS * sizeof(*SlaterElm));
  for (row = 0; row < ORBITALS; row++) {
    int column;
    for (column = row + 1; column < ORBITALS; column++) {
      const double complex value =
          (0.43 + 0.21 * (row + 1) + 0.17 * (column + 1) +
           0.031 * (row + 1) * (column + 1)) +
          (0.09 * (row + 1) - 0.057 * (column + 1) + 0.023) * I;
      SlaterElm[(size_t)row * ORBITALS + (size_t)column] = value;
      SlaterElm[(size_t)column * ORBITALS + (size_t)row] = -value;
    }
  }
}

static void setup_globals(void) {
  /* The general-pair fixtures must not depend on zero-initialized globals. */
  FlagGrandCanonical = 1;
  iFlgOrbitalGeneral = 1;
  TwoSz = -1;
  Nsite = 2;
  Nsite2 = ORBITALS;
  Nsize = 2;
  NsizeMax = ORBITALS;
  NQPFull = 1;
  NProj = 1;
  NGutzwillerIdx = 2;
  NJastrowIdx = 0;
  NSpinJastrowIdx = 0;
  NDoublonHolon2siteIdx = 0;
  NDoublonHolon4siteIdx = 0;
  LapackLWork = 128;
  NThread = omp_get_max_threads();
  initializeWorkSpaceAll();
  SlaterElm = malloc(ORBITALS * ORBITALS * sizeof(*SlaterElm));
  InvM = malloc(ORBITALS * ORBITALS * sizeof(*InvM));
  PfM = malloc(sizeof(*PfM));
  QPFullWeight = malloc(sizeof(*QPFullWeight));
  Proj = malloc(sizeof(*Proj));
  GutzwillerIdx = malloc(Nsite * sizeof(*GutzwillerIdx));
  CHECK(SlaterElm != NULL && InvM != NULL && PfM != NULL &&
            QPFullWeight != NULL && Proj != NULL && GutzwillerIdx != NULL,
        "global allocation failed");
  if (failures != 0) exit(EXIT_FAILURE);
  QPFullWeight[0] = 1.0;
  Proj[0] = -0.37;
  GutzwillerIdx[0] = 0;
  GutzwillerIdx[1] = 0;
  fill_slater();
}

static double complex set_state(const unsigned int mask, int *eleIdx,
                                int *eleCfg, int *eleNum,
                                int *eleProjCnt) {
  int count = 0;
  int rs;
  for (rs = 0; rs < ORBITALS; rs++) {
    eleCfg[rs] = -1;
    eleNum[rs] = (int)((mask >> rs) & 1U);
    eleIdx[rs] = 70 + rs;
  }
  for (rs = 0; rs < ORBITALS; rs++) {
    if (eleNum[rs] != 0) {
      eleIdx[count] = rs;
      eleCfg[rs] = count;
      count++;
    }
  }
  Ncur = count;
  MakeProjCnt(eleProjCnt, eleNum);
  memset(InvM, 0x5a, ORBITALS * ORBITALS * sizeof(*InvM));
  CHECK(CalculateMAllGC_fcmp(Ncur, eleIdx, 0, 1) == GC_MALL_OK,
        "state mask=%u rebuild", mask);
  return CalculateLogIP_fcmp(PfM, 0, 1, MPI_COMM_SELF);
}

static void snapshot(SamplerState *state, const int *eleIdx,
                     const int *eleCfg, const int *eleNum,
                     const int *eleProjCnt, const double complex logIp) {
  state->ncur = Ncur;
  memcpy(state->eleIdx, eleIdx, sizeof(state->eleIdx));
  memcpy(state->eleCfg, eleCfg, sizeof(state->eleCfg));
  memcpy(state->eleNum, eleNum, sizeof(state->eleNum));
  memcpy(state->eleProjCnt, eleProjCnt, sizeof(state->eleProjCnt));
  memcpy(state->inv, InvM, sizeof(state->inv));
  state->pf = PfM[0];
  state->logIp = logIp;
}

static void check_snapshot(const SamplerState *state, const int *eleIdx,
                           const int *eleCfg, const int *eleNum,
                           const int *eleProjCnt,
                           const double complex logIp,
                           const char *label) {
  CHECK(Ncur == state->ncur, "%s Ncur changed", label);
  CHECK(memcmp(eleIdx, state->eleIdx, sizeof(state->eleIdx)) == 0,
        "%s eleIdx changed", label);
  CHECK(memcmp(eleCfg, state->eleCfg, sizeof(state->eleCfg)) == 0,
        "%s eleCfg changed", label);
  CHECK(memcmp(eleNum, state->eleNum, sizeof(state->eleNum)) == 0,
        "%s eleNum changed", label);
  CHECK(memcmp(eleProjCnt, state->eleProjCnt,
               sizeof(state->eleProjCnt)) == 0,
        "%s eleProjCnt changed", label);
  CHECK(memcmp(InvM, state->inv, sizeof(state->inv)) == 0,
        "%s InvM changed", label);
  CHECK(PfM[0] == state->pf, "%s PfM changed", label);
  CHECK(logIp == state->logIp, "%s logIp changed", label);
}

static void check_fast_against_rebuild(const int *eleIdx,
                                       const double complex fastPf,
                                       const double complex *fastInv,
                                       const char *label) {
  const int ncur = Ncur;
  int row;
  CHECK(CalculateMAllGC_fcmp(ncur, eleIdx, 0, 1) == GC_MALL_OK,
        "%s accepted rebuild failed", label);
  CHECK(close_complex(fastPf, PfM[0]), "%s accepted Pf mismatch", label);
  for (row = 0; row < ncur; row++) {
    int column;
    for (column = 0; column < ncur; column++) {
      const size_t index = (size_t)row * ORBITALS + (size_t)column;
      CHECK(close_complex(fastInv[index], InvM[index]),
            "%s accepted inverse row=%d col=%d", label, row, column);
    }
  }
}

static void exercise_transaction(const enum GCMoveClass moveClass,
                                 const unsigned int mask, const int arg0,
                                 const int arg1, const char *label) {
  int eleIdx[ORBITALS];
  int eleCfg[ORBITALS];
  int eleNum[ORBITALS];
  int eleProjCnt[1];
  int projCntNew[1];
  double complex pfMNew[1];
  double complex logIp;
  SamplerState oldState;
  double complex fastInv[ORBITALS * ORBITALS];
  double complex fastPf;

  logIp = set_state(mask, eleIdx, eleCfg, eleNum, eleProjCnt);
  snapshot(&oldState, eleIdx, eleCfg, eleNum, eleProjCnt, INFINITY + 0.0 * I);
  logIp = oldState.logIp;
  CHECK(GCAttemptMove(moveClass, arg0, arg1, 0.5, eleIdx, eleCfg, eleNum,
                      eleProjCnt, &logIp, pfMNew, projCntNew, 0, 1,
                      MPI_COMM_SELF) == 0,
        "%s deterministic reject accepted", label);
  check_snapshot(&oldState, eleIdx, eleCfg, eleNum, eleProjCnt, logIp, label);

  logIp = set_state(mask, eleIdx, eleCfg, eleNum, eleProjCnt);
  CHECK(GCAttemptMove(moveClass, arg0, arg1, 0.0, eleIdx, eleCfg, eleNum,
                      eleProjCnt, &logIp, pfMNew, projCntNew, 0, 1,
                      MPI_COMM_SELF) == 1,
        "%s deterministic accept rejected", label);
  fastPf = PfM[0];
  memcpy(fastInv, InvM, sizeof(fastInv));
  check_fast_against_rebuild(eleIdx, fastPf, fastInv, label);
}

static void test_transactions(void) {
  exercise_transaction(GC_MOVE_HOP, 3U, 0, 2, "hop transaction");
  exercise_transaction(GC_MOVE_ADD, 0U, 0, 1, "add transaction");
  exercise_transaction(GC_MOVE_REMOVE, 15U, 0, 1,
                       "remove transaction");
}

static double pfaffian_mask(const unsigned int mask);
static void test_exact_detailed_balance(void);

/*
 * Anti-periodic-style pair matrix: parameter 1 is shared by a +1 and a -1
 * upper-triangle row and parameter 2 enters with -1 (rows are I, J,
 * parameter, sign), so F = 2 s f.
 */
static void fill_slater_signed(void) {
  static const int rows[6][4] = {
      {0, 1, 0, 1}, {2, 3, 0, 1}, {0, 2, 1, 1},
      {1, 3, 1, -1}, {0, 3, 2, -1}, {1, 2, 3, 1}};
  const double complex parameters[4] = {
      0.61 - 0.22 * I, -0.37 + 0.45 * I, 0.52 + 0.31 * I, -0.29 - 0.48 * I};
  int row;
  memset(SlaterElm, 0, ORBITALS * ORBITALS * sizeof(*SlaterElm));
  for (row = 0; row < 6; row++) {
    const double complex value =
        2.0 * (double)rows[row][3] * parameters[rows[row][2]];
    SlaterElm[(size_t)rows[row][0] * ORBITALS + (size_t)rows[row][1]] = value;
    SlaterElm[(size_t)rows[row][1] * ORBITALS + (size_t)rows[row][0]] = -value;
  }
}

static void exercise_signed_transaction(const enum GCMoveClass moveClass,
                                        const unsigned int mask,
                                        const int arg0, const int arg1,
                                        const unsigned int expectedMask,
                                        const char *label) {
  int eleIdx[ORBITALS];
  int eleCfg[ORBITALS];
  int eleNum[ORBITALS];
  int eleProjCnt[1];
  int projCntNew[1];
  double complex pfMNew[1];
  double complex logIp;
  double complex fastInv[ORBITALS * ORBITALS];
  double complex fastPf;
  double complex fastLogIp;
  SamplerState oldState;

  logIp = set_state(mask, eleIdx, eleCfg, eleNum, eleProjCnt);
  snapshot(&oldState, eleIdx, eleCfg, eleNum, eleProjCnt, logIp);
  CHECK(GCAttemptMove(moveClass, arg0, arg1, 1.0e300, eleIdx, eleCfg, eleNum,
                      eleProjCnt, &logIp, pfMNew, projCntNew, 0, 1,
                      MPI_COMM_SELF) == 0,
        "%s forced reject accepted", label);
  check_snapshot(&oldState, eleIdx, eleCfg, eleNum, eleProjCnt, logIp, label);

  logIp = set_state(mask, eleIdx, eleCfg, eleNum, eleProjCnt);
  CHECK(GCAttemptMove(moveClass, arg0, arg1, 0.0, eleIdx, eleCfg, eleNum,
                      eleProjCnt, &logIp, pfMNew, projCntNew, 0, 1,
                      MPI_COMM_SELF) == 1,
        "%s forced accept rejected", label);
  CHECK(current_mask(eleNum) == expectedMask,
        "%s produced mask %u, expected %u", label, current_mask(eleNum),
        expectedMask);
  fastPf = PfM[0];
  fastLogIp = logIp;
  memcpy(fastInv, InvM, sizeof(fastInv));
  CHECK(fabs(cabs(fastPf) - pfaffian_mask(expectedMask)) <=
            2.0e-9 * (1.0 + pfaffian_mask(expectedMask)),
        "%s |Pf| %.17g differs from independent %.17g", label, cabs(fastPf),
        pfaffian_mask(expectedMask));
  check_fast_against_rebuild(eleIdx, fastPf, fastInv, label);
  CHECK(close_complex(cexp(fastLogIp),
                      cexp(CalculateLogIP_fcmp(PfM, 0, 1, MPI_COMM_SELF))),
        "%s accepted logIp differs from rebuild", label);
}

/* Explicit hop, pair-add and pair-remove transactions on the signed matrix;
 * the original fixture is restored afterwards for the chain tests. */
static void test_signed_transactions(void) {
  exercise_signed_transaction(GC_MOVE_HOP, 9U, 1, 2, 5U,
                              "signed hop (0,3)->(0,2)");
  exercise_signed_transaction(GC_MOVE_HOP, 12U, 0, 1, 10U,
                              "signed hop (2,3)->(1,3)");
  exercise_signed_transaction(GC_MOVE_HOP, 3U, 0, 3, 10U,
                              "signed hop (0,1)->(1,3)");
  exercise_signed_transaction(GC_MOVE_ADD, 0U, 0, 3, 9U,
                              "signed add (0,3)");
  exercise_signed_transaction(GC_MOVE_ADD, 3U, 2, 3, 15U,
                              "signed add (2,3) to (0,1)");
  exercise_signed_transaction(GC_MOVE_ADD, 6U, 3, 0, 15U,
                              "signed add (3,0) to (1,2)");
  exercise_signed_transaction(GC_MOVE_REMOVE, 15U, 0, 3, 6U,
                              "signed remove slots (0,3)");
  exercise_signed_transaction(GC_MOVE_REMOVE, 15U, 1, 3, 5U,
                              "signed remove slots (1,3)");
  exercise_signed_transaction(GC_MOVE_REMOVE, 9U, 0, 1, 0U,
                              "signed remove to vacuum");
}

static void test_signed_matrix(void) {
  double complex original[ORBITALS * ORBITALS];
  memcpy(original, SlaterElm, sizeof(original));
  fill_slater_signed();
  test_signed_transactions();
  test_exact_detailed_balance();
  memcpy(SlaterElm, original, sizeof(original));
}

static double pfaffian_mask(const unsigned int mask) {
  int occupied[ORBITALS];
  int ncur = 0;
  int rs;
  for (rs = 0; rs < ORBITALS; rs++) {
    if ((mask & (1U << rs)) != 0U) occupied[ncur++] = rs;
  }
  if (ncur == 0) return 1.0;
  if (ncur == 2) {
    return cabs(-SlaterElm[(size_t)occupied[0] * ORBITALS +
                           (size_t)occupied[1]]);
  }
  {
    const double complex x01 =
        -SlaterElm[(size_t)occupied[0] * ORBITALS + occupied[1]];
    const double complex x02 =
        -SlaterElm[(size_t)occupied[0] * ORBITALS + occupied[2]];
    const double complex x03 =
        -SlaterElm[(size_t)occupied[0] * ORBITALS + occupied[3]];
    const double complex x12 =
        -SlaterElm[(size_t)occupied[1] * ORBITALS + occupied[2]];
    const double complex x13 =
        -SlaterElm[(size_t)occupied[1] * ORBITALS + occupied[3]];
    const double complex x23 =
        -SlaterElm[(size_t)occupied[2] * ORBITALS + occupied[3]];
    return cabs(x01 * x23 - x02 * x13 + x03 * x12);
  }
}

static int double_occupancy(const unsigned int mask) {
  int count = 0;
  int site;
  for (site = 0; site < Nsite; site++) {
    count += (int)(((mask >> site) & 1U) &
                   ((mask >> (site + Nsite)) & 1U));
  }
  return count;
}

static void exact_probabilities(double *probability) {
  double normalization = 0.0;
  int i;
  for (i = 0; i < STATE_COUNT; i++) {
    const double pf = pfaffian_mask(evenStates[i]);
    probability[i] =
        exp(2.0 * creal(Proj[0]) * double_occupancy(evenStates[i])) *
        pf * pf;
    normalization += probability[i];
  }
  for (i = 0; i < STATE_COUNT; i++) probability[i] /= normalization;
}

static double choose_two(const int count) {
  return 0.5 * (double)count * (double)(count - 1);
}

static double class_probability(const enum GCMoveClass moveClass,
                                const int ncur) {
  const int hop = ncur > 0 && ncur < ORBITALS;
  const int add = ORBITALS - ncur >= 2;
  const int remove = ncur >= 2;
  const double normalization = 0.5 * hop + 0.25 * add + 0.25 * remove;
  if (moveClass == GC_MOVE_HOP) return hop ? 0.5 / normalization : 0.0;
  if (moveClass == GC_MOVE_ADD) return add ? 0.25 / normalization : 0.0;
  return remove ? 0.25 / normalization : 0.0;
}

static void test_selector(void) {
  CHECK(GCSelectMoveClass(0.0, 0, ORBITALS) == GC_MOVE_ADD,
        "vacuum selector");
  CHECK(GCSelectMoveClass(0.999999, ORBITALS, ORBITALS) == GC_MOVE_REMOVE,
        "full selector");
  CHECK(GCSelectMoveClass(0.49, 2, ORBITALS) == GC_MOVE_HOP,
        "interior hop selector");
  CHECK(GCSelectMoveClass(0.50, 2, ORBITALS) == GC_MOVE_ADD,
        "interior add selector");
  CHECK(GCSelectMoveClass(0.75, 2, ORBITALS) == GC_MOVE_REMOVE,
        "interior remove selector");
}

static void test_exact_detailed_balance(void) {
  double probability[STATE_COUNT];
  int mutationKilled = 0;
  int ix;
  exact_probabilities(probability);
  for (ix = 0; ix < STATE_COUNT; ix++) {
    int iy;
    for (iy = 0; iy < STATE_COUNT; iy++) {
      const unsigned int removed = evenStates[ix] & ~evenStates[iy];
      const unsigned int added = evenStates[iy] & ~evenStates[ix];
      const int nx = popcount(evenStates[ix]);
      const int ny = popcount(evenStates[iy]);
      double qxy = 0.0;
      double qyx = 0.0;
      double proposalRatio = 1.0;
      if (nx == ny && popcount(removed) == 1 && popcount(added) == 1) {
        qxy = class_probability(GC_MOVE_HOP, nx) /
              ((double)nx * (double)(ORBITALS - nx));
        qyx = qxy;
      } else if (ny == nx + 2 && removed == 0U && popcount(added) == 2) {
        const double pAdd = class_probability(GC_MOVE_ADD, nx);
        const double pRemove = class_probability(GC_MOVE_REMOVE, ny);
        qxy = pAdd / choose_two(ORBITALS - nx);
        qyx = pRemove / choose_two(ny);
        proposalRatio = GCProposalRatioAdd(nx, ORBITALS, pAdd, pRemove);
      } else {
        continue;
      }
      {
        const double targetRatio = probability[iy] / probability[ix];
        const double acceptanceXY = fmin(1.0, targetRatio * proposalRatio);
        const double acceptanceYX =
            fmin(1.0, 1.0 / (targetRatio * proposalRatio));
        const double lhs = probability[ix] * qxy * acceptanceXY;
        const double rhs = probability[iy] * qyx * acceptanceYX;
        const double mutatedLhs =
            probability[ix] * qxy * fmin(1.0, targetRatio);
        const double mutatedRhs =
            probability[iy] * qyx * fmin(1.0, 1.0 / targetRatio);
        CHECK(fabs(lhs - rhs) < 2.0e-13,
              "detailed balance mask=%u -> %u lhs=%.17g rhs=%.17g",
              evenStates[ix], evenStates[iy], lhs, rhs);
        if (fabs(mutatedLhs - mutatedRhs) > 1.0e-8) mutationKilled = 1;
      }
    }
  }
  CHECK(mutationKilled, "proposalRatio=1 mutation survived exact gate");
}

static void test_burn_and_save(void) {
  int eleIdx[ORBITALS];
  int eleCfg[ORBITALS];
  int eleNum[ORBITALS];
  int eleProjCnt[1];
  int idxSnapshot[ORBITALS];
  int cfgSnapshot[ORBITALS];
  int numSnapshot[ORBITALS];
  int projSnapshot[1];
  double complex logIp =
      set_state(5U, eleIdx, eleCfg, eleNum, eleProjCnt);
  size_t i;
  BurnEleIdx = malloc(ORBITALS * sizeof(*BurnEleIdx));
  BurnEleCfg = malloc(ORBITALS * sizeof(*BurnEleCfg));
  BurnEleNum = malloc(ORBITALS * sizeof(*BurnEleNum));
  BurnEleProjCnt = malloc(sizeof(*BurnEleProjCnt));
  EleIdx = malloc(2U * ORBITALS * sizeof(*EleIdx));
  EleCfg = malloc(2U * ORBITALS * sizeof(*EleCfg));
  EleNum = malloc(2U * ORBITALS * sizeof(*EleNum));
  EleProjCnt = malloc(2U * sizeof(*EleProjCnt));
  EleNumSample = malloc(2U * sizeof(*EleNumSample));
  logSqPfFullSlater = malloc(2U * sizeof(*logSqPfFullSlater));
  CHECK(BurnEleIdx != NULL && BurnEleCfg != NULL && BurnEleNum != NULL &&
            BurnEleProjCnt != NULL && EleIdx != NULL && EleCfg != NULL &&
            EleNum != NULL && EleProjCnt != NULL && EleNumSample != NULL &&
            logSqPfFullSlater != NULL,
        "burn/save allocation");
  eleIdx[2] = 91;
  eleIdx[3] = 92;
  memcpy(idxSnapshot, eleIdx, sizeof(idxSnapshot));
  memcpy(cfgSnapshot, eleCfg, sizeof(cfgSnapshot));
  memcpy(numSnapshot, eleNum, sizeof(numSnapshot));
  memcpy(projSnapshot, eleProjCnt, sizeof(projSnapshot));
  copyToBurnSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt);
  memset(eleIdx, 0xa1, sizeof(eleIdx));
  memset(eleCfg, 0xa2, sizeof(eleCfg));
  memset(eleNum, 0xa3, sizeof(eleNum));
  memset(eleProjCnt, 0xa4, sizeof(eleProjCnt));
  Ncur = 0;
  copyFromBurnSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt);
  CHECK(Ncur == 2, "burn Ncur restore");
  CHECK(memcmp(eleIdx, idxSnapshot, sizeof(idxSnapshot)) == 0,
        "burn eleIdx restore");
  CHECK(memcmp(eleCfg, cfgSnapshot, sizeof(cfgSnapshot)) == 0,
        "burn eleCfg restore");
  CHECK(memcmp(eleNum, numSnapshot, sizeof(numSnapshot)) == 0,
        "burn eleNum restore");
  CHECK(memcmp(eleProjCnt, projSnapshot, sizeof(projSnapshot)) == 0,
        "burn proj restore");
  for (i = 0; i < 2U * ORBITALS; i++) {
    EleIdx[i] = -9;
    EleCfg[i] = -9;
    EleNum[i] = -9;
  }
  saveEleConfigGC(1, logIp, eleIdx, eleCfg, eleNum, eleProjCnt, Ncur);
  CHECK(memcmp(EleIdx + ORBITALS, eleIdx, sizeof(idxSnapshot)) == 0,
        "saved eleIdx capacity stride");
  CHECK(memcmp(EleCfg + ORBITALS, eleCfg, sizeof(cfgSnapshot)) == 0,
        "saved eleCfg stride");
  CHECK(memcmp(EleNum + ORBITALS, eleNum, sizeof(numSnapshot)) == 0,
        "saved eleNum stride");
  CHECK(EleNumSample[1] == Ncur, "saved ncur");
  CHECK(fabs(logSqPfFullSlater[1] -
             2.0 * (LogProjVal(eleProjCnt) + creal(logIp))) < 1.0e-14,
        "saved log weight");
}

static void test_production_chains(void) {
  double exact[STATE_COUNT];
  double frequency[CHAIN_COUNT][STATE_COUNT];
  long long totalCounts[STATE_COUNT] = {0};
  int saw02 = 0;
  int saw24 = 0;
  int chain;
  exact_probabilities(exact);
  for (chain = 0; chain < CHAIN_COUNT; chain++) {
    int eleIdx[ORBITALS];
    int eleCfg[ORBITALS];
    int eleNum[ORBITALS];
    int eleProjCnt[1];
    int projCntNew[1];
    double complex pfMNew[1];
    double complex logIp;
    long long counts[STATE_COUNT] = {0};
    int previousNcur;
    int step;
    init_gen_rand((uint32_t)(1013 + 7919 * chain));
    logIp = set_state(evenStates[chain % STATE_COUNT], eleIdx, eleCfg,
                      eleNum, eleProjCnt);
    previousNcur = Ncur;
    for (step = 0; step < CHAIN_WARMUP; step++) {
      (void)GCMakeOneStep(eleIdx, eleCfg, eleNum, eleProjCnt, &logIp,
                          pfMNew, projCntNew, 0, 1, MPI_COMM_SELF);
      if ((previousNcur == 0 && Ncur == 2) ||
          (previousNcur == 2 && Ncur == 0)) saw02 = 1;
      if ((previousNcur == 2 && Ncur == 4) ||
          (previousNcur == 4 && Ncur == 2)) saw24 = 1;
      previousNcur = Ncur;
    }
    for (step = 0; step < CHAIN_SAMPLES; step++) {
      int separation;
      for (separation = 0; separation < ORBITALS; separation++) {
        (void)GCMakeOneStep(eleIdx, eleCfg, eleNum, eleProjCnt, &logIp,
                            pfMNew, projCntNew, 0, 1, MPI_COMM_SELF);
        if ((previousNcur == 0 && Ncur == 2) ||
            (previousNcur == 2 && Ncur == 0)) saw02 = 1;
        if ((previousNcur == 2 && Ncur == 4) ||
            (previousNcur == 4 && Ncur == 2)) saw24 = 1;
        previousNcur = Ncur;
      }
      {
        const int index = state_index(current_mask(eleNum));
        CHECK(index >= 0, "chain produced odd state");
        if (index >= 0) counts[index]++;
      }
    }
    for (step = 0; step < STATE_COUNT; step++) {
      frequency[chain][step] =
          (double)counts[step] / (double)CHAIN_SAMPLES;
      totalCounts[step] += counts[step];
    }
  }
  CHECK(saw02, "production chains did not observe 0<->2 transition");
  CHECK(saw24, "production chains did not observe 2<->4 transition");
  for (chain = 0; chain < STATE_COUNT; chain++) {
    double mean = 0.0;
    double variance = 0.0;
    double standardError;
    double tolerance;
    int sample;
    for (sample = 0; sample < CHAIN_COUNT; sample++) {
      mean += frequency[sample][chain];
    }
    mean /= CHAIN_COUNT;
    for (sample = 0; sample < CHAIN_COUNT; sample++) {
      const double delta = frequency[sample][chain] - mean;
      variance += delta * delta;
    }
    variance /= CHAIN_COUNT - 1;
    standardError = sqrt(variance / CHAIN_COUNT);
    tolerance = fmax(6.0 * standardError, 5.0e-4);
    CHECK(fabs(mean - exact[chain]) <= tolerance,
          "chain state=%u mean=%.17g exact=%.17g SE=%.17g tol=%.17g",
          evenStates[chain], mean, exact[chain], standardError, tolerance);
    CHECK((double)(CHAIN_COUNT * CHAIN_SAMPLES) * exact[chain] >= 100.0,
          "chain state=%u expected count below 100", evenStates[chain]);
    CHECK(totalCounts[chain] > 0, "chain state=%u never observed",
          evenStates[chain]);
  }
}


/* ------------------------------------------------------------------------
 * Anti-parallel (Sz=0) mode: FlagGrandCanonical=1, iFlgOrbitalGeneral=0.
 * Orbitals 0,1 are up and 2,3 down; the same-spin blocks of the pair matrix
 * vanish, so only the six balanced masks carry amplitude.
 * ---------------------------------------------------------------------- */
#define ANTI_STATE_COUNT 6
static const unsigned int antiStates[ANTI_STATE_COUNT] = {0U,  5U,  6U,
                                                          9U, 10U, 15U};

static int finite_close_complex(const double complex actual,
                                const double complex expected,
                                const double tolerance) {
  return isfinite(creal(actual)) && isfinite(cimag(actual)) &&
         isfinite(creal(expected)) && isfinite(cimag(expected)) &&
         isfinite(tolerance) && tolerance >= 0.0 &&
         cabs(actual - expected) <= tolerance;
}

static void fill_slater_antiparallel(void) {
  const double complex F[2][2] = {{0.71 + 0.12 * I, 0.34 - 0.23 * I},
                                  {-0.41 + 0.26 * I, 0.62 - 0.11 * I}};
  int i;
  memset(SlaterElm, 0, ORBITALS * ORBITALS * sizeof(*SlaterElm));
  for (i = 0; i < 2; i++) {
    int j;
    for (j = 0; j < 2; j++) {
      SlaterElm[(size_t)i * ORBITALS + (size_t)(2 + j)] = F[i][j];
      SlaterElm[(size_t)(2 + j) * ORBITALS + (size_t)i] = -F[i][j];
    }
  }
}

static void enter_antiparallel_mode(void) {
  FlagGrandCanonical = 1;
  iFlgOrbitalGeneral = 0;
  TwoSz = 0;
  fill_slater_antiparallel();
}

static int anti_state_index(const unsigned int mask) {
  int i;
  for (i = 0; i < ANTI_STATE_COUNT; i++) {
    if (antiStates[i] == mask) return i;
  }
  return -1;
}

static int spin_balanced(const unsigned int mask) {
  return popcount(mask & 3U) == popcount((mask >> 2) & 3U);
}

/* Like set_state, optionally storing the electrons in reverse order. */
static double complex set_state_order(const unsigned int mask,
                                      const int reversed, int *eleIdx,
                                      int *eleCfg, int *eleNum,
                                      int *eleProjCnt) {
  double complex logIp = set_state(mask, eleIdx, eleCfg, eleNum, eleProjCnt);
  if (reversed && Ncur > 1) {
    int k;
    for (k = 0; k < Ncur / 2; k++) {
      const int other = Ncur - 1 - k;
      const int temporary = eleIdx[k];
      eleIdx[k] = eleIdx[other];
      eleIdx[other] = temporary;
    }
    for (k = 0; k < Ncur; k++) eleCfg[eleIdx[k]] = k;
    CHECK(CalculateMAllGC_fcmp(Ncur, eleIdx, 0, 1) == GC_MALL_OK,
          "reversed state mask=%u rebuild", mask);
    logIp = CalculateLogIP_fcmp(PfM, 0, 1, MPI_COMM_SELF);
  }
  return logIp;
}

static void anti_exact_probabilities(double *probability) {
  double normalization = 0.0;
  int i;
  for (i = 0; i < ANTI_STATE_COUNT; i++) {
    const double pf = pfaffian_mask(antiStates[i]);
    probability[i] =
        exp(2.0 * creal(Proj[0]) * double_occupancy(antiStates[i])) * pf * pf;
    normalization += probability[i];
  }
  for (i = 0; i < ANTI_STATE_COUNT; i++) probability[i] /= normalization;
}

/* Guards run before any Pfaffian work: the candidate scratch keeps its
 * sentinel and the state is untouched even with an always-accept draw. */
static void expect_guard_reject(const enum GCMoveClass moveClass,
                                const unsigned int mask, const int arg0Orbital,
                                const int arg1, const int argIsPosition,
                                const char *label) {
  int eleIdx[ORBITALS];
  int eleCfg[ORBITALS];
  int eleNum[ORBITALS];
  int eleProjCnt[1];
  int projCntNew[1];
  double complex pfMNew[1];
  double complex logIp;
  SamplerState before;
  int arg0;
  int second;
  logIp = set_state(mask, eleIdx, eleCfg, eleNum, eleProjCnt);
  snapshot(&before, eleIdx, eleCfg, eleNum, eleProjCnt, logIp);
  arg0 = argIsPosition ? eleCfg[arg0Orbital] : arg0Orbital;
  second = (argIsPosition && moveClass == GC_MOVE_REMOVE) ? eleCfg[arg1] : arg1;
  if (moveClass == GC_MOVE_REMOVE && arg0 > second) {
    const int temporary = arg0;
    arg0 = second;
    second = temporary;
  }
  pfMNew[0] = 123.0 + 7.0 * I;
  CHECK(GCAttemptMove(moveClass, arg0, second, 0.0, eleIdx, eleCfg, eleNum,
                      eleProjCnt, &logIp, pfMNew, projCntNew, 0, 1,
                      MPI_COMM_SELF) == 0,
        "%s was accepted", label);
  check_snapshot(&before, eleIdx, eleCfg, eleNum, eleProjCnt, logIp, label);
  CHECK(pfMNew[0] == 123.0 + 7.0 * I, "%s reached Pf calculation", label);
}

static void test_anti_guards(void) {
  /* same-spin adds from the vacuum */
  expect_guard_reject(GC_MOVE_ADD, 0U, 0, 1, 0, "same-spin up add");
  expect_guard_reject(GC_MOVE_ADD, 0U, 2, 3, 0, "same-spin down add");
  /* Non-vacuum mask 5: occupied up0/down0; spin-changing hops. */
  expect_guard_reject(GC_MOVE_HOP, 5U, 2, 1, 1, "spin-changing hop down->up");
  expect_guard_reject(GC_MOVE_HOP, 5U, 0, 3, 1, "spin-changing hop up->down");
  /* Full filling: same-spin pair removals. */
  expect_guard_reject(GC_MOVE_REMOVE, 15U, 0, 1, 1, "same-spin up remove");
  expect_guard_reject(GC_MOVE_REMOVE, 15U, 2, 3, 1, "same-spin down remove");
}

static void test_anti_initialization(void) {
  int eleIdx[ORBITALS];
  int eleCfg[ORBITALS];
  int eleNum[ORBITALS];
  int eleProjCnt[1];
  int seen[16] = {0};
  int ncur;
  for (ncur = 0; ncur <= ORBITALS; ncur += 2) {
    int seed;
    for (seed = 0; seed < 64; seed++) {
      unsigned int mask;
      init_gen_rand((uint32_t)(4099 + 31 * seed + ncur));
      Ncur = ncur;
      CHECK(makeInitialSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt, 0, 1,
                                MPI_COMM_SELF) == 0,
            "initialization ncur=%d failed", ncur);
      mask = current_mask(eleNum);
      CHECK(Ncur == ncur, "initialization changed Ncur");
      CHECK(spin_balanced(mask) && popcount(mask) == ncur,
            "initialization ncur=%d produced mask %u", ncur, mask);
      CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, Nsite, Ncur) == 1,
            "initialization ncur=%d produced an inconsistent state", ncur);
      {
        /* The spin-resolved draw stores the m up electrons first. */
        int position;
        for (position = 0; position < ncur; position++) {
          CHECK((position < ncur / 2) == (eleIdx[position] < Nsite),
                "initialization ncur=%d position %d holds orbital %d", ncur,
                position, eleIdx[position]);
        }
      }
      if (mask < 16U) seen[mask] = 1;
    }
  }
  CHECK(seen[0] && seen[5] && seen[6] && seen[9] && seen[10] && seen[15],
        "initialization did not reach every balanced state");
}

/* Every hop/add/remove candidate of the spin-balanced proposal from every
 * balanced state, in both electron orders. */
typedef struct {
  enum GCMoveClass moveClass;
  int first;  /* orbital (hop: moved electron, remove: up electron) */
  int second; /* orbital (hop: target, add: down, remove: down electron) */
} AntiMove;

static int anti_moves(const unsigned int mask, AntiMove *moves) {
  int count = 0;
  int a;
  int b;
  for (a = 0; a < ORBITALS; a++) {
    for (b = 0; b < ORBITALS; b++) {
      const int sameSpin = (a / 2) == (b / 2);
      const int aOcc = (int)((mask >> a) & 1U);
      const int bOcc = (int)((mask >> b) & 1U);
      if (sameSpin && aOcc && !bOcc) {
        moves[count].moveClass = GC_MOVE_HOP;
        moves[count].first = a;
        moves[count].second = b;
        count++;
      }
      if (a < 2 && b >= 2 && !aOcc && !bOcc) {
        moves[count].moveClass = GC_MOVE_ADD;
        moves[count].first = a;
        moves[count].second = b;
        count++;
      }
      if (a < 2 && b >= 2 && aOcc && bOcc) {
        moves[count].moveClass = GC_MOVE_REMOVE;
        moves[count].first = a;
        moves[count].second = b;
        count++;
      }
    }
  }
  return count;
}

static unsigned int apply_anti_move(const unsigned int mask,
                                    const AntiMove *move) {
  if (move->moveClass == GC_MOVE_HOP) {
    return mask ^ (1U << move->first) ^ (1U << move->second);
  }
  if (move->moveClass == GC_MOVE_ADD) {
    return mask | (1U << move->first) | (1U << move->second);
  }
  return mask ^ (1U << move->first) ^ (1U << move->second);
}

static void move_arguments(const AntiMove *move, const int *eleCfg,
                           int *arg0, int *arg1) {
  if (move->moveClass == GC_MOVE_HOP) {
    *arg0 = eleCfg[move->first];
    *arg1 = move->second;
  } else if (move->moveClass == GC_MOVE_ADD) {
    *arg0 = move->first;
    *arg1 = move->second;
  } else {
    *arg0 = eleCfg[move->first];
    *arg1 = eleCfg[move->second];
    if (*arg0 > *arg1) {
      const int temporary = *arg0;
      *arg0 = *arg1;
      *arg1 = temporary;
    }
  }
}

static double anti_proposal_ratio(const AntiMove *move, const int ncur) {
  const int m = ncur / 2;
  if (move->moveClass == GC_MOVE_ADD) {
    return GCAntiRatioAdd(m, Nsite, class_probability(GC_MOVE_ADD, ncur),
                          class_probability(GC_MOVE_REMOVE, ncur + 2));
  }
  if (move->moveClass == GC_MOVE_REMOVE) {
    return GCAntiRatioRemove(m, Nsite,
                             class_probability(GC_MOVE_REMOVE, ncur),
                             class_probability(GC_MOVE_ADD, ncur - 2));
  }
  return 1.0;
}

static void test_anti_transactions(void) {
  int stateIndex;
  int boundaryProbes = 0;
  for (stateIndex = 0; stateIndex < ANTI_STATE_COUNT; stateIndex++) {
    const unsigned int mask = antiStates[stateIndex];
    AntiMove moves[32];
    const int count = anti_moves(mask, moves);
    int k;
    for (k = 0; k < count; k++) {
      int reversed;
      for (reversed = 0; reversed < 2; reversed++) {
        int eleIdx[ORBITALS];
        int eleCfg[ORBITALS];
        int eleNum[ORBITALS];
        int eleProjCnt[1];
        int projCntNew[1];
        double complex pfMNew[1];
        double complex logIp;
        double complex rebuiltLog;
        double complex fastInv[ORBITALS * ORBITALS];
        double complex fastPf;
        SamplerState before;
        const unsigned int target = apply_anti_move(mask, &moves[k]);
        char label[128];
        int arg0;
        int arg1;
        snprintf(label, sizeof(label), "anti move class=%d %u->%u order=%d",
                 (int)moves[k].moveClass, mask, target, reversed);
        CHECK(spin_balanced(target), "%s leaves Sz=0", label);

        /* Finite forced rejection keeps every piece of state. */
        logIp = set_state_order(mask, reversed, eleIdx, eleCfg, eleNum,
                                eleProjCnt);
        logIp += 1000.0;
        snapshot(&before, eleIdx, eleCfg, eleNum, eleProjCnt, logIp);
        move_arguments(&moves[k], eleCfg, &arg0, &arg1);
        CHECK(GCAttemptMove(moves[k].moveClass, arg0, arg1, 0.5, eleIdx,
                            eleCfg, eleNum, eleProjCnt, &logIp, pfMNew,
                            projCntNew, 0, 1, MPI_COMM_SELF) == 0,
              "%s forced reject accepted", label);
        check_snapshot(&before, eleIdx, eleCfg, eleNum, eleProjCnt, logIp,
                       label);

        /* draw=0 accepts; the fast update equals a full rebuild. */
        logIp = set_state_order(mask, reversed, eleIdx, eleCfg, eleNum,
                                eleProjCnt);
        move_arguments(&moves[k], eleCfg, &arg0, &arg1);
        CHECK(GCAttemptMove(moves[k].moveClass, arg0, arg1, 0.0, eleIdx,
                            eleCfg, eleNum, eleProjCnt, &logIp, pfMNew,
                            projCntNew, 0, 1, MPI_COMM_SELF) == 1,
              "%s forced accept rejected", label);
        CHECK(current_mask(eleNum) == target, "%s produced mask %u", label,
              current_mask(eleNum));
        CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, Nsite, Ncur) == 1,
              "%s left an inconsistent state", label);
        fastPf = PfM[0];
        memcpy(fastInv, InvM, sizeof(fastInv));
        check_fast_against_rebuild(eleIdx, fastPf, fastInv, label);
        rebuiltLog = CalculateLogIP_fcmp(PfM, 0, 1, MPI_COMM_SELF);
        CHECK(finite_close_complex(cexp(logIp), cexp(rebuiltLog),
                                   2.0e-9 * (1.0 + cabs(cexp(rebuiltLog)))),
              "%s accepted log differs from rebuild", label);

        /* The Metropolis boundary uses the spin-resolved proposal ratio. */
        {
          const double pfX = pfaffian_mask(mask);
          const double pfY = pfaffian_mask(target);
          const double weight =
              (pfY * pfY) / (pfX * pfX) *
              exp(2.0 * creal(Proj[0]) *
                  (double_occupancy(target) - double_occupancy(mask))) *
              anti_proposal_ratio(&moves[k], popcount(mask));
          if (weight < 1.0 && weight > 1.0e-6) {
            int accepted;
            boundaryProbes++;
            logIp = set_state_order(mask, reversed, eleIdx, eleCfg, eleNum,
                                    eleProjCnt);
            move_arguments(&moves[k], eleCfg, &arg0, &arg1);
            accepted = GCAttemptMove(moves[k].moveClass, arg0, arg1,
                                     weight * (1.0 - 1.0e-9), eleIdx, eleCfg,
                                     eleNum, eleProjCnt, &logIp, pfMNew,
                                     projCntNew, 0, 1, MPI_COMM_SELF);
            CHECK(accepted == 1, "%s rejected just below weight %.17g", label,
                  weight);
            logIp = set_state_order(mask, reversed, eleIdx, eleCfg, eleNum,
                                    eleProjCnt);
            move_arguments(&moves[k], eleCfg, &arg0, &arg1);
            accepted = GCAttemptMove(moves[k].moveClass, arg0, arg1,
                                     weight * (1.0 + 1.0e-9), eleIdx, eleCfg,
                                     eleNum, eleProjCnt, &logIp, pfMNew,
                                     projCntNew, 0, 1, MPI_COMM_SELF);
            CHECK(accepted == 0, "%s accepted just above weight %.17g", label,
                  weight);
          }
        }
      }
    }
  }
  CHECK(boundaryProbes > 10, "too few Metropolis boundary probes (%d)",
        boundaryProbes);
}

static void test_anti_chains(void) {
  double exact[ANTI_STATE_COUNT];
  double frequency[CHAIN_COUNT][ANTI_STATE_COUNT];
  long long totalCounts[ANTI_STATE_COUNT] = {0};
  int saw02 = 0;
  int saw24 = 0;
  int imbalanced = 0;
  int chain;
  anti_exact_probabilities(exact);
  for (chain = 0; chain < CHAIN_COUNT; chain++) {
    int eleIdx[ORBITALS];
    int eleCfg[ORBITALS];
    int eleNum[ORBITALS];
    int eleProjCnt[1];
    int projCntNew[1];
    double complex pfMNew[1];
    double complex logIp;
    long long counts[ANTI_STATE_COUNT] = {0};
    int previousNcur;
    int step;
    init_gen_rand((uint32_t)(2027 + 7919 * chain));
    /* Half of the chains start with down electrons stored first. */
    logIp = set_state_order(antiStates[chain % ANTI_STATE_COUNT],
                            (chain / ANTI_STATE_COUNT) % 2, eleIdx, eleCfg,
                            eleNum, eleProjCnt);
    previousNcur = Ncur;
    for (step = 0; step < CHAIN_WARMUP + CHAIN_SAMPLES; step++) {
      int separation;
      for (separation = 0; separation < ORBITALS; separation++) {
        (void)GCMakeOneStep(eleIdx, eleCfg, eleNum, eleProjCnt, &logIp,
                            pfMNew, projCntNew, 0, 1, MPI_COMM_SELF);
        if ((previousNcur == 0 && Ncur == 2) ||
            (previousNcur == 2 && Ncur == 0)) saw02 = 1;
        if ((previousNcur == 2 && Ncur == 4) ||
            (previousNcur == 4 && Ncur == 2)) saw24 = 1;
        if (!spin_balanced(current_mask(eleNum))) imbalanced++;
        previousNcur = Ncur;
      }
      if (step >= CHAIN_WARMUP) {
        const int index = anti_state_index(current_mask(eleNum));
        if (index >= 0) counts[index]++;
      }
    }
    for (step = 0; step < ANTI_STATE_COUNT; step++) {
      frequency[chain][step] = (double)counts[step] / (double)CHAIN_SAMPLES;
      totalCounts[step] += counts[step];
    }
  }
  CHECK(imbalanced == 0, "anti chains visited %d imbalanced states",
        imbalanced);
  CHECK(saw02, "anti chains did not observe 0<->2 transition");
  CHECK(saw24, "anti chains did not observe 2<->4 transition");
  for (chain = 0; chain < ANTI_STATE_COUNT; chain++) {
    double mean = 0.0;
    double variance = 0.0;
    double standardError;
    double tolerance;
    int sample;
    for (sample = 0; sample < CHAIN_COUNT; sample++) {
      mean += frequency[sample][chain];
    }
    mean /= CHAIN_COUNT;
    for (sample = 0; sample < CHAIN_COUNT; sample++) {
      const double delta = frequency[sample][chain] - mean;
      variance += delta * delta;
    }
    variance /= CHAIN_COUNT - 1;
    standardError = sqrt(variance / CHAIN_COUNT);
    tolerance = fmax(6.0 * standardError, 5.0e-4);
    CHECK(isfinite(mean) && isfinite(tolerance) &&
              fabs(mean - exact[chain]) <= tolerance,
          "anti chain state=%u mean=%.17g exact=%.17g SE=%.17g tol=%.17g",
          antiStates[chain], mean, exact[chain], standardError, tolerance);
    CHECK(totalCounts[chain] > 0, "anti chain state=%u never observed",
          antiStates[chain]);
  }
}


/* Three sites: m and L-m both reach 1 and 2, so a biased choice among
 * several empty or occupied sites changes the stationary distribution. */
#define THREE_SITES 3
#define THREE_ORBITALS 6

static double complex det3(const double complex F[THREE_SITES][THREE_SITES],
                           const int *rows, const int *cols, const int n) {
  if (n == 0) return 1.0;
  if (n == 1) return F[rows[0]][cols[0]];
  if (n == 2) {
    return F[rows[0]][cols[0]] * F[rows[1]][cols[1]] -
           F[rows[0]][cols[1]] * F[rows[1]][cols[0]];
  }
  return F[rows[0]][cols[0]] * (F[rows[1]][cols[1]] * F[rows[2]][cols[2]] -
                                F[rows[1]][cols[2]] * F[rows[2]][cols[1]]) -
         F[rows[0]][cols[1]] * (F[rows[1]][cols[0]] * F[rows[2]][cols[2]] -
                                F[rows[1]][cols[2]] * F[rows[2]][cols[0]]) +
         F[rows[0]][cols[2]] * (F[rows[1]][cols[0]] * F[rows[2]][cols[1]] -
                                F[rows[1]][cols[1]] * F[rows[2]][cols[0]]);
}

static void test_anti_chains_three_sites(void) {
  const double complex F[THREE_SITES][THREE_SITES] = {
      {0.71 + 0.12 * I, 0.34 - 0.23 * I, -0.29 + 0.17 * I},
      {-0.41 + 0.26 * I, 0.62 - 0.11 * I, 0.27 + 0.22 * I},
      {0.23 + 0.37 * I, -0.36 + 0.14 * I, 0.58 + 0.09 * I}};
  unsigned int states[64];
  double exact[64];
  double frequency[CHAIN_COUNT][64];
  double normalization = 0.0;
  int nstates = 0;
  int imbalanced = 0;
  int chain;
  int i;
  unsigned int mask;
  Nsite = THREE_SITES;
  Nsite2 = THREE_ORBITALS;
  NsizeMax = THREE_ORBITALS;
  free(SlaterElm);
  free(InvM);
  free(GutzwillerIdx);
  SlaterElm = calloc(THREE_ORBITALS * THREE_ORBITALS, sizeof(*SlaterElm));
  InvM = calloc(THREE_ORBITALS * THREE_ORBITALS, sizeof(*InvM));
  GutzwillerIdx = calloc(THREE_SITES, sizeof(*GutzwillerIdx));
  if (SlaterElm == NULL || InvM == NULL || GutzwillerIdx == NULL) {
    fprintf(stderr, "three-site allocation failed\n");
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < THREE_SITES; i++) {
    int j;
    for (j = 0; j < THREE_SITES; j++) {
      SlaterElm[(size_t)i * THREE_ORBITALS + (size_t)(THREE_SITES + j)] =
          F[i][j];
      SlaterElm[(size_t)(THREE_SITES + j) * THREE_ORBITALS + (size_t)i] =
          -F[i][j];
    }
  }
  for (mask = 0; mask < (1U << THREE_ORBITALS); mask++) {
    int up[THREE_SITES];
    int down[THREE_SITES];
    int nu = 0;
    int nd = 0;
    int doublon = 0;
    int site;
    double complex amplitude;
    for (site = 0; site < THREE_SITES; site++) {
      const int u = (int)((mask >> site) & 1U);
      const int d = (int)((mask >> (site + THREE_SITES)) & 1U);
      if (u) up[nu++] = site;
      if (d) down[nd++] = site;
      doublon += u & d;
    }
    if (nu != nd) continue;
    amplitude = det3(F, up, down, nu);
    states[nstates] = mask;
    exact[nstates] = exp(2.0 * creal(Proj[0]) * doublon) *
                     creal(amplitude * conj(amplitude));
    normalization += exact[nstates];
    nstates++;
  }
  CHECK(nstates == 20, "three-site balanced basis has %d states", nstates);
  for (i = 0; i < nstates; i++) exact[i] /= normalization;

  for (chain = 0; chain < CHAIN_COUNT; chain++) {
    int eleIdx[THREE_ORBITALS];
    int eleCfg[THREE_ORBITALS];
    int eleNum[THREE_ORBITALS];
    int eleProjCnt[1];
    int projCntNew[1];
    double complex pfMNew[1];
    double complex logIp;
    long long counts[64] = {0};
    int step;
    init_gen_rand((uint32_t)(5153 + 6007 * chain));
    Ncur = 2 * (chain % (THREE_SITES + 1));
    CHECK(makeInitialSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt, 0, 1,
                              MPI_COMM_SELF) == 0,
          "three-site initialization failed");
    logIp = CalculateLogIP_fcmp(PfM, 0, 1, MPI_COMM_SELF);
    for (step = 0; step < CHAIN_WARMUP + CHAIN_SAMPLES; step++) {
      int separation;
      for (separation = 0; separation < THREE_ORBITALS; separation++) {
        (void)GCMakeOneStep(eleIdx, eleCfg, eleNum, eleProjCnt, &logIp,
                            pfMNew, projCntNew, 0, 1, MPI_COMM_SELF);
      }
      if (step >= CHAIN_WARMUP) {
        unsigned int current = 0U;
        int rs;
        int index = -1;
        for (rs = 0; rs < THREE_ORBITALS; rs++) {
          if (eleNum[rs] != 0) current |= 1U << rs;
        }
        for (i = 0; i < nstates; i++) {
          if (states[i] == current) index = i;
        }
        if (index < 0) {
          imbalanced++;
        } else {
          counts[index]++;
        }
      }
    }
    for (i = 0; i < nstates; i++) {
      frequency[chain][i] = (double)counts[i] / (double)CHAIN_SAMPLES;
    }
  }
  CHECK(imbalanced == 0, "three-site chains left Sz=0 %d times", imbalanced);
  for (i = 0; i < nstates; i++) {
    double mean = 0.0;
    double variance = 0.0;
    double tolerance;
    int sample;
    for (sample = 0; sample < CHAIN_COUNT; sample++) {
      mean += frequency[sample][i];
    }
    mean /= CHAIN_COUNT;
    for (sample = 0; sample < CHAIN_COUNT; sample++) {
      const double delta = frequency[sample][i] - mean;
      variance += delta * delta;
    }
    variance /= CHAIN_COUNT - 1;
    tolerance = fmax(6.0 * sqrt(variance / CHAIN_COUNT), 5.0e-4);
    CHECK(isfinite(mean) && isfinite(tolerance) &&
              fabs(mean - exact[i]) <= tolerance,
          "three-site state=%u mean=%.17g exact=%.17g tol=%.17g", states[i],
          mean, exact[i], tolerance);
  }
}

/* Two ranks split NQPFull=1, so rank 1 owns no projection; it must still take
 * part in every collective and follow the same state. */
static void test_anti_collective(void) {
  int eleIdx[ORBITALS];
  int eleCfg[ORBITALS];
  int eleNum[ORBITALS];
  int eleProjCnt[1];
  int projCntNew[1];
  double complex pfMNew[1];
  double complex logIp;
  int qpStart;
  int qpEnd;
  int rank;
  int size;
  int step;
  int accepted[3] = {0, 0, 0};
  int imbalanced = 0;
  int mismatched = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  SplitLoop(&qpStart, &qpEnd, NQPFull, rank, size);
  init_gen_rand(80713U);
  Ncur = 2;
  CHECK(makeInitialSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt, qpStart, qpEnd,
                            MPI_COMM_WORLD) == 0,
        "collective initialization failed");
  logIp = CalculateLogIP_fcmp(PfM, qpStart, qpEnd, MPI_COMM_WORLD);
  for (step = 0; step < 4000; step++) {
    const int ncurBefore = Ncur;
    const unsigned int maskBefore = current_mask(eleNum);
    unsigned int masks[2];
    unsigned int localMask;
    if (GCMakeOneStep(eleIdx, eleCfg, eleNum, eleProjCnt, &logIp, pfMNew,
                      projCntNew, qpStart, qpEnd, MPI_COMM_WORLD)) {
      if (Ncur > ncurBefore) {
        accepted[1]++;
      } else if (Ncur < ncurBefore) {
        accepted[2]++;
      } else if (current_mask(eleNum) != maskBefore) {
        accepted[0]++;
      }
    }
    localMask = current_mask(eleNum);
    if (!spin_balanced(localMask)) imbalanced++;
    MPI_Allreduce(&localMask, &masks[0], 1, MPI_UNSIGNED, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&localMask, &masks[1], 1, MPI_UNSIGNED, MPI_MAX,
                  MPI_COMM_WORLD);
    if (masks[0] != masks[1]) mismatched++;
  }
  CHECK(size == 2, "collective mode expects two ranks (got %d)", size);
  CHECK(rank == 0 ? qpEnd - qpStart == 1 : qpEnd == qpStart,
        "rank %d projection range [%d,%d)", rank, qpStart, qpEnd);
  CHECK(imbalanced == 0, "rank %d visited %d imbalanced states", rank,
        imbalanced);
  CHECK(mismatched == 0, "ranks disagreed on %d states", mismatched);
  CHECK(accepted[0] > 0 && accepted[1] > 0 && accepted[2] > 0,
        "rank %d accepted hop/add/remove = %d/%d/%d", rank, accepted[0],
        accepted[1], accepted[2]);
}

/* Fatal paths for the death-test wrapper; each must abort with a message. */
static void anti_death_initialization(void) {
  int eleIdx[ORBITALS];
  int eleCfg[ORBITALS];
  int eleNum[ORBITALS];
  int eleProjCnt[1];
  memset(SlaterElm, 0, ORBITALS * ORBITALS * sizeof(*SlaterElm));
  Ncur = 2;
  (void)makeInitialSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt, 0, 1,
                            MPI_COMM_SELF);
}

static void anti_death_nonfinite_candidate(void) {
  int eleIdx[ORBITALS];
  int eleCfg[ORBITALS];
  int eleNum[ORBITALS];
  int eleProjCnt[1];
  int projCntNew[1];
  double complex pfMNew[1];
  double complex logIp = set_state(5U, eleIdx, eleCfg, eleNum, eleProjCnt);
  /* Poison the pair (up1, down1) only after the current state is built. */
  SlaterElm[(size_t)1 * ORBITALS + 3] = NAN;
  SlaterElm[(size_t)3 * ORBITALS + 1] = NAN;
  (void)GCAttemptMove(GC_MOVE_ADD, 1, 3, 0.5, eleIdx, eleCfg, eleNum,
                      eleProjCnt, &logIp, pfMNew, projCntNew, 0, 1,
                      MPI_COMM_SELF);
}

static void anti_death_nonfinite_old_log(void) {
  int eleIdx[ORBITALS];
  int eleCfg[ORBITALS];
  int eleNum[ORBITALS];
  int eleProjCnt[1];
  int projCntNew[1];
  double complex pfMNew[1];
  double complex logIp = set_state(5U, eleIdx, eleCfg, eleNum, eleProjCnt);
  logIp = INFINITY + 0.0 * I;
  (void)GCAttemptMove(GC_MOVE_ADD, 1, 3, 0.5, eleIdx, eleCfg, eleNum,
                      eleProjCnt, &logIp, pfMNew, projCntNew, 0, 1,
                      MPI_COMM_SELF);
}

static void anti_death_burn_state(void) {
  int eleIdx[ORBITALS];
  int eleCfg[ORBITALS];
  int eleNum[ORBITALS];
  int eleProjCnt[1];
  (void)set_state(3U, eleIdx, eleCfg, eleNum, eleProjCnt);
  BurnEleIdx = malloc(ORBITALS * sizeof(*BurnEleIdx));
  BurnEleCfg = malloc(ORBITALS * sizeof(*BurnEleCfg));
  BurnEleNum = malloc(ORBITALS * sizeof(*BurnEleNum));
  BurnEleProjCnt = malloc(sizeof(*BurnEleProjCnt));
  TmpEleIdx = malloc(ORBITALS * sizeof(*TmpEleIdx));
  TmpEleCfg = malloc(ORBITALS * sizeof(*TmpEleCfg));
  TmpEleNum = malloc(ORBITALS * sizeof(*TmpEleNum));
  TmpEleProjCnt = malloc(sizeof(*TmpEleProjCnt));
  if (BurnEleIdx == NULL || BurnEleCfg == NULL || BurnEleNum == NULL ||
      BurnEleProjCnt == NULL || TmpEleIdx == NULL || TmpEleCfg == NULL ||
      TmpEleNum == NULL || TmpEleProjCnt == NULL) {
    fprintf(stderr, "burn death allocation failed\n");
    exit(EXIT_FAILURE);
  }
  /* mask 3 = up0, up1: a spin-imbalanced stored chain. */
  copyToBurnSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt);
  BurnFlag = 1;
  NVMCWarmUp = 0;
  NVMCSample = 1;
  NVMCInterval = 1;
  VMCMakeSampleGC(MPI_COMM_SELF);
}

static int run_antiparallel_death(const char *mode) {
  enter_antiparallel_mode();
  init_gen_rand(60013U);
  if (strcmp(mode, "--antiparallel-death-init") == 0) {
    anti_death_initialization();
  } else if (strcmp(mode, "--antiparallel-death-nonfinite") == 0) {
    anti_death_nonfinite_candidate();
  } else if (strcmp(mode, "--antiparallel-death-oldlog") == 0) {
    anti_death_nonfinite_old_log();
  } else if (strcmp(mode, "--antiparallel-death-burn") == 0) {
    anti_death_burn_state();
  } else {
    fprintf(stderr, "unknown death mode %s\n", mode);
    return 2;
  }
  fprintf(stderr, "death mode %s returned without aborting\n", mode);
  return 3;
}

static void cleanup_globals(void) {
  FreeWorkSpaceAll();
  free(SlaterElm);
  free(InvM);
  free(PfM);
  free(QPFullWeight);
  free(Proj);
  free(GutzwillerIdx);
  free(BurnEleIdx);
  free(BurnEleCfg);
  free(BurnEleNum);
  free(BurnEleProjCnt);
  free(EleIdx);
  free(EleCfg);
  free(EleNum);
  free(EleProjCnt);
  free(EleNumSample);
  free(logSqPfFullSlater);
}

int main(int argc, char **argv) {
  int collectiveOnly = 0;
  const char *mode = argc == 2 ? argv[1] : "";
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#endif
  collectiveOnly = strcmp(mode, "--collective-only") == 0;
  setup_globals();
  if (strncmp(mode, "--antiparallel-death", 20) == 0) {
    return run_antiparallel_death(mode);
  }
  if (strcmp(mode, "--antiparallel") == 0) {
    enter_antiparallel_mode();
    test_anti_guards();
    test_anti_initialization();
    test_anti_transactions();
    test_anti_chains();
    test_anti_chains_three_sites();
  } else if (strcmp(mode, "--antiparallel-collective") == 0) {
    enter_antiparallel_mode();
    test_anti_collective();
  } else {
    CHECK(!GCAntiEnabled(FlagGrandCanonical, iFlgOrbitalGeneral),
          "general fixture entered the anti-parallel mode");
    test_collective_rebuild_status();
    if (!collectiveOnly) {
      test_selector();
      test_exact_detailed_balance();
      test_transactions();
      test_signed_matrix();
      test_burn_and_save();
      test_production_chains();
    }
  }
  cleanup_globals();
#ifdef _mpi_use
  MPI_Finalize();
#endif
  if (failures != 0) {
    fprintf(stderr, "GCSampler_Unit: %d failure(s)\n", failures);
    return EXIT_FAILURE;
  }
  printf("GCSampler_Unit: PASS\n");
  return EXIT_SUCCESS;
}
