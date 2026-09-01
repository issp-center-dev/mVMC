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
  NThread = 1;
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
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argv;
#endif
  collectiveOnly = argc == 2 && strcmp(argv[1], "--collective-only") == 0;
  setup_globals();
  test_collective_rebuild_status();
  if (!collectiveOnly) {
    test_selector();
    test_exact_detailed_balance();
    test_transactions();
    test_burn_and_save();
    test_production_chains();
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
