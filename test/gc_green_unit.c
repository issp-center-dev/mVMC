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
extern int omp_get_num_threads(void);

#include "global.h"
#include "blas_externs.h"

#include "../src/mVMC/gc_size.c"
#include "../src/mVMC/workspace.c"
#include "../src/mVMC/projection.c"
#include "../src/mVMC/gc_config.c"
#include "../src/mVMC/matrix_gc.c"
#include "../src/mVMC/pfupdate_gc.c"

double complex CalculateIP_fcmp(double complex *const pfM, const int qpStart,
                                const int qpEnd, MPI_Comm comm) {
  double complex result = 0.0;
  int qpidx;
  (void)comm;
  for (qpidx = qpStart; qpidx < qpEnd; qpidx++) {
    result += QPFullWeight[qpidx] * pfM[qpidx - qpStart];
  }
  return result;
}

#include "../src/mVMC/locgrn_gc.c"

void StartTimer(int timer) { (void)timer; }
void StopTimer(int timer) { (void)timer; }

#include "../src/mVMC/calham_gc.c"
#include "../src/mVMC/calgrn_gc.c"

#define ORBITALS 8
#define ACTIVE 4
#define QPS 2
#define MAX_N 3
#define SMALL_COUNT 36

static int failures = 0;
static const int baseEleIdx[ACTIVE] = {0, 2, 4, 7};
static const uint64_t intGuard = UINT64_C(0x6543cdef1234abcd);
static const double complex complexGuard = 71.25 - 39.5 * I;
static const double doubleGuard = -913.125;
static int gutzStorage[4] = {0, 0, 0, 0};
static int anomalousRowsStorage[4][5] = {
    {1, 1, 0, 2, 1},
    {0, 2, 1, 1, 0},
    {0, 0, 0, 3, 1},
    {1, 3, 1, 0, 0},
};
static int *anomalousRows[4] = {
    anomalousRowsStorage[0], anomalousRowsStorage[1],
    anomalousRowsStorage[2], anomalousRowsStorage[3],
};
static double complex anomalousParameter[4] = {
    0.31 + 0.17 * I, 0.31 - 0.17 * I,
   -0.23 + 0.11 * I, -0.23 - 0.11 * I,
};

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "GCGreen_Unit FAIL: ");                                 \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, "\n");                                                  \
      failures++;                                                               \
    }                                                                           \
  } while (0)

typedef struct {
  GCGreenScratch scratch;
  int *projRaw;
  int *rsiRaw;
  int *rsjRaw;
  int *msaRaw;
  int *pfIRaw;
  double complex *pfMRaw;
  double complex *vecRaw;
  double complex *wRaw;
  double complex *smallRaw;
  double complex *pfWorkRaw;
  double *pfRRaw;
} GuardedScratch;

static int *guarded_int(const int count, int **raw) {
  *raw = malloc((size_t)(count + 4) * sizeof(**raw));
  memcpy(*raw, &intGuard, sizeof(intGuard));
  memcpy(*raw + count + 2, &intGuard, sizeof(intGuard));
  return *raw + 2;
}

static double complex *guarded_complex(const int count,
                                       double complex **raw) {
  *raw = malloc((size_t)(count + 2) * sizeof(**raw));
  (*raw)[0] = complexGuard;
  (*raw)[count + 1] = complexGuard;
  return *raw + 1;
}

static double *guarded_double(const int count, double **raw) {
  *raw = malloc((size_t)(count + 2) * sizeof(**raw));
  (*raw)[0] = doubleGuard;
  (*raw)[count + 1] = doubleGuard;
  return *raw + 1;
}

static void initialize_scratch(GuardedScratch *guarded) {
  memset(guarded, 0, sizeof(*guarded));
  guarded->scratch.maxN = MAX_N;
  guarded->scratch.pfLWork = SMALL_COUNT;
  guarded->scratch.projCntNew = guarded_int(NProj, &guarded->projRaw);
  guarded->scratch.rsi = guarded_int(MAX_N, &guarded->rsiRaw);
  guarded->scratch.rsj = guarded_int(MAX_N, &guarded->rsjRaw);
  guarded->scratch.msa = guarded_int(MAX_N, &guarded->msaRaw);
  guarded->scratch.pfIWork = guarded_int(2 * MAX_N, &guarded->pfIRaw);
  guarded->scratch.pfMNew = guarded_complex(QPS, &guarded->pfMRaw);
  guarded->scratch.vec = guarded_complex(MAX_N * ORBITALS, &guarded->vecRaw);
  guarded->scratch.w = guarded_complex(ORBITALS, &guarded->wRaw);
  guarded->scratch.smallMat =
      guarded_complex(SMALL_COUNT, &guarded->smallRaw);
  guarded->scratch.pfWork =
      guarded_complex(SMALL_COUNT, &guarded->pfWorkRaw);
  guarded->scratch.pfRWork = guarded_double(SMALL_COUNT, &guarded->pfRRaw);
}

static void check_int_guard(const int *raw, const int count,
                            const char *label) {
  uint64_t before;
  uint64_t after;
  memcpy(&before, raw, sizeof(before));
  memcpy(&after, raw + count + 2, sizeof(after));
  CHECK(before == intGuard && after == intGuard, "%s guard changed", label);
}

static void check_scratch_guards(const GuardedScratch *guarded) {
  check_int_guard(guarded->projRaw, NProj, "projection");
  check_int_guard(guarded->rsiRaw, MAX_N, "rsi");
  check_int_guard(guarded->rsjRaw, MAX_N, "rsj");
  check_int_guard(guarded->msaRaw, MAX_N, "msa");
  check_int_guard(guarded->pfIRaw, 2 * MAX_N, "pfIWork");
#define CHECK_COMPLEX_GUARD(raw, count, label)                                  \
  CHECK((raw)[0] == complexGuard && (raw)[(count) + 1] == complexGuard,         \
        "%s guard changed", (label))
  CHECK_COMPLEX_GUARD(guarded->pfMRaw, QPS, "pfMNew");
  CHECK_COMPLEX_GUARD(guarded->vecRaw, MAX_N * ORBITALS, "vec");
  CHECK_COMPLEX_GUARD(guarded->wRaw, ORBITALS, "w");
  CHECK_COMPLEX_GUARD(guarded->smallRaw, SMALL_COUNT, "smallMat");
  CHECK_COMPLEX_GUARD(guarded->pfWorkRaw, SMALL_COUNT, "pfWork");
#undef CHECK_COMPLEX_GUARD
  CHECK(guarded->pfRRaw[0] == doubleGuard &&
            guarded->pfRRaw[SMALL_COUNT + 1] == doubleGuard,
        "pfRWork guard changed");
}

static void free_scratch(GuardedScratch *guarded) {
  free(guarded->pfRRaw);
  free(guarded->pfWorkRaw);
  free(guarded->smallRaw);
  free(guarded->wRaw);
  free(guarded->vecRaw);
  free(guarded->pfMRaw);
  free(guarded->pfIRaw);
  free(guarded->msaRaw);
  free(guarded->rsjRaw);
  free(guarded->rsiRaw);
  free(guarded->projRaw);
}

static double complex pfaffian_recursive(const double complex *matrix,
                                          const int n) {
  double complex result = 0.0;
  int column;
  if (n == 0) return 1.0;
  if (n == 2) return matrix[1];
  for (column = 1; column < n; column++) {
    double complex minor[ORBITALS * ORBITALS];
    int sourceRow;
    int minorRow = 0;
    for (sourceRow = 1; sourceRow < n; sourceRow++) {
      int sourceColumn;
      int minorColumn = 0;
      if (sourceRow == column) continue;
      for (sourceColumn = 1; sourceColumn < n; sourceColumn++) {
        if (sourceColumn == column) continue;
        minor[(size_t)minorRow * (size_t)(n - 2) +
              (size_t)minorColumn++] =
            matrix[(size_t)sourceRow * (size_t)n +
                   (size_t)sourceColumn];
      }
      minorRow++;
    }
    result += ((column % 2) != 0 ? 1.0 : -1.0) * matrix[column] *
              pfaffian_recursive(minor, n - 2);
  }
  return result;
}

static double complex overlap_for_sorted_state(const int *occupied,
                                                const int count) {
  double complex overlap = 0.0;
  int qpidx;
  for (qpidx = 0; qpidx < QPS; qpidx++) {
    double complex matrix[ORBITALS * ORBITALS];
    int row;
    for (row = 0; row < count; row++) {
      int column;
      for (column = 0; column < count; column++) {
        matrix[(size_t)row * (size_t)count + (size_t)column] =
            -SlaterElm[(size_t)qpidx * ORBITALS * ORBITALS +
                       (size_t)occupied[row] * ORBITALS +
                       (size_t)occupied[column]];
      }
    }
    overlap += QPFullWeight[qpidx] * pfaffian_recursive(matrix, count);
  }
  return overlap;
}

static int parity_below(const unsigned int occupancy, const int orbital) {
  unsigned int lower = occupancy & ((1U << orbital) - 1U);
  int parity = 0;
  while (lower != 0U) {
    parity ^= (int)(lower & 1U);
    lower >>= 1;
  }
  return parity;
}

static int apply_annihilation(unsigned int *occupancy, const int orbital,
                              int *sign) {
  if (((*occupancy >> orbital) & 1U) == 0U) return 0;
  if (parity_below(*occupancy, orbital) != 0) *sign = -*sign;
  *occupancy &= ~(1U << orbital);
  return 1;
}

static int apply_creation(unsigned int *occupancy, const int orbital,
                          int *sign) {
  if (((*occupancy >> orbital) & 1U) != 0U) return 0;
  if (parity_below(*occupancy, orbital) != 0) *sign = -*sign;
  *occupancy |= 1U << orbital;
  return 1;
}

static void build_state(const int *baseIdx, const int ncur, int *eleCfg,
                        int *eleNum, int *projCnt) {
  int k;
  for (k = 0; k < ORBITALS; k++) {
    eleCfg[k] = -1;
    eleNum[k] = 0;
  }
  for (k = 0; k < ncur; k++) {
    eleCfg[baseIdx[k]] = k;
    eleNum[baseIdx[k]] = 1;
  }
  MakeProjCnt(projCnt, eleNum);
}

static double complex brute_anomalous_param(
    const int type, const int rs1, const int rs2, const int *baseIdx,
    const int ncur, const double complex baseSortedOverlap,
    const int *baseProjCnt) {
  unsigned int occupancy = 0U;
  int finalEleNum[ORBITALS] = {0};
  int finalEleIdx[ORBITALS];
  int finalProjCnt[1];
  int sign = 1;
  int k;
  int count = 0;
  for (k = 0; k < ncur; k++) occupancy |= 1U << baseIdx[k];
  if (type == 1) {
    if (!apply_creation(&occupancy, rs2, &sign) ||
        !apply_creation(&occupancy, rs1, &sign)) return 0.0;
  } else {
    if (!apply_annihilation(&occupancy, rs2, &sign) ||
        !apply_annihilation(&occupancy, rs1, &sign)) return 0.0;
  }
  for (k = 0; k < ORBITALS; k++) {
    if (((occupancy >> k) & 1U) != 0U) {
      finalEleIdx[count++] = k;
      finalEleNum[k] = 1;
    }
  }
  MakeProjCnt(finalProjCnt, finalEleNum);
  return (double)sign * conj(ProjRatio(finalProjCnt, baseProjCnt) *
                             overlap_for_sorted_state(finalEleIdx, count) /
                             baseSortedOverlap);
}

static void permute_base(const int *sorted, int *out, const int ncur,
                         const int kind) {
  int k;
  for (k = 0; k < ncur; k++) out[k] = sorted[k];
  if (ncur < 2) return;
  if (kind == 1) {
    for (k = 0; k < ncur; k++) out[k] = sorted[ncur - 1 - k];
  } else if (kind == 2) {
    for (k = 0; k < ncur; k++) out[k] = sorted[(k + 1) % ncur];
  } else if (kind == 3) {
    out[0] = sorted[1];
    out[1] = sorted[0];
  }
}

static void test_anomalous_green_exhaustive(void) {
  static const int bases[][ORBITALS] = {
      {0},
      {0, 5},
      {0, 2, 4, 7},
      {0, 1, 2, 3, 4, 5},
      {0, 1, 2, 3, 4, 5, 6, 7},
  };
  static const int ncurs[] = {0, 2, 4, 6, 8};
  int caseCount = 0;
  int nonzeroAdd = 0;
  int nonzeroRemove = 0;
  int nonzeroRemovePositiveSign = 0;
  int nonzeroRemoveNegativeSign = 0;
  int vacuumNonzeroAdd = 0;
  int vacuumNonzeroRemove = 0;
  int fullNonzeroAdd = 0;
  int fullNonzeroRemove = 0;
  int b;
  for (b = 0; b < 5; b++) {
    const int ncur = ncurs[b];
    int kind;
    for (kind = 0; kind < 4; kind++) {
      int baseIdx[ORBITALS] = {0};
      int eleCfgL[ORBITALS];
      int eleNumL[ORBITALS];
      int projCnt[1];
      int rs1;
      int rs2;
      double complex ipSorted;
      double complex ipPermuted;
      GuardedScratch guarded;
      if (ncur < 2 && kind > 0) continue;
      permute_base(bases[b], baseIdx, ncur, kind);
      build_state(baseIdx, ncur, eleCfgL, eleNumL, projCnt);
      CHECK(CalculateMAllGC_fcmp(ncur, baseIdx, 0, NQPFull) == GC_MALL_OK,
            "rebuild failed ncur=%d kind=%d", ncur, kind);
      ipPermuted = CalculateIP_fcmp(PfM, 0, NQPFull, MPI_COMM_SELF);
      ipSorted = overlap_for_sorted_state(bases[b], ncur);
      initialize_scratch(&guarded);
      for (rs1 = 0; rs1 < ORBITALS; rs1++) {
        for (rs2 = 0; rs2 < ORBITALS; rs2++) {
          double complex expectedAdd;
          double complex actualAdd;
          double complex expectedRemove;
          double complex actualRemove;
          if (rs1 == rs2) continue;
          expectedAdd = brute_anomalous_param(1, rs1, rs2, bases[b],
                                              ncur, ipSorted, projCnt);
          actualAdd = GreenFuncPairAddGC(rs1, rs2, ipPermuted, ncur,
                                         baseIdx, eleCfgL, eleNumL,
                                         projCnt, &guarded.scratch);
          CHECK(cabs(actualAdd - expectedAdd) <
                    1.0e-9 * (1.0 + cabs(expectedAdd)),
                "PairAdd ncur=%d kind=%d rs=(%d,%d)",
                ncur, kind, rs1, rs2);
          if (cabs(expectedAdd) > 1.0e-12) {
            nonzeroAdd++;
            if (ncur == 0) vacuumNonzeroAdd++;
            if (ncur == ORBITALS) fullNonzeroAdd++;
          }
          expectedRemove = brute_anomalous_param(0, rs1, rs2, bases[b],
                                                 ncur, ipSorted, projCnt);
          actualRemove = GreenFuncPairRemoveGC(rs1, rs2, ipPermuted, ncur,
                                               baseIdx, eleCfgL, eleNumL,
                                               projCnt, &guarded.scratch);
          CHECK(cabs(actualRemove - expectedRemove) <
                    1.0e-9 * (1.0 + cabs(expectedRemove)),
                "PairRemove ncur=%d kind=%d rs=(%d,%d)",
                ncur, kind, rs1, rs2);
          if (cabs(expectedRemove) > 1.0e-12) {
            nonzeroRemove++;
            if (eleCfgL[rs1] < eleCfgL[rs2]) {
              nonzeroRemovePositiveSign++;
            } else {
              nonzeroRemoveNegativeSign++;
            }
            if (ncur == 0) vacuumNonzeroRemove++;
            if (ncur == ORBITALS) fullNonzeroRemove++;
          }
          caseCount += 2;
          CHECK(eleNumL[rs1] == (eleCfgL[rs1] >= 0) &&
                    eleNumL[rs2] == (eleCfgL[rs2] >= 0),
                "eleNum not restored ncur=%d rs=(%d,%d)",
                ncur, rs1, rs2);
        }
      }
      check_scratch_guards(&guarded);
      free_scratch(&guarded);
    }
  }
  printf("anomalous oracle cases=%d nonzeroAdd=%d nonzeroRemove=%d "
         "removeSignPlus=%d removeSignMinus=%d vacuumAdd=%d fullRemove=%d\n",
         caseCount, nonzeroAdd, nonzeroRemove,
         nonzeroRemovePositiveSign, nonzeroRemoveNegativeSign,
         vacuumNonzeroAdd, fullNonzeroRemove);
  CHECK(caseCount == 1904, "anomalous oracle case inventory changed: %d",
        caseCount);
  CHECK(nonzeroAdd > 100 && nonzeroRemove > 100 &&
            nonzeroRemovePositiveSign > 0 &&
            nonzeroRemoveNegativeSign > 0,
        "anomalous oracle lacks add/remove/sign coverage");
  CHECK(vacuumNonzeroAdd > 0 && vacuumNonzeroRemove == 0 &&
            fullNonzeroAdd == 0 && fullNonzeroRemove > 0,
        "anomalous vacuum/full boundary coverage changed");
}

static void test_anomalous_omp_stress(void) {
  static const int stressBase[ORBITALS] = {0, 2, 4, 7};
  int baseIdx[ORBITALS] = {0};
  int eleCfgL[ORBITALS];
  int eleNumL[ORBITALS];
  int projCnt[1];
  double complex ip;
  double complex serialAdd;
  double complex serialRemove;
  const char *requiredText = getenv("MVMC_GC_REQUIRE_OMP_THREADS");
  const int requiredThreads = requiredText == NULL ? 0 : atoi(requiredText);
  unsigned int teamMask = 0U;
  unsigned int evaluationMask = 0U;
  int observedTeamSize = 0;
  int failuresBefore = failures;
  GuardedScratch serialGuarded;
  memcpy(baseIdx, stressBase, sizeof(baseIdx));
  build_state(baseIdx, 4, eleCfgL, eleNumL, projCnt);
  CHECK(CalculateMAllGC_fcmp(4, baseIdx, 0, NQPFull) == GC_MALL_OK,
        "OMP stress rebuild failed");
  ip = CalculateIP_fcmp(PfM, 0, NQPFull, MPI_COMM_SELF);
  initialize_scratch(&serialGuarded);
  serialAdd = GreenFuncPairAddGC(1, 6, ip, 4, baseIdx, eleCfgL, eleNumL,
                                 projCnt, &serialGuarded.scratch);
  serialRemove = GreenFuncPairRemoveGC(0, 7, ip, 4, baseIdx, eleCfgL,
                                       eleNumL, projCnt,
                                       &serialGuarded.scratch);
  CHECK(cabs(serialAdd) > 1.0e-12 && cabs(serialRemove) > 1.0e-12,
        "OMP stress pair fixture is vacuous");
  check_scratch_guards(&serialGuarded);
  free_scratch(&serialGuarded);
#pragma omp parallel shared(observedTeamSize, teamMask, evaluationMask)
  {
    GuardedScratch guarded;
    const int tid = omp_get_thread_num();
    int myEleNum[ORBITALS];
    int evaluated = 0;
    int iter;
#pragma omp critical(gc_anomalous_participation)
    {
      observedTeamSize = omp_get_num_threads();
      if (tid < (int)(8 * sizeof(teamMask))) teamMask |= 1U << tid;
    }
    initialize_scratch(&guarded);
    memcpy(myEleNum, eleNumL, sizeof(myEleNum));
#pragma omp for
    for (iter = 0; iter < 4000; iter++) {
      const double complex add = GreenFuncPairAddGC(
          1, 6, ip, 4, baseIdx, eleCfgL, myEleNum, projCnt,
          &guarded.scratch);
      const double complex remove = GreenFuncPairRemoveGC(
          0, 7, ip, 4, baseIdx, eleCfgL, myEleNum, projCnt,
          &guarded.scratch);
      evaluated = 1;
      if (cabs(add - serialAdd) > 1.0e-12 ||
          cabs(remove - serialRemove) > 1.0e-12) {
#pragma omp critical(gc_anomalous_failure)
        CHECK(0, "OMP stress mismatch iter=%d", iter);
      }
    }
#pragma omp critical(gc_anomalous_participation)
    {
      if (evaluated && tid < (int)(8 * sizeof(evaluationMask))) {
        evaluationMask |= 1U << tid;
      }
      check_scratch_guards(&guarded);
    }
    free_scratch(&guarded);
  }
  CHECK(failures == failuresBefore, "OMP stress failed");
  if (requiredThreads > 0) {
    const unsigned int requiredMask = (1U << requiredThreads) - 1U;
    CHECK(observedTeamSize == requiredThreads,
          "OMP team size got=%d expected=%d", observedTeamSize,
          requiredThreads);
    CHECK((teamMask & requiredMask) == requiredMask,
          "OMP team mask got=0x%x expected=0x%x", teamMask, requiredMask);
    CHECK((evaluationMask & requiredMask) == requiredMask,
          "OMP kernel-evaluation mask got=0x%x expected=0x%x",
          evaluationMask, requiredMask);
  }
}

static double complex brute_green(const int n, const int *rsi,
                                  const int *rsj,
                                  const double complex baseOverlap,
                                  const int *baseProjCnt) {
  unsigned int occupancy = 0U;
  int finalEleNum[ORBITALS] = {0};
  int finalEleIdx[ACTIVE];
  int finalProjCnt[1];
  int sign = 1;
  int k;
  int count = 0;
  for (k = 0; k < ACTIVE; k++) occupancy |= 1U << baseEleIdx[k];
  for (k = n - 1; k >= 0; k--) {
    if (!apply_annihilation(&occupancy, rsj[k], &sign) ||
        !apply_creation(&occupancy, rsi[k], &sign)) {
      return 0.0;
    }
  }
  for (k = 0; k < ORBITALS; k++) {
    if (((occupancy >> k) & 1U) != 0U) {
      finalEleIdx[count++] = k;
      finalEleNum[k] = 1;
    }
  }
  CHECK(count == ACTIVE, "operator changed particle count");
  MakeProjCnt(finalProjCnt, finalEleNum);
  return (double)sign * conj(ProjRatio(finalProjCnt, baseProjCnt) *
                             overlap_for_sorted_state(finalEleIdx, count) /
                             baseOverlap);
}

static void assert_state_unchanged(const int *eleIdx, const int *eleCfg,
                                   const int *eleNum, const int *projCnt,
                                   const double complex *invSnapshot,
                                   const double complex *pfSnapshot,
                                   const char *label) {
  int expectedCfg[ORBITALS];
  int expectedNum[ORBITALS] = {0};
  int expectedIdx[ORBITALS];
  int i;
  for (i = 0; i < ORBITALS; i++) {
    expectedCfg[i] = -1;
    expectedIdx[i] = -1;
  }
  for (i = 0; i < ACTIVE; i++) {
    expectedIdx[i] = baseEleIdx[i];
    expectedCfg[baseEleIdx[i]] = i;
    expectedNum[baseEleIdx[i]] = 1;
  }
  CHECK(memcmp(eleIdx, expectedIdx, sizeof(expectedIdx)) == 0,
        "%s changed eleIdx", label);
  CHECK(memcmp(eleCfg, expectedCfg, sizeof(expectedCfg)) == 0,
        "%s changed eleCfg", label);
  CHECK(memcmp(eleNum, expectedNum, sizeof(expectedNum)) == 0,
        "%s changed eleNum", label);
  {
    int expectedProj[1];
    MakeProjCnt(expectedProj, expectedNum);
    CHECK(memcmp(projCnt, expectedProj, sizeof(expectedProj)) == 0,
          "%s changed projection counts", label);
  }
  CHECK(memcmp(InvM, invSnapshot,
               QPS * ORBITALS * ORBITALS * sizeof(*InvM)) == 0,
        "%s changed shared InvM", label);
  CHECK(memcmp(PfM, pfSnapshot, QPS * sizeof(*PfM)) == 0,
        "%s changed shared PfM", label);
}

static void initialize_state(int *eleIdx, int *eleCfg, int *eleNum,
                             int *projCnt) {
  int i;
  for (i = 0; i < ORBITALS; i++) {
    eleIdx[i] = -1;
    eleCfg[i] = -1;
    eleNum[i] = 0;
  }
  for (i = 0; i < ACTIVE; i++) {
    eleIdx[i] = baseEleIdx[i];
    eleCfg[eleIdx[i]] = i;
    eleNum[eleIdx[i]] = 1;
  }
  MakeProjCnt(projCnt, eleNum);
}

static void test_green_kernels(void) {
  GuardedScratch guarded;
  int eleIdx[ORBITALS];
  int eleCfg[ORBITALS];
  int eleNum[ORBITALS];
  int projCnt[1];
  double complex invSnapshot[QPS * ORBITALS * ORBITALS];
  double complex pfSnapshot[QPS];
  double complex baseOverlap;
  int i;
  int j;
  int k;
  int l;
  initialize_scratch(&guarded);
  initialize_state(eleIdx, eleCfg, eleNum, projCnt);
  CHECK(CalculateMAllGC_fcmp(ACTIVE, eleIdx, 0, QPS) == GC_MALL_OK,
        "baseline rebuild failed");
  baseOverlap = CalculateIP_fcmp(PfM, 0, QPS, MPI_COMM_SELF);
  memcpy(invSnapshot, InvM, sizeof(invSnapshot));
  memcpy(pfSnapshot, PfM, sizeof(pfSnapshot));

  for (i = 0; i < ORBITALS; i++) {
    for (j = 0; j < ORBITALS; j++) {
      int create[1] = {i};
      int annihilate[1] = {j};
      const double complex actual = GreenFunc1GC(
          i, j, baseOverlap, ACTIVE, eleIdx, eleCfg, eleNum, projCnt,
          &guarded.scratch);
      const double complex expected =
          brute_green(1, create, annihilate, baseOverlap, projCnt);
      CHECK(cabs(actual - expected) < 3.0e-9,
            "Green1 (%d,%d) got=(%.17g,%.17g) expected=(%.17g,%.17g)",
            i, j, creal(actual), cimag(actual), creal(expected),
            cimag(expected));
      assert_state_unchanged(eleIdx, eleCfg, eleNum, projCnt, invSnapshot,
                             pfSnapshot, "Green1");
    }
  }

  for (i = 0; i < ORBITALS; i++) {
    for (j = 0; j < ORBITALS; j++) {
      for (k = 0; k < ORBITALS; k++) {
        for (l = 0; l < ORBITALS; l++) {
          int create[2] = {i, k};
          int annihilate[2] = {j, l};
          const double complex actual = GreenFunc2GC(
              i, j, k, l, baseOverlap, ACTIVE, eleIdx, eleCfg, eleNum,
              projCnt, &guarded.scratch);
          const double complex expected =
              brute_green(2, create, annihilate, baseOverlap, projCnt);
          CHECK(cabs(actual - expected) < 7.0e-9,
                "Green2 (%d,%d,%d,%d) got=(%.17g,%.17g) expected=(%.17g,%.17g)",
                i, j, k, l, creal(actual), cimag(actual), creal(expected),
                cimag(expected));
          assert_state_unchanged(eleIdx, eleCfg, eleNum, projCnt, invSnapshot,
                                 pfSnapshot, "Green2");
        }
      }
    }
  }

  {
    int create1[1] = {1};
    int annihilate1[1] = {0};
    int create2[2] = {1, 3};
    int annihilate2[2] = {0, 2};
    int create3[3] = {1, 3, 5};
    int annihilate3[3] = {0, 2, 4};
    int *creates[3] = {create1, create2, create3};
    int *annihilates[3] = {annihilate1, annihilate2, annihilate3};
    int n;
    for (n = 1; n <= 3; n++) {
      int rsi[MAX_N];
      int rsj[MAX_N];
      double complex actual;
      double complex expected;
      memcpy(rsi, creates[n - 1], (size_t)n * sizeof(*rsi));
      memcpy(rsj, annihilates[n - 1], (size_t)n * sizeof(*rsj));
      actual = GreenFuncNGC(n, rsi, rsj, baseOverlap, ACTIVE, eleIdx, eleCfg,
                           eleNum, projCnt, &guarded.scratch);
      expected = brute_green(n, creates[n - 1], annihilates[n - 1],
                             baseOverlap, projCnt);
      CHECK(cabs(actual - expected) < 7.0e-9,
            "GreenN n=%d got=(%.17g,%.17g) expected=(%.17g,%.17g)", n,
            creal(actual), cimag(actual), creal(expected), cimag(expected));
      assert_state_unchanged(eleIdx, eleCfg, eleNum, projCnt, invSnapshot,
                             pfSnapshot, "GreenN");
    }
  }
  check_scratch_guards(&guarded);
  free_scratch(&guarded);
}

static void test_hamiltonian_and_measurement(void) {
  int eleIdx[ORBITALS];
  int eleCfg[ORBITALS];
  int eleNum[ORBITALS];
  int projCnt[1];
  double complex invSnapshot[QPS * ORBITALS * ORBITALS];
  double complex pfSnapshot[QPS];
  double complex baseOverlap;
  double complex baseSortedOverlap;
  double complex expectedEnergy;
  double complex actualEnergy;
  double complex classContribution;
  double complex createContribution;
  double complex removeContribution;
  double complex anomalousPairContribution;
  double complex mutatedEnergy;
  int transferStorage[2][4] = {{1, 0, 0, 0}, {2, 1, 3, 1}};
  int *transferRows[2] = {transferStorage[0], transferStorage[1]};
  double complex transferParameter[2] = {0.41 - 0.13 * I,
                                         -0.22 + 0.09 * I};
  int interStorage[1][8] = {{1, 0, 0, 0, 2, 1, 3, 1}};
  int *interRows[1] = {interStorage[0]};
  double complex interParameter[1] = {0.17 + 0.08 * I};
  int nbodyN[1] = {3};
  int nbodyOffset[1] = {0};
  int nbodyStorage[3][4] = {
      {1, 0, 0, 0}, {3, 0, 2, 0}, {1, 1, 0, 1}};
  int *nbodyRows[3] = {nbodyStorage[0], nbodyStorage[1], nbodyStorage[2]};
  double complex nbodyParameter[1] = {-0.12 + 0.06 * I};
  int coulombIntraStorage[2] = {0, 1};
  double coulombIntraParameter[2] = {0.43, -0.19};
  int coulombInterStorage[1][2] = {{0, 2}};
  int *coulombInterRows[1] = {coulombInterStorage[0]};
  double coulombInterParameter[1] = {0.31};
  int hundStorage[1][2] = {{0, 2}};
  int *hundRows[1] = {hundStorage[0]};
  double hundParameter[1] = {0.23};
  int pairHoppingStorage[2][2] = {{0, 1}, {1, 0}};
  int *pairHoppingRows[2] = {
      pairHoppingStorage[0], pairHoppingStorage[1]};
  double pairHoppingParameter[2] = {0.29, 0.29};
  int exchangeStorage[1][2] = {{2, 3}};
  int *exchangeRows[1] = {exchangeStorage[0]};
  double exchangeParameter[1] = {-0.37};
  int cisStorage[2][4] = {{1, 0, 0, 0}, {2, 1, 3, 1}};
  int *cisRows[2] = {cisStorage[0], cisStorage[1]};
  int dcStorage[1][8] = {{1, 0, 0, 0, 2, 1, 3, 1}};
  int *dcRows[1] = {dcStorage[0]};
  int productStorage[1][2] = {{0, 1}};
  int *productRows[1] = {productStorage[0]};
  int nbodyGStorage[1] = {3};
  int nbodyGOffsetStorage[1] = {0};
  int nbodyGFactors[3][4] = {
      {1, 0, 0, 0}, {3, 0, 2, 0}, {1, 1, 0, 1}};
  int *nbodyGFactorRows[3] = {
      nbodyGFactors[0], nbodyGFactors[1], nbodyGFactors[2]};
  double complex local[2] = {0.0, 0.0};
  double complex phys1[2] = {0.0, 0.0};
  double complex phys2[1] = {0.0};
  double complex physDC[1] = {0.0};
  double complex physN[1] = {0.0};
  double complex physAnomalous[4] = {0.0, 0.0, 0.0, 0.0};
  const double measurementWeight = 0.37;
  int create1[1];
  int annihilate1[1];
  int create2[2];
  int annihilate2[2];
  int create3[3];
  int annihilate3[3];
  double complex expectedLocal0;
  double complex expectedLocal1;
  double complex expectedDC;
  double complex expectedN;
  int anomalousIdx;

  initialize_state(eleIdx, eleCfg, eleNum, projCnt);
  CHECK(CalculateMAllGC_fcmp(ACTIVE, eleIdx, 0, QPS) == GC_MALL_OK,
        "Hamiltonian baseline rebuild failed");
  baseOverlap = CalculateIP_fcmp(PfM, 0, QPS, MPI_COMM_SELF);
  baseSortedOverlap = overlap_for_sorted_state(baseEleIdx, ACTIVE);
  memcpy(invSnapshot, InvM, sizeof(invSnapshot));
  memcpy(pfSnapshot, PfM, sizeof(pfSnapshot));

  NCoulombIntra = 2;
  CoulombIntra = coulombIntraStorage;
  ParaCoulombIntra = coulombIntraParameter;
  NCoulombInter = 1;
  CoulombInter = coulombInterRows;
  ParaCoulombInter = coulombInterParameter;
  NHundCoupling = 1;
  HundCoupling = hundRows;
  ParaHundCoupling = hundParameter;
  NTransfer = 2;
  Transfer = transferRows;
  ParaTransfer = transferParameter;
  NPairHopping = 2;
  PairHopping = pairHoppingRows;
  ParaPairHopping = pairHoppingParameter;
  NExchangeCoupling = 1;
  ExchangeCoupling = exchangeRows;
  ParaExchangeCoupling = exchangeParameter;
  NInterAll = 1;
  InterAll = interRows;
  ParaInterAll = interParameter;
  NNBodyInterAll = 1;
  NBodyInterAllMaxN = 3;
  NBodyInterAllN = nbodyN;
  NBodyInterAllOffset = nbodyOffset;
  NBodyInterAllIdx = nbodyRows;
  ParaNBodyInterAll = nbodyParameter;
  NAnomalousTerm = 4;
  AnomalousTerm = anomalousRows;
  ParaAnomalousTerm = anomalousParameter;

  create1[0] = 1;
  annihilate1[0] = 0;
  classContribution = coulombIntraParameter[0] * eleNum[0] * eleNum[4] +
                      coulombIntraParameter[1] * eleNum[1] * eleNum[5];
  CHECK(cabs(classContribution) > 1.0e-6,
        "CoulombIntra fixture is vacuous");
  expectedEnergy = classContribution;
  classContribution = coulombInterParameter[0] *
                      (eleNum[0] + eleNum[4]) *
                      (eleNum[2] + eleNum[6]);
  CHECK(cabs(classContribution) > 1.0e-6,
        "CoulombInter fixture is vacuous");
  expectedEnergy += classContribution;
  classContribution = -hundParameter[0] *
                      (eleNum[0] * eleNum[2] + eleNum[4] * eleNum[6]);
  CHECK(cabs(classContribution) > 1.0e-6, "Hund fixture is vacuous");
  expectedEnergy += classContribution;
  classContribution = -transferParameter[0] *
                      brute_green(1, create1, annihilate1, baseOverlap,
                                  projCnt);
  create1[0] = 6;
  annihilate1[0] = 7;
  classContribution -= transferParameter[1] *
                       brute_green(1, create1, annihilate1, baseOverlap,
                                   projCnt);
  CHECK(cabs(classContribution) > 1.0e-6, "Transfer fixture is vacuous");
  expectedEnergy += classContribution;
  create2[0] = 0;
  annihilate2[0] = 1;
  create2[1] = 4;
  annihilate2[1] = 5;
  classContribution = pairHoppingParameter[0] *
                      brute_green(2, create2, annihilate2, baseOverlap,
                                  projCnt);
  create2[0] = 1;
  annihilate2[0] = 0;
  create2[1] = 5;
  annihilate2[1] = 4;
  classContribution += pairHoppingParameter[1] *
                       brute_green(2, create2, annihilate2, baseOverlap,
                                   projCnt);
  CHECK(cabs(classContribution) > 1.0e-6, "PairHop fixture is vacuous");
  expectedEnergy += classContribution;
  create2[0] = 2;
  annihilate2[0] = 3;
  create2[1] = 7;
  annihilate2[1] = 6;
  classContribution = exchangeParameter[0] *
                      brute_green(2, create2, annihilate2, baseOverlap,
                                  projCnt);
  create2[0] = 6;
  annihilate2[0] = 7;
  create2[1] = 3;
  annihilate2[1] = 2;
  classContribution += exchangeParameter[0] *
                       brute_green(2, create2, annihilate2, baseOverlap,
                                   projCnt);
  CHECK(cabs(classContribution) > 1.0e-6, "Exchange fixture is vacuous");
  expectedEnergy += classContribution;
  create2[0] = 1;
  annihilate2[0] = 0;
  create2[1] = 6;
  annihilate2[1] = 7;
  classContribution = interParameter[0] *
                      brute_green(2, create2, annihilate2, baseOverlap,
                                  projCnt);
  CHECK(cabs(classContribution) > 1.0e-6, "InterAll fixture is vacuous");
  expectedEnergy += classContribution;
  create3[0] = 1;
  annihilate3[0] = 0;
  create3[1] = 3;
  annihilate3[1] = 2;
  create3[2] = 5;
  annihilate3[2] = 4;
  classContribution = nbodyParameter[0] *
                      brute_green(3, create3, annihilate3, baseOverlap,
                                  projCnt);
  CHECK(cabs(classContribution) > 1.0e-6, "NBody fixture is vacuous");
  expectedEnergy += classContribution;
  createContribution =
      anomalousParameter[0] * brute_anomalous_param(
          1, 1, 6, baseEleIdx, ACTIVE, baseSortedOverlap, projCnt) +
      anomalousParameter[3] * brute_anomalous_param(
          1, 7, 0, baseEleIdx, ACTIVE, baseSortedOverlap, projCnt);
  removeContribution =
      anomalousParameter[1] * brute_anomalous_param(
          0, 6, 1, baseEleIdx, ACTIVE, baseSortedOverlap, projCnt) +
      anomalousParameter[2] * brute_anomalous_param(
          0, 0, 7, baseEleIdx, ACTIVE, baseSortedOverlap, projCnt);
  CHECK(cabs(createContribution) > 1.0e-6,
        "AnomalousTerm create fixture is vacuous");
  CHECK(cabs(removeContribution) > 1.0e-6,
        "AnomalousTerm remove fixture is vacuous");
  classContribution = createContribution + removeContribution;
  CHECK(cabs(classContribution) > 1.0e-6,
        "AnomalousTerm total fixture is vacuous");
  expectedEnergy += classContribution;
  actualEnergy = CalculateHamiltonianGC(baseOverlap, ACTIVE, eleIdx, eleCfg,
                                        eleNum, projCnt);
  CHECK(cabs(actualEnergy - expectedEnergy) < 2.0e-8,
        "Hamiltonian got=(%.17g,%.17g) expected=(%.17g,%.17g)",
        creal(actualEnergy), cimag(actualEnergy), creal(expectedEnergy),
        cimag(expectedEnergy));
  anomalousPairContribution =
      anomalousParameter[0] * brute_anomalous_param(
          1, 1, 6, baseEleIdx, ACTIVE, baseSortedOverlap, projCnt) +
      anomalousParameter[1] * brute_anomalous_param(
          0, 6, 1, baseEleIdx, ACTIVE, baseSortedOverlap, projCnt);
  anomalousParameter[0] = -anomalousParameter[0];
  anomalousParameter[1] = -anomalousParameter[1];
  mutatedEnergy = CalculateHamiltonianGC(baseOverlap, ACTIVE, eleIdx, eleCfg,
                                         eleNum, projCnt);
  CHECK(cabs((mutatedEnergy - actualEnergy) +
             2.0 * anomalousPairContribution) < 2.0e-8,
        "AnomalousTerm coefficient mutation mismatch");
  anomalousParameter[0] = -anomalousParameter[0];
  anomalousParameter[1] = -anomalousParameter[1];
  CHECK(CalculateSzGC(eleNum) == 0.0, "CalculateSzGC normalization");
  assert_state_unchanged(eleIdx, eleCfg, eleNum, projCnt, invSnapshot,
                         pfSnapshot, "Hamiltonian");

  NCisAjs = 2;
  CisAjsIdx = cisRows;
  NCisAjsCktAltDC = 1;
  CisAjsCktAltDCIdx = dcRows;
  NCisAjsCktAlt = 1;
  CisAjsCktAltIdx = productRows;
  NNBodyG = 1;
  NBodyGMaxN = 3;
  NBodyGN = nbodyGStorage;
  NBodyGOffset = nbodyGOffsetStorage;
  NBodyGIdx = nbodyGFactorRows;
  NTwist = 0;
  LocalCisAjs = local;
  PhysCisAjs = phys1;
  PhysCisAjsCktAlt = phys2;
  PhysCisAjsCktAltDC = physDC;
  PhysNBodyG = physN;
  NAnomalousG = 4;
  AnomalousG = anomalousRows;
  PhysAnomalousG = physAnomalous;
  create1[0] = 1;
  annihilate1[0] = 0;
  expectedLocal0 = brute_green(1, create1, annihilate1, baseOverlap, projCnt);
  create1[0] = 6;
  annihilate1[0] = 7;
  expectedLocal1 = brute_green(1, create1, annihilate1, baseOverlap, projCnt);
  create2[0] = 1;
  annihilate2[0] = 0;
  create2[1] = 6;
  annihilate2[1] = 7;
  expectedDC = brute_green(2, create2, annihilate2, baseOverlap, projCnt);
  create3[0] = 1;
  annihilate3[0] = 0;
  create3[1] = 3;
  annihilate3[1] = 2;
  create3[2] = 5;
  annihilate3[2] = 4;
  expectedN = brute_green(3, create3, annihilate3, baseOverlap, projCnt);
  CalculateGreenFuncGC(measurementWeight, baseOverlap, ACTIVE, eleIdx, eleCfg,
                       eleNum, projCnt);
  CHECK(cabs(local[0] - expectedLocal0) < 8.0e-9 &&
            cabs(local[1] - expectedLocal1) < 8.0e-9,
        "measurement LocalCisAjs mismatch");
  CHECK(cabs(phys1[0] - measurementWeight * expectedLocal0) < 8.0e-9 &&
            cabs(phys1[1] - measurementWeight * expectedLocal1) < 8.0e-9,
        "measurement PhysCisAjs mismatch");
  CHECK(cabs(phys2[0] - measurementWeight * expectedLocal0 *
                              conj(expectedLocal1)) < 8.0e-9,
        "measurement product mismatch");
  CHECK(cabs(physDC[0] - measurementWeight * expectedDC) < 8.0e-9,
        "measurement direct two-body mismatch");
  CHECK(cabs(physN[0] - measurementWeight * expectedN) < 8.0e-9,
        "measurement n-body mismatch");
  for (anomalousIdx = 0; anomalousIdx < NAnomalousG; anomalousIdx++) {
    const int rs1 = AnomalousG[anomalousIdx][1] +
                    AnomalousG[anomalousIdx][2] * Nsite;
    const int rs2 = AnomalousG[anomalousIdx][3] +
                    AnomalousG[anomalousIdx][4] * Nsite;
    const double complex expected = brute_anomalous_param(
        AnomalousG[anomalousIdx][0], rs1, rs2, baseEleIdx, ACTIVE,
        baseSortedOverlap, projCnt);
    CHECK(cabs(PhysAnomalousG[anomalousIdx] -
               measurementWeight * expected) < 8.0e-9,
          "AnomalousG measurement row=%d", anomalousIdx);
  }
  CHECK(cabs(physAnomalous[0]) > 1.0e-6 &&
            cabs(physAnomalous[2]) > 1.0e-6,
        "AnomalousG create/remove fixture is vacuous");
  assert_state_unchanged(eleIdx, eleCfg, eleNum, projCnt, invSnapshot,
                         pfSnapshot, "measurement");
}

static void initialize_fixture(void) {
  int qpidx;
  NThread = omp_get_max_threads();
  Nsite = 4;
  Nsite2 = ORBITALS;
  Nsize = ACTIVE;
  NsizeMax = ORBITALS;
  NQPFull = QPS;
  LapackLWork = 256;
  NProj = 1;
  NGutzwillerIdx = 1;
  NJastrowIdx = 0;
  NSpinJastrowIdx = 0;
  NDoublonHolon2siteIdx = 0;
  NDoublonHolon4siteIdx = 0;
  GutzwillerIdx = gutzStorage;
  Proj = malloc(sizeof(*Proj));
  SlaterElm = malloc(QPS * ORBITALS * ORBITALS * sizeof(*SlaterElm));
  InvM = malloc(QPS * ORBITALS * ORBITALS * sizeof(*InvM));
  PfM = malloc(QPS * sizeof(*PfM));
  QPFullWeight = malloc(QPS * sizeof(*QPFullWeight));
  CHECK(Proj != NULL && SlaterElm != NULL && InvM != NULL && PfM != NULL &&
            QPFullWeight != NULL,
        "fixture allocation failed");
  Proj[0] = 0.23;
  QPFullWeight[0] = 0.61 - 0.17 * I;
  QPFullWeight[1] = -0.29 + 0.37 * I;
  for (qpidx = 0; qpidx < QPS; qpidx++) {
    int row;
    double complex *slater =
        SlaterElm + (size_t)qpidx * ORBITALS * ORBITALS;
    memset(slater, 0, ORBITALS * ORBITALS * sizeof(*slater));
    for (row = 0; row < ORBITALS; row++) {
      int column;
      for (column = row + 1; column < ORBITALS; column++) {
        const double complex value =
            (0.41 + 0.13 * (row + 1) + 0.071 * (column + 1) +
             0.037 * qpidx * (row + column + 1)) +
            (0.19 - 0.043 * (row + 1) + 0.029 * (column + 1) +
             0.017 * qpidx) * I;
        slater[(size_t)row * ORBITALS + (size_t)column] = value;
        slater[(size_t)column * ORBITALS + (size_t)row] = -value;
      }
    }
  }
  initializeWorkSpaceAll();
}

static void free_fixture(void) {
  FreeWorkSpaceAll();
  free(QPFullWeight);
  free(PfM);
  free(InvM);
  free(SlaterElm);
  free(Proj);
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  initialize_fixture();
  test_green_kernels();
  test_anomalous_green_exhaustive();
  test_anomalous_omp_stress();
  test_hamiltonian_and_measurement();
  free_fixture();
#ifdef _mpi_use
  MPI_Finalize();
#endif
  if (failures != 0) {
    fprintf(stderr, "%d GC Green checks failed\n", failures);
    return EXIT_FAILURE;
  }
  printf("GC Green checks: PASS\n");
  return EXIT_SUCCESS;
}
