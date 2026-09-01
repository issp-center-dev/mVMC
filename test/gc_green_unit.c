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
static const int baseEleIdx[ACTIVE] = {0, 2, 5, 7};
static const uint64_t intGuard = UINT64_C(0x6543cdef1234abcd);
static const double complex complexGuard = 71.25 - 39.5 * I;
static const double doubleGuard = -913.125;
static int gutzStorage[4] = {0, 0, 0, 0};

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
    int create3[3] = {1, 3, 4};
    int annihilate3[3] = {0, 2, 5};
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
  double complex expectedEnergy;
  double complex actualEnergy;
  int transferStorage[2][4] = {{1, 0, 0, 0}, {3, 1, 3, 0}};
  int *transferRows[2] = {transferStorage[0], transferStorage[1]};
  double complex transferParameter[2] = {0.41 - 0.13 * I,
                                         -0.22 + 0.09 * I};
  int interStorage[1][8] = {{1, 0, 0, 0, 3, 1, 1, 1}};
  int *interRows[1] = {interStorage[0]};
  double complex interParameter[1] = {0.17 + 0.08 * I};
  int nbodyN[1] = {3};
  int nbodyOffset[1] = {0};
  int nbodyStorage[3][4] = {
      {1, 0, 0, 0}, {3, 0, 2, 0}, {0, 1, 1, 1}};
  int *nbodyRows[3] = {nbodyStorage[0], nbodyStorage[1], nbodyStorage[2]};
  double complex nbodyParameter[1] = {-0.12 + 0.06 * I};
  int coulombInterStorage[1][2] = {{0, 1}};
  int *coulombInterRows[1] = {coulombInterStorage[0]};
  double coulombInterParameter[1] = {0.31};
  int cisStorage[2][4] = {{1, 0, 0, 0}, {3, 1, 3, 0}};
  int *cisRows[2] = {cisStorage[0], cisStorage[1]};
  int dcStorage[1][8] = {{1, 0, 0, 0, 3, 1, 1, 1}};
  int *dcRows[1] = {dcStorage[0]};
  int productStorage[1][2] = {{0, 1}};
  int *productRows[1] = {productStorage[0]};
  int nbodyGStorage[1] = {3};
  int nbodyGOffsetStorage[1] = {0};
  int nbodyGFactors[3][4] = {
      {1, 0, 0, 0}, {3, 0, 2, 0}, {0, 1, 1, 1}};
  int *nbodyGFactorRows[3] = {
      nbodyGFactors[0], nbodyGFactors[1], nbodyGFactors[2]};
  double complex local[2] = {0.0, 0.0};
  double complex phys1[2] = {0.0, 0.0};
  double complex phys2[1] = {0.0};
  double complex physDC[1] = {0.0};
  double complex physN[1] = {0.0};
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

  initialize_state(eleIdx, eleCfg, eleNum, projCnt);
  CHECK(CalculateMAllGC_fcmp(ACTIVE, eleIdx, 0, QPS) == GC_MALL_OK,
        "Hamiltonian baseline rebuild failed");
  baseOverlap = CalculateIP_fcmp(PfM, 0, QPS, MPI_COMM_SELF);
  memcpy(invSnapshot, InvM, sizeof(invSnapshot));
  memcpy(pfSnapshot, PfM, sizeof(pfSnapshot));

  NCoulombIntra = 0;
  NCoulombInter = 1;
  CoulombInter = coulombInterRows;
  ParaCoulombInter = coulombInterParameter;
  NHundCoupling = 0;
  NTransfer = 2;
  Transfer = transferRows;
  ParaTransfer = transferParameter;
  NPairHopping = 0;
  NExchangeCoupling = 0;
  NInterAll = 1;
  InterAll = interRows;
  ParaInterAll = interParameter;
  NNBodyInterAll = 1;
  NBodyInterAllMaxN = 3;
  NBodyInterAllN = nbodyN;
  NBodyInterAllOffset = nbodyOffset;
  NBodyInterAllIdx = nbodyRows;
  ParaNBodyInterAll = nbodyParameter;

  create1[0] = 1;
  annihilate1[0] = 0;
  expectedEnergy = coulombInterParameter[0];
  expectedEnergy -= transferParameter[0] *
                    brute_green(1, create1, annihilate1, baseOverlap, projCnt);
  create1[0] = 7;
  annihilate1[0] = 3;
  expectedEnergy -= transferParameter[1] *
                    brute_green(1, create1, annihilate1, baseOverlap, projCnt);
  create2[0] = 1;
  annihilate2[0] = 0;
  create2[1] = 7;
  annihilate2[1] = 5;
  expectedEnergy += interParameter[0] *
                    brute_green(2, create2, annihilate2, baseOverlap, projCnt);
  create3[0] = 1;
  annihilate3[0] = 0;
  create3[1] = 3;
  annihilate3[1] = 2;
  create3[2] = 4;
  annihilate3[2] = 5;
  expectedEnergy += nbodyParameter[0] *
                    brute_green(3, create3, annihilate3, baseOverlap, projCnt);
  actualEnergy = CalculateHamiltonianGC(baseOverlap, ACTIVE, eleIdx, eleCfg,
                                        eleNum, projCnt);
  CHECK(cabs(actualEnergy - expectedEnergy) < 2.0e-8,
        "Hamiltonian got=(%.17g,%.17g) expected=(%.17g,%.17g)",
        creal(actualEnergy), cimag(actualEnergy), creal(expectedEnergy),
        cimag(expectedEnergy));
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
  create1[0] = 1;
  annihilate1[0] = 0;
  expectedLocal0 = brute_green(1, create1, annihilate1, baseOverlap, projCnt);
  create1[0] = 7;
  annihilate1[0] = 3;
  expectedLocal1 = brute_green(1, create1, annihilate1, baseOverlap, projCnt);
  create2[0] = 1;
  annihilate2[0] = 0;
  create2[1] = 7;
  annihilate2[1] = 5;
  expectedDC = brute_green(2, create2, annihilate2, baseOverlap, projCnt);
  create3[0] = 1;
  annihilate3[0] = 0;
  create3[1] = 3;
  annihilate3[1] = 2;
  create3[2] = 4;
  annihilate3[2] = 5;
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
