#include <complex.h>
#include <math.h>
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

/* Use the production implementations in one translation unit because
 * global.h owns the legacy process-wide storage definitions. */
#include "../src/mVMC/gc_size.c"
#include "../src/mVMC/workspace.c"
#include "../src/mVMC/gc_component.c"
#include "../src/mVMC/matrix_gc.c"

#define MAX_DIMENSION 6

static int failures = 0;

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "GCMatrix_Unit FAIL: ");                                 \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, "\n");                                                   \
      failures++;                                                               \
    }                                                                           \
  } while (0)

static void reset_component(void) {
  CHECK(GCComponentInit() == GC_MALL_OK, "component init");
}

static int popcount(unsigned int value) {
  int count = 0;
  while (value != 0U) {
    count += (int)(value & 1U);
    value >>= 1;
  }
  return count;
}

static void fill_slater(double complex *slater, const int dimension,
                        const double scale) {
  int row;
  memset(slater, 0,
         (size_t)dimension * (size_t)dimension * sizeof(*slater));
  for (row = 0; row < dimension; row++) {
    int column;
    for (column = row + 1; column < dimension; column++) {
      const double realPart =
          scale * (0.31 * (row + 1) + 0.17 * (column + 1) +
                   0.013 * (row + 1) * (column + 1));
      const double imagPart =
          scale * (0.09 * (row + 1) - 0.07 * (column + 1) + 0.023);
      const double complex value = realPart + imagPart * I;
      slater[(size_t)row * (size_t)dimension + (size_t)column] = value;
      slater[(size_t)column * (size_t)dimension + (size_t)row] = -value;
    }
  }
}

static void build_active_x(const int ncur, const int *eleIdx,
                           double complex *matrix) {
  int row;
  for (row = 0; row < ncur; row++) {
    int column;
    for (column = 0; column < ncur; column++) {
      matrix[(size_t)row * (size_t)ncur + (size_t)column] =
          -SlaterElm[(size_t)eleIdx[row] * (size_t)Nsite2 +
                     (size_t)eleIdx[column]];
    }
  }
}

/* Pf(SlaterElm_x) is the stored PfM convention (|phi_GC> = exp[sum f c+c+]|0>),
 * while InvM keeps the X = -SlaterElm convention of build_active_x. */
static void build_active_slater(const int ncur, const int *eleIdx,
                                double complex *matrix) {
  int row;
  for (row = 0; row < ncur; row++) {
    int column;
    for (column = 0; column < ncur; column++) {
      matrix[(size_t)row * (size_t)ncur + (size_t)column] =
          SlaterElm[(size_t)eleIdx[row] * (size_t)Nsite2 +
                    (size_t)eleIdx[column]];
    }
  }
}

static double complex pfaffian_recursive(const double complex *matrix,
                                          const int n) {
  double complex result = 0.0;
  int column;
  if (n == 0) return 1.0;
  if (n == 2) return matrix[1];
  for (column = 1; column < n; column++) {
    double complex minor[MAX_DIMENSION * MAX_DIMENSION];
    int sourceRow;
    int minorRow = 0;
    for (sourceRow = 1; sourceRow < n; sourceRow++) {
      int sourceColumn;
      int minorColumn = 0;
      if (sourceRow == column) continue;
      for (sourceColumn = 1; sourceColumn < n; sourceColumn++) {
        if (sourceColumn == column) continue;
        minor[(size_t)minorRow * (size_t)(n - 2) +
              (size_t)minorColumn] =
            matrix[(size_t)sourceRow * (size_t)n +
                   (size_t)sourceColumn];
        minorColumn++;
      }
      minorRow++;
    }
    result += ((column % 2) != 0 ? 1.0 : -1.0) * matrix[column] *
              pfaffian_recursive(minor, n - 2);
  }
  return result;
}

static void check_inverse(const double complex *matrix, const int ncur,
                          const char *label) {
  int row;
  for (row = 0; row < ncur; row++) {
    int column;
    for (column = 0; column < ncur; column++) {
      double complex product = 0.0;
      int inner;
      for (inner = 0; inner < ncur; inner++) {
        product += matrix[(size_t)row * (size_t)ncur + (size_t)inner] *
                   InvM[(size_t)inner * (size_t)NsizeMax +
                        (size_t)column];
      }
      CHECK(cabs(product - (row == column ? 1.0 : 0.0)) < 2.0e-10,
            "%s inverse row=%d col=%d value=(%.17g,%.17g)", label, row,
            column, creal(product), cimag(product));
    }
  }
}

static void check_inactive_sentinel(const int ncur,
                                    const double complex sentinel,
                                    const char *label) {
  int row;
  for (row = 0; row < NsizeMax; row++) {
    int column;
    for (column = 0; column < NsizeMax; column++) {
      if (row >= ncur || column >= ncur) {
        CHECK(InvM[(size_t)row * (size_t)NsizeMax + (size_t)column] ==
                  sentinel,
              "%s changed inactive row=%d col=%d", label, row, column);
      }
    }
  }
}

static void legacy_inverse(const double complex *x, const int ncur,
                           double complex *invRowMajor) {
  double complex *work =
      malloc((size_t)ncur * (size_t)ncur * sizeof(*work));
  double complex *lwork = malloc(256 * sizeof(*lwork));
  int *ipiv = malloc((size_t)ncur * sizeof(*ipiv));
  int n = ncur, lda = ncur, info = 0, lw = 256, row, col;
  for (row = 0; row < ncur; row++)
    for (col = 0; col < ncur; col++)
      work[row + (size_t)col * (size_t)ncur] =
          x[(size_t)row * (size_t)ncur + (size_t)col];
  M_ZGETRF(&n, &n, work, &lda, ipiv, &info);
  CHECK(info == 0, "legacy GETRF info=%d", info);
  M_ZGETRI(&n, work, &lda, ipiv, lwork, &lw, &info);
  CHECK(info == 0, "legacy GETRI info=%d", info);
  for (row = 0; row < ncur; row++)
    for (col = 0; col < ncur; col++)
      invRowMajor[(size_t)row * (size_t)NsizeMax + (size_t)col] =
          -work[col + (size_t)row * (size_t)ncur];
  free(ipiv);
  free(lwork);
  free(work);
}

static void check_regular_matches_legacy(const int ncur, const int *eleIdx,
                                         const char *label) {
  double complex x[MAX_DIMENSION * MAX_DIMENSION];
  double complex legacy[MAX_DIMENSION * MAX_DIMENSION];
  int row, col;
  build_active_x(ncur, eleIdx, x);
  legacy_inverse(x, ncur, legacy);
  CHECK(GCQpState[0] == GC_QP_INVERSE,
        "%s state is not INVERSE", label);
  for (row = 0; row < ncur; row++)
    for (col = 0; col < ncur; col++) {
      const size_t k =
          (size_t)row * (size_t)NsizeMax + (size_t)col;
      CHECK(cabs(InvM[k] - legacy[k]) <=
                1.0e-12 * (1.0 + cabs(legacy[k])),
            "%s inverse row=%d col=%d", label, row, col);
    }
}

static void test_dimension_patterns(const int dimension, const int ncur) {
  const double complex sentinel = 73.0 - 19.0 * I;
  unsigned int mask;
  Nsite2 = dimension;
  NsizeMax = dimension;
  NQPFull = 2;
  LapackLWork = 256;
  SlaterElm = realloc(
      SlaterElm, 2U * (size_t)dimension * (size_t)dimension *
                     sizeof(*SlaterElm));
  InvM = realloc(InvM,
                 (2U * (size_t)dimension * (size_t)dimension + 2U) *
                     sizeof(*InvM));
  CHECK(SlaterElm != NULL && InvM != NULL, "test allocation dimension=%d",
        dimension);
  if (SlaterElm == NULL || InvM == NULL) exit(EXIT_FAILURE);
  PfM = InvM + 2U * (size_t)dimension * (size_t)dimension;
  reset_component();
  fill_slater(SlaterElm, dimension, 1.0);
  fill_slater(SlaterElm + (size_t)dimension * (size_t)dimension, dimension,
              1.7);

  for (mask = 0; mask < (1U << dimension); mask++) {
    int eleIdx[MAX_DIMENSION];
    double complex matrix[MAX_DIMENSION * MAX_DIMENSION];
    double complex expectedPf;
    int count = 0;
    int rs;
    char label[80];
    if (popcount(mask) != ncur) continue;
    for (rs = 0; rs < dimension; rs++) {
      if ((mask & (1U << rs)) != 0U) eleIdx[count++] = rs;
    }
    for (rs = 0; rs < 2 * dimension * dimension; rs++) InvM[rs] = sentinel;
    PfM[0] = -301.0 + 17.0 * I;
    PfM[1] = -302.0 + 18.0 * I;
    build_active_slater(ncur, eleIdx, matrix);
    expectedPf = pfaffian_recursive(matrix, ncur);
    build_active_x(ncur, eleIdx, matrix);
    snprintf(label, sizeof(label), "M=%d n=%d mask=%u", dimension, ncur,
             mask);

    CHECK(CalculateMAllGC_fcmp(ncur, eleIdx, 0, 1) == GC_MALL_OK,
          "%s rebuild status", label);
    CHECK(cabs(PfM[0] - expectedPf) < 2.0e-10,
          "%s Pf got=(%.17g,%.17g) expected=(%.17g,%.17g)", label,
          creal(PfM[0]), cimag(PfM[0]), creal(expectedPf), cimag(expectedPf));
    CHECK(PfM[1] == -302.0 + 18.0 * I, "%s changed unused Pf slot", label);
    if (ncur > 0) {
      check_inverse(matrix, ncur, label);
      check_regular_matches_legacy(ncur, eleIdx, label);
    }
    check_inactive_sentinel(ncur, sentinel, label);
  }
}

static void test_qp_range_layout(void) {
  const int eleIdx[2] = {1, 3};
  double complex matrix[4];
  double complex expectedPf;
  const size_t block = (size_t)NsizeMax * (size_t)NsizeMax;
  size_t i;
  for (i = 0; i < 2 * block; i++) InvM[i] = 91.0 + 3.0 * I;
  PfM[0] = 41.0 + I;
  PfM[1] = 42.0 + I;
  build_active_slater(2, eleIdx, matrix);
  /* build_active_slater used QP 0; the second Slater block is scaled by 1.7. */
  expectedPf = 1.7 * pfaffian_recursive(matrix, 2);
  CHECK(CalculateMAllGC_fcmp(2, eleIdx, 1, 2) == GC_MALL_OK,
        "qpStart=1 rebuild status");
  CHECK(cabs(PfM[0] - expectedPf) < 2.0e-10,
        "rank-local Pf slot reads global qpStart");
  CHECK(PfM[1] == 42.0 + I, "rank-local rebuild changed Pf slot 1");
  for (i = block; i < 2 * block; i++) {
    CHECK(InvM[i] == 91.0 + 3.0 * I,
          "rank-local rebuild changed matrix slot 1 index=%zu", i - block);
  }
}

static void test_singular_commits_pfonly(void) {
  const int eleIdx[4] = {0, 1, 2, 3};
  const double complex invSentinel = 27.0 - 4.0 * I;
  const size_t block = (size_t)NsizeMax * (size_t)NsizeMax;
  size_t i;
  reset_component();
  memset(SlaterElm, 0, 2U * block * sizeof(*SlaterElm));
  for (i = 0; i < 2 * block; i++) InvM[i] = invSentinel;
  PfM[0] = -11.0 + 8.0 * I;
  CHECK(CalculateMAllGC_fcmp(4, eleIdx, 0, 1) == GC_MALL_OK,
        "singular matrix must classify, not fail");
  CHECK(GCQpState[0] == GC_QP_PFONLY, "singular matrix state");
  CHECK(PfM[0] == 0.0, "singular matrix PfM must be exactly 0");
  for (i = 0; i < 2 * block; i++)
    CHECK(InvM[i] == invSentinel,
          "singular matrix touched inverse %zu", i);
}

static void test_structural_zero_row(void) {
  const int with0[2] = {0, 2};
  const int without0[2] = {1, 3};
  int k;
  reset_component();
  fill_slater(SlaterElm, NsizeMax, 1.0);
  for (k = 0; k < NsizeMax; k++) {
    SlaterElm[k] = 0.0;
    SlaterElm[(size_t)k * (size_t)NsizeMax] = 0.0;
  }
  CHECK(CalculateMAllGC_fcmp(2, with0, 0, 1) == GC_MALL_OK,
        "zero-row rebuild");
  CHECK(GCQpState[0] == GC_QP_PFONLY && PfM[0] == 0.0,
        "zero-row state");
  CHECK(CalculateMAllGC_fcmp(2, without0, 0, 1) == GC_MALL_OK,
        "regular rebuild");
  CHECK(GCQpState[0] == GC_QP_INVERSE,
        "regular state after zero-row");
  check_regular_matches_legacy(2, without0, "after zero-row");
}

static void test_low_rank_is_pfonly(void) {
  const double complex u[4] = {0.7 + 0.1 * I, -0.3 + 0.4 * I,
                               0.5 - 0.2 * I, 0.9 + 0.3 * I};
  const double complex v[4] = {0.2 - 0.6 * I, 0.8 + 0.1 * I,
                               -0.4 + 0.5 * I, 0.1 + 0.7 * I};
  const int all[4] = {0, 1, 2, 3};
  const int pair[2] = {0, 1};
  int i, j;
  reset_component();
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      SlaterElm[(size_t)i * (size_t)NsizeMax + (size_t)j] =
          u[i] * v[j] - v[i] * u[j];
  CHECK(CalculateMAllGC_fcmp(4, all, 0, 1) == GC_MALL_OK,
        "low-rank rebuild");
  CHECK(GCQpState[0] == GC_QP_PFONLY, "low-rank state (got %s)",
        GCComponentStateName(GCQpState[0]));
  CHECK(cabs(PfM[0]) <= 1.0e-12, "low-rank Pf magnitude %.3e",
        cabs(PfM[0]));
  CHECK(CalculateMAllGC_fcmp(2, pair, 0, 1) == GC_MALL_OK,
        "low-rank pair");
  CHECK(GCQpState[0] == GC_QP_INVERSE,
        "two-electron low-rank is regular");
}

static void test_mixed_components(void) {
  const int eleIdx[2] = {0, 3};
  const size_t block = (size_t)NsizeMax * (size_t)NsizeMax;
  int k;
  reset_component();
  fill_slater(SlaterElm, NsizeMax, 1.0);
  fill_slater(SlaterElm + block, NsizeMax, 1.7);
  for (k = 0; k < NsizeMax; k++) {
    SlaterElm[block + (size_t)k] = 0.0;
    SlaterElm[block + (size_t)k * (size_t)NsizeMax] = 0.0;
  }
  CHECK(CalculateMAllGC_fcmp(2, eleIdx, 0, 2) == GC_MALL_OK,
        "mixed rebuild");
  CHECK(GCQpState[0] == GC_QP_INVERSE &&
            GCQpState[1] == GC_QP_PFONLY,
        "mixed states");
  CHECK(PfM[1] == 0.0, "mixed PfM[1]");
  CHECK(CalculateMAllGC_fcmp(2, eleIdx, 1, 2) == GC_MALL_OK,
        "range rebuild");
  CHECK(GCQpState[0] == GC_QP_PFONLY,
        "range slot 0 reads global qp1");
}

static void test_nonfinite_is_error(void) {
  const int eleIdx[2] = {0, 1};
  const double complex pfSentinel = 5.0 + 5.0 * I;
  reset_component();
  fill_slater(SlaterElm, NsizeMax, 1.0);
  SlaterElm[1] = NAN;
  SlaterElm[(size_t)NsizeMax] = NAN;
  PfM[0] = pfSentinel;
  GCQpState[0] = GC_QP_INVERSE;
  CHECK(CalculateMAllGC_fcmp(2, eleIdx, 0, 1) == GC_MALL_NONFINITE,
        "NaN status");
  CHECK(PfM[0] == pfSentinel, "NaN committed PfM");
}

static void test_vacuum_component(void) {
  reset_component();
  CHECK(CalculateMAllGC_fcmp(0, NULL, 0, 1) == GC_MALL_OK,
        "vacuum rebuild");
  CHECK(PfM[0] == 1.0 && GCQpState[0] == GC_QP_INVERSE,
        "vacuum Pf/state");
}

static void test_invalid_arguments(void) {
  const int valid[2] = {0, 1};
  const int bad[2] = {0, 99};
  CHECK(CalculateMAllGC_fcmp(-2, valid, 0, 1) == GC_MALL_INVALID_ARGUMENT,
        "negative ncur accepted");
  CHECK(CalculateMAllGC_fcmp(1, valid, 0, 1) == GC_MALL_INVALID_ARGUMENT,
        "odd ncur accepted");
  CHECK(CalculateMAllGC_fcmp(NsizeMax + 2, valid, 0, 1) ==
            GC_MALL_INVALID_ARGUMENT,
        "oversized ncur accepted");
  CHECK(CalculateMAllGC_fcmp(2, bad, 0, 1) == GC_MALL_INVALID_ARGUMENT,
        "bad fused orbital accepted");
  CHECK(CalculateMAllGC_fcmp(2, valid, -1, 1) == GC_MALL_INVALID_ARGUMENT,
        "negative qpStart accepted");
  CHECK(CalculateMAllGC_fcmp(2, valid, 1, 3) == GC_MALL_INVALID_ARGUMENT,
        "qpEnd beyond NQPFull accepted");
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  NThread = omp_get_max_threads();
  initializeWorkSpaceAll();
  test_dimension_patterns(4, 0);
  test_dimension_patterns(4, 2);
  test_dimension_patterns(4, 4);
  test_qp_range_layout();
  test_singular_commits_pfonly();
  test_structural_zero_row();
  test_low_rank_is_pfonly();
  test_mixed_components();
  test_nonfinite_is_error();
  test_vacuum_component();
  test_invalid_arguments();
  test_dimension_patterns(6, 2);
  test_dimension_patterns(6, 4);
  test_dimension_patterns(6, 6);
  GCComponentRelease();
  FreeWorkSpaceAll();
  free(InvM);
  free(SlaterElm);
#ifdef _mpi_use
  MPI_Finalize();
#endif
  if (failures != 0) {
    fprintf(stderr, "%d GC matrix checks failed\n", failures);
    return EXIT_FAILURE;
  }
  printf("GC matrix checks: PASS\n");
  return EXIT_SUCCESS;
}
