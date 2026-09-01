#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

extern int omp_get_thread_num(void);

#include "global.h"
#include "blas_externs.h"

/* Use the production implementations in one translation unit because
 * global.h owns the legacy process-wide storage definitions. */
#include "../src/mVMC/gc_size.c"
#include "../src/mVMC/workspace.c"
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
    build_active_x(ncur, eleIdx, matrix);
    expectedPf = pfaffian_recursive(matrix, ncur);
    snprintf(label, sizeof(label), "M=%d n=%d mask=%u", dimension, ncur,
             mask);

    CHECK(CalculateMAllGC_fcmp(ncur, eleIdx, 0, 1) == GC_MALL_OK,
          "%s rebuild status", label);
    CHECK(cabs(PfM[0] - expectedPf) < 2.0e-10,
          "%s Pf got=(%.17g,%.17g) expected=(%.17g,%.17g)", label,
          creal(PfM[0]), cimag(PfM[0]), creal(expectedPf), cimag(expectedPf));
    CHECK(PfM[1] == -302.0 + 18.0 * I, "%s changed unused Pf slot", label);
    if (ncur > 0) check_inverse(matrix, ncur, label);
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
  build_active_x(2, eleIdx, matrix);
  /* build_active_x used QP 0; the second Slater block is scaled by 1.7. */
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

static void test_failure_does_not_commit(void) {
  const int eleIdx[4] = {0, 1, 2, 3};
  const double complex invSentinel = 27.0 - 4.0 * I;
  const double complex pfSentinel = -11.0 + 8.0 * I;
  const size_t block = (size_t)NsizeMax * (size_t)NsizeMax;
  size_t i;
  memset(SlaterElm, 0, 2U * block * sizeof(*SlaterElm));
  for (i = 0; i < 2 * block; i++) InvM[i] = invSentinel;
  PfM[0] = pfSentinel;
  PfM[1] = pfSentinel;
  CHECK(CalculateMAllGC_fcmp(4, eleIdx, 0, 1) == GC_MALL_PFAFFIAN,
        "singular matrix reports Pfaffian stage");
  CHECK(PfM[0] == pfSentinel, "singular matrix committed Pfaffian");
  for (i = 0; i < 2 * block; i++) {
    CHECK(InvM[i] == invSentinel,
          "singular matrix committed inverse index=%zu", i);
  }
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
  NThread = 1;
  initializeWorkSpaceAll();
  test_dimension_patterns(4, 0);
  test_dimension_patterns(4, 2);
  test_dimension_patterns(4, 4);
  test_qp_range_layout();
  test_failure_does_not_commit();
  test_invalid_arguments();
  test_dimension_patterns(6, 2);
  test_dimension_patterns(6, 4);
  test_dimension_patterns(6, 6);
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
