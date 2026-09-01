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

#include "../src/mVMC/gc_size.c"
#include "../src/mVMC/workspace.c"
#include "../src/mVMC/matrix_gc.c"
#include "../src/mVMC/slater_gc.c"

#define ORBITALS 4
#define PARAMETERS 6

static int failures = 0;
static int orbitalIdxStorage[ORBITALS][ORBITALS];
static int orbitalSgnStorage[ORBITALS][ORBITALS];
static int *orbitalIdxRows[ORBITALS];
static int *orbitalSgnRows[ORBITALS];
static int qpTransStorage[2] = {0, 1};
static int qpTransSgnStorage[2] = {1, 1};
static int qpOptTransStorage[2] = {0, 1};
static int qpOptTransSgnStorage[2] = {1, 1};
static int *qpTransRows[1] = {qpTransStorage};
static int *qpTransSgnRows[1] = {qpTransSgnStorage};
static int *qpOptTransRows[1] = {qpOptTransStorage};
static int *qpOptTransSgnRows[1] = {qpOptTransSgnStorage};

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "GCSlaterDiff_Unit FAIL: ");                            \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, "\n");                                                  \
      failures++;                                                               \
    }                                                                           \
  } while (0)

static void refresh_slater_elements(void) {
  int rsi;
  for (rsi = 0; rsi < ORBITALS; rsi++) {
    int rsj;
    for (rsj = 0; rsj < ORBITALS; rsj++) {
      const int direct = OrbitalIdx[rsi][rsj];
      const int reverse = OrbitalIdx[rsj][rsi];
      SlaterElm[(size_t)rsi * ORBITALS + (size_t)rsj] =
          Slater[direct] * (double)OrbitalSgn[rsi][rsj] -
          Slater[reverse] * (double)OrbitalSgn[rsj][rsi];
    }
  }
}

static double complex rebuild_overlap(const int ncur, const int *eleIdx) {
  refresh_slater_elements();
  CHECK(CalculateMAllGC_fcmp(ncur, eleIdx, 0, 1) == GC_MALL_OK,
        "rebuild failed for ncur=%d", ncur);
  return QPFullWeight[0] * PfM[0];
}

static double complex finite_difference(const int ncur, const int *eleIdx,
                                        const int parameter,
                                        const int imaginary,
                                        const double epsilon,
                                        const double complex baseOverlap) {
  const double complex delta = imaginary ? I * epsilon : epsilon;
  double complex plus;
  double complex minus;
  Slater[parameter] += delta;
  plus = rebuild_overlap(ncur, eleIdx);
  Slater[parameter] -= 2.0 * delta;
  minus = rebuild_overlap(ncur, eleIdx);
  Slater[parameter] += delta;
  return (plus - minus) / (2.0 * epsilon * baseOverlap);
}

static void check_derivative(const int ncur, const int *eleIdx,
                             const int realParameter,
                             const int imaginaryParameter) {
  double complex derivative[2 * PARAMETERS];
  double complex baseOverlap = rebuild_overlap(ncur, eleIdx);
  double complex expectedReal;
  double complex expectedImag;
  memset(derivative, 0x5a, sizeof(derivative));
  SlaterElmDiffGC_fcmp(derivative, baseOverlap, eleIdx, ncur);
  expectedReal = finite_difference(ncur, eleIdx, realParameter, 0, 1.0e-5,
                                   baseOverlap);
  expectedImag = finite_difference(ncur, eleIdx, imaginaryParameter, 1, 5.0e-6,
                                   baseOverlap);
  CHECK(cabs(derivative[2 * realParameter] - expectedReal) < 2.0e-9,
        "real FD ncur=%d parameter=%d got=(%.17g,%.17g) expected=(%.17g,%.17g)",
        ncur, realParameter, creal(derivative[2 * realParameter]),
        cimag(derivative[2 * realParameter]), creal(expectedReal),
        cimag(expectedReal));
  CHECK(cabs(derivative[2 * imaginaryParameter + 1] - expectedImag) < 2.0e-9,
        "imag FD ncur=%d parameter=%d got=(%.17g,%.17g) expected=(%.17g,%.17g)",
        ncur, imaginaryParameter,
        creal(derivative[2 * imaginaryParameter + 1]),
        cimag(derivative[2 * imaginaryParameter + 1]), creal(expectedImag),
        cimag(expectedImag));
  if (ncur == 0) {
    int i;
    for (i = 0; i < 2 * PARAMETERS; i++) {
      CHECK(derivative[i] == 0.0, "vacuum derivative[%d] is nonzero", i);
    }
  }
}

static void initialize_fixture(void) {
  int row;
  int parameter = 0;
  NThread = 1;
  Nsite = 2;
  Nsite2 = ORBITALS;
  NsizeMax = ORBITALS;
  NQPFull = 1;
  NQPFix = 1;
  NMPTrans = 1;
  NSPGaussLeg = 1;
  NQPOptTrans = 1;
  NSlater = PARAMETERS;
  LapackLWork = 128;
  for (row = 0; row < ORBITALS; row++) {
    int column;
    orbitalIdxRows[row] = orbitalIdxStorage[row];
    orbitalSgnRows[row] = orbitalSgnStorage[row];
    orbitalIdxStorage[row][row] = 0;
    orbitalSgnStorage[row][row] = 0;
    for (column = row + 1; column < ORBITALS; column++) {
      orbitalIdxStorage[row][column] = parameter;
      orbitalIdxStorage[column][row] = parameter;
      orbitalSgnStorage[row][column] = 1;
      orbitalSgnStorage[column][row] = -1;
      parameter++;
    }
  }
  OrbitalIdx = orbitalIdxRows;
  OrbitalSgn = orbitalSgnRows;
  QPTrans = qpTransRows;
  QPTransSgn = qpTransSgnRows;
  QPOptTrans = qpOptTransRows;
  QPOptTransSgn = qpOptTransSgnRows;
  Slater = malloc(PARAMETERS * sizeof(*Slater));
  SlaterElm = malloc(ORBITALS * ORBITALS * sizeof(*SlaterElm));
  InvM = malloc(ORBITALS * ORBITALS * sizeof(*InvM));
  PfM = malloc(sizeof(*PfM));
  QPFullWeight = malloc(sizeof(*QPFullWeight));
  CHECK(Slater != NULL && SlaterElm != NULL && InvM != NULL && PfM != NULL &&
            QPFullWeight != NULL,
        "fixture allocation failed");
  for (parameter = 0; parameter < PARAMETERS; parameter++) {
    Slater[parameter] =
        (0.37 + 0.11 * parameter) + (0.19 - 0.047 * parameter) * I;
  }
  QPFullWeight[0] = 0.73 - 0.21 * I;
  initializeWorkSpaceAll();
}

static void free_fixture(void) {
  FreeWorkSpaceAll();
  free(QPFullWeight);
  free(PfM);
  free(InvM);
  free(SlaterElm);
  free(Slater);
}

int main(int argc, char **argv) {
  const int eleIdx2[2] = {0, 1};
  const int eleIdx4[4] = {0, 1, 2, 3};
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  initialize_fixture();
  check_derivative(0, NULL, 0, 1);
  check_derivative(2, eleIdx2, 0, 0);
  check_derivative(4, eleIdx4, 1, 4);
  free_fixture();
#ifdef _mpi_use
  MPI_Finalize();
#endif
  if (failures != 0) {
    fprintf(stderr, "%d GC Slater derivative checks failed\n", failures);
    return EXIT_FAILURE;
  }
  printf("GC Slater derivative checks: PASS\n");
  return EXIT_SUCCESS;
}
