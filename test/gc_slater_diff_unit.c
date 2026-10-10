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

#include "../src/mVMC/gc_size.c"
#include "../src/mVMC/workspace.c"
#include "../src/mVMC/matrix_gc.c"
#include "../src/mVMC/slater_gc.c"
#include "../src/mVMC/slater_fsz.c"

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

static void refresh_slater_elements(void);
/* Matrix builder used by rebuild_overlap(): the test's own expansion or the
 * production UpdateSlaterElm_fsz(). */
static void (*rebuild_matrix)(void) = refresh_slater_elements;

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
  rebuild_matrix();
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
    for (i = 0; i < 2 * NSlater; i++) {
      CHECK(derivative[i] == 0.0, "vacuum derivative[%d] is nonzero", i);
    }
  }
}

static void initialize_fixture(void) {
  int row;
  int parameter = 0;
  NThread = omp_get_max_threads();
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

/*
 * Anti-periodic-style table: parameter 1 is shared by a +1 and a -1
 * upper-triangle row, as a wrapped pair shares a translation class with an
 * unwrapped one, and parameter 2 enters with -1 only.  Two negative rows that
 * merely flip the overall sign of Pf(F) would be invisible to d ln Pf, so
 * they are placed on terms of the full Pfaffian with different weight.
 * Rows are (I, J, parameter, sign) with I<J.
 */
#define SIGNED_PARAMETERS 4
static const int signedRows[6][4] = {
    {0, 1, 0, 1}, {2, 3, 0, 1}, {0, 2, 1, 1},
    {1, 3, 1, -1}, {0, 3, 2, -1}, {1, 2, 3, 1}};

static void load_pair_table(const int rows[6][4], const int expanded) {
  int row;
  for (row = 0; row < ORBITALS; row++) {
    orbitalIdxStorage[row][row] = 0;
    orbitalSgnStorage[row][row] = 0;
  }
  for (row = 0; row < 6; row++) {
    const int first = rows[row][0];
    const int second = rows[row][1];
    const int parameter = expanded ? row : rows[row][2];
    const int sign = expanded ? 1 : rows[row][3];
    orbitalIdxStorage[first][second] = parameter;
    orbitalIdxStorage[second][first] = parameter;
    orbitalSgnStorage[first][second] = sign;
    orbitalSgnStorage[second][first] = -sign;
  }
}

static void set_signed_parameters(void) {
  int parameter;
  NSlater = SIGNED_PARAMETERS;
  load_pair_table(signedRows, 0);
  for (parameter = 0; parameter < SIGNED_PARAMETERS; parameter++) {
    Slater[parameter] =
        (0.41 - 0.13 * parameter) + (-0.23 + 0.17 * parameter) * I;
  }
}

/* Independent expansion F_IJ = q_I q_J (s_IJ f[k] - s_JI f[k]) for I != J,
 * with the reader's completion s_JI = -s_IJ and the site sign q. */
static double complex expected_element(const int first, const int second) {
  int row;
  for (row = 0; row < 6; row++) {
    const int a = signedRows[row][0];
    const int b = signedRows[row][1];
    if ((a == first && b == second) || (a == second && b == first)) {
      const double complex f = Slater[signedRows[row][2]];
      const double sign =
          (double)(a == first ? signedRows[row][3] : -signedRows[row][3]);
      const double q = (double)(qpTransSgnStorage[first % 2] *
                                qpTransSgnStorage[second % 2]);
      return q * (sign * f - (-sign) * f);
    }
  }
  return 0.0;
}

static void check_production_matrix(const char *label) {
  int first;
  UpdateSlaterElm_fsz();
  for (first = 0; first < ORBITALS; first++) {
    int second;
    for (second = 0; second < ORBITALS; second++) {
      const double complex actual =
          SlaterElm[(size_t)first * ORBITALS + (size_t)second];
      const double complex expected = expected_element(first, second);
      CHECK(actual == expected,
            "%s production F[%d][%d]=(%.17g,%.17g) expected=(%.17g,%.17g)",
            label, first, second, creal(actual), cimag(actual),
            creal(expected), cimag(expected));
    }
  }
}

/* Shared signed parameters versus one parameter per pair with the sign
 * absorbed into its value: the production matrices must be bitwise equal. */
static void check_expanded_reference(void) {
  double complex shared[ORBITALS * ORBITALS];
  double complex sharedParameters[SIGNED_PARAMETERS];
  int row;
  set_signed_parameters();
  UpdateSlaterElm_fsz();
  memcpy(shared, SlaterElm, sizeof(shared));
  memcpy(sharedParameters, Slater, sizeof(sharedParameters));
  NSlater = 6;
  load_pair_table(signedRows, 1);
  for (row = 0; row < 6; row++) {
    Slater[row] =
        (double)signedRows[row][3] * sharedParameters[signedRows[row][2]];
  }
  UpdateSlaterElm_fsz();
  CHECK(memcmp(shared, SlaterElm, sizeof(shared)) == 0,
        "expanded reference SlaterElm differs from shared signed table");
}

static void check_signed_derivatives(const char *label) {
  const int eleIdx03[2] = {0, 3};
  const int eleIdx23[2] = {2, 3};
  const int eleIdx13[2] = {1, 3};
  const int eleIdx4[4] = {0, 1, 2, 3};
  check_derivative(0, NULL, 0, 1);
  check_derivative(2, eleIdx03, 2, 2);
  check_derivative(2, eleIdx23, 0, 0);
  check_derivative(2, eleIdx13, 1, 1);
  check_derivative(4, eleIdx4, 0, 1);
  check_derivative(4, eleIdx4, 1, 0);
  check_derivative(4, eleIdx4, 2, 3);
  (void)label;
}

/* The sign must change the full-occupation derivative, otherwise the signed
 * checks above could not detect an ignored OrbitalSgn. */
static void check_sign_is_observable(void) {
  const int eleIdx4[4] = {0, 1, 2, 3};
  double complex withSign[2 * SIGNED_PARAMETERS];
  double complex withoutSign[2 * SIGNED_PARAMETERS];
  int row;
  set_signed_parameters();
  SlaterElmDiffGC_fcmp(withSign, rebuild_overlap(4, eleIdx4), eleIdx4, 4);
  for (row = 0; row < 6; row++) {
    if (signedRows[row][3] < 0) {
      orbitalSgnStorage[signedRows[row][0]][signedRows[row][1]] = 1;
      orbitalSgnStorage[signedRows[row][1]][signedRows[row][0]] = -1;
    }
  }
  SlaterElmDiffGC_fcmp(withoutSign, rebuild_overlap(4, eleIdx4), eleIdx4, 4);
  CHECK(cabs(withSign[2] - withoutSign[2]) > 1.0e-3,
        "dropping the negative OrbitalSgn does not change d/df1");
  CHECK(cabs(withSign[4] - withoutSign[4]) > 1.0e-3,
        "dropping the negative OrbitalSgn does not change d/df2");
  set_signed_parameters();
}

static void run_signed_checks(void) {
  rebuild_matrix = UpdateSlaterElm_fsz;
  set_signed_parameters();
  check_production_matrix("identity translation");
  check_signed_derivatives("identity translation");
  check_sign_is_observable();
  check_expanded_reference();

  /* One translation projection whose map is the identity but whose sign
   * depends on the site: rows and columns of site 1 flip. */
  set_signed_parameters();
  qpTransSgnStorage[1] = -1;
  check_production_matrix("site-dependent translation sign");
  check_signed_derivatives("site-dependent translation sign");
  qpTransSgnStorage[1] = 1;
  rebuild_matrix = refresh_slater_elements;
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
  run_signed_checks();
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
