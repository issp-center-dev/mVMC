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
#include "../src/mVMC/gc_component.c"
#include "../src/mVMC/matrix_gc.c"
#include "../src/mVMC/slater_gc.c"
#include "../src/mVMC/slater_fsz.c"
#include "../src/mVMC/slater.c"

#define ORBITALS 4
#define PARAMETERS 6

static int failures = 0;
static int orbitalIdxStorage[ORBITALS][ORBITALS];
static int orbitalSgnStorage[ORBITALS][ORBITALS];
static int *orbitalIdxRows[ORBITALS];
static int *orbitalSgnRows[ORBITALS];
static int qpTransStorage[2] = {0, 1};
static int qpTransSgnStorage[2] = {1, 1};
static int qpTransStorage1[2] = {1, 0};
static int qpTransSgnStorage1[2] = {1, 1};
static int qpOptTransStorage[2] = {0, 1};
static int qpOptTransSgnStorage[2] = {1, 1};
static int *qpTransRows[1] = {qpTransStorage};
static int *qpTransSgnRows[1] = {qpTransSgnStorage};
static int *qpTransRows2[2] = {qpTransStorage, qpTransStorage1};
static int *qpTransSgnRows2[2] = {qpTransSgnStorage,
                                  qpTransSgnStorage1};
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

static int finite_close(const double complex actual,
                        const double complex expected,
                        const double tolerance) {
  return isfinite(creal(actual)) && isfinite(cimag(actual)) &&
         isfinite(creal(expected)) && isfinite(cimag(expected)) &&
         isfinite(tolerance) && tolerance >= 0.0 &&
         cabs(actual - expected) <= tolerance;
}

static void refresh_slater_elements(void);
static void free_fixture(void);
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
  double complex overlap = 0.0;
  int qp;
  rebuild_matrix();
  CHECK(CalculateMAllGC_fcmp(ncur, eleIdx, 0, NQPFull) == GC_MALL_OK,
        "rebuild failed for ncur=%d", ncur);
  for (qp = 0; qp < NQPFull; qp++)
    overlap += QPFullWeight[qp] * PfM[qp];
  return overlap;
}

static double complex finite_difference(const int ncur, const int *eleIdx,
                                        const int parameter,
                                        const int imaginary,
                                        const double epsilon,
                                        const double complex baseOverlap) {
  const double complex original = Slater[parameter];
  const double complex delta = imaginary ? I * epsilon : epsilon;
  double complex plus;
  double complex minus;
  double complex restored;
  Slater[parameter] = original + delta;
  plus = rebuild_overlap(ncur, eleIdx);
  Slater[parameter] = original - delta;
  minus = rebuild_overlap(ncur, eleIdx);
  Slater[parameter] = original;
  restored = rebuild_overlap(ncur, eleIdx);
  CHECK(finite_close(restored, baseOverlap,
                     2.0e-13 * (1.0 + cabs(baseOverlap))),
        "FD baseline restore ncur=%d parameter=%d got=(%.17g,%.17g) "
        "expected=(%.17g,%.17g)",
        ncur, parameter, creal(restored), cimag(restored), creal(baseOverlap),
        cimag(baseOverlap));
  return (plus - minus) / (2.0 * epsilon * baseOverlap);
}

static void check_derivative(const int ncur, const int *eleIdx,
                             const int realParameter,
                             const int imaginaryParameter) {
  double complex *derivative =
      calloc((size_t)(2 * NSlater), sizeof(*derivative));
  double complex baseOverlap = rebuild_overlap(ncur, eleIdx);
  double complex expectedReal;
  double complex expectedImag;
  if (derivative == NULL) {
    CHECK(0, "derivative allocation failed");
    return;
  }
  memset(derivative, 0x5a, (size_t)(2 * NSlater) * sizeof(*derivative));
  SlaterElmDiffGC_fcmp(derivative, baseOverlap, eleIdx, ncur);
  expectedReal = finite_difference(ncur, eleIdx, realParameter, 0, 1.0e-5,
                                   baseOverlap);
  expectedImag = finite_difference(ncur, eleIdx, imaginaryParameter, 1, 5.0e-6,
                                   baseOverlap);
  CHECK(finite_close(derivative[2 * realParameter], expectedReal,
                     2.0e-9 * (1.0 + cabs(expectedReal))),
        "real FD ncur=%d parameter=%d got=(%.17g,%.17g) expected=(%.17g,%.17g)",
        ncur, realParameter, creal(derivative[2 * realParameter]),
        cimag(derivative[2 * realParameter]), creal(expectedReal),
        cimag(expectedReal));
  CHECK(finite_close(derivative[2 * imaginaryParameter + 1], expectedImag,
                     2.0e-9 * (1.0 + cabs(expectedImag))),
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
  free(derivative);
}

static void initialize_fixture(void) {
  int row;
  int parameter = 0;
  NThread = omp_get_max_threads();
  /* The general-pair fixture must not rely on the global's zero default. */
  iFlgOrbitalGeneral = 1;
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
  CHECK(GCComponentInit() == GC_MALL_OK, "component init");
}

static void initialize_multiqp_fixture(void) {
  initialize_fixture();
  NQPFull = 2;
  NQPFix = 2;
  NMPTrans = 2;
  QPTrans = qpTransRows2;
  QPTransSgn = qpTransSgnRows2;
  SlaterElm = realloc(
      SlaterElm, 2U * ORBITALS * ORBITALS * sizeof(*SlaterElm));
  InvM = realloc(InvM, 2U * ORBITALS * ORBITALS * sizeof(*InvM));
  PfM = realloc(PfM, 2U * sizeof(*PfM));
  QPFullWeight = realloc(QPFullWeight, 2U * sizeof(*QPFullWeight));
  CHECK(SlaterElm != NULL && InvM != NULL && PfM != NULL &&
            QPFullWeight != NULL,
        "multiqp fixture allocation");
  if (failures != 0) exit(EXIT_FAILURE);
  QPFullWeight[0] = 0.73 - 0.21 * I;
  QPFullWeight[1] = 0.37 + 0.52 * I;
  Slater[orbitalIdxStorage[0][1]] = 0.0;
  Slater[orbitalIdxStorage[0][2]] = 0.0;
  Slater[orbitalIdxStorage[0][3]] = 0.0;
  rebuild_matrix = UpdateSlaterElm_fsz;
  CHECK(GCComponentInit() == GC_MALL_OK, "multiqp component init");
}

static void check_component_states(const int ncur, const int *eleIdx,
                                   const int state0, const int state1,
                                   const char *label) {
  rebuild_matrix();
  CHECK(CalculateMAllGC_fcmp(ncur, eleIdx, 0, NQPFull) == GC_MALL_OK,
        "%s rebuild", label);
  CHECK(GCQpState[0] == state0 && GCQpState[1] == state1,
        "%s states %s/%s", label, GCComponentStateName(GCQpState[0]),
        GCComponentStateName(GCQpState[1]));
}

static void check_lowrank_component_gradient(void) {
  const double complex u[4] = {0.7 + 0.1 * I, -0.3 + 0.4 * I,
                               0.5 - 0.2 * I, 0.9 + 0.3 * I};
  const double complex v[4] = {0.2 - 0.6 * I, 0.8 + 0.1 * I,
                               -0.4 + 0.5 * I, 0.1 + 0.7 * I};
  const int eleIdx[4] = {0, 1, 2, 3};
  double complex gradient[ORBITALS * ORBITALS];
  double complex plus;
  double complex minus;
  double complex finiteDifference;
  const double epsilon = 1.0e-5;
  int i;
  NQPFull = 1;
  NQPFix = 1;
  NMPTrans = 1;
  QPTrans = qpTransRows;
  QPTransSgn = qpTransSgnRows;
  rebuild_matrix = refresh_slater_elements;
  for (i = 0; i < ORBITALS; i++) {
    int j;
    for (j = i + 1; j < ORBITALS; j++) {
      Slater[orbitalIdxStorage[i][j]] =
          0.5 * (u[i] * v[j] - v[i] * u[j]);
    }
  }
  CHECK(GCComponentInit() == GC_MALL_OK, "low-rank component init");
  rebuild_matrix();
  CHECK(CalculateMAllGC_fcmp(4, eleIdx, 0, 1) == GC_MALL_OK,
        "low-rank rebuild");
  CHECK(GCQpState[0] == GC_QP_PFONLY, "low-rank state %s",
        GCComponentStateName(GCQpState[0]));
  CHECK(GCComponentGradient(0, 0, 4, eleIdx, gradient,
                            GCComponentThreadCache()) == GC_MALL_OK,
        "gradient status");
  Slater[orbitalIdxStorage[0][1]] += epsilon;
  rebuild_matrix();
  CHECK(CalculateMAllGC_fcmp(4, eleIdx, 0, 1) == GC_MALL_OK, "plus");
  plus = PfM[0];
  Slater[orbitalIdxStorage[0][1]] -= 2.0 * epsilon;
  rebuild_matrix();
  CHECK(CalculateMAllGC_fcmp(4, eleIdx, 0, 1) == GC_MALL_OK, "minus");
  minus = PfM[0];
  Slater[orbitalIdxStorage[0][1]] += epsilon;
  finiteDifference = (plus - minus) / (2.0 * epsilon);
  CHECK(cabs(gradient[1]) > 1.0e-3,
        "low-rank cofactor is vacuous");
  CHECK(finite_close(2.0 * gradient[1], finiteDifference,
                     2.0e-9 * (1.0 + cabs(finiteDifference))),
        "low-rank dPf/df got=(%.17g,%.17g) fd=(%.17g,%.17g)",
        creal(2.0 * gradient[1]), cimag(2.0 * gradient[1]),
        creal(finiteDifference), cimag(finiteDifference));
  CHECK(gradient[NsizeMax] == -gradient[1],
        "cofactor antisymmetry");
  CHECK(gradient[0] == 0.0, "cofactor diagonal");
}

static void run_multiqp_checks(void) {
  const int occupied02[2] = {0, 2};
  const int occupied13[2] = {1, 3};
  const int occupied23[2] = {2, 3};
  initialize_multiqp_fixture();
  check_component_states(2, occupied02, GC_QP_PFONLY, GC_QP_INVERSE,
                         "{0,2}");
  check_component_states(2, occupied13, GC_QP_INVERSE, GC_QP_PFONLY,
                         "{1,3}");
  check_component_states(2, occupied23, GC_QP_INVERSE, GC_QP_INVERSE,
                         "{2,3}");
  check_derivative(2, occupied02, orbitalIdxStorage[0][2],
                   orbitalIdxStorage[1][3]);
  check_derivative(2, occupied13, orbitalIdxStorage[1][3],
                   orbitalIdxStorage[0][2]);
  check_derivative(2, occupied23, orbitalIdxStorage[2][3],
                   orbitalIdxStorage[1][2]);
  check_lowrank_component_gradient();
  free_fixture();
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

/* ------------------------------------------------------------------------
 * Anti-parallel (OrbitalAntiParallel) GC fixtures.  OrbitalIdx/OrbitalSgn are
 * Nsite x Nsite and SlaterElm is built by the production
 * UpdateSlaterElm_fcmp() with NSPGaussLeg=1.  Sizes are allocated here and
 * never shared with the static OrbitalGeneral fixture above.
 * ---------------------------------------------------------------------- */
static int **antiIdxRows = NULL;
static int **antiSgnRows = NULL;
static int *antiQPTrans = NULL;
static int *antiQPTransSgn = NULL;
static int *antiQPOpt = NULL;
static int *antiQPOptSgn = NULL;
static int *antiQPTransRows[1];
static int *antiQPTransSgnRows[1];
static int *antiQPOptRows[1];
static int *antiQPOptSgnRows[1];

static int **allocate_rows(const int count) {
  int **rows = calloc((size_t)count, sizeof(*rows));
  int row;
  if (rows == NULL) return NULL;
  for (row = 0; row < count; row++) {
    rows[row] = calloc((size_t)count, sizeof(**rows));
    if (rows[row] == NULL) return NULL;
  }
  return rows;
}

static void free_rows(int **rows, const int count) {
  int row;
  if (rows == NULL) return;
  for (row = 0; row < count; row++) free(rows[row]);
  free(rows);
}

static int anti_wraps(const int i, const int j, const int nsite) {
  return (i == 0 && j == nsite - 1) || (i == nsite - 1 && j == 0);
}

static void set_antiparallel_translation(const int permuted) {
  int site;
  for (site = 0; site < Nsite; site++) {
    antiQPTrans[site] = permuted ? (site + 1) % Nsite : site;
    antiQPTransSgn[site] = (permuted && site == Nsite - 1) ? -1 : 1;
    antiQPOpt[site] = site;
    antiQPOptSgn[site] = 1;
  }
}

/* shared: parameters indexed by (i+j)%2 for two sites and by (4i+j)%5 for
 * four sites (two classes would make the 4x4 F singular), with the (1,1)
 * entry flipped so a class enters with both signs.  ap: wrapping pairs
 * (0,L-1) and (L-1,0) carry the sign -1. */
static int anti_shared_index(const int i, const int j, const int nsite) {
  return nsite == 2 ? (i + j) % 2 : (nsite * i + j) % 5;
}

static void initialize_antiparallel_fixture(const int nsite, const int shared,
                                            const int ap) {
  int i;
  int j;
  int parameter;
  NThread = omp_get_max_threads();
  iFlgOrbitalGeneral = 0;
  Nsite = nsite;
  Nsite2 = 2 * nsite;
  NsizeMax = Nsite2;
  NQPFull = 1;
  NQPFix = 1;
  NMPTrans = 1;
  NSPGaussLeg = 1;
  NQPOptTrans = 1;
  NSlater = shared ? (nsite == 2 ? 2 : 5) : nsite * nsite;
  LapackLWork = 1024;
  antiIdxRows = allocate_rows(nsite);
  antiSgnRows = allocate_rows(nsite);
  antiQPTrans = calloc((size_t)nsite, sizeof(int));
  antiQPTransSgn = calloc((size_t)nsite, sizeof(int));
  antiQPOpt = calloc((size_t)nsite, sizeof(int));
  antiQPOptSgn = calloc((size_t)nsite, sizeof(int));
  Slater = malloc((size_t)NSlater * sizeof(*Slater));
  SlaterElm = malloc((size_t)Nsite2 * (size_t)Nsite2 * sizeof(*SlaterElm));
  InvM = malloc((size_t)NsizeMax * (size_t)NsizeMax * sizeof(*InvM));
  PfM = malloc(sizeof(*PfM));
  QPFullWeight = malloc(sizeof(*QPFullWeight));
  SPGLCosSin = malloc(sizeof(*SPGLCosSin));
  SPGLCosCos = malloc(sizeof(*SPGLCosCos));
  SPGLSinSin = malloc(sizeof(*SPGLSinSin));
  if (antiIdxRows == NULL || antiSgnRows == NULL || antiQPTrans == NULL ||
      antiQPTransSgn == NULL || antiQPOpt == NULL || antiQPOptSgn == NULL ||
      Slater == NULL || SlaterElm == NULL || InvM == NULL || PfM == NULL ||
      QPFullWeight == NULL || SPGLCosSin == NULL || SPGLCosCos == NULL ||
      SPGLSinSin == NULL) {
    fprintf(stderr, "anti-parallel fixture allocation failed\n");
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < nsite; i++) {
    for (j = 0; j < nsite; j++) {
      int sign = (ap && anti_wraps(i, j, nsite)) ? -1 : 1;
      if (shared && i == 1 && j == 1) sign = -sign;
      antiIdxRows[i][j] =
          shared ? anti_shared_index(i, j, nsite) : i * nsite + j;
      antiSgnRows[i][j] = sign;
    }
  }
  for (parameter = 0; parameter < NSlater; parameter++) {
    Slater[parameter] = (0.25 + 0.4 * cos(0.7 * parameter + 0.3)) +
                        0.3 * sin(1.3 * parameter + 0.2) * I;
  }
  OrbitalIdx = antiIdxRows;
  OrbitalSgn = antiSgnRows;
  antiQPTransRows[0] = antiQPTrans;
  antiQPTransSgnRows[0] = antiQPTransSgn;
  antiQPOptRows[0] = antiQPOpt;
  antiQPOptSgnRows[0] = antiQPOptSgn;
  QPTrans = antiQPTransRows;
  QPTransSgn = antiQPTransSgnRows;
  QPOptTrans = antiQPOptRows;
  QPOptTransSgn = antiQPOptSgnRows;
  set_antiparallel_translation(0);
  QPFullWeight[0] = 0.73 - 0.21 * I;
  SPGLCosSin[0] = 0.0;
  SPGLCosCos[0] = 1.0;
  SPGLSinSin[0] = 0.0;
  rebuild_matrix = UpdateSlaterElm_fcmp;
  initializeWorkSpaceAll();
  CHECK(GCComponentInit() == GC_MALL_OK,
        "anti-parallel component init");
}

static void free_antiparallel_fixture(void) {
  GCComponentRelease();
  FreeWorkSpaceAll();
  free(SPGLSinSin);
  free(SPGLCosCos);
  free(SPGLCosSin);
  free(QPFullWeight);
  free(PfM);
  free(InvM);
  free(SlaterElm);
  free(Slater);
  free(antiQPOptSgn);
  free(antiQPOpt);
  free(antiQPTransSgn);
  free(antiQPTrans);
  free_rows(antiSgnRows, Nsite);
  free_rows(antiIdxRows, Nsite);
  antiIdxRows = antiSgnRows = NULL;
  SPGLCosSin = SPGLCosCos = SPGLSinSin = NULL;
  rebuild_matrix = refresh_slater_elements;
}

static void run_cofactor_omp_audit(void) {
#ifdef MVMC_ENABLE_FAULT_INJECTION
  const char *requiredText = getenv("MVMC_GC_REQUIRE_OMP_THREADS");
  int eleIdx[8];
  double complex gradient[64];
  int iteration;
  int i;
  if (requiredText == NULL || *requiredText == '\0') return;
  initialize_antiparallel_fixture(4, 0, 0);
  memset(SlaterElm, 0, 64U * sizeof(*SlaterElm));
  SlaterElm[2U * 8U + 3U] = 1.0;
  SlaterElm[3U * 8U + 2U] = -1.0;
  SlaterElm[4U * 8U + 5U] = 1.0;
  SlaterElm[5U * 8U + 4U] = -1.0;
  SlaterElm[6U * 8U + 7U] = 1.0;
  SlaterElm[7U * 8U + 6U] = -1.0;
  for (i = 0; i < 8; i++) eleIdx[i] = i;
  CHECK(GCComponentRebuild(0, 0, 8, eleIdx,
                           GCComponentThreadCache()) == GC_MALL_OK,
        "cofactor audit rebuild");
  CHECK(GCQpState[0] == GC_QP_PFONLY, "cofactor audit state %s",
        GCComponentStateName(GCQpState[0]));
  for (iteration = 0; iteration < 64; iteration++) {
    CHECK(GCComponentGradient(0, 0, 8, eleIdx, gradient,
                              GCComponentThreadCache()) == GC_MALL_OK,
          "cofactor audit gradient iteration=%d", iteration);
    for (i = 0; i < 8; i++) {
      int j;
      for (j = 0; j < 8; j++) {
        const double complex expected =
            i == 0 && j == 1 ? 1.0 : i == 1 && j == 0 ? -1.0 : 0.0;
        CHECK(finite_close(gradient[(size_t)i * 8U + (size_t)j], expected,
                           1.0e-12),
              "cofactor audit gradient[%d,%d]", i, j);
      }
    }
  }
  {
    const int required = atoi(requiredText);
    GCComponentOmpAudit audit;
    GCComponentOmpAuditTotal(GC_OMP_COFACTOR, &audit);
    CHECK(audit.teamMin == required && audit.teamMax == required,
          "cofactor audit team min/max=%d/%d expected=%d", audit.teamMin,
          audit.teamMax, required);
    CHECK(audit.workers >= 2, "cofactor audit workers=%d", audit.workers);
    CHECK(audit.workCount == 64LL * 28LL,
          "cofactor audit work=%lld expected=%lld", audit.workCount,
          64LL * 28LL);
  }
  free_antiparallel_fixture();
#endif
}

static void check_all_derivatives(const int ncur, const int *eleIdx) {
  int parameter;
  for (parameter = 0; parameter < NSlater; parameter++) {
    check_derivative(ncur, eleIdx, parameter, parameter);
  }
}

/* d/df_shared must equal sum_ij sgn_ij d/dF_ij of the one-parameter-per-pair
 * expansion that represents the identical matrix. */
static void check_shared_sum(const int ncur, const int *eleIdx) {
  const int nsite = Nsite;
  const int sharedCount = NSlater;
  double complex *sharedValues = malloc((size_t)sharedCount * sizeof(*sharedValues));
  double complex *sharedDerivative =
      calloc((size_t)(2 * sharedCount), sizeof(*sharedDerivative));
  double complex *expandedDerivative =
      calloc((size_t)(2 * nsite * nsite), sizeof(*expandedDerivative));
  double complex *expandedSlater =
      malloc((size_t)nsite * (size_t)nsite * sizeof(*expandedSlater));
  int **sharedIdx = OrbitalIdx;
  int **sharedSgn = OrbitalSgn;
  int **expandedIdx = allocate_rows(nsite);
  int **expandedSgn = allocate_rows(nsite);
  double complex *sharedSlater = Slater;
  double complex sharedPf;
  double complex expandedPf;
  int i;
  int j;
  int parameter;
  if (sharedValues == NULL || sharedDerivative == NULL ||
      expandedDerivative == NULL || expandedSlater == NULL ||
      expandedIdx == NULL || expandedSgn == NULL) {
    fprintf(stderr, "shared-sum allocation failed\n");
    exit(EXIT_FAILURE);
  }
  memcpy(sharedValues, Slater, (size_t)sharedCount * sizeof(*sharedValues));
  sharedPf = rebuild_overlap(ncur, eleIdx);
  SlaterElmDiffGC_fcmp(sharedDerivative, sharedPf, eleIdx, ncur);
  for (i = 0; i < nsite; i++) {
    for (j = 0; j < nsite; j++) {
      expandedIdx[i][j] = i * nsite + j;
      expandedSgn[i][j] = 1;
      expandedSlater[i * nsite + j] =
          (double)sharedSgn[i][j] * sharedValues[sharedIdx[i][j]];
    }
  }
  OrbitalIdx = expandedIdx;
  OrbitalSgn = expandedSgn;
  Slater = expandedSlater;
  NSlater = nsite * nsite;
  expandedPf = rebuild_overlap(ncur, eleIdx);
  SlaterElmDiffGC_fcmp(expandedDerivative, expandedPf, eleIdx, ncur);
  CHECK(finite_close(expandedPf, sharedPf, 2.0e-12 * (1.0 + cabs(sharedPf))),
        "shared/expanded overlap differ ncur=%d", ncur);
  for (parameter = 0; parameter < sharedCount; parameter++) {
    int component;
    for (component = 0; component < 2; component++) {
      double complex sum = 0.0;
      for (i = 0; i < nsite; i++) {
        for (j = 0; j < nsite; j++) {
          if (sharedIdx[i][j] == parameter) {
            sum += (double)sharedSgn[i][j] *
                   expandedDerivative[2 * (i * nsite + j) + component];
          }
        }
      }
      CHECK(finite_close(sharedDerivative[2 * parameter + component], sum,
                         2.0e-12 * (1.0 + cabs(sum))),
            "shared parameter %d component %d is not the sum of its "
            "contributions (ncur=%d)",
            parameter, component, ncur);
    }
  }
  OrbitalIdx = sharedIdx;
  OrbitalSgn = sharedSgn;
  Slater = sharedSlater;
  NSlater = sharedCount;
  free_rows(expandedSgn, nsite);
  free_rows(expandedIdx, nsite);
  free(expandedSlater);
  free(expandedDerivative);
  free(sharedDerivative);
  free(sharedValues);
}

/* The same state written as OrbitalGeneral: upper-triangle parameter F/2 on
 * up-down pairs and 0 on same-spin pairs, built by UpdateSlaterElm_fsz().
 * Pfaffians agree, and d/dp = 2 d/dF by the chain rule (identity
 * translation, one parameter per anti-parallel pair). */
static void check_general_equivalence(const int ncur, const int *eleIdx) {
  const int nsite = Nsite;
  const int orbitals = 2 * nsite;
  const int generalCount = nsite * (2 * nsite - 1);
  int **antiIdx = OrbitalIdx;
  int **antiSgn = OrbitalSgn;
  double complex *antiSlater = Slater;
  const int antiCount = NSlater;
  int **generalIdx = allocate_rows(orbitals);
  int **generalSgn = allocate_rows(orbitals);
  double complex *generalSlater =
      calloc((size_t)generalCount, sizeof(*generalSlater));
  double complex *antiDerivative =
      calloc((size_t)(2 * antiCount), sizeof(*antiDerivative));
  double complex *generalDerivative =
      calloc((size_t)(2 * generalCount), sizeof(*generalDerivative));
  int pairIndex[2 * 8][2 * 8];
  double complex antiPf;
  double complex generalPf;
  int first;
  int parameter = 0;
  if (generalIdx == NULL || generalSgn == NULL || generalSlater == NULL ||
      antiDerivative == NULL || generalDerivative == NULL || orbitals > 16) {
    fprintf(stderr, "general-equivalence allocation failed\n");
    exit(EXIT_FAILURE);
  }
  antiPf = rebuild_overlap(ncur, eleIdx);
  SlaterElmDiffGC_fcmp(antiDerivative, antiPf, eleIdx, ncur);
  for (first = 0; first < orbitals; first++) {
    int second;
    for (second = first + 1; second < orbitals; second++) {
      generalIdx[first][second] = parameter;
      generalIdx[second][first] = parameter;
      generalSgn[first][second] = 1;
      generalSgn[second][first] = -1;
      pairIndex[first][second] = parameter;
      if (first < nsite && second >= nsite) {
        const int i = first;
        const int j = second - nsite;
        generalSlater[parameter] =
            0.5 * (double)antiSgn[i][j] * antiSlater[antiIdx[i][j]];
      }
      parameter++;
    }
  }
  OrbitalIdx = generalIdx;
  OrbitalSgn = generalSgn;
  Slater = generalSlater;
  NSlater = generalCount;
  iFlgOrbitalGeneral = 1;
  rebuild_matrix = UpdateSlaterElm_fsz;
  generalPf = rebuild_overlap(ncur, eleIdx);
  SlaterElmDiffGC_fcmp(generalDerivative, generalPf, eleIdx, ncur);
  CHECK(finite_close(generalPf, antiPf, 2.0e-12 * (1.0 + cabs(antiPf))),
        "General F/2 overlap differs ncur=%d anti=(%.17g,%.17g) "
        "general=(%.17g,%.17g)",
        ncur, creal(antiPf), cimag(antiPf), creal(generalPf),
        cimag(generalPf));
  {
    int i;
    for (i = 0; i < nsite; i++) {
      int j;
      for (j = 0; j < nsite; j++) {
        const int k = antiIdx[i][j];
        const int p = pairIndex[i][j + nsite];
        int component;
        for (component = 0; component < 2; component++) {
          const double complex expected =
              0.5 * (double)antiSgn[i][j] *
              generalDerivative[2 * p + component];
          CHECK(finite_close(antiDerivative[2 * k + component], expected,
                             2.0e-10 * (1.0 + cabs(expected))),
                "General chain rule pair (%d,%d) component %d ncur=%d", i, j,
                component, ncur);
        }
      }
    }
  }
  OrbitalIdx = antiIdx;
  OrbitalSgn = antiSgn;
  Slater = antiSlater;
  NSlater = antiCount;
  iFlgOrbitalGeneral = 0;
  rebuild_matrix = UpdateSlaterElm_fcmp;
  free(generalDerivative);
  free(antiDerivative);
  free(generalSlater);
  free_rows(generalSgn, orbitals);
  free_rows(generalIdx, orbitals);
}

static void run_antiparallel_checks(void) {
  /* fused index = site + spin*Nsite; electron order is arbitrary. */
  const int occupied2[] = {3, 0};
  const int occupied4[] = {2, 0, 3, 1};
  const int four2[] = {5, 2};
  const int four4[] = {6, 1, 0, 5};
  const int four6[] = {7, 3, 4, 0, 6, 2};
  const int four8[] = {4, 1, 7, 2, 0, 5, 3, 6};
  int shared;
  int ap;

  for (shared = 0; shared < 2; shared++) {
    for (ap = 0; ap < 2; ap++) {
      initialize_antiparallel_fixture(2, shared, ap);
      check_all_derivatives(0, NULL);
      check_all_derivatives(2, occupied2);
      check_all_derivatives(4, occupied4);
      if (shared) {
        check_shared_sum(2, occupied2);
        check_shared_sum(4, occupied4);
      } else {
        check_general_equivalence(2, occupied2);
        check_general_equivalence(4, occupied4);
      }
      set_antiparallel_translation(1);
      check_all_derivatives(2, occupied2);
      check_all_derivatives(4, occupied4);
      free_antiparallel_fixture();

      initialize_antiparallel_fixture(4, shared, ap);
      check_all_derivatives(0, NULL);
      check_all_derivatives(2, four2);
      check_all_derivatives(4, four4);
      check_all_derivatives(6, four6);
      check_all_derivatives(8, four8);
      if (shared) {
        check_shared_sum(4, four4);
        check_shared_sum(8, four8);
      } else {
        check_general_equivalence(2, four2);
        check_general_equivalence(4, four4);
        check_general_equivalence(6, four6);
        check_general_equivalence(8, four8);
      }
      set_antiparallel_translation(1);
      check_all_derivatives(4, four4);
      check_all_derivatives(8, four8);
      free_antiparallel_fixture();
    }
  }
}

/* ------------------------------------------------------------------------
 * Spin-projected anti-parallel fixtures (NSPGaussLeg=G, single translation).
 * The G=1 fixture above is widened to G spin-quadrature components built by
 * the production UpdateSlaterElm_fcmp() rotation, so d ln Psi/d f has to
 * include every rotated spin block, not only the up-down one.
 * ---------------------------------------------------------------------- */
static void extend_antiparallel_spin(const int gauss, const double *beta,
                                     const double *weight) {
  const size_t slaterCount =
      (size_t)gauss * (size_t)Nsite2 * (size_t)Nsite2;
  const size_t inverseCount =
      (size_t)gauss * (size_t)NsizeMax * (size_t)NsizeMax;
  int q;
  NSPGaussLeg = gauss;
  NQPFix = gauss;
  NQPFull = gauss;
  SlaterElm = realloc(SlaterElm, slaterCount * sizeof(*SlaterElm));
  InvM = realloc(InvM, inverseCount * sizeof(*InvM));
  PfM = realloc(PfM, (size_t)gauss * sizeof(*PfM));
  QPFullWeight = realloc(QPFullWeight, (size_t)gauss * sizeof(*QPFullWeight));
  SPGLCosSin = realloc(SPGLCosSin, (size_t)gauss * sizeof(*SPGLCosSin));
  SPGLCosCos = realloc(SPGLCosCos, (size_t)gauss * sizeof(*SPGLCosCos));
  SPGLSinSin = realloc(SPGLSinSin, (size_t)gauss * sizeof(*SPGLSinSin));
  if (SlaterElm == NULL || InvM == NULL || PfM == NULL ||
      QPFullWeight == NULL || SPGLCosSin == NULL || SPGLCosCos == NULL ||
      SPGLSinSin == NULL) {
    fprintf(stderr, "spin fixture allocation failed\n");
    exit(EXIT_FAILURE);
  }
  for (q = 0; q < gauss; q++) {
    const double c = cos(0.5 * beta[q]);
    const double s = sin(0.5 * beta[q]);
    SPGLCosSin[q] = c * s;
    SPGLCosCos[q] = c * c;
    SPGLSinSin[q] = s * s;
    QPFullWeight[q] = weight[q];
  }
  CHECK(GCComponentInit() == GC_MALL_OK, "spin component init");
}

/* The rotated same-spin block must be present, otherwise the derivative
 * checks below could pass with the G=1 (up-down only) chain rule. */
static void check_rotation_is_present(const int q) {
  UpdateSlaterElm_fcmp();
  CHECK(cabs(SlaterElm[(size_t)q * (size_t)Nsite2 * (size_t)Nsite2 + 1U]) >
            1.0e-3,
        "spin component %d has no up-up block", q);
}

static void run_spin_checks(void) {
  const double beta[2] = {0.37, 1.21};
  const double weight[2] = {0.7, 0.3};
  const int two2[] = {3, 0};
  const int two4[] = {2, 0, 3, 1};
  const int four2[] = {5, 2};
  const int four4[] = {6, 1, 0, 5};
  const int four6[] = {7, 3, 4, 0, 6, 2};
  const int four8[] = {4, 1, 7, 2, 0, 5, 3, 6};
  int shared;
  int ap;

  for (shared = 0; shared < 2; shared++) {
    for (ap = 0; ap < 2; ap++) {
      initialize_antiparallel_fixture(2, shared, ap);
      extend_antiparallel_spin(2, beta, weight);
      /* The shared two-site table has F01 = F10 (a pure singlet), whose
       * same-spin blocks vanish for every beta. */
      if (!shared) {
        check_rotation_is_present(0);
        check_rotation_is_present(1);
      }
      check_all_derivatives(0, NULL);
      check_all_derivatives(2, two2);
      check_all_derivatives(4, two4);
      set_antiparallel_translation(1);
      check_all_derivatives(2, two2);
      check_all_derivatives(4, two4);
      free_antiparallel_fixture();

      initialize_antiparallel_fixture(4, shared, ap);
      extend_antiparallel_spin(2, beta, weight);
      check_rotation_is_present(1);
      check_all_derivatives(0, NULL);
      check_all_derivatives(2, four2);
      check_all_derivatives(4, four4);
      check_all_derivatives(6, four6);
      check_all_derivatives(8, four8);
      set_antiparallel_translation(1);
      check_all_derivatives(4, four4);
      check_all_derivatives(8, four8);
      free_antiparallel_fixture();
    }
  }
}

/* One spin component is exactly zero (PFONLY) while the projected overlap is
 * not: its cofactor must still enter the derivative. */
static void run_spin_pfonly_checks(void) {
  const double beta[2] = {0.0, 0.5 * acos(-1.0)};
  const double weight[2] = {0.7, 0.3};
  const int occupied03[2] = {0, 3};
  double complex gradient[16];
  double complex ip;
  initialize_antiparallel_fixture(2, 0, 0);
  extend_antiparallel_spin(2, beta, weight);
  Slater[0] = 1.0;
  Slater[1] = 0.0;
  Slater[2] = 2.0 + 0.2 * I;
  Slater[3] = 3.0;
  ip = rebuild_overlap(2, occupied03);
  CHECK(GCQpState[0] == GC_QP_PFONLY, "spin zero component state %s",
        GCComponentStateName(GCQpState[0]));
  CHECK(GCQpState[1] == GC_QP_INVERSE, "spin regular component state %s",
        GCComponentStateName(GCQpState[1]));
  CHECK(PfM[0] == 0.0, "spin zero component Pf=(%.17g,%.17g)",
        creal(PfM[0]), cimag(PfM[0]));
  CHECK(cabs(ip) > 0.1, "spin composite is nonzero");
  CHECK(GCComponentGradient(0, 0, 2, occupied03, gradient,
                            GCComponentThreadCache()) == GC_MALL_OK,
        "spin zero component gradient");
  /* beta=0: A_{0up,1dn} = cc F01 with cc=1, so dPf/dF01 of the zero
   * component is the cofactor 1. */
  CHECK(finite_close(gradient[1], 1.0, 1.0e-14),
        "spin zero component dPf/dF01=(%.17g,%.17g)", creal(gradient[1]),
        cimag(gradient[1]));
  check_all_derivatives(2, occupied03);
  CHECK(GCQpState[0] == GC_QP_PFONLY,
        "spin restored zero component state %s",
        GCComponentStateName(GCQpState[0]));
  CHECK(GCQpState[1] == GC_QP_INVERSE,
        "spin restored regular component state %s",
        GCComponentStateName(GCQpState[1]));
  CHECK(PfM[0] == 0.0, "spin restored zero component Pf=(%.17g,%.17g)",
        creal(PfM[0]), cimag(PfM[0]));
  free_antiparallel_fixture();
}

/* Rotated PFONLY cofactors under OpenMP.  F = diag(0,1,2,3) has no
 * antisymmetric part and cc+ss=1, so every rotated component equals
 * [[0,D],[-D,0]] with D = F: Pf = det D = 0 at full filling and the only
 * nonzero cofactor is C[0,4] = 1*2*3 = 6.  The 28 pairs of each call are
 * distributed over the team. */
static void run_spin_cofactor_omp_audit(void) {
#ifdef MVMC_ENABLE_FAULT_INJECTION
  const char *requiredText = getenv("MVMC_GC_REQUIRE_OMP_THREADS");
  const double beta[2] = {0.37, 1.21};
  const double weight[2] = {0.7, 0.3};
  int eleIdx[8];
  double complex gradient[64];
  int iteration;
  int i;
  int q;
  if (requiredText == NULL || *requiredText == '\0') {
    CHECK(0, "--spin-omp requires MVMC_GC_REQUIRE_OMP_THREADS");
    return;
  }
  initialize_antiparallel_fixture(4, 0, 0);
  extend_antiparallel_spin(2, beta, weight);
  for (i = 0; i < NSlater; i++) Slater[i] = 0.0;
  for (i = 0; i < 4; i++) Slater[i * 4 + i] = (double)i;
  UpdateSlaterElm_fcmp();
  for (i = 0; i < 8; i++) eleIdx[i] = i;
  for (q = 0; q < 2; q++) {
    CHECK(GCComponentRebuild(q, q, 8, eleIdx, GCComponentThreadCache()) ==
              GC_MALL_OK,
          "spin cofactor audit rebuild q=%d", q);
    CHECK(GCQpState[q] == GC_QP_PFONLY,
          "rotated component %d is not PFONLY (%s)", q,
          GCComponentStateName(GCQpState[q]));
  }
  for (iteration = 0; iteration < 64; iteration++) {
    q = iteration % 2;
    CHECK(GCComponentGradient(q, q, 8, eleIdx, gradient,
                              GCComponentThreadCache()) == GC_MALL_OK,
          "spin cofactor audit gradient iteration=%d", iteration);
    for (i = 0; i < 8; i++) {
      int j;
      for (j = 0; j < 8; j++) {
        const double complex expected =
            i == 0 && j == 4 ? 6.0 : i == 4 && j == 0 ? -6.0 : 0.0;
        CHECK(finite_close(gradient[(size_t)i * 8U + (size_t)j], expected,
                           1.0e-12),
              "spin cofactor audit gradient[%d,%d] q=%d", i, j, q);
      }
    }
  }
  {
    const int required = atoi(requiredText);
    GCComponentOmpAudit audit;
    GCComponentOmpAuditTotal(GC_OMP_COFACTOR, &audit);
    CHECK(audit.teamMin == required && audit.teamMax == required,
          "spin cofactor audit team min/max=%d/%d expected=%d",
          audit.teamMin, audit.teamMax, required);
    CHECK(audit.workers >= 2, "spin cofactor audit workers=%d",
          audit.workers);
    CHECK(audit.workCount == 64LL * 28LL,
          "spin cofactor audit work=%lld expected=%lld", audit.workCount,
          64LL * 28LL);
    printf("spin cofactor audit: team=%d/%d workers=%d work=%lld\n",
           audit.teamMin, audit.teamMax, audit.workers, audit.workCount);
  }
  free_antiparallel_fixture();
#else
  CHECK(0, "--spin-omp requires MVMC_ENABLE_FAULT_INJECTION");
#endif
}

static void free_fixture(void) {
  GCComponentRelease();
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
#endif
  if (argc > 1 && strcmp(argv[1], "--antiparallel") == 0) {
    run_antiparallel_checks();
  } else if (argc > 1 && strcmp(argv[1], "--spin") == 0) {
    run_spin_checks();
  } else if (argc > 1 && strcmp(argv[1], "--spin-pfonly") == 0) {
    run_spin_pfonly_checks();
  } else if (argc > 1 && strcmp(argv[1], "--spin-omp") == 0) {
    run_spin_cofactor_omp_audit();
  } else if (argc > 1 && strcmp(argv[1], "--multiqp") == 0) {
    run_multiqp_checks();
    run_cofactor_omp_audit();
  } else {
    initialize_fixture();
    check_derivative(0, NULL, 0, 1);
    check_derivative(2, eleIdx2, 0, 0);
    check_derivative(4, eleIdx4, 1, 4);
    run_signed_checks();
    free_fixture();
  }
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
