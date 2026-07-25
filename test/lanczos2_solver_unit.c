#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "physcal_lanczos2.h"

int CalculateEne(double H1, double H2_1, double H2_2, double H3, double H4,
                 double *alpha_p, double *ene_p, double *ene_vp,
                 double *alpha_m, double *ene_m, double *ene_vm);

static int failures = 0;

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "Lanczos2Solver_Unit FAIL: ");                           \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, "\n");                                                    \
      failures++;                                                               \
    }                                                                           \
  } while (0)

static int NearlyEqual(double actual, double expected, double relative) {
  const double scale = fmax(fmax(fabs(actual), fabs(expected)), 1.0);
  return fabs(actual - expected) <= relative * scale;
}

static int NearlyEqualRelative(double actual, double expected,
                               double relative) {
  const double scale = fmax(fabs(expected), DBL_MIN);
  return fabs(actual - expected) <= relative * scale;
}

static int TokenCount(char *line) {
  int count = 0;
  char *token = strtok(line, " \t\r\n");
  while (token != NULL) {
    count++;
    token = strtok(NULL, " \t\r\n");
  }
  return count;
}

static void FillHankelMoment(const double *eigenvalue, const double *weight,
                             int count, double *moment) {
  double powerMoment[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int i;
  int power;
  int a;
  int b;
  for (i = 0; i < count; i++) {
    double value = 1.0;
    for (power = 0; power <= 6; power++) {
      powerMoment[power] += weight[i] * value;
      value *= eigenvalue[i];
    }
  }
  for (a = 0; a < LANCZOS2_POWER_COUNT; a++) {
    for (b = 0; b < LANCZOS2_POWER_COUNT; b++) {
      moment[4 * a + b] = powerMoment[a + b];
    }
  }
}

static void TestMatrixMapping(void) {
  double moment[16] = {0.0};
  double overlap[4];
  double hamiltonian[4];
  double hamiltonianSquared[4];
  double antiHermitianResidual;
  double hankelResidual;
  Lanczos2SolveStatus status;

  moment[0] = 1.0;
  moment[1] = moment[4] = 0.2;
  moment[2] = moment[8] = 0.9;
  moment[5] = 1.3;
  moment[6] = moment[9] = -0.4;
  moment[10] = 2.2;

  status = BuildLanczos2MatricesReal(
      moment, 2, overlap, hamiltonian, hamiltonianSquared,
      &antiHermitianResidual, &hankelResidual);
  CHECK(status == LANCZOS2_SOLVE_OK, "matrix mapping status=%d", status);
  CHECK(NearlyEqual(overlap[0], 1.0, 1.0e-15), "S00=%g", overlap[0]);
  CHECK(NearlyEqual(overlap[1], 0.2, 1.0e-15), "S01=%g", overlap[1]);
  CHECK(NearlyEqual(overlap[3], 1.3, 1.0e-15), "S11=%g", overlap[3]);
  CHECK(NearlyEqual(hamiltonian[0], 0.2, 1.0e-15), "H00=%g",
        hamiltonian[0]);
  CHECK(NearlyEqual(hamiltonian[1], 1.1, 1.0e-15), "H01=%g",
        hamiltonian[1]);
  CHECK(NearlyEqual(hamiltonian[3], -0.4, 1.0e-15), "H11=%g",
        hamiltonian[3]);
  CHECK(NearlyEqual(hamiltonianSquared[0], 1.3, 1.0e-15), "G00=%g",
        hamiltonianSquared[0]);
  CHECK(NearlyEqual(hamiltonianSquared[1], -0.4, 1.0e-15), "G01=%g",
        hamiltonianSquared[1]);
  CHECK(NearlyEqual(hamiltonianSquared[3], 2.2, 1.0e-15), "G11=%g",
        hamiltonianSquared[3]);
  CHECK(antiHermitianResidual == 0.0, "anti-Hermitian residual=%g",
        antiHermitianResidual);
  CHECK(hankelResidual > 0.0,
        "non-Hankel fixture unexpectedly has zero residual");
}

static void TestHankelResidualUsesAntidiagonalScale(void) {
  double moment[16] = {0.0};
  double overlap[4];
  double hamiltonian[4];
  double hamiltonianSquared[4];
  double antiHermitianResidual;
  double hankelResidual;
  Lanczos2SolveStatus status;

  moment[0] = 1.0;
  moment[1] = moment[4] = 0.2;
  moment[2] = moment[8] = 1.0;
  moment[5] = 1.1;
  moment[6] = moment[9] = 2.0;
  moment[10] = 1.0;
  moment[15] = 1.0e12;
  status = BuildLanczos2MatricesReal(
      moment, 2, overlap, hamiltonian, hamiltonianSquared,
      &antiHermitianResidual, &hankelResidual);
  CHECK(status == LANCZOS2_SOLVE_OK,
        "antidiagonal-scale fixture status=%d", status);
  CHECK(hankelResidual > 0.08,
        "low-order Hankel mismatch hidden by high-order scale: %.17g",
        hankelResidual);
}

static void TestExactThreeStateRealAndComplex(void) {
  const double eigenvalue[3] = {-2.0, 1.0, 4.0};
  const double weight[3] = {0.2, 0.3, 0.5};
  double moment[16];
  double complex complexMoment[16];
  Lanczos2Result realResult;
  Lanczos2Result complexResult;
  Lanczos2SolveStatus status;
  int i;

  FillHankelMoment(eigenvalue, weight, 3, moment);
  for (i = 0; i < 16; i++) complexMoment[i] = moment[i];

  status = SolveLanczos2Real(moment, 3, &realResult);
  CHECK(status == LANCZOS2_SOLVE_OK, "real exact solve status=%d (%s)",
        status, Lanczos2SolveError(status));
  CHECK(realResult.solveFlags == 0, "real solve flags=%d",
        realResult.solveFlags);
  CHECK(NearlyEqual(realResult.shift, 1.9, 1.0e-14), "real shift=%.17g",
        realResult.shift);
  CHECK(NearlyEqual(realResult.energy, -2.0, 1.0e-13),
        "real energy=%.17g", realResult.energy);
  CHECK(fabs(realResult.variance) <= 2.0e-13,
        "real variance=%.17g", realResult.variance);
  CHECK(realResult.antiHermitianResidual == 0.0,
        "real anti-Hermitian residual=%g",
        realResult.antiHermitianResidual);
  CHECK(realResult.hankelResidual == 0.0, "real Hankel residual=%g",
        realResult.hankelResidual);
  CHECK(realResult.eigenvalueResidual <= 2.0e-13,
        "real eigenvalue residual=%g", realResult.eigenvalueResidual);
  CHECK(NearlyEqual(creal(realResult.alpha[0]), -1.25, 2.0e-13),
        "real alpha1=%.17g", creal(realResult.alpha[0]));
  CHECK(NearlyEqual(creal(realResult.alpha[1]), 0.25, 2.0e-13),
        "real alpha2=%.17g", creal(realResult.alpha[1]));

  status = SolveLanczos2Complex(complexMoment, 3, &complexResult);
  CHECK(status == LANCZOS2_SOLVE_OK, "complex exact solve status=%d (%s)",
        status, Lanczos2SolveError(status));
  CHECK(NearlyEqual(complexResult.energy, -2.0, 1.0e-13),
        "complex energy=%.17g", complexResult.energy);
  CHECK(fabs(complexResult.variance) <= 2.0e-13,
        "complex variance=%.17g", complexResult.variance);
  CHECK(NearlyEqual(creal(complexResult.alpha[0]), -1.25, 2.0e-13),
        "complex alpha1=%.17g", creal(complexResult.alpha[0]));
  CHECK(NearlyEqual(creal(complexResult.alpha[1]), 0.25, 2.0e-13),
        "complex alpha2=%.17g", creal(complexResult.alpha[1]));
  for (i = 0; i < 3; i++) {
    CHECK(fabs(cimag(complexResult.coefficient[i])) <= 2.0e-14,
          "complex coefficient[%d] imaginary part=%.17g", i,
          cimag(complexResult.coefficient[i]));
  }
  {
    int maximumIndex = 0;
    for (i = 1; i < 3; i++) {
      if (cabs(complexResult.coefficient[i]) >
          cabs(complexResult.coefficient[maximumIndex])) {
        maximumIndex = i;
      }
    }
    CHECK(creal(complexResult.coefficient[maximumIndex]) > 0.0,
          "phase-fixed maximum coefficient is not positive");
    CHECK(fabs(cimag(complexResult.coefficient[maximumIndex])) <= 2.0e-14,
          "phase-fixed maximum coefficient is not real");
  }
}

static void TestLargeEnergyShiftAccuracy(void) {
  const double offset[6] = {0.0, 0.7, 1.6, 2.9, 4.1, 6.3};
  const double weight[6] = {0.30, 0.25, 0.18, 0.12, 0.09, 0.06};
  const unsigned char valid[6] = {1, 1, 1, 1, 1, 1};
  const double energyShift = 100.0;
  double shiftedEigenvalue[6];
  double samplePower[6 * LANCZOS2_POWER_COUNT];
  double centeredMoment[16];
  double shiftedMoment[16];
  Lanczos2Result centered;
  Lanczos2Result shifted;
  Lanczos2SolveStatus centeredStatus;
  Lanczos2SolveStatus shiftedStatus;
  double expectedEnergy;
  double expectedVarianceRatio;
  int i;
  int k;

  for (i = 0; i < 6; i++) {
    shiftedEigenvalue[i] = offset[i] - energyShift;
    for (k = 0; k < LANCZOS2_POWER_COUNT; k++) {
      samplePower[i * LANCZOS2_POWER_COUNT + k] =
          pow(shiftedEigenvalue[i], (double)k);
    }
  }
  FillHankelMoment(offset, weight, 6, centeredMoment);
  centeredStatus = SolveLanczos2Real(centeredMoment, 3, &centered);
  shiftedStatus = BuildLanczos2ShiftedMomentReal(
      shiftedMoment, samplePower, weight, valid, 6, -energyShift);
  CHECK(shiftedStatus == LANCZOS2_SOLVE_OK,
        "large-shift moment construction status=%d (%s)", shiftedStatus,
        Lanczos2SolveError(shiftedStatus));
  if (shiftedStatus == LANCZOS2_SOLVE_OK) {
    shiftedStatus = SolveLanczos2ShiftedReal(
        shiftedMoment, 3, -energyShift, &shifted);
  }
  CHECK(centeredStatus == LANCZOS2_SOLVE_OK,
        "large-shift centered solve status=%d (%s)", centeredStatus,
        Lanczos2SolveError(centeredStatus));
  CHECK(shiftedStatus == LANCZOS2_SOLVE_OK,
        "large-shift shifted solve status=%d (%s)", shiftedStatus,
        Lanczos2SolveError(shiftedStatus));
  if (centeredStatus != LANCZOS2_SOLVE_OK ||
      shiftedStatus != LANCZOS2_SOLVE_OK) {
    return;
  }
  expectedEnergy = centered.energy - energyShift;
  expectedVarianceRatio =
      centered.variance / (expectedEnergy * expectedEnergy);
  CHECK(NearlyEqualRelative(shifted.energy, expectedEnergy, 1.0e-8),
        "large-shift energy actual=%.17g expected=%.17g",
        shifted.energy, expectedEnergy);
  CHECK(NearlyEqualRelative(shifted.varianceOverEnergySquared,
                            expectedVarianceRatio, 1.0e-8),
        "large-shift variance/E2 actual=%.17g expected=%.17g",
        shifted.varianceOverEnergySquared, expectedVarianceRatio);
}

static void TestRegularizedRetry(void) {
  const double eigenvalue[2] = {-2.0, 4.0};
  const double weight[2] = {0.5, 0.5};
  double moment[16];
  Lanczos2Result result;
  const Lanczos2SolveStatus status = (
      FillHankelMoment(eigenvalue, weight, 2, moment),
      SolveLanczos2Real(moment, 3, &result));
  CHECK(status == LANCZOS2_SOLVE_OK, "regularized solve status=%d (%s)",
        status, Lanczos2SolveError(status));
  CHECK((result.solveFlags & LANCZOS2_FLAG_REGULARIZED) != 0,
        "regularized flag was not set (flags=%d)", result.solveFlags);
  CHECK(NearlyEqual(result.epsilon, 1.0e-12, 1.0e-14),
        "epsilon=%.17g", result.epsilon);
  CHECK(NearlyEqual(result.energy, -2.0, 2.0e-11),
        "regularized energy=%.17g", result.energy);
}

static void TestNonPositiveNormalizationFails(void) {
  double moment[16] = {0.0};
  Lanczos2Result result;
  Lanczos2SolveStatus status;
  const double rho = 1.0 + 5.0e-13;

  moment[0] = 1.0;
  moment[1] = moment[4] = 0.0;
  moment[2] = moment[8] = rho;
  moment[5] = 1.0;
  moment[6] = moment[9] = 0.0;
  moment[10] = 1.0;
  moment[3] = moment[12] = 2.0;
  moment[7] = moment[13] = -1.0;
  moment[11] = moment[14] = 0.0;
  moment[15] = 1.0;

  status = SolveLanczos2Real(moment, 3, &result);
  CHECK(status == LANCZOS2_SOLVE_NORMALIZATION_FAILURE,
        "non-positive normalization status=%d (%s)", status,
        Lanczos2SolveError(status));
  CHECK((result.solveFlags & LANCZOS2_FLAG_REGULARIZED) != 0,
        "non-positive normalization did not exercise retry");
}

static void TestAlphaUndefined(void) {
  double moment[16] = {0.0};
  Lanczos2Result result;
  Lanczos2SolveStatus status;

  moment[0] = 1.0;
  moment[1] = moment[4] = 0.0;
  moment[2] = moment[8] = -1.0;
  moment[5] = 1.0;
  moment[6] = moment[9] = -2.0;
  moment[10] = 5.0;

  status = SolveLanczos2Real(moment, 2, &result);
  CHECK(status == LANCZOS2_SOLVE_OK, "alpha-undefined status=%d (%s)",
        status, Lanczos2SolveError(status));
  CHECK((result.solveFlags & LANCZOS2_FLAG_ALPHA_UNDEFINED) != 0,
        "alpha-undefined flag not set (flags=%d)", result.solveFlags);
  CHECK(isnan(creal(result.alpha[0])), "undefined alpha is not NaN");
  CHECK(NearlyEqual(result.energy, -2.0, 1.0e-14),
        "alpha-undefined energy=%.17g", result.energy);
  CHECK(cabs(result.coefficient[0]) <= 1.0e-14,
        "alpha-undefined c0 magnitude=%.17g",
        cabs(result.coefficient[0]));
}

static void TestNonfiniteAndInvalidInputs(void) {
  double moment[16] = {0.0};
  Lanczos2Result result;
  Lanczos2SolveStatus status;
  moment[0] = 1.0;
  moment[6] = NAN;
  status = SolveLanczos2Real(moment, 3, &result);
  CHECK(status == LANCZOS2_SOLVE_NONFINITE_MOMENT,
        "NaN status=%d", status);

  memset(moment, 0, sizeof(moment));
  moment[0] = 1.0;
  moment[5] = 0.0;
  status = SolveLanczos2Real(moment, 2, &result);
  CHECK(status == LANCZOS2_SOLVE_NONPOSITIVE_DIAGONAL,
        "non-positive shifted diagonal status=%d", status);

  status = SolveLanczos2Real(moment, 1, &result);
  CHECK(status == LANCZOS2_SOLVE_INVALID_ARGUMENT,
        "invalid dimension status=%d", status);
  status = SolveLanczos2Real(NULL, 3, &result);
  CHECK(status == LANCZOS2_SOLVE_INVALID_ARGUMENT,
        "NULL moment status=%d", status);
}

static void TestLegacyDimensionTwoReduction(void) {
  double moment[16] = {0.0};
  Lanczos2Result result;
  double alphaPlus;
  double energyPlus;
  double variancePlus;
  double alphaMinus;
  double energyMinus;
  double varianceMinus;
  double legacyAlpha;
  double legacyEnergy;
  double legacyVariance;
  int legacyStatus;
  Lanczos2SolveStatus status;

  moment[0] = 1.0;
  moment[1] = moment[4] = 0.2;   /* QQQQ[2]: H1 */
  moment[5] = 1.3;               /* QQQQ[3]: H2_1 */
  moment[2] = moment[8] = 0.9;   /* QQQQ[10]: H2_2 */
  moment[6] = moment[9] = -0.4;  /* QQQQ[11]: H3 */
  moment[10] = 2.2;              /* QQQQ[15]: H4 */

  legacyStatus =
      CalculateEne(0.2, 1.3, 0.9, -0.4, 2.2, &alphaPlus, &energyPlus,
                   &variancePlus, &alphaMinus, &energyMinus, &varianceMinus);
  CHECK(legacyStatus == 0, "legacy CalculateEne status=%d", legacyStatus);
  if (energyPlus <= energyMinus) {
    legacyAlpha = alphaPlus;
    legacyEnergy = energyPlus;
    legacyVariance = variancePlus;
  } else {
    legacyAlpha = alphaMinus;
    legacyEnergy = energyMinus;
    legacyVariance = varianceMinus;
  }

  status = SolveLanczos2Real(moment, 2, &result);
  CHECK(status == LANCZOS2_SOLVE_OK, "d=2 solve status=%d (%s)", status,
        Lanczos2SolveError(status));
  CHECK(NearlyEqual(result.energy, legacyEnergy, 1.0e-13),
        "d=2 energy new=%.17g legacy=%.17g", result.energy, legacyEnergy);
  CHECK(NearlyEqual(result.varianceOverEnergySquared, legacyVariance,
                    1.0e-13),
        "d=2 variance/E2 new=%.17g legacy=%.17g",
        result.varianceOverEnergySquared, legacyVariance);
  CHECK(NearlyEqual(creal(result.alpha[0]), legacyAlpha, 1.0e-13),
        "d=2 alpha new=%.17g legacy=%.17g", creal(result.alpha[0]),
        legacyAlpha);
}

static void TestMomentAccumulation(void) {
  const double realPower[4] = {1.0, -2.0, 3.0, -4.0};
  const double complex complexPower[4] = {
      1.0, 1.0 + 2.0 * I, -3.0 + 0.5 * I, 2.0 - I};
  double realMoment[16] = {0.0};
  double complex complexMoment[16] = {0.0};
  Lanczos2SolveStatus status;

  status = AccumulateLanczos2MomentReal(realMoment, realPower, 0.25);
  CHECK(status == LANCZOS2_SOLVE_OK, "real accumulation status=%d", status);
  CHECK(NearlyEqual(realMoment[1 * 4 + 3], 2.0, 1.0e-15),
        "real M13=%.17g", realMoment[1 * 4 + 3]);

  status =
      AccumulateLanczos2MomentComplex(complexMoment, complexPower, 0.5);
  CHECK(status == LANCZOS2_SOLVE_OK, "complex accumulation status=%d",
        status);
  CHECK(cabs(complexMoment[1 * 4 + 2] -
             0.5 * conj(complexPower[1]) * complexPower[2]) <= 1.0e-15,
        "complex M12 has wrong conjugation");
  CHECK(cabs(complexMoment[2 * 4 + 1] -
             conj(complexMoment[1 * 4 + 2])) <= 1.0e-15,
        "complex accumulated moment is not Hermitian");

  {
    double badPower[4] = {1.0, 2.0, NAN, 4.0};
    status = AccumulateLanczos2MomentReal(realMoment, badPower, 1.0);
    CHECK(status == LANCZOS2_SOLVE_NONFINITE_MOMENT,
          "non-finite power accumulation status=%d", status);
  }
  {
    double overflowingPower[4] = {DBL_MAX, 1.0, 1.0, 1.0};
    double overflowMoment[16] = {0.0};
    double complex complexOverflowPower[4] = {
        DBL_MAX + 0.0 * I, 1.0, 1.0, 1.0};
    double complex complexOverflowMoment[16] = {0.0};
    status =
        AccumulateLanczos2MomentReal(overflowMoment, overflowingPower, 2.0);
    CHECK(status == LANCZOS2_SOLVE_NONFINITE_MOMENT,
          "real overflow accumulation status=%d", status);
    CHECK(overflowMoment[0] == 0.0,
          "real overflow partially updated the moment");
    status = AccumulateLanczos2MomentComplex(
        complexOverflowMoment, complexOverflowPower, 2.0);
    CHECK(status == LANCZOS2_SOLVE_NONFINITE_MOMENT,
          "complex overflow accumulation status=%d", status);
    CHECK(complexOverflowMoment[0] == 0.0,
          "complex overflow partially updated the moment");
  }
}

static void TestOutputLayout(void) {
  const double eigenvalue[3] = {-2.0, 1.0, 4.0};
  const double weight[3] = {0.2, 0.3, 0.5};
  double moment[16];
  double complex complexMoment[16];
  Lanczos2Result result;
  FILE *output = tmpfile();
  FILE *momentOutput = tmpfile();
  char header[2048];
  char data[4096];
  char extra[8];
  char tokenBuffer[4096];
  Lanczos2SolveStatus status;
  int i;

  CHECK(output != NULL && momentOutput != NULL, "tmpfile failed");
  if (output == NULL || momentOutput == NULL) {
    if (output != NULL) fclose(output);
    if (momentOutput != NULL) fclose(momentOutput);
    return;
  }
  FillHankelMoment(eigenvalue, weight, 3, moment);
  status = WriteLanczos2OutputReal(
      output, momentOutput, moment, 0.0, &result);
  CHECK(status == LANCZOS2_SOLVE_OK, "real output status=%d (%s)", status,
        Lanczos2SolveError(status));

  rewind(output);
  CHECK(fgets(header, sizeof(header), output) != NULL,
        "ls2_out header missing");
  CHECK(strncmp(header, "# mVMC ls2_out v1 ",
                strlen("# mVMC ls2_out v1 ")) == 0,
        "ls2_out header version mismatch: %s", header);
  CHECK(fgets(data, sizeof(data), output) != NULL, "ls2_out data missing");
  CHECK(data[strlen(data) - 1] == '\n', "ls2_out data lacks final newline");
  memcpy(tokenBuffer, data, strlen(data) + 1);
  CHECK(TokenCount(tokenBuffer) == 17, "ls2_out numeric column count != 17");
  CHECK(fgets(extra, sizeof(extra), output) == NULL,
        "ls2_out has more than two lines");

  rewind(momentOutput);
  CHECK(fgets(header, sizeof(header), momentOutput) != NULL,
        "ls2_moment header missing");
  CHECK(strncmp(header, "# mVMC ls2_moment v2 basis_shift=",
                strlen("# mVMC ls2_moment v2 basis_shift=")) == 0,
        "ls2_moment header version mismatch: %s", header);
  CHECK(fgets(data, sizeof(data), momentOutput) != NULL,
        "ls2_moment data missing");
  CHECK(data[strlen(data) - 1] == '\n',
        "ls2_moment data lacks final newline");
  memcpy(tokenBuffer, data, strlen(data) + 1);
  CHECK(TokenCount(tokenBuffer) == 32,
        "ls2_moment numeric column count != 32");
  CHECK(fgets(extra, sizeof(extra), momentOutput) == NULL,
        "ls2_moment has more than two lines");
  fclose(output);
  fclose(momentOutput);

  output = tmpfile();
  momentOutput = tmpfile();
  CHECK(output != NULL && momentOutput != NULL, "complex tmpfile failed");
  if (output == NULL || momentOutput == NULL) {
    if (output != NULL) fclose(output);
    if (momentOutput != NULL) fclose(momentOutput);
    return;
  }
  for (i = 0; i < 16; i++) complexMoment[i] = moment[i];
  status = WriteLanczos2OutputComplex(
      output, momentOutput, complexMoment, 0.0, &result);
  CHECK(status == LANCZOS2_SOLVE_OK, "complex output status=%d (%s)",
        status, Lanczos2SolveError(status));
  rewind(output);
  CHECK(fgets(header, sizeof(header), output) != NULL,
        "complex ls2_out header missing");
  CHECK(fgets(data, sizeof(data), output) != NULL,
        "complex ls2_out data missing");
  memcpy(tokenBuffer, data, strlen(data) + 1);
  CHECK(TokenCount(tokenBuffer) == 17,
        "complex ls2_out numeric column count != 17");
  rewind(momentOutput);
  CHECK(fgets(header, sizeof(header), momentOutput) != NULL,
        "complex ls2_moment header missing");
  CHECK(fgets(data, sizeof(data), momentOutput) != NULL,
        "complex ls2_moment data missing");
  memcpy(tokenBuffer, data, strlen(data) + 1);
  CHECK(TokenCount(tokenBuffer) == 32,
        "complex ls2_moment numeric column count != 32");
  fclose(output);
  fclose(momentOutput);
}

int main(void) {
  TestMatrixMapping();
  TestHankelResidualUsesAntidiagonalScale();
  TestExactThreeStateRealAndComplex();
  TestLargeEnergyShiftAccuracy();
  TestRegularizedRetry();
  TestNonPositiveNormalizationFails();
  TestAlphaUndefined();
  TestNonfiniteAndInvalidInputs();
  TestLegacyDimensionTwoReduction();
  TestMomentAccumulation();
  TestOutputLayout();

  if (failures != 0) {
    fprintf(stderr, "Lanczos2Solver_Unit: %d failure(s)\n", failures);
    return 1;
  }
  printf("Lanczos2Solver_Unit: PASS\n");
  return 0;
}
