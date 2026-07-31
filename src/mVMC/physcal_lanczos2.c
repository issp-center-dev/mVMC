#include <complex.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#include "blas_externs.h"
#include "physcal_lanczos2.h"

#define MATRIX_INDEX(row, column, leadingDimension) \
  ((row) * (leadingDimension) + (column))

enum {
  SUPPORT_VALID = 0,
  SUPPORT_WEIGHT,
  SUPPORT_M02_REAL,
  SUPPORT_M02_IMAG,
  SUPPORT_M11_REAL,
  SUPPORT_M11_IMAG,
  SUPPORT_M03_REAL,
  SUPPORT_M03_IMAG,
  SUPPORT_M12_REAL,
  SUPPORT_M12_IMAG
};

static const double POWER_LANCZOS_SUPPORT_RELATIVE_TOLERANCE = 1.0e-8;
static const double POWER_LANCZOS_SUPPORT_SCORE_TOLERANCE = 4.5;
/* A separate effect-size gate prevents a noisy, order-one disagreement
 * from being labeled pass without turning the numerical 1e-8 mismatch
 * floor into an impractical equivalence-test tolerance. */
static const double POWER_LANCZOS_SUPPORT_UNRESOLVED_RELATIVE_TOLERANCE =
    0.5;

static int ValidDimension(int dimension) {
  return dimension == 2 || dimension == 3;
}

static void ClearResult(Lanczos2Result *result, int dimension) {
  int i;
  memset(result, 0, sizeof(*result));
  result->dimension = dimension;
  result->energy = NAN;
  result->variance = NAN;
  result->varianceOverEnergySquared = NAN;
  result->eigenvalueResidual = NAN;
  for (i = 0; i < LANCZOS2_MAX_DIMENSION - 1; i++) {
    result->alpha[i] = NAN + I * NAN;
  }
}

static double RelativeResidual(double numerator, double scale) {
  return numerator / fmax(scale, DBL_MIN);
}

static void ClearSupportDiagnostic(
    PowerLanczosSupportDiagnostic *diagnostic, int hasThirdAnchor) {
  memset(diagnostic, 0, sizeof(*diagnostic));
  diagnostic->status = POWER_LANCZOS_SUPPORT_INVALID;
  diagnostic->hasThirdAnchor = hasThirdAnchor != 0;
  diagnostic->m02 = diagnostic->m11 = NAN + I * NAN;
  diagnostic->m03 = diagnostic->m12 = NAN + I * NAN;
  diagnostic->delta2 = diagnostic->delta3 = NAN + I * NAN;
  diagnostic->relativeDifference2 = NAN;
  diagnostic->relativeDifference3 = NAN;
  diagnostic->standardError2 = NAN;
  diagnostic->standardError3 = NAN;
  diagnostic->score2 = NAN;
  diagnostic->score3 = NAN;
}

Lanczos2SolveStatus RecordPowerLanczosSupportSampleComplex(
    double sample[POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH],
    const double complex localPower[LANCZOS2_POWER_COUNT], int powerCount,
    double weight) {
  double complex m02;
  double complex m11;
  double complex m03 = 0.0;
  double complex m12 = 0.0;
  int i;
  if (sample == NULL || localPower == NULL ||
      (powerCount != 3 && powerCount != 4) ||
      !isfinite(weight) || weight < 0.0) {
    return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  }
  /* A zero self-normalized reweighting contribution carries no moment and
   * remains an unmarked row; it may reduce the final valid-sample count. */
  if (weight == 0.0) {
    memset(sample, 0, POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH * sizeof(*sample));
    return LANCZOS2_SOLVE_OK;
  }
  for (i = 0; i < powerCount; i++) {
    if (!isfinite(creal(localPower[i])) ||
        !isfinite(cimag(localPower[i]))) {
      return LANCZOS2_SOLVE_NONFINITE_MOMENT;
    }
  }
  m02 = weight * conj(localPower[0]) * localPower[2];
  m11 = weight * conj(localPower[1]) * localPower[1];
  if (powerCount == 4) {
    m03 = weight * conj(localPower[0]) * localPower[3];
    m12 = weight * conj(localPower[1]) * localPower[2];
  }
  if (!isfinite(creal(m02)) || !isfinite(cimag(m02)) ||
      !isfinite(creal(m11)) || !isfinite(cimag(m11)) ||
      !isfinite(creal(m03)) || !isfinite(cimag(m03)) ||
      !isfinite(creal(m12)) || !isfinite(cimag(m12))) {
    return LANCZOS2_SOLVE_NONFINITE_MOMENT;
  }
  sample[SUPPORT_VALID] = 1.0;
  sample[SUPPORT_WEIGHT] = weight;
  sample[SUPPORT_M02_REAL] = creal(m02);
  sample[SUPPORT_M02_IMAG] = cimag(m02);
  sample[SUPPORT_M11_REAL] = creal(m11);
  sample[SUPPORT_M11_IMAG] = cimag(m11);
  sample[SUPPORT_M03_REAL] = creal(m03);
  sample[SUPPORT_M03_IMAG] = cimag(m03);
  sample[SUPPORT_M12_REAL] = creal(m12);
  sample[SUPPORT_M12_IMAG] = cimag(m12);
  return LANCZOS2_SOLVE_OK;
}

Lanczos2SolveStatus RecordPowerLanczosSupportSampleReal(
    double sample[POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH],
    const double localPower[LANCZOS2_POWER_COUNT], int powerCount,
    double weight) {
  double complex complexPower[LANCZOS2_POWER_COUNT] = {0.0, 0.0, 0.0, 0.0};
  int i;
  if (localPower == NULL) return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  for (i = 0; i < powerCount && i < LANCZOS2_POWER_COUNT; i++) {
    complexPower[i] = localPower[i];
  }
  return RecordPowerLanczosSupportSampleComplex(
      sample, complexPower, powerCount, weight);
}

static double complex SupportValue(const double *sample, int realIndex) {
  return sample[realIndex] + I * sample[realIndex + 1];
}

static int FiniteComplex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static double ComplexBlockStandardError(
    const double complex *value, int count) {
  double complex mean = 0.0;
  double sumSquared = 0.0;
  int block;
  if (count < 2) return NAN;
  for (block = 0; block < count; block++) mean += value[block];
  mean /= (double)count;
  for (block = 0; block < count; block++) {
    const double residual = cabs(value[block] - mean);
    sumSquared += residual * residual;
  }
  return sqrt(sumSquared / ((double)count * (double)(count - 1)));
}

PowerLanczosSupportStatus AnalyzePowerLanczosSupport(
    const double *sampleData, size_t sampleCapacity, int hasThirdAnchor,
    PowerLanczosSupportDiagnostic *diagnostic) {
  double complex totalM02 = 0.0;
  double complex totalM11 = 0.0;
  double complex totalM03 = 0.0;
  double complex totalM12 = 0.0;
  double blockWeight[POWER_LANCZOS_SUPPORT_MAX_BLOCKS] = {0.0};
  double complex blockM02[POWER_LANCZOS_SUPPORT_MAX_BLOCKS] = {0.0};
  double complex blockM11[POWER_LANCZOS_SUPPORT_MAX_BLOCKS] = {0.0};
  double complex blockM03[POWER_LANCZOS_SUPPORT_MAX_BLOCKS] = {0.0};
  double complex blockM12[POWER_LANCZOS_SUPPORT_MAX_BLOCKS] = {0.0};
  double complex blockDelta2[POWER_LANCZOS_SUPPORT_MAX_BLOCKS];
  double complex blockDelta3[POWER_LANCZOS_SUPPORT_MAX_BLOCKS];
  double sumWeight = 0.0;
  double sumWeightSquared = 0.0;
  size_t validCount = 0;
  size_t validIndex = 0;
  size_t sampleIndex;
  int blockCount;
  int block;
  if (diagnostic == NULL) return POWER_LANCZOS_SUPPORT_INVALID;
  ClearSupportDiagnostic(diagnostic, hasThirdAnchor);
  if (sampleData == NULL || sampleCapacity == 0) return diagnostic->status;

  for (sampleIndex = 0; sampleIndex < sampleCapacity; sampleIndex++) {
    const double *sample =
        sampleData + sampleIndex * POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH;
    int field;
    if (sample[SUPPORT_VALID] == 0.0) continue;
    if (sample[SUPPORT_VALID] != 1.0) return diagnostic->status;
    for (field = SUPPORT_WEIGHT;
         field < POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH; field++) {
      if (!isfinite(sample[field])) return diagnostic->status;
    }
    if (!(sample[SUPPORT_WEIGHT] > 0.0)) return diagnostic->status;
    validCount++;
  }
  diagnostic->sampleCount = validCount;
  if (validCount < POWER_LANCZOS_SUPPORT_MIN_SAMPLES) {
    diagnostic->status = POWER_LANCZOS_SUPPORT_INCONCLUSIVE;
    return diagnostic->status;
  }
  blockCount = (int)(validCount / 8);
  if (blockCount > POWER_LANCZOS_SUPPORT_MAX_BLOCKS) {
    blockCount = POWER_LANCZOS_SUPPORT_MAX_BLOCKS;
  }
  if (blockCount < 4) {
    diagnostic->status = POWER_LANCZOS_SUPPORT_INCONCLUSIVE;
    return diagnostic->status;
  }
  diagnostic->blockCount = blockCount;

  for (sampleIndex = 0; sampleIndex < sampleCapacity; sampleIndex++) {
    const double *sample =
        sampleData + sampleIndex * POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH;
    const double weight = sample[SUPPORT_WEIGHT];
    if (sample[SUPPORT_VALID] == 0.0) continue;
    block = (int)(validIndex * (size_t)blockCount / validCount);
    if (block >= blockCount) block = blockCount - 1;
    sumWeight += weight;
    sumWeightSquared += weight * weight;
    totalM02 += SupportValue(sample, SUPPORT_M02_REAL);
    totalM11 += SupportValue(sample, SUPPORT_M11_REAL);
    totalM03 += SupportValue(sample, SUPPORT_M03_REAL);
    totalM12 += SupportValue(sample, SUPPORT_M12_REAL);
    blockWeight[block] += weight;
    blockM02[block] += SupportValue(sample, SUPPORT_M02_REAL);
    blockM11[block] += SupportValue(sample, SUPPORT_M11_REAL);
    blockM03[block] += SupportValue(sample, SUPPORT_M03_REAL);
    blockM12[block] += SupportValue(sample, SUPPORT_M12_REAL);
    validIndex++;
  }
  if (!isfinite(sumWeight) || !isfinite(sumWeightSquared) ||
      !(sumWeight > 0.0) || !(sumWeightSquared > 0.0) ||
      !FiniteComplex(totalM02) || !FiniteComplex(totalM11) ||
      !FiniteComplex(totalM03) || !FiniteComplex(totalM12)) {
    return diagnostic->status;
  }
  diagnostic->totalWeight = sumWeight;
  diagnostic->effectiveSampleCount =
      sumWeight * sumWeight / sumWeightSquared;
  diagnostic->m02 = totalM02 / sumWeight;
  diagnostic->m11 = totalM11 / sumWeight;
  diagnostic->m03 = totalM03 / sumWeight;
  diagnostic->m12 = totalM12 / sumWeight;
  diagnostic->delta2 = diagnostic->m02 - diagnostic->m11;
  diagnostic->delta3 = diagnostic->m03 - diagnostic->m12;
  diagnostic->relativeDifference2 =
      RelativeResidual(cabs(diagnostic->delta2),
                       fmax(cabs(diagnostic->m02), cabs(diagnostic->m11)));
  diagnostic->relativeDifference3 =
      RelativeResidual(cabs(diagnostic->delta3),
                       fmax(cabs(diagnostic->m03), cabs(diagnostic->m12)));
  if (!isfinite(diagnostic->effectiveSampleCount) ||
      !FiniteComplex(diagnostic->m02) || !FiniteComplex(diagnostic->m11) ||
      !FiniteComplex(diagnostic->m03) || !FiniteComplex(diagnostic->m12) ||
      !FiniteComplex(diagnostic->delta2) ||
      !FiniteComplex(diagnostic->delta3) ||
      !isfinite(diagnostic->relativeDifference2) ||
      !isfinite(diagnostic->relativeDifference3)) {
    return diagnostic->status;
  }

  for (block = 0; block < blockCount; block++) {
    if (!(blockWeight[block] > 0.0) ||
        !FiniteComplex(blockM02[block]) ||
        !FiniteComplex(blockM11[block]) ||
        !FiniteComplex(blockM03[block]) ||
        !FiniteComplex(blockM12[block])) {
      return diagnostic->status;
    }
    blockDelta2[block] =
        (blockM02[block] - blockM11[block]) / blockWeight[block];
    blockDelta3[block] =
        (blockM03[block] - blockM12[block]) / blockWeight[block];
  }
  diagnostic->standardError2 =
      ComplexBlockStandardError(blockDelta2, blockCount);
  diagnostic->standardError3 =
      ComplexBlockStandardError(blockDelta3, blockCount);
  if (!isfinite(diagnostic->standardError2) ||
      !isfinite(diagnostic->standardError3)) {
    return diagnostic->status;
  }
  {
    const double scale2 =
        fmax(cabs(diagnostic->m02), cabs(diagnostic->m11));
    const double scale3 =
        fmax(cabs(diagnostic->m03), cabs(diagnostic->m12));
    const double floor2 = 256.0 * DBL_EPSILON * fmax(scale2, 1.0);
    const double floor3 = 256.0 * DBL_EPSILON * fmax(scale3, 1.0);
    diagnostic->score2 = cabs(diagnostic->delta2) /
                         fmax(diagnostic->standardError2, floor2);
    diagnostic->score3 = cabs(diagnostic->delta3) /
                         fmax(diagnostic->standardError3, floor3);
  }
  if (!isfinite(diagnostic->score2) || !isfinite(diagnostic->score3)) {
    return diagnostic->status;
  }
  {
    const int mismatch2 =
        diagnostic->relativeDifference2 >=
            POWER_LANCZOS_SUPPORT_RELATIVE_TOLERANCE &&
        diagnostic->score2 >= POWER_LANCZOS_SUPPORT_SCORE_TOLERANCE;
    const int mismatch3 =
        diagnostic->hasThirdAnchor &&
        diagnostic->relativeDifference3 >=
            POWER_LANCZOS_SUPPORT_RELATIVE_TOLERANCE &&
        diagnostic->score3 >= POWER_LANCZOS_SUPPORT_SCORE_TOLERANCE;
    const int unresolved2 =
        diagnostic->relativeDifference2 >=
            POWER_LANCZOS_SUPPORT_UNRESOLVED_RELATIVE_TOLERANCE &&
        diagnostic->score2 < POWER_LANCZOS_SUPPORT_SCORE_TOLERANCE;
    const int unresolved3 =
        diagnostic->hasThirdAnchor &&
        diagnostic->relativeDifference3 >=
            POWER_LANCZOS_SUPPORT_UNRESOLVED_RELATIVE_TOLERANCE &&
        diagnostic->score3 < POWER_LANCZOS_SUPPORT_SCORE_TOLERANCE;
    diagnostic->status =
        mismatch2 || mismatch3
            ? POWER_LANCZOS_SUPPORT_MISMATCH
            : (unresolved2 || unresolved3
                   ? POWER_LANCZOS_SUPPORT_INCONCLUSIVE
                   : POWER_LANCZOS_SUPPORT_PASS);
  }
  if (!diagnostic->hasThirdAnchor) {
    diagnostic->m03 = diagnostic->m12 = NAN + I * NAN;
    diagnostic->delta3 = NAN + I * NAN;
    diagnostic->relativeDifference3 = NAN;
    diagnostic->standardError3 = NAN;
    diagnostic->score3 = NAN;
  }
  return diagnostic->status;
}

const char *PowerLanczosSupportStatusName(
    PowerLanczosSupportStatus status) {
  switch (status) {
    case POWER_LANCZOS_SUPPORT_PASS:
      return "pass";
    case POWER_LANCZOS_SUPPORT_MISMATCH:
      return "mismatch";
    case POWER_LANCZOS_SUPPORT_INCONCLUSIVE:
      return "inconclusive";
    case POWER_LANCZOS_SUPPORT_INVALID:
      return "invalid";
  }
  return "invalid";
}

Lanczos2SolveStatus WritePowerLanczosSupportDiagnostic(
    FILE *output, int lanczosStep, int experimentalMode,
    const PowerLanczosSupportDiagnostic *diagnostic) {
  const char *quality;
  if (output == NULL || diagnostic == NULL ||
      (lanczosStep != 1 && lanczosStep != 2) ||
      (experimentalMode != 0 && experimentalMode != 1)) {
    return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  }
  quality = diagnostic->status == POWER_LANCZOS_SUPPORT_PASS
                ? "support-check-passed-not-proof"
                : (experimentalMode ? "biased-diagnostic-only"
                                    : "invalid-biased-estimator");
  if (fprintf(output,
              "# mVMC power_lanczos_support v1 step=%d mode=%s "
              "result=%s quality=%s scope=necessary-not-sufficient\n",
              lanczosStep, experimentalMode ? "experimental" : "strict",
              PowerLanczosSupportStatusName(diagnostic->status), quality) < 0 ||
      fprintf(output,
              "# sample_count block_count total_weight effective_sample_count "
              "M02_re M02_im M11_re M11_im delta2_re delta2_im "
              "relative_difference2 standard_error2 score2 "
              "M03_re M03_im M12_re M12_im delta3_re delta3_im "
              "relative_difference3 standard_error3 score3\n") < 0 ||
      fprintf(output,
              "%zu %d % .18e % .18e "
              "% .18e % .18e % .18e % .18e % .18e % .18e "
              "% .18e % .18e % .18e "
              "% .18e % .18e % .18e % .18e % .18e % .18e "
              "% .18e % .18e % .18e\n",
              diagnostic->sampleCount, diagnostic->blockCount,
              diagnostic->totalWeight, diagnostic->effectiveSampleCount,
              creal(diagnostic->m02), cimag(diagnostic->m02),
              creal(diagnostic->m11), cimag(diagnostic->m11),
              creal(diagnostic->delta2), cimag(diagnostic->delta2),
              diagnostic->relativeDifference2,
              diagnostic->standardError2, diagnostic->score2,
              creal(diagnostic->m03), cimag(diagnostic->m03),
              creal(diagnostic->m12), cimag(diagnostic->m12),
              creal(diagnostic->delta3), cimag(diagnostic->delta3),
              diagnostic->relativeDifference3,
              diagnostic->standardError3, diagnostic->score3) < 0 ||
      fflush(output) != 0 || ferror(output)) {
    return LANCZOS2_SOLVE_IO_FAILURE;
  }
  return LANCZOS2_SOLVE_OK;
}

static double HankelResidualReal(const double *moment) {
  double residual = 0.0;
  int order;
  int a;
  int b;
  int c;
  int d;
  for (order = 0; order <= 2 * (LANCZOS2_POWER_COUNT - 1); order++) {
    double maximum = 0.0;
    double scale = 0.0;
    for (a = 0; a < LANCZOS2_POWER_COUNT; a++) {
      b = order - a;
      if (b < 0 || b >= LANCZOS2_POWER_COUNT) continue;
      scale = fmax(scale,
                   sqrt(fabs(moment[MATRIX_INDEX(
                                 a, a, LANCZOS2_POWER_COUNT)]) *
                             fabs(moment[MATRIX_INDEX(
                                 b, b, LANCZOS2_POWER_COUNT)])));
      scale = fmax(
          scale,
          fabs(moment[MATRIX_INDEX(a, b, LANCZOS2_POWER_COUNT)]));
      for (c = 0; c < LANCZOS2_POWER_COUNT; c++) {
        d = order - c;
        if (d >= 0 && d < LANCZOS2_POWER_COUNT) {
          maximum = fmax(maximum,
                         fabs(moment[MATRIX_INDEX(
                                  a, b, LANCZOS2_POWER_COUNT)] -
                              moment[MATRIX_INDEX(
                                  c, d, LANCZOS2_POWER_COUNT)]));
        }
      }
    }
    residual = fmax(residual, RelativeResidual(maximum, scale));
  }
  return residual;
}

static double HankelResidualComplex(const double complex *moment) {
  double residual = 0.0;
  int order;
  int a;
  int b;
  int c;
  int d;
  for (order = 0; order <= 2 * (LANCZOS2_POWER_COUNT - 1); order++) {
    double maximum = 0.0;
    double scale = 0.0;
    for (a = 0; a < LANCZOS2_POWER_COUNT; a++) {
      b = order - a;
      if (b < 0 || b >= LANCZOS2_POWER_COUNT) continue;
      scale = fmax(scale,
                   sqrt(cabs(moment[MATRIX_INDEX(
                                  a, a, LANCZOS2_POWER_COUNT)]) *
                             cabs(moment[MATRIX_INDEX(
                                  b, b, LANCZOS2_POWER_COUNT)])));
      scale = fmax(
          scale,
          cabs(moment[MATRIX_INDEX(a, b, LANCZOS2_POWER_COUNT)]));
      for (c = 0; c < LANCZOS2_POWER_COUNT; c++) {
        d = order - c;
        if (d >= 0 && d < LANCZOS2_POWER_COUNT) {
          maximum = fmax(maximum,
                         cabs(moment[MATRIX_INDEX(
                                  a, b, LANCZOS2_POWER_COUNT)] -
                              moment[MATRIX_INDEX(
                                  c, d, LANCZOS2_POWER_COUNT)]));
        }
      }
    }
    residual = fmax(residual, RelativeResidual(maximum, scale));
  }
  return residual;
}

static Lanczos2SolveStatus ProjectMomentReal(const double *moment,
                                             double *projected,
                                             double *residual) {
  double maximum = 0.0;
  double scale = 0.0;
  int a;
  int b;
  for (a = 0; a < LANCZOS2_POWER_COUNT; a++) {
    for (b = 0; b < LANCZOS2_POWER_COUNT; b++) {
      const double mab =
          moment[MATRIX_INDEX(a, b, LANCZOS2_POWER_COUNT)];
      const double mba =
          moment[MATRIX_INDEX(b, a, LANCZOS2_POWER_COUNT)];
      if (!isfinite(mab)) return LANCZOS2_SOLVE_NONFINITE_MOMENT;
      maximum = fmax(maximum, fabs(mab - mba));
      scale = fmax(scale, fabs(mab));
      projected[MATRIX_INDEX(a, b, LANCZOS2_POWER_COUNT)] =
          0.5 * (mab + mba);
    }
  }
  *residual = RelativeResidual(maximum, scale);
  return LANCZOS2_SOLVE_OK;
}

static Lanczos2SolveStatus ProjectMomentComplex(
    const double complex *moment, double complex *projected,
    double *residual) {
  double maximum = 0.0;
  double scale = 0.0;
  int a;
  int b;
  for (a = 0; a < LANCZOS2_POWER_COUNT; a++) {
    for (b = 0; b < LANCZOS2_POWER_COUNT; b++) {
      const double complex mab =
          moment[MATRIX_INDEX(a, b, LANCZOS2_POWER_COUNT)];
      const double complex mba =
          conj(moment[MATRIX_INDEX(b, a, LANCZOS2_POWER_COUNT)]);
      if (!isfinite(creal(mab)) || !isfinite(cimag(mab))) {
        return LANCZOS2_SOLVE_NONFINITE_MOMENT;
      }
      maximum = fmax(maximum, cabs(mab - mba));
      scale = fmax(scale, cabs(mab));
      projected[MATRIX_INDEX(a, b, LANCZOS2_POWER_COUNT)] =
          0.5 * (mab + mba);
    }
  }
  *residual = RelativeResidual(maximum, scale);
  return LANCZOS2_SOLVE_OK;
}

static void HermitianProjectReal(double *matrix, int dimension) {
  int a;
  int b;
  for (a = 0; a < dimension; a++) {
    for (b = a + 1; b < dimension; b++) {
      const double value =
          0.5 * (matrix[MATRIX_INDEX(a, b, dimension)] +
                 matrix[MATRIX_INDEX(b, a, dimension)]);
      matrix[MATRIX_INDEX(a, b, dimension)] = value;
      matrix[MATRIX_INDEX(b, a, dimension)] = value;
    }
  }
}

static void HermitianProjectComplex(double complex *matrix, int dimension) {
  int a;
  int b;
  for (a = 0; a < dimension; a++) {
    matrix[MATRIX_INDEX(a, a, dimension)] =
        creal(matrix[MATRIX_INDEX(a, a, dimension)]);
    for (b = a + 1; b < dimension; b++) {
      const double complex value =
          0.5 * (matrix[MATRIX_INDEX(a, b, dimension)] +
                 conj(matrix[MATRIX_INDEX(b, a, dimension)]));
      matrix[MATRIX_INDEX(a, b, dimension)] = value;
      matrix[MATRIX_INDEX(b, a, dimension)] = conj(value);
    }
  }
}

static void MatricesFromProjectedMomentReal(const double *moment,
                                            int dimension, double *overlap,
                                            double *hamiltonian,
                                            double *hamiltonianSquared) {
  int a;
  int b;
  for (a = 0; a < dimension; a++) {
    for (b = 0; b < dimension; b++) {
      overlap[MATRIX_INDEX(a, b, dimension)] =
          moment[MATRIX_INDEX(a, b, LANCZOS2_POWER_COUNT)];
      hamiltonian[MATRIX_INDEX(a, b, dimension)] =
          0.5 *
          (moment[MATRIX_INDEX(a, b + 1, LANCZOS2_POWER_COUNT)] +
           moment[MATRIX_INDEX(a + 1, b, LANCZOS2_POWER_COUNT)]);
      hamiltonianSquared[MATRIX_INDEX(a, b, dimension)] =
          moment[MATRIX_INDEX(a + 1, b + 1, LANCZOS2_POWER_COUNT)];
    }
  }
  HermitianProjectReal(overlap, dimension);
  HermitianProjectReal(hamiltonian, dimension);
  HermitianProjectReal(hamiltonianSquared, dimension);
}

static void MatricesFromProjectedMomentComplex(
    const double complex *moment, int dimension, double complex *overlap,
    double complex *hamiltonian, double complex *hamiltonianSquared) {
  int a;
  int b;
  for (a = 0; a < dimension; a++) {
    for (b = 0; b < dimension; b++) {
      overlap[MATRIX_INDEX(a, b, dimension)] =
          moment[MATRIX_INDEX(a, b, LANCZOS2_POWER_COUNT)];
      hamiltonian[MATRIX_INDEX(a, b, dimension)] =
          0.5 *
          (moment[MATRIX_INDEX(a, b + 1, LANCZOS2_POWER_COUNT)] +
           moment[MATRIX_INDEX(a + 1, b, LANCZOS2_POWER_COUNT)]);
      hamiltonianSquared[MATRIX_INDEX(a, b, dimension)] =
          moment[MATRIX_INDEX(a + 1, b + 1, LANCZOS2_POWER_COUNT)];
    }
  }
  HermitianProjectComplex(overlap, dimension);
  HermitianProjectComplex(hamiltonian, dimension);
  HermitianProjectComplex(hamiltonianSquared, dimension);
}

Lanczos2SolveStatus BuildLanczos2MatricesReal(
    const double moment[LANCZOS2_MOMENT_COUNT], int dimension,
    double *overlap, double *hamiltonian, double *hamiltonianSquared,
    double *antiHermitianResidual, double *hankelResidual) {
  double projected[LANCZOS2_MOMENT_COUNT];
  Lanczos2SolveStatus status;
  if (moment == NULL || overlap == NULL || hamiltonian == NULL ||
      hamiltonianSquared == NULL || antiHermitianResidual == NULL ||
      hankelResidual == NULL || !ValidDimension(dimension)) {
    return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  }
  status = ProjectMomentReal(moment, projected, antiHermitianResidual);
  if (status != LANCZOS2_SOLVE_OK) return status;
  *hankelResidual = HankelResidualReal(projected);
  MatricesFromProjectedMomentReal(projected, dimension, overlap, hamiltonian,
                                  hamiltonianSquared);
  return LANCZOS2_SOLVE_OK;
}

Lanczos2SolveStatus BuildLanczos2MatricesComplex(
    const double complex moment[LANCZOS2_MOMENT_COUNT], int dimension,
    double complex *overlap, double complex *hamiltonian,
    double complex *hamiltonianSquared, double *antiHermitianResidual,
    double *hankelResidual) {
  double complex projected[LANCZOS2_MOMENT_COUNT];
  Lanczos2SolveStatus status;
  if (moment == NULL || overlap == NULL || hamiltonian == NULL ||
      hamiltonianSquared == NULL || antiHermitianResidual == NULL ||
      hankelResidual == NULL || !ValidDimension(dimension)) {
    return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  }
  status = ProjectMomentComplex(moment, projected, antiHermitianResidual);
  if (status != LANCZOS2_SOLVE_OK) return status;
  *hankelResidual = HankelResidualComplex(projected);
  MatricesFromProjectedMomentComplex(projected, dimension, overlap,
                                     hamiltonian, hamiltonianSquared);
  return LANCZOS2_SOLVE_OK;
}

static void BuildShiftTransform(double shift, double *transform) {
  /*
   * Rows are the shifted powers (H-shift)^k; columns are powers H^j.
   * Only the 4x4 transform is needed by the second-step moment table.
   */
  static const double binomial[4][4] = {
      {1.0, 0.0, 0.0, 0.0},
      {1.0, 1.0, 0.0, 0.0},
      {1.0, 2.0, 1.0, 0.0},
      {1.0, 3.0, 3.0, 1.0}};
  int k;
  int j;
  memset(transform, 0, LANCZOS2_MOMENT_COUNT * sizeof(*transform));
  for (k = 0; k < LANCZOS2_POWER_COUNT; k++) {
    for (j = 0; j <= k; j++) {
      transform[MATRIX_INDEX(k, j, LANCZOS2_POWER_COUNT)] =
          binomial[k][j] * pow(-shift, (double)(k - j));
    }
  }
}

static void TransformMomentReal(const double *transform, const double *moment,
                                double *shifted) {
  int a;
  int b;
  int j;
  int k;
  for (a = 0; a < LANCZOS2_POWER_COUNT; a++) {
    for (b = 0; b < LANCZOS2_POWER_COUNT; b++) {
      double value = 0.0;
      for (j = 0; j < LANCZOS2_POWER_COUNT; j++) {
        for (k = 0; k < LANCZOS2_POWER_COUNT; k++) {
          value +=
              transform[MATRIX_INDEX(a, j, LANCZOS2_POWER_COUNT)] *
              moment[MATRIX_INDEX(j, k, LANCZOS2_POWER_COUNT)] *
              transform[MATRIX_INDEX(b, k, LANCZOS2_POWER_COUNT)];
        }
      }
      shifted[MATRIX_INDEX(a, b, LANCZOS2_POWER_COUNT)] = value;
    }
  }
}

static void TransformMomentComplex(const double *transform,
                                   const double complex *moment,
                                   double complex *shifted) {
  int a;
  int b;
  int j;
  int k;
  for (a = 0; a < LANCZOS2_POWER_COUNT; a++) {
    for (b = 0; b < LANCZOS2_POWER_COUNT; b++) {
      double complex value = 0.0;
      for (j = 0; j < LANCZOS2_POWER_COUNT; j++) {
        for (k = 0; k < LANCZOS2_POWER_COUNT; k++) {
          value +=
              transform[MATRIX_INDEX(a, j, LANCZOS2_POWER_COUNT)] *
              moment[MATRIX_INDEX(j, k, LANCZOS2_POWER_COUNT)] *
              transform[MATRIX_INDEX(b, k, LANCZOS2_POWER_COUNT)];
        }
      }
      shifted[MATRIX_INDEX(a, b, LANCZOS2_POWER_COUNT)] = value;
    }
  }
}

static Lanczos2SolveStatus TransformLocalPowerReal(
    const double *transform, const double *localPower,
    double *shiftedPower) {
  int k;
  int j;
  for (k = 0; k < LANCZOS2_POWER_COUNT; k++) {
    double value = 0.0;
    for (j = 0; j <= k; j++) {
      value = fma(
          transform[MATRIX_INDEX(k, j, LANCZOS2_POWER_COUNT)],
          localPower[j], value);
    }
    if (!isfinite(value)) return LANCZOS2_SOLVE_NONFINITE_MOMENT;
    shiftedPower[k] = value;
  }
  return LANCZOS2_SOLVE_OK;
}

static Lanczos2SolveStatus TransformLocalPowerComplex(
    const double *transform, const double complex *localPower,
    double complex *shiftedPower) {
  int k;
  int j;
  for (k = 0; k < LANCZOS2_POWER_COUNT; k++) {
    double realValue = 0.0;
    double imaginaryValue = 0.0;
    for (j = 0; j <= k; j++) {
      const double coefficient =
          transform[MATRIX_INDEX(k, j, LANCZOS2_POWER_COUNT)];
      realValue = fma(coefficient, creal(localPower[j]), realValue);
      imaginaryValue =
          fma(coefficient, cimag(localPower[j]), imaginaryValue);
    }
    if (!isfinite(realValue) || !isfinite(imaginaryValue)) {
      return LANCZOS2_SOLVE_NONFINITE_MOMENT;
    }
    shiftedPower[k] = realValue + imaginaryValue * I;
  }
  return LANCZOS2_SOLVE_OK;
}

Lanczos2SolveStatus BuildLanczos2ShiftedMomentReal(
    double moment[LANCZOS2_MOMENT_COUNT], const double *samplePower,
    const double *sampleWeight, const unsigned char *sampleValid,
    size_t sampleCount, double shift) {
  double updatedMoment[LANCZOS2_MOMENT_COUNT] = {0.0};
  double transform[LANCZOS2_MOMENT_COUNT];
  double shiftedPower[LANCZOS2_POWER_COUNT];
  Lanczos2SolveStatus status;
  size_t sample;
  if (moment == NULL || samplePower == NULL || sampleWeight == NULL ||
      sampleValid == NULL || sampleCount == 0) {
    return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  }
  if (!isfinite(shift)) return LANCZOS2_SOLVE_NONFINITE_MOMENT;
  BuildShiftTransform(shift, transform);
  for (sample = 0; sample < sampleCount; sample++) {
    if (sampleValid[sample] == 0) continue;
    status = TransformLocalPowerReal(
        transform, samplePower + sample * LANCZOS2_POWER_COUNT,
        shiftedPower);
    if (status != LANCZOS2_SOLVE_OK) return status;
    status = AccumulateLanczos2MomentReal(
        updatedMoment, shiftedPower, sampleWeight[sample]);
    if (status != LANCZOS2_SOLVE_OK) return status;
  }
  memcpy(moment, updatedMoment, sizeof(updatedMoment));
  return LANCZOS2_SOLVE_OK;
}

Lanczos2SolveStatus BuildLanczos2ShiftedMomentComplex(
    double complex moment[LANCZOS2_MOMENT_COUNT],
    const double complex *samplePower, const double *sampleWeight,
    const unsigned char *sampleValid, size_t sampleCount, double shift) {
  double complex updatedMoment[LANCZOS2_MOMENT_COUNT] = {0.0};
  double transform[LANCZOS2_MOMENT_COUNT];
  double complex shiftedPower[LANCZOS2_POWER_COUNT];
  Lanczos2SolveStatus status;
  size_t sample;
  if (moment == NULL || samplePower == NULL || sampleWeight == NULL ||
      sampleValid == NULL || sampleCount == 0) {
    return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  }
  if (!isfinite(shift)) return LANCZOS2_SOLVE_NONFINITE_MOMENT;
  BuildShiftTransform(shift, transform);
  for (sample = 0; sample < sampleCount; sample++) {
    if (sampleValid[sample] == 0) continue;
    status = TransformLocalPowerComplex(
        transform, samplePower + sample * LANCZOS2_POWER_COUNT,
        shiftedPower);
    if (status != LANCZOS2_SOLVE_OK) return status;
    status = AccumulateLanczos2MomentComplex(
        updatedMoment, shiftedPower, sampleWeight[sample]);
    if (status != LANCZOS2_SOLVE_OK) return status;
  }
  memcpy(moment, updatedMoment, sizeof(updatedMoment));
  return LANCZOS2_SOLVE_OK;
}

static void RowMajorToColumnMajorReal(const double *source,
                                      double *destination, int dimension) {
  int a;
  int b;
  for (a = 0; a < dimension; a++) {
    for (b = 0; b < dimension; b++) {
      destination[a + b * dimension] =
          source[MATRIX_INDEX(a, b, dimension)];
    }
  }
}

static void RowMajorToColumnMajorComplex(const double complex *source,
                                         double complex *destination,
                                         int dimension) {
  int a;
  int b;
  for (a = 0; a < dimension; a++) {
    for (b = 0; b < dimension; b++) {
      destination[a + b * dimension] =
          source[MATRIX_INDEX(a, b, dimension)];
    }
  }
}

static int GeneralizedEigenReal(const double *hamiltonian,
                                const double *overlap, int dimension,
                                double epsilon, double *eigenvalue,
                                double *eigenvector) {
  double a[LANCZOS2_MAX_DIMENSION * LANCZOS2_MAX_DIMENSION];
  double b[LANCZOS2_MAX_DIMENSION * LANCZOS2_MAX_DIMENSION];
  double eigenvalues[LANCZOS2_MAX_DIMENSION];
  double work[64];
  const int itype = 1;
  const int leadingDimension = dimension;
  const int workSize = (int)(sizeof(work) / sizeof(work[0]));
  const char jobz = 'V';
  const char uplo = 'U';
  int i;
  int info = 0;
  RowMajorToColumnMajorReal(hamiltonian, a, dimension);
  RowMajorToColumnMajorReal(overlap, b, dimension);
  for (i = 0; i < dimension; i++) b[i + i * dimension] += epsilon;
  M_DSYGV(&itype, &jobz, &uplo, &dimension, a, &leadingDimension, b,
          &leadingDimension, eigenvalues, work, &workSize, &info);
  if (info == 0) {
    *eigenvalue = eigenvalues[0];
    for (i = 0; i < dimension; i++) eigenvector[i] = a[i];
  }
  return info;
}

static int GeneralizedEigenComplex(const double complex *hamiltonian,
                                   const double complex *overlap,
                                   int dimension, double epsilon,
                                   double *eigenvalue,
                                   double complex *eigenvector) {
  double complex a[LANCZOS2_MAX_DIMENSION * LANCZOS2_MAX_DIMENSION];
  double complex b[LANCZOS2_MAX_DIMENSION * LANCZOS2_MAX_DIMENSION];
  double eigenvalues[LANCZOS2_MAX_DIMENSION];
  double complex work[64];
  double realWork[3 * LANCZOS2_MAX_DIMENSION - 2];
  const int itype = 1;
  const int leadingDimension = dimension;
  const int workSize = (int)(sizeof(work) / sizeof(work[0]));
  const char jobz = 'V';
  const char uplo = 'U';
  int i;
  int info = 0;
  RowMajorToColumnMajorComplex(hamiltonian, a, dimension);
  RowMajorToColumnMajorComplex(overlap, b, dimension);
  for (i = 0; i < dimension; i++) b[i + i * dimension] += epsilon;
  M_ZHEGV(&itype, &jobz, &uplo, &dimension, a, &leadingDimension, b,
          &leadingDimension, eigenvalues, work, &workSize, realWork, &info);
  if (info == 0) {
    *eigenvalue = eigenvalues[0];
    for (i = 0; i < dimension; i++) eigenvector[i] = a[i];
  }
  return info;
}

static double QuadraticReal(const double *coefficient, const double *matrix,
                            int dimension) {
  double value = 0.0;
  int a;
  int b;
  for (a = 0; a < dimension; a++) {
    for (b = 0; b < dimension; b++) {
      value += coefficient[a] *
               matrix[MATRIX_INDEX(a, b, dimension)] * coefficient[b];
    }
  }
  return value;
}

static double complex QuadraticComplex(const double complex *coefficient,
                                       const double complex *matrix,
                                       int dimension) {
  double complex value = 0.0;
  int a;
  int b;
  for (a = 0; a < dimension; a++) {
    for (b = 0; b < dimension; b++) {
      value += conj(coefficient[a]) *
               matrix[MATRIX_INDEX(a, b, dimension)] * coefficient[b];
    }
  }
  return value;
}

static void FinalizeVariance(Lanczos2Result *result, double h2) {
  const double energySquared = result->energy * result->energy;
  const double scale = fmax(fmax(fabs(h2), energySquared), 1.0);
  result->variance = h2 - energySquared;
  if (result->variance < 0.0 &&
      fabs(result->variance) <= 128.0 * DBL_EPSILON * scale) {
    result->variance = 0.0;
  }
  if (energySquared > DBL_MIN) {
    result->varianceOverEnergySquared =
        result->variance / energySquared;
  } else {
    result->varianceOverEnergySquared = NAN;
  }
}

static Lanczos2SolveStatus FinalizeReal(
    const double *overlap, const double *hamiltonian,
    const double *hamiltonianSquared, const double *transform,
    const double *diagonalScale, const double *eigenvector,
    double shiftedEigenvalue, Lanczos2Result *result) {
  double shiftedCoefficient[LANCZOS2_MAX_DIMENSION] = {0.0, 0.0, 0.0};
  double coefficient[LANCZOS2_MAX_DIMENSION] = {0.0, 0.0, 0.0};
  double normalization;
  double maximum = 0.0;
  double h2;
  int maximumIndex = 0;
  int a;
  int j;
  for (a = 0; a < result->dimension; a++) {
    shiftedCoefficient[a] = diagonalScale[a] * eigenvector[a];
  }
  for (j = 0; j < result->dimension; j++) {
    for (a = j; a < result->dimension; a++) {
      coefficient[j] +=
          transform[MATRIX_INDEX(a, j, LANCZOS2_POWER_COUNT)] *
          shiftedCoefficient[a];
    }
  }
  normalization =
      QuadraticReal(coefficient, overlap, result->dimension);
  if (!isfinite(normalization) || !(normalization > 0.0)) {
    return LANCZOS2_SOLVE_NORMALIZATION_FAILURE;
  }
  normalization = sqrt(normalization);
  for (a = 0; a < result->dimension; a++) {
    coefficient[a] /= normalization;
    if (fabs(coefficient[a]) > maximum) {
      maximum = fabs(coefficient[a]);
      maximumIndex = a;
    }
  }
  if (coefficient[maximumIndex] < 0.0) {
    for (a = 0; a < result->dimension; a++) coefficient[a] = -coefficient[a];
  }
  for (a = 0; a < result->dimension; a++) {
    result->coefficient[a] = coefficient[a];
  }
  result->energy =
      QuadraticReal(coefficient, hamiltonian, result->dimension);
  h2 = QuadraticReal(coefficient, hamiltonianSquared, result->dimension);
  if (!isfinite(result->energy) || !isfinite(h2) ||
      !isfinite(shiftedEigenvalue)) {
    return LANCZOS2_SOLVE_NONFINITE_MOMENT;
  }
  FinalizeVariance(result, h2);
  if (!isfinite(result->variance)) {
    return LANCZOS2_SOLVE_NONFINITE_MOMENT;
  }
  result->eigenvalueResidual =
      RelativeResidual(fabs((shiftedEigenvalue + result->shift) -
                            result->energy),
                       fmax(fabs(result->energy), 1.0));
  if (fabs(coefficient[0]) >
      64.0 * DBL_EPSILON * fmax(maximum, DBL_MIN)) {
    for (a = 1; a < result->dimension; a++) {
      result->alpha[a - 1] = coefficient[a] / coefficient[0];
    }
  } else {
    result->solveFlags |= LANCZOS2_FLAG_ALPHA_UNDEFINED;
  }
  return LANCZOS2_SOLVE_OK;
}

static Lanczos2SolveStatus FinalizeComplex(
    const double complex *overlap, const double complex *hamiltonian,
    const double complex *hamiltonianSquared, const double *transform,
    const double *diagonalScale, const double complex *eigenvector,
    double shiftedEigenvalue, Lanczos2Result *result) {
  double complex shiftedCoefficient[LANCZOS2_MAX_DIMENSION] = {0.0, 0.0,
                                                               0.0};
  double complex coefficient[LANCZOS2_MAX_DIMENSION] = {0.0, 0.0, 0.0};
  double complex normalizationComplex;
  double normalization;
  double maximum = 0.0;
  double h2;
  double complex phase;
  int maximumIndex = 0;
  int a;
  int j;
  for (a = 0; a < result->dimension; a++) {
    shiftedCoefficient[a] = diagonalScale[a] * eigenvector[a];
  }
  for (j = 0; j < result->dimension; j++) {
    for (a = j; a < result->dimension; a++) {
      coefficient[j] +=
          transform[MATRIX_INDEX(a, j, LANCZOS2_POWER_COUNT)] *
          shiftedCoefficient[a];
    }
  }
  normalizationComplex =
      QuadraticComplex(coefficient, overlap, result->dimension);
  normalization = creal(normalizationComplex);
  if (!isfinite(normalization) || !isfinite(cimag(normalizationComplex)) ||
      !(normalization > 0.0)) {
    return LANCZOS2_SOLVE_NORMALIZATION_FAILURE;
  }
  normalization = sqrt(normalization);
  for (a = 0; a < result->dimension; a++) {
    coefficient[a] /= normalization;
    if (cabs(coefficient[a]) > maximum) {
      maximum = cabs(coefficient[a]);
      maximumIndex = a;
    }
  }
  phase = conj(coefficient[maximumIndex]) / maximum;
  for (a = 0; a < result->dimension; a++) coefficient[a] *= phase;
  for (a = 0; a < result->dimension; a++) {
    result->coefficient[a] = coefficient[a];
  }
  result->energy =
      creal(QuadraticComplex(coefficient, hamiltonian, result->dimension));
  h2 = creal(QuadraticComplex(coefficient, hamiltonianSquared,
                             result->dimension));
  if (!isfinite(result->energy) || !isfinite(h2) ||
      !isfinite(shiftedEigenvalue)) {
    return LANCZOS2_SOLVE_NONFINITE_MOMENT;
  }
  FinalizeVariance(result, h2);
  if (!isfinite(result->variance)) {
    return LANCZOS2_SOLVE_NONFINITE_MOMENT;
  }
  result->eigenvalueResidual =
      RelativeResidual(fabs((shiftedEigenvalue + result->shift) -
                            result->energy),
                       fmax(fabs(result->energy), 1.0));
  if (cabs(coefficient[0]) >
      64.0 * DBL_EPSILON * fmax(maximum, DBL_MIN)) {
    for (a = 1; a < result->dimension; a++) {
      result->alpha[a - 1] = coefficient[a] / coefficient[0];
    }
  } else {
    result->solveFlags |= LANCZOS2_FLAG_ALPHA_UNDEFINED;
  }
  return LANCZOS2_SOLVE_OK;
}

Lanczos2SolveStatus SolveLanczos2Real(
    const double moment[LANCZOS2_MOMENT_COUNT], int dimension,
    Lanczos2Result *result) {
  double projected[LANCZOS2_MOMENT_COUNT];
  double shifted[LANCZOS2_MOMENT_COUNT];
  double transform[LANCZOS2_MOMENT_COUNT];
  double overlap[9];
  double hamiltonian[9];
  double hamiltonianSquared[9];
  double shiftedOverlap[9];
  double shiftedHamiltonian[9];
  double shiftedHamiltonianSquared[9];
  double equilibratedOverlap[9];
  double equilibratedHamiltonian[9];
  double diagonalScale[LANCZOS2_MAX_DIMENSION];
  double eigenvector[LANCZOS2_MAX_DIMENSION];
  double shiftedEigenvalue = NAN;
  double maximumDiagonal = 0.0;
  Lanczos2SolveStatus status;
  int a;
  int b;
  int info;

  if (moment == NULL || result == NULL || !ValidDimension(dimension)) {
    return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  }
  ClearResult(result, dimension);
  status = ProjectMomentReal(moment, projected,
                             &result->antiHermitianResidual);
  if (status != LANCZOS2_SOLVE_OK) return status;
  result->hankelResidual = HankelResidualReal(projected);
  if (!(projected[0] > 0.0) || !isfinite(projected[0])) {
    return LANCZOS2_SOLVE_INVALID_M00;
  }
  result->shift = projected[1] / projected[0];
  if (!isfinite(result->shift)) return LANCZOS2_SOLVE_INVALID_M00;

  MatricesFromProjectedMomentReal(projected, dimension, overlap, hamiltonian,
                                  hamiltonianSquared);
  BuildShiftTransform(result->shift, transform);
  TransformMomentReal(transform, projected, shifted);
  MatricesFromProjectedMomentReal(shifted, dimension, shiftedOverlap,
                                  shiftedHamiltonian,
                                  shiftedHamiltonianSquared);

  for (a = 0; a < dimension; a++) {
    const double diagonal =
        shiftedOverlap[MATRIX_INDEX(a, a, dimension)];
    if (!isfinite(diagonal) || !(diagonal > 0.0)) {
      return LANCZOS2_SOLVE_NONPOSITIVE_DIAGONAL;
    }
    diagonalScale[a] = 1.0 / sqrt(diagonal);
  }
  for (a = 0; a < dimension; a++) {
    for (b = 0; b < dimension; b++) {
      equilibratedOverlap[MATRIX_INDEX(a, b, dimension)] =
          diagonalScale[a] *
          shiftedOverlap[MATRIX_INDEX(a, b, dimension)] *
          diagonalScale[b];
      equilibratedHamiltonian[MATRIX_INDEX(a, b, dimension)] =
          diagonalScale[a] *
          shiftedHamiltonian[MATRIX_INDEX(a, b, dimension)] *
          diagonalScale[b];
      if (!isfinite(equilibratedOverlap[MATRIX_INDEX(a, b, dimension)]) ||
          !isfinite(
              equilibratedHamiltonian[MATRIX_INDEX(a, b, dimension)])) {
        return LANCZOS2_SOLVE_NONFINITE_MOMENT;
      }
    }
    maximumDiagonal =
        fmax(maximumDiagonal,
             equilibratedOverlap[MATRIX_INDEX(a, a, dimension)]);
  }

  info = GeneralizedEigenReal(equilibratedHamiltonian,
                              equilibratedOverlap, dimension, 0.0,
                              &shiftedEigenvalue, eigenvector);
  if (info > dimension) {
    result->epsilon = fmax(maximumDiagonal, 1.0) * 1.0e-12;
    result->solveFlags |= LANCZOS2_FLAG_REGULARIZED;
    info = GeneralizedEigenReal(equilibratedHamiltonian,
                                equilibratedOverlap, dimension,
                                result->epsilon, &shiftedEigenvalue,
                                eigenvector);
  }
  result->lapackInfo = info;
  if (info != 0) return LANCZOS2_SOLVE_LAPACK_FAILURE;
  return FinalizeReal(overlap, hamiltonian, hamiltonianSquared, transform,
                      diagonalScale, eigenvector, shiftedEigenvalue, result);
}

Lanczos2SolveStatus SolveLanczos2Complex(
    const double complex moment[LANCZOS2_MOMENT_COUNT], int dimension,
    Lanczos2Result *result) {
  double complex projected[LANCZOS2_MOMENT_COUNT];
  double complex shifted[LANCZOS2_MOMENT_COUNT];
  double transform[LANCZOS2_MOMENT_COUNT];
  double complex overlap[9];
  double complex hamiltonian[9];
  double complex hamiltonianSquared[9];
  double complex shiftedOverlap[9];
  double complex shiftedHamiltonian[9];
  double complex shiftedHamiltonianSquared[9];
  double complex equilibratedOverlap[9];
  double complex equilibratedHamiltonian[9];
  double diagonalScale[LANCZOS2_MAX_DIMENSION];
  double complex eigenvector[LANCZOS2_MAX_DIMENSION];
  double shiftedEigenvalue = NAN;
  double maximumDiagonal = 0.0;
  Lanczos2SolveStatus status;
  int a;
  int b;
  int info;

  if (moment == NULL || result == NULL || !ValidDimension(dimension)) {
    return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  }
  ClearResult(result, dimension);
  status = ProjectMomentComplex(moment, projected,
                                &result->antiHermitianResidual);
  if (status != LANCZOS2_SOLVE_OK) return status;
  result->hankelResidual = HankelResidualComplex(projected);
  if (!(creal(projected[0]) > 0.0) || !isfinite(creal(projected[0]))) {
    return LANCZOS2_SOLVE_INVALID_M00;
  }
  result->shift = creal(projected[1]) / creal(projected[0]);
  if (!isfinite(result->shift)) return LANCZOS2_SOLVE_INVALID_M00;

  MatricesFromProjectedMomentComplex(projected, dimension, overlap,
                                     hamiltonian, hamiltonianSquared);
  BuildShiftTransform(result->shift, transform);
  TransformMomentComplex(transform, projected, shifted);
  MatricesFromProjectedMomentComplex(
      shifted, dimension, shiftedOverlap, shiftedHamiltonian,
      shiftedHamiltonianSquared);

  for (a = 0; a < dimension; a++) {
    const double diagonal =
        creal(shiftedOverlap[MATRIX_INDEX(a, a, dimension)]);
    if (!isfinite(diagonal) || !(diagonal > 0.0)) {
      return LANCZOS2_SOLVE_NONPOSITIVE_DIAGONAL;
    }
    diagonalScale[a] = 1.0 / sqrt(diagonal);
  }
  for (a = 0; a < dimension; a++) {
    for (b = 0; b < dimension; b++) {
      equilibratedOverlap[MATRIX_INDEX(a, b, dimension)] =
          diagonalScale[a] *
          shiftedOverlap[MATRIX_INDEX(a, b, dimension)] *
          diagonalScale[b];
      equilibratedHamiltonian[MATRIX_INDEX(a, b, dimension)] =
          diagonalScale[a] *
          shiftedHamiltonian[MATRIX_INDEX(a, b, dimension)] *
          diagonalScale[b];
      if (!isfinite(
              creal(equilibratedOverlap[MATRIX_INDEX(a, b, dimension)])) ||
          !isfinite(
              cimag(equilibratedOverlap[MATRIX_INDEX(a, b, dimension)])) ||
          !isfinite(creal(equilibratedHamiltonian[
              MATRIX_INDEX(a, b, dimension)])) ||
          !isfinite(cimag(equilibratedHamiltonian[
              MATRIX_INDEX(a, b, dimension)]))) {
        return LANCZOS2_SOLVE_NONFINITE_MOMENT;
      }
    }
    maximumDiagonal =
        fmax(maximumDiagonal,
             creal(equilibratedOverlap[MATRIX_INDEX(a, a, dimension)]));
  }

  info = GeneralizedEigenComplex(equilibratedHamiltonian,
                                 equilibratedOverlap, dimension, 0.0,
                                 &shiftedEigenvalue, eigenvector);
  if (info > dimension) {
    result->epsilon = fmax(maximumDiagonal, 1.0) * 1.0e-12;
    result->solveFlags |= LANCZOS2_FLAG_REGULARIZED;
    info = GeneralizedEigenComplex(equilibratedHamiltonian,
                                   equilibratedOverlap, dimension,
                                   result->epsilon, &shiftedEigenvalue,
                                   eigenvector);
  }
  result->lapackInfo = info;
  if (info != 0) return LANCZOS2_SOLVE_LAPACK_FAILURE;
  return FinalizeComplex(overlap, hamiltonian, hamiltonianSquared, transform,
                         diagonalScale, eigenvector, shiftedEigenvalue,
                         result);
}

static Lanczos2SolveStatus RestoreExternalShift(
    double basisShift, Lanczos2Result *result) {
  double transform[LANCZOS2_MOMENT_COUNT];
  double complex shiftedCoefficient[LANCZOS2_MAX_DIMENSION];
  double maximum = 0.0;
  double oldEnergyScale;
  double newEnergyScale;
  double complex phase;
  int maximumIndex = 0;
  int a;
  int j;
  if (result == NULL || !isfinite(basisShift)) {
    return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  }
  BuildShiftTransform(basisShift, transform);
  memcpy(shiftedCoefficient, result->coefficient,
         sizeof(shiftedCoefficient));
  memset(result->coefficient, 0, sizeof(result->coefficient));
  for (j = 0; j < result->dimension; j++) {
    for (a = j; a < result->dimension; a++) {
      result->coefficient[j] +=
          transform[MATRIX_INDEX(a, j, LANCZOS2_POWER_COUNT)] *
          shiftedCoefficient[a];
    }
  }
  for (a = 0; a < result->dimension; a++) {
    if (!isfinite(creal(result->coefficient[a])) ||
        !isfinite(cimag(result->coefficient[a]))) {
      return LANCZOS2_SOLVE_NONFINITE_MOMENT;
    }
    if (cabs(result->coefficient[a]) > maximum) {
      maximum = cabs(result->coefficient[a]);
      maximumIndex = a;
    }
  }
  if (!(maximum > 0.0) || !isfinite(maximum)) {
    return LANCZOS2_SOLVE_NORMALIZATION_FAILURE;
  }
  phase = conj(result->coefficient[maximumIndex]) / maximum;
  for (a = 0; a < result->dimension; a++) {
    result->coefficient[a] *= phase;
  }

  oldEnergyScale = fmax(fabs(result->energy), 1.0);
  result->energy += basisShift;
  result->shift += basisShift;
  if (!isfinite(result->energy) || !isfinite(result->shift)) {
    return LANCZOS2_SOLVE_NONFINITE_MOMENT;
  }
  newEnergyScale = fmax(fabs(result->energy), 1.0);
  result->eigenvalueResidual *= oldEnergyScale / newEnergyScale;
  if (result->energy * result->energy > DBL_MIN) {
    result->varianceOverEnergySquared =
        result->variance / (result->energy * result->energy);
  } else {
    result->varianceOverEnergySquared = NAN;
  }

  result->solveFlags &= ~LANCZOS2_FLAG_ALPHA_UNDEFINED;
  memset(result->alpha, 0, sizeof(result->alpha));
  if (cabs(result->coefficient[0]) >
      64.0 * DBL_EPSILON * fmax(maximum, DBL_MIN)) {
    for (a = 1; a < result->dimension; a++) {
      result->alpha[a - 1] =
          result->coefficient[a] / result->coefficient[0];
    }
  } else {
    result->solveFlags |= LANCZOS2_FLAG_ALPHA_UNDEFINED;
  }
  return LANCZOS2_SOLVE_OK;
}

Lanczos2SolveStatus SolveLanczos2ShiftedReal(
    const double moment[LANCZOS2_MOMENT_COUNT], int dimension,
    double basisShift, Lanczos2Result *result) {
  Lanczos2SolveStatus status;
  if (!isfinite(basisShift)) return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  status = SolveLanczos2Real(moment, dimension, result);
  if (status != LANCZOS2_SOLVE_OK) return status;
  return RestoreExternalShift(basisShift, result);
}

Lanczos2SolveStatus SolveLanczos2ShiftedComplex(
    const double complex moment[LANCZOS2_MOMENT_COUNT], int dimension,
    double basisShift, Lanczos2Result *result) {
  Lanczos2SolveStatus status;
  if (!isfinite(basisShift)) return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  status = SolveLanczos2Complex(moment, dimension, result);
  if (status != LANCZOS2_SOLVE_OK) return status;
  return RestoreExternalShift(basisShift, result);
}

const char *Lanczos2SolveError(Lanczos2SolveStatus status) {
  switch (status) {
    case LANCZOS2_SOLVE_OK:
      return "";
    case LANCZOS2_SOLVE_INVALID_ARGUMENT:
      return "invalid Lanczos2 solver argument";
    case LANCZOS2_SOLVE_NONFINITE_MOMENT:
      return "Lanczos2 moment or transformed matrix is not finite";
    case LANCZOS2_SOLVE_INVALID_M00:
      return "Lanczos2 M[0,0] or energy shift is invalid";
    case LANCZOS2_SOLVE_NONPOSITIVE_DIAGONAL:
      return "Lanczos2 shifted overlap has a non-positive diagonal";
    case LANCZOS2_SOLVE_LAPACK_FAILURE:
      return "Lanczos2 generalized eigensolver failed";
    case LANCZOS2_SOLVE_NORMALIZATION_FAILURE:
      return "Lanczos2 coefficient has non-positive overlap norm";
    case LANCZOS2_SOLVE_IO_FAILURE:
      return "Lanczos2 output write failed";
  }
  return "unknown Lanczos2 solver error";
}

Lanczos2SolveStatus AccumulateLanczos2MomentReal(
    double moment[LANCZOS2_MOMENT_COUNT],
    const double localPower[LANCZOS2_POWER_COUNT], double weight) {
  double updatedMoment[LANCZOS2_MOMENT_COUNT];
  int a;
  int b;
  if (moment == NULL || localPower == NULL || !isfinite(weight)) {
    return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  }
  for (a = 0; a < LANCZOS2_POWER_COUNT; a++) {
    if (!isfinite(localPower[a])) {
      return LANCZOS2_SOLVE_NONFINITE_MOMENT;
    }
  }
  for (a = 0; a < LANCZOS2_POWER_COUNT; a++) {
    for (b = 0; b < LANCZOS2_POWER_COUNT; b++) {
      const int index = MATRIX_INDEX(a, b, LANCZOS2_POWER_COUNT);
      const double term = weight * localPower[a] * localPower[b];
      updatedMoment[index] = moment[index] + term;
      if (!isfinite(moment[index]) || !isfinite(term) ||
          !isfinite(updatedMoment[index])) {
        return LANCZOS2_SOLVE_NONFINITE_MOMENT;
      }
    }
  }
  memcpy(moment, updatedMoment, sizeof(updatedMoment));
  return LANCZOS2_SOLVE_OK;
}

Lanczos2SolveStatus AccumulateLanczos2MomentComplex(
    double complex moment[LANCZOS2_MOMENT_COUNT],
    const double complex localPower[LANCZOS2_POWER_COUNT], double weight) {
  double complex updatedMoment[LANCZOS2_MOMENT_COUNT];
  int a;
  int b;
  if (moment == NULL || localPower == NULL || !isfinite(weight)) {
    return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  }
  for (a = 0; a < LANCZOS2_POWER_COUNT; a++) {
    if (!isfinite(creal(localPower[a])) ||
        !isfinite(cimag(localPower[a]))) {
      return LANCZOS2_SOLVE_NONFINITE_MOMENT;
    }
  }
  for (a = 0; a < LANCZOS2_POWER_COUNT; a++) {
    for (b = 0; b < LANCZOS2_POWER_COUNT; b++) {
      const int index = MATRIX_INDEX(a, b, LANCZOS2_POWER_COUNT);
      const double complex term =
          weight * conj(localPower[a]) * localPower[b];
      updatedMoment[index] = moment[index] + term;
      if (!isfinite(creal(moment[index])) ||
          !isfinite(cimag(moment[index])) || !isfinite(creal(term)) ||
          !isfinite(cimag(term)) || !isfinite(creal(updatedMoment[index])) ||
          !isfinite(cimag(updatedMoment[index]))) {
        return LANCZOS2_SOLVE_NONFINITE_MOMENT;
      }
    }
  }
  memcpy(moment, updatedMoment, sizeof(updatedMoment));
  return LANCZOS2_SOLVE_OK;
}

static int WriteLanczos2ResultRecord(FILE *output,
                                     const Lanczos2Result *result) {
  if (fprintf(
          output,
          "# mVMC ls2_out v1 E sigma2_over_E2 "
          "c0_re c0_im c1_re c1_im c2_re c2_im "
          "alpha1_re alpha1_im alpha2_re alpha2_im "
          "solve_flag epsilon shift antihermitian_residual "
          "hankel_residual\n") < 0) {
    return -1;
  }
  if (fprintf(
          output,
          "% .18e % .18e "
          "% .18e % .18e % .18e % .18e % .18e % .18e "
          "% .18e % .18e % .18e % .18e "
          "%d % .18e % .18e % .18e % .18e\n",
          result->energy, result->varianceOverEnergySquared,
          creal(result->coefficient[0]), cimag(result->coefficient[0]),
          creal(result->coefficient[1]), cimag(result->coefficient[1]),
          creal(result->coefficient[2]), cimag(result->coefficient[2]),
          creal(result->alpha[0]), cimag(result->alpha[0]),
          creal(result->alpha[1]), cimag(result->alpha[1]),
          result->solveFlags, result->epsilon, result->shift,
          result->antiHermitianResidual, result->hankelResidual) < 0) {
    return -1;
  }
  return ferror(output) ? -1 : 0;
}

static int WriteLanczos2MomentHeader(FILE *momentOutput,
                                     double basisShift) {
  int a;
  int b;
  if (fprintf(momentOutput,
              "# mVMC ls2_moment v2 basis_shift=% .18e",
              basisShift) < 0) {
    return -1;
  }
  for (a = 0; a < LANCZOS2_POWER_COUNT; a++) {
    for (b = 0; b < LANCZOS2_POWER_COUNT; b++) {
      if (fprintf(momentOutput, " M%d%d_re M%d%d_im", a, b, a, b) < 0) {
        return -1;
      }
    }
  }
  return fprintf(momentOutput, "\n") < 0 ? -1 : 0;
}

Lanczos2SolveStatus WriteLanczos2OutputReal(
    FILE *output, FILE *momentOutput,
    const double moment[LANCZOS2_MOMENT_COUNT], double basisShift,
    Lanczos2Result *result) {
  Lanczos2SolveStatus status;
  int i;
  if (output == NULL || momentOutput == NULL || moment == NULL ||
      result == NULL) {
    return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  }
  status = SolveLanczos2ShiftedReal(moment, 3, basisShift, result);
  if (status != LANCZOS2_SOLVE_OK) return status;
  if (WriteLanczos2ResultRecord(output, result) != 0 ||
      WriteLanczos2MomentHeader(momentOutput, basisShift) != 0) {
    return LANCZOS2_SOLVE_IO_FAILURE;
  }
  for (i = 0; i < LANCZOS2_MOMENT_COUNT; i++) {
    if (fprintf(momentOutput, "%s% .18e % .18e",
                i == 0 ? "" : " ", moment[i], 0.0) < 0) {
      return LANCZOS2_SOLVE_IO_FAILURE;
    }
  }
  if (fprintf(momentOutput, "\n") < 0 || ferror(momentOutput)) {
    return LANCZOS2_SOLVE_IO_FAILURE;
  }
  return LANCZOS2_SOLVE_OK;
}

Lanczos2SolveStatus WriteLanczos2OutputComplex(
    FILE *output, FILE *momentOutput,
    const double complex moment[LANCZOS2_MOMENT_COUNT],
    double basisShift, Lanczos2Result *result) {
  Lanczos2SolveStatus status;
  int i;
  if (output == NULL || momentOutput == NULL || moment == NULL ||
      result == NULL) {
    return LANCZOS2_SOLVE_INVALID_ARGUMENT;
  }
  status = SolveLanczos2ShiftedComplex(moment, 3, basisShift, result);
  if (status != LANCZOS2_SOLVE_OK) return status;
  if (WriteLanczos2ResultRecord(output, result) != 0 ||
      WriteLanczos2MomentHeader(momentOutput, basisShift) != 0) {
    return LANCZOS2_SOLVE_IO_FAILURE;
  }
  for (i = 0; i < LANCZOS2_MOMENT_COUNT; i++) {
    if (fprintf(momentOutput, "%s% .18e % .18e",
                i == 0 ? "" : " ", creal(moment[i]),
                cimag(moment[i])) < 0) {
      return LANCZOS2_SOLVE_IO_FAILURE;
    }
  }
  if (fprintf(momentOutput, "\n") < 0 || ferror(momentOutput)) {
    return LANCZOS2_SOLVE_IO_FAILURE;
  }
  return LANCZOS2_SOLVE_OK;
}
