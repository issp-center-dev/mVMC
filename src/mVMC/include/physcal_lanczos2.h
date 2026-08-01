#ifndef _PHYSCAL_LANCZOS2_H
#define _PHYSCAL_LANCZOS2_H

#include <complex.h>
#include <stddef.h>
#include <stdio.h>

enum {
  LANCZOS2_MAX_DIMENSION = 3,
  LANCZOS2_POWER_COUNT = 4,
  LANCZOS2_MOMENT_COUNT = 16,
  POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH = 10,
  POWER_LANCZOS_SUPPORT_MIN_SAMPLES = 32,
  POWER_LANCZOS_SUPPORT_MAX_BLOCKS = 16
};

typedef enum {
  POWER_LANCZOS_SUPPORT_PASS = 0,
  POWER_LANCZOS_SUPPORT_MISMATCH,
  POWER_LANCZOS_SUPPORT_INCONCLUSIVE,
  POWER_LANCZOS_SUPPORT_INVALID
} PowerLanczosSupportStatus;

typedef struct {
  PowerLanczosSupportStatus status;
  size_t sampleCount;
  int blockCount;
  double totalWeight;
  double effectiveSampleCount;
  int hasThirdAnchor;
  double complex m02;
  double complex m11;
  double complex m03;
  double complex m12;
  double complex delta2;
  double complex delta3;
  double relativeDifference2;
  double relativeDifference3;
  double standardError2;
  double standardError3;
  double score2;
  double score3;
} PowerLanczosSupportDiagnostic;

typedef enum {
  LANCZOS2_SOLVE_OK = 0,
  LANCZOS2_SOLVE_INVALID_ARGUMENT,
  LANCZOS2_SOLVE_NONFINITE_MOMENT,
  LANCZOS2_SOLVE_INVALID_M00,
  LANCZOS2_SOLVE_NONPOSITIVE_DIAGONAL,
  LANCZOS2_SOLVE_LAPACK_FAILURE,
  LANCZOS2_SOLVE_NORMALIZATION_FAILURE,
  LANCZOS2_SOLVE_IO_FAILURE
} Lanczos2SolveStatus;

enum {
  LANCZOS2_FLAG_REGULARIZED = 1 << 0,
  LANCZOS2_FLAG_ALPHA_UNDEFINED = 1 << 1
};

typedef struct {
  int dimension;
  int solveFlags;
  int lapackInfo;
  double energy;
  double variance;
  double varianceOverEnergySquared;
  double complex coefficient[LANCZOS2_MAX_DIMENSION];
  double complex alpha[LANCZOS2_MAX_DIMENSION - 1];
  double epsilon;
  double shift;
  double antiHermitianResidual;
  double hankelResidual;
  double eigenvalueResidual;
} Lanczos2Result;

Lanczos2SolveStatus BuildLanczos2MatricesReal(
    const double moment[LANCZOS2_MOMENT_COUNT], int dimension,
    double *overlap, double *hamiltonian, double *hamiltonianSquared,
    double *antiHermitianResidual, double *hankelResidual);

Lanczos2SolveStatus BuildLanczos2MatricesComplex(
    const double complex moment[LANCZOS2_MOMENT_COUNT], int dimension,
    double complex *overlap, double complex *hamiltonian,
    double complex *hamiltonianSquared, double *antiHermitianResidual,
    double *hankelResidual);

Lanczos2SolveStatus SolveLanczos2Real(
    const double moment[LANCZOS2_MOMENT_COUNT], int dimension,
    Lanczos2Result *result);

Lanczos2SolveStatus SolveLanczos2Complex(
    const double complex moment[LANCZOS2_MOMENT_COUNT], int dimension,
    Lanczos2Result *result);

Lanczos2SolveStatus SolveLanczos2ShiftedReal(
    const double moment[LANCZOS2_MOMENT_COUNT], int dimension,
    double basisShift, Lanczos2Result *result);

Lanczos2SolveStatus SolveLanczos2ShiftedComplex(
    const double complex moment[LANCZOS2_MOMENT_COUNT], int dimension,
    double basisShift, Lanczos2Result *result);

Lanczos2SolveStatus BuildLanczos2ShiftedMomentReal(
    double moment[LANCZOS2_MOMENT_COUNT], const double *samplePower,
    const double *sampleWeight, const unsigned char *sampleValid,
    size_t sampleCount, double shift);

Lanczos2SolveStatus BuildLanczos2ShiftedMomentComplex(
    double complex moment[LANCZOS2_MOMENT_COUNT],
    const double complex *samplePower, const double *sampleWeight,
    const unsigned char *sampleValid, size_t sampleCount, double shift);

Lanczos2SolveStatus AccumulateLanczos2MomentReal(
    double moment[LANCZOS2_MOMENT_COUNT],
    const double localPower[LANCZOS2_POWER_COUNT], double weight);

Lanczos2SolveStatus AccumulateLanczos2MomentComplex(
    double complex moment[LANCZOS2_MOMENT_COUNT],
    const double complex localPower[LANCZOS2_POWER_COUNT], double weight);

Lanczos2SolveStatus RecordPowerLanczosSupportSampleReal(
    double sample[POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH],
    const double localPower[LANCZOS2_POWER_COUNT], int powerCount,
    double weight);

Lanczos2SolveStatus RecordPowerLanczosSupportSampleComplex(
    double sample[POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH],
    const double complex localPower[LANCZOS2_POWER_COUNT], int powerCount,
    double weight);

PowerLanczosSupportStatus AnalyzePowerLanczosSupport(
    const double *sampleData, size_t sampleCapacity, int hasThirdAnchor,
    PowerLanczosSupportDiagnostic *diagnostic);

Lanczos2SolveStatus WritePowerLanczosSupportDiagnostic(
    FILE *output, int lanczosStep, int experimentalMode,
    const PowerLanczosSupportDiagnostic *diagnostic);

const char *PowerLanczosSupportStatusName(
    PowerLanczosSupportStatus status);

Lanczos2SolveStatus WriteLanczos2OutputReal(
    FILE *output, FILE *momentOutput,
    const double moment[LANCZOS2_MOMENT_COUNT], double basisShift,
    Lanczos2Result *result);

Lanczos2SolveStatus WriteLanczos2OutputComplex(
    FILE *output, FILE *momentOutput,
    const double complex moment[LANCZOS2_MOMENT_COUNT],
    double basisShift, Lanczos2Result *result);

const char *Lanczos2SolveError(Lanczos2SolveStatus status);

#endif
