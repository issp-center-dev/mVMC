#ifndef _PHYSCAL_LANCZOS2_H
#define _PHYSCAL_LANCZOS2_H

#include <complex.h>
#include <stddef.h>
#include <stdio.h>

enum {
  LANCZOS2_MAX_DIMENSION = 3,
  LANCZOS2_POWER_COUNT = 4,
  LANCZOS2_MOMENT_COUNT = 16
};

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
