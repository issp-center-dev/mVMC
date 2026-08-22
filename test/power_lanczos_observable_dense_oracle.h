#ifndef TEST_POWER_LANCZOS_OBSERVABLE_DENSE_ORACLE_H
#define TEST_POWER_LANCZOS_OBSERVABLE_DENSE_ORACLE_H

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define ORACLE_MAX_ORBITALS 6
#define ORACLE_MAX_FULL_DIM (1u << ORACLE_MAX_ORBITALS)
#define ORACLE_MAX_TERMS 32

typedef struct {
  int nsite;
  int up_electron_count;
  int down_electron_count;
  int pure_spin;
} OracleLayout;

typedef struct {
  size_t operator_count;
  size_t orbital[4];
} OracleMonomial;

typedef struct {
  double complex coefficient;
  OracleMonomial monomial;
} OracleTerm;

typedef struct {
  size_t full_dimension;
  size_t sector_dimension;
  uint64_t sector_words[ORACLE_MAX_FULL_DIM];
} OracleBasis;

int oracle_build_basis(const OracleLayout *layout, OracleBasis *basis);

int oracle_build_monomial_matrix(
    int orbital_count, const OracleMonomial *monomial,
    double complex matrix[ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM]);

int oracle_build_hamiltonian_matrix(
    int orbital_count, const OracleTerm *terms, size_t term_count,
    double complex matrix[ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM]);

void oracle_matrix_vector(
    const double complex matrix[ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM],
    size_t dimension, const double complex *vector,
    double complex *result);

int oracle_observable_matrix(
    int orbital_count, const OracleMonomial *monomial,
    const double complex *basis0, const double complex *basis1,
    double complex result[4]);

int oracle_observable_final(
    int orbital_count, const OracleMonomial *monomial,
    const double complex *basis0, const double complex *basis1,
    const double complex alpha[2], double complex *result);

#endif
