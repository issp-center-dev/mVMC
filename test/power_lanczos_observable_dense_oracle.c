#include "power_lanczos_observable_dense_oracle.h"

#include <math.h>
#include <string.h>

/*
 * Independent dense oracle: it constructs every single creation/annihilation
 * matrix on the full Fock space, multiplies those dense matrices, and only
 * then restricts vectors to the requested sector.  It deliberately imports
 * no production observable action/evaluator header or implementation.
 */

static int Popcount(uint64_t value) {
  int count = 0;
  while (value != 0) {
    value &= value - UINT64_C(1);
    ++count;
  }
  return count;
}

static int FiniteComplex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

int oracle_build_basis(const OracleLayout *layout, OracleBasis *basis) {
  uint64_t word;
  uint64_t dimension;
  if (layout == NULL || basis == NULL || layout->nsite <= 0 ||
      2 * layout->nsite > ORACLE_MAX_ORBITALS ||
      layout->up_electron_count < 0 ||
      layout->up_electron_count > layout->nsite ||
      layout->down_electron_count < 0 ||
      layout->down_electron_count > layout->nsite ||
      (layout->pure_spin != 0 && layout->pure_spin != 1) ||
      (layout->pure_spin &&
       layout->up_electron_count + layout->down_electron_count !=
           layout->nsite)) {
    return 0;
  }
  memset(basis, 0, sizeof(*basis));
  dimension = UINT64_C(1) << (2 * layout->nsite);
  basis->full_dimension = (size_t)dimension;
  for (word = 0; word < dimension; ++word) {
    int up = 0;
    int down = 0;
    int valid = 1;
    int site;
    for (site = 0; site < layout->nsite; ++site) {
      const int site_up = (int)((word >> site) & UINT64_C(1));
      const int site_down =
          (int)((word >> (site + layout->nsite)) & UINT64_C(1));
      up += site_up;
      down += site_down;
      if (layout->pure_spin && site_up + site_down != 1) valid = 0;
    }
    if (valid && up == layout->up_electron_count &&
        down == layout->down_electron_count) {
      basis->sector_words[basis->sector_dimension++] = word;
    }
  }
  return basis->sector_dimension > 0;
}

static int BuildSingleOperator(int orbital_count, size_t orbital,
                               int creation, double complex *matrix) {
  const size_t dimension = (size_t)1 << orbital_count;
  size_t source;
  if (orbital >= (size_t)orbital_count) return 0;
  memset(matrix, 0,
         ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM * sizeof(matrix[0]));
  for (source = 0; source < dimension; ++source) {
    const uint64_t mask = UINT64_C(1) << orbital;
    const int occupied = (((uint64_t)source & mask) != 0);
    uint64_t target;
    int sign;
    if ((creation && occupied) || (!creation && !occupied)) continue;
    target = creation ? (uint64_t)source | mask
                      : (uint64_t)source & ~mask;
    sign = (Popcount((uint64_t)source & (mask - UINT64_C(1))) & 1)
               ? -1
               : 1;
    matrix[(size_t)target * dimension + source] = (double)sign;
  }
  return 1;
}

static void MatrixMultiply(const double complex *left,
                           const double complex *right,
                           size_t dimension, double complex *result) {
  size_t row;
  memset(result, 0,
         ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM * sizeof(result[0]));
  for (row = 0; row < dimension; ++row) {
    size_t middle;
    for (middle = 0; middle < dimension; ++middle) {
      const double complex factor = left[row * dimension + middle];
      size_t column;
      if (creal(factor) == 0.0 && cimag(factor) == 0.0) continue;
      for (column = 0; column < dimension; ++column) {
        result[row * dimension + column] +=
            factor * right[middle * dimension + column];
      }
    }
  }
}

int oracle_build_monomial_matrix(
    int orbital_count, const OracleMonomial *monomial,
    double complex matrix[ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM]) {
  double complex current[ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM];
  double complex single[ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM];
  double complex product[ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM];
  const size_t dimension = (size_t)1 << orbital_count;
  size_t index;
  if (monomial == NULL || matrix == NULL || orbital_count <= 0 ||
      orbital_count > ORACLE_MAX_ORBITALS ||
      (monomial->operator_count != 2 && monomial->operator_count != 4)) {
    return 0;
  }
  memset(current, 0, sizeof(current));
  for (index = 0; index < dimension; ++index) {
    current[index * dimension + index] = 1.0;
  }
  for (index = 0; index < monomial->operator_count; ++index) {
    if (!BuildSingleOperator(orbital_count, monomial->orbital[index],
                             (index % 2) == 0, single)) {
      return 0;
    }
    MatrixMultiply(current, single, dimension, product);
    memcpy(current, product, sizeof(current));
  }
  memcpy(matrix, current, sizeof(current));
  return 1;
}

int oracle_build_hamiltonian_matrix(
    int orbital_count, const OracleTerm *terms, size_t term_count,
    double complex matrix[ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM]) {
  double complex monomial[ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM];
  const size_t dimension = (size_t)1 << orbital_count;
  size_t term;
  size_t entry;
  if (matrix == NULL || orbital_count <= 0 ||
      orbital_count > ORACLE_MAX_ORBITALS ||
      (term_count != 0 && terms == NULL) || term_count > ORACLE_MAX_TERMS) {
    return 0;
  }
  memset(matrix, 0,
         ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM * sizeof(matrix[0]));
  for (term = 0; term < term_count; ++term) {
    if (!FiniteComplex(terms[term].coefficient) ||
        !oracle_build_monomial_matrix(orbital_count,
                                      &terms[term].monomial, monomial)) {
      return 0;
    }
    for (entry = 0; entry < dimension * dimension; ++entry) {
      matrix[entry] += terms[term].coefficient * monomial[entry];
    }
  }
  return 1;
}

void oracle_matrix_vector(
    const double complex matrix[ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM],
    size_t dimension, const double complex *vector,
    double complex *result) {
  size_t row;
  for (row = 0; row < dimension; ++row) {
    size_t column;
    double complex value = 0.0;
    for (column = 0; column < dimension; ++column) {
      value += matrix[row * dimension + column] * vector[column];
    }
    result[row] = value;
  }
}

int oracle_observable_matrix(
    int orbital_count, const OracleMonomial *monomial,
    const double complex *basis0, const double complex *basis1,
    double complex result[4]) {
  double complex matrix[ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM];
  const double complex *basis[2] = {basis0, basis1};
  const size_t dimension = (size_t)1 << orbital_count;
  size_t row_basis;
  if (basis0 == NULL || basis1 == NULL || result == NULL ||
      !oracle_build_monomial_matrix(orbital_count, monomial, matrix)) {
    return 0;
  }
  for (row_basis = 0; row_basis < 2; ++row_basis) {
    size_t column_basis;
    for (column_basis = 0; column_basis < 2; ++column_basis) {
      size_t target;
      double complex value = 0.0;
      for (target = 0; target < dimension; ++target) {
        size_t source;
        for (source = 0; source < dimension; ++source) {
          value += conj(basis[row_basis][target]) *
                   matrix[target * dimension + source] *
                   basis[column_basis][source];
        }
      }
      result[2 * row_basis + column_basis] = value;
    }
  }
  return 1;
}

int oracle_observable_final(
    int orbital_count, const OracleMonomial *monomial,
    const double complex *basis0, const double complex *basis1,
    const double complex alpha[2], double complex *result) {
  double complex final[ORACLE_MAX_FULL_DIM];
  double complex matrix[4];
  const size_t dimension = (size_t)1 << orbital_count;
  size_t index;
  if (alpha == NULL || result == NULL) return 0;
  for (index = 0; index < dimension; ++index) {
    final[index] = alpha[0] * basis0[index] + alpha[1] * basis1[index];
  }
  if (!oracle_observable_matrix(orbital_count, monomial, final, final,
                                matrix)) {
    return 0;
  }
  *result = matrix[0];
  return 1;
}
