#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#include "blas_externs.h"
#include "krylov_gevp_solver.h"

#define MATRIX_INDEX(row, column, dimension) \
  ((size_t)(row) * (size_t)(dimension) + (size_t)(column))

enum { MVMC_KRYLOV_GEVP_WORK_SIZE = 256 };

static int finite_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static int valid_dimension(int dimension) {
  return dimension == 2 || dimension == 3;
}

static int valid_eigensolver_dimension(int dimension) {
  return dimension >= 1 && dimension <= MVMC_KRYLOV_GEVP_MAX_DIMENSION;
}

static size_t required_upper_count(int dimension) {
  return (size_t)dimension * ((size_t)dimension + 1U) / 2U;
}

static int valid_policy(const MVMCKrylovGEVPPolicy *policy) {
  return policy != NULL && policy->valid &&
         policy->policy_version == MVMC_KRYLOV_GEVP_POLICY_VERSION &&
         isfinite(policy->rank_relative_cutoff) &&
         policy->rank_relative_cutoff > 0.0 &&
         policy->rank_relative_cutoff < 1.0 &&
         isfinite(policy->overlap_negative_relative_tolerance) &&
         policy->overlap_negative_relative_tolerance >= 0.0 &&
         policy->overlap_negative_relative_tolerance < 1.0 &&
         isfinite(
             policy->hamiltonian_squared_negative_relative_tolerance) &&
         policy->hamiltonian_squared_negative_relative_tolerance >= 0.0 &&
         policy->hamiltonian_squared_negative_relative_tolerance < 1.0 &&
         isfinite(policy->maximum_input_antihermitian_effect) &&
         policy->maximum_input_antihermitian_effect >= 0.0 &&
         policy->maximum_input_antihermitian_effect < 1.0 &&
         isfinite(policy->maximum_gevp_relative_residual) &&
         policy->maximum_gevp_relative_residual > 0.0 &&
         policy->maximum_gevp_relative_residual < 1.0 &&
         isfinite(policy->negative_variance_relative_tolerance) &&
         policy->negative_variance_relative_tolerance >= 0.0 &&
         policy->negative_variance_relative_tolerance < 1.0 &&
         isfinite(policy->degenerate_root_gap_relative_threshold) &&
         policy->degenerate_root_gap_relative_threshold >= 0.0 &&
         policy->degenerate_root_gap_relative_threshold < 1.0;
}

static void invalidate_result(MVMCKrylovGEVPStatus status, int dimension,
                              MVMCKrylovGEVPResult *result) {
  int index;
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
  result->policy_version = MVMC_KRYLOV_GEVP_POLICY_VERSION;
  result->dimension = dimension;
  result->phase_pivot = -1;
  result->lapack_info_overlap_query = -1;
  result->lapack_info_overlap = -1;
  result->lapack_info_hamiltonian_squared_query = -1;
  result->lapack_info_hamiltonian_squared = -1;
  result->lapack_info_projected_query = -1;
  result->lapack_info_projected = -1;
  result->rank_relative_cutoff = NAN;
  result->condition_estimate = NAN;
  result->overlap_diagonal_imaginary_residual = NAN;
  result->hamiltonian_antihermitian_residual = NAN;
  result->hamiltonian_squared_diagonal_imaginary_residual = NAN;
  result->overlap_norm = NAN;
  result->hamiltonian_norm = NAN;
  result->hamiltonian_squared_norm = NAN;
  result->projected_energy = NAN;
  result->root_gap = NAN;
  result->normalization = NAN;
  result->energy = NAN;
  result->energy_squared = NAN;
  result->variance = NAN;
  result->gevp_relative_residual = NAN;
  for (index = 0; index < MVMC_KRYLOV_GEVP_MAX_DIMENSION; ++index) {
    result->overlap_eigenvalue[index] = NAN;
    result->hamiltonian_squared_eigenvalue[index] = NAN;
    result->coefficient[index] = NAN + I * NAN;
  }
  for (index = 0;
       index < MVMC_KRYLOV_GEVP_MAX_DIMENSION *
                   MVMC_KRYLOV_GEVP_MAX_DIMENSION;
       ++index) {
    result->root_subspace_projector[index] = NAN + I * NAN;
  }
}

static int byte_ranges_overlap(const void *left, size_t left_bytes,
                               const void *right, size_t right_bytes) {
  uintptr_t left_begin;
  uintptr_t right_begin;
  uintptr_t left_end;
  uintptr_t right_end;
  if (left == NULL || right == NULL || left_bytes == 0 || right_bytes == 0) {
    return 0;
  }
  left_begin = (uintptr_t)left;
  right_begin = (uintptr_t)right;
  if (left_begin > UINTPTR_MAX - left_bytes ||
      right_begin > UINTPTR_MAX - right_bytes) {
    return 1;
  }
  left_end = left_begin + left_bytes;
  right_end = right_begin + right_bytes;
  return left_begin < right_end && right_begin < left_end;
}

static int storage_is_disjoint(const void *overlap, const void *forward,
                               const void *reverse, const void *squared,
                               size_t input_bytes,
                               const MVMCKrylovGEVPPolicy *policy,
                               const MVMCKrylovGEVPResult *result) {
  const void *input[4] = {overlap, forward, reverse, squared};
  int left;
  int right;
  for (left = 0; left < 4; ++left) {
    if (byte_ranges_overlap(input[left], input_bytes, result,
                            sizeof(*result))) {
      return 0;
    }
    for (right = left + 1; right < 4; ++right) {
      if (byte_ranges_overlap(input[left], input_bytes, input[right],
                              input_bytes)) {
        return 0;
      }
    }
    if (byte_ranges_overlap(input[left], input_bytes, policy,
                            sizeof(*policy))) {
      return 0;
    }
  }
  return !byte_ranges_overlap(policy, sizeof(*policy), result,
                              sizeof(*result));
}

static void norm_add(double value, double *norm) {
  *norm = hypot(*norm, value);
}

static void complex_norm_add(double complex value, double *norm) {
  norm_add(creal(value), norm);
  norm_add(cimag(value), norm);
}

static double dense_frobenius_norm(const double complex *matrix,
                                   int dimension) {
  double norm = 0.0;
  int row;
  int column;
  for (row = 0; row < dimension; ++row) {
    for (column = 0; column < dimension; ++column) {
      complex_norm_add(matrix[MATRIX_INDEX(row, column, dimension)], &norm);
    }
  }
  return norm;
}

static void row_major_to_column_major(const double complex *source,
                                      double complex *destination,
                                      int dimension) {
  int row;
  int column;
  for (row = 0; row < dimension; ++row) {
    for (column = 0; column < dimension; ++column) {
      destination[row + column * dimension] =
          source[MATRIX_INDEX(row, column, dimension)];
    }
  }
}

static int hermitian_eigensolve(const double complex *matrix, int dimension,
                                double *eigenvalue,
                                double complex *eigenvector_column_major,
                                int *query_info, int *solve_info) {
  double complex query_matrix[MVMC_KRYLOV_GEVP_MAX_DIMENSION *
                              MVMC_KRYLOV_GEVP_MAX_DIMENSION];
  double complex work_query = 0.0;
  double complex work[MVMC_KRYLOV_GEVP_WORK_SIZE];
  double rwork[3 * MVMC_KRYLOV_GEVP_MAX_DIMENSION - 2];
  const char jobz = 'V';
  const char uplo = 'U';
  const int leading_dimension = dimension;
  int lwork = -1;
  int info = 0;
  int requested_work;
  if (matrix == NULL || eigenvalue == NULL ||
      eigenvector_column_major == NULL || query_info == NULL ||
      solve_info == NULL || !valid_eigensolver_dimension(dimension)) {
    return 0;
  }
  row_major_to_column_major(matrix, query_matrix, dimension);
  M_ZHEEV(&jobz, &uplo, &dimension, query_matrix, &leading_dimension,
          eigenvalue, &work_query, &lwork, rwork, &info);
  *query_info = info;
  if (info != 0 || !isfinite(creal(work_query)) ||
      creal(work_query) < 1.0 || creal(work_query) > (double)INT32_MAX) {
    *solve_info = -1;
    return 0;
  }
  requested_work = (int)ceil(creal(work_query));
  if (requested_work > MVMC_KRYLOV_GEVP_WORK_SIZE) {
    *solve_info = -1;
    return 0;
  }
  row_major_to_column_major(matrix, eigenvector_column_major, dimension);
  lwork = MVMC_KRYLOV_GEVP_WORK_SIZE;
  info = 0;
  M_ZHEEV(&jobz, &uplo, &dimension, eigenvector_column_major,
          &leading_dimension, eigenvalue, work, &lwork, rwork, &info);
  *solve_info = info;
  return info == 0;
}

static MVMCKrylovGEVPStatus build_hermitian_matrices(
    int dimension, const double complex *overlap_upper,
    const double complex *forward_upper,
    const double complex *reverse_upper,
    const double complex *squared_upper, double complex *overlap,
    double complex *hamiltonian, double complex *squared,
    MVMCKrylovGEVPResult *result) {
  double overlap_imaginary_norm = 0.0;
  double overlap_raw_norm = 0.0;
  double squared_imaginary_norm = 0.0;
  double squared_raw_norm = 0.0;
  double hamiltonian_residual_norm = 0.0;
  double hamiltonian_raw_norm = 0.0;
  size_t entry = 0;
  int row;
  int column;
  for (row = 0; row < dimension; ++row) {
    for (column = row; column < dimension; ++column, ++entry) {
      const double complex overlap_value = overlap_upper[entry];
      const double complex forward_value = forward_upper[entry];
      const double complex reverse_value = reverse_upper[entry];
      const double complex squared_value = squared_upper[entry];
      const double complex delta = forward_value - conj(reverse_value);
      double complex hermitian_value;
      if (!finite_complex(overlap_value) || !finite_complex(forward_value) ||
          !finite_complex(reverse_value) || !finite_complex(squared_value)) {
        return MVMC_KRYLOV_GEVP_NONFINITE_INPUT;
      }
      if (row == column) {
        norm_add(2.0 * cimag(overlap_value), &overlap_imaginary_norm);
        norm_add(2.0 * cimag(squared_value), &squared_imaginary_norm);
      }
      complex_norm_add(overlap_value, &overlap_raw_norm);
      complex_norm_add(squared_value, &squared_raw_norm);
      complex_norm_add(delta, &hamiltonian_residual_norm);
      complex_norm_add(forward_value, &hamiltonian_raw_norm);
      if (row != column) {
        complex_norm_add(overlap_value, &overlap_raw_norm);
        complex_norm_add(squared_value, &squared_raw_norm);
        complex_norm_add(delta, &hamiltonian_residual_norm);
        complex_norm_add(reverse_value, &hamiltonian_raw_norm);
      }
      overlap[MATRIX_INDEX(row, column, dimension)] =
          row == column ? creal(overlap_value) : overlap_value;
      squared[MATRIX_INDEX(row, column, dimension)] =
          row == column ? creal(squared_value) : squared_value;
      hermitian_value =
          0.5 * forward_value + 0.5 * conj(reverse_value);
      hamiltonian[MATRIX_INDEX(row, column, dimension)] =
          row == column ? creal(hermitian_value) : hermitian_value;
      if (row != column) {
        overlap[MATRIX_INDEX(column, row, dimension)] = conj(overlap_value);
        squared[MATRIX_INDEX(column, row, dimension)] = conj(squared_value);
        hamiltonian[MATRIX_INDEX(column, row, dimension)] =
            conj(hermitian_value);
      }
    }
  }
  result->overlap_norm = dense_frobenius_norm(overlap, dimension);
  result->hamiltonian_norm = dense_frobenius_norm(hamiltonian, dimension);
  result->hamiltonian_squared_norm =
      dense_frobenius_norm(squared, dimension);
  result->overlap_diagonal_imaginary_residual =
      overlap_raw_norm > 0.0
          ? overlap_imaginary_norm / overlap_raw_norm
          : overlap_imaginary_norm;
  result->hamiltonian_antihermitian_residual =
      hamiltonian_raw_norm > 0.0
          ? hamiltonian_residual_norm / hamiltonian_raw_norm
          : hamiltonian_residual_norm;
  result->hamiltonian_squared_diagonal_imaginary_residual =
      squared_raw_norm > 0.0
          ? squared_imaginary_norm / squared_raw_norm
          : squared_imaginary_norm;
  if (!isfinite(result->overlap_norm) ||
      !isfinite(result->hamiltonian_norm) ||
      !isfinite(result->hamiltonian_squared_norm) ||
      !isfinite(result->overlap_diagonal_imaginary_residual) ||
      !isfinite(result->hamiltonian_antihermitian_residual) ||
      !isfinite(result->hamiltonian_squared_diagonal_imaginary_residual)) {
    return MVMC_KRYLOV_GEVP_NONFINITE_INPUT;
  }
  return MVMC_KRYLOV_GEVP_OK;
}

static double complex quadratic_form(const double complex *coefficient,
                                     const double complex *matrix,
                                     int dimension) {
  double complex value = 0.0;
  int row;
  int column;
  for (row = 0; row < dimension; ++row) {
    for (column = 0; column < dimension; ++column) {
      value += conj(coefficient[row]) *
               matrix[MATRIX_INDEX(row, column, dimension)] *
               coefficient[column];
    }
  }
  return value;
}

static double vector_norm(const double complex *value, int dimension) {
  double norm = 0.0;
  int index;
  for (index = 0; index < dimension; ++index) {
    complex_norm_add(value[index], &norm);
  }
  return norm;
}

static int build_root_subspace_projector(
    int dimension, int retained_rank, int root_multiplicity,
    const double *inverse_diagonal, const double complex *whitening,
    const double complex *projected_eigenvectors,
    double complex *projector) {
  double complex basis[MVMC_KRYLOV_GEVP_MAX_DIMENSION *
                       MVMC_KRYLOV_GEVP_MAX_DIMENSION] = {0.0};
  int root;
  int row;
  int column;
  if (!valid_dimension(dimension) || retained_rank <= 0 ||
      retained_rank > dimension || root_multiplicity <= 0 ||
      root_multiplicity > retained_rank || inverse_diagonal == NULL ||
      whitening == NULL || projected_eigenvectors == NULL ||
      projector == NULL) {
    return 0;
  }
  for (root = 0; root < root_multiplicity; ++root) {
    double complex vector[MVMC_KRYLOV_GEVP_MAX_DIMENSION] = {0.0};
    double norm;
    int pass;
    int previous;
    for (row = 0; row < dimension; ++row) {
      int retained;
      for (retained = 0; retained < retained_rank; ++retained) {
        vector[row] +=
            whitening[MATRIX_INDEX(row, retained, retained_rank)] *
            projected_eigenvectors[retained + root * retained_rank];
      }
      vector[row] *= inverse_diagonal[row];
    }
    for (pass = 0; pass < 2; ++pass) {
      for (previous = 0; previous < root; ++previous) {
        double complex overlap = 0.0;
        for (row = 0; row < dimension; ++row) {
          overlap +=
              conj(basis[MATRIX_INDEX(row, previous, root_multiplicity)]) *
              vector[row];
        }
        for (row = 0; row < dimension; ++row) {
          vector[row] -=
              overlap *
              basis[MATRIX_INDEX(row, previous, root_multiplicity)];
        }
      }
    }
    norm = vector_norm(vector, dimension);
    if (!isfinite(norm) || !(norm > 0.0)) return 0;
    for (row = 0; row < dimension; ++row) {
      basis[MATRIX_INDEX(row, root, root_multiplicity)] =
          vector[row] / norm;
      if (!finite_complex(
              basis[MATRIX_INDEX(row, root, root_multiplicity)])) {
        return 0;
      }
    }
  }
  for (row = 0; row < dimension; ++row) {
    for (column = 0; column < dimension; ++column) {
      double complex value = 0.0;
      for (root = 0; root < root_multiplicity; ++root) {
        value += basis[MATRIX_INDEX(row, root, root_multiplicity)] *
                 conj(basis[MATRIX_INDEX(column, root,
                                         root_multiplicity)]);
      }
      projector[MATRIX_INDEX(row, column, dimension)] = value;
      if (!finite_complex(value)) return 0;
    }
  }
  return 1;
}

static MVMCKrylovGEVPStatus finalize_result(
    const MVMCKrylovGEVPPolicy *policy, const double complex *overlap,
    const double complex *hamiltonian, const double complex *squared,
    double complex *coefficient, MVMCKrylovGEVPResult *result) {
  double complex overlap_value;
  double complex energy_numerator;
  double complex squared_numerator;
  double complex residual[MVMC_KRYLOV_GEVP_MAX_DIMENSION] = {0.0};
  double coefficient_norm;
  double residual_norm;
  double matrix_norm_scale;
  double scaled_denominator;
  double absolute_energy;
  double variance_scale;
  double maximum = -1.0;
  double complex phase;
  int pivot = 0;
  int row;
  int column;
  coefficient_norm = vector_norm(coefficient, result->dimension);
  if (!isfinite(coefficient_norm) || !(coefficient_norm > 0.0)) {
    return MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
  }
  for (row = 0; row < result->dimension; ++row) {
    coefficient[row] /= coefficient_norm;
    if (!finite_complex(coefficient[row])) {
      return MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
    }
    if (cabs(coefficient[row]) > maximum) {
      maximum = cabs(coefficient[row]);
      pivot = row;
    }
  }
  if (!(maximum > 0.0) || !isfinite(maximum)) {
    return MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
  }
  phase = coefficient[pivot] / maximum;
  if (!finite_complex(phase) || cabs(phase) == 0.0) {
    return MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
  }
  for (row = 0; row < result->dimension; ++row) {
    coefficient[row] *= conj(phase);
    result->coefficient[row] = coefficient[row];
  }
  result->coefficient[pivot] = maximum;
  result->phase_pivot = pivot;
  overlap_value = quadratic_form(result->coefficient, overlap,
                                 result->dimension);
  energy_numerator = quadratic_form(result->coefficient, hamiltonian,
                                    result->dimension);
  squared_numerator = quadratic_form(result->coefficient, squared,
                                     result->dimension);
  if (!finite_complex(overlap_value) || !finite_complex(energy_numerator) ||
      !finite_complex(squared_numerator) ||
      !(creal(overlap_value) > 0.0) ||
      fabs(cimag(overlap_value)) >
          256.0 * DBL_EPSILON * fmax(creal(overlap_value), DBL_MIN)) {
    return MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
  }
  result->normalization = creal(overlap_value);
  result->energy = creal(energy_numerator) / result->normalization;
  result->energy_squared =
      creal(squared_numerator) / result->normalization;
  if (!isfinite(result->energy) || !isfinite(result->energy_squared)) {
    return MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
  }
  for (row = 0; row < result->dimension; ++row) {
    for (column = 0; column < result->dimension; ++column) {
      residual[row] +=
          (hamiltonian[MATRIX_INDEX(row, column, result->dimension)] -
           result->energy *
               overlap[MATRIX_INDEX(row, column, result->dimension)]) *
          result->coefficient[column];
    }
  }
  residual_norm = vector_norm(residual, result->dimension);
  coefficient_norm = vector_norm(result->coefficient, result->dimension);
  matrix_norm_scale =
      fmax(result->hamiltonian_norm, result->overlap_norm);
  absolute_energy = fabs(result->energy);
  if (!isfinite(residual_norm) || !isfinite(coefficient_norm) ||
      !(coefficient_norm > 0.0) || !isfinite(matrix_norm_scale) ||
      matrix_norm_scale < 0.0 || !isfinite(absolute_energy)) {
    return MVMC_KRYLOV_GEVP_RESIDUAL_FAILURE;
  }
  if (matrix_norm_scale == 0.0) {
    result->gevp_relative_residual = residual_norm == 0.0 ? 0.0 : NAN;
  } else if (absolute_energy <= 1.0) {
    /* Evaluate r / ((||K|| + |E| ||S||) ||c||) without forming the
       potentially overflowing dimensional denominator. */
    scaled_denominator =
        result->hamiltonian_norm / matrix_norm_scale +
        absolute_energy * (result->overlap_norm / matrix_norm_scale);
    result->gevp_relative_residual = scaled_denominator > 0.0
        ? (residual_norm / matrix_norm_scale) / scaled_denominator /
              coefficient_norm
        : (residual_norm == 0.0 ? 0.0 : NAN);
  } else {
    /* Factor |E| out as well so the dimensionless sum stays at most two. */
    scaled_denominator =
        (result->hamiltonian_norm / matrix_norm_scale) / absolute_energy +
        result->overlap_norm / matrix_norm_scale;
    result->gevp_relative_residual =
        (residual_norm / absolute_energy) / matrix_norm_scale /
        scaled_denominator / coefficient_norm;
  }
  if (!isfinite(result->gevp_relative_residual) ||
      result->gevp_relative_residual >
          policy->maximum_gevp_relative_residual) {
    return MVMC_KRYLOV_GEVP_RESIDUAL_FAILURE;
  }
  result->variance =
      result->energy_squared - result->energy * result->energy;
  variance_scale =
      fmax(fabs(result->energy_squared), result->energy * result->energy);
  if (!isfinite(result->variance) || !isfinite(variance_scale)) {
    return MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
  }
  if (result->variance < 0.0) {
    if (result->variance >=
        -policy->negative_variance_relative_tolerance * variance_scale) {
      result->variance = 0.0;
      result->variance_clamped = 1;
    } else {
      return MVMC_KRYLOV_GEVP_NEGATIVE_VARIANCE;
    }
  }
  return MVMC_KRYLOV_GEVP_OK;
}

MVMCKrylovGEVPStatus mvmc_krylov_gevp_default_policy(
    double rank_relative_cutoff, MVMCKrylovGEVPPolicy *policy) {
  MVMCKrylovGEVPPolicy candidate;
  if (policy == NULL || !isfinite(rank_relative_cutoff) ||
      rank_relative_cutoff <= 0.0 || rank_relative_cutoff >= 1.0) {
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.policy_version = MVMC_KRYLOV_GEVP_POLICY_VERSION;
  candidate.rank_relative_cutoff = rank_relative_cutoff;
  candidate.overlap_negative_relative_tolerance = 1.0e-10;
  candidate.hamiltonian_squared_negative_relative_tolerance = 1.0e-10;
  candidate.maximum_input_antihermitian_effect = 0.25;
  candidate.maximum_gevp_relative_residual = 1.0e-10;
  candidate.negative_variance_relative_tolerance = 1.0e-10;
  candidate.degenerate_root_gap_relative_threshold = 1.0e-10;
  *policy = candidate;
  return MVMC_KRYLOV_GEVP_OK;
}

MVMCKrylovGEVPStatus mvmc_krylov_gevp_solve_complex_packed(
    const MVMCKrylovGEVPPolicy *policy, int dimension,
    const double complex *overlap_upper,
    const double complex *hamiltonian_forward_upper,
    const double complex *hamiltonian_reverse_upper,
    const double complex *hamiltonian_squared_upper, size_t upper_count,
    MVMCKrylovGEVPResult *result) {
  MVMCKrylovGEVPResult candidate;
  double complex overlap[9] = {0.0};
  double complex hamiltonian[9] = {0.0};
  double complex squared[9] = {0.0};
  double complex equilibrated_overlap[9] = {0.0};
  double complex equilibrated_hamiltonian[9] = {0.0};
  double complex equilibrated_squared[9] = {0.0};
  double complex overlap_eigenvectors[9] = {0.0};
  double complex squared_eigenvectors[9] = {0.0};
  double complex whitening[9] = {0.0};
  double complex projected[9] = {0.0};
  double complex projected_eigenvectors[9] = {0.0};
  double complex coefficient[MVMC_KRYLOV_GEVP_MAX_DIMENSION] = {0.0};
  double inverse_diagonal[MVMC_KRYLOV_GEVP_MAX_DIMENSION] = {0.0};
  double projected_eigenvalue[MVMC_KRYLOV_GEVP_MAX_DIMENSION] = {0.0};
  double largest_overlap;
  double smallest_retained = INFINITY;
  double squared_scale = 0.0;
  size_t input_bytes;
  MVMCKrylovGEVPStatus status;
  int retained_index[MVMC_KRYLOV_GEVP_MAX_DIMENSION] = {0};
  int retained_rank = 0;
  int row;
  int column;
  int index;
  int left;
  int right;

  if (result == NULL || !valid_dimension(dimension) ||
      upper_count != required_upper_count(dimension) ||
      overlap_upper == NULL || hamiltonian_forward_upper == NULL ||
      hamiltonian_reverse_upper == NULL ||
      hamiltonian_squared_upper == NULL || !valid_policy(policy)) {
    if (result != NULL) {
      invalidate_result(MVMC_KRYLOV_GEVP_INVALID_ARGUMENT, dimension,
                        result);
    }
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  if (upper_count > SIZE_MAX / sizeof(*overlap_upper)) {
    invalidate_result(MVMC_KRYLOV_GEVP_INVALID_ARGUMENT, dimension, result);
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  input_bytes = upper_count * sizeof(*overlap_upper);
  if (!storage_is_disjoint(overlap_upper, hamiltonian_forward_upper,
                           hamiltonian_reverse_upper,
                           hamiltonian_squared_upper, input_bytes, policy,
                           result)) {
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  invalidate_result(MVMC_KRYLOV_GEVP_INVALID_ARGUMENT, dimension,
                    &candidate);
  candidate.rank_relative_cutoff = policy->rank_relative_cutoff;
  status = build_hermitian_matrices(
      dimension, overlap_upper, hamiltonian_forward_upper,
      hamiltonian_reverse_upper, hamiltonian_squared_upper, overlap,
      hamiltonian, squared, &candidate);
  if (status != MVMC_KRYLOV_GEVP_OK) goto fail;
  if (candidate.overlap_diagonal_imaginary_residual >
          policy->maximum_input_antihermitian_effect ||
      candidate.hamiltonian_antihermitian_residual >
          policy->maximum_input_antihermitian_effect ||
      candidate.hamiltonian_squared_diagonal_imaginary_residual >
          policy->maximum_input_antihermitian_effect) {
    status = MVMC_KRYLOV_GEVP_ANTIHERMITIAN_INPUT;
    goto fail;
  }
  for (row = 0; row < dimension; ++row) {
    const double diagonal = creal(overlap[MATRIX_INDEX(row, row, dimension)]);
    if (!isfinite(diagonal) || !(diagonal > 0.0)) {
      status = MVMC_KRYLOV_GEVP_NONPOSITIVE_OVERLAP_DIAGONAL;
      goto fail;
    }
    inverse_diagonal[row] = 1.0 / sqrt(diagonal);
    if (!isfinite(inverse_diagonal[row]) ||
        !(inverse_diagonal[row] > 0.0)) {
      status = MVMC_KRYLOV_GEVP_NONFINITE_INPUT;
      goto fail;
    }
  }
  for (row = 0; row < dimension; ++row) {
    for (column = 0; column < dimension; ++column) {
      const size_t dense_index = MATRIX_INDEX(row, column, dimension);
      equilibrated_overlap[dense_index] =
          (inverse_diagonal[row] * overlap[dense_index]) *
          inverse_diagonal[column];
      equilibrated_hamiltonian[dense_index] =
          (inverse_diagonal[row] * hamiltonian[dense_index]) *
          inverse_diagonal[column];
      equilibrated_squared[dense_index] =
          (inverse_diagonal[row] * squared[dense_index]) *
          inverse_diagonal[column];
      if (!finite_complex(equilibrated_overlap[dense_index]) ||
          !finite_complex(equilibrated_hamiltonian[dense_index]) ||
          !finite_complex(equilibrated_squared[dense_index])) {
        status = MVMC_KRYLOV_GEVP_NONFINITE_INPUT;
        goto fail;
      }
    }
  }
  if (!hermitian_eigensolve(
          equilibrated_overlap, dimension, candidate.overlap_eigenvalue,
          overlap_eigenvectors, &candidate.lapack_info_overlap_query,
          &candidate.lapack_info_overlap)) {
    status = MVMC_KRYLOV_GEVP_LAPACK_FAILURE;
    goto fail;
  }
  largest_overlap = candidate.overlap_eigenvalue[dimension - 1];
  if (!isfinite(largest_overlap) || !(largest_overlap > 0.0)) {
    status = MVMC_KRYLOV_GEVP_INDEFINITE_OVERLAP;
    goto fail;
  }
  for (index = 0; index < dimension; ++index) {
    const double eigenvalue = candidate.overlap_eigenvalue[index];
    if (!isfinite(eigenvalue) ||
        eigenvalue <
            -policy->overlap_negative_relative_tolerance * largest_overlap) {
      status = MVMC_KRYLOV_GEVP_INDEFINITE_OVERLAP;
      goto fail;
    }
    if (eigenvalue > policy->rank_relative_cutoff * largest_overlap) {
      retained_index[retained_rank++] = index;
      smallest_retained = fmin(smallest_retained, eigenvalue);
    }
  }
  if (retained_rank == 0) {
    status = MVMC_KRYLOV_GEVP_ZERO_RETAINED_RANK;
    goto fail;
  }
  candidate.retained_rank = retained_rank;
  candidate.discarded_rank = dimension - retained_rank;
  candidate.condition_estimate = largest_overlap / smallest_retained;
  if (!isfinite(candidate.condition_estimate) ||
      candidate.condition_estimate < 1.0) {
    status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
    goto fail;
  }
  if (!hermitian_eigensolve(
          equilibrated_squared, dimension,
          candidate.hamiltonian_squared_eigenvalue, squared_eigenvectors,
          &candidate.lapack_info_hamiltonian_squared_query,
          &candidate.lapack_info_hamiltonian_squared)) {
    status = MVMC_KRYLOV_GEVP_LAPACK_FAILURE;
    goto fail;
  }
  for (index = 0; index < dimension; ++index) {
    squared_scale = fmax(
        squared_scale, fabs(candidate.hamiltonian_squared_eigenvalue[index]));
  }
  if (squared_scale > 0.0 &&
      candidate.hamiltonian_squared_eigenvalue[0] <
          -policy->hamiltonian_squared_negative_relative_tolerance *
              squared_scale) {
    status = MVMC_KRYLOV_GEVP_INDEFINITE_HAMILTONIAN_SQUARED;
    goto fail;
  }
  for (row = 0; row < dimension; ++row) {
    for (index = 0; index < retained_rank; ++index) {
      const int eigen_index = retained_index[index];
      whitening[MATRIX_INDEX(row, index, retained_rank)] =
          overlap_eigenvectors[row + eigen_index * dimension] /
          sqrt(candidate.overlap_eigenvalue[eigen_index]);
    }
  }
  for (left = 0; left < retained_rank; ++left) {
    for (right = 0; right < retained_rank; ++right) {
      double complex value = 0.0;
      for (row = 0; row < dimension; ++row) {
        for (column = 0; column < dimension; ++column) {
          value += conj(whitening[MATRIX_INDEX(row, left, retained_rank)]) *
                   equilibrated_hamiltonian[
                       MATRIX_INDEX(row, column, dimension)] *
                   whitening[MATRIX_INDEX(column, right, retained_rank)];
        }
      }
      projected[MATRIX_INDEX(left, right, retained_rank)] = value;
    }
  }
  for (left = 0; left < retained_rank; ++left) {
    projected[MATRIX_INDEX(left, left, retained_rank)] =
        creal(projected[MATRIX_INDEX(left, left, retained_rank)]);
    for (right = left + 1; right < retained_rank; ++right) {
      const double complex value =
          0.5 * projected[MATRIX_INDEX(left, right, retained_rank)] +
          0.5 * conj(projected[MATRIX_INDEX(right, left, retained_rank)]);
      projected[MATRIX_INDEX(left, right, retained_rank)] = value;
      projected[MATRIX_INDEX(right, left, retained_rank)] = conj(value);
    }
  }
  if (!hermitian_eigensolve(
          projected, retained_rank, projected_eigenvalue,
          projected_eigenvectors, &candidate.lapack_info_projected_query,
          &candidate.lapack_info_projected)) {
    status = MVMC_KRYLOV_GEVP_LAPACK_FAILURE;
    goto fail;
  }
  candidate.projected_energy = projected_eigenvalue[0];
  candidate.root_gap = retained_rank > 1
                           ? projected_eigenvalue[1] - projected_eigenvalue[0]
                           : INFINITY;
  candidate.root_multiplicity = 1;
  while (candidate.root_multiplicity < retained_rank &&
         projected_eigenvalue[candidate.root_multiplicity] -
                 projected_eigenvalue[0] <=
             policy->degenerate_root_gap_relative_threshold *
                 fmax(1.0, fabs(projected_eigenvalue[0]))) {
    ++candidate.root_multiplicity;
  }
  candidate.coefficient_comparison_uses_projector =
      candidate.root_multiplicity > 1;
  if (!build_root_subspace_projector(
          dimension, retained_rank, candidate.root_multiplicity,
          inverse_diagonal, whitening, projected_eigenvectors,
          candidate.root_subspace_projector)) {
    status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
    goto fail;
  }
  for (row = 0; row < dimension; ++row) {
    double complex scaled_coefficient = 0.0;
    for (index = 0; index < retained_rank; ++index) {
      scaled_coefficient +=
          whitening[MATRIX_INDEX(row, index, retained_rank)] *
          projected_eigenvectors[index];
    }
    coefficient[row] = inverse_diagonal[row] * scaled_coefficient;
  }
  status = finalize_result(policy, overlap, hamiltonian, squared, coefficient,
                           &candidate);
  if (status != MVMC_KRYLOV_GEVP_OK) goto fail;
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_GEVP_OK;
  *result = candidate;
  return MVMC_KRYLOV_GEVP_OK;

fail:
  candidate.valid = 0;
  candidate.status = status;
  *result = candidate;
  return status;
}

MVMCKrylovGEVPStatus mvmc_krylov_gevp_solve_real_packed(
    const MVMCKrylovGEVPPolicy *policy, int dimension,
    const double *overlap_upper, const double *hamiltonian_forward_upper,
    const double *hamiltonian_reverse_upper,
    const double *hamiltonian_squared_upper, size_t upper_count,
    MVMCKrylovGEVPResult *result) {
  double complex overlap[6];
  double complex forward[6];
  double complex reverse[6];
  double complex squared[6];
  MVMCKrylovGEVPResult candidate;
  MVMCKrylovGEVPStatus status;
  size_t input_bytes;
  size_t entry;
  int index;
  if (result == NULL || !valid_dimension(dimension) ||
      upper_count != required_upper_count(dimension) ||
      overlap_upper == NULL || hamiltonian_forward_upper == NULL ||
      hamiltonian_reverse_upper == NULL ||
      hamiltonian_squared_upper == NULL || !valid_policy(policy)) {
    if (result != NULL) {
      invalidate_result(MVMC_KRYLOV_GEVP_INVALID_ARGUMENT, dimension,
                        result);
    }
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  if (upper_count > SIZE_MAX / sizeof(*overlap_upper)) {
    invalidate_result(MVMC_KRYLOV_GEVP_INVALID_ARGUMENT, dimension, result);
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  input_bytes = upper_count * sizeof(*overlap_upper);
  if (!storage_is_disjoint(overlap_upper, hamiltonian_forward_upper,
                           hamiltonian_reverse_upper,
                           hamiltonian_squared_upper, input_bytes, policy,
                           result)) {
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  for (entry = 0; entry < upper_count; ++entry) {
    overlap[entry] = overlap_upper[entry];
    forward[entry] = hamiltonian_forward_upper[entry];
    reverse[entry] = hamiltonian_reverse_upper[entry];
    squared[entry] = hamiltonian_squared_upper[entry];
  }
  status = mvmc_krylov_gevp_solve_complex_packed(
      policy, dimension, overlap, forward, reverse, squared, upper_count,
      &candidate);
  if (status == MVMC_KRYLOV_GEVP_OK) {
    for (index = 0; index < dimension; ++index) {
      if (fabs(cimag(candidate.coefficient[index])) >
          256.0 * DBL_EPSILON) {
        candidate.valid = 0;
        candidate.status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
        status = candidate.status;
        break;
      }
      candidate.coefficient[index] = creal(candidate.coefficient[index]);
    }
    for (index = 0;
         status == MVMC_KRYLOV_GEVP_OK && index < dimension * dimension;
         ++index) {
      if (fabs(cimag(candidate.root_subspace_projector[index])) >
          256.0 * DBL_EPSILON) {
        candidate.valid = 0;
        candidate.status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
        status = candidate.status;
        break;
      }
      candidate.root_subspace_projector[index] =
          creal(candidate.root_subspace_projector[index]);
    }
  }
  *result = candidate;
  return status;
}

const char *mvmc_krylov_gevp_status_string(MVMCKrylovGEVPStatus status) {
  switch (status) {
    case MVMC_KRYLOV_GEVP_OK:
      return "ok";
    case MVMC_KRYLOV_GEVP_INVALID_ARGUMENT:
      return "invalid argument";
    case MVMC_KRYLOV_GEVP_NONFINITE_INPUT:
      return "nonfinite input";
    case MVMC_KRYLOV_GEVP_ANTIHERMITIAN_INPUT:
      return "anti-Hermitian input exceeds effect-size guard";
    case MVMC_KRYLOV_GEVP_NONPOSITIVE_OVERLAP_DIAGONAL:
      return "overlap diagonal is not finite and positive";
    case MVMC_KRYLOV_GEVP_INDEFINITE_OVERLAP:
      return "overlap has a significant negative eigenvalue";
    case MVMC_KRYLOV_GEVP_ZERO_RETAINED_RANK:
      return "relative overlap cutoff retained no direction";
    case MVMC_KRYLOV_GEVP_INDEFINITE_HAMILTONIAN_SQUARED:
      return "Hamiltonian-squared matrix is significantly indefinite";
    case MVMC_KRYLOV_GEVP_LAPACK_FAILURE:
      return "Hermitian eigensolver failure";
    case MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE:
      return "coefficient normalization or Rayleigh reevaluation failure";
    case MVMC_KRYLOV_GEVP_RESIDUAL_FAILURE:
      return "full-basis scaled relative residual exceeded policy";
    case MVMC_KRYLOV_GEVP_NEGATIVE_VARIANCE:
      return "variance is negative beyond the relative roundoff tolerance";
  }
  return "unknown";
}
