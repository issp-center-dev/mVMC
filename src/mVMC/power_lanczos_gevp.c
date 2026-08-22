#include "power_lanczos_gevp.h"

#ifdef _mpi_use
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#include "blas_externs.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#define MATRIX_INDEX(row, column, dimension) \
  ((size_t)(row) * (size_t)(dimension) + (size_t)(column))

enum {
  MAX_DIMENSION = MVMC_KRYLOV_GEVP_MAX_DIMENSION,
  MAX_UPPER = 6,
  WORK_SIZE = 256
};

static int finite_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static int valid_dimension(int dimension) {
  return dimension == 2 || dimension == 3;
}

static size_t required_upper_count(int dimension) {
  return (size_t)dimension * ((size_t)dimension + 1U) / 2U;
}

static int valid_policy(const MVMCPowerLanczosGEVPPolicy *policy) {
  const double protocol_tolerance = 64.0 * DBL_EPSILON;
  const double cutoff[4] = {0x1p-48, 0x1p-40, 0x1p-32, 0x1p-24};
  return policy != NULL && policy->valid &&
         policy->policy_version == MVMC_POWER_LANCZOS_GEVP_POLICY_VERSION &&
         policy->cutoff_id >= MVMC_POWER_LANCZOS_GEVP_CUTOFF_S48 &&
         policy->cutoff_id <= MVMC_POWER_LANCZOS_GEVP_CUTOFF_S24 &&
         policy->rank_relative_cutoff == cutoff[policy->cutoff_id] &&
         policy->degenerate_root_gap_relative_threshold ==
             protocol_tolerance &&
         policy->maximum_normwise_backward_error == protocol_tolerance &&
         policy->negative_variance_relative_tolerance ==
             protocol_tolerance;
}

static void invalidate_result(MVMCKrylovGEVPStatus status, int dimension,
                              MVMCPowerLanczosGEVPResult *result) {
  int index;
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
  result->eigenspace_status = status;
  result->policy_version = MVMC_POWER_LANCZOS_GEVP_POLICY_VERSION;
  result->dimension = dimension;
  result->phase_pivot = -1;
  result->rank_relative_cutoff = NAN;
  result->condition_estimate = NAN;
  result->root_gap = NAN;
  result->normalization = NAN;
  result->energy = NAN;
  result->energy_squared = NAN;
  result->variance = NAN;
  result->raw_action_relative_residual = NAN;
  result->normwise_backward_error = NAN;
  result->relative_residual = NAN;
  result->eigenspace_relative_residual = NAN;
  result->eigenspace_normalization = NAN;
  for (index = 0; index < MAX_DIMENSION; ++index) {
    result->overlap_eigenvalue[index] = NAN;
  }
}

static int ranges_overlap(const void *left, size_t left_size,
                          const void *right, size_t right_size) {
  uintptr_t left_begin;
  uintptr_t right_begin;
  if (left == NULL || right == NULL || left_size == 0 || right_size == 0) {
    return 0;
  }
  left_begin = (uintptr_t)left;
  right_begin = (uintptr_t)right;
  if (left_begin > UINTPTR_MAX - left_size ||
      right_begin > UINTPTR_MAX - right_size) {
    return 1;
  }
  return left_begin < right_begin + right_size &&
         right_begin < left_begin + left_size;
}

static int disjoint_complex_inputs(
    const MVMCPowerLanczosGEVPPolicy *policy,
    const double complex *overlap_upper,
    const double complex *hamiltonian_forward_upper,
    const double complex *hamiltonian_reverse_upper,
    const double complex *hamiltonian_squared_upper, size_t upper_count,
    const MVMCPowerLanczosGEVPResult *result) {
  const void *inputs[4] = {overlap_upper, hamiltonian_forward_upper,
                           hamiltonian_reverse_upper,
                           hamiltonian_squared_upper};
  const size_t input_bytes = upper_count * sizeof(*overlap_upper);
  int left;
  int right;
  for (left = 0; left < 4; ++left) {
    if (ranges_overlap(inputs[left], input_bytes, policy, sizeof(*policy)) ||
        ranges_overlap(inputs[left], input_bytes, result, sizeof(*result))) {
      return 0;
    }
    for (right = left + 1; right < 4; ++right) {
      if (ranges_overlap(inputs[left], input_bytes, inputs[right],
                         input_bytes)) {
        return 0;
      }
    }
  }
  return !ranges_overlap(policy, sizeof(*policy), result, sizeof(*result));
}

static int disjoint_real_inputs(
    const MVMCPowerLanczosGEVPPolicy *policy, const double *overlap_upper,
    const double *hamiltonian_forward_upper,
    const double *hamiltonian_reverse_upper,
    const double *hamiltonian_squared_upper, size_t upper_count,
    const MVMCPowerLanczosGEVPResult *result) {
  const void *inputs[4] = {overlap_upper, hamiltonian_forward_upper,
                           hamiltonian_reverse_upper,
                           hamiltonian_squared_upper};
  const size_t input_bytes = upper_count * sizeof(*overlap_upper);
  int left;
  int right;
  for (left = 0; left < 4; ++left) {
    if (ranges_overlap(inputs[left], input_bytes, policy, sizeof(*policy)) ||
        ranges_overlap(inputs[left], input_bytes, result, sizeof(*result))) {
      return 0;
    }
    for (right = left + 1; right < 4; ++right) {
      if (ranges_overlap(inputs[left], input_bytes, inputs[right],
                         input_bytes)) {
        return 0;
      }
    }
  }
  return !ranges_overlap(policy, sizeof(*policy), result, sizeof(*result));
}

static double vector_norm(const double complex *vector, int dimension) {
  double norm = 0.0;
  int index;
  for (index = 0; index < dimension; ++index) {
    norm = hypot(norm, cabs(vector[index]));
  }
  return norm;
}

static double matrix_frobenius_norm(const double complex *matrix,
                                    int dimension) {
  double norm = 0.0;
  int row;
  int column;
  for (row = 0; row < dimension; ++row) {
    for (column = 0; column < dimension; ++column) {
      norm = hypot(norm, cabs(matrix[MATRIX_INDEX(row, column, dimension)]));
    }
  }
  return norm;
}

static void matrix_vector(const double complex *matrix,
                          const double complex *vector, int dimension,
                          double complex *product) {
  int row;
  int column;
  for (row = 0; row < dimension; ++row) {
    product[row] = 0.0;
    for (column = 0; column < dimension; ++column) {
      product[row] += matrix[MATRIX_INDEX(row, column, dimension)] *
                      vector[column];
    }
  }
}

static double complex inner_product(const double complex *left,
                                    const double complex *right,
                                    int dimension) {
  double complex value = 0.0;
  int index;
  for (index = 0; index < dimension; ++index) {
    value += conj(left[index]) * right[index];
  }
  return value;
}

static int unpack_matrices(
    int dimension, const double complex *overlap_upper,
    const double complex *hamiltonian_forward_upper,
    const double complex *hamiltonian_reverse_upper,
    const double complex *hamiltonian_squared_upper,
    double complex *overlap, double complex *hamiltonian,
    double complex *hamiltonian_squared,
    double *antihermitian_effect) {
  size_t entry = 0;
  double overlap_imaginary_norm = 0.0;
  double overlap_raw_norm = 0.0;
  double residual_norm = 0.0;
  double raw_norm = 0.0;
  double squared_imaginary_norm = 0.0;
  double squared_raw_norm = 0.0;
  int row;
  int column;
  for (row = 0; row < dimension; ++row) {
    for (column = row; column < dimension; ++column, ++entry) {
      const double complex overlap_value = overlap_upper[entry];
      const double complex squared_value =
          hamiltonian_squared_upper[entry];
      double complex hamiltonian_value =
          0.5 * hamiltonian_forward_upper[entry] +
          0.5 * conj(hamiltonian_reverse_upper[entry]);
      const double complex delta =
          hamiltonian_forward_upper[entry] -
          conj(hamiltonian_reverse_upper[entry]);
      if (!finite_complex(overlap_value) ||
          !finite_complex(hamiltonian_forward_upper[entry]) ||
          !finite_complex(hamiltonian_reverse_upper[entry]) ||
          !finite_complex(hamiltonian_value) ||
          !finite_complex(squared_value)) {
        return 0;
      }
      if (row == column) {
        overlap_imaginary_norm =
            hypot(overlap_imaginary_norm, 2.0 * cimag(overlap_value));
        squared_imaginary_norm =
            hypot(squared_imaginary_norm, 2.0 * cimag(squared_value));
      }
      overlap_raw_norm = hypot(overlap_raw_norm, cabs(overlap_value));
      squared_raw_norm = hypot(squared_raw_norm, cabs(squared_value));
      residual_norm = hypot(residual_norm, cabs(delta));
      raw_norm = hypot(raw_norm, cabs(hamiltonian_forward_upper[entry]));
      if (row == column) {
        residual_norm =
            hypot(residual_norm, 2.0 * cimag(hamiltonian_value));
      }
      if (row != column) {
        overlap_raw_norm = hypot(overlap_raw_norm, cabs(overlap_value));
        squared_raw_norm = hypot(squared_raw_norm, cabs(squared_value));
        residual_norm = hypot(residual_norm, cabs(delta));
        raw_norm = hypot(raw_norm,
                         cabs(hamiltonian_reverse_upper[entry]));
      }
      overlap[MATRIX_INDEX(row, column, dimension)] =
          row == column ? creal(overlap_value) : overlap_value;
      hamiltonian[MATRIX_INDEX(row, column, dimension)] =
          row == column ? creal(hamiltonian_value) : hamiltonian_value;
      hamiltonian_squared[MATRIX_INDEX(row, column, dimension)] =
          row == column ? creal(squared_value) : squared_value;
      if (row != column) {
        overlap[MATRIX_INDEX(column, row, dimension)] = conj(overlap_value);
        hamiltonian[MATRIX_INDEX(column, row, dimension)] =
            conj(hamiltonian_value);
        hamiltonian_squared[MATRIX_INDEX(column, row, dimension)] =
            conj(squared_value);
      }
    }
  }
  *antihermitian_effect = fmax(
      raw_norm > 0.0 ? residual_norm / raw_norm : residual_norm,
      fmax(overlap_raw_norm > 0.0
               ? overlap_imaginary_norm / overlap_raw_norm
               : overlap_imaginary_norm,
           squared_raw_norm > 0.0
               ? squared_imaginary_norm / squared_raw_norm
               : squared_imaginary_norm));
  if (!isfinite(*antihermitian_effect)) return 0;
  return 1;
}

static MVMCKrylovGEVPStatus recompute_returned_coefficient_diagnostics(
    const MVMCPowerLanczosGEVPPolicy *policy, int dimension,
    const double complex *overlap_upper,
    const double complex *hamiltonian_forward_upper,
    const double complex *hamiltonian_reverse_upper,
    const double complex *hamiltonian_squared_upper,
    MVMCPowerLanczosGEVPResult *result) {
  double complex overlap[MAX_DIMENSION * MAX_DIMENSION] = {0.0};
  double complex hamiltonian[MAX_DIMENSION * MAX_DIMENSION] = {0.0};
  double complex squared[MAX_DIMENSION * MAX_DIMENSION] = {0.0};
  double complex overlap_action[MAX_DIMENSION] = {0.0};
  double complex hamiltonian_action[MAX_DIMENSION] = {0.0};
  double complex squared_action[MAX_DIMENSION] = {0.0};
  double complex residual[MAX_DIMENSION] = {0.0};
  double complex normalization_value;
  double complex energy_value;
  double complex squared_value;
  double antihermitian_effect = 0.0;
  double normalization;
  double energy;
  double energy_squared;
  double residual_norm;
  double denominator;
  double backward_denominator;
  double variance;
  double variance_scale;
  double maximum = 0.0;
  int pivot = 0;
  int index;
  if (!unpack_matrices(dimension, overlap_upper,
                       hamiltonian_forward_upper,
                       hamiltonian_reverse_upper,
                       hamiltonian_squared_upper, overlap, hamiltonian,
                       squared, &antihermitian_effect)) {
    result->valid = 0;
    result->status = MVMC_KRYLOV_GEVP_NONFINITE_INPUT;
    return result->status;
  }
  matrix_vector(overlap, result->coefficient, dimension, overlap_action);
  normalization_value = inner_product(
      result->coefficient, overlap_action, dimension);
  if (!finite_complex(normalization_value) ||
      !(creal(normalization_value) > 0.0) ||
      fabs(cimag(normalization_value)) >
          1024.0 * DBL_EPSILON *
              fmax(1.0, fabs(creal(normalization_value)))) {
    result->valid = 0;
    result->status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
    return result->status;
  }
  normalization = sqrt(creal(normalization_value));
  for (index = 0; index < dimension; ++index) {
    result->coefficient[index] /= normalization;
    if (!finite_complex(result->coefficient[index])) {
      result->valid = 0;
      result->status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
      return result->status;
    }
    if (cabs(result->coefficient[index]) > maximum) {
      maximum = cabs(result->coefficient[index]);
      pivot = index;
    }
  }
  if (!(maximum > 0.0) || cimag(result->coefficient[pivot]) != 0.0 ||
      !(creal(result->coefficient[pivot]) > 0.0)) {
    result->valid = 0;
    result->status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
    return result->status;
  }
  result->phase_pivot = pivot;
  matrix_vector(overlap, result->coefficient, dimension, overlap_action);
  matrix_vector(hamiltonian, result->coefficient, dimension,
                hamiltonian_action);
  matrix_vector(squared, result->coefficient, dimension, squared_action);
  normalization_value = inner_product(
      result->coefficient, overlap_action, dimension);
  energy_value = inner_product(
      result->coefficient, hamiltonian_action, dimension);
  squared_value = inner_product(
      result->coefficient, squared_action, dimension);
  if (!finite_complex(normalization_value) || !finite_complex(energy_value) ||
      !finite_complex(squared_value) ||
      fabs(cimag(normalization_value)) >
          1024.0 * DBL_EPSILON *
              fmax(1.0, fabs(creal(normalization_value))) ||
      fabs(cimag(energy_value)) >
          1024.0 * DBL_EPSILON * fmax(1.0, fabs(creal(energy_value))) ||
      fabs(cimag(squared_value)) >
          1024.0 * DBL_EPSILON * fmax(1.0, fabs(creal(squared_value)))) {
    result->valid = 0;
    result->status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
    return result->status;
  }
  normalization = creal(normalization_value);
  if (!(normalization > 0.0) ||
      fabs(normalization - 1.0) > 4096.0 * DBL_EPSILON) {
    result->valid = 0;
    result->status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
    return result->status;
  }
  energy = creal(energy_value) / normalization;
  energy_squared = creal(squared_value) / normalization;
  for (index = 0; index < dimension; ++index) {
    residual[index] = hamiltonian_action[index] -
                      energy * overlap_action[index];
  }
  residual_norm = vector_norm(residual, dimension);
  denominator = fmax(
      1.0, fmax(vector_norm(hamiltonian_action, dimension),
                fabs(energy) * vector_norm(overlap_action, dimension)));
  backward_denominator = fmax(
      1.0,
      fmax(matrix_frobenius_norm(hamiltonian, dimension) *
               vector_norm(result->coefficient, dimension),
           fabs(energy) * matrix_frobenius_norm(overlap, dimension) *
               vector_norm(result->coefficient, dimension)));
  result->raw_action_relative_residual = residual_norm / denominator;
  result->normwise_backward_error = residual_norm / backward_denominator;
  result->relative_residual = result->normwise_backward_error;
  if (!isfinite(result->raw_action_relative_residual) ||
      !isfinite(result->normwise_backward_error) ||
      result->normwise_backward_error >
          policy->maximum_normwise_backward_error) {
    result->valid = 0;
    result->status = MVMC_KRYLOV_GEVP_RESIDUAL_FAILURE;
    return result->status;
  }
  variance = energy_squared - energy * energy;
  variance_scale = fmax(1.0, fmax(fabs(energy_squared), energy * energy));
  if (!isfinite(variance) || !isfinite(variance_scale)) {
    result->valid = 0;
    result->status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
    return result->status;
  }
  result->observables_valid =
      variance >= -policy->negative_variance_relative_tolerance *
                      variance_scale;
  result->variance_clamped = 0;
  result->normalization = normalization;
  result->energy = energy;
  result->energy_squared = energy_squared;
  result->variance = variance;
  result->eigenspace_status = MVMC_KRYLOV_GEVP_OK;
  result->eigenspace_relative_residual = result->normwise_backward_error;
  result->eigenspace_normalization = result->normalization;
  result->valid = 1;
  result->status = MVMC_KRYLOV_GEVP_OK;
  return result->status;
}

static int hermitian_eigensolve(const double complex *matrix, int dimension,
                                double *eigenvalue,
                                double complex *eigenvector) {
  double complex query_matrix[MAX_DIMENSION * MAX_DIMENSION];
  double complex work[WORK_SIZE];
  double rwork[3 * MAX_DIMENSION - 2];
  double complex query;
  char jobz = 'V';
  char uplo = 'U';
  int n = dimension;
  int lda = dimension;
  int lwork = -1;
  int info = 0;
  int row;
  int column;
  if (matrix == NULL || eigenvalue == NULL || eigenvector == NULL ||
      dimension < 1 || dimension > MAX_DIMENSION) {
    return 0;
  }
  for (row = 0; row < dimension; ++row) {
    for (column = 0; column < dimension; ++column) {
      query_matrix[row + column * dimension] =
          matrix[MATRIX_INDEX(row, column, dimension)];
    }
  }
  zheev_(&jobz, &uplo, &n, query_matrix, &lda, eigenvalue, &query, &lwork,
         rwork, &info);
  if (info != 0 || !finite_complex(query) || creal(query) < 1.0 ||
      creal(query) > (double)WORK_SIZE) {
    return 0;
  }
  for (row = 0; row < dimension; ++row) {
    for (column = 0; column < dimension; ++column) {
      eigenvector[row + column * dimension] =
          matrix[MATRIX_INDEX(row, column, dimension)];
    }
  }
  lwork = WORK_SIZE;
  zheev_(&jobz, &uplo, &n, eigenvector, &lda, eigenvalue, work, &lwork,
         rwork, &info);
  return info == 0;
}

static int canonical_coefficient(
    const double complex *overlap, const double complex *root_vectors,
    int dimension, int multiplicity, double complex *coefficient,
    int *phase_pivot) {
  int canonical_index;
  for (canonical_index = 0; canonical_index < dimension;
       ++canonical_index) {
    double complex projected[MAX_DIMENSION] = {0.0};
    double complex overlap_projected[MAX_DIMENSION] = {0.0};
    double complex norm_value;
    double norm;
    int root;
    int row;
    for (root = 0; root < multiplicity; ++root) {
      double complex weight = 0.0;
      for (row = 0; row < dimension; ++row) {
        weight +=
            conj(root_vectors[MATRIX_INDEX(row, root, multiplicity)]) *
            overlap[MATRIX_INDEX(row, canonical_index, dimension)];
      }
      for (row = 0; row < dimension; ++row) {
        projected[row] +=
            root_vectors[MATRIX_INDEX(row, root, multiplicity)] * weight;
      }
    }
    matrix_vector(overlap, projected, dimension, overlap_projected);
    norm_value = inner_product(projected, overlap_projected, dimension);
    if (!finite_complex(norm_value) || !(creal(norm_value) > 0.0) ||
        fabs(cimag(norm_value)) >
            1024.0 * DBL_EPSILON * fmax(1.0, creal(norm_value))) {
      continue;
    }
    norm = sqrt(creal(norm_value));
    if (!(norm > 64.0 * DBL_EPSILON) || !isfinite(norm)) continue;
    for (row = 0; row < dimension; ++row) {
      coefficient[row] = projected[row] / norm;
    }
    break;
  }
  if (canonical_index == dimension) return 0;
  {
    double maximum = -1.0;
    double complex phase;
    int pivot = 0;
    int row;
    for (row = 0; row < dimension; ++row) {
      const double magnitude = cabs(coefficient[row]);
      if (magnitude > maximum) {
        maximum = magnitude;
        pivot = row;
      }
    }
    if (!(maximum > 0.0) || !isfinite(maximum)) return 0;
    phase = coefficient[pivot] / maximum;
    if (!finite_complex(phase) || cabs(phase) == 0.0) return 0;
    for (row = 0; row < dimension; ++row) {
      coefficient[row] *= conj(phase);
      if (!finite_complex(coefficient[row])) return 0;
    }
    coefficient[pivot] = maximum;
    *phase_pivot = pivot;
  }
  return 1;
}

static int build_root_projector(const double complex *root_vectors,
                                int dimension, int multiplicity,
                                double complex *projector) {
  double complex basis[MAX_DIMENSION * MAX_DIMENSION] = {0.0};
  int root;
  int row;
  int column;
  for (root = 0; root < multiplicity; ++root) {
    double complex vector[MAX_DIMENSION] = {0.0};
    double norm;
    int pass;
    int previous;
    for (row = 0; row < dimension; ++row) {
      vector[row] = root_vectors[MATRIX_INDEX(row, root, multiplicity)];
    }
    for (pass = 0; pass < 2; ++pass) {
      for (previous = 0; previous < root; ++previous) {
        double complex value = 0.0;
        for (row = 0; row < dimension; ++row) {
          value += conj(basis[MATRIX_INDEX(
                            row, previous, multiplicity)]) * vector[row];
        }
        for (row = 0; row < dimension; ++row) {
          vector[row] -= value *
              basis[MATRIX_INDEX(row, previous, multiplicity)];
        }
      }
    }
    norm = vector_norm(vector, dimension);
    if (!isfinite(norm) || !(norm > 0.0)) return 0;
    for (row = 0; row < dimension; ++row) {
      basis[MATRIX_INDEX(row, root, multiplicity)] = vector[row] / norm;
    }
  }
  for (row = 0; row < dimension; ++row) {
    for (column = 0; column < dimension; ++column) {
      double complex value = 0.0;
      for (root = 0; root < multiplicity; ++root) {
        value += basis[MATRIX_INDEX(row, root, multiplicity)] *
                 conj(basis[MATRIX_INDEX(column, root, multiplicity)]);
      }
      projector[MATRIX_INDEX(row, column, dimension)] = value;
    }
  }
  return 1;
}

static MVMCKrylovGEVPStatus solve_complex_impl(
    const MVMCPowerLanczosGEVPPolicy *policy, int dimension,
    const double complex *overlap_upper,
    const double complex *hamiltonian_forward_upper,
    const double complex *hamiltonian_reverse_upper,
    const double complex *hamiltonian_squared_upper, size_t upper_count,
    MVMCPowerLanczosGEVPResult *result) {
  double complex overlap[MAX_DIMENSION * MAX_DIMENSION] = {0.0};
  double complex hamiltonian[MAX_DIMENSION * MAX_DIMENSION] = {0.0};
  double complex squared[MAX_DIMENSION * MAX_DIMENSION] = {0.0};
  double overlap_eigenvalue[MAX_DIMENSION] = {0.0};
  double squared_eigenvalue[MAX_DIMENSION] = {0.0};
  double root[MAX_DIMENSION] = {0.0};
  double complex overlap_eigenvector[MAX_DIMENSION * MAX_DIMENSION] = {0.0};
  double complex squared_eigenvector[MAX_DIMENSION * MAX_DIMENSION] = {0.0};
  double complex projected_eigenvector[MAX_DIMENSION * MAX_DIMENSION] = {0.0};
  double complex whitening[MAX_DIMENSION * MAX_DIMENSION] = {0.0};
  double complex projected[MAX_DIMENSION * MAX_DIMENSION] = {0.0};
  double complex root_vectors[MAX_DIMENSION * MAX_DIMENSION] = {0.0};
  double complex coefficient[MAX_DIMENSION] = {0.0};
  double complex overlap_action[MAX_DIMENSION] = {0.0};
  double complex hamiltonian_action[MAX_DIMENSION] = {0.0};
  double complex squared_action[MAX_DIMENSION] = {0.0};
  double complex residual[MAX_DIMENSION] = {0.0};
  double complex normalization_value;
  double complex energy_value;
  double complex squared_value;
  double normalization;
  double energy;
  double energy_squared;
  double variance;
  double variance_scale;
  double denominator;
  double backward_denominator;
  double antihermitian_effect = NAN;
  double largest_overlap;
  double smallest_retained = INFINITY;
  double squared_scale = 0.0;
  int retained_index[MAX_DIMENSION] = {0};
  int retained_rank = 0;
  int root_multiplicity;
  int row;
  int column;
  int index;

  invalidate_result(MVMC_KRYLOV_GEVP_INVALID_ARGUMENT, dimension, result);
  if (!valid_policy(policy) || !valid_dimension(dimension) ||
      upper_count != required_upper_count(dimension) ||
      overlap_upper == NULL || hamiltonian_forward_upper == NULL ||
      hamiltonian_reverse_upper == NULL ||
      hamiltonian_squared_upper == NULL || result == NULL) {
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  result->rank_relative_cutoff = policy->rank_relative_cutoff;
  if (!unpack_matrices(dimension, overlap_upper,
                       hamiltonian_forward_upper,
                       hamiltonian_reverse_upper,
                       hamiltonian_squared_upper, overlap, hamiltonian,
                       squared, &antihermitian_effect)) {
    result->status = MVMC_KRYLOV_GEVP_NONFINITE_INPUT;
    return result->status;
  }
  if (antihermitian_effect > 0.25) {
    result->status = MVMC_KRYLOV_GEVP_ANTIHERMITIAN_INPUT;
    return result->status;
  }
  if (!hermitian_eigensolve(overlap, dimension, overlap_eigenvalue,
                            overlap_eigenvector)) {
    result->status = MVMC_KRYLOV_GEVP_LAPACK_FAILURE;
    return result->status;
  }
  largest_overlap = overlap_eigenvalue[dimension - 1];
  if (!isfinite(largest_overlap) || !(largest_overlap > 0.0)) {
    result->status = MVMC_KRYLOV_GEVP_ZERO_RETAINED_RANK;
    return result->status;
  }
  for (index = 0; index < dimension; ++index) {
    if (!isfinite(overlap_eigenvalue[index]) ||
        overlap_eigenvalue[index] <
            -64.0 * DBL_EPSILON * largest_overlap) {
      result->status = MVMC_KRYLOV_GEVP_INDEFINITE_OVERLAP;
      return result->status;
    }
    result->overlap_eigenvalue[index] = overlap_eigenvalue[index];
  }
  for (index = dimension - 1; index >= 0; --index) {
    if (overlap_eigenvalue[index] >
        policy->rank_relative_cutoff * largest_overlap) {
      retained_index[retained_rank++] = index;
      smallest_retained = fmin(smallest_retained,
                               overlap_eigenvalue[index]);
    }
  }
  if (retained_rank == 0) {
    result->status = MVMC_KRYLOV_GEVP_ZERO_RETAINED_RANK;
    return result->status;
  }
  result->retained_rank = retained_rank;
  result->discarded_rank = dimension - retained_rank;
  result->condition_estimate = largest_overlap / smallest_retained;
  if (!isfinite(result->condition_estimate) ||
      result->condition_estimate < 1.0) {
    result->status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
    return result->status;
  }
  if (!hermitian_eigensolve(squared, dimension, squared_eigenvalue,
                            squared_eigenvector)) {
    result->status = MVMC_KRYLOV_GEVP_LAPACK_FAILURE;
    return result->status;
  }
  for (index = 0; index < dimension; ++index) {
    squared_scale = fmax(squared_scale, fabs(squared_eigenvalue[index]));
  }
  if (squared_scale > 0.0 &&
      squared_eigenvalue[0] < -64.0 * DBL_EPSILON * squared_scale) {
    result->status = MVMC_KRYLOV_GEVP_INDEFINITE_HAMILTONIAN_SQUARED;
    return result->status;
  }
  for (row = 0; row < dimension; ++row) {
    for (index = 0; index < retained_rank; ++index) {
      const int eigen_index = retained_index[index];
      whitening[MATRIX_INDEX(row, index, retained_rank)] =
          overlap_eigenvector[row + eigen_index * dimension] /
          sqrt(overlap_eigenvalue[eigen_index]);
    }
  }
  for (row = 0; row < retained_rank; ++row) {
    for (column = 0; column < retained_rank; ++column) {
      double complex value = 0.0;
      int physical_row;
      int physical_column;
      for (physical_row = 0; physical_row < dimension; ++physical_row) {
        for (physical_column = 0; physical_column < dimension;
             ++physical_column) {
          value +=
              conj(whitening[MATRIX_INDEX(
                  physical_row, row, retained_rank)]) *
              hamiltonian[MATRIX_INDEX(
                  physical_row, physical_column, dimension)] *
              whitening[MATRIX_INDEX(
                  physical_column, column, retained_rank)];
        }
      }
      projected[MATRIX_INDEX(row, column, retained_rank)] = value;
    }
  }
  for (row = 0; row < retained_rank; ++row) {
    projected[MATRIX_INDEX(row, row, retained_rank)] =
        creal(projected[MATRIX_INDEX(row, row, retained_rank)]);
    for (column = row + 1; column < retained_rank; ++column) {
      const double complex value =
          0.5 * projected[MATRIX_INDEX(row, column, retained_rank)] +
          0.5 * conj(projected[MATRIX_INDEX(
                    column, row, retained_rank)]);
      projected[MATRIX_INDEX(row, column, retained_rank)] = value;
      projected[MATRIX_INDEX(column, row, retained_rank)] = conj(value);
    }
  }
  if (!hermitian_eigensolve(projected, retained_rank, root,
                            projected_eigenvector)) {
    result->status = MVMC_KRYLOV_GEVP_LAPACK_FAILURE;
    return result->status;
  }
  root_multiplicity = 1;
  while (root_multiplicity < retained_rank &&
         root[root_multiplicity] - root[0] <=
             policy->degenerate_root_gap_relative_threshold *
                 fmax(1.0, fabs(root[0]))) {
    ++root_multiplicity;
  }
  result->root_multiplicity = root_multiplicity;
  result->root_gap = retained_rank > 1 ? root[1] - root[0] : 0.0;
  for (row = 0; row < dimension; ++row) {
    int root_index;
    for (root_index = 0; root_index < root_multiplicity; ++root_index) {
      for (index = 0; index < retained_rank; ++index) {
        root_vectors[MATRIX_INDEX(
            row, root_index, root_multiplicity)] +=
            whitening[MATRIX_INDEX(row, index, retained_rank)] *
            projected_eigenvector[index +
                                  root_index * retained_rank];
      }
    }
  }
  if (!canonical_coefficient(overlap, root_vectors, dimension,
                             root_multiplicity, coefficient,
                             &result->phase_pivot) ||
      !build_root_projector(root_vectors, dimension, root_multiplicity,
                            result->root_subspace_projector)) {
    result->status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
    return result->status;
  }
  matrix_vector(overlap, coefficient, dimension, overlap_action);
  matrix_vector(hamiltonian, coefficient, dimension, hamiltonian_action);
  matrix_vector(squared, coefficient, dimension, squared_action);
  normalization_value =
      inner_product(coefficient, overlap_action, dimension);
  energy_value = inner_product(coefficient, hamiltonian_action, dimension);
  squared_value = inner_product(coefficient, squared_action, dimension);
  if (!finite_complex(normalization_value) ||
      !finite_complex(energy_value) || !finite_complex(squared_value) ||
      !(creal(normalization_value) > 0.0) ||
      fabs(cimag(normalization_value)) >
          1024.0 * DBL_EPSILON *
              fmax(1.0, fabs(creal(normalization_value))) ||
      fabs(cimag(energy_value)) >
          1024.0 * DBL_EPSILON * fmax(1.0, fabs(creal(energy_value))) ||
      fabs(cimag(squared_value)) >
          1024.0 * DBL_EPSILON * fmax(1.0, fabs(creal(squared_value)))) {
    result->status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
    return result->status;
  }
  normalization = creal(normalization_value);
  if (fabs(normalization - 1.0) > 4096.0 * DBL_EPSILON) {
    result->status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
    return result->status;
  }
  energy = creal(energy_value) / normalization;
  energy_squared = creal(squared_value) / normalization;
  if (!isfinite(energy) || !isfinite(energy_squared)) {
    result->status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
    return result->status;
  }
  for (index = 0; index < dimension; ++index) {
    residual[index] = hamiltonian_action[index] -
                      energy * overlap_action[index];
  }
  denominator = fmax(
      1.0, fmax(vector_norm(hamiltonian_action, dimension),
                fabs(energy) * vector_norm(overlap_action, dimension)));
  result->raw_action_relative_residual =
      vector_norm(residual, dimension) / denominator;
  backward_denominator = fmax(
      1.0,
      fmax(matrix_frobenius_norm(hamiltonian, dimension) *
               vector_norm(coefficient, dimension),
           fabs(energy) * matrix_frobenius_norm(overlap, dimension) *
               vector_norm(coefficient, dimension)));
  result->normwise_backward_error =
      vector_norm(residual, dimension) / backward_denominator;
  result->relative_residual = result->normwise_backward_error;
  if (!isfinite(result->raw_action_relative_residual) ||
      !isfinite(result->normwise_backward_error) ||
      result->normwise_backward_error >
          policy->maximum_normwise_backward_error) {
    result->status = MVMC_KRYLOV_GEVP_RESIDUAL_FAILURE;
    return result->status;
  }
  variance = energy_squared - energy * energy;
  variance_scale =
      fmax(1.0, fmax(fabs(energy_squared), energy * energy));
  if (!isfinite(variance) || !isfinite(variance_scale)) {
    result->status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
    return result->status;
  }
  /* C1 release variance is recomputed from physical Psi/H Psi.  This raw
   * packed-matrix difference remains diagnostic and never selects a root. */
  if (variance >=
      -policy->negative_variance_relative_tolerance * variance_scale) {
    result->observables_valid = 1;
  }
  result->normalization = normalization;
  result->energy = energy;
  result->energy_squared = energy_squared;
  result->variance = variance;
  for (index = 0; index < dimension; ++index) {
    result->coefficient[index] = coefficient[index];
  }
  result->valid = 1;
  result->status = MVMC_KRYLOV_GEVP_OK;
  result->eigenspace_status = MVMC_KRYLOV_GEVP_OK;
  result->eigenspace_relative_residual = result->normwise_backward_error;
  result->eigenspace_normalization = result->normalization;
  return MVMC_KRYLOV_GEVP_OK;
}

MVMCKrylovGEVPStatus mvmc_power_lanczos_gevp_default_policy(
    double rank_relative_cutoff, MVMCPowerLanczosGEVPPolicy *policy) {
  MVMCPowerLanczosGEVPPolicy candidate;
  int cutoff_id = -1;
  if (rank_relative_cutoff == 0x1p-48) {
    cutoff_id = MVMC_POWER_LANCZOS_GEVP_CUTOFF_S48;
  } else if (rank_relative_cutoff == 0x1p-40) {
    cutoff_id = MVMC_POWER_LANCZOS_GEVP_CUTOFF_S40;
  } else if (rank_relative_cutoff == 0x1p-32) {
    cutoff_id = MVMC_POWER_LANCZOS_GEVP_CUTOFF_S32;
  } else if (rank_relative_cutoff == 0x1p-24) {
    cutoff_id = MVMC_POWER_LANCZOS_GEVP_CUTOFF_S24;
  }
  if (policy == NULL || cutoff_id < 0) {
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.policy_version = MVMC_POWER_LANCZOS_GEVP_POLICY_VERSION;
  candidate.cutoff_id = (MVMCPowerLanczosGEVPCutoffID)cutoff_id;
  candidate.rank_relative_cutoff = rank_relative_cutoff;
  candidate.degenerate_root_gap_relative_threshold = 64.0 * DBL_EPSILON;
  candidate.maximum_normwise_backward_error = 64.0 * DBL_EPSILON;
  candidate.negative_variance_relative_tolerance = 64.0 * DBL_EPSILON;
  *policy = candidate;
  return MVMC_KRYLOV_GEVP_OK;
}

MVMCKrylovGEVPStatus mvmc_power_lanczos_gevp_solve_complex_packed(
    const MVMCPowerLanczosGEVPPolicy *policy, int dimension,
    const double complex *overlap_upper,
    const double complex *hamiltonian_forward_upper,
    const double complex *hamiltonian_reverse_upper,
    const double complex *hamiltonian_squared_upper, size_t upper_count,
    MVMCPowerLanczosGEVPResult *result) {
  if (result == NULL) {
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  if (upper_count > SIZE_MAX / sizeof(*overlap_upper)) {
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  if (!disjoint_complex_inputs(
          policy, overlap_upper, hamiltonian_forward_upper,
          hamiltonian_reverse_upper, hamiltonian_squared_upper,
          upper_count, result)) {
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  return solve_complex_impl(policy, dimension, overlap_upper,
                            hamiltonian_forward_upper,
                            hamiltonian_reverse_upper,
                            hamiltonian_squared_upper, upper_count, result);
}

MVMCKrylovGEVPStatus mvmc_power_lanczos_gevp_solve_real_packed(
    const MVMCPowerLanczosGEVPPolicy *policy, int dimension,
    const double *overlap_upper, const double *hamiltonian_forward_upper,
    const double *hamiltonian_reverse_upper,
    const double *hamiltonian_squared_upper, size_t upper_count,
    MVMCPowerLanczosGEVPResult *result) {
  double complex overlap[MAX_UPPER] = {0.0};
  double complex forward[MAX_UPPER] = {0.0};
  double complex reverse[MAX_UPPER] = {0.0};
  double complex squared[MAX_UPPER] = {0.0};
  MVMCKrylovGEVPStatus status;
  size_t index;
  if (result == NULL) {
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  if (upper_count > MAX_UPPER ||
      upper_count > SIZE_MAX / sizeof(*overlap_upper)) {
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  if (!disjoint_real_inputs(policy, overlap_upper,
                            hamiltonian_forward_upper,
                            hamiltonian_reverse_upper,
                            hamiltonian_squared_upper, upper_count, result)) {
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  if (overlap_upper == NULL || hamiltonian_forward_upper == NULL ||
      hamiltonian_reverse_upper == NULL ||
      hamiltonian_squared_upper == NULL) {
    invalidate_result(MVMC_KRYLOV_GEVP_INVALID_ARGUMENT, dimension, result);
    return MVMC_KRYLOV_GEVP_INVALID_ARGUMENT;
  }
  for (index = 0; index < upper_count; ++index) {
    overlap[index] = overlap_upper[index];
    forward[index] = hamiltonian_forward_upper[index];
    reverse[index] = hamiltonian_reverse_upper[index];
    squared[index] = hamiltonian_squared_upper[index];
  }
  status = solve_complex_impl(policy, dimension, overlap, forward, reverse,
                              squared, upper_count, result);
  if (status == MVMC_KRYLOV_GEVP_OK) {
    for (index = 0; index < (size_t)dimension; ++index) {
      if (fabs(cimag(result->coefficient[index])) >
          1024.0 * DBL_EPSILON) {
        result->valid = 0;
        result->status = MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE;
        return result->status;
      }
      result->coefficient[index] = creal(result->coefficient[index]);
    }
    status = recompute_returned_coefficient_diagnostics(
        policy, dimension, overlap, forward, reverse, squared, result);
  }
  return status;
}

const char *mvmc_power_lanczos_gevp_policy_id(void) {
  return MVMC_POWER_LANCZOS_GEVP_POLICY_ID;
}
