#include "power_lanczos_gevp.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(condition, ...)                                             \
  do {                                                                    \
    if (!(condition)) {                                                   \
      fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__);              \
      fprintf(stderr, __VA_ARGS__);                                       \
      fputc('\n', stderr);                                                \
      ++failures;                                                         \
    }                                                                     \
  } while (0)

static double s_norm(const double complex *coefficient,
                     const double complex *overlap, int dimension) {
  double complex value = 0.0;
  int row;
  int column;
  for (row = 0; row < dimension; ++row) {
    for (column = 0; column < dimension; ++column) {
      value += conj(coefficient[row]) *
               overlap[row * dimension + column] * coefficient[column];
    }
  }
  return creal(value);
}

static double test_vector_norm(const double complex *vector, int dimension) {
  double norm = 0.0;
  int index;
  for (index = 0; index < dimension; ++index) {
    norm = hypot(norm, cabs(vector[index]));
  }
  return norm;
}

static double test_matrix_norm(const double complex *matrix, int dimension) {
  double norm = 0.0;
  int row;
  int column;
  for (row = 0; row < dimension; ++row) {
    for (column = 0; column < dimension; ++column) {
      norm = hypot(norm, cabs(matrix[row * dimension + column]));
    }
  }
  return norm;
}

static int report_close(double reported, double recomputed) {
  const double error = fabs(reported - recomputed);
  return error <= 1.0e-18 ||
         error <= 2.0e-10 * fmax(fabs(recomputed), DBL_MIN);
}

static void check_returned_real_diagnostics(
    const char *label, int dimension, const double *overlap_upper,
    const double *forward_upper, const double *reverse_upper,
    const double *squared_upper,
    const MVMCPowerLanczosGEVPResult *result) {
  double complex overlap[9] = {0.0};
  double complex hamiltonian[9] = {0.0};
  double complex squared[9] = {0.0};
  double complex overlap_action[3] = {0.0};
  double complex hamiltonian_action[3] = {0.0};
  double complex squared_action[3] = {0.0};
  double complex residual[3] = {0.0};
  double complex normalization_value = 0.0;
  double complex energy_value = 0.0;
  double complex squared_value = 0.0;
  double normalization;
  double energy;
  double energy_squared;
  double residual_norm;
  double raw_denominator;
  double backward_denominator;
  double raw_residual;
  double backward_error;
  size_t entry = 0;
  int row;
  int column;
  int index;
  for (row = 0; row < dimension; ++row) {
    for (column = row; column < dimension; ++column) {
      const double s_value = overlap_upper[entry];
      const double h_value =
          0.5 * forward_upper[entry] + 0.5 * reverse_upper[entry];
      const double b_value = squared_upper[entry];
      overlap[row * dimension + column] = s_value;
      hamiltonian[row * dimension + column] = h_value;
      squared[row * dimension + column] = b_value;
      overlap[column * dimension + row] = s_value;
      hamiltonian[column * dimension + row] = h_value;
      squared[column * dimension + row] = b_value;
      ++entry;
    }
  }
  for (row = 0; row < dimension; ++row) {
    for (column = 0; column < dimension; ++column) {
      overlap_action[row] +=
          overlap[row * dimension + column] * result->coefficient[column];
      hamiltonian_action[row] +=
          hamiltonian[row * dimension + column] * result->coefficient[column];
      squared_action[row] +=
          squared[row * dimension + column] * result->coefficient[column];
    }
  }
  for (index = 0; index < dimension; ++index) {
    normalization_value +=
        conj(result->coefficient[index]) * overlap_action[index];
    energy_value +=
        conj(result->coefficient[index]) * hamiltonian_action[index];
    squared_value +=
        conj(result->coefficient[index]) * squared_action[index];
  }
  normalization = creal(normalization_value);
  energy = creal(energy_value) / normalization;
  energy_squared = creal(squared_value) / normalization;
  for (index = 0; index < dimension; ++index) {
    residual[index] =
        hamiltonian_action[index] - energy * overlap_action[index];
  }
  residual_norm = test_vector_norm(residual, dimension);
  raw_denominator = fmax(
      1.0, fmax(test_vector_norm(hamiltonian_action, dimension),
                fabs(energy) * test_vector_norm(overlap_action, dimension)));
  backward_denominator = fmax(
      1.0,
      fmax(test_matrix_norm(hamiltonian, dimension) *
               test_vector_norm(result->coefficient, dimension),
           fabs(energy) * test_matrix_norm(overlap, dimension) *
               test_vector_norm(result->coefficient, dimension)));
  raw_residual = residual_norm / raw_denominator;
  backward_error = residual_norm / backward_denominator;
  CHECK(report_close(result->normalization, normalization),
        "%s returned normalization report %.17g recomputed %.17g", label,
        result->normalization, normalization);
  CHECK(report_close(result->energy, energy),
        "%s returned energy report %.17g recomputed %.17g", label,
        result->energy, energy);
  CHECK(report_close(result->energy_squared, energy_squared),
        "%s returned M2 report %.17g recomputed %.17g", label,
        result->energy_squared, energy_squared);
  CHECK(report_close(result->variance,
                     energy_squared - energy * energy),
        "%s returned variance report %.17g recomputed %.17g", label,
        result->variance, energy_squared - energy * energy);
  CHECK(report_close(result->raw_action_relative_residual, raw_residual),
        "%s returned raw residual %.17g recomputed %.17g", label,
        result->raw_action_relative_residual, raw_residual);
  CHECK(report_close(result->normwise_backward_error, backward_error),
        "%s returned backward error %.17g recomputed %.17g", label,
        result->normwise_backward_error, backward_error);
}

static void test_policy(void) {
  MVMCPowerLanczosGEVPPolicy policy;
  memset(&policy, 0xa5, sizeof(policy));
  CHECK(mvmc_power_lanczos_gevp_default_policy(0x1p-40, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "default policy");
  CHECK(policy.valid &&
            policy.policy_version ==
                MVMC_POWER_LANCZOS_GEVP_POLICY_VERSION &&
            policy.cutoff_id == MVMC_POWER_LANCZOS_GEVP_CUTOFF_S40 &&
            policy.rank_relative_cutoff == 0x1p-40 &&
            policy.degenerate_root_gap_relative_threshold ==
                64.0 * DBL_EPSILON &&
            policy.maximum_normwise_backward_error ==
                64.0 * DBL_EPSILON &&
            policy.negative_variance_relative_tolerance ==
                64.0 * DBL_EPSILON,
        "P6 policy fields");
  CHECK(strcmp(mvmc_power_lanczos_gevp_policy_id(),
               MVMC_POWER_LANCZOS_GEVP_POLICY_ID) == 0,
        "policy id");
  CHECK(mvmc_power_lanczos_gevp_default_policy(0.0, &policy) ==
            MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "zero cutoff rejected");
  CHECK(mvmc_power_lanczos_gevp_default_policy(1.0, &policy) ==
            MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "unit cutoff rejected");
  CHECK(mvmc_power_lanczos_gevp_default_policy(NAN, &policy) ==
            MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "nonfinite cutoff rejected");
  CHECK(mvmc_power_lanczos_gevp_default_policy(0.1, &policy) ==
            MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "off-grid cutoff rejected");
  CHECK(mvmc_power_lanczos_gevp_default_policy(0x1p-30, &policy) ==
            MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "unregistered power-of-two cutoff rejected");
}

static void test_real_unique_and_scale(void) {
  const double scale[3] = {1.0, 0.5, 2.0};
  MVMCPowerLanczosGEVPPolicy policy;
  int scale_index;
  CHECK(mvmc_power_lanczos_gevp_default_policy(0x1p-40, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "real policy");
  for (scale_index = 0; scale_index < 3; ++scale_index) {
    const double factor = scale[scale_index];
    const double overlap[3] = {factor, 0.0, factor};
    const double forward[3] = {-2.0 * factor, 0.0, factor};
    const double reverse[3] = {-2.0 * factor, 0.0, factor};
    const double squared[3] = {4.0 * factor, 0.0, factor};
    const double complex dense_overlap[4] = {factor, 0.0, 0.0, factor};
    MVMCPowerLanczosGEVPResult result;
    CHECK(mvmc_power_lanczos_gevp_solve_real_packed(
              &policy, 2, overlap, forward, reverse, squared, 3,
              &result) == MVMC_KRYLOV_GEVP_OK && result.valid,
          "real unique scale %d status=%s", scale_index,
          mvmc_krylov_gevp_status_string(result.status));
    if (!result.valid) continue;
    CHECK(result.retained_rank == 2 && result.root_multiplicity == 1,
          "real ranks scale %d", scale_index);
    CHECK(fabs(s_norm(result.coefficient, dense_overlap, 2) - 1.0) <=
              1024.0 * DBL_EPSILON,
          "real S normalization scale %d %.17g", scale_index,
          result.normalization);
    CHECK(fabs(result.energy + 2.0) <= 64.0 * DBL_EPSILON &&
              fabs(result.energy_squared - 4.0) <=
                  128.0 * DBL_EPSILON &&
              result.variance == 0.0,
          "real observables scale %d E=%.17g M2=%.17g V=%.17g",
          scale_index, result.energy, result.energy_squared,
          result.variance);
    CHECK(result.relative_residual <= 64.0 * DBL_EPSILON,
          "real P6 residual scale %d %.17g", scale_index,
          result.relative_residual);
    CHECK(cimag(result.coefficient[0]) == 0.0 &&
              cimag(result.coefficient[1]) == 0.0,
          "real coefficient is real");
    check_returned_real_diagnostics(
        "dimension-two real", 2, overlap, forward, reverse, squared,
        &result);
  }
}

static void test_complex_unique(void) {
  const double complex overlap[3] = {1.0, 0.0, 1.0};
  const double complex forward[3] = {-1.0, 0.25 * I, 2.0};
  const double complex reverse[3] = {-1.0, -0.25 * I, 2.0};
  const double complex squared[3] = {1.0625, 0.25 * I, 4.0625};
  const double complex dense_overlap[4] = {1.0, 0.0, 0.0, 1.0};
  MVMCPowerLanczosGEVPPolicy policy;
  MVMCPowerLanczosGEVPResult result;
  CHECK(mvmc_power_lanczos_gevp_default_policy(0x1p-40, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "complex policy");
  CHECK(mvmc_power_lanczos_gevp_solve_complex_packed(
            &policy, 2, overlap, forward, reverse, squared, 3,
            &result) == MVMC_KRYLOV_GEVP_OK && result.valid,
        "complex unique status=%s",
        mvmc_krylov_gevp_status_string(result.status));
  if (!result.valid) return;
  CHECK(fabs(s_norm(result.coefficient, dense_overlap, 2) - 1.0) <=
            1024.0 * DBL_EPSILON,
        "complex S normalization %.17g", result.normalization);
  CHECK(result.root_multiplicity == 1 && result.phase_pivot >= 0 &&
            cimag(result.coefficient[result.phase_pivot]) == 0.0 &&
            creal(result.coefficient[result.phase_pivot]) > 0.0,
        "complex phase/root");
  CHECK(result.relative_residual <= 64.0 * DBL_EPSILON &&
            result.variance <= 1024.0 * DBL_EPSILON,
        "complex residual/variance %.17g %.17g",
        result.relative_residual, result.variance);
}

static void test_dimension_three_unique(void) {
  const double real_s[6] = {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
  const double real_kf[6] = {-3.0, 0.0, 0.0, 0.0, 0.0, 2.0};
  const double real_kr[6] = {-3.0, 0.0, 0.0, 0.0, 0.0, 2.0};
  const double real_b[6] = {9.0, 0.0, 0.0, 0.0, 0.0, 4.0};
  const double complex complex_s[6] = {
      1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
  const double complex complex_kf[6] = {
      -2.0, 0.25 * I, 0.0, 0.5, 0.0, 2.0};
  const double complex complex_kr[6] = {
      -2.0, -0.25 * I, 0.0, 0.5, 0.0, 2.0};
  const double complex complex_b[6] = {
      4.0625, -0.375 * I, 0.0, 0.3125, 0.0, 4.0};
  MVMCPowerLanczosGEVPPolicy policy;
  MVMCPowerLanczosGEVPResult result;
  CHECK(mvmc_power_lanczos_gevp_default_policy(0x1p-40, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "dimension-three policy");
  CHECK(mvmc_power_lanczos_gevp_solve_real_packed(
            &policy, 3, real_s, real_kf, real_kr, real_b, 6, &result) ==
            MVMC_KRYLOV_GEVP_OK && result.valid &&
            result.dimension == 3 && result.retained_rank == 3 &&
            result.root_multiplicity == 1 &&
            fabs(result.energy + 3.0) <= 64.0 * DBL_EPSILON &&
            result.relative_residual <= 64.0 * DBL_EPSILON,
        "dimension-three real unique status=%s",
        mvmc_krylov_gevp_status_string(result.status));
  if (result.valid) {
    check_returned_real_diagnostics(
        "dimension-three real", 3, real_s, real_kf, real_kr, real_b,
        &result);
  }
  CHECK(mvmc_power_lanczos_gevp_solve_complex_packed(
            &policy, 3, complex_s, complex_kf, complex_kr, complex_b, 6,
            &result) == MVMC_KRYLOV_GEVP_OK && result.valid &&
            result.dimension == 3 && result.retained_rank == 3 &&
            result.root_multiplicity == 1 &&
            result.relative_residual <= 64.0 * DBL_EPSILON &&
            result.variance <= 1024.0 * DBL_EPSILON,
        "dimension-three complex unique status=%s residual=%.17g",
        mvmc_krylov_gevp_status_string(result.status),
        result.relative_residual);
}

static void test_degenerate_canonical_projection(void) {
  const double overlap[3] = {2.0, 0.5, 1.0};
  const double forward[3] = {-2.0, -0.5, -1.0};
  const double reverse[3] = {-2.0, -0.5, -1.0};
  const double squared[3] = {2.0, 0.5, 1.0};
  MVMCPowerLanczosGEVPPolicy policy;
  MVMCPowerLanczosGEVPResult result;
  CHECK(mvmc_power_lanczos_gevp_default_policy(0x1p-40, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "degenerate policy");
  CHECK(mvmc_power_lanczos_gevp_solve_real_packed(
            &policy, 2, overlap, forward, reverse, squared, 3,
            &result) == MVMC_KRYLOV_GEVP_OK && result.valid,
        "degenerate solve status=%s",
        mvmc_krylov_gevp_status_string(result.status));
  if (!result.valid) return;
  CHECK(result.root_multiplicity == 2,
        "degenerate multiplicity %d", result.root_multiplicity);
  CHECK(fabs(creal(result.coefficient[0]) - 1.0 / sqrt(2.0)) <=
            2048.0 * DBL_EPSILON &&
            cabs(result.coefficient[1]) <= 2048.0 * DBL_EPSILON,
        "canonical e0 projection %.17g%+.17gi %.17g%+.17gi",
        creal(result.coefficient[0]), cimag(result.coefficient[0]),
        creal(result.coefficient[1]), cimag(result.coefficient[1]));
  CHECK(fabs(result.normalization - 1.0) <= 1024.0 * DBL_EPSILON &&
            fabs(result.energy + 1.0) <= 1024.0 * DBL_EPSILON &&
            result.variance == 0.0,
        "degenerate observables");
}

static void test_fail_closed_boundaries(void) {
  const double zero_s[3] = {0.0, 0.0, 0.0};
  const double zero_kf[3] = {0.0, 0.0, 0.0};
  const double zero_kr[3] = {0.0, 0.0, 0.0};
  const double zero_b[3] = {0.0, 0.0, 0.0};
  const double overlap[3] = {1.0, 0.0, 1.0};
  const double forward[3] = {-2.0, 0.0, 1.0};
  const double reverse[3] = {-2.0, 0.0, 1.0};
  MVMCPowerLanczosGEVPPolicy policy;
  MVMCPowerLanczosGEVPResult result;
  CHECK(mvmc_power_lanczos_gevp_default_policy(0x1p-40, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "boundary policy");
  CHECK(mvmc_power_lanczos_gevp_solve_real_packed(
            &policy, 2, zero_s, zero_kf, zero_kr, zero_b, 3, &result) ==
            MVMC_KRYLOV_GEVP_ZERO_RETAINED_RANK &&
            !result.valid,
        "zero overlap fail-closed status=%s",
        mvmc_krylov_gevp_status_string(result.status));
  CHECK(mvmc_power_lanczos_gevp_solve_real_packed(
            &policy, 2, overlap, forward, reverse, zero_b, 3,
            &result) == MVMC_KRYLOV_GEVP_OK && result.valid &&
            !result.observables_valid && result.variance < 0.0,
        "inconsistent B leaves coefficient valid but observables invalid status=%s",
        mvmc_krylov_gevp_status_string(result.status));
  {
    const double complex invalid_s[3] = {1.0 + 1000000.0 * I, 0.0,
                                         1.0};
    const double complex valid_kf[3] = {-2.0, 0.0, 1.0};
    const double complex valid_kr[3] = {-2.0, 0.0, 1.0};
    const double complex valid_b[3] = {4.0, 0.0, 1.0};
    CHECK(mvmc_power_lanczos_gevp_solve_complex_packed(
              &policy, 2, invalid_s, valid_kf, valid_kr, valid_b, 3,
              &result) == MVMC_KRYLOV_GEVP_ANTIHERMITIAN_INPUT &&
              !result.valid,
          "overlap diagonal imaginary component rejected");
  }
  {
    const double complex valid_s[3] = {1.0, 0.0, 1.0};
    const double complex valid_kf[3] = {-2.0, 0.0, 1.0};
    const double complex valid_kr[3] = {-2.0, 0.0, 1.0};
    const double complex invalid_b[3] = {4.0 + 1000000.0 * I, 0.0,
                                         1.0};
    CHECK(mvmc_power_lanczos_gevp_solve_complex_packed(
              &policy, 2, valid_s, valid_kf, valid_kr, invalid_b, 3,
              &result) == MVMC_KRYLOV_GEVP_ANTIHERMITIAN_INPUT &&
              !result.valid,
          "squared diagonal imaginary component rejected");
  }
  {
    const double complex valid_s[3] = {1.0, 0.0, 1.0};
    const double complex invalid_kf[3] = {1000000.0 * I, 0.0, 1.0};
    const double complex invalid_kr[3] = {-1000000.0 * I, 0.0, 1.0};
    const double complex valid_b[3] = {1.0, 0.0, 1.0};
    CHECK(mvmc_power_lanczos_gevp_solve_complex_packed(
              &policy, 2, valid_s, invalid_kf, invalid_kr, valid_b, 3,
              &result) == MVMC_KRYLOV_GEVP_ANTIHERMITIAN_INPUT &&
              !result.valid,
          "Hamiltonian diagonal imaginary component rejected");
  }
}

typedef union {
  MVMCPowerLanczosGEVPPolicy policy;
  MVMCPowerLanczosGEVPResult result;
  double complex complex_value;
  double real_value;
  unsigned char bytes[sizeof(MVMCPowerLanczosGEVPResult) + 32U];
} AliasStorage;

static void test_alias_no_write(void) {
  const double overlap[3] = {1.0, 0.0, 1.0};
  const double forward[3] = {-2.0, 0.0, 1.0};
  const double reverse[3] = {-2.0, 0.0, 1.0};
  const double squared[3] = {4.0, 0.0, 1.0};
  const double complex complex_forward[3] = {-2.0, 0.0, 1.0};
  const double complex complex_reverse[3] = {-2.0, 0.0, 1.0};
  const double complex complex_squared[3] = {4.0, 0.0, 1.0};
  MVMCPowerLanczosGEVPPolicy policy;
  MVMCPowerLanczosGEVPPolicy policy_snapshot;
  MVMCPowerLanczosGEVPResult result;
  MVMCPowerLanczosGEVPResult result_snapshot;
  CHECK(mvmc_power_lanczos_gevp_default_policy(0x1p-40, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "alias policy");

  {
    double shared[4] = {1.0, 0.0, 1.0, 1.0};
    double snapshot[4];
    memcpy(snapshot, shared, sizeof(snapshot));
    memset(&result, 0x5a, sizeof(result));
    result_snapshot = result;
    policy_snapshot = policy;
    CHECK(mvmc_power_lanczos_gevp_solve_real_packed(
              &policy, 2, shared, shared, reverse, squared, 3, &result) ==
              MVMC_KRYLOV_GEVP_INVALID_ARGUMENT &&
              memcmp(shared, snapshot, sizeof(snapshot)) == 0 &&
              memcmp(&result, &result_snapshot, sizeof(result)) == 0 &&
              memcmp(&policy, &policy_snapshot, sizeof(policy)) == 0,
          "exact input/input alias rejected without write");
    CHECK(mvmc_power_lanczos_gevp_solve_real_packed(
              &policy, 2, shared, shared + 1, reverse, squared, 3,
              &result) == MVMC_KRYLOV_GEVP_INVALID_ARGUMENT &&
              memcmp(shared, snapshot, sizeof(snapshot)) == 0 &&
              memcmp(&result, &result_snapshot, sizeof(result)) == 0 &&
              memcmp(&policy, &policy_snapshot, sizeof(policy)) == 0,
          "partial input/input alias rejected without write");
  }
  {
    AliasStorage storage;
    AliasStorage snapshot;
    double *aliased_input = (double *)(void *)storage.bytes;
    MVMCPowerLanczosGEVPResult *exact_result =
        (MVMCPowerLanczosGEVPResult *)(void *)storage.bytes;
    memset(&storage, 0xa5, sizeof(storage));
    aliased_input[0] = 1.0;
    aliased_input[1] = 0.0;
    aliased_input[2] = 1.0;
    snapshot = storage;
    CHECK(mvmc_power_lanczos_gevp_solve_real_packed(
              &policy, 2, aliased_input, forward, reverse, squared, 3,
              exact_result) == MVMC_KRYLOV_GEVP_INVALID_ARGUMENT &&
              memcmp(&storage, &snapshot, sizeof(storage)) == 0,
          "exact input/result alias rejected without write");
    CHECK(mvmc_power_lanczos_gevp_solve_real_packed(
              &policy, 2, aliased_input, forward, reverse, squared, 3,
              (MVMCPowerLanczosGEVPResult *)(void *)(storage.bytes +
                                                     sizeof(double))) ==
                  MVMC_KRYLOV_GEVP_INVALID_ARGUMENT &&
              memcmp(&storage, &snapshot, sizeof(storage)) == 0,
          "partial input/result alias rejected without write");
  }
  {
    AliasStorage storage;
    AliasStorage snapshot;
    MVMCPowerLanczosGEVPPolicy *aliased_policy =
        (MVMCPowerLanczosGEVPPolicy *)(void *)storage.bytes;
    CHECK(mvmc_power_lanczos_gevp_default_policy(0x1p-40,
                                                  aliased_policy) ==
              MVMC_KRYLOV_GEVP_OK,
          "exact policy/result setup");
    snapshot = storage;
    CHECK(mvmc_power_lanczos_gevp_solve_real_packed(
              aliased_policy, 2, overlap, forward, reverse, squared, 3,
              (MVMCPowerLanczosGEVPResult *)(void *)storage.bytes) ==
                  MVMC_KRYLOV_GEVP_INVALID_ARGUMENT &&
              memcmp(&storage, &snapshot, sizeof(storage)) == 0,
          "exact policy/result alias rejected without write");
    CHECK(mvmc_power_lanczos_gevp_solve_real_packed(
              aliased_policy, 2, overlap, forward, reverse, squared, 3,
              (MVMCPowerLanczosGEVPResult *)(void *)(storage.bytes +
                                                     sizeof(double))) ==
                  MVMC_KRYLOV_GEVP_INVALID_ARGUMENT &&
              memcmp(&storage, &snapshot, sizeof(storage)) == 0,
          "partial policy/result alias rejected without write");
  }
  {
    double complex shared[4] = {1.0, 0.0, 1.0, 1.0};
    double complex snapshot[4];
    memcpy(snapshot, shared, sizeof(snapshot));
    memset(&result, 0x5a, sizeof(result));
    result_snapshot = result;
    CHECK(mvmc_power_lanczos_gevp_solve_complex_packed(
              &policy, 2, shared, shared + 1, complex_reverse,
              complex_squared, 3, &result) ==
              MVMC_KRYLOV_GEVP_INVALID_ARGUMENT &&
              memcmp(shared, snapshot, sizeof(snapshot)) == 0 &&
              memcmp(&result, &result_snapshot, sizeof(result)) == 0,
          "complex partial input/input alias rejected without write");
  }
  {
    AliasStorage storage;
    AliasStorage snapshot;
    double complex *aliased_input =
        (double complex *)(void *)storage.bytes;
    memset(&storage, 0xa5, sizeof(storage));
    aliased_input[0] = 1.0;
    aliased_input[1] = 0.0;
    aliased_input[2] = 1.0;
    snapshot = storage;
    CHECK(mvmc_power_lanczos_gevp_solve_complex_packed(
              &policy, 2, aliased_input, complex_forward,
              complex_reverse, complex_squared, 3,
              (MVMCPowerLanczosGEVPResult *)(void *)storage.bytes) ==
                  MVMC_KRYLOV_GEVP_INVALID_ARGUMENT &&
              memcmp(&storage, &snapshot, sizeof(storage)) == 0,
          "complex exact input/result alias rejected without write");
    CHECK(mvmc_power_lanczos_gevp_solve_complex_packed(
              &policy, 2, aliased_input, complex_forward,
              complex_reverse, complex_squared, SIZE_MAX,
              (MVMCPowerLanczosGEVPResult *)(void *)storage.bytes) ==
                  MVMC_KRYLOV_GEVP_INVALID_ARGUMENT &&
              memcmp(&storage, &snapshot, sizeof(storage)) == 0,
          "complex overflow/alias rejected without write");
  }
}

int main(void) {
  test_policy();
  test_real_unique_and_scale();
  test_complex_unique();
  test_dimension_three_unique();
  test_degenerate_canonical_projection();
  test_fail_closed_boundaries();
  test_alias_no_write();
  if (failures != 0) {
    fprintf(stderr, "power-Lanczos GEVP unit failures=%d\n", failures);
    return 1;
  }
  puts("power-Lanczos GEVP unit: PASS");
  return 0;
}
