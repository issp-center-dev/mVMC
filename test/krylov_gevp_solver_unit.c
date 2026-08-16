#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "krylov_gevp_solver.h"

static int failures = 0;

#define CHECK(condition, ...)                                              \
  do {                                                                     \
    if (!(condition)) {                                                    \
      fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__);               \
      fprintf(stderr, __VA_ARGS__);                                        \
      fputc('\n', stderr);                                                 \
      ++failures;                                                          \
    }                                                                      \
  } while (0)

enum { PHYSICAL_DIMENSION = 3, MAX_UPPER = 6 };

static size_t upper_index(int dimension, int row, int column) {
  size_t index = 0;
  int current_row;
  for (current_row = 0; current_row < row; ++current_row) {
    index += (size_t)(dimension - current_row);
  }
  return index + (size_t)(column - row);
}

static double complex inner_product(const double complex *left,
                                    const double complex *right) {
  double complex value = 0.0;
  int index;
  for (index = 0; index < PHYSICAL_DIMENSION; ++index) {
    value += conj(left[index]) * right[index];
  }
  return value;
}

static void matrix_vector(const double complex *matrix,
                          const double complex *vector,
                          double complex *result) {
  int row;
  int column;
  for (row = 0; row < PHYSICAL_DIMENSION; ++row) {
    result[row] = 0.0;
    for (column = 0; column < PHYSICAL_DIMENSION; ++column) {
      result[row] += matrix[row * PHYSICAL_DIMENSION + column] *
                     vector[column];
    }
  }
}

static void build_packed_matrices(
    const double complex *hamilian, const double complex *psi,
    int dimension, double common_scale, double complex *overlap,
    double complex *forward, double complex *reverse,
    double complex *squared) {
  double complex vector[4][PHYSICAL_DIMENSION] = {{0.0}};
  int order;
  int row;
  int column;
  memcpy(vector[0], psi, sizeof(vector[0]));
  for (order = 1; order < 4; ++order) {
    matrix_vector(hamilian, vector[order - 1], vector[order]);
  }
  for (row = 0; row < dimension; ++row) {
    for (column = row; column < dimension; ++column) {
      const size_t entry = upper_index(dimension, row, column);
      overlap[entry] = common_scale *
                       inner_product(vector[row], vector[column]);
      forward[entry] = common_scale *
                       inner_product(vector[row], vector[column + 1]);
      reverse[entry] = common_scale *
                       inner_product(vector[column], vector[row + 1]);
      squared[entry] = common_scale *
                       inner_product(vector[row + 1], vector[column + 1]);
    }
  }
}

static double projective_distance(const double complex *left,
                                  const double complex *right,
                                  int dimension) {
  double complex overlap = 0.0;
  double complex phase = 1.0;
  double norm = 0.0;
  int index;
  for (index = 0; index < dimension; ++index) {
    overlap += conj(left[index]) * right[index];
  }
  if (cabs(overlap) > 0.0) phase = overlap / cabs(overlap);
  for (index = 0; index < dimension; ++index) {
    norm = hypot(norm, cabs(left[index] * phase - right[index]));
  }
  return norm;
}

static double rank_one_projector_error(const MVMCKrylovGEVPResult *result) {
  double error = 0.0;
  int row;
  int column;
  for (row = 0; row < result->dimension; ++row) {
    for (column = 0; column < result->dimension; ++column) {
      const double complex expected =
          result->coefficient[row] * conj(result->coefficient[column]);
      error = hypot(
          error,
          cabs(result->root_subspace_projector[
                   row * result->dimension + column] - expected));
    }
  }
  return error;
}

static void test_policy(void) {
  MVMCKrylovGEVPPolicy policy;
  memset(&policy, 0xa5, sizeof(policy));
  CHECK(mvmc_krylov_gevp_default_policy(1.0e-10, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "default policy");
  CHECK(policy.valid &&
            policy.policy_version == MVMC_KRYLOV_GEVP_POLICY_VERSION &&
            policy.rank_relative_cutoff == 1.0e-10 &&
            policy.overlap_negative_relative_tolerance == 1.0e-10 &&
            policy.maximum_input_antihermitian_effect == 0.25 &&
            policy.maximum_gevp_relative_residual == 1.0e-10 &&
            policy.degenerate_root_gap_relative_threshold == 1.0e-10,
        "default policy fields");
  CHECK(mvmc_krylov_gevp_default_policy(0.0, &policy) ==
            MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "zero cutoff rejected");
  CHECK(mvmc_krylov_gevp_default_policy(1.0, &policy) ==
            MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "unit cutoff rejected");
  CHECK(mvmc_krylov_gevp_default_policy(NAN, &policy) ==
            MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "nonfinite cutoff rejected");
  CHECK(mvmc_krylov_gevp_default_policy(1.0e-10, NULL) ==
            MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "null policy rejected");
}

static void test_real_exact_and_common_scale(void) {
  const double complex hamiltonian[9] = {
      -2.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 3.0};
  const double complex psi[3] = {1.0, 0.8, -0.5};
  const double common_scale[3] = {
      1.0, 3.8725919148493183e-121, 2.5822498780869086e120};
  const double expected_energy[2] = {-1.8127994975654498, -2.0};
  MVMCKrylovGEVPPolicy policy;
  MVMCKrylovGEVPResult reference[2];
  int dimension;
  int scale_index;
  CHECK(mvmc_krylov_gevp_default_policy(1.0e-10, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "real policy");
  memset(reference, 0, sizeof(reference));
  for (dimension = 2; dimension <= 3; ++dimension) {
    const size_t count = (size_t)dimension * (size_t)(dimension + 1) / 2U;
    for (scale_index = 0; scale_index < 3; ++scale_index) {
      double complex overlap_complex[MAX_UPPER] = {0.0};
      double complex forward_complex[MAX_UPPER] = {0.0};
      double complex reverse_complex[MAX_UPPER] = {0.0};
      double complex squared_complex[MAX_UPPER] = {0.0};
      double overlap[MAX_UPPER] = {0.0};
      double forward[MAX_UPPER] = {0.0};
      double reverse[MAX_UPPER] = {0.0};
      double squared[MAX_UPPER] = {0.0};
      MVMCKrylovGEVPResult result;
      size_t entry;
      build_packed_matrices(hamiltonian, psi, dimension,
                            common_scale[scale_index], overlap_complex,
                            forward_complex, reverse_complex,
                            squared_complex);
      for (entry = 0; entry < count; ++entry) {
        overlap[entry] = creal(overlap_complex[entry]);
        forward[entry] = creal(forward_complex[entry]);
        reverse[entry] = creal(reverse_complex[entry]);
        squared[entry] = creal(squared_complex[entry]);
      }
      CHECK(mvmc_krylov_gevp_solve_real_packed(
                &policy, dimension, overlap, forward, reverse, squared,
                count, &result) == MVMC_KRYLOV_GEVP_OK,
            "real dimension %d scale %d solve: %s", dimension,
            scale_index, mvmc_krylov_gevp_status_string(result.status));
      if (!result.valid) continue;
      CHECK(fabs(result.energy - expected_energy[dimension - 2]) <= 2.0e-11,
            "real dimension %d energy %.17g", dimension, result.energy);
      CHECK(result.retained_rank == dimension && result.discarded_rank == 0,
            "real dimension %d rank %d", dimension, result.retained_rank);
      CHECK(result.gevp_relative_residual <= 1.0e-12,
            "real dimension %d residual %.17g", dimension,
            result.gevp_relative_residual);
      CHECK(result.hamiltonian_antihermitian_residual <= 1.0e-14,
            "real dimension %d Hermitian residual %.17g", dimension,
            result.hamiltonian_antihermitian_residual);
      CHECK(result.phase_pivot >= 0 &&
                creal(result.coefficient[result.phase_pivot]) > 0.0 &&
                cimag(result.coefficient[result.phase_pivot]) == 0.0,
            "real phase convention");
      CHECK(result.root_multiplicity == 1 &&
                !result.coefficient_comparison_uses_projector &&
                rank_one_projector_error(&result) <= 2.0e-12,
            "real nondegenerate root projector");
      for (int projector_index = 0;
           projector_index < dimension * dimension; ++projector_index) {
        CHECK(cimag(result.root_subspace_projector[projector_index]) == 0.0,
              "real projector entry %d has zero imaginary part",
              projector_index);
      }
      if (scale_index == 0) {
        reference[dimension - 2] = result;
      } else if (reference[dimension - 2].valid) {
        CHECK(fabs(result.energy - reference[dimension - 2].energy) <=
                  3.0e-12,
              "common scale energy invariance dimension %d", dimension);
        CHECK(projective_distance(result.coefficient,
                                  reference[dimension - 2].coefficient,
                                  dimension) <= 2.0e-10,
              "common scale coefficient invariance dimension %d",
              dimension);
      }
      if (dimension == 3) {
        CHECK(result.variance <= 5.0e-11,
              "full real Krylov ground variance %.17g", result.variance);
      } else {
        CHECK(result.variance > 0.0,
              "restricted real Krylov variance positive");
      }
    }
  }
}

static void test_degenerate_root_projector(void) {
  const double overlap[3] = {1.0, 0.0, 1.0};
  const double forward[3] = {-1.0, 0.0, -1.0};
  const double reverse[3] = {-1.0, 0.0, -1.0};
  const double squared[3] = {1.0, 0.0, 1.0};
  MVMCKrylovGEVPPolicy policy;
  MVMCKrylovGEVPResult result;
  int row;
  int column;
  CHECK(mvmc_krylov_gevp_default_policy(1.0e-10, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "degenerate policy");
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, overlap, forward, reverse, squared, 3,
            &result) == MVMC_KRYLOV_GEVP_OK && result.valid,
        "degenerate root solve: %s",
        mvmc_krylov_gevp_status_string(result.status));
  CHECK(result.root_multiplicity == 2 &&
            result.coefficient_comparison_uses_projector &&
            result.root_gap == 0.0,
        "degenerate root classification multiplicity=%d gap=%.17g",
        result.root_multiplicity, result.root_gap);
  for (row = 0; row < 2; ++row) {
    for (column = 0; column < 2; ++column) {
      const double complex expected = row == column ? 1.0 : 0.0;
      CHECK(cabs(result.root_subspace_projector[row * 2 + column] -
                 expected) <= 2.0e-14,
            "degenerate projector (%d,%d)", row, column);
      CHECK(cimag(result.root_subspace_projector[row * 2 + column]) == 0.0,
            "degenerate real projector (%d,%d) has zero imaginary part",
            row, column);
    }
  }
}

static void test_complex_exact(void) {
  const double complex hamiltonian[9] = {
      -1.0, 0.4 * I, 0.2 + 0.1 * I,
      -0.4 * I, 0.5, -0.3 * I,
      0.2 - 0.1 * I, 0.3 * I, 2.0};
  const double complex psi[3] = {1.0, 0.7 + 0.2 * I, -0.4 + 0.6 * I};
  const double expected_energy[2] = {
      -1.0106940604792756, -1.1079332179665404};
  MVMCKrylovGEVPPolicy policy;
  int dimension;
  CHECK(mvmc_krylov_gevp_default_policy(1.0e-10, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "complex policy");
  for (dimension = 2; dimension <= 3; ++dimension) {
    double complex overlap[MAX_UPPER] = {0.0};
    double complex forward[MAX_UPPER] = {0.0};
    double complex reverse[MAX_UPPER] = {0.0};
    double complex squared[MAX_UPPER] = {0.0};
    double complex overlap_copy[MAX_UPPER];
    double complex forward_copy[MAX_UPPER];
    double complex reverse_copy[MAX_UPPER];
    double complex squared_copy[MAX_UPPER];
    MVMCKrylovGEVPResult result;
    const size_t count = (size_t)dimension * (size_t)(dimension + 1) / 2U;
    build_packed_matrices(hamiltonian, psi, dimension, 1.0, overlap,
                          forward, reverse, squared);
    memcpy(overlap_copy, overlap, sizeof(overlap));
    memcpy(forward_copy, forward, sizeof(forward));
    memcpy(reverse_copy, reverse, sizeof(reverse));
    memcpy(squared_copy, squared, sizeof(squared));
    CHECK(mvmc_krylov_gevp_solve_complex_packed(
              &policy, dimension, overlap, forward, reverse, squared,
              count, &result) == MVMC_KRYLOV_GEVP_OK,
          "complex dimension %d solve: %s", dimension,
          mvmc_krylov_gevp_status_string(result.status));
    CHECK(result.valid &&
              fabs(result.energy - expected_energy[dimension - 2]) <=
                  3.0e-11,
          "complex dimension %d energy %.17g", dimension, result.energy);
    CHECK(result.gevp_relative_residual <= 1.0e-12,
          "complex dimension %d residual %.17g", dimension,
          result.gevp_relative_residual);
    CHECK(memcmp(overlap, overlap_copy, sizeof(overlap)) == 0 &&
              memcmp(forward, forward_copy, sizeof(forward)) == 0 &&
              memcmp(reverse, reverse_copy, sizeof(reverse)) == 0 &&
              memcmp(squared, squared_copy, sizeof(squared)) == 0,
          "complex inputs remain byte-identical");
    if (dimension == 3) {
      CHECK(result.variance <= 8.0e-11,
            "full complex Krylov ground variance %.17g", result.variance);
    }
  }
}

static void test_relative_rank_cutoff(void) {
  const double overlap[3] = {1.0, 1.0 - 1.0e-12, 1.0};
  const double forward[3] = {-1.0, -(1.0 - 1.0e-12), -1.0};
  const double reverse[3] = {-1.0, -(1.0 - 1.0e-12), -1.0};
  const double squared[3] = {1.0, 1.0 - 1.0e-12, 1.0};
  MVMCKrylovGEVPPolicy policy;
  MVMCKrylovGEVPResult result;
  CHECK(mvmc_krylov_gevp_default_policy(1.0e-10, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "rank policy");
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, overlap, forward, reverse, squared, 3, &result) ==
            MVMC_KRYLOV_GEVP_OK,
        "rank-one solve: %s", mvmc_krylov_gevp_status_string(result.status));
  CHECK(result.valid && result.retained_rank == 1 &&
            result.discarded_rank == 1 && fabs(result.energy + 1.0) < 1e-13,
        "rank-one result rank=%d energy=%.17g", result.retained_rank,
        result.energy);
  CHECK(result.gevp_relative_residual <= 1.0e-14,
        "rank-one full residual %.17g", result.gevp_relative_residual);
}

static void test_registered_cutoff_scan(void) {
  const double complex hamiltonian[9] = {
      -2.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 3.0};
  const double complex psi[3] = {1.0, 0.8, -0.5};
  const double cutoff[3] = {1.0e-12, 1.0e-10, 1.0e-8};
  double complex overlap_complex[MAX_UPPER] = {0.0};
  double complex forward_complex[MAX_UPPER] = {0.0};
  double complex reverse_complex[MAX_UPPER] = {0.0};
  double complex squared_complex[MAX_UPPER] = {0.0};
  double overlap[MAX_UPPER];
  double forward[MAX_UPPER];
  double reverse[MAX_UPPER];
  double squared[MAX_UPPER];
  MVMCKrylovGEVPPolicy policy;
  MVMCKrylovGEVPResult result;
  int scan;
  size_t entry;
  build_packed_matrices(hamiltonian, psi, 3, 1.0, overlap_complex,
                        forward_complex, reverse_complex, squared_complex);
  for (entry = 0; entry < MAX_UPPER; ++entry) {
    overlap[entry] = creal(overlap_complex[entry]);
    forward[entry] = creal(forward_complex[entry]);
    reverse[entry] = creal(reverse_complex[entry]);
    squared[entry] = creal(squared_complex[entry]);
  }
  for (scan = 0; scan < 3; ++scan) {
    CHECK(mvmc_krylov_gevp_default_policy(cutoff[scan], &policy) ==
              MVMC_KRYLOV_GEVP_OK,
          "registered cutoff policy %d", scan);
    CHECK(mvmc_krylov_gevp_solve_real_packed(
              &policy, 3, overlap, forward, reverse, squared, MAX_UPPER,
              &result) == MVMC_KRYLOV_GEVP_OK,
          "registered cutoff solve %d: %s", scan,
          mvmc_krylov_gevp_status_string(result.status));
    CHECK(result.retained_rank == 3 && fabs(result.energy + 2.0) < 3e-11,
          "registered cutoff %g rank=%d energy=%.17g", cutoff[scan],
          result.retained_rank, result.energy);
  }

  {
    const double epsilon = 1.0e-9;
    const double sensitive_overlap[3] = {
        1.0 + 0.5e-9, 1.0 - 0.5e-9, 1.0 + 0.5e-9};
    const double sensitive_forward[3] = {
        -1.0 - epsilon, -1.0 + epsilon, -1.0 - epsilon};
    const double sensitive_reverse[3] = {
        -1.0 - epsilon, -1.0 + epsilon, -1.0 - epsilon};
    const double sensitive_squared[3] = {
        1.0 + 2.0 * epsilon, 1.0 - 2.0 * epsilon,
        1.0 + 2.0 * epsilon};
    double energy[3];
    int rank[3];
    for (scan = 0; scan < 3; ++scan) {
      CHECK(mvmc_krylov_gevp_default_policy(cutoff[scan], &policy) ==
                MVMC_KRYLOV_GEVP_OK,
            "sensitivity cutoff policy");
      CHECK(mvmc_krylov_gevp_solve_real_packed(
                &policy, 2, sensitive_overlap, sensitive_forward,
                sensitive_reverse, sensitive_squared, 3, &result) ==
                MVMC_KRYLOV_GEVP_OK,
            "sensitivity cutoff solve %d: %s", scan,
            mvmc_krylov_gevp_status_string(result.status));
      energy[scan] = result.energy;
      rank[scan] = result.retained_rank;
    }
    CHECK(rank[0] == 2 && rank[1] == 2 && rank[2] == 1,
          "sensitivity rank pattern %d/%d/%d", rank[0], rank[1], rank[2]);
    CHECK(fabs(energy[0] + 2.0) < 1e-6 &&
              fabs(energy[1] + 2.0) < 1e-6 &&
              fabs(energy[2] + 1.0) < 2e-7,
          "sensitivity expected STOP energies %.17g %.17g %.17g",
          energy[0], energy[1], energy[2]);
  }
}

static void test_variance_and_discarded_residual(void) {
  const double overlap[3] = {1.0, 0.0, 1.0};
  const double forward[3] = {-1.0, 0.0, 1.0};
  const double reverse[3] = {-1.0, 0.0, 1.0};
  double squared[3] = {1.0 - 5.0e-11, 0.0, 1.0};
  MVMCKrylovGEVPPolicy policy;
  MVMCKrylovGEVPResult result;
  CHECK(mvmc_krylov_gevp_default_policy(1.0e-10, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "variance policy");
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, overlap, forward, reverse, squared, 3, &result) ==
            MVMC_KRYLOV_GEVP_OK,
        "roundoff-negative variance solve");
  CHECK(result.variance_clamped && result.variance == 0.0,
        "roundoff-negative variance clamped");
  squared[0] = 1.0 - 2.0e-10;
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, overlap, forward, reverse, squared, 3, &result) ==
            MVMC_KRYLOV_GEVP_NEGATIVE_VARIANCE,
        "material negative variance rejected: %s",
        mvmc_krylov_gevp_status_string(result.status));

  {
    const double epsilon = 1.0e-12;
    const double nearly_dependent_overlap[3] = {
        1.0 + 0.5 * epsilon, 1.0 - 0.5 * epsilon,
        1.0 + 0.5 * epsilon};
    const double coupled_hamiltonian[3] = {-0.9, -1.0, -1.1};
    const double coupled_reverse[3] = {-0.9, -1.0, -1.1};
    const double positive_squared[3] = {1.0, 0.0, 1.0};
    CHECK(mvmc_krylov_gevp_solve_real_packed(
              &policy, 2, nearly_dependent_overlap, coupled_hamiltonian,
              coupled_reverse, positive_squared, 3, &result) ==
              MVMC_KRYLOV_GEVP_RESIDUAL_FAILURE,
          "discarded coupling caught by full residual: %s residual=%.17g rank=%d eigen=%.17g/%.17g cutoff=%.17g",
          mvmc_krylov_gevp_status_string(result.status),
          result.gevp_relative_residual, result.retained_rank,
          result.overlap_eigenvalue[0], result.overlap_eigenvalue[1],
          result.rank_relative_cutoff);
  }
}

static void test_residual_denominator_overflow(void) {
  const double scale = 0.375 * DBL_MAX;
  const double coupling = 1.0e-8;
  const double overlap[3] = {scale, scale, scale};
  const double forward[3] = {
      scale * (1.0 + coupling), scale, scale * (1.0 - coupling)};
  const double reverse[3] = {
      scale * (1.0 + coupling), scale, scale * (1.0 - coupling)};
  const double squared[3] = {scale, scale, scale};
  MVMCKrylovGEVPPolicy policy;
  MVMCKrylovGEVPResult result;

  CHECK(mvmc_krylov_gevp_default_policy(1.0e-10, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "residual denominator overflow policy");
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, overlap, forward, reverse, squared, 3, &result) ==
            MVMC_KRYLOV_GEVP_RESIDUAL_FAILURE,
        "finite input with overflowing residual denominator must not pass: "
        "status=%s residual=%.17g energy=%.17g",
        mvmc_krylov_gevp_status_string(result.status),
        result.gevp_relative_residual, result.energy);
  CHECK(isfinite(result.gevp_relative_residual) &&
            result.gevp_relative_residual >
                policy.maximum_gevp_relative_residual,
        "overflow-safe residual reports the failing finite ratio");

  {
    const double larger_energy_scale = 0.2 * DBL_MAX;
    const double larger_energy = 1.5;
    const double larger_overlap[3] = {
        larger_energy_scale, larger_energy_scale, larger_energy_scale};
    const double larger_forward[3] = {
        larger_energy_scale * (larger_energy + coupling),
        larger_energy_scale * larger_energy,
        larger_energy_scale * (larger_energy - coupling)};
    const double larger_reverse[3] = {
        larger_energy_scale * (larger_energy + coupling),
        larger_energy_scale * larger_energy,
        larger_energy_scale * (larger_energy - coupling)};
    const double larger_squared[3] = {
        larger_energy_scale * larger_energy * larger_energy,
        larger_energy_scale * larger_energy * larger_energy,
        larger_energy_scale * larger_energy * larger_energy};
    CHECK(mvmc_krylov_gevp_solve_real_packed(
              &policy, 2, larger_overlap, larger_forward, larger_reverse,
              larger_squared, 3, &result) ==
              MVMC_KRYLOV_GEVP_RESIDUAL_FAILURE,
          "overflow-safe residual scaling for |energy| > 1: "
          "status=%s residual=%.17g energy=%.17g",
          mvmc_krylov_gevp_status_string(result.status),
          result.gevp_relative_residual, result.energy);
    CHECK(isfinite(result.gevp_relative_residual) &&
              result.gevp_relative_residual >
                  policy.maximum_gevp_relative_residual,
          "overflow-safe |energy| > 1 residual reports the failing ratio");
  }
}

static void test_extreme_finite_hermitian_average(void) {
  const double overlap[3] = {
      0.75 * DBL_MAX, 0.0, 0.25 * DBL_MAX};
  const double forward[3] = {
      0.75 * DBL_MAX, 0.0, 0.25 * DBL_MAX};
  const double reverse[3] = {
      0.75 * DBL_MAX, 0.0, 0.25 * DBL_MAX};
  const double squared[3] = {
      0.75 * DBL_MAX, 0.0, 0.25 * DBL_MAX};
  MVMCKrylovGEVPPolicy policy;
  MVMCKrylovGEVPResult result;

  CHECK(mvmc_krylov_gevp_default_policy(1.0e-10, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "extreme Hermitian average policy");
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, overlap, forward, reverse, squared, 3, &result) ==
            MVMC_KRYLOV_GEVP_OK,
        "finite Hermitian average must not overflow its intermediate sum: %s",
        mvmc_krylov_gevp_status_string(result.status));
  CHECK(result.valid && fabs(result.energy - 1.0) <= 64.0 * DBL_EPSILON &&
            result.variance == 0.0 && result.gevp_relative_residual == 0.0,
        "extreme Hermitian average result");
}

static void test_zero_residual_denominator(void) {
  const double overlap[3] = {1.0, 0.0, 1.0};
  const double forward[3] = {0.0, 0.0, 0.0};
  const double reverse[3] = {0.0, 0.0, 0.0};
  const double squared[3] = {0.0, 0.0, 0.0};
  MVMCKrylovGEVPPolicy policy;
  MVMCKrylovGEVPResult result;

  CHECK(mvmc_krylov_gevp_default_policy(1.0e-10, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "zero residual denominator policy");
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, overlap, forward, reverse, squared, 3, &result) ==
            MVMC_KRYLOV_GEVP_OK,
        "zero denominator and zero residual is valid: %s",
        mvmc_krylov_gevp_status_string(result.status));
  CHECK(result.valid && result.energy == 0.0 && result.energy_squared == 0.0 &&
            result.variance == 0.0 && result.gevp_relative_residual == 0.0,
        "zero residual denominator result");
}

static void test_failure_paths(void) {
  const double valid_overlap[3] = {1.0, 0.0, 1.0};
  const double valid_forward[3] = {-1.0, 0.0, 1.0};
  const double valid_reverse[3] = {-1.0, 0.0, 1.0};
  const double valid_squared[3] = {1.0, 0.0, 1.0};
  MVMCKrylovGEVPPolicy policy;
  MVMCKrylovGEVPResult result;
  double overlap[3];
  double forward[3];
  double reverse[3];
  double squared[3];
  double complex complex_overlap[3] = {1.0, 0.0, 1.0};
  double complex complex_forward[3] = {-1.0, 0.0, 1.0};
  double complex complex_reverse[3] = {-1.0, 0.0, 1.0};
  double complex complex_squared[3] = {1.0, 0.0, 1.0};
  CHECK(mvmc_krylov_gevp_default_policy(1.0e-10, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "failure policy");

  memcpy(overlap, valid_overlap, sizeof(overlap));
  overlap[1] = 2.0;
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, overlap, valid_forward, valid_reverse,
            valid_squared, 3, &result) ==
            MVMC_KRYLOV_GEVP_INDEFINITE_OVERLAP,
        "indefinite overlap rejected: %s",
        mvmc_krylov_gevp_status_string(result.status));

  memcpy(squared, valid_squared, sizeof(squared));
  squared[1] = 2.0;
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, valid_overlap, valid_forward, valid_reverse,
            squared, 3, &result) ==
            MVMC_KRYLOV_GEVP_INDEFINITE_HAMILTONIAN_SQUARED,
        "indefinite B rejected: %s",
        mvmc_krylov_gevp_status_string(result.status));

  memcpy(forward, valid_forward, sizeof(forward));
  memcpy(reverse, valid_reverse, sizeof(reverse));
  forward[1] = 1.0;
  reverse[1] = -1.0;
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, valid_overlap, forward, reverse, valid_squared, 3,
            &result) == MVMC_KRYLOV_GEVP_ANTIHERMITIAN_INPUT,
        "anti-Hermitian effect rejected: %s residual %.17g",
        mvmc_krylov_gevp_status_string(result.status),
        result.hamiltonian_antihermitian_residual);

  memcpy(overlap, valid_overlap, sizeof(overlap));
  overlap[0] = NAN;
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, overlap, valid_forward, valid_reverse,
            valid_squared, 3, &result) ==
            MVMC_KRYLOV_GEVP_NONFINITE_INPUT,
        "NaN rejected");

  memcpy(overlap, valid_overlap, sizeof(overlap));
  overlap[2] = 0.0;
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, overlap, valid_forward, valid_reverse,
            valid_squared, 3, &result) ==
            MVMC_KRYLOV_GEVP_NONPOSITIVE_OVERLAP_DIAGONAL,
        "zero diagonal rejected");

  complex_overlap[0] = I;
  CHECK(mvmc_krylov_gevp_solve_complex_packed(
            &policy, 2, complex_overlap, complex_forward, complex_reverse,
            complex_squared, 3, &result) ==
            MVMC_KRYLOV_GEVP_ANTIHERMITIAN_INPUT,
        "pure-imaginary overlap diagonal rejected by relative effect guard");

  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 1, valid_overlap, valid_forward, valid_reverse,
            valid_squared, 1, &result) ==
            MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "dimension one rejected");
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, valid_overlap, valid_forward, valid_reverse,
            valid_squared, 2, &result) ==
            MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "wrong upper count rejected");

  policy.valid = 0;
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, valid_overlap, valid_forward, valid_reverse,
            valid_squared, 3, &result) ==
            MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "invalid policy rejected");
  CHECK(mvmc_krylov_gevp_default_policy(1.0e-10, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "restore policy for degeneracy-threshold negative");
  policy.degenerate_root_gap_relative_threshold = NAN;
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, valid_overlap, valid_forward, valid_reverse,
            valid_squared, 3, &result) ==
            MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "nonfinite degeneracy threshold rejected");
}

static void test_alias_rejection(void) {
  const double valid_overlap[3] = {1.0, 0.0, 1.0};
  const double valid_reverse[3] = {-1.0, 0.0, 1.0};
  const double valid_squared[3] = {1.0, 0.0, 1.0};
  MVMCKrylovGEVPPolicy policy;
  MVMCKrylovGEVPResult result;
  double shared[3];
  double snapshot[3];
  union {
    MVMCKrylovGEVPResult result;
    double complex packed[sizeof(MVMCKrylovGEVPResult) /
                              sizeof(double complex) +
                          1];
  } output_alias;
  union {
    MVMCKrylovGEVPResult result;
    double complex packed[sizeof(MVMCKrylovGEVPResult) /
                              sizeof(double complex) +
                          1];
  } output_snapshot;
  double complex forward[3] = {-1.0, 0.0, 1.0};
  double complex reverse[3] = {-1.0, 0.0, 1.0};
  double complex squared[3] = {1.0, 0.0, 1.0};
  CHECK(mvmc_krylov_gevp_default_policy(1.0e-10, &policy) ==
            MVMC_KRYLOV_GEVP_OK,
        "alias policy");
  memcpy(shared, valid_overlap, sizeof(shared));
  memcpy(snapshot, shared, sizeof(snapshot));
  CHECK(mvmc_krylov_gevp_solve_real_packed(
            &policy, 2, shared, shared, valid_reverse, valid_squared, 3,
            &result) == MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "input alias rejected");
  CHECK(memcmp(shared, snapshot, sizeof(shared)) == 0,
        "input alias rejection is transactional");

  memset(&output_alias, 0x5a, sizeof(output_alias));
  memcpy(output_snapshot.packed, output_alias.packed, sizeof(output_alias));
  output_alias.packed[0] = 1.0;
  output_alias.packed[1] = 0.0;
  output_alias.packed[2] = 1.0;
  memcpy(&output_snapshot, &output_alias, sizeof(output_alias));
  CHECK(mvmc_krylov_gevp_solve_complex_packed(
            &policy, 2, output_alias.packed, forward, reverse, squared, 3,
            &output_alias.result) == MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
        "output/input byte-range alias rejected");
  CHECK(memcmp(&output_alias, &output_snapshot, sizeof(output_alias)) == 0,
        "output alias rejection performs no write");
}

int main(void) {
  test_policy();
  test_real_exact_and_common_scale();
  test_complex_exact();
  test_degenerate_root_projector();
  test_relative_rank_cutoff();
  test_registered_cutoff_scan();
  test_variance_and_discarded_residual();
  test_residual_denominator_overflow();
  test_extreme_finite_hermitian_average();
  test_zero_residual_denominator();
  test_failure_paths();
  test_alias_rejection();
  if (failures != 0) {
    fprintf(stderr, "krylov GEVP solver unit failures=%d\n", failures);
    return 1;
  }
  puts("krylov GEVP solver unit: PASS");
  return 0;
}
