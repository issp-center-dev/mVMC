#include "krylov_final_state_sensitivity.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failures = 0;
static int world_rank = 0;

#define CHECK(condition, ...)                                                 \
  do {                                                                        \
    if (!(condition)) {                                                       \
      fprintf(stderr,                                                         \
              "KrylovFinalStateSensitivity_Unit FAIL rank %d: ",            \
              world_rank);                                                    \
      fprintf(stderr, __VA_ARGS__);                                           \
      fprintf(stderr, " (line %d)\n", __LINE__);                            \
      ++failures;                                                             \
    }                                                                         \
  } while (0)

static int close_double(double actual, double expected, double tolerance) {
  return fabs(actual - expected) <= tolerance * fmax(1.0, fabs(expected));
}

static MVMCKrylovFinalStateQuadraticModel model2(
    const double complex *overlap, const double complex *hamiltonian,
    const double complex *squared, const double complex *observable) {
  MVMCKrylovFinalStateQuadraticModel model;
  const MVMCKrylovStatus status =
      mvmc_krylov_final_state_quadratic_model_create(
          2, overlap, hamiltonian, squared, observable, 3, &model);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && model.valid,
        "dimension-2 model creation status=%d", status);
  return model;
}

static void test_exact_metrics_and_invariance(void) {
  const double complex overlap[] = {1.0, 0.0, 1.0};
  const double complex hamiltonian[] = {0.0, 1.0, 0.0};
  const double complex squared[] = {1.0, 0.0, 1.0};
  const double complex observable[] = {1.0, 0.0, -1.0};
  const double scale = 1.0 / sqrt(2.0);
  const double complex coefficient[] = {scale, -scale};
  const double complex phase = 2.5 * (cos(0.4) + I * sin(0.4));
  const double complex rotated[] = {
      phase * coefficient[0], phase * coefficient[1]};
  const MVMCKrylovFinalStateQuadraticModel model =
      model2(overlap, hamiltonian, squared, observable);
  MVMCKrylovFinalStateQuadraticMetrics metrics;
  MVMCKrylovFinalStateQuadraticMetrics rotated_metrics;

  CHECK(mvmc_krylov_final_state_quadratic_metrics(
            &model, coefficient, 2, &metrics) == MVMC_KRYLOV_STATUS_OK &&
            metrics.valid && close_double(metrics.norm, 1.0, 2.0e-14) &&
            close_double(metrics.energy, -1.0, 2.0e-14) &&
            close_double(metrics.hamiltonian_second_moment, 1.0,
                         2.0e-14) &&
            close_double(metrics.full_support_variance, 0.0, 2.0e-14) &&
            close_double(metrics.observable, 0.0, 2.0e-14),
        "exact dimension-2 metrics");
  CHECK(mvmc_krylov_final_state_quadratic_metrics(
            &model, rotated, 2, &rotated_metrics) ==
            MVMC_KRYLOV_STATUS_OK &&
            close_double(rotated_metrics.energy, metrics.energy,
                         2.0e-14) &&
            close_double(rotated_metrics.full_support_variance,
                         metrics.full_support_variance, 2.0e-14) &&
            close_double(rotated_metrics.observable, metrics.observable,
                         2.0e-14),
        "global coefficient scale/phase invariance");
}

static void test_sensitivity_decomposition(void) {
  const double complex overlap[] = {1.0, 0.0, 1.0};
  const double complex hamiltonian[] = {0.0, 1.0, 0.0};
  const double complex squared[] = {1.0, 0.0, 1.0};
  const double complex observable[] = {1.0, 0.0, -1.0};
  const double complex candidate_overlap[] = {1.01, 0.0, 0.99};
  const double complex candidate_hamiltonian[] = {0.0, 0.99, 0.0};
  const double complex candidate_squared[] = {1.01, 0.0, 1.01};
  const double complex candidate_observable[] = {1.0, 0.01 * I, -1.0};
  const double scale = 1.0 / sqrt(2.0);
  const double complex coefficient[] = {scale, -scale};
  const double complex candidate_coefficient[] = {
      scale * 1.01, -scale * 0.99 + 0.002 * I};
  const MVMCKrylovFinalStateQuadraticModel baseline_model =
      model2(overlap, hamiltonian, squared, observable);
  const MVMCKrylovFinalStateQuadraticModel candidate_model =
      model2(candidate_overlap, candidate_hamiltonian, candidate_squared,
             candidate_observable);
  MVMCKrylovFinalStateSensitivityComparison identity;
  MVMCKrylovFinalStateSensitivityComparison coefficient_only;
  MVMCKrylovFinalStateSensitivityComparison matrix_only;
  MVMCKrylovFinalStateSensitivityComparison comparison;

  CHECK(mvmc_krylov_final_state_sensitivity_compare(
            &baseline_model, coefficient, &baseline_model, coefficient, 2,
            &identity) == MVMC_KRYLOV_STATUS_OK &&
            identity.valid &&
            close_double(identity.coefficient_projective_distance, 0.0,
                         2.0e-14) &&
            close_double(identity.overlap_relative_frobenius_error, 0.0,
                         0.0) &&
            close_double(identity.combined.energy_absolute, 0.0, 0.0) &&
            close_double(identity.combined.variance_absolute, 0.0, 0.0) &&
            close_double(identity.combined.observable_absolute, 0.0, 0.0),
        "identity sensitivity decomposition");
  CHECK(mvmc_krylov_final_state_sensitivity_compare(
            &baseline_model, coefficient, &baseline_model,
            candidate_coefficient, 2, &coefficient_only) ==
            MVMC_KRYLOV_STATUS_OK &&
            coefficient_only.valid &&
            coefficient_only.coefficient_projective_distance > 0.0 &&
            coefficient_only.overlap_relative_frobenius_error == 0.0 &&
            coefficient_only.hamiltonian_relative_frobenius_error == 0.0 &&
            coefficient_only.matrix_only.energy_absolute == 0.0 &&
            coefficient_only.matrix_only.variance_absolute == 0.0 &&
            close_double(coefficient_only.combined.energy_absolute,
                         coefficient_only.coefficient_only.energy_absolute,
                         2.0e-14) &&
            close_double(coefficient_only.combined.variance_absolute,
                         coefficient_only.coefficient_only.variance_absolute,
                         2.0e-14),
        "coefficient-only decomposition isolates matrix path");
  CHECK(mvmc_krylov_final_state_sensitivity_compare(
            &baseline_model, coefficient, &candidate_model, coefficient, 2,
            &matrix_only) == MVMC_KRYLOV_STATUS_OK &&
            matrix_only.valid &&
            matrix_only.coefficient_projective_distance == 0.0 &&
            matrix_only.coefficient_only.energy_absolute == 0.0 &&
            matrix_only.coefficient_only.variance_absolute == 0.0 &&
            matrix_only.overlap_relative_frobenius_error > 0.0 &&
            close_double(matrix_only.combined.energy_absolute,
                         matrix_only.matrix_only.energy_absolute,
                         2.0e-14) &&
            close_double(matrix_only.combined.variance_absolute,
                         matrix_only.matrix_only.variance_absolute,
                         2.0e-14),
        "matrix-only decomposition isolates coefficient path");
  CHECK(mvmc_krylov_final_state_sensitivity_compare(
            &baseline_model, coefficient, &candidate_model,
            candidate_coefficient, 2, &comparison) ==
            MVMC_KRYLOV_STATUS_OK &&
            comparison.valid &&
            comparison.coefficient_projective_distance > 0.0 &&
            comparison.overlap_relative_frobenius_error > 0.0 &&
            comparison.hamiltonian_relative_frobenius_error > 0.0 &&
            comparison.hamiltonian_squared_relative_frobenius_error > 0.0 &&
            comparison.observable_relative_frobenius_error > 0.0 &&
            comparison.coefficient_only.observable_absolute > 0.0 &&
            comparison.matrix_only.energy_absolute > 0.0 &&
            comparison.combined.variance_absolute > 0.0,
        "coefficient/matrix/combined sensitivity separation");
}

static void test_complex_dimension3(void) {
  const double complex overlap[] = {
      1.2, 0.1 + 0.05 * I, -0.03 + 0.02 * I,
      0.9, 0.04 - 0.01 * I, 1.1};
  const double complex hamiltonian[] = {
      -0.7, 0.2 + 0.1 * I, -0.05 * I,
      0.4, -0.15 + 0.03 * I, 1.3};
  const double complex squared[] = {
      1.5, -0.1 + 0.04 * I, 0.08 - 0.02 * I,
      1.1, 0.12 + 0.05 * I, 2.2};
  const double complex observable[] = {
      0.3, 0.07 - 0.02 * I, -0.04 + 0.01 * I,
      -0.5, 0.02 + 0.03 * I, 0.8};
  const double complex coefficient[] = {
      0.7 + 0.1 * I, -0.2 + 0.4 * I, 0.3 - 0.15 * I};
  MVMCKrylovFinalStateQuadraticModel model;
  MVMCKrylovFinalStateQuadraticMetrics metrics;
  CHECK(mvmc_krylov_final_state_quadratic_model_create(
            3, overlap, hamiltonian, squared, observable, 6, &model) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_final_state_quadratic_metrics(
                &model, coefficient, 3, &metrics) ==
                MVMC_KRYLOV_STATUS_OK &&
            metrics.valid && close_double(metrics.norm, 0.83605, 2.0e-14) &&
            close_double(metrics.energy, -0.23425632438251293,
                         2.0e-14) &&
            close_double(metrics.hamiltonian_second_moment,
                         1.4586448178936668, 2.0e-14) &&
            close_double(metrics.full_support_variance,
                         1.4037687923804618, 2.0e-14) &&
            close_double(metrics.observable, 0.15034985945816637,
                         2.0e-14),
        "complex dimension-3 fixed dense-oracle metrics");
}

static void test_optional_observable(void) {
  const double complex overlap[] = {1.0, 0.0, 1.0};
  const double complex hamiltonian[] = {0.0, 1.0, 0.0};
  const double complex squared[] = {1.0, 0.0, 1.0};
  const double complex coefficient[] = {1.0, -1.0};
  MVMCKrylovFinalStateQuadraticModel model;
  MVMCKrylovFinalStateQuadraticMetrics metrics;
  MVMCKrylovFinalStateSensitivityComparison comparison;

  CHECK(mvmc_krylov_final_state_quadratic_model_create(
            2, overlap, hamiltonian, squared, NULL, 3, &model) ==
            MVMC_KRYLOV_STATUS_OK &&
            model.valid && !model.has_observable &&
            mvmc_krylov_final_state_quadratic_metrics(
                &model, coefficient, 2, &metrics) ==
                MVMC_KRYLOV_STATUS_OK &&
            metrics.valid && !metrics.has_observable &&
            isnan(metrics.observable),
        "optional observable omitted from metrics");
  CHECK(mvmc_krylov_final_state_sensitivity_compare(
            &model, coefficient, &model, coefficient, 2, &comparison) ==
            MVMC_KRYLOV_STATUS_OK &&
            comparison.valid && !comparison.has_observable &&
            isnan(comparison.observable_relative_frobenius_error) &&
            isnan(comparison.combined.observable_absolute) &&
            isnan(comparison.combined.observable_scaled),
        "optional observable omitted from sensitivity comparison");
}

static void test_invalid_inputs(void) {
  const double complex overlap[] = {1.0, 0.0, 1.0};
  const double complex hamiltonian[] = {0.0, 1.0, 0.0};
  const double complex squared[] = {1.0, 0.0, 1.0};
  const double complex observable[] = {1.0, 0.0, -1.0};
  const double complex coefficient[] = {1.0, -1.0};
  MVMCKrylovFinalStateQuadraticModel model;
  MVMCKrylovFinalStateQuadraticMetrics metrics;
  MVMCKrylovFinalStateSensitivityComparison comparison;
  double complex bad_overlap[] = {1.0 + 1.0e-8 * I, 0.0, 1.0};
  double complex zero_overlap[] = {0.0, 0.0, 0.0};
  double complex bad_squared[] = {0.0, 0.0, 0.0};
  double complex bad_coefficient[] = {1.0, NAN + 0.0 * I};

  CHECK(mvmc_krylov_final_state_quadratic_model_create(
            2, bad_overlap, hamiltonian, squared, observable, 3, &model) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !model.valid,
        "non-real Hermitian diagonal rejected");
  CHECK(mvmc_krylov_final_state_quadratic_model_create(
            2, overlap, hamiltonian, squared, observable, 2, &model) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "upper count mismatch rejected");
  model = model2(zero_overlap, hamiltonian, squared, observable);
  CHECK(mvmc_krylov_final_state_quadratic_metrics(
            &model, coefficient, 2, &metrics) ==
            MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE &&
            !metrics.valid,
        "zero norm rejected");
  model = model2(overlap, hamiltonian, bad_squared, observable);
  CHECK(mvmc_krylov_final_state_quadratic_metrics(
            &model, coefficient, 2, &metrics) ==
            MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE,
        "materially negative full-support variance rejected");
  model = model2(overlap, hamiltonian, squared, observable);
  CHECK(mvmc_krylov_final_state_quadratic_metrics(
            &model, bad_coefficient, 2, &metrics) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "nonfinite coefficient rejected");
  {
    MVMCKrylovFinalStateQuadraticModel corrupt = model;
    corrupt.model_hash ^= UINT64_C(1);
    CHECK(mvmc_krylov_final_state_sensitivity_compare(
              &model, coefficient, &corrupt, coefficient, 2,
              &comparison) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              !comparison.valid,
          "tampered candidate model hash rejected");
  }
}

static void test_rank_invariance(void) {
#ifdef _mpi_use
  const double complex overlap[] = {1.0, 0.0, 1.0};
  const double complex hamiltonian[] = {0.0, 1.0, 0.0};
  const double complex squared[] = {1.0, 0.0, 1.0};
  const double complex observable[] = {1.0, 0.0, -1.0};
  const double scale = 1.0 / sqrt(2.0);
  const double complex coefficient[] = {scale, -scale};
  const MVMCKrylovFinalStateQuadraticModel model =
      model2(overlap, hamiltonian, squared, observable);
  MVMCKrylovFinalStateQuadraticMetrics metrics;
  uint64_t local[6];
  uint64_t minimum[6];
  uint64_t maximum[6];
  double fields[5];
  int index;
  CHECK(mvmc_krylov_final_state_quadratic_metrics(
            &model, coefficient, 2, &metrics) == MVMC_KRYLOV_STATUS_OK,
        "MPI sensitivity metrics");
  local[0] = model.model_hash;
  fields[0] = metrics.norm;
  fields[1] = metrics.energy;
  fields[2] = metrics.hamiltonian_second_moment;
  fields[3] = metrics.full_support_variance;
  fields[4] = metrics.observable;
  for (index = 0; index < 5; ++index) {
    memcpy(&local[index + 1], &fields[index], sizeof(local[index + 1]));
  }
  MPI_Allreduce(local, minimum, 6, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(local, maximum, 6, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  CHECK(memcmp(minimum, maximum, sizeof(minimum)) == 0,
        "sensitivity metrics differ by MPI rank");
#endif
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#else
  (void)argc;
  (void)argv;
#endif
  test_exact_metrics_and_invariance();
  test_sensitivity_decomposition();
  test_complex_dimension3();
  test_optional_observable();
  test_invalid_inputs();
  test_rank_invariance();
#ifdef _mpi_use
  {
    int global_failures = 0;
    MPI_Allreduce(&failures, &global_failures, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    failures = global_failures;
  }
#endif
  if (failures == 0 && world_rank == 0) {
    puts("krylov final-state sensitivity unit: PASS");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return failures == 0 ? 0 : 1;
}
