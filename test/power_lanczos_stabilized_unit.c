#include "power_lanczos_stabilized.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failures = 0;
static int rank = 0;
static int size = 1;

#define CHECK(condition, message)                                           \
  do {                                                                      \
    if (!(condition)) {                                                     \
      fprintf(stderr, "PowerLanczosStabilized FAIL rank %d: %s (line %d)\n", \
              rank, message, __LINE__);                                     \
      ++failures;                                                           \
    }                                                                       \
  } while (0)

static MVMCClassicPfaffianCommunicator world(void) {
#ifdef _mpi_use
  return MPI_COMM_WORLD;
#else
  return 0;
#endif
}

static void set_real_pair(double *slater, int up, int down, double value) {
  const int orbitals = 4;
  const int down_orbital = 2 + down;
  slater[down_orbital * orbitals + up] = -value;
  slater[up * orbitals + down_orbital] = value;
}

static void set_complex_pair(double complex *slater, int up, int down,
                             double complex value) {
  const int orbitals = 4;
  const int down_orbital = 2 + down;
  slater[down_orbital * orbitals + up] = -value;
  slater[up * orbitals + down_orbital] = value;
}

static void fill_real_slater(double *slater, int qp_total, double scale) {
  int qp;
  int up;
  int down;
  memset(slater, 0, (size_t)qp_total * 16 * sizeof(*slater));
  for (qp = 0; qp < qp_total; ++qp) {
    double *component = slater + (size_t)qp * 16;
    for (up = 0; up < 2; ++up) {
      for (down = 0; down < 2; ++down) {
        set_real_pair(component, up, down,
                      scale * (double)(qp + 1) *
                          (1.5 + (double)up + 0.25 * (double)down));
      }
    }
  }
}

static void fill_complex_slater(double complex *slater, int qp_total) {
  int qp;
  int up;
  int down;
  memset(slater, 0, (size_t)qp_total * 16 * sizeof(*slater));
  for (qp = 0; qp < qp_total; ++qp) {
    double complex *component = slater + (size_t)qp * 16;
    for (up = 0; up < 2; ++up) {
      for (down = 0; down < 2; ++down) {
        const double magnitude =
            (double)(qp + 1) *
            (1.5 + (double)up + 0.25 * (double)down);
        const double phase = 0.13 * (double)up - 0.21 * (double)down +
                             0.07 * (double)(up * down);
        set_complex_pair(component, up, down,
                         magnitude * (cos(phase) + I * sin(phase)));
      }
    }
  }
}

static void initialize_view(MVMCPowerLanczosClassicView *view,
                            MVMCPowerLanczosClassicArithmetic arithmetic,
                            double *slater_real,
                            double complex *slater_complex) {
  enum { qp_total = 8 };
  static int transfer_rows[4][4] = {
      {1, 0, 0, 0}, {0, 0, 1, 0},
      {1, 1, 0, 1}, {0, 1, 1, 1}};
  static int *transfer_indices[4] = {
      transfer_rows[0], transfer_rows[1],
      transfer_rows[2], transfer_rows[3]};
  static const double complex transfer_parameters[4] = {
      0.75, 0.75, 0.50, 0.50};
  static const int intra_indices[1] = {0};
  static const double intra_parameters[1] = {0.40};
  static double complex weights[qp_total];
  int qp;
  for (qp = 0; qp < qp_total; ++qp) {
    weights[qp] = 1.0 / (double)(qp + 1);
  }
  memset(view, 0, sizeof(*view));
  view->site_count = 2;
  view->up_electron_count = 1;
  view->down_electron_count = 1;
  view->arithmetic = arithmetic;
  view->transfer_count = 4;
  view->transfer_indices = transfer_indices;
  view->transfer_parameters = transfer_parameters;
  view->coulomb_intra_count = 1;
  view->coulomb_intra_indices = intra_indices;
  view->coulomb_intra_parameters = intra_parameters;
  view->qp_total = qp_total;
  view->qp_start = qp_total * rank / size;
  view->qp_end = qp_total * (rank + 1) / size;
  view->qp_weights = weights;
  view->slater_real = slater_real;
  view->slater_complex = slater_complex;
}

static MVMCPowerLanczosStabilizedResult run_case(
    int power_step, double scale,
    MVMCPowerLanczosClassicArithmetic arithmetic) {
  enum { qp_total = 8 };
  double slater_real[qp_total * 16];
  double complex slater_complex[qp_total * 16];
  MVMCPowerLanczosClassicView view;
  MVMCPowerLanczosStabilizedInput input;
  MVMCPowerLanczosStabilizedResult result;
  MVMCKrylovStatus status;
  double energy_imaginary_tolerance;
  double variance_imaginary_tolerance;
  if (arithmetic == MVMC_POWER_LANCZOS_CLASSIC_REAL) {
    fill_real_slater(slater_real, qp_total, scale);
    initialize_view(&view, arithmetic, slater_real, NULL);
  } else {
    fill_complex_slater(slater_complex, qp_total);
    initialize_view(&view, arithmetic, NULL, slater_complex);
  }
  memset(&input, 0, sizeof(input));
  input.classic_view = &view;
  input.communicator = world();
  input.power_step = power_step;
  input.seed = UINT64_C(0x53544142494c495a);
  input.warm_up = 16;
  input.sample_count = 256;
  input.interval = 1;
  status = mvmc_power_lanczos_stabilized_run(&input, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_OK, "stabilized run status");
  CHECK(result.valid && result.status == MVMC_KRYLOV_STATUS_OK,
        "stabilized result validity");
  CHECK(result.power_step == power_step &&
            result.coefficient_samples == 256 &&
            result.final_samples == 256,
        "stabilized result shape");
  CHECK(isfinite(result.energy) && isfinite(result.energy_standard_error) &&
            result.energy_standard_error >= 0.0 &&
            isfinite(result.variance) &&
            isfinite(result.variance_standard_error) &&
            result.variance_standard_error >= 0.0,
        "finite corrected energy and variance");
  CHECK(isfinite(result.condition_estimate) &&
            isfinite(result.gevp_residual) &&
            result.gevp_residual <= 1.0e-8,
        "finite GEVP diagnostics");
  energy_imaginary_tolerance =
      arithmetic == MVMC_POWER_LANCZOS_CLASSIC_REAL
          ? 1.0e-10
          : fmax(5.0e-3, 6.0 * result.energy_standard_error);
  variance_imaginary_tolerance =
      arithmetic == MVMC_POWER_LANCZOS_CLASSIC_REAL
          ? 1.0e-10
          : fmax(5.0e-3, 6.0 * result.variance_standard_error);
  if (!(fabs(result.final_energy_imaginary) <=
            energy_imaginary_tolerance &&
        fabs(result.variance_imaginary) <= variance_imaginary_tolerance &&
        isfinite(result.energy_tau_int) &&
        result.energy_tau_int >= 0.5 &&
        isfinite(result.effective_sample_count) &&
        result.effective_sample_count > 0.0)) {
    fprintf(stderr,
            "statistical diagnostics: energy_imag=%.17g "
            "variance_imag=%.17g tau=%.17g effective=%.17g\n",
            result.final_energy_imaginary, result.variance_imaginary,
            result.energy_tau_int, result.effective_sample_count);
    CHECK(0, "finite statistical diagnostics");
  }
#ifdef _mpi_use
  {
    double local[4] = {result.energy, result.variance,
                       result.energy_standard_error,
                       result.variance_standard_error};
    double minimum[4];
    double maximum[4];
    MPI_Allreduce(local, minimum, 4, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(local, maximum, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    CHECK(memcmp(minimum, maximum, sizeof(minimum)) == 0,
          "result differs between MPI ranks");
  }
#endif
  return result;
}

static void check_ed(const MVMCPowerLanczosStabilizedResult *result,
                     double expected_energy, double expected_variance,
                     const char *label) {
  const double energy_tolerance =
      fmax(5.0e-3, 6.0 * result->energy_standard_error);
  const double variance_tolerance =
      fmax(5.0e-3, 6.0 * result->variance_standard_error);
  if (fabs(result->energy - expected_energy) > energy_tolerance ||
      fabs(result->variance - expected_variance) > variance_tolerance) {
    fprintf(stderr,
            "%s ED mismatch: E=%.17g expected=%.17g tol=%.6g; "
            "variance=%.17g expected=%.17g tol=%.6g\n",
            label, result->energy, expected_energy, energy_tolerance,
            result->variance, expected_variance, variance_tolerance);
  }
  CHECK(fabs(result->energy - expected_energy) <= energy_tolerance,
        "corrected energy differs from exact diagonalization");
  CHECK(fabs(result->variance - expected_variance) <= variance_tolerance,
        "corrected variance differs from exact diagonalization");
}

int main(int argc, char **argv) {
  MVMCPowerLanczosStabilizedResult base;
  MVMCPowerLanczosStabilizedResult scaled;
  MVMCPowerLanczosStabilizedResult second_step;
  MVMCPowerLanczosStabilizedResult complex_first_step;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
  (void)argc;
  (void)argv;
#endif
  CHECK(mvmc_power_lanczos_stabilized_block_count(256) == 16,
        "valid sample count block selection");
  CHECK(mvmc_power_lanczos_stabilized_block_count(37) == 0,
        "invalid sample count block selection");
  base = run_case(1, 1.0, MVMC_POWER_LANCZOS_CLASSIC_REAL);
  check_ed(&base, -1.167342721287232, 0.0012031753879899743,
           "real first step");
  scaled = run_case(1, 1.0e150, MVMC_POWER_LANCZOS_CLASSIC_REAL);
  CHECK(fabs(base.energy - scaled.energy) <=
            1.0e-10 * fmax(1.0, fabs(base.energy)) &&
            fabs(base.variance - scaled.variance) <=
            1.0e-9 * fmax(1.0, fabs(base.variance)),
        "global amplitude scaling changes corrected observables");
  second_step = run_case(2, 1.0, MVMC_POWER_LANCZOS_CLASSIC_REAL);
  check_ed(&second_step, -1.1679611410882096, 0.000163604277214624,
           "real second step");
  complex_first_step =
      run_case(1, 1.0, MVMC_POWER_LANCZOS_CLASSIC_COMPLEX);
  check_ed(&complex_first_step, -1.1659572960326554,
           0.003155912326888499, "complex first step");
#ifdef _mpi_use
  {
    int total = 0;
    MPI_Allreduce(&failures, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    failures = total;
  }
  MPI_Finalize();
#endif
  if (failures != 0) return 1;
  if (rank == 0) printf("power-Lanczos stabilized unit: PASS\n");
  return 0;
}
