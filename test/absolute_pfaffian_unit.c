#include "absolute_pfaffian.h"

#include <complex.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failure_count = 0;

static MVMCAbsolutePfaffianResult component(MVMCPfaffianState state,
                                            double complex pfaffian);

#define CHECK(condition, message)                                      \
  do {                                                                 \
    if (!(condition)) {                                                \
      fprintf(stderr, "FAIL: %s (line %d)\n", (message), __LINE__);   \
      ++failure_count;                                                 \
    }                                                                  \
  } while (0)

static int close_complex(double complex actual, double complex expected,
                         double relative_tolerance) {
  return cabs(actual - expected) <=
         relative_tolerance * fmax(1.0, cabs(expected));
}

static void set_real_pair(double *matrix, int n, int row, int column,
                          double value) {
  matrix[row + column * n] = value;
  matrix[column + row * n] = -value;
}

static void set_complex_pair(double complex *matrix, int n, int row,
                             int column, double complex value) {
  matrix[row + column * n] = value;
  matrix[column + row * n] = -value;
}

/* Independent Laplace expansion, used only as a small dense oracle. */
static double complex reference_pfaffian(const double complex *matrix, int n) {
  double complex minor[64] = {0.0};
  double complex value = 0.0;
  int partner, old_column, old_row, new_column, new_row;

  if (n == 0) return 1.0;
  if (n == 2) return matrix[n];
  for (partner = 1; partner < n; ++partner) {
    new_column = 0;
    for (old_column = 1; old_column < n; ++old_column) {
      if (old_column == partner) continue;
      new_row = 0;
      for (old_row = 1; old_row < n; ++old_row) {
        if (old_row == partner) continue;
        minor[new_row + new_column * (n - 2)] =
            matrix[old_row + old_column * n];
        ++new_row;
      }
      ++new_column;
    }
    value += (partner % 2 == 1 ? 1.0 : -1.0) * matrix[partner * n] *
             reference_pfaffian(minor, n - 2);
  }
  return value;
}

static int real_matrix_is_value(const double *matrix, int n, double value) {
  int index;
  for (index = 0; index < n * n; ++index) {
    if (matrix[index] != value) return 0;
  }
  return 1;
}

static int complex_matrix_is_value(const double complex *matrix, int n,
                                   double complex value) {
  int index;
  for (index = 0; index < n * n; ++index) {
    if (matrix[index] != value) return 0;
  }
  return 1;
}

static void test_real_regular(void) {
  const int n = 4;
  double matrix[16] = {0.0};
  double original[16];
  double inverse[16];
  double complex oracle_matrix[16];
  MVMCAbsolutePfaffianResult result;
  MVMCPfaffianStatus status;
  int index;

  set_real_pair(matrix, n, 0, 1, 2.0);
  set_real_pair(matrix, n, 0, 2, -1.0);
  set_real_pair(matrix, n, 0, 3, 0.5);
  set_real_pair(matrix, n, 1, 2, 3.0);
  set_real_pair(matrix, n, 1, 3, 4.0);
  set_real_pair(matrix, n, 2, 3, -2.0);
  memcpy(original, matrix, sizeof(matrix));
  for (index = 0; index < 16; ++index) {
    inverse[index] = 1234.0;
    oracle_matrix[index] = matrix[index];
  }

  status = mvmc_absolute_pfaffian_real(matrix, n, n, inverse, n, 17, 0.0,
                                        0.0, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK, "real regular status");
  CHECK(result.state == MVMC_PFAFFIAN_REGULAR, "real regular state");
  CHECK(result.inverse_valid == 1, "real regular inverse valid");
  CHECK(result.factor_info == 0 && result.inverse_info == 0,
        "real regular backend info");
  CHECK(result.rebuild_generation == 17, "real generation preserved");
  CHECK(close_complex(result.pfaffian,
                      reference_pfaffian(oracle_matrix, n), 1.0e-13),
        "real Pfaffian matches independent dense oracle");
  CHECK(result.inverse_residual < 1.0e-13, "real inverse residual");
  CHECK(memcmp(matrix, original, sizeof(matrix)) == 0,
        "real source matrix remains unchanged");
}

static void test_complex_regular(void) {
  const int n = 4;
  double complex matrix[16] = {0.0};
  double complex original[16];
  double complex inverse[16];
  MVMCAbsolutePfaffianResult result;
  MVMCPfaffianStatus status;
  int index;

  set_complex_pair(matrix, n, 0, 1, 1.0 + 0.5 * I);
  set_complex_pair(matrix, n, 0, 2, -0.25 + 0.75 * I);
  set_complex_pair(matrix, n, 0, 3, 2.0 - 0.5 * I);
  set_complex_pair(matrix, n, 1, 2, 0.5 + 1.25 * I);
  set_complex_pair(matrix, n, 1, 3, -1.0 + 0.25 * I);
  set_complex_pair(matrix, n, 2, 3, 1.5 - 0.75 * I);
  memcpy(original, matrix, sizeof(matrix));
  for (index = 0; index < 16; ++index) inverse[index] = 1234.0 + I;

  status = mvmc_absolute_pfaffian_complex(
      matrix, n, n, inverse, n, 23, 0.0, 0.0, &result);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK, "complex regular status");
  CHECK(result.state == MVMC_PFAFFIAN_REGULAR, "complex regular state");
  CHECK(result.inverse_valid == 1, "complex regular inverse valid");
  CHECK(result.factor_info == 0 && result.inverse_info == 0,
        "complex regular backend info");
  CHECK(result.rebuild_generation == 23, "complex generation preserved");
  CHECK(close_complex(result.pfaffian, reference_pfaffian(matrix, n),
                      2.0e-13),
        "complex Pfaffian matches independent dense oracle");
  CHECK(result.inverse_residual < 1.0e-13, "complex inverse residual");
  CHECK(memcmp(matrix, original, sizeof(matrix)) == 0,
        "complex source matrix remains unchanged");
}

static void test_dense_six_by_six(void) {
  const int n = 6;
  double matrix[36] = {0.0};
  double inverse[36];
  double complex complex_matrix[36] = {0.0};
  double complex complex_inverse[36];
  MVMCAbsolutePfaffianResult result;
  int column, row;

  for (column = 1; column < n; ++column) {
    for (row = 0; row < column; ++row) {
      const double value = 0.17 * (double)(row + 1) +
                           0.31 * (double)(column + 1) +
                           0.07 * (double)(row * column + 1);
      set_real_pair(matrix, n, row, column, value);
      set_complex_pair(complex_matrix, n, row, column,
                       value + I * (0.11 * (double)(column - row)));
    }
  }
  CHECK(mvmc_absolute_pfaffian_real(matrix, n, n, inverse, n, 25, 0.0,
                                     0.0, &result) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "6x6 real status");
  CHECK(result.state == MVMC_PFAFFIAN_REGULAR, "6x6 real regular");
  for (column = 0; column < n * n; ++column) {
    complex_inverse[column] = matrix[column];
  }
  CHECK(close_complex(result.pfaffian,
                      reference_pfaffian(complex_inverse, n), 5.0e-12),
        "6x6 real Pfaffian matches recursive oracle");

  CHECK(mvmc_absolute_pfaffian_complex(
            complex_matrix, n, n, complex_inverse, n, 26, 0.0, 0.0,
            &result) == MVMC_PFAFFIAN_STATUS_OK,
        "6x6 complex status");
  CHECK(result.state == MVMC_PFAFFIAN_REGULAR, "6x6 complex regular");
  CHECK(close_complex(result.pfaffian,
                      reference_pfaffian(complex_matrix, n), 5.0e-12),
        "6x6 complex Pfaffian matches recursive oracle");
}

static double deterministic_value(uint64_t *state) {
  *state = *state * UINT64_C(6364136223846793005) +
           UINT64_C(1442695040888963407);
  return (double)((*state >> 11) & UINT64_C(0x1fffff)) /
             (double)UINT64_C(0x1fffff) -
         0.5;
}

static void test_dense_oracle_stress(void) {
  double matrix[64];
  double inverse[64];
  double complex complex_matrix[64];
  double complex complex_inverse[64];
  MVMCAbsolutePfaffianResult result;
  uint64_t random_state = UINT64_C(0x728f31a5c49dbe07);
  int n, sample, column, row;

  for (n = 2; n <= 8; n += 2) {
    for (sample = 0; sample < 16; ++sample) {
      memset(matrix, 0, sizeof(matrix));
      memset(complex_matrix, 0, sizeof(complex_matrix));
      for (column = 1; column < n; ++column) {
        for (row = 0; row < column; ++row) {
          double real_part = deterministic_value(&random_state);
          double imag_part = 0.5 * deterministic_value(&random_state);
          if ((row % 2) == 0 && column == row + 1) real_part += 2.0;
          set_real_pair(matrix, n, row, column, real_part);
          set_complex_pair(complex_matrix, n, row, column,
                           real_part + I * imag_part);
        }
      }

      CHECK(mvmc_absolute_pfaffian_real(
                matrix, n, n, inverse, n, (uint64_t)sample, 0.0, 0.0,
                &result) == MVMC_PFAFFIAN_STATUS_OK,
            "dense real stress status");
      CHECK(result.state == MVMC_PFAFFIAN_REGULAR,
            "dense real stress remains regular");
      for (column = 0; column < n * n; ++column) {
        complex_inverse[column] = matrix[column];
      }
      CHECK(close_complex(result.pfaffian,
                          reference_pfaffian(complex_inverse, n), 2.0e-11),
            "dense real stress Pfaffian sign matches oracle");

      CHECK(mvmc_absolute_pfaffian_complex(
                complex_matrix, n, n, complex_inverse, n,
                (uint64_t)sample, 0.0, 0.0, &result) ==
                MVMC_PFAFFIAN_STATUS_OK,
            "dense complex stress status");
      CHECK(result.state == MVMC_PFAFFIAN_REGULAR,
            "dense complex stress remains regular");
      CHECK(close_complex(result.pfaffian,
                          reference_pfaffian(complex_matrix, n), 3.0e-11),
            "dense complex stress Pfaffian phase matches oracle");
    }
  }
}

static void test_near_singular_and_scaling(void) {
  const int n = 4;
  double matrix[16] = {0.0};
  double scaled[16];
  double inverse[16];
  MVMCAbsolutePfaffianResult result, scaled_result;
  int index;

  set_real_pair(matrix, n, 0, 1, 1.0);
  set_real_pair(matrix, n, 2, 3, 1.0e-15);
  for (index = 0; index < 16; ++index) {
    scaled[index] = matrix[index] * 1.0e100;
    inverse[index] = 77.0;
  }
  CHECK(mvmc_absolute_pfaffian_real(matrix, n, n, inverse, n, 31, 0.0,
                                     0.0, &result) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "near-singular status");
  if (result.state != MVMC_PFAFFIAN_NEAR_SINGULAR) {
    fprintf(stderr,
            "near diagnostic: state=%s info=%d pivot=%.17g residual=%.17g "
            "pf=%.17g\n",
            mvmc_pfaffian_state_name(result.state), result.factor_info,
            result.scaled_min_pivot, result.inverse_residual,
            creal(result.pfaffian));
  }
  CHECK(result.state == MVMC_PFAFFIAN_NEAR_SINGULAR,
        "small scaled pivot is near singular");
  CHECK(result.inverse_valid == 0, "near-singular inverse invalid");
  CHECK(real_matrix_is_value(inverse, n, 77.0),
        "near-singular path does not touch inverse output");
  CHECK(close_complex(result.pfaffian, 1.0e-15, 1.0e-13),
        "near-singular Pfaffian remains available");

  for (index = 0; index < 16; ++index) inverse[index] = 88.0;
  CHECK(mvmc_absolute_pfaffian_real(scaled, n, n, inverse, n, 32, 0.0,
                                     0.0, &scaled_result) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "scaled near-singular status");
  if (scaled_result.state != MVMC_PFAFFIAN_NEAR_SINGULAR) {
    fprintf(stderr,
            "scaled near diagnostic: state=%s pivot=%.17g residual=%.17g "
            "pf=%.17g\n",
            mvmc_pfaffian_state_name(scaled_result.state),
            scaled_result.scaled_min_pivot, scaled_result.inverse_residual,
            creal(scaled_result.pfaffian));
  }
  CHECK(scaled_result.state == MVMC_PFAFFIAN_NEAR_SINGULAR,
        "near-singular classification is scale aware");
  CHECK(real_matrix_is_value(inverse, n, 88.0),
        "scaled near-singular path does not touch inverse output");
  CHECK(fabs(result.scaled_min_pivot - scaled_result.scaled_min_pivot) <
            1.0e-28,
        "scaled pivot metric is invariant under scalar rescaling");
}

static void test_singular_and_nonfinite(void) {
  const int n = 4;
  double matrix[16] = {0.0};
  double complex complex_matrix[16] = {0.0};
  double inverse[16];
  double complex complex_inverse[16];
  MVMCAbsolutePfaffianResult result;
  int index;

  set_real_pair(matrix, n, 0, 1, 1.0);
  for (index = 0; index < 16; ++index) inverse[index] = 99.0;
  CHECK(mvmc_absolute_pfaffian_real(matrix, n, n, inverse, n, 41, 0.0,
                                     0.0, &result) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "singular status");
  CHECK(result.state == MVMC_PFAFFIAN_SINGULAR, "singular state");
  CHECK(result.factor_info > 0, "PFAPACK singular pivot is retained");
  CHECK(result.pfaffian == 0.0, "singular Pfaffian is exactly zero");
  CHECK(result.inverse_valid == 0, "singular inverse invalid");
  CHECK(real_matrix_is_value(inverse, n, 99.0),
        "singular path performs zero inverse writes");

  set_complex_pair(complex_matrix, n, 0, 1, 1.0 + I);
  for (index = 0; index < 16; ++index) complex_inverse[index] = 6.0 + I;
  CHECK(mvmc_absolute_pfaffian_complex(
            complex_matrix, n, n, complex_inverse, n, 42, 0.0, 0.0,
            &result) == MVMC_PFAFFIAN_STATUS_OK,
        "complex singular status");
  CHECK(result.state == MVMC_PFAFFIAN_SINGULAR,
        "complex singular state");
  CHECK(result.pfaffian == 0.0 && result.inverse_valid == 0,
        "complex singular Pfaffian and inverse contract");
  CHECK(complex_matrix_is_value(complex_inverse, n, 6.0 + I),
        "complex singular path performs zero inverse writes");

  complex_matrix[0 + 1 * n] = NAN + 0.0 * I;
  complex_matrix[1 + 0 * n] = -NAN + 0.0 * I;
  for (index = 0; index < 16; ++index) complex_inverse[index] = 5.0 + I;
  CHECK(mvmc_absolute_pfaffian_complex(
            complex_matrix, n, n, complex_inverse, n, 42, 0.0, 0.0,
            &result) == MVMC_PFAFFIAN_STATUS_OK,
        "nonfinite input is a numerical classification");
  CHECK(result.state == MVMC_PFAFFIAN_NONFINITE, "nonfinite state");
  CHECK(result.inverse_valid == 0, "nonfinite inverse invalid");
  CHECK(complex_matrix_is_value(complex_inverse, n, 5.0 + I),
        "nonfinite path does not touch inverse output");
}

static void test_singular_inverse_guard_page(void) {
  const int n = 4;
  double matrix[16] = {0.0};
  MVMCAbsolutePfaffianResult result;
  const long page_size = sysconf(_SC_PAGESIZE);
  void *guard_page;

  CHECK(page_size > 0, "system page size is available");
  if (page_size <= 0) return;
  guard_page = mmap(NULL, (size_t)page_size, PROT_NONE,
                    MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  CHECK(guard_page != MAP_FAILED, "inverse guard page allocation");
  if (guard_page == MAP_FAILED) return;

  set_real_pair(matrix, n, 0, 1, 1.0);
  CHECK(mvmc_absolute_pfaffian_real(matrix, n, n, (double *)guard_page, n,
                                     43, 0.0, 0.0, &result) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "singular evaluation succeeds with inaccessible inverse output");
  CHECK(result.state == MVMC_PFAFFIAN_SINGULAR &&
            result.inverse_valid == 0,
        "singular guard-page state");
  CHECK(munmap(guard_page, (size_t)page_size) == 0,
        "inverse guard page cleanup");
}

static void test_state_transitions(void) {
  const int n = 4;
  double matrix[16] = {0.0};
  double inverse[16];
  MVMCAbsolutePfaffianResult result;
  const MVMCPfaffianState expected[] = {
      MVMC_PFAFFIAN_REGULAR, MVMC_PFAFFIAN_REGULAR,
      MVMC_PFAFFIAN_SINGULAR, MVMC_PFAFFIAN_REGULAR,
      MVMC_PFAFFIAN_SINGULAR, MVMC_PFAFFIAN_SINGULAR};
  const double second_pivot[] = {1.0, 0.5, 0.0, 0.25, 0.0, 0.0};
  int step;

  for (step = 0; step < 6; ++step) {
    memset(matrix, 0, sizeof(matrix));
    set_real_pair(matrix, n, 0, 1, 1.0);
    if (second_pivot[step] != 0.0) {
      set_real_pair(matrix, n, 2, 3, second_pivot[step]);
    }
    CHECK(mvmc_absolute_pfaffian_real(
              matrix, n, n, inverse, n, (uint64_t)(100 + step), 0.0, 0.0,
              &result) == MVMC_PFAFFIAN_STATUS_OK,
          "state transition evaluation status");
    CHECK(result.state == expected[step],
          "R-to-R, R-to-S, S-to-R, or S-to-S transition state");
    CHECK(result.rebuild_generation == (uint64_t)(100 + step),
          "transition generation increments");
  }
}

static void test_invalid_arguments(void) {
  double matrix[16] = {0.0};
  double inverse[16] = {0.0};
  MVMCAbsolutePfaffianResult result;

  CHECK(mvmc_absolute_pfaffian_real(matrix, 0, 0, inverse, 0, 0, 0.0,
                                     0.0, &result) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "zero dimension rejected");
  CHECK(mvmc_absolute_pfaffian_real(matrix, -2, 1, inverse, 1, 0, 0.0,
                                     0.0, &result) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "negative dimension rejected");
  CHECK(mvmc_absolute_pfaffian_real(matrix, 3, 3, inverse, 3, 0, 0.0,
                                     0.0, &result) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "odd dimension rejected");
  CHECK(mvmc_absolute_pfaffian_real(matrix, 4, 3, inverse, 4, 0, 0.0,
                                     0.0, &result) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "short source leading dimension rejected");
  CHECK(mvmc_absolute_pfaffian_real(matrix, INT_MAX, INT_MAX, inverse,
                                     INT_MAX, 0, 0.0, 0.0, &result) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "allocation-size overflow rejected before dereference");
  CHECK(mvmc_absolute_pfaffian_real(matrix, 4, 4, inverse, 4, 0, NAN, 0.0,
                                     &result) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "nonfinite pivot tolerance rejected");

  set_real_pair(matrix, 4, 0, 1, 1.0);
  matrix[1 + 0 * 4] = 1.0;
  CHECK(mvmc_absolute_pfaffian_real(matrix, 4, 4, inverse, 4, 0, 0.0,
                                     0.0, &result) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "non-skew source rejected");

  {
    MVMCAbsolutePfaffianResult components[2];
    double complex weights[2] = {1.0, 1.0};
    MVMCProjectedAmplitudeResult projected;
    components[0] = component(MVMC_PFAFFIAN_REGULAR, 1.0);
    components[1] = component(MVMC_PFAFFIAN_REGULAR, 2.0);
    CHECK(mvmc_projected_amplitude_slice(components, 2, weights, 0, 0, 2,
                                          &projected) ==
              MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
          "zero NQPFull rejected");
    CHECK(mvmc_projected_amplitude_slice(components, 2, weights, 2, -1, 1,
                                          &projected) ==
              MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
          "negative local QP start rejected");
    CHECK(mvmc_projected_amplitude_slice(components, 2, weights, 2, 0, 3,
                                          &projected) ==
              MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
          "local QP end beyond NQPFull rejected");
    CHECK(mvmc_projected_amplitude_slice(components, 2, weights, 2, 1, 1,
                                          &projected) ==
              MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
          "empty local QP range rejected");
    CHECK(mvmc_projected_amplitude_slice(components, 1, weights, 2, 0, 2,
                                          &projected) ==
              MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
          "inconsistent local count and QP range rejected");
    CHECK(projected.valid == 0, "invalid slice result is marked invalid");
  }
}

static MVMCAbsolutePfaffianResult component(MVMCPfaffianState state,
                                            double complex pfaffian) {
  MVMCAbsolutePfaffianResult result;
  memset(&result, 0, sizeof(result));
  result.state = state;
  result.pfaffian = pfaffian;
  return result;
}

static void test_projected_amplitude(void) {
  MVMCAbsolutePfaffianResult components[4];
  double complex weights[4];
  MVMCProjectedAmplitudeResult result;

  components[0] = component(MVMC_PFAFFIAN_REGULAR, 2.0 + 3.0 * I);
  weights[0] = 0.5 - I;
  CHECK(mvmc_projected_amplitude(components, weights, 1, &result) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "NQP=1 aggregation status");
  CHECK(result.valid == 1, "successful aggregation result is valid");
  CHECK(close_complex(result.total, 4.0 - 0.5 * I, 1.0e-15),
        "NQP=1 weighted amplitude");
  CHECK(result.regular_count == 1 && result.singular_count == 0,
        "NQP=1 state counters");

  components[0] = component(MVMC_PFAFFIAN_REGULAR, 2.0);
  components[1] = component(MVMC_PFAFFIAN_SINGULAR, 0.0);
  components[2] = component(MVMC_PFAFFIAN_NEAR_SINGULAR, -1.0 + I);
  components[3] = component(MVMC_PFAFFIAN_REGULAR, 3.0);
  weights[0] = 1.0;
  weights[1] = 7.0 - I;
  weights[2] = I;
  weights[3] = -0.5;
  CHECK(mvmc_projected_amplitude(components, weights, 4, &result) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "NQP=4 aggregation status");
  CHECK(close_complex(result.total, -0.5 - I, 1.0e-15),
        "one singular component does not zero total amplitude");
  CHECK(result.regular_count == 2 && result.near_singular_count == 1 &&
            result.singular_count == 1,
        "NQP=4 state counters");

  components[0] = component(MVMC_PFAFFIAN_REGULAR, 1.0);
  components[1] = component(MVMC_PFAFFIAN_REGULAR, -1.0);
  weights[0] = weights[1] = 1.0;
  CHECK(mvmc_projected_amplitude(components, weights, 2, &result) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "cancellation aggregation status");
  CHECK(result.total == 0.0 && result.sum_abs == 2.0 &&
            result.cancellation_ratio == 0.0,
        "exact cancellation metrics");

  components[1] =
      component(MVMC_PFAFFIAN_REGULAR, -(1.0 - 1.0e-14));
  CHECK(mvmc_projected_amplitude(components, weights, 2, &result) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "near-cancellation aggregation status");
  CHECK(cabs(result.total) < 1.1e-14 &&
            result.cancellation_ratio < 6.0e-15,
        "near-zero cancellation metric remains resolved");

  components[0] = component(MVMC_PFAFFIAN_SINGULAR, 0.0);
  components[1] = component(MVMC_PFAFFIAN_SINGULAR, 0.0);
  CHECK(mvmc_projected_amplitude(components, weights, 2, &result) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "all-zero aggregation status");
  CHECK(result.total == 0.0 && result.sum_abs == 0.0 &&
            result.singular_count == 2,
        "all-singular total amplitude is zero");

  {
    MVMCAbsolutePfaffianResult compensated_components[6];
    double complex compensated_weights[6];
    const double values[6] = {1.0e16, 1.0, 1.0, 1.0, 1.0, -1.0e16};
    int index;

    for (index = 0; index < 6; ++index) {
      compensated_components[index] =
          component(MVMC_PFAFFIAN_REGULAR, values[index]);
      compensated_weights[index] = 1.0;
    }
    CHECK(mvmc_projected_amplitude(compensated_components,
                                   compensated_weights, 6, &result) ==
              MVMC_PFAFFIAN_STATUS_OK,
          "compensated aggregation status");
    CHECK(result.valid == 1 && result.total == 4.0,
          "compensated sum recovers terms lost by plain summation");
  }

  components[0] = component(MVMC_PFAFFIAN_REGULAR, 1.0);
  components[1] = component(MVMC_PFAFFIAN_INVALID, 2.0);
  weights[0] = weights[1] = 1.0;
  CHECK(mvmc_projected_amplitude(components, weights, 2, &result) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "invalid component rejects aggregation");
  CHECK(result.valid == 0 && result.total == 0.0 &&
            result.regular_count == 0,
        "failed aggregation does not expose a partial result");

  components[0] = component(MVMC_PFAFFIAN_REGULAR, DBL_MAX);
  components[1] = component(MVMC_PFAFFIAN_REGULAR, DBL_MAX);
  weights[0] = weights[1] = 1.0;
  CHECK(mvmc_projected_amplitude(components, weights, 2, &result) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "finite terms whose aggregate overflows are rejected");
  CHECK(result.valid == 0 && result.total == 0.0,
        "overflowing aggregate does not expose a partial result");
}

static void test_rank_local_projected_amplitude(void) {
  MVMCAbsolutePfaffianResult components[4];
  double complex weights[4] = {1.0, 2.0, -I, 0.5 + 0.25 * I};
  MVMCProjectedAmplitudeResult local, expected;
  double complex reduced_total;
  double reduced_sum_abs;
  unsigned long long local_counts[3], reduced_counts[3];
  int rank = 0, size = 1;
  int qp_start, qp_end;

#ifdef _mpi_use
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  if (size != 1 && size != 2) return;

  components[0] = component(MVMC_PFAFFIAN_REGULAR, 2.0);
  components[1] = component(MVMC_PFAFFIAN_SINGULAR, 0.0);
  components[2] = component(MVMC_PFAFFIAN_NEAR_SINGULAR, -1.0 + I);
  components[3] = component(MVMC_PFAFFIAN_REGULAR, 3.0 - 0.5 * I);
  qp_start = 4 * rank / size;
  qp_end = 4 * (rank + 1) / size;
  CHECK(mvmc_projected_amplitude_slice(
            components + qp_start, (size_t)(qp_end - qp_start), weights, 4,
            qp_start, qp_end, &local) == MVMC_PFAFFIAN_STATUS_OK,
        "rank-local projected amplitude status");
  CHECK(mvmc_projected_amplitude(components, weights, 4, &expected) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "full projected amplitude reference status");

  reduced_total = local.total;
  reduced_sum_abs = local.sum_abs;
  local_counts[0] = (unsigned long long)local.regular_count;
  local_counts[1] = (unsigned long long)local.near_singular_count;
  local_counts[2] = (unsigned long long)local.singular_count;
  memcpy(reduced_counts, local_counts, sizeof(local_counts));
#ifdef _mpi_use
  if (size > 1) {
    MPI_Allreduce(&local.total, &reduced_total, 1, MPI_C_DOUBLE_COMPLEX, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&local.sum_abs, &reduced_sum_abs, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(local_counts, reduced_counts, 3, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  CHECK(close_complex(reduced_total, expected.total, 1.0e-15),
        "rank-local totals reduce to full amplitude");
  CHECK(fabs(reduced_sum_abs - expected.sum_abs) < 1.0e-14,
        "rank-local absolute sums reduce to full metric");
  CHECK(reduced_counts[0] == expected.regular_count &&
            reduced_counts[1] == expected.near_singular_count &&
            reduced_counts[2] == expected.singular_count,
        "rank-local state counters reduce to full counters");
}

static void test_value_only_api(void) {
  const int n = 4;
  double matrix[16] = {0.0};
  double original[16];
  double inverse[16];
  double complex complex_matrix[16] = {0.0};
  double complex complex_inverse[16];
  MVMCAbsolutePfaffianRealValueWorkspace *real_workspace = NULL;
  MVMCAbsolutePfaffianComplexValueWorkspace *complex_workspace = NULL;
  MVMCAbsolutePfaffianResult full;
  MVMCAbsolutePfaffianValueResult value;
  MVMCAbsolutePfaffianValueResult components[4];
  MVMCProjectedAmplitudeResult aggregate, empty;
  double complex weights[4] = {1.0, 2.0, -I, 0.5 + 0.25 * I};

  set_real_pair(matrix, n, 0, 1, 2.0);
  set_real_pair(matrix, n, 0, 2, -1.0);
  set_real_pair(matrix, n, 0, 3, 0.5);
  set_real_pair(matrix, n, 1, 2, 3.0);
  set_real_pair(matrix, n, 1, 3, 4.0);
  set_real_pair(matrix, n, 2, 3, -2.0);
  memcpy(original, matrix, sizeof(matrix));
  CHECK(mvmc_absolute_pfaffian_real(matrix, n, n, inverse, n, 501, 0.0,
                                     0.0, &full) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "value-only real full reference status");
  CHECK(mvmc_absolute_pfaffian_real_value(matrix, n, n, 0.0, &value) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "value-only real convenience status");
  CHECK(value.state == MVMC_PFAFFIAN_VALUE_WELL_PIVOTED,
        "value-only real well-pivoted classification");
  CHECK(value.factor_info == full.factor_info &&
            value.matrix_scale == full.matrix_scale &&
            value.scaled_min_pivot == full.scaled_min_pivot &&
            value.pfaffian == full.pfaffian,
        "value-only real diagnostics and Pfaffian match full P1");
  CHECK(memcmp(matrix, original, sizeof(matrix)) == 0,
        "value-only real leaves source unchanged");

  CHECK(mvmc_absolute_pfaffian_real_value_workspace_create(
            n, &real_workspace) == MVMC_PFAFFIAN_STATUS_OK,
        "value-only real workspace create");
  CHECK(mvmc_absolute_pfaffian_real_value_with_workspace(
            real_workspace, matrix, n, n, 0.0, &value) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            value.pfaffian == full.pfaffian,
        "value-only real workspace reuse first evaluation");
  CHECK(mvmc_absolute_pfaffian_real_value_with_workspace(
            real_workspace, matrix, n, n, 0.0, &value) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            value.pfaffian == full.pfaffian,
        "value-only real workspace reuse second evaluation");

  memset(matrix, 0, sizeof(matrix));
  set_real_pair(matrix, n, 0, 1, 1.0);
  set_real_pair(matrix, n, 2, 3, 1.0e-15);
  CHECK(mvmc_absolute_pfaffian_real_value_with_workspace(
            real_workspace, matrix, n, n, 0.0, &value) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            value.state == MVMC_PFAFFIAN_VALUE_NEAR_PIVOT &&
            close_complex(value.pfaffian, 1.0e-15, 1.0e-13),
        "value-only near pivot retains finite Pfaffian");
  memset(matrix, 0, sizeof(matrix));
  set_real_pair(matrix, n, 0, 1, 1.0);
  CHECK(mvmc_absolute_pfaffian_real_value_with_workspace(
            real_workspace, matrix, n, n, 0.0, &value) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            value.state == MVMC_PFAFFIAN_VALUE_SINGULAR &&
            value.factor_info > 0 && value.pfaffian == 0.0,
        "value-only singular contract");
  matrix[0] = NAN;
  CHECK(mvmc_absolute_pfaffian_real_value_with_workspace(
            real_workspace, matrix, n, n, 0.0, &value) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            value.state == MVMC_PFAFFIAN_VALUE_NONFINITE,
        "value-only nonfinite classification");
  CHECK(mvmc_absolute_pfaffian_real_value_with_workspace(
            real_workspace, matrix, n, n, NAN, &value) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            value.state == MVMC_PFAFFIAN_VALUE_INVALID,
        "value-only invalid input is atomic");

  set_complex_pair(complex_matrix, n, 0, 1, 1.0 + 0.5 * I);
  set_complex_pair(complex_matrix, n, 0, 2, -0.25 + 0.75 * I);
  set_complex_pair(complex_matrix, n, 0, 3, 2.0 - 0.5 * I);
  set_complex_pair(complex_matrix, n, 1, 2, 0.5 + 1.25 * I);
  set_complex_pair(complex_matrix, n, 1, 3, -1.0 + 0.25 * I);
  set_complex_pair(complex_matrix, n, 2, 3, 1.5 - 0.75 * I);
  CHECK(mvmc_absolute_pfaffian_complex(
            complex_matrix, n, n, complex_inverse, n, 502, 0.0, 0.0,
            &full) == MVMC_PFAFFIAN_STATUS_OK,
        "value-only complex full reference status");
  CHECK(mvmc_absolute_pfaffian_complex_value_workspace_create(
            n, &complex_workspace) == MVMC_PFAFFIAN_STATUS_OK,
        "value-only complex workspace create");
  CHECK(mvmc_absolute_pfaffian_complex_value_with_workspace(
            complex_workspace, complex_matrix, n, n, 0.0, &value) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            value.state == MVMC_PFAFFIAN_VALUE_WELL_PIVOTED &&
            close_complex(value.pfaffian, full.pfaffian, 1.0e-15),
        "value-only complex phase matches full P1");
  memset(complex_matrix, 0, sizeof(complex_matrix));
  set_complex_pair(complex_matrix, n, 0, 1, 1.0 + I);
  set_complex_pair(complex_matrix, n, 2, 3, 1.0e-15 - 0.5e-15 * I);
  CHECK(mvmc_absolute_pfaffian_complex(
            complex_matrix, n, n, complex_inverse, n, 503, 0.0, 0.0,
            &full) == MVMC_PFAFFIAN_STATUS_OK &&
            full.state == MVMC_PFAFFIAN_NEAR_SINGULAR,
        "value-only complex near full reference");
  CHECK(mvmc_absolute_pfaffian_complex_value_with_workspace(
            complex_workspace, complex_matrix, n, n, 0.0, &value) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            value.state == MVMC_PFAFFIAN_VALUE_NEAR_PIVOT &&
            value.pfaffian == full.pfaffian,
        "value-only complex near state and phase match full P1");
  memset(complex_matrix, 0, sizeof(complex_matrix));
  set_complex_pair(complex_matrix, n, 0, 1, 1.0 + I);
  CHECK(mvmc_absolute_pfaffian_complex(
            complex_matrix, n, n, complex_inverse, n, 504, 0.0, 0.0,
            &full) == MVMC_PFAFFIAN_STATUS_OK &&
            full.state == MVMC_PFAFFIAN_SINGULAR,
        "value-only complex singular full reference");
  CHECK(mvmc_absolute_pfaffian_complex_value_with_workspace(
            complex_workspace, complex_matrix, n, n, 0.0, &value) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            value.state == MVMC_PFAFFIAN_VALUE_SINGULAR &&
            value.pfaffian == full.pfaffian && value.pfaffian == 0.0,
        "value-only complex singular state matches full P1");
  complex_matrix[0] = INFINITY + 0.0 * I;
  CHECK(mvmc_absolute_pfaffian_complex_value_with_workspace(
            complex_workspace, complex_matrix, n, n, 0.0, &value) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            value.state == MVMC_PFAFFIAN_VALUE_NONFINITE,
        "value-only complex nonfinite classification");
  CHECK(mvmc_absolute_pfaffian_complex_value_with_workspace(
            complex_workspace, complex_matrix, n, n, INFINITY, &value) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            value.state == MVMC_PFAFFIAN_VALUE_INVALID,
        "value-only complex invalid input is atomic");

  components[0].state = MVMC_PFAFFIAN_VALUE_WELL_PIVOTED;
  components[0].pfaffian = 2.0;
  components[1].state = MVMC_PFAFFIAN_VALUE_SINGULAR;
  components[1].pfaffian = 0.0;
  components[2].state = MVMC_PFAFFIAN_VALUE_NEAR_PIVOT;
  components[2].pfaffian = -1.0 + I;
  components[3].state = MVMC_PFAFFIAN_VALUE_WELL_PIVOTED;
  components[3].pfaffian = 3.0 - 0.5 * I;
  CHECK(mvmc_projected_amplitude_values(
            components, weights, 4, &aggregate) == MVMC_PFAFFIAN_STATUS_OK &&
            aggregate.valid && aggregate.regular_count == 2 &&
            aggregate.near_singular_count == 1 &&
            aggregate.singular_count == 1,
        "value-only projected aggregation preserves classifications");
  CHECK(close_complex(aggregate.total, 4.625 + 1.5 * I, 1.0e-15),
        "value-only projected total matches direct oracle");
  CHECK(fabs(aggregate.sum_abs -
             (2.0 + sqrt(2.0) + cabs(1.625 + 0.5 * I))) < 1.0e-14,
        "value-only projected absolute sum matches direct oracle");
  CHECK(fabs(aggregate.cancellation_ratio -
             cabs(aggregate.total) / aggregate.sum_abs) < 1.0e-15,
        "value-only projected cancellation metric matches oracle");
  CHECK(mvmc_projected_amplitude_value_slice(
            NULL, 0, weights, 4, 2, 2, &empty) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            empty.valid && empty.total == 0.0 && empty.regular_count == 0,
        "value-only empty QP slice is valid");
  components[1].pfaffian = 1.0;
  CHECK(mvmc_projected_amplitude_values(
            components, weights, 4, &aggregate) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            aggregate.valid == 0,
        "value-only invalid singular component fails atomically");

  components[0].state = MVMC_PFAFFIAN_VALUE_SINGULAR;
  components[0].pfaffian = 0.0;
  components[1] = components[0];
  CHECK(mvmc_projected_amplitude_values(
            components, weights, 2, &aggregate) == MVMC_PFAFFIAN_STATUS_OK &&
            aggregate.valid && aggregate.total == 0.0 &&
            aggregate.sum_abs == 0.0 && aggregate.singular_count == 2,
        "all-zero value components remain a valid projected total");
  components[0].state = MVMC_PFAFFIAN_VALUE_WELL_PIVOTED;
  components[0].pfaffian = 1.0e16;
  components[1].state = MVMC_PFAFFIAN_VALUE_WELL_PIVOTED;
  components[1].pfaffian = -1.0e16;
  components[2].state = MVMC_PFAFFIAN_VALUE_WELL_PIVOTED;
  components[2].pfaffian = 1.0;
  weights[0] = weights[1] = weights[2] = 1.0;
  CHECK(mvmc_projected_amplitude_values(
            components, weights, 3, &aggregate) == MVMC_PFAFFIAN_STATUS_OK &&
            aggregate.total == 1.0,
        "value-only compensated sum retains a small cancellation remainder");

  mvmc_absolute_pfaffian_complex_value_workspace_destroy(complex_workspace);
  mvmc_absolute_pfaffian_real_value_workspace_destroy(real_workspace);
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  if (!mvmc_absolute_pfaffian_strict_fp_enabled()) {
    fputs("FAIL: absolute Pfaffian core requires strict floating-point "
          "semantics\n",
          stderr);
#ifdef _mpi_use
    MPI_Finalize();
#endif
    return 1;
  }
  test_real_regular();
  test_complex_regular();
  test_dense_six_by_six();
  test_dense_oracle_stress();
  test_near_singular_and_scaling();
  test_singular_and_nonfinite();
  test_singular_inverse_guard_page();
  test_state_transitions();
  test_invalid_arguments();
  test_projected_amplitude();
  test_rank_local_projected_amplitude();
  test_value_only_api();

  if (failure_count != 0) {
    fprintf(stderr, "%d absolute Pfaffian checks failed\n", failure_count);
#ifdef _mpi_use
    MPI_Finalize();
#endif
    return 1;
  }
  if (failure_count == 0) puts("absolute Pfaffian unit checks passed");
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return 0;
}
