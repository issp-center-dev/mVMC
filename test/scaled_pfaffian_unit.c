#include "absolute_pfaffian.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failures = 0;

#define CHECK(condition, message)                                             \
  do {                                                                        \
    if (!(condition)) {                                                       \
      fprintf(stderr, "ScaledPfaffian_Unit FAIL: %s\n", (message));         \
      ++failures;                                                             \
    }                                                                         \
  } while (0)

static int close_double(double actual, double expected, double tolerance) {
  return fabs(actual - expected) <=
         tolerance * fmax(1.0, fabs(expected));
}

static int close_complex(double complex actual, double complex expected,
                         double tolerance) {
  return cabs(actual - expected) <=
         tolerance * fmax(1.0, cabs(expected));
}

static void set_real_pair(double *matrix, int n, int left, int right,
                          double value) {
  matrix[left + (size_t)right * (size_t)n] = value;
  matrix[right + (size_t)left * (size_t)n] = -value;
}

static void set_complex_pair(double complex *matrix, int n, int left,
                             int right, double complex value) {
  matrix[left + (size_t)right * (size_t)n] = value;
  matrix[right + (size_t)left * (size_t)n] = -value;
}

static double complex dense_pfaffian_recursive(
    const double complex *matrix, int n) {
  double complex total = 0.0;
  int partner;
  if (n == 0) return 1.0;
  for (partner = 1; partner < n; ++partner) {
    double complex minor[36] = {0.0};
    int source_column, source_row, target_column = 0;
    for (source_column = 1; source_column < n; ++source_column) {
      int target_row = 0;
      if (source_column == partner) continue;
      for (source_row = 1; source_row < n; ++source_row) {
        if (source_row == partner) continue;
        minor[target_row + (size_t)target_column * (size_t)(n - 2)] =
            matrix[source_row + (size_t)source_column * (size_t)n];
        ++target_row;
      }
      ++target_column;
    }
    total += (partner % 2 == 1 ? 1.0 : -1.0) *
             matrix[(size_t)partner * (size_t)n] *
             dense_pfaffian_recursive(minor, n - 2);
  }
  return total;
}

static MVMCScaledComplex finite_value(double complex phase, double log_abs) {
  MVMCScaledComplex value;
  memset(&value, 0, sizeof(value));
  CHECK(mvmc_scaled_complex_make_finite(phase, log_abs, -INFINITY, &value) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "finite scaled fixture construction");
  return value;
}

static MVMCAbsolutePfaffianScaledValueResult scaled_component(
    MVMCPfaffianValueState factor_state, MVMCScaledComplex value) {
  MVMCAbsolutePfaffianScaledValueResult result;
  memset(&result, 0, sizeof(result));
  result.factor_state = factor_state;
  result.value = value;
  return result;
}

static void test_scaled_pure_math(void) {
  MVMCScaledComplex values[6];
  MVMCScaledComplex product;
  MVMCScaledComplex sum;
  MVMCScaledComplex sentinel;
  double complex exported = 123.0 + 7.0 * I;

  CHECK(mvmc_scaled_complex_make_finite(2.0, 1000.0, -INFINITY,
                                        &values[0]) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            values[0].state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
            close_complex(values[0].phase, 1.0, 1.0e-15) &&
            values[0].log_abs == 1000.0,
        "finite constructor normalizes phase and preserves log magnitude");
  CHECK(mvmc_scaled_complex_export_common_scale(&values[0], 1000.0,
                                                 &exported) ==
            MVMC_SCALED_EXPORT_OK &&
            close_complex(exported, 1.0, 1.0e-15),
        "common-scale export avoids overflow");
  CHECK(mvmc_scaled_complex_export_common_scale(&values[0], 0.0,
                                                 &exported) ==
            MVMC_SCALED_EXPORT_OVERFLOW,
        "overflow export is classified separately");

  values[1] = finite_value(-I, -1000.0);
  CHECK(mvmc_scaled_complex_export_common_scale(&values[1], 0.0,
                                                 &exported) ==
            MVMC_SCALED_EXPORT_UNDERFLOW &&
            exported == 0.0,
        "underflow export is not exact zero");
  CHECK(mvmc_scaled_complex_multiply(&values[0], &values[1], &product) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            product.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
            close_double(product.log_abs, 0.0, 1.0e-15) &&
            close_complex(product.phase, -I, 1.0e-15) &&
            isfinite(product.log_abs_error_bound),
        "scaled multiplication handles exponent span 2000");

  CHECK(mvmc_scaled_complex_from_raw(0.0, &values[2]) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            values[2].state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO &&
            isfinite(values[2].log_abs_error_bound),
        "raw zero imports as numeric zero");
  CHECK(mvmc_scaled_complex_from_raw(NAN + I, &values[3]) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            values[3].state == MVMC_SCALED_COMPLEX_NONFINITE,
        "raw NaN imports as nonfinite");
  CHECK(mvmc_scaled_complex_from_raw(
            DBL_MAX + I * DBL_MAX, &values[3]) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            values[3].state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
            values[3].log_abs > log(DBL_MAX) &&
            close_complex(values[3].phase, (1.0 + I) / sqrt(2.0),
                          1.0e-15),
        "finite complex raw input survives cabs-range overflow");
  CHECK(mvmc_scaled_complex_make_exact_zero(&values[4]) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            values[4].state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
        "structural exact zero constructor");
  CHECK(mvmc_scaled_complex_multiply(&values[4], &values[0], &product) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            product.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
        "exact zero product remains exact");

  CHECK(mvmc_scaled_complex_make_numeric_zero(
            -20.0, 0.0, -INFINITY, 0.0, &values[2]) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "numeric zero fixture construction");
  values[3] = finite_value(1.0, 3.0);
  values[3].log_abs_error_bound = -30.0;
  CHECK(mvmc_scaled_complex_multiply(&values[2], &values[3], &product) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            product.state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO &&
            product.log_abs_error_bound >= -17.0,
        "numeric-zero product propagates both operand error bounds");
  CHECK(mvmc_scaled_complex_sum_ordered(values + 2, 1, &sum) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            sum.state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO &&
            sum.log_abs_error_bound == -20.0,
        "all-numeric-zero sum retains its error bound");

  values[0] = finite_value(1.0, 1000.0);
  values[1] = finite_value(-1.0, 1000.0);
  CHECK(mvmc_scaled_complex_sum_ordered(values, 2, &sum) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            sum.state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO &&
            sum.max_input_log_abs == 1000.0 &&
            sum.cancellation_ratio == 0.0 &&
            isfinite(sum.log_abs_error_bound),
        "finite cancellation is numeric rather than exact zero");
  values[0] = finite_value(1.0, 1000.0);
  values[1] = finite_value(I, -1000.0);
  CHECK(mvmc_scaled_complex_sum_ordered(values, 2, &sum) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            sum.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
            close_double(sum.log_abs, 1000.0, 1.0e-15),
        "ordered sum retains dominant exponent without overflow");
  CHECK(mvmc_scaled_complex_sum_ordered(values + 4, 1, &sum) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            sum.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
        "all structural-zero sum remains exact");

  memset(&sentinel, 0x5a, sizeof(sentinel));
  values[5] = sentinel;
  CHECK(mvmc_scaled_complex_make_finite(1.0, INFINITY, -INFINITY,
                                        &sentinel) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            memcmp(&sentinel, &values[5], sizeof(sentinel)) == 0,
        "invalid pure operation preserves caller result");
  values[5] = finite_value(1.0, 0.0);
  values[5].phase = 2.0;
  memset(&sentinel, 0x5a, sizeof(sentinel));
  values[4] = sentinel;
  CHECK(mvmc_scaled_complex_multiply(&values[5], &values[0], &sentinel) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            memcmp(&sentinel, &values[4], sizeof(sentinel)) == 0,
        "non-unit input phase is rejected transactionally");
}

static void test_real_scaled_pfaffian(void) {
  double matrix[16] = {0.0};
  MVMCAbsolutePfaffianScaledValueResult result;
  MVMCAbsolutePfaffianScaledValueResult sentinel;
  MVMCAbsolutePfaffianValueResult raw_result;
  double small = exp(-500.0);
  double large = exp(500.0);
  double complex exported = 0.0;

  set_real_pair(matrix, 4, 0, 1, 2.0);
  set_real_pair(matrix, 4, 2, 3, -3.0);
  CHECK(mvmc_absolute_pfaffian_real_scaled_value(matrix, 4, 4, 0.0,
                                                 &result) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_absolute_pfaffian_real_value(matrix, 4, 4, 0.0,
                                              &raw_result) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_scaled_complex_export_common_scale(&result.value, 0.0,
                                                     &exported) ==
            MVMC_SCALED_EXPORT_OK &&
            close_complex(exported, raw_result.pfaffian, 1.0e-14),
        "representable real scaled Pfaffian matches P1 value-only API");

  memset(matrix, 0, sizeof(matrix));
  set_real_pair(matrix, 4, 0, 1, small);
  set_real_pair(matrix, 4, 2, 3, small);
  CHECK(mvmc_absolute_pfaffian_real_scaled_value(matrix, 4, 4, 0.0,
                                                 &result) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            result.factor_state == MVMC_PFAFFIAN_VALUE_WELL_PIVOTED &&
            result.value.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
            close_double(result.value.log_abs, -1000.0, 1.0e-13) &&
            close_complex(result.value.phase, 1.0, 1.0e-14),
        "real Pfaffian underflow-scale product remains finite in log form");
  CHECK(mvmc_scaled_complex_export_common_scale(&result.value, 0.0,
                                                 &exported) ==
            MVMC_SCALED_EXPORT_UNDERFLOW,
        "real underflow-scale raw export is classified");

  memset(matrix, 0, sizeof(matrix));
  set_real_pair(matrix, 4, 0, 1, large);
  set_real_pair(matrix, 4, 2, 3, -large);
  CHECK(mvmc_absolute_pfaffian_real_scaled_value(matrix, 4, 4, 0.0,
                                                 &result) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            result.value.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
            close_double(result.value.log_abs, 1000.0, 1.0e-13) &&
            close_complex(result.value.phase, -1.0, 1.0e-14),
        "real Pfaffian overflow-scale product remains finite in log form");
  CHECK(mvmc_scaled_complex_export_common_scale(&result.value, 1000.0,
                                                 &exported) ==
            MVMC_SCALED_EXPORT_OK &&
            close_complex(exported, -1.0, 1.0e-13),
        "overflow-scale Pfaffian exports at a common scale");

  memset(matrix, 0, sizeof(matrix));
  set_real_pair(matrix, 4, 0, 1, 1.0);
  set_real_pair(matrix, 4, 2, 3, 1.0e-20);
  CHECK(mvmc_absolute_pfaffian_real_scaled_value(matrix, 4, 4, 0.0,
                                                 &result) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            result.factor_state == MVMC_PFAFFIAN_VALUE_NEAR_PIVOT &&
            result.value.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO,
        "near-singular real Pfaffian retains a finite scaled value");
  memset(matrix, 0, sizeof(matrix));
  set_real_pair(matrix, 4, 0, 1, 1.0);
  CHECK(mvmc_absolute_pfaffian_real_scaled_value(matrix, 4, 4, 0.0,
                                                 &result) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            result.factor_state == MVMC_PFAFFIAN_VALUE_SINGULAR &&
            result.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
        "factorization-proven singular Pfaffian is exact zero");

  memset(&sentinel, 0x37, sizeof(sentinel));
  result = sentinel;
  CHECK(mvmc_absolute_pfaffian_real_scaled_value(NULL, 4, 4, 0.0,
                                                 &result) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            memcmp(&result, &sentinel, sizeof(result)) == 0,
        "invalid scaled Pfaffian call preserves old state");
}

static void test_complex_dense_oracle(void) {
  const int n = 6;
  double complex matrix[36] = {0.0};
  double complex expected;
  double expected_abs;
  MVMCAbsolutePfaffianScaledValueResult result;
  MVMCAbsolutePfaffianValueResult raw_result;
  double complex exported = 0.0;
  int left, right;

  for (right = 1; right < n; ++right) {
    for (left = 0; left < right; ++left) {
      double real_part = (double)(2 * left - right + 1) / 7.0;
      double imag_part = (double)(left + 3 * right + 1) / 19.0;
      set_complex_pair(matrix, n, left, right,
                       real_part + I * imag_part);
    }
  }
  expected = dense_pfaffian_recursive(matrix, n);
  expected_abs = cabs(expected);
  CHECK(expected_abs > 0.0 && isfinite(expected_abs),
        "independent dense oracle is finite");
  CHECK(mvmc_absolute_pfaffian_complex_scaled_value(matrix, n, n, 0.0,
                                                    &result) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            result.value.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
            close_double(result.value.log_abs, log(expected_abs), 2.0e-13) &&
            close_complex(result.value.phase, expected / expected_abs,
                          2.0e-13),
        "complex scaled Pfaffian matches independent dense oracle");
  CHECK(mvmc_absolute_pfaffian_complex_value(matrix, n, n, 0.0,
                                             &raw_result) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_scaled_complex_export_common_scale(&result.value, 0.0,
                                                     &exported) ==
            MVMC_SCALED_EXPORT_OK &&
            close_complex(exported, raw_result.pfaffian, 2.0e-13),
        "representable complex scaled Pfaffian matches P1 value-only API");

  matrix[0] = INFINITY;
  CHECK(mvmc_absolute_pfaffian_complex_scaled_value(matrix, n, n, 0.0,
                                                    &result) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            result.factor_state == MVMC_PFAFFIAN_VALUE_NONFINITE &&
            result.value.state == MVMC_SCALED_COMPLEX_NONFINITE,
        "complex nonfinite matrix is classified fail-loud");
}

static void test_scaled_projection(void) {
  MVMCAbsolutePfaffianScaledValueResult components[4];
  MVMCProjectedScaledAmplitudeResult aggregate;
  MVMCProjectedScaledAmplitudeResult empty;
  MVMCProjectedScaledAmplitudeResult sentinel;
  MVMCScaledComplex value;
  double complex weights[4] = {1.0, 2.0, -I, 0.5 + 0.25 * I};
  double complex exported = 0.0;
  const double complex expected = 4.625 + 1.5 * I;

  CHECK(mvmc_scaled_complex_from_raw(-2.5 + 0.75 * I, &value) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "NQP=1 component construction");
  components[0] =
      scaled_component(MVMC_PFAFFIAN_VALUE_WELL_PIVOTED, value);
  CHECK(mvmc_projected_scaled_amplitude_values(
            components, weights, 1, &aggregate) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            aggregate.valid && aggregate.well_pivoted_count == 1 &&
            aggregate.total.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO,
        "NQP=1 scaled projection remains finite");
  CHECK(mvmc_scaled_complex_export_common_scale(&aggregate.total, 0.0,
                                                 &exported) ==
            MVMC_SCALED_EXPORT_OK &&
            close_complex(exported, -2.5 + 0.75 * I, 2.0e-13),
        "NQP=1 scaled projection matches direct oracle");

  components[0] = scaled_component(MVMC_PFAFFIAN_VALUE_WELL_PIVOTED,
                                   finite_value(1.0, 0.0));
  weights[0] = DBL_MAX + I * DBL_MAX;
  CHECK(mvmc_projected_scaled_amplitude_values(
            components, weights, 1, &aggregate) ==
            MVMC_PFAFFIAN_STATUS_OK && aggregate.valid &&
            aggregate.total.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
            aggregate.total.log_abs > log(DBL_MAX) &&
            close_complex(aggregate.total.phase,
                          (1.0 + I) / sqrt(2.0), 2.0e-13),
        "finite complex weight survives raw-magnitude overflow");
  weights[0] = 1.0;

  CHECK(mvmc_scaled_complex_from_raw(2.0, &value) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "projection regular component construction");
  components[0] =
      scaled_component(MVMC_PFAFFIAN_VALUE_WELL_PIVOTED, value);
  CHECK(mvmc_scaled_complex_make_exact_zero(&value) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "projection singular component construction");
  components[1] = scaled_component(MVMC_PFAFFIAN_VALUE_SINGULAR, value);
  CHECK(mvmc_scaled_complex_from_raw(-1.0 + I, &value) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "projection near component construction");
  components[2] = scaled_component(MVMC_PFAFFIAN_VALUE_NEAR_PIVOT, value);
  CHECK(mvmc_scaled_complex_from_raw(3.0 - 0.5 * I, &value) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "projection final component construction");
  components[3] =
      scaled_component(MVMC_PFAFFIAN_VALUE_WELL_PIVOTED, value);

  CHECK(mvmc_projected_scaled_amplitude_values(
            components, weights, 4, &aggregate) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            aggregate.valid && aggregate.well_pivoted_count == 2 &&
            aggregate.near_pivot_count == 1 &&
            aggregate.exact_zero_count == 1 &&
            aggregate.numeric_zero_count == 0 &&
            aggregate.total.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO,
        "global-QP scaled projection preserves classifications");
  CHECK(mvmc_scaled_complex_export_common_scale(&aggregate.total, 0.0,
                                                 &exported) ==
            MVMC_SCALED_EXPORT_OK &&
            close_complex(exported, expected, 2.0e-13),
        "global-QP scaled projection matches direct oracle");
  CHECK(close_double(aggregate.log_sum_abs,
                     log(2.0 + sqrt(2.0) + cabs(1.625 + 0.5 * I)),
                     2.0e-13),
        "projection log absolute sum matches oracle");

  CHECK(mvmc_projected_scaled_amplitude_value_slice(
            NULL, 0, weights, 4, 2, 2, &empty) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            empty.valid &&
            empty.total.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
        "empty local QP slice is a valid exact-zero partial sum");

  memset(&sentinel, 0x2c, sizeof(sentinel));
  aggregate = sentinel;
  components[1].factor_state = MVMC_PFAFFIAN_VALUE_SINGULAR;
  components[1].value = finite_value(1.0, 0.0);
  CHECK(mvmc_projected_scaled_amplitude_values(
            components, weights, 4, &aggregate) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            memcmp(&aggregate, &sentinel, sizeof(aggregate)) == 0,
        "invalid projected component preserves old result");
  components[1] = scaled_component(MVMC_PFAFFIAN_VALUE_WELL_PIVOTED,
                                   value);
  CHECK(mvmc_scaled_complex_make_exact_zero(&components[1].value) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "inconsistent projected component fixture construction");
  CHECK(mvmc_projected_scaled_amplitude_values(
            components, weights, 4, &aggregate) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "well-pivoted component cannot carry exact-zero Pfaffian state");
  aggregate = sentinel;
  CHECK(mvmc_projected_scaled_amplitude_values(
            components, weights, SIZE_MAX / sizeof(MVMCScaledComplex) + 1,
            &aggregate) == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            memcmp(&aggregate, &sentinel, sizeof(aggregate)) == 0,
        "projected component allocation extent overflow is atomic");
}

static void test_rank_invariance(void) {
#ifdef _mpi_use
  int rank = 0;
  int size = 1;
  int qp_start;
  int qp_end;
  MVMCAbsolutePfaffianScaledValueResult components[4];
  MVMCProjectedScaledAmplitudeResult aggregate;
  MVMCProjectedScaledAmplitudeResult partial;
  MVMCScaledComplex value;
  double complex weights[4] = {1.0, -I, 0.5, 2.0 + I};
  double fields[4];
  double minimum[4];
  double maximum[4];
  int state;
  int minimum_state;
  int maximum_state;
  int index;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  for (index = 0; index < 4; ++index) {
    CHECK(mvmc_scaled_complex_make_finite(
              index % 2 == 0 ? 1.0 : I, (double)(index - 2) * 400.0,
              -INFINITY, &value) == MVMC_PFAFFIAN_STATUS_OK,
          "rank scaled component construction");
    components[index] = scaled_component(
        index == 3 ? MVMC_PFAFFIAN_VALUE_NEAR_PIVOT
                   : MVMC_PFAFFIAN_VALUE_WELL_PIVOTED,
        value);
  }
  CHECK(mvmc_projected_scaled_amplitude_values(
            components, weights, 4, &aggregate) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "rank-invariant full scaled projection");
  fields[0] = aggregate.total.log_abs;
  fields[1] = creal(aggregate.total.phase);
  fields[2] = cimag(aggregate.total.phase);
  fields[3] = aggregate.total.log_abs_error_bound;
  state = (int)aggregate.total.state;
  MPI_Allreduce(fields, minimum, 4, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(fields, maximum, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&state, &minimum_state, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&state, &maximum_state, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  CHECK(memcmp(minimum, maximum, sizeof(minimum)) == 0 &&
            minimum_state == maximum_state,
        "phase/log/status are bitwise invariant at MPI rank 2/4");

  qp_start = 4 * rank / size;
  qp_end = 4 * (rank + 1) / size;
  CHECK(mvmc_projected_scaled_amplitude_value_slice(
            components + qp_start, (size_t)(qp_end - qp_start), weights, 4,
            qp_start, qp_end, &partial) == MVMC_PFAFFIAN_STATUS_OK &&
            partial.valid,
        "rank-local scaled QP slice including empty ownership is valid");
#endif
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  test_scaled_pure_math();
  test_real_scaled_pfaffian();
  test_complex_dense_oracle();
  test_scaled_projection();
  test_rank_invariance();

  if (failures != 0) {
    fprintf(stderr, "ScaledPfaffian_Unit: %d failure(s)\n", failures);
#ifdef _mpi_use
    MPI_Finalize();
#endif
    return 1;
  }
  puts("ScaledPfaffian_Unit: PASS");
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return 0;
}
