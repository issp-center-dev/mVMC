#include "power_lanczos_numeric_evidence.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(condition, ...)                                                \
  do {                                                                       \
    if (!(condition)) {                                                      \
      fprintf(stderr, "PowerLanczosNumericEvidence_Unit FAIL: ");          \
      fprintf(stderr, __VA_ARGS__);                                          \
      fprintf(stderr, " (line %d)\n", __LINE__);                          \
      ++failures;                                                            \
    }                                                                        \
  } while (0)

static MVMCScaledComplex finite_value(double complex value,
                                      double absolute_error) {
  MVMCScaledComplex scaled;
  const double magnitude = cabs(value);
  memset(&scaled, 0, sizeof(scaled));
  CHECK(magnitude > 0.0 && isfinite(magnitude) &&
            absolute_error >= 0.0 && absolute_error < magnitude,
        "finite fixture precondition");
  CHECK(mvmc_scaled_complex_make_finite(
            value / magnitude, log(magnitude),
            absolute_error == 0.0 ? -INFINITY : log(absolute_error),
            &scaled) == MVMC_PFAFFIAN_STATUS_OK &&
            scaled.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO,
        "finite fixture construction");
  return scaled;
}

static void test_division_bound_covers_corners(void) {
  const double numerator_value = 4.0;
  const double numerator_error = 1.0e-12;
  const double denominator_value = 2.0;
  const double denominator_error = 2.0e-13;
  MVMCScaledComplex numerator =
      finite_value(numerator_value, numerator_error);
  MVMCScaledComplex denominator =
      finite_value(denominator_value, denominator_error);
  MVMCScaledComplex quotient;
  MVMCPowerLanczosNumericEvidence evidence;
  int numerator_sign;
  int denominator_sign;
  memset(&quotient, 0, sizeof(quotient));
  memset(&evidence, 0, sizeof(evidence));

  CHECK(mvmc_power_lanczos_scaled_divide(
            &numerator, &denominator, &quotient) ==
                MVMC_KRYLOV_STATUS_OK &&
            quotient.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO,
        "finite scaled division");
  CHECK(mvmc_power_lanczos_scaled_export_evidence(
            &quotient, 0.0, &evidence) == MVMC_KRYLOV_STATUS_OK &&
            evidence.valid && evidence.status == MVMC_KRYLOV_STATUS_OK &&
            evidence.version == MVMC_POWER_LANCZOS_NUMERIC_EVIDENCE_VERSION &&
            evidence.support_flags ==
                MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO &&
            isfinite(evidence.absolute_numeric_bound) &&
            evidence.absolute_numeric_bound > 0.0,
        "finite quotient export evidence");
  for (numerator_sign = -1; numerator_sign <= 1; numerator_sign += 2) {
    for (denominator_sign = -1; denominator_sign <= 1;
         denominator_sign += 2) {
      const double exact_corner =
          (numerator_value + numerator_sign * numerator_error) /
          (denominator_value + denominator_sign * denominator_error);
      CHECK(fabs(exact_corner - creal(evidence.value)) <=
                evidence.absolute_numeric_bound,
            "quotient bound covers corner (%d, %d): delta=%.17g bound=%.17g",
            numerator_sign, denominator_sign,
            fabs(exact_corner - creal(evidence.value)),
            evidence.absolute_numeric_bound);
    }
  }
}

static void test_zero_and_range_classification(void) {
  MVMCScaledComplex numerator;
  MVMCScaledComplex denominator = finite_value(2.0, 0.0);
  MVMCScaledComplex quotient;
  MVMCPowerLanczosNumericEvidence evidence;

  CHECK(mvmc_scaled_complex_make_exact_zero(&numerator) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_power_lanczos_scaled_divide(
                &numerator, &denominator, &quotient) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_scaled_export_evidence(
                &quotient, 0.0, &evidence) == MVMC_KRYLOV_STATUS_OK &&
            evidence.valid && evidence.value == 0.0 &&
            evidence.absolute_numeric_bound == 0.0 &&
            evidence.support_flags ==
                MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO,
        "exact-zero quotient remains structural exact zero");

  CHECK(mvmc_scaled_complex_make_numeric_zero(
            log(5.0e-11), 0.0, -INFINITY, 0.0, &numerator) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_power_lanczos_scaled_divide(
                &numerator, &denominator, &quotient) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_scaled_export_evidence(
                &quotient, -log(2.0), &evidence) ==
                MVMC_KRYLOV_STATUS_OK &&
            evidence.valid && evidence.value == 0.0 &&
            evidence.absolute_numeric_bound >= 5.0e-11 &&
            evidence.support_flags ==
                MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO,
        "numeric-zero quotient retains a positive propagated bound");

  numerator = finite_value(1.0, 0.0);
  numerator.log_abs = -1000.0;
  numerator.max_input_log_abs = -1000.0;
  numerator.cancellation_log_abs = -1000.0;
  CHECK(mvmc_power_lanczos_scaled_export_evidence(
            &numerator, 0.0, &evidence) == MVMC_KRYLOV_STATUS_OK &&
            evidence.valid && evidence.value == 0.0 &&
            evidence.absolute_numeric_bound > 0.0 &&
            evidence.support_flags ==
                MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO,
        "finite underflow is numeric zero rather than exact zero");

  numerator.log_abs = 1000.0;
  numerator.max_input_log_abs = 1000.0;
  numerator.cancellation_log_abs = 1000.0;
  CHECK(mvmc_power_lanczos_scaled_export_evidence(
            &numerator, 0.0, &evidence) ==
                MVMC_KRYLOV_STATUS_NONFINITE &&
            !evidence.valid && evidence.status == MVMC_KRYLOV_STATUS_NONFINITE,
        "overflow export fails closed");
}

static void test_denominator_and_recenter_boundaries(void) {
  MVMCScaledComplex numerator = finite_value(1.0, 0.0);
  MVMCScaledComplex denominator;
  MVMCScaledComplex quotient;
  MVMCPowerLanczosNumericEvidence evidence;
  MVMCPowerLanczosNumericEvidence recentered;
  const double production_value = 2.0 + 4.0 * DBL_EPSILON;
  double original_bound;

  CHECK(mvmc_scaled_complex_make_numeric_zero(
            log(1.0e-12), 0.0, -INFINITY, 0.0, &denominator) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_power_lanczos_scaled_divide(
                &numerator, &denominator, &quotient) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "numeric-zero denominator is rejected");

  quotient = finite_value(2.0, 1.0e-13);
  CHECK(mvmc_power_lanczos_scaled_export_evidence(
            &quotient, 0.0, &evidence) == MVMC_KRYLOV_STATUS_OK,
        "recenter fixture export");
  original_bound = evidence.absolute_numeric_bound;
  CHECK(mvmc_power_lanczos_numeric_evidence_recenter(
            &evidence, production_value, &recentered) ==
                MVMC_KRYLOV_STATUS_OK &&
            recentered.valid && recentered.value == production_value &&
            recentered.absolute_numeric_bound >=
                original_bound + fabs(production_value - creal(evidence.value)) &&
            fabs(2.0 - creal(recentered.value)) <=
                recentered.absolute_numeric_bound,
        "recenter widens the bound around production bytes");

  CHECK(mvmc_scaled_complex_make_exact_zero(&quotient) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_power_lanczos_scaled_export_evidence(
                &quotient, 0.0, &evidence) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_numeric_evidence_recenter(
                &evidence, DBL_TRUE_MIN, &recentered) ==
                MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE &&
            !recentered.valid,
        "exact zero cannot be recentered to a nonzero production value");
}

static void test_weighted_norm_and_guide_product(void) {
  const double central = 2.0;
  const double input_error = 1.0e-12;
  const double log_scale[1] = {0.0};
  const double guide_weight[1] = {0.5};
  MVMCScaledComplex value = finite_value(central, input_error);
  MVMCScaledComplex guide;
  MVMCPowerLanczosNumericEvidence evidence;
  int sign;

  CHECK(mvmc_power_lanczos_scaled_weighted_norm(
            &value, log_scale, guide_weight, 1, 1.0, &guide) ==
                MVMC_KRYLOV_STATUS_OK &&
            guide.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
            isfinite(guide.log_abs_error_bound),
        "weighted norm retains input and arithmetic error");
  CHECK(mvmc_power_lanczos_guide_normalized_product_evidence(
            &value, 0.0, &value, 0.0, 1, &guide, 4.0 / 3.0,
            &evidence) == MVMC_KRYLOV_STATUS_OK &&
            evidence.valid && evidence.value == 4.0 / 3.0 &&
            evidence.absolute_numeric_bound > 0.0,
        "guide-normalized product evidence");
  for (sign = -1; sign <= 1; sign += 2) {
    const double perturbed = central + sign * input_error;
    const double exact = perturbed * perturbed /
                         (1.0 + 0.5 * perturbed * perturbed);
    CHECK(fabs(exact - creal(evidence.value)) <=
              evidence.absolute_numeric_bound,
          "guide product bound covers endpoint %d: delta=%.17g bound=%.17g",
          sign, fabs(exact - creal(evidence.value)),
          evidence.absolute_numeric_bound);
  }

  CHECK(mvmc_power_lanczos_guide_normalized_product_evidence(
            &value, 0.0, &value, 0.0, -1, &guide, -4.0 / 3.0,
            &evidence) == MVMC_KRYLOV_STATUS_OK &&
            evidence.valid && evidence.value == -4.0 / 3.0,
        "negative product sign is retained");
}

int main(void) {
  test_division_bound_covers_corners();
  test_zero_and_range_classification();
  test_denominator_and_recenter_boundaries();
  test_weighted_norm_and_guide_product();
  if (failures != 0) {
    fprintf(stderr, "PowerLanczosNumericEvidence_Unit: %d failure(s)\n",
            failures);
    return 1;
  }
  printf("PowerLanczosNumericEvidence_Unit: PASS\n");
  return 0;
}
