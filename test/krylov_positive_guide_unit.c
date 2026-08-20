#include "krylov_positive_guide.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

static int failures = 0;

#define CHECK(condition, message)                                             \
  do {                                                                        \
    if (!(condition)) {                                                       \
      fprintf(stderr, "KrylovPositiveGuide_Unit FAIL: %s (line %d)\n",       \
              (message), __LINE__);                                           \
      ++failures;                                                             \
    }                                                                         \
  } while (0)

static int close_double(double actual, double expected, double tolerance) {
  return fabs(actual - expected) <= tolerance * fmax(1.0, fabs(expected));
}

static MVMCScaledComplex finite_raw(double complex raw) {
  MVMCScaledComplex value;
  memset(&value, 0, sizeof(value));
  CHECK(mvmc_scaled_complex_from_raw(raw, &value) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            value.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO,
        "finite raw scaled fixture");
  return value;
}

static MVMCScaledComplex finite_log(double log_abs) {
  MVMCScaledComplex value;
  memset(&value, 0, sizeof(value));
  CHECK(mvmc_scaled_complex_make_finite(1.0, log_abs, -INFINITY, &value) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            value.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO,
        "finite log scaled fixture");
  return value;
}

static MVMCKrylovPositiveGuidePolicy default_policy(int order) {
  MVMCKrylovPositiveGuidePolicy policy;
  int index;
  memset(&policy, 0, sizeof(policy));
  policy.order = order;
  policy.eta = 0.25;
  policy.policy_hash = UINT64_C(0x5034440000000001);
  for (index = 0; index <= MVMC_KRYLOV_MAX_ORDER; ++index) {
    policy.lambda[index] = 1.0;
  }
  return policy;
}

static void test_raw_oracle(void) {
  MVMCKrylovPositiveGuidePolicy policy = default_policy(3);
  MVMCScaledComplex values[MVMC_KRYLOV_MAX_ORDER + 1];
  MVMCKrylovPositiveGuideResult result;
  const double expected = log(0.25 + 4.0 + 2.0 * 9.0 + 0.5 * 16.0);

  policy.lambda[1] = 2.0;
  policy.lambda[2] = 7.0;
  policy.lambda[3] = 0.5;
  values[0] = finite_raw(2.0);
  values[1] = finite_raw(0.0 + 3.0 * I);
  CHECK(mvmc_scaled_complex_make_exact_zero(&values[2]) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "exact zero scaled fixture");
  values[3] = finite_raw(-4.0);

  CHECK(mvmc_krylov_positive_guide_evaluate(&policy, values, 4, &result) ==
            MVMC_KRYLOV_STATUS_OK &&
            result.valid,
        "raw oracle guide evaluation succeeds");
  CHECK(close_double(result.log_guide, expected, 1.0e-14),
        "guide log matches raw oracle");
  CHECK(result.finite_component_count == 3 &&
            result.exact_zero_component_count == 1 &&
            result.numeric_zero_component_count == 0,
        "guide records component state counts");
  CHECK(close_double(result.log_floor, log(0.25), 1.0e-15),
        "guide records floor contribution");
}

static void test_log_domain_range(void) {
  MVMCKrylovPositiveGuidePolicy policy = default_policy(1);
  MVMCScaledComplex values[2];
  MVMCKrylovPositiveGuideResult result;
  const double expected = 1400.0 + log1p(2.0 * exp(-2.0));

  policy.eta = 0.125;
  policy.lambda[0] = 1.0;
  policy.lambda[1] = 2.0;
  values[0] = finite_log(700.0);
  values[1] = finite_log(699.0);

  CHECK(mvmc_krylov_positive_guide_evaluate(&policy, values, 2, &result) ==
            MVMC_KRYLOV_STATUS_OK &&
            result.valid,
        "large log-domain guide evaluation succeeds");
  CHECK(close_double(result.log_guide, expected, 1.0e-14),
        "guide log-sum-exp avoids raw overflow");
}

static void test_basis_scale_oracle(void) {
  MVMCKrylovPositiveGuidePolicy policy = default_policy(3);
  MVMCScaledComplex values[MVMC_KRYLOV_MAX_ORDER + 1];
  MVMCKrylovPositiveGuideResult result;
  double expected;

  policy.eta = 0.125;
  policy.lambda[0] = 1.0;
  policy.lambda[1] = 0.5;
  policy.lambda[2] = 2.0;
  policy.lambda[3] = 1.0;
  policy.log_basis_scale[0] = log(3.0);
  policy.log_basis_scale[1] = log(0.25);
  policy.log_basis_scale[2] = -2.0;
  policy.log_basis_scale[3] = 0.0;
  values[0] = finite_raw(2.0);
  values[1] = finite_raw(0.0 + 4.0 * I);
  CHECK(mvmc_scaled_complex_make_exact_zero(&values[2]) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "basis-scale exact zero fixture");
  CHECK(mvmc_scaled_complex_make_exact_zero(&values[3]) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "basis-scale second exact zero fixture");
  expected = log(0.125 +
                 4.0 * exp(2.0 * policy.log_basis_scale[0]) +
                 0.5 * 16.0 *
                     exp(2.0 * policy.log_basis_scale[1]));

  CHECK(mvmc_krylov_positive_guide_evaluate(&policy, values, 4, &result) ==
            MVMC_KRYLOV_STATUS_OK &&
            result.valid,
        "basis-scaled guide evaluation succeeds");
  CHECK(close_double(result.log_guide, expected, 1.0e-14),
        "guide basis scale matches oracle");
  CHECK(result.finite_component_count == 2 &&
            result.exact_zero_component_count == 2,
        "basis scale preserves component state counts");
}

static void test_zero_floor_support(void) {
  MVMCKrylovPositiveGuidePolicy policy = default_policy(2);
  MVMCScaledComplex values[3];
  MVMCKrylovPositiveGuideResult result;

  policy.eta = 0.5;
  CHECK(mvmc_scaled_complex_make_exact_zero(&values[0]) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "exact zero fixture for floor test");
  CHECK(mvmc_scaled_complex_make_numeric_zero(-20.0, 0.0, -INFINITY, 0.0,
                                               &values[1]) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "numeric zero fixture for floor test");
  CHECK(mvmc_scaled_complex_make_exact_zero(&values[2]) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "second exact zero fixture for floor test");

  CHECK(mvmc_krylov_positive_guide_evaluate(&policy, values, 3, &result) ==
            MVMC_KRYLOV_STATUS_OK &&
            result.valid,
        "floor keeps guide positive when all components are zero-like");
  CHECK(close_double(result.log_guide, log(0.5), 1.0e-15),
        "floor-only guide value is eta");
  CHECK(isinf(result.log_signal_sum) && result.log_signal_sum < 0.0,
        "zero-like components do not create exact signal contribution");
  CHECK(result.exact_zero_component_count == 2 &&
            result.numeric_zero_component_count == 1,
        "floor guide records zero component counts");
}

static void test_rejections(void) {
  MVMCKrylovPositiveGuidePolicy policy = default_policy(1);
  MVMCScaledComplex values[2];
  MVMCKrylovPositiveGuideResult result;

  values[0] = finite_raw(1.0);
  values[1] = finite_raw(2.0);
  policy.lambda[1] = -1.0;
  CHECK(mvmc_krylov_positive_guide_evaluate(&policy, values, 2, &result) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !result.valid,
        "negative lambda is rejected");

  policy = default_policy(1);
  policy.eta = 0.0;
  CHECK(mvmc_krylov_positive_guide_evaluate(&policy, values, 2, &result) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !result.valid,
        "nonpositive eta is rejected");

  policy = default_policy(1);
  values[1].state = MVMC_SCALED_COMPLEX_NONFINITE;
  CHECK(mvmc_krylov_positive_guide_evaluate(&policy, values, 2, &result) ==
            MVMC_KRYLOV_STATUS_NONFINITE &&
            !result.valid,
        "nonfinite guide component is rejected");

  policy = default_policy(1);
  values[1] = finite_raw(2.0);
  policy.log_basis_scale[1] = INFINITY;
  CHECK(mvmc_krylov_positive_guide_evaluate(&policy, values, 2, &result) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !result.valid,
        "nonfinite basis scale is rejected");
}

static void test_acceptance(void) {
  MVMCKrylovPositiveGuidePolicy policy = default_policy(0);
  MVMCScaledComplex current_value[1];
  MVMCScaledComplex proposal_value[1];
  MVMCKrylovPositiveGuideResult current;
  MVMCKrylovPositiveGuideResult proposal;
  MVMCKrylovPositiveGuideAcceptance acceptance;
  const double expected_log_ratio = log(0.25 + 4.0) - log(0.25 + 1.0);

  current_value[0] = finite_raw(1.0);
  proposal_value[0] = finite_raw(2.0);
  CHECK(mvmc_krylov_positive_guide_evaluate(
            &policy, current_value, 1, &current) == MVMC_KRYLOV_STATUS_OK,
        "current guide for acceptance");
  CHECK(mvmc_krylov_positive_guide_evaluate(
            &policy, proposal_value, 1, &proposal) == MVMC_KRYLOV_STATUS_OK,
        "proposal guide for acceptance");
  CHECK(mvmc_krylov_positive_guide_acceptance(
            &current, &proposal, 0.0, 0.999, &acceptance) ==
            MVMC_KRYLOV_STATUS_OK &&
            acceptance.valid && acceptance.accepted &&
            close_double(acceptance.log_acceptance_ratio, expected_log_ratio,
                         1.0e-14),
        "positive guide move accepts without exponentiating the ratio");
  CHECK(mvmc_krylov_positive_guide_acceptance(
            &proposal, &current, 0.0, 0.10, &acceptance) ==
            MVMC_KRYLOV_STATUS_OK &&
            acceptance.valid && acceptance.accepted,
        "downhill guide move accepts below exp(log-ratio)");
  CHECK(mvmc_krylov_positive_guide_acceptance(
            &proposal, &current, 0.0, 0.90, &acceptance) ==
            MVMC_KRYLOV_STATUS_OK &&
            acceptance.valid && !acceptance.accepted,
        "downhill guide move rejects above exp(log-ratio)");
  CHECK(mvmc_krylov_positive_guide_acceptance(
            &proposal, &current, 0.0, 1.0, &acceptance) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !acceptance.valid,
        "uniform draw upper endpoint is rejected");
}

int main(void) {
  test_raw_oracle();
  test_log_domain_range();
  test_basis_scale_oracle();
  test_zero_floor_support();
  test_rejections();
  test_acceptance();
  if (failures == 0) {
    printf("krylov positive guide unit: PASS\n");
  }
  return failures == 0 ? 0 : 1;
}
