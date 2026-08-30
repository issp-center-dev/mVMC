#include "power_lanczos_numeric_evidence.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "power_lanczos_numeric_evidence.c requires bounded Krylov"
#endif

#include <float.h>
#include <math.h>
#include <string.h>

static int finite_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static double log_add_bound(double left, double right) {
  if (isinf(left) && left < 0.0) return right;
  if (isinf(right) && right < 0.0) return left;
  if (!isfinite(left) || !isfinite(right)) return NAN;
  if (left < right) {
    const double temporary = left;
    left = right;
    right = temporary;
  }
  return left + log1p(exp(right - left));
}

static int finite_magnitude(double complex value, double *magnitude) {
  const double real_part = fabs(creal(value));
  const double imaginary_part = fabs(cimag(value));
  const double scale = fmax(real_part, imaginary_part);
  if (magnitude == NULL || !isfinite(real_part) ||
      !isfinite(imaginary_part)) {
    return 0;
  }
  if (scale == 0.0) {
    *magnitude = 0.0;
    return 1;
  }
  *magnitude = scale * hypot(real_part / scale, imaginary_part / scale);
  return isfinite(*magnitude);
}

static int exp_upper(double log_value, double *value) {
  const double true_minimum = nextafter(0.0, 1.0);
  const double minimum_log = log(true_minimum);
  const double maximum_log = log(DBL_MAX);
  double candidate;
  if (value == NULL || !isfinite(log_value) || log_value > maximum_log) {
    return 0;
  }
  if (log_value <= minimum_log) {
    *value = true_minimum;
    return 1;
  }
  candidate = exp(log_value);
  if (!isfinite(candidate) || candidate <= 0.0) return 0;
  candidate = nextafter(candidate, INFINITY);
  if (!isfinite(candidate)) return 0;
  *value = candidate;
  return 1;
}

static int add_upper(double left, double right, double *sum) {
  double candidate;
  if (sum == NULL || !isfinite(left) || left < 0.0 ||
      !isfinite(right) || right < 0.0) {
    return 0;
  }
  candidate = left + right;
  if (!isfinite(candidate)) return 0;
  if (candidate != 0.0) candidate = nextafter(candidate, INFINITY);
  if (!isfinite(candidate)) return 0;
  *sum = candidate;
  return 1;
}

static void invalidate_evidence(MVMCKrylovStatus status,
                                MVMCPowerLanczosNumericEvidence *evidence) {
  if (evidence == NULL) return;
  memset(evidence, 0, sizeof(*evidence));
  evidence->status = status;
  evidence->value = NAN + I * NAN;
  evidence->absolute_numeric_bound = NAN;
}

static MVMCKrylovStatus pfaffian_status_to_krylov(
    MVMCPfaffianStatus status, const MVMCScaledComplex *result) {
  if (status != MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (!mvmc_scaled_complex_is_valid(result)) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  return result->state == MVMC_SCALED_COMPLEX_NONFINITE
             ? MVMC_KRYLOV_STATUS_NONFINITE
             : MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus scaled_log_scale(
    const MVMCScaledComplex *value, double log_scale,
    MVMCScaledComplex *result) {
  MVMCScaledComplex factor;
  MVMCPfaffianStatus status;
  if (result == NULL || !mvmc_scaled_complex_is_valid(value) ||
      !isfinite(log_scale)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (log_scale == 0.0) {
    *result = *value;
    return value->state == MVMC_SCALED_COMPLEX_NONFINITE
               ? MVMC_KRYLOV_STATUS_NONFINITE
               : MVMC_KRYLOV_STATUS_OK;
  }
  status = mvmc_scaled_complex_make_finite(
      1.0, log_scale, -INFINITY, &factor);
  if (status == MVMC_PFAFFIAN_STATUS_OK) {
    status = mvmc_scaled_complex_multiply(value, &factor, result);
  }
  return pfaffian_status_to_krylov(status, result);
}

static MVMCKrylovStatus scaled_conjugate(
    const MVMCScaledComplex *value, MVMCScaledComplex *result) {
  if (result == NULL || !mvmc_scaled_complex_is_valid(value)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *result = *value;
  if (result->state == MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
    result->phase = conj(result->phase);
  }
  return result->state == MVMC_SCALED_COMPLEX_NONFINITE
             ? MVMC_KRYLOV_STATUS_NONFINITE
             : MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_scaled_divide(
    const MVMCScaledComplex *numerator,
    const MVMCScaledComplex *denominator,
    MVMCScaledComplex *quotient) {
  double denominator_relative_error = 0.0;
  double denominator_lower_log;
  double quotient_log;
  double quotient_error_log = -INFINITY;
  double numerator_error_term = -INFINITY;
  double denominator_error_term = -INFINITY;
  double complex phase;
  MVMCPfaffianStatus pfaffian_status;
  if (quotient == NULL || !mvmc_scaled_complex_is_valid(numerator) ||
      !mvmc_scaled_complex_is_valid(denominator) ||
      denominator->state != MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (numerator->state == MVMC_SCALED_COMPLEX_NONFINITE) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  if (isfinite(denominator->log_abs_error_bound)) {
    denominator_relative_error =
        exp(denominator->log_abs_error_bound - denominator->log_abs);
    denominator_relative_error =
        nextafter(denominator_relative_error, INFINITY);
    if (!isfinite(denominator_relative_error) ||
        denominator_relative_error < 0.0 ||
        denominator_relative_error >= 1.0) {
      return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    }
  }
  denominator_lower_log =
      denominator->log_abs + log1p(-denominator_relative_error);
  if (!isfinite(denominator_lower_log)) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  if (numerator->state == MVMC_SCALED_COMPLEX_EXACT_ZERO) {
    pfaffian_status = mvmc_scaled_complex_make_exact_zero(quotient);
    return pfaffian_status == MVMC_PFAFFIAN_STATUS_OK
               ? MVMC_KRYLOV_STATUS_OK
               : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (numerator->state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO) {
    const double error_log =
        numerator->log_abs_error_bound - denominator_lower_log;
    if (!isfinite(error_log)) return MVMC_KRYLOV_STATUS_NONFINITE;
    pfaffian_status = mvmc_scaled_complex_make_numeric_zero(
        error_log,
        isfinite(numerator->max_input_log_abs)
            ? numerator->max_input_log_abs - denominator_lower_log
            : -INFINITY,
        isfinite(numerator->cancellation_log_abs)
            ? numerator->cancellation_log_abs - denominator_lower_log
            : -INFINITY,
        numerator->cancellation_ratio, quotient);
    return pfaffian_status == MVMC_PFAFFIAN_STATUS_OK
               ? MVMC_KRYLOV_STATUS_OK
               : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }

  quotient_log = numerator->log_abs - denominator->log_abs;
  phase = numerator->phase * conj(denominator->phase);
  if (!isfinite(quotient_log) || !finite_complex(phase)) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  if (isfinite(numerator->log_abs_error_bound)) {
    numerator_error_term = numerator->log_abs_error_bound;
  }
  if (isfinite(denominator->log_abs_error_bound)) {
    denominator_error_term =
        quotient_log + denominator->log_abs_error_bound;
  }
  if (isfinite(numerator_error_term) ||
      isfinite(denominator_error_term)) {
    quotient_error_log =
        log_add_bound(numerator_error_term, denominator_error_term) -
        denominator_lower_log;
  }
  quotient_error_log = log_add_bound(
      quotient_error_log, quotient_log + log(64.0 * DBL_EPSILON));
  if (!isfinite(quotient_error_log)) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  pfaffian_status = mvmc_scaled_complex_make_finite(
      phase, quotient_log, quotient_error_log, quotient);
  return pfaffian_status == MVMC_PFAFFIAN_STATUS_OK
             ? MVMC_KRYLOV_STATUS_OK
             : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

MVMCKrylovStatus mvmc_power_lanczos_scaled_weighted_norm(
    const MVMCScaledComplex *values, const double *log_scale,
    const double *weight, size_t value_count, double eta,
    MVMCScaledComplex *result) {
  MVMCScaledComplex terms[MVMC_KRYLOV_MAX_ORDER + 2];
  size_t index;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  MVMCPfaffianStatus pfaffian_status;
  if (result == NULL || values == NULL || log_scale == NULL ||
      weight == NULL || value_count == 0 ||
      value_count > MVMC_KRYLOV_MAX_ORDER + 1 || !isfinite(eta) ||
      eta < 0.0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (index = 0; index < value_count; ++index) {
    MVMCScaledComplex scaled;
    MVMCScaledComplex conjugated;
    MVMCScaledComplex squared;
    MVMCScaledComplex factor;
    if (!mvmc_scaled_complex_is_valid(values + index) ||
        !isfinite(log_scale[index]) || !isfinite(weight[index]) ||
        weight[index] < 0.0) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    if (weight[index] == 0.0) {
      pfaffian_status = mvmc_scaled_complex_make_exact_zero(terms + index);
      status = pfaffian_status_to_krylov(pfaffian_status, terms + index);
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      continue;
    }
    status = scaled_log_scale(values + index, log_scale[index], &scaled);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = scaled_conjugate(&scaled, &conjugated);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      pfaffian_status = mvmc_scaled_complex_multiply(
          &conjugated, &scaled, &squared);
      status = pfaffian_status_to_krylov(pfaffian_status, &squared);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      pfaffian_status = mvmc_scaled_complex_make_finite(
          1.0, log(weight[index]), -INFINITY, &factor);
      status = pfaffian_status_to_krylov(pfaffian_status, &factor);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      pfaffian_status = mvmc_scaled_complex_multiply(
          &squared, &factor, terms + index);
      status = pfaffian_status_to_krylov(pfaffian_status, terms + index);
    }
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
  }
  if (eta == 0.0) {
    pfaffian_status =
        mvmc_scaled_complex_make_exact_zero(terms + value_count);
  } else {
    pfaffian_status = mvmc_scaled_complex_make_finite(
        1.0, log(eta), -INFINITY, terms + value_count);
  }
  status = pfaffian_status_to_krylov(
      pfaffian_status, terms + value_count);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  pfaffian_status = mvmc_scaled_complex_sum_ordered(
      terms, value_count + 1, result);
  return pfaffian_status_to_krylov(pfaffian_status, result);
}

MVMCKrylovStatus mvmc_power_lanczos_guide_normalized_product_evidence(
    const MVMCScaledComplex *left, double left_log_scale,
    const MVMCScaledComplex *right, double right_log_scale, int sign,
    const MVMCScaledComplex *guide, double complex production_value,
    MVMCPowerLanczosNumericEvidence *evidence) {
  MVMCScaledComplex scaled_left;
  MVMCScaledComplex scaled_right;
  MVMCScaledComplex conjugated_left;
  MVMCScaledComplex numerator;
  MVMCScaledComplex quotient;
  MVMCPowerLanczosNumericEvidence independent;
  MVMCKrylovStatus status;
  MVMCPfaffianStatus pfaffian_status;
  if (evidence == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_evidence(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, evidence);
  if ((sign != -1 && sign != 1) ||
      !isfinite(left_log_scale) || !isfinite(right_log_scale) ||
      !finite_complex(production_value) ||
      !mvmc_scaled_complex_is_valid(left) ||
      !mvmc_scaled_complex_is_valid(right) ||
      !mvmc_scaled_complex_is_valid(guide)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = scaled_log_scale(left, left_log_scale, &scaled_left);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = scaled_log_scale(right, right_log_scale, &scaled_right);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = scaled_conjugate(&scaled_left, &conjugated_left);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    pfaffian_status = mvmc_scaled_complex_multiply(
        &conjugated_left, &scaled_right, &numerator);
    status = pfaffian_status_to_krylov(pfaffian_status, &numerator);
  }
  if (status == MVMC_KRYLOV_STATUS_OK && sign == -1 &&
      numerator.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
    numerator.phase = -numerator.phase;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_scaled_divide(
        &numerator, guide, &quotient);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_scaled_export_evidence(
        &quotient, 0.0, &independent);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_numeric_evidence_recenter(
        &independent, production_value, evidence);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_evidence(status, evidence);
  }
  return status;
}

MVMCKrylovStatus mvmc_power_lanczos_scaled_export_evidence(
    const MVMCScaledComplex *scaled_value, double common_log_scale,
    MVMCPowerLanczosNumericEvidence *evidence) {
  MVMCPowerLanczosNumericEvidence candidate;
  MVMCScaledComplexExportStatus export_status;
  double magnitude = 0.0;
  double propagated_bound = 0.0;
  double rounding_bound = 0.0;
  if (evidence == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_evidence(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, evidence);
  if (!mvmc_scaled_complex_is_valid(scaled_value) ||
      !isfinite(common_log_scale)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (scaled_value->state == MVMC_SCALED_COMPLEX_NONFINITE) {
    invalidate_evidence(MVMC_KRYLOV_STATUS_NONFINITE, evidence);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.version = MVMC_POWER_LANCZOS_NUMERIC_EVIDENCE_VERSION;
  candidate.scaled_state = scaled_value->state;
  export_status = mvmc_scaled_complex_export_common_scale(
      scaled_value, common_log_scale, &candidate.value);
  if (export_status == MVMC_SCALED_EXPORT_EXACT_ZERO) {
    candidate.support_flags =
        MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO;
  } else if (export_status == MVMC_SCALED_EXPORT_NUMERIC_ZERO) {
    if (!exp_upper(scaled_value->log_abs_error_bound - common_log_scale,
                   &candidate.absolute_numeric_bound)) {
      invalidate_evidence(MVMC_KRYLOV_STATUS_NONFINITE, evidence);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    candidate.support_flags =
        MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO;
  } else if (export_status == MVMC_SCALED_EXPORT_UNDERFLOW) {
    if (!exp_upper(scaled_value->log_abs - common_log_scale, &magnitude)) {
      invalidate_evidence(MVMC_KRYLOV_STATUS_NONFINITE, evidence);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    if (isfinite(scaled_value->log_abs_error_bound) &&
        !exp_upper(scaled_value->log_abs_error_bound - common_log_scale,
                   &propagated_bound)) {
      invalidate_evidence(MVMC_KRYLOV_STATUS_NONFINITE, evidence);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    if (!add_upper(magnitude, propagated_bound,
                   &candidate.absolute_numeric_bound)) {
      invalidate_evidence(MVMC_KRYLOV_STATUS_NONFINITE, evidence);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    candidate.support_flags =
        MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO;
  } else if (export_status == MVMC_SCALED_EXPORT_OK) {
    if (!finite_magnitude(candidate.value, &magnitude)) {
      invalidate_evidence(MVMC_KRYLOV_STATUS_NONFINITE, evidence);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    if (isfinite(scaled_value->log_abs_error_bound) &&
        !exp_upper(scaled_value->log_abs_error_bound - common_log_scale,
                   &propagated_bound)) {
      invalidate_evidence(MVMC_KRYLOV_STATUS_NONFINITE, evidence);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    rounding_bound = magnitude * (64.0 * DBL_EPSILON);
    if (rounding_bound == 0.0 && magnitude > 0.0) {
      rounding_bound = nextafter(0.0, 1.0);
    }
    if (!isfinite(rounding_bound) ||
        !add_upper(propagated_bound, rounding_bound,
                   &candidate.absolute_numeric_bound)) {
      invalidate_evidence(MVMC_KRYLOV_STATUS_NONFINITE, evidence);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    candidate.support_flags =
        MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO;
    if (candidate.absolute_numeric_bound >= magnitude) {
      candidate.support_flags |=
          MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO;
    }
  } else {
    const MVMCKrylovStatus status =
        export_status == MVMC_SCALED_EXPORT_OVERFLOW ||
                export_status == MVMC_SCALED_EXPORT_NONFINITE
            ? MVMC_KRYLOV_STATUS_NONFINITE
            : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    invalidate_evidence(status, evidence);
    return status;
  }
  if (!finite_complex(candidate.value) ||
      !isfinite(candidate.absolute_numeric_bound) ||
      candidate.absolute_numeric_bound < 0.0 ||
      ((candidate.support_flags &
        MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO) != 0 &&
       candidate.absolute_numeric_bound <= 0.0)) {
    invalidate_evidence(MVMC_KRYLOV_STATUS_NONFINITE, evidence);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  candidate.valid = 1;
  *evidence = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_numeric_evidence_recenter(
    const MVMCPowerLanczosNumericEvidence *evidence,
    double complex production_value,
    MVMCPowerLanczosNumericEvidence *recentered) {
  MVMCPowerLanczosNumericEvidence candidate;
  double distance;
  double magnitude;
  const uint8_t support_mask =
      MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO |
      MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO |
      MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO;
  if (recentered == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_evidence(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, recentered);
  if (evidence == NULL || !evidence->valid ||
      evidence->status != MVMC_KRYLOV_STATUS_OK ||
      evidence->version != MVMC_POWER_LANCZOS_NUMERIC_EVIDENCE_VERSION ||
      !finite_complex(evidence->value) ||
      !isfinite(evidence->absolute_numeric_bound) ||
      evidence->absolute_numeric_bound < 0.0 ||
      evidence->support_flags == 0 ||
      (evidence->support_flags & (uint8_t)~support_mask) != 0 ||
      !finite_complex(production_value) ||
      !finite_magnitude(production_value - evidence->value, &distance) ||
      !finite_magnitude(production_value, &magnitude)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if ((evidence->support_flags &
       MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO) != 0 &&
      (distance != 0.0 || magnitude != 0.0)) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  candidate = *evidence;
  if (!add_upper(evidence->absolute_numeric_bound, distance,
                 &candidate.absolute_numeric_bound)) {
    invalidate_evidence(MVMC_KRYLOV_STATUS_NONFINITE, recentered);
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  candidate.value = production_value;
  if (magnitude == 0.0) {
    if ((candidate.support_flags &
         MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO) == 0) {
      candidate.support_flags =
          MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO;
      if (candidate.absolute_numeric_bound == 0.0) {
        candidate.absolute_numeric_bound = nextafter(0.0, 1.0);
      }
    }
  } else {
    candidate.support_flags =
        MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO;
    if ((evidence->support_flags &
         MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO) != 0 ||
        candidate.absolute_numeric_bound >= magnitude) {
      candidate.support_flags |=
          MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO;
    }
  }
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  *recentered = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}
