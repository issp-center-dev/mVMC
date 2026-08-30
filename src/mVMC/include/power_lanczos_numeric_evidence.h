#ifndef MVMC_POWER_LANCZOS_NUMERIC_EVIDENCE_H
#define MVMC_POWER_LANCZOS_NUMERIC_EVIDENCE_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "bounded_krylov_engine.h"
#include "power_lanczos_primitive_trace.h"

#include <complex.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_NUMERIC_EVIDENCE_VERSION UINT64_C(1)

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  uint64_t version;
  MVMCScaledComplexState scaled_state;
  double complex value;
  double absolute_numeric_bound;
  uint8_t support_flags;
} MVMCPowerLanczosNumericEvidence;

/*
 * Divide two scaled values while retaining an absolute error bound.  The
 * denominator must be finite nonzero and its lower magnitude bound must be
 * strictly positive.  Exact-zero and numeric-zero numerators are preserved.
 */
MVMCKrylovStatus mvmc_power_lanczos_scaled_divide(
    const MVMCScaledComplex *numerator,
    const MVMCScaledComplex *denominator,
    MVMCScaledComplex *quotient);

/*
 * Form eta + sum_i weight[i] * |exp(log_scale[i]) * value[i]|^2 in
 * scaled form.  The nonnegative binary64 weights and eta are treated as the
 * production inputs; arithmetic error and every input value error are
 * propagated into the result.
 */
MVMCKrylovStatus mvmc_power_lanczos_scaled_weighted_norm(
    const MVMCScaledComplex *values, const double *log_scale,
    const double *weight, size_t value_count, double eta,
    MVMCScaledComplex *result);

/*
 * Certify sign * conj(exp(left_log_scale) * left) *
 * exp(right_log_scale) * right / guide.  The returned center is recentered
 * to production_value without changing that established production value.
 */
MVMCKrylovStatus mvmc_power_lanczos_guide_normalized_product_evidence(
    const MVMCScaledComplex *left, double left_log_scale,
    const MVMCScaledComplex *right, double right_log_scale, int sign,
    const MVMCScaledComplex *guide, double complex production_value,
    MVMCPowerLanczosNumericEvidence *evidence);

/*
 * Export at common_log_scale and convert the scaled error to one finite
 * absolute binary64 bound.  Finite values that underflow during export are
 * reported as NUMERIC_ZERO, never EXACT_ZERO.
 */
MVMCKrylovStatus mvmc_power_lanczos_scaled_export_evidence(
    const MVMCScaledComplex *scaled_value, double common_log_scale,
    MVMCPowerLanczosNumericEvidence *evidence);

/*
 * Recenter an evidence bound around a separately computed production value.
 * This permits evidence calculation to remain independent without changing
 * established estimator bytes.  The central-value distance is added to the
 * absolute bound with outward rounding.
 */
MVMCKrylovStatus mvmc_power_lanczos_numeric_evidence_recenter(
    const MVMCPowerLanczosNumericEvidence *evidence,
    double complex production_value,
    MVMCPowerLanczosNumericEvidence *recentered);

#endif /* bounded engine */

#endif /* MVMC_POWER_LANCZOS_NUMERIC_EVIDENCE_H */
