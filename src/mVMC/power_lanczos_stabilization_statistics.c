#include "power_lanczos_stabilization_statistics.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "stabilization statistics require the bounded power-Lanczos engine"
#endif

#include "krylov_streaming_statistics.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef MVMCKrylovStatus (*TraceExport)(
    const MVMCPowerLanczosProductionSession *, double *, size_t);

static int finite_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static int checked_multiply(size_t left, size_t right, size_t *product) {
  if (product == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *product = left * right;
  return 1;
}

static void invalidate_result(
    MVMCKrylovStatus status,
    MVMCPowerLanczosStabilizationResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
}

static MVMCKrylovStatus trace_maximum_tau(
    const MVMCPowerLanczosProductionSession *session,
    const MVMCPowerLanczosPrimitiveTraceSummary *summary,
    TraceExport export_trace, size_t first_scalar, double *maximum_tau) {
  double *scalars = NULL;
  size_t value_count = 0;
  size_t byte_count = 0;
  size_t scalar;
  size_t group;
  MVMCKrylovStatus status;
  if (session == NULL || summary == NULL || !summary->valid ||
      !summary->frozen || summary->scalar_count == 0 ||
      summary->group_count == 0 || summary->samples_per_group < 2 ||
      export_trace == NULL || first_scalar >= summary->scalar_count ||
      maximum_tau == NULL ||
      !checked_multiply(summary->scalar_count, summary->group_count,
                        &value_count) ||
      !checked_multiply(value_count, summary->samples_per_group,
                        &value_count) ||
      !checked_multiply(value_count, sizeof(*scalars), &byte_count)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  scalars = (double *)malloc(byte_count);
  if (scalars == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  status = export_trace(session, scalars, value_count);
  *maximum_tau = 1.0;
  for (scalar = first_scalar; status == MVMC_KRYLOV_STATUS_OK &&
                   scalar < summary->scalar_count;
       ++scalar) {
    for (group = 0; status == MVMC_KRYLOV_STATUS_OK &&
                    group < summary->group_count;
         ++group) {
      const size_t offset =
          (scalar * summary->group_count + group) *
          summary->samples_per_group;
      MVMCKrylovTauIntResult tau;
      status = mvmc_krylov_tau_int_geyer_initial_positive(
          scalars + offset, summary->samples_per_group, &tau);
      if (status == MVMC_KRYLOV_STATUS_OK &&
          (!tau.valid || tau.status != MVMC_KRYLOV_STATUS_OK ||
           !isfinite(tau.tau_int) || tau.tau_int < 1.0)) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      }
      if (status == MVMC_KRYLOV_STATUS_OK &&
          tau.tau_int > *maximum_tau) {
        *maximum_tau = tau.tau_int;
      }
    }
  }
  free(scalars);
  return status;
}

static int support_is_valid(uint8_t support) {
  const uint8_t allowed =
      MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO |
      MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO |
      MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO;
  return support != 0 && (support & (uint8_t)~allowed) == 0;
}

static MVMCKrylovStatus record_primitive(
    double complex value, double bound, uint8_t support,
    MVMCPowerLanczosStabilizationResult *result) {
  if (!finite_complex(value) || !isfinite(bound) || bound < 0.0 ||
      !support_is_valid(support)) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  if ((support & MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO) != 0) {
    ++result->finite_nonzero_primitive_count;
  }
  if ((support & MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO) != 0) {
    ++result->exact_zero_primitive_count;
  }
  if ((support & MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO) != 0) {
    ++result->numeric_zero_primitive_count;
  }
  if (bound > result->maximum_absolute_numeric_bound) {
    result->maximum_absolute_numeric_bound = bound;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus analyze_trace_primitives(
    const MVMCPowerLanczosProductionSession *session,
    const MVMCPowerLanczosProductionSessionSummary *session_summary,
    const MVMCPowerLanczosPrimitiveTraceSummary *coefficient,
    const MVMCPowerLanczosPrimitiveTraceSummary *final,
    MVMCPowerLanczosStabilizationResult *result) {
  double complex *values = NULL;
  double *bounds = NULL;
  uint8_t *support = NULL;
  size_t maximum_primitives;
  size_t group;
  size_t sample;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  if (session == NULL || session_summary == NULL || coefficient == NULL ||
      final == NULL || result == NULL ||
      coefficient->primitive_count != 1 + 4 * session_summary->upper_count ||
      final->primitive_count != 2 ||
      coefficient->group_count != session_summary->chain_count ||
      final->group_count != session_summary->chain_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  maximum_primitives = coefficient->primitive_count > final->primitive_count
                           ? coefficient->primitive_count
                           : final->primitive_count;
  values = (double complex *)calloc(maximum_primitives, sizeof(*values));
  bounds = (double *)calloc(maximum_primitives, sizeof(*bounds));
  support = (uint8_t *)calloc(maximum_primitives, sizeof(*support));
  if (values == NULL || bounds == NULL || support == NULL) {
    status = MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    goto done;
  }
  for (group = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  group < coefficient->group_count;
       ++group) {
    for (sample = 0; status == MVMC_KRYLOV_STATUS_OK &&
                     sample < coefficient->samples_per_group;
         ++sample) {
      double tail = NAN;
      size_t primitive;
      int outside_base_support = 0;
      status = mvmc_power_lanczos_production_session_coefficient_trace_sample(
          session, group, sample, values, coefficient->primitive_count,
          bounds, coefficient->primitive_count, support,
          coefficient->primitive_count, &tail);
      if (status == MVMC_KRYLOV_STATUS_OK && tail != 0.0) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      }
      for (primitive = 0; status == MVMC_KRYLOV_STATUS_OK &&
                          primitive < coefficient->primitive_count;
           ++primitive) {
        status = record_primitive(values[primitive], bounds[primitive],
                                  support[primitive], result);
      }
      if (status != MVMC_KRYLOV_STATUS_OK) continue;
      outside_base_support =
          (support[1] &
           (MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO |
            MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO)) != 0 &&
          (support[1] &
           MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO) == 0;
      if (outside_base_support) {
        size_t entry;
        ++result->outside_base_support_samples;
        for (entry = 1; entry < session_summary->upper_count; ++entry) {
          const size_t primitive = 1 + entry;
          if ((support[primitive] &
               MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO) != 0 &&
              cabs(values[primitive]) > bounds[primitive]) {
            ++result->nonzero_outside_support_contributions;
          }
        }
      }
    }
  }
  for (group = 0; status == MVMC_KRYLOV_STATUS_OK &&
                  group < final->group_count;
       ++group) {
    for (sample = 0; status == MVMC_KRYLOV_STATUS_OK &&
                     sample < final->samples_per_group;
         ++sample) {
      double tail = NAN;
      size_t primitive;
      status = mvmc_power_lanczos_production_session_final_trace_sample(
          session, group, sample, values, final->primitive_count,
          bounds, final->primitive_count, support,
          final->primitive_count, &tail);
      if (status == MVMC_KRYLOV_STATUS_OK &&
          (!isfinite(tail) || tail < 0.0 ||
           tail != hypot(creal(values[1]), cimag(values[1])))) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      }
      for (primitive = 0; status == MVMC_KRYLOV_STATUS_OK &&
                          primitive < final->primitive_count;
           ++primitive) {
        status = record_primitive(values[primitive], bounds[primitive],
                                  support[primitive], result);
      }
    }
  }

done:
  free(support);
  free(bounds);
  free(values);
  return status;
}

static MVMCKrylovStatus evaluate_blocks(
    const MVMCPowerLanczosProductionSession *session,
    const MVMCPowerLanczosProductionSessionSummary *summary,
    MVMCPowerLanczosStabilizationResult *result) {
  double complex *leave_one_second_moment = NULL;
  double complex *final_block_mean = NULL;
  double complex leave_one_mean = 0.0;
  double complex energy_weighted_sum = 0.0;
  double real_jackknife_sum = 0.0;
  double imag_jackknife_sum = 0.0;
  double real_energy_sum = 0.0;
  double imag_energy_sum = 0.0;
  uint64_t final_count = 0;
  uint64_t expected_block_count = 0;
  size_t block;
  double block_count;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  if (session == NULL || summary == NULL || result == NULL ||
      summary->block_count < 2) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  leave_one_second_moment = (double complex *)calloc(
      summary->block_count, sizeof(*leave_one_second_moment));
  final_block_mean = (double complex *)calloc(
      summary->block_count, sizeof(*final_block_mean));
  if (leave_one_second_moment == NULL || final_block_mean == NULL) {
    status = MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    goto done;
  }
  for (block = 0; block < summary->block_count; ++block) {
    MVMCPowerLanczosProductionCoefficientBlockResult coefficient;
    MVMCPowerLanczosProductionFinalBlockResult final;
    status = mvmc_power_lanczos_production_session_coefficient_block_result(
        session, block, &coefficient);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_power_lanczos_production_session_final_block_result(
          session, block, &final);
    }
    if (status != MVMC_KRYLOV_STATUS_OK) goto done;
    if (block == 0) {
      expected_block_count = final.sample_count;
    } else if (final.sample_count != expected_block_count) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      goto done;
    }
    if (final.sample_count > UINT64_MAX - final_count ||
        !finite_complex(coefficient.leave_one_second_moment) ||
        !finite_complex(final.energy_sum)) {
      status = MVMC_KRYLOV_STATUS_NONFINITE;
      goto done;
    }
    leave_one_second_moment[block] =
        coefficient.leave_one_second_moment;
    leave_one_mean += coefficient.leave_one_second_moment;
    final_block_mean[block] = final.energy_sum / (double)final.sample_count;
    energy_weighted_sum += final.energy_sum;
    final_count += final.sample_count;
  }
  block_count = (double)summary->block_count;
  leave_one_mean /= block_count;
  energy_weighted_sum /= (double)final_count;
  for (block = 0; block < summary->block_count; ++block) {
    const double complex leave_delta =
        leave_one_second_moment[block] - leave_one_mean;
    const double complex energy_delta =
        final_block_mean[block] - energy_weighted_sum;
    real_jackknife_sum += creal(leave_delta) * creal(leave_delta);
    imag_jackknife_sum += cimag(leave_delta) * cimag(leave_delta);
    real_energy_sum += creal(energy_delta) * creal(energy_delta);
    imag_energy_sum += cimag(energy_delta) * cimag(energy_delta);
  }
  result->energy = creal(energy_weighted_sum);
  result->energy_imaginary = cimag(energy_weighted_sum);
  result->energy_standard_error =
      sqrt(real_energy_sum / (block_count * (block_count - 1.0)));
  result->energy_imaginary_standard_error =
      sqrt(imag_energy_sum / (block_count * (block_count - 1.0)));
  result->coefficient_second_moment =
      summary->corrected_variance +
      creal(energy_weighted_sum * energy_weighted_sum);
  result->coefficient_second_moment_imaginary =
      summary->corrected_variance_imaginary +
      cimag(energy_weighted_sum * energy_weighted_sum);
  result->coefficient_second_moment_standard_error =
      sqrt((block_count - 1.0) / block_count * real_jackknife_sum);
  result->coefficient_second_moment_imaginary_standard_error =
      sqrt((block_count - 1.0) / block_count * imag_jackknife_sum);
  result->variance = summary->corrected_variance;
  result->variance_imaginary = summary->corrected_variance_imaginary;
  result->variance_standard_error = hypot(
      hypot(result->coefficient_second_moment_standard_error,
            result->coefficient_second_moment_imaginary_standard_error),
      2.0 * cabs(energy_weighted_sum) *
          hypot(result->energy_standard_error,
                result->energy_imaginary_standard_error));
  if (!isfinite(result->energy) || !isfinite(result->energy_imaginary) ||
      !isfinite(result->energy_standard_error) ||
      !isfinite(result->energy_imaginary_standard_error) ||
      !isfinite(result->coefficient_second_moment) ||
      !isfinite(result->coefficient_second_moment_imaginary) ||
      !isfinite(result->coefficient_second_moment_standard_error) ||
      !isfinite(result->coefficient_second_moment_imaginary_standard_error) ||
      !isfinite(result->variance) || !isfinite(result->variance_imaginary) ||
      !isfinite(result->variance_standard_error) ||
      fabs(result->energy - summary->final_energy) >
          64.0 * DBL_EPSILON * fmax(1.0, fabs(summary->final_energy)) ||
      fabs(result->energy_imaginary - summary->final_energy_imaginary) >
          64.0 * DBL_EPSILON *
              fmax(1.0, fabs(summary->final_energy_imaginary))) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }

done:
  free(final_block_mean);
  free(leave_one_second_moment);
  return status;
}

MVMCKrylovStatus mvmc_power_lanczos_stabilization_statistics_evaluate(
    const MVMCPowerLanczosProductionSession *session,
    int require_zero_support_evidence,
    MVMCPowerLanczosStabilizationResult *result) {
  MVMCPowerLanczosStabilizationResult candidate;
  MVMCPowerLanczosProductionSessionSummary summary;
  MVMCPowerLanczosPrimitiveTraceSummary coefficient_trace;
  MVMCPowerLanczosPrimitiveTraceSummary final_trace;
  const uint64_t inconclusive_gates =
      MVMC_POWER_LANCZOS_STABILIZATION_GATE_BLOCK_COUNT |
      MVMC_POWER_LANCZOS_STABILIZATION_GATE_COEFFICIENT_TAU |
      MVMC_POWER_LANCZOS_STABILIZATION_GATE_FINAL_TAU;
  double imaginary_numeric_tolerance;
  double energy_imaginary_tolerance;
  double variance_imaginary_tolerance;
  MVMCKrylovStatus status;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_result(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (session == NULL ||
      (require_zero_support_evidence != 0 &&
       require_zero_support_evidence != 1)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  status = mvmc_power_lanczos_production_session_summary(session, &summary);
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (!summary.valid || summary.status != MVMC_KRYLOV_STATUS_OK ||
       summary.state != MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINALIZED ||
       summary.observable_enabled || summary.block_count < 2 ||
       summary.chain_count == 0)) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_production_session_coefficient_trace_summary(
        session, &coefficient_trace);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_production_session_final_trace_summary(
        session, &final_trace);
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (coefficient_trace.samples_per_group % summary.block_count != 0 ||
       final_trace.samples_per_group % summary.block_count != 0)) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    candidate.block_count = summary.block_count;
    candidate.chain_count = summary.chain_count;
    candidate.coefficient_sample_count = summary.coefficient_sample_count;
    candidate.final_sample_count = summary.final_sample_count;
    candidate.coefficient_block_length =
        coefficient_trace.samples_per_group / summary.block_count;
    candidate.final_block_length =
        final_trace.samples_per_group / summary.block_count;
    candidate.maximum_leave_one_projective_distance =
        summary.maximum_leave_one_projective_distance;
    candidate.retained_rank = summary.coefficient_gevp.retained_rank;
    candidate.gevp_relative_residual =
        summary.coefficient_gevp.relative_residual;
    candidate.gevp_normalization_absolute_error =
        fabs(summary.coefficient_gevp.normalization - 1.0);
    candidate.session_allocated_bytes = summary.allocated_bytes;
    status = evaluate_blocks(session, &summary, &candidate);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    /* Primitive zero is the importance-normalization denominator.  Its
     * eta-scale variation can have an arbitrarily large scale-invariant tau
     * while contributing negligibly to corrected S/K/B.  The block ratios
     * retain that denominator; the mixing maximum begins at primitive one. */
    status = trace_maximum_tau(
        session, &coefficient_trace,
        mvmc_power_lanczos_production_session_coefficient_trace_export_scalars,
        2,
        &candidate.coefficient_tau_int_max);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = trace_maximum_tau(
        session, &final_trace,
        mvmc_power_lanczos_production_session_final_trace_export_scalars,
        0,
        &candidate.final_tau_int_max);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = analyze_trace_primitives(
        session, &summary, &coefficient_trace, &final_trace, &candidate);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    invalidate_result(status, result);
    return status;
  }
  candidate.coefficient_tau_to_block_length =
      candidate.coefficient_tau_int_max /
      (double)candidate.coefficient_block_length;
  candidate.final_tau_to_block_length =
      candidate.final_tau_int_max / (double)candidate.final_block_length;
  candidate.coefficient_effective_block_count =
      (double)candidate.block_count /
      fmax(1.0, candidate.coefficient_tau_to_block_length);
  candidate.final_effective_block_count =
      (double)candidate.block_count /
      fmax(1.0, candidate.final_tau_to_block_length);
  if (candidate.block_count < 32) {
    candidate.failed_gates |=
        MVMC_POWER_LANCZOS_STABILIZATION_GATE_BLOCK_COUNT;
  }
  if (candidate.coefficient_tau_to_block_length > 1.0 / 24.0) {
    candidate.failed_gates |=
        MVMC_POWER_LANCZOS_STABILIZATION_GATE_COEFFICIENT_TAU;
  }
  if (candidate.final_tau_to_block_length > 1.0 / 24.0) {
    candidate.failed_gates |=
        MVMC_POWER_LANCZOS_STABILIZATION_GATE_FINAL_TAU;
  }
  if (candidate.retained_rank < 1) {
    candidate.failed_gates |=
        MVMC_POWER_LANCZOS_STABILIZATION_GATE_GEVP_RANK;
  }
  if (!isfinite(candidate.gevp_relative_residual) ||
      candidate.gevp_relative_residual > 64.0 * DBL_EPSILON) {
    candidate.failed_gates |=
        MVMC_POWER_LANCZOS_STABILIZATION_GATE_GEVP_RESIDUAL;
  }
  if (!isfinite(candidate.gevp_normalization_absolute_error) ||
      candidate.gevp_normalization_absolute_error >
          MVMC_POWER_LANCZOS_GEVP_NORMALIZATION_TOLERANCE_FACTOR *
              DBL_EPSILON) {
    candidate.failed_gates |=
        MVMC_POWER_LANCZOS_STABILIZATION_GATE_GEVP_NORMALIZATION;
  }
  if (candidate.variance < -candidate.maximum_absolute_numeric_bound) {
    candidate.failed_gates |=
        MVMC_POWER_LANCZOS_STABILIZATION_GATE_NEGATIVE_VARIANCE;
  }
  if (require_zero_support_evidence &&
      (candidate.outside_base_support_samples == 0 ||
       candidate.nonzero_outside_support_contributions == 0)) {
    candidate.failed_gates |=
        MVMC_POWER_LANCZOS_STABILIZATION_GATE_ZERO_SUPPORT;
  }
  imaginary_numeric_tolerance = candidate.maximum_absolute_numeric_bound +
      64.0 * DBL_EPSILON *
          fmax(1.0, fmax(fabs(candidate.energy), fabs(candidate.variance)));
  energy_imaginary_tolerance = imaginary_numeric_tolerance +
      MVMC_POWER_LANCZOS_STABILIZATION_IMAGINARY_SIGMA_MULTIPLIER *
          candidate.energy_imaginary_standard_error;
  variance_imaginary_tolerance = imaginary_numeric_tolerance +
      MVMC_POWER_LANCZOS_STABILIZATION_IMAGINARY_SIGMA_MULTIPLIER *
          candidate.variance_standard_error;
  if (fabs(candidate.energy_imaginary) > energy_imaginary_tolerance ||
      fabs(candidate.variance_imaginary) > variance_imaginary_tolerance) {
    candidate.failed_gates |=
        MVMC_POWER_LANCZOS_STABILIZATION_GATE_COMPLEX_RESIDUAL;
  }
  if (candidate.failed_gates == 0) {
    candidate.decision = MVMC_POWER_LANCZOS_STABILIZATION_PASS;
  } else if ((candidate.failed_gates & ~inconclusive_gates) == 0) {
    candidate.decision = MVMC_POWER_LANCZOS_STABILIZATION_INCONCLUSIVE;
  } else {
    candidate.decision = MVMC_POWER_LANCZOS_STABILIZATION_FAIL;
  }
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.version = MVMC_POWER_LANCZOS_STABILIZATION_RESULT_VERSION;
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

const char *mvmc_power_lanczos_stabilization_decision_string(
    MVMCPowerLanczosStabilizationDecision decision) {
  switch (decision) {
    case MVMC_POWER_LANCZOS_STABILIZATION_PASS:
      return "PASS";
    case MVMC_POWER_LANCZOS_STABILIZATION_INCONCLUSIVE:
      return "INCONCLUSIVE";
    case MVMC_POWER_LANCZOS_STABILIZATION_FAIL:
      return "FAIL";
    default:
      return "INVALID";
  }
}
