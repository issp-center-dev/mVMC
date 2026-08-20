/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "krylov_positive_guide.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "krylov_positive_guide.c requires the power-Lanczos core"
#endif

#include <math.h>
#include <string.h>

static double log_add_positive(double left, double right) {
  const double maximum = fmax(left, right);
  return maximum + log1p(exp(fmin(left, right) - maximum));
}

static int validate_policy(const MVMCKrylovPositiveGuidePolicy *policy) {
  int order;
  if (policy == NULL || policy->order < 0 ||
      policy->order > MVMC_KRYLOV_MAX_ORDER || !isfinite(policy->eta) ||
      policy->eta <= 0.0) {
    return 0;
  }
  for (order = 0; order <= policy->order; ++order) {
    if (!isfinite(policy->lambda[order]) || policy->lambda[order] < 0.0) {
      return 0;
    }
    if (!isfinite(policy->log_basis_scale[order])) return 0;
  }
  return 1;
}

static void invalidate_guide_result(MVMCKrylovStatus status,
                                    MVMCKrylovPositiveGuideResult *result) {
  if (result == NULL) return;
  memset(result, 0, sizeof(*result));
  result->status = status;
  result->log_guide = NAN;
  result->log_floor = NAN;
  result->log_signal_sum = NAN;
}

MVMCKrylovStatus mvmc_krylov_positive_guide_evaluate(
    const MVMCKrylovPositiveGuidePolicy *policy,
    const MVMCScaledComplex *values, size_t value_count,
    MVMCKrylovPositiveGuideResult *result) {
  MVMCKrylovPositiveGuideResult candidate;
  int order;
  double log_guide;
  double log_signal_sum = -INFINITY;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  invalidate_guide_result(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
  if (!validate_policy(policy) || values == NULL ||
      value_count < (size_t)policy->order + 1) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }

  memset(&candidate, 0, sizeof(candidate));
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.log_floor = log(policy->eta);
  log_guide = candidate.log_floor;

  for (order = 0; order <= policy->order; ++order) {
    const MVMCScaledComplex *value = &values[order];
    double term_log;
    if (!mvmc_scaled_complex_is_valid(value)) {
      invalidate_guide_result(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    if (value->state == MVMC_SCALED_COMPLEX_NONFINITE) {
      invalidate_guide_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    if (value->state == MVMC_SCALED_COMPLEX_EXACT_ZERO) {
      ++candidate.exact_zero_component_count;
      continue;
    }
    if (value->state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO) {
      ++candidate.numeric_zero_component_count;
      continue;
    }
    ++candidate.finite_component_count;
    if (policy->lambda[order] == 0.0) continue;
    if (!isfinite(value->log_abs + policy->log_basis_scale[order])) {
      invalidate_guide_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    term_log = log(policy->lambda[order]) +
               2.0 * (value->log_abs + policy->log_basis_scale[order]);
    if (!isfinite(term_log)) {
      invalidate_guide_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    log_signal_sum = isinf(log_signal_sum) && log_signal_sum < 0.0
                         ? term_log
                         : log_add_positive(log_signal_sum, term_log);
    log_guide = log_add_positive(log_guide, term_log);
    if (!isfinite(log_signal_sum) || !isfinite(log_guide)) {
      invalidate_guide_result(MVMC_KRYLOV_STATUS_NONFINITE, result);
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }

  candidate.valid = 1;
  candidate.log_guide = log_guide;
  candidate.log_signal_sum = log_signal_sum;
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_krylov_positive_guide_acceptance(
    const MVMCKrylovPositiveGuideResult *current,
    const MVMCKrylovPositiveGuideResult *proposal,
    double log_proposal_ratio, double uniform,
    MVMCKrylovPositiveGuideAcceptance *result) {
  MVMCKrylovPositiveGuideAcceptance candidate;
  double log_acceptance_ratio;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(result, 0, sizeof(*result));
  result->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  result->log_acceptance_ratio = NAN;
  result->uniform = NAN;

  if (current == NULL || proposal == NULL || !current->valid ||
      !proposal->valid || current->status != MVMC_KRYLOV_STATUS_OK ||
      proposal->status != MVMC_KRYLOV_STATUS_OK ||
      !isfinite(current->log_guide) || !isfinite(proposal->log_guide) ||
      !isfinite(log_proposal_ratio) || !isfinite(uniform) ||
      uniform < 0.0 || uniform >= 1.0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }

  log_acceptance_ratio =
      proposal->log_guide - current->log_guide + log_proposal_ratio;
  if (!isfinite(log_acceptance_ratio)) {
    result->status = MVMC_KRYLOV_STATUS_NONFINITE;
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }

  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.status = MVMC_KRYLOV_STATUS_OK;
  candidate.log_acceptance_ratio = log_acceptance_ratio;
  candidate.uniform = uniform;
  candidate.accepted = log_acceptance_ratio >= 0.0 ||
                       uniform == 0.0 ||
                       log(uniform) < log_acceptance_ratio;
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}
