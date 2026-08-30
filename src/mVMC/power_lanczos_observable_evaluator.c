#include "power_lanczos_observable_evaluator.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "power_lanczos_observable_evaluator.c requires the bounded engine"
#endif

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

struct MVMCPowerLanczosObservableEvaluatorWorkspace {
  size_t max_requests;
  size_t word_count;
  size_t allocated_bytes;
  MVMCPowerLanczosObservableTargetGroupScratch group_scratch;
  uint64_t *target_words;
  int *request_target_index;
  int *request_sign;
  MVMCScaledComplex *target_basis;
  double complex *primitive_scratch;
  MVMCPowerLanczosNumericEvidence *evidence_scratch;
};

typedef struct {
  int zero;
  double complex phase;
  double log_abs;
} ScaledView;

static int CheckedAdd(size_t *total, size_t addition) {
  if (addition > SIZE_MAX - *total) return 0;
  *total += addition;
  return 1;
}

static int CheckedMultiply(size_t left, size_t right, size_t *product) {
  if (left != 0 && right > SIZE_MAX / left) return 0;
  *product = left * right;
  return 1;
}

static int TargetHashCapacity(size_t max_requests, size_t *capacity) {
  size_t candidate = 2;
  if (capacity == NULL || max_requests > SIZE_MAX / 2) return 0;
  while (candidate < 2 * max_requests) {
    if (candidate > SIZE_MAX / 2) return 0;
    candidate *= 2;
  }
  *capacity = candidate;
  return 1;
}

static int FiniteComplex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static MVMCKrylovStatus AllocateArray(size_t count, size_t element_size,
                                     void **allocation,
                                     size_t *allocated_bytes) {
  size_t bytes;
  if (!CheckedMultiply(count, element_size, &bytes) ||
      !CheckedAdd(allocated_bytes, bytes)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  *allocation = calloc(count, element_size);
  return *allocation != NULL ? MVMC_KRYLOV_STATUS_OK
                             : MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
}

MVMCKrylovStatus mvmc_power_lanczos_observable_evaluator_workspace_create(
    size_t max_requests, size_t word_count,
    MVMCPowerLanczosObservableEvaluatorWorkspace **workspace) {
  MVMCPowerLanczosObservableEvaluatorWorkspace *candidate;
  MVMCKrylovStatus status;
  size_t target_word_count;
  size_t target_basis_count;
  size_t primitive_count;
  size_t target_hash_capacity;
  if (workspace == NULL || *workspace != NULL || max_requests == 0 ||
      max_requests > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS ||
      word_count == 0 ||
      word_count > ((2 * MVMC_POWER_LANCZOS_OBSERVABLE_MAX_NSITE + 63) /
                    64) ||
      !CheckedMultiply(max_requests, word_count, &target_word_count) ||
      !TargetHashCapacity(max_requests, &target_hash_capacity) ||
      !CheckedMultiply(max_requests, 2, &target_basis_count) ||
      !CheckedMultiply(max_requests,
                       MVMC_POWER_LANCZOS_OBSERVABLE_MATRIX_ENTRIES,
                       &primitive_count)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  candidate = (MVMCPowerLanczosObservableEvaluatorWorkspace *)calloc(
      1, sizeof(*candidate));
  if (candidate == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  candidate->allocated_bytes = sizeof(*candidate);
  candidate->group_scratch.max_requests = max_requests;
  candidate->group_scratch.word_count = word_count;
  candidate->group_scratch.hash_capacity = target_hash_capacity;
  status = AllocateArray(
      target_word_count, sizeof(candidate->group_scratch.target_words[0]),
      (void **)&candidate->group_scratch.target_words,
      &candidate->allocated_bytes);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = AllocateArray(
        max_requests,
        sizeof(candidate->group_scratch.request_target_index[0]),
        (void **)&candidate->group_scratch.request_target_index,
        &candidate->allocated_bytes);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = AllocateArray(
        max_requests,
        sizeof(candidate->group_scratch.request_fermion_sign[0]),
        (void **)&candidate->group_scratch.request_fermion_sign,
        &candidate->allocated_bytes);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = AllocateArray(
        target_hash_capacity,
        sizeof(candidate->group_scratch.hash_slots[0]),
        (void **)&candidate->group_scratch.hash_slots,
        &candidate->allocated_bytes);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = AllocateArray(target_word_count,
                           sizeof(candidate->target_words[0]),
                           (void **)&candidate->target_words,
                           &candidate->allocated_bytes);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = AllocateArray(max_requests,
                           sizeof(candidate->request_target_index[0]),
                           (void **)&candidate->request_target_index,
                           &candidate->allocated_bytes);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = AllocateArray(max_requests, sizeof(candidate->request_sign[0]),
                           (void **)&candidate->request_sign,
                           &candidate->allocated_bytes);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = AllocateArray(target_basis_count,
                           sizeof(candidate->target_basis[0]),
                           (void **)&candidate->target_basis,
                           &candidate->allocated_bytes);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = AllocateArray(primitive_count,
                           sizeof(candidate->primitive_scratch[0]),
                           (void **)&candidate->primitive_scratch,
                           &candidate->allocated_bytes);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = AllocateArray(primitive_count,
                           sizeof(candidate->evidence_scratch[0]),
                           (void **)&candidate->evidence_scratch,
                           &candidate->allocated_bytes);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    mvmc_power_lanczos_observable_evaluator_workspace_destroy(candidate);
    return status;
  }
  candidate->max_requests = max_requests;
  candidate->word_count = word_count;
  *workspace = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_power_lanczos_observable_evaluator_workspace_destroy(
    MVMCPowerLanczosObservableEvaluatorWorkspace *workspace) {
  if (workspace == NULL) return;
  free(workspace->evidence_scratch);
  free(workspace->primitive_scratch);
  free(workspace->target_basis);
  free(workspace->request_sign);
  free(workspace->request_target_index);
  free(workspace->target_words);
  free(workspace->group_scratch.hash_slots);
  free(workspace->group_scratch.request_fermion_sign);
  free(workspace->group_scratch.request_target_index);
  free(workspace->group_scratch.target_words);
  memset(workspace, 0, sizeof(*workspace));
  free(workspace);
}

size_t mvmc_power_lanczos_observable_evaluator_workspace_bytes(
    const MVMCPowerLanczosObservableEvaluatorWorkspace *workspace) {
  return workspace != NULL ? workspace->allocated_bytes : 0;
}

static MVMCKrylovStatus ActionStatusToKrylov(
    MVMCPowerLanczosObservableActionStatus status) {
  switch (status) {
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK:
      return MVMC_KRYLOV_STATUS_OK;
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_RESOURCE_LIMIT:
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INTERNAL_INVARIANT_FAILURE:
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_ARGUMENT:
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_LAYOUT:
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_SOURCE:
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_RECORD:
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_TARGET:
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

static MVMCKrylovStatus ValidateCommon(
    MVMCPowerLanczosObservableEvaluatorWorkspace *workspace,
    MVMCKrylovBoundedWorkspace *bounded_workspace,
    const MVMCPowerLanczosObservableLayout *layout,
    const MVMCPowerLanczosObservablePlan *plan,
    const uint64_t *source_words, size_t word_count, double log_guide,
    const double basis_log_scale[2], size_t output_capacity,
    size_t output_entries_per_request) {
  size_t required_output;
  if (workspace == NULL || bounded_workspace == NULL || layout == NULL ||
      plan == NULL || source_words == NULL || basis_log_scale == NULL ||
      !mvmc_bounded_krylov_session_is_active(bounded_workspace) ||
      plan->record_count < 0 ||
      (size_t)plan->record_count > workspace->max_requests ||
      word_count != workspace->word_count || !isfinite(log_guide) ||
      !isfinite(basis_log_scale[0]) || !isfinite(basis_log_scale[1]) ||
      !CheckedMultiply((size_t)plan->record_count,
                       output_entries_per_request, &required_output) ||
      output_capacity < required_output) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus EvaluateRoot(
    MVMCKrylovBoundedWorkspace *bounded_workspace,
    const uint64_t *words, size_t word_count,
    MVMCScaledComplex basis[2]) {
  MVMCKrylovBoundedResult result;
  MVMCKrylovStatus status = mvmc_bounded_krylov_session_evaluate(
      bounded_workspace, words, word_count, &result);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (!result.valid || result.status != MVMC_KRYLOV_STATUS_OK ||
      result.evaluated_order < 1 ||
      !mvmc_scaled_complex_is_valid(&result.value[0]) ||
      !mvmc_scaled_complex_is_valid(&result.value[1])) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  basis[0] = result.value[0];
  basis[1] = result.value[1];
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus EvaluateSourceAndTargets(
    MVMCPowerLanczosObservableEvaluatorWorkspace *workspace,
    MVMCKrylovBoundedWorkspace *bounded_workspace,
    const MVMCPowerLanczosObservableLayout *layout,
    const MVMCPowerLanczosObservablePlan *plan,
    const uint64_t *source_words, size_t word_count,
    MVMCScaledComplex source_basis[2],
    MVMCPowerLanczosObservableSampleDiagnostics *diagnostics) {
  MVMCPowerLanczosObservableTargetGroupResult group;
  MVMCPowerLanczosObservableActionStatus action_status;
  MVMCKrylovStatus status;
  size_t target_index;
  action_status = mvmc_power_lanczos_observable_group_targets(
      layout, plan, &workspace->group_scratch, source_words, word_count,
      workspace->target_words,
      workspace->max_requests * workspace->word_count,
      workspace->request_target_index, workspace->max_requests,
      workspace->request_sign, workspace->max_requests, &group);
  if (action_status != MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK) {
    return ActionStatusToKrylov(action_status);
  }
  status = EvaluateRoot(bounded_workspace, source_words, word_count,
                        source_basis);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  diagnostics->engine_root_evaluations = 1;
  for (target_index = 0; target_index < group.unique_target_count;
       ++target_index) {
    const uint64_t *target =
        workspace->target_words + target_index * word_count;
    MVMCScaledComplex *target_basis =
        workspace->target_basis + 2 * target_index;
    if (memcmp(source_words, target,
               word_count * sizeof(source_words[0])) == 0) {
      target_basis[0] = source_basis[0];
      target_basis[1] = source_basis[1];
      diagnostics->source_target_reused = 1;
      continue;
    }
    status = EvaluateRoot(bounded_workspace, target, word_count,
                          target_basis);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    ++diagnostics->engine_root_evaluations;
  }
  diagnostics->request_count = group.request_count;
  diagnostics->active_request_count = group.active_request_count;
  diagnostics->exact_zero_request_count =
      group.request_count - group.active_request_count;
  diagnostics->unique_target_count = group.unique_target_count;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus MakeView(const MVMCScaledComplex *value,
                                double log_scale, ScaledView *view) {
  if (value == NULL || view == NULL || !isfinite(log_scale) ||
      !mvmc_scaled_complex_is_valid(value)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(view, 0, sizeof(*view));
  switch (value->state) {
    case MVMC_SCALED_COMPLEX_FINITE_NONZERO:
      if (!FiniteComplex(value->phase) || !isfinite(value->log_abs) ||
          !isfinite(value->log_abs + log_scale)) {
        return MVMC_KRYLOV_STATUS_NONFINITE;
      }
      view->phase = value->phase;
      view->log_abs = value->log_abs + log_scale;
      return MVMC_KRYLOV_STATUS_OK;
    case MVMC_SCALED_COMPLEX_EXACT_ZERO:
    case MVMC_SCALED_COMPLEX_NUMERIC_ZERO:
      view->zero = 1;
      return MVMC_KRYLOV_STATUS_OK;
    case MVMC_SCALED_COMPLEX_NONFINITE:
      return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

static MVMCKrylovStatus SafeExp(double log_value, double *value) {
  if (value == NULL || !isfinite(log_value)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (log_value > log(DBL_MAX)) return MVMC_KRYLOV_STATUS_NONFINITE;
  if (log_value < log(nextafter(0.0, 1.0))) {
    *value = 0.0;
    return MVMC_KRYLOV_STATUS_OK;
  }
  *value = exp(log_value);
  return isfinite(*value) ? MVMC_KRYLOV_STATUS_OK
                          : MVMC_KRYLOV_STATUS_NONFINITE;
}

static MVMCKrylovStatus GuideNormalizedProduct(
    const ScaledView *left, const ScaledView *right, int sign,
    double log_guide, double complex *result, double *log_contribution) {
  double log_magnitude;
  double magnitude;
  MVMCKrylovStatus status;
  if (left == NULL || right == NULL || result == NULL ||
      log_contribution == NULL || (sign != -1 && sign != 1) ||
      !isfinite(log_guide)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (left->zero || right->zero) {
    *result = 0.0;
    *log_contribution = -INFINITY;
    return MVMC_KRYLOV_STATUS_OK;
  }
  log_magnitude = left->log_abs + right->log_abs - log_guide;
  if (!isfinite(log_magnitude)) return MVMC_KRYLOV_STATUS_NONFINITE;
  status = SafeExp(log_magnitude, &magnitude);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  *result = (double)sign * magnitude * conj(left->phase) * right->phase;
  if (!FiniteComplex(*result)) return MVMC_KRYLOV_STATUS_NONFINITE;
  *log_contribution = log_magnitude;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus ExactZeroEvidence(
    MVMCPowerLanczosNumericEvidence *evidence) {
  MVMCScaledComplex zero;
  MVMCPfaffianStatus pfaffian_status;
  if (evidence == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  pfaffian_status = mvmc_scaled_complex_make_exact_zero(&zero);
  if (pfaffian_status != MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  return mvmc_power_lanczos_scaled_export_evidence(&zero, 0.0, evidence);
}

static void NoteLog(double value,
                    MVMCPowerLanczosObservableSampleDiagnostics *diagnostics) {
  if (value == -INFINITY) return;
  if (!isfinite(diagnostics->minimum_log_contribution) ||
      value < diagnostics->minimum_log_contribution) {
    diagnostics->minimum_log_contribution = value;
  }
  if (!isfinite(diagnostics->maximum_log_contribution) ||
      value > diagnostics->maximum_log_contribution) {
    diagnostics->maximum_log_contribution = value;
  }
}

static MVMCKrylovStatus ScaleBasisValue(
    const MVMCScaledComplex *basis, double log_scale,
    double complex coefficient, MVMCScaledComplex *result) {
  MVMCScaledComplex scale;
  MVMCScaledComplex scaled_basis;
  MVMCScaledComplex scaled_coefficient;
  MVMCPfaffianStatus pf_status;
  if (basis == NULL || result == NULL || !isfinite(log_scale) ||
      !FiniteComplex(coefficient)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (creal(coefficient) == 0.0 && cimag(coefficient) == 0.0) {
    pf_status = mvmc_scaled_complex_make_exact_zero(result);
    return pf_status == MVMC_PFAFFIAN_STATUS_OK
               ? MVMC_KRYLOV_STATUS_OK
               : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  pf_status = mvmc_scaled_complex_make_finite(
      1.0, log_scale, -INFINITY, &scale);
  if (pf_status == MVMC_PFAFFIAN_STATUS_OK) {
    pf_status = mvmc_scaled_complex_multiply(basis, &scale, &scaled_basis);
  }
  if (pf_status == MVMC_PFAFFIAN_STATUS_OK) {
    pf_status = mvmc_scaled_complex_from_raw(coefficient,
                                             &scaled_coefficient);
  }
  if (pf_status == MVMC_PFAFFIAN_STATUS_OK) {
    pf_status = mvmc_scaled_complex_multiply(
        &scaled_coefficient, &scaled_basis, result);
  }
  return pf_status == MVMC_PFAFFIAN_STATUS_OK
             ? MVMC_KRYLOV_STATUS_OK
             : MVMC_KRYLOV_STATUS_NONFINITE;
}

static MVMCKrylovStatus MakeFinalAmplitude(
    const MVMCScaledComplex basis[2], const double basis_log_scale[2],
    const double complex alpha[2], MVMCScaledComplex *result) {
  MVMCScaledComplex terms[2];
  MVMCKrylovStatus status;
  MVMCPfaffianStatus pf_status;
  int index;
  for (index = 0; index < 2; ++index) {
    status = ScaleBasisValue(&basis[index], basis_log_scale[index],
                             alpha[index], &terms[index]);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
  }
  pf_status = mvmc_scaled_complex_sum_ordered(terms, 2, result);
  return pf_status == MVMC_PFAFFIAN_STATUS_OK
             ? MVMC_KRYLOV_STATUS_OK
             : MVMC_KRYLOV_STATUS_NONFINITE;
}

MVMCKrylovStatus mvmc_power_lanczos_observable_coefficient_sample(
    MVMCPowerLanczosObservableEvaluatorWorkspace *workspace,
    MVMCKrylovBoundedWorkspace *bounded_workspace,
    const MVMCPowerLanczosObservableLayout *layout,
    const MVMCPowerLanczosObservablePlan *plan,
    const uint64_t *source_words, size_t word_count,
    double log_guide, const double basis_log_scale[2],
    double complex *matrix_entries, size_t matrix_entry_capacity,
    MVMCPowerLanczosObservableSampleDiagnostics *diagnostics) {
  return mvmc_power_lanczos_observable_coefficient_sample_with_evidence(
      workspace, bounded_workspace, layout, plan, source_words, word_count,
      log_guide, basis_log_scale, NULL, matrix_entries,
      matrix_entry_capacity, NULL, 0, diagnostics);
}

MVMCKrylovStatus
mvmc_power_lanczos_observable_coefficient_sample_with_evidence(
    MVMCPowerLanczosObservableEvaluatorWorkspace *workspace,
    MVMCKrylovBoundedWorkspace *bounded_workspace,
    const MVMCPowerLanczosObservableLayout *layout,
    const MVMCPowerLanczosObservablePlan *plan,
    const uint64_t *source_words, size_t word_count,
    double log_guide, const double basis_log_scale[2],
    const MVMCScaledComplex *guide,
    double complex *matrix_entries, size_t matrix_entry_capacity,
    MVMCPowerLanczosNumericEvidence *evidence,
    size_t evidence_capacity,
    MVMCPowerLanczosObservableSampleDiagnostics *diagnostics) {
  MVMCScaledComplex source_basis[2];
  MVMCPowerLanczosObservableSampleDiagnostics local_diagnostics;
  MVMCKrylovStatus status;
  size_t request;
  memset(&local_diagnostics, 0, sizeof(local_diagnostics));
  local_diagnostics.status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  local_diagnostics.minimum_log_contribution = NAN;
  local_diagnostics.maximum_log_contribution = NAN;
  status = ValidateCommon(workspace, bounded_workspace, layout, plan,
                          source_words, word_count, log_guide,
                          basis_log_scale, matrix_entry_capacity,
                          MVMC_POWER_LANCZOS_OBSERVABLE_MATRIX_ENTRIES);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (matrix_entries == NULL || diagnostics == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if ((evidence == NULL) != (guide == NULL) ||
      (evidence != NULL &&
       (evidence_capacity <
            (size_t)plan->record_count *
                MVMC_POWER_LANCZOS_OBSERVABLE_MATRIX_ENTRIES ||
        !mvmc_scaled_complex_is_valid(guide) ||
        guide->state != MVMC_SCALED_COMPLEX_FINITE_NONZERO))) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = EvaluateSourceAndTargets(
      workspace, bounded_workspace, layout, plan, source_words, word_count,
      source_basis, &local_diagnostics);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  for (request = 0; request < (size_t)plan->record_count; ++request) {
    size_t row;
    const int target_index = workspace->request_target_index[request];
    if (target_index == MVMC_POWER_LANCZOS_OBSERVABLE_ZERO_TARGET) {
      for (row = 0; row < 4; ++row) {
        workspace->primitive_scratch[4 * request + row] = 0.0;
        if (evidence != NULL) {
          status = ExactZeroEvidence(
              workspace->evidence_scratch + 4 * request + row);
          if (status != MVMC_KRYLOV_STATUS_OK) return status;
        }
      }
      continue;
    }
    for (row = 0; row < 2; ++row) {
      size_t column;
      ScaledView left;
      status = MakeView(&workspace->target_basis[2 * (size_t)target_index + row],
                        basis_log_scale[row], &left);
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      for (column = 0; column < 2; ++column) {
        ScaledView right;
        double log_contribution;
        status = MakeView(&source_basis[column], basis_log_scale[column],
                          &right);
        if (status != MVMC_KRYLOV_STATUS_OK) return status;
        status = GuideNormalizedProduct(
            &left, &right, workspace->request_sign[request], log_guide,
            &workspace->primitive_scratch[4 * request + 2 * row + column],
            &log_contribution);
        if (status != MVMC_KRYLOV_STATUS_OK) return status;
        if (evidence != NULL) {
          status =
              mvmc_power_lanczos_guide_normalized_product_evidence(
                  &workspace->target_basis[
                      2 * (size_t)target_index + row],
                  basis_log_scale[row], &source_basis[column],
                  basis_log_scale[column],
                  workspace->request_sign[request], guide,
                  workspace->primitive_scratch[
                      4 * request + 2 * row + column],
                  workspace->evidence_scratch +
                      4 * request + 2 * row + column);
          if (status != MVMC_KRYLOV_STATUS_OK) return status;
        }
        NoteLog(log_contribution, &local_diagnostics);
      }
    }
  }
  memcpy(matrix_entries, workspace->primitive_scratch,
         (size_t)plan->record_count * 4 * sizeof(matrix_entries[0]));
  if (evidence != NULL) {
    memcpy(evidence, workspace->evidence_scratch,
           (size_t)plan->record_count * 4 * sizeof(evidence[0]));
  }
  local_diagnostics.valid = 1;
  local_diagnostics.status = MVMC_KRYLOV_STATUS_OK;
  *diagnostics = local_diagnostics;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_observable_final_sample(
    MVMCPowerLanczosObservableEvaluatorWorkspace *workspace,
    MVMCKrylovBoundedWorkspace *bounded_workspace,
    const MVMCPowerLanczosObservableLayout *layout,
    const MVMCPowerLanczosObservablePlan *plan,
    const uint64_t *source_words, size_t word_count,
    double log_guide, const double basis_log_scale[2],
    const double complex alpha[2],
    double complex *numerators, size_t numerator_capacity,
    MVMCPowerLanczosObservableSampleDiagnostics *diagnostics) {
  return mvmc_power_lanczos_observable_final_sample_with_evidence(
      workspace, bounded_workspace, layout, plan, source_words, word_count,
      log_guide, basis_log_scale, alpha, numerators, numerator_capacity,
      NULL, 0, diagnostics);
}

MVMCKrylovStatus mvmc_power_lanczos_observable_final_sample_with_evidence(
    MVMCPowerLanczosObservableEvaluatorWorkspace *workspace,
    MVMCKrylovBoundedWorkspace *bounded_workspace,
    const MVMCPowerLanczosObservableLayout *layout,
    const MVMCPowerLanczosObservablePlan *plan,
    const uint64_t *source_words, size_t word_count,
    double log_guide, const double basis_log_scale[2],
    const double complex alpha[2],
    double complex *numerators, size_t numerator_capacity,
    MVMCPowerLanczosNumericEvidence *evidence,
    size_t evidence_capacity,
    MVMCPowerLanczosObservableSampleDiagnostics *diagnostics) {
  MVMCScaledComplex source_basis[2];
  MVMCScaledComplex source_final;
  MVMCScaledComplex evidence_guide;
  MVMCPowerLanczosObservableSampleDiagnostics local_diagnostics;
  MVMCKrylovStatus status;
  size_t request;
  memset(&local_diagnostics, 0, sizeof(local_diagnostics));
  local_diagnostics.status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  local_diagnostics.minimum_log_contribution = NAN;
  local_diagnostics.maximum_log_contribution = NAN;
  status = ValidateCommon(workspace, bounded_workspace, layout, plan,
                          source_words, word_count, log_guide,
                          basis_log_scale, numerator_capacity, 1);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (alpha == NULL || !FiniteComplex(alpha[0]) ||
      !FiniteComplex(alpha[1]) || numerators == NULL ||
      diagnostics == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (evidence != NULL &&
      evidence_capacity < (size_t)plan->record_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = EvaluateSourceAndTargets(
      workspace, bounded_workspace, layout, plan, source_words, word_count,
      source_basis, &local_diagnostics);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  status = MakeFinalAmplitude(source_basis, basis_log_scale, alpha,
                              &source_final);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (evidence != NULL) {
    const double evidence_log_scale[1] = {0.0};
    const double evidence_weight[1] = {1.0};
    status = mvmc_power_lanczos_scaled_weighted_norm(
        &source_final, evidence_log_scale, evidence_weight, 1, 0.0,
        &evidence_guide);
    if (status != MVMC_KRYLOV_STATUS_OK ||
        evidence_guide.state != MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
      return status == MVMC_KRYLOV_STATUS_OK
                 ? MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE
                 : status;
    }
  }
  for (request = 0; request < (size_t)plan->record_count; ++request) {
    const int target_index = workspace->request_target_index[request];
    if (target_index == MVMC_POWER_LANCZOS_OBSERVABLE_ZERO_TARGET) {
      workspace->primitive_scratch[request] = 0.0;
      if (evidence != NULL) {
        status = ExactZeroEvidence(workspace->evidence_scratch + request);
        if (status != MVMC_KRYLOV_STATUS_OK) return status;
      }
    } else {
      MVMCScaledComplex target_final;
      ScaledView left;
      ScaledView right;
      double log_contribution;
      status = MakeFinalAmplitude(
          &workspace->target_basis[2 * (size_t)target_index],
          basis_log_scale, alpha, &target_final);
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      status = MakeView(&target_final, 0.0, &left);
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = MakeView(&source_final, 0.0, &right);
      }
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = GuideNormalizedProduct(
            &left, &right, workspace->request_sign[request], log_guide,
            &workspace->primitive_scratch[request], &log_contribution);
      }
      if (status != MVMC_KRYLOV_STATUS_OK) return status;
      if (evidence != NULL) {
        status = mvmc_power_lanczos_guide_normalized_product_evidence(
            &target_final, 0.0, &source_final, 0.0,
            workspace->request_sign[request], &evidence_guide,
            workspace->primitive_scratch[request],
            workspace->evidence_scratch + request);
        if (status != MVMC_KRYLOV_STATUS_OK) return status;
      }
      NoteLog(log_contribution, &local_diagnostics);
    }
  }
  memcpy(numerators, workspace->primitive_scratch,
         (size_t)plan->record_count * sizeof(numerators[0]));
  if (evidence != NULL) {
    memcpy(evidence, workspace->evidence_scratch,
           (size_t)plan->record_count * sizeof(evidence[0]));
  }
  local_diagnostics.valid = 1;
  local_diagnostics.status = MVMC_KRYLOV_STATUS_OK;
  *diagnostics = local_diagnostics;
  return MVMC_KRYLOV_STATUS_OK;
}
