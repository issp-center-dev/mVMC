#include "power_lanczos_observable_evaluator.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

#define CHECK(condition, message)                       \
  do {                                                  \
    if (!(condition)) {                                 \
      fprintf(stderr, "FAIL: %s (line %d)\n", message, \
              __LINE__);                                \
      ++failures;                                       \
    }                                                   \
  } while (0)

typedef struct {
  double complex value[4];
  uint64_t calls;
  MVMCKrylovStatus forced_status;
} AmplitudeTable;

static MVMCKrylovStatus TableAmplitude(
    const uint64_t *words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
  AmplitudeTable *table = (AmplitudeTable *)context;
  const size_t index = word_count == 1 ? (size_t)(words[0] & 3) : 0;
  MVMCPfaffianStatus pf_status;
  ++table->calls;
  if (table->forced_status != MVMC_KRYLOV_STATUS_OK) {
    return table->forced_status;
  }
  memset(result, 0, sizeof(*result));
  if (creal(table->value[index]) == 0.0 &&
      cimag(table->value[index]) == 0.0) {
    pf_status = mvmc_scaled_complex_make_exact_zero(&result->value);
  } else {
    pf_status = mvmc_scaled_complex_from_raw(table->value[index],
                                             &result->value);
  }
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  if (result->value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO) {
    result->exact_zero_component_count = 1;
  } else {
    result->well_pivoted_component_count = 1;
  }
  result->local_factorization_count = 1;
  result->global_factorization_count = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static int Close(double complex actual, double complex expected) {
  return cabs(actual - expected) <= 5.0e-14 * (1.0 + cabs(expected));
}

static MVMCKrylovBoundedLimits Limits(void) {
  MVMCKrylovBoundedLimits limits;
  memset(&limits, 0, sizeof(limits));
  limits.amplitude_policy_hash = UINT64_C(0x503643324556414c);
  limits.cache_bytes = 1024 * 1024;
  limits.max_row_transitions = 32;
  limits.max_workspace_bytes = 64 * 1024 * 1024;
  limits.max_node_expansions = 1024;
  limits.max_terminal_amplitude_calls = 1024;
  limits.max_total_row_transitions = 4096;
  limits.max_order = 1;
  return limits;
}

int main(void) {
  const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0}};
  const MVMCKrylovHamiltonianTerm terms[] = {
      {1.0, 0, 2, 0, 0}, {1.0, 2, 2, 0, 1}};
  const MVMCKrylovFockModel model = {
      2, 1, 0, 0, 1, terms, 2, operators, 4};
  const MVMCPowerLanczosObservableLayout layout = {2, 1, 0, 0};
  const int raw[4] = {0, 0, 1, 0};
  MVMCPowerLanczosObservableRecord record;
  MVMCPowerLanczosObservablePlan observable_plan;
  const uint64_t source = UINT64_C(2);
  const double scales[2] = {0.0, 0.0};
  const double inverse_sqrt_two = 1.0 / sqrt(2.0);
  const double complex alpha[2] = {inverse_sqrt_two, inverse_sqrt_two};
  const MVMCKrylovBoundedLimits limits = Limits();
  MVMCKrylovBoundedPlan *bounded_plan = NULL;
  MVMCKrylovBoundedWorkspace *bounded_workspace = NULL;
  MVMCPowerLanczosObservableEvaluatorWorkspace *evaluator = NULL;
  MVMCPowerLanczosObservableSampleDiagnostics diagnostics;
  AmplitudeTable table;
  double complex matrix[4] = {99.0, 99.0, 99.0, 99.0};
  double complex numerator = 99.0;
  MVMCKrylovStatus status;
  memset(&record, 0, sizeof(record));
  record.family = MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS;
  record.row_width = 4;
  memcpy(record.raw_indices, raw, sizeof(raw));
  memset(&observable_plan, 0, sizeof(observable_plan));
  observable_plan.nsite = 2;
  observable_plan.record_count = 1;
  observable_plan.records = &record;
  memset(&table, 0, sizeof(table));
  table.value[2] = 1.0;

  status = mvmc_bounded_krylov_plan_create(&model, &limits, &bounded_plan);
  CHECK(status == MVMC_KRYLOV_STATUS_OK, "bounded plan create");
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_workspace_create(
        bounded_plan, &bounded_workspace);
  }
  CHECK(status == MVMC_KRYLOV_STATUS_OK, "bounded workspace create");
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_evaluator_workspace_create(
        1, 1, &evaluator);
  }
  CHECK(status == MVMC_KRYLOV_STATUS_OK && evaluator != NULL,
        "evaluator workspace create");
  CHECK(mvmc_power_lanczos_observable_evaluator_workspace_bytes(evaluator) >
            0,
        "evaluator workspace byte ledger");

  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_session_begin(
        bounded_workspace, TableAmplitude, &table, UINT64_C(0x1001));
  }
  CHECK(status == MVMC_KRYLOV_STATUS_OK, "coefficient session begin");
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_coefficient_sample(
        evaluator, bounded_workspace, &layout, &observable_plan, &source, 1,
        0.0, scales, matrix, 4, &diagnostics);
  }
  CHECK(status == MVMC_KRYLOV_STATUS_OK && diagnostics.valid,
        "coefficient sample");
  CHECK(Close(matrix[0], 0.0) && Close(matrix[1], 0.0) &&
            Close(matrix[2], 1.0) && Close(matrix[3], 0.0),
        "full lower-triangle matrix counterexample");
  CHECK(diagnostics.request_count == 1 &&
            diagnostics.active_request_count == 1 &&
            diagnostics.unique_target_count == 1 &&
            diagnostics.engine_root_evaluations == 2 &&
            !diagnostics.source_target_reused,
        "coefficient diagnostics");
  if (mvmc_bounded_krylov_session_is_active(bounded_workspace)) {
    CHECK(mvmc_bounded_krylov_session_end(bounded_workspace) ==
              MVMC_KRYLOV_STATUS_OK,
          "coefficient session end");
  }

  status = mvmc_bounded_krylov_session_begin(
      bounded_workspace, TableAmplitude, &table, UINT64_C(0x1002));
  CHECK(status == MVMC_KRYLOV_STATUS_OK, "final session begin");
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_final_sample(
        evaluator, bounded_workspace, &layout, &observable_plan, &source, 1,
        0.0, scales, alpha, &numerator, 1, &diagnostics);
  }
  CHECK(status == MVMC_KRYLOV_STATUS_OK && diagnostics.valid,
        "final sample");
  if (!Close(numerator, 0.5)) {
    fprintf(stderr, "direct-final actual=(%.17g,%.17g)\n",
            creal(numerator), cimag(numerator));
  }
  CHECK(Close(numerator, 0.5), "direct-final lower-triangle numerator");
  if (mvmc_bounded_krylov_session_is_active(bounded_workspace)) {
    CHECK(mvmc_bounded_krylov_session_end(bounded_workspace) ==
              MVMC_KRYLOV_STATUS_OK,
          "final session end");
  }

  matrix[0] = 123.0;
  memset(&diagnostics, 0x5a, sizeof(diagnostics));
  status = mvmc_bounded_krylov_session_begin(
      bounded_workspace, TableAmplitude, &table, UINT64_C(0x1003));
  CHECK(status == MVMC_KRYLOV_STATUS_OK, "failure session begin");
  table.forced_status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_coefficient_sample(
        evaluator, bounded_workspace, &layout, &observable_plan, &source, 1,
        0.0, scales, matrix, 4, &diagnostics);
  }
  CHECK(status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE,
        "engine failure propagation");
  CHECK(Close(matrix[0], 123.0), "failed evaluator exposed partial output");
  CHECK(!mvmc_bounded_krylov_session_is_active(bounded_workspace),
        "failed engine session remained active");

  table.forced_status = MVMC_KRYLOV_STATUS_OK;
  status = mvmc_bounded_krylov_session_begin(
      bounded_workspace, TableAmplitude, &table, UINT64_C(0x1004));
  CHECK(status == MVMC_KRYLOV_STATUS_OK, "invalid-output session begin");
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_coefficient_sample(
        evaluator, bounded_workspace, &layout, &observable_plan, &source, 1,
        0.0, scales, NULL, 4, &diagnostics);
  }
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "NULL coefficient output returned success");
  if (mvmc_bounded_krylov_session_is_active(bounded_workspace)) {
    CHECK(mvmc_bounded_krylov_session_end(bounded_workspace) ==
              MVMC_KRYLOV_STATUS_OK,
          "invalid-output session end");
  }

  mvmc_power_lanczos_observable_evaluator_workspace_destroy(evaluator);
  mvmc_bounded_krylov_workspace_destroy(bounded_workspace);
  mvmc_bounded_krylov_plan_destroy(bounded_plan);
  if (failures != 0) {
    fprintf(stderr, "%d observable evaluator assertion(s) failed\n",
            failures);
    return EXIT_FAILURE;
  }
  puts("power_lanczos_observable_evaluator_unit: PASS");
  return EXIT_SUCCESS;
}
