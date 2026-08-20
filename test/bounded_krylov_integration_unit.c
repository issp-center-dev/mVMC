#include "bounded_krylov_collective.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failures = 0;
static int world_rank = 0;
static int world_size = 1;

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "BoundedKrylovIntegration_Unit FAIL rank %d: ",         \
              world_rank);                                                      \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, " (line %d)\n", __LINE__);                              \
      ++failures;                                                               \
    }                                                                           \
  } while (0)

typedef struct {
  MVMCKrylovBoundedCollectiveWorkspace *collective;
  double complex weights[4];
  int qp_start;
  int qp_end;
  int failure_rank;
  uint64_t calls;
} CollectiveAmplitudeContext;

static MVMCAbsolutePfaffianScaledValueResult component_value(
    double complex raw) {
  MVMCAbsolutePfaffianScaledValueResult component;
  memset(&component, 0, sizeof(component));
  component.matrix_scale = 1.0;
  component.scaled_min_pivot = 0.5;
  if (creal(raw) == 0.0 && cimag(raw) == 0.0) {
    component.factor_state = MVMC_PFAFFIAN_VALUE_SINGULAR;
    component.factor_info = 1;
    (void)mvmc_scaled_complex_make_exact_zero(&component.value);
  } else {
    component.factor_state = MVMC_PFAFFIAN_VALUE_WELL_PIVOTED;
    (void)mvmc_scaled_complex_from_raw_testing(raw, &component.value);
  }
  return component;
}

static double complex terminal_component(int qp, size_t configuration) {
  static const double complex configuration_one[4] = {
      1.0 + 0.5 * I, 2.0 - I, 0.0, -0.25 + 0.75 * I};
  static const double complex configuration_two[4] = {
      -0.5 + 0.25 * I, 0.0, 3.0 + I, 0.75 - 0.5 * I};
  return configuration == 1 ? configuration_one[qp]
                            : configuration_two[qp];
}

static MVMCKrylovStatus collective_amplitude_callback(
    const uint64_t *words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
  CollectiveAmplitudeContext *callback =
      (CollectiveAmplitudeContext *)context;
  MVMCAbsolutePfaffianScaledValueResult local[4];
  MVMCKrylovBoundedCollectiveResult collective_result;
  MVMCKrylovStatus local_status = MVMC_KRYLOV_STATUS_OK;
  const size_t configuration =
      word_count == 1 ? (size_t)(words[0] & UINT64_C(63)) : 63;
  int qp;
  ++callback->calls;
  for (qp = callback->qp_start; qp < callback->qp_end; ++qp) {
    local[qp - callback->qp_start] =
        component_value(terminal_component(qp, configuration));
  }
  if (world_rank == callback->failure_rank && callback->calls == 1) {
    local_status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  return mvmc_bounded_krylov_collective_gather_projected_amplitude(
      callback->collective, local_status, local,
      (size_t)(callback->qp_end - callback->qp_start), callback->weights,
      4, callback->qp_start, callback->qp_end, result,
      &collective_result);
}

static MVMCKrylovFockModel hopping_model(
    const MVMCKrylovHamiltonianTerm *terms,
    const MVMCKrylovFermionOperator *operators) {
  MVMCKrylovFockModel model;
  memset(&model, 0, sizeof(model));
  model.site_count = 2;
  model.up_electron_count = 1;
  model.hermitian = 1;
  model.terms = terms;
  model.term_count = 2;
  model.operators = operators;
  model.operator_count = 4;
  return model;
}

static MVMCKrylovBoundedLimits bounded_limits(void) {
  MVMCKrylovBoundedLimits limits;
  limits.amplitude_policy_hash = UINT64_C(0x5150312d34504f4c);
  limits.cache_bytes = 4096;
  limits.max_row_transitions = 8;
  limits.max_workspace_bytes = (size_t)1024 * 1024;
  limits.max_node_expansions = 1000;
  limits.max_terminal_amplitude_calls = 1000;
  limits.max_total_row_transitions = 10000;
  limits.max_order = 3;
  return limits;
}

static int scaled_to_raw(const MVMCScaledComplex *value,
                         double complex *raw) {
  if (value->state == MVMC_SCALED_COMPLEX_EXACT_ZERO) {
    *raw = 0.0;
    return 1;
  }
  return mvmc_scaled_complex_export_common_scale(value, 0.0, raw) ==
         MVMC_SCALED_EXPORT_OK;
}

static int close_complex(double complex actual, double complex expected,
                         double tolerance) {
  return cabs(actual - expected) <= tolerance * (1.0 + cabs(expected));
}

static double complex projected_terminal(
    const CollectiveAmplitudeContext *context, size_t configuration) {
  MVMCAbsolutePfaffianScaledValueResult components[4];
  MVMCProjectedScaledAmplitudeResult projected;
  double complex raw = NAN + I * NAN;
  int qp;
  for (qp = 0; qp < 4; ++qp) {
    components[qp] = component_value(terminal_component(qp, configuration));
  }
  CHECK(mvmc_projected_scaled_amplitude_values(
            components, context->weights, 4, &projected) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            scaled_to_raw(&projected.total, &raw),
        "integration projection oracle failed");
  return raw;
}

static void check_rank_invariance(const MVMCKrylovBoundedResult *result,
                                  uint64_t calls, uint64_t plan_hash) {
#ifdef _mpi_use
  uint64_t local[7];
  uint64_t minimum[7];
  uint64_t maximum[7];
  int order;
  local[0] = plan_hash;
  local[1] = calls;
  local[2] = result->statistics.trace_hash;
  local[3] = result->statistics.node_expansions;
  local[4] = result->statistics.raw_row_transitions;
  local[5] = result->statistics.terminal_amplitude_calls;
  local[6] = (uint64_t)(unsigned int)result->status;
  MPI_Allreduce(local, minimum, 7, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(local, maximum, 7, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  CHECK(memcmp(minimum, maximum, sizeof(minimum)) == 0,
        "integrated plan/ROW/callback trace differs by rank");
  for (order = 0; order <= 3; ++order) {
    double fields[6] = {
        creal(result->value[order].phase),
        cimag(result->value[order].phase),
        result->value[order].log_abs,
        result->value[order].log_abs_error_bound,
        result->value[order].max_input_log_abs,
        result->value[order].cancellation_ratio,
    };
    double field_min[6];
    double field_max[6];
    int state = (int)result->value[order].state;
    int state_min;
    int state_max;
    MPI_Allreduce(fields, field_min, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(fields, field_max, 6, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&state, &state_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&state, &state_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    CHECK(memcmp(field_min, field_max, sizeof(field_min)) == 0 &&
              state_min == state_max,
          "integrated scaled value differs by rank at order %d", order);
  }
#else
  (void)result;
  (void)calls;
  (void)plan_hash;
#endif
}

static void run_integration(
    MVMCKrylovBoundedCollectiveWorkspace *collective) {
  const double complex hopping = 1.25 + 0.5 * I;
  const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
  };
  const MVMCKrylovHamiltonianTerm terms[] = {
      {hopping, 0, 2, 1, 0}, {conj(hopping), 2, 2, 1, 1}};
  const MVMCKrylovFockModel model = hopping_model(terms, operators);
  const MVMCKrylovBoundedLimits limits = bounded_limits();
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  MVMCKrylovBoundedCollectiveResult preflight;
  MVMCKrylovBoundedResult result;
  CollectiveAmplitudeContext context;
  MVMCKrylovStatus local_status;
  MVMCKrylovStatus status;
  const uint64_t root = UINT64_C(1);
  uint64_t plan_hash = 0;
  double complex terminal_one;
  double complex terminal_two;
  double complex expected[4];
  int order;
  memset(&context, 0, sizeof(context));
  context.collective = collective;
  context.weights[0] = 1.0;
  context.weights[1] = -0.5 + 0.25 * I;
  context.weights[2] = 2.0;
  context.weights[3] = -I;
  context.qp_start = 4 * world_rank / world_size;
  context.qp_end = 4 * (world_rank + 1) / world_size;
  context.failure_rank = -1;

  local_status =
      mvmc_bounded_krylov_plan_create(&model, &limits, &plan);
  if (local_status == MVMC_KRYLOV_STATUS_OK) {
    plan_hash = mvmc_bounded_krylov_plan_hash(plan);
    local_status = mvmc_bounded_krylov_workspace_create(plan, &workspace);
  }
  status = mvmc_bounded_krylov_collective_synchronize(
      collective, local_status, plan_hash, workspace != NULL, &preflight);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && preflight.valid,
        "integrated plan/allocation preflight failed");
  if (status != MVMC_KRYLOV_STATUS_OK || workspace == NULL) {
    mvmc_bounded_krylov_workspace_destroy(workspace);
    mvmc_bounded_krylov_plan_destroy(plan);
    return;
  }
  status = mvmc_bounded_krylov_evaluate(
      workspace, &root, 1, collective_amplitude_callback, &context, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid &&
            result.statistics.engine_heap_allocations == 0,
        "integrated allocation-free bounded evaluation failed");
  terminal_one = projected_terminal(&context, 1);
  terminal_two = projected_terminal(&context, 2);
  expected[0] = terminal_one;
  expected[1] = hopping * terminal_two;
  expected[2] = hopping * conj(hopping) * terminal_one;
  expected[3] = hopping * conj(hopping) * hopping * terminal_two;
  for (order = 0; order <= 3; ++order) {
    double complex raw = NAN + I * NAN;
    CHECK(scaled_to_raw(&result.value[order], &raw) &&
              close_complex(raw, expected[order], 5.0e-13),
          "integrated QP=4 result mismatch at order %d", order);
  }
  check_rank_invariance(&result, context.calls, plan_hash);

  context.calls = 0;
  context.failure_rank = world_size - 1;
  status = mvmc_bounded_krylov_evaluate(
      workspace, &root, 1, collective_amplitude_callback, &context, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE && !result.valid &&
            result.evaluated_order == -1,
        "one-rank integrated callback failure did not fail atomically");
  check_rank_invariance(&result, context.calls, plan_hash);
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
}

int main(int argc, char **argv) {
  MVMCKrylovBoundedCollectiveWorkspace *collective = NULL;
  MVMCKrylovStatus status;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  status = mvmc_bounded_krylov_collective_workspace_create(
      MPI_COMM_WORLD, 4, 4, &collective);
#else
  (void)argc;
  (void)argv;
  status = mvmc_bounded_krylov_collective_workspace_create(0, 4, 4,
                                                            &collective);
#endif
  CHECK(status == MVMC_KRYLOV_STATUS_OK && collective != NULL,
        "integration collective workspace creation failed");
  if (collective != NULL) run_integration(collective);
  mvmc_bounded_krylov_collective_workspace_destroy(collective);
#ifdef _mpi_use
  {
    int global_failures = 0;
    MPI_Allreduce(&failures, &global_failures, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    failures = global_failures;
  }
#endif
  if (failures != 0) {
    if (world_rank == 0) {
      fprintf(stderr, "BoundedKrylovIntegration_Unit: %d failure(s)\n",
              failures);
    }
#ifdef _mpi_use
    MPI_Finalize();
#endif
    return EXIT_FAILURE;
  }
  if (world_rank == 0) {
    printf("BoundedKrylovIntegration_Unit: PASS (%d rank%s)\n", world_size,
           world_size == 1 ? "" : "s");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}
