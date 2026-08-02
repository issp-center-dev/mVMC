#include "bounded_krylov_engine.h"

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

_Static_assert(sizeof(double) == sizeof(uint64_t),
               "bitwise scaled-value tests require 64-bit double storage");

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "BoundedKrylovEngine_Unit FAIL rank %d: ",              \
              world_rank);                                                      \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, " (line %d)\n", __LINE__);                              \
      ++failures;                                                               \
    }                                                                           \
  } while (0)

typedef struct {
  double complex values[64];
  uint64_t calls;
  uint64_t trace_hash;
  MVMCKrylovStatus forced_status;
  int return_nonfinite;
} TableAmplitude;

static uint64_t hash_u64_value(uint64_t hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    hash ^= value & UINT64_C(0xff);
    hash *= UINT64_C(1099511628211);
    value >>= 8;
  }
  return hash;
}

static void scaled_value_field_bits(const MVMCScaledComplex *value,
                                    uint64_t bits[8]) {
  double fields[7];
  size_t field;
  bits[0] = (uint64_t)(unsigned int)value->state;
  fields[0] = creal(value->phase);
  fields[1] = cimag(value->phase);
  fields[2] = value->log_abs;
  fields[3] = value->log_abs_error_bound;
  fields[4] = value->max_input_log_abs;
  fields[5] = value->cancellation_log_abs;
  fields[6] = value->cancellation_ratio;
  for (field = 0; field < 7; ++field) {
    memcpy(&bits[field + 1], &fields[field], sizeof(bits[field + 1]));
  }
}

static int scaled_values_bitwise_equal(const MVMCScaledComplex *left,
                                       const MVMCScaledComplex *right) {
  uint64_t left_bits[8];
  uint64_t right_bits[8];
  scaled_value_field_bits(left, left_bits);
  scaled_value_field_bits(right, right_bits);
  return memcmp(left_bits, right_bits, sizeof(left_bits)) == 0;
}

static size_t small_configuration_index(const uint64_t *words,
                                        size_t word_count) {
  return word_count == 1 ? (size_t)(words[0] & UINT64_C(63)) : 63;
}

static MVMCKrylovStatus scaled_table_callback(
    const uint64_t *words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
  TableAmplitude *table = (TableAmplitude *)context;
  const size_t index = small_configuration_index(words, word_count);
  ++table->calls;
  table->trace_hash = hash_u64_value(table->trace_hash, (uint64_t)index);
  if (table->forced_status != MVMC_KRYLOV_STATUS_OK) {
    return table->forced_status;
  }
  memset(result, 0, sizeof(*result));
  if (table->return_nonfinite) {
    memset(&result->value, 0, sizeof(result->value));
    result->value.state = MVMC_SCALED_COMPLEX_NONFINITE;
    return MVMC_KRYLOV_STATUS_OK;
  }
  if (creal(table->values[index]) == 0.0 &&
      cimag(table->values[index]) == 0.0) {
    if (mvmc_scaled_complex_make_exact_zero(&result->value) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    result->exact_zero_component_count = 1;
  } else {
    if (mvmc_scaled_complex_from_raw_testing(table->values[index],
                                              &result->value) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    }
    result->well_pivoted_component_count = 1;
  }
  result->local_factorization_count = 1;
  result->global_factorization_count = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus raw_table_callback(
    const uint64_t *words, size_t word_count, void *context,
    MVMCKrylovAmplitudeResult *result) {
  TableAmplitude *table = (TableAmplitude *)context;
  const size_t index = small_configuration_index(words, word_count);
  ++table->calls;
  memset(result, 0, sizeof(*result));
  if (table->forced_status != MVMC_KRYLOV_STATUS_OK) {
    return table->forced_status;
  }
  result->value = table->return_nonfinite ? NAN + 0.0 * I
                                           : table->values[index];
  result->total_zero = creal(result->value) == 0.0 &&
                       cimag(result->value) == 0.0;
  result->regular_component_count = result->total_zero ? 0 : 1;
  result->local_factorization_count = 1;
  result->global_factorization_count = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovFockModel electronic_model(
    size_t site_count, size_t up_count, size_t down_count,
    const MVMCKrylovHamiltonianTerm *terms, size_t term_count,
    const MVMCKrylovFermionOperator *operators, size_t operator_count) {
  MVMCKrylovFockModel model;
  memset(&model, 0, sizeof(model));
  model.site_count = site_count;
  model.up_electron_count = up_count;
  model.down_electron_count = down_count;
  model.hermitian = 1;
  model.terms = terms;
  model.term_count = term_count;
  model.operators = operators;
  model.operator_count = operator_count;
  return model;
}

static MVMCKrylovBoundedLimits bounded_limits(int order,
                                               size_t cache_bytes) {
  MVMCKrylovBoundedLimits limits;
  limits.amplitude_policy_hash = UINT64_C(0x504f4c4943590001);
  limits.cache_bytes = cache_bytes;
  limits.max_row_transitions = 4096;
  limits.max_workspace_bytes = (size_t)512 * 1024 * 1024;
  limits.max_node_expansions = UINT64_C(1000000);
  limits.max_terminal_amplitude_calls = UINT64_C(1000000);
  limits.max_total_row_transitions = UINT64_C(10000000);
  limits.max_order = order;
  return limits;
}

static MVMCKrylovLimits reference_limits(int order) {
  MVMCKrylovLimits limits;
  limits.max_states = 4096;
  limits.max_transitions = 65536;
  limits.max_amplitude_evaluations = 4096;
  limits.max_bytes = (size_t)64 * 1024 * 1024;
  limits.max_order = order;
  return limits;
}

static MVMCKrylovStatus evaluate_bounded(
    const MVMCKrylovFockModel *model, const uint64_t *root,
    size_t word_count, const MVMCKrylovBoundedLimits *limits,
    TableAmplitude *amplitude, MVMCKrylovBoundedResult *result,
    uint64_t *plan_hash, size_t *workspace_bytes, size_t *cache_set_bytes,
    size_t *cache_allocated_bytes) {
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  MVMCKrylovStatus status =
      mvmc_bounded_krylov_plan_create(model, limits, &plan);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (plan_hash != NULL) *plan_hash = mvmc_bounded_krylov_plan_hash(plan);
  if (cache_set_bytes != NULL) {
    *cache_set_bytes = mvmc_bounded_krylov_cache_set_bytes(plan);
  }
  status = mvmc_bounded_krylov_workspace_create(plan, &workspace);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    if (workspace_bytes != NULL) {
      *workspace_bytes = mvmc_bounded_krylov_workspace_bytes(workspace);
    }
    status = mvmc_bounded_krylov_evaluate(
        workspace, root, word_count, scaled_table_callback, amplitude, result);
    if (cache_allocated_bytes != NULL) {
      *cache_allocated_bytes = result->statistics.cache_allocated_bytes;
    }
  }
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
  return status;
}

static MVMCKrylovStatus evaluate_reference(
    const MVMCKrylovFockModel *model, const uint64_t *root,
    size_t word_count, int order, TableAmplitude *amplitude,
    MVMCKrylovResult *result) {
  MVMCKrylovStatus status;
  const MVMCKrylovLimits limits = reference_limits(order);
  MVMCKrylovWorkspace *workspace =
      mvmc_krylov_workspace_create(model->site_count, &limits, &status);
  if (workspace == NULL) return status;
  status = mvmc_krylov_evaluate(workspace, model, root, word_count,
                                raw_table_callback, amplitude, result);
  mvmc_krylov_workspace_destroy(workspace);
  return status;
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

static int finite_nonnegative_timings(
    const MVMCKrylovBoundedStatistics *statistics) {
  int depth;
  const double values[] = {
      statistics->reset_seconds,
      statistics->total_seconds,
      statistics->evaluation_wall_seconds,
      statistics->amplitude_wall_seconds,
      statistics->connectivity_wall_seconds,
      statistics->cache_wall_seconds,
  };
  size_t index;
  for (index = 0; index < sizeof(values) / sizeof(values[0]); ++index) {
    if (!isfinite(values[index]) || values[index] < 0.0) return 0;
  }
  for (depth = 0; depth <= MVMC_KRYLOV_MAX_ORDER; ++depth) {
    if (!isfinite(statistics->depth_wall_seconds[depth]) ||
        statistics->depth_wall_seconds[depth] < 0.0) return 0;
  }
  return statistics->total_seconds + 1.0e-12 >=
         statistics->reset_seconds;
}

static void check_atomic_failure(const MVMCKrylovBoundedResult *result,
                                 MVMCKrylovStatus status,
                                 const char *label) {
  int order;
  CHECK(!result->valid && result->status == status,
        "%s did not return invalid status %d", label, (int)status);
  CHECK(result->evaluated_order == -1 &&
            result->statistics.completed_order == -1,
        "%s exposed a completed partial order", label);
  for (order = 0; order <= MVMC_KRYLOV_MAX_ORDER; ++order) {
    CHECK(result->value[order].state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
          "%s exposed partial value[%d]", label, order);
  }
}

static void test_dense_ed_and_p3(void) {
  const double complex a = 0.7 + 0.2 * I;
  const double complex b = -0.3 + 0.4 * I;
  const double diagonal[3] = {0.25, -0.5, 0.75};
  const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 2},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
      {MVMC_KRYLOV_FERMION_CREATE, 2},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
  };
  const MVMCKrylovHamiltonianTerm terms[] = {
      {diagonal[0], 0, 2, 0, 0}, {diagonal[1], 2, 2, 0, 1},
      {diagonal[2], 4, 2, 0, 2}, {a, 6, 2, 1, 0},
      {conj(a), 8, 2, 1, 1},     {b, 10, 2, 1, 2},
      {conj(b), 12, 2, 1, 3},
  };
  const MVMCKrylovFockModel model = electronic_model(
      3, 1, 0, terms, sizeof(terms) / sizeof(terms[0]), operators,
      sizeof(operators) / sizeof(operators[0]));
  const double complex matrix[3][3] = {
      {diagonal[0], a, 0.0},
      {conj(a), diagonal[1], b},
      {0.0, conj(b), diagonal[2]},
  };
  double complex expected[4][3];
  TableAmplitude fixture;
  int order;
  int row;
  memset(&fixture, 0, sizeof(fixture));
  fixture.trace_hash = UINT64_C(1469598103934665603);
  expected[0][0] = 0.2 - 0.1 * I;
  expected[0][1] = -0.3 + 0.8 * I;
  expected[0][2] = 0.6 + 0.4 * I;
  fixture.values[1] = expected[0][0];
  fixture.values[2] = expected[0][1];
  fixture.values[4] = expected[0][2];
  for (order = 1; order <= 3; ++order) {
    for (row = 0; row < 3; ++row) {
      int column;
      expected[order][row] = 0.0;
      for (column = 0; column < 3; ++column) {
        expected[order][row] +=
            matrix[row][column] * expected[order - 1][column];
      }
    }
  }
  for (row = 0; row < 3; ++row) {
    const uint64_t root = UINT64_C(1) << row;
    const MVMCKrylovBoundedLimits limits =
        bounded_limits(3, (size_t)1024 * 1024);
    MVMCKrylovBoundedResult bounded;
    MVMCKrylovResult reference;
    TableAmplitude bounded_amplitude = fixture;
    TableAmplitude reference_amplitude = fixture;
    CHECK(evaluate_bounded(&model, &root, 1, &limits, &bounded_amplitude,
                           &bounded, NULL, NULL, NULL, NULL) ==
              MVMC_KRYLOV_STATUS_OK &&
              bounded.valid,
          "dense ED bounded root %d failed", row);
    CHECK(evaluate_reference(&model, &root, 1, 3, &reference_amplitude,
                             &reference) == MVMC_KRYLOV_STATUS_OK &&
              reference.valid,
          "dense ED P3 root %d failed", row);
    for (order = 0; order <= 3; ++order) {
      double complex actual = NAN + I * NAN;
      CHECK(scaled_to_raw(&bounded.value[order], &actual),
            "dense root %d order %d did not export", row, order);
      CHECK(close_complex(actual, expected[order][row], 3.0e-13),
            "dense ED mismatch root=%d order=%d", row, order);
      CHECK(close_complex(actual, reference.value[order], 3.0e-13),
            "P3 A/B mismatch root=%d order=%d", row, order);
    }
    CHECK(bounded.statistics.engine_heap_allocations == 0,
          "evaluate hot path reported a heap allocation");
    CHECK(bounded.statistics.plan_bytes >= bounded.statistics.model_bytes &&
              bounded.statistics.workspace_bytes >=
                  bounded.statistics.frame_bytes +
                      bounded.statistics.scratch_bytes,
          "byte breakdown is inconsistent");
    CHECK(finite_nonnegative_timings(&bounded.statistics),
          "bounded timing fields are not finite/nonnegative");
  }
}

static void test_zero_bridge_and_repeated_reset(void) {
  const double complex hopping = 1.0 + 2.0 * I;
  const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
  };
  const MVMCKrylovHamiltonianTerm terms[] = {
      {hopping, 0, 2, 1, 0}, {conj(hopping), 2, 2, 1, 1}};
  const MVMCKrylovFockModel model =
      electronic_model(2, 1, 0, terms, 2, operators, 4);
  MVMCKrylovBoundedLimits limits = bounded_limits(3, 4096);
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  MVMCKrylovBoundedResult first;
  MVMCKrylovBoundedResult second;
  TableAmplitude amplitude;
  const uint64_t root = UINT64_C(1);
  double complex raw;
  memset(&amplitude, 0, sizeof(amplitude));
  amplitude.trace_hash = UINT64_C(1469598103934665603);
  amplitude.values[2] = 3.0 - I;
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_bounded_krylov_workspace_create(plan, &workspace) ==
                MVMC_KRYLOV_STATUS_OK,
        "zero-bridge workspace creation failed");
  if (workspace == NULL) {
    mvmc_bounded_krylov_plan_destroy(plan);
    return;
  }
  CHECK(mvmc_bounded_krylov_evaluate(workspace, &root, 1,
                                      scaled_table_callback, &amplitude,
                                      &first) == MVMC_KRYLOV_STATUS_OK,
        "zero-bridge evaluation failed");
  CHECK(first.value[0].state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
        "zero bridge v0 was not exact zero");
  CHECK(scaled_to_raw(&first.value[1], &raw) &&
            close_complex(raw, hopping * amplitude.values[2], 1.0e-14),
        "zero bridge was truncated at depth 1");
  CHECK(first.value[2].state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
        "zero bridge v2 was not exact zero");
  CHECK(scaled_to_raw(&first.value[3], &raw) &&
            close_complex(raw,
                          hopping * amplitude.values[2] *
                              creal(hopping * conj(hopping)),
                          2.0e-13),
        "zero bridge depth-3 value mismatch");
  amplitude.values[2] = -0.75 + 0.25 * I;
  amplitude.calls = 0;
  amplitude.trace_hash = UINT64_C(1469598103934665603);
  CHECK(mvmc_bounded_krylov_evaluate(workspace, &root, 1,
                                      scaled_table_callback, &amplitude,
                                      &second) == MVMC_KRYLOV_STATUS_OK,
        "repeated evaluation failed");
  CHECK(scaled_to_raw(&second.value[1], &raw) &&
            close_complex(raw, hopping * amplitude.values[2], 1.0e-14),
        "logical reset reused a stale root cache");
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
}

static void test_scaled_multipath_cancellation(void) {
  const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
      {MVMC_KRYLOV_FERMION_CREATE, 2},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
  };
  const MVMCKrylovHamiltonianTerm terms[] = {
      {1.0, 0, 2, 1, 0}, {1.0, 2, 2, 1, 1},
      {1.0, 4, 2, 1, 2}, {1.0, 6, 2, 1, 3},
  };
  const MVMCKrylovFockModel model =
      electronic_model(3, 1, 0, terms, 4, operators, 8);
  const MVMCKrylovBoundedLimits limits = bounded_limits(1, 0);
  const uint64_t root = UINT64_C(2);
  TableAmplitude amplitude;
  MVMCKrylovBoundedResult result;
  memset(&amplitude, 0, sizeof(amplitude));
  amplitude.values[1] = 1.0;
  amplitude.values[4] = -1.0;
  CHECK(evaluate_bounded(&model, &root, 1, &limits, &amplitude, &result,
                         NULL, NULL, NULL, NULL) ==
            MVMC_KRYLOV_STATUS_OK,
        "scaled multipath cancellation evaluation failed");
  CHECK(result.value[1].state == MVMC_SCALED_COMPLEX_NUMERIC_ZERO &&
            result.statistics.computed_numeric_zero_values > 0,
        "finite multipath cancellation was not diagnosed as NUMERIC_ZERO");
}

static void test_duplicate_permutation_and_cache_grid(void) {
  const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
  };
  MVMCKrylovHamiltonianTerm terms_a[] = {
      {1.0e16, 0, 2, 1, 10}, {1.0, 0, 2, 1, 12},
      {-1.0e16, 0, 2, 1, 11}, {1.0e16, 2, 2, 1, 10},
      {1.0, 2, 2, 1, 12}, {-1.0e16, 2, 2, 1, 11},
  };
  const MVMCKrylovHamiltonianTerm terms_b[] = {
      {-1.0e16, 2, 2, 1, 11}, {1.0, 0, 2, 1, 12},
      {1.0e16, 2, 2, 1, 10},  {-1.0e16, 0, 2, 1, 11},
      {1.0, 2, 2, 1, 12},     {1.0e16, 0, 2, 1, 10},
  };
  MVMCKrylovFockModel model_a =
      electronic_model(2, 1, 0, terms_a, 6, operators, 4);
  const MVMCKrylovFockModel model_b =
      electronic_model(2, 1, 0, terms_b, 6, operators, 4);
  const uint64_t root = UINT64_C(1);
  TableAmplitude fixture;
  MVMCKrylovBoundedLimits hash_limits = bounded_limits(3, 0);
  MVMCKrylovBoundedPlan *plan_a = NULL;
  MVMCKrylovBoundedPlan *plan_b = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  MVMCKrylovBoundedResult deep_copy_result;
  uint64_t canonical_hash = 0;
  size_t set_bytes = 0;
  size_t capacity_index;
  size_t capacities[4];
  MVMCScaledComplex baseline[MVMC_KRYLOV_MAX_ORDER + 1];
  int have_baseline = 0;
  memset(&fixture, 0, sizeof(fixture));
  fixture.trace_hash = UINT64_C(1469598103934665603);
  fixture.values[1] = 3.0;
  fixture.values[2] = -2.0;
  CHECK(((1.0 + 1.0e16) + -1.0e16) == 0.0,
        "dynamic-range fixture no longer defeats naive summation");
  CHECK(mvmc_bounded_krylov_plan_create(&model_a, &hash_limits, &plan_a) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_bounded_krylov_plan_create(&model_b, &hash_limits, &plan_b) ==
                MVMC_KRYLOV_STATUS_OK,
        "permutation plans failed");
  if (plan_a != NULL && plan_b != NULL) {
    canonical_hash = mvmc_bounded_krylov_plan_hash(plan_a);
    CHECK(canonical_hash == mvmc_bounded_krylov_plan_hash(plan_b),
          "term permutation changed canonical plan hash");
    set_bytes = mvmc_bounded_krylov_cache_set_bytes(plan_a);
  }
  CHECK(set_bytes > 0, "cache set byte formula is zero");
  CHECK(mvmc_bounded_krylov_workspace_create(plan_a, &workspace) ==
            MVMC_KRYLOV_STATUS_OK,
        "deep-copy workspace creation failed");
  terms_a[0].coefficient = 17.0;
  model_a.hermitian = 0;
  if (workspace != NULL) {
    TableAmplitude amplitude = fixture;
    CHECK(mvmc_bounded_krylov_evaluate(workspace, &root, 1,
                                        scaled_table_callback, &amplitude,
                                        &deep_copy_result) ==
              MVMC_KRYLOV_STATUS_OK,
          "immutable plan did not survive source mutation");
    {
      double complex raw;
      CHECK(scaled_to_raw(&deep_copy_result.value[1], &raw) &&
                close_complex(raw, -2.0, 1.0e-14),
            "Neumaier duplicate residual was lost");
    }
  }
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan_b);
  mvmc_bounded_krylov_plan_destroy(plan_a);

  terms_a[0].coefficient = 1.0e16;
  model_a.hermitian = 1;
  capacities[0] = 0;
  capacities[1] = set_bytes + set_bytes / 2;
  capacities[2] = (size_t)64 * 1024 * 1024;
  capacities[3] = (size_t)256 * 1024 * 1024;
  for (capacity_index = 0;
       capacity_index < (world_size == 1 ? 4U : 2U); ++capacity_index) {
    MVMCKrylovBoundedLimits limits =
        bounded_limits(3, capacities[capacity_index]);
    MVMCKrylovBoundedResult result;
    TableAmplitude amplitude = fixture;
    size_t actual_cache_bytes = 0;
    CHECK(evaluate_bounded(&model_a, &root, 1, &limits, &amplitude, &result,
                           NULL, NULL, NULL, &actual_cache_bytes) ==
              MVMC_KRYLOV_STATUS_OK,
          "cache grid capacity %zu failed", capacities[capacity_index]);
    if (!have_baseline) {
      memcpy(baseline, result.value, sizeof(baseline));
      have_baseline = 1;
    } else {
      int order;
      for (order = 0; order <= 3; ++order) {
        CHECK(scaled_values_bitwise_equal(&baseline[order],
                                          &result.value[order]),
              "cache capacity changed scaled-value fields at order %d",
              order);
      }
    }
    CHECK(actual_cache_bytes <= capacities[capacity_index] &&
              (set_bytes == 0 || actual_cache_bytes % set_bytes == 0),
          "cache request was not truncated to complete sets");
    if (capacity_index == 0) {
      CHECK(actual_cache_bytes == 0 && result.statistics.cache_insertions == 0 &&
                result.statistics.cache_hits[0] == 0,
            "zero-capacity cache was not disabled");
    }
    if (capacity_index == 1) {
      CHECK(actual_cache_bytes == set_bytes &&
                result.statistics.cache_evictions > 0,
            "one-set collision/eviction fixture did not evict");
    }
  }
  {
    MVMCKrylovBoundedLimits changed = hash_limits;
    MVMCKrylovBoundedPlan *changed_plan = NULL;
    changed.amplitude_policy_hash ^= UINT64_C(1);
    CHECK(mvmc_bounded_krylov_plan_create(&model_a, &changed,
                                           &changed_plan) ==
              MVMC_KRYLOV_STATUS_OK &&
              mvmc_bounded_krylov_plan_hash(changed_plan) != canonical_hash,
          "amplitude-policy change did not change plan hash");
    mvmc_bounded_krylov_plan_destroy(changed_plan);
  }
}

static void test_fermion_sign_pure_spin_and_multiword(void) {
  {
    const MVMCKrylovFermionOperator operators[] = {
        {MVMC_KRYLOV_FERMION_CREATE, 0},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
        {MVMC_KRYLOV_FERMION_CREATE, 1},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
    };
    const MVMCKrylovHamiltonianTerm terms[] = {
        {1.0, 0, 2, 1, 0}, {1.0, 2, 2, 1, 1}};
    const MVMCKrylovFockModel model =
        electronic_model(2, 1, 0, terms, 2, operators, 4);
    const MVMCKrylovBoundedLimits limits = bounded_limits(1, 0);
    const uint64_t root = UINT64_C(1);
    TableAmplitude amplitude;
    MVMCKrylovBoundedResult result;
    double complex raw;
    memset(&amplitude, 0, sizeof(amplitude));
    amplitude.values[2] = -2.0;
    CHECK(evaluate_bounded(&model, &root, 1, &limits, &amplitude, &result,
                           NULL, NULL, NULL, NULL) ==
              MVMC_KRYLOV_STATUS_OK &&
              scaled_to_raw(&result.value[1], &raw) &&
              close_complex(raw, -2.0, 1.0e-14),
          "electronic fermion sign/orientation mismatch");
  }
  {
    const MVMCKrylovFermionOperator operators[] = {
        {MVMC_KRYLOV_FERMION_CREATE, 0},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
        {MVMC_KRYLOV_FERMION_CREATE, 3},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
        {MVMC_KRYLOV_FERMION_CREATE, 1},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
        {MVMC_KRYLOV_FERMION_CREATE, 2},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 3},
        {MVMC_KRYLOV_FERMION_CREATE, 1},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
        {MVMC_KRYLOV_FERMION_CREATE, 2},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
    };
    const MVMCKrylovHamiltonianTerm terms[] = {
        {0.75, 0, 4, 2, 0}, {0.75, 4, 4, 2, 1},
        {-0.25, 8, 4, 3, 0}};
    MVMCKrylovFockModel model =
        electronic_model(2, 1, 1, terms, 3, operators, 12);
    const MVMCKrylovBoundedLimits limits = bounded_limits(1, 0);
    const uint64_t root = UINT64_C(6);
    TableAmplitude amplitude;
    MVMCKrylovBoundedResult result;
    double complex raw;
    model.pure_spin = 1;
    memset(&amplitude, 0, sizeof(amplitude));
    amplitude.values[6] = 2.0;
    amplitude.values[9] = 2.0;
    CHECK(evaluate_bounded(&model, &root, 1, &limits, &amplitude, &result,
                           NULL, NULL, NULL, NULL) ==
              MVMC_KRYLOV_STATUS_OK &&
              scaled_to_raw(&result.value[1], &raw) &&
              close_complex(raw, 1.0, 1.0e-14),
          "pure-spin Exchange/XXZ sign or diagonal mismatch");
  }
  {
    const MVMCKrylovFockModel model =
        electronic_model(33, 1, 1, NULL, 0, NULL, 0);
    const MVMCKrylovBoundedLimits limits = bounded_limits(0, 0);
    uint64_t root[2] = {UINT64_C(1), UINT64_C(2)};
    TableAmplitude amplitude;
    MVMCKrylovBoundedResult result;
    memset(&amplitude, 0, sizeof(amplitude));
    amplitude.values[63] = 1.0;
    CHECK(evaluate_bounded(&model, root, 2, &limits, &amplitude, &result,
                           NULL, NULL, NULL, NULL) ==
              MVMC_KRYLOV_STATUS_OK,
          "valid multiword configuration failed");
    root[1] |= UINT64_C(4);
    CHECK(evaluate_bounded(&model, root, 2, &limits, &amplitude, &result,
                           NULL, NULL, NULL, NULL) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "multiword padding bit was accepted");
  }
}

static void test_limits_wrap_and_atomic_failures(void) {
  const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
  };
  const MVMCKrylovHamiltonianTerm terms[] = {
      {1.0, 0, 2, 1, 0}, {1.0, 2, 2, 1, 1}};
  const MVMCKrylovFockModel model =
      electronic_model(2, 1, 0, terms, 2, operators, 4);
  const uint64_t root = UINT64_C(1);
  TableAmplitude fixture;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  MVMCKrylovBoundedResult result;
  MVMCKrylovBoundedLimits limits;
  size_t set_bytes;
  memset(&fixture, 0, sizeof(fixture));
  fixture.values[1] = 1.0;
  fixture.values[2] = 2.0;

  limits = bounded_limits(-1, 0);
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            plan == NULL,
        "negative max order was accepted");
  limits = bounded_limits(MVMC_KRYLOV_MAX_ORDER + 1, 0);
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            plan == NULL,
        "max order above the supported domain was accepted");
  limits = bounded_limits(1, 0);
  limits.max_node_expansions = 0;
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            plan == NULL,
        "zero node limit was accepted");
  {
    MVMCKrylovFockModel invalid_model = model;
    invalid_model.hermitian = 0;
    limits = bounded_limits(1, 0);
    CHECK(mvmc_bounded_krylov_plan_create(&invalid_model, &limits, &plan) ==
              MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL &&
              plan == NULL,
          "non-Hermitian model was accepted");
  }
  {
    MVMCKrylovHamiltonianTerm invalid_terms[2];
    MVMCKrylovFockModel invalid_model = model;
    memcpy(invalid_terms, terms, sizeof(invalid_terms));
    invalid_terms[0].coefficient = NAN + 0.0 * I;
    invalid_model.terms = invalid_terms;
    limits = bounded_limits(1, 0);
    CHECK(mvmc_bounded_krylov_plan_create(&invalid_model, &limits, &plan) ==
              MVMC_KRYLOV_STATUS_NONFINITE &&
              plan == NULL,
          "nonfinite model coefficient was accepted");
  }
  {
    MVMCKrylovFermionOperator invalid_operators[4];
    MVMCKrylovFockModel invalid_model = model;
    memcpy(invalid_operators, operators, sizeof(invalid_operators));
    invalid_operators[0].orbital = 4;
    invalid_model.operators = invalid_operators;
    limits = bounded_limits(1, 0);
    CHECK(mvmc_bounded_krylov_plan_create(&invalid_model, &limits, &plan) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              plan == NULL,
          "out-of-range model orbital was accepted");
  }

  limits = bounded_limits(2, 0);
  limits.max_row_transitions = 1;
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
            plan == NULL,
        "insufficient max-row limit was accepted");

  limits = bounded_limits(2, 0);
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_OK,
        "limit fixture plan creation failed");
  if (plan == NULL) return;
  set_bytes = mvmc_bounded_krylov_cache_set_bytes(plan);
  mvmc_bounded_krylov_plan_destroy(plan);
  plan = NULL;

  limits = bounded_limits(3, SIZE_MAX);
  limits.max_workspace_bytes = SIZE_MAX;
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_bounded_krylov_workspace_create(plan, &workspace) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            workspace == NULL,
        "overflowing cache capacity was accepted");
  mvmc_bounded_krylov_plan_destroy(plan);
  plan = NULL;

  limits = bounded_limits(2, 0);
  limits.max_workspace_bytes = 1;
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_bounded_krylov_workspace_create(plan, &workspace) ==
                MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
            workspace == NULL,
        "workspace byte limit did not fire");
  mvmc_bounded_krylov_plan_destroy(plan);
  plan = NULL;

  limits = bounded_limits(2, 0);
  limits.max_node_expansions = 1;
  {
    TableAmplitude amplitude = fixture;
    CHECK(evaluate_bounded(&model, &root, 1, &limits, &amplitude, &result,
                           NULL, NULL, NULL, NULL) ==
              MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
          "node expansion limit did not fire");
    check_atomic_failure(&result, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                         "node limit");
  }

  limits = bounded_limits(1, 0);
  limits.max_terminal_amplitude_calls = 1;
  {
    TableAmplitude amplitude = fixture;
    CHECK(evaluate_bounded(&model, &root, 1, &limits, &amplitude, &result,
                           NULL, NULL, NULL, NULL) ==
              MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
          "terminal call limit did not fire");
    check_atomic_failure(&result, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                         "terminal limit");
  }

  limits = bounded_limits(2, 0);
  limits.max_total_row_transitions = 1;
  {
    TableAmplitude amplitude = fixture;
    CHECK(evaluate_bounded(&model, &root, 1, &limits, &amplitude, &result,
                           NULL, NULL, NULL, NULL) ==
              MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
          "total row-transition limit did not fire");
    check_atomic_failure(&result, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                         "row-transition limit");
  }

  limits = bounded_limits(0, 0);
  {
    TableAmplitude amplitude = fixture;
    amplitude.forced_status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    CHECK(evaluate_bounded(&model, &root, 1, &limits, &amplitude, &result,
                           NULL, NULL, NULL, NULL) ==
              MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE,
          "callback failure was not propagated");
    check_atomic_failure(&result, MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE,
                         "callback failure");
  }
  {
    TableAmplitude amplitude = fixture;
    amplitude.return_nonfinite = 1;
    CHECK(evaluate_bounded(&model, &root, 1, &limits, &amplitude, &result,
                           NULL, NULL, NULL, NULL) ==
              MVMC_KRYLOV_STATUS_NONFINITE,
          "nonfinite callback value was accepted");
    check_atomic_failure(&result, MVMC_KRYLOV_STATUS_NONFINITE,
                         "nonfinite callback");
  }

  limits = bounded_limits(0, set_bytes);
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_bounded_krylov_workspace_create(plan, &workspace) ==
                MVMC_KRYLOV_STATUS_OK,
        "cache wrap fixture creation failed");
  if (workspace != NULL) {
    TableAmplitude amplitude = fixture;
    CHECK(mvmc_bounded_krylov_testing_force_cache_counters(
              workspace, UINT64_MAX, 0) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_bounded_krylov_evaluate(
                  workspace, &root, 1, scaled_table_callback, &amplitude,
                  &result) == MVMC_KRYLOV_STATUS_OK &&
              result.statistics.cache_epoch_full_clears == 1,
          "epoch wrap did not perform a full logical clear");
    amplitude = fixture;
    CHECK(mvmc_bounded_krylov_testing_force_cache_counters(
              workspace, 7, UINT64_MAX) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_bounded_krylov_evaluate(
                  workspace, &root, 1, scaled_table_callback, &amplitude,
                  &result) == MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
          "access-counter overflow did not fail loud");
    check_atomic_failure(&result, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                         "cache access overflow");
  }
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
}

static void test_rank_invariance(void) {
#ifdef _mpi_use
  const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
  };
  const MVMCKrylovHamiltonianTerm terms[] = {
      {0.5 + 0.25 * I, 0, 2, 1, 0},
      {0.5 - 0.25 * I, 2, 2, 1, 1},
  };
  const MVMCKrylovFockModel model =
      electronic_model(2, 1, 0, terms, 2, operators, 4);
  const MVMCKrylovBoundedLimits limits = bounded_limits(3, 4096);
  const uint64_t root = UINT64_C(1);
  TableAmplitude amplitude;
  MVMCKrylovBoundedResult result;
  uint64_t plan_hash = 0;
  uint64_t fields[11];
  uint64_t minimum[11];
  uint64_t maximum[11];
  int order;
  memset(&amplitude, 0, sizeof(amplitude));
  amplitude.trace_hash = UINT64_C(1469598103934665603);
  amplitude.values[1] = 0.25 + 0.5 * I;
  amplitude.values[2] = -0.75 + 0.125 * I;
  CHECK(evaluate_bounded(&model, &root, 1, &limits, &amplitude, &result,
                         &plan_hash, NULL, NULL, NULL) ==
            MVMC_KRYLOV_STATUS_OK,
        "rank-invariance evaluation failed");
  fields[0] = plan_hash;
  fields[1] = result.statistics.trace_hash;
  fields[2] = amplitude.trace_hash;
  fields[3] = amplitude.calls;
  fields[4] = result.statistics.node_expansions;
  fields[5] = result.statistics.raw_row_transitions;
  fields[6] = result.statistics.cache_insertions;
  fields[7] = result.statistics.cache_evictions;
  fields[8] = result.statistics.row_transition_peak;
  fields[9] = result.statistics.terminal_amplitude_calls;
  fields[10] = (uint64_t)(unsigned int)result.status;
  MPI_Allreduce(fields, minimum, 11, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(fields, maximum, 11, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  CHECK(memcmp(minimum, maximum, sizeof(minimum)) == 0,
        "plan/ROW/cache/callback trace differs by MPI rank");
  for (order = 0; order <= 3; ++order) {
    uint64_t local[8];
    uint64_t gathered_min[8];
    uint64_t gathered_max[8];
    scaled_value_field_bits(&result.value[order], local);
    MPI_Allreduce(local, gathered_min, 8, MPI_UINT64_T, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(local, gathered_max, 8, MPI_UINT64_T, MPI_MAX,
                  MPI_COMM_WORLD);
    CHECK(memcmp(gathered_min, gathered_max, sizeof(gathered_min)) == 0,
          "scaled-value fields differ by rank at order %d", order);
  }
#endif
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#else
  (void)argc;
  (void)argv;
#endif
  test_dense_ed_and_p3();
  test_zero_bridge_and_repeated_reset();
  test_scaled_multipath_cancellation();
  test_duplicate_permutation_and_cache_grid();
  test_fermion_sign_pure_spin_and_multiword();
  test_limits_wrap_and_atomic_failures();
  test_rank_invariance();

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
      fprintf(stderr, "BoundedKrylovEngine_Unit: %d failure(s)\n", failures);
    }
#ifdef _mpi_use
    MPI_Finalize();
#endif
    return 1;
  }
  if (world_rank == 0) puts("BoundedKrylovEngine_Unit: PASS");
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return 0;
}
