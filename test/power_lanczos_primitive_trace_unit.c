#include "power_lanczos_primitive_trace.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

#define CHECK(condition, ...)                                                \
  do {                                                                       \
    if (!(condition)) {                                                      \
      fprintf(stderr, "PowerLanczosPrimitiveTrace_Unit FAIL: ");           \
      fprintf(stderr, __VA_ARGS__);                                          \
      fprintf(stderr, " (line %d)\n", __LINE__);                           \
      ++failures;                                                            \
    }                                                                        \
  } while (0)

static MVMCPowerLanczosPrimitiveTrace *create_trace(size_t primitives,
                                                     size_t groups,
                                                     size_t samples) {
  MVMCPowerLanczosPrimitiveTrace *trace = NULL;
  CHECK(mvmc_power_lanczos_primitive_trace_create(
            MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_COEFFICIENT,
            primitives, groups, samples, 1024 * 1024, &trace) ==
                MVMC_KRYLOV_STATUS_OK &&
            trace != NULL,
        "trace create");
  return trace;
}

static void test_complete_trace(void) {
  MVMCPowerLanczosPrimitiveTrace *trace = create_trace(2, 2, 3);
  MVMCPowerLanczosPrimitiveTraceSummary summary;
  MVMCPowerLanczosPrimitiveTraceGroup group_record;
  double scalars[24];
  size_t group;
  size_t sample;
  CHECK(mvmc_power_lanczos_primitive_trace_allocated_bytes(trace) > 0,
        "allocated byte getter before freeze");
  CHECK(mvmc_power_lanczos_primitive_trace_group(trace, 0,
                                                  &group_record) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !group_record.valid,
        "getter before freeze");
  CHECK(mvmc_power_lanczos_primitive_trace_register_group(trace, 0, 0, 2) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_primitive_trace_register_group(trace, 1, 2,
                                                               1) ==
                MVMC_KRYLOV_STATUS_OK,
        "canonical group registration");
  for (sample = 0; sample < 3; ++sample) {
    for (group = 0; group < 2; ++group) {
      double complex values[2];
      double bounds[2] = {0.0, 0.01};
      uint8_t support[2] = {
          MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO,
          MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO |
              MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO};
      values[0] = (double)(10 * group + sample) +
                  I * (double)(100 + 10 * group + sample);
      values[1] = (double)(20 + 10 * group + sample) -
                  I * (double)(30 + 10 * group + sample);
      if (group == 1 && sample == 2) {
        values[1] = 0.0;
        bounds[1] = 0.0;
        support[1] = MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO;
      }
      CHECK(mvmc_power_lanczos_primitive_trace_append(
                trace, group, sample, values, bounds, support, 2,
                (double)(10 * group + sample)) == MVMC_KRYLOV_STATUS_OK,
            "append group %zu sample %zu", group, sample);
    }
  }
  CHECK(mvmc_power_lanczos_primitive_trace_freeze(trace) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_primitive_trace_summary(trace, &summary) ==
                MVMC_KRYLOV_STATUS_OK &&
            summary.valid && summary.frozen && summary.version ==
                MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_VERSION &&
            summary.stage ==
                MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_COEFFICIENT &&
            summary.primitive_count == 2 && summary.scalar_count == 4 &&
            summary.group_count == 2 && summary.samples_per_group == 3 &&
            summary.appended_sample_count == 6 &&
            summary.allocated_bytes > sizeof(summary),
        "freeze/summary");
  for (group = 0; group < 2; ++group) {
    CHECK(mvmc_power_lanczos_primitive_trace_group(trace, group,
                                                    &group_record) ==
                  MVMC_KRYLOV_STATUS_OK &&
              group_record.valid && group_record.group_ordinal == group &&
              group_record.leader_world_rank == 2 * group &&
              group_record.chain_size == (group == 0 ? 2 : 1) &&
              group_record.sample_count == 3,
          "group getter %zu", group);
  }
  {
    double complex values[2];
    double bounds[2];
    uint8_t support[2];
    double tail = -1.0;
    CHECK(mvmc_power_lanczos_primitive_trace_sample(
              trace, 1, 2, values, 2, bounds, 2, support, 2, &tail) ==
                  MVMC_KRYLOV_STATUS_OK &&
              creal(values[0]) == 12.0 && cimag(values[0]) == 112.0 &&
              values[1] == 0.0 && bounds[1] == 0.0 &&
              support[1] ==
                  MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO &&
              tail == 12.0,
          "sample getter");
    values[0] = 1.0;
    bounds[0] = 1.0;
    support[0] = 255;
    tail = 1.0;
    CHECK(mvmc_power_lanczos_primitive_trace_sample(
              trace, 0, 0, values, 1, bounds, 1, support, 1, &tail) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              values[0] == 0.0 && bounds[0] == 0.0 && support[0] == 0 &&
              tail == 0.0,
          "undersized getter zeroes outputs");
  }
  CHECK(mvmc_power_lanczos_primitive_trace_export_scalars(
            trace, scalars, sizeof(scalars) / sizeof(scalars[0])) ==
                MVMC_KRYLOV_STATUS_OK &&
            scalars[0] == 0.0 && scalars[3] == 10.0 &&
            scalars[6] == 100.0 && scalars[9] == 110.0 &&
            scalars[12] == 20.0 && scalars[15] == 30.0 &&
            scalars[18] == -30.0 && scalars[21] == -40.0 &&
            scalars[23] == 0.0,
        "scalar-major export");
  {
    const double complex value = 1.0;
    const double bound = 0.0;
    const uint8_t support =
        MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO;
    CHECK(mvmc_power_lanczos_primitive_trace_append(
              trace, 0, 0, &value, &bound, &support, 1, 0.0) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              mvmc_power_lanczos_primitive_trace_summary(trace, &summary) ==
                  MVMC_KRYLOV_STATUS_OK &&
              summary.valid && summary.frozen,
          "frozen trace immutable");
  }
  mvmc_power_lanczos_primitive_trace_destroy(trace);
}

static void test_failure_boundaries(void) {
  MVMCPowerLanczosPrimitiveTrace *trace = NULL;
  MVMCPowerLanczosPrimitiveTraceSummary summary;
  CHECK(mvmc_power_lanczos_primitive_trace_create(
            MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_FINAL, 1, 1, 1, 1,
            &trace) == MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
            trace == NULL,
        "allocation cap");
  CHECK(mvmc_power_lanczos_primitive_trace_create(
            MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_FINAL, SIZE_MAX, 2, 2,
            SIZE_MAX, &trace) == MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
            trace == NULL,
        "shape overflow");

  trace = create_trace(1, 1, 1);
  CHECK(mvmc_power_lanczos_primitive_trace_register_group(trace, 0, 0, 1) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_primitive_trace_freeze(trace) ==
                MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE &&
            mvmc_power_lanczos_primitive_trace_summary(trace, &summary) ==
                MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE &&
            !summary.valid && !summary.frozen,
        "incomplete trace fail-closed");
  mvmc_power_lanczos_primitive_trace_destroy(trace);

  trace = create_trace(1, 2, 1);
  CHECK(mvmc_power_lanczos_primitive_trace_register_group(trace, 0, 0, 1) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_primitive_trace_register_group(trace, 1, 0,
                                                               1) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "duplicate leader rejected");
  mvmc_power_lanczos_primitive_trace_destroy(trace);

  trace = create_trace(1, 1, 2);
  CHECK(mvmc_power_lanczos_primitive_trace_register_group(trace, 0, 0, 1) ==
                MVMC_KRYLOV_STATUS_OK,
        "append-order group");
  {
    const double complex value = 1.0;
    const double bound = 0.0;
    const uint8_t support =
        MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO;
    CHECK(mvmc_power_lanczos_primitive_trace_append(
              trace, 0, 1, &value, &bound, &support, 1, 0.0) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "nonsequential append rejected");
  }
  mvmc_power_lanczos_primitive_trace_destroy(trace);
}

static void test_sample_mutations(void) {
  size_t mutation;
  for (mutation = 0; mutation < 7; ++mutation) {
    MVMCPowerLanczosPrimitiveTrace *trace = create_trace(1, 1, 1);
    double complex value = 1.0;
    double bound = 0.0;
    uint8_t support =
        MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO;
    double tail = 0.0;
    MVMCKrylovStatus expected = MVMC_KRYLOV_STATUS_NONFINITE;
    CHECK(mvmc_power_lanczos_primitive_trace_register_group(trace, 0, 0, 1) ==
                  MVMC_KRYLOV_STATUS_OK,
          "mutation group %zu", mutation);
    switch (mutation) {
      case 0:
        value = NAN;
        break;
      case 1:
        bound = -1.0;
        break;
      case 2:
        support = 0;
        break;
      case 3:
        support = UINT8_C(128);
        break;
      case 4:
        support = MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO;
        value = 0.0;
        break;
      case 5:
        support = MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO;
        break;
      case 6:
        tail = INFINITY;
        expected = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
        break;
    }
    CHECK(mvmc_power_lanczos_primitive_trace_append(
              trace, 0, 0, &value, &bound, &support, 1, tail) == expected,
          "sample mutation %zu", mutation);
    mvmc_power_lanczos_primitive_trace_destroy(trace);
  }
}

int main(void) {
  test_complete_trace();
  test_failure_boundaries();
  test_sample_mutations();
  if (failures != 0) {
    fprintf(stderr, "%d primitive trace checks failed\n", failures);
    return 1;
  }
  puts("power-Lanczos primitive trace checks passed");
  return 0;
}
