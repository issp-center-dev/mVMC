#ifndef MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_H
#define MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "bounded_krylov_engine.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_VERSION UINT64_C(1)

typedef enum {
  MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_COEFFICIENT = 1,
  MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_FINAL = 2
} MVMCPowerLanczosPrimitiveTraceStage;

typedef enum {
  MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO = 1u << 0,
  MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO = 1u << 1,
  MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO = 1u << 2
} MVMCPowerLanczosPrimitiveSupportFlag;

typedef struct MVMCPowerLanczosPrimitiveTrace
    MVMCPowerLanczosPrimitiveTrace;

typedef struct {
  int valid;
  MVMCPowerLanczosPrimitiveTraceStage stage;
  size_t group_ordinal;
  size_t leader_world_rank;
  size_t chain_size;
  size_t sample_count;
} MVMCPowerLanczosPrimitiveTraceGroup;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  uint64_t version;
  MVMCPowerLanczosPrimitiveTraceStage stage;
  int frozen;
  size_t primitive_count;
  size_t scalar_count;
  size_t group_count;
  size_t samples_per_group;
  size_t appended_sample_count;
  size_t allocated_bytes;
} MVMCPowerLanczosPrimitiveTraceSummary;

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_create(
    MVMCPowerLanczosPrimitiveTraceStage stage, size_t primitive_count,
    size_t group_count, size_t samples_per_group,
    size_t maximum_allocated_bytes,
    MVMCPowerLanczosPrimitiveTrace **trace);

void mvmc_power_lanczos_primitive_trace_destroy(
    MVMCPowerLanczosPrimitiveTrace *trace);

size_t mvmc_power_lanczos_primitive_trace_allocated_bytes(
    const MVMCPowerLanczosPrimitiveTrace *trace);

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_register_group(
    MVMCPowerLanczosPrimitiveTrace *trace, size_t group_ordinal,
    size_t leader_world_rank, size_t chain_size);

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_append(
    MVMCPowerLanczosPrimitiveTrace *trace, size_t group_ordinal,
    size_t sample_index, const double complex *values,
    const double *absolute_numeric_bounds, const uint8_t *support_flags,
    size_t primitive_count, double tail_magnitude);

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_freeze(
    MVMCPowerLanczosPrimitiveTrace *trace);

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_summary(
    const MVMCPowerLanczosPrimitiveTrace *trace,
    MVMCPowerLanczosPrimitiveTraceSummary *summary);

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_group(
    const MVMCPowerLanczosPrimitiveTrace *trace, size_t group_ordinal,
    MVMCPowerLanczosPrimitiveTraceGroup *group);

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_sample(
    const MVMCPowerLanczosPrimitiveTrace *trace, size_t group_ordinal,
    size_t sample_index, double complex *values, size_t value_capacity,
    double *absolute_numeric_bounds, size_t bound_capacity,
    uint8_t *support_flags, size_t support_capacity,
    double *tail_magnitude);

/*
 * Export scalar-major, group-major, sample-major real/imag coordinates.
 * Scalar 2*p is real(values[p]); scalar 2*p+1 is imag(values[p]).
 */
MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_export_scalars(
    const MVMCPowerLanczosPrimitiveTrace *trace, double *scalars,
    size_t scalar_capacity);

#endif /* bounded engine */

#endif /* MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_H */
