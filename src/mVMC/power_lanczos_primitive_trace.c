#include "power_lanczos_primitive_trace.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "power_lanczos_primitive_trace.c requires bounded Krylov"
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>

enum {
  TRACE_BUILDING = 1,
  TRACE_FROZEN = 2,
  TRACE_FAILED = 3,
  SUPPORT_MASK = MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO |
                 MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO |
                 MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO
};

struct MVMCPowerLanczosPrimitiveTrace {
  int valid;
  int state;
  MVMCKrylovStatus status;
  MVMCPowerLanczosPrimitiveTraceStage stage;
  size_t primitive_count;
  size_t scalar_count;
  size_t group_count;
  size_t samples_per_group;
  size_t sample_slot_count;
  size_t value_count;
  size_t appended_sample_count;
  size_t allocated_bytes;
  MVMCPowerLanczosPrimitiveTraceGroup *groups;
  size_t *next_sample;
  double complex *values;
  double *bounds;
  uint8_t *support;
  double *tails;
};

static int checked_add(size_t left, size_t right, size_t *result) {
  if (result == NULL || left > SIZE_MAX - right) return 0;
  *result = left + right;
  return 1;
}

static int checked_multiply(size_t left, size_t right, size_t *result) {
  if (result == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *result = left * right;
  return 1;
}

static int reserve(size_t *offset, size_t count, size_t item_size,
                   size_t alignment, size_t *start) {
  size_t aligned;
  size_t bytes;
  if (offset == NULL || start == NULL || alignment == 0 ||
      (alignment & (alignment - 1)) != 0 ||
      *offset > SIZE_MAX - (alignment - 1)) {
    return 0;
  }
  aligned = (*offset + alignment - 1) & ~(alignment - 1);
  if (!checked_multiply(count, item_size, &bytes) ||
      !checked_add(aligned, bytes, offset)) {
    return 0;
  }
  *start = aligned;
  return 1;
}

static int valid_stage(MVMCPowerLanczosPrimitiveTraceStage stage) {
  return stage == MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_COEFFICIENT ||
         stage == MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_FINAL;
}

static MVMCKrylovStatus fail_trace(MVMCPowerLanczosPrimitiveTrace *trace,
                                   MVMCKrylovStatus status) {
  if (trace != NULL && trace->state == TRACE_BUILDING) {
    trace->valid = 0;
    trace->state = TRACE_FAILED;
    trace->status = status;
  }
  return status;
}

static size_t value_offset(const MVMCPowerLanczosPrimitiveTrace *trace,
                           size_t primitive, size_t group, size_t sample) {
  return (primitive * trace->group_count + group) *
             trace->samples_per_group +
         sample;
}

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_create(
    MVMCPowerLanczosPrimitiveTraceStage stage, size_t primitive_count,
    size_t group_count, size_t samples_per_group,
    size_t maximum_allocated_bytes,
    MVMCPowerLanczosPrimitiveTrace **trace) {
  MVMCPowerLanczosPrimitiveTrace *created;
  size_t sample_slots;
  size_t value_count;
  size_t scalar_count;
  size_t offset = sizeof(*created);
  size_t group_offset;
  size_t next_offset;
  size_t value_offset_bytes;
  size_t bound_offset;
  size_t support_offset;
  size_t tail_offset;
  if (trace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *trace = NULL;
  if (!valid_stage(stage) || primitive_count == 0 || group_count == 0 ||
      samples_per_group == 0 || maximum_allocated_bytes == 0 ||
      !checked_multiply(group_count, samples_per_group, &sample_slots) ||
      !checked_multiply(primitive_count, sample_slots, &value_count) ||
      !checked_multiply(primitive_count, 2, &scalar_count) ||
      !reserve(&offset, group_count, sizeof(MVMCPowerLanczosPrimitiveTraceGroup),
               _Alignof(MVMCPowerLanczosPrimitiveTraceGroup), &group_offset) ||
      !reserve(&offset, group_count, sizeof(size_t), _Alignof(size_t),
               &next_offset) ||
      !reserve(&offset, value_count, sizeof(double complex),
               _Alignof(double complex), &value_offset_bytes) ||
      !reserve(&offset, value_count, sizeof(double), _Alignof(double),
               &bound_offset) ||
      !reserve(&offset, value_count, sizeof(uint8_t), _Alignof(uint8_t),
               &support_offset) ||
      !reserve(&offset, sample_slots, sizeof(double), _Alignof(double),
               &tail_offset)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  if (offset > maximum_allocated_bytes) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  created = (MVMCPowerLanczosPrimitiveTrace *)calloc(1, offset);
  if (created == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  created->valid = 1;
  created->state = TRACE_BUILDING;
  created->status = MVMC_KRYLOV_STATUS_OK;
  created->stage = stage;
  created->primitive_count = primitive_count;
  created->scalar_count = scalar_count;
  created->group_count = group_count;
  created->samples_per_group = samples_per_group;
  created->sample_slot_count = sample_slots;
  created->value_count = value_count;
  created->allocated_bytes = offset;
  created->groups = (MVMCPowerLanczosPrimitiveTraceGroup *)(
      (unsigned char *)created + group_offset);
  created->next_sample =
      (size_t *)((unsigned char *)created + next_offset);
  created->values =
      (double complex *)((unsigned char *)created + value_offset_bytes);
  created->bounds = (double *)((unsigned char *)created + bound_offset);
  created->support = (uint8_t *)((unsigned char *)created + support_offset);
  created->tails = (double *)((unsigned char *)created + tail_offset);
  *trace = created;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_power_lanczos_primitive_trace_destroy(
    MVMCPowerLanczosPrimitiveTrace *trace) {
  size_t bytes;
  if (trace == NULL) return;
  bytes = trace->allocated_bytes;
  if (bytes >= sizeof(*trace)) memset(trace, 0, bytes);
  free(trace);
}

size_t mvmc_power_lanczos_primitive_trace_allocated_bytes(
    const MVMCPowerLanczosPrimitiveTrace *trace) {
  return trace != NULL ? trace->allocated_bytes : 0;
}

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_register_group(
    MVMCPowerLanczosPrimitiveTrace *trace, size_t group_ordinal,
    size_t leader_world_rank, size_t chain_size) {
  MVMCPowerLanczosPrimitiveTraceGroup candidate;
  size_t group;
  if (trace == NULL || !trace->valid || trace->state != TRACE_BUILDING ||
      group_ordinal >= trace->group_count || chain_size == 0) {
    return fail_trace(trace, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
  }
  if (trace->groups[group_ordinal].valid) {
    return fail_trace(trace, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
  }
  for (group = 0; group < trace->group_count; ++group) {
    if (trace->groups[group].valid &&
        trace->groups[group].leader_world_rank == leader_world_rank) {
      return fail_trace(trace, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
    }
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.stage = trace->stage;
  candidate.group_ordinal = group_ordinal;
  candidate.leader_world_rank = leader_world_rank;
  candidate.chain_size = chain_size;
  candidate.sample_count = trace->samples_per_group;
  trace->groups[group_ordinal] = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_append(
    MVMCPowerLanczosPrimitiveTrace *trace, size_t group_ordinal,
    size_t sample_index, const double complex *values,
    const double *absolute_numeric_bounds, const uint8_t *support_flags,
    size_t primitive_count, double tail_magnitude) {
  size_t primitive;
  if (trace == NULL || !trace->valid || trace->state != TRACE_BUILDING ||
      group_ordinal >= trace->group_count ||
      !trace->groups[group_ordinal].valid || values == NULL ||
      absolute_numeric_bounds == NULL || support_flags == NULL ||
      primitive_count != trace->primitive_count ||
      sample_index != trace->next_sample[group_ordinal] ||
      sample_index >= trace->samples_per_group ||
      !isfinite(tail_magnitude) || tail_magnitude < 0.0) {
    return fail_trace(trace, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
  }
  for (primitive = 0; primitive < primitive_count; ++primitive) {
    const uint8_t flags = support_flags[primitive];
    const double bound = absolute_numeric_bounds[primitive];
    const double complex value = values[primitive];
    if (!isfinite(creal(value)) || !isfinite(cimag(value)) ||
        !isfinite(bound) || bound < 0.0 || flags == 0 ||
        (flags & (uint8_t)~SUPPORT_MASK) != 0 ||
        ((flags & MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO) != 0 &&
         bound <= 0.0) ||
        (flags == MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_EXACT_ZERO &&
         (creal(value) != 0.0 || cimag(value) != 0.0 || bound != 0.0)) ||
        (flags == MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_NUMERIC_ZERO &&
         (creal(value) != 0.0 || cimag(value) != 0.0))) {
      return fail_trace(trace, MVMC_KRYLOV_STATUS_NONFINITE);
    }
  }
  for (primitive = 0; primitive < primitive_count; ++primitive) {
    const size_t offset = value_offset(trace, primitive, group_ordinal,
                                       sample_index);
    trace->values[offset] = values[primitive];
    trace->bounds[offset] = absolute_numeric_bounds[primitive];
    trace->support[offset] = support_flags[primitive];
  }
  trace->tails[group_ordinal * trace->samples_per_group + sample_index] =
      tail_magnitude;
  ++trace->next_sample[group_ordinal];
  ++trace->appended_sample_count;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_freeze(
    MVMCPowerLanczosPrimitiveTrace *trace) {
  size_t group;
  if (trace == NULL || !trace->valid || trace->state != TRACE_BUILDING) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (trace->appended_sample_count != trace->sample_slot_count) {
    return fail_trace(trace, MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE);
  }
  for (group = 0; group < trace->group_count; ++group) {
    if (!trace->groups[group].valid ||
        trace->groups[group].group_ordinal != group ||
        trace->groups[group].stage != trace->stage ||
        trace->groups[group].sample_count != trace->samples_per_group ||
        trace->next_sample[group] != trace->samples_per_group ||
        (group != 0 &&
         trace->groups[group - 1].leader_world_rank >=
             trace->groups[group].leader_world_rank)) {
      return fail_trace(trace,
                        MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE);
    }
  }
  trace->state = TRACE_FROZEN;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_summary(
    const MVMCPowerLanczosPrimitiveTrace *trace,
    MVMCPowerLanczosPrimitiveTraceSummary *summary) {
  MVMCPowerLanczosPrimitiveTraceSummary candidate;
  if (summary == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(summary, 0, sizeof(*summary));
  summary->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (trace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = trace->valid;
  candidate.status = trace->status;
  candidate.version = MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_VERSION;
  candidate.stage = trace->stage;
  candidate.frozen = trace->state == TRACE_FROZEN;
  candidate.primitive_count = trace->primitive_count;
  candidate.scalar_count = trace->scalar_count;
  candidate.group_count = trace->group_count;
  candidate.samples_per_group = trace->samples_per_group;
  candidate.appended_sample_count = trace->appended_sample_count;
  candidate.allocated_bytes = trace->allocated_bytes;
  *summary = candidate;
  return candidate.status;
}

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_group(
    const MVMCPowerLanczosPrimitiveTrace *trace, size_t group_ordinal,
    MVMCPowerLanczosPrimitiveTraceGroup *group) {
  if (group == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(group, 0, sizeof(*group));
  if (trace == NULL || !trace->valid || trace->state != TRACE_FROZEN ||
      group_ordinal >= trace->group_count ||
      !trace->groups[group_ordinal].valid) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *group = trace->groups[group_ordinal];
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_sample(
    const MVMCPowerLanczosPrimitiveTrace *trace, size_t group_ordinal,
    size_t sample_index, double complex *values, size_t value_capacity,
    double *absolute_numeric_bounds, size_t bound_capacity,
    uint8_t *support_flags, size_t support_capacity,
    double *tail_magnitude) {
  size_t primitive;
  if (values != NULL && value_capacity != 0) {
    memset(values, 0, value_capacity * sizeof(*values));
  }
  if (absolute_numeric_bounds != NULL && bound_capacity != 0) {
    memset(absolute_numeric_bounds, 0,
           bound_capacity * sizeof(*absolute_numeric_bounds));
  }
  if (support_flags != NULL && support_capacity != 0) {
    memset(support_flags, 0, support_capacity * sizeof(*support_flags));
  }
  if (tail_magnitude != NULL) *tail_magnitude = 0.0;
  if (trace == NULL || !trace->valid || trace->state != TRACE_FROZEN ||
      group_ordinal >= trace->group_count ||
      sample_index >= trace->samples_per_group || values == NULL ||
      absolute_numeric_bounds == NULL || support_flags == NULL ||
      tail_magnitude == NULL || value_capacity < trace->primitive_count ||
      bound_capacity < trace->primitive_count ||
      support_capacity < trace->primitive_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (primitive = 0; primitive < trace->primitive_count; ++primitive) {
    const size_t offset = value_offset(trace, primitive, group_ordinal,
                                       sample_index);
    values[primitive] = trace->values[offset];
    absolute_numeric_bounds[primitive] = trace->bounds[offset];
    support_flags[primitive] = trace->support[offset];
  }
  *tail_magnitude =
      trace->tails[group_ordinal * trace->samples_per_group + sample_index];
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_primitive_trace_export_scalars(
    const MVMCPowerLanczosPrimitiveTrace *trace, double *scalars,
    size_t scalar_capacity) {
  size_t required;
  size_t primitive;
  size_t group;
  size_t sample;
  if (trace == NULL || !trace->valid || trace->state != TRACE_FROZEN ||
      scalars == NULL ||
      !checked_multiply(trace->scalar_count, trace->sample_slot_count,
                        &required) ||
      scalar_capacity < required) {
    if (scalars != NULL && scalar_capacity != 0) {
      memset(scalars, 0, scalar_capacity * sizeof(*scalars));
    }
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(scalars, 0, required * sizeof(*scalars));
  for (primitive = 0; primitive < trace->primitive_count; ++primitive) {
    for (group = 0; group < trace->group_count; ++group) {
      for (sample = 0; sample < trace->samples_per_group; ++sample) {
        const size_t input = value_offset(trace, primitive, group, sample);
        const size_t slot = group * trace->samples_per_group + sample;
        scalars[(2 * primitive) * trace->sample_slot_count + slot] =
            creal(trace->values[input]);
        scalars[(2 * primitive + 1) * trace->sample_slot_count + slot] =
            cimag(trace->values[input]);
      }
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}
