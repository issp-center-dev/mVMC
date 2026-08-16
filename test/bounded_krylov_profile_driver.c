/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "absolute_pfaffian.h"
#include "bounded_krylov_collective.h"
#include "classic_krylov_model.h"
#include "classic_pfaffian_matrix.h"

#include <complex.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

enum {
  PROFILE_SCHEMA_VERSION = 1,
  PROFILE_MAX_SITE_COUNT = 16,
  PROFILE_DEPTH_COUNT = MVMC_KRYLOV_MAX_ORDER + 1
};

typedef struct {
  MVMCClassicKrylovRawModel raw;
  MVMCClassicKrylovTransfer *transfers;
  MVMCClassicKrylovSiteCoupling *intra;
  MVMCClassicKrylovPairCoupling inter[2];
  MVMCClassicKrylovPairCoupling hund[2];
} ProfileModel;

typedef struct {
  size_t site_count;
  size_t up_electron_count;
  size_t down_electron_count;
  size_t orbital_count;
  size_t word_count;
  int nsize;
  int qp_total;
  int qp_start;
  int qp_end;
  size_t local_qp_count;
  size_t slater_matrix_count;
  size_t matrix_count;
  double scaled_pivot_tolerance;
  double projection_parameter;
  double complex global_weights[4];
  double *slater;
  double *matrix;
  int *ele_idx;
  int *ele_num;
  MVMCAbsolutePfaffianRealValueWorkspace *value_workspace;
  MVMCAbsolutePfaffianScaledValueResult local_components[4];
  MVMCKrylovBoundedCollectiveWorkspace *collective;
  uint64_t callback_calls;
} ProfileScaledAmplitude;

typedef enum {
  COUNTER_ROOT_EVALUATIONS = 0,
  COUNTER_NODE_EXPANSIONS,
  COUNTER_RECURSION_0,
  COUNTER_RECURSION_1,
  COUNTER_RECURSION_2,
  COUNTER_RECURSION_3,
  COUNTER_CACHE_HITS_0,
  COUNTER_CACHE_HITS_1,
  COUNTER_CACHE_HITS_2,
  COUNTER_CACHE_HITS_3,
  COUNTER_CACHE_MISSES_0,
  COUNTER_CACHE_MISSES_1,
  COUNTER_CACHE_MISSES_2,
  COUNTER_CACHE_MISSES_3,
  COUNTER_CACHE_INSERTIONS,
  COUNTER_CACHE_EVICTIONS,
  COUNTER_CACHE_EPOCH_CLEARS,
  COUNTER_RAW_TRANSITIONS,
  COUNTER_MERGED_DUPLICATES,
  COUNTER_CANCELLED_ZERO,
  COUNTER_TERMINAL_CALLS,
  COUNTER_WELL_PIVOTED_COMPONENTS,
  COUNTER_NEAR_PIVOT_COMPONENTS,
  COUNTER_EXACT_ZERO_COMPONENTS,
  COUNTER_NUMERIC_ZERO_COMPONENTS,
  COUNTER_LOCAL_FACTORIZATIONS,
  COUNTER_GLOBAL_FACTORIZATIONS,
  COUNTER_COMPUTED_EXACT_ZERO,
  COUNTER_COMPUTED_NUMERIC_ZERO,
  COUNTER_TRACE_HASH,
  COUNTER_ROW_TRANSITION_PEAK,
  COUNTER_CACHE_ENTRIES_PEAK,
  COUNTER_PLAN_BYTES,
  COUNTER_MODEL_BYTES,
  COUNTER_ENGINE_WORKSPACE_BYTES,
  COUNTER_FRAME_BYTES,
  COUNTER_SCRATCH_BYTES,
  COUNTER_CACHE_REQUESTED_BYTES,
  COUNTER_CACHE_ALLOCATED_BYTES,
  COUNTER_CACHE_SET_COUNT,
  COUNTER_COLLECTIVE_WORKSPACE_BYTES,
  COUNTER_AMPLITUDE_WORKSPACE_BYTES,
  COUNTER_ENGINE_HEAP_ALLOCATIONS,
  COUNTER_RSS_BYTES,
  PROFILE_COUNTER_COUNT
} ProfileCounter;

typedef enum {
  TIME_RESET = 0,
  TIME_EVALUATION,
  TIME_DEPTH_0,
  TIME_DEPTH_1,
  TIME_DEPTH_2,
  TIME_DEPTH_3,
  TIME_AMPLITUDE,
  TIME_CONNECTIVITY,
  TIME_CACHE,
  TIME_TOTAL,
  PROFILE_TIME_COUNT
} ProfileTime;

static int world_rank(void) {
#ifdef _mpi_use
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
#else
  return 0;
#endif
}

static int world_size(void) {
#ifdef _mpi_use
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
#else
  return 1;
#endif
}

static MVMCKrylovBoundedCommunicator world_communicator(void) {
#ifdef _mpi_use
  return MPI_COMM_WORLD;
#else
  return 0;
#endif
}

static int collective_all_ready(int local_ready) {
#ifdef _mpi_use
  int global_ready = 0;
  if (MPI_Allreduce(&local_ready, &global_ready, 1, MPI_INT, MPI_MIN,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    return 0;
  }
  return global_ready;
#else
  return local_ready;
#endif
}

static uint64_t resident_set_bytes(void) {
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) != 0 || usage.ru_maxrss < 0) return 0;
#if defined(__APPLE__)
  return (uint64_t)usage.ru_maxrss;
#else
  if ((uint64_t)usage.ru_maxrss > UINT64_MAX / UINT64_C(1024)) {
    return UINT64_MAX;
  }
  return (uint64_t)usage.ru_maxrss * UINT64_C(1024);
#endif
}

static int parse_int(const char *text, int minimum, int maximum, int *value) {
  char *end = NULL;
  long parsed;
  if (text == NULL || value == NULL || *text == '\0') return 0;
  errno = 0;
  parsed = strtol(text, &end, 10);
  if (errno != 0 || end == text || *end != '\0' || parsed < minimum ||
      parsed > maximum) {
    return 0;
  }
  *value = (int)parsed;
  return 1;
}

static int parse_size_arg(const char *text, size_t maximum, size_t *value) {
  char *end = NULL;
  unsigned long long parsed;
  if (text == NULL || value == NULL || *text == '\0') return 0;
  errno = 0;
  parsed = strtoull(text, &end, 10);
  if (errno != 0 || end == text || *end != '\0' || parsed > maximum ||
      parsed > (unsigned long long)SIZE_MAX) {
    return 0;
  }
  *value = (size_t)parsed;
  return 1;
}

static int popcount_u32(unsigned int value) {
  int count = 0;
  while (value != 0U) {
    value &= value - 1U;
    ++count;
  }
  return count;
}

static size_t fixed_masks(int site_count, unsigned int **masks) {
  const unsigned int end = 1U << site_count;
  const int electron_count = site_count / 2;
  unsigned int *created;
  unsigned int mask;
  size_t count = 0;
  size_t cursor = 0;

  if (masks == NULL) return 0;
  *masks = NULL;
  for (mask = 0; mask < end; ++mask) {
    if (popcount_u32(mask) == electron_count) ++count;
  }
  created = (unsigned int *)malloc(count * sizeof(*created));
  if (created == NULL) return 0;
  for (mask = 0; mask < end; ++mask) {
    if (popcount_u32(mask) == electron_count) created[cursor++] = mask;
  }
  *masks = created;
  return count;
}

static int checked_multiply_size(size_t left, size_t right, size_t *value) {
  if (value == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *value = left * right;
  return 1;
}

static int checked_add_size(size_t left, size_t right, size_t *value) {
  if (value == NULL || right > SIZE_MAX - left) return 0;
  *value = left + right;
  return 1;
}

static int initialize_model(ProfileModel *fixture, int site_count) {
  const int transfer_count = 4 * site_count + 1;
  int cursor = 0;
  int bond;
  int spin;

  if (fixture == NULL) return 0;
  memset(fixture, 0, sizeof(*fixture));
  fixture->transfers = (MVMCClassicKrylovTransfer *)calloc(
      (size_t)transfer_count, sizeof(*fixture->transfers));
  fixture->intra = (MVMCClassicKrylovSiteCoupling *)calloc(
      (size_t)site_count, sizeof(*fixture->intra));
  if (fixture->transfers == NULL || fixture->intra == NULL) return 0;

  for (bond = 0; bond < site_count; ++bond) {
    const int left = bond;
    const int right = (bond + 1) % site_count;
    const double hopping = 0.75 + 0.015625 * (double)(bond % 5);
    for (spin = 0; spin < 2; ++spin) {
      fixture->transfers[cursor++] =
          (MVMCClassicKrylovTransfer){left, spin, right, spin, hopping};
      fixture->transfers[cursor++] =
          (MVMCClassicKrylovTransfer){right, spin, left, spin, hopping};
    }
  }
  fixture->transfers[cursor++] =
      (MVMCClassicKrylovTransfer){1, 0, 1, 0, 0.21875};
  if (cursor != transfer_count) return 0;
  for (bond = 0; bond < site_count; ++bond) {
    const double coefficient =
        (bond & 1) ? -0.0625 * (double)(1 + bond % 3)
                   : 0.09375 * (double)(1 + bond % 3);
    fixture->intra[bond] =
        (MVMCClassicKrylovSiteCoupling){bond, coefficient};
  }
  fixture->inter[0] =
      (MVMCClassicKrylovPairCoupling){0, site_count / 2, 0.078125};
  fixture->inter[1] = (MVMCClassicKrylovPairCoupling){
      1, (1 + site_count / 2) % site_count, -0.046875};
  fixture->hund[0] =
      (MVMCClassicKrylovPairCoupling){0, 1, 0.0546875};
  fixture->hund[1] = (MVMCClassicKrylovPairCoupling){
      site_count / 2, site_count / 2 + 1, -0.0390625};

  fixture->raw.site_count = site_count;
  fixture->raw.up_electron_count = site_count / 2;
  fixture->raw.down_electron_count = site_count / 2;
  fixture->raw.transfer_count = transfer_count;
  fixture->raw.transfers = fixture->transfers;
  fixture->raw.coulomb_intra_count = site_count;
  fixture->raw.coulomb_intra = fixture->intra;
  fixture->raw.coulomb_inter_count = 2;
  fixture->raw.coulomb_inter = fixture->inter;
  fixture->raw.hund_count = 2;
  fixture->raw.hund = fixture->hund;
  return 1;
}

static void destroy_model(ProfileModel *fixture) {
  if (fixture == NULL) return;
  free(fixture->intra);
  free(fixture->transfers);
  memset(fixture, 0, sizeof(*fixture));
}

static int initialize_slater(int site_count, int qp_total, double **slater,
                             double complex **weights) {
  const size_t dimension = (size_t)2 * (size_t)site_count;
  const size_t matrix_count = dimension * dimension;
  double *created_slater;
  double complex *created_weights;
  int qp;
  int up;
  int down;

  if (slater == NULL || weights == NULL) return 0;
  *slater = NULL;
  *weights = NULL;
  created_slater =
      (double *)calloc((size_t)qp_total * matrix_count, sizeof(double));
  created_weights =
      (double complex *)calloc((size_t)qp_total, sizeof(double complex));
  if (created_slater == NULL || created_weights == NULL) {
    free(created_weights);
    free(created_slater);
    return 0;
  }
  for (qp = 0; qp < qp_total; ++qp) {
    static const double weight_values[4] = {1.0, 0.5, -0.25, 0.125};
    created_weights[qp] = weight_values[qp];
    for (up = 0; up < site_count; ++up) {
      for (down = 0; down < site_count; ++down) {
        const size_t base = (size_t)qp * matrix_count;
        const size_t up_orbital = (size_t)up;
        const size_t down_orbital = (size_t)site_count + (size_t)down;
        const double value =
            sin(0.31 * (double)(qp + 1) * (double)(up + 1) +
                0.17 * (double)(down + 1)) +
            0.5 * cos(0.11 * (double)(qp + 2) * (double)(up + 1) *
                      (double)(down + 1)) +
            (up == down ? 0.25 : 0.0);
        created_slater[base + down_orbital * dimension + up_orbital] =
            -value;
        created_slater[base + up_orbital * dimension + down_orbital] =
            value;
      }
    }
  }
  *slater = created_slater;
  *weights = created_weights;
  return 1;
}

static int build_occupation(ProfileScaledAmplitude *workspace,
                            const uint64_t *words, size_t word_count) {
  size_t orbital;
  size_t up_count = 0;
  size_t down_count = 0;

  if (workspace == NULL || words == NULL ||
      word_count != workspace->word_count) {
    return 0;
  }
  if ((workspace->orbital_count % 64) != 0 &&
      (words[word_count - 1] >> (workspace->orbital_count % 64)) != 0) {
    return 0;
  }
  memset(workspace->ele_num, 0,
         workspace->orbital_count * sizeof(*workspace->ele_num));
  for (orbital = 0; orbital < workspace->orbital_count; ++orbital) {
    const int occupied =
        (int)((words[orbital / 64] >> (orbital % 64)) & UINT64_C(1));
    workspace->ele_num[orbital] = occupied;
    if (!occupied) continue;
    if (orbital < workspace->site_count) {
      if (up_count >= workspace->up_electron_count) return 0;
      workspace->ele_idx[up_count++] = (int)orbital;
    } else {
      if (down_count >= workspace->down_electron_count) return 0;
      workspace->ele_idx[workspace->up_electron_count + down_count++] =
          (int)(orbital - workspace->site_count);
    }
  }
  return up_count == workspace->up_electron_count &&
         down_count == workspace->down_electron_count;
}

static double projection_log_value(ProfileScaledAmplitude *workspace) {
  size_t site;
  size_t double_occupancy = 0;

  for (site = 0; site < workspace->site_count; ++site) {
    double_occupancy +=
        (size_t)(workspace->ele_num[site] *
                 workspace->ele_num[workspace->site_count + site]);
  }
  return workspace->projection_parameter * (double)double_occupancy;
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

static MVMCKrylovStatus profile_scaled_amplitude(
    const uint64_t *configuration_words, size_t word_count,
    void *context, MVMCKrylovScaledAmplitudeResult *result) {
  ProfileScaledAmplitude *workspace = (ProfileScaledAmplitude *)context;
  MVMCKrylovScaledAmplitudeResult projected;
  MVMCKrylovBoundedCollectiveResult collective_result;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  MVMCPfaffianStatus pf_status = MVMC_PFAFFIAN_STATUS_OK;
  double log_projection;
  size_t local_qp;
  int all_valid;

  if (workspace == NULL || result == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(result, 0, sizeof(*result));
  ++workspace->callback_calls;
  all_valid = build_occupation(workspace, configuration_words, word_count);
  if (!collective_all_ready(all_valid)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  log_projection = projection_log_value(workspace);
  if (!isfinite(log_projection)) status = MVMC_KRYLOV_STATUS_NONFINITE;
  for (local_qp = 0; status == MVMC_KRYLOV_STATUS_OK &&
                     local_qp < workspace->local_qp_count; ++local_qp) {
    const int global_qp = workspace->qp_start + (int)local_qp;
    pf_status = mvmc_classic_pfaffian_build_real_matrix(
        workspace->slater +
            (size_t)global_qp * workspace->slater_matrix_count,
        (int)workspace->site_count, workspace->nsize, workspace->ele_idx,
        workspace->matrix);
    if (pf_status == MVMC_PFAFFIAN_STATUS_OK) {
      pf_status = mvmc_absolute_pfaffian_real_scaled_value_with_workspace(
          workspace->value_workspace, workspace->matrix, workspace->nsize,
          workspace->nsize, workspace->scaled_pivot_tolerance,
          &workspace->local_components[local_qp]);
    }
    if (pf_status != MVMC_PFAFFIAN_STATUS_OK) {
      status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    } else if (workspace->local_components[local_qp].factor_state ==
               MVMC_PFAFFIAN_VALUE_NONFINITE) {
      status = MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  status = mvmc_bounded_krylov_collective_gather_projected_amplitude(
      workspace->collective, status, workspace->local_components,
      workspace->local_qp_count, workspace->global_weights,
      workspace->qp_total, workspace->qp_start, workspace->qp_end,
      &projected, &collective_result);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (log_projection != 0.0) {
    MVMCScaledComplex projection;
    MVMCScaledComplex scaled;
    if (mvmc_scaled_complex_make_finite(1.0 + 0.0 * I, log_projection,
                                        -INFINITY, &projection) !=
            MVMC_PFAFFIAN_STATUS_OK ||
        mvmc_scaled_complex_multiply(&projection, &projected.value,
                                     &scaled) != MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    projected.value = scaled;
  }
  *result = projected;
  return MVMC_KRYLOV_STATUS_OK;
}

static int create_profile_amplitude(
    int site_count, int qp_total, const double *slater,
    const double complex *weights, double projection_parameter,
    MVMCKrylovBoundedCollectiveWorkspace *collective,
    ProfileScaledAmplitude **workspace) {
  ProfileScaledAmplitude *created;
  size_t dimension = (size_t)2 * (size_t)site_count;
  size_t slater_matrix_count;
  size_t slater_count;
  size_t slater_bytes;
  int qp;

  if (workspace == NULL || slater == NULL || weights == NULL ||
      collective == NULL) {
    return 0;
  }
  *workspace = NULL;
  created = (ProfileScaledAmplitude *)calloc(1, sizeof(*created));
  if (created == NULL) return 0;
  created->site_count = (size_t)site_count;
  created->up_electron_count = (size_t)(site_count / 2);
  created->down_electron_count = (size_t)(site_count / 2);
  created->orbital_count = dimension;
  created->word_count = (dimension + 63) / 64;
  created->nsize = site_count;
  created->qp_total = qp_total;
  created->qp_start = world_rank() * qp_total / world_size();
  created->qp_end = (world_rank() + 1) * qp_total / world_size();
  created->local_qp_count =
      (size_t)(created->qp_end - created->qp_start);
  created->matrix_count = (size_t)created->nsize * (size_t)created->nsize;
  created->scaled_pivot_tolerance = 0.0;
  created->projection_parameter = projection_parameter;
  created->collective = collective;
  for (qp = 0; qp < qp_total; ++qp) created->global_weights[qp] = weights[qp];
  if (!checked_multiply_size(dimension, dimension, &slater_matrix_count) ||
      !checked_multiply_size((size_t)qp_total, slater_matrix_count,
                             &slater_count) ||
      !checked_multiply_size(slater_count, sizeof(*created->slater),
                             &slater_bytes)) {
    free(created);
    return 0;
  }
  created->slater_matrix_count = slater_matrix_count;
  created->slater = (double *)malloc(slater_bytes);
  created->matrix =
      (double *)malloc(created->matrix_count * sizeof(*created->matrix));
  created->ele_idx =
      (int *)calloc((size_t)created->nsize, sizeof(*created->ele_idx));
  created->ele_num =
      (int *)calloc(created->orbital_count, sizeof(*created->ele_num));
  if (created->slater == NULL || created->matrix == NULL ||
      created->ele_idx == NULL || created->ele_num == NULL ||
      mvmc_absolute_pfaffian_real_value_workspace_create(
          created->nsize, &created->value_workspace) !=
          MVMC_PFAFFIAN_STATUS_OK) {
    mvmc_absolute_pfaffian_real_value_workspace_destroy(
        created->value_workspace);
    free(created->ele_num);
    free(created->ele_idx);
    free(created->matrix);
    free(created->slater);
    free(created);
    return 0;
  }
  memcpy(created->slater, slater, slater_bytes);
  *workspace = created;
  return 1;
}

static void destroy_profile_amplitude(ProfileScaledAmplitude *workspace) {
  if (workspace == NULL) return;
  mvmc_absolute_pfaffian_real_value_workspace_destroy(
      workspace->value_workspace);
  free(workspace->ele_num);
  free(workspace->ele_idx);
  free(workspace->matrix);
  free(workspace->slater);
  free(workspace);
}

static size_t profile_amplitude_workspace_bytes(
    const ProfileScaledAmplitude *workspace) {
  size_t bytes;
  size_t item;
  if (workspace == NULL) return 0;
  bytes = sizeof(*workspace);
#define ADD_BYTES(count, type)                                             \
  do {                                                                     \
    if (!checked_multiply_size((count), sizeof(type), &item) ||            \
        !checked_add_size(bytes, item, &bytes)) return 0;                  \
  } while (0)
  ADD_BYTES((size_t)workspace->nsize, int);
  ADD_BYTES(workspace->orbital_count, int);
  ADD_BYTES((size_t)workspace->qp_total * workspace->slater_matrix_count,
            double);
  ADD_BYTES(workspace->matrix_count, double);
#undef ADD_BYTES
  if (!checked_add_size(
          bytes,
          mvmc_absolute_pfaffian_real_value_workspace_bytes(
              workspace->value_workspace),
          &bytes)) {
    return 0;
  }
  return bytes;
}

static int profile_limits(size_t cache_bytes,
                          const MVMCKrylovFockModel *model,
                          MVMCKrylovBoundedLimits *limits) {
  if (limits == NULL || model == NULL) return 0;
  memset(limits, 0, sizeof(*limits));
  limits->amplitude_policy_hash = UINT64_C(0x5034435f5245414c);
  limits->cache_bytes = cache_bytes;
  limits->max_row_transitions = 4096;
  limits->max_workspace_bytes = (size_t)4 * 1024 * 1024 * 1024;
  if (limits->max_workspace_bytes < cache_bytes) {
    limits->max_workspace_bytes = SIZE_MAX;
  }
  limits->max_node_expansions = UINT64_C(100000000000);
  limits->max_terminal_amplitude_calls = UINT64_C(100000000000);
  limits->max_total_row_transitions = UINT64_C(100000000000);
  limits->max_order = MVMC_KRYLOV_MAX_ORDER;
  (void)model;
  return 1;
}

static void add_statistics(uint64_t *counters, double *times,
                           const MVMCKrylovBoundedStatistics *statistics) {
  int depth;
  counters[COUNTER_ROOT_EVALUATIONS] += statistics->root_evaluations;
  counters[COUNTER_NODE_EXPANSIONS] += statistics->node_expansions;
  for (depth = 0; depth < PROFILE_DEPTH_COUNT; ++depth) {
    counters[COUNTER_RECURSION_0 + depth] +=
        statistics->recursion_calls[depth];
    counters[COUNTER_CACHE_HITS_0 + depth] += statistics->cache_hits[depth];
    counters[COUNTER_CACHE_MISSES_0 + depth] +=
        statistics->cache_misses[depth];
    times[TIME_DEPTH_0 + depth] += statistics->depth_wall_seconds[depth];
  }
  counters[COUNTER_CACHE_INSERTIONS] += statistics->cache_insertions;
  counters[COUNTER_CACHE_EVICTIONS] += statistics->cache_evictions;
  counters[COUNTER_CACHE_EPOCH_CLEARS] +=
      statistics->cache_epoch_full_clears;
  counters[COUNTER_RAW_TRANSITIONS] += statistics->raw_row_transitions;
  counters[COUNTER_MERGED_DUPLICATES] +=
      statistics->merged_duplicate_transitions;
  counters[COUNTER_CANCELLED_ZERO] +=
      statistics->cancelled_zero_transitions;
  counters[COUNTER_TERMINAL_CALLS] += statistics->terminal_amplitude_calls;
  counters[COUNTER_WELL_PIVOTED_COMPONENTS] +=
      statistics->well_pivoted_component_count;
  counters[COUNTER_NEAR_PIVOT_COMPONENTS] +=
      statistics->near_pivot_component_count;
  counters[COUNTER_EXACT_ZERO_COMPONENTS] +=
      statistics->exact_zero_component_count;
  counters[COUNTER_NUMERIC_ZERO_COMPONENTS] +=
      statistics->numeric_zero_component_count;
  counters[COUNTER_LOCAL_FACTORIZATIONS] +=
      statistics->local_factorization_count;
  counters[COUNTER_GLOBAL_FACTORIZATIONS] +=
      statistics->global_factorization_count;
  counters[COUNTER_COMPUTED_EXACT_ZERO] +=
      statistics->computed_exact_zero_values;
  counters[COUNTER_COMPUTED_NUMERIC_ZERO] +=
      statistics->computed_numeric_zero_values;
  counters[COUNTER_TRACE_HASH] = statistics->trace_hash;
  if ((uint64_t)statistics->row_transition_peak >
      counters[COUNTER_ROW_TRANSITION_PEAK]) {
    counters[COUNTER_ROW_TRANSITION_PEAK] =
        (uint64_t)statistics->row_transition_peak;
  }
  if ((uint64_t)statistics->cache_entries_peak >
      counters[COUNTER_CACHE_ENTRIES_PEAK]) {
    counters[COUNTER_CACHE_ENTRIES_PEAK] =
        (uint64_t)statistics->cache_entries_peak;
  }
  counters[COUNTER_PLAN_BYTES] = (uint64_t)statistics->plan_bytes;
  counters[COUNTER_MODEL_BYTES] = (uint64_t)statistics->model_bytes;
  counters[COUNTER_ENGINE_WORKSPACE_BYTES] =
      (uint64_t)statistics->workspace_bytes;
  counters[COUNTER_FRAME_BYTES] = (uint64_t)statistics->frame_bytes;
  counters[COUNTER_SCRATCH_BYTES] = (uint64_t)statistics->scratch_bytes;
  counters[COUNTER_CACHE_REQUESTED_BYTES] =
      (uint64_t)statistics->cache_requested_bytes;
  counters[COUNTER_CACHE_ALLOCATED_BYTES] =
      (uint64_t)statistics->cache_allocated_bytes;
  counters[COUNTER_CACHE_SET_COUNT] = (uint64_t)statistics->cache_set_count;
  counters[COUNTER_ENGINE_HEAP_ALLOCATIONS] +=
      statistics->engine_heap_allocations;
  times[TIME_RESET] += statistics->reset_seconds;
  times[TIME_EVALUATION] += statistics->evaluation_wall_seconds;
  times[TIME_AMPLITUDE] += statistics->amplitude_wall_seconds;
  times[TIME_CONNECTIVITY] += statistics->connectivity_wall_seconds;
  times[TIME_CACHE] += statistics->cache_wall_seconds;
  times[TIME_TOTAL] += statistics->total_seconds;
}

static void reduce_profile(const uint64_t *local_counters,
                           const double *local_times,
                           uint64_t *sum_counters, uint64_t *max_counters,
                           double *sum_times, double *max_times) {
#ifdef _mpi_use
  MPI_Reduce(local_counters, sum_counters, PROFILE_COUNTER_COUNT,
             MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(local_counters, max_counters, PROFILE_COUNTER_COUNT,
             MPI_UINT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(local_times, sum_times, PROFILE_TIME_COUNT, MPI_DOUBLE, MPI_SUM,
             0, MPI_COMM_WORLD);
  MPI_Reduce(local_times, max_times, PROFILE_TIME_COUNT, MPI_DOUBLE, MPI_MAX,
             0, MPI_COMM_WORLD);
#else
  memcpy(sum_counters, local_counters,
         PROFILE_COUNTER_COUNT * sizeof(*sum_counters));
  memcpy(max_counters, local_counters,
         PROFILE_COUNTER_COUNT * sizeof(*max_counters));
  memcpy(sum_times, local_times, PROFILE_TIME_COUNT * sizeof(*sum_times));
  memcpy(max_times, local_times, PROFILE_TIME_COUNT * sizeof(*max_times));
#endif
}

static void print_resource(const char *scope, const uint64_t *counter,
                           const double *time) {
  printf("RESOURCE scope=%s roots=%" PRIu64, scope,
         counter[COUNTER_ROOT_EVALUATIONS]);
  printf(" node_expansions=%" PRIu64,
         counter[COUNTER_NODE_EXPANSIONS]);
  printf(" recursion=%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64,
         counter[COUNTER_RECURSION_0], counter[COUNTER_RECURSION_1],
         counter[COUNTER_RECURSION_2], counter[COUNTER_RECURSION_3]);
  printf(" cache_hits=%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64,
         counter[COUNTER_CACHE_HITS_0], counter[COUNTER_CACHE_HITS_1],
         counter[COUNTER_CACHE_HITS_2], counter[COUNTER_CACHE_HITS_3]);
  printf(" cache_misses=%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64,
         counter[COUNTER_CACHE_MISSES_0], counter[COUNTER_CACHE_MISSES_1],
         counter[COUNTER_CACHE_MISSES_2], counter[COUNTER_CACHE_MISSES_3]);
  printf(" cache_insertions=%" PRIu64 " cache_evictions=%" PRIu64,
         counter[COUNTER_CACHE_INSERTIONS],
         counter[COUNTER_CACHE_EVICTIONS]);
  printf(" cache_epoch_clears=%" PRIu64,
         counter[COUNTER_CACHE_EPOCH_CLEARS]);
  printf(" raw_transitions=%" PRIu64 " merged_duplicates=%" PRIu64,
         counter[COUNTER_RAW_TRANSITIONS],
         counter[COUNTER_MERGED_DUPLICATES]);
  printf(" cancelled_zero=%" PRIu64 " terminal_calls=%" PRIu64,
         counter[COUNTER_CANCELLED_ZERO], counter[COUNTER_TERMINAL_CALLS]);
  printf(" well_pivoted=%" PRIu64 " near_pivot=%" PRIu64,
         counter[COUNTER_WELL_PIVOTED_COMPONENTS],
         counter[COUNTER_NEAR_PIVOT_COMPONENTS]);
  printf(" exact_zero_components=%" PRIu64,
         counter[COUNTER_EXACT_ZERO_COMPONENTS]);
  printf(" numeric_zero_components=%" PRIu64,
         counter[COUNTER_NUMERIC_ZERO_COMPONENTS]);
  printf(" local_factorizations=%" PRIu64,
         counter[COUNTER_LOCAL_FACTORIZATIONS]);
  printf(" global_factorizations=%" PRIu64,
         counter[COUNTER_GLOBAL_FACTORIZATIONS]);
  printf(" computed_exact_zero=%" PRIu64,
         counter[COUNTER_COMPUTED_EXACT_ZERO]);
  printf(" computed_numeric_zero=%" PRIu64,
         counter[COUNTER_COMPUTED_NUMERIC_ZERO]);
  printf(" trace_hash=%" PRIu64, counter[COUNTER_TRACE_HASH]);
  printf(" row_peak=%" PRIu64 " cache_entries_peak=%" PRIu64,
         counter[COUNTER_ROW_TRANSITION_PEAK],
         counter[COUNTER_CACHE_ENTRIES_PEAK]);
  printf(" plan_bytes=%" PRIu64 " model_bytes=%" PRIu64,
         counter[COUNTER_PLAN_BYTES], counter[COUNTER_MODEL_BYTES]);
  printf(" engine_workspace_bytes=%" PRIu64,
         counter[COUNTER_ENGINE_WORKSPACE_BYTES]);
  printf(" frame_bytes=%" PRIu64 " scratch_bytes=%" PRIu64,
         counter[COUNTER_FRAME_BYTES], counter[COUNTER_SCRATCH_BYTES]);
  printf(" cache_requested_bytes=%" PRIu64,
         counter[COUNTER_CACHE_REQUESTED_BYTES]);
  printf(" cache_allocated_bytes=%" PRIu64,
         counter[COUNTER_CACHE_ALLOCATED_BYTES]);
  printf(" cache_set_count=%" PRIu64,
         counter[COUNTER_CACHE_SET_COUNT]);
  printf(" collective_workspace_bytes=%" PRIu64,
         counter[COUNTER_COLLECTIVE_WORKSPACE_BYTES]);
  printf(" amplitude_workspace_bytes=%" PRIu64,
         counter[COUNTER_AMPLITUDE_WORKSPACE_BYTES]);
  printf(" engine_heap_allocations=%" PRIu64,
         counter[COUNTER_ENGINE_HEAP_ALLOCATIONS]);
  printf(" rss_bytes=%" PRIu64, counter[COUNTER_RSS_BYTES]);
  printf(" reset_seconds=%.17g evaluation_seconds=%.17g",
         time[TIME_RESET], time[TIME_EVALUATION]);
  printf(" depth_seconds=%.17g,%.17g,%.17g,%.17g",
         time[TIME_DEPTH_0], time[TIME_DEPTH_1], time[TIME_DEPTH_2],
         time[TIME_DEPTH_3]);
  printf(" amplitude_seconds=%.17g connectivity_seconds=%.17g",
         time[TIME_AMPLITUDE], time[TIME_CONNECTIVITY]);
  printf(" cache_seconds=%.17g total_seconds=%.17g\n",
         time[TIME_CACHE], time[TIME_TOTAL]);
}

#if defined(MVMC_BOUNDED_KRYLOV_PROFILE_DRIVER_NO_MAIN) &&                   \
    defined(__GNUC__)
#define MVMC_PROFILE_DRIVER_MAYBE_UNUSED __attribute__((unused))
#else
#define MVMC_PROFILE_DRIVER_MAYBE_UNUSED
#endif

static int MVMC_PROFILE_DRIVER_MAYBE_UNUSED run_profile(
    int site_count, int qp_total, int requested_samples, size_t cache_bytes,
    int audit, int persistent_session) {
  ProfileModel fixture;
  MVMCClassicKrylovModelWorkspace *model_workspace = NULL;
  ProfileScaledAmplitude *amplitude_workspace = NULL;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *krylov_workspace = NULL;
  MVMCKrylovBoundedCollectiveWorkspace *collective_workspace = NULL;
  const MVMCKrylovFockModel *model = NULL;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  MVMCKrylovBoundedCollectiveResult collective_result;
  unsigned int *masks = NULL;
  size_t mask_count = 0;
  size_t sector_size = 0;
  size_t sample_count = 0;
  double *slater = NULL;
  double complex *weights = NULL;
  const double projection_parameter = -0.27;
  uint64_t local_counters[PROFILE_COUNTER_COUNT] = {0};
  uint64_t sum_counters[PROFILE_COUNTER_COUNT] = {0};
  uint64_t max_counters[PROFILE_COUNTER_COUNT] = {0};
  double local_times[PROFILE_TIME_COUNT] = {0.0};
  double sum_times[PROFILE_TIME_COUNT] = {0.0};
  double max_times[PROFILE_TIME_COUNT] = {0.0};
  const int rank = world_rank();
  const int size = world_size();
  uint64_t plan_hash = 0;
  size_t sample;
  int local_ready = 1;
  int return_code = 1;

  memset(&fixture, 0, sizeof(fixture));
  mask_count = fixed_masks(site_count, &masks);
  if (mask_count == 0 || !checked_multiply_size(mask_count, mask_count,
                                                &sector_size)) {
    local_ready = 0;
  }
  if (local_ready) {
    sample_count = requested_samples == 0 ? sector_size
                                          : (size_t)requested_samples;
    if (sample_count == 0 || sample_count > sector_size ||
        !initialize_model(&fixture, site_count) ||
        !initialize_slater(site_count, qp_total, &slater, &weights)) {
      local_ready = 0;
    }
  }
  if (!collective_all_ready(local_ready)) {
    if (rank == 0) {
      fprintf(stderr, "bounded profile fixture allocation/preflight failed\n");
    }
    goto cleanup;
  }
  status = mvmc_classic_krylov_model_workspace_create_from_root(
      rank == 0 ? &fixture.raw : NULL, world_communicator(),
      &model_workspace);
  if (!collective_all_ready(status == MVMC_KRYLOV_STATUS_OK)) {
    if (rank == 0) {
      fprintf(stderr, "model create failed collectively: %s\n",
              mvmc_krylov_status_string(status));
    }
    goto cleanup;
  }
  model = mvmc_classic_krylov_model(model_workspace);
  local_ready = model != NULL && profile_limits(cache_bytes, model, &limits);
  if (local_ready) {
    status = mvmc_bounded_krylov_plan_create(model, &limits, &plan);
    local_ready = status == MVMC_KRYLOV_STATUS_OK && plan != NULL;
  }
  if (local_ready) {
    plan_hash = mvmc_bounded_krylov_plan_hash(plan);
    status = mvmc_bounded_krylov_workspace_create(plan, &krylov_workspace);
    local_ready = status == MVMC_KRYLOV_STATUS_OK &&
                  krylov_workspace != NULL;
  }
  if (local_ready) {
    status = mvmc_bounded_krylov_collective_workspace_create(
        world_communicator(), (size_t)qp_total, PROFILE_DEPTH_COUNT,
        &collective_workspace);
    local_ready = status == MVMC_KRYLOV_STATUS_OK &&
                  collective_workspace != NULL;
  }
  if (local_ready) {
    local_ready = create_profile_amplitude(
        site_count, qp_total, slater, weights, projection_parameter,
        collective_workspace, &amplitude_workspace);
    status = local_ready ? MVMC_KRYLOV_STATUS_OK
                         : MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
  status = mvmc_bounded_krylov_collective_synchronize(
      collective_workspace, status, plan_hash, local_ready,
      &collective_result);
  if (!collective_all_ready(status == MVMC_KRYLOV_STATUS_OK &&
                            collective_result.valid && local_ready)) {
    if (rank == 0) {
      fprintf(stderr, "bounded engine preflight failed collectively: %s\n",
              mvmc_krylov_status_string(status));
    }
    goto cleanup;
  }
  local_counters[COUNTER_COLLECTIVE_WORKSPACE_BYTES] =
      (uint64_t)mvmc_bounded_krylov_collective_workspace_bytes(
          collective_workspace);
  local_counters[COUNTER_AMPLITUDE_WORKSPACE_BYTES] =
      (uint64_t)profile_amplitude_workspace_bytes(amplitude_workspace);
  if (rank == 0) {
    printf("%s schema=%d fixture=%s",
           persistent_session ? "SESSION_PROFILE" : "PROFILE",
           persistent_session ? 2 : PROFILE_SCHEMA_VERSION,
           persistent_session ? "p4s9_target_session_spread_roots"
                              : "p4c_bounded_electronic_real");
    printf(" site_count=%d qp_total=%d sample_count=%zu sector_size=%zu",
           site_count, qp_total, sample_count, sector_size);
    printf(" rank_count=%d cache_requested_bytes=%zu audit=%d\n",
           size, cache_bytes, audit);
  }
  if (persistent_session) {
    uint64_t generation_hash = plan_hash ^ limits.amplitude_policy_hash ^
                               UINT64_C(0x5034533953505246);
    if (generation_hash == 0) generation_hash = UINT64_C(1);
    status = mvmc_bounded_krylov_session_begin(
        krylov_workspace, profile_scaled_amplitude, amplitude_workspace,
        generation_hash);
    if (!collective_all_ready(status == MVMC_KRYLOV_STATUS_OK)) {
      if (rank == 0) {
        fprintf(stderr, "bounded session profile begin failed: %s\n",
                mvmc_krylov_status_string(status));
      }
      goto cleanup;
    }
  }
  for (sample = 0; sample < sample_count; ++sample) {
    const size_t sector_index =
        sample_count == 1
            ? (sector_size - 1) / 2
            : sample * (sector_size - 1) / (sample_count - 1);
    const size_t up_index = sector_index / mask_count;
    const size_t down_index = sector_index % mask_count;
    const uint64_t configuration =
        (uint64_t)masks[up_index] |
        ((uint64_t)masks[down_index] << (unsigned int)site_count);
    MVMCKrylovBoundedResult result;
    int order;
    status = persistent_session
                 ? mvmc_bounded_krylov_session_evaluate(
                       krylov_workspace, &configuration, 1, &result)
                 : mvmc_bounded_krylov_evaluate(
                       krylov_workspace, &configuration, 1,
                       profile_scaled_amplitude, amplitude_workspace,
                       &result);
    local_ready = status == MVMC_KRYLOV_STATUS_OK && result.valid &&
                  result.evaluated_order == MVMC_KRYLOV_MAX_ORDER;
    if (!collective_all_ready(local_ready)) {
      if (rank == 0) {
        fprintf(stderr,
                "bounded profile state=%" PRIu64
                " failed collectively: %s\n",
                configuration, mvmc_krylov_status_string(status));
      }
      goto cleanup;
    }
    add_statistics(local_counters, local_times, &result.statistics);
    if (rank == 0) {
      printf("ROW %zu %zu %" PRIu64, sample, sector_index, configuration);
      for (order = 0; order < PROFILE_DEPTH_COUNT; ++order) {
        double complex raw = NAN + I * NAN;
        if (!scaled_to_raw(&result.value[order], &raw)) {
          fprintf(stderr, "bounded profile value export failed at order %d\n",
                  order);
          goto cleanup;
        }
        printf(" %.17g %.17g", creal(raw), cimag(raw));
      }
      putchar('\n');
      if (audit) {
        printf("STAT %zu", sample);
        printf(" reset_seconds=%.17g total_seconds=%.17g",
               result.statistics.reset_seconds,
               result.statistics.total_seconds);
        printf(" depth_seconds=%.17g,%.17g,%.17g,%.17g",
               result.statistics.depth_wall_seconds[0],
               result.statistics.depth_wall_seconds[1],
               result.statistics.depth_wall_seconds[2],
               result.statistics.depth_wall_seconds[3]);
        printf(" terminal_calls=%" PRIu64,
               result.statistics.terminal_amplitude_calls);
        printf(" raw_transitions=%" PRIu64,
               result.statistics.raw_row_transitions);
        printf(" cache_hits=%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64,
               result.statistics.cache_hits[0],
               result.statistics.cache_hits[1],
               result.statistics.cache_hits[2],
               result.statistics.cache_hits[3]);
        printf(" cache_misses=%" PRIu64 ",%" PRIu64 ",%" PRIu64
               ",%" PRIu64,
               result.statistics.cache_misses[0],
               result.statistics.cache_misses[1],
               result.statistics.cache_misses[2],
               result.statistics.cache_misses[3]);
        printf(" row_peak=%zu cache_entries_peak=%zu",
               result.statistics.row_transition_peak,
               result.statistics.cache_entries_peak);
        printf(" session_active=%d session_generation=%" PRIu64,
               result.statistics.persistent_session_active,
               result.statistics.amplitude_generation_hash);
        printf(" session_root_evaluation=%" PRIu64,
               result.statistics.session_root_evaluation);
        printf(" cache_reset_performed=%" PRIu64,
               result.statistics.cache_reset_performed);
        printf(" cache_resident_start=%zu cache_resident_end=%zu",
               result.statistics.cache_entries_resident_start,
               result.statistics.cache_entries_resident_end);
        printf(" workspace_bytes=%zu cache_allocated_bytes=%zu\n",
               result.statistics.workspace_bytes,
               result.statistics.cache_allocated_bytes);
      }
    }
  }
  if (persistent_session) {
    status = mvmc_bounded_krylov_session_end(krylov_workspace);
    if (!collective_all_ready(status == MVMC_KRYLOV_STATUS_OK)) {
      if (rank == 0) {
        fprintf(stderr, "bounded session profile end failed: %s\n",
                mvmc_krylov_status_string(status));
      }
      goto cleanup;
    }
  }
  local_counters[COUNTER_RSS_BYTES] = resident_set_bytes();
  reduce_profile(local_counters, local_times, sum_counters, max_counters,
                 sum_times, max_times);
  if (rank == 0 && audit) {
    print_resource("rank_sum", sum_counters, sum_times);
    print_resource("rank_max", max_counters, max_times);
  }
  return_code = 0;

cleanup:
  destroy_profile_amplitude(amplitude_workspace);
  mvmc_bounded_krylov_collective_workspace_destroy(collective_workspace);
  mvmc_bounded_krylov_workspace_destroy(krylov_workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
  mvmc_classic_krylov_model_workspace_destroy(model_workspace);
  free(weights);
  free(slater);
  destroy_model(&fixture);
  free(masks);
  return return_code;
}

#undef MVMC_PROFILE_DRIVER_MAYBE_UNUSED

#ifndef MVMC_BOUNDED_KRYLOV_PROFILE_DRIVER_NO_MAIN
int main(int argc, char **argv) {
  int site_count;
  int qp_total;
  int sample_count;
  int audit;
  int persistent_session = 0;
  size_t cache_bytes;
  int result;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#endif
  if ((argc != 6 && argc != 7) ||
      !parse_int(argv[1], 4, PROFILE_MAX_SITE_COUNT, &site_count) ||
      (site_count & 1) != 0 || !parse_int(argv[2], 1, 4, &qp_total) ||
      (qp_total != 1 && qp_total != 4) ||
      !parse_int(argv[3], 0, 100000000, &sample_count) ||
      !parse_size_arg(argv[4], (size_t)4 * 1024 * 1024 * 1024,
                      &cache_bytes) ||
      !parse_int(argv[5], 0, 1, &audit) ||
      (argc == 7 && strcmp(argv[6], "session-profile") != 0)) {
    if (world_rank() == 0) {
      fprintf(stderr,
              "usage: %s EVEN_SITE_COUNT(4..16) QP_TOTAL(1|4) "
              "SAMPLES(0=all) CACHE_BYTES AUDIT(0|1) [session-profile]\n",
              argv[0]);
    }
    result = EXIT_FAILURE;
  } else {
    persistent_session = argc == 7;
    result = run_profile(site_count, qp_total, sample_count, cache_bytes,
                         audit, persistent_session) == 0
                 ? EXIT_SUCCESS
                 : EXIT_FAILURE;
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return result;
}
#endif
