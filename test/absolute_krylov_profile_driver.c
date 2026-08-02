/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "classic_krylov_amplitude.h"
#include "classic_krylov_model.h"
#include "krylov_fock_reference.h"

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

typedef enum {
  COUNTER_ROOT_EVALUATIONS = 0,
  COUNTER_RECURSION_0,
  COUNTER_RECURSION_1,
  COUNTER_RECURSION_2,
  COUNTER_RECURSION_3,
  COUNTER_UNIQUE_0,
  COUNTER_UNIQUE_1,
  COUNTER_UNIQUE_2,
  COUNTER_UNIQUE_3,
  COUNTER_MEMO_HITS_0,
  COUNTER_MEMO_HITS_1,
  COUNTER_MEMO_HITS_2,
  COUNTER_MEMO_HITS_3,
  COUNTER_MEMO_MISSES_0,
  COUNTER_MEMO_MISSES_1,
  COUNTER_MEMO_MISSES_2,
  COUNTER_MEMO_MISSES_3,
  COUNTER_RAW_TRANSITIONS,
  COUNTER_MERGED_DUPLICATES,
  COUNTER_CANCELLED_ZERO,
  COUNTER_TERMINAL_REQUESTS,
  COUNTER_TERMINAL_CACHE_HITS,
  COUNTER_REGULAR_COMPONENTS,
  COUNTER_NEAR_PIVOT_COMPONENTS,
  COUNTER_SINGULAR_COMPONENTS,
  COUNTER_TOTAL_ZERO,
  COUNTER_LOCAL_FACTORIZATIONS,
  COUNTER_GLOBAL_FACTORIZATIONS,
  COUNTER_FRONTIER_PEAK,
  COUNTER_MEMO_PEAK,
  COUNTER_KRYLOV_WORKSPACE_BYTES,
  COUNTER_MODEL_WORKSPACE_BYTES,
  COUNTER_AMPLITUDE_WORKSPACE_BYTES,
  COUNTER_RSS_BYTES,
  PROFILE_COUNTER_COUNT
} ProfileCounter;

typedef enum {
  TIME_DEPTH_0 = 0,
  TIME_DEPTH_1,
  TIME_DEPTH_2,
  TIME_DEPTH_3,
  TIME_AMPLITUDE,
  TIME_CONNECTIVITY,
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

static MVMCClassicPfaffianCommunicator world_communicator(void) {
#ifdef _mpi_use
  return MPI_COMM_WORLD;
#else
  return 0;
#endif
}

static double monotonic_seconds(void) {
  struct timespec now;
  if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return -1.0;
  return (double)now.tv_sec + 1.0e-9 * (double)now.tv_nsec;
}

static double elapsed_seconds(double start) {
  const double end = monotonic_seconds();
  return start >= 0.0 && end >= start ? end - start : 0.0;
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

static MVMCClassicKrylovAmplitudeLayout amplitude_layout(
    int site_count, int qp_total, const int *gutzwiller_indices,
    const double complex *projection_parameter) {
  MVMCClassicKrylovAmplitudeLayout layout;
  const int rank = world_rank();
  const int size = world_size();
  memset(&layout, 0, sizeof(layout));
  layout.site_count = (size_t)site_count;
  layout.up_electron_count = (size_t)(site_count / 2);
  layout.down_electron_count = (size_t)(site_count / 2);
  layout.qp_total = qp_total;
  layout.qp_start = rank * qp_total / size;
  layout.qp_end = (rank + 1) * qp_total / size;
  layout.scaled_pivot_tolerance = 0.0;
  layout.nproj = 1;
  layout.ngutzwiller_idx = 1;
  layout.gutzwiller_idx = gutzwiller_indices;
  layout.projection_parameters = projection_parameter;
  layout.communicator = world_communicator();
  return layout;
}

static int checked_multiply_size(size_t left, size_t right, size_t *value) {
  if (value == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *value = left * right;
  return 1;
}

static int profile_limits(size_t sector_size, size_t term_count,
                          MVMCKrylovLimits *limits) {
  size_t max_states;
  size_t max_transitions;
  if (limits == NULL || sector_size > (SIZE_MAX - 64) / 4) return 0;
  max_states = 4 * sector_size + 64;
  if (!checked_multiply_size(max_states, term_count, &max_transitions)) {
    return 0;
  }
  limits->max_states = max_states;
  limits->max_transitions = max_transitions;
  limits->max_amplitude_evaluations = sector_size;
  limits->max_bytes = (size_t)1024 * 1024 * 1024;
  limits->max_order = MVMC_KRYLOV_MAX_ORDER;
  return 1;
}

static void add_statistics(uint64_t *counters, double *times,
                           const MVMCKrylovStatistics *statistics) {
  int depth;
  counters[COUNTER_ROOT_EVALUATIONS] += statistics->root_evaluations;
  for (depth = 0; depth < PROFILE_DEPTH_COUNT; ++depth) {
    counters[COUNTER_RECURSION_0 + depth] +=
        statistics->recursion_calls[depth];
    counters[COUNTER_UNIQUE_0 + depth] += statistics->unique_states[depth];
    counters[COUNTER_MEMO_HITS_0 + depth] += statistics->memo_hits[depth];
    counters[COUNTER_MEMO_MISSES_0 + depth] += statistics->memo_misses[depth];
    times[TIME_DEPTH_0 + depth] += statistics->depth_wall_seconds[depth];
  }
  counters[COUNTER_RAW_TRANSITIONS] += statistics->raw_transitions;
  counters[COUNTER_MERGED_DUPLICATES] +=
      statistics->merged_duplicate_transitions;
  counters[COUNTER_CANCELLED_ZERO] +=
      statistics->cancelled_zero_transitions;
  counters[COUNTER_TERMINAL_REQUESTS] +=
      statistics->terminal_amplitude_requests;
  counters[COUNTER_TERMINAL_CACHE_HITS] +=
      statistics->terminal_cache_hits;
  counters[COUNTER_REGULAR_COMPONENTS] +=
      statistics->regular_component_count;
  counters[COUNTER_NEAR_PIVOT_COMPONENTS] +=
      statistics->near_pivot_component_count;
  counters[COUNTER_SINGULAR_COMPONENTS] +=
      statistics->singular_component_count;
  counters[COUNTER_TOTAL_ZERO] += statistics->total_zero_count;
  counters[COUNTER_LOCAL_FACTORIZATIONS] +=
      statistics->local_factorization_count;
  counters[COUNTER_GLOBAL_FACTORIZATIONS] +=
      statistics->global_factorization_count;
  if ((uint64_t)statistics->frontier_peak >
      counters[COUNTER_FRONTIER_PEAK]) {
    counters[COUNTER_FRONTIER_PEAK] = (uint64_t)statistics->frontier_peak;
  }
  if ((uint64_t)statistics->memo_entries_peak >
      counters[COUNTER_MEMO_PEAK]) {
    counters[COUNTER_MEMO_PEAK] =
        (uint64_t)statistics->memo_entries_peak;
  }
  times[TIME_AMPLITUDE] += statistics->amplitude_wall_seconds;
  times[TIME_CONNECTIVITY] += statistics->connectivity_wall_seconds;
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
  if (world_rank() == 0) {
    sum_counters[COUNTER_REGULAR_COMPONENTS] =
        max_counters[COUNTER_REGULAR_COMPONENTS];
    sum_counters[COUNTER_NEAR_PIVOT_COMPONENTS] =
        max_counters[COUNTER_NEAR_PIVOT_COMPONENTS];
    sum_counters[COUNTER_SINGULAR_COMPONENTS] =
        max_counters[COUNTER_SINGULAR_COMPONENTS];
    sum_counters[COUNTER_TOTAL_ZERO] = max_counters[COUNTER_TOTAL_ZERO];
    sum_counters[COUNTER_GLOBAL_FACTORIZATIONS] =
        max_counters[COUNTER_GLOBAL_FACTORIZATIONS];
  }
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
  printf(" recursion=%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64,
         counter[COUNTER_RECURSION_0], counter[COUNTER_RECURSION_1],
         counter[COUNTER_RECURSION_2], counter[COUNTER_RECURSION_3]);
  printf(" unique=%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64,
         counter[COUNTER_UNIQUE_0], counter[COUNTER_UNIQUE_1],
         counter[COUNTER_UNIQUE_2], counter[COUNTER_UNIQUE_3]);
  printf(" memo_hits=%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64,
         counter[COUNTER_MEMO_HITS_0], counter[COUNTER_MEMO_HITS_1],
         counter[COUNTER_MEMO_HITS_2], counter[COUNTER_MEMO_HITS_3]);
  printf(" memo_misses=%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64,
         counter[COUNTER_MEMO_MISSES_0], counter[COUNTER_MEMO_MISSES_1],
         counter[COUNTER_MEMO_MISSES_2], counter[COUNTER_MEMO_MISSES_3]);
  printf(" raw_transitions=%" PRIu64 " merged_duplicates=%" PRIu64,
         counter[COUNTER_RAW_TRANSITIONS],
         counter[COUNTER_MERGED_DUPLICATES]);
  printf(" cancelled_zero=%" PRIu64 " terminal_requests=%" PRIu64,
         counter[COUNTER_CANCELLED_ZERO],
         counter[COUNTER_TERMINAL_REQUESTS]);
  printf(" terminal_cache_hits=%" PRIu64,
         counter[COUNTER_TERMINAL_CACHE_HITS]);
  printf(" regular_logical=%" PRIu64 " near_pivot_logical=%" PRIu64,
         counter[COUNTER_REGULAR_COMPONENTS],
         counter[COUNTER_NEAR_PIVOT_COMPONENTS]);
  printf(" singular_logical=%" PRIu64 " total_zero_logical=%" PRIu64,
         counter[COUNTER_SINGULAR_COMPONENTS], counter[COUNTER_TOTAL_ZERO]);
  printf(" local_factorizations=%" PRIu64,
         counter[COUNTER_LOCAL_FACTORIZATIONS]);
  printf(" global_factorizations_logical=%" PRIu64,
         counter[COUNTER_GLOBAL_FACTORIZATIONS]);
  printf(" frontier_peak=%" PRIu64 " memo_peak=%" PRIu64,
         counter[COUNTER_FRONTIER_PEAK], counter[COUNTER_MEMO_PEAK]);
  printf(" krylov_workspace_bytes=%" PRIu64,
         counter[COUNTER_KRYLOV_WORKSPACE_BYTES]);
  printf(" model_workspace_bytes=%" PRIu64,
         counter[COUNTER_MODEL_WORKSPACE_BYTES]);
  printf(" amplitude_workspace_bytes=%" PRIu64,
         counter[COUNTER_AMPLITUDE_WORKSPACE_BYTES]);
  printf(" rss_bytes=%" PRIu64, counter[COUNTER_RSS_BYTES]);
  printf(" depth_seconds=%.17g,%.17g,%.17g,%.17g",
         time[TIME_DEPTH_0], time[TIME_DEPTH_1], time[TIME_DEPTH_2],
         time[TIME_DEPTH_3]);
  printf(" amplitude_seconds=%.17g connectivity_seconds=%.17g",
         time[TIME_AMPLITUDE], time[TIME_CONNECTIVITY]);
  printf(" total_seconds=%.17g\n", time[TIME_TOTAL]);
}

static int run_profile(int site_count, int qp_total, int requested_samples,
                       int audit) {
  ProfileModel fixture;
  MVMCClassicKrylovModelWorkspace *model_workspace = NULL;
  MVMCClassicKrylovRealAmplitudeWorkspace *amplitude_workspace = NULL;
  MVMCKrylovWorkspace *krylov_workspace = NULL;
  const MVMCKrylovFockModel *model = NULL;
  MVMCKrylovLimits limits;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  unsigned int *masks = NULL;
  size_t mask_count = 0;
  size_t sector_size = 0;
  size_t sample_count = 0;
  double *slater = NULL;
  double complex *weights = NULL;
  int *gutzwiller_indices = NULL;
  const double complex projection_parameter[1] = {-0.27};
  MVMCClassicKrylovAmplitudeLayout layout;
  uint64_t local_counters[PROFILE_COUNTER_COUNT] = {0};
  uint64_t sum_counters[PROFILE_COUNTER_COUNT] = {0};
  uint64_t max_counters[PROFILE_COUNTER_COUNT] = {0};
  double local_times[PROFILE_TIME_COUNT] = {0.0};
  double sum_times[PROFILE_TIME_COUNT] = {0.0};
  double max_times[PROFILE_TIME_COUNT] = {0.0};
  const int rank = world_rank();
  const int size = world_size();
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
    sample_count = requested_samples == 0 ? sector_size :
                                               (size_t)requested_samples;
    if (sample_count == 0 || sample_count > sector_size ||
        !initialize_model(&fixture, site_count) ||
        !initialize_slater(site_count, qp_total, &slater, &weights)) {
      local_ready = 0;
    }
  }
  if (local_ready) {
    gutzwiller_indices =
        (int *)calloc((size_t)site_count, sizeof(*gutzwiller_indices));
    if (gutzwiller_indices == NULL) local_ready = 0;
  }
  if (!collective_all_ready(local_ready)) {
    if (rank == 0) {
      fprintf(stderr, "profile fixture allocation/preflight failed\n");
    }
    goto cleanup;
  }
  layout = amplitude_layout(site_count, qp_total, gutzwiller_indices,
                            projection_parameter);
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
  status = mvmc_classic_krylov_real_amplitude_workspace_create(
      &layout, slater, weights, &amplitude_workspace);
  if (!collective_all_ready(status == MVMC_KRYLOV_STATUS_OK)) {
    if (rank == 0) {
      fprintf(stderr, "amplitude create failed collectively: %s\n",
              mvmc_krylov_status_string(status));
    }
    goto cleanup;
  }
  model = mvmc_classic_krylov_model(model_workspace);
  local_ready = model != NULL;
  if (local_ready &&
      !profile_limits(sector_size, model->term_count, &limits)) {
    local_ready = 0;
  }
  if (local_ready) {
    krylov_workspace =
        mvmc_krylov_workspace_create((size_t)site_count, &limits, &status);
    local_ready = krylov_workspace != NULL &&
                  status == MVMC_KRYLOV_STATUS_OK;
  }
  if (!collective_all_ready(local_ready)) {
    if (rank == 0) {
      fprintf(stderr, "Krylov workspace create failed collectively: %s\n",
              mvmc_krylov_status_string(status));
    }
    goto cleanup;
  }
  local_counters[COUNTER_KRYLOV_WORKSPACE_BYTES] =
      (uint64_t)mvmc_krylov_workspace_bytes(krylov_workspace);
  local_counters[COUNTER_MODEL_WORKSPACE_BYTES] =
      (uint64_t)mvmc_classic_krylov_model_workspace_bytes(model_workspace);
  local_counters[COUNTER_AMPLITUDE_WORKSPACE_BYTES] =
      (uint64_t)mvmc_classic_krylov_real_amplitude_workspace_bytes(
          amplitude_workspace);
  if (rank == 0) {
    printf("PROFILE schema=%d fixture=p3_scaling_electronic_real",
           PROFILE_SCHEMA_VERSION);
    printf(" site_count=%d qp_total=%d sample_count=%zu sector_size=%zu",
           site_count, qp_total, sample_count, sector_size);
    printf(" rank_count=%d audit=%d\n", size, audit);
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
    MVMCKrylovResult result;
    const double evaluation_start = monotonic_seconds();
    int order;
    status = mvmc_krylov_evaluate(
        krylov_workspace, model, &configuration, 1,
        mvmc_classic_krylov_real_amplitude, amplitude_workspace, &result);
    local_ready = status == MVMC_KRYLOV_STATUS_OK && result.valid &&
                  result.evaluated_order == MVMC_KRYLOV_MAX_ORDER;
    if (!collective_all_ready(local_ready)) {
      if (rank == 0) {
        fprintf(stderr,
                "profile state=%" PRIu64 " failed collectively: %s\n",
                configuration, mvmc_krylov_status_string(status));
      }
      goto cleanup;
    }
    local_times[TIME_TOTAL] += elapsed_seconds(evaluation_start);
    add_statistics(local_counters, local_times, &result.statistics);
    if (rank == 0) {
      printf("ROW %zu %zu %" PRIu64, sample, sector_index, configuration);
      for (order = 0; order < PROFILE_DEPTH_COUNT; ++order) {
        printf(" %.17g %.17g", creal(result.value[order]),
               cimag(result.value[order]));
      }
      putchar('\n');
      if (audit) {
        printf("STAT %zu", sample);
        printf(" depth_seconds=%.17g,%.17g,%.17g,%.17g",
               result.statistics.depth_wall_seconds[0],
               result.statistics.depth_wall_seconds[1],
               result.statistics.depth_wall_seconds[2],
               result.statistics.depth_wall_seconds[3]);
        printf(" terminal_requests=%" PRIu64,
               result.statistics.terminal_amplitude_requests);
        printf(" raw_transitions=%" PRIu64,
               result.statistics.raw_transitions);
        printf(" memo_peak=%zu workspace_bytes=%zu\n",
               result.statistics.memo_entries_peak,
               result.statistics.workspace_bytes);
      }
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
  mvmc_krylov_workspace_destroy(krylov_workspace);
  mvmc_classic_krylov_real_amplitude_workspace_destroy(amplitude_workspace);
  mvmc_classic_krylov_model_workspace_destroy(model_workspace);
  free(gutzwiller_indices);
  free(weights);
  free(slater);
  destroy_model(&fixture);
  free(masks);
  return return_code;
}

int main(int argc, char **argv) {
  int site_count;
  int qp_total;
  int sample_count;
  int audit;
  int result;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#endif
  if (argc != 5 ||
      !parse_int(argv[1], 4, PROFILE_MAX_SITE_COUNT, &site_count) ||
      (site_count & 1) != 0 || !parse_int(argv[2], 1, 4, &qp_total) ||
      (qp_total != 1 && qp_total != 4) ||
      !parse_int(argv[3], 0, 100000000, &sample_count) ||
      !parse_int(argv[4], 0, 1, &audit)) {
    if (world_rank() == 0) {
      fprintf(stderr,
              "usage: %s EVEN_SITE_COUNT(4..16) QP_TOTAL(1|4) "
              "SAMPLES(0=all) AUDIT(0|1)\n",
              argv[0]);
    }
    result = EXIT_FAILURE;
  } else {
    result = run_profile(site_count, qp_total, sample_count, audit) == 0
                 ? EXIT_SUCCESS
                 : EXIT_FAILURE;
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return result;
}
