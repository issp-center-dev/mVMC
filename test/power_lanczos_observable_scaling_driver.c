/*
 * Testing-only P6-C2 operator-count scaling driver.  This executable calls
 * the production observable action/evaluator/block implementations directly;
 * it does not relax the corrected production selector gate.
 */

#define MVMC_BOUNDED_KRYLOV_PROFILE_DRIVER_NO_MAIN
#include "bounded_krylov_profile_driver.c"

#include "power_lanczos_observable_blocks.h"
#include "power_lanczos_observable_census.h"
#include "power_lanczos_observable_evaluator.h"

#include <float.h>
#include <limits.h>

enum {
  SCALING_SCHEMA_VERSION = 2,
  SCALING_NSITE = 16,
  SCALING_NQP = 8,
  SCALING_MAX_SOURCES = 64,
  SCALING_HASH_HEX = 64
};

typedef struct {
  const char *id;
  const char *model;
  const char *arithmetic;
  int pure_spin;
  int required_family;
  int low_reuse;
} ScalingStratum;

typedef struct {
  MVMCClassicKrylovRawModel raw;
  MVMCClassicKrylovPairCoupling exchange[SCALING_NSITE];
} ScalingPureModel;

typedef struct {
  double setup_seconds;
  double active_seconds;
  uint64_t peak_rss_bytes;
  uint64_t exact_allocated_bytes;
  uint64_t semantic_hash;
  size_t minimum_unique_targets;
  size_t maximum_unique_targets;
} ScalingReducedMetrics;

static const ScalingStratum ScalingStrata[] = {
    {"SC-ONEBODY", "electronic", "real", 0,
     MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, 0},
    {"SC-QUARTIC", "pure_spin", "real", 1,
     MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX, 0},
    {"SC-MIXED", "electronic", "complex", 0, -1, 0},
    {"SC-LOW-REUSE", "electronic", "complex", 0, -1, 1},
};

static double ScalingWallSeconds(void) {
  struct timespec now;
  if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return -1.0;
  return (double)now.tv_sec + 1.0e-9 * (double)now.tv_nsec;
}

static double ScalingElapsedSeconds(double start) {
  const double end = ScalingWallSeconds();
  return start >= 0.0 && end >= start ? end - start : -1.0;
}

static const ScalingStratum *FindStratum(const char *id) {
  size_t index;
  for (index = 0; index < sizeof(ScalingStrata) / sizeof(ScalingStrata[0]);
       ++index) {
    if (strcmp(ScalingStrata[index].id, id) == 0) {
      return &ScalingStrata[index];
    }
  }
  return NULL;
}

static int IsLowerHexSha256(const char *text) {
  size_t index;
  if (text == NULL || strlen(text) != SCALING_HASH_HEX) return 0;
  for (index = 0; index < SCALING_HASH_HEX; ++index) {
    if (!((text[index] >= '0' && text[index] <= '9') ||
          (text[index] >= 'a' && text[index] <= 'f'))) {
      return 0;
    }
  }
  return 1;
}

static uint64_t SeedPrefix(const char *sha256) {
  uint64_t value = 0;
  int index;
  for (index = 0; index < 16; ++index) {
    const char byte = sha256[index];
    value <<= 4;
    value |= (uint64_t)(byte <= '9' ? byte - '0' : byte - 'a' + 10);
  }
  return value == 0 ? UINT64_C(1) : value;
}

static int InitializePureModel(ScalingPureModel *fixture) {
  int bond;
  if (fixture == NULL) return 0;
  memset(fixture, 0, sizeof(*fixture));
  for (bond = 0; bond < SCALING_NSITE; ++bond) {
    fixture->exchange[bond] = (MVMCClassicKrylovPairCoupling){
        bond, (bond + 1) % SCALING_NSITE,
        0.5 + 0.015625 * (double)(bond % 7)};
  }
  fixture->raw.site_count = SCALING_NSITE;
  fixture->raw.up_electron_count = SCALING_NSITE / 2;
  fixture->raw.down_electron_count = SCALING_NSITE / 2;
  fixture->raw.pure_spin = 1;
  fixture->raw.exchange_count = SCALING_NSITE;
  fixture->raw.exchange = fixture->exchange;
  return 1;
}

static int BuildSourceBatch(const ScalingStratum *stratum,
                            const char *timing_seed_sha256,
                            size_t source_count, uint64_t *sources,
                            unsigned int **all_masks,
                            size_t *all_mask_count) {
  size_t source;
  uint64_t seed;
  if (stratum == NULL || sources == NULL || all_masks == NULL ||
      all_mask_count == NULL || source_count == 0 ||
      source_count > SCALING_MAX_SOURCES) {
    return 0;
  }
  *all_masks = NULL;
  *all_mask_count = 0;
  seed = SeedPrefix(timing_seed_sha256);
  if (stratum->low_reuse) {
    unsigned int reservoir[20];
    size_t reservoir_count = 0;
    unsigned int mask;
    const unsigned int fixed = UINT32_C(0x001f);
    for (mask = 0; mask < 64; ++mask) {
      if (popcount_u32(mask) == 3) reservoir[reservoir_count++] = mask;
    }
    if (reservoir_count != 20) return 0;
    for (source = 0; source < source_count; ++source) {
      const size_t ordinal = (source + (size_t)(seed % 400)) % 400;
      const unsigned int up =
          fixed | (reservoir[ordinal / 20] << 10);
      const unsigned int down =
          fixed | (reservoir[ordinal % 20] << 10);
      sources[source] = (uint64_t)up | ((uint64_t)down << SCALING_NSITE);
    }
  } else {
    size_t mask_count = fixed_masks(SCALING_NSITE, all_masks);
    if (mask_count == 0) return 0;
    *all_mask_count = mask_count;
    for (source = 0; source < source_count; ++source) {
      const size_t up_index =
          ((size_t)(seed % mask_count) + source * 131) % mask_count;
      size_t down_index;
      unsigned int down;
      if (stratum->pure_spin) {
        down = (~(*all_masks)[up_index]) & UINT32_C(0xffff);
      } else {
        down_index =
            ((size_t)((seed >> 17) % mask_count) + source * 313) %
            mask_count;
        down = (*all_masks)[down_index];
      }
      sources[source] = (uint64_t)(*all_masks)[up_index] |
                        ((uint64_t)down << SCALING_NSITE);
    }
  }
  for (source = 0; source < source_count; ++source) {
    size_t prior;
    for (prior = 0; prior < source; ++prior) {
      if (sources[prior] == sources[source]) return 0;
    }
  }
  return 1;
}

static MVMCPowerLanczosObservableCensusStatus BuildObservablePlan(
    const char *one_path, const char *quartic_ex_path,
    const char *quartic_path, MVMCPowerLanczosObservablePlan *plan,
    char *diagnostic, size_t diagnostic_capacity) {
  const char *paths[MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT] = {
      strcmp(one_path, "-") == 0 ? NULL : one_path,
      strcmp(quartic_ex_path, "-") == 0 ? NULL : quartic_ex_path,
      strcmp(quartic_path, "-") == 0 ? NULL : quartic_path};
  MVMCPowerLanczosObservableCensusStatus status =
      MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
#ifdef _mpi_use
  unsigned char *wire = NULL;
  size_t wire_size = 0;
  uint64_t wire_size_u64 = 0;
  int status_wire;
  if (world_rank() == 0) {
    status = mvmc_power_lanczos_observable_plan_build_from_files(
        SCALING_NSITE, SCALING_NSITE, paths, plan, diagnostic,
        diagnostic_capacity);
  }
  status_wire = (int)status;
  MPI_Bcast(&status_wire, 1, MPI_INT, 0, MPI_COMM_WORLD);
  status = (MVMCPowerLanczosObservableCensusStatus)status_wire;
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  if (world_rank() == 0) {
    status = mvmc_power_lanczos_observable_plan_wire_size(plan, &wire_size);
    if (status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
        wire_size > (size_t)INT_MAX) {
      status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_RESOURCE_LIMIT;
    }
    if (status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) {
      wire = (unsigned char *)malloc(wire_size);
      if (wire == NULL) {
        status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ALLOCATION_FAILURE;
      }
    }
    if (status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) {
      size_t packed = 0;
      status = mvmc_power_lanczos_observable_plan_pack(
          plan, wire, wire_size, &packed);
      if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK ||
          packed != wire_size) {
        status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
      }
    }
    wire_size_u64 = (uint64_t)wire_size;
  }
  status_wire = (int)status;
  MPI_Bcast(&status_wire, 1, MPI_INT, 0, MPI_COMM_WORLD);
  status = (MVMCPowerLanczosObservableCensusStatus)status_wire;
  if (status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) {
    MPI_Bcast(&wire_size_u64, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    if (wire_size_u64 == 0 || wire_size_u64 > (uint64_t)SIZE_MAX ||
        wire_size_u64 > (uint64_t)INT_MAX) {
      status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_RESOURCE_LIMIT;
    } else if (world_rank() != 0) {
      wire_size = (size_t)wire_size_u64;
      wire = (unsigned char *)malloc(wire_size);
      if (wire == NULL) {
        status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ALLOCATION_FAILURE;
      }
    }
  }
  if (!collective_all_ready(
          status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK)) {
    free(wire);
    return status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK
               ? status
               : MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ALLOCATION_FAILURE;
  }
  MPI_Bcast(wire, (int)wire_size, MPI_BYTE, 0, MPI_COMM_WORLD);
  if (world_rank() != 0) {
    status = mvmc_power_lanczos_observable_plan_unpack(
        wire, wire_size, plan, diagnostic, diagnostic_capacity);
  }
  free(wire);
  return status;
#else
  status = mvmc_power_lanczos_observable_plan_build_from_files(
      SCALING_NSITE, SCALING_NSITE, paths, plan, diagnostic,
      diagnostic_capacity);
  return status;
#endif
}

static int ValidatePlanProfile(const ScalingStratum *stratum,
                               const MVMCPowerLanczosObservablePlan *plan,
                               size_t operator_count) {
  int family;
  if (stratum == NULL || plan == NULL || plan->record_count < 0 ||
      (size_t)plan->record_count != operator_count) {
    return 0;
  }
  if (stratum->required_family >= 0) {
    for (family = 0; family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT;
         ++family) {
      const int expected = family == stratum->required_family
                               ? (int)operator_count
                               : 0;
      if (plan->family_count[family] != expected) return 0;
    }
  } else if (operator_count >= 3 &&
             (plan->family_count[0] == 0 || plan->family_count[1] == 0 ||
              plan->family_count[2] == 0)) {
    return 0;
  }
  return 1;
}

static uint64_t FnvBytes(uint64_t hash, const void *data, size_t size) {
  const unsigned char *bytes = (const unsigned char *)data;
  size_t index;
  for (index = 0; index < size; ++index) {
    hash ^= bytes[index];
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

static uint64_t HashPayload(const double complex *coefficient,
                            size_t coefficient_count,
                            const double complex *final_values,
                            size_t final_count,
                            const uint64_t *block_counts,
                            size_t block_count) {
  uint64_t hash = UINT64_C(1469598103934665603);
  hash = FnvBytes(hash, coefficient,
                  coefficient_count * sizeof(coefficient[0]));
  hash = FnvBytes(hash, final_values,
                  final_count * sizeof(final_values[0]));
  return FnvBytes(hash, block_counts,
                  block_count * sizeof(block_counts[0]));
}

static int WriteU64Little(FILE *stream, uint64_t value) {
  unsigned char bytes[8];
  int index;
  for (index = 0; index < 8; ++index) {
    bytes[index] = (unsigned char)(value >> (8 * index));
  }
  return fwrite(bytes, 1, sizeof(bytes), stream) == sizeof(bytes);
}

static int WriteComplexPayload(const char *path,
                               const double complex *values,
                               size_t value_count) {
  FILE *stream = fopen(path, "wb");
  size_t index;
  if (stream == NULL) return 0;
  for (index = 0; index < value_count; ++index) {
    uint64_t real_bits;
    uint64_t imaginary_bits;
    const double real_value = creal(values[index]);
    const double imaginary_value = cimag(values[index]);
    memcpy(&real_bits, &real_value, sizeof(real_bits));
    memcpy(&imaginary_bits, &imaginary_value, sizeof(imaginary_bits));
    if (!WriteU64Little(stream, real_bits) ||
        !WriteU64Little(stream, imaginary_bits)) {
      fclose(stream);
      return 0;
    }
  }
  return fclose(stream) == 0;
}

static int WriteMetrics(
    const char *path, const ScalingStratum *stratum,
    const MVMCPowerLanczosObservablePlan *observable_plan,
    size_t operator_count, size_t source_count, int replicate,
    const char *run_kind, const char *timing_seed_sha256,
    size_t block_count, size_t coefficient_payload_bytes,
    size_t final_payload_bytes, uint64_t callback_calls,
    const ScalingReducedMetrics *metrics) {
  FILE *stream = fopen(path, "wb");
  int result;
  if (stream == NULL) return 0;
  result = fprintf(
      stream,
      "{\n"
      "  \"schema_version\": %d,\n"
      "  \"timer_id\": \"p6c2_observable_source_batch_v2\",\n"
      "  \"stratum_id\": \"%s\",\n"
      "  \"operator_count\": %zu,\n"
      "  \"saved_source_count\": %zu,\n"
      "  \"replicate_ordinal\": %d,\n"
      "  \"run_kind\": \"%s\",\n"
      "  \"model\": \"%s\",\n"
      "  \"arithmetic\": \"%s\",\n"
      "  \"Nsite\": %d,\n"
      "  \"NQPFull\": %d,\n"
      "  \"family_counts\": {\"cisajs\": %d, \"cisajscktaltex\": %d, \"cisajscktalt\": %d},\n"
      "  \"raw_observable_census_sha256\": \"%s\",\n"
      "  \"operator_census_sha256\": \"%s\",\n"
      "  \"timing_seed_sha256\": \"%s\",\n"
      "  \"mpi_world_size\": %d,\n"
      "  \"setup_seconds\": %.17g,\n"
      "  \"active_seconds\": %.17g,\n"
      "  \"wall_seconds\": %.17g,\n"
      "  \"peak_RSS_bytes\": %" PRIu64 ",\n"
      "  \"exact_allocated_bytes\": %" PRIu64 ",\n"
      "  \"minimum_unique_target_count\": %zu,\n"
      "  \"maximum_unique_target_count\": %zu,\n"
      "  \"block_count\": %zu,\n"
      "  \"coefficient_payload_bytes\": %zu,\n"
      "  \"final_payload_bytes\": %zu,\n"
      "  \"artifact_payload_bytes\": %zu,\n"
      "  \"artifact_files\": 2,\n"
      "  \"amplitude_callback_calls\": %" PRIu64 ",\n"
      "  \"normalized_semantic_payload_fnv1a64\": \"%016" PRIx64 "\",\n"
      "  \"status\": \"PASS\"\n"
      "}\n",
      SCALING_SCHEMA_VERSION, stratum->id, operator_count, source_count,
      replicate, run_kind, stratum->model, stratum->arithmetic,
      SCALING_NSITE, SCALING_NQP, observable_plan->family_count[0],
      observable_plan->family_count[1], observable_plan->family_count[2],
      observable_plan->raw_observable_census_sha256,
      observable_plan->semantic_observable_census_sha256,
      timing_seed_sha256, world_size(), metrics->setup_seconds,
      metrics->active_seconds,
      metrics->setup_seconds + metrics->active_seconds,
      metrics->peak_rss_bytes, metrics->exact_allocated_bytes,
      metrics->minimum_unique_targets, metrics->maximum_unique_targets,
      block_count, coefficient_payload_bytes, final_payload_bytes,
      coefficient_payload_bytes + final_payload_bytes, callback_calls,
      metrics->semantic_hash);
  return result > 0 && fclose(stream) == 0;
}

static int AddSizeToU64(uint64_t *total, size_t value) {
  if (total == NULL || (uint64_t)value > UINT64_MAX - *total) return 0;
  *total += (uint64_t)value;
  return 1;
}

static int RunScaling(
    const ScalingStratum *stratum, size_t operator_count,
    size_t source_count, int replicate, const char *run_kind,
    const char *timing_seed_sha256, const char *one_path,
    const char *quartic_ex_path, const char *quartic_path,
    const char *coefficient_path, const char *final_path,
    const char *metrics_path) {
  const size_t block_length = source_count >= 32 ? source_count / 32 : 1;
  const size_t block_count = source_count / block_length;
  const size_t coefficient_entries = 4 * operator_count;
  const double scales[2] = {0.0, 0.0};
  const double inverse_sqrt_two = 1.0 / sqrt(2.0);
  double complex alpha[2] = {inverse_sqrt_two, inverse_sqrt_two};
  ProfileModel electronic_fixture;
  ScalingPureModel pure_fixture;
  MVMCClassicKrylovModelWorkspace *model_workspace = NULL;
  MVMCKrylovBoundedPlan *bounded_plan = NULL;
  MVMCKrylovBoundedWorkspace *bounded_workspace = NULL;
  MVMCKrylovBoundedCollectiveWorkspace *collective_workspace = NULL;
  ProfileScaledAmplitude *amplitude_workspace = NULL;
  MVMCPowerLanczosObservableEvaluatorWorkspace *evaluator = NULL;
  MVMCPowerLanczosObservableBlockAccumulator *coefficient_blocks = NULL;
  MVMCPowerLanczosObservableBlockAccumulator *final_blocks = NULL;
  MVMCPowerLanczosObservableBlockSummary coefficient_summary;
  MVMCPowerLanczosObservableBlockSummary final_summary;
  MVMCPowerLanczosObservablePlan observable_plan;
  MVMCKrylovBoundedLimits limits;
  const MVMCKrylovFockModel *model = NULL;
  uint64_t sources[SCALING_MAX_SOURCES];
  size_t coefficient_unique_targets[SCALING_MAX_SOURCES];
  unsigned int *all_masks = NULL;
  size_t all_mask_count = 0;
  double *slater = NULL;
  double complex *weights = NULL;
  double complex *coefficient_sample = NULL;
  double complex *final_sample = NULL;
  double complex *coefficient_export = NULL;
  double complex *final_export = NULL;
  uint64_t *coefficient_counts = NULL;
  uint64_t *final_counts = NULL;
  char diagnostic[512] = {0};
  MVMCPowerLanczosObservableCensusStatus census_status;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  ScalingReducedMetrics local_metrics;
  ScalingReducedMetrics reduced_metrics;
  uint64_t generation;
  uint64_t exact_allocated = 0;
  uint64_t callback_calls = 0;
  uint64_t reduced_callback_calls = 0;
  size_t source;
  size_t coefficient_payload_bytes = 0;
  size_t final_payload_bytes = 0;
  double setup_start = ScalingWallSeconds();
  double active_start;
  int local_ready = 1;
  int result = 1;

  memset(&electronic_fixture, 0, sizeof(electronic_fixture));
  memset(&pure_fixture, 0, sizeof(pure_fixture));
  memset(&observable_plan, 0, sizeof(observable_plan));
  memset(&coefficient_summary, 0, sizeof(coefficient_summary));
  memset(&final_summary, 0, sizeof(final_summary));
  memset(&local_metrics, 0, sizeof(local_metrics));
  memset(&reduced_metrics, 0, sizeof(reduced_metrics));
  memset(coefficient_unique_targets, 0, sizeof(coefficient_unique_targets));
  local_metrics.minimum_unique_targets = operator_count;
  if (strcmp(stratum->arithmetic, "complex") == 0) {
    alpha[1] *= 1.0 + 0.125 * I;
  }

  mvmc_power_lanczos_observable_plan_init(&observable_plan);
  census_status = BuildObservablePlan(
      one_path, quartic_ex_path, quartic_path, &observable_plan,
      diagnostic, sizeof(diagnostic));
  local_ready = census_status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
                ValidatePlanProfile(stratum, &observable_plan,
                                    operator_count) &&
                BuildSourceBatch(stratum, timing_seed_sha256, source_count,
                                 sources, &all_masks, &all_mask_count);
  if (!collective_all_ready(local_ready)) {
    if (world_rank() == 0) {
      fprintf(stderr, "P6-C2 scaling input failed: census=%s detail=%s\n",
              mvmc_power_lanczos_observable_census_status_string(
                  census_status),
              diagnostic);
    }
    goto cleanup;
  }
  if (stratum->pure_spin) {
    local_ready = InitializePureModel(&pure_fixture);
  } else {
    local_ready = initialize_model(&electronic_fixture, SCALING_NSITE);
  }
  if (!collective_all_ready(local_ready)) goto cleanup;
  status = mvmc_classic_krylov_model_workspace_create_from_root(
      world_rank() == 0
          ? (stratum->pure_spin ? &pure_fixture.raw
                                : &electronic_fixture.raw)
          : NULL,
      world_communicator(), &model_workspace);
  local_ready = status == MVMC_KRYLOV_STATUS_OK && model_workspace != NULL;
  if (!collective_all_ready(local_ready)) goto cleanup;
  model = mvmc_classic_krylov_model(model_workspace);
  local_ready = model != NULL && profile_limits(64 * 1024 * 1024, model,
                                                &limits);
  if (local_ready) {
    limits.max_order = 1;
    status = mvmc_bounded_krylov_plan_create(model, &limits, &bounded_plan);
    local_ready = status == MVMC_KRYLOV_STATUS_OK && bounded_plan != NULL;
  }
  if (local_ready) {
    status = mvmc_bounded_krylov_workspace_create(bounded_plan,
                                                   &bounded_workspace);
    local_ready = status == MVMC_KRYLOV_STATUS_OK &&
                  bounded_workspace != NULL;
  }
  if (local_ready) {
    status = mvmc_bounded_krylov_collective_workspace_create(
        world_communicator(), SCALING_NQP, 2, &collective_workspace);
    local_ready = status == MVMC_KRYLOV_STATUS_OK &&
                  collective_workspace != NULL;
  }
  if (local_ready) {
    local_ready = initialize_slater(SCALING_NSITE, SCALING_NQP, &slater,
                                    &weights);
  }
  if (local_ready && strcmp(stratum->arithmetic, "complex") == 0) {
    int qp;
    for (qp = 0; qp < SCALING_NQP; ++qp) {
      const double angle = 0.125 * (double)(qp + 1);
      weights[qp] *= cos(angle) + I * sin(angle);
    }
  }
  if (local_ready) {
    local_ready = create_profile_amplitude(
        SCALING_NSITE, SCALING_NQP, slater, weights, -0.27,
        collective_workspace, &amplitude_workspace);
  }
  if (local_ready) {
    status = mvmc_power_lanczos_observable_evaluator_workspace_create(
        operator_count, 1, &evaluator);
    local_ready = status == MVMC_KRYLOV_STATUS_OK && evaluator != NULL;
  }
  if (local_ready) {
    status = mvmc_power_lanczos_observable_block_create(
        MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT, operator_count,
        block_count, block_length, &coefficient_blocks);
    local_ready = status == MVMC_KRYLOV_STATUS_OK;
  }
  if (local_ready) {
    status = mvmc_power_lanczos_observable_block_create(
        MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_FINAL, operator_count,
        block_count, block_length, &final_blocks);
    local_ready = status == MVMC_KRYLOV_STATUS_OK;
  }
  coefficient_sample =
      (double complex *)calloc(coefficient_entries,
                               sizeof(coefficient_sample[0]));
  final_sample =
      (double complex *)calloc(operator_count, sizeof(final_sample[0]));
  coefficient_export = (double complex *)calloc(
      block_count * coefficient_entries, sizeof(coefficient_export[0]));
  final_export = (double complex *)calloc(
      block_count * operator_count, sizeof(final_export[0]));
  coefficient_counts =
      (uint64_t *)calloc(block_count, sizeof(coefficient_counts[0]));
  final_counts = (uint64_t *)calloc(block_count, sizeof(final_counts[0]));
  local_ready = local_ready && coefficient_sample != NULL &&
                final_sample != NULL && coefficient_export != NULL &&
                final_export != NULL && coefficient_counts != NULL &&
                final_counts != NULL;
  if (!collective_all_ready(local_ready)) goto cleanup;

  status = mvmc_power_lanczos_observable_block_summary(
      coefficient_blocks, &coefficient_summary);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_block_summary(final_blocks,
                                                         &final_summary);
  }
  local_ready = status == MVMC_KRYLOV_STATUS_OK;
  if (!collective_all_ready(local_ready)) goto cleanup;
  if (!AddSizeToU64(&exact_allocated, sizeof(observable_plan)) ||
      !AddSizeToU64(&exact_allocated,
                    operator_count * sizeof(observable_plan.records[0])) ||
      !AddSizeToU64(&exact_allocated,
                    mvmc_classic_krylov_model_workspace_bytes(
                        model_workspace)) ||
      !AddSizeToU64(&exact_allocated,
                    mvmc_bounded_krylov_plan_bytes(bounded_plan)) ||
      !AddSizeToU64(&exact_allocated,
                    mvmc_bounded_krylov_workspace_bytes(
                        bounded_workspace)) ||
      !AddSizeToU64(&exact_allocated,
                    mvmc_bounded_krylov_collective_workspace_bytes(
                        collective_workspace)) ||
      !AddSizeToU64(&exact_allocated,
                    profile_amplitude_workspace_bytes(amplitude_workspace)) ||
      !AddSizeToU64(&exact_allocated,
                    mvmc_power_lanczos_observable_evaluator_workspace_bytes(
                        evaluator)) ||
      !AddSizeToU64(&exact_allocated, coefficient_summary.allocated_bytes) ||
      !AddSizeToU64(&exact_allocated, final_summary.allocated_bytes) ||
      !AddSizeToU64(&exact_allocated,
                    (coefficient_entries + operator_count +
                     block_count * coefficient_entries +
                     block_count * operator_count) *
                        sizeof(double complex)) ||
      !AddSizeToU64(&exact_allocated,
                    2 * block_count * sizeof(uint64_t)) ||
      !AddSizeToU64(&exact_allocated, all_mask_count * sizeof(unsigned int)) ||
      !AddSizeToU64(&exact_allocated, sizeof(sources)) ||
      !AddSizeToU64(&exact_allocated,
                    sizeof(coefficient_unique_targets))) {
    goto cleanup;
  }
  generation = SeedPrefix(timing_seed_sha256) ^
               mvmc_bounded_krylov_plan_hash(bounded_plan) ^
               UINT64_C(0x503643325343414c);
  if (generation == 0) generation = UINT64_C(1);
  local_metrics.setup_seconds = ScalingElapsedSeconds(setup_start);

  active_start = ScalingWallSeconds();
  status = mvmc_bounded_krylov_session_begin(
      bounded_workspace, profile_scaled_amplitude, amplitude_workspace,
      generation);
  for (source = 0; status == MVMC_KRYLOV_STATUS_OK && source < source_count;
       ++source) {
    MVMCPowerLanczosObservableLayout layout = {
        SCALING_NSITE, SCALING_NSITE / 2, SCALING_NSITE / 2,
        stratum->pure_spin};
    MVMCPowerLanczosObservableSampleDiagnostics diagnostics;
    status = mvmc_power_lanczos_observable_coefficient_sample(
        evaluator, bounded_workspace, &layout, &observable_plan,
        &sources[source], 1, 0.0, scales, coefficient_sample,
        coefficient_entries, &diagnostics);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      if (!diagnostics.valid ||
          diagnostics.unique_target_count > operator_count) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      } else {
        coefficient_unique_targets[source] =
            diagnostics.unique_target_count;
        if (diagnostics.unique_target_count <
            local_metrics.minimum_unique_targets) {
          local_metrics.minimum_unique_targets =
              diagnostics.unique_target_count;
        }
        if (diagnostics.unique_target_count >
            local_metrics.maximum_unique_targets) {
          local_metrics.maximum_unique_targets =
              diagnostics.unique_target_count;
        }
        status = mvmc_power_lanczos_observable_block_add_sample(
            coefficient_blocks, coefficient_sample, coefficient_entries);
      }
    }
  }
  if (mvmc_bounded_krylov_session_is_active(bounded_workspace)) {
    const MVMCKrylovStatus end_status =
        mvmc_bounded_krylov_session_end(bounded_workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) status = end_status;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_session_begin(
        bounded_workspace, profile_scaled_amplitude, amplitude_workspace,
        generation ^ UINT64_C(0x8000000000000000));
  }
  for (source = 0; status == MVMC_KRYLOV_STATUS_OK && source < source_count;
       ++source) {
    MVMCPowerLanczosObservableLayout layout = {
        SCALING_NSITE, SCALING_NSITE / 2, SCALING_NSITE / 2,
        stratum->pure_spin};
    MVMCPowerLanczosObservableSampleDiagnostics diagnostics;
    status = mvmc_power_lanczos_observable_final_sample(
        evaluator, bounded_workspace, &layout, &observable_plan,
        &sources[source], 1, 0.0, scales, alpha, final_sample,
        operator_count, &diagnostics);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      if (!diagnostics.valid ||
          diagnostics.unique_target_count !=
              coefficient_unique_targets[source]) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      } else {
        status = mvmc_power_lanczos_observable_block_add_sample(
            final_blocks, final_sample, operator_count);
      }
    }
  }
  if (mvmc_bounded_krylov_session_is_active(bounded_workspace)) {
    const MVMCKrylovStatus end_status =
        mvmc_bounded_krylov_session_end(bounded_workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) status = end_status;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_block_export(
        coefficient_blocks, coefficient_export,
        block_count * coefficient_entries, coefficient_counts,
        block_count);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_block_export(
        final_blocks, final_export, block_count * operator_count,
        final_counts, block_count);
  }
  local_metrics.active_seconds = ScalingElapsedSeconds(active_start);
  local_ready = status == MVMC_KRYLOV_STATUS_OK;
  if (local_ready) {
    size_t block;
    for (block = 0; block < block_count; ++block) {
      if (coefficient_counts[block] != block_length ||
          final_counts[block] != block_length) {
        local_ready = 0;
      }
    }
  }
  if (local_ready && stratum->low_reuse &&
      10 * local_metrics.minimum_unique_targets < 9 * operator_count) {
    local_ready = 0;
  }
  local_metrics.peak_rss_bytes = resident_set_bytes();
  local_metrics.exact_allocated_bytes = exact_allocated;
  local_metrics.semantic_hash = HashPayload(
      coefficient_export, block_count * coefficient_entries, final_export,
      block_count * operator_count, coefficient_counts, block_count);
  callback_calls = amplitude_workspace->callback_calls;
  if (!collective_all_ready(local_ready)) {
    if (world_rank() == 0) {
      fprintf(stderr,
              "P6-C2 scaling active batch failed: status=%s min_U=%zu R=%zu\n",
              mvmc_krylov_status_string(status),
              local_metrics.minimum_unique_targets, operator_count);
    }
    goto cleanup;
  }

#ifdef _mpi_use
  {
    uint64_t hash_min = 0;
    uint64_t hash_max = 0;
    uint64_t minimum_u = (uint64_t)local_metrics.minimum_unique_targets;
    uint64_t maximum_u = (uint64_t)local_metrics.maximum_unique_targets;
    uint64_t reduced_minimum_u = 0;
    uint64_t reduced_maximum_u = 0;
    MPI_Reduce(&local_metrics.setup_seconds, &reduced_metrics.setup_seconds,
               1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_metrics.active_seconds, &reduced_metrics.active_seconds,
               1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_metrics.peak_rss_bytes,
               &reduced_metrics.peak_rss_bytes, 1, MPI_UINT64_T, MPI_MAX, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&local_metrics.exact_allocated_bytes,
               &reduced_metrics.exact_allocated_bytes, 1, MPI_UINT64_T,
               MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&minimum_u, &reduced_minimum_u, 1, MPI_UINT64_T, MPI_MIN, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&maximum_u, &reduced_maximum_u, 1, MPI_UINT64_T, MPI_MAX, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&local_metrics.semantic_hash, &hash_min, 1, MPI_UINT64_T,
               MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_metrics.semantic_hash, &hash_max, 1, MPI_UINT64_T,
               MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&callback_calls, &reduced_callback_calls, 1, MPI_UINT64_T,
               MPI_SUM, 0, MPI_COMM_WORLD);
    if (world_rank() == 0) {
      reduced_metrics.minimum_unique_targets = (size_t)reduced_minimum_u;
      reduced_metrics.maximum_unique_targets = (size_t)reduced_maximum_u;
      reduced_metrics.semantic_hash = hash_min;
      local_ready = hash_min == hash_max;
    }
    MPI_Bcast(&local_ready, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
#else
  reduced_metrics = local_metrics;
  reduced_callback_calls = callback_calls;
#endif
  if (!local_ready) {
    if (world_rank() == 0) {
      fprintf(stderr, "P6-C2 scaling MPI semantic payload mismatch\n");
    }
    goto cleanup;
  }
  if (mvmc_power_lanczos_observable_block_payload_bytes(
          MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT, operator_count,
          block_count, &coefficient_payload_bytes) !=
          MVMC_KRYLOV_STATUS_OK ||
      mvmc_power_lanczos_observable_block_payload_bytes(
          MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_FINAL, operator_count,
          block_count, &final_payload_bytes) != MVMC_KRYLOV_STATUS_OK) {
    goto cleanup;
  }
  if (world_rank() == 0) {
    local_ready = WriteComplexPayload(
                      coefficient_path, coefficient_export,
                      block_count * coefficient_entries) &&
                  WriteComplexPayload(final_path, final_export,
                                      block_count * operator_count) &&
                  WriteMetrics(
                      metrics_path, stratum, &observable_plan,
                      operator_count, source_count, replicate, run_kind,
                      timing_seed_sha256, block_count,
                      coefficient_payload_bytes, final_payload_bytes,
                      reduced_callback_calls, &reduced_metrics);
  }
#ifdef _mpi_use
  MPI_Bcast(&local_ready, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  if (!local_ready) goto cleanup;
  result = 0;

cleanup:
  if (mvmc_bounded_krylov_session_is_active(bounded_workspace)) {
    (void)mvmc_bounded_krylov_session_end(bounded_workspace);
  }
  free(final_counts);
  free(coefficient_counts);
  free(final_export);
  free(coefficient_export);
  free(final_sample);
  free(coefficient_sample);
  mvmc_power_lanczos_observable_block_destroy(final_blocks);
  mvmc_power_lanczos_observable_block_destroy(coefficient_blocks);
  mvmc_power_lanczos_observable_evaluator_workspace_destroy(evaluator);
  destroy_profile_amplitude(amplitude_workspace);
  free(weights);
  free(slater);
  mvmc_bounded_krylov_collective_workspace_destroy(collective_workspace);
  mvmc_bounded_krylov_workspace_destroy(bounded_workspace);
  mvmc_bounded_krylov_plan_destroy(bounded_plan);
  mvmc_classic_krylov_model_workspace_destroy(model_workspace);
  destroy_model(&electronic_fixture);
  free(all_masks);
  mvmc_power_lanczos_observable_plan_destroy(&observable_plan);
  return result;
}

int main(int argc, char **argv) {
  const ScalingStratum *stratum = NULL;
  size_t operator_count = 0;
  size_t source_count = 0;
  int replicate = 0;
  int result = EXIT_FAILURE;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#endif
  if (argc != 13 || (stratum = FindStratum(argv[1])) == NULL ||
      !parse_size_arg(argv[2], 512, &operator_count) ||
      operator_count == 0 ||
      !parse_size_arg(argv[3], SCALING_MAX_SOURCES, &source_count) ||
      source_count == 0 || !parse_int(argv[4], 0, 6, &replicate) ||
      (strcmp(argv[5], "warmup") != 0 &&
       strcmp(argv[5], "measured") != 0 &&
       strcmp(argv[5], "parity") != 0 &&
       strcmp(argv[5], "selftest") != 0) ||
      !IsLowerHexSha256(argv[6])) {
    if (world_rank() == 0) {
      fprintf(
          stderr,
          "usage: %s STRATUM R SOURCES REPLICATE "
          "KIND(warmup|measured|parity|selftest) TIMING_SEED_SHA256 "
          "ONE_DEF|- QUARTIC_EX_DEF|- QUARTIC_DEF|- "
          "COEFFICIENT_BIN FINAL_BIN METRICS_JSON\n",
          argv[0]);
    }
  } else {
    result = RunScaling(stratum, operator_count, source_count, replicate,
                        argv[5], argv[6], argv[7], argv[8], argv[9],
                        argv[10], argv[11], argv[12]) == 0
                 ? EXIT_SUCCESS
                 : EXIT_FAILURE;
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return result;
}
