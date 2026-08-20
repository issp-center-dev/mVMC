/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#define MVMC_BOUNDED_KRYLOV_PROFILE_DRIVER_NO_MAIN
#include "bounded_krylov_profile_driver.c"

#include "krylov_matrix_measurement.h"
#include "krylov_positive_sampler.h"

#include <float.h>
#include <limits.h>

enum {
  MARKOV_SCHEMA_VERSION = 2,
  MARKOV_SESSION_SCHEMA_VERSION = 4,
  MARKOV_ORDER = MVMC_KRYLOV_MAX_ORDER,
  MARKOV_DIMENSION = MVMC_KRYLOV_MAX_ORDER,
  MARKOV_FAMILY_COUNT = 3,
  MARKOV_OFFICIAL_BLOCKS = MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS,
  MARKOV_DIAGNOSTIC_BLOCKS = MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS
};

typedef enum {
  MARKOV_FAMILY_S = 0,
  MARKOV_FAMILY_K,
  MARKOV_FAMILY_B
} MarkovMatrixFamily;

typedef struct {
  double complex exact_sum;
  double unsigned_sum;
  double complex exact_theta;
  double budget;
} MarkovEntryBudget;

typedef struct {
  MarkovEntryBudget entry[MARKOV_FAMILY_COUNT][6];
  double target_sum;
} MarkovExactBudget;

static void markov_hash_byte(uint64_t *hash, unsigned char value) {
  *hash ^= (uint64_t)value;
  *hash *= UINT64_C(1099511628211);
}

static void markov_hash_u64(uint64_t *hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    markov_hash_byte(hash, (unsigned char)(value & UINT64_C(0xff)));
    value >>= 8;
  }
}

static uint64_t markov_double_bits(double value) {
  uint64_t bits = 0;
  if (sizeof(bits) == sizeof(value)) memcpy(&bits, &value, sizeof(bits));
  return bits;
}

static int parse_markov_double(const char *text, double minimum,
                               double maximum, double *value) {
  char *end = NULL;
  double parsed;
  if (text == NULL || value == NULL || *text == '\0') return 0;
  errno = 0;
  parsed = strtod(text, &end);
  if (errno != 0 || end == text || *end != '\0' || !isfinite(parsed) ||
      parsed < minimum || parsed > maximum) {
    return 0;
  }
  *value = parsed;
  return 1;
}

static int parse_markov_u64(const char *text, uint64_t *value) {
  char *end = NULL;
  unsigned long long parsed;
  if (text == NULL || value == NULL || *text == '\0') return 0;
  errno = 0;
  parsed = strtoull(text, &end, 0);
  if (errno != 0 || end == text || *end != '\0') return 0;
  *value = (uint64_t)parsed;
  return 1;
}

static int markov_reduce_max_times(const double *local, double *maximum,
                                   size_t count) {
  if (local == NULL || maximum == NULL || count == 0 || count > INT_MAX) {
    return 0;
  }
#ifdef _mpi_use
  return MPI_Reduce(local, maximum, (int)count, MPI_DOUBLE, MPI_MAX, 0,
                    world_communicator()) == MPI_SUCCESS;
#else
  memcpy(maximum, local, count * sizeof(local[0]));
  return 1;
#endif
}

static const char *markov_family_name(MarkovMatrixFamily family) {
  switch (family) {
    case MARKOV_FAMILY_S:
      return "S";
    case MARKOV_FAMILY_K:
      return "K";
    case MARKOV_FAMILY_B:
      return "B";
  }
  return "unknown";
}

static double complex *markov_family_values(
    MarkovMatrixFamily family, double complex *overlap,
    double complex *hamiltonian, double complex *hamiltonian_squared) {
  switch (family) {
    case MARKOV_FAMILY_S:
      return overlap;
    case MARKOV_FAMILY_K:
      return hamiltonian;
    case MARKOV_FAMILY_B:
      return hamiltonian_squared;
  }
  return NULL;
}

static MVMCKrylovMatrixKind markov_matrix_kind(MarkovMatrixFamily family) {
  switch (family) {
    case MARKOV_FAMILY_S:
      return MVMC_KRYLOV_MATRIX_OVERLAP;
    case MARKOV_FAMILY_K:
      return MVMC_KRYLOV_MATRIX_HAMILTONIAN;
    case MARKOV_FAMILY_B:
      return MVMC_KRYLOV_MATRIX_HAMILTONIAN_SQUARED;
  }
  return MVMC_KRYLOV_MATRIX_OVERLAP;
}

static size_t markov_series_index(size_t sample_count,
                                  MarkovMatrixFamily family, size_t entry,
                                  size_t upper_count, size_t sample) {
  return (((size_t)family * upper_count + entry) * sample_count) + sample;
}

static int markov_configuration_at(size_t sector_index,
                                   const unsigned int *masks,
                                   size_t mask_count,
                                   int site_count,
                                   uint64_t *configuration) {
  const size_t up_index = sector_index / mask_count;
  const size_t down_index = sector_index % mask_count;
  if (masks == NULL || configuration == NULL || up_index >= mask_count ||
      down_index >= mask_count) {
    return 0;
  }
  *configuration =
      (uint64_t)masks[up_index] |
      ((uint64_t)masks[down_index] << (unsigned int)site_count);
  return 1;
}

static uint64_t markov_policy_hash(int site_count, int qp_total,
                                   size_t cache_bytes, double rho,
                                   double eta,
                                   const double *log_basis_scale) {
  uint64_t hash = UINT64_C(1469598103934665603);
  int order;
  markov_hash_u64(&hash, UINT64_C(0x5034534d41524b4f));
  markov_hash_u64(&hash, (uint64_t)(unsigned int)site_count);
  markov_hash_u64(&hash, (uint64_t)(unsigned int)qp_total);
  markov_hash_u64(&hash, (uint64_t)cache_bytes);
  markov_hash_u64(&hash, markov_double_bits(rho));
  markov_hash_u64(&hash, markov_double_bits(eta));
  for (order = 0; order <= MARKOV_ORDER; ++order) {
    markov_hash_u64(&hash, markov_double_bits(log_basis_scale[order]));
  }
  return hash;
}

static int markov_init_guide_policy(
    int site_count, int qp_total, size_t cache_bytes, double rho, double eta,
    const double *log_basis_scale,
    MVMCKrylovPositiveGuidePolicy *policy) {
  int order;
  if (policy == NULL || log_basis_scale == NULL || !isfinite(eta) ||
      eta <= 0.0) {
    return 0;
  }
  memset(policy, 0, sizeof(*policy));
  policy->order = MARKOV_ORDER;
  policy->eta = eta;
  for (order = 0; order <= MARKOV_ORDER; ++order) {
    if (!isfinite(log_basis_scale[order])) return 0;
    policy->lambda[order] = 1.0;
    policy->log_basis_scale[order] = log_basis_scale[order];
  }
  policy->policy_hash =
      markov_policy_hash(site_count, qp_total, cache_bytes, rho, eta,
                         log_basis_scale);
  return 1;
}

static int markov_init_measurement_policy(
    double eta, const double *log_basis_scale,
    MVMCKrylovMatrixMeasurementPolicy *policy) {
  int order;
  if (policy == NULL || log_basis_scale == NULL || !isfinite(eta) ||
      eta <= 0.0) {
    return 0;
  }
  memset(policy, 0, sizeof(*policy));
  policy->order = MARKOV_ORDER;
  policy->eta = eta;
  for (order = 0; order <= MARKOV_ORDER; ++order) {
    if (!isfinite(log_basis_scale[order])) return 0;
    policy->guide_lambda[order] = 1.0;
    policy->target_weight[order] = 1.0;
    policy->log_basis_scale[order] = log_basis_scale[order];
  }
  return 1;
}

static int markov_evaluate_exact_sector(
    MVMCKrylovBoundedWorkspace *krylov_workspace,
    ProfileScaledAmplitude *amplitude_workspace, const unsigned int *masks,
    size_t mask_count, size_t sector_size, int site_count,
    MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT],
    double complex (*raw_values)[PROFILE_DEPTH_COUNT],
    uint64_t *configurations, double norms[PROFILE_DEPTH_COUNT]) {
  size_t sector_index;
  int order;
  memset(norms, 0, PROFILE_DEPTH_COUNT * sizeof(norms[0]));
  for (sector_index = 0; sector_index < sector_size; ++sector_index) {
    uint64_t configuration = 0;
    MVMCKrylovBoundedResult result;
    MVMCKrylovStatus status;
    if (!markov_configuration_at(sector_index, masks, mask_count,
                                 site_count, &configuration)) {
      return 0;
    }
    status = mvmc_bounded_krylov_evaluate(
        krylov_workspace, &configuration, 1, profile_scaled_amplitude,
        amplitude_workspace, &result);
    if (!collective_all_ready(status == MVMC_KRYLOV_STATUS_OK &&
                              result.valid &&
                              result.evaluated_order == MARKOV_ORDER)) {
      if (world_rank() == 0) {
        fprintf(stderr,
                "MARKOV exact sector state=%zu failed collectively: %s\n",
                sector_index, mvmc_krylov_status_string(status));
      }
      return 0;
    }
    configurations[sector_index] = configuration;
    for (order = 0; order <= MARKOV_ORDER; ++order) {
      double complex raw = 0.0;
      exact_values[sector_index][order] = result.value[order];
      if (!scaled_to_raw(&result.value[order], &raw) ||
          !isfinite(creal(raw)) || !isfinite(cimag(raw))) {
        return 0;
      }
      raw_values[sector_index][order] = raw;
      norms[order] += creal(conj(raw) * raw);
      if (!isfinite(norms[order])) return 0;
    }
  }
  for (order = 0; order <= MARKOV_ORDER; ++order) {
    if (!isfinite(norms[order]) || norms[order] <= 0.0) return 0;
  }
  return 1;
}

static int markov_compute_scales_and_eta(
    size_t sector_size,
    const double complex (*raw_values)[PROFILE_DEPTH_COUNT],
    const double norms[PROFILE_DEPTH_COUNT], double rho,
    double log_basis_scale[PROFILE_DEPTH_COUNT], double *eta) {
  size_t sector_index;
  int order;
  double target_sum = 0.0;
  if (raw_values == NULL || norms == NULL || log_basis_scale == NULL ||
      eta == NULL || !isfinite(rho) || rho <= 0.0 ||
      !isfinite(norms[0]) || norms[0] <= 0.0) {
    return 0;
  }
  for (order = 0; order <= MARKOV_ORDER; ++order) {
    if (!isfinite(norms[order]) || norms[order] <= 0.0) return 0;
    log_basis_scale[order] =
        order == 0 ? 0.0 : 0.5 * (log(norms[0]) - log(norms[order]));
    if (!isfinite(log_basis_scale[order])) return 0;
  }
  for (sector_index = 0; sector_index < sector_size; ++sector_index) {
    double target = 0.0;
    for (order = 0; order <= MARKOV_ORDER; ++order) {
      const double scaled_abs =
          cabs(raw_values[sector_index][order]) *
          exp(log_basis_scale[order]);
      target += scaled_abs * scaled_abs;
    }
    if (!isfinite(target) || target <= 0.0) return 0;
    target_sum += target;
    if (!isfinite(target_sum)) return 0;
  }
  *eta = rho * (target_sum / (double)sector_size);
  return isfinite(*eta) && *eta > 0.0;
}

static int markov_compute_exact_budget(
    const MVMCKrylovMatrixMeasurementPolicy *measurement_policy,
    const MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT],
    size_t sector_size, size_t upper_count, MarkovExactBudget *budget,
    double *guide_weights) {
  size_t sector_index;
  MarkovMatrixFamily family;
  if (measurement_policy == NULL || exact_values == NULL ||
      budget == NULL || guide_weights == NULL) {
    return 0;
  }
  memset(budget, 0, sizeof(*budget));
  for (sector_index = 0; sector_index < sector_size; ++sector_index) {
    double complex overlap[6];
    double complex hamiltonian[6];
    double complex hamiltonian_adjoint[6];
    double complex hamiltonian_squared[6];
    MVMCKrylovMatrixMeasurementSampleDiagnostics diagnostics;
    double guide;
    MVMCKrylovStatus status =
        mvmc_krylov_matrix_measurement_sample_with_adjoint(
            measurement_policy, exact_values[sector_index],
            PROFILE_DEPTH_COUNT, overlap, hamiltonian,
            hamiltonian_adjoint, hamiltonian_squared, upper_count,
            &diagnostics);
    if (status != MVMC_KRYLOV_STATUS_OK || !diagnostics.valid ||
        !isfinite(diagnostics.log_guide)) {
      return 0;
    }
    guide = exp(diagnostics.log_guide);
    if (!isfinite(guide) || guide <= 0.0) return 0;
    guide_weights[sector_index] = guide;
    budget->target_sum += diagnostics.denominator * guide;
    if (!isfinite(budget->target_sum)) return 0;
    for (family = MARKOV_FAMILY_S; family <= MARKOV_FAMILY_B; ++family) {
      double complex *values = markov_family_values(
          family, overlap, hamiltonian, hamiltonian_squared);
      size_t entry;
      if (values == NULL) return 0;
      for (entry = 0; entry < upper_count; ++entry) {
        const double complex raw = values[entry] * guide;
        MarkovEntryBudget *item = &budget->entry[family][entry];
        if (!isfinite(creal(raw)) || !isfinite(cimag(raw))) return 0;
        item->exact_sum += raw;
        item->unsigned_sum += cabs(raw);
        if (!isfinite(creal(item->exact_sum)) ||
            !isfinite(cimag(item->exact_sum)) ||
            !isfinite(item->unsigned_sum)) {
          return 0;
        }
      }
    }
  }
  if (!isfinite(budget->target_sum) || budget->target_sum <= 0.0) {
    return 0;
  }
  for (family = MARKOV_FAMILY_S; family <= MARKOV_FAMILY_B; ++family) {
    size_t entry;
    for (entry = 0; entry < upper_count; ++entry) {
      MarkovEntryBudget *item = &budget->entry[family][entry];
      const double unsigned_scale = item->unsigned_sum / budget->target_sum;
      item->exact_theta = item->exact_sum / budget->target_sum;
      item->budget = fmax(1.0e-12, 0.02 * unsigned_scale);
      if (!isfinite(creal(item->exact_theta)) ||
          !isfinite(cimag(item->exact_theta)) ||
          !isfinite(item->budget) || item->budget <= 0.0) {
        return 0;
      }
    }
  }
  return 1;
}

static size_t markov_inverse_cdf(double uniform, const double *weights,
                                 size_t count, double *total_weight) {
  double total = 0.0;
  double threshold;
  double cumulative = 0.0;
  size_t index;
  for (index = 0; index < count; ++index) {
    total += weights[index];
  }
  threshold = uniform * total;
  for (index = 0; index < count; ++index) {
    cumulative += weights[index];
    if (threshold <= cumulative || index + 1 == count) {
      if (total_weight != NULL) *total_weight = total;
      return index;
    }
  }
  if (total_weight != NULL) *total_weight = total;
  return count - 1;
}

static int markov_accumulator_storage_size(size_t block_count,
                                           size_t upper_count,
                                           size_t *value) {
  return checked_multiply_size(block_count, upper_count, value);
}

static int markov_init_block_accumulator(
    size_t block_count, uint64_t block_length, size_t upper_count,
    MVMCKrylovMatrixMeasurementBlockAccumulator *accumulator,
    MVMCKrylovMatrixMeasurementAccumulator **blocks,
    MVMCKrylovStreamingComplexSum **overlap,
    MVMCKrylovStreamingComplexSum **hamiltonian,
    MVMCKrylovStreamingComplexSum **hamiltonian_adjoint,
    MVMCKrylovStreamingComplexSum **hamiltonian_squared) {
  size_t storage_count = 0;
  if (accumulator == NULL || blocks == NULL || overlap == NULL ||
      hamiltonian == NULL || hamiltonian_adjoint == NULL ||
      hamiltonian_squared == NULL ||
      !markov_accumulator_storage_size(block_count, upper_count,
                                       &storage_count)) {
    return 0;
  }
  *blocks = (MVMCKrylovMatrixMeasurementAccumulator *)calloc(
      block_count, sizeof(**blocks));
  *overlap = (MVMCKrylovStreamingComplexSum *)calloc(
      storage_count, sizeof(**overlap));
  *hamiltonian = (MVMCKrylovStreamingComplexSum *)calloc(
      storage_count, sizeof(**hamiltonian));
  *hamiltonian_adjoint = (MVMCKrylovStreamingComplexSum *)calloc(
      storage_count, sizeof(**hamiltonian_adjoint));
  *hamiltonian_squared = (MVMCKrylovStreamingComplexSum *)calloc(
      storage_count, sizeof(**hamiltonian_squared));
  if (*blocks == NULL || *overlap == NULL || *hamiltonian == NULL ||
      *hamiltonian_adjoint == NULL || *hamiltonian_squared == NULL) {
    return 0;
  }
  return mvmc_krylov_matrix_measurement_block_accumulator_init(
             MARKOV_ORDER, block_count, block_length, *blocks, *overlap,
             *hamiltonian, *hamiltonian_adjoint, *hamiltonian_squared,
             storage_count, accumulator) == MVMC_KRYLOV_STATUS_OK;
}

static void markov_free_block_storage(
    MVMCKrylovMatrixMeasurementAccumulator *blocks,
    MVMCKrylovStreamingComplexSum *overlap,
    MVMCKrylovStreamingComplexSum *hamiltonian,
    MVMCKrylovStreamingComplexSum *hamiltonian_adjoint,
    MVMCKrylovStreamingComplexSum *hamiltonian_squared) {
  free(hamiltonian_squared);
  free(hamiltonian_adjoint);
  free(hamiltonian);
  free(overlap);
  free(blocks);
}

static int markov_print_tau(const char *series_name, const char *component,
                            const double *values, size_t sample_count,
                            double block_length, double maximum_tau_int,
                            double block_multiplier, int required,
                            int require_nonzero_variance,
                            double *maximum_seen, int *all_pass) {
  MVMCKrylovTauIntResult tau;
  MVMCKrylovTauIntGateResult gate;
  int effective_required = required;
  MVMCKrylovStatus status =
      mvmc_krylov_tau_int_geyer_initial_positive(values, sample_count, &tau);
  if (status != MVMC_KRYLOV_STATUS_OK || !tau.valid) return 0;
  if (require_nonzero_variance && tau.variance == 0.0) {
    effective_required = 0;
  }
  status = mvmc_krylov_tau_int_gate_check(
      tau.tau_int, block_length, maximum_tau_int, block_multiplier, &gate);
  if (status != MVMC_KRYLOV_STATUS_OK || !gate.valid) return 0;
  if (maximum_seen != NULL && tau.tau_int > *maximum_seen) {
    *maximum_seen = tau.tau_int;
  }
  if (all_pass != NULL && effective_required && !gate.passed) *all_pass = 0;
  printf("TAU series=%s component=%s sample_count=%zu variance=%.17g",
         series_name, component, sample_count, tau.variance);
  printf(" positive_pair_count=%zu tau_int=%.17g", tau.positive_pair_count,
         tau.tau_int);
  printf(" block_length=%.17g block_multiplier=%.17g pass=%d required=%d\n",
         block_length, gate.actual_block_multiplier, gate.passed,
         effective_required);
  return 1;
}

static int markov_entry_indices(size_t upper_count, size_t entry,
                                size_t *row, size_t *column) {
  size_t r;
  size_t cursor = 0;
  for (r = 0; r < MARKOV_DIMENSION; ++r) {
    size_t c;
    for (c = r; c < MARKOV_DIMENSION; ++c) {
      if (cursor == entry) {
        *row = r;
        *column = c;
        return 1;
      }
      ++cursor;
    }
  }
  (void)upper_count;
  return 0;
}

static int run_markov(int site_count, int qp_total, int sample_count,
                      size_t cache_bytes, double rho, uint64_t seed,
                      const MVMCKrylovPositiveSamplerProposalPolicy
                          *requested_proposal_policy,
                      uint64_t rng_stream, int persistent_session) {
  ProfileModel fixture;
  MVMCClassicKrylovModelWorkspace *model_workspace = NULL;
  ProfileScaledAmplitude *amplitude_workspace = NULL;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *krylov_workspace = NULL;
  MVMCKrylovBoundedCollectiveWorkspace *collective_workspace = NULL;
  const MVMCKrylovFockModel *model = NULL;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovPositiveGuidePolicy guide_policy;
  MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
  MVMCKrylovMatrixMeasurementPolicy measurement_policy;
  MVMCKrylovPositiveSamplerSnapshot current;
  MVMCKrylovPositiveSamplerRng rng;
  MVMCKrylovPositiveSamplerManifest final_manifest;
  MVMCKrylovPositiveSamplerTraceStatistics trace_statistics;
  MVMCKrylovBoundedStatistics final_session_statistics;
  MVMCKrylovMatrixMeasurementAccumulator total_accumulator;
  MVMCKrylovMatrixMeasurementBlockAccumulator official_accumulator;
  MVMCKrylovMatrixMeasurementBlockAccumulator diagnostic_accumulator;
  MVMCKrylovMatrixMeasurementDiagnosticsSummary diagnostics_summary;
  MVMCKrylovBoundedCollectiveResult collective_result;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  MVMCKrylovStreamingComplexSum *total_overlap = NULL;
  MVMCKrylovStreamingComplexSum *total_hamiltonian = NULL;
  MVMCKrylovStreamingComplexSum *total_hamiltonian_adjoint = NULL;
  MVMCKrylovStreamingComplexSum *total_hamiltonian_squared = NULL;
  MVMCKrylovMatrixMeasurementAccumulator *official_blocks = NULL;
  MVMCKrylovStreamingComplexSum *official_overlap = NULL;
  MVMCKrylovStreamingComplexSum *official_hamiltonian = NULL;
  MVMCKrylovStreamingComplexSum *official_hamiltonian_adjoint = NULL;
  MVMCKrylovStreamingComplexSum *official_hamiltonian_squared = NULL;
  MVMCKrylovMatrixMeasurementAccumulator *diagnostic_blocks = NULL;
  MVMCKrylovStreamingComplexSum *diagnostic_overlap = NULL;
  MVMCKrylovStreamingComplexSum *diagnostic_hamiltonian = NULL;
  MVMCKrylovStreamingComplexSum *diagnostic_hamiltonian_adjoint = NULL;
  MVMCKrylovStreamingComplexSum *diagnostic_hamiltonian_squared = NULL;
  unsigned int *masks = NULL;
  MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT] = NULL;
  double complex (*raw_values)[PROFILE_DEPTH_COUNT] = NULL;
  uint64_t *configurations = NULL;
  double *guide_weights = NULL;
  double *w_series = NULL;
  double complex *entry_series = NULL;
  double *scratch_series = NULL;
  double *slater = NULL;
  double complex *weights = NULL;
  double norms[PROFILE_DEPTH_COUNT];
  double log_basis_scale[PROFILE_DEPTH_COUNT];
  double eta = NAN;
  double initial_uniform = NAN;
  double total_guide_weight = NAN;
  MarkovExactBudget exact_budget;
  const double projection_parameter = -0.27;
  const double minimum_leave_one_denominator = 1.0e-12;
  const double minimum_abs_denominator_mean = 1.0e-12;
  const double maximum_denominator_relative_se = 1.0;
  const double maximum_tau_int = 16.0;
  const double block_length_multiplier = 16.0;
  const double maximum_block_stability_ratio = 1.25;
  const double maximum_conservative_se_budget_ratio = 1.0;
  const size_t official_block_count = MARKOV_OFFICIAL_BLOCKS;
  const size_t diagnostic_block_count = MARKOV_DIAGNOSTIC_BLOCKS;
  size_t official_block_length = 0;
  size_t diagnostic_block_length = 0;
  size_t mask_count = 0;
  size_t sector_size = 0;
  size_t upper_count = 0;
  size_t dimension = 0;
  size_t series_count = 0;
  size_t initial_sector_index = 0;
  uint64_t current_words[1] = {0};
  uint64_t proposal_words[1] = {0};
  uint64_t plan_hash = 0;
  uint64_t amplitude_generation_hash = 0;
  int local_ready = 1;
  int return_code = 1;
  int session_started = 0;
  int p4s_pass = 1;
  int tau_pass = 1;
  int block_stability_pass = 1;
  int block_pathology_pass = 1;
  int budget_pass = 1;
  int denominator_pass = 1;
  double maximum_tau_seen = 0.0;
  double maximum_official_se_budget_ratio = 0.0;
  double maximum_se_budget_ratio = 0.0;
  double maximum_stability_ratio = 0.0;
  double local_timing[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double rank_max_timing[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int rank = world_rank();
  size_t sample;
  int order;
  MarkovMatrixFamily family;

  memset(&fixture, 0, sizeof(fixture));
  memset(&current, 0, sizeof(current));
  memset(&rng, 0, sizeof(rng));
  memset(&proposal_policy, 0, sizeof(proposal_policy));
  memset(&final_manifest, 0, sizeof(final_manifest));
  memset(&trace_statistics, 0, sizeof(trace_statistics));
  memset(&final_session_statistics, 0, sizeof(final_session_statistics));
  if (requested_proposal_policy == NULL ||
      !requested_proposal_policy->valid ||
      requested_proposal_policy->status != MVMC_KRYLOV_STATUS_OK ||
      sample_count <= 0 ||
      ((size_t)sample_count % official_block_count) != 0 ||
      ((size_t)sample_count % diagnostic_block_count) != 0) {
    local_ready = 0;
  }
  if (local_ready) {
    official_block_length = (size_t)sample_count / official_block_count;
    diagnostic_block_length = (size_t)sample_count / diagnostic_block_count;
    mask_count = fixed_masks(site_count, &masks);
    local_ready = mask_count != 0 &&
                  checked_multiply_size(mask_count, mask_count,
                                        &sector_size) &&
                  initialize_model(&fixture, site_count) &&
                  initialize_slater(site_count, qp_total, &slater, &weights);
  }
  if (local_ready) {
    exact_values = (MVMCScaledComplex (*)[PROFILE_DEPTH_COUNT])calloc(
        sector_size, sizeof(*exact_values));
    raw_values = (double complex (*)[PROFILE_DEPTH_COUNT])calloc(
        sector_size, sizeof(*raw_values));
    configurations = (uint64_t *)calloc(sector_size, sizeof(*configurations));
    guide_weights = (double *)calloc(sector_size, sizeof(*guide_weights));
    w_series = (double *)calloc((size_t)sample_count, sizeof(*w_series));
    scratch_series =
        (double *)calloc((size_t)sample_count, sizeof(*scratch_series));
    local_ready = exact_values != NULL && raw_values != NULL &&
                  configurations != NULL && guide_weights != NULL &&
                  w_series != NULL && scratch_series != NULL;
  }
  if (!collective_all_ready(local_ready)) {
    if (rank == 0) {
      fprintf(stderr, "bounded Markov fixture allocation/preflight failed\n");
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
      fprintf(stderr, "bounded Markov preflight failed collectively: %s\n",
              mvmc_krylov_status_string(status));
    }
    goto cleanup;
  }
  if (!markov_evaluate_exact_sector(
          krylov_workspace, amplitude_workspace, masks, mask_count,
          sector_size, site_count, exact_values, raw_values,
          configurations, norms) ||
      !markov_compute_scales_and_eta(
          sector_size,
          (const double complex(*)[PROFILE_DEPTH_COUNT])raw_values, norms,
          rho, log_basis_scale, &eta) ||
      !markov_init_guide_policy(site_count, qp_total, cache_bytes, rho, eta,
                                log_basis_scale, &guide_policy) ||
      !markov_init_measurement_policy(eta, log_basis_scale,
                                      &measurement_policy) ||
      mvmc_krylov_matrix_measurement_dimension(
          MARKOV_ORDER, &dimension, &upper_count) != MVMC_KRYLOV_STATUS_OK ||
      upper_count != 6 || dimension != MARKOV_DIMENSION ||
      !markov_compute_exact_budget(
          &measurement_policy,
          (const MVMCScaledComplex(*)[PROFILE_DEPTH_COUNT])exact_values,
          sector_size, upper_count, &exact_budget, guide_weights)) {
    if (rank == 0) {
      fprintf(stderr, "bounded Markov exact calibration failed\n");
    }
    goto cleanup;
  }
  proposal_policy = *requested_proposal_policy;
  if (!checked_multiply_size((size_t)sample_count,
                             MARKOV_FAMILY_COUNT * upper_count,
                             &series_count)) {
    goto cleanup;
  }
  entry_series =
      (double complex *)calloc(series_count, sizeof(*entry_series));
  if (entry_series == NULL) goto cleanup;

  mvmc_bounded_krylov_workspace_destroy(krylov_workspace);
  krylov_workspace = NULL;
  status = mvmc_bounded_krylov_workspace_create(plan, &krylov_workspace);
  if (!collective_all_ready(status == MVMC_KRYLOV_STATUS_OK &&
                            krylov_workspace != NULL)) {
    if (rank == 0) {
      fprintf(stderr, "bounded Markov measurement workspace reset failed\n");
    }
    goto cleanup;
  }

  if (persistent_session) {
    amplitude_generation_hash =
        plan_hash ^ limits.amplitude_policy_hash ^
        UINT64_C(0x503453394f464643);
    if (amplitude_generation_hash == 0) {
      amplitude_generation_hash = UINT64_C(1);
    }
    status = mvmc_bounded_krylov_session_begin(
        krylov_workspace, profile_scaled_amplitude, amplitude_workspace,
        amplitude_generation_hash);
    session_started = status == MVMC_KRYLOV_STATUS_OK;
    if (!collective_all_ready(session_started)) {
      if (rank == 0) {
        fprintf(stderr, "bounded Markov session begin failed: %s\n",
                mvmc_krylov_status_string(status));
      }
      goto cleanup;
    }
  }

  status = mvmc_krylov_positive_sampler_rng_seed(seed, rng_stream, &rng);
  if (status != MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_rng_draw_uniform(
          &rng, &initial_uniform) != MVMC_KRYLOV_STATUS_OK) {
    goto cleanup;
  }
  initial_sector_index = markov_inverse_cdf(
      initial_uniform, guide_weights, sector_size, &total_guide_weight);
  current_words[0] = configurations[initial_sector_index];
  status = mvmc_krylov_positive_sampler_initialize(
      krylov_workspace, &guide_policy, current_words, 1,
      profile_scaled_amplitude, amplitude_workspace, &current);
  if (!collective_all_ready(status == MVMC_KRYLOV_STATUS_OK &&
                            current.valid)) {
    if (rank == 0) {
      fprintf(stderr, "bounded Markov sampler initialize failed: %s\n",
              mvmc_krylov_status_string(status));
    }
    goto cleanup;
  }
  if (persistent_session &&
      !collective_all_ready(
          current.krylov.statistics.persistent_session_active &&
          current.krylov.statistics.amplitude_generation_hash ==
              amplitude_generation_hash &&
          current.krylov.statistics.session_root_evaluation == UINT64_C(1) &&
          current.krylov.statistics.cache_reset_performed == UINT64_C(1) &&
          current.krylov.statistics.engine_heap_allocations == UINT64_C(0))) {
    if (rank == 0) {
      fprintf(stderr, "bounded Markov session initialization audit failed\n");
    }
    goto cleanup;
  }

  if (rank == 0) {
    total_overlap = (MVMCKrylovStreamingComplexSum *)calloc(
        upper_count, sizeof(*total_overlap));
    total_hamiltonian = (MVMCKrylovStreamingComplexSum *)calloc(
        upper_count, sizeof(*total_hamiltonian));
    total_hamiltonian_adjoint = (MVMCKrylovStreamingComplexSum *)calloc(
        upper_count, sizeof(*total_hamiltonian_adjoint));
    total_hamiltonian_squared = (MVMCKrylovStreamingComplexSum *)calloc(
        upper_count, sizeof(*total_hamiltonian_squared));
    if (total_overlap == NULL || total_hamiltonian == NULL ||
        total_hamiltonian_adjoint == NULL ||
        total_hamiltonian_squared == NULL ||
        mvmc_krylov_matrix_measurement_accumulator_init_with_diagnostics(
            MARKOV_ORDER, total_overlap, total_hamiltonian,
            total_hamiltonian_adjoint, total_hamiltonian_squared,
            upper_count, &total_accumulator) != MVMC_KRYLOV_STATUS_OK ||
        !markov_init_block_accumulator(
            official_block_count, official_block_length, upper_count,
            &official_accumulator, &official_blocks, &official_overlap,
            &official_hamiltonian, &official_hamiltonian_adjoint,
            &official_hamiltonian_squared) ||
        !markov_init_block_accumulator(
            diagnostic_block_count, diagnostic_block_length, upper_count,
            &diagnostic_accumulator, &diagnostic_blocks,
            &diagnostic_overlap, &diagnostic_hamiltonian,
            &diagnostic_hamiltonian_adjoint,
            &diagnostic_hamiltonian_squared) ||
        mvmc_krylov_positive_sampler_trace_statistics_reset(
            &trace_statistics) != MVMC_KRYLOV_STATUS_OK) {
      fprintf(stderr, "bounded Markov accumulator initialization failed\n");
      goto cleanup;
    }

    printf("MARKOV schema=%d fixture=%s",
           persistent_session ? MARKOV_SESSION_SCHEMA_VERSION
                              : MARKOV_SCHEMA_VERSION,
           persistent_session ? "p4s9_long_direct_session_official"
                              : "p4s_bounded_markov_real");
    printf(" site_count=%d qp_total=%d sample_count=%d sector_size=%zu",
           site_count, qp_total, sample_count, sector_size);
    printf(" rank_count=%d cache_requested_bytes=%zu rho=%.17g eta=%.17g",
           world_size(), cache_bytes, rho, eta);
    printf(" seed=%" PRIu64 " seed_hex=0x%016" PRIx64,
           seed, seed);
    printf(" rng_stream=%" PRIu64, rng_stream);
    printf(" persistent_session=%d amplitude_generation_hash=%" PRIu64,
           persistent_session, amplitude_generation_hash);
    printf(" proposal_kernel=%" PRIu64,
           proposal_policy.kernel_id);
    printf(" global_numerator=%zu global_denominator=%zu",
           proposal_policy.global_numerator,
           proposal_policy.global_denominator);
    printf(" neighbor_numerator=%zu neighbor_denominator=%zu",
           proposal_policy.neighbor_numerator,
           proposal_policy.neighbor_denominator);
    printf(" distance_numerator=%zu distance_denominator=%zu",
           proposal_policy.distance_numerator,
           proposal_policy.distance_denominator);
    printf(" distance_rounding_rule=%" PRIu64,
           proposal_policy.distance_rounding_rule);
    printf(" proposal_policy_hash=%" PRIu64,
           proposal_policy.policy_hash);
    printf(" initial_uniform=%.17g initial_sector_index=%zu",
           initial_uniform, initial_sector_index);
    printf(" initial_configuration=%" PRIu64, current_words[0]);
    printf(" total_guide_weight=%.17g", total_guide_weight);
    printf(" official_block_count=%zu official_block_length=%zu",
           official_block_count, official_block_length);
    printf(" diagnostic_block_count=%zu diagnostic_block_length=%zu",
           diagnostic_block_count, diagnostic_block_length);
    printf(" minimum_leave_one_denominator=%.17g",
           minimum_leave_one_denominator);
    printf(" minimum_abs_denominator_mean=%.17g",
           minimum_abs_denominator_mean);
    printf(" maximum_denominator_relative_se=%.17g",
           maximum_denominator_relative_se);
    printf(" maximum_tau_int=%.17g block_length_multiplier=%.17g",
           maximum_tau_int, block_length_multiplier);
    printf(" maximum_block_stability_ratio=%.17g",
           maximum_block_stability_ratio);
    printf(" maximum_conservative_se_budget_ratio=%.17g\n",
           maximum_conservative_se_budget_ratio);
    for (order = 0; order <= MARKOV_ORDER; ++order) {
      printf("SCALE order=%d norm=%.17g log_basis_scale=%.17g\n",
             order, norms[order], log_basis_scale[order]);
    }
  }

  for (sample = 0; sample < (size_t)sample_count; ++sample) {
    MVMCKrylovPositiveSamplerProposalStepResult step;
    MVMCKrylovMatrixMeasurementSampleDiagnostics diagnostics;
    double complex overlap[6];
    double complex hamiltonian[6];
    double complex hamiltonian_adjoint[6];
    double complex hamiltonian_squared[6];
    MarkovMatrixFamily family;
    status = mvmc_krylov_positive_sampler_step_mixture_rng(
        krylov_workspace, &guide_policy, model, &proposal_policy,
        current_words, 1, &current, &rng, profile_scaled_amplitude,
        amplitude_workspace, proposal_words, 1, &step);
    if (!collective_all_ready(status == MVMC_KRYLOV_STATUS_OK &&
                              step.valid && current.valid &&
                              isfinite(step.proposal_seconds) &&
                              step.proposal_seconds >= 0.0 &&
                              isfinite(step.component_selection_seconds) &&
                              step.component_selection_seconds >= 0.0 &&
                              isfinite(step.global_subset_seconds) &&
                              step.global_subset_seconds >= 0.0 &&
                              isfinite(step.shell_generation_seconds) &&
                              step.shell_generation_seconds >= 0.0 &&
                              isfinite(step.bounded_evaluation_seconds) &&
                              step.bounded_evaluation_seconds >= 0.0 &&
                              isfinite(step.total_step_seconds) &&
                              step.total_step_seconds >= 0.0)) {
      if (rank == 0) {
        fprintf(stderr, "bounded Markov step %zu failed: %s\n", sample,
                mvmc_krylov_status_string(status));
      }
      goto cleanup;
    }
    if (persistent_session &&
        !collective_all_ready(
            step.step.proposal_krylov.statistics.persistent_session_active &&
            step.step.proposal_krylov.statistics.amplitude_generation_hash ==
                amplitude_generation_hash &&
            step.step.proposal_krylov.statistics.session_root_evaluation ==
                (uint64_t)sample + UINT64_C(2) &&
            step.step.proposal_krylov.statistics.cache_reset_performed ==
                UINT64_C(0) &&
            step.step.proposal_krylov.statistics.engine_heap_allocations ==
                UINT64_C(0))) {
      if (rank == 0) {
        fprintf(stderr, "bounded Markov session step audit failed at %zu\n",
                sample);
      }
      goto cleanup;
    }
    if (persistent_session) {
      final_session_statistics = step.step.proposal_krylov.statistics;
    }
    local_timing[0] += step.proposal_seconds;
    local_timing[1] += step.component_selection_seconds;
    local_timing[2] += step.global_subset_seconds;
    local_timing[3] += step.shell_generation_seconds;
    local_timing[4] += step.bounded_evaluation_seconds;
    local_timing[5] += step.total_step_seconds;
    if (rank != 0) continue;

    status = mvmc_krylov_matrix_measurement_sample_with_adjoint(
        &measurement_policy, current.krylov.value, PROFILE_DEPTH_COUNT,
        overlap, hamiltonian, hamiltonian_adjoint, hamiltonian_squared,
        upper_count, &diagnostics);
    if (status != MVMC_KRYLOV_STATUS_OK || !diagnostics.valid ||
        mvmc_krylov_matrix_measurement_accumulator_add_sample(
            &total_accumulator, &measurement_policy, current.krylov.value,
            PROFILE_DEPTH_COUNT, NULL) != MVMC_KRYLOV_STATUS_OK ||
        mvmc_krylov_matrix_measurement_block_accumulator_add_sample(
            &official_accumulator, &measurement_policy, current.krylov.value,
            PROFILE_DEPTH_COUNT, NULL) != MVMC_KRYLOV_STATUS_OK ||
        mvmc_krylov_matrix_measurement_block_accumulator_add_sample(
            &diagnostic_accumulator, &measurement_policy,
            current.krylov.value, PROFILE_DEPTH_COUNT,
            NULL) != MVMC_KRYLOV_STATUS_OK ||
        mvmc_krylov_positive_sampler_trace_statistics_record_step(
            &step, current_words, 1, &trace_statistics) !=
            MVMC_KRYLOV_STATUS_OK) {
      fprintf(stderr, "bounded Markov measurement add failed at sample %zu\n",
              sample);
      goto cleanup;
    }

    w_series[sample] = diagnostics.denominator;
    for (family = MARKOV_FAMILY_S; family <= MARKOV_FAMILY_B; ++family) {
      double complex *values = markov_family_values(
          family, overlap, hamiltonian, hamiltonian_squared);
      size_t entry;
      if (values == NULL) goto cleanup;
      for (entry = 0; entry < upper_count; ++entry) {
        entry_series[markov_series_index((size_t)sample_count, family,
                                         entry, upper_count, sample)] =
            values[entry];
      }
    }
    printf("SAMPLE sample=%zu configuration=%" PRIu64, sample,
           current_words[0]);
    printf(" accepted=%d component=%d self_proposal=%d",
           step.step.accepted, (int)step.component, step.self_proposal);
    printf(" configuration_changed=%d selected_neighbor=%zu",
           step.configuration_changed, step.selected_neighbor_index);
    printf(" component_draw_count=%zu global_subset_draw_count=%zu",
           step.component_draw_count, step.global_subset_draw_count);
    printf(" shell_draw_count=%zu shell_max_distance=%zu",
           step.shell_draw_count, step.shell_max_distance);
    printf(" shell_distance=%zu shell_up_distance=%zu",
           step.shell_distance, step.shell_up_distance);
    printf(" shell_down_distance=%zu shell_count=%zu",
           step.shell_down_distance, step.shell_count);
    printf(" uniform=%.17g", step.uniform);
    printf(" log_acceptance_ratio=%.17g denominator=%.17g",
           step.step.acceptance.log_acceptance_ratio,
           diagnostics.denominator);
    printf(" zero_target=%d accepted_generation=%" PRIu64 "\n",
           diagnostics.zero_target_sample, current.accepted_generation);
  }

  if (persistent_session) {
    status = mvmc_bounded_krylov_session_end(krylov_workspace);
    if (!collective_all_ready(status == MVMC_KRYLOV_STATUS_OK)) {
      if (rank == 0) {
        fprintf(stderr, "bounded Markov session end failed: %s\n",
                mvmc_krylov_status_string(status));
      }
      goto cleanup;
    }
    session_started = 0;
  }

  if (!markov_reduce_max_times(local_timing, rank_max_timing, 6)) {
    goto cleanup;
  }

  if (rank == 0) {
    printf("TIMING sample_count=%d", sample_count);
    printf(" proposal_seconds=%.17g component_selection_seconds=%.17g",
           rank_max_timing[0], rank_max_timing[1]);
    printf(" global_subset_seconds=%.17g bounded_evaluation_seconds=%.17g",
           rank_max_timing[2], rank_max_timing[4]);
    printf(" shell_generation_seconds=%.17g", rank_max_timing[3]);
    printf(" total_step_seconds=%.17g total_seconds_per_step=%.17g\n",
           rank_max_timing[5], rank_max_timing[5] / (double)sample_count);
    if (persistent_session) {
      printf("SESSION schema=%d amplitude_generation_hash=%" PRIu64,
             MARKOV_SESSION_SCHEMA_VERSION, amplitude_generation_hash);
      printf(" session_root_evaluations=%" PRIu64,
             final_session_statistics.session_root_evaluation);
      printf(" cache_reset_count=1");
      printf(" cache_resident_end=%zu cache_entries_peak=%zu",
             final_session_statistics.cache_entries_resident_end,
             final_session_statistics.cache_entries_peak);
      printf(" cache_hits=%" PRIu64 " cache_misses=%" PRIu64,
             final_session_statistics.cache_hits[0] +
                 final_session_statistics.cache_hits[1] +
                 final_session_statistics.cache_hits[2] +
                 final_session_statistics.cache_hits[3],
             final_session_statistics.cache_misses[0] +
                 final_session_statistics.cache_misses[1] +
                 final_session_statistics.cache_misses[2] +
                 final_session_statistics.cache_misses[3]);
      printf(" engine_heap_allocations=%" PRIu64 " session_end_pass=1\n",
             final_session_statistics.engine_heap_allocations);
    }
    printf("TRACE attempted=%" PRIu64 " completed=%" PRIu64,
           trace_statistics.attempted_steps,
           trace_statistics.completed_steps);
    printf(" accepted=%" PRIu64 " rejected=%" PRIu64,
           trace_statistics.accepted_steps, trace_statistics.rejected_steps);
    printf(" neighbor_attempted=%" PRIu64 " neighbor_accepted=%" PRIu64,
           trace_statistics.neighbor_attempted_steps,
           trace_statistics.neighbor_accepted_steps);
    printf(" neighbor_rejected=%" PRIu64,
           trace_statistics.neighbor_rejected_steps);
    printf(" global_attempted=%" PRIu64 " global_accepted=%" PRIu64,
           trace_statistics.global_attempted_steps,
           trace_statistics.global_accepted_steps);
    printf(" global_rejected=%" PRIu64 " global_self=%" PRIu64,
           trace_statistics.global_rejected_steps,
           trace_statistics.global_self_proposals);
    printf(" shell_attempted=%" PRIu64 " shell_accepted=%" PRIu64,
           trace_statistics.shell_attempted_steps,
           trace_statistics.shell_accepted_steps);
    printf(" shell_rejected=%" PRIu64,
           trace_statistics.shell_rejected_steps);
    printf(" configuration_changing_accepted=%" PRIu64,
           trace_statistics.configuration_changing_accepted_moves);
    printf(" positive_support=%" PRIu64 " support_violation=%" PRIu64,
           trace_statistics.positive_support_steps,
           trace_statistics.support_violation_steps);
    printf(" floor_supported_zero=%" PRIu64,
           trace_statistics.floor_supported_zero_steps);
    printf(" finite_components=%" PRIu64 " exact_zero_components=%" PRIu64,
           trace_statistics.finite_proposal_components,
           trace_statistics.exact_zero_proposal_components);
    printf(" numeric_zero_components=%" PRIu64,
           trace_statistics.numeric_zero_proposal_components);
    printf(" terminal_amplitude_calls=%" PRIu64,
           trace_statistics.terminal_amplitude_calls);
    printf(" min_log_acceptance_ratio=%.17g",
           trace_statistics.min_log_acceptance_ratio);
    printf(" max_log_acceptance_ratio=%.17g",
           trace_statistics.max_log_acceptance_ratio);
    printf(" sum_log_acceptance_ratio=%.17g",
           trace_statistics.sum_log_acceptance_ratio);
    printf(" min_proposal_log_guide=%.17g",
           trace_statistics.min_proposal_log_guide);
    printf(" max_proposal_log_guide=%.17g",
           trace_statistics.max_proposal_log_guide);
    printf(" trace_hash=%" PRIu64 "\n", trace_statistics.trace_hash);

    status = mvmc_krylov_matrix_measurement_diagnostics_summary(
        &total_accumulator, minimum_abs_denominator_mean,
        maximum_denominator_relative_se, &diagnostics_summary);
    if (status != MVMC_KRYLOV_STATUS_OK || !diagnostics_summary.valid) {
      fprintf(stderr, "bounded Markov diagnostics summary failed\n");
      goto cleanup;
    }
    denominator_pass = diagnostics_summary.denominator_stable;
    printf("SUMMARY sample_count=%" PRIu64,
           diagnostics_summary.sample_count);
    printf(" zero_target_sample_count=%" PRIu64,
           diagnostics_summary.zero_target_sample_count);
    printf(" denominator_sum=%.17g denominator_mean=%.17g",
           diagnostics_summary.denominator_sum,
           diagnostics_summary.denominator_mean);
    printf(" denominator_relative_se=%.17g denominator_stable=%d",
           diagnostics_summary.denominator_relative_se,
           diagnostics_summary.denominator_stable);
    printf(" effective_sample_fraction=%.17g",
           diagnostics_summary.effective_sample_fraction);
    printf(" zero_target_sample_fraction=%.17g",
           diagnostics_summary.zero_target_sample_fraction);
    printf(" minimum_denominator=%.17g maximum_denominator=%.17g",
           diagnostics_summary.minimum_denominator,
           diagnostics_summary.maximum_denominator);
    printf(" denominator_tail_ratio=%.17g log_contribution_span=%.17g",
           diagnostics_summary.denominator_tail_ratio,
           diagnostics_summary.log_contribution_span);
    printf(" hamiltonian_antihermitian_residual=%.17g",
           diagnostics_summary.hamiltonian_antihermitian_residual);
    printf(" hamiltonian_norm=%.17g\n",
           diagnostics_summary.hamiltonian_norm);

    if (!markov_print_tau("W", "real", w_series, (size_t)sample_count,
                          (double)official_block_length, maximum_tau_int,
                          block_length_multiplier, 1, 0,
                          &maximum_tau_seen, &tau_pass)) {
      goto cleanup;
    }

    for (family = MARKOV_FAMILY_S; family <= MARKOV_FAMILY_B; ++family) {
      size_t entry;
      for (entry = 0; entry < upper_count; ++entry) {
        size_t row = 0;
        size_t column = 0;
        MVMCKrylovJackknifeBlock official_blocks_for_entry
            [MARKOV_OFFICIAL_BLOCKS];
        MVMCKrylovJackknifeBlock diagnostic_blocks_for_entry
            [MARKOV_DIAGNOSTIC_BLOCKS];
        MVMCKrylovJackknifeResult official;
        MVMCKrylovJackknifeResult diagnostic;
        MVMCKrylovBlockStabilityResult stability;
        const MarkovEntryBudget *entry_budget =
            &exact_budget.entry[family][entry];
        double official_se_budget_ratio;
        double diagnostic_se_budget_ratio;
        double conservative_se_budget_ratio;
        double smaller_se;
        double larger_se;
        int pathology_pass;
        double complex theta;
        size_t index;
        if (!markov_entry_indices(upper_count, entry, &row, &column)) {
          goto cleanup;
        }
        status =
            mvmc_krylov_matrix_measurement_block_accumulator_extract_blocks(
                &official_accumulator, markov_matrix_kind(family), row,
                column, official_blocks_for_entry, official_block_count);
        if (status != MVMC_KRYLOV_STATUS_OK ||
            mvmc_krylov_matrix_measurement_block_accumulator_extract_blocks(
                &diagnostic_accumulator, markov_matrix_kind(family), row,
                column, diagnostic_blocks_for_entry,
                diagnostic_block_count) != MVMC_KRYLOV_STATUS_OK ||
            mvmc_krylov_jackknife_ratio(
                official_blocks_for_entry, official_block_count,
                minimum_leave_one_denominator, &official) !=
                MVMC_KRYLOV_STATUS_OK ||
            mvmc_krylov_jackknife_ratio(
                diagnostic_blocks_for_entry, diagnostic_block_count,
                minimum_leave_one_denominator, &diagnostic) !=
                MVMC_KRYLOV_STATUS_OK ||
            mvmc_krylov_block_stability_check(
                &official, &diagnostic, maximum_block_stability_ratio,
                &stability) != MVMC_KRYLOV_STATUS_OK) {
          goto cleanup;
        }
        official_se_budget_ratio = official.complex_se / entry_budget->budget;
        diagnostic_se_budget_ratio =
            diagnostic.complex_se / entry_budget->budget;
        conservative_se_budget_ratio =
            fmax(official.complex_se, diagnostic.complex_se) /
            entry_budget->budget;
        if (!isfinite(official_se_budget_ratio) ||
            !isfinite(diagnostic_se_budget_ratio) ||
            !isfinite(conservative_se_budget_ratio)) {
          goto cleanup;
        }
        if (official_se_budget_ratio > maximum_official_se_budget_ratio) {
          maximum_official_se_budget_ratio = official_se_budget_ratio;
        }
        if (conservative_se_budget_ratio > maximum_se_budget_ratio) {
          maximum_se_budget_ratio = conservative_se_budget_ratio;
        }
        if (stability.symmetric_se_ratio > maximum_stability_ratio) {
          maximum_stability_ratio = stability.symmetric_se_ratio;
        }
        smaller_se = fmin(official.complex_se, diagnostic.complex_se);
        larger_se = fmax(official.complex_se, diagnostic.complex_se);
        pathology_pass = !(larger_se > 0.0 && smaller_se == 0.0);
        if (!official.denominator_stable || !stability.passed ||
            !pathology_pass) {
          if (!official.denominator_stable) denominator_pass = 0;
          if (!stability.passed) block_stability_pass = 0;
          if (!pathology_pass) block_pathology_pass = 0;
        }
        if (conservative_se_budget_ratio >
            maximum_conservative_se_budget_ratio) {
          budget_pass = 0;
        }
        theta = official.estimate;
        for (index = 0; index < (size_t)sample_count; ++index) {
          const double complex numerator =
              entry_series[markov_series_index((size_t)sample_count,
                                               family, entry, upper_count,
                                               index)];
          const double complex delta = numerator - theta * w_series[index];
          scratch_series[index] = creal(delta);
        }
        {
          char name[32];
          snprintf(name, sizeof(name), "%s_%zu%zu",
                   markov_family_name(family), row, column);
          if (!markov_print_tau(name, "real", scratch_series,
                                (size_t)sample_count,
                                (double)official_block_length,
                                maximum_tau_int, block_length_multiplier,
                                1, 1, &maximum_tau_seen, &tau_pass)) {
            goto cleanup;
          }
          for (index = 0; index < (size_t)sample_count; ++index) {
            const double complex numerator =
                entry_series[markov_series_index((size_t)sample_count,
                                                 family, entry, upper_count,
                                                 index)];
            const double complex delta = numerator - theta * w_series[index];
            scratch_series[index] = cimag(delta);
          }
          if (!markov_print_tau(name, "imag", scratch_series,
                                (size_t)sample_count,
                                (double)official_block_length,
                                maximum_tau_int, block_length_multiplier,
                                1, 1, &maximum_tau_seen, &tau_pass)) {
            goto cleanup;
          }
        }
        printf("ENTRY family=%s row=%zu column=%zu",
               markov_family_name(family), row, column);
        printf(" exact_re=%.17g exact_im=%.17g",
               creal(entry_budget->exact_theta),
               cimag(entry_budget->exact_theta));
        printf(" estimate_re=%.17g estimate_im=%.17g",
               creal(official.estimate), cimag(official.estimate));
        printf(" official_se=%.17g diagnostic_se=%.17g",
               official.complex_se, diagnostic.complex_se);
        printf(" official_se_real=%.17g official_se_imag=%.17g",
               official.se_real, official.se_imag);
        printf(" covariance_real_imag=%.17g",
               official.covariance_real_imag);
        printf(" stability_ratio=%.17g stability_pass=%d",
               stability.symmetric_se_ratio, stability.passed);
        printf(" pathology_pass=%d", pathology_pass);
        printf(" denominator_stable=%d unstable_leave_one_blocks=%zu",
               official.denominator_stable,
               official.unstable_leave_one_blocks);
        printf(" budget=%.17g se_budget_ratio=%.17g",
               entry_budget->budget, official_se_budget_ratio);
        printf(" diagnostic_se_budget_ratio=%.17g",
               diagnostic_se_budget_ratio);
        printf(" conservative_se_budget_ratio=%.17g budget_pass=%d\n",
               conservative_se_budget_ratio,
               conservative_se_budget_ratio <=
                   maximum_conservative_se_budget_ratio);
      }
    }
    p4s_pass = trace_statistics.completed_steps == (uint64_t)sample_count &&
               trace_statistics.attempted_steps ==
               trace_statistics.completed_steps &&
               trace_statistics.support_violation_steps == 0 &&
               tau_pass && block_pathology_pass && budget_pass &&
               denominator_pass;
    status = mvmc_krylov_positive_sampler_manifest_create(
        &guide_policy, &proposal_policy, model, &limits, plan_hash,
        current_words, 1, &current, &rng, &final_manifest);
    if (status != MVMC_KRYLOV_STATUS_OK || !final_manifest.valid) {
      goto cleanup;
    }
    printf("MANIFEST rng_algorithm=%" PRIu64 " rng_state=%" PRIu64,
           final_manifest.rng_algorithm, final_manifest.rng_state);
    printf(" rng_stream=%" PRIu64 " rng_draws=%" PRIu64,
           final_manifest.rng_stream, final_manifest.rng_draws);
    printf(" policy_hash=%" PRIu64 " guide_shape_hash=%" PRIu64,
           final_manifest.policy_hash, final_manifest.guide_shape_hash);
    printf(" proposal_policy_hash=%" PRIu64,
           final_manifest.proposal_policy_hash);
    printf(" proposal_model_hash=%" PRIu64,
           final_manifest.proposal_model_hash);
    printf(" bounded_plan_hash=%" PRIu64,
           final_manifest.bounded_plan_hash);
    printf(" amplitude_policy_hash=%" PRIu64,
           final_manifest.amplitude_policy_hash);
    printf(" current_configuration_hash=%" PRIu64,
           final_manifest.current_configuration_hash);
    printf(" current_scale_hash=%" PRIu64,
           final_manifest.current_scale_hash);
    printf(" accepted_generation=%" PRIu64 "\n",
           final_manifest.accepted_generation);
    printf("DECISION p4s_decision=%s support_pass=%d tau_pass=%d",
           p4s_pass ? "GO" : "STOP",
           trace_statistics.support_violation_steps == 0, tau_pass);
    printf(" block_stability_pass=%d budget_pass=%d denominator_pass=%d",
           block_stability_pass, budget_pass, denominator_pass);
    printf(" block_pathology_pass=%d", block_pathology_pass);
    printf(" maximum_tau_int=%.17g maximum_se_budget_ratio=%.17g",
           maximum_tau_seen, maximum_se_budget_ratio);
    printf(" maximum_official_se_budget_ratio=%.17g",
           maximum_official_se_budget_ratio);
    printf(" maximum_block_stability_ratio=%.17g\n",
           maximum_stability_ratio);
  }
  return_code = 0;

cleanup:
  if (session_started && krylov_workspace != NULL) {
    (void)mvmc_bounded_krylov_session_end(krylov_workspace);
  }
  markov_free_block_storage(diagnostic_blocks, diagnostic_overlap,
                            diagnostic_hamiltonian,
                            diagnostic_hamiltonian_adjoint,
                            diagnostic_hamiltonian_squared);
  markov_free_block_storage(official_blocks, official_overlap,
                            official_hamiltonian,
                            official_hamiltonian_adjoint,
                            official_hamiltonian_squared);
  free(total_hamiltonian_squared);
  free(total_hamiltonian_adjoint);
  free(total_hamiltonian);
  free(total_overlap);
  free(entry_series);
  free(scratch_series);
  free(w_series);
  free(guide_weights);
  free(configurations);
  free(raw_values);
  free(exact_values);
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

#ifndef MVMC_BOUNDED_KRYLOV_MARKOV_DRIVER_NO_MAIN
int main(int argc, char **argv) {
  int site_count;
  int qp_total;
  int sample_count;
  size_t cache_bytes;
  double rho;
  uint64_t seed;
  MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
  uint64_t rng_stream = 0;
  int persistent_session = 0;
  int arguments_valid;
  int result;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#endif
  memset(&proposal_policy, 0, sizeof(proposal_policy));
  arguments_valid =
      argc >= 7 &&
      parse_int(argv[1], 4, PROFILE_MAX_SITE_COUNT, &site_count) &&
      (site_count & 1) == 0 && parse_int(argv[2], 1, 4, &qp_total) &&
      (qp_total == 1 || qp_total == 4) &&
      parse_int(argv[3], 1, 100000000, &sample_count) &&
      parse_size_arg(argv[4], (size_t)4 * 1024 * 1024 * 1024,
                     &cache_bytes) &&
      parse_markov_double(argv[5], DBL_MIN, DBL_MAX, &rho) &&
      parse_markov_u64(argv[6], &seed);
  if (arguments_valid && argc == 7) {
    arguments_valid =
        mvmc_krylov_positive_sampler_proposal_policy_create(
            0, 1, &proposal_policy) == MVMC_KRYLOV_STATUS_OK;
  } else if (arguments_valid &&
             (argc == 9 || argc == 10 || argc == 11)) {
    size_t global_numerator = 0;
    size_t global_denominator = 0;
    arguments_valid =
        parse_size_arg(argv[7], UINT64_MAX, &global_numerator) &&
        parse_size_arg(argv[8], UINT64_MAX, &global_denominator) &&
        mvmc_krylov_positive_sampler_proposal_policy_create(
            global_numerator, global_denominator, &proposal_policy) ==
            MVMC_KRYLOV_STATUS_OK &&
        (argc < 10 || parse_markov_u64(argv[9], &rng_stream)) &&
        (argc != 11 || strcmp(argv[10], "session") == 0);
    persistent_session = arguments_valid && argc == 11;
  } else if (arguments_valid && (argc == 12 || argc == 13) &&
             strcmp(argv[7], "shell") == 0) {
    size_t neighbor_numerator = 0;
    size_t neighbor_denominator = 0;
    size_t distance_numerator = 0;
    size_t distance_denominator = 0;
    arguments_valid =
        parse_size_arg(argv[8], UINT64_MAX, &neighbor_numerator) &&
        parse_size_arg(argv[9], UINT64_MAX, &neighbor_denominator) &&
        parse_size_arg(argv[10], UINT64_MAX, &distance_numerator) &&
        parse_size_arg(argv[11], UINT64_MAX, &distance_denominator) &&
        mvmc_krylov_positive_sampler_shell_policy_create(
            neighbor_numerator, neighbor_denominator, distance_numerator,
            distance_denominator, &proposal_policy) ==
            MVMC_KRYLOV_STATUS_OK &&
        (argc != 13 || parse_markov_u64(argv[12], &rng_stream));
  } else {
    arguments_valid = 0;
  }
  if (!arguments_valid) {
    if (world_rank() == 0) {
      fprintf(stderr,
              "usage: %s EVEN_SITE_COUNT(4..16) QP_TOTAL(1|4) "
              "SAMPLES CACHE_BYTES RHO SEED [GLOBAL_NUMERATOR "
              "GLOBAL_DENOMINATOR [RNG_STREAM [session]]]\n"
              "       %s EVEN_SITE_COUNT(4..16) QP_TOTAL(1|4) "
              "SAMPLES CACHE_BYTES RHO SEED shell NEIGHBOR_NUMERATOR "
              "NEIGHBOR_DENOMINATOR DISTANCE_NUMERATOR "
              "DISTANCE_DENOMINATOR [RNG_STREAM]\n",
              argv[0],
              argv[0]);
    }
    result = EXIT_FAILURE;
  } else {
    result = run_markov(site_count, qp_total, sample_count, cache_bytes,
                        rho, seed, &proposal_policy, rng_stream,
                        persistent_session) == 0
                 ? EXIT_SUCCESS
                 : EXIT_FAILURE;
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return result;
}
#endif
