/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
at your option any later version.
*/

#define MVMC_BOUNDED_KRYLOV_MARKOV_DRIVER_NO_MAIN
#include "bounded_krylov_markov_driver.c"

enum {
  PILOT_SCHEMA_VERSION = 1,
  PILOT_SHELL_SCHEMA_VERSION = 2,
  PILOT_FLAT_SCHEMA_VERSION = 3,
  PILOT_SURROGATE_SCHEMA_VERSION = 4,
  PILOT_PARTIAL_SCHEMA_VERSION = 5,
  PILOT_PARTIAL_NEIGHBOR_SCHEMA_VERSION = 6,
  PILOT_PARTIAL_CALLBACK_SCHEMA_VERSION = 7,
  PILOT_LONG_DIRECT_SCHEMA_VERSION = 8,
  PILOT_SEED_COUNT = 4,
  PILOT_COMPONENT_COUNT = 2,
  PILOT_SERIES_COUNT =
      MARKOV_FAMILY_COUNT * 6 * PILOT_COMPONENT_COUNT
};

typedef enum {
  PILOT_MODE_GLOBAL = 0,
  PILOT_MODE_SHELL,
  PILOT_MODE_FLAT_GUIDE,
  PILOT_MODE_SURROGATE,
  PILOT_MODE_PARTIAL_GUIDE,
  PILOT_MODE_PARTIAL_NEIGHBOR,
  PILOT_MODE_PARTIAL_CALLBACK,
  PILOT_MODE_LONG_DIRECT
} PilotMode;

typedef struct {
  MVMCKrylovBoundedWorkspace *workspace;
  ProfileScaledAmplitude *amplitude_workspace;
  const MVMCKrylovPositiveGuidePolicy *guide_policy;
  const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy;
  uint64_t evaluation_count;
  uint64_t terminal_amplitude_calls;
  double evaluation_seconds;
  double depth_seconds[MVMC_KRYLOV_MAX_ORDER + 1];
} PilotPartialCallbackContext;

typedef struct {
  double *guide_weights;
  MVMCKrylovPositiveGuideResult *guides;
  double *influence;
  double iid_variance[PILOT_SERIES_COUNT];
  double mean_denominator;
  double minimum_guide;
  double maximum_guide;
  double mean_guide;
  double guide_maximum_minimum_ratio;
  double uniform_sector_self_probability;
  double importance_weight_maximum_mean_ratio;
  double importance_weight_ess_fraction;
  MarkovExactBudget budget;
  uint64_t table_hash;
} PilotCandidateTable;

typedef struct {
  double *weights;
  double *log_weights;
  double *cumulative_weights;
  double total_weight;
  double minimum_weight;
  double maximum_weight;
  double mean_weight;
  double minimum_target_ratio;
  double maximum_target_ratio;
  double target_ratio_maximum_mean_ratio;
  double target_ratio_ess_fraction;
  double ideal_independence_acceptance;
  double ideal_independence_configuration_change;
  uint64_t table_hash;
} PilotSurrogateTable;

typedef struct {
  size_t current_index;
  uint64_t current_words;
  MVMCKrylovPositiveSamplerRng rng;
  uint64_t trace_hash;
  uint64_t accepted;
  uint64_t rejected;
  uint64_t neighbor_attempted;
  uint64_t global_attempted;
  uint64_t global_self;
  uint64_t shell_attempted;
  uint64_t surrogate_inner_attempted;
  uint64_t surrogate_inner_accepted;
  uint64_t surrogate_inner_rejected;
  uint64_t surrogate_final_changed;
  uint64_t surrogate_final_self;
  double surrogate_log_weight;
} PilotChain;

static size_t pilot_series_index(MarkovMatrixFamily family, size_t entry,
                                 int component) {
  return (((size_t)family * 6 + entry) * PILOT_COMPONENT_COUNT) +
         (size_t)component;
}

static uint64_t pilot_guide_shape_hash(
    const MVMCKrylovPositiveGuidePolicy *policy) {
  uint64_t hash = UINT64_C(1469598103934665603);
  int order;
  markov_hash_u64(&hash, policy->policy_hash);
  markov_hash_u64(&hash, (uint64_t)(unsigned int)policy->order);
  markov_hash_u64(&hash, markov_double_bits(policy->eta));
  for (order = 0; order <= policy->order; ++order) {
    markov_hash_u64(&hash, markov_double_bits(policy->lambda[order]));
    markov_hash_u64(&hash,
                    markov_double_bits(policy->log_basis_scale[order]));
  }
  return hash;
}

static void pilot_table_destroy(PilotCandidateTable *table) {
  if (table == NULL) return;
  free(table->influence);
  free(table->guides);
  free(table->guide_weights);
  memset(table, 0, sizeof(*table));
}

static void pilot_surrogate_table_destroy(PilotSurrogateTable *table) {
  if (table == NULL) return;
  free(table->log_weights);
  free(table->weights);
  free(table->cumulative_weights);
  memset(table, 0, sizeof(*table));
}

static MVMCKrylovStatus pilot_partial_log_weight_callback(
    const uint64_t *configuration_words, size_t word_count,
    void *context, double *log_weight) {
  PilotPartialCallbackContext *callback =
      (PilotPartialCallbackContext *)context;
  MVMCKrylovBoundedResult bounded;
  MVMCKrylovPositiveSamplerSurrogateWeightResult weight;
  MVMCKrylovStatus status;
  int order;
  if (configuration_words == NULL || word_count == 0 || callback == NULL ||
      callback->workspace == NULL || callback->amplitude_workspace == NULL ||
      callback->guide_policy == NULL ||
      callback->surrogate_policy == NULL || log_weight == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_bounded_krylov_evaluate(
      callback->workspace, configuration_words, word_count,
      profile_scaled_amplitude, callback->amplitude_workspace, &bounded);
  if (status != MVMC_KRYLOV_STATUS_OK || !bounded.valid ||
      bounded.evaluated_order != callback->surrogate_policy->partial_order) {
    return status == MVMC_KRYLOV_STATUS_OK
               ? MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE
               : status;
  }
  status = mvmc_krylov_positive_sampler_surrogate_weight_partial(
      callback->surrogate_policy, callback->guide_policy, bounded.value,
      (size_t)bounded.evaluated_order + 1, &weight);
  if (status != MVMC_KRYLOV_STATUS_OK || !weight.valid) return status;
  if (callback->evaluation_count == UINT64_MAX ||
      bounded.statistics.terminal_amplitude_calls >
          UINT64_MAX - callback->terminal_amplitude_calls) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  if (!isfinite(bounded.statistics.evaluation_wall_seconds) ||
      bounded.statistics.evaluation_wall_seconds < 0.0) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  for (order = 0; order <= callback->surrogate_policy->partial_order;
       ++order) {
    if (!isfinite(bounded.statistics.depth_wall_seconds[order]) ||
        bounded.statistics.depth_wall_seconds[order] < 0.0) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
  }
  ++callback->evaluation_count;
  callback->terminal_amplitude_calls +=
      bounded.statistics.terminal_amplitude_calls;
  callback->evaluation_seconds +=
      bounded.statistics.evaluation_wall_seconds;
  for (order = 0; order <= callback->surrogate_policy->partial_order;
       ++order) {
    callback->depth_seconds[order] +=
        bounded.statistics.depth_wall_seconds[order];
  }
  *log_weight = weight.log_weight;
  return MVMC_KRYLOV_STATUS_OK;
}

static int pilot_partial_callback_anchor(
    PilotPartialCallbackContext *callback,
    const uint64_t *configurations, size_t sector_size,
    const PilotSurrogateTable *partial_table,
    double *maximum_log_residual, double *maximum_weight_relative_residual,
    uint64_t *callback_table_hash) {
  size_t state;
  if (callback == NULL || configurations == NULL || sector_size == 0 ||
      partial_table == NULL || partial_table->log_weights == NULL ||
      maximum_log_residual == NULL ||
      maximum_weight_relative_residual == NULL ||
      callback_table_hash == NULL) {
    return 0;
  }
  *maximum_log_residual = 0.0;
  *maximum_weight_relative_residual = 0.0;
  *callback_table_hash = UINT64_C(1469598103934665603);
  for (state = 0; state < sector_size; ++state) {
    double log_weight = NAN;
    double log_residual;
    double relative_residual;
    if (pilot_partial_log_weight_callback(
            &configurations[state], 1, callback, &log_weight) !=
            MVMC_KRYLOV_STATUS_OK ||
        !isfinite(log_weight)) {
      return 0;
    }
    log_residual = fabs(log_weight - partial_table->log_weights[state]);
    relative_residual = fabs(expm1(
        log_weight - partial_table->log_weights[state]));
    if (!isfinite(log_residual) || !isfinite(relative_residual)) return 0;
    *maximum_log_residual = fmax(*maximum_log_residual, log_residual);
    *maximum_weight_relative_residual =
        fmax(*maximum_weight_relative_residual, relative_residual);
    markov_hash_u64(callback_table_hash, markov_double_bits(log_weight));
  }
  return callback->evaluation_count == (uint64_t)sector_size;
}

static int pilot_table_create(
    const MVMCKrylovMatrixMeasurementPolicy *measurement_policy,
    const MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT],
    size_t sector_size, PilotCandidateTable *table) {
  double *mean = NULL;
  double *second_moment = NULL;
  double total_guide = 0.0;
  double weighted_denominator_second_moment = 0.0;
  double maximum_denominator = 0.0;
  size_t state;
  size_t series;
  if (measurement_policy == NULL || exact_values == NULL || table == NULL ||
      sector_size == 0 ||
      sector_size > SIZE_MAX / PILOT_SERIES_COUNT ||
      sector_size * PILOT_SERIES_COUNT > SIZE_MAX / sizeof(double)) {
    return 0;
  }
  memset(table, 0, sizeof(*table));
  table->minimum_guide = DBL_MAX;
  table->guide_weights =
      (double *)calloc(sector_size, sizeof(*table->guide_weights));
  table->guides = (MVMCKrylovPositiveGuideResult *)calloc(
      sector_size, sizeof(*table->guides));
  table->influence = (double *)calloc(
      sector_size * PILOT_SERIES_COUNT, sizeof(*table->influence));
  mean = (double *)calloc(PILOT_SERIES_COUNT, sizeof(*mean));
  second_moment =
      (double *)calloc(PILOT_SERIES_COUNT, sizeof(*second_moment));
  if (table->guide_weights == NULL || table->guides == NULL ||
      table->influence == NULL || mean == NULL || second_moment == NULL ||
      !markov_compute_exact_budget(
          measurement_policy, exact_values, sector_size, 6,
          &table->budget, table->guide_weights)) {
    free(second_moment);
    free(mean);
    pilot_table_destroy(table);
    return 0;
  }
  table->table_hash = UINT64_C(1469598103934665603);
  for (state = 0; state < sector_size; ++state) {
    double complex overlap[6];
    double complex hamiltonian[6];
    double complex hamiltonian_adjoint[6];
    double complex hamiltonian_squared[6];
    MVMCKrylovMatrixMeasurementSampleDiagnostics diagnostics;
    MarkovMatrixFamily family;
    const double guide = table->guide_weights[state];
    if (mvmc_krylov_matrix_measurement_sample_with_adjoint(
            measurement_policy, exact_values[state], PROFILE_DEPTH_COUNT,
            overlap, hamiltonian, hamiltonian_adjoint,
            hamiltonian_squared, 6, &diagnostics) !=
            MVMC_KRYLOV_STATUS_OK ||
        !diagnostics.valid || !isfinite(guide) || guide <= 0.0 ||
        !isfinite(diagnostics.denominator) || diagnostics.denominator <= 0.0) {
      free(second_moment);
      free(mean);
      pilot_table_destroy(table);
      return 0;
    }
    table->guides[state].valid = 1;
    table->guides[state].status = MVMC_KRYLOV_STATUS_OK;
    table->guides[state].log_guide = log(guide);
    total_guide += guide;
    table->minimum_guide = fmin(table->minimum_guide, guide);
    table->maximum_guide = fmax(table->maximum_guide, guide);
    maximum_denominator = fmax(maximum_denominator,
                               diagnostics.denominator);
    weighted_denominator_second_moment +=
        guide * diagnostics.denominator * diagnostics.denominator;
    for (family = MARKOV_FAMILY_S; family <= MARKOV_FAMILY_B; ++family) {
      double complex *values = markov_family_values(
          family, overlap, hamiltonian, hamiltonian_squared);
      size_t entry;
      if (values == NULL) {
        free(second_moment);
        free(mean);
        pilot_table_destroy(table);
        return 0;
      }
      for (entry = 0; entry < 6; ++entry) {
        const double complex delta =
            values[entry] -
            table->budget.entry[family][entry].exact_theta *
                diagnostics.denominator;
        const size_t real_series = pilot_series_index(family, entry, 0);
        const size_t imag_series = pilot_series_index(family, entry, 1);
        table->influence[state * PILOT_SERIES_COUNT + real_series] =
            creal(delta);
        table->influence[state * PILOT_SERIES_COUNT + imag_series] =
            cimag(delta);
      }
    }
    markov_hash_u64(&table->table_hash, markov_double_bits(guide));
    markov_hash_u64(&table->table_hash,
                    markov_double_bits(diagnostics.denominator));
  }
  if (!isfinite(total_guide) || total_guide <= 0.0 ||
      !isfinite(weighted_denominator_second_moment) ||
      weighted_denominator_second_moment <= 0.0 ||
      !isfinite(table->minimum_guide) || table->minimum_guide <= 0.0 ||
      !isfinite(table->maximum_guide) || table->maximum_guide <= 0.0) {
    free(second_moment);
    free(mean);
    pilot_table_destroy(table);
    return 0;
  }
  table->mean_denominator = table->budget.target_sum / total_guide;
  table->mean_guide = total_guide / (double)sector_size;
  table->guide_maximum_minimum_ratio =
      table->maximum_guide / table->minimum_guide;
  table->uniform_sector_self_probability = 1.0 / (double)sector_size;
  table->importance_weight_maximum_mean_ratio =
      maximum_denominator / table->mean_denominator;
  table->importance_weight_ess_fraction =
      (table->budget.target_sum * table->budget.target_sum) /
      (total_guide * weighted_denominator_second_moment);
  for (state = 0; state < sector_size; ++state) {
    const double probability = table->guide_weights[state] / total_guide;
    for (series = 0; series < PILOT_SERIES_COUNT; ++series) {
      const double value =
          table->influence[state * PILOT_SERIES_COUNT + series];
      mean[series] += probability * value;
      second_moment[series] += probability * value * value;
    }
  }
  for (series = 0; series < PILOT_SERIES_COUNT; ++series) {
    double variance = second_moment[series] - mean[series] * mean[series];
    const double scale = fabs(second_moment[series]) +
                         fabs(mean[series] * mean[series]);
    if (variance < 0.0 && fabs(variance) <= 256.0 * DBL_EPSILON * scale) {
      variance = 0.0;
    }
    if (!isfinite(variance) || variance < 0.0) {
      free(second_moment);
      free(mean);
      pilot_table_destroy(table);
      return 0;
    }
    table->iid_variance[series] = variance;
    markov_hash_u64(&table->table_hash, markov_double_bits(variance));
  }
  free(second_moment);
  free(mean);
  return isfinite(table->mean_denominator) &&
         table->mean_denominator > 0.0 && isfinite(table->mean_guide) &&
         table->mean_guide > 0.0 &&
         isfinite(table->guide_maximum_minimum_ratio) &&
         table->guide_maximum_minimum_ratio >= 1.0 &&
         isfinite(table->importance_weight_maximum_mean_ratio) &&
         table->importance_weight_maximum_mean_ratio >= 1.0 &&
         isfinite(table->importance_weight_ess_fraction) &&
         table->importance_weight_ess_fraction > 0.0 &&
         table->importance_weight_ess_fraction <= 1.0 + 256.0 * DBL_EPSILON;
}

static int pilot_iid_summary(const PilotCandidateTable *table,
                             double *maximum_ratio,
                             double *required_tau,
                             MarkovMatrixFamily *limiting_family,
                             size_t *limiting_entry) {
  MarkovMatrixFamily family;
  if (table == NULL || maximum_ratio == NULL || required_tau == NULL ||
      limiting_family == NULL || limiting_entry == NULL ||
      !isfinite(table->mean_denominator) || table->mean_denominator <= 0.0) {
    return 0;
  }
  *maximum_ratio = -1.0;
  *required_tau = NAN;
  *limiting_family = MARKOV_FAMILY_S;
  *limiting_entry = 0;
  for (family = MARKOV_FAMILY_S; family <= MARKOV_FAMILY_B; ++family) {
    size_t entry;
    for (entry = 0; entry < 6; ++entry) {
      const double variance =
          table->iid_variance[pilot_series_index(family, entry, 0)] +
          table->iid_variance[pilot_series_index(family, entry, 1)];
      const double budget = table->budget.entry[family][entry].budget;
      const double ratio =
          sqrt(variance /
               (4096.0 * table->mean_denominator *
                table->mean_denominator)) /
          budget;
      if (!isfinite(ratio) || ratio < 0.0 || !isfinite(budget) ||
          budget <= 0.0) {
        return 0;
      }
      if (ratio > *maximum_ratio) {
        *maximum_ratio = ratio;
        *limiting_family = family;
        *limiting_entry = entry;
      }
    }
  }
  if (*maximum_ratio <= 0.0) return 0;
  *required_tau = (0.90 / *maximum_ratio) * (0.90 / *maximum_ratio);
  return isfinite(*required_tau) && *required_tau > 0.0;
}

static int pilot_surrogate_table_create(
    const MVMCKrylovPositiveGuidePolicy *guide_policy,
    const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy,
    const MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT],
    const PilotCandidateTable *target_table, size_t sector_size,
    int compute_ideal_independence, PilotSurrogateTable *table) {
  double total_weight = 0.0;
  double total_target = 0.0;
  double weighted_ratio_second = 0.0;
  double maximum_ratio = 0.0;
  size_t state;
  if (guide_policy == NULL || surrogate_policy == NULL ||
      exact_values == NULL || target_table == NULL || table == NULL ||
      sector_size == 0) {
    return 0;
  }
  memset(table, 0, sizeof(*table));
  table->minimum_weight = DBL_MAX;
  table->minimum_target_ratio = DBL_MAX;
  table->weights =
      (double *)calloc(sector_size, sizeof(*table->weights));
  table->log_weights =
      (double *)calloc(sector_size, sizeof(*table->log_weights));
  table->cumulative_weights =
      (double *)calloc(sector_size, sizeof(*table->cumulative_weights));
  if (table->weights == NULL || table->log_weights == NULL ||
      table->cumulative_weights == NULL) {
    pilot_surrogate_table_destroy(table);
    return 0;
  }
  table->table_hash = UINT64_C(1469598103934665603);
  for (state = 0; state < sector_size; ++state) {
    MVMCKrylovPositiveSamplerSurrogateWeightResult result;
    double weight;
    double ratio;
    if (mvmc_krylov_positive_sampler_surrogate_weight_zeroth(
            surrogate_policy, guide_policy, &exact_values[state][0],
            &result) != MVMC_KRYLOV_STATUS_OK || !result.valid) {
      pilot_surrogate_table_destroy(table);
      return 0;
    }
    weight = exp(result.log_weight);
    ratio = target_table->guide_weights[state] / weight;
    if (!isfinite(weight) || weight <= 0.0 || !isfinite(ratio) ||
        ratio <= 0.0) {
      pilot_surrogate_table_destroy(table);
      return 0;
    }
    table->weights[state] = weight;
    table->log_weights[state] = result.log_weight;
    total_weight += weight;
    table->cumulative_weights[state] = total_weight;
    total_target += target_table->guide_weights[state];
    weighted_ratio_second += weight * ratio * ratio;
    table->minimum_weight = fmin(table->minimum_weight, weight);
    table->maximum_weight = fmax(table->maximum_weight, weight);
    table->minimum_target_ratio =
        fmin(table->minimum_target_ratio, ratio);
    table->maximum_target_ratio =
        fmax(table->maximum_target_ratio, ratio);
    maximum_ratio = fmax(maximum_ratio, ratio);
    markov_hash_u64(&table->table_hash, markov_double_bits(weight));
    markov_hash_u64(&table->table_hash, markov_double_bits(ratio));
  }
  if (!isfinite(total_weight) || total_weight <= 0.0 ||
      !isfinite(total_target) || total_target <= 0.0 ||
      !isfinite(weighted_ratio_second) || weighted_ratio_second <= 0.0) {
    pilot_surrogate_table_destroy(table);
    return 0;
  }
  table->mean_weight = total_weight / (double)sector_size;
  table->total_weight = total_weight;
  table->target_ratio_maximum_mean_ratio =
      maximum_ratio / (total_target / total_weight);
  table->target_ratio_ess_fraction =
      (total_target * total_target) /
      (total_weight * weighted_ratio_second);
  if (compute_ideal_independence) {
    size_t x;
    for (x = 0; x < sector_size; ++x) {
      const double target_x = target_table->guide_weights[x];
      const double surrogate_x = table->weights[x];
      const double pi_x = target_x / total_target;
      size_t y;
      for (y = 0; y < sector_size; ++y) {
        const double target_y = target_table->guide_weights[y];
        const double surrogate_y = table->weights[y];
        const double q_y = surrogate_y / total_weight;
        const double acceptance =
            fmin(1.0, (target_y * surrogate_x) /
                           (target_x * surrogate_y));
        const double contribution = pi_x * q_y * acceptance;
        table->ideal_independence_acceptance += contribution;
        if (x != y) {
          table->ideal_independence_configuration_change += contribution;
        }
      }
    }
  } else {
    table->ideal_independence_acceptance = NAN;
    table->ideal_independence_configuration_change = NAN;
  }
  return isfinite(table->minimum_weight) && table->minimum_weight > 0.0 &&
         isfinite(table->maximum_weight) &&
         table->maximum_weight >= table->minimum_weight &&
         isfinite(table->minimum_target_ratio) &&
         table->minimum_target_ratio > 0.0 &&
         isfinite(table->maximum_target_ratio) &&
         table->maximum_target_ratio >= table->minimum_target_ratio &&
         isfinite(table->target_ratio_maximum_mean_ratio) &&
         table->target_ratio_maximum_mean_ratio >= 1.0 &&
         isfinite(table->target_ratio_ess_fraction) &&
         table->target_ratio_ess_fraction > 0.0 &&
         table->target_ratio_ess_fraction <= 1.0 + 256.0 * DBL_EPSILON &&
         (!compute_ideal_independence ||
          (isfinite(table->ideal_independence_acceptance) &&
           table->ideal_independence_acceptance > 0.0 &&
           table->ideal_independence_acceptance <=
               1.0 + 256.0 * DBL_EPSILON &&
           isfinite(table->ideal_independence_configuration_change) &&
           table->ideal_independence_configuration_change >= 0.0 &&
           table->ideal_independence_configuration_change <=
               table->ideal_independence_acceptance +
                   256.0 * DBL_EPSILON));
}

static uint64_t pilot_partial_policy_hash(
    const MVMCKrylovPositiveGuidePolicy *target_policy,
    int partial_order, size_t floor_multiplier) {
  uint64_t hash = UINT64_C(1469598103934665603);
  markov_hash_u64(&hash, UINT64_C(0x5034533750415254));
  markov_hash_u64(&hash, pilot_guide_shape_hash(target_policy));
  markov_hash_u64(&hash, (uint64_t)(unsigned int)partial_order);
  markov_hash_u64(&hash, (uint64_t)floor_multiplier);
  markov_hash_u64(&hash, UINT64_C(1));
  return hash;
}

static uint64_t pilot_partial_neighbor_policy_hash(
    uint64_t partial_policy_hash, size_t step_count) {
  uint64_t hash = UINT64_C(1469598103934665603);
  markov_hash_u64(&hash, UINT64_C(0x50345337504b4d48));
  markov_hash_u64(&hash, partial_policy_hash);
  markov_hash_u64(&hash, (uint64_t)step_count);
  markov_hash_u64(&hash, UINT64_C(1));
  return hash;
}

static int pilot_partial_table_create(
    const MVMCKrylovPositiveGuidePolicy *target_policy,
    int partial_order, size_t floor_multiplier,
    const MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT],
    const double complex (*raw_values)[PROFILE_DEPTH_COUNT],
    const PilotCandidateTable *target_table, size_t sector_size,
    PilotSurrogateTable *table, double *maximum_anchor_relative_residual,
    uint64_t *partial_policy_hash) {
  MVMCKrylovPositiveGuidePolicy partial_policy;
  double total_weight = 0.0;
  double total_target = 0.0;
  double weighted_ratio_second = 0.0;
  double maximum_ratio = 0.0;
  size_t state;
  int order;
  if (target_policy == NULL || exact_values == NULL || raw_values == NULL ||
      target_table == NULL || table == NULL ||
      maximum_anchor_relative_residual == NULL ||
      partial_policy_hash == NULL || partial_order < 1 ||
      partial_order >= target_policy->order || floor_multiplier == 0 ||
      sector_size == 0) {
    return 0;
  }
  memset(&partial_policy, 0, sizeof(partial_policy));
  partial_policy.order = partial_order;
  partial_policy.eta =
      (double)floor_multiplier * target_policy->eta;
  for (order = 0; order <= partial_order; ++order) {
    partial_policy.lambda[order] = 1.0;
    partial_policy.log_basis_scale[order] =
        target_policy->log_basis_scale[order];
  }
  partial_policy.policy_hash =
      pilot_partial_policy_hash(target_policy, partial_order,
                                floor_multiplier);
  if (!isfinite(partial_policy.eta) || partial_policy.eta <= 0.0 ||
      partial_policy.policy_hash == 0) {
    return 0;
  }
  memset(table, 0, sizeof(*table));
  table->minimum_weight = DBL_MAX;
  table->minimum_target_ratio = DBL_MAX;
  table->weights =
      (double *)calloc(sector_size, sizeof(*table->weights));
  table->log_weights =
      (double *)calloc(sector_size, sizeof(*table->log_weights));
  table->cumulative_weights =
      (double *)calloc(sector_size, sizeof(*table->cumulative_weights));
  if (table->weights == NULL || table->log_weights == NULL ||
      table->cumulative_weights == NULL) {
    pilot_surrogate_table_destroy(table);
    return 0;
  }
  table->table_hash = UINT64_C(1469598103934665603);
  *maximum_anchor_relative_residual = 0.0;
  *partial_policy_hash = partial_policy.policy_hash;
  for (state = 0; state < sector_size; ++state) {
    MVMCKrylovPositiveGuideResult result;
    double raw_weight = partial_policy.eta;
    double weight;
    double ratio;
    if (mvmc_krylov_positive_guide_evaluate(
            &partial_policy, exact_values[state], PROFILE_DEPTH_COUNT,
            &result) != MVMC_KRYLOV_STATUS_OK || !result.valid) {
      pilot_surrogate_table_destroy(table);
      return 0;
    }
    for (order = 0; order <= partial_order; ++order) {
      const double scaled_abs =
          cabs(raw_values[state][order]) *
          exp(target_policy->log_basis_scale[order]);
      raw_weight += scaled_abs * scaled_abs;
    }
    weight = exp(result.log_guide);
    ratio = target_table->guide_weights[state] / weight;
    if (!isfinite(raw_weight) || raw_weight <= 0.0 || !isfinite(weight) ||
        weight <= 0.0 || !isfinite(ratio) || ratio <= 0.0) {
      pilot_surrogate_table_destroy(table);
      return 0;
    }
    *maximum_anchor_relative_residual = fmax(
        *maximum_anchor_relative_residual,
        fabs(weight - raw_weight) / fmax(weight, raw_weight));
    table->weights[state] = weight;
    table->log_weights[state] = result.log_guide;
    total_weight += weight;
    table->cumulative_weights[state] = total_weight;
    total_target += target_table->guide_weights[state];
    weighted_ratio_second += weight * ratio * ratio;
    table->minimum_weight = fmin(table->minimum_weight, weight);
    table->maximum_weight = fmax(table->maximum_weight, weight);
    table->minimum_target_ratio =
        fmin(table->minimum_target_ratio, ratio);
    table->maximum_target_ratio =
        fmax(table->maximum_target_ratio, ratio);
    maximum_ratio = fmax(maximum_ratio, ratio);
    markov_hash_u64(&table->table_hash, markov_double_bits(weight));
    markov_hash_u64(&table->table_hash, markov_double_bits(ratio));
  }
  if (!isfinite(total_weight) || total_weight <= 0.0 ||
      !isfinite(total_target) || total_target <= 0.0 ||
      !isfinite(weighted_ratio_second) || weighted_ratio_second <= 0.0) {
    pilot_surrogate_table_destroy(table);
    return 0;
  }
  table->mean_weight = total_weight / (double)sector_size;
  table->total_weight = total_weight;
  table->target_ratio_maximum_mean_ratio =
      maximum_ratio / (total_target / total_weight);
  table->target_ratio_ess_fraction =
      (total_target * total_target) /
      (total_weight * weighted_ratio_second);
  for (state = 0; state < sector_size; ++state) {
    const double target_x = target_table->guide_weights[state];
    const double surrogate_x = table->weights[state];
    const double pi_x = target_x / total_target;
    size_t proposal;
    for (proposal = 0; proposal < sector_size; ++proposal) {
      const double target_y = target_table->guide_weights[proposal];
      const double surrogate_y = table->weights[proposal];
      const double q_y = surrogate_y / total_weight;
      const double acceptance =
          fmin(1.0, (target_y * surrogate_x) /
                         (target_x * surrogate_y));
      const double contribution = pi_x * q_y * acceptance;
      table->ideal_independence_acceptance += contribution;
      if (state != proposal) {
        table->ideal_independence_configuration_change += contribution;
      }
    }
  }
  return isfinite(table->minimum_weight) && table->minimum_weight > 0.0 &&
         isfinite(table->maximum_weight) &&
         table->maximum_weight >= table->minimum_weight &&
         isfinite(table->target_ratio_maximum_mean_ratio) &&
         table->target_ratio_maximum_mean_ratio >= 1.0 &&
         isfinite(table->target_ratio_ess_fraction) &&
         table->target_ratio_ess_fraction > 0.0 &&
         table->target_ratio_ess_fraction <= 1.0 + 256.0 * DBL_EPSILON &&
         isfinite(table->ideal_independence_acceptance) &&
         table->ideal_independence_acceptance > 0.0 &&
         table->ideal_independence_acceptance <=
             1.0 + 256.0 * DBL_EPSILON &&
         isfinite(table->ideal_independence_configuration_change) &&
         table->ideal_independence_configuration_change >= 0.0 &&
         table->ideal_independence_configuration_change <=
             table->ideal_independence_acceptance +
                 256.0 * DBL_EPSILON &&
         isfinite(*maximum_anchor_relative_residual) &&
         *maximum_anchor_relative_residual <= 1024.0 * DBL_EPSILON;
}

static int pilot_lookup_create(const uint64_t *configurations,
                               size_t sector_size, int site_count,
                               size_t **lookup, size_t *lookup_count) {
  size_t state;
  size_t count;
  if (configurations == NULL || lookup == NULL || lookup_count == NULL ||
      site_count <= 0 || site_count > 8) {
    return 0;
  }
  count = (size_t)1 << (2 * (unsigned int)site_count);
  *lookup = (size_t *)malloc(count * sizeof(**lookup));
  if (*lookup == NULL) return 0;
  for (state = 0; state < count; ++state) (*lookup)[state] = SIZE_MAX;
  for (state = 0; state < sector_size; ++state) {
    if (configurations[state] >= count ||
        (*lookup)[configurations[state]] != SIZE_MAX) {
      free(*lookup);
      *lookup = NULL;
      return 0;
    }
    (*lookup)[configurations[state]] = state;
  }
  *lookup_count = count;
  return 1;
}

static int pilot_chain_initialize(uint64_t seed,
                                  const double *guide_weights,
                                  size_t sector_size,
                                  const uint64_t *configurations,
                                  PilotChain *chain) {
  double uniform = NAN;
  if (guide_weights == NULL || configurations == NULL || chain == NULL) {
    return 0;
  }
  memset(chain, 0, sizeof(*chain));
  if (mvmc_krylov_positive_sampler_rng_seed(seed, UINT64_C(0),
                                             &chain->rng) !=
          MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_rng_draw_uniform(
          &chain->rng, &uniform) != MVMC_KRYLOV_STATUS_OK) {
    return 0;
  }
  chain->current_index =
      markov_inverse_cdf(uniform, guide_weights, sector_size, NULL);
  chain->current_words = configurations[chain->current_index];
  chain->trace_hash = UINT64_C(1469598103934665603);
  chain->surrogate_log_weight = NAN;
  return 1;
}

static int pilot_chain_advance(
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const PilotCandidateTable *table, const size_t *lookup,
    size_t lookup_count, size_t sample_count,
    double *series, PilotChain *chain) {
  size_t sample;
  for (sample = 0; sample < sample_count; ++sample) {
    MVMCKrylovPositiveSamplerProposalDrawResult draw;
    MVMCKrylovPositiveGuideAcceptance acceptance;
    uint64_t proposal_words = 0;
    size_t proposal_index;
    size_t series_index;
    MVMCKrylovStatus status =
        mvmc_krylov_positive_sampler_draw_mixture_rng(
            model, proposal_policy, &chain->current_words, 1,
            &chain->rng, &proposal_words, 1, &draw);
    if (status != MVMC_KRYLOV_STATUS_OK || !draw.valid ||
        proposal_words >= lookup_count ||
        lookup[proposal_words] == SIZE_MAX) {
      return 0;
    }
    proposal_index = lookup[proposal_words];
    status = mvmc_krylov_positive_guide_acceptance(
        &table->guides[chain->current_index],
        &table->guides[proposal_index], draw.log_proposal_ratio,
        draw.uniform, &acceptance);
    if (status != MVMC_KRYLOV_STATUS_OK || !acceptance.valid) return 0;
    if (draw.component == MVMC_KRYLOV_PROPOSAL_COMPONENT_GLOBAL) {
      ++chain->global_attempted;
      if (draw.self_proposal) ++chain->global_self;
    } else if (draw.component ==
               MVMC_KRYLOV_PROPOSAL_COMPONENT_NEIGHBOR) {
      ++chain->neighbor_attempted;
    } else if (draw.component == MVMC_KRYLOV_PROPOSAL_COMPONENT_SHELL) {
      if (draw.self_proposal || draw.shell_distance == 0 ||
          draw.shell_count == 0) {
        return 0;
      }
      ++chain->shell_attempted;
    } else {
      return 0;
    }
    if (acceptance.accepted) {
      chain->current_index = proposal_index;
      chain->current_words = proposal_words;
      ++chain->accepted;
    } else {
      ++chain->rejected;
    }
    chain->rng = draw.rng_after;
    markov_hash_u64(&chain->trace_hash, chain->current_words);
    markov_hash_u64(&chain->trace_hash,
                    (uint64_t)(unsigned int)draw.component);
    markov_hash_u64(&chain->trace_hash,
                    (uint64_t)(unsigned int)acceptance.accepted);
    markov_hash_u64(&chain->trace_hash, chain->rng.state);
    if (series != NULL) {
      for (series_index = 0; series_index < PILOT_SERIES_COUNT;
           ++series_index) {
        series[series_index * sample_count + sample] =
            table->influence[chain->current_index * PILOT_SERIES_COUNT +
                             series_index];
      }
    }
  }
  return 1;
}

static int pilot_chain_advance_surrogate(
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy,
    const PilotCandidateTable *target_table,
    const PilotSurrogateTable *surrogate_table, const size_t *lookup,
    size_t lookup_count, size_t sample_count, double *series,
    PilotChain *chain) {
  size_t sample;
  if (model == NULL || surrogate_policy == NULL || target_table == NULL ||
      surrogate_table == NULL || lookup == NULL || chain == NULL) {
    return 0;
  }
  for (sample = 0; sample < sample_count; ++sample) {
    const size_t outer_start_index = chain->current_index;
    size_t scratch_index = outer_start_index;
    uint64_t scratch_words = chain->current_words;
    size_t inner_step;
    double outer_uniform = NAN;
    MVMCKrylovPositiveGuideAcceptance outer_acceptance;
    for (inner_step = 0; inner_step < surrogate_policy->step_count;
         ++inner_step) {
      uint64_t proposal_words = 0;
      size_t proposal_index;
      size_t neighbor_count = 0;
      size_t selected_neighbor_index = 0;
      double neighbor_log_proposal_ratio = NAN;
      double uniform = NAN;
      int accepted;
      if (mvmc_krylov_fock_proposal_count_neighbors(
              model, &scratch_words, 1, &neighbor_count) !=
              MVMC_KRYLOV_STATUS_OK ||
          neighbor_count == 0 ||
          mvmc_krylov_positive_sampler_rng_draw_bounded(
              &chain->rng, neighbor_count, &selected_neighbor_index) !=
              MVMC_KRYLOV_STATUS_OK ||
          mvmc_krylov_fock_proposal_select_neighbor(
              model, &scratch_words, 1, selected_neighbor_index,
              &proposal_words, 1, &neighbor_count) !=
              MVMC_KRYLOV_STATUS_OK ||
          mvmc_krylov_fock_proposal_log_ratio(
              model, &scratch_words, 1, &proposal_words, 1,
              &neighbor_log_proposal_ratio) != MVMC_KRYLOV_STATUS_OK ||
          neighbor_log_proposal_ratio != 0.0 ||
          proposal_words >= lookup_count ||
          lookup[proposal_words] == SIZE_MAX ||
          mvmc_krylov_positive_sampler_rng_draw_uniform(
              &chain->rng, &uniform) != MVMC_KRYLOV_STATUS_OK) {
        return 0;
      }
      proposal_index = lookup[proposal_words];
      accepted =
          surrogate_table->log_weights[proposal_index] >=
              surrogate_table->log_weights[scratch_index] ||
          uniform == 0.0 ||
          log(uniform) < surrogate_table->log_weights[proposal_index] -
                             surrogate_table->log_weights[scratch_index];
      if (chain->surrogate_inner_attempted == UINT64_MAX ||
          (accepted && chain->surrogate_inner_accepted == UINT64_MAX) ||
          (!accepted && chain->surrogate_inner_rejected == UINT64_MAX)) {
        return 0;
      }
      ++chain->surrogate_inner_attempted;
      if (accepted) {
        scratch_index = proposal_index;
        scratch_words = proposal_words;
        ++chain->surrogate_inner_accepted;
      } else {
        ++chain->surrogate_inner_rejected;
      }
    }
    if (mvmc_krylov_positive_sampler_rng_draw_uniform(
            &chain->rng, &outer_uniform) != MVMC_KRYLOV_STATUS_OK ||
        mvmc_krylov_positive_guide_acceptance(
            &target_table->guides[outer_start_index],
            &target_table->guides[scratch_index],
            surrogate_table->log_weights[outer_start_index] -
                surrogate_table->log_weights[scratch_index],
            outer_uniform, &outer_acceptance) != MVMC_KRYLOV_STATUS_OK ||
        !outer_acceptance.valid) {
      return 0;
    }
    if ((scratch_index == outer_start_index &&
         chain->surrogate_final_self == UINT64_MAX) ||
        (scratch_index != outer_start_index &&
         chain->surrogate_final_changed == UINT64_MAX) ||
        (outer_acceptance.accepted && chain->accepted == UINT64_MAX) ||
        (!outer_acceptance.accepted && chain->rejected == UINT64_MAX)) {
      return 0;
    }
    if (scratch_index == outer_start_index) {
      ++chain->surrogate_final_self;
    } else {
      ++chain->surrogate_final_changed;
    }
    if (outer_acceptance.accepted) {
      chain->current_index = scratch_index;
      chain->current_words = scratch_words;
      ++chain->accepted;
    } else {
      ++chain->rejected;
    }
    markov_hash_u64(&chain->trace_hash, chain->current_words);
    markov_hash_u64(&chain->trace_hash, surrogate_policy->policy_hash);
    markov_hash_u64(&chain->trace_hash,
                    (uint64_t)(unsigned int)outer_acceptance.accepted);
    markov_hash_u64(&chain->trace_hash, chain->rng.state);
    if (series != NULL) {
      size_t series_index;
      for (series_index = 0; series_index < PILOT_SERIES_COUNT;
           ++series_index) {
        series[series_index * sample_count + sample] =
            target_table->influence[
                chain->current_index * PILOT_SERIES_COUNT + series_index];
      }
    }
  }
  return 1;
}

static int pilot_chain_advance_surrogate_callback(
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy,
    const PilotCandidateTable *target_table,
    PilotPartialCallbackContext *callback, const size_t *lookup,
    size_t lookup_count, size_t sample_count, double *series,
    PilotChain *chain) {
  size_t sample;
  if (model == NULL || surrogate_policy == NULL || target_table == NULL ||
      callback == NULL || lookup == NULL || chain == NULL) {
    return 0;
  }
  if (!isfinite(chain->surrogate_log_weight) &&
      pilot_partial_log_weight_callback(
          &chain->current_words, 1, callback,
          &chain->surrogate_log_weight) != MVMC_KRYLOV_STATUS_OK) {
    return 0;
  }
  for (sample = 0; sample < sample_count; ++sample) {
    MVMCKrylovPositiveSamplerSurrogateDrawResult draw;
    MVMCKrylovPositiveGuideAcceptance outer_acceptance;
    uint64_t proposal_words = 0;
    size_t proposal_index;
    double outer_uniform = NAN;
    const size_t outer_start_index = chain->current_index;
    MVMCKrylovStatus status =
        mvmc_krylov_positive_sampler_draw_surrogate_rng(
            model, surrogate_policy, &chain->current_words, 1,
            chain->surrogate_log_weight, &chain->rng,
            pilot_partial_log_weight_callback, callback, &proposal_words, 1,
            &draw);
    if (status != MVMC_KRYLOV_STATUS_OK || !draw.valid ||
        proposal_words >= lookup_count || lookup[proposal_words] == SIZE_MAX) {
      return 0;
    }
    proposal_index = lookup[proposal_words];
    chain->rng = draw.rng_after;
    if (mvmc_krylov_positive_sampler_rng_draw_uniform(
            &chain->rng, &outer_uniform) != MVMC_KRYLOV_STATUS_OK ||
        mvmc_krylov_positive_guide_acceptance(
            &target_table->guides[outer_start_index],
            &target_table->guides[proposal_index], draw.log_proposal_ratio,
            outer_uniform, &outer_acceptance) != MVMC_KRYLOV_STATUS_OK ||
        !outer_acceptance.valid ||
        draw.surrogate_evaluation_count >
            UINT64_MAX - chain->surrogate_inner_attempted ||
        draw.inner_accepted >
            UINT64_MAX - chain->surrogate_inner_accepted ||
        draw.inner_rejected >
            UINT64_MAX - chain->surrogate_inner_rejected ||
        (draw.final_configuration_changed &&
         chain->surrogate_final_changed == UINT64_MAX) ||
        (!draw.final_configuration_changed &&
         chain->surrogate_final_self == UINT64_MAX) ||
        (outer_acceptance.accepted && chain->accepted == UINT64_MAX) ||
        (!outer_acceptance.accepted && chain->rejected == UINT64_MAX)) {
      return 0;
    }
    chain->surrogate_inner_attempted += draw.surrogate_evaluation_count;
    chain->surrogate_inner_accepted += draw.inner_accepted;
    chain->surrogate_inner_rejected += draw.inner_rejected;
    if (draw.final_configuration_changed) {
      ++chain->surrogate_final_changed;
    } else {
      ++chain->surrogate_final_self;
    }
    if (outer_acceptance.accepted) {
      chain->current_index = proposal_index;
      chain->current_words = proposal_words;
      chain->surrogate_log_weight = draw.final_log_weight;
      ++chain->accepted;
    } else {
      ++chain->rejected;
    }
    markov_hash_u64(&chain->trace_hash, chain->current_words);
    markov_hash_u64(&chain->trace_hash, surrogate_policy->policy_hash);
    markov_hash_u64(&chain->trace_hash,
                    (uint64_t)(unsigned int)outer_acceptance.accepted);
    markov_hash_u64(&chain->trace_hash, chain->rng.state);
    if (series != NULL) {
      size_t series_index;
      for (series_index = 0; series_index < PILOT_SERIES_COUNT;
           ++series_index) {
        series[series_index * sample_count + sample] =
            target_table->influence[
                chain->current_index * PILOT_SERIES_COUNT + series_index];
      }
    }
  }
  return 1;
}

static size_t pilot_weighted_index(
    double uniform, const PilotSurrogateTable *table, size_t sector_size) {
  const double threshold = uniform * table->total_weight;
  size_t lower = 0;
  size_t upper = sector_size;
  while (lower + 1 < upper) {
    const size_t middle = lower + (upper - lower) / 2;
    if (threshold < table->cumulative_weights[middle]) {
      upper = middle;
    } else {
      lower = middle;
    }
  }
  if (threshold < table->cumulative_weights[lower]) return lower;
  return upper < sector_size ? upper : sector_size - 1;
}

static int pilot_chain_advance_partial_independence(
    const PilotCandidateTable *target_table,
    const PilotSurrogateTable *partial_table,
    const uint64_t *configurations, size_t sector_size,
    uint64_t partial_policy_hash, size_t sample_count,
    double *series, PilotChain *chain) {
  size_t sample;
  if (target_table == NULL || partial_table == NULL ||
      configurations == NULL || sector_size == 0 ||
      partial_policy_hash == 0 || chain == NULL ||
      partial_table->cumulative_weights == NULL ||
      !isfinite(partial_table->total_weight) ||
      partial_table->total_weight <= 0.0) {
    return 0;
  }
  for (sample = 0; sample < sample_count; ++sample) {
    double proposal_uniform = NAN;
    double acceptance_uniform = NAN;
    size_t proposal_index;
    MVMCKrylovPositiveGuideAcceptance acceptance;
    if (mvmc_krylov_positive_sampler_rng_draw_uniform(
            &chain->rng, &proposal_uniform) != MVMC_KRYLOV_STATUS_OK ||
        mvmc_krylov_positive_sampler_rng_draw_uniform(
            &chain->rng, &acceptance_uniform) != MVMC_KRYLOV_STATUS_OK) {
      return 0;
    }
    proposal_index = pilot_weighted_index(
        proposal_uniform, partial_table, sector_size);
    if (proposal_index >= sector_size ||
        mvmc_krylov_positive_guide_acceptance(
            &target_table->guides[chain->current_index],
            &target_table->guides[proposal_index],
            partial_table->log_weights[chain->current_index] -
                partial_table->log_weights[proposal_index],
            acceptance_uniform, &acceptance) != MVMC_KRYLOV_STATUS_OK ||
        !acceptance.valid || chain->global_attempted == UINT64_MAX ||
        (proposal_index == chain->current_index &&
         chain->global_self == UINT64_MAX) ||
        (acceptance.accepted && chain->accepted == UINT64_MAX) ||
        (!acceptance.accepted && chain->rejected == UINT64_MAX)) {
      return 0;
    }
    ++chain->global_attempted;
    if (proposal_index == chain->current_index) ++chain->global_self;
    if (acceptance.accepted) {
      chain->current_index = proposal_index;
      chain->current_words = configurations[proposal_index];
      ++chain->accepted;
    } else {
      ++chain->rejected;
    }
    markov_hash_u64(&chain->trace_hash, chain->current_words);
    markov_hash_u64(&chain->trace_hash, partial_policy_hash);
    markov_hash_u64(&chain->trace_hash,
                    (uint64_t)(unsigned int)acceptance.accepted);
    markov_hash_u64(&chain->trace_hash, chain->rng.state);
    if (series != NULL) {
      size_t series_index;
      for (series_index = 0; series_index < PILOT_SERIES_COUNT;
           ++series_index) {
        series[series_index * sample_count + sample] =
            target_table->influence[
                chain->current_index * PILOT_SERIES_COUNT + series_index];
      }
    }
  }
  return 1;
}

static int pilot_configuration_distance(
    const MVMCKrylovFockModel *model, uint64_t current,
    uint64_t proposal, size_t *distance) {
  size_t site;
  size_t up_removed = 0;
  size_t down_removed = 0;
  if (model == NULL || distance == NULL) return 0;
  for (site = 0; site < model->site_count; ++site) {
    const uint64_t up_mask = UINT64_C(1) << site;
    const uint64_t down_mask = UINT64_C(1) << (model->site_count + site);
    if ((current & up_mask) != 0 && (proposal & up_mask) == 0) {
      ++up_removed;
    }
    if ((current & down_mask) != 0 && (proposal & down_mask) == 0) {
      ++down_removed;
    }
  }
  if (model->pure_spin) {
    if (up_removed != down_removed) return 0;
    *distance = up_removed;
  } else {
    if (up_removed > SIZE_MAX - down_removed) return 0;
    *distance = up_removed + down_removed;
  }
  return 1;
}

static int pilot_exact_l4_balance(
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const uint64_t *configurations, const double *guide_weights,
    size_t sector_size, const size_t *lookup, size_t lookup_count,
    double *maximum_proposal_row_residual,
    double *maximum_proposal_symmetry_residual,
    double *maximum_db_residual, double *maximum_stationary_residual,
    int *shell_cardinality_pass) {
  double *transition = NULL;
  double total_guide = 0.0;
  double maximum_flux = 0.0;
  double maximum_probability = 0.0;
  double maximum_proposal = 0.0;
  size_t resolved_maximum = 0;
  size_t resolved_distance = 0;
  size_t production_shell_count = 0;
  const int shell_kernel =
      proposal_policy->kernel_id ==
      MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_SHELL_MIXTURE;
  size_t x;
  (void)lookup;
  (void)lookup_count;
  if (maximum_proposal_row_residual == NULL ||
      maximum_proposal_symmetry_residual == NULL ||
      maximum_db_residual == NULL ||
      maximum_stationary_residual == NULL ||
      shell_cardinality_pass == NULL) {
    return 0;
  }
  if (model->site_count != 4) {
    *maximum_proposal_row_residual = 0.0;
    *maximum_proposal_symmetry_residual = 0.0;
    *maximum_db_residual = 0.0;
    *maximum_stationary_residual = 0.0;
    *shell_cardinality_pass = 1;
    return 1;
  }
  if (sector_size > SIZE_MAX / sector_size ||
      sector_size * sector_size > SIZE_MAX / sizeof(*transition)) {
    return 0;
  }
  transition =
      (double *)calloc(sector_size * sector_size, sizeof(*transition));
  if (transition == NULL) return 0;
  *maximum_proposal_row_residual = 0.0;
  *maximum_proposal_symmetry_residual = 0.0;
  *shell_cardinality_pass = 1;
  if (shell_kernel &&
      (mvmc_krylov_fock_proposal_resolve_shell_distance(
           model, proposal_policy->distance_numerator,
           proposal_policy->distance_denominator, &resolved_maximum,
           &resolved_distance) != MVMC_KRYLOV_STATUS_OK ||
       mvmc_krylov_fock_proposal_count_shell(
           model, resolved_distance, &production_shell_count) !=
           MVMC_KRYLOV_STATUS_OK ||
       resolved_maximum == 0 || resolved_distance == 0 ||
       production_shell_count == 0)) {
    free(transition);
    return 0;
  }
  for (x = 0; x < sector_size; ++x) total_guide += guide_weights[x];
  for (x = 0; x < sector_size; ++x) {
    size_t neighbor_count = 0;
    size_t y;
    size_t observed_neighbor_count = 0;
    size_t observed_shell_count = 0;
    double proposal_row_sum = 0.0;
    if (mvmc_krylov_fock_proposal_count_neighbors(
            model, &configurations[x], 1, &neighbor_count) !=
            MVMC_KRYLOV_STATUS_OK ||
        neighbor_count == 0) {
      free(transition);
      return 0;
    }
    for (y = 0; y < sector_size; ++y) {
      size_t distance = 0;
      if (!pilot_configuration_distance(
              model, configurations[x], configurations[y], &distance)) {
        free(transition);
        return 0;
      }
      if (distance == 1) ++observed_neighbor_count;
      if (shell_kernel && distance == resolved_distance) {
        ++observed_shell_count;
      }
    }
    if (observed_neighbor_count != neighbor_count ||
        (shell_kernel && observed_shell_count != production_shell_count)) {
      *shell_cardinality_pass = 0;
    }
    for (y = 0; y < sector_size; ++y) {
      size_t distance = 0;
      double probability = 0.0;
      if (!pilot_configuration_distance(
              model, configurations[x], configurations[y], &distance)) {
        free(transition);
        return 0;
      }
      if (proposal_policy->kernel_id ==
          MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_GLOBAL_MIXTURE) {
        const double global_weight =
            (double)proposal_policy->global_numerator /
            (double)proposal_policy->global_denominator;
        probability += global_weight / (double)sector_size;
        if (distance == 1) {
          probability += (1.0 - global_weight) / (double)neighbor_count;
        }
      } else if (shell_kernel) {
        const double neighbor_weight =
            (double)proposal_policy->neighbor_numerator /
            (double)proposal_policy->neighbor_denominator;
        if (distance == 1) {
          probability += neighbor_weight / (double)observed_neighbor_count;
        }
        if (distance == resolved_distance) {
          probability +=
              (1.0 - neighbor_weight) / (double)observed_shell_count;
        }
      } else {
        free(transition);
        return 0;
      }
      transition[x * sector_size + y] = probability;
      proposal_row_sum += probability;
      maximum_proposal = fmax(maximum_proposal, probability);
    }
    *maximum_proposal_row_residual =
        fmax(*maximum_proposal_row_residual, fabs(proposal_row_sum - 1.0));
  }
  for (x = 0; x < sector_size; ++x) {
    size_t y;
    for (y = 0; y < sector_size; ++y) {
      *maximum_proposal_symmetry_residual = fmax(
          *maximum_proposal_symmetry_residual,
          fabs(transition[x * sector_size + y] -
               transition[y * sector_size + x]));
    }
  }
  if (maximum_proposal <= 0.0) {
    free(transition);
    return 0;
  }
  *maximum_proposal_row_residual /= maximum_proposal;
  *maximum_proposal_symmetry_residual /= maximum_proposal;
  for (x = 0; x < sector_size; ++x) {
    size_t y;
    double off_diagonal_sum = 0.0;
    for (y = 0; y < sector_size; ++y) {
      if (y != x) {
        const double acceptance =
            fmin(1.0, guide_weights[y] / guide_weights[x]);
        transition[x * sector_size + y] *= acceptance;
        off_diagonal_sum += transition[x * sector_size + y];
      }
    }
    transition[x * sector_size + x] = 1.0 - off_diagonal_sum;
  }
  *maximum_db_residual = 0.0;
  *maximum_stationary_residual = 0.0;
  for (x = 0; x < sector_size; ++x) {
    const double qx = guide_weights[x] / total_guide;
    double stationary_sum = 0.0;
    size_t y;
    maximum_probability = fmax(maximum_probability, qx);
    for (y = 0; y < sector_size; ++y) {
      const double qy = guide_weights[y] / total_guide;
      const double forward = qx * transition[x * sector_size + y];
      const double reverse = qy * transition[y * sector_size + x];
      maximum_flux = fmax(maximum_flux, fmax(fabs(forward), fabs(reverse)));
      *maximum_db_residual =
          fmax(*maximum_db_residual, fabs(forward - reverse));
      stationary_sum +=
          qy * transition[y * sector_size + x];
    }
    *maximum_stationary_residual =
        fmax(*maximum_stationary_residual, fabs(stationary_sum - qx));
  }
  free(transition);
  if (maximum_flux <= 0.0 || maximum_probability <= 0.0) return 0;
  *maximum_db_residual /= maximum_flux;
  *maximum_stationary_residual /= maximum_probability;
  return isfinite(*maximum_db_residual) &&
         isfinite(*maximum_stationary_residual) &&
         isfinite(*maximum_proposal_row_residual) &&
         isfinite(*maximum_proposal_symmetry_residual) &&
         *shell_cardinality_pass &&
         *maximum_proposal_row_residual <= 256.0 * DBL_EPSILON &&
         *maximum_proposal_symmetry_residual <= 256.0 * DBL_EPSILON &&
         *maximum_db_residual <= 256.0 * DBL_EPSILON &&
         *maximum_stationary_residual <= 256.0 * DBL_EPSILON;
}

static void pilot_dense_transition_multiply(const double *left,
                                            const double *right,
                                            size_t count,
                                            double *product) {
  size_t row;
  for (row = 0; row < count; ++row) {
    size_t column;
    for (column = 0; column < count; ++column) {
      double value = 0.0;
      size_t inner;
      for (inner = 0; inner < count; ++inner) {
        value += left[row * count + inner] *
                 right[inner * count + column];
      }
      product[row * count + column] = value;
    }
  }
}

static int pilot_exact_l4_surrogate_balance(
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy,
    const uint64_t *configurations, const double *target_weights,
    const double *surrogate_weights, size_t sector_size,
    double *inner_row_residual, double *inner_db_residual,
    double *powered_row_residual, double *powered_db_residual,
    double *powered_ratio_residual, double *outer_db_residual,
    double *outer_stationary_residual) {
  double *inner = NULL;
  double *power = NULL;
  double *temporary = NULL;
  double *outer = NULL;
  double total_target = 0.0;
  double total_surrogate = 0.0;
  double maximum_inner_flux = 0.0;
  double maximum_powered_flux = 0.0;
  double maximum_outer_flux = 0.0;
  double maximum_target_probability = 0.0;
  size_t x;
  if (inner_row_residual == NULL || inner_db_residual == NULL ||
      powered_row_residual == NULL || powered_db_residual == NULL ||
      powered_ratio_residual == NULL || outer_db_residual == NULL ||
      outer_stationary_residual == NULL || model == NULL ||
      surrogate_policy == NULL || configurations == NULL ||
      target_weights == NULL || surrogate_weights == NULL) {
    return 0;
  }
  *inner_row_residual = 0.0;
  *inner_db_residual = 0.0;
  *powered_row_residual = 0.0;
  *powered_db_residual = 0.0;
  *powered_ratio_residual = 0.0;
  *outer_db_residual = 0.0;
  *outer_stationary_residual = 0.0;
  if (model->site_count != 4) return 1;
  if (sector_size > SIZE_MAX / sector_size ||
      sector_size * sector_size > SIZE_MAX / sizeof(*inner)) {
    return 0;
  }
  inner = (double *)calloc(sector_size * sector_size, sizeof(*inner));
  power = (double *)calloc(sector_size * sector_size, sizeof(*power));
  temporary =
      (double *)calloc(sector_size * sector_size, sizeof(*temporary));
  outer = (double *)calloc(sector_size * sector_size, sizeof(*outer));
  if (inner == NULL || power == NULL || temporary == NULL || outer == NULL) {
    free(outer);
    free(temporary);
    free(power);
    free(inner);
    return 0;
  }
  for (x = 0; x < sector_size; ++x) {
    size_t neighbor_count = 0;
    size_t y;
    double off_diagonal = 0.0;
    if (mvmc_krylov_fock_proposal_count_neighbors(
            model, &configurations[x], 1, &neighbor_count) !=
            MVMC_KRYLOV_STATUS_OK ||
        neighbor_count == 0) {
      goto fail;
    }
    for (y = 0; y < sector_size; ++y) {
      size_t distance = 0;
      if (!pilot_configuration_distance(model, configurations[x],
                                        configurations[y], &distance)) {
        goto fail;
      }
      if (distance == 1) {
        const double acceptance =
            fmin(1.0, surrogate_weights[y] / surrogate_weights[x]);
        inner[x * sector_size + y] =
            acceptance / (double)neighbor_count;
        off_diagonal += inner[x * sector_size + y];
      }
    }
    inner[x * sector_size + x] = 1.0 - off_diagonal;
    power[x * sector_size + x] = 1.0;
    total_target += target_weights[x];
    total_surrogate += surrogate_weights[x];
  }
  for (x = 0; x < sector_size; ++x) {
    size_t y;
    double row_sum = 0.0;
    for (y = 0; y < sector_size; ++y) {
      const double forward =
          surrogate_weights[x] * inner[x * sector_size + y];
      const double reverse =
          surrogate_weights[y] * inner[y * sector_size + x];
      row_sum += inner[x * sector_size + y];
      maximum_inner_flux =
          fmax(maximum_inner_flux, fmax(fabs(forward), fabs(reverse)));
      *inner_db_residual =
          fmax(*inner_db_residual, fabs(forward - reverse));
    }
    *inner_row_residual =
        fmax(*inner_row_residual, fabs(row_sum - 1.0));
  }
  if (maximum_inner_flux <= 0.0 || total_target <= 0.0 ||
      total_surrogate <= 0.0) {
    goto fail;
  }
  *inner_db_residual /= maximum_inner_flux;
  {
    size_t step;
    for (step = 0; step < surrogate_policy->step_count; ++step) {
      double *swap;
      pilot_dense_transition_multiply(power, inner, sector_size, temporary);
      swap = power;
      power = temporary;
      temporary = swap;
    }
  }
  for (x = 0; x < sector_size; ++x) {
    size_t y;
    double row_sum = 0.0;
    double outer_off_diagonal = 0.0;
    for (y = 0; y < sector_size; ++y) {
      const double forward =
          surrogate_weights[x] * power[x * sector_size + y];
      const double reverse =
          surrogate_weights[y] * power[y * sector_size + x];
      row_sum += power[x * sector_size + y];
      maximum_powered_flux =
          fmax(maximum_powered_flux, fmax(fabs(forward), fabs(reverse)));
      *powered_db_residual =
          fmax(*powered_db_residual, fabs(forward - reverse));
      if (x != y && power[x * sector_size + y] > 0.0 &&
          power[y * sector_size + x] > 0.0) {
        const double observed =
            log(power[y * sector_size + x]) -
            log(power[x * sector_size + y]);
        const double expected =
            log(surrogate_weights[x]) - log(surrogate_weights[y]);
        *powered_ratio_residual =
            fmax(*powered_ratio_residual, fabs(observed - expected));
      }
      if (x != y) {
        const double acceptance =
            fmin(1.0, (target_weights[y] * surrogate_weights[x]) /
                           (target_weights[x] * surrogate_weights[y]));
        outer[x * sector_size + y] =
            power[x * sector_size + y] * acceptance;
        outer_off_diagonal += outer[x * sector_size + y];
      }
    }
    *powered_row_residual =
        fmax(*powered_row_residual, fabs(row_sum - 1.0));
    outer[x * sector_size + x] = 1.0 - outer_off_diagonal;
  }
  if (maximum_powered_flux <= 0.0) goto fail;
  *powered_db_residual /= maximum_powered_flux;
  for (x = 0; x < sector_size; ++x) {
    const double pi_x = target_weights[x] / total_target;
    double stationary = 0.0;
    size_t y;
    maximum_target_probability = fmax(maximum_target_probability, pi_x);
    for (y = 0; y < sector_size; ++y) {
      const double pi_y = target_weights[y] / total_target;
      const double forward = pi_x * outer[x * sector_size + y];
      const double reverse = pi_y * outer[y * sector_size + x];
      maximum_outer_flux =
          fmax(maximum_outer_flux, fmax(fabs(forward), fabs(reverse)));
      *outer_db_residual =
          fmax(*outer_db_residual, fabs(forward - reverse));
      stationary += pi_y * outer[y * sector_size + x];
    }
    *outer_stationary_residual =
        fmax(*outer_stationary_residual, fabs(stationary - pi_x));
  }
  if (maximum_outer_flux <= 0.0 || maximum_target_probability <= 0.0) {
    goto fail;
  }
  *outer_db_residual /= maximum_outer_flux;
  *outer_stationary_residual /= maximum_target_probability;
  free(outer);
  free(temporary);
  free(power);
  free(inner);
  return isfinite(*inner_row_residual) &&
         isfinite(*inner_db_residual) &&
         isfinite(*powered_row_residual) &&
         isfinite(*powered_db_residual) &&
         isfinite(*powered_ratio_residual) &&
         isfinite(*outer_db_residual) &&
         isfinite(*outer_stationary_residual) &&
         *inner_row_residual <= 256.0 * DBL_EPSILON &&
         *inner_db_residual <= 256.0 * DBL_EPSILON &&
         *powered_row_residual <= 256.0 * DBL_EPSILON &&
         *powered_db_residual <= 256.0 * DBL_EPSILON &&
         *powered_ratio_residual <= 1024.0 * DBL_EPSILON &&
         *outer_db_residual <= 256.0 * DBL_EPSILON &&
         *outer_stationary_residual <= 256.0 * DBL_EPSILON;

fail:
  free(outer);
  free(temporary);
  free(power);
  free(inner);
  return 0;
}

static int pilot_exact_l4_partial_independence_balance(
    int site_count, const double *target_weights,
    const PilotSurrogateTable *partial_table, size_t sector_size,
    double *row_residual, double *db_residual,
    double *stationary_residual) {
  double *transition = NULL;
  double total_target = 0.0;
  double maximum_flux = 0.0;
  double maximum_probability = 0.0;
  size_t x;
  if (target_weights == NULL || partial_table == NULL ||
      row_residual == NULL || db_residual == NULL ||
      stationary_residual == NULL) {
    return 0;
  }
  *row_residual = 0.0;
  *db_residual = 0.0;
  *stationary_residual = 0.0;
  if (site_count != 4) return 1;
  if (sector_size == 0 || partial_table->weights == NULL ||
      !isfinite(partial_table->total_weight) ||
      partial_table->total_weight <= 0.0 ||
      sector_size > SIZE_MAX / sector_size ||
      sector_size * sector_size > SIZE_MAX / sizeof(*transition)) {
    return 0;
  }
  transition =
      (double *)calloc(sector_size * sector_size, sizeof(*transition));
  if (transition == NULL) return 0;
  for (x = 0; x < sector_size; ++x) total_target += target_weights[x];
  if (!isfinite(total_target) || total_target <= 0.0) goto fail_partial;
  for (x = 0; x < sector_size; ++x) {
    double accepted_sum = 0.0;
    size_t y;
    for (y = 0; y < sector_size; ++y) {
      const double q_y =
          partial_table->weights[y] / partial_table->total_weight;
      const double acceptance =
          fmin(1.0, (target_weights[y] * partial_table->weights[x]) /
                         (target_weights[x] * partial_table->weights[y]));
      transition[x * sector_size + y] = q_y * acceptance;
      accepted_sum += transition[x * sector_size + y];
    }
    if (!isfinite(accepted_sum) || accepted_sum <= 0.0 ||
        accepted_sum > 1.0 + 256.0 * DBL_EPSILON) {
      goto fail_partial;
    }
    transition[x * sector_size + x] += 1.0 - accepted_sum;
  }
  for (x = 0; x < sector_size; ++x) {
    const double pi_x = target_weights[x] / total_target;
    double row_sum = 0.0;
    double stationary = 0.0;
    size_t y;
    maximum_probability = fmax(maximum_probability, pi_x);
    for (y = 0; y < sector_size; ++y) {
      const double pi_y = target_weights[y] / total_target;
      const double forward = pi_x * transition[x * sector_size + y];
      const double reverse = pi_y * transition[y * sector_size + x];
      row_sum += transition[x * sector_size + y];
      stationary += pi_y * transition[y * sector_size + x];
      maximum_flux =
          fmax(maximum_flux, fmax(fabs(forward), fabs(reverse)));
      *db_residual = fmax(*db_residual, fabs(forward - reverse));
    }
    *row_residual = fmax(*row_residual, fabs(row_sum - 1.0));
    *stationary_residual =
        fmax(*stationary_residual, fabs(stationary - pi_x));
  }
  if (maximum_flux <= 0.0 || maximum_probability <= 0.0) {
    goto fail_partial;
  }
  *db_residual /= maximum_flux;
  *stationary_residual /= maximum_probability;
  free(transition);
  return isfinite(*row_residual) && isfinite(*db_residual) &&
         isfinite(*stationary_residual) &&
         *row_residual <= 256.0 * DBL_EPSILON &&
         *db_residual <= 256.0 * DBL_EPSILON &&
         *stationary_residual <= 256.0 * DBL_EPSILON;

fail_partial:
  free(transition);
  return 0;
}

static int pilot_restart_check(
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const PilotCandidateTable *table, const uint64_t *configurations,
    size_t sector_size, const size_t *lookup, size_t lookup_count,
    uint64_t seed) {
  PilotChain full;
  PilotChain split;
  if (!pilot_chain_initialize(seed, table->guide_weights, sector_size,
                              configurations, &full) ||
      !pilot_chain_initialize(seed, table->guide_weights, sector_size,
                              configurations, &split) ||
      !pilot_chain_advance(model, proposal_policy, table, lookup,
                           lookup_count, 257, NULL, &full) ||
      !pilot_chain_advance(model, proposal_policy, table, lookup,
                           lookup_count, 101, NULL, &split) ||
      !pilot_chain_advance(model, proposal_policy, table, lookup,
                           lookup_count, 156, NULL, &split)) {
    return 0;
  }
  return full.current_index == split.current_index &&
         full.current_words == split.current_words &&
         full.rng.state == split.rng.state &&
         full.rng.draws == split.rng.draws &&
         full.trace_hash == split.trace_hash &&
         full.accepted == split.accepted &&
         full.rejected == split.rejected &&
         full.neighbor_attempted == split.neighbor_attempted &&
         full.global_attempted == split.global_attempted &&
         full.global_self == split.global_self &&
         full.shell_attempted == split.shell_attempted;
}

static int pilot_surrogate_restart_check(
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy,
    const PilotCandidateTable *target_table,
    const PilotSurrogateTable *surrogate_table,
    const uint64_t *configurations, size_t sector_size,
    const size_t *lookup, size_t lookup_count, uint64_t seed) {
  PilotChain full;
  PilotChain split;
  if (!pilot_chain_initialize(seed, target_table->guide_weights, sector_size,
                              configurations, &full) ||
      !pilot_chain_initialize(seed, target_table->guide_weights, sector_size,
                              configurations, &split) ||
      !pilot_chain_advance_surrogate(
          model, surrogate_policy, target_table, surrogate_table, lookup,
          lookup_count, 257, NULL, &full) ||
      !pilot_chain_advance_surrogate(
          model, surrogate_policy, target_table, surrogate_table, lookup,
          lookup_count, 101, NULL, &split) ||
      !pilot_chain_advance_surrogate(
          model, surrogate_policy, target_table, surrogate_table, lookup,
          lookup_count, 156, NULL, &split)) {
    return 0;
  }
  return full.current_index == split.current_index &&
         full.current_words == split.current_words &&
         full.rng.state == split.rng.state &&
         full.rng.draws == split.rng.draws &&
         full.trace_hash == split.trace_hash &&
         full.accepted == split.accepted &&
         full.rejected == split.rejected &&
         full.surrogate_inner_attempted ==
             split.surrogate_inner_attempted &&
         full.surrogate_inner_accepted == split.surrogate_inner_accepted &&
         full.surrogate_inner_rejected == split.surrogate_inner_rejected &&
         full.surrogate_final_changed == split.surrogate_final_changed &&
         full.surrogate_final_self == split.surrogate_final_self;
}

static int pilot_surrogate_chains_equal(
    const PilotChain *left, const PilotChain *right) {
  return left != NULL && right != NULL &&
         left->current_index == right->current_index &&
         left->current_words == right->current_words &&
         left->rng.state == right->rng.state &&
         left->rng.draws == right->rng.draws &&
         left->trace_hash == right->trace_hash &&
         left->accepted == right->accepted &&
         left->rejected == right->rejected &&
         left->surrogate_inner_attempted ==
             right->surrogate_inner_attempted &&
         left->surrogate_inner_accepted == right->surrogate_inner_accepted &&
         left->surrogate_inner_rejected == right->surrogate_inner_rejected &&
         left->surrogate_final_changed == right->surrogate_final_changed &&
         left->surrogate_final_self == right->surrogate_final_self;
}

static int pilot_partial_callback_restart_check(
    const MVMCKrylovFockModel *model,
    const MVMCKrylovPositiveSamplerSurrogatePolicy *surrogate_policy,
    const PilotCandidateTable *target_table,
    PilotPartialCallbackContext *callback,
    const uint64_t *configurations, size_t sector_size,
    const size_t *lookup, size_t lookup_count, uint64_t seed) {
  PilotChain full;
  PilotChain split;
  if (!pilot_chain_initialize(seed, target_table->guide_weights, sector_size,
                              configurations, &full) ||
      !pilot_chain_initialize(seed, target_table->guide_weights, sector_size,
                              configurations, &split) ||
      !pilot_chain_advance_surrogate_callback(
          model, surrogate_policy, target_table, callback, lookup,
          lookup_count, 257, NULL, &full) ||
      !pilot_chain_advance_surrogate_callback(
          model, surrogate_policy, target_table, callback, lookup,
          lookup_count, 101, NULL, &split) ||
      !pilot_chain_advance_surrogate_callback(
          model, surrogate_policy, target_table, callback, lookup,
          lookup_count, 156, NULL, &split)) {
    return 0;
  }
  return pilot_surrogate_chains_equal(&full, &split) &&
         isfinite(full.surrogate_log_weight) &&
         isfinite(split.surrogate_log_weight) &&
         full.surrogate_log_weight == split.surrogate_log_weight;
}

static int pilot_partial_restart_check(
    const PilotCandidateTable *target_table,
    const PilotSurrogateTable *partial_table,
    const uint64_t *configurations, size_t sector_size,
    uint64_t partial_policy_hash, uint64_t seed) {
  PilotChain full;
  PilotChain split;
  if (!pilot_chain_initialize(seed, target_table->guide_weights, sector_size,
                              configurations, &full) ||
      !pilot_chain_initialize(seed, target_table->guide_weights, sector_size,
                              configurations, &split) ||
      !pilot_chain_advance_partial_independence(
          target_table, partial_table, configurations, sector_size,
          partial_policy_hash, 257, NULL, &full) ||
      !pilot_chain_advance_partial_independence(
          target_table, partial_table, configurations, sector_size,
          partial_policy_hash, 101, NULL, &split) ||
      !pilot_chain_advance_partial_independence(
          target_table, partial_table, configurations, sector_size,
          partial_policy_hash, 156, NULL, &split)) {
    return 0;
  }
  return full.current_index == split.current_index &&
         full.current_words == split.current_words &&
         full.rng.state == split.rng.state &&
         full.rng.draws == split.rng.draws &&
         full.trace_hash == split.trace_hash &&
         full.accepted == split.accepted &&
         full.rejected == split.rejected &&
         full.global_attempted == split.global_attempted &&
         full.global_self == split.global_self;
}

static int pilot_emit_flat_control(
    int control_index, const char *control_id, int uniform_guide,
    int site_count, int qp_total, size_t cache_bytes, double rho,
    const MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT],
    const double complex (*raw_values)[PROFILE_DEPTH_COUNT],
    const double norms[PROFILE_DEPTH_COUNT], size_t sector_size) {
  MVMCKrylovPositiveGuidePolicy guide_policy;
  MVMCKrylovMatrixMeasurementPolicy measurement_policy;
  PilotCandidateTable table;
  double log_basis_scale[PROFILE_DEPTH_COUNT];
  double eta = NAN;
  double maximum_iid_ratio = NAN;
  double required_tau = NAN;
  MarkovMatrixFamily limiting_family = MARKOV_FAMILY_S;
  size_t limiting_entry = 0;
  int order;

  memset(&table, 0, sizeof(table));
  if (control_id == NULL ||
      !markov_compute_scales_and_eta(sector_size, raw_values, norms,
                                     uniform_guide ? 1.0 : rho,
                                     log_basis_scale, &eta)) {
    return 0;
  }
  if (uniform_guide) {
    memset(&guide_policy, 0, sizeof(guide_policy));
    memset(&measurement_policy, 0, sizeof(measurement_policy));
    eta = 1.0;
    guide_policy.order = MARKOV_ORDER;
    guide_policy.eta = eta;
    guide_policy.policy_hash = UINT64_C(0x50345335554e4946);
    measurement_policy.order = MARKOV_ORDER;
    measurement_policy.eta = eta;
    for (order = 0; order <= MARKOV_ORDER; ++order) {
      guide_policy.lambda[order] = 0.0;
      guide_policy.log_basis_scale[order] = log_basis_scale[order];
      measurement_policy.guide_lambda[order] = 0.0;
      measurement_policy.target_weight[order] = 1.0;
      measurement_policy.log_basis_scale[order] = log_basis_scale[order];
    }
  } else if (!markov_init_guide_policy(
                 site_count, qp_total, cache_bytes, rho, eta,
                 log_basis_scale, &guide_policy) ||
             !markov_init_measurement_policy(
                 eta, log_basis_scale, &measurement_policy)) {
    return 0;
  }
  if (!pilot_table_create(&measurement_policy, exact_values, sector_size,
                          &table) ||
      !pilot_iid_summary(&table, &maximum_iid_ratio, &required_tau,
                         &limiting_family, &limiting_entry)) {
    pilot_table_destroy(&table);
    return 0;
  }
  if (world_rank() == 0) {
    printf("PILOT_CONTROL schema=%d control_index=%d control_id=%s",
           PILOT_FLAT_SCHEMA_VERSION, control_index, control_id);
    printf(" promotable=0 uniform_guide=%d rho=%.17g eta=%.17g",
           uniform_guide, uniform_guide ? 0.0 : rho, eta);
    printf(" guide_policy_hash=%" PRIu64 " guide_shape_hash=%" PRIu64,
           guide_policy.policy_hash, pilot_guide_shape_hash(&guide_policy));
    printf(" table_hash=%" PRIu64 " anchor_count=%zu",
           table.table_hash, sector_size);
    printf(" guide_minimum=%.17g guide_maximum=%.17g guide_mean=%.17g",
           table.minimum_guide, table.maximum_guide, table.mean_guide);
    printf(" guide_maximum_minimum_ratio=%.17g",
           table.guide_maximum_minimum_ratio);
    printf(" uniform_sector_self_probability=%.17g",
           table.uniform_sector_self_probability);
    printf(" importance_weight_maximum_mean_ratio=%.17g",
           table.importance_weight_maximum_mean_ratio);
    printf(" importance_weight_ess_fraction=%.17g",
           table.importance_weight_ess_fraction);
    printf(" maximum_iid_se_budget_ratio=%.17g required_tau_for_0p90=%.17g",
           maximum_iid_ratio, required_tau);
    printf(" limiting_family=%s limiting_entry=%zu\n",
           markov_family_name(limiting_family), limiting_entry);
  }
  pilot_table_destroy(&table);
  return 1;
}

static int pilot_emit_surrogate_control(
    int control_index, size_t floor_multiplier,
    int site_count, int qp_total, size_t cache_bytes,
    const MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT],
    const double complex (*raw_values)[PROFILE_DEPTH_COUNT],
    const double norms[PROFILE_DEPTH_COUNT], size_t sector_size) {
  MVMCKrylovPositiveGuidePolicy guide_policy;
  MVMCKrylovMatrixMeasurementPolicy measurement_policy;
  MVMCKrylovPositiveSamplerSurrogatePolicy surrogate_policy;
  PilotCandidateTable target_table;
  PilotSurrogateTable surrogate_table;
  double log_basis_scale[PROFILE_DEPTH_COUNT];
  double eta = NAN;
  double maximum_iid_ratio = NAN;
  double required_tau = NAN;
  MarkovMatrixFamily limiting_family = MARKOV_FAMILY_S;
  size_t limiting_entry = 0;
  memset(&target_table, 0, sizeof(target_table));
  memset(&surrogate_table, 0, sizeof(surrogate_table));
  if (!markov_compute_scales_and_eta(sector_size, raw_values, norms, 1.0e-2,
                                     log_basis_scale, &eta) ||
      !markov_init_guide_policy(site_count, qp_total, cache_bytes, 1.0e-2,
                                eta, log_basis_scale, &guide_policy) ||
      !markov_init_measurement_policy(eta, log_basis_scale,
                                      &measurement_policy) ||
      mvmc_krylov_positive_sampler_surrogate_policy_create(
          1, floor_multiplier, &guide_policy, &surrogate_policy) !=
          MVMC_KRYLOV_STATUS_OK ||
      !pilot_table_create(&measurement_policy, exact_values, sector_size,
                          &target_table) ||
      !pilot_iid_summary(&target_table, &maximum_iid_ratio, &required_tau,
                         &limiting_family, &limiting_entry) ||
      !pilot_surrogate_table_create(
          &guide_policy, &surrogate_policy, exact_values, &target_table,
          sector_size, 1, &surrogate_table)) {
    pilot_surrogate_table_destroy(&surrogate_table);
    pilot_table_destroy(&target_table);
    return 0;
  }
  if (world_rank() == 0) {
    printf("PILOT_SURROGATE_CONTROL schema=%d control_index=%d",
           PILOT_SURROGATE_SCHEMA_VERSION, control_index);
    printf(" control_id=exact_q_independence promotable=0 alpha=%zu",
           floor_multiplier);
    printf(" rho=0.01 eta=%.17g guide_policy_hash=%" PRIu64,
           eta, guide_policy.policy_hash);
    printf(" guide_shape_hash=%" PRIu64,
           pilot_guide_shape_hash(&guide_policy));
    printf(" surrogate_policy_hash=%" PRIu64,
           surrogate_policy.policy_hash);
    printf(" target_table_hash=%" PRIu64,
           target_table.table_hash);
    printf(" surrogate_table_hash=%" PRIu64 " anchor_count=%zu",
           surrogate_table.table_hash, sector_size);
    printf(" surrogate_minimum=%.17g surrogate_maximum=%.17g",
           surrogate_table.minimum_weight, surrogate_table.maximum_weight);
    printf(" surrogate_mean=%.17g target_ratio_minimum=%.17g",
           surrogate_table.mean_weight,
           surrogate_table.minimum_target_ratio);
    printf(" target_ratio_maximum=%.17g",
           surrogate_table.maximum_target_ratio);
    printf(" target_ratio_maximum_mean_ratio=%.17g",
           surrogate_table.target_ratio_maximum_mean_ratio);
    printf(" target_ratio_ess_fraction=%.17g",
           surrogate_table.target_ratio_ess_fraction);
    printf(" ideal_acceptance=%.17g ideal_configuration_change=%.17g",
           surrogate_table.ideal_independence_acceptance,
           surrogate_table.ideal_independence_configuration_change);
    printf(" maximum_iid_se_budget_ratio=%.17g required_tau_for_0p90=%.17g",
           maximum_iid_ratio, required_tau);
    printf(" limiting_family=%s limiting_entry=%zu\n",
           markov_family_name(limiting_family), limiting_entry);
  }
  pilot_surrogate_table_destroy(&surrogate_table);
  pilot_table_destroy(&target_table);
  return 1;
}

static int pilot_run_candidate(
    int candidate_index, int site_count, int qp_total, int sample_count,
    size_t cache_bytes, double rho,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const uint64_t seeds[PILOT_SEED_COUNT], int pilot_schema_version,
    size_t predicted_sample_count,
    double maximum_predicted_ratio_threshold,
    double maximum_tau_threshold,
    const MVMCKrylovFockModel *model,
    const MVMCKrylovBoundedLimits *limits,
    const MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT],
    const double complex (*raw_values)[PROFILE_DEPTH_COUNT],
    const double norms[PROFILE_DEPTH_COUNT], const uint64_t *configurations,
    size_t sector_size, const size_t *lookup, size_t lookup_count) {
  MVMCKrylovPositiveGuidePolicy guide_policy;
  MVMCKrylovMatrixMeasurementPolicy measurement_policy;
  PilotCandidateTable table;
  double log_basis_scale[PROFILE_DEPTH_COUNT];
  double eta = NAN;
  double maximum_predicted_ratio = 0.0;
  double maximum_tau = 0.0;
  double maximum_iid_ratio = NAN;
  double required_tau = NAN;
  double db_residual = NAN;
  double stationary_residual = NAN;
  double proposal_row_residual = NAN;
  double proposal_symmetry_residual = NAN;
  size_t shell_max_distance = 0;
  size_t shell_distance = 0;
  size_t shell_count = 0;
  uint64_t proposal_model_hash = 0;
  MarkovMatrixFamily iid_limiting_family = MARKOV_FAMILY_S;
  size_t iid_limiting_entry = 0;
  int table_pass;
  int balance_pass;
  int restart_pass;
  int shell_cardinality_pass = 1;
  int seed_index;

  memset(&table, 0, sizeof(table));
  if (proposal_policy == NULL || seeds == NULL ||
      predicted_sample_count == 0 ||
      !isfinite(maximum_predicted_ratio_threshold) ||
      maximum_predicted_ratio_threshold <= 0.0 ||
      !isfinite(maximum_tau_threshold) || maximum_tau_threshold <= 0.0 ||
      !markov_compute_scales_and_eta(sector_size, raw_values, norms, rho,
                                     log_basis_scale, &eta) ||
      !markov_init_guide_policy(site_count, qp_total, cache_bytes, rho, eta,
                                log_basis_scale, &guide_policy) ||
      !markov_init_measurement_policy(eta, log_basis_scale,
                                      &measurement_policy) ||
      mvmc_krylov_positive_sampler_proposal_model_policy_hash(
          proposal_policy, model, &configurations[0], 1,
          &proposal_model_hash) !=
          MVMC_KRYLOV_STATUS_OK) {
    return 0;
  }
  if (proposal_policy->kernel_id ==
      MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_SHELL_MIXTURE) {
    if (mvmc_krylov_fock_proposal_resolve_shell_distance(
            model, proposal_policy->distance_numerator,
            proposal_policy->distance_denominator, &shell_max_distance,
            &shell_distance) != MVMC_KRYLOV_STATUS_OK ||
        mvmc_krylov_fock_proposal_count_shell(
            model, shell_distance, &shell_count) !=
            MVMC_KRYLOV_STATUS_OK) {
      return 0;
    }
  }
  table_pass = pilot_table_create(&measurement_policy, exact_values,
                                  sector_size, &table);
  table_pass = table_pass &&
               pilot_iid_summary(&table, &maximum_iid_ratio, &required_tau,
                                 &iid_limiting_family,
                                 &iid_limiting_entry);
  balance_pass = table_pass && pilot_exact_l4_balance(
                                   model, proposal_policy, configurations,
                                   table.guide_weights, sector_size, lookup,
                                   lookup_count, &proposal_row_residual,
                                   &proposal_symmetry_residual, &db_residual,
                                   &stationary_residual,
                                   &shell_cardinality_pass);
  restart_pass = table_pass && pilot_restart_check(
                                  model, proposal_policy, &table,
                                  configurations, sector_size, lookup,
                                  lookup_count, seeds[0]);
  if (!table_pass || !balance_pass || !restart_pass) {
    pilot_table_destroy(&table);
    return 0;
  }

  if (world_rank() == 0) {
    printf("PILOT_CANDIDATE schema=%d candidate_index=%d",
           pilot_schema_version, candidate_index);
    printf(" site_count=%d qp_total=%d sample_count_per_seed=%d",
           site_count, qp_total, sample_count);
    printf(" cache_requested_bytes=%zu rho=%.17g eta=%.17g",
           cache_bytes, rho, eta);
    printf(" predicted_official_sample_count=%zu",
           predicted_sample_count);
    printf(" maximum_predicted_ratio_threshold=%.17g",
           maximum_predicted_ratio_threshold);
    printf(" maximum_tau_threshold=%.17g", maximum_tau_threshold);
    printf(" proposal_kernel=%" PRIu64, proposal_policy->kernel_id);
    printf(" global_numerator=%zu global_denominator=%zu",
           proposal_policy->global_numerator,
           proposal_policy->global_denominator);
    printf(" neighbor_numerator=%zu neighbor_denominator=%zu",
           proposal_policy->neighbor_numerator,
           proposal_policy->neighbor_denominator);
    printf(" distance_numerator=%zu distance_denominator=%zu",
           proposal_policy->distance_numerator,
           proposal_policy->distance_denominator);
    printf(" distance_rounding_rule=%" PRIu64,
           proposal_policy->distance_rounding_rule);
    printf(" shell_max_distance=%zu shell_distance=%zu shell_count=%zu",
           shell_max_distance, shell_distance, shell_count);
    printf(" proposal_policy_hash=%" PRIu64,
           proposal_policy->policy_hash);
    printf(" proposal_model_hash=%" PRIu64,
           proposal_model_hash);
    printf(" guide_policy_hash=%" PRIu64 " guide_shape_hash=%" PRIu64,
           guide_policy.policy_hash, pilot_guide_shape_hash(&guide_policy));
    printf(" table_hash=%" PRIu64 " anchor_count=%zu",
           table.table_hash, sector_size);
    printf(" guide_minimum=%.17g guide_maximum=%.17g guide_mean=%.17g",
           table.minimum_guide, table.maximum_guide, table.mean_guide);
    printf(" guide_maximum_minimum_ratio=%.17g",
           table.guide_maximum_minimum_ratio);
    printf(" uniform_sector_self_probability=%.17g",
           table.uniform_sector_self_probability);
    printf(" importance_weight_maximum_mean_ratio=%.17g",
           table.importance_weight_maximum_mean_ratio);
    printf(" importance_weight_ess_fraction=%.17g",
           table.importance_weight_ess_fraction);
    printf(" maximum_iid_se_budget_ratio=%.17g required_tau_for_0p90=%.17g",
           maximum_iid_ratio, required_tau);
    printf(" iid_limiting_family=%s iid_limiting_entry=%zu",
           markov_family_name(iid_limiting_family), iid_limiting_entry);
    printf(" exact_balance_pass=%d db_residual=%.17g",
           balance_pass, db_residual);
    printf(" stationary_residual=%.17g", stationary_residual);
    printf(" proposal_row_residual=%.17g", proposal_row_residual);
    printf(" proposal_symmetry_residual=%.17g",
           proposal_symmetry_residual);
    printf(" shell_cardinality_pass=%d restart_pass=%d\n",
           shell_cardinality_pass, restart_pass);
  }

  for (seed_index = 0; seed_index < PILOT_SEED_COUNT; ++seed_index) {
    PilotChain chain;
    double *series = NULL;
    double seed_maximum_tau = 0.0;
    double seed_maximum_ratio = 0.0;
    MarkovMatrixFamily family;
    if ((size_t)sample_count >
        SIZE_MAX / (PILOT_SERIES_COUNT * sizeof(*series))) {
      pilot_table_destroy(&table);
      return 0;
    }
    series = (double *)calloc(
        (size_t)sample_count * PILOT_SERIES_COUNT, sizeof(*series));
    if (series == NULL ||
        !pilot_chain_initialize(seeds[seed_index], table.guide_weights,
                                sector_size, configurations, &chain) ||
        !pilot_chain_advance(model, proposal_policy, &table,
                             lookup, lookup_count,
                             (size_t)sample_count, series, &chain)) {
      free(series);
      pilot_table_destroy(&table);
      return 0;
    }
    for (family = MARKOV_FAMILY_S; family <= MARKOV_FAMILY_B; ++family) {
      size_t entry;
      for (entry = 0; entry < 6; ++entry) {
        double component_tau[2];
        int component;
        double predicted_se;
        double predicted_ratio;
        for (component = 0; component < 2; ++component) {
          const size_t index = pilot_series_index(family, entry, component);
          MVMCKrylovTauIntResult tau;
          if (mvmc_krylov_tau_int_geyer_initial_positive(
                  &series[index * (size_t)sample_count],
                  (size_t)sample_count, &tau) != MVMC_KRYLOV_STATUS_OK ||
              !tau.valid) {
            free(series);
            pilot_table_destroy(&table);
            return 0;
          }
          component_tau[component] = tau.tau_int;
          seed_maximum_tau = fmax(seed_maximum_tau, tau.tau_int);
          if (world_rank() == 0) {
            printf("PILOT_SERIES candidate_index=%d seed_index=%d",
                   candidate_index, seed_index);
            printf(" family=%s entry=%zu component=%s",
                   markov_family_name(family), entry,
                   component == 0 ? "real" : "imag");
            printf(" iid_variance=%.17g tau_int=%.17g\n",
                   table.iid_variance[index], tau.tau_int);
          }
        }
        predicted_se = sqrt(
            (table.iid_variance[pilot_series_index(family, entry, 0)] *
                 component_tau[0] +
             table.iid_variance[pilot_series_index(family, entry, 1)] *
                 component_tau[1]) /
            ((double)predicted_sample_count * table.mean_denominator *
             table.mean_denominator));
        predicted_ratio =
            predicted_se / table.budget.entry[family][entry].budget;
        seed_maximum_ratio = fmax(seed_maximum_ratio, predicted_ratio);
        if (world_rank() == 0) {
          printf("PILOT_ENTRY candidate_index=%d seed_index=%d",
                 candidate_index, seed_index);
          printf(" family=%s entry=%zu predicted_se=%.17g",
                 markov_family_name(family), entry, predicted_se);
          printf(" budget=%.17g predicted_se_budget_ratio=%.17g\n",
                 table.budget.entry[family][entry].budget,
                 predicted_ratio);
        }
      }
    }
    maximum_tau = fmax(maximum_tau, seed_maximum_tau);
    maximum_predicted_ratio =
        fmax(maximum_predicted_ratio, seed_maximum_ratio);
    if (world_rank() == 0) {
      printf("PILOT_SEED candidate_index=%d seed_index=%d",
             candidate_index, seed_index);
      printf(" seed=%" PRIu64 " seed_hex=0x%016" PRIx64,
             seeds[seed_index], seeds[seed_index]);
      printf(" maximum_tau_int=%.17g maximum_predicted_se_budget_ratio=%.17g",
             seed_maximum_tau, seed_maximum_ratio);
      printf(" accepted=%" PRIu64 " rejected=%" PRIu64,
             chain.accepted, chain.rejected);
      printf(" neighbor_attempted=%" PRIu64 " global_attempted=%" PRIu64,
             chain.neighbor_attempted, chain.global_attempted);
      printf(" global_self=%" PRIu64 " shell_attempted=%" PRIu64,
             chain.global_self, chain.shell_attempted);
      printf(" trace_hash=%" PRIu64 "\n", chain.trace_hash);
    }
    free(series);
  }
  if (world_rank() == 0) {
    printf("PILOT_DECISION candidate_index=%d", candidate_index);
    printf(" physical_case_eligible=%d",
           maximum_predicted_ratio <= maximum_predicted_ratio_threshold &&
               maximum_tau <= maximum_tau_threshold && balance_pass &&
               restart_pass);
    printf(" maximum_predicted_se_budget_ratio=%.17g",
           maximum_predicted_ratio);
    printf(" maximum_tau_int=%.17g table_pass=%d balance_pass=%d",
           maximum_tau, table_pass, balance_pass);
    printf(" restart_pass=%d resource_limits_hash=%" PRIu64 "\n",
           restart_pass, limits->amplitude_policy_hash);
  }
  pilot_table_destroy(&table);
  return 1;
}

static int pilot_run_surrogate_candidate(
    int candidate_index, size_t surrogate_step_count,
    size_t floor_multiplier, int site_count, int qp_total,
    int sample_count, size_t cache_bytes,
    const uint64_t seeds[PILOT_SEED_COUNT],
    const MVMCKrylovFockModel *model,
    const MVMCKrylovBoundedLimits *limits,
    const MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT],
    const double complex (*raw_values)[PROFILE_DEPTH_COUNT],
    const double norms[PROFILE_DEPTH_COUNT],
    const uint64_t *configurations, size_t sector_size,
    const size_t *lookup, size_t lookup_count) {
  MVMCKrylovPositiveGuidePolicy guide_policy;
  MVMCKrylovMatrixMeasurementPolicy measurement_policy;
  MVMCKrylovPositiveSamplerProposalPolicy neighbor_policy;
  MVMCKrylovPositiveSamplerSurrogatePolicy surrogate_policy;
  PilotCandidateTable target_table;
  PilotSurrogateTable surrogate_table;
  double log_basis_scale[PROFILE_DEPTH_COUNT];
  double eta = NAN;
  double maximum_predicted_ratio = 0.0;
  double maximum_tau = 0.0;
  double maximum_iid_ratio = NAN;
  double required_tau = NAN;
  double inner_row_residual = NAN;
  double inner_db_residual = NAN;
  double powered_row_residual = NAN;
  double powered_db_residual = NAN;
  double powered_ratio_residual = NAN;
  double outer_db_residual = NAN;
  double outer_stationary_residual = NAN;
  MarkovMatrixFamily iid_limiting_family = MARKOV_FAMILY_S;
  size_t iid_limiting_entry = 0;
  uint64_t proposal_model_hash = 0;
  int table_pass;
  int balance_pass;
  int restart_pass;
  int seed_index;

  memset(&target_table, 0, sizeof(target_table));
  memset(&surrogate_table, 0, sizeof(surrogate_table));
  if (!markov_compute_scales_and_eta(sector_size, raw_values, norms, 1.0e-2,
                                     log_basis_scale, &eta) ||
      !markov_init_guide_policy(site_count, qp_total, cache_bytes, 1.0e-2,
                                eta, log_basis_scale, &guide_policy) ||
      !markov_init_measurement_policy(eta, log_basis_scale,
                                      &measurement_policy) ||
      mvmc_krylov_positive_sampler_proposal_policy_create(
          0, 1, &neighbor_policy) != MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_surrogate_policy_create(
          surrogate_step_count, floor_multiplier, &guide_policy,
          &surrogate_policy) != MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_proposal_model_policy_hash(
          &neighbor_policy, model, &configurations[0], 1,
          &proposal_model_hash) != MVMC_KRYLOV_STATUS_OK) {
    return 0;
  }
  table_pass =
      pilot_table_create(&measurement_policy, exact_values, sector_size,
                         &target_table) &&
      pilot_iid_summary(&target_table, &maximum_iid_ratio, &required_tau,
                        &iid_limiting_family, &iid_limiting_entry) &&
      pilot_surrogate_table_create(
          &guide_policy, &surrogate_policy, exact_values, &target_table,
          sector_size, 0, &surrogate_table);
  balance_pass =
      table_pass && pilot_exact_l4_surrogate_balance(
                        model, &surrogate_policy, configurations,
                        target_table.guide_weights, surrogate_table.weights,
                        sector_size, &inner_row_residual, &inner_db_residual,
                        &powered_row_residual, &powered_db_residual,
                        &powered_ratio_residual, &outer_db_residual,
                        &outer_stationary_residual);
  restart_pass =
      table_pass && pilot_surrogate_restart_check(
                        model, &surrogate_policy, &target_table,
                        &surrogate_table, configurations, sector_size, lookup,
                        lookup_count, seeds[0]);
  if (!table_pass || !balance_pass || !restart_pass) {
    pilot_surrogate_table_destroy(&surrogate_table);
    pilot_table_destroy(&target_table);
    return 0;
  }
  if (world_rank() == 0) {
    printf("PILOT_CANDIDATE schema=%d candidate_index=%d",
           PILOT_SURROGATE_SCHEMA_VERSION, candidate_index);
    printf(" site_count=%d qp_total=%d sample_count_per_seed=%d",
           site_count, qp_total, sample_count);
    printf(" cache_requested_bytes=%zu rho=0.01 eta=%.17g",
           cache_bytes, eta);
    printf(" surrogate_step_count=%zu floor_multiplier=%zu",
           surrogate_step_count, floor_multiplier);
    printf(" inner_kernel=neighbor_only proposal_policy_hash=%" PRIu64,
           neighbor_policy.policy_hash);
    printf(" proposal_model_hash=%" PRIu64,
           proposal_model_hash);
    printf(" guide_policy_hash=%" PRIu64 " guide_shape_hash=%" PRIu64,
           guide_policy.policy_hash, pilot_guide_shape_hash(&guide_policy));
    printf(" surrogate_policy_hash=%" PRIu64,
           surrogate_policy.policy_hash);
    printf(" target_table_hash=%" PRIu64,
           target_table.table_hash);
    printf(" surrogate_table_hash=%" PRIu64 " anchor_count=%zu",
           surrogate_table.table_hash, sector_size);
    printf(" surrogate_minimum=%.17g surrogate_maximum=%.17g",
           surrogate_table.minimum_weight,
           surrogate_table.maximum_weight);
    printf(" surrogate_mean=%.17g target_ratio_minimum=%.17g",
           surrogate_table.mean_weight,
           surrogate_table.minimum_target_ratio);
    printf(" target_ratio_maximum=%.17g",
           surrogate_table.maximum_target_ratio);
    printf(" target_ratio_maximum_mean_ratio=%.17g",
           surrogate_table.target_ratio_maximum_mean_ratio);
    printf(" target_ratio_ess_fraction=%.17g",
           surrogate_table.target_ratio_ess_fraction);
    printf(" maximum_iid_se_budget_ratio=%.17g required_tau_for_0p90=%.17g",
           maximum_iid_ratio, required_tau);
    printf(" iid_limiting_family=%s iid_limiting_entry=%zu",
           markov_family_name(iid_limiting_family), iid_limiting_entry);
    printf(" exact_balance_pass=%d inner_row_residual=%.17g",
           balance_pass, inner_row_residual);
    printf(" inner_db_residual=%.17g powered_row_residual=%.17g",
           inner_db_residual, powered_row_residual);
    printf(" powered_db_residual=%.17g powered_ratio_residual=%.17g",
           powered_db_residual, powered_ratio_residual);
    printf(" outer_db_residual=%.17g outer_stationary_residual=%.17g",
           outer_db_residual, outer_stationary_residual);
    printf(" restart_pass=%d\n", restart_pass);
  }

  for (seed_index = 0; seed_index < PILOT_SEED_COUNT; ++seed_index) {
    PilotChain chain;
    double *series = NULL;
    double seed_maximum_tau = 0.0;
    double seed_maximum_ratio = 0.0;
    MarkovMatrixFamily family;
    if ((size_t)sample_count >
        SIZE_MAX / (PILOT_SERIES_COUNT * sizeof(*series))) {
      pilot_surrogate_table_destroy(&surrogate_table);
      pilot_table_destroy(&target_table);
      return 0;
    }
    series = (double *)calloc(
        (size_t)sample_count * PILOT_SERIES_COUNT, sizeof(*series));
    if (series == NULL ||
        !pilot_chain_initialize(seeds[seed_index], target_table.guide_weights,
                                sector_size, configurations, &chain) ||
        !pilot_chain_advance_surrogate(
            model, &surrogate_policy, &target_table, &surrogate_table, lookup,
            lookup_count, (size_t)sample_count, series, &chain)) {
      free(series);
      pilot_surrogate_table_destroy(&surrogate_table);
      pilot_table_destroy(&target_table);
      return 0;
    }
    for (family = MARKOV_FAMILY_S; family <= MARKOV_FAMILY_B; ++family) {
      size_t entry;
      for (entry = 0; entry < 6; ++entry) {
        double component_tau[2];
        double predicted_se;
        double predicted_ratio;
        int component;
        for (component = 0; component < 2; ++component) {
          const size_t index = pilot_series_index(family, entry, component);
          MVMCKrylovTauIntResult tau;
          if (mvmc_krylov_tau_int_geyer_initial_positive(
                  &series[index * (size_t)sample_count],
                  (size_t)sample_count, &tau) != MVMC_KRYLOV_STATUS_OK ||
              !tau.valid) {
            free(series);
            pilot_surrogate_table_destroy(&surrogate_table);
            pilot_table_destroy(&target_table);
            return 0;
          }
          component_tau[component] = tau.tau_int;
          seed_maximum_tau = fmax(seed_maximum_tau, tau.tau_int);
          if (world_rank() == 0) {
            printf("PILOT_SERIES candidate_index=%d seed_index=%d",
                   candidate_index, seed_index);
            printf(" family=%s entry=%zu component=%s",
                   markov_family_name(family), entry,
                   component == 0 ? "real" : "imag");
            printf(" iid_variance=%.17g tau_int=%.17g\n",
                   target_table.iid_variance[index], tau.tau_int);
          }
        }
        predicted_se = sqrt(
            (target_table.iid_variance[
                 pilot_series_index(family, entry, 0)] * component_tau[0] +
             target_table.iid_variance[
                 pilot_series_index(family, entry, 1)] * component_tau[1]) /
            (4096.0 * target_table.mean_denominator *
             target_table.mean_denominator));
        predicted_ratio =
            predicted_se / target_table.budget.entry[family][entry].budget;
        seed_maximum_ratio = fmax(seed_maximum_ratio, predicted_ratio);
        if (world_rank() == 0) {
          printf("PILOT_ENTRY candidate_index=%d seed_index=%d",
                 candidate_index, seed_index);
          printf(" family=%s entry=%zu predicted_se=%.17g",
                 markov_family_name(family), entry, predicted_se);
          printf(" budget=%.17g predicted_se_budget_ratio=%.17g\n",
                 target_table.budget.entry[family][entry].budget,
                 predicted_ratio);
        }
      }
    }
    maximum_tau = fmax(maximum_tau, seed_maximum_tau);
    maximum_predicted_ratio =
        fmax(maximum_predicted_ratio, seed_maximum_ratio);
    if (world_rank() == 0) {
      printf("PILOT_SEED candidate_index=%d seed_index=%d",
             candidate_index, seed_index);
      printf(" seed=%" PRIu64 " seed_hex=0x%016" PRIx64,
             seeds[seed_index], seeds[seed_index]);
      printf(" maximum_tau_int=%.17g maximum_predicted_se_budget_ratio=%.17g",
             seed_maximum_tau, seed_maximum_ratio);
      printf(" accepted=%" PRIu64 " rejected=%" PRIu64,
             chain.accepted, chain.rejected);
      printf(" surrogate_inner_attempted=%" PRIu64,
             chain.surrogate_inner_attempted);
      printf(" surrogate_inner_accepted=%" PRIu64,
             chain.surrogate_inner_accepted);
      printf(" surrogate_inner_rejected=%" PRIu64,
             chain.surrogate_inner_rejected);
      printf(" surrogate_final_changed=%" PRIu64,
             chain.surrogate_final_changed);
      printf(" surrogate_final_self=%" PRIu64,
             chain.surrogate_final_self);
      printf(" trace_hash=%" PRIu64 "\n", chain.trace_hash);
    }
    free(series);
  }
  if (world_rank() == 0) {
    printf("PILOT_DECISION candidate_index=%d", candidate_index);
    printf(" physical_case_eligible=%d",
           maximum_predicted_ratio <= 0.90 && maximum_tau <= 12.0 &&
               balance_pass && restart_pass);
    printf(" maximum_predicted_se_budget_ratio=%.17g",
           maximum_predicted_ratio);
    printf(" maximum_tau_int=%.17g table_pass=%d balance_pass=%d",
           maximum_tau, table_pass, balance_pass);
    printf(" restart_pass=%d resource_limits_hash=%" PRIu64 "\n",
           restart_pass, limits->amplitude_policy_hash);
  }
  pilot_surrogate_table_destroy(&surrogate_table);
  pilot_table_destroy(&target_table);
  return 1;
}

static int pilot_run_partial_candidate(
    int candidate_index, int partial_order, size_t floor_multiplier,
    int site_count, int qp_total, int sample_count, size_t cache_bytes,
    const uint64_t seeds[PILOT_SEED_COUNT],
    const MVMCKrylovBoundedLimits *limits,
    const MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT],
    const double complex (*raw_values)[PROFILE_DEPTH_COUNT],
    const double norms[PROFILE_DEPTH_COUNT],
    const uint64_t *configurations, size_t sector_size) {
  MVMCKrylovPositiveGuidePolicy guide_policy;
  MVMCKrylovMatrixMeasurementPolicy measurement_policy;
  PilotCandidateTable target_table;
  PilotSurrogateTable partial_table;
  double log_basis_scale[PROFILE_DEPTH_COUNT];
  double eta = NAN;
  double maximum_anchor_relative_residual = NAN;
  double maximum_predicted_ratio = 0.0;
  double maximum_tau = 0.0;
  double maximum_iid_ratio = NAN;
  double required_tau = NAN;
  double row_residual = NAN;
  double db_residual = NAN;
  double stationary_residual = NAN;
  MarkovMatrixFamily iid_limiting_family = MARKOV_FAMILY_S;
  size_t iid_limiting_entry = 0;
  uint64_t partial_policy_hash = 0;
  int table_pass;
  int balance_pass;
  int restart_pass;
  int seed_index;

  memset(&target_table, 0, sizeof(target_table));
  memset(&partial_table, 0, sizeof(partial_table));
  if (limits == NULL ||
      !markov_compute_scales_and_eta(sector_size, raw_values, norms, 1.0e-2,
                                     log_basis_scale, &eta) ||
      !markov_init_guide_policy(site_count, qp_total, cache_bytes, 1.0e-2,
                                eta, log_basis_scale, &guide_policy) ||
      !markov_init_measurement_policy(eta, log_basis_scale,
                                      &measurement_policy)) {
    return 0;
  }
  table_pass =
      pilot_table_create(&measurement_policy, exact_values, sector_size,
                         &target_table) &&
      pilot_iid_summary(&target_table, &maximum_iid_ratio, &required_tau,
                        &iid_limiting_family, &iid_limiting_entry) &&
      pilot_partial_table_create(
          &guide_policy, partial_order, floor_multiplier, exact_values,
          raw_values, &target_table, sector_size, &partial_table,
          &maximum_anchor_relative_residual, &partial_policy_hash);
  balance_pass =
      table_pass && pilot_exact_l4_partial_independence_balance(
                        site_count, target_table.guide_weights,
                        &partial_table, sector_size, &row_residual,
                        &db_residual, &stationary_residual);
  restart_pass =
      table_pass && pilot_partial_restart_check(
                        &target_table, &partial_table, configurations,
                        sector_size, partial_policy_hash, seeds[0]);
  if (!table_pass || !balance_pass || !restart_pass) {
    pilot_surrogate_table_destroy(&partial_table);
    pilot_table_destroy(&target_table);
    return 0;
  }
  if (world_rank() == 0) {
    printf("PILOT_PARTIAL_CANDIDATE schema=%d candidate_index=%d",
           PILOT_PARTIAL_SCHEMA_VERSION, candidate_index);
    printf(" site_count=%d qp_total=%d sample_count_per_seed=%d",
           site_count, qp_total, sample_count);
    printf(" cache_requested_bytes=%zu rho=0.01 eta=%.17g",
           cache_bytes, eta);
    printf(" partial_order=%d floor_multiplier=%zu",
           partial_order, floor_multiplier);
    printf(" proposal=exact_q_independence guide_policy_hash=%" PRIu64,
           guide_policy.policy_hash);
    printf(" guide_shape_hash=%" PRIu64,
           pilot_guide_shape_hash(&guide_policy));
    printf(" partial_policy_hash=%" PRIu64,
           partial_policy_hash);
    printf(" target_table_hash=%" PRIu64,
           target_table.table_hash);
    printf(" partial_table_hash=%" PRIu64 " anchor_count=%zu",
           partial_table.table_hash, sector_size);
    printf(" maximum_anchor_relative_residual=%.17g",
           maximum_anchor_relative_residual);
    printf(" partial_minimum=%.17g partial_maximum=%.17g",
           partial_table.minimum_weight, partial_table.maximum_weight);
    printf(" partial_mean=%.17g target_ratio_minimum=%.17g",
           partial_table.mean_weight,
           partial_table.minimum_target_ratio);
    printf(" target_ratio_maximum=%.17g",
           partial_table.maximum_target_ratio);
    printf(" target_ratio_maximum_mean_ratio=%.17g",
           partial_table.target_ratio_maximum_mean_ratio);
    printf(" target_ratio_ess_fraction=%.17g",
           partial_table.target_ratio_ess_fraction);
    printf(" ideal_acceptance=%.17g ideal_configuration_change=%.17g",
           partial_table.ideal_independence_acceptance,
           partial_table.ideal_independence_configuration_change);
    printf(" maximum_iid_se_budget_ratio=%.17g required_tau_for_0p90=%.17g",
           maximum_iid_ratio, required_tau);
    printf(" iid_limiting_family=%s iid_limiting_entry=%zu",
           markov_family_name(iid_limiting_family), iid_limiting_entry);
    printf(" exact_balance_pass=%d row_residual=%.17g",
           balance_pass, row_residual);
    printf(" db_residual=%.17g stationary_residual=%.17g",
           db_residual, stationary_residual);
    printf(" restart_pass=%d\n", restart_pass);
  }

  for (seed_index = 0; seed_index < PILOT_SEED_COUNT; ++seed_index) {
    PilotChain chain;
    double *series = NULL;
    double seed_maximum_tau = 0.0;
    double seed_maximum_ratio = 0.0;
    MarkovMatrixFamily family;
    if ((size_t)sample_count >
        SIZE_MAX / (PILOT_SERIES_COUNT * sizeof(*series))) {
      pilot_surrogate_table_destroy(&partial_table);
      pilot_table_destroy(&target_table);
      return 0;
    }
    series = (double *)calloc(
        (size_t)sample_count * PILOT_SERIES_COUNT, sizeof(*series));
    if (series == NULL ||
        !pilot_chain_initialize(seeds[seed_index], target_table.guide_weights,
                                sector_size, configurations, &chain) ||
        !pilot_chain_advance_partial_independence(
            &target_table, &partial_table, configurations, sector_size,
            partial_policy_hash, (size_t)sample_count, series, &chain)) {
      free(series);
      pilot_surrogate_table_destroy(&partial_table);
      pilot_table_destroy(&target_table);
      return 0;
    }
    for (family = MARKOV_FAMILY_S; family <= MARKOV_FAMILY_B; ++family) {
      size_t entry;
      for (entry = 0; entry < 6; ++entry) {
        double component_tau[2];
        double predicted_se;
        double predicted_ratio;
        int component;
        for (component = 0; component < 2; ++component) {
          const size_t index = pilot_series_index(family, entry, component);
          MVMCKrylovTauIntResult tau;
          if (mvmc_krylov_tau_int_geyer_initial_positive(
                  &series[index * (size_t)sample_count],
                  (size_t)sample_count, &tau) != MVMC_KRYLOV_STATUS_OK ||
              !tau.valid) {
            free(series);
            pilot_surrogate_table_destroy(&partial_table);
            pilot_table_destroy(&target_table);
            return 0;
          }
          component_tau[component] = tau.tau_int;
          seed_maximum_tau = fmax(seed_maximum_tau, tau.tau_int);
          if (world_rank() == 0) {
            printf("PILOT_SERIES candidate_index=%d seed_index=%d",
                   candidate_index, seed_index);
            printf(" family=%s entry=%zu component=%s",
                   markov_family_name(family), entry,
                   component == 0 ? "real" : "imag");
            printf(" iid_variance=%.17g tau_int=%.17g\n",
                   target_table.iid_variance[index], tau.tau_int);
          }
        }
        predicted_se = sqrt(
            (target_table.iid_variance[
                 pilot_series_index(family, entry, 0)] * component_tau[0] +
             target_table.iid_variance[
                 pilot_series_index(family, entry, 1)] * component_tau[1]) /
            (4096.0 * target_table.mean_denominator *
             target_table.mean_denominator));
        predicted_ratio =
            predicted_se / target_table.budget.entry[family][entry].budget;
        seed_maximum_ratio = fmax(seed_maximum_ratio, predicted_ratio);
        if (world_rank() == 0) {
          printf("PILOT_ENTRY candidate_index=%d seed_index=%d",
                 candidate_index, seed_index);
          printf(" family=%s entry=%zu predicted_se=%.17g",
                 markov_family_name(family), entry, predicted_se);
          printf(" budget=%.17g predicted_se_budget_ratio=%.17g\n",
                 target_table.budget.entry[family][entry].budget,
                 predicted_ratio);
        }
      }
    }
    maximum_tau = fmax(maximum_tau, seed_maximum_tau);
    maximum_predicted_ratio =
        fmax(maximum_predicted_ratio, seed_maximum_ratio);
    if (world_rank() == 0) {
      printf("PILOT_SEED candidate_index=%d seed_index=%d",
             candidate_index, seed_index);
      printf(" seed=%" PRIu64 " seed_hex=0x%016" PRIx64,
             seeds[seed_index], seeds[seed_index]);
      printf(" maximum_tau_int=%.17g maximum_predicted_se_budget_ratio=%.17g",
             seed_maximum_tau, seed_maximum_ratio);
      printf(" accepted=%" PRIu64 " rejected=%" PRIu64,
             chain.accepted, chain.rejected);
      printf(" independence_attempted=%" PRIu64,
             chain.global_attempted);
      printf(" independence_self=%" PRIu64,
             chain.global_self);
      printf(" trace_hash=%" PRIu64 "\n", chain.trace_hash);
    }
    free(series);
  }
  if (world_rank() == 0) {
    printf("PILOT_DECISION candidate_index=%d", candidate_index);
    printf(" physical_case_eligible=%d",
           maximum_predicted_ratio <= 0.90 && maximum_tau <= 12.0 &&
               table_pass && balance_pass && restart_pass);
    printf(" maximum_predicted_se_budget_ratio=%.17g",
           maximum_predicted_ratio);
    printf(" maximum_tau_int=%.17g table_pass=%d balance_pass=%d",
           maximum_tau, table_pass, balance_pass);
    printf(" restart_pass=%d resource_limits_hash=%" PRIu64 "\n",
           restart_pass, limits->amplitude_policy_hash);
  }
  pilot_surrogate_table_destroy(&partial_table);
  pilot_table_destroy(&target_table);
  return 1;
}

static int pilot_run_partial_neighbor_candidate(
    int candidate_index, int partial_order, size_t floor_multiplier,
    size_t step_count, int site_count, int qp_total, int sample_count,
    size_t cache_bytes, const uint64_t seeds[PILOT_SEED_COUNT],
    const MVMCKrylovFockModel *model,
    const MVMCKrylovBoundedLimits *limits,
    const MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT],
    const double complex (*raw_values)[PROFILE_DEPTH_COUNT],
    const double norms[PROFILE_DEPTH_COUNT],
    const uint64_t *configurations, size_t sector_size,
    const size_t *lookup, size_t lookup_count) {
  MVMCKrylovPositiveGuidePolicy guide_policy;
  MVMCKrylovMatrixMeasurementPolicy measurement_policy;
  MVMCKrylovPositiveSamplerProposalPolicy neighbor_policy;
  MVMCKrylovPositiveSamplerSurrogatePolicy pilot_policy;
  PilotCandidateTable target_table;
  PilotSurrogateTable partial_table;
  double log_basis_scale[PROFILE_DEPTH_COUNT];
  double eta = NAN;
  double maximum_anchor_relative_residual = NAN;
  double maximum_predicted_ratio = 0.0;
  double maximum_tau = 0.0;
  double maximum_iid_ratio = NAN;
  double required_tau = NAN;
  double inner_row_residual = NAN;
  double inner_db_residual = NAN;
  double powered_row_residual = NAN;
  double powered_db_residual = NAN;
  double powered_ratio_residual = NAN;
  double outer_db_residual = NAN;
  double outer_stationary_residual = NAN;
  MarkovMatrixFamily iid_limiting_family = MARKOV_FAMILY_S;
  size_t iid_limiting_entry = 0;
  uint64_t proposal_model_hash = 0;
  uint64_t partial_policy_hash = 0;
  int table_pass;
  int balance_pass;
  int restart_pass;
  int seed_index;

  memset(&target_table, 0, sizeof(target_table));
  memset(&partial_table, 0, sizeof(partial_table));
  memset(&pilot_policy, 0, sizeof(pilot_policy));
  if (model == NULL || limits == NULL || lookup == NULL ||
      !markov_compute_scales_and_eta(sector_size, raw_values, norms, 1.0e-2,
                                     log_basis_scale, &eta) ||
      !markov_init_guide_policy(site_count, qp_total, cache_bytes, 1.0e-2,
                                eta, log_basis_scale, &guide_policy) ||
      !markov_init_measurement_policy(eta, log_basis_scale,
                                      &measurement_policy) ||
      mvmc_krylov_positive_sampler_proposal_policy_create(
          0, 1, &neighbor_policy) != MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_proposal_model_policy_hash(
          &neighbor_policy, model, &configurations[0], 1,
          &proposal_model_hash) != MVMC_KRYLOV_STATUS_OK) {
    return 0;
  }
  table_pass =
      pilot_table_create(&measurement_policy, exact_values, sector_size,
                         &target_table) &&
      pilot_iid_summary(&target_table, &maximum_iid_ratio, &required_tau,
                        &iid_limiting_family, &iid_limiting_entry) &&
      pilot_partial_table_create(
          &guide_policy, partial_order, floor_multiplier, exact_values,
          raw_values, &target_table, sector_size, &partial_table,
          &maximum_anchor_relative_residual, &partial_policy_hash);
  pilot_policy.step_count = step_count;
  pilot_policy.policy_hash =
      pilot_partial_neighbor_policy_hash(partial_policy_hash, step_count);
  if (pilot_policy.policy_hash == 0) table_pass = 0;
  balance_pass =
      table_pass && pilot_exact_l4_surrogate_balance(
                        model, &pilot_policy, configurations,
                        target_table.guide_weights, partial_table.weights,
                        sector_size, &inner_row_residual, &inner_db_residual,
                        &powered_row_residual, &powered_db_residual,
                        &powered_ratio_residual, &outer_db_residual,
                        &outer_stationary_residual);
  restart_pass =
      table_pass && pilot_surrogate_restart_check(
                        model, &pilot_policy, &target_table, &partial_table,
                        configurations, sector_size, lookup, lookup_count,
                        seeds[0]);
  if (!table_pass || !balance_pass || !restart_pass) {
    pilot_surrogate_table_destroy(&partial_table);
    pilot_table_destroy(&target_table);
    return 0;
  }
  if (world_rank() == 0) {
    printf("PILOT_PARTIAL_NEIGHBOR_CANDIDATE schema=%d candidate_index=%d",
           PILOT_PARTIAL_NEIGHBOR_SCHEMA_VERSION, candidate_index);
    printf(" site_count=%d qp_total=%d sample_count_per_seed=%d",
           site_count, qp_total, sample_count);
    printf(" cache_requested_bytes=%zu rho=0.01 eta=%.17g",
           cache_bytes, eta);
    printf(" partial_order=%d floor_multiplier=%zu step_count=%zu",
           partial_order, floor_multiplier, step_count);
    printf(" inner_kernel=neighbor_only proposal_policy_hash=%" PRIu64,
           neighbor_policy.policy_hash);
    printf(" proposal_model_hash=%" PRIu64, proposal_model_hash);
    printf(" guide_policy_hash=%" PRIu64 " guide_shape_hash=%" PRIu64,
           guide_policy.policy_hash, pilot_guide_shape_hash(&guide_policy));
    printf(" partial_policy_hash=%" PRIu64,
           partial_policy_hash);
    printf(" partial_neighbor_policy_hash=%" PRIu64,
           pilot_policy.policy_hash);
    printf(" target_table_hash=%" PRIu64,
           target_table.table_hash);
    printf(" partial_table_hash=%" PRIu64 " anchor_count=%zu",
           partial_table.table_hash, sector_size);
    printf(" maximum_anchor_relative_residual=%.17g",
           maximum_anchor_relative_residual);
    printf(" partial_minimum=%.17g partial_maximum=%.17g",
           partial_table.minimum_weight, partial_table.maximum_weight);
    printf(" partial_mean=%.17g target_ratio_minimum=%.17g",
           partial_table.mean_weight, partial_table.minimum_target_ratio);
    printf(" target_ratio_maximum=%.17g",
           partial_table.maximum_target_ratio);
    printf(" target_ratio_maximum_mean_ratio=%.17g",
           partial_table.target_ratio_maximum_mean_ratio);
    printf(" target_ratio_ess_fraction=%.17g",
           partial_table.target_ratio_ess_fraction);
    printf(" ideal_acceptance=%.17g ideal_configuration_change=%.17g",
           partial_table.ideal_independence_acceptance,
           partial_table.ideal_independence_configuration_change);
    printf(" maximum_iid_se_budget_ratio=%.17g required_tau_for_0p90=%.17g",
           maximum_iid_ratio, required_tau);
    printf(" iid_limiting_family=%s iid_limiting_entry=%zu",
           markov_family_name(iid_limiting_family), iid_limiting_entry);
    printf(" exact_balance_pass=%d inner_row_residual=%.17g",
           balance_pass, inner_row_residual);
    printf(" inner_db_residual=%.17g powered_row_residual=%.17g",
           inner_db_residual, powered_row_residual);
    printf(" powered_db_residual=%.17g powered_ratio_residual=%.17g",
           powered_db_residual, powered_ratio_residual);
    printf(" outer_db_residual=%.17g outer_stationary_residual=%.17g",
           outer_db_residual, outer_stationary_residual);
    printf(" restart_pass=%d\n", restart_pass);
  }
  for (seed_index = 0; seed_index < PILOT_SEED_COUNT; ++seed_index) {
    PilotChain chain;
    double *series = NULL;
    double seed_maximum_tau = 0.0;
    double seed_maximum_ratio = 0.0;
    MarkovMatrixFamily family;
    if ((size_t)sample_count >
        SIZE_MAX / (PILOT_SERIES_COUNT * sizeof(*series))) {
      pilot_surrogate_table_destroy(&partial_table);
      pilot_table_destroy(&target_table);
      return 0;
    }
    series = (double *)calloc(
        (size_t)sample_count * PILOT_SERIES_COUNT, sizeof(*series));
    if (series == NULL ||
        !pilot_chain_initialize(seeds[seed_index], target_table.guide_weights,
                                sector_size, configurations, &chain) ||
        !pilot_chain_advance_surrogate(
            model, &pilot_policy, &target_table, &partial_table, lookup,
            lookup_count, (size_t)sample_count, series, &chain)) {
      free(series);
      pilot_surrogate_table_destroy(&partial_table);
      pilot_table_destroy(&target_table);
      return 0;
    }
    for (family = MARKOV_FAMILY_S; family <= MARKOV_FAMILY_B; ++family) {
      size_t entry;
      for (entry = 0; entry < 6; ++entry) {
        double component_tau[2];
        double predicted_se;
        double predicted_ratio;
        int component;
        for (component = 0; component < 2; ++component) {
          const size_t series_index =
              pilot_series_index(family, entry, component);
          MVMCKrylovTauIntResult tau;
          if (mvmc_krylov_tau_int_geyer_initial_positive(
                  &series[series_index * (size_t)sample_count],
                  (size_t)sample_count, &tau) != MVMC_KRYLOV_STATUS_OK ||
              !tau.valid) {
            free(series);
            pilot_surrogate_table_destroy(&partial_table);
            pilot_table_destroy(&target_table);
            return 0;
          }
          component_tau[component] = tau.tau_int;
          seed_maximum_tau = fmax(seed_maximum_tau, tau.tau_int);
          if (world_rank() == 0) {
            printf("PILOT_SERIES candidate_index=%d seed_index=%d",
                   candidate_index, seed_index);
            printf(" family=%s entry=%zu component=%s",
                   markov_family_name(family), entry,
                   component == 0 ? "real" : "imag");
            printf(" iid_variance=%.17g tau_int=%.17g\n",
                   target_table.iid_variance[series_index], tau.tau_int);
          }
        }
        predicted_se = sqrt(
            (target_table.iid_variance[
                 pilot_series_index(family, entry, 0)] * component_tau[0] +
             target_table.iid_variance[
                 pilot_series_index(family, entry, 1)] * component_tau[1]) /
            (4096.0 * target_table.mean_denominator *
             target_table.mean_denominator));
        predicted_ratio =
            predicted_se / target_table.budget.entry[family][entry].budget;
        seed_maximum_ratio = fmax(seed_maximum_ratio, predicted_ratio);
        if (world_rank() == 0) {
          printf("PILOT_ENTRY candidate_index=%d seed_index=%d",
                 candidate_index, seed_index);
          printf(" family=%s entry=%zu predicted_se=%.17g",
                 markov_family_name(family), entry, predicted_se);
          printf(" budget=%.17g predicted_se_budget_ratio=%.17g\n",
                 target_table.budget.entry[family][entry].budget,
                 predicted_ratio);
        }
      }
    }
    maximum_tau = fmax(maximum_tau, seed_maximum_tau);
    maximum_predicted_ratio =
        fmax(maximum_predicted_ratio, seed_maximum_ratio);
    if (world_rank() == 0) {
      printf("PILOT_SEED candidate_index=%d seed_index=%d",
             candidate_index, seed_index);
      printf(" seed=%" PRIu64 " seed_hex=0x%016" PRIx64,
             seeds[seed_index], seeds[seed_index]);
      printf(" maximum_tau_int=%.17g maximum_predicted_se_budget_ratio=%.17g",
             seed_maximum_tau, seed_maximum_ratio);
      printf(" accepted=%" PRIu64 " rejected=%" PRIu64,
             chain.accepted, chain.rejected);
      printf(" surrogate_inner_attempted=%" PRIu64,
             chain.surrogate_inner_attempted);
      printf(" surrogate_inner_accepted=%" PRIu64,
             chain.surrogate_inner_accepted);
      printf(" surrogate_inner_rejected=%" PRIu64,
             chain.surrogate_inner_rejected);
      printf(" surrogate_final_changed=%" PRIu64,
             chain.surrogate_final_changed);
      printf(" surrogate_final_self=%" PRIu64,
             chain.surrogate_final_self);
      printf(" trace_hash=%" PRIu64 "\n", chain.trace_hash);
    }
    free(series);
  }
  if (world_rank() == 0) {
    printf("PILOT_DECISION candidate_index=%d", candidate_index);
    printf(" physical_case_eligible=%d",
           maximum_predicted_ratio <= 0.90 && maximum_tau <= 12.0 &&
               table_pass && balance_pass && restart_pass);
    printf(" maximum_predicted_se_budget_ratio=%.17g",
           maximum_predicted_ratio);
    printf(" maximum_tau_int=%.17g table_pass=%d balance_pass=%d",
           maximum_tau, table_pass, balance_pass);
    printf(" restart_pass=%d resource_limits_hash=%" PRIu64 "\n",
           restart_pass, limits->amplitude_policy_hash);
  }
  pilot_surrogate_table_destroy(&partial_table);
  pilot_table_destroy(&target_table);
  return 1;
}

static int pilot_run_partial_callback_candidate(
    int candidate_index, int partial_order, size_t floor_multiplier,
    size_t step_count, int site_count, int qp_total, int sample_count,
    size_t cache_bytes, const uint64_t seeds[PILOT_SEED_COUNT],
    const MVMCKrylovFockModel *model,
    const MVMCKrylovBoundedLimits *limits,
    MVMCKrylovBoundedWorkspace *partial_workspace,
    uint64_t partial_plan_hash,
    ProfileScaledAmplitude *amplitude_workspace,
    const MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT],
    const double complex (*raw_values)[PROFILE_DEPTH_COUNT],
    const double norms[PROFILE_DEPTH_COUNT],
    const uint64_t *configurations, size_t sector_size,
    const size_t *lookup, size_t lookup_count) {
  MVMCKrylovPositiveGuidePolicy guide_policy;
  MVMCKrylovMatrixMeasurementPolicy measurement_policy;
  MVMCKrylovPositiveSamplerProposalPolicy neighbor_policy;
  MVMCKrylovPositiveSamplerSurrogatePolicy callback_policy;
  PilotPartialCallbackContext callback;
  PilotCandidateTable target_table;
  PilotSurrogateTable partial_table;
  double log_basis_scale[PROFILE_DEPTH_COUNT];
  double eta = NAN;
  double maximum_anchor_relative_residual = NAN;
  double maximum_callback_log_residual = NAN;
  double maximum_callback_weight_relative_residual = NAN;
  double inner_row_residual = NAN;
  double inner_db_residual = NAN;
  double powered_row_residual = NAN;
  double powered_db_residual = NAN;
  double powered_ratio_residual = NAN;
  double outer_db_residual = NAN;
  double outer_stationary_residual = NAN;
  uint64_t proposal_model_hash = 0;
  uint64_t partial_policy_hash = 0;
  uint64_t callback_table_hash = 0;
  int table_pass;
  int balance_pass;
  int anchor_pass;
  int restart_pass;
  int replay_pass = 1;
  int seed_index;

  memset(&target_table, 0, sizeof(target_table));
  memset(&partial_table, 0, sizeof(partial_table));
  memset(&callback, 0, sizeof(callback));
  if (model == NULL || limits == NULL || partial_workspace == NULL ||
      partial_plan_hash == 0 || amplitude_workspace == NULL ||
      lookup == NULL || sample_count != 256 || partial_order != 2 ||
      step_count != 32 ||
      !markov_compute_scales_and_eta(sector_size, raw_values, norms, 1.0e-2,
                                     log_basis_scale, &eta) ||
      !markov_init_guide_policy(site_count, qp_total, cache_bytes, 1.0e-2,
                                eta, log_basis_scale, &guide_policy) ||
      !markov_init_measurement_policy(eta, log_basis_scale,
                                      &measurement_policy) ||
      mvmc_krylov_positive_sampler_proposal_policy_create(
          0, 1, &neighbor_policy) != MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_proposal_model_policy_hash(
          &neighbor_policy, model, &configurations[0], 1,
          &proposal_model_hash) != MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_surrogate_partial_policy_create(
          step_count, partial_order, floor_multiplier, &guide_policy,
          &callback_policy) != MVMC_KRYLOV_STATUS_OK) {
    return 0;
  }
  table_pass =
      pilot_table_create(&measurement_policy, exact_values, sector_size,
                         &target_table) &&
      pilot_partial_table_create(
          &guide_policy, partial_order, floor_multiplier, exact_values,
          raw_values, &target_table, sector_size, &partial_table,
          &maximum_anchor_relative_residual, &partial_policy_hash);
  balance_pass =
      table_pass && pilot_exact_l4_surrogate_balance(
                        model, &callback_policy, configurations,
                        target_table.guide_weights, partial_table.weights,
                        sector_size, &inner_row_residual, &inner_db_residual,
                        &powered_row_residual, &powered_db_residual,
                        &powered_ratio_residual, &outer_db_residual,
                        &outer_stationary_residual);
  callback.workspace = partial_workspace;
  callback.amplitude_workspace = amplitude_workspace;
  callback.guide_policy = &guide_policy;
  callback.surrogate_policy = &callback_policy;
  anchor_pass = table_pass && pilot_partial_callback_anchor(
                                  &callback, configurations, sector_size,
                                  &partial_table,
                                  &maximum_callback_log_residual,
                                  &maximum_callback_weight_relative_residual,
                                  &callback_table_hash) &&
                maximum_callback_log_residual <= 1024.0 * DBL_EPSILON &&
                maximum_callback_weight_relative_residual <=
                    1024.0 * DBL_EPSILON;
  callback.evaluation_count = 0;
  callback.terminal_amplitude_calls = 0;
  restart_pass =
      anchor_pass && pilot_partial_callback_restart_check(
                         model, &callback_policy, &target_table, &callback,
                         configurations, sector_size, lookup, lookup_count,
                         seeds[0]);
  if (!table_pass || !balance_pass || !anchor_pass || !restart_pass) {
    pilot_surrogate_table_destroy(&partial_table);
    pilot_table_destroy(&target_table);
    return 0;
  }
  if (world_rank() == 0) {
    printf("PILOT_PARTIAL_CALLBACK_CANDIDATE schema=%d candidate_index=%d",
           PILOT_PARTIAL_CALLBACK_SCHEMA_VERSION, candidate_index);
    printf(" site_count=%d qp_total=%d sample_count_per_seed=%d",
           site_count, qp_total, sample_count);
    printf(" cache_requested_bytes=%zu rho=0.01 eta=%.17g",
           cache_bytes, eta);
    printf(" partial_order=%d floor_multiplier=%zu step_count=%zu",
           partial_order, floor_multiplier, step_count);
    printf(" callback=max_order_2_bounded inner_kernel=neighbor_only");
    printf(" proposal_policy_hash=%" PRIu64, neighbor_policy.policy_hash);
    printf(" proposal_model_hash=%" PRIu64, proposal_model_hash);
    printf(" guide_policy_hash=%" PRIu64 " guide_shape_hash=%" PRIu64,
           guide_policy.policy_hash, pilot_guide_shape_hash(&guide_policy));
    printf(" partial_policy_hash=%" PRIu64,
           partial_policy_hash);
    printf(" callback_policy_hash=%" PRIu64,
           callback_policy.policy_hash);
    printf(" target_table_hash=%" PRIu64,
           target_table.table_hash);
    printf(" partial_table_hash=%" PRIu64,
           partial_table.table_hash);
    printf(" callback_table_hash=%" PRIu64,
           callback_table_hash);
    printf(" partial_plan_hash=%" PRIu64 " anchor_count=%zu",
           partial_plan_hash, sector_size);
    printf(" maximum_anchor_relative_residual=%.17g",
           maximum_anchor_relative_residual);
    printf(" maximum_callback_log_residual=%.17g",
           maximum_callback_log_residual);
    printf(" maximum_callback_weight_relative_residual=%.17g",
           maximum_callback_weight_relative_residual);
    printf(" exact_balance_pass=%d inner_row_residual=%.17g",
           balance_pass, inner_row_residual);
    printf(" inner_db_residual=%.17g powered_row_residual=%.17g",
           inner_db_residual, powered_row_residual);
    printf(" powered_db_residual=%.17g powered_ratio_residual=%.17g",
           powered_db_residual, powered_ratio_residual);
    printf(" outer_db_residual=%.17g outer_stationary_residual=%.17g",
           outer_db_residual, outer_stationary_residual);
    printf(" anchor_pass=%d restart_pass=%d\n",
           anchor_pass, restart_pass);
  }
  callback.evaluation_count = 0;
  callback.terminal_amplitude_calls = 0;
  for (seed_index = 0; seed_index < PILOT_SEED_COUNT; ++seed_index) {
    PilotChain table_chain;
    PilotChain callback_chain;
    double *table_series = NULL;
    double *callback_series = NULL;
    const size_t series_value_count =
        (size_t)sample_count * PILOT_SERIES_COUNT;
    uint64_t evaluations_before = callback.evaluation_count;
    uint64_t terminal_before = callback.terminal_amplitude_calls;
    int seed_replay_pass;
    memset(&table_chain, 0, sizeof(table_chain));
    memset(&callback_chain, 0, sizeof(callback_chain));
    if ((size_t)sample_count >
        SIZE_MAX / (PILOT_SERIES_COUNT * sizeof(*table_series))) {
      replay_pass = 0;
      break;
    }
    table_series = (double *)calloc(series_value_count,
                                    sizeof(*table_series));
    callback_series = (double *)calloc(series_value_count,
                                       sizeof(*callback_series));
    seed_replay_pass =
        table_series != NULL && callback_series != NULL &&
        pilot_chain_initialize(seeds[seed_index], target_table.guide_weights,
                               sector_size, configurations, &table_chain) &&
        pilot_chain_initialize(seeds[seed_index], target_table.guide_weights,
                               sector_size, configurations, &callback_chain) &&
        pilot_chain_advance_surrogate(
            model, &callback_policy, &target_table, &partial_table, lookup,
            lookup_count, (size_t)sample_count, table_series, &table_chain) &&
        pilot_chain_advance_surrogate_callback(
            model, &callback_policy, &target_table, &callback, lookup,
            lookup_count, (size_t)sample_count, callback_series,
            &callback_chain) &&
        pilot_surrogate_chains_equal(&table_chain, &callback_chain) &&
        memcmp(table_series, callback_series,
               series_value_count * sizeof(*table_series)) == 0 &&
        isfinite(callback_chain.surrogate_log_weight) &&
        fabs(callback_chain.surrogate_log_weight -
             partial_table.log_weights[callback_chain.current_index]) <=
            1024.0 * DBL_EPSILON;
    if (world_rank() == 0) {
      printf("PILOT_CALLBACK_SEED candidate_index=%d seed_index=%d",
             candidate_index, seed_index);
      printf(" seed=%" PRIu64 " seed_hex=0x%016" PRIx64,
             seeds[seed_index], seeds[seed_index]);
      printf(" replay_pass=%d accepted=%" PRIu64 " rejected=%" PRIu64,
             seed_replay_pass, callback_chain.accepted,
             callback_chain.rejected);
      printf(" surrogate_inner_attempted=%" PRIu64,
             callback_chain.surrogate_inner_attempted);
      printf(" surrogate_inner_accepted=%" PRIu64,
             callback_chain.surrogate_inner_accepted);
      printf(" surrogate_inner_rejected=%" PRIu64,
             callback_chain.surrogate_inner_rejected);
      printf(" surrogate_final_changed=%" PRIu64,
             callback_chain.surrogate_final_changed);
      printf(" surrogate_final_self=%" PRIu64,
             callback_chain.surrogate_final_self);
      printf(" callback_evaluations=%" PRIu64,
             callback.evaluation_count - evaluations_before);
      printf(" callback_terminal_amplitude_calls=%" PRIu64,
             callback.terminal_amplitude_calls - terminal_before);
      printf(" trace_hash=%" PRIu64 "\n", callback_chain.trace_hash);
    }
    free(callback_series);
    free(table_series);
    if (!seed_replay_pass) replay_pass = 0;
  }
  if (world_rank() == 0) {
    printf("PILOT_CALLBACK_DECISION candidate_index=%d callback_eligible=%d",
           candidate_index, replay_pass && anchor_pass && balance_pass &&
                                restart_pass);
    printf(" table_pass=%d balance_pass=%d anchor_pass=%d",
           table_pass, balance_pass, anchor_pass);
    printf(" restart_pass=%d replay_pass=%d resource_limits_hash=%" PRIu64,
           restart_pass, replay_pass, limits->amplitude_policy_hash);
    printf("\n");
  }
  pilot_surrogate_table_destroy(&partial_table);
  pilot_table_destroy(&target_table);
  return replay_pass;
}

static int run_pilot(int site_count, int qp_total, int sample_count,
                     size_t cache_bytes, PilotMode mode) {
  static const size_t numerators[4] = {1, 1, 1, 1};
  static const size_t denominators[4] = {8, 4, 2, 1};
  static const double rho_values[3] = {1.0e-2, 1.0e-4, 1.0e-6};
  static const size_t shell_distance_numerators[3] = {1, 1, 2};
  static const size_t shell_distance_denominators[3] = {3, 2, 3};
  static const double flat_rho_values[7] = {
      3.0e-2, 1.0e-1, 3.0e-1, 1.0, 3.0, 10.0, 100.0};
  static const size_t surrogate_step_counts[3] = {4, 16, 32};
  static const size_t surrogate_floor_multipliers[3] = {1, 10, 100};
  static const uint64_t global_seeds[PILOT_SEED_COUNT] = {
      UINT64_C(0x5034533350494c31), UINT64_C(0x5034533350494c32),
      UINT64_C(0x5034533350494c33), UINT64_C(0x5034533350494c34)};
  static const uint64_t shell_seeds[PILOT_SEED_COUNT] = {
      UINT64_C(0x5034533450494c31), UINT64_C(0x5034533450494c32),
      UINT64_C(0x5034533450494c33), UINT64_C(0x5034533450494c34)};
  static const uint64_t flat_seeds[PILOT_SEED_COUNT] = {
      UINT64_C(0x5034533550494c31), UINT64_C(0x5034533550494c32),
      UINT64_C(0x5034533550494c33), UINT64_C(0x5034533550494c34)};
  static const uint64_t surrogate_seeds[PILOT_SEED_COUNT] = {
      UINT64_C(0x5034533650494c31), UINT64_C(0x5034533650494c32),
      UINT64_C(0x5034533650494c33), UINT64_C(0x5034533650494c34)};
  static const uint64_t partial_seeds[PILOT_SEED_COUNT] = {
      UINT64_C(0x5034533750494c31), UINT64_C(0x5034533750494c32),
      UINT64_C(0x5034533750494c33), UINT64_C(0x5034533750494c34)};
  static const uint64_t partial_neighbor_seeds[PILOT_SEED_COUNT] = {
      UINT64_C(0x5034533850494c31), UINT64_C(0x5034533850494c32),
      UINT64_C(0x5034533850494c33), UINT64_C(0x5034533850494c34)};
  static const uint64_t long_direct_seeds[PILOT_SEED_COUNT] = {
      UINT64_C(0x5034533950494c31), UINT64_C(0x5034533950494c32),
      UINT64_C(0x5034533950494c33), UINT64_C(0x5034533950494c34)};
  ProfileModel fixture;
  MVMCClassicKrylovModelWorkspace *model_workspace = NULL;
  ProfileScaledAmplitude *amplitude_workspace = NULL;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedPlan *partial_plan = NULL;
  MVMCKrylovBoundedWorkspace *krylov_workspace = NULL;
  MVMCKrylovBoundedWorkspace *partial_workspace = NULL;
  MVMCKrylovBoundedCollectiveWorkspace *collective_workspace = NULL;
  const MVMCKrylovFockModel *model = NULL;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovBoundedLimits partial_limits;
  MVMCKrylovBoundedCollectiveResult collective_result;
  unsigned int *masks = NULL;
  MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT] = NULL;
  double complex (*raw_values)[PROFILE_DEPTH_COUNT] = NULL;
  uint64_t *configurations = NULL;
  size_t *lookup = NULL;
  double *slater = NULL;
  double complex *weights = NULL;
  double norms[PROFILE_DEPTH_COUNT];
  const double projection_parameter = -0.27;
  size_t mask_count = 0;
  size_t sector_size = 0;
  size_t lookup_count = 0;
  uint64_t plan_hash = 0;
  uint64_t partial_plan_hash = 0;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  int local_ready = 1;
  int return_code = 1;
  int fraction_index;
  int rank = world_rank();

  memset(&fixture, 0, sizeof(fixture));
  mask_count = fixed_masks(site_count, &masks);
  local_ready = mask_count != 0 &&
                checked_multiply_size(mask_count, mask_count, &sector_size) &&
                initialize_model(&fixture, site_count) &&
                initialize_slater(site_count, qp_total, &slater, &weights);
  if (local_ready) {
    exact_values = (MVMCScaledComplex (*)[PROFILE_DEPTH_COUNT])calloc(
        sector_size, sizeof(*exact_values));
    raw_values = (double complex (*)[PROFILE_DEPTH_COUNT])calloc(
        sector_size, sizeof(*raw_values));
    configurations =
        (uint64_t *)calloc(sector_size, sizeof(*configurations));
    local_ready = exact_values != NULL && raw_values != NULL &&
                  configurations != NULL;
  }
  if (!collective_all_ready(local_ready)) goto cleanup;
  status = mvmc_classic_krylov_model_workspace_create_from_root(
      rank == 0 ? &fixture.raw : NULL, world_communicator(),
      &model_workspace);
  if (!collective_all_ready(status == MVMC_KRYLOV_STATUS_OK)) goto cleanup;
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
  if (local_ready && mode == PILOT_MODE_PARTIAL_CALLBACK) {
    partial_limits = limits;
    partial_limits.max_order = 2;
    status = mvmc_bounded_krylov_plan_create(
        model, &partial_limits, &partial_plan);
    local_ready = status == MVMC_KRYLOV_STATUS_OK && partial_plan != NULL;
  }
  if (local_ready && mode == PILOT_MODE_PARTIAL_CALLBACK) {
    partial_plan_hash = mvmc_bounded_krylov_plan_hash(partial_plan);
    status = mvmc_bounded_krylov_workspace_create(
        partial_plan, &partial_workspace);
    local_ready = status == MVMC_KRYLOV_STATUS_OK &&
                  partial_workspace != NULL && partial_plan_hash != 0;
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
    goto cleanup;
  }
  if (!markov_evaluate_exact_sector(
          krylov_workspace, amplitude_workspace, masks, mask_count,
          sector_size, site_count, exact_values, raw_values,
          configurations, norms) ||
      !pilot_lookup_create(configurations, sector_size, site_count, &lookup,
                           &lookup_count)) {
    goto cleanup;
  }
  if (rank == 0) {
    const int schema =
        mode == PILOT_MODE_LONG_DIRECT
            ? PILOT_LONG_DIRECT_SCHEMA_VERSION
            : (mode == PILOT_MODE_PARTIAL_CALLBACK
            ? PILOT_PARTIAL_CALLBACK_SCHEMA_VERSION
            : (mode == PILOT_MODE_PARTIAL_NEIGHBOR
            ? PILOT_PARTIAL_NEIGHBOR_SCHEMA_VERSION
            : (mode == PILOT_MODE_PARTIAL_GUIDE
                   ? PILOT_PARTIAL_SCHEMA_VERSION
            : (mode == PILOT_MODE_SURROGATE
                   ? PILOT_SURROGATE_SCHEMA_VERSION
            : (mode == PILOT_MODE_FLAT_GUIDE
                   ? PILOT_FLAT_SCHEMA_VERSION
                   : (mode == PILOT_MODE_SHELL
                          ? PILOT_SHELL_SCHEMA_VERSION
                          : PILOT_SCHEMA_VERSION))))));
    const char *fixture_name =
        mode == PILOT_MODE_LONG_DIRECT
            ? "p4s9_long_direct_exact_table"
            : (mode == PILOT_MODE_PARTIAL_CALLBACK
            ? "p4s7_partial_guide_real_callback_stage_b2"
            : (mode == PILOT_MODE_PARTIAL_NEIGHBOR
            ? "p4s7_partial_guide_neighbor_stage_b1"
            : (mode == PILOT_MODE_PARTIAL_GUIDE
                   ? "p4s7_partial_guide_exact_q_prefilter"
            : (mode == PILOT_MODE_SURROGATE
                   ? "p4s6_surrogate_exact_table"
            : (mode == PILOT_MODE_FLAT_GUIDE
                   ? "p4s5_flat_guide_exact_table"
                   : (mode == PILOT_MODE_SHELL
                          ? "p4s4_fixed_distance_exact_table"
                          : "p4s3_exact_table"))))));
    printf("PILOT_RUN schema=%d fixture=%s",
           schema, fixture_name);
    printf(" site_count=%d qp_total=%d rank_count=%d",
           site_count, qp_total, world_size());
    printf(" sample_count_per_seed=%d cache_requested_bytes=%zu",
           sample_count, cache_bytes);
    printf(" sector_size=%zu anchor_count=%zu plan_hash=%" PRIu64 "\n",
           sector_size, sector_size, plan_hash);
  }
  if (mode == PILOT_MODE_LONG_DIRECT) {
    MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
    if (mvmc_krylov_positive_sampler_proposal_policy_create(
            0, 1, &proposal_policy) != MVMC_KRYLOV_STATUS_OK ||
        !pilot_run_candidate(
            0, site_count, qp_total, sample_count, cache_bytes, 1.0e-2,
            &proposal_policy, long_direct_seeds,
            PILOT_LONG_DIRECT_SCHEMA_VERSION, 32768, 0.80, 12.0,
            model, &limits, exact_values, raw_values, norms, configurations,
            sector_size, lookup, lookup_count)) {
      goto cleanup;
    }
  } else if (mode == PILOT_MODE_PARTIAL_CALLBACK) {
    static const size_t stage_b2_alphas[2] = {1, 10};
    int alpha_index;
    for (alpha_index = 0; alpha_index < 2; ++alpha_index) {
      if (!pilot_run_partial_callback_candidate(
              alpha_index, 2, stage_b2_alphas[alpha_index], 32,
              site_count, qp_total, sample_count, cache_bytes,
              partial_neighbor_seeds, model, &partial_limits,
              partial_workspace, partial_plan_hash, amplitude_workspace,
              exact_values, raw_values, norms, configurations, sector_size,
              lookup, lookup_count)) {
        goto cleanup;
      }
    }
  } else if (mode == PILOT_MODE_PARTIAL_NEIGHBOR) {
    static const size_t stage_b_alphas[2] = {1, 10};
    static const size_t stage_b_steps[2] = {16, 32};
    int alpha_index;
    for (alpha_index = 0; alpha_index < 2; ++alpha_index) {
      int step_index;
      for (step_index = 0; step_index < 2; ++step_index) {
        const int candidate_index = alpha_index * 2 + step_index;
        if (!pilot_run_partial_neighbor_candidate(
                candidate_index, 2, stage_b_alphas[alpha_index],
                stage_b_steps[step_index], site_count, qp_total,
                sample_count, cache_bytes, partial_neighbor_seeds, model,
                &limits, exact_values, raw_values, norms, configurations,
                sector_size, lookup, lookup_count)) {
          goto cleanup;
        }
      }
    }
  } else if (mode == PILOT_MODE_PARTIAL_GUIDE) {
    int partial_order;
    for (partial_order = 1; partial_order <= 2; ++partial_order) {
      int alpha_index;
      for (alpha_index = 0; alpha_index < 3; ++alpha_index) {
        const int candidate_index = (partial_order - 1) * 3 + alpha_index;
        if (!pilot_run_partial_candidate(
                candidate_index, partial_order,
                surrogate_floor_multipliers[alpha_index], site_count,
                qp_total, sample_count, cache_bytes, partial_seeds, &limits,
                exact_values, raw_values, norms, configurations,
                sector_size)) {
          goto cleanup;
        }
      }
    }
  } else if (mode == PILOT_MODE_SURROGATE) {
    int step_index;
    for (fraction_index = 0; fraction_index < 3; ++fraction_index) {
      if (!pilot_emit_surrogate_control(
              fraction_index, surrogate_floor_multipliers[fraction_index],
              site_count, qp_total, cache_bytes, exact_values, raw_values,
              norms, sector_size)) {
        goto cleanup;
      }
    }
    for (step_index = 0; step_index < 3; ++step_index) {
      int alpha_index;
      for (alpha_index = 0; alpha_index < 3; ++alpha_index) {
        const int candidate_index = step_index * 3 + alpha_index;
        if (!pilot_run_surrogate_candidate(
                candidate_index, surrogate_step_counts[step_index],
                surrogate_floor_multipliers[alpha_index], site_count,
                qp_total, sample_count, cache_bytes, surrogate_seeds, model,
                &limits, exact_values, raw_values, norms, configurations,
                sector_size, lookup, lookup_count)) {
          goto cleanup;
        }
      }
    }
  } else if (mode == PILOT_MODE_FLAT_GUIDE) {
    MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
    if (!pilot_emit_flat_control(
            0, "rho_0p01", 0, site_count, qp_total, cache_bytes, 1.0e-2,
            exact_values, raw_values, norms, sector_size) ||
        !pilot_emit_flat_control(
            1, "uniform_guide_limit", 1, site_count, qp_total, cache_bytes,
            0.0, exact_values, raw_values, norms, sector_size) ||
        mvmc_krylov_positive_sampler_proposal_policy_create(
            1, 1, &proposal_policy) != MVMC_KRYLOV_STATUS_OK) {
      goto cleanup;
    }
    for (fraction_index = 0; fraction_index < 7; ++fraction_index) {
      if (!pilot_run_candidate(
              fraction_index, site_count, qp_total, sample_count,
              cache_bytes, flat_rho_values[fraction_index], &proposal_policy,
              flat_seeds, PILOT_FLAT_SCHEMA_VERSION, 4096, 0.90, 12.0,
              model, &limits,
              exact_values, raw_values, norms, configurations, sector_size,
              lookup, lookup_count)) {
        goto cleanup;
      }
    }
  } else if (mode == PILOT_MODE_SHELL) {
    for (fraction_index = 0; fraction_index < 3; ++fraction_index) {
      MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
      if (mvmc_krylov_positive_sampler_shell_policy_create(
              1, 8, shell_distance_numerators[fraction_index],
              shell_distance_denominators[fraction_index],
              &proposal_policy) != MVMC_KRYLOV_STATUS_OK ||
          !pilot_run_candidate(
              fraction_index, site_count, qp_total, sample_count,
              cache_bytes, 1.0e-2, &proposal_policy, shell_seeds,
              PILOT_SHELL_SCHEMA_VERSION, 4096, 0.90, 12.0,
              model, &limits, exact_values, raw_values, norms,
              configurations, sector_size, lookup, lookup_count)) {
        goto cleanup;
      }
    }
  } else {
    for (fraction_index = 0; fraction_index < 4; ++fraction_index) {
      int rho_index;
      for (rho_index = 0; rho_index < 3; ++rho_index) {
        const int candidate_index = fraction_index * 3 + rho_index;
        MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
        if (mvmc_krylov_positive_sampler_proposal_policy_create(
                numerators[fraction_index], denominators[fraction_index],
                &proposal_policy) != MVMC_KRYLOV_STATUS_OK ||
            !pilot_run_candidate(
                candidate_index, site_count, qp_total, sample_count,
                cache_bytes, rho_values[rho_index], &proposal_policy,
                global_seeds, PILOT_SCHEMA_VERSION, 4096, 0.90, 12.0,
                model, &limits,
                exact_values, raw_values, norms, configurations,
                sector_size, lookup, lookup_count)) {
          goto cleanup;
        }
      }
    }
  }
  return_code = 0;

cleanup:
  free(lookup);
  free(configurations);
  free(raw_values);
  free(exact_values);
  destroy_profile_amplitude(amplitude_workspace);
  mvmc_bounded_krylov_collective_workspace_destroy(collective_workspace);
  mvmc_bounded_krylov_workspace_destroy(krylov_workspace);
  mvmc_bounded_krylov_workspace_destroy(partial_workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
  mvmc_bounded_krylov_plan_destroy(partial_plan);
  mvmc_classic_krylov_model_workspace_destroy(model_workspace);
  free(weights);
  free(slater);
  destroy_model(&fixture);
  free(masks);
  return return_code;
}

#ifndef MVMC_BOUNDED_KRYLOV_MARKOV_PILOT_DRIVER_NO_MAIN
int main(int argc, char **argv) {
  int site_count;
  int qp_total;
  int sample_count;
  size_t cache_bytes;
  PilotMode mode = PILOT_MODE_GLOBAL;
  int result;
  (void)&run_markov;
  (void)&parse_markov_double;
  (void)&parse_markov_u64;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#endif
  if ((argc != 5 && argc != 6) ||
      !parse_int(argv[1], 4, 8, &site_count) || (site_count & 1) != 0 ||
      !parse_int(argv[2], 1, 4, &qp_total) ||
      (qp_total != 1 && qp_total != 4) ||
      !parse_int(argv[3], 32, 100000000, &sample_count) ||
      !parse_size_arg(argv[4], (size_t)4 * 1024 * 1024 * 1024,
                      &cache_bytes) ||
      (argc == 6 && strcmp(argv[5], "shell") != 0 &&
       strcmp(argv[5], "flat") != 0 &&
       strcmp(argv[5], "surrogate") != 0 &&
       strcmp(argv[5], "partial") != 0 &&
       strcmp(argv[5], "partial-neighbor") != 0 &&
       strcmp(argv[5], "partial-callback") != 0 &&
       strcmp(argv[5], "long-direct") != 0)) {
    if (world_rank() == 0) {
      fprintf(stderr,
              "usage: %s EVEN_SITE_COUNT(4..8) QP_TOTAL(1|4) "
              "SAMPLES_PER_SEED CACHE_BYTES "
              "[shell|flat|surrogate|partial|partial-neighbor|"
              "partial-callback|long-direct]\n",
              argv[0]);
    }
    result = EXIT_FAILURE;
  } else {
    if (argc == 6) {
      mode = strcmp(argv[5], "long-direct") == 0
                 ? PILOT_MODE_LONG_DIRECT
                 : (strcmp(argv[5], "shell") == 0
                 ? PILOT_MODE_SHELL
                 : (strcmp(argv[5], "flat") == 0
                        ? PILOT_MODE_FLAT_GUIDE
                        : (strcmp(argv[5], "surrogate") == 0
                               ? PILOT_MODE_SURROGATE
                               : (strcmp(argv[5], "partial") == 0
                                      ? PILOT_MODE_PARTIAL_GUIDE
                                      : (strcmp(argv[5],
                                                "partial-neighbor") == 0
                                             ? PILOT_MODE_PARTIAL_NEIGHBOR
                                             : PILOT_MODE_PARTIAL_CALLBACK)))));
    }
    result = run_pilot(site_count, qp_total, sample_count, cache_bytes,
                       mode) == 0
                 ? EXIT_SUCCESS
                 : EXIT_FAILURE;
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return result;
}
#endif
