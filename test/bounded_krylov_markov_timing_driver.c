/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
at your option any later version.
*/

#define MVMC_BOUNDED_KRYLOV_MARKOV_PILOT_DRIVER_NO_MAIN
#include "bounded_krylov_markov_pilot_driver.c"

#include <time.h>

enum {
  TIMING_SCHEMA_VERSION = 1,
  TIMING_FRACTION_COUNT = 5,
  TIMING_RHO_COUNT = 3,
  TIMING_REPEAT_COUNT = 7,
  TIMING_COMPONENT_COUNT = 5
};

enum {
  TIMING_PARTIAL_SCHEMA_VERSION = 2,
  TIMING_PARTIAL_COMPONENT_COUNT = 7
};

enum {
  TIMING_SESSION_SCHEMA_VERSION = 3,
  TIMING_SESSION_COMPONENT_COUNT = 8
};

static double timing_wall_seconds(void) {
  struct timespec now;
  if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) return -1.0;
  return (double)now.tv_sec + 1.0e-9 * (double)now.tv_nsec;
}

static double timing_elapsed_seconds(double start) {
  const double end = timing_wall_seconds();
  return start >= 0.0 && end >= start ? end - start : -1.0;
}

static uint64_t timing_amplitude_generation_hash(
    const ProfileScaledAmplitude *workspace) {
  uint64_t hash = UINT64_C(1469598103934665603);
  size_t index;
  size_t slater_count;
  if (workspace == NULL ||
      !checked_multiply_size((size_t)workspace->qp_total,
                             workspace->slater_matrix_count,
                             &slater_count)) {
    return 0;
  }
  markov_hash_u64(&hash, UINT64_C(0x5034533953455353));
  markov_hash_u64(&hash, (uint64_t)workspace->site_count);
  markov_hash_u64(&hash, (uint64_t)workspace->up_electron_count);
  markov_hash_u64(&hash, (uint64_t)workspace->down_electron_count);
  markov_hash_u64(&hash, (uint64_t)(unsigned int)workspace->qp_total);
  markov_hash_u64(&hash, markov_double_bits(workspace->projection_parameter));
  markov_hash_u64(&hash,
                  markov_double_bits(workspace->scaled_pivot_tolerance));
  for (index = 0; index < (size_t)workspace->qp_total; ++index) {
    markov_hash_u64(&hash,
                    markov_double_bits(creal(workspace->global_weights[index])));
    markov_hash_u64(&hash,
                    markov_double_bits(cimag(workspace->global_weights[index])));
  }
  for (index = 0; index < slater_count; ++index) {
    markov_hash_u64(&hash, markov_double_bits(workspace->slater[index]));
  }
  return hash == 0 ? UINT64_C(1) : hash;
}

static int timing_guide_equal(
    const MVMCKrylovPositiveGuideResult *left,
    const MVMCKrylovPositiveGuideResult *right) {
  return left != NULL && right != NULL && left->valid == right->valid &&
         left->status == right->status &&
         fabs(left->log_guide - right->log_guide) <=
             256.0 * DBL_EPSILON * fmax(1.0, fabs(right->log_guide));
}

static int timing_acceptance_equal(
    const MVMCKrylovPositiveGuideAcceptance *left,
    const MVMCKrylovPositiveGuideAcceptance *right) {
  return left != NULL && right != NULL && left->valid == right->valid &&
         left->status == right->status &&
         fabs(left->log_acceptance_ratio - right->log_acceptance_ratio) <=
             512.0 * DBL_EPSILON *
                 fmax(1.0, fabs(right->log_acceptance_ratio)) &&
         left->uniform == right->uniform &&
         left->accepted == right->accepted;
}

static int timing_run_repeat(
    int site_count, int qp_total, int sample_count, size_t cache_bytes,
    int fraction_index, int rho_index, uint64_t stream,
    const MVMCKrylovFockModel *model, const MVMCKrylovBoundedPlan *plan,
    ProfileScaledAmplitude *amplitude_workspace,
    const MVMCKrylovPositiveGuidePolicy *guide_policy,
    const PilotCandidateTable *table, const uint64_t *configurations,
    size_t sector_size, const size_t *lookup, size_t lookup_count,
    size_t global_numerator,
    size_t global_denominator, double rho, int persistent_session,
    uint64_t amplitude_generation_hash) {
  const uint64_t seed = persistent_session
                            ? UINT64_C(0x5034533954494d45)
                            : UINT64_C(0x5034533352305631);
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
  MVMCKrylovPositiveSamplerSnapshot current;
  MVMCKrylovPositiveSamplerRng rng;
  double initial_uniform = NAN;
  double local_timing[TIMING_SESSION_COMPONENT_COUNT] = {
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double rank_max_timing[TIMING_SESSION_COMPONENT_COUNT] = {
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  uint64_t current_words = 0;
  uint64_t accepted = 0;
  uint64_t rejected = 0;
  uint64_t neighbor_attempted = 0;
  uint64_t global_attempted = 0;
  uint64_t global_self = 0;
  uint64_t proposal_model_hash = 0;
  uint64_t session_cache_hits = 0;
  uint64_t session_cache_misses = 0;
  uint64_t session_cache_insertions = 0;
  uint64_t session_cache_evictions = 0;
  uint64_t session_cache_resets = 0;
  uint64_t callback_calls_start = 0;
  uint64_t allocated_capacity_bytes = 0;
  uint64_t rank_max_allocated_capacity_bytes = 0;
  uint64_t rank_max_peak_rss_bytes = 0;
  uint64_t statistical_trace_hash = UINT64_C(1469598103934665603);
  size_t session_resident_peak = 0;
  size_t initial_index;
  size_t exact_current_index;
  double chain_start = -1.0;
  double component_start = -1.0;
  int sample;

  if (mvmc_krylov_positive_sampler_proposal_policy_create(
          global_numerator, global_denominator, &proposal_policy) !=
          MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_proposal_model_hash(
          model, &configurations[0], 1, &proposal_model_hash) !=
          MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_rng_seed(seed, stream, &rng) !=
          MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_rng_draw_uniform(
          &rng, &initial_uniform) != MVMC_KRYLOV_STATUS_OK) {
    return 0;
  }
  initial_index = markov_inverse_cdf(
      initial_uniform, table->guide_weights, sector_size, NULL);
  exact_current_index = initial_index;
  current_words = configurations[initial_index];
  if (mvmc_bounded_krylov_workspace_create(plan, &workspace) !=
          MVMC_KRYLOV_STATUS_OK ||
      workspace == NULL) {
    mvmc_bounded_krylov_workspace_destroy(workspace);
    return 0;
  }

  chain_start = timing_wall_seconds();
  callback_calls_start = amplitude_workspace->callback_calls;
  if (persistent_session) {
    component_start = timing_wall_seconds();
    if (amplitude_generation_hash == 0 ||
        mvmc_bounded_krylov_session_begin(
            workspace, profile_scaled_amplitude, amplitude_workspace,
            amplitude_generation_hash) != MVMC_KRYLOV_STATUS_OK) {
      mvmc_bounded_krylov_workspace_destroy(workspace);
      return 0;
    }
    local_timing[5] = timing_elapsed_seconds(component_start);
  }
  component_start = timing_wall_seconds();
  if (mvmc_krylov_positive_sampler_initialize(
          workspace, guide_policy, &current_words, 1,
          profile_scaled_amplitude, amplitude_workspace, &current) !=
          MVMC_KRYLOV_STATUS_OK ||
      !current.valid ||
      fabs(current.guide.log_guide - table->guides[initial_index].log_guide) >
          256.0 * DBL_EPSILON *
              fmax(1.0, fabs(table->guides[initial_index].log_guide))) {
    mvmc_bounded_krylov_workspace_destroy(workspace);
    return 0;
  }
  if (persistent_session &&
      !timing_guide_equal(&current.guide, &table->guides[initial_index])) {
    fprintf(stderr,
            "P4-S9 session initial exact-guide mismatch stream=%" PRIu64
            " index=%zu actual=%.17g expected=%.17g\n",
            stream, initial_index, current.guide.log_guide,
            table->guides[initial_index].log_guide);
    mvmc_bounded_krylov_workspace_destroy(workspace);
    return 0;
  }
  local_timing[6] = timing_elapsed_seconds(component_start);
  if (persistent_session) {
    int depth;
    const MVMCKrylovBoundedStatistics *statistics = &current.krylov.statistics;
    size_t allocated = 0;
    if (!statistics->persistent_session_active ||
        statistics->amplitude_generation_hash != amplitude_generation_hash ||
        statistics->session_root_evaluation != 1 ||
        statistics->cache_reset_performed != 1) {
      mvmc_bounded_krylov_workspace_destroy(workspace);
      return 0;
    }
    if (!checked_add_size(statistics->plan_bytes, statistics->model_bytes,
                          &allocated) ||
        !checked_add_size(allocated, statistics->workspace_bytes,
                          &allocated) ||
        !checked_add_size(
            allocated,
            mvmc_bounded_krylov_collective_workspace_bytes(
                amplitude_workspace->collective),
            &allocated) ||
        !checked_add_size(allocated,
                          profile_amplitude_workspace_bytes(
                              amplitude_workspace),
                          &allocated)) {
      mvmc_bounded_krylov_workspace_destroy(workspace);
      return 0;
    }
    allocated_capacity_bytes = (uint64_t)allocated;
    session_cache_resets += statistics->cache_reset_performed;
    session_cache_insertions += statistics->cache_insertions;
    session_cache_evictions += statistics->cache_evictions;
    session_resident_peak = statistics->cache_entries_peak;
    for (depth = 0; depth <= MVMC_KRYLOV_MAX_ORDER; ++depth) {
      session_cache_hits += statistics->cache_hits[depth];
      session_cache_misses += statistics->cache_misses[depth];
    }
  }
  for (sample = 0; sample < sample_count; ++sample) {
    MVMCKrylovPositiveSamplerProposalStepResult step;
    MVMCKrylovPositiveGuideAcceptance exact_acceptance;
    uint64_t proposal_words = 0;
    size_t proposal_index = SIZE_MAX;
    MVMCKrylovStatus status =
        mvmc_krylov_positive_sampler_step_mixture_rng(
            workspace, guide_policy, model, &proposal_policy,
            &current_words, 1, &current, &rng, profile_scaled_amplitude,
            amplitude_workspace, &proposal_words, 1, &step);
    if (!collective_all_ready(
            status == MVMC_KRYLOV_STATUS_OK && step.valid &&
            isfinite(step.proposal_seconds) &&
            step.proposal_seconds >= 0.0 &&
            isfinite(step.component_selection_seconds) &&
            step.component_selection_seconds >= 0.0 &&
            isfinite(step.global_subset_seconds) &&
            step.global_subset_seconds >= 0.0 &&
            isfinite(step.bounded_evaluation_seconds) &&
            step.bounded_evaluation_seconds >= 0.0 &&
            isfinite(step.total_step_seconds) &&
            step.total_step_seconds >= 0.0)) {
      mvmc_bounded_krylov_workspace_destroy(workspace);
      return 0;
    }
    if (persistent_session) {
      if (lookup == NULL || proposal_words >= lookup_count ||
          lookup[proposal_words] == SIZE_MAX) {
        mvmc_bounded_krylov_workspace_destroy(workspace);
        return 0;
      }
      proposal_index = lookup[proposal_words];
      status = mvmc_krylov_positive_guide_acceptance(
          &table->guides[exact_current_index],
          &table->guides[proposal_index], step.log_proposal_ratio,
          step.uniform, &exact_acceptance);
      if (!timing_guide_equal(&step.step.proposal_guide,
                              &table->guides[proposal_index]) ||
          status != MVMC_KRYLOV_STATUS_OK ||
          !timing_acceptance_equal(&step.step.acceptance,
                                   &exact_acceptance)) {
        fprintf(stderr,
                "P4-S9 session exact-table mismatch stream=%" PRIu64
                " sample=%d current=%zu proposal=%zu guide=%.17g/%.17g"
                " accepted=%d/%d status=%d\n",
                stream, sample, exact_current_index, proposal_index,
                step.step.proposal_guide.log_guide,
                table->guides[proposal_index].log_guide,
                step.step.acceptance.accepted, exact_acceptance.accepted,
                (int)status);
        mvmc_bounded_krylov_workspace_destroy(workspace);
        return 0;
      }
      if (exact_acceptance.accepted) exact_current_index = proposal_index;
      if (current_words != configurations[exact_current_index]) {
        mvmc_bounded_krylov_workspace_destroy(workspace);
        return 0;
      }
      markov_hash_u64(&statistical_trace_hash, current_words);
      markov_hash_u64(&statistical_trace_hash,
                      (uint64_t)(unsigned int)step.component);
      markov_hash_u64(&statistical_trace_hash,
                      (uint64_t)(unsigned int)exact_acceptance.accepted);
      markov_hash_u64(&statistical_trace_hash, rng.state);
    }
    if (persistent_session) {
      int depth;
      const MVMCKrylovBoundedStatistics *statistics =
          &step.step.proposal_krylov.statistics;
      if (!statistics->persistent_session_active ||
          statistics->amplitude_generation_hash !=
              amplitude_generation_hash ||
          statistics->session_root_evaluation != (uint64_t)sample + 2 ||
          statistics->cache_reset_performed != 0) {
        mvmc_bounded_krylov_workspace_destroy(workspace);
        return 0;
      }
      if (statistics->cache_insertions >
              UINT64_MAX - session_cache_insertions ||
          statistics->cache_evictions >
              UINT64_MAX - session_cache_evictions ||
          statistics->cache_reset_performed >
              UINT64_MAX - session_cache_resets) {
        mvmc_bounded_krylov_workspace_destroy(workspace);
        return 0;
      }
      session_cache_insertions += statistics->cache_insertions;
      session_cache_evictions += statistics->cache_evictions;
      session_cache_resets += statistics->cache_reset_performed;
      if (statistics->cache_entries_peak > session_resident_peak) {
        session_resident_peak = statistics->cache_entries_peak;
      }
      for (depth = 0; depth <= MVMC_KRYLOV_MAX_ORDER; ++depth) {
        if (statistics->cache_hits[depth] >
                UINT64_MAX - session_cache_hits ||
            statistics->cache_misses[depth] >
                UINT64_MAX - session_cache_misses) {
          mvmc_bounded_krylov_workspace_destroy(workspace);
          return 0;
        }
        session_cache_hits += statistics->cache_hits[depth];
        session_cache_misses += statistics->cache_misses[depth];
      }
    }
    local_timing[0] += step.proposal_seconds;
    local_timing[1] += step.component_selection_seconds;
    local_timing[2] += step.global_subset_seconds;
    local_timing[3] += step.bounded_evaluation_seconds;
    local_timing[4] += step.total_step_seconds;
    if (step.step.accepted) {
      ++accepted;
    } else {
      ++rejected;
    }
    if (step.component == MVMC_KRYLOV_PROPOSAL_COMPONENT_GLOBAL) {
      ++global_attempted;
      if (step.self_proposal) ++global_self;
    } else {
      ++neighbor_attempted;
    }
  }
  if (persistent_session &&
      mvmc_bounded_krylov_session_end(workspace) !=
          MVMC_KRYLOV_STATUS_OK) {
    mvmc_bounded_krylov_workspace_destroy(workspace);
    return 0;
  }
  local_timing[7] = timing_elapsed_seconds(chain_start);
  if (persistent_session &&
      (mvmc_bounded_krylov_collective_max_u64(
           amplitude_workspace->collective, allocated_capacity_bytes,
           &rank_max_allocated_capacity_bytes) != MVMC_KRYLOV_STATUS_OK ||
       mvmc_bounded_krylov_collective_max_u64(
           amplitude_workspace->collective, resident_set_bytes(),
           &rank_max_peak_rss_bytes) != MVMC_KRYLOV_STATUS_OK)) {
    mvmc_bounded_krylov_workspace_destroy(workspace);
    return 0;
  }
  if (!markov_reduce_max_times(local_timing, rank_max_timing,
                               persistent_session
                                   ? TIMING_SESSION_COMPONENT_COUNT
                                   : TIMING_COMPONENT_COUNT)) {
    mvmc_bounded_krylov_workspace_destroy(workspace);
    return 0;
  }
  if (world_rank() == 0) {
    printf("%s schema=%d site_count=%d qp_total=%d",
           persistent_session ? "TIMING_SESSION_REPEAT" : "TIMING_REPEAT",
           persistent_session ? TIMING_SESSION_SCHEMA_VERSION
                              : TIMING_SCHEMA_VERSION,
           site_count, qp_total);
    printf(" sample_count=%d cache_requested_bytes=%zu", sample_count,
           cache_bytes);
    printf(" fraction_index=%d rho_index=%d rho=%.17g",
           fraction_index, rho_index, rho);
    printf(" global_numerator=%zu global_denominator=%zu",
           global_numerator, global_denominator);
    printf(" seed=%" PRIu64 " seed_hex=0x%016" PRIx64,
           seed, seed);
    printf(" rng_stream=%" PRIu64, stream);
    printf(" proposal_policy_hash=%" PRIu64,
           proposal_policy.policy_hash);
    printf(" proposal_model_hash=%" PRIu64,
           proposal_model_hash);
    printf(" proposal_seconds=%.17g component_selection_seconds=%.17g",
           rank_max_timing[0], rank_max_timing[1]);
    printf(" global_subset_seconds=%.17g bounded_evaluation_seconds=%.17g",
           rank_max_timing[2], rank_max_timing[3]);
    printf(" total_step_seconds=%.17g total_seconds_per_step=%.17g",
           rank_max_timing[4], rank_max_timing[4] / (double)sample_count);
    printf(" accepted=%" PRIu64 " rejected=%" PRIu64,
           accepted, rejected);
    printf(" neighbor_attempted=%" PRIu64 " global_attempted=%" PRIu64,
           neighbor_attempted, global_attempted);
    printf(" global_self=%" PRIu64 " final_rng_draws=%" PRIu64 "\n",
           global_self, rng.draws);
    if (persistent_session) {
      printf("TIMING_SESSION_RESOURCE schema=%d site_count=%d qp_total=%d",
             TIMING_SESSION_SCHEMA_VERSION, site_count, qp_total);
      printf(" sample_count=%d cache_requested_bytes=%zu", sample_count,
             cache_bytes);
      printf(" rng_stream=%" PRIu64, stream);
      printf(" amplitude_generation_hash=%" PRIu64,
             amplitude_generation_hash);
      printf(" session_begin_seconds=%.17g initialization_seconds=%.17g",
             rank_max_timing[5], rank_max_timing[6]);
      printf(" combined_session_seconds=%.17g seconds_per_configuration=%.17g",
             rank_max_timing[7], rank_max_timing[7] / (double)sample_count);
      printf(" session_root_evaluations=%d cache_resets=%" PRIu64,
             sample_count + 1, session_cache_resets);
      printf(" cache_hits=%" PRIu64 " cache_misses=%" PRIu64,
             session_cache_hits, session_cache_misses);
      printf(" cache_insertions=%" PRIu64 " cache_evictions=%" PRIu64,
             session_cache_insertions, session_cache_evictions);
      printf(" cache_resident_peak=%zu callback_calls=%" PRIu64,
             session_resident_peak,
             amplitude_workspace->callback_calls - callback_calls_start);
      printf(" allocated_capacity_bytes=%" PRIu64 " peak_rss_bytes=%" PRIu64
             "\n",
             rank_max_allocated_capacity_bytes, rank_max_peak_rss_bytes);
      printf("TIMING_SESSION_CORRECTNESS schema=%d site_count=%d qp_total=%d",
             TIMING_SESSION_SCHEMA_VERSION, site_count, qp_total);
      printf(" rng_stream=%" PRIu64 " exact_table_match=1", stream);
      printf(" final_sector_index=%zu final_configuration=%" PRIu64,
             exact_current_index, current_words);
      printf(" statistical_trace_hash=%" PRIu64 "\n",
             statistical_trace_hash);
    }
  }
  mvmc_bounded_krylov_workspace_destroy(workspace);
  return 1;
}

static int timing_run_partial_callback_repeat(
    int candidate_index, size_t floor_multiplier, int site_count,
    int qp_total, int sample_count, size_t cache_bytes, uint64_t stream,
    const MVMCKrylovFockModel *model,
    const MVMCKrylovBoundedPlan *full_plan,
    const MVMCKrylovBoundedPlan *partial_plan,
    uint64_t partial_plan_hash,
    ProfileScaledAmplitude *amplitude_workspace,
    const MVMCKrylovPositiveGuidePolicy *guide_policy,
    const PilotCandidateTable *table, const uint64_t *configurations,
    size_t sector_size) {
  const uint64_t seed = UINT64_C(0x5034533352305631);
  MVMCKrylovBoundedWorkspace *full_workspace = NULL;
  MVMCKrylovBoundedWorkspace *partial_workspace = NULL;
  MVMCKrylovPositiveSamplerProposalPolicy neighbor_policy;
  MVMCKrylovPositiveSamplerSurrogatePolicy surrogate_policy;
  MVMCKrylovPositiveSamplerSnapshot current;
  MVMCKrylovPositiveSamplerRng rng;
  PilotPartialCallbackContext callback;
  double initial_uniform = NAN;
  double current_log_weight = NAN;
  double local_timing[TIMING_PARTIAL_COMPONENT_COUNT] = {
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double rank_max_timing[TIMING_PARTIAL_COMPONENT_COUNT] = {
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  uint64_t current_words = 0;
  uint64_t accepted = 0;
  uint64_t rejected = 0;
  uint64_t final_changed = 0;
  uint64_t skipped_full_evaluations = 0;
  uint64_t proposal_model_hash = 0;
  size_t initial_index;
  int sample;

  memset(&callback, 0, sizeof(callback));
  if (model == NULL || full_plan == NULL || partial_plan == NULL ||
      partial_plan_hash == 0 || amplitude_workspace == NULL ||
      guide_policy == NULL || table == NULL || configurations == NULL ||
      sample_count <= 0 ||
      mvmc_krylov_positive_sampler_proposal_policy_create(
          0, 1, &neighbor_policy) != MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_surrogate_partial_policy_create(
          32, 2, floor_multiplier, guide_policy, &surrogate_policy) !=
          MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_proposal_model_policy_hash(
          &neighbor_policy, model, &configurations[0], 1,
          &proposal_model_hash) != MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_rng_seed(seed, stream, &rng) !=
          MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_rng_draw_uniform(
          &rng, &initial_uniform) != MVMC_KRYLOV_STATUS_OK) {
    return 0;
  }
  initial_index = markov_inverse_cdf(
      initial_uniform, table->guide_weights, sector_size, NULL);
  current_words = configurations[initial_index];
  if (mvmc_bounded_krylov_workspace_create(
          full_plan, &full_workspace) != MVMC_KRYLOV_STATUS_OK ||
      full_workspace == NULL ||
      mvmc_bounded_krylov_workspace_create(
          partial_plan, &partial_workspace) != MVMC_KRYLOV_STATUS_OK ||
      partial_workspace == NULL ||
      mvmc_krylov_positive_sampler_initialize(
          full_workspace, guide_policy, &current_words, 1,
          profile_scaled_amplitude, amplitude_workspace, &current) !=
          MVMC_KRYLOV_STATUS_OK ||
      !current.valid ||
      fabs(current.guide.log_guide -
           table->guides[initial_index].log_guide) >
          256.0 * DBL_EPSILON *
              fmax(1.0, fabs(table->guides[initial_index].log_guide))) {
    mvmc_bounded_krylov_workspace_destroy(partial_workspace);
    mvmc_bounded_krylov_workspace_destroy(full_workspace);
    return 0;
  }
  callback.workspace = partial_workspace;
  callback.amplitude_workspace = amplitude_workspace;
  callback.guide_policy = guide_policy;
  callback.surrogate_policy = &surrogate_policy;
  if (pilot_partial_log_weight_callback(
          &current_words, 1, &callback,
          &current_log_weight) != MVMC_KRYLOV_STATUS_OK) {
    mvmc_bounded_krylov_workspace_destroy(partial_workspace);
    mvmc_bounded_krylov_workspace_destroy(full_workspace);
    return 0;
  }
  memset(callback.depth_seconds, 0, sizeof(callback.depth_seconds));
  callback.evaluation_count = 0;
  callback.terminal_amplitude_calls = 0;
  callback.evaluation_seconds = 0.0;
  for (sample = 0; sample < sample_count; ++sample) {
      MVMCKrylovPositiveSamplerSurrogateDrawResult draw;
      MVMCKrylovPositiveSamplerStepResult step;
      MVMCKrylovPositiveGuideAcceptance self_acceptance;
      uint64_t proposal_words = 0;
      double outer_uniform = NAN;
      double final_start = -1.0;
      double final_seconds = 0.0;
      double evaluation_before = callback.evaluation_seconds;
      double total_start = timing_wall_seconds();
      MVMCKrylovStatus status =
          mvmc_krylov_positive_sampler_draw_surrogate_rng(
              model, &surrogate_policy, &current_words, 1,
              current_log_weight,
              &rng, pilot_partial_log_weight_callback, &callback,
              &proposal_words, 1, &draw);
      if (status != MVMC_KRYLOV_STATUS_OK || !draw.valid) {
        mvmc_bounded_krylov_workspace_destroy(partial_workspace);
        mvmc_bounded_krylov_workspace_destroy(full_workspace);
        return 0;
      }
      rng = draw.rng_after;
      if (mvmc_krylov_positive_sampler_rng_draw_uniform(
              &rng, &outer_uniform) != MVMC_KRYLOV_STATUS_OK) {
        mvmc_bounded_krylov_workspace_destroy(partial_workspace);
        mvmc_bounded_krylov_workspace_destroy(full_workspace);
        return 0;
      }
      if (proposal_words == current_words) {
        status = mvmc_krylov_positive_guide_acceptance(
            &current.guide, &current.guide, draw.log_proposal_ratio,
            outer_uniform, &self_acceptance);
        if (status != MVMC_KRYLOV_STATUS_OK || !self_acceptance.valid ||
            !self_acceptance.accepted || draw.final_configuration_changed) {
          mvmc_bounded_krylov_workspace_destroy(partial_workspace);
          mvmc_bounded_krylov_workspace_destroy(full_workspace);
          return 0;
        }
        ++accepted;
        ++skipped_full_evaluations;
      } else {
        final_start = timing_wall_seconds();
        status = mvmc_krylov_positive_sampler_step(
            full_workspace, guide_policy, &current_words, 1, &current,
            &proposal_words, 1, draw.log_proposal_ratio, outer_uniform,
            profile_scaled_amplitude, amplitude_workspace, &step);
        final_seconds = timing_elapsed_seconds(final_start);
        if (status != MVMC_KRYLOV_STATUS_OK || !step.valid ||
            !isfinite(final_seconds) || final_seconds < 0.0) {
          mvmc_bounded_krylov_workspace_destroy(partial_workspace);
          mvmc_bounded_krylov_workspace_destroy(full_workspace);
          return 0;
        }
        ++final_changed;
        if (step.accepted) {
          current_log_weight = draw.final_log_weight;
          ++accepted;
        } else {
          ++rejected;
        }
      }
      local_timing[4] += fmax(
          0.0, draw.proposal_seconds -
                   (callback.evaluation_seconds - evaluation_before));
      local_timing[5] += final_seconds;
      {
        const double total_seconds = timing_elapsed_seconds(total_start);
        if (!isfinite(total_seconds) || total_seconds < 0.0) {
          mvmc_bounded_krylov_workspace_destroy(partial_workspace);
          mvmc_bounded_krylov_workspace_destroy(full_workspace);
          return 0;
        }
        local_timing[6] += total_seconds;
      }
  }
  local_timing[0] = callback.depth_seconds[0];
  local_timing[1] = callback.depth_seconds[1];
  local_timing[2] = callback.depth_seconds[2];
  local_timing[3] = callback.evaluation_seconds;
  if (!markov_reduce_max_times(local_timing, rank_max_timing,
                               TIMING_PARTIAL_COMPONENT_COUNT)) {
    mvmc_bounded_krylov_workspace_destroy(partial_workspace);
    mvmc_bounded_krylov_workspace_destroy(full_workspace);
    return 0;
  }
  if (world_rank() == 0) {
    printf("TIMING_PARTIAL_CALLBACK_REPEAT schema=%d candidate_index=%d",
           TIMING_PARTIAL_SCHEMA_VERSION, candidate_index);
    printf(" site_count=%d qp_total=%d sample_count=%d",
           site_count, qp_total, sample_count);
    printf(" cache_requested_bytes=%zu partial_order=2", cache_bytes);
    printf(" floor_multiplier=%zu step_count=32", floor_multiplier);
    printf(" seed=%" PRIu64 " seed_hex=0x%016" PRIx64, seed, seed);
    printf(" rng_stream=%" PRIu64, stream);
    printf(" proposal_policy_hash=%" PRIu64,
           neighbor_policy.policy_hash);
    printf(" proposal_model_hash=%" PRIu64,
           proposal_model_hash);
    printf(" callback_policy_hash=%" PRIu64,
           surrogate_policy.policy_hash);
    printf(" partial_plan_hash=%" PRIu64, partial_plan_hash);
    printf(" partial_depth0_seconds=%.17g partial_depth1_seconds=%.17g",
           rank_max_timing[0], rank_max_timing[1]);
    printf(" partial_depth2_seconds=%.17g partial_evaluation_seconds=%.17g",
           rank_max_timing[2], rank_max_timing[3]);
    printf(" inner_overhead_seconds=%.17g final_bounded_seconds=%.17g",
           rank_max_timing[4], rank_max_timing[5]);
    printf(" total_step_seconds=%.17g total_seconds_per_step=%.17g",
           rank_max_timing[6], rank_max_timing[6] / (double)sample_count);
    printf(" accepted=%" PRIu64 " rejected=%" PRIu64,
           accepted, rejected);
    printf(" final_changed=%" PRIu64 " skipped_full_evaluations=%" PRIu64,
           final_changed, skipped_full_evaluations);
    printf(" callback_evaluations=%" PRIu64,
           callback.evaluation_count);
    printf(" callback_terminal_amplitude_calls=%" PRIu64,
           callback.terminal_amplitude_calls);
    printf(" final_rng_draws=%" PRIu64 "\n", rng.draws);
  }
  mvmc_bounded_krylov_workspace_destroy(partial_workspace);
  mvmc_bounded_krylov_workspace_destroy(full_workspace);
  return 1;
}

static int run_timing(int site_count, int qp_total, int sample_count,
                      size_t cache_bytes, int persistent_session) {
  static const size_t numerators[TIMING_FRACTION_COUNT] = {0, 1, 1, 1, 1};
  static const size_t denominators[TIMING_FRACTION_COUNT] = {1, 8, 4, 2, 1};
  static const double rho_values[TIMING_RHO_COUNT] = {
      1.0e-2, 1.0e-4, 1.0e-6};
  ProfileModel fixture;
  MVMCClassicKrylovModelWorkspace *model_workspace = NULL;
  ProfileScaledAmplitude *amplitude_workspace = NULL;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *calibration_workspace = NULL;
  MVMCKrylovBoundedCollectiveWorkspace *collective_workspace = NULL;
  const MVMCKrylovFockModel *model = NULL;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovBoundedCollectiveResult collective_result;
  unsigned int *masks = NULL;
  MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT] = NULL;
  double complex (*raw_values)[PROFILE_DEPTH_COUNT] = NULL;
  const MVMCScaledComplex (*exact_values_readonly)[PROFILE_DEPTH_COUNT] =
      NULL;
  const double complex (*raw_values_readonly)[PROFILE_DEPTH_COUNT] = NULL;
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
  uint64_t amplitude_generation_hash = 0;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  int local_ready = 1;
  int return_code = 1;
  int rank = world_rank();
  int rho_index;

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
    exact_values_readonly =
        (const MVMCScaledComplex (*)[PROFILE_DEPTH_COUNT])exact_values;
    raw_values_readonly =
        (const double complex (*)[PROFILE_DEPTH_COUNT])raw_values;
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
    status =
        mvmc_bounded_krylov_workspace_create(plan, &calibration_workspace);
    local_ready = status == MVMC_KRYLOV_STATUS_OK &&
                  calibration_workspace != NULL;
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
  if (local_ready && persistent_session) {
    amplitude_generation_hash =
        timing_amplitude_generation_hash(amplitude_workspace);
    local_ready = amplitude_generation_hash != 0;
    status = local_ready ? MVMC_KRYLOV_STATUS_OK
                         : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  status = mvmc_bounded_krylov_collective_synchronize(
      collective_workspace, status, plan_hash, local_ready,
      &collective_result);
  if (!collective_all_ready(status == MVMC_KRYLOV_STATUS_OK &&
                            collective_result.valid && local_ready)) {
    goto cleanup;
  }
  if (!markov_evaluate_exact_sector(
          calibration_workspace, amplitude_workspace, masks, mask_count,
          sector_size, site_count, exact_values, raw_values,
          configurations, norms) ||
      !pilot_lookup_create(configurations, sector_size, site_count, &lookup,
                           &lookup_count)) {
    goto cleanup;
  }
  mvmc_bounded_krylov_workspace_destroy(calibration_workspace);
  calibration_workspace = NULL;
  if (rank == 0) {
    printf("%s schema=%d fixture=%s",
           persistent_session ? "TIMING_SESSION_RUN" : "TIMING_RUN",
           persistent_session ? TIMING_SESSION_SCHEMA_VERSION
                              : TIMING_SCHEMA_VERSION,
           persistent_session ? "p4s9_long_direct_session_timing"
                              : "p4s3_matched_timing");
    printf(" site_count=%d qp_total=%d rank_count=%d", site_count,
           qp_total, world_size());
    printf(" sample_count=%d repeat_count=%d cache_requested_bytes=%zu",
           sample_count, TIMING_REPEAT_COUNT, cache_bytes);
    printf(" sector_size=%zu plan_hash=%" PRIu64 "\n", sector_size,
           plan_hash);
  }
  for (rho_index = 0;
       rho_index < (persistent_session ? 1 : TIMING_RHO_COUNT);
       ++rho_index) {
    double log_basis_scale[PROFILE_DEPTH_COUNT];
    double eta = NAN;
    MVMCKrylovMatrixMeasurementPolicy measurement_policy;
    PilotCandidateTable table;
    int fraction_index;
    memset(&table, 0, sizeof(table));
    if (!markov_compute_scales_and_eta(
            sector_size, raw_values_readonly, norms, rho_values[rho_index],
            log_basis_scale, &eta) ||
        !markov_init_measurement_policy(eta, log_basis_scale,
                                        &measurement_policy) ||
        !pilot_table_create(&measurement_policy, exact_values_readonly,
                            sector_size, &table)) {
      pilot_table_destroy(&table);
      goto cleanup;
    }
    for (fraction_index = 0;
         fraction_index <
             (persistent_session ? 1 : TIMING_FRACTION_COUNT);
         ++fraction_index) {
      uint64_t stream;
      for (stream = 0; stream < TIMING_REPEAT_COUNT; ++stream) {
        MVMCKrylovPositiveGuidePolicy guide_policy;
        int order;
        memset(&guide_policy, 0, sizeof(guide_policy));
        guide_policy.order = MARKOV_ORDER;
        guide_policy.eta = eta;
        for (order = 0; order <= MARKOV_ORDER; ++order) {
          guide_policy.lambda[order] = 1.0;
          guide_policy.log_basis_scale[order] = log_basis_scale[order];
        }
        guide_policy.policy_hash = markov_policy_hash(
            site_count, qp_total, cache_bytes, rho_values[rho_index], eta,
            log_basis_scale);
        (void)guide_policy;
        if (!timing_run_repeat(
                site_count, qp_total, sample_count, cache_bytes,
                fraction_index, rho_index, stream, model, plan,
                amplitude_workspace, &guide_policy, &table,
                configurations, sector_size, lookup, lookup_count,
                numerators[fraction_index], denominators[fraction_index],
                rho_values[rho_index], persistent_session,
                amplitude_generation_hash)) {
          pilot_table_destroy(&table);
          goto cleanup;
        }
      }
    }
    pilot_table_destroy(&table);
  }
  return_code = 0;

cleanup:
  free(lookup);
  free(configurations);
  free(raw_values);
  free(exact_values);
  destroy_profile_amplitude(amplitude_workspace);
  mvmc_bounded_krylov_collective_workspace_destroy(collective_workspace);
  mvmc_bounded_krylov_workspace_destroy(calibration_workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
  mvmc_classic_krylov_model_workspace_destroy(model_workspace);
  free(weights);
  free(slater);
  destroy_model(&fixture);
  free(masks);
  return return_code;
}

static int run_partial_callback_timing(
    int site_count, int qp_total, int sample_count, size_t cache_bytes) {
  static const size_t floor_multipliers[2] = {1, 10};
  ProfileModel fixture;
  MVMCClassicKrylovModelWorkspace *model_workspace = NULL;
  ProfileScaledAmplitude *amplitude_workspace = NULL;
  MVMCKrylovBoundedPlan *full_plan = NULL;
  MVMCKrylovBoundedPlan *partial_plan = NULL;
  MVMCKrylovBoundedWorkspace *calibration_workspace = NULL;
  MVMCKrylovBoundedCollectiveWorkspace *collective_workspace = NULL;
  const MVMCKrylovFockModel *model = NULL;
  MVMCKrylovBoundedLimits full_limits;
  MVMCKrylovBoundedLimits partial_limits;
  MVMCKrylovBoundedCollectiveResult collective_result;
  MVMCKrylovPositiveGuidePolicy guide_policy;
  MVMCKrylovMatrixMeasurementPolicy measurement_policy;
  PilotCandidateTable table;
  unsigned int *masks = NULL;
  MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT] = NULL;
  double complex (*raw_values)[PROFILE_DEPTH_COUNT] = NULL;
  const MVMCScaledComplex (*exact_values_readonly)[PROFILE_DEPTH_COUNT] =
      NULL;
  const double complex (*raw_values_readonly)[PROFILE_DEPTH_COUNT] = NULL;
  uint64_t *configurations = NULL;
  double *slater = NULL;
  double complex *weights = NULL;
  double norms[PROFILE_DEPTH_COUNT];
  double log_basis_scale[PROFILE_DEPTH_COUNT];
  double eta = NAN;
  const double projection_parameter = -0.27;
  size_t mask_count = 0;
  size_t sector_size = 0;
  uint64_t full_plan_hash = 0;
  uint64_t partial_plan_hash = 0;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  int local_ready = 1;
  int return_code = 1;
  int rank = world_rank();
  uint64_t stream;

  memset(&fixture, 0, sizeof(fixture));
  memset(&table, 0, sizeof(table));
  mask_count = fixed_masks(site_count, &masks);
  local_ready = sample_count == 128 && mask_count != 0 &&
                checked_multiply_size(mask_count, mask_count, &sector_size) &&
                initialize_model(&fixture, site_count) &&
                initialize_slater(site_count, qp_total, &slater, &weights);
  if (local_ready) {
    exact_values = (MVMCScaledComplex (*)[PROFILE_DEPTH_COUNT])calloc(
        sector_size, sizeof(*exact_values));
    raw_values = (double complex (*)[PROFILE_DEPTH_COUNT])calloc(
        sector_size, sizeof(*raw_values));
    exact_values_readonly =
        (const MVMCScaledComplex (*)[PROFILE_DEPTH_COUNT])exact_values;
    raw_values_readonly =
        (const double complex (*)[PROFILE_DEPTH_COUNT])raw_values;
    configurations =
        (uint64_t *)calloc(sector_size, sizeof(*configurations));
    local_ready = exact_values != NULL && raw_values != NULL &&
                  configurations != NULL;
  }
  if (!collective_all_ready(local_ready)) goto cleanup_partial;
  status = mvmc_classic_krylov_model_workspace_create_from_root(
      rank == 0 ? &fixture.raw : NULL, world_communicator(),
      &model_workspace);
  if (!collective_all_ready(status == MVMC_KRYLOV_STATUS_OK)) {
    goto cleanup_partial;
  }
  model = mvmc_classic_krylov_model(model_workspace);
  local_ready = model != NULL &&
                profile_limits(cache_bytes, model, &full_limits);
  if (local_ready) {
    partial_limits = full_limits;
    partial_limits.max_order = 2;
    status = mvmc_bounded_krylov_plan_create(
        model, &full_limits, &full_plan);
    local_ready = status == MVMC_KRYLOV_STATUS_OK && full_plan != NULL;
  }
  if (local_ready) {
    status = mvmc_bounded_krylov_plan_create(
        model, &partial_limits, &partial_plan);
    local_ready = status == MVMC_KRYLOV_STATUS_OK && partial_plan != NULL;
  }
  if (local_ready) {
    full_plan_hash = mvmc_bounded_krylov_plan_hash(full_plan);
    partial_plan_hash = mvmc_bounded_krylov_plan_hash(partial_plan);
    status = mvmc_bounded_krylov_workspace_create(
        full_plan, &calibration_workspace);
    local_ready = status == MVMC_KRYLOV_STATUS_OK &&
                  calibration_workspace != NULL && full_plan_hash != 0 &&
                  partial_plan_hash != 0;
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
      collective_workspace, status, full_plan_hash, local_ready,
      &collective_result);
  if (!collective_all_ready(status == MVMC_KRYLOV_STATUS_OK &&
                            collective_result.valid && local_ready)) {
    goto cleanup_partial;
  }
  if (!markov_evaluate_exact_sector(
          calibration_workspace, amplitude_workspace, masks, mask_count,
          sector_size, site_count, exact_values, raw_values,
          configurations, norms) ||
      !markov_compute_scales_and_eta(
          sector_size, raw_values_readonly, norms, 1.0e-2,
          log_basis_scale, &eta) ||
      !markov_init_guide_policy(site_count, qp_total, cache_bytes, 1.0e-2,
                                eta, log_basis_scale, &guide_policy) ||
      !markov_init_measurement_policy(eta, log_basis_scale,
                                      &measurement_policy) ||
      !pilot_table_create(&measurement_policy, exact_values_readonly,
                          sector_size, &table)) {
    goto cleanup_partial;
  }
  mvmc_bounded_krylov_workspace_destroy(calibration_workspace);
  calibration_workspace = NULL;
  if (rank == 0) {
    printf("TIMING_PARTIAL_RUN schema=%d ",
           TIMING_PARTIAL_SCHEMA_VERSION);
    printf("fixture=p4s7_partial_callback_matched_timing");
    printf(" site_count=%d qp_total=%d rank_count=%d",
           site_count, qp_total, world_size());
    printf(" sample_count=%d repeat_count=%d cache_requested_bytes=%zu",
           sample_count, TIMING_REPEAT_COUNT, cache_bytes);
    printf(" sector_size=%zu full_plan_hash=%" PRIu64,
           sector_size, full_plan_hash);
    printf(" partial_plan_hash=%" PRIu64 "\n", partial_plan_hash);
  }
  for (stream = 0; stream < TIMING_REPEAT_COUNT; ++stream) {
    int candidate_index;
    if (!timing_run_repeat(
            site_count, qp_total, sample_count, cache_bytes,
            4, 0, stream, model, full_plan, amplitude_workspace,
            &guide_policy, &table, configurations, sector_size,
            NULL, 0, 1, 1, 1.0e-2, 0, 0)) {
      goto cleanup_partial;
    }
    for (candidate_index = 0; candidate_index < 2; ++candidate_index) {
      if (!timing_run_partial_callback_repeat(
              candidate_index, floor_multipliers[candidate_index],
              site_count, qp_total, sample_count, cache_bytes, stream,
              model, full_plan, partial_plan, partial_plan_hash,
              amplitude_workspace, &guide_policy, &table,
              configurations, sector_size)) {
        goto cleanup_partial;
      }
    }
  }
  return_code = 0;

cleanup_partial:
  pilot_table_destroy(&table);
  free(configurations);
  free(raw_values);
  free(exact_values);
  destroy_profile_amplitude(amplitude_workspace);
  mvmc_bounded_krylov_collective_workspace_destroy(collective_workspace);
  mvmc_bounded_krylov_workspace_destroy(calibration_workspace);
  mvmc_bounded_krylov_plan_destroy(partial_plan);
  mvmc_bounded_krylov_plan_destroy(full_plan);
  mvmc_classic_krylov_model_workspace_destroy(model_workspace);
  free(weights);
  free(slater);
  destroy_model(&fixture);
  free(masks);
  return return_code;
}

int main(int argc, char **argv) {
  int site_count;
  int qp_total;
  int sample_count = 128;
  size_t cache_bytes;
  int partial_callback_mode = 0;
  int persistent_session_mode = 0;
  int result;
  (void)&run_pilot;
  (void)&run_markov;
  (void)&parse_markov_double;
  (void)&parse_markov_u64;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#endif
  if ((argc != 4 && argc != 5 && argc != 6) ||
      !parse_int(argv[1], 4, 8, &site_count) || (site_count & 1) != 0 ||
      !parse_int(argv[2], 1, 4, &qp_total) ||
      (qp_total != 1 && qp_total != 4) ||
      !parse_size_arg(argv[3], (size_t)4 * 1024 * 1024 * 1024,
                      &cache_bytes) ||
      ((argc == 5 || argc == 6) &&
       !parse_int(argv[4], 1, 32768, &sample_count)) ||
      (argc == 6 && strcmp(argv[5], "partial-callback") != 0 &&
       strcmp(argv[5], "long-direct-session") != 0) ||
      (argc != 6 && sample_count > 128) ||
      (argc == 6 && strcmp(argv[5], "partial-callback") == 0 &&
       sample_count != 128)) {
    if (world_rank() == 0) {
      fprintf(stderr,
              "usage: %s EVEN_SITE_COUNT(4..8) QP_TOTAL(1|4) "
              "CACHE_BYTES [SAMPLE_COUNT "
              "[partial-callback|long-direct-session]]\n",
              argv[0]);
    }
    result = EXIT_FAILURE;
  } else {
    partial_callback_mode =
        argc == 6 && strcmp(argv[5], "partial-callback") == 0;
    persistent_session_mode =
        argc == 6 && strcmp(argv[5], "long-direct-session") == 0;
    result = (partial_callback_mode
                  ? run_partial_callback_timing(
                        site_count, qp_total, sample_count, cache_bytes)
                  : run_timing(site_count, qp_total, sample_count,
                               cache_bytes, persistent_session_mode)) == 0
                 ? EXIT_SUCCESS
                 : EXIT_FAILURE;
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return result;
}
