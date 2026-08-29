#include "krylov_positive_sampler.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(condition, ...)                                                 \
  do {                                                                        \
    if (!(condition)) {                                                       \
      fprintf(stderr, "KrylovPositiveSampler_Unit FAIL: ");                  \
      fprintf(stderr, __VA_ARGS__);                                           \
      fprintf(stderr, " (line %d)\n", __LINE__);                            \
      ++failures;                                                             \
    }                                                                         \
  } while (0)

typedef struct {
  double complex values[64];
  MVMCKrylovStatus forced_status;
  uint64_t calls;
} TableAmplitude;

typedef struct {
  double weights[64];
  MVMCKrylovStatus forced_status;
  size_t fail_on_call;
  size_t calls;
} SurrogateTable;

typedef struct {
  uint64_t trace_hash;
  uint64_t final_words;
  uint64_t accepted;
  uint64_t rejected;
  MVMCKrylovPositiveSamplerRng rng;
  MVMCKrylovPositiveSamplerSnapshot current;
} SessionTrace;

static int close_double(double actual, double expected, double tolerance) {
  return fabs(actual - expected) <= tolerance * fmax(1.0, fabs(expected));
}

static uint64_t hash_u64_value(uint64_t hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    hash ^= value & UINT64_C(0xff);
    hash *= UINT64_C(1099511628211);
    value >>= 8;
  }
  return hash;
}

static size_t replacement_distance_word(uint64_t current,
                                        uint64_t proposal,
                                        size_t site_count) {
  const uint64_t mask =
      (UINT64_C(1) << site_count) - UINT64_C(1);
  uint64_t removed = (current & ~proposal & mask) |
                     (((current >> site_count) &
                       ~(proposal >> site_count) & mask)
                      << site_count);
  size_t distance = 0;
  while (removed != 0) {
    distance += (size_t)(removed & UINT64_C(1));
    removed >>= 1;
  }
  return distance;
}

static MVMCKrylovStatus table_callback(
    const uint64_t *words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
  TableAmplitude *table = (TableAmplitude *)context;
  const size_t index = word_count == 1 ? (size_t)(words[0] & UINT64_C(63))
                                      : 63;
  ++table->calls;
  if (table->forced_status != MVMC_KRYLOV_STATUS_OK) {
    return table->forced_status;
  }
  memset(result, 0, sizeof(*result));
  if (creal(table->values[index]) == 0.0 &&
      cimag(table->values[index]) == 0.0) {
    if (mvmc_scaled_complex_make_exact_zero(&result->value) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    result->exact_zero_component_count = 1;
  } else if (mvmc_scaled_complex_from_raw_testing(
                 table->values[index], &result->value) !=
             MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  } else {
    result->well_pivoted_component_count = 1;
  }
  result->local_factorization_count = 1;
  result->global_factorization_count = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus surrogate_table_callback(
    const uint64_t *words, size_t word_count, void *context,
    double *log_weight) {
  SurrogateTable *table = (SurrogateTable *)context;
  size_t index;
  if (words == NULL || word_count == 0 || table == NULL ||
      log_weight == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  ++table->calls;
  if (table->forced_status != MVMC_KRYLOV_STATUS_OK) {
    return table->forced_status;
  }
  if (table->fail_on_call != 0 && table->calls == table->fail_on_call) {
    return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  index = (size_t)(words[0] & UINT64_C(63));
  if (!isfinite(table->weights[index]) || table->weights[index] <= 0.0) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  *log_weight = log(table->weights[index]);
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovFockModel diagonal_model(
    const MVMCKrylovHamiltonianTerm *terms,
    const MVMCKrylovFermionOperator *operators) {
  MVMCKrylovFockModel model;
  memset(&model, 0, sizeof(model));
  model.site_count = 2;
  model.up_electron_count = 1;
  model.hermitian = 1;
  model.terms = terms;
  model.term_count = 2;
  model.operators = operators;
  model.operator_count = 4;
  return model;
}

static MVMCKrylovFockModel proposal_model(void) {
  MVMCKrylovFockModel model;
  memset(&model, 0, sizeof(model));
  model.site_count = 2;
  model.up_electron_count = 1;
  model.hermitian = 1;
  return model;
}

static MVMCKrylovBoundedLimits bounded_limits(void) {
  MVMCKrylovBoundedLimits limits;
  limits.amplitude_policy_hash = UINT64_C(0x50344453414d504c);
  limits.cache_bytes = 4096;
  limits.max_row_transitions = 16;
  limits.max_workspace_bytes = (size_t)1024 * 1024;
  limits.max_node_expansions = 1000;
  limits.max_terminal_amplitude_calls = 1000;
  limits.max_total_row_transitions = 10000;
  limits.max_order = 1;
  return limits;
}

static MVMCKrylovPositiveGuidePolicy guide_policy(void) {
  MVMCKrylovPositiveGuidePolicy policy;
  int order;
  memset(&policy, 0, sizeof(policy));
  policy.order = 1;
  policy.eta = 0.25;
  policy.policy_hash = UINT64_C(0x5034444755494445);
  for (order = 0; order <= MVMC_KRYLOV_MAX_ORDER; ++order) {
    policy.lambda[order] = 1.0;
  }
  return policy;
}

static void with_workspace(
    void (*body)(MVMCKrylovBoundedWorkspace *, TableAmplitude *)) {
  const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
  };
  const MVMCKrylovHamiltonianTerm terms[] = {
      {2.0, 0, 2, 0, 0},
      {0.5, 2, 2, 0, 1},
  };
  const MVMCKrylovFockModel model = diagonal_model(terms, operators);
  const MVMCKrylovBoundedLimits limits = bounded_limits();
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  TableAmplitude amplitude;
  memset(&amplitude, 0, sizeof(amplitude));
  amplitude.values[1] = 1.0;
  amplitude.values[2] = 4.0;
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_OK,
        "sampler plan create");
  CHECK(plan != NULL &&
            mvmc_bounded_krylov_workspace_create(plan, &workspace) ==
                MVMC_KRYLOV_STATUS_OK &&
            workspace != NULL,
        "sampler workspace create");
  if (workspace != NULL) {
    body(workspace, &amplitude);
  }
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
}

static void with_workspace_and_plan(
    void (*body)(MVMCKrylovBoundedWorkspace *, TableAmplitude *,
                 uint64_t, const MVMCKrylovBoundedLimits *)) {
  const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
  };
  const MVMCKrylovHamiltonianTerm terms[] = {
      {2.0, 0, 2, 0, 0},
      {0.5, 2, 2, 0, 1},
  };
  const MVMCKrylovFockModel model = diagonal_model(terms, operators);
  const MVMCKrylovBoundedLimits limits = bounded_limits();
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  TableAmplitude amplitude;
  memset(&amplitude, 0, sizeof(amplitude));
  amplitude.values[1] = 1.0;
  amplitude.values[2] = 4.0;
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_OK,
        "manifest plan create");
  CHECK(plan != NULL &&
            mvmc_bounded_krylov_workspace_create(plan, &workspace) ==
                MVMC_KRYLOV_STATUS_OK &&
            workspace != NULL,
        "manifest workspace create");
  if (workspace != NULL && plan != NULL) {
    body(workspace, &amplitude, mvmc_bounded_krylov_plan_hash(plan),
         &limits);
  }
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
}

static int run_session_trace(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude,
    const MVMCKrylovPositiveGuidePolicy *policy, uint64_t generation,
    size_t step_count, size_t split_after, SessionTrace *trace) {
  const MVMCKrylovFockModel model = proposal_model();
  MVMCKrylovPositiveSamplerSnapshot current;
  MVMCKrylovPositiveSamplerRng rng;
  uint64_t current_words = UINT64_C(1);
  uint64_t hash = UINT64_C(1469598103934665603);
  uint64_t accepted = 0;
  uint64_t rejected = 0;
  size_t step_index;
  if (workspace == NULL || amplitude == NULL || policy == NULL ||
      trace == NULL || split_after > step_count ||
      mvmc_bounded_krylov_session_begin(
          workspace, table_callback, amplitude, generation) !=
          MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_initialize(
          workspace, policy, &current_words, 1, table_callback, amplitude,
          &current) != MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_rng_seed(
          UINT64_C(0x5034533952535452), UINT64_C(7), &rng) !=
          MVMC_KRYLOV_STATUS_OK) {
    return 0;
  }
  for (step_index = 0; step_index < step_count; ++step_index) {
    MVMCKrylovPositiveSamplerProposalStepResult step;
    uint64_t proposal_words = 0;
    if (step_index == split_after && split_after != step_count) {
      if (mvmc_bounded_krylov_session_end(workspace) !=
              MVMC_KRYLOV_STATUS_OK ||
          mvmc_bounded_krylov_session_begin(
              workspace, table_callback, amplitude, generation) !=
              MVMC_KRYLOV_STATUS_OK) {
        return 0;
      }
    }
    if (mvmc_krylov_positive_sampler_step_rng(
            workspace, policy, &model, &current_words, 1, &current, &rng,
            table_callback, amplitude, &proposal_words, 1, &step) !=
            MVMC_KRYLOV_STATUS_OK ||
        !step.valid) {
      return 0;
    }
    if (step.step.accepted) {
      ++accepted;
    } else {
      ++rejected;
    }
    hash = hash_u64_value(hash, current_words);
    hash = hash_u64_value(
        hash, (uint64_t)(unsigned int)step.step.accepted);
    hash = hash_u64_value(hash, rng.state);
  }
  if (mvmc_bounded_krylov_session_end(workspace) !=
      MVMC_KRYLOV_STATUS_OK) {
    return 0;
  }
  memset(trace, 0, sizeof(*trace));
  trace->trace_hash = hash;
  trace->final_words = current_words;
  trace->accepted = accepted;
  trace->rejected = rejected;
  trace->rng = rng;
  trace->current = current;
  return 1;
}

static void test_accept_reject_body(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude) {
  const MVMCKrylovPositiveGuidePolicy policy = guide_policy();
  MVMCKrylovPositiveSamplerSnapshot current;
  MVMCKrylovPositiveSamplerStepResult step;
  uint64_t current_words[1] = {UINT64_C(1)};
  uint64_t high_words[1] = {UINT64_C(2)};
  uint64_t low_words[1] = {UINT64_C(1)};
  const double low_log_guide = log(0.25 + 1.0 + 4.0);
  const double high_log_guide = log(0.25 + 16.0 + 4.0);

  CHECK(mvmc_krylov_positive_sampler_initialize(
            workspace, &policy, current_words, 1, table_callback,
            amplitude, &current) == MVMC_KRYLOV_STATUS_OK &&
            current.valid,
        "sampler initialize current");
  CHECK(close_double(current.guide.log_guide, low_log_guide, 1.0e-14) &&
            current.accepted_generation == 0,
        "initial current guide/generation");

  CHECK(mvmc_krylov_positive_sampler_step(
            workspace, &policy, current_words, 1, &current, high_words, 1,
            0.0, 0.999, table_callback, amplitude, &step) ==
            MVMC_KRYLOV_STATUS_OK &&
            step.valid && step.accepted && current_words[0] == high_words[0],
        "higher-guide proposal commits");
  CHECK(current.accepted_generation == 1 && step.accepted_generation == 1 &&
            close_double(current.guide.log_guide, high_log_guide, 1.0e-14),
        "accepted proposal updates current snapshot");

  CHECK(mvmc_krylov_positive_sampler_step(
            workspace, &policy, current_words, 1, &current, low_words, 1,
            0.0, 0.999, table_callback, amplitude, &step) ==
            MVMC_KRYLOV_STATUS_OK &&
            step.valid && !step.accepted &&
            current_words[0] == high_words[0],
        "lower-guide proposal rejects without changing current words");
  CHECK(current.accepted_generation == 1 && step.accepted_generation == 1 &&
            close_double(current.guide.log_guide, high_log_guide, 1.0e-14) &&
            close_double(step.proposal_guide.log_guide, low_log_guide,
                         1.0e-14),
        "rejected proposal leaves current snapshot intact");
}

static void test_persistent_session_body(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude) {
  const MVMCKrylovPositiveGuidePolicy policy = guide_policy();
  const uint64_t generation = UINT64_C(0x5034533953414d50);
  uint64_t legacy_words[1] = {UINT64_C(1)};
  uint64_t session_words[1] = {UINT64_C(1)};
  const uint64_t proposal_words[1] = {UINT64_C(2)};
  MVMCKrylovPositiveSamplerSnapshot legacy_current;
  MVMCKrylovPositiveSamplerSnapshot session_current;
  MVMCKrylovPositiveSamplerSnapshot invalid_snapshot;
  MVMCKrylovPositiveSamplerStepResult legacy_step;
  MVMCKrylovPositiveSamplerStepResult session_step;
  TableAmplitude substitute = *amplitude;
  SessionTrace uninterrupted;
  SessionTrace split;
  int order;

  CHECK(mvmc_krylov_positive_sampler_initialize(
            workspace, &policy, legacy_words, 1, table_callback, amplitude,
            &legacy_current) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_step(
                workspace, &policy, legacy_words, 1, &legacy_current,
                proposal_words, 1, 0.0, 0.0, table_callback, amplitude,
                &legacy_step) == MVMC_KRYLOV_STATUS_OK,
        "legacy reference for persistent sampler failed");

  amplitude->calls = 0;
  CHECK(mvmc_bounded_krylov_session_begin(
            workspace, table_callback, amplitude, generation) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_initialize(
                workspace, &policy, session_words, 1, table_callback,
                amplitude, &session_current) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_step(
                workspace, &policy, session_words, 1, &session_current,
                proposal_words, 1, 0.0, 0.0, table_callback, amplitude,
                &session_step) == MVMC_KRYLOV_STATUS_OK,
        "persistent sampler initialize/step failed");
  CHECK(session_current.krylov.statistics.persistent_session_active == 1 &&
            session_current.krylov.statistics.session_root_evaluation == 2 &&
            session_step.proposal_krylov.statistics
                    .persistent_session_active == 1 &&
            session_step.proposal_krylov.statistics
                    .session_root_evaluation == 2 &&
            session_step.proposal_krylov.statistics.cache_reset_performed ==
                0 &&
            session_step.proposal_krylov.statistics
                    .cache_entries_resident_start > 0,
        "persistent sampler did not use the active engine session");
  CHECK(session_step.accepted == legacy_step.accepted &&
            session_words[0] == legacy_words[0] &&
            session_current.guide.log_guide ==
                legacy_current.guide.log_guide &&
            session_step.proposal_guide.log_guide ==
                legacy_step.proposal_guide.log_guide,
        "persistent sampler changed guide/acceptance bytes");
  for (order = 0; order <= policy.order; ++order) {
    CHECK(memcmp(&session_step.proposal_krylov.value[order],
                 &legacy_step.proposal_krylov.value[order],
                 sizeof(session_step.proposal_krylov.value[order])) == 0,
          "persistent sampler changed proposal value[%d] bytes", order);
  }
  CHECK(mvmc_bounded_krylov_session_end(workspace) ==
            MVMC_KRYLOV_STATUS_OK,
        "persistent sampler session end failed");

  CHECK(mvmc_bounded_krylov_session_begin(
            workspace, table_callback, amplitude, generation + UINT64_C(1)) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_initialize(
                workspace, &policy, session_words, 1, table_callback,
                &substitute, &invalid_snapshot) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !mvmc_bounded_krylov_session_is_active(workspace) &&
            !invalid_snapshot.valid,
        "sampler callback-context substitution did not fail closed");

  CHECK(run_session_trace(
            workspace, amplitude, &policy, generation + UINT64_C(2), 64,
            64, &uninterrupted) &&
            run_session_trace(
                workspace, amplitude, &policy, generation + UINT64_C(2),
                64, 23, &split),
        "persistent sampler uninterrupted/split trace failed");
  CHECK(uninterrupted.trace_hash == split.trace_hash &&
            uninterrupted.final_words == split.final_words &&
            uninterrupted.accepted == split.accepted &&
            uninterrupted.rejected == split.rejected &&
            memcmp(&uninterrupted.rng, &split.rng,
                   sizeof(uninterrupted.rng)) == 0 &&
            uninterrupted.current.accepted_generation ==
                split.current.accepted_generation &&
            uninterrupted.current.guide.log_guide ==
                split.current.guide.log_guide,
        "cold split restart changed session statistical bytes");
}

static void test_failure_transaction_body(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude) {
  const MVMCKrylovPositiveGuidePolicy policy = guide_policy();
  MVMCKrylovPositiveSamplerSnapshot current;
  MVMCKrylovPositiveSamplerStepResult step;
  uint64_t current_words[1] = {UINT64_C(2)};
  uint64_t proposal_words[1] = {UINT64_C(1)};
  const double current_log_guide = log(0.25 + 16.0 + 4.0);

  CHECK(mvmc_krylov_positive_sampler_initialize(
            workspace, &policy, current_words, 1, table_callback,
            amplitude, &current) == MVMC_KRYLOV_STATUS_OK &&
            current.valid,
        "failure test initialize current");
  amplitude->forced_status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  CHECK(mvmc_krylov_positive_sampler_step(
            workspace, &policy, current_words, 1, &current, proposal_words,
            1, 0.0, 0.0, table_callback, amplitude, &step) ==
            MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
            !step.valid && current_words[0] == UINT64_C(2),
        "proposal evaluation failure does not commit");
  CHECK(current.accepted_generation == 0 &&
            close_double(current.guide.log_guide, current_log_guide,
                         1.0e-14),
        "failure leaves current snapshot intact");
}

static void test_policy_mismatch_body(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude) {
  MVMCKrylovPositiveGuidePolicy policy = guide_policy();
  MVMCKrylovPositiveSamplerSnapshot current;
  MVMCKrylovPositiveSamplerStepResult step;
  uint64_t current_words[1] = {UINT64_C(1)};
  uint64_t proposal_words[1] = {UINT64_C(2)};

  CHECK(mvmc_krylov_positive_sampler_initialize(
            workspace, &policy, current_words, 1, table_callback,
            amplitude, &current) == MVMC_KRYLOV_STATUS_OK &&
            current.valid,
        "policy mismatch initialize current");
  policy.policy_hash ^= UINT64_C(1);
  CHECK(mvmc_krylov_positive_sampler_step(
            workspace, &policy, current_words, 1, &current, proposal_words,
            1, 0.0, 0.5, table_callback, amplitude, &step) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !step.valid && current_words[0] == UINT64_C(1),
        "policy hash mismatch is rejected transactionally");
}

static void test_selected_neighbor_trace_body(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude) {
  const MVMCKrylovPositiveGuidePolicy policy = guide_policy();
  const MVMCKrylovFockModel proposal = proposal_model();
  MVMCKrylovPositiveSamplerSnapshot current;
  MVMCKrylovPositiveSamplerProposalStepResult step;
  uint64_t current_words[1] = {UINT64_C(1)};
  uint64_t proposal_words[1] = {0};

  CHECK(mvmc_krylov_positive_sampler_initialize(
            workspace, &policy, current_words, 1, table_callback,
            amplitude, &current) == MVMC_KRYLOV_STATUS_OK &&
            current.valid,
        "selected-neighbor trace initialize current");
  CHECK(mvmc_krylov_positive_sampler_step_selected_neighbor(
            workspace, &policy, &proposal, current_words, 1, &current,
            0, 0.999, table_callback, amplitude, proposal_words, 1,
            &step) == MVMC_KRYLOV_STATUS_OK &&
            step.valid && step.neighbor_count == 1 &&
            close_double(step.log_proposal_ratio, 0.0, 1.0e-14) &&
            step.step.accepted && current_words[0] == UINT64_C(2) &&
            proposal_words[0] == UINT64_C(2),
        "selected-neighbor accepted trace step");
  CHECK(current.accepted_generation == 1 &&
            step.step.accepted_generation == 1,
        "selected-neighbor accepted generation");

  proposal_words[0] = 0;
  CHECK(mvmc_krylov_positive_sampler_step_selected_neighbor(
            workspace, &policy, &proposal, current_words, 1, &current,
            0, 0.999, table_callback, amplitude, proposal_words, 1,
            &step) == MVMC_KRYLOV_STATUS_OK &&
            step.valid && !step.step.accepted &&
            current_words[0] == UINT64_C(2) &&
            proposal_words[0] == UINT64_C(1),
        "selected-neighbor rejected trace step");
  CHECK(current.accepted_generation == 1 &&
            step.step.accepted_generation == 1,
        "selected-neighbor reject keeps generation");
}

static uint64_t run_rng_steps(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *policy,
    const MVMCKrylovFockModel *proposal,
    uint64_t *words,
    MVMCKrylovPositiveSamplerSnapshot *current,
    MVMCKrylovPositiveSamplerRng *rng,
    TableAmplitude *amplitude,
    size_t step_count,
    MVMCKrylovPositiveSamplerTraceStatistics *statistics) {
  uint64_t trace_hash = UINT64_C(1469598103934665603);
  size_t step_index;
  for (step_index = 0; step_index < step_count; ++step_index) {
    uint64_t proposal_words[1] = {0};
    MVMCKrylovPositiveSamplerProposalStepResult step;
    const MVMCKrylovStatus status =
        mvmc_krylov_positive_sampler_step_rng(
            workspace, policy, proposal, words, 1, current, rng,
            table_callback, amplitude, proposal_words, 1, &step);
    CHECK(status == MVMC_KRYLOV_STATUS_OK && step.valid,
          "rng step %zu failed: %s", step_index,
          mvmc_krylov_status_string(status));
    if (status != MVMC_KRYLOV_STATUS_OK || !step.valid) break;
    CHECK(step.rng_after.draws == rng->draws &&
              step.rng_after.state == rng->state,
          "rng step %zu did not publish committed rng", step_index);
    if (statistics != NULL) {
      CHECK(mvmc_krylov_positive_sampler_trace_statistics_record_step(
                &step, words, 1, statistics) == MVMC_KRYLOV_STATUS_OK,
            "trace statistics record step %zu", step_index);
    }
    trace_hash = hash_u64_value(trace_hash, words[0]);
    trace_hash = hash_u64_value(
        trace_hash, (uint64_t)step.step.accepted);
    trace_hash = hash_u64_value(trace_hash, rng->draws);
  }
  return trace_hash;
}

static uint64_t run_mixture_steps(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveGuidePolicy *policy,
    const MVMCKrylovFockModel *proposal,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    uint64_t *words, MVMCKrylovPositiveSamplerSnapshot *current,
    MVMCKrylovPositiveSamplerRng *rng, TableAmplitude *amplitude,
    size_t step_count,
    MVMCKrylovPositiveSamplerTraceStatistics *statistics) {
  uint64_t trace_hash = UINT64_C(1469598103934665603);
  size_t step_index;
  for (step_index = 0; step_index < step_count; ++step_index) {
    uint64_t proposal_words[1] = {0};
    MVMCKrylovPositiveSamplerProposalStepResult step;
    const MVMCKrylovStatus status =
        mvmc_krylov_positive_sampler_step_mixture_rng(
            workspace, policy, proposal, proposal_policy, words, 1,
            current, rng, table_callback, amplitude, proposal_words, 1,
            &step);
    CHECK(status == MVMC_KRYLOV_STATUS_OK && step.valid,
          "mixture step %zu failed: %s", step_index,
          mvmc_krylov_status_string(status));
    if (status != MVMC_KRYLOV_STATUS_OK || !step.valid) break;
    CHECK(step.rng_after.draws == rng->draws &&
              step.rng_after.state == rng->state &&
              step.proposal_policy_hash == proposal_policy->policy_hash &&
              step.proposal_model_hash != 0,
          "mixture step %zu metadata mismatch", step_index);
    if (statistics != NULL) {
      CHECK(mvmc_krylov_positive_sampler_trace_statistics_record_step(
                &step, words, 1, statistics) == MVMC_KRYLOV_STATUS_OK,
            "mixture statistics record step %zu", step_index);
    }
    trace_hash = hash_u64_value(trace_hash, words[0]);
    trace_hash = hash_u64_value(trace_hash, (uint64_t)step.component);
    trace_hash = hash_u64_value(trace_hash, (uint64_t)step.self_proposal);
    trace_hash = hash_u64_value(trace_hash, rng->draws);
  }
  return trace_hash;
}

static int trace_statistics_equal(
    const MVMCKrylovPositiveSamplerTraceStatistics *left,
    const MVMCKrylovPositiveSamplerTraceStatistics *right) {
  return left->valid == right->valid &&
         left->status == right->status &&
         left->attempted_steps == right->attempted_steps &&
         left->completed_steps == right->completed_steps &&
         left->accepted_steps == right->accepted_steps &&
         left->rejected_steps == right->rejected_steps &&
         left->neighbor_attempted_steps == right->neighbor_attempted_steps &&
         left->neighbor_accepted_steps == right->neighbor_accepted_steps &&
         left->neighbor_rejected_steps == right->neighbor_rejected_steps &&
         left->global_attempted_steps == right->global_attempted_steps &&
         left->global_accepted_steps == right->global_accepted_steps &&
         left->global_rejected_steps == right->global_rejected_steps &&
         left->global_self_proposals == right->global_self_proposals &&
         left->shell_attempted_steps == right->shell_attempted_steps &&
         left->shell_accepted_steps == right->shell_accepted_steps &&
         left->shell_rejected_steps == right->shell_rejected_steps &&
         left->configuration_changing_accepted_moves ==
             right->configuration_changing_accepted_moves &&
         left->positive_support_steps == right->positive_support_steps &&
         left->support_violation_steps == right->support_violation_steps &&
         left->floor_supported_zero_steps ==
             right->floor_supported_zero_steps &&
         left->finite_proposal_components ==
             right->finite_proposal_components &&
         left->exact_zero_proposal_components ==
             right->exact_zero_proposal_components &&
         left->numeric_zero_proposal_components ==
             right->numeric_zero_proposal_components &&
         left->terminal_amplitude_calls == right->terminal_amplitude_calls &&
         left->trace_hash == right->trace_hash &&
         close_double(left->min_log_acceptance_ratio,
                      right->min_log_acceptance_ratio, 0.0) &&
         close_double(left->max_log_acceptance_ratio,
                      right->max_log_acceptance_ratio, 0.0) &&
         close_double(left->sum_log_acceptance_ratio,
                      right->sum_log_acceptance_ratio, 0.0) &&
         close_double(left->min_proposal_log_guide,
                      right->min_proposal_log_guide, 0.0) &&
         close_double(left->max_proposal_log_guide,
                      right->max_proposal_log_guide, 0.0);
}

static void test_rng_manifest_restart_body(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude,
    uint64_t plan_hash, const MVMCKrylovBoundedLimits *limits) {
  const MVMCKrylovPositiveGuidePolicy policy = guide_policy();
  MVMCKrylovPositiveGuidePolicy alternate_guide_policy = guide_policy();
  const MVMCKrylovFockModel proposal = proposal_model();
  MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
  MVMCKrylovPositiveSamplerProposalPolicy alternate_proposal_policy;
  MVMCKrylovPositiveSamplerSurrogatePolicy surrogate_k4_alpha1;
  MVMCKrylovPositiveSamplerSurrogatePolicy surrogate_k16_alpha1;
  MVMCKrylovPositiveSamplerSurrogatePolicy surrogate_k4_alpha10;
  MVMCKrylovFockModel alternate_proposal = proposal_model();
  MVMCKrylovPositiveSamplerSnapshot full_current;
  MVMCKrylovPositiveSamplerSnapshot split_current;
  MVMCKrylovPositiveSamplerSnapshot resumed_current;
  MVMCKrylovPositiveSamplerSnapshot alternate_guide_current;
  MVMCKrylovPositiveSamplerRng full_rng;
  MVMCKrylovPositiveSamplerRng split_rng;
  MVMCKrylovPositiveSamplerRng resumed_rng;
  MVMCKrylovPositiveSamplerManifest split_manifest;
  MVMCKrylovPositiveSamplerManifest full_manifest;
  MVMCKrylovPositiveSamplerManifest resumed_manifest;
  MVMCKrylovPositiveSamplerManifest alternate_guide_manifest;
  MVMCKrylovPositiveSamplerManifest surrogate_manifest;
  uint64_t full_words[1] = {UINT64_C(1)};
  uint64_t split_words[1] = {UINT64_C(1)};
  uint64_t resumed_words[1] = {0};
  uint64_t full_hash;
  uint64_t split_hash;
  uint64_t resumed_hash;
  int matches = 0;

  CHECK(mvmc_krylov_positive_sampler_proposal_policy_create(
            0, 1, &proposal_policy) == MVMC_KRYLOV_STATUS_OK,
        "restart proposal policy create");
  CHECK(mvmc_krylov_positive_sampler_proposal_policy_create(
            1, 2, &alternate_proposal_policy) == MVMC_KRYLOV_STATUS_OK,
        "alternate proposal policy create");
  alternate_proposal.site_count = 3;

  CHECK(mvmc_krylov_positive_sampler_rng_seed(
            UINT64_C(0x123456789abcdef0), UINT64_C(7), &full_rng) ==
            MVMC_KRYLOV_STATUS_OK &&
            full_rng.valid,
        "full rng seed");
  CHECK(mvmc_krylov_positive_sampler_rng_seed(
            UINT64_C(0x123456789abcdef0), UINT64_C(7), &split_rng) ==
            MVMC_KRYLOV_STATUS_OK &&
            split_rng.valid,
        "split rng seed");
  CHECK(mvmc_krylov_positive_sampler_initialize(
            workspace, &policy, full_words, 1, table_callback,
            amplitude, &full_current) == MVMC_KRYLOV_STATUS_OK &&
            full_current.valid,
        "full rng initialize");
  CHECK(mvmc_krylov_positive_sampler_initialize(
            workspace, &policy, split_words, 1, table_callback,
            amplitude, &split_current) == MVMC_KRYLOV_STATUS_OK &&
            split_current.valid,
        "split rng initialize");

  full_hash = run_rng_steps(workspace, &policy, &proposal, full_words,
                            &full_current, &full_rng, amplitude, 4, NULL);
  split_hash = run_rng_steps(workspace, &policy, &proposal, split_words,
                             &split_current, &split_rng, amplitude, 2, NULL);
  CHECK(mvmc_krylov_positive_sampler_manifest_create(
            &policy, &proposal_policy, &proposal, limits, plan_hash,
            split_words, 1, &split_current, &split_rng,
            &split_manifest) == MVMC_KRYLOV_STATUS_OK &&
            split_manifest.valid,
        "split manifest create");

  resumed_words[0] = split_words[0];
  resumed_current = split_current;
  resumed_rng = split_rng;
  CHECK(mvmc_krylov_positive_sampler_manifest_matches(
            &split_manifest, &policy, &proposal_policy, &proposal, limits,
            plan_hash, resumed_words, 1, &resumed_current, &resumed_rng,
            &matches) ==
            MVMC_KRYLOV_STATUS_OK &&
            matches,
        "restored manifest does not match");
  CHECK(mvmc_krylov_positive_sampler_surrogate_policy_create(
            4, 1, &policy, &surrogate_k4_alpha1) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_surrogate_policy_create(
                16, 1, &policy, &surrogate_k16_alpha1) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_surrogate_policy_create(
                4, 10, &policy, &surrogate_k4_alpha10) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_surrogate_manifest_create(
                &policy, &proposal_policy, &surrogate_k4_alpha1, &proposal,
                limits, plan_hash, resumed_words, 1, &resumed_current,
                &resumed_rng, &surrogate_manifest) ==
                MVMC_KRYLOV_STATUS_OK &&
            surrogate_manifest.valid &&
            surrogate_manifest.surrogate_policy_hash ==
                surrogate_k4_alpha1.policy_hash,
        "surrogate manifest create");
  CHECK(mvmc_krylov_positive_sampler_surrogate_manifest_matches(
            &surrogate_manifest, &policy, &proposal_policy,
            &surrogate_k4_alpha1, &proposal, limits, plan_hash,
            resumed_words, 1, &resumed_current, &resumed_rng, &matches) ==
            MVMC_KRYLOV_STATUS_OK &&
            matches,
        "surrogate manifest exact match");
  CHECK(mvmc_krylov_positive_sampler_surrogate_manifest_matches(
            &surrogate_manifest, &policy, &proposal_policy,
            &surrogate_k16_alpha1, &proposal, limits, plan_hash,
            resumed_words, 1, &resumed_current, &resumed_rng, &matches) ==
            MVMC_KRYLOV_STATUS_OK &&
            !matches,
        "surrogate cross-K mismatch");
  CHECK(mvmc_krylov_positive_sampler_surrogate_manifest_matches(
            &surrogate_manifest, &policy, &proposal_policy,
            &surrogate_k4_alpha10, &proposal, limits, plan_hash,
            resumed_words, 1, &resumed_current, &resumed_rng, &matches) ==
            MVMC_KRYLOV_STATUS_OK &&
            !matches,
        "surrogate cross-alpha mismatch");
  alternate_guide_policy.eta = 1.0;
  alternate_guide_policy.policy_hash ^= UINT64_C(0x503453354752484f);
  alternate_guide_current = split_current;
  alternate_guide_current.policy_hash = alternate_guide_policy.policy_hash;
  CHECK(mvmc_krylov_positive_guide_evaluate(
            &alternate_guide_policy, alternate_guide_current.krylov.value,
            (size_t)alternate_guide_current.krylov.evaluated_order + 1,
            &alternate_guide_current.guide) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_manifest_create(
                &alternate_guide_policy, &proposal_policy, &proposal, limits,
                plan_hash, split_words, 1, &alternate_guide_current,
                &split_rng, &alternate_guide_manifest) ==
                MVMC_KRYLOV_STATUS_OK &&
            alternate_guide_manifest.valid &&
            alternate_guide_manifest.guide_shape_hash !=
                split_manifest.guide_shape_hash &&
            alternate_guide_manifest.proposal_policy_hash ==
                split_manifest.proposal_policy_hash &&
            alternate_guide_manifest.proposal_model_hash ==
                split_manifest.proposal_model_hash,
        "guide-rho manifest identity");
  CHECK(mvmc_krylov_positive_sampler_manifest_matches(
            &split_manifest, &alternate_guide_policy, &proposal_policy,
            &proposal, limits, plan_hash, split_words, 1,
            &alternate_guide_current, &split_rng, &matches) ==
            MVMC_KRYLOV_STATUS_OK &&
            !matches,
        "guide-rho mismatch not detected");
  CHECK(mvmc_krylov_positive_sampler_manifest_matches(
            &split_manifest, &policy, &alternate_proposal_policy,
            &proposal, limits, plan_hash, resumed_words, 1,
            &resumed_current, &resumed_rng, &matches) ==
            MVMC_KRYLOV_STATUS_OK &&
            !matches,
        "proposal-policy mismatch not detected");
  CHECK(mvmc_krylov_positive_sampler_manifest_matches(
            &split_manifest, &policy, &proposal_policy,
            &alternate_proposal, limits, plan_hash, resumed_words, 1,
            &resumed_current, &resumed_rng, &matches) ==
            MVMC_KRYLOV_STATUS_OK &&
            !matches,
        "proposal-model mismatch not detected");
  resumed_rng.draws += 2;
  CHECK(mvmc_krylov_positive_sampler_manifest_matches(
            &split_manifest, &policy, &proposal_policy, &proposal, limits,
            plan_hash, resumed_words, 1, &resumed_current, &resumed_rng,
            &matches) ==
            MVMC_KRYLOV_STATUS_OK &&
            !matches,
        "rng mismatch not detected");
  resumed_rng = split_rng;

  resumed_hash = run_rng_steps(workspace, &policy, &proposal, resumed_words,
                               &resumed_current, &resumed_rng, amplitude, 2,
                               NULL);
  CHECK(full_words[0] == resumed_words[0] &&
            full_current.accepted_generation ==
                resumed_current.accepted_generation &&
            full_rng.state == resumed_rng.state &&
            full_rng.draws == resumed_rng.draws &&
            close_double(full_current.guide.log_guide,
                         resumed_current.guide.log_guide, 1.0e-14),
        "full and resumed chains diverged");
  CHECK(mvmc_krylov_positive_sampler_manifest_create(
            &policy, &proposal_policy, &proposal, limits, plan_hash,
            full_words, 1, &full_current, &full_rng,
            &full_manifest) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_manifest_create(
                &policy, &proposal_policy, &proposal, limits, plan_hash,
                resumed_words, 1, &resumed_current, &resumed_rng,
                &resumed_manifest) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_manifest_matches(
                &full_manifest, &policy, &proposal_policy, &proposal,
                limits, plan_hash, resumed_words, 1, &resumed_current,
                &resumed_rng, &matches) ==
                MVMC_KRYLOV_STATUS_OK &&
            matches &&
            full_manifest.current_scale_hash ==
                resumed_manifest.current_scale_hash,
        "full/resumed manifest mismatch");
  CHECK(full_hash != 0 && split_hash != 0 && resumed_hash != 0,
        "trace hashes were not updated");
}

static void test_full_markov_trace_statistics_body(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude) {
  const MVMCKrylovPositiveGuidePolicy policy = guide_policy();
  const MVMCKrylovFockModel proposal = proposal_model();
  MVMCKrylovPositiveSamplerSnapshot current;
  MVMCKrylovPositiveSamplerSnapshot repeat_current;
  MVMCKrylovPositiveSamplerRng rng;
  MVMCKrylovPositiveSamplerRng repeat_rng;
  MVMCKrylovPositiveSamplerTraceStatistics statistics;
  MVMCKrylovPositiveSamplerTraceStatistics repeat_statistics;
  uint64_t words[1] = {UINT64_C(1)};
  uint64_t repeat_words[1] = {UINT64_C(1)};
  const size_t step_count = 16;

  CHECK(mvmc_krylov_positive_sampler_rng_seed(
            UINT64_C(0xf00dcafebeef1234), UINT64_C(3), &rng) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_rng_seed(
                UINT64_C(0xf00dcafebeef1234), UINT64_C(3),
                &repeat_rng) == MVMC_KRYLOV_STATUS_OK,
        "trace rng seed");
  CHECK(mvmc_krylov_positive_sampler_trace_statistics_reset(
            &statistics) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_trace_statistics_reset(
                &repeat_statistics) == MVMC_KRYLOV_STATUS_OK,
        "trace statistics reset");
  CHECK(mvmc_krylov_positive_sampler_initialize(
            workspace, &policy, words, 1, table_callback,
            amplitude, &current) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_initialize(
                workspace, &policy, repeat_words, 1, table_callback,
                amplitude, &repeat_current) == MVMC_KRYLOV_STATUS_OK,
        "trace initialize");

  (void)run_rng_steps(workspace, &policy, &proposal, words, &current,
                      &rng, amplitude, step_count, &statistics);
  (void)run_rng_steps(workspace, &policy, &proposal, repeat_words,
                      &repeat_current, &repeat_rng, amplitude,
                      step_count, &repeat_statistics);

  CHECK(statistics.completed_steps == (uint64_t)step_count &&
            statistics.attempted_steps == statistics.completed_steps &&
            statistics.accepted_steps + statistics.rejected_steps ==
                statistics.completed_steps &&
            statistics.accepted_steps > 0 &&
            statistics.rejected_steps > 0,
        "acceptance evidence counts");
  CHECK(statistics.positive_support_steps == statistics.completed_steps &&
            statistics.support_violation_steps == 0 &&
            statistics.floor_supported_zero_steps == 0 &&
            statistics.exact_zero_proposal_components == 0 &&
            statistics.numeric_zero_proposal_components == 0 &&
            statistics.finite_proposal_components ==
                2 * statistics.completed_steps,
        "nonzero-support trace evidence counts");
  CHECK(statistics.min_log_acceptance_ratio < 0.0 &&
            statistics.max_log_acceptance_ratio > 0.0 &&
            statistics.min_proposal_log_guide <=
                statistics.max_proposal_log_guide &&
            statistics.terminal_amplitude_calls > 0 &&
            statistics.trace_hash != UINT64_C(1469598103934665603),
        "trace statistic ranges/hash");
  CHECK(trace_statistics_equal(&statistics, &repeat_statistics),
        "repeat trace statistics mismatch");
}

static void test_floor_supported_zero_evidence_body(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude) {
  const MVMCKrylovPositiveGuidePolicy policy = guide_policy();
  const MVMCKrylovFockModel proposal = proposal_model();
  MVMCKrylovPositiveSamplerSnapshot current;
  MVMCKrylovPositiveSamplerProposalStepResult step;
  MVMCKrylovPositiveSamplerTraceStatistics statistics;
  uint64_t current_words[1] = {UINT64_C(2)};
  uint64_t proposal_words[1] = {0};

  amplitude->values[1] = 0.0;
  amplitude->values[2] = 4.0;
  CHECK(mvmc_krylov_positive_sampler_trace_statistics_reset(
            &statistics) == MVMC_KRYLOV_STATUS_OK,
        "zero-support statistics reset");
  CHECK(mvmc_krylov_positive_sampler_initialize(
            workspace, &policy, current_words, 1, table_callback,
            amplitude, &current) == MVMC_KRYLOV_STATUS_OK &&
            current.valid,
        "zero-support initialize");
  CHECK(mvmc_krylov_positive_sampler_step_selected_neighbor(
            workspace, &policy, &proposal, current_words, 1, &current,
            0, 0.0, table_callback, amplitude, proposal_words, 1,
            &step) == MVMC_KRYLOV_STATUS_OK &&
            step.valid,
        "zero-support selected step");
  CHECK(mvmc_krylov_positive_sampler_trace_statistics_record_step(
            &step, current_words, 1, &statistics) ==
            MVMC_KRYLOV_STATUS_OK,
        "zero-support statistics record");
  CHECK(statistics.completed_steps == 1 &&
            statistics.accepted_steps == 1 &&
            statistics.positive_support_steps == 1 &&
            statistics.support_violation_steps == 0 &&
            statistics.floor_supported_zero_steps == 1 &&
            statistics.exact_zero_proposal_components > 0 &&
            statistics.numeric_zero_proposal_components == 0,
        "floor-supported zero evidence counts");
}

static void test_zero_fraction_preserves_legacy_body(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude) {
  const MVMCKrylovPositiveGuidePolicy policy = guide_policy();
  const MVMCKrylovFockModel proposal = proposal_model();
  MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
  MVMCKrylovPositiveSamplerSnapshot legacy_current;
  MVMCKrylovPositiveSamplerSnapshot mixture_current;
  MVMCKrylovPositiveSamplerRng legacy_rng;
  MVMCKrylovPositiveSamplerRng mixture_rng;
  uint64_t legacy_words[1] = {UINT64_C(1)};
  uint64_t mixture_words[1] = {UINT64_C(1)};
  size_t step_index;

  CHECK(mvmc_krylov_positive_sampler_proposal_policy_create(
            0, 1, &proposal_policy) == MVMC_KRYLOV_STATUS_OK,
        "zero-fraction policy create");
  CHECK(mvmc_krylov_positive_sampler_rng_seed(
            UINT64_C(0x503453335a45524f), UINT64_C(11), &legacy_rng) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_rng_seed(
                UINT64_C(0x503453335a45524f), UINT64_C(11),
                &mixture_rng) == MVMC_KRYLOV_STATUS_OK,
        "zero-fraction rng seed");
  CHECK(mvmc_krylov_positive_sampler_initialize(
            workspace, &policy, legacy_words, 1, table_callback,
            amplitude, &legacy_current) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_initialize(
                workspace, &policy, mixture_words, 1, table_callback,
                amplitude, &mixture_current) == MVMC_KRYLOV_STATUS_OK,
        "zero-fraction initialize");

  for (step_index = 0; step_index < 32; ++step_index) {
    MVMCKrylovPositiveSamplerProposalStepResult legacy_step;
    MVMCKrylovPositiveSamplerProposalStepResult mixture_step;
    uint64_t legacy_proposal[1] = {0};
    uint64_t mixture_proposal[1] = {0};
    const MVMCKrylovStatus legacy_status =
        mvmc_krylov_positive_sampler_step_rng(
            workspace, &policy, &proposal, legacy_words, 1,
            &legacy_current, &legacy_rng, table_callback, amplitude,
            legacy_proposal, 1, &legacy_step);
    const MVMCKrylovStatus mixture_status =
        mvmc_krylov_positive_sampler_step_mixture_rng(
            workspace, &policy, &proposal, &proposal_policy,
            mixture_words, 1, &mixture_current, &mixture_rng,
            table_callback, amplitude, mixture_proposal, 1, &mixture_step);
    CHECK(legacy_status == MVMC_KRYLOV_STATUS_OK &&
              mixture_status == MVMC_KRYLOV_STATUS_OK &&
              legacy_step.valid && mixture_step.valid,
          "zero-fraction step %zu status", step_index);
    CHECK(legacy_words[0] == mixture_words[0] &&
              legacy_proposal[0] == mixture_proposal[0] &&
              legacy_rng.state == mixture_rng.state &&
              legacy_rng.draws == mixture_rng.draws &&
              legacy_step.selected_neighbor_index ==
                  mixture_step.selected_neighbor_index &&
              legacy_step.neighbor_count == mixture_step.neighbor_count &&
              legacy_step.uniform == mixture_step.uniform &&
              legacy_step.step.accepted == mixture_step.step.accepted &&
              legacy_current.accepted_generation ==
                  mixture_current.accepted_generation &&
              legacy_current.guide.log_guide ==
                  mixture_current.guide.log_guide,
          "zero-fraction step %zu changed legacy sequence", step_index);
  }
}

static void test_mixture_counters_body(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude) {
  const MVMCKrylovPositiveGuidePolicy policy = guide_policy();
  const MVMCKrylovFockModel proposal = proposal_model();
  MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
  MVMCKrylovPositiveSamplerSnapshot current;
  MVMCKrylovPositiveSamplerRng rng;
  MVMCKrylovPositiveSamplerTraceStatistics statistics;
  uint64_t words[1] = {UINT64_C(1)};

  CHECK(mvmc_krylov_positive_sampler_proposal_policy_create(
            1, 2, &proposal_policy) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_rng_seed(
                UINT64_C(0x503453334d495831), UINT64_C(5), &rng) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_trace_statistics_reset(
                &statistics) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_initialize(
                workspace, &policy, words, 1, table_callback, amplitude,
                &current) == MVMC_KRYLOV_STATUS_OK,
        "mixture counter setup");
  (void)run_mixture_steps(workspace, &policy, &proposal,
                          &proposal_policy, words, &current, &rng,
                          amplitude, 128, &statistics);
  CHECK(statistics.completed_steps == 128 &&
            statistics.neighbor_attempted_steps +
                    statistics.global_attempted_steps ==
                statistics.completed_steps &&
            statistics.neighbor_attempted_steps > 0 &&
            statistics.global_attempted_steps > 0 &&
            statistics.neighbor_accepted_steps +
                    statistics.neighbor_rejected_steps ==
                statistics.neighbor_attempted_steps &&
            statistics.global_accepted_steps +
                    statistics.global_rejected_steps ==
                statistics.global_attempted_steps &&
            statistics.global_self_proposals > 0 &&
            statistics.configuration_changing_accepted_moves <=
                statistics.accepted_steps,
        "mixture component/self/mobility counters");
}

static void test_mixture_degenerate_and_failure_body(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude) {
  const MVMCKrylovPositiveGuidePolicy policy = guide_policy();
  MVMCKrylovFockModel one_state = proposal_model();
  const MVMCKrylovFockModel regular = proposal_model();
  MVMCKrylovPositiveSamplerProposalPolicy global_only;
  MVMCKrylovPositiveSamplerProposalPolicy mixed;
  MVMCKrylovPositiveSamplerSnapshot current;
  MVMCKrylovPositiveSamplerProposalStepResult step;
  MVMCKrylovPositiveSamplerRng rng;
  MVMCKrylovPositiveSamplerRng saved_rng;
  uint64_t words[1] = {UINT64_C(1)};
  uint64_t proposal_words[1] = {UINT64_MAX};

  one_state.site_count = 1;
  one_state.up_electron_count = 1;
  CHECK(mvmc_krylov_positive_sampler_proposal_policy_create(
            1, 1, &global_only) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_proposal_policy_create(
                1, 2, &mixed) == MVMC_KRYLOV_STATUS_OK,
        "degenerate policies create");
  CHECK(mvmc_krylov_positive_sampler_rng_seed(
            UINT64_C(0x5034533344454731), UINT64_C(2), &rng) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_initialize(
                workspace, &policy, words, 1, table_callback, amplitude,
                &current) == MVMC_KRYLOV_STATUS_OK,
        "degenerate setup");
  CHECK(mvmc_krylov_positive_sampler_step_mixture_rng(
            workspace, &policy, &one_state, &global_only, words, 1,
            &current, &rng, table_callback, amplitude, proposal_words, 1,
            &step) == MVMC_KRYLOV_STATUS_OK &&
            step.valid &&
            step.component == MVMC_KRYLOV_PROPOSAL_COMPONENT_GLOBAL &&
            step.self_proposal && !step.configuration_changed &&
            step.step.accepted && step.neighbor_count == 0 &&
            step.global_subset_draw_count == 0 && words[0] == UINT64_C(1),
        "global-only one-state self proposal");

  saved_rng = rng;
  proposal_words[0] = UINT64_MAX;
  CHECK(mvmc_krylov_positive_sampler_step_mixture_rng(
            workspace, &policy, &one_state, &mixed, words, 1, &current,
            &rng, table_callback, amplitude, proposal_words, 1, &step) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !step.valid && rng.state == saved_rng.state &&
            rng.draws == saved_rng.draws && proposal_words[0] == 0 &&
            words[0] == UINT64_C(1),
        "zero-neighbor mixed policy fails without renormalization");

  amplitude->forced_status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  proposal_words[0] = UINT64_MAX;
  CHECK(mvmc_krylov_positive_sampler_step_mixture_rng(
            workspace, &policy, &regular, &mixed, words, 1, &current,
            &rng, table_callback, amplitude, proposal_words, 1, &step) ==
            MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
            !step.valid && rng.state == saved_rng.state &&
            rng.draws == saved_rng.draws && proposal_words[0] == 0 &&
            words[0] == UINT64_C(1) && current.accepted_generation == 1,
        "mixture evaluation failure is transactional");
}

static void test_mixture_restart_boundaries_body(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude,
    uint64_t plan_hash, const MVMCKrylovBoundedLimits *limits) {
  const MVMCKrylovPositiveGuidePolicy policy = guide_policy();
  const MVMCKrylovFockModel proposal = proposal_model();
  const size_t numerators[] = {0, 1, 1, 1, 1};
  const size_t denominators[] = {1, 8, 4, 2, 1};
  size_t candidate_index;

  for (candidate_index = 0; candidate_index < 5; ++candidate_index) {
    MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
    MVMCKrylovPositiveSamplerSnapshot full_current;
    MVMCKrylovPositiveSamplerSnapshot resumed_current;
    MVMCKrylovPositiveSamplerRng full_rng;
    MVMCKrylovPositiveSamplerRng resumed_rng;
    MVMCKrylovPositiveSamplerManifest split_manifest;
    MVMCKrylovPositiveSamplerTraceStatistics full_statistics;
    MVMCKrylovPositiveSamplerTraceStatistics resumed_statistics;
    uint64_t full_words[1] = {UINT64_C(1)};
    uint64_t resumed_words[1] = {UINT64_C(1)};
    int matches = 0;

    CHECK(mvmc_krylov_positive_sampler_proposal_policy_create(
              numerators[candidate_index], denominators[candidate_index],
              &proposal_policy) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_rng_seed(
                  UINT64_C(0x5034533352535431),
                  (uint64_t)candidate_index, &full_rng) ==
                  MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_rng_seed(
                  UINT64_C(0x5034533352535431),
                  (uint64_t)candidate_index, &resumed_rng) ==
                  MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_trace_statistics_reset(
                  &full_statistics) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_trace_statistics_reset(
                  &resumed_statistics) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_initialize(
                  workspace, &policy, full_words, 1, table_callback,
                  amplitude, &full_current) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_initialize(
                  workspace, &policy, resumed_words, 1, table_callback,
                  amplitude, &resumed_current) == MVMC_KRYLOV_STATUS_OK,
          "mixture restart setup %zu", candidate_index);
    (void)run_mixture_steps(
        workspace, &policy, &proposal, &proposal_policy, full_words,
        &full_current, &full_rng, amplitude, 40, &full_statistics);
    (void)run_mixture_steps(
        workspace, &policy, &proposal, &proposal_policy, resumed_words,
        &resumed_current, &resumed_rng, amplitude, 17,
        &resumed_statistics);
    CHECK(mvmc_krylov_positive_sampler_manifest_create(
              &policy, &proposal_policy, &proposal, limits, plan_hash,
              resumed_words, 1, &resumed_current, &resumed_rng,
              &split_manifest) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_manifest_matches(
                  &split_manifest, &policy, &proposal_policy, &proposal,
                  limits, plan_hash, resumed_words, 1, &resumed_current,
                  &resumed_rng, &matches) == MVMC_KRYLOV_STATUS_OK &&
              matches,
          "mixture split manifest %zu", candidate_index);
    (void)run_mixture_steps(
        workspace, &policy, &proposal, &proposal_policy, resumed_words,
        &resumed_current, &resumed_rng, amplitude, 23,
        &resumed_statistics);
    CHECK(full_words[0] == resumed_words[0] &&
              full_rng.state == resumed_rng.state &&
              full_rng.draws == resumed_rng.draws &&
              full_current.accepted_generation ==
                  resumed_current.accepted_generation &&
              full_current.guide.log_guide ==
                  resumed_current.guide.log_guide &&
              trace_statistics_equal(&full_statistics,
                                     &resumed_statistics),
          "mixture full/restart mismatch %zu", candidate_index);
  }
}

static void test_shell_policy_and_draw(void) {
  MVMCKrylovFockModel proposal;
  MVMCKrylovPositiveSamplerProposalPolicy one_third;
  MVMCKrylovPositiveSamplerProposalPolicy one_half;
  MVMCKrylovPositiveSamplerProposalPolicy invalid_policy;
  MVMCKrylovPositiveSamplerRng rng;
  uint64_t current = UINT64_C(0x33);
  uint64_t proposal_words = 0;
  uint64_t model_hash_one_third = 0;
  uint64_t model_hash_one_half = 0;
  size_t neighbor_draws = 0;
  size_t shell_draws = 0;
  size_t draw_index;

  memset(&proposal, 0, sizeof(proposal));
  proposal.site_count = 4;
  proposal.up_electron_count = 2;
  proposal.down_electron_count = 2;
  proposal.hermitian = 1;
  CHECK(mvmc_krylov_positive_sampler_shell_policy_create(
            1, 8, 1, 3, &one_third) == MVMC_KRYLOV_STATUS_OK &&
            one_third.valid &&
            one_third.kernel_id ==
                MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_SHELL_MIXTURE &&
            one_third.neighbor_numerator == 1 &&
            one_third.neighbor_denominator == 8 &&
            one_third.distance_numerator == 1 &&
            one_third.distance_denominator == 3 &&
            one_third.distance_rounding_rule ==
                MVMC_KRYLOV_PROPOSAL_DISTANCE_ROUND_HALF_UP_V1,
        "one-third shell policy create");
  CHECK(mvmc_krylov_positive_sampler_shell_policy_create(
            1, 8, 1, 2, &one_half) == MVMC_KRYLOV_STATUS_OK &&
            one_half.policy_hash != one_third.policy_hash,
        "one-half shell policy create/hash");
  CHECK(mvmc_krylov_positive_sampler_shell_policy_create(
            0, 8, 1, 2, &invalid_policy) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            mvmc_krylov_positive_sampler_shell_policy_create(
                8, 8, 1, 2, &invalid_policy) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            mvmc_krylov_positive_sampler_shell_policy_create(
                1, 8, 0, 2, &invalid_policy) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            mvmc_krylov_positive_sampler_shell_policy_create(
                1, 8, 3, 2, &invalid_policy) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "invalid shell policy accepted");
  CHECK(mvmc_krylov_positive_sampler_proposal_model_policy_hash(
            &one_third, &proposal, &current, 1,
            &model_hash_one_third) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_proposal_model_policy_hash(
                &one_half, &proposal, &current, 1,
                &model_hash_one_half) == MVMC_KRYLOV_STATUS_OK &&
            model_hash_one_third != 0 && model_hash_one_half != 0 &&
            model_hash_one_third != model_hash_one_half,
        "resolved shell model hashes");
  CHECK(mvmc_krylov_positive_sampler_rng_seed(
            UINT64_C(0x503453345348454c), UINT64_C(9), &rng) ==
            MVMC_KRYLOV_STATUS_OK,
        "shell draw rng seed");

  for (draw_index = 0; draw_index < 4096; ++draw_index) {
    MVMCKrylovPositiveSamplerProposalDrawResult draw;
    const MVMCKrylovStatus status =
        mvmc_krylov_positive_sampler_draw_mixture_rng(
            &proposal, &one_half, &current, 1, &rng, &proposal_words, 1,
            &draw);
    CHECK(status == MVMC_KRYLOV_STATUS_OK && draw.valid &&
              draw.component_draw_count == 1 &&
              draw.proposal_policy_hash == one_half.policy_hash &&
              draw.proposal_model_hash == model_hash_one_half &&
              draw.log_proposal_ratio == 0.0 &&
              isfinite(draw.uniform) && draw.uniform >= 0.0 &&
              draw.uniform < 1.0 && proposal_words != current,
          "shell mixture draw %zu", draw_index);
    if (status != MVMC_KRYLOV_STATUS_OK || !draw.valid) break;
    if (draw.component == MVMC_KRYLOV_PROPOSAL_COMPONENT_NEIGHBOR) {
      ++neighbor_draws;
      CHECK(draw.selected_neighbor_index < draw.neighbor_count &&
                draw.shell_draw_count == 0 &&
                draw.shell_distance == 0 && draw.shell_count == 0 &&
                replacement_distance_word(current, proposal_words, 4) == 1,
            "neighbor component metadata %zu", draw_index);
    } else if (draw.component == MVMC_KRYLOV_PROPOSAL_COMPONENT_SHELL) {
      ++shell_draws;
      CHECK(draw.selected_neighbor_index == SIZE_MAX &&
                draw.shell_max_distance == 4 &&
                draw.shell_distance == 2 &&
                draw.shell_up_distance + draw.shell_down_distance == 2 &&
                draw.shell_count == 18 &&
                replacement_distance_word(current, proposal_words, 4) == 2 &&
                isfinite(draw.shell_generation_seconds) &&
                draw.shell_generation_seconds >= 0.0,
            "shell component metadata %zu", draw_index);
    } else {
      CHECK(0, "unexpected shell policy component %d",
            (int)draw.component);
    }
    rng = draw.rng_after;
  }
  CHECK(neighbor_draws + shell_draws == 4096 &&
            neighbor_draws >= 430 && neighbor_draws <= 600 &&
            shell_draws >= 3496 && shell_draws <= 3666,
        "shell mixture component census neighbor=%zu shell=%zu",
        neighbor_draws, shell_draws);

  {
    MVMCKrylovFockModel pure_spin = proposal;
    uint64_t pure_current = UINT64_C(0xc3);
    uint64_t pure_proposal = 0;
    size_t attempt;
    int observed_shell = 0;
    pure_spin.pure_spin = 1;
    for (attempt = 0; attempt < 64 && !observed_shell; ++attempt) {
      MVMCKrylovPositiveSamplerProposalDrawResult draw;
      const MVMCKrylovStatus status =
          mvmc_krylov_positive_sampler_draw_mixture_rng(
              &pure_spin, &one_half, &pure_current, 1, &rng,
              &pure_proposal, 1, &draw);
      CHECK(status == MVMC_KRYLOV_STATUS_OK && draw.valid &&
                pure_proposal != pure_current,
            "pure-spin shell mixture draw %zu", attempt);
      if (status != MVMC_KRYLOV_STATUS_OK || !draw.valid) break;
      rng = draw.rng_after;
      if (draw.component == MVMC_KRYLOV_PROPOSAL_COMPONENT_SHELL) {
        observed_shell = 1;
        CHECK(draw.shell_max_distance == 2 && draw.shell_distance == 1 &&
                  draw.shell_up_distance == 1 &&
                  draw.shell_down_distance == 1 &&
                  draw.shell_count == 4 &&
                  ((pure_proposal >> 4) & UINT64_C(15)) ==
                      ((~pure_proposal) & UINT64_C(15)),
              "pure-spin shell metadata/complement");
      }
    }
    CHECK(observed_shell, "pure-spin shell component not observed");
  }

  {
    MVMCKrylovPositiveSamplerProposalDrawResult draw;
    const uint64_t saved_current = current;
    CHECK(mvmc_krylov_positive_sampler_draw_mixture_rng(
              &proposal, &one_half, &current, 1, &rng,
              &current, 1, &draw) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              current == saved_current && !draw.valid,
          "overlapping sampler buffers changed current");
  }
  invalid_policy = one_half;
  invalid_policy.distance_numerator = 1;
  invalid_policy.distance_denominator = 3;
  proposal_words = UINT64_MAX;
  CHECK(mvmc_krylov_positive_sampler_draw_mixture_rng(
            &proposal, &invalid_policy, &current, 1, &rng, &proposal_words,
            1, &(MVMCKrylovPositiveSamplerProposalDrawResult){0}) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && proposal_words == 0,
        "tampered shell policy hash accepted");
}

static void test_shell_restart_trace_and_failure_body(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude,
    uint64_t plan_hash, const MVMCKrylovBoundedLimits *limits) {
  const MVMCKrylovPositiveGuidePolicy policy = guide_policy();
  const MVMCKrylovFockModel proposal = proposal_model();
  const size_t distance_numerators[] = {1, 1, 2};
  const size_t distance_denominators[] = {3, 2, 3};
  size_t candidate_index;

  for (candidate_index = 0; candidate_index < 3; ++candidate_index) {
    MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
    MVMCKrylovPositiveSamplerSnapshot full_current;
    MVMCKrylovPositiveSamplerSnapshot resumed_current;
    MVMCKrylovPositiveSamplerRng full_rng;
    MVMCKrylovPositiveSamplerRng resumed_rng;
    MVMCKrylovPositiveSamplerManifest split_manifest;
    MVMCKrylovPositiveSamplerTraceStatistics full_statistics;
    MVMCKrylovPositiveSamplerTraceStatistics resumed_statistics;
    uint64_t full_words[1] = {UINT64_C(1)};
    uint64_t resumed_words[1] = {UINT64_C(1)};
    int matches = 0;

    CHECK(mvmc_krylov_positive_sampler_shell_policy_create(
              1, 8, distance_numerators[candidate_index],
              distance_denominators[candidate_index], &proposal_policy) ==
              MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_rng_seed(
                  UINT64_C(0x5034533452535431),
                  (uint64_t)candidate_index, &full_rng) ==
                  MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_rng_seed(
                  UINT64_C(0x5034533452535431),
                  (uint64_t)candidate_index, &resumed_rng) ==
                  MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_trace_statistics_reset(
                  &full_statistics) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_trace_statistics_reset(
                  &resumed_statistics) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_initialize(
                  workspace, &policy, full_words, 1, table_callback,
                  amplitude, &full_current) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_initialize(
                  workspace, &policy, resumed_words, 1, table_callback,
                  amplitude, &resumed_current) == MVMC_KRYLOV_STATUS_OK,
          "shell restart setup %zu", candidate_index);
    (void)run_mixture_steps(
        workspace, &policy, &proposal, &proposal_policy, full_words,
        &full_current, &full_rng, amplitude, 40, &full_statistics);
    (void)run_mixture_steps(
        workspace, &policy, &proposal, &proposal_policy, resumed_words,
        &resumed_current, &resumed_rng, amplitude, 17,
        &resumed_statistics);
    CHECK(mvmc_krylov_positive_sampler_manifest_create(
              &policy, &proposal_policy, &proposal, limits, plan_hash,
              resumed_words, 1, &resumed_current, &resumed_rng,
              &split_manifest) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_manifest_matches(
                  &split_manifest, &policy, &proposal_policy, &proposal,
                  limits, plan_hash, resumed_words, 1, &resumed_current,
                  &resumed_rng, &matches) == MVMC_KRYLOV_STATUS_OK &&
              matches,
          "shell split manifest %zu", candidate_index);
    (void)run_mixture_steps(
        workspace, &policy, &proposal, &proposal_policy, resumed_words,
        &resumed_current, &resumed_rng, amplitude, 23,
        &resumed_statistics);
    CHECK(full_words[0] == resumed_words[0] &&
              full_rng.state == resumed_rng.state &&
              full_rng.draws == resumed_rng.draws &&
              full_current.accepted_generation ==
                  resumed_current.accepted_generation &&
              trace_statistics_equal(&full_statistics,
                                     &resumed_statistics) &&
              full_statistics.neighbor_attempted_steps +
                      full_statistics.shell_attempted_steps ==
                  full_statistics.completed_steps &&
              full_statistics.shell_attempted_steps > 0 &&
              full_statistics.global_attempted_steps == 0 &&
              full_statistics.shell_accepted_steps +
                      full_statistics.shell_rejected_steps ==
                  full_statistics.shell_attempted_steps,
          "shell full/restart trace mismatch %zu", candidate_index);
  }

  {
    MVMCKrylovPositiveSamplerProposalPolicy proposal_policy = {0};
    MVMCKrylovPositiveSamplerSnapshot current = {0};
    MVMCKrylovPositiveSamplerProposalStepResult step;
    MVMCKrylovPositiveSamplerRng rng = {0};
    MVMCKrylovPositiveSamplerRng saved_rng;
    uint64_t words[1] = {UINT64_C(1)};
    uint64_t proposal_words[1] = {UINT64_MAX};
    CHECK(mvmc_krylov_positive_sampler_shell_policy_create(
              1, 8, 1, 2, &proposal_policy) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_rng_seed(
                  UINT64_C(0x503453344641494c), UINT64_C(4), &rng) ==
                  MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_positive_sampler_initialize(
                  workspace, &policy, words, 1, table_callback, amplitude,
                  &current) == MVMC_KRYLOV_STATUS_OK,
          "shell failure setup");
    saved_rng = rng;
    amplitude->forced_status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    CHECK(mvmc_krylov_positive_sampler_step_mixture_rng(
              workspace, &policy, &proposal, &proposal_policy, words, 1,
              &current, &rng, table_callback, amplitude, proposal_words, 1,
              &step) == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
              !step.valid && words[0] == UINT64_C(1) &&
              current.accepted_generation == 0 &&
              rng.state == saved_rng.state && rng.draws == saved_rng.draws &&
              proposal_words[0] == 0,
          "shell evaluation failure is not transactional");
    amplitude->forced_status = MVMC_KRYLOV_STATUS_OK;
  }
}

static void test_small_sector_exact_balance_body(
    MVMCKrylovBoundedWorkspace *workspace, TableAmplitude *amplitude) {
  const MVMCKrylovPositiveGuidePolicy policy = guide_policy();
  MVMCKrylovPositiveSamplerSnapshot state[2];
  MVMCKrylovPositiveGuideAcceptance forward;
  MVMCKrylovPositiveGuideAcceptance reverse;
  uint64_t words[2] = {UINT64_C(1), UINT64_C(2)};
  double guide_weight[2];
  double stationary[2];
  double transition[2][2];
  double forward_acceptance;
  double reverse_acceptance;
  const double other_proposal_probability = 0.75;
  const double tolerance = 256.0 * DBL_EPSILON;

  CHECK(mvmc_krylov_positive_sampler_initialize(
            workspace, &policy, &words[0], 1, table_callback, amplitude,
            &state[0]) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_initialize(
                workspace, &policy, &words[1], 1, table_callback,
                amplitude, &state[1]) == MVMC_KRYLOV_STATUS_OK,
        "small-sector exact state evaluation");
  CHECK(mvmc_krylov_positive_guide_acceptance(
            &state[0].guide, &state[1].guide, 0.0, 0.0,
            &forward) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_guide_acceptance(
                &state[1].guide, &state[0].guide, 0.0, 0.0,
                &reverse) == MVMC_KRYLOV_STATUS_OK,
        "small-sector acceptance ratios");
  guide_weight[0] = exp(state[0].guide.log_guide);
  guide_weight[1] = exp(state[1].guide.log_guide);
  stationary[0] = guide_weight[0] / (guide_weight[0] + guide_weight[1]);
  stationary[1] = guide_weight[1] / (guide_weight[0] + guide_weight[1]);
  forward_acceptance = exp(fmin(0.0, forward.log_acceptance_ratio));
  reverse_acceptance = exp(fmin(0.0, reverse.log_acceptance_ratio));
  transition[0][1] = other_proposal_probability * forward_acceptance;
  transition[0][0] = 1.0 - transition[0][1];
  transition[1][0] = other_proposal_probability * reverse_acceptance;
  transition[1][1] = 1.0 - transition[1][0];
  CHECK(close_double(transition[0][0] + transition[0][1], 1.0,
                     tolerance) &&
            close_double(transition[1][0] + transition[1][1], 1.0,
                         tolerance),
        "small-sector transition row sums");
  CHECK(close_double(stationary[0] * transition[0][1],
                     stationary[1] * transition[1][0], tolerance),
        "small-sector detailed balance");
  CHECK(close_double(stationary[0] * transition[0][0] +
                         stationary[1] * transition[1][0],
                     stationary[0], tolerance) &&
            close_double(stationary[0] * transition[0][1] +
                             stationary[1] * transition[1][1],
                         stationary[1], tolerance),
        "small-sector stationarity");
}

static void test_surrogate_policy_weight_and_draw(void) {
  MVMCKrylovPositiveGuidePolicy guide = guide_policy();
  MVMCKrylovPositiveGuidePolicy other_guide = guide_policy();
  MVMCKrylovPositiveGuidePolicy partial_guide = guide_policy();
  MVMCKrylovPositiveSamplerSurrogatePolicy k4_alpha10;
  MVMCKrylovPositiveSamplerSurrogatePolicy k16_alpha10;
  MVMCKrylovPositiveSamplerSurrogatePolicy k4_alpha100;
  MVMCKrylovPositiveSamplerSurrogatePolicy k32_alpha100;
  MVMCKrylovPositiveSamplerSurrogatePolicy partial_m1;
  MVMCKrylovPositiveSamplerSurrogatePolicy partial_m2;
  MVMCKrylovPositiveSamplerSurrogatePolicy invalid_surrogate;
  MVMCKrylovPositiveSamplerSurrogateWeightResult weight;
  MVMCScaledComplex zeroth;
  MVMCScaledComplex partial_values[3];
  MVMCKrylovFockModel model = proposal_model();
  MVMCKrylovPositiveSamplerRng rng;
  MVMCKrylovPositiveSamplerRng repeat_rng;
  MVMCKrylovPositiveSamplerSurrogateDrawResult draw;
  MVMCKrylovPositiveSamplerSurrogateDrawResult repeat_draw;
  SurrogateTable table;
  uint64_t current = UINT64_C(1);
  uint64_t proposal = UINT64_C(0xaaaaaaaaaaaaaaaa);
  uint64_t repeat_proposal = 0;
  size_t index;

  memset(&table, 0, sizeof(table));
  for (index = 0; index < 64; ++index) {
    table.weights[index] = 1.0 + (double)index;
  }
  guide.log_basis_scale[0] = 0.0;
  other_guide.log_basis_scale[0] = 0.0;
  other_guide.eta = 1.0;
  other_guide.policy_hash ^= UINT64_C(0x503453365152484f);
  partial_guide.order = 3;
  partial_guide.policy_hash ^= UINT64_C(0x5034533750415254);
  CHECK(mvmc_krylov_positive_sampler_surrogate_policy_create(
            4, 10, &guide, &k4_alpha10) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_surrogate_policy_create(
                16, 10, &guide, &k16_alpha10) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_surrogate_policy_create(
                4, 100, &guide, &k4_alpha100) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_surrogate_policy_create(
                32, 100, &guide, &k32_alpha100) == MVMC_KRYLOV_STATUS_OK &&
            k4_alpha10.valid &&
            k4_alpha10.policy_hash != k16_alpha10.policy_hash &&
            k4_alpha10.policy_hash != k4_alpha100.policy_hash,
        "surrogate policy identity");
  CHECK(mvmc_krylov_positive_sampler_surrogate_partial_policy_create(
            32, 1, 1, &partial_guide, &partial_m1) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_surrogate_partial_policy_create(
                32, 2, 1, &partial_guide, &partial_m2) ==
                MVMC_KRYLOV_STATUS_OK &&
            partial_m1.valid && partial_m2.valid &&
            partial_m1.partial_order == 1 && partial_m2.partial_order == 2 &&
            partial_m1.policy_hash != partial_m2.policy_hash,
        "partial surrogate policy identity");
  CHECK(mvmc_krylov_positive_sampler_surrogate_partial_policy_create(
            32, 4, 1, &partial_guide, &invalid_surrogate) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "partial surrogate order bounds");
  CHECK(mvmc_krylov_positive_sampler_surrogate_policy_create(
            SIZE_MAX, 1, &guide, &invalid_surrogate) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "surrogate step counter overflow rejected");
  CHECK(mvmc_scaled_complex_from_raw_testing(2.0, &zeroth) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_krylov_positive_sampler_surrogate_weight_zeroth(
                &k4_alpha10, &guide, &zeroth, &weight) ==
                MVMC_KRYLOV_STATUS_OK &&
            weight.valid && !weight.floor_only &&
            close_double(exp(weight.log_weight), 6.5, 1.0e-14),
        "surrogate finite zeroth weight");
  CHECK(mvmc_scaled_complex_make_exact_zero(&zeroth) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_krylov_positive_sampler_surrogate_weight_zeroth(
                &k4_alpha10, &guide, &zeroth, &weight) ==
                MVMC_KRYLOV_STATUS_OK &&
            weight.valid && weight.floor_only &&
            close_double(exp(weight.log_weight), 2.5, 1.0e-14),
        "surrogate floor-only zeroth weight");
  CHECK(mvmc_scaled_complex_make_numeric_zero(
            -20.0, 0.0, -INFINITY, 0.0, &zeroth) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_krylov_positive_sampler_surrogate_weight_zeroth(
                &k4_alpha10, &guide, &zeroth, &weight) ==
                MVMC_KRYLOV_STATUS_OK &&
            weight.valid && weight.floor_only &&
            close_double(exp(weight.log_weight), 2.5, 1.0e-14),
        "surrogate numeric-zero floor weight");
  CHECK(mvmc_scaled_complex_from_raw_testing(NAN + I, &zeroth) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_krylov_positive_sampler_surrogate_weight_zeroth(
                &k4_alpha10, &guide, &zeroth, &weight) ==
                MVMC_KRYLOV_STATUS_NONFINITE && !weight.valid,
        "surrogate nonfinite zeroth rejected");
  CHECK(mvmc_krylov_positive_sampler_surrogate_weight_zeroth(
            &k4_alpha10, &other_guide, &zeroth, &weight) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && !weight.valid,
        "surrogate guide mismatch rejected");

  CHECK(mvmc_scaled_complex_from_raw_testing(2.0, &partial_values[0]) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_scaled_complex_from_raw_testing(3.0, &partial_values[1]) ==
                MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_scaled_complex_make_exact_zero(&partial_values[2]) ==
                MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_krylov_positive_sampler_surrogate_weight_partial(
                &partial_m2, &partial_guide, partial_values, 3, &weight) ==
                MVMC_KRYLOV_STATUS_OK &&
            weight.valid && !weight.floor_only &&
            close_double(exp(weight.log_weight), 13.25, 1.0e-14),
        "partial surrogate finite and exact-zero weight");
  CHECK(mvmc_scaled_complex_make_numeric_zero(
            -20.0, 0.0, -INFINITY, 0.0, &partial_values[0]) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_scaled_complex_make_exact_zero(&partial_values[1]) ==
                MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_krylov_positive_sampler_surrogate_weight_partial(
                &partial_m1, &partial_guide, partial_values, 2, &weight) ==
                MVMC_KRYLOV_STATUS_OK &&
            weight.valid && weight.floor_only &&
            close_double(exp(weight.log_weight), 0.25, 1.0e-14),
        "partial surrogate all-zero floor weight");
  CHECK(mvmc_scaled_complex_from_raw_testing(NAN + I, &partial_values[1]) ==
            MVMC_PFAFFIAN_STATUS_OK &&
            mvmc_krylov_positive_sampler_surrogate_weight_partial(
                &partial_m1, &partial_guide, partial_values, 2, &weight) ==
                MVMC_KRYLOV_STATUS_NONFINITE && !weight.valid,
        "partial surrogate nonfinite rejected");
  CHECK(mvmc_krylov_positive_sampler_surrogate_weight_partial(
            &partial_m2, &partial_guide, partial_values, 2, &weight) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && !weight.valid &&
            mvmc_krylov_positive_sampler_surrogate_weight_zeroth(
                &partial_m2, &partial_guide, partial_values, &weight) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && !weight.valid,
        "partial surrogate value census and zeroth API separation");

  CHECK(mvmc_krylov_positive_sampler_rng_seed(
            UINT64_C(0x5034533655524e47), UINT64_C(3), &rng) ==
            MVMC_KRYLOV_STATUS_OK,
        "surrogate rng seed");
  repeat_rng = rng;
  CHECK(mvmc_krylov_positive_sampler_draw_surrogate_rng(
            &model, &k4_alpha10, &current, 1, log(table.weights[current]),
            &rng, surrogate_table_callback, &table, &proposal, 1, &draw) ==
            MVMC_KRYLOV_STATUS_OK &&
            draw.valid && draw.step_count == 4 &&
            draw.surrogate_evaluation_count == 4 &&
            draw.inner_accepted + draw.inner_rejected == 4 &&
            draw.inner_configuration_changes == draw.inner_accepted &&
            draw.surrogate_policy_hash == k4_alpha10.policy_hash &&
            draw.proposal_model_hash != 0 &&
            draw.final_configuration_hash != 0 &&
            close_double(draw.log_proposal_ratio,
                         log(table.weights[current]) -
                             log(table.weights[proposal]),
                         1.0e-14),
        "surrogate K-step draw");
  table.calls = 0;
  CHECK(mvmc_krylov_positive_sampler_draw_surrogate_rng(
            &model, &k4_alpha10, &current, 1, log(table.weights[current]),
            &repeat_rng, surrogate_table_callback, &table, &repeat_proposal,
            1, &repeat_draw) == MVMC_KRYLOV_STATUS_OK &&
            repeat_draw.valid && proposal == repeat_proposal &&
            draw.rng_after.state == repeat_draw.rng_after.state &&
            draw.rng_after.draws == repeat_draw.rng_after.draws &&
            draw.inner_accepted == repeat_draw.inner_accepted &&
            draw.inner_rejected == repeat_draw.inner_rejected,
        "surrogate deterministic replay");

  table.calls = 0;
  table.fail_on_call = 2;
  proposal = UINT64_C(0xaaaaaaaaaaaaaaaa);
  CHECK(mvmc_krylov_positive_sampler_draw_surrogate_rng(
            &model, &k4_alpha10, &current, 1, log(table.weights[current]),
            &rng, surrogate_table_callback, &table, &proposal, 1, &draw) ==
            MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
            !draw.valid && proposal == UINT64_C(0xaaaaaaaaaaaaaaaa) &&
            current == UINT64_C(1),
        "surrogate callback failure transaction");
  table.fail_on_call = 0;
  CHECK(mvmc_krylov_positive_sampler_draw_surrogate_rng(
            &model, &k4_alpha10, &current, 1, log(table.weights[current]),
            &rng, surrogate_table_callback, &table, &current, 1, &draw) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && current == UINT64_C(1),
        "surrogate overlapping buffers rejected");

  for (index = 0; index < 64; ++index) {
    table.weights[index] = index == current ? 1.0 : DBL_MIN;
  }
  table.calls = 0;
  CHECK(mvmc_krylov_positive_sampler_draw_surrogate_rng(
            &model, &k4_alpha10, &current, 1, 0.0, &rng,
            surrogate_table_callback, &table, &proposal, 1, &draw) ==
            MVMC_KRYLOV_STATUS_OK &&
            draw.valid && !draw.final_configuration_changed &&
            proposal == current && draw.inner_accepted == 0 &&
            draw.inner_rejected == 4,
        "surrogate final-self transaction");
  for (index = 0; index < 64; ++index) table.weights[index] = NAN;
  table.weights[current] = 1.0;
  proposal = UINT64_C(0xaaaaaaaaaaaaaaaa);
  table.calls = 0;
  CHECK(mvmc_krylov_positive_sampler_draw_surrogate_rng(
            &model, &k4_alpha10, &current, 1, 0.0, &rng,
            surrogate_table_callback, &table, &proposal, 1, &draw) ==
            MVMC_KRYLOV_STATUS_NONFINITE && !draw.valid &&
            proposal == UINT64_C(0xaaaaaaaaaaaaaaaa) &&
            current == UINT64_C(1),
        "surrogate nonfinite callback transaction");
  for (index = 0; index < 64; ++index) {
    table.weights[index] = 1.0 + (double)index;
  }

  {
    const size_t site_counts[] = {4, 6, 8};
    const MVMCKrylovPositiveSamplerSurrogatePolicy *policies[] = {
        &k4_alpha10, &k16_alpha10, &k32_alpha100};
    size_t case_index;
    for (case_index = 0; case_index < 3; ++case_index) {
      MVMCKrylovFockModel pure_spin = {0};
      MVMCKrylovStatus status;
      const size_t site_count = site_counts[case_index];
      const size_t half = site_count / 2;
      const uint64_t site_mask =
          (UINT64_C(1) << site_count) - UINT64_C(1);
      const uint64_t up = (UINT64_C(1) << half) - UINT64_C(1);
      uint64_t pure_current = up | ((site_mask ^ up) << site_count);
      uint64_t pure_proposal = 0;
      pure_spin.site_count = site_count;
      pure_spin.up_electron_count = half;
      pure_spin.down_electron_count = half;
      pure_spin.pure_spin = 1;
      pure_spin.hermitian = 1;
      table.calls = 0;
      status = mvmc_krylov_positive_sampler_draw_surrogate_rng(
          &pure_spin, policies[case_index], &pure_current, 1,
          log(table.weights[pure_current & UINT64_C(63)]), &rng,
          surrogate_table_callback, &table, &pure_proposal, 1, &draw);
      CHECK(status == MVMC_KRYLOV_STATUS_OK && draw.valid &&
                ((pure_proposal >> site_count) & site_mask) ==
                    ((~pure_proposal) & site_mask) &&
                mvmc_krylov_fock_validate(&pure_spin, &pure_proposal, 1) ==
                    MVMC_KRYLOV_STATUS_OK,
            "surrogate pure-spin L%zu draw status=%d", site_count,
            (int)status);
    }
  }
  {
    MVMCKrylovFockModel multiword = {0};
    MVMCKrylovStatus status;
    uint64_t multi_current[2] = {UINT64_C(1), UINT64_C(0)};
    uint64_t multi_proposal[2] = {UINT64_MAX, UINT64_MAX};
    multiword.site_count = 40;
    multiword.up_electron_count = 1;
    multiword.hermitian = 1;
    table.calls = 0;
    status = mvmc_krylov_positive_sampler_draw_surrogate_rng(
        &multiword, &k4_alpha10, multi_current, 2,
        log(table.weights[1]), &rng, surrogate_table_callback, &table,
        multi_proposal, 2, &draw);
    CHECK(status == MVMC_KRYLOV_STATUS_OK &&
              draw.valid && draw.surrogate_evaluation_count == 4,
          "surrogate multiword draw status=%d", (int)status);
  }
}

static void test_accept_reject(void) {
  with_workspace(test_accept_reject_body);
}

static void test_persistent_session(void) {
  with_workspace(test_persistent_session_body);
}

static void test_failure_transaction(void) {
  with_workspace(test_failure_transaction_body);
}

static void test_policy_mismatch(void) {
  with_workspace(test_policy_mismatch_body);
}

static void test_selected_neighbor_trace(void) {
  with_workspace(test_selected_neighbor_trace_body);
}

static void test_rng_manifest_restart(void) {
  with_workspace_and_plan(test_rng_manifest_restart_body);
}

static void test_full_markov_trace_statistics(void) {
  with_workspace(test_full_markov_trace_statistics_body);
}

static void test_floor_supported_zero_evidence(void) {
  with_workspace(test_floor_supported_zero_evidence_body);
}

static void test_zero_fraction_preserves_legacy(void) {
  with_workspace(test_zero_fraction_preserves_legacy_body);
}

static void test_mixture_counters(void) {
  with_workspace(test_mixture_counters_body);
}

static void test_mixture_degenerate_and_failure(void) {
  with_workspace(test_mixture_degenerate_and_failure_body);
}

static void test_mixture_restart_boundaries(void) {
  with_workspace_and_plan(test_mixture_restart_boundaries_body);
}

static void test_small_sector_exact_balance(void) {
  with_workspace(test_small_sector_exact_balance_body);
}

static void test_shell_restart_trace_and_failure(void) {
  with_workspace_and_plan(test_shell_restart_trace_and_failure_body);
}

int main(void) {
  test_accept_reject();
  test_persistent_session();
  test_failure_transaction();
  test_policy_mismatch();
  test_selected_neighbor_trace();
  test_rng_manifest_restart();
  test_full_markov_trace_statistics();
  test_floor_supported_zero_evidence();
  test_zero_fraction_preserves_legacy();
  test_mixture_counters();
  test_mixture_degenerate_and_failure();
  test_mixture_restart_boundaries();
  test_shell_policy_and_draw();
  test_shell_restart_trace_and_failure();
  test_surrogate_policy_weight_and_draw();
  test_small_sector_exact_balance();
  if (failures == 0) {
    printf("krylov positive sampler unit: PASS\n");
  }
  return failures == 0 ? 0 : 1;
}
