#include "krylov_final_state_chain.h"
#include "krylov_streaming_statistics.h"

#include <complex.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

#define DRIVER_SCHEMA 2
#define MAX_STATES 4

typedef struct {
  double complex values[64];
} TableAmplitude;

typedef struct {
  const char *id;
  MVMCKrylovFockModel model;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovFinalStatePolicy policy;
  TableAmplitude amplitude;
  uint64_t states[MAX_STATES];
  size_t state_count;
  uint64_t amplitude_generation_hash;
} Fixture;

typedef struct {
  uint64_t attempted;
  uint64_t accepted;
  uint64_t changed;
  uint64_t exact_zero_proposals;
  uint64_t trace_hash;
  double *energy_real;
  double *energy_imag;
  double *second;
  size_t sample_count;
  MVMCKrylovFinalStateChainManifest final_manifest;
  uint64_t final_words;
  MVMCKrylovPositiveSamplerRng final_rng;
  MVMCKrylovFinalStateSnapshot final_snapshot;
} Trace;

static int world_rank = 0;
static int world_size = 1;

static void hash_u64(uint64_t *hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    *hash ^= value & UINT64_C(0xff);
    *hash *= UINT64_C(1099511628211);
    value >>= 8;
  }
}

static uint64_t double_bits(double value) {
  uint64_t bits = 0;
  memcpy(&bits, &value, sizeof(bits));
  return bits;
}

static int parse_size(const char *text, size_t *value) {
  char *end = NULL;
  uintmax_t parsed;
  if (text == NULL || value == NULL || text[0] == '-') return 0;
  parsed = strtoumax(text, &end, 0);
  if (end == text || *end != '\0' || parsed > SIZE_MAX) return 0;
  *value = (size_t)parsed;
  return 1;
}

static int parse_u64(const char *text, uint64_t *value) {
  char *end = NULL;
  uintmax_t parsed;
  if (text == NULL || value == NULL || text[0] == '-') return 0;
  parsed = strtoumax(text, &end, 0);
  if (end == text || *end != '\0' || parsed > UINT64_MAX) return 0;
  *value = (uint64_t)parsed;
  return 1;
}

static MVMCKrylovStatus table_callback(
    const uint64_t *words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
  const TableAmplitude *table = (const TableAmplitude *)context;
  const size_t index = word_count == 1
                           ? (size_t)(words[0] & UINT64_C(63))
                           : 63;
  if (table == NULL || result == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(result, 0, sizeof(*result));
  if (creal(table->values[index]) == 0.0 &&
      cimag(table->values[index]) == 0.0) {
    if (mvmc_scaled_complex_make_exact_zero(&result->value) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    result->exact_zero_component_count = 1;
  } else {
    if (mvmc_scaled_complex_from_raw_testing(
            table->values[index], &result->value) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    }
    result->well_pivoted_component_count = 1;
  }
  result->local_factorization_count = 1;
  result->global_factorization_count = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovBoundedLimits bounded_limits(uint64_t hash) {
  MVMCKrylovBoundedLimits limits;
  memset(&limits, 0, sizeof(limits));
  limits.amplitude_policy_hash = hash;
  limits.cache_bytes = 16384;
  limits.max_row_transitions = 64;
  limits.max_workspace_bytes = (size_t)4 * 1024 * 1024;
  limits.max_node_expansions = 10000;
  limits.max_terminal_amplitude_calls = 10000;
  limits.max_total_row_transitions = 100000;
  limits.max_order = 2;
  return limits;
}

static MVMCKrylovFockModel two_state_model(void) {
  static const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
  };
  static const MVMCKrylovHamiltonianTerm terms[] = {
      {1.0, 0, 2, 1, 0},
      {1.0, 2, 2, 1, 1},
  };
  MVMCKrylovFockModel model;
  memset(&model, 0, sizeof(model));
  model.site_count = 2;
  model.up_electron_count = 1;
  model.hermitian = 1;
  model.terms = terms;
  model.term_count = sizeof(terms) / sizeof(terms[0]);
  model.operators = operators;
  model.operator_count = sizeof(operators) / sizeof(operators[0]);
  return model;
}

static MVMCKrylovFockModel ring4_model(void) {
  static const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
      {MVMC_KRYLOV_FERMION_CREATE, 2},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 2},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 3},
      {MVMC_KRYLOV_FERMION_CREATE, 3},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
      {MVMC_KRYLOV_FERMION_CREATE, 3},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 3},
  };
  static const MVMCKrylovHamiltonianTerm terms[] = {
      {1.0, 0, 2, 1, 0},   {1.0, 2, 2, 1, 1},
      {1.0, 4, 2, 1, 2},   {1.0, 6, 2, 1, 3},
      {1.0, 8, 2, 1, 4},   {1.0, 10, 2, 1, 5},
      {1.0, 12, 2, 1, 6},  {1.0, 14, 2, 1, 7},
  };
  MVMCKrylovFockModel model;
  memset(&model, 0, sizeof(model));
  model.site_count = 4;
  model.up_electron_count = 1;
  model.hermitian = 1;
  model.terms = terms;
  model.term_count = sizeof(terms) / sizeof(terms[0]);
  model.operators = operators;
  model.operator_count = sizeof(operators) / sizeof(operators[0]);
  return model;
}

static int fixture_create(const char *id, Fixture *fixture) {
  double complex coefficient[2];
  MVMCKrylovStatus status;
  memset(fixture, 0, sizeof(*fixture));
  fixture->id = id;
  if (strcmp(id, "two_state_support_bridge") == 0 ||
      strcmp(id, "two_state_variance_gap") == 0) {
    fixture->model = two_state_model();
    fixture->limits = bounded_limits(UINT64_C(0x5035433253544154));
    fixture->states[0] = UINT64_C(1);
    fixture->states[1] = UINT64_C(2);
    fixture->state_count = 2;
    fixture->amplitude.values[1] = 1.0;
    coefficient[0] = 1.0;
    coefficient[1] = strcmp(id, "two_state_support_bridge") == 0
                         ? 2.0
                         : 0.0;
    fixture->amplitude_generation_hash = UINT64_C(0x5035433253544745);
  } else if (strcmp(id, "complex_ring4_near_zero") == 0) {
    fixture->model = ring4_model();
    fixture->limits = bounded_limits(UINT64_C(0x50354352494e4734));
    fixture->states[0] = UINT64_C(1);
    fixture->states[1] = UINT64_C(2);
    fixture->states[2] = UINT64_C(4);
    fixture->states[3] = UINT64_C(8);
    fixture->state_count = 4;
    fixture->amplitude.values[1] = 1.0;
    fixture->amplitude.values[2] = 0.0;
    fixture->amplitude.values[4] = 0.05001 - 0.025 * I;
    fixture->amplitude.values[8] = -0.2 + 0.1 * I;
    coefficient[0] = 1.0;
    coefficient[1] = 0.25;
    fixture->amplitude_generation_hash = UINT64_C(0x503543523447454e);
  } else {
    return 0;
  }
  status = mvmc_krylov_final_state_policy_create(
      1, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC,
      fixture->amplitude_generation_hash, coefficient, 2,
      &fixture->policy);
  return status == MVMC_KRYLOV_STATUS_OK;
}

static int export_scaled(const MVMCScaledComplex *value,
                         double complex *raw) {
  const MVMCScaledComplexExportStatus status =
      mvmc_scaled_complex_export_common_scale(value, 0.0, raw);
  return status == MVMC_SCALED_EXPORT_OK ||
         status == MVMC_SCALED_EXPORT_EXACT_ZERO;
}

static int begin_workspace(
    Fixture *fixture, MVMCKrylovBoundedPlan **plan,
    MVMCKrylovBoundedWorkspace **workspace) {
  MVMCKrylovStatus status = mvmc_bounded_krylov_plan_create(
      &fixture->model, &fixture->limits, plan);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_workspace_create(*plan, workspace);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_session_begin(
        *workspace, table_callback, &fixture->amplitude,
        fixture->amplitude_generation_hash);
  }
  return status == MVMC_KRYLOV_STATUS_OK;
}

static void end_workspace(MVMCKrylovBoundedPlan *plan,
                          MVMCKrylovBoundedWorkspace *workspace) {
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
}

static int print_anchors(Fixture *fixture,
                         MVMCKrylovBoundedWorkspace *workspace) {
  size_t index;
  size_t support_count = 0;
  uint64_t anchor_hash = UINT64_C(1469598103934665603);
  for (index = 0; index < fixture->state_count; ++index) {
    MVMCKrylovBoundedResult krylov;
    MVMCKrylovFinalStateEvaluation final_state;
    double complex raw[3];
    double complex psi = NAN + I * NAN;
    double complex hpsi = NAN + I * NAN;
    int order;
    if (mvmc_bounded_krylov_session_evaluate_bound(
            workspace, &fixture->states[index], 1, table_callback,
            &fixture->amplitude, &krylov) != MVMC_KRYLOV_STATUS_OK ||
        mvmc_krylov_final_state_evaluate(
            &fixture->policy, &krylov, &final_state) !=
            MVMC_KRYLOV_STATUS_OK) {
      return 0;
    }
    for (order = 0; order <= 2; ++order) {
      if (!export_scaled(&krylov.value[order], &raw[order])) return 0;
    }
    if (!export_scaled(&final_state.amplitude, &psi) ||
        !export_scaled(&final_state.hamiltonian_amplitude, &hpsi)) {
      return 0;
    }
    if (final_state.sampleable) ++support_count;
    hash_u64(&anchor_hash, fixture->states[index]);
    for (order = 0; order <= 2; ++order) {
      hash_u64(&anchor_hash, double_bits(creal(raw[order])));
      hash_u64(&anchor_hash, double_bits(cimag(raw[order])));
    }
    hash_u64(&anchor_hash, double_bits(creal(psi)));
    hash_u64(&anchor_hash, double_bits(cimag(psi)));
    hash_u64(&anchor_hash, double_bits(creal(hpsi)));
    hash_u64(&anchor_hash, double_bits(cimag(hpsi)));
    hash_u64(&anchor_hash,
             (uint64_t)(unsigned int)final_state.sampleable);
    if (world_rank == 0) {
      printf("ANCHOR index=%zu words=%" PRIu64
             " v0_real=%.17g v0_imag=%.17g"
             " v1_real=%.17g v1_imag=%.17g"
             " v2_real=%.17g v2_imag=%.17g"
             " psi_real=%.17g psi_imag=%.17g"
             " hpsi_real=%.17g hpsi_imag=%.17g"
             " sampleable=%d psi0_state=%d psi_state=%d hpsi_state=%d\n",
             index, fixture->states[index], creal(raw[0]), cimag(raw[0]),
             creal(raw[1]), cimag(raw[1]), creal(raw[2]), cimag(raw[2]),
             creal(psi), cimag(psi), creal(hpsi), cimag(hpsi),
             final_state.sampleable, (int)krylov.value[0].state,
             (int)final_state.amplitude.state,
             (int)final_state.hamiltonian_amplitude.state);
    }
  }
#ifdef _mpi_use
  {
    uint64_t local[] = {anchor_hash, (uint64_t)support_count};
    uint64_t minimum[2];
    uint64_t maximum[2];
    MPI_Allreduce(local, minimum, 2, MPI_UINT64_T, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(local, maximum, 2, MPI_UINT64_T, MPI_MAX,
                  MPI_COMM_WORLD);
    if (memcmp(minimum, maximum, sizeof(minimum)) != 0) return 0;
  }
#endif
  if (world_rank == 0) {
    printf("CONNECTIVITY sector_count=%zu target_support_count=%zu"
           " positive_pair_count=%zu connected_by_policy=%d"
           " pairwise_enumerated=0 proposal=global_uniform\n",
           fixture->state_count, support_count,
           support_count * support_count, support_count > 0);
  }
  return support_count > 0;
}

static int trace_rank_invariant(const Trace *trace) {
#ifdef _mpi_use
  uint64_t local[] = {
      trace->trace_hash,
      trace->final_words,
      trace->final_rng.state,
      trace->final_rng.draws,
      trace->final_snapshot.accepted_generation,
      trace->final_manifest.component_state_hash,
  };
  uint64_t minimum[6];
  uint64_t maximum[6];
  MPI_Allreduce(local, minimum, 6, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(local, maximum, 6, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  return memcmp(minimum, maximum, sizeof(minimum)) == 0;
#else
  (void)trace;
  return 1;
#endif
}

static void trace_destroy(Trace *trace) {
  free(trace->energy_real);
  free(trace->energy_imag);
  free(trace->second);
  memset(trace, 0, sizeof(*trace));
}

static int trace_allocate(size_t sample_count, Trace *trace) {
  memset(trace, 0, sizeof(*trace));
  trace->trace_hash = UINT64_C(1469598103934665603);
  trace->sample_count = sample_count;
  if (sample_count == 0) return 1;
  trace->energy_real = (double *)calloc(sample_count, sizeof(double));
  trace->energy_imag = (double *)calloc(sample_count, sizeof(double));
  trace->second = (double *)calloc(sample_count, sizeof(double));
  return trace->energy_real != NULL && trace->energy_imag != NULL &&
         trace->second != NULL;
}

static int run_range(
    Fixture *fixture, MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    size_t burn_in, size_t total_samples, size_t begin_step,
    size_t end_step, uint64_t *words, MVMCKrylovFinalStateSnapshot *current,
    MVMCKrylovPositiveSamplerRng *rng, Trace *trace) {
  size_t step_index;
  for (step_index = begin_step; step_index < end_step; ++step_index) {
    MVMCKrylovFinalStateChainStepResult step;
    uint64_t proposal_words = 0;
    double complex energy = NAN + I * NAN;
    double complex second = NAN + I * NAN;
    const MVMCKrylovStatus status =
        mvmc_krylov_final_state_chain_step_mixture_rng(
            workspace, &fixture->policy, &fixture->model, proposal_policy,
            words, 1, current, rng, table_callback, &fixture->amplitude,
            &proposal_words, 1, &step);
    if (status != MVMC_KRYLOV_STATUS_OK || !step.valid) return 0;
    ++trace->attempted;
    trace->accepted += (uint64_t)(unsigned int)step.step.accepted;
    trace->changed += (uint64_t)(unsigned int)step.configuration_changed;
    trace->exact_zero_proposals +=
        (uint64_t)(unsigned int)step.step.acceptance.exact_zero_proposal;
    hash_u64(&trace->trace_hash, *words);
    hash_u64(&trace->trace_hash,
             (uint64_t)(unsigned int)step.step.accepted);
    hash_u64(&trace->trace_hash, rng->draws);
    if (step_index >= burn_in) {
      const size_t sample = step_index - burn_in;
      if (sample >= total_samples ||
          mvmc_krylov_final_state_observable_measure(
              MVMC_KRYLOV_FINAL_OBSERVABLE_ENERGY_COMPLEX,
              &current->final_state, &energy) != MVMC_KRYLOV_STATUS_OK ||
          mvmc_krylov_final_state_observable_measure(
              MVMC_KRYLOV_FINAL_OBSERVABLE_LOCAL_ENERGY_ABS_SQUARED_DIAGNOSTIC,
              &current->final_state, &second) != MVMC_KRYLOV_STATUS_OK) {
        return 0;
      }
      trace->energy_real[sample] = creal(energy);
      trace->energy_imag[sample] = cimag(energy);
      trace->second[sample] = creal(second);
      hash_u64(&trace->trace_hash, double_bits(creal(energy)));
      hash_u64(&trace->trace_hash, double_bits(cimag(energy)));
      hash_u64(&trace->trace_hash, double_bits(creal(second)));
    }
  }
  return 1;
}

static int run_trace(Fixture *fixture, size_t burn_in, size_t sample_count,
                     size_t restart_split, uint64_t seed, uint64_t stream,
                     int restart_mode, Trace *trace) {
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
  MVMCKrylovFinalStateSnapshot current;
  MVMCKrylovPositiveSamplerRng rng;
  uint64_t words = fixture->states[0];
  const size_t total_steps = burn_in + sample_count;
  int ok = 0;
  if (!trace_allocate(sample_count, trace) ||
      !begin_workspace(fixture, &plan, &workspace) ||
      mvmc_krylov_positive_sampler_proposal_policy_create(
          1, 1, &proposal_policy) != MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_final_state_sampler_initialize(
          workspace, &fixture->policy, &words, 1, table_callback,
          &fixture->amplitude, &current) != MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_positive_sampler_rng_seed(seed, stream, &rng) !=
          MVMC_KRYLOV_STATUS_OK) {
    goto cleanup;
  }
  if (!restart_mode) {
    if (!run_range(fixture, workspace, &proposal_policy, burn_in,
                   sample_count, 0, total_steps, &words, &current, &rng,
                   trace)) {
      goto cleanup;
    }
  } else {
    MVMCKrylovFinalStateChainManifest restart_manifest;
    if (restart_split == 0 || restart_split >= total_steps ||
        !run_range(fixture, workspace, &proposal_policy, burn_in,
                   sample_count, 0, restart_split, &words, &current, &rng,
                   trace) ||
        mvmc_krylov_final_state_chain_manifest_create(
            &fixture->policy, &proposal_policy, &fixture->model,
            &fixture->limits, mvmc_bounded_krylov_plan_hash(plan), &words,
            1, &current, &rng, &restart_manifest) !=
            MVMC_KRYLOV_STATUS_OK) {
      goto cleanup;
    }
    end_workspace(plan, workspace);
    plan = NULL;
    workspace = NULL;
    if (!begin_workspace(fixture, &plan, &workspace)) goto cleanup;
    mvmc_bounded_krylov_session_end(workspace);
    memset(&current, 0, sizeof(current));
    memset(&rng, 0, sizeof(rng));
    if (mvmc_krylov_final_state_chain_restart_restore(
            workspace, &restart_manifest, &fixture->policy,
            &proposal_policy, &fixture->model, &fixture->limits,
            mvmc_bounded_krylov_plan_hash(plan), &words, 1,
            table_callback, &fixture->amplitude, &current, &rng) !=
            MVMC_KRYLOV_STATUS_OK ||
        !run_range(fixture, workspace, &proposal_policy, burn_in,
                   sample_count, restart_split, total_steps, &words,
                   &current, &rng, trace)) {
      goto cleanup;
    }
  }
  if (mvmc_krylov_final_state_chain_manifest_create(
          &fixture->policy, &proposal_policy, &fixture->model,
          &fixture->limits, mvmc_bounded_krylov_plan_hash(plan), &words, 1,
          &current, &rng, &trace->final_manifest) !=
      MVMC_KRYLOV_STATUS_OK) {
    goto cleanup;
  }
  trace->final_words = words;
  trace->final_rng = rng;
  trace->final_snapshot = current;
  ok = 1;
cleanup:
  if (workspace != NULL || plan != NULL) end_workspace(plan, workspace);
  if (!ok) trace_destroy(trace);
  return ok;
}

static int traces_equal(const Trace *left, const Trace *right) {
  return left->attempted == right->attempted &&
         left->accepted == right->accepted && left->changed == right->changed &&
         left->exact_zero_proposals == right->exact_zero_proposals &&
         left->trace_hash == right->trace_hash &&
         left->sample_count == right->sample_count &&
         left->final_words == right->final_words &&
         memcmp(&left->final_rng, &right->final_rng,
                sizeof(left->final_rng)) == 0 &&
         left->final_snapshot.accepted_generation ==
             right->final_snapshot.accepted_generation &&
         memcmp(left->energy_real, right->energy_real,
                left->sample_count * sizeof(double)) == 0 &&
         memcmp(left->energy_imag, right->energy_imag,
                left->sample_count * sizeof(double)) == 0 &&
         memcmp(left->second, right->second,
                left->sample_count * sizeof(double)) == 0;
}

static double mean_range(const double *values, size_t begin, size_t end) {
  double sum = 0.0;
  double compensation = 0.0;
  size_t index;
  for (index = begin; index < end; ++index) {
    const double corrected = values[index] - compensation;
    const double updated = sum + corrected;
    compensation = (updated - sum) - corrected;
    sum = updated;
  }
  return sum / (double)(end - begin);
}

static int compare_double(const void *left, const void *right) {
  const double a = *(const double *)left;
  const double b = *(const double *)right;
  return (a > b) - (a < b);
}

static double quantile_sorted(const double *values, size_t count,
                              double probability) {
  const double position = probability * (double)(count - 1);
  const size_t lower = (size_t)floor(position);
  const size_t upper = (size_t)ceil(position);
  const double fraction = position - (double)lower;
  return values[lower] + fraction * (values[upper] - values[lower]);
}

static int print_trace(const Trace *trace, size_t burn_in,
                       size_t block_length, int restart_pass) {
  MVMCKrylovTauIntResult tau_real;
  MVMCKrylovTauIntResult tau_imag;
  MVMCKrylovTauIntResult tau_second;
  double *absolute = NULL;
  double *deviation = NULL;
  double median;
  double mad;
  size_t block;
  size_t block_count;
  size_t index;
  if (trace->sample_count == 0 || block_length == 0 ||
      trace->sample_count % block_length != 0 ||
      mvmc_krylov_tau_int_geyer_initial_positive(
          trace->energy_real, trace->sample_count, &tau_real) !=
          MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_tau_int_geyer_initial_positive(
          trace->energy_imag, trace->sample_count, &tau_imag) !=
          MVMC_KRYLOV_STATUS_OK ||
      mvmc_krylov_tau_int_geyer_initial_positive(
          trace->second, trace->sample_count, &tau_second) !=
          MVMC_KRYLOV_STATUS_OK) {
    return 0;
  }
  absolute = (double *)malloc(trace->sample_count * sizeof(double));
  deviation = (double *)malloc(trace->sample_count * sizeof(double));
  if (absolute == NULL || deviation == NULL) {
    free(absolute);
    free(deviation);
    return 0;
  }
  for (index = 0; index < trace->sample_count; ++index) {
    absolute[index] = hypot(trace->energy_real[index],
                            trace->energy_imag[index]);
  }
  qsort(absolute, trace->sample_count, sizeof(double), compare_double);
  median = quantile_sorted(absolute, trace->sample_count, 0.5);
  for (index = 0; index < trace->sample_count; ++index) {
    deviation[index] = fabs(absolute[index] - median);
  }
  qsort(deviation, trace->sample_count, sizeof(double), compare_double);
  mad = quantile_sorted(deviation, trace->sample_count, 0.5);
  block_count = trace->sample_count / block_length;
  if (world_rank == 0) {
    printf("TRACE attempted=%" PRIu64 " accepted=%" PRIu64
           " changed=%" PRIu64 " exact_zero_proposals=%" PRIu64
           " acceptance_rate=%.17g trace_hash=%" PRIu64
           " final_words=%" PRIu64 " rng_draws=%" PRIu64
           " accepted_generation=%" PRIu64 " restart_pass=%d"
           " burn_in=%zu sample_count=%zu block_length=%zu block_count=%zu\n",
           trace->attempted, trace->accepted, trace->changed,
           trace->exact_zero_proposals,
           (double)trace->accepted / (double)trace->attempted,
           trace->trace_hash, trace->final_words, trace->final_rng.draws,
           trace->final_snapshot.accepted_generation, restart_pass,
           burn_in, trace->sample_count, block_length, block_count);
    printf("TAU energy_real=%.17g energy_imag=%.17g second=%.17g\n",
           tau_real.tau_int, tau_imag.tau_int, tau_second.tau_int);
    printf("TAIL median=%.17g mad=%.17g q90=%.17g q99=%.17g"
           " q999=%.17g maximum=%.17g\n",
           median, mad, quantile_sorted(absolute, trace->sample_count, 0.9),
           quantile_sorted(absolute, trace->sample_count, 0.99),
           quantile_sorted(absolute, trace->sample_count, 0.999),
           absolute[trace->sample_count - 1]);
    for (block = 0; block < block_count; ++block) {
      const size_t begin = block * block_length;
      const size_t end = begin + block_length;
      printf("BLOCK index=%zu energy_real=%.17g energy_imag=%.17g"
             " second=%.17g\n",
             block, mean_range(trace->energy_real, begin, end),
             mean_range(trace->energy_imag, begin, end),
             mean_range(trace->second, begin, end));
    }
    printf("MANIFEST state_version=%" PRIu64
           " final_policy_hash=%" PRIu64
           " coefficient_source=%" PRIu64
           " coefficient_provenance_hash=%" PRIu64
           " proposal_policy_hash=%" PRIu64
           " proposal_model_hash=%" PRIu64
           " bounded_plan_hash=%" PRIu64
           " amplitude_policy_hash=%" PRIu64
           " amplitude_generation_hash=%" PRIu64
           " component_state_hash=%" PRIu64 "\n",
           trace->final_manifest.state_version,
           trace->final_manifest.final_policy_hash,
           trace->final_manifest.coefficient_source,
           trace->final_manifest.coefficient_provenance_hash,
           trace->final_manifest.proposal_policy_hash,
           trace->final_manifest.proposal_model_hash,
           trace->final_manifest.bounded_plan_hash,
           trace->final_manifest.amplitude_policy_hash,
           trace->final_manifest.amplitude_generation_hash,
           trace->final_manifest.component_state_hash);
  }
  free(absolute);
  free(deviation);
  return 1;
}

int main(int argc, char **argv) {
  Fixture fixture;
  MVMCKrylovBoundedPlan *anchor_plan = NULL;
  MVMCKrylovBoundedWorkspace *anchor_workspace = NULL;
  Trace full;
  Trace restarted;
  size_t sample_count;
  size_t burn_in;
  size_t block_length;
  size_t restart_split;
  uint64_t seed;
  uint64_t stream;
  int local_ok = 1;
  int global_ok = 1;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif
  memset(&full, 0, sizeof(full));
  memset(&restarted, 0, sizeof(restarted));
  if (argc != 8 || !fixture_create(argv[1], &fixture) ||
      !parse_size(argv[2], &sample_count) ||
      !parse_size(argv[3], &burn_in) ||
      !parse_size(argv[4], &block_length) ||
      !parse_u64(argv[5], &seed) || !parse_u64(argv[6], &stream) ||
      !parse_size(argv[7], &restart_split) ||
      (sample_count != 0 &&
       (block_length == 0 || sample_count % block_length != 0 ||
        restart_split == 0 || restart_split >= burn_in + sample_count))) {
    if (world_rank == 0) {
      fprintf(stderr,
              "usage: %s fixture sample_count burn_in block_length"
              " seed stream restart_split\n",
              argv[0]);
    }
    local_ok = 0;
    goto finish;
  }
  if (world_rank == 0) {
    printf("RUN schema=%d fixture=%s rank_count=%d sample_count=%zu"
           " burn_in=%zu block_length=%zu seed=%" PRIu64
           " stream=%" PRIu64 " restart_split=%zu\n",
           DRIVER_SCHEMA, fixture.id, world_size, sample_count, burn_in,
           block_length, seed, stream, restart_split);
  }
  if (!begin_workspace(&fixture, &anchor_plan, &anchor_workspace) ||
      !print_anchors(&fixture, anchor_workspace)) {
    local_ok = 0;
    goto finish;
  }
  end_workspace(anchor_plan, anchor_workspace);
  anchor_plan = NULL;
  anchor_workspace = NULL;
  if (sample_count == 0) {
    if (world_rank == 0) puts("DECISION deterministic=PASS stochastic=NA");
    goto finish;
  }
  if (!run_trace(&fixture, burn_in, sample_count, restart_split, seed,
                 stream, 0, &full) ||
      !run_trace(&fixture, burn_in, sample_count, restart_split, seed,
                 stream, 1, &restarted) ||
      !traces_equal(&full, &restarted) ||
      !trace_rank_invariant(&full) ||
      !print_trace(&full, burn_in, block_length, 1)) {
    local_ok = 0;
    goto finish;
  }
  if (world_rank == 0) puts("DECISION deterministic=PASS stochastic=RECORDED");
finish:
  if (anchor_workspace != NULL || anchor_plan != NULL) {
    end_workspace(anchor_plan, anchor_workspace);
  }
  trace_destroy(&full);
  trace_destroy(&restarted);
#ifdef _mpi_use
  MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Finalize();
#else
  global_ok = local_ok;
#endif
  return global_ok ? 0 : 1;
}
