#include "power_lanczos_corrected_dispatch.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "corrected dispatch requires the bounded power-Lanczos engine"
#endif

#include "power_lanczos_production_session.h"
#include "power_lanczos_stabilization_statistics.h"
#include "krylov_fock_proposal.h"

#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SCALE_PILOT_WARM_UP ((size_t)4096)
#define SCALE_PILOT_SAMPLE_COUNT ((size_t)4096)
#define DEFAULT_STAGE_WARM_UP ((size_t)4096)
#define DEFAULT_STAGE_SAMPLE_COUNT ((size_t)16384)
#define DEFAULT_STAGE_INTERVAL ((size_t)1)
#define STABILIZATION_BLOCK_COUNT ((size_t)32)
#define OUTPUT_CAPACITY ((size_t)262144)
#define CORRECTED_BOOTSTRAP_STREAM UINT64_C(0x503643424f4f5431)

static int checked_add_u64(uint64_t left, uint64_t right, uint64_t *sum) {
  if (sum == NULL || right > UINT64_MAX - left) return 0;
  *sum = left + right;
  return 1;
}

static uint64_t hash_u64(uint64_t hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    hash ^= (value >> (8 * byte)) & UINT64_C(0xff);
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

static MVMCKrylovStatus synchronize_status(
    MVMCKrylovStatus local,
    MVMCClassicPfaffianCommunicator communicator) {
#ifdef _mpi_use
  int encoded = (int)local;
  int reduced = 0;
  if (MPI_Allreduce(&encoded, &reduced, 1, MPI_INT, MPI_MAX,
                    communicator) != MPI_SUCCESS ||
      reduced < (int)MVMC_KRYLOV_STATUS_OK ||
      reduced > (int)MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  return (MVMCKrylovStatus)reduced;
#else
  (void)communicator;
  return local;
#endif
}

static size_t resolved_control(size_t value, size_t fallback) {
  return value == 0 ? fallback : value;
}

static MVMCKrylovStatus bootstrap_draw_bounded(
    void *context, size_t bound, size_t *value) {
  return mvmc_krylov_positive_sampler_rng_draw_bounded(
      (MVMCKrylovPositiveSamplerRng *)context, bound, value);
}

static MVMCKrylovStatus corrected_uniform_sector_bootstrap(
    const MVMCPowerLanczosCorrectedDispatchInput *input,
    const MVMCPowerLanczosClassicBridge *bridge,
    uint64_t **root_words, size_t *root_word_count,
    uint64_t *proposal_counter) {
  const MVMCKrylovFockModel *model;
  MVMCKrylovPositiveSamplerRng rng;
  MVMCKrylovFockUniformProposalResult proposal;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  size_t word_count;
  uint64_t *words = NULL;
  if (input == NULL || bridge == NULL || root_words == NULL ||
      root_word_count == NULL || proposal_counter == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *root_words = NULL;
  *root_word_count = 0;
  *proposal_counter = 0;
  model = mvmc_power_lanczos_classic_bridge_model(bridge);
  if (model == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  word_count = mvmc_krylov_fock_word_count(model->site_count);
  if (word_count == 0 || word_count > SIZE_MAX / sizeof(*words)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  if (input->mpi_world_rank == 0) {
    words = (uint64_t *)calloc(word_count, sizeof(*words));
    if (words == NULL) {
      status = MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_krylov_positive_sampler_rng_seed(
          input->resolved_base_seed ^ input->base_generation,
          CORRECTED_BOOTSTRAP_STREAM ^ input->run_index, &rng);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      memset(&proposal, 0, sizeof(proposal));
      status = mvmc_krylov_fock_proposal_draw_uniform_sector(
          model, bootstrap_draw_bounded, &rng, words, word_count, &proposal);
      if (status == MVMC_KRYLOV_STATUS_OK &&
          (!proposal.valid || proposal.status != MVMC_KRYLOV_STATUS_OK ||
           proposal.word_count != word_count)) {
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
      }
    }
  }
  status = synchronize_status(status, input->world_communicator);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    free(words);
    return status;
  }
  if (input->mpi_world_rank == 0) {
    *root_words = words;
    *root_word_count = word_count;
  }
  /* One corrected bootstrap proposal was completed.  The production session
   * subsequently evaluates this sector-valid state with the bounded absolute
   * evaluator and its strictly positive guide floor. */
  *proposal_counter = UINT64_C(1);
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus validate_controls(
    const MVMCPowerLanczosCorrectedChainControls *controls,
    size_t *coefficient_warm_up, size_t *coefficient_sample_count,
    size_t *coefficient_interval, size_t *final_warm_up,
    size_t *final_sample_count, size_t *final_interval) {
  if (controls == NULL || coefficient_warm_up == NULL ||
      coefficient_sample_count == NULL || coefficient_interval == NULL ||
      final_warm_up == NULL || final_sample_count == NULL ||
      final_interval == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  *coefficient_warm_up =
      resolved_control(controls->coefficient_warm_up, DEFAULT_STAGE_WARM_UP);
  *coefficient_sample_count = resolved_control(
      controls->coefficient_sample_count, DEFAULT_STAGE_SAMPLE_COUNT);
  *coefficient_interval =
      resolved_control(controls->coefficient_interval, DEFAULT_STAGE_INTERVAL);
  *final_warm_up =
      resolved_control(controls->final_warm_up, DEFAULT_STAGE_WARM_UP);
  *final_sample_count =
      resolved_control(controls->final_sample_count, DEFAULT_STAGE_SAMPLE_COUNT);
  *final_interval =
      resolved_control(controls->final_interval, DEFAULT_STAGE_INTERVAL);
  if (*coefficient_sample_count < STABILIZATION_BLOCK_COUNT ||
      *final_sample_count < STABILIZATION_BLOCK_COUNT ||
      *coefficient_sample_count % STABILIZATION_BLOCK_COUNT != 0 ||
      *final_sample_count % STABILIZATION_BLOCK_COUNT != 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus create_limits(
    const MVMCKrylovFockModel *model, int power_step,
    uint64_t generation_hash, MVMCKrylovBoundedLimits *limits) {
  uint64_t frontier = 1;
  uint64_t nodes = 1;
  uint64_t transitions = 0;
  int depth;
  if (model == NULL || limits == NULL || model->term_count == 0 ||
      power_step < 1 || power_step > 2 ||
      model->term_count > SIZE_MAX / sizeof(MVMCKrylovHamiltonianTerm)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (depth = 0; depth < power_step + 1; ++depth) {
    uint64_t next;
    if (model->term_count > UINT64_MAX / frontier) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    next = frontier * (uint64_t)model->term_count;
    if (!checked_add_u64(transitions, next, &transitions) ||
        !checked_add_u64(nodes, next, &nodes)) {
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    frontier = next;
  }
  if (nodes > UINT64_C(10000000) ||
      transitions > UINT64_C(100000000)) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  memset(limits, 0, sizeof(*limits));
  limits->amplitude_policy_hash = generation_hash;
  /* This cache is an optimization only.  Keep the corrected route's fixed
   * per-rank footprint small; the bounded engine remains correct on misses. */
  limits->cache_bytes = (size_t)4 * 1024 * 1024;
  limits->max_row_transitions = model->term_count;
  limits->max_workspace_bytes = (size_t)512 * 1024 * 1024;
  limits->max_node_expansions = nodes;
  limits->max_terminal_amplitude_calls = frontier;
  limits->max_total_row_transitions = transitions;
  limits->max_order = power_step + 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus build_execution(
    const MVMCPowerLanczosCorrectedDispatchInput *input,
    const MVMCPowerLanczosClassicBridge *bridge,
    MVMCPowerLanczosProductionExecution *execution) {
  const MVMCKrylovFockModel *model =
      mvmc_power_lanczos_classic_bridge_model(bridge);
  const size_t basis_count = (size_t)input->power_step + 1;
  uint64_t policy_hash = UINT64_C(1469598103934665603);
  size_t index;
  MVMCKrylovStatus status;
  memset(execution, 0, sizeof(*execution));
  execution->world_communicator = input->world_communicator;
  execution->chain_communicator = input->chain_communicator;
  execution->proposal_model = model;
  execution->amplitude =
      mvmc_power_lanczos_classic_bridge_amplitude(bridge);
  execution->amplitude_context =
      mvmc_power_lanczos_classic_bridge_amplitude_context(bridge);
  execution->amplitude_generation_hash =
      mvmc_power_lanczos_classic_bridge_generation_hash(bridge);
  status = create_limits(model, input->power_step,
                         execution->amplitude_generation_hash,
                         &execution->bounded_limits);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  execution->coefficient_guide_policy.order = (int)basis_count;
  execution->coefficient_guide_policy.eta = 0x1p-40;
  execution->matrix_policy.order = (int)basis_count;
  execution->matrix_policy.eta = 0x1p-40;
  for (index = 0; index <= basis_count; ++index) {
    const double reference_weight = index == 0 ? 0x1p-40 : 1.0;
    execution->coefficient_guide_policy.lambda[index] = reference_weight;
    execution->matrix_policy.guide_lambda[index] = reference_weight;
    execution->matrix_policy.target_weight[index] = 1.0;
  }
  policy_hash = hash_u64(policy_hash, UINT64_C(0x53544142494c495a));
  policy_hash = hash_u64(policy_hash, (uint64_t)(unsigned)input->power_step);
  policy_hash = hash_u64(policy_hash,
                         execution->amplitude_generation_hash);
  execution->coefficient_guide_policy.policy_hash =
      policy_hash == 0 ? UINT64_C(1) : policy_hash;
  status = mvmc_krylov_positive_sampler_proposal_policy_create(
      1, 16, &execution->proposal_policy);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (mvmc_power_lanczos_gevp_default_policy(
          0x1p-40, &execution->gevp_policy) != MVMC_KRYLOV_GEVP_OK) {
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  execution->block_count = STABILIZATION_BLOCK_COUNT;
  execution->scale_pilot_enabled = 1;
  execution->scale_pilot_warm_up = SCALE_PILOT_WARM_UP;
  execution->scale_pilot_sample_count = SCALE_PILOT_SAMPLE_COUNT;
  execution->eta_relative = 0x1p-40;
  status = validate_controls(
      &input->controls, &execution->coefficient_warm_up,
      &execution->coefficient_sample_count,
      &execution->coefficient_interval, &execution->final_warm_up,
      &execution->final_sample_count, &execution->final_interval);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  execution->maximum_leave_one_projective_distance = 1.0;
  execution->observable_enabled = 0;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus terminal_sampling_support(
    const MVMCPowerLanczosCorrectedDispatchInput *input,
    const MVMCPowerLanczosClassicBridge *bridge, double *log_support) {
  MVMCKrylovScaledAmplitudeResult amplitude;
  uint64_t *words = NULL;
  const size_t expected_words =
      mvmc_krylov_fock_word_count((size_t)input->classic_view->site_count);
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  if (log_support == NULL || expected_words == 0 ||
      expected_words > SIZE_MAX / sizeof(*words)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  words = (uint64_t *)calloc(expected_words, sizeof(*words));
  if (words == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  if (input->mpi_world_rank == 0) {
    if (input->root_terminal_words == NULL ||
        input->root_terminal_word_count != expected_words) {
      status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    } else {
      memcpy(words, input->root_terminal_words,
             expected_words * sizeof(*words));
    }
  }
  status = synchronize_status(status, input->world_communicator);
#ifdef _mpi_use
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (expected_words > (size_t)INT_MAX ||
       MPI_Bcast(words, (int)expected_words, MPI_UINT64_T, 0,
                 input->world_communicator) != MPI_SUCCESS)) {
    status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#endif
  memset(&amplitude, 0, sizeof(amplitude));
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_classic_bridge_amplitude(bridge)(
        words, expected_words,
        mvmc_power_lanczos_classic_bridge_amplitude_context(bridge),
        &amplitude);
  }
  status = synchronize_status(status, input->world_communicator);
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (!mvmc_scaled_complex_is_valid(&amplitude.value) ||
       amplitude.value.state != MVMC_SCALED_COMPLEX_FINITE_NONZERO ||
       !isfinite(amplitude.value.log_abs))) {
    status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  status = synchronize_status(status, input->world_communicator);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    *log_support = 2.0 * amplitude.value.log_abs;
    if (!isfinite(*log_support)) status = MVMC_KRYLOV_STATUS_NONFINITE;
  }
  free(words);
  return synchronize_status(status, input->world_communicator);
}

static MVMCKrylovStatus root_render_and_commit(
    const MVMCPowerLanczosCorrectedDispatchInput *input,
    const MVMCPowerLanczosProductionSession *session,
    MVMCPowerLanczosCorrectedDispatchResult *result) {
  MVMCPowerLanczosStabilizationResult statistics;
  MVMCPowerLanczosOutputFile file;
  MVMCPowerLanczosOutputTransactionConfig transaction;
  MVMCPowerLanczosOutputTransactionResult committed;
  char *output = NULL;
  size_t output_size = 0;
  MVMCKrylovStatus status;
  memset(&statistics, 0, sizeof(statistics));
  status = mvmc_power_lanczos_stabilization_statistics_evaluate(
      session, 0, &statistics);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  output = (char *)malloc(OUTPUT_CAPACITY);
  if (output == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  status = mvmc_power_lanczos_stabilization_output_render(
      session, &statistics, input->root_identity, output, OUTPUT_CAPACITY,
      &output_size);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    memset(&file, 0, sizeof(file));
    file.target_basename = input->root_output_basename;
    file.contents = (const unsigned char *)output;
    file.content_size = output_size;
    memset(&transaction, 0, sizeof(transaction));
    transaction.directory = input->root_output_directory;
    transaction.files = &file;
    transaction.file_count = 1;
    transaction.manifest_index = 0;
    result->output_status = mvmc_power_lanczos_output_transaction_commit(
        &transaction, &committed);
    if (result->output_status !=
            MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_OK ||
        !committed.valid || !committed.committed ||
        !committed.rollback_complete) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    } else {
      result->decision = statistics.decision;
      result->failed_gates = statistics.failed_gates;
      result->energy = statistics.energy;
      result->energy_standard_error = statistics.energy_standard_error;
      result->variance = statistics.variance;
      result->variance_standard_error = statistics.variance_standard_error;
      result->output_size = output_size;
      memcpy(result->output_sha256, committed.sha256[0],
             sizeof(result->output_sha256));
    }
  }
  free(output);
  return status;
}

MVMCKrylovStatus mvmc_power_lanczos_corrected_dispatch_run(
    const MVMCPowerLanczosCorrectedDispatchInput *input,
    MVMCPowerLanczosCorrectedDispatchResult *result) {
  MVMCPowerLanczosClassicBridge *bridge = NULL;
  MVMCPowerLanczosProductionSession *session = NULL;
  MVMCPowerLanczosProductionSessionConfig config;
  MVMCPowerLanczosProductionExecution execution;
  MVMCPowerLanczosTerminalSnapshotInput terminal;
  const MVMCPowerLanczosTerminalSnapshotInput *root_terminal = NULL;
  uint64_t *bootstrap_words = NULL;
  const uint64_t *selected_terminal_words = NULL;
  size_t selected_terminal_word_count = 0;
  uint64_t selected_terminal_proposal_counter = 0;
  double log_sampling_support = NAN;
  const char *stage = "input-validation";
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  int root_status;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(result, 0, sizeof(*result));
  result->version = MVMC_POWER_LANCZOS_CORRECTED_DISPATCH_VERSION;
  result->bootstrap_mode = input == NULL
                               ? MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_PROVIDED_TERMINAL
                               : input->bootstrap_mode;
  result->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  result->output_status =
      MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_INVALID_ARGUMENT;
  if (input == NULL || input->classic_view == NULL ||
      (input->power_step != 1 && input->power_step != 2) ||
      input->resolved_base_seed == 0 || input->base_generation == 0 ||
      (input->bootstrap_mode !=
           MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_PROVIDED_TERMINAL &&
       input->bootstrap_mode !=
           MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_UNIFORM_SECTOR) ||
      (input->bootstrap_mode ==
           MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_PROVIDED_TERMINAL &&
       input->terminal_proposal_counter == 0) ||
      (input->bootstrap_mode ==
           MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_UNIFORM_SECTOR &&
       input->terminal_proposal_counter != 0) ||
      input->mpi_world_size == 0 ||
      input->mpi_world_rank >= input->mpi_world_size ||
      input->split_size == 0 ||
      input->split_size > input->mpi_world_size) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (input->mpi_world_rank == 0 &&
      (input->root_identity == NULL ||
       input->root_output_directory == NULL ||
       !mvmc_power_lanczos_output_basename_valid(
           input->root_output_basename) ||
       (input->bootstrap_mode ==
            MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_PROVIDED_TERMINAL &&
        (input->root_terminal_words == NULL ||
         input->root_terminal_word_count == 0)) ||
       (input->bootstrap_mode ==
            MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_UNIFORM_SECTOR &&
        (input->root_terminal_words != NULL ||
         input->root_terminal_word_count != 0)))) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = synchronize_status(status, input->world_communicator);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    result->status = status;
    return status;
  }
  stage = "classic-bridge";
  status = mvmc_power_lanczos_classic_bridge_create(
      input->classic_view, input->world_communicator,
      input->chain_communicator, &bridge);
  if (status == MVMC_KRYLOV_STATUS_OK &&
      input->bootstrap_mode ==
          MVMC_POWER_LANCZOS_CORRECTED_BOOTSTRAP_PROVIDED_TERMINAL) {
    stage = "terminal-support";
    status = terminal_sampling_support(input, bridge, &log_sampling_support);
    selected_terminal_words = input->root_terminal_words;
    selected_terminal_word_count = input->root_terminal_word_count;
    selected_terminal_proposal_counter = input->terminal_proposal_counter;
  } else if (status == MVMC_KRYLOV_STATUS_OK) {
    stage = "corrected-bootstrap";
    status = corrected_uniform_sector_bootstrap(
        input, bridge, &bootstrap_words, &selected_terminal_word_count,
        &selected_terminal_proposal_counter);
    selected_terminal_words = bootstrap_words;
    log_sampling_support = 0.0;
  }
  result->bootstrap_proposal_counter = selected_terminal_proposal_counter;
  if (status == MVMC_KRYLOV_STATUS_OK) {
    stage = "session-create";
    memset(&config, 0, sizeof(config));
    config.power_step = input->power_step;
    config.resolved_base_seed = input->resolved_base_seed;
    config.run_index = input->run_index;
    config.mpi_world_rank = input->mpi_world_rank;
    config.mpi_world_size = input->mpi_world_size;
    config.split_size = input->split_size;
    config.base_generation = input->base_generation;
    config.configuration_bit_count =
        (size_t)input->classic_view->site_count * 2;
    memset(&terminal, 0, sizeof(terminal));
    if (input->mpi_world_rank == 0) {
      terminal.canonical_words = selected_terminal_words;
      terminal.word_count = selected_terminal_word_count;
      terminal.configuration_bit_count = config.configuration_bit_count;
      terminal.generation = input->base_generation;
      terminal.terminal_proposal_counter =
          selected_terminal_proposal_counter;
      terminal.expected_particle_count =
          (size_t)input->classic_view->up_electron_count +
          (size_t)input->classic_view->down_electron_count;
      terminal.log_sampling_support = log_sampling_support;
      root_terminal = &terminal;
    }
    status = mvmc_power_lanczos_production_session_create(
        &config, input->world_communicator, root_terminal, &session);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    stage = "session-verify";
    status = mvmc_power_lanczos_production_session_verify(session);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    stage = "execution-build";
    status = build_execution(input, bridge, &execution);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    stage = "execution-prepare";
    status = mvmc_power_lanczos_production_session_prepare_execution(
        session, &execution);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    stage = "coefficient-run";
    status = mvmc_power_lanczos_production_session_run_coefficient(session);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    stage = "coefficient-freeze";
    status = mvmc_power_lanczos_production_session_freeze_coefficient(session);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    stage = "final-run";
    status = mvmc_power_lanczos_production_session_run_final(session);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    stage = "session-finalize";
    status = mvmc_power_lanczos_production_session_finalize(session);
  }
  root_status = (int)status;
  if (status == MVMC_KRYLOV_STATUS_OK && input->mpi_world_rank == 0) {
    stage = "output-commit";
    status = root_render_and_commit(input, session, result);
    root_status = (int)status;
  }
#ifdef _mpi_use
  if (MPI_Bcast(&root_status, 1, MPI_INT, 0,
                input->world_communicator) != MPI_SUCCESS) {
    status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  } else if (root_status < (int)MVMC_KRYLOV_STATUS_OK ||
             root_status >
                 (int)MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE) {
    status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  } else {
    status = (MVMCKrylovStatus)root_status;
  }
#else
  status = (MVMCKrylovStatus)root_status;
#endif
  if (status == MVMC_KRYLOV_STATUS_OK) {
#ifdef _mpi_use
    unsigned long long packed[4];
    double values[4];
    if (input->mpi_world_rank == 0) {
      packed[0] = (unsigned long long)result->decision;
      packed[1] = (unsigned long long)result->failed_gates;
      packed[2] = (unsigned long long)result->output_size;
      packed[3] = (unsigned long long)result->output_status;
      values[0] = result->energy;
      values[1] = result->energy_standard_error;
      values[2] = result->variance;
      values[3] = result->variance_standard_error;
    }
    if (MPI_Bcast(packed, 4, MPI_UNSIGNED_LONG_LONG, 0,
                  input->world_communicator) != MPI_SUCCESS ||
        MPI_Bcast(values, 4, MPI_DOUBLE, 0,
                  input->world_communicator) != MPI_SUCCESS ||
        MPI_Bcast(result->output_sha256,
                  MVMC_POWER_LANCZOS_SHA256_HEX_CAPACITY, MPI_CHAR, 0,
                  input->world_communicator) != MPI_SUCCESS) {
      status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    } else if (input->mpi_world_rank != 0) {
      result->decision = (MVMCPowerLanczosStabilizationDecision)packed[0];
      result->failed_gates = (uint64_t)packed[1];
      result->output_size = (size_t)packed[2];
      result->output_status =
          (MVMCPowerLanczosOutputTransactionStatus)packed[3];
      result->energy = values[0];
      result->energy_standard_error = values[1];
      result->variance = values[2];
      result->variance_standard_error = values[3];
    }
#endif
  }
  mvmc_power_lanczos_production_session_destroy(session);
  mvmc_power_lanczos_classic_bridge_destroy(bridge);
  free(bootstrap_words);
  if (status != MVMC_KRYLOV_STATUS_OK && input->mpi_world_rank == 0) {
    fprintf(stderr, "P6 corrected dispatch stage %s failed: %s\n", stage,
            mvmc_krylov_status_string(status));
  }
  result->status = status;
  result->valid = status == MVMC_KRYLOV_STATUS_OK;
  return status;
}

#undef SCALE_PILOT_WARM_UP
#undef SCALE_PILOT_SAMPLE_COUNT
#undef DEFAULT_STAGE_WARM_UP
#undef DEFAULT_STAGE_SAMPLE_COUNT
#undef DEFAULT_STAGE_INTERVAL
#undef STABILIZATION_BLOCK_COUNT
#undef OUTPUT_CAPACITY
#undef CORRECTED_BOOTSTRAP_STREAM
