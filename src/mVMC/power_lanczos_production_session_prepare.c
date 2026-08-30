#include "power_lanczos_production_session_internal.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "power_lanczos_production_session_prepare.c requires bounded Krylov"
#endif

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

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

static int checked_proposal_count(size_t warm_up, size_t sample_count,
                                  size_t interval) {
  uint64_t samples;
  uint64_t stride;
  uint64_t warm;
#if SIZE_MAX > UINT64_MAX
  if (warm_up > UINT64_MAX || sample_count > UINT64_MAX ||
      interval > UINT64_MAX) {
    return 0;
  }
#endif
  samples = (uint64_t)sample_count;
  stride = (uint64_t)interval;
  warm = (uint64_t)warm_up;
  return stride != 0 &&
         (samples == 0 || stride <= (UINT64_MAX - warm) / samples);
}

static uint64_t double_bits(double value) {
  uint64_t bits = 0;
  memcpy(&bits, &value, sizeof(bits));
  return bits;
}

static void hash_u64(uint64_t *hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    *hash ^= (value >> (8 * byte)) & UINT64_C(0xff);
    *hash *= UINT64_C(1099511628211);
  }
}

static void hash_bytes(uint64_t *hash, const unsigned char *bytes,
                       size_t count) {
  size_t index;
  for (index = 0; index < count; ++index) {
    hash_u64(hash, (uint64_t)bytes[index]);
  }
}

static MVMCKrylovStatus confirmation_status(
    MVMCPowerLanczosObservableConfirmationStatus status) {
  switch (status) {
    case MVMC_POWER_LANCZOS_CONFIRMATION_OK:
      return MVMC_KRYLOV_STATUS_OK;
    case MVMC_POWER_LANCZOS_CONFIRMATION_RESOURCE_LIMIT:
      return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    case MVMC_POWER_LANCZOS_CONFIRMATION_ALLOCATION_FAILURE:
      return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    case MVMC_POWER_LANCZOS_CONFIRMATION_NONFINITE:
      return MVMC_KRYLOV_STATUS_NONFINITE;
    case MVMC_POWER_LANCZOS_CONFIRMATION_GEVP_FAILURE:
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    default:
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
}

static int valid_status(MVMCKrylovStatus status) {
  return status >= MVMC_KRYLOV_STATUS_OK &&
         status <= MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

static MVMCKrylovStatus synchronize_world_status(
    MVMCPowerLanczosProductionSession *session,
    MVMCKrylovStatus local_status) {
  MVMCKrylovStatus effective = valid_status(local_status)
                                      ? local_status
                                      : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
#ifdef _mpi_use
  struct {
    int value;
    int index;
  } local, global;
  local.value = (int)effective;
  local.index = (int)session->config.mpi_world_rank;
  if (MPI_Allreduce(&local, &global, 1, MPI_2INT, MPI_MAXLOC,
                    session->world_communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  effective = (MVMCKrylovStatus)global.value;
#else
  (void)session;
#endif
  return effective;
}

static int valid_guide_policy(const MVMCKrylovPositiveGuidePolicy *policy,
                              size_t basis_count) {
  size_t index;
  if (policy == NULL || policy->order != (int)basis_count ||
      !isfinite(policy->eta) || policy->eta <= 0.0 ||
      policy->policy_hash == 0) {
    return 0;
  }
  for (index = 0; index <= basis_count; ++index) {
    if (!isfinite(policy->lambda[index]) || policy->lambda[index] < 0.0 ||
        !isfinite(policy->log_basis_scale[index])) {
      return 0;
    }
  }
  return 1;
}

static int valid_matrix_policy(
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    size_t basis_count) {
  size_t index;
  int has_target = 0;
  if (policy == NULL || policy->order != (int)basis_count ||
      !isfinite(policy->eta) || policy->eta <= 0.0) {
    return 0;
  }
  for (index = 0; index <= basis_count; ++index) {
    if (!isfinite(policy->guide_lambda[index]) ||
        policy->guide_lambda[index] < 0.0 ||
        !isfinite(policy->target_weight[index]) ||
        policy->target_weight[index] < 0.0 ||
        !isfinite(policy->log_basis_scale[index])) {
      return 0;
    }
    has_target = has_target || policy->target_weight[index] > 0.0;
  }
  return has_target;
}

static int guide_policies_match(
    const MVMCKrylovPositiveGuidePolicy *guide,
    const MVMCKrylovMatrixMeasurementPolicy *matrix, size_t basis_count) {
  size_t index;
  if (guide == NULL || matrix == NULL || guide->eta != matrix->eta) {
    return 0;
  }
  for (index = 0; index <= basis_count; ++index) {
    if (guide->lambda[index] != matrix->guide_lambda[index] ||
        guide->log_basis_scale[index] != matrix->log_basis_scale[index]) {
      return 0;
    }
  }
  return 1;
}

static int valid_scale_pilot_contract(
    const MVMCPowerLanczosProductionExecution *execution,
    size_t basis_count) {
  size_t index;
  if (execution == NULL) return 0;
  if (!execution->scale_pilot_enabled) {
    return execution->scale_pilot_warm_up == 0 &&
           execution->scale_pilot_sample_count == 0 &&
           execution->eta_relative == 0.0;
  }
  if (execution->scale_pilot_enabled != 1 ||
      execution->scale_pilot_warm_up == 0 ||
      execution->scale_pilot_sample_count == 0 ||
      !isfinite(execution->eta_relative) ||
      execution->eta_relative <= 0.0 ||
      execution->coefficient_guide_policy.eta != 0x1p-40 ||
      execution->matrix_policy.eta != 0x1p-40 ||
      !checked_proposal_count(execution->scale_pilot_warm_up,
                              execution->scale_pilot_sample_count, 1)) {
    return 0;
  }
  for (index = 0; index <= basis_count; ++index) {
    const double reference_weight = index == 0 ? 0x1p-40 : 1.0;
    if (execution->coefficient_guide_policy.lambda[index] !=
            reference_weight ||
        execution->matrix_policy.guide_lambda[index] != reference_weight ||
        execution->matrix_policy.target_weight[index] != 1.0 ||
        execution->coefficient_guide_policy.log_basis_scale[index] != 0.0 ||
        execution->matrix_policy.log_basis_scale[index] != 0.0) {
      return 0;
    }
  }
  return 1;
}

static int valid_proposal_policy(
    const MVMCKrylovPositiveSamplerProposalPolicy *policy) {
  MVMCKrylovPositiveSamplerProposalPolicy expected;
  MVMCKrylovStatus status;
  if (policy == NULL || !policy->valid ||
      policy->status != MVMC_KRYLOV_STATUS_OK ||
      policy->version != MVMC_KRYLOV_PROPOSAL_POLICY_VERSION ||
      policy->policy_hash == 0 ||
      policy->subset_order != MVMC_KRYLOV_PROPOSAL_SUBSET_CANONICAL_V1 ||
      policy->self_report_rule != MVMC_KRYLOV_PROPOSAL_SELF_REPORT_V1 ||
      policy->rng_algorithm != MVMC_KRYLOV_POSITIVE_SAMPLER_RNG_SPLITMIX64) {
    return 0;
  }
  if (policy->kernel_id ==
      MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_GLOBAL_MIXTURE) {
    status = mvmc_krylov_positive_sampler_proposal_policy_create(
        policy->global_numerator, policy->global_denominator, &expected);
    return status == MVMC_KRYLOV_STATUS_OK &&
           policy->neighbor_numerator == 0 &&
           policy->neighbor_denominator == 0 &&
           policy->distance_numerator == 0 &&
           policy->distance_denominator == 0 &&
           policy->distance_rounding_rule == 0 &&
           policy->policy_hash == expected.policy_hash;
  }
  if (policy->kernel_id ==
      MVMC_KRYLOV_PROPOSAL_KERNEL_NEIGHBOR_SHELL_MIXTURE) {
    status = mvmc_krylov_positive_sampler_shell_policy_create(
        policy->neighbor_numerator, policy->neighbor_denominator,
        policy->distance_numerator, policy->distance_denominator, &expected);
    return status == MVMC_KRYLOV_STATUS_OK &&
           policy->global_numerator == 0 &&
           policy->global_denominator == 0 &&
           policy->policy_hash == expected.policy_hash;
  }
  return 0;
}

static int valid_gevp_policy(const MVMCPowerLanczosGEVPPolicy *policy) {
  const double cutoff[4] = {0x1p-48, 0x1p-40, 0x1p-32, 0x1p-24};
  const double tolerance = 64.0 * DBL_EPSILON;
  return policy != NULL && policy->valid &&
         policy->policy_version == MVMC_POWER_LANCZOS_GEVP_POLICY_VERSION &&
         policy->cutoff_id >= MVMC_POWER_LANCZOS_GEVP_CUTOFF_S48 &&
         policy->cutoff_id <= MVMC_POWER_LANCZOS_GEVP_CUTOFF_S24 &&
         policy->rank_relative_cutoff == cutoff[policy->cutoff_id] &&
         policy->degenerate_root_gap_relative_threshold == tolerance &&
         policy->maximum_normwise_backward_error == tolerance &&
         policy->negative_variance_relative_tolerance == tolerance;
}

static int valid_model_shape(
    const MVMCPowerLanczosProductionSession *session,
    const MVMCKrylovFockModel *model) {
  size_t index;
  size_t bit_count;
  if (model == NULL || model->site_count == 0 ||
      model->site_count > SIZE_MAX / 2 ||
      !checked_multiply(model->site_count, 2, &bit_count) ||
      bit_count != session->config.configuration_bit_count ||
      (model->term_count != 0 && model->terms == NULL) ||
      (model->operator_count != 0 && model->operators == NULL)) {
    return 0;
  }
  for (index = 0; index < model->term_count; ++index) {
    const MVMCKrylovHamiltonianTerm *term = &model->terms[index];
    if (term->operator_offset > model->operator_count ||
        term->operator_count > model->operator_count - term->operator_offset ||
        !isfinite(creal(term->coefficient)) ||
        !isfinite(cimag(term->coefficient))) {
      return 0;
    }
  }
  return 1;
}

static int valid_execution_shape(
    const MVMCPowerLanczosProductionSession *session,
    const MVMCPowerLanczosProductionExecution *execution,
    size_t *maximum_entries_per_block) {
  size_t coefficient_entries;
  size_t final_entries = 2;
  size_t observable_entries;
  size_t request_count = 0;
  size_t coefficient_block_length;
  size_t final_block_length;
  size_t global_block_length;
  if (execution == NULL || maximum_entries_per_block == NULL ||
      execution->amplitude == NULL || execution->amplitude_generation_hash == 0 ||
      !valid_model_shape(session, execution->proposal_model) ||
      execution->bounded_limits.max_order != (int)session->basis_count ||
      !valid_guide_policy(&execution->coefficient_guide_policy,
                          session->basis_count) ||
      !valid_proposal_policy(&execution->proposal_policy) ||
      !valid_matrix_policy(&execution->matrix_policy,
                           session->basis_count) ||
      !guide_policies_match(&execution->coefficient_guide_policy,
                            &execution->matrix_policy,
                            session->basis_count) ||
      !valid_scale_pilot_contract(execution, session->basis_count) ||
      !valid_gevp_policy(&execution->gevp_policy) ||
      execution->block_count < 2 || execution->coefficient_sample_count == 0 ||
      execution->final_sample_count == 0 || execution->coefficient_interval == 0 ||
      execution->final_interval == 0 ||
      execution->coefficient_sample_count % execution->block_count != 0 ||
      execution->final_sample_count % execution->block_count != 0 ||
      execution->coefficient_sample_count >
          MVMC_POWER_LANCZOS_BLOCK_MAX_EXACT_SAMPLE_COUNT ||
      execution->final_sample_count >
          MVMC_POWER_LANCZOS_BLOCK_MAX_EXACT_SAMPLE_COUNT ||
      !isfinite(execution->maximum_leave_one_projective_distance) ||
      execution->maximum_leave_one_projective_distance < 0.0 ||
      !checked_proposal_count(execution->coefficient_warm_up,
                              execution->coefficient_sample_count,
                              execution->coefficient_interval) ||
      !checked_proposal_count(execution->final_warm_up,
                              execution->final_sample_count,
                              execution->final_interval) ||
      !checked_multiply(session->upper_count, 4, &coefficient_entries)) {
    return 0;
  }
  coefficient_block_length =
      execution->coefficient_sample_count / execution->block_count;
  final_block_length = execution->final_sample_count / execution->block_count;
  if (coefficient_block_length == 0 || final_block_length == 0 ||
      session->chain_count == 0 ||
      !checked_multiply(execution->coefficient_sample_count,
                        session->chain_count, &global_block_length) ||
      global_block_length >
          (size_t)MVMC_POWER_LANCZOS_BLOCK_MAX_EXACT_SAMPLE_COUNT ||
      !checked_multiply(execution->final_sample_count,
                        session->chain_count, &global_block_length) ||
      global_block_length >
          (size_t)MVMC_POWER_LANCZOS_BLOCK_MAX_EXACT_SAMPLE_COUNT ||
      !checked_multiply(coefficient_block_length, session->chain_count,
                        &global_block_length) ||
      global_block_length >
          (size_t)MVMC_POWER_LANCZOS_BLOCK_MAX_EXACT_SAMPLE_COUNT ||
      !checked_multiply(final_block_length, session->chain_count,
                        &global_block_length) ||
      global_block_length >
          (size_t)MVMC_POWER_LANCZOS_BLOCK_MAX_EXACT_SAMPLE_COUNT) {
    return 0;
  }
  if (execution->observable_enabled) {
    const MVMCPowerLanczosObservablePlan *plan = execution->observable_plan;
    const MVMCPowerLanczosObservableLayout *layout =
        &execution->observable_layout;
    if (session->config.power_step != 1 || plan == NULL ||
        plan->record_count <= 0 || plan->records == NULL ||
        execution->block_count <
            MVMC_POWER_LANCZOS_CONFIRMATION_MIN_BLOCK_COUNT ||
        coefficient_block_length != final_block_length ||
        plan->nsite != layout->nsite ||
        layout->nsite <= 0 || (size_t)layout->nsite !=
                                  execution->proposal_model->site_count ||
        layout->up_electron_count < 0 ||
        (size_t)layout->up_electron_count !=
            execution->proposal_model->up_electron_count ||
        layout->down_electron_count < 0 ||
        (size_t)layout->down_electron_count !=
            execution->proposal_model->down_electron_count ||
        layout->pure_spin != execution->proposal_model->pure_spin ||
        mvmc_power_lanczos_observable_word_count(layout) == 0 ||
        mvmc_power_lanczos_observable_plan_rehash(plan, NULL, 0) !=
            MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK ||
        mvmc_power_lanczos_observable_plan_semantic_rehash(plan, NULL, 0) !=
            MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) {
      return 0;
    }
    request_count = (size_t)plan->record_count;
    if (!checked_multiply(request_count,
                          MVMC_POWER_LANCZOS_OBSERVABLE_MATRIX_ENTRIES,
                          &observable_entries) ||
        !checked_add(coefficient_entries, observable_entries,
                     &coefficient_entries) ||
        !checked_add(final_entries, request_count, &final_entries)) {
      return 0;
    }
  } else if (execution->observable_plan != NULL) {
    return 0;
  }
  *maximum_entries_per_block =
      coefficient_entries > final_entries ? coefficient_entries
                                           : final_entries;
  return 1;
}

static uint64_t execution_fingerprint(
    const MVMCPowerLanczosProductionSession *session,
    const MVMCPowerLanczosProductionExecution *execution) {
  const MVMCKrylovFockModel *model = execution->proposal_model;
  uint64_t hash = UINT64_C(1469598103934665603);
  size_t index;
  hash_u64(&hash, MVMC_POWER_LANCZOS_PRODUCTION_SESSION_VERSION);
  hash_u64(&hash, (uint64_t)(unsigned int)session->config.power_step);
  hash_u64(&hash, session->config.resolved_base_seed);
  hash_u64(&hash, session->config.run_index);
  hash_u64(&hash, (uint64_t)session->config.mpi_world_size);
  hash_u64(&hash, (uint64_t)session->config.split_size);
  hash_u64(&hash, session->config.base_generation);
  hash_u64(&hash, (uint64_t)session->config.configuration_bit_count);
  hash_u64(&hash, session->snapshot_selection.version);
  hash_u64(&hash,
           (uint64_t)(unsigned int)session->snapshot_selection.selection_kind);
  hash_u64(&hash,
           (uint64_t)session->snapshot_selection.owner_world_rank);
  hash_u64(&hash, session->snapshot_selection.terminal_proposal_counter);
  hash_u64(&hash, session->snapshot_selection.generation);
  hash_u64(&hash, (uint64_t)session->snapshot_selection.word_count);
  hash_u64(&hash,
           (uint64_t)session->snapshot_selection.configuration_bit_count);
  hash_u64(&hash,
           (uint64_t)session->snapshot_selection.expected_particle_count);
  hash_bytes(&hash, session->snapshot_selection.snapshot_sha256,
             MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES);
  hash_u64(&hash, execution->amplitude_generation_hash);
  hash_u64(&hash, (uint64_t)model->site_count);
  hash_u64(&hash, (uint64_t)model->up_electron_count);
  hash_u64(&hash, (uint64_t)model->down_electron_count);
  hash_u64(&hash, (uint64_t)(unsigned int)model->pure_spin);
  hash_u64(&hash, (uint64_t)(unsigned int)model->hermitian);
  hash_u64(&hash, (uint64_t)model->term_count);
  hash_u64(&hash, (uint64_t)model->operator_count);
  for (index = 0; index < model->term_count; ++index) {
    const MVMCKrylovHamiltonianTerm *term = &model->terms[index];
    hash_u64(&hash, double_bits(creal(term->coefficient)));
    hash_u64(&hash, double_bits(cimag(term->coefficient)));
    hash_u64(&hash, (uint64_t)term->operator_offset);
    hash_u64(&hash, (uint64_t)term->operator_count);
    hash_u64(&hash, (uint64_t)(uint32_t)term->source_kind);
    hash_u64(&hash, (uint64_t)term->source_index);
  }
  for (index = 0; index < model->operator_count; ++index) {
    hash_u64(&hash, (uint64_t)(unsigned int)model->operators[index].kind);
    hash_u64(&hash, (uint64_t)model->operators[index].orbital);
  }
  hash_u64(&hash, execution->bounded_limits.amplitude_policy_hash);
  hash_u64(&hash, (uint64_t)execution->bounded_limits.cache_bytes);
  hash_u64(&hash, (uint64_t)execution->bounded_limits.max_row_transitions);
  hash_u64(&hash, (uint64_t)execution->bounded_limits.max_workspace_bytes);
  hash_u64(&hash, execution->bounded_limits.max_node_expansions);
  hash_u64(&hash, execution->bounded_limits.max_terminal_amplitude_calls);
  hash_u64(&hash, execution->bounded_limits.max_total_row_transitions);
  hash_u64(&hash, (uint64_t)(unsigned int)execution->bounded_limits.max_order);
  hash_u64(&hash, execution->coefficient_guide_policy.policy_hash);
  hash_u64(&hash, double_bits(execution->coefficient_guide_policy.eta));
  hash_u64(&hash, double_bits(execution->matrix_policy.eta));
  for (index = 0; index <= session->basis_count; ++index) {
    hash_u64(&hash,
             double_bits(execution->coefficient_guide_policy.lambda[index]));
    hash_u64(&hash, double_bits(
                        execution->coefficient_guide_policy.log_basis_scale[index]));
    hash_u64(&hash, double_bits(execution->matrix_policy.guide_lambda[index]));
    hash_u64(&hash, double_bits(execution->matrix_policy.target_weight[index]));
    hash_u64(&hash,
             double_bits(execution->matrix_policy.log_basis_scale[index]));
  }
  hash_u64(&hash, execution->proposal_policy.policy_hash);
  hash_u64(&hash, (uint64_t)execution->block_count);
  hash_u64(&hash, (uint64_t)(unsigned int)execution->scale_pilot_enabled);
  hash_u64(&hash, (uint64_t)execution->scale_pilot_warm_up);
  hash_u64(&hash, (uint64_t)execution->scale_pilot_sample_count);
  hash_u64(&hash, double_bits(execution->eta_relative));
  hash_u64(&hash, (uint64_t)execution->coefficient_warm_up);
  hash_u64(&hash, (uint64_t)execution->coefficient_sample_count);
  hash_u64(&hash, (uint64_t)execution->coefficient_interval);
  hash_u64(&hash, (uint64_t)execution->final_warm_up);
  hash_u64(&hash, (uint64_t)execution->final_sample_count);
  hash_u64(&hash, (uint64_t)execution->final_interval);
  hash_u64(&hash, (uint64_t)(unsigned int)execution->gevp_policy.cutoff_id);
  hash_u64(&hash, double_bits(execution->gevp_policy.rank_relative_cutoff));
  hash_u64(&hash, double_bits(
                      execution->maximum_leave_one_projective_distance));
  hash_u64(&hash, (uint64_t)(unsigned int)execution->observable_enabled);
  if (execution->observable_enabled) {
    const MVMCPowerLanczosObservablePlan *plan = execution->observable_plan;
    hash_u64(&hash, (uint64_t)(unsigned int)
                        execution->observable_layout.nsite);
    hash_u64(&hash, (uint64_t)(unsigned int)
                        execution->observable_layout.up_electron_count);
    hash_u64(&hash, (uint64_t)(unsigned int)
                        execution->observable_layout.down_electron_count);
    hash_u64(&hash, (uint64_t)(unsigned int)
                        execution->observable_layout.pure_spin);
    hash_u64(&hash, (uint64_t)(unsigned int)plan->record_count);
    hash_bytes(&hash,
               (const unsigned char *)plan->raw_observable_census_sha256,
               MVMC_POWER_LANCZOS_OBSERVABLE_SHA256_HEX_CAPACITY);
    hash_bytes(
        &hash,
        (const unsigned char *)plan->semantic_observable_census_sha256,
        MVMC_POWER_LANCZOS_OBSERVABLE_SHA256_HEX_CAPACITY);
  }
  return hash == 0 ? UINT64_C(1) : hash;
}

static void *allocate_owned(size_t count, size_t item_size,
                            size_t *owned_bytes,
                            MVMCKrylovStatus *status) {
  size_t bytes = 0;
  size_t next = 0;
  void *memory;
  if (*status != MVMC_KRYLOV_STATUS_OK || count == 0) return NULL;
  if (!checked_multiply(count, item_size, &bytes) ||
      !checked_add(*owned_bytes, bytes, &next)) {
    *status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    return NULL;
  }
  memory = calloc(count, item_size);
  if (memory == NULL) {
    *status = MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    return NULL;
  }
  *owned_bytes = next;
  return memory;
}

MVMCKrylovStatus mvmc_power_lanczos_production_session_prepare_execution(
    MVMCPowerLanczosProductionSession *session,
    const MVMCPowerLanczosProductionExecution *execution) {
  MVMCPowerLanczosSnapshotSummary snapshot_summary;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  size_t term_bytes = 0;
  size_t operator_bytes = 0;
  size_t block_storage_count = 0;
  size_t final_storage_count = 0;
  size_t coefficient_evidence_sum_count = 0;
  size_t block_sum_count = 0;
  size_t leave_one_count = 0;
  size_t digest_count = 0;
  size_t provenance_digest_count = 0;
  size_t observable_coefficient_sample_count = 0;
  size_t observable_coefficient_block_count = 0;
  size_t observable_final_block_count = 0;
  size_t observable_result_block_count = 0;
  size_t observable_leave_one_coefficient_count = 0;
  size_t coefficient_primitive_count = 0;
  size_t final_primitive_count = 0;
  size_t maximum_primitive_count = 0;
  size_t maximum_entries_per_block = 0;
  size_t owned_bytes = 0;
  size_t resource_bytes = 0;
  uint64_t fingerprint = 0;
#ifdef _mpi_use
  int initialized = 0;
  int world_rank = -1;
  int world_size = 0;
  int chain_rank = -1;
  int chain_size = 0;
  int leader_rank = -1;
  int leader_size = 0;
  uint64_t minimum_fingerprint = 0;
  uint64_t maximum_fingerprint = 0;
#endif
  if (session == NULL || execution == NULL || !session->valid ||
      session->status != MVMC_KRYLOV_STATUS_OK ||
      session->state !=
          MVMC_POWER_LANCZOS_PRODUCTION_SESSION_SNAPSHOTS_VERIFIED ||
      session->execution_prepared) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#ifdef _mpi_use
  if (execution->world_communicator == MPI_COMM_NULL ||
      execution->chain_communicator == MPI_COMM_NULL ||
      MPI_Initialized(&initialized) != MPI_SUCCESS || !initialized ||
      MPI_Comm_dup(execution->world_communicator,
                   &session->world_communicator) != MPI_SUCCESS) {
    return mvmc_power_lanczos_production_session_fail_internal(
        session, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
  }
  if (MPI_Comm_set_errhandler(session->world_communicator,
                              MPI_ERRORS_RETURN) != MPI_SUCCESS ||
      MPI_Comm_rank(session->world_communicator, &world_rank) != MPI_SUCCESS ||
      MPI_Comm_size(session->world_communicator, &world_size) != MPI_SUCCESS ||
      world_rank < 0 || world_size <= 0 ||
      (size_t)world_rank != session->config.mpi_world_rank ||
      (size_t)world_size != session->config.mpi_world_size) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#else
  (void)execution->world_communicator;
  (void)execution->chain_communicator;
  if (session->config.mpi_world_rank != 0 ||
      session->config.mpi_world_size != 1 || session->config.split_size != 1) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#endif
  status = synchronize_world_status(session, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
#ifdef _mpi_use
  if (MPI_Comm_dup(execution->chain_communicator,
                   &session->chain_communicator) != MPI_SUCCESS ||
      MPI_Comm_set_errhandler(session->chain_communicator,
                              MPI_ERRORS_RETURN) != MPI_SUCCESS ||
      MPI_Comm_rank(session->chain_communicator, &chain_rank) != MPI_SUCCESS ||
      MPI_Comm_size(session->chain_communicator, &chain_size) != MPI_SUCCESS ||
      chain_rank < 0 || chain_size <= 0) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    const size_t group_base =
        session->config.mpi_world_rank / session->config.split_size *
        session->config.split_size;
    size_t expected_size = session->config.mpi_world_size - group_base;
    if (expected_size > session->config.split_size) {
      expected_size = session->config.split_size;
    }
    if ((size_t)chain_rank != session->config.mpi_world_rank - group_base ||
        (size_t)chain_size != expected_size) {
      status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  status = synchronize_world_status(session, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  session->chain_rank = (size_t)chain_rank;
  session->chain_size = (size_t)chain_size;
#else
  session->chain_rank = 0;
  session->chain_size = 1;
#endif
  session->chain_count = session->config.mpi_world_size /
                         session->config.split_size;
  if (session->config.mpi_world_size % session->config.split_size != 0) {
    ++session->chain_count;
  }
  status = valid_execution_shape(session, execution,
                                 &maximum_entries_per_block)
               ? MVMC_KRYLOV_STATUS_OK
               : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_snapshot_summary(session->snapshot,
                                                 &snapshot_summary);
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (snapshot_summary.word_count !=
           mvmc_krylov_fock_word_count(execution->proposal_model->site_count) ||
       snapshot_summary.configuration_bit_count !=
           session->config.configuration_bit_count ||
       session->snapshot_selection.expected_particle_count !=
           execution->proposal_model->up_electron_count +
               execution->proposal_model->down_electron_count)) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  fingerprint = status == MVMC_KRYLOV_STATUS_OK
                    ? execution_fingerprint(session, execution)
                    : 0;
#ifdef _mpi_use
  if (MPI_Allreduce(&fingerprint, &minimum_fingerprint, 1, MPI_UINT64_T,
                    MPI_MIN, session->world_communicator) != MPI_SUCCESS ||
      MPI_Allreduce(&fingerprint, &maximum_fingerprint, 1, MPI_UINT64_T,
                    MPI_MAX, session->world_communicator) != MPI_SUCCESS) {
    status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  } else if (status == MVMC_KRYLOV_STATUS_OK &&
             (minimum_fingerprint == 0 ||
              minimum_fingerprint != maximum_fingerprint)) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
#endif
  status = synchronize_world_status(session, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
#ifdef _mpi_use
  if (MPI_Comm_split(session->world_communicator,
                     session->chain_rank == 0 ? 0 : MPI_UNDEFINED,
                     world_rank, &session->leader_communicator) != MPI_SUCCESS) {
    status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  } else if (session->chain_rank == 0 &&
             (session->leader_communicator == MPI_COMM_NULL ||
              MPI_Comm_set_errhandler(session->leader_communicator,
                                      MPI_ERRORS_RETURN) != MPI_SUCCESS ||
              MPI_Comm_rank(session->leader_communicator,
                            &leader_rank) != MPI_SUCCESS ||
              MPI_Comm_size(session->leader_communicator,
                            &leader_size) != MPI_SUCCESS ||
              leader_rank < 0 || leader_size <= 0 ||
              (size_t)leader_rank !=
                  session->config.mpi_world_rank /
                      session->config.split_size ||
              (size_t)leader_size != session->chain_count)) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  } else if (session->chain_rank != 0 &&
             session->leader_communicator != MPI_COMM_NULL) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  session->leader_rank = leader_rank < 0 ? SIZE_MAX : (size_t)leader_rank;
#else
  session->leader_rank = 0;
#endif
  status = synchronize_world_status(session, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  if (!checked_multiply(execution->proposal_model->term_count,
                        sizeof(*session->model_terms), &term_bytes) ||
      !checked_multiply(execution->proposal_model->operator_count,
                        sizeof(*session->model_operators), &operator_bytes)) {
    status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  if (status == MVMC_KRYLOV_STATUS_OK && term_bytes != 0) {
    session->model_terms =
        (MVMCKrylovHamiltonianTerm *)allocate_owned(
            execution->proposal_model->term_count,
            sizeof(*session->model_terms), &owned_bytes, &status);
  }
  if (status == MVMC_KRYLOV_STATUS_OK && operator_bytes != 0) {
    session->model_operators =
        (MVMCKrylovFermionOperator *)allocate_owned(
            execution->proposal_model->operator_count,
            sizeof(*session->model_operators), &owned_bytes, &status);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    if (term_bytes != 0) {
      memcpy(session->model_terms, execution->proposal_model->terms,
             term_bytes);
    }
    if (operator_bytes != 0) {
      memcpy(session->model_operators, execution->proposal_model->operators,
             operator_bytes);
    }
    session->model = *execution->proposal_model;
    session->model.terms = session->model_terms;
    session->model.operators = session->model_operators;
  }
  if (status == MVMC_KRYLOV_STATUS_OK && execution->observable_enabled) {
    MVMCPowerLanczosObservablePlan copied = *execution->observable_plan;
    MVMCPowerLanczosObservableRecord *records =
        (MVMCPowerLanczosObservableRecord *)allocate_owned(
            (size_t)execution->observable_plan->record_count,
            sizeof(*records), &owned_bytes, &status);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      memcpy(records, execution->observable_plan->records,
             (size_t)execution->observable_plan->record_count *
                 sizeof(*records));
      copied.records = records;
      session->observable_layout = execution->observable_layout;
      session->observable_plan = copied;
      session->observable_request_count =
          (size_t)execution->observable_plan->record_count;
    }
  }
  status = synchronize_world_status(session, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  if (!checked_multiply(session->upper_count, 4,
                        &session->coefficient_entries_per_block)) {
    status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  session->final_entries_per_block = 2;
  if (status == MVMC_KRYLOV_STATUS_OK && execution->observable_enabled) {
    size_t observable_matrix_entries = 0;
    if (!checked_multiply(session->observable_request_count,
                          MVMC_POWER_LANCZOS_OBSERVABLE_MATRIX_ENTRIES,
                          &observable_matrix_entries) ||
        !checked_add(session->coefficient_entries_per_block,
                     observable_matrix_entries,
                     &session->coefficient_entries_per_block) ||
        !checked_add(session->final_entries_per_block,
                     session->observable_request_count,
                     &session->final_entries_per_block)) {
      status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (!checked_add(session->coefficient_entries_per_block, 1,
                    &coefficient_primitive_count) ||
       !checked_add(session->final_entries_per_block, 0,
                    &final_primitive_count))) {
    status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  maximum_primitive_count =
      coefficient_primitive_count > final_primitive_count
          ? coefficient_primitive_count
          : final_primitive_count;
  status = synchronize_world_status(session, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  status = mvmc_bounded_krylov_plan_create(
      &session->model, &execution->bounded_limits, &session->bounded_plan);
  status = synchronize_world_status(session, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  resource_bytes = mvmc_bounded_krylov_plan_bytes(session->bounded_plan);
  if (!checked_add(owned_bytes, resource_bytes, &owned_bytes)) {
    status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  status = synchronize_world_status(session, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  status = mvmc_bounded_krylov_workspace_create(
      session->bounded_plan, &session->bounded_workspace);
  status = synchronize_world_status(session, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  resource_bytes =
      mvmc_bounded_krylov_workspace_bytes(session->bounded_workspace);
  if (!checked_add(owned_bytes, resource_bytes, &owned_bytes)) {
    status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  status = synchronize_world_status(session, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  if (execution->observable_enabled) {
    MVMCPowerLanczosObservableBlockSummary block_summary;
    status = mvmc_power_lanczos_observable_evaluator_workspace_create(
        session->observable_request_count, snapshot_summary.word_count,
        &session->observable_evaluator);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      resource_bytes =
          mvmc_power_lanczos_observable_evaluator_workspace_bytes(
              session->observable_evaluator);
      if (!checked_add(owned_bytes, resource_bytes, &owned_bytes)) {
        status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      }
    }
    status = synchronize_world_status(session, status);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_power_lanczos_observable_block_create(
          MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT,
          session->observable_request_count, execution->block_count,
          (uint64_t)(execution->coefficient_sample_count /
                     execution->block_count),
          &session->observable_coefficient_blocks);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_power_lanczos_observable_block_summary(
          session->observable_coefficient_blocks, &block_summary);
    }
    if (status == MVMC_KRYLOV_STATUS_OK &&
        !checked_add(owned_bytes, block_summary.allocated_bytes,
                     &owned_bytes)) {
      status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    status = synchronize_world_status(session, status);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_power_lanczos_observable_block_create(
          MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_FINAL,
          session->observable_request_count, execution->block_count,
          (uint64_t)(execution->final_sample_count /
                     execution->block_count),
          &session->observable_final_blocks);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_power_lanczos_observable_block_summary(
          session->observable_final_blocks, &block_summary);
    }
    if (status == MVMC_KRYLOV_STATUS_OK &&
        !checked_add(owned_bytes, block_summary.allocated_bytes,
                     &owned_bytes)) {
      status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    }
    status = synchronize_world_status(session, status);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = confirmation_status(
          mvmc_power_lanczos_observable_confirmation_create(
              session->observable_request_count, execution->block_count,
              &execution->gevp_policy, &session->observable_confirmation));
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      resource_bytes =
          mvmc_power_lanczos_observable_confirmation_allocated_bytes(
              session->observable_confirmation);
      if (!checked_add(owned_bytes, resource_bytes, &owned_bytes)) {
        status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      }
    }
    status = synchronize_world_status(session, status);
    if (status != MVMC_KRYLOV_STATUS_OK) {
      return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                  status);
    }
  }
  status = mvmc_power_lanczos_chain_rng_create(
#ifdef _mpi_use
      session->chain_communicator, session->config.mpi_world_rank,
#else
      execution->chain_communicator, session->config.mpi_world_rank,
#endif
      session->config.mpi_world_size, session->config.split_size,
      &session->chain_rng_workspace);
  status = synchronize_world_status(session, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  resource_bytes = mvmc_power_lanczos_chain_rng_allocated_bytes(
      session->chain_rng_workspace);
  if (!checked_add(owned_bytes, resource_bytes, &owned_bytes)) {
    status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  status = synchronize_world_status(session, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
#ifdef _mpi_use
  if (status == MVMC_KRYLOV_STATUS_OK && session->chain_rank == 0) {
    status = mvmc_power_lanczos_block_collective_create(
        session->leader_communicator, execution->block_count,
        maximum_entries_per_block, &session->block_collective_workspace);
  }
#else
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_block_collective_create(
        0, execution->block_count, maximum_entries_per_block,
        &session->block_collective_workspace);
  }
#endif
  status = synchronize_world_status(session, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  resource_bytes = mvmc_power_lanczos_block_collective_allocated_bytes(
      session->block_collective_workspace);
  if (!checked_add(owned_bytes, resource_bytes, &owned_bytes) ||
      !checked_multiply(execution->block_count, session->upper_count,
                        &block_storage_count) ||
      !checked_multiply(block_storage_count, 4,
                        &coefficient_evidence_sum_count) ||
      !checked_multiply(execution->block_count, 2, &final_storage_count) ||
      !checked_multiply(execution->block_count,
                        maximum_entries_per_block, &block_sum_count) ||
      !checked_multiply(execution->block_count, session->basis_count,
                        &leave_one_count) ||
      !checked_multiply(session->chain_size,
                        MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES,
                        &digest_count) ||
      !checked_multiply(session->chain_count,
                        MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES,
                        &provenance_digest_count)) {
    status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  if (status == MVMC_KRYLOV_STATUS_OK && execution->observable_enabled &&
      (!checked_multiply(
           session->observable_request_count,
           MVMC_POWER_LANCZOS_OBSERVABLE_MATRIX_ENTRIES,
           &observable_coefficient_sample_count) ||
       !checked_multiply(execution->block_count,
                         observable_coefficient_sample_count,
                         &observable_coefficient_block_count) ||
       !checked_multiply(execution->block_count,
                         session->observable_request_count,
                         &observable_final_block_count) ||
       !checked_multiply(execution->block_count,
                         session->observable_request_count,
                         &observable_result_block_count) ||
       !checked_multiply(
           execution->block_count,
           MVMC_POWER_LANCZOS_CONFIRMATION_DIMENSION,
           &observable_leave_one_coefficient_count))) {
    status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  session->coefficient_blocks =
      (MVMCKrylovMatrixMeasurementAccumulator *)allocate_owned(
          execution->block_count, sizeof(*session->coefficient_blocks),
          &owned_bytes, &status);
  session->coefficient_overlap =
      (MVMCKrylovStreamingComplexSum *)allocate_owned(
          block_storage_count, sizeof(*session->coefficient_overlap),
          &owned_bytes, &status);
  session->coefficient_hamiltonian =
      (MVMCKrylovStreamingComplexSum *)allocate_owned(
          block_storage_count, sizeof(*session->coefficient_hamiltonian),
          &owned_bytes, &status);
  session->coefficient_hamiltonian_adjoint =
      (MVMCKrylovStreamingComplexSum *)allocate_owned(
          block_storage_count,
          sizeof(*session->coefficient_hamiltonian_adjoint), &owned_bytes,
          &status);
  session->coefficient_hamiltonian_squared =
      (MVMCKrylovStreamingComplexSum *)allocate_owned(
          block_storage_count,
          sizeof(*session->coefficient_hamiltonian_squared), &owned_bytes,
          &status);
  session->final_local_accumulators =
      (MVMCKrylovStreamingComplexSum *)allocate_owned(
          final_storage_count, sizeof(*session->final_local_accumulators),
          &owned_bytes, &status);
  session->observable_coefficient_sample = (double complex *)allocate_owned(
      observable_coefficient_sample_count,
      sizeof(*session->observable_coefficient_sample), &owned_bytes, &status);
  session->observable_final_sample = (double complex *)allocate_owned(
      session->observable_request_count,
      sizeof(*session->observable_final_sample), &owned_bytes, &status);
  session->observable_numeric_evidence_sample =
      (MVMCPowerLanczosNumericEvidence *)allocate_owned(
          observable_coefficient_sample_count,
          sizeof(*session->observable_numeric_evidence_sample),
          &owned_bytes, &status);
  session->primitive_trace_values = (double complex *)allocate_owned(
      session->chain_rank == 0 ? maximum_primitive_count : 0,
      sizeof(*session->primitive_trace_values), &owned_bytes, &status);
  session->primitive_trace_bounds = (double *)allocate_owned(
      session->chain_rank == 0 ? maximum_primitive_count : 0,
      sizeof(*session->primitive_trace_bounds), &owned_bytes, &status);
  session->primitive_trace_support = (uint8_t *)allocate_owned(
      session->chain_rank == 0 ? maximum_primitive_count : 0,
      sizeof(*session->primitive_trace_support), &owned_bytes, &status);
  session->observable_coefficient_block_sums =
      (double complex *)allocate_owned(
          observable_coefficient_block_count,
          sizeof(*session->observable_coefficient_block_sums), &owned_bytes,
          &status);
  session->observable_final_block_sums = (double complex *)allocate_owned(
      observable_final_block_count,
      sizeof(*session->observable_final_block_sums), &owned_bytes, &status);
  session->observable_block_counts = (uint64_t *)allocate_owned(
      execution->observable_enabled ? execution->block_count : 0,
      sizeof(*session->observable_block_counts), &owned_bytes, &status);
  session->observable_coefficient_estimates =
      (double complex *)allocate_owned(
          session->observable_request_count,
          sizeof(*session->observable_coefficient_estimates), &owned_bytes,
          &status);
  session->observable_final_estimates = (double complex *)allocate_owned(
      session->observable_request_count,
      sizeof(*session->observable_final_estimates), &owned_bytes, &status);
  session->observable_leave_one_estimates =
      (double complex *)allocate_owned(
          observable_result_block_count,
          sizeof(*session->observable_leave_one_estimates), &owned_bytes,
          &status);
  session->observable_final_block_means =
      (double complex *)allocate_owned(
          observable_result_block_count,
          sizeof(*session->observable_final_block_means), &owned_bytes,
          &status);
  session->observable_leave_one_coefficients =
      (double complex *)allocate_owned(
          observable_leave_one_coefficient_count,
          sizeof(*session->observable_leave_one_coefficients), &owned_bytes,
          &status);
  session->observable_leave_one_projective_distances =
      (double *)allocate_owned(
          execution->observable_enabled ? execution->block_count : 0,
          sizeof(*session->observable_leave_one_projective_distances),
          &owned_bytes, &status);
  session->proposal_words = (uint64_t *)allocate_owned(
      snapshot_summary.word_count, sizeof(*session->proposal_words),
      &owned_bytes, &status);
  session->local_block_counts = (uint64_t *)allocate_owned(
      execution->block_count, sizeof(*session->local_block_counts),
      &owned_bytes, &status);
  session->global_block_counts = (uint64_t *)allocate_owned(
      execution->block_count, sizeof(*session->global_block_counts),
      &owned_bytes, &status);
  session->local_block_sums = (double complex *)allocate_owned(
      block_sum_count, sizeof(*session->local_block_sums), &owned_bytes,
      &status);
  session->global_block_sums = (double complex *)allocate_owned(
      block_sum_count, sizeof(*session->global_block_sums), &owned_bytes,
      &status);
  session->coefficient_evidence_block_counts =
      (uint64_t *)allocate_owned(
          execution->block_count,
          sizeof(*session->coefficient_evidence_block_counts), &owned_bytes,
          &status);
  session->coefficient_evidence_block_sums =
      (double complex *)allocate_owned(
          coefficient_evidence_sum_count,
          sizeof(*session->coefficient_evidence_block_sums), &owned_bytes,
          &status);
  session->final_evidence_block_counts = (uint64_t *)allocate_owned(
      execution->block_count, sizeof(*session->final_evidence_block_counts),
      &owned_bytes, &status);
  session->final_evidence_block_sums = (double complex *)allocate_owned(
      final_storage_count, sizeof(*session->final_evidence_block_sums),
      &owned_bytes, &status);
  session->leave_one_coefficients = (double complex *)allocate_owned(
      leave_one_count, sizeof(*session->leave_one_coefficients), &owned_bytes,
      &status);
  session->leave_one_energies = (double *)allocate_owned(
      execution->block_count, sizeof(*session->leave_one_energies),
      &owned_bytes, &status);
  session->leave_one_second_moments = (double complex *)allocate_owned(
      execution->block_count, sizeof(*session->leave_one_second_moments),
      &owned_bytes, &status);
  session->leave_one_projective_distances = (double *)allocate_owned(
      execution->block_count,
      sizeof(*session->leave_one_projective_distances), &owned_bytes,
      &status);
  session->chain_digest_buffer = (unsigned char *)allocate_owned(
      digest_count, sizeof(*session->chain_digest_buffer), &owned_bytes,
      &status);
  session->chain_provenance_digests = (unsigned char *)allocate_owned(
      provenance_digest_count,
      sizeof(*session->chain_provenance_digests), &owned_bytes, &status);
  if (status == MVMC_KRYLOV_STATUS_OK && session->chain_rank == 0) {
    status = mvmc_power_lanczos_primitive_trace_create(
        MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_COEFFICIENT,
        coefficient_primitive_count, 1,
        execution->coefficient_sample_count,
        execution->bounded_limits.max_workspace_bytes,
        &session->coefficient_local_trace);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_power_lanczos_primitive_trace_register_group(
          session->coefficient_local_trace, 0,
          session->config.mpi_world_rank, session->chain_size);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      resource_bytes =
          mvmc_power_lanczos_primitive_trace_allocated_bytes(
              session->coefficient_local_trace);
      if (!checked_add(owned_bytes, resource_bytes, &owned_bytes)) {
        status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      }
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK && session->chain_rank == 0) {
    status = mvmc_power_lanczos_primitive_trace_create(
        MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_FINAL,
        final_primitive_count, 1, execution->final_sample_count,
        execution->bounded_limits.max_workspace_bytes,
        &session->final_local_trace);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_power_lanczos_primitive_trace_register_group(
          session->final_local_trace, 0,
          session->config.mpi_world_rank, session->chain_size);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      resource_bytes =
          mvmc_power_lanczos_primitive_trace_allocated_bytes(
              session->final_local_trace);
      if (!checked_add(owned_bytes, resource_bytes, &owned_bytes)) {
        status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
      }
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_krylov_matrix_measurement_block_accumulator_init(
        (int)session->basis_count, execution->block_count,
        (uint64_t)(execution->coefficient_sample_count /
                   execution->block_count),
        session->coefficient_blocks, session->coefficient_overlap,
        session->coefficient_hamiltonian,
        session->coefficient_hamiltonian_adjoint,
        session->coefficient_hamiltonian_squared, block_storage_count,
        &session->coefficient_accumulator);
  }
  if (status == MVMC_KRYLOV_STATUS_OK &&
      !checked_add(session->allocated_bytes, owned_bytes,
                   &resource_bytes)) {
    status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  status = synchronize_world_status(session, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    return mvmc_power_lanczos_production_session_fail_internal(session,
                                                                status);
  }
  session->execution = *execution;
  session->execution_fingerprint = fingerprint;
  session->execution.proposal_model = &session->model;
  session->execution.observable_plan = execution->observable_enabled
                                           ? &session->observable_plan
                                           : NULL;
#ifdef _mpi_use
  session->execution.world_communicator = session->world_communicator;
  session->execution.chain_communicator = session->chain_communicator;
#endif
  session->block_entry_count = maximum_entries_per_block;
  session->coefficient_primitive_count = coefficient_primitive_count;
  session->final_primitive_count = final_primitive_count;
  session->maximum_primitive_count = maximum_primitive_count;
  session->allocated_bytes = resource_bytes;
  session->execution_prepared = 1;
  return MVMC_KRYLOV_STATUS_OK;
}
