#ifndef MVMC_POWER_LANCZOS_PRODUCTION_SESSION_H
#define MVMC_POWER_LANCZOS_PRODUCTION_SESSION_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "krylov_matrix_measurement.h"
#include "power_lanczos_block_collective.h"
#include "power_lanczos_chain_rng.h"
#include "power_lanczos_gevp.h"
#include "power_lanczos_observable_action.h"
#include "power_lanczos_observable_confirmation.h"
#include "power_lanczos_primitive_trace.h"
#include "power_lanczos_rng.h"
#include "power_lanczos_snapshot.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_PRODUCTION_SESSION_VERSION UINT64_C(6)
#define MVMC_POWER_LANCZOS_TERMINAL_SNAPSHOT_SELECTION_VERSION UINT64_C(1)

typedef enum {
  MVMC_POWER_LANCZOS_TERMINAL_SNAPSHOT_WORLD_RANK_0_BASE_TERMINAL_AFTER_LAST_PROPOSAL = 1
} MVMCPowerLanczosTerminalSnapshotSelectionKind;

typedef struct {
  const uint64_t *canonical_words;
  size_t word_count;
  size_t configuration_bit_count;
  uint64_t generation;
  uint64_t terminal_proposal_counter;
  size_t expected_particle_count;
  double log_sampling_support;
} MVMCPowerLanczosTerminalSnapshotInput;

typedef struct {
  int valid;
  uint64_t version;
  MVMCPowerLanczosTerminalSnapshotSelectionKind selection_kind;
  size_t owner_world_rank;
  uint64_t terminal_proposal_counter;
  uint64_t generation;
  size_t word_count;
  size_t configuration_bit_count;
  size_t expected_particle_count;
  unsigned char snapshot_sha256[
      MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES];
} MVMCPowerLanczosTerminalSnapshotSelection;

typedef enum {
  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_CREATED = 0,
  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_SNAPSHOTS_VERIFIED,
  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_RUNNING,
  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_READY,
  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_FROZEN,
  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINAL_RUNNING,
  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINAL_READY,
  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINALIZED,
  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FAILED
} MVMCPowerLanczosProductionSessionState;

typedef struct {
  int power_step;
  uint64_t resolved_base_seed;
  uint64_t run_index;
  size_t mpi_world_rank;
  size_t mpi_world_size;
  size_t split_size;
  uint64_t base_generation;
  size_t configuration_bit_count;
} MVMCPowerLanczosProductionSessionConfig;

typedef struct MVMCPowerLanczosProductionSession
    MVMCPowerLanczosProductionSession;

typedef struct {
  MVMCKrylovBoundedCommunicator world_communicator;
  MVMCKrylovBoundedCommunicator chain_communicator;
  const MVMCKrylovFockModel *proposal_model;
  MVMCKrylovScaledAmplitudeCallback amplitude;
  void *amplitude_context;
  uint64_t amplitude_generation_hash;
  MVMCKrylovBoundedLimits bounded_limits;
  MVMCKrylovPositiveGuidePolicy coefficient_guide_policy;
  MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
  MVMCKrylovMatrixMeasurementPolicy matrix_policy;
  MVMCPowerLanczosGEVPPolicy gevp_policy;
  size_t block_count;
  int scale_pilot_enabled;
  size_t scale_pilot_warm_up;
  size_t scale_pilot_sample_count;
  double eta_relative;
  size_t coefficient_warm_up;
  size_t coefficient_sample_count;
  size_t coefficient_interval;
  size_t final_warm_up;
  size_t final_sample_count;
  size_t final_interval;
  double maximum_leave_one_projective_distance;
  int observable_enabled;
  MVMCPowerLanczosObservableLayout observable_layout;
  const MVMCPowerLanczosObservablePlan *observable_plan;
} MVMCPowerLanczosProductionExecution;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  uint64_t version;
  MVMCPowerLanczosProductionSessionState state;
  int power_step;
  size_t basis_count;
  size_t upper_count;
  size_t allocated_bytes;
  int execution_prepared;
  size_t chain_rank;
  size_t chain_size;
  size_t chain_count;
  size_t block_count;
  uint64_t execution_fingerprint;
  int scale_pilot_enabled;
  uint64_t scale_pilot_warm_up;
  uint64_t scale_pilot_sample_count_per_chain;
  uint64_t scale_pilot_proposals;
  uint64_t scale_pilot_accepted_steps;
  uint64_t scale_pilot_sample_count;
  double eta_relative;
  double resolved_eta;
  double log_basis_scale[MVMC_KRYLOV_MAX_ORDER + 1];
  uint64_t coefficient_proposals;
  uint64_t coefficient_accepted_steps;
  uint64_t final_proposals;
  uint64_t final_accepted_steps;
  uint64_t coefficient_sample_count;
  uint64_t final_sample_count;
  int observable_enabled;
  size_t observable_request_count;
  char observable_census_sha256[
      MVMC_POWER_LANCZOS_OBSERVABLE_SHA256_HEX_CAPACITY];
  char observable_semantic_sha256[
      MVMC_POWER_LANCZOS_OBSERVABLE_SHA256_HEX_CAPACITY];
  double maximum_leave_one_projective_distance;
  double final_energy;
  double final_energy_imaginary;
  double corrected_variance;
  double corrected_variance_imaginary;
  double final_local_energy_abs_squared_diagnostic;
  uint64_t coefficient_provenance_hash;
  uint64_t final_policy_hash;
  unsigned char coefficient_provenance_sha256[
      MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES];
  MVMCPowerLanczosGEVPResult coefficient_gevp;
  MVMCPowerLanczosTerminalSnapshotSelection snapshot_selection;
  MVMCPowerLanczosSnapshotSummary snapshot;
  MVMCPowerLanczosRngDomain base_rng;
  MVMCPowerLanczosRngDomain coefficient_rng;
  MVMCPowerLanczosRngDomain final_rng;
} MVMCPowerLanczosProductionSessionSummary;

typedef struct {
  int valid;
  size_t request_index;
  int raw_ordinal;
  MVMCPowerLanczosObservableFamily family;
  char canonical_operator_id[
      MVMC_POWER_LANCZOS_OBSERVABLE_OPERATOR_ID_CAPACITY];
  double complex coefficient_estimate;
  double complex final_estimate;
} MVMCPowerLanczosProductionObservableResult;

typedef struct {
  int valid;
  size_t block_index;
  size_t request_index;
  double complex coefficient_leave_one_estimate;
  double complex final_block_mean;
  double complex leave_one_coefficient[
      MVMC_POWER_LANCZOS_CONFIRMATION_DIMENSION];
  double leave_one_projective_distance;
} MVMCPowerLanczosProductionObservableBlockResult;

/*
 * Immutable packed-upper coefficient-stage matrix evidence.  Only the first
 * upper_count entries in each array are populated; the fixed capacity covers
 * the supported P6 dimensions (step 1 and step 2).
 */
typedef struct {
  int valid;
  size_t basis_count;
  size_t upper_count;
  uint64_t sample_count;
  double complex overlap[6];
  double complex hamiltonian[6];
  double complex hamiltonian_adjoint[6];
  double complex hamiltonian_squared[6];
} MVMCPowerLanczosProductionMatrixResult;

typedef struct {
  int valid;
  size_t block_index;
  size_t basis_count;
  size_t upper_count;
  uint64_t sample_count;
  uint64_t leave_one_sample_count;
  double complex overlap_sum[6];
  double complex hamiltonian_sum[6];
  double complex hamiltonian_adjoint_sum[6];
  double complex hamiltonian_squared_sum[6];
  double complex leave_one_coefficient[3];
  double leave_one_energy;
  double complex leave_one_second_moment;
  double leave_one_projective_distance;
} MVMCPowerLanczosProductionCoefficientBlockResult;

typedef struct {
  int valid;
  size_t block_index;
  uint64_t sample_count;
  double complex energy_sum;
  double complex local_energy_abs_squared_sum;
} MVMCPowerLanczosProductionFinalBlockResult;

MVMCKrylovStatus mvmc_power_lanczos_production_dimension(
    int power_step, size_t *basis_count, size_t *upper_count);

MVMCKrylovStatus mvmc_power_lanczos_production_session_create(
    const MVMCPowerLanczosProductionSessionConfig *config,
    MVMCKrylovBoundedCommunicator world_communicator,
    const MVMCPowerLanczosTerminalSnapshotInput *rank_0_terminal,
    MVMCPowerLanczosProductionSession **session);

#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
MVMCKrylovStatus
mvmc_power_lanczos_production_session_testing_corrupt_terminal_after_broadcast(
    size_t target_world_rank, size_t word_index, uint64_t xor_mask);
#endif

void mvmc_power_lanczos_production_session_destroy(
    MVMCPowerLanczosProductionSession *session);

MVMCKrylovStatus mvmc_power_lanczos_production_session_verify(
    MVMCPowerLanczosProductionSession *session);

MVMCKrylovStatus mvmc_power_lanczos_production_session_prepare_execution(
    MVMCPowerLanczosProductionSession *session,
    const MVMCPowerLanczosProductionExecution *execution);

MVMCKrylovStatus mvmc_power_lanczos_production_session_run_coefficient(
    MVMCPowerLanczosProductionSession *session);

MVMCKrylovStatus mvmc_power_lanczos_production_session_freeze_coefficient(
    MVMCPowerLanczosProductionSession *session);

MVMCKrylovStatus mvmc_power_lanczos_production_session_run_final(
    MVMCPowerLanczosProductionSession *session);

MVMCKrylovStatus mvmc_power_lanczos_production_session_finalize(
    MVMCPowerLanczosProductionSession *session);

MVMCKrylovStatus mvmc_power_lanczos_production_session_summary(
    const MVMCPowerLanczosProductionSession *session,
    MVMCPowerLanczosProductionSessionSummary *summary);

MVMCKrylovStatus mvmc_power_lanczos_production_session_observable_result(
    const MVMCPowerLanczosProductionSession *session, size_t request_index,
    MVMCPowerLanczosProductionObservableResult *result);

MVMCKrylovStatus
mvmc_power_lanczos_production_session_observable_block_result(
    const MVMCPowerLanczosProductionSession *session, size_t block_index,
    size_t request_index,
    MVMCPowerLanczosProductionObservableBlockResult *result);

MVMCKrylovStatus mvmc_power_lanczos_production_session_matrix_result(
    const MVMCPowerLanczosProductionSession *session,
    MVMCPowerLanczosProductionMatrixResult *result);

MVMCKrylovStatus
mvmc_power_lanczos_production_session_coefficient_block_result(
    const MVMCPowerLanczosProductionSession *session, size_t block_index,
    MVMCPowerLanczosProductionCoefficientBlockResult *result);

MVMCKrylovStatus mvmc_power_lanczos_production_session_final_block_result(
    const MVMCPowerLanczosProductionSession *session, size_t block_index,
    MVMCPowerLanczosProductionFinalBlockResult *result);

MVMCKrylovStatus
mvmc_power_lanczos_production_session_coefficient_trace_summary(
    const MVMCPowerLanczosProductionSession *session,
    MVMCPowerLanczosPrimitiveTraceSummary *summary);

MVMCKrylovStatus mvmc_power_lanczos_production_session_final_trace_summary(
    const MVMCPowerLanczosProductionSession *session,
    MVMCPowerLanczosPrimitiveTraceSummary *summary);

MVMCKrylovStatus
mvmc_power_lanczos_production_session_coefficient_trace_group(
    const MVMCPowerLanczosProductionSession *session,
    size_t group_ordinal, MVMCPowerLanczosPrimitiveTraceGroup *group);

MVMCKrylovStatus mvmc_power_lanczos_production_session_final_trace_group(
    const MVMCPowerLanczosProductionSession *session,
    size_t group_ordinal, MVMCPowerLanczosPrimitiveTraceGroup *group);

MVMCKrylovStatus
mvmc_power_lanczos_production_session_coefficient_trace_sample(
    const MVMCPowerLanczosProductionSession *session,
    size_t group_ordinal, size_t sample_index, double complex *values,
    size_t value_capacity, double *absolute_numeric_bounds,
    size_t bound_capacity, uint8_t *support_flags,
    size_t support_capacity, double *tail_magnitude);

MVMCKrylovStatus mvmc_power_lanczos_production_session_final_trace_sample(
    const MVMCPowerLanczosProductionSession *session,
    size_t group_ordinal, size_t sample_index, double complex *values,
    size_t value_capacity, double *absolute_numeric_bounds,
    size_t bound_capacity, uint8_t *support_flags,
    size_t support_capacity, double *tail_magnitude);

MVMCKrylovStatus
mvmc_power_lanczos_production_session_coefficient_trace_export_scalars(
    const MVMCPowerLanczosProductionSession *session, double *scalars,
    size_t scalar_capacity);

MVMCKrylovStatus
mvmc_power_lanczos_production_session_final_trace_export_scalars(
    const MVMCPowerLanczosProductionSession *session, double *scalars,
    size_t scalar_capacity);

const char *mvmc_power_lanczos_production_session_state_string(
    MVMCPowerLanczosProductionSessionState state);

#endif

#endif
