#ifndef MVMC_POWER_LANCZOS_PRODUCTION_SESSION_INTERNAL_H
#define MVMC_POWER_LANCZOS_PRODUCTION_SESSION_INTERNAL_H

#include "krylov_final_state.h"
#include "power_lanczos_observable_blocks.h"
#include "power_lanczos_observable_evaluator.h"
#include "power_lanczos_primitive_trace.h"
#include "power_lanczos_production_session.h"

struct MVMCPowerLanczosProductionSession {
  int valid;
  int execution_prepared;
  MVMCKrylovStatus status;
  MVMCPowerLanczosProductionSessionState state;
  MVMCPowerLanczosProductionSessionConfig config;
  MVMCPowerLanczosProductionExecution execution;
  size_t basis_count;
  size_t upper_count;
  size_t allocated_bytes;
  size_t chain_rank;
  size_t chain_size;
  size_t chain_count;
  size_t leader_rank;
  size_t block_entry_count;
  size_t coefficient_entries_per_block;
  size_t final_entries_per_block;
  size_t coefficient_primitive_count;
  size_t final_primitive_count;
  size_t maximum_primitive_count;
  size_t observable_request_count;
  uint64_t execution_fingerprint;
  MVMCPowerLanczosTerminalSnapshotSelection snapshot_selection;
  MVMCPowerLanczosSnapshot *snapshot;
  MVMCPowerLanczosRngDomain base_rng;
  MVMCPowerLanczosRngDomain coefficient_rng;
  MVMCPowerLanczosRngDomain final_rng;
  MVMCKrylovFockModel model;
  MVMCKrylovHamiltonianTerm *model_terms;
  MVMCKrylovFermionOperator *model_operators;
  MVMCKrylovBoundedPlan *bounded_plan;
  MVMCKrylovBoundedWorkspace *bounded_workspace;
  MVMCPowerLanczosChainRngWorkspace *chain_rng_workspace;
  MVMCPowerLanczosBlockCollectiveWorkspace *block_collective_workspace;
  MVMCPowerLanczosObservableEvaluatorWorkspace *observable_evaluator;
  MVMCPowerLanczosObservableBlockAccumulator *observable_coefficient_blocks;
  MVMCPowerLanczosObservableBlockAccumulator *observable_final_blocks;
  MVMCPowerLanczosObservableConfirmationSession *observable_confirmation;
  MVMCPowerLanczosPrimitiveTrace *coefficient_local_trace;
  MVMCPowerLanczosPrimitiveTrace *final_local_trace;
  MVMCPowerLanczosPrimitiveTrace *coefficient_trace;
  MVMCPowerLanczosPrimitiveTrace *final_trace;
  MVMCPowerLanczosObservableLayout observable_layout;
  MVMCPowerLanczosObservablePlan observable_plan;
  MVMCKrylovMatrixMeasurementBlockAccumulator coefficient_accumulator;
  MVMCKrylovMatrixMeasurementAccumulator *coefficient_blocks;
  MVMCKrylovStreamingComplexSum *coefficient_overlap;
  MVMCKrylovStreamingComplexSum *coefficient_hamiltonian;
  MVMCKrylovStreamingComplexSum *coefficient_hamiltonian_adjoint;
  MVMCKrylovStreamingComplexSum *coefficient_hamiltonian_squared;
  MVMCKrylovStreamingComplexSum *final_local_accumulators;
  double complex *observable_coefficient_sample;
  double complex *observable_final_sample;
  MVMCPowerLanczosNumericEvidence *observable_numeric_evidence_sample;
  double complex *primitive_trace_values;
  double *primitive_trace_bounds;
  uint8_t *primitive_trace_support;
  double complex *observable_coefficient_block_sums;
  double complex *observable_final_block_sums;
  uint64_t *observable_block_counts;
  double complex *observable_coefficient_estimates;
  double complex *observable_final_estimates;
  double complex *observable_leave_one_estimates;
  double complex *observable_final_block_means;
  double complex *observable_leave_one_coefficients;
  double *observable_leave_one_projective_distances;
  uint64_t *proposal_words;
  uint64_t *local_block_counts;
  uint64_t *global_block_counts;
  double complex *local_block_sums;
  double complex *global_block_sums;
  uint64_t *coefficient_evidence_block_counts;
  double complex *coefficient_evidence_block_sums;
  uint64_t *final_evidence_block_counts;
  double complex *final_evidence_block_sums;
  int coefficient_evidence_valid;
  int final_evidence_valid;
  double complex *leave_one_coefficients;
  double *leave_one_energies;
  double complex *leave_one_second_moments;
  double *leave_one_projective_distances;
  unsigned char *chain_digest_buffer;
  unsigned char *chain_provenance_digests;
  uint64_t *coefficient_words;
  uint64_t *final_words;
  size_t word_count;
  MVMCKrylovPositiveSamplerRng coefficient_rank_rng;
  MVMCKrylovPositiveSamplerRng final_rank_rng;
  MVMCKrylovPositiveSamplerSnapshot coefficient_sampler;
  MVMCKrylovFinalStateSnapshot final_sampler;
  uint64_t coefficient_proposals;
  uint64_t coefficient_accepted_steps;
  uint64_t scale_pilot_proposals;
  uint64_t scale_pilot_accepted_steps;
  uint64_t scale_pilot_sample_count;
  uint64_t final_proposals;
  uint64_t final_accepted_steps;
  uint64_t coefficient_sample_count;
  uint64_t final_sample_count;
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
  double complex coefficient_overlap_mean[6];
  double complex coefficient_hamiltonian_mean[6];
  double complex coefficient_hamiltonian_adjoint_mean[6];
  double complex coefficient_hamiltonian_squared_mean[6];
  MVMCKrylovFinalStatePolicy final_policy;
#ifdef _mpi_use
  MPI_Comm world_communicator;
  MPI_Comm chain_communicator;
  MPI_Comm leader_communicator;
#endif
};

MVMCKrylovStatus mvmc_power_lanczos_production_session_fail_internal(
    MVMCPowerLanczosProductionSession *session,
    MVMCKrylovStatus status);

#endif
