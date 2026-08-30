#ifndef MVMC_POWER_LANCZOS_OBSERVABLE_CONFIRMATION_H
#define MVMC_POWER_LANCZOS_OBSERVABLE_CONFIRMATION_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "power_lanczos_gevp.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_CONFIRMATION_DIMENSION 2
#define MVMC_POWER_LANCZOS_CONFIRMATION_UPPER_COUNT 3
#define MVMC_POWER_LANCZOS_CONFIRMATION_MATRIX_FAMILY_COUNT 4
#define MVMC_POWER_LANCZOS_CONFIRMATION_MIN_BLOCK_COUNT 32
#define MVMC_POWER_LANCZOS_CONFIRMATION_MAX_EXACT_SAMPLE_COUNT \
  UINT64_C(9007199254740992)

typedef enum {
  MVMC_POWER_LANCZOS_CONFIRMATION_OK = 0,
  MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_ARGUMENT,
  MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_STATE,
  MVMC_POWER_LANCZOS_CONFIRMATION_BLOCK_CONTRACT,
  MVMC_POWER_LANCZOS_CONFIRMATION_RESOURCE_LIMIT,
  MVMC_POWER_LANCZOS_CONFIRMATION_ALLOCATION_FAILURE,
  MVMC_POWER_LANCZOS_CONFIRMATION_NONFINITE,
  MVMC_POWER_LANCZOS_CONFIRMATION_GEVP_FAILURE,
  MVMC_POWER_LANCZOS_CONFIRMATION_INTERNAL_FAILURE
} MVMCPowerLanczosObservableConfirmationStatus;

typedef enum {
  MVMC_POWER_LANCZOS_CONFIRMATION_CONFIGURED = 0,
  MVMC_POWER_LANCZOS_CONFIRMATION_COEFFICIENT_SAMPLING,
  MVMC_POWER_LANCZOS_CONFIRMATION_COEFFICIENT_READY,
  MVMC_POWER_LANCZOS_CONFIRMATION_COEFFICIENT_FROZEN,
  MVMC_POWER_LANCZOS_CONFIRMATION_FINAL_READY,
  MVMC_POWER_LANCZOS_CONFIRMATION_FINALIZED,
  MVMC_POWER_LANCZOS_CONFIRMATION_FAILED
} MVMCPowerLanczosObservableConfirmationState;

typedef struct MVMCPowerLanczosObservableConfirmationSession
    MVMCPowerLanczosObservableConfirmationSession;

typedef struct {
  int valid;
  MVMCPowerLanczosObservableConfirmationStatus status;
  MVMCPowerLanczosObservableConfirmationState state;
  size_t request_count;
  size_t block_count;
  size_t coefficient_blocks_added;
  size_t final_blocks_added;
  uint64_t coefficient_block_length;
  uint64_t final_block_length;
  uint64_t coefficient_sample_count;
  uint64_t final_sample_count;
  double complex coefficient[MVMC_POWER_LANCZOS_CONFIRMATION_DIMENSION];
  double energy;
  double maximum_leave_one_projective_distance;
} MVMCPowerLanczosObservableConfirmationSummary;

/*
 * The session consumes rank-reduced nonoverlap block sums.  Both stages must
 * provide exactly block_count blocks, and every block in both stages must use
 * the same nonzero sample_count no greater than 2^53.  S, K_forward,
 * K_reverse and B use upper-packed row-then-column order (00,01,11).
 * observable_matrix_sums use request-major full row order (00,01,10,11).
 */
MVMCPowerLanczosObservableConfirmationStatus
mvmc_power_lanczos_observable_confirmation_create(
    size_t request_count, size_t block_count,
    const MVMCPowerLanczosGEVPPolicy *gevp_policy,
    MVMCPowerLanczosObservableConfirmationSession **session);

void mvmc_power_lanczos_observable_confirmation_destroy(
    MVMCPowerLanczosObservableConfirmationSession *session);

MVMCPowerLanczosObservableConfirmationStatus
mvmc_power_lanczos_observable_confirmation_add_coefficient_block(
    MVMCPowerLanczosObservableConfirmationSession *session,
    uint64_t sample_count,
    const double complex overlap_upper_sums[3],
    const double complex hamiltonian_forward_upper_sums[3],
    const double complex hamiltonian_reverse_upper_sums[3],
    const double complex hamiltonian_squared_upper_sums[3],
    const double complex *observable_matrix_sums,
    size_t observable_matrix_count);

MVMCPowerLanczosObservableConfirmationStatus
mvmc_power_lanczos_observable_confirmation_freeze_coefficient(
    MVMCPowerLanczosObservableConfirmationSession *session);

MVMCPowerLanczosObservableConfirmationStatus
mvmc_power_lanczos_observable_confirmation_add_final_block(
    MVMCPowerLanczosObservableConfirmationSession *session,
    uint64_t sample_count, const double complex *observable_sums,
    size_t observable_count);

/*
 * Finalize and atomically export all requested arrays.  leave_one_estimates
 * and final_block_means are block-major/request-major.  leave_one_coefficient
 * is block-major with two complex values per block.  All six output ranges
 * must be pairwise disjoint; failures leave every output untouched.
 */
MVMCPowerLanczosObservableConfirmationStatus
mvmc_power_lanczos_observable_confirmation_finalize(
    MVMCPowerLanczosObservableConfirmationSession *session,
    double complex *coefficient_estimates, size_t coefficient_capacity,
    double complex *final_estimates, size_t final_capacity,
    double complex *leave_one_estimates, size_t leave_one_capacity,
    double complex *final_block_means, size_t final_block_capacity,
    double complex *leave_one_coefficient,
    size_t leave_one_coefficient_capacity,
    double *leave_one_projective_distance,
    size_t leave_one_projective_capacity);

MVMCPowerLanczosObservableConfirmationStatus
mvmc_power_lanczos_observable_confirmation_summary(
    const MVMCPowerLanczosObservableConfirmationSession *session,
    MVMCPowerLanczosObservableConfirmationSummary *summary);

size_t mvmc_power_lanczos_observable_confirmation_allocated_bytes(
    const MVMCPowerLanczosObservableConfirmationSession *session);

const char *mvmc_power_lanczos_observable_confirmation_status_string(
    MVMCPowerLanczosObservableConfirmationStatus status);

#endif

#endif
