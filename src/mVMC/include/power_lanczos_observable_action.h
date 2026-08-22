#ifndef MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_H
#define MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_H

#include "power_lanczos_observable_census.h"

#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_OBSERVABLE_ZERO_TARGET (-1)

typedef enum {
  MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK = 0,
  MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_ARGUMENT,
  MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_LAYOUT,
  MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_SOURCE,
  MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_RECORD,
  MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_TARGET,
  MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_RESOURCE_LIMIT,
  MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INTERNAL_INVARIANT_FAILURE
} MVMCPowerLanczosObservableActionStatus;

typedef struct {
  int nsite;
  int up_electron_count;
  int down_electron_count;
  int pure_spin;
} MVMCPowerLanczosObservableLayout;

typedef struct {
  int valid;
  MVMCPowerLanczosObservableActionStatus status;
  int applied;
  int fermion_sign;
} MVMCPowerLanczosObservableActionResult;

typedef struct {
  int valid;
  MVMCPowerLanczosObservableActionStatus status;
  size_t request_count;
  size_t active_request_count;
  size_t unique_target_count;
  int source_target_present;
} MVMCPowerLanczosObservableTargetGroupResult;

/*
 * Caller-owned transactional scratch for target grouping.  Production keeps
 * this storage in the persistent evaluator workspace so the per-sample path
 * does not consume a large stack frame or allocate dynamically.
 */
typedef struct {
  size_t max_requests;
  size_t word_count;
  uint64_t *target_words;
  int *request_target_index;
  int *request_fermion_sign;
  int *hash_slots;
  size_t hash_capacity;
} MVMCPowerLanczosObservableTargetGroupScratch;

/* Returns zero for an invalid layout. */
size_t mvmc_power_lanczos_observable_word_count(
    const MVMCPowerLanczosObservableLayout *layout);

MVMCPowerLanczosObservableActionStatus
mvmc_power_lanczos_observable_validate_configuration(
    const MVMCPowerLanczosObservableLayout *layout,
    const uint64_t *configuration_words, size_t word_count);

/*
 * Apply one raw monomial right-to-left in canonical spin-orbital order
 * p = site + spin * nsite.  On an error, or when applied is zero, target_words
 * is unchanged.  A successful Pauli-forbidden action, or an action projected
 * outside the physical pure-spin sector, has status OK, applied=0, and
 * fermion_sign=0.
 */
MVMCPowerLanczosObservableActionStatus
mvmc_power_lanczos_observable_apply(
    const MVMCPowerLanczosObservableLayout *layout,
    const MVMCPowerLanczosObservableRecord *record,
    const uint64_t *source_words, size_t word_count,
    uint64_t *target_words, size_t target_word_capacity,
    MVMCPowerLanczosObservableActionResult *result);

/*
 * Apply every ordered request and group byte-identical nonzero targets by
 * first appearance.  request_target_index uses
 * MVMC_POWER_LANCZOS_OBSERVABLE_ZERO_TARGET for exact-zero actions.
 * All caller-owned output arrays and result remain unmodified on failure.
 */
MVMCPowerLanczosObservableActionStatus
mvmc_power_lanczos_observable_group_targets(
    const MVMCPowerLanczosObservableLayout *layout,
    const MVMCPowerLanczosObservablePlan *plan,
    MVMCPowerLanczosObservableTargetGroupScratch *scratch,
    const uint64_t *source_words, size_t word_count,
    uint64_t *unique_target_words, size_t unique_target_word_capacity,
    int *request_target_index, size_t request_target_capacity,
    int *request_fermion_sign, size_t request_sign_capacity,
    MVMCPowerLanczosObservableTargetGroupResult *result);

const char *mvmc_power_lanczos_observable_action_status_string(
    MVMCPowerLanczosObservableActionStatus status);

#endif
