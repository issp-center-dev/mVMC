#include "power_lanczos_observable_action.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

static void Require(int condition, const char *message) {
  if (!condition) {
    fprintf(stderr, "FAIL: %s\n", message);
    ++failures;
  }
}

static int InitializeGroupScratch(
    size_t max_requests, size_t word_count,
    MVMCPowerLanczosObservableTargetGroupScratch *scratch) {
  size_t target_word_count;
  size_t hash_capacity = 2;
  if (scratch == NULL || max_requests == 0 || word_count == 0 ||
      max_requests > SIZE_MAX / word_count || max_requests > SIZE_MAX / 2) {
    return 0;
  }
  while (hash_capacity < 2 * max_requests) hash_capacity *= 2;
  target_word_count = max_requests * word_count;
  memset(scratch, 0, sizeof(*scratch));
  scratch->max_requests = max_requests;
  scratch->word_count = word_count;
  scratch->hash_capacity = hash_capacity;
  scratch->target_words =
      (uint64_t *)calloc(target_word_count, sizeof(scratch->target_words[0]));
  scratch->request_target_index = (int *)calloc(
      max_requests, sizeof(scratch->request_target_index[0]));
  scratch->request_fermion_sign = (int *)calloc(
      max_requests, sizeof(scratch->request_fermion_sign[0]));
  scratch->hash_slots =
      (int *)calloc(hash_capacity, sizeof(scratch->hash_slots[0]));
  if (scratch->target_words == NULL ||
      scratch->request_target_index == NULL ||
      scratch->request_fermion_sign == NULL || scratch->hash_slots == NULL) {
    free(scratch->hash_slots);
    free(scratch->request_fermion_sign);
    free(scratch->request_target_index);
    free(scratch->target_words);
    memset(scratch, 0, sizeof(*scratch));
    return 0;
  }
  return 1;
}

static void DestroyGroupScratch(
    MVMCPowerLanczosObservableTargetGroupScratch *scratch) {
  if (scratch == NULL) return;
  free(scratch->hash_slots);
  free(scratch->request_fermion_sign);
  free(scratch->request_target_index);
  free(scratch->target_words);
  memset(scratch, 0, sizeof(*scratch));
}

static MVMCPowerLanczosObservableRecord Record(
    MVMCPowerLanczosObservableFamily family, const int *indices) {
  MVMCPowerLanczosObservableRecord record;
  int index;
  memset(&record, 0, sizeof(record));
  record.family = family;
  record.row_width =
      family == MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS ? 4 : 8;
  for (index = 0; index < record.row_width; ++index) {
    record.raw_indices[index] = indices[index];
  }
  return record;
}

static void CheckAction(
    const MVMCPowerLanczosObservableLayout *layout,
    const MVMCPowerLanczosObservableRecord *record, uint64_t source,
    int expected_applied, int expected_sign, uint64_t expected_target,
    const char *message) {
  uint64_t target = UINT64_C(0xfeedface);
  MVMCPowerLanczosObservableActionResult result;
  const MVMCPowerLanczosObservableActionStatus status =
      mvmc_power_lanczos_observable_apply(layout, record, &source, 1,
                                          &target, 1, &result);
  Require(status == MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK, message);
  Require(result.valid && result.status == status, message);
  Require(result.applied == expected_applied, message);
  Require(result.fermion_sign == expected_sign, message);
  if (expected_applied) {
    Require(target == expected_target, message);
  } else {
    Require(target == UINT64_C(0xfeedface),
            "exact-zero action modified target");
  }
}

static void TestOneBody(void) {
  const MVMCPowerLanczosObservableLayout layout = {4, 2, 1, 0};
  const uint64_t source = (UINT64_C(1) << 1) | (UINT64_C(1) << 3) |
                          (UINT64_C(1) << 4);
  const int onsite_indices[4] = {1, 0, 1, 0};
  const int negative_indices[4] = {0, 0, 3, 0};
  const int positive_indices[4] = {2, 0, 1, 0};
  const int pauli_indices[4] = {1, 0, 3, 0};
  const int empty_indices[4] = {2, 1, 2, 1};
  MVMCPowerLanczosObservableRecord onsite =
      Record(MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, onsite_indices);
  MVMCPowerLanczosObservableRecord negative =
      Record(MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, negative_indices);
  MVMCPowerLanczosObservableRecord positive =
      Record(MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, positive_indices);
  MVMCPowerLanczosObservableRecord pauli =
      Record(MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, pauli_indices);
  MVMCPowerLanczosObservableRecord empty =
      Record(MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, empty_indices);
  CheckAction(&layout, &onsite, source, 1, 1, source,
              "onsite density action");
  CheckAction(&layout, &negative, source, 1, -1,
              (UINT64_C(1) << 0) | (UINT64_C(1) << 1) |
                  (UINT64_C(1) << 4),
              "negative-parity hopping action");
  CheckAction(&layout, &positive, source, 1, 1,
              (UINT64_C(1) << 2) | (UINT64_C(1) << 3) |
                  (UINT64_C(1) << 4),
              "positive-parity hopping action");
  CheckAction(&layout, &pauli, source, 0, 0, 0,
              "Pauli-forbidden occupied create");
  CheckAction(&layout, &empty, source, 0, 0, 0,
              "Pauli-forbidden empty annihilate");
}

static void TestQuarticAndPureSpin(void) {
  const MVMCPowerLanczosObservableLayout electronic = {4, 2, 1, 0};
  const uint64_t source = (UINT64_C(1) << 1) | (UINT64_C(1) << 3) |
                          (UINT64_C(1) << 4);
  const int density_indices[8] = {1, 0, 1, 0, 0, 1, 0, 1};
  const int distinct_indices[8] = {0, 0, 3, 0, 2, 1, 0, 1};
  const int repeated_indices[8] = {0, 0, 3, 0, 0, 0, 1, 0};
  MVMCPowerLanczosObservableRecord density = Record(
      MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX, density_indices);
  MVMCPowerLanczosObservableRecord distinct = Record(
      MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX, distinct_indices);
  MVMCPowerLanczosObservableRecord repeated = Record(
      MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX, repeated_indices);
  const MVMCPowerLanczosObservableLayout pure_spin = {2, 1, 1, 1};
  const uint64_t spin_source =
      (UINT64_C(1) << 1) | (UINT64_C(1) << 2);
  const int exchange_indices[8] = {0, 0, 1, 0, 1, 1, 0, 1};
  const int projected_zero_indices[4] = {0, 0, 1, 0};
  MVMCPowerLanczosObservableRecord exchange = Record(
      MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALT, exchange_indices);
  MVMCPowerLanczosObservableRecord projected_zero = Record(
      MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, projected_zero_indices);
  CheckAction(&electronic, &density, source, 1, 1, source,
              "quartic density-density action");
  CheckAction(&electronic, &distinct, source, 1, -1,
              (UINT64_C(1) << 0) | (UINT64_C(1) << 1) |
                  (UINT64_C(1) << 6),
              "four-distinct-orbital action");
  CheckAction(&electronic, &repeated, source, 0, 0, 0,
              "quartic repeated-orbital Pauli zero");
  CheckAction(&pure_spin, &exchange, spin_source, 1, 1,
              (UINT64_C(1) << 0) | (UINT64_C(1) << 3),
              "pure-spin transverse exchange action");
  CheckAction(&pure_spin, &projected_zero, spin_source, 0, 0, 0,
              "pure-spin out-of-sector one-body action");
  {
    MVMCPowerLanczosObservablePlan plan;
    MVMCPowerLanczosObservableTargetGroupResult result;
    uint64_t target = UINT64_C(0xdecafbad);
    int target_index = 99;
    int fermion_sign = 99;
    MVMCPowerLanczosObservableTargetGroupScratch scratch;
    MVMCPowerLanczosObservableActionStatus status;
    memset(&plan, 0, sizeof(plan));
    plan.nsite = pure_spin.nsite;
    plan.record_count = 1;
    plan.records = &projected_zero;
    Require(InitializeGroupScratch(1, 1, &scratch),
            "allocate pure-spin grouping scratch");
    status = mvmc_power_lanczos_observable_group_targets(
        &pure_spin, &plan, &scratch, &spin_source, 1, &target, 1,
        &target_index, 1, &fermion_sign, 1, &result);
    Require(status == MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK &&
                result.valid,
            "pure-spin projected-zero grouping failed");
    Require(result.request_count == 1 &&
                result.active_request_count == 0 &&
                result.unique_target_count == 0 &&
                !result.source_target_present,
            "pure-spin projected-zero grouping counts mismatch");
    Require(target_index == MVMC_POWER_LANCZOS_OBSERVABLE_ZERO_TARGET &&
                fermion_sign == 0,
            "pure-spin projected-zero grouping result mismatch");
    Require(target == UINT64_C(0xdecafbad),
            "pure-spin projected-zero grouping modified target");
    DestroyGroupScratch(&scratch);
  }
}

static void TestValidationAndTransaction(void) {
  MVMCPowerLanczosObservableLayout layout = {4, 2, 1, 0};
  const uint64_t source = (UINT64_C(1) << 1) | (UINT64_C(1) << 3) |
                          (UINT64_C(1) << 4);
  const int indices[4] = {0, 0, 3, 0};
  MVMCPowerLanczosObservableRecord record =
      Record(MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, indices);
  MVMCPowerLanczosObservableActionResult result;
  uint64_t target = UINT64_C(0x12345678);
  MVMCPowerLanczosObservableActionStatus status;
  status = mvmc_power_lanczos_observable_apply(
      &layout, &record, &source, 1, &target, 0, &result);
  Require(status == MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_ARGUMENT,
          "short target capacity accepted");
  Require(target == UINT64_C(0x12345678),
          "failed action modified target");
  record.raw_indices[1] = 1;
  status = mvmc_power_lanczos_observable_apply(
      &layout, &record, &source, 1, &target, 1, &result);
  Require(status == MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_RECORD,
          "spin-nonconserving record accepted");
  Require(target == UINT64_C(0x12345678),
          "invalid record modified target");
  layout.pure_spin = 1;
  Require(mvmc_power_lanczos_observable_word_count(&layout) == 0,
          "invalid pure-spin layout accepted");
  layout = (MVMCPowerLanczosObservableLayout){4, 2, 1, 0};
  {
    const uint64_t invalid_source = source | (UINT64_C(1) << 63);
    Require(mvmc_power_lanczos_observable_validate_configuration(
                &layout, &invalid_source, 1) ==
                MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_SOURCE,
            "noncanonical high bit accepted");
  }
  Require(strcmp(mvmc_power_lanczos_observable_action_status_string(
                     MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_RESOURCE_LIMIT),
                 "resource_limit") == 0,
          "status string mismatch");
}

static void TestTargetGrouping(void) {
  const MVMCPowerLanczosObservableLayout layout = {4, 2, 1, 0};
  const uint64_t source = (UINT64_C(1) << 1) | (UINT64_C(1) << 3) |
                          (UINT64_C(1) << 4);
  const int rows[5][8] = {
      {0, 0, 3, 0, 0, 0, 0, 0},
      {2, 0, 1, 0, 0, 0, 0, 0},
      {1, 0, 1, 0, 0, 1, 0, 1},
      {1, 0, 1, 0, 0, 1, 0, 1},
      {1, 0, 3, 0, 0, 0, 0, 0}};
  MVMCPowerLanczosObservableRecord records[5];
  MVMCPowerLanczosObservablePlan plan;
  uint64_t targets[5] = {UINT64_C(0xaaaaaaaa), UINT64_C(0xaaaaaaaa),
                         UINT64_C(0xaaaaaaaa), UINT64_C(0xaaaaaaaa),
                         UINT64_C(0xaaaaaaaa)};
  int indices[5] = {99, 99, 99, 99, 99};
  int signs[5] = {99, 99, 99, 99, 99};
  MVMCPowerLanczosObservableTargetGroupResult result;
  MVMCPowerLanczosObservableTargetGroupScratch scratch;
  MVMCPowerLanczosObservableActionStatus status;
  memset(&plan, 0, sizeof(plan));
  records[0] = Record(MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, rows[0]);
  records[1] = Record(MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, rows[1]);
  records[2] = Record(MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX,
                      rows[2]);
  records[3] = Record(MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALT, rows[3]);
  records[4] = Record(MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, rows[4]);
  plan.nsite = 4;
  plan.record_count = 5;
  plan.records = records;
  Require(InitializeGroupScratch(5, 1, &scratch),
          "allocate grouping scratch");
  status = mvmc_power_lanczos_observable_group_targets(
      &layout, &plan, &scratch, &source, 1, targets, 5, indices, 5, signs, 5,
      &result);
  Require(status == MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK && result.valid,
          "target grouping failed");
  Require(result.request_count == 5 && result.active_request_count == 4 &&
              result.unique_target_count == 3 &&
              result.source_target_present,
          "target grouping counts mismatch");
  Require(indices[0] == 0 && indices[1] == 1 && indices[2] == 2 &&
              indices[3] == 2 &&
              indices[4] == MVMC_POWER_LANCZOS_OBSERVABLE_ZERO_TARGET,
          "target grouping stable indices mismatch");
  Require(signs[0] == -1 && signs[1] == 1 && signs[2] == 1 &&
              signs[3] == 1 && signs[4] == 0,
          "target grouping signs mismatch");
  Require(targets[2] == source,
          "cross-family identical action did not reuse target");
  indices[0] = 77;
  signs[0] = 88;
  targets[0] = UINT64_C(0xbbbbbbbb);
  status = mvmc_power_lanczos_observable_group_targets(
      &layout, &plan, &scratch, &source, 1, targets, 4, indices, 5, signs, 5,
      &result);
  Require(status == MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_RESOURCE_LIMIT,
          "short grouping capacity accepted");
  Require(indices[0] == 77 && signs[0] == 88 &&
              targets[0] == UINT64_C(0xbbbbbbbb),
          "failed grouping modified outputs");
  scratch.max_requests = 4;
  status = mvmc_power_lanczos_observable_group_targets(
      &layout, &plan, &scratch, &source, 1, targets, 5, indices, 5, signs, 5,
      &result);
  Require(status == MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_ARGUMENT,
          "undersized grouping scratch accepted");
  Require(indices[0] == 77 && signs[0] == 88 &&
              targets[0] == UINT64_C(0xbbbbbbbb),
          "invalid grouping scratch modified outputs");
  DestroyGroupScratch(&scratch);
}

int main(void) {
  TestOneBody();
  TestQuarticAndPureSpin();
  TestValidationAndTransaction();
  TestTargetGrouping();
  if (failures != 0) {
    fprintf(stderr, "%d observable action assertion(s) failed\n", failures);
    return EXIT_FAILURE;
  }
  puts("power_lanczos_observable_action_unit: PASS");
  return EXIT_SUCCESS;
}
