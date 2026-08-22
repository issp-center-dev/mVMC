#include "power_lanczos_observable_action.h"

#include <limits.h>
#include <stdint.h>
#include <string.h>

#define MVMC_OBSERVABLE_MAX_ORBITALS \
  (2 * MVMC_POWER_LANCZOS_OBSERVABLE_MAX_NSITE)
#define MVMC_OBSERVABLE_MAX_WORDS \
  ((MVMC_OBSERVABLE_MAX_ORBITALS + 63) / 64)

static uint64_t PopcountWord(uint64_t value) {
  uint64_t count = 0;
  while (value != 0) {
    value &= value - UINT64_C(1);
    ++count;
  }
  return count;
}

static int ConfigurationBit(const uint64_t *words, size_t orbital) {
  return (int)((words[orbital / 64] >> (orbital % 64)) & UINT64_C(1));
}

static void SetConfigurationBit(uint64_t *words, size_t orbital,
                                int value) {
  const uint64_t mask = UINT64_C(1) << (orbital % 64);
  if (value) {
    words[orbital / 64] |= mask;
  } else {
    words[orbital / 64] &= ~mask;
  }
}

static uint64_t CountOccupiedBefore(const uint64_t *words,
                                    size_t orbital) {
  const size_t full_words = orbital / 64;
  const size_t partial_bits = orbital % 64;
  uint64_t count = 0;
  size_t word;
  for (word = 0; word < full_words; ++word) {
    count += PopcountWord(words[word]);
  }
  if (partial_bits != 0) {
    const uint64_t mask =
        (UINT64_C(1) << partial_bits) - UINT64_C(1);
    count += PopcountWord(words[full_words] & mask);
  }
  return count;
}

static int ValidLayout(const MVMCPowerLanczosObservableLayout *layout) {
  if (layout == NULL || layout->nsite <= 0 ||
      layout->nsite > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_NSITE ||
      layout->up_electron_count < 0 ||
      layout->up_electron_count > layout->nsite ||
      layout->down_electron_count < 0 ||
      layout->down_electron_count > layout->nsite ||
      (layout->pure_spin != 0 && layout->pure_spin != 1)) {
    return 0;
  }
  return !layout->pure_spin ||
         layout->up_electron_count + layout->down_electron_count ==
             layout->nsite;
}

size_t mvmc_power_lanczos_observable_word_count(
    const MVMCPowerLanczosObservableLayout *layout) {
  if (!ValidLayout(layout)) return 0;
  return ((size_t)2 * (size_t)layout->nsite + 63) / 64;
}

MVMCPowerLanczosObservableActionStatus
mvmc_power_lanczos_observable_validate_configuration(
    const MVMCPowerLanczosObservableLayout *layout,
    const uint64_t *configuration_words, size_t word_count) {
  size_t site;
  size_t up_count = 0;
  size_t down_count = 0;
  size_t expected_word_count;
  size_t used_bits;
  if (!ValidLayout(layout)) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_LAYOUT;
  }
  expected_word_count = mvmc_power_lanczos_observable_word_count(layout);
  if (configuration_words == NULL || word_count != expected_word_count) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_ARGUMENT;
  }
  used_bits = ((size_t)2 * (size_t)layout->nsite) % 64;
  if (used_bits != 0) {
    const uint64_t used_mask =
        (UINT64_C(1) << used_bits) - UINT64_C(1);
    if ((configuration_words[word_count - 1] & ~used_mask) != 0) {
      return MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_SOURCE;
    }
  }
  for (site = 0; site < (size_t)layout->nsite; ++site) {
    const int up = ConfigurationBit(configuration_words, site);
    const int down = ConfigurationBit(
        configuration_words, site + (size_t)layout->nsite);
    up_count += (size_t)up;
    down_count += (size_t)down;
    if (layout->pure_spin && up + down != 1) {
      return MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_SOURCE;
    }
  }
  if (up_count != (size_t)layout->up_electron_count ||
      down_count != (size_t)layout->down_electron_count) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_SOURCE;
  }
  return MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK;
}

static int ValidRecord(const MVMCPowerLanczosObservableLayout *layout,
                       const MVMCPowerLanczosObservableRecord *record) {
  int pair;
  int expected_width;
  if (record == NULL || record->family < 0 ||
      record->family >= MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT) {
    return 0;
  }
  expected_width =
      record->family == MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS ? 4 : 8;
  if (record->row_width != expected_width) return 0;
  for (pair = 0; pair < expected_width / 2; ++pair) {
    const int site = record->raw_indices[2 * pair];
    const int spin = record->raw_indices[2 * pair + 1];
    if (site < 0 || site >= layout->nsite || spin < 0 || spin > 1) {
      return 0;
    }
  }
  if (record->raw_indices[1] != record->raw_indices[3]) return 0;
  return expected_width == 4 ||
         record->raw_indices[5] == record->raw_indices[7];
}

static MVMCPowerLanczosObservableActionStatus ValidateTarget(
    const MVMCPowerLanczosObservableLayout *layout,
    const uint64_t *configuration_words, size_t word_count) {
  const MVMCPowerLanczosObservableActionStatus status =
      mvmc_power_lanczos_observable_validate_configuration(
          layout, configuration_words, word_count);
  if (status == MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK) return status;
  if (status == MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_SOURCE) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_TARGET;
  }
  return status;
}

MVMCPowerLanczosObservableActionStatus
mvmc_power_lanczos_observable_apply(
    const MVMCPowerLanczosObservableLayout *layout,
    const MVMCPowerLanczosObservableRecord *record,
    const uint64_t *source_words, size_t word_count,
    uint64_t *target_words, size_t target_word_capacity,
    MVMCPowerLanczosObservableActionResult *result) {
  uint64_t scratch[MVMC_OBSERVABLE_MAX_WORDS];
  MVMCPowerLanczosObservableActionResult local_result;
  MVMCPowerLanczosObservableActionStatus status;
  size_t reverse_index;
  size_t operator_count;
  int sign = 1;
  memset(&local_result, 0, sizeof(local_result));
  local_result.status =
      MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_ARGUMENT;
  if (result == NULL || target_words == NULL ||
      target_word_capacity < word_count) {
    return local_result.status;
  }
  status = mvmc_power_lanczos_observable_validate_configuration(
      layout, source_words, word_count);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK) {
    return status;
  }
  if (!ValidRecord(layout, record)) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_RECORD;
  }
  if (word_count > MVMC_OBSERVABLE_MAX_WORDS) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_RESOURCE_LIMIT;
  }
  memcpy(scratch, source_words, word_count * sizeof(scratch[0]));
  operator_count = (size_t)record->row_width / 2;
  for (reverse_index = operator_count; reverse_index > 0; --reverse_index) {
    const size_t pair = reverse_index - 1;
    const int site = record->raw_indices[2 * pair];
    const int spin = record->raw_indices[2 * pair + 1];
    const size_t orbital =
        (size_t)site + (size_t)spin * (size_t)layout->nsite;
    const int creation = (pair % 2) == 0;
    const int occupied = ConfigurationBit(scratch, orbital);
    if ((creation && occupied) || (!creation && !occupied)) {
      local_result.valid = 1;
      local_result.status = MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK;
      local_result.applied = 0;
      local_result.fermion_sign = 0;
      *result = local_result;
      return local_result.status;
    }
    if ((CountOccupiedBefore(scratch, orbital) & UINT64_C(1)) != 0) {
      sign = -sign;
    }
    SetConfigurationBit(scratch, orbital, creation);
  }
  status = ValidateTarget(layout, scratch, word_count);
  if (status == MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_TARGET &&
      layout->pure_spin) {
    local_result.valid = 1;
    local_result.status = MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK;
    local_result.applied = 0;
    local_result.fermion_sign = 0;
    *result = local_result;
    return local_result.status;
  }
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK) return status;
  memcpy(target_words, scratch, word_count * sizeof(target_words[0]));
  local_result.valid = 1;
  local_result.status = MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK;
  local_result.applied = 1;
  local_result.fermion_sign = sign;
  *result = local_result;
  return local_result.status;
}

static uint64_t ConfigurationHash(const uint64_t *words,
                                  size_t word_count) {
  uint64_t hash = UINT64_C(1469598103934665603);
  size_t word;
  for (word = 0; word < word_count; ++word) {
    uint64_t value = words[word];
    int byte;
    for (byte = 0; byte < 8; ++byte) {
      hash ^= value & UINT64_C(0xff);
      hash *= UINT64_C(1099511628211);
      value >>= 8;
    }
  }
  return hash;
}

MVMCPowerLanczosObservableActionStatus
mvmc_power_lanczos_observable_group_targets(
    const MVMCPowerLanczosObservableLayout *layout,
    const MVMCPowerLanczosObservablePlan *plan,
    MVMCPowerLanczosObservableTargetGroupScratch *scratch,
    const uint64_t *source_words, size_t word_count,
    uint64_t *unique_target_words, size_t unique_target_word_capacity,
    int *request_target_index, size_t request_target_capacity,
    int *request_fermion_sign, size_t request_sign_capacity,
    MVMCPowerLanczosObservableTargetGroupResult *result) {
  uint64_t *local_targets;
  int *local_indices;
  int *local_signs;
  int *hash_slots;
  MVMCPowerLanczosObservableTargetGroupResult local_result;
  MVMCPowerLanczosObservableActionStatus status;
  size_t request_count;
  size_t request;
  size_t required_target_words;
  memset(&local_result, 0, sizeof(local_result));
  local_result.status =
      MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_ARGUMENT;
  if (result == NULL || plan == NULL || scratch == NULL ||
      unique_target_words == NULL ||
      request_target_index == NULL || request_fermion_sign == NULL) {
    return local_result.status;
  }
  status = mvmc_power_lanczos_observable_validate_configuration(
      layout, source_words, word_count);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK) return status;
  if (plan->nsite != layout->nsite || plan->record_count < 0 ||
      plan->record_count > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS ||
      (plan->record_count > 0 && plan->records == NULL)) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_ARGUMENT;
  }
  request_count = (size_t)plan->record_count;
  if (scratch->max_requests < request_count ||
      scratch->word_count != word_count || scratch->target_words == NULL ||
      scratch->request_target_index == NULL ||
      scratch->request_fermion_sign == NULL || scratch->hash_slots == NULL ||
      scratch->hash_capacity == 0 ||
      (scratch->hash_capacity & (scratch->hash_capacity - 1)) != 0 ||
      request_count > scratch->hash_capacity / 2) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_ARGUMENT;
  }
  if (request_count > SIZE_MAX / word_count) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_RESOURCE_LIMIT;
  }
  required_target_words = request_count * word_count;
  if (unique_target_word_capacity < required_target_words ||
      request_target_capacity < request_count ||
      request_sign_capacity < request_count) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_RESOURCE_LIMIT;
  }
  local_targets = scratch->target_words;
  local_indices = scratch->request_target_index;
  local_signs = scratch->request_fermion_sign;
  hash_slots = scratch->hash_slots;
  for (request = 0; request < scratch->hash_capacity; ++request) {
    hash_slots[request] = -1;
  }
  for (request = 0; request < request_count; ++request) {
    uint64_t target[MVMC_OBSERVABLE_MAX_WORDS];
    MVMCPowerLanczosObservableActionResult action;
    size_t slot;
    status = mvmc_power_lanczos_observable_apply(
        layout, &plan->records[request], source_words, word_count,
        target, MVMC_OBSERVABLE_MAX_WORDS, &action);
    if (status != MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK ||
        !action.valid) {
      return status != MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK
                 ? status
                 : MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INTERNAL_INVARIANT_FAILURE;
    }
    if (!action.applied) {
      local_indices[request] =
          MVMC_POWER_LANCZOS_OBSERVABLE_ZERO_TARGET;
      local_signs[request] = 0;
      continue;
    }
    ++local_result.active_request_count;
    slot = (size_t)(ConfigurationHash(target, word_count) &
                    (scratch->hash_capacity - 1));
    while (hash_slots[slot] >= 0) {
      const size_t target_index = (size_t)hash_slots[slot];
      if (memcmp(local_targets + target_index * word_count, target,
                 word_count * sizeof(target[0])) == 0) {
        break;
      }
      slot = (slot + 1) & (scratch->hash_capacity - 1);
    }
    if (hash_slots[slot] < 0) {
      const size_t target_index = local_result.unique_target_count;
      if (target_index >= request_count) {
        return MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INTERNAL_INVARIANT_FAILURE;
      }
      memcpy(local_targets + target_index * word_count, target,
             word_count * sizeof(target[0]));
      hash_slots[slot] = (int)target_index;
      ++local_result.unique_target_count;
      if (memcmp(source_words, target,
                 word_count * sizeof(target[0])) == 0) {
        local_result.source_target_present = 1;
      }
    }
    local_indices[request] = hash_slots[slot];
    local_signs[request] = action.fermion_sign;
  }
  memcpy(unique_target_words, local_targets,
         local_result.unique_target_count * word_count *
             sizeof(unique_target_words[0]));
  memcpy(request_target_index, local_indices,
         request_count * sizeof(request_target_index[0]));
  memcpy(request_fermion_sign, local_signs,
         request_count * sizeof(request_fermion_sign[0]));
  local_result.valid = 1;
  local_result.status = MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK;
  local_result.request_count = request_count;
  *result = local_result;
  return local_result.status;
}

const char *mvmc_power_lanczos_observable_action_status_string(
    MVMCPowerLanczosObservableActionStatus status) {
  switch (status) {
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_OK:
      return "ok";
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_ARGUMENT:
      return "invalid_argument";
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_LAYOUT:
      return "invalid_layout";
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_SOURCE:
      return "invalid_source";
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_RECORD:
      return "invalid_record";
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INVALID_TARGET:
      return "invalid_target";
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_RESOURCE_LIMIT:
      return "resource_limit";
    case MVMC_POWER_LANCZOS_OBSERVABLE_ACTION_INTERNAL_INVARIANT_FAILURE:
      return "internal_invariant_failure";
  }
  return "unknown";
}
