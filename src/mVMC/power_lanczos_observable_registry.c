#include "power_lanczos_observable_registry.h"

#include <string.h>

static MVMCPowerLanczosObservablePlan RawPlan;
static MVMCPowerLanczosLegacyAugmentedPlan LegacyPlan;
static int RawPlanPublished;
static int LegacyPlanPublished;

void mvmc_power_lanczos_observable_registry_reset(void) {
  mvmc_power_lanczos_observable_plan_destroy(&RawPlan);
  mvmc_power_lanczos_legacy_augmented_plan_destroy(&LegacyPlan);
  RawPlanPublished = 0;
  LegacyPlanPublished = 0;
}

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_registry_publish_raw(
    MVMCPowerLanczosObservablePlan *source) {
  MVMCPowerLanczosObservableCensusStatus status;
  if (source == NULL) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
  }
  status = mvmc_power_lanczos_observable_plan_rehash(source, NULL, 0);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  mvmc_power_lanczos_observable_plan_destroy(&RawPlan);
  RawPlan = *source;
  memset(source, 0, sizeof(*source));
  RawPlanPublished = 1;
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
}

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_registry_publish_legacy(
    MVMCPowerLanczosLegacyAugmentedPlan *source) {
  MVMCPowerLanczosObservableCensusStatus status;
  if (source == NULL) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
  }
  status = mvmc_power_lanczos_legacy_augmented_plan_rehash(source, NULL, 0);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  mvmc_power_lanczos_legacy_augmented_plan_destroy(&LegacyPlan);
  LegacyPlan = *source;
  memset(source, 0, sizeof(*source));
  LegacyPlanPublished = 1;
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
}

const MVMCPowerLanczosObservablePlan *
mvmc_power_lanczos_observable_registry_raw(void) {
  return RawPlanPublished ? &RawPlan : NULL;
}

const MVMCPowerLanczosLegacyAugmentedPlan *
mvmc_power_lanczos_observable_registry_legacy(void) {
  return LegacyPlanPublished ? &LegacyPlan : NULL;
}
