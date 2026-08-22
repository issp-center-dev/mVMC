#ifndef MVMC_POWER_LANCZOS_OBSERVABLE_REGISTRY_H
#define MVMC_POWER_LANCZOS_OBSERVABLE_REGISTRY_H

#include "power_lanczos_observable_census.h"

/*
 * Process-local ownership boundary for the two intentionally separate plans.
 * Publication validates the stored digest and moves the source object only
 * after validation succeeds.  The caller must initialize source before use.
 */
void mvmc_power_lanczos_observable_registry_reset(void);

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_registry_publish_raw(
    MVMCPowerLanczosObservablePlan *source);

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_registry_publish_legacy(
    MVMCPowerLanczosLegacyAugmentedPlan *source);

const MVMCPowerLanczosObservablePlan *
mvmc_power_lanczos_observable_registry_raw(void);

const MVMCPowerLanczosLegacyAugmentedPlan *
mvmc_power_lanczos_observable_registry_legacy(void);

#endif
