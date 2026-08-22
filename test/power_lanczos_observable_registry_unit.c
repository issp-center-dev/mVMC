#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "power_lanczos_observable_registry.h"

static int failures;

#define CHECK(condition, message)                                             \
  do {                                                                        \
    if (!(condition)) {                                                       \
      fprintf(stderr,                                                         \
              "PowerLanczosObservableRegistry_Unit FAIL: %s\n", message);  \
      ++failures;                                                             \
    }                                                                         \
  } while (0)

static void InitializeMap(int storage[4][4], int *rows[4]) {
  int creation;
  int annihilation;
  for (creation = 0; creation < 4; ++creation) {
    rows[creation] = storage[creation];
    for (annihilation = 0; annihilation < 4; ++annihilation) {
      storage[creation][annihilation] = -1;
    }
  }
}

int main(void) {
  const char *const empty_paths[MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT] = {
      NULL, NULL, NULL};
  int storage[4][4];
  int *rows[4];
  char diagnostic[512];
  char published_sha[65];
  MVMCPowerLanczosObservablePlan raw;
  MVMCPowerLanczosObservablePlan replacement;
  MVMCPowerLanczosLegacyAugmentedPlan legacy;
  MVMCPowerLanczosObservableCensusStatus status;

  mvmc_power_lanczos_observable_registry_reset();
  CHECK(mvmc_power_lanczos_observable_registry_raw() == NULL &&
            mvmc_power_lanczos_observable_registry_legacy() == NULL,
        "reset must clear both ownership domains");

  mvmc_power_lanczos_observable_plan_init(&raw);
  status = mvmc_power_lanczos_observable_plan_build_from_files(
      2, 2, empty_paths, &raw, diagnostic, sizeof(diagnostic));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK,
        "build raw source plan");
  status = mvmc_power_lanczos_observable_registry_publish_raw(&raw);
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
            raw.records == NULL && raw.schema_version == 0 &&
            mvmc_power_lanczos_observable_registry_raw() != NULL,
        "raw publication must move ownership");
  if (mvmc_power_lanczos_observable_registry_raw() != NULL) {
    strcpy(published_sha,
           mvmc_power_lanczos_observable_registry_raw()
               ->raw_observable_census_sha256);
  } else {
    published_sha[0] = '\0';
  }

  mvmc_power_lanczos_observable_plan_init(&replacement);
  status = mvmc_power_lanczos_observable_plan_build_from_files(
      3, 1, empty_paths, &replacement, diagnostic, sizeof(diagnostic));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK,
        "build replacement raw source plan");
  replacement.raw_observable_census_sha256[0] =
      replacement.raw_observable_census_sha256[0] == '0' ? '1' : '0';
  status = mvmc_power_lanczos_observable_registry_publish_raw(&replacement);
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DIGEST_MISMATCH &&
            replacement.schema_version != 0 &&
            mvmc_power_lanczos_observable_registry_raw() != NULL &&
            strcmp(mvmc_power_lanczos_observable_registry_raw()
                       ->raw_observable_census_sha256,
                   published_sha) == 0,
        "failed raw publication must preserve source and registry");
  mvmc_power_lanczos_observable_plan_destroy(&replacement);

  mvmc_power_lanczos_observable_plan_init(&replacement);
  status = mvmc_power_lanczos_observable_plan_build_from_files(
      3, 1, empty_paths, &replacement, diagnostic, sizeof(diagnostic));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK,
        "build empty-digest replacement raw source plan");
  replacement.raw_observable_census_sha256[0] = '\0';
  status = mvmc_power_lanczos_observable_registry_publish_raw(&replacement);
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DIGEST_MISMATCH &&
            replacement.schema_version != 0 &&
            mvmc_power_lanczos_observable_registry_raw() != NULL &&
            strcmp(mvmc_power_lanczos_observable_registry_raw()
                       ->raw_observable_census_sha256,
                   published_sha) == 0,
        "empty raw digest must preserve source and registry");
  mvmc_power_lanczos_observable_plan_destroy(&replacement);

  InitializeMap(storage, rows);
  rows[0][1] = 0;
  rows[3][2] = 1;
  mvmc_power_lanczos_legacy_augmented_plan_init(&legacy);
  status = mvmc_power_lanczos_legacy_augmented_plan_build(
      2, 2, rows, &legacy, diagnostic, sizeof(diagnostic));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK,
        "build legacy source plan");
  status = mvmc_power_lanczos_observable_registry_publish_legacy(&legacy);
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
            legacy.orbital_pairs == NULL && legacy.schema_version == 0 &&
            mvmc_power_lanczos_observable_registry_legacy() != NULL,
        "legacy publication must move ownership");
  CHECK(mvmc_power_lanczos_observable_registry_raw() != NULL,
        "legacy publication must not replace raw ownership");

  mvmc_power_lanczos_observable_registry_reset();
  CHECK(mvmc_power_lanczos_observable_registry_raw() == NULL &&
            mvmc_power_lanczos_observable_registry_legacy() == NULL,
        "final reset must destroy both plans");
  if (failures != 0) return EXIT_FAILURE;
  printf("PowerLanczosObservableRegistry_Unit: PASS\n");
  return EXIT_SUCCESS;
}
