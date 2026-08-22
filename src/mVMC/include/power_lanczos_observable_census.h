#ifndef MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_H
#define MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_H

#include <stddef.h>
#include <stdio.h>

#define MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_VERSION 2
#define MVMC_POWER_LANCZOS_LEGACY_AUGMENTED_PLAN_VERSION 1
#define MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT 3
#define MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS 3584
#define MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ONE_BODY 512
#define MVMC_POWER_LANCZOS_OBSERVABLE_MAX_QUARTIC 1536
#define MVMC_POWER_LANCZOS_OBSERVABLE_MAX_NSITE 16
#define MVMC_POWER_LANCZOS_OBSERVABLE_MAX_NSITE_UC 16
#define MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ROW_WIDTH 8
#define MVMC_POWER_LANCZOS_OBSERVABLE_BASENAME_CAPACITY 256
#define MVMC_POWER_LANCZOS_OBSERVABLE_SHA256_HEX_CAPACITY 65
#define MVMC_POWER_LANCZOS_OBSERVABLE_OPERATOR_ID_CAPACITY 128

typedef enum {
  MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS = 0,
  MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX = 1,
  MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALT = 2
} MVMCPowerLanczosObservableFamily;

typedef enum {
  MVMC_POWER_LANCZOS_OBSERVABLE_SELF_ADJOINT = 0,
  MVMC_POWER_LANCZOS_OBSERVABLE_ADJOINT_REQUESTED,
  MVMC_POWER_LANCZOS_OBSERVABLE_ADJOINT_UNREQUESTED
} MVMCPowerLanczosObservableAdjointClass;

typedef enum {
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK = 0,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_IO_ERROR,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_MALFORMED_COUNT,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_NEGATIVE_COUNT,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_RESOURCE_LIMIT,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_COUNT_MISMATCH,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ROW_WIDTH,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SITE_RANGE,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SPIN_RANGE,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SPIN_NONCONSERVING,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DUPLICATE,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SIZE_OVERFLOW,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ALLOCATION_FAILURE,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SOURCE_CHANGED,
  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DIGEST_MISMATCH
} MVMCPowerLanczosObservableCensusStatus;

typedef struct {
  char raw_file_basename[MVMC_POWER_LANCZOS_OBSERVABLE_BASENAME_CAPACITY];
  char raw_file_sha256[MVMC_POWER_LANCZOS_OBSERVABLE_SHA256_HEX_CAPACITY];
  int signed_count;
} MVMCPowerLanczosObservableFile;

typedef struct {
  int raw_ordinal; /* Zero-based and stable within the raw family file. */
  MVMCPowerLanczosObservableFamily family;
  int row_width;
  int raw_indices[MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ROW_WIDTH];
  char canonical_operator_id[
      MVMC_POWER_LANCZOS_OBSERVABLE_OPERATOR_ID_CAPACITY];
  MVMCPowerLanczosObservableAdjointClass adjoint_class;
} MVMCPowerLanczosObservableRecord;

typedef struct {
  int schema_version;
  int nsite;
  int nsite_uc;
  int family_count[MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT];
  int record_count;
  MVMCPowerLanczosObservableFile
      files[MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT];
  MVMCPowerLanczosObservableRecord *records;
  char raw_observable_census_sha256[
      MVMC_POWER_LANCZOS_OBSERVABLE_SHA256_HEX_CAPACITY];
  /* File-name/byte-format independent identity for result certificates. */
  char semantic_observable_census_sha256[
      MVMC_POWER_LANCZOS_OBSERVABLE_SHA256_HEX_CAPACITY];
} MVMCPowerLanczosObservablePlan;

/*
 * Legacy indirect-evaluator identity.  This is deliberately a different
 * object and digest schema from the ordered raw observable plan above.
 */
typedef struct {
  int schema_version;
  int nsite;
  int pair_count;
  int *orbital_pairs; /* ordinal-major [pair_count][creation, annihilation] */
  char legacy_augmented_one_body_sha256[
      MVMC_POWER_LANCZOS_OBSERVABLE_SHA256_HEX_CAPACITY];
} MVMCPowerLanczosLegacyAugmentedPlan;

/* Initialize once before build/unpack, and destroy once after final use. */
void mvmc_power_lanczos_observable_plan_init(
    MVMCPowerLanczosObservablePlan *plan);

void mvmc_power_lanczos_observable_plan_destroy(
    MVMCPowerLanczosObservablePlan *plan);

void mvmc_power_lanczos_legacy_augmented_plan_init(
    MVMCPowerLanczosLegacyAugmentedPlan *plan);

void mvmc_power_lanczos_legacy_augmented_plan_destroy(
    MVMCPowerLanczosLegacyAugmentedPlan *plan);

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_legacy_augmented_plan_build(
    int nsite, int pair_count, int *const *pair_to_ordinal,
    MVMCPowerLanczosLegacyAugmentedPlan *plan, char *diagnostic,
    size_t diagnostic_capacity);

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_legacy_augmented_plan_rehash(
    const MVMCPowerLanczosLegacyAugmentedPlan *plan, char *actual_sha256,
    size_t actual_sha256_capacity);

/*
 * The paths use family order cisajs, cisajscktaltex, cisajscktalt.
 * A NULL or empty path denotes an empty family.  All raw validation completes
 * before the returned record array is allocated.
 */
MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_plan_build_from_files(
    int nsite, int nsite_uc,
    const char *const paths[MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT],
    MVMCPowerLanczosObservablePlan *plan, char *diagnostic,
    size_t diagnostic_capacity);

/* Rebuild canonical JSON and verify the stored SHA-256. */
MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_plan_rehash(
    const MVMCPowerLanczosObservablePlan *plan, char *actual_sha256,
    size_t actual_sha256_capacity);

/* Rebuild and verify the file-name independent certificate identity. */
MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_plan_semantic_rehash(
    const MVMCPowerLanczosObservablePlan *plan, char *actual_sha256,
    size_t actual_sha256_capacity);

/* Write the exact canonical JSON byte sequence used by plan_rehash(). */
MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_plan_write_canonical_json(
    const MVMCPowerLanczosObservablePlan *plan, FILE *stream);

/* Architecture-independent wire form used for rank-0 broadcast. */
MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_plan_wire_size(
    const MVMCPowerLanczosObservablePlan *plan, size_t *wire_size);

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_plan_pack(
    const MVMCPowerLanczosObservablePlan *plan, void *wire,
    size_t wire_capacity, size_t *wire_size);

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_plan_unpack(
    const void *wire, size_t wire_size,
    MVMCPowerLanczosObservablePlan *plan, char *diagnostic,
    size_t diagnostic_capacity);

const char *mvmc_power_lanczos_observable_family_name(
    MVMCPowerLanczosObservableFamily family);

const char *mvmc_power_lanczos_observable_adjoint_class_name(
    MVMCPowerLanczosObservableAdjointClass adjoint_class);

const char *mvmc_power_lanczos_observable_census_status_string(
    MVMCPowerLanczosObservableCensusStatus status);

#endif
