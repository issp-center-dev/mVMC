#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "power_lanczos_observable_census.h"

static int failures = 0;

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "PowerLanczosObservableCensus_Unit FAIL: ");            \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, "\n");                                                  \
      ++failures;                                                               \
    }                                                                           \
  } while (0)

static int WriteText(const char *path, const char *text) {
  FILE *stream = fopen(path, "wb");
  size_t length = strlen(text);
  if (stream == NULL) return 0;
  if (fwrite(text, 1, length, stream) != length) {
    fclose(stream);
    return 0;
  }
  return fclose(stream) == 0;
}

static unsigned int ReadBigEndianU32(const unsigned char *bytes) {
  return ((unsigned int)bytes[0] << 24) |
         ((unsigned int)bytes[1] << 16) |
         ((unsigned int)bytes[2] << 8) | (unsigned int)bytes[3];
}

static size_t FirstRecordOffset(const unsigned char *wire, size_t wire_size) {
  size_t offset = 8u + 7u * 4u;
  int family;
  for (family = 0; family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT;
       ++family) {
    unsigned int basename_length;
    unsigned int digest_length;
    if (offset + 8u > wire_size) return wire_size;
    offset += 4u;
    basename_length = ReadBigEndianU32(wire + offset);
    offset += 4u;
    if ((size_t)basename_length > wire_size - offset) return wire_size;
    offset += basename_length;
    if (offset + 4u > wire_size) return wire_size;
    digest_length = ReadBigEndianU32(wire + offset);
    offset += 4u;
    if ((size_t)digest_length > wire_size - offset) return wire_size;
    offset += digest_length;
  }
  return offset;
}

static int WriteObservable(const char *path, const char *label, int count,
                           const char *rows) {
  FILE *stream = fopen(path, "wb");
  if (stream == NULL) return 0;
  if (fprintf(stream,
              "====================\n%s %d\n====================\n"
              "observable rows\n====================\n%s",
              label, count, rows) < 0) {
    fclose(stream);
    return 0;
  }
  return fclose(stream) == 0;
}

static MVMCPowerLanczosObservableCensusStatus Build(
    int nsite, int nsite_uc, const char *one, const char *quartic_ex,
    const char *quartic, MVMCPowerLanczosObservablePlan *plan,
    char diagnostic[512]) {
  const char *paths[MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT] = {
      one, quartic_ex, quartic};
  return mvmc_power_lanczos_observable_plan_build_from_files(
      nsite, nsite_uc, paths, plan, diagnostic, 512);
}

static void CheckRejected(
    const char *label, int nsite, int nsite_uc, const char *one,
    const char *quartic_ex, const char *quartic,
    MVMCPowerLanczosObservableCensusStatus expected) {
  MVMCPowerLanczosObservablePlan plan;
  MVMCPowerLanczosObservableCensusStatus status;
  char diagnostic[512];
  mvmc_power_lanczos_observable_plan_init(&plan);
  status = Build(nsite, nsite_uc, one, quartic_ex, quartic, &plan,
                 diagnostic);
  CHECK(status == expected, "%s: status=%s expected=%s diagnostic=%s", label,
        mvmc_power_lanczos_observable_census_status_string(status),
        mvmc_power_lanczos_observable_census_status_string(expected),
        diagnostic);
  CHECK(plan.records == NULL && plan.record_count == 0,
        "%s: rejected input published a partial plan", label);
  mvmc_power_lanczos_observable_plan_destroy(&plan);
}

static void TestPositive(const char *directory) {
  char one[512];
  char quartic_ex[512];
  char quartic[512];
  char renamed_one[512];
  char renamed_quartic_ex[512];
  char renamed_quartic[512];
  char diagnostic[512];
  char rehash[65];
  MVMCPowerLanczosObservablePlan plan;
  MVMCPowerLanczosObservablePlan repeat;
  MVMCPowerLanczosObservablePlan renamed;
  MVMCPowerLanczosObservablePlan received;
  MVMCPowerLanczosObservableCensusStatus status;
  unsigned char *wire = NULL;
  size_t wire_size = 0;
  size_t packed_size = 0;
  FILE *json;
  char json_bytes[8192];
  size_t json_size;
  char original_raw_digest_first;
  char original_semantic_digest_first;
  snprintf(one, sizeof(one), "%s/greenone.def", directory);
  snprintf(quartic_ex, sizeof(quartic_ex), "%s/greentwoex.def", directory);
  snprintf(quartic, sizeof(quartic), "%s/greentwo.def", directory);
  snprintf(renamed_one, sizeof(renamed_one), "%s/one-renamed.def", directory);
  snprintf(renamed_quartic_ex, sizeof(renamed_quartic_ex),
           "%s/quartic-ex-renamed.def", directory);
  snprintf(renamed_quartic, sizeof(renamed_quartic),
           "%s/quartic-renamed.def", directory);
  CHECK(WriteObservable(one, "NCisAjs", 3,
                        "0 0 0 0\n0 0 1 0\n1 0 0 0\n"),
        "write positive one-body input");
  CHECK(WriteObservable(quartic_ex, "NCisAjsCktAlt", 2,
                        "0 0 0 0 1 1 1 1\n"
                        "0 0 1 0 2 1 3 1\n"),
        "write positive quartic-ex input");
  CHECK(WriteObservable(quartic, "NCisAjsCktAltDC", 1,
                        "0 0 0 0 1 1 1 1\n"),
        "write positive quartic input");

  mvmc_power_lanczos_observable_plan_init(&plan);
  mvmc_power_lanczos_observable_plan_init(&repeat);
  mvmc_power_lanczos_observable_plan_init(&renamed);
  mvmc_power_lanczos_observable_plan_init(&received);
  status = Build(4, 2, one, quartic_ex, quartic, &plan, diagnostic);
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK,
        "positive status=%s diagnostic=%s",
        mvmc_power_lanczos_observable_census_status_string(status),
        diagnostic);
  CHECK(plan.schema_version == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_VERSION &&
            plan.nsite == 4 && plan.nsite_uc == 2,
        "positive plan metadata");
  CHECK(plan.record_count == 6 && plan.family_count[0] == 3 &&
            plan.family_count[1] == 2 && plan.family_count[2] == 1,
        "positive family counts");
  CHECK(strcmp(plan.files[0].raw_file_sha256,
               "01353d6829f5082146a9dc45d7813dbd70ad5a77fa38cebe4e8d8a92f1f50226") ==
            0 &&
            strcmp(plan.files[1].raw_file_sha256,
                   "5532056934c96a983ec26725a7aab79bd923cbd3375427a62600624c3ef1f234") ==
                0 &&
            strcmp(plan.files[2].raw_file_sha256,
                   "604cb18854abb2c1a39ebb8ff3c844a6bacad9eaeb86d8a07ee2c1cb4e712339") ==
                0,
        "raw input byte SHA-256 values");
  CHECK(plan.records != NULL, "positive record allocation");
  if (plan.records != NULL) {
    CHECK(strcmp(plan.records[0].canonical_operator_id, "cdag_0_c_0") == 0 &&
              plan.records[0].adjoint_class ==
                  MVMC_POWER_LANCZOS_OBSERVABLE_SELF_ADJOINT,
          "onsite self-adjoint classification");
    CHECK(plan.records[1].adjoint_class ==
              MVMC_POWER_LANCZOS_OBSERVABLE_ADJOINT_REQUESTED &&
              plan.records[2].adjoint_class ==
                  MVMC_POWER_LANCZOS_OBSERVABLE_ADJOINT_REQUESTED,
          "one-body requested adjoint pair");
    CHECK(plan.records[3].family ==
              MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX &&
              plan.records[5].family ==
                  MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALT &&
              strcmp(plan.records[3].canonical_operator_id,
                     plan.records[5].canonical_operator_id) == 0 &&
              plan.records[3].adjoint_class ==
                  MVMC_POWER_LANCZOS_OBSERVABLE_SELF_ADJOINT,
          "quartic family identity is separate from action identity");
    plan.records[0].raw_indices[4] = 1;
    status = mvmc_power_lanczos_observable_plan_wire_size(&plan, &wire_size);
    CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT,
          "nonzero in-memory index padding must be rejected");
    plan.records[0].raw_indices[4] = 0;
  }
  CHECK(strcmp(plan.raw_observable_census_sha256,
               "c1bccf0f44655c77feb9398debf4c2068658337e26632d473efc9e731856e95f") ==
            0,
        "canonical census SHA-256 exact value");
  CHECK(strcmp(plan.semantic_observable_census_sha256,
               "a6366568b59651701ebb2cbc1c4a51914b6a58356acae7d8c794a2d72abfaf22") ==
            0,
        "semantic census SHA-256 exact value");
  status = mvmc_power_lanczos_observable_plan_rehash(&plan, rehash,
                                                     sizeof(rehash));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
            strcmp(rehash, plan.raw_observable_census_sha256) == 0,
        "canonical census rehash");
  status = mvmc_power_lanczos_observable_plan_semantic_rehash(
      &plan, rehash, sizeof(rehash));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
            strcmp(rehash, plan.semantic_observable_census_sha256) == 0,
        "semantic census rehash");

  status = Build(4, 2, one, quartic_ex, quartic, &repeat, diagnostic);
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
            strcmp(repeat.raw_observable_census_sha256,
                   plan.raw_observable_census_sha256) == 0,
        "repeated parse is deterministic");

  CHECK(WriteObservable(renamed_one, "NCisAjs", 3,
                        "0 0 0 0\n0 0 1 0\n1 0 0 0\n"),
        "write renamed one-body input");
  CHECK(WriteObservable(renamed_quartic_ex, "NCisAjsCktAlt", 2,
                        "0 0 0 0 1 1 1 1\n"
                        "0 0 1 0 2 1 3 1\n"),
        "write renamed quartic-ex input");
  CHECK(WriteObservable(renamed_quartic, "NCisAjsCktAltDC", 1,
                        "0 0 0 0 1 1 1 1\n"),
        "write renamed quartic input");
  status = Build(4, 2, renamed_one, renamed_quartic_ex, renamed_quartic,
                 &renamed, diagnostic);
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
            strcmp(renamed.raw_observable_census_sha256,
                   plan.raw_observable_census_sha256) != 0 &&
            strcmp(renamed.semantic_observable_census_sha256,
                   plan.semantic_observable_census_sha256) == 0,
        "semantic census identity must ignore raw file basenames");

  status = mvmc_power_lanczos_observable_plan_wire_size(&plan, &wire_size);
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK && wire_size > 0,
        "wire size status=%s",
        mvmc_power_lanczos_observable_census_status_string(status));
  wire = (unsigned char *)malloc(wire_size);
  CHECK(wire != NULL, "allocate wire buffer");
  if (wire != NULL) {
    status = mvmc_power_lanczos_observable_plan_pack(
        &plan, wire, wire_size, &packed_size);
    CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
              packed_size == wire_size,
          "pack wire status=%s size=%zu expected=%zu",
          mvmc_power_lanczos_observable_census_status_string(status),
          packed_size, wire_size);
    status = mvmc_power_lanczos_observable_plan_unpack(
        wire, wire_size, &received, diagnostic, sizeof(diagnostic));
    CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
              received.record_count == plan.record_count &&
              strcmp(received.raw_observable_census_sha256,
                     plan.raw_observable_census_sha256) == 0 &&
              strcmp(received.semantic_observable_census_sha256,
                     plan.semantic_observable_census_sha256) == 0,
          "wire round trip status=%s diagnostic=%s",
          mvmc_power_lanczos_observable_census_status_string(status),
          diagnostic);
    mvmc_power_lanczos_observable_plan_destroy(&received);

    {
      const size_t first_record = FirstRecordOffset(wire, wire_size);
      const size_t first_padding_byte = first_record + 7u * 4u + 3u;
      CHECK(first_record < wire_size && first_padding_byte < wire_size,
            "locate first one-body wire padding");
      if (first_padding_byte < wire_size) {
        wire[first_padding_byte] = 1;
        status = mvmc_power_lanczos_observable_plan_unpack(
            wire, wire_size, &received, diagnostic, sizeof(diagnostic));
        CHECK(status ==
                  MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT &&
                  received.records == NULL,
              "nonzero wire index padding status=%s diagnostic=%s",
              mvmc_power_lanczos_observable_census_status_string(status),
              diagnostic);
        wire[first_padding_byte] = 0;
      }
    }

    wire[wire_size - 1] = wire[wire_size - 1] == '0' ? '1' : '0';
    status = mvmc_power_lanczos_observable_plan_unpack(
        wire, wire_size, &received, diagnostic, sizeof(diagnostic));
    CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DIGEST_MISMATCH &&
              received.records == NULL,
          "wire digest mutation status=%s diagnostic=%s",
          mvmc_power_lanczos_observable_census_status_string(status),
          diagnostic);
    status = mvmc_power_lanczos_observable_plan_unpack(
        wire, wire_size - 1, &received, diagnostic, sizeof(diagnostic));
    CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT &&
              received.records == NULL,
          "truncated wire status=%s diagnostic=%s",
          mvmc_power_lanczos_observable_census_status_string(status),
          diagnostic);
    free(wire);
  }

  json = tmpfile();
  CHECK(json != NULL, "create canonical JSON scratch stream");
  if (json != NULL) {
    status = mvmc_power_lanczos_observable_plan_write_canonical_json(&plan,
                                                                    json);
    CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK,
          "write canonical JSON status=%s",
          mvmc_power_lanczos_observable_census_status_string(status));
    CHECK(fseek(json, 0, SEEK_SET) == 0, "rewind canonical JSON");
    json_size = fread(json_bytes, 1, sizeof(json_bytes) - 1, json);
    json_bytes[json_size] = '\0';
    CHECK(strstr(json_bytes,
                 "{\"raw_file_basename\":\"greenone.def\","
                 "\"raw_file_sha256\":") == json_bytes + 1,
          "canonical field order and basename");
    CHECK(strstr(json_bytes, "\"raw_ordinal\":0,\"family\":\"cisajs\","
                             "\"raw_indices\":[0,0,0,0]") != NULL,
          "canonical ordinal/family/raw-index encoding");
    fclose(json);
  }

  original_raw_digest_first = plan.raw_observable_census_sha256[0];
  original_semantic_digest_first =
      plan.semantic_observable_census_sha256[0];
  plan.raw_observable_census_sha256[0] =
      plan.raw_observable_census_sha256[0] == '0' ? '1' : '0';
  status = mvmc_power_lanczos_observable_plan_rehash(&plan, rehash,
                                                     sizeof(rehash));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DIGEST_MISMATCH,
        "tampered digest must fail binding");
  plan.raw_observable_census_sha256[0] = original_raw_digest_first;
  plan.semantic_observable_census_sha256[0] =
      plan.semantic_observable_census_sha256[0] == '0' ? '1' : '0';
  status = mvmc_power_lanczos_observable_plan_rehash(&plan, rehash,
                                                     sizeof(rehash));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DIGEST_MISMATCH,
        "tampered semantic digest must fail complete binding");
  status = mvmc_power_lanczos_observable_plan_semantic_rehash(
      &plan, rehash, sizeof(rehash));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DIGEST_MISMATCH,
        "tampered semantic digest must fail semantic binding");
  plan.semantic_observable_census_sha256[0] =
      original_semantic_digest_first;
  mvmc_power_lanczos_observable_plan_destroy(&renamed);
  mvmc_power_lanczos_observable_plan_destroy(&repeat);
  mvmc_power_lanczos_observable_plan_destroy(&received);
  mvmc_power_lanczos_observable_plan_destroy(&plan);
}

static void InitializeLegacyMap(int storage[8][8], int *rows[8]) {
  int creation;
  int annihilation;
  for (creation = 0; creation < 8; ++creation) {
    rows[creation] = storage[creation];
    for (annihilation = 0; annihilation < 8; ++annihilation) {
      storage[creation][annihilation] = -1;
    }
  }
}

static void TestLegacyAugmentedPlan(void) {
  int storage[8][8];
  int repeat_storage[8][8];
  int *rows[8];
  int *repeat_rows[8];
  char diagnostic[512];
  char rehash[65];
  char original_sha[65];
  int *original_pairs;
  MVMCPowerLanczosLegacyAugmentedPlan plan;
  MVMCPowerLanczosLegacyAugmentedPlan repeat;
  MVMCPowerLanczosObservableCensusStatus status;

  InitializeLegacyMap(storage, rows);
  InitializeLegacyMap(repeat_storage, repeat_rows);
  rows[2][3] = 0;
  rows[7][1] = 1;
  rows[0][0] = 2;
  repeat_rows[2][3] = 0;
  repeat_rows[7][1] = 1;
  repeat_rows[0][0] = 2;
  mvmc_power_lanczos_legacy_augmented_plan_init(&plan);
  mvmc_power_lanczos_legacy_augmented_plan_init(&repeat);

  status = mvmc_power_lanczos_legacy_augmented_plan_build(
      4, 3, rows, &plan, diagnostic, sizeof(diagnostic));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK,
        "legacy augmented positive status=%s diagnostic=%s",
        mvmc_power_lanczos_observable_census_status_string(status),
        diagnostic);
  CHECK(plan.schema_version ==
            MVMC_POWER_LANCZOS_LEGACY_AUGMENTED_PLAN_VERSION &&
            plan.nsite == 4 && plan.pair_count == 3 &&
            plan.orbital_pairs != NULL,
        "legacy augmented metadata");
  if (plan.orbital_pairs != NULL) {
    CHECK(plan.orbital_pairs[0] == 2 && plan.orbital_pairs[1] == 3 &&
              plan.orbital_pairs[2] == 7 && plan.orbital_pairs[3] == 1 &&
              plan.orbital_pairs[4] == 0 && plan.orbital_pairs[5] == 0,
          "legacy augmented ordinal-major ordering");
  }
  CHECK(strcmp(plan.legacy_augmented_one_body_sha256,
               "7c70a24526992bb7a9355e048ce172c555bdba96ff4998b3b85131424af1b62d") ==
            0,
        "legacy augmented canonical SHA-256 exact value");
  status = mvmc_power_lanczos_legacy_augmented_plan_rehash(
      &plan, rehash, sizeof(rehash));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
            strcmp(rehash, plan.legacy_augmented_one_body_sha256) == 0,
        "legacy augmented rehash");

  status = mvmc_power_lanczos_legacy_augmented_plan_build(
      4, 3, repeat_rows, &repeat, diagnostic, sizeof(diagnostic));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
            strcmp(repeat.legacy_augmented_one_body_sha256,
                   plan.legacy_augmented_one_body_sha256) == 0,
        "legacy augmented repeated build is deterministic");

  strcpy(original_sha, plan.legacy_augmented_one_body_sha256);
  original_pairs = plan.orbital_pairs;
  InitializeLegacyMap(storage, rows);
  rows[0][0] = 0;
  rows[1][1] = 0;
  status = mvmc_power_lanczos_legacy_augmented_plan_build(
      4, 2, rows, &plan, diagnostic, sizeof(diagnostic));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DUPLICATE &&
            plan.orbital_pairs == original_pairs &&
            strcmp(plan.legacy_augmented_one_body_sha256, original_sha) == 0,
        "legacy duplicate rejection is transactional");

  InitializeLegacyMap(storage, rows);
  rows[0][0] = 0;
  rows[1][1] = 1;
  status = mvmc_power_lanczos_legacy_augmented_plan_build(
      4, 3, rows, &plan, diagnostic, sizeof(diagnostic));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_COUNT_MISMATCH &&
            plan.orbital_pairs == original_pairs,
        "legacy missing ordinal rejection is transactional");

  InitializeLegacyMap(storage, rows);
  rows[0][0] = 3;
  status = mvmc_power_lanczos_legacy_augmented_plan_build(
      4, 3, rows, &plan, diagnostic, sizeof(diagnostic));
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT &&
            plan.orbital_pairs == original_pairs,
        "legacy out-of-range ordinal rejection is transactional");

  if (plan.orbital_pairs != NULL) {
    plan.orbital_pairs[0] ^= 1;
    status = mvmc_power_lanczos_legacy_augmented_plan_rehash(
        &plan, rehash, sizeof(rehash));
    CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DIGEST_MISMATCH,
          "legacy pair mutation must fail digest binding");
    plan.orbital_pairs[0] ^= 1;
  }
  mvmc_power_lanczos_legacy_augmented_plan_destroy(&repeat);
  mvmc_power_lanczos_legacy_augmented_plan_destroy(&plan);
}

static void TestNegative(const char *directory) {
  char path[512];
  char diagnostic[512];
  MVMCPowerLanczosObservablePlan empty_plan;
  MVMCPowerLanczosObservableCensusStatus status;
  snprintf(path, sizeof(path), "%s/negative.def", directory);

  mvmc_power_lanczos_observable_plan_init(&empty_plan);
  CHECK(WriteObservable(path, "NCisAjs", 0, ""),
        "write valid empty family");
  status = Build(4, 2, path, NULL, NULL, &empty_plan, diagnostic);
  CHECK(status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
            empty_plan.record_count == 0 && empty_plan.records == NULL &&
            strcmp(empty_plan.raw_observable_census_sha256,
                   "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945") ==
                0,
        "valid empty census status=%s diagnostic=%s",
        mvmc_power_lanczos_observable_census_status_string(status),
        diagnostic);
  mvmc_power_lanczos_observable_plan_destroy(&empty_plan);

  CHECK(WriteObservable(path, "NCisAjs", -1, ""), "write negative count");
  CheckRejected("negative count", 4, 2, path, NULL, NULL,
                MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_NEGATIVE_COUNT);

  CHECK(WriteText(path,
                  "====\nNCisAjs nope\n====\nrows\n====\n0 0 0 0\n"),
        "write malformed count");
  CheckRejected("malformed count", 4, 2, path, NULL, NULL,
                MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_MALFORMED_COUNT);

  CHECK(WriteObservable(path, "NCisAjs", 1, "0 0 0\n"),
        "write short row");
  CheckRejected("short row", 4, 2, path, NULL, NULL,
                MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ROW_WIDTH);

  CHECK(WriteObservable(path, "NCisAjs", 1, "0 0 0 0 7\n"),
        "write wide row");
  CheckRejected("wide row", 4, 2, path, NULL, NULL,
                MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ROW_WIDTH);

  CHECK(WriteObservable(path, "NCisAjs", 1, "4 0 0 0\n"),
        "write site range row");
  CheckRejected("site range", 4, 2, path, NULL, NULL,
                MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SITE_RANGE);

  CHECK(WriteObservable(path, "NCisAjs", 1, "0 2 0 2\n"),
        "write spin range row");
  CheckRejected("spin range", 4, 2, path, NULL, NULL,
                MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SPIN_RANGE);

  CHECK(WriteObservable(path, "NCisAjs", 1, "0 0 1 1\n"),
        "write spin-changing row");
  CheckRejected("spin nonconserving", 4, 2, path, NULL, NULL,
                MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SPIN_NONCONSERVING);

  CHECK(WriteObservable(path, "NCisAjs", 2,
                        "0 0 1 0\n0 0 1 0\n"),
        "write duplicate rows");
  CheckRejected("ordered duplicate", 4, 2, path, NULL, NULL,
                MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DUPLICATE);

  CHECK(WriteObservable(path, "NCisAjs", 2, "0 0 1 0\n"),
        "write missing row");
  CheckRejected("declared row missing", 4, 2, path, NULL, NULL,
                MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_COUNT_MISMATCH);

  CHECK(WriteObservable(path, "NCisAjs", 1,
                        "0 0 1 0\n1 0 0 0\n"),
        "write extra row");
  CheckRejected("undeclared extra row", 4, 2, path, NULL, NULL,
                MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_COUNT_MISMATCH);

  CHECK(WriteObservable(path, "NCisAjs", 0, "0 0 1 0\n"),
        "write zero-count file with data");
  CheckRejected("zero count with data", 4, 2, path, NULL, NULL,
                MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_COUNT_MISMATCH);

  CHECK(WriteObservable(path, "NCisAjs", 9, ""),
        "write formula-limit header");
  CheckRejected("family formula cap", 2, 2, path, NULL, NULL,
                MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_RESOURCE_LIMIT);

  CheckRejected("Nsite scope cap", 17, 1, NULL, NULL, NULL,
                MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_RESOURCE_LIMIT);
  CheckRejected("NsiteUC relation", 4, 5, NULL, NULL, NULL,
                MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_RESOURCE_LIMIT);
}

int main(void) {
  char directory[] = "power_lanczos_observable_census.XXXXXX";
  char path[512];
  if (mkdtemp(directory) == NULL) {
    fprintf(stderr,
            "PowerLanczosObservableCensus_Unit FAIL: cannot create scratch directory\n");
    return EXIT_FAILURE;
  }
  TestPositive(directory);
  TestNegative(directory);
  TestLegacyAugmentedPlan();

  snprintf(path, sizeof(path), "%s/greenone.def", directory);
  unlink(path);
  snprintf(path, sizeof(path), "%s/greentwoex.def", directory);
  unlink(path);
  snprintf(path, sizeof(path), "%s/greentwo.def", directory);
  unlink(path);
  snprintf(path, sizeof(path), "%s/negative.def", directory);
  unlink(path);
  snprintf(path, sizeof(path), "%s/one-renamed.def", directory);
  unlink(path);
  snprintf(path, sizeof(path), "%s/quartic-ex-renamed.def", directory);
  unlink(path);
  snprintf(path, sizeof(path), "%s/quartic-renamed.def", directory);
  unlink(path);
  rmdir(directory);

  if (failures != 0) {
    fprintf(stderr, "PowerLanczosObservableCensus_Unit: %d failure(s)\n",
            failures);
    return EXIT_FAILURE;
  }
  printf("PowerLanczosObservableCensus_Unit: PASS\n");
  return EXIT_SUCCESS;
}
