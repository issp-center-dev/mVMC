#include "power_lanczos_json_writer.h"

#include <locale.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(condition, ...)                                                \
  do {                                                                       \
    if (!(condition)) {                                                      \
      fprintf(stderr, "PowerLanczosJsonWriter_Unit FAIL: ");               \
      fprintf(stderr, __VA_ARGS__);                                          \
      fprintf(stderr, " (line %d)\n", __LINE__);                           \
      ++failures;                                                            \
    }                                                                        \
  } while (0)

static void test_canonical_record(void) {
  static const char expected[] =
      "{\"a\":[null,true,1,1.25,\"x\\\"y\\\\z\"],"
      "\"b\":{\"a\":-2,\"z\":false}}\n";
  char buffer[256];
  MVMCPowerLanczosJsonWriter writer;
  size_t size = 0;
  CHECK(mvmc_power_lanczos_json_writer_initialize(
            buffer, sizeof(buffer), &writer) == MVMC_KRYLOV_STATUS_OK &&
            writer.valid && writer.version ==
                MVMC_POWER_LANCZOS_JSON_WRITER_VERSION,
        "writer initialize");
  CHECK(mvmc_power_lanczos_json_object_begin(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_key(&writer, "a") ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_array_begin(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_null(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_boolean(&writer, 1) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_uint64(&writer, UINT64_C(1)) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_double(&writer, 1.25) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_string(&writer, "x\"y\\z") ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_array_end(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_key(&writer, "b") ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_object_begin(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_key(&writer, "a") ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_int64(&writer, INT64_C(-2)) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_key(&writer, "z") ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_boolean(&writer, 0) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_object_end(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_object_end(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_finish(&writer, &size) ==
                MVMC_KRYLOV_STATUS_OK &&
            size == sizeof(expected) - 1 && strcmp(buffer, expected) == 0,
        "canonical JSON bytes: %s", buffer);
}

static void test_number_boundaries(void) {
  char first[256];
  char second[256];
  size_t first_size = 0;
  size_t second_size = 0;
  MVMCPowerLanczosJsonWriter writer;
  CHECK(mvmc_power_lanczos_json_writer_initialize(
            first, sizeof(first), &writer) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_array_begin(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_double(&writer, -0.0) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_double(&writer, 0x1.0p-1022) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_double(&writer, 0x1.fffffffffffffp+1023) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_uint64(&writer, UINT64_MAX) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_int64(&writer, INT64_MIN) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_array_end(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_finish(&writer, &first_size) ==
                MVMC_KRYLOV_STATUS_OK,
        "finite number boundaries");
  CHECK(strncmp(first, "[0,", 3) == 0 && strchr(first, ',') != NULL,
        "negative zero normalization");
  CHECK(mvmc_power_lanczos_json_writer_initialize(
            second, sizeof(second), &writer) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_array_begin(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_double(&writer, -0.0) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_double(&writer, 0x1.0p-1022) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_double(&writer, 0x1.fffffffffffffp+1023) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_uint64(&writer, UINT64_MAX) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_int64(&writer, INT64_MIN) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_array_end(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_finish(&writer, &second_size) ==
                MVMC_KRYLOV_STATUS_OK &&
            first_size == second_size &&
            memcmp(first, second, first_size + 1) == 0,
        "deterministic finite rendering");
}

static void test_locale_independence(void) {
  static const char *candidates[] = {"de_DE.UTF-8", "de_DE.utf8",
                                     "fr_FR.UTF-8", "fr_FR.utf8", NULL};
  char original[128];
  char buffer[64];
  const char *current;
  const char *selected = NULL;
  size_t index;
  size_t size = 0;
  MVMCPowerLanczosJsonWriter writer;

  current = setlocale(LC_NUMERIC, NULL);
  CHECK(current != NULL && strlen(current) < sizeof(original),
        "query original numeric locale");
  if (current == NULL || strlen(current) >= sizeof(original)) return;
  memcpy(original, current, strlen(current) + 1);
  for (index = 0; candidates[index] != NULL; ++index) {
    if (setlocale(LC_NUMERIC, candidates[index]) != NULL) {
      selected = candidates[index];
      break;
    }
  }
  if (selected != NULL) {
    CHECK(mvmc_power_lanczos_json_writer_initialize(
              buffer, sizeof(buffer), &writer) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_json_double(&writer, 1.25) ==
                  MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_json_finish(&writer, &size) ==
                  MVMC_KRYLOV_STATUS_OK &&
              size == 5 && strcmp(buffer, "1.25\n") == 0,
          "locale-independent double under %s: %s", selected, buffer);
  }
  CHECK(setlocale(LC_NUMERIC, original) != NULL,
        "restore original numeric locale");
}

static void test_public_strings(void) {
  static const char *invalid[] = {
      "/Users/person/build", "/home/person/build", "~/private",
      "C:\\private\\build", "prefix D:/private/build", "\\\\server\\share",
      "prefix \\Users\\person\\build", "file://private", "Dropbox/project",
      "token=abc", "SECRET:value", "Bearer abc", "line\nvalue", NULL};
  size_t index;
  CHECK(mvmc_power_lanczos_json_public_string_valid("") &&
            mvmc_power_lanczos_json_public_string_valid(
                "power-lanczos-p6a-production-contract-v3") &&
            mvmc_power_lanczos_json_public_string_valid(
                "2026-08-28T01:20:00Z"),
        "valid public strings");
  for (index = 0; invalid[index] != NULL; ++index) {
    CHECK(!mvmc_power_lanczos_json_public_string_valid(invalid[index]),
          "unsafe public string accepted: %s", invalid[index]);
  }
}

static void test_fail_closed(void) {
  char buffer[64];
  MVMCPowerLanczosJsonWriter writer;
  size_t size = 99;
  CHECK(mvmc_power_lanczos_json_writer_initialize(
            buffer, sizeof(buffer), &writer) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_object_begin(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_key(&writer, "b") ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_null(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_key(&writer, "a") ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !writer.valid &&
            mvmc_power_lanczos_json_finish(&writer, &size) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            size == 0,
        "out-of-order key fail-closed");

  CHECK(mvmc_power_lanczos_json_writer_initialize(
            buffer, sizeof(buffer), &writer) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_array_begin(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_double(&writer, NAN) ==
                MVMC_KRYLOV_STATUS_NONFINITE &&
            !writer.valid &&
            mvmc_power_lanczos_json_finish(&writer, &size) ==
                MVMC_KRYLOV_STATUS_NONFINITE,
        "nonfinite number rejected");

  CHECK(mvmc_power_lanczos_json_writer_initialize(
            buffer, 8, &writer) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_string(&writer, "too-long") ==
                MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
            !writer.valid,
        "capacity failure");

  CHECK(mvmc_power_lanczos_json_writer_initialize(
            buffer, sizeof(buffer), &writer) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_object_begin(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_key(&writer, "a") ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_object_end(&writer) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "missing object value rejected");

  CHECK(mvmc_power_lanczos_json_writer_initialize(
            buffer, sizeof(buffer), &writer) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_null(&writer) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_json_finish(&writer, &size) ==
                MVMC_KRYLOV_STATUS_OK &&
            writer.finished && strcmp(buffer, "null\n") == 0 &&
            mvmc_power_lanczos_json_finish(&writer, &size) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            size == 0 && strcmp(buffer, "null\n") == 0,
        "double finish rejected without byte mutation");
}

static void test_depth_limit(void) {
  char buffer[256];
  MVMCPowerLanczosJsonWriter writer;
  size_t depth;
  CHECK(mvmc_power_lanczos_json_writer_initialize(
            buffer, sizeof(buffer), &writer) == MVMC_KRYLOV_STATUS_OK,
        "depth writer initialize");
  for (depth = 0; depth < MVMC_POWER_LANCZOS_JSON_MAX_DEPTH; ++depth) {
    CHECK(mvmc_power_lanczos_json_array_begin(&writer) ==
              MVMC_KRYLOV_STATUS_OK,
          "depth %zu", depth);
  }
  CHECK(mvmc_power_lanczos_json_array_begin(&writer) ==
            MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
        "maximum depth rejected");
}

int main(void) {
  test_canonical_record();
  test_number_boundaries();
  test_locale_independence();
  test_public_strings();
  test_fail_closed();
  test_depth_limit();
  if (failures != 0) {
    fprintf(stderr, "%d JSON writer checks failed\n", failures);
    return 1;
  }
  puts("power-Lanczos JSON writer checks passed");
  return 0;
}
