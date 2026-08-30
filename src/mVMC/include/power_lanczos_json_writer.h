#ifndef MVMC_POWER_LANCZOS_JSON_WRITER_H
#define MVMC_POWER_LANCZOS_JSON_WRITER_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "bounded_krylov_engine.h"

#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_JSON_WRITER_VERSION UINT64_C(1)
#define MVMC_POWER_LANCZOS_JSON_MAX_DEPTH 32
#define MVMC_POWER_LANCZOS_JSON_KEY_CAPACITY 128

typedef enum {
  MVMC_POWER_LANCZOS_JSON_OBJECT = 1,
  MVMC_POWER_LANCZOS_JSON_ARRAY = 2
} MVMCPowerLanczosJsonContainer;

typedef struct {
  MVMCPowerLanczosJsonContainer container;
  size_t value_count;
  int expecting_value;
  char last_key[MVMC_POWER_LANCZOS_JSON_KEY_CAPACITY];
} MVMCPowerLanczosJsonFrame;

typedef struct {
  int valid;
  int finished;
  MVMCKrylovStatus status;
  uint64_t version;
  char *buffer;
  size_t capacity;
  size_t used;
  size_t root_value_count;
  size_t depth;
  MVMCPowerLanczosJsonFrame frames[MVMC_POWER_LANCZOS_JSON_MAX_DEPTH];
} MVMCPowerLanczosJsonWriter;

int mvmc_power_lanczos_json_public_string_valid(const char *value);

MVMCKrylovStatus mvmc_power_lanczos_json_writer_initialize(
    char *buffer, size_t capacity, MVMCPowerLanczosJsonWriter *writer);

MVMCKrylovStatus mvmc_power_lanczos_json_object_begin(
    MVMCPowerLanczosJsonWriter *writer);
MVMCKrylovStatus mvmc_power_lanczos_json_object_end(
    MVMCPowerLanczosJsonWriter *writer);
MVMCKrylovStatus mvmc_power_lanczos_json_array_begin(
    MVMCPowerLanczosJsonWriter *writer);
MVMCKrylovStatus mvmc_power_lanczos_json_array_end(
    MVMCPowerLanczosJsonWriter *writer);
MVMCKrylovStatus mvmc_power_lanczos_json_key(
    MVMCPowerLanczosJsonWriter *writer, const char *key);
MVMCKrylovStatus mvmc_power_lanczos_json_null(
    MVMCPowerLanczosJsonWriter *writer);
MVMCKrylovStatus mvmc_power_lanczos_json_boolean(
    MVMCPowerLanczosJsonWriter *writer, int value);
MVMCKrylovStatus mvmc_power_lanczos_json_uint64(
    MVMCPowerLanczosJsonWriter *writer, uint64_t value);
MVMCKrylovStatus mvmc_power_lanczos_json_int64(
    MVMCPowerLanczosJsonWriter *writer, int64_t value);
MVMCKrylovStatus mvmc_power_lanczos_json_double(
    MVMCPowerLanczosJsonWriter *writer, double value);
MVMCKrylovStatus mvmc_power_lanczos_json_string(
    MVMCPowerLanczosJsonWriter *writer, const char *value);

/* Finish appends one newline and a trailing NUL. output_size includes newline. */
MVMCKrylovStatus mvmc_power_lanczos_json_finish(
    MVMCPowerLanczosJsonWriter *writer, size_t *output_size);

#endif /* bounded engine */

#endif /* MVMC_POWER_LANCZOS_JSON_WRITER_H */
