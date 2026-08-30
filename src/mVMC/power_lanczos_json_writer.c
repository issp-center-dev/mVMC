#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200809L
#endif

#include "power_lanczos_json_writer.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "power_lanczos_json_writer.c requires bounded Krylov"
#endif

#include <inttypes.h>
#include <locale.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#if defined(__APPLE__)
#include <xlocale.h>
#endif

static MVMCKrylovStatus writer_fail(MVMCPowerLanczosJsonWriter *writer,
                                    MVMCKrylovStatus status) {
  if (writer != NULL) {
    if (!writer->valid && writer->status != MVMC_KRYLOV_STATUS_OK) {
      return writer->status;
    }
    writer->valid = 0;
    writer->status = status;
  }
  return status;
}

static int writer_ready(const MVMCPowerLanczosJsonWriter *writer) {
  return writer != NULL && writer->valid && !writer->finished &&
         writer->status == MVMC_KRYLOV_STATUS_OK && writer->buffer != NULL &&
         writer->capacity > writer->used;
}

static MVMCKrylovStatus append_bytes(MVMCPowerLanczosJsonWriter *writer,
                                     const char *bytes, size_t length) {
  if (!writer_ready(writer) || bytes == NULL ||
      length > writer->capacity - writer->used - 1) {
    return writer_fail(writer, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
  }
  if (length != 0) memcpy(writer->buffer + writer->used, bytes, length);
  writer->used += length;
  writer->buffer[writer->used] = '\0';
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus append_character(MVMCPowerLanczosJsonWriter *writer,
                                         char character) {
  return append_bytes(writer, &character, 1);
}

static int ascii_printable_string(const char *value, int allow_empty) {
  const unsigned char *cursor = (const unsigned char *)value;
  if (value == NULL || (!allow_empty && *value == '\0')) return 0;
  while (*cursor != '\0') {
    if (*cursor < 0x20 || *cursor > 0x7e) return 0;
    ++cursor;
  }
  return 1;
}

static int ascii_equal_fold(unsigned char left, unsigned char right) {
  if (left >= 'A' && left <= 'Z') left = (unsigned char)(left + 'a' - 'A');
  if (right >= 'A' && right <= 'Z') {
    right = (unsigned char)(right + 'a' - 'A');
  }
  return left == right;
}

static int contains_fold(const char *value, const char *needle) {
  size_t value_length;
  size_t needle_length;
  size_t offset;
  size_t index;
  if (value == NULL || needle == NULL) return 0;
  value_length = strlen(value);
  needle_length = strlen(needle);
  if (needle_length == 0 || needle_length > value_length) return 0;
  for (offset = 0; offset <= value_length - needle_length; ++offset) {
    for (index = 0; index < needle_length; ++index) {
      if (!ascii_equal_fold((unsigned char)value[offset + index],
                            (unsigned char)needle[index])) {
        break;
      }
    }
    if (index == needle_length) return 1;
  }
  return 0;
}

static int contains_windows_absolute_path(const char *value) {
  size_t index;
  const size_t length = value == NULL ? 0 : strlen(value);
  for (index = 0; index + 2 < length; ++index) {
    const unsigned char drive = (unsigned char)value[index];
    if (((drive >= 'A' && drive <= 'Z') ||
         (drive >= 'a' && drive <= 'z')) &&
        value[index + 1] == ':' &&
        (value[index + 2] == '/' || value[index + 2] == '\\')) {
      return 1;
    }
  }
  return 0;
}

int mvmc_power_lanczos_json_public_string_valid(const char *value) {
  static const char *forbidden[] = {
      "/users/",  "/home/",     "\\users\\",  "\\home\\",
      "dropbox",  "file://",    "token=",      "token:",
      "secret=",  "secret:",    "credential=", "credential:",
      "bearer ",  NULL};
  size_t index;
  if (!ascii_printable_string(value, 1)) return 0;
  if (value[0] == '/' || value[0] == '~' ||
      (value[0] == '\\' && value[1] == '\\') ||
      contains_windows_absolute_path(value)) {
    return 0;
  }
  for (index = 0; forbidden[index] != NULL; ++index) {
    if (contains_fold(value, forbidden[index])) return 0;
  }
  return 1;
}

static MVMCKrylovStatus append_quoted(MVMCPowerLanczosJsonWriter *writer,
                                      const char *value) {
  const char *cursor;
  MVMCKrylovStatus status;
  status = append_character(writer, '"');
  for (cursor = value;
       status == MVMC_KRYLOV_STATUS_OK && *cursor != '\0'; ++cursor) {
    if (*cursor == '"' || *cursor == '\\') {
      status = append_character(writer, '\\');
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = append_character(writer, *cursor);
    }
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = append_character(writer, '"');
  }
  return status;
}

static MVMCKrylovStatus before_value(MVMCPowerLanczosJsonWriter *writer) {
  MVMCPowerLanczosJsonFrame *frame;
  if (!writer_ready(writer)) {
    return writer == NULL ? MVMC_KRYLOV_STATUS_INVALID_ARGUMENT
                          : writer->status;
  }
  if (writer->depth == 0) {
    if (writer->root_value_count != 0) {
      return writer_fail(writer, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
    }
    ++writer->root_value_count;
    return MVMC_KRYLOV_STATUS_OK;
  }
  frame = &writer->frames[writer->depth - 1];
  if (frame->container == MVMC_POWER_LANCZOS_JSON_OBJECT) {
    if (!frame->expecting_value) {
      return writer_fail(writer, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
    }
    frame->expecting_value = 0;
    ++frame->value_count;
    return MVMC_KRYLOV_STATUS_OK;
  }
  if (frame->container != MVMC_POWER_LANCZOS_JSON_ARRAY) {
    return writer_fail(writer,
                       MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE);
  }
  if (frame->value_count != 0 &&
      append_character(writer, ',') != MVMC_KRYLOV_STATUS_OK) {
    return writer->status;
  }
  ++frame->value_count;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus container_begin(
    MVMCPowerLanczosJsonWriter *writer,
    MVMCPowerLanczosJsonContainer container, char opening) {
  MVMCPowerLanczosJsonFrame *frame;
  MVMCKrylovStatus status = before_value(writer);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (writer->depth >= MVMC_POWER_LANCZOS_JSON_MAX_DEPTH) {
    return writer_fail(writer, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
  }
  status = append_character(writer, opening);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  frame = &writer->frames[writer->depth++];
  memset(frame, 0, sizeof(*frame));
  frame->container = container;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus container_end(
    MVMCPowerLanczosJsonWriter *writer,
    MVMCPowerLanczosJsonContainer container, char closing) {
  MVMCPowerLanczosJsonFrame *frame;
  if (!writer_ready(writer) || writer->depth == 0) {
    return writer_fail(writer, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
  }
  frame = &writer->frames[writer->depth - 1];
  if (frame->container != container || frame->expecting_value) {
    return writer_fail(writer, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
  }
  --writer->depth;
  return append_character(writer, closing);
}

MVMCKrylovStatus mvmc_power_lanczos_json_writer_initialize(
    char *buffer, size_t capacity, MVMCPowerLanczosJsonWriter *writer) {
  if (writer == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(writer, 0, sizeof(*writer));
  writer->status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  writer->version = MVMC_POWER_LANCZOS_JSON_WRITER_VERSION;
  if (buffer == NULL || capacity < 2) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  buffer[0] = '\0';
  writer->valid = 1;
  writer->status = MVMC_KRYLOV_STATUS_OK;
  writer->buffer = buffer;
  writer->capacity = capacity;
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_json_object_begin(
    MVMCPowerLanczosJsonWriter *writer) {
  return container_begin(writer, MVMC_POWER_LANCZOS_JSON_OBJECT, '{');
}

MVMCKrylovStatus mvmc_power_lanczos_json_object_end(
    MVMCPowerLanczosJsonWriter *writer) {
  return container_end(writer, MVMC_POWER_LANCZOS_JSON_OBJECT, '}');
}

MVMCKrylovStatus mvmc_power_lanczos_json_array_begin(
    MVMCPowerLanczosJsonWriter *writer) {
  return container_begin(writer, MVMC_POWER_LANCZOS_JSON_ARRAY, '[');
}

MVMCKrylovStatus mvmc_power_lanczos_json_array_end(
    MVMCPowerLanczosJsonWriter *writer) {
  return container_end(writer, MVMC_POWER_LANCZOS_JSON_ARRAY, ']');
}

MVMCKrylovStatus mvmc_power_lanczos_json_key(
    MVMCPowerLanczosJsonWriter *writer, const char *key) {
  MVMCPowerLanczosJsonFrame *frame;
  size_t length;
  MVMCKrylovStatus status;
  if (!writer_ready(writer) || writer->depth == 0 ||
      !ascii_printable_string(key, 0)) {
    return writer_fail(writer, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
  }
  length = strlen(key);
  if (length >= MVMC_POWER_LANCZOS_JSON_KEY_CAPACITY) {
    return writer_fail(writer, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
  }
  frame = &writer->frames[writer->depth - 1];
  if (frame->container != MVMC_POWER_LANCZOS_JSON_OBJECT ||
      frame->expecting_value ||
      (frame->value_count != 0 && strcmp(frame->last_key, key) >= 0)) {
    return writer_fail(writer, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
  }
  if (frame->value_count != 0) {
    status = append_character(writer, ',');
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
  }
  status = append_quoted(writer, key);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = append_character(writer, ':');
  }
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  memcpy(frame->last_key, key, length + 1);
  frame->expecting_value = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus literal(MVMCPowerLanczosJsonWriter *writer,
                                const char *value) {
  MVMCKrylovStatus status = before_value(writer);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  return append_bytes(writer, value, strlen(value));
}

MVMCKrylovStatus mvmc_power_lanczos_json_null(
    MVMCPowerLanczosJsonWriter *writer) {
  return literal(writer, "null");
}

MVMCKrylovStatus mvmc_power_lanczos_json_boolean(
    MVMCPowerLanczosJsonWriter *writer, int value) {
  if (value != 0 && value != 1) {
    return writer_fail(writer, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
  }
  return literal(writer, value ? "true" : "false");
}

MVMCKrylovStatus mvmc_power_lanczos_json_uint64(
    MVMCPowerLanczosJsonWriter *writer, uint64_t value) {
  char buffer[32];
  const int length = snprintf(buffer, sizeof(buffer), "%" PRIu64, value);
  if (length <= 0 || (size_t)length >= sizeof(buffer)) {
    return writer_fail(writer,
                       MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE);
  }
  return literal(writer, buffer);
}

MVMCKrylovStatus mvmc_power_lanczos_json_int64(
    MVMCPowerLanczosJsonWriter *writer, int64_t value) {
  char buffer[32];
  const int length = snprintf(buffer, sizeof(buffer), "%" PRId64, value);
  if (length <= 0 || (size_t)length >= sizeof(buffer)) {
    return writer_fail(writer,
                       MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE);
  }
  return literal(writer, buffer);
}

MVMCKrylovStatus mvmc_power_lanczos_json_double(
    MVMCPowerLanczosJsonWriter *writer, double value) {
  char buffer[64];
  int length;
  locale_t c_locale;
  locale_t previous;
  if (!isfinite(value)) {
    return writer_fail(writer, MVMC_KRYLOV_STATUS_NONFINITE);
  }
  if (value == 0.0) return literal(writer, "0");
  c_locale = newlocale(LC_NUMERIC_MASK, "C", (locale_t)0);
  if (c_locale == (locale_t)0) {
    return writer_fail(writer,
                       MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE);
  }
  previous = uselocale(c_locale);
  if (previous == (locale_t)0) {
    freelocale(c_locale);
    return writer_fail(writer,
                       MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE);
  }
  length = snprintf(buffer, sizeof(buffer), "%.17g", value);
  if (uselocale(previous) == (locale_t)0) {
    /* c_locale is still active for this thread; freeing it would dangle. */
    return writer_fail(writer,
                       MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE);
  }
  freelocale(c_locale);
  if (length <= 0 || (size_t)length >= sizeof(buffer) ||
      strchr(buffer, ',') != NULL) {
    return writer_fail(writer,
                       MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE);
  }
  return literal(writer, buffer);
}

MVMCKrylovStatus mvmc_power_lanczos_json_string(
    MVMCPowerLanczosJsonWriter *writer, const char *value) {
  MVMCKrylovStatus status;
  if (!mvmc_power_lanczos_json_public_string_valid(value)) {
    return writer_fail(writer, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
  }
  status = before_value(writer);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  return append_quoted(writer, value);
}

MVMCKrylovStatus mvmc_power_lanczos_json_finish(
    MVMCPowerLanczosJsonWriter *writer, size_t *output_size) {
  MVMCKrylovStatus status;
  if (output_size != NULL) *output_size = 0;
  if (!writer_ready(writer) || output_size == NULL || writer->depth != 0 ||
      writer->root_value_count != 1) {
    return writer_fail(writer, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
  }
  status = append_character(writer, '\n');
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  writer->finished = 1;
  *output_size = writer->used;
  return MVMC_KRYLOV_STATUS_OK;
}
