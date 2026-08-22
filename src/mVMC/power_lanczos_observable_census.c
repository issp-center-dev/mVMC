#include "power_lanczos_observable_census.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

enum { MVMC_OBSERVABLE_INPUT_LINE_CAPACITY = 4096 };

typedef struct {
  uint32_t state[8];
  uint64_t bit_count;
  unsigned char block[64];
  size_t block_size;
} SHA256Context;

typedef struct {
  MVMCPowerLanczosObservableFamily family;
  int width;
  int indices[MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ROW_WIDTH];
} ParsedRecord;

typedef struct {
  char *bytes;
  size_t size;
  size_t capacity;
} JsonBuilder;

static const unsigned char ObservableWireMagic[8] = {
    'P', '6', 'O', 'B', 'S', 'C', '1', '\0'};

static const uint32_t SHA256Constants[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u,
    0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
    0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
    0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
    0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
    0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u,
    0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
    0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
    0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
    0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u,
    0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u,
    0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
    0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
    0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};

static uint32_t RotateRight(uint32_t value, unsigned int shift) {
  return (value >> shift) | (value << (32u - shift));
}

static void SHA256Transform(SHA256Context *context,
                            const unsigned char block[64]) {
  uint32_t words[64];
  uint32_t a, b, c, d, e, f, g, h;
  int index;

  for (index = 0; index < 16; ++index) {
    const int offset = 4 * index;
    words[index] = ((uint32_t)block[offset] << 24) |
                   ((uint32_t)block[offset + 1] << 16) |
                   ((uint32_t)block[offset + 2] << 8) |
                   (uint32_t)block[offset + 3];
  }
  for (index = 16; index < 64; ++index) {
    const uint32_t s0 = RotateRight(words[index - 15], 7) ^
                        RotateRight(words[index - 15], 18) ^
                        (words[index - 15] >> 3);
    const uint32_t s1 = RotateRight(words[index - 2], 17) ^
                        RotateRight(words[index - 2], 19) ^
                        (words[index - 2] >> 10);
    words[index] = words[index - 16] + s0 + words[index - 7] + s1;
  }

  a = context->state[0];
  b = context->state[1];
  c = context->state[2];
  d = context->state[3];
  e = context->state[4];
  f = context->state[5];
  g = context->state[6];
  h = context->state[7];
  for (index = 0; index < 64; ++index) {
    const uint32_t sum0 = RotateRight(a, 2) ^ RotateRight(a, 13) ^
                          RotateRight(a, 22);
    const uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
    const uint32_t sum1 = RotateRight(e, 6) ^ RotateRight(e, 11) ^
                          RotateRight(e, 25);
    const uint32_t choice = (e & f) ^ ((~e) & g);
    const uint32_t temporary1 =
        h + sum1 + choice + SHA256Constants[index] + words[index];
    const uint32_t temporary2 = sum0 + majority;
    h = g;
    g = f;
    f = e;
    e = d + temporary1;
    d = c;
    c = b;
    b = a;
    a = temporary1 + temporary2;
  }
  context->state[0] += a;
  context->state[1] += b;
  context->state[2] += c;
  context->state[3] += d;
  context->state[4] += e;
  context->state[5] += f;
  context->state[6] += g;
  context->state[7] += h;
}

static void SHA256Init(SHA256Context *context) {
  static const uint32_t initial[8] = {
      0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
      0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u};
  memcpy(context->state, initial, sizeof(initial));
  context->bit_count = 0;
  context->block_size = 0;
}

static void SHA256Update(SHA256Context *context, const void *input,
                         size_t input_size) {
  const unsigned char *bytes = (const unsigned char *)input;
  while (input_size > 0) {
    size_t available = sizeof(context->block) - context->block_size;
    size_t take = input_size < available ? input_size : available;
    memcpy(context->block + context->block_size, bytes, take);
    context->block_size += take;
    context->bit_count += (uint64_t)take * 8u;
    bytes += take;
    input_size -= take;
    if (context->block_size == sizeof(context->block)) {
      SHA256Transform(context, context->block);
      context->block_size = 0;
    }
  }
}

static void SHA256Final(SHA256Context *context, unsigned char digest[32]) {
  uint64_t bit_count = context->bit_count;
  int index;
  context->block[context->block_size++] = 0x80u;
  if (context->block_size > 56) {
    while (context->block_size < 64) context->block[context->block_size++] = 0;
    SHA256Transform(context, context->block);
    context->block_size = 0;
  }
  while (context->block_size < 56) context->block[context->block_size++] = 0;
  for (index = 7; index >= 0; --index) {
    context->block[context->block_size++] =
        (unsigned char)(bit_count >> (8 * index));
  }
  SHA256Transform(context, context->block);
  for (index = 0; index < 8; ++index) {
    digest[4 * index] = (unsigned char)(context->state[index] >> 24);
    digest[4 * index + 1] = (unsigned char)(context->state[index] >> 16);
    digest[4 * index + 2] = (unsigned char)(context->state[index] >> 8);
    digest[4 * index + 3] = (unsigned char)context->state[index];
  }
}

static void DigestToHex(const unsigned char digest[32], char output[65]) {
  static const char hex[] = "0123456789abcdef";
  int index;
  for (index = 0; index < 32; ++index) {
    output[2 * index] = hex[digest[index] >> 4];
    output[2 * index + 1] = hex[digest[index] & 15u];
  }
  output[64] = '\0';
}

static MVMCPowerLanczosObservableCensusStatus SetDiagnostic(
    MVMCPowerLanczosObservableCensusStatus status, char *diagnostic,
    size_t diagnostic_capacity, const char *format, ...) {
  if (diagnostic != NULL && diagnostic_capacity > 0) {
    va_list arguments;
    va_start(arguments, format);
    (void)vsnprintf(diagnostic, diagnostic_capacity, format, arguments);
    va_end(arguments);
  }
  return status;
}

static int IsWhitespaceOnly(const char *text) {
  const unsigned char *cursor = (const unsigned char *)text;
  while (*cursor != '\0') {
    if (!isspace(*cursor)) return 0;
    ++cursor;
  }
  return 1;
}

static int ReadBoundedLine(FILE *stream, char *line, size_t capacity) {
  size_t length;
  if (fgets(line, (int)capacity, stream) == NULL) return feof(stream) ? 1 : -1;
  length = strlen(line);
  if (length > 0 && line[length - 1] == '\n') return 0;
  if (feof(stream)) return 0;
  return -2;
}

static int ParseCountLine(const char *line, int *count) {
  const unsigned char *cursor = (const unsigned char *)line;
  char *end;
  long value;
  if (line == NULL || count == NULL) return 0;
  while (isspace(*cursor)) ++cursor;
  if (*cursor == '\0') return 0;
  while (*cursor != '\0' && !isspace(*cursor)) ++cursor;
  if (*cursor == '\0') return 0;
  while (isspace(*cursor)) ++cursor;
  errno = 0;
  value = strtol((const char *)cursor, &end, 10);
  if (errno == ERANGE || end == (const char *)cursor || value < INT_MIN ||
      value > INT_MAX) {
    return 0;
  }
  while (isspace((unsigned char)*end)) ++end;
  if (*end != '\0') return 0;
  *count = (int)value;
  return 1;
}

static int ParseIntegerRow(const char *line, int width, int indices[8]) {
  const char *cursor = line;
  int index;
  for (index = 0; index < width; ++index) {
    char *end;
    long value;
    while (isspace((unsigned char)*cursor)) ++cursor;
    if (*cursor == '\0') return 0;
    errno = 0;
    value = strtol(cursor, &end, 10);
    if (errno == ERANGE || end == cursor || value < INT_MIN ||
        value > INT_MAX) {
      return 0;
    }
    indices[index] = (int)value;
    cursor = end;
  }
  while (isspace((unsigned char)*cursor)) ++cursor;
  return *cursor == '\0';
}

static const char *PathBasename(const char *path) {
  const char *slash = strrchr(path, '/');
  const char *backslash = strrchr(path, '\\');
  const char *separator = slash;
  if (backslash != NULL && (separator == NULL || backslash > separator)) {
    separator = backslash;
  }
  return separator == NULL ? path : separator + 1;
}

static MVMCPowerLanczosObservableCensusStatus HashFile(
    const char *path, char output[65], int *contains_nul, char *diagnostic,
    size_t diagnostic_capacity) {
  unsigned char buffer[8192];
  unsigned char digest[32];
  SHA256Context context;
  FILE *stream = fopen(path, "rb");
  size_t count;
  if (stream == NULL) {
    return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_IO_ERROR,
                         diagnostic, diagnostic_capacity,
                         "cannot open observable file '%s'", path);
  }
  SHA256Init(&context);
  *contains_nul = 0;
  while ((count = fread(buffer, 1, sizeof(buffer), stream)) > 0) {
    if (memchr(buffer, '\0', count) != NULL) *contains_nul = 1;
    SHA256Update(&context, buffer, count);
  }
  if (ferror(stream)) {
    fclose(stream);
    return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_IO_ERROR,
                         diagnostic, diagnostic_capacity,
                         "cannot read observable file '%s'", path);
  }
  fclose(stream);
  SHA256Final(&context, digest);
  DigestToHex(digest, output);
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
}

static MVMCPowerLanczosObservableCensusStatus InspectFile(
    const char *path, MVMCPowerLanczosObservableFile *metadata,
    char *diagnostic, size_t diagnostic_capacity) {
  char line[MVMC_OBSERVABLE_INPUT_LINE_CAPACITY];
  char first_sha[65];
  const char *basename;
  FILE *stream;
  int header_index;
  int contains_nul = 0;
  MVMCPowerLanczosObservableCensusStatus status;
  memset(metadata, 0, sizeof(*metadata));
  if (path == NULL || path[0] == '\0') {
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
  }
  basename = PathBasename(path);
  if (basename[0] == '\0' ||
      strlen(basename) >= sizeof(metadata->raw_file_basename)) {
    return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT,
                         diagnostic, diagnostic_capacity,
                         "observable basename is empty or too long");
  }
  status = HashFile(path, first_sha, &contains_nul, diagnostic,
                    diagnostic_capacity);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  if (contains_nul) {
    return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ROW_WIDTH,
                         diagnostic, diagnostic_capacity,
                         "observable file '%s' contains a NUL byte", basename);
  }
  stream = fopen(path, "rb");
  if (stream == NULL) {
    return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_IO_ERROR,
                         diagnostic, diagnostic_capacity,
                         "cannot reopen observable file '%s'", path);
  }
  for (header_index = 0; header_index < 5; ++header_index) {
    int line_status = ReadBoundedLine(stream, line, sizeof(line));
    if (line_status != 0) {
      fclose(stream);
      return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_MALFORMED_COUNT,
                           diagnostic, diagnostic_capacity,
                           "observable file '%s' has an incomplete header",
                           basename);
    }
    if (header_index == 1 && !ParseCountLine(line, &metadata->signed_count)) {
      fclose(stream);
      return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_MALFORMED_COUNT,
                           diagnostic, diagnostic_capacity,
                           "observable file '%s' has a malformed signed count",
                           basename);
    }
  }
  fclose(stream);
  if (metadata->signed_count < 0) {
    return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_NEGATIVE_COUNT,
                         diagnostic, diagnostic_capacity,
                         "observable file '%s' has negative count %d", basename,
                         metadata->signed_count);
  }
  strcpy(metadata->raw_file_basename, basename);
  strcpy(metadata->raw_file_sha256, first_sha);
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
}

static int SameParsedRecord(const ParsedRecord *left,
                            const ParsedRecord *right) {
  return left->family == right->family && left->width == right->width &&
         memcmp(left->indices, right->indices,
                (size_t)left->width * sizeof(left->indices[0])) == 0;
}

static MVMCPowerLanczosObservableCensusStatus ParseRows(
    const char *path, MVMCPowerLanczosObservableFamily family, int count,
    int nsite, ParsedRecord *records, int record_offset, char *diagnostic,
    size_t diagnostic_capacity) {
  char line[MVMC_OBSERVABLE_INPUT_LINE_CAPACITY];
  FILE *stream;
  int header_index;
  int ordinal;
  int width = family == MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS ? 4 : 8;
  if (count == 0 && (path == NULL || path[0] == '\0')) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
  }
  stream = fopen(path, "rb");
  if (stream == NULL) {
    return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_IO_ERROR,
                         diagnostic, diagnostic_capacity,
                         "cannot reopen observable file '%s'", path);
  }
  for (header_index = 0; header_index < 5; ++header_index) {
    if (ReadBoundedLine(stream, line, sizeof(line)) != 0) {
      fclose(stream);
      return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_COUNT_MISMATCH,
                           diagnostic, diagnostic_capacity,
                           "observable file '%s' lost its header", path);
    }
  }
  for (ordinal = 0; ordinal < count; ++ordinal) {
    ParsedRecord *record = records + record_offset + ordinal;
    int index;
    int prior;
    int line_status = ReadBoundedLine(stream, line, sizeof(line));
    memset(record, 0, sizeof(*record));
    record->family = family;
    record->width = width;
    if (line_status == 1) {
      fclose(stream);
      return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_COUNT_MISMATCH,
                           diagnostic, diagnostic_capacity,
                           "%s row count is smaller than declared count %d",
                           mvmc_power_lanczos_observable_family_name(family),
                           count);
    }
    if (line_status != 0 ||
        !ParseIntegerRow(line, width, record->indices)) {
      fclose(stream);
      return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ROW_WIDTH,
                           diagnostic, diagnostic_capacity,
                           "%s raw ordinal %d does not contain exactly %d integers",
                           mvmc_power_lanczos_observable_family_name(family),
                           ordinal, width);
    }
    for (index = 0; index < width; index += 2) {
      if (record->indices[index] < 0 || record->indices[index] >= nsite) {
        fclose(stream);
        return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SITE_RANGE,
                             diagnostic, diagnostic_capacity,
                             "%s raw ordinal %d has site %d outside [0,%d)",
                             mvmc_power_lanczos_observable_family_name(family),
                             ordinal, record->indices[index], nsite);
      }
      if (record->indices[index + 1] < 0 ||
          record->indices[index + 1] > 1) {
        fclose(stream);
        return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SPIN_RANGE,
                             diagnostic, diagnostic_capacity,
                             "%s raw ordinal %d has spin %d outside {0,1}",
                             mvmc_power_lanczos_observable_family_name(family),
                             ordinal, record->indices[index + 1]);
      }
    }
    if (record->indices[1] != record->indices[3] ||
        (width == 8 && record->indices[5] != record->indices[7])) {
      fclose(stream);
      return SetDiagnostic(
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SPIN_NONCONSERVING,
          diagnostic, diagnostic_capacity,
          "%s raw ordinal %d is not bilinear-spin-conserving",
          mvmc_power_lanczos_observable_family_name(family), ordinal);
    }
    for (prior = record_offset; prior < record_offset + ordinal; ++prior) {
      if (SameParsedRecord(record, records + prior)) {
        fclose(stream);
        return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DUPLICATE,
                             diagnostic, diagnostic_capacity,
                             "%s raw ordinal %d duplicates ordinal %d",
                             mvmc_power_lanczos_observable_family_name(family),
                             ordinal, prior - record_offset);
      }
    }
  }
  while (1) {
    int line_status = ReadBoundedLine(stream, line, sizeof(line));
    if (line_status == 1) break;
    if (line_status != 0 || !IsWhitespaceOnly(line)) {
      fclose(stream);
      return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_COUNT_MISMATCH,
                           diagnostic, diagnostic_capacity,
                           "%s contains data after its declared %d rows",
                           mvmc_power_lanczos_observable_family_name(family),
                           count);
    }
  }
  fclose(stream);
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
}

static void MakeAdjoint(const ParsedRecord *source, ParsedRecord *adjoint) {
  memset(adjoint, 0, sizeof(*adjoint));
  adjoint->family = source->family;
  adjoint->width = source->width;
  if (source->width == 4) {
    adjoint->indices[0] = source->indices[2];
    adjoint->indices[1] = source->indices[3];
    adjoint->indices[2] = source->indices[0];
    adjoint->indices[3] = source->indices[1];
  } else {
    int pair;
    for (pair = 0; pair < 4; ++pair) {
      adjoint->indices[2 * pair] = source->indices[6 - 2 * pair];
      adjoint->indices[2 * pair + 1] = source->indices[7 - 2 * pair];
    }
  }
}

static MVMCPowerLanczosObservableAdjointClass ClassifyAdjoint(
    const ParsedRecord *records, int record_count, int record_index) {
  const ParsedRecord *source = records + record_index;
  ParsedRecord adjoint;
  int index;
  MakeAdjoint(source, &adjoint);
  if (SameParsedRecord(source, &adjoint) ||
      (source->width == 8 &&
       ((source->indices[0] == source->indices[2] &&
         source->indices[1] == source->indices[3] &&
         source->indices[4] == source->indices[6] &&
         source->indices[5] == source->indices[7]) ||
        (source->indices[4] == source->indices[2] &&
         source->indices[5] == source->indices[3] &&
         source->indices[6] == source->indices[0] &&
         source->indices[7] == source->indices[1])))) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_SELF_ADJOINT;
  }
  for (index = 0; index < record_count; ++index) {
    if (SameParsedRecord(records + index, &adjoint)) {
      return MVMC_POWER_LANCZOS_OBSERVABLE_ADJOINT_REQUESTED;
    }
  }
  return MVMC_POWER_LANCZOS_OBSERVABLE_ADJOINT_UNREQUESTED;
}

static int MakeOperatorId(const ParsedRecord *record, int nsite, char *output,
                          size_t output_capacity) {
  int orbitals[4] = {0, 0, 0, 0};
  int pair;
  if (record == NULL || output == NULL ||
      (record->width != 4 && record->width != 8)) {
    return -1;
  }
  for (pair = 0; pair < record->width / 2; ++pair) {
    orbitals[pair] =
        record->indices[2 * pair] + record->indices[2 * pair + 1] * nsite;
  }
  if (record->width == 4) {
    return snprintf(output, output_capacity, "cdag_%d_c_%d", orbitals[0],
                    orbitals[1]);
  }
  return snprintf(output, output_capacity, "cdag_%d_c_%d_cdag_%d_c_%d",
                  orbitals[0], orbitals[1], orbitals[2], orbitals[3]);
}

static int JsonReserve(JsonBuilder *builder, size_t addition) {
  size_t required;
  size_t capacity;
  char *replacement;
  if (addition > SIZE_MAX - builder->size - 1) return 0;
  required = builder->size + addition + 1;
  if (required <= builder->capacity) return 1;
  capacity = builder->capacity == 0 ? 1024 : builder->capacity;
  while (capacity < required) {
    if (capacity > SIZE_MAX / 2) {
      capacity = required;
      break;
    }
    capacity *= 2;
  }
  replacement = (char *)realloc(builder->bytes, capacity);
  if (replacement == NULL) return 0;
  builder->bytes = replacement;
  builder->capacity = capacity;
  return 1;
}

static int JsonAppendBytes(JsonBuilder *builder, const char *bytes,
                           size_t size) {
  if (!JsonReserve(builder, size)) return 0;
  memcpy(builder->bytes + builder->size, bytes, size);
  builder->size += size;
  builder->bytes[builder->size] = '\0';
  return 1;
}

static int JsonAppendLiteral(JsonBuilder *builder, const char *literal) {
  return JsonAppendBytes(builder, literal, strlen(literal));
}

static int JsonAppendInteger(JsonBuilder *builder, int value) {
  char buffer[64];
  int length = snprintf(buffer, sizeof(buffer), "%d", value);
  return length > 0 && (size_t)length < sizeof(buffer) &&
         JsonAppendBytes(builder, buffer, (size_t)length);
}

static int JsonAppendEscaped(JsonBuilder *builder, const char *value) {
  const unsigned char *cursor = (const unsigned char *)value;
  if (!JsonAppendLiteral(builder, "\"")) return 0;
  while (*cursor != '\0') {
    char escaped[7];
    if (*cursor == '"' || *cursor == '\\') {
      escaped[0] = '\\';
      escaped[1] = (char)*cursor;
      if (!JsonAppendBytes(builder, escaped, 2)) return 0;
    } else if (*cursor < 0x20u) {
      static const char hex[] = "0123456789abcdef";
      escaped[0] = '\\';
      escaped[1] = 'u';
      escaped[2] = '0';
      escaped[3] = '0';
      escaped[4] = hex[*cursor >> 4];
      escaped[5] = hex[*cursor & 15u];
      if (!JsonAppendBytes(builder, escaped, 6)) return 0;
    } else if (!JsonAppendBytes(builder, (const char *)cursor, 1)) {
      return 0;
    }
    ++cursor;
  }
  return JsonAppendLiteral(builder, "\"");
}

static int CheckedSizeAdd(size_t *total, size_t addition) {
  if (addition > SIZE_MAX - *total) return 0;
  *total += addition;
  return 1;
}

static int BoundedStringLength(const char *value, size_t capacity,
                               size_t *length) {
  size_t index;
  for (index = 0; index < capacity; ++index) {
    if (value[index] == '\0') {
      *length = index;
      return 1;
    }
  }
  return 0;
}

static void WriteU32(unsigned char **cursor, uint32_t value) {
  (*cursor)[0] = (unsigned char)(value >> 24);
  (*cursor)[1] = (unsigned char)(value >> 16);
  (*cursor)[2] = (unsigned char)(value >> 8);
  (*cursor)[3] = (unsigned char)value;
  *cursor += 4;
}

static int ReadU32(const unsigned char **cursor, const unsigned char *end,
                   uint32_t *value) {
  if ((size_t)(end - *cursor) < 4) return 0;
  *value = ((uint32_t)(*cursor)[0] << 24) |
           ((uint32_t)(*cursor)[1] << 16) |
           ((uint32_t)(*cursor)[2] << 8) | (uint32_t)(*cursor)[3];
  *cursor += 4;
  return 1;
}

static int SkipWireBytes(const unsigned char **cursor,
                         const unsigned char *end, size_t count) {
  if (count > (size_t)(end - *cursor)) return 0;
  *cursor += count;
  return 1;
}

static MVMCPowerLanczosObservableCensusStatus ValidatePlan(
    const MVMCPowerLanczosObservablePlan *plan) {
  ParsedRecord parsed[MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS];
  size_t census_digest_length;
  size_t semantic_digest_length;
  int family;
  int expected_count = 0;
  int record_index = 0;
  if (plan == NULL || plan->schema_version !=
                          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_VERSION ||
      plan->nsite < 1 || plan->nsite > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_NSITE ||
      plan->nsite_uc < 1 || plan->nsite_uc > plan->nsite ||
      plan->nsite_uc > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_NSITE_UC) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
  }
  if (!BoundedStringLength(plan->raw_observable_census_sha256,
                           sizeof(plan->raw_observable_census_sha256),
                           &census_digest_length) ||
      !BoundedStringLength(plan->semantic_observable_census_sha256,
                           sizeof(plan->semantic_observable_census_sha256),
                           &semantic_digest_length) ||
      (census_digest_length != 0 && census_digest_length != 64) ||
      (semantic_digest_length != 0 && semantic_digest_length != 64)) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
  }
  for (family = 0; family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT;
       ++family) {
    int count = plan->family_count[family];
    int ordinal;
    int formula_factor = family == MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS ? 2 : 6;
    int hard_cap = family == MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS
                       ? MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ONE_BODY
                       : MVMC_POWER_LANCZOS_OBSERVABLE_MAX_QUARTIC;
    int formula_cap = formula_factor * plan->nsite_uc * plan->nsite;
    size_t basename_length;
    size_t file_digest_length;
    size_t digest_index;
    if (count < 0 || count > formula_cap || count > hard_cap ||
        plan->files[family].signed_count != count) {
      return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_RESOURCE_LIMIT;
    }
    if (!BoundedStringLength(plan->files[family].raw_file_basename,
                             sizeof(plan->files[family].raw_file_basename),
                             &basename_length) ||
        !BoundedStringLength(plan->files[family].raw_file_sha256,
                             sizeof(plan->files[family].raw_file_sha256),
                             &file_digest_length) ||
        ((basename_length == 0) != (file_digest_length == 0)) ||
        (basename_length > 0 && file_digest_length != 64) ||
        (count > 0 && basename_length == 0)) {
      return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
    }
    if (file_digest_length == 64) {
      for (digest_index = 0; digest_index < 64; ++digest_index) {
        const char digit = plan->files[family].raw_file_sha256[digest_index];
        if (!((digit >= '0' && digit <= '9') ||
              (digit >= 'a' && digit <= 'f'))) {
          return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
        }
      }
    }
    if (count > INT_MAX - expected_count) {
      return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SIZE_OVERFLOW;
    }
    expected_count += count;
    for (ordinal = 0; ordinal < count; ++ordinal, ++record_index) {
      const MVMCPowerLanczosObservableRecord *record;
      ParsedRecord *parsed_record;
      char operator_id[MVMC_POWER_LANCZOS_OBSERVABLE_OPERATOR_ID_CAPACITY];
      size_t record_operator_id_length;
      int index;
      if (plan->records == NULL) {
        return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
      }
      record = plan->records + record_index;
      parsed_record = parsed + record_index;
      if (record->family != (MVMCPowerLanczosObservableFamily)family ||
          record->raw_ordinal != ordinal ||
          record->row_width != (family == 0 ? 4 : 8)) {
        return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
      }
      memset(parsed_record, 0, sizeof(*parsed_record));
      parsed_record->family = record->family;
      parsed_record->width = record->row_width;
      memcpy(parsed_record->indices, record->raw_indices,
             (size_t)parsed_record->width *
                 sizeof(parsed_record->indices[0]));
      for (index = parsed_record->width;
           index < MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ROW_WIDTH; ++index) {
        if (record->raw_indices[index] != 0) {
          return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
        }
      }
      for (index = 0; index < parsed_record->width; index += 2) {
        if (parsed_record->indices[index] < 0 ||
            parsed_record->indices[index] >= plan->nsite) {
          return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SITE_RANGE;
        }
        if (parsed_record->indices[index + 1] < 0 ||
            parsed_record->indices[index + 1] > 1) {
          return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SPIN_RANGE;
        }
      }
      if (parsed_record->indices[1] != parsed_record->indices[3] ||
          (parsed_record->width == 8 &&
           parsed_record->indices[5] != parsed_record->indices[7])) {
        return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SPIN_NONCONSERVING;
      }
      if (!BoundedStringLength(record->canonical_operator_id,
                               sizeof(record->canonical_operator_id),
                               &record_operator_id_length) ||
          MakeOperatorId(parsed_record, plan->nsite, operator_id,
                         sizeof(operator_id)) <= 0 ||
          strlen(operator_id) != record_operator_id_length ||
          memcmp(operator_id, record->canonical_operator_id,
                 record_operator_id_length) != 0) {
        return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
      }
    }
  }
  if (expected_count != plan->record_count ||
      expected_count > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_COUNT_MISMATCH;
  }
  for (record_index = 0; record_index < expected_count; ++record_index) {
    int prior;
    if (ClassifyAdjoint(parsed, expected_count, record_index) !=
        plan->records[record_index].adjoint_class) {
      return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
    }
    for (prior = 0; prior < record_index; ++prior) {
      if (SameParsedRecord(parsed + prior, parsed + record_index)) {
        return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DUPLICATE;
      }
    }
  }
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
}

static MVMCPowerLanczosObservableCensusStatus BuildCanonicalJson(
    const MVMCPowerLanczosObservablePlan *plan, JsonBuilder *builder) {
  MVMCPowerLanczosObservableCensusStatus status = ValidatePlan(plan);
  int record_index;
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  memset(builder, 0, sizeof(*builder));
  if (!JsonAppendLiteral(builder, "[")) goto allocation_failure;
  for (record_index = 0; record_index < plan->record_count; ++record_index) {
    const MVMCPowerLanczosObservableRecord *record =
        plan->records + record_index;
    const MVMCPowerLanczosObservableFile *file = plan->files + record->family;
    int index;
    if (record_index > 0 && !JsonAppendLiteral(builder, ",")) {
      goto allocation_failure;
    }
    if (!JsonAppendLiteral(builder, "{\"raw_file_basename\":" ) ||
        !JsonAppendEscaped(builder, file->raw_file_basename) ||
        !JsonAppendLiteral(builder, ",\"raw_file_sha256\":" ) ||
        !JsonAppendEscaped(builder, file->raw_file_sha256) ||
        !JsonAppendLiteral(builder, ",\"raw_ordinal\":" ) ||
        !JsonAppendInteger(builder, record->raw_ordinal) ||
        !JsonAppendLiteral(builder, ",\"family\":" ) ||
        !JsonAppendEscaped(
            builder, mvmc_power_lanczos_observable_family_name(record->family)) ||
        !JsonAppendLiteral(builder, ",\"raw_indices\":[")) {
      goto allocation_failure;
    }
    for (index = 0; index < record->row_width; ++index) {
      if ((index > 0 && !JsonAppendLiteral(builder, ",")) ||
          !JsonAppendInteger(builder, record->raw_indices[index])) {
        goto allocation_failure;
      }
    }
    if (!JsonAppendLiteral(builder, "],\"canonical_operator_id\":" ) ||
        !JsonAppendEscaped(builder, record->canonical_operator_id) ||
        !JsonAppendLiteral(builder, ",\"adjoint_class\":" ) ||
        !JsonAppendEscaped(
            builder,
            mvmc_power_lanczos_observable_adjoint_class_name(
                record->adjoint_class)) ||
        !JsonAppendLiteral(builder, "}")) {
      goto allocation_failure;
    }
  }
  if (!JsonAppendLiteral(builder, "]")) goto allocation_failure;
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;

allocation_failure:
  free(builder->bytes);
  memset(builder, 0, sizeof(*builder));
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ALLOCATION_FAILURE;
}

static MVMCPowerLanczosObservableCensusStatus ComputePlanHash(
    const MVMCPowerLanczosObservablePlan *plan, char output[65]) {
  JsonBuilder builder;
  SHA256Context context;
  unsigned char digest[32];
  MVMCPowerLanczosObservableCensusStatus status =
      BuildCanonicalJson(plan, &builder);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  SHA256Init(&context);
  SHA256Update(&context, builder.bytes, builder.size);
  SHA256Final(&context, digest);
  DigestToHex(digest, output);
  free(builder.bytes);
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
}

static MVMCPowerLanczosObservableCensusStatus BuildSemanticCanonicalJson(
    const MVMCPowerLanczosObservablePlan *plan, JsonBuilder *builder) {
  MVMCPowerLanczosObservableCensusStatus status = ValidatePlan(plan);
  int family;
  int record_index;
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  memset(builder, 0, sizeof(*builder));
  if (!JsonAppendLiteral(
          builder,
          "{\"schema\":\"raw_observable_semantic_v1\",\"nsite\":" ) ||
      !JsonAppendInteger(builder, plan->nsite) ||
      !JsonAppendLiteral(builder, ",\"nsite_uc\":" ) ||
      !JsonAppendInteger(builder, plan->nsite_uc) ||
      !JsonAppendLiteral(builder, ",\"family_count\":[")) {
    goto allocation_failure;
  }
  for (family = 0; family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT;
       ++family) {
    if ((family > 0 && !JsonAppendLiteral(builder, ",")) ||
        !JsonAppendInteger(builder, plan->family_count[family])) {
      goto allocation_failure;
    }
  }
  if (!JsonAppendLiteral(builder, "],\"records\":[")) {
    goto allocation_failure;
  }
  for (record_index = 0; record_index < plan->record_count; ++record_index) {
    const MVMCPowerLanczosObservableRecord *record =
        plan->records + record_index;
    int index;
    if ((record_index > 0 && !JsonAppendLiteral(builder, ",")) ||
        !JsonAppendLiteral(builder, "{\"raw_ordinal\":" ) ||
        !JsonAppendInteger(builder, record->raw_ordinal) ||
        !JsonAppendLiteral(builder, ",\"family\":" ) ||
        !JsonAppendEscaped(
            builder, mvmc_power_lanczos_observable_family_name(record->family)) ||
        !JsonAppendLiteral(builder, ",\"raw_indices\":[")) {
      goto allocation_failure;
    }
    for (index = 0; index < record->row_width; ++index) {
      if ((index > 0 && !JsonAppendLiteral(builder, ",")) ||
          !JsonAppendInteger(builder, record->raw_indices[index])) {
        goto allocation_failure;
      }
    }
    if (!JsonAppendLiteral(builder, "],\"canonical_operator_id\":" ) ||
        !JsonAppendEscaped(builder, record->canonical_operator_id) ||
        !JsonAppendLiteral(builder, ",\"adjoint_class\":" ) ||
        !JsonAppendEscaped(
            builder,
            mvmc_power_lanczos_observable_adjoint_class_name(
                record->adjoint_class)) ||
        !JsonAppendLiteral(builder, "}")) {
      goto allocation_failure;
    }
  }
  if (!JsonAppendLiteral(builder, "]}")) goto allocation_failure;
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;

allocation_failure:
  free(builder->bytes);
  memset(builder, 0, sizeof(*builder));
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ALLOCATION_FAILURE;
}

static MVMCPowerLanczosObservableCensusStatus ComputePlanSemanticHash(
    const MVMCPowerLanczosObservablePlan *plan, char output[65]) {
  JsonBuilder builder;
  SHA256Context context;
  unsigned char digest[32];
  MVMCPowerLanczosObservableCensusStatus status =
      BuildSemanticCanonicalJson(plan, &builder);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  SHA256Init(&context);
  SHA256Update(&context, builder.bytes, builder.size);
  SHA256Final(&context, digest);
  DigestToHex(digest, output);
  free(builder.bytes);
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
}

static MVMCPowerLanczosObservableCensusStatus
ComputeLegacyAugmentedPlanHash(
    const MVMCPowerLanczosLegacyAugmentedPlan *plan, char output[65]) {
  JsonBuilder builder;
  SHA256Context context;
  unsigned char digest[32];
  int ordinal;
  if (plan == NULL ||
      plan->schema_version != MVMC_POWER_LANCZOS_LEGACY_AUGMENTED_PLAN_VERSION ||
      plan->nsite <= 0 ||
      plan->nsite > INT_MAX / 2 || plan->pair_count < 0 ||
      (plan->pair_count > 0 && plan->orbital_pairs == NULL)) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
  }
  memset(&builder, 0, sizeof(builder));
  if (!JsonAppendLiteral(&builder,
                         "{\"schema\":\"legacy_augmented_one_body_v1\","
                         "\"nsite\":") ||
      !JsonAppendInteger(&builder, plan->nsite) ||
      !JsonAppendLiteral(&builder, ",\"pairs\":[")) {
    goto allocation_failure;
  }
  for (ordinal = 0; ordinal < plan->pair_count; ++ordinal) {
    const int creation = plan->orbital_pairs[2 * ordinal];
    const int annihilation = plan->orbital_pairs[2 * ordinal + 1];
    if (creation < 0 || creation >= 2 * plan->nsite || annihilation < 0 ||
        annihilation >= 2 * plan->nsite) {
      free(builder.bytes);
      return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
    }
    if ((ordinal > 0 && !JsonAppendLiteral(&builder, ",")) ||
        !JsonAppendLiteral(&builder, "[") ||
        !JsonAppendInteger(&builder, creation) ||
        !JsonAppendLiteral(&builder, ",") ||
        !JsonAppendInteger(&builder, annihilation) ||
        !JsonAppendLiteral(&builder, "]")) {
      goto allocation_failure;
    }
  }
  if (!JsonAppendLiteral(&builder, "]}")) goto allocation_failure;
  SHA256Init(&context);
  SHA256Update(&context, builder.bytes, builder.size);
  SHA256Final(&context, digest);
  DigestToHex(digest, output);
  free(builder.bytes);
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;

allocation_failure:
  free(builder.bytes);
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ALLOCATION_FAILURE;
}

void mvmc_power_lanczos_observable_plan_init(
    MVMCPowerLanczosObservablePlan *plan) {
  if (plan != NULL) memset(plan, 0, sizeof(*plan));
}

void mvmc_power_lanczos_observable_plan_destroy(
    MVMCPowerLanczosObservablePlan *plan) {
  if (plan == NULL) return;
  free(plan->records);
  memset(plan, 0, sizeof(*plan));
}

void mvmc_power_lanczos_legacy_augmented_plan_init(
    MVMCPowerLanczosLegacyAugmentedPlan *plan) {
  if (plan != NULL) memset(plan, 0, sizeof(*plan));
}

void mvmc_power_lanczos_legacy_augmented_plan_destroy(
    MVMCPowerLanczosLegacyAugmentedPlan *plan) {
  if (plan == NULL) return;
  free(plan->orbital_pairs);
  memset(plan, 0, sizeof(*plan));
}

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_legacy_augmented_plan_build(
    int nsite, int pair_count, int *const *pair_to_ordinal,
    MVMCPowerLanczosLegacyAugmentedPlan *plan, char *diagnostic,
    size_t diagnostic_capacity) {
  MVMCPowerLanczosLegacyAugmentedPlan candidate;
  unsigned char *seen = NULL;
  size_t pair_slots;
  int orbital_count;
  int creation;
  int observed = 0;
  MVMCPowerLanczosObservableCensusStatus status;
  if (diagnostic != NULL && diagnostic_capacity > 0) diagnostic[0] = '\0';
  if (plan == NULL || pair_to_ordinal == NULL || nsite <= 0 ||
      nsite > INT_MAX / 2 || pair_count < 0) {
    return SetDiagnostic(
        MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT, diagnostic,
        diagnostic_capacity, "invalid legacy augmented plan arguments");
  }
  orbital_count = 2 * nsite;
  if ((long long)pair_count >
          (long long)orbital_count * (long long)orbital_count ||
      (size_t)pair_count > SIZE_MAX / (2u * sizeof(int))) {
    return SetDiagnostic(
        MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SIZE_OVERFLOW, diagnostic,
        diagnostic_capacity, "legacy augmented plan size overflow");
  }
  mvmc_power_lanczos_legacy_augmented_plan_init(&candidate);
  candidate.schema_version =
      MVMC_POWER_LANCZOS_LEGACY_AUGMENTED_PLAN_VERSION;
  candidate.nsite = nsite;
  candidate.pair_count = pair_count;
  pair_slots = (size_t)pair_count * 2u;
  if (pair_count > 0) {
    candidate.orbital_pairs = (int *)calloc(pair_slots, sizeof(int));
    seen = (unsigned char *)calloc((size_t)pair_count, 1);
    if (candidate.orbital_pairs == NULL || seen == NULL) {
      free(seen);
      mvmc_power_lanczos_legacy_augmented_plan_destroy(&candidate);
      return SetDiagnostic(
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ALLOCATION_FAILURE,
          diagnostic, diagnostic_capacity,
          "cannot allocate legacy augmented plan");
    }
  }
  for (creation = 0; creation < orbital_count; ++creation) {
    int annihilation;
    if (pair_to_ordinal[creation] == NULL) {
      status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
      goto reject;
    }
    for (annihilation = 0; annihilation < orbital_count; ++annihilation) {
      const int ordinal = pair_to_ordinal[creation][annihilation];
      if (ordinal == -1) continue;
      if (ordinal < 0 || ordinal >= pair_count) {
        status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
        goto reject;
      }
      if (seen[ordinal]) {
        status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DUPLICATE;
        goto reject;
      }
      seen[ordinal] = 1;
      candidate.orbital_pairs[2 * ordinal] = creation;
      candidate.orbital_pairs[2 * ordinal + 1] = annihilation;
      ++observed;
    }
  }
  if (observed != pair_count) {
    status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_COUNT_MISMATCH;
    goto reject;
  }
  status = ComputeLegacyAugmentedPlanHash(
      &candidate, candidate.legacy_augmented_one_body_sha256);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) goto reject;
  free(seen);
  mvmc_power_lanczos_legacy_augmented_plan_destroy(plan);
  *plan = candidate;
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;

reject:
  free(seen);
  mvmc_power_lanczos_legacy_augmented_plan_destroy(&candidate);
  return SetDiagnostic(status, diagnostic, diagnostic_capacity,
                       "invalid legacy augmented one-body index map");
}

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_legacy_augmented_plan_rehash(
    const MVMCPowerLanczosLegacyAugmentedPlan *plan, char *actual_sha256,
    size_t actual_sha256_capacity) {
  char actual[65];
  MVMCPowerLanczosObservableCensusStatus status =
      ComputeLegacyAugmentedPlanHash(plan, actual);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  if (actual_sha256 != NULL) {
    if (actual_sha256_capacity < sizeof(actual)) {
      return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
    }
    memcpy(actual_sha256, actual, sizeof(actual));
  }
  if (strcmp(actual, plan->legacy_augmented_one_body_sha256) != 0) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DIGEST_MISMATCH;
  }
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
}

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_plan_build_from_files(
    int nsite, int nsite_uc,
    const char *const paths[MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT],
    MVMCPowerLanczosObservablePlan *plan, char *diagnostic,
    size_t diagnostic_capacity) {
  ParsedRecord parsed[MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS];
  MVMCPowerLanczosObservablePlan candidate;
  int family;
  int total = 0;
  int offset = 0;
  MVMCPowerLanczosObservableCensusStatus status;
  if (diagnostic != NULL && diagnostic_capacity > 0) diagnostic[0] = '\0';
  if (paths == NULL || plan == NULL || nsite < 1 ||
      nsite > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_NSITE || nsite_uc < 1 ||
      nsite_uc > nsite ||
      nsite_uc > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_NSITE_UC) {
    return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_RESOURCE_LIMIT,
                         diagnostic, diagnostic_capacity,
                         "P6 observable scope requires 1 <= NsiteUC <= Nsite <= 16");
  }
  mvmc_power_lanczos_observable_plan_init(&candidate);
  candidate.schema_version = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_VERSION;
  candidate.nsite = nsite;
  candidate.nsite_uc = nsite_uc;
  for (family = 0; family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT;
       ++family) {
    int count;
    int factor = family == MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS ? 2 : 6;
    int hard_cap = family == MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS
                       ? MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ONE_BODY
                       : MVMC_POWER_LANCZOS_OBSERVABLE_MAX_QUARTIC;
    int formula_cap = factor * nsite_uc * nsite;
    status = InspectFile(paths[family], candidate.files + family, diagnostic,
                         diagnostic_capacity);
    if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
    count = candidate.files[family].signed_count;
    if (count > formula_cap || count > hard_cap) {
      return SetDiagnostic(
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_RESOURCE_LIMIT, diagnostic,
          diagnostic_capacity,
          "%s count %d exceeds formula cap %d or hard cap %d",
          mvmc_power_lanczos_observable_family_name(
              (MVMCPowerLanczosObservableFamily)family),
          count, formula_cap, hard_cap);
    }
    if (count > INT_MAX - total) {
      return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SIZE_OVERFLOW,
                           diagnostic, diagnostic_capacity,
                           "observable family count sum overflows int");
    }
    candidate.family_count[family] = count;
    total += count;
  }
  if (total > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS) {
    return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_RESOURCE_LIMIT,
                         diagnostic, diagnostic_capacity,
                         "combined observable count %d exceeds hard cap %d",
                         total, MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS);
  }
  memset(parsed, 0, (size_t)total * sizeof(parsed[0]));
  for (family = 0; family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT;
       ++family) {
    status = ParseRows(paths[family],
                       (MVMCPowerLanczosObservableFamily)family,
                       candidate.family_count[family], nsite, parsed, offset,
                       diagnostic, diagnostic_capacity);
    if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
    {
      char current_sha[65];
      int contains_nul = 0;
      if (paths[family] != NULL && paths[family][0] != '\0') {
        status = HashFile(paths[family], current_sha, &contains_nul, diagnostic,
                          diagnostic_capacity);
        if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
        if (contains_nul || strcmp(current_sha,
                                   candidate.files[family].raw_file_sha256) != 0) {
          return SetDiagnostic(
              MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SOURCE_CHANGED, diagnostic,
              diagnostic_capacity, "%s changed during strict census parsing",
              candidate.files[family].raw_file_basename);
        }
      }
    }
    offset += candidate.family_count[family];
  }
  candidate.record_count = total;
  if (total > 0) {
    candidate.records = (MVMCPowerLanczosObservableRecord *)calloc(
        (size_t)total, sizeof(*candidate.records));
    if (candidate.records == NULL) {
      return SetDiagnostic(
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ALLOCATION_FAILURE, diagnostic,
          diagnostic_capacity, "cannot allocate validated observable plan");
    }
  }
  for (offset = 0; offset < total; ++offset) {
    MVMCPowerLanczosObservableRecord *record = candidate.records + offset;
    int family_start = 0;
    int id_length;
    for (family = 0; family < (int)parsed[offset].family; ++family) {
      family_start += candidate.family_count[family];
    }
    record->raw_ordinal = offset - family_start;
    record->family = parsed[offset].family;
    record->row_width = parsed[offset].width;
    memcpy(record->raw_indices, parsed[offset].indices,
           (size_t)record->row_width * sizeof(record->raw_indices[0]));
    id_length = MakeOperatorId(parsed + offset, nsite,
                               record->canonical_operator_id,
                               sizeof(record->canonical_operator_id));
    if (id_length <= 0 ||
        (size_t)id_length >= sizeof(record->canonical_operator_id)) {
      mvmc_power_lanczos_observable_plan_destroy(&candidate);
      return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SIZE_OVERFLOW,
                           diagnostic, diagnostic_capacity,
                           "canonical operator ID exceeds its fixed capacity");
    }
    record->adjoint_class = ClassifyAdjoint(parsed, total, offset);
  }
  status = ComputePlanSemanticHash(
      &candidate, candidate.semantic_observable_census_sha256);
  if (status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) {
    status = ComputePlanHash(&candidate,
                             candidate.raw_observable_census_sha256);
  }
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) {
    mvmc_power_lanczos_observable_plan_destroy(&candidate);
    return SetDiagnostic(status, diagnostic, diagnostic_capacity,
                         "cannot hash canonical observable census");
  }
  mvmc_power_lanczos_observable_plan_destroy(plan);
  *plan = candidate;
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
}

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_plan_rehash(
    const MVMCPowerLanczosObservablePlan *plan, char *actual_sha256,
    size_t actual_sha256_capacity) {
  char actual[65];
  char semantic_actual[65];
  MVMCPowerLanczosObservableCensusStatus status = ComputePlanHash(plan, actual);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  status = ComputePlanSemanticHash(plan, semantic_actual);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  if (actual_sha256 != NULL) {
    if (actual_sha256_capacity < sizeof(actual)) {
      return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
    }
    strcpy(actual_sha256, actual);
  }
  if (strcmp(plan->raw_observable_census_sha256, actual) != 0 ||
      strcmp(plan->semantic_observable_census_sha256, semantic_actual) != 0) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DIGEST_MISMATCH;
  }
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
}

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_plan_semantic_rehash(
    const MVMCPowerLanczosObservablePlan *plan, char *actual_sha256,
    size_t actual_sha256_capacity) {
  char actual[65];
  MVMCPowerLanczosObservableCensusStatus status =
      ComputePlanSemanticHash(plan, actual);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  if (actual_sha256 != NULL) {
    if (actual_sha256_capacity < sizeof(actual)) {
      return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
    }
    memcpy(actual_sha256, actual, sizeof(actual));
  }
  return strcmp(plan->semantic_observable_census_sha256, actual) == 0
             ? MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK
             : MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DIGEST_MISMATCH;
}

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_plan_write_canonical_json(
    const MVMCPowerLanczosObservablePlan *plan, FILE *stream) {
  JsonBuilder builder;
  MVMCPowerLanczosObservableCensusStatus status;
  if (stream == NULL) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
  }
  status = BuildCanonicalJson(plan, &builder);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  if (builder.size > 0 &&
      fwrite(builder.bytes, 1, builder.size, stream) != builder.size) {
    free(builder.bytes);
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_IO_ERROR;
  }
  free(builder.bytes);
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
}

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_plan_wire_size(
    const MVMCPowerLanczosObservablePlan *plan, size_t *wire_size) {
  size_t total = sizeof(ObservableWireMagic) + 7u * sizeof(uint32_t) + 128u;
  char rehash[65];
  int family;
  int record_index;
  MVMCPowerLanczosObservableCensusStatus status;
  if (wire_size == NULL) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
  }
  status = mvmc_power_lanczos_observable_plan_rehash(plan, rehash,
                                                     sizeof(rehash));
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  for (family = 0; family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT;
       ++family) {
    size_t basename_length;
    size_t digest_length;
    if (!BoundedStringLength(plan->files[family].raw_file_basename,
                             sizeof(plan->files[family].raw_file_basename),
                             &basename_length) ||
        !BoundedStringLength(plan->files[family].raw_file_sha256,
                             sizeof(plan->files[family].raw_file_sha256),
                             &digest_length) ||
        !CheckedSizeAdd(&total, 3u * sizeof(uint32_t)) ||
        !CheckedSizeAdd(&total, basename_length) ||
        !CheckedSizeAdd(&total, digest_length)) {
      return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SIZE_OVERFLOW;
    }
  }
  for (record_index = 0; record_index < plan->record_count; ++record_index) {
    size_t operator_id_length;
    if (!BoundedStringLength(
            plan->records[record_index].canonical_operator_id,
            sizeof(plan->records[record_index].canonical_operator_id),
            &operator_id_length) ||
        !CheckedSizeAdd(&total, 13u * sizeof(uint32_t)) ||
        !CheckedSizeAdd(&total, operator_id_length)) {
      return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SIZE_OVERFLOW;
    }
  }
  *wire_size = total;
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
}

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_plan_pack(
    const MVMCPowerLanczosObservablePlan *plan, void *wire,
    size_t wire_capacity, size_t *wire_size) {
  unsigned char *cursor = (unsigned char *)wire;
  size_t required;
  int family;
  int record_index;
  MVMCPowerLanczosObservableCensusStatus status =
      mvmc_power_lanczos_observable_plan_wire_size(plan, &required);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) return status;
  if (wire_size != NULL) *wire_size = required;
  if (wire == NULL || wire_capacity < required) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT;
  }
  memcpy(cursor, ObservableWireMagic, sizeof(ObservableWireMagic));
  cursor += sizeof(ObservableWireMagic);
  WriteU32(&cursor, (uint32_t)plan->schema_version);
  WriteU32(&cursor, (uint32_t)plan->nsite);
  WriteU32(&cursor, (uint32_t)plan->nsite_uc);
  WriteU32(&cursor, (uint32_t)plan->record_count);
  for (family = 0; family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT;
       ++family) {
    WriteU32(&cursor, (uint32_t)plan->family_count[family]);
  }
  for (family = 0; family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT;
       ++family) {
    size_t basename_length = strlen(plan->files[family].raw_file_basename);
    size_t digest_length = strlen(plan->files[family].raw_file_sha256);
    WriteU32(&cursor, (uint32_t)plan->files[family].signed_count);
    WriteU32(&cursor, (uint32_t)basename_length);
    memcpy(cursor, plan->files[family].raw_file_basename, basename_length);
    cursor += basename_length;
    WriteU32(&cursor, (uint32_t)digest_length);
    memcpy(cursor, plan->files[family].raw_file_sha256, digest_length);
    cursor += digest_length;
  }
  for (record_index = 0; record_index < plan->record_count; ++record_index) {
    const MVMCPowerLanczosObservableRecord *record =
        plan->records + record_index;
    size_t operator_id_length = strlen(record->canonical_operator_id);
    int index;
    WriteU32(&cursor, (uint32_t)record->raw_ordinal);
    WriteU32(&cursor, (uint32_t)record->family);
    WriteU32(&cursor, (uint32_t)record->row_width);
    for (index = 0; index < MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ROW_WIDTH;
         ++index) {
      WriteU32(&cursor, (uint32_t)record->raw_indices[index]);
    }
    WriteU32(&cursor, (uint32_t)operator_id_length);
    memcpy(cursor, record->canonical_operator_id, operator_id_length);
    cursor += operator_id_length;
    WriteU32(&cursor, (uint32_t)record->adjoint_class);
  }
  memcpy(cursor, plan->raw_observable_census_sha256, 64);
  cursor += 64;
  memcpy(cursor, plan->semantic_observable_census_sha256, 64);
  cursor += 64;
  if ((size_t)(cursor - (unsigned char *)wire) != required) {
    return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SIZE_OVERFLOW;
  }
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
}

MVMCPowerLanczosObservableCensusStatus
mvmc_power_lanczos_observable_plan_unpack(
    const void *wire, size_t wire_size, MVMCPowerLanczosObservablePlan *plan,
    char *diagnostic, size_t diagnostic_capacity) {
  const unsigned char *begin = (const unsigned char *)wire;
  const unsigned char *cursor = begin;
  const unsigned char *end;
  const unsigned char *record_begin;
  const unsigned char *digest_begin;
  MVMCPowerLanczosObservablePlan candidate;
  ParsedRecord parsed[MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS];
  MVMCPowerLanczosObservableAdjointClass
      adjoint_classes[MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS];
  uint32_t value;
  int family;
  int record_index;
  int count_sum = 0;
  MVMCPowerLanczosObservableCensusStatus status;
  if (diagnostic != NULL && diagnostic_capacity > 0) diagnostic[0] = '\0';
  if (wire == NULL || plan == NULL || wire_size < sizeof(ObservableWireMagic) ||
      wire_size > (size_t)INT_MAX) {
    return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT,
                         diagnostic, diagnostic_capacity,
                         "invalid observable census wire arguments");
  }
  end = begin + wire_size;
  if (memcmp(cursor, ObservableWireMagic, sizeof(ObservableWireMagic)) != 0) {
    return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT,
                         diagnostic, diagnostic_capacity,
                         "observable census wire magic mismatch");
  }
  cursor += sizeof(ObservableWireMagic);
  mvmc_power_lanczos_observable_plan_init(&candidate);
  if (!ReadU32(&cursor, end, &value)) goto truncated;
  candidate.schema_version = (int)value;
  if (!ReadU32(&cursor, end, &value) || value > INT_MAX) goto truncated;
  candidate.nsite = (int)value;
  if (!ReadU32(&cursor, end, &value) || value > INT_MAX) goto truncated;
  candidate.nsite_uc = (int)value;
  if (!ReadU32(&cursor, end, &value) || value > INT_MAX) goto truncated;
  candidate.record_count = (int)value;
  for (family = 0; family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT;
       ++family) {
    if (!ReadU32(&cursor, end, &value) || value > INT_MAX) goto truncated;
    candidate.family_count[family] = (int)value;
    if (candidate.family_count[family] > INT_MAX - count_sum) {
      return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SIZE_OVERFLOW,
                           diagnostic, diagnostic_capacity,
                           "observable wire family sum overflows int");
    }
    count_sum += candidate.family_count[family];
  }
  if (candidate.schema_version !=
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_VERSION ||
      candidate.nsite < 1 ||
      candidate.nsite > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_NSITE ||
      candidate.nsite_uc < 1 || candidate.nsite_uc > candidate.nsite ||
      candidate.record_count != count_sum ||
      candidate.record_count > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS) {
    return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_RESOURCE_LIMIT,
                         diagnostic, diagnostic_capacity,
                         "observable wire count or scope metadata is invalid");
  }
  for (family = 0; family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT;
       ++family) {
    int factor = family == MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS ? 2 : 6;
    int hard_cap = family == MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS
                       ? MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ONE_BODY
                       : MVMC_POWER_LANCZOS_OBSERVABLE_MAX_QUARTIC;
    int formula_cap = factor * candidate.nsite_uc * candidate.nsite;
    if (candidate.family_count[family] > formula_cap ||
        candidate.family_count[family] > hard_cap) {
      return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_RESOURCE_LIMIT,
                           diagnostic, diagnostic_capacity,
                           "observable wire family %d exceeds its count cap",
                           family);
    }
  }
  for (family = 0; family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT;
       ++family) {
    uint32_t basename_length;
    uint32_t digest_length;
    const unsigned char *basename_bytes;
    const unsigned char *digest_bytes;
    uint32_t digest_index;
    if (!ReadU32(&cursor, end, &value) || value > INT_MAX ||
        !ReadU32(&cursor, end, &basename_length) ||
        basename_length >= sizeof(candidate.files[family].raw_file_basename)) {
      goto truncated;
    }
    basename_bytes = cursor;
    if (!SkipWireBytes(&cursor, end, basename_length) ||
        !ReadU32(&cursor, end, &digest_length) ||
        digest_length >= sizeof(candidate.files[family].raw_file_sha256)) {
      goto truncated;
    }
    digest_bytes = cursor;
    if (!SkipWireBytes(&cursor, end, digest_length)) goto truncated;
    candidate.files[family].signed_count = (int)value;
    if (candidate.files[family].signed_count !=
            candidate.family_count[family] ||
        ((basename_length == 0) != (digest_length == 0)) ||
        (basename_length > 0 && digest_length != 64) ||
        (candidate.family_count[family] > 0 && basename_length == 0)) {
      return SetDiagnostic(
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT, diagnostic,
          diagnostic_capacity,
          "observable wire family %d metadata does not match its count",
          family);
    }
    if (memchr(basename_bytes, '\0', basename_length) != NULL) {
      return SetDiagnostic(
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT, diagnostic,
          diagnostic_capacity,
          "observable wire family %d basename contains a NUL byte", family);
    }
    for (digest_index = 0; digest_index < digest_length; ++digest_index) {
      const unsigned char digit = digest_bytes[digest_index];
      if (!((digit >= '0' && digit <= '9') ||
            (digit >= 'a' && digit <= 'f'))) {
        return SetDiagnostic(
            MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT,
            diagnostic, diagnostic_capacity,
            "observable wire family %d SHA-256 is not lowercase hex", family);
      }
    }
  }
  record_begin = cursor;
  memset(parsed, 0, (size_t)candidate.record_count * sizeof(parsed[0]));
  for (record_index = 0; record_index < candidate.record_count;
       ++record_index) {
    ParsedRecord *parsed_record = parsed + record_index;
    uint32_t ordinal;
    uint32_t record_family;
    uint32_t row_width;
    uint32_t operator_id_length;
    uint32_t adjoint_class;
    const unsigned char *operator_id;
    char expected_operator_id[
        MVMC_POWER_LANCZOS_OBSERVABLE_OPERATOR_ID_CAPACITY];
    int expected_family = 0;
    int expected_family_start = 0;
    int index;
    if (!ReadU32(&cursor, end, &ordinal) ||
        !ReadU32(&cursor, end, &record_family) ||
        !ReadU32(&cursor, end, &row_width)) {
      goto truncated;
    }
    while (expected_family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT &&
           record_index >=
               expected_family_start +
                   candidate.family_count[expected_family]) {
      expected_family_start += candidate.family_count[expected_family];
      ++expected_family;
    }
    if (expected_family >= MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT ||
        record_family != (uint32_t)expected_family ||
        ordinal != (uint32_t)(record_index - expected_family_start) ||
        row_width != (uint32_t)(expected_family == 0 ? 4 : 8)) {
      return SetDiagnostic(
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT, diagnostic,
          diagnostic_capacity,
          "observable wire record %d breaks family/ordinal ordering",
          record_index);
    }
    parsed_record->family = (MVMCPowerLanczosObservableFamily)record_family;
    parsed_record->width = (int)row_width;
    for (index = 0; index < MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ROW_WIDTH;
         ++index) {
      if (!ReadU32(&cursor, end, &value)) goto truncated;
      if (value > INT_MAX) {
        return SetDiagnostic(
            MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT,
            diagnostic, diagnostic_capacity,
            "observable wire record %d has a negative/oversized index",
            record_index);
      }
      parsed_record->indices[index] = (int)value;
      if (index >= parsed_record->width && value != 0) {
        return SetDiagnostic(
            MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT,
            diagnostic, diagnostic_capacity,
            "observable wire record %d has nonzero index padding",
            record_index);
      }
    }
    if (!ReadU32(&cursor, end, &operator_id_length) ||
        operator_id_length == 0 ||
        operator_id_length >=
            MVMC_POWER_LANCZOS_OBSERVABLE_OPERATOR_ID_CAPACITY ||
        operator_id_length > (uint32_t)(end - cursor)) {
      goto truncated;
    }
    operator_id = cursor;
    if (!SkipWireBytes(&cursor, end, operator_id_length) ||
        !ReadU32(&cursor, end, &adjoint_class)) {
      goto truncated;
    }
    if (ordinal > INT_MAX ||
        record_family >= MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT ||
        (row_width != 4 && row_width != 8) ||
        adjoint_class > MVMC_POWER_LANCZOS_OBSERVABLE_ADJOINT_UNREQUESTED) {
      return SetDiagnostic(
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT, diagnostic,
          diagnostic_capacity, "observable wire record %d is invalid",
          record_index);
    }
    for (index = 0; index < parsed_record->width; index += 2) {
      if (parsed_record->indices[index] < 0 ||
          parsed_record->indices[index] >= candidate.nsite) {
        return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SITE_RANGE,
                             diagnostic, diagnostic_capacity,
                             "observable wire record %d has invalid site",
                             record_index);
      }
      if (parsed_record->indices[index + 1] < 0 ||
          parsed_record->indices[index + 1] > 1) {
        return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SPIN_RANGE,
                             diagnostic, diagnostic_capacity,
                             "observable wire record %d has invalid spin",
                             record_index);
      }
    }
    if (parsed_record->indices[1] != parsed_record->indices[3] ||
        (parsed_record->width == 8 &&
         parsed_record->indices[5] != parsed_record->indices[7])) {
      return SetDiagnostic(
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SPIN_NONCONSERVING,
          diagnostic, diagnostic_capacity,
          "observable wire record %d is not spin-conserving", record_index);
    }
    if (MakeOperatorId(parsed_record, candidate.nsite, expected_operator_id,
                       sizeof(expected_operator_id)) <= 0 ||
        strlen(expected_operator_id) != operator_id_length ||
        memcmp(expected_operator_id, operator_id, operator_id_length) != 0) {
      return SetDiagnostic(
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT, diagnostic,
          diagnostic_capacity,
          "observable wire record %d has a stale operator ID", record_index);
    }
    adjoint_classes[record_index] =
        (MVMCPowerLanczosObservableAdjointClass)adjoint_class;
  }
  digest_begin = cursor;
  if ((size_t)(end - cursor) != 128) goto truncated;
  memcpy(candidate.raw_observable_census_sha256, cursor, 64);
  candidate.raw_observable_census_sha256[64] = '\0';
  for (record_index = 0; record_index < 64; ++record_index) {
    const char digit = candidate.raw_observable_census_sha256[record_index];
    if (!((digit >= '0' && digit <= '9') ||
          (digit >= 'a' && digit <= 'f'))) {
      return SetDiagnostic(
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT, diagnostic,
          diagnostic_capacity, "observable wire digest is not lowercase hex");
    }
  }
  cursor += 64;
  memcpy(candidate.semantic_observable_census_sha256, cursor, 64);
  candidate.semantic_observable_census_sha256[64] = '\0';
  for (record_index = 0; record_index < 64; ++record_index) {
    const char digit =
        candidate.semantic_observable_census_sha256[record_index];
    if (!((digit >= '0' && digit <= '9') ||
          (digit >= 'a' && digit <= 'f'))) {
      return SetDiagnostic(
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT, diagnostic,
          diagnostic_capacity,
          "observable wire semantic digest is not lowercase hex");
    }
  }
  cursor += 64;
  if (cursor != end) goto truncated;
  for (record_index = 0; record_index < candidate.record_count;
       ++record_index) {
    int prior;
    if (ClassifyAdjoint(parsed, candidate.record_count, record_index) !=
        adjoint_classes[record_index]) {
      return SetDiagnostic(
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT, diagnostic,
          diagnostic_capacity,
          "observable wire record %d has a stale adjoint class", record_index);
    }
    for (prior = 0; prior < record_index; ++prior) {
      if (SameParsedRecord(parsed + prior, parsed + record_index)) {
        return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DUPLICATE,
                             diagnostic, diagnostic_capacity,
                             "observable wire record %d is duplicated",
                             record_index);
      }
    }
  }

  cursor = begin + sizeof(ObservableWireMagic) + 7u * sizeof(uint32_t);
  for (family = 0; family < MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT;
       ++family) {
    uint32_t basename_length;
    uint32_t digest_length;
    const unsigned char *basename_bytes;
    const unsigned char *digest_bytes;
    if (!ReadU32(&cursor, end, &value)) goto truncated;
    candidate.files[family].signed_count = (int)value;
    if (!ReadU32(&cursor, end, &basename_length) ||
        basename_length >=
            sizeof(candidate.files[family].raw_file_basename)) {
      goto truncated;
    }
    basename_bytes = cursor;
    if (!SkipWireBytes(&cursor, end, basename_length)) goto truncated;
    memcpy(candidate.files[family].raw_file_basename, basename_bytes,
           basename_length);
    candidate.files[family].raw_file_basename[basename_length] = '\0';
    if (!ReadU32(&cursor, end, &digest_length) ||
        digest_length >= sizeof(candidate.files[family].raw_file_sha256)) {
      goto truncated;
    }
    digest_bytes = cursor;
    if (!SkipWireBytes(&cursor, end, digest_length)) goto truncated;
    memcpy(candidate.files[family].raw_file_sha256, digest_bytes,
           digest_length);
    candidate.files[family].raw_file_sha256[digest_length] = '\0';
  }
  if (cursor != record_begin) goto truncated;
  if (candidate.record_count > 0) {
    candidate.records = (MVMCPowerLanczosObservableRecord *)calloc(
        (size_t)candidate.record_count, sizeof(*candidate.records));
    if (candidate.records == NULL) {
      return SetDiagnostic(
          MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ALLOCATION_FAILURE, diagnostic,
          diagnostic_capacity, "cannot allocate received observable plan");
    }
  }
  for (record_index = 0; record_index < candidate.record_count;
       ++record_index) {
    MVMCPowerLanczosObservableRecord *record =
        candidate.records + record_index;
    uint32_t operator_id_length;
    const unsigned char *operator_id_bytes;
    int index;
    if (!ReadU32(&cursor, end, &value)) goto truncated;
    record->raw_ordinal = (int)value;
    if (!ReadU32(&cursor, end, &value)) goto truncated;
    record->family = (MVMCPowerLanczosObservableFamily)value;
    if (!ReadU32(&cursor, end, &value)) goto truncated;
    record->row_width = (int)value;
    for (index = 0; index < MVMC_POWER_LANCZOS_OBSERVABLE_MAX_ROW_WIDTH;
         ++index) {
      if (!ReadU32(&cursor, end, &value)) goto truncated;
      record->raw_indices[index] = (int)(uint32_t)value;
    }
    if (!ReadU32(&cursor, end, &operator_id_length) ||
        operator_id_length >= sizeof(record->canonical_operator_id)) {
      goto truncated;
    }
    operator_id_bytes = cursor;
    if (!SkipWireBytes(&cursor, end, operator_id_length)) goto truncated;
    memcpy(record->canonical_operator_id, operator_id_bytes,
           operator_id_length);
    record->canonical_operator_id[operator_id_length] = '\0';
    if (!ReadU32(&cursor, end, &value)) goto truncated;
    record->adjoint_class = (MVMCPowerLanczosObservableAdjointClass)value;
  }
  if (cursor != digest_begin) goto truncated;
  status = mvmc_power_lanczos_observable_plan_rehash(&candidate, NULL, 0);
  if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) {
    mvmc_power_lanczos_observable_plan_destroy(&candidate);
    return SetDiagnostic(status, diagnostic, diagnostic_capacity,
                         "received observable plan failed canonical rehash");
  }
  mvmc_power_lanczos_observable_plan_destroy(plan);
  *plan = candidate;
  return MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;

truncated:
  mvmc_power_lanczos_observable_plan_destroy(&candidate);
  return SetDiagnostic(MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT,
                       diagnostic, diagnostic_capacity,
                       "observable census wire is truncated or oversized");
}

const char *mvmc_power_lanczos_observable_family_name(
    MVMCPowerLanczosObservableFamily family) {
  switch (family) {
    case MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS:
      return "cisajs";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX:
      return "cisajscktaltex";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALT:
      return "cisajscktalt";
  }
  return "invalid";
}

const char *mvmc_power_lanczos_observable_adjoint_class_name(
    MVMCPowerLanczosObservableAdjointClass adjoint_class) {
  switch (adjoint_class) {
    case MVMC_POWER_LANCZOS_OBSERVABLE_SELF_ADJOINT:
      return "self_adjoint";
    case MVMC_POWER_LANCZOS_OBSERVABLE_ADJOINT_REQUESTED:
      return "has_requested_adjoint";
    case MVMC_POWER_LANCZOS_OBSERVABLE_ADJOINT_UNREQUESTED:
      return "nonhermitian_unpaired";
  }
  return "invalid";
}

const char *mvmc_power_lanczos_observable_census_status_string(
    MVMCPowerLanczosObservableCensusStatus status) {
  switch (status) {
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK:
      return "ok";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_INVALID_ARGUMENT:
      return "invalid argument";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_IO_ERROR:
      return "I/O error";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_MALFORMED_COUNT:
      return "malformed signed count";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_NEGATIVE_COUNT:
      return "negative signed count";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_RESOURCE_LIMIT:
      return "P6 resource limit";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_COUNT_MISMATCH:
      return "raw row count mismatch";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ROW_WIDTH:
      return "raw row width";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SITE_RANGE:
      return "site range";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SPIN_RANGE:
      return "spin range";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SPIN_NONCONSERVING:
      return "bilinear spin nonconserving";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DUPLICATE:
      return "duplicate ordered raw row";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SIZE_OVERFLOW:
      return "size overflow";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_ALLOCATION_FAILURE:
      return "allocation failure";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SOURCE_CHANGED:
      return "raw source changed while parsing";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_DIGEST_MISMATCH:
      return "canonical census digest mismatch";
  }
  return "unknown observable census status";
}
