#include "power_lanczos_sha256.h"

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include <stdint.h>
#include <string.h>

typedef struct {
  uint32_t state[8];
  uint64_t bit_count;
  unsigned char block[64];
  size_t block_size;
} SHA256Context;

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

static uint32_t rotate_right(uint32_t value, unsigned int shift) {
  return (value >> shift) | (value << (32u - shift));
}

static void transform(SHA256Context *context,
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
    const uint32_t s0 = rotate_right(words[index - 15], 7) ^
                        rotate_right(words[index - 15], 18) ^
                        (words[index - 15] >> 3);
    const uint32_t s1 = rotate_right(words[index - 2], 17) ^
                        rotate_right(words[index - 2], 19) ^
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
    const uint32_t sum0 = rotate_right(a, 2) ^ rotate_right(a, 13) ^
                          rotate_right(a, 22);
    const uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
    const uint32_t sum1 = rotate_right(e, 6) ^ rotate_right(e, 11) ^
                          rotate_right(e, 25);
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

static void initialize(SHA256Context *context) {
  static const uint32_t initial[8] = {
      0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
      0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u};
  memcpy(context->state, initial, sizeof(initial));
  context->bit_count = 0;
  context->block_size = 0;
}

static void update(SHA256Context *context, const unsigned char *bytes,
                   size_t input_size) {
  while (input_size > 0) {
    const size_t available = sizeof(context->block) - context->block_size;
    const size_t take = input_size < available ? input_size : available;
    memcpy(context->block + context->block_size, bytes, take);
    context->block_size += take;
    context->bit_count += (uint64_t)take * UINT64_C(8);
    bytes += take;
    input_size -= take;
    if (context->block_size == sizeof(context->block)) {
      transform(context, context->block);
      context->block_size = 0;
    }
  }
}

static void finalize(SHA256Context *context, unsigned char digest[32]) {
  const uint64_t bit_count = context->bit_count;
  int index;
  context->block[context->block_size++] = 0x80u;
  if (context->block_size > 56) {
    while (context->block_size < 64) {
      context->block[context->block_size++] = 0;
    }
    transform(context, context->block);
    context->block_size = 0;
  }
  while (context->block_size < 56) {
    context->block[context->block_size++] = 0;
  }
  for (index = 7; index >= 0; --index) {
    context->block[context->block_size++] =
        (unsigned char)(bit_count >> (8 * index));
  }
  transform(context, context->block);
  for (index = 0; index < 8; ++index) {
    digest[4 * index] = (unsigned char)(context->state[index] >> 24);
    digest[4 * index + 1] =
        (unsigned char)(context->state[index] >> 16);
    digest[4 * index + 2] =
        (unsigned char)(context->state[index] >> 8);
    digest[4 * index + 3] = (unsigned char)context->state[index];
  }
  memset(context, 0, sizeof(*context));
}

int mvmc_power_lanczos_sha256_bytes(
    const void *input, size_t input_size, unsigned char digest[32]) {
  SHA256Context context;
  if (digest == NULL || (input == NULL && input_size != 0) ||
      input_size > UINT64_MAX / UINT64_C(8)) {
    return 0;
  }
  initialize(&context);
  if (input_size != 0) {
    update(&context, (const unsigned char *)input, input_size);
  }
  finalize(&context, digest);
  return 1;
}

int mvmc_power_lanczos_sha256_hex(const void *input, size_t input_size,
                                  char hex[65]) {
  static const char digits[] = "0123456789abcdef";
  unsigned char digest[32];
  size_t index;
  if (hex == NULL ||
      !mvmc_power_lanczos_sha256_bytes(input, input_size, digest)) {
    return 0;
  }
  for (index = 0; index < sizeof(digest); ++index) {
    hex[2 * index] = digits[digest[index] >> 4];
    hex[2 * index + 1] = digits[digest[index] & 15u];
  }
  hex[64] = '\0';
  memset(digest, 0, sizeof(digest));
  return 1;
}

#endif
