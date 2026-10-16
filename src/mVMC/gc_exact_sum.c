#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>

#include "gc_exact_sum.h"

/* Each finite binary64 is an integer multiple of 2^-1074. Its magnitude
 * occupies at most 2098 bits; summing at most INT_MAX values needs at most
 * another sizeof(int)*CHAR_BIT bits. Separate positive/negative base-2^32
 * expansions avoid intermediate overflow even for DBL_MAX cancellation.
 * Integer-only accumulation also makes fast-math reassociation irrelevant.
 * This helper is compiled separately from the production amalgamation. */
enum { GC_EXACT_LIMBS = (2098 + sizeof(int) * CHAR_BIT + 31) / 32 };

static int GCExactAdd(uint32_t *sum, unsigned int word, uint64_t value) {
  while (value != 0) {
    if (word >= GC_EXACT_LIMBS) return -1;
    value += sum[word];
    sum[word++] = (uint32_t)value;
    value >>= 32;
  }
  return 0;
}

int GCExactComplexSumIsZero(const double complex *values, const int count,
                          int *isZero) {
  int component, i, allZero = 1;
  if (isZero == NULL) return -1;
  *isZero = 0;
  if (values == NULL || count <= 0) return -1;
  if (FLT_RADIX != 2 || DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024 ||
      DBL_MIN_EXP != -1021 || sizeof(double) != sizeof(uint64_t) ||
      sizeof(double complex) != 2 * sizeof(uint64_t)) return -1;
  for (component = 0; component < 2; ++component) {
    uint32_t sums[2][GC_EXACT_LIMBS] = {{0}};
    for (i = 0; i < count; ++i) {
      uint64_t bits[2], significand;
      unsigned int exponent, shift, word, offset, sign;
      /* C complex has the representation of {real, imaginary} doubles. */
      memcpy(bits, &values[i], sizeof bits);
      exponent = (unsigned int)((bits[component] >> 52) & 0x7ffU);
      if (exponent == 0x7ffU) return -1;
      sign = (unsigned int)(bits[component] >> 63);
      significand = bits[component] & UINT64_C(0x000fffffffffffff);
      shift = exponent == 0 ? 0 : exponent - 1;
      if (exponent != 0) significand |= UINT64_C(0x0010000000000000);
      word = shift / 32;
      offset = shift % 32;
      if (GCExactAdd(sums[sign], word,
                     (uint64_t)(uint32_t)significand << offset) != 0 ||
          GCExactAdd(sums[sign], word + 1,
                     (significand >> 32) << offset) != 0) return -1;
    }
    if (memcmp(sums[0], sums[1], sizeof sums[0]) != 0) allZero = 0;
  }
  *isZero = allZero;
  return 0;
}

int GCExactComplexSumSelfTest(void) {
  const double one = 1.0;
  uint64_t bits;
  const double complex zero[] = {0.25, -0.25 * I, -0.25, 0.25 * I};
  const double complex nonzero[] = {1.0, 0x1p-53, -1.0};
  const double complex imaginary[] = {I, 0x1p-53 * I, -I};
  const uint64_t invalidBits[2] = {UINT64_C(0x7ff0000000000000), 0};
  double complex invalid;
  int isZero;
  if (sizeof one != sizeof bits || sizeof invalid != sizeof invalidBits)
    return -1;
  memcpy(&bits, &one, sizeof bits);
  if (bits != UINT64_C(0x3ff0000000000000)) return -1;
  memcpy(&invalid, invalidBits, sizeof invalid);
  if (GCExactComplexSumIsZero(zero, 4, &isZero) != 0 || !isZero ||
      GCExactComplexSumIsZero(nonzero, 3, &isZero) != 0 || isZero ||
      GCExactComplexSumIsZero(imaginary, 3, &isZero) != 0 || isZero ||
      GCExactComplexSumIsZero(&invalid, 1, &isZero) == 0) return -1;
  return 0;
}
