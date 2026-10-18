/* Independent support, proposal distribution and detailed-balance checks. */
#include "gc_antiparallel.h"
#include "gc_spin_support.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures;
#define CHECK(x) do { if (!(x)) { \
  fprintf(stderr, "GC spin support failure at line %d\n", __LINE__); \
  ++failures; } } while (0)

static uint32_t random_state = 1729;
static uint32_t random_word(void) {
  random_state ^= random_state << 13;
  random_state ^= random_state >> 17;
  random_state ^= random_state << 5;
  return random_state;
}

static int bits(unsigned x) {
  int n = 0;
  for (; x; x >>= 1) n += (int)(x & 1U);
  return n;
}

static double choose2(int n) { return (double)n * (n - 1) / 2.0; }

static void symmetry(void) {
  int idx[2][2] = {{0, 1}, {2, 3}}, signs[2][2] = {{1, 1}, {1, 1}};
  int *ir[2] = {idx[0], idx[1]}, *sr[2] = {signs[0], signs[1]};
  double complex f[4] = {0, 0.5 + 0.25 * I, -0.5 - 0.25 * I, 0};
  int k;
  CHECK(GCSpinPureTriplet(2, 4, ir, sr, f) == 1);
  f[3] = 1e-100; /* No tolerance-based rounding to the triplet manifold. */
  CHECK(GCSpinPureTriplet(2, 4, ir, sr, f) == 0);
  f[3] = 0;
  CHECK(GCSpinPureTriplet(2, 4, ir, sr, f) == 1);
  idx[1][0] = 1; signs[1][0] = -1;
  CHECK(GCSpinPureTriplet(2, 4, ir, sr, f) == 1);
  signs[1][0] = 1;
  CHECK(GCSpinPureTriplet(2, 4, ir, sr, f) == 0);
  f[3] = NAN;
  CHECK(GCSpinPureTriplet(2, 4, ir, sr, f) == -1);
  f[3] = 0; idx[1][1] = -1;
  CHECK(GCSpinPureTriplet(2, 4, ir, sr, f) == -1);
  for (k = 0; k <= 12; k += 2) {
    CHECK(GCSpinForbiddenSector(1, 0, 2, 0, 1, k) == (k % 4 == 2));
    CHECK(!GCSpinForbiddenSector(1, 0, 1, 0, 1, k));
    CHECK(!GCSpinForbiddenSector(0, 0, 2, 0, 1, k));
    CHECK(!GCSpinForbiddenSector(1, 1, 2, 0, 1, k));
    CHECK(!GCSpinForbiddenSector(1, 0, 2, 0, 0, k));
    CHECK(!GCSpinForbiddenSector(1, 0, 2, 0, -1, k));
  }
}

static void balance(void) {
  int l, m;
  for (l = 2; l <= 16; ++l) {
    for (m = 0; m + 2 <= l; ++m) {
      const double forward = 0.125 / pow(choose2(l - m), 2);
      const double reverse = 0.125 / pow(choose2(m + 2), 2);
      const double ratio = GCAntiRatioTwoPairs(m, l, 1);
      const double back = GCAntiRatioTwoPairs(m + 2, l, 0);
      /* A nontrivial pair of target weights, independent of proposal code. */
      const double wx = exp(-0.31 * m), wy = exp(-0.31 * (m + 2));
      CHECK(fabs(ratio - reverse / forward) < 1e-10 * ratio);
      CHECK(fabs(ratio * back - 1) < 1e-14);
      CHECK(fabs(wx * forward * fmin(1, wy / wx * ratio) -
                 wy * reverse * fmin(1, wx / wy * back)) < 1e-14);
    }
    CHECK(GCAntiRatioTwoPairs(l, l, 1) == 0);
    CHECK(GCAntiRatioTwoPairs(0, l, 0) == 0);
  }
  CHECK(GCAntiSelectBridge(0) == GC_ANTI_GLOBAL);
  CHECK(GCAntiSelectBridge(1.0 / 64) == GC_ANTI_ADD_TWO);
  CHECK(GCAntiSelectBridge(9.0 / 64) == GC_ANTI_REMOVE_TWO);
  CHECK(GCAntiSelectBridge(17.0 / 64) == -1);
}

static void draws(void) {
  int histogram[256] = {0}, before[8] = {0}, after[8];
  int iteration, mask, supported = 0;
  double ratio;
  for (iteration = 0; iteration < 140000; ++iteration) {
    int i;
    CHECK(GCAntiDrawBridge(GC_ANTI_GLOBAL, before, 4, 0, after,
                           random_word, &ratio) == 1);
    CHECK(ratio == 1);
    for (i = 0, mask = 0; i < 8; ++i) mask |= after[i] << i;
    CHECK(bits((unsigned)mask & 15) == bits((unsigned)mask >> 4));
    ++histogram[mask];
  }
  for (mask = 0; mask < 256; ++mask) {
    if (bits((unsigned)mask & 15) == bits((unsigned)mask >> 4)) {
      ++supported;
      CHECK(abs(histogram[mask] - 2000) < 230);
    } else CHECK(histogram[mask] == 0);
  }
  CHECK(supported == 70);
  memset(histogram, 0, sizeof histogram);
  for (iteration = 0; iteration < 72000; ++iteration) {
    int i, restored[8];
    CHECK(GCAntiDrawBridge(GC_ANTI_ADD_TWO, before, 4, 0, after,
                           random_word, &ratio) == 1);
    CHECK(ratio == 36);
    for (i = 0, mask = 0; i < 8; ++i) mask |= after[i] << i;
    CHECK(bits((unsigned)mask & 15) == 2 && bits((unsigned)mask >> 4) == 2);
    ++histogram[mask];
    CHECK(GCAntiDrawBridge(GC_ANTI_REMOVE_TWO, after, 4, 2, restored,
                           random_word, &ratio) == 1);
    CHECK(memcmp(restored, before, sizeof before) == 0);
  }
  for (mask = 0; mask < 256; ++mask)
    if (bits((unsigned)mask & 15) == 2 && bits((unsigned)mask >> 4) == 2)
      CHECK(abs(histogram[mask] - 2000) < 230);
}

int main(void) {
  symmetry(); balance(); draws();
  if (failures) return EXIT_FAILURE;
  puts("GC spin support and bridge proposals: PASS");
  return EXIT_SUCCESS;
}
