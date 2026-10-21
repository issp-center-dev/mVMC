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

/* ---- Finite total spin S>0 pure helpers (Task 1) ---- */

static void finite_symmetry(void) {
  int idx[3][3], signs[3][3], i, j;
  int *ir[3] = {idx[0], idx[1], idx[2]}, *sr[3] = {signs[0], signs[1], signs[2]};
  double complex f[9];
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j) { idx[i][j] = 3 * i + j; signs[i][j] = 1; }
  for (i = 0; i < 9; ++i) f[i] = 0;
  /* The zero table is both; S>0 treats it as a pure singlet. */
  CHECK(GCSpinPureSinglet(3, 9, ir, sr, f) == 1);
  CHECK(GCSpinPureTriplet(3, 9, ir, sr, f) == 1);
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j) f[3 * i + j] = (i + 1) * (j + 1) + 0.25 * I;
  CHECK(GCSpinPureSinglet(3, 9, ir, sr, f) == 1);
  CHECK(GCSpinPureTriplet(3, 9, ir, sr, f) == 0);
  f[5] = nextafter(creal(f[5]), 0.0) + cimag(f[5]) * I; /* F12 != F21 */
  CHECK(GCSpinPureSinglet(3, 9, ir, sr, f) == 0);
  f[5] = f[7];
  CHECK(GCSpinPureSinglet(3, 9, ir, sr, f) == 1);
  /* Shared parameter with sign -1 encodes an antisymmetric pair. */
  idx[2][1] = 5; signs[2][1] = -1;
  CHECK(GCSpinPureSinglet(3, 9, ir, sr, f) == 0);
  signs[2][1] = 1;
  CHECK(GCSpinPureSinglet(3, 9, ir, sr, f) == 1);
  f[0] = INFINITY;
  CHECK(GCSpinPureSinglet(3, 9, ir, sr, f) == -1);
  f[0] = 1; idx[0][0] = 9;
  CHECK(GCSpinPureSinglet(3, 9, ir, sr, f) == -1);
  idx[0][0] = 0; signs[0][0] = 0;
  CHECK(GCSpinPureSinglet(3, 9, ir, sr, f) == -1);
  /* Pair parity: m-S even; S=0 reduces to the Phase 3a N%4==2 rule. */
  for (i = 0; i <= 24; ++i) {
    int spin;
    CHECK(GCSpinForbiddenSector(1, 0, 2, 0, 1, i) ==
          (i % 2 == 0 && i % 4 == 2));
    for (spin = 0; spin <= 4; ++spin) {
      CHECK(GCSpinForbiddenSector(1, 0, 3, spin, 1, i) ==
            (i % 2 == 0 && (((i / 2 - spin) % 2) + 2) % 2 == 1));
      CHECK(!GCSpinForbiddenSector(1, 0, 1, spin, 1, i));
      CHECK(!GCSpinForbiddenSector(1, 0, 3, spin, 0, i));
      CHECK(!GCSpinForbiddenSector(1, 1, 3, spin, 1, i));
    }
  }
}

/* Brute force over every M=0 occupation of L<=6 sites. */
static void finite_support(void) {
  int nsite;
  for (nsite = 1; nsite <= 6; ++nsite) {
    const unsigned limit = 1U << (2 * nsite);
    int spin;
    for (spin = 0; spin <= nsite; ++spin) {
      int electrons;
      for (electrons = 0; electrons <= 2 * nsite; electrons += 2) {
        int seenLow = 1000, seenHigh = -1, low, high, triplet;
        unsigned mask;
        for (mask = 0; mask < limit; ++mask) {
          int occupation[12], rs, singles, doubles = 0;
          const unsigned up = mask & ((1U << nsite) - 1U);
          const unsigned down = mask >> nsite;
          if (bits(mask) != electrons || bits(up) != bits(down)) continue;
          for (rs = 0; rs < 2 * nsite; ++rs)
            occupation[rs] = (int)((mask >> rs) & 1U);
          singles = GCSpinCountSingles(occupation, nsite);
          CHECK(singles == bits(up ^ down));
          CHECK(GCSpinConfigurationZero(spin, singles) ==
                (spin > 0 && bits(up ^ down) < 2 * spin));
          doubles = bits(up & down);
          if (!GCSpinConfigurationZero(spin, singles)) {
            if (doubles < seenLow) seenLow = doubles;
            if (doubles > seenHigh) seenHigh = doubles;
          }
        }
        /* The closed form covers exactly the observed doublon counts. */
        CHECK(GCSpinDoublonRange(nsite, spin, electrons, &low, &high) ==
              (seenHigh >= 0));
        if (seenHigh >= 0) CHECK(low == seenLow && high == seenHigh);
        for (triplet = -1; triplet <= 1; ++triplet)
          CHECK(GCSpinSectorPossible(nsite, spin, triplet, electrons) ==
                (seenHigh >= 0 &&
                 !(triplet == 1 && (electrons / 2 - spin) % 2 != 0)));
      }
      CHECK(GCSpinDoublonRange(nsite, spin, 1, NULL, NULL) == -1);
      CHECK(!GCSpinSectorPossible(nsite, spin, 0, 2 * nsite + 2));
    }
  }
  {
    int bad[4] = {1, 2, 0, 1};
    CHECK(GCSpinCountSingles(bad, 2) == -1);
    CHECK(GCSpinCountSingles(NULL, 2) == -1);
  }
}

/* Structured initial draws: always valid, D uniform in its range, and the
 * configuration uniform for a given D (L=4: S=1 at N=4, S=2 at N=4). */
static void finite_initial(void) {
  int spin;
  for (spin = 1; spin <= 2; ++spin) {
    long histogram[256] = {0};
    long byDoubles[3] = {0, 0, 0};
    int perDoubles[3] = {0, 0, 0};
    int occupation[8], iteration, mask;
    const int draws = 120000;
    for (iteration = 0; iteration < draws; ++iteration) {
      int rs;
      CHECK(GCAntiDrawSpinConfiguration(4, 2, spin, random_word,
                                        occupation) == 1);
      for (rs = 0, mask = 0; rs < 8; ++rs) mask |= occupation[rs] << rs;
      CHECK(bits((unsigned)mask & 15) == 2 && bits((unsigned)mask >> 4) == 2);
      CHECK(bits(((unsigned)mask & 15) ^ ((unsigned)mask >> 4)) >= 2 * spin);
      ++histogram[mask];
    }
    for (mask = 0; mask < 256; ++mask) {
      const unsigned up = (unsigned)mask & 15, down = (unsigned)mask >> 4;
      if (bits(up) == 2 && bits(down) == 2 && bits(up ^ down) >= 2 * spin) {
        perDoubles[bits(up & down)]++;
        byDoubles[bits(up & down)] += histogram[mask];
      } else {
        CHECK(histogram[mask] == 0);
      }
    }
    if (spin == 1) {
      /* D=0 (6 states) and D=1 (24 states), each D with probability 1/2. */
      CHECK(perDoubles[0] == 6 && perDoubles[1] == 24 && perDoubles[2] == 0);
      CHECK(labs(byDoubles[0] - draws / 2) < 1000);
      for (mask = 0; mask < 256; ++mask) {
        const unsigned up = (unsigned)mask & 15, down = (unsigned)mask >> 4;
        if (bits(up) != 2 || bits(down) != 2 || bits(up ^ down) < 2) continue;
        if (bits(up & down) == 0)
          CHECK(labs(histogram[mask] - draws / 12) < 600);
        else
          CHECK(labs(histogram[mask] - draws / 48) < 330);
      }
    } else {
      CHECK(perDoubles[0] == 6 && byDoubles[0] == draws);
      for (mask = 0; mask < 256; ++mask)
        if (histogram[mask])
          CHECK(labs(histogram[mask] - draws / 6) < 900);
    }
  }
  {
    int occupation[8];
    CHECK(GCAntiDrawSpinConfiguration(4, 1, 2, random_word, occupation) == 0);
    CHECK(GCAntiDrawSpinConfiguration(4, 4, 1, random_word, occupation) == 0);
    CHECK(GCAntiDrawSpinConfiguration(4, 5, 0, random_word, occupation) == -1);
    CHECK(GCAntiDrawSpinConfiguration(4, 2, -1, random_word, occupation) ==
          -1);
    CHECK(GCAntiDrawSpinConfiguration(4, 2, 1, NULL, occupation) == -1);
  }
}

/* Spin swap: uniform over (up-only, down-only) pairs, counts unchanged, the
 * reverse proposal has the same probability, impossible swaps self-loop. */
static void finite_swap(void) {
  unsigned mask;
  for (mask = 0; mask < 256; ++mask) {
    const unsigned up = mask & 15, down = mask >> 4;
    const int upOnly = bits(up & ~down), downOnly = bits(down & ~up);
    int current[8], candidate[8], rs, iteration;
    long histogram[256] = {0};
    double ratio = 0;
    if (bits(up) != bits(down)) continue;
    for (rs = 0; rs < 8; ++rs) current[rs] = (int)((mask >> rs) & 1U);
    if (upOnly == 0 || downOnly == 0) {
      CHECK(GCAntiDrawBridge(GC_ANTI_SWAP, current, 4, bits(up), candidate,
                             random_word, &ratio) == 0);
      continue;
    }
    for (iteration = 0; iteration < 1600 * upOnly * downOnly; ++iteration) {
      unsigned target = 0;
      CHECK(GCAntiDrawBridge(GC_ANTI_SWAP, current, 4, bits(up), candidate,
                             random_word, &ratio) == 1);
      CHECK(ratio == 1.0);
      for (rs = 0; rs < 8; ++rs) target |= (unsigned)candidate[rs] << rs;
      ++histogram[target];
    }
    for (rs = 0; rs < 256; ++rs) {
      const unsigned t = (unsigned)rs, tu = t & 15, td = t >> 4;
      const unsigned changed = (up ^ tu) | (down ^ td);
      /* Valid swap: one up-only site i becomes down-only and one down-only
       * site j becomes up-only; doubles and empties are untouched. */
      const int valid = bits(up ^ tu) == 2 && (up ^ tu) == (down ^ td) &&
                        bits(up & ~down & ~tu & td) == 1 &&
                        bits(down & ~up & ~td & tu) == 1 &&
                        (changed & (up & down)) == 0;
      if (valid) {
        CHECK(labs(histogram[rs] - 1600) < 200);
        CHECK(bits(tu & ~td) == upOnly && bits(td & ~tu) == downOnly);
      } else {
        CHECK(histogram[rs] == 0);
      }
    }
  }
  /* Finite-S mixture: 1/64 global, 8/64 add, 8/64 remove, 8/64 swap. */
  {
    int k, counts[5] = {0, 0, 0, 0, 0};
    for (k = 0; k < 64; ++k) {
      const int move = GCAntiSelectFiniteSpinBridge((k + 0.5) / 64.0);
      counts[move < 0 ? 4 : move]++;
    }
    CHECK(counts[GC_ANTI_GLOBAL] == 1 && counts[GC_ANTI_ADD_TWO] == 8 &&
          counts[GC_ANTI_REMOVE_TWO] == 8 && counts[GC_ANTI_SWAP] == 8 &&
          counts[4] == 39);
  }
}

int main(void) {
  symmetry(); balance(); draws();
  finite_symmetry(); finite_support(); finite_initial(); finite_swap();
  if (failures) return EXIT_FAILURE;
  puts("GC spin support and bridge proposals: PASS");
  return EXIT_SUCCESS;
}
