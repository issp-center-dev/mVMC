#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gc_antiparallel.h"

#define MAX_SITES 4
#define MAX_ORBITALS (2 * MAX_SITES)
#define MAX_STATES (1 << MAX_ORBITALS)

static int failures = 0;

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "GC_AntiPair_Proposal_Unit FAIL: ");                     \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, "\n");                                                   \
      failures++;                                                               \
    }                                                                           \
  } while (0)

static int finite_close(const double actual, const double expected,
                        const double tolerance) {
  return isfinite(actual) && isfinite(expected) && isfinite(tolerance) &&
         tolerance >= 0.0 && fabs(actual - expected) <= tolerance;
}

/* Fake 32-bit generator returning a fixed script. */
static const uint32_t *fakeScript = NULL;
static int fakeLength = 0;
static int fakePosition = 0;

static uint32_t fake_next32(void) {
  if (fakePosition >= fakeLength) {
    fprintf(stderr, "GC_AntiPair_Proposal_Unit FAIL: fake RNG exhausted\n");
    exit(EXIT_FAILURE);
  }
  return fakeScript[fakePosition++];
}

static void set_script(const uint32_t *script, const int length) {
  fakeScript = script;
  fakeLength = length;
  fakePosition = 0;
}

static void test_mode_and_ratio(void) {
  CHECK(GCAntiEnabled(1, 0) && !GCAntiEnabled(0, 0) && !GCAntiEnabled(1, 1),
        "mode predicate");
  CHECK(GCAntiEnabled(2, 0) && !GCAntiEnabled(0, 1), "nonzero GC flag");
  CHECK(fabs(GCAntiRatioAdd(0, 4, 1.0, .25) - 4.0) < 1e-14, "add m=0");
  CHECK(fabs(GCAntiRatioRemove(1, 4, .25, 1.0) - .25) < 1e-14, "remove m=1");
  CHECK(fabs(GCAntiRatioAdd(3, 4, .25, 1.0) - .25) < 1e-14, "add m=L-1");
  CHECK(fabs(GCAntiRatioRemove(4, 4, 1.0, .25) - 4.0) < 1e-14, "remove m=L");
  CHECK(fabs(GCAntiRatioAdd(1, 3, .25, .25) - 1.0) < 1e-14, "add m=1 L=3");
  /* Invalid sector, site count or class probability returns 0. */
  CHECK(GCAntiRatioAdd(-1, 4, .25, .25) == 0.0, "add negative m");
  CHECK(GCAntiRatioAdd(4, 4, .25, .25) == 0.0, "add from full filling");
  CHECK(GCAntiRatioAdd(0, 0, .25, .25) == 0.0, "add nsite=0");
  CHECK(GCAntiRatioAdd(0, -4, .25, .25) == 0.0, "add negative nsite");
  CHECK(GCAntiRatioAdd(0, 4, 0.0, .25) == 0.0, "add pAdd=0");
  CHECK(GCAntiRatioAdd(0, 4, -.25, .25) == 0.0, "add pAdd<0");
  CHECK(GCAntiRatioAdd(0, 4, .25, -.25) == 0.0, "add pRemove<0");
  CHECK(GCAntiRatioAdd(0, 4, NAN, .25) == 0.0, "add pAdd NaN");
  CHECK(GCAntiRatioAdd(0, 4, .25, INFINITY) == 0.0, "add pRemove Inf");
  CHECK(GCAntiRatioRemove(0, 4, .25, .25) == 0.0, "remove from vacuum");
  CHECK(GCAntiRatioRemove(5, 4, .25, .25) == 0.0, "remove m>L");
  CHECK(GCAntiRatioRemove(-1, 4, .25, .25) == 0.0, "remove negative m");
  CHECK(GCAntiRatioRemove(1, 0, .25, .25) == 0.0, "remove nsite=0");
  CHECK(GCAntiRatioRemove(1, 4, 0.0, .25) == 0.0, "remove pRemove=0");
  CHECK(GCAntiRatioRemove(1, 4, .25, -1.0) == 0.0, "remove pAdd<0");
  CHECK(GCAntiRatioRemove(1, 4, NAN, .25) == 0.0, "remove pRemove NaN");
  CHECK(GCAntiRatioRemove(1, 4, .25, -INFINITY) == 0.0, "remove pAdd -Inf");
}

static void test_accept_log(void) {
  int accepted = -1;
  CHECK(GCAntiAcceptLog(0, 500, 0, 1, .9, &accepted) == 0 && accepted,
        "large ratio must not overflow into rejection");
  accepted = -1;
  CHECK(GCAntiAcceptLog(0, -INFINITY + 0.0 * I, 0, 1, 0, &accepted) == 0 &&
            !accepted,
        "zero candidate is a valid rejection");
  accepted = -1;
  CHECK(GCAntiAcceptLog(0, NAN + 0.0 * I, 0, 1, .5, &accepted) != 0 &&
            !accepted,
        "NaN candidate is an error");
  CHECK(GCAntiAcceptLog(0, -500, 0, 1, 0, &accepted) == 0 && accepted,
        "draw=0 accepts any finite ratio");
  CHECK(GCAntiAcceptLog(0, -500, 0, 1, 1e-300, &accepted) == 0 && !accepted,
        "tiny ratio is rejected for a positive draw");
  CHECK(GCAntiAcceptLog(0, 0, 0, 1, .999999, &accepted) == 0 && accepted,
        "ratio 1 accepts every draw below 1");
  /* exp(2*(0.1+0.2))*0.5 = 0.91105940...; probe both sides. */
  CHECK(GCAntiAcceptLog(0.3, 0.5, 0.1, 0.5, 0.911, &accepted) == 0 && accepted,
        "draw below the Metropolis weight");
  CHECK(GCAntiAcceptLog(0.3, 0.5, 0.1, 0.5, 0.9111, &accepted) == 0 &&
            !accepted,
        "draw above the Metropolis weight");
  CHECK(GCAntiAcceptLog(0, 1, 0, 1, .5, NULL) != 0, "NULL result");
  CHECK(GCAntiAcceptLog(NAN, 0, 0, 1, .5, &accepted) != 0, "old log NaN");
  CHECK(GCAntiAcceptLog(INFINITY, 0, 0, 1, .5, &accepted) != 0,
        "old log +Inf");
  CHECK(GCAntiAcceptLog(-INFINITY, 0, 0, 1, .5, &accepted) != 0,
        "old log -Inf");
  CHECK(GCAntiAcceptLog(0.0 + NAN * I, 0, 0, 1, .5, &accepted) != 0,
        "old log imaginary NaN");
  CHECK(GCAntiAcceptLog(0, INFINITY, 0, 1, .5, &accepted) != 0 && !accepted,
        "new log +Inf");
  CHECK(GCAntiAcceptLog(0, 0.0 + INFINITY * I, 0, 1, .5, &accepted) != 0,
        "new log imaginary Inf");
  CHECK(GCAntiAcceptLog(0, 0, NAN, 1, .5, &accepted) != 0,
        "projection log NaN");
  CHECK(GCAntiAcceptLog(0, 0, INFINITY, 1, .5, &accepted) != 0,
        "projection log Inf");
  CHECK(GCAntiAcceptLog(0, 0, 0, 0, .5, &accepted) != 0, "ratio=0");
  CHECK(GCAntiAcceptLog(0, 0, 0, -1, .5, &accepted) != 0, "ratio<0");
  CHECK(GCAntiAcceptLog(0, 0, 0, NAN, .5, &accepted) != 0, "ratio NaN");
  CHECK(GCAntiAcceptLog(0, 0, 0, INFINITY, .5, &accepted) != 0, "ratio Inf");
  CHECK(GCAntiAcceptLog(0, 0, 0, 1, -0.1, &accepted) != 0, "draw<0");
  CHECK(GCAntiAcceptLog(0, 0, 0, 1, 1.0, &accepted) != 0, "draw=1");
  CHECK(GCAntiAcceptLog(0, 0, 0, 1, NAN, &accepted) != 0, "draw NaN");
  CHECK(GCAntiAcceptLog(0, 0, 0, 1, 1e300, &accepted) != 0, "draw huge");
}

static void test_draw_below(void) {
  static const uint32_t script[] = {0, 1, 2, 3};
  static const uint32_t high[] = {UINT32_MAX, 7};
  static const uint32_t zero[] = {0};
  uint32_t value;
  set_script(script, 4);
  value = GCAntiDrawBelow(3, fake_next32);
  /* 2^32 mod 3 = 1, so the first word 0 is discarded. */
  CHECK(value == 1 && fakePosition == 2, "bound 3 rejection (value=%u)", value);
  set_script(zero, 1);
  CHECK(GCAntiDrawBelow(1, fake_next32) == 0 && fakePosition == 1, "bound 1");
  set_script(high, 2);
  CHECK(GCAntiDrawBelow(3, fake_next32) == UINT32_MAX % 3, "largest word");
  set_script(script, 4);
  CHECK(GCAntiDrawBelow(0, fake_next32) == UINT32_MAX && fakePosition == 0,
        "bound 0 consumes nothing");
  CHECK(GCAntiDrawBelow(5, NULL) == UINT32_MAX, "NULL generator");
  set_script(script, 4);
  CHECK(GCAntiDrawBelow(4, fake_next32) == 0 && fakePosition == 1,
        "power of two bound has no rejection");
}

static void fill_state(const unsigned mask, const int nsite, int *eleIdx,
                       int *eleCfg, int *eleNum, int *ncur) {
  int rs;
  int position = 0;
  for (rs = 0; rs < 2 * nsite; rs++) {
    eleNum[rs] = (int)((mask >> rs) & 1U);
    eleCfg[rs] = -1;
  }
  for (rs = 0; rs < 2 * MAX_SITES; rs++) eleIdx[rs] = 1000 + rs;
  /* Interleave spins so the electron order is not up-block first. */
  for (rs = 2 * nsite - 1; rs >= 0; rs--) {
    if (eleNum[rs]) {
      eleIdx[position] = rs;
      eleCfg[rs] = position;
      position++;
    }
  }
  *ncur = position;
}

static void test_find_and_validate(void) {
  int eleIdx[MAX_ORBITALS];
  int eleCfg[MAX_ORBITALS];
  int eleNum[MAX_ORBITALS];
  int ncur;
  /* nsite=4: up {1,3}, down {0,2} -> fused 1,3,4,6 */
  const unsigned mask = (1U << 1) | (1U << 3) | (1U << 4) | (1U << 6);
  fill_state(mask, 4, eleIdx, eleCfg, eleNum, &ncur);
  CHECK(ncur == 4, "fixture particle count");
  CHECK(GCAntiFindOccupied(eleNum, 4, 0, 0) == 1, "up k=0");
  CHECK(GCAntiFindOccupied(eleNum, 4, 0, 1) == 3, "up k=1");
  CHECK(GCAntiFindOccupied(eleNum, 4, 1, 0) == 4, "down k=0");
  CHECK(GCAntiFindOccupied(eleNum, 4, 1, 1) == 6, "down k=1");
  CHECK(GCAntiFindOccupied(eleNum, 4, 0, 2) == -1, "up k out of range");
  CHECK(GCAntiFindOccupied(eleNum, 4, 0, -1) == -1, "negative k");
  CHECK(GCAntiFindOccupied(eleNum, 4, 2, 0) == -1, "spin 2");
  CHECK(GCAntiFindOccupied(eleNum, 4, -1, 0) == -1, "spin -1");
  CHECK(GCAntiFindOccupied(NULL, 4, 0, 0) == -1, "NULL occupation");
  CHECK(GCAntiFindOccupied(eleNum, 0, 0, 0) == -1, "nsite 0");

  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 4) == 1,
        "valid interleaved state");
  /* Scratch beyond ncur is unconstrained. */
  eleIdx[5] = -77;
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 4) == 1,
        "scratch tail ignored");
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 2) == 0,
        "ncur mismatch");
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 3) == 0, "odd ncur");
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 10) == 0,
        "ncur above 2*nsite");
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, -2) == 0,
        "negative ncur");
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 0, 0) == 0, "nsite 0");
  CHECK(GCAntiValidateConfig(NULL, eleCfg, eleNum, 4, 4) == 0, "NULL eleIdx");
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 0x7fffffff, 0) == 0,
        "2*nsite overflow");

  /* Spin imbalance: move down electron 6 to up orbital 2. */
  fill_state((1U << 1) | (1U << 2) | (1U << 3) | (1U << 4), 4, eleIdx, eleCfg,
             eleNum, &ncur);
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 4) == 0,
        "spin imbalance");
  /* Duplicate electron in eleIdx. */
  fill_state(mask, 4, eleIdx, eleCfg, eleNum, &ncur);
  eleIdx[1] = eleIdx[0];
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 4) == 0,
        "duplicate electron");
  /* eleCfg does not point back. */
  fill_state(mask, 4, eleIdx, eleCfg, eleNum, &ncur);
  eleCfg[eleIdx[0]] = 1;
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 4) == 0,
        "eleCfg mismatch");
  /* Out-of-range orbital index. */
  fill_state(mask, 4, eleIdx, eleCfg, eleNum, &ncur);
  eleIdx[2] = 8;
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 4) == 0,
        "orbital index out of range");
  fill_state(mask, 4, eleIdx, eleCfg, eleNum, &ncur);
  eleIdx[2] = -1;
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 4) == 0,
        "negative orbital index");
  /* Occupation outside {0,1}. */
  fill_state(mask, 4, eleIdx, eleCfg, eleNum, &ncur);
  eleNum[0] = 2;
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 4) == 0,
        "occupation 2");
  /* Empty orbital with a stale eleCfg. */
  fill_state(mask, 4, eleIdx, eleCfg, eleNum, &ncur);
  eleCfg[0] = 0;
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 4) == 0,
        "stale eleCfg on an empty orbital");
  /* Occupied orbital not listed in eleIdx. */
  fill_state(mask, 4, eleIdx, eleCfg, eleNum, &ncur);
  eleNum[0] = 1;
  eleNum[5] = 1;
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 4) == 0,
        "occupied orbitals not listed");
  /* Vacuum and full filling are valid. */
  fill_state(0U, 4, eleIdx, eleCfg, eleNum, &ncur);
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 0) == 1, "vacuum");
  fill_state(0xffU, 4, eleIdx, eleCfg, eleNum, &ncur);
  CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, 4, 8) == 1,
        "full filling");
}

/* Reweighted class probabilities used by the GC sampler. */
static void class_probabilities(const int m, const int nsite, double *pHop,
                                double *pAdd, double *pRemove) {
  const int hop = m > 0 && m < nsite;
  const int add = m < nsite;
  const int remove = m > 0;
  const double norm = 0.5 * hop + 0.25 * add + 0.25 * remove;
  *pHop = hop ? 0.5 / norm : 0.0;
  *pAdd = add ? 0.25 / norm : 0.0;
  *pRemove = remove ? 0.25 / norm : 0.0;
}

static int spin_count(const unsigned mask, const int nsite, const int spin) {
  int count = 0;
  int site;
  for (site = 0; site < nsite; site++) {
    count += (int)((mask >> (site + spin * nsite)) & 1U);
  }
  return count;
}

/* q[x][y] for the spin-balanced proposal; returns 0 when consistent. */
static void build_proposal(const int nsite, double q[MAX_STATES][MAX_STATES],
                           int *balanced, int *nbalanced) {
  const unsigned states = 1U << (2 * nsite);
  unsigned x;
  *nbalanced = 0;
  memset(q, 0, sizeof(double) * MAX_STATES * MAX_STATES);
  for (x = 0; x < states; x++) {
    const int m = spin_count(x, nsite, 0);
    double pHop, pAdd, pRemove;
    int s, a, b;
    if (m != spin_count(x, nsite, 1)) continue;
    balanced[(*nbalanced)++] = (int)x;
    class_probabilities(m, nsite, &pHop, &pAdd, &pRemove);
    if (pHop > 0) {
      for (s = 0; s < 2; s++) {
        for (a = 0; a < nsite; a++) {
          for (b = 0; b < nsite; b++) {
            const unsigned from = 1U << (a + s * nsite);
            const unsigned to = 1U << (b + s * nsite);
            if ((x & from) && !(x & to)) {
              q[x][x ^ from ^ to] += pHop / (2.0 * m * (nsite - m));
            }
          }
        }
      }
    }
    for (a = 0; a < nsite; a++) {
      for (b = 0; b < nsite; b++) {
        const unsigned up = 1U << a;
        const unsigned down = 1U << (b + nsite);
        if (pAdd > 0 && !(x & up) && !(x & down)) {
          q[x][x | up | down] += pAdd / ((double)(nsite - m) * (nsite - m));
        }
        if (pRemove > 0 && (x & up) && (x & down)) {
          q[x][x ^ up ^ down] += pRemove / ((double)m * m);
        }
      }
    }
  }
}

static double choose_two(const int n) { return 0.5 * n * (n - 1); }

static void test_detailed_balance(const int nsite) {
  static double q[MAX_STATES][MAX_STATES];
  int balanced[MAX_STATES];
  int nbalanced;
  int ix, iy;
  int choose_two_broken = 0;
  build_proposal(nsite, q, balanced, &nbalanced);
  for (ix = 0; ix < nbalanced; ix++) {
    double total = 0.0;
    const int x = balanced[ix];
    for (iy = 0; iy < nbalanced; iy++) total += q[x][balanced[iy]];
    CHECK(finite_close(total, 1.0, 2e-13), "L=%d state %d proposal sum %.17g",
          nsite, x, total);
  }
  for (ix = 0; ix < nbalanced; ix++) {
    for (iy = 0; iy < nbalanced; iy++) {
      const int x = balanced[ix];
      const int y = balanced[iy];
      const int mx = spin_count((unsigned)x, nsite, 0);
      const int my = spin_count((unsigned)y, nsite, 0);
      double pHopX, pAddX, pRemoveX, pHopY, pAddY, pRemoveY;
      double ratio;
      double pix = 1.0 + x;
      double piy = 1.0 + y;
      double ax, ay;
      int accepted;
      if (x == y || q[x][y] == 0.0) continue;
      CHECK(q[y][x] > 0.0, "L=%d reverse proposal missing %d->%d", nsite, y,
            x);
      class_probabilities(mx, nsite, &pHopX, &pAddX, &pRemoveX);
      class_probabilities(my, nsite, &pHopY, &pAddY, &pRemoveY);
      if (my == mx + 1) {
        ratio = GCAntiRatioAdd(mx, nsite, pAddX, pRemoveY);
        if (fabs((pRemoveY / pAddX) * choose_two(2 * nsite - 2 * mx) /
                     choose_two(2 * mx + 2) -
                 ratio) > 1e-12) {
          choose_two_broken = 1;
        }
      } else if (my == mx - 1) {
        ratio = GCAntiRatioRemove(mx, nsite, pRemoveX, pAddY);
      } else {
        ratio = 1.0;
      }
      CHECK(finite_close(ratio, q[y][x] / q[x][y],
                         2e-13 * (1.0 + q[y][x] / q[x][y])),
            "L=%d ratio %d->%d helper=%.17g enumerated=%.17g", nsite, x, y,
            ratio, q[y][x] / q[x][y]);
      ax = fmin(1.0, piy / pix * ratio);
      ay = fmin(1.0, pix / piy / ratio);
      CHECK(finite_close(pix * q[x][y] * ax, piy * q[y][x] * ay,
                         2e-13 * (1.0 + pix * q[x][y] * ax)),
            "L=%d detailed balance %d<->%d", nsite, x, y);
      /* The log acceptance must realize min(1, weight). */
      if (ax < 1.0) {
        const double oldLog = 0.5 * log(pix);
        const double newLog = 0.5 * log(piy);
        CHECK(GCAntiAcceptLog(oldLog, newLog, 0.0, ratio, ax * (1 - 1e-9),
                              &accepted) == 0 && accepted,
              "L=%d accept below weight %d->%d", nsite, x, y);
        CHECK(GCAntiAcceptLog(oldLog, newLog, 0.0, ratio, ax * (1 + 1e-9),
                              &accepted) == 0 && !accepted,
              "L=%d reject above weight %d->%d", nsite, x, y);
      }
    }
  }
  /* The old all-orbital choose-two count cannot replace the spin count. */
  CHECK(choose_two_broken, "L=%d choose-two ratio unexpectedly agrees", nsite);
}

int main(void) {
  int nsite;
  test_mode_and_ratio();
  test_accept_log();
  test_draw_below();
  test_find_and_validate();
  for (nsite = 2; nsite <= 4; nsite++) test_detailed_balance(nsite);
  if (failures != 0) {
    fprintf(stderr, "GC_AntiPair_Proposal_Unit: %d failure(s)\n", failures);
    return EXIT_FAILURE;
  }
  printf("GC_AntiPair_Proposal_Unit passed\n");
  return EXIT_SUCCESS;
}
