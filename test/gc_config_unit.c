#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gc_config.h"

#define MAX_ORBITALS 8

static int failures = 0;

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "GCConfig_Unit FAIL: ");                                 \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, "\n");                                                   \
      failures++;                                                               \
    }                                                                           \
  } while (0)

static void swap_int(int *left, int *right) {
  const int temporary = *left;
  *left = *right;
  *right = temporary;
}

static int permutation_sign_by_cycles(const int *permutation, const int n) {
  int visited[MAX_ORBITALS] = {0};
  int i;
  int sign = 1;
  for (i = 0; i < n; i++) {
    int cycleLength = 0;
    int j = i;
    if (visited[i]) continue;
    while (!visited[j]) {
      visited[j] = 1;
      j = permutation[j];
      cycleLength++;
    }
    if (cycleLength % 2 == 0) sign = -sign;
  }
  return sign;
}

static int parity_oracle(const int pos0, const int pos1, const int ncurOld) {
  int permutation[MAX_ORBITALS];
  int i;
  for (i = 0; i < ncurOld; i++) permutation[i] = i;
  if (pos1 != ncurOld - 1) {
    swap_int(&permutation[pos1], &permutation[ncurOld - 1]);
  }
  if (pos0 != ncurOld - 2) {
    swap_int(&permutation[pos0], &permutation[ncurOld - 2]);
  }
  return permutation_sign_by_cycles(permutation, ncurOld);
}

static void initialize_ordered_configuration(const int ncur, int *eleIdx,
                                             int *eleCfg, int *eleNum) {
  static const int ordering[MAX_ORBITALS] = {7, 0, 5, 2, 6, 1, 4, 3};
  int i;
  for (i = 0; i < MAX_ORBITALS; i++) {
    eleIdx[i] = ordering[i];
    eleCfg[i] = -1;
    eleNum[i] = 0;
  }
  for (i = 0; i < ncur; i++) {
    eleCfg[eleIdx[i]] = i;
    eleNum[eleIdx[i]] = 1;
  }
}

static void test_remove_parity_and_restore(void) {
  int ncurOld;
  for (ncurOld = 2; ncurOld <= MAX_ORBITALS; ncurOld += 2) {
    int pos0;
    for (pos0 = 0; pos0 < ncurOld; pos0++) {
      int pos1;
      for (pos1 = pos0 + 1; pos1 < ncurOld; pos1++) {
        int eleIdx[MAX_ORBITALS];
        int eleCfg[MAX_ORBITALS];
        int eleNum[MAX_ORBITALS];
        int originalIdx[MAX_ORBITALS];
        int originalCfg[MAX_ORBITALS];
        int originalNum[MAX_ORBITALS];
        int expectedIdx[MAX_ORBITALS];
        int expectedCfg[MAX_ORBITALS];
        int expectedNum[MAX_ORBITALS];
        int ncur = ncurOld;
        int i;
        int sign;

        initialize_ordered_configuration(ncurOld, eleIdx, eleCfg, eleNum);
        memcpy(originalIdx, eleIdx, sizeof(eleIdx));
        memcpy(originalCfg, eleCfg, sizeof(eleCfg));
        memcpy(originalNum, eleNum, sizeof(eleNum));
        memcpy(expectedIdx, eleIdx, sizeof(eleIdx));
        if (pos1 != ncurOld - 1) {
          swap_int(&expectedIdx[pos1], &expectedIdx[ncurOld - 1]);
        }
        if (pos0 != ncurOld - 2) {
          swap_int(&expectedIdx[pos0], &expectedIdx[ncurOld - 2]);
        }
        for (i = 0; i < MAX_ORBITALS; i++) {
          expectedCfg[i] = -1;
          expectedNum[i] = 0;
        }
        for (i = 0; i < ncurOld - 2; i++) {
          expectedCfg[expectedIdx[i]] = i;
          expectedNum[expectedIdx[i]] = 1;
        }

        sign = GCRemovePair(pos0, pos1, eleIdx, eleCfg, eleNum, &ncur);
        CHECK(sign == parity_oracle(pos0, pos1, ncurOld),
              "parity n=%d positions=(%d,%d): got %d expected %d", ncurOld,
              pos0, pos1, sign, parity_oracle(pos0, pos1, ncurOld));
        CHECK(GCRemoveParitySign(pos0, pos1, ncurOld) ==
                  parity_oracle(pos0, pos1, ncurOld),
              "standalone parity n=%d positions=(%d,%d)", ncurOld, pos0,
              pos1);
        CHECK(ncur == ncurOld - 2, "remove count n=%d got %d", ncurOld,
              ncur);
        CHECK(memcmp(eleIdx, expectedIdx, sizeof(eleIdx)) == 0,
              "remove permutation n=%d positions=(%d,%d)", ncurOld, pos0,
              pos1);
        CHECK(memcmp(eleCfg, expectedCfg, sizeof(eleCfg)) == 0,
              "remove inverse map n=%d positions=(%d,%d)", ncurOld, pos0,
              pos1);
        CHECK(memcmp(eleNum, expectedNum, sizeof(eleNum)) == 0,
              "remove occupation n=%d positions=(%d,%d)", ncurOld, pos0,
              pos1);

        GCRestoreRemovedPair(pos0, pos1, eleIdx, eleCfg, eleNum, &ncur);
        CHECK(ncur == ncurOld, "restore count n=%d got %d", ncurOld, ncur);
        CHECK(memcmp(eleIdx, originalIdx, sizeof(eleIdx)) == 0,
              "restore list n=%d positions=(%d,%d)", ncurOld, pos0, pos1);
        CHECK(memcmp(eleCfg, originalCfg, sizeof(eleCfg)) == 0,
              "restore inverse map n=%d positions=(%d,%d)", ncurOld, pos0,
              pos1);
        CHECK(memcmp(eleNum, originalNum, sizeof(eleNum)) == 0,
              "restore occupation n=%d positions=(%d,%d)", ncurOld, pos0,
              pos1);
      }
    }
  }
}

static void test_add_remove_round_trip(void) {
  int eleIdx[MAX_ORBITALS];
  int eleCfg[MAX_ORBITALS];
  int eleNum[MAX_ORBITALS];
  int originalIdx[MAX_ORBITALS];
  int originalCfg[MAX_ORBITALS];
  int originalNum[MAX_ORBITALS];
  int ncur = 4;
  int i;

  initialize_ordered_configuration(ncur, eleIdx, eleCfg, eleNum);
  memcpy(originalIdx, eleIdx, sizeof(eleIdx));
  memcpy(originalCfg, eleCfg, sizeof(eleCfg));
  memcpy(originalNum, eleNum, sizeof(eleNum));
  GCAddPair(6, 3, eleIdx, eleCfg, eleNum, &ncur);
  CHECK(ncur == 6, "add count");
  CHECK(eleIdx[4] == 6 && eleIdx[5] == 3, "add append ordering");
  CHECK(eleCfg[6] == 4 && eleCfg[3] == 5, "add inverse map");
  CHECK(eleNum[6] == 1 && eleNum[3] == 1, "add occupation");
  CHECK(GCRemovePair(4, 5, eleIdx, eleCfg, eleNum, &ncur) == 1,
        "removing appended tail has positive parity");
  CHECK(ncur == 4, "add/remove count round trip");
  CHECK(memcmp(eleIdx, originalIdx, (size_t)ncur * sizeof(int)) == 0,
        "add/remove active list round trip");
  for (i = 0; i < MAX_ORBITALS; i++) {
    CHECK(eleCfg[i] == originalCfg[i], "add/remove eleCfg rs=%d", i);
    CHECK(eleNum[i] == originalNum[i], "add/remove eleNum rs=%d", i);
  }
}

static int popcount_int(const unsigned int value) {
  unsigned int work = value;
  int count = 0;
  while (work != 0) {
    count += (int)(work & 1U);
    work >>= 1;
  }
  return count;
}

static double choose_two(const int count) {
  return 0.5 * (double)count * (double)(count - 1);
}

static double class_probability(const int ncur, const int nsite2,
                                const int addClass) {
  const int hopPossible = ncur > 0 && ncur < nsite2;
  const int addPossible = nsite2 - ncur >= 2;
  const int removePossible = ncur >= 2;
  double normalization = 0.0;
  if (hopPossible) normalization += 0.50;
  if (addPossible) normalization += 0.25;
  if (removePossible) normalization += 0.25;
  if (normalization == 0.0) return 0.0;
  if (addClass) return addPossible ? 0.25 / normalization : 0.0;
  return removePossible ? 0.25 / normalization : 0.0;
}

static void test_kth_empty(void) {
  unsigned int mask;
  for (mask = 0; mask < (1U << MAX_ORBITALS); mask++) {
    int eleNum[MAX_ORBITALS];
    int expected[MAX_ORBITALS];
    int emptyCount = 0;
    int rs;
    for (rs = 0; rs < MAX_ORBITALS; rs++) {
      eleNum[rs] = (int)((mask >> rs) & 1U);
      if (eleNum[rs] == 0) expected[emptyCount++] = rs;
    }
    for (rs = 0; rs < emptyCount; rs++) {
      CHECK(GCFindKthEmpty(eleNum, MAX_ORBITALS, rs) == expected[rs],
            "k-th empty mask=%u k=%d", mask, rs);
    }
    CHECK(GCFindKthEmpty(eleNum, MAX_ORBITALS, -1) == -1,
          "negative empty index mask=%u", mask);
    CHECK(GCFindKthEmpty(eleNum, MAX_ORBITALS, emptyCount) == -1,
          "out-of-range empty index mask=%u", mask);
  }
}

static void test_physical_state_detailed_balance(void) {
  static const unsigned int states[8] = {0U, 3U, 5U, 6U, 9U, 10U, 12U, 15U};
  int sawMutationFailure = 0;
  int ix;
  for (ix = 0; ix < 8; ix++) {
    int iy;
    const int nx = popcount_int(states[ix]);
    const double pix = 1.0 + 0.17 * (double)(ix + 1);
    for (iy = 0; iy < 8; iy++) {
      const int ny = popcount_int(states[iy]);
      const double piy = 1.0 + 0.17 * (double)(iy + 1);
      const unsigned int added = states[iy] & ~states[ix];
      const unsigned int removed = states[ix] & ~states[iy];
      if (ny == nx + 2 && popcount_int(added) == 2 && removed == 0U) {
        const double pAddX = class_probability(nx, 4, 1);
        const double pRemoveY = class_probability(ny, 4, 0);
        const double qForward = pAddX / choose_two(4 - nx);
        const double qReverse = pRemoveY / choose_two(ny);
        const double ratio = GCProposalRatioAdd(nx, 4, pAddX, pRemoveY);
        const double reverseRatio =
            GCProposalRatioRemove(ny, 4, pRemoveY, pAddX);
        const double acceptanceForward = fmin(1.0, piy / pix * ratio);
        const double acceptanceReverse =
            fmin(1.0, pix / piy * reverseRatio);
        const double lhs = pix * qForward * acceptanceForward;
        const double rhs = piy * qReverse * acceptanceReverse;
        const double mutatedLhs = pix * qForward * fmin(1.0, piy / pix);
        const double mutatedRhs = piy * qReverse * fmin(1.0, pix / piy);

        CHECK(fabs(ratio - qReverse / qForward) < 1.0e-14,
              "add proposal ratio states %u -> %u", states[ix], states[iy]);
        CHECK(fabs(reverseRatio - qForward / qReverse) < 1.0e-14,
              "remove proposal ratio states %u -> %u", states[iy],
              states[ix]);
        CHECK(fabs(ratio * reverseRatio - 1.0) < 1.0e-14,
              "reciprocal proposal ratios states %u <-> %u", states[ix],
              states[iy]);
        CHECK(fabs(lhs - rhs) < 1.0e-13,
              "detailed balance states %u <-> %u: %.17g vs %.17g",
              states[ix], states[iy], lhs, rhs);
        if (fabs(mutatedLhs - mutatedRhs) > 1.0e-8) sawMutationFailure = 1;
      }
    }
  }
  CHECK(sawMutationFailure,
        "proposalRatio=1 mutation unexpectedly preserves all detailed balance edges");
}

int main(void) {
  test_remove_parity_and_restore();
  test_add_remove_round_trip();
  test_kth_empty();
  test_physical_state_detailed_balance();
  if (failures != 0) {
    fprintf(stderr, "%d GC configuration checks failed\n", failures);
    return EXIT_FAILURE;
  }
  printf("GC configuration checks: PASS\n");
  return EXIT_SUCCESS;
}
