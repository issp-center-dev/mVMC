/*
 * GC Lanczos spin support descriptor (gc_krylov_spin_support.h): canonical
 * inactive form, S = 0 folding (ZERO->TRIPLET, SINGLET->GENERAL, singlet
 * not read), S > 0 four classes, rejection of every non-canonical or
 * out-of-range value without silent folding, and the per-configuration
 * zero reasons against an independent bit-mask reference for every
 * balanced configuration of L = 1..5.  Failures leave the result unchanged.
 */
#include "gc_krylov_spin_support.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int failure_count = 0;

#define CHECK(condition, message)                                    \
  do {                                                               \
    if (!(condition)) {                                              \
      fprintf(stderr, "FAIL: %s (line %d)\n", (message), __LINE__); \
      ++failure_count;                                               \
    }                                                                \
  } while (0)

#define MAX_L 5

static MVMCGCSpinSupport make(int active, int gauss, int spin, int symmetry) {
  MVMCGCSpinSupport support;
  support.active = active;
  support.gauss_count = gauss;
  support.total_spin = spin;
  support.pair_symmetry = (MVMCGCPairSymmetry)symmetry;
  return support;
}

static int validate(MVMCGCSpinSupport support, int site_count,
                    MVMCKrylovConfigurationDomain domain) {
  return mvmc_gc_spin_support_validate(&support, site_count, domain) ==
         MVMC_KRYLOV_STATUS_OK;
}

static int same(MVMCGCSpinSupport left, MVMCGCSpinSupport right) {
  return left.active == right.active &&
         left.gauss_count == right.gauss_count &&
         left.total_spin == right.total_spin &&
         left.pair_symmetry == right.pair_symmetry;
}

static void test_validate(void) {
  const MVMCKrylovConfigurationDomain B = MVMC_KRYLOV_DOMAIN_GC_BALANCED;
  const MVMCKrylovConfigurationDomain E = MVMC_KRYLOV_DOMAIN_GC_EVEN;
  int s, sym;
  CHECK(same(mvmc_gc_spin_support_inactive(), make(0, 1, 0, 0)),
        "canonical inactive value");
  CHECK(validate(make(0, 1, 0, 0), 4, B) && validate(make(0, 1, 0, 0), 4, E),
        "inactive canonical in both GC domains");
  CHECK(!validate(make(0, 0, 0, 0), 4, B) && !validate(make(0, 2, 0, 0), 4, B) &&
            !validate(make(0, 1, 1, 0), 4, B) &&
            !validate(make(0, 1, 0, MVMC_GC_PAIR_TRIPLET), 4, B) &&
            !validate(make(0, -1, 0, 0), 4, B),
        "non-canonical inactive values are rejected, not folded");
  CHECK(!validate(make(0, 1, 0, 0), 4, MVMC_KRYLOV_DOMAIN_FIXED),
        "fixed domain rejected");
  for (s = 0; s <= 2; ++s) {
    for (sym = 0; sym < 4; ++sym) {
      const int expected = s > 0 || sym == MVMC_GC_PAIR_GENERAL ||
                           sym == MVMC_GC_PAIR_TRIPLET;
      CHECK(validate(make(1, 2, s, sym), 4, B) == expected,
            "active S/symmetry combinations (S=0: GENERAL/TRIPLET only)");
      CHECK(!validate(make(1, 2, s, sym), 4, E),
            "active projection requires the balanced domain");
    }
  }
  CHECK(!validate(make(2, 2, 1, 0), 4, B) && !validate(make(-1, 2, 1, 0), 4, B),
        "active flag must be 0 or 1");
  CHECK(!validate(make(1, 1, 1, 0), 4, B) && !validate(make(1, 0, 1, 0), 4, B) &&
            !validate(make(1, -3, 1, 0), 4, B) &&
            validate(make(1, INT_MAX, 1, 0), 4, B),
        "active requires G >= 2 (no upper bound here)");
  CHECK(!validate(make(1, 2, -1, 0), 4, B) && !validate(make(1, 2, 3, 0), 4, B) &&
            !validate(make(1, 2, 3, 0), 5, B) && validate(make(1, 2, 2, 0), 5, B) &&
            !validate(make(1, 2, INT_MAX, 0), 4, B) &&
            !validate(make(1, 2, INT_MIN, 0), 4, B),
        "0 <= S <= floor(L/2)");
  CHECK(!validate(make(1, 2, 1, 4), 4, B) && !validate(make(1, 2, 1, -1), 4, B),
        "symmetry outside the enum rejected");
  CHECK(!validate(make(0, 1, 0, 0), 0, B) && !validate(make(0, 1, 0, 0), -1, B) &&
            !validate(make(0, 1, 0, 0), INT_MAX / 2 + 1, B) &&
            validate(make(0, 1, 0, 0), INT_MAX / 2, B) &&
            !validate(make(0, 1, 0, 0), INT_MIN, B),
        "1 <= L <= INT_MAX/2 checked first");
  CHECK(mvmc_gc_spin_support_validate(NULL, 4, B) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "NULL descriptor");
}

static void test_classify(void) {
  MVMCGCSpinSupport out;
  int triplet, singlet;
  for (triplet = 0; triplet < 2; ++triplet) {
    for (singlet = 0; singlet < 2; ++singlet) {
      CHECK(mvmc_gc_spin_support_classify(1, 4, 0, triplet, singlet, &out) ==
                    MVMC_KRYLOV_STATUS_OK &&
                same(out, make(1, 4, 0,
                               triplet ? MVMC_GC_PAIR_TRIPLET
                                       : MVMC_GC_PAIR_GENERAL)),
            "S=0 folds ZERO->TRIPLET and SINGLET->GENERAL");
      CHECK(mvmc_gc_spin_support_classify(1, 4, 1, triplet, singlet, &out) ==
                    MVMC_KRYLOV_STATUS_OK &&
                (int)out.pair_symmetry ==
                    (triplet && singlet ? MVMC_GC_PAIR_ZERO
                     : triplet          ? MVMC_GC_PAIR_TRIPLET
                     : singlet          ? MVMC_GC_PAIR_SINGLET
                                        : MVMC_GC_PAIR_GENERAL),
            "S>0 uses the four exact classes");
    }
  }
  /* S=0 production reads GCSpinTriplet only: an invalid singlet is ignored. */
  CHECK(mvmc_gc_spin_support_classify(1, 2, 0, 1, -1, &out) ==
                MVMC_KRYLOV_STATUS_OK &&
            same(out, make(1, 2, 0, MVMC_GC_PAIR_TRIPLET)),
        "S=0 never reads the singlet classification");
  out = make(7, 7, 7, 1);
  CHECK(mvmc_gc_spin_support_classify(1, 2, 0, -1, 0, &out) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            same(out, make(7, 7, 7, 1)) &&
            mvmc_gc_spin_support_classify(1, 2, 1, 1, -1, &out) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            mvmc_gc_spin_support_classify(1, 2, 1, 2, 0, &out) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            mvmc_gc_spin_support_classify(2, 2, 1, 1, 0, &out) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            same(out, make(7, 7, 7, 1)),
        "invalid classifications fail without writing");
  CHECK(mvmc_gc_spin_support_classify(0, 1, 3, -1, -1, &out) ==
                MVMC_KRYLOV_STATUS_OK &&
            same(out, mvmc_gc_spin_support_inactive()),
        "G=1 is the canonical inactive value (S is not used)");
  CHECK(mvmc_gc_spin_support_classify(0, 0, 0, 0, 0, &out) ==
                MVMC_KRYLOV_STATUS_OK &&
            !validate(out, 4, MVMC_KRYLOV_DOMAIN_GC_BALANCED),
        "an invalid inactive G is kept for validate to reject");
  CHECK(mvmc_gc_spin_support_classify(1, 0, 1, 0, 0, &out) ==
                MVMC_KRYLOV_STATUS_OK &&
            !validate(out, 4, MVMC_KRYLOV_DOMAIN_GC_BALANCED),
        "an invalid active G is not hidden as inactive");
}

static int popcount(unsigned value) {
  int count = 0;
  while (value) {
    count += (int)(value & 1u);
    value >>= 1;
  }
  return count;
}

/* Independent reference from the mask bits. */
static unsigned reference_reasons(unsigned mask, int L, int S, int symmetry) {
  const unsigned low = (1u << L) - 1u;
  const unsigned up = mask & low, down = (mask >> L) & low;
  const int n = popcount(mask);
  const int singles = popcount(up ^ down);
  unsigned reasons = 0;
  if (S > 0 && singles < 2 * S) reasons |= 1u;
  if ((symmetry == MVMC_GC_PAIR_TRIPLET || symmetry == MVMC_GC_PAIR_ZERO) &&
      ((n / 2 - S) % 2 + 2) % 2 == 1)
    reasons |= 2u;
  if (S > 0 &&
      (symmetry == MVMC_GC_PAIR_SINGLET || symmetry == MVMC_GC_PAIR_ZERO))
    reasons |= 4u;
  return reasons;
}

static void occupation_of(unsigned mask, int L, int *occupation) {
  int orbital;
  for (orbital = 0; orbital < 2 * L; ++orbital) {
    occupation[orbital] = (int)((mask >> orbital) & 1u);
  }
}

static void test_evaluate_exhaustive(void) {
  int occupation[2 * MAX_L];
  int L, S, symmetry;
  unsigned long cases = 0, zeros = 0, distinct_singles = 0;
  for (L = 1; L <= MAX_L; ++L) {
    for (S = 0; S <= L / 2; ++S) {
      for (symmetry = 0; symmetry < 4; ++symmetry) {
        const MVMCGCSpinSupport support = make(1, 3, S, symmetry);
        unsigned mask;
        int seen[2 * MAX_L + 1][MAX_L + 1];
        if (!validate(support, L, MVMC_KRYLOV_DOMAIN_GC_BALANCED)) continue;
        memset(seen, 0, sizeof(seen));
        for (mask = 0; mask < (1u << (2 * L)); ++mask) {
          MVMCGCSpinSupportResult result;
          const unsigned low = (1u << L) - 1u;
          const int balanced =
              popcount(mask & low) == popcount((mask >> L) & low);
          MVMCKrylovStatus status;
          occupation_of(mask, L, occupation);
          memset(&result, 0x5a, sizeof(result));
          status = mvmc_gc_spin_support_evaluate(
              &support, L, MVMC_KRYLOV_DOMAIN_GC_BALANCED, occupation,
              2 * (size_t)L, &result);
          if (!balanced) {
            MVMCGCSpinSupportResult sentinel;
            memset(&sentinel, 0x5a, sizeof(sentinel));
            CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                      memcmp(&result, &sentinel, sizeof(result)) == 0,
                  "unbalanced configuration rejected, result unchanged");
            continue;
          }
          ++cases;
          CHECK(status == MVMC_KRYLOV_STATUS_OK &&
                    result.electron_count == popcount(mask) &&
                    result.single_count ==
                        popcount((mask & low) ^ ((mask >> L) & low)) &&
                    result.zero_reasons ==
                        reference_reasons(mask, L, S, symmetry),
                "zero reasons equal the independent mask reference");
          zeros += result.zero_reasons != 0;
          if (!seen[result.electron_count][result.single_count]++) {
            ++distinct_singles;
          }
        }
      }
    }
  }
  CHECK(cases > 1000 && zeros > 100 && distinct_singles > 50,
        "the enumeration covers zero and nonzero reasons and same-N singles");
  printf("support evaluate: %lu configurations, %lu structural zeros\n", cases,
         zeros);
}

static void test_evaluate_inactive_and_invalid(void) {
  const MVMCGCSpinSupport inactive = mvmc_gc_spin_support_inactive();
  const MVMCGCSpinSupport active = make(1, 2, 1, MVMC_GC_PAIR_TRIPLET);
  MVMCGCSpinSupportResult result, sentinel;
  int occupation[8];
  memset(&sentinel, 0x5a, sizeof(sentinel));
  /* L4 N4 with two doublons: for S=1 triplet both a configuration zero
   * (no singles) and a parity zero (N/2-S=1). */
  occupation_of(0x33u, 4, occupation);
  CHECK(mvmc_gc_spin_support_evaluate(&inactive, 4,
                                      MVMC_KRYLOV_DOMAIN_GC_BALANCED,
                                      occupation, 8, &result) ==
                MVMC_KRYLOV_STATUS_OK &&
            result.zero_reasons == 0 && result.electron_count == 4 &&
            result.single_count == 0,
        "inactive support never masks but still counts");
  CHECK(mvmc_gc_spin_support_evaluate(&active, 4,
                                      MVMC_KRYLOV_DOMAIN_GC_BALANCED,
                                      occupation, 8, &result) ==
                MVMC_KRYLOV_STATUS_OK &&
            result.zero_reasons == (MVMC_GC_SPIN_ZERO_CONFIGURATION |
                                    MVMC_GC_SPIN_ZERO_PARITY),
        "reasons are ORed (configuration and parity)");
  /* GC_EVEN with odd N or a 2 entry is invalid even when inactive. */
  occupation_of(0x1u, 4, occupation);
  result = sentinel;
  CHECK(mvmc_gc_spin_support_evaluate(&inactive, 4, MVMC_KRYLOV_DOMAIN_GC_EVEN,
                                      occupation, 8, &result) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            memcmp(&result, &sentinel, sizeof(result)) == 0,
        "odd N rejected in GC_EVEN, result unchanged");
  occupation_of(0x3u, 4, occupation);
  CHECK(mvmc_gc_spin_support_evaluate(&inactive, 4, MVMC_KRYLOV_DOMAIN_GC_EVEN,
                                      occupation, 8, &result) ==
                MVMC_KRYLOV_STATUS_OK && result.single_count == 2,
        "GC_EVEN accepts an unbalanced even configuration when inactive");
  occupation[5] = 2;
  result = sentinel;
  CHECK(mvmc_gc_spin_support_evaluate(&inactive, 4, MVMC_KRYLOV_DOMAIN_GC_EVEN,
                                      occupation, 8, &result) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            memcmp(&result, &sentinel, sizeof(result)) == 0,
        "occupation other than 0/1 rejected even when inactive");
  occupation[5] = -1;
  CHECK(mvmc_gc_spin_support_evaluate(&inactive, 4, MVMC_KRYLOV_DOMAIN_GC_EVEN,
                                      occupation, 8, &result) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "negative occupation rejected");
  occupation_of(0x11u, 4, occupation);
  CHECK(mvmc_gc_spin_support_evaluate(&active, 4,
                                      MVMC_KRYLOV_DOMAIN_GC_BALANCED,
                                      occupation, 7, &result) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            mvmc_gc_spin_support_evaluate(&active, 4,
                                          MVMC_KRYLOV_DOMAIN_GC_BALANCED,
                                          NULL, 8, &result) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            mvmc_gc_spin_support_evaluate(&active, 4,
                                          MVMC_KRYLOV_DOMAIN_GC_BALANCED,
                                          occupation, 8, NULL) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "occupation count and NULL arguments");
  {
    const MVMCGCSpinSupport bad = make(1, 2, 0, MVMC_GC_PAIR_ZERO);
    result = sentinel;
    CHECK(mvmc_gc_spin_support_evaluate(&bad, 4,
                                        MVMC_KRYLOV_DOMAIN_GC_BALANCED,
                                        occupation, 8, &result) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              memcmp(&result, &sentinel, sizeof(result)) == 0,
          "an unfolded S=0 ZERO descriptor is rejected by evaluate");
  }
  {
    const MVMCGCSpinSupport singlet = make(1, 2, 1, MVMC_GC_PAIR_SINGLET);
    occupation_of(0x21u, 4, occupation); /* up0, down1: two singles */
    CHECK(mvmc_gc_spin_support_all_zero(&singlet) &&
              !mvmc_gc_spin_support_all_zero(&active) &&
              mvmc_gc_spin_support_evaluate(&singlet, 4,
                                            MVMC_KRYLOV_DOMAIN_GC_BALANCED,
                                            occupation, 8, &result) ==
                  MVMC_KRYLOV_STATUS_OK &&
              result.zero_reasons == MVMC_GC_SPIN_ZERO_SINGLET,
          "an S>0 singlet is a zero of every configuration");
  }
}

int main(void) {
  test_validate();
  test_classify();
  test_evaluate_exhaustive();
  test_evaluate_inactive_and_invalid();
  printf(failure_count == 0 ? "gc_krylov_spin_support_unit: PASS\n"
                            : "gc_krylov_spin_support_unit: FAIL\n");
  return failure_count == 0 ? 0 : 1;
}
