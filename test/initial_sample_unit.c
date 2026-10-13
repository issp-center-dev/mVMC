#include <complex.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

extern int omp_get_thread_num(void);
extern int omp_get_max_threads(void);

#include "global.h"
#include "SFMT.h"

#include "../src/mVMC/initial_sample.c"

#define SITE_MAX 64
#define SEED_COUNT 200

static int failures = 0;

#define CHECK(condition, ...)                                                  \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "InitialSample_Unit FAIL: ");                          \
      fprintf(stderr, __VA_ARGS__);                                            \
      fprintf(stderr, "\n");                                                   \
      failures++;                                                              \
    }                                                                          \
  } while (0)

static int locSpn[SITE_MAX];

/* Sites 0, 2, 4, ... have a local spin, up to nLoc of them. */
static void setup(const int nsite, const int nLoc, const int ne,
                  const int twoSz) {
  int ri, count = 0;
  Nsite = nsite;
  Nsite2 = 2 * nsite;
  Ne = ne;
  Nsize = 2 * ne;
  TwoSz = twoSz;
  LocSpn = locSpn;
  for (ri = 0; ri < nsite; ri++) {
    locSpn[ri] = (ri % 2 == 0 && count < nLoc) ? 1 : 0;
    count += locSpn[ri];
  }
}

static int up_total(const int fsz) {
  return fsz ? Ne + ((TwoSz == -1) ? 0 : TwoSz / 2) : Ne;
}

/* The local spins as they were placed before the fix */
static void legacy(const int fsz, int *eleIdx, int *eleCfg, int *eleSpn) {
  int ri, mi, si, msi, rsi;
  for (msi = 0; msi < Nsize; msi++) eleIdx[msi] = -1;
  for (rsi = 0; rsi < Nsite2; rsi++) eleCfg[rsi] = -1;
  if (fsz) {
    for (msi = 0; msi < Nsize; msi++) {
      eleSpn[msi] = (msi < up_total(fsz)) ? 0 : 1;
    }
  }
  for (ri = 0; ri < Nsite; ri++) {
    if (LocSpn[ri] == 1) {
      if (fsz) {
        do {
          msi = gen_rand32() % Nsize;
          si = eleSpn[msi];
        } while (eleIdx[msi] != -1);
        eleCfg[ri + si * Nsite] = msi;
        eleIdx[msi] = ri;
      } else {
        do {
          mi = gen_rand32() % Ne;
          si = (genrand_real2() < 0.5) ? 0 : 1;
        } while (eleIdx[mi + si * Ne] != -1);
        eleCfg[ri + si * Nsite] = mi;
        eleIdx[mi + si * Ne] = ri;
      }
    }
  }
}

/* Number of the electrons of spin si which are not placed yet */
static int rest(const int fsz, const int *eleIdx, const int si) {
  const int nUp = up_total(fsz);
  int msi, count = 0;
  for (msi = 0; msi < Nsize; msi++) {
    if (eleIdx[msi] == -1 && ((msi < nUp) ? 0 : 1) == si) count++;
  }
  return count;
}

static int fits(const int fsz, const int *eleIdx) {
  int ri, nFree = 0;
  for (ri = 0; ri < Nsite; ri++) {
    if (LocSpn[ri] != 1) nFree++;
  }
  return rest(fsz, eleIdx, 0) <= nFree && rest(fsz, eleIdx, 1) <= nFree;
}

/* Every local spin has one electron, the other sites have none, and eleIdx
   and eleCfg agree with each other. */
static void check_valid(const char *name, const int seed, const int fsz,
                        const int *eleIdx, const int *eleCfg,
                        const int *eleSpn) {
  const int nUp = up_total(fsz);
  int ri, si, msi, placed = 0;
  for (ri = 0; ri < Nsite; ri++) {
    const int up = eleCfg[ri];
    const int down = eleCfg[ri + Nsite];
    if (LocSpn[ri] == 1) {
      CHECK((up == -1) != (down == -1),
            "%s seed %d: site %d has %s", name, seed, ri,
            (up == -1) ? "no electron" : "two electrons");
    } else {
      CHECK(up == -1 && down == -1,
            "%s seed %d: site %d without local spin is occupied", name, seed,
            ri);
    }
    for (si = 0; si < 2; si++) {
      const int label = eleCfg[ri + si * Nsite];
      if (label == -1) continue;
      msi = fsz ? label : label + si * Ne;
      CHECK(label >= 0 && msi < Nsize && eleIdx[msi] == ri,
            "%s seed %d: eleCfg and eleIdx disagree at site %d", name, seed,
            ri);
      CHECK(((msi < nUp) ? 0 : 1) == si,
            "%s seed %d: electron %d has the wrong spin", name, seed, msi);
      placed++;
    }
  }
  for (msi = 0; msi < Nsize; msi++) {
    if (eleIdx[msi] != -1) placed--;
    if (fsz) {
      CHECK(eleSpn[msi] == ((msi < nUp) ? 0 : 1),
            "%s seed %d: eleSpn[%d] is wrong", name, seed, msi);
    }
  }
  CHECK(placed == 0, "%s seed %d: eleIdx has electrons which eleCfg has not",
        name, seed);
  CHECK(fits(fsz, eleIdx),
        "%s seed %d: %d up and %d down itinerant electrons do not fit", name,
        seed, rest(fsz, eleIdx, 0), rest(fsz, eleIdx, 1));
}

/* Returns the number of seeds with which the draw before the fix did not
   fit, i.e. with which vmc.out did not end. */
static int test_case(const char *name, const int fsz, const int nTryMax) {
  int seed, status, nNotFit = 0;
  for (seed = 1; seed <= SEED_COUNT; seed++) {
    int oldIdx[2 * SITE_MAX], oldCfg[2 * SITE_MAX], oldSpn[2 * SITE_MAX];
    int newIdx[2 * SITE_MAX], newCfg[2 * SITE_MAX], newSpn[2 * SITE_MAX];
    uint32_t oldNext, newNext;

    init_gen_rand(seed);
    legacy(fsz, oldIdx, oldCfg, oldSpn);
    oldNext = gen_rand32();

    init_gen_rand(seed);
    status = InitialSamplePlaceLocalSpinTry(newIdx, newCfg,
                                            fsz ? newSpn : NULL, nTryMax);
    newNext = gen_rand32();

    CHECK(status == 0, "%s seed %d: status %d", name, seed, status);
    if (status != 0) continue;
    check_valid(name, seed, fsz, newIdx, newCfg, newSpn);

    if (!fits(fsz, oldIdx)) {
      nNotFit++;
    } else if (nTryMax > 0) {
      /* nothing changes if the draw before the fix fits */
      CHECK(memcmp(oldIdx, newIdx, sizeof(int) * Nsize) == 0 &&
                memcmp(oldCfg, newCfg, sizeof(int) * Nsite2) == 0 &&
                (!fsz || memcmp(oldSpn, newSpn, sizeof(int) * Nsize) == 0),
            "%s seed %d: the sample differs from that before the fix", name,
            seed);
      CHECK(oldNext == newNext,
            "%s seed %d: random numbers are consumed differently", name,
            seed);
    }
  }
  return nNotFit;
}

static void test_unchanged(void) {
  /* Kondo chain, 4 local spins and 4 conduction sites, half filling and
     below: the draw always fits */
  setup(8, 4, 4, 0);
  CHECK(test_case("half_filling", 0, 100) == 0, "half_filling: no fit");
  CHECK(test_case("half_filling_fsz", 1, 100) == 0, "half_filling_fsz: no fit");
  setup(8, 4, 3, 0);
  CHECK(test_case("below_half_filling", 0, 100) == 0,
        "below_half_filling: no fit");
  /* Heisenberg chain */
  setup(8, 8, 4, 0);
  CHECK(test_case("spin", 0, 100) == 0, "spin: no fit");
  CHECK(test_case("spin_fsz", 1, 100) == 0, "spin_fsz: no fit");
  setup(8, 8, 4, -1);
  CHECK(test_case("spin_fsz_sz_free", 1, 100) == 0, "spin_fsz_sz_free: no fit");
  /* Hubbard chain: no local spin, no random number is used */
  setup(6, 0, 3, 0);
  CHECK(test_case("hubbard", 0, 100) == 0, "hubbard: no fit");
  CHECK(test_case("hubbard_fsz", 1, 100) == 0, "hubbard_fsz: no fit");
}

static void test_above_half_filling(void) {
  /* These are the inputs with which vmc.out did not end for some seeds. */
  setup(8, 4, 5, 0); /* Ncond = 6 */
  CHECK(test_case("ncond6", 0, 100) > 0, "ncond6: the draw always fits");
  CHECK(test_case("ncond6_fsz", 1, 100) > 0,
        "ncond6_fsz: the draw always fits");
  setup(8, 4, 6, 0); /* Ncond = 8, the conduction band is full */
  CHECK(test_case("ncond8", 0, 100) > 0, "ncond8: the draw always fits");
  CHECK(test_case("ncond8_fsz", 1, 100) > 0,
        "ncond8_fsz: the draw always fits");
  setup(8, 4, 6, 4); /* all the local spins must be up */
  CHECK(test_case("ncond8_fsz_2sz4", 1, 100) > 0,
        "ncond8_fsz_2sz4: the draw always fits");
}

static void test_restricted_draw(void) {
  /* nTryMax = 0: the local spins are placed with the restricted number of
     up spins from the beginning */
  setup(8, 4, 5, 0);
  test_case("restricted_ncond6", 0, 0);
  test_case("restricted_ncond6_fsz", 1, 0);
  setup(8, 4, 6, 0);
  test_case("restricted_ncond8", 0, 0);
  test_case("restricted_ncond8_fsz", 1, 0);
  setup(8, 4, 6, 4);
  test_case("restricted_ncond8_fsz_2sz4", 1, 0);
  setup(8, 8, 4, 0);
  test_case("restricted_spin", 0, 0);
  setup(6, 0, 3, 0);
  test_case("restricted_hubbard", 0, 0);
  /* 24 local spins which must be up altogether: the draw at random fits
     with the probability C(48,24)/C(72,24) = 2e-6, so that the restricted
     draw is used after 100 draws */
  setup(48, 24, 36, 24);
  CHECK(test_case("large_fsz", 1, 100) == SEED_COUNT,
        "large_fsz: the draw at random fits");
}

static void test_no_fit(void) {
  int eleIdx[2 * SITE_MAX], eleCfg[2 * SITE_MAX], eleSpn[2 * SITE_MAX];
  /* Ncond = Nsite with local spins on all the sites */
  setup(4, 4, 4, 0);
  CHECK(InitialSamplePlaceLocalSpin(eleIdx, eleCfg, NULL) == 1,
        "spin with Ncond = Nsite is accepted");
  CHECK(InitialSamplePlaceLocalSpin(eleIdx, eleCfg, eleSpn) == 1,
        "spin with Ncond = Nsite is accepted (FSZ)");
  /* more electrons than twice the sites */
  setup(6, 0, 7, 0);
  CHECK(InitialSamplePlaceLocalSpin(eleIdx, eleCfg, NULL) == 1,
        "Ne > Nsite is accepted");
  /* 8 up electrons in 6 sites */
  setup(6, 0, 5, 6);
  CHECK(InitialSamplePlaceLocalSpin(eleIdx, eleCfg, eleSpn) == 1,
        "Ne + Sz > Nsite is accepted (FSZ)");
  /* the conduction band is full and the local spins are not */
  setup(8, 4, 6, 6);
  CHECK(InitialSamplePlaceLocalSpin(eleIdx, eleCfg, eleSpn) == 1,
        "2Sz > NLocalSpin with the full conduction band is accepted (FSZ)");
}

int main(void) {
  test_unchanged();
  test_above_half_filling();
  test_restricted_draw();
  test_no_fit();
  if (failures != 0) {
    fprintf(stderr, "InitialSample_Unit: %d failure(s)\n", failures);
    return EXIT_FAILURE;
  }
  printf("InitialSample_Unit: PASS\n");
  return EXIT_SUCCESS;
}
