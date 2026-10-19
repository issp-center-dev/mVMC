/*
 * Finite total-spin (S>0) GC projection, Task 1/2 of the finite-S plan.
 * The public S>0 input gate stays closed; this Testing-only binary sets the
 * globals directly and drives the production amplitude, Green, derivative,
 * sampler and initialization code.
 *
 *   --probe <fixture> <out>  production values for the independent oracle
 *   --chain <fixture> <seed> <sweeps> <out>
 *                            production chain histogram (oracle compares)
 *   --mask | --transactions | --init | --restart | --support
 *                            self-contained unit checks
 *   --death-*                intended aborts (runtest_gc_antiparallel_death)
 *
 * Fixture text (written by runtest_gc_finite_spin.py):
 *   GC_FINITE_SPIN_FIXTURE_V1
 *   L S G gutzwiller(%a)
 *   L*L lines "re im" (%a), row-major F_ij
 */
#include <complex.h>
#include <math.h>
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
#include "blas_externs.h"
#include "SFMT.h"

#include "../src/mVMC/gc_size.c"
#include "../src/mVMC/workspace.c"
#include "../src/mVMC/gc_component.c"
#include "../src/mVMC/projection.c"
#include "../src/mVMC/gc_config.c"
#include "../src/mVMC/matrix_gc.c"
#include "../src/mVMC/pfupdate_gc.c"
#include "../src/mVMC/splitloop.c"
#include "../src/mVMC/gauleg.c"
#include "../src/mVMC/legendrepoly.c"
#include "../src/mVMC/qp.c"
#include "../src/mVMC/gc_antiparallel.c"
#include "../src/mVMC/vmcmake_gc.c"
#include "../src/mVMC/slater.c"
#include "../src/mVMC/slater_gc.c"
#include "../src/mVMC/locgrn_gc.c"

#define MAX_SITES 8
#define MAX_ORBITALS (2 * MAX_SITES)

static int failures = 0;

#define CHECK(condition, ...)                                                  \
  do {                                                                         \
    if (!(condition)) {                                                        \
      fprintf(stderr, "GC_FS_Unit FAIL: ");                                  \
      fprintf(stderr, __VA_ARGS__);                                            \
      fprintf(stderr, "\n");                                                 \
      failures++;                                                              \
    }                                                                          \
  } while (0)

/* ------------------------------------------------------------------------ */
/* Fixture: identity orbital table, one identity translation, G spin nodes. */
/* ------------------------------------------------------------------------ */
static int *idxStorage, *sgnStorage;
static int *idxRows[MAX_SITES], *sgnRows[MAX_SITES];
static int transRow[MAX_SITES], transSgnRow[MAX_SITES];
static int *transRows[1] = {transRow}, *transSgnRows[1] = {transSgnRow};
static double complex paraTrans[1];

static void release_fixture(void) {
  GCComponentRelease();
  FreeWorkSpaceAll();
  free(idxStorage);
  free(sgnStorage);
  free(Slater);
  free(SlaterElm);
  free(InvM);
  free(PfM);
  free(QPFullWeight);
  free(Proj);
  free(GutzwillerIdx);
  idxStorage = sgnStorage = NULL;
  Slater = NULL;
  SlaterElm = InvM = PfM = QPFullWeight = NULL;
  Proj = NULL;
  GutzwillerIdx = NULL;
}

static void setup_fixture(const int nsite, const int spin, const int gauss,
                          const double complex *pair,
                          const double gutzwiller) {
  int i, j;
  size_t qpCount;
  if (nsite < 1 || nsite > MAX_SITES || gauss < 1) {
    fprintf(stderr, "unsupported fixture L=%d G=%d\n", nsite, gauss);
    exit(EXIT_FAILURE);
  }
  FlagGrandCanonical = 1;
  iFlgOrbitalGeneral = 0;
  TwoSz = 0;
  Nsite = nsite;
  Nsite2 = 2 * nsite;
  NsizeMax = Nsite2;
  Nsize = Nsite2;
  Ne = nsite;
  NSPGaussLeg = gauss;
  NSPStot = spin;
  NMPTrans = 1;
  NQPOptTrans = 1;
  NQPFix = gauss;
  NQPFull = gauss;
  FlagOptTrans = 0;
  NOptTrans = 0;
  NSlater = nsite * nsite;
  NProj = 1;
  NGutzwillerIdx = 1;
  NJastrowIdx = 0;
  NSpinJastrowIdx = 0;
  NDoublonHolon2siteIdx = 0;
  NDoublonHolon4siteIdx = 0;
  NGCInitNelec = -1;
  LapackLWork = 4096;
  NThread = omp_get_max_threads();
  idxStorage = malloc((size_t)NSlater * sizeof(*idxStorage));
  sgnStorage = malloc((size_t)NSlater * sizeof(*sgnStorage));
  Slater = malloc((size_t)NSlater * sizeof(*Slater));
  SlaterElm = malloc((size_t)NQPFull * Nsite2 * Nsite2 * sizeof(*SlaterElm));
  InvM = malloc((size_t)NQPFull * NsizeMax * NsizeMax * sizeof(*InvM));
  PfM = malloc((size_t)NQPFull * sizeof(*PfM));
  qpCount = (size_t)NQPFull + (size_t)NQPFix + 5U * (size_t)gauss;
  QPFullWeight = malloc(qpCount * sizeof(*QPFullWeight));
  Proj = malloc(sizeof(*Proj));
  GutzwillerIdx = calloc((size_t)nsite, sizeof(*GutzwillerIdx));
  if (idxStorage == NULL || sgnStorage == NULL || Slater == NULL ||
      SlaterElm == NULL || InvM == NULL || PfM == NULL ||
      QPFullWeight == NULL || Proj == NULL || GutzwillerIdx == NULL) {
    fprintf(stderr, "fixture allocation failed\n");
    exit(EXIT_FAILURE);
  }
  QPFixWeight = QPFullWeight + NQPFull;
  SPGLCos = QPFixWeight + NQPFix;
  SPGLSin = SPGLCos + gauss;
  SPGLCosSin = SPGLSin + gauss;
  SPGLCosCos = SPGLCosSin + gauss;
  SPGLSinSin = SPGLCosCos + gauss;
  paraTrans[0] = 1.0;
  ParaQPTrans = paraTrans;
  for (i = 0; i < nsite; ++i) {
    transRow[i] = i;
    transSgnRow[i] = 1;
    idxRows[i] = idxStorage + i * nsite;
    sgnRows[i] = sgnStorage + i * nsite;
    for (j = 0; j < nsite; ++j) {
      idxRows[i][j] = i * nsite + j;
      sgnRows[i][j] = 1;
      Slater[i * nsite + j] = pair[i * nsite + j];
    }
  }
  OrbitalIdx = idxRows;
  OrbitalSgn = sgnRows;
  QPTrans = transRows;
  QPTransSgn = transSgnRows;
  QPOptTrans = transRows;
  QPOptTransSgn = transSgnRows;
  Proj[0] = gutzwiller;
  initializeWorkSpaceAll();
  InitQPWeight();
  UpdateSlaterElm_fcmp();
  if (GCComponentInit() != GC_MALL_OK) {
    fprintf(stderr, "component init failed\n");
    exit(EXIT_FAILURE);
  }
}

static double complex *read_fixture(const char *path, int *nsite, int *spin,
                                    int *gauss, double *gutzwiller) {
  char magic[64], a[64], b[64], g[64];
  double complex *pair;
  int k;
  FILE *fp = fopen(path, "r");
  if (fp == NULL || fscanf(fp, "%63s", magic) != 1 ||
      strcmp(magic, "GC_FINITE_SPIN_FIXTURE_V1") != 0 ||
      fscanf(fp, "%d %d %d %63s", nsite, spin, gauss, g) != 4 ||
      *nsite < 1 || *nsite > MAX_SITES || *spin < 0 || *gauss < 1) {
    fprintf(stderr, "bad fixture %s\n", path);
    exit(EXIT_FAILURE);
  }
  *gutzwiller = strtod(g, NULL);
  pair = malloc((size_t)(*nsite) * (size_t)(*nsite) * sizeof(*pair));
  for (k = 0; pair != NULL && k < (*nsite) * (*nsite); ++k) {
    if (fscanf(fp, "%63s %63s", a, b) != 2) {
      fprintf(stderr, "short fixture %s\n", path);
      exit(EXIT_FAILURE);
    }
    pair[k] = strtod(a, NULL) + strtod(b, NULL) * I;
  }
  if (pair == NULL || fscanf(fp, "%63s", magic) != EOF) {
    fprintf(stderr, "trailing or unreadable fixture %s\n", path);
    exit(EXIT_FAILURE);
  }
  fclose(fp);
  return pair;
}

/* Built-in fixtures for self-contained modes (values from the oracle). */
static void builtin_fixture(const char *name, const int gauss) {
  static const double complex l2Triplet[4] = {0, 1, -1, 0};
  static const double complex l2Singlet[4] = {0.5, 0.25, 0.25, -0.5};
  static const double complex l4Max[16] = {0, 1, 2, 3, -2, 1, -1, 1,
                                           1, 3, 0, -2, 2, -1, 4, 1};
  static const double complex l4Triplet[16] = {0, 3, 1, 1, -3, 0, -4, 2,
                                               -1, 4, 0, -6, -1, -2, 6, 0};
  if (strcmp(name, "L2_S1") == 0) {
    setup_fixture(2, 1, gauss, l2Triplet, 0.0);
  } else if (strcmp(name, "L2_S1_singlet") == 0) {
    setup_fixture(2, 1, gauss, l2Singlet, 0.0);
  } else if (strcmp(name, "L4_S2") == 0) {
    setup_fixture(4, 2, gauss, l4Max, -0.3);
  } else if (strcmp(name, "L4_S1_triplet") == 0) {
    setup_fixture(4, 1, gauss, l4Triplet, 0.0);
  } else if (strcmp(name, "L4_S1_general") == 0) {
    /* finite_spin_fixture("L4_S1_general") of the oracle (exact binary
     * fractions), with a Gutzwiller factor. */
    double complex general[16];
    int i, j;
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
        general[4 * i + j] =
            (11 + 3 * i - 2 * j + 17 * (i == j) + 5 * (i * j % 3)) / 32.0 +
            I * ((1 + 2 * i + j - 3 * (i == j)) / 64.0);
    general[1] += 0.75;
    general[8] -= 0.25 * I;
    general[13] += 0.375;
    setup_fixture(4, 1, gauss, general, 0.4);
  } else {
    fprintf(stderr, "unknown builtin fixture %s\n", name);
    exit(EXIT_FAILURE);
  }
}

/* ------------------------------------------------------------------------ */
/* State helpers                                                            */
/* ------------------------------------------------------------------------ */
typedef struct {
  int ncur;
  int eleIdx[MAX_ORBITALS];
  int eleCfg[MAX_ORBITALS];
  int eleNum[MAX_ORBITALS];
  int eleProjCnt[1];
  double complex logIp;
  double complex *inv;
  double complex *pf;
  int *state;
} Snapshot;

static int popcount(unsigned value) {
  int count = 0;
  for (; value != 0U; value >>= 1) count += (int)(value & 1U);
  return count;
}

static int balanced(const unsigned mask) {
  const unsigned low = (1U << Nsite) - 1U;
  return popcount(mask & low) == popcount(mask >> Nsite);
}

static unsigned mask_of(const int *eleNum) {
  unsigned mask = 0U;
  int rs;
  for (rs = 0; rs < Nsite2; ++rs)
    if (eleNum[rs]) mask |= 1U << rs;
  return mask;
}

/* Occupy mask in ascending orbital order and rebuild every component. */
static int load_state(const unsigned mask, int *eleIdx, int *eleCfg,
                      int *eleNum, int *eleProjCnt) {
  int rs, count = 0;
  for (rs = 0; rs < NsizeMax; ++rs) eleIdx[rs] = -1;
  for (rs = 0; rs < Nsite2; ++rs) {
    eleNum[rs] = (int)((mask >> rs) & 1U);
    eleCfg[rs] = -1;
    if (eleNum[rs]) {
      eleIdx[count] = rs;
      eleCfg[rs] = count++;
    }
  }
  Ncur = count;
  MakeProjCnt(eleProjCnt, eleNum);
  return CalculateMAllGC_fcmp(Ncur, eleIdx, 0, NQPFull);
}

static void take_snapshot(Snapshot *s, const int *eleIdx, const int *eleCfg,
                          const int *eleNum, const int *eleProjCnt,
                          const double complex logIp) {
  const size_t inv = (size_t)NQPFull * NsizeMax * NsizeMax;
  s->ncur = Ncur;
  memcpy(s->eleIdx, eleIdx, (size_t)NsizeMax * sizeof(int));
  memcpy(s->eleCfg, eleCfg, (size_t)Nsite2 * sizeof(int));
  memcpy(s->eleNum, eleNum, (size_t)Nsite2 * sizeof(int));
  s->eleProjCnt[0] = eleProjCnt[0];
  s->logIp = logIp;
  s->inv = malloc(inv * sizeof(*s->inv));
  s->pf = malloc((size_t)NQPFull * sizeof(*s->pf));
  s->state = malloc((size_t)NQPFull * sizeof(*s->state));
  if (s->inv == NULL || s->pf == NULL || s->state == NULL) exit(EXIT_FAILURE);
  memcpy(s->inv, InvM, inv * sizeof(*InvM));
  memcpy(s->pf, PfM, (size_t)NQPFull * sizeof(*PfM));
  memcpy(s->state, GCQpState, (size_t)NQPFull * sizeof(*GCQpState));
}

static int same_snapshot(const Snapshot *s, const int *eleIdx,
                         const int *eleCfg, const int *eleNum,
                         const int *eleProjCnt, const double complex logIp) {
  const size_t inv = (size_t)NQPFull * NsizeMax * NsizeMax;
  return s->ncur == Ncur &&
         memcmp(s->eleIdx, eleIdx, (size_t)NsizeMax * sizeof(int)) == 0 &&
         memcmp(s->eleCfg, eleCfg, (size_t)Nsite2 * sizeof(int)) == 0 &&
         memcmp(s->eleNum, eleNum, (size_t)Nsite2 * sizeof(int)) == 0 &&
         s->eleProjCnt[0] == eleProjCnt[0] &&
         memcmp(&s->logIp, &logIp, sizeof(logIp)) == 0 &&
         memcmp(s->inv, InvM, inv * sizeof(*InvM)) == 0 &&
         memcmp(s->pf, PfM, (size_t)NQPFull * sizeof(*PfM)) == 0 &&
         memcmp(s->state, GCQpState, (size_t)NQPFull * sizeof(int)) == 0;
}

static void free_snapshot(Snapshot *s) {
  free(s->inv);
  free(s->pf);
  free(s->state);
}

static void allocate_green_scratch(GCGreenScratch *scratch) {
  const int maxN = 2;
  const int small = 4 * maxN * maxN;
  scratch->maxN = maxN;
  scratch->pfLWork = small;
  scratch->projCntNew = malloc((size_t)NProj * sizeof(int));
  scratch->rsi = malloc((size_t)maxN * sizeof(int));
  scratch->rsj = malloc((size_t)maxN * sizeof(int));
  scratch->msa = malloc((size_t)maxN * sizeof(int));
  scratch->pfIWork = malloc((size_t)(2 * maxN) * sizeof(int));
  scratch->pfMNew = malloc((size_t)NQPFull * sizeof(double complex));
  scratch->vec =
      malloc((size_t)maxN * (size_t)NsizeMax * sizeof(double complex));
  scratch->w = malloc((size_t)NsizeMax * sizeof(double complex));
  scratch->smallMat = malloc((size_t)small * sizeof(double complex));
  scratch->pfWork = malloc((size_t)small * sizeof(double complex));
  scratch->pfRWork = malloc((size_t)small * sizeof(double));
  if (scratch->projCntNew == NULL || scratch->rsi == NULL ||
      scratch->rsj == NULL || scratch->msa == NULL ||
      scratch->pfIWork == NULL || scratch->pfMNew == NULL ||
      scratch->vec == NULL || scratch->w == NULL ||
      scratch->smallMat == NULL || scratch->pfWork == NULL ||
      scratch->pfRWork == NULL)
    exit(EXIT_FAILURE);
}

static void free_green_scratch(GCGreenScratch *scratch) {
  free(scratch->projCntNew);
  free(scratch->rsi);
  free(scratch->rsj);
  free(scratch->msa);
  free(scratch->pfIWork);
  free(scratch->pfMNew);
  free(scratch->vec);
  free(scratch->w);
  free(scratch->smallMat);
  free(scratch->pfWork);
  free(scratch->pfRWork);
}

static void put_complex(FILE *fp, const double complex value) {
  fprintf(fp, " %a %a", creal(value), cimag(value));
}

/* ------------------------------------------------------------------------ */
/* --probe: every balanced configuration of the fixture                    */
/* ------------------------------------------------------------------------ */
static int run_probe(const char *fixturePath, const char *outPath) {
  int nsite, spin, gauss, q;
  double gutzwiller;
  double complex *pair =
      read_fixture(fixturePath, &nsite, &spin, &gauss, &gutzwiller);
  FILE *fp;
  unsigned mask;
  int eleIdx[MAX_ORBITALS], eleCfg[MAX_ORBITALS], eleNum[MAX_ORBITALS];
  int eleProjCnt[1];
  double complex *srOptO;
  GCGreenScratch scratch;
  setup_fixture(nsite, spin, gauss, pair, gutzwiller);
  free(pair);
  srOptO = malloc((size_t)(2 * NSlater) * sizeof(*srOptO));
  fp = fopen(outPath, "w");
  if (fp == NULL || srOptO == NULL) return EXIT_FAILURE;
  allocate_green_scratch(&scratch);
  fprintf(fp, "GC_FINITE_SPIN_PROBE_V1 %d %d %d %d %d %d\n", Nsite,
          NSPStot, NSPGaussLeg, GCFiniteSpinActive(), GCSpinTriplet,
          GCSpinSinglet);
  for (q = 0; q < NQPFull; ++q) {
    fprintf(fp, "Q %d %a %a %a", q, creal(SPGLCosSin[q]),
            creal(SPGLCosCos[q]), creal(SPGLSinSin[q]));
    put_complex(fp, QPFullWeight[q]);
    fprintf(fp, "\n");
  }
  for (mask = 0U; mask < (1U << Nsite2); ++mask) {
    int status, singles, k;
    double complex ip;
    if (!balanced(mask)) continue;
    status = load_state(mask, eleIdx, eleCfg, eleNum, eleProjCnt);
    CHECK(status == GC_MALL_OK, "probe rebuild mask=%u", mask);
    singles = GCSpinCountSingles(eleNum, Nsite);
    ip = CalculateIPGC_fcmp(PfM, Ncur, eleNum, 0, NQPFull, MPI_COMM_SELF);
    fprintf(fp, "A %u %d %d", mask, Ncur, singles);
    put_complex(fp, ip);
    fprintf(fp, "\n");
    CHECK(isfinite(creal(ip)) && isfinite(cimag(ip)),
          "nonfinite probe amplitude mask=%u", mask);
    if (ip != 0.0) {
      int i, j, a, b;
      SlaterElmDiffGC_fcmp(srOptO, ip, eleIdx, Ncur);
      for (k = 0; k < NSlater; ++k) {
        fprintf(fp, "D %u %d", mask, k);
        put_complex(fp, srOptO[2 * k]);
        put_complex(fp, srOptO[2 * k + 1]);
        fprintf(fp, "\n");
      }
      for (i = 0; i < Nsite2; ++i) {
        for (j = 0; j < Nsite2; ++j) {
          fprintf(fp, "G1 %u %d %d", mask, i, j);
          put_complex(fp, GreenFunc1GC(i, j, ip, Ncur, eleIdx, eleCfg,
                                       eleNum, eleProjCnt, &scratch));
          fprintf(fp, "\nPA %u %d %d", mask, i, j);
          put_complex(fp, GreenFuncPairAddGC(i, j, ip, Ncur, eleIdx, eleCfg,
                                             eleNum, eleProjCnt, &scratch));
          fprintf(fp, "\nPR %u %d %d", mask, i, j);
          put_complex(fp,
                      GreenFuncPairRemoveGC(i, j, ip, Ncur, eleIdx, eleCfg,
                                            eleNum, eleProjCnt, &scratch));
          fprintf(fp, "\n");
        }
      }
      /* Sz-conserving two-body products: up-up x down-down and the
       * spin-flip pairs c+_{a up} c_{b dn} c+_{i dn} c_{j up}. */
      for (a = 0; a < Nsite; ++a)
        for (b = 0; b < Nsite; ++b)
          for (i = 0; i < Nsite; ++i)
            for (j = 0; j < Nsite; ++j) {
              const int set[2][4] = {
                  {a, b, i + Nsite, j + Nsite},
                  {a, b + Nsite, i + Nsite, j}};
              int kind;
              for (kind = 0; kind < 2; ++kind) {
                fprintf(fp, "G2 %u %d %d %d %d", mask, set[kind][0],
                        set[kind][1], set[kind][2], set[kind][3]);
                put_complex(fp, GreenFunc2GC(set[kind][0], set[kind][1],
                                             set[kind][2], set[kind][3], ip,
                                             Ncur, eleIdx, eleCfg, eleNum,
                                             eleProjCnt, &scratch));
                fprintf(fp, "\n");
              }
            }
      /* The Green functions restore the configuration exactly. */
      CHECK(mask_of(eleNum) == mask, "Green changed the configuration");
    } else if (GCSpinTriplet == 1 && singles >= 2 * NSPStot) {
      /* Pure-triplet pair-parity zero: the raw, unmasked absolute
       * derivative (ip = 1) is the derivative along symmetry-breaking
       * directions; it must not be multiplied by the parameter flag. */
      SlaterElmDiffGC_fcmp(srOptO, 1.0, eleIdx, Ncur);
      for (k = 0; k < NSlater; ++k) {
        fprintf(fp, "Z %u %d", mask, k);
        put_complex(fp, srOptO[2 * k]);
        fprintf(fp, "\n");
      }
    }
  }
  free_green_scratch(&scratch);
  free(srOptO);
  if (fclose(fp) != 0) return EXIT_FAILURE;
  release_fixture();
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}

/* ------------------------------------------------------------------------ */
/* --chain: production sampler histogram                                   */
/* ------------------------------------------------------------------------ */
static int run_chain(const char *fixturePath, const char *seedText,
                     const char *sweepText, const char *outPath) {
  int nsite, spin, gauss, sweep, step;
  double gutzwiller;
  double complex *pair =
      read_fixture(fixturePath, &nsite, &spin, &gauss, &gutzwiller);
  const unsigned seed = (unsigned)strtoul(seedText, NULL, 10);
  const long sweeps = strtol(sweepText, NULL, 10);
  int eleIdx[MAX_ORBITALS], eleCfg[MAX_ORBITALS], eleNum[MAX_ORBITALS];
  int eleProjCnt[1], projCntNew[1];
  double complex *pfMNew, logIp;
  long *histogram;
  long accepted = 0;
  unsigned mask;
  FILE *fp;
  setup_fixture(nsite, spin, gauss, pair, gutzwiller);
  free(pair);
  pfMNew = malloc((size_t)NQPFull * sizeof(*pfMNew));
  histogram = calloc((size_t)1 << Nsite2, sizeof(*histogram));
  if (pfMNew == NULL || histogram == NULL || sweeps <= 0) return EXIT_FAILURE;
  init_gen_rand(seed);
  Ncur = Nsite;
  (void)makeInitialSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt, 0, NQPFull,
                            MPI_COMM_SELF);
  logIp = CalculateLogIPGC_fcmp(PfM, Ncur, eleNum, 0, NQPFull, MPI_COMM_SELF);
  for (sweep = 0; sweep < sweeps; ++sweep) {
    for (step = 0; step < Nsite; ++step) {
      accepted += GCMakeOneStep(eleIdx, eleCfg, eleNum, eleProjCnt, &logIp,
                                pfMNew, projCntNew, 0, NQPFull,
                                MPI_COMM_SELF);
    }
    CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, Nsite, Ncur),
          "chain state invalid");
    CHECK(!GCSpinConfigurationZero(NSPStot,
                                   GCSpinCountSingles(eleNum, Nsite)),
          "chain visited a structural zero");
    histogram[mask_of(eleNum)]++;
  }
  fp = fopen(outPath, "w");
  if (fp == NULL) return EXIT_FAILURE;
  fprintf(fp, "GC_FINITE_SPIN_CHAIN_V1 %u %ld %ld\n", seed, sweeps, accepted);
  for (mask = 0U; mask < (1U << Nsite2); ++mask)
    if (histogram[mask] != 0) fprintf(fp, "H %u %ld\n", mask, histogram[mask]);
  if (fclose(fp) != 0) return EXIT_FAILURE;
  free(histogram);
  free(pfMNew);
  release_fixture();
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}

/* ------------------------------------------------------------------------ */
/* --mask: wrapper contract                                                 */
/* ------------------------------------------------------------------------ */
static void test_mask(void) {
  int eleIdx[MAX_ORBITALS], eleCfg[MAX_ORBITALS], eleNum[MAX_ORBITALS];
  int eleProjCnt[1];
  double complex values[64];
  unsigned mask;
  int q;
  builtin_fixture("L4_S1_triplet", 4);
  CHECK(GCFiniteSpinActive() && GCSpinTriplet == 1 && GCSpinSinglet == 0,
        "triplet fixture classification %d %d", GCSpinTriplet,
        GCSpinSinglet);
  for (q = 0; q < NQPFull; ++q) values[q] = 1.0e-16 * (q + 1);
  for (mask = 0U; mask < (1U << Nsite2); ++mask) {
    int n, singles, zero;
    double complex ip;
    if (!balanced(mask)) continue;
    (void)load_state(mask, eleIdx, eleCfg, eleNum, eleProjCnt);
    n = Ncur;
    singles = GCSpinCountSingles(eleNum, Nsite);
    zero = singles < 2 || (n / 2 - 1) % 2 != 0;
    CHECK(GCProjectionSectorZero(n, eleNum) == zero,
          "sector zero mask=%u N=%d singles=%d", mask, n, singles);
    ip = CalculateIPGC_fcmp(values, n, eleNum, 0, NQPFull, MPI_COMM_SELF);
    CHECK(zero ? ip == 0.0 : ip != 0.0,
          "proven zero must remove only finite residue mask=%u", mask);
    values[0] = NAN;
    ip = CalculateIPGC_fcmp(values, n, eleNum, 0, NQPFull, MPI_COMM_SELF);
    CHECK(!isfinite(creal(ip)), "nonfinite hidden by the mask at %u", mask);
    values[0] = 1.0e-16;
    /* Invalid candidates: a particle count different from the occupation,
     * a missing occupation and a non-0/1 entry are errors, not zeros. */
    CHECK(GCProjectionSectorZero(n + 2, eleNum) == -1, "N mismatch");
    CHECK(GCProjectionSectorZero(n, NULL) == -1, "NULL occupation");
    CHECK(isnan(creal(CalculateIPGC_fcmp(values, n + 2, eleNum, 0, NQPFull,
                                         MPI_COMM_SELF))),
          "N mismatch did not give NaN");
  }
  eleNum[0] = 2;
  CHECK(GCProjectionSectorZero(Ncur, eleNum) == -1, "occupation 2 accepted");
  release_fixture();

  /* S=0 keeps the Phase 3a contract and needs no occupation. */
  builtin_fixture("L4_S1_triplet", 4);
  NSPStot = 0;
  UpdateSlaterElm_fcmp();
  CHECK(!GCFiniteSpinActive() && GCSpinSinglet == 0 && GCSpinTriplet == 1,
        "S=0 classification");
  for (q = 0; q <= 8; q += 2)
    CHECK(GCProjectionSectorZero(q, NULL) == (q % 4 == 2), "S=0 N=%d", q);
  release_fixture();

  /* G=1 is no projection: S>0 has no effect and the occupation is unused. */
  builtin_fixture("L4_S1_triplet", 1);
  CHECK(!GCFiniteSpinActive() && GCSpinTriplet == 0 && GCSpinSinglet == 0,
        "G=1 must not classify");
  for (q = 0; q <= 8; q += 2)
    CHECK(GCProjectionSectorZero(q, NULL) == 0, "G=1 N=%d", q);
  release_fixture();
}

/* ------------------------------------------------------------------------ */
/* --transactions: rejected bridges leave every array; accepted ones equal */
/* a full rebuild of the new configuration.                                */
/* ------------------------------------------------------------------------ */
static void check_rebuilt(const int *eleIdx, const int *eleNum,
                          const double complex logIp, const char *label) {
  const size_t pfBytes = (size_t)NQPFull * sizeof(*PfM);
  double complex *pf = malloc(pfBytes);
  double complex again;
  int q;
  if (pf == NULL) exit(EXIT_FAILURE);
  memcpy(pf, PfM, pfBytes);
  CHECK(CalculateMAllGC_fcmp(Ncur, eleIdx, 0, NQPFull) == GC_MALL_OK,
        "%s rebuild", label);
  for (q = 0; q < NQPFull; ++q)
    CHECK(cabs(pf[q] - PfM[q]) <= 1.0e-10 * (1.0 + cabs(PfM[q])),
          "%s Pfaffian differs from a full rebuild at qp=%d", label, q);
  again = CalculateLogIPGC_fcmp(PfM, Ncur, eleNum, 0, NQPFull, MPI_COMM_SELF);
  CHECK(cabs(cexp(again - logIp) - 1.0) <= 1.0e-10,
        "%s amplitude differs from a full rebuild", label);
  free(pf);
}

static void test_transactions(const char *fixture) {
  int eleIdx[MAX_ORBITALS], eleCfg[MAX_ORBITALS], eleNum[MAX_ORBITALS];
  int eleProjCnt[1], projCntNew[1];
  double complex *pfMNew, logIp;
  long counts[4][2] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}};
  int kind, trial;
  builtin_fixture(fixture, 6);
  pfMNew = malloc((size_t)NQPFull * sizeof(*pfMNew));
  if (pfMNew == NULL) exit(EXIT_FAILURE);
  init_gen_rand(4049U);
  Ncur = Nsite;
  (void)makeInitialSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt, 0, NQPFull,
                            MPI_COMM_SELF);
  logIp = CalculateLogIPGC_fcmp(PfM, Ncur, eleNum, 0, NQPFull, MPI_COMM_SELF);
  for (trial = 0; trial < 400; ++trial) {
    for (kind = 0; kind < 4; ++kind) {
      static const enum GCAntiBridgeMove moves[4] = {
          GC_ANTI_SWAP, GC_ANTI_ADD_TWO, GC_ANTI_REMOVE_TWO, GC_ANTI_GLOBAL};
      Snapshot before;
      int accepted;
      take_snapshot(&before, eleIdx, eleCfg, eleNum, eleProjCnt, logIp);
      accepted = GCAttemptBridge(moves[kind], eleIdx, eleCfg, eleNum,
                                 eleProjCnt, &logIp, pfMNew, projCntNew, 0,
                                 NQPFull, MPI_COMM_SELF);
      counts[kind][accepted != 0]++;
      if (accepted) {
        CHECK(GCAntiValidateConfig(eleIdx, eleCfg, eleNum, Nsite, Ncur),
              "accepted bridge state invalid");
        CHECK(!GCSpinConfigurationZero(
                  NSPStot, GCSpinCountSingles(eleNum, Nsite)),
              "accepted a structural zero");
        if (moves[kind] == GC_ANTI_SWAP)
          CHECK(Ncur == before.ncur &&
                    GCSpinCountSingles(eleNum, Nsite) ==
                        GCSpinCountSingles(before.eleNum, Nsite),
                "swap changed N or singles");
        check_rebuilt(eleIdx, eleNum, logIp, "accepted bridge");
      } else {
        CHECK(same_snapshot(&before, eleIdx, eleCfg, eleNum, eleProjCnt,
                            logIp),
              "rejected bridge kind=%d changed the sampler state", kind);
      }
      free_snapshot(&before);
    }
  }
  /* Equal-magnitude spin states (maximal S, a two-single triplet) accept
   * every swap; the general fixture with four singles must also reject. */
  CHECK(counts[0][1] > 0, "no accepted swap");
  if (strcmp(fixture, "L4_S1_general") == 0)
    CHECK(counts[0][0] > 0, "general swap never rejected");
  if (strcmp(fixture, "L4_S2") != 0) {
    CHECK(counts[1][1] + counts[2][1] > 0, "no accepted N+-4 move");
    CHECK(counts[3][1] > 0, "no accepted global move");
  }
  printf("transactions %s swap %ld/%ld add %ld/%ld remove %ld/%ld "
         "global %ld/%ld (rejected/accepted)\n",
         fixture, counts[0][0], counts[0][1], counts[1][0], counts[1][1],
         counts[2][0], counts[2][1], counts[3][0], counts[3][1]);
  free(pfMNew);
  release_fixture();
}

/* ------------------------------------------------------------------------ */
/* --init / --restart                                                       */
/* ------------------------------------------------------------------------ */
static void test_init(void) {
  int eleIdx[MAX_ORBITALS], eleCfg[MAX_ORBITALS], eleNum[MAX_ORBITALS];
  int eleProjCnt[1], trial, n;
  long seen[256] = {0};
  /* L4 S=2: only N=4 with four singles; the default N=Nsite=4 is valid. */
  builtin_fixture("L4_S2", 4);
  init_gen_rand(77U);
  for (trial = 0; trial < 3000; ++trial) {
    Ncur = 4;
    NGCInitNelec = trial % 2 ? 4 : -1;
    CHECK(makeInitialSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt, 0, NQPFull,
                              MPI_COMM_SELF) == 0, "L4 S2 init");
    CHECK(Ncur == 4 && GCSpinCountSingles(eleNum, Nsite) == 4,
          "L4 S2 init N=%d", Ncur);
    seen[mask_of(eleNum)]++;
  }
  {
    static const unsigned support[6] = {60, 90, 105, 150, 165, 195};
    long total = 0;
    int k;
    for (k = 0; k < 6; ++k) {
      /* 3000 draws over 6 equally likely states: 500 +- ~5 sigma. */
      CHECK(labs(seen[support[k]] - 500) < 110, "init histogram %u: %ld",
            support[k], seen[support[k]]);
      total += seen[support[k]];
    }
    CHECK(total == 3000, "init left the maximal-S support");
  }
  release_fixture();

  /* L4 triplet S=1: N=4 is forbidden. Unspecified initial N moves to the
   * nearest permitted sector (tie -> the smaller N=2). */
  builtin_fixture("L4_S1_triplet", 4);
  init_gen_rand(91U);
  Ncur = 4;
  NGCInitNelec = -1;
  CHECK(makeInitialSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt, 0, NQPFull,
                            MPI_COMM_SELF) == 0, "auto init");
  CHECK(Ncur == 2, "auto initial N %d (expected 2)", Ncur);
  for (n = 2; n <= 6; n += 4) {
    Ncur = n;
    NGCInitNelec = n;
    CHECK(makeInitialSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt, 0, NQPFull,
                              MPI_COMM_SELF) == 0 && Ncur == n,
          "explicit initial N=%d", n);
  }
  /* Ordering helper: distance from the center, ties to the smaller N. */
  {
    int order[16], count = GCFiniteSpinAppendByDistance(4, order, 0);
    CHECK(count == 2 && order[0] == 2 && order[1] == 6,
          "distance order %d: %d %d", count, order[0], order[1]);
    order[0] = 6;
    count = GCFiniteSpinAppendByDistance(4, order, 1);
    CHECK(count == 2 && order[0] == 6 && order[1] == 2,
          "listed sectors must not repeat");
  }
  release_fixture();
}

/* A parameter generation that turns the saved N into a forbidden sector:
 * restore from a burn sample at N=4 of a general F, then make F an exact
 * triplet (N=4 is now zero for S=1). The sampler restarts in a permitted N
 * with the full warm-up and records only permitted samples. */
static void test_restart(void) {
  static const double complex general[16] = {
      0.5, 3, 1, 1, -3, 0.25, -4, 2, -1, 4, 0, -6, -1, -2, 6, 0.125};
  static const double complex triplet[16] = {0, 3, 1, 1, -3, 0, -4, 2,
                                             -1, 4, 0, -6, -1, -2, 6, 0};
  int eleIdx[MAX_ORBITALS], eleCfg[MAX_ORBITALS], eleNum[MAX_ORBITALS];
  int eleProjCnt[1], k;
  setup_fixture(4, 1, 4, general, 0.0);
  CHECK(GCSpinTriplet == 0, "general fixture classified as triplet");
  init_gen_rand(313U);
  /* mask 0x96 = up1,up2 + down0,down3: N=4 with four singles. */
  CHECK(load_state(0x96U, eleIdx, eleCfg, eleNum, eleProjCnt) == GC_MALL_OK,
        "restart state");
  CHECK(isfinite(creal(CalculateLogIPGC_fcmp(PfM, Ncur, eleNum, 0, NQPFull,
                                             MPI_COMM_SELF))),
        "general N=4 state is zero");
  BurnEleIdx = malloc(MAX_ORBITALS * sizeof(int));
  BurnEleCfg = malloc(MAX_ORBITALS * sizeof(int));
  BurnEleNum = malloc(MAX_ORBITALS * sizeof(int));
  BurnEleProjCnt = malloc(sizeof(int));
  TmpEleIdx = malloc(MAX_ORBITALS * sizeof(int));
  TmpEleCfg = malloc(MAX_ORBITALS * sizeof(int));
  TmpEleNum = malloc(MAX_ORBITALS * sizeof(int));
  TmpEleProjCnt = malloc(sizeof(int));
  EleIdx = malloc(8 * MAX_ORBITALS * sizeof(int));
  EleCfg = malloc(8 * MAX_ORBITALS * sizeof(int));
  EleNum = malloc(8 * MAX_ORBITALS * sizeof(int));
  EleProjCnt = malloc(8 * sizeof(int));
  EleNumSample = malloc(8 * sizeof(int));
  logSqPfFullSlater = malloc(8 * sizeof(double));
  copyToBurnSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt);
  memcpy(Slater, triplet, sizeof(triplet));
  UpdateSlaterElm_fcmp();
  CHECK(GCSpinTriplet == 1, "triplet generation not detected");
  BurnFlag = 1;
  NGCInitNelec = 6;
  NVMCWarmUp = 50;
  NVMCInterval = 1;
  NVMCSample = 8;
  VMCMakeSampleGC(MPI_COMM_SELF);
  for (k = 0; k < 8; ++k)
    CHECK(EleNumSample[k] == 2 || EleNumSample[k] == 6,
          "restart sampled forbidden N=%d", EleNumSample[k]);
  CHECK(BurnFlag == 1, "burn flag after the restarted run");
  free(BurnEleIdx);
  free(BurnEleCfg);
  free(BurnEleNum);
  free(BurnEleProjCnt);
  free(TmpEleIdx);
  free(TmpEleCfg);
  free(TmpEleNum);
  free(TmpEleProjCnt);
  free(EleIdx);
  free(EleCfg);
  free(EleNum);
  free(EleProjCnt);
  free(EleNumSample);
  free(logSqPfFullSlater);
  release_fixture();
}

/* ------------------------------------------------------------------------ */
/* --support: production proposal kernels and connectivity mutations        */
/* ------------------------------------------------------------------------ */
static uint32_t support_state = 2463534242U;
static uint32_t support_word(void) {
  support_state ^= support_state << 13;
  support_state ^= support_state >> 17;
  support_state ^= support_state << 5;
  return support_state;
}

static int to_occupation(unsigned mask, int nsite, int *occupation) {
  int rs;
  for (rs = 0; rs < 2 * nsite; ++rs) occupation[rs] = (int)((mask >> rs) & 1U);
  return popcount(mask & ((1U << nsite) - 1U));
}

/* Candidates reached from mask by the local hop/N+-2 kernel (same-spin hop,
 * up-down pair add/remove), enumerated directly. */
static int local_neighbour(unsigned from, unsigned to, int nsite) {
  const unsigned diff = from ^ to;
  const unsigned low = (1U << nsite) - 1U;
  const int changed = popcount(diff);
  if (changed == 2) {
    const unsigned up = diff & low, down = diff >> nsite;
    if (popcount(up) == 2 || popcount(down) == 2)
      return popcount(from & diff) == 1; /* same-spin hop */
    return popcount(up) == 1 && popcount(down) == 1 &&
           (popcount(from & diff) == 0 || popcount(from & diff) == 2);
  }
  return 0;
}

#define KIND_LOCAL 1
#define KIND_SWAP 2
#define KIND_TWO 4
#define KIND_GLOBAL 8

/* Production draws (GCAntiDrawBridge) define the swap/two-pair/global edges;
 * enough draws per state make every proposal edge appear. */
static int components(const unsigned *support, const int count,
                      const int nsite, const int kinds) {
  int parent[64], i, j, comps = count;
  int current[MAX_ORBITALS], candidate[MAX_ORBITALS];
  for (i = 0; i < count; ++i) parent[i] = i;
  for (i = 0; i < count; ++i) {
    const int m = to_occupation(support[i], nsite, current);
    int kind, draw;
    for (j = 0; j < count; ++j) {
      if (j != i && (kinds & KIND_LOCAL) &&
          local_neighbour(support[i], support[j], nsite)) {
        int a = i, b = j;
        while (parent[a] != a) a = parent[a];
        while (parent[b] != b) b = parent[b];
        if (a != b) { parent[a] = b; --comps; }
      }
    }
    for (kind = 0; kind < 4; ++kind) {
      static const enum GCAntiBridgeMove moves[4] = {
          GC_ANTI_SWAP, GC_ANTI_ADD_TWO, GC_ANTI_REMOVE_TWO, GC_ANTI_GLOBAL};
      static const int bits[4] = {KIND_SWAP, KIND_TWO, KIND_TWO, KIND_GLOBAL};
      if (!(kinds & bits[kind])) continue;
      for (draw = 0; draw < 20000; ++draw) {
        double ratio;
        unsigned target = 0U;
        int rs;
        if (GCAntiDrawBridge(moves[kind], current, nsite, m, candidate,
                             support_word, &ratio) != 1)
          break;
        for (rs = 0; rs < 2 * nsite; ++rs)
          if (candidate[rs]) target |= 1U << rs;
        for (j = 0; j < count; ++j) {
          if (support[j] == target) {
            int a = i, b = j;
            while (parent[a] != a) a = parent[a];
            while (parent[b] != b) b = parent[b];
            if (a != b) { parent[a] = b; --comps; }
          }
        }
      }
    }
  }
  return comps;
}

static void test_support(void) {
  /* (a) maximal S: no local edge; the swap connects (L=2 S=1, L=4 S=2). */
  static const unsigned l2[2] = {6, 9};
  static const unsigned l4[6] = {60, 90, 105, 150, 165, 195};
  /* (b) L=4 S=1 triplet support: N=2 and N=6 with at least two singles. */
  unsigned triplet[64];
  /* (c) abstract L=6 target of plan section 6.3. */
  static const unsigned abstract[2] = {129, 3966};
  int count = 0, all = KIND_LOCAL | KIND_SWAP | KIND_TWO | KIND_GLOBAL;
  unsigned mask;
  int occupation[MAX_ORBITALS];
  Nsite = 4;
  for (mask = 0U; mask < 256U; ++mask) {
    const int n = popcount(mask);
    if (!balanced(mask) || (n != 2 && n != 6)) continue;
    (void)to_occupation(mask, 4, occupation);
    if (GCSpinCountSingles(occupation, 4) >= 2) triplet[count++] = mask;
  }
  CHECK(count == 24, "triplet support size %d", count);
  Nsite = 2;
  CHECK(components(l2, 2, 2, KIND_LOCAL | KIND_TWO) == 2,
        "L2 maximal S connected without swap (mutation a)");
  CHECK(components(l2, 2, 2, KIND_LOCAL | KIND_SWAP | KIND_TWO) == 1,
        "L2 swap does not connect");
  Nsite = 4;
  CHECK(components(l4, 6, 4, KIND_LOCAL | KIND_TWO) == 6,
        "L4 maximal S has local edges (mutation a)");
  CHECK(components(l4, 6, 4, KIND_LOCAL | KIND_SWAP | KIND_TWO) == 1,
        "L4 swap does not connect");
  CHECK(components(triplet, count, 4, KIND_LOCAL | KIND_SWAP) == 2,
        "triplet N=2/6 joined without two-pair move (mutation b)");
  CHECK(components(triplet, count, 4, KIND_LOCAL | KIND_SWAP | KIND_TWO) == 1,
        "two-pair move does not join N=2/6");
  Nsite = 6;
  CHECK(components(abstract, 2, 6, KIND_LOCAL | KIND_SWAP | KIND_TWO) == 2,
        "abstract target joined without global (mutation c)");
  CHECK(components(abstract, 2, 6, all) == 1, "global does not join");

  /* Production mixture: fixed probabilities summing to one. */
  CHECK(GCAntiSelectFiniteSpinBridge(0.0) == GC_ANTI_GLOBAL &&
            GCAntiSelectFiniteSpinBridge(nextafter(1.0 / 64, 0)) ==
                GC_ANTI_GLOBAL &&
            GCAntiSelectFiniteSpinBridge(1.0 / 64) == GC_ANTI_ADD_TWO &&
            GCAntiSelectFiniteSpinBridge(9.0 / 64) == GC_ANTI_REMOVE_TWO &&
            GCAntiSelectFiniteSpinBridge(17.0 / 64) == GC_ANTI_SWAP &&
            GCAntiSelectFiniteSpinBridge(nextafter(25.0 / 64, 0)) ==
                GC_ANTI_SWAP &&
            GCAntiSelectFiniteSpinBridge(25.0 / 64) == -1 &&
            GCAntiSelectFiniteSpinBridge(nextafter(1.0, 0)) == -1,
        "finite-S mixture boundaries");
  /* S=0 mixture unchanged (Phase 3 random-number contract). */
  CHECK(GCAntiSelectBridge(17.0 / 64) == -1 &&
            GCAntiSelectBridge(nextafter(17.0 / 64, 0)) ==
                GC_ANTI_REMOVE_TWO,
        "S=0 mixture changed");
}

/* ------------------------------------------------------------------------ */
/* --death-*                                                                */
/* ------------------------------------------------------------------------ */
static int run_death(const char *mode) {
  int eleIdx[MAX_ORBITALS], eleCfg[MAX_ORBITALS], eleNum[MAX_ORBITALS];
  int eleProjCnt[1], projCntNew[1], rank = 0, size = 1;
  double complex pfMNew[64], logIp;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  init_gen_rand(60013U);
  if (strcmp(mode, "--death-singlet") == 0) {
    builtin_fixture("L2_S1_singlet", 4);
    GCRequireSpinSupport(MPI_COMM_WORLD);
  } else if (strcmp(mode, "--death-explicit-n") == 0) {
    builtin_fixture("L4_S1_triplet", 4);
    Ncur = 4;
    NGCInitNelec = 4;
    (void)makeInitialSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt, 0, NQPFull,
                              MPI_COMM_WORLD);
  } else if (strcmp(mode, "--death-explicit-range") == 0) {
    builtin_fixture("L4_S2", 4);
    Ncur = 2;
    NGCInitNelec = 2;
    (void)makeInitialSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt, 0, NQPFull,
                              MPI_COMM_WORLD);
  } else if (strcmp(mode, "--death-rank-singlet") == 0) {
    builtin_fixture("L4_S1_triplet", 4);
    GCSpinSinglet = rank == 0 ? 0 : 1;
    GCRequireSpinSupport(MPI_COMM_WORLD);
  } else if (strcmp(mode, "--death-swap-nonfinite") == 0) {
    int start, end;
    int k;
    /* G=2 leaves an empty QP range on the third rank of np3. */
    builtin_fixture("L4_S2", 2);
    SplitLoop(&start, &end, NQPFull, rank, size);
    (void)load_state(60U, eleIdx, eleCfg, eleNum, eleProjCnt);
    CHECK(CalculateMAllGC_fcmp(Ncur, eleIdx, start, end) == GC_MALL_OK,
          "distributed rebuild");
    logIp = CalculateLogIPGC_fcmp(PfM, Ncur, eleNum, start, end,
                                  MPI_COMM_WORLD);
    /* A nonfinite component on one rank, including an empty-QP rank, must
     * reach every rank before the projected sum is evaluated. */
    if (rank == (size > 1 ? 1 : 0) && start < end)
      for (k = 0; k < Nsite2 * Nsite2; ++k)
        SlaterElm[(size_t)start * Nsite2 * Nsite2 + k] = NAN;
    (void)GCAttemptBridge(GC_ANTI_SWAP, eleIdx, eleCfg, eleNum, eleProjCnt,
                          &logIp, pfMNew, projCntNew, start, end,
                          MPI_COMM_WORLD);
  } else if (strcmp(mode, "--death-restart-same-generation") == 0) {
    /* A stored state that is zero although the pair function did not
     * change (mask 51 has two doublons, impossible at maximal S) must abort
     * instead of being hidden by the generation-change restart. */
    static int burnIdx[MAX_ORBITALS], burnCfg[MAX_ORBITALS],
        burnNum[MAX_ORBITALS], burnProj[1], tmpIdx[MAX_ORBITALS],
        tmpCfg[MAX_ORBITALS], tmpNum[MAX_ORBITALS], tmpProj[1];
    builtin_fixture("L4_S2", 4);
    BurnEleIdx = burnIdx; BurnEleCfg = burnCfg; BurnEleNum = burnNum;
    BurnEleProjCnt = burnProj; TmpEleIdx = tmpIdx; TmpEleCfg = tmpCfg;
    TmpEleNum = tmpNum; TmpEleProjCnt = tmpProj;
    (void)load_state(51U, eleIdx, eleCfg, eleNum, eleProjCnt);
    copyToBurnSampleGC(eleIdx, eleCfg, eleNum, eleProjCnt);
    BurnFlag = 1;
    NVMCWarmUp = 0;
    NVMCSample = 1;
    NVMCInterval = 1;
    VMCMakeSampleGC(MPI_COMM_WORLD);
  } else {
    fprintf(stderr, "unknown death mode %s\n", mode);
    return 2;
  }
  fprintf(stderr, "death mode %s returned without aborting\n", mode);
  return 3;
}

int main(int argc, char **argv) {
  const char *mode = argc >= 2 ? argv[1] : "";
  int status = EXIT_SUCCESS;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#endif
  if (strncmp(mode, "--death-", 8) == 0) return run_death(mode);
  if (strcmp(mode, "--probe") == 0 && argc == 4) {
    status = run_probe(argv[2], argv[3]);
  } else if (strcmp(mode, "--chain") == 0 && argc == 6) {
    status = run_chain(argv[2], argv[3], argv[4], argv[5]);
  } else if (strcmp(mode, "--mask") == 0) {
    test_mask();
  } else if (strcmp(mode, "--transactions") == 0) {
    test_transactions("L4_S2");
    test_transactions("L4_S1_triplet");
    test_transactions("L4_S1_general");
  } else if (strcmp(mode, "--init") == 0) {
    test_init();
  } else if (strcmp(mode, "--restart") == 0) {
    test_restart();
  } else if (strcmp(mode, "--support") == 0) {
    test_support();
  } else {
    fprintf(stderr, "usage: %s --probe|--chain|--mask|--transactions|"
            "--init|--restart|--support|--death-*\n", argv[0]);
    status = 2;
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  if (status != EXIT_SUCCESS) return status;
  if (failures != 0) {
    fprintf(stderr, "GC_FS_Unit %s: %d failure(s)\n", mode, failures);
    return EXIT_FAILURE;
  }
  printf("GC_FS_Unit %s: PASS\n", mode);
  return EXIT_SUCCESS;
}
