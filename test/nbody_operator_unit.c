#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "backflow_nbody.h"
#include "nbody_operator.h"

int NExUpdatePath = 0;
int Nsite = 0;

#include "sector_projection.h"

static int failures = 0;
static unsigned long kindCounts[4] = {0, 0, 0, 0};

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "NBodyOperator_Unit FAIL: ");                             \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, "\n");                                                    \
      failures++;                                                               \
    }                                                                           \
  } while (0)

static int PopcountBelow(uint64_t state, int orbital) {
  uint64_t mask;
  if (orbital <= 0) return 0;
  mask = (UINT64_C(1) << orbital) - UINT64_C(1);
  return __builtin_popcountll(state & mask);
}

static int ApplyOpsBitstring(uint64_t state, int n,
                             const int *rsi, const int *rsj,
                             uint64_t *target, int *sign) {
  int k;
  int value = 1;

  if (n <= 0 || rsi == NULL || rsj == NULL ||
      target == NULL || sign == NULL) {
    return 0;
  }

  for (k = n - 1; k >= 0; k--) {
    const uint64_t annihilationBit = UINT64_C(1) << rsj[k];
    const uint64_t creationBit = UINT64_C(1) << rsi[k];

    if ((state & annihilationBit) == 0) return 0;
    if (PopcountBelow(state, rsj[k]) & 1) value = -value;
    state &= ~annihilationBit;

    if ((state & creationBit) != 0) return 0;
    if (PopcountBelow(state, rsi[k]) & 1) value = -value;
    state |= creationBit;
  }

  *target = state;
  *sign = value;
  return 1;
}

static void FillOccupation(uint64_t state, int nOrbitals, int *occupation) {
  int orbital;
  for (orbital = 0; orbital < nOrbitals; orbital++) {
    occupation[orbital] = (int)((state >> orbital) & UINT64_C(1));
  }
}

static void CheckOneReduction(uint64_t state, int n, int nOrbitals,
                              const int *rsi, const int *rsj,
                              const char *label) {
  int occupation[8];
  int rsiInput[8];
  int rsjInput[8];
  int rsiWork[8];
  int rsjWork[8];
  uint64_t originalTarget = 0;
  uint64_t reducedTarget = 0;
  int originalSign = 0;
  int reducedSign = 0;
  int originalNonzero;
  int reducedNonzero;
  int status;
  NBodyReduction reduction;

  CHECK(n <= 8 && nOrbitals <= 8, "%s: test dimensions exceed fixture", label);
  if (n > 8 || nOrbitals > 8) return;

  memcpy(rsiInput, rsi, (size_t)n * sizeof(int));
  memcpy(rsjInput, rsj, (size_t)n * sizeof(int));
  FillOccupation(state, nOrbitals, occupation);
  originalNonzero = ApplyOpsBitstring(state, n, rsi, rsj,
                                      &originalTarget, &originalSign);

  memset(rsiWork, 0x5a, sizeof(rsiWork));
  memset(rsjWork, 0xa5, sizeof(rsjWork));
  status = ReduceNBodyTerm(n, rsi, rsj, occupation, nOrbitals,
                           rsiWork, rsjWork, n, &reduction);
  CHECK(status == 0, "%s: reducer returned status %d", label, status);
  if (status != 0) return;

  CHECK(memcmp(rsiInput, rsi, (size_t)n * sizeof(int)) == 0,
        "%s: reducer modified rsi input", label);
  CHECK(memcmp(rsjInput, rsj, (size_t)n * sizeof(int)) == 0,
        "%s: reducer modified rsj input", label);
  CHECK(reduction.kind >= NBODY_REDUCED_ZERO &&
            reduction.kind <= NBODY_REDUCED_INVALID,
        "%s: invalid reduction kind %d", label, (int)reduction.kind);
  CHECK(reduction.sign == 1 || reduction.sign == -1,
        "%s: invalid reduction sign %d", label, reduction.sign);
  if (reduction.kind >= NBODY_REDUCED_ZERO &&
      reduction.kind <= NBODY_REDUCED_INVALID) {
    kindCounts[reduction.kind]++;
  }

  if (!originalNonzero) {
    CHECK(reduction.kind == NBODY_REDUCED_ZERO,
          "%s: oracle zero but reducer kind=%d order=%d sign=%d",
          label, (int)reduction.kind, reduction.order, reduction.sign);
    CHECK(reduction.order == 0, "%s: zero reduction order=%d",
          label, reduction.order);
    return;
  }

  CHECK(reduction.kind != NBODY_REDUCED_ZERO,
        "%s: reducer returned zero for nonzero oracle", label);
  CHECK(reduction.kind != NBODY_REDUCED_INVALID,
        "%s: reducer returned invalid for valid oracle", label);

  if (reduction.kind == NBODY_REDUCED_SCALAR) {
    CHECK(reduction.order == 0, "%s: scalar order=%d", label, reduction.order);
    CHECK(originalTarget == state,
          "%s: scalar changed state 0x%llx -> 0x%llx", label,
          (unsigned long long)state, (unsigned long long)originalTarget);
    CHECK(reduction.sign == originalSign,
          "%s: scalar sign=%d oracle=%d", label,
          reduction.sign, originalSign);
    return;
  }

  CHECK(reduction.kind == NBODY_REDUCED_HOPS,
        "%s: expected hops, kind=%d", label, (int)reduction.kind);
  CHECK(reduction.order > 0 && reduction.order <= n,
        "%s: invalid hop order=%d", label, reduction.order);
  if (reduction.order <= 0 || reduction.order > n) return;

  reducedNonzero = ApplyOpsBitstring(state, reduction.order,
                                     rsiWork, rsjWork,
                                     &reducedTarget, &reducedSign);
  CHECK(reducedNonzero, "%s: reduced hops evaluate to zero", label);
  if (!reducedNonzero) return;
  CHECK(reducedTarget == originalTarget,
        "%s: target mismatch original=0x%llx reduced=0x%llx", label,
        (unsigned long long)originalTarget,
        (unsigned long long)reducedTarget);
  CHECK(reduction.sign * reducedSign == originalSign,
        "%s: sign mismatch reduction=%d reduced=%d original=%d", label,
        reduction.sign, reducedSign, originalSign);
}

static void TestInvalidArguments(void) {
  int rsi[2] = {1, 2};
  int rsj[2] = {0, 1};
  int occupation[3] = {1, 1, 0};
  int rsiBefore[2];
  int rsjBefore[2];
  int occupationBefore[3];
  int badOccupation[3] = {1, 2, 0};
  int rsiWork[2] = {12345, 12345};
  int rsjWork[2] = {-12345, -12345};
  int sharedWork[2] = {777, 777};
  int badRsi[2] = {1, 3};
  union {
    int work[4];
    NBodyReduction reduction;
  } outputAlias;
  unsigned char outputAliasBefore[sizeof(outputAlias)];
  NBodyReduction reduction;

  CHECK(ReduceNBodyTerm(0, rsi, rsj, occupation, 3,
                        rsiWork, rsjWork, 2, &reduction) != 0,
        "n<=0 was accepted");
  CHECK(ReduceNBodyTerm(2, rsi, rsj, occupation, 0,
                        rsiWork, rsjWork, 2, &reduction) != 0,
        "nOrbitals<=0 was accepted");
  CHECK(ReduceNBodyTerm(2, NULL, rsj, occupation, 3,
                        rsiWork, rsjWork, 2, &reduction) != 0,
        "NULL rsi was accepted");
  CHECK(ReduceNBodyTerm(2, rsi, NULL, occupation, 3,
                        rsiWork, rsjWork, 2, &reduction) != 0,
        "NULL rsj was accepted");
  CHECK(ReduceNBodyTerm(2, rsi, rsj, NULL, 3,
                        rsiWork, rsjWork, 2, &reduction) != 0,
        "NULL occupation was accepted");
  CHECK(ReduceNBodyTerm(2, rsi, rsj, occupation, 3,
                        NULL, rsjWork, 2, &reduction) != 0,
        "NULL rsiWork was accepted");
  CHECK(ReduceNBodyTerm(2, rsi, rsj, occupation, 3,
                        rsiWork, NULL, 2, &reduction) != 0,
        "NULL rsjWork was accepted");
  CHECK(ReduceNBodyTerm(2, rsi, rsj, occupation, 3,
                        rsiWork, rsjWork, 2, NULL) != 0,
        "NULL reduction was accepted");
  CHECK(ReduceNBodyTerm(2, badRsi, rsj, occupation, 3,
                        rsiWork, rsjWork, 2, &reduction) != 0,
        "out-of-range orbital was accepted");
  CHECK(ReduceNBodyTerm(2, rsi, rsj, badOccupation, 3,
                        rsiWork, rsjWork, 2, &reduction) != 0,
        "non-binary occupation was accepted");
  CHECK(ReduceNBodyTerm(2, rsi, rsj, occupation, 3,
                        sharedWork, sharedWork, 2, &reduction) != 0,
        "aliased work arrays were accepted");

  memcpy(rsiBefore, rsi, sizeof(rsi));
  memcpy(rsjBefore, rsj, sizeof(rsj));
  memcpy(occupationBefore, occupation, sizeof(occupation));
  CHECK(ReduceNBodyTerm(2, rsi, rsj, occupation, 3,
                        rsj, rsiWork, 2, &reduction) != 0,
        "rsiWork aliasing rsj input was accepted");
  CHECK(ReduceNBodyTerm(2, rsi, rsj, occupation, 3,
                        rsiWork, rsi, 2, &reduction) != 0,
        "rsjWork aliasing rsi input was accepted");
  CHECK(ReduceNBodyTerm(2, rsi, rsj, occupation, 3,
                        occupation, rsjWork, 2, &reduction) != 0,
        "rsiWork aliasing occupation input was accepted");
  CHECK(memcmp(rsiBefore, rsi, sizeof(rsi)) == 0
            && memcmp(rsjBefore, rsj, sizeof(rsj)) == 0
            && memcmp(occupationBefore, occupation, sizeof(occupation)) == 0,
        "rejected cross-alias call modified an input array");

  memset(&outputAlias, 0x5a, sizeof(outputAlias));
  memcpy(outputAliasBefore, &outputAlias, sizeof(outputAlias));
  CHECK(ReduceNBodyTerm(2, rsi, rsj, occupation, 3,
                        outputAlias.work, rsjWork, 2,
                        &outputAlias.reduction) != 0,
        "reduction aliasing a work array was accepted");
  CHECK(memcmp(outputAliasBefore, &outputAlias, sizeof(outputAlias)) == 0,
        "rejected reduction/work alias modified shared storage");

  rsiWork[0] = rsiWork[1] = 12345;
  rsjWork[0] = rsjWork[1] = -12345;
  CHECK(ReduceNBodyTerm(2, rsi, rsj, occupation, 3,
                        rsiWork, rsjWork, 1, &reduction) != 0,
        "undersized work arrays were accepted");
  CHECK(rsiWork[0] == 12345 && rsiWork[1] == 12345 &&
            rsjWork[0] == -12345 && rsjWork[1] == -12345,
        "undersized work arrays were modified");
  CHECK(reduction.kind == NBODY_REDUCED_INVALID &&
            reduction.order == 0 && reduction.sign == 1,
        "invalid result was not fully initialized");
}

static void TestNamedSixOrbitalCases(void) {
  {
    const int rsi[3] = {5, 2, 4};
    const int rsj[3] = {2, 1, 0};
    CheckOneReduction(UINT64_C(0x03), 3, 6, rsi, rsj,
                      "six-orbital annihilation-creation contraction");
  }
  {
    const int rsi[3] = {2, 4, 5};
    const int rsj[3] = {1, 2, 5};
    CheckOneReduction(UINT64_C(0x26), 3, 6, rsi, rsj,
                      "six-orbital creation-annihilation minus contraction");
  }
  {
    const int rsi[3] = {0, 1, 2};
    const int rsj[3] = {0, 1, 2};
    CheckOneReduction(UINT64_C(0x07), 3, 6, rsi, rsj,
                      "six-orbital complete scalar reduction");
  }
  {
    const int rsi[3] = {3, 4, 5};
    const int rsj[3] = {0, 0, 1};
    CheckOneReduction(UINT64_C(0x03), 3, 6, rsi, rsj,
                      "six-orbital repeated annihilation zero");
  }
}

static void TestExhaustiveFourOrbital(void) {
  int n;
  uint64_t state;
  unsigned long caseCount = 0;

  for (n = 1; n <= 4; n++) {
    uint64_t termCount = UINT64_C(1) << (4 * n);
    uint64_t termCode;
    for (termCode = 0; termCode < termCount; termCode++) {
      uint64_t code = termCode;
      int rsi[4];
      int rsj[4];
      int k;
      for (k = 0; k < n; k++) {
        int pair = (int)(code & UINT64_C(0x0f));
        rsi[k] = pair >> 2;
        rsj[k] = pair & 3;
        code >>= 4;
      }
      for (state = 0; state < UINT64_C(16); state++) {
        CheckOneReduction(state, n, 4, rsi, rsj, "exhaustive-four-orbital");
        caseCount++;
      }
    }
  }

  CHECK(caseCount == 1118464UL,
        "unexpected exhaustive case count %lu", caseCount);
}

static void TestSectorHelpers(void) {
  int tJAllowed[4] = {1, 0, 0, 1};
  int tJRejected[4] = {1, 0, 1, 0};
  int tJPreAllowed[4] = {1, 0, 0, 0};
  int tJPreRejected[4] = {1, 0, 0, 1};
  int doublonAllowed[4] = {1, 0, 1, 0};
  int doublonRejected[4] = {1, 0, 0, 1};
  int unrestricted[4] = {1, 1, 1, 1};
  int before[4];

  Nsite = 2;
  NExUpdatePath = 4;
  CHECK(IsSectorStateAllowed(tJAllowed) == 1, "t-J allowed state rejected");
  CHECK(IsSectorStateAllowed(tJRejected) == 0, "t-J doublon accepted");
  CHECK(IsSectorPreserved_1hopPre(1, 0, 0, tJPreAllowed) == 1,
        "t-J allowed pre-hop rejected");
  CHECK(IsSectorPreserved_1hopPre(1, 0, 0, tJPreRejected) == 0,
        "t-J doublon-producing pre-hop accepted");

  NExUpdatePath = 5;
  CHECK(IsSectorStateAllowed(tJAllowed) == 1, "t-J path 5 allowed state rejected");
  CHECK(IsSectorStateAllowed(tJRejected) == 0, "t-J path 5 doublon accepted");

  NExUpdatePath = 6;
  CHECK(IsSectorStateAllowed(doublonAllowed) == 1,
        "doublon-only allowed state rejected");
  CHECK(IsSectorStateAllowed(doublonRejected) == 0,
        "doublon-only single occupancy accepted");

  NExUpdatePath = 0;
  CHECK(IsSectorStateAllowed(unrestricted) == 1,
        "unrestricted sector rejected");

  NExUpdatePath = 4;
  memcpy(before, tJAllowed, sizeof(before));
  (void)WouldPreserve_2hopPre(1, 0, 2, 3, tJAllowed);
  CHECK(memcmp(before, tJAllowed, sizeof(before)) == 0,
        "WouldPreserve_2hopPre did not restore occupation");
}

static BFNBodyScratchDimensions ScratchDimensions(int useFsz) {
  BFNBodyScratchDimensions dimensions;
  memset(&dimensions, 0, sizeof(dimensions));
  dimensions.maxOrder = 4;
  dimensions.useFsz = useFsz;
  dimensions.nsize = 4;
  dimensions.nsite2 = 8;
  dimensions.nproj = 3;
  dimensions.nsite = 4;
  dimensions.nrange = 2;
  dimensions.nqpFull = 2;
  dimensions.lapackLWork = 64;
  if (useFsz) {
    dimensions.bfFszGreenBufferCount = 96;
    dimensions.pfUpdateIntCount = 24;
    dimensions.pfUpdateDoubleCount = 80;
    dimensions.bfHopIntSize = 40;
  }
  return dimensions;
}

static size_t SumIntSlices(const BFNBodyScratchSizes *sizes) {
  return sizes->rsiCount + sizes->rsjCount + sizes->movedCount
      + sizes->eleIdxCount + sizes->eleCfgCount + sizes->eleNumCount
      + sizes->eleSpnCount + sizes->projCntCount + sizes->projBFCntCount
      + sizes->pfIWorkCount + sizes->affectedCount
      + sizes->hopIntWorkCount;
}

static size_t SumComplexSlices(const BFNBodyScratchSizes *sizes) {
  return sizes->greenBufferCount + sizes->slaterCount
      + sizes->candidatePfCount + sizes->pfBufMCount
      + sizes->pfWorkCount;
}

static void CheckScratchPointerLayout(const BFNBodyScratchSizes *sizes,
                                      const BFNBodyScratch *scratch,
                                      int *intBase,
                                      double complex *complexBase,
                                      double *doubleBase) {
  size_t offset = 0;
#define CHECK_INT_SLICE(field, countField)                                      \
  do {                                                                          \
    if (sizes->countField == 0) {                                                \
      CHECK(scratch->field == NULL, #field " zero slice is not NULL");          \
    } else {                                                                     \
      CHECK(scratch->field == intBase + offset, #field " offset mismatch");     \
    }                                                                            \
    offset += sizes->countField;                                                 \
  } while (0)
  CHECK_INT_SLICE(rsi, rsiCount);
  CHECK_INT_SLICE(rsj, rsjCount);
  CHECK_INT_SLICE(moved, movedCount);
  CHECK_INT_SLICE(eleIdx, eleIdxCount);
  CHECK_INT_SLICE(eleCfg, eleCfgCount);
  CHECK_INT_SLICE(eleNum, eleNumCount);
  CHECK_INT_SLICE(eleSpn, eleSpnCount);
  CHECK_INT_SLICE(projCnt, projCntCount);
  CHECK_INT_SLICE(projBFCnt, projBFCntCount);
  CHECK_INT_SLICE(pfIWork, pfIWorkCount);
  CHECK_INT_SLICE(affected, affectedCount);
  CHECK_INT_SLICE(hopIntWork, hopIntWorkCount);
#undef CHECK_INT_SLICE
  CHECK(offset == sizes->intCount, "integer slice end mismatch");

  offset = 0;
#define CHECK_COMPLEX_SLICE(field, countField)                                  \
  do {                                                                          \
    CHECK(scratch->field == complexBase + offset, #field " offset mismatch");   \
    offset += sizes->countField;                                                 \
  } while (0)
  CHECK_COMPLEX_SLICE(greenBuffer, greenBufferCount);
  CHECK_COMPLEX_SLICE(slater, slaterCount);
  CHECK_COMPLEX_SLICE(candidatePf, candidatePfCount);
  CHECK_COMPLEX_SLICE(pfBufM, pfBufMCount);
  CHECK_COMPLEX_SLICE(pfWork, pfWorkCount);
#undef CHECK_COMPLEX_SLICE
  CHECK(offset == sizes->complexCount, "complex slice end mismatch");

  CHECK(scratch->pfRWork == doubleBase, "pfRWork offset mismatch");
  CHECK(sizes->pfRWorkCount == sizes->doubleCount,
        "double slice total mismatch");
}

static void CheckExactScratchBinding(const BFNBodyScratchSizes *sizes,
                                     const char *label) {
  BFNBodyScratch scratch;
  int bindStatus;
  int *intBase = (int *)malloc(sizes->intCount * sizeof(int));
  double complex *complexBase = (double complex *)malloc(
      sizes->complexCount * sizeof(double complex));
  double *doubleBase = (double *)malloc(sizes->doubleCount * sizeof(double));

  CHECK(intBase != NULL && complexBase != NULL && doubleBase != NULL,
        "%s binding fixture allocation failed", label);
  if (intBase == NULL || complexBase == NULL || doubleBase == NULL) {
    free(intBase);
    free(complexBase);
    free(doubleBase);
    return;
  }

  memset(&scratch, 0xa5, sizeof(scratch));
  bindStatus = BindBFNBodyScratch(sizes, intBase, sizes->intCount,
                                  complexBase, sizes->complexCount,
                                  doubleBase, sizes->doubleCount, &scratch);
  CHECK(bindStatus == BF_NBODY_OK, "%s exact-size binding failed", label);
  if (bindStatus == BF_NBODY_OK) {
    CHECK(scratch.maxOrder == sizes->maxOrder
              && scratch.useFsz == sizes->useFsz
              && memcmp(&scratch.sizes, sizes, sizeof(*sizes)) == 0,
          "%s binding did not retain authoritative metadata", label);
    CheckScratchPointerLayout(sizes, &scratch, intBase, complexBase,
                              doubleBase);
  }

  free(intBase);
  free(complexBase);
  free(doubleBase);
}

static void TestScratchSizes(void) {
  BFNBodyScratchDimensions dimensions = ScratchDimensions(0);
  BFNBodyScratchDimensions invalid[10];
  BFNBodyScratchDimensions overflow;
  BFNBodyScratchSizes nonFsz;
  BFNBodyScratchSizes fsz;
  BFNBodyScratchSizes zeroNproj;
  BFNBodyScratchSizes sizes;
  BFNBodyScratch scratch;
  BFNBodyScratch before;
  int *intBase;
  double complex *complexBase;
  double *doubleBase;
  int i;

  for (i = 0; i < 10; i++) invalid[i] = dimensions;
  invalid[0].maxOrder = 0;
  invalid[1].nsize = 0;
  invalid[2].nsite2 = 0;
  invalid[3].nproj = -1;
  invalid[4].nsite = 0;
  invalid[5].nrange = 0;
  invalid[6].nqpFull = 0;
  invalid[7].lapackLWork = 0;
  invalid[8].nsize = -1;
  invalid[9].nrange = -1;
  for (i = 0; i < 10; i++) {
    memset(&sizes, 0x5a, sizeof(sizes));
    CHECK(ComputeBFNBodyScratchSizes(&invalid[i], &sizes)
              == BF_NBODY_WORKSPACE_ERROR,
          "invalid scratch dimensions case %d accepted", i);
  }

  dimensions = ScratchDimensions(0);
  dimensions.nproj = 0;
  CHECK(ComputeBFNBodyScratchSizes(&dimensions, &zeroNproj) == BF_NBODY_OK
            && zeroNproj.projCntCount == 0,
        "zero-length optional NProj slice rejected");
  dimensions = ScratchDimensions(0);
  dimensions.nsite2 = 6;
  CHECK(ComputeBFNBodyScratchSizes(&dimensions, &sizes)
            == BF_NBODY_WORKSPACE_ERROR,
        "nsite2 != 2*nsite accepted");
  dimensions = ScratchDimensions(0);
  dimensions.useFsz = 2;
  CHECK(ComputeBFNBodyScratchSizes(&dimensions, &sizes)
            == BF_NBODY_WORKSPACE_ERROR,
        "invalid useFsz accepted");

  overflow = ScratchDimensions(0);
  overflow.nsite = INT_MAX/2;
  overflow.nsite2 = INT_MAX-1;
  overflow.nrange = INT_MAX;
  CHECK(ComputeBFNBodyScratchSizes(&overflow, &sizes)
            == BF_NBODY_WORKSPACE_ERROR,
        "overflowing scratch multiplication accepted");

  CHECK(ComputeBFNBodyScratchSizes(NULL, &sizes)
            == BF_NBODY_WORKSPACE_ERROR,
        "NULL dimensions accepted");
  CHECK(ComputeBFNBodyScratchSizes(&dimensions, NULL)
            == BF_NBODY_WORKSPACE_ERROR,
        "NULL sizes output accepted");

  dimensions = ScratchDimensions(0);
  CHECK(ComputeBFNBodyScratchSizes(&dimensions, &nonFsz) == BF_NBODY_OK,
        "non-FSZ scratch sizing failed");
  dimensions = ScratchDimensions(1);
  CHECK(ComputeBFNBodyScratchSizes(&dimensions, &fsz) == BF_NBODY_OK,
        "FSZ scratch sizing failed");

  CHECK(nonFsz.intCount == SumIntSlices(&nonFsz),
        "non-FSZ integer total is not the slice sum");
  CHECK(nonFsz.complexCount == SumComplexSlices(&nonFsz),
        "non-FSZ complex total is not the slice sum");
  CHECK(nonFsz.doubleCount == nonFsz.pfRWorkCount,
        "non-FSZ double total is not the slice sum");
  CHECK(fsz.intCount == SumIntSlices(&fsz),
        "FSZ integer total is not the slice sum");
  CHECK(fsz.complexCount == SumComplexSlices(&fsz),
        "FSZ complex total is not the slice sum");
  CHECK(fsz.doubleCount == fsz.pfRWorkCount,
        "FSZ double total is not the slice sum");
  CHECK(nonFsz.affectedCount == 0 && nonFsz.hopIntWorkCount == 0,
        "non-FSZ scratch contains FSZ-only slices");
  CHECK(fsz.affectedCount == (size_t)dimensions.nsize
            && fsz.hopIntWorkCount == (size_t)dimensions.bfHopIntSize,
        "FSZ-only scratch slices have wrong sizes");
  CHECK(nonFsz.greenBufferCount == 10
            && fsz.greenBufferCount
                == dimensions.bfFszGreenBufferCount
            && nonFsz.pfIWorkCount == 4
            && fsz.pfIWorkCount == dimensions.pfUpdateIntCount
            && nonFsz.pfRWorkCount == 64
            && fsz.pfRWorkCount == dimensions.pfUpdateDoubleCount,
        "FSZ update-work maxima have wrong sizes");
  CHECK(nonFsz.rsiCount == fsz.rsiCount
            && nonFsz.rsjCount == fsz.rsjCount
            && nonFsz.movedCount == fsz.movedCount
            && nonFsz.eleIdxCount == fsz.eleIdxCount
            && nonFsz.eleCfgCount == fsz.eleCfgCount
            && nonFsz.eleNumCount == fsz.eleNumCount
            && nonFsz.eleSpnCount == fsz.eleSpnCount
            && nonFsz.projCntCount == fsz.projCntCount
            && nonFsz.projBFCntCount == fsz.projBFCntCount
            && nonFsz.slaterCount == fsz.slaterCount
            && nonFsz.candidatePfCount == fsz.candidatePfCount
            && nonFsz.pfBufMCount == fsz.pfBufMCount
            && nonFsz.pfWorkCount == fsz.pfWorkCount,
        "FSZ changed a non-update scratch slice");

  sizes = fsz;
  intBase = (int *)malloc(sizes.intCount*sizeof(int));
  complexBase = (double complex *)malloc(
      sizes.complexCount*sizeof(double complex));
  doubleBase = (double *)malloc(sizes.doubleCount*sizeof(double));
  CHECK(intBase != NULL && complexBase != NULL && doubleBase != NULL,
        "scratch fixture allocation failed");
  if (intBase == NULL || complexBase == NULL || doubleBase == NULL) {
    free(intBase);
    free(complexBase);
    free(doubleBase);
    return;
  }

  memset(&scratch, 0xa5, sizeof(scratch));
  before = scratch;
  CHECK(BindBFNBodyScratch(&sizes, intBase, sizes.intCount-1,
                           complexBase, sizes.complexCount,
                           doubleBase, sizes.doubleCount, &scratch)
            == BF_NBODY_WORKSPACE_ERROR,
        "undersized integer base accepted");
  CHECK(memcmp(&scratch, &before, sizeof(scratch)) == 0,
        "failed integer binding modified scratch");
  CHECK(BindBFNBodyScratch(&sizes, intBase, sizes.intCount,
                           complexBase, sizes.complexCount-1,
                           doubleBase, sizes.doubleCount, &scratch)
            == BF_NBODY_WORKSPACE_ERROR,
        "undersized complex base accepted");
  CHECK(memcmp(&scratch, &before, sizeof(scratch)) == 0,
        "failed complex binding modified scratch");
  CHECK(BindBFNBodyScratch(&sizes, intBase, sizes.intCount,
                           complexBase, sizes.complexCount,
                           doubleBase, sizes.doubleCount-1, &scratch)
            == BF_NBODY_WORKSPACE_ERROR,
        "undersized double base accepted");
  CHECK(memcmp(&scratch, &before, sizeof(scratch)) == 0,
        "failed double binding modified scratch");

  CHECK(BindBFNBodyScratch(&sizes, intBase, sizes.intCount,
                           complexBase, sizes.complexCount,
                           doubleBase, sizes.doubleCount, &scratch)
            == BF_NBODY_OK,
        "exact-size scratch binding failed");
  CHECK(scratch.maxOrder == sizes.maxOrder
            && scratch.useFsz == sizes.useFsz,
        "scratch metadata did not come from sizes");
  CHECK(memcmp(&scratch.sizes, &sizes, sizeof(sizes)) == 0,
        "bound scratch did not retain authoritative sizes");
  CheckScratchPointerLayout(&sizes, &scratch, intBase, complexBase,
                            doubleBase);

  free(intBase);
  free(complexBase);
  free(doubleBase);

  CheckExactScratchBinding(&nonFsz, "non-FSZ");
  CheckExactScratchBinding(&zeroNproj, "zero-NProj");
}

int main(void) {
  TestInvalidArguments();
  TestNamedSixOrbitalCases();
  TestExhaustiveFourOrbital();
  TestSectorHelpers();
  TestScratchSizes();

  CHECK(kindCounts[NBODY_REDUCED_ZERO] > 0, "zero class was not covered");
  CHECK(kindCounts[NBODY_REDUCED_SCALAR] > 0, "scalar class was not covered");
  CHECK(kindCounts[NBODY_REDUCED_HOPS] > 0, "hops class was not covered");

  if (failures != 0) {
    fprintf(stderr, "NBodyOperator_Unit: %d failure(s)\n", failures);
    return 1;
  }

  printf("NBodyOperator_Unit PASS zero=%lu scalar=%lu hops=%lu\n",
         kindCounts[NBODY_REDUCED_ZERO],
         kindCounts[NBODY_REDUCED_SCALAR],
         kindCounts[NBODY_REDUCED_HOPS]);
  return 0;
}
