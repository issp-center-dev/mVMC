#include <complex.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

#include "include/absolute_pfaffian.h"
#include "include/gc_component.h"
#include "include/gc_size.h"
#include "include/global.h"
#include "include/matrix_gc.h"

#ifndef _SRC_GC_COMPONENT
#define _SRC_GC_COMPONENT

#define GC_COMPONENT_SLOTS 4
/* Same factor as effective_residual_tolerance() in absolute_pfaffian.c. */
#define GC_COMPONENT_RESIDUAL_FACTOR 256.0

typedef struct {
  int n;
  unsigned long long lastUse;
  MVMCAbsolutePfaffianComplexWorkspace *full;
  MVMCAbsolutePfaffianComplexValueWorkspace *value;
} GCComponentSlot;

struct GCComponentCache {
  GCComponentSlot slots[GC_COMPONENT_SLOTS];
  unsigned long long clock;
  unsigned int probeSeed;
  int capacity;
  int *idx;
  double complex *x;
  double complex *inverse;
  double complex *minor;
  double complex *vector;
  GCComponentCounters counters;
};

static GCComponentCache **gcCaches = NULL;
static int gcCacheCount = 0;
static int gcCapacity = 0;
static int gcQpCapacity = 0;

const char *GCComponentStateName(const int state) {
  return state == GC_QP_INVERSE ? "INVERSE"
         : state == GC_QP_PFONLY ? "PFONLY" : "unknown";
}

static int GCComponentStatusFromBackend(const MVMCPfaffianStatus status) {
  switch (status) {
    case MVMC_PFAFFIAN_STATUS_OK:
      return GC_MALL_OK;
    case MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE:
      return GC_MALL_FP_MODE;
    case MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE:
      return GC_MALL_ALLOCATION;
    case MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE:
      return GC_MALL_PFAFFIAN;
    default:
      return GC_MALL_INVALID_ARGUMENT;
  }
}

static void GCComponentSlotClear(GCComponentSlot *slot) {
  if (slot->full != NULL)
    mvmc_absolute_pfaffian_complex_workspace_destroy(slot->full);
  if (slot->value != NULL)
    mvmc_absolute_pfaffian_complex_value_workspace_destroy(slot->value);
  slot->full = NULL;
  slot->value = NULL;
  slot->n = 0;
  slot->lastUse = 0;
}

static int GCComponentSlotFor(GCComponentCache *cache, const int n,
                              const int needFull, GCComponentSlot **out) {
  GCComponentSlot *slot = NULL;
  MVMCPfaffianStatus status;
  int i;
  *out = NULL;
  for (i = 0; i < GC_COMPONENT_SLOTS && slot == NULL; i++)
    if (cache->slots[i].n == n) slot = &cache->slots[i];
  for (i = 0; i < GC_COMPONENT_SLOTS && slot == NULL; i++)
    if (cache->slots[i].n == 0) slot = &cache->slots[i];
  if (slot == NULL) {
    slot = &cache->slots[0];
    for (i = 1; i < GC_COMPONENT_SLOTS; i++)
      if (cache->slots[i].lastUse < slot->lastUse) slot = &cache->slots[i];
    GCComponentSlotClear(slot);
  }
  if (slot->n == 0) {
    status = mvmc_absolute_pfaffian_complex_value_workspace_create(
        n, &slot->value);
    if (status != MVMC_PFAFFIAN_STATUS_OK) {
      GCComponentSlotClear(slot);
      return GCComponentStatusFromBackend(status);
    }
    slot->n = n;
  }
  if (needFull && slot->full == NULL) {
    status = mvmc_absolute_pfaffian_complex_workspace_create(n, &slot->full);
    if (status != MVMC_PFAFFIAN_STATUS_OK) {
      GCComponentSlotClear(slot);
      return GCComponentStatusFromBackend(status);
    }
  }
  slot->lastUse = ++cache->clock;
  *out = slot;
  return GC_MALL_OK;
}

void GCComponentRelease(void) {
  int t;
  if (gcCaches != NULL) {
    for (t = 0; t < gcCacheCount; t++) {
      int s;
      for (s = 0; s < GC_COMPONENT_SLOTS; s++)
        GCComponentSlotClear(&gcCaches[t]->slots[s]);
      free(gcCaches[t]->vector);
      free(gcCaches[t]->minor);
      free(gcCaches[t]->inverse);
      free(gcCaches[t]->x);
      free(gcCaches[t]->idx);
      free(gcCaches[t]);
    }
    free(gcCaches);
  }
  gcCaches = NULL;
  gcCacheCount = 0;
  gcCapacity = 0;
  gcQpCapacity = 0;
  free(GCQpState);
  GCQpState = NULL;
}

int GCComponentInit(void) {
  size_t matrixCount;
  int t;
  if (NThread <= 0 || NsizeMax <= 0 || NQPFull <= 0)
    return GC_MALL_INVALID_ARGUMENT;
  if (gcCaches != NULL && gcCacheCount == NThread &&
      gcCapacity == NsizeMax && gcQpCapacity == NQPFull)
    return GC_MALL_OK;
  matrixCount = GCCheckedMulSize((size_t)NsizeMax, (size_t)NsizeMax);
  GCComponentRelease();
  GCQpState = (int *)GCCheckedMalloc(
      GCCheckedMulSize(sizeof(int), (size_t)NQPFull), "GCQpState");
  for (t = 0; t < NQPFull; t++) GCQpState[t] = GC_QP_INVERSE;
  gcCaches = (GCComponentCache **)GCCheckedMalloc(
      GCCheckedMulSize(sizeof(*gcCaches), (size_t)NThread),
      "GC component caches");
  for (t = 0; t < NThread; t++) {
    GCComponentCache *cache = (GCComponentCache *)GCCheckedMalloc(
        sizeof(*cache), "GC component cache");
    memset(cache, 0, sizeof(*cache));
    cache->capacity = NsizeMax;
    cache->probeSeed = 0x9e3779b9u + (unsigned int)t;
    cache->idx = (int *)GCCheckedMalloc(
        GCCheckedMulSize(sizeof(int), (size_t)NsizeMax),
        "GC component index buffer");
    cache->x = (double complex *)GCCheckedMalloc(
        GCCheckedMulSize(sizeof(double complex), matrixCount),
        "GC component X");
    cache->inverse = (double complex *)GCCheckedMalloc(
        GCCheckedMulSize(sizeof(double complex), matrixCount),
        "GC component inverse");
    cache->minor = (double complex *)GCCheckedMalloc(
        GCCheckedMulSize(sizeof(double complex), matrixCount),
        "GC component minor");
    cache->vector = (double complex *)GCCheckedMalloc(
        GCCheckedMulSize(sizeof(double complex),
                         GCCheckedMulSize(2, (size_t)NsizeMax)),
        "GC component probe vectors");
    gcCaches[t] = cache;
  }
  gcCacheCount = NThread;
  gcCapacity = NsizeMax;
  gcQpCapacity = NQPFull;
  return GC_MALL_OK;
}

GCComponentCache *GCComponentThreadCache(void) {
  const int thread = omp_get_thread_num();
  if (gcCaches == NULL || thread < 0 || thread >= gcCacheCount) {
    fprintf(stderr,
            "Error: GC component cache is not initialized for thread %d.\n",
            thread);
    exit(EXIT_FAILURE);
  }
  if (gcCapacity < NsizeMax || gcQpCapacity < NQPFull) {
    fprintf(stderr,
            "Error: GC component cache capacity (NsizeMax=%d, NQPFull=%d) "
            "is smaller than the current dimensions (%d, %d); call "
            "GCComponentInit after changing them.\n",
            gcCapacity, gcQpCapacity, NsizeMax, NQPFull);
    exit(EXIT_FAILURE);
  }
  return gcCaches[thread];
}

GCComponentCounters *GCComponentCountersOf(GCComponentCache *cache) {
  return &cache->counters;
}

void GCComponentCountersTotal(GCComponentCounters *total) {
  int t;
  memset(total, 0, sizeof(*total));
  for (t = 0; t < gcCacheCount; t++) {
    const GCComponentCounters *c = &gcCaches[t]->counters;
    total->rebuild += c->rebuild;
    total->toPfonly += c->toPfonly;
    total->recover += c->recover;
    total->probeRebuild += c->probeRebuild;
    total->debugNearSingular += c->debugNearSingular;
    total->directValue += c->directValue;
    total->cofactorGradient += c->cofactorGradient;
    total->compositeZero += c->compositeZero;
  }
}

void GCComponentFinalize(const GCComponentComm comm) {
  GCComponentCounters total;
  long long local[7], global[7];
  int rank = 0;
  if (gcCaches == NULL) return;
  GCComponentCountersTotal(&total);
  local[0] = total.rebuild;
  local[1] = total.toPfonly;
  local[2] = total.recover;
  local[3] = total.probeRebuild;
  local[4] = total.debugNearSingular;
  local[5] = total.directValue;
  local[6] = total.cofactorGradient;
#ifdef _mpi_use
  MPI_Reduce(local, global, 7, MPI_LONG_LONG, MPI_SUM, 0, comm);
  MPI_Comm_rank(comm, &rank);
#else
  (void)comm;
  memcpy(global, local, sizeof(global));
#endif
  if (rank == 0) {
    printf("GC component stats: thread_capacity=%d rebuild=%lld "
           "to_pfonly=%lld recover=%lld probe_rebuild=%lld "
           "debug_near_singular=%lld direct_value=%lld "
           "cofactor_gradient=%lld composite_zero=%lld\n",
           gcCacheCount, global[0], global[1], global[2], global[3],
           global[4], global[5], global[6], total.compositeZero);
  }
  GCComponentRelease();
}

static void GCComponentBuildX(const double complex *slater, const int ncur,
                              const int *idx, double complex *x) {
  int row, col;
  for (col = 0; col < ncur; col++)
    for (row = 0; row < ncur; row++)
      x[row + (size_t)col * (size_t)ncur] =
          -slater[(size_t)idx[row] * (size_t)Nsite2 + (size_t)idx[col]];
}

static int GCComponentArgumentsValid(const int qpGlobal, const int ncur,
                                     const int *eleIdx) {
  int i;
  if (qpGlobal < 0 || qpGlobal >= NQPFull || ncur < 0 ||
      ncur > NsizeMax || (ncur % 2) != 0 ||
      (ncur > 0 && eleIdx == NULL) || SlaterElm == NULL)
    return 0;
  for (i = 0; i < ncur; i++)
    if (eleIdx[i] < 0 || eleIdx[i] >= Nsite2) return 0;
  return 1;
}

int GCComponentRebuild(const int qpGlobal, const int qpLocal,
                       const int ncur, const int *eleIdx,
                       GCComponentCache *cache) {
  const size_t slaterStride =
      GCCheckedMulSize((size_t)Nsite2, (size_t)Nsite2);
  const size_t invStride =
      GCCheckedMulSize((size_t)NsizeMax, (size_t)NsizeMax);
  MVMCAbsolutePfaffianResult result;
  MVMCPfaffianStatus status;
  GCComponentSlot *slot;
  double complex *invM;
  int row, col, slotStatus;
  if (cache == NULL ||
      !GCComponentArgumentsValid(qpGlobal, ncur, eleIdx) ||
      qpLocal < 0 || qpLocal >= NQPFull || GCQpState == NULL ||
      InvM == NULL || PfM == NULL)
    return GC_MALL_INVALID_ARGUMENT;
  cache->counters.rebuild++;
  if (ncur == 0) {
    PfM[qpLocal] = 1.0;
    GCQpState[qpLocal] = GC_QP_INVERSE;
    return GC_MALL_OK;
  }
  slotStatus = GCComponentSlotFor(cache, ncur, 1, &slot);
  if (slotStatus != GC_MALL_OK) return slotStatus;
  GCComponentBuildX(SlaterElm + (size_t)qpGlobal * slaterStride, ncur,
                    eleIdx, cache->x);
  status = mvmc_absolute_pfaffian_complex_with_workspace(
      slot->full, cache->x, ncur, ncur, cache->inverse, ncur,
      (uint64_t)cache->counters.rebuild, 0.0, 0.0, &result);
  if (status != MVMC_PFAFFIAN_STATUS_OK)
    return GCComponentStatusFromBackend(status);
  if (result.state == MVMC_PFAFFIAN_NONFINITE) return GC_MALL_NONFINITE;
  if (result.state == MVMC_PFAFFIAN_INVALID)
    return GC_MALL_INVALID_ARGUMENT;
  PfM[qpLocal] =
      ((ncur / 2) % 2 == 0) ? result.pfaffian : -result.pfaffian;
  if (result.state == MVMC_PFAFFIAN_REGULAR) {
    invM = InvM + (size_t)qpLocal * invStride;
    for (row = 0; row < ncur; row++)
      for (col = 0; col < ncur; col++)
        invM[(size_t)row * (size_t)NsizeMax + (size_t)col] =
            cache->inverse[row + (size_t)col * (size_t)ncur];
    GCQpState[qpLocal] = GC_QP_INVERSE;
  } else {
    GCQpState[qpLocal] = GC_QP_PFONLY;
  }
  return GC_MALL_OK;
}

int GCComponentValue(const int qpGlobal, const int ncur,
                     const int *eleIdx, GCComponentCache *cache,
                     double complex *pf) {
  const size_t slaterStride =
      GCCheckedMulSize((size_t)Nsite2, (size_t)Nsite2);
  MVMCAbsolutePfaffianValueResult result;
  MVMCPfaffianStatus status;
  GCComponentSlot *slot;
  int slotStatus;
  if (cache == NULL || pf == NULL ||
      !GCComponentArgumentsValid(qpGlobal, ncur, eleIdx))
    return GC_MALL_INVALID_ARGUMENT;
  cache->counters.directValue++;
  if (ncur == 0) {
    *pf = 1.0;
    return GC_MALL_OK;
  }
  slotStatus = GCComponentSlotFor(cache, ncur, 0, &slot);
  if (slotStatus != GC_MALL_OK) return slotStatus;
  GCComponentBuildX(SlaterElm + (size_t)qpGlobal * slaterStride, ncur,
                    eleIdx, cache->x);
  status = mvmc_absolute_pfaffian_complex_value_with_workspace(
      slot->value, cache->x, ncur, ncur, 0.0, &result);
  if (status != MVMC_PFAFFIAN_STATUS_OK)
    return GCComponentStatusFromBackend(status);
  if (result.state == MVMC_PFAFFIAN_VALUE_NONFINITE)
    return GC_MALL_NONFINITE;
  if (result.state == MVMC_PFAFFIAN_VALUE_INVALID)
    return GC_MALL_INVALID_ARGUMENT;
  *pf = result.state == MVMC_PFAFFIAN_VALUE_SINGULAR
            ? 0.0
            : (((ncur / 2) % 2 == 0) ? result.pfaffian
                                     : -result.pfaffian);
  return GC_MALL_OK;
}

int GCComponentValueAppended(const int qpGlobal, const int ncur,
                             const int *eleIdx, const int rsa,
                             const int rsb, GCComponentCache *cache,
                             double complex *pf) {
  int i;
  if (cache == NULL ||
      !GCComponentArgumentsValid(qpGlobal, ncur, eleIdx) ||
      ncur + 2 > NsizeMax || rsa < 0 || rsa >= Nsite2 || rsb < 0 ||
      rsb >= Nsite2 || rsa == rsb)
    return GC_MALL_INVALID_ARGUMENT;
  for (i = 0; i < ncur; i++) {
    if (eleIdx[i] == rsa || eleIdx[i] == rsb)
      return GC_MALL_INVALID_ARGUMENT;
    cache->idx[i] = eleIdx[i];
  }
  cache->idx[ncur] = rsa;
  cache->idx[ncur + 1] = rsb;
  return GCComponentValue(qpGlobal, ncur + 2, cache->idx, cache, pf);
}

int GCComponentValueExcluding(const int qpGlobal, const int ncur,
                              const int *eleIdx, const int pos0,
                              const int pos1, GCComponentCache *cache,
                              double complex *pf) {
  int i, count = 0;
  if (cache == NULL ||
      !GCComponentArgumentsValid(qpGlobal, ncur, eleIdx) || ncur < 2 ||
      pos0 < 0 || pos0 >= ncur || pos1 < 0 || pos1 >= ncur ||
      pos0 == pos1)
    return GC_MALL_INVALID_ARGUMENT;
  for (i = 0; i < ncur; i++)
    if (i != pos0 && i != pos1) cache->idx[count++] = eleIdx[i];
  return GCComponentValue(qpGlobal, ncur - 2, cache->idx, cache, pf);
}

#endif
