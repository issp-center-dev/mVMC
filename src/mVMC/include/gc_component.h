#ifndef _GC_COMPONENT_H
#define _GC_COMPONENT_H

#include <complex.h>
#ifdef _mpi_use
#include <mpi.h>
typedef MPI_Comm GCComponentComm;
#else
typedef int GCComponentComm;
#endif

#define GC_QP_INVERSE 0
#define GC_QP_PFONLY 1

typedef struct GCComponentCache GCComponentCache;

#ifdef MVMC_ENABLE_FAULT_INJECTION
enum GCComponentOmpPath {
  GC_OMP_GREEN,
  GC_OMP_COFACTOR,
  GC_OMP_PATHS
};

typedef struct {
  int teamMin;
  int teamMax;
  int workers;
  long long workCount;
} GCComponentOmpAudit;
#endif

typedef struct {
  long long rebuild;
  long long toPfonly;
  long long recover;
  long long probeRebuild;
  long long debugNearSingular;
  long long directValue;
  long long cofactorGradient;
  long long compositeZero;
} GCComponentCounters;

int GCComponentInit(void);
void GCComponentRelease(void);
void GCComponentFinalize(GCComponentComm comm);
GCComponentCache *GCComponentThreadCache(void);
GCComponentCounters *GCComponentCountersOf(GCComponentCache *cache);
void GCComponentCountersTotal(GCComponentCounters *total);
const char *GCComponentStateName(int state);

int GCComponentRebuild(int qpGlobal, int qpLocal, int ncur,
                       const int *eleIdx, GCComponentCache *cache);
int GCComponentValue(int qpGlobal, int ncur, const int *eleIdx,
                     GCComponentCache *cache, double complex *pf);
int GCComponentValueAppended(int qpGlobal, int ncur, const int *eleIdx,
                             int rsa, int rsb, GCComponentCache *cache,
                             double complex *pf);
int GCComponentValueExcluding(int qpGlobal, int ncur, const int *eleIdx,
                              int pos0, int pos1, GCComponentCache *cache,
                              double complex *pf);
int GCComponentProbe(int qpGlobal, int qpLocal, int ncur,
                     const int *eleIdx, GCComponentCache *cache);
int GCComponentGradient(int qpGlobal, int qpLocal, int ncur,
                        const int *eleIdx, double complex *gradient,
                        GCComponentCache *cache);
#ifdef MVMC_ENABLE_FAULT_INJECTION
void GCComponentRecordOmpWork(enum GCComponentOmpPath path);
void GCComponentOmpAuditTotal(enum GCComponentOmpPath path,
                              GCComponentOmpAudit *out);
#endif

#endif
