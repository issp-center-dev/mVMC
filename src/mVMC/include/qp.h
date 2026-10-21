#ifndef _QP
#define _QP

#include "mvmc_mpi.h"

void InitQPWeight();
/* Proven support zeros only; finite component validation still precedes
 * applying this mask to any amplitude. eleNum is the candidate occupation;
 * it is required (non-NULL) only for an active finite S>0 projection, which
 * also masks configurations with fewer than 2S singles. Returns 1 for a
 * proven zero, 0 otherwise and -1 for an invalid candidate occupation. */
int GCProjectionSectorZero(int newN, const int *eleNum);
/* Active anti-parallel GC projection onto a finite total spin S>0. */
int GCFiniteSpinActive(void);

double complex CalculateLogIP_fcmp(double complex * const pfM, const int qpStart, const int qpEnd, MPI_Comm comm);
double complex CalculateIP_fcmp(double complex * const pfM, const int qpStart, const int qpEnd, MPI_Comm comm);
/* newN is the candidate particle count, including pair-add/remove changes;
 * eleNum the candidate occupation (see GCProjectionSectorZero). An invalid
 * candidate occupation gives a NaN amplitude after the collective sum. */
double complex CalculateIPGC_fcmp(double complex *pfM, int newN,
                                  const int *eleNum, int qpStart, int qpEnd,
                                  MPI_Comm comm);
double complex CalculateLogIPGC_fcmp(double complex *pfM, int newN,
                                     const int *eleNum, int qpStart,
                                     int qpEnd, MPI_Comm comm);
void UpdateQPWeight();
#ifdef MVMC_ENABLE_FAULT_INJECTION
/* Testing only: rank 0 writes the spin quadrature coefficients and full qp
 * weights to MVMC_GC_SPIN_QP_DUMP (lossless %a); collective over comm. */
void GCTestDumpSpinQP(MPI_Comm comm);
#endif

#endif
