#ifndef _QP
#define _QP

void InitQPWeight();

double complex CalculateLogIP_fcmp(double complex * const pfM, const int qpStart, const int qpEnd, MPI_Comm comm);
double complex CalculateIP_fcmp(double complex * const pfM, const int qpStart, const int qpEnd, MPI_Comm comm);
/* newN is the candidate particle count, including pair-add/remove changes. */
double complex CalculateIPGC_fcmp(double complex *pfM, int newN,
                                  int qpStart, int qpEnd, MPI_Comm comm);
double complex CalculateLogIPGC_fcmp(double complex *pfM, int newN,
                                     int qpStart, int qpEnd, MPI_Comm comm);
void UpdateQPWeight();
#ifdef MVMC_ENABLE_FAULT_INJECTION
/* Testing only: rank 0 writes the spin quadrature coefficients and full qp
 * weights to MVMC_GC_SPIN_QP_DUMP (lossless %a); collective over comm. */
void GCTestDumpSpinQP(MPI_Comm comm);
#endif

#endif
