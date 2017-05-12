#ifndef _QP_REAL
#define _QP_REAL

double CalculateLogIP_real(double* const pfM, const int qpStart, const int qpEnd, MPI_Comm comm);
double CalculateIP_real(double* const pfM, const int qpStart, const int qpEnd, MPI_Comm comm);

#endif
