#ifndef _QP
#define _QP

void InitQPWeight();

double complex CalculateLogIP_fcmp(double complex * const pfM, const int qpStart, const int qpEnd, MPI_Comm comm);
double complex CalculateIP_fcmp(double complex * const pfM, const int qpStart, const int qpEnd, MPI_Comm comm);
void UpdateQPWeight();

#endif
