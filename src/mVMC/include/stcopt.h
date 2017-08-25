#ifndef _STCOPT_HEADER
#define _STCOPT_HEADER
#include <mpi.h>

int StochasticOpt(MPI_Comm comm);
void stcOptInit(double *const s, double *const g, const int nSmat, const int *const smatToParaIdx);
int stcOptMain(double *g, const int nSmat, const int *smatToParaIdx, MPI_Comm);

#endif
