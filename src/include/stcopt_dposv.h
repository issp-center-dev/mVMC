#ifndef _STCOPT_PDPOSV
#ifndef _STCOPT_DPOSV
#define _STCOPT_DPOSV
#include <mpi.h>

int StochasticOpt(MPI_Comm comm);
void stcOptInit(double *const s, double *const g, const int nSmat, int *const smatToParaIdx);
int stcOptMain(double *const s, double *const g, const int nSmat);


#endif
#endif
