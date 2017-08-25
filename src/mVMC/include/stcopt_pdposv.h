#ifndef _STCOPT_DPOSV
#ifndef _STCOPT_PDPOSV
#define _STCOPT_PDPOSV

int StochasticOptDiag(MPI_Comm comm);
int stcOptMainDiag(double *const r, int const nSmat, int *const smatToParaIdx,
               MPI_Comm comm, int const optNum);

#endif
#endif

