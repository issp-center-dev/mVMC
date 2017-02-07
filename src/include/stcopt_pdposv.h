#ifndef _STCOPT_PDPOSV
#define _STCOPT_PDPOSV

int StochasticOpt(MPI_Comm comm);
int stcOptMain(double *r, const int nSmat, const int *smatToParaIdx, MPI_Comm comm);
int StochasticOptDiag(MPI_Comm comm);
int stcOptMainDiag(double *const r, int const nSmat, int *const smatToParaIdx,
               MPI_Comm comm, int const optNum);

#endif
