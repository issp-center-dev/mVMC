#ifndef _VMCMAKE_REAL
#define _VMCMAKE_REAL
#include <complex.h>
#include <mpi.h>

void VMCMakeSample_real(MPI_Comm comm);
void VMC_BF_MakeSample_real(MPI_Comm comm);

int makeInitialSample_real(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,
                           const int qpStart, const int qpEnd, MPI_Comm comm);

int makeInitialSampleBF_real(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt, int *eleProjBFCnt,
                             const int qpStart, const int qpEnd, MPI_Comm comm);

#endif


