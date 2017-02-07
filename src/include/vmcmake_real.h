#ifndef _VMCMAKE_REAL
#define _VMCMAKE_REAL
#include <complex.h>
#include <mpi.h>

void VMCMakeSample_real(MPI_Comm comm);
void VMC_BF_MakeSample_real(MPI_Comm comm);

void copyFromBurnSampleBF(int *eleIdx);
void copyToBurnSampleBF(const int *eleIdx);
void saveEleConfigBF(const int sample, const double logIp,
                     const int *eleIdx, const int *eleCfg, const int *eleNum, const int *eleProjCnt, const int *eleProjBFCnt);


int makeInitialSample_real(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,
                           const int qpStart, const int qpEnd, MPI_Comm comm);

int makeInitialSampleBF_real(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt, int *eleProjBFCnt,
                             const int qpStart, const int qpEnd, MPI_Comm comm);

#endif


