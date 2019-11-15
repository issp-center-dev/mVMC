#ifndef _VMCMAKE_FSZ_REAL
#define _VMCMAKE_FSZ_REAL
#include <mpi.h>

void VMCMakeSample_fsz_real(MPI_Comm comm);

int makeInitialSample_fsz_real(int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,int *eleSpn,
                      const int qpStart, const int qpEnd, MPI_Comm comm);

#endif


