#ifndef _AVERAGE
#define _AVERAGE
#include <mpi.h>
void WeightAverageWE(MPI_Comm comm);
void WeightAverageSROpt(MPI_Comm comm);
void WeightAverageSROpt_real(MPI_Comm comm);
void WeightAverageGreenFunc(MPI_Comm comm);
#endif
