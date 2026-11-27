#ifndef _AVERAGE
#define _AVERAGE
#include <stddef.h>
#include "mvmc_mpi.h"
void WeightAverageWE(MPI_Comm comm);
void WeightAverageSROpt(MPI_Comm comm);
void WeightAverageSROpt_real(MPI_Comm comm);
/* Returns the bin's Lanczos support sample count on rank zero. */
size_t WeightAverageGreenFunc(MPI_Comm comm);
void WeightAverageNearZeroGuide(MPI_Comm comm);
#endif
