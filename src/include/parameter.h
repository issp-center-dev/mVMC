#ifndef _VMC_INCLUDE_PARAMETER
#define _VMC_INCLUDE_PARAMETER

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "sfmt.h"
#include "global.h"

#ifdef _mpi_use
#include <mpi.h>
#else
typedef int MPI_Comm;
MPI_Comm MPI_COMM_WORLD=0;
inline void MPI_Init(int argc, char* argv[]) {return;}
inline void MPI_Finalize() {return;}
inline void MPI_Abort(MPI_Comm comm, int errorcode) {exit(errorcode); return;}
inline void MPI_Barrier(MPI_Comm comm) {return;}
inline void MPI_Comm_size(MPI_Comm comm, int *size) {*size = 1; return;}
inline void MPI_Comm_rank(MPI_Comm comm, int *rank) {*rank = 0; return;}
#endif /* _mpi_use */

void InitParameter();
int ReadInitParameter(char *initFile);
void SyncModifiedParameter(MPI_Comm comm);
void SetFlagShift();

void shiftGJ();
double shiftDH2();
double shiftDH4();

#endif
