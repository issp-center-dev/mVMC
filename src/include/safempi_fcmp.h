/*-------------------------------------------------------------
 * Variational Monte Carlo
 * safe version of mpi functions
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#pragma once
#include <complex.h>
#define D_MpiSendMax  1048576 /* 2^20 */
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

void SafeMpiReduce_fcmp(double complex *send, double complex *recv, int nData, MPI_Comm comm);
void SafeMpiAllReduce_fcmp(double complex *send, double complex *recv, int nData, MPI_Comm comm);
void SafeMpiBcast_fcmp(double complex *buff, int nData, MPI_Comm comm);
