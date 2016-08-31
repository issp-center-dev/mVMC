/*-------------------------------------------------------------
 * Variational Monte Carlo
 * safe version of mpi functions
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#ifndef _VMC_INCLUDE_SAFEMPI
#define _VMC_INCLUDE_SAFEMPI

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

#define D_MpiSendMax  1048576 /* 2^20 */

void SafeMpiReduce(double *send, double *recv, int nData, MPI_Comm comm);
void SafeMpiAllReduce(double *send, double *recv, int nData, MPI_Comm comm);
void SafeMpiBcast(double *buff, int nData, MPI_Comm comm);
void SafeMpiBcastInt(int *buff, int nData, MPI_Comm comm);

#endif
