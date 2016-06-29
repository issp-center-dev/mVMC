/*-------------------------------------------------------------
 * Variational Monte Carlo
 * safe version of mpi functions
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

#define D_MpiSendMax  1048576 /* 2^20 */

void SafeMpiReduce_fcmp(double complex *send, double complex *recv, int nData, MPI_Comm comm);
void SafeMpiAllReduce_fcmp(double complex *send, double complex *recv, int nData, MPI_Comm comm);
void SafeMpiBcast_fcmp(double complex *buff, int nData, MPI_Comm comm);

void SafeMpiReduce_fcmp(double complex *send, double complex *recv, int nData, MPI_Comm comm) {
  #ifdef _mpi_use
  int nSend = D_MpiSendMax; /* defined in global.h */
  int idx = 0;

  while(idx<nData) {
    if(idx+nSend > nData) nSend = nData-idx;
    MPI_Reduce(send+idx,recv+idx,nSend,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm);
    idx += nSend;
  }

  #endif
  return;
}

void SafeMpiAllReduce_fcmp(double complex *send, double complex *recv, int nData, MPI_Comm comm) {
  #ifdef _mpi_use
  int nSend = D_MpiSendMax; /* defined in global.h */
  int idx = 0;

  while(idx<nData) {
    if(idx+nSend > nData) nSend = nData-idx;
    MPI_Allreduce(send+idx,recv+idx,nSend,MPI_DOUBLE_COMPLEX,MPI_SUM,comm);
    idx += nSend;
  }

  #endif
  return;
}

void SafeMpiBcast_fcmp(double complex *buff, int nData, MPI_Comm comm) {
  #ifdef _mpi_use
  int nSend = D_MpiSendMax; /* defined in global.h */
  int idx = 0;

  while(idx<nData) {
    if(idx+nSend > nData) nSend = nData-idx;
    MPI_Bcast(buff+idx,nSend,MPI_DOUBLE_COMPLEX,0,comm);
    idx += nSend;
  }

  #endif
  return;
}
