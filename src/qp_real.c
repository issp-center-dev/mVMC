/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Quantum Projection
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

double CalculateLogIP_real(double* const pfM, const int qpStart, const int qpEnd, MPI_Comm comm);
double CalculateIP_real(double* const pfM, const int qpStart, const int qpEnd, MPI_Comm comm);


/* Calculate logarithm of inner product <phi|L|x> */
double CalculateLogIP_real(double* const pfM, const int qpStart, const int qpEnd, MPI_Comm comm) {
  const int qpNum = qpEnd-qpStart;
  double  ip=0.0+0.0*I, ip2;
  int qpidx;
  int size;
  MPI_Comm_size(comm,&size);

  #pragma loop noalias
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    ip += creal(QPFullWeight[qpidx+qpStart]) * pfM[qpidx];
  }
  if(size>1) {
    MPI_Allreduce(&ip, &ip2, 1, MPI_DOUBLE, MPI_SUM, comm);//TBC
    ip = ip2;
  }
  return clog(ip);
}

/* Calculate inner product <phi|L|x> */
double  CalculateIP_real(double* const pfM, const int qpStart, const int qpEnd, MPI_Comm comm) {
  const int qpNum = qpEnd-qpStart;
  double ip=0.0+0.0*I, ip2;
  int qpidx;
  int size;
  MPI_Comm_size(comm,&size);

  #pragma loop noalias
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    ip += creal(QPFullWeight[qpidx+qpStart]) * pfM[qpidx];
    //printf("DEBUG: qpidx =%d: %lf %lf   \n",qpidx,creal( QPFullWeight[qpidx+qpStart] ),creal(pfM[qpidx]));
  }
  if(size>1) {
    MPI_Allreduce(&ip, &ip2, 1, MPI_DOUBLE, MPI_SUM, comm);
    ip = ip2;
  }
  return ip;
}
