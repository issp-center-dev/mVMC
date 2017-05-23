/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is developed based on the mVMC-mini program
(https://github.com/fiber-miniapp/mVMC-mini)
which follows "The BSD 3-Clause License".

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details. 

You should have received a copy of the GNU General Public License 
along with this program. If not, see http://www.gnu.org/licenses/. 
*/
/*-------------------------------------------------------------
 * Variational Monte Carlo
 * Quantum Projection
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "qp_real.h"

#ifndef _SRC_QP_REAL
#define _SRC_QP_REAL

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
#endif

