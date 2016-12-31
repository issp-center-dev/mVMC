/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Mitsuaki Kawamura, Takeo Kato, Masatoshi Imada.

his program is developed based on the mVMC-mini program
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
 * safe version of mpi functions
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#pragma once
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
