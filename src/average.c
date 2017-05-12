/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Mitsuaki Kawamura, Takeo Kato, Masatoshi Imada.

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
 * calculate weighted averages of physical quantities
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include <complex.h>
#include "global.h"
#include "average.h"
#ifndef _SRC_AVERAGE
#define _SRC_AVERAGE

void weightAverageReduce(int n, double *vec, MPI_Comm comm);
void weightAverageReduce_fcmp(int n, double complex *vec, MPI_Comm comm);
void weightAverageReduce_real(int n, double *vec, MPI_Comm comm);


/* calculate average of Wc, Etot and Etot2 ;Sztot for fsz*/
/* All processes will have the result */
void WeightAverageWE(MPI_Comm comm) {
  const int n=4;//fsz
  double complex invW;
  int rank,size;
  double complex send[n], recv[n];
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  /* Wc, Etot and Etot2 */
  if(size>1) {
    send[0] = Wc;
    send[1] = Etot;
    send[2] = Etot2;
    send[3] = Sztot;

    SafeMpiAllReduce_fcmp(send,recv,n,comm);

    Wc    = recv[0];
    invW  = 1.0/Wc;
    Etot  = recv[1]*invW;
    Etot2 = recv[2]*invW;
    Sztot = recv[3]*invW;
  } else {
    invW  = 1.0/Wc;
    Etot  *= invW;
    Etot2 *= invW;
    Sztot *= invW;
  }

  return;
}

/* calculate average of SROptOO and SROptHO */
/* All processes will have the result */
void WeightAverageSROpt(MPI_Comm comm) {
  int i,n;
  double invW = 1.0/Wc;
  double complex *vec,*buf;
  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  /* SROptOO and SROptHO */ // except for SROptO 
  if(NSRCG == 0){
    n = 2*SROptSize*(2*SROptSize+1);
  }else{
    n = 2*SROptSize*3;
  }
  vec = SROptOO;
  if(size>1) {
    RequestWorkSpaceComplex(n);
    buf = GetWorkSpaceComplex(n);

    SafeMpiAllReduce_fcmp(vec,buf,n,comm);

    #pragma omp parallel for default(shared) private(i)
    #pragma loop noalias
    for(i=0;i<n;i++) vec[i] = buf[i] * invW;

    ReleaseWorkSpaceComplex();
 } else {
    #pragma omp parallel for default(shared) private(i)
    #pragma loop noalias
    for(i=0;i<n;i++) vec[i] *= invW;
  }

  return;
}

/* calculate average of SROptOO_real and SROptHO_real */
/* All processes will have the result */
void WeightAverageSROpt_real(MPI_Comm comm) {
  int i,n;
  double invW = 1.0/Wc;
  double *vec,*buf;
  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  /* SROptOO and SROptHO */ // except for SROptO 
  if(NSRCG == 0){
    n = SROptSize*(SROptSize+1);
  }else{
    n = SROptSize*3;
  }
  vec = SROptOO_real;
  if(size>1) {
    RequestWorkSpaceDouble(n);
    buf = GetWorkSpaceDouble(n);

    SafeMpiAllReduce(vec,buf,n,comm);

    #pragma omp parallel for default(shared) private(i)
    #pragma loop noalias
    for(i=0;i<n;i++) vec[i] = buf[i] * invW;

    ReleaseWorkSpaceDouble();
 } else {
    #pragma omp parallel for default(shared) private(i)
    #pragma loop noalias
    for(i=0;i<n;i++) vec[i] *= invW;
  }
  return;
}



/* calculate average of SROptOO and SROptHO */
/* All processes will have the result */
/*
void WeightAverageSROpt_real(MPI_Comm comm) {
  int i,n,j,int_x,int_y;
  double invW = 1.0/Wc;
  double *vec,*buf;
  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  n   = SROptSize*(SROptSize+1);
  vec = (double*)malloc(sizeof(double)*n);
  j   = 0;
  #pragma omp parallel for default(shared) private(i,int_x,int_y,j)
  #pragma loop noalias
  for(i=0;i<4*n;i++){
    int_x  = i%(2*SROptSize);
    int_y  = (i-int_x)/(2*SROptSize);
    if(int_x%2==0 && int_y%2==0){
      j      = int_x/2+(int_y/2)*SROptSize;
      vec[j] = creal(SROptOO[i]);// only real part TBC
    }
  }
  if(size>1) {
    RequestWorkSpaceDouble(n);
    buf = GetWorkSpaceDouble(n);

    SafeMpiAllReduce(vec,buf,n,comm);

    #pragma omp parallel for default(shared) private(i)
    #pragma loop noalias
    for(i=0;i<n;i++) vec[i] = buf[i] * invW;

    ReleaseWorkSpaceDouble();
 } else {
    #pragma omp parallel for default(shared) private(i)
    #pragma loop noalias
    for(i=0;i<n;i++) vec[i] *= invW;
  }
  #pragma omp parallel for default(shared) private(i,int_x,int_y,j)
  #pragma loop noalias
  for(i=0;i<4*n;i++){
    int_x  = i%(2*SROptSize);
    int_y  = (i-int_x)/(2*SROptSize);
    if(int_x%2==0 && int_y%2==0){
      j          = int_x/2+(int_y/2)*SROptSize;
      SROptOO[i] = vec[j];// only real part TBC
    }
  }
  free(vec);

  return;
}
*/

/* calculate average of Green functions */
/* Only rank=0 process will have the result */
void WeightAverageGreenFunc(MPI_Comm comm) {
  int n;
  double complex *vec;
    double *vec_real;
  /* Green functions */
  /* CisAjs, CisAjsCktAlt and CisAjsCktAltDC */
  n = NCisAjs+NCisAjsCktAlt+NCisAjsCktAltDC;
  vec = PhysCisAjs;
  weightAverageReduce_fcmp(n,vec,comm);
  
  if(NLanczosMode>0){
    /* QQQQ */
    n = NLSHam*NLSHam*NLSHam*NLSHam;
      if(AllComplexFlag==0 && iFlgOrbitalGeneral==0){
        vec_real=QQQQ_real;
        weightAverageReduce_real(n,vec_real,comm);
      }
      else{
        vec = QQQQ;
        weightAverageReduce_fcmp(n,vec,comm);
      }
    if(NLanczosMode>1){
      /* QCisAjsQ and QCisAjsCktAltQ */
      n = NLSHam*NLSHam*NCisAjs + NLSHam*NLSHam*NCisAjsCktAltDC;
        if(AllComplexFlag==0){
            vec_real=QCisAjsQ_real;
            weightAverageReduce_real(n,vec_real,comm);
        }
        else{
            vec = QCisAjsQ;
            weightAverageReduce_fcmp(n, vec, comm);
        }
    }
  }
  return;
}

void weightAverageReduce(int n, double *vec, MPI_Comm comm) {
  int i;
  const double invW = 1.0/Wc;
  double *buf;
  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  if(size>1) {
    RequestWorkSpaceDouble(n);
    buf = GetWorkSpaceDouble(n);

    SafeMpiReduce(vec,buf,n,comm);
    if(rank==0) {
      #pragma omp parallel for default(shared) private(i)
      #pragma loop noalias
      for(i=0;i<n;i++) vec[i] = buf[i] * invW;
    }

    ReleaseWorkSpaceDouble();
  } else {
    #pragma omp parallel for default(shared) private(i)
    #pragma loop noalias
    for(i=0;i<n;i++) vec[i] *= invW;
  }

  return;
}

void weightAverageReduce_fcmp(int n, double  complex *vec, MPI_Comm comm) {
  int i;
  const double complex invW = 1.0/Wc;
  double complex *buf;
  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  if(size>1) {
    RequestWorkSpaceComplex(n);
    buf = GetWorkSpaceComplex(n);

    SafeMpiReduce_fcmp(vec,buf,n,comm);
    if(rank==0) {
      #pragma omp parallel for default(shared) private(i)
      #pragma loop noalias
      for(i=0;i<n;i++) vec[i] = buf[i] * invW;
    }

    ReleaseWorkSpaceComplex();
  } else {
    #pragma omp parallel for default(shared) private(i)
    #pragma loop noalias
    for(i=0;i<n;i++) vec[i] *= invW;
  }

  return;
}

void weightAverageReduce_real(int n, double *vec, MPI_Comm comm) {
    int i;
    const double invW = 1.0/Wc;
    double *buf;
    int rank,size;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);

    if(size>1) {
        RequestWorkSpaceDouble(n);
        buf = GetWorkSpaceDouble(n);

        SafeMpiReduce(vec,buf,n,comm);
        if(rank==0) {
#pragma omp parallel for default(shared) private(i)
#pragma loop noalias
            for(i=0;i<n;i++) vec[i] = buf[i] * invW;
        }

        ReleaseWorkSpaceDouble();
    } else {
#pragma omp parallel for default(shared) private(i)
#pragma loop noalias
        for(i=0;i<n;i++) vec[i] *= invW;
    }

    return;
}

#endif
