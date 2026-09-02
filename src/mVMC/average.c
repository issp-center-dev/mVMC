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
 * calculate weighted averages of physical quantities
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include <complex.h>
#include <stdint.h>
#include "global.h"
#include "average.h"
#include "physcal_lanczos2.h"
#ifndef _SRC_AVERAGE
#define _SRC_AVERAGE

void weightAverageReduce(int n, double *vec, MPI_Comm comm);
void weightAverageReduce_fcmp(int n, double complex *vec, MPI_Comm comm);
void weightAverageReduce_real(int n, double *vec, MPI_Comm comm);

static void CheckPowerLanczosSupport(MPI_Comm comm) {
  PowerLanczosSupportDiagnostic diagnostic;
  PowerLanczosSupportStatus supportStatus = POWER_LANCZOS_SUPPORT_INVALID;
  Lanczos2SolveStatus outputStatus = LANCZOS2_SOLVE_OK;
  double *reduced = PowerLanczosSupportSampleData;
  size_t analyzedCapacity = (size_t)PowerLanczosSupportSampleCapacity;
  const int count =
      PowerLanczosSupportSampleCapacity *
      POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH;
  int statusCode = POWER_LANCZOS_SUPPORT_INVALID;
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  if(size>1) {
    /*
     * Each rank stores only its local NVMCSample rows.  A root-only gather
     * preserves the previous rank-major sample order without replicating a
     * world-sized audit buffer, or its reduction buffer, on every rank.
     */
    if(analyzedCapacity > SIZE_MAX / (size_t)size) {
      if(rank == 0) {
        fprintf(stderr,
                "Error: power-Lanczos gathered support capacity is too "
                "large.\n");
      }
      MPI_Abort(comm,EXIT_FAILURE);
    }
    analyzedCapacity *= (size_t)size;
    if(rank == 0) {
      if(analyzedCapacity >
             SIZE_MAX / POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH ||
         analyzedCapacity * POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH >
             SIZE_MAX / sizeof(double)) {
        fprintf(stderr,
                "Error: power-Lanczos gathered support buffer is too "
                "large.\n");
        MPI_Abort(comm,EXIT_FAILURE);
      }
      reduced = (double*)malloc(
          analyzedCapacity * POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH *
          sizeof(double));
      if(reduced == NULL) {
        fprintf(stderr,
                "Error: failed to allocate power-Lanczos support gather "
                "buffer.\n");
        MPI_Abort(comm,EXIT_FAILURE);
      }
    }
#ifdef _mpi_use
    MPI_Gather(PowerLanczosSupportSampleData,count,MPI_DOUBLE,
               reduced,count,MPI_DOUBLE,0,comm);
#endif
  }
  if(rank == 0) {
    supportStatus = AnalyzePowerLanczosSupport(
        reduced,analyzedCapacity,
        NLanczosStep==2,&diagnostic);
    outputStatus = WritePowerLanczosSupportDiagnostic(
        FileLSSupport,NLanczosStep,NLanczosSupportMode,&diagnostic);
    if(outputStatus != LANCZOS2_SOLVE_OK) {
      fprintf(stderr,
              "Error: power-Lanczos support diagnostic output failed: %s.\n",
              Lanczos2SolveError(outputStatus));
      supportStatus = POWER_LANCZOS_SUPPORT_INVALID;
    }
    statusCode = (int)supportStatus;
    if(supportStatus != POWER_LANCZOS_SUPPORT_PASS) {
      fprintf(stderr,
              "%s: power-Lanczos support %s: "
              "M02-M11=(%.6e,%.6e), score=%.6e, "
              "M03-M12=(%.6e,%.6e), score=%.6e, samples=%zu.\n",
              NLanczosSupportMode == 0 ? "Error" : "Warning",
              PowerLanczosSupportStatusName(supportStatus),
              creal(diagnostic.delta2),cimag(diagnostic.delta2),
              diagnostic.score2,creal(diagnostic.delta3),
              cimag(diagnostic.delta3),diagnostic.score3,
              diagnostic.sampleCount);
    }
  }
#ifdef _mpi_use
  MPI_Bcast(&statusCode,1,MPI_INT,0,comm);
#endif
  if(size>1 && rank == 0) free(reduced);
  if(NLanczosSupportMode == 0 &&
     statusCode != POWER_LANCZOS_SUPPORT_PASS) {
    MPI_Abort(comm,EXIT_FAILURE);
  }
}


/* calculate average of Wc, Etot and Etot2 ;Sztot,Sztot2 for fsz*/
/* All processes will have the result */
void WeightAverageWE(MPI_Comm comm) {
  if(FlagGrandCanonical != 0) {
    const int nGC=7;
    double complex invWGC;
    int rankGC,sizeGC;
    double complex sendGC[nGC],recvGC[nGC];
    MPI_Comm_rank(comm,&rankGC);
    MPI_Comm_size(comm,&sizeGC);
    if(sizeGC>1) {
      sendGC[0] = Wc;
      sendGC[1] = Etot;
      sendGC[2] = Etot2;
      sendGC[3] = Sztot;
      sendGC[4] = Sztot2;
      sendGC[5] = Ntot;
      sendGC[6] = Ntot2;
      SafeMpiAllReduce_fcmp(sendGC,recvGC,nGC,comm);
      Wc     = recvGC[0];
      invWGC = 1.0/Wc;
      Etot   = recvGC[1]*invWGC;
      Etot2  = recvGC[2]*invWGC;
      Sztot  = recvGC[3]*invWGC;
      Sztot2 = recvGC[4]*invWGC;
      Ntot   = recvGC[5]*invWGC;
      Ntot2  = recvGC[6]*invWGC;
    } else {
      invWGC = 1.0/Wc;
      Etot   *= invWGC;
      Etot2  *= invWGC;
      Sztot  *= invWGC;
      Sztot2 *= invWGC;
      Ntot   *= invWGC;
      Ntot2  *= invWGC;
    }
    return;
  }
  const int n=5;//fsz
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
    send[4] = Sztot2;

    SafeMpiAllReduce_fcmp(send,recv,n,comm);

    Wc     = recv[0];
    invW   = 1.0/Wc;
    Etot   = recv[1]*invW;
    Etot2  = recv[2]*invW;
    Sztot  = recv[3]*invW;
    Sztot2 = recv[4]*invW;
  } else {
    invW  = 1.0/Wc;
    Etot  *= invW;
    Etot2 *= invW;
    Sztot *= invW;
    Sztot2 *= invW;
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
  //n = NCisAjs+NCisAjsCktAlt+NCisAjsCktAltDC;
  n = NCisAjs+NCisAjsCktAlt+NCisAjsCktAltDC+NTwist;
  vec = PhysCisAjs;
  weightAverageReduce_fcmp(n,vec,comm);

  if (NNBodyG > 0) {
    vec = PhysNBodyG;
    weightAverageReduce_fcmp(NNBodyG,vec,comm);
  }

  if (NAnomalousG > 0) {
    weightAverageReduce_fcmp(NAnomalousG, PhysAnomalousG, comm);
  }

  if(NLanczosMode>0 && NLanczosStep==1){
    /* QQQQ */
    n = NLSHam*NLSHam*NLSHam*NLSHam;
    //if(AllComplexFlag==0 && iFlgOrbitalGeneral==0){
    if(AllComplexFlag==0){
      vec_real=QQQQ_real;
      weightAverageReduce_real(n,vec_real,comm);
    }
    else{
      vec = QQQQ;
      weightAverageReduce_fcmp(n,vec,comm);
    }
    if(NLanczosMode>1){
      /* QCisAjsQ and QCisAjsCktAltQ */
      n = NLSHam*NLSHam*NCisAjs + NLSHam*NLSHam*(NCisAjsCktAltDC+NCisAjsCktAlt);
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
  if(NLanczosMode>0 && NLanczosStep==2){
    Lanczos2SolveStatus lanczos2Status;
    int rank = 0;
    MPI_Comm_rank(comm,&rank);
    LS2MomentBasisShift = creal(Etot);
    /*
     * Every rank must center its local moments in the same basis before
     * their reduction.  WeightAverageWE already all-reduces Etot; keep
     * that invariant local to this use site as defensive hardening.
     */
    SafeMpiBcast(&LS2MomentBasisShift,1,comm);
    if(AllComplexFlag==0){
      lanczos2Status = BuildLanczos2ShiftedMomentReal(
          LS2Moment_real,LS2SamplePower_real,LS2SampleWeight,
          LS2SampleValid,(size_t)NVMCSample,LS2MomentBasisShift);
      if(lanczos2Status != LANCZOS2_SOLVE_OK) {
        fprintf(stderr,
                "Error: Lanczos2 shifted moment construction failed on "
                "rank %d: %s.\n",
                rank,Lanczos2SolveError(lanczos2Status));
        MPI_Abort(comm,EXIT_FAILURE);
      }
      weightAverageReduce_real(LANCZOS2_MOMENT_COUNT,LS2Moment_real,comm);
    }else{
      lanczos2Status = BuildLanczos2ShiftedMomentComplex(
          LS2Moment,LS2SamplePower,LS2SampleWeight,LS2SampleValid,
          (size_t)NVMCSample,LS2MomentBasisShift);
      if(lanczos2Status != LANCZOS2_SOLVE_OK) {
        fprintf(stderr,
                "Error: Lanczos2 shifted moment construction failed on "
                "rank %d: %s.\n",
                rank,Lanczos2SolveError(lanczos2Status));
        MPI_Abort(comm,EXIT_FAILURE);
      }
      weightAverageReduce_fcmp(LANCZOS2_MOMENT_COUNT,LS2Moment,comm);
    }
  }
  if(NLanczosMode>0){
    CheckPowerLanczosSupport(comm);
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
