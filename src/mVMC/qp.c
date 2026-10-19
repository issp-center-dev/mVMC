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
#include "qp.h"
#ifndef _QP_SRC
#define _QP_SRC
#include "global.h"
#include <complex.h>
#include "legendrepoly.h"
#include "workspace.h"
#include "gauleg.h"
#include "gc_exact_sum.h"
#include "gc_spin_support.h"

/* Set once from the distributed input coefficients; read-only during a run. */
static int GCTransWeightSumExactlyZero = 0;

/* calculate SPGLCos, SPGLSin and QPFullWeight */
void InitQPWeight() {
  int i,j,idx;
  double *beta;
  double *weight;
  double w;

  GCTransWeightSumExactlyZero = 0;
  if (FlagGrandCanonical &&
      (GCExactComplexSumSelfTest() != 0 ||
       GCExactComplexSumIsZero(ParaQPTrans, NMPTrans,
                               &GCTransWeightSumExactlyZero) != 0)) {
    fprintf(stderr, "Error: GC translation-weight exact-sum initialization failed.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  RequestWorkSpaceDouble(2*NSPGaussLeg);
  beta   = GetWorkSpaceDouble(NSPGaussLeg);
  weight = GetWorkSpaceDouble(NSPGaussLeg);

  if(NSPGaussLeg==1) {
    beta[0] = 0.0;
    weight[0] = 1.0;
    SPGLCos[0] = 1.0;
    SPGLSin[0] = 0.0;
    SPGLCosSin[0] = 0.0;
    SPGLCosCos[0] = 1.0;
    SPGLSinSin[0] = 0.0;

    for(j=0;j<NMPTrans;j++) {
      //printf("XDEBUG: j=%d %lf %lf\n",j, creal(ParaQPTrans[j]), cimag(ParaQPTrans[j]));
      QPFixWeight[j] = ParaQPTrans[j];
    }

  } else {
    GaussLeg(0, M_PI, beta, weight, NSPGaussLeg);

    #pragma omp parallel for default(shared) private(i,j,w,idx)
    for(i=0;i<NSPGaussLeg;i++) {
      SPGLCos[i]    = cos(0.5 * beta[i])+0*I;    //TBC
      SPGLSin[i]    = sin(0.5 * beta[i])+0*I;    //TBC
      SPGLCosSin[i] = SPGLCos[i]*SPGLSin[i]+0*I; //TBC
      SPGLCosCos[i] = SPGLCos[i]*SPGLCos[i]+0*I; //TBC
      SPGLSinSin[i] = SPGLSin[i]*SPGLSin[i]+0*I; //TBC
      
      w = 0.5*sin(beta[i])*weight[i]*LegendrePoly(cos(beta[i]), NSPStot);
      
      for(j=0;j<NMPTrans;j++) {
        idx = i + j*NSPGaussLeg;
        QPFixWeight[idx] = w * ParaQPTrans[j];
      }
    }
  }

  UpdateQPWeight();

  ReleaseWorkSpaceDouble();
  return;
}

#ifdef MVMC_ENABLE_FAULT_INJECTION
/* Testing only: the production quadrature coefficients, so that independent
 * oracles can use exactly the binary64 values of this run.  Header
 * "GC_SPIN_QP_V1 G <G> Q <NQPFull>", then "q spidx cs cc ss w_re w_im". */
void GCTestDumpSpinQP(MPI_Comm comm) {
  int rank;
  int failed = 0;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0 && FlagGrandCanonical) {
    const char *path = getenv("MVMC_GC_SPIN_QP_DUMP");
    if (path != NULL && *path != '\0') {
      FILE *fp = fopen(path, "w");
      if (fp == NULL) {
        failed = 1;
      } else {
        int q;
        if (fprintf(fp, "GC_SPIN_QP_V1 G %d Q %d\n", NSPGaussLeg, NQPFull) <
            0)
          failed = 1;
        for (q = 0; q < NQPFull; q++) {
          const int s = q % NSPGaussLeg;
          if (fprintf(fp, "%d %d %a %a %a %a %a\n", q, s,
                      creal(SPGLCosSin[s]), creal(SPGLCosCos[s]),
                      creal(SPGLSinSin[s]), creal(QPFullWeight[q]),
                      cimag(QPFullWeight[q])) < 0)
            failed = 1;
        }
        if (fclose(fp) != 0) failed = 1;
      }
    }
  }
  /* No-MPI builds have no MPI_Bcast/MPI_INT stub. */
#ifdef _mpi_use
  MPI_Bcast(&failed, 1, MPI_INT, 0, comm);
#endif
  if (failed) {
    if (rank == 0) fprintf(stderr, "Error: GC spin qp dump failed.\n");
    MPI_Abort(comm, EXIT_FAILURE);
  }
}
#endif


/* Calculate logarithm of inner product <phi|L|x> */
double complex CalculateLogIP_fcmp(double complex * const pfM, const int qpStart, const int qpEnd, MPI_Comm comm) {
  const int qpNum = qpEnd-qpStart;
  double complex ip=0.0+0.0*I, ip2;
  int qpidx;
  int size;
  MPI_Comm_size(comm,&size);

  #pragma loop noalias
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    ip += QPFullWeight[qpidx+qpStart] * pfM[qpidx];
  }
  if(size>1) {
    MPI_Allreduce(&ip, &ip2, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);//TBC
    ip = ip2;
  }
  //printf("QDEBUG: ip= %lf %lf\n",creal(clog(ip)),cimag(clog(ip)));
  return clog(ip);
}

/* Calculate inner product <phi|L|x> */
double complex CalculateIP_fcmp(double complex * const pfM, const int qpStart, const int qpEnd, MPI_Comm comm) {
  const int qpNum = qpEnd-qpStart;
  double complex ip=0.0+0.0*I, ip2;
  int qpidx;
  int size;
  MPI_Comm_size(comm,&size);

  #pragma loop noalias
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    ip += QPFullWeight[qpidx+qpStart] * pfM[qpidx];
    //printf("DEBUG: qpidx =%d: %lf %lf   \n",qpidx,creal( QPFullWeight[qpidx+qpStart] ),creal(pfM[qpidx]));
  }
  if(size>1) {
    MPI_Allreduce(&ip, &ip2, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);
    ip = ip2;
  }
  return ip;
}

/* Vacuum/full filling are translation invariant; pure-triplet S=0 states
 * also exclude odd pair number. Remove only finite roundoff in those proven
 * sectors, keeping the QP reduction order and nonfinite propagation. */
int GCProjectionSectorZero(const int newN) {
  return (FlagGrandCanonical && GCTransWeightSumExactlyZero &&
          (newN == 0 || newN == Nsite2)) ||
         GCSpinForbiddenSector(FlagGrandCanonical, iFlgOrbitalGeneral,
                               NSPGaussLeg, NSPStot, GCSpinTriplet, newN);
}

double complex CalculateIPGC_fcmp(double complex *const pfM, const int newN,
                                  const int qpStart, const int qpEnd,
                                  MPI_Comm comm) {
  const double complex ip = CalculateIP_fcmp(pfM, qpStart, qpEnd, comm);
  if (GCProjectionSectorZero(newN) &&
      isfinite(creal(ip)) && isfinite(cimag(ip)))
    return 0.0 + 0.0 * I;
  return ip;
}

double complex CalculateLogIPGC_fcmp(double complex *const pfM, const int newN,
                                     const int qpStart, const int qpEnd,
                                     MPI_Comm comm) {
  return clog(CalculateIPGC_fcmp(pfM, newN, qpStart, qpEnd, comm));
}

void UpdateQPWeight() {
  int i,j,offset;
  double complex tmp; //TBC

  if(FlagOptTrans>0) {
    for(i=0;i<NOptTrans;i++) {
      offset = i*NQPFix;
      tmp = OptTrans[i];
      for(j=0;j<NQPFix;j++) {
        QPFullWeight[offset+j] = tmp * QPFixWeight[j];
      }
    }
  } else {
    for(j=0;j<NQPFix;j++) {
      QPFullWeight[j] = QPFixWeight[j];
    }
  }

  return;
}

#endif
