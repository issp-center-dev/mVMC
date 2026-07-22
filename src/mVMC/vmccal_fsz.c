/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

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
 * calculate physical quantities
 *-------------------------------------------------------------
 * by Satoshi Morita 
 *-------------------------------------------------------------*/
#ifndef _SRC_VMCCAL_FSZ
#define _SRC_VMCCAL_FSZ

#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vmccal_fsz.h"
#include "matrix.h"
#include "calham_fsz_real.h"
#include "calham_fsz.h"
#include "calgrn_fsz.h"
#include "lslocgrn_fsz.h"
#include "lslocgrn_bf_fsz.h"

static int isFiniteLSLQFSZComplex(const double complex *lslq) {
  int i;
  for(i=0; i<NLSHam*NLSHam; i++) {
    if(!isfinite(creal(lslq[i])) || !isfinite(cimag(lslq[i]))) return 0;
  }
  return 1;
}

static int isFiniteLSLQFSZReal(const double *lslq) {
  int i;
  for(i=0; i<NLSHam*NLSHam; i++) {
    if(!isfinite(lslq[i])) return 0;
  }
  return 1;
}

static long getLanczosNonfiniteInjectionSample(void) {
  const char *value = getenv("MVMC_LANCZOS_TEST_NONFINITE_SAMPLE");
  char *end = NULL;
  long sample;
  if(value == NULL || value[0] == '\0') return -1;
  errno = 0;
  sample = strtol(value, &end, 10);
  if(errno != 0 || end == value || *end != '\0' || sample < 0) return -2;
  return sample;
}

static void finalizeLanczosFSZAccounting(MPI_Comm comm, int rank,
                                         long long checkedLocal,
                                         long long rejectedLocal) {
  long long local[2] = {checkedLocal, rejectedLocal};
  long long global[2] = {checkedLocal, rejectedLocal};
#ifdef _mpi_use
  MPI_Allreduce(local, global, 2, MPI_LONG_LONG, MPI_SUM, comm);
#endif
  if(global[0] <= 0) {
    if(rank == 0) {
      fprintf(stderr, "Error: FSZ Lanczos checked zero samples.\n");
    }
    MPI_Abort(comm, EXIT_FAILURE);
  }
  if(global[1] > 0 && global[1]*100 <= global[0]) {
    if(rank == 0) {
      fprintf(stderr,
              "warning: FSZ Lanczos rejected %lld/%lld samples due to "
              "non-finite local vectors or numerical candidate rebuild failures.\n",
              global[1], global[0]);
    }
  } else if(global[1]*100 > global[0]) {
    if(rank == 0) {
      fprintf(stderr,
              "Error: FSZ Lanczos rejected %lld/%lld samples (>1%%); biased output will not be written.\n",
              global[1], global[0]);
    }
    MPI_Abort(comm, EXIT_FAILURE);
  }
}

static void clearStoredOSampleRange_fsz(const int sampleStart, const int sampleEnd) {
  const int sampleSize = sampleEnd - sampleStart;

  if(sampleSize <= 0) return;

  if(AllComplexFlag == 0) {
    const size_t offset = (size_t)sampleStart * (size_t)SROptSize;
    const size_t nStore = (size_t)sampleSize * (size_t)SROptSize;
    size_t i;
    for(i=0; i<nStore; i++) SROptO_Store_real[offset+i] = 0.0;
  } else {
    const size_t stride = 2 * (size_t)SROptSize;
    const size_t offset = (size_t)sampleStart * stride;
    const size_t nStore = (size_t)sampleSize * stride;
    size_t i;
    for(i=0; i<nStore; i++) SROptO_Store[offset+i] = 0.0 + 0.0*I;
  }
}

static int calculateBFFSZIPForFD(int *eleIdx, int *eleSpn, int *eleNum,
                                 const int *eleProjBFCnt,
                                 const int qpStart, const int qpEnd,
                                 double complex *ip) {
  int info;

  MakeSlaterElmBF_fsz(eleNum, eleProjBFCnt);
  info = CalculateMAll_BF_fsz(eleIdx, eleSpn, qpStart, qpEnd);
  *ip = CalculateIP_fcmp(PfM, qpStart, qpEnd, MPI_COMM_SELF);
  return info;
}

static void dumpBFFSZSRDiffCheck(const char *path, int *eleIdx, int *eleSpn,
                                 int *eleNum, const int *eleProjBFCnt,
                                 const int qpStart, const int qpEnd,
                                 const int sample) {
  FILE *fp;
  double complex *diffNoBF;
  double complex *analyticPacked;
  double complex *analyticProjBF;
  double complex *analyticSlater;
  double complex *projBFStore;
  double complex *slaterStore;
  double complex ipNoBF = 0.0 + 0.0*I;
  double complex ipBF = 0.0 + 0.0*I;
  double complex ipPlus = 0.0 + 0.0*I;
  double complex ipMinus = 0.0 + 0.0*I;
  double complex fd;
  double complex projBFAnalyticAtMax = 0.0 + 0.0*I;
  double complex projBFFDAtMax = 0.0 + 0.0*I;
  double complex orbitalAnalyticAtMax = 0.0 + 0.0*I;
  double complex orbitalFDAtMax = 0.0 + 0.0*I;
  const double h = 1.0e-6;
  double maxIdentityOrbitalDiff = 0.0;
  double maxProjBFRealDiff = 0.0;
  double maxProjBFImagDiff = 0.0;
  double maxProjBFDiff = 0.0;
  double maxOrbitalRealDiff = 0.0;
  double maxOrbitalImagDiff = 0.0;
  double maxOrbitalDiff = 0.0;
  int maxIdentityOrbitalIdx = -1;
  int maxProjBFRealIdx = -1;
  int maxProjBFImagIdx = -1;
  int maxProjBFIdx = -1;
  int maxProjBFIsImag = 0;
  int maxOrbitalRealIdx = -1;
  int maxOrbitalImagIdx = -1;
  int maxOrbitalIdx = -1;
  int maxOrbitalIsImag = 0;
  int nonzeroProjBFFDCount = 0;
  int nonzeroOrbitalFDCount = 0;
  int nanCount = 0;
  int fdFailCount = 0;
  int infoNoBF;
  int infoBF;
  int infoPlus;
  int infoMinus;
  int noBFReady;
  int bfReady;
  int idx;

  if(path == NULL || path[0] == '\0') return;
  if(NProjBF <= 0 || ProjBF == NULL || NSlater <= 0) return;

  diffNoBF = (double complex *)calloc((size_t)2 * (size_t)NSlater, sizeof(double complex));
  analyticPacked = (double complex *)calloc((size_t)2 * (size_t)(NProjBF + NSlater),
                                            sizeof(double complex));
  projBFStore = (double complex *)calloc((size_t)NProjBF, sizeof(double complex));
  slaterStore = (double complex *)calloc((size_t)NSlater, sizeof(double complex));
  if(diffNoBF == NULL || analyticPacked == NULL ||
      projBFStore == NULL || slaterStore == NULL) {
    fprintf(stderr, "Error: memory allocation failed for BF-FSZ SR diff dump.\n");
    free(diffNoBF);
    free(analyticPacked);
    free(projBFStore);
    free(slaterStore);
    return;
  }
  analyticProjBF = analyticPacked;
  analyticSlater = analyticPacked + 2*NProjBF;

  for(idx=0;idx<NProjBF;idx++) projBFStore[idx] = ProjBF[idx];
  for(idx=0;idx<NSlater;idx++) slaterStore[idx] = Slater[idx];

  infoNoBF = CalculateMAll_fsz(eleIdx, eleSpn, qpStart, qpEnd);
  ipNoBF = CalculateIP_fcmp(PfM, qpStart, qpEnd, MPI_COMM_SELF);
  noBFReady = (infoNoBF == 0 && isfinite(creal(ipNoBF)) &&
               isfinite(cimag(ipNoBF)) && cabs(ipNoBF) > 0.0);
  if(noBFReady) {
    SlaterElmDiff_fsz(diffNoBF, ipNoBF, eleIdx, eleSpn);
  } else {
    fdFailCount++;
  }

  infoBF = calculateBFFSZIPForFD(eleIdx, eleSpn, eleNum, eleProjBFCnt,
                                 qpStart, qpEnd, &ipBF);
  bfReady = (infoBF == 0 && isfinite(creal(ipBF)) &&
             isfinite(cimag(ipBF)) && cabs(ipBF) > 0.0);
  if(bfReady) {
    BackFlowDiff_fsz(analyticProjBF, ipBF, eleIdx, eleSpn, eleNum, eleProjBFCnt);
    SlaterElmBFDiff_fsz(analyticSlater, ipBF, eleIdx, eleSpn, eleNum, eleProjBFCnt);
  } else {
    fdFailCount++;
  }

  if(noBFReady && bfReady) {
    for(idx=0;idx<2*NSlater;idx++) {
      const double diff = cabs(analyticSlater[idx] - diffNoBF[idx]);
      if(!isfinite(diff)) {
        nanCount++;
      } else if(diff > maxIdentityOrbitalDiff) {
        maxIdentityOrbitalDiff = diff;
        maxIdentityOrbitalIdx = idx;
      }
    }
  }

  if(bfReady) {
    for(idx=0;idx<NProjBF;idx++) {
      double diff;

      ProjBF[idx] = projBFStore[idx] + h;
      infoPlus = calculateBFFSZIPForFD(eleIdx, eleSpn, eleNum, eleProjBFCnt,
                                       qpStart, qpEnd, &ipPlus);
      ProjBF[idx] = projBFStore[idx] - h;
      infoMinus = calculateBFFSZIPForFD(eleIdx, eleSpn, eleNum, eleProjBFCnt,
                                        qpStart, qpEnd, &ipMinus);
      ProjBF[idx] = projBFStore[idx];

      if(infoPlus != 0 || infoMinus != 0 ||
          !isfinite(creal(ipPlus)) || !isfinite(cimag(ipPlus)) ||
          !isfinite(creal(ipMinus)) || !isfinite(cimag(ipMinus))) {
        fdFailCount++;
      } else {
        fd = ((ipPlus - ipMinus) / (2.0*h)) / ipBF;
        diff = cabs(fd - analyticProjBF[2*idx]);
        if(!isfinite(diff) || !isfinite(cabs(fd))) {
          nanCount++;
        } else {
          if(cabs(fd) > 1.0e-12) nonzeroProjBFFDCount++;
          if(diff > maxProjBFRealDiff) {
            maxProjBFRealDiff = diff;
            maxProjBFRealIdx = idx;
          }
          if(diff > maxProjBFDiff) {
            maxProjBFDiff = diff;
            maxProjBFIdx = idx;
            maxProjBFIsImag = 0;
            projBFAnalyticAtMax = analyticProjBF[2*idx];
            projBFFDAtMax = fd;
          }
        }
      }

      ProjBF[idx] = projBFStore[idx] + I*h;
      infoPlus = calculateBFFSZIPForFD(eleIdx, eleSpn, eleNum, eleProjBFCnt,
                                       qpStart, qpEnd, &ipPlus);
      ProjBF[idx] = projBFStore[idx] - I*h;
      infoMinus = calculateBFFSZIPForFD(eleIdx, eleSpn, eleNum, eleProjBFCnt,
                                        qpStart, qpEnd, &ipMinus);
      ProjBF[idx] = projBFStore[idx];

      if(infoPlus != 0 || infoMinus != 0 ||
          !isfinite(creal(ipPlus)) || !isfinite(cimag(ipPlus)) ||
          !isfinite(creal(ipMinus)) || !isfinite(cimag(ipMinus))) {
        fdFailCount++;
      } else {
        fd = ((ipPlus - ipMinus) / (2.0*h)) / ipBF;
        diff = cabs(fd - analyticProjBF[2*idx+1]);
        if(!isfinite(diff) || !isfinite(cabs(fd))) {
          nanCount++;
        } else {
          if(cabs(fd) > 1.0e-12) nonzeroProjBFFDCount++;
          if(diff > maxProjBFImagDiff) {
            maxProjBFImagDiff = diff;
            maxProjBFImagIdx = idx;
          }
          if(diff > maxProjBFDiff) {
            maxProjBFDiff = diff;
            maxProjBFIdx = idx;
            maxProjBFIsImag = 1;
            projBFAnalyticAtMax = analyticProjBF[2*idx+1];
            projBFFDAtMax = fd;
          }
        }
      }
    }

    for(idx=0;idx<NSlater;idx++) {
      double diff;

      Slater[idx] = slaterStore[idx] + h;
      infoPlus = calculateBFFSZIPForFD(eleIdx, eleSpn, eleNum, eleProjBFCnt,
                                       qpStart, qpEnd, &ipPlus);
      Slater[idx] = slaterStore[idx] - h;
      infoMinus = calculateBFFSZIPForFD(eleIdx, eleSpn, eleNum, eleProjBFCnt,
                                        qpStart, qpEnd, &ipMinus);
      Slater[idx] = slaterStore[idx];

      if(infoPlus != 0 || infoMinus != 0 ||
          !isfinite(creal(ipPlus)) || !isfinite(cimag(ipPlus)) ||
          !isfinite(creal(ipMinus)) || !isfinite(cimag(ipMinus))) {
        fdFailCount++;
      } else {
        fd = ((ipPlus - ipMinus) / (2.0*h)) / ipBF;
        diff = cabs(fd - analyticSlater[2*idx]);
        if(!isfinite(diff) || !isfinite(cabs(fd))) {
          nanCount++;
        } else {
          if(cabs(fd) > 1.0e-12) nonzeroOrbitalFDCount++;
          if(diff > maxOrbitalRealDiff) {
            maxOrbitalRealDiff = diff;
            maxOrbitalRealIdx = idx;
          }
          if(diff > maxOrbitalDiff) {
            maxOrbitalDiff = diff;
            maxOrbitalIdx = idx;
            maxOrbitalIsImag = 0;
            orbitalAnalyticAtMax = analyticSlater[2*idx];
            orbitalFDAtMax = fd;
          }
        }
      }

      Slater[idx] = slaterStore[idx] + I*h;
      infoPlus = calculateBFFSZIPForFD(eleIdx, eleSpn, eleNum, eleProjBFCnt,
                                       qpStart, qpEnd, &ipPlus);
      Slater[idx] = slaterStore[idx] - I*h;
      infoMinus = calculateBFFSZIPForFD(eleIdx, eleSpn, eleNum, eleProjBFCnt,
                                        qpStart, qpEnd, &ipMinus);
      Slater[idx] = slaterStore[idx];

      if(infoPlus != 0 || infoMinus != 0 ||
          !isfinite(creal(ipPlus)) || !isfinite(cimag(ipPlus)) ||
          !isfinite(creal(ipMinus)) || !isfinite(cimag(ipMinus))) {
        fdFailCount++;
      } else {
        fd = ((ipPlus - ipMinus) / (2.0*h)) / ipBF;
        diff = cabs(fd - analyticSlater[2*idx+1]);
        if(!isfinite(diff) || !isfinite(cabs(fd))) {
          nanCount++;
        } else {
          if(cabs(fd) > 1.0e-12) nonzeroOrbitalFDCount++;
          if(diff > maxOrbitalImagDiff) {
            maxOrbitalImagDiff = diff;
            maxOrbitalImagIdx = idx;
          }
          if(diff > maxOrbitalDiff) {
            maxOrbitalDiff = diff;
            maxOrbitalIdx = idx;
            maxOrbitalIsImag = 1;
            orbitalAnalyticAtMax = analyticSlater[2*idx+1];
            orbitalFDAtMax = fd;
          }
        }
      }
    }
  }

  for(idx=0;idx<NProjBF;idx++) ProjBF[idx] = projBFStore[idx];
  for(idx=0;idx<NSlater;idx++) Slater[idx] = slaterStore[idx];
  MakeSlaterElmBF_fsz(eleNum, eleProjBFCnt);

  fp = fopen(path, "w");
  if(fp == NULL) {
    fprintf(stderr, "Error: failed to open BF-FSZ SR diff dump file: %s\n", path);
  } else {
    fprintf(fp, "sample %d\n", sample);
    fprintf(fp, "info_no_bf %d\n", infoNoBF);
    fprintf(fp, "info_bf %d\n", infoBF);
    fprintf(fp, "nprojbf %d\n", NProjBF);
    fprintf(fp, "nslater %d\n", NSlater);
    fprintf(fp, "step %.17e\n", h);
    fprintf(fp, "fd_fail_count %d\n", fdFailCount);
    fprintf(fp, "nan_count %d\n", nanCount);
    fprintf(fp, "nonzero_projbf_fd_count %d\n", nonzeroProjBFFDCount);
    fprintf(fp, "nonzero_orbital_fd_count %d\n", nonzeroOrbitalFDCount);
    fprintf(fp, "max_abs_identity_orbital_diff %.17e\n", maxIdentityOrbitalDiff);
    fprintf(fp, "max_identity_orbital_idx %d\n", maxIdentityOrbitalIdx);
    fprintf(fp, "max_abs_projbf_fd_real %.17e\n", maxProjBFRealDiff);
    fprintf(fp, "max_projbf_real_idx %d\n", maxProjBFRealIdx);
    fprintf(fp, "max_abs_projbf_fd_imag %.17e\n", maxProjBFImagDiff);
    fprintf(fp, "max_projbf_imag_idx %d\n", maxProjBFImagIdx);
    fprintf(fp, "max_abs_projbf_fd_diff %.17e\n", maxProjBFDiff);
    fprintf(fp, "max_projbf_idx %d\n", maxProjBFIdx);
    fprintf(fp, "max_projbf_is_imag %d\n", maxProjBFIsImag);
    fprintf(fp, "max_abs_orbital_fd_real %.17e\n", maxOrbitalRealDiff);
    fprintf(fp, "max_orbital_real_idx %d\n", maxOrbitalRealIdx);
    fprintf(fp, "max_abs_orbital_fd_imag %.17e\n", maxOrbitalImagDiff);
    fprintf(fp, "max_orbital_imag_idx %d\n", maxOrbitalImagIdx);
    fprintf(fp, "max_abs_orbital_fd_diff %.17e\n", maxOrbitalDiff);
    fprintf(fp, "max_orbital_idx %d\n", maxOrbitalIdx);
    fprintf(fp, "max_orbital_is_imag %d\n", maxOrbitalIsImag);
    fprintf(fp, "projbf_analytic_at_max %.17e %.17e\n",
            creal(projBFAnalyticAtMax), cimag(projBFAnalyticAtMax));
    fprintf(fp, "projbf_fd_at_max %.17e %.17e\n",
            creal(projBFFDAtMax), cimag(projBFFDAtMax));
    fprintf(fp, "orbital_analytic_at_max %.17e %.17e\n",
            creal(orbitalAnalyticAtMax), cimag(orbitalAnalyticAtMax));
    fprintf(fp, "orbital_fd_at_max %.17e %.17e\n",
            creal(orbitalFDAtMax), cimag(orbitalFDAtMax));
    fclose(fp);
  }

  free(diffNoBF);
  free(analyticPacked);
  free(projBFStore);
  free(slaterStore);
}

void VMCMainCal_fsz(MPI_Comm comm_parent, MPI_Comm comm) {
  int *eleIdx,*eleCfg,*eleNum,*eleProjCnt,*eleSpn; //fsz
  double complex e,ip;
  double w,x;
  double sqrtw;
  double complex we;
  double Sz;

  const int qpStart=0;
  const int qpEnd=NQPFull;
  int sample,sampleStart,sampleEnd,sampleSize;
  int i,info;

  /* optimazation for Kei */
  const int nProj=NProj;
  double complex *srOptO = SROptO;
  double         *srOptO_real = SROptO_real;
  int tmp_i;

  int rank,size,parentRank,parentSize,int_i;
  LSLanczosFSZScratch lanczosScratch;
  int lanczosScratchReady = 0;
  int lanczosWarnings = 0;
  long lanczosInjectSample = -1;
  long lanczosInjectParentRank = -1;
  long long lanczosCheckedSamples = 0;
  long long lanczosRejectedSamples = 0;
  FILE *lanczosOracleDump = NULL;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm_parent,&parentSize);
  MPI_Comm_rank(comm_parent,&parentRank);
#ifdef __DEBUG_DETAILDETAIL
  printf("  Debug: SplitLoop\n");
#endif
  SplitLoop(&sampleStart,&sampleEnd,NVMCSample,rank,size);

  /* initialization */
  StartTimer(24);
  clearPhysQuantity();
  StopTimer(24);

  memset(&lanczosScratch, 0, sizeof(lanczosScratch));
  if(NVMCCalMode == 1 && NLanczosMode == 1) {
    const char *dumpValue = getenv("MVMC_LANCZOS_ORACLE_DUMP");
    lanczosInjectSample = getLanczosNonfiniteInjectionSample();
    lanczosInjectParentRank = LSLanczosTestNonfiniteParentRank();
    if(lanczosInjectSample == -2) {
      if(rank == 0) {
        fprintf(stderr,
                "Error: MVMC_LANCZOS_TEST_NONFINITE_SAMPLE must be a non-negative integer.\n");
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
    if(lanczosInjectParentRank == -2 ||
       lanczosInjectParentRank >= parentSize) {
      if(parentRank == 0) {
        fprintf(stderr,
                "Error: MVMC_LANCZOS_TEST_NONFINITE_PARENT_RANK must identify a parent communicator rank.\n");
      }
      MPI_Abort(comm_parent, EXIT_FAILURE);
    }
    if(lanczosInjectParentRank >= 0 &&
       lanczosInjectParentRank != parentRank) {
      lanczosInjectSample = -1;
    }
    if(LSLanczosFSZScratchInit(&lanczosScratch, AllComplexFlag == 0) != 0) {
      if(rank == 0) {
        fprintf(stderr, "Error: failed to allocate FSZ Lanczos scratch.\n");
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
    lanczosScratchReady = 1;
    if(dumpValue != NULL && dumpValue[0] != '\0' && strcmp(dumpValue, "0") != 0) {
      char dumpPath[1024];
      int pathLength;
      const char *basePath = strcmp(dumpValue, "1") == 0
                           ? "lanczos_oracle_fsz.dat" : dumpValue;
      if(parentSize > 1) {
        pathLength = snprintf(dumpPath, sizeof(dumpPath), "%s.rank%04d",
                              basePath, parentRank);
      } else {
        pathLength = snprintf(dumpPath, sizeof(dumpPath), "%s", basePath);
      }
      if(pathLength < 0 || (size_t)pathLength >= sizeof(dumpPath)) {
        fprintf(stderr,
                "Error: FSZ Lanczos oracle dump path is too long on rank %d.\n",
                parentRank);
        MPI_Abort(comm_parent, EXIT_FAILURE);
      }
      lanczosOracleDump = fopen(dumpPath, "w");
      if(lanczosOracleDump == NULL) {
        fprintf(stderr,
                "Error: failed to open FSZ Lanczos oracle dump '%s' on rank %d.\n",
                dumpPath, parentRank);
        MPI_Abort(comm_parent, EXIT_FAILURE);
      }
    }
  }

  if(NVMCCalMode==0 && NStoreO!=0 && NSRCG==0) {
    clearStoredOSampleRange_fsz(sampleStart, sampleEnd);
  }

  for(sample=sampleStart;sample<sampleEnd;sample++) {
    eleIdx = EleIdx + sample*Nsize;
    eleCfg = EleCfg + sample*Nsite2;
    eleNum = EleNum + sample*Nsite2;
    eleProjCnt = EleProjCnt + sample*NProj;
    eleSpn     = EleSpn + sample*Nsize; //fsz

    StartTimer(40);
#ifdef _DEBUG_DETAIL
    printf("  Debug: sample=%d: CalculateMAll \n",sample);
#endif
    //info = CalculateMAll_fsz(eleIdx,eleSpn,qpStart,qpEnd);//info = CalculateMAll_fcmp(eleIdx,qpStart,qpEnd); // InvM,PfM will change
    if(AllComplexFlag==0){
      info = CalculateMAll_fsz_real(eleIdx,eleSpn,qpStart,qpEnd); // InvM_real,PfM_real will change
#pragma omp parallel for default(shared) private(tmp_i)
      for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++)  InvM[tmp_i]= InvM_real[tmp_i]; // InvM will be used in  SlaterElmDiff_fcmp
    }else{
      info = CalculateMAll_fsz(eleIdx,eleSpn,qpStart,qpEnd); // InvM,PfM will change
    }
    StopTimer(40);

    if(info!=0) {
      fprintf(stderr,"warning: VMCMainCal rank:%d sample:%d info:%d (CalculateMAll)\n",rank,sample,info);
      continue;
    }
#ifdef _DEBUG_DETAIL
    printf("  Debug: sample=%d: CalculateIP \n",sample);
#endif
    //ip = CalculateIP_fcmp(PfM,qpStart,qpEnd,MPI_COMM_SELF);
    if(AllComplexFlag==0){
      ip = CalculateIP_real(PfM_real,qpStart,qpEnd,MPI_COMM_SELF);
    }else{
      ip = CalculateIP_fcmp(PfM,qpStart,qpEnd,MPI_COMM_SELF);
    } 
#ifdef _DEBUG_DETAIL
    printf("  Debug: sample=%d: LogProjVal \n",sample);
#endif
    /* calculate reweight */
    x = LogProjVal(eleProjCnt);
    if (reweight==1){
       w = exp(2.0*(log(cabs(ip))+x) - logSqPfFullSlater[sample]);
    }else{
       w =1.0;
    }
    //LogProjVal(eleProjCnt);
    //w =1.0;
#ifdef _DEBUG_DETAIL
    printf("  Debug: sample=%d: isfinite \n",sample);
#endif
    if( !isfinite(w) ) {
      fprintf(stderr,"warning: VMCMainCal rank:%d sample:%d w=%e\n",rank,sample,w);
      continue;
    }

    StartTimer(41);
    /* calculate energy */
#ifdef _DEBUG_DETAIL
    printf("  Debug: sample=%d: calculateHam \n",sample);
#endif
#ifdef _DEBUG_DETAIL
    printf("  Debug: sample=%d: calculateHam_cmp \n",sample);
#endif
    //e  = CalculateHamiltonian_fsz(ip,eleIdx,eleCfg,eleNum,eleProjCnt,eleSpn);//fsz
    if(AllComplexFlag==0){
#ifdef _DEBUG_VMCCAL
      printf("  Debug: sample=%d: calculateHam_real \n",sample);
#endif
      e = CalculateHamiltonian_fsz_real(creal(ip),eleIdx,eleCfg,eleNum,eleProjCnt,eleSpn);
    }else{
#ifdef _DEBUG_VMCCAL
      printf("  Debug: sample=%d: calculateHam_cmp \n",sample);
#endif
      e = CalculateHamiltonian_fsz(ip,eleIdx,eleCfg,eleNum,eleProjCnt,eleSpn);
    }
    Sz = CalculateSz_fsz(ip,eleIdx,eleCfg,eleNum,eleProjCnt,eleSpn);//fsz
    //printf("MDEBUG: Sz=%lf \n",Sz);
    //printf("MDEBUG: e= %lf %lf ip= %lf %lf \n",creal(e),cimag(e),creal(ip),cimag(ip));
    StopTimer(41);
    if( !isfinite(creal(e) + cimag(e)) ) {
      fprintf(stderr,"warning: VMCMainCal rank:%d sample:%d e=%e\n",rank,sample,creal(e)); //TBC
      continue;
    }

    if(NVMCCalMode == 1 && NLanczosMode == 1) {
      int lanczosInfo;
      int lanczosFinite;
      lanczosCheckedSamples++;
      StartTimer(43);
      if(AllComplexFlag == 0) {
        lanczosInfo = LSLocalQ_fsz_real(creal(e), creal(ip), eleIdx,
                                        eleCfg, eleNum, eleProjCnt, eleSpn,
                                        &lanczosScratch, LSLQ_real);
      } else {
        lanczosInfo = LSLocalQ_fsz(e, ip, eleIdx, eleCfg, eleNum,
                                   eleProjCnt, eleSpn, &lanczosScratch,
                                   LSLQ);
      }
      StopTimer(43);
      if(lanczosInfo == LSLANCZOS_NUMERIC_REJECT) {
        lanczosRejectedSamples++;
        if(lanczosWarnings < 3) {
          fprintf(stderr,
                  "warning: FSZ Lanczos rejected numerical candidate "
                  "rebuild rank:%d sample:%d.\n", rank, sample);
          lanczosWarnings++;
        }
        continue;
      }
      if(lanczosInfo != LSLANCZOS_OK) {
        fprintf(stderr,
                "Error: FSZ Lanczos state or contract failure rank:%d "
                "sample:%d status:%d.\n", rank, sample, lanczosInfo);
        MPI_Abort(comm, EXIT_FAILURE);
      }
      if(AllComplexFlag == 0) {
        if((long)sample == lanczosInjectSample) LSLQ_real[3] = NAN;
        lanczosFinite = isFiniteLSLQFSZReal(LSLQ_real);
      } else {
        if((long)sample == lanczosInjectSample) LSLQ[3] = NAN + 0.0*I;
        lanczosFinite = isFiniteLSLQFSZComplex(LSLQ);
      }
      if(!lanczosFinite) {
        lanczosRejectedSamples++;
        if(lanczosWarnings < 3) {
          fprintf(stderr,
                  "warning: FSZ Lanczos rejected non-finite local vector rank:%d sample:%d.\n",
                  rank, sample);
          lanczosWarnings++;
        }
        continue;
      }
      if(lanczosOracleDump != NULL) {
        fprintf(lanczosOracleDump, "sample %d occ", sample);
        for(i=0; i<Nsite2; i++) fprintf(lanczosOracleDump, " %d", eleNum[i]);
        fprintf(lanczosOracleDump, " lslq");
        if(AllComplexFlag == 0) {
          for(i=0; i<NLSHam*NLSHam; i++) {
            fprintf(lanczosOracleDump, " %.17e %.17e", LSLQ_real[i], 0.0);
          }
        } else {
          for(i=0; i<NLSHam*NLSHam; i++) {
            fprintf(lanczosOracleDump, " %.17e %.17e",
                    creal(LSLQ[i]), cimag(LSLQ[i]));
          }
        }
        fprintf(lanczosOracleDump, "\n");
      }
    }

    Wc    += w;
    Etot  += w * e;
    Sztot += w * Sz;
    Sztot2 += w * Sz*Sz;
    Etot2 += w * conj(e) * e;
#ifdef _DEBUG_DETAIL
    printf("  Debug: sample=%d: calculateOpt \n",sample);
#endif
    if(NVMCCalMode==0) {
      /* Calculate O for correlation fauctors */
      srOptO[0] = 1.0+0.0*I;//   real 
      srOptO[1] = 0.0+0.0*I;//   real 
      #pragma loop noalias
      for(i=0;i<nProj;i++){ 
        srOptO[(i+1)*2]     = (double)(eleProjCnt[i]); // even real
        srOptO[(i+1)*2+1]   = 0.0+0.0*I;               // odd  comp
      }

      StartTimer(42);
      /* SlaterElmDiff */
      SlaterElmDiff_fsz(SROptO+2*NProj+2,ip,eleIdx,eleSpn) ;//SlaterElmDiff_fcmp(SROptO+2*NProj+2,ip,eleIdx); //TBC: using InvM not InvM_real
      StopTimer(42);
      
      if(FlagOptTrans>0) { // this part will be not used
        calculateOptTransDiff(SROptO+2*NProj+2*NSlater+2, ip); //TBC
      }
      //[s] this part will be used for real varaibles
      if(AllComplexFlag==0){
#pragma loop noalias
        for(i=0;i<SROptSize;i++){ 
          srOptO_real[i] = creal(srOptO[2*i]);       
        }
      }
      //[e]
      
      StartTimer(43);
      /* Calculate OO and HO */
      if(NSRCG==0 && NStoreO==0){
        //calculateOO(SROptOO,SROptHO,SROptO,w,e,SROptSize);
        if(AllComplexFlag==0){
          calculateOO_real(SROptOO_real,SROptHO_real,SROptO_real,w,creal(e),SROptSize);
        }else{
          calculateOO(SROptOO,SROptHO,SROptO,w,e,SROptSize);
        } 
      }else{
        we    = w*e;
        sqrtw = sqrt(w); 
        /*#pragma omp parallel for default(shared) private(int_i)
        for(int_i=0;int_i<SROptSize*2;int_i++){
        // SROptO_Store for fortran
          SROptO_Store[int_i+sample*(2*SROptSize)]  = sqrtw*SROptO[int_i];
          SROptHO[int_i]                           += we*SROptO[int_i]; 
        }*/
        if(AllComplexFlag==0){
          #pragma omp parallel for default(shared) private(int_i)
          for(int_i=0;int_i<SROptSize;int_i++){
            // SROptO_Store for fortran
            SROptO_Store_real[int_i+sample*SROptSize]  = sqrtw*SROptO_real[int_i];
            SROptHO_real[int_i]                       += creal(we)*SROptO_real[int_i]; 
          }
        }else{
          #pragma omp parallel for default(shared) private(int_i)
          for(int_i=0;int_i<SROptSize*2;int_i++){
            // SROptO_Store for fortran
            SROptO_Store[int_i+sample*(2*SROptSize)]  = sqrtw*SROptO[int_i];
            SROptHO[int_i]                           += we*SROptO[int_i]; 
          }
        }
      } 
      StopTimer(43);
    } else if(NVMCCalMode==1) {
      StartTimer(42);
      /* Calculate Green Function */
      CalculateGreenFunc_fsz(w,ip,eleIdx,eleCfg,eleNum,eleSpn,eleProjCnt);
      StopTimer(42);

      if(NLanczosMode>0){
        if(AllComplexFlag == 0) {
          calculateQQQQ_real(QQQQ_real, LSLQ_real, w, NLSHam);
        } else {
          calculateQQQQ(QQQQ, LSLQ, w, NLSHam);
        }
      }
    }
  } /* end of for(sample) */

  if(NVMCCalMode == 1 && NLanczosMode == 1) {
    finalizeLanczosFSZAccounting(comm_parent, parentRank, lanczosCheckedSamples,
                                 lanczosRejectedSamples);
  }
  if(lanczosOracleDump != NULL) fclose(lanczosOracleDump);
  if(lanczosScratchReady) LSLanczosFSZScratchFree(&lanczosScratch);

// calculate OO and HO at NVMCCalMode==0
  if(NVMCCalMode==0){
    if(NStoreO!=0 || NSRCG!=0){
      sampleSize=sampleEnd-sampleStart;
      if(NSRCG!=0 || sampleSize>0){
        /*StartTimer(45);
        calculateOO_Store(SROptOO,SROptHO,SROptO_Store,w,e,2*SROptSize,sampleSize);
        StopTimer(45);*/
        if(AllComplexFlag==0){
          double *srOptO_Store_ptr = SROptO_Store_real;
          if(NSRCG==0) {
            srOptO_Store_ptr += (size_t)sampleStart * (size_t)SROptSize;
          }
          StartTimer(45);
          calculateOO_Store_real(SROptOO_real,SROptHO_real,srOptO_Store_ptr,creal(w),creal(e),SROptSize,sampleSize);
          StopTimer(45);
        }else{
          double complex *srOptO_Store_ptr = SROptO_Store;
          if(NSRCG==0) {
            srOptO_Store_ptr += (size_t)sampleStart * (2 * (size_t)SROptSize);
          }
          StartTimer(45);
          calculateOO_Store(SROptOO,SROptHO,srOptO_Store_ptr,w,e,2*SROptSize,sampleSize);
          StopTimer(45);
        }
      }
    }
  }
  return;
}

void VMC_BF_MainCal_fsz(MPI_Comm comm_parent, MPI_Comm comm) {
  int *eleIdx,*eleCfg,*eleNum,*eleProjCnt,*eleSpn,*eleProjBFCnt;
  double complex e,ip;
  double w,x;
  double sqrtw;
  double complex we;
  double Sz;
  double c2SzStart = 0.0;
  const char *bfFSZSRDiffDumpPath = getenv("MVMC_BF_FSZ_SR_DIFF_DUMP");

  const int qpStart=0;
  const int qpEnd=NQPFull;
  int sample,sampleStart,sampleEnd,sampleSize;
  int i,info;
  int int_i;
  const int nProj=NProj;
  double complex *srOptO = SROptO;
  double         *srOptO_real = SROptO_real;

  int rank,size,parentRank,parentSize;
  LSLanczosBFFSZScratch lanczosScratch;
  int lanczosScratchReady = 0;
  int lanczosWarnings = 0;
  long lanczosInjectSample = -1;
  long lanczosInjectParentRank = -1;
  long long lanczosCheckedSamples = 0;
  long long lanczosRejectedSamples = 0;
  FILE *lanczosOracleDump = NULL;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm_parent,&parentSize);
  MPI_Comm_rank(comm_parent,&parentRank);
  SplitLoop(&sampleStart,&sampleEnd,NVMCSample,rank,size);

  StartTimer(24);
  clearPhysQuantity();
  StopTimer(24);

  memset(&lanczosScratch, 0, sizeof(lanczosScratch));
  if(NVMCCalMode == 1 && NLanczosMode == 1) {
    const char *dumpValue = getenv("MVMC_LANCZOS_ORACLE_DUMP");
    lanczosInjectSample = getLanczosNonfiniteInjectionSample();
    lanczosInjectParentRank = LSLanczosTestNonfiniteParentRank();
    if(lanczosInjectSample == -2) {
      if(rank == 0) {
        fprintf(stderr,
                "Error: MVMC_LANCZOS_TEST_NONFINITE_SAMPLE must be a non-negative integer.\n");
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
    if(lanczosInjectParentRank == -2 ||
       lanczosInjectParentRank >= parentSize) {
      if(parentRank == 0) {
        fprintf(stderr,
                "Error: MVMC_LANCZOS_TEST_NONFINITE_PARENT_RANK must identify a parent communicator rank.\n");
      }
      MPI_Abort(comm_parent, EXIT_FAILURE);
    }
    if(lanczosInjectParentRank >= 0 &&
       lanczosInjectParentRank != parentRank) {
      lanczosInjectSample = -1;
    }
    if(LSLanczosBFFSZScratchInit(&lanczosScratch) != 0) {
      if(rank == 0) {
        fprintf(stderr, "Error: failed to allocate BF-FSZ Lanczos scratch.\n");
      }
      MPI_Abort(comm, EXIT_FAILURE);
    }
    lanczosScratchReady = 1;
    if(dumpValue != NULL && dumpValue[0] != '\0' &&
       strcmp(dumpValue, "0") != 0) {
      char dumpPath[1024];
      int pathLength;
      const char *basePath = strcmp(dumpValue, "1") == 0
                           ? "lanczos_oracle_bf_fsz.dat" : dumpValue;
      if(parentSize > 1) {
        pathLength = snprintf(dumpPath, sizeof(dumpPath), "%s.rank%04d",
                              basePath, parentRank);
      } else {
        pathLength = snprintf(dumpPath, sizeof(dumpPath), "%s", basePath);
      }
      if(pathLength < 0 || (size_t)pathLength >= sizeof(dumpPath)) {
        fprintf(stderr,
                "Error: BF-FSZ Lanczos oracle dump path is too long on rank %d.\n",
                parentRank);
        MPI_Abort(comm_parent, EXIT_FAILURE);
      }
      lanczosOracleDump = fopen(dumpPath, "w");
      if(lanczosOracleDump == NULL) {
        fprintf(stderr,
                "Error: failed to open BF-FSZ Lanczos oracle dump '%s' on rank %d.\n",
                dumpPath, parentRank);
        MPI_Abort(comm_parent, EXIT_FAILURE);
      }
    }
  }

  if(NVMCCalMode==0 && NStoreO!=0 && NSRCG==0) {
    clearStoredOSampleRange_fsz(sampleStart, sampleEnd);
  }

  for(sample=sampleStart;sample<sampleEnd;sample++) {
    eleIdx = EleIdx + sample*Nsize;
    eleCfg = EleCfg + sample*Nsite2;
    eleNum = EleNum + sample*Nsite2;
    eleProjCnt = EleProjCnt + sample*NProj;
    eleProjBFCnt = EleProjBFCnt + sample*16*Nsite*Nrange;
    eleSpn = EleSpn + sample*Nsize;

    StartTimer(40);
    MakeSlaterElmBF_fsz(eleNum,eleProjBFCnt);
    if(NVMCCalMode == 1 && rank == 0 && sample == sampleStart &&
        bfFSZSRDiffDumpPath != NULL && bfFSZSRDiffDumpPath[0] != '\0') {
      dumpBFFSZSRDiffCheck(bfFSZSRDiffDumpPath, eleIdx, eleSpn, eleNum,
                           eleProjBFCnt, qpStart, qpEnd, sample);
    }
    info = CalculateMAll_BF_fsz(eleIdx,eleSpn,qpStart,qpEnd);
    StopTimer(40);

    if(info!=0) {
      fprintf(stderr,"warning: VMC_BF_MainCal_fsz rank:%d sample:%d info:%d (CalculateMAll)\n",rank,sample,info);
      continue;
    }

    ip = CalculateIP_fcmp(PfM,qpStart,qpEnd,MPI_COMM_SELF);
    if(!isfinite(creal(ip)) || !isfinite(cimag(ip))) {
      fprintf(stderr,"warning: VMC_BF_MainCal_fsz rank:%d sample:%d ip=%e %e\n",
              rank,sample,creal(ip),cimag(ip));
      continue;
    }

    x = LogProjVal(eleProjCnt);
    if(reweight==1){
      w = exp(2.0*(log(cabs(ip))+x) - logSqPfFullSlater[sample]);
    }else{
      w = 1.0;
    }

    if(!isfinite(w)) {
      fprintf(stderr,"warning: VMC_BF_MainCal_fsz rank:%d sample:%d w=%e\n",rank,sample,w);
      continue;
    }

    StartTimer(41);
    e = CalculateHamiltonianBF_fsz(ip,eleIdx,eleCfg,eleNum,eleProjCnt,eleSpn,eleProjBFCnt);
    if(BFFSZC2DetailProfileEnabled) {
      c2SzStart = BFFSZC2DetailMonotonicSeconds();
    }
    Sz = CalculateSz_fsz(ip,eleIdx,eleCfg,eleNum,eleProjCnt,eleSpn);
    if(BFFSZC2DetailProfileEnabled) {
      BFFSZC2DetailSzSecondsRank0 += BFFSZC2DetailMonotonicSeconds()-c2SzStart;
    }
    StopTimer(41);

    if(!isfinite(creal(e) + cimag(e))) {
      fprintf(stderr,"warning: VMC_BF_MainCal_fsz rank:%d sample:%d e=%e\n",rank,sample,creal(e));
      continue;
    }

    if(NVMCCalMode == 1 && NLanczosMode == 1) {
      int lanczosInfo;
      int lanczosFinite;
      lanczosCheckedSamples++;
      StartTimer(43);
      if(AllComplexFlag == 0) {
        lanczosInfo = LSLocalQBF_fsz_real(
            e, ip, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn,
            eleProjBFCnt, &lanczosScratch, LSLQ_real);
      } else {
        lanczosInfo = LSLocalQBF_fsz(
            e, ip, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn,
            eleProjBFCnt, &lanczosScratch, LSLQ);
      }
      StopTimer(43);
      if(lanczosInfo == LSLANCZOS_REAL_GATE_FAILURE) {
        fprintf(stderr,
                "Error: BF-FSZ real Lanczos vector has a significant imaginary "
                "component rank:%d sample:%d index:%d value:(%.17e,%.17e) "
                "tolerance:%.17e.\n", rank, sample,
                lanczosScratch.realGateIndex,
                creal(lanczosScratch.realGateValue),
                cimag(lanczosScratch.realGateValue),
                lanczosScratch.realGateTolerance);
        MPI_Abort(comm, EXIT_FAILURE);
      }
      if(lanczosInfo == LSLANCZOS_NUMERIC_REJECT) {
        lanczosRejectedSamples++;
        if(lanczosWarnings < 3) {
          fprintf(stderr,
                  "warning: BF-FSZ Lanczos rejected numerical candidate "
                  "rebuild rank:%d sample:%d.\n", rank, sample);
          lanczosWarnings++;
        }
        continue;
      }
      if(lanczosInfo != LSLANCZOS_OK) {
        fprintf(stderr,
                "Error: BF-FSZ Lanczos state or contract failure rank:%d "
                "sample:%d status:%d.\n", rank, sample, lanczosInfo);
        MPI_Abort(comm, EXIT_FAILURE);
      }
      if((long)sample == lanczosInjectSample) {
        if(AllComplexFlag == 0) {
          LSLQ_real[3] = NAN;
        } else {
          LSLQ[3] = NAN + 0.0*I;
        }
      }
      lanczosFinite = AllComplexFlag == 0
                    ? isFiniteLSLQFSZReal(LSLQ_real)
                    : isFiniteLSLQFSZComplex(LSLQ);
      if(!lanczosFinite) {
        lanczosRejectedSamples++;
        if(lanczosWarnings < 3) {
          fprintf(stderr,
                  "warning: BF-FSZ Lanczos rejected non-finite local vector "
                  "rank:%d sample:%d.\n", rank, sample);
          lanczosWarnings++;
        }
        continue;
      }
      if(lanczosOracleDump != NULL) {
        fprintf(lanczosOracleDump, "sample %d occ", sample);
        for(i=0; i<Nsite2; i++) fprintf(lanczosOracleDump, " %d", eleNum[i]);
        fprintf(lanczosOracleDump, " lslq");
        if(AllComplexFlag == 0) {
          for(i=0; i<NLSHam*NLSHam; i++) {
            fprintf(lanczosOracleDump, " %.17e %.17e", LSLQ_real[i], 0.0);
          }
        } else {
          for(i=0; i<NLSHam*NLSHam; i++) {
            fprintf(lanczosOracleDump, " %.17e %.17e",
                    creal(LSLQ[i]), cimag(LSLQ[i]));
          }
        }
        fprintf(lanczosOracleDump, "\n");
      }
    }

    Wc    += w;
    Etot  += w * e;
    Sztot += w * Sz;
    Sztot2 += w * Sz*Sz;
    Etot2 += w * conj(e) * e;

    if(NVMCCalMode==0) {
      srOptO[0] = 1.0+0.0*I;
      srOptO[1] = 0.0+0.0*I;
#pragma loop noalias
      for(i=0;i<nProj;i++){
        srOptO[(i+1)*2]   = (double)(eleProjCnt[i]);
        srOptO[(i+1)*2+1] = 0.0+0.0*I;
      }

      BackFlowDiff_fsz(SROptO+2*NProj+2,ip,eleIdx,eleSpn,eleNum,eleProjBFCnt);

      StartTimer(42);
      SlaterElmBFDiff_fsz(SROptO+2*NProj+2*NProjBF+2,ip,eleIdx,eleSpn,eleNum,eleProjBFCnt);
      StopTimer(42);

      if(FlagOptTrans>0) {
        calculateOptTransDiff(SROptO+2*NProj+2*NProjBF+2*NSlater+2, ip);
      }

      if(AllComplexFlag==0){
#pragma loop noalias
        for(i=0;i<SROptSize;i++){
          srOptO_real[i] = creal(srOptO[2*i]);
        }
      }

      StartTimer(43);
      if(NSRCG==0 && NStoreO==0){
        if(AllComplexFlag==0){
          calculateOO_real(SROptOO_real,SROptHO_real,SROptO_real,w,creal(e),SROptSize);
        }else{
          calculateOO(SROptOO,SROptHO,SROptO,w,e,SROptSize);
        }
      }else{
        we    = w*e;
        sqrtw = sqrt(w);
        if(AllComplexFlag==0){
#pragma omp parallel for default(shared) private(int_i)
          for(int_i=0;int_i<SROptSize;int_i++){
            SROptO_Store_real[int_i+sample*SROptSize] = sqrtw*SROptO_real[int_i];
            SROptHO_real[int_i]                      += creal(we)*SROptO_real[int_i];
          }
        }else{
#pragma omp parallel for default(shared) private(int_i)
          for(int_i=0;int_i<SROptSize*2;int_i++){
            SROptO_Store[int_i+sample*(2*SROptSize)] = sqrtw*SROptO[int_i];
            SROptHO[int_i]                          += we*SROptO[int_i];
          }
        }
      }
      StopTimer(43);
    } else if(NVMCCalMode==1) {
      StartTimer(42);
      CalculateGreenFuncBF_fsz(w,ip,eleIdx,eleCfg,eleNum,eleSpn,eleProjCnt,eleProjBFCnt);
      StopTimer(42);

      if(NLanczosMode>0) {
        if(AllComplexFlag == 0) {
          calculateQQQQ_real(QQQQ_real, LSLQ_real, w, NLSHam);
        } else {
          calculateQQQQ(QQQQ, LSLQ, w, NLSHam);
        }
      }
    }
  }

  if(NVMCCalMode == 1 && NLanczosMode == 1) {
    finalizeLanczosFSZAccounting(comm_parent, parentRank, lanczosCheckedSamples,
                                 lanczosRejectedSamples);
  }
  if(lanczosOracleDump != NULL) fclose(lanczosOracleDump);
  if(lanczosScratchReady) LSLanczosBFFSZScratchFree(&lanczosScratch);

  if(NVMCCalMode==0){
    if(NStoreO!=0 || NSRCG!=0){
      sampleSize=sampleEnd-sampleStart;
      if(NSRCG!=0 || sampleSize>0){
        if(AllComplexFlag==0){
          double *srOptO_Store_ptr = SROptO_Store_real;
          if(NSRCG==0) {
            srOptO_Store_ptr += (size_t)sampleStart * (size_t)SROptSize;
          }
          StartTimer(45);
          calculateOO_Store_real(SROptOO_real,SROptHO_real,srOptO_Store_ptr,creal(w),creal(e),SROptSize,sampleSize);
          StopTimer(45);
        }else{
          double complex *srOptO_Store_ptr = SROptO_Store;
          if(NSRCG==0) {
            srOptO_Store_ptr += (size_t)sampleStart * (2 * (size_t)SROptSize);
          }
          StartTimer(45);
          calculateOO_Store(SROptOO,SROptHO,srOptO_Store_ptr,w,e,2*SROptSize,sampleSize);
          StopTimer(45);
        }
      }
    }
  }
  return;
}

#endif
