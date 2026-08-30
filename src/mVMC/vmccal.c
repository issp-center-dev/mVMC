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
 * calculate physical quantities
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#ifndef _SRC_VMCCAL
#define _SRC_VMCCAL

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "vmccal.h"
#include "matrix.h"
#include "calham_real.h"
#include "calham.h"
#include "lslocgrn_real.c"
#include "lslocgrn.c"
#include "lslocgrn_heisenberg.c"
#include "calgrn.c"
#include "physcal_lanczos2.h"

//#define _DEBUG_VMCCAL
//#define _DEBUG_VMCCAL_DETAIL

void clearPhysQuantity();

void calculateOptTransDiff(double complex *srOptO, const double complex ipAll);
void calculateOO_matvec(double complex *srOptOO, double complex *srOptHO, const double complex *srOptO,
                 const double complex w, const double complex e, const int srOptSize);
void calculateOO(double complex *srOptOO, double complex *srOptHO, const double complex *srOptO,
                 const double  w, const double complex e, const int srOptSize);
void calculateOO_real(double *srOptOO, double *srOptHO, const double *srOptO,
                 const double w, const double e, const int srOptSize);
void calculateOO_Store_real(double *srOptOO_real, double *srOptHO_real,  double *srOptO_real,
                 const double w, const double e,  int srOptSize, int sampleSize);
void calculateOO_Store(double complex *srOptOO, double complex *srOptHO,  double complex *srOptO,
                 const double w, const double complex e,  int srOptSize, int sampleSize);

void calculateQQQQ_real(double *qqqq, const double *lslq, const double w, const int nLSHam);

void calculateQQQQ(double complex *qqqq, const double complex*lslq, const double w, const int nLSHam);
void calculateQCAQ(double complex *qcaq, const double complex*lslca, const double complex*lslq,
                   const double w, const int nLSHam, const int nCA);
void calculateQCACAQ(double complex *qcacaq, const double complex*lslca, const double w,
                     const int nLSHam, const int nCA, const int nCACA,
                     int **cacaIdx);

void calculateQCAQ_real(double *qcaq, const double *lslca, const double *lslq,
                   const double w, const int nLSHam, const int nCA);

void calculateQCACAQ_real(double *qcacaq, const double *lslca, const double w,
                          const int nLSHam, const int nCA, const int nCACA,
                          int **cacaIdx);

void calculateQCACAQDC_real(double *qcacaq, const double *lslq, const double w,
                            const int nLSHam, const int nCA, const int nCACA,
                            int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,
                            const double h1, const double ip);

void calculateQCACAQDC(double complex *qcacaq, const double complex *lslq, const double w,
                       const int nLSHam, const int nCA, const int nCACA,
                       int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,
                       const double complex h1, const double complex ip,double complex *rbmCnt);

static void RecordPowerLanczosSupportLSLQSample(
    const double *lslqReal, const double complex *lslqComplex,
    double weight, int parentRank, int sample, MPI_Comm comm) {
  Lanczos2SolveStatus status;
  const size_t offset =
      (size_t)sample * POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH;
  if(AllComplexFlag == 0) {
    const double power[LANCZOS2_POWER_COUNT] = {
        lslqReal[0],lslqReal[1],lslqReal[3],0.0};
    status = RecordPowerLanczosSupportSampleReal(
        PowerLanczosSupportSampleData + offset,power,3,weight);
  }else{
    const double complex power[LANCZOS2_POWER_COUNT] = {
        lslqComplex[0],lslqComplex[1],lslqComplex[3],0.0};
    status = RecordPowerLanczosSupportSampleComplex(
        PowerLanczosSupportSampleData + offset,power,3,weight);
  }
  if(status != LANCZOS2_SOLVE_OK) {
    fprintf(stderr,
            "Error: power-Lanczos support sample failure on rank %d, "
            "sample %d: %s.\n",
            parentRank,sample,Lanczos2SolveError(status));
    MPI_Abort(comm,EXIT_FAILURE);
  }
}

static void clearStoredOSampleRange(const int sampleStart, const int sampleEnd) {
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

typedef struct {
  int enabled;
  int *discrete;
  double *matrixReal;
  double complex *matrixComplex;
  size_t discreteCount;
  size_t matrixCount;
} Lanczos2StateSnapshot;

static int Lanczos2StateSnapshotInit(Lanczos2StateSnapshot *snapshot) {
#ifndef MVMC_ENABLE_FAULT_INJECTION
  memset(snapshot, 0, sizeof(*snapshot));
  return 0;
#else
  const char *value = getenv("MVMC_LANCZOS2_STATE_CHECK");
  size_t matrixStride;
  size_t matrixSize;
  memset(snapshot, 0, sizeof(*snapshot));
  if (value == NULL || value[0] == '\0' || strcmp(value, "0") == 0) {
    return 0;
  }
  snapshot->enabled = 1;
  if (Nsize < 0 || Nsite2 < 0 || NProj < 0 || NQPFull < 0 ||
      (size_t)Nsite2 >
          (SIZE_MAX - (size_t)Nsize - (size_t)NProj) / 2) {
    return -1;
  }
  snapshot->discreteCount =
      (size_t)Nsize + 2 * (size_t)Nsite2 + (size_t)NProj;
  if (Nsize != 0 &&
      (size_t)Nsize > SIZE_MAX / (size_t)Nsize) return -1;
  matrixStride = (size_t)Nsize * (size_t)Nsize + 1;
  if (matrixStride != 0 &&
      (size_t)NQPFull > SIZE_MAX / matrixStride) return -1;
  snapshot->matrixCount = (size_t)NQPFull * matrixStride;
  if (snapshot->discreteCount > SIZE_MAX / sizeof(*snapshot->discrete)) {
    return -1;
  }
  snapshot->discrete =
      (int *)malloc(snapshot->discreteCount * sizeof(*snapshot->discrete));
  matrixSize = AllComplexFlag == 0 ? sizeof(*snapshot->matrixReal)
                                   : sizeof(*snapshot->matrixComplex);
  if (snapshot->matrixCount > SIZE_MAX / matrixSize) return -1;
  if (AllComplexFlag == 0) {
    snapshot->matrixReal =
        (double *)malloc(snapshot->matrixCount * matrixSize);
  } else {
    snapshot->matrixComplex =
        (double complex *)malloc(snapshot->matrixCount * matrixSize);
  }
  if (snapshot->discrete == NULL ||
      (AllComplexFlag == 0 && snapshot->matrixReal == NULL) ||
      (AllComplexFlag != 0 && snapshot->matrixComplex == NULL)) {
    free(snapshot->discrete);
    free(snapshot->matrixReal);
    free(snapshot->matrixComplex);
    memset(snapshot, 0, sizeof(*snapshot));
    return -1;
  }
  return 0;
#endif
}

static void Lanczos2StateSnapshotCapture(
    Lanczos2StateSnapshot *snapshot, const int *eleIdx, const int *eleCfg,
    const int *eleNum, const int *eleProjCnt) {
  int *destination;
  if (!snapshot->enabled) return;
  destination = snapshot->discrete;
  memcpy(destination, eleIdx, (size_t)Nsize * sizeof(*destination));
  destination += Nsize;
  memcpy(destination, eleCfg, (size_t)Nsite2 * sizeof(*destination));
  destination += Nsite2;
  memcpy(destination, eleNum, (size_t)Nsite2 * sizeof(*destination));
  destination += Nsite2;
  if (NProj > 0) {
    memcpy(destination, eleProjCnt, (size_t)NProj * sizeof(*destination));
  }
  if (AllComplexFlag == 0) {
    memcpy(snapshot->matrixReal, InvM_real,
           snapshot->matrixCount * sizeof(*snapshot->matrixReal));
  } else {
    memcpy(snapshot->matrixComplex, InvM,
           snapshot->matrixCount * sizeof(*snapshot->matrixComplex));
  }
}

static int Lanczos2StateSnapshotMatches(
    const Lanczos2StateSnapshot *snapshot, const int *eleIdx,
    const int *eleCfg, const int *eleNum, const int *eleProjCnt) {
  const int *source;
  if (!snapshot->enabled) return 1;
  source = snapshot->discrete;
  if (memcmp(source, eleIdx, (size_t)Nsize * sizeof(*source)) != 0) return 0;
  source += Nsize;
  if (memcmp(source, eleCfg, (size_t)Nsite2 * sizeof(*source)) != 0) return 0;
  source += Nsite2;
  if (memcmp(source, eleNum, (size_t)Nsite2 * sizeof(*source)) != 0) return 0;
  source += Nsite2;
  if (NProj > 0 &&
      memcmp(source, eleProjCnt, (size_t)NProj * sizeof(*source)) != 0) {
    return 0;
  }
  if (AllComplexFlag == 0) {
    return memcmp(snapshot->matrixReal, InvM_real,
                  snapshot->matrixCount *
                      sizeof(*snapshot->matrixReal)) == 0;
  }
  return memcmp(snapshot->matrixComplex, InvM,
                snapshot->matrixCount *
                    sizeof(*snapshot->matrixComplex)) == 0;
}

static void Lanczos2StateSnapshotFree(Lanczos2StateSnapshot *snapshot) {
  free(snapshot->discrete);
  free(snapshot->matrixReal);
  free(snapshot->matrixComplex);
  memset(snapshot, 0, sizeof(*snapshot));
}

void VMCMainCal(MPI_Comm comm_parent, MPI_Comm comm) {
  int *eleIdx,*eleCfg,*eleNum,*eleProjCnt;
  double complex e,ip;
  double x,w;
  double sqrtw;
  double complex we;

  const int qpStart=0;
  const int qpEnd=NQPFull;
  int sample,sampleStart,sampleEnd,sampleSize;
  int i,info,tmp_i;

  /* optimazation for Kei */
  const int nProj=NProj;
  double complex *srOptO = SROptO;
  double         *srOptO_real = SROptO_real;

  double complex *rbmCnt;
  const int nSizeRBM=NRBM_PhysLayerIdx + Nneuron;

  int rank,size,parentRank,parentSize,int_i;
  FILE *lanczosOracleDump = NULL;
  Lanczos2StateSnapshot lanczos2StateSnapshot;
  Lanczos2HeisenbergScratch lanczos2HeisenbergScratch;
#ifdef MVMC_ENABLE_FAULT_INJECTION
  int lanczos2InjectNonfinite = 0;
  int lanczos2CalHCAGuardAudit = 0;
#endif
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm_parent,&parentSize);
  MPI_Comm_rank(comm_parent,&parentRank);
#ifdef _DEBUG_VMCCAL
  printf("  Debug: SplitLoop\n");
#endif
  SplitLoop(&sampleStart,&sampleEnd,NVMCSample,rank,size);

  /* initialization */
  StartTimer(24);
  clearPhysQuantity();
  StopTimer(24);

  memset(&lanczos2StateSnapshot, 0, sizeof(lanczos2StateSnapshot));
  memset(&lanczos2HeisenbergScratch, 0, sizeof(lanczos2HeisenbergScratch));
#ifdef MVMC_ENABLE_FAULT_INJECTION
  /*
   * Test-only hook for exercising the collective non-finite failure path.
   * Production builds (Testing=OFF) do not compile this mutation.
   */
  {
    const char *injectValue =
        getenv("MVMC_LANCZOS2_TEST_NONFINITE_SAMPLE");
    lanczos2InjectNonfinite =
        injectValue != NULL && injectValue[0] != '\0' &&
        strcmp(injectValue, "0") != 0;
  }
  {
    const char *auditValue =
        getenv("MVMC_LANCZOS2_TEST_CALHCA_GUARD_AUDIT");
    lanczos2CalHCAGuardAudit =
        auditValue != NULL && auditValue[0] != '\0' &&
        strcmp(auditValue, "0") != 0;
    if(lanczos2CalHCAGuardAudit &&
       (parentSize != 1 || NVMCCalMode != 1 || NLanczosMode <= 0 ||
        NLanczosEstimatorMode != 1 ||
        NLanczosStep != 2 || NQPFull <= 1)) {
      fprintf(stderr,
              "Error: Lanczos2 calHCA guard audit requires single-rank "
              "step=2 physical calculation with NQPFull>1.\n");
      MPI_Abort(comm_parent, EXIT_FAILURE);
    }
    Lanczos2CalHCAGuardAuditRealDirectCount = 0;
    Lanczos2CalHCAGuardAuditRealZeroComponentCount = 0;
    Lanczos2CalHCAGuardAuditComplexDirectCount = 0;
    Lanczos2CalHCAGuardAuditComplexZeroComponentCount = 0;
    Lanczos2CalHCAGuardAuditRealEnabled = lanczos2CalHCAGuardAudit;
    Lanczos2CalHCAGuardAuditComplexEnabled = lanczos2CalHCAGuardAudit;
  }
#endif
  if(NVMCCalMode == 1 && NLanczosMode > 0 &&
     NLanczosEstimatorMode == 1 && NLanczosStep == 2 &&
     Lanczos2StateSnapshotInit(&lanczos2StateSnapshot) != 0) {
    fprintf(stderr,
            "Error: failed to allocate Lanczos2 state snapshot on rank %d.\n",
            parentRank);
    MPI_Abort(comm_parent, EXIT_FAILURE);
  }
  if(NVMCCalMode == 1 && NLanczosMode > 0 &&
     NLanczosEstimatorMode == 1 && NLanczosStep == 2 &&
     NExUpdatePath == 2 &&
     Lanczos2HeisenbergScratchInit(
         &lanczos2HeisenbergScratch, AllComplexFlag != 0) != 0) {
    fprintf(stderr,
            "Error: failed to allocate Lanczos2 Heisenberg scratch "
            "on rank %d.\n",
            parentRank);
    MPI_Abort(comm_parent, EXIT_FAILURE);
  }
#ifdef MVMC_ENABLE_FAULT_INJECTION
  if(lanczos2HeisenbergScratch.auditEnabled) {
    if(parentSize != 1 || NVMCCalMode != 1 || NLanczosMode != 1 ||
       NLanczosEstimatorMode != 1 ||
       NLanczosStep != 2 || NExUpdatePath != 2 || NQPFull != 1 ||
       NExchangeCoupling <= 0) {
      fprintf(stderr,
              "Error: Lanczos2 Heisenberg audit requires a single-rank "
              "pure-spin Exchange calculation with NQPFull=1.\n");
      MPI_Abort(comm_parent, EXIT_FAILURE);
    }
    Lanczos2HeisenbergZeroOverlapAuditRealCount = 0;
    Lanczos2HeisenbergZeroOverlapAuditComplexCount = 0;
    Lanczos2HeisenbergZeroOverlapAuditRealEnabled =
        AllComplexFlag == 0;
    Lanczos2HeisenbergZeroOverlapAuditComplexEnabled =
        AllComplexFlag != 0;
  }
#endif

  if(NVMCCalMode == 1 && NLanczosMode > 0 &&
     NLanczosEstimatorMode == 1) {
#ifdef MVMC_ENABLE_FAULT_INJECTION
    const char *dumpValue = getenv("MVMC_LANCZOS_ORACLE_DUMP");
    if(dumpValue != NULL && dumpValue[0] != '\0' && strcmp(dumpValue, "0") != 0) {
      char dumpPath[1024];
      int pathLength;
      const char *basePath = strcmp(dumpValue, "1") == 0
                           ? "lanczos_oracle_nonfsz.dat" : dumpValue;
      if(parentSize > 1) {
        pathLength = snprintf(dumpPath, sizeof(dumpPath), "%s.rank%04d",
                              basePath, parentRank);
      } else {
        pathLength = snprintf(dumpPath, sizeof(dumpPath), "%s", basePath);
      }
      if(pathLength < 0 || (size_t)pathLength >= sizeof(dumpPath)) {
        fprintf(stderr,
                "Error: non-FSZ Lanczos oracle dump path is too long on rank %d.\n",
                parentRank);
        MPI_Abort(comm_parent, EXIT_FAILURE);
      }
      lanczosOracleDump = fopen(dumpPath, "w");
      if(lanczosOracleDump == NULL) {
        fprintf(stderr,
                "Error: failed to open non-FSZ Lanczos oracle dump '%s' on rank %d.\n",
                dumpPath, parentRank);
        MPI_Abort(comm_parent, EXIT_FAILURE);
      }
    }
#endif
  }

  if(NVMCCalMode==0 && NStoreO!=0 && NSRCG==0) {
    clearStoredOSampleRange(sampleStart, sampleEnd);
  }

  for(sample=sampleStart;sample<sampleEnd;sample++) {

    eleIdx = EleIdx + sample*Nsize;
    eleCfg = EleCfg + sample*Nsite2;
    eleNum = EleNum + sample*Nsite2;
    eleProjCnt = NProj > 0 ? EleProjCnt + sample*NProj : NULL;
    rbmCnt = FlagRBM ? RBMCnt + sample*nSizeRBM : NULL;

    StartTimer(40);
#ifdef _DEBUG_VMCCAL
    printf("  Debug: sample=%d: CalculateMAll \n",sample);
#endif
    if(AllComplexFlag==0){
      info = CalculateMAll_real(eleIdx,qpStart,qpEnd); // InvM_real,PfM_real will change
#pragma omp parallel for default(shared) private(tmp_i)
      for(tmp_i=0;tmp_i<NQPFull*(Nsize*Nsize+1);tmp_i++)  InvM[tmp_i]= InvM_real[tmp_i]; // InvM will be used in  SlaterElmDiff_fcmp
    }else{
      info = CalculateMAll_fcmp(eleIdx,qpStart,qpEnd); // InvM,PfM will change
    }
    StopTimer(40);

    if(info!=0) {
      fprintf(stderr,"warning: VMCMainCal rank:%d sample:%d info:%d (CalculateMAll)\n",rank,sample,info);
      continue;
    }
#ifdef _DEBUG_VMCCAL
    printf("  Debug: sample=%d: CalculateIP \n",sample);
#endif
    if(AllComplexFlag==0){
      ip = CalculateIP_real(PfM_real,qpStart,qpEnd,MPI_COMM_SELF);
    }else{
      ip = CalculateIP_fcmp(PfM,qpStart,qpEnd,MPI_COMM_SELF);
    }

#ifdef _DEBUG_VMCCAL
    printf("  Debug: sample=%d: LogProjVal \n",sample);
#endif
    x = LogProjVal(eleProjCnt);
    /* calculate reweight */
    if (reweight==1){
       w = 2.0*(log(cabs(ip))+x);
       if (FlagRBM) {
         w += 2.0*creal(LogWeightRBM(rbmCnt));
       }
       w = exp(w - logSqPfFullSlater[sample]);
    }else{
       w =1.0;
    }

#ifdef _DEBUG_VMCCAL
    printf("  Debug: sample=%d: isfinite \n",sample);
#endif
    if( !isfinite(w) ) {
      fprintf(stderr,"warning: VMCMainCal rank:%d sample:%d w=%e\n",rank,sample,w);
      continue;
    }

    StartTimer(41);
    /* calculate energy */
#ifdef _DEBUG_VMCCAL
    printf("  Debug: sample=%d: calculateHam \n",sample);
#endif
    if(AllComplexFlag==0){
#ifdef _DEBUG_VMCCAL
      printf("  Debug: sample=%d: calculateHam_real \n",sample);
#endif
      e = CalculateHamiltonian_real(creal(ip),eleIdx,eleCfg,eleNum,eleProjCnt);
    }else{
#ifdef _DEBUG_VMCCAL
      printf("  Debug: sample=%d: calculateHam_cmp \n",sample);
#endif
      e = CalculateHamiltonian(ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
    }
    //printf("MDEBUG: %lf %lf \n",creal(e),cimag(e));
    StopTimer(41);

#ifdef _DEBUG_VMCCAL
    printf("  Debug: sample=%d: e = %lf %lf \n",sample, creal(e), cimag(e));
#endif
    if( !isfinite(creal(e) + cimag(e)) ) {
      fprintf(stderr,"warning: VMCMainCal rank:%d sample:%d e=%e\n",rank,sample,creal(e)); //TBC
      continue;
    }

    Wc += w;
    Etot  += w * e;
    Etot2 += w * conj(e) * e;
#ifdef _DEBUG_VMCCAL
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

      if (FlagRBM) {
        RBMDiff(SROptO+2*NProj+2,rbmCnt,eleNum);
      }

      StartTimer(42);
      /* SlaterElmDiff */
      SlaterElmDiff_fcmp(SROptO+2*NProj+2*NRBM+2,ip,eleIdx); //TBC: using InvM not InvM_real
      StopTimer(42);

      if(FlagOptTrans>0) { // this part will be not used
        calculateOptTransDiff(SROptO+2*NProj+2*FlagRBM*NRBM+2*NSlater+2, ip); //TBC
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
#ifdef _DEBUG_VMCCAL
      fprintf(stdout, "Debug: Start: CalcGreenFunc\n");
#endif
      CalculateGreenFunc(w,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
      StopTimer(42);

      if(NLanczosMode>0 && NLanczosEstimatorMode==1){
#ifdef _DEBUG_VMCCAL
  fprintf(stdout, "Debug: Start: Lanczos\n");
#endif
        if(NLanczosStep==1){
        // ignoring Lanczos: to be added
        /* Calculate local QQQQ */
        StartTimer(43);
        if(AllComplexFlag==0) {
          LSLocalQ_real(creal(e),creal(ip),eleIdx,eleCfg,eleNum,eleProjCnt, LSLQ_real);
          RecordPowerLanczosSupportLSLQSample(
              LSLQ_real,NULL,w,parentRank,sample,comm_parent);
          calculateQQQQ_real(QQQQ_real,LSLQ_real,w,NLSHam);
        }else{
          LSLocalQ(e,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt, LSLQ);
          RecordPowerLanczosSupportLSLQSample(
              NULL,LSLQ,w,parentRank,sample,comm_parent);
          calculateQQQQ(QQQQ,LSLQ,w,NLSHam);
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
        StopTimer(43);

        // LanczosGreen
        if(NLanczosMode>1){
          // Calculate local QcisAjsQ
          StartTimer(44);
          if(AllComplexFlag==0) {

            LSLocalCisAjs_real(creal(e),creal(ip),eleIdx,eleCfg,eleNum,eleProjCnt);
            calculateQCAQ_real(QCisAjsQ_real,LSLCisAjs_real,LSLQ_real,w,NLSHam,NCisAjs);
            calculateQCACAQ_real(QCisAjsCktAltQ_real,LSLCisAjs_real,w,NLSHam,NCisAjs,
                                 NCisAjsCktAlt, CisAjsCktAltIdx);
            calculateQCACAQDC_real(QCisAjsCktAltQDC_real,LSLQ_real,w,NLSHam,NCisAjs,
                                   NCisAjsCktAltDC,eleIdx,eleCfg,eleNum,eleProjCnt,creal(e),creal(ip));

          }
          else{
            LSLocalCisAjs(e,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
            calculateQCAQ(QCisAjsQ,LSLCisAjs,LSLQ,w,NLSHam,NCisAjs);
            calculateQCACAQ(QCisAjsCktAltQ,LSLCisAjs,w,NLSHam,NCisAjs,
                            NCisAjsCktAlt,CisAjsCktAltIdx);
            calculateQCACAQDC(QCisAjsCktAltQDC,LSLQ,w,NLSHam,NCisAjs,
                              NCisAjsCktAltDC,eleIdx,eleCfg,eleNum,eleProjCnt,e,ip,rbmCnt);
          }
          StopTimer(44);
        }
        }else{
          Lanczos2SolveStatus lanczos2Status;
          StartTimer(95);
          Lanczos2StateSnapshotCapture(
              &lanczos2StateSnapshot,eleIdx,eleCfg,eleNum,eleProjCnt);
          if(AllComplexFlag==0) {
            if(NExUpdatePath == 2) {
              lanczos2Status =
                  CalculateLS2HeisenbergLocalPowerReal(
                      creal(e),creal(ip),eleIdx,eleCfg,eleNum,eleProjCnt,
                      &lanczos2HeisenbergScratch,LS2LocalPower_real) == 0
                      ? LANCZOS2_SOLVE_OK
                      : LANCZOS2_SOLVE_INVALID_ARGUMENT;
            }else{
              CalculateLS2LocalPower_real(
                  creal(e),creal(ip),eleIdx,eleCfg,eleNum,eleProjCnt,
                  LS2LocalPower_real);
              lanczos2Status = LANCZOS2_SOLVE_OK;
            }
#ifdef MVMC_ENABLE_FAULT_INJECTION
            if(lanczos2InjectNonfinite && parentRank == 0 &&
               sample == sampleStart) {
              LS2LocalPower_real[LANCZOS2_POWER_COUNT-1] = NAN;
            }
#endif
            for(i=0;
                lanczos2Status == LANCZOS2_SOLVE_OK &&
                i<LANCZOS2_POWER_COUNT;
                i++) {
              if(!isfinite(LS2LocalPower_real[i])) {
                lanczos2Status = LANCZOS2_SOLVE_NONFINITE_MOMENT;
                break;
              }
            }
            if(lanczos2Status == LANCZOS2_SOLVE_OK) {
              memcpy(LS2SamplePower_real +
                         (size_t)sample * LANCZOS2_POWER_COUNT,
                     LS2LocalPower_real,
                     LANCZOS2_POWER_COUNT * sizeof(double));
            }
          }else{
            if(NExUpdatePath == 2) {
              lanczos2Status =
                  CalculateLS2HeisenbergLocalPowerComplex(
                      e,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt,
                      &lanczos2HeisenbergScratch,LS2LocalPower) == 0
                      ? LANCZOS2_SOLVE_OK
                      : LANCZOS2_SOLVE_INVALID_ARGUMENT;
            }else{
              CalculateLS2LocalPower(
                  e,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt,
                  LS2LocalPower);
              lanczos2Status = LANCZOS2_SOLVE_OK;
            }
#ifdef MVMC_ENABLE_FAULT_INJECTION
            if(lanczos2InjectNonfinite && parentRank == 0 &&
               sample == sampleStart) {
              LS2LocalPower[LANCZOS2_POWER_COUNT-1] = NAN + 0.0*I;
            }
#endif
            for(i=0;
                lanczos2Status == LANCZOS2_SOLVE_OK &&
                i<LANCZOS2_POWER_COUNT;
                i++) {
              if(!isfinite(creal(LS2LocalPower[i])) ||
                 !isfinite(cimag(LS2LocalPower[i]))) {
                lanczos2Status = LANCZOS2_SOLVE_NONFINITE_MOMENT;
                break;
              }
            }
            if(lanczos2Status == LANCZOS2_SOLVE_OK) {
              memcpy(LS2SamplePower +
                         (size_t)sample * LANCZOS2_POWER_COUNT,
                     LS2LocalPower,
                     LANCZOS2_POWER_COUNT * sizeof(double complex));
            }
          }
          if(!Lanczos2StateSnapshotMatches(
                 &lanczos2StateSnapshot,eleIdx,eleCfg,eleNum,eleProjCnt)) {
            fprintf(stderr,
                    "Error: Lanczos2 state restoration failure on rank %d, "
                    "sample %d.\n",
                    parentRank,sample);
            MPI_Abort(comm_parent,EXIT_FAILURE);
          }
          if(lanczos2Status != LANCZOS2_SOLVE_OK) {
            fprintf(stderr,
                    "Error: Lanczos2 local power failure on rank %d, "
                    "sample %d, stage accumulation: %s.\n",
                    parentRank,sample,Lanczos2SolveError(lanczos2Status));
            MPI_Abort(comm_parent,EXIT_FAILURE);
          }
          LS2SampleWeight[sample] = w;
          LS2SampleValid[sample] = 1;
          if(AllComplexFlag==0) {
            lanczos2Status = RecordPowerLanczosSupportSampleReal(
                PowerLanczosSupportSampleData +
                    (size_t)sample *
                    POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH,
                LS2LocalPower_real,LANCZOS2_POWER_COUNT,w);
          }else{
            lanczos2Status = RecordPowerLanczosSupportSampleComplex(
                PowerLanczosSupportSampleData +
                    (size_t)sample *
                    POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH,
                LS2LocalPower,LANCZOS2_POWER_COUNT,w);
          }
          if(lanczos2Status != LANCZOS2_SOLVE_OK) {
            fprintf(stderr,
                    "Error: power-Lanczos support sample failure on rank %d, "
                    "sample %d: %s.\n",
                    parentRank,sample,Lanczos2SolveError(lanczos2Status));
            MPI_Abort(comm_parent,EXIT_FAILURE);
          }
          if(lanczosOracleDump != NULL) {
            fprintf(lanczosOracleDump, "sample %d occ", sample);
            for(i=0; i<Nsite2; i++) {
              fprintf(lanczosOracleDump, " %d", eleNum[i]);
            }
            fprintf(lanczosOracleDump, " ls2power");
            if(AllComplexFlag == 0) {
              for(i=0; i<LANCZOS2_POWER_COUNT; i++) {
                fprintf(lanczosOracleDump, " %.17e %.17e",
                        LS2LocalPower_real[i],0.0);
              }
            }else{
              for(i=0; i<LANCZOS2_POWER_COUNT; i++) {
                fprintf(lanczosOracleDump, " %.17e %.17e",
                        creal(LS2LocalPower[i]),cimag(LS2LocalPower[i]));
              }
            }
            fprintf(lanczosOracleDump, "\n");
          }
          StopTimer(95);
        }
      }
    }
  } /* end of for(sample) */

#ifdef MVMC_ENABLE_FAULT_INJECTION
  if(lanczos2CalHCAGuardAudit) {
    const long long directCount =
        AllComplexFlag == 0
            ? Lanczos2CalHCAGuardAuditRealDirectCount
            : Lanczos2CalHCAGuardAuditComplexDirectCount;
    const long long zeroComponentCount =
        AllComplexFlag == 0
            ? Lanczos2CalHCAGuardAuditRealZeroComponentCount
            : Lanczos2CalHCAGuardAuditComplexZeroComponentCount;
    if(directCount <= 0 || zeroComponentCount <= 0) {
      fprintf(stderr,
              "Error: Lanczos2 calHCA guard audit did not exercise a direct "
              "multi-QP zero-component branch "
              "(direct=%lld, zero_component=%lld).\n",
              directCount,zeroComponentCount);
      MPI_Abort(comm_parent, EXIT_FAILURE);
    }
    printf("Lanczos2 calHCA guard audit: %s direct=%lld "
           "zero_component=%lld\n",
           AllComplexFlag == 0 ? "real" : "complex",
           directCount,zeroComponentCount);
  }
#endif

#ifdef MVMC_ENABLE_FAULT_INJECTION
  if(lanczos2HeisenbergScratch.auditEnabled) {
    const long long zeroOverlapCount =
        AllComplexFlag == 0
            ? Lanczos2HeisenbergZeroOverlapAuditRealCount
            : Lanczos2HeisenbergZeroOverlapAuditComplexCount;
    if(lanczos2HeisenbergScratch.auditedSamples <= 0 ||
       zeroOverlapCount <= 0) {
      fprintf(stderr,
              "Error: Lanczos2 Heisenberg audit did not exercise all "
              "required paths (samples=%lld, zero_overlap=%lld).\n",
              lanczos2HeisenbergScratch.auditedSamples,zeroOverlapCount);
      MPI_Abort(comm_parent, EXIT_FAILURE);
    }
    printf("Lanczos2 Heisenberg audit: %s samples=%lld "
           "zero_overlap=%lld green2=%lld green4=%lld green6=%lld\n",
           AllComplexFlag == 0 ? "real" : "complex",
           lanczos2HeisenbergScratch.auditedSamples,zeroOverlapCount,
           lanczos2HeisenbergScratch.greenCalls[0],
           lanczos2HeisenbergScratch.greenCalls[1],
           lanczos2HeisenbergScratch.greenCalls[2]);
  }
#endif

  if(lanczosOracleDump != NULL) fclose(lanczosOracleDump);
  Lanczos2StateSnapshotFree(&lanczos2StateSnapshot);
  Lanczos2HeisenbergScratchFree(&lanczos2HeisenbergScratch);

// calculate OO and HO at NVMCCalMode==0
  if(NVMCCalMode==0){
    if(NSRCG!=0 || NStoreO!=0){
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

void clearPhysQuantity(){
  int i,n;
  double complex *vec;
  double  *vec_real;
//[s] MERGE BY TM
  Wc = Etot = Etot2 = Sztot=Sztot2 =0.0;//fsz
  //Wc = Etot = Etot2 = 0.0;
  Dbtot = Dbtot2 = 0.0;
//[e] MERGE BY TM
  if(NVMCCalMode==0) {
    /* SROptOO, SROptHO, SROptO */
    if(NSRCG!=0){
      n = (2*SROptSize)*4; // TBC
    }else{
      n = (2*SROptSize)*(2*SROptSize+2); // TBC
    }
    vec = SROptOO;
    #pragma omp parallel for default(shared) private(i)
    for(i=0;i<n;i++) vec[i] = 0.0+0.0*I;

    if (AllComplexFlag == 0){
    // only for real variables
      if(NSRCG!=0){
        n = (SROptSize)*4; // TBC
      }else{
        n = (SROptSize)*(SROptSize+2); // TBC
      }
      vec_real = SROptOO_real;
      #pragma omp parallel for default(shared) private(i)
      for(i=0;i<n;i++) vec_real[i] = 0.0;
    }
  } else if(NVMCCalMode==1) {
    /* CisAjs, CisAjsCktAlt, CisAjsCktAltDC */
    //n = NCisAjs+NCisAjsCktAlt+NCisAjsCktAltDC;
    n = NCisAjs+NCisAjsCktAlt+NCisAjsCktAltDC + NTwist;
    vec = PhysCisAjs;
#pragma omp parallel for default(shared) private(i)
    for(i=0;i<n;i++) vec[i] = 0.0+0.0*I;

    vec = PhysNBodyG;
#pragma omp parallel for default(shared) private(i)
    for(i=0;i<NNBodyG;i++) vec[i] = 0.0+0.0*I;

    if(NLanczosMode>0 && NLanczosEstimatorMode==1 && NLanczosStep==1) {
      /* QQQQ, LSLQ */
        //[TODO]: Check the value n
      n = NLSHam*NLSHam*NLSHam*NLSHam + NLSHam*NLSHam;
      vec = QQQQ;
#pragma omp parallel for default(shared) private(i)
      for(i=0;i<n;i++) vec[i] = 0.0+0.0*I;

      n = NLSHam*NLSHam*NLSHam*NLSHam + NLSHam*NLSHam;
      vec_real = QQQQ_real;
#pragma omp parallel for default(shared) private(i)
      for(i=0;i<n;i++) vec_real[i] = 0.0;

      if(NLanczosMode>1) {
        /* QCisAjsQ, QCisAjsCktAltQ, LSLCisAjs */
        //[TODO]: Check the value n
        n = NLSHam*NLSHam*NCisAjs + NLSHam*NLSHam*(NCisAjsCktAltDC+NCisAjsCktAlt)
          + NLSHam*NCisAjs;
        vec = QCisAjsQ;
#pragma omp parallel for default(shared) private(i)
        for(i=0;i<n;i++) vec[i] = 0.0+0.0*I;

        n = NLSHam*NLSHam*NCisAjs + NLSHam*NLSHam*(NCisAjsCktAltDC+NCisAjsCktAlt)
          + NLSHam*NCisAjs;
        vec_real = QCisAjsQ_real;
#pragma omp parallel for default(shared) private(i)
        for(i=0;i<n;i++) vec_real[i] = 0.0;
      }
    }
    if(NLanczosMode>0 && NLanczosEstimatorMode==1 && NLanczosStep==2) {
      if(AllComplexFlag==0) {
        for(i=0;i<LANCZOS2_MOMENT_COUNT+LANCZOS2_POWER_COUNT;i++) {
          LS2Moment_real[i] = 0.0;
        }
        for(i=0;i<NVMCSample*LANCZOS2_POWER_COUNT;i++) {
          LS2SamplePower_real[i] = 0.0;
        }
      }else{
        for(i=0;i<LANCZOS2_MOMENT_COUNT+LANCZOS2_POWER_COUNT;i++) {
          LS2Moment[i] = 0.0;
        }
        for(i=0;i<NVMCSample*LANCZOS2_POWER_COUNT;i++) {
          LS2SamplePower[i] = 0.0;
        }
      }
      for(i=0;i<NVMCSample;i++) {
        LS2SampleWeight[i] = 0.0;
        LS2SampleValid[i] = 0;
      }
      LS2MomentBasisShift = 0.0;
    }
    if(NLanczosMode>0 && NLanczosEstimatorMode==1) {
      for(i=0;
          i<PowerLanczosSupportSampleCapacity *
                POWER_LANCZOS_SUPPORT_SAMPLE_WIDTH;
          i++) {
        PowerLanczosSupportSampleData[i] = 0.0;
      }
    }
  }
  return;
}

void calculateOptTransDiff(double complex *srOptO, const double complex ipAll) {
  int i,j;
  double complex ip;
  double complex *pfM;

  for(i=0;i<NQPOptTrans;++i) {
    ip = 0.0;
    pfM = PfM + i*NQPFix;
    for(j=0;j<NQPFix;++j) {
      ip += QPFixWeight[j] * pfM[j];
    }
    srOptO[i] = ip/ipAll;
  }

  return;
}

void calculateOO_Store_real(double *srOptOO_real, double *srOptHO_real, double *srOptO_Store_real,
                 const double w, const double e, int srOptSize, int sampleSize) {

  int i,j;
  char jobz, uplo;
  double alpha,beta,o;

  alpha = 1.0;
  beta  = 0.0;

  jobz = 'N';
  uplo = 'T';
  if(NSRCG==0){
    M_DGEMM(&jobz,&uplo,&srOptSize,&srOptSize,&sampleSize,&alpha,srOptO_Store_real,&srOptSize,srOptO_Store_real,&srOptSize,&beta,srOptOO_real,&srOptSize);
  }else{
#pragma omp parallel for default(shared) private(i)
#pragma loop noalias
    for(i=0; i<srOptSize; ++i){
      srOptOO_real[i] = 0.0;
      srOptOO_real[i+srOptSize] = 0.0;
    }
    if(reweight==1){
    for(j=0; j<sampleSize; ++j){
      /* Store holds sqrt(w_j)*O; row 0 (O_0 = 1) gives sqrt(w_j),
         so multiplying once more by it makes <O> w-weighted,
         consistent with the dposv/GEMM paths. With reweight off
         (w_j = 1) the original loop below is kept bit-identical. */
      const double sqrtw = srOptO_Store_real[j*srOptSize];
#pragma omp parallel for default(shared) private(i,o)
#pragma loop noalias
    for(i=0; i<srOptSize; ++i){
      o = srOptO_Store_real[i+j*srOptSize];
      srOptOO_real[i] += sqrtw*o;
      srOptOO_real[i+srOptSize] += o*o;
    }}
    }else{
    for(j=0; j<sampleSize; ++j){
#pragma omp parallel for default(shared) private(i,o)
#pragma loop noalias
    for(i=0; i<srOptSize; ++i){
      o = srOptO_Store_real[i+j*srOptSize];
      srOptOO_real[i] += o;
      srOptOO_real[i+srOptSize] += o*o;
    }}
    }
  }

  return;
}



void calculateOO_Store(double complex *srOptOO, double complex *srOptHO, double complex *srOptO_Store,
                 const double w, const double complex e, int srOptSize, int sampleSize) {

  int i,j;
  char jobz, uplo;
  double complex alpha,beta,o;

  alpha = 1.0;
  beta  = 0.0;

  jobz = 'N';
  uplo = 'C';
  if(NSRCG==0){
    M_ZGEMM(&jobz,&uplo,&srOptSize,&srOptSize,&sampleSize,&alpha,srOptO_Store,&srOptSize,srOptO_Store,&srOptSize,&beta,srOptOO,&srOptSize);
  }else{
#pragma omp parallel for default(shared) private(i)
#pragma loop noalias
    for(i=0; i<srOptSize; ++i){
      srOptOO[i] = 0.0;
      srOptOO[i+srOptSize] = 0.0;
    }
    if(reweight==1){
    for(j=0; j<sampleSize; ++j){
      /* Store holds sqrt(w_j)*O; row 0 (O_0 = 1) gives sqrt(w_j),
         so multiplying once more by it makes <O> w-weighted,
         consistent with the dposv/GEMM paths. With reweight off
         (w_j = 1) the original loop below is kept bit-identical. */
      const double sqrtw = creal(srOptO_Store[j*srOptSize]);
#pragma omp parallel for default(shared) private(i,o)
#pragma loop noalias
    for(i=0; i<srOptSize; ++i){
      o = srOptO_Store[i+j*srOptSize];
      srOptOO[i] += sqrtw*o;
      srOptOO[i+srOptSize] += creal(o)*creal(o)+cimag(o)*cimag(o);
    } }
    }else{
    for(j=0; j<sampleSize; ++j){
#pragma omp parallel for default(shared) private(i,o)
#pragma loop noalias
    for(i=0; i<srOptSize; ++i){
      o = srOptO_Store[i+j*srOptSize];
      srOptOO[i] += o;
      srOptOO[i+srOptSize] += creal(o)*creal(o)+cimag(o)*cimag(o);
    } }
    }
  }

  return;
}




//void calculateOO(double complex *srOptOO, double complex *srOptHO, const double complex *srOptO,
//                 const double w, const double complex e, const int srOptSize) {
//  double we=w*e;
//
//  #define M_DAXPY daxpy_
//  #define M_DGER dger_
//
//  extern int M_DAXPY(const int *n, const double *alpha, const double *x, const int *incx,
//                     double *y, const int *incy);
//  extern int M_DGER(const int *m, const int *n, const double *alpha,
//                    const double *x, const int *incx, const double *y, const int *incy,
//                    double *a, const int *lda);
//  int m,n,incx,incy,lda;
//  m=n=lda=srOptSize;
//  incx=incy=1;
//
//  /* OO[i][j] += w*O[i]*O[j] */
//  M_DGER(&m, &n, &w, srOptO, &incx, srOptO, &incy, srOptOO, &lda);
//
//  /* HO[i] += w*e*O[i] */
//  M_DAXPY(&n, &we, srOptO, &incx, srOptHO, &incy);
//
//  return;
//}

void calculateOO_matvec(double complex *srOptOO, double complex *srOptHO, const double complex *srOptO,
                 const double complex w, const double complex e, const int srOptSize) {
  double complex we=w*e;

  int m,n,incx,incy,lda;
  m=n=lda=2*srOptSize;
  incx=incy=1;

//   OO[i][j] += w*O[i]*O[j]
  M_ZGERC(&m, &n, &w, srOptO, &incx, srOptO, &incy, srOptOO, &lda);
//   HO[i] += w*e*O[i]
  M_ZAXPY(&n, &we, srOptO, &incx, srOptHO, &incy);
  return;
}

void calculateOO(double complex *srOptOO, double complex *srOptHO, const double complex *srOptO,
                 const double w, const double complex e, const int srOptSize){
  int i,j;
  double complex tmp;
  #pragma omp parallel for default(shared) private(j,tmp)
  //    private(i,j,tmp,srOptOO)
#pragma loop noalias
  for(j=0;j<2*srOptSize;j++) {
    tmp                            = w * srOptO[j];
    srOptOO[0*(2*srOptSize)+j]    += tmp;      // update O
    srOptOO[1*(2*srOptSize)+j]    += 0.0;      // update
    srOptHO[j]                    += e * tmp;  // update HO
  }

  #pragma omp parallel for default(shared) private(i,j,tmp)
#pragma loop noalias
  for(i=2;i<2*srOptSize;i++) {
    tmp            = w * srOptO[i];
    for(j=0;j<2*srOptSize;j++) {
      srOptOO[i*(2*srOptSize)+j] += w*(srOptO[j])*conj(srOptO[i]); // TBC
      //srOptOO[j+i*(2*srOptSize)] += w*(srOptO[j])*(srOptO[i]); // TBC
    }
  }

  return;
}

void calculateOO_real(double *srOptOO, double *srOptHO, const double *srOptO,
                 const double w, const double e, const int srOptSize) {
  double we=w*e;
  int m,n,incx,incy,lda;
  m=n=lda=srOptSize;
  incx=incy=1;

  /* OO[i][j] += w*O[i]*O[j] */
  M_DGER(&m, &n, &w, srOptO, &incx, srOptO, &incy, srOptOO, &lda);

  /* HO[i] += w*e*O[i] */
  M_DAXPY(&n, &we, srOptO, &incx, srOptHO, &incy);

  return;
}

void calculateQQQQ_real(double *qqqq, const double *lslq, const double w, const int nLSHam) {
  const int n=nLSHam*nLSHam*nLSHam*nLSHam;
  int rq,rp,ri,rj;
  int i,tmp;

  /* QQQQ[rq][rp][ri][rj] += w * LSLQ[rq][ri] * LSLQ[rp][rj] */
# pragma omp parallel for default(shared) private(i,tmp,rq,rp,ri,rj)
  for(i=0;i<n;++i) {
    rj = i%nLSHam;   tmp=i/nLSHam;
    ri = tmp%nLSHam; tmp=tmp/nLSHam;
    rp = tmp%nLSHam; tmp=tmp/nLSHam;
    rq = tmp%nLSHam;

    qqqq[i] += w * lslq[rq*nLSHam+ri] * lslq[rp*nLSHam+rj];
  }

  return;
}


void calculateQQQQ(double complex *qqqq, const double complex*lslq, const double w, const int nLSHam) {
  const int n=nLSHam*nLSHam*nLSHam*nLSHam;
  int rq,rp,ri,rj;
  int i,tmp;

  /* QQQQ[rq][rp][ri][rj] += w * LSLQ[rq][ri] * LSLQ[rp][rj] */
  # pragma omp parallel for default(shared) private(i,tmp,rq,rp,ri,rj)
  for(i=0;i<n;++i) {
    rj = i%nLSHam;   tmp=i/nLSHam;
    ri = tmp%nLSHam; tmp=tmp/nLSHam;
    rp = tmp%nLSHam; tmp=tmp/nLSHam;
    rq = tmp%nLSHam;

    qqqq[i] += w * conj(lslq[rq*nLSHam+ri]) * lslq[rp*nLSHam+rj];
   //fprintf(stdout, "Debug: qqqq[%d]= %lf, %lf \n", i, creal(qqqq[i]),cimag(qqqq[i]));
  }

  return;
}

void calculateQCAQ(double complex*qcaq, const double complex*lslca, const double complex*lslq,
                   const double w, const int nLSHam, const int nCA) {
  const int n=nLSHam*nLSHam*nCA;
  int rq,rp,idx;
  int i,tmp;

  /* QCisAjsQ[rq][rp][idx] += w * LSLCisAjs[rq][idx] * LSLQ[rp][0] */
# pragma omp parallel for default(shared) private(i,tmp,idx,rp,rq)
  for(i=0;i<n;++i) {
    idx = i%nCA;     tmp = i/nCA;
    rp = tmp%nLSHam; tmp = tmp/nLSHam;
    rq = tmp%nLSHam;

    qcaq[i] += w * conj(lslca[rq*nCA+idx]) * lslq[rp*nLSHam];
  }

  return;
}

void calculateQCACAQ(double complex *qcacaq, const double complex *lslca, const double w,
                     const int nLSHam, const int nCA, const int nCACA,
                     int **cacaIdx) {
  const int n=nLSHam*nLSHam*nCACA;
  int rq,rp,ri,rj,idx;
  int i,tmp;

  /* QCisAjsCktAltQ[rq][rp][idx] += w * LSLCisAjs[rq][ri] * LSLCisAjs[rp][rj] */
# pragma omp parallel for default(shared) private(i,tmp,idx,rp,rq,ri,rj)
  for(i=0;i<n;++i) {
    idx = i%nCACA;   tmp = i/nCACA;
    rp = tmp%nLSHam; tmp = tmp/nLSHam;
    rq = tmp%nLSHam;

    ri = cacaIdx[idx][0];
    rj = cacaIdx[idx][1];
    //fprintf(stdout, "Debug: ri= %d, rj=%d, rp=%d, rq=%d \n", ri, rj, rp, rq);

    qcacaq[i] += w * conj(lslca[rq*nCA+ri]) * lslca[rp*nCA+rj];
  }

  return;
}

void calculateQCAQ_real(double *qcaq, const double *lslca, const double *lslq,
                   const double w, const int nLSHam, const int nCA) {
  const int n=nLSHam*nLSHam*nCA;
  int rq,rp,idx;
  int i,tmp;

  /* QCisAjsQ[rq][rp][idx] += w * LSLCisAjs[rq][idx] * LSLQ[rp][0] */
# pragma omp parallel for default(shared) private(i,tmp,idx,rp,rq)
  for(i=0;i<n;++i) {
    idx = i%nCA;     tmp = i/nCA;
    rp = tmp%nLSHam; tmp = tmp/nLSHam;
    rq = tmp%nLSHam;

    qcaq[i] += w * lslca[rq*nCA+idx] * lslq[rp*nLSHam];
    //fprintf(stdout, "Debug: qcaq[%d]= %lf, %lf \n", i, creal(qcaq[i]),cimag(qcaq[i]));

  }

  return;
}

void calculateQCACAQ_real(double *qcacaq, const double *lslca, const double w,
                     const int nLSHam, const int nCA, const int nCACA,
                     int **cacaIdx) {
  const int n=nLSHam*nLSHam*nCACA;
  int rq,rp,ri,rj,idx;
  int i,tmp;


  /* QCisAjsCktAltQ[rq][rp][idx] += w * LSLCisAjs[rq][ri] * LSLCisAjs[rp][rj] */
# pragma omp parallel for default(shared) private(i,tmp,idx,rp,rq,ri,rj)
  for(i=0;i<n;++i) {
    idx = i%nCACA;   tmp = i/nCACA;
    rp = tmp%nLSHam; tmp = tmp/nLSHam;
    rq = tmp%nLSHam;

    ri = cacaIdx[idx][0];
    rj = cacaIdx[idx][1];
    qcacaq[i] += w * lslca[rq*nCA+ri] * lslca[rp*nCA+rj];

  }
  /*
  for(i=0;i<n;++i) {
    fprintf(stdout, "Debug: qcacaq[%d]= %lf, %lf \n", i, creal(qcacaq[i]),cimag(qcacaq[i]));
    //fprintf(stdout, "Debug: lslca[%d]= %lf, %lf \n", i, creal(lslca[i]),cimag(lslca[i]));
  }
  */

  return;

}

void calculateQCACAQDC_real(double *qcacaq, const double *lslq, const double w,
                            const int nLSHam, const int nCA, const int nCACA,
                            int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,
                            const double h1, const double ip) {
  const int n=nLSHam*nLSHam*nCACA;
  int rq,rp,ri,rj,rk,rl,s,t,idx;
  int i,tmp;

  for(i=0;i<n;++i) {
    idx = i%nCACA;   tmp = i/nCACA;
    rp = tmp%nLSHam; tmp = tmp/nLSHam;
    rq = tmp%nLSHam;

    if (rq) {
      // Calls 2-body routines to evaluate HCACA.
      ri = CisAjsCktAltDCIdx[idx][0];
      rj = CisAjsCktAltDCIdx[idx][2];
      s  = CisAjsCktAltDCIdx[idx][1];
      rk = CisAjsCktAltDCIdx[idx][4];
      rl = CisAjsCktAltDCIdx[idx][6];
      t  = CisAjsCktAltDCIdx[idx][5];
      qcacaq[i] += w * lslq[rp*nLSHam] * calHCACA_real(ri,rj,rk,rl,s,t,h1,ip,
                                                       eleIdx,eleCfg,eleNum,eleProjCnt);
    } else
      qcacaq[i] += w * lslq[rp*nLSHam] * LocalCisAjsCktAltDC[idx];
  }

  return;
}

void calculateQCACAQDC(double complex *qcacaq, const double complex *lslq, const double w,
                       const int nLSHam, const int nCA, const int nCACA,
                       int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt,
                       const double complex h1, const double complex ip,double complex *rbmCnt) {
  const int n=nLSHam*nLSHam*nCACA;
  int rq,rp,ri,rj,rk,rl,s,t,idx;
  int i,tmp;

  for(i=0;i<n;++i) {
    idx = i%nCACA;   tmp = i/nCACA;
    rp = tmp%nLSHam; tmp = tmp/nLSHam;
    rq = tmp%nLSHam;

    if (rq) {
      // Calls 2-body routines to evaluate HCACA.
      ri = CisAjsCktAltDCIdx[idx][0];
      rj = CisAjsCktAltDCIdx[idx][2];
      s  = CisAjsCktAltDCIdx[idx][1];
      rk = CisAjsCktAltDCIdx[idx][4];
      rl = CisAjsCktAltDCIdx[idx][6];
      t  = CisAjsCktAltDCIdx[idx][5];
      qcacaq[i] += w * lslq[rp*nLSHam] * calHCACA(ri,rj,rk,rl,s,t,h1,ip,
                                                  eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
    } else
      qcacaq[i] += w * lslq[rp*nLSHam] * LocalCisAjsCktAltDC[idx];
  }

  return;
}
#endif
