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
 * calculate physical quantities
 *-------------------------------------------------------------
 * by Satoshi Morita 
 *-------------------------------------------------------------*/
#ifndef _SRC_VMCCAL
#define _SRC_VMCCAL
#include "vmccal.h"
#include "matrix.h"
#include "calham_real.h"
#include "calham.h"
#include "lslocgrn_real.c"
#include "lslocgrn.c"
#include "calgrn.c"

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

void VMCMainCal(MPI_Comm comm) {
  int *eleIdx,*eleCfg,*eleNum,*eleProjCnt;
  double complex e,ip;
  double w;
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

  int rank,size,int_i;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
#ifdef _DEBUG
  printf("  Debug: SplitLoop\n");
#endif
  SplitLoop(&sampleStart,&sampleEnd,NVMCSample,rank,size);

  /* initialization */
  StartTimer(24);
  clearPhysQuantity();
  StopTimer(24);
  for(sample=sampleStart;sample<sampleEnd;sample++) {
    eleIdx = EleIdx + sample*Nsize;
    eleCfg = EleCfg + sample*Nsite2;
    eleNum = EleNum + sample*Nsite2;
    eleProjCnt = EleProjCnt + sample*NProj;
//DEBUG
    /* for(i=0;i<Nsite;i++) {
      printf("Debug: sample=%d: i=%d  up=%d down =%d \n",sample,i,eleCfg[i+0*Nsite],eleCfg[i+1*Nsite]);
      }*/
//DEBUG

    StartTimer(40);
#ifdef _DEBUG
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
#ifdef _DEBUG
    printf("  Debug: sample=%d: CalculateIP \n",sample);
#endif
    if(AllComplexFlag==0){
      ip = CalculateIP_real(PfM_real,qpStart,qpEnd,MPI_COMM_SELF);
    }else{
      ip = CalculateIP_fcmp(PfM,qpStart,qpEnd,MPI_COMM_SELF);
    } 

    //x = LogProjVal(eleProjCnt);
#ifdef _DEBUG
    printf("  Debug: sample=%d: LogProjVal \n",sample);
#endif
    LogProjVal(eleProjCnt);
    /* calculate reweight */
    //w = exp(2.0*(log(fabs(ip))+x) - logSqPfFullSlater[sample]);
    w =1.0;
#ifdef _DEBUG
    printf("  Debug: sample=%d: isfinite \n",sample);
#endif
    if( !isfinite(w) ) {
      fprintf(stderr,"warning: VMCMainCal rank:%d sample:%d w=%e\n",rank,sample,w);
      continue;
    }

    StartTimer(41);
    /* calculate energy */
#ifdef _DEBUG
    printf("  Debug: sample=%d: calculateHam \n",sample);
#endif
    if(AllComplexFlag==0){
#ifdef _DEBUG
      printf("  Debug: sample=%d: calculateHam_real \n",sample);
#endif
      e = CalculateHamiltonian_real(creal(ip),eleIdx,eleCfg,eleNum,eleProjCnt);
    }else{
#ifdef _DEBUG
      printf("  Debug: sample=%d: calculateHam_cmp \n",sample);
#endif
      e = CalculateHamiltonian(ip,eleIdx,eleCfg,eleNum,eleProjCnt);
    }
    //printf("DEBUG: rank=%d: sample=%d ip= %lf %lf\n",rank,sample,creal(ip),cimag(ip));
    StopTimer(41);
    if( !isfinite(creal(e) + cimag(e)) ) {
      fprintf(stderr,"warning: VMCMainCal rank:%d sample:%d e=%e\n",rank,sample,creal(e)); //TBC
      continue;
    }

    Wc += w;
    Etot  += w * e;
    Etot2 += w * conj(e) * e;
#ifdef _DEBUG
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
      SlaterElmDiff_fcmp(SROptO+2*NProj+2,ip,eleIdx); //TBC: using InvM not InvM_real
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
      if(NStoreO==0){
        //calculateOO_matvec(SROptOO,SROptHO,SROptO,w,e,SROptSize);
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
      CalculateGreenFunc(w,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
      StopTimer(42);

      if(NLanczosMode>0){
        // ignoring Lanczos: to be added
        /* Calculate local QQQQ */
        StartTimer(43);
        if(AllComplexFlag==0) {
          LSLocalQ_real(creal(e),creal(ip),eleIdx,eleCfg,eleNum,eleProjCnt, LSLQ_real);
          calculateQQQQ_real(QQQQ_real,LSLQ_real,w,NLSHam);
        }else{
            LSLocalQ(e,ip,eleIdx,eleCfg,eleNum,eleProjCnt, LSLQ);
            calculateQQQQ(QQQQ,LSLQ,w,NLSHam);
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
                            NCisAjsCktAltDC, CisAjsCktAltLzIdx);
            
          }
          else{

            LSLocalCisAjs(e,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
            calculateQCAQ(QCisAjsQ,LSLCisAjs,LSLQ,w,NLSHam,NCisAjs);
            calculateQCACAQ(QCisAjsCktAltQ,LSLCisAjs,w,NLSHam,NCisAjs,
                            NCisAjsCktAltDC,CisAjsCktAltLzIdx);
             
          }
          StopTimer(44);
        }

      }
    }
  } /* end of for(sample) */

// calculate OO and HO at NVMCCalMode==0
  if(NStoreO!=0 && NVMCCalMode==0){
    sampleSize=sampleEnd-sampleStart;
    if(AllComplexFlag==0){
      StartTimer(45);
      calculateOO_Store_real(SROptOO_real,SROptHO_real,SROptO_Store_real,creal(w),creal(e),SROptSize,sampleSize);
      StopTimer(45);
    }else{
      StartTimer(45);
      calculateOO_Store(SROptOO,SROptHO,SROptO_Store,w,e,2*SROptSize,sampleSize);
      StopTimer(45);
    }
  }
  return;
}

void VMC_BF_MainCal(MPI_Comm comm) {
    int *eleIdx, *eleCfg, *eleNum, *eleProjCnt, *eleProjBFCnt;
    double complex e, ip; //db is double?
    double x, w, db;
    double we, sqrtw;
    int int_i, sampleSize, tmp_i;
    const int qpStart = 0;
    const int qpEnd = NQPFull;
    int sample, sampleStart, sampleEnd;
    int i, info;
    double complex *InvM_Moto, *PfM_Moto;
    double *InvM_real_Moto, *PfM_real_Moto;

    /* optimazation for Kei */
    const int nProj = NProj;
    double complex *srOptO = SROptO;
    double *srOptO_real = SROptO_real;

    double rtmp;

    int rank, size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    SplitLoop(&sampleStart, &sampleEnd, NVMCSample, rank, size);
    //SplitLoop(&sampleStart,&sampleEnd,NExactSample,rank,size);

    /* initialization */
    StartTimer(24);
    clearPhysQuantity();
    StopTimer(24);

    InvM_Moto = InvM;
    PfM_Moto = PfM;

    InvM_real_Moto = InvM_real;
    PfM_real_Moto = PfM_real;

    for (sample = sampleStart; sample < sampleEnd; sample++) {
        eleIdx = EleIdx + sample * Nsize;
        eleCfg = EleCfg + sample * Nsite2;
        eleNum = EleNum + sample * Nsite2;
        eleProjCnt = EleProjCnt + sample * NProj;
        eleProjBFCnt = EleProjBFCnt + sample * 16 * Nsite * Nrange;

        StartTimer(45);
        MakeSlaterElmBF_fcmp(eleNum, eleProjBFCnt);
        StopTimer(45);

        StartTimer(40);
        if (AllComplexFlag == 0) {
#pragma omp parallel for default(shared) private(tmp_i)
            for (tmp_i = 0; tmp_i < NQPFull * (2 * Nsite) * (2 * Nsite); tmp_i++)
                SlaterElm_real[tmp_i] = creal(SlaterElm[tmp_i]);

            info = CalculateMAll_BF_real(eleIdx, qpStart, qpEnd);  // InvM_real,PfM_real will change
#pragma omp parallel for default(shared) private(tmp_i)
            for (tmp_i = 0; tmp_i < NQPFull * (Nsize * Nsize + 1); tmp_i++)
                InvM[tmp_i] = InvM_real[tmp_i]; // InvM will be used in  SlaterElmDiff_fcmp
        } else {//complex
            info = CalculateMAll_BF_fcmp(eleIdx,qpStart,qpEnd);
        }
        StopTimer(40);

        if (info != 0) {
            fprintf(stderr, "waring: VMCMainCal rank:%d sample:%d info:%d (CalculateMAll)\n", rank, sample, info);
            continue;
        }

        if (AllComplexFlag == 0) {
            ip = CalculateIP_real(PfM_real, qpStart, qpEnd, MPI_COMM_SELF);
        } else {
            ip = CalculateIP_fcmp(PfM, qpStart, qpEnd, MPI_COMM_SELF);
        }

        LogProjVal(eleProjCnt);
        /* calculate reweight */
        // TODO: check ~ Is it OK to fix w ?
        //w = exp(2.0*(log(fabs(ip))+x) - logSqPfFullSlater[sample]);
        w = 1.0;
        /*
        if(log(fabs(1.0-w)) > -5) {
            printf("w=%.3e\n",w);
        }
        */
        if (!isfinite(w)) {
            fprintf(stderr, "waring: VMCMainCal rank:%d sample:%d w=%e\n", rank, sample, w);
            continue;
        }

        StartTimer(41);
        if (AllComplexFlag == 0) {
            e = CalculateHamiltonianBF_real(creal(ip), eleIdx, eleCfg, eleNum, eleProjCnt, eleProjBFCnt);
        } else {/* calculate energy */
            e = CalculateHamiltonianBF_fcmp(ip, eleIdx, eleCfg, eleNum, eleProjCnt, eleProjBFCnt);
        }

        /* calculate double occupation D */
        db = CalculateDoubleOccupation(eleIdx, eleCfg, eleNum, eleProjCnt);
        StopTimer(41);
        if (!isfinite(e)) {
            fprintf(stderr, "waring: VMCMainCal rank:%d sample:%d e=%e\n", rank, sample, creal(e));
            continue;
        }

        Wc += w;
        Etot += w * e;
        Etot2 += w * e * e;
        Dbtot += w * db;
        Dbtot2 += w * db * db;

        if (NVMCCalMode == 0) {
            /* Calculate O for correlation fauctors */
            /*SROptO[0] = 1.0;
#pragma loop noalias
            for(i=0;i<nProj;i++) srOptO[i+1] = (double)(eleProjCnt[i]);
*/
            srOptO[0] = 1.0 + 0.0 * I;//   real
            srOptO[1] = 0.0 + 0.0 * I;//   real
#pragma loop noalias
            for (i = 0; i < nProj; i++) {
                srOptO[(i + 1) * 2] = (double) (eleProjCnt[i]); // even real
                srOptO[(i + 1) * 2 + 1] = 0.0 + 0.0 * I;               // odd  comp
            }

            //StartTimer(44);
            /* BackflowDiff */
            //BackFlowDiff(SROptO+NProj+1,ip,eleIdx,eleNum,eleProjCnt,eleProjBFCnt);
            BackFlowDiff_fcmp(SROptO + 2 * NProj + 2, ip, eleIdx, eleNum, eleProjCnt, eleProjBFCnt);
            //StopTimer(44);

            StartTimer(42);
            /* SlaterElmDiff */
            //SlaterElmBFDiff_fcmp(SROptO+NProj+NProjBF+1,ip,eleIdx,eleNum,eleCfg,eleProjCnt,eleProjBFCnt);
            SlaterElmBFDiff_fcmp(SROptO + 2 * NProj + 2 * NProjBF + 2, ip, eleIdx, eleNum, eleCfg, eleProjCnt,
                                 eleProjBFCnt);
            StopTimer(42);

            if (FlagOptTrans > 0) {
                calculateOptTransDiff(SROptO + 2 * NProj + 2 * NProjBF + 2 * NSlater + 2, ip);
            }

            //[s] this part will be used for real varaibles
            if (AllComplexFlag == 0) {
#pragma loop noalias
                for (i = 0; i < SROptSize; i++) {
                    srOptO_real[i] = creal(srOptO[2 * i]);
                }
            }
            //[e]

            StartTimer(43);
            /* Calculate OO and HO */
            if (NStoreO == 0) {
                if (AllComplexFlag == 0) {
                    calculateOO_real(SROptOO_real, SROptHO_real, SROptO_real, w, creal(e), SROptSize);
                } else {
                    calculateOO(SROptOO, SROptHO, SROptO, w, e, SROptSize);
                }
            } else {
                we = w * e;
                sqrtw = sqrt(w);
                if (AllComplexFlag == 0) {
#pragma omp parallel for default(shared) private(int_i)
                    for (int_i = 0; int_i < SROptSize; int_i++) {
                        // SROptO_Store for fortran
                        SROptO_Store_real[int_i + sample * SROptSize] = sqrtw * SROptO_real[int_i];
                        SROptHO_real[int_i] += creal(we) * SROptO_real[int_i];
                    }
                } else {
#pragma omp parallel for default(shared) private(int_i)
                    for (int_i = 0; int_i < SROptSize * 2; int_i++) {
                        // SROptO_Store for fortran
                        SROptO_Store[int_i + sample * (2 * SROptSize)] = sqrtw * SROptO[int_i];
                        SROptHO[int_i] += we * SROptO[int_i];
                    }
                }
            }
            StopTimer(43);

        } else if (NVMCCalMode == 1) {
            StartTimer(42);
            /* Calculate Green Function */
            CalculateGreenFuncBF(w, ip, eleIdx, eleCfg, eleNum, eleProjCnt, eleProjBFCnt);
            StopTimer(42);

            if (NLanczosMode > 0) {
                // ignoring Lanczos: to be added
                /* Calculate local QQQQ */
                StartTimer(43);
                if (AllComplexFlag == 0) {
                    LSLocalQ_real(creal(e), creal(ip), eleIdx, eleCfg, eleNum, eleProjCnt, LSLQ_real);
                    calculateQQQQ_real(QQQQ_real, LSLQ_real, w, NLSHam);
                } else {
                    LSLocalQ(e, ip, eleIdx, eleCfg, eleNum, eleProjCnt, LSLQ);
                    calculateQQQQ(QQQQ, LSLQ, w, NLSHam);
                    return;
                }
                StopTimer(43);
                if (NLanczosMode > 1) {
                    /* Calculate local QcisAjsQ */
                    StartTimer(44);
                    if (AllComplexFlag == 0) {
                        LSLocalCisAjs_real(creal(e), creal(ip), eleIdx, eleCfg, eleNum, eleProjCnt);
                        calculateQCAQ_real(QCisAjsQ_real, LSLCisAjs_real, LSLQ_real, w, NLSHam, NCisAjs);
                        calculateQCACAQ_real(QCisAjsCktAltQ_real, LSLCisAjs_real, w, NLSHam, NCisAjs,
                                             NCisAjsCktAltDC, CisAjsCktAltLzIdx);
                    } else {
                        LSLocalCisAjs(e, ip, eleIdx, eleCfg, eleNum, eleProjCnt);
                        calculateQCAQ(QCisAjsQ, LSLCisAjs, LSLQ, w, NLSHam, NCisAjs);
                        calculateQCACAQ(QCisAjsCktAltQ, LSLCisAjs, w, NLSHam, NCisAjs,
                                        NCisAjsCktAltDC, CisAjsCktAltLzIdx);
                        return;
                    }
                    StopTimer(44);
                }
            }
        }
    } /* end of for(sample) */

    // calculate OO and HO at NVMCCalMode==0
    if(NStoreO!=0 && NVMCCalMode==0){
        sampleSize=sampleEnd-sampleStart;
        if(AllComplexFlag==0){
            StartTimer(45);
            calculateOO_Store_real(SROptOO_real,SROptHO_real,SROptO_Store_real,creal(w),creal(e),SROptSize,sampleSize);
            StopTimer(45);
        }else{
            StartTimer(45);
            calculateOO_Store(SROptOO,SROptHO,SROptO_Store,w,e,2*SROptSize,sampleSize);
            StopTimer(45);
        }
    }

    InvM = InvM_Moto;
    PfM = PfM_Moto;

    return;
}


void clearPhysQuantity(){
  int i,n;
  double complex *vec;
  double  *vec_real;
//[s] MERGE BY TM
  Wc = Etot = Etot2 = Sztot =0.0;//fsz
  //Wc = Etot = Etot2 = 0.0;
  Dbtot = Dbtot2 = 0.0;
//[e] MERGE BY TM
  if(NVMCCalMode==0) {
    /* SROptOO, SROptHO, SROptO */
    n = (2*SROptSize)*(2*SROptSize+2); // TBC
    vec = SROptOO;
    #pragma omp parallel for default(shared) private(i)
    for(i=0;i<n;i++) vec[i] = 0.0+0.0*I;

// only for real variables
    n = (SROptSize)*(SROptSize+2); // TBC
    vec_real = SROptOO_real;
    #pragma omp parallel for default(shared) private(i)
    for(i=0;i<n;i++) vec_real[i] = 0.0;
  } else if(NVMCCalMode==1) {
    /* CisAjs, CisAjsCktAlt, CisAjsCktAltDC */
    n = 2*NCisAjs+NCisAjsCktAltDC+NCisAjsCktAltDC;

    vec = PhysCisAjs;
    #pragma omp parallel for default(shared) private(i)
    for(i=0;i<n;i++) vec[i] = 0.0+0.0*I;

    if(NLanczosMode>0) {
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
         n = NLSHam*NLSHam*NCisAjs + NLSHam*NLSHam*NCisAjsCktAltDC
          + NLSHam*NCisAjs;
        vec = QCisAjsQ;
#pragma omp parallel for default(shared) private(i)
        for(i=0;i<n;i++) vec[i] = 0.0+0.0*I;

        n = NLSHam*NLSHam*NCisAjs + NLSHam*NLSHam*NCisAjsCktAltDC
          + NLSHam*NCisAjs;
        vec_real = QCisAjsQ_real;
#pragma omp parallel for default(shared) private(i)
        for(i=0;i<n;i++) vec_real[i] = 0.0;

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

//#define M_DGEM dgemm_

extern int
dgemm_(char *jobz, char *uplo, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb,
       double *beta, double *c, int *ldc);

  char jobz, uplo;
  double alpha,beta;
  
  alpha = 1.0;
  beta  = 0.0;
  
  jobz = 'N';
  uplo = 'T';
  dgemm_(&jobz,&uplo,&srOptSize,&srOptSize,&sampleSize,&alpha,srOptO_Store_real,&srOptSize,srOptO_Store_real,&srOptSize,&beta,srOptOO_real,&srOptSize);

  return;
}



void calculateOO_Store(double complex *srOptOO, double complex *srOptHO, double complex *srOptO_Store,
                 const double w, const double complex e, int srOptSize, int sampleSize) {

  //#define M_DGEM dgemm_

  extern int zgemm_(char *jobz, char *uplo, int *m,int *n,int *k,double complex *alpha,  double complex *a, int *lda, double complex *b, int *ldb,
                    double complex *beta,double complex *c,int *ldc);

  char jobz, uplo;
  double complex alpha,beta;
  
  alpha = 1.0;
  beta  = 0.0;
  
  jobz = 'N';
  uplo = 'C';
  zgemm_(&jobz,&uplo,&srOptSize,&srOptSize,&sampleSize,&alpha,srOptO_Store,&srOptSize,srOptO_Store,&srOptSize,&beta,srOptOO,&srOptSize);

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

  #define M_ZAXPY zaxpy_
  #define M_ZGERC zgerc_

  extern int M_ZAXPY(const int *n, const double complex *alpha, const double complex *x, const int *incx,
                     double complex *y, const int *incy);
  extern int M_ZGERC(const int *m, const int *n, const double complex *alpha,
                    const double complex *x, const int *incx, const double complex *y, const int *incy, 
                    double complex *a, const int *lda);
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

  #define M_DAXPY daxpy_
  #define M_DGER dger_

  extern int M_DAXPY(const int *n, const double *alpha, const double *x, const int *incx,
                     double *y, const int *incy);
  extern int M_DGER(const int *m, const int *n, const double *alpha,
                    const double *x, const int *incx, const double *y, const int *incy, 
                    double *a, const int *lda);
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
#endif
