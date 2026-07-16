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
 * calculate Hamiltonian
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "calham_fsz.h"

#pragma once

double complex GreenFuncN_fsz(const int n, int *rsi, int *rsj,
                  const double complex ip, int *eleIdx, const int *eleCfg,
                  int *eleNum, const int *eleProjCnt, int *eleSpn,
                  double complex *buffer, int *bufferInt, double *rwork);

double complex GreenFunc1BF_fsz(const int ri, const int rj, const int s, const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer, int *affected,
                  int *hopIntWork, int hopIntWorkSize,
                  int *pfIWork, double *pfRWork,
                  double complex *pfBufM, double complex *pfWork);

double CalculateSz_fsz(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt,int *eleSpn) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  double  Sz=0.0;
  int ri;
    
  Sz = 0.0;
//#pragma omp for default(none) firstprivate(Nsite,n0,n1) private(idx,ri) reduction(+:Sz)
  for(ri=0;ri<Nsite;ri++) {
    Sz += n0[ri]-n1[ri];
  }
  return Sz;
}



double complex CalculateHamiltonian_fsz(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt,int *eleSpn) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  double complex e=0.0, tmp;
  int idx;
  int ri,rj,s,rk,rl,t,u,v;
  int k, offset, nbody;
  int *myEleIdx, *myEleNum, *myProjCntNew,*myEleSpn;
  int *myNBodyRsi, *myNBodyRsj;
  double complex *myBuffer;
  double *myRWork;
  double complex myEnergy;
  int maxNBodyInterAll = (NBodyInterAllMaxN > 2) ? NBodyInterAllMaxN : 2;
  int needNBodyInterAllRWork = (NNBodyInterAll > 0 && NBodyInterAllMaxN > 2);

  RequestWorkSpaceThreadInt(Nsize+Nsize+Nsite2+NProj+2*maxNBodyInterAll);
  RequestWorkSpaceThreadComplex(NQPFull+maxNBodyInterAll*Nsize);
  if (needNBodyInterAllRWork) RequestWorkSpaceThreadDouble(LapackLWork);
  /* GreenFunc1: NQPFull, GreenFunc2/NBodyInterAll: NQPFull+maxNBodyInterAll*Nsize */

  /*
#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myProjCntNew,myBuffer,myEnergy,idx)  \
  reduction(+:e)
  */
#pragma omp parallel default(none)                                      \
  private(myEleIdx,myEleSpn,myEleNum,myProjCntNew,myNBodyRsi,myNBodyRsj,myBuffer,myRWork,myEnergy, idx, ri, rj, rk, rl, s, t,u,v) \
  firstprivate(ip, Nsize, Nsite2, NProj, NQPFull, NCoulombIntra, CoulombIntra, ParaCoulombIntra, \
               NCoulombInter, CoulombInter, ParaCoulombInter, NHundCoupling, HundCoupling, ParaHundCoupling, \
               NTransfer, Transfer, ParaTransfer, NPairHopping, PairHopping, ParaPairHopping, \
               NExchangeCoupling, ExchangeCoupling, ParaExchangeCoupling, NInterAll, InterAll, ParaInterAll, \
               Nsite, maxNBodyInterAll, needNBodyInterAllRWork, LapackLWork, \
               NNBodyInterAll, NBodyInterAllN, NBodyInterAllOffset, NBodyInterAllIdx, ParaNBodyInterAll, n0, n1) \
  shared(eleCfg, eleProjCnt, eleIdx, eleNum,eleSpn) reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleSpn = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myNBodyRsi = GetWorkSpaceThreadInt(maxNBodyInterAll);
    myNBodyRsj = GetWorkSpaceThreadInt(maxNBodyInterAll);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull+maxNBodyInterAll*Nsize);
    myRWork = needNBodyInterAllRWork ? GetWorkSpaceThreadDouble(LapackLWork) : NULL;

    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleSpn[idx] = eleSpn[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
    #pragma omp barrier
    
    myEnergy = 0.0;

    #pragma omp master
    {StartTimer(70);}

    /* CoulombIntra */
    #pragma omp for private(idx,ri) nowait
    for(idx=0;idx<NCoulombIntra;idx++) {
      ri = CoulombIntra[idx];
      myEnergy += ParaCoulombIntra[idx] * n0[ri] * n1[ri];
    }
  
    /* CoulombInter */
    #pragma omp for private(idx,ri,rj) nowait
    for(idx=0;idx<NCoulombInter;idx++) {
      ri = CoulombInter[idx][0];
      rj = CoulombInter[idx][1];
      myEnergy += ParaCoulombInter[idx] * (n0[ri]+n1[ri]) * (n0[rj]+n1[rj]);
    }

    /* HundCoupling */
    #pragma omp for private(idx,ri,rj) nowait
    for(idx=0;idx<NHundCoupling;idx++) {
      ri = HundCoupling[idx][0];
      rj = HundCoupling[idx][1];
      myEnergy -= ParaHundCoupling[idx] * (n0[ri]*n0[rj] + n1[ri]*n1[rj]);
      /* Caution: negative sign */
    }

    #pragma omp master
    {StopTimer(70);StartTimer(71);}

    /* Transfer */
    #pragma omp for private(idx,ri,rj,s,t) schedule(dynamic) nowait
    for(idx=0;idx<NTransfer;idx++) {
      ri = Transfer[idx][0];
      s  = Transfer[idx][1];
      rj = Transfer[idx][2];
      t  = Transfer[idx][3];
      if(s==t){
        myEnergy -= ParaTransfer[idx]
         * GreenFunc1_fsz(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      }else{ 
        myEnergy -= ParaTransfer[idx]
        * GreenFunc1_fsz2(ri,rj,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      }
      /* Caution: negative sign */
    }

    #pragma omp master
    {StopTimer(71);StartTimer(72);}

    /* Pair Hopping */
    #pragma omp for private(idx,ri,rj) schedule(dynamic) nowait
    for(idx=0;idx<NPairHopping;idx++) {
      ri = PairHopping[idx][0];
      rj = PairHopping[idx][1];
    
      myEnergy += ParaPairHopping[idx]
        * GreenFunc2_fsz(ri,rj,ri,rj,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
    }

    /* Exchange Coupling */
    #pragma omp for private(idx,ri,rj,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NExchangeCoupling;idx++) {
      ri = ExchangeCoupling[idx][0];
      rj = ExchangeCoupling[idx][1];
    
      tmp =  GreenFunc2_fsz(ri,rj,rj,ri,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      //printf("idx=%d ri=%d rj=%d:  tmp=%lf \n",idx,ri,rj,creal(tmp));
      tmp += GreenFunc2_fsz(ri,rj,rj,ri,1,0,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      //tmp =  GreenFunc2_fsz2(ri,rj,rj,ri,0,0,1,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      //tmp += GreenFunc2_fsz2(ri,rj,rj,ri,1,1,0,0,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      myEnergy += ParaExchangeCoupling[idx] * tmp;
    }

    /* Inter All */
    #pragma omp for private(idx,ri,rj,s,rk,rl,t) schedule(dynamic) nowait
    for(idx=0;idx<NInterAll;idx++) {
      ri = InterAll[idx][0];
      s  = InterAll[idx][1];
      rj = InterAll[idx][2];
      t  = InterAll[idx][3];
      rk = InterAll[idx][4];
      u  = InterAll[idx][5];
      rl = InterAll[idx][6];
      v  = InterAll[idx][7];
      
      if(s==t && u==v){
        myEnergy += ParaInterAll[idx]
          * GreenFunc2_fsz(ri,rj,rk,rl,s,u,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      }else{
        myEnergy += ParaInterAll[idx]
          * GreenFunc2_fsz2(ri,rj,rk,rl,s,t,u,v,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      } 
    }

    /* NBodyInterAll */
    #pragma omp for private(idx,k,offset,nbody,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NNBodyInterAll;idx++) {
      nbody = NBodyInterAllN[idx];
      offset = NBodyInterAllOffset[idx];
      for(k=0;k<nbody;k++) {
        myNBodyRsi[k] = NBodyInterAllIdx[offset+k][0]
                      + NBodyInterAllIdx[offset+k][1]*Nsite;
        myNBodyRsj[k] = NBodyInterAllIdx[offset+k][2]
                      + NBodyInterAllIdx[offset+k][3]*Nsite;
      }
      tmp = GreenFuncN_fsz(nbody,myNBodyRsi,myNBodyRsj,ip,myEleIdx,eleCfg,
                           myEleNum,eleProjCnt,myEleSpn,myBuffer,
                           myProjCntNew,myRWork);
      myEnergy += ParaNBodyInterAll[idx] * tmp;
    }

    #pragma omp master
    {StopTimer(72);}

    e += myEnergy;
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  if (needNBodyInterAllRWork) ReleaseWorkSpaceThreadDouble();
  return e;
}

double complex CalculateHamiltonianBF_fsz(const double complex ip, int *eleIdx, int *eleCfg,
                             int *eleNum, const int *eleProjCnt, int *eleSpn,
                             const int *eleProjBFCnt) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  double complex e=0.0;
  int idx;
  int ri,rj,s,t;
  int *myEleIdx, *myEleCfg, *myEleNum, *myEleSpn;
  int *myProjCntNew, *myProjBFCntNew, *myAffected, *myHopIntWork;
  int *myPfIWork;
  double complex *myBuffer, *myPfBufM, *myPfWork;
  double *myPfRWork;
  const size_t bfSlaterSize = (size_t)NQPFull*(size_t)Nsite2*(size_t)Nsite2;
  size_t pfUpdateComplexCount,pfUpdateIntCount,pfUpdateDoubleCount;
  size_t bfBufferSize,bfComplexSize,nsizeSquared;
  long long bfBaseIntSizeLL,bfIntSizeLL;
  int bfBaseIntSize, bfHopIntSize, bfIntSize, bfPfIntSize, bfPfDoubleSize;

  if(NPairHopping > 0 || NExchangeCoupling > 0 || NInterAll > 0 || NNBodyInterAll > 0) {
    fprintf(stderr, "Error: CalculateHamiltonianBF_fsz supports only density and one-body transfer terms.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(GetSlaterElmBF_fsz_hop_int_work_size(&bfHopIntSize) != 0) {
    fprintf(stderr, "Error: invalid BF-FSZ Hamiltonian hop integer workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(GetCalculateNewPfMBF_fsz_rows_work_size(&pfUpdateComplexCount,
      &pfUpdateIntCount,&pfUpdateDoubleCount) != 0
      || pfUpdateIntCount > INT_MAX || pfUpdateDoubleCount > INT_MAX
      || bfSlaterSize > SIZE_MAX-2*(size_t)NQPFull
      || bfSlaterSize+2*(size_t)NQPFull > SIZE_MAX-pfUpdateComplexCount) {
    fprintf(stderr,"Error: invalid BF-FSZ Hamiltonian Pfaffian workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  nsizeSquared = (size_t)Nsize*(size_t)Nsize;
  bfBufferSize = 2*(size_t)NQPFull + bfSlaterSize + pfUpdateComplexCount;
  if(nsizeSquared > SIZE_MAX-(size_t)LapackLWork
      || bfBufferSize > SIZE_MAX-nsizeSquared-(size_t)LapackLWork) {
    fprintf(stderr,"Error: BF-FSZ Hamiltonian complex workspace size overflow.\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  bfComplexSize = bfBufferSize+nsizeSquared+(size_t)LapackLWork;
  if(bfComplexSize > INT_MAX) {
    fprintf(stderr,"Error: BF-FSZ Hamiltonian complex workspace is too large.\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  bfPfIntSize = (Nsize > (int)pfUpdateIntCount)
      ? Nsize : (int)pfUpdateIntCount;
  bfPfDoubleSize = (LapackLWork > (int)pfUpdateDoubleCount)
      ? LapackLWork : (int)pfUpdateDoubleCount;
  bfBaseIntSizeLL = 2LL*Nsize + 2LL*Nsite2 + (long long)NProj
      + 16LL*Nsite*Nrange;
  if(bfBaseIntSizeLL < 0 || bfBaseIntSizeLL > INT_MAX) {
    fprintf(stderr, "Error: invalid BF-FSZ Hamiltonian base integer workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  bfBaseIntSize = (int)bfBaseIntSizeLL;
  bfIntSizeLL = (long long)bfBaseIntSize+Nsize+bfHopIntSize+bfPfIntSize;
  if(bfIntSizeLL < 0 || bfIntSizeLL > INT_MAX) {
    fprintf(stderr, "Error: invalid BF-FSZ Hamiltonian integer workspace size.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  bfIntSize = (int)bfIntSizeLL;

  RequestWorkSpaceInt(bfIntSize);
  RequestWorkSpaceComplex((int)bfComplexSize);
  RequestWorkSpaceDouble(bfPfDoubleSize);

  myEleIdx = GetWorkSpaceInt(Nsize);
  myEleCfg = GetWorkSpaceInt(Nsite2);
  myEleNum = GetWorkSpaceInt(Nsite2);
  myEleSpn = GetWorkSpaceInt(Nsize);
  myProjCntNew = GetWorkSpaceInt(NProj);
  myProjBFCntNew = GetWorkSpaceInt(16*Nsite*Nrange);
  myAffected = GetWorkSpaceInt(Nsize);
  myHopIntWork = GetWorkSpaceInt(bfHopIntSize);
  myPfIWork = GetWorkSpaceInt(bfPfIntSize);
  myBuffer = GetWorkSpaceComplex((int)bfBufferSize);
  myPfBufM = GetWorkSpaceComplex(Nsize*Nsize);
  myPfWork = GetWorkSpaceComplex(LapackLWork);
  myPfRWork = GetWorkSpaceDouble(bfPfDoubleSize);

  for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
  for(idx=0;idx<Nsite2;idx++) myEleCfg[idx] = eleCfg[idx];
  for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
  for(idx=0;idx<Nsize;idx++) myEleSpn[idx] = eleSpn[idx];

  for(idx=0;idx<NCoulombIntra;idx++) {
    ri = CoulombIntra[idx];
    e += ParaCoulombIntra[idx] * n0[ri] * n1[ri];
  }

  for(idx=0;idx<NCoulombInter;idx++) {
    ri = CoulombInter[idx][0];
    rj = CoulombInter[idx][1];
    e += ParaCoulombInter[idx] * (n0[ri]+n1[ri]) * (n0[rj]+n1[rj]);
  }

  for(idx=0;idx<NHundCoupling;idx++) {
    ri = HundCoupling[idx][0];
    rj = HundCoupling[idx][1];
    e -= ParaHundCoupling[idx] * (n0[ri]*n0[rj] + n1[ri]*n1[rj]);
  }

  for(idx=0;idx<NTransfer;idx++) {
    ri = Transfer[idx][0];
    s  = Transfer[idx][1];
    rj = Transfer[idx][2];
    t  = Transfer[idx][3];
    if(s != t) {
      fprintf(stderr, "Error: BackFlow FSZ does not support spin-changing transfer terms yet.\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    e -= ParaTransfer[idx]
      * GreenFunc1BF_fsz(ri,rj,s,ip,myEleIdx,myEleCfg,myEleNum,eleProjCnt,myEleSpn,
                         myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer,
                         myAffected,myHopIntWork,bfHopIntSize,
                         myPfIWork,myPfRWork,myPfBufM,myPfWork);
  }

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceComplex();
  ReleaseWorkSpaceDouble();
  return e;
}

/* Calculate the CoulombIntra, CoulombInter, Hund terms, */
/* which can be calculated by number operators. */
/* This function will be used in the Lanczos mode */
double complex CalculateHamiltonian0_fsz(const int *eleNum) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  double complex e=0.0;
  int idx;
  int ri,rj;
  double complex myEnergy;

#pragma omp parallel default(shared)\
  private(myEnergy) reduction(+:e)
  {
    myEnergy = 0.0;

    /* CoulombIntra */
    #pragma omp for private(idx,ri)
    for(idx=0;idx<NCoulombIntra;idx++) {
      ri = CoulombIntra[idx];
      myEnergy += ParaCoulombIntra[idx] * n0[ri] * n1[ri];
    }
  
    /* CoulombInter */
    #pragma omp for private(idx,ri,rj)
    for(idx=0;idx<NCoulombInter;idx++) {
      ri = CoulombInter[idx][0];
      rj = CoulombInter[idx][1];
      myEnergy += ParaCoulombInter[idx] * (n0[ri]+n1[ri]) * (n0[rj]+n1[rj]);
    }

    /* HundCoupling */
    #pragma omp for private(idx,ri,rj)
    for(idx=0;idx<NHundCoupling;idx++) {
      ri = HundCoupling[idx][0];
      rj = HundCoupling[idx][1];
      myEnergy -= ParaHundCoupling[idx] * (n0[ri]*n0[rj] + n1[ri]*n1[rj]);
      /* Caution: negative sign */
    }

    e += myEnergy;
  }

  return e;
}

/* Calculate the transfer terms, */
/* which can be calculated by 1-body Green function. */
/* This function will be used in the Lanczos mode */
double complex CalculateHamiltonian1_fsz(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt,int *eleSpn) {
  double complex e=0.0;
  int idx;
  int ri,rj,s;
  int *myEleIdx, *myEleNum, *myProjCntNew,*myEleSpn;
  double complex *myBuffer;
  double complex myEnergy;

  RequestWorkSpaceThreadInt(Nsize+Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadComplex(NQPFull);
  /* GreenFunc1: NQPFull */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleSpn,myEleNum,myProjCntNew,myBuffer,myEnergy,idx)  \
  reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleSpn = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull);

    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleSpn[idx] = eleSpn[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];

    myEnergy = 0.0;

    /* Transfer */
    #pragma omp for private(idx,ri,rj,s) schedule(dynamic) nowait
    for(idx=0;idx<NTransfer;idx++) {
      ri = Transfer[idx][0];
      s  = Transfer[idx][1];
      rj = Transfer[idx][2];
      
      myEnergy -= ParaTransfer[idx]
        * GreenFunc1_fsz(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      /* Caution: negative sign */
    }

    e += myEnergy;
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  return e;
}

/* Calculate the exchange coupling, pair hopping, interAll terms, */
/* which can be calculated by 2-body Green function. */
/* This function will be used in the Lanczos mode */
double complex CalculateHamiltonian2_fsz(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt,int *eleSpn) {
  double e=0.0, tmp;
  int idx;
  int ri,rj,s,rk,rl,t;
  int *myEleIdx, *myEleNum, *myProjCntNew,*myEleSpn;
  double complex *myBuffer;
  double complex myEnergy;

  RequestWorkSpaceThreadInt(Nsize+Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadComplex(NQPFull+2*Nsize);
  /* GreenFunc2: NQPFull+2*Nsize */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleSpn,myEleNum,myProjCntNew,myBuffer,myEnergy,idx)  \
  reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleSpn = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull+2*Nsize);

    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleSpn[idx] = eleSpn[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];

    myEnergy = 0.0;

    /* Pair Hopping */
    #pragma omp for private(idx,ri,rj) schedule(dynamic) nowait
    for(idx=0;idx<NPairHopping;idx++) {
      ri = PairHopping[idx][0];
      rj = PairHopping[idx][1];
    
      myEnergy += ParaPairHopping[idx]
        * GreenFunc2_fsz(ri,rj,ri,rj,0,1,ip,myEleIdx,eleCfg,myEleNum,
                     eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
    }

    /* Exchange Coupling */
    #pragma omp for private(idx,ri,rj,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NExchangeCoupling;idx++) {
      ri = ExchangeCoupling[idx][0];
      rj = ExchangeCoupling[idx][1];
    
      tmp =  GreenFunc2_fsz(ri,rj,rj,ri,0,1,ip,myEleIdx,eleCfg,myEleNum,
                        eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      tmp += GreenFunc2_fsz(ri,rj,rj,ri,1,0,ip,myEleIdx,eleCfg,myEleNum,
                        eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      myEnergy += ParaExchangeCoupling[idx] * tmp;
    }

    /* Inter All */
    #pragma omp for private(idx,ri,rj,s,rk,rl,t) schedule(dynamic) nowait
    for(idx=0;idx<NInterAll;idx++) {
      ri = InterAll[idx][0];
      rj = InterAll[idx][2];
      s  = InterAll[idx][3];
      rk = InterAll[idx][4];
      rl = InterAll[idx][6];
      t  = InterAll[idx][7];
      
      myEnergy += ParaInterAll[idx]
        * GreenFunc2_fsz(ri,rj,rk,rl,s,t,ip,myEleIdx,eleCfg,myEleNum,
                     eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
    }

    e += myEnergy;
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  return e;
}
