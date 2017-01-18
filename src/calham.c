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
 * calculate Hamiltonian
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#pragma once

double complex CalculateHamiltonian(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt);
double complex CalculateHamiltonian0(const int *eleNum);
double complex CalculateHamiltonian1(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt);
double complex CalculateHamiltonian2(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt);

double CalculateDoubleOccupation(int *eleIdx, const int *eleCfg,
                                 int *eleNum, const int *eleProjCnt) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  double db=0.0;
  int ri;

#pragma omp master
  {StartTimer(70);}

  /* CoulombIntra */
#pragma omp parallel for private(ri) reduction(+:db)
  for(ri=0;ri<Nsite;ri++) {
    db += (double)(n0[ri] * n1[ri]);
  }

#pragma omp master
  {StopTimer(70);}

  return db;
}

double complex CalculateHamiltonian(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  double complex e=0.0, tmp;
  int idx;
  int ri,rj,s,rk,rl,t;
  int *myEleIdx, *myEleNum, *myProjCntNew;
  double complex *myBuffer;
  double complex myEnergy;

  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadComplex(NQPFull+2*Nsize);
  /* GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize */

  /*
#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myProjCntNew,myBuffer,myEnergy,idx)  \
  reduction(+:e)
  */
#pragma omp parallel default(none)                                      \
  private(myEleIdx,myEleNum,myProjCntNew,myBuffer,myEnergy, idx, ri, rj, rk, rl, s, t) \
  firstprivate(ip, Nsize, Nsite2, NProj, NQPFull, NCoulombIntra, CoulombIntra, ParaCoulombIntra, \
               NCoulombInter, CoulombInter, ParaCoulombInter, NHundCoupling, HundCoupling, ParaHundCoupling, \
               NTransfer, Transfer, ParaTransfer, NPairHopping, PairHopping, ParaPairHopping, \
               NExchangeCoupling, ExchangeCoupling, ParaExchangeCoupling, NInterAll, InterAll, ParaInterAll, n0, n1) \
  shared(eleCfg, eleProjCnt, eleIdx, eleNum) reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull+2*Nsize);

    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
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
    #pragma omp for private(idx,ri,rj,s) schedule(dynamic) nowait
    for(idx=0;idx<NTransfer;idx++) {
      ri = Transfer[idx][0];
      rj = Transfer[idx][2];
      s  = Transfer[idx][3];
      
      myEnergy -= ParaTransfer[idx]
        * GreenFunc1(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
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
        * GreenFunc2(ri,rj,ri,rj,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
    }

    /* Exchange Coupling */
    #pragma omp for private(idx,ri,rj,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NExchangeCoupling;idx++) {
      ri = ExchangeCoupling[idx][0];
      rj = ExchangeCoupling[idx][1];
    
      tmp =  GreenFunc2(ri,rj,rj,ri,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
      tmp += GreenFunc2(ri,rj,rj,ri,1,0,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
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
        * GreenFunc2(ri,rj,rk,rl,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
    }

    #pragma omp master
    {StopTimer(72);}

    e += myEnergy;
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  return e;
}

/* Calculate the CoulombIntra, CoulombInter, Hund terms, */
/* which can be calculated by number operators. */
/* This function will be used in the Lanczos mode */
double complex CalculateHamiltonian0(const int *eleNum) {
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
double complex CalculateHamiltonian1(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt) {
  double complex e=0.0;
  int idx;
  int ri,rj,s;
  int *myEleIdx, *myEleNum, *myProjCntNew;
  double complex *myBuffer;
  double complex myEnergy;

  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadComplex(NQPFull);
  /* GreenFunc1: NQPFull */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myProjCntNew,myBuffer,myEnergy,idx)  \
  reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull);

    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];

    myEnergy = 0.0;

    /* Transfer */
    #pragma omp for private(idx,ri,rj,s) schedule(dynamic) nowait
    for(idx=0;idx<NTransfer;idx++) {
      ri = Transfer[idx][0];
      rj = Transfer[idx][2];
      s  = Transfer[idx][3];
      
      myEnergy -= ParaTransfer[idx]
        * GreenFunc1(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
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
double complex CalculateHamiltonian2(const double complex ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt) {
  double e=0.0, tmp;
  int idx;
  int ri,rj,s,rk,rl,t;
  int *myEleIdx, *myEleNum, *myProjCntNew;
  double complex *myBuffer;
  double complex myEnergy;

  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadComplex(NQPFull+2*Nsize);
  /* GreenFunc2: NQPFull+2*Nsize */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myProjCntNew,myBuffer,myEnergy,idx)  \
  reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull+2*Nsize);

    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];

    myEnergy = 0.0;

    /* Pair Hopping */
    #pragma omp for private(idx,ri,rj) schedule(dynamic) nowait
    for(idx=0;idx<NPairHopping;idx++) {
      ri = PairHopping[idx][0];
      rj = PairHopping[idx][1];
    
      myEnergy += ParaPairHopping[idx]
        * GreenFunc2(ri,rj,ri,rj,0,1,ip,myEleIdx,eleCfg,myEleNum,
                     eleProjCnt,myProjCntNew,myBuffer);
    }

    /* Exchange Coupling */
    #pragma omp for private(idx,ri,rj,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NExchangeCoupling;idx++) {
      ri = ExchangeCoupling[idx][0];
      rj = ExchangeCoupling[idx][1];
    
      tmp =  GreenFunc2(ri,rj,rj,ri,0,1,ip,myEleIdx,eleCfg,myEleNum,
                        eleProjCnt,myProjCntNew,myBuffer);
      tmp += GreenFunc2(ri,rj,rj,ri,1,0,ip,myEleIdx,eleCfg,myEleNum,
                        eleProjCnt,myProjCntNew,myBuffer);
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
        * GreenFunc2(ri,rj,rk,rl,s,t,ip,myEleIdx,eleCfg,myEleNum,
                     eleProjCnt,myProjCntNew,myBuffer);
    }

    e += myEnergy;
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  return e;
}
