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
#ifndef _SRC_CALHAMREAL
#define _SRC_CALHAMREAL
#include "workspace.c"
#include "complex.h"
#include "slater.h"
#include "locgrn_real.h"

double CalculateHamiltonian_real(const double ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt);

double CalculateHamiltonianBF_real(const double ip, int *eleIdx, const int *eleCfg,
                              int *eleNum, const int *eleProjCnt, const int *eleProjBFCnt);


double CalculateHamiltonian0_real(const int *eleNum);

///
/// \param ip
/// \param eleIdx
/// \param eleCfg
/// \param eleNum
/// \param eleProjCnt
/// \return myEnergy
double CalculateHamiltonian_real(const double ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  double e=0.0, tmp;
  int idx;
  int ri,rj,s,rk,rl,t;
  int *myEleIdx, *myEleNum, *myProjCntNew;
  double  *myBuffer;
  double  myEnergy;

  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadDouble(NQPFull+2*Nsize);
  /* GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize */

  #pragma omp parallel default(none) \
    private(myEleIdx,myEleNum,myProjCntNew,myBuffer,myEnergy, idx, ri, rj, rk, rl, s, t) \
    firstprivate(ip, Nsize, Nsite2, NProj, NQPFull, NCoulombIntra, CoulombIntra, ParaCoulombIntra,   \
    NCoulombInter, CoulombInter, ParaCoulombInter, NHundCoupling, HundCoupling, ParaHundCoupling,    \
    NTransfer, Transfer, ParaTransfer, NPairHopping, PairHopping, ParaPairHopping,    \
    NExchangeCoupling, ExchangeCoupling, ParaExchangeCoupling, NInterAll, InterAll, ParaInterAll, n0, n1)\
    shared(eleCfg, eleProjCnt, eleIdx, eleNum) reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myBuffer = GetWorkSpaceThreadDouble(NQPFull+2*Nsize);

    #pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
    #pragma omp barrier
    
    myEnergy = 0.0;

    #pragma omp master
    {StartTimer(70);}
#ifdef _DEBUG
#pragma omp master
    printf("    Debug: CoulombIntra\n");
#endif
    /* CoulombIntra */
    #pragma omp for private(idx,ri) nowait
    for(idx=0;idx<NCoulombIntra;idx++) {
      ri = CoulombIntra[idx];
      myEnergy += ParaCoulombIntra[idx] * n0[ri] * n1[ri];
    }

#ifdef _DEBUG
#pragma omp master
    printf("    Debug: CoulombInter\n");
#endif
    /* CoulombInter */
    #pragma omp for private(idx,ri,rj) nowait
    for(idx=0;idx<NCoulombInter;idx++) {
      ri = CoulombInter[idx][0];
      rj = CoulombInter[idx][1];
      myEnergy += ParaCoulombInter[idx] * (n0[ri]+n1[ri]) * (n0[rj]+n1[rj]);
    }

#ifdef _DEBUG
#pragma omp master
    printf("    Debug: HundCoupling\n");
#endif
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

#ifdef _DEBUG
#pragma omp master
    printf("    Debug: Transfer\n");
#endif
    /* Transfer */
#pragma omp for private(idx,ri,rj,s) schedule(dynamic) nowait
    for(idx=0;idx<NTransfer;idx++) {
      ri = Transfer[idx][0];
      rj = Transfer[idx][2];
      s  = Transfer[idx][3];
      
      myEnergy -= creal(ParaTransfer[idx])
        * GreenFunc1_real(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
      /* Caution: negative sign */
    }

    #pragma omp master
    {StopTimer(71);StartTimer(72);}

#ifdef _DEBUG
#pragma omp master
    printf("    Debug: PairHopping\n");
#endif
    /* Pair Hopping */
    #pragma omp for private(idx,ri,rj) schedule(dynamic) nowait
    for(idx=0;idx<NPairHopping;idx++) {
      ri = PairHopping[idx][0];
      rj = PairHopping[idx][1];
    
      myEnergy += ParaPairHopping[idx]
        * GreenFunc2_real(ri,rj,ri,rj,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
    }

#ifdef _DEBUG
#pragma omp master
    printf("    Debug: ExchangeCoupling, NExchangeCoupling=%d\n", NExchangeCoulpling);
#endif
    /* Exchange Coupling */
    #pragma omp for private(idx,ri,rj,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NExchangeCoupling;idx++) {
      ri = ExchangeCoupling[idx][0];
      rj = ExchangeCoupling[idx][1];
    
      tmp =  GreenFunc2_real(ri,rj,rj,ri,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
      tmp += GreenFunc2_real(ri,rj,rj,ri,1,0,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
      myEnergy += ParaExchangeCoupling[idx] * tmp;
    }
    
#ifdef _DEBUG
#pragma omp master
    printf("    Debug: InterAll, NInterAll=%d\n", NInterAll);
#endif
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
        * GreenFunc2_real(ri,rj,rk,rl,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
    }

    #pragma omp master
    {StopTimer(72);}
#ifdef _DEBUG    
    printf("    Debug: myEnergy=%lf\n", myEnergy);
#endif
    e += myEnergy;
  }
#ifdef _DEBUG
#pragma omp master
  printf("    Debug: Release\n");
#endif
  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadDouble();
#ifdef _DEBUG
#pragma omp master
  printf("    Debug: HamRealFinish\n", NInterAll);
#endif
  return e;
}

double CalculateHamiltonianBF_real(const double ip, int *eleIdx, const int *eleCfg,
                              int *eleNum, const int *eleProjCnt, const int *eleProjBFCnt) {
    const int *n0 = eleNum;
    const int *n1 = eleNum + Nsite;
    double e=0.0, tmp;
    int idx;
    int ri,rj,s,rk,rl,t;
    int *myEleIdx, *myEleNum, *myEleCfg, *myProjCntNew, *myProjBFCntNew;
    //double sltTmp[NThread*NQPFull*Nsite2*Nsite2];
    double *mySltBFTmp;
    double *myBuffer;
    double myEnergy;

    RequestWorkSpaceThreadInt(Nsize+2*Nsite2+NProj+16*Nsite*Nrange);
    RequestWorkSpaceThreadDouble(NQPFull+2*Nsize+NQPFull*Nsite2*Nsite2);
    /* GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myEleCfg,myProjCntNew,myProjBFCntNew,myBuffer,mySltBFTmp,myEnergy)\
  reduction(+:e)
    {
        myEleIdx = GetWorkSpaceThreadInt(Nsize);
        myEleNum = GetWorkSpaceThreadInt(Nsite2);
        myEleCfg = GetWorkSpaceThreadInt(Nsite2);
        myProjCntNew   = GetWorkSpaceThreadInt(NProj);
        myProjBFCntNew = GetWorkSpaceThreadInt(16*Nsite*Nrange);
        myBuffer   = GetWorkSpaceThreadDouble(NQPFull+2*Nsize);
        mySltBFTmp = GetWorkSpaceThreadDouble(NQPFull*Nsite2*Nsite2);

#pragma loop noalias
        for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
#pragma loop noalias
        for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
#pragma loop noalias
        for(idx=0;idx<Nsite2;idx++) myEleCfg[idx] = eleCfg[idx];

        StoreSlaterElmBF_real(mySltBFTmp);
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

            myEnergy -= creal(ParaTransfer[idx])
                        //* GreenFunc1(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
                        * GreenFunc1BF_real(ri,rj,s,ip,mySltBFTmp,myEleIdx,myEleCfg,myEleNum,eleProjCnt,myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer);
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
                        * GreenFunc2_real(ri,rj,ri,rj,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
        }

        /* Exchange Coupling */
#pragma omp for private(idx,ri,rj,tmp) schedule(dynamic) nowait
        for(idx=0;idx<NExchangeCoupling;idx++) {
            ri = ExchangeCoupling[idx][0];
            rj = ExchangeCoupling[idx][1];

            tmp =  GreenFunc2_real(ri,rj,rj,ri,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
            tmp += GreenFunc2_real(ri,rj,rj,ri,1,0,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
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
            myEnergy += creal(ParaInterAll[idx])
                        * GreenFunc2_real(ri,rj,rk,rl,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
        }

#pragma omp master
        {StopTimer(72);}

        e += myEnergy;
    }

    ReleaseWorkSpaceThreadInt();
    ReleaseWorkSpaceThreadDouble();
    return e;
}

/* Calculate the CoulombIntra, CoulombInter, Hund terms, */
/* which can be calculated by number operators. */
///
/// \param eleNum [in]
/// \return myEnergy
/// \version 1.0
double CalculateHamiltonian0_real(const int *eleNum)
{
    const int *n0 = eleNum;
    const int *n1 = eleNum + Nsite;
    double e=0.0;
    int idx;
    int ri,rj;
    double myEnergy;

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

#endif
