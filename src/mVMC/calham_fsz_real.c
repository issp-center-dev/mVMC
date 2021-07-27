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
 *-------------------------------------------------------------*/

#include "calham_fsz_real.h"

#pragma once

double CalculateHamiltonian_fsz_real(const double ip, int *eleIdx, const int *eleCfg,
                             int *eleNum, const int *eleProjCnt,int *eleSpn) {
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  double e=0.0, tmp;
  int idx;
  int ri,rj,s,rk,rl,t,u,v;
  int *myEleIdx, *myEleNum, *myProjCntNew,*myEleSpn;
  double *myBuffer;
  double  myEnergy;

  RequestWorkSpaceThreadInt(Nsize+Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadDouble(NQPFull+2*Nsize);
  /* GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize */

  /*
#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myProjCntNew,myBuffer,myEnergy,idx)  \
  reduction(+:e)
  */
#pragma omp parallel default(none)                                      \
  private(myEleIdx,myEleSpn,myEleNum,myProjCntNew,myBuffer,myEnergy, idx, ri, rj, rk, rl, s, t,u,v) \
  firstprivate(ip, Nsize, Nsite2, NProj, NQPFull, NCoulombIntra, CoulombIntra, ParaCoulombIntra, \
               NCoulombInter, CoulombInter, ParaCoulombInter, NHundCoupling, HundCoupling, ParaHundCoupling, \
               NTransfer, Transfer, ParaTransfer, NPairHopping, PairHopping, ParaPairHopping, \
               NExchangeCoupling, ExchangeCoupling, ParaExchangeCoupling, NInterAll, InterAll, ParaInterAll, n0, n1, ip) \
  shared(eleCfg, eleProjCnt, eleIdx, eleNum,eleSpn) reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleSpn = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myBuffer = GetWorkSpaceThreadDouble(NQPFull+2*Nsize);

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
         * GreenFunc1_fsz_real(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      }else{ 
        myEnergy -= ParaTransfer[idx]
        * GreenFunc1_fsz2_real(ri,rj,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
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
        * GreenFunc2_fsz_real(ri,rj,ri,rj,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
    }

    /* Exchange Coupling */
    #pragma omp for private(idx,ri,rj,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NExchangeCoupling;idx++) {
      ri = ExchangeCoupling[idx][0];
      rj = ExchangeCoupling[idx][1];
    
      tmp =  GreenFunc2_fsz_real(ri,rj,rj,ri,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      tmp += GreenFunc2_fsz_real(ri,rj,rj,ri,1,0,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
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
          * GreenFunc2_fsz_real(ri,rj,rk,rl,s,u,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      }else{
        myEnergy += ParaInterAll[idx]
          * GreenFunc2_fsz2_real(ri,rj,rk,rl,s,t,u,v,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,myProjCntNew,myBuffer);
      } 
    }

    #pragma omp master
    {StopTimer(72);}

    e += myEnergy;
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadDouble();
  return e;
}


