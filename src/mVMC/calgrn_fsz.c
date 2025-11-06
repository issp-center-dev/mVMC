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
 * Cauculate Green Functions
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "calgrn_fsz.h"
#ifndef _CALGRN_FSZ_SRC
#define _CALGRN_FSZ_SRC


void CalculateGreenFunc_fsz(const double w, const double complex ip, int *eleIdx, int *eleCfg,
                         int *eleNum, int *eleSpn,int *eleProjCnt) {

  int idx,idx0,idx1;
  int ri,rj,s,rk,rl,t,u,v;
  double complex tmp;
  int *myEleIdx, *myEleNum, *myProjCntNew,*myEleSpn;
  double complex *myBuffer;
  double two_pi = 2.0*acos(-1.0);
  int site_idx, spin_idx, rsi_idx, rsi;
  int x, y, z;
  double Wx, Wy, Wz; 
  double weight;

  RequestWorkSpaceThreadInt(Nsize+Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadComplex(NQPFull+2*Nsize);
  /* GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize */
  
  for(idx=0;idx<NTwist;idx++) {
    weight = 0.0;
    for(site_idx=0;site_idx<2*Nsite;site_idx++) {
      ri  = TwistIdx[idx][2*site_idx];
      s   = TwistIdx[idx][2*site_idx+1];
      rsi = ri + s*Nsite;
      
      x = LatticeIdx[ri][0];
      y = LatticeIdx[ri][1];
      z = LatticeIdx[ri][2];
      
      Wx = ParaTwist[idx][3*site_idx + 0];
      Wy = ParaTwist[idx][3*site_idx + 1];
      Wz = ParaTwist[idx][3*site_idx + 2];
      tmp = x*Wx + y*Wy + z*Wz; 
      
      weight += tmp*(double)eleNum[rsi];
      //printf("weight = %lf x=%d y=%d z=%d Wx=%d Wy=%d Wz=%d eleNum=%d\n",
      //       weight,x,y,z,Wx,Wy,Wz,eleNum[rsi]);

    }
    PhysTwist[idx] += w*cexp(I*two_pi*weight);
  }

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myEleSpn,myProjCntNew,myBuffer,idx)
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

    #pragma omp master
    {StartTimer(50);}

    #pragma omp for private(idx,ri,rj,s,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NCisAjs;idx++) {
      ri = CisAjsIdx[idx][0];
      s  = CisAjsIdx[idx][1];
      rj = CisAjsIdx[idx][2];
      t  = CisAjsIdx[idx][3];
      if(s==t){
        tmp = GreenFunc1_fsz(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,
                       myProjCntNew,myBuffer);
      }else{
        tmp = GreenFunc1_fsz2(ri,rj,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,
                       myProjCntNew,myBuffer);
      }     
      LocalCisAjs[idx] = tmp;
    }

    #pragma omp master
    {StopTimer(50);StartTimer(51);}
    
    #pragma omp for private(idx,ri,rj,s,rk,rl,t,tmp) schedule(dynamic)
    for(idx=0;idx<NCisAjsCktAltDC;idx++) {
      /*
      ri = CisAjsCktAltDCIdx[idx][0];
      rj = CisAjsCktAltDCIdx[idx][1];
      s  = CisAjsCktAltDCIdx[idx][2];
      rk = CisAjsCktAltDCIdx[idx][3];
      rl = CisAjsCktAltDCIdx[idx][4];
      t  = CisAjsCktAltDCIdx[idx][5];
      */
      ri = CisAjsCktAltDCIdx[idx][0];
      s  = CisAjsCktAltDCIdx[idx][1];
      rj = CisAjsCktAltDCIdx[idx][2];
      t  = CisAjsCktAltDCIdx[idx][3];
//
      rk = CisAjsCktAltDCIdx[idx][4];
      u  = CisAjsCktAltDCIdx[idx][5];
      rl = CisAjsCktAltDCIdx[idx][6];
      v  = CisAjsCktAltDCIdx[idx][7];

      if(s==t && u==v){
        tmp = GreenFunc2_fsz(ri,rj,rk,rl,s,u,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,
                       myProjCntNew,myBuffer);
      }else{
        tmp = GreenFunc2_fsz2(ri,rj,rk,rl,s,t,u,v,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myEleSpn,
                       myProjCntNew,myBuffer);
      }
      PhysCisAjsCktAltDC[idx] += w*tmp;
    }
    
    #pragma omp master
    {StopTimer(51);StartTimer(52);}

    #pragma omp for private(idx) nowait
    for(idx=0;idx<NCisAjs;idx++) {
      PhysCisAjs[idx] += w*LocalCisAjs[idx];
    }
    
    #pragma omp master
    {StopTimer(52);StartTimer(53);}

    #pragma omp for private(idx,idx0,idx1) nowait
    for(idx=0;idx<NCisAjsCktAlt;idx++) {
      idx0 = CisAjsCktAltIdx[idx][0];
      idx1 = CisAjsCktAltIdx[idx][1];
      PhysCisAjsCktAlt[idx] += w*LocalCisAjs[idx0]*conj(LocalCisAjs[idx1]);// TBC conj ok?
    }

    #pragma omp master
    {StopTimer(53);}
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  return;
}
#endif
