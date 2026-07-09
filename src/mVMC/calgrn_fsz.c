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

double complex GreenFunc1BF_fsz(const int ri, const int rj, const int s, const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer);

double complex GreenFunc2BF_fsz(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double complex ip,
                  int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, const int *eleProjBFCnt, int *projBFCntNew,
                  double complex *buffer);


void CalculateGreenFunc_fsz(const double w, const double complex ip, int *eleIdx, int *eleCfg,
                         int *eleNum, int *eleSpn,int *eleProjCnt) {

  int idx,idx0,idx1;
  int k, offset, nbody;
  int ri,rj,s,rk,rl,t,u,v;
  double complex tmp;
  int *myEleIdx, *myEleNum, *myProjCntNew,*myEleSpn, *myRsi, *myRsj;
  double complex *myBuffer;
  double *myRWork;
  double two_pi = 2.0*acos(-1.0);
  int site_idx, spin_idx, rsi_idx, rsi;
  int x, y, z;
  double Wx, Wy, Wz;
  double weight;
  int maxNBodyG = (NBodyGMaxN > 2) ? NBodyGMaxN : 2;
  int needNBodyRWork = (NNBodyG > 0 && NBodyGMaxN > 2);

  RequestWorkSpaceThreadInt(Nsize+Nsize+Nsite2+NProj+2*maxNBodyG);
  RequestWorkSpaceThreadComplex(NQPFull+maxNBodyG*Nsize);
  if (needNBodyRWork) RequestWorkSpaceThreadDouble(LapackLWork);
  /* GreenFunc1: NQPFull, GreenFunc2/NBodyG: NQPFull+maxNBodyG*Nsize */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myEleSpn,myProjCntNew,myRsi,myRsj,myBuffer,myRWork,idx)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleSpn = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myRsi = GetWorkSpaceThreadInt(maxNBodyG);
    myRsj = GetWorkSpaceThreadInt(maxNBodyG);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull+maxNBodyG*Nsize);
    myRWork = needNBodyRWork ? GetWorkSpaceThreadDouble(LapackLWork) : NULL;

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
    {StopTimer(53);StartTimer(54);}

    for(idx=0;idx<NTwist;idx++) {
      #pragma omp single
      weight = 0.0;
      #pragma omp for private(site_idx, ri, s, rsi, x,y,z, Wx, Wy,Wz,tmp) reduction(+:weight)
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

        weight += tmp*(double)myEleNum[rsi];
      }
      #pragma omp single
      PhysTwist[idx] += w*cexp(I*two_pi*weight);
    }

    #pragma omp for private(idx,k,offset,nbody,tmp) schedule(dynamic)
    for(idx=0;idx<NNBodyG;idx++) {
      nbody = NBodyGN[idx];
      offset = NBodyGOffset[idx];
      for(k=0;k<nbody;k++) {
        myRsi[k] = NBodyGIdx[offset+k][0] + NBodyGIdx[offset+k][1]*Nsite;
        myRsj[k] = NBodyGIdx[offset+k][2] + NBodyGIdx[offset+k][3]*Nsite;
      }
      tmp = GreenFuncN_fsz(nbody,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,
                           eleProjCnt,myEleSpn,myBuffer,myProjCntNew,myRWork);
      PhysNBodyG[idx] += w*tmp;
    }

    #pragma omp master
    {StopTimer(54);}
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  if (needNBodyRWork) ReleaseWorkSpaceThreadDouble();
  return;
}

void CalculateGreenFuncBF_fsz(const double w, const double complex ip, int *eleIdx, int *eleCfg,
                         int *eleNum, int *eleSpn, int *eleProjCnt,
                         const int *eleProjBFCnt) {
  int idx,idx0,idx1;
  int ri,rj,s,t,rk,rl,u,v;
  double complex tmp;
  int *myEleIdx, *myEleCfg, *myEleNum, *myEleSpn;
  int *myProjCntNew, *myProjBFCntNew;
  double complex *myBuffer;
  const int bfSlaterSize = NQPFull*Nsite2*Nsite2;

  if(NTwist > 0 || NNBodyG > 0) {
    fprintf(stderr, "Error: CalculateGreenFuncBF_fsz does not support Twist or NBodyG yet.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  RequestWorkSpaceInt(Nsize+Nsite2+Nsite2+Nsize+NProj+16*Nsite*Nrange);
  RequestWorkSpaceComplex(NQPFull + bfSlaterSize);

  myEleIdx = GetWorkSpaceInt(Nsize);
  myEleCfg = GetWorkSpaceInt(Nsite2);
  myEleNum = GetWorkSpaceInt(Nsite2);
  myEleSpn = GetWorkSpaceInt(Nsize);
  myProjCntNew = GetWorkSpaceInt(NProj);
  myProjBFCntNew = GetWorkSpaceInt(16*Nsite*Nrange);
  myBuffer = GetWorkSpaceComplex(NQPFull + bfSlaterSize);

  for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
  for(idx=0;idx<Nsite2;idx++) myEleCfg[idx] = eleCfg[idx];
  for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
  for(idx=0;idx<Nsize;idx++) myEleSpn[idx] = eleSpn[idx];

  StartTimer(50);
  for(idx=0;idx<NCisAjs;idx++) {
    ri = CisAjsIdx[idx][0];
    s  = CisAjsIdx[idx][1];
    rj = CisAjsIdx[idx][2];
    t  = CisAjsIdx[idx][3];
    if(s != t) {
      fprintf(stderr, "Error: BackFlow FSZ does not support spin-changing OneBodyG yet.\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    tmp = GreenFunc1BF_fsz(ri,rj,s,ip,myEleIdx,myEleCfg,myEleNum,eleProjCnt,myEleSpn,
                           myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer);
    LocalCisAjs[idx] = tmp;
  }
  StopTimer(50);

  StartTimer(51);
  for(idx=0;idx<NCisAjsCktAltDC;idx++) {
    ri = CisAjsCktAltDCIdx[idx][0];
    s  = CisAjsCktAltDCIdx[idx][1];
    rj = CisAjsCktAltDCIdx[idx][2];
    t  = CisAjsCktAltDCIdx[idx][3];
    rk = CisAjsCktAltDCIdx[idx][4];
    u  = CisAjsCktAltDCIdx[idx][5];
    rl = CisAjsCktAltDCIdx[idx][6];
    v  = CisAjsCktAltDCIdx[idx][7];

    if(s != t || u != v) {
      fprintf(stderr, "Error: BackFlow FSZ does not support spin-changing TwoBodyG yet.\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    tmp = GreenFunc2BF_fsz(ri,rj,rk,rl,s,u,ip,myEleIdx,myEleCfg,myEleNum,eleProjCnt,myEleSpn,
                           myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer);
    PhysCisAjsCktAltDC[idx] += w*tmp;
  }
  StopTimer(51);

  StartTimer(52);
  for(idx=0;idx<NCisAjs;idx++) {
    PhysCisAjs[idx] += w*LocalCisAjs[idx];
  }
  StopTimer(52);

  StartTimer(53);
  for(idx=0;idx<NCisAjsCktAlt;idx++) {
    idx0 = CisAjsCktAltIdx[idx][0];
    idx1 = CisAjsCktAltIdx[idx][1];
    PhysCisAjsCktAlt[idx] += w*LocalCisAjs[idx0]*conj(LocalCisAjs[idx1]);
  }
  StopTimer(53);

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceComplex();
  return;
}
#endif
