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
 * local Q and Green Functions with Lanczos step
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "lslocgrn.h"

#ifndef _SRC_LSLOCGRN_CMP
#define _SRC_LSLOCGRN_CMP

#include "locgrn.h"
#include "workspace.h"
#include "vmccal.h"
#include "calham.h"
#include "qp.h"
#include "pfupdate.h"
#include "pfupdate_two_fcmp.h"
#include "projection.h"

double complex calculateHK(const double complex h1, const double complex ip, int *eleIdx, int *eleCfg,
                   int *eleNum, int *eleProjCnt, double complex *rbmCnt);
double complex calculateHW(const double complex h1, const double complex ip, int *eleIdx, int *eleCfg,
                   int *eleNum, int *eleProjCnt, double complex *rbmCnt);

double complex calHCA(const int ri, const int rj, const int s,
              const double complex h1, const double complex ip, int *eleIdx, int *eleCfg,
              int *eleNum, int *eleProjCnt, double complex *rbmCnt);
double complex calHCACA(const int ri, const int rj, const int rk, const int rl,
                const int si,const int sk,
                const double complex h1, const double complex ip, int *eleIdx, int *eleCfg,
                int *eleNum, int *eleProjCnt, double complex *rbmCnt);

double complex checkGF1(const int ri, const int rj, const int s, const double complex ip,
                int *eleIdx, const int *eleCfg, int *eleNum);
double complex calHCA1(const int ri, const int rj, const int s,
               const double complex ip, int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt, double complex *rbmCnt);
double complex calHCA2(const int ri, const int rj, const int s,
               const double complex ip, int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt, double complex *rbmCnt);

double complex checkGF2(const int ri, const int rj, const int rk, const int rl,
                const int s, const int t, const double complex ip,
                int *eleIdx, const int *eleCfg, int *eleNum);
double complex calHCACA1(const int ri, const int rj, const int rk, const int rl,
                 const int si,const int sk,
                 const double complex ip, int *eleIdx, int *eleCfg,
                 int *eleNum, int *eleProjCnt, double complex *rbmCnt);
double complex calHCACA2(const int ri, const int rj, const int rk, const int rl,
                 const int si,const int sk,
                 const double complex ip, int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt, double complex *rbmCnt);

void copyMAll(double complex *invM_from, double complex *pfM_from, double complex *invM_to, double complex *pfM_to);

/* Calculate <psi|QQ|x>/<psi|x> */
void LSLocalQ(const double complex h1, const double complex ip, int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt, double complex *rbmCnt, double complex *_LSLQ)
{
  double complex e0,h2;

  e0 = CalculateHamiltonian0(eleNum); /* V */

  h2 = h1*e0; /* HV = (V+K+W)V */
  h2 += calculateHK(h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
  h2 += calculateHW(h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);

  /* calculate local Q (IQ) */
  _LSLQ[0] = 1.0; /* I */
  _LSLQ[1] = h1;  /* H = V+K+W */

  /* calculate local Q (KQ) */
  _LSLQ[2] = h1;  /* H */
  _LSLQ[3] = h2;  /* H*H */

  return;
}

/* Calculate <psi|QCisAjs|x>/<psi|x> */
void LSLocalCisAjs(const double complex h1, const double complex ip, int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt, double complex *rbmCnt) {
  const int nCisAjs=NCisAjs;
  double complex*lsLCisAjs = LSLCisAjs;
  double complex*localCisAjs = LocalCisAjs;
  int ri,rj,s;
  int idx;

  /* copy local ICisAjs */
  #pragma loop noalias
  for(idx=0;idx<nCisAjs;idx++){
    lsLCisAjs[idx] = localCisAjs[idx];
  }

  for(idx=0;idx<nCisAjs;idx++){
    ri = CisAjsIdx[idx][0];
    rj = CisAjsIdx[idx][2];
    s  = CisAjsIdx[idx][3];

    /* calculate local HCisAjs */
    LSLCisAjs[idx+nCisAjs] = calHCA(ri,rj,s,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
  }
  return;
}

double complex calculateHK(const double complex h1, const double complex ip, int *eleIdx, int *eleCfg,
                   int *eleNum, int *eleProjCnt, double complex *rbmCnt) {
  int idx,ri,rj,s;
  double complex val=0.0;

  for(idx=0;idx<NTransfer;idx++) {
    ri = Transfer[idx][0];
    rj = Transfer[idx][2];
    s  = Transfer[idx][3];
  
    val -= ParaTransfer[idx] * calHCA(ri,rj,s,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
    /* Caution: negative sign */
  }

  return val;
}

double complex calculateHW(const double complex h1, const double complex ip, int *eleIdx, int *eleCfg,
                   int *eleNum, int *eleProjCnt, double complex *rbmCnt) {
  int idx,ri,rj,s,rk,rl,t;
  double complex val=0.0,tmp;

  /* Pair Hopping */
  for(idx=0;idx<NPairHopping;idx++) {
    ri = PairHopping[idx][0];
    rj = PairHopping[idx][1];
    
    val += ParaPairHopping[idx]
      * calHCACA(ri,rj,ri,rj,0,1,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
  }

  /* Exchange Coupling */
  for(idx=0;idx<NExchangeCoupling;idx++) {
    ri = ExchangeCoupling[idx][0];
    rj = ExchangeCoupling[idx][1];
   
    tmp =  calHCACA(ri,rj,rj,ri,0,1,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
    tmp += calHCACA(ri,rj,rj,ri,1,0,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
    val += ParaExchangeCoupling[idx] * tmp;
  }

  /* Inter All */
  for(idx=0;idx<NInterAll;idx++) {
    ri = InterAll[idx][0];
    rj = InterAll[idx][2];
    s  = InterAll[idx][3];
    rk = InterAll[idx][4];
    rl = InterAll[idx][6];
    t  = InterAll[idx][7];
      
    val += ParaInterAll[idx]
      * calHCACA(ri,rj,rk,rl,s,t,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
  }

  return val;
}

/* calculate <psi| H C_is A_js |x>/<psi|x> */
double complex calHCA(const int ri, const int rj, const int s,
              const double complex h1, const double complex ip, int *eleIdx, int *eleCfg,
              int *eleNum, int *eleProjCnt, double complex *rbmCnt) {
  int rsi=ri+s*Nsite;
  int rsj=rj+s*Nsite;
  double complex val;
  double complex g;

  /* check */
  if(rsi==rsj) {
    if(eleNum[rsi]==1) return h1;
    else return 0.0;
  } else {
    if(eleNum[rsj]==0) return 0.0;
    if(eleNum[rsi]==1) return 0.0;
  }

  g = checkGF1(ri,rj,s,ip,eleIdx,eleCfg,eleNum);
  if(cabs(g)>1.0e-12) {
    val = calHCA1(ri,rj,s,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
  } else {
    val = calHCA2(ri,rj,s,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
  }

  return val;
}

double complex checkGF1(const int ri, const int rj, const int s, const double complex ip,
                int *eleIdx, const int *eleCfg, int *eleNum) {
  double complex z;
  int mj,msj,rsi,rsj;
  double complex pfMNew[NQPFull];

  mj = eleCfg[rj+s*Nsite];
  msj = mj + s*Ne;
  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;

  /* hopping */
  eleIdx[msj] = ri;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;

  /* calculate Pfaffian */
    CalculateNewPfM(mj, s, pfMNew, eleIdx, 0, NQPFull);
    z = CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  return z/ip;
}

/* calculate <psi| H C_is A_js |x>/<psi|x> = <psi|x'>/<psi|x> * <psi|H|x'>/<psi|x'> */
double complex calHCA1(const int ri, const int rj, const int s,
               const double complex ip, int *eleIdx, int *eleCfg,
               int *eleNum, int *eleProjCnt, double complex *rbmCnt) {
  complex double *oldInvM; /* [NQPFull*Nsize*Nsize;] */
  complex double *oldPfM;  /* [NQPFull] */
  int *projCntNew;
  double complex *rbmCntNew;

  int rsi=ri+s*Nsite;
  int rsj=rj+s*Nsite;
  int mj;
  double complex ipNew,z,e;

  RequestWorkSpaceInt(NProj);
  RequestWorkSpaceComplex(NQPFull*(Nsize*Nsize+1) + FlagRBM*(NRBM_PhysLayerIdx+Nneuron));

  projCntNew = GetWorkSpaceInt(NProj);
  oldInvM = GetWorkSpaceComplex(NQPFull*Nsize*Nsize);
  oldPfM  = GetWorkSpaceComplex(NQPFull);
  if (FlagRBM) {
    rbmCntNew = GetWorkSpaceComplex(NRBM_PhysLayerIdx + Nneuron);
  }

  /* copy InvM and PfM */
  copyMAll(InvM,PfM,oldInvM,oldPfM);

  /* The mj-th electron with spin s hops to site ri */
  mj = eleCfg[rsj];
  eleIdx[mj+s*Ne] = ri;
  eleCfg[rsj] = -1;
  eleCfg[rsi] = mj;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;

  UpdateProjCnt(rj, ri, s, projCntNew, eleProjCnt, eleNum);
  if (FlagRBM) {
    UpdateRBMCnt(rj, ri, s, rbmCntNew, rbmCnt, eleNum);
  }
  z = ProjRatio(projCntNew,eleProjCnt);
  if (FlagRBM) {
    z *= RBMRatio(rbmCntNew,rbmCnt);
  }

  UpdateMAll(mj,s,eleIdx,0,NQPFull);
  ipNew = CalculateIP_fcmp(PfM,0,NQPFull,MPI_COMM_SELF);

  e = CalculateHamiltonian(ipNew,eleIdx,eleCfg,eleNum,projCntNew,rbmCntNew);

  /* revert hopping */
  eleIdx[mj+s*Ne] = rj;
  eleCfg[rsj] = mj;
  eleCfg[rsi] = -1;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  /* restore InvM and PfM */
  copyMAll(oldInvM,oldPfM,InvM,PfM);

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceComplex();
  return e*conj(z*ipNew/ip);
}

/* calculate <psi| H C_is A_js |x>/<psi|x> for <psi|CA|x>/<psi|x>=0 */
/* Assuming ri!=rj, eleNum[rsi]=1, eleNum[rsj]=0 */
double complex calHCA2(const int ri, const int rj, const int s,
               const double complex ip, int *eleIdx, int *eleCfg,
               int *eleNum, int *eleProjCnt, double complex *rbmCnt) {
  const int nsize=Nsize;
  const int nsite2=Nsite2;

  int idx;
  int rk,rl,sk;
  double complex val=0.0;
  int rsi = ri+s*Nsite;
  int rsj = rj+s*Nsite;

  double complex g;

  double complex *buffer;
  int *bufferInt;

  int *myEleIdx, *myEleNum, *myBufferInt, *myRsi, *myRsj;
  double complex *myBuffer;
  double complex myValue=0;
  double complex v=0.0;
  double complex *rbmCntNew, *myRBMCntNew;

  RequestWorkSpaceInt(NProj);      /* for GreenFunc1 */
  RequestWorkSpaceComplex(NQPFull + FlagRBM*(NRBM_PhysLayerIdx+Nneuron)); /* for GreenFunc1 */
  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj+6);
  RequestWorkSpaceThreadComplex(NQPFull+3*Nsize + FlagRBM*(NRBM_PhysLayerIdx+Nneuron));

  bufferInt = GetWorkSpaceInt(NProj);
  buffer = GetWorkSpaceComplex(NQPFull);
  if (FlagRBM) {
    rbmCntNew = GetWorkSpaceComplex(NRBM_PhysLayerIdx+Nneuron);
  }

  /* H0 term */
  /* <psi|H0 CA|x>/<psi|x> = H0(x') <psi|CA|x>/<psi|x> */
  g = GreenFunc1(ri,rj,s,ip,eleIdx,eleCfg,eleNum,eleProjCnt,bufferInt,rbmCnt,rbmCntNew,buffer);

  /* hopping */
  eleNum[rsi] = 1;
  eleNum[rsj] = 0;

  val = g*CalculateHamiltonian0(eleNum);

  /* revert hopping */
  eleNum[rsi] = 0;
  eleNum[rsj] = 1;

  /* end of H0 term */

#pragma omp parallel default(shared)\
  private(myRBMCntNew,myEleIdx,myEleNum,myBufferInt,myBuffer,myValue,myRsi,myRsj)  \
  reduction(+:v)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myBufferInt = GetWorkSpaceThreadInt(NProj);
    myRsi = GetWorkSpaceThreadInt(3);
    myRsj = GetWorkSpaceThreadInt(3);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull+3*Nsize);
    if (FlagRBM) {
      myRBMCntNew  = GetWorkSpaceThreadComplex(NRBM_PhysLayerIdx+Nneuron);
    }

    #pragma loop noalias
    for(idx=0;idx<nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<nsite2;idx++) myEleNum[idx] = eleNum[idx];
    
    myValue = 0.0;

    /* Transfer */
    #pragma omp for private(idx,rk,rl,sk) schedule(dynamic) nowait
    for(idx=0;idx<NTransfer;idx++) {
      rk = Transfer[idx][0];
      rl = Transfer[idx][2];
      sk = Transfer[idx][3];
      
      myValue -= ParaTransfer[idx]
        * GreenFunc2(rk,rl,ri,rj,sk,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBufferInt,rbmCnt,myRBMCntNew,myBuffer);
      /* Caution: negative sign */
    }

    /* Pair Hopping */
    #pragma omp for private(idx,rk,rl) schedule(dynamic) nowait
    for(idx=0;idx<NPairHopping;idx++) {
      rk = PairHopping[idx][0];
      rl = PairHopping[idx][1];
      myRsi[0] = rk; /* s=0 */
      myRsj[0] = rl; /* s=0 */
      myRsi[1] = rk+Nsite; /* s=1 */
      myRsj[1] = rl+Nsite; /* s=1 */
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      
      myValue += ParaPairHopping[idx]
        * GreenFuncN(3,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,rbmCnt,myRBMCntNew,myBuffer,myBufferInt);
    }

    /* Exchange Coupling */
    #pragma omp for private(idx,rk,rl) schedule(dynamic) nowait
    for(idx=0;idx<NExchangeCoupling;idx++) {
      rk = ExchangeCoupling[idx][0];
      rl = ExchangeCoupling[idx][1];
      myRsi[0] = rk; /* s=0 */
      myRsj[0] = rl; /* s=0 */
      myRsi[1] = rl+Nsite; /* s=1 */
      myRsj[1] = rk+Nsite; /* s=1 */
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      myValue += ParaExchangeCoupling[idx]
        * GreenFuncN(3,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,rbmCnt,myRBMCntNew,myBuffer,myBufferInt);

      myRsi[0] = rk+Nsite; /* s=1 */
      myRsj[0] = rl+Nsite; /* s=1 */
      myRsi[1] = rl; /* s=0 */
      myRsj[1] = rk; /* s=0 */
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      myValue += ParaExchangeCoupling[idx]
        * GreenFuncN(3,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,rbmCnt,myRBMCntNew,myBuffer,myBufferInt);
    }
    
    /* Inter All */
    #pragma omp for private(idx) schedule(dynamic) nowait
    for(idx=0;idx<NInterAll;idx++) {
      myRsi[0] = InterAll[idx][0] + InterAll[idx][3]*Nsite;
      myRsj[0] = InterAll[idx][2] + InterAll[idx][3]*Nsite;
      myRsi[1] = InterAll[idx][4] + InterAll[idx][7]*Nsite;
      myRsj[1] = InterAll[idx][6] + InterAll[idx][7]*Nsite;
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      myValue += ParaInterAll[idx]
        * GreenFuncN(3,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,rbmCnt,myRBMCntNew,myBuffer,myBufferInt);
    }
    
    v += myValue;
  }
  val += v;

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceComplex();
  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();

  return val;
}

double complex calHCACA(const int ri, const int rj, const int rk, const int rl,
                const int si,const int sk,
                const double complex h1, const double complex ip, int *eleIdx, int *eleCfg,
                int *eleNum, int *eleProjCnt, double complex *rbmCnt) {
  int rsi=ri+si*Nsite;
  int rsj=rj+si*Nsite;
  int rsk=rk+sk*Nsite;
  int rsl=rl+sk*Nsite;

  double complex val;
  double complex g;

  /* check */
  if(rsk==rsl) {
    if(eleNum[rsk]==1) {
      return calHCA(ri,rj,si,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
    } else return 0;
  } else if(rsj==rsk) {
    if(eleNum[rsj]==1) return 0;
    else {
      return calHCA(ri,rl,si,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
    }
  } else if(rsj==rsl) {
    return 0;
  } else if(rsi==rsj) {
    if(eleNum[rsi]==1) {
      return calHCA(rk,rl,sk,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
    } else return 0;
  } else if(rsi==rsk) {
    return 0;
  } else if(rsi==rsl) {
    if(eleNum[rsi]==1) {
      return -calHCA(rk,rj,sk,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
    } else return 0;
  } else {
    if(eleNum[rsl]==0) return 0.0;
    if(eleNum[rsk]==1) return 0.0;
    if(eleNum[rsj]==0) return 0.0;
    if(eleNum[rsi]==1) return 0.0;
  }

  g = checkGF2(ri,rj,rk,rl,si,sk,ip,eleIdx,eleCfg,eleNum);
  if(cabs(g)>1.0e-12) {
    val = calHCACA1(ri,rj,rk,rl,si,sk,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
  } else {
    val = calHCACA2(ri,rj,rk,rl,si,sk,ip,eleIdx,eleCfg,eleNum,eleProjCnt,rbmCnt);
  }

  return val;
}

double complex checkGF2(const int ri, const int rj, const int rk, const int rl,
                const int s, const int t, const double complex ip,
                int *eleIdx, const int *eleCfg, int *eleNum) {
  double complex z;
  int mj,msj,ml,mtl;
  int rsi,rsj,rtk,rtl;
  double complex *pfMNew;
  double complex *buffer;

  RequestWorkSpaceComplex(NQPFull+2*Nsize);
  pfMNew = GetWorkSpaceComplex(NQPFull);
  buffer = GetWorkSpaceComplex(2*Nsize);

  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;
  rtk = rk + t*Nsite;
  rtl = rl + t*Nsite;

  mj = eleCfg[rj+s*Nsite];
  ml = eleCfg[rl+t*Nsite];
  msj = mj + s*Ne;
  mtl = ml + t*Ne;

  /* hopping */
  eleIdx[mtl] = rk;
  eleNum[rtl] = 0;
  eleNum[rtk] = 1;

  eleIdx[msj] = ri;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;

  /* calculate Pfaffian */
  CalculateNewPfMTwo_fcmp(ml, t, mj, s, pfMNew, eleIdx, 0, NQPFull, buffer);
  z = CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[mtl] = rl;
  eleNum[rtl] = 1;
  eleNum[rtk] = 0;
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  ReleaseWorkSpaceComplex();
  return z/ip;
}

double complex calHCACA1(const int ri, const int rj, const int rk, const int rl,
                 const int si,const int sk,
                 const double complex ip, int *eleIdx, int *eleCfg,
                 int *eleNum, int *eleProjCnt, double complex *rbmCnt) {
  double complex *oldInvM; /* [NQPFull*Nsize*Nsize;] */
  double complex *oldPfM;  /* [NQPFull] */
  int *projCntNew;
  double complex *rbmCntNew;

  int rsi=ri+si*Nsite;
  int rsj=rj+si*Nsite;
  int rsk=rk+sk*Nsite;
  int rsl=rl+sk*Nsite;
  int mj,ml;
  double complex ipNew,z,e;

  RequestWorkSpaceInt(NProj);
  RequestWorkSpaceComplex(NQPFull*(Nsize*Nsize+1) + FlagRBM*(NRBM_PhysLayerIdx+Nneuron));

  projCntNew = GetWorkSpaceInt(NProj);
  oldInvM = GetWorkSpaceComplex(NQPFull*Nsize*Nsize);
  oldPfM  = GetWorkSpaceComplex(NQPFull);
  if (FlagRBM) {
    rbmCntNew = GetWorkSpaceComplex(NRBM_PhysLayerIdx + Nneuron);
  }

  /* copy InvM and PfM */
  copyMAll(InvM,PfM,oldInvM,oldPfM);

  /* The ml-th electron with spin sk hops from rl to rk */
  ml = eleCfg[rsl];
  eleIdx[ml+sk*Ne] = rk;
  eleCfg[rsl] = -1;
  eleCfg[rsk] = ml;
  eleNum[rsl] = 0;
  eleNum[rsk] = 1;
  UpdateProjCnt(rl, rk, sk, projCntNew, eleProjCnt, eleNum);
  if (FlagRBM) {
    UpdateRBMCnt(rl, rk, sk, rbmCntNew, rbmCnt, eleNum);
  }

  /* The mj-th electron with spin si hops from rj to ri */
  mj = eleCfg[rsj];
  eleIdx[mj+si*Ne] = ri;
  eleCfg[rsj] = -1;
  eleCfg[rsi] = mj;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt(rj, ri, si, projCntNew, projCntNew, eleNum);
  if (FlagRBM) {
    UpdateRBMCnt(rj, ri, si, rbmCntNew, rbmCntNew, eleNum);
  }

  z = ProjRatio(projCntNew,eleProjCnt);
  if (FlagRBM) {
    z *= RBMRatio(rbmCntNew,rbmCnt);
  }

  UpdateMAllTwo_fcmp(ml, sk, mj, si, rl, rj, eleIdx, 0, NQPFull);
  ipNew = CalculateIP_fcmp(PfM,0,NQPFull,MPI_COMM_SELF);

  e = CalculateHamiltonian(ipNew,eleIdx,eleCfg,eleNum,projCntNew,rbmCnt);

  /* revert hopping */
  eleIdx[mj+si*Ne] = rj;
  eleCfg[rsj] = mj;
  eleCfg[rsi] = -1;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  eleIdx[ml+sk*Ne] = rl;
  eleCfg[rsl] = ml;
  eleCfg[rsk] = -1;
  eleNum[rsl] = 1;
  eleNum[rsk] = 0;

  /* restore InvM and PfM */
  copyMAll(oldInvM,oldPfM,InvM,PfM);

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceComplex();
  if (FlagRBM) {
    return e*conj(z*ipNew/ip);
  }
  else {
    return e*z*ipNew/ip;
  }
}

/* calculate <psi| H C_is A_js C_kt A_lt|x>/<psi|x> for <psi|CACA|x>/<psi|x>=0 */
/* Assuming ri,rj,rk,rl are different, eleNum[rsi]=1, eleNum[rsj]=0, eleNum[rsk]=1, eleNum[rsl]=0  */
double complex calHCACA2(const int ri, const int rj, const int rk, const int rl,
                 const int si,const int sk,
                 const double complex ip, int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt, double complex *rbmCnt) {
  const int nsize=Nsize;
  const int nsite2=Nsite2;

  const int rsi = ri+si*Nsite;
  const int rsj = rj+si*Nsite;
  const int rsk = rk+sk*Nsite;
  const int rsl = rl+sk*Nsite;

  int idx,r0,r1;
  double complex val=0.0;

  double complex g;

  double complex *buffer;
  int *bufferInt;

  int *myEleIdx, *myEleNum, *myBufferInt, *myRsi, *myRsj;
  double complex *myBuffer;
  double complex myValue=0.0;
  double complex v=0.0;
  double complex *rbmCntNew, *myRBMCntNew;

  RequestWorkSpaceInt(NProj);      /* for GreenFunc2 */
  RequestWorkSpaceComplex(NQPFull+2*Nsize + FlagRBM * (NRBM_PhysLayerIdx+Nneuron)); /* for GreenFunc1 */
  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj+8);
  RequestWorkSpaceThreadComplex(NQPFull+3*Nsize + FlagRBM * (NRBM_PhysLayerIdx+Nneuron));

  bufferInt = GetWorkSpaceInt(NProj);
  buffer = GetWorkSpaceComplex(NQPFull+2*Nsize);
  if (FlagRBM) {
    rbmCntNew = GetWorkSpaceComplex(NRBM_PhysLayerIdx+Nneuron);
  }

  /* H0 term */
  /* <psi|H0 CACA|x>/<psi|x> = H0(x') <psi|CACA|x>/<psi|x> */
  g = GreenFunc2(ri,rj,rk,rl,si,sk,ip,
                 eleIdx,eleCfg,eleNum,eleProjCnt,bufferInt,rbmCnt,myRBMCntNew,buffer);

  /* hopping */
  eleNum[rsi] = 1;
  eleNum[rsj] = 0;
  eleNum[rsk] = 1;
  eleNum[rsl] = 0;

  val = g*CalculateHamiltonian0(eleNum);

  /* revert hopping */
  eleNum[rsi] = 0;
  eleNum[rsj] = 1;
  eleNum[rsk] = 0;
  eleNum[rsl] = 1;

  /* end of H0 term */

#pragma omp parallel default(shared)\
  private(myRBMCntNew,myEleIdx,myEleNum,myBufferInt,myBuffer,myValue,myRsi,myRsj)  \
  reduction(+:v)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myBufferInt = GetWorkSpaceThreadInt(NProj);
    myRsi = GetWorkSpaceThreadInt(4);
    myRsj = GetWorkSpaceThreadInt(4);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull+4*Nsize);
    if (FlagRBM) {
      myRBMCntNew  = GetWorkSpaceThreadComplex(NRBM_PhysLayerIdx+Nneuron);
    }

    #pragma loop noalias
    for(idx=0;idx<nsize;idx++) myEleIdx[idx] = eleIdx[idx];
    #pragma loop noalias
    for(idx=0;idx<nsite2;idx++) myEleNum[idx] = eleNum[idx];
    
    myValue = 0.0;

    /* Transfer */
    #pragma omp for private(idx) schedule(dynamic) nowait
    for(idx=0;idx<NTransfer;idx++) {
      myRsi[0] = Transfer[idx][0]+Transfer[idx][1]*Nsite;
      myRsj[0] = Transfer[idx][2]+Transfer[idx][3]*Nsite;
      myRsi[1] = rsi;
      myRsj[1] = rsj;
      myRsi[2] = rsk;
      myRsj[2] = rsl;
      
      myValue -= ParaTransfer[idx]
        * GreenFuncN(3,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,rbmCnt,myRBMCntNew,myBuffer,myBufferInt);
      /* Caution: negative sign */
    }

    /* Pair Hopping */
    #pragma omp for private(idx,r0,r1) schedule(dynamic) nowait
    for(idx=0;idx<NPairHopping;idx++) {
      r0 = PairHopping[idx][0];
      r1 = PairHopping[idx][1];
      myRsi[0] = r0; /* s=0 */
      myRsj[0] = r1; /* s=0 */
      myRsi[1] = r0+Nsite; /* s=1 */
      myRsj[1] = r1+Nsite; /* s=1 */
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      myRsi[3] = rsk;
      myRsj[3] = rsl;
      
      myValue += ParaPairHopping[idx]
        * GreenFuncN(4,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,rbmCnt,myRBMCntNew,myBuffer,myBufferInt);
    }

    /* Exchange Coupling */
    #pragma omp for private(idx,r0,r1) schedule(dynamic) nowait
    for(idx=0;idx<NExchangeCoupling;idx++) {
      r0 = ExchangeCoupling[idx][0];
      r1 = ExchangeCoupling[idx][1];
      myRsi[0] = r0; /* s=0 */
      myRsj[0] = r1; /* s=0 */
      myRsi[1] = r1+Nsite; /* s=1 */
      myRsj[1] = r0+Nsite; /* s=1 */
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      myRsi[3] = rsk;
      myRsj[3] = rsl;
      myValue += ParaExchangeCoupling[idx]
        * GreenFuncN(4,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,rbmCnt,myRBMCntNew,myBuffer,myBufferInt);

      myRsi[0] = r0+Nsite; /* s=1 */
      myRsj[0] = r1+Nsite; /* s=1 */
      myRsi[1] = r1; /* s=0 */
      myRsj[1] = r0; /* s=0 */
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      myRsi[3] = rsk;
      myRsj[3] = rsl;
      myValue += ParaExchangeCoupling[idx]
        * GreenFuncN(4,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,rbmCnt,myRBMCntNew,myBuffer,myBufferInt);
    }
    
    /* Inter All */
    #pragma omp for private(idx) schedule(dynamic) nowait
    for(idx=0;idx<NInterAll;idx++) {
      myRsi[0] = InterAll[idx][0] + InterAll[idx][3]*Nsite;
      myRsj[0] = InterAll[idx][2] + InterAll[idx][3]*Nsite;
      myRsi[1] = InterAll[idx][4] + InterAll[idx][7]*Nsite;
      myRsj[1] = InterAll[idx][6] + InterAll[idx][7]*Nsite;
      myRsi[2] = rsi;
      myRsj[2] = rsj;
      myRsi[3] = rsk;
      myRsj[3] = rsl;
      myValue += ParaInterAll[idx]
        * GreenFuncN(4,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,rbmCnt,myRBMCntNew,myBuffer,myBufferInt);
    }
    
    v += myValue;
  }
  val += v;

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceComplex();
  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();

  return val;
}


/* copy invM and pfM */
void copyMAll(complex double *invM_from, complex double *pfM_from, complex double *invM_to, complex double *pfM_to) {
  int i,n;

  n = NQPFull*Nsize*Nsize;
  #pragma loop noalias
  for(i=0;i<n;i++) invM_to[i] = invM_from[i];

  n = NQPFull;
  #pragma loop noalias
  for(i=0;i<n;i++) pfM_to[i] = pfM_from[i];

  return;
}

#endif
