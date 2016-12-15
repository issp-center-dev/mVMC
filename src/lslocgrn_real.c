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
 * local Q and Green Functions with Lanczos step
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#ifndef _SRC_LSLOCGRN_REAL
#define _SRC_LSLOCGRN_REAL
#include "lslocgrn_real.h"

/* Calculate <psi|QQ|x>/<psi|x> */
void LSLocalQ_real(const double h1, const double ip, int *eleIdx, int *eleCfg, int *eleNum, int *eleProjCnt)
{
  double complex e0,h2;

  e0 = CalculateHamiltonian0_real(eleNum); /* V */

  h2 = h1*e0; /* HV = (V+K+W)V */
  h2 += calculateHK(h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
  h2 += calculateHW(h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);

  /* calculate local Q (IQ) */
  LSLQ[0] = 1.0; /* I */
  LSLQ[1] = h1;  /* H = V+K+W */

  /* calculate local Q (KQ) */
  LSLQ[2] = h1;  /* H */
  LSLQ[3] = h2;  /* H*H */

  return;
}

///
/// \param h1
/// \param ip
/// \param eleIdx
/// \param eleCfg
/// \param eleNum
/// \param eleProjCnt
/// \return val
/// \version 1.0
double calculateHK_real(const double h1, const double ip, int *eleIdx, int *eleCfg,
                           int *eleNum, int *eleProjCnt) {
  int idx,ri,rj,s;
  double val=0.0;

  for(idx=0;idx<NTransfer;idx++) {
    ri = Transfer[idx][0];
    rj = Transfer[idx][2];
    s  = Transfer[idx][3];

    val -= creal(ParaTransfer[idx]) * calHCA_real(ri,rj,s,h1,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
    /* Caution: negative sign */
  }

  return val;
}

/// calculate <psi| H C_is A_js |x>/<psi|x>
/// \param ri
/// \param rj
/// \param s
/// \param h1
/// \param ip
/// \param eleIdx
/// \param eleCfg
/// \param eleNum
/// \param eleProjCnt
/// \return val
/// \version 1.0
double calHCA_real(const int ri, const int rj, const int s,
                      const double h1, const double ip, int *eleIdx, int *eleCfg,
                      int *eleNum, int *eleProjCnt) {
  int rsi=ri+s*Nsite;
  int rsj=rj+s*Nsite;
  double val;
  double g;

  /* check */
  if(rsi==rsj) {
    if(eleNum[rsi]==1) return h1;
    else return 0.0;
  } else {
    if(eleNum[rsj]==0) return 0.0;
    if(eleNum[rsi]==1) return 0.0;
  }

  g = checkGF1_real(ri,rj,s,ip,eleIdx,eleCfg,eleNum);
  if(fabs(g)>1.0e-12) {
    val = calHCA1_real(ri,rj,s,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
  } else {
    val = calHCA2_real(ri,rj,s,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
  }

  return val;
}

///
/// \param ri
/// \param rj
/// \param s
/// \param ip
/// \param eleIdx
/// \param eleCfg
/// \param eleNum
/// \return z/ip
/// \version 1.0
double checkGF1_real(const int ri, const int rj, const int s, const double ip,
                int *eleIdx, const int *eleCfg, int *eleNum) {
  double z;
  int mj,msj,rsi,rsj;
  double pfMNew[NQPFull];

  mj = eleCfg[rj+s*Nsite];
  msj = mj + s*Ne;
  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;

  /* hopping */
  eleIdx[msj] = ri;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;

  /* calculate Pfaffian */
  CalculateNewPfM_real(mj, s, pfMNew, eleIdx, 0, NQPFull);
  z = CalculateIP_real(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  return z/ip;
}

/// calculate <psi| H C_is A_js |x>/<psi|x> = <psi|x'>/<psi|x> * <psi|H|x'>/<psi|x'>
/// \param ri
/// \param rj
/// \param s
/// \param ip
/// \param eleIdx
/// \param eleCfg
/// \param eleNum
/// \param eleProjCnt
/// \return
/// \version 1.0
double calHCA1_real(const int ri, const int rj, const int s,
               const double ip, int *eleIdx, int *eleCfg,
               int *eleNum, int *eleProjCnt) {
    double *oldInvM; /* [NQPFull*Nsize*Nsize;] */
    double *oldPfM;  /* [NQPFull] */
    int *projCntNew;

    int rsi=ri+s*Nsite;
    int rsj=rj+s*Nsite;
    int mj;
    double ipNew,z,e;

    RequestWorkSpaceInt(NProj);
    RequestWorkSpaceDouble(NQPFull*(Nsize*Nsize+1));

    projCntNew = GetWorkSpaceInt(NProj);
    oldInvM = GetWorkSpaceDouble(NQPFull*Nsize*Nsize);
    oldPfM  = GetWorkSpaceDouble(NQPFull);

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
    z = ProjRatio(projCntNew,eleProjCnt);

    UpdateMAll(mj,s,eleIdx,0,NQPFull);
    ipNew = CalculateIP(PfM,0,NQPFull,MPI_COMM_SELF);

    e = CalculateHamiltonian(ipNew,eleIdx,eleCfg,eleNum,projCntNew);

    /* revert hopping */
    eleIdx[mj+s*Ne] = rj;
    eleCfg[rsj] = mj;
    eleCfg[rsi] = -1;
    eleNum[rsj] = 1;
    eleNum[rsi] = 0;

    /* restore InvM and PfM */
    copyMAll(oldInvM,oldPfM,InvM,PfM);

    ReleaseWorkSpaceInt();
    ReleaseWorkSpaceDouble();
    return e*z*ipNew/ip;
}

/* calculate <psi| H C_is A_js |x>/<psi|x> for <psi|CA|x>/<psi|x>=0 */
/* Assuming ri!=rj, eleNum[rsi]=1, eleNum[rsj]=0 */
double complex calHCA2_real(const int ri, const int rj, const int s,
                       const double ip, int *eleIdx, int *eleCfg,
                       int *eleNum, int *eleProjCnt) {
    const int nsize=Nsize;
    const int nsite2=Nsite2;

    int idx;
    int rk,rl,sk;
    double complex val=0.0;
    int rsi = ri+s*Nsite;
    int rsj = rj+s*Nsite;

    double g;

    double *buffer;
    int *bufferInt;

    int *myEleIdx, *myEleNum, *myBufferInt, *myRsi, *myRsj;
    double *myBuffer;
    double complex myValue=0;
    double complex v=0.0;

    RequestWorkSpaceInt(NProj);      /* for GreenFunc1 */
    RequestWorkSpaceDouble(NQPFull); /* for GreenFunc1 */
    RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj+6);
    RequestWorkSpaceThreadDouble(NQPFull+3*Nsize);

    bufferInt = GetWorkSpaceInt(NProj);
    buffer = GetWorkSpaceDouble(NQPFull);

    /* H0 term */
    /* <psi|H0 CA|x>/<psi|x> = H0(x') <psi|CA|x>/<psi|x> */
    g = GreenFunc1(ri,rj,s,ip,eleIdx,eleCfg,eleNum,eleProjCnt,bufferInt,buffer);

    /* hopping */
    eleNum[rsi] = 1;
    eleNum[rsj] = 0;

    val = g*CalculateHamiltonian0(eleNum);

    /* revert hopping */
    eleNum[rsi] = 0;
    eleNum[rsj] = 1;

    /* end of H0 term */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myBufferInt,myBuffer,myValue,myRsi,myRsj)  \
  reduction(+:v)
    {
        myEleIdx = GetWorkSpaceThreadInt(Nsize);
        myEleNum = GetWorkSpaceThreadInt(Nsite2);
        myBufferInt = GetWorkSpaceThreadInt(NProj);
        myRsi = GetWorkSpaceThreadInt(3);
        myRsj = GetWorkSpaceThreadInt(3);
        myBuffer = GetWorkSpaceThreadDouble(NQPFull+3*Nsize);

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
                       * GreenFunc2(rk,rl,ri,rj,sk,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBufferInt,myBuffer);
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
                       * GreenFuncN(3,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBuffer,myBufferInt);
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
                       * GreenFuncN(3,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBuffer,myBufferInt);

            myRsi[0] = rk+Nsite; /* s=1 */
            myRsj[0] = rl+Nsite; /* s=1 */
            myRsi[1] = rl; /* s=0 */
            myRsj[1] = rk; /* s=0 */
            myRsi[2] = rsi;
            myRsj[2] = rsj;
            myValue += ParaExchangeCoupling[idx]
                       * GreenFuncN(3,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBuffer,myBufferInt);
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
                       * GreenFuncN(3,myRsi,myRsj,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myBuffer,myBufferInt);
        }

        v += myValue;
    }
    val += v;

    ReleaseWorkSpaceInt();
    ReleaseWorkSpaceDouble();
    ReleaseWorkSpaceThreadInt();
    ReleaseWorkSpaceThreadDouble();

    return val;
}

#endif