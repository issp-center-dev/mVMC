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
 * calculate Hamiltonian
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/
#include "calham.h"

#pragma once

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

  const int nHamiltonianTwo = NPairHopping + NExchangeCoupling * 2 + NInterAll;
  int *lazy_info = malloc(    sizeof(int) * nHamiltonianTwo * 2);
  int *lazy_rsi  = malloc(2 * sizeof(int) * nHamiltonianTwo);
  int *lazy_msj  = malloc(2 * sizeof(int) * nHamiltonianTwo);
  double complex *lazy_ip  = malloc(sizeof(double complex) * nHamiltonianTwo);
  double complex *lazy_pfa = malloc(sizeof(double complex) * nHamiltonianTwo * NQPFull);
  memset(lazy_info,                    0, sizeof(int) * nHamiltonianTwo);
  memset(lazy_info + nHamiltonianTwo, -1, sizeof(int) * nHamiltonianTwo);

  for (int mi=0; mi<Ne;  mi++) EleSpn[mi] = 0;
  for (int mi=Ne;mi<Ne*2;mi++) EleSpn[mi] = 1;

  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadComplex(NQPFull+2*Nsize);
  /* GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize */

  // shared variables: eleCfg, eleProjCnt, eleIdx, eleNum
  // implicitly shared variables: [ip].
  // OpenMP =3.X does NOT allow specifying implicitly shared variable [ip]
  // as shared(), while in OpenMP >4.0 [ip] MUST be specified.
  // Hence we have to fall back to default(shared) here.
#pragma omp parallel default(shared)                                      \
  private(myEleIdx,myEleNum,myProjCntNew,myBuffer,myEnergy, idx, ri, rj, rk, rl, s, t) \
  firstprivate(Nsize, Nsite2, NProj, NQPFull, NCoulombIntra, CoulombIntra, ParaCoulombIntra, \
               NCoulombInter, CoulombInter, ParaCoulombInter, NHundCoupling, HundCoupling, ParaHundCoupling, \
               NTransfer, Transfer, ParaTransfer, NPairHopping, PairHopping, ParaPairHopping, \
               NExchangeCoupling, ExchangeCoupling, ParaExchangeCoupling, NInterAll, InterAll, ParaInterAll, n0, n1) \
  shared(eleCfg, eleProjCnt, eleIdx, eleNum, lazy_info) reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myBuffer = GetWorkSpaceThreadComplex(NQPFull+2*Nsize);

    void *pfOrbital[NQPFull];
    void *pfUpdator[NQPFull];
    void *pfMat[NQPFull];
    void *pfMap[NQPFull];

    // Attaching thread-private objects to thread-shared InvM.
    // These objects no long need mutating states in this use. Just functor-like stuff.
    updated_tdi_v_seq_init_precomp_z(NQPFull, Nsite, Nsite2, Nsize,
                                     SlaterElm, Nsite2*Nsite2,
                                     InvM, Nsize*Nsize,
                                     eleIdx, EleSpn,
                                     2 /* GF @ measure: 2 at max. */, PfM,
                                     pfUpdator, pfOrbital, pfMat, pfMap);

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

    /****************/
    /* 2-body terms */
    /****************/
    int noffset_lazy = 0, noffset_rsij = 0;

    /* Pair Hopping */
    #pragma omp for private(idx,ri,rj) schedule(dynamic) nowait
    for(idx=0;idx<NPairHopping;idx++) {
      ri = PairHopping[idx][0];
      rj = PairHopping[idx][1];

      int *lazy_info_loc = lazy_info + noffset_lazy + idx;
      int *lazy_rsi_loc  = lazy_rsi  + noffset_rsij + idx * 2;
      int *lazy_msj_loc  = lazy_msj  + noffset_rsij + idx * 2;
      double complex *lazy_ip_loc = lazy_ip + noffset_lazy + idx;

      *lazy_ip_loc = ParaPairHopping[idx] *
        GreenFunc2_(ri,rj,ri,rj,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer,
                    lazy_info_loc, lazy_rsi_loc, lazy_msj_loc);
      if ( !*lazy_info_loc ) myEnergy += *lazy_ip_loc;
    }
    noffset_lazy += NPairHopping;
    noffset_rsij += NPairHopping * 2;

    /* Exchange Coupling */
    #pragma omp for private(idx,ri,rj,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NExchangeCoupling;idx++) {
      ri = ExchangeCoupling[idx][0];
      rj = ExchangeCoupling[idx][1];
    
      int *lazy_info_even = lazy_info + noffset_lazy +  idx * 2;
      int *lazy_info_odd  = lazy_info + noffset_lazy +  idx * 2 + 1;
      int *lazy_rsi_even  = lazy_rsi  + noffset_rsij + (idx * 2    ) * 2;
      int *lazy_rsi_odd   = lazy_rsi  + noffset_rsij + (idx * 2 + 1) * 2;
      int *lazy_msj_even  = lazy_msj  + noffset_rsij + (idx * 2    ) * 2;
      int *lazy_msj_odd   = lazy_msj  + noffset_rsij + (idx * 2 + 1) * 2;
      double complex *lazy_ip_even = lazy_ip + noffset_lazy + idx * 2;
      double complex *lazy_ip_odd  = lazy_ip + noffset_lazy + idx * 2 + 1;

      *lazy_ip_even = ParaExchangeCoupling[idx] * GreenFunc2_(ri,rj,rj,ri,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer, lazy_info_even, lazy_rsi_even, lazy_msj_even);
      *lazy_ip_odd  = ParaExchangeCoupling[idx] * GreenFunc2_(ri,rj,rj,ri,1,0,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer, lazy_info_odd, lazy_rsi_odd, lazy_msj_odd);
      if ( !*lazy_info_even ) myEnergy += *lazy_ip_even;
      if ( !*lazy_info_odd  ) myEnergy += *lazy_ip_odd;
    }
    noffset_lazy += 2 * NExchangeCoupling;
    noffset_rsij += 2 * NExchangeCoupling * 2;

    /* Inter All */
    #pragma omp for private(idx,ri,rj,s,rk,rl,t) schedule(dynamic) nowait
    for(idx=0;idx<NInterAll;idx++) {
      ri = InterAll[idx][0];
      rj = InterAll[idx][2];
      s  = InterAll[idx][3];
      rk = InterAll[idx][4];
      rl = InterAll[idx][6];
      t  = InterAll[idx][7];

      int *lazy_info_loc = lazy_info + noffset_lazy + idx;
      int *lazy_rsi_loc  = lazy_rsi  + noffset_rsij + idx * 2;
      int *lazy_msj_loc  = lazy_msj  + noffset_rsij + idx * 2;
      double complex *lazy_ip_loc = lazy_ip + noffset_lazy + idx;

      *lazy_ip_loc = ParaInterAll[idx] *
        GreenFunc2_(ri,rj,rk,rl,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer,
                    lazy_info_loc, lazy_rsi_loc, lazy_msj_loc);
      if ( !*lazy_info_loc ) myEnergy += *lazy_ip_loc;
    }
    noffset_lazy += NInterAll;
    noffset_rsij += NInterAll * 2;

    /* Batch-compute 2-body Green's functions. */
    #pragma omp barrier
#if 0
    int num_qp_var0 = 0;
    // Pack lazy info.
    for (idx=0; idx<nHamiltonianTwo; ++idx)
      if (lazy_info[idx]) {
        if (omp_get_thread_num() == 0) {
          lazy_rsi[num_qp_var0 * 2    ] = lazy_rsi[idx * 2]; \
          lazy_rsi[num_qp_var0 * 2 + 1] = lazy_rsi[idx * 2 + 1]; \
          lazy_msj[num_qp_var0 * 2    ] = lazy_msj[idx * 2]; \
          lazy_msj[num_qp_var0 * 2 + 1] = lazy_msj[idx * 2 + 1]; \
          lazy_info[nHamiltonianTwo + num_qp_var0] = idx; \
        }
        num_qp_var0++;
      }
    #pragma omp barrier
    updated_tdi_v_omp_var0_proc_batch_greentwo_z(NQPFull, num_qp_var0,
                                                 NULL, lazy_info + nHamiltonianTwo,
                                                 lazy_rsi, lazy_msj,
                                                 lazy_pfa,
                                                 pfUpdator, pfOrbital, pfMat, pfMap);
#else
    updated_tdi_v_omp_var1_proc_batch_greentwo_z(NQPFull, nHamiltonianTwo,
                                                 lazy_info, lazy_info + nHamiltonianTwo,
                                                 lazy_rsi, lazy_msj,
                                                 lazy_pfa,
                                                 pfUpdator, pfOrbital, pfMat, pfMap);
#endif
    #pragma omp barrier
    #pragma omp for private(idx,ri,rj,tmp) schedule(dynamic) nowait
    for(idx=0;idx<nHamiltonianTwo;idx++)
      if ( lazy_info[idx] )
        myEnergy += conj(CalculateIP_fcmp(lazy_pfa + idx * NQPFull, 0, NQPFull, MPI_COMM_SELF)) * lazy_ip[idx];

    #pragma omp master
    {StopTimer(72);}

    e += myEnergy;

    // Freeing within OMP region is thread-private.
    updated_tdi_v_free_z(NQPFull, pfUpdator, pfOrbital, pfMat, pfMap);
  }
  free(lazy_info);
  free(lazy_rsi);
  free(lazy_msj);
  free(lazy_ip);
  free(lazy_pfa);

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

double complex CalculateHamiltonianBF_fcmp(const double complex ip, int *eleIdx, const int *eleCfg,
                              int *eleNum, const int *eleProjCnt, const int *eleProjBFCnt)
{
  const int *n0 = eleNum;
  const int *n1 = eleNum + Nsite;
  double complex e=0.0, tmp;
  int idx;
  int ri,rj,s,rk,rl,t;
  int *myEleIdx, *myEleNum, *myEleCfg, *myProjCntNew, *myProjBFCntNew;
  //double sltTmp[NThread*NQPFull*Nsite2*Nsite2];
  double complex *mySltBFTmp;
  double complex *myBuffer;
  double complex myEnergy;

  RequestWorkSpaceThreadInt(Nsize+2*Nsite2+NProj+16*Nsite*Nrange);
  RequestWorkSpaceThreadComplex(NQPFull+2*Nsize+NQPFull*Nsite2*Nsite2);
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
    myBuffer   = GetWorkSpaceThreadComplex(NQPFull+2*Nsize);
    mySltBFTmp = GetWorkSpaceThreadComplex(NQPFull*Nsite2*Nsite2);

#pragma loop noalias
    for(idx=0;idx<Nsize;idx++) myEleIdx[idx] = eleIdx[idx];
#pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleNum[idx] = eleNum[idx];
#pragma loop noalias
    for(idx=0;idx<Nsite2;idx++) myEleCfg[idx] = eleCfg[idx];

    StoreSlaterElmBF_fcmp(mySltBFTmp);
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
                  *GreenFunc1BF(ri,rj,s,ip,mySltBFTmp,myEleIdx,myEleCfg,myEleNum,eleProjCnt,myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer);
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
