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
#ifndef _SRC_CALHAMREAL
#define _SRC_CALHAMREAL
#include "workspace.h"
#include "complex.h"
#include "slater.h"
#include "locgrn_real.h"
#include "calham_real.h"

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

  const int nHamiltonianTwo = NPairHopping + NExchangeCoupling * 2 + NInterAll;
  int *lazy_info = malloc(    sizeof(int) * nHamiltonianTwo * 2);
  int *lazy_rsi  = malloc(2 * sizeof(int) * nHamiltonianTwo);
  int *lazy_msj  = malloc(2 * sizeof(int) * nHamiltonianTwo);
  double *lazy_ip  = malloc(sizeof(double) * nHamiltonianTwo);
  double *lazy_pfa = malloc(sizeof(double) * nHamiltonianTwo * NQPFull);
  memset(lazy_info,                    0, sizeof(int) * nHamiltonianTwo);
  memset(lazy_info + nHamiltonianTwo, -1, sizeof(int) * nHamiltonianTwo);

  for (int mi=0; mi<Ne;  mi++) EleSpn[mi] = 0;
  for (int mi=Ne;mi<Ne*2;mi++) EleSpn[mi] = 1;

  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadDouble(NQPFull+2*Nsize);
  /* GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize */

  // shared variables: eleCfg, eleProjCnt, eleIdx, eleNum
  // implicitly shared variables: [ip].
  // OpenMP =3.X does NOT allow specifying implicitly shared variable [ip]
  // as shared(), while in OpenMP >4.0 [ip] MUST be specified.
  // Hence we have to fall back to default(shared) here.
  #pragma omp parallel default(shared) \
    private(myEleIdx,myEleNum,myProjCntNew,myBuffer,myEnergy, idx, ri, rj, rk, rl, s, t) \
    firstprivate(Nsize, Nsite2, NProj, NQPFull, NCoulombIntra, CoulombIntra, ParaCoulombIntra,   \
    NCoulombInter, CoulombInter, ParaCoulombInter, NHundCoupling, HundCoupling, ParaHundCoupling,    \
    NTransfer, Transfer, ParaTransfer, NPairHopping, PairHopping, ParaPairHopping,    \
    NExchangeCoupling, ExchangeCoupling, ParaExchangeCoupling, NInterAll, InterAll, ParaInterAll, n0, n1)\
    shared(eleCfg, eleProjCnt, eleIdx, eleNum, lazy_info) reduction(+:e)
  {
    myEleIdx = GetWorkSpaceThreadInt(Nsize);
    myEleNum = GetWorkSpaceThreadInt(Nsite2);
    myProjCntNew = GetWorkSpaceThreadInt(NProj);
    myBuffer = GetWorkSpaceThreadDouble(NQPFull+2*Nsize);

    void *pfOrbital[NQPFull];
    void *pfUpdator[NQPFull];
    void *pfMat[NQPFull];
    void *pfMap[NQPFull];

    // Attaching thread-private objects to thread-shared InvM.
    // These objects no long need mutating states in this use. Just functor-like stuff.
    updated_tdi_v_seq_init_precomp_d(NQPFull, Nsite, Nsite2, Nsize,
                                     SlaterElm_real, Nsite2*Nsite2,
                                     InvM_real, Nsize*Nsize,
                                     eleIdx, EleSpn,
                                     2 /* GF @ measure: 2 at max. */, PfM_real,
                                     pfUpdator, pfOrbital, pfMat, pfMap);

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
      double *lazy_ip_loc = lazy_ip + noffset_lazy + idx;

      *lazy_ip_loc = ParaPairHopping[idx]
        * GreenFunc2_real_(ri,rj,ri,rj,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer,
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
      double *lazy_ip_even = lazy_ip + noffset_lazy + idx * 2;
      double *lazy_ip_odd  = lazy_ip + noffset_lazy + idx * 2 + 1;

      *lazy_ip_even = ParaExchangeCoupling[idx] * GreenFunc2_real_(ri,rj,rj,ri,0,1,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer, lazy_info_even, lazy_rsi_even, lazy_msj_even);
      *lazy_ip_odd  = ParaExchangeCoupling[idx] * GreenFunc2_real_(ri,rj,rj,ri,1,0,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer, lazy_info_odd, lazy_rsi_odd, lazy_msj_odd);
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
      double *lazy_ip_loc = lazy_ip + noffset_lazy + idx;

      *lazy_ip_loc = ParaInterAll[idx]
        * GreenFunc2_real_(ri,rj,rk,rl,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer,
                           lazy_info_loc, lazy_rsi_loc, lazy_msj_loc);
      if ( !*lazy_info_loc ) myEnergy += *lazy_ip_loc;
    }
    noffset_lazy += NInterAll;
    noffset_rsij += NInterAll * 2;

    /* Batch-compute 2-body Green's functions. */
    #pragma omp barrier
    // if ( Nsize <= 100 ) // Heuristics: Huge Nelec seems to cause parallelize-over-nGF spill L2.
    if ( 1 )
    {
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
      fprintf(stdout, "GF2: %d / %d\n", num_qp_var0, nHamiltonianTwo);
      #pragma omp barrier
      updated_tdi_v_omp_var0_proc_batch_greentwo_d(NQPFull, num_qp_var0,
                                                   NULL, lazy_info + nHamiltonianTwo,
                                                   lazy_rsi, lazy_msj,
                                                   lazy_pfa,
                                                   pfUpdator, pfOrbital, pfMat, pfMap);
    } else {
      updated_tdi_v_omp_var1_proc_batch_greentwo_d(NQPFull, nHamiltonianTwo,
                                                   lazy_info, lazy_info + nHamiltonianTwo,
                                                   lazy_rsi, lazy_msj,
                                                   lazy_pfa,
                                                   pfUpdator, pfOrbital, pfMat, pfMap);
    }
    #pragma omp barrier
    #pragma omp for private(idx,ri,rj,tmp) schedule(dynamic) nowait
    for(idx=0;idx<nHamiltonianTwo;idx++)
      if ( lazy_info[idx] )
        myEnergy += CalculateIP_real(lazy_pfa + idx * NQPFull, 0, NQPFull, MPI_COMM_SELF) * lazy_ip[idx];

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
        // * GreenFunc1(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer);
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
