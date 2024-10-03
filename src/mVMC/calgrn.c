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
#include "calgrn.h"
#ifndef _CALGRN_SRC
#define _CALGRN_SRC

void CalculateGreenFunc(const double w, const double complex ip, int *eleIdx, int *eleCfg,
                        int *eleNum, int *eleProjCnt) {

  int idx,idx0,idx1;
  int ri,rj,s,rk,rl,t;
  double complex tmp;
  int *myEleIdx, *myEleNum, *myProjCntNew;
  double complex *myBuffer;

  int *lazy_info = malloc(    sizeof(int) * NCisAjsCktAltDC * 2);
  int *lazy_rsi  = malloc(2 * sizeof(int) * NCisAjsCktAltDC);
  int *lazy_msj  = malloc(2 * sizeof(int) * NCisAjsCktAltDC);
  double complex *lazy_ip  = malloc(sizeof(double complex) * NCisAjsCktAltDC);
  double complex *lazy_pfa = malloc(sizeof(double complex) * NCisAjsCktAltDC * NQPFull);
  memset(lazy_info,                    0, sizeof(int) * NCisAjsCktAltDC);
  memset(lazy_info + NCisAjsCktAltDC, -1, sizeof(int) * NCisAjsCktAltDC);

  for (int mi=0; mi<Ne;  mi++) EleSpn[mi] = 0;
  for (int mi=Ne;mi<Ne*2;mi++) EleSpn[mi] = 1;

  RequestWorkSpaceThreadInt(Nsize+Nsite2+NProj);
  RequestWorkSpaceThreadComplex(NQPFull+2*Nsize);
  /* GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize */

  #pragma omp parallel default(shared)		\
  private(myEleIdx,myEleNum,myProjCntNew,myBuffer,idx)
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

    #pragma omp master
    {StartTimer(50);}

    #pragma omp for private(idx,ri,rj,s,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NCisAjs;idx++) {
      ri = CisAjsIdx[idx][0];
      rj = CisAjsIdx[idx][2];
      s  = CisAjsIdx[idx][3];
      tmp = GreenFunc1(ri,rj,s,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,
                       myProjCntNew,myBuffer);
      LocalCisAjs[idx] = tmp;
    }
    #pragma omp master
    {StopTimer(50);StartTimer(51);}
    
    #pragma omp for private(idx,ri,rj,s,rk,rl,t) schedule(dynamic)
    for(idx=0;idx<NCisAjsCktAltDC;idx++) {
      ri = CisAjsCktAltDCIdx[idx][0];
      rj = CisAjsCktAltDCIdx[idx][2];
      s  = CisAjsCktAltDCIdx[idx][1];
      rk = CisAjsCktAltDCIdx[idx][4];
      rl = CisAjsCktAltDCIdx[idx][6];
      t  = CisAjsCktAltDCIdx[idx][5];

      int *lazy_info_loc = lazy_info + idx;
      int *lazy_rsi_loc  = lazy_rsi  + idx * 2;
      int *lazy_msj_loc  = lazy_msj  + idx * 2;
      double complex *lazy_ip_loc = lazy_ip + idx;

      *lazy_ip_loc = w *
        GreenFunc2_(ri,rj,rk,rl,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,myProjCntNew,myBuffer,
          lazy_info_loc, lazy_rsi_loc, lazy_msj_loc);
      if ( !*lazy_info_loc ) {
        PhysCisAjsCktAltDC[idx] += *lazy_ip_loc;
      } /* else {
        // Reverse lookup for the idx
        *lazy_info_loc = idx + 1;
      } */
    }

    /* Batch-compute 2-body Green's functions. */
    #pragma omp barrier
    if ( Nsize <= 200 ) { // Heuristics: Huge Nelec seems to cause parallelize-over-nGF spill L2.
      int num_qp_var0 = 0;
      // Pack lazy info.
      for (idx=0; idx<NCisAjsCktAltDC; ++idx)
        if (lazy_info[idx]) {
          if (omp_get_thread_num() == 0) {
            lazy_rsi[num_qp_var0 * 2    ] = lazy_rsi[idx * 2]; \
            lazy_rsi[num_qp_var0 * 2 + 1] = lazy_rsi[idx * 2 + 1]; \
            lazy_msj[num_qp_var0 * 2    ] = lazy_msj[idx * 2]; \
            lazy_msj[num_qp_var0 * 2 + 1] = lazy_msj[idx * 2 + 1]; \
            lazy_info[NCisAjsCktAltDC + num_qp_var0] = idx; \
          }
          num_qp_var0++;
        }
      #pragma omp barrier
      updated_tdi_v_omp_var0_proc_batch_greentwo_z(NQPFull, num_qp_var0,
                                                   NULL, lazy_info + NCisAjsCktAltDC,
                                                   lazy_rsi, lazy_msj,
                                                   lazy_pfa,
                                                   pfUpdator, pfOrbital, pfMat, pfMap);
    } else {
      updated_tdi_v_omp_var1_proc_batch_greentwo_z(NQPFull, NCisAjsCktAltDC,
                                                   lazy_info, lazy_info + NCisAjsCktAltDC,
                                                   lazy_rsi, lazy_msj,
                                                   lazy_pfa,
                                                   pfUpdator, pfOrbital, pfMat, pfMap);
    }
    #pragma omp barrier
    #pragma omp for private(idx) schedule(dynamic) nowait
    for(idx=0;idx<NCisAjsCktAltDC;idx++)
      if ( lazy_info[idx] )
        PhysCisAjsCktAltDC[idx] += conj(CalculateIP_fcmp(lazy_pfa + idx * NQPFull, 0, NQPFull, MPI_COMM_SELF)) * lazy_ip[idx];

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
  free(lazy_info);
  free(lazy_rsi);
  free(lazy_msj);
  free(lazy_ip);
  free(lazy_pfa);

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  return;
}


void CalculateGreenFuncBF(const double w, const double ip, int *eleIdx, int *eleCfg,
                          int *eleNum, int *eleProjCnt, const int *eleProjBFCnt) {

  int idx,idx0,idx1;
  int ri,rj,s,rk,rl,t;
  double complex tmp;
  int *myEleIdx, *myEleNum, *myEleCfg, *myProjCntNew, *myProjBFCntNew;
  double complex* mySltBFTmp;
  double complex* myBuffer;

  RequestWorkSpaceThreadInt(Nsize+2*Nsite2+NProj+16*Nsite*Nrange);
  RequestWorkSpaceThreadComplex(NQPFull+2*Nsize+NQPFull*Nsite2*Nsite2);
  /* GreenFunc1: NQPFull, GreenFunc2: NQPFull+2*Nsize */

#pragma omp parallel default(shared)\
  private(myEleIdx,myEleNum,myEleCfg,myProjCntNew,myProjBFCntNew,myBuffer,mySltBFTmp)
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
#pragma omp master
    {StartTimer(50);}

#pragma omp for private(idx,ri,rj,s,tmp) schedule(dynamic) nowait
    for(idx=0;idx<NCisAjs;idx++) {
      ri = CisAjsIdx[idx][0];
      rj = CisAjsIdx[idx][2];
      s  = CisAjsIdx[idx][3];
      tmp = GreenFunc1BF(ri,rj,s,ip,mySltBFTmp,myEleIdx,myEleCfg,myEleNum,eleProjCnt,myProjCntNew,eleProjBFCnt,myProjBFCntNew,myBuffer);
      LocalCisAjs[idx] = tmp;
    }

#pragma omp master
    {StopTimer(50);StartTimer(51);}

#pragma omp for private(idx,ri,rj,s,rk,rl,t,tmp) schedule(dynamic)
    for(idx=0;idx<NCisAjsCktAltDC;idx++) {
      ri = CisAjsCktAltDCIdx[idx][0];
      rj = CisAjsCktAltDCIdx[idx][2];
      s  = CisAjsCktAltDCIdx[idx][1];
      rk = CisAjsCktAltDCIdx[idx][4];
      rl = CisAjsCktAltDCIdx[idx][6];
      t  = CisAjsCktAltDCIdx[idx][5];
      tmp = GreenFunc2(ri,rj,rk,rl,s,t,ip,myEleIdx,eleCfg,myEleNum,eleProjCnt,
                       myProjCntNew,myBuffer);
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
      PhysCisAjsCktAlt[idx] += w*LocalCisAjs[idx0]*LocalCisAjs[idx1];
    }

#pragma omp master
    {StopTimer(53);}
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  return;
}
#endif
