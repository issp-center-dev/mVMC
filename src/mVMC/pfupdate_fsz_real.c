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
 * fast Pfaffian update
 *-------------------------------------------------------------
 *-------------------------------------------------------------*/

#include "pfupdate_fsz_real.h"
#ifndef _PFUDATE_FSZ_REAL_SRC
#define _PFUDATE_FSZ_REAL_SRC

/* Calculate new pfaffian. The ma-th electron with spin s hops. */
void CalculateNewPfM_fsz_real(const int ma, const int s, double *pfMNew, const int *eleIdx,const int *eleSpn,
                     const int qpStart, const int qpEnd) {
  #pragma procedure serial
  const int qpNum = qpEnd-qpStart;
  //const int msa = ma+s*Ne;//fsz
  const int msa = ma;
  const int rsa = eleIdx[msa] + s*Nsite;

  int qpidx;
  int msj,rsj;
  const double *sltE_a; /* update elements of msa-th row */
  const double *invM_a;
  double ratio;

  /* optimization for Kei */
  const int nsize = Nsize;

  #pragma loop noalias
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    sltE_a = SlaterElm_real + (qpidx+qpStart)*Nsite2*Nsite2 + rsa*Nsite2;
    invM_a = InvM_real + qpidx*Nsize*Nsize + msa*Nsize;

    ratio = 0.0;
    for(msj=0;msj<nsize;msj++) { //fsz
      rsj = eleIdx[msj]+eleSpn[msj]*Nsite;//fsz
      ratio += invM_a[msj] * sltE_a[rsj];
    }
    
    pfMNew[qpidx] = -ratio*PfM_real[qpidx];
  }

  return;
}

/* thread parallel version of CalculateNewPfM_fsz */
void CalculateNewPfM2_fsz_real(const int ma, const int s, double *pfMNew, const int *eleIdx,const int *eleSpn,
                     const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  //const int msa = ma+s*Ne;
  const int msa = ma;
  const int rsa = eleIdx[msa] + s*Nsite;

  int qpidx;
  int msj,rsj;
  const double *sltE_a; /* update elements of msa-th row */
  const double *invM_a;
  double ratio;

  /* optimization for Kei */
  const int nsize = Nsize;

  #pragma omp parallel for default(shared)        \
    private(qpidx,msj,sltE_a,invM_a,ratio,rsj)
  #pragma loop noalias
  for(qpidx=0;qpidx<qpNum;qpidx++) {
    sltE_a = SlaterElm_real + (qpidx+qpStart)*Nsite2*Nsite2 + rsa*Nsite2;
    invM_a = InvM_real + qpidx*Nsize*Nsize + msa*Nsize;

    ratio = 0.0;
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj]+eleSpn[msj]*Nsite;//fsz
      ratio += invM_a[msj] * sltE_a[rsj];
    }

    pfMNew[qpidx] = -ratio*PfM_real[qpidx];
  }

  return;
}

/* Update PfM and InvM. The ma-th electron with spin s hops to site ra=eleIdx[msi] */
void UpdateMAll_fsz_real(const int ma, const int s, const int *eleIdx,const int *eleSpn,
                const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  int qpidx;
  double *vec1,*vec2;

  RequestWorkSpaceThreadDouble(2*Nsize);

  #pragma omp parallel default(shared) private(vec1,vec2)
  {
    vec1 = GetWorkSpaceThreadDouble(Nsize);
    vec2 = GetWorkSpaceThreadDouble(Nsize);
   
    #pragma omp for private(qpidx)
    #pragma loop nounroll
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      updateMAll_child_fsz_real(ma, s, eleIdx,eleSpn, qpStart, qpEnd, qpidx, vec1, vec2);
    }
  }

  ReleaseWorkSpaceThreadDouble();
  return;
}

void updateMAll_child_fsz_real(const int ma, const int s, const int *eleIdx,const int *eleSpn,
                      const int qpStart, const int qpEnd, const int qpidx,
                      double *vec1, double *vec2) {
  #pragma procedure serial
  /* const int qpNum = qpEnd-qpStart; */
  //const int msa = ma+s*Ne;
  const int msa = ma;
  const int rsa = eleIdx[msa] + s*Nsite;
  const int nsize = Nsize; /* optimization for Kei */

  int msi,msj,rsj;

  const double *sltE_a; /* update elements of msa-th row */
  double sltE_aj;
  double *invM;
  double *invM_i,*invM_j,*invM_a;

  double vec1_i,vec2_i;
  double invVec1_a;
  double tmp;

  sltE_a = SlaterElm_real + (qpidx+qpStart)*Nsite2*Nsite2 + rsa*Nsite2;

  invM = InvM_real + qpidx*Nsize*Nsize;
  invM_a = invM + msa*Nsize;

  for(msi=0;msi<nsize;msi++) vec1[msi] = 0.0+0.0*I; //TBC

  /* Calculate vec1[i] = sum_j invM[i][j] sltE[a][j] */
  /* Note that invM[i][j] = -invM[j][i] */
  #pragma loop noalias
  for(msj=0;msj<nsize;msj++) {
    //rsj = eleIdx[msj] + (msj/Ne)*Nsite;
    rsj = eleIdx[msj] + eleSpn[msj]*Nsite; //fsz
    sltE_aj = sltE_a[rsj];
    invM_j = invM + msj*Nsize;

    for(msi=0;msi<nsize;msi++) {
      vec1[msi] += -invM_j[msi] * sltE_aj;
    }
  }

  /* Update Pfaffian */
  /* Calculate -1.0/bufV_a to reduce devision */
  tmp = vec1[msa];
  PfM_real[qpidx] *= -tmp;
  invVec1_a = -1.0/tmp;

  /* Calculate vec2[i] = -InvM[a][i]/vec1[a] */
  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    vec2[msi] = invM_a[msi] * invVec1_a;
  }

  /* Update InvM */
  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    invM_i = invM + msi*Nsize;
    vec1_i = vec1[msi];
    vec2_i = vec2[msi];

    for(msj=0;msj<nsize;msj++) {
      invM_i[msj] += vec1_i * vec2[msj] - vec1[msj] * vec2_i;
    }

    invM_i[msa] -= vec2_i;
  }

  #pragma loop noalias
  for(msj=0;msj<nsize;msj++) {
    invM_a[msj] += vec2[msj];
  }
  /* end of update invM */

  return;
}

#endif
