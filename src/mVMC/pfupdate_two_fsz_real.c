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
 * fast Pfaffian update for two-electron hopping
 *-------------------------------------------------------------
 *-------------------------------------------------------------*/

#include "pfupdate_two_fsz_real.h"
#ifndef _PFUDATE_TWO_FSZ_REAL_SRC
#define _PFUDATE_TWO_FSZ_REAL_SRC

/* Calculate new pfaffian. 
   The ma-th electron with spin s hops
   and then the mb-th electron with spin t hops. */
/* buffer size = 2*Nsize */
void CalculateNewPfMTwo_fsz_real(const int ma, const int s, const int mb, const int t,
                        double *pfMNew, const int *eleIdx,const int *eleSpn,
                        const int qpStart, const int qpEnd, double *buffer) {
  #pragma procedure serial
  const int qpNum = qpEnd-qpStart;
  const int msa = ma;//+s*Ne; fsz
  const int msb = mb;//+t*Ne; fsz
  int qpidx;
  double *vec_a = buffer;
  double *vec_b = buffer + Nsize;

  if(msa==msb) {
    CalculateNewPfM_fsz_real(mb, t, pfMNew, eleIdx,eleSpn, qpStart, qpEnd);
    return;
  }

  for(qpidx=0;qpidx<qpNum;qpidx++) {
    calculateNewPfMTwo_child_fsz_real(ma, s, mb, t, pfMNew, eleIdx,eleSpn,
                             qpStart, qpEnd, qpidx, vec_a, vec_b);
  }

  return;
}

/* Calculate new pfaffian. 
   The ma-th electron with spin s hops
   and then the mb-th electron with spin t hops. */
/* thread parallel version of CalculateNewPfMTwo */
void CalculateNewPfMTwo2_fsz_real(const int ma, const int s, const int mb, const int t,
                         double *pfMNew, const int *eleIdx,const int *eleSpn,
                         const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  const int msa = ma;//+s*Ne; fsz
  const int msb = mb;//+t*Ne; fsz
  int qpidx;
  double *vec_a,*vec_b;

  if(msa==msb) {
    CalculateNewPfM2_fsz_real(mb, t, pfMNew, eleIdx,eleSpn, qpStart, qpEnd);
    return;
  }

  RequestWorkSpaceThreadDouble(2*Nsize);

  #pragma omp parallel default(shared) private(vec_a,vec_b)
  {
    vec_a = GetWorkSpaceThreadDouble(Nsize);
    vec_b = GetWorkSpaceThreadDouble(Nsize);

    #pragma omp for private(qpidx)
    #pragma loop nounroll
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      calculateNewPfMTwo_child_fsz_real(ma, s, mb, t, pfMNew, eleIdx,eleSpn,
                               qpStart, qpEnd, qpidx, vec_a, vec_b);
    }
  }
  
  ReleaseWorkSpaceThreadDouble();
  return;
}


void calculateNewPfMTwo_child_fsz_real(const int ma, const int s, const int mb, const int t,
                              double *pfMNew, const int *eleIdx,const int *eleSpn,
                              const int qpStart, const int qpEnd, const int qpidx,
                              double *vec_a, double *vec_b) {
  #pragma procedure serial
  /* const int qpNum = qpEnd-qpStart; */
  const int msa = ma;//+s*Ne; fsz
  const int msb = mb;//+t*Ne; fsz
  const int rsa = eleIdx[msa] + s*Nsite;
  const int rsb = eleIdx[msb] + t*Nsite;
  const int nsize = Nsize;

  int msi,msj;
  int rsi;

  double p_a,p_b,q_a,q_b,bMa;
  double ratio,tmp;

  const double *sltE;
  const double *sltE_a; /* update elements of msa-th row */
  const double *sltE_b; /* update elements of msb-th row */

  double *invM;
  const double *invM_a, *invM_b, *invM_i;
  double invM_ab,invM_ai,invM_bi;

  double vec_ba,vec_ai,vec_bi;

  sltE = SlaterElm_real + (qpidx+qpStart)*Nsite2*Nsite2;
  sltE_a = sltE + rsa*Nsite2;
  sltE_b = sltE + rsb*Nsite2;

  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    rsi = eleIdx[msi] + eleSpn[msi]*Nsite;// fsz(msi/Ne)*Nsite;
    vec_a[msi] = sltE_a[rsi];
    vec_b[msi] = sltE_b[rsi];
  }
  vec_ba = vec_b[msa];

  invM = InvM_real + qpidx*Nsize*Nsize;
  invM_a = invM + msa*Nsize;
  invM_b = invM + msb*Nsize;
  invM_ab = invM_a[msb];

  p_a = p_b = q_a = q_b = bMa = 0.0;

  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    vec_ai = vec_a[msi];
    vec_bi = vec_b[msi];
    invM_ai = invM_a[msi];
    invM_bi = invM_b[msi];

    p_a += invM_ai * vec_ai;
    p_b += invM_bi * vec_ai;
    q_a += invM_ai * vec_bi;
    q_b += invM_bi * vec_bi;
  }

  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    invM_i = invM + msi*Nsize;
    tmp = 0.0;
    for(msj=0;msj<nsize;msj++) {
      tmp += invM_i[msj] * vec_a[msj];
    }
    bMa += vec_b[msi] * tmp;
  }

  /* Calculate ratio PfMNew/PfMOld */
  ratio = invM_ab*vec_ba + invM_ab*bMa + p_a*q_b - p_b*q_a;

  /* Update pfMNew */
  pfMNew[qpidx] = ratio*PfM_real[qpidx];

  return;
}

// s comp
/* Update PfM and InvM. The ma-th electron with spin s hops from raOld to site ra=eleIdx[msa],
   and then the mb-th electron with spin t hops from rbOld to site rb=eleIdx[msb] */
void UpdateMAllTwo_fsz_real(const int ma, const int s, const int mb, const int t,
                   const int raOld, const int rbOld,
                   const int *eleIdx,const int *eleSpn, const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  int qpidx;
  double *vec1,*vec2,*vec3,*vec4;

  RequestWorkSpaceThreadDouble(4*Nsize);

#pragma omp parallel default(shared) private(vec1,vec2,vec3,vec4,qpidx)
  {
    vec1 = GetWorkSpaceThreadDouble(Nsize);
    vec2 = GetWorkSpaceThreadDouble(Nsize);
    vec3 = GetWorkSpaceThreadDouble(Nsize);
    vec4 = GetWorkSpaceThreadDouble(Nsize);
   
    #pragma omp for
    #pragma loop nounroll
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      updateMAllTwo_child_fsz_real(ma, s, mb, t, raOld, rbOld, eleIdx,eleSpn, qpStart, qpEnd, qpidx,
                          vec1, vec2, vec3, vec4);
    }
  }

  ReleaseWorkSpaceThreadDouble();
  return;
}

void updateMAllTwo_child_fsz_real(const int ma, const int s, const int mb, const int t,
                         const int raOld, const int rbOld,
                         const int *eleIdx,const int *eleSpn, const int qpStart, const int qpEnd, const int qpidx,
                         double *vecP, double *vecQ, double *vecS, double *vecT) {
  #pragma procedure serial
  const int msa = ma;//+s*Ne; fsz
  const int msb = mb;//+t*Ne; fsz
  const int rsa = eleIdx[msa] + s*Nsite;
  const int rsb = eleIdx[msb] + t*Nsite;
  const int rsaOld = raOld + s*Nsite;
  const int rsbOld = raOld + t*Nsite;
  const int nsize = Nsize;

  const double *sltE = SlaterElm_real + (qpidx+qpStart)*Nsite2*Nsite2;;
  const double *sltE_a = sltE + rsa*Nsite2; /* update elements of msa-th row */
  const double *sltE_b = sltE + rsb*Nsite2; /* update elements of msb-th row */
  const double mOld_ab = sltE[rsaOld*Nsite2 + rsbOld];

  double *invM = InvM_real + qpidx*Nsize*Nsize;
  double *invM_a = invM + msa*Nsize;
  double *invM_b = invM + msb*Nsize;
  double *invM_i;
  double invM_ab = invM_a[msb];

  double ratio,det,invDet,bMa;
  double a,b,c,d,e,f;
  double p_i,p_j,q_i,q_j,s_i,s_j,t_i,t_j;

  int msi,msj;
  int rsi;

  /* initialize vecP[i], vecQ[i] */
  /* vecS[i], vecT[i] are temporally used as
     vecS[i] = sltE[a][j], vecT[i] = sltE[b][j]. */
  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    vecP[msi]=0.0;
    vecQ[msi]=0.0;
    rsi = eleIdx[msi] + eleSpn[msi]*Nsite;//(msi/Ne)*Nsite; //fsz
    vecS[msi]=sltE_a[rsi];
    vecT[msi]=sltE_b[rsi];
  }
  /* Set vecS[b] = mOld_ab, which is (a,b)-elements of the old M. */
  vecS[msb] = mOld_ab;

  /* Calculate vecP[i], vecQ[i] */
  /* vecP[i]= sum_j invM[i][j]*sltE[a][j] */
  /* vecQ[i]= sum_j invM[i][j]*sltE[b][j] */
  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    invM_i = invM + msi*Nsize;
    for(msj=0;msj<nsize;msj++) {
      vecP[msi] += invM_i[msj]*vecS[msj];
      vecQ[msi] += invM_i[msj]*vecT[msj];
    }
  }

  /* Update Paffian */
  bMa = 0.0;
  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    bMa += vecT[msi] * vecP[msi];
  }
  ratio = invM_ab*vecT[msa] + invM_ab*bMa
        + vecP[msa]*vecQ[msb] - vecP[msb]*vecQ[msa];
  PfM_real[qpidx] *= ratio;

  /* Set coefficients */
  a = -vecP[msa];
  b = vecP[msb];
  c = vecQ[msa];
  d = -vecQ[msb];
  e = -bMa - vecT[msa];
  f = invM_a[msb];

  det = a*d - b*c - e*f;
  invDet = 1.0/det;

  /* Calculate vecS[i], vecT[i] */
  /* vecS[i]= invM[a][i]/D, vecT[i]= invM[b][i]/D */
  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    vecS[msi] = invDet * invM_a[msi];
    vecT[msi] = invDet * invM_b[msi];
  }  

  /* Update InvM */
  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    invM_i = invM + msi*Nsize;
    p_i = vecP[msi];
    q_i = vecQ[msi];
    s_i = vecS[msi];
    t_i = vecT[msi];

    for(msj=0;msj<nsize;msj++) {
      p_j = vecP[msj];
      q_j = vecQ[msj];
      s_j = vecS[msj];
      t_j = vecT[msj];

      invM_i[msj] += a*(q_i * t_j - q_j * t_i)
        + b*(q_i * s_j - q_j * s_i)
        + c*(p_i * t_j - p_j * t_i)
        + d*(p_i * s_j - p_j * s_i)
        + e*det*(s_i * t_j - s_j * t_i)
        + f*invDet*(p_i * q_j - q_i * p_j);
    }
    invM_i[msa] += -c*t_i -d*s_i - f*invDet*q_i;
    invM_i[msb] += -a*t_i -b*s_i + f*invDet*p_i;
  }

  #pragma loop noalias
  for(msj=0;msj<nsize;msj++) {
    p_j = vecP[msj];
    q_j = vecQ[msj];
    s_j = vecS[msj];
    t_j = vecT[msj];

    invM_a[msj] += c*t_j + d*s_j + f*invDet*q_j;
    invM_b[msj] += a*t_j + b*s_j - f*invDet*p_j;
  }

  invM_a[msb] += f*invDet;
  invM_b[msa] -= f*invDet;

  return;
}

//e comp
#endif
