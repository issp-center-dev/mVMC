/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Mitsuaki Kawamura, Takeo Kato, Masatoshi Imada.

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
 * local Green Functions
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

#include "locgrn_real.h"
#include "pfupdate_real.h"
#include "pfupdate_two_real.h"

#ifndef _SRC_LOCGRN_REAL
#define _SRC_LOCGRN_REAL

double calculateNewPfMN_real_child(const int qpidx, const int n, const int *msa, const int *rsa,
                                   const int *eleIdx, double *buffer);

/* Calculate 1-body Green function <CisAjs> */
/* buffer size = NQPFull */
double  GreenFunc1_real(const int ri, const int rj, const int s, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  int *projCntNew, double *buffer) {
  double  z;
  int mj,msj,rsi,rsj;
  double  *pfMNew_real = buffer; /* NQPFull */

  if(ri==rj) return eleNum[ri+s*Nsite];
  if(eleNum[ri+s*Nsite]==1 || eleNum[rj+s*Nsite]==0) return 0.0;

  mj = eleCfg[rj+s*Nsite];
  msj = mj + s*Ne;
  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;

  /* hopping */
  eleIdx[msj] = ri;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt(rj, ri, s, projCntNew, eleProjCnt, eleNum);
  z = ProjRatio(projCntNew,eleProjCnt);

  /* calculate Pfaffian */
  CalculateNewPfM_real(mj, s, pfMNew_real, eleIdx, 0, NQPFull);
  z *= CalculateIP_real(pfMNew_real, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  return z/ip;//TBC
}

/* Calculate 2-body Green function <psi|CisAjsCktAlt|x>/<psi|x> */
/* buffer size = NQPFull+2*Nsize */
double GreenFunc2_real(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  int *projCntNew, double *buffer) {
  double z;
  int mj,msj,ml,mtl;
  int rsi,rsj,rtk,rtl;
  double *pfMNew_real = buffer; /* [NQPFull] */
  double *bufV   = buffer+NQPFull; /* 2*Nsize */

  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;
  rtk = rk + t*Nsite;
  rtl = rl + t*Nsite;

  if(s==t) {
    if(rk==rl) { /* CisAjsNks */
      if(eleNum[rtk]==0) return 0.0;
      else return GreenFunc1_real(ri,rj,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,projCntNew,buffer); /* CisAjs */
    }else if(rj==rl) {
      return 0.0; /* CisAjsCksAjs (j!=k) */
    }else if(ri==rl) { /* AjsCksNis */
      if(eleNum[rsi]==0) return 0.0;
      else if(rj==rk) return 1.0-eleNum[rsj];
      else return -GreenFunc1_real(rk,rj,s,ip,eleIdx,eleCfg,eleNum,
                              eleProjCnt,projCntNew,buffer); /* -CksAjs */
    }else if(rj==rk) { /* CisAls(1-Njs) */
      if(eleNum[rsj]==1) return 0.0;
      else if(ri==rl) return eleNum[rsi];
      else return GreenFunc1_real(ri,rl,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,projCntNew,buffer); /* CisAls */
    }else if(ri==rk) {
      return 0.0; /* CisAjsCisAls (i!=j) */
    }else if(ri==rj) { /* NisCksAls (i!=k,l) */
      if(eleNum[rsi]==0) return 0.0;
      else return GreenFunc1_real(rk,rl,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,projCntNew,buffer); /* CksAls */
    }
  } if(s!=t) {
    if(rk==rl) { /* CisAjsNkt */
      if(eleNum[rtk]==0) return 0.0;
      else if(ri==rj) return eleNum[rsi];
      else return GreenFunc1_real(ri,rj,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,projCntNew,buffer); /* CisAjs */
    }else if(ri==rj) { /* NisCktAlt */
      if(eleNum[rsi]==0) return 0.0;
      else return GreenFunc1_real(rk,rl,t,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,projCntNew,buffer); /* CktAlt */
    }
  }

  if(eleNum[rsi]==1 || eleNum[rsj]==0 || eleNum[rtk]==1 || eleNum[rtl]==0) return 0.0;

  mj = eleCfg[rj+s*Nsite];
  ml = eleCfg[rl+t*Nsite];
  msj = mj + s*Ne;
  mtl = ml + t*Ne;

  /* hopping */
  eleIdx[mtl] = rk;
  eleNum[rtl] = 0;
  eleNum[rtk] = 1;
  UpdateProjCnt(rl, rk, t, projCntNew, eleProjCnt, eleNum);
  eleIdx[msj] = ri;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt(rj, ri, s, projCntNew, projCntNew, eleNum);

  z = ProjRatio(projCntNew,eleProjCnt);

  /* calculate Pfaffian */
  CalculateNewPfMTwo_real(ml, t, mj, s, pfMNew_real, eleIdx, 0, NQPFull, bufV);
  z *= CalculateIP_real(pfMNew_real, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[mtl] = rl;
  eleNum[rtl] = 1;
  eleNum[rtk] = 0;
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  return z/ip;//TBC
}


/// Calculate n-body Green function
/// <phi| c1 a1 c2 a2 ... cn an |x>
/// c_k = c_rsi[k], a_k = a_rsj[k]
/// buffer size    = NQPFull + n*Nsize
/// bufferInt size = NProj
/// \param n
/// \param rsi the array of indices of creation operators
/// \param rsj the array of indices of annihilation operators
/// \param ip
/// \param eleIdx
/// \param eleCfg
/// \param eleNum
/// \param eleProjCnt
/// \param buffer
/// \param bufferInt
/// \return
/// \version 1.0
double GreenFuncN_real(const int n, int *rsi, int *rsj, const double  ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  double *buffer, int *bufferInt){
  int ri,rj,rk,rl,si,sj,sk,mj;
  int k,l,m,rsk;
  double z,x;
  int qpidx;

  int *projCntNew = bufferInt; /* [NProj] */

  int msj[n];
  double *pfMNew = buffer; /* [NQPFull] */
  double *bufV = buffer+NQPFull; /* [n*Nsize] */

  for(k=0;k<n;k++) {
    si = rsi[k]/Nsite;
    sj = rsj[k]/Nsite;
    if(si!=sj) return 0;
  }

  if(n<=0) return 0;
  else if(n==1) {
    ri = rsi[0]%Nsite;
    rj = rsj[0]%Nsite;
    si = rsi[0]/Nsite;
    return GreenFunc1_real(ri,rj,si,ip,eleIdx,eleCfg,eleNum,eleProjCnt,projCntNew,buffer);
  } else if(n==2) {
    ri = rsi[0]%Nsite;
    rj = rsj[0]%Nsite;
    si = rsi[0]/Nsite;
    rk = rsi[1]%Nsite;
    rl = rsj[1]%Nsite;
    sk = rsi[1]/Nsite;
    return GreenFunc2_real(ri,rj,rk,rl,si,sk,ip,eleIdx,eleCfg,eleNum,eleProjCnt,projCntNew,buffer);
  }

  /* reduction */
  for(k=n-1;k>=0;k--) {
    /* ** check for an annihilation operator at rsj[k] ** */
    rsk = rsj[k];
    for(l=k+1;l<n;l++) {
      /* rsj[k] == rsi[l] */
      if(rsk==rsi[l]) {
        rsj[k] = rsj[l];
        for(m=l;m<n-1;m++) { /* shift */
          rsi[m] = rsi[m+1];
          rsj[m] = rsj[m+1];
        }
        return GreenFuncN_real(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,eleProjCnt,buffer,bufferInt);
      }
      /* rsj[k] == rsj[l] */
      if(rsk==rsj[l]) return 0;
    }
    /* check electron number */
    if(eleNum[rsk]==0) return 0;

    /* ** check for a creation operator at rsi[k] ** */
    rsk = rsi[k];
    /* rsi[k] == rsj[k] */
    if(rsk==rsj[k]) {
      for(m=k;m<n-1;m++) { /* shift */
        rsi[m] = rsi[m+1];
        rsj[m] = rsj[m+1];
      }
      return GreenFuncN_real(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,eleProjCnt,buffer,bufferInt);
    }
    for(l=k+1;l<n;l++) {
      /* rsi[k] == rsi[l] */
      if(rsk==rsi[l]) return 0;
      /* rsi[k] == rsj[l] (k<l) */
      if(rsk==rsj[l]) {
        rsi[k] = rsi[l];
        for(m=l;m<n-1;m++) { /* shift */
          rsi[m] = rsi[m+1];
          rsj[m] = rsj[m+1];
        }
        return (-1.0)*GreenFuncN_real(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,eleProjCnt,buffer,bufferInt);
      }
    }
    /* check electron number */
    if(eleNum[rsk]==1) return 0;
  }

  /* hopping */
  #pragma loop noalias
  for(k=0;k<n;k++) {
    ri = rsi[k]%Nsite;
    rj = rsj[k]%Nsite;
    sj = rsj[k]/Nsite;
    mj = eleCfg[rsj[k]];
    msj[k] = mj + sj*Ne;

    eleIdx[msj[k]] = ri;
    eleNum[rsj[k]] = 0;
    eleNum[rsi[k]] = 1;
  }

  MakeProjCnt(projCntNew,eleNum);
  x = ProjRatio(projCntNew,eleProjCnt);

  /* calculateNewPfM */
  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    pfMNew[qpidx] = calculateNewPfMN_real_child(qpidx,n,msj,rsj,eleIdx,bufV);
  }
  z = CalculateIP_real(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hoppint */
  #pragma loop noalias
  for(k=0;k<n;k++) {
    rj = rsj[k]%Nsite;
    eleIdx[msj[k]] = rj;
    eleNum[rsj[k]] = 1;
    eleNum[rsi[k]] = 0;
  }

  return x*z/ip;
}

///* msa[k]-th electron hops from rsa[k] to eleIdx[msa[k]] */
///* buffer size = n*Nsize */
double calculateNewPfMN_real_child(const int qpidx, const int n, const int *msa, const int *rsa,
                              const int *eleIdx, double *buffer) {
  const int nsize = Nsize;
  const int n2 = 2*n;
  const double  *sltE;
  const double  *sltE_k;
  const double  *invM;
  const double  *invM_i, *invM_k, *invM_l;

  double  *vec; /* vec[n][nsize] */
  double  *vec_k, *vec_l;
  double  mat[n2*n2]; /* mat[n2][n2] */
  double  *mat_k;
  double sgn;

  int rsi,rsk,msi,msj,k,l;
  double val,tmp;

  /* for DSKPFA */
  char uplo='U', mthd='P';
  int nn,lda,info=0;
  double pfaff;
  int iwork[n2];
  double work[n2*n2]; /* [n2][n2] */
  int lwork = n2*n2;
  nn=lda=n2;

  sltE = SlaterElm_real + qpidx*Nsite2*Nsite2;
  invM = InvM_real + qpidx*Nsize*Nsize;

  vec = buffer; /* n*nsize */

  #pragma loop noalias
  for(k=0;k<n;k++) {
    rsk = eleIdx[msa[k]] + (msa[k]/Ne)*Nsite;
    sltE_k = sltE + rsk*Nsite2;
    vec_k = vec + k*nsize;
    #pragma loop norecurrence
    for(msi=0;msi<nsize;msi++) {
      rsi = eleIdx[msi] + (msi/Ne)*Nsite;
      vec_k[msi] = sltE_k[rsi];
    }
  }

  /* X_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k;
    vec_k = vec + k*nsize;
    for(l=k+1;l<n;l++) {
      vec_l = vec + l*nsize;
      val = 0.0;
      for(msi=0;msi<nsize;msi++) {
        invM_i = invM + msi*nsize;
        tmp = 0.0;
        for(msj=0;msj<nsize;msj++) {
          tmp += invM_i[msj] * vec_l[msj];
        }
        val += tmp * vec_k[msi];
      }
      mat_k[l] = val + vec_k[msa[l]];
    }
  }

  /* Y_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*k + n;
    vec_k = vec + k*nsize;
    for(l=0;l<n;l++) {
      invM_l = invM + msa[l]*nsize;
      val = 0.0;
      for(msi=0;msi<nsize;msi++) {
        val += vec_k[msi] * invM_l[msi];
      }
      mat_k[l] = val;
    }
  }

  /* Z_kl */
  for(k=0;k<n;k++) {
    mat_k = mat + n2*(k+n) + n;
    invM_k = invM + msa[k]*nsize;
    for(l=k+1;l<n;l++) {
      mat_k[l] = invM_k[msa[l]];
    }
  }

  #pragma loop noalias
  for(k=0;k<n2;k++) {
    #pragma loop norecurrence
    for(l=0;l<k;l++) {
      mat[n2*k + l] = -mat[n2*l + k]; /* transpose */
    }
    mat[n2*k + k] = 0.0; /* diagonal elements */
  }

  /* calculate Pf M */
  M_DSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, &info);

  sgn = ( (n*(n-1)/2)%2==0 ) ? 1.0 : -1.0;

  return sgn * pfaff * PfM_real[qpidx];
}


/* Calculate 1-body Green function <CisAjs> */
/* buffer size = NQPFull */
double GreenFunc1BF_real(const int ri, const int rj, const int s, const double ip, double *bufM,
                    int *eleIdx, int *eleCfg, int *eleNum, const int *eleProjCnt,
                    int *projCntNew, const int *eleProjBFCnt,int *projBFCntNew, double *buffer) {
  double z,logz;
  int msi,mj,msj,rsi,rsj;
  int qpidx,i;
  const int nsize=Nsize;
  double *sltE;
  double *sltE_i;
  double *bufM_i, *bufM_i2;
  //double complex bufM[NQPFull*Nsize*Nsize];
  double *pfMNew_real = buffer; /* NQPFull */
  int msaTmp[NQPFull*Nsite],icount[NQPFull];

  if(ri==rj) return eleNum[ri+s*Nsite];
  if(eleNum[ri+s*Nsite]==1 || eleNum[rj+s*Nsite]==0) return 0.0;

  /* store SlaterElmBF before hopping */
  mj = eleCfg[rj+s*Nsite];
  msj = mj + s*Ne;
  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;

  /* hopping */
  eleCfg[rsj] = -1;
  eleCfg[rsi] = mj;
  eleIdx[msj] = ri;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt(rj, ri, s, projCntNew, eleProjCnt, eleNum);
  z = ProjRatio(projCntNew,eleProjCnt);

  /* calculate Pfaffian */
  //printf("1");
  StartTimer(81);
  MakeProjBFCnt(projBFCntNew, eleNum);
  StopTimer(81);
  StartTimer(82);
  UpdateSlaterElmBFGrn_real(mj, rj, ri, s, eleCfg, eleNum, projBFCntNew, msaTmp, icount, bufM);
  StopTimer(82);
  StartTimer(83);
  CalculateNewPfMBF_real(icount, msaTmp, pfMNew_real, eleIdx, 0, NQPFull, bufM);
  StopTimer(83);
  z *= CalculateIP_real(pfMNew_real, 0, NQPFull, MPI_COMM_SELF);
  //printf("1");
  //MakeSlaterElmBF(eleNum,projBFCntNew);
  //UpdateSlaterElmBF2(eleNum,projBFCntNew);
  //printf("1");
  //CalculateMAllBF_NoThread(eleIdx,0,NQPFull,bufM);
  //CalculateMAllBF_NoThread(eleIdx,0,NQPFull,SlaterElmBF);
  //printf("1\n");
  //z *= CalculateIP(PfM,0,NQPFull,MPI_COMM_SELF);

  /* revert hopping */
  eleCfg[rsj] = mj;
  eleCfg[rsi] = -1;
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  //UpdateSlaterElmBFTmp5(mj, ri, rj, s, eleCfg, eleNum, eleProjBFCnt, msaTmp, icount, bufM);
  StoreSlaterElmBF_real(bufM);
  //MakeProjBFCnt(projBFCntNew, eleNum);
  //MakeSlaterElmBF(eleNum,projBFCntNew);
  //CalculateMAllBF_NoThread(eleIdx,0,NQPFull,bufM);
  //CalculateMAllBF_NoThread(eleIdx,0,NQPFull,SlaterElmBF);

  return z/ip;
}


#endif
