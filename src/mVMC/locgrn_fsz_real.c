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
 * local Green Functions
 *-------------------------------------------------------------
 *-------------------------------------------------------------*/
#include "locgrn_fsz_real.h"
#include "projection.h"
#include "pfupdate_fsz_real.h"
#include "pfupdate_two_fsz_real.h"
#include "qp_real.h"

#pragma once

/* Calculate 1-body Green function <CisAjs> */
/* buffer size = NQPFull */
double GreenFunc1_fsz_real(const int ri, const int rj, const int s, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double *buffer) {
  double z;
  int mj,msj,rsi,rsj;
  double *pfMNew = buffer; /* NQPFull */

  if(ri==rj) return eleNum[ri+s*Nsite];
  if(eleNum[ri+s*Nsite]==1 || eleNum[rj+s*Nsite]==0) return 0.0;

  mj  = eleCfg[rj+s*Nsite];
  msj = mj;// + s*Ne;
  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;

  /* hopping */
  eleIdx[msj] = ri;
  eleSpn[msj] = s;//fsz
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt(rj, ri, s, projCntNew, eleProjCnt, eleNum);
  z = ProjRatio(projCntNew,eleProjCnt);

  /* calculate Pfaffian */
  CalculateNewPfM_fsz_real(mj, s, pfMNew, eleIdx,eleSpn, 0, NQPFull);//fsz
  z *= CalculateIP_real(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[msj] = rj;
  eleSpn[msj] = s; //fsz
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  return z/ip;//TBC
}

/* Calculate 1-body Green function <CisAjt> */ // s!=t
/* buffer size = NQPFull */
double GreenFunc1_fsz2_real(const int ri, const int rj, const int s,const int t, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double *buffer) {
  double z;
  int mj,rsi,rtj;//msj
  double *pfMNew = buffer; /* NQPFull */

  //if(ri==rj) return eleNum[ri+s*Nsite]; //fsz
  if(eleNum[ri+s*Nsite]==1 || eleNum[rj+t*Nsite]==0) return 0.0;

  mj  = eleCfg[rj+t*Nsite];
  //mtj = mj;// + s*Ne;
  rsi = ri + s*Nsite;
  rtj = rj + t*Nsite;

  /* hopping */
  // (j,t) -> (i,s)
  eleIdx[mj] = ri;
  eleSpn[mj] = s;//fsz
  eleNum[rtj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt_fsz(rj, ri, t,s, projCntNew, eleProjCnt, eleNum);
  z = ProjRatio(projCntNew,eleProjCnt);

  /* calculate Pfaffian */
  CalculateNewPfM_fsz_real(mj, s, pfMNew, eleIdx,eleSpn, 0, NQPFull);//fsz: note EleSpn[mj]=s
  z *= CalculateIP_real(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[mj] = rj;
  eleSpn[mj] = t; //fsz
  eleNum[rtj] = 1;
  eleNum[rsi] = 0;

  return z/ip;//TBC
}


/* Calculate 2-body Green function <psi|CisAjsCktAlt|x>/<psi|x> */
/* buffer size = NQPFull+2*Nsize */
double GreenFunc2_fsz_real(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double *buffer) {
  double z;
  int mj,msj,ml,mtl;
  int rsi,rsj,rtk,rtl;
  double *pfMNew = buffer; /* [NQPFull] */
  double *bufV   = buffer+NQPFull; /* 2*Nsize */

  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;
  rtk = rk + t*Nsite;
  rtl = rl + t*Nsite;

  if(s==t) {
    if(rk==rl) { /* CisAjsNks */
      if(eleNum[rtk]==0) return 0.0;
      else return GreenFunc1_fsz_real(ri,rj,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CisAjs */
    }else if(rj==rl) {
      return 0.0; /* CisAjsCksAjs (j!=k) */
    }else if(ri==rl) { /* AjsCksNis */
      if(eleNum[rsi]==0) return 0.0;
      else if(rj==rk) return 1.0-eleNum[rsj];
      else return -GreenFunc1_fsz_real(rk,rj,s,ip,eleIdx,eleCfg,eleNum,
                              eleProjCnt,eleSpn,projCntNew,buffer); /* -CksAjs */
    }else if(rj==rk) { /* CisAls(1-Njs) */
      if(eleNum[rsj]==1) return 0.0;
      else if(ri==rl) return eleNum[rsi];
      else return GreenFunc1_fsz_real(ri,rl,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CisAls */
    }else if(ri==rk) {
      return 0.0; /* CisAjsCisAls (i!=j) */
    }else if(ri==rj) { /* NisCksAls (i!=k,l) */
      if(eleNum[rsi]==0) return 0.0;
      else return GreenFunc1_fsz_real(rk,rl,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CksAls */
    }
  }else{
    if(rk==rl) { /* CisAjsNkt */
      if(eleNum[rtk]==0) return 0.0;
      else if(ri==rj) return eleNum[rsi];
      else return GreenFunc1_fsz_real(ri,rj,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CisAjs */
    }else if(ri==rj) { /* NisCktAlt */
      if(eleNum[rsi]==0) return 0.0;
      else return GreenFunc1_fsz_real(rk,rl,t,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CktAlt */
    }
  }

  if(eleNum[rsi]==1 || eleNum[rsj]==0 || eleNum[rtk]==1 || eleNum[rtl]==0) return 0.0;

  mj = eleCfg[rj+s*Nsite];
  ml = eleCfg[rl+t*Nsite];
  msj = mj;// + s*Ne;
  mtl = ml;// + t*Ne;

  /* hopping */
  eleIdx[mtl] = rk;
  eleSpn[mtl] = t;
  eleNum[rtl] = 0;
  eleNum[rtk] = 1;
  UpdateProjCnt(rl, rk, t, projCntNew, eleProjCnt, eleNum);
  eleIdx[msj] = ri;
  eleSpn[msj] = s;
  eleNum[rsj] = 0;
  eleNum[rsi] = 1;
  UpdateProjCnt(rj, ri, s, projCntNew, projCntNew, eleNum);

  z = ProjRatio(projCntNew,eleProjCnt);

  /* calculate Pfaffian */
  CalculateNewPfMTwo_fsz_real(ml, t, mj, s, pfMNew, eleIdx,eleSpn, 0, NQPFull, bufV);
  z *= CalculateIP_real(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[mtl] = rl;
  eleSpn[mtl] = t;
  eleNum[rtl] = 1;
  eleNum[rtk] = 0;
  eleIdx[msj] = rj;
  eleSpn[msj] = s;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  return z/ip;//TBC
}

/* Calculate 2-body Green function <psi|CisAjtCkuAlv|x>/<psi|x> */
/* buffer size = NQPFull+2*Nsize */
double GreenFunc2_fsz2_real(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t,const int u,const int v, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double *buffer) {
  double z;
  //int mj,msj,ml,mtl;
  int mj,ml; // mj: (rj,t) -> (ri,s) ml:(rl,v) -> (rk,u)
  //int rsi,rtj,ruk,rvl;
  int XI,XJ,XK,XL; //fsz
  double *pfMNew = buffer; /* [NQPFull] */
  double *bufV   = buffer+NQPFull; /* 2*Nsize */

  XI = ri + s*Nsite;
  XJ = rj + t*Nsite;
  XK = rk + u*Nsite;
  XL = rl + v*Nsite;


  if(XI == XJ) {
    if(XJ == XK) {
      if(XK == XL) {
        // I=J=K=L #1
        return eleNum[XI];
      }else{
        // I=J=K !=L #2
        if(eleNum[XI]==1) return 0.0; // 
        else return GreenFunc1_fsz2_real(rk,rl,u,v,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CkuAlv */
      }      
    }else if(XJ == XL){
      // I=J=L !=K #3
      return 0.0;
    }else if(XK == XL){
      // I=J ! K=L #4
      return eleNum[XI]*eleNum[XK];
    }else {
      // I=J   K!=L #5
      if(eleNum[XI]==0) return 0.0; // 
      else return GreenFunc1_fsz2_real(rk,rl,u,v,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer); /* CkuAlv */
    } 
  }else if(XI == XK){
    if(XJ == XL){
      // I=K != J=L  #6
      return 0.0;
    }else if(XK == XL){
      // I=K=L !=J #7
      return 0.0;
    }else{
      // I=K != L!=J #8
      return 0.0;
    }
  }else if(XI == XL){
    if(XJ == XK){
      // I=L != J==K #9
      return eleNum[XI]*(1-eleNum[XJ]);
    }else{
      // I=L != J!=K #10
      if(eleNum[XI]==0) return 0.0; // 
      else return -1.0*GreenFunc1_fsz2_real(rk,rj,u,t,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer); /* CkuAjt */
    }
  }else if(XJ == XK){
   if(XK == XL){
     // I != J=K=L  #11
     if(eleNum[XJ]==0) return 0.0; // 
     else return GreenFunc1_fsz2_real(ri,rj,s,t,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer); /* CisAjt */
   }else{
     // I != J=K !=L  #12
     if(eleNum[XJ]==1) return 0.0; // 
     else return GreenFunc1_fsz2_real(ri,rl,s,v,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer); /* CisAlv */
   } 
  }else if(XJ == XL){
    // I != J=L !=K  #13
    return 0.0; // 
  }else if(XK == XL){
    // I != J != K =L  #14
    if(eleNum[XK]==0) return 0.0; // 
    else return GreenFunc1_fsz2_real(ri,rj,s,t,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer); /* CisAjt */
  }
     
  // from here, no pair exists 
  if(eleNum[XI]==1 || eleNum[XJ]==0 || eleNum[XK]==1 || eleNum[XL]==0) return 0.0;


  mj = eleCfg[XJ]; // electron exists
  ml = eleCfg[XL]; // electron exists
  //msj = mj;// + s*Ne;
  //mtl = ml;// + t*Ne;

  /* hopping */
  eleIdx[ml] = rk; // ml : (rk,u)
  eleSpn[ml] = u;
  eleNum[XL] = 0;
  eleNum[XK] = 1;
  UpdateProjCnt_fsz(rl, rk, v,u, projCntNew, eleProjCnt, eleNum); // (rl,v) -> (rk,u) v!=u
  eleIdx[mj]  = ri; // mj : (ri,s)
  eleSpn[mj]  = s;
  eleNum[XJ] = 0;
  eleNum[XI] = 1;
  UpdateProjCnt_fsz(rj,ri, t,s, projCntNew, projCntNew, eleNum); // (rj,t) -> (ri,s) t!=s

  z = ProjRatio(projCntNew,eleProjCnt);

  /* calculate Pfaffian */
  CalculateNewPfMTwo_fsz_real(ml, u, mj, s, pfMNew, eleIdx,eleSpn, 0, NQPFull, bufV); // ml -> rk,u; mj -> ri,s
  z *= CalculateIP_real(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[ml]  = rl;  // ml : (rl,v)
  eleSpn[ml]  = v;
  eleNum[XL] = 1;
  eleNum[XK] = 0;
  eleIdx[mj]  = rj;  // mj : (rj,t)
  eleSpn[mj]  = t;
  eleNum[XJ] = 1;
  eleNum[XI] = 0;

  return z/ip;//TBC
}


