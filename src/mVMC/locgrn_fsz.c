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
 * by Satoshi Morita
 *-------------------------------------------------------------*/

double complex GreenFunc1_fsz(const int ri, const int rj, const int s, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double complex *buffer);

double complex GreenFunc1_fsz2(const int ri, const int rj, const int s,const int t, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double complex *buffer);

double complex GreenFunc2_fsz(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double complex  ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double complex *buffer);
/*
double complex GreenFuncN(const int n, int *rsi, int *rsj, const double complex  ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  double complex *buffer, int *bufferInt);
double complex calculateNewPfMN_child(const int qpidx, const int n, const int *msa, const int *rsa,
                              const int *eleIdx, double complex *buffer);
*/
/* Calculate 1-body Green function <CisAjs> */
/* buffer size = NQPFull */
double complex GreenFunc1_fsz(const int ri, const int rj, const int s, const double complex  ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, complex double *buffer) {
  double complex z;
  int mj,msj,rsi,rsj;
  double complex *pfMNew = buffer; /* NQPFull */

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
  CalculateNewPfM_fsz(mj, s, pfMNew, eleIdx,eleSpn, 0, NQPFull);//fsz
  z *= CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[msj] = rj;
  eleSpn[msj] = s; //fsz
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  return conj(z/ip);//TBC
}

/* Calculate 1-body Green function <CisAjt> */ // s!=t
/* buffer size = NQPFull */
double complex GreenFunc1_fsz2(const int ri, const int rj, const int s,const int t, const double complex  ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, complex double *buffer) {
  double complex z;
  int mj,rsi,rtj;//msj
  double complex *pfMNew = buffer; /* NQPFull */

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
  CalculateNewPfM_fsz(mj, s, pfMNew, eleIdx,eleSpn, 0, NQPFull);//fsz: note EleSpn[mj]=s
  z *= CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[mj] = rj;
  eleSpn[mj] = t; //fsz
  eleNum[rtj] = 1;
  eleNum[rsi] = 0;

  return conj(z/ip);//TBC
}


/* Calculate 2-body Green function <psi|CisAjsCktAlt|x>/<psi|x> */
/* buffer size = NQPFull+2*Nsize */
double complex GreenFunc2_fsz_(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double complex *buffer,
                  int *lazy_info, int *lazy_rsi, int *lazy_msj) {
  double complex z;
  int mj,msj,ml,mtl;
  int rsi,rsj,rtk,rtl;
  double complex *pfMNew = buffer; /* [NQPFull] */
  double complex *bufV   = buffer+NQPFull; /* 2*Nsize */

  rsi = ri + s*Nsite;
  rsj = rj + s*Nsite;
  rtk = rk + t*Nsite;
  rtl = rl + t*Nsite;

  if(s==t) {
    if(rk==rl) { /* CisAjsNks */
      if(eleNum[rtk]==0) return 0.0;
      else return GreenFunc1_fsz(ri,rj,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CisAjs */
    }else if(rj==rl) {
      return 0.0; /* CisAjsCksAjs (j!=k) */
    }else if(ri==rl) { /* AjsCksNis */
      if(eleNum[rsi]==0) return 0.0;
      else if(rj==rk) return 1.0-eleNum[rsj];
      else return -GreenFunc1_fsz(rk,rj,s,ip,eleIdx,eleCfg,eleNum,
                              eleProjCnt,eleSpn,projCntNew,buffer); /* -CksAjs */
    }else if(rj==rk) { /* CisAls(1-Njs) */
      if(eleNum[rsj]==1) return 0.0;
      else if(ri==rl) return eleNum[rsi];
      else return GreenFunc1_fsz(ri,rl,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CisAls */
    }else if(ri==rk) {
      return 0.0; /* CisAjsCisAls (i!=j) */
    }else if(ri==rj) { /* NisCksAls (i!=k,l) */
      if(eleNum[rsi]==0) return 0.0;
      else return GreenFunc1_fsz(rk,rl,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CksAls */
    }
  }else{
    if(rk==rl) { /* CisAjsNkt */
      if(eleNum[rtk]==0) return 0.0;
      else if(ri==rj) return eleNum[rsi];
      else return GreenFunc1_fsz(ri,rj,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,eleSpn,projCntNew,buffer); /* CisAjs */
    }else if(ri==rj) { /* NisCktAlt */
      if(eleNum[rsi]==0) return 0.0;
      else return GreenFunc1_fsz(rk,rl,t,ip,eleIdx,eleCfg,eleNum,
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
  if (!lazy_info) {
    CalculateNewPfMTwo_fsz(ml, t, mj, s, pfMNew, eleIdx,eleSpn, 0, NQPFull, bufV);
    z *= CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);
  } else {
    // Save invocation hook instead of calling.
    lazy_msj[0] = msj; //< to edit to ri, revert to rj
    lazy_msj[1] = mtl; //< to edit to rk, revert to rl
    lazy_rsi[0] = rsi;
    lazy_rsi[1] = rtk;
    lazy_info[0] = 1; //< "is delayed" flag..
  }

  /* revert hopping */
  eleIdx[mtl] = rl;
  eleSpn[mtl] = t;
  eleNum[rtl] = 1;
  eleNum[rtk] = 0;
  eleIdx[msj] = rj;
  eleSpn[msj] = s;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  return conj(z/ip);//TBC
}

double complex GreenFunc2_fsz(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double complex *buffer) {
  return GreenFunc2_fsz_(ri, rj, rk, rl, s, t, ip, eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn,
                         projCntNew, buffer, 0, 0, 0);
}

/* Calculate 2-body Green function <psi|CisAjtCkuAlv|x>/<psi|x> */
/* buffer size = NQPFull+2*Nsize */
double complex GreenFunc2_fsz2_(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t,const int u,const int v, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double complex *buffer,
                  int *lazy_info, int *lazy_rsi, int *lazy_msj) {
  double complex z;
  //int mj,msj,ml,mtl;
  int mj,ml; // mj: (rj,t) -> (ri,s) ml:(rl,v) -> (rk,u)
  //int rsi,rtj,ruk,rvl;
  int XI,XJ,XK,XL; //fsz
  double complex *pfMNew = buffer; /* [NQPFull] */
  double complex *bufV   = buffer+NQPFull; /* 2*Nsize */

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
        else return GreenFunc1_fsz2(rk,rl,u,v,ip,eleIdx,eleCfg,eleNum,
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
      else return GreenFunc1_fsz2(rk,rl,u,v,ip,eleIdx,eleCfg,eleNum,
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
      else return -1.0*GreenFunc1_fsz2(rk,rj,u,t,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer); /* CkuAjt */
    }
  }else if(XJ == XK){
   if(XK == XL){
     // I != J=K=L  #11
     if(eleNum[XJ]==0) return 0.0; // 
     else return GreenFunc1_fsz2(ri,rj,s,t,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer); /* CisAjt */
   }else{
     // I != J=K !=L  #12
     if(eleNum[XJ]==1) return 0.0; // 
     else return GreenFunc1_fsz2(ri,rl,s,v,ip,eleIdx,eleCfg,eleNum,
                            eleProjCnt,eleSpn,projCntNew,buffer); /* CisAlv */
   } 
  }else if(XJ == XL){
    // I != J=L !=K  #13
    return 0.0; // 
  }else if(XK == XL){
    // I != J != K =L  #14
    if(eleNum[XK]==0) return 0.0; // 
    else return GreenFunc1_fsz2(ri,rj,s,t,ip,eleIdx,eleCfg,eleNum,
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
  if (!lazy_info) {
    CalculateNewPfMTwo_fsz(ml, u, mj, s, pfMNew, eleIdx,eleSpn, 0, NQPFull, bufV); // ml -> rk,u; mj -> ri,s
    z *= CalculateIP_fcmp(pfMNew, 0, NQPFull, MPI_COMM_SELF);
  } else {
    lazy_msj[0] = mj;
    lazy_msj[1] = ml;
    lazy_rsi[0] = XI;
    lazy_rsi[1] = XK;
    lazy_info[0] = 1;
  }

  /* revert hopping */
  eleIdx[ml]  = rl;  // ml : (rl,v)
  eleSpn[ml]  = v;
  eleNum[XL] = 1;
  eleNum[XK] = 0;
  eleIdx[mj]  = rj;  // mj : (rj,t)
  eleSpn[mj]  = t;
  eleNum[XJ] = 1;
  eleNum[XI] = 0;

  return conj(z/ip);//TBC
}


double complex GreenFunc2_fsz2(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t,const int u,const int v, const double complex ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,int *eleSpn,
                  int *projCntNew, double complex *buffer) {
  return GreenFunc2_fsz2_(ri, rj, rk, rl, s, t, u, v, ip,
                          eleIdx, eleCfg, eleNum, eleProjCnt, eleSpn, projCntNew, buffer, 0, 0, 0);
}

// ignore GreenFuncN: to be added

/* Calculate n-body Green function */
/* <phi| c1 a1 c2 a2 ... cn an |x> */
/* c_k = c_rsi[k], a_k = a_rsj[k] */
/* rsi is the array of indices of creation operators */
/* rsj is the array of indices of annihilation operators */
/* buffer size    = NQPFull + n*Nsize */
/* bufferInt size = NProj */

//double complex GreenFuncN(const int n, int *rsi, int *rsj, const double complex ip,
//                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
//                  double complex *buffer, int *bufferInt){
//  int ri,rj,rk,rl,si,sj,sk,mj;
//  int k,l,m,rsk;
//  double complex z,x;
//  int qpidx;
//
//  int *projCntNew = bufferInt; /* [NProj] */
//
//  int msj[n];
//  double complex *pfMNew = buffer; /* [NQPFull] */
//  double complex *bufV = buffer+NQPFull; /* [n*Nsize] */
//
//  for(k=0;k<n;k++) {
//    si = rsi[k]/Nsite;
//    sj = rsj[k]/Nsite;
//    if(si!=sj) return 0;
//  }
//
//  if(n<=0) return 0;
//  else if(n==1) {
//    ri = rsi[0]%Nsite;
//    rj = rsj[0]%Nsite;
//    si = rsi[0]/Nsite;
//    return GreenFunc1(ri,rj,si,ip,eleIdx,eleCfg,eleNum,eleProjCnt,projCntNew,buffer);
//  } else if(n==2) {
//    ri = rsi[0]%Nsite;
//    rj = rsj[0]%Nsite;
//    si = rsi[0]/Nsite;
//    rk = rsi[1]%Nsite;
//    rl = rsj[1]%Nsite;
//    sk = rsi[1]/Nsite;
//    return GreenFunc2(ri,rj,rk,rl,si,sk,ip,eleIdx,eleCfg,eleNum,eleProjCnt,projCntNew,buffer);
//  }
//
//  /* reduction */
//  for(k=n-1;k>=0;k--) {
//    /* ** check for an annihilation operator at rsj[k] ** */
//    rsk = rsj[k];
//    for(l=k+1;l<n;l++) {
//      /* rsj[k] == rsi[l] */
//      if(rsk==rsi[l]) {
//        rsj[k] = rsj[l];
//        for(m=l;m<n-1;m++) { /* shift */
//          rsi[m] = rsi[m+1];
//          rsj[m] = rsj[m+1];
//        }
//        return GreenFuncN(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,eleProjCnt,buffer,bufferInt);
//      }
//      /* rsj[k] == rsj[l] */
//      if(rsk==rsj[l]) return 0;
//    }
//    /* check electron number */
//    if(eleNum[rsk]==0) return 0;
//
//    /* ** check for a creation operator at rsi[k] ** */
//    rsk = rsi[k];
//    /* rsi[k] == rsj[k] */
//    if(rsk==rsj[k]) {
//      for(m=k;m<n-1;m++) { /* shift */
//        rsi[m] = rsi[m+1];
//        rsj[m] = rsj[m+1];
//      }
//      return GreenFuncN(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,eleProjCnt,buffer,bufferInt);
//    }
//    for(l=k+1;l<n;l++) {
//      /* rsi[k] == rsi[l] */
//      if(rsk==rsi[l]) return 0;
//      /* rsi[k] == rsj[l] (k<l) */
//      if(rsk==rsj[l]) {
//        rsi[k] = rsi[l];
//        for(m=l;m<n-1;m++) { /* shift */
//          rsi[m] = rsi[m+1];
//          rsj[m] = rsj[m+1];
//        }
//        return (-1.0)*GreenFuncN(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,eleProjCnt,buffer,bufferInt);
//      }
//    }
//    /* check electron number */
//    if(eleNum[rsk]==1) return 0;
//  }
//
//  /* hopping */
//  #pragma loop noalias
//  for(k=0;k<n;k++) {
//    ri = rsi[k]%Nsite;
//    rj = rsj[k]%Nsite;
//    sj = rsj[k]/Nsite;
//    mj = eleCfg[rsj[k]];
//    msj[k] = mj + sj*Ne;
//    
//    eleIdx[msj[k]] = ri;
//    eleNum[rsj[k]] = 0;
//    eleNum[rsi[k]] = 1;
//  }
//
//  MakeProjCnt(projCntNew,eleNum);
//  x = ProjRatio(projCntNew,eleProjCnt);
//
//  /* calculateNewPfM */
//  for(qpidx=0;qpidx<NQPFull;qpidx++) {
//    pfMNew[qpidx] = calculateNewPfMN_child(qpidx,n,msj,rsj,eleIdx,bufV);
//  }
//  z = CalculateIP(pfMNew, 0, NQPFull, MPI_COMM_SELF);
//
//  /* revert hoppint */
//  #pragma loop noalias
//  for(k=0;k<n;k++) {
//    rj = rsj[k]%Nsite;    
//    eleIdx[msj[k]] = rj;
//    eleNum[rsj[k]] = 1;
//    eleNum[rsi[k]] = 0;
//  }
//
//  return x*z/ip;
//}
//
///* msa[k]-th electron hops from rsa[k] to eleIdx[msa[k]] */
///* buffer size = n*Nsize */
//double complex calculateNewPfMN_child(const int qpidx, const int n, const int *msa, const int *rsa,
//                              const int *eleIdx, double complex *buffer) {
//  const int nsize = Nsize;
//  const int n2 = 2*n;
//  const double complex *sltE;
//  const double complex *sltE_k;
//  const double complex *invM;
//  const double complex *invM_i, *invM_k, *invM_l;
//
//  double complex *vec; /* vec[n][nsize] */
//  double complex *vec_k, *vec_l;
//  double complex mat[n2*n2]; /* mat[n2][n2] */
//  double complex *mat_k;
//  double sgn;
//
//  int rsi,rsk,msi,msj,k,l;
//  double complex val,tmp;
//
//  /* for DSKPFA */
//  char uplo='U', mthd='P';
//  int nn,lda,info=0;
//  double complex pfaff;
//  int iwork[n2];
//  double complex work[n2*n2]; /* [n2][n2] */
//  int lwork = n2*n2;
//  nn=lda=n2;
//
//  sltE = SlaterElm + qpidx*Nsite2*Nsite2;
//  invM = InvM + qpidx*Nsize*Nsize;
//
//  vec = buffer; /* n*nsize */
//
//  #pragma loop noalias
//  for(k=0;k<n;k++) {
//    rsk = eleIdx[msa[k]] + (msa[k]/Ne)*Nsite;
//    sltE_k = sltE + rsk*Nsite2;
//    vec_k = vec + k*nsize;
//    #pragma loop norecurrence
//    for(msi=0;msi<nsize;msi++) {
//      rsi = eleIdx[msi] + (msi/Ne)*Nsite;
//      vec_k[msi] = sltE_k[rsi];
//    }
//  }
//
//  /* X_kl */
//  for(k=0;k<n;k++) {
//    mat_k = mat + n2*k;
//    vec_k = vec + k*nsize;
//    for(l=k+1;l<n;l++) {
//      vec_l = vec + l*nsize;
//      val = 0.0;
//      for(msi=0;msi<nsize;msi++) {
//        invM_i = invM + msi*nsize;
//        tmp = 0.0;
//        for(msj=0;msj<nsize;msj++) {
//          tmp += invM_i[msj] * vec_l[msj];
//        }
//        val += tmp * vec_k[msi];
//      }
//      mat_k[l] = val + vec_k[msa[l]];
//    }
//  }
//
//  /* Y_kl */
//  for(k=0;k<n;k++) {
//    mat_k = mat + n2*k + n;
//    vec_k = vec + k*nsize;
//    for(l=0;l<n;l++) {
//      invM_l = invM + msa[l]*nsize;
//      val = 0.0;
//      for(msi=0;msi<nsize;msi++) {
//        val += vec_k[msi] * invM_l[msi];
//      }
//      mat_k[l] = val;
//    }
//  }
//  
//  /* Z_kl */
//  for(k=0;k<n;k++) {
//    mat_k = mat + n2*(k+n) + n;
//    invM_k = invM + msa[k]*nsize;
//    for(l=k+1;l<n;l++) {
//      mat_k[l] = invM_k[msa[l]];
//    }
//  }
//
//  #pragma loop noalias
//  for(k=0;k<n2;k++) {
//    #pragma loop norecurrence
//    for(l=0;l<k;l++) {
//      mat[n2*k + l] = -mat[n2*l + k]; /* transpose */
//    }
//    mat[n2*k + k] = 0.0; /* diagonal elements */
//  }
//
//  /* calculate Pf M */
//  // Coverage-0 and commented out,
//  // but should be useful sometimes.
//  //M_DSKPFA(&uplo, &mthd, &nn, mat, &lda, &pfaff, iwork, work, &lwork, &info);
//  info = 1; // Skipping inverse.
//  M_ZSKPFA(&uplo, &mthd, &n, bufM, &lda, &pfaff, iwork, work, &lwork/*, rwork*/, &info);
//
//  sgn = ( (n*(n-1)/2)%2==0 ) ? 1.0 : -1.0;
//
//  return sgn * pfaff * PfM[qpidx];
//}
