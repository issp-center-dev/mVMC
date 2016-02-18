/*-------------------------------------------------------------
 * Variational Monte Carlo
 * local Green Functions
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

double GreenFunc1(const int ri, const int rj, const int s, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  int *projCntNew, double *buffer);
double GreenFunc2(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  int *projCntNew, double *buffer);
double GreenFuncN(const int n, int *rsi, int *rsj, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  double *buffer, int *bufferInt);
double calculateNewPfMN_child(const int qpidx, const int n, const int *msa, const int *rsa,
                              const int *eleIdx, double *buffer);

/* Calculate 1-body Green function <CisAjs> */
/* buffer size = NQPFull */
double GreenFunc1(const int ri, const int rj, const int s, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
                  int *projCntNew, double *buffer) {
  double z;
  int mj,msj,rsi,rsj;
  double *pfMNew = buffer; /* NQPFull */

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
  CalculateNewPfM(mj, s, pfMNew, eleIdx, 0, NQPFull);
  z *= CalculateIP(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  return z/ip;
}

/* Calculate 2-body Green function <psi|CisAjsCktAlt|x>/<psi|x> */
/* buffer size = NQPFull+2*Nsize */
double GreenFunc2(const int ri, const int rj, const int rk, const int rl,
                  const int s, const int t, const double ip,
                  int *eleIdx, const int *eleCfg, int *eleNum, const int *eleProjCnt,
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
      else return GreenFunc1(ri,rj,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,projCntNew,buffer); /* CisAjs */
    }else if(rj==rl) {
      return 0.0; /* CisAjsCksAjs (j!=k) */
    }else if(ri==rl) { /* AjsCksNis */
      if(eleNum[rsi]==0) return 0.0;
      else if(rj==rk) return 1.0-eleNum[rsj];
      else return -GreenFunc1(rk,rj,s,ip,eleIdx,eleCfg,eleNum,
                              eleProjCnt,projCntNew,buffer); /* -CksAjs */
    }else if(rj==rk) { /* CisAls(1-Njs) */
      if(eleNum[rsj]==1) return 0.0;
      else if(ri==rl) return eleNum[rsi];
      else return GreenFunc1(ri,rl,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,projCntNew,buffer); /* CisAls */
    }else if(ri==rk) {
      return 0.0; /* CisAjsCisAls (i!=j) */
    }else if(ri==rj) { /* NisCksAls (i!=k,l) */
      if(eleNum[rsi]==0) return 0.0;
      else return GreenFunc1(rk,rl,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,projCntNew,buffer); /* CksAls */
    }
  } if(s!=t) {
    if(rk==rl) { /* CisAjsNkt */
      if(eleNum[rtk]==0) return 0.0;
      else if(ri==rj) return eleNum[rsi];
      else return GreenFunc1(ri,rj,s,ip,eleIdx,eleCfg,eleNum,
                             eleProjCnt,projCntNew,buffer); /* CisAjs */
    }else if(ri==rj) { /* NisCktAlt */
      if(eleNum[rsi]==0) return 0.0;
      else return GreenFunc1(rk,rl,t,ip,eleIdx,eleCfg,eleNum,
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
  CalculateNewPfMTwo(ml, t, mj, s, pfMNew, eleIdx, 0, NQPFull, bufV);
  z *= CalculateIP(pfMNew, 0, NQPFull, MPI_COMM_SELF);

  /* revert hopping */
  eleIdx[mtl] = rl;
  eleNum[rtl] = 1;
  eleNum[rtk] = 0;
  eleIdx[msj] = rj;
  eleNum[rsj] = 1;
  eleNum[rsi] = 0;

  return z/ip;
}

/* Calculate n-body Green function */
/* <phi| c1 a1 c2 a2 ... cn an |x> */
/* c_k = c_rsi[k], a_k = a_rsj[k] */
/* rsi is the array of indices of creation operators */
/* rsj is the array of indices of annihilation operators */
/* buffer size    = NQPFull + n*Nsize */
/* bufferInt size = NProj */
double GreenFuncN(const int n, int *rsi, int *rsj, const double ip,
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
    return GreenFunc1(ri,rj,si,ip,eleIdx,eleCfg,eleNum,eleProjCnt,projCntNew,buffer);
  } else if(n==2) {
    ri = rsi[0]%Nsite;
    rj = rsj[0]%Nsite;
    si = rsi[0]/Nsite;
    rk = rsi[1]%Nsite;
    rl = rsj[1]%Nsite;
    sk = rsi[1]/Nsite;
    return GreenFunc2(ri,rj,rk,rl,si,sk,ip,eleIdx,eleCfg,eleNum,eleProjCnt,projCntNew,buffer);
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
        return GreenFuncN(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,eleProjCnt,buffer,bufferInt);
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
      return GreenFuncN(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,eleProjCnt,buffer,bufferInt);
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
        return (-1.0)*GreenFuncN(n-1,rsi,rsj,ip,eleIdx,eleCfg,eleNum,eleProjCnt,buffer,bufferInt);
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
    pfMNew[qpidx] = calculateNewPfMN_child(qpidx,n,msj,rsj,eleIdx,bufV);
  }
  z = CalculateIP(pfMNew, 0, NQPFull, MPI_COMM_SELF);

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

/* msa[k]-th electron hops from rsa[k] to eleIdx[msa[k]] */
/* buffer size = n*Nsize */
double calculateNewPfMN_child(const int qpidx, const int n, const int *msa, const int *rsa,
                              const int *eleIdx, double *buffer) {
  const int nsize = Nsize;
  const int n2 = 2*n;
  const double *sltE;
  const double *sltE_k;
  const double *invM;
  const double *invM_i, *invM_k, *invM_l;

  double *vec; /* vec[n][nsize] */
  double *vec_k, *vec_l;
  double mat[n2*n2]; /* mat[n2][n2] */
  double *mat_k;
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

  sltE = SlaterElm + qpidx*Nsite2*Nsite2;
  invM = InvM + qpidx*Nsize*Nsize;

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

  return sgn * pfaff * PfM[qpidx];
}
