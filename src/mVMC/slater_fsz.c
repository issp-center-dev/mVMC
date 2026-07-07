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
 * slater elements
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

void UpdateSlaterElm_fsz();
void MakeSlaterElmBF_fsz(const int *eleNum, const int *eleProjBFCnt);
void SlaterElmDiff_fsz(double complex *srOptO, const double complex ip, int *eleIdx,int *eleSpn);

static inline int BFThetaCount_fsz(const int mu, const int centerSite,
                                   const int spin, const int rangeSlot,
                                   const int *eleProjBFCnt) {
  const int nSiteRange = Nsite*Nrange;
  const int channelOffset = (spin == 0) ? 0 : 4*nSiteRange;
  if(mu < 0 || mu >= 4) return 0;
  return eleProjBFCnt[channelOffset + mu*nSiteRange + centerSite*Nrange + rangeSlot];
}

static int BFHasActiveTheta_fsz(const int centerSite, const int spin,
                                const int *eleProjBFCnt) {
  int rangeSlot;
  /* Match the fcmp BackFlow eta gate: only the doublon-holon theta channel activates eta. */
  for(rangeSlot=0;rangeSlot<Nrange;rangeSlot++) {
    if(BFThetaCount_fsz(1, centerSite, spin, rangeSlot, eleProjBFCnt) != 0) {
      return 1;
    }
  }
  return 0;
}

static inline double complex SlaterOrbital_fsz(const int ri, const int si,
                                               const int rj, const int sj) {
  const int rsi = ri + si*Nsite;
  const int rsj = rj + sj*Nsite;
  return Slater[OrbitalIdx[rsi][rsj]] * (double)OrbitalSgn[rsi][rsj];
}

static int BFSubIndexFromMuRange_fsz(const int mu, const int centerSite,
                                     const int rangeSite) {
  int d, xtmp, idx;
  d = RangeIdx[centerSite][rangeSite];
  xtmp = 4*d + mu;
  idx = xtmp - 3 - d;
  if(xtmp%4 == 0) idx = -1;
  if(xtmp == 0) idx = 0;
  return idx;
}

static double complex SubSlaterElmBF_fsz(const int ri, const int si,
                                         const int rj, const int sj,
                                         const int *eleProjBFCnt) {
  int mu, nu;
  int k, l;
  int rk, rl;
  int nidx, midx, bfidx;
  int cntI, cntJ;
  double complex slt = 0.0 + 0.0*I;
  double eta;

  for(mu=0;mu<4;mu++) {
    for(nu=0;nu<4;nu++) {
      if(mu == 0 && nu == 0) continue;
      for(k=0;k<Nrange;k++) {
        rk = PosBF[ri][k];
        nidx = BFSubIndexFromMuRange_fsz(mu, ri, rk);
        if(nidx < 0) continue;
        cntI = BFThetaCount_fsz(mu, ri, si, k, eleProjBFCnt);
        if(cntI == 0) continue;

        for(l=0;l<Nrange;l++) {
          rl = PosBF[rj][l];
          midx = BFSubIndexFromMuRange_fsz(nu, rj, rl);
          if(midx < 0) continue;
          cntJ = BFThetaCount_fsz(nu, rj, sj, l, eleProjBFCnt);
          if(cntJ == 0) continue;

          bfidx = BFSubIdx[nidx][midx];
          slt += -ProjBF[bfidx] * (double)(cntI*cntJ)
                 * SlaterOrbital_fsz(rk, si, rl, sj);
        }
      }
    }
  }

  eta = (BFHasActiveTheta_fsz(ri, si, eleProjBFCnt) ||
         BFHasActiveTheta_fsz(rj, sj, eleProjBFCnt))
      ? creal(ProjBF[0]) : 1.0;
  slt += eta * SlaterOrbital_fsz(ri, si, rj, sj);

  return slt;
}

void UpdateSlaterElm_fsz() {
  int ri,ori,tri,sgni,rsi0,rsi1;
  int rj,orj,trj,sgnj,rsj0,rsj1;
  int qpidx,mpidx,optidx;
  int tri0,tri1,trj0,trj1; //fsz
  double complex slt_i0j0,slt_j0i0;
  double complex slt_i0j1,slt_j1i0;
  double complex slt_i1j0,slt_j0i1;
  double complex slt_i1j1,slt_j1i1;
  int *xqp, *xqpSgn, *xqpOpt, *xqpOptSgn;
  double complex *sltE,*sltE_i0,*sltE_i1;

  #pragma omp parallel for default(shared)        \
    private(qpidx,optidx,mpidx,                   \
            xqpOpt,xqpOptSgn,xqp,xqpSgn,sltE,     \
            ri,ori,tri,sgni,rsi0,rsi1,sltE_i0,sltE_i1,      \
            rj,orj,trj,sgnj,rsj0,rsj1,slt_i0j0,slt_j0i0,slt_i0j1,slt_j1i0,slt_i1j0,slt_j0i1,slt_i1j1,slt_j1i1,\
            tri0,tri1,trj0,trj1)
  #pragma loop noalias
  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    //printf("qpidx=%d \n",qpidx);
    // note: NQPFix = NSPGaussLeg * NMPTrans, NQPFull = NQPFix * NQPOptTrans(=1);
    // qpidx  = optidx*NQPFix+NSPGaussLeg*mpidx+spidx 
    // optidx will not be used
    optidx    = qpidx / NQPFix;              // optidx =0:  
    mpidx     = (qpidx%NQPFix) / NSPGaussLeg;// mpidx !=0: momentum projection
//  spidx     = qpidx % NSPGaussLeg;         // spidx  =0: spin projection

    xqpOpt    = QPOptTrans[optidx];          // xqpOpta   = 1
    xqpOptSgn = QPOptTransSgn[optidx];       // xqpOptSgn = 1
    xqp       = QPTrans[mpidx];              // xqp[*]=QPTrans[mpidx][*] 
    xqpSgn    = QPTransSgn[mpidx];           // xqpSgn[*]=QPTransSgn[mpidx][*] 
//    cs        = SPGLCosSin[spidx]; // spin projection is not used fsz
//    cc        = SPGLCosCos[spidx];
//    ss        = SPGLSinSin[spidx];
    
    sltE      = SlaterElm + qpidx*Nsite2*Nsite2;
    // for fsz: only for translational projection
    // spin projection is not implemented in ver.1.0
    for(ri=0;ri<Nsite;ri++) {
      ori     = xqpOpt[ri]; // ori = ri
      tri     = xqp[ori];   // ori -> tri momentum projection
      tri0    = tri;        // fsz up
      tri1    = tri0+Nsite; // fsz down
      sgni    = xqpSgn[ori]*xqpOptSgn[ri];
      rsi0    = ri;       // up
      rsi1    = ri+Nsite; // down
      sltE_i0 = sltE + rsi0*Nsite2; // sltE_i0[*]=sltE[rsi0*Nsite2][*]
      sltE_i1 = sltE + rsi1*Nsite2; // sltE_i1[*]=sltE[rsi1*Nsite2][*]
      
      for(rj=0;rj<Nsite;rj++) {
        //
        orj  = xqpOpt[rj];
        trj  = xqp[orj];
        //
        trj0 = trj;       // fsz up
        trj1 = trj+Nsite; // fsz down
        sgnj = xqpSgn[orj]*xqpOptSgn[rj];
        //
        rsj0 = rj;
        rsj1 = rj+Nsite;
        //
// for fsz        
        slt_i0j0        = Slater[ OrbitalIdx[tri0][trj0] ] * (double)(OrbitalSgn[tri0][trj0]*sgni*sgnj);
        slt_j0i0        = Slater[ OrbitalIdx[trj0][tri0] ] * (double)(OrbitalSgn[trj0][tri0]*sgni*sgnj);

        slt_i0j1        = Slater[ OrbitalIdx[tri0][trj1] ] * (double)(OrbitalSgn[tri0][trj1]*sgni*sgnj);
        slt_j1i0        = Slater[ OrbitalIdx[trj1][tri0] ] * (double)(OrbitalSgn[trj1][tri0]*sgni*sgnj);

        slt_i1j0        = Slater[ OrbitalIdx[tri1][trj0] ] * (double)(OrbitalSgn[tri1][trj0]*sgni*sgnj);
        slt_j0i1        = Slater[ OrbitalIdx[trj0][tri1] ] * (double)(OrbitalSgn[trj0][tri1]*sgni*sgnj);

        slt_i1j1        = Slater[ OrbitalIdx[tri1][trj1] ] * (double)(OrbitalSgn[tri1][trj1]*sgni*sgnj);
        slt_j1i1        = Slater[ OrbitalIdx[trj1][tri1] ] * (double)(OrbitalSgn[trj1][tri1]*sgni*sgnj);

        // F_{IJ}-F_{JI}
        sltE_i0[rsj0] =  slt_i0j0-slt_j0i0;   // up   - up
        sltE_i0[rsj1] =  slt_i0j1-slt_j1i0;   // up   - down
        sltE_i1[rsj0] =  slt_i1j0-slt_j0i1;   // down - up
        sltE_i1[rsj1] =  slt_i1j1-slt_j1i1;   // down - down 
      }// for rj 
    }// for ri
  }//for qpidx
  return;
}

void MakeSlaterElmBF_fsz(const int *eleNum, const int *eleProjBFCnt) {
  int ri,ori,tri,sgni,rsi;
  int rj,orj,trj,sgnj,rsj;
  int si,sj;
  int qpidx,mpidx,optidx;
  int *xqp, *xqpSgn, *xqpOpt, *xqpOptSgn;
  double complex slt_ij, slt_ji;
  double complex *sltE, *sltE_i;

  (void)eleNum;

  #pragma omp parallel for default(shared)        \
    private(qpidx,optidx,mpidx,                   \
            xqpOpt,xqpOptSgn,xqp,xqpSgn,sltE,     \
            ri,ori,tri,sgni,rsi,sltE_i,           \
            rj,orj,trj,sgnj,rsj,si,sj,slt_ij,slt_ji)
  #pragma loop noalias
  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    optidx    = qpidx / NQPFix;
    mpidx     = (qpidx%NQPFix) / NSPGaussLeg;

    xqpOpt    = QPOptTrans[optidx];
    xqpOptSgn = QPOptTransSgn[optidx];
    xqp       = QPTrans[mpidx];
    xqpSgn    = QPTransSgn[mpidx];

    sltE      = SlaterElmBF + qpidx*Nsite2*Nsite2;

    for(ri=0;ri<Nsite;ri++) {
      ori  = xqpOpt[ri];
      tri  = xqp[ori];
      sgni = xqpSgn[ori]*xqpOptSgn[ri];

      for(si=0;si<2;si++) {
        rsi = ri + si*Nsite;
        sltE_i = sltE + rsi*Nsite2;

        for(rj=0;rj<Nsite;rj++) {
          orj  = xqpOpt[rj];
          trj  = xqp[orj];
          sgnj = xqpSgn[orj]*xqpOptSgn[rj];

          for(sj=0;sj<2;sj++) {
            rsj = rj + sj*Nsite;
            slt_ij = SubSlaterElmBF_fsz(tri, si, trj, sj, eleProjBFCnt)
                     * (double)(sgni*sgnj);
            slt_ji = SubSlaterElmBF_fsz(trj, sj, tri, si, eleProjBFCnt)
                     * (double)(sgni*sgnj);
            sltE_i[rsj] = slt_ij - slt_ji;
          }
        }
      }
    }
  }
  return;
}

// Calculating Tr[Inv[M]*D_k(X)]
void SlaterElmDiff_fsz(double complex *srOptO, const double complex ip, int *eleIdx,int *eleSpn) {
  const int nBuf=NSlater*NQPFull;
  const int nsize = Nsize;       // Nsize=2*Ne
  const int nQPFull = NQPFull;   // note: NQPFix = NSPGaussLeg * NMPTrans, NQPFull = NQPFix * NQPOptTrans(=1);
  const int nMPTrans = NMPTrans; // number of translational operators
  const int nSlater = NSlater;   // number of slater
  const int nTrans = NMPTrans * NQPOptTrans; //usually NQPOptTrans=1

  const double complex invIP = 1.0/ip;
  int msi,msj,ri,rj,ori,orj,tri,trj,sgni,sgnj;
  int mpidx,orbidx,qpidx,optidx,i;
  int si,sj;//fsz
  double complex cs; // including Pf
  int *xqp,*xqpSgn,*xqpOpt,*xqpOptSgn;
  double complex *invM,*invM_i;

  int *orbitalIdx_i;
  int *transOrbIdx; /* transOrbIdx[mpidx][msi][msj] */
  int *transOrbSgn; /* transOrbSgn[mpidx][msi][msj] */
  int *tOrbIdx,*tOrbIdx_i;
  int *tOrbSgn,*tOrbSgn_i;
  double complex *buf, *buffer;
  double complex tmp;

  RequestWorkSpaceInt(2*nTrans*Nsize*Nsize);
  RequestWorkSpaceComplex(NQPFull*NSlater);

  transOrbIdx = GetWorkSpaceInt(nTrans*Nsize*Nsize); /* transOrbIdx[mpidx][msi][msj] */
  transOrbSgn = GetWorkSpaceInt(nTrans*Nsize*Nsize); /* transOrbSgn[mpidx][msi][msj] */
  buffer = GetWorkSpaceComplex(NQPFull*NSlater);

  for(i=0;i<nBuf;i++) buffer[i]=0.0;

  #pragma omp parallel for default(shared)                        \
    private(qpidx,optidx,mpidx,msi,msj,xqp,xqpSgn,xqpOpt,xqpOptSgn,\
            ri,ori,tri,sgni,rj,orj,trj,sgnj,si,sj,                       \
            tOrbIdx,tOrbIdx_i,tOrbSgn,tOrbSgn_i,orbitalIdx_i)
  #pragma loop noalias
// this loop for making transObrIdx,transOrbSgn
// transOrbIdx (NMPTrans * NQPOptTrans)*(nsize*nsize) nsize=2Ne
// transOrbSgn (NMPTrans * NQPOptTrans)*(nsize*nsize) nsize=2Ne
  for(qpidx=0;qpidx<nTrans;qpidx++) {//only for trans op. // nTran=nMPTrans*(NQPOptTrans=1)=nMPTrans
    optidx    = qpidx / nMPTrans;                // qpidx=mpidx+nMPTrans*optidx
    mpidx     = qpidx % nMPTrans; 
                                                 // QPOptTrans:     NQPOptTrans* Nsite matrix :  usually NQPOptTrans = 1
                                                 // QPOptTransSgn:  NQPOptTrans* Nsite matrix :  usually NQPOptTrans = 1 
    xqpOpt    = QPOptTrans[optidx];              // xqpOpt[]    = QPOptTrans[optidx][]    : usually optidx=0 and not be used
    xqpOptSgn = QPOptTransSgn[optidx];           // xqpOptSgn[] = QPOptTransSgn[optidx][] : usually optidx=0 and not be used
    xqp       = QPTrans[mpidx];                  // xqp    =  QPTrans[mpidx][]
    xqpSgn    = QPTransSgn[mpidx];               // xqpSgn =  QPTransSgn[mpidx][]
    tOrbIdx   = transOrbIdx + qpidx*nsize*nsize; // tOrbIdx[msi][msj]=OrbitalIdx[tri][trj]
    tOrbSgn   = transOrbSgn + qpidx*nsize*nsize; // tOrbSgn : sign of f_ij
    for(msi=0;msi<nsize;msi++) {                 // nsize=2*Ne including spin degrees of freedom
      ri           = eleIdx[msi];                //  ri  : postion where the msi-th electron exists  
      si           = eleSpn[msi];                //  si  : spin of the msi-th electron  :fsz
      ori          = xqpOpt[ri];                 // ori  : ri -> ori -> tri
      tri          = xqp[ori]+si*Nsite;          // tri  : fsz
      sgni         = xqpSgn[ori]*xqpOptSgn[ri];  // sgni :
      tOrbIdx_i    = tOrbIdx + msi*nsize;        // 
      tOrbSgn_i    = tOrbSgn + msi*nsize;        //
      orbitalIdx_i = OrbitalIdx[tri];            //
      for(msj=0;msj<nsize;msj++) { //nsize=2*Ne //
        rj             = eleIdx[msj];           //
        sj             = eleSpn[msj];        // sj spin of the msj-th electron  :fsz
        orj            = xqpOpt[rj];         
        trj            = xqp[orj]+sj*Nsite;       //fsz
        sgnj           = xqpSgn[orj]*xqpOptSgn[rj]; //xqpSgn
        tOrbIdx_i[msj] = orbitalIdx_i[trj];
        tOrbSgn_i[msj] = sgni*sgnj*OrbitalSgn[tri][trj];
      }
    }
  }
// calculating Tr(X^{-1}*dX/df_{msi,msj})=-2*alpha(sigma(msi),sigma(msj))(X^{-1})_{msi,msj}
  #pragma omp parallel for default(shared)        \
    private(qpidx,mpidx,cs,                   \
            tOrbIdx,tOrbSgn,invM,buf,msi,msj,             \
            tOrbIdx_i,tOrbSgn_i,invM_i,orbidx)
  #pragma loop noalias
  for(qpidx=0;qpidx<nQPFull;qpidx++) { // nQPFull = NQPFix * NQPOptTrans: usually = NSPGaussLeg * NMPTrans
    mpidx = qpidx / NSPGaussLeg;       // qpidx   = NSPGaussLeg*mpidx + spidx

    cs = PfM[qpidx];// * SPGLCosSin[spidx]; //  PfM

    tOrbIdx = transOrbIdx + mpidx*nsize*nsize; // tOrbIdx[msi][msj] = transOrbIdx[mpidx][msi][msj]
    tOrbSgn = transOrbSgn + mpidx*nsize*nsize; // tOrbSgn[msi][msj] = transOrbSgn[mpidx][msi][msj]
    invM    = InvM        + qpidx*Nsize*Nsize; // invM  M^-1 and PfM, InvM: NQPFull*(Nsize*Nsize)+PfM
    buf     = buffer      + qpidx*NSlater;     // buf   fij, buffer: NSlater*NQPFull

    #pragma loop norecurrence
    for(msi=0;msi<nsize;msi++) { //fsz
      tOrbIdx_i = tOrbIdx + msi*nsize;         // tOrbIdx_i[*] = tOrbIdx[msi][*]
      tOrbSgn_i = tOrbSgn + msi*nsize;         // tOrbSgn_i[*] = tOrbSgn[msi][*]
      invM_i    = invM + msi*nsize;            // invM[]       = invM[msi][]
      for(msj=0;msj<nsize;msj++) {             // 
        orbidx       = tOrbIdx_i[msj];         // 
        buf[orbidx] += -1.0*invM_i[msj]*cs*tOrbSgn_i[msj]; // invM[msi][msj]// 2*(1/2)=1
      }
    }
  }

  /* store SROptO[] */
  for(orbidx=0;orbidx<nSlater;orbidx++) {
    srOptO[2*orbidx]   = 0.0+0.0*I; // 0
    srOptO[2*orbidx+1] = 0.0+0.0*I; // 0
  }
  #pragma loop noalias
  for(qpidx=0;qpidx<nQPFull;qpidx++) {
    tmp = QPFullWeight[qpidx];
    buf = buffer + qpidx*nSlater;
    for(orbidx=0;orbidx<nSlater;orbidx++) {
      srOptO[2*orbidx]   += tmp * buf[orbidx];   //real      TBC
      srOptO[2*orbidx+1] += tmp * buf[orbidx]*I; //imaginary TBC
      //printf("Re DEBUG: tmp=%lf :orbidx=%d srOptO=%lf %lf invIP=%lf %lf \n",tmp,orbidx,creal(srOptO[2*orbidx]),cimag(srOptO[2*orbidx]),creal(invIP),cimag(invIP));
      //printf("Im DEBUG:orbidx=%d srOptO=%lf %lf invIP=%lf %lf \n",orbidx,creal(srOptO[2*orbidx+1]),cimag(srOptO[2*orbidx+1]),creal(invIP),cimag(invIP));
    }
  }
  for(orbidx=0;orbidx<nSlater;orbidx++) {
    srOptO[2*orbidx]   *= invIP;
    srOptO[2*orbidx+1] *= invIP;
  }

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceComplex();
  return;
}
