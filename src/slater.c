/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Mitsuaki Kawamura, Takeo Kato, Masatoshi Imada.

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

void UpdateSlaterElm_fcmp();
void SlaterElmDiff_fcmp(double complex *srOptO, const double complex ip, int *eleIdx);

void UpdateSlaterElm_fcmp() {
  int ri,ori,tri,sgni,rsi0,rsi1;
  int rj,orj,trj,sgnj,rsj0,rsj1;
  int qpidx,mpidx,spidx,optidx;
  double cs,cc,ss;
  double complex slt_ij,slt_ji;
  int *xqp, *xqpSgn, *xqpOpt, *xqpOptSgn;
  double complex *sltE,*sltE_i0,*sltE_i1;

  #pragma omp parallel for default(shared)        \
    private(qpidx,optidx,mpidx,spidx,                      \
            xqpOpt,xqpOptSgn,xqp,xqpSgn,cs,cc,ss,sltE,     \
            ri,ori,tri,sgni,rsi0,rsi1,sltE_i0,sltE_i1,      \
            rj,orj,trj,sgnj,rsj0,rsj1,slt_ij,slt_ji)
  #pragma loop noalias
  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    // qpidx  = optidx*NQPFix+NSPGaussLeg*mpidx+spidx 
    // optidx will not be used
    optidx    = qpidx / NQPFix;   
    mpidx     = (qpidx%NQPFix) / NSPGaussLeg;
    spidx     = qpidx % NSPGaussLeg;

    xqpOpt    = QPOptTrans[optidx];
    xqpOptSgn = QPOptTransSgn[optidx];
    xqp       = QPTrans[mpidx];
    xqpSgn    = QPTransSgn[mpidx];
    cs        = SPGLCosSin[spidx];
    cc        = SPGLCosCos[spidx];
    ss        = SPGLSinSin[spidx];
    
    sltE      = SlaterElm + qpidx*Nsite2*Nsite2;
    
    for(ri=0;ri<Nsite;ri++) {
      ori     = xqpOpt[ri];
      tri     = xqp[ori];
      sgni    = xqpSgn[ori]*xqpOptSgn[ri];
      rsi0    = ri;
      rsi1    = ri+Nsite;
      sltE_i0 = sltE + rsi0*Nsite2;
      sltE_i1 = sltE + rsi1*Nsite2;
      
      for(rj=0;rj<Nsite;rj++) {
        orj = xqpOpt[rj];
        trj = xqp[orj];
        sgnj = xqpSgn[orj]*xqpOptSgn[rj];
        rsj0 = rj;
        rsj1 = rj+Nsite;
        
        slt_ij = Slater[ OrbitalIdx[tri][trj] ] * (double)(OrbitalSgn[tri][trj]*sgni*sgnj);
        slt_ji = Slater[ OrbitalIdx[trj][tri] ] * (double)(OrbitalSgn[trj][tri]*sgni*sgnj);
        
        sltE_i0[rsj0] = -(slt_ij - slt_ji)*cs;   // up   - up
        sltE_i0[rsj1] =   slt_ij*cc + slt_ji*ss; // up   - down
        sltE_i1[rsj0] = -slt_ij*ss - slt_ji*cc;  // down - up
        sltE_i1[rsj1] =  (slt_ij - slt_ji)*cs;   // down - down 
      }
    }
  }

  return;
}

// Calculating Tr[Inv[M]*D_k(X)]
void SlaterElmDiff_fcmp(double complex *srOptO, const double complex ip, int *eleIdx) {
  const int nBuf=NSlater*NQPFull;
  const int nsize = Nsize;
  const int ne = Ne;
  const int nQPFull = NQPFull;
  const int nMPTrans = NMPTrans; // number of translational operators
  const int nSlater = NSlater;
  const int nTrans = NMPTrans * NQPOptTrans; //usually NQPOptTrans=1

  const double complex invIP = 1.0/ip;
  int msi,msj,ri,rj,ori,orj,tri,trj,sgni,sgnj;
  int mpidx,spidx,orbidx,qpidx,optidx,i;
  double complex cs,cc,ss; // including Pf
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
            ri,ori,tri,sgni,rj,orj,trj,sgnj,                       \
            tOrbIdx,tOrbIdx_i,tOrbSgn,tOrbSgn_i,orbitalIdx_i)
  #pragma loop noalias
  for(qpidx=0;qpidx<nTrans;qpidx++) {            // nTran=nMPTrans*NQPOptTrans: usually nMPTrans
    optidx    = qpidx / nMPTrans;                // qpidx=mpidx+nMPTrans*optidx
    mpidx     = qpidx % nMPTrans; 
                                                 // QPOptTrans:     NQPOptTrans* Nsite matrix :  usually NQPOptTrans = 1
                                                 // QPOptTransSgn:  NQPOptTrans* Nsite matrix :  usually NQPOptTrans = 1 
    xqpOpt    = QPOptTrans[optidx];              // xqpOpt[]    = QPOptTrans[optidx][]    : usually optidx=0 and not be used
    xqpOptSgn = QPOptTransSgn[optidx];           // xqpOptSgn[] = QPOptTransSgn[optidx][] : usually optidx=0 and not be used
    xqp       = QPTrans[mpidx];                  // xqp    =  QPTrans[mpidx][]
    xqpSgn    = QPTransSgn[mpidx];               // xqpSgn =  QPTransSgn[mpidx][]
    tOrbIdx   = transOrbIdx + qpidx*nsize*nsize; // tOrbIdx : f_ij
    tOrbSgn   = transOrbSgn + qpidx*nsize*nsize; // tOrbSgn : sign of f_ij
    for(msi=0;msi<nsize;msi++) {                 // nsize=2*Ne
      ri           = eleIdx[msi];                //  ri  : postion where the msi-th electron exists  
      ori          = xqpOpt[ri];                 // ori  : 
      tri          = xqp[ori];                   // tri  : 
      sgni         = xqpSgn[ori]*xqpOptSgn[ri];  // sgni :
      tOrbIdx_i    = tOrbIdx + msi*nsize;        // 
      tOrbSgn_i    = tOrbSgn + msi*nsize;        //
      orbitalIdx_i = OrbitalIdx[tri];            //
      for(msj=0;msj<nsize;msj++) { //nsize=2*Ne //
        rj             = eleIdx[msj];           //
        orj            = xqpOpt[rj];         
        trj            = xqp[orj];
        sgnj           = xqpSgn[orj]*xqpOptSgn[rj];
        tOrbIdx_i[msj] = orbitalIdx_i[trj];
        tOrbSgn_i[msj] = sgni*sgnj*OrbitalSgn[tri][trj];
      }
    }
  }
// calculating Tr(X^{-1}*dX/df_{msi,msj})=-2*alpha(sigma(msi),sigma(msj))(X^{-1})_{msi,msj}
  #pragma omp parallel for default(shared)        \
    private(qpidx,mpidx,spidx,cs,cc,ss,                   \
            tOrbIdx,tOrbSgn,invM,buf,msi,msj,             \
            tOrbIdx_i,tOrbSgn_i,invM_i,orbidx)
  #pragma loop noalias
  for(qpidx=0;qpidx<nQPFull;qpidx++) { // nQPFull = NQPFix * NQPOptTrans: usually = NSPGaussLeg * NMPTrans
    mpidx = qpidx / NSPGaussLeg;       // qpidx   = NSPGaussLeg*mpidx + spidx
    spidx = qpidx % NSPGaussLeg;

    cs = PfM[qpidx] * SPGLCosSin[spidx]; // spin rotation + PfM
    cc = PfM[qpidx] * SPGLCosCos[spidx]; // spin rotation + PfM
    ss = PfM[qpidx] * SPGLSinSin[spidx]; // spin rotation + PfM

    tOrbIdx = transOrbIdx + mpidx*nsize*nsize; // tOrbIdx[msi][msj] = transOrbIdx[mpidx][msi][msj]
    tOrbSgn = transOrbSgn + mpidx*nsize*nsize; // tOrbSgn[msi][msj] = transOrbSgn[mpidx][msi][msj]
    invM    = InvM        + qpidx*Nsize*Nsize; // invM  M^-1 and PfM, InvM: NQPFull*(Nsize*Nsize)+PfM
    buf     = buffer      + qpidx*NSlater;     // buf   fij, buffer: NSlater*NQPFull

    #pragma loop norecurrence
    for(msi=0;msi<ne;msi++) {
      tOrbIdx_i = tOrbIdx + msi*nsize;         // tOrbIdx_i[] = tOrbIdx[msi][msj]
      tOrbSgn_i = tOrbSgn + msi*nsize;         // tOrbSgn_i[] = tOrbSgn[msi][msj]
      invM_i    = invM + msi*nsize;            // invM[]      = invM[msi][]
      for(msj=0;msj<ne;msj++) {                // up-up
        /* si=0 sj=0*/
        orbidx       = tOrbIdx_i[msj];         // 
        buf[orbidx] += invM_i[msj]*cs*tOrbSgn_i[msj]; // invM[msi][msj]
      }
      for(msj=ne;msj<nsize;msj++) {            // up-down
        /* si=0 sj=1*/
        orbidx       = tOrbIdx_i[msj];
        buf[orbidx] -= invM_i[msj]*cc*tOrbSgn_i[msj];
      }
    }
    #pragma loop norecurrence
    for(msi=ne;msi<nsize;msi++) { 
      tOrbIdx_i = tOrbIdx + msi*nsize;
      tOrbSgn_i = tOrbSgn + msi*nsize;
      invM_i = invM + msi*nsize;
      for(msj=0;msj<ne;msj++) {    // down-up
        /* si=1 sj=0*/
        orbidx = tOrbIdx_i[msj];
        buf[orbidx] += invM_i[msj]*ss*tOrbSgn_i[msj];
      }
      for(msj=ne;msj<nsize;msj++) {// down-down
        /* si=1 sj=1*/
        orbidx = tOrbIdx_i[msj];
        buf[orbidx] -= invM_i[msj]*cs*tOrbSgn_i[msj];
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
