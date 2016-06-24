/*-------------------------------------------------------------
 * Variational Monte Carlo
 * slater elements
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

void UpdateSlaterElm();
void SlaterElmDiff(double *srOptO, const double ip, int *eleIdx);

void UpdateSlaterElm() {
  int ri,ori,tri,sgni,rsi0,rsi1;
  int rj,orj,trj,sgnj,rsj0,rsj1;
  int qpidx,mpidx,spidx,optidx;
  double cs,cc,ss;
  double slt_ij,slt_ji;
  int *xqp, *xqpSgn, *xqpOpt, *xqpOptSgn;
  double *sltE,*sltE_i0,*sltE_i1;

  #pragma omp parallel for default(shared)        \
    private(qpidx,optidx,mpidx,spidx,                      \
            xqpOpt,xqpOptSgn,xqp,xqpSgn,cs,cc,ss,sltE,     \
            ri,ori,tri,sgni,rsi0,rsi1,sltE_i0,sltE_i1,      \
            rj,orj,trj,sgnj,rsj0,rsj1,slt_ij,slt_ji)
  #pragma loop noalias
  for(qpidx=0;qpidx<NQPFull;qpidx++) {
    optidx = qpidx / NQPFix;
    mpidx = (qpidx%NQPFix) / NSPGaussLeg;
    spidx = qpidx % NSPGaussLeg;

    xqpOpt = QPOptTrans[optidx];
    xqpOptSgn = QPOptTransSgn[optidx];
    xqp = QPTrans[mpidx];
    xqpSgn = QPTransSgn[mpidx];
    cs = SPGLCosSin[spidx];
    cc = SPGLCosCos[spidx];
    ss = SPGLSinSin[spidx];
    
    sltE = SlaterElm + qpidx*Nsite2*Nsite2;
    
    for(ri=0;ri<Nsite;ri++) {
      ori = xqpOpt[ri];
      tri = xqp[ori];
      sgni = xqpSgn[ori]*xqpOptSgn[ri];
      rsi0 = ri;
      rsi1 = ri+Nsite;
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
        
        sltE_i0[rsj0] = -(slt_ij - slt_ji)*cs;
        sltE_i0[rsj1] = slt_ij*cc + slt_ji*ss;
        sltE_i1[rsj0] = -slt_ij*ss - slt_ji*cc;
        sltE_i1[rsj1] = (slt_ij - slt_ji)*cs;
      }
    }
  }

  return;
}

void SlaterElmDiff(double *srOptO, const double ip, int *eleIdx) {
  const int nBuf=NSlater*NQPFull;
  const int nsize = Nsize;
  const int ne = Ne;
  const int nQPFull = NQPFull;
  const int nMPTrans = NMPTrans;
  const int nSlater = NSlater;
  const int nTrans = NMPTrans * NQPOptTrans;

  const double invIP = 1.0/ip;
  int msi,msj,ri,rj,ori,orj,tri,trj,sgni,sgnj;
  int mpidx,spidx,orbidx,qpidx,optidx,i;
  double cs,cc,ss;
  int *xqp,*xqpSgn,*xqpOpt,*xqpOptSgn;
  double *invM,*invM_i;

  int *orbitalIdx_i;
  int *transOrbIdx; /* transOrbIdx[mpidx][msi][msj] */
  int *transOrbSgn; /* transOrbSgn[mpidx][msi][msj] */
  int *tOrbIdx,*tOrbIdx_i;
  int *tOrbSgn,*tOrbSgn_i;
  double *buf, *buffer;
  double tmp;

  RequestWorkSpaceInt(2*nTrans*Nsize*Nsize);
  RequestWorkSpaceDouble(NQPFull*NSlater);

  transOrbIdx = GetWorkSpaceInt(nTrans*Nsize*Nsize); /* transOrbIdx[mpidx][msi][msj] */
  transOrbSgn = GetWorkSpaceInt(nTrans*Nsize*Nsize); /* transOrbSgn[mpidx][msi][msj] */
  buffer = GetWorkSpaceDouble(NQPFull*NSlater);

  for(i=0;i<nBuf;i++) buffer[i]=0.0;

  #pragma omp parallel for default(shared)                        \
    private(qpidx,optidx,mpidx,msi,msj,xqp,xqpSgn,xqpOpt,xqpOptSgn,\
            ri,ori,tri,sgni,rj,orj,trj,sgnj,                       \
            tOrbIdx,tOrbIdx_i,tOrbSgn,tOrbSgn_i,orbitalIdx_i)
  #pragma loop noalias
  for(qpidx=0;qpidx<nTrans;qpidx++) {
    optidx = qpidx / nMPTrans;
    mpidx = qpidx % nMPTrans;

    xqpOpt = QPOptTrans[optidx];
    xqpOptSgn = QPOptTransSgn[optidx];
    xqp = QPTrans[mpidx];
    xqpSgn = QPTransSgn[mpidx];
    tOrbIdx = transOrbIdx + qpidx*nsize*nsize;
    tOrbSgn = transOrbSgn + qpidx*nsize*nsize;
    for(msi=0;msi<nsize;msi++) {
      ri = eleIdx[msi];
      ori = xqpOpt[ri];
      tri = xqp[ori];
      sgni = xqpSgn[ori]*xqpOptSgn[ri];
      tOrbIdx_i = tOrbIdx + msi*nsize;
      tOrbSgn_i = tOrbSgn + msi*nsize;
      orbitalIdx_i = OrbitalIdx[tri];
      for(msj=0;msj<nsize;msj++) {
        rj = eleIdx[msj];
        orj = xqpOpt[rj];
        trj = xqp[orj];
        sgnj = xqpSgn[orj]*xqpOptSgn[rj];
        tOrbIdx_i[msj] = orbitalIdx_i[trj];
        tOrbSgn_i[msj] = sgni*sgnj*OrbitalSgn[tri][trj];
      }
    }
  }

  #pragma omp parallel for default(shared)        \
    private(qpidx,mpidx,spidx,cs,cc,ss,                   \
            tOrbIdx,tOrbSgn,invM,buf,msi,msj,             \
            tOrbIdx_i,tOrbSgn_i,invM_i,orbidx)
  #pragma loop noalias
  for(qpidx=0;qpidx<nQPFull;qpidx++) {
    mpidx = qpidx / NSPGaussLeg;
    spidx = qpidx % NSPGaussLeg;

    cs = PfM[qpidx] * SPGLCosSin[spidx];
    cc = PfM[qpidx] * SPGLCosCos[spidx];
    ss = PfM[qpidx] * SPGLSinSin[spidx];

    tOrbIdx = transOrbIdx + mpidx*nsize*nsize;
    tOrbSgn = transOrbSgn + mpidx*nsize*nsize;
    invM = InvM + qpidx*Nsize*Nsize;
    buf = buffer + qpidx*NSlater;

    #pragma loop norecurrence
    for(msi=0;msi<ne;msi++) {
      tOrbIdx_i = tOrbIdx + msi*nsize;
      tOrbSgn_i = tOrbSgn + msi*nsize;
      invM_i = invM + msi*nsize;
      for(msj=0;msj<ne;msj++) {
        /* si=0 sj=0*/
        orbidx = tOrbIdx_i[msj];
        buf[orbidx] += invM_i[msj]*cs*tOrbSgn_i[msj];
      }
      for(msj=ne;msj<nsize;msj++) {
        /* si=0 sj=1*/
        orbidx = tOrbIdx_i[msj];
        buf[orbidx] -= invM_i[msj]*cc*tOrbSgn_i[msj];
      }
    }
    #pragma loop norecurrence
    for(msi=ne;msi<nsize;msi++) {
      tOrbIdx_i = tOrbIdx + msi*nsize;
      tOrbSgn_i = tOrbSgn + msi*nsize;
      invM_i = invM + msi*nsize;
      for(msj=0;msj<ne;msj++) {
        /* si=1 sj=0*/
        orbidx = tOrbIdx_i[msj];
        buf[orbidx] += invM_i[msj]*ss*tOrbSgn_i[msj];
      }
      for(msj=ne;msj<nsize;msj++) {
        /* si=1 sj=1*/
        orbidx = tOrbIdx_i[msj];
        buf[orbidx] -= invM_i[msj]*cs*tOrbSgn_i[msj];
      }
    }
  }

  /* store SROptO[] */
  for(orbidx=0;orbidx<nSlater;orbidx++) {
    srOptO[orbidx] = 0.0;
  }
  #pragma loop noalias
  for(qpidx=0;qpidx<nQPFull;qpidx++) {
    tmp = QPFullWeight[qpidx];
    buf = buffer + qpidx*nSlater;
    for(orbidx=0;orbidx<nSlater;orbidx++) {
      srOptO[orbidx] += tmp * buf[orbidx];
    }
  }
  for(orbidx=0;orbidx<nSlater;orbidx++) {
    srOptO[orbidx] *= invIP;
    //printf("DEBUG:orbidx=%d srOptO=%lf %lf\n",orbidx,srOptO[orbidx],invIP);
  }

  ReleaseWorkSpaceInt();
  ReleaseWorkSpaceDouble();
  return;
}
