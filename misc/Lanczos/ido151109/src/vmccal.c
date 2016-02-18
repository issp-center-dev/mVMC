/*-------------------------------------------------------------
 * Variational Monte Carlo
 * calculate physical quantities
 *-------------------------------------------------------------
 * by Satoshi Morita and Ryui Kaneko
 *-------------------------------------------------------------*/

void VMCMainCal(MPI_Comm comm);
void clearPhysQuantity();
void calculateOptTransDiff(double *srOptO, const double ipAll);
void calculateOO(double *srOptOO, double *srOptHO, const double *srOptO,
                 const double w, const double e, const int srOptSize);
void calculateOO_Store(double *srOptOO, double *srOptO,
                 const double w, const double e,  int srOptSize, int sampleSize);
void calculateHHHH_Store(double *lsHHHH, double *lsHH_Store, int lsSize, int sampleSize);
void calculateQQQQ_Store(double *qqqq, double *lsHHHH, const int nLSHam, const int nHHHH);
void calculateQQQQ(double *qqqq, const double *lslq, const double w, const int nLSHam,const int size);
void calculateQCAQ(double *qcaq, const double *lslca, const double *lslq,
                   const double w, const int nLSHam, const int nCA);
void calculateQCACAQ(double *qcacaq, const double *lslca, const double w,
                     const int nLSHam, const int nCA, const int nCACA,
                     int **cacaIdx);

void VMCMainCal(MPI_Comm comm) {
  int *eleIdx,*eleCfg,*eleNum,*eleProjCnt;
  double e,x,w,ip;
  double we,sqrtw;
  double wec,wek,ec,ek;

  const int qpStart=0;
  const int qpEnd=NQPFull;
  int sample,sampleStart,sampleEnd,sampleSize;
  int transferStart,transferEnd;
  int dcStart,dcEnd;
  int i,info;
  double *lsHH, *lsW, *lsH,*lsHcH, *lsHkH;
  double recv[NVMCSample*(NSplitSize+3)];
  //double *recv;

  /* optimazation for Kei */
  const int nProj=NProj;
  double *srOptO = SROptO;

  int rank,size,int_i;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  sampleStart=0;
  sampleEnd=NVMCSample;
  //SplitLoop(&sampleStart,&sampleEnd,NVMCSample,rank,size);
  SplitLoop(&transferStart,&transferEnd,NTransfer,rank,size);
  //SplitLoop(&dcStart,&dcEnd,NCisAjsCktAltDC,rank,size);
//printf("tStart=%d,End=%d\n",transferStart,transferEnd);
  /* initialization */
  clearPhysQuantity();

  for(sample=sampleStart;sample<sampleEnd;sample++) {
    eleIdx = EleIdx + sample*Nsize;
    eleCfg = EleCfg + sample*Nsite2;
    eleNum = EleNum + sample*Nsite2;
    eleProjCnt = EleProjCnt + sample*NProj;

    StartTimer(40);
    info = CalculateMAll(eleIdx,qpStart,qpEnd);
    StopTimer(40);

    if(info!=0) {
      fprintf(stderr,"waring: VMCMainCal rank:%d sample:%d info:%d (CalculateMAll)\n",rank,sample,info);
      continue;
    }

    ip = CalculateIP(PfM,qpStart,qpEnd,MPI_COMM_SELF);
    x = LogProjVal(eleProjCnt);
    /* calculate reweight */
    w = exp(2.0*(log(fabs(ip))+x) - logSqPfFullSlater[sample]);
    if( !isfinite(w) ) {
      fprintf(stderr,"waring: VMCMainCal rank:%d sample:%d w=%e\n",rank,sample,w);
      continue;
    }

    StartTimer(41);
    /* calculate energy */
    e  = CalculateHamiltonian(ip,eleIdx,eleCfg,eleNum,eleProjCnt);
    //ec  = CalculateHamiltonian0(eleNum);
    //ek  = CalculateHamiltonian1(ip,eleIdx,eleCfg,eleNum,eleProjCnt,transferStart,transferEnd);
   // e = ec + ek;
    StopTimer(41);
    if( !isfinite(e) ) {
      fprintf(stderr,"waring: VMCMainCal rank:%d sample:%d e=%e\n",rank,sample,e);
      continue;
    }

    Wc += w/(double)size;
    Etot  += w * e;
    //Ec += w*ec/(double)size;
    //Ec2 += w*ec*ec/(double)size;
    //Ek += w*ek;
    //Eck += w*ec*ek;
    //Ek2 += w*ek*ek;//need ek x ek' in other rank->Store ek for Etot2
    //printf("ec=%.2e, ek=%.2e\n",ec,ek);
    Etot2 += w * e * e;

    if(NVMCCalMode==0) {
      /* Calculate O for correlation fauctors */
      SROptO[0] = 1.0;
      #pragma loop noalias
      for(i=0;i<nProj;i++) srOptO[i+1] = (double)(eleProjCnt[i]);

      StartTimer(42);
      /* SlaterElmDiff */
      SlaterElmDiff(SROptO+NProj+1,ip,eleIdx);
      StopTimer(42);
      
      if(FlagOptTrans>0) {
        calculateOptTransDiff(SROptO+NProj+NSlater+1, ip);
      }

      StartTimer(43);
      /* Calculate OO and HO */
      if(NStoreO==0){
        calculateOO(SROptOO,SROptHO,SROptO,w,e,SROptSize);
      }else{
        //we    = w*e;
        wec    = w*ec/(double)size;
        wek    = w*ek;
        sqrtw = sqrt(w/(double)size); 
        #pragma omp parallel for default(shared) private(int_i)
        for(int_i=0;int_i<SROptSize;int_i++){
          // SROptO_Store for fortran
          SROptO_Store[int_i+sample*SROptSize]  = sqrtw*SROptO[int_i];
          //SROptHO[int_i]                       += we*SROptO[int_i]; 
          SROptHcO[int_i]                       += wec*SROptO[int_i];
          SROptHkO[int_i]                       += wek*SROptO[int_i]; 
        }
      }
      StopTimer(43);
    } else if(NVMCCalMode==1) {
      StartTimer(42);
      /* Calculate Green Function */
      CalculateGreenFunc(w,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
      StopTimer(42);

      if(NLanczosMode>0){
        /* Calculate local QQQQ */
        StartTimer(43);
        if(NStoreO==0){
          LSLocalQ(e,ip,eleIdx,eleCfg,eleNum,eleProjCnt,transferStart,transferEnd,comm);
          calculateQQQQ(QQQQ,LSLQ,w,NLSHam,size);
        }else{
          sqrtw = sqrt(w);
          lsHH = LSHH_Store  + sample*(3+NSplitSize);
        
          LSLocalQ_Store(sqrtw, e,ip,eleIdx,eleCfg,eleNum,eleProjCnt,
                       transferStart,transferEnd,lsHH,rank);
        }
        StopTimer(43);
        if(NLanczosMode>1){
          /* Calculate local QcisAjsQ */
          StartTimer(44);
          LSLocalCisAjs(e,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
          calculateQCAQ(QCisAjsQ,LSLCisAjs,LSLQ,w,NLSHam,NCisAjs);
          calculateQCACAQ(QCisAjsCktAltQ,LSLCisAjs,w,NLSHam,NCisAjs,
                          NCisAjsCktAlt,CisAjsCktAltIdx);
          StopTimer(44);
        }
      }
    }
  } /* end of for(sample) */

// calculate OO at NVMCCalMode==0
  if(NStoreO!=0){
    sampleSize=sampleEnd-sampleStart;
    if(NVMCCalMode==0){
      StartTimer(45);
      calculateOO_Store(SROptOO,SROptO_Store,w,e,SROptSize,sampleSize);
      StopTimer(45);
    }else if(NLanczosMode>0){
      StartTimer(45);
      if(size>1) {
        //SafeMpiReduce(LSHH_Store, recv, (NSplitSize+3)*sampleSize, comm);
        //if(rank==0){
        //  #pragma loop noalias
        //  for(i=0;i<(NSplitSize+3)*sampleSize;i++) LSHH_Store[i] = recv[i];
        //}
        SafeMpiAllReduce(LSHH_Store, recv, (NSplitSize+3)*sampleSize, comm);
        #pragma loop noalias
        for(i=0;i<(NSplitSize+3)*sampleSize;i++) LSHH_Store[i] = recv[i]/sqrt((double)NSplitSize);
      }
      //if(rank==0){
        calculateHHHH_Store(LSHHHH,LSHH_Store,3+NSplitSize,sampleSize);
        calculateQQQQ_Store(QQQQ,LSHHHH,NLSHam,3+NSplitSize);
      //}
      StopTimer(45);
    }
  }
  return;
}

void clearPhysQuantity(){
  int i,n;
  double *vec;
  Wc = Etot = Etot2 = 0.0;
  Ek = Ek2 = Ec = Ec2 = Eck = 0.0;
  if(NVMCCalMode==0) {
    /* SROptOO, SROptHcO, SROptHkO, SROptHO, SROptO */
    n = SROptSize*(SROptSize+4);
    vec = SROptOO;
    for(i=0;i<n;i++) vec[i] = 0.0;
  } else if(NVMCCalMode==1) {
    /* CisAjs, CisAjsCktAlt, CisAjsCktAltDC */
    n = 2*NCisAjs+NCisAjsCktAlt+NCisAjsCktAltDC;
    vec = PhysCisAjs;
    for(i=0;i<n;i++) vec[i] = 0.0;
    if(NLanczosMode>0) {
      /* QQQQ, LSLQ */
      n = NLSHam*NLSHam*NLSHam*NLSHam + NLSHam*NLSHam;
      vec = QQQQ;
      for(i=0;i<n;i++) vec[i] = 0.0;
      if(NStoreO!=0){
        n = (3+NSplitSize)*NVMCSample;
        vec = LSHH_Store;
        for(i=0;i<n;i++) vec[i] = 0.0;
        n = (3+NSplitSize)*(3+NSplitSize);
        vec = LSHHHH;
        for(i=0;i<n;i++) vec[i] = 0.0;
      }
      if(NLanczosMode>1) {
        /* QCisAjsQ, QCisAjsCktAltQ, LSLCisAjs */
        n = NLSHam*NLSHam*NCisAjs + NLSHam*NLSHam*NCisAjsCktAlt
          + NLSHam*NCisAjs;
        vec = QCisAjsQ;
        for(i=0;i<n;i++) vec[i] = 0.0;
      }
    }
  }
  return;
}

void calculateOptTransDiff(double *srOptO, const double ipAll) {
  int i,j;
  double ip;
  double *pfM;

  for(i=0;i<NQPOptTrans;++i) {
    ip = 0.0;
    pfM = PfM + i*NQPFix;
    for(j=0;j<NQPFix;++j) {
      ip += QPFixWeight[j] * pfM[j];
    }
    srOptO[i] = ip/ipAll;
  }

  return;
}

void calculateOO_Store(double *srOptOO, double *srOptO_Store,
                 const double w, const double e, int srOptSize, int sampleSize) {

  //#define M_DGEM dgemm_

  extern int dgemm_(char *jobz, char *uplo, int *m,int *n,int *k,double *alpha,  double *a, int *lda, double *b, int *ldb,
                    double *beta,double *c,int *ldc);

  char jobz, uplo;
  double alpha,beta;
  
  alpha = 1.0;
  beta  = 0.0;
  
  jobz = 'N';
  uplo = 'T';
  dgemm_(&jobz,&uplo,&srOptSize,&srOptSize,&sampleSize,&alpha,srOptO_Store,&srOptSize,srOptO_Store,&srOptSize,&beta,srOptOO,&srOptSize);

  return;
}

void calculateHHHH_Store(double *lsHHHH, double *lsHH_Store, int lsSize, int sampleSize) {

  //#define M_DGEM dgemm_

  extern int dgemm_(char *jobz, char *uplo, int *m,int *n,int *k,double *alpha,  double *a, int *lda, double *b, int *ldb,
                    double *beta,double *c,int *ldc);

  char jobz, uplo;
  double alpha,beta;
  
  alpha = 1.0;
  beta  = 0.0;
  
  jobz = 'N';
  uplo = 'T';
  dgemm_(&jobz,&uplo,&lsSize,&lsSize,&sampleSize,&alpha,lsHH_Store,&lsSize,lsHH_Store,&lsSize,&beta,lsHHHH,&lsSize);

  return;
}



void calculateOO(double *srOptOO, double *srOptHO, const double *srOptO,
                 const double w, const double e, const int srOptSize) {
  double we=w*e;

  #define M_DAXPY daxpy_
  #define M_DGER dger_

  extern int M_DAXPY(const int *n, const double *alpha, const double *x, const int *incx,
                     double *y, const int *incy);
  extern int M_DGER(const int *m, const int *n, const double *alpha,
                    const double *x, const int *incx, const double *y, const int *incy, 
                    double *a, const int *lda);
  int m,n,incx,incy,lda;
  m=n=lda=srOptSize;
  incx=incy=1;

  /* OO[i][j] += w*O[i]*O[j] */
  M_DGER(&m, &n, &w, srOptO, &incx, srOptO, &incy, srOptOO, &lda);

  /* HO[i] += w*e*O[i] */
  M_DAXPY(&n, &we, srOptO, &incx, srOptHO, &incy);

  return;
}

void calculateQQQQ_Store(double *qqqq, double *lsHHHH, const int nLSHam, const int nHHHH) {
  int i,j;
  double *lsHHHH_i;
  double tmp2,tmp4,tmp5;
  tmp2=tmp4=tmp5=0.0;

  /* QQQQ[rq][rp][ri][rj] += w * LSLQ[rq][ri] * LSLQ[rp][rj] */
  qqqq[0] = lsHHHH[0];
  qqqq[1] = lsHHHH[1];
  # pragma omp parallel for default(shared) private(i,j,lsHHHH_i)\
  reduction(+:tmp2,tmp4,tmp5)
  for(i=0;i<nHHHH;++i) {
    lsHHHH_i = lsHHHH + i*nHHHH;
    for(j=i;j<nHHHH;++j) {
      if(i==0){
        if(j==i){
          qqqq[0] = lsHHHH_i[j];
        }else if(j==i+1){
          qqqq[1] = lsHHHH_i[j];
        }else{
          tmp2 += lsHHHH_i[j];
        }
      }else if(i==1){
        if(j==i){
          qqqq[3] = lsHHHH_i[j];
        }else{
          tmp4 += lsHHHH_i[j];
        }
      }else{
        if(j==i){
          tmp5 += lsHHHH_i[j];
        }else{
          tmp5 += 2.0*lsHHHH_i[j];
        }
      }
    }
  }
  qqqq[2] = tmp2;
  qqqq[4] = tmp4;
  qqqq[5] = tmp5;

  return;
}

void calculateQQQQ(double *qqqq, const double *lslq, const double w, const int nLSHam, const int size) {
  const int n=nLSHam*nLSHam*nLSHam*nLSHam;
  int rq,rp,ri,rj;
  double ratio = 1.0/(double)size;
  int i,tmp;

  /* QQQQ[rq][rp][ri][rj] += w * LSLQ[rq][ri] * LSLQ[rp][rj] */
  # pragma omp parallel for default(shared) private(i,tmp,rq,rp,ri,rj)
  for(i=0;i<n;++i) {
    rj = i%nLSHam;   tmp=i/nLSHam;
    ri = tmp%nLSHam; tmp=tmp/nLSHam;
    rp = tmp%nLSHam; tmp=tmp/nLSHam;
    rq = tmp%nLSHam;

    //if( (rq==1 && ri==1) || (rp==1 && rj ==1)){
    //  ratio = 1.0;
    //}else{
    //}

    qqqq[i] += w * lslq[rq*nLSHam+ri] * lslq[rp*nLSHam+rj] * ratio;
  }

  return;
}

void calculateQCAQ(double *qcaq, const double *lslca, const double *lslq,
                   const double w, const int nLSHam, const int nCA) {
  const int n=nLSHam*nLSHam*nCA;
  int rq,rp,idx;
  int i,tmp;

  /* QCisAjsQ[rq][rp][idx] += w * LSLCisAjs[rq][idx] * LSLQ[rp][0] */
# pragma omp parallel for default(shared) private(i,tmp,idx,rp,rq)
  for(i=0;i<n;++i) {
    idx = i%nCA;     tmp = i/nCA;
    rp = tmp%nLSHam; tmp = tmp/nLSHam;
    rq = tmp%nLSHam;

    qcaq[i] += w * lslca[rq*nCA+idx] * lslq[rp*nLSHam];
  }

  return;
}

void calculateQCACAQ(double *qcacaq, const double *lslca, const double w,
                     const int nLSHam, const int nCA, const int nCACA,
                     int **cacaIdx) {
  const int n=nLSHam*nLSHam*nCACA;
  int rq,rp,ri,rj,idx;
  int i,tmp;

  /* QCisAjsCktAltQ[rq][rp][idx] += w * LSLCisAjs[rq][ri] * LSLCisAjs[rp][rj] */
# pragma omp parallel for default(shared) private(i,tmp,idx,rp,rq,ri,rj)
  for(i=0;i<n;++i) {
    idx = i%nCACA;   tmp = i/nCACA;
    rp = tmp%nLSHam; tmp = tmp/nLSHam;
    rq = tmp%nLSHam;

    ri = cacaIdx[idx][0];
    rj = cacaIdx[idx][1];

    qcacaq[i] += w * lslca[rq*nCA+ri] * lslca[rp*nCA+rj];
  }

  return;
}

