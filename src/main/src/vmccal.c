/*-------------------------------------------------------------
 * Variational Monte Carlo
 * calculate physical quantities
 *-------------------------------------------------------------
 * by Satoshi Morita 
 *-------------------------------------------------------------*/

void VMCMainCal(MPI_Comm comm);
void clearPhysQuantity();
void calculateOptTransDiff(double complex *srOptO, const double complex ipAll);
void calculateOO(double complex *srOptOO, double complex *srOptHO, const double complex *srOptO,
                 const double w, const double complex e, const int srOptSize);
void calculateOO_Store(double complex *srOptOO, double complex *srOptHO,  double complex *srOptO,
                 const double w, const double complex e,  int srOptSize, int sampleSize);
void calculateQQQQ(double *qqqq, const double *lslq, const double w, const int nLSHam);
void calculateQCAQ(double *qcaq, const double *lslca, const double *lslq,
                   const double w, const int nLSHam, const int nCA);
void calculateQCACAQ(double *qcacaq, const double *lslca, const double w,
                     const int nLSHam, const int nCA, const int nCACA,
                     int **cacaIdx);

void VMCMainCal(MPI_Comm comm) {
  int *eleIdx,*eleCfg,*eleNum,*eleProjCnt;
  double complex e,x,ip;
  double w;
  double we,sqrtw;

  const int qpStart=0;
  const int qpEnd=NQPFull;
  int sample,sampleStart,sampleEnd,sampleSize;
  int i,info;

  /* optimazation for Kei */
  const int nProj=NProj;
  double complex *srOptO = SROptO;

  int rank,size,int_i;
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  SplitLoop(&sampleStart,&sampleEnd,NVMCSample,rank,size);

  /* initialization */
  clearPhysQuantity();

  for(sample=sampleStart;sample<sampleEnd;sample++) {
    eleIdx = EleIdx + sample*Nsize;
    eleCfg = EleCfg + sample*Nsite2;
    eleNum = EleNum + sample*Nsite2;
    eleProjCnt = EleProjCnt + sample*NProj;

    StartTimer(40);
    info = CalculateMAll_fcmp(eleIdx,qpStart,qpEnd);
    StopTimer(40);

    if(info!=0) {
      fprintf(stderr,"waring: VMCMainCal rank:%d sample:%d info:%d (CalculateMAll)\n",rank,sample,info);
      continue;
    }

    ip = CalculateIP_fcmp(PfM,qpStart,qpEnd,MPI_COMM_SELF);
   // printf("DEBUG:  ip= %lf %lf\n",creal(ip),cimag(ip));
    x = LogProjVal(eleProjCnt);
    /* calculate reweight */
    //w = exp(2.0*(log(fabs(ip))+x) - logSqPfFullSlater[sample]);
    w =1.0;
    if( !isfinite(w) ) {
      fprintf(stderr,"waring: VMCMainCal rank:%d sample:%d w=%e\n",rank,sample,w);
      continue;
    }

    StartTimer(41);
    /* calculate energy */
    e = CalculateHamiltonian(ip,eleIdx,eleCfg,eleNum,eleProjCnt);
    //printf("DEBUG:  e= %lf %lf\n",creal(e),cimag(e));
    StopTimer(41);
    if( !isfinite(e) ) {
      fprintf(stderr,"waring: VMCMainCal rank:%d sample:%d e=%e\n",rank,sample,creal(e)); //TBC
      continue;
    }

    Wc += w;
    Etot  += w * e;
    Etot2 += w * conj(e) * e;

    if(NVMCCalMode==0) {
      /* Calculate O for correlation fauctors */
      SROptO[0] = 1.0+0.0*I;//   real 
      SROptO[1] = 0.0+0.0*I;//   real 
      #pragma loop noalias
      for(i=0;i<nProj;i++){ 
        srOptO[(i+1)*2]     = (double)(eleProjCnt[i]); // even real
        srOptO[(i+1)*2+1]   = 0.0+0.0*I;               // odd  comp
      }

      StartTimer(42);
      /* SlaterElmDiff */
      SlaterElmDiff_fcmp(SROptO+2*NProj+2,ip,eleIdx); //TBC
      //printf("DEBUG:  NPara=%d NProj=%d NSlater=%d \n",NPara,NProj,NSlater);
      //for(i=0;i<NPara+1;i++){ 
      //  printf("DEBUG:  i=%d %lf %lf\n",i,creal(SROptO[i*2]),cimag(SROptO[i*2+1]));
      //}
      StopTimer(42);
      
      if(FlagOptTrans>0) { // this part will be not used
        calculateOptTransDiff(SROptO+2*NProj+2*NSlater+2, ip); //TBC
      }

      //for(i=0;i<NPara*2+2;i++){ 
      //  printf("DEBUG: i=%d %lf %lf \n",i,creal(srOptO[i]),cimag(srOptO[i]));
      //}
      StartTimer(43);
      /* Calculate OO and HO */
      if(NStoreO==0){
        calculateOO(SROptOO,SROptHO,SROptO,w,e,SROptSize);
      }else{
        we    = w*e;
        sqrtw = sqrt(w); 
        #pragma omp parallel for default(shared) private(int_i)
        for(int_i=0;int_i<SROptSize;int_i++){
          // SROptO_Store for fortran
          SROptO_Store[int_i+sample*SROptSize]  = sqrtw*SROptO[int_i];
          SROptHO[int_i]                       += we*SROptO[int_i]; 
        }
      } 
      StopTimer(43);

    } else if(NVMCCalMode==1) {
      StartTimer(42);
      /* Calculate Green Function */
      CalculateGreenFunc(w,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
      StopTimer(42);

      if(NLanczosMode>0){
        // ignoring Lanczos: to be added
        /* Calculate local QQQQ */
        //StartTimer(43);
        //LSLocalQ(e,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
        //calculateQQQQ(QQQQ,LSLQ,w,NLSHam);
        //StopTimer(43);
        //if(NLanczosMode>1){
          /* Calculate local QcisAjsQ */
          //StartTimer(44);
          //LSLocalCisAjs(e,ip,eleIdx,eleCfg,eleNum,eleProjCnt);
          //calculateQCAQ(QCisAjsQ,LSLCisAjs,LSLQ,w,NLSHam,NCisAjs);
          //calculateQCACAQ(QCisAjsCktAltQ,LSLCisAjs,w,NLSHam,NCisAjs,
          //                NCisAjsCktAlt,CisAjsCktAltIdx);
          //StopTimer(44);
        //}
      }
    }
  } /* end of for(sample) */

// calculate OO and HO at NVMCCalMode==0
  if(NStoreO!=0 && NVMCCalMode==0){
    sampleSize=sampleEnd-sampleStart;
    StartTimer(45);
    calculateOO_Store(SROptOO,SROptHO,SROptO_Store,w,e,SROptSize,sampleSize);
    StopTimer(45);
  }
  return;
}

void clearPhysQuantity(){
  int i,n;
  double complex *vec;
  Wc = Etot = Etot2 = 0.0;
  if(NVMCCalMode==0) {
    /* SROptOO, SROptHO, SROptO */
    n = SROptSize*(SROptSize+2);
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

void calculateOptTransDiff(double complex *srOptO, const double complex ipAll) {
  int i,j;
  double complex ip;
  double complex *pfM;

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

void calculateOO_Store(double complex *srOptOO, double complex *srOptHO, double complex *srOptO_Store,
                 const double w, const double complex e, int srOptSize, int sampleSize) {

  //#define M_DGEM dgemm_

  extern int dgemm_(char *jobz, char *uplo, int *m,int *n,int *k,double *alpha,  double *a, int *lda, double *b, int *ldb,
                    double *beta,double *c,int *ldc);

  char jobz, uplo;
  double alpha,beta;
  
  alpha = 1.0;
  beta  = 0.0;
  
  jobz = 'N';
  uplo = 'T';
  //dgemm_(&jobz,&uplo,&srOptSize,&srOptSize,&sampleSize,&alpha,srOptO_Store,&srOptSize,srOptO_Store,&srOptSize,&beta,srOptOO,&srOptSize);

  return;
}




//void calculateOO(double complex *srOptOO, double complex *srOptHO, const double complex *srOptO,
//                 const double w, const double complex e, const int srOptSize) {
//  double we=w*e;
//
//  #define M_DAXPY daxpy_
//  #define M_DGER dger_
//
//  extern int M_DAXPY(const int *n, const double *alpha, const double *x, const int *incx,
//                     double *y, const int *incy);
//  extern int M_DGER(const int *m, const int *n, const double *alpha,
//                    const double *x, const int *incx, const double *y, const int *incy, 
//                    double *a, const int *lda);
//  int m,n,incx,incy,lda;
//  m=n=lda=srOptSize;
//  incx=incy=1;
//
//  /* OO[i][j] += w*O[i]*O[j] */
//  M_DGER(&m, &n, &w, srOptO, &incx, srOptO, &incy, srOptOO, &lda);
//
//  /* HO[i] += w*e*O[i] */
//  M_DAXPY(&n, &we, srOptO, &incx, srOptHO, &incy);
//
//  return;
//}

void calculateOO(double complex *srOptOO, double complex *srOptHO, const double complex *srOptO,
                 const double w, const double complex e, const int srOptSize){
  int i,j;
  double complex tmp;
  //#pragma omp parallel for default(shared)        \
    private(i,j,tmp,srOptOO)
  //#pragma loop noalias
  for(i=0;i<2*srOptSize;i++) {
    tmp         = w * srOptO[i];
    //printf("i=%d %lf \n",i,creal(tmp));
    srOptHO[i] += e * tmp;
    for(j=0;j<2*srOptSize;j++) {
      srOptOO[j+i*(2*srOptSize)] += w*srOptO[j]*(srOptO[i]); // TBC
    }
    /* HO[i] += w*e*O[i] */
  }

  return;
}

void calculateQQQQ(double *qqqq, const double *lslq, const double w, const int nLSHam) {
  const int n=nLSHam*nLSHam*nLSHam*nLSHam;
  int rq,rp,ri,rj;
  int i,tmp;

  /* QQQQ[rq][rp][ri][rj] += w * LSLQ[rq][ri] * LSLQ[rp][rj] */
  # pragma omp parallel for default(shared) private(i,tmp,rq,rp,ri,rj)
  for(i=0;i<n;++i) {
    rj = i%nLSHam;   tmp=i/nLSHam;
    ri = tmp%nLSHam; tmp=tmp/nLSHam;
    rp = tmp%nLSHam; tmp=tmp/nLSHam;
    rq = tmp%nLSHam;

    qqqq[i] += w * lslq[rq*nLSHam+ri] * lslq[rp*nLSHam+rj];
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

