/*-------------------------------------------------------------
 * Variational Monte Carlo
 * matrix Package (LAPACK and Pfapack)
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/


int CalculateMAll_fcmp(const int *eleIdx, const int qpStart, const int qpEnd);
int calculateMAll_child_fcmp(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
                        double complex *bufM, int *iwork, double complex *work, int lwork,double *rwork);

#ifdef _SYSTEM_A
  #define M_DGETRF DGETRF
  #define M_DGETRI DGETRI
  #define M_DSKPFA DSKPFA
  #define M_ZGETRF ZGETRF
  #define M_ZGETRI ZGETRI
  #define M_ZSKPFA ZSKPFA
#elif _lapack_small_nounderscore
  #define M_DGETRF dgetrf
  #define M_DGETRI dgetri
  #define M_DSKPFA dskpfa
  #define M_ZGETRF zgetrf
  #define M_ZGETRI zgetri
  #define M_ZSKPFA zskpfa
#else
  #define M_DGETRF dgetrf_
  #define M_DGETRI dgetri_
  #define M_DSKPFA dskpfa_
  #define M_ZGETRF zgetrf_
  #define M_ZGETRI zgetri_
  #define M_ZSKPFA zskpfa_
#endif

#define D_PfLimit 1.0e-100

int M_DGETRF(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
int M_DGETRI(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
int M_DSKPFA(const char *uplo, const char *mthd, const int *n,
             double *a, const int *lda, double *pfaff, int *iwork,
             double *work, const int *lwork, int *info);
int M_ZGETRF(int *m, int *n, double complex *a, int *lda, int *ipiv, int *info);
int M_ZGETRI(int *n, double complex *a, int *lda, int *ipiv, double complex *work, int *lwork, int *info);
int M_ZSKPFA(const char *uplo, const char *mthd, const int *n,
             double complex *a, const int *lda, double complex *pfaff, int *iwork,
             double complex *work, const int *lwork, double *rwork, int *info);


int getLWork() {
  char uplo='U', mthd='P';
  int n,lda,lwork,info=0;
  double pfaff;
  int iwork;
  double a;
  double optSize1,optSize2;

  /* ask the optimal size of work */
  n=lda=Nsize;
  lwork=-1;
  M_DGETRI(&n, &a, &lda, &iwork, &optSize1, &lwork, &info);
  lwork=-1;
  M_DSKPFA(&uplo, &mthd, &n, &a, &lda, &pfaff, &iwork, &optSize2, &lwork, &info);

  lwork = (optSize1>optSize2) ? (int)optSize1 : (int)optSize2;
  return lwork;
}

int getLWork_fcmp() {
  char uplo='U', mthd='P';
  int n,lda,lwork,info=0;
  double rwork;
  double complex pfaff;
  int iwork;
  double complex a;
  double complex optSize1,optSize2;

  /* ask the optimal size of work */
  n=lda=Ne;
  lwork=-1;
  M_ZGETRI(&n, &a, &lda, &iwork, &optSize1, &lwork, &info);
  lwork=-1;
  M_ZSKPFA(&uplo, &mthd, &n, &a, &lda, &pfaff, &iwork, &optSize2, &lwork, &rwork, &info);

  lwork = (creal(optSize1)>creal(optSize2)) ? (int)creal(optSize1) : (int)creal(optSize2);
  return lwork;
}


//==============s fcmp =============//
/* Calculate PfM and InvM from qpidx=qpStart to qpEnd */
int CalculateMAll_fcmp(const int *eleIdx, const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  int qpidx;

  int info = 0;

  double complex *myBufM;
  double complex *myWork;
  int            *myIWork;
  int             myInfo;
  double         *myRWork;

  RequestWorkSpaceThreadInt(Nsize);         //int

  RequestWorkSpaceThreadComplex(Nsize*Nsize+LapackLWork);

  RequestWorkSpaceThreadDouble(LapackLWork); // TBC for rwork

  #pragma omp parallel default(shared)              \
    private(myIWork,myWork,myInfo,myBufM)
  {
    myIWork = GetWorkSpaceThreadInt(Nsize); // int

    myBufM  = GetWorkSpaceThreadComplex(Nsize*Nsize); //comp
    myWork  = GetWorkSpaceThreadComplex(LapackLWork); // comp

    myRWork = GetWorkSpaceThreadDouble(LapackLWork); //TBC for rwork

    #pragma omp for private(qpidx)
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      if(info!=0) continue;
      
      myInfo = calculateMAll_child_fcmp(eleIdx, qpStart, qpEnd, qpidx,
                                   myBufM, myIWork, myWork, LapackLWork,myRWork);
      if(myInfo!=0) {
        #pragma omp critical
        info=myInfo;
      }
    }
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadComplex();
  ReleaseWorkSpaceThreadDouble();
  return info;
}

int calculateMAll_child_fcmp(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
                        double complex *bufM, int *iwork, double complex *work, int lwork,double *rwork) {
  #pragma procedure serial
  /* const int qpNum = qpEnd-qpStart; */
  int msi,msj;
  int rsi,rsj;

  char uplo='U', mthd='P';
  int m,n,lda,info=0;
  double complex pfaff;

  /* optimization for Kei */
  const int nsize = Nsize;

  const double complex *sltE = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2;
  const double complex *sltE_i;

  double complex *invM = InvM + qpidx*Nsize*Nsize;
  double complex *invM_i;

  double complex *bufM_i, *bufM_i2;

  m=n=lda=Nsize;

  /* store bufM */
  /* Note that bufM is column-major and skew-symmetric. */
  /* bufM[msj][msi] = -sltE[rsi][rsj] */
  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    rsi = eleIdx[msi] + (msi/Ne)*Nsite;
    bufM_i = bufM + msi*Nsize;
    sltE_i = sltE + rsi*Nsite2;
    #pragma loop norecurrence
    for(msj=0;msj<nsize;msj++) {
      rsj = eleIdx[msj] + (msj/Ne)*Nsite;
      bufM_i[msj] = -sltE_i[rsj];
//      printf("DEBUG: msi=%d msj=%d bufM=%lf %lf \n",msi,msj,creal(bufM_i[msj]),cimag(bufM_i[msj]));
    }
  }

  /* copy bufM to invM */
  /* For Pfaffian calculation, invM is used as second buffer */
  #pragma loop noalias
  for(msi=0;msi<nsize*nsize;msi++) {
    invM[msi] = bufM[msi];
  }
  /* calculate Pf M */
  //printf("DEBUG: n=%d \n",n);
  //for(msi=0;msi<nsize;msi++){
  //  for(msj=0;msj<nsize;msj++){
  //    printf("DEBUG: msi=%d msj=%d bufM %lf %lf \n",msi,msj,creal(bufM[msi+msj*n]),cimag(bufM[msi+msj*n]));
  //  }
  //}
  M_ZSKPFA(&uplo, &mthd, &n, invM, &lda, &pfaff, iwork, work, &lwork, rwork, &info); //TBC
  //printf("DEBUG: pfaff=%lf %lf\n",creal(pfaff),cimag(pfaff));
  if(info!=0) return info;
  if(!isfinite(pfaff)) return qpidx+1;
  PfM[qpidx] = pfaff;

  /* DInv */
  M_ZGETRF(&m, &n, bufM, &lda, iwork, &info); /* ipiv = iwork */
  if(info!=0) return info;
  //for(msi=0;msi<nsize*nsize;msi++) {
  //  printf("DEBUG: M_ZGETRF %d %lf %lf \n",msi,creal(bufM[msi]),cimag(bufM[msi]));
  //}
  
  M_ZGETRI(&n, bufM, &lda, iwork, work, &lwork, &info);
  if(info!=0) return info;
  
  /* store InvM */
  /* BufM is column-major, InvM is row-major */
  #pragma loop noalias
  for(msi=0;msi<nsize;msi++) {
    invM_i = invM + msi*Nsize;
    bufM_i = bufM + msi*Nsize;
    bufM_i2 = bufM + msi;
    for(msj=0;msj<nsize;msj++) {
      invM_i[msj] = 0.5*(bufM_i2[msj*nsize] - bufM_i[msj]);
      //printf("DEBUG: msj=%d invM=%lf %lf \n",msj,creal(invM_i[msj]),cimag(invM_i[msj]));
      /* invM[i][j] = 0.5*(bufM[i][j]-bufM[j][i]) */
    }
  }

  return info;
}


//==============e fcmp =============//
