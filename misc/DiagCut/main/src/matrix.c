/*-------------------------------------------------------------
 * Variational Monte Carlo
 * matrix Package (LAPACK and Pfapack)
 *-------------------------------------------------------------
 * by Satoshi Morita
 *-------------------------------------------------------------*/

int CalculateMAll(const int *eleIdx, const int qpStart, const int qpEnd);
int calculateMAll_child(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
                        double *bufM, int *iwork, double *work, int lwork);

#ifdef _SYSTEM_A
  #define M_DGETRF DGETRF
  #define M_DGETRI DGETRI
  #define M_DSKPFA DSKPFA
#elif _lapack_small_nounderscore
  #define M_DGETRF dgetrf
  #define M_DGETRI dgetri
  #define M_DSKPFA dskpfa
#else
  #define M_DGETRF dgetrf_
  #define M_DGETRI dgetri_
  #define M_DSKPFA dskpfa_
#endif

#define D_PfLimit 1.0e-100

int M_DGETRF(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
int M_DGETRI(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
int M_DSKPFA(const char *uplo, const char *mthd, const int *n,
             double *a, const int *lda, double *pfaff, int *iwork,
             double *work, const int *lwork, int *info);


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

/* Calculate PfM and InvM from qpidx=qpStart to qpEnd */
int CalculateMAll(const int *eleIdx, const int qpStart, const int qpEnd) {
  const int qpNum = qpEnd-qpStart;
  int qpidx;

  int info = 0;

  double *myBufM;
  double *myWork;
  int *myIWork;
  int myInfo;

  RequestWorkSpaceThreadInt(Nsize);
  RequestWorkSpaceThreadDouble(Nsize*Nsize+LapackLWork);

  #pragma omp parallel default(shared)              \
    private(myIWork,myWork,myInfo,myBufM)
  {
    myIWork = GetWorkSpaceThreadInt(Nsize);
    myBufM = GetWorkSpaceThreadDouble(Nsize*Nsize);
    myWork = GetWorkSpaceThreadDouble(LapackLWork);

    #pragma omp for private(qpidx)
    for(qpidx=0;qpidx<qpNum;qpidx++) {
      if(info!=0) continue;
      
      myInfo = calculateMAll_child(eleIdx, qpStart, qpEnd, qpidx,
                                   myBufM, myIWork, myWork, LapackLWork);
      if(myInfo!=0) {
        #pragma omp critical
        info=myInfo;
      }
    }
  }

  ReleaseWorkSpaceThreadInt();
  ReleaseWorkSpaceThreadDouble();
  return info;
}

int calculateMAll_child(const int *eleIdx, const int qpStart, const int qpEnd, const int qpidx,
                        double *bufM, int *iwork, double *work, int lwork) {
  #pragma procedure serial
  /* const int qpNum = qpEnd-qpStart; */
  int msi,msj;
  int rsi,rsj;

  char uplo='U', mthd='P';
  int m,n,lda,info=0;
  double pfaff;

  /* optimization for Kei */
  const int nsize = Nsize;

  const double *sltE = SlaterElm + (qpidx+qpStart)*Nsite2*Nsite2;
  const double *sltE_i;

  double *invM = InvM + qpidx*Nsize*Nsize;
  double *invM_i;

  double *bufM_i, *bufM_i2;

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
    }
  }

  /* copy bufM to invM */
  /* For Pfaffian calculation, invM is used as second buffer */
  #pragma loop noalias
  for(msi=0;msi<nsize*nsize;msi++) {
    invM[msi] = bufM[msi];
  }
  /* calculate Pf M */
  M_DSKPFA(&uplo, &mthd, &n, invM, &lda, &pfaff, iwork, work, &lwork, &info);
  if(info!=0) return info;
  if(!isfinite(pfaff)) return qpidx+1;
  PfM[qpidx] = pfaff;

  /* DInv */
  M_DGETRF(&m, &n, bufM, &lda, iwork, &info); /* ipiv = iwork */
  if(info!=0) return info;
  
  M_DGETRI(&n, bufM, &lda, iwork, work, &lwork, &info);
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
      /* invM[i][j] = 0.5*(bufM[i][j]-bufM[j][i]) */
    }
  }

  return info;
}
