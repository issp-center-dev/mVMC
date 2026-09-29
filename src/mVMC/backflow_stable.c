/* Non-FSZ BackFlow numerical kernels (NSPGaussLeg == 1).
 * The multi-row ratio formula subtracts O(||A^-1||) terms even for a
 * no-op; its error can grow as eps*kappa(A)^2. Evaluate the occupied
 * Pfaffian directly instead. Keep the row-local Slater construction.
 * A = [0 F; -F^T 0], A^-1 = [0 -F^-T; F^-1 0]. Invert F once and
 * mirror it, rather than separately rounding the reciprocal blocks.
 * The lower Slater triangle is authoritative, as in the full Pfaffian
 * evaluation with UPLO='U' and column-major -Slater^T. No conjugation.
 * Included by backflow_pf.c; all workspaces are local to the caller.
 */

static int BFStableBuild_real(const double *slater, const int *eleIdx,
    int globalQp, int n, const int *msa, const double *vec, double *a) {
  if(!slater || !eleIdx || !a || Nsize <= 0 || Ne <= 0 || Nsize != 2LL*Ne
     || Nsite <= 0 || Nsite2 != 2LL*Nsite || globalQp < 0 || globalQp >= NQPFull
     || n < 0 || n > Nsize || (n && (!msa || !vec))) return BF_PF_INVALID_ARGUMENT;
  int *rowMap = malloc((size_t)Nsize*sizeof(int));
  if(!rowMap) return BF_PF_INVALID_ARGUMENT;
  for(int i=0;i<Nsize;i++) {
    if(eleIdx[i] < 0 || eleIdx[i] >= Nsite) { free(rowMap); return BF_PF_INVALID_ARGUMENT; }
    rowMap[i]=-1;
  }
  for(int k=0;k<n;k++) {
    if(msa[k]<0 || msa[k]>=Nsize || rowMap[msa[k]]!=-1) {
      free(rowMap); return BF_PF_INVALID_ARGUMENT;
    }
    rowMap[msa[k]]=k;
  }
  const double *slt=slater+(size_t)globalQp*Nsite2*Nsite2;
  int status=BF_PF_OK;
  for(int i=0;i<Nsize;i++) {
    a[(size_t)i*Nsize+i]=0.0;
    for(int j=0;j<i;j++) {
      /* Unchanged entries use the base configuration; changed rows already
       * contain the candidate entries, including intersections of two rows. */
      double value;
      if(rowMap[i]>=0) value=vec[(size_t)rowMap[i]*Nsize+j];
      else if(rowMap[j]>=0) value=-vec[(size_t)rowMap[j]*Nsize+i];
      else value=slt[(size_t)(eleIdx[i]+(i/Ne)*Nsite)*Nsite2+eleIdx[j]+(j/Ne)*Nsite];
      if(!isfinite(value)) status=BF_PF_NONFINITE;
      if(i/Ne==j/Ne && value!=0.0) status=BF_PF_INVALID_ARGUMENT;
      a[(size_t)i*Nsize+j]=-value; /* column-major A */
      a[(size_t)j*Nsize+i]=value;
    }
  }
  free(rowMap);
  return status;
}

static int BFStableFactor_real(double *a, double *pf, double *inverse,
    int *pivot, double *work, int lwork, int globalQp, int stage) {
  int dim=Nsize, block=Ne, info=0, status=BF_PF_OK;
  char uplo='U', method='P';
  double *f=NULL;
  if(inverse) {
    f=malloc((size_t)block*block*sizeof(*f));
    if(!f) return BF_PF_INVALID_ARGUMENT;
    for(int j=0;j<block;j++) for(int i=0;i<block;i++)
      f[(size_t)j*block+i]=a[(size_t)(block+j)*dim+i];
  }
  M_DSKPFA(&uplo,&method,&dim,a,&dim,pf,pivot,work,&lwork,&info);
  if(stage) {
    info=BFInjectNumericalInfo(stage,globalQp,info);
    if(BFInjectNumericalNan(stage,globalQp)) *pf=NAN;
  }
  status=BFStatusFromLapack(info);
  const double value=*pf;
  if(status==BF_PF_OK && !isfinite(value)) status=BF_PF_NONFINITE;
  if(status==BF_PF_OK && inverse) {
    M_DGETRF(&block,&block,f,&block,pivot,&info);
    if(stage==BF_FAIL_ACCEPT_PF) info=BFInjectNumericalInfo(BF_FAIL_ACCEPT_LU,globalQp,info);
    status=BFStatusFromLapack(info);
    if(status==BF_PF_OK) {
      M_DGETRI(&block,f,&block,pivot,work,&lwork,&info);
      if(stage==BF_FAIL_ACCEPT_PF) info=BFInjectNumericalInfo(BF_FAIL_ACCEPT_INVERSE,globalQp,info);
      status=BFStatusFromLapack(info);
    }
    if(status==BF_PF_OK) {
      memset(inverse,0,(size_t)dim*dim*sizeof(*inverse));
      for(int j=0;j<block;j++) for(int i=0;i<block;i++) {
        const double value=f[(size_t)j*block+i];
        if(!isfinite(value)) status=BF_PF_NONFINITE;
        inverse[(size_t)(block+i)*dim+j]=value;
        inverse[(size_t)j*dim+block+i]=-value;
      }
    }
  }
  free(f);
  return status;
}

static int BFStableEvaluate_real(const double *slater, const int *eleIdx,
    int globalQp, int n, const int *msa, const double *vec,
    double *pf, double *inverse, int stage) {
  if(Nsize<=0 || LapackLWork<Nsize || !pf
     || (size_t)Nsize>SIZE_MAX/(size_t)Nsize/sizeof(double)
     || (size_t)LapackLWork>SIZE_MAX/sizeof(double)) return BF_PF_INVALID_ARGUMENT;
  double *a=malloc((size_t)Nsize*Nsize*sizeof(*a));
  double *work=malloc((size_t)LapackLWork*sizeof(*work));
  int *pivot=malloc((size_t)Nsize*sizeof(*pivot));

  int status=BF_PF_INVALID_ARGUMENT;
  if(a && work && pivot) {
    status=BFStableBuild_real(slater,eleIdx,globalQp,n,msa,vec,a);
    if(status==BF_PF_OK) status=BFStableFactor_real(a,pf,inverse,pivot,work,LapackLWork,globalQp,stage);
  }
  free(a); free(work); free(pivot);
  return status;
}

static int BFStableBuild_fcmp(const double complex *slater, const int *eleIdx,
    int globalQp, int n, const int *msa, const double complex *vec, double complex *a) {
  if(!slater || !eleIdx || !a || Nsize <= 0 || Ne <= 0 || Nsize != 2LL*Ne
     || Nsite <= 0 || Nsite2 != 2LL*Nsite || globalQp < 0 || globalQp >= NQPFull
     || n < 0 || n > Nsize || (n && (!msa || !vec))) return BF_PF_INVALID_ARGUMENT;
  int *rowMap = malloc((size_t)Nsize*sizeof(int));
  if(!rowMap) return BF_PF_INVALID_ARGUMENT;
  for(int i=0;i<Nsize;i++) {
    if(eleIdx[i] < 0 || eleIdx[i] >= Nsite) { free(rowMap); return BF_PF_INVALID_ARGUMENT; }
    rowMap[i]=-1;
  }
  for(int k=0;k<n;k++) {
    if(msa[k]<0 || msa[k]>=Nsize || rowMap[msa[k]]!=-1) {
      free(rowMap); return BF_PF_INVALID_ARGUMENT;
    }
    rowMap[msa[k]]=k;
  }
  const double complex *slt=slater+(size_t)globalQp*Nsite2*Nsite2;
  int status=BF_PF_OK;
  for(int i=0;i<Nsize;i++) {
    a[(size_t)i*Nsize+i]=0.0;
    for(int j=0;j<i;j++) {
      /* Unchanged entries use the base configuration; changed rows already
       * contain the candidate entries, including intersections of two rows. */
      double complex value;
      if(rowMap[i]>=0) value=vec[(size_t)rowMap[i]*Nsize+j];
      else if(rowMap[j]>=0) value=-vec[(size_t)rowMap[j]*Nsize+i];
      else value=slt[(size_t)(eleIdx[i]+(i/Ne)*Nsite)*Nsite2+eleIdx[j]+(j/Ne)*Nsite];
      if(!(isfinite(creal(value)) && isfinite(cimag(value)))) status=BF_PF_NONFINITE;
      if(i/Ne==j/Ne && value!=0.0) status=BF_PF_INVALID_ARGUMENT;
      a[(size_t)i*Nsize+j]=-value; /* column-major A */
      a[(size_t)j*Nsize+i]=value;
    }
  }
  free(rowMap);
  return status;
}

static int BFStableFactor_fcmp(double complex *a, double complex *pf, double complex *inverse,
    int *pivot, double complex *work, int lwork, double *rwork, int globalQp, int stage) {
  int dim=Nsize, block=Ne, info=0, status=BF_PF_OK;
  char uplo='U', method='P';
  double complex *f=NULL;
  if(inverse) {
    f=malloc((size_t)block*block*sizeof(*f));
    if(!f) return BF_PF_INVALID_ARGUMENT;
    for(int j=0;j<block;j++) for(int i=0;i<block;i++)
      f[(size_t)j*block+i]=a[(size_t)(block+j)*dim+i];
  }
  M_ZSKPFA(&uplo,&method,&dim,a,&dim,pf,pivot,work,&lwork, rwork,&info);
  if(stage) {
    info=BFInjectNumericalInfo(stage,globalQp,info);
    if(BFInjectNumericalNan(stage,globalQp)) *pf=NAN;
  }
  status=BFStatusFromLapack(info);
  const double complex value=*pf;
  if(status==BF_PF_OK && !(isfinite(creal(value)) && isfinite(cimag(value)))) status=BF_PF_NONFINITE;
  if(status==BF_PF_OK && inverse) {
    M_ZGETRF(&block,&block,f,&block,pivot,&info);
    if(stage==BF_FAIL_ACCEPT_PF) info=BFInjectNumericalInfo(BF_FAIL_ACCEPT_LU,globalQp,info);
    status=BFStatusFromLapack(info);
    if(status==BF_PF_OK) {
      M_ZGETRI(&block,f,&block,pivot,work,&lwork,&info);
      if(stage==BF_FAIL_ACCEPT_PF) info=BFInjectNumericalInfo(BF_FAIL_ACCEPT_INVERSE,globalQp,info);
      status=BFStatusFromLapack(info);
    }
    if(status==BF_PF_OK) {
      memset(inverse,0,(size_t)dim*dim*sizeof(*inverse));
      for(int j=0;j<block;j++) for(int i=0;i<block;i++) {
        const double complex value=f[(size_t)j*block+i];
        if(!(isfinite(creal(value)) && isfinite(cimag(value)))) status=BF_PF_NONFINITE;
        inverse[(size_t)(block+i)*dim+j]=value;
        inverse[(size_t)j*dim+block+i]=-value;
      }
    }
  }
  free(f);
  return status;
}

static int BFStableEvaluate_fcmp(const double complex *slater, const int *eleIdx,
    int globalQp, int n, const int *msa, const double complex *vec,
    double complex *pf, double complex *inverse, int stage) {
  if(Nsize<=0 || LapackLWork<Nsize || !pf
     || (size_t)Nsize>SIZE_MAX/(size_t)Nsize/sizeof(double complex)
     || (size_t)LapackLWork>SIZE_MAX/sizeof(double complex)) return BF_PF_INVALID_ARGUMENT;
  double complex *a=malloc((size_t)Nsize*Nsize*sizeof(*a));
  double complex *work=malloc((size_t)LapackLWork*sizeof(*work));
  int *pivot=malloc((size_t)Nsize*sizeof(*pivot));
  double *rwork=malloc((size_t)LapackLWork*sizeof(*rwork));
  int status=BF_PF_INVALID_ARGUMENT;
  if(a && work && pivot && rwork) {
    status=BFStableBuild_fcmp(slater,eleIdx,globalQp,n,msa,vec,a);
    if(status==BF_PF_OK) status=BFStableFactor_fcmp(a,pf,inverse,pivot,work,LapackLWork, rwork,globalQp,stage);
  }
  free(a); free(work); free(pivot); free(rwork);
  return status;
}
