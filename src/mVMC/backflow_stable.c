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
#include "include/backflow_stable.h"

static int BFStableBuild_real(const double *slater, const int *eleIdx,
    int globalQp, int n, const int *msa, const double *vec, double *f,
    int *rowMap) {
  if(!slater || !eleIdx || !f || !rowMap
     || Nsize <= 0 || Ne <= 0 || Nsize != 2LL*Ne
     || Nsite <= 0 || Nsite2 != 2LL*Nsite || globalQp < 0 || globalQp >= NQPFull
     || n < 0 || n > Nsize || (n && (!msa || !vec))) return BF_PF_INVALID_ARGUMENT;
  for(int i=0;i<Nsize;i++) {
    if(eleIdx[i] < 0 || eleIdx[i] >= Nsite) return BF_PF_INVALID_ARGUMENT;
    rowMap[i]=-1;
  }
  for(int k=0;k<n;k++) {
    if(msa[k]<0 || msa[k]>=Nsize || rowMap[msa[k]]!=-1)
      return BF_PF_INVALID_ARGUMENT;
    rowMap[msa[k]]=k;
  }
  const double *slt=slater+(size_t)globalQp*Nsite2*Nsite2;
  int status=BF_PF_OK;
  /* Preserve the exact zero-block contract even though only F is factored. */
  for(int spin=0;spin<2;spin++) for(int i=1;i<Ne;i++) for(int j=0;j<i;j++) {
      const int ii=spin*Ne+i, jj=spin*Ne+j;
      /* Unchanged entries use the base configuration; changed rows already
       * contain the candidate entries, including intersections of two rows. */
      double value;
      if(rowMap[ii]>=0) value=vec[(size_t)rowMap[ii]*Nsize+jj];
      else if(rowMap[jj]>=0) value=-vec[(size_t)rowMap[jj]*Nsize+ii];
      else value=slt[(size_t)(eleIdx[ii]+spin*Nsite)*Nsite2+eleIdx[jj]+spin*Nsite];
      if(!isfinite(value)) status=BF_PF_NONFINITE;
      if(value!=0.0) status=BF_PF_INVALID_ARGUMENT;
  }
  for(int j=0;j<Ne;j++) for(int i=0;i<Ne;i++) {
    const int down=Ne+j, up=i;
    double value;
    if(rowMap[down]>=0) value=vec[(size_t)rowMap[down]*Nsize+up];
    else if(rowMap[up]>=0) value=-vec[(size_t)rowMap[up]*Nsize+down];
    else value=slt[(size_t)(eleIdx[down]+Nsite)*Nsite2+eleIdx[up]];
    if(!isfinite(value)) status=BF_PF_NONFINITE;
    f[(size_t)j*Ne+i]=-value; /* column-major F from A=[0 F;-F^T 0] */
  }
  return status;
}

static int BFStableFactor_real(double *f, double *pf, double *inverse,
    int *pivot, double *work, int lwork, int globalQp, int stage) {
  int dim=Nsize, block=Ne, info=0, status=BF_PF_OK;
  M_DGETRF(&block,&block,f,&block,pivot,&info);
  if(info<0) status=BFStatusFromLapack(info);
  if(status==BF_PF_OK && info>0) *pf=0.0;
  if(status==BF_PF_OK && info==0) {
    double value=((block&3)==2 || (block&3)==3) ? -1.0 : 1.0;
    for(int i=0;i<block;i++) {
      if(pivot[i]!=i+1) value=-value;
      value*=f[(size_t)i*block+i];
    }
    *pf=value;
  }
  if(stage) {
    const int pfInfo=BFInjectNumericalInfo(stage,globalQp,0);
    if(pfInfo) status=BFStatusFromLapack(pfInfo);
    if(BFInjectNumericalNan(stage,globalQp)) *pf=NAN;
  }
  if(status==BF_PF_OK && !isfinite(*pf)) status=BF_PF_NONFINITE;
  if(status==BF_PF_OK && inverse) {
    /* A singular F still represents a valid zero Pfaffian, but it cannot
     * supply the inverse required by a committed state. */
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
  return status;
}

static int BFStableEvaluate_real(const double *slater, const int *eleIdx,
    int globalQp, int n, const int *msa, const double *vec,
    double *pf, double *inverse, int stage) {
  if(Nsize<=0 || Ne<=0 || LapackLWork<Nsize || !pf
     || (size_t)Ne>SIZE_MAX/(size_t)Ne/sizeof(double)
     || (size_t)LapackLWork>SIZE_MAX/sizeof(double)) return BF_PF_INVALID_ARGUMENT;
  double *a=malloc((size_t)Ne*Ne*sizeof(*a));
  double *work=malloc((size_t)LapackLWork*sizeof(*work));
  int *pivot=malloc((size_t)Nsize*sizeof(*pivot));

  int status=BF_PF_INVALID_ARGUMENT;
  if(a && work && pivot) {
    status=BFStableBuild_real(slater,eleIdx,globalQp,n,msa,vec,a,pivot);
    if(status==BF_PF_OK) status=BFStableFactor_real(a,pf,inverse,pivot,work,LapackLWork,globalQp,stage);
  }
  free(a); free(work); free(pivot);
  return status;
}

static int BFStableEvaluate_real_workspace(const double *slater,
    const int *eleIdx, int globalQp, int n, const int *msa,
    const double *vec, double *pf, double *inverse, int stage,
    BFStableWorkspaceReal *scratch) {
  if(!scratch || !scratch->factor || !scratch->iwork || !pf
     || (inverse && (!scratch->work || scratch->lwork<Ne)))
    return BF_PF_INVALID_ARGUMENT;
  int status=BFStableBuild_real(slater,eleIdx,globalQp,n,msa,vec,
      scratch->factor,scratch->iwork);
  if(status==BF_PF_OK) status=BFStableFactor_real(scratch->factor,pf,inverse,
      scratch->iwork,scratch->work,scratch->lwork,globalQp,stage);
  return status;
}

static int BFStableBuild_fcmp(const double complex *slater, const int *eleIdx,
    int globalQp, int n, const int *msa, const double complex *vec,
    double complex *f, int *rowMap) {
  if(!slater || !eleIdx || !f || !rowMap
     || Nsize <= 0 || Ne <= 0 || Nsize != 2LL*Ne
     || Nsite <= 0 || Nsite2 != 2LL*Nsite || globalQp < 0 || globalQp >= NQPFull
     || n < 0 || n > Nsize || (n && (!msa || !vec))) return BF_PF_INVALID_ARGUMENT;
  for(int i=0;i<Nsize;i++) {
    if(eleIdx[i] < 0 || eleIdx[i] >= Nsite) return BF_PF_INVALID_ARGUMENT;
    rowMap[i]=-1;
  }
  for(int k=0;k<n;k++) {
    if(msa[k]<0 || msa[k]>=Nsize || rowMap[msa[k]]!=-1)
      return BF_PF_INVALID_ARGUMENT;
    rowMap[msa[k]]=k;
  }
  const double complex *slt=slater+(size_t)globalQp*Nsite2*Nsite2;
  int status=BF_PF_OK;
  for(int spin=0;spin<2;spin++) for(int i=1;i<Ne;i++) for(int j=0;j<i;j++) {
      const int ii=spin*Ne+i, jj=spin*Ne+j;
      /* Unchanged entries use the base configuration; changed rows already
       * contain the candidate entries, including intersections of two rows. */
      double complex value;
      if(rowMap[ii]>=0) value=vec[(size_t)rowMap[ii]*Nsize+jj];
      else if(rowMap[jj]>=0) value=-vec[(size_t)rowMap[jj]*Nsize+ii];
      else value=slt[(size_t)(eleIdx[ii]+spin*Nsite)*Nsite2+eleIdx[jj]+spin*Nsite];
      if(!(isfinite(creal(value)) && isfinite(cimag(value)))) status=BF_PF_NONFINITE;
      if(value!=0.0) status=BF_PF_INVALID_ARGUMENT;
  }
  for(int j=0;j<Ne;j++) for(int i=0;i<Ne;i++) {
    const int down=Ne+j, up=i;
    double complex value;
    if(rowMap[down]>=0) value=vec[(size_t)rowMap[down]*Nsize+up];
    else if(rowMap[up]>=0) value=-vec[(size_t)rowMap[up]*Nsize+down];
    else value=slt[(size_t)(eleIdx[down]+Nsite)*Nsite2+eleIdx[up]];
    if(!(isfinite(creal(value)) && isfinite(cimag(value)))) status=BF_PF_NONFINITE;
    f[(size_t)j*Ne+i]=-value;
  }
  return status;
}

static int BFStableFactor_fcmp(double complex *f, double complex *pf, double complex *inverse,
    int *pivot, double complex *work, int lwork, double *rwork, int globalQp, int stage) {
  int dim=Nsize, block=Ne, info=0, status=BF_PF_OK;
  (void)rwork;
  M_ZGETRF(&block,&block,f,&block,pivot,&info);
  if(info<0) status=BFStatusFromLapack(info);
  if(status==BF_PF_OK && info>0) *pf=0.0;
  if(status==BF_PF_OK && info==0) {
    double complex value=((block&3)==2 || (block&3)==3) ? -1.0 : 1.0;
    for(int i=0;i<block;i++) {
      if(pivot[i]!=i+1) value=-value;
      value*=f[(size_t)i*block+i];
    }
    *pf=value;
  }
  if(stage) {
    const int pfInfo=BFInjectNumericalInfo(stage,globalQp,0);
    if(pfInfo) status=BFStatusFromLapack(pfInfo);
    if(BFInjectNumericalNan(stage,globalQp)) *pf=NAN;
  }
  if(status==BF_PF_OK && !(isfinite(creal(*pf)) && isfinite(cimag(*pf)))) status=BF_PF_NONFINITE;
  if(status==BF_PF_OK && inverse) {
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
  return status;
}

static int BFStableEvaluate_fcmp(const double complex *slater, const int *eleIdx,
    int globalQp, int n, const int *msa, const double complex *vec,
    double complex *pf, double complex *inverse, int stage) {
  if(Nsize<=0 || Ne<=0 || LapackLWork<Nsize || !pf
     || (size_t)Ne>SIZE_MAX/(size_t)Ne/sizeof(double complex)
     || (size_t)LapackLWork>SIZE_MAX/sizeof(double complex)) return BF_PF_INVALID_ARGUMENT;
  double complex *a=malloc((size_t)Ne*Ne*sizeof(*a));
  double complex *work=malloc((size_t)LapackLWork*sizeof(*work));
  int *pivot=malloc((size_t)Nsize*sizeof(*pivot));
  int status=BF_PF_INVALID_ARGUMENT;
  if(a && work && pivot) {
    status=BFStableBuild_fcmp(slater,eleIdx,globalQp,n,msa,vec,a,pivot);
    if(status==BF_PF_OK) status=BFStableFactor_fcmp(a,pf,inverse,pivot,work,LapackLWork,NULL,globalQp,stage);
  }
  free(a); free(work); free(pivot);
  return status;
}

static int BFStableEvaluate_fcmp_workspace(const double complex *slater,
    const int *eleIdx, int globalQp, int n, const int *msa,
    const double complex *vec, double complex *pf, double complex *inverse,
    int stage, BFStableWorkspaceFcmp *scratch) {
  if(!scratch || !scratch->factor || !scratch->iwork || !pf
     || (inverse && (!scratch->work || scratch->lwork<Ne)))
    return BF_PF_INVALID_ARGUMENT;
  int status=BFStableBuild_fcmp(slater,eleIdx,globalQp,n,msa,vec,
      scratch->factor,scratch->iwork);
  if(status==BF_PF_OK) status=BFStableFactor_fcmp(scratch->factor,pf,inverse,
      scratch->iwork,scratch->work,scratch->lwork,NULL,globalQp,stage);
  return status;
}
