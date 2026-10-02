#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _mpi_use
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif
#include "matrix.h"
#include "blas_externs.h"

static int Nsize, Ne, Nsite, Nsite2, NQPFull=1, LapackLWork=1024;
enum { BF_FAIL_PROPOSAL=1, BF_FAIL_ACCEPT_LU, BF_FAIL_ACCEPT_INVERSE, BF_FAIL_ACCEPT_PF };
#define BFInjectNumericalInfo(stage,qp,info) (info)
#define BFInjectNumericalNan(stage,qp) 0
static int BFStatusFromLapack(int info) {
  return info<0 ? BF_PF_INVALID_ARGUMENT : (info>0 ? BF_PF_LAPACK_FAILURE : BF_PF_OK);
}
#include "../src/mVMC/backflow_stable.c"

#define CHECK(c) do { if(!(c)) { fprintf(stderr,"line %d: %s\n",__LINE__,#c); return 1; } } while(0)

/* Orthogonal Hadamard factors give independently known determinant and inverse.
 * kappa(F)=2^exponent, including the no-op counterexample's ill-conditioned range.
 * Complex row/column phases exercise transpose (not conjugate transpose). */
static int exercise(int block, int exponent, int complexMode) {
  double complex a[256]={0}, d[256], expected[256]={0}, pf, initial, vec[256];
  double ar[256], dr[256], vr[256], pr;
  double complex factorC[64], workC[1024];
  double factorR[64], workR[1024];
  int iworkC[16], iworkR[16];
  BFStableWorkspaceFcmp scratchC={factorC,workC,iworkC,1024};
  BFStableWorkspaceReal scratchR={factorR,workR,iworkR,1024};
  int labels[16], rows[16];
  Nsize=Nsite2=2*block; Ne=Nsite=block;
  double complex determinant=1.0;
  for(int k=0;k<block;k++) determinant*=ldexp(1.0,-exponent*k/(block-1));
  for(int i=0;i<block;i++) {
    double complex phase=complexMode ? cexp(I*(0.13*i+0.07*i)) : 1.0;
    determinant*=phase;
    for(int j=0;j<block;j++) {
      double f=0.0, inverse=0.0;
      for(int k=0;k<block;k++) {
        int bits=(i&k)^(j&k), parity=0;
        while(bits) { parity^=bits&1; bits>>=1; }
        double sign=parity ? -1.0 : 1.0;
        f+=sign*ldexp(1.0,-exponent*k/(block-1))/block;
        inverse+=sign*ldexp(1.0,exponent*k/(block-1))/block;
      }
      double complex p=complexMode ? cexp(I*(0.13*i+0.07*j)) : 1.0;
      a[i*Nsize+block+j]=f*p;
      a[(block+j)*Nsize+i]=-f*p;
      /* F^-1[j,i] carries the reciprocal phases of F[i,j]. */
      expected[(block+j)*Nsize+i]=inverse/p;
      expected[i*Nsize+block+j]=-inverse/p;
    }
  }
  if((block*(block-1)/2)%2) determinant=-determinant;
  for(int i=0;i<Nsize;i++) { labels[i]=i%block; rows[i]=Nsize-1-i; }
  for(int i=0;i<Nsize*Nsize;i++) ar[i]=creal(a[i]);
  double maxinv=0.0;
  for(int i=0;i<Nsize*Nsize;i++) maxinv=fmax(maxinv,cabs(expected[i]));
  for(int iteration=0;iteration<32;iteration++) {
    int count=iteration%Nsize+1;
    for(int k=0;k<count;k++) for(int j=0;j<Nsize;j++) {
      vec[k*Nsize+j]=a[rows[k]*Nsize+j]; vr[k*Nsize+j]=creal(vec[k*Nsize+j]);
    }
    if(complexMode) {
      CHECK(BFStableEvaluate_fcmp_workspace(a,labels,0,count,rows,vec,
            &pf,d,0,&scratchC)==0);
    } else {
      CHECK(BFStableEvaluate_real_workspace(ar,labels,0,count,rows,vr,
            &pr,dr,0,&scratchR)==0);
      pf=pr;
      for(int i=0;i<Nsize*Nsize;i++) d[i]=dr[i];
    }
    if(iteration==0) initial=pf;
    CHECK(pf==initial); /* no-op ratio is exactly one, with no chain drift */
    CHECK(cabs(pf/determinant-1.0)<2e-8);
    for(int i=0;i<Nsize;i++) for(int j=0;j<Nsize;j++) {
      CHECK(d[i*Nsize+j]==-d[j*Nsize+i]);
      CHECK(cabs(d[i*Nsize+j]-expected[i*Nsize+j])/maxinv<2e-8);
      double complex residual=-(i==j);
      for(int k=0;k<Nsize;k++) residual+=a[i*Nsize+k]*d[k*Nsize+j];
      CHECK(cabs(residual)/maxinv<2e-14);
    }
  }
  /* A real change of two rows, including their intersection. A diagonal
   * congruence has a known Pfaffian multiplier and inverse transformation. */
  rows[0]=1; rows[1]=block;
  for(int k=0;k<2;k++) for(int j=0;j<Nsize;j++) {
    double scale=(k==0 ? 1.25 : 0.75)*(j==rows[0] ? 1.25 : j==rows[1] ? 0.75 : 1.0);
    vec[k*Nsize+j]=scale*a[rows[k]*Nsize+j];vr[k*Nsize+j]=creal(vec[k*Nsize+j]);
  }
  if(complexMode) CHECK(BFStableEvaluate_fcmp(a,labels,0,2,rows,vec,&pf,d,0)==0);
  else { CHECK(BFStableEvaluate_real(ar,labels,0,2,rows,vr,&pr,dr,0)==0); pf=pr; }
  CHECK(cabs(pf/determinant-1.25*0.75)<2e-8);
  rows[1]=rows[0];
  CHECK(BFStableEvaluate_real(ar,labels,0,2,rows,vr,&pr,dr,0)==BF_PF_INVALID_ARGUMENT);
  rows[1]=Nsize;
  CHECK(BFStableEvaluate_real(ar,labels,0,2,rows,vr,&pr,dr,0)==BF_PF_INVALID_ARGUMENT);
  CHECK(BFStableEvaluate_real(ar,labels,0,-1,rows,vr,&pr,dr,0)==BF_PF_INVALID_ARGUMENT);
  return 0;
}

/* Cover every Pfaffian prefactor parity and force LU row pivots with a known
 * row permutation of a lower triangular F. */
static int exercise_pivot_sign(int block) {
  const int dim=2*block;
  double complex *slater=calloc((size_t)dim*dim,sizeof(*slater));
  double *slaterR=calloc((size_t)dim*dim,sizeof(*slaterR));
  double complex *factorC=malloc((size_t)block*block*sizeof(*factorC));
  double *factorR=malloc((size_t)block*block*sizeof(*factorR));
  double complex workC[64], pfC;
  double workR[64], pfR;
  int *labels=malloc((size_t)dim*sizeof(*labels));
  int *iworkC=malloc((size_t)dim*sizeof(*iworkC));
  int *iworkR=malloc((size_t)dim*sizeof(*iworkR));
  CHECK(slater && slaterR && factorC && factorR && labels && iworkC && iworkR);
  Nsize=Nsite2=dim;Ne=Nsite=block;
  double determinant=1.0;
  for(int i=0;i<block;i++) determinant*=1.0+0.125*i;
  if(block>1) determinant=-determinant;
  if((block&3)==2 || (block&3)==3) determinant=-determinant;
  for(int i=0;i<dim;i++) labels[i]=i%block;
  for(int i=0;i<block;i++) for(int j=0;j<block;j++) {
    int source=i;
    if(block>1) source=(i==0 ? block-1 : (i==block-1 ? 0 : i));
    const double f=source<j ? 0.0 : (source==j ? 1.0+0.125*j
        : 0.03125*(1+source+j));
    slater[(size_t)(block+j)*dim+i]=-f;
    slaterR[(size_t)(block+j)*dim+i]=-f;
  }
  BFStableWorkspaceFcmp scratchC={factorC,workC,iworkC,64};
  BFStableWorkspaceReal scratchR={factorR,workR,iworkR,64};
  CHECK(BFStableEvaluate_fcmp_workspace(slater,labels,0,0,NULL,NULL,
        &pfC,NULL,0,&scratchC)==BF_PF_OK);
  CHECK(BFStableEvaluate_real_workspace(slaterR,labels,0,0,NULL,NULL,
        &pfR,NULL,0,&scratchR)==BF_PF_OK);
  CHECK(cabs(pfC/determinant-1.0)<2e-14);
  CHECK(fabs(pfR/determinant-1.0)<2e-14);
  free(slater);free(slaterR);free(factorC);free(factorR);
  free(labels);free(iworkC);free(iworkR);
  return 0;
}

int main(void) {
  for(int block=2;block<=8;block*=2) for(int e=0;e<=24;e+=8) for(int z=0;z<2;z++)
    CHECK(exercise(block,e,z)==0);
  for(int block=1;block<=9;block++) CHECK(exercise_pivot_sign(block)==0);
  Nsize=Nsite2=4;Ne=Nsite=2;
  double zero[16]={0}, pf, inv[16];int labels[4]={0,1,0,1};
  CHECK(BFStableEvaluate_real(zero,labels,0,0,NULL,NULL,&pf,NULL,0)==0 && pf==0.0);
  CHECK(BFStableEvaluate_real(zero,labels,0,0,NULL,NULL,&pf,inv,0)==BF_PF_LAPACK_FAILURE);
  zero[8]=NAN;
  CHECK(BFStableEvaluate_real(zero,labels,0,0,NULL,NULL,&pf,NULL,0)==BF_PF_NONFINITE);
  zero[8]=0.0;zero[4]=1.0;
  CHECK(BFStableEvaluate_real(zero,labels,0,0,NULL,NULL,&pf,NULL,0)==BF_PF_INVALID_ARGUMENT);
  zero[4]=NAN;
  CHECK(BFStableEvaluate_real(zero,labels,0,0,NULL,NULL,&pf,NULL,0)==BF_PF_INVALID_ARGUMENT);
  puts("BackFlow stable occupied kernels: real/complex, no-op chains, two-row changes, inverse, invalid and singular cases passed");
  return 0;
}
