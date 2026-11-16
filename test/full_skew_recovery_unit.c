/* Exercise the production matrix wrappers. The spies delegate to the actual
 * PFAPACK/ltl2inv routines except in explicitly labelled error-contract tests. */
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _mpi_use
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif
extern int omp_get_thread_num(void);
extern int omp_get_max_threads(void);
#include "global.h"
#include "blas_externs.h"

static int fault, inversions, last_info, failures, cases, rank_id, rank_count=1;
static int snapshot;
static int spy_dsktrf(const char *u,const char *m,const int *n,double *a,
    const int *lda,int *p,double *w,const int *lw,int *info) {
  int rc=M_DSKTRF(u,m,n,a,lda,p,w,lw,info);
  last_info=*info;
  if(fault==1) *info=-4;
  if(fault>=2 && fault<=4) { *info=2; a[*lda]=fault==2?0:fault==3?NAN:INFINITY; }
  return rc;
}
static int spy_zsktrf(const char *u,const char *m,const int *n,double complex *a,
    const int *lda,int *p,double complex *w,const int *lw,int *info) {
  int rc=M_ZSKTRF(u,m,n,a,lda,p,w,lw,info);
  last_info=*info;
  if(fault==1) *info=-4;
  if(fault>=2 && fault<=4) { *info=2; a[*lda]=fault==2?0:fault==3?NAN:INFINITY; }
  if(fault==8) {
    double parts[2]={1,NAN}; *info=2; memcpy(&a[*lda],parts,sizeof(parts));
  }
  return rc;
}
static void spy_inv_d(int n,double *a,int ld,int *p,double *v,double *w,int lw) {
  ++inversions; utu2inv_d(n,a,ld,p,v,w,lw);
  if(fault==5) a[0]=INFINITY;
}
static void spy_inv_z(int n,double complex *a,int ld,int *p,double complex *v,double complex *w,int lw) {
  ++inversions; utu2inv_z(n,a,ld,p,v,w,lw);
  if(fault==5) a[0]=INFINITY;
}
static void spy_pf_d(int n,double *a,int ld,int *p,double *pf) {
  utu2pfa_d(n,a,ld,p,pf); if(fault==6) *pf=INFINITY;
}
static void spy_pf_z(int n,double complex *a,int ld,int *p,double complex *pf) {
  utu2pfa_z(n,a,ld,p,pf);
  if(fault==6) *pf=INFINITY;
  if(fault==7) *pf=DBL_MAX+I*DBL_MAX;
}
#undef M_DSKTRF
#undef M_ZSKTRF
#define M_DSKTRF spy_dsktrf
#define M_ZSKTRF spy_zsktrf
#define utu2inv_d spy_inv_d
#define utu2inv_z spy_inv_z
#define utu2pfa_d spy_pf_d
#define utu2pfa_z spy_pf_z
#include "../src/mVMC/matrix.c"

#define CHECK(x,msg) do { if(!(x)) { \
  fprintf(stderr,"rank=%d n=%d fault=%d: %s (line %d)\n",rank_id,Nsize,fault,msg,__LINE__); ++failures; } } while(0)
static void *allocate(size_t n,size_t s) {
  void *p=calloc(n,s); if(!p) {fprintf(stderr,"allocation failed\n");exit(2);} return p;
}
static double complex direct_pf(const double complex *a,int n) {
  if(n==0) return 1;
  double complex sum=0,sub[64];
  for(int j=1;j<n;j++) {
    int r=0;
    for(int i=1;i<n;i++) if(i!=j) {
      int c=0;
      for(int k=1;k<n;k++) if(k!=j) sub[r*(n-2)+c++]=a[i*n+k];
      ++r;
    }
    sum+=(j%2?1:-1)*a[j]*direct_pf(sub,n-2);
  }
  return sum;
}
/* Analytic block/tridiagonal Pfaffians, with an optional odd permutation. */
static double complex make_matrix(double complex *a,int n,int kind,int cx,double scale) {
  double complex expected=1;
  for(int i=0;i<n;i+=2) {
    double complex v=scale*(1.5+(i/2)%3)*(cx?1+0.25*I:1);
    if(kind==3 && i==0) v=0;
    if(kind==4 && i==2*(n/4)) v=0;
    if(kind==5 && i==n-2) v=0;
    if(kind==6) v=0;
    a[i*n+i+1]=v; a[(i+1)*n+i]=-v; expected*=v;
  }
  if(kind==1) for(int i=1;i<n-1;i+=2) {
    double complex v=scale*(cx?0.125+0.0625*I:0.125);
    a[i*n+i+1]=v; a[(i+1)*n+i]=-v;
  }
  if(kind==2) {
    double complex *copy=allocate(n*n,sizeof(*copy));memcpy(copy,a,n*n*sizeof(*a));
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) a[i*n+j]=copy[((i+1)%n)*n+(j+1)%n];
    expected=-expected; free(copy);
  }
  return expected;
}
static int rebuild(int cx,int fsz,int *idx,int *spn,int start,int end) {
  return cx ? (fsz?CalculateMAll_fsz(idx,spn,start,end):CalculateMAll_fcmp(idx,start,end))
            : (fsz?CalculateMAll_fsz_real(idx,spn,start,end):CalculateMAll_real(idx,start,end));
}
static void run_case(int n,int kind,int cx,int fsz,double scale,int expect_ok) {
  const int start=2+2*rank_id, end=start+2, nsq=n*n;
  int *idx=allocate(n,sizeof(*idx)), *spn=allocate(n,sizeof(*spn));
  double complex *a=allocate(nsq,sizeof(*a));
  Nsize=n;Ne=n/2;Nsite=n+1;Nsite2=2*Nsite;NQPFull=2+2*rank_count;
  LapackLWork=n*n+64*n;
  const int stride=Nsite2*Nsite2;
  SlaterElm=allocate((size_t)NQPFull*stride,sizeof(*SlaterElm));
  SlaterElm_real=allocate((size_t)NQPFull*stride,sizeof(*SlaterElm_real));
  InvM=allocate(3*nsq,sizeof(*InvM));InvM_real=allocate(3*nsq,sizeof(*InvM_real));
  PfM=allocate(3,sizeof(*PfM));PfM_real=allocate(3,sizeof(*PfM_real));
  for(int i=0;i<n;i++) {idx[i]=(i+1)%Nsite;spn[i]=fsz?i%2:i/Ne;}
  double complex expected=make_matrix(a,n,kind>=7?0:kind,cx,scale);
  if(n<=8 && fabs(scale)<1e10) {
    double complex oracle=direct_pf(a,n);
    CHECK(cabs(oracle-expected)<=1e-12*fmax(1,cabs(expected)),"independent Pfaffian oracle");
  }
  for(int qp=0;qp<NQPFull;qp++) for(int i=0;i<n;i++) for(int j=0;j<n;j++) {
    int off=qp*stride+(idx[i]+spn[i]*Nsite)*Nsite2+idx[j]+spn[j]*Nsite;
    double complex value=(1+0.125*qp)*a[i*n+j];
    if(((kind==7 && qp==start) || (kind==8 && qp==start+1)) &&
       ((i==0 && j==1) || (i==1 && j==0))) value=0;
    SlaterElm[off]=value;SlaterElm_real[off]=creal(value);
  }
  for(int i=0;i<3;i++) PfM[i]=PfM_real[i]=1234567;
  for(int i=0;i<3*nsq;i++) InvM[i]=InvM_real[i]=7654321;
  inversions=0;last_info=0;
  int status=rebuild(cx,fsz,idx,spn,start,end);
  if(n==6 && !fault && (kind==0 || kind==3 || kind==4 || kind==5))
    CHECK(last_info==(kind==5?5:4),"U-order hidden singular INFO vector");
  ++cases;
  if(expect_ok) {
    CHECK(status==0,"nonsingular wrapper succeeds");
    CHECK(inversions==2,"both rank-local QPs inverted");
    for(int qp=0;qp<2 && status==0;qp++) {
      double complex pf=cx?PfM[qp]:PfM_real[qp];
      double qpScale=1+0.125*(start+qp);
      double complex qpExpected=expected*pow(qpScale,n/2);
      if(fault==7) CHECK(creal(pf)==DBL_MAX && cimag(pf)==DBL_MAX,"finite components whose sum overflows");
      else CHECK(cabs(pf-qpExpected)<=1e-10*fmax(DBL_MIN,cabs(qpExpected)),"Pfaffian sign/phase/value");
      double residual=0,skew=0;
      for(int i=0;i<n;i++) for(int j=0;j<n;j++) {
        double complex left=0,right=0;
        for(int k=0;k<n;k++) {
          left+=qpScale*a[i*n+k]*(cx?InvM[qp*nsq+k*n+j]:InvM_real[qp*nsq+k*n+j]);
          right+=(cx?InvM[qp*nsq+i*n+k]:InvM_real[qp*nsq+i*n+k])*qpScale*a[k*n+j];
        }
        CHECK(isfinite(creal(left)) && isfinite(cimag(left)) &&
              isfinite(creal(right)) && isfinite(cimag(right)),"finite matrix products");
        residual=fmax(residual,fmax(cabs(left-(i==j)),cabs(right-(i==j))));
        double complex v=cx?InvM[qp*nsq+i*n+j]:InvM_real[qp*nsq+i*n+j];
        double complex t=cx?InvM[qp*nsq+j*n+i]:InvM_real[qp*nsq+j*n+i];
        CHECK(isfinite(creal(v)) && isfinite(cimag(v)),"finite inverse");
        skew=fmax(skew,cabs(v+t));
      }
      CHECK(residual<1e-10,"left/right inverse with mVMC row-major sign");
      CHECK(skew<1e-10,"skew inverse");
      if(snapshot) {
        printf("n%d cx%d fsz%d qp%d pf=%a,%a",n,cx,fsz,qp,creal(pf),cimag(pf));
        for(int i=0;i<nsq;i++) {
          double complex v=cx?InvM[qp*nsq+i]:InvM_real[qp*nsq+i];
          printf(" %a,%a",creal(v),cimag(v));
        }
        puts("");
      }
    }
    if(kind==1 && !fault) CHECK(last_info==0,"normal-path fixture has INFO=0");
  } else {
    CHECK(fault==1?status==-4:status>0,"failure status is preserved or positive");
    if(kind==8) {
      double complex first=cx?PfM[0]:PfM_real[0];
      double complex firstExpected=expected*pow(1+0.125*start,n/2);
      CHECK(cabs(first-firstExpected)<=1e-10*fmax(DBL_MIN,cabs(firstExpected)),"regular QP before singular QP");
      CHECK(PfM[1]==1234567 && PfM_real[1]==1234567,"failed QP does not publish Pfaffian");
      CHECK(inversions==1,"only the regular QP is inverted");
    } else CHECK(PfM[0]==1234567 && PfM_real[0]==1234567,"failure does not publish Pfaffian");
    if(kind!=8 && fault!=5 && scale!=1e-310) CHECK(inversions==0,"invalid factors/Pf do not reach inverse");
  }
  CHECK(PfM[2]==1234567 && PfM_real[2]==1234567,"QP output boundary");
  for(int i=2*nsq;i<3*nsq;i++) CHECK(InvM[i]==7654321 && InvM_real[i]==7654321,"inverse output boundary");
  CHECK(rebuild(cx,fsz,idx,spn,start,start)==0,"empty QP slice succeeds");
  free(a);free(idx);free(spn);free(SlaterElm);free(SlaterElm_real);free(InvM);free(InvM_real);free(PfM);free(PfM_real);
}
int main(int argc,char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&rank_id);MPI_Comm_size(MPI_COMM_WORLD,&rank_count);
#endif
  snapshot=argc>1 && strcmp(argv[1],"--normal-snapshot")==0;
  NThread=omp_get_max_threads();initializeWorkSpaceAll();
  const int sizes[]={2,4,6,8,130};
  for(int cx=0;cx<2;cx++) for(int fsz=0;fsz<2;fsz++) {
    for(size_t s=0;s<sizeof(sizes)/sizeof(sizes[0]);s++) {
      for(int kind=0;kind<9;kind++) if(!snapshot || kind==1)
        run_case(sizes[s],kind,cx,fsz,1,kind<3);
    }
    if(!snapshot) {
      run_case(4,0,cx,fsz,1e-12,1);
      run_case(4,0,cx,fsz,1e-180,1); /* finite inverse; Pf underflow is not a new cutoff */
      run_case(4,0,cx,fsz,1e200,0); /* finite factors, nonfinite Pf */
      run_case(4,0,cx,fsz,1e-310,0); /* finite nonzero pivots, nonfinite inverse */
      for(fault=1;fault<=8;fault++) if(cx || (fault!=7 && fault!=8))
        run_case(4,0,cx,fsz,1,fault==7);
      fault=0;
    }
  }
  FreeWorkSpaceAll();
  int total=failures;
#ifdef _mpi_use
  MPI_Allreduce(&failures,&total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  if(rank_id==0) fprintf(stderr,"full skew recovery: cases/rank=%d ranks=%d failures=%d\n",cases,rank_count,total);
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return total?1:0;
}
