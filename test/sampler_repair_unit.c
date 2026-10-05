/* Exercise the production repair state machine with a deterministic factorizer.
 * The solver MPI smoke separately covers real factorization and wiring. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MPI_COMM_WORLD 0
#define MPI_INT 0
static int MPI_Abort(int comm, int code) { (void)comm; exit(code); }
static int MPI_Comm_rank(int comm, int *rank) { (void)comm; *rank=0; return 0; }
static int MPI_Bcast(void *p,int n,int type,int root,int comm) {
  (void)p;(void)n;(void)type;(void)root;(void)comm;return 0;
}
static int NVMCCalMode=1, NLanczosMode=1, NLanczosStep=1, NLanczosEstimatorMode=0;
static int AllComplexFlag=0, iFlgOrbitalGeneral=0, NBackFlowIdx=0;
static int FlagGrandCanonical=0, FlagRBM=0, NSplitSize=1, NExUpdatePath=1;
static int NQPFull=3,Nsize=2,LapackLWork=8;
static double pf[3],inverse[12], *PfM_real=pf,*InvM_real=inverse;
static int factor_calls, factor_fail;
static double refreshed_value=0.08;
static int calculateMAll_child_real(const int *idx,int start,int end,int q,
    double *buf,int *iw,double *w,int nw,double *p,double *inv) {
  (void)idx;(void)buf;(void)iw;(void)w;(void)nw;
  assert(start==0 && end==3);factor_calls++;
  if(factor_fail)return 1;
  p[q]=refreshed_value;inv[4*q]=123;return 0;
}
#include "../src/mVMC/sampler_repair.c"
static void accept(double q1) {
  double proposal[3]={1,q1,1e-9};
  SamplerRepairAcceptPre(proposal);memcpy(pf,proposal,sizeof pf);
  SamplerRepairAcceptPost(NULL,0,3);
}
int main(int argc,char **argv) {
  setenv("MVMC_SAMPLER_DRIFT_LOG","0",1);
  unsetenv("MVMC_SAMPLER_REPAIR");
  if(argc>1 && !strcmp(argv[1],"invalid")) {
    setenv("MVMC_SAMPLER_REPAIR","junk",1);SamplerRepairInit();return 0;
  }
  if(argc>1 && !strcmp(argv[1],"unsupported")) {
    NBackFlowIdx=1;setenv("MVMC_SAMPLER_REPAIR","1",1);SamplerRepairInit();return 0;
  }
  SamplerRepairInit();assert(Sr.enabled && !Sr.logging);
  pf[0]=1;pf[1]=1e-5;pf[2]=1e-9;SamplerRepairRecompute();
  accept(0.1);assert(factor_calls==1 && inverse[4]==123);
  assert(Sr.reference[1]==0.1 && Sr.floor[1]==refreshed_value && Sr.total==1);
  SamplerRepairFinalize();assert(!Sr.enabled && Sr.reference==NULL);
  // An initially small component that fell, then grew gradually, selects min
  // history even when neither last-refresh nor single-step growth reaches G.
  factor_calls=0;SamplerRepairInit();pf[0]=1;pf[1]=0.01;pf[2]=1e-9;
  SamplerRepairRecompute();accept(1e-9);accept(1e-7);
  assert(factor_calls==0 && Sr.floor[1]==1e-9);
  accept(2e-6);assert(factor_calls==1 && Sr.floor[1]==refreshed_value);
  pf[1]=0.02;SamplerRepairRecompute();assert(Sr.reference[1]==0.02 && Sr.floor[1]==0.02);
  if(argc>1 && !strcmp(argv[1],"factor-failure")) {
    factor_fail=1;accept(0.000001);accept(0.1);return 0;
  }
  // A zero reference is rebuilt, not silently clipped or excluded.
  pf[1]=0;SamplerRepairRecompute();accept(0);assert(factor_calls==2);
  // Logging summarizes the existing measurement and does not alter state.
  Sr.logging=1;SamplerRepairSetBin(0);
  SamplerRepairMeasure(1,0,0);SamplerRepairMeasure(NAN,0,0);
  assert(Sr.samples==2 && Sr.nonfinite==1 && Sr.gt1==1 && Sr.max_d==2);
  SamplerRepairSetBin(1);assert(Sr.samples==0 && Sr.max_d==0);
  SamplerRepairFinalize();
  // Explicit off retains the sampler; unsupported paths are auto-disabled.
  setenv("MVMC_SAMPLER_REPAIR","0",1);SamplerRepairInit();assert(!Sr.enabled && !Sr.reference);
  SamplerRepairAcceptPre(NULL);assert(!SamplerRepairAcceptPost(NULL,0,3));SamplerRepairFinalize();
  unsetenv("MVMC_SAMPLER_REPAIR");NBackFlowIdx=1;SamplerRepairInit();assert(!Sr.enabled);SamplerRepairFinalize();
  NBackFlowIdx=0;NLanczosMode=0;SamplerRepairInit();assert(!Sr.enabled);SamplerRepairFinalize();
  // Counter-only volume grows with bins; detail has a per-kind run budget.
  NLanczosMode=1;SamplerRepairInit();Sr.logging=1;Sr.file=tmpfile();assert(Sr.file);
  SamplerRepairSetBin(0);
  for(int i=0;i<40;i++)SamplerRepairProposal(1,i,1,0,-INFINITY,0,0);
  SamplerRepairProposal(2,0,1,0,-1000,0,0);
  SamplerRepairProposal(2,1,1,0,1000,0,INFINITY);
  for(int i=0;i<40;i++)SamplerRepairCurrent(3,i,1,0,-INFINITY,"continue-after-repair");
  for(int i=0;i<40;i++) {
    SamplerRepairMeasure(-INFINITY,0,0);
    SamplerRepairMeasureFinish(i,0,NAN,2);
  }
  // Each factorization receives a zero current component; record actual output.
  refreshed_value=0.08;
  for(int i=0;i<40;i++) {
    pf[0]=1;pf[1]=0;pf[2]=1e-9;SamplerRepairRecompute();accept(0);
  }
  assert(Sr.proposal_zero==40 && Sr.ratio_underflow==1 && Sr.ratio_nonfinite==1);
  assert(Sr.current_nonfinite==40 && Sr.measurement_zero==40);
  assert(Sr.measure_weight_zero==40 && Sr.skip_energy==40);
  assert(Sr.component_zero_pre==40 && Sr.component_zero_post==0);
  assert(Sr.events_written==16 && Sr.events_suppressed==146);
  SamplerRepairSetBin(1);assert(Sr.events_written==16 && Sr.proposal_zero==0);
  SamplerRepairProposal(4,0,1,0,-INFINITY,0,0);assert(Sr.events_written==16);
  SamplerRepairFlush();
  fflush(Sr.file);long bytes=ftell(Sr.file);assert(bytes>0 && bytes<8192);
  rewind(Sr.file);char line[512];int events=0,counts=0;
  while(fgets(line,sizeof line,Sr.file)) {
    if(line[0]=='E') {events++;assert(strchr(line,'\n'));}
    if(line[0]=='Z')counts++;
  }
  assert(events==16 && counts==2);
  printf("bounded log: %ld bytes, %d details, %ld suppressed\n",bytes,events,Sr.events_suppressed);
  fseek(Sr.file,0,SEEK_END);
  SamplerRepairFinalize();
  puts("sampler repair unit PASS");return 0;
}
