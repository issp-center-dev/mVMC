/* Testing-only measurement replay. Sampler trajectory regression is separate:
 * every measurement route receives the same saved particle labels and counts. */
#ifdef MVMC_ENABLE_FAULT_INJECTION
#include <limits.h>
static int BFReplayReadInt(FILE *fp, int *value) {
  char token[128], *end;
  long parsed;
  if(fscanf(fp,"%127s",token) != 1) return 1;
  errno=0;
  parsed=strtol(token,&end,10);
  if(errno || *end || parsed < INT_MIN || parsed > INT_MAX) return 1;
  *value=(int)parsed;
  return 0;
}

static int BFReplayMeasurementConfigurations(MPI_Comm comm) {
  int rank, enabled=0, invalid=0;
  MPI_Comm_rank(comm,&rank);
  if(rank == 0) {
    const char *path=getenv("MVMC_BF_TEST_REPLAY_CONFIG");
    const char *record=getenv("MVMC_BF_TEST_RECORD_CONFIG");
    if(path && *path && record && *record) {
      fprintf(stderr,"Error: invalid BackFlow replay: record and replay are exclusive.\n");
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    if(record && *record) {
      FILE *fp=fopen(record,"w");
      int failed=!fp;
      if(fp) {
        failed=fprintf(fp,"%d %d\n",NVMCSample,Nsize) < 0;
        for(int sample=0;sample<NVMCSample;sample++)
          for(int i=0;i<Nsize;i++)
            if(fprintf(fp,"%d%c",EleIdx[(size_t)sample*Nsize+i],i+1==Nsize?'\n':' ') < 0) failed=1;
        if(fclose(fp)) failed=1;
      }
      if(failed) {
        fprintf(stderr,"Error: failed to record BackFlow replay configurations.\n");
        MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
      }
    }
    enabled=path && *path;
    if(enabled) {
      FILE *fp=fopen(path,"r");
      int samples=0, particles=0;
      invalid=!fp || Nsize <= 0 || Ne <= 0 || Nsize != 2LL*Ne
        || NVMCSample <= 0 || NVMCSample > INT_MAX/Nsize;
      if(!invalid) invalid=BFReplayReadInt(fp,&samples) || BFReplayReadInt(fp,&particles)
        || samples != NVMCSample || particles != Nsize;
      for(int i=0;!invalid && i<NVMCSample*Nsize;i++)
        invalid=BFReplayReadInt(fp,EleIdx+i);
      if(!invalid) {
        char extra[2];
        invalid=fscanf(fp,"%1s",extra) != EOF || ferror(fp);
      }
      if(fp) fclose(fp);
    }
  }
  MPI_Bcast(&enabled,1,MPI_INT,0,comm);
  if(!enabled) return 0;
  MPI_Bcast(&invalid,1,MPI_INT,0,comm);
  if(invalid) {
    if(rank == 0) fprintf(stderr,"Error: invalid BackFlow measurement replay file.\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  MPI_Bcast(EleIdx,NVMCSample*Nsize,MPI_INT,0,comm);
  for(int sample=0;sample<NVMCSample;sample++) {
    int *idx=EleIdx+(size_t)sample*Nsize;
    int *cfg=EleCfg+(size_t)sample*Nsite2;
    int *num=EleNum+(size_t)sample*Nsite2;
    for(int i=0;i<Nsite2;i++) {cfg[i]=-1;num[i]=0;}
    for(int i=0;i<Nsize;i++) {
      if(idx[i] < 0 || idx[i] >= Nsite || i/Ne > 1) {invalid=1;break;}
      int site=idx[i]+(i/Ne)*Nsite;
      if(num[site]) {invalid=1;break;}
      cfg[site]=i%Ne;num[site]=1;
    }
    for(int site=0;site<Nsite;site++)
      if(LocSpn[site] == 1 && num[site]+num[site+Nsite] != 1) invalid=1;
    if(invalid) {
      if(rank == 0) fprintf(stderr,"Error: invalid BackFlow replay configuration at sample %d.\n",sample);
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
    }
    MakeProjCnt(EleProjCnt+(size_t)sample*NProj,num);
    MakeProjBFCnt(EleProjBFCnt+(size_t)sample*16*Nsite*Nrange,num);
  }
  if(rank == 0) fprintf(stderr,"BackFlow measurement replay: samples=%d particles=%d\n",NVMCSample,Nsize);
  return 1;
}
#else
#define BFReplayMeasurementConfigurations(comm) 0
#endif
