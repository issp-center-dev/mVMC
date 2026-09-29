/* Included by backflow_vmc.c. Non-FSZ candidate transactions share a sampler
 * communicator: status agreement always precedes a collective amplitude. */
#include "include/backflow_acceptance.h"
typedef struct {
  MPI_Comm comm;
  int real, rank, qpStart, qpEnd;
  int proposalRecovered, preparedFull;
  int started;
  int lastDecision;
  int lastAcceptRecovered;
  size_t slaterCount, invCount;
  double complex *slater, *pf, *inv;
  double *slaterR, *pfR, *invR;
  int *proposalIdx;
  long long proposalRecovery, acceptRecovery, periodicRecovery, zeroProposal;
} BFSampleTransaction;

#ifdef MVMC_ENABLE_FAULT_INJECTION
static void BFTransactionTrace(BFSampleTransaction *tx, int kind, double u,
    double complex oldLog, double complex newLog) {
  if(getenv("MVMC_BF_TEST_TRACE") == NULL) return;
  int worldRank;
  char path[80];
  MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
  snprintf(path,sizeof(path),"bf_transaction_trace.rank%d.jsonl",worldRank);
  FILE *fp=fopen(path,"a");
  if(!fp) MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  fprintf(fp,"{\"kind\":%d,\"accept\":%d,\"u\":%.17e,\"log_committed\":\"%.17e\","
             "\"log_candidate\":\"%.17e\",\"proposal_recovery\":%d,\"accept_recovery\":%d,\"qp_start\":%d,\"qp_end\":%d,",
          kind,tx->lastDecision,u,creal(oldLog),creal(newLog),tx->proposalRecovered,tx->lastAcceptRecovered,
          tx->qpStart,tx->qpEnd);
  fprintf(fp,"\"ele_idx\":[");
  for(int i=0;i<Nsize;i++) fprintf(fp,"%s%d",i?",":"",TmpEleIdx[i]);
  fprintf(fp,"],\"candidate_ele_idx\":[");
  for(int i=0;i<Nsize;i++) fprintf(fp,"%s%d",i?",":"",tx->proposalIdx[i]);
  fprintf(fp,"],\"count\":[");
  for(int i=0;i<NProj;i++) fprintf(fp,"%s%d",i?",":"",TmpEleProjCnt[i]);
  fprintf(fp,"],\"bf_count\":[");
  for(int i=0;i<16*Nsite*Nrange;i++) fprintf(fp,"%s%d",i?",":"",TmpEleProjBFCnt[i]);
  fprintf(fp,"],\"pf\":[");
  for(int q=0;q<tx->qpEnd-tx->qpStart;q++) {
    double complex value=tx->real ? PfM_real[q] : PfM[q];
    fprintf(fp,"%s[%.17e,%.17e]",q?",":"",creal(value),cimag(value));
  }
  fprintf(fp,"],\"inverse\":[");
  for(size_t i=0;i<tx->invCount;i++) {
    double complex value=tx->real ? InvM_real[i] : InvM[i];
    fprintf(fp,"%s[%.17e,%.17e]",i?",":"",creal(value),cimag(value));
  }
  fprintf(fp,"],\"slater\":[");
  for(size_t i=0;i<tx->slaterCount;i++) {
    double complex value=tx->real ? SlaterElmBF_real[i] : SlaterElmBF[i];
    fprintf(fp,"%s[%.17e,%.17e]",i?",":"",creal(value),cimag(value));
  }
  fprintf(fp,"]}\n");
  fclose(fp);
}
#else
#define BFTransactionTrace(tx,kind,u,oldLog,newLog) ((void)0)
#endif

static void BFTransactionInit(BFSampleTransaction *tx) {
#ifdef MVMC_ENABLE_FAULT_INJECTION
  unsigned masks[4]={0,0,0,0};
  int invalid=0;
  if(tx->rank == 0) {
    const char *raw=getenv("MVMC_BF_TEST_FAILURE");
    if(raw && *raw) {
      char *copy=strdup(raw), *token;
      if(!copy) invalid=1;
      else for(token=strtok(copy,",");token;token=strtok(NULL,",")) {
        int stage=0;
        if(strcmp(token,"accept_reject")==0) {masks[3]=1;continue;}
        const char *names[]={"", "proposal", "accept_lu", "accept_inverse", "accept_pf",
                             "full_pf", "full_inverse", "periodic"};
        for(int i=1;i<=BF_FAIL_PERIODIC;i++) if(strcmp(token,names[i])==0) stage=i;
        if(strcmp(token,"proposal_negative")==0) {stage=BF_FAIL_PROPOSAL;masks[1]|=1u<<stage;}
        if(strcmp(token,"proposal_nan")==0) {stage=BF_FAIL_PROPOSAL;masks[2]|=1u<<stage;}
        if(!stage) invalid=1;
        else masks[0]|=1u<<stage;
      }
      free(copy);
    }
  }
  MPI_Bcast(&invalid,1,MPI_INT,0,tx->comm);
  MPI_Bcast(masks,4,MPI_UNSIGNED,0,tx->comm);
  if(invalid) {
    if(tx->rank == 0) fprintf(stderr,"Error: invalid MVMC_BF_TEST_FAILURE.\n");
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  BFFailureMask=masks[0];BFFailureNegativeMask=masks[1];BFFailureNanMask=masks[2];
  BFFailureConsumed=0;
  BFFailureForceAccept=(int)masks[3];
#else
  (void)tx;
#endif
}

static int BFTransactionInject(BFSampleTransaction *tx, int stage, int status) {
  if(tx->qpStart == 0 && tx->qpEnd > 0 && status == BF_PF_OK)
    return BFStatusFromLapack(BFInjectNumericalInfo(stage,0,0));
  return status;
}

static int BFTransactionAgree(BFSampleTransaction *tx, int status, const char *stage) {
  int global;
  if(status < BF_PF_OK || status > BF_PF_INVALID_ARGUMENT)
    status = BF_PF_INVALID_ARGUMENT;
  if(status != BF_PF_OK)
    fprintf(stderr,"BackFlow numerical status: stage=%s rank=%d qp=[%d,%d) status=%d\n",
            stage,tx->rank,tx->qpStart,tx->qpEnd,status);
  MPI_Allreduce(&status,&global,1,MPI_INT,MPI_MAX,tx->comm);
  if(global == BF_PF_INVALID_ARGUMENT) {
    if(tx->rank == 0) fprintf(stderr,"Error: BackFlow invalid numerical arguments at %s.\n",stage);
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  return global;
}

static void BFTransactionRequire(BFSampleTransaction *tx, int status, const char *stage) {
  if(BFTransactionAgree(tx,status,stage) != BF_PF_OK) {
    if(tx->rank == 0) fprintf(stderr,"Error: BackFlow checked recovery failed at %s.\n",stage);
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
}

static int BFTransactionFinitePf(BFSampleTransaction *tx) {
  for(int q=0;q<tx->qpEnd-tx->qpStart;q++) {
    if(tx->real) { if(!isfinite(tx->pfR[q])) return BF_PF_NONFINITE; }
    else if(!(isfinite(creal(tx->pf[q])) && isfinite(cimag(tx->pf[q])))) return BF_PF_NONFINITE;
  }
  return BF_PF_OK;
}

/* Pfaffian-only recovery: a singular inverse is not evidence of zero projected
 * amplitude. Empty-QP ranks participate in every subsequent agreement. */
static int BFTransactionFullPf(BFSampleTransaction *tx, const int *counts) {
  if(tx->qpStart == tx->qpEnd) return BF_PF_OK;
  return tx->real
    ? CalculateBFCanonicalPf_real(TmpEleIdx,TmpEleNum,counts,tx->qpStart,tx->qpEnd,tx->pfR)
    : CalculateBFCanonicalPf_fcmp(TmpEleIdx,TmpEleNum,counts,tx->qpStart,tx->qpEnd,tx->pf);
}

static int BFTransactionRecoveryPf(BFSampleTransaction *tx, const int *counts) {
  return BFTransactionInject(tx,BF_FAIL_FULL_PF,BFTransactionFullPf(tx,counts));
}

static double complex BFTransactionLogPf(BFSampleTransaction *tx) {
  if(tx->real) return CalculateLogIP_real(tx->pfR,tx->qpStart,tx->qpEnd,tx->comm);
  return CalculateLogIP_fcmp(tx->pf,tx->qpStart,tx->qpEnd,tx->comm);
}

static int BFTransactionValidLog(double complex value) {
  return isfinite(cimag(value)) && (isfinite(creal(value)) || creal(value) == -INFINITY);
}

static double complex BFTransactionProposal(BFSampleTransaction *tx, const int *counts, int status) {
  double complex logValue;
  tx->proposalRecovered = tx->preparedFull = 0;
  memcpy(tx->proposalIdx,TmpEleIdx,(size_t)Nsize*sizeof(int));
  if(BFUseCanonicalNonFszPath()) {
    status=BFTransactionInject(tx,BF_FAIL_PROPOSAL,status);
    if(tx->qpStart == 0 && tx->qpEnd > 0 && BFInjectNumericalNan(BF_FAIL_PROPOSAL,0)) {
      if(tx->real) tx->pfR[0]=NAN; else tx->pf[0]=NAN;
    }
  }
  if(status == BF_PF_OK) status = BFTransactionFinitePf(tx);
  status = BFTransactionAgree(tx,status,"proposal");
  if(status != BF_PF_OK) {
    tx->proposalRecovered = 1;
    tx->proposalRecovery++;
    BFTransactionRequire(tx,BFTransactionRecoveryPf(tx,counts),"proposal full Pfaffian");
    BFTransactionRequire(tx,BFTransactionFinitePf(tx),"proposal recovered values");
  }
  logValue = BFTransactionLogPf(tx);
  if(!BFTransactionValidLog(logValue)) {
    /* Even finite QP contributions may overflow their projected sum. Retry the
     * same configuration once; never silently reject an unresolved amplitude. */
    tx->proposalRecovered = 1;
    tx->proposalRecovery++;
    BFTransactionRequire(tx,BFTransactionRecoveryPf(tx,counts),"nonfinite amplitude recovery");
    BFTransactionRequire(tx,BFTransactionFinitePf(tx),"nonfinite recovered values");
    logValue = BFTransactionLogPf(tx);
    BFTransactionRequire(tx,BFTransactionValidLog(logValue) ? BF_PF_OK : BF_PF_NONFINITE,
                         "recovered log amplitude");
  }
  if(creal(logValue) == -INFINITY) tx->zeroProposal++;
  return logValue;
}

/* The long-double expression avoids overflow even when all input logarithms
 * are finite. u is drawn once by the caller and reused after recovery. */
static int BFTransactionDecision(BFSampleTransaction *tx, double projection,
    double complex newLog, double complex oldLog, double u) {
  int accept=0;
  int status=(isfinite(cimag(newLog)) && isfinite(cimag(oldLog))
      && BFLogAcceptance(projection,creal(newLog),creal(oldLog),u,&accept))
      ? BF_PF_OK : BF_PF_NONFINITE;
  BFTransactionRequire(tx,status,"log acceptance");
  return accept;
}

static int BFTransactionFullState(BFSampleTransaction *tx, const int *counts) {
  tx->preparedFull = 1;
  return tx->real
    ? RebuildSlaterMAllBF_real(TmpEleIdx,TmpEleNum,counts,tx->qpStart,tx->qpEnd,
        tx->slater,tx->pf,tx->inv,tx->slaterR,tx->pfR,tx->invR)
    : RebuildSlaterMAllBF_fcmp(TmpEleIdx,TmpEleNum,counts,tx->qpStart,tx->qpEnd,
        tx->slater,tx->pf,tx->inv);
}

static int BFTransactionAccept(BFSampleTransaction *tx, const int *counts,
    const int *icount, const int *msa, double projection, double complex oldLog,
    double complex *newLog, double u) {
  int status;
  tx->lastDecision=tx->lastAcceptRecovered=0;
  int provisional=BFTransactionDecision(tx,projection,*newLog,oldLog,u);
#ifdef MVMC_ENABLE_FAULT_INJECTION
  if(BFFailureForceAccept && !provisional && isfinite(creal(*newLog))) {
    /* Model a finite but wrong proposal value followed by an accept-stage
     * failure. Recovery must undo this provisional acceptance with the same u. */
    BFFailureForceAccept=0;
    for(int q=0;q<tx->qpEnd-tx->qpStart;q++) {
      if(tx->real) tx->pfR[q]*=exp(50.0); else tx->pf[q]*=exp(50.0);
    }
    *newLog+=50.0;
    BFFailureMask|=1u<<BF_FAIL_ACCEPT_LU;
    BFFailureConsumed&=~(1u<<BF_FAIL_ACCEPT_LU);
    provisional=BFTransactionDecision(tx,projection,*newLog,oldLog,u);
  }
#endif
  if(!provisional) return 0;
#ifdef MVMC_ENABLE_FAULT_INJECTION
  if(getenv("MVMC_BF_TEST_TRACE")) {
    for(int q=0;q<tx->qpEnd-tx->qpStart;q++) {
      double complex value=tx->real ? tx->pfR[q] : tx->pf[q];
      if(value == 0.0)
        fprintf(stderr,"BackFlow zero QP proposal: qp=%d projected_log=%.17e accept=1\n",
                tx->qpStart+q,creal(*newLog));
    }
  }
#endif
  if(BFUseCanonicalNonFszPath() || tx->proposalRecovered) {
    status = BFTransactionFullState(tx,counts);
    status = BFTransactionInject(tx,BF_FAIL_ACCEPT_LU,status);
    status = BFTransactionInject(tx,BF_FAIL_ACCEPT_INVERSE,status);
    status = BFTransactionInject(tx,BF_FAIL_ACCEPT_PF,status);
  }
  else {
    if(tx->qpEnd > tx->qpStart) AddBFProfileCounter(BFPROF_LEGACY_NONFSZ_ACCEPT,1);
    if(tx->real)
      status = UpdateMAll_BF_real(icount,msa,tx->pfR,TmpEleIdx,tx->qpStart,tx->qpEnd,tx->invR);
    else
      status = UpdateMAll_BF_fcmp(icount,msa,tx->pf,TmpEleIdx,tx->qpStart,tx->qpEnd,tx->inv);
  }
  status = BFTransactionAgree(tx,status,"accept prepare");
  if(status != BF_PF_OK) {
    tx->acceptRecovery++;
    tx->lastAcceptRecovered=1;
    BFTransactionRequire(tx,BFTransactionRecoveryPf(tx,counts),"accept full Pfaffian");
    BFTransactionRequire(tx,BFTransactionFinitePf(tx),"accept recovered values");
    *newLog = BFTransactionLogPf(tx);
    if(!BFTransactionDecision(tx,projection,*newLog,oldLog,u)) return 0;
    BFTransactionRequire(tx,BFTransactionInject(tx,BF_FAIL_FULL_INVERSE,
        BFTransactionFullState(tx,counts)),"accept full inverse");
  }
  BFTransactionRequire(tx,BFTransactionFinitePf(tx),"prepared candidate values");
  *newLog = BFTransactionLogPf(tx);
  /* PfM from prepare/recovery is authoritative. The same u decides again
   * before any live PfM/InvM or acceptance counter is committed. */
  tx->lastDecision=BFTransactionDecision(tx,projection,*newLog,oldLog,u);
  return tx->lastDecision;
}

static void BFTransactionCommit(BFSampleTransaction *tx, const int *counts) {
  const size_t np = (size_t)(tx->qpEnd-tx->qpStart);
  if(tx->preparedFull) {
    /* Also refresh the legacy eta tables associated with this configuration. */
    MakeSlaterElmBF_fcmp(TmpEleNum,counts);
    if(tx->real) memcpy(SlaterElmBF_real,tx->slaterR,tx->slaterCount*sizeof(double));
  }
  if(tx->real) {
    memcpy(PfM_real,tx->pfR,np*sizeof(double));
    memcpy(InvM_real,tx->invR,tx->invCount*sizeof(double));
    if(BFUseCanonicalNonFszPath()) {
      memcpy(PfM,tx->pf,np*sizeof(double complex));
      memcpy(InvM,tx->inv,tx->invCount*sizeof(double complex));
    }
  } else {
    memcpy(PfM,tx->pf,np*sizeof(double complex));
    memcpy(InvM,tx->inv,tx->invCount*sizeof(double complex));
  }
}

/* Preserve the normal periodic arithmetic by rebuilding from the existing
 * Slater table into temporary outputs; independently rebuild only on failure. */
static int BFTransactionCurrentState(BFSampleTransaction *tx) {
  int status=BF_PF_OK;
  int *iw=malloc((size_t)Nsize*sizeof(int));
  double *rw=malloc((size_t)LapackLWork*sizeof(double));
  double complex *cm=malloc((size_t)Nsize*Nsize*sizeof(double complex));
  double complex *cw=malloc((size_t)LapackLWork*sizeof(double complex));
  if(!iw || !rw || !cm || !cw) status=BF_PF_INVALID_ARGUMENT;
  tx->preparedFull=0;
  if(status == BF_PF_OK) for(int q=0;q<tx->qpEnd-tx->qpStart;q++) {
    int local = tx->real
      ? calculateMAll_BF_real_child(TmpEleIdx,tx->qpStart,tx->qpEnd,q,
          (double *)cm,iw,(double *)cw,LapackLWork,tx->pfR,tx->invR)
      : calculateMAll_BF_fcmp_child(TmpEleIdx,tx->qpStart,tx->qpEnd,q,
          cm,iw,cw,LapackLWork,rw,tx->pf,tx->inv);
    if(local > status) status=local;
    if(tx->real && BFUseCanonicalNonFszPath()) {
      local=calculateMAll_BF_fcmp_child(TmpEleIdx,tx->qpStart,tx->qpEnd,q,
          cm,iw,cw,LapackLWork,rw,tx->pf,tx->inv);
      if(local > status) status=local;
    }
  }
  free(iw);free(rw);free(cm);free(cw);
  return status;
}

static double complex BFTransactionRefresh(BFSampleTransaction *tx, const int *counts) {
  int local=BFTransactionCurrentState(tx);
  if(tx->started) local=BFTransactionInject(tx,BF_FAIL_PERIODIC,local);
  int status=BFTransactionAgree(tx,local,"periodic prepare");
  if(status != BF_PF_OK) {
    tx->periodicRecovery++;
    BFTransactionRequire(tx,BFTransactionInject(tx,BF_FAIL_FULL_INVERSE,
        BFTransactionFullState(tx,counts)),"periodic full inverse");
  }
  BFTransactionRequire(tx,BFTransactionFinitePf(tx),"periodic values");
  double complex value=BFTransactionLogPf(tx);
  BFTransactionRequire(tx,isfinite(creal(value)) && isfinite(cimag(value))
      ? BF_PF_OK : BF_PF_NONFINITE,"periodic log amplitude");
  BFTransactionCommit(tx,counts);
  return value;
}

static void BFTransactionFinish(void) {
#ifdef MVMC_ENABLE_FAULT_INJECTION
  /* A zero-QP rank may not have consumed a kernel hook. Do not carry sampler
   * injections into its subsequent Green-function measurement work. */
  BFFailureMask=BFFailureNegativeMask=BFFailureNanMask=0;
  BFFailureForceAccept=0;
#endif
}
