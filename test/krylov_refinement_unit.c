/* Independent analytic Pfaffian, Metropolis boundary, node-rescue, and
 * exponential tests. No threshold turns a small amplitude into zero. */
#include "krylov_refinement.h"
#include "gc_krylov_amplitude.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef _mpi_use
#include <mpi.h>
#endif

static int failures;
#define CHECK(c) do { if (!(c)) { ++failures; \
  fprintf(stderr, "refinement FAIL line %d: %s\n", __LINE__, #c); } } while (0)

typedef struct {
  double complex slater[16], weight, exponent;
  int failure;
} Fixture;

static MVMCKrylovStatus base(const uint64_t *words, size_t count,
                            void *context, MVMCRefinedBase *value) {
  Fixture *f = context;
  int occupied[4], n = 0, k;
  int64_t exponent_count = words[0] != 0;
  if (f->failure) return MVMC_KRYLOV_STATUS_NONFINITE;
  CHECK(count == 1);
  for (k = 0; k < 4; ++k) if ((words[0] >> k) & 1) occupied[n++] = k;
  return mvmc_refined_base_build(value, f->slater, 4, occupied, (size_t)n,
      &f->weight, 1, 0, &f->exponent, &exponent_count, 1);
}
static void entry(Fixture *f, int i, int j, double complex z) {
  f->slater[i*4+j] = z; f->slater[j*4+i] = -z;
}
static const MVMCKrylovFermionOperator identity_ops[4] = {
  {MVMC_KRYLOV_FERMION_CREATE,0},{MVMC_KRYLOV_FERMION_ANNIHILATE,0},
  {MVMC_KRYLOV_FERMION_ANNIHILATE,0},{MVMC_KRYLOV_FERMION_CREATE,0}};
static const MVMCKrylovHamiltonianTerm identity_terms[2] = {
  {2,0,2,0,0},{2,2,2,0,0}};
static void analytic_cases(void) {
  Fixture f = {{0}, 1.0, 0.0, 0};
  MVMCKrylovFockModel model = {0};
  MVMCKrylovFinalStatePolicy policy;
  MVMCKrylovRefinementResult result, sentinel;
  MVMCKrylovRefinementLimits limits = {16,128,128,4096,1000,512};
  double complex c[2] = {1.0,0.0}, energy;
  uint64_t vacuum = 0, full = 15, pair = 3;
  MVMCKrylovStatus status;
  model.site_count=2; model.hermitian=1;
  model.domain=MVMC_KRYLOV_DOMAIN_GC_EVEN;
  model.terms=identity_terms; model.term_count=2;
  model.operators=identity_ops; model.operator_count=4;
  CHECK(mvmc_krylov_final_state_policy_create(1,
      MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT, 19, c, 2, &policy)==0);
  /* Pf = 3*a-1, exactly -2^-54 despite a singular PFAPACK result. */
  entry(&f,0,1,0x1.5555555555555p-2);
  entry(&f,0,2,1); entry(&f,1,3,1); entry(&f,2,3,3);
  status=mvmc_krylov_refine_final(&model,&policy,&vacuum,&full,1,0.5,
      base,&f,NULL,&result);
  CHECK(status==0 && result.valid && !result.accepted &&
        !result.exact_zero_proposal);
  status=mvmc_krylov_refine_final(&model,&policy,&vacuum,&full,1,0,
      base,&f,NULL,&result);
  CHECK(status==0 && result.accepted);
  if(status==0 && result.accepted) {
    CHECK(fabs(result.evaluation.amplitude.log_abs + 54*log(2))<1e-13);
    CHECK(creal(result.evaluation.amplitude.phase)<0);
    CHECK(mvmc_krylov_final_state_local_energy_export(&result.evaluation,
                                                    &energy)==0);
    CHECK(cabs(energy-2)<1e-13);
  }
  /* Pf = +2^-53, including the previously understated error case. */
  entry(&f,0,1,0x1.5555555555556p-2);
  status=mvmc_krylov_refine_final(&model,&policy,&vacuum,&full,1,0,
      base,&f,NULL,&result);
  CHECK(status==0 && result.accepted);
  if(status==0) {
    CHECK(fabs(result.evaluation.amplitude.log_abs+53*log(2))<1e-13);
    CHECK(creal(result.evaluation.amplitude.phase)>0);
  }
  /* Exact node and exact strict boundary u == |Phi_y/Phi_x|^2. */
  memset(f.slater,0,sizeof(f.slater));
  entry(&f,0,1,1); entry(&f,0,2,1); entry(&f,1,3,1); entry(&f,2,3,1);
  status=mvmc_krylov_refine_final(&model,&policy,&vacuum,&full,1,0,
      base,&f,NULL,&result);
  CHECK(status==0 && !result.accepted && result.exact_zero_proposal);
  entry(&f,0,1,0.5);
  status=mvmc_krylov_refine_final(&model,&policy,&vacuum,&pair,1,0.25,
      base,&f,NULL,&result);
  CHECK(status==0 && !result.accepted && !result.exact_zero_proposal);
  status=mvmc_krylov_refine_final(&model,&policy,&vacuum,&pair,1,
      nextafter(0.25,0),base,&f,NULL,&result);
  CHECK(status==0 && result.accepted);
  /* Complex phase, cancellation in frozen coefficients, exp overflow. */
  entry(&f,0,1,0.25+0.5*I);
  c[1]=-0.5+0x1p-54;
  CHECK(mvmc_krylov_final_state_policy_create(1,
      MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT,20,c,2,&policy)==0);
  f.exponent=1000;
  status=mvmc_krylov_refine_final(&model,&policy,&vacuum,&pair,1,0.9,
      base,&f,NULL,&result);
  CHECK(status==0 && result.accepted);
  if(status==0) {
    CHECK(mvmc_krylov_final_state_local_energy_export(&result.evaluation,
                                                    &energy)==0);
    CHECK(cabs(energy-2)<1e-13);
    CHECK(result.evaluation.amplitude.log_abs>900);
  }
  f.exponent=-1000;
  status=mvmc_krylov_refine_final(&model,&policy,&vacuum,&pair,1,0.9,
      base,&f,NULL,&result);
  CHECK(status==0 && !result.accepted);
  /* Failure is transactional; resource limits never produce acceptance. */
  memset(&sentinel,0x5a,sizeof(sentinel)); result=sentinel;
  limits.max_node_calls=1;
  CHECK(mvmc_krylov_refine_final(&model,&policy,&vacuum,&pair,1,0.9,
      base,&f,&limits,&result)==MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
  CHECK(memcmp(&result,&sentinel,sizeof(result))==0);
  f.failure=1;
  CHECK(mvmc_krylov_refine_final(&model,&policy,&vacuum,&pair,1,0.9,
      base,&f,NULL,&result)==MVMC_KRYLOV_STATUS_NONFINITE);
  CHECK(memcmp(&result,&sentinel,sizeof(result))==0);
  f.failure=0; f.exponent=0;
  /* These enter the compiled C++ backend, not a mock finite guard.
   * GMP must never receive NaN/Inf, and failures must leave output intact. */
  {
    const double bad[] = {NAN, INFINITY, -INFINITY};
    size_t i;
    for (i=0; i<sizeof(bad)/sizeof(bad[0]); ++i) {
      result=sentinel;
      CHECK(mvmc_krylov_refine_final(&model,&policy,&vacuum,&pair,1,bad[i],
          base,&f,NULL,&result)==MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
      CHECK(memcmp(&result,&sentinel,sizeof(result))==0);
      f.weight=bad[i];
      CHECK(mvmc_krylov_refine_final(&model,&policy,&vacuum,&pair,1,0.5,
          base,&f,NULL,&result)==MVMC_KRYLOV_STATUS_NONFINITE);
      CHECK(memcmp(&result,&sentinel,sizeof(result))==0);
      f.weight=1; f.exponent=bad[i];
      CHECK(mvmc_krylov_refine_final(&model,&policy,&vacuum,&pair,1,0.5,
          base,&f,NULL,&result)==MVMC_KRYLOV_STATUS_NONFINITE);
      CHECK(memcmp(&result,&sentinel,sizeof(result))==0);
      f.exponent=0;
    }
  }
  /* Psi(pair)=0 but H Psi(pair)=1: no base support mask on H Psi/Phi. */
  {
    MVMCKrylovFermionOperator ops[4] = {
      {MVMC_KRYLOV_FERMION_CREATE,0},{MVMC_KRYLOV_FERMION_CREATE,1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE,1},{MVMC_KRYLOV_FERMION_ANNIHILATE,0}};
    MVMCKrylovHamiltonianTerm terms[2]={{1,0,2,0,0},{1,2,2,0,0}};
    memset(f.slater,0,sizeof(f.slater));
    model.terms=terms; model.term_count=2; model.operators=ops; model.operator_count=4;
    c[0]=1; c[1]=0;
    CHECK(mvmc_krylov_final_state_policy_create(1,
        MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT,23,c,2,&policy)==0);
    status=mvmc_krylov_refine_final(&model,&policy,&vacuum,&pair,1,0,
        base,&f,NULL,&result);
    CHECK(status==0 && !result.accepted && result.exact_zero_proposal);
    c[0]=1; c[1]=1;
    CHECK(mvmc_krylov_final_state_policy_create(1,
        MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT,21,c,2,&policy)==0);
    status=mvmc_krylov_refine_final(&model,&policy,&vacuum,&pair,1,0.75,
        base,&f,NULL,&result);
    CHECK(status==0 && result.accepted);
    if(status==0) {
      CHECK(result.krylov.value[0].state==MVMC_SCALED_COMPLEX_EXACT_ZERO);
      CHECK(result.krylov.value[1].state==MVMC_SCALED_COMPLEX_FINITE_NONZERO);
      CHECK(mvmc_krylov_final_state_local_energy_export(&result.evaluation,
                                                      &energy)==0);
      CHECK(cabs(energy-1)<1e-13);
    }
    /* 1-exp(2^-160) needs more than 128 bits.  Distinct exponents
     * cannot be cancelled by rational symbolic grouping. */
    entry(&f,0,1,-1);
    f.exponent=0x1p-160;
    status=mvmc_krylov_refine_final(&model,&policy,&vacuum,&pair,1,0.5,
        base,&f,NULL,&result);
    CHECK(status==0 && result.accepted && result.precision>=256);
    if(status==0) {
      CHECK(fabs(result.evaluation.amplitude.log_abs+160*log(2))<1e-12);
      CHECK(mvmc_krylov_final_state_local_energy_export(&result.evaluation,
                                                      &energy)==0);
      CHECK(cabs(energy-1)<1e-13);
    }
    limits.max_node_calls=1000; limits.max_precision=128;
    result=sentinel;
    CHECK(mvmc_krylov_refine_final(&model,&policy,&vacuum,&pair,1,0.5,
        base,&f,&limits,&result)==MVMC_KRYLOV_STATUS_RESOURCE_LIMIT);
    CHECK(memcmp(&result,&sentinel,sizeof(result))==0);
    /* A zero Hamiltonian gives an exactly zero, exportable local energy. */
    f.exponent=0;
    model.terms=NULL; model.term_count=0;
    model.operators=NULL; model.operator_count=0;
    status=mvmc_krylov_refine_final(&model,&policy,&vacuum,&pair,1,0.5,
        base,&f,NULL,&result);
    CHECK(status==0 && result.accepted);
    if(status==0) {
      MVMCKrylovFinalStateSnapshot snapshot;
      result.krylov.statistics.amplitude_generation_hash=123;
      CHECK(result.evaluation.local_energy.state==MVMC_SCALED_COMPLEX_EXACT_ZERO);
      CHECK(mvmc_krylov_final_state_snapshot_from_evaluation(&policy,&pair,1,
          &result.krylov,&result.evaluation,&snapshot)==0);
      CHECK(snapshot.amplitude_generation_hash==123 && snapshot.integrity_hash!=0);
      CHECK(snapshot.integrity_hash==
          mvmc_krylov_final_state_snapshot_integrity_hash(&snapshot));
    }
  }
}

/* Replicated full QP data, including an MPI rank owning no QP component. */
static void collective_case(int rank,int size) {
  MVMCGCKrylovAmplitudeLayout layout = {0};
  MVMCGCKrylovAmplitudeWorkspace *workspace = NULL;
  MVMCKrylovFockModel model = {0};
  MVMCKrylovFinalStatePolicy policy;
  MVMCKrylovRefinementResult result;
  double complex slater[16]={0}, weight=1, c[2]={1,0};
  uint64_t vacuum=0,pair=3;
  layout.site_count=2; layout.domain=MVMC_KRYLOV_DOMAIN_GC_EVEN;
  layout.qp_total=1; layout.qp_start=rank/size; layout.qp_end=(rank+1)/size;
  layout.translation_count=1; layout.translation_weights=&weight;
  layout.projection.site_count=2;
#ifdef _mpi_use
  layout.communicator=MPI_COMM_WORLD;
#endif
  slater[1]=0.5; slater[4]=-0.5;
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(&layout,slater,&weight,
                                                 &workspace)==0);
  model.site_count=2; model.hermitian=1;
  model.domain=MVMC_KRYLOV_DOMAIN_GC_EVEN;
  model.terms=identity_terms; model.term_count=2;
  model.operators=identity_ops; model.operator_count=4;
  CHECK(mvmc_krylov_final_state_policy_create(1,
      MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT,22,c,2,&policy)==0);
  if(workspace) {
    CHECK(mvmc_gc_krylov_refine_final(workspace,&model,&policy,&vacuum,
                                     &pair,1,0.1,&result)==0);
    CHECK(result.accepted && result.valid);
    if(size>1) {
      CHECK(mvmc_gc_krylov_refine_final(workspace,&model,&policy,&vacuum,
          &pair,1,rank==1 ? 0.2 : 0.1,&result)==
          MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE);
      {
        MVMCKrylovHamiltonianTerm bad_terms[2];
        memcpy(bad_terms,identity_terms,sizeof(bad_terms));
        if(rank==1) bad_terms[0].coefficient=NAN;
        model.terms=bad_terms;
        CHECK(mvmc_gc_krylov_refine_final(workspace,&model,&policy,&vacuum,
            &pair,1,0.1,&result)==MVMC_KRYLOV_STATUS_NONFINITE);
      }
    }
    mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
  }
}
int main(int argc,char **argv) {
  int rank=0,size=1,total;
#ifdef _mpi_use
  MPI_Init(&argc,&argv); MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
#else
  (void)argc; (void)argv;
#endif
  analytic_cases(); collective_case(rank,size); total=failures;
#ifdef _mpi_use
  MPI_Allreduce(&failures,&total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  if(rank==0) printf("Krylov refinement: %s\n",total ? "FAIL" : "PASS");
  return total ? 1 : 0;
}
