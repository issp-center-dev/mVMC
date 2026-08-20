#include "krylov_final_state_chain.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failures = 0;
static int world_rank = 0;

#define CHECK(condition, ...)                                                 \
  do {                                                                        \
    if (!(condition)) {                                                       \
      fprintf(stderr, "KrylovFinalStateChain_Unit FAIL rank %d: ",          \
              world_rank);                                                    \
      fprintf(stderr, __VA_ARGS__);                                           \
      fprintf(stderr, " (line %d)\n", __LINE__);                            \
      ++failures;                                                             \
    }                                                                         \
  } while (0)

typedef struct {
  double complex values[64];
  MVMCKrylovStatus forced_status;
} TableAmplitude;

typedef struct {
  uint64_t words;
  MVMCKrylovFinalStateSnapshot current;
  MVMCKrylovPositiveSamplerRng rng;
  uint64_t trace_hash;
} ChainRun;

static int rng_equal(const MVMCKrylovPositiveSamplerRng *left,
                     const MVMCKrylovPositiveSamplerRng *right) {
  return left->valid == right->valid &&
         left->state_version == right->state_version &&
         left->algorithm == right->algorithm && left->state == right->state &&
         left->stream == right->stream && left->draws == right->draws;
}

static void hash_u64(uint64_t *hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    *hash ^= value & UINT64_C(0xff);
    *hash *= UINT64_C(1099511628211);
    value >>= 8;
  }
}

static MVMCKrylovStatus table_callback(
    const uint64_t *words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
  TableAmplitude *table = (TableAmplitude *)context;
  const size_t index = word_count == 1
                           ? (size_t)(words[0] & UINT64_C(63))
                           : 63;
  if (table == NULL || result == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (table->forced_status != MVMC_KRYLOV_STATUS_OK) {
    return table->forced_status;
  }
  memset(result, 0, sizeof(*result));
  if (creal(table->values[index]) == 0.0 &&
      cimag(table->values[index]) == 0.0) {
    if (mvmc_scaled_complex_make_exact_zero(&result->value) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    result->exact_zero_component_count = 1;
  } else {
    if (mvmc_scaled_complex_from_raw(
            table->values[index], &result->value) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    }
    result->well_pivoted_component_count = 1;
  }
  result->local_factorization_count = 1;
  result->global_factorization_count = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovFockModel two_state_model(void) {
  static const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
  };
  static const MVMCKrylovHamiltonianTerm terms[] = {
      {1.0, 0, 2, 1, 0},
      {1.0, 2, 2, 1, 1},
  };
  MVMCKrylovFockModel model;
  memset(&model, 0, sizeof(model));
  model.site_count = 2;
  model.up_electron_count = 1;
  model.hermitian = 1;
  model.terms = terms;
  model.term_count = 2;
  model.operators = operators;
  model.operator_count = 4;
  return model;
}

static MVMCKrylovBoundedLimits bounded_limits(void) {
  MVMCKrylovBoundedLimits limits;
  memset(&limits, 0, sizeof(limits));
  limits.amplitude_policy_hash = UINT64_C(0x503543434841494e);
  limits.cache_bytes = 4096;
  limits.max_row_transitions = 16;
  limits.max_workspace_bytes = (size_t)1024 * 1024;
  limits.max_node_expansions = 1000;
  limits.max_terminal_amplitude_calls = 1000;
  limits.max_total_row_transitions = 10000;
  limits.max_order = 2;
  return limits;
}

static MVMCKrylovBoundedWorkspace *create_workspace(
    const MVMCKrylovFockModel *model, const MVMCKrylovBoundedLimits *limits,
    MVMCKrylovBoundedPlan **plan_out) {
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  CHECK(mvmc_bounded_krylov_plan_create(model, limits, &plan) ==
            MVMC_KRYLOV_STATUS_OK &&
            plan != NULL,
        "chain plan creation");
  if (plan != NULL) {
    CHECK(mvmc_bounded_krylov_workspace_create(plan, &workspace) ==
              MVMC_KRYLOV_STATUS_OK &&
              workspace != NULL,
          "chain workspace creation");
  }
  *plan_out = plan;
  return workspace;
}

static void destroy_workspace(MVMCKrylovBoundedWorkspace *workspace,
                              MVMCKrylovBoundedPlan *plan) {
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
}

static MVMCKrylovFinalStatePolicy final_policy(double second_coefficient) {
  const double complex coefficient[] = {1.0, second_coefficient};
  MVMCKrylovFinalStatePolicy policy;
  CHECK(mvmc_krylov_final_state_policy_create(
            1, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC,
            UINT64_C(0x503543434f454646), coefficient, 2, &policy) ==
            MVMC_KRYLOV_STATUS_OK,
        "chain final policy");
  return policy;
}

static MVMCKrylovPositiveSamplerProposalPolicy proposal_policy(
    size_t global_numerator, size_t global_denominator) {
  MVMCKrylovPositiveSamplerProposalPolicy policy;
  CHECK(mvmc_krylov_positive_sampler_proposal_policy_create(
            global_numerator, global_denominator, &policy) ==
            MVMC_KRYLOV_STATUS_OK,
        "chain proposal policy");
  return policy;
}

static void initialize_run(MVMCKrylovBoundedWorkspace *workspace,
                           const MVMCKrylovFinalStatePolicy *policy,
                           TableAmplitude *amplitude, ChainRun *run) {
  memset(run, 0, sizeof(*run));
  run->words = UINT64_C(1);
  run->trace_hash = UINT64_C(1469598103934665603);
  CHECK(mvmc_bounded_krylov_session_begin(
            workspace, table_callback, amplitude,
            UINT64_C(0x503543414d504745)) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_final_state_sampler_initialize(
                workspace, policy, &run->words, 1, table_callback,
                amplitude, &run->current) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_positive_sampler_rng_seed(
                UINT64_C(0x503543554e495431),
                UINT64_C(0x5035435354524541), &run->rng) ==
                MVMC_KRYLOV_STATUS_OK,
        "chain run initialization");
  CHECK(run->current.integrity_hash != 0 &&
            run->current.integrity_hash ==
                mvmc_krylov_final_state_snapshot_integrity_hash(
                    &run->current),
        "chain snapshot integrity initialized");
}

static void run_steps(MVMCKrylovBoundedWorkspace *workspace,
                      const MVMCKrylovFinalStatePolicy *policy,
                      const MVMCKrylovFockModel *model,
                      const MVMCKrylovPositiveSamplerProposalPolicy *proposal,
                      TableAmplitude *amplitude, size_t count, ChainRun *run) {
  size_t step_index;
  for (step_index = 0; step_index < count; ++step_index) {
    const uint64_t draws_before = run->rng.draws;
    uint64_t proposal_words = 0;
    MVMCKrylovFinalStateChainStepResult step;
    const MVMCKrylovStatus status =
        mvmc_krylov_final_state_chain_step_mixture_rng(
            workspace, policy, model, proposal, &run->words, 1,
            &run->current, &run->rng, table_callback, amplitude,
            &proposal_words, 1, &step);
    CHECK(status == MVMC_KRYLOV_STATUS_OK && step.valid &&
              step.rng_after.state == run->rng.state &&
              step.rng_after.draws == run->rng.draws &&
              run->rng.draws > draws_before &&
              step.proposal.proposal_policy_hash == proposal->policy_hash,
          "chain step %zu", step_index);
    if (status != MVMC_KRYLOV_STATUS_OK || !step.valid) return;
    hash_u64(&run->trace_hash, run->words);
    hash_u64(&run->trace_hash, (uint64_t)(unsigned int)step.step.accepted);
    hash_u64(&run->trace_hash,
             (uint64_t)(unsigned int)step.configuration_changed);
    hash_u64(&run->trace_hash, run->rng.draws);
  }
}

static void test_observable_inventory(void) {
  const double complex values[] = {1.0, 2.0, 3.0};
  MVMCKrylovBoundedResult krylov;
  MVMCKrylovFinalStateEvaluation evaluation;
  const MVMCKrylovFinalStatePolicy policy = final_policy(2.0);
  double complex measured = NAN + I * NAN;
  int order;
  memset(&krylov, 0, sizeof(krylov));
  krylov.valid = 1;
  krylov.status = MVMC_KRYLOV_STATUS_OK;
  krylov.evaluated_order = 2;
  for (order = 0; order <= 2; ++order) {
    CHECK(mvmc_scaled_complex_from_raw(
              values[order], &krylov.value[order]) ==
              MVMC_PFAFFIAN_STATUS_OK,
          "observable raw value %d", order);
  }
  CHECK(mvmc_krylov_final_state_evaluate(
            &policy, &krylov, &evaluation) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_final_state_observable_measure(
                MVMC_KRYLOV_FINAL_OBSERVABLE_ENERGY_COMPLEX,
                &evaluation, &measured) == MVMC_KRYLOV_STATUS_OK &&
            cabs(measured - 8.0 / 5.0) <= 2.0e-14,
        "energy observable measurement");
  CHECK(mvmc_krylov_final_state_observable_measure(
            MVMC_KRYLOV_FINAL_OBSERVABLE_LOCAL_ENERGY_ABS_SQUARED_DIAGNOSTIC,
            &evaluation, &measured) == MVMC_KRYLOV_STATUS_OK &&
            fabs(creal(measured) - 64.0 / 25.0) <= 2.0e-14 &&
            cimag(measured) == 0.0,
        "direct second-moment diagnostic measurement");
  measured = 0.0;
  CHECK(mvmc_krylov_final_state_observable_measure(
            MVMC_KRYLOV_FINAL_OBSERVABLE_ARBITRARY_UNSUPPORTED,
            &evaluation, &measured) ==
            MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL &&
            isnan(creal(measured)) && isnan(cimag(measured)),
        "unsupported observable fails fast");
}

static void test_transaction_restart_and_manifest(void) {
  const MVMCKrylovFockModel model = two_state_model();
  const MVMCKrylovBoundedLimits limits = bounded_limits();
  const MVMCKrylovFinalStatePolicy policy = final_policy(2.0);
  const MVMCKrylovFinalStatePolicy different_policy = final_policy(2.001);
  const MVMCKrylovPositiveSamplerProposalPolicy proposal =
      proposal_policy(1, 1);
  const MVMCKrylovPositiveSamplerProposalPolicy different_proposal =
      proposal_policy(1, 2);
  MVMCKrylovBoundedPlan *full_plan = NULL;
  MVMCKrylovBoundedPlan *split_plan = NULL;
  MVMCKrylovBoundedWorkspace *full_workspace =
      create_workspace(&model, &limits, &full_plan);
  MVMCKrylovBoundedWorkspace *split_workspace =
      create_workspace(&model, &limits, &split_plan);
  TableAmplitude amplitude;
  ChainRun full;
  ChainRun split;
  MVMCKrylovFinalStateChainManifest manifest;
  MVMCKrylovFinalStateChainManifest final_manifest;
  uint64_t restart_words = 0;
  int matches = 0;

  memset(&amplitude, 0, sizeof(amplitude));
  amplitude.values[1] = 1.0;
  initialize_run(full_workspace, &policy, &amplitude, &full);
  initialize_run(split_workspace, &policy, &amplitude, &split);
  run_steps(full_workspace, &policy, &model, &proposal, &amplitude, 20,
            &full);
  run_steps(split_workspace, &policy, &model, &proposal, &amplitude, 7,
            &split);
  CHECK(mvmc_krylov_final_state_chain_manifest_create(
            &policy, &proposal, &model, &limits,
            mvmc_bounded_krylov_plan_hash(split_plan), &split.words, 1,
            &split.current, &split.rng, &manifest) ==
            MVMC_KRYLOV_STATUS_OK &&
            manifest.valid && manifest.component_state_hash != 0 &&
            manifest.final_policy_hash == policy.policy_hash &&
            manifest.coefficient_provenance_hash == policy.provenance_hash,
        "split restart manifest creation");
  restart_words = split.words;
  CHECK(mvmc_krylov_final_state_chain_manifest_matches(
            &manifest, &policy, &proposal, &model, &limits,
            mvmc_bounded_krylov_plan_hash(split_plan), &split.words, 1,
            &split.current, &split.rng, &matches) ==
            MVMC_KRYLOV_STATUS_OK &&
            matches,
        "live manifest match");
  {
    MVMCKrylovFinalStateSnapshot corrupt = split.current;
    MVMCKrylovFinalStateChainManifest rejected;
    corrupt.final_state.log_weight += 1.0e-6;
    CHECK(mvmc_krylov_final_state_chain_manifest_create(
              &policy, &proposal, &model, &limits,
              mvmc_bounded_krylov_plan_hash(split_plan), &split.words, 1,
              &corrupt, &split.rng, &rejected) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              !rejected.valid,
          "corrupt final evaluation rejected by manifest");
  }
  {
    MVMCKrylovFinalStateSnapshot corrupt = split.current;
    MVMCKrylovFinalStateChainManifest rejected;
    corrupt.krylov.value[1].log_abs += 1.0e-6;
    CHECK(mvmc_krylov_final_state_evaluate(
              &policy, &corrupt.krylov, &corrupt.final_state) ==
              MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_final_state_chain_manifest_create(
                  &policy, &proposal, &model, &limits,
                  mvmc_bounded_krylov_plan_hash(split_plan), &split.words, 1,
                  &corrupt, &split.rng, &rejected) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              !rejected.valid,
          "self-consistent corrupt Krylov snapshot rejected by integrity");
  }
  matches = 1;
  CHECK(mvmc_krylov_final_state_chain_manifest_matches(
            &manifest, &policy, &different_proposal, &model, &limits,
            mvmc_bounded_krylov_plan_hash(split_plan), &split.words, 1,
            &split.current, &split.rng, &matches) ==
            MVMC_KRYLOV_STATUS_OK &&
            !matches,
        "proposal-policy mismatch rejected");
  matches = 1;
  {
    MVMCKrylovPositiveSamplerRng wrong_rng = split.rng;
    ++wrong_rng.draws;
    CHECK(mvmc_krylov_final_state_chain_manifest_matches(
              &manifest, &policy, &proposal, &model, &limits,
              mvmc_bounded_krylov_plan_hash(split_plan), &split.words, 1,
              &split.current, &wrong_rng, &matches) ==
              MVMC_KRYLOV_STATUS_OK &&
              !matches,
          "RNG draw-count mismatch rejected");
  }

  destroy_workspace(split_workspace, split_plan);
  split_workspace = create_workspace(&model, &limits, &split_plan);
  memset(&split.current, 0, sizeof(split.current));
  memset(&split.rng, 0, sizeof(split.rng));
  CHECK(mvmc_krylov_final_state_chain_restart_restore(
            split_workspace, &manifest, &policy, &proposal, &model,
            &limits, mvmc_bounded_krylov_plan_hash(split_plan),
            &split.words, 1, table_callback, &amplitude, &split.current,
            &split.rng) == MVMC_KRYLOV_STATUS_OK,
        "cold restart restore");
  run_steps(split_workspace, &policy, &model, &proposal, &amplitude, 13,
            &split);
  CHECK(full.words == split.words &&
            full.trace_hash == split.trace_hash &&
            rng_equal(&full.rng, &split.rng) &&
            full.current.accepted_generation ==
                split.current.accepted_generation,
        "full/split trace identity");
  CHECK(mvmc_krylov_final_state_chain_manifest_create(
            &policy, &proposal, &model, &limits,
            mvmc_bounded_krylov_plan_hash(split_plan), &split.words, 1,
            &split.current, &split.rng, &final_manifest) ==
            MVMC_KRYLOV_STATUS_OK &&
            final_manifest.rng_draws == split.rng.draws,
        "final manifest creation");

  {
    MVMCKrylovFinalStateSnapshot rejected_current;
    MVMCKrylovPositiveSamplerRng rejected_rng;
    memset(&rejected_current, 0, sizeof(rejected_current));
    memset(&rejected_rng, 0, sizeof(rejected_rng));
    CHECK(mvmc_krylov_final_state_chain_restart_restore(
              split_workspace, &manifest, &different_policy, &proposal,
              &model, &limits,
              mvmc_bounded_krylov_plan_hash(split_plan), &restart_words,
              1, table_callback, &amplitude, &rejected_current,
              &rejected_rng) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              !rejected_current.valid && !rejected_rng.valid,
          "coefficient-policy restart mismatch rejected");
  }
  {
    MVMCKrylovFinalStateChainManifest zero_generation = manifest;
    MVMCKrylovFinalStateSnapshot rejected_current;
    MVMCKrylovPositiveSamplerRng rejected_rng;
    zero_generation.amplitude_generation_hash = 0;
    memset(&rejected_current, 0, sizeof(rejected_current));
    memset(&rejected_rng, 0, sizeof(rejected_rng));
    CHECK(mvmc_krylov_final_state_chain_restart_restore(
              split_workspace, &zero_generation, &policy, &proposal, &model,
              &limits, mvmc_bounded_krylov_plan_hash(split_plan),
              &restart_words, 1, table_callback, &amplitude,
              &rejected_current, &rejected_rng) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              !rejected_current.valid && !rejected_rng.valid,
          "zero-generation restart manifest rejected");
  }

  {
    uint64_t proposal_words = 0;
    const uint64_t saved_words = split.words;
    const MVMCKrylovPositiveSamplerRng saved_rng = split.rng;
    const MVMCKrylovFinalStateSnapshot saved_current = split.current;
    MVMCKrylovFinalStatePolicy invalid_policy = policy;
    MVMCKrylovFinalStateChainStepResult step;
    invalid_policy.policy_hash ^= UINT64_C(1);
    CHECK(mvmc_krylov_final_state_chain_step_mixture_rng(
              split_workspace, &invalid_policy, &model, &proposal,
              &split.words, 1, &split.current, &split.rng, table_callback,
              &amplitude, &proposal_words, 1, &step) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              !step.valid && split.words == saved_words &&
              rng_equal(&split.rng, &saved_rng) &&
              memcmp(&split.current, &saved_current,
                     sizeof(saved_current)) == 0,
          "failed step leaves chain and RNG unchanged");
  }

  {
    uint64_t proposal_words = 0;
    const uint64_t saved_words = split.words;
    const MVMCKrylovPositiveSamplerRng saved_rng = split.rng;
    const MVMCKrylovFinalStateSnapshot saved_current = split.current;
    MVMCKrylovFinalStateChainStepResult step;
    CHECK(mvmc_bounded_krylov_session_end(split_workspace) ==
              MVMC_KRYLOV_STATUS_OK &&
              mvmc_bounded_krylov_session_begin(
                  split_workspace, table_callback, &amplitude,
                  split.current.amplitude_generation_hash) ==
              MVMC_KRYLOV_STATUS_OK,
          "empty-cache session restart before callback failure");
    amplitude.forced_status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    CHECK(mvmc_krylov_final_state_chain_step_mixture_rng(
              split_workspace, &policy, &model, &proposal, &split.words, 1,
              &split.current, &split.rng, table_callback, &amplitude,
              &proposal_words, 1, &step) ==
              MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
              !step.valid && split.words == saved_words &&
              rng_equal(&split.rng, &saved_rng) &&
              memcmp(&split.current, &saved_current,
                     sizeof(saved_current)) == 0 &&
              mvmc_bounded_krylov_session_is_active(split_workspace),
          "callback failure is transactional and session is recovered");
    amplitude.forced_status = MVMC_KRYLOV_STATUS_OK;
    CHECK(mvmc_krylov_final_state_chain_step_mixture_rng(
              split_workspace, &policy, &model, &proposal, &split.words, 1,
              &split.current, &split.rng, table_callback, &amplitude,
              &proposal_words, 1, &step) == MVMC_KRYLOV_STATUS_OK &&
              step.valid && split.rng.draws > saved_rng.draws,
          "chain continues after transient callback failure");
  }

#ifdef _mpi_use
  {
    uint64_t local[] = {final_manifest.component_state_hash,
                        final_manifest.proposal_model_hash,
                        final_manifest.rng_state,
                        final_manifest.rng_draws,
                        split.trace_hash};
    uint64_t minimum[5];
    uint64_t maximum[5];
    MPI_Allreduce(local, minimum, 5, MPI_UINT64_T, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(local, maximum, 5, MPI_UINT64_T, MPI_MAX,
                  MPI_COMM_WORLD);
    CHECK(memcmp(minimum, maximum, sizeof(minimum)) == 0,
          "chain manifest/trace differs by MPI rank");
  }
#endif

  destroy_workspace(full_workspace, full_plan);
  destroy_workspace(split_workspace, split_plan);
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#else
  (void)argc;
  (void)argv;
#endif
  test_observable_inventory();
  test_transaction_restart_and_manifest();
#ifdef _mpi_use
  {
    int global_failures = 0;
    MPI_Allreduce(&failures, &global_failures, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    failures = global_failures;
  }
#endif
  if (failures == 0 && world_rank == 0) {
    puts("krylov final-state chain unit: PASS");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return failures == 0 ? 0 : 1;
}
