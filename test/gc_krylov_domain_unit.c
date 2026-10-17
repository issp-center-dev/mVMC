/*
 * Configuration-domain contract for the grand-canonical first-step
 * stabilized Lanczos: fixed-N, GC-balanced (Nup=Ndown, all even N) and
 * GC-even (all even N) validation in the reference and bounded engines,
 * plan-hash identity, fixed-sector proposal rejection and the classic model
 * domain creator (rank-0 metadata, collective mismatch rejection).
 */
#include "bounded_krylov_engine.h"
#include "classic_krylov_model.h"
#include "krylov_fock_proposal.h"
#include "krylov_positive_sampler.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failure_count = 0;

#define CHECK(condition, message)                                    \
  do {                                                               \
    if (!(condition)) {                                              \
      fprintf(stderr, "FAIL: %s (line %d)\n", (message), __LINE__); \
      ++failure_count;                                               \
    }                                                                \
  } while (0)

/* Base-commit 2f394cac plan hash of fixed_hash_model(): fixed domain keeps it. */
#define FIXED_PLAN_HASH UINT64_C(0xadb00ea4e6234af6)

static int world_rank(void) {
  int rank = 0;
#ifdef _mpi_use
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  return rank;
}

static int world_size(void) {
  int size = 1;
#ifdef _mpi_use
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  return size;
}

static MVMCClassicPfaffianCommunicator world_communicator(void) {
#ifdef _mpi_use
  return MPI_COMM_WORLD;
#else
  return 0;
#endif
}

static int popcount_mask(unsigned mask) {
  int count = 0;
  while (mask != 0) {
    count += (int)(mask & 1u);
    mask >>= 1;
  }
  return count;
}

/* L=2: hopping -1 on both spins, U=3/2, chemical potential -1/4 n. */
static const MVMCKrylovHamiltonianTerm kTerms[10] = {
    {-1.0 + 0.0 * I, 0, 2, 1, 0},  {-1.0 + 0.0 * I, 2, 2, 1, 1},
    {-1.0 + 0.0 * I, 4, 2, 1, 2},  {-1.0 + 0.0 * I, 6, 2, 1, 3},
    {1.5 + 0.0 * I, 8, 4, 2, 0},   {1.5 + 0.0 * I, 12, 4, 2, 1},
    {-0.25 + 0.0 * I, 16, 2, 1, 4}, {-0.25 + 0.0 * I, 18, 2, 1, 5},
    {-0.25 + 0.0 * I, 20, 2, 1, 6}, {-0.25 + 0.0 * I, 22, 2, 1, 7}};
static const MVMCKrylovFermionOperator kOperators[24] = {
    {MVMC_KRYLOV_FERMION_CREATE, 0}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
    {MVMC_KRYLOV_FERMION_CREATE, 1}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
    {MVMC_KRYLOV_FERMION_CREATE, 2}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 3},
    {MVMC_KRYLOV_FERMION_CREATE, 3}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
    {MVMC_KRYLOV_FERMION_CREATE, 0}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
    {MVMC_KRYLOV_FERMION_CREATE, 2}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
    {MVMC_KRYLOV_FERMION_CREATE, 1}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
    {MVMC_KRYLOV_FERMION_CREATE, 3}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 3},
    {MVMC_KRYLOV_FERMION_CREATE, 0}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
    {MVMC_KRYLOV_FERMION_CREATE, 1}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
    {MVMC_KRYLOV_FERMION_CREATE, 2}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
    {MVMC_KRYLOV_FERMION_CREATE, 3}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 3}};

static MVMCKrylovFockModel two_site_model(MVMCKrylovConfigurationDomain domain,
                                          size_t up, size_t down) {
  MVMCKrylovFockModel model;
  memset(&model, 0, sizeof(model));
  model.site_count = 2;
  model.up_electron_count = up;
  model.down_electron_count = down;
  model.hermitian = 1;
  model.terms = kTerms;
  model.term_count = 10;
  model.operators = kOperators;
  model.operator_count = 24;
  model.domain = domain;
  return model;
}

/* Same model as the base-commit hash probe (tmp/.../base_plan_hash_probe.c). */
static const MVMCKrylovHamiltonianTerm kHashTerms[3] = {
    {-1.0 + 0.0 * I, 0, 2, 1, 0},
    {-1.0 + 0.0 * I, 2, 2, 1, 1},
    {1.5 + 0.0 * I, 4, 4, 2, 0}};
static const MVMCKrylovFermionOperator kHashOperators[8] = {
    {MVMC_KRYLOV_FERMION_CREATE, 0}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
    {MVMC_KRYLOV_FERMION_CREATE, 1}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
    {MVMC_KRYLOV_FERMION_CREATE, 0}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
    {MVMC_KRYLOV_FERMION_CREATE, 2}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 2}};

static MVMCKrylovFockModel fixed_hash_model(void) {
  MVMCKrylovFockModel model;
  memset(&model, 0, sizeof(model));
  model.site_count = 2;
  model.up_electron_count = 1;
  model.down_electron_count = 1;
  model.hermitian = 1;
  model.terms = kHashTerms;
  model.term_count = 3;
  model.operators = kHashOperators;
  model.operator_count = 8;
  model.domain = MVMC_KRYLOV_DOMAIN_FIXED;
  return model;
}

static MVMCKrylovBoundedLimits hash_limits(void) {
  MVMCKrylovBoundedLimits limits;
  memset(&limits, 0, sizeof(limits));
  limits.amplitude_policy_hash = UINT64_C(0x1234);
  limits.cache_bytes = 4096;
  limits.max_row_transitions = 8;
  limits.max_workspace_bytes = 1u << 20;
  limits.max_node_expansions = 64;
  limits.max_terminal_amplitude_calls = 64;
  limits.max_total_row_transitions = 256;
  limits.max_order = 2;
  return limits;
}

static MVMCKrylovBoundedLimits engine_limits(void) {
  MVMCKrylovBoundedLimits limits;
  memset(&limits, 0, sizeof(limits));
  limits.amplitude_policy_hash = UINT64_C(0x474344304d41494e);
  limits.cache_bytes = 1u << 16;
  limits.max_row_transitions = 64;
  limits.max_workspace_bytes = (size_t)64 * 1024 * 1024;
  limits.max_node_expansions = UINT64_C(100000);
  limits.max_terminal_amplitude_calls = UINT64_C(100000);
  limits.max_total_row_transitions = UINT64_C(1000000);
  limits.max_order = 2;
  return limits;
}

static double complex table_value(uint64_t mask) {
  return (1.0 + (double)mask) + 0.25 * (double)mask * I;
}

static MVMCKrylovStatus scaled_callback(const uint64_t *words,
                                        size_t word_count, void *context,
                                        MVMCKrylovScaledAmplitudeResult *result) {
  (void)context;
  if (words == NULL || word_count != 1 || result == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(result, 0, sizeof(*result));
  if (mvmc_scaled_complex_from_raw_testing(table_value(words[0]),
                                            &result->value) !=
      MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  result->well_pivoted_component_count = 1;
  result->local_factorization_count = 1;
  result->global_factorization_count = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus raw_callback(const uint64_t *words, size_t word_count,
                                     void *context,
                                     MVMCKrylovAmplitudeResult *result) {
  (void)context;
  if (words == NULL || word_count != 1 || result == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(result, 0, sizeof(*result));
  result->value = table_value(words[0]);
  result->regular_component_count = 1;
  result->local_factorization_count = 1;
  result->global_factorization_count = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus bounded_values(MVMCKrylovBoundedWorkspace *workspace,
                                       uint64_t mask, int session,
                                       double complex values[3]) {
  MVMCKrylovBoundedResult result;
  MVMCKrylovStatus status;
  int order;
  memset(&result, 0, sizeof(result));
  status = session
               ? mvmc_bounded_krylov_session_evaluate(workspace, &mask, 1,
                                                      &result)
               : mvmc_bounded_krylov_evaluate(workspace, &mask, 1,
                                              scaled_callback, NULL, &result);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  for (order = 0; order <= 2; ++order) {
    MVMCScaledComplexExportStatus exported =
        mvmc_scaled_complex_export_common_scale(&result.value[order], 0.0,
                                                &values[order]);
    if (exported == MVMC_SCALED_EXPORT_EXACT_ZERO) {
      values[order] = 0.0;
    } else if (exported != MVMC_SCALED_EXPORT_OK) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static int expected_accept(MVMCKrylovConfigurationDomain domain,
                           unsigned mask) {
  const int up = popcount_mask(mask & 3u);
  const int down = popcount_mask((mask >> 2) & 3u);
  switch (domain) {
    case MVMC_KRYLOV_DOMAIN_FIXED:
      return up == 1 && down == 1;
    case MVMC_KRYLOV_DOMAIN_GC_BALANCED:
      return up == down;
    case MVMC_KRYLOV_DOMAIN_GC_EVEN:
      return (up + down) % 2 == 0;
  }
  return 0;
}

static void test_plan_example(void) {
  MVMCKrylovFockModel model;
  uint64_t words = UINT64_C(0x5); /* up0 + down0 */
  memset(&model, 0, sizeof(model));
  model.site_count = 2;
  model.hermitian = 1;
  model.domain = MVMC_KRYLOV_DOMAIN_GC_BALANCED;
  CHECK(mvmc_krylov_fock_validate(&model, &words, 1) == MVMC_KRYLOV_STATUS_OK,
        "balanced accepts up0+down0");
  words = UINT64_C(0x3); /* two up */
  CHECK(mvmc_krylov_fock_validate(&model, &words, 1) != MVMC_KRYLOV_STATUS_OK,
        "balanced rejects two up");
  model.domain = MVMC_KRYLOV_DOMAIN_GC_EVEN;
  CHECK(mvmc_krylov_fock_validate(&model, &words, 1) == MVMC_KRYLOV_STATUS_OK,
        "even accepts two up");
}

static void test_mask_sets(void) {
  const MVMCKrylovConfigurationDomain domains[3] = {
      MVMC_KRYLOV_DOMAIN_FIXED, MVMC_KRYLOV_DOMAIN_GC_BALANCED,
      MVMC_KRYLOV_DOMAIN_GC_EVEN};
  const int expected_counts[3] = {4, 6, 8};
  int d;
  for (d = 0; d < 3; ++d) {
    const int fixed = domains[d] == MVMC_KRYLOV_DOMAIN_FIXED;
    MVMCKrylovFockModel model =
        two_site_model(domains[d], fixed ? 1u : 0u, fixed ? 1u : 0u);
    MVMCKrylovBoundedLimits limits = engine_limits();
    MVMCKrylovBoundedPlan *plan = NULL;
    MVMCKrylovBoundedWorkspace *workspace = NULL;
    unsigned mask;
    int reference_count = 0;
    int bounded_count = 0;
    limits.max_order = 0;
    CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
              MVMC_KRYLOV_STATUS_OK,
          "plan for every domain");
    CHECK(plan != NULL && mvmc_bounded_krylov_workspace_create(
                              plan, &workspace) == MVMC_KRYLOV_STATUS_OK,
          "workspace for every domain");
    for (mask = 0; mask < 16u; ++mask) {
      uint64_t words = mask;
      MVMCKrylovBoundedResult result;
      const int expected = expected_accept(domains[d], mask);
      const MVMCKrylovStatus reference =
          mvmc_krylov_fock_validate(&model, &words, 1);
      MVMCKrylovStatus bounded = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      if (workspace != NULL) {
        bounded = mvmc_bounded_krylov_evaluate(workspace, &words, 1,
                                               scaled_callback, NULL, &result);
      }
      CHECK((reference == MVMC_KRYLOV_STATUS_OK) == expected,
            "reference domain acceptance");
      CHECK(expected || reference == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
            "reference rejection is INVALID_ARGUMENT");
      CHECK((bounded == MVMC_KRYLOV_STATUS_OK) == expected,
            "bounded domain acceptance");
      CHECK(expected || bounded == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
            "bounded rejection is INVALID_ARGUMENT");
      reference_count += reference == MVMC_KRYLOV_STATUS_OK;
      bounded_count += bounded == MVMC_KRYLOV_STATUS_OK;
    }
    CHECK(reference_count == expected_counts[d], "reference accepted count");
    CHECK(bounded_count == expected_counts[d], "bounded accepted count");
    {
      uint64_t excess = UINT64_C(0x15); /* bit 4 is outside 2L=4 */
      CHECK(mvmc_krylov_fock_validate(&model, &excess, 1) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
            "excess bit rejected");
    }
    mvmc_bounded_krylov_workspace_destroy(workspace);
    mvmc_bounded_krylov_plan_destroy(plan);
  }
}

static void test_model_validation(void) {
  MVMCKrylovBoundedLimits limits = engine_limits();
  MVMCKrylovBoundedPlan *plan = NULL;
  uint64_t vacuum = 0;
  MVMCKrylovFockModel model;

  model = two_site_model(MVMC_KRYLOV_DOMAIN_GC_BALANCED, 1, 1);
  CHECK(mvmc_krylov_fock_validate(&model, &vacuum, 1) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "GC model with fixed counts rejected by reference");
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "GC model with fixed counts rejected by bounded");
  CHECK(plan == NULL, "rejected plan stays NULL");

  model = two_site_model(MVMC_KRYLOV_DOMAIN_GC_EVEN, 0, 0);
  model.pure_spin = 1;
  CHECK(mvmc_krylov_fock_validate(&model, &vacuum, 1) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "GC pure_spin rejected by reference");
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "GC pure_spin rejected by bounded");

  model = two_site_model((MVMCKrylovConfigurationDomain)7, 0, 0);
  CHECK(mvmc_krylov_fock_validate(&model, &vacuum, 1) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "invalid domain rejected by reference");
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "invalid domain rejected by bounded");

  /* Fixed pure-spin regression: one electron per site. */
  model = two_site_model(MVMC_KRYLOV_DOMAIN_FIXED, 1, 1);
  model.pure_spin = 1;
  {
    uint64_t one_each = UINT64_C(0x9); /* up0 + down1 */
    uint64_t doublon = UINT64_C(0x5);  /* up0 + down0 */
    CHECK(mvmc_krylov_fock_validate(&model, &one_each, 1) ==
              MVMC_KRYLOV_STATUS_OK,
          "fixed pure-spin accepts singly occupied sites");
    CHECK(mvmc_krylov_fock_validate(&model, &doublon, 1) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "fixed pure-spin rejects doublon");
  }
  model = two_site_model(MVMC_KRYLOV_DOMAIN_FIXED, 0, 0);
  CHECK(mvmc_krylov_fock_validate(&model, &vacuum, 1) == MVMC_KRYLOV_STATUS_OK,
        "fixed vacuum sector keeps its previous contract");
}

static void test_multiword(void) {
  MVMCKrylovFockModel model;
  uint64_t words[2];
  memset(&model, 0, sizeof(model));
  model.site_count = 33;
  model.hermitian = 1;
  model.domain = MVMC_KRYLOV_DOMAIN_GC_BALANCED;
  words[0] = UINT64_C(1) << 32; /* up site 32 */
  words[1] = UINT64_C(1) << 1;  /* orbital 65 = down site 32 */
  CHECK(mvmc_krylov_fock_validate(&model, words, 2) == MVMC_KRYLOV_STATUS_OK,
        "multiword balanced pair across words");
  words[1] |= UINT64_C(1) << 2; /* orbital 66 does not exist */
  CHECK(mvmc_krylov_fock_validate(&model, words, 2) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "multiword excess bit");
  words[0] = 0;
  words[1] = UINT64_C(0x3); /* down sites 31 and 32 */
  CHECK(mvmc_krylov_fock_validate(&model, words, 2) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "multiword balanced rejects two down");
  model.domain = MVMC_KRYLOV_DOMAIN_GC_EVEN;
  CHECK(mvmc_krylov_fock_validate(&model, words, 2) == MVMC_KRYLOV_STATUS_OK,
        "multiword even accepts two down");
  words[1] = UINT64_C(0x1);
  CHECK(mvmc_krylov_fock_validate(&model, words, 2) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "multiword even rejects odd N");
  CHECK(mvmc_krylov_fock_validate(&model, words, 1) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "multiword word count");
}

static MVMCKrylovStatus zero_draw(void *context, size_t bound,
                                  size_t *value) {
  (void)context;
  if (bound == 0 || value == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *value = 0;
  return MVMC_KRYLOV_STATUS_OK;
}

static void test_fixed_sector_proposal_rejects_gc(void) {
  const MVMCKrylovConfigurationDomain domains[2] = {
      MVMC_KRYLOV_DOMAIN_GC_BALANCED, MVMC_KRYLOV_DOMAIN_GC_EVEN};
  int d;
  for (d = 0; d < 2; ++d) {
    const MVMCKrylovFockModel model = two_site_model(domains[d], 0, 0);
    uint64_t words = UINT64_C(0x5);
    uint64_t proposal = UINT64_C(0xdead);
    size_t count = 99;
    uint64_t hash = 1;
    MVMCKrylovFockUniformProposalResult uniform;
    CHECK(mvmc_krylov_fock_proposal_count_neighbors(&model, &words, 1,
                                                    &count) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "fixed-sector neighbor count rejects GC");
    CHECK(mvmc_krylov_fock_proposal_draw_uniform_sector(
              &model, zero_draw, NULL, &proposal, 1, &uniform) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "fixed-sector uniform draw rejects GC");
    CHECK(mvmc_krylov_positive_sampler_proposal_model_hash(&model, &words, 1,
                                                           &hash) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "proposal model hash rejects GC");
    CHECK(hash == 0, "rejected proposal hash is cleared");
  }
  {
    const MVMCKrylovFockModel fixed =
        two_site_model(MVMC_KRYLOV_DOMAIN_FIXED, 1, 1);
    uint64_t words = UINT64_C(0x5);
    size_t count = 0;
    uint64_t hash = 0;
    CHECK(mvmc_krylov_fock_proposal_count_neighbors(&fixed, &words, 1,
                                                    &count) ==
              MVMC_KRYLOV_STATUS_OK && count > 0,
          "fixed-sector proposal still accepts fixed domain");
    CHECK(mvmc_krylov_positive_sampler_proposal_model_hash(&fixed, &words, 1,
                                                           &hash) ==
              MVMC_KRYLOV_STATUS_OK && hash != 0,
          "proposal model hash for fixed domain");
  }
}

static uint64_t plan_hash_for(const MVMCKrylovFockModel *model) {
  MVMCKrylovBoundedLimits limits = hash_limits();
  MVMCKrylovBoundedPlan *plan = NULL;
  uint64_t hash = 0;
  if (mvmc_bounded_krylov_plan_create(model, &limits, &plan) ==
      MVMC_KRYLOV_STATUS_OK) {
    hash = mvmc_bounded_krylov_plan_hash(plan);
  }
  mvmc_bounded_krylov_plan_destroy(plan);
  return hash;
}

static void test_plan_hash(void) {
  MVMCKrylovFockModel model = fixed_hash_model();
  uint64_t fixed;
  uint64_t balanced;
  uint64_t even;
  fixed = plan_hash_for(&model);
  CHECK(fixed == FIXED_PLAN_HASH, "fixed-domain plan hash unchanged");
  model.up_electron_count = 0;
  model.down_electron_count = 0;
  model.domain = MVMC_KRYLOV_DOMAIN_GC_BALANCED;
  balanced = plan_hash_for(&model);
  model.domain = MVMC_KRYLOV_DOMAIN_GC_EVEN;
  even = plan_hash_for(&model);
  CHECK(balanced != 0 && even != 0, "GC plan hashes exist");
  CHECK(balanced != fixed && even != fixed && balanced != even,
        "domain is part of the plan identity");
  CHECK(plan_hash_for(&model) == even, "GC plan hash is deterministic");
}

static int close_complex(double complex actual, double complex expected) {
  return cabs(actual - expected) <= 1e-12 * (1.0 + cabs(expected));
}

static void test_bounded_values_across_sectors(void) {
  const MVMCKrylovFockModel model =
      two_site_model(MVMC_KRYLOV_DOMAIN_GC_EVEN, 0, 0);
  const uint64_t roots[6] = {UINT64_C(0xf), UINT64_C(0x5), UINT64_C(0x0),
                             UINT64_C(0x3), UINT64_C(0x5), UINT64_C(0x9)};
  MVMCKrylovBoundedLimits limits = engine_limits();
  MVMCKrylovLimits reference_limits;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  MVMCKrylovWorkspace *reference = NULL;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  size_t root;

  memset(&reference_limits, 0, sizeof(reference_limits));
  reference_limits.max_states = 4096;
  reference_limits.max_transitions = 65536;
  reference_limits.max_amplitude_evaluations = 4096;
  reference_limits.max_bytes = (size_t)64 * 1024 * 1024;
  reference_limits.max_order = 2;
  reference = mvmc_krylov_workspace_create(2, &reference_limits, &status);
  CHECK(reference != NULL && status == MVMC_KRYLOV_STATUS_OK,
        "reference workspace");
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_OK,
        "GC-even plan");
  CHECK(plan != NULL && mvmc_bounded_krylov_workspace_create(
                            plan, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "GC-even workspace");
  if (reference == NULL || workspace == NULL) {
    mvmc_bounded_krylov_workspace_destroy(workspace);
    mvmc_bounded_krylov_plan_destroy(plan);
    mvmc_krylov_workspace_destroy(reference);
    return;
  }
  CHECK(mvmc_bounded_krylov_session_begin(workspace, scaled_callback, NULL,
                                          UINT64_C(0x474344304b455931)) ==
            MVMC_KRYLOV_STATUS_OK,
        "session begin");
  for (root = 0; root < 6; ++root) {
    double complex session_values[3] = {0, 0, 0};
    double complex cold_values[3] = {0, 0, 0};
    MVMCKrylovResult expected;
    uint64_t words = roots[root];
    int order;
    CHECK(bounded_values(workspace, roots[root], 1, session_values) ==
              MVMC_KRYLOV_STATUS_OK,
          "session evaluation across N");
    CHECK(mvmc_krylov_evaluate(reference, &model, &words, 1, raw_callback,
                               NULL, &expected) == MVMC_KRYLOV_STATUS_OK,
          "reference evaluation across N");
    for (order = 0; order <= 2; ++order) {
      CHECK(close_complex(session_values[order], expected.value[order]),
            "session value equals reference");
    }
    if (root == 0) {
      /* Full filling is an eigenstate: E = 2U - 4 mu = 3 - 1. */
      CHECK(close_complex(session_values[1], 2.0 * table_value(0xf)),
            "full filling q1 by hand");
      CHECK(close_complex(session_values[2], 4.0 * table_value(0xf)),
            "full filling q2 by hand");
    }
    if (root == 2) {
      CHECK(close_complex(session_values[0], table_value(0)) &&
                session_values[1] == 0.0 && session_values[2] == 0.0,
            "vacuum is annihilated by H");
    }
    if (root == 5) {
      CHECK(mvmc_bounded_krylov_session_end(workspace) ==
                MVMC_KRYLOV_STATUS_OK,
            "session end");
      CHECK(bounded_values(workspace, roots[root], 0, cold_values) ==
                MVMC_KRYLOV_STATUS_OK,
            "cold evaluation after session");
      for (order = 0; order <= 2; ++order) {
        CHECK(close_complex(cold_values[order], session_values[order]),
              "cold equals cached session value");
      }
    }
  }
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
  mvmc_krylov_workspace_destroy(reference);
}

static const MVMCClassicKrylovTransfer kRawTransfers[6] = {
    {0, 0, 1, 0, 1.0 + 0.0 * I}, {1, 0, 0, 0, 1.0 + 0.0 * I},
    {0, 1, 1, 1, 1.0 + 0.0 * I}, {1, 1, 0, 1, 1.0 + 0.0 * I},
    {0, 0, 0, 0, 0.25 + 0.0 * I}, {1, 1, 1, 1, 0.25 + 0.0 * I}};
static const MVMCClassicKrylovSiteCoupling kRawIntra[2] = {{0, 1.5},
                                                           {1, 1.5}};
static const MVMCClassicKrylovPairCoupling kRawInter[1] = {{0, 1, 0.125}};
static const MVMCClassicKrylovPairCoupling kRawHund[1] = {{0, 1, 0.0625}};
static const MVMCClassicKrylovPairCoupling kRawExchange[1] = {{0, 1, 0.5}};

static MVMCClassicKrylovRawModel raw_model(int up, int down) {
  MVMCClassicKrylovRawModel raw;
  memset(&raw, 0, sizeof(raw));
  raw.site_count = 2;
  raw.up_electron_count = up;
  raw.down_electron_count = down;
  raw.transfer_count = 6;
  raw.transfers = kRawTransfers;
  raw.coulomb_intra_count = 2;
  raw.coulomb_intra = kRawIntra;
  raw.coulomb_inter_count = 1;
  raw.coulomb_inter = kRawInter;
  raw.hund_count = 1;
  raw.hund = kRawHund;
  return raw;
}

static uint64_t model_fingerprint(const MVMCKrylovFockModel *model) {
  uint64_t hash = UINT64_C(1469598103934665603);
  size_t index;
  const unsigned char *bytes;
#define MIX(value)                                       \
  do {                                                   \
    hash ^= (uint64_t)(value);                           \
    hash *= UINT64_C(1099511628211);                     \
  } while (0)
  MIX(model->site_count);
  MIX(model->up_electron_count);
  MIX(model->down_electron_count);
  MIX((unsigned)model->pure_spin);
  MIX((unsigned)model->domain);
  MIX(model->term_count);
  MIX(model->operator_count);
  for (index = 0; index < model->term_count; ++index) {
    size_t byte;
    bytes = (const unsigned char *)&model->terms[index].coefficient;
    for (byte = 0; byte < sizeof(double complex); ++byte) MIX(bytes[byte]);
    MIX(model->terms[index].operator_offset);
    MIX(model->terms[index].operator_count);
  }
  for (index = 0; index < model->operator_count; ++index) {
    MIX((unsigned)model->operators[index].kind);
    MIX(model->operators[index].orbital);
  }
#undef MIX
  return hash;
}

static int same_terms(const MVMCKrylovFockModel *left,
                      const MVMCKrylovFockModel *right) {
  size_t index;
  if (left->term_count != right->term_count ||
      left->operator_count != right->operator_count) {
    return 0;
  }
  for (index = 0; index < left->term_count; ++index) {
    if (left->terms[index].coefficient != right->terms[index].coefficient ||
        left->terms[index].operator_offset !=
            right->terms[index].operator_offset ||
        left->terms[index].operator_count !=
            right->terms[index].operator_count ||
        left->terms[index].source_kind != right->terms[index].source_kind ||
        left->terms[index].source_index != right->terms[index].source_index) {
      return 0;
    }
  }
  for (index = 0; index < left->operator_count; ++index) {
    if (left->operators[index].kind != right->operators[index].kind ||
        left->operators[index].orbital != right->operators[index].orbital) {
      return 0;
    }
  }
  return 1;
}

static void test_classic_model_domain(void) {
  const int rank = world_rank();
  MVMCClassicKrylovRawModel gc_raw = raw_model(0, 0);
  MVMCClassicKrylovRawModel fixed_raw = raw_model(1, 1);
  MVMCClassicKrylovModelWorkspace *gc = NULL;
  MVMCClassicKrylovModelWorkspace *fixed = NULL;
  MVMCClassicKrylovModelWorkspace *rejected = NULL;
  MVMCKrylovStatus status;

  /* Non-root ranks may pass raw == NULL. */
  status = mvmc_classic_krylov_model_workspace_create_domain_from_root(
      rank == 0 ? &gc_raw : NULL, MVMC_KRYLOV_DOMAIN_GC_BALANCED,
      world_communicator(), &gc);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && gc != NULL,
        "GC-balanced classic model");
  status = mvmc_classic_krylov_model_workspace_create_from_root(
      rank == 0 ? &fixed_raw : NULL, world_communicator(), &fixed);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && fixed != NULL,
        "legacy creator remains fixed-N");
  if (gc != NULL && fixed != NULL) {
    const MVMCKrylovFockModel *gc_model = mvmc_classic_krylov_model(gc);
    const MVMCKrylovFockModel *fixed_model = mvmc_classic_krylov_model(fixed);
    CHECK(gc_model->domain == MVMC_KRYLOV_DOMAIN_GC_BALANCED &&
              gc_model->up_electron_count == 0 &&
              gc_model->down_electron_count == 0 && gc_model->pure_spin == 0,
          "GC model metadata");
    CHECK(fixed_model->domain == MVMC_KRYLOV_DOMAIN_FIXED &&
              fixed_model->up_electron_count == 1,
          "fixed model metadata");
    CHECK(same_terms(gc_model, fixed_model),
          "GC and fixed share the Hamiltonian conversion");
#ifdef _mpi_use
    {
      unsigned long long local = (unsigned long long)model_fingerprint(gc_model);
      unsigned long long minimum = 0;
      unsigned long long maximum = 0;
      MPI_Allreduce(&local, &minimum, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN,
                    MPI_COMM_WORLD);
      MPI_Allreduce(&local, &maximum, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX,
                    MPI_COMM_WORLD);
      CHECK(minimum == maximum, "GC model identical on every rank");
    }
#else
    (void)model_fingerprint;
#endif
  }
  mvmc_classic_krylov_model_workspace_destroy(gc);
  mvmc_classic_krylov_model_workspace_destroy(fixed);

  status = mvmc_classic_krylov_model_workspace_create_domain_from_root(
      rank == 0 ? &fixed_raw : NULL, MVMC_KRYLOV_DOMAIN_GC_EVEN,
      world_communicator(), &rejected);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && rejected == NULL,
        "GC domain requires zero counts");
  status = mvmc_classic_krylov_model_workspace_create_domain_from_root(
      rank == 0 ? &gc_raw : NULL, MVMC_KRYLOV_DOMAIN_FIXED,
      world_communicator(), &rejected);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && rejected == NULL,
        "fixed domain keeps positive counts");
  status = mvmc_classic_krylov_model_workspace_create_domain_from_root(
      rank == 0 ? &gc_raw : NULL, (MVMCKrylovConfigurationDomain)9,
      world_communicator(), &rejected);
  CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && rejected == NULL,
        "invalid classic domain");
  {
    MVMCClassicKrylovRawModel spin = gc_raw;
    spin.pure_spin = 1;
    status = mvmc_classic_krylov_model_workspace_create_domain_from_root(
        rank == 0 ? &spin : NULL, MVMC_KRYLOV_DOMAIN_GC_BALANCED,
        world_communicator(), &rejected);
    CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && rejected == NULL,
          "GC rejects pure_spin");
  }
  {
    MVMCClassicKrylovRawModel exchange = gc_raw;
    exchange.exchange_count = 1;
    exchange.exchange = kRawExchange;
    status = mvmc_classic_krylov_model_workspace_create_domain_from_root(
        rank == 0 ? &exchange : NULL, MVMC_KRYLOV_DOMAIN_GC_EVEN,
        world_communicator(), &rejected);
    CHECK(status == MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL && rejected == NULL,
          "GC rejects Exchange");
  }
  if (world_size() > 1) {
    /* Ranks disagree on the domain: collective rejection, no hang. */
    status = mvmc_classic_krylov_model_workspace_create_domain_from_root(
        rank == 0 ? &gc_raw : NULL,
        rank == 0 ? MVMC_KRYLOV_DOMAIN_GC_BALANCED : MVMC_KRYLOV_DOMAIN_GC_EVEN,
        world_communicator(), &rejected);
    CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && rejected == NULL,
          "domain mismatch rejected collectively");
  }
}

int main(int argc, char **argv) {
  int total = 0;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  test_plan_example();
  test_mask_sets();
  test_model_validation();
  test_multiword();
  test_fixed_sector_proposal_rejects_gc();
  test_plan_hash();
  test_bounded_values_across_sectors();
  test_classic_model_domain();
#ifdef _mpi_use
  MPI_Allreduce(&failure_count, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  total = failure_count;
#endif
  if (world_rank() == 0) {
    printf(total == 0 ? "gc_krylov_domain_unit: PASS\n"
                      : "gc_krylov_domain_unit: FAIL\n");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return total == 0 ? 0 : 1;
}
