#include "power_lanczos_production_session.h"
#include "power_lanczos_sha256.h"
#include "power_lanczos_snapshot.h"
#include "power_lanczos_stabilization_output.h"
#include "power_lanczos_stabilization_statistics.h"

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failures = 0;
static int world_rank = 0;
static int world_size = 1;

#define CHECK(condition, ...)                                                 \
  do {                                                                        \
    if (!(condition)) {                                                       \
      fprintf(stderr, "PowerLanczosProductionSession_Unit FAIL rank %d: ",  \
              world_rank);                                                    \
      fprintf(stderr, __VA_ARGS__);                                           \
      fprintf(stderr, " (line %d)\n", __LINE__);                           \
      ++failures;                                                             \
    }                                                                         \
  } while (0)

static void test_sha256(void) {
  char digest[MVMC_POWER_LANCZOS_SHA256_HEX_CAPACITY];
  CHECK(mvmc_power_lanczos_sha256_hex(NULL, 0, digest) &&
            strcmp(digest,
                   "e3b0c44298fc1c149afbf4c8996fb924"
                   "27ae41e4649b934ca495991b7852b855") == 0,
        "empty SHA-256 vector");
  CHECK(mvmc_power_lanczos_sha256_hex("abc", 3, digest) &&
            strcmp(digest,
                   "ba7816bf8f01cfea414140de5dae2223"
                   "b00361a396177a9cb410ff61f20015ad") == 0,
        "abc SHA-256 vector");
  CHECK(!mvmc_power_lanczos_sha256_hex(NULL, 1, digest),
        "NULL nonempty SHA input rejected");
}

static void test_rng(void) {
  static const uint64_t expected[4][3] = {
      {UINT64_C(0xa1b9d5b340250210), UINT64_C(0x53504556c629bf73),
       UINT64_C(0xd70aa66b0f815041)},
      {UINT64_C(0xc43e59da24342f89), UINT64_C(0xa53492842b302584),
       UINT64_C(0xd6df8738705275b2)},
      {UINT64_C(0xa9a711b1d121d94c), UINT64_C(0x7fab24b597acc5a9),
       UINT64_C(0x68b77e2814fcb2ba)},
      {UINT64_C(0xb66d436255bc1080), UINT64_C(0xf1687c09b9a0166d),
       UINT64_C(0xdceef04ae2455571)}};
  const MVMCPowerLanczosRngStage stages[3] = {
      MVMC_POWER_LANCZOS_RNG_STAGE_BASE,
      MVMC_POWER_LANCZOS_RNG_STAGE_COEFFICIENT,
      MVMC_POWER_LANCZOS_RNG_STAGE_FINAL};
  MVMCPowerLanczosRngDomain domain;
  size_t rank;
  size_t stage;
  for (rank = 0; rank < 4; ++rank) {
    for (stage = 0; stage < 3; ++stage) {
      int matches = 0;
      CHECK(mvmc_power_lanczos_rng_derive(
                UINT64_C(0x503641434f4e4631), UINT64_C(7), stages[stage],
                rank, 4, 2, &domain) == MVMC_KRYLOV_STATUS_OK &&
                domain.valid && domain.split_rank == rank % 2 &&
                domain.derived_seed == expected[rank][stage],
            "frozen RNG vector rank=%zu stage=%zu actual=0x%016" PRIx64,
            rank, stage, domain.derived_seed);
      CHECK(mvmc_power_lanczos_rng_domain_matches(
                &domain, UINT64_C(0x503641434f4e4631), UINT64_C(7),
                stages[stage], rank, 4, 2, &matches) ==
                    MVMC_KRYLOV_STATUS_OK &&
                matches,
            "RNG domain self match rank=%zu stage=%zu", rank, stage);
      CHECK(mvmc_power_lanczos_rng_domain_matches(
                &domain, UINT64_C(0x503641434f4e4631), UINT64_C(8),
                stages[stage], rank, 4, 2, &matches) ==
                    MVMC_KRYLOV_STATUS_OK &&
                !matches,
            "RNG run-index mutation rank=%zu stage=%zu", rank, stage);
    }
  }
  CHECK(mvmc_power_lanczos_rng_derive(
            1, 0, MVMC_POWER_LANCZOS_RNG_STAGE_BASE, 0, 0, 1,
            &domain) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "zero world size rejected");
  CHECK(mvmc_power_lanczos_rng_derive(
            1, 0, MVMC_POWER_LANCZOS_RNG_STAGE_BASE, 2, 2, 1,
            &domain) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "out-of-range rank rejected");
  CHECK(mvmc_power_lanczos_rng_derive(
            1, 0, (MVMCPowerLanczosRngStage)99, 0, 1, 1,
            &domain) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "unknown RNG stage rejected");
}

static void test_snapshot_lifecycle(void) {
  static const unsigned char expected_sha256[32] = {
      0x48, 0x8a, 0x3d, 0x42, 0x08, 0x5b, 0xbd, 0x51,
      0xec, 0x9e, 0x71, 0x75, 0x88, 0x34, 0x45, 0x68,
      0xdc, 0xc1, 0x96, 0x7b, 0xca, 0x3a, 0xb0, 0xc3,
      0xa2, 0xb1, 0x9a, 0xfe, 0x1a, 0x6d, 0x48, 0x50};
  uint64_t source[2] = {UINT64_C(0x0123456789abcdef), UINT64_C(0x15)};
  const uint64_t original[2] = {source[0], source[1]};
  MVMCPowerLanczosSnapshot *snapshot = NULL;
  MVMCPowerLanczosSnapshotSummary initial;
  MVMCPowerLanczosSnapshotSummary completed;
  uint64_t *coefficient_words = NULL;
  uint64_t *final_words = NULL;
  size_t word_count = 0;
  CHECK(mvmc_power_lanczos_snapshot_create(
            source, 2, 69, UINT64_C(17), &snapshot) ==
            MVMC_KRYLOV_STATUS_OK,
        "snapshot create");
  source[0] = 0;
  source[1] = 0;
  CHECK(mvmc_power_lanczos_snapshot_summary(snapshot, &initial) ==
            MVMC_KRYLOV_STATUS_OK &&
            initial.valid && initial.word_count == 2 &&
            initial.configuration_bit_count == 69 &&
            memcmp(initial.base_sha256, expected_sha256, 32) == 0 &&
            memcmp(initial.base_sha256, initial.coefficient_sha256, 32) == 0 &&
            memcmp(initial.base_sha256, initial.final_sha256, 32) == 0,
        "initial snapshot summary");
  CHECK(mvmc_power_lanczos_snapshot_coefficient_begin(
            snapshot, &coefficient_words, &word_count) ==
                MVMC_KRYLOV_STATUS_OK &&
            coefficient_words != source && word_count == 2 &&
            coefficient_words[0] == original[0] &&
            coefficient_words[1] == original[1],
        "coefficient deep clone begin");
  coefficient_words[0] ^= UINT64_C(0x80);
  CHECK(mvmc_power_lanczos_snapshot_verify(snapshot) ==
            MVMC_KRYLOV_STATUS_OK,
        "active coefficient clone may mutate");
  CHECK(mvmc_power_lanczos_snapshot_coefficient_complete(snapshot) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_snapshot_summary(snapshot, &completed) ==
                MVMC_KRYLOV_STATUS_OK &&
            memcmp(completed.base_sha256, completed.coefficient_sha256,
                   32) != 0 &&
            memcmp(completed.base_sha256, completed.final_sha256, 32) == 0,
        "coefficient completion freezes clone A");
  CHECK(mvmc_power_lanczos_snapshot_final_begin(
            snapshot, &final_words, &word_count) ==
                MVMC_KRYLOV_STATUS_OK &&
            final_words != source && final_words != coefficient_words &&
            final_words[0] == original[0] && final_words[1] == original[1],
        "final clone B preserved and disjoint");
  final_words[0] ^= UINT64_C(0x40);
  CHECK(mvmc_power_lanczos_snapshot_verify(snapshot) ==
            MVMC_KRYLOV_STATUS_OK,
        "active final clone may mutate");
  CHECK(mvmc_power_lanczos_snapshot_final_complete(snapshot) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_snapshot_summary(snapshot, &completed) ==
                MVMC_KRYLOV_STATUS_OK &&
            completed.state ==
                MVMC_POWER_LANCZOS_SNAPSHOT_FINAL_COMPLETE &&
            memcmp(completed.base_sha256, completed.final_sha256, 32) != 0 &&
            mvmc_power_lanczos_snapshot_verify(snapshot) ==
                MVMC_KRYLOV_STATUS_OK,
        "final completion freezes clone B");
  final_words[0] ^= UINT64_C(0x10);
  CHECK(mvmc_power_lanczos_snapshot_verify(snapshot) ==
            MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE,
        "frozen final tamper fails loud");
  mvmc_power_lanczos_snapshot_destroy(snapshot);

  CHECK(mvmc_power_lanczos_snapshot_create(
            original, 2, 69, UINT64_C(17), &snapshot) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_snapshot_coefficient_begin(
                snapshot, &coefficient_words, &word_count) ==
                MVMC_KRYLOV_STATUS_OK,
        "coefficient tamper fixture recreate");
  coefficient_words[0] ^= UINT64_C(0x80);
  CHECK(mvmc_power_lanczos_snapshot_coefficient_complete(snapshot) ==
            MVMC_KRYLOV_STATUS_OK,
        "coefficient tamper fixture freeze");
  coefficient_words[0] ^= UINT64_C(0x20);
  CHECK(mvmc_power_lanczos_snapshot_verify(snapshot) ==
            MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE,
        "frozen coefficient tamper fails loud");
  mvmc_power_lanczos_snapshot_destroy(snapshot);
}

static void test_snapshot_failures(void) {
  const uint64_t valid[2] = {UINT64_C(3), UINT64_C(1)};
  const uint64_t bad_padding[2] = {UINT64_C(3), UINT64_C(2)};
  MVMCPowerLanczosSnapshot *snapshot = NULL;
  uint64_t *coefficient_words = NULL;
  size_t word_count = 0;
  CHECK(mvmc_power_lanczos_snapshot_create(
            bad_padding, 2, 65, 1, &snapshot) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "noncanonical padding rejected");
  CHECK(mvmc_power_lanczos_snapshot_create(
            valid, 2, 0, 1, &snapshot) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "zero bit count rejected");
  CHECK(mvmc_power_lanczos_snapshot_create(
            valid, 2, 65, 0, &snapshot) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "zero generation rejected");
  CHECK(mvmc_power_lanczos_snapshot_create(
            valid, 2, 65, 1, &snapshot) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_snapshot_coefficient_begin(
                snapshot, &coefficient_words, &word_count) ==
                MVMC_KRYLOV_STATUS_OK,
        "tamper fixture setup");
  CHECK(mvmc_power_lanczos_snapshot_testing_corrupt(
            snapshot, MVMC_POWER_LANCZOS_SNAPSHOT_TEST_FINAL, 0,
            UINT64_C(1)) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_power_lanczos_snapshot_coefficient_complete(snapshot) ==
                MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE,
        "clone B tamper detected before final");
  mvmc_power_lanczos_snapshot_destroy(snapshot);
}

static MVMCPowerLanczosProductionSessionConfig session_config(int step) {
  MVMCPowerLanczosProductionSessionConfig config;
  memset(&config, 0, sizeof(config));
  config.power_step = step;
  config.resolved_base_seed = UINT64_C(0x503641434f4e4631);
  config.run_index = UINT64_C(7);
  config.mpi_world_rank = (size_t)world_rank;
  config.mpi_world_size = (size_t)world_size;
  config.split_size = 2;
  config.base_generation = UINT64_C(19);
  config.configuration_bit_count = 65;
  return config;
}

static MVMCKrylovBoundedCommunicator test_world_communicator(void) {
#ifdef _mpi_use
  return MPI_COMM_WORLD;
#else
  return 0;
#endif
}

static size_t test_particle_count(const uint64_t *words, size_t word_count) {
  size_t total = 0;
  size_t index;
  if (words == NULL) return 0;
  for (index = 0; index < word_count; ++index) {
    uint64_t word = words[index];
    while (word != 0) {
      word &= word - UINT64_C(1);
      ++total;
    }
  }
  return total;
}

static MVMCKrylovStatus create_test_session(
    const MVMCPowerLanczosProductionSessionConfig *config,
    const uint64_t *canonical_words, size_t word_count,
    MVMCPowerLanczosProductionSession **session) {
  MVMCPowerLanczosTerminalSnapshotInput terminal;
  const MVMCPowerLanczosTerminalSnapshotInput *root_terminal = NULL;
  memset(&terminal, 0, sizeof(terminal));
  if (config != NULL && config->mpi_world_rank == 0) {
    terminal.canonical_words = canonical_words;
    terminal.word_count = word_count;
    terminal.configuration_bit_count = config->configuration_bit_count;
    terminal.generation = config->base_generation;
    terminal.terminal_proposal_counter = UINT64_C(4096);
    terminal.expected_particle_count =
        test_particle_count(canonical_words, word_count);
    terminal.log_sampling_support = 0.0;
    root_terminal = &terminal;
  }
  return mvmc_power_lanczos_production_session_create(
      config, test_world_communicator(), root_terminal, session);
}

static void test_session(void) {
  const uint64_t words[2] = {UINT64_C(0x55aa), UINT64_C(1)};
  int step;
  for (step = 1; step <= 2; ++step) {
    MVMCPowerLanczosProductionSessionConfig config = session_config(step);
    MVMCPowerLanczosProductionSession *session = NULL;
    MVMCPowerLanczosProductionSessionSummary summary;
    CHECK(create_test_session(
              &config, words, 2, &session) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_production_session_verify(session) ==
                  MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_production_session_summary(
                  session, &summary) == MVMC_KRYLOV_STATUS_OK &&
              summary.valid &&
              summary.state ==
                  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_SNAPSHOTS_VERIFIED &&
              summary.basis_count == (size_t)step + 1 &&
              summary.upper_count == (step == 1 ? 3u : 6u) &&
              summary.snapshot_selection.valid &&
              summary.snapshot_selection.version ==
                  MVMC_POWER_LANCZOS_TERMINAL_SNAPSHOT_SELECTION_VERSION &&
              summary.snapshot_selection.selection_kind ==
                  MVMC_POWER_LANCZOS_TERMINAL_SNAPSHOT_WORLD_RANK_0_BASE_TERMINAL_AFTER_LAST_PROPOSAL &&
              summary.snapshot_selection.owner_world_rank == 0 &&
              summary.snapshot_selection.terminal_proposal_counter ==
                  UINT64_C(4096) &&
              summary.snapshot_selection.generation ==
                  summary.snapshot.generation &&
              summary.snapshot_selection.word_count ==
                  summary.snapshot.word_count &&
              summary.snapshot_selection.configuration_bit_count ==
                  summary.snapshot.configuration_bit_count &&
              summary.snapshot_selection.expected_particle_count == 9 &&
              memcmp(summary.snapshot_selection.snapshot_sha256,
                     summary.snapshot.base_sha256,
                     MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES) == 0 &&
              summary.base_rng.derived_seed !=
                  summary.coefficient_rng.derived_seed &&
              summary.base_rng.derived_seed != summary.final_rng.derived_seed &&
              summary.coefficient_rng.derived_seed !=
                  summary.final_rng.derived_seed,
          "session step %d create/verify", step);
    mvmc_power_lanczos_production_session_destroy(session);
  }
  {
    MVMCPowerLanczosProductionSessionConfig config = session_config(3);
    MVMCPowerLanczosProductionSession *session = NULL;
    CHECK(create_test_session(
              &config, words, 2, &session) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              session == NULL,
          "step 3 rejected");
  }
}

static void test_terminal_snapshot_failures(void) {
  const uint64_t words[2] = {UINT64_C(0x55aa), UINT64_C(1)};
  MVMCPowerLanczosProductionSessionConfig config = session_config(1);
  MVMCPowerLanczosTerminalSnapshotInput terminal;
  const MVMCPowerLanczosTerminalSnapshotInput *root_terminal;
  MVMCPowerLanczosProductionSession *session = NULL;
  memset(&terminal, 0, sizeof(terminal));
  terminal.canonical_words = words;
  terminal.word_count = 2;
  terminal.configuration_bit_count = config.configuration_bit_count;
  terminal.generation = config.base_generation;
  terminal.terminal_proposal_counter = UINT64_C(4096);
  terminal.expected_particle_count = 10;
  terminal.log_sampling_support = 0.0;
  root_terminal = world_rank == 0 ? &terminal : NULL;
  CHECK(mvmc_power_lanczos_production_session_create(
            &config, test_world_communicator(), root_terminal, &session) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            session == NULL,
        "root particle-count mismatch rejected collectively");

  terminal.expected_particle_count = 9;
  terminal.log_sampling_support = NAN;
  CHECK(mvmc_power_lanczos_production_session_create(
            &config, test_world_communicator(), root_terminal, &session) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            session == NULL,
        "nonfinite root sampling support rejected collectively");

  terminal.log_sampling_support = 0.0;
  terminal.terminal_proposal_counter = 0;
  CHECK(mvmc_power_lanczos_production_session_create(
            &config, test_world_communicator(), root_terminal, &session) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            session == NULL,
        "zero terminal proposal counter rejected collectively");

  terminal.terminal_proposal_counter = UINT64_C(4096);
  terminal.generation = config.base_generation + UINT64_C(1);
  CHECK(mvmc_power_lanczos_production_session_create(
            &config, test_world_communicator(), root_terminal, &session) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            session == NULL,
        "root generation mismatch rejected collectively");

  terminal.generation = config.base_generation;
  if (world_rank == world_size - 1) {
    config.mpi_world_rank = config.mpi_world_size;
  }
  CHECK(mvmc_power_lanczos_production_session_create(
            &config, test_world_communicator(), root_terminal, &session) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            session == NULL,
        "rank-local communicator identity mismatch rejected collectively");
  config = session_config(1);

#ifdef _mpi_use
  CHECK(mvmc_power_lanczos_production_session_create(
            &config, MPI_COMM_NULL, root_terminal, &session) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            session == NULL,
        "null world communicator rejected before collective entry");
  if (world_size > 1) {
    CHECK(mvmc_power_lanczos_production_session_create(
              &config, test_world_communicator(), &terminal, &session) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              session == NULL,
          "non-root candidate input rejected collectively");
    CHECK(mvmc_power_lanczos_production_session_testing_corrupt_terminal_after_broadcast(
              1, 0, UINT64_C(1)) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_production_session_create(
                  &config, test_world_communicator(), root_terminal,
                  &session) ==
                  MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE &&
              session == NULL,
          "rank-local post-broadcast substitution rejected collectively");
  }
#endif
}

static MVMCKrylovStatus unused_amplitude(
    const uint64_t *configuration_words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
  (void)configuration_words;
  (void)word_count;
  (void)context;
  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(result, 0, sizeof(*result));
  return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
}

static MVMCPowerLanczosProductionExecution execution_fixture(
    int step, const MVMCKrylovFockModel *model,
    MVMCKrylovBoundedCommunicator chain_communicator) {
  MVMCPowerLanczosProductionExecution execution;
  size_t basis_count = (size_t)step + 1;
  size_t index;
  memset(&execution, 0, sizeof(execution));
#ifdef _mpi_use
  execution.world_communicator = MPI_COMM_WORLD;
#else
  execution.world_communicator = 0;
#endif
  execution.chain_communicator = chain_communicator;
  execution.proposal_model = model;
  execution.amplitude = unused_amplitude;
  execution.amplitude_generation_hash = UINT64_C(0x5036444558454341);
  execution.bounded_limits.amplitude_policy_hash =
      UINT64_C(0x503644424f554e44);
  execution.bounded_limits.cache_bytes = 4096;
  execution.bounded_limits.max_row_transitions = 16;
  execution.bounded_limits.max_workspace_bytes = (size_t)1024 * 1024;
  execution.bounded_limits.max_node_expansions = UINT64_C(10000);
  execution.bounded_limits.max_terminal_amplitude_calls = UINT64_C(10000);
  execution.bounded_limits.max_total_row_transitions = UINT64_C(100000);
  execution.bounded_limits.max_order = (int)basis_count;
  execution.coefficient_guide_policy.order = (int)basis_count;
  execution.coefficient_guide_policy.eta = 0.25;
  execution.coefficient_guide_policy.policy_hash =
      UINT64_C(0x5036444755494445) ^ (uint64_t)step;
  execution.matrix_policy.order = (int)basis_count;
  execution.matrix_policy.eta = 0.25;
  for (index = 0; index <= basis_count; ++index) {
    execution.coefficient_guide_policy.lambda[index] = 1.0;
    execution.matrix_policy.guide_lambda[index] = 1.0;
    execution.matrix_policy.target_weight[index] = 1.0;
  }
  CHECK(mvmc_krylov_positive_sampler_proposal_policy_create(
            1, 2, &execution.proposal_policy) == MVMC_KRYLOV_STATUS_OK,
        "execution proposal policy fixture");
  CHECK(mvmc_power_lanczos_gevp_default_policy(
            0x1p-40, &execution.gevp_policy) == MVMC_KRYLOV_GEVP_OK,
        "execution GEVP policy fixture");
  execution.block_count = 2;
  execution.coefficient_sample_count = 4;
  execution.coefficient_interval = 1;
  execution.final_sample_count = 4;
  execution.final_interval = 1;
  execution.maximum_leave_one_projective_distance = 0.25;
  return execution;
}

static void test_execution_prepare(void) {
  static const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0}};
  static const MVMCKrylovHamiltonianTerm terms[] = {
      {1.0, 0, 2, 1, 0}, {1.0, 2, 2, 1, 1}};
  MVMCKrylovFockModel model;
  const uint64_t words[1] = {UINT64_C(2)};
  size_t split_size = world_size > 1 ? 2u : 1u;
  int step;
#ifdef _mpi_use
  MPI_Comm chain_communicator = MPI_COMM_NULL;
  CHECK(MPI_Comm_split(MPI_COMM_WORLD,
                       world_rank / (int)split_size, world_rank,
                       &chain_communicator) == MPI_SUCCESS,
        "execution chain communicator fixture");
#else
  MVMCKrylovBoundedCommunicator chain_communicator = 0;
#endif
  memset(&model, 0, sizeof(model));
  model.site_count = 2;
  model.up_electron_count = 1;
  model.hermitian = 1;
  model.terms = terms;
  model.term_count = 2;
  model.operators = operators;
  model.operator_count = 4;
  for (step = 1; step <= 2; ++step) {
    MVMCPowerLanczosProductionSessionConfig config = session_config(step);
    MVMCPowerLanczosProductionExecution execution =
        execution_fixture(step, &model, chain_communicator);
    MVMCPowerLanczosProductionSession *session = NULL;
    MVMCPowerLanczosProductionSessionSummary before;
    MVMCPowerLanczosProductionSessionSummary after;
    size_t expected_chain_size = split_size;
    config.split_size = split_size;
    config.configuration_bit_count = 4;
    if ((size_t)world_rank / split_size * split_size + expected_chain_size >
        (size_t)world_size) {
      expected_chain_size =
          (size_t)world_size - (size_t)world_rank / split_size * split_size;
    }
    CHECK(create_test_session(
              &config, words, 1, &session) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_production_session_summary(
                  session, &before) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_production_session_prepare_execution(
                  session, &execution) == MVMC_KRYLOV_STATUS_OK,
          "execution prepare step %d", step);
    execution.block_count = 99;
    CHECK(mvmc_power_lanczos_production_session_summary(session, &after) ==
                  MVMC_KRYLOV_STATUS_OK &&
              after.valid && after.execution_prepared &&
              after.state ==
                  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_SNAPSHOTS_VERIFIED &&
              after.chain_rank == (size_t)world_rank % split_size &&
              after.chain_size == expected_chain_size &&
              after.chain_count ==
                  ((size_t)world_size + split_size - 1) / split_size &&
              after.block_count == 2 &&
              after.allocated_bytes > before.allocated_bytes &&
              isnan(after.maximum_leave_one_projective_distance) &&
              isnan(after.final_energy) &&
              isnan(after.final_energy_imaginary) &&
              isnan(after.corrected_variance) &&
              isnan(after.corrected_variance_imaginary) &&
              isnan(after.final_local_energy_abs_squared_diagnostic),
          "execution prepared summary step %d", step);
    CHECK(mvmc_power_lanczos_production_session_prepare_execution(
              session, &execution) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "execution prepare re-entry rejected step %d", step);
    mvmc_power_lanczos_production_session_destroy(session);
  }
  {
    MVMCPowerLanczosProductionSessionConfig config = session_config(1);
    MVMCPowerLanczosProductionExecution execution =
        execution_fixture(1, &model, chain_communicator);
    MVMCPowerLanczosProductionSession *session = NULL;
    MVMCPowerLanczosProductionSessionSummary summary;
    config.split_size = split_size;
    config.configuration_bit_count = 4;
    execution.coefficient_sample_count = 3;
    CHECK(create_test_session(
              &config, words, 1, &session) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_production_session_prepare_execution(
                  session, &execution) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              mvmc_power_lanczos_production_session_summary(
                  session, &summary) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              !summary.valid &&
              summary.state == MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FAILED,
          "invalid block partition fails session");
    mvmc_power_lanczos_production_session_destroy(session);
  }
  {
    MVMCPowerLanczosProductionSessionConfig config = session_config(1);
    MVMCPowerLanczosProductionExecution execution =
        execution_fixture(1, &model, chain_communicator);
    MVMCPowerLanczosProductionSession *session = NULL;
    MVMCPowerLanczosProductionSessionSummary summary;
    config.split_size = split_size;
    config.configuration_bit_count = 4;
    model.down_electron_count = 1;
    CHECK(create_test_session(
              &config, words, 1, &session) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_production_session_prepare_execution(
                  session, &execution) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              mvmc_power_lanczos_production_session_summary(
                  session, &summary) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              !summary.valid &&
              summary.state == MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FAILED,
          "terminal particle sector mismatch fails collectively");
    model.down_electron_count = 0;
    mvmc_power_lanczos_production_session_destroy(session);
  }
#ifdef _mpi_use
  if (world_size > 1) {
    MVMCPowerLanczosProductionSessionConfig config = session_config(1);
    MVMCPowerLanczosProductionExecution execution =
        execution_fixture(1, &model, chain_communicator);
    MVMCPowerLanczosProductionSession *session = NULL;
    MVMCPowerLanczosProductionSessionSummary summary;
    config.split_size = split_size;
    config.configuration_bit_count = 4;
    if (world_rank == world_size - 1) {
      execution.maximum_leave_one_projective_distance = 0.5;
    }
    CHECK(create_test_session(
              &config, words, 1, &session) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_production_session_prepare_execution(
                  session, &execution) ==
                  MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE &&
              mvmc_power_lanczos_production_session_summary(
                  session, &summary) ==
                  MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE &&
              !summary.valid &&
              summary.state == MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FAILED,
          "one-rank execution fingerprint mismatch fails all ranks");
    mvmc_power_lanczos_production_session_destroy(session);
  }
  if (chain_communicator != MPI_COMM_NULL) {
    MPI_Comm_free(&chain_communicator);
  }
#endif
}

static MVMCKrylovStatus three_site_amplitude(
    const uint64_t *configuration_words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
  static const double complex value[3] = {
      1.0 + 0.0 * I, 0.7 + 0.2 * I, 1.2 - 0.1 * I};
  uint64_t occupied;
  size_t site;
  (void)context;
  if (configuration_words == NULL || word_count != 1 || result == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  occupied = configuration_words[0] & UINT64_C(7);
  if ((configuration_words[0] & ~UINT64_C(63)) != 0 || occupied == 0 ||
      (occupied & (occupied - 1)) != 0 ||
      (configuration_words[0] & UINT64_C(56)) != 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  site = occupied == UINT64_C(1) ? 0u
         : occupied == UINT64_C(2) ? 1u
                                  : 2u;
  memset(result, 0, sizeof(*result));
  if (mvmc_scaled_complex_from_raw(value[site], &result->value) !=
      MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  result->well_pivoted_component_count = 1;
  result->local_factorization_count = 1;
  result->global_factorization_count = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

typedef struct {
  int failure_rank;
  uint64_t failure_call;
  uint64_t calls;
} AmplitudeFailureContext;

static MVMCKrylovStatus failing_three_site_amplitude(
    const uint64_t *configuration_words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
  AmplitudeFailureContext *failure = (AmplitudeFailureContext *)context;
  if (failure == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  ++failure->calls;
  if (world_rank == failure->failure_rank &&
      failure->calls == failure->failure_call) {
    if (result != NULL) memset(result, 0, sizeof(*result));
    return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  return three_site_amplitude(configuration_words, word_count, NULL, result);
}

static MVMCKrylovFockModel three_site_model(void) {
  static const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 2},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
      {MVMC_KRYLOV_FERMION_CREATE, 2},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1}};
  static const MVMCKrylovHamiltonianTerm terms[] = {
      {0.2, 0, 2, 1, 0},  {0.7, 2, 2, 1, 1},
      {1.3, 4, 2, 1, 2},  {-0.4, 6, 2, 1, 3},
      {-0.4, 8, 2, 1, 4}, {-0.3, 10, 2, 1, 5},
      {-0.3, 12, 2, 1, 6}};
  MVMCKrylovFockModel model;
  memset(&model, 0, sizeof(model));
  model.site_count = 3;
  model.up_electron_count = 1;
  model.hermitian = 1;
  model.terms = terms;
  model.term_count = sizeof(terms) / sizeof(terms[0]);
  model.operators = operators;
  model.operator_count = sizeof(operators) / sizeof(operators[0]);
  return model;
}

static void test_execution_full_lifecycle(void) {
  const MVMCKrylovFockModel model = three_site_model();
  const uint64_t words[1] = {UINT64_C(1)};
  const size_t split_size = world_size > 1 ? 2u : 1u;
  int step;
#ifdef _mpi_use
  MPI_Comm chain_communicator = MPI_COMM_NULL;
  CHECK(MPI_Comm_split(MPI_COMM_WORLD,
                       world_rank / (int)split_size, world_rank,
                       &chain_communicator) == MPI_SUCCESS,
        "full lifecycle chain communicator fixture");
#else
  MVMCKrylovBoundedCommunicator chain_communicator = 0;
#endif
  for (step = 1; step <= 2; ++step) {
    MVMCPowerLanczosProductionSessionConfig config = session_config(step);
    MVMCPowerLanczosProductionExecution execution =
        execution_fixture(step, &model, chain_communicator);
    MVMCPowerLanczosProductionSession *session = NULL;
    MVMCPowerLanczosProductionSessionSummary summary;
    const uint64_t coefficient_proposals = UINT64_C(8 + 16 + 16 + 512);
    const uint64_t final_proposals = UINT64_C(16 + 1024);
    size_t scale_index;
    MVMCKrylovStatus final_status;
    MVMCKrylovStatus finalize_status;
    MVMCKrylovStatus summary_status;
    config.split_size = split_size;
    config.configuration_bit_count = 6;
    execution.amplitude = three_site_amplitude;
    execution.bounded_limits.max_row_transitions = 32;
    execution.bounded_limits.max_node_expansions = UINT64_C(100000);
    execution.bounded_limits.max_terminal_amplitude_calls = UINT64_C(100000);
    execution.bounded_limits.max_total_row_transitions = UINT64_C(1000000);
    execution.scale_pilot_enabled = 1;
    execution.scale_pilot_warm_up = 8;
    execution.scale_pilot_sample_count = 16;
    execution.eta_relative = 0x1p-40;
    execution.coefficient_guide_policy.eta = 0x1p-40;
    execution.matrix_policy.eta = 0x1p-40;
    for (scale_index = 0; scale_index <= (size_t)step + 1;
         ++scale_index) {
      execution.coefficient_guide_policy.lambda[scale_index] =
          scale_index == 0 ? 0x1p-40 : 1.0;
      execution.coefficient_guide_policy.log_basis_scale[scale_index] = 0.0;
      execution.matrix_policy.guide_lambda[scale_index] =
          scale_index == 0 ? 0x1p-40 : 1.0;
      execution.matrix_policy.target_weight[scale_index] = 1.0;
      execution.matrix_policy.log_basis_scale[scale_index] = 0.0;
    }
    execution.coefficient_warm_up = 16;
    execution.coefficient_sample_count = 512;
    execution.final_warm_up = 16;
    execution.final_sample_count = 1024;
    execution.block_count = 4;
    execution.maximum_leave_one_projective_distance = 1.0;
    CHECK(create_test_session(
              &config, words, 1, &session) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_production_session_prepare_execution(
                  session, &execution) == MVMC_KRYLOV_STATUS_OK,
          "full lifecycle prepare step %d", step);
    CHECK(mvmc_power_lanczos_production_session_freeze_coefficient(session) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              mvmc_power_lanczos_production_session_run_final(session) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              mvmc_power_lanczos_production_session_finalize(session) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "full lifecycle stage skip rejected step %d", step);
    CHECK(mvmc_power_lanczos_production_session_run_coefficient(session) ==
                  MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_production_session_summary(
                  session, &summary) == MVMC_KRYLOV_STATUS_OK &&
              summary.state ==
                  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_READY &&
              summary.scale_pilot_enabled == 1 &&
              summary.scale_pilot_warm_up == UINT64_C(8) &&
              summary.scale_pilot_sample_count_per_chain == UINT64_C(16) &&
              summary.scale_pilot_proposals == UINT64_C(24) &&
              summary.scale_pilot_sample_count ==
                  UINT64_C(16) * (uint64_t)summary.chain_count &&
              summary.eta_relative == 0x1p-40 &&
              summary.resolved_eta ==
                  0x1p-40 * (double)((size_t)step + 2) &&
              summary.coefficient_proposals == coefficient_proposals &&
              summary.coefficient_sample_count ==
                  UINT64_C(512) * (uint64_t)summary.chain_count,
          "full lifecycle coefficient step %d status=%s", step,
          session == NULL ? "null" : "created");
    {
      MVMCPowerLanczosPrimitiveTraceSummary trace_summary;
      const MVMCKrylovStatus trace_status =
          mvmc_power_lanczos_production_session_coefficient_trace_summary(
              session, &trace_summary);
      if (world_rank == 0) {
        double complex values[25];
        double bounds[25];
        uint8_t support[25];
        double tail = NAN;
        MVMCPowerLanczosPrimitiveTraceGroup group;
        CHECK(trace_status == MVMC_KRYLOV_STATUS_OK &&
                  trace_summary.valid && trace_summary.frozen &&
                  trace_summary.stage ==
                      MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_COEFFICIENT &&
                  trace_summary.group_count == summary.chain_count &&
                  trace_summary.samples_per_group == 512 &&
                  trace_summary.primitive_count ==
                      1 + 4 * summary.upper_count,
              "coefficient trace summary step %d", step);
        CHECK(mvmc_power_lanczos_production_session_coefficient_trace_group(
                  session, summary.chain_count - 1, &group) ==
                      MVMC_KRYLOV_STATUS_OK &&
                  group.valid &&
                  group.group_ordinal == summary.chain_count - 1 &&
                  group.leader_world_rank ==
                      (summary.chain_count - 1) * split_size &&
                  group.chain_size ==
                      (size_t)world_size - group.leader_world_rank &&
                  group.sample_count == 512,
              "coefficient trace group identity step %d", step);
        CHECK(mvmc_power_lanczos_production_session_coefficient_trace_sample(
                  session, summary.chain_count - 1, 511, values, 25,
                  bounds, 25, support, 25, &tail) ==
                      MVMC_KRYLOV_STATUS_OK &&
                  isfinite(creal(values[0])) &&
                  isfinite(cimag(values[0])) && bounds[0] >= 0.0 &&
                  support[0] != 0 && tail == 0.0,
              "coefficient trace sample step %d", step);
      } else {
        CHECK(trace_status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
              "coefficient trace is root-owned step %d", step);
      }
    }
    CHECK(mvmc_power_lanczos_production_session_run_coefficient(session) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              mvmc_power_lanczos_production_session_run_final(session) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "coefficient re-entry/final skip rejected step %d", step);
    CHECK(mvmc_power_lanczos_production_session_freeze_coefficient(session) ==
                  MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_production_session_summary(
                  session, &summary) == MVMC_KRYLOV_STATUS_OK &&
              summary.state ==
                  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_COEFFICIENT_FROZEN &&
              summary.coefficient_provenance_hash != 0 &&
              summary.final_policy_hash != 0 &&
              summary.coefficient_gevp.valid &&
              summary.coefficient_gevp.retained_rank == step + 1 &&
              summary.maximum_leave_one_projective_distance >= 0.0 &&
              summary.maximum_leave_one_projective_distance <= 1.0,
          "full lifecycle freeze step %d", step);
    CHECK(mvmc_power_lanczos_production_session_freeze_coefficient(session) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "freeze re-entry rejected step %d", step);
    {
      MVMCPowerLanczosProductionMatrixResult frozen_matrix;
      MVMCPowerLanczosProductionCoefficientBlockResult frozen_block;
      MVMCPowerLanczosProductionFinalBlockResult premature_final_block;
      memset(&premature_final_block, 0xff,
             sizeof(premature_final_block));
      CHECK(mvmc_power_lanczos_production_session_matrix_result(
                session, &frozen_matrix) == MVMC_KRYLOV_STATUS_OK &&
                frozen_matrix.valid &&
                frozen_matrix.sample_count ==
                    summary.coefficient_sample_count,
            "coefficient-frozen matrix evidence step %d", step);
      CHECK(mvmc_power_lanczos_production_session_coefficient_block_result(
                session, 0, &frozen_block) == MVMC_KRYLOV_STATUS_OK &&
                frozen_block.valid && frozen_block.sample_count > 0,
            "coefficient-frozen block evidence step %d", step);
      CHECK(mvmc_power_lanczos_production_session_final_block_result(
                session, 0, &premature_final_block) ==
                    MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                !premature_final_block.valid,
            "final evidence unavailable before final stage step %d", step);
    }
    final_status =
        mvmc_power_lanczos_production_session_run_final(session);
    if (final_status == MVMC_KRYLOV_STATUS_OK) {
      MVMCPowerLanczosProductionFinalBlockResult ready_final_block;
      CHECK(mvmc_power_lanczos_production_session_final_block_result(
                session, 0, &ready_final_block) == MVMC_KRYLOV_STATUS_OK &&
                ready_final_block.valid &&
                ready_final_block.sample_count > 0,
            "final-ready block evidence step %d", step);
      {
        MVMCPowerLanczosPrimitiveTraceSummary trace_summary;
        const MVMCKrylovStatus trace_status =
            mvmc_power_lanczos_production_session_final_trace_summary(
                session, &trace_summary);
        if (world_rank == 0) {
          double complex values[2];
          double bounds[2];
          uint8_t support[2];
          double tail = NAN;
          MVMCPowerLanczosPrimitiveTraceGroup group;
          CHECK(trace_status == MVMC_KRYLOV_STATUS_OK &&
                    trace_summary.valid && trace_summary.frozen &&
                    trace_summary.stage ==
                        MVMC_POWER_LANCZOS_PRIMITIVE_TRACE_FINAL &&
                    trace_summary.group_count == summary.chain_count &&
                    trace_summary.samples_per_group == 1024 &&
                    trace_summary.primitive_count == 2,
                "final trace summary step %d", step);
          CHECK(mvmc_power_lanczos_production_session_final_trace_group(
                    session, summary.chain_count - 1, &group) ==
                        MVMC_KRYLOV_STATUS_OK &&
                    group.valid &&
                    group.group_ordinal == summary.chain_count - 1 &&
                    group.leader_world_rank ==
                        (summary.chain_count - 1) * split_size &&
                    group.chain_size ==
                        (size_t)world_size - group.leader_world_rank &&
                    group.sample_count == 1024,
                "final trace group identity step %d", step);
          CHECK(mvmc_power_lanczos_production_session_final_trace_sample(
                    session, summary.chain_count - 1, 1023, values, 2,
                    bounds, 2, support, 2, &tail) ==
                        MVMC_KRYLOV_STATUS_OK &&
                    values[0] == 1.0 && bounds[0] == 0.0 &&
                    support[0] ==
                        MVMC_POWER_LANCZOS_PRIMITIVE_SUPPORT_FINITE_NONZERO &&
                    isfinite(creal(values[1])) &&
                    isfinite(cimag(values[1])) && bounds[1] >= 0.0 &&
                    tail == hypot(creal(values[1]), cimag(values[1])),
                "final trace sample step %d", step);
        } else {
          CHECK(trace_status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "final trace is root-owned step %d", step);
        }
      }
    }
    finalize_status = final_status == MVMC_KRYLOV_STATUS_OK
                          ? mvmc_power_lanczos_production_session_finalize(
                                session)
                          : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    summary_status =
        mvmc_power_lanczos_production_session_summary(session, &summary);
    CHECK(final_status == MVMC_KRYLOV_STATUS_OK &&
              finalize_status == MVMC_KRYLOV_STATUS_OK &&
              summary_status == MVMC_KRYLOV_STATUS_OK &&
              summary.state ==
                  MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINALIZED &&
              summary.final_proposals == final_proposals &&
              summary.scale_pilot_accepted_steps <=
                  summary.scale_pilot_proposals * summary.chain_count &&
              summary.coefficient_accepted_steps <=
                  (summary.coefficient_proposals -
                   summary.scale_pilot_proposals) *
                      summary.chain_count &&
              summary.final_accepted_steps <=
                  summary.final_proposals * summary.chain_count &&
              summary.final_sample_count ==
                  UINT64_C(1024) * (uint64_t)summary.chain_count &&
              isfinite(summary.final_energy) &&
              isfinite(summary.final_energy_imaginary) &&
              isfinite(summary.corrected_variance) &&
              isfinite(summary.corrected_variance_imaginary) &&
              isfinite(summary.final_local_energy_abs_squared_diagnostic),
          "full lifecycle final step %d run=%s finalize=%s summary=%s "
          "state=%s E=%.17g var=%.17g diag=%.17g",
          step, mvmc_krylov_status_string(final_status),
          mvmc_krylov_status_string(finalize_status),
          mvmc_krylov_status_string(summary_status),
          mvmc_power_lanczos_production_session_state_string(summary.state),
          summary.final_energy, summary.corrected_variance,
          summary.final_local_energy_abs_squared_diagnostic);
    CHECK(fabs(summary.final_energy - summary.coefficient_gevp.energy) <=
              1.0e-2,
          "scaled GEVP/final-chain energy mismatch step %d gevp=%.17g "
          "final=%.17g",
          step, summary.coefficient_gevp.energy, summary.final_energy);
    if (finalize_status == MVMC_KRYLOV_STATUS_OK) {
      MVMCPowerLanczosProductionMatrixResult matrix;
      MVMCPowerLanczosProductionCoefficientBlockResult coefficient_block;
      MVMCPowerLanczosProductionFinalBlockResult final_block;
      MVMCPowerLanczosStabilizationResult stabilization;
      CHECK(summary.version ==
                MVMC_POWER_LANCZOS_PRODUCTION_SESSION_VERSION &&
                summary.version == UINT64_C(6),
            "immutable evidence session version step %d", step);
      for (scale_index = 0; scale_index <= (size_t)step + 1;
           ++scale_index) {
        CHECK(isfinite(summary.log_basis_scale[scale_index]),
              "scale pilot finite basis scale step %d index %zu", step,
              scale_index);
      }
      CHECK(fabs(summary.log_basis_scale[0] -
                 summary.log_basis_scale[1]) > 1.0e-6,
            "scale pilot must exercise nonuniform basis coordinates step %d",
            step);
      CHECK(mvmc_power_lanczos_production_session_matrix_result(
                session, &matrix) == MVMC_KRYLOV_STATUS_OK &&
                matrix.valid && matrix.basis_count == (size_t)step + 1 &&
                matrix.upper_count == (step == 1 ? 3u : 6u) &&
                matrix.sample_count == summary.coefficient_sample_count,
            "immutable matrix result step %d", step);
      CHECK(mvmc_power_lanczos_production_session_coefficient_block_result(
                session, 0, &coefficient_block) ==
                    MVMC_KRYLOV_STATUS_OK &&
                coefficient_block.valid &&
                coefficient_block.sample_count ==
                    UINT64_C(128) * (uint64_t)summary.chain_count &&
                coefficient_block.leave_one_sample_count ==
                    summary.coefficient_sample_count -
                        coefficient_block.sample_count &&
                isfinite(coefficient_block.leave_one_energy) &&
                isfinite(creal(
                    coefficient_block.leave_one_second_moment)) &&
                isfinite(cimag(
                    coefficient_block.leave_one_second_moment)) &&
                coefficient_block.leave_one_projective_distance >= 0.0 &&
                coefficient_block.leave_one_projective_distance <= 1.0,
            "immutable coefficient block step %d", step);
      CHECK(mvmc_power_lanczos_production_session_final_block_result(
                session, execution.block_count - 1, &final_block) ==
                    MVMC_KRYLOV_STATUS_OK &&
                final_block.valid &&
                final_block.sample_count ==
                    UINT64_C(256) * (uint64_t)summary.chain_count &&
                isfinite(creal(final_block.energy_sum)) &&
                isfinite(cimag(final_block.energy_sum)),
            "immutable final block step %d", step);
      if (world_rank == 0) {
        static const MVMCPowerLanczosStabilizationOutputIdentity identity = {
            "s2-focused-session",
            "0123456789012345678901234567890123456789",
            "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
            "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
            "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
            "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
            "appleclang-focused",
            "STABILIZATION_SEED_1"};
        char first_output[32768];
        char second_output[32768];
        char too_small[8];
        size_t first_size = 0;
        size_t second_size = 0;
        size_t too_small_size = 99;
        MVMCKrylovStatus first_render_status;
        MVMCKrylovStatus second_render_status;
        CHECK(mvmc_power_lanczos_stabilization_statistics_evaluate(
                  session, 0, &stabilization) == MVMC_KRYLOV_STATUS_OK &&
                  stabilization.valid &&
                  stabilization.status == MVMC_KRYLOV_STATUS_OK &&
                  stabilization.version ==
                      MVMC_POWER_LANCZOS_STABILIZATION_RESULT_VERSION &&
                  stabilization.block_count == execution.block_count &&
                  stabilization.chain_count == summary.chain_count &&
                  stabilization.coefficient_sample_count ==
                      summary.coefficient_sample_count &&
                  stabilization.final_sample_count ==
                      summary.final_sample_count &&
                  stabilization.coefficient_block_length == 128 &&
                  stabilization.final_block_length == 256 &&
                  isfinite(stabilization.coefficient_tau_int_max) &&
                  stabilization.coefficient_tau_int_max >= 1.0 &&
                  isfinite(stabilization.coefficient_effective_block_count) &&
                  stabilization.coefficient_effective_block_count > 0.0 &&
                  isfinite(stabilization.final_tau_int_max) &&
                  stabilization.final_tau_int_max >= 1.0 &&
                  isfinite(stabilization.final_effective_block_count) &&
                  stabilization.final_effective_block_count > 0.0 &&
                  isfinite(stabilization.energy) &&
                  isfinite(stabilization.variance) &&
                  isfinite(stabilization.variance_standard_error) &&
                  stabilization.maximum_absolute_numeric_bound >= 0.0 &&
                  (stabilization.failed_gates &
                   MVMC_POWER_LANCZOS_STABILIZATION_GATE_BLOCK_COUNT) != 0 &&
                  (stabilization.failed_gates &
                   MVMC_POWER_LANCZOS_STABILIZATION_GATE_NEGATIVE_VARIANCE) ==
                      0 &&
                  stabilization.decision !=
                      MVMC_POWER_LANCZOS_STABILIZATION_PASS,
              "lean stabilization evidence step %d status=%s decision=%s "
              "failed=0x%016" PRIx64,
              step, mvmc_krylov_status_string(stabilization.status),
              mvmc_power_lanczos_stabilization_decision_string(
                  stabilization.decision),
              stabilization.failed_gates);
        first_render_status =
            mvmc_power_lanczos_stabilization_output_render(
                session, &stabilization, &identity, first_output,
                sizeof(first_output), &first_size);
        second_render_status =
            mvmc_power_lanczos_stabilization_output_render(
                session, &stabilization, &identity, second_output,
                sizeof(second_output), &second_size);
        CHECK(first_render_status == MVMC_KRYLOV_STATUS_OK &&
                  second_render_status == MVMC_KRYLOV_STATUS_OK &&
                  first_size == second_size && first_size > 0 &&
                  first_size < sizeof(first_output) &&
                  memcmp(first_output, second_output, first_size + 1) == 0 &&
                  strstr(first_output, "\"blocks\"") != NULL &&
                  strstr(first_output, "\"schema_version\":1") != NULL &&
                  strstr(first_output, "\"log_basis_scale\"") != NULL &&
                  strstr(first_output, "\"resolved_eta\"") != NULL &&
                  strstr(first_output,
                         "\"coefficient_accepted_steps_total\"") != NULL &&
                  strstr(first_output,
                         "\"final_accepted_steps_total\"") != NULL &&
                  strstr(first_output,
                         "\"scale_pilot_sample_count_per_chain\":16") !=
                      NULL &&
                  strstr(first_output, "primitive") == NULL &&
                  strstr(first_output, "observable") == NULL &&
                  strstr(first_output, "bootstrap") == NULL,
              "deterministic compact output step %d status=%s/%s size=%zu",
              step, mvmc_krylov_status_string(first_render_status),
              mvmc_krylov_status_string(second_render_status), first_size);
        too_small[0] = 'x';
        CHECK(mvmc_power_lanczos_stabilization_output_render(
                  session, &stabilization, &identity, too_small,
                  sizeof(too_small), &too_small_size) ==
                      MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
                  too_small_size == 0 && too_small[0] == '\0',
              "compact output capacity failure is fail-closed step %d",
              step);
      } else {
        memset(&stabilization, 0xff, sizeof(stabilization));
        CHECK(mvmc_power_lanczos_stabilization_statistics_evaluate(
                  session, 0, &stabilization) ==
                      MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                  !stabilization.valid,
              "lean stabilization evidence remains root-owned step %d",
              step);
      }
    }
    mvmc_power_lanczos_production_session_destroy(session);
  }
  {
    MVMCPowerLanczosProductionSessionConfig config = session_config(1);
    MVMCPowerLanczosProductionExecution execution =
        execution_fixture(1, &model, chain_communicator);
    MVMCPowerLanczosProductionSession *session = NULL;
    MVMCPowerLanczosProductionSessionSummary summary;
    AmplitudeFailureContext failure;
    memset(&failure, 0, sizeof(failure));
    failure.failure_rank = world_size - 1;
    failure.failure_call = 2;
    config.split_size = split_size;
    config.configuration_bit_count = 6;
    execution.amplitude = failing_three_site_amplitude;
    execution.amplitude_context = &failure;
    execution.bounded_limits.max_row_transitions = 32;
    execution.coefficient_sample_count = 32;
    execution.final_sample_count = 32;
    execution.block_count = 4;
    execution.maximum_leave_one_projective_distance = 1.0;
    CHECK(create_test_session(
              &config, words, 1, &session) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_production_session_prepare_execution(
                  session, &execution) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_power_lanczos_production_session_run_coefficient(
                  session) == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
              mvmc_power_lanczos_production_session_summary(
                  session, &summary) ==
                  MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
              !summary.valid &&
              summary.state == MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FAILED,
          "rank-local amplitude failure converges globally");
    mvmc_power_lanczos_production_session_destroy(session);
  }
#ifdef _mpi_use
  if (chain_communicator != MPI_COMM_NULL) {
    MPI_Comm_free(&chain_communicator);
  }
#endif
}

static void test_mpi_domain_census(void) {
#ifdef _mpi_use
  MVMCPowerLanczosRngDomain coefficient;
  uint64_t local;
  uint64_t gathered[4] = {0, 0, 0, 0};
  int left;
  int right;
  CHECK(world_size <= 4 &&
            mvmc_power_lanczos_rng_derive(
                UINT64_C(0x503641434f4e4631), UINT64_C(7),
                MVMC_POWER_LANCZOS_RNG_STAGE_COEFFICIENT,
                (size_t)world_rank, (size_t)world_size, 2,
                &coefficient) == MVMC_KRYLOV_STATUS_OK,
        "MPI coefficient domain");
  local = coefficient.derived_seed;
  MPI_Allgather(&local, 1, MPI_UINT64_T, gathered, 1, MPI_UINT64_T,
                MPI_COMM_WORLD);
  for (left = 0; left < world_size; ++left) {
    for (right = left + 1; right < world_size; ++right) {
      CHECK(gathered[left] != gathered[right],
            "MPI rank seed collision %d/%d", left, right);
    }
  }
#endif
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#else
  (void)argc;
  (void)argv;
#endif
  test_sha256();
  test_rng();
  test_snapshot_lifecycle();
  test_snapshot_failures();
  test_session();
  test_terminal_snapshot_failures();
  test_execution_prepare();
  test_execution_full_lifecycle();
  test_mpi_domain_census();
#ifdef _mpi_use
  {
    int global_failures = 0;
    MPI_Allreduce(&failures, &global_failures, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    failures = global_failures;
  }
#endif
  if (failures == 0 && world_rank == 0) {
    puts("power-Lanczos production session unit: PASS");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return failures == 0 ? 0 : 1;
}
