#include "classic_krylov_model.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failure_count = 0;

#define CHECK(condition, ...)                                 \
  do {                                                        \
    if (!(condition)) {                                       \
      fprintf(stderr, "FAIL: ");                             \
      fprintf(stderr, __VA_ARGS__);                           \
      fprintf(stderr, " (line %d)\n", __LINE__);             \
      ++failure_count;                                        \
    }                                                         \
  } while (0)

static void rank_and_size(int *rank, int *size) {
#ifdef _mpi_use
  MPI_Comm_rank(MPI_COMM_WORLD, rank);
  MPI_Comm_size(MPI_COMM_WORLD, size);
#else
  *rank = 0;
  *size = 1;
#endif
}

static MVMCClassicPfaffianCommunicator world_communicator(void) {
#ifdef _mpi_use
  return MPI_COMM_WORLD;
#else
  return 0;
#endif
}

static MVMCClassicKrylovRawModel electronic_raw(
    const MVMCClassicKrylovTransfer *transfers, int transfer_count,
    const MVMCClassicKrylovSiteCoupling *intra, int intra_count,
    const MVMCClassicKrylovPairCoupling *inter, int inter_count,
    const MVMCClassicKrylovPairCoupling *hund, int hund_count) {
  MVMCClassicKrylovRawModel raw;
  memset(&raw, 0, sizeof(raw));
  raw.site_count = 2;
  raw.up_electron_count = 1;
  raw.down_electron_count = 1;
  raw.transfer_count = transfer_count;
  raw.transfers = transfers;
  raw.coulomb_intra_count = intra_count;
  raw.coulomb_intra = intra;
  raw.coulomb_inter_count = inter_count;
  raw.coulomb_inter = inter;
  raw.hund_count = hund_count;
  raw.hund = hund;
  return raw;
}

static MVMCKrylovStatus create_root_model(
    const MVMCClassicKrylovRawModel *raw,
    MVMCClassicKrylovModelWorkspace **workspace) {
  int rank, size;
  rank_and_size(&rank, &size);
  (void)size;
  return mvmc_classic_krylov_model_workspace_create_from_root(
      rank == 0 ? raw : NULL, world_communicator(), workspace);
}

typedef struct {
  double complex values[16];
} PureSpinOracle;

static MVMCKrylovStatus pure_spin_oracle_amplitude(
    const uint64_t *configuration_words, size_t word_count,
    void *context, MVMCKrylovAmplitudeResult *result) {
  PureSpinOracle *oracle = (PureSpinOracle *)context;
  uint64_t configuration;
  if (configuration_words == NULL || word_count != 1 || oracle == NULL ||
      result == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  configuration = configuration_words[0];
  if (configuration >= 16) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(result, 0, sizeof(*result));
  result->value = oracle->values[configuration];
  result->total_zero = result->value == 0.0;
  result->regular_component_count = result->total_zero ? 0 : 1;
  result->local_factorization_count = 1;
  result->global_factorization_count = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static int evaluate_first_order(
    const MVMCKrylovFockModel *model, uint64_t root,
    PureSpinOracle *oracle, double complex *value) {
  MVMCKrylovLimits limits;
  MVMCKrylovWorkspace *workspace;
  MVMCKrylovResult result;
  MVMCKrylovStatus status;
  memset(&limits, 0, sizeof(limits));
  limits.max_states = 16;
  limits.max_transitions = 256;
  limits.max_amplitude_evaluations = 16;
  limits.max_bytes = (size_t)1024 * 1024;
  limits.max_order = 1;
  workspace = mvmc_krylov_workspace_create(model->site_count, &limits,
                                            &status);
  if (workspace == NULL) return 0;
  status = mvmc_krylov_evaluate(
      workspace, model, &root, 1, pure_spin_oracle_amplitude, oracle,
      &result);
  mvmc_krylov_workspace_destroy(workspace);
  if (status != MVMC_KRYLOV_STATUS_OK || !result.valid ||
      result.evaluated_order != 1) {
    return 0;
  }
  *value = result.value[1];
  return 1;
}

static void test_electronic_mapping_and_permutation(void) {
  const MVMCClassicKrylovTransfer transfers_a[] = {
      {0, 0, 1, 0, 0.5 + 0.25 * I},
      {1, 0, 0, 0, 0.7 - 0.25 * I},
      {0, 0, 1, 0, 0.2},
      {0, 1, 0, 1, 0.3},
  };
  const MVMCClassicKrylovTransfer transfers_b[] = {
      {0, 1, 0, 1, 0.3},
      {0, 0, 1, 0, 0.2},
      {0, 0, 1, 0, 0.5 + 0.25 * I},
      {1, 0, 0, 0, 0.7 - 0.25 * I},
  };
  const MVMCClassicKrylovSiteCoupling intra_a[] = {
      {1, 0.7}, {1, -0.2}};
  const MVMCClassicKrylovSiteCoupling intra_b[] = {
      {1, -0.2}, {1, 0.7}};
  const MVMCClassicKrylovPairCoupling inter_a[] = {
      {1, 0, 0.4}, {0, 1, 0.1}};
  const MVMCClassicKrylovPairCoupling inter_b[] = {
      {0, 1, 0.1}, {1, 0, 0.4}};
  const MVMCClassicKrylovPairCoupling hund[] = {{0, 1, 0.6}};
  MVMCClassicKrylovRawModel raw_a = electronic_raw(
      transfers_a, 4, intra_a, 2, inter_a, 2, hund, 1);
  MVMCClassicKrylovRawModel raw_b = electronic_raw(
      transfers_b, 4, intra_b, 2, inter_b, 2, hund, 1);
  MVMCClassicKrylovModelWorkspace *workspace_a = NULL;
  MVMCClassicKrylovModelWorkspace *workspace_b = NULL;
  const MVMCKrylovFockModel *model_a;
  const MVMCKrylovFockModel *model_b;

  CHECK(create_root_model(&raw_a, &workspace_a) == MVMC_KRYLOV_STATUS_OK,
        "electronic model A create");
  CHECK(create_root_model(&raw_b, &workspace_b) == MVMC_KRYLOV_STATUS_OK,
        "electronic model B create");
  if (workspace_a == NULL || workspace_b == NULL) goto cleanup;
  model_a = mvmc_classic_krylov_model(workspace_a);
  model_b = mvmc_classic_krylov_model(workspace_b);
  CHECK(model_a != NULL && model_b != NULL,
        "model accessor returned NULL");
  if (model_a == NULL || model_b == NULL) goto cleanup;
  CHECK(model_a->term_count == 10 && model_a->operator_count == 34,
        "electronic family expansion count");
  CHECK(model_a->term_count == model_b->term_count &&
            model_a->operator_count == model_b->operator_count &&
            memcmp(model_a->terms, model_b->terms,
                   model_a->term_count * sizeof(*model_a->terms)) == 0 &&
            memcmp(model_a->operators, model_b->operators,
                   model_a->operator_count * sizeof(*model_a->operators)) ==
                0,
        "raw permutation changed canonical term list");
  CHECK(model_a->terms[0].source_kind ==
            MVMC_CLASSIC_KRYLOV_SOURCE_TRANSFER &&
            model_a->terms[0].coefficient == -(0.7 + 0.25 * I),
        "duplicate Transfer coefficient/sign mapping");
  CHECK(model_a->terms[3].source_kind ==
            MVMC_CLASSIC_KRYLOV_SOURCE_COULOMB_INTRA &&
            cabs(model_a->terms[3].coefficient - 0.5) < 2.0e-15,
        "CoulombIntra merge/mapping");
  CHECK(model_a->terms[8].source_kind ==
            MVMC_CLASSIC_KRYLOV_SOURCE_HUND &&
            model_a->terms[8].coefficient == -0.6,
        "Hund coefficient sign mapping");
  CHECK(mvmc_classic_krylov_model_workspace_bytes(workspace_a) >
            sizeof(*model_a),
        "model workspace byte accounting");

cleanup:
  mvmc_classic_krylov_model_workspace_destroy(workspace_b);
  mvmc_classic_krylov_model_workspace_destroy(workspace_a);
}

static void test_pure_spin_mapping(void) {
  const MVMCClassicKrylovPairCoupling inter[] = {{0, 1, -0.25}};
  const MVMCClassicKrylovPairCoupling hund[] = {{0, 1, -0.5}};
  const MVMCClassicKrylovPairCoupling exchange[] = {{1, 0, -0.5}};
  MVMCClassicKrylovRawModel raw;
  MVMCClassicKrylovModelWorkspace *workspace = NULL;
  const MVMCKrylovFockModel *model;
  PureSpinOracle oracle;
  double complex value = NAN + I * NAN;
  memset(&raw, 0, sizeof(raw));
  raw.site_count = 2;
  raw.up_electron_count = 1;
  raw.down_electron_count = 1;
  raw.pure_spin = 1;
  raw.coulomb_inter_count = 1;
  raw.coulomb_inter = inter;
  raw.hund_count = 1;
  raw.hund = hund;
  raw.exchange_count = 1;
  raw.exchange = exchange;
  CHECK(create_root_model(&raw, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "pure-spin model create");
  if (workspace == NULL) return;
  model = mvmc_classic_krylov_model(workspace);
  CHECK(model != NULL && model->pure_spin == 1 && model->term_count == 8 &&
            model->operator_count == 32,
        "pure-spin family expansion count");
  if (model != NULL) {
    CHECK(model->terms[0].source_kind ==
              MVMC_CLASSIC_KRYLOV_SOURCE_COULOMB_INTER &&
              model->terms[0].coefficient == -0.25 &&
              model->terms[4].source_kind ==
                  MVMC_CLASSIC_KRYLOV_SOURCE_HUND &&
              model->terms[4].coefficient == 0.5 &&
              model->terms[6].source_kind ==
                  MVMC_CLASSIC_KRYLOV_SOURCE_EXCHANGE &&
              model->terms[6].coefficient == -0.5,
          "StdFace CoulombInter/Hund/Exchange coefficient mapping");
    CHECK(model->operators[24].orbital == 0 &&
              model->operators[25].orbital == 1 &&
              model->operators[26].orbital == 3 &&
              model->operators[27].orbital == 2,
          "Exchange forward orientation");
    /* The up-orbitals-before-down Fock convention reverses the coordinate
       phase of the two physical spin basis states; the eigenvalue set is
       nevertheless the Heisenberg {-3/4, +1/4} spectrum. */
    memset(&oracle, 0, sizeof(oracle));
    oracle.values[UINT64_C(6)] = 1.0;
    oracle.values[UINT64_C(9)] = 1.0;
    CHECK(evaluate_first_order(model, UINT64_C(6), &oracle, &value) &&
              cabs(value + 0.75) < 2.0e-15,
          "StdFace Heisenberg symmetric-Fock eigenvalue actual=(%.17g,%.17g)",
          creal(value), cimag(value));
    oracle.values[UINT64_C(9)] = -1.0;
    CHECK(evaluate_first_order(model, UINT64_C(6), &oracle, &value) &&
              cabs(value - 0.25) < 2.0e-15,
          "StdFace Heisenberg antisymmetric-Fock eigenvalue actual=(%.17g,%.17g)",
          creal(value), cimag(value));
  }
  mvmc_classic_krylov_model_workspace_destroy(workspace);
}

static void expect_status(MVMCClassicKrylovRawModel *raw,
                          MVMCKrylovStatus expected, const char *message) {
  MVMCClassicKrylovModelWorkspace *workspace = NULL;
  CHECK(create_root_model(raw, &workspace) == expected && workspace == NULL,
        "%s", message);
  mvmc_classic_krylov_model_workspace_destroy(workspace);
}

static void test_failure_contract(void) {
  MVMCClassicKrylovTransfer transfers[2] = {
      {0, 0, 1, 0, 1.0}, {1, 0, 0, 0, 1.0}};
  MVMCClassicKrylovRawModel raw =
      electronic_raw(transfers, 2, NULL, 0, NULL, 0, NULL, 0);

  raw.transfer_count = -1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "negative raw count rejected");
  raw.transfer_count = 2;
  raw.coulomb_intra_count = -1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "negative CoulombIntra count rejected");
  raw.coulomb_intra_count = 0;
  raw.coulomb_inter_count = -1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "negative CoulombInter count rejected");
  raw.coulomb_inter_count = 0;
  raw.hund_count = -1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "negative Hund count rejected");
  raw.hund_count = 0;
  raw.exchange_count = -1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "negative Exchange count rejected");
  raw.exchange_count = 0;
  raw.pair_hopping_count = -1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "negative PairHopping count rejected");
  raw.pair_hopping_count = 0;
  raw.inter_all_count = -1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "negative InterAll count rejected");
  raw.inter_all_count = 0;
  raw.nbody_inter_all_count = -1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "negative NBodyInterAll count rejected");
  raw.nbody_inter_all_count = 0;
  raw.hund_count = 1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "positive count with NULL binding rejected");
  raw.hund_count = 0;
  transfers[0].output_site = 2;
  expect_status(&raw, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "out-of-range Transfer site rejected");
  transfers[0].output_site = 0;
  transfers[0].input_spin = 1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL,
                "spin-changing Transfer rejected as unsupported");
  transfers[0].input_spin = 0;
  transfers[0].coefficient = NAN;
  expect_status(&raw, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "nonfinite Transfer rejected");
  transfers[0].coefficient = 1.0;
  raw.transfer_count = 1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL,
                "missing Transfer adjoint rejected");
  raw.transfer_count = 2;
  transfers[1].coefficient = 0.5;
  expect_status(&raw, MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL,
                "mismatched Transfer adjoint rejected");
  transfers[1].coefficient = nextafter(1.0, 2.0);
  expect_status(&raw, MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL,
                "one-ulp Transfer adjoint mismatch rejected");
  transfers[1].coefficient = 1.0;
  transfers[0].output_site = transfers[0].input_site = 0;
  transfers[0].coefficient = 1.0 + I;
  raw.transfer_count = 1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL,
                "complex diagonal Transfer rejected");

  transfers[0] = (MVMCClassicKrylovTransfer){0, 0, 1, 0, 1.0};
  transfers[1] = (MVMCClassicKrylovTransfer){1, 0, 0, 0, 1.0};
  raw.transfer_count = 2;
  raw.pair_hopping_count = 1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL,
                "PairHopping rejected");
  raw.pair_hopping_count = 0;
  raw.inter_all_count = 1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL,
                "InterAll rejected");
  raw.inter_all_count = 0;
  raw.nbody_inter_all_count = 1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL,
                "NBodyInterAll rejected");
  raw.nbody_inter_all_count = 0;

  {
    MVMCClassicKrylovPairCoupling exchange = {0, 1, 1.0};
    raw.exchange_count = 1;
    raw.exchange = &exchange;
    expect_status(&raw, MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL,
                  "electronic Exchange rejected");
    exchange.first_site = 2;
    expect_status(&raw, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                  "malformed electronic Exchange is invalid before unsupported");
    exchange.first_site = 0;
    exchange.coefficient = NAN;
    expect_status(&raw, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                  "nonfinite electronic Exchange is invalid before unsupported");
    raw.exchange_count = 0;
    raw.exchange = NULL;
  }
  transfers[1].input_site = 2;
  transfers[0].input_spin = 1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "malformed Transfer is invalid before spin-change unsupported");
  transfers[0].input_spin = 0;
  transfers[1].input_site = 0;
  raw.pure_spin = 1;
  expect_status(&raw, MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL,
                "pure-spin Transfer rejected");
}

static void test_cancelled_duplicate(void) {
  const MVMCClassicKrylovTransfer transfers[] = {
      {0, 0, 1, 0, 1.0}, {0, 0, 1, 0, -1.0},
      {1, 0, 0, 0, 1.0}, {1, 0, 0, 0, -1.0}};
  MVMCClassicKrylovRawModel raw =
      electronic_raw(transfers, 4, NULL, 0, NULL, 0, NULL, 0);
  MVMCClassicKrylovModelWorkspace *workspace = NULL;
  const MVMCKrylovFockModel *model;
  CHECK(create_root_model(&raw, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "cancelled duplicate model create");
  model = mvmc_classic_krylov_model(workspace);
  CHECK(model != NULL && model->term_count == 0 &&
            model->operator_count == 0,
        "exactly cancelled raw terms were not removed");
  mvmc_classic_krylov_model_workspace_destroy(workspace);
}

static void test_split_communicator(void) {
#ifdef _mpi_use
  const MVMCClassicKrylovTransfer transfers[] = {
      {0, 0, 1, 0, 1.0}, {1, 0, 0, 0, 1.0}};
  MVMCClassicKrylovRawModel raw =
      electronic_raw(transfers, 2, NULL, 0, NULL, 0, NULL, 0);
  MVMCClassicKrylovModelWorkspace *workspace = NULL;
  MPI_Comm split = MPI_COMM_NULL;
  int world_rank, world_size, split_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size < 2) return;
  MPI_Comm_split(MPI_COMM_WORLD, world_rank % 2, world_rank, &split);
  MPI_Comm_rank(split, &split_rank);
  CHECK(mvmc_classic_krylov_model_workspace_create_from_root(
            split_rank == 0 ? &raw : NULL, split, &workspace) ==
            MVMC_KRYLOV_STATUS_OK,
        "split communicator root-only model create");
  mvmc_classic_krylov_model_workspace_destroy(workspace);
  MPI_Comm_free(&split);
#endif
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  test_electronic_mapping_and_permutation();
  test_pure_spin_mapping();
  test_failure_contract();
  test_cancelled_duplicate();
  test_split_communicator();

  if (failure_count != 0) {
    fprintf(stderr, "%d classic Krylov model checks failed\n", failure_count);
#ifdef _mpi_use
    MPI_Finalize();
#endif
    return 1;
  }
  puts("classic Krylov model unit checks passed");
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return 0;
}
