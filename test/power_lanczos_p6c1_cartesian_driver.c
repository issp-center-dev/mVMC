#include "classic_krylov_amplitude.h"
#include "classic_krylov_model.h"
#include "bounded_krylov_engine.h"
#include "krylov_fock_reference.h"
#include "power_lanczos_gevp.h"

#include <complex.h>
#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

enum {
  SITE_COUNT = 4,
  ORBITAL_COUNT = 8,
  MAX_QP = 8,
  MAX_PROJECTION = 32,
  MAX_STATES = 36,
  MAX_VECTORS = 4,
  MAX_UPPER = 6
};

typedef struct {
  MVMCClassicKrylovRawModel raw;
  MVMCClassicKrylovTransfer transfers[16];
  MVMCClassicKrylovSiteCoupling intra[SITE_COUNT];
  MVMCClassicKrylovPairCoupling inter[SITE_COUNT];
  MVMCClassicKrylovPairCoupling hund[SITE_COUNT];
  MVMCClassicKrylovPairCoupling exchange[SITE_COUNT];
} FixtureModel;

typedef struct {
  char fixture_id[16];
  char logical_sha[65];
  char rendered_sha[65];
  int model;
  int order;
  int arithmetic;
  int nqp;
  int profile;
  int support;
  int nproj;
  double complex pair[SITE_COUNT * SITE_COUNT];
  double complex weights[MAX_QP];
  int mapping[MAX_QP][SITE_COUNT];
  double complex projection[MAX_PROJECTION];
  uint64_t proof_x;
  uint64_t proof_y;
} FixtureInput;

typedef struct {
  uint64_t source;
} DeltaAmplitude;

static const double CUTOFFS[4] = {
    0x1p-48, 0x1p-40, 0x1p-32, 0x1p-24};
static const char *const CUTOFF_IDS[4] = {"S48", "S40", "S32", "S24"};
static const char FROZEN_MANIFEST_SHA[65] =
    "78e9d30c6603d1529b46da2e2d763fbe52b74be76528f8a6709e7e5f0882d34d";
static const char FROZEN_REGISTRY_SHA[65] =
    "ae1e4eb7a1e0766c5c8727c4ae66d0e056f884dd0017cd02ff99d8176d7488fb";
static const char FROZEN_MODULE_SHA[65] =
    "2ad8053037e9c974b90d374ca69324ac317aaad3f5641a05265c2cf257345481";
static const char NUMERICAL_POLICY_SHA[65] =
    "550e488543a88ba50d91294b66262924a1845a6a91fc5b53351d45acfd256389";

#ifdef MVMC_P6C1_PREFLIGHT
enum { EXPECTED_CASE_COUNT = 9 };
static const char INPUT_HEADER[] = "P6C1_CARTESIAN_PREFLIGHT_INPUT_V1";
static const char OUTPUT_HEADER[] = "P6C1_CARTESIAN_PREFLIGHT_OUTPUT_V1";
static const char *const PREFLIGHT_FIXTURE_IDS[EXPECTED_CASE_COUNT] = {
    "F0001", "F1765", "F0254", "F1095", "F0423",
    "F0670", "F2014", "F1007", "F1679"};
#else
enum { EXPECTED_CASE_COUNT = 2016 };
static const char INPUT_HEADER[] = "P6C1_CARTESIAN_INPUT_V1";
static const char OUTPUT_HEADER[] = "P6C1_CARTESIAN_OUTPUT_V1";
#endif

static int world_rank(void) {
#ifdef _mpi_use
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
#else
  return 0;
#endif
}

static int world_size(void) {
#ifdef _mpi_use
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
#else
  return 1;
#endif
}

static MVMCClassicPfaffianCommunicator world_communicator(void) {
#ifdef _mpi_use
  return MPI_COMM_WORLD;
#else
  return 0;
#endif
}

static int collective_failure(int local_failure) {
#ifdef _mpi_use
  int global_failure = 0;
  if (MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    return 1;
  }
  return global_failure;
#else
  return local_failure;
#endif
}

static int bit_count(unsigned int value) {
  int count = 0;
  while (value != 0U) {
    count += (int)(value & 1U);
    value >>= 1;
  }
  return count;
}

static int canonical_sha256(const char *value) {
  int index;
  if (value == NULL || strlen(value) != 64) return 0;
  for (index = 0; index < 64; ++index) {
    if (!((value[index] >= '0' && value[index] <= '9') ||
          (value[index] >= 'a' && value[index] <= 'f'))) {
      return 0;
    }
  }
  return 1;
}

static int valid_configuration(uint64_t bits, int pure_spin) {
  const unsigned int up = (unsigned int)(bits & 0x0fU);
  const unsigned int down = (unsigned int)((bits >> SITE_COUNT) & 0x0fU);
  return (bits & ~UINT64_C(0xff)) == 0 && bit_count(up) == 2 &&
         bit_count(down) == 2 &&
         (!pure_spin || down == ((~up) & 0x0fU));
}

static int configuration_list(int pure_spin, uint64_t *states) {
  int count = 0;
  uint64_t bits;
  for (bits = 0; bits < UINT64_C(256); ++bits) {
    if (valid_configuration(bits, pure_spin)) states[count++] = bits;
  }
  return count;
}

static void print_real(double value) {
  if (value == 0.0) value = 0.0;
  printf(" %a", value);
}

static void print_complex(double complex value) {
  print_real(creal(value));
  print_real(cimag(value));
}

static void print_error_bound(double value) {
  if (isinf(value) && value < 0.0) {
    printf(" NINF");
  } else {
    print_real(value);
  }
}

static double complex inner_product(const double complex *left,
                                    const double complex *right,
                                    int count) {
  double complex value = 0.0;
  int index;
  for (index = 0; index < count; ++index) {
    value += conj(left[index]) * right[index];
  }
  return value;
}

static int observable_consistent(double left, double right) {
  const double difference = fabs(left - right);
  return difference <= 1.0e-10 ||
         difference <= 1.0e-11 * fmax(1.0, fabs(left));
}

static size_t upper_index(int dimension, int row, int column) {
  size_t index = 0;
  int current;
  for (current = 0; current < row; ++current) {
    index += (size_t)(dimension - current);
  }
  return index + (size_t)(column - row);
}

static MVMCKrylovBoundedLimits bounded_limits(int maximum_order) {
  MVMCKrylovBoundedLimits limits;
  limits.amplitude_policy_hash = UINT64_C(0x50364331414d5001);
  limits.cache_bytes = 1024 * 1024;
  limits.max_row_transitions = 4096;
  limits.max_workspace_bytes = 64 * 1024 * 1024;
  limits.max_node_expansions = UINT64_C(1000000);
  limits.max_terminal_amplitude_calls = UINT64_C(1000000);
  limits.max_total_row_transitions = UINT64_C(10000000);
  limits.max_order = maximum_order;
  return limits;
}

static int scaled_to_raw(const MVMCScaledComplex *value,
                         double complex *raw) {
  MVMCScaledComplexExportStatus status;
  if (value == NULL || raw == NULL) return 0;
  status = mvmc_scaled_complex_export_common_scale(value, 0.0, raw);
  return status == MVMC_SCALED_EXPORT_OK ||
         status == MVMC_SCALED_EXPORT_EXACT_ZERO ||
         status == MVMC_SCALED_EXPORT_NUMERIC_ZERO;
}

static void initialize_electronic_model(FixtureModel *fixture,
                                        int complex_case) {
  static const int bonds[SITE_COUNT][2] = {
      {0, 1}, {1, 2}, {2, 3}, {3, 0}};
  static const double imaginary[SITE_COUNT] = {
      1.0 / 8.0, 1.0 / 16.0, -1.0 / 8.0, -1.0 / 16.0};
  static const double intra[SITE_COUNT] = {4.0, 17.0 / 4.0, 9.0 / 2.0,
                                           19.0 / 4.0};
  static const double inter[SITE_COUNT] = {1.0, 9.0 / 8.0, 5.0 / 4.0,
                                           11.0 / 8.0};
  static const double hund[SITE_COUNT] = {1.0 / 4.0, 5.0 / 16.0,
                                          3.0 / 8.0, 7.0 / 16.0};
  int transfer_count = 0;
  int bond;
  int spin;
  memset(fixture, 0, sizeof(*fixture));
  for (bond = 0; bond < SITE_COUNT; ++bond) {
    const double complex value =
        -1.0 + (complex_case ? imaginary[bond] : 0.0) * I;
    for (spin = 0; spin < 2; ++spin) {
      fixture->transfers[transfer_count++] =
          (MVMCClassicKrylovTransfer){
              bonds[bond][0], spin, bonds[bond][1], spin, value};
      fixture->transfers[transfer_count++] =
          (MVMCClassicKrylovTransfer){
              bonds[bond][1], spin, bonds[bond][0], spin, conj(value)};
    }
    fixture->intra[bond] =
        (MVMCClassicKrylovSiteCoupling){bond, intra[bond]};
    fixture->inter[bond] =
        (MVMCClassicKrylovPairCoupling){
            bonds[bond][0], bonds[bond][1], inter[bond]};
    fixture->hund[bond] =
        (MVMCClassicKrylovPairCoupling){
            bonds[bond][0], bonds[bond][1], hund[bond]};
  }
  fixture->raw.site_count = SITE_COUNT;
  fixture->raw.up_electron_count = 2;
  fixture->raw.down_electron_count = 2;
  fixture->raw.transfer_count = transfer_count;
  fixture->raw.transfers = fixture->transfers;
  fixture->raw.coulomb_intra_count = SITE_COUNT;
  fixture->raw.coulomb_intra = fixture->intra;
  fixture->raw.coulomb_inter_count = SITE_COUNT;
  fixture->raw.coulomb_inter = fixture->inter;
  fixture->raw.hund_count = SITE_COUNT;
  fixture->raw.hund = fixture->hund;
}

static void initialize_spin_model(FixtureModel *fixture) {
  static const int bonds[SITE_COUNT][2] = {
      {0, 1}, {1, 2}, {2, 3}, {3, 0}};
  static const double jxy[SITE_COUNT] = {1.0, 9.0 / 8.0, 5.0 / 4.0,
                                         11.0 / 8.0};
  static const double jz[SITE_COUNT] = {3.0 / 4.0, 13.0 / 16.0,
                                        7.0 / 8.0, 15.0 / 16.0};
  int bond;
  memset(fixture, 0, sizeof(*fixture));
  for (bond = 0; bond < SITE_COUNT; ++bond) {
    fixture->hund[bond] =
        (MVMCClassicKrylovPairCoupling){
            bonds[bond][0], bonds[bond][1], -0.5 * jz[bond]};
    fixture->exchange[bond] =
        (MVMCClassicKrylovPairCoupling){
            bonds[bond][0], bonds[bond][1], 0.5 * jxy[bond]};
  }
  fixture->raw.site_count = SITE_COUNT;
  fixture->raw.up_electron_count = 2;
  fixture->raw.down_electron_count = 2;
  fixture->raw.pure_spin = 1;
  fixture->raw.hund_count = SITE_COUNT;
  fixture->raw.hund = fixture->hund;
  fixture->raw.exchange_count = SITE_COUNT;
  fixture->raw.exchange = fixture->exchange;
}

static int read_fixture(FILE *input, int expected_ordinal,
                        FixtureInput *fixture) {
  char token[64];
  char expected_id[16];
  int pair_index;
  int q;
  int projection_index;
  unsigned long long proof_x;
  unsigned long long proof_y;
  memset(fixture, 0, sizeof(*fixture));
#ifdef MVMC_P6C1_PREFLIGHT
  if (expected_ordinal < 1 || expected_ordinal > EXPECTED_CASE_COUNT) {
    return 0;
  }
  snprintf(expected_id, sizeof(expected_id), "%s",
           PREFLIGHT_FIXTURE_IDS[expected_ordinal - 1]);
#else
  snprintf(expected_id, sizeof(expected_id), "F%04d", expected_ordinal);
#endif
  if (fscanf(input, "%63s %15s %d %d %d %d %d %d %64s %64s %d",
             token, fixture->fixture_id, &fixture->model, &fixture->order,
             &fixture->arithmetic, &fixture->nqp, &fixture->profile,
             &fixture->support, fixture->logical_sha,
             fixture->rendered_sha, &fixture->nproj) != 11 ||
      strcmp(token, "CASE") != 0 ||
      strcmp(fixture->fixture_id, expected_id) != 0 ||
      fixture->model < 0 || fixture->model > 1 ||
      fixture->order < 1 || fixture->order > 2 ||
      fixture->arithmetic < 0 || fixture->arithmetic > 1 ||
      (fixture->nqp != 1 && fixture->nqp != 4 && fixture->nqp != 8) ||
      fixture->profile < 0 || fixture->profile > 6 ||
      fixture->support < 0 || fixture->support > 2 ||
      fixture->nproj < 0 || fixture->nproj > MAX_PROJECTION ||
      !canonical_sha256(fixture->logical_sha) ||
      !canonical_sha256(fixture->rendered_sha)) {
    return 0;
  }
  if (fscanf(input, "%63s", token) != 1 || strcmp(token, "PAIR") != 0) {
    return 0;
  }
  for (pair_index = 0; pair_index < SITE_COUNT * SITE_COUNT;
       ++pair_index) {
    double real;
    double imaginary;
    if (fscanf(input, "%la %la", &real, &imaginary) != 2) return 0;
    fixture->pair[pair_index] = real + I * imaginary;
  }
  if (fscanf(input, "%63s", token) != 1 || strcmp(token, "QP") != 0) {
    return 0;
  }
  for (q = 0; q < fixture->nqp; ++q) {
    int input_q;
    double real;
    double imaginary;
    if (fscanf(input, "%d %la %la %d %d %d %d", &input_q, &real,
               &imaginary, &fixture->mapping[q][0],
               &fixture->mapping[q][1], &fixture->mapping[q][2],
               &fixture->mapping[q][3]) != 7 || input_q != q) {
      return 0;
    }
    fixture->weights[q] = real + I * imaginary;
    for (int site = 0; site < SITE_COUNT; ++site) {
      if (fixture->mapping[q][site] < 0 ||
          fixture->mapping[q][site] >= SITE_COUNT) {
        return 0;
      }
    }
  }
  if (fscanf(input, "%63s", token) != 1 || strcmp(token, "PROJ") != 0) {
    return 0;
  }
  for (projection_index = 0; projection_index < fixture->nproj;
       ++projection_index) {
    double value;
    if (fscanf(input, "%la", &value) != 1) return 0;
    fixture->projection[projection_index] = value;
  }
  if (fscanf(input, "%63s %llu %llu", token, &proof_x, &proof_y) != 3 ||
      strcmp(token, "PROOF") != 0) {
    return 0;
  }
  fixture->proof_x = (uint64_t)proof_x;
  fixture->proof_y = (uint64_t)proof_y;
  if (fscanf(input, "%63s", token) != 1 || strcmp(token, "END") != 0) {
    return 0;
  }
  if (fixture->arithmetic == 0) {
    for (pair_index = 0; pair_index < SITE_COUNT * SITE_COUNT;
         ++pair_index) {
      if (cimag(fixture->pair[pair_index]) != 0.0) return 0;
    }
    for (q = 0; q < fixture->nqp; ++q) {
      if (cimag(fixture->weights[q]) != 0.0) return 0;
    }
  }
  return 1;
}

static int profile_projection_count(int profile) {
  static const int counts[7] = {0, 4, 6, 6, 6, 10, 32};
  return profile >= 0 && profile < 7 ? counts[profile] : -1;
}

static MVMCClassicKrylovAmplitudeLayout amplitude_layout(
    const FixtureInput *fixture) {
  static const int gutzwiller[SITE_COUNT] = {0, 1, 2, 3};
  static const int jastrow_flat[SITE_COUNT * SITE_COUNT] = {
      -1, 0, 1, 2,
       0,-1, 3, 4,
       1, 3,-1, 5,
       2, 4, 5,-1};
  static const int *const jastrow_rows[SITE_COUNT] = {
      jastrow_flat, jastrow_flat + 4, jastrow_flat + 8,
      jastrow_flat + 12};
  static const int dh2_flat[2 * SITE_COUNT] = {
      3, 1, 0, 2, 1, 3, 2, 0};
  static const int *const dh2_rows[1] = {dh2_flat};
  static const int dh4_flat[4 * SITE_COUNT] = {
      2, 3, 1, 2,
      3, 0, 2, 3,
      0, 1, 3, 0,
      1, 2, 0, 1};
  static const int *const dh4_rows[1] = {dh4_flat};
  MVMCClassicKrylovAmplitudeLayout layout;
  const int rank = world_rank();
  const int size = world_size();
  memset(&layout, 0, sizeof(layout));
  layout.site_count = SITE_COUNT;
  layout.up_electron_count = 2;
  layout.down_electron_count = 2;
  layout.pure_spin = fixture->model;
  layout.qp_total = fixture->nqp;
  layout.qp_start = fixture->nqp * rank / size;
  layout.qp_end = fixture->nqp * (rank + 1) / size;
  layout.scaled_pivot_tolerance = 0.0;
  layout.nproj = fixture->nproj;
  if (fixture->profile == 1 || fixture->profile == 6) {
    layout.ngutzwiller_idx = 4;
    layout.gutzwiller_idx = gutzwiller;
  }
  if (fixture->profile == 2 || fixture->profile == 6) {
    layout.njastrow_idx = 6;
    layout.jastrow_idx = jastrow_rows;
  }
  if (fixture->profile == 3 || fixture->profile == 6) {
    layout.nspin_jastrow_idx = 6;
    layout.spin_jastrow_idx = jastrow_rows;
  }
  if (fixture->profile == 4 || fixture->profile == 6) {
    layout.ndoublon_holon_2site_idx = 1;
    layout.doublon_holon_2site_idx = dh2_rows;
  }
  if (fixture->profile == 5 || fixture->profile == 6) {
    layout.ndoublon_holon_4site_idx = 1;
    layout.doublon_holon_4site_idx = dh4_rows;
  }
  layout.projection_parameters = fixture->projection;
  layout.communicator = world_communicator();
  return layout;
}

static void build_slater(const FixtureInput *fixture,
                         double *real_slater,
                         double complex *complex_slater) {
  int q;
  memset(real_slater, 0,
         MAX_QP * ORBITAL_COUNT * ORBITAL_COUNT * sizeof(*real_slater));
  memset(complex_slater, 0,
         MAX_QP * ORBITAL_COUNT * ORBITAL_COUNT *
             sizeof(*complex_slater));
  for (q = 0; q < fixture->nqp; ++q) {
    int up;
    int down;
    for (up = 0; up < SITE_COUNT; ++up) {
      for (down = 0; down < SITE_COUNT; ++down) {
        const int transformed_up = fixture->mapping[q][up];
        const int transformed_down = fixture->mapping[q][down];
        const int down_orbital = SITE_COUNT + down;
        const double complex value =
            fixture->pair[transformed_up * SITE_COUNT + transformed_down];
        const size_t offset =
            (size_t)q * ORBITAL_COUNT * ORBITAL_COUNT;
        complex_slater[offset +
                       (size_t)up * ORBITAL_COUNT + down_orbital] = value;
        complex_slater[offset +
                       (size_t)down_orbital * ORBITAL_COUNT + up] = -value;
        real_slater[offset +
                    (size_t)up * ORBITAL_COUNT + down_orbital] =
            creal(value);
        real_slater[offset +
                    (size_t)down_orbital * ORBITAL_COUNT + up] =
            -creal(value);
      }
    }
  }
}

static MVMCKrylovStatus delta_scaled_amplitude(
    const uint64_t *configuration_words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
  const DeltaAmplitude *delta = (const DeltaAmplitude *)context;
  if (configuration_words == NULL || word_count != 1 || delta == NULL ||
      result == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(result, 0, sizeof(*result));
  if (configuration_words[0] == delta->source) {
    if (mvmc_scaled_complex_from_raw(1.0, &result->value) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    }
    result->well_pivoted_component_count = 1;
  } else {
    if (mvmc_scaled_complex_make_exact_zero(&result->value) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    }
    result->exact_zero_component_count = 1;
  }
  result->local_factorization_count = 1;
  result->global_factorization_count = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static int emit_full_hamiltonian(int model_id, int arithmetic,
                                 const MVMCKrylovFockModel *model) {
  uint64_t states[MAX_STATES];
  const int state_count = configuration_list(model_id, states);
  MVMCKrylovBoundedLimits limits = bounded_limits(1);
  MVMCKrylovStatus status;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  int source;
  int target;
  status = mvmc_bounded_krylov_plan_create(model, &limits, &plan);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_workspace_create(plan, &workspace);
  }
  if (workspace == NULL || status != MVMC_KRYLOV_STATUS_OK) {
    mvmc_bounded_krylov_plan_destroy(plan);
    return 0;
  }
  if (world_rank() == 0) {
    printf("MODEL_H %d %d %d %zu %zu\n", model_id, arithmetic,
           state_count, model->term_count, model->operator_count);
  }
  for (target = 0; target < state_count; ++target) {
    for (source = 0; source < state_count; ++source) {
      DeltaAmplitude delta;
      MVMCKrylovBoundedResult result;
      double complex value;
      delta.source = states[source];
      status = mvmc_bounded_krylov_evaluate(
          workspace, &states[target], 1, delta_scaled_amplitude, &delta,
          &result);
      if (status != MVMC_KRYLOV_STATUS_OK || !result.valid ||
          result.evaluated_order != 1 ||
          !scaled_to_raw(&result.value[1], &value)) {
        mvmc_bounded_krylov_workspace_destroy(workspace);
        mvmc_bounded_krylov_plan_destroy(plan);
        return 0;
      }
      if (world_rank() == 0) {
        printf("HROW %" PRIu64 " %" PRIu64, states[target],
               states[source]);
        print_complex(value);
        putchar('\n');
      }
    }
  }
  if (world_rank() == 0) printf("END_MODEL_H %d %d\n", model_id, arithmetic);
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
  return 1;
}

static void print_pack(const char *name, const double complex *values,
                       size_t count) {
  size_t index;
  printf("PACK %s %zu", name, count);
  for (index = 0; index < count; ++index) print_complex(values[index]);
  putchar('\n');
}

static int evaluate_fixture(const FixtureInput *fixture,
                            int *model_h_emitted) {
  FixtureModel fixture_model;
  MVMCClassicKrylovModelWorkspace *model_workspace = NULL;
  MVMCClassicKrylovRealAmplitudeWorkspace *real_workspace = NULL;
  MVMCClassicKrylovComplexAmplitudeWorkspace *complex_workspace = NULL;
  MVMCKrylovBoundedPlan *krylov_plan = NULL;
  MVMCKrylovBoundedWorkspace *krylov_workspace = NULL;
  const MVMCKrylovFockModel *model;
  MVMCKrylovScaledAmplitudeCallback amplitude;
  void *amplitude_context = NULL;
  MVMCClassicKrylovAmplitudeLayout layout = amplitude_layout(fixture);
  MVMCKrylovBoundedLimits limits = bounded_limits(fixture->order + 1);
  MVMCKrylovStatus status;
  double real_slater[MAX_QP * ORBITAL_COUNT * ORBITAL_COUNT];
  double complex complex_slater[
      MAX_QP * ORBITAL_COUNT * ORBITAL_COUNT];
  uint64_t states[MAX_STATES];
  double complex vectors[MAX_VECTORS][MAX_STATES] = {{0.0}};
  double complex overlap[MAX_UPPER] = {0.0};
  double complex forward[MAX_UPPER] = {0.0};
  double complex reverse[MAX_UPPER] = {0.0};
  double complex squared[MAX_UPPER] = {0.0};
  const int state_count = configuration_list(fixture->model, states);
  const int dimension = fixture->order + 1;
  const size_t upper_count =
      (size_t)dimension * (size_t)(dimension + 1) / 2U;
  const int rank = world_rank();
  const int size = world_size();
  uint64_t generation_hash = 0;
  uint64_t plan_hash = 0;
  int state_index;
  int row;
  int column;
  int cutoff_index;
  int all_amplitudes_zero = 1;
  int all_vectors_zero = 1;
  int setup_failure = 0;
  int return_code = 0;

  if (fixture->nproj != profile_projection_count(fixture->profile)) {
    if (rank == 0) fprintf(stderr, "%s projection count mismatch\n",
                           fixture->fixture_id);
    return 0;
  }
  if (fixture->model == 0) {
    initialize_electronic_model(&fixture_model, fixture->arithmetic);
  } else {
    initialize_spin_model(&fixture_model);
  }
  status = mvmc_classic_krylov_model_workspace_create_from_root(
      rank == 0 ? &fixture_model.raw : NULL, world_communicator(),
      &model_workspace);
  if (status != MVMC_KRYLOV_STATUS_OK) goto cleanup;
  model = mvmc_classic_krylov_model(model_workspace);
  {
    const int key = fixture->model * 2 + fixture->arithmetic;
    if (!model_h_emitted[key]) {
      setup_failure = !emit_full_hamiltonian(
          fixture->model, fixture->arithmetic, model);
      if (collective_failure(setup_failure)) goto cleanup;
      model_h_emitted[key] = 1;
    }
  }
  build_slater(fixture, real_slater, complex_slater);
  if (fixture->arithmetic == 0) {
    status = mvmc_classic_krylov_real_amplitude_workspace_create(
        &layout, real_slater, fixture->weights, &real_workspace);
    amplitude = mvmc_classic_krylov_real_scaled_amplitude;
    amplitude_context = real_workspace;
    generation_hash =
        mvmc_classic_krylov_real_amplitude_generation_hash(real_workspace);
  } else {
    status = mvmc_classic_krylov_complex_amplitude_workspace_create(
        &layout, complex_slater, fixture->weights, &complex_workspace);
    amplitude = mvmc_classic_krylov_complex_scaled_amplitude;
    amplitude_context = complex_workspace;
    generation_hash =
        mvmc_classic_krylov_complex_amplitude_generation_hash(
            complex_workspace);
  }
  setup_failure = status != MVMC_KRYLOV_STATUS_OK ||
                  amplitude_context == NULL || generation_hash == 0;
  if (collective_failure(setup_failure)) goto cleanup;
  status = mvmc_bounded_krylov_plan_create(model, &limits, &krylov_plan);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    plan_hash = mvmc_bounded_krylov_plan_hash(krylov_plan);
    status = mvmc_bounded_krylov_workspace_create(
        krylov_plan, &krylov_workspace);
  }
  setup_failure = krylov_workspace == NULL ||
                  status != MVMC_KRYLOV_STATUS_OK || plan_hash == 0;
  if (collective_failure(setup_failure)) goto cleanup;
  if (rank == 0) {
    printf("BEGIN %s %d %d %d %d %d %d %s %s %" PRIu64
           " %" PRIu64 " %d %zu %zu\n",
           fixture->fixture_id, fixture->model, fixture->order,
           fixture->arithmetic, fixture->nqp, fixture->profile,
           fixture->support, fixture->logical_sha, fixture->rendered_sha,
           generation_hash, plan_hash, state_count, model->term_count,
           model->operator_count);
    printf("PARTITIONS %d", size);
    for (int partition_rank = 0; partition_rank < size; ++partition_rank) {
      printf(" %d %d %d", partition_rank,
             fixture->nqp * partition_rank / size,
             fixture->nqp * (partition_rank + 1) / size);
    }
    putchar('\n');
    printf("PIVOT 0x0p+0 0x1p-44\n");
    printf("PROOF %" PRIu64 " %" PRIu64 "\n", fixture->proof_x,
           fixture->proof_y);
  }
  for (state_index = 0; state_index < state_count; ++state_index) {
    MVMCKrylovScaledAmplitudeResult amplitude_result;
    MVMCKrylovBoundedResult result;
    double complex amplitude_value;
    int vector_index;
    status = amplitude(&states[state_index], 1, amplitude_context,
                       &amplitude_result);
    if (collective_failure(status != MVMC_KRYLOV_STATUS_OK)) goto cleanup;
    setup_failure = !scaled_to_raw(&amplitude_result.value, &amplitude_value);
    if (collective_failure(setup_failure)) goto cleanup;
    if (amplitude_value != 0.0 ||
        (amplitude_result.value.state != MVMC_SCALED_COMPLEX_EXACT_ZERO &&
         amplitude_result.value.state != MVMC_SCALED_COMPLEX_NUMERIC_ZERO)) {
      all_amplitudes_zero = 0;
    }
    status = mvmc_bounded_krylov_evaluate(
        krylov_workspace, &states[state_index], 1, amplitude,
        amplitude_context, &result);
    setup_failure = status != MVMC_KRYLOV_STATUS_OK || !result.valid ||
                    result.evaluated_order != fixture->order + 1;
    if (collective_failure(setup_failure)) goto cleanup;
    for (vector_index = 0; vector_index < fixture->order + 2;
         ++vector_index) {
      setup_failure = !scaled_to_raw(
          &result.value[vector_index], &vectors[vector_index][state_index]);
      if (collective_failure(setup_failure)) goto cleanup;
      if (vectors[vector_index][state_index] != 0.0) {
        all_vectors_zero = 0;
      }
    }
    if (rank == 0) {
      printf("ROW %" PRIu64 " %d %" PRIu64 " %" PRIu64 " %" PRIu64
             " %" PRIu64 " %" PRIu64 " %" PRIu64,
             states[state_index], (int)amplitude_result.value.state,
             amplitude_result.well_pivoted_component_count,
             amplitude_result.near_pivot_component_count,
             amplitude_result.exact_zero_component_count,
             amplitude_result.numeric_zero_component_count,
             amplitude_result.local_factorization_count,
             amplitude_result.global_factorization_count);
      print_error_bound(amplitude_result.value.log_abs_error_bound);
      print_complex(amplitude_value);
      printf(" %d", fixture->order + 2);
      for (vector_index = 0; vector_index < fixture->order + 2;
           ++vector_index) {
        print_complex(vectors[vector_index][state_index]);
      }
      putchar('\n');
      printf("VSTATE %d", fixture->order + 2);
      for (vector_index = 0; vector_index < fixture->order + 2;
           ++vector_index) {
        printf(" %d", (int)result.value[vector_index].state);
        print_error_bound(result.value[vector_index].log_abs_error_bound);
      }
      putchar('\n');
    }
  }
  for (row = 0; row < dimension; ++row) {
    for (column = row; column < dimension; ++column) {
      const size_t entry = upper_index(dimension, row, column);
      overlap[entry] = inner_product(vectors[row], vectors[column],
                                     state_count);
      forward[entry] = inner_product(vectors[row], vectors[column + 1],
                                     state_count);
      reverse[entry] = inner_product(vectors[column], vectors[row + 1],
                                     state_count);
      squared[entry] = inner_product(vectors[row + 1],
                                     vectors[column + 1], state_count);
    }
  }
  if (fixture->support == 2) {
    int packed_zero = 1;
    size_t entry;
    for (entry = 0; entry < upper_count; ++entry) {
      if (overlap[entry] != 0.0 || forward[entry] != 0.0 ||
          reverse[entry] != 0.0 || squared[entry] != 0.0) {
        packed_zero = 0;
      }
    }
    if (!all_amplitudes_zero || !all_vectors_zero || !packed_zero) {
      if (rank == 0) {
        fprintf(stderr, "%s singular runtime predicate mismatch\n",
                fixture->fixture_id);
      }
      goto cleanup;
    }
  }
  if (rank == 0) {
    print_pack("S", overlap, upper_count);
    print_pack("KF", forward, upper_count);
    print_pack("KR", reverse, upper_count);
    print_pack("B", squared, upper_count);
  }
  for (cutoff_index = 0; cutoff_index < 4; ++cutoff_index) {
    if (fixture->support == 2) {
      if (rank == 0) {
        printf("SOLVE %s NQP_SINGULAR NOT_EVALUATED\n",
               CUTOFF_IDS[cutoff_index]);
      }
    } else {
      MVMCPowerLanczosGEVPPolicy policy;
      MVMCPowerLanczosGEVPResult result;
      MVMCKrylovGEVPStatus solve_status;
      if (mvmc_power_lanczos_gevp_default_policy(
              CUTOFFS[cutoff_index], &policy) != MVMC_KRYLOV_GEVP_OK) {
        goto cleanup;
      }
      if (fixture->arithmetic == 0) {
        double real_overlap[MAX_UPPER] = {0.0};
        double real_forward[MAX_UPPER] = {0.0};
        double real_reverse[MAX_UPPER] = {0.0};
        double real_squared[MAX_UPPER] = {0.0};
        size_t entry;
        for (entry = 0; entry < upper_count; ++entry) {
          if (cimag(overlap[entry]) != 0.0 ||
              cimag(forward[entry]) != 0.0 ||
              cimag(reverse[entry]) != 0.0 ||
              cimag(squared[entry]) != 0.0) {
            goto cleanup;
          }
          real_overlap[entry] = creal(overlap[entry]);
          real_forward[entry] = creal(forward[entry]);
          real_reverse[entry] = creal(reverse[entry]);
          real_squared[entry] = creal(squared[entry]);
        }
        solve_status = mvmc_power_lanczos_gevp_solve_real_packed(
            &policy, dimension, real_overlap, real_forward, real_reverse,
            real_squared, upper_count, &result);
      } else {
        solve_status = mvmc_power_lanczos_gevp_solve_complex_packed(
            &policy, dimension, overlap, forward, reverse, squared,
            upper_count, &result);
      }
      if (solve_status != MVMC_KRYLOV_GEVP_OK || !result.valid) {
        if (rank == 0) {
          fprintf(stderr, "%s %s P6 GEVP failed: %s\n",
                  fixture->fixture_id, CUTOFF_IDS[cutoff_index],
                  mvmc_krylov_gevp_status_string(solve_status));
        }
        goto cleanup;
      }
      {
        double complex optimized[MAX_STATES] = {0.0};
        double complex optimized_h[MAX_STATES] = {0.0};
        double complex energy_numerator;
        double state_norm;
        double state_energy;
        double state_m2;
        double state_variance = 0.0;
        int coefficient_index;
        for (state_index = 0; state_index < state_count; ++state_index) {
          for (coefficient_index = 0; coefficient_index < dimension;
               ++coefficient_index) {
            optimized[state_index] +=
                result.coefficient[coefficient_index] *
                vectors[coefficient_index][state_index];
            optimized_h[state_index] +=
                result.coefficient[coefficient_index] *
                vectors[coefficient_index + 1][state_index];
          }
        }
        state_norm = creal(inner_product(optimized, optimized, state_count));
        energy_numerator = inner_product(optimized, optimized_h, state_count);
        if (!(state_norm > 0.0) || !isfinite(state_norm) ||
            fabs(cimag(energy_numerator)) >
                1024.0 * DBL_EPSILON *
                    fmax(1.0, fabs(creal(energy_numerator)))) {
          goto cleanup;
        }
        state_energy = creal(energy_numerator) / state_norm;
        state_m2 = creal(inner_product(optimized_h, optimized_h,
                                      state_count)) / state_norm;
        for (state_index = 0; state_index < state_count; ++state_index) {
          const double complex residual =
              optimized_h[state_index] -
              state_energy * optimized[state_index];
          state_variance += creal(conj(residual) * residual);
        }
        state_variance /= state_norm;
        if (!isfinite(state_energy) || !isfinite(state_m2) ||
            !isfinite(state_variance) || state_variance < 0.0 ||
            !observable_consistent(state_energy, result.energy) ||
            !observable_consistent(state_m2, result.energy_squared)) {
          goto cleanup;
        }
      if (rank == 0) {
        printf("SOLVE %s OK %d %d %d %d %d",
               CUTOFF_IDS[cutoff_index], result.dimension,
               result.retained_rank, result.discarded_rank,
               result.root_multiplicity, result.phase_pivot);
        print_real(result.normalization);
        print_real(state_energy);
        print_real(state_m2);
        print_real(state_variance);
        print_real(result.energy);
        print_real(result.energy_squared);
        print_real(result.variance);
        print_real(result.normwise_backward_error);
        print_real(result.raw_action_relative_residual);
        print_real(result.root_gap);
        print_real(result.condition_estimate);
        printf(" %d %d", result.variance_clamped,
               result.observables_valid);
        putchar('\n');
        printf("COEFF %s %d", CUTOFF_IDS[cutoff_index], dimension);
        for (coefficient_index = 0; coefficient_index < dimension;
             ++coefficient_index) {
          print_complex(result.coefficient[coefficient_index]);
        }
        putchar('\n');
        for (state_index = 0; state_index < state_count; ++state_index) {
          printf("STATE %s %" PRIu64, CUTOFF_IDS[cutoff_index],
                 states[state_index]);
          print_complex(optimized[state_index]);
          print_complex(optimized_h[state_index]);
          putchar('\n');
        }
      }
      }
    }
  }
  if (rank == 0) printf("END_CASE %s\n", fixture->fixture_id);
  return_code = 1;

cleanup:
  mvmc_bounded_krylov_workspace_destroy(krylov_workspace);
  mvmc_bounded_krylov_plan_destroy(krylov_plan);
  mvmc_classic_krylov_complex_amplitude_workspace_destroy(complex_workspace);
  mvmc_classic_krylov_real_amplitude_workspace_destroy(real_workspace);
  mvmc_classic_krylov_model_workspace_destroy(model_workspace);
  return return_code;
}

static int emit_zero_packed_negative_probe(void) {
  double overlap[3] = {0.0, 0.0, 0.0};
  double forward[3] = {0.0, 0.0, 0.0};
  double reverse[3] = {0.0, 0.0, 0.0};
  double squared[3] = {0.0, 0.0, 0.0};
  MVMCKrylovGEVPPolicy policy;
  MVMCKrylovGEVPResult result;
  MVMCKrylovGEVPStatus status;
  if (mvmc_krylov_gevp_default_policy(0x1p-40, &policy) !=
      MVMC_KRYLOV_GEVP_OK) {
    return 0;
  }
  status = mvmc_krylov_gevp_solve_real_packed(
      &policy, 2, overlap, forward, reverse, squared, 3, &result);
  if (world_rank() == 0) {
    printf("NEGATIVE ZERO_PACKED_P4 %d %d\n", (int)status, result.valid);
  }
  return status == MVMC_KRYLOV_GEVP_NONPOSITIVE_OVERLAP_DIAGONAL &&
         !result.valid;
}

int main(int argc, char **argv) {
  FILE *input = NULL;
  char header[64];
  char manifest_sha[65];
  char registry_sha[65];
  char module_sha[65];
  char numerical_policy_sha[65];
  char final_token[64];
  int case_count = 0;
  int model_h_emitted[4] = {0, 0, 0, 0};
  int failure = 0;
  int ordinal;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#endif
  if (argc != 2) {
    if (world_rank() == 0) {
      fprintf(stderr, "usage: %s <canonical-input>\n", argv[0]);
    }
    failure = 1;
  } else {
    input = fopen(argv[1], "rb");
    failure = input == NULL;
  }
  if (collective_failure(failure)) goto cleanup;
  if (fscanf(input, "%63s %d %64s %64s %64s %64s", header, &case_count,
             manifest_sha, registry_sha, module_sha,
             numerical_policy_sha) != 6 ||
      strcmp(header, INPUT_HEADER) != 0 ||
      case_count != EXPECTED_CASE_COUNT ||
      strcmp(manifest_sha, FROZEN_MANIFEST_SHA) != 0 ||
      strcmp(registry_sha, FROZEN_REGISTRY_SHA) != 0 ||
      strcmp(module_sha, FROZEN_MODULE_SHA) != 0 ||
      strcmp(numerical_policy_sha, NUMERICAL_POLICY_SHA) != 0) {
    failure = 1;
  }
  if (collective_failure(failure)) goto cleanup;
  if (world_rank() == 0) {
    printf("%s %d %d %s %s %s %s\n", OUTPUT_HEADER, world_size(),
           case_count, manifest_sha, registry_sha, module_sha,
           numerical_policy_sha);
  }
  for (ordinal = 1; ordinal <= case_count; ++ordinal) {
    FixtureInput fixture;
    failure = !read_fixture(input, ordinal, &fixture);
    if (collective_failure(failure)) goto cleanup;
    failure = !evaluate_fixture(&fixture, model_h_emitted);
    if (collective_failure(failure)) goto cleanup;
  }
  failure = fscanf(input, "%63s", final_token) != 1 ||
            strcmp(final_token, "END_INPUT") != 0 ||
            fscanf(input, "%63s", final_token) == 1;
  if (collective_failure(failure)) goto cleanup;
  failure = !emit_zero_packed_negative_probe();
  if (collective_failure(failure)) goto cleanup;
  if (world_rank() == 0) printf("END_OUTPUT %d\n", case_count);

cleanup:
  if (input != NULL) fclose(input);
  if (failure && world_rank() == 0) {
    fprintf(stderr, "P6-C1 Cartesian driver failed\n");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return failure ? EXIT_FAILURE : EXIT_SUCCESS;
}
