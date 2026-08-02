#include "classic_krylov_amplitude.h"
#include "classic_krylov_model.h"
#include "krylov_fock_reference.h"

#include <complex.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

enum { SITE_COUNT = 4, ORBITAL_COUNT = 8 };

typedef struct {
  MVMCClassicKrylovRawModel raw;
  MVMCClassicKrylovTransfer transfers[20];
  MVMCClassicKrylovSiteCoupling intra[SITE_COUNT];
  MVMCClassicKrylovPairCoupling inter[2];
  MVMCClassicKrylovPairCoupling hund[4];
  MVMCClassicKrylovPairCoupling exchange[4];
} FixtureModel;

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

static int bit_count(unsigned int value) {
  int count = 0;
  while (value != 0U) {
    count += (int)(value & 1U);
    value >>= 1;
  }
  return count;
}

static void initialize_electronic_model(FixtureModel *fixture,
                                        int complex_case) {
  static const int bonds[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  static const double intra_values[SITE_COUNT] = {0.25, -0.125, 0.375,
                                                   -0.3125};
  int transfer_count = 0;
  int bond;
  int spin;

  memset(fixture, 0, sizeof(*fixture));
  for (bond = 0; bond < 4; ++bond) {
    const int left = bonds[bond][0];
    const int right = bonds[bond][1];
    const double complex hopping =
        complex_case ? (0.75 + (0.125 + 0.0625 * (double)bond) * I)
                     : (0.75 + 0.0 * I);
    for (spin = 0; spin < 2; ++spin) {
      fixture->transfers[transfer_count++] =
          (MVMCClassicKrylovTransfer){left, spin, right, spin, hopping};
      fixture->transfers[transfer_count++] =
          (MVMCClassicKrylovTransfer){right, spin, left, spin,
                                      conj(hopping)};
    }
  }
  fixture->transfers[transfer_count++] =
      (MVMCClassicKrylovTransfer){1, 0, 1, 0, 0.4375};
  for (spin = 0; spin < SITE_COUNT; ++spin) {
    fixture->intra[spin] =
        (MVMCClassicKrylovSiteCoupling){spin, intra_values[spin]};
  }
  fixture->inter[0] = (MVMCClassicKrylovPairCoupling){0, 2, 0.1875};
  fixture->inter[1] = (MVMCClassicKrylovPairCoupling){1, 3, -0.15625};
  fixture->hund[0] = (MVMCClassicKrylovPairCoupling){0, 1, 0.21875};
  fixture->hund[1] = (MVMCClassicKrylovPairCoupling){2, 3, -0.09375};

  fixture->raw.site_count = SITE_COUNT;
  fixture->raw.up_electron_count = 2;
  fixture->raw.down_electron_count = 2;
  fixture->raw.transfer_count = transfer_count;
  fixture->raw.transfers = fixture->transfers;
  fixture->raw.coulomb_intra_count = SITE_COUNT;
  fixture->raw.coulomb_intra = fixture->intra;
  fixture->raw.coulomb_inter_count = 2;
  fixture->raw.coulomb_inter = fixture->inter;
  fixture->raw.hund_count = 2;
  fixture->raw.hund = fixture->hund;
}

static void initialize_spin_model(FixtureModel *fixture) {
  static const double hund_values[4] = {0.375, -0.25, 0.1875, -0.125};
  static const double exchange_values[4] = {0.625, -0.3125, 0.4375,
                                             -0.5625};
  int bond;

  memset(fixture, 0, sizeof(*fixture));
  for (bond = 0; bond < 4; ++bond) {
    fixture->hund[bond] =
        (MVMCClassicKrylovPairCoupling){bond, (bond + 1) % 4,
                                        hund_values[bond]};
    fixture->exchange[bond] =
        (MVMCClassicKrylovPairCoupling){bond, (bond + 1) % 4,
                                        exchange_values[bond]};
  }
  fixture->raw.site_count = SITE_COUNT;
  fixture->raw.up_electron_count = 2;
  fixture->raw.down_electron_count = 2;
  fixture->raw.pure_spin = 1;
  fixture->raw.hund_count = 4;
  fixture->raw.hund = fixture->hund;
  fixture->raw.exchange_count = 4;
  fixture->raw.exchange = fixture->exchange;
}

static void set_real_pair(double *slater, int up, int down, double value) {
  const int down_orbital = SITE_COUNT + down;
  slater[(size_t)down_orbital * ORBITAL_COUNT + (size_t)up] = -value;
  slater[(size_t)up * ORBITAL_COUNT + (size_t)down_orbital] = value;
}

static void set_complex_pair(double complex *slater, int up, int down,
                             double complex value) {
  const int down_orbital = SITE_COUNT + down;
  slater[(size_t)down_orbital * ORBITAL_COUNT + (size_t)up] = -value;
  slater[(size_t)up * ORBITAL_COUNT + (size_t)down_orbital] = value;
}

static void initialize_real_slater(double *slater) {
  static const double pairs[SITE_COUNT][SITE_COUNT] = {
      {1.5, -0.5, 1.5, 1.0},
      {-1.5, 0.5, -1.5, -1.0},
      {-0.5, 1.5, 1.5, 2.0},
      {-2.0, -2.0, -1.0, 0.5}};
  int up;
  int down;
  memset(slater, 0, ORBITAL_COUNT * ORBITAL_COUNT * sizeof(*slater));
  for (up = 0; up < SITE_COUNT; ++up) {
    for (down = 0; down < SITE_COUNT; ++down) {
      set_real_pair(slater, up, down, pairs[up][down]);
    }
  }
}

static void initialize_complex_slater(double complex *slater) {
  double complex pairs[SITE_COUNT][SITE_COUNT] = {
      {2.0 + I, 1.0 + 1.5 * I, 0.5 - 0.5 * I, 0.5 - 1.5 * I},
      {-2.0 - I, -1.0 - 1.5 * I, -0.5 + 0.5 * I, -0.5 + 1.5 * I},
      {0.5 + I, 0.5 + 2.0 * I, -0.5 - 2.0 * I, -0.5 - 1.5 * I},
      {1.5 + 1.5 * I, 0.5 + 0.5 * I, 1.5 - 1.5 * I,
       1.5 + 1.5 * I}};
  int up;
  int down;
  memset(slater, 0, ORBITAL_COUNT * ORBITAL_COUNT * sizeof(*slater));
  for (up = 0; up < SITE_COUNT; ++up) {
    for (down = 0; down < SITE_COUNT; ++down) {
      set_complex_pair(slater, up, down, pairs[up][down]);
    }
  }
}

static MVMCClassicKrylovAmplitudeLayout amplitude_layout(int pure_spin) {
  static const int gutzwiller_indices[SITE_COUNT] = {0, 0, 0, 0};
  static const double complex projection_parameter[1] = {-0.27};
  MVMCClassicKrylovAmplitudeLayout layout;
  const int rank = world_rank();
  const int size = world_size();
  memset(&layout, 0, sizeof(layout));
  layout.site_count = SITE_COUNT;
  layout.up_electron_count = 2;
  layout.down_electron_count = 2;
  layout.pure_spin = pure_spin;
  layout.qp_total = 1;
  layout.qp_start = rank / size;
  layout.qp_end = (rank + 1) / size;
  layout.scaled_pivot_tolerance = 64.0;
  if (!pure_spin) {
    layout.nproj = 1;
    layout.ngutzwiller_idx = 1;
    layout.gutzwiller_idx = gutzwiller_indices;
    layout.projection_parameters = projection_parameter;
  }
  layout.communicator = world_communicator();
  return layout;
}

static MVMCKrylovLimits fixture_limits(void) {
  MVMCKrylovLimits limits;
  limits.max_states = 4096;
  limits.max_transitions = 262144;
  limits.max_amplitude_evaluations = 4096;
  limits.max_bytes = 64 * 1024 * 1024;
  limits.max_order = 3;
  return limits;
}

static int evaluate_case(const char *name, int pure_spin, int complex_case) {
  FixtureModel fixture;
  MVMCClassicKrylovModelWorkspace *model_workspace = NULL;
  MVMCClassicKrylovRealAmplitudeWorkspace *real_workspace = NULL;
  MVMCClassicKrylovComplexAmplitudeWorkspace *complex_workspace = NULL;
  MVMCKrylovWorkspace *krylov_workspace = NULL;
  const MVMCKrylovFockModel *model;
  MVMCKrylovAmplitudeCallback callback;
  void *amplitude_context;
  MVMCKrylovLimits limits = fixture_limits();
  MVMCKrylovStatus status;
  MVMCClassicKrylovAmplitudeLayout layout = amplitude_layout(pure_spin);
  double real_slater[ORBITAL_COUNT * ORBITAL_COUNT];
  double complex complex_slater[ORBITAL_COUNT * ORBITAL_COUNT];
  const double complex weights[1] = {1.0};
  unsigned int up_mask;
  unsigned int down_mask;
  const int rank = world_rank();
  const int size = world_size();
  int return_code = 1;

  if (pure_spin) {
    initialize_spin_model(&fixture);
  } else {
    initialize_electronic_model(&fixture, complex_case);
  }
  status = mvmc_classic_krylov_model_workspace_create_from_root(
      rank == 0 ? &fixture.raw : NULL, world_communicator(),
      &model_workspace);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    fprintf(stderr, "%s model create failed on rank %d: %s\n", name, rank,
            mvmc_krylov_status_string(status));
    goto cleanup;
  }
  if (complex_case) {
    initialize_complex_slater(complex_slater);
    status = mvmc_classic_krylov_complex_amplitude_workspace_create(
        &layout, complex_slater, weights, &complex_workspace);
    callback = mvmc_classic_krylov_complex_amplitude;
    amplitude_context = complex_workspace;
  } else {
    initialize_real_slater(real_slater);
    status = mvmc_classic_krylov_real_amplitude_workspace_create(
        &layout, real_slater, weights, &real_workspace);
    callback = mvmc_classic_krylov_real_amplitude;
    amplitude_context = real_workspace;
  }
  if (status != MVMC_KRYLOV_STATUS_OK || amplitude_context == NULL) {
    fprintf(stderr, "%s amplitude create failed on rank %d: %s\n", name,
            rank, mvmc_krylov_status_string(status));
    goto cleanup;
  }
  model = mvmc_classic_krylov_model(model_workspace);
  krylov_workspace =
      mvmc_krylov_workspace_create(SITE_COUNT, &limits, &status);
  if (krylov_workspace == NULL || status != MVMC_KRYLOV_STATUS_OK) {
    fprintf(stderr, "%s Krylov workspace create failed: %s\n", name,
            mvmc_krylov_status_string(status));
    goto cleanup;
  }
  if (rank == 0) {
    printf("BEGIN %s rank_count %d empty_qp_ranks %d\n", name, size,
           size - 1);
  }
  for (up_mask = 0; up_mask < (1U << SITE_COUNT); ++up_mask) {
    if (bit_count(up_mask) != 2) continue;
    for (down_mask = 0; down_mask < (1U << SITE_COUNT); ++down_mask) {
      uint64_t configuration;
      MVMCKrylovResult result;
      int order;
      if (bit_count(down_mask) != 2) continue;
      if (pure_spin && down_mask != ((~up_mask) & 0x0fU)) continue;
      configuration = (uint64_t)up_mask |
                      ((uint64_t)down_mask << SITE_COUNT);
      status = mvmc_krylov_evaluate(
          krylov_workspace, model, &configuration, 1, callback,
          amplitude_context, &result);
      if (status != MVMC_KRYLOV_STATUS_OK || !result.valid ||
          result.evaluated_order != 3) {
        fprintf(stderr, "%s state=%" PRIu64 " failed on rank %d: %s\n",
                name, configuration, rank,
                mvmc_krylov_status_string(status));
        goto cleanup;
      }
      if (rank == 0) {
        printf("ROW %" PRIu64, configuration);
        for (order = 0; order <= 3; ++order) {
          printf(" %.17g %.17g", creal(result.value[order]),
                 cimag(result.value[order]));
        }
        putchar('\n');
      }
    }
  }
  if (rank == 0) printf("END %s\n", name);
  return_code = 0;

cleanup:
  mvmc_krylov_workspace_destroy(krylov_workspace);
  mvmc_classic_krylov_complex_amplitude_workspace_destroy(complex_workspace);
  mvmc_classic_krylov_real_amplitude_workspace_destroy(real_workspace);
  mvmc_classic_krylov_model_workspace_destroy(model_workspace);
  return return_code;
}

int main(int argc, char **argv) {
  int failed = 0;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  failed |= evaluate_case("electronic_real", 0, 0);
  failed |= evaluate_case("electronic_complex", 0, 1);
  failed |= evaluate_case("spin_real", 1, 0);
  failed |= evaluate_case("spin_complex", 1, 1);
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
