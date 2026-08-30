#include "power_lanczos_classic_bridge.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "power_lanczos_classic_bridge.c requires bounded Krylov"
#endif

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  MVMCClassicKrylovTransfer *transfers;
  MVMCClassicKrylovSiteCoupling *coulomb_intra;
  MVMCClassicKrylovPairCoupling *coulomb_inter;
  MVMCClassicKrylovPairCoupling *hund;
  MVMCClassicKrylovPairCoupling *exchange;
} MVMCPowerLanczosClassicRawCopies;

struct MVMCPowerLanczosClassicBridge {
  int valid;
  MVMCPowerLanczosClassicArithmetic arithmetic;
  MVMCClassicKrylovModelWorkspace *model_workspace;
  MVMCClassicKrylovRealAmplitudeWorkspace *real_amplitude_workspace;
  MVMCClassicKrylovComplexAmplitudeWorkspace *complex_amplitude_workspace;
  size_t model_bytes;
  size_t amplitude_bytes;
  size_t allocated_bytes;
  uint64_t amplitude_generation_hash;
};

static int checked_multiply(size_t left, size_t right, size_t *product) {
  if (product == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *product = left * right;
  return 1;
}

static int checked_add(size_t left, size_t right, size_t *sum) {
  if (sum == NULL || right > SIZE_MAX - left) return 0;
  *sum = left + right;
  return 1;
}

static int valid_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static int valid_site(int site, int site_count) {
  return site >= 0 && site < site_count;
}

static int valid_model_rows(
    const MVMCPowerLanczosClassicView *view, int require_pointers) {
  int index;
  if (!require_pointers) return 1;
  if ((view->transfer_count != 0 &&
       (view->transfer_indices == NULL ||
        view->transfer_parameters == NULL)) ||
      (view->coulomb_intra_count != 0 &&
       (view->coulomb_intra_indices == NULL ||
        view->coulomb_intra_parameters == NULL)) ||
      (view->coulomb_inter_count != 0 &&
       (view->coulomb_inter_indices == NULL ||
        view->coulomb_inter_parameters == NULL)) ||
      (view->hund_count != 0 &&
       (view->hund_indices == NULL || view->hund_parameters == NULL)) ||
      (view->exchange_count != 0 &&
       (view->exchange_indices == NULL ||
        view->exchange_parameters == NULL))) {
    return 0;
  }
  for (index = 0; index < view->transfer_count; ++index) {
    const int *row = view->transfer_indices[index];
    if (row == NULL || !valid_site(row[0], view->site_count) ||
        (row[1] != 0 && row[1] != 1) ||
        !valid_site(row[2], view->site_count) ||
        (row[3] != 0 && row[3] != 1) ||
        !valid_complex(view->transfer_parameters[index])) {
      return 0;
    }
  }
  for (index = 0; index < view->coulomb_intra_count; ++index) {
    if (!valid_site(view->coulomb_intra_indices[index], view->site_count) ||
        !isfinite(view->coulomb_intra_parameters[index])) {
      return 0;
    }
  }
#define VALIDATE_PAIR_ROWS(rows, values, count)                            \
  do {                                                                     \
    for (index = 0; index < (count); ++index) {                            \
      if ((rows)[index] == NULL ||                                         \
          !valid_site((rows)[index][0], view->site_count) ||               \
          !valid_site((rows)[index][1], view->site_count) ||               \
          !isfinite((values)[index])) {                                    \
        return 0;                                                          \
      }                                                                    \
    }                                                                      \
  } while (0)
  VALIDATE_PAIR_ROWS(view->coulomb_inter_indices,
                     view->coulomb_inter_parameters,
                     view->coulomb_inter_count);
  VALIDATE_PAIR_ROWS(view->hund_indices, view->hund_parameters,
                     view->hund_count);
  VALIDATE_PAIR_ROWS(view->exchange_indices, view->exchange_parameters,
                     view->exchange_count);
#undef VALIDATE_PAIR_ROWS
  return 1;
}

static MVMCKrylovStatus validate_view_shape(
    const MVMCPowerLanczosClassicView *view, int require_model_pointers,
    int chain_rank, int chain_size) {
  long long expected_start;
  long long expected_end;
  if (view == NULL || view->site_count <= 0 ||
      view->site_count > INT_MAX / 2 ||
      view->up_electron_count <= 0 ||
      view->up_electron_count != view->down_electron_count ||
      view->up_electron_count > view->site_count ||
      (view->pure_spin != 0 && view->pure_spin != 1) ||
      (view->pure_spin &&
       view->up_electron_count + view->down_electron_count !=
           view->site_count) ||
      (view->arithmetic != MVMC_POWER_LANCZOS_CLASSIC_REAL &&
       view->arithmetic != MVMC_POWER_LANCZOS_CLASSIC_COMPLEX) ||
      view->transfer_count < 0 || view->coulomb_intra_count < 0 ||
      view->coulomb_inter_count < 0 || view->hund_count < 0 ||
      view->exchange_count < 0 || view->pair_hopping_count < 0 ||
      view->inter_all_count < 0 || view->nbody_inter_all_count < 0 ||
      view->nbody_g_count < 0 || view->qp_total <= 0 ||
      chain_rank < 0 || chain_size <= 0 || chain_rank >= chain_size ||
      view->nproj < 0 || view->ngutzwiller_idx < 0 ||
      view->njastrow_idx < 0 || view->nspin_jastrow_idx < 0 ||
      view->ndoublon_holon_2site_idx < 0 ||
      view->ndoublon_holon_4site_idx < 0 ||
      !isfinite(view->scaled_pivot_tolerance) ||
      view->qp_weights == NULL ||
      (view->arithmetic == MVMC_POWER_LANCZOS_CLASSIC_REAL &&
       view->slater_real == NULL) ||
      (view->arithmetic == MVMC_POWER_LANCZOS_CLASSIC_COMPLEX &&
       view->slater_complex == NULL)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (view->unsupported_amplitude_features != 0 ||
      view->pair_hopping_count != 0 || view->inter_all_count != 0 ||
      view->nbody_inter_all_count != 0 || view->nbody_g_count != 0 ||
      (view->pure_spin &&
       (view->transfer_count != 0 || view->coulomb_intra_count != 0)) ||
      (!view->pure_spin && view->exchange_count != 0)) {
    return MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL;
  }
  expected_start =
      (long long)view->qp_total * (long long)chain_rank / chain_size;
  expected_end =
      (long long)view->qp_total * (long long)(chain_rank + 1) / chain_size;
  if (expected_start < 0 || expected_end < expected_start ||
      expected_end > INT_MAX || view->qp_start != (int)expected_start ||
      view->qp_end != (int)expected_end ||
      !valid_model_rows(view, require_model_pointers)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static uint64_t hash_u64(uint64_t hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    hash ^= (value >> (8 * byte)) & UINT64_C(0xff);
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

static uint64_t view_metadata_hash(
    const MVMCPowerLanczosClassicView *view) {
  const int values[] = {
      view->site_count,
      view->up_electron_count,
      view->down_electron_count,
      view->pure_spin,
      (int)view->arithmetic,
      view->transfer_count,
      view->coulomb_intra_count,
      view->coulomb_inter_count,
      view->hund_count,
      view->exchange_count,
      view->pair_hopping_count,
      view->inter_all_count,
      view->nbody_inter_all_count,
      view->nbody_g_count,
      view->qp_total,
      view->nproj,
      view->ngutzwiller_idx,
      view->njastrow_idx,
      view->nspin_jastrow_idx,
      view->ndoublon_holon_2site_idx,
      view->ndoublon_holon_4site_idx};
  uint64_t hash = UINT64_C(1469598103934665603);
  uint64_t pivot = 0;
  size_t index;
  memcpy(&pivot, &view->scaled_pivot_tolerance, sizeof(pivot));
  hash = hash_u64(hash, MVMC_POWER_LANCZOS_CLASSIC_BRIDGE_VERSION);
  for (index = 0; index < sizeof(values) / sizeof(values[0]); ++index) {
    hash = hash_u64(hash, (uint64_t)(uint32_t)values[index]);
  }
  hash = hash_u64(hash, pivot);
  hash = hash_u64(hash, (uint64_t)view->unsupported_amplitude_features);
  return hash == 0 ? UINT64_C(1) : hash;
}

static MVMCKrylovStatus synchronize_status(
    MVMCKrylovStatus local_status,
    MVMCClassicPfaffianCommunicator communicator) {
#ifdef _mpi_use
  int local = (int)local_status;
  int global = 0;
  if (MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX,
                    communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (global < (int)MVMC_KRYLOV_STATUS_OK ||
      global > (int)MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  return (MVMCKrylovStatus)global;
#else
  (void)communicator;
  return local_status;
#endif
}

static MVMCKrylovStatus audit_world_metadata(
    const MVMCPowerLanczosClassicView *view,
    MVMCClassicPfaffianCommunicator world_communicator) {
#ifdef _mpi_use
  const uint64_t local = view_metadata_hash(view);
  uint64_t minimum = 0;
  uint64_t maximum = 0;
  if (MPI_Allreduce(&local, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                    world_communicator) != MPI_SUCCESS ||
      MPI_Allreduce(&local, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                    world_communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  return minimum == maximum ? MVMC_KRYLOV_STATUS_OK
                            : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
#else
  (void)view;
  (void)world_communicator;
  return MVMC_KRYLOV_STATUS_OK;
#endif
}

static int allocate_array(int count, size_t element_size, void **output) {
  size_t bytes;
  if (output == NULL) return 0;
  *output = NULL;
  if (count == 0) return 1;
  if (count < 0 || !checked_multiply((size_t)count, element_size, &bytes)) {
    return 0;
  }
  *output = calloc(1, bytes);
  return *output != NULL;
}

static void free_raw_copies(MVMCPowerLanczosClassicRawCopies *copies) {
  if (copies == NULL) return;
  free(copies->exchange);
  free(copies->hund);
  free(copies->coulomb_inter);
  free(copies->coulomb_intra);
  free(copies->transfers);
  memset(copies, 0, sizeof(*copies));
}

static int build_root_raw(
    const MVMCPowerLanczosClassicView *view,
    MVMCPowerLanczosClassicRawCopies *copies,
    MVMCClassicKrylovRawModel *raw) {
  int index;
  memset(copies, 0, sizeof(*copies));
  memset(raw, 0, sizeof(*raw));
  if (!allocate_array(view->transfer_count, sizeof(*copies->transfers),
                      (void **)&copies->transfers) ||
      !allocate_array(view->coulomb_intra_count,
                      sizeof(*copies->coulomb_intra),
                      (void **)&copies->coulomb_intra) ||
      !allocate_array(view->coulomb_inter_count,
                      sizeof(*copies->coulomb_inter),
                      (void **)&copies->coulomb_inter) ||
      !allocate_array(view->hund_count, sizeof(*copies->hund),
                      (void **)&copies->hund) ||
      !allocate_array(view->exchange_count, sizeof(*copies->exchange),
                      (void **)&copies->exchange)) {
    free_raw_copies(copies);
    return 0;
  }
  for (index = 0; index < view->transfer_count; ++index) {
    const int *row = view->transfer_indices[index];
    copies->transfers[index].output_site = row[0];
    copies->transfers[index].output_spin = row[1];
    copies->transfers[index].input_site = row[2];
    copies->transfers[index].input_spin = row[3];
    copies->transfers[index].coefficient =
        view->transfer_parameters[index];
  }
  for (index = 0; index < view->coulomb_intra_count; ++index) {
    copies->coulomb_intra[index].site =
        view->coulomb_intra_indices[index];
    copies->coulomb_intra[index].coefficient =
        view->coulomb_intra_parameters[index];
  }
#define COPY_PAIR_ROWS(destination, rows, values, count)                  \
  do {                                                                     \
    for (index = 0; index < (count); ++index) {                            \
      (destination)[index].first_site = (rows)[index][0];                  \
      (destination)[index].second_site = (rows)[index][1];                 \
      (destination)[index].coefficient = (values)[index];                  \
    }                                                                      \
  } while (0)
  COPY_PAIR_ROWS(copies->coulomb_inter, view->coulomb_inter_indices,
                 view->coulomb_inter_parameters, view->coulomb_inter_count);
  COPY_PAIR_ROWS(copies->hund, view->hund_indices, view->hund_parameters,
                 view->hund_count);
  COPY_PAIR_ROWS(copies->exchange, view->exchange_indices,
                 view->exchange_parameters, view->exchange_count);
#undef COPY_PAIR_ROWS
  raw->site_count = view->site_count;
  raw->up_electron_count = view->up_electron_count;
  raw->down_electron_count = view->down_electron_count;
  raw->pure_spin = view->pure_spin;
  raw->transfer_count = view->transfer_count;
  raw->transfers = copies->transfers;
  raw->coulomb_intra_count = view->coulomb_intra_count;
  raw->coulomb_intra = copies->coulomb_intra;
  raw->coulomb_inter_count = view->coulomb_inter_count;
  raw->coulomb_inter = copies->coulomb_inter;
  raw->hund_count = view->hund_count;
  raw->hund = copies->hund;
  raw->exchange_count = view->exchange_count;
  raw->exchange = copies->exchange;
  raw->pair_hopping_count = view->pair_hopping_count;
  raw->inter_all_count = view->inter_all_count;
  raw->nbody_inter_all_count = view->nbody_inter_all_count;
  return 1;
}

static MVMCClassicKrylovAmplitudeLayout amplitude_layout(
    const MVMCPowerLanczosClassicView *view,
    MVMCClassicPfaffianCommunicator chain_communicator) {
  MVMCClassicKrylovAmplitudeLayout layout;
  memset(&layout, 0, sizeof(layout));
  layout.site_count = (size_t)view->site_count;
  layout.up_electron_count = (size_t)view->up_electron_count;
  layout.down_electron_count = (size_t)view->down_electron_count;
  layout.pure_spin = view->pure_spin;
  layout.qp_total = view->qp_total;
  layout.qp_start = view->qp_start;
  layout.qp_end = view->qp_end;
  layout.scaled_pivot_tolerance = view->scaled_pivot_tolerance;
  layout.nproj = view->nproj;
  layout.ngutzwiller_idx = view->ngutzwiller_idx;
  layout.njastrow_idx = view->njastrow_idx;
  layout.nspin_jastrow_idx = view->nspin_jastrow_idx;
  layout.ndoublon_holon_2site_idx = view->ndoublon_holon_2site_idx;
  layout.ndoublon_holon_4site_idx = view->ndoublon_holon_4site_idx;
  layout.gutzwiller_idx = view->gutzwiller_idx;
  layout.jastrow_idx = (const int *const *)view->jastrow_idx;
  layout.spin_jastrow_idx =
      (const int *const *)view->spin_jastrow_idx;
  layout.doublon_holon_2site_idx =
      (const int *const *)view->doublon_holon_2site_idx;
  layout.doublon_holon_4site_idx =
      (const int *const *)view->doublon_holon_4site_idx;
  layout.projection_parameters = view->projection_parameters;
  layout.unsupported_features = view->unsupported_amplitude_features;
  layout.communicator = chain_communicator;
  return layout;
}

MVMCKrylovStatus mvmc_power_lanczos_classic_bridge_create(
    const MVMCPowerLanczosClassicView *view,
    MVMCClassicPfaffianCommunicator world_communicator,
    MVMCClassicPfaffianCommunicator chain_communicator,
    MVMCPowerLanczosClassicBridge **bridge) {
  MVMCPowerLanczosClassicBridge *created = NULL;
  MVMCPowerLanczosClassicRawCopies copies;
  MVMCClassicKrylovRawModel raw;
  MVMCClassicKrylovAmplitudeLayout layout;
  const MVMCKrylovFockModel *model = NULL;
  MVMCKrylovStatus status;
  int world_rank = 0;
  int chain_rank = 0;
  int chain_size = 1;
  int root_raw_ready = 1;
  uint64_t minimum_generation = 0;
  uint64_t maximum_generation = 0;
  memset(&copies, 0, sizeof(copies));
  memset(&raw, 0, sizeof(raw));
  if (bridge == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *bridge = NULL;
#ifdef _mpi_use
  if (world_communicator == MPI_COMM_NULL ||
      chain_communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(world_communicator, &world_rank) != MPI_SUCCESS ||
      MPI_Comm_rank(chain_communicator, &chain_rank) != MPI_SUCCESS ||
      MPI_Comm_size(chain_communicator, &chain_size) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#else
  (void)world_communicator;
  (void)chain_communicator;
#endif
  status = validate_view_shape(view, world_rank == 0,
                               chain_rank, chain_size);
  status = synchronize_status(status, world_communicator);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  status = audit_world_metadata(view, world_communicator);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (world_rank == 0) root_raw_ready = build_root_raw(view, &copies, &raw);
  status = synchronize_status(
      root_raw_ready ? MVMC_KRYLOV_STATUS_OK
                     : MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE,
      world_communicator);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    free_raw_copies(&copies);
    return status;
  }
  created = (MVMCPowerLanczosClassicBridge *)calloc(1, sizeof(*created));
  status = synchronize_status(
      created == NULL ? MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE
                      : MVMC_KRYLOV_STATUS_OK,
      world_communicator);
  if (status != MVMC_KRYLOV_STATUS_OK || created == NULL) {
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    goto fail;
  }
  status = mvmc_classic_krylov_model_workspace_create_from_root(
      world_rank == 0 ? &raw : NULL, world_communicator,
      &created->model_workspace);
  free_raw_copies(&copies);
  if (status != MVMC_KRYLOV_STATUS_OK) goto fail;
  model = mvmc_classic_krylov_model(created->model_workspace);
  if (model == NULL) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    goto fail;
  }
  layout = amplitude_layout(view, chain_communicator);
  if (view->arithmetic == MVMC_POWER_LANCZOS_CLASSIC_REAL) {
    status = mvmc_classic_krylov_real_amplitude_workspace_create(
        &layout, view->slater_real, view->qp_weights,
        &created->real_amplitude_workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      created->amplitude_generation_hash =
          mvmc_classic_krylov_real_amplitude_generation_hash(
              created->real_amplitude_workspace);
      created->amplitude_bytes =
          mvmc_classic_krylov_real_amplitude_workspace_bytes(
              created->real_amplitude_workspace);
    }
  } else {
    status = mvmc_classic_krylov_complex_amplitude_workspace_create(
        &layout, view->slater_complex, view->qp_weights,
        &created->complex_amplitude_workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      created->amplitude_generation_hash =
          mvmc_classic_krylov_complex_amplitude_generation_hash(
              created->complex_amplitude_workspace);
      created->amplitude_bytes =
          mvmc_classic_krylov_complex_amplitude_workspace_bytes(
              created->complex_amplitude_workspace);
    }
  }
  status = synchronize_status(status, world_communicator);
  if (status != MVMC_KRYLOV_STATUS_OK) goto fail;
#ifdef _mpi_use
  if (MPI_Allreduce(&created->amplitude_generation_hash,
                    &minimum_generation, 1, MPI_UINT64_T, MPI_MIN,
                    world_communicator) != MPI_SUCCESS ||
      MPI_Allreduce(&created->amplitude_generation_hash,
                    &maximum_generation, 1, MPI_UINT64_T, MPI_MAX,
                    world_communicator) != MPI_SUCCESS) {
    status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    goto fail;
  }
#else
  minimum_generation = created->amplitude_generation_hash;
  maximum_generation = created->amplitude_generation_hash;
#endif
  if (minimum_generation == 0 ||
      minimum_generation != maximum_generation) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    goto fail;
  }
  created->model_bytes =
      mvmc_classic_krylov_model_workspace_bytes(created->model_workspace);
  if (created->model_bytes == 0 || created->amplitude_bytes == 0 ||
      !checked_add(sizeof(*created), created->model_bytes,
                   &created->allocated_bytes) ||
      !checked_add(created->allocated_bytes, created->amplitude_bytes,
                   &created->allocated_bytes)) {
    status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
    goto fail;
  }
  created->arithmetic = view->arithmetic;
  created->valid = 1;
  *bridge = created;
  return MVMC_KRYLOV_STATUS_OK;

fail:
  free_raw_copies(&copies);
  mvmc_power_lanczos_classic_bridge_destroy(created);
  return status;
}

void mvmc_power_lanczos_classic_bridge_destroy(
    MVMCPowerLanczosClassicBridge *bridge) {
  if (bridge == NULL) return;
  mvmc_classic_krylov_complex_amplitude_workspace_destroy(
      bridge->complex_amplitude_workspace);
  mvmc_classic_krylov_real_amplitude_workspace_destroy(
      bridge->real_amplitude_workspace);
  mvmc_classic_krylov_model_workspace_destroy(bridge->model_workspace);
  memset(bridge, 0, sizeof(*bridge));
  free(bridge);
}

const MVMCKrylovFockModel *mvmc_power_lanczos_classic_bridge_model(
    const MVMCPowerLanczosClassicBridge *bridge) {
  if (bridge == NULL || !bridge->valid) return NULL;
  return mvmc_classic_krylov_model(bridge->model_workspace);
}

MVMCKrylovScaledAmplitudeCallback
mvmc_power_lanczos_classic_bridge_amplitude(
    const MVMCPowerLanczosClassicBridge *bridge) {
  if (bridge == NULL || !bridge->valid) return NULL;
  return bridge->arithmetic == MVMC_POWER_LANCZOS_CLASSIC_REAL
             ? mvmc_classic_krylov_real_scaled_amplitude
             : mvmc_classic_krylov_complex_scaled_amplitude;
}

void *mvmc_power_lanczos_classic_bridge_amplitude_context(
    const MVMCPowerLanczosClassicBridge *bridge) {
  if (bridge == NULL || !bridge->valid) return NULL;
  return bridge->arithmetic == MVMC_POWER_LANCZOS_CLASSIC_REAL
             ? (void *)bridge->real_amplitude_workspace
             : (void *)bridge->complex_amplitude_workspace;
}

uint64_t mvmc_power_lanczos_classic_bridge_generation_hash(
    const MVMCPowerLanczosClassicBridge *bridge) {
  return bridge == NULL || !bridge->valid
             ? 0
             : bridge->amplitude_generation_hash;
}

MVMCKrylovStatus mvmc_power_lanczos_classic_bridge_summary(
    const MVMCPowerLanczosClassicBridge *bridge,
    MVMCPowerLanczosClassicBridgeSummary *summary) {
  MVMCPowerLanczosClassicBridgeSummary candidate;
  const MVMCKrylovFockModel *model;
  if (summary == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(summary, 0, sizeof(*summary));
  if (bridge == NULL || !bridge->valid) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  model = mvmc_power_lanczos_classic_bridge_model(bridge);
  if (model == NULL) return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.version = MVMC_POWER_LANCZOS_CLASSIC_BRIDGE_VERSION;
  candidate.arithmetic = bridge->arithmetic;
  candidate.model_bytes = bridge->model_bytes;
  candidate.amplitude_bytes = bridge->amplitude_bytes;
  candidate.allocated_bytes = bridge->allocated_bytes;
  candidate.term_count = model->term_count;
  candidate.operator_count = model->operator_count;
  candidate.amplitude_generation_hash =
      bridge->amplitude_generation_hash;
  *summary = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}
