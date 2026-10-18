#include "power_lanczos_gc_bridge.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "power_lanczos_gc_bridge.c requires bounded Krylov"
#endif

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct MVMCPowerLanczosGCBridge {
  int valid;
  MVMCClassicKrylovModelWorkspace *model_workspace;
  MVMCGCKrylovAmplitudeWorkspace *amplitude_workspace;
  uint64_t amplitude_generation_hash;
};

static MVMCKrylovStatus validate_gc_view_shape(
    const MVMCPowerLanczosGCView *view, int chain_rank, int chain_size) {
  const MVMCPowerLanczosClassicView *data;
  long long expected_start;
  long long expected_end;
  int index;
  if (view == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  data = &view->data;
  if (data->site_count <= 0 || data->site_count > INT_MAX / 2 ||
      data->up_electron_count != 0 || data->down_electron_count != 0 ||
      data->pure_spin != 0 ||
      data->arithmetic != MVMC_POWER_LANCZOS_CLASSIC_COMPLEX ||
      (view->domain != MVMC_KRYLOV_DOMAIN_GC_BALANCED &&
       view->domain != MVMC_KRYLOV_DOMAIN_GC_EVEN) ||
      view->translation_count <= 0 || view->translation_weights == NULL ||
      (view->spin_triplet_s0 != 0 && view->spin_triplet_s0 != 1) ||
      (view->spin_triplet_s0 &&
       view->domain != MVMC_KRYLOV_DOMAIN_GC_BALANCED) ||
      data->transfer_count < 0 || data->coulomb_intra_count < 0 ||
      data->coulomb_inter_count < 0 || data->hund_count < 0 ||
      data->exchange_count < 0 || data->pair_hopping_count < 0 ||
      data->inter_all_count < 0 || data->nbody_inter_all_count < 0 ||
      data->nbody_g_count < 0 || data->qp_total <= 0 ||
      chain_rank < 0 || chain_size <= 0 || chain_rank >= chain_size ||
      data->nproj < 0 || data->ngutzwiller_idx < 0 ||
      data->njastrow_idx < 0 || data->nspin_jastrow_idx < 0 ||
      data->ndoublon_holon_2site_idx < 0 ||
      data->ndoublon_holon_4site_idx < 0 ||
      !isfinite(data->scaled_pivot_tolerance) || data->qp_weights == NULL ||
      data->slater_complex == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (index = 0; index < view->translation_count; ++index) {
    if (!isfinite(creal(view->translation_weights[index])) ||
        !isfinite(cimag(view->translation_weights[index]))) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  if (data->unsupported_amplitude_features != 0 ||
      data->pair_hopping_count != 0 || data->inter_all_count != 0 ||
      data->nbody_inter_all_count != 0 || data->nbody_g_count != 0 ||
      data->exchange_count != 0) {
    return MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL;
  }
  expected_start =
      (long long)data->qp_total * (long long)chain_rank / chain_size;
  expected_end =
      (long long)data->qp_total * (long long)(chain_rank + 1) / chain_size;
  if (expected_start < 0 || expected_end < expected_start ||
      expected_end > INT_MAX || data->qp_start != (int)expected_start ||
      data->qp_end != (int)expected_end) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus synchronize_status(
    MVMCKrylovStatus local_status,
    MVMCClassicPfaffianCommunicator communicator) {
#ifdef _mpi_use
  int local = (int)local_status;
  int global = 0;
  if (MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, communicator) !=
          MPI_SUCCESS ||
      global < (int)MVMC_KRYLOV_STATUS_OK ||
      global > (int)MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  return (MVMCKrylovStatus)global;
#else
  (void)communicator;
  return local_status;
#endif
}

#ifdef _mpi_use
static uint64_t hash_u64(uint64_t hash, uint64_t value) {
  int byte;
  for (byte = 0; byte < 8; ++byte) {
    hash ^= (value >> (8 * byte)) & UINT64_C(0xff);
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

static uint64_t gc_view_metadata_hash(const MVMCPowerLanczosGCView *view) {
  const MVMCPowerLanczosClassicView *data = &view->data;
  const int values[] = {
      data->site_count, (int)view->domain, view->translation_count,
      view->spin_triplet_s0,
      data->transfer_count, data->coulomb_intra_count,
      data->coulomb_inter_count, data->hund_count, data->qp_total,
      data->nproj, data->ngutzwiller_idx, data->njastrow_idx,
      data->nspin_jastrow_idx, data->ndoublon_holon_2site_idx,
      data->ndoublon_holon_4site_idx};
  uint64_t hash = UINT64_C(1469598103934665603);
  uint64_t bits = 0;
  size_t index;
  hash = hash_u64(hash, MVMC_POWER_LANCZOS_GC_BRIDGE_VERSION);
  for (index = 0; index < sizeof(values) / sizeof(values[0]); ++index) {
    hash = hash_u64(hash, (uint64_t)(uint32_t)values[index]);
  }
  memcpy(&bits, &data->scaled_pivot_tolerance, sizeof(bits));
  hash = hash_u64(hash, bits);
  for (index = 0; index < (size_t)view->translation_count; ++index) {
    double parts[2];
    parts[0] = creal(view->translation_weights[index]);
    parts[1] = cimag(view->translation_weights[index]);
    memcpy(&bits, parts, sizeof(bits));
    hash = hash_u64(hash, bits);
    memcpy(&bits, parts + 1, sizeof(bits));
    hash = hash_u64(hash, bits);
  }
  return hash == 0 ? UINT64_C(1) : hash;
}
#endif

static MVMCKrylovStatus audit_world_metadata(
    const MVMCPowerLanczosGCView *view,
    MVMCClassicPfaffianCommunicator world_communicator) {
#ifdef _mpi_use
  const uint64_t local = gc_view_metadata_hash(view);
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

static MVMCGCKrylovAmplitudeLayout amplitude_layout(
    const MVMCPowerLanczosGCView *view,
    MVMCClassicPfaffianCommunicator chain_communicator) {
  const MVMCPowerLanczosClassicView *data = &view->data;
  MVMCGCKrylovAmplitudeLayout layout;
  memset(&layout, 0, sizeof(layout));
  layout.site_count = data->site_count;
  layout.domain = view->domain;
  layout.qp_total = data->qp_total;
  layout.qp_start = data->qp_start;
  layout.qp_end = data->qp_end;
  layout.translation_count = view->translation_count;
  layout.translation_weights = view->translation_weights;
  layout.spin_triplet_s0 = view->spin_triplet_s0;
  layout.scaled_pivot_tolerance = data->scaled_pivot_tolerance;
  layout.projection.site_count = (size_t)data->site_count;
  layout.projection.nproj = data->nproj;
  layout.projection.ngutzwiller_idx = data->ngutzwiller_idx;
  layout.projection.njastrow_idx = data->njastrow_idx;
  layout.projection.nspin_jastrow_idx = data->nspin_jastrow_idx;
  layout.projection.ndoublon_holon_2site_idx =
      data->ndoublon_holon_2site_idx;
  layout.projection.ndoublon_holon_4site_idx =
      data->ndoublon_holon_4site_idx;
  layout.projection.gutzwiller_idx = data->gutzwiller_idx;
  layout.projection.jastrow_idx = (const int *const *)data->jastrow_idx;
  layout.projection.spin_jastrow_idx =
      (const int *const *)data->spin_jastrow_idx;
  layout.projection.doublon_holon_2site_idx =
      (const int *const *)data->doublon_holon_2site_idx;
  layout.projection.doublon_holon_4site_idx =
      (const int *const *)data->doublon_holon_4site_idx;
  layout.projection.parameters = data->projection_parameters;
  layout.communicator = chain_communicator;
  return layout;
}

MVMCKrylovStatus mvmc_power_lanczos_gc_bridge_create(
    const MVMCPowerLanczosGCView *view,
    MVMCClassicPfaffianCommunicator world_comm,
    MVMCClassicPfaffianCommunicator chain_comm,
    MVMCPowerLanczosGCBridge **bridge) {
  MVMCPowerLanczosGCBridge *created = NULL;
  MVMCGCKrylovAmplitudeLayout layout;
  MVMCKrylovStatus status;
  int chain_rank = 0;
  int chain_size = 1;
  uint64_t minimum_generation = 0;
  uint64_t maximum_generation = 0;
  if (bridge == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *bridge = NULL;
#ifdef _mpi_use
  if (world_comm == MPI_COMM_NULL || chain_comm == MPI_COMM_NULL ||
      MPI_Comm_rank(chain_comm, &chain_rank) != MPI_SUCCESS ||
      MPI_Comm_size(chain_comm, &chain_size) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#else
  (void)world_comm;
  (void)chain_comm;
#endif
  status = validate_gc_view_shape(view, chain_rank, chain_size);
  status = synchronize_status(status, world_comm);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  status = audit_world_metadata(view, world_comm);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  created = (MVMCPowerLanczosGCBridge *)calloc(1, sizeof(*created));
  status = synchronize_status(created == NULL
                                  ? MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE
                                  : MVMC_KRYLOV_STATUS_OK,
                              world_comm);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    free(created);
    return status;
  }
  status = mvmc_power_lanczos_model_create_from_view(
      &view->data, view->domain, world_comm, &created->model_workspace);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    layout = amplitude_layout(view, chain_comm);
    status = mvmc_gc_krylov_amplitude_workspace_create(
        &layout, view->data.slater_complex, view->data.qp_weights,
        &created->amplitude_workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      created->amplitude_generation_hash =
          mvmc_gc_krylov_amplitude_generation_hash(
              created->amplitude_workspace);
    }
    status = synchronize_status(status, world_comm);
  }
#ifdef _mpi_use
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (MPI_Allreduce(&created->amplitude_generation_hash,
                     &minimum_generation, 1, MPI_UINT64_T, MPI_MIN,
                     world_comm) != MPI_SUCCESS ||
       MPI_Allreduce(&created->amplitude_generation_hash,
                     &maximum_generation, 1, MPI_UINT64_T, MPI_MAX,
                     world_comm) != MPI_SUCCESS)) {
    status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#else
  minimum_generation = created->amplitude_generation_hash;
  maximum_generation = created->amplitude_generation_hash;
#endif
  if (status == MVMC_KRYLOV_STATUS_OK &&
      (minimum_generation == 0 || minimum_generation != maximum_generation)) {
    status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    mvmc_power_lanczos_gc_bridge_destroy(created);
    return status;
  }
  created->valid = 1;
  *bridge = created;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_power_lanczos_gc_bridge_destroy(MVMCPowerLanczosGCBridge *bridge) {
  if (bridge == NULL) return;
  mvmc_gc_krylov_amplitude_workspace_destroy(bridge->amplitude_workspace);
  mvmc_classic_krylov_model_workspace_destroy(bridge->model_workspace);
  memset(bridge, 0, sizeof(*bridge));
  free(bridge);
}

const MVMCKrylovFockModel *mvmc_power_lanczos_gc_bridge_model(
    const MVMCPowerLanczosGCBridge *bridge) {
  if (bridge == NULL || !bridge->valid) return NULL;
  return mvmc_classic_krylov_model(bridge->model_workspace);
}

MVMCKrylovScaledAmplitudeCallback mvmc_power_lanczos_gc_bridge_amplitude(
    const MVMCPowerLanczosGCBridge *bridge) {
  if (bridge == NULL || !bridge->valid) return NULL;
  return mvmc_gc_krylov_scaled_amplitude;
}

void *mvmc_power_lanczos_gc_bridge_amplitude_context(
    const MVMCPowerLanczosGCBridge *bridge) {
  if (bridge == NULL || !bridge->valid) return NULL;
  return bridge->amplitude_workspace;
}

uint64_t mvmc_power_lanczos_gc_bridge_generation_hash(
    const MVMCPowerLanczosGCBridge *bridge) {
  return bridge == NULL || !bridge->valid ? 0
                                          : bridge->amplitude_generation_hash;
}
