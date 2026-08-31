#ifndef MVMC_POWER_LANCZOS_CLASSIC_BRIDGE_H
#define MVMC_POWER_LANCZOS_CLASSIC_BRIDGE_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "classic_krylov_amplitude.h"
#include "classic_krylov_model.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_CLASSIC_BRIDGE_VERSION UINT64_C(1)

typedef enum {
  MVMC_POWER_LANCZOS_CLASSIC_REAL = 0,
  MVMC_POWER_LANCZOS_CLASSIC_COMPLEX = 1
} MVMCPowerLanczosClassicArithmetic;

/*
 * A read-only, production-shaped view of synchronized classic mVMC data.
 * Counts remain signed so malformed inputs are rejected before conversion to
 * size_t or allocation.  Model rows need to be present only on world rank 0;
 * amplitude and projection data are required on every rank because amplitude
 * workspace creation is collective on each chain communicator.
 */
typedef struct {
  int site_count;
  int up_electron_count;
  int down_electron_count;
  int pure_spin;
  MVMCPowerLanczosClassicArithmetic arithmetic;

  int transfer_count;
  int *const *transfer_indices;
  const double complex *transfer_parameters;
  int coulomb_intra_count;
  const int *coulomb_intra_indices;
  const double *coulomb_intra_parameters;
  int coulomb_inter_count;
  int *const *coulomb_inter_indices;
  const double *coulomb_inter_parameters;
  int hund_count;
  int *const *hund_indices;
  const double *hund_parameters;
  int exchange_count;
  int *const *exchange_indices;
  const double *exchange_parameters;
  int pair_hopping_count;
  int inter_all_count;
  int nbody_inter_all_count;
  int nbody_g_count;

  int qp_total;
  int qp_start;
  int qp_end;
  double scaled_pivot_tolerance;
  int nproj;
  int ngutzwiller_idx;
  int njastrow_idx;
  int nspin_jastrow_idx;
  int ndoublon_holon_2site_idx;
  int ndoublon_holon_4site_idx;
  const int *gutzwiller_idx;
  int *const *jastrow_idx;
  int *const *spin_jastrow_idx;
  int *const *doublon_holon_2site_idx;
  int *const *doublon_holon_4site_idx;
  const double complex *projection_parameters;
  const double complex *qp_weights;
  const double *slater_real;
  const double complex *slater_complex;
  uint32_t unsupported_amplitude_features;
} MVMCPowerLanczosClassicView;

typedef struct MVMCPowerLanczosClassicBridge
    MVMCPowerLanczosClassicBridge;

MVMCKrylovStatus mvmc_power_lanczos_classic_bridge_create(
    const MVMCPowerLanczosClassicView *view,
    MVMCClassicPfaffianCommunicator world_communicator,
    MVMCClassicPfaffianCommunicator chain_communicator,
    MVMCPowerLanczosClassicBridge **bridge);

void mvmc_power_lanczos_classic_bridge_destroy(
    MVMCPowerLanczosClassicBridge *bridge);

const MVMCKrylovFockModel *mvmc_power_lanczos_classic_bridge_model(
    const MVMCPowerLanczosClassicBridge *bridge);

MVMCKrylovScaledAmplitudeCallback
mvmc_power_lanczos_classic_bridge_amplitude(
    const MVMCPowerLanczosClassicBridge *bridge);

void *mvmc_power_lanczos_classic_bridge_amplitude_context(
    const MVMCPowerLanczosClassicBridge *bridge);

uint64_t mvmc_power_lanczos_classic_bridge_generation_hash(
    const MVMCPowerLanczosClassicBridge *bridge);

#endif /* bounded engine */

#endif /* MVMC_POWER_LANCZOS_CLASSIC_BRIDGE_H */
