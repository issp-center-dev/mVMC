#ifndef MVMC_POWER_LANCZOS_GC_BRIDGE_H
#define MVMC_POWER_LANCZOS_GC_BRIDGE_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "gc_krylov_amplitude.h"
#include "power_lanczos_classic_bridge.h"

#include <complex.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_GC_BRIDGE_VERSION UINT64_C(1)

/*
 * Production-shaped grand-canonical view.  data.up/down_electron_count and
 * data.pure_spin are 0, data.arithmetic is COMPLEX and data.slater_complex is
 * the full QP SlaterElm.  General pairing is expressed by the domain
 * (GC_EVEN); the canonical UNSUPPORTED_GENERAL_ORBITAL bit is not used and
 * any unsupported amplitude feature is rejected.  translation_weights are the
 * raw input translation weights (ParaQPTrans) used by the boundary rule.
 * Model rows are read on world rank 0 only; amplitude data is replicated on
 * every rank of each chain communicator.
 */
typedef struct {
  MVMCPowerLanczosClassicView data;
  MVMCKrylovConfigurationDomain domain;
  int translation_count;
  const double complex *translation_weights;
} MVMCPowerLanczosGCView;

typedef struct MVMCPowerLanczosGCBridge MVMCPowerLanczosGCBridge;

MVMCKrylovStatus mvmc_power_lanczos_gc_bridge_create(
    const MVMCPowerLanczosGCView *view,
    MVMCClassicPfaffianCommunicator world_comm,
    MVMCClassicPfaffianCommunicator chain_comm,
    MVMCPowerLanczosGCBridge **bridge);
void mvmc_power_lanczos_gc_bridge_destroy(MVMCPowerLanczosGCBridge *bridge);
const MVMCKrylovFockModel *mvmc_power_lanczos_gc_bridge_model(
    const MVMCPowerLanczosGCBridge *bridge);
MVMCKrylovScaledAmplitudeCallback mvmc_power_lanczos_gc_bridge_amplitude(
    const MVMCPowerLanczosGCBridge *bridge);
void *mvmc_power_lanczos_gc_bridge_amplitude_context(
    const MVMCPowerLanczosGCBridge *bridge);
uint64_t mvmc_power_lanczos_gc_bridge_generation_hash(
    const MVMCPowerLanczosGCBridge *bridge);

#endif /* MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE */

#endif /* MVMC_POWER_LANCZOS_GC_BRIDGE_H */
