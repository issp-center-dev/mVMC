/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_CLASSIC_KRYLOV_AMPLITUDE_H
#define MVMC_CLASSIC_KRYLOV_AMPLITUDE_H

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                     \
    defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "classic_pfaffian_collective.h"
#include "krylov_fock_reference.h"
#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#include "bounded_krylov_engine.h"
#endif

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define MVMC_CLASSIC_KRYLOV_UNSUPPORTED_RBM UINT32_C(1)
#define MVMC_CLASSIC_KRYLOV_UNSUPPORTED_BACKFLOW UINT32_C(2)
#define MVMC_CLASSIC_KRYLOV_UNSUPPORTED_GENERAL_ORBITAL UINT32_C(4)
#define MVMC_CLASSIC_KRYLOV_UNSUPPORTED_BLOCK_UPDATE UINT32_C(8)

typedef struct {
  size_t site_count;
  size_t up_electron_count;
  size_t down_electron_count;
  int pure_spin;
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
  /*
   * Production-shaped index bindings.  Jastrow tables are site_count-square;
   * DH rows contain 2*site_count and 4*site_count site indices respectively.
   */
  const int *const *jastrow_idx;
  const int *const *spin_jastrow_idx;
  const int *const *doublon_holon_2site_idx;
  const int *const *doublon_holon_4site_idx;
  const double complex *projection_parameters;
  uint32_t unsupported_features;
  MVMCClassicPfaffianCommunicator communicator;
} MVMCClassicKrylovAmplitudeLayout;

typedef struct MVMCClassicKrylovRealAmplitudeWorkspace
    MVMCClassicKrylovRealAmplitudeWorkspace;
typedef struct MVMCClassicKrylovComplexAmplitudeWorkspace
    MVMCClassicKrylovComplexAmplitudeWorkspace;

/*
 * Creation is collective on layout->communicator.  It verifies replicated
 * Slater data, QP weights, and every classic projection binding exactly,
 * then keeps an immutable private copy.  Evaluation performs no allocation
 * and never writes production PfM/InvM, accepted sampler state, projection
 * counts, caller buffers, or the input configuration.  Projected values are
 * gathered and summed in global QP index order, independent of MPI rank
 * count and partitioning.
 */
MVMCKrylovStatus mvmc_classic_krylov_real_amplitude_workspace_create(
    const MVMCClassicKrylovAmplitudeLayout *layout,
    const double *slater_elm, const double complex *global_weights,
    MVMCClassicKrylovRealAmplitudeWorkspace **workspace);

void mvmc_classic_krylov_real_amplitude_workspace_destroy(
    MVMCClassicKrylovRealAmplitudeWorkspace *workspace);
size_t mvmc_classic_krylov_real_amplitude_workspace_bytes(
    const MVMCClassicKrylovRealAmplitudeWorkspace *workspace);

MVMCKrylovStatus mvmc_classic_krylov_complex_amplitude_workspace_create(
    const MVMCClassicKrylovAmplitudeLayout *layout,
    const double complex *slater_elm,
    const double complex *global_weights,
    MVMCClassicKrylovComplexAmplitudeWorkspace **workspace);

void mvmc_classic_krylov_complex_amplitude_workspace_destroy(
    MVMCClassicKrylovComplexAmplitudeWorkspace *workspace);
size_t mvmc_classic_krylov_complex_amplitude_workspace_bytes(
    const MVMCClassicKrylovComplexAmplitudeWorkspace *workspace);

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)
MVMCKrylovStatus mvmc_classic_krylov_real_amplitude(
    const uint64_t *configuration_words, size_t word_count,
    void *context, MVMCKrylovAmplitudeResult *result);

MVMCKrylovStatus mvmc_classic_krylov_complex_amplitude(
    const uint64_t *configuration_words, size_t word_count,
    void *context, MVMCKrylovAmplitudeResult *result);
#endif

uint64_t mvmc_classic_krylov_real_amplitude_generation_hash(
    const MVMCClassicKrylovRealAmplitudeWorkspace *workspace);
uint64_t mvmc_classic_krylov_complex_amplitude_generation_hash(
    const MVMCClassicKrylovComplexAmplitudeWorkspace *workspace);

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
MVMCKrylovStatus mvmc_classic_krylov_real_scaled_amplitude(
    const uint64_t *configuration_words, size_t word_count,
    void *context, MVMCKrylovScaledAmplitudeResult *result);

MVMCKrylovStatus mvmc_classic_krylov_complex_scaled_amplitude(
    const uint64_t *configuration_words, size_t word_count,
    void *context, MVMCKrylovScaledAmplitudeResult *result);
#endif

#endif /* reference or bounded engine */

#endif /* MVMC_CLASSIC_KRYLOV_AMPLITUDE_H */
