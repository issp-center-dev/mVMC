/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_GC_KRYLOV_AMPLITUDE_H
#define MVMC_GC_KRYLOV_AMPLITUDE_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "bounded_krylov_engine.h"
#include "classic_krylov_projection.h"
#include "classic_pfaffian_collective.h"
#include "krylov_fock_reference.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

/*
 * Private absolute amplitude of the grand-canonical wave function for every
 * configuration of a GC domain.  For the ascending occupied orbitals x of a
 * configuration with N electrons the QP component is the production
 * GCComponentValue (-1)^(N/2) Pf(-Slater[x,x]); N = 0 is the empty Pfaffian
 * 1.  Components are summed in global QP order, multiplied by the Re(Proj)
 * correlation factor, and the Phase 3b boundary rule (N = 0 or 2L with an
 * exactly zero sum of the raw translation weights) gives EXACT_ZERO after the
 * components were checked.  Nonfinite components are never hidden.
 *
 * Slater (qp_total * (2L)^2, row-major as SlaterElm), QP weights, translation
 * weights and the projection binding are replicated on every rank of the
 * communicator and copied privately; the production GC sampler state is
 * neither read nor written.
 */
typedef struct {
  int site_count;
  MVMCKrylovConfigurationDomain domain;
  int qp_total, qp_start, qp_end;
  int translation_count;
  const double complex *translation_weights;
  /* Proven by exact comparison of the unrotated parameter values for this
   * generation: active S=0 projection and F^T=-F. Balanced domain only. */
  int spin_triplet_s0;
  double scaled_pivot_tolerance;
  MVMCClassicKrylovProjectionLayout projection;
  MVMCClassicPfaffianCommunicator communicator;
} MVMCGCKrylovAmplitudeLayout;

typedef struct MVMCGCKrylovAmplitudeWorkspace MVMCGCKrylovAmplitudeWorkspace;

/* Collective on layout->communicator. */
MVMCKrylovStatus mvmc_gc_krylov_amplitude_workspace_create(
    const MVMCGCKrylovAmplitudeLayout *layout,
    const double complex *slater, const double complex *qp_weights,
    MVMCGCKrylovAmplitudeWorkspace **workspace);
void mvmc_gc_krylov_amplitude_workspace_destroy(
    MVMCGCKrylovAmplitudeWorkspace *workspace);
size_t mvmc_gc_krylov_amplitude_workspace_bytes(
    const MVMCGCKrylovAmplitudeWorkspace *workspace);
uint64_t mvmc_gc_krylov_amplitude_generation_hash(
    const MVMCGCKrylovAmplitudeWorkspace *workspace);

/* Scaled amplitude callback; collective on the workspace communicator. */
MVMCKrylovStatus mvmc_gc_krylov_scaled_amplitude(
    const uint64_t *words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result);

#if defined(MVMC_ENABLE_FAULT_INJECTION)
typedef enum {
  MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONE = 0,
  /* Workspace creation reports an allocation failure on the target rank. */
  MVMC_GC_KRYLOV_AMPLITUDE_FAULT_ALLOCATION,
  /* Every evaluation reports a nonfinite component on the target rank. */
  MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONFINITE
} MVMCGCKrylovAmplitudeFault;

/* Testing only; world_rank selects the MPI_COMM_WORLD rank (0 serial). */
void mvmc_gc_krylov_amplitude_testing_set_fault(
    MVMCGCKrylovAmplitudeFault fault, int world_rank);
#endif

#endif /* MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE */

#endif /* MVMC_GC_KRYLOV_AMPLITUDE_H */
