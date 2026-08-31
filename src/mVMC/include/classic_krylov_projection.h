/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_CLASSIC_KRYLOV_PROJECTION_H
#define MVMC_CLASSIC_KRYLOV_PROJECTION_H

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                     \
    defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "krylov_fock_reference.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

typedef struct {
  size_t site_count;
  int nproj;
  int ngutzwiller_idx;
  int njastrow_idx;
  int nspin_jastrow_idx;
  int ndoublon_holon_2site_idx;
  int ndoublon_holon_4site_idx;
  const int *gutzwiller_idx;
  const int *const *jastrow_idx;
  const int *const *spin_jastrow_idx;
  const int *const *doublon_holon_2site_idx;
  const int *const *doublon_holon_4site_idx;
  const double complex *parameters;
} MVMCClassicKrylovProjectionLayout;

typedef struct MVMCClassicKrylovProjectionWorkspace
    MVMCClassicKrylovProjectionWorkspace;

/*
 * The workspace owns an immutable copy of every classic projection binding.
 * Jastrow tables are site_count-square.  Each DH2/DH4 shell row contains
 * 2*site_count or 4*site_count neighbor indices respectively.  Parameters
 * follow the production MakeProjCnt order:
 * Gutzwiller, Jastrow, spin-Jastrow, 6*DH2, then 10*DH4.
 */
MVMCKrylovStatus mvmc_classic_krylov_projection_workspace_create(
    const MVMCClassicKrylovProjectionLayout *layout,
    MVMCClassicKrylovProjectionWorkspace **workspace);

void mvmc_classic_krylov_projection_workspace_destroy(
    MVMCClassicKrylovProjectionWorkspace *workspace);

size_t mvmc_classic_krylov_projection_workspace_bytes(
    const MVMCClassicKrylovProjectionWorkspace *workspace);

/* Canonical byte record used for exact replicated-binding MPI audits. */
const unsigned char *mvmc_classic_krylov_projection_binding_bytes(
    const MVMCClassicKrylovProjectionWorkspace *workspace,
    size_t *byte_count);

/*
 * Evaluation is allocation-free and transactional.  counts may be NULL when
 * count_capacity is zero.  Otherwise it must hold at least nproj entries.
 * Neither output is modified on error.  The workspace is not concurrently
 * reentrant because it owns the count scratch buffer.
 */
MVMCKrylovStatus mvmc_classic_krylov_projection_evaluate(
    MVMCClassicKrylovProjectionWorkspace *workspace,
    const int *ele_num, size_t orbital_count,
    int64_t *counts, size_t count_capacity, double *log_factor);

#endif /* reference or bounded engine */

#endif /* MVMC_CLASSIC_KRYLOV_PROJECTION_H */
