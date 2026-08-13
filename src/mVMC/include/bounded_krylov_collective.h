/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_BOUNDED_KRYLOV_COLLECTIVE_H
#define MVMC_BOUNDED_KRYLOV_COLLECTIVE_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "bounded_krylov_engine.h"

#ifdef _mpi_use
#include <mpi.h>
typedef MPI_Comm MVMCKrylovBoundedCommunicator;
#else
typedef int MVMCKrylovBoundedCommunicator;
#endif

typedef struct MVMCKrylovBoundedCollectiveWorkspace
    MVMCKrylovBoundedCollectiveWorkspace;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int failure_rank;
  uint64_t plan_hash;
} MVMCKrylovBoundedCollectiveResult;

/*
 * Creation is collective.  The communicator is duplicated and all buffers
 * for QP gathers and rank-ordered accumulator merges are allocated here.
 * Calls using the workspace allocate no memory.  Every later collective call
 * requires a valid workspace on every rank and the same call order.
 */
MVMCKrylovStatus mvmc_bounded_krylov_collective_workspace_create(
    MVMCKrylovBoundedCommunicator communicator,
    size_t max_qp_components, size_t max_accumulator_count,
    MVMCKrylovBoundedCollectiveWorkspace **workspace);
void mvmc_bounded_krylov_collective_workspace_destroy(
    MVMCKrylovBoundedCollectiveWorkspace *workspace);
size_t mvmc_bounded_krylov_collective_workspace_bytes(
    const MVMCKrylovBoundedCollectiveWorkspace *workspace);

/*
 * Synchronize plan identity, allocation readiness, and a local evaluation or
 * callback status.  The highest enum severity wins; ties use the lowest rank.
 * A plan mismatch is INTERNAL_INVARIANT_FAILURE.  All ranks receive the same
 * status and failure rank.
 */
MVMCKrylovStatus mvmc_bounded_krylov_collective_synchronize(
    MVMCKrylovBoundedCollectiveWorkspace *workspace,
    MVMCKrylovStatus local_status, uint64_t local_plan_hash,
    int local_allocation_ok, MVMCKrylovBoundedCollectiveResult *result);

MVMCKrylovStatus mvmc_bounded_krylov_collective_max_u64(
    MVMCKrylovBoundedCollectiveWorkspace *workspace,
    uint64_t local_value, uint64_t *global_maximum);

/*
 * Gather a rank-ordered contiguous partition [0, qp_total) into global QP
 * order.  Empty local slices are valid.  global_components is committed only
 * after every rank, range, and scaled component has passed validation.
 */
MVMCKrylovStatus mvmc_bounded_krylov_collective_gather_scaled_components(
    MVMCKrylovBoundedCollectiveWorkspace *workspace,
    MVMCKrylovStatus local_status,
    const MVMCAbsolutePfaffianScaledValueResult *local_components,
    size_t local_component_count, int qp_total, int qp_start, int qp_end,
    MVMCAbsolutePfaffianScaledValueResult *global_components,
    MVMCKrylovBoundedCollectiveResult *result);

/*
 * Allocation-free global-QP projection for terminal-amplitude callbacks.
 * The weight array must be replicated; an exact stable-hash audit rejects a
 * mismatch before gathering.  The returned amplitude is committed only on
 * collective success.
 */
MVMCKrylovStatus
mvmc_bounded_krylov_collective_gather_projected_amplitude(
    MVMCKrylovBoundedCollectiveWorkspace *workspace,
    MVMCKrylovStatus local_status,
    const MVMCAbsolutePfaffianScaledValueResult *local_components,
    size_t local_component_count, const double complex *global_weights,
    int qp_total, int qp_start, int qp_end,
    MVMCKrylovScaledAmplitudeResult *amplitude,
    MVMCKrylovBoundedCollectiveResult *result);

/*
 * Merge one equal-length scaled accumulator array per rank.  For every array
 * index, operands are summed in increasing communicator-rank order.  Output
 * is transactional and identical on all ranks.
 */
MVMCKrylovStatus mvmc_bounded_krylov_collective_merge_scaled_accumulators(
    MVMCKrylovBoundedCollectiveWorkspace *workspace,
    MVMCKrylovStatus local_status,
    const MVMCScaledComplex *local_values, size_t value_count,
    MVMCScaledComplex *global_values,
    MVMCKrylovBoundedCollectiveResult *result);

#endif /* MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE */

#endif /* MVMC_BOUNDED_KRYLOV_COLLECTIVE_H */
