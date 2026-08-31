/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_CLASSIC_PFAFFIAN_COLLECTIVE_H
#define MVMC_CLASSIC_PFAFFIAN_COLLECTIVE_H

#include "classic_pfaffian_state.h"

#include <stddef.h>
#include <stdint.h>

#ifdef _mpi_use
#include <mpi.h>
typedef MPI_Comm MVMCClassicPfaffianCommunicator;
#else
typedef int MVMCClassicPfaffianCommunicator;
#endif

typedef struct MVMCClassicPfaffianCollectiveWorkspace
    MVMCClassicPfaffianCollectiveWorkspace;

typedef struct {
  int valid;
  MVMCPfaffianStatus status;
  int failure_rank;
  int failure_local_qp;
  int failure_global_qp;
  MVMCProjectedAmplitudeResult aggregate;
} MVMCClassicPfaffianCollectiveResult;

typedef struct {
  int valid;
  MVMCPfaffianStatus status;
  int accepted;
  double log_acceptance_ratio;
} MVMCClassicPfaffianMetropolisResult;

/*
 * Creation and destruction are collective on communicator in MPI builds.
 * The workspace duplicates the passed communicator and preallocates its
 * ownership audit buffer; aggregate/prepare calls allocate no memory.
 * MPI communicators must retain the default MPI_ERRORS_ARE_FATAL handler.
 * Primitive MPI failures cannot be recovered collectively after only a
 * subset of ranks returns an error.
 */
MVMCPfaffianStatus mvmc_classic_pfaffian_collective_workspace_create(
    MVMCClassicPfaffianCommunicator communicator,
    MVMCClassicPfaffianCollectiveWorkspace **workspace);
void mvmc_classic_pfaffian_collective_workspace_destroy(
    MVMCClassicPfaffianCollectiveWorkspace *workspace);
#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                     \
    defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
size_t mvmc_classic_pfaffian_collective_workspace_bytes(
    const MVMCClassicPfaffianCollectiveWorkspace *workspace);
#endif

/*
 * All ranks must call in the same order with a workspace from the same
 * communicator.  Local QP ranges must form one rank-ordered contiguous
 * partition [0, qp_total); empty ranges are valid.  A local failure still
 * participates in every collective and is reported identically on all
 * ranks.  Severity is fixed as factorization > allocation > unsupported-FP
 * > invalid-input > OK, with the lowest communicator rank breaking ties.
 * Failure leaves result->aggregate zeroed and invalid.
 * A null result on one rank is propagated as invalid input without allowing
 * that rank to skip the collective sequence.
 * The workspace itself must be valid on every rank; without its communicator
 * a rank cannot participate in failure propagation.
 */
MVMCPfaffianStatus mvmc_classic_pfaffian_collective_aggregate(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    MVMCPfaffianStatus local_status,
    const MVMCProjectedAmplitudeResult *local_aggregate,
    int qp_total, int qp_start, int qp_end,
    int failure_local_qp, int failure_global_qp,
    MVMCClassicPfaffianCollectiveResult *result);

/*
 * Collective absolute proposals.  Each rank first evaluates its bound local
 * QP slice, then all ranks enter the same status/aggregate collective.  If
 * any rank fails, every local candidate is discarded and accepted state is
 * unchanged.  Successful calls leave local candidates ready for publish.
 */
MVMCPfaffianStatus mvmc_classic_pfaffian_real_prepare_collective(
    MVMCClassicPfaffianRealWorkspace *state_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const double *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCClassicPfaffianCollectiveResult *result);

MVMCPfaffianStatus mvmc_classic_pfaffian_complex_prepare_collective(
    MVMCClassicPfaffianComplexWorkspace *state_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const double complex *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance,
    MVMCClassicPfaffianCollectiveResult *result);

/*
 * Make one communicator-wide Metropolis decision from authoritative global
 * totals.  Inputs must be identical on every rank.  A zero proposal total is
 * a valid rejection; a zero current total is invalid.  Infinite log factors
 * are supported, NaN is not.  The caller owns the single uniform draw and
 * must pass a value in [0,1).
 */
MVMCPfaffianStatus mvmc_classic_pfaffian_collective_metropolis(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    double complex current_total, double complex proposal_total,
    double log_non_pfaffian_ratio, double uniform,
    MVMCClassicPfaffianMetropolisResult *result);

/* Collective guard for a branch that must precede proposal evaluation. */
MVMCPfaffianStatus mvmc_classic_pfaffian_collective_preflight(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    int local_valid, int branch_kind);

/*
 * Return the communicator-wide logical AND without treating false as an
 * error.  This is used to select an optional branch collectively.  local_true
 * must be exactly 0 or 1 on every rank; invalid values are propagated only
 * after every rank has participated in the reduction.
 */
MVMCPfaffianStatus mvmc_classic_pfaffian_collective_all_true(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    int local_true, int *all_true);

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                     \
    defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
/* Creation-time exact replicated-data audit.  This routine may allocate. */
MVMCPfaffianStatus mvmc_classic_pfaffian_collective_all_equal_bytes(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    const void *local_data, size_t byte_count, int *all_equal);

/* Allocation-free exact audit for evaluation-time configuration words. */
MVMCPfaffianStatus mvmc_classic_pfaffian_collective_all_equal_u64(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    const uint64_t *local_values, size_t value_count, int *all_equal);

MVMCPfaffianStatus mvmc_classic_pfaffian_collective_sum_u64(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    uint64_t local_value, uint64_t *global_sum);

MVMCPfaffianStatus mvmc_classic_pfaffian_collective_max_int(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    int local_value, int *global_max);

/*
 * Gather a validated contiguous QP partition into global QP index order.
 * No allocation is performed; global_components has qp_total elements on
 * every rank and local_components has qp_end-qp_start elements.
 */
MVMCPfaffianStatus
mvmc_classic_pfaffian_collective_gather_value_components(
    MVMCClassicPfaffianCollectiveWorkspace *workspace,
    const MVMCAbsolutePfaffianValueResult *local_components,
    size_t local_component_count, int qp_total, int qp_start, int qp_end,
    MVMCAbsolutePfaffianValueResult *global_components);
#endif /* reference or bounded engine */

#endif /* MVMC_CLASSIC_PFAFFIAN_COLLECTIVE_H */
