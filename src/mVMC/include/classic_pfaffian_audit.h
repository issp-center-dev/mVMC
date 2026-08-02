/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_CLASSIC_PFAFFIAN_AUDIT_H
#define MVMC_CLASSIC_PFAFFIAN_AUDIT_H

#include "classic_pfaffian_collective.h"

#if defined(MVMC_ENABLE_PFAFFIAN_STATE_AUDIT)

typedef enum {
  MVMC_CLASSIC_AUDIT_OFF = 0,
  MVMC_CLASSIC_AUDIT_GUARDED,
  MVMC_CLASSIC_AUDIT_ALWAYS
} MVMCClassicPfaffianAuditMode;

typedef struct {
  int valid;
  MVMCPfaffianStatus status;
  MVMCClassicPfaffianAuditMode requested_mode;
  MVMCClassicPfaffianAuditMode effective_mode;
  int executed;
  int fallback;
  int mismatch;
  int failure_rank;
  int failure_local_qp;
  int failure_global_qp;
  MVMCProjectedAmplitudeResult observed_aggregate;
} MVMCClassicPfaffianAuditReport;

typedef struct {
  int valid;
  size_t component_count;
  size_t matrix_element_count;
  MVMCAbsolutePfaffianResult *components;
  double *inverse;
} MVMCClassicPfaffianRealObservation;

typedef struct {
  int valid;
  size_t component_count;
  size_t matrix_element_count;
  MVMCAbsolutePfaffianResult *components;
  double complex *inverse;
} MVMCClassicPfaffianComplexObservation;

typedef MVMCPfaffianStatus (*MVMCClassicPfaffianRealObserver)(
    void *context, const MVMCClassicPfaffianState *absolute_candidate,
    const double *absolute_inverse,
    MVMCClassicPfaffianRealObservation *observation);

typedef MVMCPfaffianStatus (*MVMCClassicPfaffianComplexObserver)(
    void *context, const MVMCClassicPfaffianState *absolute_candidate,
    const double complex *absolute_inverse,
    MVMCClassicPfaffianComplexObservation *observation);

typedef struct MVMCClassicPfaffianRealAuditWorkspace
    MVMCClassicPfaffianRealAuditWorkspace;
typedef struct MVMCClassicPfaffianComplexAuditWorkspace
    MVMCClassicPfaffianComplexAuditWorkspace;

/* Creation is allowed to allocate; every later audit call is allocation-free. */
MVMCPfaffianStatus mvmc_classic_pfaffian_real_audit_workspace_create(
    const MVMCClassicPfaffianRealWorkspace *state_workspace,
    MVMCClassicPfaffianAuditMode mode,
    double absolute_tolerance, double relative_tolerance,
    MVMCClassicPfaffianRealAuditWorkspace **workspace);
void mvmc_classic_pfaffian_real_audit_workspace_destroy(
    MVMCClassicPfaffianRealAuditWorkspace *workspace);

MVMCPfaffianStatus mvmc_classic_pfaffian_complex_audit_workspace_create(
    const MVMCClassicPfaffianComplexWorkspace *state_workspace,
    MVMCClassicPfaffianAuditMode mode,
    double absolute_tolerance, double relative_tolerance,
    MVMCClassicPfaffianComplexAuditWorkspace **workspace);
void mvmc_classic_pfaffian_complex_audit_workspace_destroy(
    MVMCClassicPfaffianComplexAuditWorkspace *workspace);

MVMCClassicPfaffianAuditMode mvmc_classic_pfaffian_real_audit_mode(
    const MVMCClassicPfaffianRealAuditWorkspace *workspace);
MVMCClassicPfaffianAuditMode mvmc_classic_pfaffian_complex_audit_mode(
    const MVMCClassicPfaffianComplexAuditWorkspace *workspace);

int mvmc_classic_pfaffian_real_audit_compatible(
    const MVMCClassicPfaffianRealAuditWorkspace *audit_workspace,
    const MVMCClassicPfaffianRealWorkspace *state_workspace);
int mvmc_classic_pfaffian_complex_audit_compatible(
    const MVMCClassicPfaffianComplexAuditWorkspace *audit_workspace,
    const MVMCClassicPfaffianComplexWorkspace *state_workspace);

/*
 * Observe and compare one proposal.  In GUARDED mode every rank first checks
 * the accepted state; one ineligible rank selects absolute-only fallback for
 * all ranks.  The observer writes only to the supplied preallocated output.
 * These are collective operations: every rank in collective_workspace must
 * call with a valid, layout-compatible workspace, state, weight array, and
 * observer.  Sampler entry points enforce that precondition collectively.
 */
MVMCPfaffianStatus mvmc_classic_pfaffian_real_audit_candidate(
    MVMCClassicPfaffianRealAuditWorkspace *audit_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const MVMCClassicPfaffianState *accepted_state,
    const MVMCClassicPfaffianState *absolute_candidate,
    const double *absolute_candidate_inverse,
    const MVMCProjectedAmplitudeResult *absolute_global_aggregate,
    const double complex *global_weights,
    MVMCClassicPfaffianRealObserver observer, void *observer_context,
    MVMCClassicPfaffianAuditReport *report);

MVMCPfaffianStatus mvmc_classic_pfaffian_complex_audit_candidate(
    MVMCClassicPfaffianComplexAuditWorkspace *audit_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const MVMCClassicPfaffianState *accepted_state,
    const MVMCClassicPfaffianState *absolute_candidate,
    const double complex *absolute_candidate_inverse,
    const MVMCProjectedAmplitudeResult *absolute_global_aggregate,
    const double complex *global_weights,
    MVMCClassicPfaffianComplexObserver observer, void *observer_context,
    MVMCClassicPfaffianAuditReport *report);

/*
 * Force a state comparison, used by periodic rebuild regardless of mode.
 * These are collective operations with the same all-rank validity contract.
 * A forced comparison reports effective_mode=ALWAYS even if requested_mode=OFF.
 */
MVMCPfaffianStatus mvmc_classic_pfaffian_real_audit_compare(
    MVMCClassicPfaffianRealAuditWorkspace *audit_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const MVMCClassicPfaffianState *reference_state,
    const double *reference_inverse,
    const MVMCProjectedAmplitudeResult *reference_global_aggregate,
    const MVMCClassicPfaffianState *observed_state,
    const double *observed_inverse,
    MVMCClassicPfaffianAuditReport *report);

MVMCPfaffianStatus mvmc_classic_pfaffian_complex_audit_compare(
    MVMCClassicPfaffianComplexAuditWorkspace *audit_workspace,
    MVMCClassicPfaffianCollectiveWorkspace *collective_workspace,
    const MVMCClassicPfaffianState *reference_state,
    const double complex *reference_inverse,
    const MVMCProjectedAmplitudeResult *reference_global_aggregate,
    const MVMCClassicPfaffianState *observed_state,
    const double complex *observed_inverse,
    MVMCClassicPfaffianAuditReport *report);

#endif /* MVMC_ENABLE_PFAFFIAN_STATE_AUDIT */

#endif /* MVMC_CLASSIC_PFAFFIAN_AUDIT_H */
