/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_CLASSIC_PFAFFIAN_STATE_H
#define MVMC_CLASSIC_PFAFFIAN_STATE_H

#include "absolute_pfaffian.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

typedef struct {
  int valid;
  uint64_t generation;
  const MVMCAbsolutePfaffianResult *components;
  MVMCProjectedAmplitudeResult local_aggregate;
} MVMCClassicPfaffianState;

typedef struct MVMCClassicPfaffianRealWorkspace
    MVMCClassicPfaffianRealWorkspace;
typedef struct MVMCClassicPfaffianComplexWorkspace
    MVMCClassicPfaffianComplexWorkspace;

/*
 * Classic non-FSZ layout: nsize is even, the first/second nsize/2 electron
 * indices are spin up/down, and Slater matrices are global-QP-major followed
 * by row-major (2*nsite)-square matrices.  The caller must provide at least
 * qp_total*(2*nsite)*(2*nsite) Slater elements and nsize electron indices.
 * QP weights are global and replicated.  qp_start == qp_end is a valid empty
 * rank-local slice.
 */
MVMCPfaffianStatus mvmc_classic_pfaffian_real_workspace_create(
    int nsite, int nsize, int qp_total, int qp_start, int qp_end,
    MVMCClassicPfaffianRealWorkspace **workspace);
void mvmc_classic_pfaffian_real_workspace_destroy(
    MVMCClassicPfaffianRealWorkspace *workspace);

MVMCPfaffianStatus mvmc_classic_pfaffian_complex_workspace_create(
    int nsite, int nsize, int qp_total, int qp_start, int qp_end,
    MVMCClassicPfaffianComplexWorkspace **workspace);
void mvmc_classic_pfaffian_complex_workspace_destroy(
    MVMCClassicPfaffianComplexWorkspace *workspace);

/*
 * Prepare evaluates every local component out of place.  It never modifies
 * accepted state or a legacy mirror and performs no allocation after create.
 */
MVMCPfaffianStatus mvmc_classic_pfaffian_real_prepare(
    MVMCClassicPfaffianRealWorkspace *workspace,
    const double *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance);

MVMCPfaffianStatus mvmc_classic_pfaffian_complex_prepare(
    MVMCClassicPfaffianComplexWorkspace *workspace,
    const double complex *slater_elm, const int *ele_idx,
    const double complex *global_weights,
    double scaled_pivot_tolerance, double residual_tolerance);

/*
 * Publish validates the complete candidate before writing.  legacy_storage
 * is the existing full allocation whose first qp_total*nsize*nsize elements
 * are InvM and whose tail qp_total elements are aliased as PfM.  Local QP
 * component k is mirrored at local slot k, matching CalculateMAll*.
 */
MVMCPfaffianStatus mvmc_classic_pfaffian_real_publish(
    MVMCClassicPfaffianRealWorkspace *workspace,
    double *legacy_storage, size_t legacy_element_count);

MVMCPfaffianStatus mvmc_classic_pfaffian_complex_publish(
    MVMCClassicPfaffianComplexWorkspace *workspace,
    double complex *legacy_storage, size_t legacy_element_count);

void mvmc_classic_pfaffian_real_discard_candidate(
    MVMCClassicPfaffianRealWorkspace *workspace);
void mvmc_classic_pfaffian_complex_discard_candidate(
    MVMCClassicPfaffianComplexWorkspace *workspace);

const MVMCClassicPfaffianState *mvmc_classic_pfaffian_real_accepted(
    const MVMCClassicPfaffianRealWorkspace *workspace);
const MVMCClassicPfaffianState *mvmc_classic_pfaffian_real_candidate(
    const MVMCClassicPfaffianRealWorkspace *workspace);
const double *mvmc_classic_pfaffian_real_accepted_inverse(
    const MVMCClassicPfaffianRealWorkspace *workspace);
const double *mvmc_classic_pfaffian_real_candidate_inverse(
    const MVMCClassicPfaffianRealWorkspace *workspace);

const MVMCClassicPfaffianState *mvmc_classic_pfaffian_complex_accepted(
    const MVMCClassicPfaffianComplexWorkspace *workspace);
const MVMCClassicPfaffianState *mvmc_classic_pfaffian_complex_candidate(
    const MVMCClassicPfaffianComplexWorkspace *workspace);
const double complex *mvmc_classic_pfaffian_complex_accepted_inverse(
    const MVMCClassicPfaffianComplexWorkspace *workspace);
const double complex *mvmc_classic_pfaffian_complex_candidate_inverse(
    const MVMCClassicPfaffianComplexWorkspace *workspace);

int mvmc_classic_pfaffian_real_failure_local_qp(
    const MVMCClassicPfaffianRealWorkspace *workspace);
int mvmc_classic_pfaffian_real_failure_global_qp(
    const MVMCClassicPfaffianRealWorkspace *workspace);
int mvmc_classic_pfaffian_complex_failure_local_qp(
    const MVMCClassicPfaffianComplexWorkspace *workspace);
int mvmc_classic_pfaffian_complex_failure_global_qp(
    const MVMCClassicPfaffianComplexWorkspace *workspace);

size_t mvmc_classic_pfaffian_real_legacy_element_count(
    const MVMCClassicPfaffianRealWorkspace *workspace);
size_t mvmc_classic_pfaffian_complex_legacy_element_count(
    const MVMCClassicPfaffianComplexWorkspace *workspace);

/*
 * Return the immutable rank-local QP ownership bound at workspace creation.
 * All output pointers are required; false is returned for invalid arguments.
 */
int mvmc_classic_pfaffian_real_qp_range(
    const MVMCClassicPfaffianRealWorkspace *workspace,
    int *qp_total, int *qp_start, int *qp_end);
int mvmc_classic_pfaffian_complex_qp_range(
    const MVMCClassicPfaffianComplexWorkspace *workspace,
    int *qp_total, int *qp_start, int *qp_end);

/* Read-only layout metadata used to bind preallocated audit workspaces. */
int mvmc_classic_pfaffian_real_layout(
    const MVMCClassicPfaffianRealWorkspace *workspace,
    int *nsize, int *qp_total, int *qp_start, int *qp_end);
int mvmc_classic_pfaffian_complex_layout(
    const MVMCClassicPfaffianComplexWorkspace *workspace,
    int *nsize, int *qp_total, int *qp_start, int *qp_end);

#endif /* MVMC_CLASSIC_PFAFFIAN_STATE_H */
