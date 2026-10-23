/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_KRYLOV_FINAL_STATE_H
#define MVMC_KRYLOV_FINAL_STATE_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_P5_CORE) ||                             \
    defined(MVMC_ENABLE_POWER_LANCZOS_P5_TESTING)

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "P5 final-state core requires the bounded Krylov engine"
#endif

#include "bounded_krylov_engine.h"
#include "krylov_fock_proposal.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define MVMC_KRYLOV_FINAL_STATE_POLICY_VERSION UINT64_C(2)

typedef enum {
  MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT = 1,
  MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC = 2,
  MVMC_KRYLOV_FINAL_COEFFICIENT_PERTURBED_TESTING = 3,
  MVMC_KRYLOV_FINAL_COEFFICIENT_PRODUCTION_NOISY = 4
} MVMCKrylovFinalCoefficientSource;

typedef struct {
  int order;
  MVMCKrylovFinalCoefficientSource source;
  uint64_t provenance_hash;
  uint64_t policy_hash;
  double complex coefficient[MVMC_KRYLOV_MAX_ORDER + 1];
} MVMCKrylovFinalStatePolicy;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int order;
  uint64_t policy_hash;
  int sampleable;
  int local_energy_available;
  int coefficient_zero_count;
  int basis_finite_count;
  int basis_exact_zero_count;
  int basis_numeric_zero_count;
  MVMCScaledComplex amplitude;
  MVMCScaledComplex hamiltonian_amplitude;
  MVMCScaledComplex local_energy;
  double log_weight;
} MVMCKrylovFinalStateEvaluation;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int accepted;
  int exact_zero_proposal;
  double uniform;
  double log_target_ratio;
  double log_proposal_ratio;
  double log_acceptance_ratio;
} MVMCKrylovFinalStateAcceptance;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  size_t word_count;
  uint64_t policy_hash;
  uint64_t configuration_hash;
  uint64_t amplitude_generation_hash;
  uint64_t integrity_hash;
  uint64_t accepted_generation;
  MVMCKrylovBoundedResult krylov;
  MVMCKrylovFinalStateEvaluation final_state;
} MVMCKrylovFinalStateSnapshot;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int accepted;
  uint64_t accepted_generation;
  MVMCKrylovBoundedResult proposal_krylov;
  MVMCKrylovFinalStateEvaluation proposal_final_state;
  MVMCKrylovFinalStateAcceptance acceptance;
} MVMCKrylovFinalStateStepResult;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  int configuration_changed;
  size_t neighbor_count;
  size_t selected_neighbor_index;
  double log_proposal_ratio;
  MVMCKrylovFinalStateStepResult step;
} MVMCKrylovFinalStateProposalStepResult;

MVMCKrylovStatus mvmc_krylov_final_state_policy_create(
    int order, MVMCKrylovFinalCoefficientSource source,
    uint64_t provenance_hash, const double complex *coefficient,
    size_t coefficient_count, MVMCKrylovFinalStatePolicy *policy);

MVMCKrylovStatus mvmc_krylov_final_state_policy_create_scaled_basis(
    int order, MVMCKrylovFinalCoefficientSource source,
    uint64_t provenance_hash, const double complex *scaled_coefficient,
    const double *log_basis_scale, size_t coefficient_count,
    MVMCKrylovFinalStatePolicy *policy);

uint64_t mvmc_krylov_final_state_policy_hash(
    const MVMCKrylovFinalStatePolicy *policy);

uint64_t mvmc_krylov_final_state_snapshot_integrity_hash(
    const MVMCKrylovFinalStateSnapshot *snapshot);

/* Construct a checked snapshot from an independently refined evaluation.
 * No configuration or existing snapshot is changed on failure. */
MVMCKrylovStatus mvmc_krylov_final_state_snapshot_from_evaluation(
    const MVMCKrylovFinalStatePolicy *policy, const uint64_t *words,
    size_t word_count, const MVMCKrylovBoundedResult *krylov,
    const MVMCKrylovFinalStateEvaluation *evaluation,
    MVMCKrylovFinalStateSnapshot *snapshot);

MVMCKrylovStatus mvmc_krylov_final_state_evaluate(
    const MVMCKrylovFinalStatePolicy *policy,
    const MVMCKrylovBoundedResult *krylov,
    MVMCKrylovFinalStateEvaluation *result);

MVMCScaledComplexExportStatus mvmc_krylov_final_state_local_energy_export(
    const MVMCKrylovFinalStateEvaluation *evaluation,
    double complex *local_energy);

MVMCKrylovStatus mvmc_krylov_final_state_acceptance(
    const MVMCKrylovFinalStateEvaluation *current,
    const MVMCKrylovFinalStateEvaluation *proposal,
    double log_proposal_ratio, double uniform,
    MVMCKrylovFinalStateAcceptance *result);

MVMCKrylovStatus mvmc_krylov_final_state_sampler_initialize(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy,
    const uint64_t *current_words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovFinalStateSnapshot *snapshot);

MVMCKrylovStatus mvmc_krylov_final_state_sampler_step(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy,
    uint64_t *current_words, size_t current_word_count,
    MVMCKrylovFinalStateSnapshot *current,
    const uint64_t *proposal_words, size_t proposal_word_count,
    double log_proposal_ratio, double uniform,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCKrylovFinalStateStepResult *result);

/*
 * Why an initialization candidate may be drawn again.  Only a completed
 * evaluation whose Phi is EXACT_ZERO or NUMERIC_ZERO, or whose finite
 * H Phi / Phi cannot be formed safely (local energy unavailable), is
 * recoverable; a failing amplitude callback, nonfinite values and every
 * other status are NO_RETRY.
 */
typedef enum {
  MVMC_KRYLOV_BOOTSTRAP_NO_RETRY = 0,
  MVMC_KRYLOV_BOOTSTRAP_EXACT_ZERO,
  MVMC_KRYLOV_BOOTSTRAP_NUMERIC_ZERO,
  MVMC_KRYLOV_BOOTSTRAP_LOCAL_ENERGY
} MVMCKrylovBootstrapReason;

/*
 * Diagnostic copy of the last final-state evaluation.  evaluation_available
 * means that evaluation.amplitude (Phi) and evaluation.hamiltonian_amplitude
 * (H Phi) are meaningful; it never makes evaluation.valid or sampleable true
 * and never validates a failed snapshot or step result.
 */
typedef struct {
  int evaluation_available;
  MVMCKrylovBootstrapReason retry_reason;
  MVMCKrylovStatus status;
  MVMCKrylovFinalStateEvaluation evaluation;
} MVMCKrylovFinalDiagnostic;

/*
 * Diagnosed forms of initialize/step.  diagnostic may be NULL; the legacy
 * entries are these with NULL and keep their snapshot/result contracts.
 * A NUMERIC_ZERO proposal or an unavailable proposal local energy is still
 * an error of the step; EXACT_ZERO proposals are still rejected.
 */
MVMCKrylovStatus mvmc_krylov_final_state_sampler_initialize_diagnosed(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy,
    const uint64_t *words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *context,
    MVMCKrylovFinalStateSnapshot *snapshot,
    MVMCKrylovFinalDiagnostic *diagnostic);

MVMCKrylovStatus mvmc_krylov_final_state_sampler_step_diagnosed(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy,
    uint64_t *words, size_t word_count,
    MVMCKrylovFinalStateSnapshot *current,
    const uint64_t *proposal_words, size_t proposal_word_count,
    double log_proposal_ratio, double uniform,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *context,
    MVMCKrylovFinalStateStepResult *result,
    MVMCKrylovFinalDiagnostic *diagnostic);

MVMCKrylovStatus mvmc_krylov_final_state_sampler_step_selected_neighbor(
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFinalStatePolicy *policy,
    const MVMCKrylovFockModel *proposal_model,
    uint64_t *current_words, size_t current_word_count,
    MVMCKrylovFinalStateSnapshot *current,
    size_t selected_neighbor_index, double uniform,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words, size_t proposal_word_count,
    MVMCKrylovFinalStateProposalStepResult *result);

#endif /* P5 core or testing */

#endif /* MVMC_KRYLOV_FINAL_STATE_H */
