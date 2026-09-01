/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_POWER_LANCZOS_INDEPENDENT_H
#define MVMC_POWER_LANCZOS_INDEPENDENT_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#include "bounded_krylov_engine.h"
#include "power_lanczos_classic_bridge.h"

#define MVMC_POWER_LANCZOS_INDEPENDENT_VERSION UINT64_C(2)

typedef struct MVMCPowerLanczosIndependentModel
    MVMCPowerLanczosIndependentModel;
typedef struct MVMCPowerLanczosIndependentStreamingModel
    MVMCPowerLanczosIndependentStreamingModel;

typedef struct {
    MVMCScaledComplex q[2];
    MVMCScaledComplex r[2];
    uint64_t physical_row_terms;
    uint64_t amplitude_calls;
} MVMCPowerLanczosIndependentSnapshot;

typedef struct {
    const MVMCPowerLanczosClassicView *classic_view;
    MVMCClassicPfaffianCommunicator world_communicator;
    MVMCClassicPfaffianCommunicator chain_communicator;
    size_t physical_transfer_count;
    int *const *physical_transfer_indices;
    const double complex *physical_transfer_parameters;
    size_t physical_inter_all_count;
    int *const *physical_inter_all_indices;
    const double complex *physical_inter_all_parameters;
    size_t hprime_transfer_count;
    int *const *hprime_transfer_indices;
    const double complex *hprime_transfer_parameters;
    size_t hprime_inter_all_count;
    int *const *hprime_inter_all_indices;
    const double complex *hprime_inter_all_parameters;
    uint64_t seed;
    size_t warm_up;
    size_t sample_count;
    size_t interval;
} MVMCPowerLanczosIndependentInput;

typedef struct {
    int valid;
    MVMCKrylovStatus status;
    uint64_t version;
    uint64_t coefficient_samples;
    uint64_t final_samples;
    uint64_t physical_row_terms;
    uint64_t sampling_chains;
    size_t block_count;
    int retained_rank;
    double condition_estimate;
    double gevp_residual;
    double energy;
    double energy_standard_error;
    double variance;
    double variance_standard_error;
    double final_energy_imaginary;
    double variance_imaginary;
    double energy_tau_int;
    double effective_sample_count;
} MVMCPowerLanczosIndependentResult;

/* Lower Trans (-ParaTransfer) and written-order CACA InterAll rows. */
MVMCKrylovStatus mvmc_power_lanczos_independent_model_create(
    int site_count, int up_electron_count, int down_electron_count,
    size_t transfer_count, int *const *transfer_indices,
    const double complex *transfer_parameters, size_t inter_all_count,
    int *const *inter_all_indices, const double complex *inter_all_parameters,
    MVMCPowerLanczosIndependentModel **result);

void mvmc_power_lanczos_independent_model_destroy(
    MVMCPowerLanczosIndependentModel *model);

const MVMCKrylovFockModel *mvmc_power_lanczos_independent_model_view(
    const MVMCPowerLanczosIndependentModel *model);

/* Validate and retain a non-owning view of the physical Trans/InterAll rows. */
MVMCKrylovStatus mvmc_power_lanczos_independent_streaming_model_create(
    int site_count, int up_electron_count, int down_electron_count,
    size_t transfer_count, int *const *transfer_indices,
    const double complex *transfer_parameters, size_t inter_all_count,
    int *const *inter_all_indices, const double complex *inter_all_parameters,
    MVMCPowerLanczosIndependentStreamingModel **result);

void mvmc_power_lanczos_independent_streaming_model_destroy(
    MVMCPowerLanczosIndependentStreamingModel *model);

/* Exact depth-one H action with O(configuration + term) storage. */
MVMCKrylovStatus mvmc_power_lanczos_independent_stream_apply(
    const MVMCKrylovFockModel *model, const uint64_t *configuration_words,
    size_t word_count, MVMCKrylovScaledAmplitudeCallback child_amplitude,
    void *child_context, MVMCScaledComplex *result, uint64_t *applied_terms,
    uint64_t *amplitude_calls);

/* Same action without copying the physical Hamiltonian into term objects. */
MVMCKrylovStatus mvmc_power_lanczos_independent_stream_apply_view(
    const MVMCPowerLanczosIndependentStreamingModel *model,
    const uint64_t *configuration_words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback child_amplitude, void *child_context,
    MVMCScaledComplex *result, uint64_t *applied_terms,
    uint64_t *amplitude_calls);

/* q=(Psi,H'Psi), r=(HPsi,HH'Psi); H' workspace has max_order=1. */
MVMCKrylovStatus mvmc_power_lanczos_independent_evaluate(
    const MVMCKrylovFockModel *physical_model,
    MVMCKrylovBoundedWorkspace *hprime_workspace,
    const uint64_t *configuration_words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCPowerLanczosIndependentSnapshot *result);

MVMCKrylovStatus mvmc_power_lanczos_independent_evaluate_view(
    const MVMCPowerLanczosIndependentStreamingModel *physical_model,
    MVMCKrylovBoundedWorkspace *hprime_workspace,
    const uint64_t *configuration_words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCPowerLanczosIndependentSnapshot *result);

/* One common-normalization contribution for the 2x2 GEVP. */
MVMCKrylovStatus mvmc_power_lanczos_independent_matrix_sample(
    const MVMCPowerLanczosIndependentSnapshot *snapshot, double eta,
    double complex overlap_upper[3],
    double complex hamiltonian_forward_upper[3],
    double complex hamiltonian_reverse_upper[3],
    double complex hamiltonian_squared_upper[3], double *guide_denominator);

/* Phi=c.q and HPhi=c.r in scaled/log form. */
MVMCKrylovStatus mvmc_power_lanczos_independent_final_evaluate(
    const MVMCPowerLanczosIndependentSnapshot *snapshot,
    const double complex coefficient[2], MVMCScaledComplex *phi,
    MVMCScaledComplex *h_phi, MVMCScaledComplex *local_energy,
    double *log_weight);

MVMCKrylovStatus mvmc_power_lanczos_independent_run(
    const MVMCPowerLanczosIndependentInput *input,
    MVMCPowerLanczosIndependentResult *result);

#endif

#endif
