#ifndef MVMC_POWER_LANCZOS_OBSERVABLE_EVALUATOR_H
#define MVMC_POWER_LANCZOS_OBSERVABLE_EVALUATOR_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "bounded_krylov_engine.h"
#include "power_lanczos_observable_action.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_OBSERVABLE_MATRIX_ENTRIES 4

typedef struct MVMCPowerLanczosObservableEvaluatorWorkspace
    MVMCPowerLanczosObservableEvaluatorWorkspace;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  size_t request_count;
  size_t active_request_count;
  size_t exact_zero_request_count;
  size_t unique_target_count;
  uint64_t engine_root_evaluations;
  int source_target_reused;
  double minimum_log_contribution;
  double maximum_log_contribution;
} MVMCPowerLanczosObservableSampleDiagnostics;

MVMCKrylovStatus mvmc_power_lanczos_observable_evaluator_workspace_create(
    size_t max_requests, size_t word_count,
    MVMCPowerLanczosObservableEvaluatorWorkspace **workspace);

void mvmc_power_lanczos_observable_evaluator_workspace_destroy(
    MVMCPowerLanczosObservableEvaluatorWorkspace *workspace);

size_t mvmc_power_lanczos_observable_evaluator_workspace_bytes(
    const MVMCPowerLanczosObservableEvaluatorWorkspace *workspace);

/*
 * Evaluate one coefficient-chain sample.  The matrix output is ordered as
 * request-major (0,0),(0,1),(1,0),(1,1).  basis_log_scale defines the
 * represented basis w_i=exp(basis_log_scale[i])*v_i.  log_guide is log q(y).
 * The bounded workspace must already own an active persistent session.
 */
MVMCKrylovStatus mvmc_power_lanczos_observable_coefficient_sample(
    MVMCPowerLanczosObservableEvaluatorWorkspace *workspace,
    MVMCKrylovBoundedWorkspace *bounded_workspace,
    const MVMCPowerLanczosObservableLayout *layout,
    const MVMCPowerLanczosObservablePlan *plan,
    const uint64_t *source_words, size_t word_count,
    double log_guide, const double basis_log_scale[2],
    double complex *matrix_entries, size_t matrix_entry_capacity,
    MVMCPowerLanczosObservableSampleDiagnostics *diagnostics);

/*
 * Evaluate the independent final-chain primitive using
 * A=alpha[0]*w_0+alpha[1]*w_1 directly.  This path never contracts the
 * coefficient-stage matrix.
 */
MVMCKrylovStatus mvmc_power_lanczos_observable_final_sample(
    MVMCPowerLanczosObservableEvaluatorWorkspace *workspace,
    MVMCKrylovBoundedWorkspace *bounded_workspace,
    const MVMCPowerLanczosObservableLayout *layout,
    const MVMCPowerLanczosObservablePlan *plan,
    const uint64_t *source_words, size_t word_count,
    double log_guide, const double basis_log_scale[2],
    const double complex alpha[2],
    double complex *numerators, size_t numerator_capacity,
    MVMCPowerLanczosObservableSampleDiagnostics *diagnostics);

#endif

#endif
