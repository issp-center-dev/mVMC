#ifndef MVMC_POWER_LANCZOS_GEVP_H
#define MVMC_POWER_LANCZOS_GEVP_H

#include "krylov_gevp_solver.h"

#include <complex.h>
#include <stddef.h>

#define MVMC_POWER_LANCZOS_GEVP_POLICY_VERSION 2
#define MVMC_POWER_LANCZOS_GEVP_POLICY_ID \
  "power_lanczos_p6c1_exact_cartesian_normwise_gevp_v2"

typedef enum {
  MVMC_POWER_LANCZOS_GEVP_CUTOFF_S48 = 0,
  MVMC_POWER_LANCZOS_GEVP_CUTOFF_S40,
  MVMC_POWER_LANCZOS_GEVP_CUTOFF_S32,
  MVMC_POWER_LANCZOS_GEVP_CUTOFF_S24
} MVMCPowerLanczosGEVPCutoffID;

typedef struct {
  int valid;
  int policy_version;
  MVMCPowerLanczosGEVPCutoffID cutoff_id;
  double rank_relative_cutoff;
  double degenerate_root_gap_relative_threshold;
  double maximum_normwise_backward_error;
  double negative_variance_relative_tolerance;
} MVMCPowerLanczosGEVPPolicy;

typedef struct {
  int valid;
  int observables_valid;
  MVMCKrylovGEVPStatus status;
  MVMCKrylovGEVPStatus eigenspace_status;
  int policy_version;
  int dimension;
  int retained_rank;
  int discarded_rank;
  int root_multiplicity;
  int phase_pivot;
  int variance_clamped;
  double rank_relative_cutoff;
  double overlap_eigenvalue[MVMC_KRYLOV_GEVP_MAX_DIMENSION];
  double condition_estimate;
  double root_gap;
  double normalization;
  double energy;
  double energy_squared;
  double variance;
  double raw_action_relative_residual;
  double normwise_backward_error;
  double relative_residual;
  double eigenspace_relative_residual;
  double eigenspace_normalization;
  double complex coefficient[MVMC_KRYLOV_GEVP_MAX_DIMENSION];
  double complex root_subspace_projector[
      MVMC_KRYLOV_GEVP_MAX_DIMENSION * MVMC_KRYLOV_GEVP_MAX_DIMENSION];
} MVMCPowerLanczosGEVPResult;

MVMCKrylovGEVPStatus mvmc_power_lanczos_gevp_default_policy(
    double rank_relative_cutoff, MVMCPowerLanczosGEVPPolicy *policy);

/*
 * Inputs use the row-then-column upper-packed convention of the P4
 * packed-matrix engine.  The P6-C1 result is independently canonicalized by
 * projecting e0,e1,... into the lowest-root subspace with the S inner
 * product, S-normalizing, phase fixing, and recomputing the P6 residual.
 * All input arrays, policy, and result must be mutually disjoint.  An alias
 * failure performs no write.  Version 2 is release-qualified only for the
 * registered P6-C1 exact Cartesian cells, whose accepted positive fixtures
 * have full retained rank and a unique lowest root at every cutoff.  The
 * acceptance residual is the normwise GEVP backward error; the raw action
 * residual is retained as a conditioning diagnostic.  Matrix-difference
 * variance is diagnostic only because C1 validates release variance from
 * the corrected physical state.  observables_valid is false when that raw
 * matrix-difference diagnostic is materially negative; valid still denotes
 * a usable, independently certified coefficient.
 */
MVMCKrylovGEVPStatus mvmc_power_lanczos_gevp_solve_complex_packed(
    const MVMCPowerLanczosGEVPPolicy *policy, int dimension,
    const double complex *overlap_upper,
    const double complex *hamiltonian_forward_upper,
    const double complex *hamiltonian_reverse_upper,
    const double complex *hamiltonian_squared_upper, size_t upper_count,
    MVMCPowerLanczosGEVPResult *result);

MVMCKrylovGEVPStatus mvmc_power_lanczos_gevp_solve_real_packed(
    const MVMCPowerLanczosGEVPPolicy *policy, int dimension,
    const double *overlap_upper, const double *hamiltonian_forward_upper,
    const double *hamiltonian_reverse_upper,
    const double *hamiltonian_squared_upper, size_t upper_count,
    MVMCPowerLanczosGEVPResult *result);

const char *mvmc_power_lanczos_gevp_policy_id(void);

#endif
