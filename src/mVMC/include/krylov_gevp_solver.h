#ifndef MVMC_KRYLOV_GEVP_SOLVER_H
#define MVMC_KRYLOV_GEVP_SOLVER_H

#include <complex.h>
#include <stddef.h>

#define MVMC_KRYLOV_GEVP_POLICY_VERSION 2
#define MVMC_KRYLOV_GEVP_POLICY_ID \
  "power_lanczos_zero_support_p4f_solver_v2"

enum { MVMC_KRYLOV_GEVP_MAX_DIMENSION = 3 };

typedef enum {
  MVMC_KRYLOV_GEVP_OK = 0,
  MVMC_KRYLOV_GEVP_INVALID_ARGUMENT,
  MVMC_KRYLOV_GEVP_NONFINITE_INPUT,
  MVMC_KRYLOV_GEVP_ANTIHERMITIAN_INPUT,
  MVMC_KRYLOV_GEVP_NONPOSITIVE_OVERLAP_DIAGONAL,
  MVMC_KRYLOV_GEVP_INDEFINITE_OVERLAP,
  MVMC_KRYLOV_GEVP_ZERO_RETAINED_RANK,
  MVMC_KRYLOV_GEVP_INDEFINITE_HAMILTONIAN_SQUARED,
  MVMC_KRYLOV_GEVP_LAPACK_FAILURE,
  MVMC_KRYLOV_GEVP_NORMALIZATION_FAILURE,
  MVMC_KRYLOV_GEVP_RESIDUAL_FAILURE,
  MVMC_KRYLOV_GEVP_NEGATIVE_VARIANCE
} MVMCKrylovGEVPStatus;

typedef struct {
  int valid;
  int policy_version;
  double rank_relative_cutoff;
  double overlap_negative_relative_tolerance;
  double hamiltonian_squared_negative_relative_tolerance;
  double maximum_input_antihermitian_effect;
  double maximum_gevp_relative_residual;
  double negative_variance_relative_tolerance;
  double degenerate_root_gap_relative_threshold;
} MVMCKrylovGEVPPolicy;

typedef struct {
  int valid;
  MVMCKrylovGEVPStatus status;
  int policy_version;
  int dimension;
  int retained_rank;
  int discarded_rank;
  int root_multiplicity;
  int coefficient_comparison_uses_projector;
  int phase_pivot;
  int variance_clamped;
  int lapack_info_overlap_query;
  int lapack_info_overlap;
  int lapack_info_hamiltonian_squared_query;
  int lapack_info_hamiltonian_squared;
  int lapack_info_projected_query;
  int lapack_info_projected;
  double rank_relative_cutoff;
  double overlap_eigenvalue[MVMC_KRYLOV_GEVP_MAX_DIMENSION];
  double hamiltonian_squared_eigenvalue[MVMC_KRYLOV_GEVP_MAX_DIMENSION];
  double condition_estimate;
  double overlap_diagonal_imaginary_residual;
  double hamiltonian_antihermitian_residual;
  double hamiltonian_squared_diagonal_imaginary_residual;
  double overlap_norm;
  double hamiltonian_norm;
  double hamiltonian_squared_norm;
  double projected_energy;
  double root_gap;
  double normalization;
  double energy;
  double energy_squared;
  double variance;
  double gevp_relative_residual;
  double complex coefficient[MVMC_KRYLOV_GEVP_MAX_DIMENSION];
  double complex root_subspace_projector[
      MVMC_KRYLOV_GEVP_MAX_DIMENSION * MVMC_KRYLOV_GEVP_MAX_DIMENSION];
} MVMCKrylovGEVPResult;

MVMCKrylovGEVPStatus mvmc_krylov_gevp_default_policy(
    double rank_relative_cutoff, MVMCKrylovGEVPPolicy *policy);

/* Arrays use the same row-then-column upper-packed order as
 * mvmc_krylov_streaming_upper_index().  For i <= j, forward stores K_ij and
 * reverse stores the independently measured raw lower entry K_ji.  The
 * solver forms (forward_ij + conj(reverse_ij))/2.  All four arrays must have
 * one common positive normalization and must not overlap each other, policy,
 * or result.  An alias failure performs no write; otherwise valid==0 makes
 * coefficients and observables unusable even when diagnostics are present. */
MVMCKrylovGEVPStatus mvmc_krylov_gevp_solve_complex_packed(
    const MVMCKrylovGEVPPolicy *policy, int dimension,
    const double complex *overlap_upper,
    const double complex *hamiltonian_forward_upper,
    const double complex *hamiltonian_reverse_upper,
    const double complex *hamiltonian_squared_upper, size_t upper_count,
    MVMCKrylovGEVPResult *result);

MVMCKrylovGEVPStatus mvmc_krylov_gevp_solve_real_packed(
    const MVMCKrylovGEVPPolicy *policy, int dimension,
    const double *overlap_upper, const double *hamiltonian_forward_upper,
    const double *hamiltonian_reverse_upper,
    const double *hamiltonian_squared_upper, size_t upper_count,
    MVMCKrylovGEVPResult *result);

const char *mvmc_krylov_gevp_status_string(MVMCKrylovGEVPStatus status);

#endif
