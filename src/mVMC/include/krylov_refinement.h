#ifndef MVMC_KRYLOV_REFINEMENT_H
#define MVMC_KRYLOV_REFINEMENT_H

#include "krylov_final_state.h"

/* The opaque base value is q * exp(g), with q complex rational and g real
 * rational.  The builder imports finite binary64 inputs exactly. */
typedef struct MVMCRefinedBase MVMCRefinedBase;

typedef MVMCKrylovStatus (*MVMCRefinedBaseCallback)(
    const uint64_t *words, size_t word_count, void *context,
    MVMCRefinedBase *value);

typedef struct {
  size_t max_matrix_entries;
  size_t max_cache_entries;
  size_t max_expression_terms;
  size_t max_rational_bits;
  uint64_t max_node_calls;
  unsigned max_precision;
} MVMCKrylovRefinementLimits;

typedef struct {
  int valid;
  int accepted;
  int exact_zero_proposal;
  unsigned precision;
  uint64_t node_calls;
  MVMCKrylovBoundedResult krylov;
  MVMCKrylovFinalStateEvaluation evaluation;
} MVMCKrylovRefinementResult;

/* Local, allocation-owning refinement.  It makes no MPI calls or RNG draws.
 * Only symmetric proposals are supported.  Failure never modifies *result.
 * Rejection does not require a representable local energy.  Acceptance
 * supplies a refined evaluation and Krylov values for snapshot replacement.
 * Resource exhaustion and unresolved intervals remain explicit failures. */
MVMCKrylovStatus mvmc_krylov_refine_final(
    const MVMCKrylovFockModel *model,
    const MVMCKrylovFinalStatePolicy *policy,
    const uint64_t *current_words, const uint64_t *proposal_words,
    size_t word_count, double uniform,
    MVMCRefinedBaseCallback amplitude, void *context,
    const MVMCKrylovRefinementLimits *limits,
    MVMCKrylovRefinementResult *result);

/* Pf(Slater[x,x]) equals the production (-1)^(N/2) Pf(-Slater[x,x]).
 * Slater is QP-major, each full matrix row-major.  Parameters enter through
 * exp(sum Re(parameter[i])*counts[i]).  A proven support zero is applied
 * only after all supplied components and parameters have been checked. */
MVMCKrylovStatus mvmc_refined_base_build(
    MVMCRefinedBase *result, const double complex *slater,
    size_t orbital_count, const int *occupied, size_t electron_count,
    const double complex *weights, size_t qp_count, int support_zero,
    const double complex *parameters, const int64_t *counts, size_t nproj);

#endif
