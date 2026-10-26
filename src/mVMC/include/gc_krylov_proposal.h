/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_GC_KRYLOV_PROPOSAL_H
#define MVMC_GC_KRYLOV_PROPOSAL_H

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                         \
    defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "gc_krylov_spin_support.h"
#include "krylov_fock_proposal.h"
#include "krylov_fock_reference.h"

#include <stddef.h>
#include <stdint.h>

/*
 * Symmetric configuration proposals for the grand-canonical domains, used
 * by the chain root only (no MPI, no Hamiltonian, no allocation).  Every
 * random integer comes from draw(context, bound) in [0, bound); a callback
 * error is returned and an out-of-range value is an invariant failure.  On
 * any failure *out is not a usable candidate.  out must not alias current.
 *
 * uniform: every configuration of the domain with equal probability, by a
 *   bijective construction (GC_EVEN: 2L-1 free bits plus the parity bit;
 *   GC_BALANCED: an L-element subset of the labels U_i, Dbar_i selected
 *   sequentially; U_i selected means up_i occupied, Dbar_i not selected
 *   means down_i occupied).  No combinatorial counts or floating point.
 * initial: uniform among the configurations with initial_n electrons
 *   (GC_BALANCED: initial_n/2 up and initial_n/2 down sites).
 * local: GC_EVEN flips two distinct orbitals chosen as an ordered pair
 *   (pair addition, removal or hop); GC_BALANCED with probability 1/2 flips
 *   (up_i, down_j) when both are empty or both occupied, otherwise with
 *   probability 1/2 exchanges two distinct same-spin sites of one spin when
 *   their occupations differ; every other draw (and L = 1 hops) is a self
 *   move.  Forward and reverse probabilities are equal, so
 *   log(reverse/forward) = 0; nothing is renormalized by availability.
 */
MVMCKrylovStatus mvmc_gc_krylov_draw_uniform(
    const MVMCKrylovFockModel *model,
    MVMCKrylovFockProposalBoundedDraw draw, void *context,
    uint64_t *out, size_t word_count);
MVMCKrylovStatus mvmc_gc_krylov_draw_initial(
    const MVMCKrylovFockModel *model, int initial_n,
    MVMCKrylovFockProposalBoundedDraw draw, void *context,
    uint64_t *out, size_t word_count);
MVMCKrylovStatus mvmc_gc_krylov_draw_local(
    const MVMCKrylovFockModel *model, const uint64_t *current,
    MVMCKrylovFockProposalBoundedDraw draw, void *context,
    uint64_t *out, size_t word_count);

/* Balanced-domain N +/- 4 proposal: choose two distinct sites of each spin
 * uniformly from all sites, then flip all four only if their occupations
 * agree. Otherwise self-loop. No availability conditioning: Q is symmetric. */
MVMCKrylovStatus mvmc_gc_krylov_draw_two_pairs(
    const MVMCKrylovFockModel *model, const uint64_t *current,
    MVMCKrylovFockProposalBoundedDraw draw, void *context,
    uint64_t *out, size_t word_count);

/*
 * Finite-S (active, S > 0, GC_BALANCED) start sectors.  N = 2m is allowed
 * when max(0, 2m - L) <= m - S and, for a TRIPLET pair function, m - S is
 * even.  explicit_n = 1 keeps requested_n alone (a forbidden N is
 * INVALID_ARGUMENT); explicit_n = 0 lists every allowed N by
 * (|N - requested_n|, N).  An S > 0 SINGLET/ZERO support is an
 * AMPLITUDE_FAILURE.  requested_n must be even in [0, 2L] (never rounded);
 * capacity must be at least L + 1.  No allocation, RNG or MPI; on failure
 * sectors and count are unchanged.  The integer test does not prove that a
 * general pair function is nonzero there.
 */
MVMCKrylovStatus mvmc_gc_krylov_initial_sectors(
    const MVMCKrylovFockModel *model, const MVMCGCSpinSupport *support,
    int requested_n, int explicit_n,
    int *sectors, size_t capacity, size_t *count);

/*
 * Structured finite-S start configuration with initial_n (an allowed
 * sector) electrons, for the chain root only (no allocation, no MPI).  D is
 * uniform in [max(0, 2m - L), m - S]; then D doubly occupied sites, m - D
 * up-only and m - D down-only sites are chosen as uniform subsets of the
 * still unused sites in ascending order.  A choice with a single option
 * (one D value, zero or all remaining sites) draws nothing, so every draw
 * has bound >= 2.  This is not uniform over all supported configurations
 * and is used for initialization only.  out must not alias a current
 * configuration; on failure it is not a usable candidate.
 */
MVMCKrylovStatus mvmc_gc_krylov_draw_spin_initial(
    const MVMCKrylovFockModel *model, const MVMCGCSpinSupport *support,
    int initial_n, MVMCKrylovFockProposalBoundedDraw draw, void *context,
    uint64_t *out, size_t word_count);

#endif /* reference or bounded engine */

#endif /* MVMC_GC_KRYLOV_PROPOSAL_H */
