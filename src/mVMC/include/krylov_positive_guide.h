/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_KRYLOV_POSITIVE_GUIDE_H
#define MVMC_KRYLOV_POSITIVE_GUIDE_H

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) &&                         \
    defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "bounded_krylov_engine.h"

typedef struct {
  int order;
  double eta;
  double lambda[MVMC_KRYLOV_MAX_ORDER + 1];
  double log_basis_scale[MVMC_KRYLOV_MAX_ORDER + 1];
  uint64_t policy_hash;
} MVMCKrylovPositiveGuidePolicy;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  double log_guide;
  double log_floor;
  double log_signal_sum;
  int finite_component_count;
  int exact_zero_component_count;
  int numeric_zero_component_count;
} MVMCKrylovPositiveGuideResult;

typedef struct {
  int valid;
  MVMCKrylovStatus status;
  double log_acceptance_ratio;
  double uniform;
  int accepted;
} MVMCKrylovPositiveGuideAcceptance;

MVMCKrylovStatus mvmc_krylov_positive_guide_evaluate(
    const MVMCKrylovPositiveGuidePolicy *policy,
    const MVMCScaledComplex *values, size_t value_count,
    MVMCKrylovPositiveGuideResult *result);

MVMCKrylovStatus mvmc_krylov_positive_guide_acceptance(
    const MVMCKrylovPositiveGuideResult *current,
    const MVMCKrylovPositiveGuideResult *proposal,
    double log_proposal_ratio, double uniform,
    MVMCKrylovPositiveGuideAcceptance *result);

#endif

#endif /* MVMC_KRYLOV_POSITIVE_GUIDE_H */
