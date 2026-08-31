/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_CLASSIC_KRYLOV_CONFIGURATION_H
#define MVMC_CLASSIC_KRYLOV_CONFIGURATION_H

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                     \
    defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "krylov_fock_reference.h"

#include <stddef.h>
#include <stdint.h>

typedef struct {
  int site_count;
  int up_electron_count;
  int down_electron_count;
  int pure_spin;
} MVMCClassicKrylovConfigurationLayout;

/*
 * Canonical spin-orbital order is every up-spin site followed by every
 * down-spin site.  The classic non-FSZ P6 contract is deliberately narrower:
 * both spin populations must be positive and equal (2Sz=0).  Pure-spin
 * configurations additionally contain exactly one electron per site.
 */
size_t mvmc_classic_krylov_configuration_word_count(
    const MVMCClassicKrylovConfigurationLayout *layout);

/*
 * Decode canonical words into the classic layout.  Particle labels in eleCfg
 * are assigned in ascending occupied-site order independently for each spin;
 * eleIdx stores up labels first and down labels second.  Error returns do not
 * modify any output buffer.
 */
MVMCKrylovStatus mvmc_classic_krylov_configuration_decode(
    const MVMCClassicKrylovConfigurationLayout *layout,
    const uint64_t *configuration_words, size_t word_count,
    int *ele_idx, int *ele_cfg, int *ele_num);

/*
 * Validate all three classic views against one another and encode canonical
 * words transactionally.  This is intentionally stricter than consulting
 * eleNum alone so stale particle labels cannot cross the production adapter.
 */
MVMCKrylovStatus mvmc_classic_krylov_configuration_encode(
    const MVMCClassicKrylovConfigurationLayout *layout,
    const int *ele_idx, const int *ele_cfg, const int *ele_num,
    uint64_t *configuration_words, size_t word_count);

#endif /* reference or bounded engine */

#endif /* MVMC_CLASSIC_KRYLOV_CONFIGURATION_H */
