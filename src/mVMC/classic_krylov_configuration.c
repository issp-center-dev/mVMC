/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "classic_krylov_configuration.h"

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                     \
    defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include <limits.h>
#include <string.h>

static int valid_layout(const MVMCClassicKrylovConfigurationLayout *layout,
                        size_t *orbital_count, size_t *word_count) {
  size_t orbitals;

  if (layout == NULL || orbital_count == NULL || word_count == NULL ||
      layout->site_count <= 0 || layout->site_count > INT_MAX / 2 ||
      layout->up_electron_count <= 0 ||
      layout->up_electron_count != layout->down_electron_count ||
      layout->up_electron_count > layout->site_count ||
      (layout->pure_spin != 0 && layout->pure_spin != 1) ||
      (layout->pure_spin &&
       layout->up_electron_count + layout->down_electron_count !=
           layout->site_count)) {
    return 0;
  }
  orbitals = (size_t)layout->site_count * 2u;
  if (orbitals > SIZE_MAX - 63u) return 0;
  *orbital_count = orbitals;
  *word_count = (orbitals + 63u) / 64u;
  return *word_count != 0;
}

size_t mvmc_classic_krylov_configuration_word_count(
    const MVMCClassicKrylovConfigurationLayout *layout) {
  size_t orbital_count = 0;
  size_t word_count = 0;
  return valid_layout(layout, &orbital_count, &word_count) ? word_count : 0;
}

static int valid_word_padding(const uint64_t *words, size_t word_count,
                              size_t orbital_count) {
  const size_t used = orbital_count % 64u;
  return used == 0u || (words[word_count - 1u] >> used) == 0u;
}

static int occupied(const uint64_t *words, size_t orbital) {
  return (int)((words[orbital / 64u] >> (orbital % 64u)) & UINT64_C(1));
}

MVMCKrylovStatus mvmc_classic_krylov_configuration_decode(
    const MVMCClassicKrylovConfigurationLayout *layout,
    const uint64_t *configuration_words, size_t word_count,
    int *ele_idx, int *ele_cfg, int *ele_num) {
  size_t orbital_count = 0;
  size_t expected_words = 0;
  int up_count = 0;
  int down_count = 0;
  size_t orbital;

  if (!valid_layout(layout, &orbital_count, &expected_words) ||
      configuration_words == NULL || word_count != expected_words ||
      ele_idx == NULL || ele_cfg == NULL || ele_num == NULL ||
      !valid_word_padding(configuration_words, word_count, orbital_count)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (orbital = 0; orbital < orbital_count; ++orbital) {
    if (!occupied(configuration_words, orbital)) continue;
    if (orbital < (size_t)layout->site_count) {
      if (up_count == layout->up_electron_count) {
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      }
      ++up_count;
    } else {
      if (down_count == layout->down_electron_count) {
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      }
      ++down_count;
    }
  }
  if (up_count != layout->up_electron_count ||
      down_count != layout->down_electron_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (layout->pure_spin) {
    int site;
    for (site = 0; site < layout->site_count; ++site) {
      if (occupied(configuration_words, (size_t)site) +
              occupied(configuration_words,
                       (size_t)layout->site_count + (size_t)site) !=
          1) {
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
      }
    }
  }

  memset(ele_idx, 0,
         (size_t)(layout->up_electron_count +
                  layout->down_electron_count) * sizeof(*ele_idx));
  for (orbital = 0; orbital < orbital_count; ++orbital) {
    ele_num[orbital] = occupied(configuration_words, orbital);
    ele_cfg[orbital] = -1;
  }
  up_count = 0;
  down_count = 0;
  for (orbital = 0; orbital < orbital_count; ++orbital) {
    const int site = (int)(orbital % (size_t)layout->site_count);
    if (!ele_num[orbital]) continue;
    if (orbital < (size_t)layout->site_count) {
      ele_idx[up_count] = site;
      ele_cfg[orbital] = up_count;
      ++up_count;
    } else {
      ele_idx[layout->up_electron_count + down_count] = site;
      ele_cfg[orbital] = down_count;
      ++down_count;
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static int valid_classic_views(
    const MVMCClassicKrylovConfigurationLayout *layout,
    const int *ele_idx, const int *ele_cfg, const int *ele_num) {
  int spin;
  int site;
  int label;

  for (spin = 0; spin < 2; ++spin) {
    const int count = spin == 0 ? layout->up_electron_count
                                : layout->down_electron_count;
    const int offset = spin == 0 ? 0 : layout->up_electron_count;
    int occupied_count = 0;
    for (site = 0; site < layout->site_count; ++site) {
      const int orbital = site + spin * layout->site_count;
      if ((ele_num[orbital] != 0 && ele_num[orbital] != 1) ||
          (ele_num[orbital] == 0 && ele_cfg[orbital] != -1) ||
          (ele_num[orbital] == 1 &&
           (ele_cfg[orbital] < 0 || ele_cfg[orbital] >= count))) {
        return 0;
      }
      if (ele_num[orbital] == 1) ++occupied_count;
    }
    if (occupied_count != count) return 0;
    for (label = 0; label < count; ++label) {
      const int particle_site = ele_idx[offset + label];
      const int orbital = particle_site + spin * layout->site_count;
      if (particle_site < 0 || particle_site >= layout->site_count ||
          ele_num[orbital] != 1 || ele_cfg[orbital] != label) {
        return 0;
      }
    }
  }
  if (layout->pure_spin) {
    for (site = 0; site < layout->site_count; ++site) {
      if (ele_num[site] + ele_num[layout->site_count + site] != 1) return 0;
    }
  }
  return 1;
}

MVMCKrylovStatus mvmc_classic_krylov_configuration_encode(
    const MVMCClassicKrylovConfigurationLayout *layout,
    const int *ele_idx, const int *ele_cfg, const int *ele_num,
    uint64_t *configuration_words, size_t word_count) {
  size_t orbital_count = 0;
  size_t expected_words = 0;
  size_t orbital;

  if (!valid_layout(layout, &orbital_count, &expected_words) ||
      ele_idx == NULL || ele_cfg == NULL || ele_num == NULL ||
      configuration_words == NULL || word_count != expected_words ||
      !valid_classic_views(layout, ele_idx, ele_cfg, ele_num)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  memset(configuration_words, 0,
         word_count * sizeof(*configuration_words));
  for (orbital = 0; orbital < orbital_count; ++orbital) {
    if (ele_num[orbital]) {
      configuration_words[orbital / 64u] |=
          UINT64_C(1) << (orbital % 64u);
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

#endif /* reference or bounded engine */
