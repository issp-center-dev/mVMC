/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "classic_krylov_projection.h"

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE) ||                     \
    defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct MVMCClassicKrylovProjectionWorkspace {
  size_t site_count;
  size_t orbital_count;
  size_t nproj;
  int ngutzwiller_idx;
  int njastrow_idx;
  int nspin_jastrow_idx;
  int ndoublon_holon_2site_idx;
  int ndoublon_holon_4site_idx;
  int *gutzwiller_idx;
  int *jastrow_idx;
  int *spin_jastrow_idx;
  int *doublon_holon_2site_idx;
  int *doublon_holon_4site_idx;
  double complex *parameters;
  int64_t *scratch_counts;
  unsigned char *binding_bytes;
  size_t binding_byte_count;
};

static int checked_multiply(size_t left, size_t right, size_t *product) {
  if (product == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *product = left * right;
  return 1;
}

static int checked_add(size_t left, size_t right, size_t *sum) {
  if (sum == NULL || right > SIZE_MAX - left) return 0;
  *sum = left + right;
  return 1;
}

static int append_bytes(unsigned char **cursor, size_t *remaining,
                        const void *source, size_t byte_count) {
  if (cursor == NULL || *cursor == NULL || remaining == NULL ||
      byte_count > *remaining || (byte_count != 0 && source == NULL)) {
    return 0;
  }
  if (byte_count != 0) memcpy(*cursor, source, byte_count);
  *cursor += byte_count;
  *remaining -= byte_count;
  return 1;
}

static int expected_projection_count(
    const MVMCClassicKrylovProjectionLayout *layout, size_t *expected) {
  size_t value;

  if (layout == NULL || expected == NULL || layout->site_count == 0 ||
      layout->site_count > (size_t)INT_MAX ||
      layout->site_count > SIZE_MAX / 4u || layout->nproj < 0 ||
      layout->ngutzwiller_idx < 0 || layout->njastrow_idx < 0 ||
      layout->nspin_jastrow_idx < 0 ||
      layout->ndoublon_holon_2site_idx < 0 ||
      layout->ndoublon_holon_4site_idx < 0) {
    return 0;
  }
  value = (size_t)layout->ngutzwiller_idx;
  if (!checked_add(value, (size_t)layout->njastrow_idx, &value) ||
      !checked_add(value, (size_t)layout->nspin_jastrow_idx, &value) ||
      (size_t)layout->ndoublon_holon_2site_idx >
          (SIZE_MAX - value) / 6u) {
    return 0;
  }
  value += 6u * (size_t)layout->ndoublon_holon_2site_idx;
  if ((size_t)layout->ndoublon_holon_4site_idx >
      (SIZE_MAX - value) / 10u) {
    return 0;
  }
  value += 10u * (size_t)layout->ndoublon_holon_4site_idx;
  if (value > (size_t)INT_MAX || value != (size_t)layout->nproj) return 0;
  *expected = value;
  return 1;
}

static int valid_jastrow_binding(const int *const *binding,
                                 size_t site_count, int class_count) {
  size_t row;
  size_t column;

  if (class_count == 0) return 1;
  if (binding == NULL) return 0;
  for (row = 0; row < site_count; ++row) {
    if (binding[row] == NULL) return 0;
    for (column = 0; column < site_count; ++column) {
      const int value = binding[row][column];
      if (row == column) {
        if (value != -1) return 0;
      } else if (value < 0 || value >= class_count ||
                 binding[column] == NULL ||
                 binding[column][row] != value) {
        return 0;
      }
    }
  }
  return 1;
}

static int valid_dh_binding(const int *const *binding, size_t site_count,
                            int shell_count, size_t neighbor_count) {
  size_t row_count;
  int shell;
  size_t index;

  if (shell_count == 0) return 1;
  if (binding == NULL ||
      !checked_multiply(site_count, neighbor_count, &row_count)) {
    return 0;
  }
  for (shell = 0; shell < shell_count; ++shell) {
    if (binding[shell] == NULL) return 0;
    for (index = 0; index < row_count; ++index) {
      if (binding[shell][index] < 0 ||
          (size_t)binding[shell][index] >= site_count) {
        return 0;
      }
    }
  }
  return 1;
}

static int copy_flat_rows(int **destination, const int *const *source,
                          size_t row_count, size_t column_count) {
  size_t element_count;
  size_t byte_count;
  size_t row;

  if (destination == NULL) return 0;
  *destination = NULL;
  if (row_count == 0 || column_count == 0) return 1;
  if (source == NULL || !checked_multiply(row_count, column_count,
                                          &element_count) ||
      !checked_multiply(element_count, sizeof(**destination), &byte_count)) {
    return 0;
  }
  *destination = (int *)malloc(byte_count);
  if (*destination == NULL) return 0;
  for (row = 0; row < row_count; ++row) {
    memcpy(*destination + row * column_count, source[row],
           column_count * sizeof(**destination));
  }
  return 1;
}

static int allocate_and_copy(
    const MVMCClassicKrylovProjectionLayout *layout,
    MVMCClassicKrylovProjectionWorkspace *workspace) {
  size_t byte_count;

  if (layout->ngutzwiller_idx != 0) {
    if (!checked_multiply(layout->site_count,
                          sizeof(*workspace->gutzwiller_idx), &byte_count)) {
      return 0;
    }
    workspace->gutzwiller_idx = (int *)malloc(byte_count);
    if (workspace->gutzwiller_idx == NULL) return 0;
    memcpy(workspace->gutzwiller_idx, layout->gutzwiller_idx, byte_count);
  }
  if (!copy_flat_rows(&workspace->jastrow_idx, layout->jastrow_idx,
                      layout->njastrow_idx == 0 ? 0 : layout->site_count,
                      layout->site_count) ||
      !copy_flat_rows(&workspace->spin_jastrow_idx,
                      layout->spin_jastrow_idx,
                      layout->nspin_jastrow_idx == 0 ? 0 : layout->site_count,
                      layout->site_count) ||
      !copy_flat_rows(&workspace->doublon_holon_2site_idx,
                      layout->doublon_holon_2site_idx,
                      (size_t)layout->ndoublon_holon_2site_idx,
                      2u * layout->site_count) ||
      !copy_flat_rows(&workspace->doublon_holon_4site_idx,
                      layout->doublon_holon_4site_idx,
                      (size_t)layout->ndoublon_holon_4site_idx,
                      4u * layout->site_count)) {
    return 0;
  }
  if (workspace->nproj != 0) {
    if (!checked_multiply(workspace->nproj,
                          sizeof(*workspace->parameters), &byte_count)) {
      return 0;
    }
    workspace->parameters = (double complex *)malloc(byte_count);
    workspace->scratch_counts = (int64_t *)calloc(
        workspace->nproj, sizeof(*workspace->scratch_counts));
    if (workspace->parameters == NULL || workspace->scratch_counts == NULL) {
      return 0;
    }
    memcpy(workspace->parameters, layout->parameters, byte_count);
  }
  return 1;
}

static int add_array_bytes(size_t count, size_t element_size,
                           size_t *total) {
  size_t bytes;
  return checked_multiply(count, element_size, &bytes) &&
         checked_add(*total, bytes, total);
}

static int create_binding_record(
    MVMCClassicKrylovProjectionWorkspace *workspace) {
  uint64_t metadata[7];
  size_t jastrow_count;
  size_t dh2_count;
  size_t dh4_count;
  size_t byte_count = sizeof(metadata);
  unsigned char *cursor;
  size_t remaining;

  if (!checked_multiply(workspace->site_count, workspace->site_count,
                        &jastrow_count) ||
      !checked_multiply((size_t)workspace->ndoublon_holon_2site_idx,
                        2u * workspace->site_count, &dh2_count) ||
      !checked_multiply((size_t)workspace->ndoublon_holon_4site_idx,
                        4u * workspace->site_count, &dh4_count) ||
      !add_array_bytes(workspace->ngutzwiller_idx == 0
                           ? 0 : workspace->site_count,
                       sizeof(*workspace->gutzwiller_idx), &byte_count) ||
      !add_array_bytes(workspace->njastrow_idx == 0 ? 0 : jastrow_count,
                       sizeof(*workspace->jastrow_idx), &byte_count) ||
      !add_array_bytes(workspace->nspin_jastrow_idx == 0 ? 0 : jastrow_count,
                       sizeof(*workspace->spin_jastrow_idx), &byte_count) ||
      !add_array_bytes(dh2_count,
                       sizeof(*workspace->doublon_holon_2site_idx),
                       &byte_count) ||
      !add_array_bytes(dh4_count,
                       sizeof(*workspace->doublon_holon_4site_idx),
                       &byte_count) ||
      !add_array_bytes(workspace->nproj, sizeof(*workspace->parameters),
                       &byte_count)) {
    return 0;
  }
  workspace->binding_bytes = (unsigned char *)malloc(byte_count);
  if (workspace->binding_bytes == NULL) return 0;
  workspace->binding_byte_count = byte_count;
  metadata[0] = (uint64_t)workspace->site_count;
  metadata[1] = (uint64_t)workspace->nproj;
  metadata[2] = (uint64_t)workspace->ngutzwiller_idx;
  metadata[3] = (uint64_t)workspace->njastrow_idx;
  metadata[4] = (uint64_t)workspace->nspin_jastrow_idx;
  metadata[5] = (uint64_t)workspace->ndoublon_holon_2site_idx;
  metadata[6] = (uint64_t)workspace->ndoublon_holon_4site_idx;
  cursor = workspace->binding_bytes;
  remaining = byte_count;
  return append_bytes(&cursor, &remaining, metadata, sizeof(metadata)) &&
         append_bytes(&cursor, &remaining, workspace->gutzwiller_idx,
                      (workspace->ngutzwiller_idx == 0
                           ? 0 : workspace->site_count) *
                          sizeof(*workspace->gutzwiller_idx)) &&
         append_bytes(&cursor, &remaining, workspace->jastrow_idx,
                      (workspace->njastrow_idx == 0 ? 0 : jastrow_count) *
                          sizeof(*workspace->jastrow_idx)) &&
         append_bytes(&cursor, &remaining, workspace->spin_jastrow_idx,
                      (workspace->nspin_jastrow_idx == 0
                           ? 0 : jastrow_count) *
                          sizeof(*workspace->spin_jastrow_idx)) &&
         append_bytes(&cursor, &remaining,
                      workspace->doublon_holon_2site_idx,
                      dh2_count *
                          sizeof(*workspace->doublon_holon_2site_idx)) &&
         append_bytes(&cursor, &remaining,
                      workspace->doublon_holon_4site_idx,
                      dh4_count *
                          sizeof(*workspace->doublon_holon_4site_idx)) &&
         append_bytes(&cursor, &remaining, workspace->parameters,
                      workspace->nproj * sizeof(*workspace->parameters)) &&
         remaining == 0;
}

MVMCKrylovStatus mvmc_classic_krylov_projection_workspace_create(
    const MVMCClassicKrylovProjectionLayout *layout,
    MVMCClassicKrylovProjectionWorkspace **workspace) {
  MVMCClassicKrylovProjectionWorkspace *created;
  size_t expected = 0;
  size_t orbital_count = 0;
  size_t matrix_count = 0;
  size_t dh2_row_count = 0;
  size_t dh4_row_count = 0;
  size_t index;

  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (!expected_projection_count(layout, &expected) ||
      !checked_multiply(layout->site_count, 2u, &orbital_count) ||
      !checked_multiply(layout->site_count, layout->site_count,
                        &matrix_count) ||
      !checked_multiply(layout->site_count, 2u, &dh2_row_count) ||
      !checked_multiply(layout->site_count, 4u, &dh4_row_count) ||
      orbital_count == 0 || matrix_count == 0 || dh2_row_count == 0 ||
      dh4_row_count == 0 ||
      (layout->ngutzwiller_idx != 0 && layout->gutzwiller_idx == NULL) ||
      (expected != 0 && layout->parameters == NULL) ||
      !valid_jastrow_binding(layout->jastrow_idx, layout->site_count,
                              layout->njastrow_idx) ||
      !valid_jastrow_binding(layout->spin_jastrow_idx, layout->site_count,
                              layout->nspin_jastrow_idx) ||
      !valid_dh_binding(layout->doublon_holon_2site_idx,
                        layout->site_count,
                        layout->ndoublon_holon_2site_idx, 2u) ||
      !valid_dh_binding(layout->doublon_holon_4site_idx,
                        layout->site_count,
                        layout->ndoublon_holon_4site_idx, 4u)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (index = 0; index < layout->site_count; ++index) {
    if (layout->ngutzwiller_idx != 0 &&
        (layout->gutzwiller_idx[index] < 0 ||
         layout->gutzwiller_idx[index] >= layout->ngutzwiller_idx)) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  for (index = 0; index < expected; ++index) {
    if (!isfinite(creal(layout->parameters[index])) ||
        !isfinite(cimag(layout->parameters[index])) ||
        cimag(layout->parameters[index]) != 0.0) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  created = (MVMCClassicKrylovProjectionWorkspace *)calloc(
      1, sizeof(*created));
  if (created == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  created->site_count = layout->site_count;
  created->orbital_count = orbital_count;
  created->nproj = expected;
  created->ngutzwiller_idx = layout->ngutzwiller_idx;
  created->njastrow_idx = layout->njastrow_idx;
  created->nspin_jastrow_idx = layout->nspin_jastrow_idx;
  created->ndoublon_holon_2site_idx =
      layout->ndoublon_holon_2site_idx;
  created->ndoublon_holon_4site_idx =
      layout->ndoublon_holon_4site_idx;
  if (!allocate_and_copy(layout, created) ||
      !create_binding_record(created)) {
    mvmc_classic_krylov_projection_workspace_destroy(created);
    return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
  *workspace = created;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_classic_krylov_projection_workspace_destroy(
    MVMCClassicKrylovProjectionWorkspace *workspace) {
  if (workspace == NULL) return;
  free(workspace->binding_bytes);
  free(workspace->scratch_counts);
  free(workspace->parameters);
  free(workspace->doublon_holon_4site_idx);
  free(workspace->doublon_holon_2site_idx);
  free(workspace->spin_jastrow_idx);
  free(workspace->jastrow_idx);
  free(workspace->gutzwiller_idx);
  free(workspace);
}

size_t mvmc_classic_krylov_projection_workspace_bytes(
    const MVMCClassicKrylovProjectionWorkspace *workspace) {
  size_t bytes;
  size_t jastrow_count;
  size_t dh2_count;
  size_t dh4_count;

  if (workspace == NULL ||
      !checked_multiply(workspace->site_count, workspace->site_count,
                        &jastrow_count) ||
      !checked_multiply((size_t)workspace->ndoublon_holon_2site_idx,
                        2u * workspace->site_count, &dh2_count) ||
      !checked_multiply((size_t)workspace->ndoublon_holon_4site_idx,
                        4u * workspace->site_count, &dh4_count)) {
    return 0;
  }
  bytes = sizeof(*workspace);
  if (!add_array_bytes(workspace->ngutzwiller_idx == 0
                           ? 0 : workspace->site_count,
                       sizeof(*workspace->gutzwiller_idx), &bytes) ||
      !add_array_bytes(workspace->njastrow_idx == 0 ? 0 : jastrow_count,
                       sizeof(*workspace->jastrow_idx), &bytes) ||
      !add_array_bytes(workspace->nspin_jastrow_idx == 0 ? 0 : jastrow_count,
                       sizeof(*workspace->spin_jastrow_idx), &bytes) ||
      !add_array_bytes(dh2_count,
                       sizeof(*workspace->doublon_holon_2site_idx), &bytes) ||
      !add_array_bytes(dh4_count,
                       sizeof(*workspace->doublon_holon_4site_idx), &bytes) ||
      !add_array_bytes(workspace->nproj, sizeof(*workspace->parameters),
                       &bytes) ||
      !add_array_bytes(workspace->nproj, sizeof(*workspace->scratch_counts),
                       &bytes) ||
      !checked_add(bytes, workspace->binding_byte_count, &bytes)) {
    return 0;
  }
  return bytes;
}

const unsigned char *mvmc_classic_krylov_projection_binding_bytes(
    const MVMCClassicKrylovProjectionWorkspace *workspace,
    size_t *byte_count) {
  if (byte_count == NULL) return NULL;
  *byte_count = workspace == NULL ? 0 : workspace->binding_byte_count;
  return workspace == NULL ? NULL : workspace->binding_bytes;
}

MVMCKrylovStatus mvmc_classic_krylov_projection_evaluate(
    MVMCClassicKrylovProjectionWorkspace *workspace,
    const int *ele_num, size_t orbital_count,
    int64_t *counts, size_t count_capacity, double *log_factor) {
  const int *up;
  const int *down;
  size_t site;
  size_t other;
  size_t offset;
  size_t parameter;
  double value = 0.0;
  int shell;

  if (workspace == NULL || ele_num == NULL || log_factor == NULL ||
      orbital_count != workspace->orbital_count ||
      (counts == NULL && count_capacity != 0) ||
      (counts != NULL && count_capacity < workspace->nproj)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (site = 0; site < orbital_count; ++site) {
    if (ele_num[site] != 0 && ele_num[site] != 1) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  if (workspace->nproj != 0) {
    memset(workspace->scratch_counts, 0,
           workspace->nproj * sizeof(*workspace->scratch_counts));
  }
  up = ele_num;
  down = ele_num + workspace->site_count;
  if (workspace->ngutzwiller_idx != 0) {
    for (site = 0; site < workspace->site_count; ++site) {
      workspace->scratch_counts[workspace->gutzwiller_idx[site]] +=
          (int64_t)(up[site] * down[site]);
    }
  }
  offset = (size_t)workspace->ngutzwiller_idx;
  if (workspace->njastrow_idx != 0) {
    for (site = 0; site < workspace->site_count; ++site) {
      const int charge = up[site] + down[site] - 1;
      if (charge == 0) continue;
      for (other = site + 1u; other < workspace->site_count; ++other) {
        const int other_charge = up[other] + down[other] - 1;
        const int index = workspace->jastrow_idx[
            site * workspace->site_count + other];
        workspace->scratch_counts[offset + (size_t)index] +=
            (int64_t)(charge * other_charge);
      }
    }
  }
  offset += (size_t)workspace->njastrow_idx;
  if (workspace->nspin_jastrow_idx != 0) {
    for (site = 0; site < workspace->site_count; ++site) {
      const int spin = up[site] - down[site];
      if (spin == 0) continue;
      for (other = site + 1u; other < workspace->site_count; ++other) {
        const int other_spin = up[other] - down[other];
        const int index = workspace->spin_jastrow_idx[
            site * workspace->site_count + other];
        workspace->scratch_counts[offset + (size_t)index] +=
            (int64_t)(spin * other_spin);
      }
    }
  }
  offset += (size_t)workspace->nspin_jastrow_idx;
  for (shell = 0; shell < workspace->ndoublon_holon_2site_idx; ++shell) {
    const int *neighbors = workspace->doublon_holon_2site_idx +
        (size_t)shell * 2u * workspace->site_count;
    for (site = 0; site < workspace->site_count; ++site) {
      const int occupation = up[site] + down[site];
      const int xi = occupation / 2;
      int xm = 0;
      size_t neighbor;
      if (occupation == 1) continue;
      for (neighbor = 0; neighbor < 2u; ++neighbor) {
        const size_t neighbor_site =
            (size_t)neighbors[2u * site + neighbor];
        xm += xi == 0
                  ? up[neighbor_site] * down[neighbor_site]
                  : (1 - up[neighbor_site]) * (1 - down[neighbor_site]);
      }
      parameter = offset + (size_t)shell +
          (size_t)(xi + 2 * xm) *
              (size_t)workspace->ndoublon_holon_2site_idx;
      ++workspace->scratch_counts[parameter];
    }
  }
  offset += 6u * (size_t)workspace->ndoublon_holon_2site_idx;
  for (shell = 0; shell < workspace->ndoublon_holon_4site_idx; ++shell) {
    const int *neighbors = workspace->doublon_holon_4site_idx +
        (size_t)shell * 4u * workspace->site_count;
    for (site = 0; site < workspace->site_count; ++site) {
      const int occupation = up[site] + down[site];
      const int xi = occupation / 2;
      int xm = 0;
      size_t neighbor;
      if (occupation == 1) continue;
      for (neighbor = 0; neighbor < 4u; ++neighbor) {
        const size_t neighbor_site =
            (size_t)neighbors[4u * site + neighbor];
        xm += xi == 0
                  ? up[neighbor_site] * down[neighbor_site]
                  : (1 - up[neighbor_site]) * (1 - down[neighbor_site]);
      }
      parameter = offset + (size_t)shell +
          (size_t)(xi + 2 * xm) *
              (size_t)workspace->ndoublon_holon_4site_idx;
      ++workspace->scratch_counts[parameter];
    }
  }
  for (parameter = 0; parameter < workspace->nproj; ++parameter) {
    const double term = creal(workspace->parameters[parameter]) *
        (double)workspace->scratch_counts[parameter];
    if (!isfinite(term) || !isfinite(value + term)) {
      return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    value += term;
  }
  if (counts != NULL && workspace->nproj != 0) {
    memcpy(counts, workspace->scratch_counts,
           workspace->nproj * sizeof(*counts));
  }
  *log_factor = value;
  return MVMC_KRYLOV_STATUS_OK;
}

#endif /* reference or bounded engine */
