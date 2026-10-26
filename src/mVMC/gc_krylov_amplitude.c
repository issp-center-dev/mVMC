/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "gc_krylov_amplitude.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "gc_krylov_amplitude.c requires the bounded Krylov engine"
#endif

#include "absolute_pfaffian.h"
#include "bounded_krylov_collective.h"
#include "gc_exact_sum.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct MVMCGCKrylovAmplitudeWorkspace {
  size_t site_count;
  size_t orbital_count;
  size_t word_count;
  MVMCKrylovConfigurationDomain domain;
  int qp_total;
  int qp_start;
  int qp_end;
  size_t local_qp_count;
  size_t slater_matrix_count;
  double scaled_pivot_tolerance;
  int boundary_exact_zero;
  MVMCGCSpinSupport spin_support;
  uint64_t generation_hash;
  double complex *slater;
  double complex *global_weights;
  double complex *projection_parameters;
  size_t projection_parameter_count;
  double complex *matrix;
  int *ele_idx;
  int *ele_num;
  MVMCAbsolutePfaffianScaledValueResult *scaled_components;
  MVMCAbsolutePfaffianComplexCapacityWorkspace *pfaffian;
  MVMCClassicKrylovProjectionWorkspace *projection_workspace;
  MVMCClassicPfaffianCollectiveWorkspace *collective_workspace;
  MVMCKrylovBoundedCollectiveWorkspace *bounded_collective_workspace;
};

#if defined(MVMC_ENABLE_FAULT_INJECTION)
static MVMCGCKrylovAmplitudeFault testing_fault =
    MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONE;
static int testing_fault_rank = -1;

void mvmc_gc_krylov_amplitude_testing_set_fault(
    MVMCGCKrylovAmplitudeFault fault, int world_rank) {
  testing_fault = fault;
  testing_fault_rank = world_rank;
}

static int testing_fault_active(MVMCGCKrylovAmplitudeFault fault) {
  int rank = 0;
#ifdef _mpi_use
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  return testing_fault == fault && testing_fault_rank == rank;
}
#endif

static int checked_multiply(size_t left, size_t right, size_t *product) {
  if (product == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *product = left * right;
  return 1;
}

static uint64_t hash_bytes(uint64_t hash, const void *data,
                           size_t byte_count) {
  const unsigned char *bytes = (const unsigned char *)data;
  size_t index;
  for (index = 0; index < byte_count; ++index) {
    hash ^= (uint64_t)bytes[index];
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

static int finite_complex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static MVMCKrylovStatus map_pfaffian_status(MVMCPfaffianStatus status) {
  switch (status) {
    case MVMC_PFAFFIAN_STATUS_OK:
      return MVMC_KRYLOV_STATUS_OK;
    case MVMC_PFAFFIAN_STATUS_ALLOCATION_FAILURE:
      return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    case MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT:
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    case MVMC_PFAFFIAN_STATUS_UNSUPPORTED_FP_MODE:
    case MVMC_PFAFFIAN_STATUS_FACTORIZATION_FAILURE:
      return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
}

static MVMCKrylovStatus validate_layout(
    const MVMCGCKrylovAmplitudeLayout *layout,
    const double complex *slater, const double complex *qp_weights) {
  int index;
  if (layout == NULL || slater == NULL || qp_weights == NULL ||
      layout->site_count <= 0 || layout->site_count > INT_MAX / 2 ||
      (layout->domain != MVMC_KRYLOV_DOMAIN_GC_BALANCED &&
       layout->domain != MVMC_KRYLOV_DOMAIN_GC_EVEN) ||
      layout->qp_total <= 0 || layout->qp_start < 0 ||
      layout->qp_end < layout->qp_start ||
      layout->qp_end > layout->qp_total ||
      layout->translation_count <= 0 ||
      layout->translation_weights == NULL ||
      mvmc_gc_spin_support_validate(&layout->spin_support,
                                    layout->site_count, layout->domain) !=
          MVMC_KRYLOV_STATUS_OK ||
      !isfinite(layout->scaled_pivot_tolerance) ||
      layout->projection.site_count != (size_t)layout->site_count ||
      layout->projection.nproj < 0 ||
      (layout->projection.nproj != 0 &&
       layout->projection.parameters == NULL)) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  for (index = 0; index < layout->qp_total; ++index) {
    if (!finite_complex(qp_weights[index])) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  for (index = 0; index < layout->translation_count; ++index) {
    if (!finite_complex(layout->translation_weights[index])) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static int allocate_buffers(MVMCGCKrylovAmplitudeWorkspace *workspace,
                            const MVMCGCKrylovAmplitudeLayout *layout) {
  size_t slater_count;
  size_t bytes;
  const size_t nproj = (size_t)layout->projection.nproj;
#if defined(MVMC_ENABLE_FAULT_INJECTION)
  if (testing_fault_active(MVMC_GC_KRYLOV_AMPLITUDE_FAULT_ALLOCATION)) {
    return 0;
  }
#endif
  if (!checked_multiply((size_t)layout->qp_total,
                        workspace->slater_matrix_count, &slater_count) ||
      !checked_multiply(slater_count, sizeof(*workspace->slater), &bytes) ||
      !checked_multiply(workspace->slater_matrix_count,
                        sizeof(*workspace->matrix), &bytes) ||
      !checked_multiply((size_t)layout->qp_total,
                        sizeof(*workspace->global_weights), &bytes) ||
      !checked_multiply(workspace->orbital_count,
                        sizeof(*workspace->ele_idx), &bytes) ||
      !checked_multiply(workspace->local_qp_count,
                        sizeof(*workspace->scaled_components), &bytes)) {
    return 0;
  }
  workspace->slater =
      (double complex *)malloc(slater_count * sizeof(*workspace->slater));
  workspace->matrix = (double complex *)malloc(
      workspace->slater_matrix_count * sizeof(*workspace->matrix));
  workspace->global_weights = (double complex *)malloc(
      (size_t)layout->qp_total * sizeof(*workspace->global_weights));
  workspace->ele_idx =
      (int *)calloc(workspace->orbital_count, sizeof(*workspace->ele_idx));
  workspace->ele_num =
      (int *)calloc(workspace->orbital_count, sizeof(*workspace->ele_num));
  if (nproj != 0) {
    workspace->projection_parameters = (double complex *)calloc(
        nproj, sizeof(*workspace->projection_parameters));
  }
  if (workspace->local_qp_count != 0) {
    workspace->scaled_components =
        (MVMCAbsolutePfaffianScaledValueResult *)calloc(
            workspace->local_qp_count,
            sizeof(*workspace->scaled_components));
  }
  return workspace->slater != NULL && workspace->matrix != NULL &&
         workspace->global_weights != NULL && workspace->ele_idx != NULL &&
         workspace->ele_num != NULL &&
         (nproj == 0 || workspace->projection_parameters != NULL) &&
         (workspace->local_qp_count == 0 ||
          workspace->scaled_components != NULL);
}

/* Exact replicated-data audit on the communicator and generation hash. */
static MVMCPfaffianStatus audit_binding(
    MVMCGCKrylovAmplitudeWorkspace *workspace,
    const MVMCGCKrylovAmplitudeLayout *layout,
    const double complex *slater, const double complex *qp_weights) {
  uint64_t metadata[12];
  uint64_t tolerance_bits = 0;
  const unsigned char *projection_binding;
  size_t projection_binding_bytes = 0;
  const size_t slater_bytes = (size_t)workspace->qp_total *
                              workspace->slater_matrix_count *
                              sizeof(*slater);
  const size_t translation_bytes =
      (size_t)layout->translation_count *
      sizeof(*layout->translation_weights);
  MVMCPfaffianStatus status;
  int all_equal = 0;
  uint64_t hash;

  memcpy(&tolerance_bits, &workspace->scaled_pivot_tolerance,
         sizeof(tolerance_bits));
  metadata[0] = UINT64_C(0x4743414d504c4933); /* GC amplitude v3 */
  metadata[1] = (uint64_t)workspace->site_count;
  metadata[2] = (uint64_t)(unsigned)workspace->domain;
  metadata[3] = (uint64_t)workspace->qp_total;
  metadata[4] = (uint64_t)(unsigned)layout->translation_count;
  metadata[5] = (uint64_t)(unsigned)workspace->boundary_exact_zero;
  metadata[6] = tolerance_bits;
  /* Support descriptor schema v1 as explicit integers (no padding). */
  metadata[7] = (uint64_t)MVMC_GC_SPIN_SUPPORT_SCHEMA;
  metadata[8] = (uint64_t)(unsigned)workspace->spin_support.active;
  metadata[9] = (uint64_t)(unsigned)workspace->spin_support.gauss_count;
  metadata[10] = (uint64_t)(unsigned)workspace->spin_support.total_spin;
  metadata[11] = (uint64_t)(unsigned)workspace->spin_support.pair_symmetry;
  status = mvmc_classic_pfaffian_collective_all_equal_u64(
      workspace->collective_workspace, metadata,
      sizeof(metadata) / sizeof(metadata[0]), &all_equal);
  if (status == MVMC_PFAFFIAN_STATUS_OK && all_equal) {
    projection_binding = mvmc_classic_krylov_projection_binding_bytes(
        workspace->projection_workspace, &projection_binding_bytes);
    if (projection_binding == NULL || projection_binding_bytes == 0) {
      return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
    }
    status = mvmc_classic_pfaffian_collective_all_equal_bytes(
        workspace->collective_workspace, projection_binding,
        projection_binding_bytes, &all_equal);
  } else {
    projection_binding = NULL;
  }
  if (status == MVMC_PFAFFIAN_STATUS_OK && all_equal) {
    status = mvmc_classic_pfaffian_collective_all_equal_bytes(
        workspace->collective_workspace, layout->translation_weights,
        translation_bytes, &all_equal);
  }
  if (status == MVMC_PFAFFIAN_STATUS_OK && all_equal) {
    status = mvmc_classic_pfaffian_collective_all_equal_bytes(
        workspace->collective_workspace, qp_weights,
        (size_t)workspace->qp_total * sizeof(*qp_weights), &all_equal);
  }
  if (status == MVMC_PFAFFIAN_STATUS_OK && all_equal) {
    status = mvmc_classic_pfaffian_collective_all_equal_bytes(
        workspace->collective_workspace, slater, slater_bytes, &all_equal);
  }
  if (status != MVMC_PFAFFIAN_STATUS_OK || !all_equal) {
    return status == MVMC_PFAFFIAN_STATUS_OK
               ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT
               : status;
  }
  hash = hash_bytes(UINT64_C(1469598103934665603), metadata,
                    sizeof(metadata));
  hash = hash_bytes(hash, projection_binding, projection_binding_bytes);
  hash = hash_bytes(hash, layout->translation_weights, translation_bytes);
  hash = hash_bytes(hash, qp_weights,
                    (size_t)workspace->qp_total * sizeof(*qp_weights));
  hash = hash_bytes(hash, slater, slater_bytes);
  workspace->generation_hash = hash == 0 ? UINT64_C(1) : hash;
  return MVMC_PFAFFIAN_STATUS_OK;
}

static MVMCPfaffianStatus audit_qp_partition(
    MVMCGCKrylovAmplitudeWorkspace *workspace) {
  MVMCProjectedAmplitudeResult local;
  MVMCClassicPfaffianCollectiveResult collective;
  memset(&local, 0, sizeof(local));
  local.valid = 1;
  local.regular_count = workspace->local_qp_count;
  return mvmc_classic_pfaffian_collective_aggregate(
      workspace->collective_workspace, MVMC_PFAFFIAN_STATUS_OK, &local,
      workspace->qp_total, workspace->qp_start, workspace->qp_end, -1, -1,
      &collective);
}

static MVMCKrylovStatus synchronize(
    MVMCGCKrylovAmplitudeWorkspace *workspace, MVMCKrylovStatus local) {
  int global_status = (int)MVMC_KRYLOV_STATUS_OK;
  if (mvmc_classic_pfaffian_collective_max_int(
          workspace->collective_workspace, (int)local, &global_status) !=
      MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  return (MVMCKrylovStatus)global_status;
}

MVMCKrylovStatus mvmc_gc_krylov_amplitude_workspace_create(
    const MVMCGCKrylovAmplitudeLayout *layout,
    const double complex *slater, const double complex *qp_weights,
    MVMCGCKrylovAmplitudeWorkspace **workspace) {
  MVMCGCKrylovAmplitudeWorkspace *created = NULL;
  MVMCClassicPfaffianCollectiveWorkspace *collective = NULL;
  MVMCClassicKrylovProjectionLayout projection;
  MVMCKrylovStatus status;
  MVMCPfaffianStatus pf_status;
  int all_ready = 0;
  int boundary = 0;
  int index;

  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  *workspace = NULL;
  if (layout == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  pf_status = mvmc_classic_pfaffian_collective_workspace_create(
      layout->communicator, &collective);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK) {
    return map_pfaffian_status(pf_status);
  }
  created = (MVMCGCKrylovAmplitudeWorkspace *)calloc(1, sizeof(*created));
  pf_status = mvmc_classic_pfaffian_collective_all_true(
      collective, created != NULL, &all_ready);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK || !all_ready) {
    mvmc_classic_pfaffian_collective_workspace_destroy(collective);
    free(created);
    return pf_status == MVMC_PFAFFIAN_STATUS_OK
               ? MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE
               : map_pfaffian_status(pf_status);
  }
  created->collective_workspace = collective;

  /* Signed shapes and finite inputs first, then the exact boundary flag. */
  status = validate_layout(layout, slater, qp_weights);
  if (status == MVMC_KRYLOV_STATUS_OK &&
      GCExactComplexSumIsZero(layout->translation_weights,
                              layout->translation_count, &boundary) != 0) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = synchronize(created, status);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    mvmc_gc_krylov_amplitude_workspace_destroy(created);
    return status;
  }
  created->site_count = (size_t)layout->site_count;
  created->orbital_count = 2 * (size_t)layout->site_count;
  created->word_count = mvmc_krylov_fock_word_count(created->site_count);
  created->domain = layout->domain;
  created->qp_total = layout->qp_total;
  created->qp_start = layout->qp_start;
  created->qp_end = layout->qp_end;
  created->local_qp_count = (size_t)(layout->qp_end - layout->qp_start);
  created->scaled_pivot_tolerance = layout->scaled_pivot_tolerance;
  created->boundary_exact_zero = boundary;
  created->spin_support = layout->spin_support;
  created->projection_parameter_count = (size_t)layout->projection.nproj;
  status = checked_multiply(created->orbital_count, created->orbital_count,
                            &created->slater_matrix_count) &&
                   created->word_count != 0 &&
                   allocate_buffers(created, layout)
               ? MVMC_KRYLOV_STATUS_OK
               : MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = map_pfaffian_status(
        mvmc_absolute_pfaffian_complex_capacity_workspace_create(
            (int)created->orbital_count, &created->pfaffian));
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    /* The correlation factor uses Re(Proj), as the GC LogProjVal does. */
    for (index = 0; index < layout->projection.nproj; ++index) {
      created->projection_parameters[index] =
          creal(layout->projection.parameters[index]) + 0.0 * I;
    }
    projection = layout->projection;
    projection.parameters = created->projection_parameters;
    status = mvmc_classic_krylov_projection_workspace_create(
        &projection, &created->projection_workspace);
  }
  status = synchronize(created, status);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_collective_workspace_create(
        layout->communicator, (size_t)created->qp_total, 1,
        &created->bounded_collective_workspace);
    status = synchronize(created, status);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    mvmc_gc_krylov_amplitude_workspace_destroy(created);
    return status;
  }
  pf_status = audit_binding(created, layout, slater, qp_weights);
  if (pf_status == MVMC_PFAFFIAN_STATUS_OK) {
    pf_status = audit_qp_partition(created);
  }
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK) {
    mvmc_gc_krylov_amplitude_workspace_destroy(created);
    return pf_status == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT
               ? MVMC_KRYLOV_STATUS_INVALID_ARGUMENT
               : MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  memcpy(created->slater, slater,
         (size_t)created->qp_total * created->slater_matrix_count *
             sizeof(*created->slater));
  memcpy(created->global_weights, qp_weights,
         (size_t)created->qp_total * sizeof(*created->global_weights));
  *workspace = created;
  return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_gc_krylov_amplitude_workspace_destroy(
    MVMCGCKrylovAmplitudeWorkspace *workspace) {
  if (workspace == NULL) return;
  mvmc_bounded_krylov_collective_workspace_destroy(
      workspace->bounded_collective_workspace);
  mvmc_classic_krylov_projection_workspace_destroy(
      workspace->projection_workspace);
  mvmc_absolute_pfaffian_complex_capacity_workspace_destroy(
      workspace->pfaffian);
  mvmc_classic_pfaffian_collective_workspace_destroy(
      workspace->collective_workspace);
  free(workspace->scaled_components);
  free(workspace->ele_num);
  free(workspace->ele_idx);
  free(workspace->matrix);
  free(workspace->projection_parameters);
  free(workspace->global_weights);
  free(workspace->slater);
  free(workspace);
}

size_t mvmc_gc_krylov_amplitude_workspace_bytes(
    const MVMCGCKrylovAmplitudeWorkspace *workspace) {
  if (workspace == NULL) return 0;
  return sizeof(*workspace) +
         (size_t)workspace->qp_total * workspace->slater_matrix_count *
             sizeof(*workspace->slater) +
         workspace->slater_matrix_count * sizeof(*workspace->matrix) +
         (size_t)workspace->qp_total * sizeof(*workspace->global_weights) +
         2 * workspace->orbital_count * sizeof(*workspace->ele_idx) +
         workspace->local_qp_count * sizeof(*workspace->scaled_components) +
         mvmc_absolute_pfaffian_complex_capacity_workspace_bytes(
             workspace->pfaffian) +
         mvmc_classic_krylov_projection_workspace_bytes(
             workspace->projection_workspace) +
         mvmc_bounded_krylov_collective_workspace_bytes(
             workspace->bounded_collective_workspace);
}

uint64_t mvmc_gc_krylov_amplitude_generation_hash(
    const MVMCGCKrylovAmplitudeWorkspace *workspace) {
  return workspace == NULL ? 0 : workspace->generation_hash;
}

/* Ascending occupied orbitals and occupation numbers of one configuration. */
static int decode_configuration(MVMCGCKrylovAmplitudeWorkspace *workspace,
                                const uint64_t *words, size_t word_count,
                                int *electron_count) {
  const size_t used_bits = workspace->orbital_count % 64;
  size_t orbital;
  int up = 0;
  int down = 0;
  int count = 0;
  if (words == NULL || word_count != workspace->word_count) return 0;
  if (used_bits != 0 &&
      (words[word_count - 1] & ~((UINT64_C(1) << used_bits) - 1)) != 0) {
    return 0;
  }
  for (orbital = 0; orbital < workspace->orbital_count; ++orbital) {
    const int occupied = (int)((words[orbital / 64] >> (orbital % 64)) & 1);
    workspace->ele_num[orbital] = occupied;
    if (occupied) {
      workspace->ele_idx[count++] = (int)orbital;
      if (orbital < workspace->site_count) {
        ++up;
      } else {
        ++down;
      }
    }
  }
  if ((workspace->domain == MVMC_KRYLOV_DOMAIN_GC_BALANCED && up != down) ||
      (count % 2) != 0) {
    return 0;
  }
  *electron_count = count;
  return 1;
}

/* Spin support of the configuration last decoded into ele_num. */
static MVMCKrylovStatus support_of(
    const MVMCGCKrylovAmplitudeWorkspace *workspace,
    MVMCGCSpinSupportResult *support) {
#if defined(MVMC_ENABLE_FAULT_INJECTION)
  if (testing_fault_active(MVMC_GC_KRYLOV_AMPLITUDE_FAULT_SUPPORT))
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
#endif
  return mvmc_gc_spin_support_evaluate(
      &workspace->spin_support, (int)workspace->site_count, workspace->domain,
      workspace->ele_num, workspace->orbital_count, support);
}

/* X[r,c] = -Slater[x_r, x_c] (column-major), as GCComponentBuildX. */
static void build_component_matrix(
    const MVMCGCKrylovAmplitudeWorkspace *workspace, int global_qp,
    int electron_count) {
  const double complex *slater =
      workspace->slater +
      (size_t)global_qp * workspace->slater_matrix_count;
  const size_t stride = workspace->orbital_count;
  int row, column;
  for (column = 0; column < electron_count; ++column) {
    for (row = 0; row < electron_count; ++row) {
      workspace->matrix[row + (size_t)column * (size_t)electron_count] =
          -slater[(size_t)workspace->ele_idx[row] * stride +
                  (size_t)workspace->ele_idx[column]];
    }
  }
}

static MVMCKrylovStatus apply_projection(
    double log_projection, MVMCKrylovScaledAmplitudeResult *amplitude) {
  MVMCScaledComplex factor;
  MVMCScaledComplex projected;
  if (log_projection == 0.0) return MVMC_KRYLOV_STATUS_OK;
  if (mvmc_scaled_complex_make_finite(1.0 + 0.0 * I, log_projection,
                                      -INFINITY, &factor) !=
          MVMC_PFAFFIAN_STATUS_OK ||
      mvmc_scaled_complex_multiply(&factor, &amplitude->value, &projected) !=
          MVMC_PFAFFIAN_STATUS_OK ||
      projected.state == MVMC_SCALED_COMPLEX_NONFINITE) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  amplitude->value = projected;
  return MVMC_KRYLOV_STATUS_OK;
}

#if defined(MVMC_HAVE_MPFR_REFINEMENT)
static MVMCKrylovStatus refined_base(
    const uint64_t *words, size_t word_count, void *context,
    MVMCRefinedBase *result) {
  MVMCGCKrylovAmplitudeWorkspace *w = context;
  MVMCGCSpinSupportResult support;
  int n = 0;
  int64_t *counts = NULL;
  double log_projection = 0.0;
  int support_zero;
  MVMCKrylovStatus status;
  if (!decode_configuration(w, words, word_count, &n))
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  /* Same pure helper as the ordinary callback.  No MPI here: the wrapper's
   * synchronize after the engine returns propagates this failure. */
  if (support_of(w, &support) != MVMC_KRYLOV_STATUS_OK)
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  support_zero =
      (w->boundary_exact_zero && (n == 0 || (size_t)n == w->orbital_count)) ||
      support.zero_reasons != 0;
#if defined(MVMC_ENABLE_FAULT_INJECTION)
  if (testing_fault_active(MVMC_GC_KRYLOV_AMPLITUDE_FAULT_REFINED_UNMASKED))
    support_zero =
        w->boundary_exact_zero && (n == 0 || (size_t)n == w->orbital_count);
#endif
  if (w->projection_parameter_count) {
    counts = calloc(w->projection_parameter_count, sizeof(*counts));
    if (counts == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
  }
  status = mvmc_classic_krylov_projection_evaluate(
      w->projection_workspace, w->ele_num, w->orbital_count, counts,
      w->projection_parameter_count, &log_projection);
  if (status == MVMC_KRYLOV_STATUS_OK)
    status = mvmc_refined_base_build(
        result, w->slater, w->orbital_count, w->ele_idx, (size_t)n,
        w->global_weights, (size_t)w->qp_total, support_zero,
        w->projection_parameters, counts, w->projection_parameter_count);
  free(counts);
  return status;
}

MVMCKrylovStatus mvmc_gc_krylov_refine_final(
    MVMCGCKrylovAmplitudeWorkspace *w, const MVMCKrylovFockModel *model,
    const MVMCKrylovFinalStatePolicy *policy, const uint64_t *current,
    const uint64_t *proposal, size_t word_count, double uniform,
    MVMCKrylovRefinementResult *result) {
  MVMCKrylovRefinementResult candidate;
  MVMCKrylovFinalStateSnapshot snapshot;
  MVMCKrylovStatus status;
  uint64_t agreement[5];
  int equal = 0;
  if (w == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  status = synchronize(w, result != NULL && model != NULL && policy != NULL &&
      current != NULL && proposal != NULL && word_count == w->word_count
      ? MVMC_KRYLOV_STATUS_OK : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (mvmc_classic_pfaffian_collective_all_equal_u64(
          w->collective_workspace, current, word_count, &equal) !=
          MVMC_PFAFFIAN_STATUS_OK || !equal ||
      mvmc_classic_pfaffian_collective_all_equal_u64(
          w->collective_workspace, proposal, word_count, &equal) !=
          MVMC_PFAFFIAN_STATUS_OK || !equal)
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  agreement[0] = policy->policy_hash;
  memcpy(agreement + 1, &uniform, sizeof(uniform));
  if (mvmc_classic_pfaffian_collective_all_equal_u64(
          w->collective_workspace, agreement, 2, &equal) !=
          MVMC_PFAFFIAN_STATUS_OK || !equal)
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  memset(&candidate, 0, sizeof(candidate));
  status = mvmc_krylov_refine_final(model, policy, current, proposal,
      word_count, uniform, refined_base, w, NULL, &candidate);
  status = synchronize(w, status);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  candidate.krylov.statistics.amplitude_generation_hash = w->generation_hash;
  agreement[0] = (uint64_t)candidate.accepted;
  agreement[1] = (uint64_t)candidate.exact_zero_proposal;
  agreement[2] = (uint64_t)candidate.precision;
  agreement[3] = candidate.node_calls;
  agreement[4] = 0;
  if (candidate.accepted) {
    status = mvmc_krylov_final_state_snapshot_from_evaluation(policy, proposal,
        word_count, &candidate.krylov, &candidate.evaluation, &snapshot);
    if (status == MVMC_KRYLOV_STATUS_OK) agreement[4] = snapshot.integrity_hash;
  }
  status = synchronize(w, status);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (mvmc_classic_pfaffian_collective_all_equal_u64(
          w->collective_workspace, agreement, 5, &equal) !=
          MVMC_PFAFFIAN_STATUS_OK || !equal)
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  *result = candidate;
  return MVMC_KRYLOV_STATUS_OK;
}
#endif

MVMCKrylovStatus mvmc_gc_krylov_scaled_amplitude(
    const uint64_t *words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
  MVMCGCKrylovAmplitudeWorkspace *workspace =
      (MVMCGCKrylovAmplitudeWorkspace *)context;
  MVMCKrylovScaledAmplitudeResult projected;
  MVMCKrylovBoundedCollectiveResult collective;
  MVMCGCSpinSupportResult support;
  MVMCKrylovStatus status;
  MVMCPfaffianStatus pf_status;
  double log_projection = 0.0;
  int electron_count = 0;
  int support_zero = 0;
  int all_equal = 0;
  int all_valid = 0;
  size_t local_qp;

  if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  memset(result, 0, sizeof(*result));
  if (workspace == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;

  /* Every rank evaluates the same configuration of the domain. */
  pf_status = mvmc_classic_pfaffian_collective_all_true(
      workspace->collective_workspace,
      words != NULL && word_count == workspace->word_count, &all_valid);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (!all_valid) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  pf_status = mvmc_classic_pfaffian_collective_all_equal_u64(
      workspace->collective_workspace, words, word_count, &all_equal);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK || !all_equal) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  pf_status = mvmc_classic_pfaffian_collective_all_true(
      workspace->collective_workspace,
      decode_configuration(workspace, words, word_count, &electron_count),
      &all_valid);
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  if (!all_valid) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;

  /* Support of this configuration, then the correlation factor from the
   * private occupation numbers.  A helper failure is an invariant failure
   * agreed below; no collective is skipped. */
  memset(&support, 0, sizeof(support));
  status = support_of(workspace, &support) == MVMC_KRYLOV_STATUS_OK
               ? MVMC_KRYLOV_STATUS_OK
               : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_classic_krylov_projection_evaluate(
        workspace->projection_workspace, workspace->ele_num,
        workspace->orbital_count, NULL, 0, &log_projection);
  }
  if (status == MVMC_KRYLOV_STATUS_OK && !isfinite(log_projection)) {
    status = MVMC_KRYLOV_STATUS_NONFINITE;
  }
  status = synchronize(workspace, status);
  support_zero = support.zero_reasons != 0;
#if defined(MVMC_ENABLE_FAULT_INJECTION)
  if (testing_fault_active(MVMC_GC_KRYLOV_AMPLITUDE_FAULT_ORDINARY_UNMASKED))
    support_zero = 0;
#endif

  /* Local QP components (-1)^(N/2) Pf(-Slater[x,x]) with capacity scratch. */
  for (local_qp = 0; status == MVMC_KRYLOV_STATUS_OK &&
                     local_qp < workspace->local_qp_count; ++local_qp) {
    MVMCAbsolutePfaffianScaledValueResult *component =
        workspace->scaled_components + local_qp;
    build_component_matrix(workspace, workspace->qp_start + (int)local_qp,
                           electron_count);
    pf_status = mvmc_absolute_pfaffian_complex_scaled_value_with_capacity(
        workspace->pfaffian, electron_count == 0 ? NULL : workspace->matrix,
        electron_count, electron_count, workspace->scaled_pivot_tolerance,
        component);
    if (pf_status != MVMC_PFAFFIAN_STATUS_OK) {
      status = map_pfaffian_status(pf_status);
    } else if (component->factor_state == MVMC_PFAFFIAN_VALUE_NONFINITE) {
      status = MVMC_KRYLOV_STATUS_NONFINITE;
    } else if (component->factor_state == MVMC_PFAFFIAN_VALUE_INVALID) {
      status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    } else if ((electron_count / 2) % 2 != 0 &&
               component->value.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
      component->value.phase = -component->value.phase;
    }
  }
#if defined(MVMC_ENABLE_FAULT_INJECTION)
  if (status == MVMC_KRYLOV_STATUS_OK &&
      testing_fault_active(MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONFINITE)) {
    status = MVMC_KRYLOV_STATUS_NONFINITE;
  }
#endif
  status = mvmc_bounded_krylov_collective_gather_projected_amplitude(
      workspace->bounded_collective_workspace, status,
      workspace->scaled_components, workspace->local_qp_count,
      workspace->global_weights, workspace->qp_total, workspace->qp_start,
      workspace->qp_end, &projected, &collective);
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (projected.value.state == MVMC_SCALED_COMPLEX_NONFINITE) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  /* Proven symmetry zeros, only after all ranks checked their components.
   * Do not convert arbitrary numerical zeros or hide nonfinite failures. */
  if ((workspace->boundary_exact_zero &&
      (electron_count == 0 ||
       (size_t)electron_count == workspace->orbital_count)) ||
      support_zero) {
    if (mvmc_scaled_complex_make_exact_zero(&projected.value) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
  } else {
    status = apply_projection(log_projection, &projected);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
  }
  *result = projected;
  return MVMC_KRYLOV_STATUS_OK;
}
