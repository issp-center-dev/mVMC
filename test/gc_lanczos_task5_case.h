/*
 * Test-only builders for the frozen Task 5 fixtures (gc_lanczos_task5_fixture.h):
 * exact binary64 copies of the saved Slater/QP weights, the Fock model of
 * the listed Hamiltonian, production-shaped Transfer rows (H = -sum t c^+ c)
 * for the GC bridge, the spin support descriptor and the snapshot
 * comparison rule of the Task 5 specification (k * eps phase/log).
 */
#ifndef GC_LANCZOS_TASK5_CASE_H
#define GC_LANCZOS_TASK5_CASE_H

#include "gc_lanczos_task5_fixture.h"
#include "gc_krylov_amplitude.h"
#include "power_lanczos_gc_bridge.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  const Task5Case *source;
  int site_count;
  size_t orbital_count;
  double complex *slater;
  double complex *weight;
  double complex translation[1];
  int *gutzwiller_idx;
  double complex projection[1];
  int nproj;
  double complex coefficient[2];
  MVMCGCSpinSupport support;
  MVMCKrylovHamiltonianTerm *terms;
  MVMCKrylovFermionOperator *operators;
  MVMCKrylovFockModel model;
  int (*transfer_rows)[4];
  int **transfer_index;
  double complex *transfer_value;
} Task5Built;

/* Exact copy: double complex has the layout of double[2]. */
static inline double complex task5_complex(const double value[2]) {
  double complex z;
  memcpy(&z, value, sizeof(z));
  return z;
}

static inline void task5_free(Task5Built *built) {
  if (built == NULL) return;
  free(built->slater);
  free(built->weight);
  free(built->gutzwiller_idx);
  free(built->terms);
  free(built->operators);
  free(built->transfer_rows);
  free(built->transfer_index);
  free(built->transfer_value);
  memset(built, 0, sizeof(*built));
}

static inline int task5_build(const Task5Case *source, Task5Built *built) {
  const size_t norb = 2 * (size_t)source->site_count;
  const size_t matrix = norb * norb;
  size_t k;
  memset(built, 0, sizeof(*built));
  built->source = source;
  built->site_count = source->site_count;
  built->orbital_count = norb;
  built->slater = (double complex *)calloc(
      (size_t)source->qp_total * matrix, sizeof(*built->slater));
  built->weight = (double complex *)calloc((size_t)source->qp_total,
                                           sizeof(*built->weight));
  built->gutzwiller_idx =
      (int *)calloc((size_t)source->site_count, sizeof(int));
  built->terms = (MVMCKrylovHamiltonianTerm *)calloc(
      source->term_count, sizeof(*built->terms));
  built->operators = (MVMCKrylovFermionOperator *)calloc(
      2 * source->term_count, sizeof(*built->operators));
  built->transfer_rows = (int(*)[4])calloc(source->term_count,
                                           sizeof(*built->transfer_rows));
  built->transfer_index =
      (int **)calloc(source->term_count, sizeof(*built->transfer_index));
  built->transfer_value = (double complex *)calloc(
      source->term_count, sizeof(*built->transfer_value));
  if (built->slater == NULL || built->weight == NULL ||
      built->gutzwiller_idx == NULL || built->terms == NULL ||
      built->operators == NULL || built->transfer_rows == NULL ||
      built->transfer_index == NULL || built->transfer_value == NULL) {
    task5_free(built);
    return 0;
  }
  for (k = 0; k < (size_t)source->qp_total * matrix; ++k) {
    built->slater[k] = task5_complex(source->slater[k]);
  }
  for (k = 0; k < (size_t)source->qp_total; ++k) {
    built->weight[k] = task5_complex(source->weight[k]);
  }
  built->translation[0] = 1.0;
  built->projection[0] = task5_complex(source->gutzwiller);
  built->nproj = built->projection[0] != 0.0 ? 1 : 0;
  built->coefficient[0] = task5_complex(source->coefficient[0]);
  built->coefficient[1] = task5_complex(source->coefficient[1]);
  built->support.active = 1;
  built->support.gauss_count = source->gauss_count;
  built->support.total_spin = source->total_spin;
  built->support.pair_symmetry = (MVMCGCPairSymmetry)source->symmetry;
  for (k = 0; k < source->term_count; ++k) {
    const Task5Term *term = source->terms + k;
    const int out = term->op[0][0], in = term->op[1][0];
    built->terms[k].coefficient = task5_complex(term->coefficient);
    built->terms[k].operator_offset = 2 * k;
    built->terms[k].operator_count = 2;
    built->terms[k].source_kind = 0;
    built->terms[k].source_index = k;
    built->operators[2 * k].kind = term->op[0][1]
                                       ? MVMC_KRYLOV_FERMION_CREATE
                                       : MVMC_KRYLOV_FERMION_ANNIHILATE;
    built->operators[2 * k].orbital = (size_t)out;
    built->operators[2 * k + 1].kind = term->op[1][1]
                                           ? MVMC_KRYLOV_FERMION_CREATE
                                           : MVMC_KRYLOV_FERMION_ANNIHILATE;
    built->operators[2 * k + 1].orbital = (size_t)in;
    /* Only c^+_out c_in rows exist in the fixtures. */
    built->transfer_rows[k][0] = out % source->site_count;
    built->transfer_rows[k][1] = out / source->site_count;
    built->transfer_rows[k][2] = in % source->site_count;
    built->transfer_rows[k][3] = in / source->site_count;
    built->transfer_index[k] = built->transfer_rows[k];
    built->transfer_value[k] = -built->terms[k].coefficient;
  }
  built->model.site_count = (size_t)source->site_count;
  built->model.hermitian = 1;
  built->model.terms = built->terms;
  built->model.term_count = source->term_count;
  built->model.operators = built->operators;
  built->model.operator_count = 2 * source->term_count;
  built->model.domain = MVMC_KRYLOV_DOMAIN_GC_BALANCED;
  return 1;
}

/* Amplitude layout on comm with the QP partition of (rank, size). */
static inline void task5_layout(const Task5Built *built,
                         MVMCClassicPfaffianCommunicator comm, int rank,
                         int size, MVMCGCKrylovAmplitudeLayout *layout) {
  const int qp = built->source->qp_total;
  memset(layout, 0, sizeof(*layout));
  layout->site_count = built->site_count;
  layout->domain = MVMC_KRYLOV_DOMAIN_GC_BALANCED;
  layout->qp_total = qp;
  layout->qp_start = qp * rank / size;
  layout->qp_end = qp * (rank + 1) / size;
  layout->translation_count = 1;
  layout->translation_weights = built->translation;
  layout->spin_support = built->support;
  layout->scaled_pivot_tolerance = 0.0;
  layout->projection.site_count = (size_t)built->site_count;
  if (built->nproj) {
    layout->projection.nproj = 1;
    layout->projection.ngutzwiller_idx = 1;
    layout->projection.gutzwiller_idx = built->gutzwiller_idx;
    layout->projection.parameters = built->projection;
  }
  layout->communicator = comm;
}

/* Production-shaped GC view; model rows on world rank 0 only. */
static inline void task5_view(const Task5Built *built, int world_rank,
                       int chain_rank, int chain_size,
                       MVMCPowerLanczosGCView *view) {
  const int qp = built->source->qp_total;
  memset(view, 0, sizeof(*view));
  view->data.site_count = built->site_count;
  view->data.arithmetic = MVMC_POWER_LANCZOS_CLASSIC_COMPLEX;
  view->data.transfer_count = (int)built->source->term_count;
  if (world_rank == 0) {
    view->data.transfer_indices = built->transfer_index;
    view->data.transfer_parameters = built->transfer_value;
  }
  view->data.qp_total = qp;
  view->data.qp_start = qp * chain_rank / chain_size;
  view->data.qp_end = qp * (chain_rank + 1) / chain_size;
  if (built->nproj) {
    view->data.nproj = 1;
    view->data.ngutzwiller_idx = 1;
    view->data.gutzwiller_idx = built->gutzwiller_idx;
    view->data.projection_parameters = built->projection;
  }
  view->data.qp_weights = built->weight;
  view->data.slater_complex = built->slater;
  view->domain = MVMC_KRYLOV_DOMAIN_GC_BALANCED;
  view->translation_count = 1;
  view->translation_weights = built->translation;
  view->spin_support = built->support;
}

/* Specification rule: state exact; |phase| parts <= k eps and
 * |log_abs| <= k eps (1 + |ref|), k = 16 refined, 64 ordinary. */
static inline int task5_scaled_matches(const MVMCScaledComplex *value,
                                const Task5Scaled *reference, double k) {
  const double eps = 0x1p-52;
  if (!reference->finite) {
    return value->state == MVMC_SCALED_COMPLEX_EXACT_ZERO;
  }
  return value->state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
         fabs(creal(value->phase) - reference->phase[0]) <= k * eps &&
         fabs(cimag(value->phase) - reference->phase[1]) <= k * eps &&
         fabs(value->log_abs - reference->log_abs) <=
             k * eps * (1.0 + fabs(reference->log_abs));
}

static inline const Task5Scaled *task5_reference(const Task5Scaled *table,
                                          size_t count, unsigned mask,
                                          int value) {
  size_t k;
  for (k = 0; k < count; ++k) {
    if (table[k].configuration == mask && table[k].value == value) {
      return table + k;
    }
  }
  return NULL;
}

#endif /* GC_LANCZOS_TASK5_CASE_H */
