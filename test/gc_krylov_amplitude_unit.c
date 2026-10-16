/*
 * Grand-canonical absolute amplitude and bridge contract: the production sign
 * (-1)^(N/2) Pf(-Slater[x,x]) equals Pf(Slater[x,x]) over ascending orbitals
 * (independent recursive Pfaffian), Re(Proj) correlation factor, Phase 3b
 * boundary zero, node configurations, eleIdx permutation parity, caller
 * buffer invariance, input validation, MPI QP partitions with empty ranks,
 * rank-1 faults and the production-shaped GC bridge.
 */
#include "bounded_krylov_engine.h"
#include "gc_krylov_amplitude.h"
#include "power_lanczos_gc_bridge.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failure_count = 0;

#define CHECK(condition, message)                                    \
  do {                                                               \
    if (!(condition)) {                                              \
      fprintf(stderr, "FAIL: %s (line %d)\n", (message), __LINE__); \
      ++failure_count;                                               \
    }                                                                \
  } while (0)

#define L2 2
#define ORB2 4
#define MAX_ORB 8

static int world_rank(void) {
  int rank = 0;
#ifdef _mpi_use
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  return rank;
}

static int world_size(void) {
  int size = 1;
#ifdef _mpi_use
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  return size;
}

static MVMCClassicPfaffianCommunicator world_comm(void) {
#ifdef _mpi_use
  return MPI_COMM_WORLD;
#else
  return 0;
#endif
}

static MVMCClassicPfaffianCommunicator self_comm(void) {
#ifdef _mpi_use
  return MPI_COMM_SELF;
#else
  return 0;
#endif
}

/* Independent recursive Pfaffian (expansion along the first row). */
static double complex pfaffian(const double complex *a, int n, int lda) {
  double complex total = 0.0;
  double complex minor[MAX_ORB * MAX_ORB];
  int column;
  if (n == 0) return 1.0;
  if (n % 2) return 0.0;
  for (column = 1; column < n; ++column) {
    int keep[MAX_ORB];
    int count = 0;
    int r, c, index;
    for (index = 1; index < n; ++index) {
      if (index != column) keep[count++] = index;
    }
    for (r = 0; r < count; ++r) {
      for (c = 0; c < count; ++c) {
        minor[r * MAX_ORB + c] = a[keep[r] * lda + keep[c]];
      }
    }
    total += ((column % 2) ? 1.0 : -1.0) * a[0 * lda + column] *
             pfaffian(minor, count, MAX_ORB);
  }
  return total;
}

/* Row-major antisymmetric Slater with the anti-parallel pair F (L=2). */
static void antiparallel_slater(const double complex f[L2][L2],
                                double complex *slater) {
  int i, j;
  memset(slater, 0, sizeof(double complex) * ORB2 * ORB2);
  for (i = 0; i < L2; ++i) {
    for (j = 0; j < L2; ++j) {
      slater[i * ORB2 + (j + L2)] = f[i][j];
      slater[(j + L2) * ORB2 + i] = -f[i][j];
    }
  }
}

static double complex slater_amplitude(const double complex *slater,
                                       int orbital_count, uint64_t mask) {
  double complex sub[MAX_ORB * MAX_ORB];
  int occupied[MAX_ORB];
  int count = 0;
  int orbital, r, c;
  for (orbital = 0; orbital < orbital_count; ++orbital) {
    if ((mask >> orbital) & 1) occupied[count++] = orbital;
  }
  for (r = 0; r < count; ++r) {
    for (c = 0; c < count; ++c) {
      sub[r * MAX_ORB + c] =
          slater[occupied[r] * orbital_count + occupied[c]];
    }
  }
  return pfaffian(sub, count, MAX_ORB);
}

static int popcount64(uint64_t value) {
  int count = 0;
  while (value) {
    count += (int)(value & 1);
    value >>= 1;
  }
  return count;
}

static int balanced_mask(uint64_t mask, int site_count) {
  const uint64_t low = (UINT64_C(1) << site_count) - 1;
  return popcount64(mask & low) == popcount64(mask >> site_count);
}

typedef struct {
  MVMCGCKrylovAmplitudeLayout layout;
  int gutzwiller_idx[L2];
  int jastrow_rows[L2][L2];
  const int *jastrow_idx[L2];
  double complex projection[2];
} AmplitudeFixture;

static void fixture_layout(AmplitudeFixture *fixture,
                           MVMCKrylovConfigurationDomain domain,
                           int qp_total, const double complex *translation,
                           int translation_count, int with_projection,
                           MVMCClassicPfaffianCommunicator comm,
                           int partitioned) {
  int rank = 0;
  int size = 1;
#ifdef _mpi_use
  if (partitioned) {
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
  }
#else
  (void)partitioned;
#endif
  memset(fixture, 0, sizeof(*fixture));
  fixture->layout.site_count = L2;
  fixture->layout.domain = domain;
  fixture->layout.qp_total = qp_total;
  fixture->layout.qp_start = qp_total * rank / size;
  fixture->layout.qp_end = qp_total * (rank + 1) / size;
  fixture->layout.translation_count = translation_count;
  fixture->layout.translation_weights = translation;
  fixture->layout.scaled_pivot_tolerance = 0.0;
  fixture->layout.projection.site_count = L2;
  fixture->layout.communicator = comm;
  if (with_projection) {
    fixture->gutzwiller_idx[0] = 0;
    fixture->gutzwiller_idx[1] = 0;
    fixture->jastrow_rows[0][0] = -1;
    fixture->jastrow_rows[0][1] = 0;
    fixture->jastrow_rows[1][0] = 0;
    fixture->jastrow_rows[1][1] = -1;
    fixture->jastrow_idx[0] = fixture->jastrow_rows[0];
    fixture->jastrow_idx[1] = fixture->jastrow_rows[1];
    /* Imaginary parts must not enter the GC correlation factor. */
    fixture->projection[0] = -0.375 + 0.75 * I;
    fixture->projection[1] = 0.1875 - 0.5 * I;
    fixture->layout.projection.nproj = 2;
    fixture->layout.projection.ngutzwiller_idx = 1;
    fixture->layout.projection.njastrow_idx = 1;
    fixture->layout.projection.gutzwiller_idx = fixture->gutzwiller_idx;
    fixture->layout.projection.jastrow_idx = fixture->jastrow_idx;
    fixture->layout.projection.parameters = fixture->projection;
  }
}

static double reference_log_projection(uint64_t mask) {
  /* Gutzwiller n_up n_dn (idx 0) and Jastrow (n0-1)(n1-1) (idx 1). */
  int n_up[L2], n_dn[L2], site, count_g = 0, count_j;
  for (site = 0; site < L2; ++site) {
    n_up[site] = (int)((mask >> site) & 1);
    n_dn[site] = (int)((mask >> (site + L2)) & 1);
    count_g += n_up[site] * n_dn[site];
  }
  count_j = (n_up[0] + n_dn[0] - 1) * (n_up[1] + n_dn[1] - 1);
  return -0.375 * count_g + 0.1875 * count_j;
}

static MVMCKrylovStatus evaluate(MVMCGCKrylovAmplitudeWorkspace *workspace,
                                 uint64_t mask,
                                 MVMCKrylovScaledAmplitudeResult *result) {
  uint64_t words = mask;
  return mvmc_gc_krylov_scaled_amplitude(&words, 1, workspace, result);
}

static int matches(const MVMCKrylovScaledAmplitudeResult *result,
                   double complex expected) {
  if (expected == 0.0) {
    return result->value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO;
  }
  return result->value.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
         cabs(result->value.phase - expected / cabs(expected)) <= 3e-12 &&
         fabs(result->value.log_abs - log(cabs(expected))) <= 3e-11;
}

static const double complex kF2[L2][L2] = {
    {1.0, (1.0 + 2.0 * I) / 8.0}, {(3.0 - 1.0 * I) / 16.0, 0.75}};

static void test_single_qp_amplitudes(void) {
  double complex slater[ORB2 * ORB2];
  double complex slater_copy[ORB2 * ORB2];
  const double complex translation[1] = {1.0};
  const double complex weight[1] = {1.0};
  const MVMCKrylovConfigurationDomain domains[2] = {
      MVMC_KRYLOV_DOMAIN_GC_BALANCED, MVMC_KRYLOV_DOMAIN_GC_EVEN};
  int d, with_projection;
  antiparallel_slater(kF2, slater);
  /* General pairing: add same-spin components (G2 fixture). */
  for (d = 0; d < 2; ++d) {
    for (with_projection = 0; with_projection < 2; ++with_projection) {
      AmplitudeFixture fixture;
      AmplitudeFixture fixture_copy;
      MVMCGCKrylovAmplitudeWorkspace *workspace = NULL;
      uint64_t mask;
      if (domains[d] == MVMC_KRYLOV_DOMAIN_GC_EVEN) {
        slater[0 * ORB2 + 1] = (1.0 + 1.0 * I) / 4.0;
        slater[1 * ORB2 + 0] = -(1.0 + 1.0 * I) / 4.0;
        slater[2 * ORB2 + 3] = (2.0 - 1.0 * I) / 8.0;
        slater[3 * ORB2 + 2] = -(2.0 - 1.0 * I) / 8.0;
      }
      memcpy(slater_copy, slater, sizeof(slater));
      fixture_layout(&fixture, domains[d], 1, translation, 1,
                     with_projection, world_comm(), 1);
      fixture_copy = fixture;
      {
        const MVMCKrylovStatus created =
            mvmc_gc_krylov_amplitude_workspace_create(
                &fixture.layout, slater, weight, &workspace);
        if (created != MVMC_KRYLOV_STATUS_OK) {
          fprintf(stderr, "create domain=%d proj=%d status=%s\n",
                  (int)domains[d], with_projection,
                  mvmc_krylov_status_string(created));
        }
        CHECK(created == MVMC_KRYLOV_STATUS_OK && workspace != NULL,
              "single-QP GC amplitude workspace");
      }
      if (workspace == NULL) continue;
      CHECK(mvmc_gc_krylov_amplitude_generation_hash(workspace) != 0 &&
                mvmc_gc_krylov_amplitude_workspace_bytes(workspace) > 0,
            "generation hash and bytes");
      for (mask = 0; mask < 16; ++mask) {
        MVMCKrylovScaledAmplitudeResult result;
        const int allowed =
            domains[d] == MVMC_KRYLOV_DOMAIN_GC_BALANCED
                ? balanced_mask(mask, L2)
                : popcount64(mask) % 2 == 0;
        const MVMCKrylovStatus status = evaluate(workspace, mask, &result);
        if (!allowed) {
          CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                "configuration outside the domain is rejected");
          continue;
        }
        if (status != MVMC_KRYLOV_STATUS_OK) {
          fprintf(stderr, "domain=%d proj=%d mask=0x%llx status=%s\n",
                  (int)domains[d], with_projection,
                  (unsigned long long)mask,
                  mvmc_krylov_status_string(status));
        }
        CHECK(status == MVMC_KRYLOV_STATUS_OK, "allowed configuration");
        if (status == MVMC_KRYLOV_STATUS_OK) {
          double complex expected = slater_amplitude(slater, ORB2, mask);
          if (with_projection) expected *= exp(reference_log_projection(mask));
          CHECK(matches(&result, expected),
                "amplitude equals Pf(Slater[x,x]) times Re(Proj) factor");
        }
      }
      {
        uint64_t excess[1] = {UINT64_C(0x15)};
        MVMCKrylovScaledAmplitudeResult result;
        CHECK(mvmc_gc_krylov_scaled_amplitude(excess, 1, workspace,
                                              &result) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
              "excess configuration bit rejected");
        CHECK(mvmc_gc_krylov_scaled_amplitude(excess, 2, workspace,
                                              &result) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
              "word count rejected");
      }
      CHECK(memcmp(slater, slater_copy, sizeof(slater)) == 0 &&
                memcmp(fixture.projection, fixture_copy.projection,
                       sizeof(fixture.projection)) == 0 &&
                memcmp(fixture.jastrow_rows, fixture_copy.jastrow_rows,
                       sizeof(fixture.jastrow_rows)) == 0 &&
                memcmp(fixture.gutzwiller_idx, fixture_copy.gutzwiller_idx,
                       sizeof(fixture.gutzwiller_idx)) == 0 &&
                fixture.layout.qp_start == fixture_copy.layout.qp_start &&
                fixture.layout.qp_end == fixture_copy.layout.qp_end,
            "caller buffers are not modified");
      mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
    }
  }
}

static void test_permutation_parity(void) {
  /* Pf over a permuted occupied order times its parity is the same value. */
  double complex slater[ORB2 * ORB2];
  double complex sub[MAX_ORB * MAX_ORB];
  const int order[4] = {2, 0, 3, 1};
  int parity = 1, i, j, r, c;
  antiparallel_slater(kF2, slater);
  for (i = 0; i < 4; ++i) {
    for (j = i + 1; j < 4; ++j) {
      if (order[i] > order[j]) parity = -parity;
    }
  }
  for (r = 0; r < 4; ++r) {
    for (c = 0; c < 4; ++c) {
      sub[r * MAX_ORB + c] = -slater[order[r] * ORB2 + order[c]];
    }
  }
  /* (-1)^(N/2) Pf(-S[p(x)]) parity(p) = Pf(S[x]) for x ascending. */
  CHECK(cabs((double)parity * pfaffian(sub, 4, MAX_ORB) -
             slater_amplitude(slater, ORB2, UINT64_C(0xf))) <= 1e-14,
        "eleIdx permutation parity and the (-1)^(N/2) sign");
}

static void test_boundary_rule(void) {
  double complex slater[2 * ORB2 * ORB2];
  const double complex exact[2] = {1.0, -1.0};
  const double complex tiny[3] = {1.0, 0x1p-53, -1.0};
  double complex weights3[3];
  AmplitudeFixture fixture;
  MVMCGCKrylovAmplitudeWorkspace *workspace = NULL;
  MVMCKrylovScaledAmplitudeResult result;
  int index;
  antiparallel_slater(kF2, slater);
  for (index = 0; index < ORB2 * ORB2; ++index) {
    slater[ORB2 * ORB2 + index] = 2.0 * slater[index];
  }
  fixture_layout(&fixture, MVMC_KRYLOV_DOMAIN_GC_BALANCED, 2, exact, 2, 0,
                 world_comm(), 1);
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &fixture.layout, slater, exact, &workspace) ==
            MVMC_KRYLOV_STATUS_OK,
        "boundary workspace");
  if (workspace != NULL) {
    CHECK(evaluate(workspace, 0, &result) == MVMC_KRYLOV_STATUS_OK &&
              result.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
          "vacuum is the exact boundary zero");
    CHECK(evaluate(workspace, 0xf, &result) == MVMC_KRYLOV_STATUS_OK &&
              result.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
          "full filling is the exact boundary zero");
    CHECK(evaluate(workspace, 0x5, &result) == MVMC_KRYLOV_STATUS_OK &&
              matches(&result, -slater_amplitude(slater, ORB2, 0x5)),
          "N=2 keeps 1 - 2 = -1 times Pf");
    mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
    workspace = NULL;
  }
  /* 1 + 2^-53 - 1 is not exactly zero: no override. */
  memcpy(slater + ORB2 * ORB2, slater, sizeof(double complex) * ORB2 * ORB2);
  {
    double complex slater3[3 * ORB2 * ORB2];
    for (index = 0; index < 3; ++index) {
      memcpy(slater3 + index * ORB2 * ORB2, slater,
             sizeof(double complex) * ORB2 * ORB2);
      weights3[index] = tiny[index];
    }
    fixture_layout(&fixture, MVMC_KRYLOV_DOMAIN_GC_BALANCED, 3, tiny, 3, 0,
                   world_comm(), 1);
    CHECK(mvmc_gc_krylov_amplitude_workspace_create(
              &fixture.layout, slater3, weights3, &workspace) ==
              MVMC_KRYLOV_STATUS_OK,
          "near-boundary workspace");
    if (workspace != NULL) {
      CHECK(evaluate(workspace, 0, &result) == MVMC_KRYLOV_STATUS_OK &&
                result.value.state != MVMC_SCALED_COMPLEX_EXACT_ZERO,
            "near-zero translation sum is not overridden");
      mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
      workspace = NULL;
    }
  }
  /* A nonfinite component is reported, never hidden by the boundary zero. */
  slater[ORB2 * ORB2 + 0 * ORB2 + 2] = NAN;
  fixture_layout(&fixture, MVMC_KRYLOV_DOMAIN_GC_BALANCED, 2, exact, 2, 0,
                 world_comm(), 1);
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &fixture.layout, slater, exact, &workspace) ==
            MVMC_KRYLOV_STATUS_OK,
        "nonfinite Slater is accepted at creation");
  if (workspace != NULL) {
    CHECK(evaluate(workspace, 0xf, &result) == MVMC_KRYLOV_STATUS_NONFINITE,
          "nonfinite full-filling component is not hidden");
    CHECK(evaluate(workspace, 0, &result) == MVMC_KRYLOV_STATUS_OK &&
              result.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
          "vacuum does not read the nonfinite entry");
    mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
  }
}

static void test_node_and_vacuum(void) {
  /* Node fixture F=[[1,0],[1,1]]: q0(0x9)=0, q1(0x9)=-2 with hopping -1. */
  static const double complex node[L2][L2] = {{1.0, 0.0}, {1.0, 1.0}};
  static const MVMCKrylovHamiltonianTerm terms[4] = {
      {-1.0 + 0.0 * I, 0, 2, 1, 0}, {-1.0 + 0.0 * I, 2, 2, 1, 1},
      {-1.0 + 0.0 * I, 4, 2, 1, 2}, {-1.0 + 0.0 * I, 6, 2, 1, 3}};
  static const MVMCKrylovFermionOperator operators[8] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
      {MVMC_KRYLOV_FERMION_CREATE, 2}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 3},
      {MVMC_KRYLOV_FERMION_CREATE, 3}, {MVMC_KRYLOV_FERMION_ANNIHILATE, 2}};
  double complex slater[ORB2 * ORB2];
  const double complex translation[1] = {1.0};
  const double complex weight[1] = {0.5};
  AmplitudeFixture fixture;
  MVMCGCKrylovAmplitudeWorkspace *workspace = NULL;
  MVMCKrylovFockModel model;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *engine = NULL;
  MVMCKrylovBoundedResult values;
  MVMCKrylovScaledAmplitudeResult result;
  uint64_t words = UINT64_C(0x9);
  antiparallel_slater(node, slater);
  fixture_layout(&fixture, MVMC_KRYLOV_DOMAIN_GC_BALANCED, 1, translation, 1,
                 0, world_comm(), 1);
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &fixture.layout, slater, weight, &workspace) ==
            MVMC_KRYLOV_STATUS_OK,
        "node workspace");
  if (workspace == NULL) return;
  CHECK(evaluate(workspace, 0, &result) == MVMC_KRYLOV_STATUS_OK &&
            matches(&result, 0.5),
        "vacuum amplitude is the QP weight");
  CHECK(evaluate(workspace, 0x9, &result) == MVMC_KRYLOV_STATUS_OK &&
            result.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
        "node configuration q0 is an exact zero");
  memset(&model, 0, sizeof(model));
  model.site_count = L2;
  model.hermitian = 1;
  model.terms = terms;
  model.term_count = 4;
  model.operators = operators;
  model.operator_count = 8;
  model.domain = MVMC_KRYLOV_DOMAIN_GC_BALANCED;
  memset(&limits, 0, sizeof(limits));
  limits.amplitude_policy_hash =
      mvmc_gc_krylov_amplitude_generation_hash(workspace);
  limits.cache_bytes = 1u << 16;
  limits.max_row_transitions = 4;
  limits.max_workspace_bytes = (size_t)64 * 1024 * 1024;
  limits.max_node_expansions = 1000;
  limits.max_terminal_amplitude_calls = 1000;
  limits.max_total_row_transitions = 10000;
  limits.max_order = 2;
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_bounded_krylov_workspace_create(plan, &engine) ==
                MVMC_KRYLOV_STATUS_OK,
        "bounded engine for the node fixture");
  if (engine != NULL) {
    double complex q1 = 0.0;
    CHECK(mvmc_bounded_krylov_evaluate(engine, &words, 1,
                                       mvmc_gc_krylov_scaled_amplitude,
                                       workspace, &values) ==
              MVMC_KRYLOV_STATUS_OK &&
              values.value[0].state == MVMC_SCALED_COMPLEX_EXACT_ZERO &&
              mvmc_scaled_complex_export_common_scale(&values.value[1], 0.0,
                                                      &q1) ==
                  MVMC_SCALED_EXPORT_OK,
          "node q0/q1 through the bounded engine");
    /* Weight 0.5 scales the fixture's q1 = -2 to -1. */
    CHECK(cabs(q1 - (-1.0)) <= 1e-13, "node q1 = 0.5 * (-2)");
  }
  mvmc_bounded_krylov_workspace_destroy(engine);
  mvmc_bounded_krylov_plan_destroy(plan);
  mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
}

static void test_invalid_layouts(void) {
  double complex slater[ORB2 * ORB2];
  const double complex translation[1] = {1.0};
  const double complex nan_translation[1] = {NAN};
  const double complex weight[1] = {1.0};
  const double complex nan_weight[1] = {NAN};
  AmplitudeFixture fixture;
  MVMCGCKrylovAmplitudeWorkspace *workspace = NULL;
  antiparallel_slater(kF2, slater);
  fixture_layout(&fixture, MVMC_KRYLOV_DOMAIN_FIXED, 1, translation, 1, 0,
                 world_comm(), 1);
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &fixture.layout, slater, weight, &workspace) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            workspace == NULL,
        "fixed domain rejected");
  fixture_layout(&fixture, MVMC_KRYLOV_DOMAIN_GC_EVEN, 1, translation, 1, 0,
                 world_comm(), 1);
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &fixture.layout, slater, nan_weight, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "nonfinite QP weight rejected");
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &fixture.layout, NULL, weight, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "missing Slater rejected");
  fixture.layout.translation_weights = nan_translation;
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &fixture.layout, slater, weight, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "nonfinite translation weight rejected");
  fixture_layout(&fixture, MVMC_KRYLOV_DOMAIN_GC_EVEN, 1, translation, 1, 0,
                 world_comm(), 1);
  fixture.layout.projection.site_count = 3;
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &fixture.layout, slater, weight, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "projection site mismatch rejected");
  fixture_layout(&fixture, MVMC_KRYLOV_DOMAIN_GC_EVEN, 0, translation, 1, 0,
                 world_comm(), 1);
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &fixture.layout, slater, weight, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "qp_total=0 rejected");
  fixture_layout(&fixture, MVMC_KRYLOV_DOMAIN_GC_EVEN, 1, translation, 1, 0,
                 world_comm(), 1);
  fixture.layout.site_count = -1;
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &fixture.layout, slater, weight, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "negative site count rejected");
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(NULL, slater, weight,
                                                  &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "NULL layout rejected");
}

/* QP partition of the chain communicator against a COMM_SELF reference. */
static void test_partition(void) {
  double complex slater[3 * ORB2 * ORB2];
  const double complex translation[3] = {0.5, 0.25 + 0.25 * I, -0.125};
  AmplitudeFixture partitioned;
  AmplitudeFixture serial;
  MVMCGCKrylovAmplitudeWorkspace *split = NULL;
  MVMCGCKrylovAmplitudeWorkspace *whole = NULL;
  uint64_t mask;
  int qp, index;
  for (qp = 0; qp < 3; ++qp) {
    antiparallel_slater(kF2, slater + qp * ORB2 * ORB2);
    for (index = 0; index < ORB2 * ORB2; ++index) {
      slater[qp * ORB2 * ORB2 + index] *= 1.0 + 0.5 * qp;
    }
  }
  fixture_layout(&partitioned, MVMC_KRYLOV_DOMAIN_GC_BALANCED, 3, translation,
                 3, 1, world_comm(), 1);
  fixture_layout(&serial, MVMC_KRYLOV_DOMAIN_GC_BALANCED, 3, translation, 3,
                 1, self_comm(), 0);
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &partitioned.layout, slater, translation, &split) ==
                MVMC_KRYLOV_STATUS_OK &&
            mvmc_gc_krylov_amplitude_workspace_create(
                &serial.layout, slater, translation, &whole) ==
                MVMC_KRYLOV_STATUS_OK,
        "partitioned and serial workspaces");
  if (split != NULL && whole != NULL) {
    CHECK(mvmc_gc_krylov_amplitude_generation_hash(split) ==
              mvmc_gc_krylov_amplitude_generation_hash(whole),
          "generation hash independent of the QP partition");
    for (mask = 0; mask < 16; ++mask) {
      MVMCKrylovScaledAmplitudeResult left, right;
      if (!balanced_mask(mask, L2)) continue;
      CHECK(evaluate(split, mask, &left) == MVMC_KRYLOV_STATUS_OK &&
                evaluate(whole, mask, &right) == MVMC_KRYLOV_STATUS_OK,
            "partitioned evaluation");
      CHECK(left.value.state == right.value.state &&
                (left.value.state != MVMC_SCALED_COMPLEX_FINITE_NONZERO ||
                 (cabs(left.value.phase - right.value.phase) <= 3e-12 &&
                  fabs(left.value.log_abs - right.value.log_abs) <= 3e-11)),
            "QP partition does not change the amplitude");
    }
  }
  mvmc_gc_krylov_amplitude_workspace_destroy(split);
  mvmc_gc_krylov_amplitude_workspace_destroy(whole);
}

static void test_rank_faults(void) {
  double complex slater[2 * ORB2 * ORB2];
  const double complex translation[2] = {0.5, 0.5};
  AmplitudeFixture fixture;
  MVMCGCKrylovAmplitudeWorkspace *workspace = NULL;
  MVMCKrylovScaledAmplitudeResult result;
  const int target = world_size() > 1 ? 1 : 0;
  antiparallel_slater(kF2, slater);
  antiparallel_slater(kF2, slater + ORB2 * ORB2);
  if (world_size() > 1) {
    /* Replicated-binding mismatch on rank 1. */
    if (world_rank() == 1) slater[ORB2 * ORB2 + 2] *= 2.0;
    fixture_layout(&fixture, MVMC_KRYLOV_DOMAIN_GC_BALANCED, 2, translation,
                   2, 0, world_comm(), 1);
    CHECK(mvmc_gc_krylov_amplitude_workspace_create(
              &fixture.layout, slater, translation, &workspace) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              workspace == NULL,
          "rank-1 binding mismatch rejected on every rank");
    antiparallel_slater(kF2, slater + ORB2 * ORB2);
  }
#if defined(MVMC_ENABLE_FAULT_INJECTION)
  fixture_layout(&fixture, MVMC_KRYLOV_DOMAIN_GC_BALANCED, 2, translation, 2,
                 0, world_comm(), 1);
  mvmc_gc_krylov_amplitude_testing_set_fault(
      MVMC_GC_KRYLOV_AMPLITUDE_FAULT_ALLOCATION, target);
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &fixture.layout, slater, translation, &workspace) ==
                MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE &&
            workspace == NULL,
        "allocation fault reported on every rank");
  mvmc_gc_krylov_amplitude_testing_set_fault(
      MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONE, -1);
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &fixture.layout, slater, translation, &workspace) ==
            MVMC_KRYLOV_STATUS_OK,
        "workspace after clearing the fault");
  if (workspace != NULL) {
    mvmc_gc_krylov_amplitude_testing_set_fault(
        MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONFINITE, target);
    CHECK(evaluate(workspace, 0x5, &result) == MVMC_KRYLOV_STATUS_NONFINITE,
          "nonfinite fault reported on every rank");
    mvmc_gc_krylov_amplitude_testing_set_fault(
        MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONE, -1);
    CHECK(evaluate(workspace, 0x5, &result) == MVMC_KRYLOV_STATUS_OK,
          "evaluation after clearing the fault");
    mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
  }
#else
  (void)target;
  (void)result;
  (void)workspace;
  (void)fixture;
#endif
}

static void test_gc_bridge(void) {
  /* Production-shaped rows: t=1 hopping, mu=1/4, U=3/2 on both sites. */
  int transfer_rows[6][4] = {{0, 0, 1, 0}, {1, 0, 0, 0}, {0, 1, 1, 1},
                             {1, 1, 0, 1}, {0, 0, 0, 0}, {1, 1, 1, 1}};
  int *transfer_index[6];
  const double complex transfer_value[6] = {1.0, 1.0, 1.0, 1.0, 0.25, 0.25};
  const int intra_index[2] = {0, 1};
  const double intra_value[2] = {1.5, 1.5};
  const double complex translation[1] = {1.0};
  double complex slater[ORB2 * ORB2];
  double complex weight[1] = {1.0};
  MVMCPowerLanczosGCView view;
  MVMCPowerLanczosGCBridge *bridge = NULL;
  int index;
  int chain_rank = 0;
  int chain_size = 1;
#ifdef _mpi_use
  MPI_Comm_rank(MPI_COMM_WORLD, &chain_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &chain_size);
#endif
  antiparallel_slater(kF2, slater);
  for (index = 0; index < 6; ++index) transfer_index[index] = transfer_rows[index];
  memset(&view, 0, sizeof(view));
  view.data.site_count = L2;
  view.data.arithmetic = MVMC_POWER_LANCZOS_CLASSIC_COMPLEX;
  view.data.transfer_count = 6;
  view.data.transfer_indices = world_rank() == 0 ? transfer_index : NULL;
  view.data.transfer_parameters = world_rank() == 0 ? transfer_value : NULL;
  view.data.coulomb_intra_count = 2;
  view.data.coulomb_intra_indices = world_rank() == 0 ? intra_index : NULL;
  view.data.coulomb_intra_parameters =
      world_rank() == 0 ? intra_value : NULL;
  view.data.qp_total = 1;
  view.data.qp_start = chain_rank * 1 / chain_size;
  view.data.qp_end = (chain_rank + 1) * 1 / chain_size;
  view.data.qp_weights = weight;
  view.data.slater_complex = slater;
  view.domain = MVMC_KRYLOV_DOMAIN_GC_BALANCED;
  view.translation_count = 1;
  view.translation_weights = translation;
  CHECK(mvmc_power_lanczos_gc_bridge_create(&view, world_comm(), world_comm(),
                                            &bridge) ==
                MVMC_KRYLOV_STATUS_OK &&
            bridge != NULL,
        "GC bridge from a production-shaped view");
  if (bridge != NULL) {
    const MVMCKrylovFockModel *model =
        mvmc_power_lanczos_gc_bridge_model(bridge);
    MVMCKrylovScaledAmplitudeResult result;
    uint64_t words = UINT64_C(0x5);
    CHECK(model != NULL && model->domain == MVMC_KRYLOV_DOMAIN_GC_BALANCED &&
              model->up_electron_count == 0 && model->term_count == 8,
          "bridge model has the GC domain and the converted rows");
    CHECK(mvmc_power_lanczos_gc_bridge_amplitude(bridge) ==
                  mvmc_gc_krylov_scaled_amplitude &&
              mvmc_power_lanczos_gc_bridge_amplitude_context(bridge) != NULL &&
              mvmc_power_lanczos_gc_bridge_generation_hash(bridge) != 0,
          "bridge amplitude callback and generation");
    CHECK(mvmc_power_lanczos_gc_bridge_amplitude(bridge)(
              &words, 1, mvmc_power_lanczos_gc_bridge_amplitude_context(bridge),
              &result) == MVMC_KRYLOV_STATUS_OK &&
              matches(&result, slater_amplitude(slater, ORB2, words)),
          "bridge amplitude");
    mvmc_power_lanczos_gc_bridge_destroy(bridge);
    bridge = NULL;
  }
  view.data.up_electron_count = 1;
  CHECK(mvmc_power_lanczos_gc_bridge_create(&view, world_comm(), world_comm(),
                                            &bridge) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "GC view rejects fixed counts");
  view.data.up_electron_count = 0;
  view.data.unsupported_amplitude_features = 4u;
  CHECK(mvmc_power_lanczos_gc_bridge_create(&view, world_comm(), world_comm(),
                                            &bridge) ==
            MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL,
        "GC view rejects the canonical general-orbital bit");
  view.data.unsupported_amplitude_features = 0u;
  view.data.exchange_count = 1;
  CHECK(mvmc_power_lanczos_gc_bridge_create(&view, world_comm(), world_comm(),
                                            &bridge) ==
            MVMC_KRYLOV_STATUS_UNSUPPORTED_MODEL,
        "GC view rejects Exchange");
  view.data.exchange_count = 0;
  view.data.arithmetic = MVMC_POWER_LANCZOS_CLASSIC_REAL;
  CHECK(mvmc_power_lanczos_gc_bridge_create(&view, world_comm(), world_comm(),
                                            &bridge) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "GC view requires complex arithmetic");
  view.data.arithmetic = MVMC_POWER_LANCZOS_CLASSIC_COMPLEX;
  view.domain = MVMC_KRYLOV_DOMAIN_FIXED;
  CHECK(mvmc_power_lanczos_gc_bridge_create(&view, world_comm(), world_comm(),
                                            &bridge) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            bridge == NULL,
        "GC view requires a GC domain");
}

int main(int argc, char **argv) {
  int total = 0;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#else
  (void)argc;
  (void)argv;
#endif
  test_single_qp_amplitudes();
  test_permutation_parity();
  test_boundary_rule();
  test_node_and_vacuum();
  test_invalid_layouts();
  test_partition();
  test_rank_faults();
  test_gc_bridge();
#ifdef _mpi_use
  MPI_Allreduce(&failure_count, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  total = failure_count;
#endif
  if (world_rank() == 0) {
    printf(total == 0 ? "gc_krylov_amplitude_unit: PASS\n"
                      : "gc_krylov_amplitude_unit: FAIL\n");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return total == 0 ? 0 : 1;
}
