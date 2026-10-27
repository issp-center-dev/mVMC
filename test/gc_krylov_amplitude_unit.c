/*
 * Grand-canonical absolute amplitude and bridge contract: the production sign
 * (-1)^(N/2) Pf(-Slater[x,x]) equals Pf(Slater[x,x]) over ascending orbitals
 * (independent recursive Pfaffian), Re(Proj) correlation factor, Phase 3b
 * boundary zero, node configurations, eleIdx permutation parity, caller
 * buffer invariance, input validation, MPI QP partitions with empty ranks,
 * rank-1 faults and the production-shaped GC bridge.  Task 5: the spin
 * support descriptor (S0 canonical hash, S>0 frozen R1/E1/R1_G2/S2_G2 base
 * states, H Psi/Phi never masked, singlet/near-triplet, NaN under a mask,
 * private snapshot, ordinary/refined mask mutations on S2_G2, metadata
 * audit and hash coverage, support-helper failure agreement).
 */
#include "bounded_krylov_engine.h"
#include "gc_krylov_amplitude.h"
#include "gc_lanczos_task5_case.h"
#include "krylov_final_state.h"
#include "power_lanczos_gc_bridge.h"

#include <complex.h>
#include <limits.h>
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
  fixture->layout.spin_support = mvmc_gc_spin_support_inactive();
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

static void test_spin_parity_rule(void) {
  double complex slater[2 * ORB2 * ORB2];
  const double complex translation[1] = {1};
  const double complex weights[2] = {0.5, 0.5};
  AmplitudeFixture fixture;
  MVMCGCKrylovAmplitudeWorkspace *masked = NULL, *raw = NULL;
  MVMCKrylovScaledAmplitudeResult result;
  int q;
  for (q = 0; q < 2; ++q) {
    const double beta = q == 0 ? 0.7 : acos(-1.0)-0.7;
    const double c = cos(beta/2), s = sin(beta/2);
    double complex *a = slater + q*ORB2*ORB2;
    memset(a, 0, ORB2*ORB2*sizeof(*a));
    a[1] = -2*c*s; a[4] = -a[1];
    a[11] = 2*c*s; a[14] = -a[11];
    a[3] = c*c-s*s; a[12] = -a[3];
    a[6] = -a[3]; a[9] = -a[6];
  }
  fixture_layout(&fixture, MVMC_KRYLOV_DOMAIN_GC_BALANCED, 2, translation,
                 1, 0, world_comm(), 1);
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &fixture.layout, slater, weights, &raw) == MVMC_KRYLOV_STATUS_OK,
        "raw spin workspace");
  CHECK(mvmc_gc_spin_support_classify(1, 2, 0, 1, 1,
                                      &fixture.layout.spin_support) ==
                MVMC_KRYLOV_STATUS_OK &&
            fixture.layout.spin_support.pair_symmetry == MVMC_GC_PAIR_TRIPLET,
        "S=0 triplet descriptor (singlet not read)");
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(
            &fixture.layout, slater, weights, &masked) == MVMC_KRYLOV_STATUS_OK,
        "triplet spin workspace");
  if (raw && masked) {
    CHECK(mvmc_gc_krylov_amplitude_generation_hash(raw) !=
              mvmc_gc_krylov_amplitude_generation_hash(masked),
          "support proof is part of generation binding");
    CHECK(evaluate(masked, 0x9, &result) == MVMC_KRYLOV_STATUS_OK &&
              result.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
          "triplet N=2 is exact zero");
    CHECK(evaluate(masked, 0, &result) == MVMC_KRYLOV_STATUS_OK &&
              matches(&result, 1), "triplet vacuum survives");
    CHECK(evaluate(masked, 0xf, &result) == MVMC_KRYLOV_STATUS_OK &&
              matches(&result, -1), "triplet full survives");
    /* Mutating the caller's proof cannot change an existing generation. */
    fixture.layout.spin_support = mvmc_gc_spin_support_inactive();
    CHECK(evaluate(masked, 0x9, &result) == MVMC_KRYLOV_STATUS_OK &&
              result.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
          "immutable support snapshot");
    mvmc_gc_krylov_amplitude_testing_set_fault(
        MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONFINITE,
        world_size() > 1 ? 1 : 0);
    CHECK(evaluate(masked, 0x9, &result) == MVMC_KRYLOV_STATUS_NONFINITE,
          "rank fault not hidden by triplet mask");
    mvmc_gc_krylov_amplitude_testing_set_fault(
        MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONE, -1);
  }
  mvmc_gc_krylov_amplitude_workspace_destroy(raw);
  mvmc_gc_krylov_amplitude_workspace_destroy(masked);
  masked = NULL;
  if (world_size() > 1) {
    fixture.layout.spin_support.active = 1;
    fixture.layout.spin_support.gauss_count = 2;
    fixture.layout.spin_support.pair_symmetry =
        world_rank() == 0 ? MVMC_GC_PAIR_TRIPLET : MVMC_GC_PAIR_GENERAL;
    CHECK(mvmc_gc_krylov_amplitude_workspace_create(
              &fixture.layout, slater, weights, &masked) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && masked == NULL,
          "rank-mismatched support proof rejected");
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
  view.spin_support = mvmc_gc_spin_support_inactive();
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

/* ---- Task 5: spin support descriptor ------------------------------------ */

static void fault_everywhere(MVMCGCKrylovAmplitudeFault fault) {
  mvmc_gc_krylov_amplitude_testing_set_fault(
      fault, fault == MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONE ? -1 : world_rank());
}

static MVMCGCKrylovAmplitudeWorkspace *task5_workspace(
    const Task5Built *built, MVMCKrylovStatus *status) {
  MVMCGCKrylovAmplitudeLayout layout;
  MVMCGCKrylovAmplitudeWorkspace *workspace = NULL;
  task5_layout(built, world_comm(), world_rank(), world_size(), &layout);
  *status = mvmc_gc_krylov_amplitude_workspace_create(
      &layout, built->slater, built->weight, &workspace);
  return workspace;
}

static int task5_engine(const Task5Built *built,
                        MVMCGCKrylovAmplitudeWorkspace *workspace,
                        MVMCKrylovBoundedPlan **plan,
                        MVMCKrylovBoundedWorkspace **engine) {
  MVMCKrylovBoundedLimits limits;
  memset(&limits, 0, sizeof(limits));
  limits.amplitude_policy_hash =
      mvmc_gc_krylov_amplitude_generation_hash(workspace);
  limits.cache_bytes = (size_t)1 << 20;
  limits.max_row_transitions = built->model.term_count;
  limits.max_workspace_bytes = (size_t)64 * 1024 * 1024;
  limits.max_node_expansions = 100000;
  limits.max_terminal_amplitude_calls = 100000;
  limits.max_total_row_transitions = 1000000;
  limits.max_order = 2;
  *plan = NULL;
  *engine = NULL;
  return mvmc_bounded_krylov_plan_create(&built->model, &limits, plan) ==
             MVMC_KRYLOV_STATUS_OK &&
         mvmc_bounded_krylov_workspace_create(*plan, engine) ==
             MVMC_KRYLOV_STATUS_OK;
}

/* Psi from the callback, Psi/H Psi/H^2 Psi from the bounded engine and
 * Phi/H Phi from the raw policy, each against the frozen 64-eps rule. */
static void task5_check_values(const char *label, const Task5Built *built,
                               const Task5Scaled *table, size_t count,
                               const double complex coefficient[2]) {
  MVMCKrylovStatus status;
  MVMCGCKrylovAmplitudeWorkspace *workspace = task5_workspace(built, &status);
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *engine = NULL;
  MVMCKrylovFinalStatePolicy policy;
  size_t k;
  int compared = 0;
  CHECK(status == MVMC_KRYLOV_STATUS_OK && workspace != NULL, label);
  if (workspace == NULL) return;
  CHECK(task5_engine(built, workspace, &plan, &engine), "Task 5 engine");
  CHECK(mvmc_krylov_final_state_policy_create(
            1, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC, 31, coefficient, 2,
            &policy) == MVMC_KRYLOV_STATUS_OK,
        "raw synthetic policy");
  for (k = 0; engine != NULL && k < count; ++k) {
    uint64_t words = table[k].configuration;
    MVMCKrylovScaledAmplitudeResult base;
    MVMCKrylovBoundedResult values;
    MVMCKrylovFinalStateEvaluation final_state;
    const Task5Scaled *phi, *hphi;
    int v;
    if (table[k].value != 0) continue;
    CHECK(evaluate(workspace, words, &base) == MVMC_KRYLOV_STATUS_OK &&
              task5_scaled_matches(&base.value, table + k, 64.0),
          "base Psi state and value (64 eps)");
    status = mvmc_bounded_krylov_evaluate(engine, &words, 1,
                                          mvmc_gc_krylov_scaled_amplitude,
                                          workspace, &values);
    CHECK(status == MVMC_KRYLOV_STATUS_OK, "bounded Krylov values");
    if (status != MVMC_KRYLOV_STATUS_OK) continue;
    for (v = 0; v < 3; ++v) {
      const Task5Scaled *ref = task5_reference(table, count, (unsigned)words, v);
      if (ref == NULL || !task5_scaled_matches(values.value + v, ref, 64.0)) {
        fprintf(stderr, "%s: mask %u q%d state %d log %.17g\n", label,
                (unsigned)words, v, (int)values.value[v].state,
                values.value[v].log_abs);
      }
      CHECK(ref != NULL && task5_scaled_matches(values.value + v, ref, 64.0),
            "Psi/H Psi/H^2 Psi from the masked base only");
    }
    phi = task5_reference(table, count, (unsigned)words, 3);
    hphi = task5_reference(table, count, (unsigned)words, 4);
    status = mvmc_krylov_final_state_evaluate(&policy, &values, &final_state);
    CHECK(status == MVMC_KRYLOV_STATUS_OK && phi != NULL && hphi != NULL &&
              task5_scaled_matches(&final_state.amplitude, phi, 64.0) &&
              task5_scaled_matches(&final_state.hamiltonian_amplitude, hphi,
                                   64.0),
          "Phi and H Phi are not masked outside the base");
    ++compared;
  }
  CHECK(compared > 0, "reference configurations compared");
  mvmc_bounded_krylov_workspace_destroy(engine);
  mvmc_bounded_krylov_plan_destroy(plan);
  mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
}

static void test_task5_frozen_values(void) {
  Task5Built built;
  double complex c[2];
  CHECK(task5_build(&task5_r1, &built), "R1 build");
  c[0] = built.coefficient[0];
  c[1] = built.coefficient[1];
  task5_check_values("R1 S1 G8 triplet", &built, task5_r1.scaled,
                     task5_r1.scaled_count, c);
  c[1] = task5_complex(task5_r1c_c1);
  task5_check_values("R1C complex c1", &built, task5_r1c_scaled,
                     sizeof(task5_r1c_scaled) / sizeof(task5_r1c_scaled[0]), c);
  task5_free(&built);
  CHECK(task5_build(&task5_e1, &built), "E1 build");
  task5_check_values("E1 S1 G8 general exp(D/8)", &built, task5_e1.scaled,
                     task5_e1.scaled_count, built.coefficient);
  task5_free(&built);
  CHECK(task5_build(&task5_s2_g2, &built), "S2_G2 build");
  task5_check_values("S2_G2 S2 G2 general", &built, task5_s2_g2.scaled,
                     task5_s2_g2.scaled_count, built.coefficient);
  task5_free(&built);
  /* R1_G2 (2 QPs, empty QP rank at np3): states and the 6/9 symmetry. */
  CHECK(task5_build(&task5_r1_g2, &built), "R1_G2 build");
  {
    MVMCKrylovStatus status;
    MVMCGCKrylovAmplitudeWorkspace *workspace = task5_workspace(&built, &status);
    MVMCKrylovScaledAmplitudeResult six, nine, five, vacuum;
    CHECK(status == MVMC_KRYLOV_STATUS_OK && workspace != NULL, "R1_G2 workspace");
    if (workspace != NULL) {
      CHECK(evaluate(workspace, 6, &six) == MVMC_KRYLOV_STATUS_OK &&
                evaluate(workspace, 9, &nine) == MVMC_KRYLOV_STATUS_OK &&
                evaluate(workspace, 5, &five) == MVMC_KRYLOV_STATUS_OK &&
                evaluate(workspace, 0, &vacuum) == MVMC_KRYLOV_STATUS_OK,
            "R1_G2 evaluations");
      CHECK(six.value.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
                nine.value.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
                cabs(six.value.phase + nine.value.phase) <= 64 * 0x1p-52 &&
                fabs(six.value.log_abs - nine.value.log_abs) <=
                    64 * 0x1p-52 * (1.0 + fabs(six.value.log_abs)) &&
                five.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO &&
                vacuum.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
            "R1_G2: N2 singles survive (no S=0 rule), doublon/vacuum masked");
      mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
    }
  }
  task5_free(&built);
}

static void test_task5_mask_mutations(void) {
  Task5Built built;
  MVMCKrylovStatus status;
  MVMCGCKrylovAmplitudeWorkspace *workspace;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *engine = NULL;
  MVMCKrylovFinalStatePolicy policy;
  MVMCKrylovScaledAmplitudeResult result;
  const uint64_t current = 60, proposal = 53;
  int mutation;
  CHECK(task5_build(&task5_s2_g2, &built), "S2_G2 build");
  workspace = task5_workspace(&built, &status);
  CHECK(workspace != NULL, "S2_G2 workspace");
  if (workspace == NULL) {
    task5_free(&built);
    return;
  }
  CHECK(task5_engine(&built, workspace, &plan, &engine), "S2_G2 engine");
  CHECK(mvmc_krylov_final_state_policy_create(
            1, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC, 37, built.coefficient,
            2, &policy) == MVMC_KRYLOV_STATUS_OK,
        "S2_G2 raw policy");
  /* Same N4: 60 has four singles, 53 two (S=2 needs four). */
  CHECK(evaluate(workspace, proposal, &result) == MVMC_KRYLOV_STATUS_OK &&
            result.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
        "S2_G2 Psi(53) is an exact structural zero");
  fault_everywhere(MVMC_GC_KRYLOV_AMPLITUDE_FAULT_ORDINARY_UNMASKED);
  CHECK(evaluate(workspace, proposal, &result) == MVMC_KRYLOV_STATUS_OK &&
            result.value.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO,
        "ordinary mask removal is visible in the base state");
  fault_everywhere(MVMC_GC_KRYLOV_AMPLITUDE_FAULT_REFINED_UNMASKED);
  CHECK(evaluate(workspace, proposal, &result) == MVMC_KRYLOV_STATUS_OK &&
            result.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
        "refined-only mutation leaves the ordinary base masked");
  fault_everywhere(MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONE);
  for (mutation = 0; engine != NULL && mutation < 2; ++mutation) {
    MVMCKrylovFinalStateSnapshot snapshot;
    MVMCKrylovFinalStateStepResult step;
    uint64_t words = current;
    CHECK(mvmc_krylov_final_state_sampler_initialize(
              engine, &policy, &words, 1, mvmc_gc_krylov_scaled_amplitude,
              workspace, &snapshot) == MVMC_KRYLOV_STATUS_OK,
          "S2_G2 current 60");
    if (mutation)
      fault_everywhere(MVMC_GC_KRYLOV_AMPLITUDE_FAULT_ORDINARY_UNMASKED);
    status = mvmc_krylov_final_state_sampler_step(
        engine, &policy, &words, 1, &snapshot, &proposal, 1, 0.0, 0.5,
        mvmc_gc_krylov_scaled_amplitude, workspace, &step);
    fault_everywhere(MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONE);
    if (!mutation) {
      CHECK(status == MVMC_KRYLOV_STATUS_OK && !step.accepted &&
                step.acceptance.exact_zero_proposal && words == current,
            "ordinary step: masked 53 is an exact-zero proposal");
    } else {
      CHECK(status == MVMC_KRYLOV_STATUS_OK && step.accepted &&
                words == proposal,
            "ordinary mask mutation turns the rejection into acceptance");
    }
  }
#if defined(MVMC_HAVE_MPFR_REFINEMENT)
  {
    MVMCKrylovRefinementResult refined;
    const MVMCGCKrylovAmplitudeFault faults[3] = {
        MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONE,
        MVMC_GC_KRYLOV_AMPLITUDE_FAULT_ORDINARY_UNMASKED,
        MVMC_GC_KRYLOV_AMPLITUDE_FAULT_REFINED_UNMASKED};
    int f;
    for (f = 0; f < 3; ++f) {
      fault_everywhere(faults[f]);
      memset(&refined, 0, sizeof(refined));
      status = mvmc_gc_krylov_refine_final(workspace, &built.model, &policy,
                                           &current, &proposal, 1, 0.5,
                                           &refined);
      fault_everywhere(MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONE);
      if (f < 2) {
        CHECK(status == MVMC_KRYLOV_STATUS_OK && refined.valid &&
                  !refined.accepted && refined.exact_zero_proposal,
              "refined S2_G2 60->53 u=1/2 rejects the masked zero");
      } else {
        CHECK(status == MVMC_KRYLOV_STATUS_OK && refined.valid &&
                  refined.accepted,
              "refined mask mutation turns the rejection into acceptance");
      }
    }
  }
#endif
  mvmc_bounded_krylov_workspace_destroy(engine);
  mvmc_bounded_krylov_plan_destroy(plan);
  mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
  task5_free(&built);
}

static void test_task5_private_snapshot(void) {
  Task5Built built;
  MVMCGCKrylovAmplitudeLayout layout;
  MVMCGCKrylovAmplitudeWorkspace *workspace = NULL;
  MVMCKrylovScaledAmplitudeResult before60, before53, after60, after53;
  uint64_t hash;
  size_t k;
  CHECK(task5_build(&task5_s2_g2, &built), "S2_G2 build");
  task5_layout(&built, world_comm(), world_rank(), world_size(), &layout);
  CHECK(mvmc_gc_krylov_amplitude_workspace_create(&layout, built.slater,
                                                  built.weight, &workspace) ==
            MVMC_KRYLOV_STATUS_OK,
        "S2_G2 workspace");
  if (workspace != NULL) {
    CHECK(evaluate(workspace, 60, &before60) == MVMC_KRYLOV_STATUS_OK &&
              evaluate(workspace, 53, &before53) == MVMC_KRYLOV_STATUS_OK,
          "evaluation before caller changes");
    hash = mvmc_gc_krylov_amplitude_generation_hash(workspace);
    /* Caller mutates every input it owns after creation. */
    layout.spin_support = mvmc_gc_spin_support_inactive();
    built.support.total_spin = 0;
    for (k = 0; k < (size_t)task5_s2_g2.qp_total * built.orbital_count *
                        built.orbital_count;
         ++k) {
      built.slater[k] *= 2.0;
    }
    built.weight[0] = -built.weight[0];
    built.translation[0] = 0.0;
    CHECK(evaluate(workspace, 60, &after60) == MVMC_KRYLOV_STATUS_OK &&
              evaluate(workspace, 53, &after53) == MVMC_KRYLOV_STATUS_OK &&
              memcmp(&after60.value, &before60.value, sizeof(after60.value)) ==
                  0 &&
              after53.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO &&
              before53.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO &&
              mvmc_gc_krylov_amplitude_generation_hash(workspace) == hash,
          "private descriptor/Slater/weight snapshot: values and hash unchanged");
    mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
  }
  task5_free(&built);
}

/* S0 canonical form, near/pure triplet, and an S>0 singlet. */
static void test_task5_classification(void) {
  double complex slater[2 * ORB2 * ORB2];
  const double complex translation[1] = {1};
  const double complex weights[2] = {0.5, 0.5};
  AmplitudeFixture fixture;
  MVMCGCKrylovAmplitudeWorkspace *a = NULL, *b = NULL;
  MVMCKrylovScaledAmplitudeResult result;
  double complex F[L2 * L2] = {0.0, 1.0, -1.0, 0.0};
  int idx_rows[L2][L2] = {{0, 1}, {2, 3}};
  int sgn_rows[L2][L2] = {{1, 1}, {1, 1}};
  int *index[L2] = {idx_rows[0], idx_rows[1]};
  int *sign[L2] = {sgn_rows[0], sgn_rows[1]};
  int triplet, singlet, q;
  for (q = 0; q < 2; ++q) {
    const double beta = q == 0 ? 0.7 : acos(-1.0) - 0.7;
    const double c = cos(beta / 2), s = sin(beta / 2);
    double complex *m = slater + q * ORB2 * ORB2;
    memset(m, 0, ORB2 * ORB2 * sizeof(*m));
    m[1] = -2 * c * s; m[4] = -m[1];
    m[11] = 2 * c * s; m[14] = -m[11];
    m[3] = c * c - s * s; m[12] = -m[3];
    m[6] = -m[3]; m[9] = -m[6];
  }
  fixture_layout(&fixture, MVMC_KRYLOV_DOMAIN_GC_BALANCED, 2, translation, 1,
                 0, world_comm(), 1);
  /* Real-F classification of an exact triplet (singlet 0) and of the zero
   * F (both 1) versus the production global path (singlet never read). */
  triplet = GCSpinPureTriplet(L2, 4, index, sign, F);
  singlet = GCSpinPureSinglet(L2, 4, index, sign, F);
  CHECK(triplet == 1 && singlet == 0, "exact triplet F classified");
  CHECK(mvmc_gc_spin_support_classify(1, 2, 0, triplet, singlet,
                                      &fixture.layout.spin_support) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_gc_krylov_amplitude_workspace_create(
                &fixture.layout, slater, weights, &a) == MVMC_KRYLOV_STATUS_OK,
        "real-F S0 descriptor");
  {
    double complex zero[L2 * L2] = {0.0, 0.0, 0.0, 0.0};
    const int zt = GCSpinPureTriplet(L2, 4, index, sign, zero);
    const int zs = GCSpinPureSinglet(L2, 4, index, sign, zero);
    CHECK(zt == 1 && zs == 1, "zero F is both triplet and singlet");
    CHECK(mvmc_gc_spin_support_classify(1, 2, 0, zt, zs,
                                        &fixture.layout.spin_support) ==
                  MVMC_KRYLOV_STATUS_OK &&
              fixture.layout.spin_support.pair_symmetry ==
                  MVMC_GC_PAIR_TRIPLET &&
              mvmc_gc_krylov_amplitude_workspace_create(
                  &fixture.layout, slater, weights, &b) ==
                  MVMC_KRYLOV_STATUS_OK,
          "S0 ZERO folds to TRIPLET");
  }
  if (a != NULL && b != NULL) {
    CHECK(mvmc_gc_krylov_amplitude_generation_hash(a) ==
              mvmc_gc_krylov_amplitude_generation_hash(b),
          "real-F and production-global S0 descriptors hash identically");
    CHECK(evaluate(a, 0x9, &result) == MVMC_KRYLOV_STATUS_OK &&
              result.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
          "S0 triplet parity zero (N%4==2)");
  }
  mvmc_gc_krylov_amplitude_workspace_destroy(a);
  mvmc_gc_krylov_amplitude_workspace_destroy(b);
  a = b = NULL;
  /* A near-triplet F is GENERAL: no parity mask from a tolerance. */
  F[1] = 1.0 + 0x1p-52;
  triplet = GCSpinPureTriplet(L2, 4, index, sign, F);
  CHECK(triplet == 0, "near-triplet F is not a triplet");
  CHECK(mvmc_gc_spin_support_classify(1, 2, 0, triplet, 0,
                                      &fixture.layout.spin_support) ==
                MVMC_KRYLOV_STATUS_OK &&
            fixture.layout.spin_support.pair_symmetry ==
                MVMC_GC_PAIR_GENERAL &&
            mvmc_gc_krylov_amplitude_workspace_create(
                &fixture.layout, slater, weights, &a) == MVMC_KRYLOV_STATUS_OK,
        "near-triplet descriptor");
  if (a != NULL) {
    CHECK(evaluate(a, 0x9, &result) == MVMC_KRYLOV_STATUS_OK &&
              result.value.state != MVMC_SCALED_COMPLEX_EXACT_ZERO,
          "GENERAL S0 keeps the N%4==2 configuration");
    mvmc_gc_krylov_amplitude_workspace_destroy(a);
    a = NULL;
  }
  /* S>0 pure singlet: the helper zeros everything; the bridge refuses. */
  {
    Task5Built built;
    MVMCKrylovStatus status;
    MVMCPowerLanczosGCView view;
    MVMCPowerLanczosGCBridge *bridge = NULL;
    int chain_rank = world_rank(), chain_size = world_size();
    int sym;
    CHECK(task5_build(&task5_r1, &built), "R1 build");
    built.support.pair_symmetry = MVMC_GC_PAIR_SINGLET;
    a = task5_workspace(&built, &status);
    CHECK(status == MVMC_KRYLOV_STATUS_OK && a != NULL &&
              evaluate(a, 6, &result) == MVMC_KRYLOV_STATUS_OK &&
              result.value.state == MVMC_SCALED_COMPLEX_EXACT_ZERO,
          "S>0 singlet descriptor zeros Psi(6)");
    mvmc_gc_krylov_amplitude_workspace_destroy(a);
    for (sym = MVMC_GC_PAIR_SINGLET; sym <= MVMC_GC_PAIR_ZERO;
         sym += MVMC_GC_PAIR_ZERO - MVMC_GC_PAIR_SINGLET) {
      built.support.pair_symmetry = (MVMCGCPairSymmetry)sym;
      task5_view(&built, world_rank(), chain_rank, chain_size, &view);
      CHECK(mvmc_power_lanczos_gc_bridge_create(&view, world_comm(),
                                                world_comm(), &bridge) ==
                    MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
                bridge == NULL,
            "S>0 singlet/zero bridge is an amplitude failure");
    }
    built.support.pair_symmetry = MVMC_GC_PAIR_TRIPLET;
    task5_view(&built, world_rank(), chain_rank, chain_size, &view);
    CHECK(mvmc_power_lanczos_gc_bridge_create(&view, world_comm(),
                                              world_comm(), &bridge) ==
                  MVMC_KRYLOV_STATUS_OK &&
              bridge != NULL,
          "R1 triplet bridge");
    mvmc_power_lanczos_gc_bridge_destroy(bridge);
    task5_free(&built);
  }
}

/* A nonfinite component under a structural zero is never hidden. */
static void test_task5_nonfinite_under_mask(void) {
  Task5Built built;
  MVMCKrylovStatus status;
  MVMCGCKrylovAmplitudeWorkspace *workspace;
  MVMCKrylovScaledAmplitudeResult result;
  CHECK(task5_build(&task5_s2_g2, &built), "S2_G2 build");
  /* 53 occupies orbitals {0,2,4,5}; entry (0,2) is used only there. */
  built.slater[0 * built.orbital_count + 2] = NAN;
  built.slater[2 * built.orbital_count + 0] = NAN;
  workspace = task5_workspace(&built, &status);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && workspace != NULL,
        "NaN Slater accepted at creation");
  if (workspace != NULL) {
    CHECK(evaluate(workspace, 53, &result) == MVMC_KRYLOV_STATUS_NONFINITE,
          "masked configuration still reports the NaN component");
    CHECK(evaluate(workspace, 60, &result) == MVMC_KRYLOV_STATUS_OK &&
              result.value.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO,
          "an unaffected configuration is evaluated");
#if defined(MVMC_HAVE_MPFR_REFINEMENT)
    {
      MVMCKrylovFinalStatePolicy policy;
      MVMCKrylovRefinementResult refined;
      const uint64_t current = 60, proposal = 53;
      CHECK(mvmc_krylov_final_state_policy_create(
                1, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC, 41,
                built.coefficient, 2, &policy) == MVMC_KRYLOV_STATUS_OK &&
                mvmc_gc_krylov_refine_final(workspace, &built.model, &policy,
                                            &current, &proposal, 1, 0.5,
                                            &refined) ==
                    MVMC_KRYLOV_STATUS_NONFINITE,
            "refined builder validates NaN even for a support zero");
    }
#endif
    mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
  }
  task5_free(&built);
}

/* Rank-agreed metadata, invalid values and hash coverage of every field. */
static void test_task5_metadata(void) {
  Task5Built built;
  MVMCKrylovStatus status;
  MVMCGCKrylovAmplitudeWorkspace *workspace;
  MVMCGCSpinSupport variants[6];
  uint64_t hashes[6];
  int field, k, j;
  CHECK(task5_build(&task5_r1_g2, &built), "R1_G2 build");
  /* Every invalid value on every rank. */
  {
    const MVMCGCSpinSupport bad[] = {
        {2, 2, 1, MVMC_GC_PAIR_TRIPLET}, {1, 0, 1, MVMC_GC_PAIR_TRIPLET},
        {1, -1, 1, MVMC_GC_PAIR_TRIPLET}, {1, 1, 1, MVMC_GC_PAIR_TRIPLET},
        {1, 2, -1, MVMC_GC_PAIR_TRIPLET}, {1, 2, 2, MVMC_GC_PAIR_TRIPLET},
        {1, 2, INT_MAX, MVMC_GC_PAIR_TRIPLET},
        {1, 2, 0, MVMC_GC_PAIR_ZERO}, {1, 2, 0, MVMC_GC_PAIR_SINGLET},
        {1, 2, 1, (MVMCGCPairSymmetry)7}, {0, 2, 0, MVMC_GC_PAIR_GENERAL},
        {0, 1, 1, MVMC_GC_PAIR_GENERAL}};
    for (k = 0; k < (int)(sizeof(bad) / sizeof(bad[0])); ++k) {
      built.support = bad[k];
      workspace = task5_workspace(&built, &status);
      CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
            "invalid descriptor rejected on every rank");
      mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
    }
  }
  {
    MVMCGCKrylovAmplitudeLayout layout;
    built.support = mvmc_gc_spin_support_inactive();
    task5_layout(&built, world_comm(), world_rank(), world_size(), &layout);
    layout.site_count = INT_MAX / 2 + 1;
    layout.projection.site_count = (size_t)layout.site_count;
    workspace = NULL;
    CHECK(mvmc_gc_krylov_amplitude_workspace_create(
              &layout, built.slater, built.weight, &workspace) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              workspace == NULL,
          "L overflow rejected before any allocation");
  }
  /* Hash covers active, G, S and the symmetry separately. */
  variants[0] = built.support = (MVMCGCSpinSupport){1, 2, 1, MVMC_GC_PAIR_GENERAL};
  variants[1] = (MVMCGCSpinSupport){1, 3, 1, MVMC_GC_PAIR_GENERAL};
  variants[2] = (MVMCGCSpinSupport){1, 2, 0, MVMC_GC_PAIR_GENERAL};
  variants[3] = (MVMCGCSpinSupport){1, 2, 1, MVMC_GC_PAIR_TRIPLET};
  variants[4] = (MVMCGCSpinSupport){1, 2, 1, MVMC_GC_PAIR_SINGLET};
  variants[5] = mvmc_gc_spin_support_inactive();
  for (k = 0; k < 6; ++k) {
    built.support = variants[k];
    workspace = task5_workspace(&built, &status);
    CHECK(status == MVMC_KRYLOV_STATUS_OK && workspace != NULL,
          "valid descriptor variant");
    hashes[k] = mvmc_gc_krylov_amplitude_generation_hash(workspace);
    mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
  }
  for (k = 0; k < 6; ++k) {
    for (j = k + 1; j < 6; ++j) {
      CHECK(hashes[k] != hashes[j], "every descriptor field enters the hash");
    }
  }
  /* One rank differs in one field: same error on every rank, no hang. */
  if (world_size() > 1) {
    for (field = 0; field < 4; ++field) {
      MVMCPowerLanczosGCView view;
      MVMCPowerLanczosGCBridge *bridge = NULL;
      built.support = (MVMCGCSpinSupport){1, 2, 1, MVMC_GC_PAIR_TRIPLET};
      if (world_rank() == 1) {
        if (field == 0) built.support = mvmc_gc_spin_support_inactive();
        if (field == 1) built.support.gauss_count = 3;
        if (field == 2) built.support.total_spin = 0;
        if (field == 3) built.support.pair_symmetry = MVMC_GC_PAIR_GENERAL;
      }
      workspace = task5_workspace(&built, &status);
      CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
            "one-rank descriptor mismatch rejected by the chain audit");
      mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
      task5_view(&built, world_rank(), world_rank(), world_size(), &view);
      CHECK(mvmc_power_lanczos_gc_bridge_create(&view, world_comm(),
                                                world_comm(), &bridge) ==
                    MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                bridge == NULL,
            "one-rank descriptor mismatch rejected by the world audit");
      /* Separate chains still share the world audit. */
      task5_view(&built, world_rank(), 0, 1, &view);
      CHECK(mvmc_power_lanczos_gc_bridge_create(&view, world_comm(),
                                                self_comm(), &bridge) ==
                    MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                bridge == NULL,
            "mismatch across single-rank chains rejected by the world audit");
    }
    built.support = (MVMCGCSpinSupport){1, 2, 1, MVMC_GC_PAIR_TRIPLET};
    if (world_rank() == 1) built.support.gauss_count = 0;
    workspace = task5_workspace(&built, &status);
    CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
          "one-rank invalid value agreed as INVALID_ARGUMENT");
  }
  CHECK(MVMC_POWER_LANCZOS_GC_BRIDGE_VERSION == UINT64_C(3),
        "bridge version 3");
  /* Support helper failure: invariant failure on every rank, both paths. */
  built.support = (MVMCGCSpinSupport){1, 2, 1, MVMC_GC_PAIR_TRIPLET};
  workspace = task5_workspace(&built, &status);
  CHECK(status == MVMC_KRYLOV_STATUS_OK && workspace != NULL, "R1_G2 workspace");
  if (workspace != NULL) {
    MVMCKrylovScaledAmplitudeResult result;
    const int target = world_size() > 1 ? world_size() - 1 : 0;
    mvmc_gc_krylov_amplitude_testing_set_fault(
        MVMC_GC_KRYLOV_AMPLITUDE_FAULT_SUPPORT, target);
    CHECK(evaluate(workspace, 6, &result) ==
              MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE,
          "support helper failure agreed on every rank (ordinary)");
#if defined(MVMC_HAVE_MPFR_REFINEMENT)
    {
      MVMCKrylovFinalStatePolicy policy;
      MVMCKrylovRefinementResult refined;
      const uint64_t current = 6, proposal = 5;
      CHECK(mvmc_krylov_final_state_policy_create(
                1, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC, 43,
                built.coefficient, 2, &policy) == MVMC_KRYLOV_STATUS_OK &&
                mvmc_gc_krylov_refine_final(workspace, &built.model, &policy,
                                            &current, &proposal, 1, 0.125,
                                            &refined) ==
                    MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE,
            "support helper failure agreed on every rank (refined)");
    }
#endif
    mvmc_gc_krylov_amplitude_testing_set_fault(
        MVMC_GC_KRYLOV_AMPLITUDE_FAULT_NONE, -1);
    CHECK(evaluate(workspace, 6, &result) == MVMC_KRYLOV_STATUS_OK &&
              result.value.state == MVMC_SCALED_COMPLEX_FINITE_NONZERO,
          "evaluation after clearing the support fault");
    mvmc_gc_krylov_amplitude_workspace_destroy(workspace);
  }
  task5_free(&built);
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
  test_spin_parity_rule();
  test_node_and_vacuum();
  test_invalid_layouts();
  test_partition();
  test_rank_faults();
  test_gc_bridge();
  test_task5_frozen_values();
  test_task5_mask_mutations();
  test_task5_private_snapshot();
  test_task5_classification();
  test_task5_nonfinite_under_mask();
  test_task5_metadata();
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
