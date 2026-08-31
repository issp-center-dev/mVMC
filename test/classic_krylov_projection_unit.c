#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "classic_krylov_projection.h"

static int failures = 0;

#define CHECK(condition, message)                                         \
  do {                                                                    \
    if (!(condition)) {                                                   \
      fprintf(stderr, "ClassicKrylovProjection_Unit FAIL: %s\n",        \
              (message));                                                 \
      ++failures;                                                         \
    }                                                                     \
  } while (0)

typedef struct {
  int gutzwiller[4];
  int jastrow[4][4];
  int spin_jastrow[4][4];
  int dh2[8];
  int dh4[16];
  const int *jastrow_rows[4];
  const int *spin_jastrow_rows[4];
  const int *dh2_rows[1];
  const int *dh4_rows[1];
  double complex parameters[32];
  MVMCClassicKrylovProjectionLayout layout;
} ProjectionFixture;

static void initialize_fixture(ProjectionFixture *fixture) {
  static const int pair_map[4][4] = {
      {-1, 0, 1, 2}, {0, -1, 3, 4},
      {1, 3, -1, 5}, {2, 4, 5, -1}};
  static const int dh2[8] = {3, 1, 0, 2, 1, 3, 2, 0};
  static const int dh4[16] = {
      2, 3, 1, 2, 3, 0, 2, 3,
      0, 1, 3, 0, 1, 2, 0, 1};
  static const double gutzwiller[4] = {1.0 / 8.0, 1.0 / 4.0,
                                        3.0 / 8.0, 1.0 / 2.0};
  static const int pair_sites[6][2] = {
      {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  int index;

  memset(fixture, 0, sizeof(*fixture));
  for (index = 0; index < 4; ++index) {
    fixture->gutzwiller[index] = index;
    memcpy(fixture->jastrow[index], pair_map[index],
           sizeof(fixture->jastrow[index]));
    memcpy(fixture->spin_jastrow[index], pair_map[index],
           sizeof(fixture->spin_jastrow[index]));
    fixture->jastrow_rows[index] = fixture->jastrow[index];
    fixture->spin_jastrow_rows[index] = fixture->spin_jastrow[index];
    fixture->parameters[index] = gutzwiller[index];
  }
  memcpy(fixture->dh2, dh2, sizeof(dh2));
  memcpy(fixture->dh4, dh4, sizeof(dh4));
  fixture->dh2_rows[0] = fixture->dh2;
  fixture->dh4_rows[0] = fixture->dh4;
  for (index = 0; index < 6; ++index) {
    const int left = pair_sites[index][0];
    const int right = pair_sites[index][1];
    const int forward = (right - left + 4) % 4;
    const int backward = (left - right + 4) % 4;
    const int distance = forward < backward ? forward : backward;
    fixture->parameters[4 + index] = (1.0 + distance) / 32.0;
    fixture->parameters[10 + index] =
        ((left + right) % 2 == 0 ? 1.0 : -1.0) *
        (1.0 + distance) / 48.0;
  }
  for (index = 0; index < 6; ++index) {
    fixture->parameters[16 + index] =
        (index % 2 == 0 ? 1.0 : -1.0) * (index + 1) / 128.0;
  }
  for (index = 0; index < 10; ++index) {
    fixture->parameters[22 + index] =
        (index % 2 == 0 ? 1.0 : -1.0) * (index + 1) / 256.0;
  }
  fixture->layout.site_count = 4;
  fixture->layout.nproj = 32;
  fixture->layout.ngutzwiller_idx = 4;
  fixture->layout.njastrow_idx = 6;
  fixture->layout.nspin_jastrow_idx = 6;
  fixture->layout.ndoublon_holon_2site_idx = 1;
  fixture->layout.ndoublon_holon_4site_idx = 1;
  fixture->layout.gutzwiller_idx = fixture->gutzwiller;
  fixture->layout.jastrow_idx = fixture->jastrow_rows;
  fixture->layout.spin_jastrow_idx = fixture->spin_jastrow_rows;
  fixture->layout.doublon_holon_2site_idx = fixture->dh2_rows;
  fixture->layout.doublon_holon_4site_idx = fixture->dh4_rows;
  fixture->layout.parameters = fixture->parameters;
}

static void set_occupation(uint64_t bits, int ele_num[8]) {
  int orbital;
  for (orbital = 0; orbital < 8; ++orbital) {
    ele_num[orbital] = (int)((bits >> orbital) & UINT64_C(1));
  }
}

static void test_combined_semantics_and_deep_copy(void) {
  static const int64_t expected_double_holon[32] = {
      1, 1, 0, 0, 1, -1, -1, -1, -1, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0,
      0, 0, 0, 0, 0, 0, 2, 2, 0, 0};
  static const int64_t expected_spin[32] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      -1, 1, -1, -1, 1, -1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  ProjectionFixture fixture;
  MVMCClassicKrylovProjectionWorkspace *workspace = NULL;
  int ele_num[8];
  int ele_copy[8];
  int64_t counts[32];
  double log_factor = -9.0;
  const unsigned char *binding;
  unsigned char *binding_copy;
  size_t binding_size = 0;

  initialize_fixture(&fixture);
  CHECK(mvmc_classic_krylov_projection_workspace_create(
            &fixture.layout, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "combined workspace create");
  CHECK(mvmc_classic_krylov_projection_workspace_bytes(workspace) >
            sizeof(fixture),
        "workspace byte accounting");
  binding = mvmc_classic_krylov_projection_binding_bytes(
      workspace, &binding_size);
  binding_copy = (unsigned char *)malloc(binding_size);
  CHECK(binding != NULL && binding_size != 0 && binding_copy != NULL,
        "binding record available");
  if (binding_copy != NULL) memcpy(binding_copy, binding, binding_size);

  set_occupation(UINT64_C(0x33), ele_num);
  memcpy(ele_copy, ele_num, sizeof(ele_num));
  CHECK(mvmc_classic_krylov_projection_evaluate(
            workspace, ele_num, 8, counts, 32, &log_factor) ==
            MVMC_KRYLOV_STATUS_OK,
        "double-holon combined evaluation");
  CHECK(memcmp(counts, expected_double_holon, sizeof(counts)) == 0,
        "all combined family counts match P6-A semantics");
  CHECK(fabs(log_factor - 21.0 / 128.0) < 1.0e-15,
        "double-holon exact rational log factor");
  CHECK(memcmp(ele_num, ele_copy, sizeof(ele_num)) == 0,
        "evaluation preserves occupation input bytes");
  counts[0] = 333;
  log_factor = 19.0;
  CHECK(mvmc_classic_krylov_projection_evaluate(
            workspace, ele_num, 8, counts, 31, &log_factor) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            counts[0] == 333 && log_factor == 19.0,
        "short count buffer preserves outputs");

  set_occupation(UINT64_C(0xa5), ele_num);
  CHECK(mvmc_classic_krylov_projection_evaluate(
            workspace, ele_num, 8, counts, 32, &log_factor) ==
            MVMC_KRYLOV_STATUS_OK,
        "spin-Jastrow evaluation");
  CHECK(memcmp(counts, expected_spin, sizeof(counts)) == 0,
        "spin-Jastrow signed counts");
  CHECK(fabs(log_factor - 7.0 / 24.0) < 1.0e-15,
        "spin-Jastrow exact rational log factor");

  fixture.parameters[0] = 99.0;
  fixture.gutzwiller[0] = 3;
  fixture.jastrow[0][1] = 5;
  fixture.dh2[0] = 0;
  set_occupation(UINT64_C(0x33), ele_num);
  CHECK(mvmc_classic_krylov_projection_evaluate(
            workspace, ele_num, 8, counts, 32, &log_factor) ==
            MVMC_KRYLOV_STATUS_OK &&
            memcmp(counts, expected_double_holon, sizeof(counts)) == 0 &&
            fabs(log_factor - 21.0 / 128.0) < 1.0e-15,
        "workspace is independent of caller binding mutation");
  binding = mvmc_classic_krylov_projection_binding_bytes(
      workspace, &binding_size);
  CHECK(binding_copy != NULL &&
            memcmp(binding, binding_copy, binding_size) == 0,
        "immutable canonical binding record");
  free(binding_copy);
  mvmc_classic_krylov_projection_workspace_destroy(workspace);
}

static void test_empty_and_transactional_errors(void) {
  MVMCClassicKrylovProjectionLayout layout;
  MVMCClassicKrylovProjectionWorkspace *workspace = NULL;
  int ele_num[8];
  int64_t sentinel[2] = {111, 222};
  double log_factor = 17.0;

  memset(&layout, 0, sizeof(layout));
  layout.site_count = 4;
  set_occupation(UINT64_C(0xa5), ele_num);
  CHECK(mvmc_classic_krylov_projection_workspace_create(
            &layout, &workspace) == MVMC_KRYLOV_STATUS_OK,
        "empty projection workspace create");
  CHECK(mvmc_classic_krylov_projection_evaluate(
            workspace, ele_num, 8, NULL, 0, &log_factor) ==
            MVMC_KRYLOV_STATUS_OK && log_factor == 0.0,
        "empty projection factor is identity");
  log_factor = 17.0;
  ele_num[2] = 2;
  CHECK(mvmc_classic_krylov_projection_evaluate(
            workspace, ele_num, 8, sentinel, 2, &log_factor) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            sentinel[0] == 111 && sentinel[1] == 222 &&
            log_factor == 17.0,
        "malformed occupation preserves outputs");
  mvmc_classic_krylov_projection_workspace_destroy(workspace);
}

static void test_invalid_bindings(void) {
  ProjectionFixture fixture;
  MVMCClassicKrylovProjectionWorkspace *workspace =
      (MVMCClassicKrylovProjectionWorkspace *)(uintptr_t)1;

  initialize_fixture(&fixture);
  fixture.layout.nproj = 31;
  CHECK(mvmc_classic_krylov_projection_workspace_create(
            &fixture.layout, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && workspace == NULL,
        "projection count mismatch rejected");

  initialize_fixture(&fixture);
  fixture.jastrow[0][1] = 5;
  CHECK(mvmc_classic_krylov_projection_workspace_create(
            &fixture.layout, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "asymmetric Jastrow map rejected");

  initialize_fixture(&fixture);
  fixture.spin_jastrow[0][0] = 0;
  CHECK(mvmc_classic_krylov_projection_workspace_create(
            &fixture.layout, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "non-sentinel Jastrow diagonal rejected");

  initialize_fixture(&fixture);
  fixture.dh4[3] = 4;
  CHECK(mvmc_classic_krylov_projection_workspace_create(
            &fixture.layout, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "out-of-range DH neighbor rejected");

  initialize_fixture(&fixture);
  fixture.parameters[7] = 1.0 + I;
  CHECK(mvmc_classic_krylov_projection_workspace_create(
            &fixture.layout, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "complex classic projection parameter rejected");

  initialize_fixture(&fixture);
  fixture.parameters[7] = NAN;
  CHECK(mvmc_classic_krylov_projection_workspace_create(
            &fixture.layout, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "nonfinite projection parameter rejected");

  initialize_fixture(&fixture);
  fixture.layout.ngutzwiller_idx = -1;
  CHECK(mvmc_classic_krylov_projection_workspace_create(
            &fixture.layout, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "negative correlator count rejected");

  initialize_fixture(&fixture);
  fixture.gutzwiller[2] = 4;
  CHECK(mvmc_classic_krylov_projection_workspace_create(
            &fixture.layout, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "out-of-range Gutzwiller index rejected");

  initialize_fixture(&fixture);
  fixture.layout.parameters = NULL;
  CHECK(mvmc_classic_krylov_projection_workspace_create(
            &fixture.layout, &workspace) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "missing parameter binding rejected");
}

int main(void) {
  test_combined_semantics_and_deep_copy();
  test_empty_and_transactional_errors();
  test_invalid_bindings();
  if (failures != 0) return 1;
  printf("ClassicKrylovProjection_Unit: PASS\n");
  return 0;
}
