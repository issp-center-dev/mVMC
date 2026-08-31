#include <stdio.h>
#include <string.h>

#include "classic_krylov_configuration.h"

static int failures = 0;

#define CHECK(condition, message)                                         \
  do {                                                                    \
    if (!(condition)) {                                                   \
      fprintf(stderr, "ClassicKrylovConfiguration_Unit FAIL: %s\n",     \
              (message));                                                 \
      ++failures;                                                         \
    }                                                                     \
  } while (0)

static MVMCClassicKrylovConfigurationLayout electronic_layout(void) {
  MVMCClassicKrylovConfigurationLayout layout = {4, 2, 2, 0};
  return layout;
}

static void test_electronic_round_trip(void) {
  const MVMCClassicKrylovConfigurationLayout layout = electronic_layout();
  const uint64_t input = UINT64_C(0xa5);
  uint64_t output = UINT64_C(0xffffffffffffffff);
  int ele_idx[4];
  int ele_cfg[8];
  int ele_num[8];

  CHECK(mvmc_classic_krylov_configuration_word_count(&layout) == 1,
        "electronic word count");
  CHECK(mvmc_classic_krylov_configuration_decode(
            &layout, &input, 1, ele_idx, ele_cfg, ele_num) ==
            MVMC_KRYLOV_STATUS_OK,
        "electronic decode");
  CHECK(ele_idx[0] == 0 && ele_idx[1] == 2 &&
            ele_idx[2] == 1 && ele_idx[3] == 3,
        "ascending per-spin particle order");
  CHECK(ele_cfg[0] == 0 && ele_cfg[2] == 1 &&
            ele_cfg[5] == 0 && ele_cfg[7] == 1,
        "classic particle labels");
  CHECK(mvmc_classic_krylov_configuration_encode(
            &layout, ele_idx, ele_cfg, ele_num, &output, 1) ==
            MVMC_KRYLOV_STATUS_OK && output == input,
        "electronic round trip");
}

static void test_pure_spin_contract(void) {
  MVMCClassicKrylovConfigurationLayout layout = electronic_layout();
  uint64_t valid = UINT64_C(0x5a);
  uint64_t invalid = UINT64_C(0x99);
  int ele_idx[4];
  int ele_cfg[8];
  int ele_num[8];
  layout.pure_spin = 1;
  CHECK(mvmc_classic_krylov_configuration_decode(
            &layout, &valid, 1, ele_idx, ele_cfg, ele_num) ==
            MVMC_KRYLOV_STATUS_OK,
        "pure-spin valid decode");
  CHECK(mvmc_classic_krylov_configuration_decode(
            &layout, &invalid, 1, ele_idx, ele_cfg, ele_num) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "pure-spin double/empty sites rejected");
}

static void test_transactional_failures(void) {
  MVMCClassicKrylovConfigurationLayout layout = electronic_layout();
  uint64_t words = UINT64_C(0x5a);
  uint64_t sentinel = UINT64_C(0x123456789abcdef0);
  int ele_idx[4] = {1, 3, 0, 2};
  int ele_cfg[8] = {-1, 0, -1, 1, 0, -1, 1, -1};
  int ele_num[8] = {0, 1, 0, 1, 1, 0, 1, 0};
  int out_idx[4] = {7, 7, 7, 7};
  int out_cfg[8] = {7, 7, 7, 7, 7, 7, 7, 7};
  int out_num[8] = {7, 7, 7, 7, 7, 7, 7, 7};

  ele_cfg[6] = 0;
  CHECK(mvmc_classic_krylov_configuration_encode(
            &layout, ele_idx, ele_cfg, ele_num, &sentinel, 1) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            sentinel == UINT64_C(0x123456789abcdef0),
        "malformed classic views preserve output");
  CHECK(mvmc_classic_krylov_configuration_decode(
            &layout, &words, 0, out_idx, out_cfg, out_num) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            out_idx[0] == 7 && out_cfg[0] == 7 && out_num[0] == 7,
        "decode failure preserves outputs");

  layout.down_electron_count = 1;
  CHECK(mvmc_classic_krylov_configuration_word_count(&layout) == 0,
        "nonzero 2Sz layout rejected");
  layout = electronic_layout();
  layout.site_count = -1;
  CHECK(mvmc_classic_krylov_configuration_word_count(&layout) == 0,
        "negative site count rejected");
  layout = electronic_layout();
  layout.up_electron_count = -1;
  layout.down_electron_count = -1;
  CHECK(mvmc_classic_krylov_configuration_word_count(&layout) == 0,
        "negative electron counts rejected");
  layout = electronic_layout();
  layout.up_electron_count = 0;
  layout.down_electron_count = 0;
  CHECK(mvmc_classic_krylov_configuration_word_count(&layout) == 0,
        "zero electron counts rejected");
}

static void test_multiword_padding(void) {
  MVMCClassicKrylovConfigurationLayout layout = {40, 1, 1, 0};
  uint64_t words[2] = {UINT64_C(1), UINT64_C(1) << 15};
  int ele_idx[2];
  int ele_cfg[80];
  int ele_num[80];

  CHECK(mvmc_classic_krylov_configuration_word_count(&layout) == 2,
        "multiword count");
  CHECK(mvmc_classic_krylov_configuration_decode(
            &layout, words, 2, ele_idx, ele_cfg, ele_num) ==
            MVMC_KRYLOV_STATUS_OK && ele_idx[0] == 0 && ele_idx[1] == 39,
        "multiword decode");
  words[1] |= UINT64_C(1) << 16;
  CHECK(mvmc_classic_krylov_configuration_decode(
            &layout, words, 2, ele_idx, ele_cfg, ele_num) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "high padding bit rejected");
}

int main(void) {
  test_electronic_round_trip();
  test_pure_spin_contract();
  test_transactional_failures();
  test_multiword_padding();
  if (failures != 0) return 1;
  printf("ClassicKrylovConfiguration_Unit: PASS\n");
  return 0;
}
