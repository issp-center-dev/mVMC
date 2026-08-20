#include "krylov_fock_proposal.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(condition, ...)                                                 \
  do {                                                                        \
    if (!(condition)) {                                                       \
      fprintf(stderr, "KrylovFockProposal_Unit FAIL: ");                    \
      fprintf(stderr, __VA_ARGS__);                                           \
      fprintf(stderr, " (line %d)\n", __LINE__);                            \
      ++failures;                                                             \
    }                                                                         \
  } while (0)

static MVMCKrylovFockModel model_with_sector(
    size_t site_count, size_t up_count, size_t down_count, int pure_spin) {
  MVMCKrylovFockModel model;
  memset(&model, 0, sizeof(model));
  model.site_count = site_count;
  model.up_electron_count = up_count;
  model.down_electron_count = down_count;
  model.pure_spin = pure_spin;
  model.hermitian = 1;
  return model;
}

static int close_double(double actual, double expected) {
  return fabs(actual - expected) <= 1.0e-14 * fmax(1.0, fabs(expected));
}

typedef struct {
  uint64_t state;
  size_t calls;
  size_t fail_at;
  int out_of_range;
} TestBoundedDraw;

static uint64_t test_splitmix64(uint64_t value) {
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

static MVMCKrylovStatus test_draw_bounded(
    void *context, size_t bound, size_t *value) {
  TestBoundedDraw *draw = (TestBoundedDraw *)context;
  uint64_t uint_bound;
  uint64_t threshold;
  uint64_t generated;
  if (draw == NULL || value == NULL || bound == 0 ||
      (size_t)(uint64_t)bound != bound) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (draw->calls == draw->fail_at) {
    return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  ++draw->calls;
  if (draw->out_of_range) {
    *value = bound;
    return MVMC_KRYLOV_STATUS_OK;
  }
  uint_bound = (uint64_t)bound;
  threshold = (uint64_t)(-uint_bound) % uint_bound;
  do {
    draw->state += UINT64_C(0x9e3779b97f4a7c15);
    generated = test_splitmix64(draw->state);
  } while (generated < threshold);
  *value = (size_t)(generated % uint_bound);
  return MVMC_KRYLOV_STATUS_OK;
}

static size_t census_state_index(
    uint64_t value, uint64_t *states, size_t *state_count,
    size_t state_capacity) {
  size_t index;
  for (index = 0; index < *state_count; ++index) {
    if (states[index] == value) return index;
  }
  if (*state_count >= state_capacity) return state_capacity;
  states[*state_count] = value;
  ++(*state_count);
  return *state_count - 1;
}

static size_t popcount_u64(uint64_t value) {
  size_t count = 0;
  while (value != 0) {
    count += (size_t)(value & UINT64_C(1));
    value >>= 1;
  }
  return count;
}

static size_t electronic_replacement_distance(uint64_t current,
                                              uint64_t proposal,
                                              size_t site_count) {
  const uint64_t site_mask =
      (UINT64_C(1) << site_count) - UINT64_C(1);
  const uint64_t current_up = current & site_mask;
  const uint64_t proposal_up = proposal & site_mask;
  const uint64_t current_down = (current >> site_count) & site_mask;
  const uint64_t proposal_down = (proposal >> site_count) & site_mask;
  return popcount_u64(current_up & ~proposal_up & site_mask) +
         popcount_u64(current_down & ~proposal_down & site_mask);
}

static size_t pure_spin_replacement_distance(uint64_t current,
                                             uint64_t proposal,
                                             size_t site_count) {
  const uint64_t site_mask =
      (UINT64_C(1) << site_count) - UINT64_C(1);
  return popcount_u64((current & site_mask) &
                      ~(proposal & site_mask) & site_mask);
}

static void test_electronic_neighbor_selection(void) {
  const MVMCKrylovFockModel model = model_with_sector(4, 2, 1, 0);
  const uint64_t current = UINT64_C(1) | UINT64_C(4) | UINT64_C(32);
  uint64_t proposal = 0;
  size_t neighbor_count = 0;

  CHECK(mvmc_krylov_fock_proposal_count_neighbors(
            &model, &current, 1, &neighbor_count) ==
            MVMC_KRYLOV_STATUS_OK &&
            neighbor_count == 7,
        "electronic neighbor count");
  CHECK(mvmc_krylov_fock_proposal_select_neighbor(
            &model, &current, 1, 2, &proposal, 1, &neighbor_count) ==
            MVMC_KRYLOV_STATUS_OK &&
            proposal == (UINT64_C(1) | UINT64_C(2) | UINT64_C(32)) &&
            neighbor_count == 7,
        "electronic selected up hop");
  CHECK(mvmc_krylov_fock_proposal_select_neighbor(
            &model, &current, 1, 7, &proposal, 1, &neighbor_count) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "out-of-range neighbor accepted");
}

static void test_pure_spin_exchange_selection(void) {
  const MVMCKrylovFockModel model = model_with_sector(4, 2, 2, 1);
  const uint64_t current = UINT64_C(5) | (UINT64_C(10) << 4);
  uint64_t proposal = 0;
  size_t neighbor_count = 0;

  CHECK(mvmc_krylov_fock_proposal_count_neighbors(
            &model, &current, 1, &neighbor_count) ==
            MVMC_KRYLOV_STATUS_OK &&
            neighbor_count == 4,
        "pure-spin exchange count");
  CHECK(mvmc_krylov_fock_proposal_select_neighbor(
            &model, &current, 1, 1, &proposal, 1, &neighbor_count) ==
            MVMC_KRYLOV_STATUS_OK &&
            proposal == (UINT64_C(12) | (UINT64_C(3) << 4)),
        "pure-spin selected exchange");
  CHECK(mvmc_krylov_fock_validate(&model, &proposal, 1) ==
            MVMC_KRYLOV_STATUS_OK,
        "pure-spin proposal is outside sector");
}

static void test_log_ratio_and_neighbor_contract(void) {
  const MVMCKrylovFockModel model = model_with_sector(4, 2, 1, 0);
  const uint64_t current = UINT64_C(1) | UINT64_C(4) | UINT64_C(32);
  uint64_t proposal = 0;
  double log_ratio = 99.0;
  size_t neighbor_count = 0;

  CHECK(mvmc_krylov_fock_proposal_select_neighbor(
            &model, &current, 1, 0, &proposal, 1, &neighbor_count) ==
            MVMC_KRYLOV_STATUS_OK,
        "proposal for log-ratio fixture");
  CHECK(mvmc_krylov_fock_proposal_log_ratio(
            &model, &current, 1, &proposal, 1, &log_ratio) ==
            MVMC_KRYLOV_STATUS_OK &&
            close_double(log_ratio, 0.0),
        "uniform-sector log proposal ratio");
  CHECK(mvmc_krylov_fock_proposal_log_ratio(
            &model, &current, 1, &current, 1, &log_ratio) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "self proposal accepted as neighbor");
}

static void test_multiword_padding_contract(void) {
  const MVMCKrylovFockModel model = model_with_sector(33, 1, 1, 0);
  uint64_t words[2] = {UINT64_C(1), UINT64_C(2)};
  size_t neighbor_count = 0;

  CHECK(mvmc_krylov_fock_proposal_count_neighbors(
            &model, words, 2, &neighbor_count) ==
            MVMC_KRYLOV_STATUS_OK &&
            neighbor_count == 64,
        "multiword proposal count");
  words[1] |= UINT64_C(4);
  CHECK(mvmc_krylov_fock_proposal_count_neighbors(
            &model, words, 2, &neighbor_count) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "multiword padding bit accepted");
}

static void test_uniform_electronic_sector(void) {
  const MVMCKrylovFockModel model = model_with_sector(4, 2, 1, 0);
  const size_t sample_count = 131072;
  const size_t expected_state_count = 24;
  const size_t expected_frequency = sample_count / expected_state_count;
  TestBoundedDraw draw = {
      UINT64_C(0x123456789abcdef0), 0, SIZE_MAX, 0};
  uint64_t states[24] = {0};
  size_t frequencies[24] = {0};
  size_t state_count = 0;
  size_t sample;

  for (sample = 0; sample < sample_count; ++sample) {
    uint64_t proposal = UINT64_MAX;
    MVMCKrylovFockUniformProposalResult result;
    size_t index;
    CHECK(mvmc_krylov_fock_proposal_draw_uniform_sector(
              &model, test_draw_bounded, &draw, &proposal, 1, &result) ==
              MVMC_KRYLOV_STATUS_OK &&
              result.valid && result.word_count == 1 &&
              result.active_orbital_count == 8 && result.draw_count <= 8 &&
              mvmc_krylov_fock_validate(&model, &proposal, 1) ==
                  MVMC_KRYLOV_STATUS_OK,
          "electronic uniform proposal sample %zu", sample);
    index = census_state_index(proposal, states, &state_count, 24);
    CHECK(index < 24, "electronic uniform census overflow");
    if (index < 24) ++frequencies[index];
  }
  CHECK(state_count == expected_state_count,
        "electronic uniform state count: %zu", state_count);
  for (sample = 0; sample < state_count; ++sample) {
    const size_t actual = frequencies[sample];
    const size_t difference = actual > expected_frequency
                                  ? actual - expected_frequency
                                  : expected_frequency - actual;
    CHECK(difference <= 700,
          "electronic state %zu frequency %zu outside fixed census",
          sample, actual);
  }
}

static void test_uniform_pure_spin_sector(void) {
  const MVMCKrylovFockModel model = model_with_sector(4, 2, 2, 1);
  const size_t sample_count = 65536;
  const size_t expected_state_count = 6;
  const size_t expected_frequency = sample_count / expected_state_count;
  TestBoundedDraw draw = {
      UINT64_C(0xfedcba9876543210), 0, SIZE_MAX, 0};
  uint64_t states[6] = {0};
  size_t frequencies[6] = {0};
  size_t state_count = 0;
  size_t sample;

  for (sample = 0; sample < sample_count; ++sample) {
    uint64_t proposal = UINT64_MAX;
    MVMCKrylovFockUniformProposalResult result;
    size_t index;
    CHECK(mvmc_krylov_fock_proposal_draw_uniform_sector(
              &model, test_draw_bounded, &draw, &proposal, 1, &result) ==
              MVMC_KRYLOV_STATUS_OK && result.valid &&
              mvmc_krylov_fock_validate(&model, &proposal, 1) ==
                  MVMC_KRYLOV_STATUS_OK,
          "pure-spin uniform proposal sample %zu", sample);
    CHECK(((proposal >> 4) & UINT64_C(0xf)) ==
              ((~proposal) & UINT64_C(0xf)),
          "pure-spin down subset is not the up complement");
    index = census_state_index(proposal, states, &state_count, 6);
    CHECK(index < 6, "pure-spin uniform census overflow");
    if (index < 6) ++frequencies[index];
  }
  CHECK(state_count == expected_state_count,
        "pure-spin uniform state count: %zu", state_count);
  for (sample = 0; sample < state_count; ++sample) {
    const size_t actual = frequencies[sample];
    const size_t difference = actual > expected_frequency
                                  ? actual - expected_frequency
                                  : expected_frequency - actual;
    CHECK(difference <= 800,
          "pure-spin state %zu frequency %zu outside fixed census",
          sample, actual);
  }
}

static void test_uniform_determinism_multiword_and_degenerate(void) {
  const MVMCKrylovFockModel multiword = model_with_sector(33, 1, 1, 0);
  const MVMCKrylovFockModel full_down = model_with_sector(4, 0, 4, 0);
  const MVMCKrylovFockModel one_spin = model_with_sector(1, 1, 0, 1);
  TestBoundedDraw first = {UINT64_C(0xabcdef), 0, SIZE_MAX, 0};
  TestBoundedDraw second = first;
  uint64_t first_words[2] = {UINT64_MAX, UINT64_MAX};
  uint64_t second_words[2] = {UINT64_MAX, UINT64_MAX};
  uint64_t degenerate = UINT64_MAX;
  MVMCKrylovFockUniformProposalResult first_result;
  MVMCKrylovFockUniformProposalResult second_result;

  CHECK(mvmc_krylov_fock_proposal_draw_uniform_sector(
            &multiword, test_draw_bounded, &first, first_words, 2,
            &first_result) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_fock_proposal_draw_uniform_sector(
                &multiword, test_draw_bounded, &second, second_words, 2,
                &second_result) == MVMC_KRYLOV_STATUS_OK &&
            memcmp(first_words, second_words, sizeof(first_words)) == 0 &&
            first_result.draw_count == second_result.draw_count &&
            first.calls == second.calls &&
            mvmc_krylov_fock_validate(&multiword, first_words, 2) ==
                MVMC_KRYLOV_STATUS_OK &&
            (first_words[1] & ~UINT64_C(3)) == 0,
        "multiword deterministic uniform proposal");

  first.state = UINT64_C(1);
  first.calls = 0;
  CHECK(mvmc_krylov_fock_proposal_draw_uniform_sector(
            &full_down, test_draw_bounded, &first, &degenerate, 1,
            &first_result) == MVMC_KRYLOV_STATUS_OK &&
            degenerate == UINT64_C(0xf0) && first_result.draw_count == 0,
        "full/empty electronic degenerate proposal");
  degenerate = UINT64_MAX;
  CHECK(mvmc_krylov_fock_proposal_draw_uniform_sector(
            &one_spin, test_draw_bounded, &first, &degenerate, 1,
            &first_result) == MVMC_KRYLOV_STATUS_OK &&
            degenerate == UINT64_C(1) && first_result.draw_count == 0,
        "one-state pure-spin degenerate proposal");
}

static void test_uniform_failure_contract(void) {
  MVMCKrylovFockModel model = model_with_sector(4, 2, 1, 0);
  TestBoundedDraw draw = {UINT64_C(7), 0, 1, 0};
  uint64_t proposal = UINT64_MAX;
  MVMCKrylovFockUniformProposalResult result;

  CHECK(mvmc_krylov_fock_proposal_draw_uniform_sector(
            &model, test_draw_bounded, &draw, &proposal, 1, &result) ==
            MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
            !result.valid && proposal == 0,
        "callback failure did not invalidate proposal");
  draw.calls = 0;
  draw.fail_at = SIZE_MAX;
  draw.out_of_range = 1;
  proposal = UINT64_MAX;
  CHECK(mvmc_krylov_fock_proposal_draw_uniform_sector(
            &model, test_draw_bounded, &draw, &proposal, 1, &result) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !result.valid && proposal == 0,
        "out-of-range callback value accepted");
  draw.out_of_range = 0;
  proposal = UINT64_MAX;
  CHECK(mvmc_krylov_fock_proposal_draw_uniform_sector(
            &model, test_draw_bounded, &draw, &proposal, 0, &result) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !result.valid,
        "wrong proposal word count accepted");
  model.up_electron_count = 5;
  proposal = UINT64_MAX;
  CHECK(mvmc_krylov_fock_proposal_draw_uniform_sector(
            &model, test_draw_bounded, &draw, &proposal, 1, &result) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !result.valid && proposal == 0,
        "invalid particle count accepted");
  model = model_with_sector(4, 2, 1, 1);
  proposal = UINT64_MAX;
  CHECK(mvmc_krylov_fock_proposal_draw_uniform_sector(
            &model, test_draw_bounded, &draw, &proposal, 1, &result) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !result.valid && proposal == 0,
        "invalid pure-spin counts accepted");
  model = model_with_sector(4, 2, 1, 0);
  model.hermitian = 0;
  proposal = UINT64_MAX;
  CHECK(mvmc_krylov_fock_proposal_draw_uniform_sector(
            &model, test_draw_bounded, &draw, &proposal, 1, &result) ==
            MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL &&
            !result.valid && proposal == 0,
        "non-Hermitian model accepted");
}

static void test_shell_distance_resolution_and_cardinality(void) {
  const MVMCKrylovFockModel electronic = model_with_sector(4, 2, 1, 0);
  const MVMCKrylovFockModel half4 = model_with_sector(4, 2, 2, 0);
  const MVMCKrylovFockModel half6 = model_with_sector(6, 3, 3, 0);
  const MVMCKrylovFockModel half8 = model_with_sector(8, 4, 4, 0);
  const MVMCKrylovFockModel pure4 = model_with_sector(4, 2, 2, 1);
  const MVMCKrylovFockModel frozen = model_with_sector(4, 4, 4, 0);
  size_t maximum = 0;
  size_t distance = 0;
  size_t count = 0;

  CHECK(mvmc_krylov_fock_proposal_max_shell_distance(
            &electronic, &maximum) == MVMC_KRYLOV_STATUS_OK &&
            maximum == 3,
        "electronic maximum shell distance");
  CHECK(mvmc_krylov_fock_proposal_count_shell(&electronic, 1, &count) ==
            MVMC_KRYLOV_STATUS_OK && count == 7,
        "electronic d=1 cardinality");
  CHECK(mvmc_krylov_fock_proposal_count_shell(&electronic, 2, &count) ==
            MVMC_KRYLOV_STATUS_OK && count == 13,
        "electronic d=2 cardinality");
  CHECK(mvmc_krylov_fock_proposal_count_shell(&electronic, 3, &count) ==
            MVMC_KRYLOV_STATUS_OK && count == 3,
        "electronic d=3 cardinality");
  CHECK(mvmc_krylov_fock_proposal_count_shell(&pure4, 1, &count) ==
            MVMC_KRYLOV_STATUS_OK && count == 4,
        "pure-spin d=1 cardinality");
  CHECK(mvmc_krylov_fock_proposal_count_shell(&pure4, 2, &count) ==
            MVMC_KRYLOV_STATUS_OK && count == 1,
        "pure-spin d=2 cardinality");

  CHECK(mvmc_krylov_fock_proposal_resolve_shell_distance(
            &half4, 1, 3, &maximum, &distance) ==
            MVMC_KRYLOV_STATUS_OK && maximum == 4 && distance == 1,
        "L4 one-third distance");
  CHECK(mvmc_krylov_fock_proposal_resolve_shell_distance(
            &half4, 1, 2, &maximum, &distance) ==
            MVMC_KRYLOV_STATUS_OK && distance == 2,
        "L4 one-half distance");
  CHECK(mvmc_krylov_fock_proposal_resolve_shell_distance(
            &half4, 2, 3, &maximum, &distance) ==
            MVMC_KRYLOV_STATUS_OK && distance == 3,
        "L4 two-thirds distance");
  CHECK(mvmc_krylov_fock_proposal_resolve_shell_distance(
            &half6, 1, 3, &maximum, &distance) ==
            MVMC_KRYLOV_STATUS_OK && maximum == 6 && distance == 2,
        "L6 one-third distance");
  CHECK(mvmc_krylov_fock_proposal_resolve_shell_distance(
            &half6, 1, 2, &maximum, &distance) ==
            MVMC_KRYLOV_STATUS_OK && distance == 3,
        "L6 one-half distance");
  CHECK(mvmc_krylov_fock_proposal_resolve_shell_distance(
            &half6, 2, 3, &maximum, &distance) ==
            MVMC_KRYLOV_STATUS_OK && distance == 4,
        "L6 two-thirds distance");
  CHECK(mvmc_krylov_fock_proposal_resolve_shell_distance(
            &half8, 1, 3, &maximum, &distance) ==
            MVMC_KRYLOV_STATUS_OK && maximum == 8 && distance == 3,
        "L8 one-third distance");
  CHECK(mvmc_krylov_fock_proposal_resolve_shell_distance(
            &half8, 1, 2, &maximum, &distance) ==
            MVMC_KRYLOV_STATUS_OK && distance == 4,
        "L8 one-half distance");
  CHECK(mvmc_krylov_fock_proposal_resolve_shell_distance(
            &half8, 2, 3, &maximum, &distance) ==
            MVMC_KRYLOV_STATUS_OK && distance == 5,
        "L8 two-thirds distance");
  CHECK(mvmc_krylov_fock_proposal_resolve_shell_distance(
            &pure4, 1, 3, &maximum, &distance) ==
            MVMC_KRYLOV_STATUS_OK && maximum == 2 && distance == 1 &&
            mvmc_krylov_fock_proposal_resolve_shell_distance(
                &pure4, 1, 2, &maximum, &count) ==
                MVMC_KRYLOV_STATUS_OK && count == 1 &&
            mvmc_krylov_fock_proposal_resolve_shell_distance(
                &pure4, 2, 3, &maximum, &count) ==
                MVMC_KRYLOV_STATUS_OK && count == 1,
        "pure-spin duplicate resolved kernels");
  CHECK(mvmc_krylov_fock_proposal_resolve_shell_distance(
            &frozen, 1, 2, &maximum, &distance) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT && maximum == 0 &&
            distance == 0,
        "Dmax=0 distance accepted");
  CHECK(mvmc_krylov_fock_proposal_resolve_shell_distance(
            &half4, 0, 3, &maximum, &distance) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            mvmc_krylov_fock_proposal_resolve_shell_distance(
                &half4, 4, 3, &maximum, &distance) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "invalid shell distance fraction accepted");
}

static void test_shell_exhaustive_reverse_pair_census(void) {
  const MVMCKrylovFockModel electronic = model_with_sector(4, 2, 1, 0);
  const MVMCKrylovFockModel pure_spin = model_with_sector(4, 2, 2, 1);
  uint64_t electronic_states[24];
  uint64_t pure_states[6];
  size_t electronic_count = 0;
  size_t pure_count = 0;
  uint64_t up;

  for (up = 0; up < UINT64_C(16); ++up) {
    uint64_t down;
    if (popcount_u64(up) != 2) continue;
    for (down = 0; down < UINT64_C(16); ++down) {
      if (popcount_u64(down) == 1) {
        electronic_states[electronic_count++] = up | (down << 4);
      }
    }
    pure_states[pure_count++] = up | (((~up) & UINT64_C(15)) << 4);
  }
  CHECK(electronic_count == 24 && pure_count == 6,
        "L4 sector enumeration counts");
  {
    size_t source;
    for (source = 0; source < electronic_count; ++source) {
      size_t distance;
      for (distance = 1; distance <= 3; ++distance) {
        size_t target;
        size_t observed = 0;
        size_t expected = 0;
        CHECK(mvmc_krylov_fock_proposal_count_shell(
                  &electronic, distance, &expected) ==
                  MVMC_KRYLOV_STATUS_OK,
              "electronic shell count for reverse census");
        for (target = 0; target < electronic_count; ++target) {
          const size_t forward = electronic_replacement_distance(
              electronic_states[source], electronic_states[target], 4);
          const size_t reverse = electronic_replacement_distance(
              electronic_states[target], electronic_states[source], 4);
          CHECK(forward == reverse,
                "electronic reverse distance mismatch %zu -> %zu", source,
                target);
          if (forward == distance) ++observed;
        }
        CHECK(observed == expected,
              "electronic shell census source=%zu d=%zu: %zu != %zu",
              source, distance, observed, expected);
      }
    }
  }
  {
    size_t source;
    for (source = 0; source < pure_count; ++source) {
      size_t distance;
      for (distance = 1; distance <= 2; ++distance) {
        size_t target;
        size_t observed = 0;
        size_t expected = 0;
        CHECK(mvmc_krylov_fock_proposal_count_shell(
                  &pure_spin, distance, &expected) ==
                  MVMC_KRYLOV_STATUS_OK,
              "pure-spin shell count for reverse census");
        for (target = 0; target < pure_count; ++target) {
          const size_t forward = pure_spin_replacement_distance(
              pure_states[source], pure_states[target], 4);
          const size_t reverse = pure_spin_replacement_distance(
              pure_states[target], pure_states[source], 4);
          CHECK(forward == reverse,
                "pure-spin reverse distance mismatch %zu -> %zu", source,
                target);
          if (forward == distance) ++observed;
        }
        CHECK(observed == expected,
              "pure-spin shell census source=%zu d=%zu: %zu != %zu",
              source, distance, observed, expected);
      }
    }
  }
}

static void test_pure_spin_l6_l8_shell_anchors(void) {
  const size_t sites[] = {6, 8};
  const size_t expected_sector_counts[] = {20, 70};
  const size_t numerators[] = {1, 1, 2};
  const size_t denominators[] = {3, 2, 3};
  size_t case_index;
  for (case_index = 0; case_index < 2; ++case_index) {
    const size_t site_count = sites[case_index];
    const size_t half = site_count / 2;
    const uint64_t site_mask =
        (UINT64_C(1) << site_count) - UINT64_C(1);
    const MVMCKrylovFockModel model =
        model_with_sector(site_count, half, half, 1);
    uint64_t states[70] = {0};
    size_t state_count = 0;
    uint64_t up_mask;
    size_t candidate;
    MVMCKrylovFockProposalConnectivity connectivity;

    for (up_mask = 0; up_mask <= site_mask; ++up_mask) {
      if (popcount_u64(up_mask) == half) {
        states[state_count++] =
            up_mask | ((site_mask ^ up_mask) << site_count);
      }
    }
    CHECK(state_count == expected_sector_counts[case_index],
          "pure-spin L%zu sector count %zu", site_count, state_count);
    CHECK(mvmc_krylov_fock_proposal_check_connectivity(
              &model, 100, &connectivity) == MVMC_KRYLOV_STATUS_OK &&
              connectivity.valid && connectivity.connected &&
              connectivity.sector_count == state_count &&
              connectivity.visited_count == state_count,
          "pure-spin L%zu neighbor connectivity", site_count);

    for (candidate = 0; candidate < 3; ++candidate) {
      size_t maximum = 0;
      size_t distance = 0;
      size_t expected = 0;
      size_t source;
      CHECK(mvmc_krylov_fock_proposal_resolve_shell_distance(
                &model, numerators[candidate], denominators[candidate],
                &maximum, &distance) == MVMC_KRYLOV_STATUS_OK &&
                maximum == half && distance >= 1 && distance <= maximum &&
                mvmc_krylov_fock_proposal_count_shell(
                    &model, distance, &expected) == MVMC_KRYLOV_STATUS_OK,
            "pure-spin L%zu candidate %zu distance/cardinality", site_count,
            candidate);
      for (source = 0; source < state_count; ++source) {
        size_t target;
        size_t observed = 0;
        for (target = 0; target < state_count; ++target) {
          const size_t forward = pure_spin_replacement_distance(
              states[source], states[target], site_count);
          const size_t reverse = pure_spin_replacement_distance(
              states[target], states[source], site_count);
          CHECK(forward == reverse,
                "pure-spin L%zu reverse distance %zu -> %zu", site_count,
                source, target);
          if (forward == distance) ++observed;
        }
        CHECK(observed == expected,
              "pure-spin L%zu shell source=%zu d=%zu: %zu != %zu",
              site_count, source, distance, observed, expected);
      }
    }
  }
}

static void test_shell_uniform_draws(void) {
  const MVMCKrylovFockModel electronic = model_with_sector(4, 2, 1, 0);
  const MVMCKrylovFockModel pure_spin = model_with_sector(4, 2, 2, 1);
  const uint64_t electronic_current = UINT64_C(0x13);
  const uint64_t pure_current = UINT64_C(0xc3);
  const size_t electronic_samples = 131072;
  const size_t pure_samples = 65536;
  TestBoundedDraw electronic_draw = {
      UINT64_C(0x735337656c656331), 0, SIZE_MAX, 0};
  TestBoundedDraw pure_draw = {
      UINT64_C(0x7353377075726531), 0, SIZE_MAX, 0};
  uint64_t electronic_states[13] = {0};
  size_t electronic_frequencies[13] = {0};
  uint64_t pure_states[4] = {0};
  size_t pure_frequencies[4] = {0};
  size_t electronic_state_count = 0;
  size_t pure_state_count = 0;
  size_t sample;

  for (sample = 0; sample < electronic_samples; ++sample) {
    uint64_t proposal = UINT64_MAX;
    MVMCKrylovFockShellProposalResult result;
    const MVMCKrylovStatus status = mvmc_krylov_fock_proposal_draw_shell(
        &electronic, &electronic_current, 1, 2, test_draw_bounded,
        &electronic_draw, &proposal, 1, &result);
    size_t index;
    CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid &&
              result.max_distance == 3 && result.distance == 2 &&
              result.up_distance + result.down_distance == 2 &&
              result.shell_count == 13 &&
              electronic_replacement_distance(electronic_current, proposal,
                                              4) == 2 &&
              mvmc_krylov_fock_validate(&electronic, &proposal, 1) ==
                  MVMC_KRYLOV_STATUS_OK,
          "electronic shell draw %zu", sample);
    index = census_state_index(proposal, electronic_states,
                               &electronic_state_count, 13);
    CHECK(index < 13, "electronic shell census overflow");
    if (index < 13) ++electronic_frequencies[index];
  }
  CHECK(electronic_state_count == 13,
        "electronic shell observed state count %zu", electronic_state_count);
  for (sample = 0; sample < electronic_state_count; ++sample) {
    const size_t expected = electronic_samples / 13;
    const size_t actual = electronic_frequencies[sample];
    const size_t difference =
        actual > expected ? actual - expected : expected - actual;
    CHECK(difference <= 700,
          "electronic shell state %zu frequency %zu", sample, actual);
  }

  for (sample = 0; sample < pure_samples; ++sample) {
    uint64_t proposal = UINT64_MAX;
    MVMCKrylovFockShellProposalResult result;
    const MVMCKrylovStatus status = mvmc_krylov_fock_proposal_draw_shell(
        &pure_spin, &pure_current, 1, 1, test_draw_bounded, &pure_draw,
        &proposal, 1, &result);
    size_t index;
    CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid &&
              result.max_distance == 2 && result.distance == 1 &&
              result.up_distance == 1 && result.down_distance == 1 &&
              result.shell_count == 4 &&
              pure_spin_replacement_distance(pure_current, proposal, 4) ==
                  1 &&
              ((proposal >> 4) & UINT64_C(15)) ==
                  ((~proposal) & UINT64_C(15)) &&
              mvmc_krylov_fock_validate(&pure_spin, &proposal, 1) ==
                  MVMC_KRYLOV_STATUS_OK,
          "pure-spin shell draw %zu", sample);
    index = census_state_index(proposal, pure_states, &pure_state_count, 4);
    CHECK(index < 4, "pure-spin shell census overflow");
    if (index < 4) ++pure_frequencies[index];
  }
  CHECK(pure_state_count == 4, "pure shell observed state count %zu",
        pure_state_count);
  for (sample = 0; sample < pure_state_count; ++sample) {
    const size_t expected = pure_samples / 4;
    const size_t actual = pure_frequencies[sample];
    const size_t difference =
        actual > expected ? actual - expected : expected - actual;
    CHECK(difference <= 600, "pure shell state %zu frequency %zu", sample,
          actual);
  }
}

static void test_shell_multiword_failure_and_overflow(void) {
  const MVMCKrylovFockModel multiword = model_with_sector(33, 1, 1, 0);
  MVMCKrylovFockModel huge =
      model_with_sector(SIZE_MAX / 2, SIZE_MAX / 4, 0, 0);
  const uint64_t current[2] = {UINT64_C(1) | (UINT64_C(1) << 33), 0};
  uint64_t proposal[2] = {UINT64_MAX, UINT64_MAX};
  uint64_t overlap[3] = {current[0], current[1], UINT64_C(0x55)};
  TestBoundedDraw draw = {UINT64_C(0x7353346661696c31), 0, SIZE_MAX, 0};
  MVMCKrylovFockShellProposalResult result;
  size_t count = 0;
  size_t maximum = 0;
  size_t distance = 0;

  CHECK(mvmc_krylov_fock_proposal_draw_shell(
            &multiword, current, 2, 2, test_draw_bounded, &draw, proposal, 2,
            &result) == MVMC_KRYLOV_STATUS_OK && result.valid &&
            result.shell_count == 1024 && result.distance == 2 &&
            result.up_distance == 1 && result.down_distance == 1 &&
            mvmc_krylov_fock_validate(&multiword, proposal, 2) ==
                MVMC_KRYLOV_STATUS_OK,
        "multiword shell draw");

  draw.calls = 0;
  draw.fail_at = 0;
  proposal[0] = UINT64_MAX;
  proposal[1] = UINT64_MAX;
  CHECK(mvmc_krylov_fock_proposal_draw_shell(
            &multiword, current, 2, 2, test_draw_bounded, &draw, proposal, 2,
            &result) == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
            !result.valid && proposal[0] == 0 && proposal[1] == 0,
        "shell callback failure did not clear output");
  draw.calls = 0;
  draw.fail_at = SIZE_MAX;
  draw.out_of_range = 1;
  proposal[0] = UINT64_MAX;
  proposal[1] = UINT64_MAX;
  CHECK(mvmc_krylov_fock_proposal_draw_shell(
            &multiword, current, 2, 2, test_draw_bounded, &draw, proposal, 2,
            &result) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            !result.valid && proposal[0] == 0 && proposal[1] == 0,
        "shell out-of-range callback accepted");
  draw.out_of_range = 0;
  proposal[0] = UINT64_MAX;
  proposal[1] = UINT64_MAX;
  CHECK(mvmc_krylov_fock_proposal_draw_shell(
            &multiword, current, 2, 0, test_draw_bounded, &draw, proposal, 2,
            &result) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            proposal[0] == 0 && proposal[1] == 0,
        "d=0 shell accepted");

  CHECK(mvmc_krylov_fock_proposal_draw_shell(
            &multiword, overlap, 2, 1, test_draw_bounded, &draw, overlap + 1,
            2, &result) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            overlap[0] == current[0] && overlap[1] == current[1] &&
            overlap[2] == UINT64_C(0x55),
        "overlapping shell buffers changed input");
  CHECK(mvmc_krylov_fock_proposal_count_shell(&huge, 2, &count) ==
            MVMC_KRYLOV_STATUS_RESOURCE_LIMIT && count == 0,
        "overflowing shell cardinality accepted");
  CHECK(mvmc_krylov_fock_proposal_resolve_shell_distance(
            &multiword, SIZE_MAX, SIZE_MAX, &maximum, &distance) ==
            MVMC_KRYLOV_STATUS_RESOURCE_LIMIT && maximum == 0 &&
            distance == 0,
        "overflowing rational distance arithmetic accepted");
}

static void test_connectivity(void) {
  const MVMCKrylovFockModel electronic = model_with_sector(4, 2, 1, 0);
  const MVMCKrylovFockModel pure_spin = model_with_sector(4, 2, 2, 1);
  const MVMCKrylovFockModel too_large = model_with_sector(33, 1, 1, 0);
  MVMCKrylovFockProposalConnectivity result;

  CHECK(mvmc_krylov_fock_proposal_check_connectivity(
            &electronic, 64, &result) == MVMC_KRYLOV_STATUS_OK &&
            result.valid && result.connected && result.sector_count == 24 &&
            result.visited_count == 24,
        "electronic sector connectivity");
  CHECK(mvmc_krylov_fock_proposal_check_connectivity(
            &pure_spin, 16, &result) == MVMC_KRYLOV_STATUS_OK &&
            result.valid && result.connected && result.sector_count == 6 &&
            result.visited_count == 6,
        "pure-spin sector connectivity");
  CHECK(mvmc_krylov_fock_proposal_check_connectivity(
            &electronic, 8, &result) ==
            MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
            !result.valid,
        "connectivity max-state gate ignored");
  CHECK(mvmc_krylov_fock_proposal_check_connectivity(
            &too_large, 128, &result) ==
            MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
            !result.valid,
        "multiword exact connectivity gate ignored");
}

int main(void) {
  test_electronic_neighbor_selection();
  test_pure_spin_exchange_selection();
  test_log_ratio_and_neighbor_contract();
  test_multiword_padding_contract();
  test_uniform_electronic_sector();
  test_uniform_pure_spin_sector();
  test_uniform_determinism_multiword_and_degenerate();
  test_uniform_failure_contract();
  test_shell_distance_resolution_and_cardinality();
  test_shell_exhaustive_reverse_pair_census();
  test_pure_spin_l6_l8_shell_anchors();
  test_shell_uniform_draws();
  test_shell_multiword_failure_and_overflow();
  test_connectivity();
  return failures == 0 ? 0 : 1;
}
