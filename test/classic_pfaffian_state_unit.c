#include "classic_pfaffian_state.h"

#include <complex.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static int failure_count = 0;

#define CHECK(condition, message)                                      \
  do {                                                                 \
    if (!(condition)) {                                                \
      fprintf(stderr, "FAIL: %s (line %d)\n", (message), __LINE__);   \
      ++failure_count;                                                 \
    }                                                                  \
  } while (0)

static int close_complex(double complex actual, double complex expected,
                         double tolerance) {
  return cabs(actual - expected) <= tolerance * fmax(1.0, cabs(expected));
}

static void set_real_pair(double *matrix, int n, int row, int column,
                          double value) {
  matrix[(size_t)row * (size_t)n + (size_t)column] = value;
  matrix[(size_t)column * (size_t)n + (size_t)row] = -value;
}

static void set_complex_pair(double complex *matrix, int n, int row,
                             int column, double complex value) {
  matrix[(size_t)row * (size_t)n + (size_t)column] = value;
  matrix[(size_t)column * (size_t)n + (size_t)row] = -value;
}

static void fill_real_slater(double *matrix, double scale) {
  memset(matrix, 0, 16 * sizeof(*matrix));
  set_real_pair(matrix, 4, 0, 1, 2.0 * scale);
  set_real_pair(matrix, 4, 0, 2, -1.0 * scale);
  set_real_pair(matrix, 4, 0, 3, 0.5 * scale);
  set_real_pair(matrix, 4, 1, 2, 3.0 * scale);
  set_real_pair(matrix, 4, 1, 3, 4.0 * scale);
  set_real_pair(matrix, 4, 2, 3, -2.0 * scale);
}

static void fill_complex_slater(double complex *matrix,
                                double complex scale) {
  memset(matrix, 0, 16 * sizeof(*matrix));
  set_complex_pair(matrix, 4, 0, 1, (1.0 + 0.5 * I) * scale);
  set_complex_pair(matrix, 4, 0, 2, (-0.25 + 0.75 * I) * scale);
  set_complex_pair(matrix, 4, 0, 3, (2.0 - 0.5 * I) * scale);
  set_complex_pair(matrix, 4, 1, 2, (0.5 + 1.25 * I) * scale);
  set_complex_pair(matrix, 4, 1, 3, (-1.0 + 0.25 * I) * scale);
  set_complex_pair(matrix, 4, 2, 3, (1.5 - 0.75 * I) * scale);
}

static void selected_real_matrix(const double *slater, const int *ele_idx,
                                 double *matrix) {
  int row, column;
  for (row = 0; row < 4; ++row) {
    const int slater_row = ele_idx[row] + (row / 2) * 2;
    for (column = 0; column < 4; ++column) {
      const int slater_column = ele_idx[column] + (column / 2) * 2;
      matrix[row * 4 + column] = slater[slater_row * 4 + slater_column];
    }
  }
}

static void selected_complex_matrix(const double complex *slater,
                                    const int *ele_idx,
                                    double complex *matrix) {
  int row, column;
  for (row = 0; row < 4; ++row) {
    const int slater_row = ele_idx[row] + (row / 2) * 2;
    for (column = 0; column < 4; ++column) {
      const int slater_column = ele_idx[column] + (column / 2) * 2;
      matrix[row * 4 + column] = slater[slater_row * 4 + slater_column];
    }
  }
}

static double complex pfaffian4(const double complex *matrix) {
  return matrix[1] * matrix[11] - matrix[2] * matrix[7] +
         matrix[3] * matrix[6];
}

static double real_pfaffian4(const double *matrix) {
  return matrix[1] * matrix[11] - matrix[2] * matrix[7] +
         matrix[3] * matrix[6];
}

static void check_real_inverse(const double *matrix,
                               const double *legacy_inverse,
                               const char *message) {
  int row, column, inner;
  for (row = 0; row < 4; ++row) {
    for (column = 0; column < 4; ++column) {
      double product = 0.0;
      for (inner = 0; inner < 4; ++inner) {
        product += matrix[row * 4 + inner] *
                   legacy_inverse[inner * 4 + column];
      }
      if (fabs(product - (row == column ? 1.0 : 0.0)) > 2.0e-12) {
        CHECK(0, message);
        return;
      }
    }
  }
}

static void check_complex_inverse(const double complex *matrix,
                                  const double complex *legacy_inverse,
                                  const char *message) {
  int row, column, inner;
  for (row = 0; row < 4; ++row) {
    for (column = 0; column < 4; ++column) {
      double complex product = 0.0;
      for (inner = 0; inner < 4; ++inner) {
        product += matrix[row * 4 + inner] *
                   legacy_inverse[inner * 4 + column];
      }
      if (cabs(product - (row == column ? 1.0 : 0.0)) > 3.0e-12) {
        CHECK(0, message);
        return;
      }
    }
  }
}

static void test_real_transaction_and_transitions(void) {
  const int ele_idx[4] = {1, 0, 0, 1};
  const double complex weights[1] = {1.0};
  double slater[16];
  double regular_slater[16];
  double selected[16];
  double legacy[17];
  double legacy_snapshot[17];
  MVMCClassicPfaffianRealWorkspace *workspace = NULL;
  const MVMCClassicPfaffianState *state;
  const double *raw_inverse;
  size_t index;

  fill_real_slater(slater, 1.0);
  memcpy(regular_slater, slater, sizeof(slater));
  selected_real_matrix(slater, ele_idx, selected);
  for (index = 0; index < 17; ++index) legacy[index] = 321.0;
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 4, 1, 0, 1, &workspace) == MVMC_PFAFFIAN_STATUS_OK,
        "create real workspace");
  CHECK(mvmc_classic_pfaffian_real_legacy_element_count(workspace) == 17,
        "real legacy element count includes PfM tail");
  memcpy(legacy_snapshot, legacy, sizeof(legacy));
  CHECK(mvmc_classic_pfaffian_real_publish(workspace, legacy, 17) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "real publish before prepare is rejected");
  CHECK(memcmp(legacy, legacy_snapshot, sizeof(legacy)) == 0 &&
            !mvmc_classic_pfaffian_real_accepted(workspace)->valid,
        "real publish before prepare preserves mirror and accepted state");
  CHECK(mvmc_classic_pfaffian_real_prepare(
            workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "prepare initial real candidate");
  state = mvmc_classic_pfaffian_real_candidate(workspace);
  CHECK(state != NULL && state->valid && state->generation == 1,
        "candidate generation one before publish");
  CHECK(!mvmc_classic_pfaffian_real_accepted(workspace)->valid,
        "prepare leaves accepted state invalid");
  CHECK(close_complex(state->components[0].pfaffian,
                      real_pfaffian4(selected), 2.0e-13),
        "real adapter matrix Pfaffian matches independent selection");
  CHECK(legacy[0] == 321.0 && legacy[16] == 321.0,
        "prepare does not touch legacy mirror");
  memcpy(legacy_snapshot, legacy, sizeof(legacy));
  CHECK(mvmc_classic_pfaffian_real_publish(workspace, legacy, 16) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "short real legacy allocation rejected");
  CHECK(memcmp(legacy, legacy_snapshot, sizeof(legacy)) == 0,
        "short publish leaves legacy bytes unchanged");
  CHECK(!mvmc_classic_pfaffian_real_accepted(workspace)->valid,
        "short publish leaves accepted state unchanged");
  CHECK(mvmc_classic_pfaffian_real_publish(workspace, legacy, 17) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "publish initial real candidate");
  state = mvmc_classic_pfaffian_real_accepted(workspace);
  raw_inverse = mvmc_classic_pfaffian_real_accepted_inverse(workspace);
  CHECK(state->valid && state->generation == 1,
        "real accepted generation one");
  CHECK(close_complex(legacy[16], state->components[0].pfaffian, 1.0e-14),
        "real PfM tail mirrors accepted Pfaffian");
  for (index = 0; index < 16; ++index) {
    CHECK(fabs(legacy[index] + raw_inverse[index]) < 2.0e-13,
          "real legacy flat inverse is negative raw LAPACK inverse");
  }
  check_real_inverse(selected, legacy,
                     "real row-major legacy mirror is matrix inverse");
  memcpy(legacy_snapshot, legacy, sizeof(legacy));
  CHECK(mvmc_classic_pfaffian_real_publish(workspace, legacy, 17) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "real duplicate publish is rejected");
  CHECK(memcmp(legacy, legacy_snapshot, sizeof(legacy)) == 0 &&
            mvmc_classic_pfaffian_real_accepted(workspace)->generation == 1,
        "real duplicate publish preserves mirror and accepted generation");

  fill_real_slater(slater, 0.75);
  selected_real_matrix(slater, ele_idx, selected);
  CHECK(mvmc_classic_pfaffian_real_prepare(
            workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "prepare regular to regular transition");
  CHECK(mvmc_classic_pfaffian_real_publish(workspace, legacy, 17) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "publish regular to regular transition");
  CHECK(mvmc_classic_pfaffian_real_accepted(workspace)->generation == 2 &&
            mvmc_classic_pfaffian_real_accepted(workspace)
                    ->components[0]
                    .state == MVMC_PFAFFIAN_REGULAR,
        "real regular to regular publish increments generation");
  check_real_inverse(selected, legacy,
                     "real regular to regular mirror is usable");

  memcpy(legacy_snapshot, legacy, sizeof(legacy));
  slater[1] = NAN;
  slater[4] = NAN;
  CHECK(mvmc_classic_pfaffian_real_prepare(
            workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "nonfinite real candidate rejected");
  CHECK(mvmc_classic_pfaffian_real_failure_local_qp(workspace) == 0 &&
            mvmc_classic_pfaffian_real_failure_global_qp(workspace) == 0,
        "real candidate failure reports local and global QP");
  CHECK(mvmc_classic_pfaffian_real_candidate(workspace) == NULL,
        "failed real candidate is not publishable");
  CHECK(mvmc_classic_pfaffian_real_accepted(workspace)->generation == 2,
        "failed real candidate preserves generation");
  CHECK(memcmp(legacy, legacy_snapshot, sizeof(legacy)) == 0,
        "failed real candidate preserves legacy mirror");

  memset(slater, 0, sizeof(slater));
  CHECK(mvmc_classic_pfaffian_real_prepare(
            workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "prepare regular to singular transition");
  CHECK(mvmc_classic_pfaffian_real_candidate(workspace)
                ->components[0]
                .state == MVMC_PFAFFIAN_SINGULAR,
        "real singular component classified");
  CHECK(mvmc_classic_pfaffian_real_publish(workspace, legacy, 17) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "publish real singular state");
  CHECK(mvmc_classic_pfaffian_real_accepted(workspace)->generation == 3,
        "real singular publish increments generation");
  CHECK(legacy[16] == 0.0, "real singular PfM mirror is exact zero");
  for (index = 0; index < 16; ++index) {
    CHECK(legacy[index] == 0.0, "real singular inverse mirror is zero");
  }

  CHECK(mvmc_classic_pfaffian_real_prepare(
            workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "prepare singular to singular transition");
  CHECK(mvmc_classic_pfaffian_real_publish(workspace, legacy, 17) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "publish singular to singular transition");
  CHECK(mvmc_classic_pfaffian_real_accepted(workspace)->generation == 4 &&
            mvmc_classic_pfaffian_real_accepted(workspace)
                    ->components[0]
                    .state == MVMC_PFAFFIAN_SINGULAR,
        "real singular to singular publish remains inverse-free");

  memcpy(slater, regular_slater, sizeof(slater));
  selected_real_matrix(slater, ele_idx, selected);
  CHECK(mvmc_classic_pfaffian_real_prepare(
            workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "prepare singular to regular transition");
  CHECK(mvmc_classic_pfaffian_real_publish(workspace, legacy, 17) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "publish singular to regular transition");
  CHECK(mvmc_classic_pfaffian_real_accepted(workspace)->generation == 5 &&
            mvmc_classic_pfaffian_real_accepted(workspace)
                    ->components[0]
                    .state == MVMC_PFAFFIAN_REGULAR,
        "real singular to regular rebuild restores inverse state");
  check_real_inverse(selected, legacy,
                     "real singular to regular mirror is usable");

  CHECK(mvmc_classic_pfaffian_real_prepare(
            workspace, slater, ele_idx, weights, 2.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "prepare forced near-singular state");
  CHECK(mvmc_classic_pfaffian_real_candidate(workspace)
                ->components[0]
                .state == MVMC_PFAFFIAN_NEAR_SINGULAR,
        "large pivot tolerance produces near-singular state");
  mvmc_classic_pfaffian_real_discard_candidate(workspace);
  CHECK(mvmc_classic_pfaffian_real_candidate(workspace) == NULL &&
            mvmc_classic_pfaffian_real_accepted(workspace)->generation == 5,
        "discard preserves accepted real state and generation");

  CHECK(mvmc_classic_pfaffian_real_prepare(
            workspace, slater, ele_idx, weights, 2.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "prepare near-singular candidate for publish");
  state = mvmc_classic_pfaffian_real_candidate(workspace);
  CHECK(state != NULL &&
            state->components[0].state == MVMC_PFAFFIAN_NEAR_SINGULAR &&
            isfinite(creal(state->components[0].pfaffian)) &&
            state->components[0].pfaffian != 0.0,
        "near-singular publish candidate has finite nonzero Pfaffian");
  CHECK(mvmc_classic_pfaffian_real_publish(workspace, legacy, 17) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "publish regular to near-singular transition");
  state = mvmc_classic_pfaffian_real_accepted(workspace);
  CHECK(state->generation == 6 &&
            state->components[0].state == MVMC_PFAFFIAN_NEAR_SINGULAR &&
            close_complex(legacy[16], state->components[0].pfaffian,
                          2.0e-13),
        "near-singular publish preserves finite Pfaffian and generation");
  for (index = 0; index < 16; ++index) {
    CHECK(legacy[index] == 0.0,
          "near-singular publish zeroes legacy inverse mirror");
  }

  mvmc_classic_pfaffian_real_workspace_destroy(workspace);
}

static void test_real_multiqup_and_alias_boundary(void) {
  const int ele_idx[4] = {1, 0, 0, 1};
  const double complex weights[4] = {1.0, 2.0, 3.0, 4.0};
  double slater[64];
  double selected[16];
  double legacy[68];
  double complex expected_total = 0.0;
  MVMCClassicPfaffianRealWorkspace *workspace = NULL;
  const MVMCClassicPfaffianState *candidate;
  int qp;
  size_t index;

  for (qp = 0; qp < 4; ++qp) {
    if (qp == 2) {
      memset(slater + qp * 16, 0, 16 * sizeof(*slater));
    } else {
      fill_real_slater(slater + qp * 16, (double)(qp + 1));
    }
  }
  for (index = 0; index < 68; ++index) legacy[index] = 777.0;
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 4, 4, 1, 4, &workspace) == MVMC_PFAFFIAN_STATUS_OK,
        "create real NQP4 local slice");
  CHECK(mvmc_classic_pfaffian_real_prepare(
            workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "prepare real NQP4 local slice");
  candidate = mvmc_classic_pfaffian_real_candidate(workspace);
  for (qp = 1; qp < 4; ++qp) {
    expected_total += weights[qp] * candidate->components[qp - 1].pfaffian;
  }
  CHECK(close_complex(candidate->local_aggregate.total, expected_total,
                      2.0e-13),
        "real slice aggregate uses global QP weights");
  CHECK(candidate->local_aggregate.singular_count == 1 &&
            cabs(candidate->local_aggregate.total) > 0.0,
        "one singular component preserves finite nonzero total");
  CHECK(mvmc_classic_pfaffian_real_publish(workspace, legacy, 68) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "publish real NQP4 local slice");
  for (index = 48; index < 64; ++index) {
    CHECK(legacy[index] == 777.0,
          "real publish leaves unused inverse capacity unchanged");
  }
  CHECK(legacy[67] == 777.0,
        "real publish leaves unused PfM tail capacity unchanged");
  for (qp = 0; qp < 3; ++qp) {
    if (qp == 1) {
      for (index = 0; index < 16; ++index) {
        CHECK(legacy[qp * 16 + index] == 0.0,
              "singular local QP inverse mirror is zero");
      }
    } else {
      selected_real_matrix(slater + (qp + 1) * 16, ele_idx, selected);
      check_real_inverse(selected, legacy + qp * 16,
                         "real local QP mirror uses matching global Slater");
    }
    CHECK(close_complex(legacy[64 + qp],
                        mvmc_classic_pfaffian_real_accepted(workspace)
                            ->components[qp]
                            .pfaffian,
                        2.0e-13),
          "real local PfM occupies local tail slot");
  }
  mvmc_classic_pfaffian_real_workspace_destroy(workspace);
}

static void test_complex_multiqup(void) {
  const int ele_idx[4] = {1, 0, 0, 1};
  const double complex weights[4] = {
      1.0, 0.5 + 0.25 * I, -0.75 + 0.5 * I, 1.25 - 0.5 * I};
  double complex slater[64];
  double complex selected[16];
  double complex legacy[68];
  double complex legacy_snapshot[68];
  double complex expected_total = 0.0;
  MVMCClassicPfaffianComplexWorkspace *workspace = NULL;
  const MVMCClassicPfaffianState *candidate;
  const double complex *raw_inverse;
  int qp;
  size_t index;

  for (qp = 0; qp < 4; ++qp) {
    fill_complex_slater(slater + qp * 16,
                        (double)(qp + 1) + 0.1 * (double)qp * I);
  }
  for (index = 0; index < 68; ++index) legacy[index] = 456.0 - 2.0 * I;
  CHECK(mvmc_classic_pfaffian_complex_workspace_create(
            2, 4, 4, 0, 4, &workspace) == MVMC_PFAFFIAN_STATUS_OK,
        "create complex NQP4 workspace");
  memcpy(legacy_snapshot, legacy, sizeof(legacy));
  CHECK(mvmc_classic_pfaffian_complex_publish(workspace, legacy, 68) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "complex publish before prepare is rejected");
  CHECK(memcmp(legacy, legacy_snapshot, sizeof(legacy)) == 0 &&
            !mvmc_classic_pfaffian_complex_accepted(workspace)->valid,
        "complex publish before prepare preserves mirror and accepted state");
  CHECK(mvmc_classic_pfaffian_complex_prepare(
            workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "prepare complex NQP4 candidate");
  candidate = mvmc_classic_pfaffian_complex_candidate(workspace);
  for (qp = 0; qp < 4; ++qp) {
    selected_complex_matrix(slater + qp * 16, ele_idx, selected);
    CHECK(close_complex(candidate->components[qp].pfaffian,
                        pfaffian4(selected), 3.0e-12),
          "complex adapter Pfaffian matches independent matrix selection");
    expected_total += weights[qp] * candidate->components[qp].pfaffian;
  }
  CHECK(close_complex(candidate->local_aggregate.total, expected_total,
                      3.0e-12),
        "complex aggregate includes phase and global weights");
  CHECK(mvmc_classic_pfaffian_complex_publish(workspace, legacy, 68) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "publish complex NQP4 candidate");
  raw_inverse = mvmc_classic_pfaffian_complex_accepted_inverse(workspace);
  for (qp = 0; qp < 4; ++qp) {
    selected_complex_matrix(slater + qp * 16, ele_idx, selected);
    check_complex_inverse(selected, legacy + qp * 16,
                          "complex row-major legacy mirror is inverse");
    for (index = 0; index < 16; ++index) {
      CHECK(cabs(legacy[qp * 16 + index] + raw_inverse[qp * 16 + index]) <
                3.0e-12,
            "complex legacy flat inverse is negative raw LAPACK inverse");
    }
    CHECK(close_complex(legacy[64 + qp],
                        mvmc_classic_pfaffian_complex_accepted(workspace)
                            ->components[qp]
                            .pfaffian,
                        3.0e-12),
          "complex PfM tail mirrors accepted phase");
  }
  memcpy(legacy_snapshot, legacy, sizeof(legacy));
  CHECK(mvmc_classic_pfaffian_complex_publish(workspace, legacy, 68) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "complex duplicate publish is rejected");
  CHECK(memcmp(legacy, legacy_snapshot, sizeof(legacy)) == 0 &&
            mvmc_classic_pfaffian_complex_accepted(workspace)->generation == 1,
        "complex duplicate publish preserves mirror and accepted generation");

  memcpy(legacy_snapshot, legacy, sizeof(legacy));
  slater[2 * 16 + 1] = NAN + I;
  CHECK(mvmc_classic_pfaffian_complex_prepare(
            workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "complex nonfinite candidate rejected");
  CHECK(mvmc_classic_pfaffian_complex_failure_local_qp(workspace) == 2 &&
            mvmc_classic_pfaffian_complex_failure_global_qp(workspace) == 2,
        "complex failure reports local and global QP");
  CHECK(mvmc_classic_pfaffian_complex_accepted(workspace)->generation == 1 &&
            memcmp(legacy, legacy_snapshot, sizeof(legacy)) == 0,
        "complex failure preserves accepted generation and mirror");

  memset(slater, 0, sizeof(slater));
  CHECK(mvmc_classic_pfaffian_complex_prepare(
            workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "prepare complex regular to singular transition");
  CHECK(mvmc_classic_pfaffian_complex_candidate(workspace)
                ->local_aggregate.singular_count == 4,
        "all complex components become singular");
  CHECK(mvmc_classic_pfaffian_complex_publish(workspace, legacy, 68) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "publish complex singular state");
  CHECK(mvmc_classic_pfaffian_complex_accepted(workspace)->generation == 2,
        "complex singular publish increments generation");
  for (index = 0; index < 68; ++index) {
    CHECK(legacy[index] == 0.0,
          "complex singular state zeroes local mirrors and PfM tail");
  }

  for (qp = 0; qp < 4; ++qp) {
    fill_complex_slater(slater + qp * 16,
                        (double)(qp + 1) + 0.1 * (double)qp * I);
  }
  CHECK(mvmc_classic_pfaffian_complex_prepare(
            workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "prepare complex singular to regular transition");
  CHECK(mvmc_classic_pfaffian_complex_publish(workspace, legacy, 68) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "publish complex singular to regular transition");
  CHECK(mvmc_classic_pfaffian_complex_accepted(workspace)->generation == 3 &&
            mvmc_classic_pfaffian_complex_accepted(workspace)
                    ->local_aggregate.regular_count == 4,
        "complex singular to regular rebuild restores all components");
  mvmc_classic_pfaffian_complex_workspace_destroy(workspace);
}

static void test_exact_total_cancellation(void) {
  const int ele_idx[4] = {1, 0, 0, 1};
  const double complex weights[2] = {1.0, -1.0};
  double slater[32];
  MVMCClassicPfaffianRealWorkspace *workspace = NULL;
  const MVMCClassicPfaffianState *candidate;

  fill_real_slater(slater, 1.0);
  fill_real_slater(slater + 16, 1.0);
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 4, 2, 0, 2, &workspace) == MVMC_PFAFFIAN_STATUS_OK,
        "create exact-cancellation workspace");
  CHECK(mvmc_classic_pfaffian_real_prepare(
            workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "exact-cancellation candidate remains numerically valid");
  candidate = mvmc_classic_pfaffian_real_candidate(workspace);
  CHECK(candidate != NULL && candidate->local_aggregate.valid &&
            candidate->local_aggregate.total == 0.0 &&
            candidate->local_aggregate.sum_abs > 0.0 &&
            candidate->local_aggregate.regular_count == 2,
        "regular components can produce exact-zero projected total");
  mvmc_classic_pfaffian_real_discard_candidate(workspace);
  CHECK(!mvmc_classic_pfaffian_real_accepted(workspace)->valid,
        "discarded zero-total proposal is not published");
  mvmc_classic_pfaffian_real_workspace_destroy(workspace);
}

static void test_empty_slice_and_invalid_inputs(void) {
  const int ele_idx[4] = {0, 1, 0, 1};
  const int negative_ele_idx[4] = {-1, 1, 0, 1};
  const int high_ele_idx[4] = {2, 1, 0, 1};
  double slater[16];
  double complex complex_slater[16];
  double complex weights[1] = {1.0};
  double legacy[17];
  double legacy_snapshot[17];
  MVMCClassicPfaffianRealWorkspace *real_workspace = NULL;
  MVMCClassicPfaffianComplexWorkspace *complex_workspace = NULL;
  size_t index;

  fill_real_slater(slater, 1.0);
  fill_complex_slater(complex_slater, 1.0);
  for (index = 0; index < 17; ++index) legacy[index] = 99.0;
  memcpy(legacy_snapshot, legacy, sizeof(legacy));
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 4, 1, 1, 1, &real_workspace) == MVMC_PFAFFIAN_STATUS_OK,
        "empty local slice workspace is valid");
  CHECK(mvmc_classic_pfaffian_real_prepare(
            real_workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "empty local slice prepare succeeds");
  CHECK(mvmc_classic_pfaffian_real_candidate(real_workspace)
                ->local_aggregate.valid &&
            mvmc_classic_pfaffian_real_candidate(real_workspace)
                    ->local_aggregate.total == 0.0,
        "empty slice contributes a valid neutral aggregate");
  CHECK(mvmc_classic_pfaffian_real_publish(real_workspace, legacy, 17) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "empty local slice publish succeeds");
  CHECK(memcmp(legacy, legacy_snapshot, sizeof(legacy)) == 0,
        "empty slice publish does not touch legacy allocation");
  CHECK(mvmc_classic_pfaffian_real_accepted(real_workspace)->generation == 1,
        "empty rank follows global generation");
  mvmc_classic_pfaffian_real_workspace_destroy(real_workspace);

  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 3, 1, 0, 1, &real_workspace) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "odd nsize rejected");
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 4, 1, 1, 0, &real_workspace) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "reversed QP range rejected");
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            INT_MAX / 2 + 1, 4, 1, 0, 1, &real_workspace) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "nsite doubling overflow rejected");
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 4, 1, 0, 1, &real_workspace) == MVMC_PFAFFIAN_STATUS_OK,
        "real invalid-weight workspace create");
  CHECK(mvmc_classic_pfaffian_real_prepare(
            real_workspace, slater, negative_ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "negative electron index rejected before Slater access");
  CHECK(mvmc_classic_pfaffian_real_prepare(
            real_workspace, slater, high_ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "electron index at nsite rejected before Slater access");
  CHECK(!mvmc_classic_pfaffian_real_accepted(real_workspace)->valid &&
            mvmc_classic_pfaffian_real_candidate(real_workspace) == NULL,
        "invalid electron indices leave real state unpublished");
  weights[0] = 1.0 + 1.0e-30 * I;
  CHECK(mvmc_classic_pfaffian_real_prepare(
            real_workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "real path rejects imaginary QP weight");
  weights[0] = DBL_MAX;
  CHECK(mvmc_classic_pfaffian_real_prepare(
            real_workspace, slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "real path rejects a nonfinite weighted component");
  CHECK(mvmc_classic_pfaffian_real_failure_local_qp(real_workspace) == 0 &&
            mvmc_classic_pfaffian_real_failure_global_qp(real_workspace) == 0,
        "weighted-component overflow reports local and global QP");
  mvmc_classic_pfaffian_real_workspace_destroy(real_workspace);

  weights[0] = NAN + I;
  CHECK(mvmc_classic_pfaffian_complex_workspace_create(
            2, 4, 1, 0, 1, &complex_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "complex invalid-weight workspace create");
  CHECK(mvmc_classic_pfaffian_complex_prepare(
            complex_workspace, complex_slater, ele_idx, weights, 0.0, 0.0) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "complex path rejects nonfinite QP weight");
  CHECK(mvmc_classic_pfaffian_complex_failure_global_qp(
            complex_workspace) == 0,
        "complex invalid weight reports global QP");
  mvmc_classic_pfaffian_complex_workspace_destroy(complex_workspace);
}

int main(void) {
  test_real_transaction_and_transitions();
  test_real_multiqup_and_alias_boundary();
  test_complex_multiqup();
  test_exact_total_cancellation();
  test_empty_slice_and_invalid_inputs();

  if (failure_count != 0) {
    fprintf(stderr, "classic Pfaffian state unit: %d failure(s)\n",
            failure_count);
    return 1;
  }
  printf("classic Pfaffian state unit checks passed\n");
  return 0;
}
