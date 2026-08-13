#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "krylov_fock_reference.h"

static int failures = 0;

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "AbsoluteKrylovReference_Unit FAIL: ");                  \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, "\n");                                                   \
      ++failures;                                                               \
    }                                                                           \
  } while (0)

typedef struct {
  double complex values[64];
  size_t calls;
  size_t trace[64];
  size_t trace_count;
  MVMCKrylovStatus forced_status;
  int return_nonfinite;
  int mismatch_zero_classification;
} TableAmplitude;

static size_t SmallConfigurationIndex(const uint64_t *words,
                                      size_t word_count) {
  if (word_count != 1) return 63;
  return (size_t)(words[0] & UINT64_C(63));
}

static MVMCKrylovStatus TableAmplitudeCallback(
    const uint64_t *words, size_t word_count, void *context,
    MVMCKrylovAmplitudeResult *result) {
  TableAmplitude *table = (TableAmplitude *)context;
  const size_t index = SmallConfigurationIndex(words, word_count);
  ++table->calls;
  if (table->trace_count < sizeof(table->trace) / sizeof(table->trace[0])) {
    table->trace[table->trace_count++] = index;
  }
  if (table->forced_status != MVMC_KRYLOV_STATUS_OK) {
    return table->forced_status;
  }
  memset(result, 0, sizeof(*result));
  result->value = table->return_nonfinite ? NAN + 0.0 * I : table->values[index];
  result->total_zero =
      creal(result->value) == 0.0 && cimag(result->value) == 0.0;
  if (table->mismatch_zero_classification) {
    result->total_zero = !result->total_zero;
  }
  result->regular_component_count = 1;
  result->local_factorization_count = 1;
  result->global_factorization_count = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovLimits DefaultLimits(int order) {
  MVMCKrylovLimits limits;
  limits.max_states = 512;
  limits.max_transitions = 4096;
  limits.max_amplitude_evaluations = 512;
  limits.max_bytes = 16 * 1024 * 1024;
  limits.max_order = order;
  return limits;
}

static int ComplexClose(double complex actual, double complex expected,
                        double tolerance) {
  const double scale = 1.0 + cabs(expected);
  return cabs(actual - expected) <= tolerance * scale;
}

static int DeterministicStatisticsEqual(const MVMCKrylovStatistics *lhs,
                                        const MVMCKrylovStatistics *rhs) {
  return lhs->root_evaluations == rhs->root_evaluations &&
         lhs->requested_order == rhs->requested_order &&
         lhs->completed_order == rhs->completed_order &&
         memcmp(lhs->recursion_calls, rhs->recursion_calls,
                sizeof(lhs->recursion_calls)) == 0 &&
         memcmp(lhs->unique_states, rhs->unique_states,
                sizeof(lhs->unique_states)) == 0 &&
         memcmp(lhs->memo_hits, rhs->memo_hits, sizeof(lhs->memo_hits)) == 0 &&
         memcmp(lhs->memo_misses, rhs->memo_misses,
                sizeof(lhs->memo_misses)) == 0 &&
         lhs->raw_transitions == rhs->raw_transitions &&
         lhs->merged_duplicate_transitions ==
             rhs->merged_duplicate_transitions &&
         lhs->cancelled_zero_transitions ==
             rhs->cancelled_zero_transitions &&
         lhs->terminal_amplitude_requests ==
             rhs->terminal_amplitude_requests &&
         lhs->terminal_cache_hits == rhs->terminal_cache_hits &&
         lhs->regular_component_count == rhs->regular_component_count &&
         lhs->near_pivot_component_count ==
             rhs->near_pivot_component_count &&
         lhs->singular_component_count == rhs->singular_component_count &&
         lhs->total_zero_count == rhs->total_zero_count &&
         lhs->local_factorization_count == rhs->local_factorization_count &&
         lhs->global_factorization_count == rhs->global_factorization_count &&
         lhs->frontier_peak == rhs->frontier_peak &&
         lhs->memo_entries_peak == rhs->memo_entries_peak &&
         lhs->workspace_bytes == rhs->workspace_bytes;
}

static int FiniteNonnegativeTimings(
    const MVMCKrylovStatistics *statistics) {
  int depth;
  if (statistics == NULL ||
      !isfinite(statistics->amplitude_wall_seconds) ||
      statistics->amplitude_wall_seconds < 0.0 ||
      !isfinite(statistics->connectivity_wall_seconds) ||
      statistics->connectivity_wall_seconds < 0.0) {
    return 0;
  }
  for (depth = 0; depth <= MVMC_KRYLOV_MAX_ORDER; ++depth) {
    if (!isfinite(statistics->depth_wall_seconds[depth]) ||
        statistics->depth_wall_seconds[depth] < 0.0) {
      return 0;
    }
  }
  return 1;
}

static MVMCKrylovFockModel ElectronicModel(
    size_t site_count, size_t up_count, size_t down_count,
    const MVMCKrylovHamiltonianTerm *terms, size_t term_count,
    const MVMCKrylovFermionOperator *operators, size_t operator_count) {
  MVMCKrylovFockModel model;
  memset(&model, 0, sizeof(model));
  model.site_count = site_count;
  model.up_electron_count = up_count;
  model.down_electron_count = down_count;
  model.hermitian = 1;
  model.terms = terms;
  model.term_count = term_count;
  model.operators = operators;
  model.operator_count = operator_count;
  return model;
}

static MVMCKrylovStatus Evaluate(
    const MVMCKrylovFockModel *model, uint64_t root,
    const MVMCKrylovLimits *limits, TableAmplitude *amplitude,
    MVMCKrylovResult *result) {
  MVMCKrylovStatus status;
  MVMCKrylovWorkspace *workspace =
      mvmc_krylov_workspace_create(model->site_count, limits, &status);
  if (workspace == NULL) return status;
  status = mvmc_krylov_evaluate(workspace, model, &root, 1,
                                TableAmplitudeCallback, amplitude, result);
  mvmc_krylov_workspace_destroy(workspace);
  return status;
}

static void TestConfigurationContract(void) {
  MVMCKrylovFockModel model = ElectronicModel(33, 1, 1, NULL, 0, NULL, 0);
  uint64_t words[2] = {UINT64_C(1), UINT64_C(2)};
  CHECK(mvmc_krylov_fock_word_count(33) == 2,
        "33-site configuration must use two words");
  CHECK(mvmc_krylov_fock_validate(&model, words, 2) ==
            MVMC_KRYLOV_STATUS_OK,
        "multiword electronic configuration rejected");
  words[1] |= UINT64_C(4);
  CHECK(mvmc_krylov_fock_validate(&model, words, 2) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "nonzero padding bit accepted");
  words[1] &= ~UINT64_C(4);
  CHECK(mvmc_krylov_fock_validate(&model, words, 1) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "wrong word count accepted");

  model = ElectronicModel(2, 1, 1, NULL, 0, NULL, 0);
  model.pure_spin = 1;
  {
    const uint64_t valid = UINT64_C(1) | (UINT64_C(1) << 3);
    const uint64_t double_occupied = UINT64_C(1) | (UINT64_C(1) << 2);
    CHECK(mvmc_krylov_fock_validate(&model, &valid, 1) ==
              MVMC_KRYLOV_STATUS_OK,
          "valid pure-spin configuration rejected");
    CHECK(mvmc_krylov_fock_validate(&model, &double_occupied, 1) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "pure-spin double occupation/empty site accepted");
  }
}

static void TestDenseThreeStateED(void) {
  const double complex a = 0.7 + 0.2 * I;
  const double complex b = -0.3 + 0.4 * I;
  const double diagonal[3] = {0.25, -0.5, 0.75};
  const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 2},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
      {MVMC_KRYLOV_FERMION_CREATE, 2},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
  };
  const MVMCKrylovHamiltonianTerm terms[] = {
      {diagonal[0], 0, 2, 0, 0}, {diagonal[1], 2, 2, 0, 1},
      {diagonal[2], 4, 2, 0, 2}, {a, 6, 2, 1, 0},
      {conj(a), 8, 2, 1, 1},     {b, 10, 2, 1, 2},
      {conj(b), 12, 2, 1, 3},
  };
  const MVMCKrylovFockModel model = ElectronicModel(
      3, 1, 0, terms, sizeof(terms) / sizeof(terms[0]), operators,
      sizeof(operators) / sizeof(operators[0]));
  const double complex matrix[3][3] = {
      {diagonal[0], a, 0.0},
      {conj(a), diagonal[1], b},
      {0.0, conj(b), diagonal[2]},
  };
  double complex expected[4][3];
  TableAmplitude amplitude;
  MVMCKrylovLimits limits = DefaultLimits(3);
  int order;
  int row;
  memset(&amplitude, 0, sizeof(amplitude));
  expected[0][0] = 0.2 - 0.1 * I;
  expected[0][1] = -0.3 + 0.8 * I;
  expected[0][2] = 0.6 + 0.4 * I;
  amplitude.values[1] = expected[0][0];
  amplitude.values[2] = expected[0][1];
  amplitude.values[4] = expected[0][2];
  for (order = 1; order <= 3; ++order) {
    for (row = 0; row < 3; ++row) {
      int column;
      expected[order][row] = 0.0;
      for (column = 0; column < 3; ++column) {
        expected[order][row] +=
            matrix[row][column] * expected[order - 1][column];
      }
    }
  }
  for (row = 0; row < 3; ++row) {
    MVMCKrylovResult result;
    MVMCKrylovStatus status =
        Evaluate(&model, UINT64_C(1) << row, &limits, &amplitude, &result);
    CHECK(status == MVMC_KRYLOV_STATUS_OK && result.valid,
          "dense ED root %d failed: %s", row,
          mvmc_krylov_status_string(status));
    for (order = 0; order <= 3; ++order) {
      CHECK(ComplexClose(result.value[order], expected[order][row], 2e-13),
            "dense ED root=%d order=%d got=(%.17g,%.17g) expected=(%.17g,%.17g)",
            row, order, creal(result.value[order]), cimag(result.value[order]),
            creal(expected[order][row]), cimag(expected[order][row]));
    }
  }
}

static void TestZeroBridgesAndMemoization(void) {
  const double complex a = 1.0 + 2.0 * I;
  const MVMCKrylovFermionOperator two_site_operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
  };
  const MVMCKrylovHamiltonianTerm two_site_terms[] = {
      {a, 0, 2, 1, 0}, {conj(a), 2, 2, 1, 1}};
  const MVMCKrylovFockModel two_site_model = ElectronicModel(
      2, 1, 0, two_site_terms, 2, two_site_operators, 4);
  TableAmplitude amplitude;
  MVMCKrylovLimits limits = DefaultLimits(3);
  MVMCKrylovResult result;
  const double complex terminal = 3.0 - I;
  memset(&amplitude, 0, sizeof(amplitude));
  amplitude.values[1] = 0.0;
  amplitude.values[2] = terminal;
  CHECK(Evaluate(&two_site_model, UINT64_C(1), &limits, &amplitude, &result) ==
            MVMC_KRYLOV_STATUS_OK,
        "one-step zero bridge evaluation failed");
  CHECK(result.value[0] == 0.0,
        "zero bridge v0 must remain exact zero");
  CHECK(ComplexClose(result.value[1], a * terminal, 1e-14),
        "one-step zero bridge was truncated");
  CHECK(result.value[2] == 0.0,
        "two-site zero bridge v2 must return to zero amplitude");
  CHECK(ComplexClose(result.value[3], a * terminal * (creal(a * conj(a))),
                     1e-14),
        "one-step zero bridge v3 mismatch");
  CHECK(result.statistics.terminal_amplitude_requests == 2,
        "terminal amplitude memoization count=%llu expected=2",
        (unsigned long long)result.statistics.terminal_amplitude_requests);
  CHECK(result.statistics.memo_hits[0] > 0,
        "depth-0 memo was not reused");
  CHECK(FiniteNonnegativeTimings(&result.statistics),
        "profile timings must be finite and nonnegative");

  {
    const MVMCKrylovFermionOperator chain_operators[] = {
        {MVMC_KRYLOV_FERMION_CREATE, 0},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
        {MVMC_KRYLOV_FERMION_CREATE, 1},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
        {MVMC_KRYLOV_FERMION_CREATE, 1},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
        {MVMC_KRYLOV_FERMION_CREATE, 2},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
    };
    const MVMCKrylovHamiltonianTerm chain_terms[] = {
        {2.0, 0, 2, 1, 0}, {2.0, 2, 2, 1, 1},
        {-0.5, 4, 2, 1, 2}, {-0.5, 6, 2, 1, 3},
    };
    const MVMCKrylovFockModel chain_model = ElectronicModel(
        3, 1, 0, chain_terms, 4, chain_operators, 8);
    memset(&amplitude, 0, sizeof(amplitude));
    amplitude.values[4] = 5.0;
    limits.max_order = 2;
    CHECK(Evaluate(&chain_model, UINT64_C(1), &limits, &amplitude, &result) ==
              MVMC_KRYLOV_STATUS_OK,
          "two-step bridge evaluation failed");
    CHECK(result.value[0] == 0.0 && result.value[1] == 0.0,
          "two-step bridge lower orders must be zero");
    CHECK(ComplexClose(result.value[2], -5.0, 1e-14),
          "two-step bridge v2 mismatch: %.17g%+.17gi",
          creal(result.value[2]), cimag(result.value[2]));
  }
}

static void TestFermionSignAndExchange(void) {
  const MVMCKrylovFermionOperator hop_operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 2},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
  };
  const MVMCKrylovHamiltonianTerm hop_terms[] = {
      {1.0, 0, 2, 1, 0}, {1.0, 2, 2, 1, 1}};
  const MVMCKrylovFockModel hop_model =
      ElectronicModel(3, 2, 0, hop_terms, 2, hop_operators, 4);
  TableAmplitude amplitude;
  MVMCKrylovLimits limits = DefaultLimits(1);
  MVMCKrylovResult result;
  memset(&amplitude, 0, sizeof(amplitude));
  amplitude.values[UINT64_C(6)] = 1.0;
  CHECK(Evaluate(&hop_model, UINT64_C(3), &limits, &amplitude, &result) ==
            MVMC_KRYLOV_STATUS_OK,
        "fermion sign evaluation failed");
  CHECK(ComplexClose(result.value[1], -1.0, 1e-14),
        "same-spin occupied-orbital parity sign mismatch");

  {
    const MVMCKrylovFermionOperator down_operators[] = {
        {MVMC_KRYLOV_FERMION_CREATE, 5},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 3},
        {MVMC_KRYLOV_FERMION_CREATE, 3},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 5},
    };
    const MVMCKrylovHamiltonianTerm down_terms[] = {
        {1.0, 0, 2, 1, 0}, {1.0, 2, 2, 1, 1}};
    const MVMCKrylovFockModel down_model =
        ElectronicModel(3, 1, 2, down_terms, 2, down_operators, 4);
    memset(&amplitude, 0, sizeof(amplitude));
    amplitude.values[UINT64_C(49)] = 1.0;
    CHECK(Evaluate(&down_model, UINT64_C(25), &limits, &amplitude,
                   &result) == MVMC_KRYLOV_STATUS_OK,
          "down-block boundary parity evaluation failed");
    CHECK(ComplexClose(result.value[1], -1.0, 1e-14),
          "up/down block boundary parity sign mismatch");
  }

  {
    const MVMCKrylovFermionOperator exchange_operators[] = {
        {MVMC_KRYLOV_FERMION_CREATE, 0},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
        {MVMC_KRYLOV_FERMION_CREATE, 3},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
        {MVMC_KRYLOV_FERMION_CREATE, 1},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
        {MVMC_KRYLOV_FERMION_CREATE, 2},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 3},
    };
    const MVMCKrylovHamiltonianTerm exchange_terms[] = {
        {0.75, 0, 4, 2, 0}, {0.75, 4, 4, 2, 1}};
    MVMCKrylovFockModel exchange_model = ElectronicModel(
        2, 1, 1, exchange_terms, 2, exchange_operators, 8);
    exchange_model.pure_spin = 1;
    memset(&amplitude, 0, sizeof(amplitude));
    amplitude.values[UINT64_C(9)] = 2.0;
    CHECK(Evaluate(&exchange_model, UINT64_C(6), &limits, &amplitude,
                   &result) == MVMC_KRYLOV_STATUS_OK,
          "Exchange orientation evaluation failed");
    CHECK(ComplexClose(result.value[1], 1.5, 1e-14),
          "Exchange orientation/sign mismatch");
    memset(&amplitude, 0, sizeof(amplitude));
    amplitude.values[UINT64_C(6)] = 2.0;
    CHECK(Evaluate(&exchange_model, UINT64_C(9), &limits, &amplitude,
                   &result) == MVMC_KRYLOV_STATUS_OK,
          "reverse Exchange orientation evaluation failed");
    CHECK(ComplexClose(result.value[1], 1.5, 1e-14),
          "reverse Exchange orientation/sign mismatch");
  }
}

static void TestDuplicateMergeAndPermutation(void) {
  const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
  };
  const MVMCKrylovHamiltonianTerm terms_a[] = {
      {1.0e16, 0, 2, 1, 10}, {1.0, 0, 2, 1, 12},
      {-1.0e16, 0, 2, 1, 11}, {1.0e16, 2, 2, 1, 10},
      {1.0, 2, 2, 1, 12}, {-1.0e16, 2, 2, 1, 11},
  };
  const MVMCKrylovHamiltonianTerm terms_b[] = {
      {-1.0e16, 2, 2, 1, 11}, {1.0, 0, 2, 1, 12},
      {1.0e16, 2, 2, 1, 10},  {-1.0e16, 0, 2, 1, 11},
      {1.0, 2, 2, 1, 12},     {1.0e16, 0, 2, 1, 10},
  };
  const MVMCKrylovFockModel model_a =
      ElectronicModel(2, 1, 0, terms_a, 6, operators, 4);
  const MVMCKrylovFockModel model_b =
      ElectronicModel(2, 1, 0, terms_b, 6, operators, 4);
  TableAmplitude amplitude_a;
  TableAmplitude amplitude_b;
  MVMCKrylovLimits limits = DefaultLimits(3);
  MVMCKrylovResult result_a;
  MVMCKrylovResult result_b;
  memset(&amplitude_a, 0, sizeof(amplitude_a));
  memset(&amplitude_b, 0, sizeof(amplitude_b));
  amplitude_a.values[1] = amplitude_b.values[1] = 3.0;
  amplitude_a.values[2] = amplitude_b.values[2] = -2.0;
  CHECK(Evaluate(&model_a, UINT64_C(1), &limits, &amplitude_a, &result_a) ==
            MVMC_KRYLOV_STATUS_OK,
        "duplicate model A failed");
  CHECK(Evaluate(&model_b, UINT64_C(1), &limits, &amplitude_b, &result_b) ==
            MVMC_KRYLOV_STATUS_OK,
        "duplicate model B failed");
  CHECK(memcmp(result_a.value, result_b.value, sizeof(result_a.value)) == 0,
        "term permutation changed bitwise Krylov values");
  CHECK(amplitude_a.trace_count == amplitude_b.trace_count &&
            memcmp(amplitude_a.trace, amplitude_b.trace,
                   amplitude_a.trace_count * sizeof(amplitude_a.trace[0])) ==
                0,
        "term permutation changed terminal callback order");
  CHECK(DeterministicStatisticsEqual(&result_a.statistics,
                                     &result_b.statistics),
        "term permutation changed deterministic statistics");
  CHECK(ComplexClose(result_a.value[1], -2.0, 1e-14),
        "compensated duplicate merge lost residual coefficient");
  CHECK(result_a.statistics.merged_duplicate_transitions > 0,
        "duplicate transition counter was not incremented");

  {
    const double complex first = 0.75 + 1.25 * I;
    const double complex second = -0.5 + 0.125 * I;
    const MVMCKrylovHamiltonianTerm complex_terms[] = {
        {first, 0, 2, 1, 0},       {second, 0, 2, 1, 1},
        {conj(first), 2, 2, 1, 0}, {conj(second), 2, 2, 1, 1},
    };
    const MVMCKrylovFockModel complex_model =
        ElectronicModel(2, 1, 0, complex_terms, 4, operators, 4);
    TableAmplitude complex_amplitude;
    MVMCKrylovResult complex_result;
    memset(&complex_amplitude, 0, sizeof(complex_amplitude));
    complex_amplitude.values[2] = 2.0 - I;
    CHECK(Evaluate(&complex_model, UINT64_C(1), &limits,
                   &complex_amplitude, &complex_result) ==
              MVMC_KRYLOV_STATUS_OK,
          "complex constructive duplicate model failed");
    CHECK(ComplexClose(complex_result.value[1],
                       (first + second) * (2.0 - I), 1e-14),
          "complex constructive duplicate merge mismatch");
  }

  {
    const double complex cancel = 1.0 + 2.0 * I;
    const MVMCKrylovHamiltonianTerm cancel_terms[] = {
        {cancel, 0, 2, 1, 0},       {-cancel, 0, 2, 1, 1},
        {conj(cancel), 2, 2, 1, 0}, {-conj(cancel), 2, 2, 1, 1},
    };
    const MVMCKrylovFockModel cancel_model =
        ElectronicModel(2, 1, 0, cancel_terms, 4, operators, 4);
    TableAmplitude cancel_amplitude;
    MVMCKrylovResult cancel_result;
    memset(&cancel_amplitude, 0, sizeof(cancel_amplitude));
    cancel_amplitude.values[1] = 1.0;
    cancel_amplitude.values[2] = 9.0;
    CHECK(Evaluate(&cancel_model, UINT64_C(1), &limits, &cancel_amplitude,
                   &cancel_result) == MVMC_KRYLOV_STATUS_OK,
          "destructive duplicate model failed");
    CHECK(cancel_result.value[1] == 0.0 && cancel_result.value[2] == 0.0 &&
              cancel_result.value[3] == 0.0,
          "exactly cancelled connection was retained");
    CHECK(cancel_result.statistics.cancelled_zero_transitions > 0,
          "cancelled transition counter was not incremented");
    CHECK(cancel_result.statistics.terminal_amplitude_requests == 1,
          "cancelled neighbor still evaluated terminal amplitude");
  }
}

static void CheckAtomicFailure(const MVMCKrylovResult *result,
                               MVMCKrylovStatus expected,
                               const char *label) {
  int order;
  CHECK(!result->valid && result->status == expected,
        "%s did not return atomic invalid result", label);
  CHECK(result->evaluated_order == -1,
        "%s exposed a partial evaluated order", label);
  for (order = 0; order <= MVMC_KRYLOV_MAX_ORDER; ++order) {
    CHECK(result->value[order] == 0.0,
          "%s exposed partial value[%d]", label, order);
  }
}

static void TestFailuresAndLimits(void) {
  const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
  };
  const MVMCKrylovHamiltonianTerm terms[] = {
      {1.0, 0, 2, 1, 0}, {1.0, 2, 2, 1, 1}};
  MVMCKrylovFockModel model =
      ElectronicModel(2, 1, 0, terms, 2, operators, 4);
  TableAmplitude amplitude;
  MVMCKrylovLimits limits = DefaultLimits(2);
  MVMCKrylovResult result;
  MVMCKrylovStatus status;
  MVMCKrylovWorkspace *workspace;
  memset(&amplitude, 0, sizeof(amplitude));
  amplitude.values[1] = 1.0;
  amplitude.values[2] = 2.0;

  model.hermitian = 0;
  status = Evaluate(&model, UINT64_C(1), &limits, &amplitude, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_NON_HERMITIAN_MODEL,
        "non-Hermitian model status=%d", (int)status);
  model.hermitian = 1;

  limits.max_states = 1;
  status = Evaluate(&model, UINT64_C(1), &limits, &amplitude, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
        "state limit did not fire");
  CheckAtomicFailure(&result, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                     "state limit");

  limits = DefaultLimits(2);
  limits.max_transitions = 1;
  status = Evaluate(&model, UINT64_C(1), &limits, &amplitude, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
        "transition limit did not fire");
  CheckAtomicFailure(&result, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                     "transition limit");

  limits = DefaultLimits(1);
  limits.max_amplitude_evaluations = 1;
  status = Evaluate(&model, UINT64_C(1), &limits, &amplitude, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
        "amplitude limit did not fire");
  CheckAtomicFailure(&result, MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
                     "amplitude limit");

  limits = DefaultLimits(0);
  limits.max_bytes = 1;
  workspace = mvmc_krylov_workspace_create(2, &limits, &status);
  CHECK(workspace == NULL && status == MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
        "byte limit did not fail workspace creation");
  mvmc_krylov_workspace_destroy(workspace);

  limits = DefaultLimits(0);
  limits.max_states = 0;
  workspace = mvmc_krylov_workspace_create(2, &limits, &status);
  CHECK(workspace == NULL && status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "zero state limit was not rejected");

  limits = DefaultLimits(0);
  limits.max_transitions = 0;
  workspace = mvmc_krylov_workspace_create(2, &limits, &status);
  CHECK(workspace == NULL && status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "zero transition limit was not rejected");

  limits = DefaultLimits(0);
  limits.max_amplitude_evaluations = 0;
  workspace = mvmc_krylov_workspace_create(2, &limits, &status);
  CHECK(workspace == NULL && status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "zero amplitude-evaluation limit was not rejected");

  limits = DefaultLimits(0);
  limits.max_bytes = 0;
  workspace = mvmc_krylov_workspace_create(2, &limits, &status);
  CHECK(workspace == NULL && status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "zero byte limit was not rejected");

  limits = DefaultLimits(-1);
  workspace = mvmc_krylov_workspace_create(2, &limits, &status);
  CHECK(workspace == NULL && status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "negative max order was not rejected");

  limits = DefaultLimits(MVMC_KRYLOV_MAX_ORDER + 1);
  workspace = mvmc_krylov_workspace_create(2, &limits, &status);
  CHECK(workspace == NULL && status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "max order above supported range was not rejected");

  limits = DefaultLimits(0);
  limits.max_states = SIZE_MAX;
  workspace = mvmc_krylov_workspace_create(2, &limits, &status);
  CHECK(workspace == NULL && status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "overflowing state limit was not rejected");

  limits = DefaultLimits(0);
  limits.max_transitions = SIZE_MAX;
  workspace = mvmc_krylov_workspace_create(33, &limits, &status);
  CHECK(workspace == NULL && status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "overflowing transition-word allocation was not rejected");

  workspace = mvmc_krylov_workspace_create(SIZE_MAX, &limits, &status);
  CHECK(workspace == NULL && status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "overflowing site count was not rejected");

  limits = DefaultLimits(0);
  memset(&amplitude, 0, sizeof(amplitude));
  amplitude.return_nonfinite = 1;
  status = Evaluate(&model, UINT64_C(1), &limits, &amplitude, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_NONFINITE,
        "nonfinite callback value status=%d", (int)status);
  CheckAtomicFailure(&result, MVMC_KRYLOV_STATUS_NONFINITE,
                     "nonfinite amplitude");

  memset(&amplitude, 0, sizeof(amplitude));
  amplitude.forced_status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  status = Evaluate(&model, UINT64_C(1), &limits, &amplitude, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE,
        "callback failure status=%d", (int)status);
  CheckAtomicFailure(&result, MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE,
                     "callback failure");

  memset(&amplitude, 0, sizeof(amplitude));
  amplitude.values[1] = 1.0;
  amplitude.mismatch_zero_classification = 1;
  status = Evaluate(&model, UINT64_C(1), &limits, &amplitude, &result);
  CHECK(status == MVMC_KRYLOV_STATUS_NONFINITE,
        "inconsistent zero classification accepted");

  {
    MVMCKrylovFermionOperator invalid_operators[4];
    memcpy(invalid_operators, operators, sizeof(invalid_operators));
    invalid_operators[0].orbital = 4;
    model.operators = invalid_operators;
    memset(&amplitude, 0, sizeof(amplitude));
    status = Evaluate(&model, UINT64_C(1), &limits, &amplitude, &result);
    CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "out-of-range orbital accepted");
    model.operators = operators;
  }

  {
    MVMCKrylovFermionOperator invalid_operators[4];
    memcpy(invalid_operators, operators, sizeof(invalid_operators));
    invalid_operators[0].kind = (MVMCKrylovFermionOperatorKind)99;
    model.operators = invalid_operators;
    status = Evaluate(&model, UINT64_C(1), &limits, &amplitude, &result);
    CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "illegal operator kind accepted");
    model.operators = operators;
  }

  {
    MVMCKrylovHamiltonianTerm invalid_terms[2];
    memcpy(invalid_terms, terms, sizeof(invalid_terms));
    invalid_terms[0].operator_count = 3;
    model.terms = invalid_terms;
    status = Evaluate(&model, UINT64_C(1), &limits, &amplitude, &result);
    CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "odd fermion operator count accepted");
    CheckAtomicFailure(&result, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                       "odd operator count");
    model.terms = terms;
  }

  {
    MVMCKrylovHamiltonianTerm invalid_terms[2];
    memcpy(invalid_terms, terms, sizeof(invalid_terms));
    invalid_terms[0].coefficient = NAN + 0.0 * I;
    model.terms = invalid_terms;
    status = Evaluate(&model, UINT64_C(1), &limits, &amplitude, &result);
    CHECK(status == MVMC_KRYLOV_STATUS_NONFINITE,
          "nonfinite model coefficient status=%d", (int)status);
    CheckAtomicFailure(&result, MVMC_KRYLOV_STATUS_NONFINITE,
                       "nonfinite model coefficient");
    model.terms = terms;
  }

  {
    const MVMCKrylovHamiltonianTerm overflow_terms[] = {
        {DBL_MAX, 0, 2, 1, 0}, {DBL_MAX, 2, 2, 1, 1}};
    const MVMCKrylovFockModel overflow_model =
        ElectronicModel(2, 1, 0, overflow_terms, 2, operators, 4);
    limits = DefaultLimits(1);
    memset(&amplitude, 0, sizeof(amplitude));
    amplitude.values[2] = 2.0;
    status = Evaluate(&overflow_model, UINT64_C(1), &limits, &amplitude,
                      &result);
    CHECK(status == MVMC_KRYLOV_STATUS_NONFINITE,
          "nonfinite coefficient-amplitude product status=%d", (int)status);
    CheckAtomicFailure(&result, MVMC_KRYLOV_STATUS_NONFINITE,
                       "nonfinite product");
  }

  {
    const MVMCKrylovHamiltonianTerm sum_overflow_terms[] = {
        {DBL_MAX, 0, 2, 1, 0}, {DBL_MAX, 0, 2, 1, 1},
        {DBL_MAX, 2, 2, 1, 0}, {DBL_MAX, 2, 2, 1, 1},
    };
    const MVMCKrylovFockModel sum_overflow_model = ElectronicModel(
        2, 1, 0, sum_overflow_terms, 4, operators, 4);
    limits = DefaultLimits(1);
    memset(&amplitude, 0, sizeof(amplitude));
    amplitude.values[2] = 1.0;
    status = Evaluate(&sum_overflow_model, UINT64_C(1), &limits, &amplitude,
                      &result);
    CHECK(status == MVMC_KRYLOV_STATUS_NONFINITE,
          "nonfinite duplicate sum status=%d", (int)status);
    CheckAtomicFailure(&result, MVMC_KRYLOV_STATUS_NONFINITE,
                       "nonfinite duplicate sum");
  }

  {
    const MVMCKrylovFermionOperator spin_flip_operators[] = {
        {MVMC_KRYLOV_FERMION_CREATE, 2},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
        {MVMC_KRYLOV_FERMION_CREATE, 0},
        {MVMC_KRYLOV_FERMION_ANNIHILATE, 2},
    };
    const MVMCKrylovHamiltonianTerm spin_flip_terms[] = {
        {1.0, 0, 2, 1, 0}, {1.0, 2, 2, 1, 1}};
    const MVMCKrylovFockModel spin_flip_model = ElectronicModel(
        2, 1, 0, spin_flip_terms, 2, spin_flip_operators, 4);
    limits = DefaultLimits(1);
    memset(&amplitude, 0, sizeof(amplitude));
    amplitude.values[1] = 1.0;
    status = Evaluate(&spin_flip_model, UINT64_C(1), &limits, &amplitude,
                      &result);
    CHECK(status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "sector-changing terminal status=%d", (int)status);
    CheckAtomicFailure(&result, MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
                       "sector-changing terminal");
  }

  {
    const MVMCKrylovFockModel filled_model =
        ElectronicModel(2, 2, 0, terms, 2, operators, 4);
    limits = DefaultLimits(1);
    memset(&amplitude, 0, sizeof(amplitude));
    amplitude.values[3] = 2.0;
    status = Evaluate(&filled_model, UINT64_C(3), &limits, &amplitude,
                      &result);
    CHECK(status == MVMC_KRYLOV_STATUS_OK && result.value[1] == 0.0,
          "Pauli-blocked hopping was not a normal zero transition");
    CHECK(result.statistics.raw_transitions == 0,
          "Pauli-blocked term incremented transition counter");
  }
}

static void TestRepeatedEvaluationAndOrderDomain(void) {
  const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
  };
  const MVMCKrylovHamiltonianTerm terms[] = {
      {1.0, 0, 2, 1, 0}, {1.0, 2, 2, 1, 1}};
  const MVMCKrylovFockModel model =
      ElectronicModel(2, 1, 0, terms, 2, operators, 4);
  int order;
  for (order = 0; order <= 3; ++order) {
    MVMCKrylovLimits limits = DefaultLimits(order);
    MVMCKrylovStatus status;
    MVMCKrylovWorkspace *workspace =
        mvmc_krylov_workspace_create(2, &limits, &status);
    TableAmplitude amplitude;
    MVMCKrylovResult first;
    MVMCKrylovResult second;
    const uint64_t root = UINT64_C(1);
    memset(&amplitude, 0, sizeof(amplitude));
    amplitude.values[1] = 2.0;
    amplitude.values[2] = 5.0;
    CHECK(workspace != NULL && status == MVMC_KRYLOV_STATUS_OK,
          "workspace create failed for order %d", order);
    if (workspace == NULL) continue;
    CHECK(mvmc_krylov_evaluate(workspace, &model, &root, 1,
                               TableAmplitudeCallback, &amplitude,
                               &first) == MVMC_KRYLOV_STATUS_OK,
          "first evaluation failed for order %d", order);
    CHECK(mvmc_krylov_evaluate(workspace, &model, &root, 1,
                               TableAmplitudeCallback, &amplitude,
                               &second) == MVMC_KRYLOV_STATUS_OK,
          "second evaluation failed for order %d", order);
    CHECK(memcmp(first.value, second.value, sizeof(first.value)) == 0,
          "repeated evaluation retained stale cache for order %d", order);
    CHECK(first.evaluated_order == order && second.evaluated_order == order,
          "evaluated order mismatch for %d", order);
    mvmc_krylov_workspace_destroy(workspace);
  }
}

int main(void) {
  TestConfigurationContract();
  TestDenseThreeStateED();
  TestZeroBridgesAndMemoization();
  TestFermionSignAndExchange();
  TestDuplicateMergeAndPermutation();
  TestFailuresAndLimits();
  TestRepeatedEvaluationAndOrderDomain();

  if (failures != 0) {
    fprintf(stderr, "AbsoluteKrylovReference_Unit: %d failure(s)\n", failures);
    return 1;
  }
  printf("AbsoluteKrylovReference_Unit: PASS\n");
  return 0;
}
