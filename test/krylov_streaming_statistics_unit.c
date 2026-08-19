#include "krylov_streaming_statistics.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(condition, ...)                                                 \
  do {                                                                        \
    if (!(condition)) {                                                       \
      fprintf(stderr, "KrylovStreamingStatistics_Unit FAIL: ");              \
      fprintf(stderr, __VA_ARGS__);                                           \
      fprintf(stderr, " (line %d)\n", __LINE__);                            \
      ++failures;                                                             \
    }                                                                         \
  } while (0)

static int close_double(double actual, double expected, double tolerance) {
  return fabs(actual - expected) <= tolerance * fmax(1.0, fabs(expected));
}

static int close_complex(
    double complex actual, double complex expected, double tolerance) {
  return close_double(creal(actual), creal(expected), tolerance) &&
         close_double(cimag(actual), cimag(expected), tolerance);
}

static void test_upper_triangle_accumulator(void) {
  MVMCKrylovStreamingComplexSum direct_entries[6];
  MVMCKrylovStreamingComplexSum left_entries[6];
  MVMCKrylovStreamingComplexSum right_entries[6];
  MVMCKrylovStreamingMatrixAccumulator direct;
  MVMCKrylovStreamingMatrixAccumulator left;
  MVMCKrylovStreamingMatrixAccumulator right;
  MVMCKrylovStreamingMatrixSummary summary;
  const double complex sample_a[6] = {
      1.0 + 0.5 * I, 2.0 + 1.0 * I, 3.0 + 1.5 * I,
      4.0 + 2.0 * I, 5.0 + 2.5 * I, 6.0 + 3.0 * I};
  const double complex sample_b[6] = {
      -2.0 + 1.0 * I, -1.0 + 2.0 * I, 0.0 + 3.0 * I,
      1.0 + 4.0 * I, 2.0 + 5.0 * I, 3.0 + 6.0 * I};
  double complex direct_mean[6];
  double complex merged_mean[6];
  size_t upper_count = 0;
  size_t index = 999;
  size_t entry;

  CHECK(mvmc_krylov_streaming_upper_count(3, &upper_count) ==
            MVMC_KRYLOV_STATUS_OK &&
            upper_count == 6,
        "upper-triangle count for dim=3");
  CHECK(mvmc_krylov_streaming_upper_index(3, 0, 2, &index) ==
            MVMC_KRYLOV_STATUS_OK &&
            index == 2,
        "upper index row-major upper");
  CHECK(mvmc_krylov_streaming_upper_index(3, 2, 1, &index) ==
            MVMC_KRYLOV_STATUS_OK &&
            index == 4,
        "upper index canonicalizes lower input");
  CHECK(mvmc_krylov_streaming_upper_count(0, &upper_count) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "zero dimension rejected");

  CHECK(mvmc_krylov_streaming_matrix_accumulator_init(
            3, direct_entries, 6, &direct) == MVMC_KRYLOV_STATUS_OK &&
            direct.valid,
        "direct matrix accumulator init");
  CHECK(mvmc_krylov_streaming_matrix_accumulator_init(
            3, left_entries, 6, &left) == MVMC_KRYLOV_STATUS_OK &&
            left.valid,
        "left matrix accumulator init");
  CHECK(mvmc_krylov_streaming_matrix_accumulator_init(
            3, right_entries, 6, &right) == MVMC_KRYLOV_STATUS_OK &&
            right.valid,
        "right matrix accumulator init");

  CHECK(mvmc_krylov_streaming_matrix_accumulator_add_sample(
            &direct, sample_a, 6, 1.0) == MVMC_KRYLOV_STATUS_OK,
        "direct sample A");
  CHECK(mvmc_krylov_streaming_matrix_accumulator_add_sample(
            &direct, sample_b, 6, 1.0) == MVMC_KRYLOV_STATUS_OK,
        "direct sample B");
  CHECK(mvmc_krylov_streaming_matrix_accumulator_add_sample(
            &left, sample_a, 6, 1.0) == MVMC_KRYLOV_STATUS_OK,
        "left sample A");
  CHECK(mvmc_krylov_streaming_matrix_accumulator_add_sample(
            &right, sample_b, 6, 1.0) == MVMC_KRYLOV_STATUS_OK,
        "right sample B");
  CHECK(mvmc_krylov_streaming_matrix_accumulator_merge(&left, &right) ==
            MVMC_KRYLOV_STATUS_OK &&
            left.sample_count == 2,
        "deterministic partition merge");
  CHECK(mvmc_krylov_streaming_matrix_accumulator_mean(
            &direct, direct_mean, 6) == MVMC_KRYLOV_STATUS_OK,
        "direct matrix mean");
  CHECK(mvmc_krylov_streaming_matrix_accumulator_mean(
            &left, merged_mean, 6) == MVMC_KRYLOV_STATUS_OK,
        "merged matrix mean");
  for (entry = 0; entry < 6; ++entry) {
    CHECK(close_complex(merged_mean[entry], direct_mean[entry], 1.0e-15),
          "merged mean entry %zu", entry);
  }
  CHECK(mvmc_krylov_streaming_matrix_accumulator_summary(
            &left, &summary) == MVMC_KRYLOV_STATUS_OK &&
            summary.valid && summary.dimension == 3 &&
            summary.upper_count == 6 && summary.sample_count == 2 &&
            summary.storage_bytes ==
                6 * sizeof(MVMCKrylovStreamingComplexSum),
        "matrix accumulator storage summary");
}

typedef struct {
  double complex estimate;
  double denominator_relative_se;
  double variance_real;
  double variance_imag;
  double covariance_real_imag;
} RatioOracle;

static RatioOracle ratio_oracle(
    const double complex *numerator, const double *denominator,
    size_t sample_count) {
  RatioOracle oracle;
  double complex numerator_sum = 0.0;
  double denominator_sum = 0.0;
  double denominator_mean;
  double denominator_ss = 0.0;
  double residual_real_ss = 0.0;
  double residual_imag_ss = 0.0;
  double residual_real_imag_ss = 0.0;
  const double count = (double)sample_count;
  size_t sample;

  memset(&oracle, 0, sizeof(oracle));
  for (sample = 0; sample < sample_count; ++sample) {
    numerator_sum += numerator[sample];
    denominator_sum += denominator[sample];
  }
  denominator_mean = denominator_sum / count;
  oracle.estimate = numerator_sum / denominator_sum;
  for (sample = 0; sample < sample_count; ++sample) {
    const double denominator_delta = denominator[sample] - denominator_mean;
    const double residual_real =
        creal(numerator[sample]) - creal(oracle.estimate) *
                                     denominator[sample];
    const double residual_imag =
        cimag(numerator[sample]) - cimag(oracle.estimate) *
                                     denominator[sample];
    denominator_ss += denominator_delta * denominator_delta;
    residual_real_ss += residual_real * residual_real;
    residual_imag_ss += residual_imag * residual_imag;
    residual_real_imag_ss += residual_real * residual_imag;
  }
  oracle.denominator_relative_se =
      sqrt(denominator_ss / ((count - 1.0) * count)) /
      fabs(denominator_mean);
  oracle.variance_real =
      residual_real_ss / ((count - 1.0) * count *
                          denominator_mean * denominator_mean);
  oracle.variance_imag =
      residual_imag_ss / ((count - 1.0) * count *
                          denominator_mean * denominator_mean);
  oracle.covariance_real_imag =
      residual_real_imag_ss / ((count - 1.0) * count *
                               denominator_mean * denominator_mean);
  return oracle;
}

static void test_self_normalized_ratio(void) {
  const double denominator[4] = {2.0, 4.0, 5.0, 3.0};
  const double complex numerator[4] = {
      4.0 + 1.0 * I, 9.0 + 0.0 * I, 10.0 + 5.0 * I,
      7.0 + 3.0 * I};
  const RatioOracle oracle = ratio_oracle(numerator, denominator, 4);
  MVMCKrylovStreamingRatioAccumulator accumulator;
  MVMCKrylovStreamingRatioResult result;
  size_t sample;

  CHECK(mvmc_krylov_streaming_ratio_accumulator_reset(&accumulator) ==
            MVMC_KRYLOV_STATUS_OK &&
            accumulator.valid,
        "ratio accumulator reset");
  for (sample = 0; sample < 4; ++sample) {
    CHECK(mvmc_krylov_streaming_ratio_accumulator_add_sample(
              &accumulator, numerator[sample], denominator[sample]) ==
              MVMC_KRYLOV_STATUS_OK,
          "ratio sample %zu", sample);
  }
  CHECK(mvmc_krylov_streaming_ratio_accumulator_finalize(
            &accumulator, 1.0e-12, 1.0, &result) ==
            MVMC_KRYLOV_STATUS_OK &&
            result.valid && result.denominator_stable,
        "ratio finalize");
  CHECK(close_complex(result.estimate, oracle.estimate, 1.0e-14),
        "self-normalized estimate");
  CHECK(close_double(result.denominator_relative_se,
                     oracle.denominator_relative_se, 1.0e-14),
        "denominator relative SE");
  CHECK(close_double(result.variance_real, oracle.variance_real, 1.0e-14),
        "ratio real variance");
  CHECK(close_double(result.variance_imag, oracle.variance_imag, 1.0e-14),
        "ratio imag variance");
  CHECK(close_double(result.covariance_real_imag,
                     oracle.covariance_real_imag, 1.0e-14),
        "ratio real/imag covariance");
  CHECK(mvmc_krylov_streaming_ratio_accumulator_finalize(
            &accumulator, 1.0e-12, 1.0e-6, &result) ==
            MVMC_KRYLOV_STATUS_OK &&
            result.valid && !result.denominator_stable,
        "denominator stability gate failure");
  CHECK(mvmc_krylov_streaming_ratio_accumulator_add_sample(
            &accumulator, 1.0, 0.0) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "non-positive denominator rejected");
}

static MVMCKrylovJackknifeResult jackknife_oracle(
    const MVMCKrylovJackknifeBlock *blocks, size_t block_count) {
  MVMCKrylovJackknifeResult result;
  double complex total_numerator = 0.0;
  double denominator_total = 0.0;
  double complex leave_one_sum = 0.0;
  double complex mean;
  double block_count_double = (double)block_count;
  size_t block;

  memset(&result, 0, sizeof(result));
  for (block = 0; block < block_count; ++block) {
    total_numerator += blocks[block].numerator;
    denominator_total += blocks[block].denominator;
    result.sample_count += blocks[block].sample_count;
  }
  result.estimate = total_numerator / denominator_total;
  for (block = 0; block < block_count; ++block) {
    leave_one_sum +=
        (total_numerator - blocks[block].numerator) /
        (denominator_total - blocks[block].denominator);
  }
  mean = leave_one_sum / block_count_double;
  result.jackknife_mean = mean;
  for (block = 0; block < block_count; ++block) {
    const double complex leave =
        (total_numerator - blocks[block].numerator) /
        (denominator_total - blocks[block].denominator);
    const double complex delta = leave - mean;
    result.variance_real += creal(delta) * creal(delta);
    result.variance_imag += cimag(delta) * cimag(delta);
    result.covariance_real_imag += creal(delta) * cimag(delta);
  }
  result.variance_real *= (block_count_double - 1.0) / block_count_double;
  result.variance_imag *= (block_count_double - 1.0) / block_count_double;
  result.covariance_real_imag *=
      (block_count_double - 1.0) / block_count_double;
  result.complex_se = sqrt(result.variance_real + result.variance_imag);
  return result;
}

static void fill_jackknife_blocks(
    MVMCKrylovJackknifeBlock *blocks, size_t block_count,
    uint64_t samples_per_block, double denominator_base,
    double denominator_step, double complex center) {
  size_t block;
  const double midpoint = 0.5 * (double)(block_count - 1);
  for (block = 0; block < block_count; ++block) {
    const double offset = (double)block - midpoint;
    const double denominator =
        denominator_base + denominator_step * (double)(block % 5);
    blocks[block].denominator = denominator;
    blocks[block].numerator =
        denominator * center + 0.025 * offset -
        0.015 * offset * I;
    blocks[block].sample_count = samples_per_block;
  }
}

static void split_diagnostic_blocks(
    const MVMCKrylovJackknifeBlock *official_blocks,
    MVMCKrylovJackknifeBlock *diagnostic_blocks) {
  size_t block;
  for (block = 0; block < MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS; ++block) {
    const double signed_offset = (double)((int)(block % 3) - 1);
    const double complex split_delta = 0.08 * signed_offset +
                                       0.04 * signed_offset * I;
    diagnostic_blocks[2 * block].denominator =
        0.5 * official_blocks[block].denominator;
    diagnostic_blocks[2 * block].numerator =
        0.5 * official_blocks[block].numerator + split_delta;
    diagnostic_blocks[2 * block].sample_count =
        MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCK_LENGTH;
    diagnostic_blocks[2 * block + 1].denominator =
        0.5 * official_blocks[block].denominator;
    diagnostic_blocks[2 * block + 1].numerator =
        0.5 * official_blocks[block].numerator - split_delta;
    diagnostic_blocks[2 * block + 1].sample_count =
        MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCK_LENGTH;
  }
}

static void test_jackknife_and_trace_gates(void) {
  MVMCKrylovJackknifeBlock official_blocks[
      MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS];
  MVMCKrylovJackknifeBlock diagnostic_blocks[
      MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS];
  MVMCKrylovJackknifeResult official;
  MVMCKrylovJackknifeResult diagnostic;
  MVMCKrylovJackknifeResult oracle;
  MVMCKrylovJackknifeResult failing_diagnostic;
  MVMCKrylovBlockStabilityResult stability;
  MVMCKrylovTauIntResult tau;
  MVMCKrylovTauIntGateResult tau_gate;
  double tau_values[MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS *
                    MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCK_LENGTH];
  size_t sample;

  fill_jackknife_blocks(
      official_blocks, MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS,
      MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCK_LENGTH, 256.0, 0.25,
      1.2 + 0.4 * I);
  split_diagnostic_blocks(official_blocks, diagnostic_blocks);

  CHECK(mvmc_krylov_jackknife_ratio(
            official_blocks, MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS,
            1.0e-12, &official) == MVMC_KRYLOV_STATUS_OK &&
            official.valid && official.block_count == 16 &&
            official.sample_count == 4096 &&
            official.denominator_stable,
        "official 16x256 jackknife");
  CHECK(mvmc_krylov_jackknife_ratio(
            diagnostic_blocks, MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS,
            1.0e-12, &diagnostic) == MVMC_KRYLOV_STATUS_OK &&
            diagnostic.valid && diagnostic.block_count == 32 &&
            diagnostic.sample_count == 4096 &&
            diagnostic.denominator_stable,
        "diagnostic 32x128 jackknife");
  oracle = jackknife_oracle(
      official_blocks, MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS);
  CHECK(close_complex(official.estimate, oracle.estimate, 1.0e-14),
        "jackknife full estimate oracle");
  CHECK(close_complex(official.jackknife_mean, oracle.jackknife_mean,
                      1.0e-14),
        "jackknife leave-one mean oracle");
  CHECK(close_double(official.variance_real, oracle.variance_real, 1.0e-14),
        "jackknife real variance oracle");
  CHECK(close_double(official.variance_imag, oracle.variance_imag, 1.0e-14),
        "jackknife imag variance oracle");
  CHECK(close_double(official.covariance_real_imag,
                     oracle.covariance_real_imag, 1.0e-14),
        "jackknife covariance oracle");
  {
    MVMCKrylovJackknifeBlock zero_blocks[3] = {
        {2.0 + 0.5 * I, 4.0, 2},
        {0.0, 0.0, 2},
        {1.0 - 0.25 * I, 3.0, 2}};
    MVMCKrylovJackknifeResult zero_jackknife;
    MVMCKrylovJackknifeResult zero_oracle =
        jackknife_oracle(zero_blocks, 3);
    CHECK(mvmc_krylov_jackknife_ratio(
              zero_blocks, 3, 1.0e-12, &zero_jackknife) ==
              MVMC_KRYLOV_STATUS_OK &&
              zero_jackknife.valid && zero_jackknife.denominator_stable,
          "zero-denominator block jackknife");
    CHECK(close_complex(zero_jackknife.estimate,
                        zero_oracle.estimate, 1.0e-14),
          "zero-denominator block estimate oracle");
    CHECK(close_complex(zero_jackknife.jackknife_mean,
                        zero_oracle.jackknife_mean, 1.0e-14),
          "zero-denominator block leave-one oracle");
    zero_blocks[1].numerator = 1.0e-6;
    CHECK(mvmc_krylov_jackknife_ratio(
              zero_blocks, 3, 1.0e-12, &zero_jackknife) ==
              MVMC_KRYLOV_STATUS_NONFINITE,
          "nonzero numerator with zero denominator rejected");
  }

  CHECK(mvmc_krylov_block_stability_check(
            &official, &diagnostic, 1.25, &stability) ==
            MVMC_KRYLOV_STATUS_OK &&
            stability.valid && stability.passed &&
            stability.symmetric_se_ratio <= 1.25 &&
            stability.official_block_count == 16 &&
            stability.diagnostic_block_count == 32,
        "16x256 vs 32x128 block stability pass");
  failing_diagnostic = diagnostic;
  failing_diagnostic.complex_se = official.complex_se * 2.0;
  CHECK(mvmc_krylov_block_stability_check(
            &official, &failing_diagnostic, 1.25, &stability) ==
            MVMC_KRYLOV_STATUS_OK &&
            stability.valid && !stability.passed,
        "block stability fail");

  for (sample = 0; sample < sizeof(tau_values) / sizeof(tau_values[0]);
       ++sample) {
    tau_values[sample] = sample % 2 == 0 ? 1.0 : -1.0;
  }
  CHECK(mvmc_krylov_tau_int_geyer_initial_positive(
            tau_values, sizeof(tau_values) / sizeof(tau_values[0]),
            &tau) == MVMC_KRYLOV_STATUS_OK &&
            tau.valid && tau.tau_int >= 1.0 && tau.tau_int <= 16.0,
        "Geyer IPS tau_int diagnostic");
  CHECK(mvmc_krylov_tau_int_gate_check(
            tau.tau_int, MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCK_LENGTH,
            16.0, 16.0, &tau_gate) == MVMC_KRYLOV_STATUS_OK &&
            tau_gate.valid && tau_gate.passed,
        "official block length covers tau_int");
  CHECK(mvmc_krylov_tau_int_gate_check(
            20.0, MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCK_LENGTH, 16.0,
            16.0, &tau_gate) == MVMC_KRYLOV_STATUS_OK &&
            tau_gate.valid && !tau_gate.passed,
        "maximum tau_int gate fails");
  CHECK(mvmc_krylov_tau_int_gate_check(
            10.0, 100.0, 16.0, 16.0, &tau_gate) ==
            MVMC_KRYLOV_STATUS_OK &&
            tau_gate.valid && !tau_gate.passed,
        "minimum block length gate fails");
}

int main(void) {
  test_upper_triangle_accumulator();
  test_self_normalized_ratio();
  test_jackknife_and_trace_gates();
  return failures == 0 ? 0 : 1;
}
