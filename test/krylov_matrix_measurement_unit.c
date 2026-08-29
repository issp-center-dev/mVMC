#include "krylov_matrix_measurement.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(condition, ...)                                                 \
  do {                                                                        \
    if (!(condition)) {                                                       \
      fprintf(stderr, "KrylovMatrixMeasurement_Unit FAIL: ");                \
      fprintf(stderr, __VA_ARGS__);                                           \
      fprintf(stderr, " (line %d)\n", __LINE__);                            \
      ++failures;                                                             \
    }                                                                         \
  } while (0)

enum {
  ORDER = 3,
  DIMENSION = 3,
  UPPER_COUNT = 6,
  SAMPLE_COUNT = 4,
  OFFICIAL_TRACE_SAMPLES = MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS *
                           MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCK_LENGTH,
  DIAGNOSTIC_TRACE_SAMPLES = MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS *
                             MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCK_LENGTH,
  EXACT_SUPPORT_STATES = 4,
  EXACT_MOMENT_COUNT = 7
};

typedef struct {
  double complex overlap[UPPER_COUNT];
  double complex hamiltonian[UPPER_COUNT];
  double complex hamiltonian_adjoint[UPPER_COUNT];
  double complex hamiltonian_squared[UPPER_COUNT];
  double denominator;
} RawOracleSample;

static int close_double(double actual, double expected, double tolerance) {
  return fabs(actual - expected) <= tolerance * fmax(1.0, fabs(expected));
}

static int close_complex(
    double complex actual, double complex expected, double tolerance) {
  return cabs(actual - expected) <= tolerance * fmax(1.0, cabs(expected));
}

static size_t upper_index(size_t row, size_t column) {
  size_t index = 0;
  CHECK(mvmc_krylov_streaming_upper_index(
            DIMENSION, row, column, &index) == MVMC_KRYLOV_STATUS_OK,
        "oracle upper index");
  return index;
}

static MVMCKrylovMatrixMeasurementPolicy policy_fixture(void) {
  MVMCKrylovMatrixMeasurementPolicy policy;
  memset(&policy, 0, sizeof(policy));
  policy.order = ORDER;
  policy.eta = 0.125;
  policy.guide_lambda[0] = 1.0;
  policy.guide_lambda[1] = 0.5;
  policy.guide_lambda[2] = 2.0;
  policy.guide_lambda[3] = 1.5;
  policy.target_weight[0] = 1.0;
  policy.target_weight[1] = 1.0;
  policy.target_weight[2] = 1.0;
  policy.target_weight[3] = 1.0;
  policy.log_basis_scale[0] = 700.0;
  policy.log_basis_scale[1] = 699.25;
  policy.log_basis_scale[2] = 700.5;
  policy.log_basis_scale[3] = 699.75;
  return policy;
}

static MVMCKrylovMatrixMeasurementPolicy one_sided_policy_fixture(void) {
  MVMCKrylovMatrixMeasurementPolicy policy;
  memset(&policy, 0, sizeof(policy));
  policy.order = ORDER;
  policy.eta = 0.25;
  policy.guide_lambda[0] = 0.75;
  policy.target_weight[0] = 1.0;
  return policy;
}

static void fill_raw_dimensionless(double complex raw[SAMPLE_COUNT][4]) {
  raw[0][0] = 1.0 + 0.5 * I;
  raw[0][1] = -0.25 + 0.75 * I;
  raw[0][2] = 0.4 - 0.2 * I;
  raw[0][3] = 1.2 + 0.1 * I;

  raw[1][0] = 0.8 - 0.1 * I;
  raw[1][1] = 0.0;
  raw[1][2] = -0.6 + 0.3 * I;
  raw[1][3] = 0.5 + 0.2 * I;

  raw[2][0] = -0.3 + 0.9 * I;
  raw[2][1] = 0.7 + 0.4 * I;
  raw[2][2] = 0.2 + 0.8 * I;
  raw[2][3] = -0.4 + 0.6 * I;

  raw[3][0] = 1.4 - 0.2 * I;
  raw[3][1] = -0.5 - 0.3 * I;
  raw[3][2] = 0.0;
  raw[3][3] = 0.25 + 0.75 * I;
}

static void fill_trace_raw_dimensionless(
    size_t sample, double complex raw[4]) {
  const double a = (double)(sample % 257);
  const double b = (double)((37 * sample + 11) % 251);
  const double c = (double)((17 * sample + 23) % 263);
  raw[0] = 1.0 + 1.0e-3 * a + (0.2 + 7.0e-4 * b) * I;
  raw[1] = -0.4 + 8.0e-4 * b + (0.7 - 5.0e-4 * c) * I;
  raw[2] = 0.3 - 6.0e-4 * c + (-0.25 + 9.0e-4 * a) * I;
  raw[3] = 0.85 + 4.0e-4 * (a - b) +
           (0.1 + 3.0e-4 * (b - c)) * I;
}

static void make_scaled_values(
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    const double complex raw[4], MVMCScaledComplex scaled[4]) {
  int index;
  for (index = 0; index <= ORDER; ++index) {
    if (creal(raw[index]) == 0.0 && cimag(raw[index]) == 0.0) {
      CHECK(mvmc_scaled_complex_make_exact_zero(&scaled[index]) ==
                MVMC_PFAFFIAN_STATUS_OK,
            "make exact-zero scaled fixture");
    } else {
      const double magnitude = cabs(raw[index]);
      CHECK(mvmc_scaled_complex_make_finite(
                raw[index] / magnitude,
                log(magnitude) - policy->log_basis_scale[index],
                -INFINITY, &scaled[index]) == MVMC_PFAFFIAN_STATUS_OK,
            "make finite scaled fixture");
    }
  }
}

static void matrix_vector_product4(
    const double complex matrix[EXACT_SUPPORT_STATES][EXACT_SUPPORT_STATES],
    const double complex input[EXACT_SUPPORT_STATES],
    double complex output[EXACT_SUPPORT_STATES]) {
  int row;
  for (row = 0; row < EXACT_SUPPORT_STATES; ++row) {
    int column;
    output[row] = 0.0;
    for (column = 0; column < EXACT_SUPPORT_STATES; ++column) {
      output[row] += matrix[row][column] * input[column];
    }
  }
}

static void fill_one_sided_exact_support(
    double complex powers[EXACT_MOMENT_COUNT][EXACT_SUPPORT_STATES],
    double complex moments[EXACT_MOMENT_COUNT]) {
  const double complex h[EXACT_SUPPORT_STATES][EXACT_SUPPORT_STATES] = {
      {0.5, 0.3 + 0.2 * I, 0.0, -0.1 + 0.05 * I},
      {0.3 - 0.2 * I, -0.25, 0.4 - 0.1 * I, 0.2},
      {0.0, 0.4 + 0.1 * I, 0.75, -0.35 + 0.15 * I},
      {-0.1 - 0.05 * I, 0.2, -0.35 - 0.15 * I, -0.4},
  };
  int moment;
  int state;
  powers[0][0] = 1.0;
  powers[0][1] = I;
  powers[0][2] = -1.0;
  powers[0][3] = -I;
  for (moment = 1; moment < EXACT_MOMENT_COUNT; ++moment) {
    matrix_vector_product4(h, powers[moment - 1], powers[moment]);
  }
  for (moment = 0; moment < EXACT_MOMENT_COUNT; ++moment) {
    moments[moment] = 0.0;
    for (state = 0; state < EXACT_SUPPORT_STATES; ++state) {
      moments[moment] += conj(powers[0][state]) * powers[moment][state];
    }
  }
}

static double guide_value(
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    const double complex raw[4]) {
  double guide = policy->eta;
  int index;
  for (index = 0; index <= ORDER; ++index) {
    guide += policy->guide_lambda[index] *
             creal(conj(raw[index]) * raw[index]);
  }
  return guide;
}

static RawOracleSample raw_oracle_sample(
    const MVMCKrylovMatrixMeasurementPolicy *policy,
    const double complex raw[4]) {
  RawOracleSample oracle;
  const double guide = guide_value(policy, raw);
  size_t row;
  size_t column;
  int index;
  memset(&oracle, 0, sizeof(oracle));
  for (index = 0; index <= ORDER; ++index) {
    oracle.denominator += policy->target_weight[index] *
                          creal(conj(raw[index]) * raw[index]) / guide;
  }
  for (row = 0; row < DIMENSION; ++row) {
    for (column = row; column < DIMENSION; ++column) {
      const size_t entry = upper_index(row, column);
      const double right_ratio = exp(
          policy->log_basis_scale[column] -
          policy->log_basis_scale[column + 1]);
      const double left_ratio = exp(
          policy->log_basis_scale[row] -
          policy->log_basis_scale[row + 1]);
      oracle.overlap[entry] =
          conj(raw[row]) * raw[column] / guide;
      oracle.hamiltonian[entry] =
          right_ratio * conj(raw[row]) * raw[column + 1] / guide;
      oracle.hamiltonian_adjoint[entry] =
          left_ratio * conj(raw[column]) * raw[row + 1] / guide;
      oracle.hamiltonian_squared[entry] =
          left_ratio * right_ratio * conj(raw[row + 1]) *
          raw[column + 1] / guide;
    }
  }
  return oracle;
}

static void add_oracle_sample(
    RawOracleSample *sum, const RawOracleSample *sample) {
  size_t entry;
  sum->denominator += sample->denominator;
  for (entry = 0; entry < UPPER_COUNT; ++entry) {
    sum->overlap[entry] += sample->overlap[entry];
    sum->hamiltonian[entry] += sample->hamiltonian[entry];
    sum->hamiltonian_adjoint[entry] += sample->hamiltonian_adjoint[entry];
    sum->hamiltonian_squared[entry] += sample->hamiltonian_squared[entry];
  }
}

static void scale_oracle(RawOracleSample *sum, double factor) {
  size_t entry;
  sum->denominator *= factor;
  for (entry = 0; entry < UPPER_COUNT; ++entry) {
    sum->overlap[entry] *= factor;
    sum->hamiltonian[entry] *= factor;
    sum->hamiltonian_adjoint[entry] *= factor;
    sum->hamiltonian_squared[entry] *= factor;
  }
}

static void test_single_sample_against_raw_oracle(void) {
  const MVMCKrylovMatrixMeasurementPolicy policy = policy_fixture();
  double complex raw[SAMPLE_COUNT][4];
  MVMCScaledComplex scaled[4];
  RawOracleSample oracle;
  MVMCKrylovMatrixMeasurementSampleDiagnostics diagnostics;
  double complex overlap[UPPER_COUNT];
  double complex hamiltonian[UPPER_COUNT];
  double complex hamiltonian_adjoint[UPPER_COUNT];
  double complex hamiltonian_squared[UPPER_COUNT];
  size_t dimension = 0;
  size_t upper_count = 0;
  size_t entry;

  fill_raw_dimensionless(raw);
  make_scaled_values(&policy, raw[0], scaled);
  oracle = raw_oracle_sample(&policy, raw[0]);

  CHECK(mvmc_krylov_matrix_measurement_dimension(
            ORDER, &dimension, &upper_count) == MVMC_KRYLOV_STATUS_OK &&
            dimension == DIMENSION && upper_count == UPPER_COUNT,
        "matrix measurement dimension");
  CHECK(mvmc_krylov_matrix_measurement_sample_with_adjoint(
            &policy, scaled, 4, overlap, hamiltonian,
            hamiltonian_adjoint, hamiltonian_squared, UPPER_COUNT,
            &diagnostics) ==
            MVMC_KRYLOV_STATUS_OK &&
            diagnostics.valid && diagnostics.finite_component_count == 4 &&
            diagnostics.zero_component_count == 0,
        "single sample measurement");
  CHECK(close_double(exp(diagnostics.log_guide),
                     guide_value(&policy, raw[0]), 1.0e-13),
        "log guide raw oracle");
  CHECK(close_double(diagnostics.denominator, oracle.denominator, 1.0e-13),
        "sample denominator oracle");
  for (entry = 0; entry < UPPER_COUNT; ++entry) {
    CHECK(close_complex(overlap[entry], oracle.overlap[entry], 1.0e-13),
          "sample S entry %zu", entry);
    CHECK(close_complex(hamiltonian[entry], oracle.hamiltonian[entry],
                        1.0e-13),
          "sample K entry %zu", entry);
    CHECK(close_complex(hamiltonian_adjoint[entry],
                        oracle.hamiltonian_adjoint[entry], 1.0e-13),
          "sample K-adjoint entry %zu", entry);
    CHECK(close_complex(hamiltonian_squared[entry],
                        oracle.hamiltonian_squared[entry], 1.0e-13),
          "sample B entry %zu", entry);
  }
}

static void test_accumulator_block_totals(void) {
  const MVMCKrylovMatrixMeasurementPolicy policy = policy_fixture();
  double complex raw[SAMPLE_COUNT][4];
  MVMCKrylovStreamingComplexSum overlap_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_squared_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum left_overlap_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum left_hamiltonian_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum left_hamiltonian_squared_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum right_overlap_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum right_hamiltonian_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum right_hamiltonian_squared_entries[UPPER_COUNT];
  MVMCKrylovMatrixMeasurementAccumulator accumulator;
  MVMCKrylovMatrixMeasurementAccumulator left;
  MVMCKrylovMatrixMeasurementAccumulator right;
  RawOracleSample oracle_sum;
  MVMCKrylovJackknifeBlock blocks[2];
  MVMCKrylovJackknifeResult jackknife;
  double complex overlap_mean[UPPER_COUNT];
  double complex hamiltonian_mean[UPPER_COUNT];
  double complex hamiltonian_squared_mean[UPPER_COUNT];
  double denominator_mean = NAN;
  size_t entry;
  size_t sample;

  memset(&oracle_sum, 0, sizeof(oracle_sum));
  fill_raw_dimensionless(raw);
  CHECK(mvmc_krylov_matrix_measurement_accumulator_init(
            ORDER, overlap_entries, hamiltonian_entries,
            hamiltonian_squared_entries, UPPER_COUNT, &accumulator) ==
            MVMC_KRYLOV_STATUS_OK &&
            accumulator.valid,
        "full measurement accumulator init");
  CHECK(mvmc_krylov_matrix_measurement_accumulator_init(
            ORDER, left_overlap_entries, left_hamiltonian_entries,
            left_hamiltonian_squared_entries, UPPER_COUNT, &left) ==
            MVMC_KRYLOV_STATUS_OK &&
            left.valid,
        "left measurement accumulator init");
  CHECK(mvmc_krylov_matrix_measurement_accumulator_init(
            ORDER, right_overlap_entries, right_hamiltonian_entries,
            right_hamiltonian_squared_entries, UPPER_COUNT, &right) ==
            MVMC_KRYLOV_STATUS_OK &&
            right.valid,
        "right measurement accumulator init");

  for (sample = 0; sample < SAMPLE_COUNT; ++sample) {
    MVMCScaledComplex scaled[4];
    const RawOracleSample oracle = raw_oracle_sample(&policy, raw[sample]);
    add_oracle_sample(&oracle_sum, &oracle);
    make_scaled_values(&policy, raw[sample], scaled);
    CHECK(mvmc_krylov_matrix_measurement_accumulator_add_sample(
              &accumulator, &policy, scaled, 4, NULL) ==
              MVMC_KRYLOV_STATUS_OK,
          "full accumulator sample %zu", sample);
    CHECK(mvmc_krylov_matrix_measurement_accumulator_add_sample(
              sample < 2 ? &left : &right, &policy, scaled, 4, NULL) ==
              MVMC_KRYLOV_STATUS_OK,
          "partition accumulator sample %zu", sample);
  }
  scale_oracle(&oracle_sum, 1.0 / (double)SAMPLE_COUNT);

  CHECK(mvmc_krylov_matrix_measurement_accumulator_mean(
            &accumulator, overlap_mean, hamiltonian_mean,
            hamiltonian_squared_mean, UPPER_COUNT, &denominator_mean) ==
            MVMC_KRYLOV_STATUS_OK,
        "matrix accumulator mean");
  CHECK(close_double(denominator_mean, oracle_sum.denominator, 1.0e-13),
        "mean denominator oracle");
  for (entry = 0; entry < UPPER_COUNT; ++entry) {
    CHECK(close_complex(overlap_mean[entry], oracle_sum.overlap[entry],
                        1.0e-13),
          "mean S entry %zu", entry);
    CHECK(close_complex(hamiltonian_mean[entry],
                        oracle_sum.hamiltonian[entry], 1.0e-13),
          "mean K entry %zu", entry);
    CHECK(close_complex(hamiltonian_squared_mean[entry],
                        oracle_sum.hamiltonian_squared[entry], 1.0e-13),
          "mean B entry %zu", entry);
  }

  CHECK(mvmc_krylov_matrix_measurement_extract_block(
            &left, MVMC_KRYLOV_MATRIX_HAMILTONIAN, 0, 1,
            &blocks[0]) == MVMC_KRYLOV_STATUS_OK &&
            blocks[0].sample_count == 2,
        "left block extraction");
  CHECK(mvmc_krylov_matrix_measurement_extract_block(
            &right, MVMC_KRYLOV_MATRIX_HAMILTONIAN, 0, 1,
            &blocks[1]) == MVMC_KRYLOV_STATUS_OK &&
            blocks[1].sample_count == 2,
        "right block extraction");
  CHECK(mvmc_krylov_matrix_measurement_extract_block(
            &left, MVMC_KRYLOV_MATRIX_HAMILTONIAN_ADJOINT, 0, 1,
            &blocks[0]) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "K-adjoint extraction requires diagnostic storage");
  CHECK(blocks[0].sample_count == 0 && blocks[0].numerator == 0.0 &&
            blocks[0].denominator == 0.0,
        "failed K-adjoint extraction clears output");
  CHECK(mvmc_krylov_matrix_measurement_extract_block(
            &left, MVMC_KRYLOV_MATRIX_HAMILTONIAN, 0, 1,
            &blocks[0]) == MVMC_KRYLOV_STATUS_OK,
        "left block extraction restored after negative check");
  CHECK(mvmc_krylov_jackknife_ratio(
            blocks, 2, 1.0e-12, &jackknife) ==
            MVMC_KRYLOV_STATUS_OK &&
            jackknife.valid,
        "block totals feed jackknife ratio");
  entry = upper_index(0, 1);
  CHECK(close_complex(
            jackknife.estimate,
            (blocks[0].numerator + blocks[1].numerator) /
                (blocks[0].denominator + blocks[1].denominator),
            1.0e-15),
        "jackknife estimate from extracted blocks");
  CHECK(close_complex(
            jackknife.estimate,
            accumulator.hamiltonian.entries[entry].real.sum /
                    (blocks[0].denominator + blocks[1].denominator) +
                I * accumulator.hamiltonian.entries[entry].imag.sum /
                    (blocks[0].denominator + blocks[1].denominator),
            1.0e-13),
        "full accumulator numerator matches block ratio");
}

static void test_continuous_trace_block_accumulators(void) {
  const MVMCKrylovMatrixMeasurementPolicy policy = policy_fixture();
  MVMCKrylovMatrixMeasurementAccumulator full_accumulator;
  MVMCKrylovStreamingComplexSum full_overlap_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum full_hamiltonian_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum full_hamiltonian_squared_entries[UPPER_COUNT];
  MVMCKrylovMatrixMeasurementAccumulator
      official_block_accumulators[MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS];
  MVMCKrylovStreamingComplexSum official_overlap_entries[
      MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS * UPPER_COUNT];
  MVMCKrylovStreamingComplexSum official_hamiltonian_entries[
      MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS * UPPER_COUNT];
  MVMCKrylovStreamingComplexSum official_hamiltonian_squared_entries[
      MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS * UPPER_COUNT];
  MVMCKrylovMatrixMeasurementBlockAccumulator official;
  MVMCKrylovJackknifeBlock
      official_blocks[MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS];
  MVMCKrylovJackknifeResult official_jackknife;
  MVMCKrylovMatrixMeasurementAccumulator
      diagnostic_block_accumulators[MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS];
  MVMCKrylovStreamingComplexSum diagnostic_overlap_entries[
      MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS * UPPER_COUNT];
  MVMCKrylovStreamingComplexSum diagnostic_hamiltonian_entries[
      MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS * UPPER_COUNT];
  MVMCKrylovStreamingComplexSum diagnostic_hamiltonian_squared_entries[
      MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS * UPPER_COUNT];
  MVMCKrylovMatrixMeasurementBlockAccumulator diagnostic;
  MVMCKrylovJackknifeBlock
      diagnostic_blocks[MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS];
  MVMCKrylovJackknifeResult diagnostic_jackknife;
  MVMCKrylovBlockStabilityResult stability;
  MVMCKrylovJackknifeBlock full_block;
  size_t sample;
  size_t block;

  CHECK(OFFICIAL_TRACE_SAMPLES == DIAGNOSTIC_TRACE_SAMPLES,
        "official and diagnostic traces must cover the same sample count");
  CHECK(mvmc_krylov_matrix_measurement_accumulator_init(
            ORDER, full_overlap_entries, full_hamiltonian_entries,
            full_hamiltonian_squared_entries, UPPER_COUNT,
            &full_accumulator) == MVMC_KRYLOV_STATUS_OK &&
            full_accumulator.valid,
        "full trace accumulator init");
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_init(
            ORDER, MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS,
            MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCK_LENGTH,
            official_block_accumulators, official_overlap_entries,
            official_hamiltonian_entries, NULL,
            official_hamiltonian_squared_entries,
            MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS * UPPER_COUNT,
            &official) == MVMC_KRYLOV_STATUS_OK &&
            official.valid &&
            official.block_count == MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS,
        "official block accumulator init");
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_init(
            ORDER, MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS,
            MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCK_LENGTH,
            diagnostic_block_accumulators, diagnostic_overlap_entries,
            diagnostic_hamiltonian_entries, NULL,
            diagnostic_hamiltonian_squared_entries,
            MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS * UPPER_COUNT,
            &diagnostic) == MVMC_KRYLOV_STATUS_OK &&
            diagnostic.valid &&
            diagnostic.block_count == MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS,
        "diagnostic block accumulator init");

  for (sample = 0; sample < OFFICIAL_TRACE_SAMPLES; ++sample) {
    double complex raw[4];
    MVMCScaledComplex scaled[4];
    fill_trace_raw_dimensionless(sample, raw);
    make_scaled_values(&policy, raw, scaled);
    CHECK(mvmc_krylov_matrix_measurement_accumulator_add_sample(
              &full_accumulator, &policy, scaled, 4, NULL) ==
              MVMC_KRYLOV_STATUS_OK,
          "full trace sample %zu", sample);
    CHECK(mvmc_krylov_matrix_measurement_block_accumulator_add_sample(
              &official, &policy, scaled, 4, NULL) ==
              MVMC_KRYLOV_STATUS_OK,
          "official trace sample %zu", sample);
    CHECK(mvmc_krylov_matrix_measurement_block_accumulator_add_sample(
              &diagnostic, &policy, scaled, 4, NULL) ==
              MVMC_KRYLOV_STATUS_OK,
          "diagnostic trace sample %zu", sample);
  }

  CHECK(official.sample_count == OFFICIAL_TRACE_SAMPLES,
        "official block sample count");
  CHECK(diagnostic.sample_count == DIAGNOSTIC_TRACE_SAMPLES,
        "diagnostic block sample count");
  for (block = 0; block < MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS; ++block) {
    CHECK(official_block_accumulators[block].sample_count ==
              MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCK_LENGTH,
          "official block length %zu", block);
  }
  for (block = 0; block < MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS; ++block) {
    CHECK(diagnostic_block_accumulators[block].sample_count ==
              MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCK_LENGTH,
          "diagnostic block length %zu", block);
  }

  CHECK(mvmc_krylov_matrix_measurement_extract_block(
            &full_accumulator, MVMC_KRYLOV_MATRIX_HAMILTONIAN,
            1, 2, &full_block) == MVMC_KRYLOV_STATUS_OK,
        "full trace block extraction");
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_extract_blocks(
            &official, MVMC_KRYLOV_MATRIX_HAMILTONIAN, 1, 2,
            official_blocks, MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS) ==
            MVMC_KRYLOV_STATUS_OK,
        "official block extraction");
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_extract_blocks(
            &diagnostic, MVMC_KRYLOV_MATRIX_HAMILTONIAN, 1, 2,
            diagnostic_blocks, MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS) ==
            MVMC_KRYLOV_STATUS_OK,
        "diagnostic block extraction");
  CHECK(mvmc_krylov_jackknife_ratio(
            official_blocks, MVMC_KRYLOV_OFFICIAL_JACKKNIFE_BLOCKS,
            1.0e-12, &official_jackknife) == MVMC_KRYLOV_STATUS_OK &&
            official_jackknife.valid,
        "official jackknife ratio");
  CHECK(mvmc_krylov_jackknife_ratio(
            diagnostic_blocks, MVMC_KRYLOV_DIAGNOSTIC_JACKKNIFE_BLOCKS,
            1.0e-12, &diagnostic_jackknife) == MVMC_KRYLOV_STATUS_OK &&
            diagnostic_jackknife.valid,
        "diagnostic jackknife ratio");
  CHECK(close_complex(official_jackknife.estimate,
                      full_block.numerator / full_block.denominator,
                      1.0e-13),
        "official estimate matches full trace");
  CHECK(close_complex(diagnostic_jackknife.estimate,
                      full_block.numerator / full_block.denominator,
                      1.0e-13),
        "diagnostic estimate matches full trace");
  /* This deterministic trace has partition aliasing; P4-S owns the physics gate. */
  CHECK(mvmc_krylov_block_stability_check(
            &official_jackknife, &diagnostic_jackknife, 32.0,
            &stability) == MVMC_KRYLOV_STATUS_OK &&
            stability.valid && stability.passed,
        "official/diagnostic stability check");

  {
    double complex raw[4];
    MVMCScaledComplex scaled[4];
    fill_trace_raw_dimensionless(OFFICIAL_TRACE_SAMPLES, raw);
    make_scaled_values(&policy, raw, scaled);
    CHECK(mvmc_krylov_matrix_measurement_block_accumulator_add_sample(
              &official, &policy, scaled, 4, NULL) ==
              MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
          "official block accumulator rejects overflow");
  }
}

static void test_one_sided_moment_oracle(void) {
  const MVMCKrylovMatrixMeasurementPolicy policy =
      one_sided_policy_fixture();
  double complex powers[EXACT_MOMENT_COUNT][EXACT_SUPPORT_STATES];
  double complex moments[EXACT_MOMENT_COUNT];
  MVMCKrylovStreamingComplexSum overlap_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_squared_entries[UPPER_COUNT];
  MVMCKrylovMatrixMeasurementAccumulator accumulator;
  double complex overlap_mean[UPPER_COUNT];
  double complex hamiltonian_mean[UPPER_COUNT];
  double complex hamiltonian_squared_mean[UPPER_COUNT];
  double denominator_mean = NAN;
  int state;
  size_t row;
  size_t column;

  fill_one_sided_exact_support(powers, moments);
  CHECK(close_complex(moments[0], (double)EXACT_SUPPORT_STATES,
                      1.0e-13),
        "mu_0 exact support norm");
  CHECK(mvmc_krylov_matrix_measurement_accumulator_init(
            ORDER, overlap_entries, hamiltonian_entries,
            hamiltonian_squared_entries, UPPER_COUNT,
            &accumulator) == MVMC_KRYLOV_STATUS_OK &&
            accumulator.valid,
        "one-sided accumulator init");
  for (state = 0; state < EXACT_SUPPORT_STATES; ++state) {
    double complex raw[4];
    MVMCScaledComplex scaled[4];
    raw[0] = powers[0][state];
    raw[1] = powers[1][state];
    raw[2] = powers[2][state];
    raw[3] = powers[3][state];
    make_scaled_values(&policy, raw, scaled);
    CHECK(mvmc_krylov_matrix_measurement_accumulator_add_sample(
              &accumulator, &policy, scaled, 4, NULL) ==
              MVMC_KRYLOV_STATUS_OK,
          "one-sided support state %d", state);
  }
  CHECK(mvmc_krylov_matrix_measurement_accumulator_mean(
            &accumulator, overlap_mean, hamiltonian_mean,
            hamiltonian_squared_mean, UPPER_COUNT, &denominator_mean) ==
            MVMC_KRYLOV_STATUS_OK,
        "one-sided accumulator mean");
  CHECK(close_double(denominator_mean, 1.0, 1.0e-13),
        "one-sided denominator mean");
  for (row = 0; row < DIMENSION; ++row) {
    for (column = row; column < DIMENSION; ++column) {
      const size_t entry = upper_index(row, column);
      CHECK(close_complex(overlap_mean[entry] / denominator_mean,
                          moments[row + column] / moments[0],
                          2.0e-13),
            "S_%zu%zu one-sided mu", row, column);
      CHECK(close_complex(hamiltonian_mean[entry] / denominator_mean,
                          moments[row + column + 1] / moments[0],
                          2.0e-13),
            "K_%zu%zu one-sided mu", row, column);
      CHECK(close_complex(hamiltonian_squared_mean[entry] /
                              denominator_mean,
                          moments[row + column + 2] / moments[0],
                          2.0e-13),
            "B_%zu%zu one-sided mu", row, column);
    }
  }
}

static double raw_hamiltonian_residual(const RawOracleSample *mean,
                                       double *norm) {
  double complex dense[DIMENSION][DIMENSION] = {{0.0}};
  double residual_squared = 0.0;
  double norm_squared = 0.0;
  size_t row;
  size_t column;
  for (row = 0; row < DIMENSION; ++row) {
    for (column = row; column < DIMENSION; ++column) {
      const size_t entry = upper_index(row, column);
      dense[row][column] = mean->hamiltonian[entry];
      if (row != column) {
        dense[column][row] = mean->hamiltonian_adjoint[entry];
      }
    }
  }
  for (row = 0; row < DIMENSION; ++row) {
    for (column = 0; column < DIMENSION; ++column) {
      const double complex value = dense[row][column];
      const double complex delta = value - conj(dense[column][row]);
      residual_squared += creal(delta) * creal(delta) +
                          cimag(delta) * cimag(delta);
      norm_squared += creal(value) * creal(value) +
                      cimag(value) * cimag(value);
    }
  }
  *norm = sqrt(norm_squared);
  return *norm > 0.0 ? sqrt(residual_squared) / *norm
                     : sqrt(residual_squared);
}

static void test_offdiagonal_antihermitian_frobenius(void) {
  MVMCKrylovStreamingComplexSum overlap_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_adjoint_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_squared_entries[UPPER_COUNT];
  MVMCKrylovMatrixMeasurementAccumulator accumulator;
  MVMCKrylovMatrixMeasurementDiagnosticsSummary summary;
  const size_t entry = upper_index(0, 1);

  CHECK(mvmc_krylov_matrix_measurement_accumulator_init_with_diagnostics(
            ORDER, overlap_entries, hamiltonian_entries,
            hamiltonian_adjoint_entries, hamiltonian_squared_entries,
            UPPER_COUNT, &accumulator) == MVMC_KRYLOV_STATUS_OK,
        "off-diagonal residual accumulator init");
  accumulator.sample_count = 2;
  accumulator.overlap.sample_count = 2;
  accumulator.hamiltonian.sample_count = 2;
  accumulator.hamiltonian_adjoint.sample_count = 2;
  accumulator.hamiltonian_squared.sample_count = 2;
  accumulator.hamiltonian.entries[entry].real.sum = 2.0;
  accumulator.denominator.sum = 2.0;
  accumulator.denominator_squared.sum = 2.0;
  accumulator.minimum_denominator = 1.0;
  accumulator.maximum_denominator = 1.0;
  accumulator.minimum_log_contribution = 0.0;
  accumulator.maximum_log_contribution = 0.0;

  CHECK(mvmc_krylov_matrix_measurement_diagnostics_summary(
            &accumulator, 1.0e-12, 1.0, &summary) ==
            MVMC_KRYLOV_STATUS_OK &&
            summary.valid && summary.hamiltonian_residual_available,
        "off-diagonal residual summary");
  CHECK(close_double(summary.hamiltonian_norm, 1.0, 1.0e-15),
        "off-diagonal full-matrix norm");
  CHECK(close_double(summary.hamiltonian_antihermitian_residual,
                     sqrt(2.0), 1.0e-15),
        "off-diagonal Frobenius residual counts both matrix halves");
}

static void test_zero_target_samples_are_accepted(void) {
  const MVMCKrylovMatrixMeasurementPolicy policy = policy_fixture();
  double complex raw[SAMPLE_COUNT][4];
  double complex zero_raw[4] = {0.0, 0.0, 0.0, 0.0};
  MVMCScaledComplex zero_scaled[4];
  MVMCKrylovMatrixMeasurementSampleDiagnostics diagnostics;
  double complex overlap[UPPER_COUNT];
  double complex hamiltonian[UPPER_COUNT];
  double complex hamiltonian_adjoint[UPPER_COUNT];
  double complex hamiltonian_squared[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum overlap_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_adjoint_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_squared_entries[UPPER_COUNT];
  MVMCKrylovMatrixMeasurementAccumulator accumulator;
  MVMCKrylovMatrixMeasurementDiagnosticsSummary summary;
  RawOracleSample oracle_sum;
  RawOracleSample oracle_mean;
  double complex overlap_mean[UPPER_COUNT];
  double complex hamiltonian_mean[UPPER_COUNT];
  double complex hamiltonian_squared_mean[UPPER_COUNT];
  double denominator_mean = NAN;
  double expected_norm = NAN;
  double expected_residual;
  size_t entry;

  fill_raw_dimensionless(raw);
  make_scaled_values(&policy, zero_raw, zero_scaled);
  CHECK(mvmc_krylov_matrix_measurement_sample_with_adjoint(
            &policy, zero_scaled, 4, overlap, hamiltonian,
            hamiltonian_adjoint, hamiltonian_squared, UPPER_COUNT,
            &diagnostics) == MVMC_KRYLOV_STATUS_OK &&
            diagnostics.valid && diagnostics.zero_target_sample &&
            diagnostics.denominator == 0.0 &&
            diagnostics.target_over_guide == 0.0 &&
            diagnostics.guide_signal_over_guide == 0.0 &&
            diagnostics.finite_component_count == 0 &&
            diagnostics.zero_component_count == 4,
        "zero-target sample accepted");
  CHECK(close_double(exp(diagnostics.log_guide), policy.eta, 1.0e-13),
        "zero-target guide floor");
  for (entry = 0; entry < UPPER_COUNT; ++entry) {
    CHECK(overlap[entry] == 0.0 && hamiltonian[entry] == 0.0 &&
              hamiltonian_adjoint[entry] == 0.0 &&
              hamiltonian_squared[entry] == 0.0,
          "zero-target contribution entry %zu", entry);
  }

  CHECK(mvmc_krylov_matrix_measurement_accumulator_init_with_diagnostics(
            ORDER, overlap_entries, hamiltonian_entries,
            hamiltonian_adjoint_entries, hamiltonian_squared_entries,
            UPPER_COUNT, &accumulator) == MVMC_KRYLOV_STATUS_OK &&
            accumulator.valid,
        "zero-target accumulator init");
  memset(&oracle_sum, 0, sizeof(oracle_sum));
  {
    MVMCScaledComplex scaled[4];
    const RawOracleSample first = raw_oracle_sample(&policy, raw[0]);
    make_scaled_values(&policy, raw[0], scaled);
    add_oracle_sample(&oracle_sum, &first);
    CHECK(mvmc_krylov_matrix_measurement_accumulator_add_sample(
              &accumulator, &policy, scaled, 4, NULL) ==
              MVMC_KRYLOV_STATUS_OK,
          "normal sample before zero-target");
  }
  CHECK(mvmc_krylov_matrix_measurement_accumulator_add_sample(
            &accumulator, &policy, zero_scaled, 4, &diagnostics) ==
            MVMC_KRYLOV_STATUS_OK &&
            diagnostics.valid && diagnostics.zero_target_sample,
        "zero-target accumulator sample");
  {
    MVMCScaledComplex scaled[4];
    const RawOracleSample second = raw_oracle_sample(&policy, raw[1]);
    make_scaled_values(&policy, raw[1], scaled);
    add_oracle_sample(&oracle_sum, &second);
    CHECK(mvmc_krylov_matrix_measurement_accumulator_add_sample(
              &accumulator, &policy, scaled, 4, NULL) ==
              MVMC_KRYLOV_STATUS_OK,
          "normal sample after zero-target");
  }

  oracle_mean = oracle_sum;
  scale_oracle(&oracle_mean, 1.0 / 3.0);
  CHECK(accumulator.sample_count == 3 &&
            accumulator.zero_target_sample_count == 1 &&
            accumulator.status == MVMC_KRYLOV_STATUS_OK,
        "zero-target accumulator counts");
  CHECK(mvmc_krylov_matrix_measurement_accumulator_mean(
            &accumulator, overlap_mean, hamiltonian_mean,
            hamiltonian_squared_mean, UPPER_COUNT, &denominator_mean) ==
            MVMC_KRYLOV_STATUS_OK,
        "zero-target accumulator mean");
  CHECK(close_double(denominator_mean, oracle_mean.denominator, 1.0e-13),
        "zero-target denominator mean");
  for (entry = 0; entry < UPPER_COUNT; ++entry) {
    CHECK(close_complex(overlap_mean[entry], oracle_mean.overlap[entry],
                        1.0e-13),
          "zero-target mean S entry %zu", entry);
    CHECK(close_complex(hamiltonian_mean[entry],
                        oracle_mean.hamiltonian[entry], 1.0e-13),
          "zero-target mean K entry %zu", entry);
    CHECK(close_complex(hamiltonian_squared_mean[entry],
                        oracle_mean.hamiltonian_squared[entry], 1.0e-13),
          "zero-target mean B entry %zu", entry);
  }

  expected_residual =
      raw_hamiltonian_residual(&oracle_mean, &expected_norm);
  CHECK(mvmc_krylov_matrix_measurement_diagnostics_summary(
            &accumulator, 1.0e-12, 1.0, &summary) ==
            MVMC_KRYLOV_STATUS_OK &&
            summary.valid && summary.zero_target_sample_count == 1 &&
            close_double(summary.zero_target_sample_fraction, 1.0 / 3.0,
                         1.0e-13) &&
            summary.minimum_denominator == 0.0,
        "zero-target diagnostics summary");
  CHECK(close_double(summary.hamiltonian_antihermitian_residual,
                     expected_residual, 1.0e-13),
        "zero-target anti-Hermitian residual");
}

static void test_zero_target_block_jackknife(void) {
  const MVMCKrylovMatrixMeasurementPolicy policy = policy_fixture();
  enum { ZERO_BLOCK_COUNT = 3, ZERO_BLOCK_LENGTH = 2 };
  double complex raw[SAMPLE_COUNT][4];
  double complex zero_raw[4] = {0.0, 0.0, 0.0, 0.0};
  MVMCScaledComplex zero_scaled[4];
  MVMCKrylovMatrixMeasurementAccumulator
      block_accumulators[ZERO_BLOCK_COUNT];
  MVMCKrylovStreamingComplexSum overlap_entries[
      ZERO_BLOCK_COUNT * UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_entries[
      ZERO_BLOCK_COUNT * UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_squared_entries[
      ZERO_BLOCK_COUNT * UPPER_COUNT];
  MVMCKrylovMatrixMeasurementBlockAccumulator accumulator;
  MVMCKrylovJackknifeBlock blocks[ZERO_BLOCK_COUNT];
  MVMCKrylovJackknifeResult jackknife;
  double complex total_numerator = 0.0;
  double denominator_total = 0.0;
  size_t block;

  fill_raw_dimensionless(raw);
  make_scaled_values(&policy, zero_raw, zero_scaled);
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_init(
            ORDER, ZERO_BLOCK_COUNT, ZERO_BLOCK_LENGTH,
            block_accumulators, overlap_entries, hamiltonian_entries, NULL,
            hamiltonian_squared_entries, ZERO_BLOCK_COUNT * UPPER_COUNT,
            &accumulator) == MVMC_KRYLOV_STATUS_OK &&
            accumulator.valid,
        "zero-target block accumulator init");
  {
    size_t sample;
    for (sample = 0; sample < 2; ++sample) {
      MVMCScaledComplex scaled[4];
      make_scaled_values(&policy, raw[sample], scaled);
      CHECK(mvmc_krylov_matrix_measurement_block_accumulator_add_sample(
                &accumulator, &policy, scaled, 4, NULL) ==
                MVMC_KRYLOV_STATUS_OK,
            "positive block sample %zu", sample);
    }
  }
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_add_sample(
            &accumulator, &policy, zero_scaled, 4, NULL) ==
            MVMC_KRYLOV_STATUS_OK,
        "first zero-target block sample");
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_add_sample(
            &accumulator, &policy, zero_scaled, 4, NULL) ==
            MVMC_KRYLOV_STATUS_OK,
        "second zero-target block sample");
  {
    size_t sample;
    for (sample = 2; sample < 4; ++sample) {
      MVMCScaledComplex scaled[4];
      make_scaled_values(&policy, raw[sample], scaled);
      CHECK(mvmc_krylov_matrix_measurement_block_accumulator_add_sample(
                &accumulator, &policy, scaled, 4, NULL) ==
                MVMC_KRYLOV_STATUS_OK,
            "trailing positive block sample %zu", sample);
    }
  }
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_extract_blocks(
            &accumulator, MVMC_KRYLOV_MATRIX_HAMILTONIAN, 0, 1,
            blocks, ZERO_BLOCK_COUNT) == MVMC_KRYLOV_STATUS_OK,
        "zero-target block extraction");
  CHECK(blocks[1].sample_count == ZERO_BLOCK_LENGTH &&
            blocks[1].denominator == 0.0 && blocks[1].numerator == 0.0,
        "zero-target block total");
  for (block = 0; block < ZERO_BLOCK_COUNT; ++block) {
    total_numerator += blocks[block].numerator;
    denominator_total += blocks[block].denominator;
  }
  CHECK(mvmc_krylov_jackknife_ratio(
            blocks, ZERO_BLOCK_COUNT, 1.0e-12, &jackknife) ==
            MVMC_KRYLOV_STATUS_OK &&
            jackknife.valid && jackknife.denominator_stable,
        "zero-target block jackknife");
  CHECK(close_complex(jackknife.estimate,
                      total_numerator / denominator_total, 1.0e-13),
        "zero-target jackknife estimate");
}

static void test_diagnostics_summary(void) {
  const MVMCKrylovMatrixMeasurementPolicy policy = policy_fixture();
  double complex raw[SAMPLE_COUNT][4];
  MVMCKrylovStreamingComplexSum overlap_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_adjoint_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_squared_entries[UPPER_COUNT];
  MVMCKrylovMatrixMeasurementAccumulator accumulator;
  MVMCKrylovMatrixMeasurementDiagnosticsSummary summary;
  MVMCKrylovJackknifeBlock adjoint_block;
  RawOracleSample oracle_sum;
  RawOracleSample oracle_mean;
  double denominator_sum = 0.0;
  double denominator_squared_sum = 0.0;
  double minimum_denominator = INFINITY;
  double maximum_denominator = 0.0;
  double expected_norm = NAN;
  double expected_residual;
  size_t sample;

  memset(&oracle_sum, 0, sizeof(oracle_sum));
  fill_raw_dimensionless(raw);
  CHECK(mvmc_krylov_matrix_measurement_accumulator_init_with_diagnostics(
            ORDER, overlap_entries, hamiltonian_entries,
            hamiltonian_adjoint_entries, hamiltonian_squared_entries,
            UPPER_COUNT, &accumulator) == MVMC_KRYLOV_STATUS_OK &&
            accumulator.valid,
        "diagnostics accumulator init");
  for (sample = 0; sample < SAMPLE_COUNT; ++sample) {
    MVMCScaledComplex scaled[4];
    MVMCKrylovMatrixMeasurementSampleDiagnostics diagnostics;
    const RawOracleSample oracle = raw_oracle_sample(&policy, raw[sample]);
    add_oracle_sample(&oracle_sum, &oracle);
    denominator_sum += oracle.denominator;
    denominator_squared_sum += oracle.denominator * oracle.denominator;
    minimum_denominator = fmin(minimum_denominator, oracle.denominator);
    maximum_denominator = fmax(maximum_denominator, oracle.denominator);
    make_scaled_values(&policy, raw[sample], scaled);
    CHECK(mvmc_krylov_matrix_measurement_accumulator_add_sample(
              &accumulator, &policy, scaled, 4, &diagnostics) ==
              MVMC_KRYLOV_STATUS_OK &&
              diagnostics.valid,
          "diagnostics accumulator sample %zu", sample);
  }
  oracle_mean = oracle_sum;
  scale_oracle(&oracle_mean, 1.0 / (double)SAMPLE_COUNT);
  expected_residual =
      raw_hamiltonian_residual(&oracle_mean, &expected_norm);

  CHECK(mvmc_krylov_matrix_measurement_extract_block(
            &accumulator, MVMC_KRYLOV_MATRIX_HAMILTONIAN_ADJOINT, 0, 1,
            &adjoint_block) == MVMC_KRYLOV_STATUS_OK &&
            adjoint_block.sample_count == SAMPLE_COUNT,
        "K-adjoint diagnostic block extraction");
  CHECK(close_complex(adjoint_block.numerator,
                      oracle_sum.hamiltonian_adjoint[upper_index(0, 1)],
                      1.0e-13),
        "K-adjoint block numerator oracle");
  CHECK(close_double(adjoint_block.denominator, denominator_sum, 1.0e-13),
        "K-adjoint block denominator oracle");

  CHECK(mvmc_krylov_matrix_measurement_diagnostics_summary(
            &accumulator, 1.0e-12, 1.0, &summary) ==
            MVMC_KRYLOV_STATUS_OK &&
            summary.valid && summary.denominator_stable &&
            summary.hamiltonian_residual_available,
        "diagnostics summary");
  CHECK(close_double(summary.denominator_sum, denominator_sum, 1.0e-13),
        "denominator sum summary");
  CHECK(close_double(summary.denominator_mean,
                     denominator_sum / (double)SAMPLE_COUNT, 1.0e-13),
        "denominator mean summary");
  CHECK(close_double(summary.effective_sample_count,
                     denominator_sum * denominator_sum /
                         denominator_squared_sum,
                     1.0e-13),
        "ESS summary");
  CHECK(close_double(summary.effective_sample_fraction,
                     denominator_sum * denominator_sum /
                         (denominator_squared_sum * (double)SAMPLE_COUNT),
                     1.0e-13),
        "ESS fraction summary");
  CHECK(close_double(summary.minimum_denominator, minimum_denominator,
                     1.0e-13),
        "minimum denominator summary");
  CHECK(close_double(summary.maximum_denominator, maximum_denominator,
                     1.0e-13),
        "maximum denominator summary");
  CHECK(close_double(summary.denominator_tail_ratio,
                     maximum_denominator / summary.denominator_mean,
                     1.0e-13),
        "tail ratio summary");
  CHECK(summary.log_contribution_span > 0.0,
        "log contribution span not recorded");
  CHECK(close_double(summary.hamiltonian_norm, expected_norm, 1.0e-13),
        "Hamiltonian norm summary");
  CHECK(close_double(summary.hamiltonian_antihermitian_residual,
                     expected_residual, 1.0e-13),
        "anti-Hermitian residual summary");
  if (summary.denominator_relative_se > 0.0) {
    CHECK(mvmc_krylov_matrix_measurement_diagnostics_summary(
              &accumulator, 1.0e-12,
              0.5 * summary.denominator_relative_se, &summary) ==
              MVMC_KRYLOV_STATUS_OK &&
              summary.valid && !summary.denominator_stable,
          "denominator stability failure");
  }
}

static void test_single_sample_denominator_is_not_stable(void) {
  const MVMCKrylovMatrixMeasurementPolicy policy = policy_fixture();
  double complex raw[SAMPLE_COUNT][4];
  MVMCScaledComplex scaled[4];
  MVMCKrylovStreamingComplexSum overlap_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_adjoint_entries[UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_squared_entries[UPPER_COUNT];
  MVMCKrylovMatrixMeasurementAccumulator accumulator;
  MVMCKrylovMatrixMeasurementDiagnosticsSummary summary;

  fill_raw_dimensionless(raw);
  make_scaled_values(&policy, raw[0], scaled);
  CHECK(mvmc_krylov_matrix_measurement_accumulator_init_with_diagnostics(
            ORDER, overlap_entries, hamiltonian_entries,
            hamiltonian_adjoint_entries, hamiltonian_squared_entries,
            UPPER_COUNT, &accumulator) == MVMC_KRYLOV_STATUS_OK,
        "single-sample diagnostics accumulator init");
  CHECK(mvmc_krylov_matrix_measurement_accumulator_add_sample(
            &accumulator, &policy, scaled, 4, NULL) ==
            MVMC_KRYLOV_STATUS_OK,
        "single-sample diagnostics add");
  CHECK(mvmc_krylov_matrix_measurement_diagnostics_summary(
            &accumulator, 0.0, 0.0, &summary) ==
            MVMC_KRYLOV_STATUS_OK &&
            summary.valid && summary.sample_count == 1 &&
            summary.denominator_relative_se == 0.0 &&
            !summary.denominator_stable,
        "single sample cannot establish denominator stability");
}

static void test_block_accumulator_invalid_inputs(void) {
  const MVMCKrylovMatrixMeasurementPolicy policy = policy_fixture();
  enum { BLOCK_COUNT = 2, BLOCK_LENGTH = 2 };
  MVMCKrylovMatrixMeasurementAccumulator block_accumulators[BLOCK_COUNT];
  MVMCKrylovStreamingComplexSum overlap_entries[BLOCK_COUNT * UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_entries[
      BLOCK_COUNT * UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_adjoint_entries[
      BLOCK_COUNT * UPPER_COUNT];
  MVMCKrylovStreamingComplexSum hamiltonian_squared_entries[
      BLOCK_COUNT * UPPER_COUNT];
  MVMCKrylovMatrixMeasurementBlockAccumulator accumulator;
  MVMCKrylovMatrixMeasurementAccumulator single;
  MVMCKrylovJackknifeBlock extracted[BLOCK_COUNT];
  double complex raw[SAMPLE_COUNT][4];
  size_t sample;

  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_init(
            ORDER, 0, BLOCK_LENGTH, block_accumulators, overlap_entries,
            hamiltonian_entries, hamiltonian_adjoint_entries,
            hamiltonian_squared_entries, BLOCK_COUNT * UPPER_COUNT,
            &accumulator) == MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
        "zero block count rejected");
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_init(
            ORDER, BLOCK_COUNT, BLOCK_LENGTH, NULL, overlap_entries,
            hamiltonian_entries, hamiltonian_adjoint_entries,
            hamiltonian_squared_entries, BLOCK_COUNT * UPPER_COUNT,
            &accumulator) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "null block storage rejected");
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_init(
            ORDER, BLOCK_COUNT, BLOCK_LENGTH, block_accumulators,
            overlap_entries, hamiltonian_entries,
            hamiltonian_adjoint_entries, hamiltonian_squared_entries,
            BLOCK_COUNT * UPPER_COUNT - 1, &accumulator) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "block storage count mismatch rejected");
  CHECK(mvmc_krylov_matrix_measurement_accumulator_init(
            ORDER, overlap_entries, overlap_entries,
            hamiltonian_squared_entries, UPPER_COUNT, &single) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "single accumulator aliased storage rejected");
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_init(
            ORDER, BLOCK_COUNT, BLOCK_LENGTH, block_accumulators,
            overlap_entries, overlap_entries,
            hamiltonian_adjoint_entries, hamiltonian_squared_entries,
            BLOCK_COUNT * UPPER_COUNT, &accumulator) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "block accumulator aliased required storage rejected");
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_init(
            ORDER, BLOCK_COUNT, BLOCK_LENGTH, block_accumulators,
            overlap_entries, hamiltonian_entries, overlap_entries,
            hamiltonian_squared_entries, BLOCK_COUNT * UPPER_COUNT,
            &accumulator) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "block accumulator aliased optional storage rejected");

  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_init(
            ORDER, BLOCK_COUNT, BLOCK_LENGTH, block_accumulators,
            overlap_entries, hamiltonian_entries,
            hamiltonian_adjoint_entries, hamiltonian_squared_entries,
            BLOCK_COUNT * UPPER_COUNT, &accumulator) ==
            MVMC_KRYLOV_STATUS_OK && accumulator.valid,
        "block accumulator valid setup");
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_extract_blocks(
            &accumulator, MVMC_KRYLOV_MATRIX_HAMILTONIAN, 0, 1,
            extracted, BLOCK_COUNT) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "partially filled block partition rejected");

  fill_raw_dimensionless(raw);
  for (sample = 0; sample < SAMPLE_COUNT; ++sample) {
    MVMCScaledComplex scaled[4];
    make_scaled_values(&policy, raw[sample], scaled);
    CHECK(mvmc_krylov_matrix_measurement_block_accumulator_add_sample(
              &accumulator, &policy, scaled, 4, NULL) ==
              MVMC_KRYLOV_STATUS_OK,
          "invalid-input fixture sample %zu", sample);
  }
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_extract_blocks(
            &accumulator, MVMC_KRYLOV_MATRIX_HAMILTONIAN, 0, 1,
            extracted, BLOCK_COUNT - 1) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "block extraction count mismatch rejected");
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_extract_blocks(
            &accumulator, MVMC_KRYLOV_MATRIX_HAMILTONIAN, 0, 1,
            NULL, BLOCK_COUNT) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "null block extraction output rejected");
  CHECK(mvmc_krylov_matrix_measurement_block_accumulator_extract_blocks(
            &accumulator, MVMC_KRYLOV_MATRIX_HAMILTONIAN, 0, 1,
            extracted, BLOCK_COUNT) == MVMC_KRYLOV_STATUS_OK &&
            extracted[0].sample_count == BLOCK_LENGTH &&
            extracted[1].sample_count == BLOCK_LENGTH,
        "filled block partition remains extractable");
}

static void test_invalid_inputs(void) {
  MVMCKrylovMatrixMeasurementPolicy policy = policy_fixture();
  MVMCScaledComplex scaled[4];
  double complex raw[SAMPLE_COUNT][4];
  double complex overlap[UPPER_COUNT];
  double complex hamiltonian[UPPER_COUNT];
  double complex hamiltonian_squared[UPPER_COUNT];
  MVMCKrylovMatrixMeasurementSampleDiagnostics diagnostics;

  fill_raw_dimensionless(raw);
  make_scaled_values(&policy, raw[0], scaled);
  policy.order = 0;
  CHECK(mvmc_krylov_matrix_measurement_sample(
            &policy, scaled, 4, overlap, hamiltonian,
            hamiltonian_squared, UPPER_COUNT, &diagnostics) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "invalid order rejected");
  policy = policy_fixture();
  policy.target_weight[0] = 0.0;
  policy.target_weight[1] = 0.0;
  policy.target_weight[2] = 0.0;
  policy.target_weight[3] = 0.0;
  CHECK(mvmc_krylov_matrix_measurement_sample(
            &policy, scaled, 4, overlap, hamiltonian,
            hamiltonian_squared, UPPER_COUNT, &diagnostics) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "zero target rejected");
  policy = policy_fixture();
  scaled[1].state = MVMC_SCALED_COMPLEX_NONFINITE;
  CHECK(mvmc_krylov_matrix_measurement_sample(
            &policy, scaled, 4, overlap, hamiltonian,
            hamiltonian_squared, UPPER_COUNT, &diagnostics) ==
            MVMC_KRYLOV_STATUS_NONFINITE,
        "nonfinite scaled value rejected");
  {
    double complex target_mismatch_raw[4] = {
        0.0, 1.0 + 0.25 * I, 0.0, 0.0};
    policy = one_sided_policy_fixture();
    make_scaled_values(&policy, target_mismatch_raw, scaled);
    CHECK(mvmc_krylov_matrix_measurement_sample(
              &policy, scaled, 4, overlap, hamiltonian,
              hamiltonian_squared, UPPER_COUNT, &diagnostics) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "zero target with nonzero measurement support rejected");
  }
}

static void test_uniform_floor_only_guide(void) {
  MVMCKrylovMatrixMeasurementPolicy policy = policy_fixture();
  MVMCScaledComplex scaled[4];
  double complex raw[SAMPLE_COUNT][4];
  double complex overlap[UPPER_COUNT];
  double complex hamiltonian[UPPER_COUNT];
  double complex hamiltonian_squared[UPPER_COUNT];
  MVMCKrylovMatrixMeasurementSampleDiagnostics diagnostics;
  double expected_denominator = 0.0;
  int index;

  fill_raw_dimensionless(raw);
  policy.eta = 2.0;
  for (index = 0; index <= ORDER; ++index) {
    policy.guide_lambda[index] = 0.0;
    expected_denominator +=
        policy.target_weight[index] *
        creal(conj(raw[0][index]) * raw[0][index]) / policy.eta;
  }
  make_scaled_values(&policy, raw[0], scaled);
  CHECK(mvmc_krylov_matrix_measurement_sample(
            &policy, scaled, 4, overlap, hamiltonian,
            hamiltonian_squared, UPPER_COUNT, &diagnostics) ==
            MVMC_KRYLOV_STATUS_OK &&
            diagnostics.valid,
        "uniform floor-only guide accepted");
  CHECK(close_double(diagnostics.log_guide, log(policy.eta), 1.0e-14) &&
            close_double(diagnostics.denominator, expected_denominator,
                         1.0e-13),
        "uniform floor-only guide value");
}

int main(void) {
  test_single_sample_against_raw_oracle();
  test_accumulator_block_totals();
  test_continuous_trace_block_accumulators();
  test_one_sided_moment_oracle();
  test_offdiagonal_antihermitian_frobenius();
  test_zero_target_samples_are_accepted();
  test_zero_target_block_jackknife();
  test_diagnostics_summary();
  test_single_sample_denominator_is_not_stable();
  test_block_accumulator_invalid_inputs();
  test_uniform_floor_only_guide();
  test_invalid_inputs();
  return failures == 0 ? 0 : 1;
}
