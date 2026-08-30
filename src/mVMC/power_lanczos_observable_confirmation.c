#include "power_lanczos_observable_confirmation.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "power_lanczos_observable_confirmation.c requires the bounded engine"
#endif

#include "power_lanczos_observable_census.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum {
  CONFIRMATION_S = 0,
  CONFIRMATION_K_FORWARD,
  CONFIRMATION_K_REVERSE,
  CONFIRMATION_B
};

struct MVMCPowerLanczosObservableConfirmationSession {
  MVMCPowerLanczosObservableConfirmationStatus status;
  MVMCPowerLanczosObservableConfirmationState state;
  size_t request_count;
  size_t block_count;
  size_t coefficient_blocks_added;
  size_t final_blocks_added;
  uint64_t coefficient_block_length;
  uint64_t final_block_length;
  uint64_t coefficient_sample_count;
  uint64_t final_sample_count;
  size_t observable_matrix_stride;
  size_t allocated_bytes;
  MVMCPowerLanczosGEVPPolicy gevp_policy;
  MVMCPowerLanczosGEVPResult full_gevp;
  uint64_t *coefficient_counts;
  uint64_t *final_counts;
  double complex *matrix_blocks;
  double complex *observable_matrix_blocks;
  double complex *final_blocks;
  double complex *leave_one_coefficients;
  double *leave_one_projective_distances;
};

static int CheckedAddSize(size_t *total, size_t addition) {
  if (total == NULL || addition > SIZE_MAX - *total) return 0;
  *total += addition;
  return 1;
}

static int CheckedMultiplySize(size_t left, size_t right, size_t *product) {
  if (product == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
  *product = left * right;
  return 1;
}

static int CheckedAddU64(uint64_t *total, uint64_t addition) {
  if (total == NULL || addition > UINT64_MAX - *total) return 0;
  *total += addition;
  return 1;
}

static int FiniteComplex(double complex value) {
  return isfinite(creal(value)) && isfinite(cimag(value));
}

static int FiniteArray(const double complex *values, size_t count) {
  size_t index;
  if (values == NULL) return 0;
  for (index = 0; index < count; ++index) {
    if (!FiniteComplex(values[index])) return 0;
  }
  return 1;
}

typedef struct {
  uintptr_t begin;
  uintptr_t end;
} OutputRange;

static int MakeOutputRange(const void *pointer, size_t count,
                           size_t element_size, OutputRange *range) {
  size_t bytes;
  uintptr_t begin;
  if (pointer == NULL || range == NULL ||
      !CheckedMultiplySize(count, element_size, &bytes)) {
    return 0;
  }
  begin = (uintptr_t)pointer;
  if (bytes == 0 || bytes > UINTPTR_MAX - begin) return 0;
  range->begin = begin;
  range->end = begin + bytes;
  return 1;
}

static int OutputRangesDisjoint(const OutputRange *ranges, size_t count) {
  size_t left;
  size_t right;
  if (ranges == NULL) return 0;
  for (left = 0; left < count; ++left) {
    for (right = left + 1; right < count; ++right) {
      if (ranges[left].begin < ranges[right].end &&
          ranges[right].begin < ranges[left].end) {
        return 0;
      }
    }
  }
  return 1;
}

static MVMCPowerLanczosObservableConfirmationStatus Fail(
    MVMCPowerLanczosObservableConfirmationSession *session,
    MVMCPowerLanczosObservableConfirmationStatus status) {
  if (session != NULL) {
    session->status = status;
    session->state = MVMC_POWER_LANCZOS_CONFIRMATION_FAILED;
  }
  return status;
}

static size_t MatrixIndex(size_t block, size_t family, size_t entry) {
  return ((block * MVMC_POWER_LANCZOS_CONFIRMATION_MATRIX_FAMILY_COUNT +
           family) *
              MVMC_POWER_LANCZOS_CONFIRMATION_UPPER_COUNT +
          entry);
}

static double complex ContractFull(const double complex alpha[2],
                                   const double complex matrix[4]) {
  return conj(alpha[0]) *
             (matrix[0] * alpha[0] + matrix[1] * alpha[1]) +
         conj(alpha[1]) *
             (matrix[2] * alpha[0] + matrix[3] * alpha[1]);
}

static int Solve(const MVMCPowerLanczosObservableConfirmationSession *session,
                 size_t omitted_block,
                 MVMCPowerLanczosGEVPResult *result,
                 uint64_t *sample_count) {
  double complex matrices[
      MVMC_POWER_LANCZOS_CONFIRMATION_MATRIX_FAMILY_COUNT]
                         [MVMC_POWER_LANCZOS_CONFIRMATION_UPPER_COUNT];
  size_t block;
  size_t family;
  size_t entry;
  uint64_t count = 0;
  memset(matrices, 0, sizeof(matrices));
  for (block = 0; block < session->block_count; ++block) {
    if (block == omitted_block) continue;
    if (!CheckedAddU64(&count, session->coefficient_counts[block])) return 0;
    for (family = 0;
         family < MVMC_POWER_LANCZOS_CONFIRMATION_MATRIX_FAMILY_COUNT;
         ++family) {
      for (entry = 0;
           entry < MVMC_POWER_LANCZOS_CONFIRMATION_UPPER_COUNT; ++entry) {
        matrices[family][entry] +=
            session->matrix_blocks[MatrixIndex(block, family, entry)];
        if (!FiniteComplex(matrices[family][entry])) return 0;
      }
    }
  }
  if (count == 0) return 0;
  for (family = 0;
       family < MVMC_POWER_LANCZOS_CONFIRMATION_MATRIX_FAMILY_COUNT;
       ++family) {
    for (entry = 0;
         entry < MVMC_POWER_LANCZOS_CONFIRMATION_UPPER_COUNT; ++entry) {
      matrices[family][entry] /= (double)count;
    }
  }
  memset(result, 0, sizeof(*result));
  if (mvmc_power_lanczos_gevp_solve_complex_packed(
          &session->gevp_policy, MVMC_POWER_LANCZOS_CONFIRMATION_DIMENSION,
          matrices[CONFIRMATION_S], matrices[CONFIRMATION_K_FORWARD],
          matrices[CONFIRMATION_K_REVERSE], matrices[CONFIRMATION_B],
          MVMC_POWER_LANCZOS_CONFIRMATION_UPPER_COUNT,
          result) != MVMC_KRYLOV_GEVP_OK ||
      !result->valid || !result->observables_valid ||
      !FiniteComplex(result->coefficient[0]) ||
      !FiniteComplex(result->coefficient[1]) || !isfinite(result->energy)) {
    return 0;
  }
  if (sample_count != NULL) *sample_count = count;
  return 1;
}

static int ObservableMean(
    const MVMCPowerLanczosObservableConfirmationSession *session,
    size_t omitted_block, size_t request, double complex matrix[4]) {
  uint64_t count = 0;
  size_t block;
  size_t entry;
  memset(matrix, 0, 4 * sizeof(matrix[0]));
  for (block = 0; block < session->block_count; ++block) {
    const double complex *source;
    if (block == omitted_block) continue;
    if (!CheckedAddU64(&count, session->coefficient_counts[block])) return 0;
    source = session->observable_matrix_blocks +
             block * session->observable_matrix_stride + 4 * request;
    for (entry = 0; entry < 4; ++entry) {
      matrix[entry] += source[entry];
      if (!FiniteComplex(matrix[entry])) return 0;
    }
  }
  if (count == 0) return 0;
  for (entry = 0; entry < 4; ++entry) matrix[entry] /= (double)count;
  return 1;
}

static int ProjectiveDistance(
    const MVMCPowerLanczosObservableConfirmationSession *session,
    const double complex left[2], const double complex right[2],
    double *distance) {
  double complex overlap_upper[3] = {0.0, 0.0, 0.0};
  double complex slr;
  double complex sll;
  double complex srr;
  double value;
  uint64_t count = 0;
  size_t block;
  size_t entry;
  if (distance == NULL) return 0;
  for (block = 0; block < session->block_count; ++block) {
    if (!CheckedAddU64(&count, session->coefficient_counts[block])) return 0;
    for (entry = 0; entry < 3; ++entry) {
      overlap_upper[entry] +=
          session->matrix_blocks[MatrixIndex(block, CONFIRMATION_S, entry)];
    }
  }
  if (count == 0) return 0;
  for (entry = 0; entry < 3; ++entry) overlap_upper[entry] /= (double)count;
  slr = conj(left[0]) *
            (overlap_upper[0] * right[0] +
             overlap_upper[1] * right[1]) +
        conj(left[1]) *
            (conj(overlap_upper[1]) * right[0] +
             overlap_upper[2] * right[1]);
  sll = conj(left[0]) *
            (overlap_upper[0] * left[0] + overlap_upper[1] * left[1]) +
        conj(left[1]) *
            (conj(overlap_upper[1]) * left[0] +
             overlap_upper[2] * left[1]);
  srr = conj(right[0]) *
            (overlap_upper[0] * right[0] +
             overlap_upper[1] * right[1]) +
        conj(right[1]) *
            (conj(overlap_upper[1]) * right[0] +
             overlap_upper[2] * right[1]);
  if (!FiniteComplex(slr) || !FiniteComplex(sll) || !FiniteComplex(srr) ||
      fabs(cimag(sll)) > 1.0e-10 * fmax(1.0, fabs(creal(sll))) ||
      fabs(cimag(srr)) > 1.0e-10 * fmax(1.0, fabs(creal(srr))) ||
      creal(sll) <= 0.0 || creal(srr) <= 0.0) {
    return 0;
  }
  value = 1.0 - (cabs(slr) * cabs(slr)) / (creal(sll) * creal(srr));
  if (value < 0.0 && value > -1.0e-12) value = 0.0;
  if (!isfinite(value) || value < 0.0 || value > 1.0 + 1.0e-10) return 0;
  if (value > 1.0) value = 1.0;
  *distance = sqrt(value);
  return isfinite(*distance);
}

MVMCPowerLanczosObservableConfirmationStatus
mvmc_power_lanczos_observable_confirmation_create(
    size_t request_count, size_t block_count,
    const MVMCPowerLanczosGEVPPolicy *gevp_policy,
    MVMCPowerLanczosObservableConfirmationSession **session) {
  MVMCPowerLanczosObservableConfirmationSession *candidate;
  size_t matrix_count;
  size_t observable_count;
  size_t final_count;
  size_t bytes;
  if (session == NULL || *session != NULL || gevp_policy == NULL ||
      !gevp_policy->valid || request_count == 0 ||
      request_count > MVMC_POWER_LANCZOS_OBSERVABLE_MAX_RECORDS ||
      block_count < MVMC_POWER_LANCZOS_CONFIRMATION_MIN_BLOCK_COUNT) {
    return MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_ARGUMENT;
  }
  if (!CheckedMultiplySize(
          block_count,
          MVMC_POWER_LANCZOS_CONFIRMATION_MATRIX_FAMILY_COUNT *
              MVMC_POWER_LANCZOS_CONFIRMATION_UPPER_COUNT,
          &matrix_count) ||
      !CheckedMultiplySize(request_count, (size_t)4, &observable_count) ||
      !CheckedMultiplySize(block_count, observable_count,
                           &observable_count) ||
      !CheckedMultiplySize(block_count, request_count, &final_count)) {
    return MVMC_POWER_LANCZOS_CONFIRMATION_RESOURCE_LIMIT;
  }
  candidate = (MVMCPowerLanczosObservableConfirmationSession *)calloc(
      1, sizeof(*candidate));
  if (candidate == NULL) {
    return MVMC_POWER_LANCZOS_CONFIRMATION_ALLOCATION_FAILURE;
  }
  candidate->status = MVMC_POWER_LANCZOS_CONFIRMATION_OK;
  candidate->state = MVMC_POWER_LANCZOS_CONFIRMATION_CONFIGURED;
  candidate->request_count = request_count;
  candidate->block_count = block_count;
  candidate->observable_matrix_stride = 4 * request_count;
  candidate->gevp_policy = *gevp_policy;
  candidate->allocated_bytes = sizeof(*candidate);
#define ALLOCATE_FIELD(field, count, type)                                  \
  do {                                                                      \
    if (!CheckedMultiplySize((count), sizeof(type), &bytes) ||              \
        !CheckedAddSize(&candidate->allocated_bytes, bytes)) {              \
      mvmc_power_lanczos_observable_confirmation_destroy(candidate);        \
      return MVMC_POWER_LANCZOS_CONFIRMATION_RESOURCE_LIMIT;                \
    }                                                                       \
    candidate->field = (type *)calloc((count), sizeof(type));               \
    if (candidate->field == NULL) {                                         \
      mvmc_power_lanczos_observable_confirmation_destroy(candidate);        \
      return MVMC_POWER_LANCZOS_CONFIRMATION_ALLOCATION_FAILURE;            \
    }                                                                       \
  } while (0)
  ALLOCATE_FIELD(coefficient_counts, block_count, uint64_t);
  ALLOCATE_FIELD(final_counts, block_count, uint64_t);
  ALLOCATE_FIELD(matrix_blocks, matrix_count, double complex);
  ALLOCATE_FIELD(observable_matrix_blocks, observable_count, double complex);
  ALLOCATE_FIELD(final_blocks, final_count, double complex);
  ALLOCATE_FIELD(leave_one_coefficients, 2 * block_count, double complex);
  ALLOCATE_FIELD(leave_one_projective_distances, block_count, double);
#undef ALLOCATE_FIELD
  candidate->state = MVMC_POWER_LANCZOS_CONFIRMATION_COEFFICIENT_SAMPLING;
  *session = candidate;
  return MVMC_POWER_LANCZOS_CONFIRMATION_OK;
}

void mvmc_power_lanczos_observable_confirmation_destroy(
    MVMCPowerLanczosObservableConfirmationSession *session) {
  if (session == NULL) return;
  free(session->leave_one_projective_distances);
  free(session->leave_one_coefficients);
  free(session->final_blocks);
  free(session->observable_matrix_blocks);
  free(session->matrix_blocks);
  free(session->final_counts);
  free(session->coefficient_counts);
  memset(session, 0, sizeof(*session));
  free(session);
}

MVMCPowerLanczosObservableConfirmationStatus
mvmc_power_lanczos_observable_confirmation_add_coefficient_block(
    MVMCPowerLanczosObservableConfirmationSession *session,
    uint64_t sample_count,
    const double complex overlap_upper_sums[3],
    const double complex hamiltonian_forward_upper_sums[3],
    const double complex hamiltonian_reverse_upper_sums[3],
    const double complex hamiltonian_squared_upper_sums[3],
    const double complex *observable_matrix_sums,
    size_t observable_matrix_count) {
  const double complex *families[4] = {
      overlap_upper_sums, hamiltonian_forward_upper_sums,
      hamiltonian_reverse_upper_sums, hamiltonian_squared_upper_sums};
  const size_t block = session == NULL ? 0 : session->coefficient_blocks_added;
  uint64_t total;
  size_t family;
  if (session == NULL) {
    return MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_ARGUMENT;
  }
  if (session->state != MVMC_POWER_LANCZOS_CONFIRMATION_COEFFICIENT_SAMPLING) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_STATE);
  }
  if (block >= session->block_count) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_INTERNAL_FAILURE);
  }
  if (sample_count == 0 ||
      observable_matrix_count != session->observable_matrix_stride ||
      observable_matrix_sums == NULL) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_ARGUMENT);
  }
  if (!FiniteArray(observable_matrix_sums, observable_matrix_count)) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_NONFINITE);
  }
  if (sample_count > MVMC_POWER_LANCZOS_CONFIRMATION_MAX_EXACT_SAMPLE_COUNT) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_RESOURCE_LIMIT);
  }
  if (block != 0 && sample_count != session->coefficient_block_length) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_BLOCK_CONTRACT);
  }
  for (family = 0; family < 4; ++family) {
    if (!FiniteArray(families[family], 3)) {
      return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_NONFINITE);
    }
  }
  total = session->coefficient_sample_count;
  if (!CheckedAddU64(&total, sample_count) ||
      total > MVMC_POWER_LANCZOS_CONFIRMATION_MAX_EXACT_SAMPLE_COUNT) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_RESOURCE_LIMIT);
  }
  for (family = 0; family < 4; ++family) {
    memcpy(session->matrix_blocks + MatrixIndex(block, family, 0),
           families[family], 3 * sizeof(families[family][0]));
  }
  memcpy(session->observable_matrix_blocks +
             block * session->observable_matrix_stride,
         observable_matrix_sums,
         observable_matrix_count * sizeof(observable_matrix_sums[0]));
  session->coefficient_counts[block] = sample_count;
  if (block == 0) session->coefficient_block_length = sample_count;
  session->coefficient_sample_count = total;
  ++session->coefficient_blocks_added;
  if (session->coefficient_blocks_added == session->block_count) {
    session->state = MVMC_POWER_LANCZOS_CONFIRMATION_COEFFICIENT_READY;
  }
  return MVMC_POWER_LANCZOS_CONFIRMATION_OK;
}

MVMCPowerLanczosObservableConfirmationStatus
mvmc_power_lanczos_observable_confirmation_freeze_coefficient(
    MVMCPowerLanczosObservableConfirmationSession *session) {
  size_t block;
  if (session == NULL ||
      session->state != MVMC_POWER_LANCZOS_CONFIRMATION_COEFFICIENT_READY) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_STATE);
  }
  if (!Solve(session, SIZE_MAX, &session->full_gevp, NULL)) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_GEVP_FAILURE);
  }
  for (block = 0; block < session->block_count; ++block) {
    MVMCPowerLanczosGEVPResult leave_one;
    if (!Solve(session, block, &leave_one, NULL) ||
        !ProjectiveDistance(session, session->full_gevp.coefficient,
                            leave_one.coefficient,
                            &session->leave_one_projective_distances[block])) {
      return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_GEVP_FAILURE);
    }
    session->leave_one_coefficients[2 * block] = leave_one.coefficient[0];
    session->leave_one_coefficients[2 * block + 1] =
        leave_one.coefficient[1];
  }
  session->state = MVMC_POWER_LANCZOS_CONFIRMATION_COEFFICIENT_FROZEN;
  return MVMC_POWER_LANCZOS_CONFIRMATION_OK;
}

MVMCPowerLanczosObservableConfirmationStatus
mvmc_power_lanczos_observable_confirmation_add_final_block(
    MVMCPowerLanczosObservableConfirmationSession *session,
    uint64_t sample_count, const double complex *observable_sums,
    size_t observable_count) {
  const size_t block = session == NULL ? 0 : session->final_blocks_added;
  uint64_t total;
  if (session == NULL) {
    return MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_ARGUMENT;
  }
  if (session->state != MVMC_POWER_LANCZOS_CONFIRMATION_COEFFICIENT_FROZEN) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_STATE);
  }
  if (block >= session->block_count) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_INTERNAL_FAILURE);
  }
  if (sample_count == 0 ||
      observable_count != session->request_count || observable_sums == NULL) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_ARGUMENT);
  }
  if (!FiniteArray(observable_sums, observable_count)) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_NONFINITE);
  }
  if (sample_count > MVMC_POWER_LANCZOS_CONFIRMATION_MAX_EXACT_SAMPLE_COUNT) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_RESOURCE_LIMIT);
  }
  if (sample_count != session->coefficient_block_length ||
      (block != 0 && sample_count != session->final_block_length)) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_BLOCK_CONTRACT);
  }
  total = session->final_sample_count;
  if (!CheckedAddU64(&total, sample_count) ||
      total > MVMC_POWER_LANCZOS_CONFIRMATION_MAX_EXACT_SAMPLE_COUNT) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_RESOURCE_LIMIT);
  }
  memcpy(session->final_blocks + block * session->request_count,
         observable_sums, observable_count * sizeof(observable_sums[0]));
  session->final_counts[block] = sample_count;
  if (block == 0) session->final_block_length = sample_count;
  session->final_sample_count = total;
  ++session->final_blocks_added;
  if (session->final_blocks_added == session->block_count) {
    session->state = MVMC_POWER_LANCZOS_CONFIRMATION_FINAL_READY;
  }
  return MVMC_POWER_LANCZOS_CONFIRMATION_OK;
}

MVMCPowerLanczosObservableConfirmationStatus
mvmc_power_lanczos_observable_confirmation_finalize(
    MVMCPowerLanczosObservableConfirmationSession *session,
    double complex *coefficient_estimates, size_t coefficient_capacity,
    double complex *final_estimates, size_t final_capacity,
    double complex *leave_one_estimates, size_t leave_one_capacity,
    double complex *final_block_means, size_t final_block_capacity,
    double complex *leave_one_coefficient,
    size_t leave_one_coefficient_capacity,
    double *leave_one_projective_distance,
    size_t leave_one_projective_capacity) {
  double complex *candidate_coefficient = NULL;
  double complex *candidate_final = NULL;
  double complex *candidate_leave_one = NULL;
  double complex *candidate_final_blocks = NULL;
  size_t total_values;
  size_t coefficient_values;
  OutputRange output_ranges[6];
  size_t block;
  size_t request;
  MVMCPowerLanczosObservableConfirmationStatus status =
      MVMC_POWER_LANCZOS_CONFIRMATION_OK;
  if (session == NULL) {
    return MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_ARGUMENT;
  }
  if (session->state != MVMC_POWER_LANCZOS_CONFIRMATION_FINAL_READY) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_STATE);
  }
  if (!CheckedMultiplySize(session->block_count, session->request_count,
                           &total_values) ||
      !CheckedMultiplySize(session->block_count, (size_t)2,
                           &coefficient_values)) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_RESOURCE_LIMIT);
  }
  if (coefficient_estimates == NULL || final_estimates == NULL ||
      leave_one_estimates == NULL || final_block_means == NULL ||
      leave_one_coefficient == NULL ||
      leave_one_projective_distance == NULL ||
      coefficient_capacity < session->request_count ||
      final_capacity < session->request_count ||
      leave_one_capacity < total_values ||
      final_block_capacity < total_values ||
      leave_one_coefficient_capacity < coefficient_values ||
      leave_one_projective_capacity < session->block_count) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_ARGUMENT);
  }
  if (!MakeOutputRange(coefficient_estimates, session->request_count,
                       sizeof(coefficient_estimates[0]), &output_ranges[0]) ||
      !MakeOutputRange(final_estimates, session->request_count,
                       sizeof(final_estimates[0]), &output_ranges[1]) ||
      !MakeOutputRange(leave_one_estimates, total_values,
                       sizeof(leave_one_estimates[0]), &output_ranges[2]) ||
      !MakeOutputRange(final_block_means, total_values,
                       sizeof(final_block_means[0]), &output_ranges[3]) ||
      !MakeOutputRange(leave_one_coefficient, coefficient_values,
                       sizeof(leave_one_coefficient[0]), &output_ranges[4]) ||
      !MakeOutputRange(leave_one_projective_distance, session->block_count,
                       sizeof(leave_one_projective_distance[0]),
                       &output_ranges[5]) ||
      !OutputRangesDisjoint(output_ranges, 6)) {
    return Fail(session, MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_ARGUMENT);
  }
  candidate_coefficient =
      (double complex *)calloc(session->request_count,
                               sizeof(candidate_coefficient[0]));
  candidate_final =
      (double complex *)calloc(session->request_count,
                               sizeof(candidate_final[0]));
  candidate_leave_one =
      (double complex *)calloc(total_values, sizeof(candidate_leave_one[0]));
  candidate_final_blocks =
      (double complex *)calloc(total_values,
                               sizeof(candidate_final_blocks[0]));
  if (candidate_coefficient == NULL || candidate_final == NULL ||
      candidate_leave_one == NULL || candidate_final_blocks == NULL) {
    status = MVMC_POWER_LANCZOS_CONFIRMATION_ALLOCATION_FAILURE;
    goto cleanup;
  }
  for (request = 0; request < session->request_count; ++request) {
    double complex matrix[4];
    if (!ObservableMean(session, SIZE_MAX, request, matrix)) {
      status = MVMC_POWER_LANCZOS_CONFIRMATION_NONFINITE;
      goto cleanup;
    }
    candidate_coefficient[request] =
        ContractFull(session->full_gevp.coefficient, matrix);
  }
  for (block = 0; block < session->block_count; ++block) {
    const double complex *alpha =
        session->leave_one_coefficients + 2 * block;
    for (request = 0; request < session->request_count; ++request) {
      double complex matrix[4];
      const size_t index = block * session->request_count + request;
      if (!ObservableMean(session, block, request, matrix) ||
          session->final_counts[block] == 0) {
        status = MVMC_POWER_LANCZOS_CONFIRMATION_NONFINITE;
        goto cleanup;
      }
      candidate_leave_one[index] = ContractFull(alpha, matrix);
      candidate_final_blocks[index] =
          session->final_blocks[index] /
          (double)session->final_counts[block];
      candidate_final[request] += session->final_blocks[index];
      if (!FiniteComplex(candidate_leave_one[index]) ||
          !FiniteComplex(candidate_final_blocks[index]) ||
          !FiniteComplex(candidate_final[request])) {
        status = MVMC_POWER_LANCZOS_CONFIRMATION_NONFINITE;
        goto cleanup;
      }
    }
  }
  for (request = 0; request < session->request_count; ++request) {
    candidate_final[request] /= (double)session->final_sample_count;
    if (!FiniteComplex(candidate_final[request])) {
      status = MVMC_POWER_LANCZOS_CONFIRMATION_NONFINITE;
      goto cleanup;
    }
  }
  memcpy(coefficient_estimates, candidate_coefficient,
         session->request_count * sizeof(coefficient_estimates[0]));
  memcpy(final_estimates, candidate_final,
         session->request_count * sizeof(final_estimates[0]));
  memcpy(leave_one_estimates, candidate_leave_one,
         total_values * sizeof(leave_one_estimates[0]));
  memcpy(final_block_means, candidate_final_blocks,
         total_values * sizeof(final_block_means[0]));
  memcpy(leave_one_coefficient, session->leave_one_coefficients,
         coefficient_values * sizeof(leave_one_coefficient[0]));
  memcpy(leave_one_projective_distance,
         session->leave_one_projective_distances,
         session->block_count * sizeof(leave_one_projective_distance[0]));
  session->state = MVMC_POWER_LANCZOS_CONFIRMATION_FINALIZED;

cleanup:
  free(candidate_final_blocks);
  free(candidate_leave_one);
  free(candidate_final);
  free(candidate_coefficient);
  if (status != MVMC_POWER_LANCZOS_CONFIRMATION_OK) {
    return Fail(session, status);
  }
  return MVMC_POWER_LANCZOS_CONFIRMATION_OK;
}

MVMCPowerLanczosObservableConfirmationStatus
mvmc_power_lanczos_observable_confirmation_summary(
    const MVMCPowerLanczosObservableConfirmationSession *session,
    MVMCPowerLanczosObservableConfirmationSummary *summary) {
  MVMCPowerLanczosObservableConfirmationSummary candidate;
  size_t block;
  if (session == NULL || summary == NULL) {
    return MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_ARGUMENT;
  }
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = session->state != MVMC_POWER_LANCZOS_CONFIRMATION_FAILED;
  candidate.status = session->status;
  candidate.state = session->state;
  candidate.request_count = session->request_count;
  candidate.block_count = session->block_count;
  candidate.coefficient_blocks_added = session->coefficient_blocks_added;
  candidate.final_blocks_added = session->final_blocks_added;
  candidate.coefficient_block_length = session->coefficient_block_length;
  candidate.final_block_length = session->final_block_length;
  candidate.coefficient_sample_count = session->coefficient_sample_count;
  candidate.final_sample_count = session->final_sample_count;
  if (session->state >= MVMC_POWER_LANCZOS_CONFIRMATION_COEFFICIENT_FROZEN &&
      session->state != MVMC_POWER_LANCZOS_CONFIRMATION_FAILED) {
    candidate.coefficient[0] = session->full_gevp.coefficient[0];
    candidate.coefficient[1] = session->full_gevp.coefficient[1];
    candidate.energy = session->full_gevp.energy;
    for (block = 0; block < session->block_count; ++block) {
      if (session->leave_one_projective_distances[block] >
          candidate.maximum_leave_one_projective_distance) {
        candidate.maximum_leave_one_projective_distance =
            session->leave_one_projective_distances[block];
      }
    }
  }
  *summary = candidate;
  return MVMC_POWER_LANCZOS_CONFIRMATION_OK;
}

size_t mvmc_power_lanczos_observable_confirmation_allocated_bytes(
    const MVMCPowerLanczosObservableConfirmationSession *session) {
  return session == NULL ? 0 : session->allocated_bytes;
}

const char *mvmc_power_lanczos_observable_confirmation_status_string(
    MVMCPowerLanczosObservableConfirmationStatus status) {
  switch (status) {
    case MVMC_POWER_LANCZOS_CONFIRMATION_OK:
      return "OK";
    case MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_ARGUMENT:
      return "INVALID_ARGUMENT";
    case MVMC_POWER_LANCZOS_CONFIRMATION_INVALID_STATE:
      return "INVALID_STATE";
    case MVMC_POWER_LANCZOS_CONFIRMATION_BLOCK_CONTRACT:
      return "BLOCK_CONTRACT";
    case MVMC_POWER_LANCZOS_CONFIRMATION_RESOURCE_LIMIT:
      return "RESOURCE_LIMIT";
    case MVMC_POWER_LANCZOS_CONFIRMATION_ALLOCATION_FAILURE:
      return "ALLOCATION_FAILURE";
    case MVMC_POWER_LANCZOS_CONFIRMATION_NONFINITE:
      return "NONFINITE";
    case MVMC_POWER_LANCZOS_CONFIRMATION_GEVP_FAILURE:
      return "GEVP_FAILURE";
    case MVMC_POWER_LANCZOS_CONFIRMATION_INTERNAL_FAILURE:
      return "INTERNAL_FAILURE";
  }
  return "UNKNOWN";
}
