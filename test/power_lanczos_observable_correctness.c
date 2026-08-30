#include "power_lanczos_observable_blocks.h"
#include "power_lanczos_observable_dense_oracle.h"
#include "power_lanczos_observable_evaluator.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

#define CASE_MAX_RECORDS 64
#define CASE_MAX_OPERATORS (4 * ORACLE_MAX_TERMS)

typedef struct {
  const char *name;
  MVMCPowerLanczosObservableLayout layout;
  MVMCKrylovFermionOperator operators[CASE_MAX_OPERATORS];
  MVMCKrylovHamiltonianTerm terms[ORACLE_MAX_TERMS];
  OracleTerm oracle_terms[ORACLE_MAX_TERMS];
  size_t operator_count;
  size_t term_count;
  MVMCPowerLanczosObservableRecord records[CASE_MAX_RECORDS];
  OracleMonomial oracle_records[CASE_MAX_RECORDS];
  size_t record_count;
  double complex psi[ORACLE_MAX_FULL_DIM];
  double basis_log_scale[2];
  double complex alpha[2];
  int density_pair_record[ORACLE_MAX_ORBITALS][ORACLE_MAX_ORBITALS];
} CorrectnessCase;

typedef struct {
  const double complex *psi;
  uint64_t calls;
} AmplitudeContext;

typedef struct {
  double complex matrix[CASE_MAX_RECORDS * 4];
  double complex final[CASE_MAX_RECORDS];
  double complex oracle_matrix[CASE_MAX_RECORDS * 4];
  double complex oracle_final[CASE_MAX_RECORDS];
  uint64_t semantic_hash;
} CaseResult;

static int failures = 0;
static int world_rank = 0;

#define CHECK(condition, ...)                                           \
  do {                                                                  \
    if (!(condition)) {                                                 \
      fprintf(stderr, "ObservableCorrectness FAIL rank %d: ",          \
              world_rank);                                              \
      fprintf(stderr, __VA_ARGS__);                                     \
      fprintf(stderr, " (line %d)\n", __LINE__);                       \
      ++failures;                                                       \
    }                                                                   \
  } while (0)

static int CloseComplex(double complex actual, double complex expected) {
  const double error = cabs(actual - expected);
  const double scale = fmax(cabs(actual), cabs(expected));
  return error <= fmax(5.0e-14, 5.0e-15 * scale);
}

static uint64_t HashBytes(uint64_t hash, const void *bytes, size_t size) {
  const unsigned char *cursor = (const unsigned char *)bytes;
  size_t index;
  for (index = 0; index < size; ++index) {
    hash ^= (uint64_t)cursor[index];
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

static MVMCKrylovStatus TableAmplitude(
    const uint64_t *words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
  AmplitudeContext *amplitude = (AmplitudeContext *)context;
  MVMCPfaffianStatus pf_status;
  if (words == NULL || word_count != 1 || amplitude == NULL ||
      result == NULL || words[0] >= ORACLE_MAX_FULL_DIM) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  ++amplitude->calls;
  memset(result, 0, sizeof(*result));
  if (creal(amplitude->psi[words[0]]) == 0.0 &&
      cimag(amplitude->psi[words[0]]) == 0.0) {
    pf_status = mvmc_scaled_complex_make_exact_zero(&result->value);
    result->exact_zero_component_count = 1;
  } else {
    pf_status = mvmc_scaled_complex_from_raw(amplitude->psi[words[0]],
                                             &result->value);
    result->well_pivoted_component_count = 1;
  }
  if (pf_status != MVMC_PFAFFIAN_STATUS_OK) {
    return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
  }
  result->local_factorization_count = 1;
  result->global_factorization_count = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static int AddTerm(CorrectnessCase *test_case, double complex coefficient,
                   size_t operator_count, const size_t orbital[4]) {
  MVMCKrylovHamiltonianTerm *term;
  OracleTerm *oracle_term;
  size_t index;
  if (test_case == NULL ||
      (operator_count != 2 && operator_count != 4) ||
      test_case->term_count >= ORACLE_MAX_TERMS ||
      test_case->operator_count + operator_count > CASE_MAX_OPERATORS) {
    return 0;
  }
  term = &test_case->terms[test_case->term_count];
  oracle_term = &test_case->oracle_terms[test_case->term_count];
  memset(term, 0, sizeof(*term));
  memset(oracle_term, 0, sizeof(*oracle_term));
  term->coefficient = coefficient;
  term->operator_offset = test_case->operator_count;
  term->operator_count = operator_count;
  term->source_kind = 7;
  term->source_index = test_case->term_count;
  oracle_term->coefficient = coefficient;
  oracle_term->monomial.operator_count = operator_count;
  for (index = 0; index < operator_count; ++index) {
    test_case->operators[test_case->operator_count + index].kind =
        (index % 2) == 0 ? MVMC_KRYLOV_FERMION_CREATE
                         : MVMC_KRYLOV_FERMION_ANNIHILATE;
    test_case->operators[test_case->operator_count + index].orbital =
        orbital[index];
    oracle_term->monomial.orbital[index] = orbital[index];
  }
  test_case->operator_count += operator_count;
  ++test_case->term_count;
  return 1;
}

static int AddRecord(CorrectnessCase *test_case,
                     MVMCPowerLanczosObservableFamily family,
                     size_t operator_count, const size_t orbital[4]) {
  MVMCPowerLanczosObservableRecord *record;
  OracleMonomial *oracle;
  size_t index;
  if (test_case == NULL || test_case->record_count >= CASE_MAX_RECORDS ||
      (family == MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS &&
       operator_count != 2) ||
      (family != MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS &&
       operator_count != 4)) {
    return 0;
  }
  record = &test_case->records[test_case->record_count];
  oracle = &test_case->oracle_records[test_case->record_count];
  memset(record, 0, sizeof(*record));
  memset(oracle, 0, sizeof(*oracle));
  record->family = family;
  record->row_width = (int)(2 * operator_count);
  oracle->operator_count = operator_count;
  for (index = 0; index < operator_count; ++index) {
    record->raw_indices[2 * index] =
        (int)(orbital[index] % (size_t)test_case->layout.nsite);
    record->raw_indices[2 * index + 1] =
        (int)(orbital[index] / (size_t)test_case->layout.nsite);
    oracle->orbital[index] = orbital[index];
  }
  ++test_case->record_count;
  return 1;
}

static MVMCKrylovBoundedLimits Limits(size_t term_count) {
  MVMCKrylovBoundedLimits limits;
  memset(&limits, 0, sizeof(limits));
  limits.amplitude_policy_hash = UINT64_C(0x503643324f524143);
  limits.cache_bytes = 4 * 1024 * 1024;
  limits.max_row_transitions = term_count + 8;
  limits.max_workspace_bytes = 128 * 1024 * 1024;
  limits.max_node_expansions = UINT64_C(100000);
  limits.max_terminal_amplitude_calls = UINT64_C(100000);
  limits.max_total_row_transitions = UINT64_C(1000000);
  limits.max_order = 1;
  return limits;
}

static int RunCase(const CorrectnessCase *test_case, uint64_t generation,
                   CaseResult *case_result) {
  OracleLayout oracle_layout;
  OracleBasis basis;
  double complex hamiltonian[ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM];
  double complex hpsi[ORACLE_MAX_FULL_DIM];
  double complex scaled0[ORACLE_MAX_FULL_DIM];
  double complex scaled1[ORACLE_MAX_FULL_DIM];
  double complex coefficient_sample[CASE_MAX_RECORDS * 4];
  double complex final_sample[CASE_MAX_RECORDS];
  double complex coefficient_blocks[CASE_MAX_RECORDS * 4];
  double complex final_blocks[CASE_MAX_RECORDS];
  uint64_t coefficient_counts[1];
  uint64_t final_counts[1];
  MVMCKrylovFockModel model;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovBoundedPlan *bounded_plan = NULL;
  MVMCKrylovBoundedWorkspace *bounded_workspace = NULL;
  MVMCPowerLanczosObservableEvaluatorWorkspace *evaluator = NULL;
  MVMCPowerLanczosObservableBlockAccumulator *coefficient_accumulator = NULL;
  MVMCPowerLanczosObservableBlockAccumulator *final_accumulator = NULL;
  MVMCPowerLanczosObservablePlan observable_plan;
  MVMCPowerLanczosObservableSampleDiagnostics diagnostics;
  AmplitudeContext amplitude;
  MVMCKrylovStatus status;
  size_t sector_index;
  size_t record;
  double actual_log_guide;
  int ok = 1;
  memset(case_result, 0, sizeof(*case_result));
  oracle_layout.nsite = test_case->layout.nsite;
  oracle_layout.up_electron_count = test_case->layout.up_electron_count;
  oracle_layout.down_electron_count = test_case->layout.down_electron_count;
  oracle_layout.pure_spin = test_case->layout.pure_spin;
  if (!oracle_build_basis(&oracle_layout, &basis) ||
      !oracle_build_hamiltonian_matrix(
          2 * test_case->layout.nsite, test_case->oracle_terms,
          test_case->term_count, hamiltonian)) {
    CHECK(0, "%s oracle setup", test_case->name);
    return 0;
  }
  actual_log_guide = -log((double)basis.sector_dimension);
  oracle_matrix_vector(hamiltonian, basis.full_dimension, test_case->psi,
                       hpsi);
  if (strcmp(test_case->name, "electronic_complex_zero_support") == 0) {
    double complex bridge_operator[
        ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM];
    CHECK(oracle_build_monomial_matrix(
              2 * test_case->layout.nsite,
              &test_case->oracle_records[1], bridge_operator) &&
              cabs(test_case->psi[13]) == 0.0 && cabs(hpsi[13]) > 1.0e-12 &&
              cabs(test_case->psi[14]) > 1.0e-12 &&
              cabs(bridge_operator[13 * basis.full_dimension + 14]) == 1.0,
          "%s explicit zero-support bridge", test_case->name);
  }
  if (strcmp(test_case->name, "pure_spin_transverse") == 0) {
    double complex bridge_operator[
        ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM];
    CHECK(oracle_build_monomial_matrix(
              2 * test_case->layout.nsite,
              &test_case->oracle_records[0], bridge_operator) &&
              cabs(test_case->psi[9]) == 0.0 && cabs(hpsi[9]) > 1.0e-12 &&
              cabs(test_case->psi[6]) > 1.0e-12 &&
              cabs(bridge_operator[9 * basis.full_dimension + 6]) == 1.0,
          "%s pure-spin zero-support bridge", test_case->name);
  }
  for (sector_index = 0; sector_index < basis.full_dimension;
       ++sector_index) {
    scaled0[sector_index] =
        exp(test_case->basis_log_scale[0]) * test_case->psi[sector_index];
    scaled1[sector_index] =
        exp(test_case->basis_log_scale[1]) * hpsi[sector_index];
  }

  memset(&model, 0, sizeof(model));
  model.site_count = (size_t)test_case->layout.nsite;
  model.up_electron_count =
      (size_t)test_case->layout.up_electron_count;
  model.down_electron_count =
      (size_t)test_case->layout.down_electron_count;
  model.pure_spin = test_case->layout.pure_spin;
  model.hermitian = 1;
  model.terms = test_case->terms;
  model.term_count = test_case->term_count;
  model.operators = test_case->operators;
  model.operator_count = test_case->operator_count;
  limits = Limits(test_case->term_count);
  memset(&observable_plan, 0, sizeof(observable_plan));
  observable_plan.nsite = test_case->layout.nsite;
  observable_plan.record_count = (int)test_case->record_count;
  observable_plan.records =
      (MVMCPowerLanczosObservableRecord *)test_case->records;
  amplitude.psi = test_case->psi;
  amplitude.calls = 0;

  status = mvmc_bounded_krylov_plan_create(&model, &limits, &bounded_plan);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_workspace_create(
        bounded_plan, &bounded_workspace);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_evaluator_workspace_create(
        test_case->record_count, 1, &evaluator);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_block_create(
        MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT,
        test_case->record_count, 1, (uint64_t)basis.sector_dimension,
        &coefficient_accumulator);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_block_create(
        MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_FINAL,
        test_case->record_count, 1, (uint64_t)basis.sector_dimension,
        &final_accumulator);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    CHECK(0, "%s production workspace setup status=%d", test_case->name,
          (int)status);
    ok = 0;
    goto cleanup;
  }

  status = mvmc_bounded_krylov_session_begin(
      bounded_workspace, TableAmplitude, &amplitude, generation);
  for (sector_index = 0;
       status == MVMC_KRYLOV_STATUS_OK &&
       sector_index < basis.sector_dimension;
       ++sector_index) {
    const uint64_t source = basis.sector_words[sector_index];
    status = mvmc_power_lanczos_observable_coefficient_sample(
        evaluator, bounded_workspace, &test_case->layout, &observable_plan,
        &source, 1, actual_log_guide, test_case->basis_log_scale,
        coefficient_sample, test_case->record_count * 4, &diagnostics);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      CHECK(diagnostics.engine_root_evaluations <=
                1 + diagnostics.unique_target_count,
            "%s target grouping evaluation bound", test_case->name);
      status = mvmc_power_lanczos_observable_block_add_sample(
          coefficient_accumulator, coefficient_sample,
          test_case->record_count * 4);
    }
  }
  if (mvmc_bounded_krylov_session_is_active(bounded_workspace)) {
    const MVMCKrylovStatus end_status =
        mvmc_bounded_krylov_session_end(bounded_workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) status = end_status;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_block_export(
        coefficient_accumulator, coefficient_blocks,
        test_case->record_count * 4, coefficient_counts, 1);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    CHECK(0, "%s coefficient chain status=%d", test_case->name,
          (int)status);
    ok = 0;
    goto cleanup;
  }

  status = mvmc_bounded_krylov_session_begin(
      bounded_workspace, TableAmplitude, &amplitude,
      generation ^ UINT64_C(0x8000000000000000));
  for (sector_index = 0;
       status == MVMC_KRYLOV_STATUS_OK &&
       sector_index < basis.sector_dimension;
       ++sector_index) {
    const uint64_t source = basis.sector_words[sector_index];
    status = mvmc_power_lanczos_observable_final_sample(
        evaluator, bounded_workspace, &test_case->layout, &observable_plan,
        &source, 1, actual_log_guide, test_case->basis_log_scale,
        test_case->alpha, final_sample, test_case->record_count,
        &diagnostics);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_power_lanczos_observable_block_add_sample(
          final_accumulator, final_sample, test_case->record_count);
    }
  }
  if (mvmc_bounded_krylov_session_is_active(bounded_workspace)) {
    const MVMCKrylovStatus end_status =
        mvmc_bounded_krylov_session_end(bounded_workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) status = end_status;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_block_export(
        final_accumulator, final_blocks, test_case->record_count,
        final_counts, 1);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) {
    CHECK(0, "%s final chain status=%d", test_case->name, (int)status);
    ok = 0;
    goto cleanup;
  }
  CHECK(coefficient_counts[0] == basis.sector_dimension &&
            final_counts[0] == basis.sector_dimension,
        "%s block sample count", test_case->name);

  for (record = 0; record < test_case->record_count; ++record) {
    double complex oracle_matrix[4];
    double complex oracle_final;
    double complex quadratic = 0.0;
    size_t entry;
    CHECK(oracle_observable_matrix(
              2 * test_case->layout.nsite,
              &test_case->oracle_records[record], scaled0, scaled1,
              oracle_matrix),
          "%s oracle matrix record %zu", test_case->name, record);
    CHECK(oracle_observable_final(
              2 * test_case->layout.nsite,
              &test_case->oracle_records[record], scaled0, scaled1,
              test_case->alpha, &oracle_final),
          "%s oracle final record %zu", test_case->name, record);
    for (entry = 0; entry < 4; ++entry) {
      const double complex production =
          coefficient_blocks[4 * record + entry] /
          (double)coefficient_counts[0];
      case_result->matrix[4 * record + entry] = production;
      case_result->oracle_matrix[4 * record + entry] =
          oracle_matrix[entry];
      CHECK(CloseComplex(production, oracle_matrix[entry]),
            "%s matrix record %zu entry %zu actual=(%.17g,%.17g) "
            "oracle=(%.17g,%.17g)",
            test_case->name, record, entry, creal(production),
            cimag(production), creal(oracle_matrix[entry]),
            cimag(oracle_matrix[entry]));
    }
    quadratic = conj(test_case->alpha[0]) *
                    (oracle_matrix[0] * test_case->alpha[0] +
                     oracle_matrix[1] * test_case->alpha[1]) +
                conj(test_case->alpha[1]) *
                    (oracle_matrix[2] * test_case->alpha[0] +
                     oracle_matrix[3] * test_case->alpha[1]);
    case_result->final[record] =
        final_blocks[record] / (double)final_counts[0];
    case_result->oracle_final[record] = oracle_final;
    CHECK(CloseComplex(case_result->final[record], oracle_final),
          "%s direct final record %zu", test_case->name, record);
    CHECK(CloseComplex(quadratic, oracle_final),
          "%s oracle quadratic/direct record %zu", test_case->name,
          record);
  }
  case_result->semantic_hash = UINT64_C(1469598103934665603);
  case_result->semantic_hash = HashBytes(
      case_result->semantic_hash, case_result->matrix,
      test_case->record_count * 4 * sizeof(case_result->matrix[0]));
  case_result->semantic_hash = HashBytes(
      case_result->semantic_hash, case_result->final,
      test_case->record_count * sizeof(case_result->final[0]));

cleanup:
  mvmc_power_lanczos_observable_block_destroy(final_accumulator);
  mvmc_power_lanczos_observable_block_destroy(coefficient_accumulator);
  mvmc_power_lanczos_observable_evaluator_workspace_destroy(evaluator);
  mvmc_bounded_krylov_workspace_destroy(bounded_workspace);
  mvmc_bounded_krylov_plan_destroy(bounded_plan);
  return ok;
}

static int RunDistributedBlockCase(const CorrectnessCase *test_case,
                                   uint64_t generation, int world_size) {
  OracleLayout oracle_layout;
  OracleBasis basis;
  double complex hamiltonian[ORACLE_MAX_FULL_DIM * ORACLE_MAX_FULL_DIM];
  double complex hpsi[ORACLE_MAX_FULL_DIM];
  double complex scaled0[ORACLE_MAX_FULL_DIM];
  double complex scaled1[ORACLE_MAX_FULL_DIM];
  MVMCKrylovFockModel model;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovBoundedPlan *bounded_plan = NULL;
  MVMCKrylovBoundedWorkspace *bounded_workspace = NULL;
  MVMCPowerLanczosObservableEvaluatorWorkspace *evaluator = NULL;
  MVMCPowerLanczosObservableBlockAccumulator *local_accumulator = NULL;
  MVMCPowerLanczosObservablePlan observable_plan;
  MVMCPowerLanczosObservableSampleDiagnostics diagnostics;
  AmplitudeContext amplitude;
  MVMCKrylovStatus status;
  double complex *local_samples = NULL;
  double complex *local_block = NULL;
  double complex *gathered_samples = NULL;
  double complex *gathered_blocks = NULL;
  uint64_t local_count_wire[1] = {0};
  uint64_t *gathered_counts = NULL;
  size_t local_sector_count;
  size_t entry_count;
  size_t local_index;
  double actual_log_guide;
  int local_ready = 1;
  int collective_ready = 1;
  int ok = 1;

  oracle_layout.nsite = test_case->layout.nsite;
  oracle_layout.up_electron_count = test_case->layout.up_electron_count;
  oracle_layout.down_electron_count = test_case->layout.down_electron_count;
  oracle_layout.pure_spin = test_case->layout.pure_spin;
  if (world_size <= 0 ||
      !oracle_build_basis(&oracle_layout, &basis) ||
      basis.sector_dimension % (size_t)world_size != 0 ||
      !oracle_build_hamiltonian_matrix(
          2 * test_case->layout.nsite, test_case->oracle_terms,
          test_case->term_count, hamiltonian)) {
    CHECK(0, "%s distributed sector setup for mpi=%d", test_case->name,
          world_size);
    return 0;
  }
  local_sector_count = basis.sector_dimension / (size_t)world_size;
  entry_count = test_case->record_count * 4;
  actual_log_guide = -log((double)basis.sector_dimension);
  oracle_matrix_vector(hamiltonian, basis.full_dimension, test_case->psi,
                       hpsi);
  for (local_index = 0; local_index < basis.full_dimension; ++local_index) {
    scaled0[local_index] =
        exp(test_case->basis_log_scale[0]) * test_case->psi[local_index];
    scaled1[local_index] =
        exp(test_case->basis_log_scale[1]) * hpsi[local_index];
  }

  memset(&model, 0, sizeof(model));
  model.site_count = (size_t)test_case->layout.nsite;
  model.up_electron_count =
      (size_t)test_case->layout.up_electron_count;
  model.down_electron_count =
      (size_t)test_case->layout.down_electron_count;
  model.pure_spin = test_case->layout.pure_spin;
  model.hermitian = 1;
  model.terms = test_case->terms;
  model.term_count = test_case->term_count;
  model.operators = test_case->operators;
  model.operator_count = test_case->operator_count;
  limits = Limits(test_case->term_count);
  memset(&observable_plan, 0, sizeof(observable_plan));
  observable_plan.nsite = test_case->layout.nsite;
  observable_plan.record_count = (int)test_case->record_count;
  observable_plan.records =
      (MVMCPowerLanczosObservableRecord *)test_case->records;
  amplitude.psi = test_case->psi;
  amplitude.calls = 0;

  local_samples = (double complex *)calloc(
      local_sector_count * entry_count, sizeof(local_samples[0]));
  local_block =
      (double complex *)calloc(entry_count, sizeof(local_block[0]));
  status = mvmc_bounded_krylov_plan_create(&model, &limits, &bounded_plan);
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_workspace_create(
        bounded_plan, &bounded_workspace);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_evaluator_workspace_create(
        test_case->record_count, 1, &evaluator);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_block_create(
        MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT,
        test_case->record_count, 1, (uint64_t)local_sector_count,
        &local_accumulator);
  }
  local_ready = local_samples != NULL && local_block != NULL &&
                status == MVMC_KRYLOV_STATUS_OK;
#ifdef _mpi_use
  MPI_Allreduce(&local_ready, &collective_ready, 1, MPI_INT, MPI_MIN,
                MPI_COMM_WORLD);
#else
  collective_ready = local_ready;
#endif
  if (!collective_ready) {
    CHECK(0, "%s distributed workspace setup status=%d", test_case->name,
          (int)status);
    ok = 0;
    goto cleanup;
  }

  status = mvmc_bounded_krylov_session_begin(
      bounded_workspace, TableAmplitude, &amplitude, generation);
  for (local_index = 0;
       status == MVMC_KRYLOV_STATUS_OK &&
       local_index < local_sector_count;
       ++local_index) {
    const size_t sector_index =
        (size_t)world_rank * local_sector_count + local_index;
    const uint64_t source = basis.sector_words[sector_index];
    double complex *sample =
        local_samples + local_index * entry_count;
    status = mvmc_power_lanczos_observable_coefficient_sample(
        evaluator, bounded_workspace, &test_case->layout, &observable_plan,
        &source, 1, actual_log_guide, test_case->basis_log_scale,
        sample, entry_count, &diagnostics);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_power_lanczos_observable_block_add_sample(
          local_accumulator, sample, entry_count);
    }
  }
  if (mvmc_bounded_krylov_session_is_active(bounded_workspace)) {
    const MVMCKrylovStatus end_status =
        mvmc_bounded_krylov_session_end(bounded_workspace);
    if (status == MVMC_KRYLOV_STATUS_OK) status = end_status;
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_power_lanczos_observable_block_export(
        local_accumulator, local_block, entry_count, local_count_wire, 1);
  }
  local_ready = status == MVMC_KRYLOV_STATUS_OK &&
                local_count_wire[0] == local_sector_count;
#ifdef _mpi_use
  MPI_Allreduce(&local_ready, &collective_ready, 1, MPI_INT, MPI_MIN,
                MPI_COMM_WORLD);
#else
  collective_ready = local_ready;
#endif
  if (!collective_ready) {
    CHECK(0, "%s rank-local partial status=%d count=%llu",
          test_case->name, (int)status,
          (unsigned long long)local_count_wire[0]);
    ok = 0;
    goto cleanup;
  }

#ifdef _mpi_use
  if (world_rank == 0) {
    gathered_samples = (double complex *)calloc(
        basis.sector_dimension * entry_count,
        sizeof(gathered_samples[0]));
    gathered_blocks = (double complex *)calloc(
        (size_t)world_size * entry_count, sizeof(gathered_blocks[0]));
    gathered_counts = (uint64_t *)calloc(
        (size_t)world_size, sizeof(gathered_counts[0]));
    collective_ready = gathered_samples != NULL &&
                       gathered_blocks != NULL && gathered_counts != NULL;
  }
  MPI_Bcast(&collective_ready, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (!collective_ready) {
    CHECK(0, "%s MPI gather allocation", test_case->name);
    ok = 0;
    goto cleanup;
  }
  MPI_Gather(local_samples, (int)(local_sector_count * entry_count),
             MPI_C_DOUBLE_COMPLEX, gathered_samples,
             (int)(local_sector_count * entry_count), MPI_C_DOUBLE_COMPLEX,
             0, MPI_COMM_WORLD);
  MPI_Gather(local_block, (int)entry_count, MPI_C_DOUBLE_COMPLEX,
             gathered_blocks, (int)entry_count, MPI_C_DOUBLE_COMPLEX,
             0, MPI_COMM_WORLD);
  MPI_Gather(local_count_wire, 1, MPI_UINT64_T, gathered_counts, 1,
             MPI_UINT64_T, 0, MPI_COMM_WORLD);
#else
  gathered_samples = local_samples;
  gathered_blocks = local_block;
  gathered_counts = local_count_wire;
#endif

  if (world_rank == 0) {
    MVMCPowerLanczosObservableBlockAccumulator **rank_accumulators = NULL;
    MVMCPowerLanczosObservableBlockAccumulator *merged = NULL;
    MVMCPowerLanczosObservableBlockAccumulator *serial = NULL;
    MVMCPowerLanczosObservableBlockSummary merged_summary;
    double complex *reconstructed_block = NULL;
    double complex *merged_block = NULL;
    double complex *serial_block = NULL;
    uint64_t reconstructed_count[1] = {0};
    uint64_t merged_count[1] = {0};
    uint64_t serial_count[1] = {0};
    size_t rank_index;
    size_t record;

    memset(&merged_summary, 0, sizeof(merged_summary));
    rank_accumulators =
        (MVMCPowerLanczosObservableBlockAccumulator **)calloc(
            (size_t)world_size, sizeof(rank_accumulators[0]));
    reconstructed_block = (double complex *)calloc(
        entry_count, sizeof(reconstructed_block[0]));
    merged_block =
        (double complex *)calloc(entry_count, sizeof(merged_block[0]));
    serial_block =
        (double complex *)calloc(entry_count, sizeof(serial_block[0]));
    if (rank_accumulators == NULL || reconstructed_block == NULL ||
        merged_block == NULL || serial_block == NULL) {
      CHECK(0, "%s reducer reconstruction allocation", test_case->name);
      ok = 0;
      goto root_cleanup;
    }
    status = mvmc_power_lanczos_observable_block_create(
        MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT,
        test_case->record_count, 1, (uint64_t)local_sector_count, &merged);
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_power_lanczos_observable_block_create(
          MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT,
          test_case->record_count, 1, (uint64_t)basis.sector_dimension,
          &serial);
    }
    for (rank_index = 0;
         status == MVMC_KRYLOV_STATUS_OK &&
         rank_index < (size_t)world_size;
         ++rank_index) {
      status = mvmc_power_lanczos_observable_block_create(
          MVMC_POWER_LANCZOS_OBSERVABLE_BLOCK_COEFFICIENT,
          test_case->record_count, 1, (uint64_t)local_sector_count,
          &rank_accumulators[rank_index]);
    }
    for (rank_index = 0;
         status == MVMC_KRYLOV_STATUS_OK &&
         rank_index < (size_t)world_size;
         ++rank_index) {
      for (local_index = 0;
           status == MVMC_KRYLOV_STATUS_OK &&
           local_index < local_sector_count;
           ++local_index) {
        const double complex *sample =
            gathered_samples +
            (rank_index * local_sector_count + local_index) * entry_count;
        status = mvmc_power_lanczos_observable_block_add_sample(
            rank_accumulators[rank_index], sample, entry_count);
        if (status == MVMC_KRYLOV_STATUS_OK) {
          status = mvmc_power_lanczos_observable_block_add_sample(
              serial, sample, entry_count);
        }
      }
      if (status == MVMC_KRYLOV_STATUS_OK) {
        status = mvmc_power_lanczos_observable_block_export(
            rank_accumulators[rank_index], reconstructed_block,
            entry_count, reconstructed_count, 1);
      }
      CHECK(status == MVMC_KRYLOV_STATUS_OK &&
                reconstructed_count[0] == gathered_counts[rank_index],
            "%s gathered partial count rank %zu", test_case->name,
            rank_index);
      for (local_index = 0;
           status == MVMC_KRYLOV_STATUS_OK && local_index < entry_count;
           ++local_index) {
        CHECK(CloseComplex(
                  reconstructed_block[local_index],
                  gathered_blocks[rank_index * entry_count + local_index]),
              "%s gathered partial entry rank %zu entry %zu",
              test_case->name, rank_index, local_index);
      }
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_power_lanczos_observable_block_reduce_rank_ordered(
          (const MVMCPowerLanczosObservableBlockAccumulator *const *)
              rank_accumulators,
          (size_t)world_size, merged);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_power_lanczos_observable_block_export(
          merged, merged_block, entry_count, merged_count, 1);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_power_lanczos_observable_block_export(
          serial, serial_block, entry_count, serial_count, 1);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
      status = mvmc_power_lanczos_observable_block_summary(
          merged, &merged_summary);
    }
    CHECK(status == MVMC_KRYLOV_STATUS_OK,
          "%s rank-order reduction status=%d", test_case->name,
          (int)status);
    CHECK(merged_count[0] == basis.sector_dimension &&
              serial_count[0] == basis.sector_dimension &&
              merged_summary.completed_sample_count ==
                  basis.sector_dimension,
          "%s merged/serial sample count", test_case->name);
    for (local_index = 0; local_index < entry_count; ++local_index) {
      CHECK(CloseComplex(merged_block[local_index],
                         serial_block[local_index]),
            "%s merged/serial entry %zu mpi=%d", test_case->name,
            local_index, world_size);
    }
    for (record = 0; record < test_case->record_count; ++record) {
      double complex oracle_matrix[4];
      size_t entry;
      CHECK(oracle_observable_matrix(
                2 * test_case->layout.nsite,
                &test_case->oracle_records[record], scaled0, scaled1,
                oracle_matrix),
            "%s distributed oracle record %zu", test_case->name, record);
      for (entry = 0; entry < 4; ++entry) {
        const size_t index = 4 * record + entry;
        CHECK(CloseComplex(
                  serial_block[index] / (double)serial_count[0],
                  oracle_matrix[entry]),
              "%s distributed serial/oracle record %zu entry %zu",
              test_case->name, record, entry);
      }
    }

root_cleanup:
    if (rank_accumulators != NULL) {
      for (rank_index = 0; rank_index < (size_t)world_size; ++rank_index) {
        mvmc_power_lanczos_observable_block_destroy(
            rank_accumulators[rank_index]);
      }
    }
    free(serial_block);
    free(merged_block);
    free(reconstructed_block);
    mvmc_power_lanczos_observable_block_destroy(serial);
    mvmc_power_lanczos_observable_block_destroy(merged);
    free(rank_accumulators);
  }

cleanup:
#ifdef _mpi_use
  free(gathered_counts);
  free(gathered_blocks);
  free(gathered_samples);
#endif
  mvmc_power_lanczos_observable_block_destroy(local_accumulator);
  mvmc_power_lanczos_observable_evaluator_workspace_destroy(evaluator);
  mvmc_bounded_krylov_workspace_destroy(bounded_workspace);
  mvmc_bounded_krylov_plan_destroy(bounded_plan);
  free(local_block);
  free(local_samples);
  return ok;
}

static void BuildLowerCase(CorrectnessCase *test_case) {
  const size_t forward[4] = {0, 1, 0, 0};
  const size_t adjoint[4] = {1, 0, 0, 0};
  memset(test_case, 0, sizeof(*test_case));
  test_case->name = "lower_triangle";
  test_case->layout = (MVMCPowerLanczosObservableLayout){2, 1, 0, 0};
  CHECK(AddTerm(test_case, 1.0, 2, forward), "lower forward term");
  CHECK(AddTerm(test_case, 1.0, 2, adjoint), "lower adjoint term");
  CHECK(AddRecord(test_case, MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS,
                  2, forward),
        "lower observable");
  test_case->psi[2] = 1.0;
  test_case->alpha[0] = 1.0 / sqrt(2.0);
  test_case->alpha[1] = 1.0 / sqrt(2.0);
}

static void BuildElectronicCase(CorrectnessCase *test_case) {
  const double diagonal[6] = {0.2, -0.1, 0.35, -0.4, 0.15, 0.5};
  const double complex up_hop = 0.7 + 0.2 * I;
  const double complex down_hop = -0.3 + 0.4 * I;
  const size_t hop_up[4] = {0, 1, 0, 0};
  const size_t hop_up_adj[4] = {1, 0, 0, 0};
  const size_t hop_down[4] = {3, 5, 0, 0};
  const size_t hop_down_adj[4] = {5, 3, 0, 0};
  const size_t density_density[4] = {0, 0, 3, 3};
  const size_t density_up0[4] = {0, 0, 0, 0};
  const size_t parity_hop[4] = {0, 2, 0, 0};
  const size_t distinct[4] = {0, 2, 5, 3};
  const size_t distinct_adjoint[4] = {3, 5, 2, 0};
  const size_t repeated[4] = {0, 1, 0, 2};
  OracleLayout layout;
  OracleBasis basis;
  size_t orbital;
  size_t index;
  size_t left_orbital;
  size_t right_orbital;
  memset(test_case, 0, sizeof(*test_case));
  for (left_orbital = 0; left_orbital < ORACLE_MAX_ORBITALS;
       ++left_orbital) {
    for (right_orbital = 0; right_orbital < ORACLE_MAX_ORBITALS;
         ++right_orbital) {
      test_case->density_pair_record[left_orbital][right_orbital] = -1;
    }
  }
  test_case->name = "electronic_complex_zero_support";
  test_case->layout = (MVMCPowerLanczosObservableLayout){3, 2, 1, 0};
  for (orbital = 0; orbital < 6; ++orbital) {
    const size_t density[4] = {orbital, orbital, 0, 0};
    CHECK(AddTerm(test_case, diagonal[orbital], 2, density),
          "electronic diagonal term");
  }
  CHECK(AddTerm(test_case, up_hop, 2, hop_up), "electronic up hop");
  CHECK(AddTerm(test_case, conj(up_hop), 2, hop_up_adj),
        "electronic up hop adjoint");
  CHECK(AddTerm(test_case, down_hop, 2, hop_down),
        "electronic down hop");
  CHECK(AddTerm(test_case, conj(down_hop), 2, hop_down_adj),
        "electronic down hop adjoint");
  CHECK(AddTerm(test_case, 0.6, 4, density_density),
        "electronic interaction");
  CHECK(AddRecord(test_case, MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, 2,
                  density_up0),
        "density record");
  CHECK(AddRecord(test_case, MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, 2,
                  hop_up),
        "hop record");
  CHECK(AddRecord(test_case, MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, 2,
                  hop_up_adj),
        "hop adjoint record");
  CHECK(AddRecord(test_case, MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, 2,
                  parity_hop),
        "parity record");
  CHECK(AddRecord(test_case, MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, 2,
                  hop_down),
        "down hop record");
  CHECK(AddRecord(test_case,
                  MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX, 4,
                  density_density),
        "density-density record");
  test_case->density_pair_record[0][3] = 5;
  CHECK(AddRecord(test_case,
                  MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX, 4,
                  distinct),
        "distinct record");
  CHECK(AddRecord(test_case,
                  MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX, 4,
                  repeated),
        "repeated Pauli record");
  CHECK(AddRecord(test_case,
                  MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALT, 4,
                  distinct),
        "cross-family record");
  CHECK(AddRecord(test_case,
                  MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX, 4,
                  distinct_adjoint),
        "distinct adjoint record");
  CHECK(AddRecord(test_case,
                  MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX, 4,
                  distinct),
        "same-family cancellation partner record");
  for (left_orbital = 0; left_orbital < 6; ++left_orbital) {
    for (right_orbital = 0; right_orbital < 6; ++right_orbital) {
      const size_t pair[4] = {left_orbital, left_orbital,
                              right_orbital, right_orbital};
      if (test_case->density_pair_record[left_orbital][right_orbital] >= 0) {
        continue;
      }
      test_case->density_pair_record[left_orbital][right_orbital] =
          (int)test_case->record_count;
      CHECK(AddRecord(test_case,
                      MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX,
                      4, pair),
            "structure-factor density pair record");
    }
  }
  layout = (OracleLayout){3, 2, 1, 0};
  CHECK(oracle_build_basis(&layout, &basis), "electronic basis");
  for (index = 0; index < basis.sector_dimension; ++index) {
    const uint64_t word = basis.sector_words[index];
    test_case->psi[word] =
        (0.17 * (double)(index + 1) - 0.4) +
        I * (0.11 * (double)((index * 3 + 1) % 7) - 0.25);
  }
  test_case->psi[UINT64_C(13)] = 0.0;
  test_case->basis_log_scale[0] = 0.125;
  test_case->basis_log_scale[1] = -0.2;
  test_case->alpha[0] = 0.6 + 0.2 * I;
  test_case->alpha[1] = -0.3 + 0.7 * I;
}

static void BuildPureSpinCase(CorrectnessCase *test_case) {
  const size_t exchange[4] = {0, 1, 3, 2};
  const size_t exchange_adjoint[4] = {2, 3, 1, 0};
  const size_t density[4] = {0, 0, 0, 0};
  size_t left;
  size_t right;
  memset(test_case, 0, sizeof(*test_case));
  for (left = 0; left < ORACLE_MAX_ORBITALS; ++left) {
    for (right = 0; right < ORACLE_MAX_ORBITALS; ++right) {
      test_case->density_pair_record[left][right] = -1;
    }
  }
  test_case->name = "pure_spin_transverse";
  test_case->layout = (MVMCPowerLanczosObservableLayout){2, 1, 1, 1};
  CHECK(AddTerm(test_case, 0.8 + 0.3 * I, 4, exchange),
        "pure exchange term");
  CHECK(AddTerm(test_case, 0.8 - 0.3 * I, 4, exchange_adjoint),
        "pure exchange adjoint term");
  CHECK(AddTerm(test_case, 0.2, 2, density), "pure density term");
  CHECK(AddRecord(test_case,
                  MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALT, 4,
                  exchange),
        "pure exchange record");
  CHECK(AddRecord(test_case,
                  MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALT, 4,
                  exchange_adjoint),
        "pure exchange adjoint record");
  CHECK(AddRecord(test_case, MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, 2,
                  density),
        "pure density record");
  for (left = 0; left < 4; ++left) {
    for (right = 0; right < 4; ++right) {
      const size_t pair[4] = {left, left, right, right};
      test_case->density_pair_record[left][right] =
          (int)test_case->record_count;
      CHECK(AddRecord(test_case,
                      MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALT,
                      4, pair),
            "pure-spin density pair record");
    }
  }
  test_case->psi[6] = 1.0 + 0.2 * I;
  test_case->psi[9] = 0.0;
  test_case->basis_log_scale[0] = -0.1;
  test_case->basis_log_scale[1] = 0.15;
  test_case->alpha[0] = 0.4 - 0.3 * I;
  test_case->alpha[1] = 0.8 + 0.1 * I;
}

static void BuildParityPositiveCase(CorrectnessCase *test_case) {
  const size_t hop[4] = {0, 2, 0, 0};
  const size_t adjoint[4] = {2, 0, 0, 0};
  memset(test_case, 0, sizeof(*test_case));
  test_case->name = "one_body_parity_positive";
  test_case->layout = (MVMCPowerLanczosObservableLayout){3, 1, 0, 0};
  CHECK(AddTerm(test_case, 1.0, 2, hop), "positive parity hop");
  CHECK(AddTerm(test_case, 1.0, 2, adjoint),
        "positive parity hop adjoint");
  CHECK(AddRecord(test_case, MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS,
                  2, hop),
        "positive parity observable");
  test_case->psi[UINT64_C(4)] = 1.0;
  test_case->alpha[0] = 1.0 / sqrt(2.0);
  test_case->alpha[1] = 1.0 / sqrt(2.0);
}

static void BuildParityNegativeCase(CorrectnessCase *test_case) {
  const size_t hop[4] = {0, 2, 0, 0};
  const size_t adjoint[4] = {2, 0, 0, 0};
  memset(test_case, 0, sizeof(*test_case));
  test_case->name = "one_body_parity_negative";
  test_case->layout = (MVMCPowerLanczosObservableLayout){3, 2, 0, 0};
  CHECK(AddTerm(test_case, -1.0, 2, hop), "negative parity hop");
  CHECK(AddTerm(test_case, -1.0, 2, adjoint),
        "negative parity hop adjoint");
  CHECK(AddRecord(test_case, MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS,
                  2, hop),
        "negative parity observable");
  test_case->psi[UINT64_C(6)] = 1.0;
  test_case->alpha[0] = 1.0 / sqrt(2.0);
  test_case->alpha[1] = 1.0 / sqrt(2.0);
}

static void BuildOneBodyPauliZeroCase(CorrectnessCase *test_case) {
  const size_t hop[4] = {0, 1, 0, 0};
  const size_t density[4] = {0, 0, 0, 0};
  memset(test_case, 0, sizeof(*test_case));
  test_case->name = "one_body_pauli_zero";
  test_case->layout = (MVMCPowerLanczosObservableLayout){2, 2, 0, 0};
  CHECK(AddTerm(test_case, 1.0, 2, density), "Pauli-zero density term");
  CHECK(AddRecord(test_case, MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS,
                  2, hop),
        "Pauli-zero observable");
  test_case->psi[UINT64_C(3)] = 1.0;
  test_case->alpha[0] = 1.0 / sqrt(2.0);
  test_case->alpha[1] = 1.0 / sqrt(2.0);
}

static void BuildDistributedCase(CorrectnessCase *test_case) {
  const size_t density_up0[4] = {0, 0, 0, 0};
  const size_t hop_up[4] = {0, 1, 0, 0};
  const size_t hop_up_adjoint[4] = {1, 0, 0, 0};
  const size_t hop_down[4] = {2, 3, 0, 0};
  const size_t hop_down_adjoint[4] = {3, 2, 0, 0};
  OracleLayout layout;
  OracleBasis basis;
  size_t index;
  memset(test_case, 0, sizeof(*test_case));
  test_case->name = "distributed_electronic_block";
  test_case->layout = (MVMCPowerLanczosObservableLayout){2, 1, 1, 0};
  CHECK(AddTerm(test_case, 0.25, 2, density_up0),
        "distributed density term");
  CHECK(AddTerm(test_case, 0.6 + 0.2 * I, 2, hop_up),
        "distributed up hop");
  CHECK(AddTerm(test_case, 0.6 - 0.2 * I, 2, hop_up_adjoint),
        "distributed up hop adjoint");
  CHECK(AddTerm(test_case, -0.35 + 0.15 * I, 2, hop_down),
        "distributed down hop");
  CHECK(AddTerm(test_case, -0.35 - 0.15 * I, 2, hop_down_adjoint),
        "distributed down hop adjoint");
  CHECK(AddRecord(test_case, MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, 2,
                  density_up0),
        "distributed density record");
  CHECK(AddRecord(test_case, MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS, 2,
                  hop_up),
        "distributed hop record");
  layout = (OracleLayout){2, 1, 1, 0};
  CHECK(oracle_build_basis(&layout, &basis), "distributed basis");
  CHECK(basis.sector_dimension == 4, "distributed basis dimension");
  for (index = 0; index < basis.sector_dimension; ++index) {
    const uint64_t word = basis.sector_words[index];
    test_case->psi[word] =
        (0.09 * (double)(index + 1) - 0.2) +
        I * (0.07 * (double)((index * 5 + 2) % 11) - 0.3);
  }
  test_case->basis_log_scale[0] = 0.08;
  test_case->basis_log_scale[1] = -0.12;
}

static void CheckCombination(const char *label, const CaseResult *result,
                             const int *records,
                             const double complex *coefficients,
                             size_t count) {
  double complex production_final = 0.0;
  double complex oracle_final = 0.0;
  size_t entry;
  size_t term;
  for (term = 0; term < count; ++term) {
    production_final += coefficients[term] * result->final[records[term]];
    oracle_final +=
        coefficients[term] * result->oracle_final[records[term]];
  }
  CHECK(CloseComplex(production_final, oracle_final),
        "%s final representative", label);
  for (entry = 0; entry < 4; ++entry) {
    double complex production_matrix = 0.0;
    double complex oracle_matrix = 0.0;
    for (term = 0; term < count; ++term) {
      production_matrix += coefficients[term] *
                           result->matrix[4 * records[term] + entry];
      oracle_matrix += coefficients[term] *
                       result->oracle_matrix[4 * records[term] + entry];
    }
    CHECK(CloseComplex(production_matrix, oracle_matrix),
          "%s matrix entry %zu representative", label, entry);
  }
}

static void CheckRepresentativeRegistry(
    const CorrectnessCase *electronic_case,
    const CaseResult *electronic,
    const CorrectnessCase *pure_case, const CaseResult *pure) {
  const int density_record[1] = {0};
  const double complex unit[1] = {1.0};
  const int bond_records[2] = {1, 2};
  const double complex kinetic_coefficients[2] = {1.0, 1.0};
  const double complex current_coefficients[2] = {-I, I};
  const int density_density_record[1] = {5};
  const int transverse_records[2] = {0, 1};
  int longitudinal_records[4];
  const double complex longitudinal_coefficients[4] = {
      0.25, -0.25, -0.25, 0.25};
  int charge_structure_records[36];
  double complex charge_coefficients[36];
  int spin_structure_records[16];
  double complex spin_coefficients[16];
  const double charge_momentum = 2.0 * acos(-1.0) / 3.0;
  const double spin_momentum = acos(-1.0);
  size_t cursor = 0;
  size_t left;
  size_t right;
  CheckCombination("site/spin density", electronic, density_record, unit, 1);
  CheckCombination("bond kinetic energy", electronic, bond_records,
                   kinetic_coefficients, 2);
  CheckCombination("complex bond current", electronic, bond_records,
                   current_coefficients, 2);
  CheckCombination("density-density correlation", electronic,
                   density_density_record, unit, 1);
  CheckCombination("transverse spin exchange", pure, transverse_records,
                   kinetic_coefficients, 2);
  longitudinal_records[0] = pure_case->density_pair_record[0][1];
  longitudinal_records[1] = pure_case->density_pair_record[0][3];
  longitudinal_records[2] = pure_case->density_pair_record[2][1];
  longitudinal_records[3] = pure_case->density_pair_record[2][3];
  CheckCombination("longitudinal spin correlation", pure,
                   longitudinal_records, longitudinal_coefficients, 4);
  for (left = 0; left < 6; ++left) {
    for (right = 0; right < 6; ++right) {
      const int left_site = (int)(left % 3);
      const int right_site = (int)(right % 3);
      const double complex phase =
          cexp(I * charge_momentum *
               (double)(left_site - right_site)) / 3.0;
      charge_structure_records[cursor] =
          electronic_case->density_pair_record[left][right];
      charge_coefficients[cursor] = phase;
      ++cursor;
    }
  }
  CheckCombination("fixed-momentum charge structure factor", electronic,
                   charge_structure_records, charge_coefficients, cursor);
  cursor = 0;
  for (left = 0; left < 4; ++left) {
    for (right = 0; right < 4; ++right) {
      const int left_site = (int)(left % 2);
      const int right_site = (int)(right % 2);
      const int left_spin_sign = left < 2 ? 1 : -1;
      const int right_spin_sign = right < 2 ? 1 : -1;
      const double complex phase =
          cexp(I * spin_momentum *
               (double)(left_site - right_site)) / 2.0;
      spin_structure_records[cursor] =
          pure_case->density_pair_record[left][right];
      spin_coefficients[cursor] =
          0.25 * (double)(left_spin_sign * right_spin_sign) * phase;
      ++cursor;
    }
  }
  CheckCombination("fixed-momentum spin structure factor", pure,
                   spin_structure_records, spin_coefficients, cursor);
}

static void CheckControls(const CorrectnessCase *electronic_case,
                          const CorrectnessCase *pure_case,
                          const CaseResult *lower,
                          const CaseResult *electronic,
                          const CaseResult *pure,
                          const CaseResult *parity_positive,
                          const CaseResult *parity_negative,
                          const CaseResult *pauli_zero) {
  size_t entry;
  CHECK(CloseComplex(lower->matrix[0], 0.0) &&
            CloseComplex(lower->matrix[1], 0.0) &&
            CloseComplex(lower->matrix[2], 1.0) &&
            CloseComplex(lower->matrix[3], 0.0),
        "mandatory lower-triangle counterexample");
  CHECK(CloseComplex(lower->final[0], 0.5),
        "mandatory lower-triangle direct final");
  for (entry = 0; entry < 4; ++entry) {
    CHECK(CloseComplex(electronic->matrix[4 * 8 + entry],
                       electronic->matrix[4 * 6 + entry]),
          "quartic cross-family action entry %zu", entry);
    CHECK(CloseComplex(electronic->matrix[4 * 9 + 2 * (entry % 2) +
                                             entry / 2],
                       conj(electronic->matrix[4 * 6 + entry])),
          "quartic adjoint relation entry %zu", entry);
  }
  for (entry = 0; entry < 4; ++entry) {
    CHECK(CloseComplex(electronic->matrix[4 * 6 + entry] -
                           electronic->matrix[4 * 10 + entry],
                       0.0),
          "same-family quartic exact cancellation entry %zu", entry);
  }
  CHECK(CloseComplex(electronic->matrix[0],
                     conj(electronic->matrix[0])) &&
            CloseComplex(electronic->matrix[1],
                         conj(electronic->matrix[2])) &&
            CloseComplex(electronic->matrix[3],
                         conj(electronic->matrix[3])),
        "self-adjoint density Hermiticity");
  CHECK(fabs(cimag(electronic->matrix[4])) > 1.0e-12,
        "non-Hermitian hopping lost complex diagonal");
  CHECK(CloseComplex(electronic->matrix[4 * 7 + 0], 0.0) &&
            CloseComplex(electronic->matrix[4 * 7 + 1], 0.0) &&
            CloseComplex(electronic->matrix[4 * 7 + 2], 0.0) &&
            CloseComplex(electronic->matrix[4 * 7 + 3], 0.0),
        "repeated-orbital Pauli zero");
  CHECK(CloseComplex(pure->matrix[4], conj(pure->matrix[0])) &&
            CloseComplex(pure->matrix[5], conj(pure->matrix[2])) &&
            CloseComplex(pure->matrix[6], conj(pure->matrix[1])) &&
            CloseComplex(pure->matrix[7], conj(pure->matrix[3])),
        "pure-spin exchange adjoint relation");
  CHECK(CloseComplex(parity_positive->matrix[0], 0.0) &&
            CloseComplex(parity_positive->matrix[1], 0.0) &&
            CloseComplex(parity_positive->matrix[2], 1.0) &&
            CloseComplex(parity_positive->matrix[3], 0.0),
        "one-body positive fermion parity");
  CHECK(CloseComplex(parity_negative->matrix[0], 0.0) &&
            CloseComplex(parity_negative->matrix[1], 0.0) &&
            CloseComplex(parity_negative->matrix[2], -1.0) &&
            CloseComplex(parity_negative->matrix[3], 0.0),
        "one-body negative fermion parity");
  CHECK(CloseComplex(pauli_zero->matrix[0], 0.0) &&
            CloseComplex(pauli_zero->matrix[1], 0.0) &&
            CloseComplex(pauli_zero->matrix[2], 0.0) &&
            CloseComplex(pauli_zero->matrix[3], 0.0) &&
            CloseComplex(pauli_zero->final[0], 0.0),
        "one-body occupied-target Pauli zero");
  CheckRepresentativeRegistry(electronic_case, electronic,
                              pure_case, pure);
}

static const char *FamilyName(MVMCPowerLanczosObservableFamily family) {
  switch (family) {
    case MVMC_POWER_LANCZOS_OBSERVABLE_CISAJS:
      return "cisajs";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALTEX:
      return "cisajscktaltex";
    case MVMC_POWER_LANCZOS_OBSERVABLE_CISAJSCKTALT:
      return "cisajscktalt";
  }
  return "invalid";
}

static void PrintComplexJson(double complex value) {
  printf("[%.17g,%.17g]", creal(value), cimag(value));
}

static void EmitGroupJson(const CorrectnessCase *test_case,
                          const CaseResult *result, const char *model,
                          const char *arithmetic, const char *support) {
  size_t record;
  printf("{\"name\":\"%s\",\"model\":\"%s\","
         "\"arithmetic\":\"%s\",\"support\":\"%s\","
         "\"layout\":{\"nsite\":%d,\"up_electron_count\":%d,"
         "\"down_electron_count\":%d,\"pure_spin\":%s},"
         "\"alpha\":[",
         test_case->name, model, arithmetic, support,
         test_case->layout.nsite, test_case->layout.up_electron_count,
         test_case->layout.down_electron_count,
         test_case->layout.pure_spin ? "true" : "false");
  PrintComplexJson(test_case->alpha[0]);
  printf(",");
  PrintComplexJson(test_case->alpha[1]);
  printf("],\"record_count\":%zu,\"records\":[",
         test_case->record_count);
  for (record = 0; record < test_case->record_count; ++record) {
    const MVMCPowerLanczosObservableRecord *observable =
        &test_case->records[record];
    size_t entry;
    int raw_index;
    if (record != 0) printf(",");
    printf("{\"ordinal\":%zu,\"family\":\"%s\","
           "\"row_width\":%d,\"raw_indices\":[",
           record, FamilyName(observable->family), observable->row_width);
    for (raw_index = 0; raw_index < observable->row_width; ++raw_index) {
      if (raw_index != 0) printf(",");
      printf("%d", observable->raw_indices[raw_index]);
    }
    printf("],\"production_matrix\":[");
    for (entry = 0; entry < 4; ++entry) {
      if (entry != 0) printf(",");
      PrintComplexJson(result->matrix[4 * record + entry]);
    }
    printf("],\"oracle_matrix\":[");
    for (entry = 0; entry < 4; ++entry) {
      if (entry != 0) printf(",");
      PrintComplexJson(result->oracle_matrix[4 * record + entry]);
    }
    printf("],\"production_direct_final\":");
    PrintComplexJson(result->final[record]);
    printf(",\"oracle_direct_final\":");
    PrintComplexJson(result->oracle_final[record]);
    printf("}");
  }
  printf("]}");
}

static void EmitSnapshotJson(
    int world_size, uint64_t semantic_hash,
    const CorrectnessCase *lower_case, const CaseResult *lower_result,
    const CorrectnessCase *electronic_case,
    const CaseResult *electronic_result,
    const CorrectnessCase *pure_case, const CaseResult *pure_result,
    const CorrectnessCase *parity_positive_case,
    const CaseResult *parity_positive_result,
    const CorrectnessCase *parity_negative_case,
    const CaseResult *parity_negative_result,
    const CorrectnessCase *pauli_zero_case,
    const CaseResult *pauli_zero_result) {
  printf("{\"schema_version\":2,"
         "\"protocol_id\":\"p6c2-representative-correctness-v1\","
         "\"world_size\":%d,\"semantic_fnv64\":\"%016llx\","
         "\"passed\":true,\"groups\":[",
         world_size, (unsigned long long)semantic_hash);
  EmitGroupJson(lower_case, lower_result, "electronic", "real", "regular");
  printf(",");
  EmitGroupJson(electronic_case, electronic_result, "electronic",
                "complex", "zero_support");
  printf(",");
  EmitGroupJson(pure_case, pure_result, "pure_spin", "complex",
                "zero_support");
  printf(",");
  EmitGroupJson(parity_positive_case, parity_positive_result, "electronic",
                "real", "regular");
  printf(",");
  EmitGroupJson(parity_negative_case, parity_negative_result, "electronic",
                "real", "regular");
  printf(",");
  EmitGroupJson(pauli_zero_case, pauli_zero_result, "electronic", "real",
                "regular");
  printf("]}\n");
}

int main(int argc, char **argv) {
  CorrectnessCase lower_case;
  CorrectnessCase electronic_case;
  CorrectnessCase pure_case;
  CorrectnessCase parity_positive_case;
  CorrectnessCase parity_negative_case;
  CorrectnessCase pauli_zero_case;
  CorrectnessCase distributed_case;
  CaseResult lower_result;
  CaseResult electronic_result;
  CaseResult pure_result;
  CaseResult parity_positive_result;
  CaseResult parity_negative_result;
  CaseResult pauli_zero_result;
  uint64_t semantic_hash;
  int emit_json = 0;
  int world_size = 1;
  int argument;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif
  for (argument = 1; argument < argc; ++argument) {
    if (strcmp(argv[argument], "--emit-json") == 0) {
      emit_json = 1;
    } else {
      CHECK(0, "unknown argument: %s", argv[argument]);
    }
  }
  BuildLowerCase(&lower_case);
  BuildElectronicCase(&electronic_case);
  BuildPureSpinCase(&pure_case);
  BuildParityPositiveCase(&parity_positive_case);
  BuildParityNegativeCase(&parity_negative_case);
  BuildOneBodyPauliZeroCase(&pauli_zero_case);
  BuildDistributedCase(&distributed_case);
  CHECK(RunCase(&lower_case, UINT64_C(0x2001), &lower_result),
        "lower case execution");
  CHECK(RunCase(&electronic_case, UINT64_C(0x2002), &electronic_result),
        "electronic case execution");
  CHECK(RunCase(&pure_case, UINT64_C(0x2003), &pure_result),
        "pure case execution");
  CHECK(RunCase(&parity_positive_case, UINT64_C(0x2005),
                &parity_positive_result),
        "positive parity case execution");
  CHECK(RunCase(&parity_negative_case, UINT64_C(0x2006),
                &parity_negative_result),
        "negative parity case execution");
  CHECK(RunCase(&pauli_zero_case, UINT64_C(0x2007),
                &pauli_zero_result),
        "one-body Pauli-zero case execution");
  CHECK(RunDistributedBlockCase(&distributed_case, UINT64_C(0x2004),
                                world_size),
        "distributed block case execution");
  CheckControls(&electronic_case, &pure_case, &lower_result,
                &electronic_result, &pure_result,
                &parity_positive_result, &parity_negative_result,
                &pauli_zero_result);
  semantic_hash = lower_result.semantic_hash;
  semantic_hash = HashBytes(semantic_hash, &electronic_result.semantic_hash,
                            sizeof(electronic_result.semantic_hash));
  semantic_hash = HashBytes(semantic_hash, &pure_result.semantic_hash,
                            sizeof(pure_result.semantic_hash));
  semantic_hash = HashBytes(semantic_hash,
                            &parity_positive_result.semantic_hash,
                            sizeof(parity_positive_result.semantic_hash));
  semantic_hash = HashBytes(semantic_hash,
                            &parity_negative_result.semantic_hash,
                            sizeof(parity_negative_result.semantic_hash));
  semantic_hash = HashBytes(semantic_hash,
                            &pauli_zero_result.semantic_hash,
                            sizeof(pauli_zero_result.semantic_hash));
#ifdef _mpi_use
  {
    uint64_t minimum_hash = 0;
    uint64_t maximum_hash = 0;
    int global_failures = 0;
    MPI_Allreduce(&semantic_hash, &minimum_hash, 1, MPI_UINT64_T, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&semantic_hash, &maximum_hash, 1, MPI_UINT64_T, MPI_MAX,
                  MPI_COMM_WORLD);
    CHECK(minimum_hash == maximum_hash,
          "MPI semantic payload differs across ranks");
    MPI_Allreduce(&failures, &global_failures, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    failures = global_failures;
  }
#endif
  if (world_rank == 0 && failures == 0) {
    if (emit_json) {
      EmitSnapshotJson(
          world_size, semantic_hash, &lower_case, &lower_result,
          &electronic_case, &electronic_result, &pure_case, &pure_result,
          &parity_positive_case, &parity_positive_result,
          &parity_negative_case, &parity_negative_result,
          &pauli_zero_case, &pauli_zero_result);
    } else {
      printf("power_lanczos_observable_correctness: PASS cases=7 "
             "records=71 representatives=8 distributed_sector=4 "
             "mpi=%d hash=%016llx\n",
             world_size, (unsigned long long)semantic_hash);
    }
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
