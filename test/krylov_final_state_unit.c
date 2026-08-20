#include "krylov_final_state.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

static int failures = 0;
static int world_rank = 0;

#define CHECK(condition, ...)                                                 \
  do {                                                                        \
    if (!(condition)) {                                                       \
      fprintf(stderr, "KrylovFinalState_Unit FAIL rank %d: ", world_rank);   \
      fprintf(stderr, __VA_ARGS__);                                           \
      fprintf(stderr, " (line %d)\n", __LINE__);                            \
      ++failures;                                                             \
    }                                                                         \
  } while (0)

typedef struct {
  double complex values[64];
  uint64_t numeric_zero_mask;
  MVMCKrylovStatus forced_status;
} TableAmplitude;

static int close_double(double actual, double expected, double tolerance) {
  return fabs(actual - expected) <= tolerance * fmax(1.0, fabs(expected));
}

static int close_complex(double complex actual, double complex expected,
                         double tolerance) {
  return cabs(actual - expected) <= tolerance * fmax(1.0, cabs(expected));
}

static MVMCKrylovStatus table_callback(
    const uint64_t *words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
  TableAmplitude *table = (TableAmplitude *)context;
  const size_t index = word_count == 1
                           ? (size_t)(words[0] & UINT64_C(63))
                           : 63;
  if (table == NULL || result == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  if (table->forced_status != MVMC_KRYLOV_STATUS_OK) {
    return table->forced_status;
  }
  memset(result, 0, sizeof(*result));
  if ((table->numeric_zero_mask & (UINT64_C(1) << index)) != 0) {
    if (mvmc_scaled_complex_make_numeric_zero(
            log(1.0e-12), 0.0, -INFINITY, 0.0, &result->value) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    result->numeric_zero_component_count = 1;
  } else if (creal(table->values[index]) == 0.0 &&
             cimag(table->values[index]) == 0.0) {
    if (mvmc_scaled_complex_make_exact_zero(&result->value) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    result->exact_zero_component_count = 1;
  } else {
    if (mvmc_scaled_complex_from_raw_testing(
            table->values[index], &result->value) !=
        MVMC_PFAFFIAN_STATUS_OK) {
      return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    }
    result->well_pivoted_component_count = 1;
  }
  result->local_factorization_count = 1;
  result->global_factorization_count = 1;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovFockModel two_state_model(
    const MVMCKrylovHamiltonianTerm *terms,
    const MVMCKrylovFermionOperator *operators) {
  MVMCKrylovFockModel model;
  memset(&model, 0, sizeof(model));
  model.site_count = 2;
  model.up_electron_count = 1;
  model.hermitian = 1;
  model.terms = terms;
  model.term_count = 2;
  model.operators = operators;
  model.operator_count = 4;
  return model;
}

static MVMCKrylovBoundedLimits bounded_limits(int max_order) {
  MVMCKrylovBoundedLimits limits;
  memset(&limits, 0, sizeof(limits));
  limits.amplitude_policy_hash = UINT64_C(0x5035415441424c45);
  limits.cache_bytes = 4096;
  limits.max_row_transitions = 16;
  limits.max_workspace_bytes = (size_t)1024 * 1024;
  limits.max_node_expansions = 1000;
  limits.max_terminal_amplitude_calls = 1000;
  limits.max_total_row_transitions = 10000;
  limits.max_order = max_order;
  return limits;
}

static MVMCKrylovBoundedWorkspace *create_workspace(
    int max_order, MVMCKrylovBoundedPlan **plan_out) {
  static const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
  };
  static const MVMCKrylovHamiltonianTerm terms[] = {
      {1.0, 0, 2, 1, 0},
      {1.0, 2, 2, 1, 1},
  };
  const MVMCKrylovFockModel model = two_state_model(terms, operators);
  const MVMCKrylovBoundedLimits limits = bounded_limits(max_order);
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = NULL;
  CHECK(mvmc_bounded_krylov_plan_create(&model, &limits, &plan) ==
            MVMC_KRYLOV_STATUS_OK &&
            plan != NULL,
        "two-state plan creation");
  if (plan != NULL) {
    CHECK(mvmc_bounded_krylov_workspace_create(plan, &workspace) ==
              MVMC_KRYLOV_STATUS_OK &&
              workspace != NULL,
          "two-state workspace creation");
  }
  *plan_out = plan;
  return workspace;
}

static void destroy_workspace(MVMCKrylovBoundedWorkspace *workspace,
                              MVMCKrylovBoundedPlan *plan) {
  mvmc_bounded_krylov_workspace_destroy(workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
}

static int export_scaled(const MVMCScaledComplex *value,
                         double complex *raw) {
  const MVMCScaledComplexExportStatus status =
      mvmc_scaled_complex_export_common_scale(value, 0.0, raw);
  return status == MVMC_SCALED_EXPORT_OK ||
         status == MVMC_SCALED_EXPORT_EXACT_ZERO;
}

static MVMCKrylovBoundedResult raw_krylov_result(
    const double complex *values, int evaluated_order) {
  MVMCKrylovBoundedResult result;
  int order;
  memset(&result, 0, sizeof(result));
  result.valid = 1;
  result.status = MVMC_KRYLOV_STATUS_OK;
  result.evaluated_order = evaluated_order;
  for (order = 0; order <= evaluated_order; ++order) {
    MVMCPfaffianStatus status;
    if (creal(values[order]) == 0.0 && cimag(values[order]) == 0.0) {
      status = mvmc_scaled_complex_make_exact_zero(&result.value[order]);
    } else {
      status = mvmc_scaled_complex_from_raw_testing(
          values[order], &result.value[order]);
    }
    CHECK(status == MVMC_PFAFFIAN_STATUS_OK,
          "raw Krylov conversion order %d", order);
  }
  return result;
}

static void test_policy_contract(void) {
  const double complex p1[] = {1.0, -1.0};
  const double complex p2[] = {1.0, 0.25 * I, -0.5};
  const double complex zero[] = {0.0, 0.0};
  const double complex nonfinite[] = {1.0, NAN + 0.0 * I};
  const double complex infinite[] = {1.0, INFINITY + 0.0 * I};
  MVMCKrylovFinalStatePolicy exact;
  MVMCKrylovFinalStatePolicy synthetic;
  MVMCKrylovFinalStatePolicy perturbed;
  MVMCKrylovFinalStatePolicy invalid;

  CHECK(mvmc_krylov_final_state_policy_create(
            1, MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT,
            UINT64_C(0x4558414354503101), p1, 2, &exact) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_final_state_policy_hash(&exact) != 0,
        "exact policy creation");
  CHECK(mvmc_krylov_final_state_policy_create(
            2, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC,
            UINT64_C(0x53594e5448503201), p2, 3, &synthetic) ==
            MVMC_KRYLOV_STATUS_OK &&
            exact.policy_hash != synthetic.policy_hash,
        "synthetic p2 policy creation");
  CHECK(mvmc_krylov_final_state_policy_create(
            1, MVMC_KRYLOV_FINAL_COEFFICIENT_PERTURBED_TESTING,
            UINT64_C(0x5045525455524201), p1, 2, &perturbed) ==
            MVMC_KRYLOV_STATUS_OK &&
            perturbed.policy_hash != exact.policy_hash &&
            perturbed.policy_hash != synthetic.policy_hash,
        "perturbed Testing-only policy creation");
  CHECK(mvmc_krylov_final_state_policy_create(
            -1, MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT, 1, p1, 2,
            &invalid) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "negative order rejected");
  CHECK(mvmc_krylov_final_state_policy_create(
            0, MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT, 1, p1, 1,
            &invalid) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "order zero rejected");
  CHECK(mvmc_krylov_final_state_policy_create(
            3, MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT, 1, p2, 3,
            &invalid) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "order three rejected because H Psi needs v4");
  CHECK(mvmc_krylov_final_state_policy_create(
            1, (MVMCKrylovFinalCoefficientSource)99, 1, p1, 2,
            &invalid) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "unknown coefficient source rejected");
  CHECK(mvmc_krylov_final_state_policy_create(
            1, MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT, 0, p1, 2,
            &invalid) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "zero provenance rejected");
  CHECK(mvmc_krylov_final_state_policy_create(
            1, MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT, 1, p1, 1,
            &invalid) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "coefficient count mismatch rejected");
  CHECK(mvmc_krylov_final_state_policy_create(
            1, MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT, 1, zero, 2,
            &invalid) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "all-zero coefficients rejected");
  CHECK(mvmc_krylov_final_state_policy_create(
            1, MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT, 1, nonfinite, 2,
            &invalid) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "nonfinite coefficients rejected");
  CHECK(mvmc_krylov_final_state_policy_create(
            1, MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT, 1, infinite, 2,
            &invalid) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "infinite coefficients rejected");
  exact.policy_hash ^= UINT64_C(1);
  CHECK(mvmc_krylov_final_state_policy_hash(&exact) == 0,
        "corrupt policy hash rejected");
}

static void test_exact_support_bridge(void) {
  const double scale = 1.0 / sqrt(2.0);
  const double complex coefficient[] = {scale, -scale};
  const double complex restricted[] = {1.0, 0.0};
  const uint64_t states[] = {UINT64_C(1), UINT64_C(2)};
  MVMCKrylovFinalStatePolicy policy;
  MVMCKrylovFinalStatePolicy restricted_policy;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = create_workspace(2, &plan);
  TableAmplitude amplitude;
  double norm = 0.0;
  double energy = 0.0;
  double energy2 = 0.0;
  double restricted_norm = 0.0;
  double complex restricted_energy_numerator = 0.0;
  double restricted_hamiltonian_norm = 0.0;
  double restricted_direct_energy = 0.0;
  double restricted_direct_energy2 = 0.0;
  size_t state;

  memset(&amplitude, 0, sizeof(amplitude));
  amplitude.values[1] = 1.0;
  CHECK(mvmc_krylov_final_state_policy_create(
            1, MVMC_KRYLOV_FINAL_COEFFICIENT_EXACT,
            UINT64_C(0x3253544154454745), coefficient, 2, &policy) ==
            MVMC_KRYLOV_STATUS_OK,
        "support bridge exact policy");
  CHECK(mvmc_krylov_final_state_policy_create(
            1, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC,
            UINT64_C(0x5245535452494354), restricted, 2,
            &restricted_policy) == MVMC_KRYLOV_STATUS_OK,
        "restricted synthetic policy");
  for (state = 0; state < 2 && workspace != NULL; ++state) {
    MVMCKrylovBoundedResult krylov;
    MVMCKrylovFinalStateEvaluation final_state;
    double complex local_energy = NAN + I * NAN;
    const uint64_t words = states[state];
    double weight;
    CHECK(mvmc_bounded_krylov_evaluate(
              workspace, &words, 1, table_callback, &amplitude, &krylov) ==
              MVMC_KRYLOV_STATUS_OK,
          "support bridge bounded evaluation state %zu", state);
    CHECK(mvmc_krylov_final_state_evaluate(
              &policy, &krylov, &final_state) == MVMC_KRYLOV_STATUS_OK &&
              final_state.valid && final_state.sampleable &&
              final_state.local_energy_available,
          "support bridge final evaluation state %zu", state);
    CHECK(mvmc_krylov_final_state_local_energy_export(
              &final_state, &local_energy) == MVMC_SCALED_EXPORT_OK &&
              close_complex(local_energy, -1.0, 2.0e-13),
          "support bridge local energy state %zu=(%.17g,%.17g)", state,
          creal(local_energy), cimag(local_energy));
    if (state == 1) {
      CHECK(krylov.value[0].state == MVMC_SCALED_COMPLEX_EXACT_ZERO &&
                final_state.amplitude.state ==
                    MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
                final_state.basis_exact_zero_count == 2,
            "Psi0 zero/Psi1 nonzero support bridge");
    }
    weight = exp(final_state.log_weight);
    norm += weight;
    energy += weight * creal(local_energy);
    energy2 += weight * creal(local_energy * conj(local_energy));

    CHECK(mvmc_krylov_final_state_evaluate(
              &restricted_policy, &krylov, &final_state) ==
              MVMC_KRYLOV_STATUS_OK,
          "restricted final evaluation state %zu", state);
    CHECK(final_state.coefficient_zero_count == 1 &&
              final_state.basis_numeric_zero_count == 0,
          "restricted coefficient/basis zero counts state %zu", state);
    {
      double complex psi_value = NAN + I * NAN;
      double complex hpsi_value = NAN + I * NAN;
      CHECK(export_scaled(&final_state.amplitude, &psi_value) &&
                export_scaled(&final_state.hamiltonian_amplitude,
                              &hpsi_value),
            "restricted amplitude export state %zu", state);
      restricted_norm += creal(conj(psi_value) * psi_value);
      restricted_energy_numerator += conj(psi_value) * hpsi_value;
      restricted_hamiltonian_norm +=
          creal(conj(hpsi_value) * hpsi_value);
    }
    if (state == 0) {
      CHECK(final_state.sampleable &&
                mvmc_krylov_final_state_local_energy_export(
                    &final_state, &local_energy) ==
                    MVMC_SCALED_EXPORT_EXACT_ZERO &&
                creal(local_energy) == 0.0 && cimag(local_energy) == 0.0,
            "restricted support local energy zero");
      restricted_direct_energy += creal(local_energy);
      restricted_direct_energy2 +=
          creal(local_energy * conj(local_energy));
    } else {
      CHECK(final_state.amplitude.state == MVMC_SCALED_COMPLEX_EXACT_ZERO &&
                final_state.hamiltonian_amplitude.state ==
                    MVMC_SCALED_COMPLEX_FINITE_NONZERO &&
                !final_state.sampleable &&
                !final_state.local_energy_available,
            "Psi_p zero/H Psi_p nonzero classified separately");
    }
  }
  CHECK(close_double(norm, 1.0, 2.0e-13) &&
            close_double(energy / norm, -1.0, 2.0e-13) &&
            close_double(energy2 / norm - (energy / norm) *
                                           (energy / norm),
                         0.0, 2.0e-13),
        "exact direct enumeration energy/variance");
  CHECK(close_double(restricted_norm, 1.0, 2.0e-13) &&
            close_complex(restricted_energy_numerator, 0.0, 2.0e-13) &&
            close_double(restricted_hamiltonian_norm / restricted_norm -
                             pow(cabs(restricted_energy_numerator /
                                      restricted_norm),
                                 2.0),
                         1.0, 2.0e-13) &&
            close_double(restricted_direct_energy2 -
                             restricted_direct_energy *
                                 restricted_direct_energy,
                         0.0, 2.0e-13),
        "full-support variance one/direct restricted variance zero");
  destroy_workspace(workspace, plan);
}

static void test_synthetic_proposal_transaction(void) {
  static const MVMCKrylovFermionOperator operators[] = {
      {MVMC_KRYLOV_FERMION_CREATE, 0},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 1},
      {MVMC_KRYLOV_FERMION_CREATE, 1},
      {MVMC_KRYLOV_FERMION_ANNIHILATE, 0},
  };
  static const MVMCKrylovHamiltonianTerm terms[] = {
      {1.0, 0, 2, 1, 0},
      {1.0, 2, 2, 1, 1},
  };
  const MVMCKrylovFockModel model = two_state_model(terms, operators);
  const double complex coefficient[] = {1.0, 2.0};
  const double complex zero_target_coefficient[] = {1.0, 0.0};
  MVMCKrylovFinalStatePolicy policy;
  MVMCKrylovFinalStatePolicy zero_target_policy;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *workspace = create_workspace(2, &plan);
  TableAmplitude amplitude;
  MVMCKrylovFinalStateSnapshot current;
  MVMCKrylovFinalStateProposalStepResult step;
  uint64_t current_words[1] = {UINT64_C(1)};
  uint64_t proposal_words[1] = {0};

  memset(&amplitude, 0, sizeof(amplitude));
  amplitude.values[1] = 1.0;
  CHECK(mvmc_krylov_final_state_policy_create(
            1, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC,
            UINT64_C(0x5745494748543134), coefficient, 2, &policy) ==
            MVMC_KRYLOV_STATUS_OK,
        "proposal policy");
  CHECK(mvmc_krylov_final_state_sampler_initialize(
            workspace, &policy, current_words, 1, table_callback,
            &amplitude, &current) == MVMC_KRYLOV_STATUS_OK,
        "proposal sampler initialize");
  CHECK(mvmc_krylov_final_state_sampler_step_selected_neighbor(
            workspace, &policy, &model, current_words, 1, &current, 0,
            0.999, table_callback, &amplitude, proposal_words, 1,
            &step) == MVMC_KRYLOV_STATUS_OK &&
            step.valid && step.configuration_changed &&
            step.step.accepted && current_words[0] == UINT64_C(2) &&
            current.accepted_generation == 1 &&
            current.integrity_hash != 0 &&
            memcmp(&current.krylov, &step.step.proposal_krylov,
                   sizeof(current.krylov)) == 0 &&
            memcmp(&current.final_state, &step.step.proposal_final_state,
                   sizeof(current.final_state)) == 0 &&
            close_double(step.step.acceptance.log_target_ratio, log(4.0),
                         2.0e-13),
        "weight 1:4 forward proposal accepted");
  CHECK(mvmc_krylov_final_state_sampler_step_selected_neighbor(
            workspace, &policy, &model, current_words, 1, &current, 0,
            0.30, table_callback, &amplitude, proposal_words, 1,
            &step) == MVMC_KRYLOV_STATUS_OK &&
            step.valid && !step.configuration_changed &&
            !step.step.accepted && current_words[0] == UINT64_C(2) &&
            current.accepted_generation == 1 &&
            close_double(step.step.acceptance.log_target_ratio, -log(4.0),
                         2.0e-13),
        "reverse probability 1/4 rejected above threshold");

  {
    const uint64_t saved_generation = current.accepted_generation;
    MVMCKrylovFinalStateSnapshot saved;
    current.accepted_generation = UINT64_MAX;
    saved = current;
    CHECK(mvmc_krylov_final_state_sampler_step_selected_neighbor(
              workspace, &policy, &model, current_words, 1, &current, 0,
              0.0, table_callback, &amplitude, proposal_words, 1,
              &step) == MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
              !step.valid && current_words[0] == UINT64_C(2) &&
              memcmp(&current, &saved, sizeof(current)) == 0,
          "accepted-generation overflow is transactional");
    current.accepted_generation = saved_generation;
  }

  {
    const uint64_t saved_hash = current.configuration_hash;
    const uint64_t saved_generation = current.accepted_generation;
    const MVMCKrylovFinalStateSnapshot saved = current;
    amplitude.forced_status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    CHECK(mvmc_krylov_final_state_sampler_step_selected_neighbor(
              workspace, &policy, &model, current_words, 1, &current, 0,
              0.0, table_callback, &amplitude, proposal_words, 1,
              &step) == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
              !step.valid && current_words[0] == UINT64_C(2) &&
              current.configuration_hash == saved_hash &&
              current.accepted_generation == saved_generation &&
              memcmp(&current, &saved, sizeof(current)) == 0,
          "callback failure leaves accepted state unchanged");
    amplitude.forced_status = MVMC_KRYLOV_STATUS_OK;
  }
  destroy_workspace(workspace, plan);

  workspace = create_workspace(2, &plan);
  current_words[0] = UINT64_C(1);
  CHECK(mvmc_krylov_final_state_policy_create(
            1, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC,
            UINT64_C(0x5a45524f54415247), zero_target_coefficient, 2,
            &zero_target_policy) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_final_state_sampler_initialize(
                workspace, &zero_target_policy, current_words, 1,
                table_callback, &amplitude, &current) ==
                MVMC_KRYLOV_STATUS_OK,
        "zero-target sampler initialize");
  CHECK(mvmc_krylov_final_state_sampler_step_selected_neighbor(
            workspace, &zero_target_policy, &model, current_words, 1,
            &current, 0, 0.0, table_callback, &amplitude, proposal_words, 1,
            &step) == MVMC_KRYLOV_STATUS_OK &&
            step.valid && !step.configuration_changed &&
            step.step.acceptance.exact_zero_proposal &&
            current_words[0] == UINT64_C(1),
        "exact-zero target proposal is a deterministic reject");

  amplitude.numeric_zero_mask = UINT64_C(1) << 2;
  CHECK(mvmc_krylov_final_state_sampler_step_selected_neighbor(
            workspace, &zero_target_policy, &model, current_words, 1,
            &current, 0, 0.0, table_callback, &amplitude, proposal_words, 1,
            &step) == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE &&
            !step.valid && current_words[0] == UINT64_C(1),
        "numeric-zero target proposal fails loud");
  amplitude.numeric_zero_mask = 0;
  destroy_workspace(workspace, plan);
}

static void test_complex_p2_and_failures(void) {
  const double complex values[] = {
      1.0 + 0.2 * I, -0.5 + 0.3 * I, 0.7 - 0.1 * I, 0.2 + 0.9 * I};
  const double complex coefficient[] = {
      0.8 + 0.1 * I, -0.2 + 0.4 * I, 0.3 - 0.5 * I};
  double complex rotated[3];
  const double complex phase = cos(0.7) + I * sin(0.7);
  MVMCKrylovBoundedResult krylov = raw_krylov_result(values, 3);
  MVMCKrylovFinalStatePolicy policy;
  MVMCKrylovFinalStatePolicy rotated_policy;
  MVMCKrylovFinalStateEvaluation evaluation;
  MVMCKrylovFinalStateEvaluation rotated_evaluation;
  double complex expected_amplitude = 0.0;
  double complex expected_hamiltonian = 0.0;
  double complex actual_amplitude = NAN + I * NAN;
  double complex local_energy = NAN + I * NAN;
  int index;

  for (index = 0; index < 3; ++index) {
    expected_amplitude += coefficient[index] * values[index];
    expected_hamiltonian += coefficient[index] * values[index + 1];
    rotated[index] = phase * coefficient[index];
  }
  CHECK(mvmc_krylov_final_state_policy_create(
            2, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC,
            UINT64_C(0x434f4d504c455832), coefficient, 3, &policy) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_final_state_evaluate(
                &policy, &krylov, &evaluation) == MVMC_KRYLOV_STATUS_OK &&
            evaluation.sampleable &&
            export_scaled(&evaluation.amplitude, &actual_amplitude) &&
            mvmc_krylov_final_state_local_energy_export(
                &evaluation, &local_energy) == MVMC_SCALED_EXPORT_OK,
        "complex p2 evaluation");
  CHECK(close_complex(actual_amplitude, expected_amplitude, 5.0e-13) &&
            close_complex(local_energy,
                          expected_hamiltonian / expected_amplitude,
                          5.0e-13),
        "complex p2 raw oracle");
  CHECK(mvmc_krylov_final_state_policy_create(
            2, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC,
            UINT64_C(0x434f4d504c455833), rotated, 3,
            &rotated_policy) == MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_final_state_evaluate(
                &rotated_policy, &krylov, &rotated_evaluation) ==
                MVMC_KRYLOV_STATUS_OK,
        "rotated complex policy evaluation");
  {
    double complex rotated_local = NAN + I * NAN;
    CHECK(close_double(rotated_evaluation.log_weight,
                       evaluation.log_weight, 5.0e-13) &&
              mvmc_krylov_final_state_local_energy_export(
                  &rotated_evaluation, &rotated_local) ==
                  MVMC_SCALED_EXPORT_OK &&
              close_complex(rotated_local, local_energy, 5.0e-13),
          "global coefficient phase invariance");
  }
  {
    MVMCKrylovFinalStateAcceptance acceptance;
    MVMCKrylovFinalStateEvaluation mismatched = rotated_evaluation;
    mismatched.policy_hash ^= UINT64_C(1);
    CHECK(mvmc_krylov_final_state_acceptance(
              &evaluation, &mismatched, 0.0, 0.5, &acceptance) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              !acceptance.valid,
          "proposal from a different coefficient policy rejected");
  }

  {
    const double complex cancellation_values[] = {1.0, -1.0, 0.25};
    const double complex cancellation_coefficient[] = {1.0, 1.0};
    MVMCKrylovBoundedResult cancellation =
        raw_krylov_result(cancellation_values, 2);
    MVMCKrylovFinalStatePolicy cancellation_policy;
    CHECK(mvmc_krylov_final_state_policy_create(
              1, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC,
              UINT64_C(0x43414e43454c4c31), cancellation_coefficient, 2,
              &cancellation_policy) == MVMC_KRYLOV_STATUS_OK &&
              mvmc_krylov_final_state_evaluate(
                  &cancellation_policy, &cancellation, &evaluation) ==
                  MVMC_KRYLOV_STATUS_OK &&
              evaluation.valid &&
              evaluation.amplitude.state ==
                  MVMC_SCALED_COMPLEX_NUMERIC_ZERO &&
              !evaluation.sampleable,
          "finite-term cancellation is numeric-zero, not exact-zero");
  }
  {
    MVMCKrylovBoundedResult invalid = krylov;
    invalid.evaluated_order = 2;
    CHECK(mvmc_krylov_final_state_evaluate(
              &policy, &invalid, &evaluation) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              !evaluation.valid,
          "insufficient Krylov order rejected");
    invalid = krylov;
    invalid.value[2].state = MVMC_SCALED_COMPLEX_NONFINITE;
    CHECK(mvmc_krylov_final_state_evaluate(
              &policy, &invalid, &evaluation) ==
              MVMC_KRYLOV_STATUS_NONFINITE &&
              !evaluation.valid,
          "nonfinite basis value rejected");
  }
}

static void test_rank_invariance(void) {
#ifdef _mpi_use
  const double complex values[] = {
      1.0 + 0.2 * I, -0.5 + 0.3 * I, 0.7 - 0.1 * I, 0.2 + 0.9 * I};
  const double complex coefficient[] = {
      0.8 + 0.1 * I, -0.2 + 0.4 * I, 0.3 - 0.5 * I};
  const MVMCKrylovBoundedResult krylov = raw_krylov_result(values, 3);
  MVMCKrylovFinalStatePolicy policy;
  MVMCKrylovFinalStateEvaluation evaluation;
  uint64_t local[9];
  uint64_t minimum[9];
  uint64_t maximum[9];
  double fields[5];
  int index;
  CHECK(mvmc_krylov_final_state_policy_create(
            2, MVMC_KRYLOV_FINAL_COEFFICIENT_SYNTHETIC,
            UINT64_C(0x4d504952414e4b32), coefficient, 3, &policy) ==
            MVMC_KRYLOV_STATUS_OK &&
            mvmc_krylov_final_state_evaluate(
                &policy, &krylov, &evaluation) == MVMC_KRYLOV_STATUS_OK,
        "MPI rank final-state evaluation");
  local[0] = policy.policy_hash;
  local[1] = (uint64_t)(unsigned int)evaluation.amplitude.state;
  local[2] = (uint64_t)(unsigned int)evaluation.local_energy.state;
  local[3] = (uint64_t)(unsigned int)evaluation.sampleable;
  fields[0] = creal(evaluation.amplitude.phase);
  fields[1] = cimag(evaluation.amplitude.phase);
  fields[2] = evaluation.amplitude.log_abs;
  fields[3] = evaluation.local_energy.log_abs;
  fields[4] = evaluation.log_weight;
  for (index = 0; index < 5; ++index) {
    memcpy(&local[index + 4], &fields[index], sizeof(local[index + 4]));
  }
  MPI_Allreduce(local, minimum, 9, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(local, maximum, 9, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  CHECK(memcmp(minimum, maximum, sizeof(minimum)) == 0,
        "P5 policy/evaluation differs by MPI rank");
#endif
}

int main(int argc, char **argv) {
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#else
  (void)argc;
  (void)argv;
#endif
  test_policy_contract();
  test_exact_support_bridge();
  test_synthetic_proposal_transaction();
  test_complex_p2_and_failures();
  test_rank_invariance();

#ifdef _mpi_use
  {
    int global_failures = 0;
    MPI_Allreduce(&failures, &global_failures, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    failures = global_failures;
  }
#endif
  if (failures == 0 && world_rank == 0) {
    puts("krylov final state unit: PASS");
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return failures == 0 ? 0 : 1;
}
