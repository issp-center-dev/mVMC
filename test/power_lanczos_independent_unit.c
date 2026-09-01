#include <math.h>
#include <stdio.h>
#include <string.h>

#include "krylov_gevp_solver.h"
#include "power_lanczos_independent.h"

typedef struct {
    double complex amplitude[4];
} DenseAmplitude;

static int basis_index(uint64_t words) {
    static const uint64_t basis[4] = {UINT64_C(5), UINT64_C(9), UINT64_C(6),
                                      UINT64_C(10)};
    int index;
    for (index = 0; index < 4; ++index)
        if (words == basis[index]) return index;
    return -1;
}

static MVMCKrylovStatus amplitude_callback(
    const uint64_t *words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
    DenseAmplitude *dense = (DenseAmplitude *)context;
    const int index = word_count == 1 ? basis_index(words[0]) : -1;
    if (dense == NULL || result == NULL || index < 0)
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    memset(result, 0, sizeof(*result));
    if (creal(dense->amplitude[index]) == 0.0 &&
        cimag(dense->amplitude[index]) == 0.0)
        return mvmc_scaled_complex_make_exact_zero(&result->value) ==
                       MVMC_PFAFFIAN_STATUS_OK
                   ? MVMC_KRYLOV_STATUS_OK
                   : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    return mvmc_scaled_complex_from_raw(dense->amplitude[index],
                                        &result->value) ==
                   MVMC_PFAFFIAN_STATUS_OK
               ? MVMC_KRYLOV_STATUS_OK
               : MVMC_KRYLOV_STATUS_NONFINITE;
}

static double complex raw(const MVMCScaledComplex *value) {
    double complex result = NAN + I * NAN;
    const MVMCScaledComplexExportStatus status =
        mvmc_scaled_complex_export_common_scale(value, 0.0, &result);
    if (status == MVMC_SCALED_EXPORT_EXACT_ZERO ||
        status == MVMC_SCALED_EXPORT_NUMERIC_ZERO)
        return 0.0;
    return status == MVMC_SCALED_EXPORT_OK ? result : NAN + I * NAN;
}

static int close_complex(double complex actual, double complex expected,
                         double tolerance, const char *label, int index) {
    const double error = cabs(actual - expected);
    const double scale = fmax(1.0, cabs(expected));
    if (!isfinite(error) || error > tolerance * scale) {
        fprintf(stderr,
                "%s[%d] actual=(%.17g,%.17g) expected=(%.17g,%.17g) "
                "error=%.3e\n",
                label, index, creal(actual), cimag(actual), creal(expected),
                cimag(expected), error);
        return 0;
    }
    return 1;
}

static void dense_apply_h(const double complex input[4],
                          double complex output[4]) {
    const double t = 0.5;
    output[0] = 0.7 * input[0] + t * input[1] + t * input[2];
    output[1] = t * input[0] + t * input[3];
    output[2] = t * input[0] + t * input[3];
    output[3] = t * input[1] + t * input[2];
}

static void dense_apply_hprime(const double complex input[4],
                               double complex output[4]) {
    static const double diagonal[4] = {0.2, -0.2, 0.0, -0.4};
    int index;
    for (index = 0; index < 4; ++index)
        output[index] = diagonal[index] * input[index];
}

static int run_oracle(void) {
    int transfer_rows_storage[4][4] = {
        {0, 0, 1, 0}, {1, 0, 0, 0}, {0, 1, 1, 1}, {1, 1, 0, 1}};
    int *transfer_rows[4] = {transfer_rows_storage[0], transfer_rows_storage[1],
                             transfer_rows_storage[2],
                             transfer_rows_storage[3]};
    double complex transfer_parameters[4] = {-0.5, -0.5, -0.5, -0.5};
    int inter_row_storage[1][8] = {{0, 0, 0, 0, 0, 1, 0, 1}};
    int *inter_rows[1] = {inter_row_storage[0]};
    double complex inter_parameters[1] = {0.7};
    int hprime_rows_storage[2][4] = {{0, 0, 0, 0}, {1, 1, 1, 1}};
    int *hprime_rows[2] = {hprime_rows_storage[0], hprime_rows_storage[1]};
    double complex hprime_parameters[2] = {-0.2, 0.4};
    const uint64_t basis[4] = {UINT64_C(5), UINT64_C(9), UINT64_C(6),
                               UINT64_C(10)};
    DenseAmplitude dense = {{1.0, 0.0, 2.0, -0.5}};
    double complex expected_q1[4], expected_r0[4], expected_r1[4];
    double complex expected_h2[4];
    double complex exact_s[3] = {0}, exact_k[3] = {0};
    double complex exact_kr[3] = {0}, exact_b[3] = {0};
    double complex sampled_s[3] = {0}, sampled_k[3] = {0};
    double complex sampled_kr[3] = {0}, sampled_b[3] = {0};
    MVMCPowerLanczosIndependentModel *physical = NULL, *hprime = NULL;
    MVMCPowerLanczosIndependentStreamingModel *physical_view = NULL;
    MVMCKrylovBoundedPlan *plan = NULL;
    MVMCKrylovBoundedWorkspace *workspace = NULL;
    MVMCKrylovBoundedLimits limits;
    MVMCKrylovStatus status;
    int ok = 1;
    int x, row, column, packed;

    dense_apply_hprime(dense.amplitude, expected_q1);
    dense_apply_h(dense.amplitude, expected_r0);
    dense_apply_h(expected_q1, expected_r1);
    dense_apply_h(expected_r0, expected_h2);
    status = mvmc_power_lanczos_independent_model_create(
        2, 1, 1, 4, transfer_rows, transfer_parameters, 1, inter_rows,
        inter_parameters, &physical);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_power_lanczos_independent_streaming_model_create(
            2, 1, 1, 4, transfer_rows, transfer_parameters, 1, inter_rows,
            inter_parameters, &physical_view);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_power_lanczos_independent_model_create(
            2, 1, 1, 2, hprime_rows, hprime_parameters, 0, NULL, NULL, &hprime);
    memset(&limits, 0, sizeof(limits));
    limits.amplitude_policy_hash = UINT64_C(0x494e445054455354);
    limits.cache_bytes = 4096;
    limits.max_row_transitions = 2;
    limits.max_workspace_bytes = 1024 * 1024;
    limits.max_node_expansions = 3;
    limits.max_terminal_amplitude_calls = 2;
    limits.max_total_row_transitions = 2;
    limits.max_order = 1;
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_bounded_krylov_plan_create(
            mvmc_power_lanczos_independent_model_view(hprime), &limits, &plan);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_bounded_krylov_workspace_create(plan, &workspace);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_bounded_krylov_session_begin(
            workspace, amplitude_callback, &dense,
            limits.amplitude_policy_hash);

    for (x = 0; status == MVMC_KRYLOV_STATUS_OK && x < 4; ++x) {
        MVMCPowerLanczosIndependentSnapshot snapshot;
        double complex s[3], k[3], kr[3], b[3];
        double guide_floor_fraction = NAN;
        double guide;
        status = mvmc_power_lanczos_independent_evaluate(
            mvmc_power_lanczos_independent_model_view(physical), workspace,
            basis + x, 1, amplitude_callback, &dense, &snapshot);
        if (status != MVMC_KRYLOV_STATUS_OK) break;
        ok &= close_complex(raw(snapshot.q + 0), dense.amplitude[x], 1e-13,
                            "q0", x);
        ok &=
            close_complex(raw(snapshot.q + 1), expected_q1[x], 1e-13, "q1", x);
        ok &=
            close_complex(raw(snapshot.r + 0), expected_r0[x], 1e-13, "r0", x);
        ok &=
            close_complex(raw(snapshot.r + 1), expected_r1[x], 1e-13, "r1", x);
        if (x == 0) {
            const double complex base_only[2] = {1.0, 0.0};
            MVMCScaledComplex phi, hphi, local;
            double log_weight = NAN;
            status = mvmc_power_lanczos_independent_final_evaluate(
                &snapshot, base_only, &phi, &hphi, &local, &log_weight);
            if (status == MVMC_KRYLOV_STATUS_OK) {
                ok &= close_complex(raw(&phi), dense.amplitude[x], 1e-13,
                                    "rank-one Phi", x);
                ok &= close_complex(raw(&hphi), expected_r0[x], 1e-13,
                                    "rank-one HPhi", x);
                ok &= close_complex(raw(&local),
                                    expected_r0[x] / dense.amplitude[x], 1e-13,
                                    "rank-one local", x);
            }
        }
        if (status == MVMC_KRYLOV_STATUS_OK) {
            MVMCPowerLanczosIndependentSnapshot view_snapshot;
            status = mvmc_power_lanczos_independent_evaluate_view(
                physical_view, workspace, basis + x, 1, amplitude_callback,
                &dense, &view_snapshot);
            if (status == MVMC_KRYLOV_STATUS_OK) {
                ok &= close_complex(raw(view_snapshot.q + 0),
                                    dense.amplitude[x], 1e-13, "view q0", x);
                ok &= close_complex(raw(view_snapshot.q + 1), expected_q1[x],
                                    1e-13, "view q1", x);
                ok &= close_complex(raw(view_snapshot.r + 0), expected_r0[x],
                                    1e-13, "view r0", x);
                ok &= close_complex(raw(view_snapshot.r + 1), expected_r1[x],
                                    1e-13, "view r1", x);
            }
        }
        status = mvmc_power_lanczos_independent_matrix_sample(
            &snapshot, 0x1p-40, s, k, kr, b, &guide_floor_fraction);
        guide = 0x1p-40 + pow(cabs(dense.amplitude[x]), 2.0) +
                pow(cabs(expected_q1[x]), 2.0);
        if (status != MVMC_KRYLOV_STATUS_OK ||
            fabs(guide_floor_fraction - 0x1p-40 / guide) > 1e-13)
            ok = 0;
        packed = 0;
        for (row = 0; row < 2; ++row) {
            const double complex qrow =
                row == 0 ? dense.amplitude[x] : expected_q1[x];
            const double complex rrow =
                row == 0 ? expected_r0[x] : expected_r1[x];
            for (column = row; column < 2; ++column, ++packed) {
                const double complex qcolumn =
                    column == 0 ? dense.amplitude[x] : expected_q1[x];
                const double complex rcolumn =
                    column == 0 ? expected_r0[x] : expected_r1[x];
                sampled_s[packed] += guide * s[packed];
                sampled_k[packed] += guide * k[packed];
                sampled_kr[packed] += guide * kr[packed];
                sampled_b[packed] += guide * b[packed];
                exact_s[packed] += conj(qrow) * qcolumn;
                exact_k[packed] += conj(qrow) * rcolumn;
                exact_kr[packed] += conj(qcolumn) * rrow;
                exact_b[packed] += conj(rrow) * rcolumn;
            }
        }
    }
    for (packed = 0; packed < 3; ++packed) {
        ok &= close_complex(sampled_s[packed], exact_s[packed], 2e-13, "S",
                            packed);
        ok &= close_complex(sampled_k[packed], exact_k[packed], 2e-13, "K",
                            packed);
        ok &= close_complex(sampled_kr[packed], exact_kr[packed], 2e-13, "Kr",
                            packed);
        ok &= close_complex(sampled_b[packed], exact_b[packed], 2e-13, "B",
                            packed);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
        MVMCKrylovGEVPPolicy policy;
        MVMCKrylovGEVPResult gevp;
        const double s00 = creal(exact_s[0]);
        const double s01 = creal(exact_s[1]);
        const double s11 = creal(exact_s[2]);
        const double k00 = 0.5 * creal(exact_k[0] + exact_kr[0]);
        const double k01 = 0.5 * creal(exact_k[1] + exact_kr[1]);
        const double k11 = 0.5 * creal(exact_k[2] + exact_kr[2]);
        const double quadratic = s00 * s11 - s01 * s01;
        const double linear = -(k00 * s11 + k11 * s00 - 2.0 * k01 * s01);
        const double constant = k00 * k11 - k01 * k01;
        const double discriminant =
            linear * linear - 4.0 * quadratic * constant;
        const double dense_energy =
            (-linear - sqrt(discriminant)) / (2.0 * quadratic);
        double complex phi[4], hphi[4];
        double norm = 0.0, first = 0.0, second = 0.0;
        status = mvmc_krylov_gevp_default_policy(1e-12, &policy) ==
                             MVMC_KRYLOV_GEVP_OK &&
                         mvmc_krylov_gevp_solve_complex_packed(
                             &policy, 2, exact_s, exact_k, exact_kr, exact_b, 3,
                             &gevp) == MVMC_KRYLOV_GEVP_OK
                     ? MVMC_KRYLOV_STATUS_OK
                     : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
        if (status == MVMC_KRYLOV_STATUS_OK) {
            ok &= close_complex(gevp.energy, dense_energy, 3e-13, "GEVP energy",
                                0);
            for (x = 0; x < 4; ++x) {
                MVMCPowerLanczosIndependentSnapshot snapshot;
                MVMCScaledComplex scaled_phi, scaled_hphi, local;
                double log_weight = NAN;
                phi[x] = gevp.coefficient[0] * dense.amplitude[x] +
                         gevp.coefficient[1] * expected_q1[x];
                hphi[x] = gevp.coefficient[0] * expected_r0[x] +
                          gevp.coefficient[1] * expected_r1[x];
                status = mvmc_power_lanczos_independent_evaluate_view(
                    physical_view, workspace, basis + x, 1, amplitude_callback,
                    &dense, &snapshot);
                if (status == MVMC_KRYLOV_STATUS_OK)
                    status = mvmc_power_lanczos_independent_final_evaluate(
                        &snapshot, gevp.coefficient, &scaled_phi, &scaled_hphi,
                        &local, &log_weight);
                if (status == MVMC_KRYLOV_STATUS_OK) {
                    ok &= close_complex(raw(&scaled_phi), phi[x], 4e-13,
                                        "final Phi", x);
                    ok &= close_complex(raw(&scaled_hphi), hphi[x], 4e-13,
                                        "final HPhi", x);
                }
                if (status != MVMC_KRYLOV_STATUS_OK) break;
                norm += pow(cabs(phi[x]), 2.0);
                first += creal(conj(phi[x]) * hphi[x]);
                second += pow(cabs(hphi[x]), 2.0);
            }
            if (status == MVMC_KRYLOV_STATUS_OK) {
                const double final_energy = first / norm;
                const double final_variance =
                    second / norm - final_energy * final_energy;
                ok &= close_complex(final_energy, dense_energy, 4e-13,
                                    "final energy", 0);
                ok &= close_complex(gevp.energy_squared, second / norm, 4e-13,
                                    "final second moment", 0);
                ok &= close_complex(gevp.variance, final_variance, 4e-13,
                                    "final variance", 0);
            }
            if (status == MVMC_KRYLOV_STATUS_OK) {
                const double amplitude_scale = 1e150;
                MVMCPowerLanczosIndependentSnapshot scaled_snapshot;
                MVMCScaledComplex scaled_phi, scaled_hphi, scaled_local;
                double scaled_log_weight = NAN;
                status = mvmc_bounded_krylov_session_end(workspace);
                for (x = 0; x < 4; ++x) dense.amplitude[x] *= amplitude_scale;
                if (status == MVMC_KRYLOV_STATUS_OK)
                    status = mvmc_bounded_krylov_session_begin(
                        workspace, amplitude_callback, &dense,
                        limits.amplitude_policy_hash);
                if (status == MVMC_KRYLOV_STATUS_OK)
                    status = mvmc_power_lanczos_independent_evaluate_view(
                        physical_view, workspace, basis, 1, amplitude_callback,
                        &dense, &scaled_snapshot);
                if (status == MVMC_KRYLOV_STATUS_OK)
                    status = mvmc_power_lanczos_independent_final_evaluate(
                        &scaled_snapshot, gevp.coefficient, &scaled_phi,
                        &scaled_hphi, &scaled_local, &scaled_log_weight);
                if (status == MVMC_KRYLOV_STATUS_OK) {
                    ok &= close_complex(raw(&scaled_local), hphi[0] / phi[0],
                                        8e-13, "scaled local energy", 0);
                    ok &= fabs(scaled_log_weight - 2.0 * (log(amplitude_scale) +
                                                          log(cabs(phi[0])))) <
                          8e-13;
                }
                for (x = 0; x < 4; ++x) dense.amplitude[x] /= amplitude_scale;
            }
        }
    }
    if (workspace != NULL && mvmc_bounded_krylov_session_is_active(workspace))
        (void)mvmc_bounded_krylov_session_end(workspace);
    mvmc_bounded_krylov_workspace_destroy(workspace);
    mvmc_bounded_krylov_plan_destroy(plan);
    workspace = NULL;
    plan = NULL;
    limits.cache_bytes = 0;
    limits.max_row_transitions = 5;
    limits.max_node_expansions = 6;
    limits.max_terminal_amplitude_calls = 5;
    limits.max_total_row_transitions = 5;
    status = mvmc_bounded_krylov_plan_create(
        mvmc_power_lanczos_independent_model_view(physical), &limits, &plan);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_bounded_krylov_workspace_create(plan, &workspace);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_bounded_krylov_session_begin(
            workspace, amplitude_callback, &dense,
            limits.amplitude_policy_hash);
    for (x = 0; status == MVMC_KRYLOV_STATUS_OK && x < 4; ++x) {
        MVMCPowerLanczosIndependentSnapshot equal_snapshot;
        status = mvmc_power_lanczos_independent_evaluate(
            mvmc_power_lanczos_independent_model_view(physical), workspace,
            basis + x, 1, amplitude_callback, &dense, &equal_snapshot);
        if (status == MVMC_KRYLOV_STATUS_OK) {
            ok &= close_complex(raw(equal_snapshot.q + 1), expected_r0[x],
                                2e-13, "H'=H q1", x);
            ok &= close_complex(raw(equal_snapshot.r + 1), expected_h2[x],
                                3e-13, "H'=H r1", x);
        }
    }
    if (workspace != NULL && mvmc_bounded_krylov_session_is_active(workspace))
        (void)mvmc_bounded_krylov_session_end(workspace);
    mvmc_bounded_krylov_workspace_destroy(workspace);
    mvmc_bounded_krylov_plan_destroy(plan);
    mvmc_power_lanczos_independent_model_destroy(hprime);
    mvmc_power_lanczos_independent_model_destroy(physical);
    mvmc_power_lanczos_independent_streaming_model_destroy(physical_view);
    if (status != MVMC_KRYLOV_STATUS_OK) {
        fprintf(stderr, "oracle status=%s\n",
                mvmc_krylov_status_string(status));
        ok = 0;
    }
    return ok;
}

static int run_negative(void) {
    int bad_storage[1][4] = {{0, 0, 1, 1}};
    int *bad_rows[1] = {bad_storage[0]};
    double complex parameter[1] = {-1.0};
    MVMCPowerLanczosIndependentModel *model = NULL;
    return mvmc_power_lanczos_independent_model_create(
               2, 1, 1, 1, bad_rows, parameter, 0, NULL, NULL, &model) ==
               MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
           model == NULL;
}

int main(void) {
    if (!run_oracle() || !run_negative()) return 1;
    puts("power_lanczos_independent_unit: PASS");
    return 0;
}
