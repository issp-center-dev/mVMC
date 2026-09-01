#include "power_lanczos_independent.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "power_lanczos_independent.c requires the bounded Krylov engine"
#endif

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct MVMCPowerLanczosIndependentModel {
    MVMCKrylovFockModel model;
    MVMCKrylovHamiltonianTerm *terms;
    MVMCKrylovFermionOperator *operators;
};

struct MVMCPowerLanczosIndependentStreamingModel {
    size_t site_count;
    size_t up_electron_count;
    size_t down_electron_count;
    size_t transfer_count;
    int *const *transfer_indices;
    const double complex *transfer_parameters;
    size_t inter_all_count;
    int *const *inter_all_indices;
    const double complex *inter_all_parameters;
};

typedef struct {
    MVMCKrylovBoundedWorkspace *workspace;
    MVMCKrylovScaledAmplitudeCallback amplitude;
    void *amplitude_context;
} HPrimeAmplitudeContext;

typedef struct {
    int empty;
    double log_scale;
    double log_error_bound;
    double max_input_log_abs;
    double complex sum;
    double complex compensation;
} ScaledAccumulator;

static double log_add(double left, double right);

static int finite_complex(double complex value) {
    return isfinite(creal(value)) && isfinite(cimag(value));
}

static int checked_add(size_t left, size_t right, size_t *sum) {
    if (sum == NULL || right > SIZE_MAX - left) return 0;
    *sum = left + right;
    return 1;
}

static int checked_multiply(size_t left, size_t right, size_t *product) {
    if (product == NULL || (left != 0 && right > SIZE_MAX / left)) return 0;
    *product = left * right;
    return 1;
}

static int valid_transfer_row(int site_count, const int *row) {
    return row != NULL && row[0] >= 0 && row[0] < site_count && row[1] >= 0 &&
           row[1] <= 1 && row[2] >= 0 && row[2] < site_count && row[3] >= 0 &&
           row[3] <= 1 && row[1] == row[3];
}

static int valid_inter_all_row(int site_count, const int *row) {
    return row != NULL && row[0] >= 0 && row[0] < site_count && row[1] >= 0 &&
           row[1] <= 1 && row[2] >= 0 && row[2] < site_count && row[3] >= 0 &&
           row[3] <= 1 && row[4] >= 0 && row[4] < site_count && row[5] >= 0 &&
           row[5] <= 1 && row[6] >= 0 && row[6] < site_count && row[7] >= 0 &&
           row[7] <= 1 && row[1] == row[3] && row[5] == row[7];
}

static int bit_value(const uint64_t *words, size_t orbital) {
    return (int)((words[orbital / 64] >> (orbital % 64)) & UINT64_C(1));
}

static void set_bit(uint64_t *words, size_t orbital, int value) {
    const uint64_t mask = UINT64_C(1) << (orbital % 64);
    if (value)
        words[orbital / 64] |= mask;
    else
        words[orbital / 64] &= ~mask;
}

static uint64_t occupied_before(const uint64_t *words, size_t orbital) {
    const size_t full_words = orbital / 64;
    const size_t partial_bits = orbital % 64;
    uint64_t count = 0;
    size_t index;
    for (index = 0; index < full_words; ++index) {
#if defined(__GNUC__) || defined(__clang__)
        count += (uint64_t)__builtin_popcountll(words[index]);
#else
        uint64_t value = words[index];
        while (value != 0) {
            value &= value - UINT64_C(1);
            ++count;
        }
#endif
    }
    if (partial_bits != 0) {
        const uint64_t mask = (UINT64_C(1) << partial_bits) - UINT64_C(1);
#if defined(__GNUC__) || defined(__clang__)
        count += (uint64_t)__builtin_popcountll(words[full_words] & mask);
#else
        uint64_t value = words[full_words] & mask;
        while (value != 0) {
            value &= value - UINT64_C(1);
            ++count;
        }
#endif
    }
    return count;
}

static MVMCKrylovStatus apply_term(const MVMCKrylovFockModel *model,
                                   const MVMCKrylovHamiltonianTerm *term,
                                   const uint64_t *input, size_t word_count,
                                   uint64_t *output, int *applied, int *sign) {
    size_t reverse;
    int phase = 1;
    memcpy(output, input, word_count * sizeof(*output));
    *applied = 0;
    *sign = 1;
    for (reverse = term->operator_count; reverse > 0; --reverse) {
        const MVMCKrylovFermionOperator *item =
            model->operators + term->operator_offset + reverse - 1;
        const int occupied = bit_value(output, item->orbital);
        if ((item->kind == MVMC_KRYLOV_FERMION_CREATE && occupied) ||
            (item->kind == MVMC_KRYLOV_FERMION_ANNIHILATE && !occupied)) {
            return MVMC_KRYLOV_STATUS_OK;
        }
        if ((occupied_before(output, item->orbital) & UINT64_C(1)) != 0)
            phase = -phase;
        set_bit(output, item->orbital,
                item->kind == MVMC_KRYLOV_FERMION_CREATE);
    }
    *applied = 1;
    *sign = phase;
    return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus accumulator_add(ScaledAccumulator *accumulator,
                                        const MVMCScaledComplex *value) {
    double complex term;
    double complex updated;
    if (accumulator == NULL || !mvmc_scaled_complex_is_valid(value) ||
        value->state == MVMC_SCALED_COMPLEX_NONFINITE)
        return MVMC_KRYLOV_STATUS_NONFINITE;
    accumulator->log_error_bound =
        log_add(accumulator->log_error_bound, value->log_abs_error_bound);
    accumulator->max_input_log_abs =
        fmax(accumulator->max_input_log_abs, value->max_input_log_abs);
    if (value->state != MVMC_SCALED_COMPLEX_FINITE_NONZERO)
        return MVMC_KRYLOV_STATUS_OK;
    if (accumulator->empty) {
        accumulator->empty = 0;
        accumulator->log_scale = value->log_abs;
        accumulator->sum = value->phase;
        accumulator->compensation = 0.0;
        return MVMC_KRYLOV_STATUS_OK;
    }
    if (value->log_abs > accumulator->log_scale) {
        const double factor = exp(accumulator->log_scale - value->log_abs);
        if (!isfinite(factor)) return MVMC_KRYLOV_STATUS_NONFINITE;
        accumulator->sum *= factor;
        accumulator->compensation *= factor;
        accumulator->log_scale = value->log_abs;
    }
    term = value->phase * exp(value->log_abs - accumulator->log_scale);
    updated = accumulator->sum + term;
    if (cabs(accumulator->sum) >= cabs(term))
        accumulator->compensation += (accumulator->sum - updated) + term;
    else
        accumulator->compensation += (term - updated) + accumulator->sum;
    accumulator->sum = updated;
    return finite_complex(accumulator->sum) &&
                   finite_complex(accumulator->compensation)
               ? MVMC_KRYLOV_STATUS_OK
               : MVMC_KRYLOV_STATUS_NONFINITE;
}

static MVMCKrylovStatus accumulator_finish(const ScaledAccumulator *accumulator,
                                           MVMCScaledComplex *result) {
    const double complex raw = accumulator->sum + accumulator->compensation;
    const double magnitude = cabs(raw);
    if (accumulator->empty || magnitude == 0.0) {
        if (isfinite(accumulator->log_error_bound))
            return mvmc_scaled_complex_make_numeric_zero(
                       accumulator->log_error_bound,
                       accumulator->max_input_log_abs, -INFINITY, 0.0,
                       result) == MVMC_PFAFFIAN_STATUS_OK
                       ? MVMC_KRYLOV_STATUS_OK
                       : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
        return mvmc_scaled_complex_make_exact_zero(result) ==
                       MVMC_PFAFFIAN_STATUS_OK
                   ? MVMC_KRYLOV_STATUS_OK
                   : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    if (!finite_complex(raw) || !isfinite(magnitude) ||
        mvmc_scaled_complex_make_finite(
            raw / magnitude, accumulator->log_scale + log(magnitude),
            log_add(accumulator->log_error_bound, accumulator->log_scale +
                                                      log(magnitude) +
                                                      log(128.0 * DBL_EPSILON)),
            result) != MVMC_PFAFFIAN_STATUS_OK)
        return MVMC_KRYLOV_STATUS_NONFINITE;
    return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus append_term(MVMCPowerLanczosIndependentModel *owned,
                                    size_t term_index, size_t operator_offset,
                                    double complex coefficient,
                                    const size_t *orbitals,
                                    size_t orbital_count) {
    size_t index;
    MVMCKrylovHamiltonianTerm *term = owned->terms + term_index;
    if (!finite_complex(coefficient) ||
        (orbital_count != 2 && orbital_count != 4))
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    term->coefficient = coefficient;
    term->operator_offset = operator_offset;
    term->operator_count = orbital_count;
    term->source_kind = orbital_count == 2 ? 100 : 101;
    term->source_index = term_index;
    for (index = 0; index < orbital_count; ++index) {
        owned->operators[operator_offset + index].kind =
            index % 2 == 0 ? MVMC_KRYLOV_FERMION_CREATE
                           : MVMC_KRYLOV_FERMION_ANNIHILATE;
        owned->operators[operator_offset + index].orbital = orbitals[index];
    }
    return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_independent_model_create(
    int site_count, int up_electron_count, int down_electron_count,
    size_t transfer_count, int *const *transfer_indices,
    const double complex *transfer_parameters, size_t inter_all_count,
    int *const *inter_all_indices, const double complex *inter_all_parameters,
    MVMCPowerLanczosIndependentModel **result) {
    MVMCPowerLanczosIndependentModel *owned = NULL;
    size_t term_count = 0, operator_count = 0;
    size_t term_bytes = 0, operator_bytes = 0;
    size_t term_index = 0, operator_index = 0, index;
    MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
    if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    *result = NULL;
    if (site_count <= 0 || up_electron_count < 0 || down_electron_count < 0 ||
        up_electron_count > site_count || down_electron_count > site_count ||
        (transfer_count != 0 &&
         (transfer_indices == NULL || transfer_parameters == NULL)) ||
        (inter_all_count != 0 &&
         (inter_all_indices == NULL || inter_all_parameters == NULL)) ||
        !checked_add(transfer_count, inter_all_count, &term_count) ||
        !checked_multiply(transfer_count, 2, &operator_count) ||
        inter_all_count > (SIZE_MAX - operator_count) / 4)
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    operator_count += 4 * inter_all_count;
    if (!checked_multiply(term_count, sizeof(*owned->terms), &term_bytes) ||
        !checked_multiply(operator_count, sizeof(*owned->operators),
                          &operator_bytes))
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    owned = (MVMCPowerLanczosIndependentModel *)calloc(1, sizeof(*owned));
    if (owned != NULL && term_count != 0)
        owned->terms = (MVMCKrylovHamiltonianTerm *)calloc(1, term_bytes);
    if (owned != NULL && operator_count != 0)
        owned->operators =
            (MVMCKrylovFermionOperator *)calloc(1, operator_bytes);
    if (owned == NULL || (term_count != 0 && owned->terms == NULL) ||
        (operator_count != 0 && owned->operators == NULL)) {
        mvmc_power_lanczos_independent_model_destroy(owned);
        return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    }
    for (index = 0; status == MVMC_KRYLOV_STATUS_OK && index < transfer_count;
         ++index) {
        const int *row = transfer_indices[index];
        size_t orbitals[2];
        if (!valid_transfer_row(site_count, row)) {
            status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
            break;
        }
        orbitals[0] = (size_t)row[0] + (size_t)row[1] * (size_t)site_count;
        orbitals[1] = (size_t)row[2] + (size_t)row[3] * (size_t)site_count;
        status = append_term(owned, term_index++, operator_index,
                             -transfer_parameters[index], orbitals, 2);
        operator_index += 2;
    }
    for (index = 0; status == MVMC_KRYLOV_STATUS_OK && index < inter_all_count;
         ++index) {
        const int *row = inter_all_indices[index];
        size_t orbitals[4];
        if (!valid_inter_all_row(site_count, row)) {
            status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
            break;
        }
        orbitals[0] = (size_t)row[0] + (size_t)row[1] * (size_t)site_count;
        orbitals[1] = (size_t)row[2] + (size_t)row[3] * (size_t)site_count;
        orbitals[2] = (size_t)row[4] + (size_t)row[5] * (size_t)site_count;
        orbitals[3] = (size_t)row[6] + (size_t)row[7] * (size_t)site_count;
        status = append_term(owned, term_index++, operator_index,
                             inter_all_parameters[index], orbitals, 4);
        operator_index += 4;
    }
    if (status != MVMC_KRYLOV_STATUS_OK || term_index != term_count ||
        operator_index != operator_count) {
        mvmc_power_lanczos_independent_model_destroy(owned);
        return status == MVMC_KRYLOV_STATUS_OK
                   ? MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE
                   : status;
    }
    owned->model.site_count = (size_t)site_count;
    owned->model.up_electron_count = (size_t)up_electron_count;
    owned->model.down_electron_count = (size_t)down_electron_count;
    owned->model.pure_spin = 0;
    owned->model.hermitian = 1;
    owned->model.terms = owned->terms;
    owned->model.term_count = term_count;
    owned->model.operators = owned->operators;
    owned->model.operator_count = operator_count;
    *result = owned;
    return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_power_lanczos_independent_model_destroy(
    MVMCPowerLanczosIndependentModel *model) {
    if (model == NULL) return;
    free(model->operators);
    free(model->terms);
    free(model);
}

const MVMCKrylovFockModel *mvmc_power_lanczos_independent_model_view(
    const MVMCPowerLanczosIndependentModel *model) {
    return model == NULL ? NULL : &model->model;
}

MVMCKrylovStatus mvmc_power_lanczos_independent_streaming_model_create(
    int site_count, int up_electron_count, int down_electron_count,
    size_t transfer_count, int *const *transfer_indices,
    const double complex *transfer_parameters, size_t inter_all_count,
    int *const *inter_all_indices, const double complex *inter_all_parameters,
    MVMCPowerLanczosIndependentStreamingModel **result) {
    MVMCPowerLanczosIndependentStreamingModel *created;
    size_t index;
    if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    *result = NULL;
    if (site_count <= 0 || up_electron_count < 0 || down_electron_count < 0 ||
        up_electron_count > site_count || down_electron_count > site_count ||
        (transfer_count != 0 &&
         (transfer_indices == NULL || transfer_parameters == NULL)) ||
        (inter_all_count != 0 &&
         (inter_all_indices == NULL || inter_all_parameters == NULL)))
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    for (index = 0; index < transfer_count; ++index)
        if (!valid_transfer_row(site_count, transfer_indices[index]) ||
            !finite_complex(transfer_parameters[index]))
            return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    for (index = 0; index < inter_all_count; ++index)
        if (!valid_inter_all_row(site_count, inter_all_indices[index]) ||
            !finite_complex(inter_all_parameters[index]))
            return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    created = (MVMCPowerLanczosIndependentStreamingModel *)calloc(
        1, sizeof(*created));
    if (created == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    created->site_count = (size_t)site_count;
    created->up_electron_count = (size_t)up_electron_count;
    created->down_electron_count = (size_t)down_electron_count;
    created->transfer_count = transfer_count;
    created->transfer_indices = transfer_indices;
    created->transfer_parameters = transfer_parameters;
    created->inter_all_count = inter_all_count;
    created->inter_all_indices = inter_all_indices;
    created->inter_all_parameters = inter_all_parameters;
    *result = created;
    return MVMC_KRYLOV_STATUS_OK;
}

void mvmc_power_lanczos_independent_streaming_model_destroy(
    MVMCPowerLanczosIndependentStreamingModel *model) {
    free(model);
}

MVMCKrylovStatus mvmc_power_lanczos_independent_stream_apply(
    const MVMCKrylovFockModel *model, const uint64_t *configuration_words,
    size_t word_count, MVMCKrylovScaledAmplitudeCallback child_amplitude,
    void *child_context, MVMCScaledComplex *result, uint64_t *applied_terms,
    uint64_t *amplitude_calls) {
    uint64_t *neighbor = NULL;
    ScaledAccumulator accumulator;
    MVMCKrylovStatus status;
    size_t index;
    if (result == NULL || model == NULL || configuration_words == NULL ||
        child_amplitude == NULL || word_count == 0 ||
        word_count > SIZE_MAX / sizeof(*neighbor))
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    status = mvmc_krylov_fock_validate(model, configuration_words, word_count);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    neighbor = (uint64_t *)malloc(word_count * sizeof(*neighbor));
    if (neighbor == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    memset(&accumulator, 0, sizeof(accumulator));
    accumulator.empty = 1;
    accumulator.log_error_bound = -INFINITY;
    accumulator.max_input_log_abs = -INFINITY;
    if (applied_terms != NULL) *applied_terms = 0;
    if (amplitude_calls != NULL) *amplitude_calls = 0;
    for (index = 0; index < model->term_count; ++index) {
        const MVMCKrylovHamiltonianTerm *term = model->terms + index;
        MVMCKrylovScaledAmplitudeResult child;
        MVMCScaledComplex weight, contribution;
        int applied = 0, sign = 1;
        if (creal(term->coefficient) == 0.0 && cimag(term->coefficient) == 0.0)
            continue;
        status = apply_term(model, term, configuration_words, word_count,
                            neighbor, &applied, &sign);
        if (status != MVMC_KRYLOV_STATUS_OK || !applied) {
            if (status != MVMC_KRYLOV_STATUS_OK) break;
            continue;
        }
        memset(&child, 0, sizeof(child));
        status = child_amplitude(neighbor, word_count, child_context, &child);
        if (status != MVMC_KRYLOV_STATUS_OK) break;
        if (mvmc_scaled_complex_from_raw(conj(term->coefficient * (double)sign),
                                         &weight) != MVMC_PFAFFIAN_STATUS_OK ||
            mvmc_scaled_complex_multiply(&weight, &child.value,
                                         &contribution) !=
                MVMC_PFAFFIAN_STATUS_OK) {
            status = MVMC_KRYLOV_STATUS_NONFINITE;
            break;
        }
        status = accumulator_add(&accumulator, &contribution);
        if (status != MVMC_KRYLOV_STATUS_OK) break;
        if ((applied_terms != NULL && *applied_terms == UINT64_MAX) ||
            (amplitude_calls != NULL && *amplitude_calls == UINT64_MAX)) {
            status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
            break;
        }
        if (applied_terms != NULL) ++*applied_terms;
        if (amplitude_calls != NULL) ++*amplitude_calls;
    }
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = accumulator_finish(&accumulator, result);
    free(neighbor);
    return status;
}

static MVMCKrylovStatus validate_streaming_configuration(
    const MVMCPowerLanczosIndependentStreamingModel *model,
    const uint64_t *words, size_t word_count) {
    const size_t orbital_count = 2 * model->site_count;
    const size_t used_bits = orbital_count % 64;
    size_t site, up_count = 0, down_count = 0;
    if (words == NULL ||
        word_count != mvmc_krylov_fock_word_count(model->site_count))
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    if (used_bits != 0) {
        const uint64_t used_mask = (UINT64_C(1) << used_bits) - UINT64_C(1);
        if ((words[word_count - 1] & ~used_mask) != 0)
            return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    for (site = 0; site < model->site_count; ++site) {
        up_count += (size_t)bit_value(words, site);
        down_count += (size_t)bit_value(words, site + model->site_count);
    }
    return up_count == model->up_electron_count &&
                   down_count == model->down_electron_count
               ? MVMC_KRYLOV_STATUS_OK
               : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
}

static MVMCKrylovStatus stream_one_term(
    const MVMCKrylovFockModel *scratch_model,
    const MVMCKrylovHamiltonianTerm *term, const uint64_t *configuration_words,
    size_t word_count, uint64_t *neighbor,
    MVMCKrylovScaledAmplitudeCallback child_amplitude, void *child_context,
    ScaledAccumulator *accumulator, uint64_t *applied_terms,
    uint64_t *amplitude_calls) {
    MVMCKrylovScaledAmplitudeResult child;
    MVMCScaledComplex weight, contribution;
    MVMCKrylovStatus status;
    int applied = 0, sign = 1;
    if (creal(term->coefficient) == 0.0 && cimag(term->coefficient) == 0.0)
        return MVMC_KRYLOV_STATUS_OK;
    status = apply_term(scratch_model, term, configuration_words, word_count,
                        neighbor, &applied, &sign);
    if (status != MVMC_KRYLOV_STATUS_OK || !applied) return status;
    memset(&child, 0, sizeof(child));
    status = child_amplitude(neighbor, word_count, child_context, &child);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    if (mvmc_scaled_complex_from_raw(conj(term->coefficient * (double)sign),
                                     &weight) != MVMC_PFAFFIAN_STATUS_OK ||
        mvmc_scaled_complex_multiply(&weight, &child.value, &contribution) !=
            MVMC_PFAFFIAN_STATUS_OK)
        return MVMC_KRYLOV_STATUS_NONFINITE;
    status = accumulator_add(accumulator, &contribution);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    if (applied_terms != NULL) {
        if (*applied_terms == UINT64_MAX)
            return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
        ++*applied_terms;
    }
    if (amplitude_calls != NULL) {
        if (*amplitude_calls == UINT64_MAX)
            return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
        ++*amplitude_calls;
    }
    return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_independent_stream_apply_view(
    const MVMCPowerLanczosIndependentStreamingModel *model,
    const uint64_t *configuration_words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback child_amplitude, void *child_context,
    MVMCScaledComplex *result, uint64_t *applied_terms,
    uint64_t *amplitude_calls) {
    MVMCKrylovFockModel scratch_model;
    MVMCKrylovHamiltonianTerm term;
    MVMCKrylovFermionOperator operators[4];
    ScaledAccumulator accumulator;
    uint64_t *neighbor = NULL;
    MVMCKrylovStatus status;
    size_t index, operator_index;
    if (result == NULL || model == NULL || configuration_words == NULL ||
        child_amplitude == NULL || word_count == 0 ||
        word_count > SIZE_MAX / sizeof(*neighbor))
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    status = validate_streaming_configuration(model, configuration_words,
                                              word_count);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    neighbor = (uint64_t *)malloc(word_count * sizeof(*neighbor));
    if (neighbor == NULL) return MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    memset(&scratch_model, 0, sizeof(scratch_model));
    scratch_model.operators = operators;
    memset(&term, 0, sizeof(term));
    term.operator_offset = 0;
    memset(&accumulator, 0, sizeof(accumulator));
    accumulator.empty = 1;
    accumulator.log_error_bound = -INFINITY;
    accumulator.max_input_log_abs = -INFINITY;
    if (applied_terms != NULL) *applied_terms = 0;
    if (amplitude_calls != NULL) *amplitude_calls = 0;
    for (index = 0;
         status == MVMC_KRYLOV_STATUS_OK && index < model->transfer_count;
         ++index) {
        const int *row = model->transfer_indices[index];
        term.coefficient = -model->transfer_parameters[index];
        term.operator_count = 2;
        for (operator_index = 0; operator_index < 2; ++operator_index) {
            operators[operator_index].kind =
                operator_index == 0 ? MVMC_KRYLOV_FERMION_CREATE
                                    : MVMC_KRYLOV_FERMION_ANNIHILATE;
            operators[operator_index].orbital =
                (size_t)row[2 * operator_index] +
                (size_t)row[2 * operator_index + 1] * model->site_count;
        }
        status = stream_one_term(&scratch_model, &term, configuration_words,
                                 word_count, neighbor, child_amplitude,
                                 child_context, &accumulator, applied_terms,
                                 amplitude_calls);
    }
    for (index = 0;
         status == MVMC_KRYLOV_STATUS_OK && index < model->inter_all_count;
         ++index) {
        const int *row = model->inter_all_indices[index];
        term.coefficient = model->inter_all_parameters[index];
        term.operator_count = 4;
        for (operator_index = 0; operator_index < 4; ++operator_index) {
            operators[operator_index].kind =
                operator_index % 2 == 0 ? MVMC_KRYLOV_FERMION_CREATE
                                        : MVMC_KRYLOV_FERMION_ANNIHILATE;
            operators[operator_index].orbital =
                (size_t)row[2 * operator_index] +
                (size_t)row[2 * operator_index + 1] * model->site_count;
        }
        status = stream_one_term(&scratch_model, &term, configuration_words,
                                 word_count, neighbor, child_amplitude,
                                 child_context, &accumulator, applied_terms,
                                 amplitude_calls);
    }
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = accumulator_finish(&accumulator, result);
    free(neighbor);
    return status;
}

static MVMCKrylovStatus hprime_amplitude(
    const uint64_t *configuration_words, size_t word_count, void *context,
    MVMCKrylovScaledAmplitudeResult *result) {
    HPrimeAmplitudeContext *hprime = (HPrimeAmplitudeContext *)context;
    MVMCKrylovBoundedResult evaluated;
    MVMCKrylovStatus status;
    if (hprime == NULL || result == NULL)
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    memset(result, 0, sizeof(*result));
    status = mvmc_bounded_krylov_session_evaluate_bound(
        hprime->workspace, configuration_words, word_count, hprime->amplitude,
        hprime->amplitude_context, &evaluated);
    if (status == MVMC_KRYLOV_STATUS_OK &&
        (!evaluated.valid || evaluated.evaluated_order < 1))
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    if (status == MVMC_KRYLOV_STATUS_OK) result->value = evaluated.value[1];
    return status;
}

typedef MVMCKrylovStatus (*IndependentStreamApply)(
    const void *physical_model, const uint64_t *configuration_words,
    size_t word_count, MVMCKrylovScaledAmplitudeCallback child_amplitude,
    void *child_context, MVMCScaledComplex *result, uint64_t *applied_terms,
    uint64_t *amplitude_calls);

static MVMCKrylovStatus stream_owned_adapter(
    const void *physical_model, const uint64_t *configuration_words,
    size_t word_count, MVMCKrylovScaledAmplitudeCallback child_amplitude,
    void *child_context, MVMCScaledComplex *result, uint64_t *applied_terms,
    uint64_t *amplitude_calls) {
    return mvmc_power_lanczos_independent_stream_apply(
        (const MVMCKrylovFockModel *)physical_model, configuration_words,
        word_count, child_amplitude, child_context, result, applied_terms,
        amplitude_calls);
}

static MVMCKrylovStatus stream_view_adapter(
    const void *physical_model, const uint64_t *configuration_words,
    size_t word_count, MVMCKrylovScaledAmplitudeCallback child_amplitude,
    void *child_context, MVMCScaledComplex *result, uint64_t *applied_terms,
    uint64_t *amplitude_calls) {
    return mvmc_power_lanczos_independent_stream_apply_view(
        (const MVMCPowerLanczosIndependentStreamingModel *)physical_model,
        configuration_words, word_count, child_amplitude, child_context, result,
        applied_terms, amplitude_calls);
}

static MVMCKrylovStatus evaluate_with_stream(
    const void *physical_model, IndependentStreamApply stream_apply,
    MVMCKrylovBoundedWorkspace *hprime_workspace,
    const uint64_t *configuration_words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCPowerLanczosIndependentSnapshot *result) {
    MVMCKrylovBoundedResult q;
    HPrimeAmplitudeContext hprime;
    uint64_t terms = 0, calls = 0;
    MVMCKrylovStatus status;
    if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    memset(result, 0, sizeof(*result));
    if (physical_model == NULL || stream_apply == NULL ||
        hprime_workspace == NULL || configuration_words == NULL ||
        amplitude == NULL ||
        !mvmc_bounded_krylov_session_is_active(hprime_workspace))
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    status = mvmc_bounded_krylov_session_evaluate_bound(
        hprime_workspace, configuration_words, word_count, amplitude,
        amplitude_context, &q);
    if (status == MVMC_KRYLOV_STATUS_OK && (!q.valid || q.evaluated_order < 1))
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    result->q[0] = q.value[0];
    result->q[1] = q.value[1];
    status =
        stream_apply(physical_model, configuration_words, word_count, amplitude,
                     amplitude_context, &result->r[0], &terms, &calls);
    hprime.workspace = hprime_workspace;
    hprime.amplitude = amplitude;
    hprime.amplitude_context = amplitude_context;
    if (status == MVMC_KRYLOV_STATUS_OK) {
        uint64_t second_terms = 0, second_calls = 0;
        status = stream_apply(physical_model, configuration_words, word_count,
                              hprime_amplitude, &hprime, &result->r[1],
                              &second_terms, &second_calls);
        if (UINT64_MAX - terms < second_terms ||
            UINT64_MAX - calls < second_calls)
            status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
        else {
            terms += second_terms;
            calls += second_calls;
        }
    }
    result->physical_row_terms = terms;
    result->amplitude_calls = calls;
    return status;
}

MVMCKrylovStatus mvmc_power_lanczos_independent_evaluate(
    const MVMCKrylovFockModel *physical_model,
    MVMCKrylovBoundedWorkspace *hprime_workspace,
    const uint64_t *configuration_words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCPowerLanczosIndependentSnapshot *result) {
    return evaluate_with_stream(
        physical_model, stream_owned_adapter, hprime_workspace,
        configuration_words, word_count, amplitude, amplitude_context, result);
}

MVMCKrylovStatus mvmc_power_lanczos_independent_evaluate_view(
    const MVMCPowerLanczosIndependentStreamingModel *physical_model,
    MVMCKrylovBoundedWorkspace *hprime_workspace,
    const uint64_t *configuration_words, size_t word_count,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    MVMCPowerLanczosIndependentSnapshot *result) {
    return evaluate_with_stream(
        physical_model, stream_view_adapter, hprime_workspace,
        configuration_words, word_count, amplitude, amplitude_context, result);
}

static MVMCKrylovStatus normalized_product(const MVMCScaledComplex *left,
                                           const MVMCScaledComplex *right,
                                           double log_guide,
                                           double complex *result) {
    double magnitude;
    if (!mvmc_scaled_complex_is_valid(left) ||
        !mvmc_scaled_complex_is_valid(right) || result == NULL ||
        !isfinite(log_guide))
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    if (left->state == MVMC_SCALED_COMPLEX_NONFINITE ||
        right->state == MVMC_SCALED_COMPLEX_NONFINITE)
        return MVMC_KRYLOV_STATUS_NONFINITE;
    if (left->state != MVMC_SCALED_COMPLEX_FINITE_NONZERO ||
        right->state != MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
        *result = 0.0;
        return MVMC_KRYLOV_STATUS_OK;
    }
    magnitude = exp(left->log_abs + right->log_abs - log_guide);
    *result = magnitude * conj(left->phase) * right->phase;
    return isfinite(magnitude) && finite_complex(*result)
               ? MVMC_KRYLOV_STATUS_OK
               : MVMC_KRYLOV_STATUS_NONFINITE;
}

static double log_add(double left, double right) {
    double maximum;
    if (isinf(left) && left < 0.0) return right;
    if (isinf(right) && right < 0.0) return left;
    maximum = fmax(left, right);
    return maximum + log1p(exp(fmin(left, right) - maximum));
}

MVMCKrylovStatus mvmc_power_lanczos_independent_matrix_sample(
    const MVMCPowerLanczosIndependentSnapshot *snapshot, double eta,
    double complex overlap_upper[3],
    double complex hamiltonian_forward_upper[3],
    double complex hamiltonian_reverse_upper[3],
    double complex hamiltonian_squared_upper[3], double *guide_denominator) {
    double log_guide;
    size_t row, column, packed = 0;
    MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
    if (snapshot == NULL || overlap_upper == NULL ||
        hamiltonian_forward_upper == NULL ||
        hamiltonian_reverse_upper == NULL ||
        hamiltonian_squared_upper == NULL || guide_denominator == NULL ||
        !isfinite(eta) || eta <= 0.0)
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    log_guide = log(eta);
    for (row = 0; row < 2; ++row) {
        if (!mvmc_scaled_complex_is_valid(snapshot->q + row) ||
            !mvmc_scaled_complex_is_valid(snapshot->r + row))
            return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
        if (snapshot->q[row].state == MVMC_SCALED_COMPLEX_FINITE_NONZERO)
            log_guide = log_add(log_guide, 2.0 * snapshot->q[row].log_abs);
    }
    if (!isfinite(log_guide)) return MVMC_KRYLOV_STATUS_NONFINITE;
    *guide_denominator = exp(log(eta) - log_guide);
    for (row = 0; status == MVMC_KRYLOV_STATUS_OK && row < 2; ++row) {
        for (column = row; status == MVMC_KRYLOV_STATUS_OK && column < 2;
             ++column, ++packed) {
            status = normalized_product(snapshot->q + row, snapshot->q + column,
                                        log_guide, overlap_upper + packed);
            if (status == MVMC_KRYLOV_STATUS_OK)
                status = normalized_product(snapshot->q + row,
                                            snapshot->r + column, log_guide,
                                            hamiltonian_forward_upper + packed);
            if (status == MVMC_KRYLOV_STATUS_OK)
                status = normalized_product(snapshot->q + column,
                                            snapshot->r + row, log_guide,
                                            hamiltonian_reverse_upper + packed);
            if (status == MVMC_KRYLOV_STATUS_OK)
                status = normalized_product(snapshot->r + row,
                                            snapshot->r + column, log_guide,
                                            hamiltonian_squared_upper + packed);
        }
    }
    return status;
}

static MVMCKrylovStatus linear_combination(const MVMCScaledComplex value[2],
                                           const double complex coefficient[2],
                                           MVMCScaledComplex *result) {
    MVMCScaledComplex terms[2];
    size_t index;
    for (index = 0; index < 2; ++index) {
        MVMCScaledComplex factor;
        if (!finite_complex(coefficient[index]))
            return MVMC_KRYLOV_STATUS_NONFINITE;
        if (creal(coefficient[index]) == 0.0 &&
            cimag(coefficient[index]) == 0.0) {
            if (mvmc_scaled_complex_make_exact_zero(terms + index) !=
                MVMC_PFAFFIAN_STATUS_OK)
                return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
            continue;
        }
        if (mvmc_scaled_complex_from_raw(coefficient[index], &factor) !=
                MVMC_PFAFFIAN_STATUS_OK ||
            mvmc_scaled_complex_multiply(&factor, value + index,
                                         terms + index) !=
                MVMC_PFAFFIAN_STATUS_OK)
            return MVMC_KRYLOV_STATUS_NONFINITE;
    }
    return mvmc_scaled_complex_sum_ordered(terms, 2, result) ==
                   MVMC_PFAFFIAN_STATUS_OK
               ? MVMC_KRYLOV_STATUS_OK
               : MVMC_KRYLOV_STATUS_NONFINITE;
}

MVMCKrylovStatus mvmc_power_lanczos_independent_final_evaluate(
    const MVMCPowerLanczosIndependentSnapshot *snapshot,
    const double complex coefficient[2], MVMCScaledComplex *phi,
    MVMCScaledComplex *h_phi, MVMCScaledComplex *local_energy,
    double *log_weight) {
    MVMCKrylovStatus status;
    double complex quotient_phase;
    double quotient_log;
    if (snapshot == NULL || coefficient == NULL || phi == NULL ||
        h_phi == NULL || local_energy == NULL || log_weight == NULL)
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    status = linear_combination(snapshot->q, coefficient, phi);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = linear_combination(snapshot->r, coefficient, h_phi);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    if (phi->state != MVMC_SCALED_COMPLEX_FINITE_NONZERO) {
        *log_weight =
            phi->state == MVMC_SCALED_COMPLEX_EXACT_ZERO ? -INFINITY : NAN;
        (void)mvmc_scaled_complex_make_exact_zero(local_energy);
        return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    }
    *log_weight = 2.0 * phi->log_abs;
    if (h_phi->state == MVMC_SCALED_COMPLEX_EXACT_ZERO)
        return mvmc_scaled_complex_make_exact_zero(local_energy) ==
                       MVMC_PFAFFIAN_STATUS_OK
                   ? MVMC_KRYLOV_STATUS_OK
                   : MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    if (h_phi->state != MVMC_SCALED_COMPLEX_FINITE_NONZERO)
        return MVMC_KRYLOV_STATUS_NONFINITE;
    quotient_phase = h_phi->phase * conj(phi->phase);
    quotient_log = h_phi->log_abs - phi->log_abs;
    if (!isfinite(*log_weight) || !isfinite(quotient_log) ||
        mvmc_scaled_complex_make_finite(quotient_phase, quotient_log,
                                        quotient_log + log(128.0 * DBL_EPSILON),
                                        local_energy) !=
            MVMC_PFAFFIAN_STATUS_OK)
        return MVMC_KRYLOV_STATUS_NONFINITE;
    return MVMC_KRYLOV_STATUS_OK;
}
