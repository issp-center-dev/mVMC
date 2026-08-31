#include "power_lanczos_independent.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE) || \
    !defined(MVMC_ENABLE_POWER_LANCZOS_P5_CORE)
#error "independent power-Lanczos driver requires the production core"
#endif

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "krylov_gevp_solver.h"
#include "krylov_positive_guide.h"
#include "krylov_positive_sampler.h"
#include "krylov_streaming_statistics.h"

enum { INDEPENDENT_PACKED = 3, INDEPENDENT_MAX_BLOCKS = 16 };

typedef struct {
    double complex sum;
    double complex compensation;
} ComplexSum;

typedef struct {
    uint64_t count;
    ComplexSum overlap[INDEPENDENT_PACKED];
    ComplexSum forward[INDEPENDENT_PACKED];
    ComplexSum reverse[INDEPENDENT_PACKED];
    ComplexSum squared[INDEPENDENT_PACKED];
} MatrixBlock;

typedef struct {
    MVMCPowerLanczosIndependentSnapshot amplitude;
    MVMCKrylovPositiveGuideResult guide;
} CoefficientSnapshot;

typedef struct {
    MVMCPowerLanczosIndependentSnapshot amplitude;
    MVMCScaledComplex phi;
    MVMCScaledComplex h_phi;
    MVMCScaledComplex local_energy;
    double log_weight;
} FinalSnapshot;

static int finite_complex(double complex value) {
    return isfinite(creal(value)) && isfinite(cimag(value));
}

static void invalidate_result(MVMCKrylovStatus status,
                              MVMCPowerLanczosIndependentResult *result) {
    if (result == NULL) return;
    memset(result, 0, sizeof(*result));
    result->version = MVMC_POWER_LANCZOS_INDEPENDENT_VERSION;
    result->status = status;
    result->energy = NAN;
    result->energy_standard_error = NAN;
    result->variance = NAN;
    result->variance_standard_error = NAN;
    result->final_energy_imaginary = NAN;
    result->variance_imaginary = NAN;
    result->energy_tau_int = NAN;
    result->effective_sample_count = NAN;
    result->condition_estimate = NAN;
    result->gevp_residual = NAN;
}

static MVMCKrylovStatus synchronize_status(
    MVMCKrylovStatus local, MVMCClassicPfaffianCommunicator communicator) {
#ifdef _mpi_use
    int encoded = (int)local;
    int reduced = 0;
    if (MPI_Allreduce(&encoded, &reduced, 1, MPI_INT, MPI_MAX, communicator) !=
            MPI_SUCCESS ||
        reduced < (int)MVMC_KRYLOV_STATUS_OK ||
        reduced > (int)MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE)
        return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    return (MVMCKrylovStatus)reduced;
#else
    (void)communicator;
    return local;
#endif
}

static uint64_t mix_u64(uint64_t value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

static uint64_t domain_seed(uint64_t seed, uint64_t domain) {
    const uint64_t result = mix_u64(seed ^ domain);
    return result == 0 ? UINT64_C(1) : result;
}

static size_t block_count_for(size_t sample_count) {
    size_t blocks;
    for (blocks = INDEPENDENT_MAX_BLOCKS; blocks >= 4; --blocks) {
        if (sample_count % blocks == 0 && sample_count / blocks >= 2)
            return blocks;
    }
    return 0;
}

static void complex_sum_add(ComplexSum *sum, double complex value) {
    const double complex updated = sum->sum + value;
    if (cabs(sum->sum) >= cabs(value))
        sum->compensation += (sum->sum - updated) + value;
    else
        sum->compensation += (value - updated) + sum->sum;
    sum->sum = updated;
}

static double complex complex_sum_value(const ComplexSum *sum) {
    return sum->sum + sum->compensation;
}

static void matrix_block_add_sample(MatrixBlock *block,
                                    const double complex overlap[3],
                                    const double complex forward[3],
                                    const double complex reverse[3],
                                    const double complex squared[3]) {
    size_t index;
    ++block->count;
    for (index = 0; index < INDEPENDENT_PACKED; ++index) {
        complex_sum_add(block->overlap + index, overlap[index]);
        complex_sum_add(block->forward + index, forward[index]);
        complex_sum_add(block->reverse + index, reverse[index]);
        complex_sum_add(block->squared + index, squared[index]);
    }
}

static void matrix_blocks_combine(const MatrixBlock *blocks, size_t block_count,
                                  size_t omitted, double complex overlap[3],
                                  double complex forward[3],
                                  double complex reverse[3],
                                  double complex squared[3], uint64_t *count) {
    size_t block, index;
    memset(overlap, 0, 3 * sizeof(*overlap));
    memset(forward, 0, 3 * sizeof(*forward));
    memset(reverse, 0, 3 * sizeof(*reverse));
    memset(squared, 0, 3 * sizeof(*squared));
    *count = 0;
    for (block = 0; block < block_count; ++block) {
        if (block == omitted) continue;
        *count += blocks[block].count;
        for (index = 0; index < INDEPENDENT_PACKED; ++index) {
            overlap[index] += complex_sum_value(blocks[block].overlap + index);
            forward[index] += complex_sum_value(blocks[block].forward + index);
            reverse[index] += complex_sum_value(blocks[block].reverse + index);
            squared[index] += complex_sum_value(blocks[block].squared + index);
        }
    }
}

static MVMCKrylovStatus map_gevp(MVMCKrylovGEVPStatus status) {
    if (status == MVMC_KRYLOV_GEVP_OK) return MVMC_KRYLOV_STATUS_OK;
    if (status == MVMC_KRYLOV_GEVP_INVALID_ARGUMENT)
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    if (status == MVMC_KRYLOV_GEVP_NONFINITE_INPUT)
        return MVMC_KRYLOV_STATUS_NONFINITE;
    return MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
}

static MVMCKrylovStatus solve_matrix(const MVMCKrylovGEVPPolicy *policy,
                                     const MatrixBlock *blocks,
                                     size_t block_count, size_t omitted,
                                     MVMCKrylovGEVPResult *result,
                                     double complex *second_moment) {
    double complex overlap[3], forward[3], reverse[3], squared[3];
    uint64_t count = 0;
    double complex normalization = 0.0;
    double complex second = 0.0;
    MVMCKrylovGEVPStatus gevp_status;
    int row, column;
    size_t packed = 0;
    matrix_blocks_combine(blocks, block_count, omitted, overlap, forward,
                          reverse, squared, &count);
    if (count == 0) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    gevp_status = mvmc_krylov_gevp_solve_complex_packed(
        policy, 2, overlap, forward, reverse, squared, 3, result);
    if (gevp_status != MVMC_KRYLOV_GEVP_OK) return map_gevp(gevp_status);
    for (row = 0; row < 2; ++row) {
        for (column = row; column < 2; ++column, ++packed) {
            normalization += conj(result->coefficient[row]) * overlap[packed] *
                             result->coefficient[column];
            second += conj(result->coefficient[row]) * squared[packed] *
                      result->coefficient[column];
            if (row != column) {
                normalization += conj(result->coefficient[column]) *
                                 conj(overlap[packed]) *
                                 result->coefficient[row];
                second += conj(result->coefficient[column]) *
                          conj(squared[packed]) * result->coefficient[row];
            }
        }
    }
    if (!finite_complex(normalization) || !finite_complex(second) ||
        !(creal(normalization) > 0.0))
        return MVMC_KRYLOV_STATUS_NONFINITE;
    *second_moment = second / normalization;
    return finite_complex(*second_moment) ? MVMC_KRYLOV_STATUS_OK
                                          : MVMC_KRYLOV_STATUS_NONFINITE;
}

static MVMCKrylovStatus draw_bounded(void *context, size_t bound,
                                     size_t *value) {
    return mvmc_krylov_positive_sampler_rng_draw_bounded(
        (MVMCKrylovPositiveSamplerRng *)context, bound, value);
}

static MVMCKrylovStatus draw_uniform(
    const MVMCKrylovFockModel *model,
    MVMCClassicPfaffianCommunicator communicator,
    MVMCKrylovPositiveSamplerRng *rng, uint64_t *words, size_t word_count) {
    MVMCKrylovFockUniformProposalResult drawn;
    MVMCKrylovStatus status = mvmc_krylov_fock_proposal_draw_uniform_sector(
        model, draw_bounded, rng, words, word_count, &drawn);
    if (status == MVMC_KRYLOV_STATUS_OK && !drawn.valid)
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
#ifdef _mpi_use
    if (status == MVMC_KRYLOV_STATUS_OK &&
        (word_count > (size_t)INT_MAX ||
         MPI_Bcast(words, (int)word_count, MPI_UINT64_T, 0, communicator) !=
             MPI_SUCCESS))
        status = MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
#endif
    return synchronize_status(status, communicator);
}

static MVMCKrylovStatus coefficient_evaluate(
    const MVMCPowerLanczosIndependentStreamingModel *physical,
    MVMCKrylovBoundedWorkspace *workspace, const uint64_t *words,
    size_t word_count, MVMCKrylovScaledAmplitudeCallback amplitude,
    void *amplitude_context, const MVMCKrylovPositiveGuidePolicy *guide,
    CoefficientSnapshot *snapshot) {
    MVMCKrylovStatus status = mvmc_power_lanczos_independent_evaluate_view(
        physical, workspace, words, word_count, amplitude, amplitude_context,
        &snapshot->amplitude);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_krylov_positive_guide_evaluate(
            guide, snapshot->amplitude.q, 2, &snapshot->guide);
    return status;
}

static MVMCKrylovStatus coefficient_step(
    const MVMCPowerLanczosIndependentStreamingModel *physical,
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const MVMCKrylovPositiveGuidePolicy *guide, uint64_t *words,
    size_t word_count, CoefficientSnapshot *current,
    MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words) {
    MVMCKrylovPositiveSamplerProposalDrawResult draw;
    MVMCKrylovPositiveGuideAcceptance acceptance;
    CoefficientSnapshot proposal;
    MVMCKrylovStatus status = mvmc_krylov_positive_sampler_draw_mixture_rng(
        proposal_model, proposal_policy, words, word_count, rng, proposal_words,
        word_count, &draw);
    if (status == MVMC_KRYLOV_STATUS_OK) {
        *rng = draw.rng_after;
        status = coefficient_evaluate(physical, workspace, proposal_words,
                                      word_count, amplitude, amplitude_context,
                                      guide, &proposal);
    }
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_krylov_positive_guide_acceptance(
            &current->guide, &proposal.guide, draw.log_proposal_ratio,
            draw.uniform, &acceptance);
    if (status == MVMC_KRYLOV_STATUS_OK && acceptance.accepted) {
        memcpy(words, proposal_words, word_count * sizeof(*words));
        *current = proposal;
    }
    return status;
}

static MVMCKrylovStatus final_evaluate(
    const MVMCPowerLanczosIndependentStreamingModel *physical,
    MVMCKrylovBoundedWorkspace *workspace, const uint64_t *words,
    size_t word_count, MVMCKrylovScaledAmplitudeCallback amplitude,
    void *amplitude_context, const double complex coefficient[2],
    FinalSnapshot *snapshot, int *is_node) {
    MVMCKrylovStatus status;
    if (is_node == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    *is_node = 0;
    status = mvmc_power_lanczos_independent_evaluate_view(
        physical, workspace, words, word_count, amplitude, amplitude_context,
        &snapshot->amplitude);
    if (status != MVMC_KRYLOV_STATUS_OK) return status;
    status = mvmc_power_lanczos_independent_final_evaluate(
        &snapshot->amplitude, coefficient, &snapshot->phi, &snapshot->h_phi,
        &snapshot->local_energy, &snapshot->log_weight);
    if (status == MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE) {
        *is_node = 1;
        status = MVMC_KRYLOV_STATUS_OK;
    }
    return status;
}

static MVMCKrylovStatus final_step(
    const MVMCPowerLanczosIndependentStreamingModel *physical,
    MVMCKrylovBoundedWorkspace *workspace,
    const MVMCKrylovFockModel *proposal_model,
    const MVMCKrylovPositiveSamplerProposalPolicy *proposal_policy,
    const double complex coefficient[2], uint64_t *words, size_t word_count,
    FinalSnapshot *current, MVMCKrylovPositiveSamplerRng *rng,
    MVMCKrylovScaledAmplitudeCallback amplitude, void *amplitude_context,
    uint64_t *proposal_words) {
    MVMCKrylovPositiveSamplerProposalDrawResult draw;
    FinalSnapshot proposal;
    int proposal_is_node = 0;
    MVMCKrylovStatus status = mvmc_krylov_positive_sampler_draw_mixture_rng(
        proposal_model, proposal_policy, words, word_count, rng, proposal_words,
        word_count, &draw);
    if (status == MVMC_KRYLOV_STATUS_OK) {
        *rng = draw.rng_after;
        status = final_evaluate(physical, workspace, proposal_words, word_count,
                                amplitude, amplitude_context, coefficient,
                                &proposal, &proposal_is_node);
    }
    if (status == MVMC_KRYLOV_STATUS_OK && proposal_is_node)
        return MVMC_KRYLOV_STATUS_OK;
    if (status == MVMC_KRYLOV_STATUS_OK) {
        const double log_ratio =
            proposal.log_weight - current->log_weight + draw.log_proposal_ratio;
        if (!isfinite(log_ratio)) return MVMC_KRYLOV_STATUS_NONFINITE;
        if (log_ratio >= 0.0 || draw.uniform == 0.0 ||
            log(draw.uniform) < log_ratio) {
            memcpy(words, proposal_words, word_count * sizeof(*words));
            *current = proposal;
        }
    }
    return status;
}

static double standard_error(const double complex *block_sum,
                             const uint64_t *block_count, size_t blocks) {
    double mean = 0.0, variance = 0.0;
    size_t block;
    for (block = 0; block < blocks; ++block)
        mean += creal(block_sum[block]) / (double)block_count[block];
    mean /= (double)blocks;
    for (block = 0; block < blocks; ++block) {
        const double value =
            creal(block_sum[block]) / (double)block_count[block];
        variance += (value - mean) * (value - mean);
    }
    return sqrt(variance / ((double)blocks * ((double)blocks - 1.0)));
}

static double jackknife_se(const double *leave_one, size_t blocks) {
    double mean = 0.0, variance = 0.0;
    size_t block;
    for (block = 0; block < blocks; ++block) mean += leave_one[block];
    mean /= (double)blocks;
    for (block = 0; block < blocks; ++block)
        variance += (leave_one[block] - mean) * (leave_one[block] - mean);
    return sqrt(((double)blocks - 1.0) / (double)blocks * variance);
}

MVMCKrylovStatus mvmc_power_lanczos_independent_run(
    const MVMCPowerLanczosIndependentInput *input,
    MVMCPowerLanczosIndependentResult *result) {
    MVMCPowerLanczosClassicBridge *bridge = NULL;
    MVMCPowerLanczosIndependentStreamingModel *physical = NULL;
    MVMCPowerLanczosIndependentModel *hprime_owned = NULL;
    const MVMCKrylovFockModel *hprime = NULL;
    const MVMCKrylovFockModel *proposal_model = NULL;
    MVMCKrylovScaledAmplitudeCallback amplitude = NULL;
    void *amplitude_context = NULL;
    MVMCKrylovBoundedPlan *plan = NULL;
    MVMCKrylovBoundedWorkspace *workspace = NULL;
    MVMCKrylovBoundedLimits limits;
    MVMCKrylovPositiveSamplerProposalPolicy proposal_policy;
    MVMCKrylovPositiveGuidePolicy guide;
    MVMCKrylovPositiveSamplerRng bootstrap_rng, coefficient_rng, final_rng;
    CoefficientSnapshot coefficient_snapshot;
    FinalSnapshot final_snapshot;
    MatrixBlock blocks[INDEPENDENT_MAX_BLOCKS];
    MVMCKrylovGEVPPolicy gevp_policy;
    MVMCKrylovGEVPResult gevp;
    MVMCKrylovTauIntResult tau;
    double complex final_block_sum[INDEPENDENT_MAX_BLOCKS];
    uint64_t final_block_count[INDEPENDENT_MAX_BLOCKS];
    double leave_second[INDEPENDENT_MAX_BLOCKS];
    double *energy_series = NULL;
    uint64_t *words = NULL, *proposal_words = NULL;
    size_t word_count = 0, block_count = 0;
    double complex energy_sum = 0.0, energy = NAN + I * NAN;
    double complex second = NAN + I * NAN, variance = NAN + I * NAN;
    double energy_se = NAN, second_se = NAN;
    uint64_t generation_hash = 0, physical_row_terms = 0;
    MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
    size_t step, sample, block;

    if (result == NULL) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    invalidate_result(MVMC_KRYLOV_STATUS_INVALID_ARGUMENT, result);
    if (input == NULL || input->classic_view == NULL || input->seed == 0 ||
        input->sample_count < 8 || input->interval == 0 ||
        (input->hprime_transfer_count == 0 &&
         input->hprime_inter_all_count == 0) ||
        (input->physical_transfer_count == 0 &&
         input->physical_inter_all_count == 0) ||
        input->sample_count > SIZE_MAX / sizeof(*energy_series))
        return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    block_count = block_count_for(input->sample_count);
    if (block_count == 0) return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    status = mvmc_power_lanczos_classic_bridge_create(
        input->classic_view, input->communicator, input->communicator, &bridge);
    if (status == MVMC_KRYLOV_STATUS_OK) {
        amplitude = mvmc_power_lanczos_classic_bridge_amplitude(bridge);
        amplitude_context =
            mvmc_power_lanczos_classic_bridge_amplitude_context(bridge);
        proposal_model = mvmc_power_lanczos_classic_bridge_model(bridge);
        generation_hash =
            mvmc_power_lanczos_classic_bridge_generation_hash(bridge);
        if (amplitude == NULL || amplitude_context == NULL ||
            proposal_model == NULL || generation_hash == 0)
            status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_power_lanczos_independent_streaming_model_create(
            input->classic_view->site_count,
            input->classic_view->up_electron_count,
            input->classic_view->down_electron_count,
            input->physical_transfer_count, input->physical_transfer_indices,
            input->physical_transfer_parameters,
            input->physical_inter_all_count, input->physical_inter_all_indices,
            input->physical_inter_all_parameters, &physical);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_power_lanczos_independent_model_create(
            input->classic_view->site_count,
            input->classic_view->up_electron_count,
            input->classic_view->down_electron_count,
            input->hprime_transfer_count, input->hprime_transfer_indices,
            input->hprime_transfer_parameters, input->hprime_inter_all_count,
            input->hprime_inter_all_indices, input->hprime_inter_all_parameters,
            &hprime_owned);
    hprime = mvmc_power_lanczos_independent_model_view(hprime_owned);
    if (status == MVMC_KRYLOV_STATUS_OK) {
        word_count = mvmc_krylov_fock_word_count(
            (size_t)input->classic_view->site_count);
        words = (uint64_t *)calloc(word_count, sizeof(*words));
        proposal_words =
            (uint64_t *)calloc(word_count, sizeof(*proposal_words));
        energy_series =
            (double *)calloc(input->sample_count, sizeof(*energy_series));
        if (words == NULL || proposal_words == NULL || energy_series == NULL)
            status = MVMC_KRYLOV_STATUS_ALLOCATION_FAILURE;
    }
    memset(&limits, 0, sizeof(limits));
    if (status == MVMC_KRYLOV_STATUS_OK) {
        const uint64_t terms = (uint64_t)hprime->term_count;
        if (terms > UINT64_MAX - UINT64_C(2))
            status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
        else {
            limits.amplitude_policy_hash = generation_hash;
            limits.cache_bytes = (size_t)16 * 1024 * 1024;
            limits.max_row_transitions = hprime->term_count;
            limits.max_workspace_bytes = (size_t)512 * 1024 * 1024;
            limits.max_node_expansions = terms + UINT64_C(2);
            limits.max_terminal_amplitude_calls = terms + UINT64_C(1);
            limits.max_total_row_transitions = terms;
            limits.max_order = 1;
        }
    }
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_bounded_krylov_plan_create(hprime, &limits, &plan);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_bounded_krylov_workspace_create(plan, &workspace);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_krylov_positive_sampler_proposal_policy_create(
            1, 16, &proposal_policy);
    memset(&guide, 0, sizeof(guide));
    guide.order = 1;
    guide.eta = 0x1p-40;
    guide.lambda[0] = guide.lambda[1] = 1.0;
    guide.policy_hash = UINT64_C(0x494e4447504f4c31);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_krylov_positive_sampler_rng_seed(
            domain_seed(input->seed, UINT64_C(0x494e44424f4f5431)),
            UINT64_C(0x494e44424f4f5431), &bootstrap_rng);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = draw_uniform(proposal_model, input->communicator,
                              &bootstrap_rng, words, word_count);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_bounded_krylov_session_begin(
            workspace, amplitude, amplitude_context, generation_hash);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = coefficient_evaluate(physical, workspace, words, word_count,
                                      amplitude, amplitude_context, &guide,
                                      &coefficient_snapshot);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_krylov_positive_sampler_rng_seed(
            domain_seed(input->seed, UINT64_C(0x494e44434f454631)),
            UINT64_C(0x494e44434f454631), &coefficient_rng);
    for (step = 0; status == MVMC_KRYLOV_STATUS_OK && step < input->warm_up;
         ++step)
        status = coefficient_step(physical, workspace, proposal_model,
                                  &proposal_policy, &guide, words, word_count,
                                  &coefficient_snapshot, &coefficient_rng,
                                  amplitude, amplitude_context, proposal_words);
    memset(blocks, 0, sizeof(blocks));
    for (sample = 0;
         status == MVMC_KRYLOV_STATUS_OK && sample < input->sample_count;
         ++sample) {
        double complex overlap[3], forward[3], reverse[3], squared[3];
        double denominator;
        for (step = 0;
             status == MVMC_KRYLOV_STATUS_OK && step < input->interval; ++step)
            status = coefficient_step(
                physical, workspace, proposal_model, &proposal_policy, &guide,
                words, word_count, &coefficient_snapshot, &coefficient_rng,
                amplitude, amplitude_context, proposal_words);
        if (status == MVMC_KRYLOV_STATUS_OK)
            status = mvmc_power_lanczos_independent_matrix_sample(
                &coefficient_snapshot.amplitude, guide.eta, overlap, forward,
                reverse, squared, &denominator);
        (void)denominator;
        if (status == MVMC_KRYLOV_STATUS_OK) {
            block = sample * block_count / input->sample_count;
            matrix_block_add_sample(blocks + block, overlap, forward, reverse,
                                    squared);
            if (UINT64_MAX - physical_row_terms <
                coefficient_snapshot.amplitude.physical_row_terms)
                status = MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
            else
                physical_row_terms +=
                    coefficient_snapshot.amplitude.physical_row_terms;
        }
    }
    status = synchronize_status(status, input->communicator);
    if (status == MVMC_KRYLOV_STATUS_OK &&
        mvmc_krylov_gevp_default_policy(0x1p-40, &gevp_policy) !=
            MVMC_KRYLOV_GEVP_OK)
        status = MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = solve_matrix(&gevp_policy, blocks, block_count, SIZE_MAX,
                              &gevp, &second);
    if (status == MVMC_KRYLOV_STATUS_OK)
        status = mvmc_krylov_positive_sampler_rng_seed(
            domain_seed(input->seed, UINT64_C(0x494e4446494e4131)),
            UINT64_C(0x494e4446494e4131), &final_rng);
    if (status == MVMC_KRYLOV_STATUS_OK) {
        int final_origin_found = 0;
        for (step = 0; step < 1024 && !final_origin_found; ++step) {
            int final_origin_is_node = 0;
            status = draw_uniform(proposal_model, input->communicator,
                                  &final_rng, words, word_count);
            if (status != MVMC_KRYLOV_STATUS_OK) break;
            status =
                final_evaluate(physical, workspace, words, word_count,
                               amplitude, amplitude_context, gevp.coefficient,
                               &final_snapshot, &final_origin_is_node);
            if (status == MVMC_KRYLOV_STATUS_OK && !final_origin_is_node)
                final_origin_found = 1;
        }
        if (status == MVMC_KRYLOV_STATUS_OK && !final_origin_found)
            status = MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE;
    }
    for (step = 0; status == MVMC_KRYLOV_STATUS_OK && step < input->warm_up;
         ++step)
        status = final_step(physical, workspace, proposal_model,
                            &proposal_policy, gevp.coefficient, words,
                            word_count, &final_snapshot, &final_rng, amplitude,
                            amplitude_context, proposal_words);
    memset(final_block_sum, 0, sizeof(final_block_sum));
    memset(final_block_count, 0, sizeof(final_block_count));
    for (sample = 0;
         status == MVMC_KRYLOV_STATUS_OK && sample < input->sample_count;
         ++sample) {
        double complex local_energy = NAN + I * NAN;
        for (step = 0;
             status == MVMC_KRYLOV_STATUS_OK && step < input->interval; ++step)
            status = final_step(physical, workspace, proposal_model,
                                &proposal_policy, gevp.coefficient, words,
                                word_count, &final_snapshot, &final_rng,
                                amplitude, amplitude_context, proposal_words);
        if (status == MVMC_KRYLOV_STATUS_OK) {
            const MVMCScaledComplexExportStatus export_status =
                mvmc_scaled_complex_export_common_scale(
                    &final_snapshot.local_energy, 0.0, &local_energy);
            if (export_status == MVMC_SCALED_EXPORT_EXACT_ZERO)
                local_energy = 0.0;
            else if (export_status != MVMC_SCALED_EXPORT_OK ||
                     !finite_complex(local_energy))
                status = MVMC_KRYLOV_STATUS_NONFINITE;
        }
        if (status == MVMC_KRYLOV_STATUS_OK) {
            block = sample * block_count / input->sample_count;
            energy_series[sample] = creal(local_energy);
            energy_sum += local_energy;
            final_block_sum[block] += local_energy;
            ++final_block_count[block];
        }
    }
    if (workspace != NULL && mvmc_bounded_krylov_session_is_active(workspace)) {
        const MVMCKrylovStatus end = mvmc_bounded_krylov_session_end(workspace);
        if (status == MVMC_KRYLOV_STATUS_OK) status = end;
    }
    status = synchronize_status(status, input->communicator);
    if (status == MVMC_KRYLOV_STATUS_OK) {
        energy = energy_sum / (double)input->sample_count;
        variance = second - energy * energy;
        energy_se =
            standard_error(final_block_sum, final_block_count, block_count);
    }
    for (block = 0; status == MVMC_KRYLOV_STATUS_OK && block < block_count;
         ++block) {
        MVMCKrylovGEVPResult leave_gevp;
        double complex leave = NAN + I * NAN;
        status = solve_matrix(&gevp_policy, blocks, block_count, block,
                              &leave_gevp, &leave);
        leave_second[block] = creal(leave);
    }
    if (status == MVMC_KRYLOV_STATUS_OK) {
        second_se = jackknife_se(leave_second, block_count);
        status = mvmc_krylov_tau_int_geyer_initial_positive(
            energy_series, input->sample_count, &tau);
    }
    if (status == MVMC_KRYLOV_STATUS_OK &&
        (!finite_complex(energy) || !finite_complex(variance) ||
         !isfinite(energy_se) || !isfinite(second_se)))
        status = MVMC_KRYLOV_STATUS_NONFINITE;
    if (status == MVMC_KRYLOV_STATUS_OK) {
        result->valid = 1;
        result->status = MVMC_KRYLOV_STATUS_OK;
        result->version = MVMC_POWER_LANCZOS_INDEPENDENT_VERSION;
        result->coefficient_samples = (uint64_t)input->sample_count;
        result->final_samples = (uint64_t)input->sample_count;
        result->physical_row_terms = physical_row_terms;
        result->block_count = block_count;
        result->retained_rank = gevp.retained_rank;
        result->condition_estimate = gevp.condition_estimate;
        result->gevp_residual = gevp.gevp_relative_residual;
        result->energy = creal(energy);
        result->energy_standard_error = energy_se;
        result->variance = creal(variance);
        result->variance_standard_error =
            hypot(second_se, 2.0 * cabs(energy) * energy_se);
        result->final_energy_imaginary = cimag(energy);
        result->variance_imaginary = cimag(variance);
        result->energy_tau_int = tau.tau_int;
        result->effective_sample_count =
            (double)input->sample_count / (2.0 * tau.tau_int);
    } else {
        invalidate_result(status, result);
    }
    free(energy_series);
    free(proposal_words);
    free(words);
    mvmc_bounded_krylov_workspace_destroy(workspace);
    mvmc_bounded_krylov_plan_destroy(plan);
    mvmc_power_lanczos_independent_model_destroy(hprime_owned);
    mvmc_power_lanczos_independent_streaming_model_destroy(physical);
    mvmc_power_lanczos_classic_bridge_destroy(bridge);
    return status;
}
