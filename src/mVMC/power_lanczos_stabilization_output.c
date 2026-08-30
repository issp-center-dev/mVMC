#include "power_lanczos_stabilization_output.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "stabilization output requires the bounded power-Lanczos engine"
#endif

#include "power_lanczos_json_writer.h"

#include <complex.h>
#include <math.h>
#include <string.h>

#define RETURN_IF_ERROR(expression)                \
  do {                                             \
    const MVMCKrylovStatus status_ = (expression); \
    if (status_ != MVMC_KRYLOV_STATUS_OK) {        \
      return status_;                              \
    }                                              \
  } while (0)

static int lowercase_hex(const char *value, size_t length) {
  size_t index;
  if (value == NULL || strlen(value) != length) return 0;
  for (index = 0; index < length; ++index) {
    const unsigned char c = (unsigned char)value[index];
    if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) return 0;
  }
  return 1;
}

static int identity_valid(
    const MVMCPowerLanczosStabilizationOutputIdentity *identity) {
  size_t source_length;
  if (identity == NULL || identity->run_id == NULL ||
      identity->environment_id == NULL || identity->seed_id == NULL ||
      identity->source_commit == NULL || identity->input_sha256 == NULL ||
      identity->binary_sha256 == NULL || identity->run_id[0] == '\0' ||
      identity->environment_id[0] == '\0' || identity->seed_id[0] == '\0' ||
      !mvmc_power_lanczos_json_public_string_valid(identity->run_id) ||
      !mvmc_power_lanczos_json_public_string_valid(
          identity->environment_id) ||
      !mvmc_power_lanczos_json_public_string_valid(identity->seed_id)) {
    return 0;
  }
  source_length = strlen(identity->source_commit);
  return (source_length == 40 || source_length == 64) &&
         lowercase_hex(identity->source_commit, source_length) &&
         lowercase_hex(identity->input_sha256, 64) &&
         lowercase_hex(identity->binary_sha256, 64);
}

static void digest_hex(const unsigned char digest[32], char output[65]) {
  static const char digits[] = "0123456789abcdef";
  size_t index;
  for (index = 0; index < 32; ++index) {
    output[2 * index] = digits[digest[index] >> 4];
    output[2 * index + 1] = digits[digest[index] & 0x0f];
  }
  output[64] = '\0';
}

static MVMCKrylovStatus key_double(
    MVMCPowerLanczosJsonWriter *writer, const char *key, double value) {
  RETURN_IF_ERROR(mvmc_power_lanczos_json_key(writer, key));
  return mvmc_power_lanczos_json_double(writer, value);
}

static MVMCKrylovStatus key_int64(
    MVMCPowerLanczosJsonWriter *writer, const char *key, int64_t value) {
  RETURN_IF_ERROR(mvmc_power_lanczos_json_key(writer, key));
  return mvmc_power_lanczos_json_int64(writer, value);
}

static MVMCKrylovStatus key_uint64(
    MVMCPowerLanczosJsonWriter *writer, const char *key, uint64_t value) {
  RETURN_IF_ERROR(mvmc_power_lanczos_json_key(writer, key));
  return mvmc_power_lanczos_json_uint64(writer, value);
}

static MVMCKrylovStatus key_string(
    MVMCPowerLanczosJsonWriter *writer, const char *key,
    const char *value) {
  RETURN_IF_ERROR(mvmc_power_lanczos_json_key(writer, key));
  return mvmc_power_lanczos_json_string(writer, value);
}

static MVMCKrylovStatus write_complex(
    MVMCPowerLanczosJsonWriter *writer, double complex value) {
  if (!isfinite(creal(value)) || !isfinite(cimag(value))) {
    return MVMC_KRYLOV_STATUS_NONFINITE;
  }
  RETURN_IF_ERROR(mvmc_power_lanczos_json_array_begin(writer));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_double(writer, creal(value)));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_double(writer, cimag(value)));
  return mvmc_power_lanczos_json_array_end(writer);
}

static MVMCKrylovStatus write_blocks(
    MVMCPowerLanczosJsonWriter *writer,
    const MVMCPowerLanczosProductionSession *session, size_t block_count) {
  size_t block;
  RETURN_IF_ERROR(mvmc_power_lanczos_json_array_begin(writer));
  for (block = 0; block < block_count; ++block) {
    MVMCPowerLanczosProductionCoefficientBlockResult coefficient;
    MVMCPowerLanczosProductionFinalBlockResult final;
    RETURN_IF_ERROR(
        mvmc_power_lanczos_production_session_coefficient_block_result(
            session, block, &coefficient));
    RETURN_IF_ERROR(mvmc_power_lanczos_production_session_final_block_result(
        session, block, &final));
    if (!coefficient.valid || !final.valid) {
      return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
    }
    RETURN_IF_ERROR(mvmc_power_lanczos_json_object_begin(writer));
    RETURN_IF_ERROR(mvmc_power_lanczos_json_key(writer, "coefficient"));
    RETURN_IF_ERROR(mvmc_power_lanczos_json_object_begin(writer));
    RETURN_IF_ERROR(key_double(writer, "leave_one_projective_distance",
                               coefficient.leave_one_projective_distance));
    RETURN_IF_ERROR(key_uint64(writer, "leave_one_sample_count",
                               coefficient.leave_one_sample_count));
    RETURN_IF_ERROR(
        mvmc_power_lanczos_json_key(writer, "leave_one_second_moment"));
    RETURN_IF_ERROR(write_complex(writer, coefficient.leave_one_second_moment));
    RETURN_IF_ERROR(
        key_uint64(writer, "sample_count", coefficient.sample_count));
    RETURN_IF_ERROR(mvmc_power_lanczos_json_object_end(writer));
    RETURN_IF_ERROR(mvmc_power_lanczos_json_key(writer, "final"));
    RETURN_IF_ERROR(mvmc_power_lanczos_json_object_begin(writer));
    RETURN_IF_ERROR(mvmc_power_lanczos_json_key(writer, "energy_sum"));
    RETURN_IF_ERROR(write_complex(writer, final.energy_sum));
    RETURN_IF_ERROR(mvmc_power_lanczos_json_key(
        writer, "local_energy_abs_squared_sum"));
    RETURN_IF_ERROR(
        write_complex(writer, final.local_energy_abs_squared_sum));
    RETURN_IF_ERROR(key_uint64(writer, "sample_count", final.sample_count));
    RETURN_IF_ERROR(mvmc_power_lanczos_json_object_end(writer));
    RETURN_IF_ERROR(key_uint64(writer, "index", (uint64_t)block));
    RETURN_IF_ERROR(mvmc_power_lanczos_json_object_end(writer));
  }
  return mvmc_power_lanczos_json_array_end(writer);
}

static MVMCKrylovStatus write_gevp(
    MVMCPowerLanczosJsonWriter *writer,
    const MVMCPowerLanczosGEVPResult *gevp) {
  int index;
  if (gevp == NULL || !gevp->valid || gevp->dimension < 1 ||
      gevp->dimension > MVMC_KRYLOV_GEVP_MAX_DIMENSION) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  RETURN_IF_ERROR(mvmc_power_lanczos_json_object_begin(writer));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_key(writer, "coefficient"));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_array_begin(writer));
  for (index = 0; index < gevp->dimension; ++index) {
    RETURN_IF_ERROR(write_complex(writer, gevp->coefficient[index]));
  }
  RETURN_IF_ERROR(mvmc_power_lanczos_json_array_end(writer));
  RETURN_IF_ERROR(
      key_double(writer, "condition_estimate", gevp->condition_estimate));
  RETURN_IF_ERROR(
      key_int64(writer, "discarded_rank", gevp->discarded_rank));
  RETURN_IF_ERROR(key_double(writer, "energy", gevp->energy));
  RETURN_IF_ERROR(key_double(writer, "normalization", gevp->normalization));
  RETURN_IF_ERROR(key_double(writer, "normwise_backward_error",
                             gevp->normwise_backward_error));
  RETURN_IF_ERROR(
      mvmc_power_lanczos_json_key(writer, "overlap_eigenvalue"));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_array_begin(writer));
  for (index = 0; index < gevp->dimension; ++index) {
    RETURN_IF_ERROR(mvmc_power_lanczos_json_double(
        writer, gevp->overlap_eigenvalue[index]));
  }
  RETURN_IF_ERROR(mvmc_power_lanczos_json_array_end(writer));
  RETURN_IF_ERROR(key_int64(writer, "phase_pivot", gevp->phase_pivot));
  RETURN_IF_ERROR(key_double(writer, "rank_relative_cutoff",
                             gevp->rank_relative_cutoff));
  RETURN_IF_ERROR(
      key_double(writer, "relative_residual", gevp->relative_residual));
  RETURN_IF_ERROR(key_int64(writer, "retained_rank", gevp->retained_rank));
  RETURN_IF_ERROR(key_double(writer, "root_gap", gevp->root_gap));
  RETURN_IF_ERROR(
      key_int64(writer, "root_multiplicity", gevp->root_multiplicity));
  RETURN_IF_ERROR(key_double(writer, "variance", gevp->variance));
  return mvmc_power_lanczos_json_object_end(writer);
}

static MVMCKrylovStatus write_identity(
    MVMCPowerLanczosJsonWriter *writer,
    const MVMCPowerLanczosStabilizationOutputIdentity *identity) {
  RETURN_IF_ERROR(mvmc_power_lanczos_json_object_begin(writer));
  RETURN_IF_ERROR(
      key_string(writer, "binary_sha256", identity->binary_sha256));
  RETURN_IF_ERROR(
      key_string(writer, "environment_id", identity->environment_id));
  RETURN_IF_ERROR(
      key_string(writer, "input_sha256", identity->input_sha256));
  RETURN_IF_ERROR(key_string(writer, "run_id", identity->run_id));
  RETURN_IF_ERROR(key_string(writer, "seed_id", identity->seed_id));
  RETURN_IF_ERROR(
      key_string(writer, "source_commit", identity->source_commit));
  return mvmc_power_lanczos_json_object_end(writer);
}

static MVMCKrylovStatus write_policy(
    MVMCPowerLanczosJsonWriter *writer,
    const MVMCPowerLanczosProductionSessionSummary *summary) {
  size_t index;
  RETURN_IF_ERROR(mvmc_power_lanczos_json_object_begin(writer));
  RETURN_IF_ERROR(key_string(writer, "coefficient_guide_id",
                             summary->power_step == 1 ? "L1_FLAT"
                                                      : "L2_FLAT"));
  RETURN_IF_ERROR(key_double(writer, "eta_relative", summary->eta_relative));
  RETURN_IF_ERROR(key_string(writer, "gevp_policy_id",
                             mvmc_power_lanczos_gevp_policy_id()));
  RETURN_IF_ERROR(
      key_double(writer, "global_proposal_mixture", 1.0 / 16.0));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_key(writer, "log_basis_scale"));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_array_begin(writer));
  for (index = 0; index <= summary->basis_count; ++index) {
    RETURN_IF_ERROR(mvmc_power_lanczos_json_double(
        writer, summary->log_basis_scale[index]));
  }
  RETURN_IF_ERROR(mvmc_power_lanczos_json_array_end(writer));
  RETURN_IF_ERROR(key_double(writer, "near_zero_omission", 0.0));
  RETURN_IF_ERROR(key_string(writer, "rank_cutoff_id", "S40"));
  RETURN_IF_ERROR(key_double(writer, "resolved_eta", summary->resolved_eta));
  RETURN_IF_ERROR(key_uint64(writer, "scale_pilot_sample_count",
                             summary->scale_pilot_sample_count));
  RETURN_IF_ERROR(key_uint64(writer, "scale_pilot_warm_up",
                             summary->scale_pilot_warm_up));
  return mvmc_power_lanczos_json_object_end(writer);
}

static MVMCKrylovStatus write_statistics(
    MVMCPowerLanczosJsonWriter *writer,
    const MVMCPowerLanczosStabilizationResult *statistics) {
  RETURN_IF_ERROR(mvmc_power_lanczos_json_object_begin(writer));
  RETURN_IF_ERROR(key_uint64(writer, "coefficient_block_length",
                             (uint64_t)statistics->coefficient_block_length));
  RETURN_IF_ERROR(key_double(writer, "coefficient_effective_block_count",
                             statistics->coefficient_effective_block_count));
  RETURN_IF_ERROR(key_double(writer, "coefficient_second_moment",
                             statistics->coefficient_second_moment));
  RETURN_IF_ERROR(key_double(
      writer, "coefficient_second_moment_imaginary",
      statistics->coefficient_second_moment_imaginary));
  RETURN_IF_ERROR(key_double(
      writer, "coefficient_second_moment_imaginary_standard_error",
      statistics->coefficient_second_moment_imaginary_standard_error));
  RETURN_IF_ERROR(key_double(
      writer, "coefficient_second_moment_standard_error",
      statistics->coefficient_second_moment_standard_error));
  RETURN_IF_ERROR(key_double(writer, "coefficient_tau_int_max",
                             statistics->coefficient_tau_int_max));
  RETURN_IF_ERROR(key_double(writer, "coefficient_tau_to_block_length",
                             statistics->coefficient_tau_to_block_length));
  RETURN_IF_ERROR(key_double(writer, "energy", statistics->energy));
  RETURN_IF_ERROR(
      key_double(writer, "energy_imaginary", statistics->energy_imaginary));
  RETURN_IF_ERROR(key_double(
      writer, "energy_imaginary_standard_error",
      statistics->energy_imaginary_standard_error));
  RETURN_IF_ERROR(key_double(writer, "energy_standard_error",
                             statistics->energy_standard_error));
  RETURN_IF_ERROR(key_uint64(writer, "final_block_length",
                             (uint64_t)statistics->final_block_length));
  RETURN_IF_ERROR(key_double(writer, "final_effective_block_count",
                             statistics->final_effective_block_count));
  RETURN_IF_ERROR(key_double(writer, "final_tau_int_max",
                             statistics->final_tau_int_max));
  RETURN_IF_ERROR(key_double(writer, "final_tau_to_block_length",
                             statistics->final_tau_to_block_length));
  RETURN_IF_ERROR(key_double(
      writer, "maximum_leave_one_projective_distance",
      statistics->maximum_leave_one_projective_distance));
  RETURN_IF_ERROR(key_double(writer, "variance", statistics->variance));
  RETURN_IF_ERROR(key_double(writer, "variance_imaginary",
                             statistics->variance_imaginary));
  RETURN_IF_ERROR(key_double(writer, "variance_standard_error",
                             statistics->variance_standard_error));
  return mvmc_power_lanczos_json_object_end(writer);
}

static MVMCKrylovStatus render_record(
    const MVMCPowerLanczosProductionSession *session,
    const MVMCPowerLanczosProductionSessionSummary *summary,
    const MVMCPowerLanczosStabilizationResult *statistics,
    const MVMCPowerLanczosStabilizationOutputIdentity *identity,
    char *output, size_t output_capacity, size_t *output_size) {
  MVMCPowerLanczosJsonWriter writer;
  char snapshot_sha256[65];
  digest_hex(summary->snapshot_selection.snapshot_sha256, snapshot_sha256);
  RETURN_IF_ERROR(mvmc_power_lanczos_json_writer_initialize(
      output, output_capacity, &writer));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_object_begin(&writer));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_key(&writer, "blocks"));
  RETURN_IF_ERROR(write_blocks(&writer, session, summary->block_count));
  RETURN_IF_ERROR(key_string(
      &writer, "decision", mvmc_power_lanczos_stabilization_decision_string(
                                statistics->decision)));
  RETURN_IF_ERROR(
      key_uint64(&writer, "failed_gates", statistics->failed_gates));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_key(&writer, "gevp"));
  RETURN_IF_ERROR(write_gevp(&writer, &summary->coefficient_gevp));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_key(&writer, "identity"));
  RETURN_IF_ERROR(write_identity(&writer, identity));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_key(&writer, "policy"));
  RETURN_IF_ERROR(write_policy(&writer, summary));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_key(&writer, "resources"));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_object_begin(&writer));
  RETURN_IF_ERROR(key_uint64(&writer, "allocated_bytes",
                             (uint64_t)statistics->session_allocated_bytes));
  RETURN_IF_ERROR(key_uint64(&writer, "block_count",
                             (uint64_t)statistics->block_count));
  RETURN_IF_ERROR(key_uint64(&writer, "chain_count",
                             (uint64_t)statistics->chain_count));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_object_end(&writer));
  RETURN_IF_ERROR(key_uint64(
      &writer, "schema_version",
      MVMC_POWER_LANCZOS_STABILIZATION_OUTPUT_SCHEMA_VERSION));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_key(&writer, "session"));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_object_begin(&writer));
  RETURN_IF_ERROR(key_uint64(&writer, "basis_count",
                             (uint64_t)summary->basis_count));
  RETURN_IF_ERROR(key_uint64(&writer, "coefficient_accepted_steps_total",
                             summary->coefficient_accepted_steps));
  RETURN_IF_ERROR(key_uint64(&writer, "coefficient_proposals",
                             summary->coefficient_proposals));
  RETURN_IF_ERROR(key_uint64(&writer, "coefficient_sample_count",
                             summary->coefficient_sample_count));
  RETURN_IF_ERROR(key_uint64(
      &writer, "coefficient_scored_proposals_total",
      (summary->coefficient_proposals - summary->scale_pilot_proposals) *
          (uint64_t)summary->chain_count));
  RETURN_IF_ERROR(key_uint64(&writer, "execution_fingerprint",
                             summary->execution_fingerprint));
  RETURN_IF_ERROR(key_uint64(&writer, "final_accepted_steps_total",
                             summary->final_accepted_steps));
  RETURN_IF_ERROR(key_uint64(&writer, "final_policy_hash",
                             summary->final_policy_hash));
  RETURN_IF_ERROR(
      key_uint64(&writer, "final_proposals", summary->final_proposals));
  RETURN_IF_ERROR(key_uint64(&writer, "final_proposals_total",
                             summary->final_proposals *
                                 (uint64_t)summary->chain_count));
  RETURN_IF_ERROR(key_uint64(&writer, "final_sample_count",
                             summary->final_sample_count));
  RETURN_IF_ERROR(
      key_int64(&writer, "power_step", summary->power_step));
  RETURN_IF_ERROR(key_uint64(&writer, "scale_pilot_accepted_steps_total",
                             summary->scale_pilot_accepted_steps));
  RETURN_IF_ERROR(key_uint64(&writer, "scale_pilot_proposals",
                             summary->scale_pilot_proposals));
  RETURN_IF_ERROR(key_uint64(&writer, "scale_pilot_proposals_total",
                             summary->scale_pilot_proposals *
                                 (uint64_t)summary->chain_count));
  RETURN_IF_ERROR(key_uint64(&writer, "scale_pilot_sample_count_per_chain",
                             summary->scale_pilot_sample_count_per_chain));
  RETURN_IF_ERROR(
      key_string(&writer, "snapshot_sha256", snapshot_sha256));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_object_end(&writer));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_key(&writer, "statistics"));
  RETURN_IF_ERROR(write_statistics(&writer, statistics));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_key(&writer, "support"));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_object_begin(&writer));
  RETURN_IF_ERROR(key_uint64(&writer, "exact_zero_count",
                             statistics->exact_zero_primitive_count));
  RETURN_IF_ERROR(key_uint64(&writer, "finite_nonzero_count",
                             statistics->finite_nonzero_primitive_count));
  RETURN_IF_ERROR(key_double(&writer, "maximum_absolute_numeric_bound",
                             statistics->maximum_absolute_numeric_bound));
  RETURN_IF_ERROR(key_uint64(
      &writer, "nonzero_outside_support_contributions",
      statistics->nonzero_outside_support_contributions));
  RETURN_IF_ERROR(key_uint64(&writer, "numeric_zero_count",
                             statistics->numeric_zero_primitive_count));
  RETURN_IF_ERROR(key_uint64(&writer, "outside_base_support_samples",
                             statistics->outside_base_support_samples));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_object_end(&writer));
  RETURN_IF_ERROR(mvmc_power_lanczos_json_object_end(&writer));
  return mvmc_power_lanczos_json_finish(&writer, output_size);
}

MVMCKrylovStatus mvmc_power_lanczos_stabilization_output_render(
    const MVMCPowerLanczosProductionSession *session,
    const MVMCPowerLanczosStabilizationResult *statistics,
    const MVMCPowerLanczosStabilizationOutputIdentity *identity,
    char *output, size_t output_capacity, size_t *output_size) {
  MVMCPowerLanczosProductionSessionSummary summary;
  MVMCKrylovStatus status;
  if (output_size != NULL) *output_size = 0;
  if (output != NULL && output_capacity != 0) output[0] = '\0';
  if (session == NULL || statistics == NULL || !statistics->valid ||
      statistics->status != MVMC_KRYLOV_STATUS_OK ||
      statistics->version != MVMC_POWER_LANCZOS_STABILIZATION_RESULT_VERSION ||
      statistics->decision < MVMC_POWER_LANCZOS_STABILIZATION_PASS ||
      statistics->decision > MVMC_POWER_LANCZOS_STABILIZATION_FAIL ||
      !identity_valid(identity) || output == NULL || output_capacity < 2 ||
      output_size == NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = mvmc_power_lanczos_production_session_summary(session, &summary);
  if (status != MVMC_KRYLOV_STATUS_OK || !summary.valid ||
      summary.status != MVMC_KRYLOV_STATUS_OK ||
      summary.state != MVMC_POWER_LANCZOS_PRODUCTION_SESSION_FINALIZED ||
      summary.observable_enabled || !summary.scale_pilot_enabled ||
      summary.scale_pilot_warm_up == 0 ||
      summary.scale_pilot_sample_count_per_chain == 0 ||
      summary.scale_pilot_sample_count == 0 ||
      summary.eta_relative != 0x1p-40 || !isfinite(summary.resolved_eta) ||
      summary.resolved_eta <= 0.0 || summary.power_step < 1 ||
      summary.power_step > 2 || statistics->block_count != summary.block_count ||
      statistics->chain_count != summary.chain_count ||
      statistics->coefficient_sample_count != summary.coefficient_sample_count ||
      statistics->final_sample_count != summary.final_sample_count) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  status = render_record(session, &summary, statistics, identity, output,
                         output_capacity, output_size);
  if (status != MVMC_KRYLOV_STATUS_OK) {
    output[0] = '\0';
    *output_size = 0;
  }
  return status;
}

#undef RETURN_IF_ERROR
