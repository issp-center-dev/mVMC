/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#define MVMC_BOUNDED_KRYLOV_MARKOV_DRIVER_NO_MAIN
#include "bounded_krylov_markov_driver.c"

#include "krylov_gevp_solver.h"

#include <float.h>

enum {
  P4F_REPARSE_SCHEMA_VERSION = 1,
  P4F_FAMILY_COUNT = 4,
  P4F_UPPER_COUNT = 6,
  P4F_LINE_CAPACITY = 8192
};

typedef enum {
  P4F_FAMILY_S = 0,
  P4F_FAMILY_K,
  P4F_FAMILY_KR,
  P4F_FAMILY_B
} P4FMatrixFamily;

typedef struct {
  int schema;
  int site_count;
  int qp_total;
  size_t sample_count;
  size_t sector_size;
  int rank_count;
  size_t cache_bytes;
  double rho;
  double eta;
  size_t official_block_count;
  size_t official_block_length;
  char fixture[96];
} P4FTraceHeader;

typedef struct {
  P4FTraceHeader header;
  double norm[PROFILE_DEPTH_COUNT];
  double log_basis_scale[PROFILE_DEPTH_COUNT];
  int scale_seen[PROFILE_DEPTH_COUNT];
  uint64_t *configuration;
  double *denominator;
  size_t samples_seen;
  uint64_t sequence_fnv1a64;
  double complex original_exact[MARKOV_FAMILY_COUNT][P4F_UPPER_COUNT];
  double complex original_estimate[MARKOV_FAMILY_COUNT][P4F_UPPER_COUNT];
  int original_entry_seen[MARKOV_FAMILY_COUNT][P4F_UPPER_COUNT];
  double original_antihermitian_residual;
  int summary_seen;
} P4FTrace;

typedef struct {
  uint64_t sample_count;
  double denominator;
  double complex matrix[P4F_FAMILY_COUNT][P4F_UPPER_COUNT];
} P4FMatrixTotals;

static int p4f_line_complete(const char *line, FILE *stream) {
  const size_t length = strlen(line);
  return length > 0 &&
         (line[length - 1] == '\n' || feof(stream));
}

static int p4f_field(const char *line, const char *key, char *value,
                     size_t capacity) {
  const size_t key_length = strlen(key);
  const char *cursor = line;
  if (line == NULL || key == NULL || value == NULL || capacity == 0) {
    return 0;
  }
  while ((cursor = strstr(cursor, key)) != NULL) {
    const int left_boundary = cursor == line || cursor[-1] == ' ';
    const char *begin = cursor + key_length;
    size_t length = 0;
    if (left_boundary && *begin == '=') {
      ++begin;
      while (begin[length] != '\0' && begin[length] != '\n' &&
             begin[length] != ' ' && begin[length] != '\t' &&
             begin[length] != '\r') {
        ++length;
      }
      if (length == 0 || length >= capacity) return 0;
      memcpy(value, begin, length);
      value[length] = '\0';
      return 1;
    }
    cursor += key_length;
  }
  return 0;
}

static int p4f_parse_u64_field(const char *line, const char *key,
                                uint64_t *value) {
  char text[128];
  return p4f_field(line, key, text, sizeof(text)) &&
         parse_markov_u64(text, value);
}

static int p4f_parse_size_field(const char *line, const char *key,
                                size_t *value) {
  uint64_t parsed = 0;
  if (!p4f_parse_u64_field(line, key, &parsed) || parsed > SIZE_MAX) {
    return 0;
  }
  *value = (size_t)parsed;
  return 1;
}

static int p4f_parse_int_field(const char *line, const char *key,
                               int minimum, int maximum, int *value) {
  char text[128];
  return p4f_field(line, key, text, sizeof(text)) &&
         parse_int(text, minimum, maximum, value);
}

static int p4f_parse_double_field(const char *line, const char *key,
                                  double minimum, double maximum,
                                  double *value) {
  char text[128];
  return p4f_field(line, key, text, sizeof(text)) &&
         parse_markov_double(text, minimum, maximum, value);
}

static int p4f_close(double left, double right, double multiplier) {
  const double scale = fmax(1.0, fmax(fabs(left), fabs(right)));
  return isfinite(left) && isfinite(right) &&
         fabs(left - right) <= multiplier * DBL_EPSILON * scale;
}

static int p4f_close_complex(double complex left, double complex right,
                             double multiplier) {
  return p4f_close(creal(left), creal(right), multiplier) &&
         p4f_close(cimag(left), cimag(right), multiplier);
}

static int p4f_entry_index(size_t row, size_t column, size_t *entry) {
  return mvmc_krylov_streaming_upper_index(
             MARKOV_DIMENSION, row, column, entry) ==
         MVMC_KRYLOV_STATUS_OK;
}

static int p4f_parse_header(const char *line, P4FTraceHeader *header) {
  size_t expected_sector_size = 0;
  unsigned int *masks = NULL;
  size_t mask_count;
  if (strncmp(line, "MARKOV ", 7) != 0 || header == NULL) return 0;
  memset(header, 0, sizeof(*header));
  if (!p4f_field(line, "fixture", header->fixture,
                  sizeof(header->fixture)) ||
      (strcmp(header->fixture, "p4s9_long_direct_session_official") != 0 &&
       strcmp(header->fixture, "p4s_bounded_markov_real") != 0) ||
      !p4f_parse_int_field(line, "schema", 2, 4, &header->schema) ||
      !p4f_parse_int_field(line, "site_count", 4, 8,
                           &header->site_count) ||
      (header->site_count & 1) != 0 ||
      !p4f_parse_int_field(line, "qp_total", 1, 4,
                           &header->qp_total) ||
      (header->qp_total != 1 && header->qp_total != 4) ||
      !p4f_parse_size_field(line, "sample_count", &header->sample_count) ||
      header->sample_count == 0 ||
      !p4f_parse_size_field(line, "sector_size", &header->sector_size) ||
      !p4f_parse_int_field(line, "rank_count", 1, INT_MAX,
                           &header->rank_count) ||
      !p4f_parse_size_field(line, "cache_requested_bytes",
                            &header->cache_bytes) ||
      !p4f_parse_double_field(line, "rho", DBL_MIN, DBL_MAX,
                              &header->rho) ||
      !p4f_parse_double_field(line, "eta", DBL_MIN, DBL_MAX,
                              &header->eta) ||
      !p4f_parse_size_field(line, "official_block_count",
                            &header->official_block_count) ||
      !p4f_parse_size_field(line, "official_block_length",
                            &header->official_block_length) ||
      header->official_block_count != MARKOV_OFFICIAL_BLOCKS ||
      header->official_block_length == 0 ||
      header->official_block_count >
          SIZE_MAX / header->official_block_length ||
      header->sample_count !=
          header->official_block_count * header->official_block_length) {
    return 0;
  }
  mask_count = fixed_masks(header->site_count, &masks);
  if (mask_count != 0 &&
      checked_multiply_size(mask_count, mask_count,
                            &expected_sector_size) &&
      expected_sector_size == header->sector_size) {
    free(masks);
    return 1;
  }
  free(masks);
  return 0;
}

static int p4f_original_family(const char *name,
                               MarkovMatrixFamily *family) {
  if (strcmp(name, "S") == 0) {
    *family = MARKOV_FAMILY_S;
  } else if (strcmp(name, "K") == 0) {
    *family = MARKOV_FAMILY_K;
  } else if (strcmp(name, "B") == 0) {
    *family = MARKOV_FAMILY_B;
  } else {
    return 0;
  }
  return 1;
}

static int p4f_parse_trace(FILE *stream, P4FTrace *trace) {
  char line[P4F_LINE_CAPACITY];
  size_t scale_count = 0;
  if (stream == NULL || trace == NULL ||
      fgets(line, sizeof(line), stream) == NULL ||
      !p4f_line_complete(line, stream) ||
      !p4f_parse_header(line, &trace->header)) {
    return 0;
  }
  trace->configuration =
      (uint64_t *)calloc(trace->header.sample_count,
                         sizeof(*trace->configuration));
  trace->denominator =
      (double *)calloc(trace->header.sample_count,
                       sizeof(*trace->denominator));
  if (trace->configuration == NULL || trace->denominator == NULL) return 0;
  trace->sequence_fnv1a64 = UINT64_C(1469598103934665603);
  trace->original_antihermitian_residual = NAN;
  while (fgets(line, sizeof(line), stream) != NULL) {
    if (!p4f_line_complete(line, stream)) return 0;
    if (strncmp(line, "SCALE ", 6) == 0) {
      int order = -1;
      if (!p4f_parse_int_field(line, "order", 0, MARKOV_ORDER, &order) ||
          trace->scale_seen[order] ||
          !p4f_parse_double_field(line, "norm", DBL_MIN, DBL_MAX,
                                  &trace->norm[order]) ||
          !p4f_parse_double_field(line, "log_basis_scale", -DBL_MAX,
                                  DBL_MAX,
                                  &trace->log_basis_scale[order])) {
        return 0;
      }
      trace->scale_seen[order] = 1;
      ++scale_count;
    } else if (strncmp(line, "SAMPLE ", 7) == 0) {
      size_t sample = 0;
      uint64_t configuration = 0;
      double denominator = NAN;
      if (!p4f_parse_size_field(line, "sample", &sample) ||
          sample != trace->samples_seen ||
          sample >= trace->header.sample_count ||
          !p4f_parse_u64_field(line, "configuration", &configuration) ||
          !p4f_parse_double_field(line, "denominator", 0.0, DBL_MAX,
                                  &denominator)) {
        return 0;
      }
      trace->configuration[sample] = configuration;
      trace->denominator[sample] = denominator;
      markov_hash_u64(&trace->sequence_fnv1a64, (uint64_t)sample);
      markov_hash_u64(&trace->sequence_fnv1a64, configuration);
      markov_hash_u64(&trace->sequence_fnv1a64,
                      markov_double_bits(denominator));
      ++trace->samples_seen;
    } else if (strncmp(line, "SUMMARY ", 8) == 0) {
      if (trace->summary_seen ||
          !p4f_parse_double_field(
              line, "hamiltonian_antihermitian_residual", 0.0, DBL_MAX,
              &trace->original_antihermitian_residual)) {
        return 0;
      }
      trace->summary_seen = 1;
    } else if (strncmp(line, "ENTRY ", 6) == 0) {
      char name[16];
      MarkovMatrixFamily family;
      size_t row = 0;
      size_t column = 0;
      size_t entry = 0;
      double exact_re;
      double exact_im;
      double estimate_re;
      double estimate_im;
      if (!p4f_field(line, "family", name, sizeof(name)) ||
          !p4f_original_family(name, &family) ||
          !p4f_parse_size_field(line, "row", &row) ||
          !p4f_parse_size_field(line, "column", &column) ||
          row > column || column >= MARKOV_DIMENSION ||
          !p4f_entry_index(row, column, &entry) ||
          trace->original_entry_seen[family][entry] ||
          !p4f_parse_double_field(line, "exact_re", -DBL_MAX, DBL_MAX,
                                  &exact_re) ||
          !p4f_parse_double_field(line, "exact_im", -DBL_MAX, DBL_MAX,
                                  &exact_im) ||
          !p4f_parse_double_field(line, "estimate_re", -DBL_MAX, DBL_MAX,
                                  &estimate_re) ||
          !p4f_parse_double_field(line, "estimate_im", -DBL_MAX, DBL_MAX,
                                  &estimate_im)) {
        return 0;
      }
      trace->original_exact[family][entry] = exact_re + I * exact_im;
      trace->original_estimate[family][entry] =
          estimate_re + I * estimate_im;
      trace->original_entry_seen[family][entry] = 1;
    }
  }
  if (ferror(stream) || scale_count != PROFILE_DEPTH_COUNT ||
      trace->samples_seen != trace->header.sample_count ||
      !trace->summary_seen) {
    return 0;
  }
  {
    MarkovMatrixFamily family;
    for (family = MARKOV_FAMILY_S; family <= MARKOV_FAMILY_B; ++family) {
      size_t entry;
      for (entry = 0; entry < P4F_UPPER_COUNT; ++entry) {
        if (!trace->original_entry_seen[family][entry]) return 0;
      }
    }
  }
  return 1;
}

static const char *p4f_family_name(P4FMatrixFamily family) {
  switch (family) {
    case P4F_FAMILY_S:
      return "S";
    case P4F_FAMILY_K:
      return "K";
    case P4F_FAMILY_KR:
      return "KR";
    case P4F_FAMILY_B:
      return "B";
  }
  return "unknown";
}

static MVMCKrylovMatrixKind p4f_matrix_kind(P4FMatrixFamily family) {
  switch (family) {
    case P4F_FAMILY_S:
      return MVMC_KRYLOV_MATRIX_OVERLAP;
    case P4F_FAMILY_K:
      return MVMC_KRYLOV_MATRIX_HAMILTONIAN;
    case P4F_FAMILY_KR:
      return MVMC_KRYLOV_MATRIX_HAMILTONIAN_ADJOINT;
    case P4F_FAMILY_B:
      return MVMC_KRYLOV_MATRIX_HAMILTONIAN_SQUARED;
  }
  return MVMC_KRYLOV_MATRIX_OVERLAP;
}

static double complex *p4f_sample_family(
    P4FMatrixFamily family, double complex *overlap,
    double complex *hamiltonian, double complex *hamiltonian_reverse,
    double complex *hamiltonian_squared) {
  switch (family) {
    case P4F_FAMILY_S:
      return overlap;
    case P4F_FAMILY_K:
      return hamiltonian;
    case P4F_FAMILY_KR:
      return hamiltonian_reverse;
    case P4F_FAMILY_B:
      return hamiltonian_squared;
  }
  return NULL;
}

static int p4f_extract_totals(
    const MVMCKrylovMatrixMeasurementAccumulator *accumulator,
    P4FMatrixTotals *totals) {
  P4FMatrixFamily family;
  if (accumulator == NULL || totals == NULL) return 0;
  memset(totals, 0, sizeof(*totals));
  for (family = P4F_FAMILY_S; family <= P4F_FAMILY_B; ++family) {
    size_t entry;
    for (entry = 0; entry < P4F_UPPER_COUNT; ++entry) {
      size_t row = 0;
      size_t column = 0;
      MVMCKrylovJackknifeBlock block;
      if (!markov_entry_indices(P4F_UPPER_COUNT, entry, &row, &column) ||
          mvmc_krylov_matrix_measurement_extract_block(
              accumulator, p4f_matrix_kind(family), row, column,
              &block) != MVMC_KRYLOV_STATUS_OK) {
        return 0;
      }
      if (family == P4F_FAMILY_S && entry == 0) {
        totals->sample_count = block.sample_count;
        totals->denominator = block.denominator;
      } else if (block.sample_count != totals->sample_count ||
                 !p4f_close(block.denominator, totals->denominator,
                            32.0)) {
        return 0;
      }
      totals->matrix[family][entry] = block.numerator;
    }
  }
  return totals->sample_count != 0 && isfinite(totals->denominator) &&
         totals->denominator > 0.0;
}

static void p4f_project_upper(const double complex *source, int dimension,
                              double complex *destination) {
  size_t target = 0;
  int row;
  for (row = 0; row < dimension; ++row) {
    int column;
    for (column = row; column < dimension; ++column) {
      size_t source_entry = 0;
      (void)p4f_entry_index((size_t)row, (size_t)column, &source_entry);
      destination[target++] = source[source_entry];
    }
  }
}

static int p4f_solve_and_print(const char *scope, int dimension,
                               double cutoff,
                               const P4FMatrixTotals *totals) {
  double complex overlap[P4F_UPPER_COUNT] = {0.0};
  double complex forward[P4F_UPPER_COUNT] = {0.0};
  double complex reverse[P4F_UPPER_COUNT] = {0.0};
  double complex squared[P4F_UPPER_COUNT] = {0.0};
  MVMCKrylovGEVPPolicy policy;
  MVMCKrylovGEVPResult result;
  MVMCKrylovGEVPStatus status;
  const size_t count =
      (size_t)dimension * ((size_t)dimension + 1U) / 2U;
  int coefficient;
  p4f_project_upper(totals->matrix[P4F_FAMILY_S], dimension, overlap);
  p4f_project_upper(totals->matrix[P4F_FAMILY_K], dimension, forward);
  p4f_project_upper(totals->matrix[P4F_FAMILY_KR], dimension, reverse);
  p4f_project_upper(totals->matrix[P4F_FAMILY_B], dimension, squared);
  if (mvmc_krylov_gevp_default_policy(cutoff, &policy) !=
      MVMC_KRYLOV_GEVP_OK) {
    return 0;
  }
  status = mvmc_krylov_gevp_solve_complex_packed(
      &policy, dimension, overlap, forward, reverse, squared, count,
      &result);
  printf("GEVP scope=%s dimension=%d cutoff=%.17g status=%s valid=%d",
         scope, dimension, cutoff, mvmc_krylov_gevp_status_string(status),
         result.valid);
  printf(" retained_rank=%d discarded_rank=%d phase_pivot=%d",
         result.retained_rank, result.discarded_rank, result.phase_pivot);
  printf(" root_multiplicity=%d comparison_uses_projector=%d",
         result.root_multiplicity,
         result.coefficient_comparison_uses_projector);
  printf(" energy=%.17g energy_squared=%.17g variance=%.17g",
         result.energy, result.energy_squared, result.variance);
  printf(" residual=%.17g antihermitian_residual=%.17g",
         result.gevp_relative_residual,
         result.hamiltonian_antihermitian_residual);
  printf(" condition_estimate=%.17g root_gap=%.17g",
         result.condition_estimate, result.root_gap);
  for (coefficient = 0; coefficient < dimension; ++coefficient) {
    printf(" coefficient%d_re=%.17g coefficient%d_im=%.17g",
           coefficient, creal(result.coefficient[coefficient]),
           coefficient, cimag(result.coefficient[coefficient]));
  }
  for (coefficient = 0; coefficient < dimension * dimension;
       ++coefficient) {
    printf(" projector%d_re=%.17g projector%d_im=%.17g",
           coefficient,
           creal(result.root_subspace_projector[coefficient]),
           coefficient,
           cimag(result.root_subspace_projector[coefficient]));
  }
  fputc('\n', stdout);
  return status == MVMC_KRYLOV_GEVP_OK && result.valid;
}

static void p4f_add_totals(P4FMatrixTotals *destination,
                           const P4FMatrixTotals *source) {
  P4FMatrixFamily family;
  destination->sample_count += source->sample_count;
  destination->denominator += source->denominator;
  for (family = P4F_FAMILY_S; family <= P4F_FAMILY_B; ++family) {
    size_t entry;
    for (entry = 0; entry < P4F_UPPER_COUNT; ++entry) {
      destination->matrix[family][entry] +=
          source->matrix[family][entry];
    }
  }
}

static void p4f_subtract_totals(P4FMatrixTotals *destination,
                                const P4FMatrixTotals *source) {
  P4FMatrixFamily family;
  destination->sample_count -= source->sample_count;
  destination->denominator -= source->denominator;
  for (family = P4F_FAMILY_S; family <= P4F_FAMILY_B; ++family) {
    size_t entry;
    for (entry = 0; entry < P4F_UPPER_COUNT; ++entry) {
      destination->matrix[family][entry] -=
          source->matrix[family][entry];
    }
  }
}

/* The frozen P4-S9 schema-3 logs predate the P4-E full-Frobenius fix and
 * count each off-diagonal anti-Hermitian difference once.  Reconstruct that
 * historical diagnostic only for source-log identity checking; every new
 * P4-F diagnostic and gate uses the corrected full Frobenius convention. */
static double p4f_legacy_antihermitian_residual(
    const P4FMatrixTotals *totals) {
  double residual_norm = 0.0;
  double hamiltonian_norm = 0.0;
  size_t entry;
  for (entry = 0; entry < P4F_UPPER_COUNT; ++entry) {
    size_t row = 0;
    size_t column = 0;
    const double complex forward = totals->matrix[P4F_FAMILY_K][entry];
    const double complex reverse = totals->matrix[P4F_FAMILY_KR][entry];
    const double complex delta = forward - conj(reverse);
    (void)markov_entry_indices(P4F_UPPER_COUNT, entry, &row, &column);
    residual_norm = hypot(residual_norm, cabs(delta));
    hamiltonian_norm = hypot(hamiltonian_norm, cabs(forward));
    if (row != column) {
      hamiltonian_norm = hypot(hamiltonian_norm, cabs(reverse));
    }
  }
  return hamiltonian_norm > 0.0
             ? residual_norm / hamiltonian_norm
             : residual_norm;
}

static int p4f_emit_solver_results(const P4FMatrixTotals *full,
                                   const P4FMatrixTotals *exact,
                                   const P4FMatrixTotals *blocks,
                                   size_t block_count) {
  const double cutoff[3] = {1.0e-12, 1.0e-10, 1.0e-8};
  size_t cutoff_index;
  int dimension;
  for (dimension = 2; dimension <= 3; ++dimension) {
    for (cutoff_index = 0; cutoff_index < 3; ++cutoff_index) {
      if (!p4f_solve_and_print("full", dimension, cutoff[cutoff_index],
                               full) ||
          !p4f_solve_and_print("exact", dimension, cutoff[cutoff_index],
                               exact)) {
        return 0;
      }
    }
  }
  {
    const size_t prefix_blocks[3] = {4, 8, 16};
    size_t prefix_index;
    for (prefix_index = 0; prefix_index < 3; ++prefix_index) {
      const size_t count = prefix_blocks[prefix_index];
      P4FMatrixTotals prefix;
      size_t block;
      char scope[96];
      if (count > block_count) return 0;
      memset(&prefix, 0, sizeof(prefix));
      for (block = 0; block < count; ++block) {
        p4f_add_totals(&prefix, &blocks[block]);
      }
      for (dimension = 2; dimension <= 3; ++dimension) {
        snprintf(scope, sizeof(scope), "prefix%" PRIu64,
                 prefix.sample_count);
        if (!p4f_solve_and_print(scope, dimension, 1.0e-10, &prefix)) {
          return 0;
        }
        for (block = 0; block < count; ++block) {
          P4FMatrixTotals leave_one = prefix;
          p4f_subtract_totals(&leave_one, &blocks[block]);
          snprintf(scope, sizeof(scope), "prefix%" PRIu64 "_leave%zu",
                   prefix.sample_count, block);
          if (!p4f_solve_and_print(scope, dimension, 1.0e-10,
                                   &leave_one)) {
            return 0;
          }
        }
      }
    }
  }
  return 1;
}

static int p4f_reparse(FILE *stream) {
  P4FTrace trace;
  ProfileModel fixture;
  MVMCClassicKrylovModelWorkspace *model_workspace = NULL;
  ProfileScaledAmplitude *amplitude_workspace = NULL;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *krylov_workspace = NULL;
  MVMCKrylovBoundedCollectiveWorkspace *collective_workspace = NULL;
  const MVMCKrylovFockModel *model = NULL;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovMatrixMeasurementPolicy measurement_policy;
  MVMCKrylovMatrixMeasurementAccumulator total_accumulator;
  MVMCKrylovMatrixMeasurementBlockAccumulator block_accumulator;
  MVMCKrylovMatrixMeasurementDiagnosticsSummary diagnostics_summary;
  MVMCKrylovBoundedCollectiveResult collective_result;
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  MVMCKrylovStreamingComplexSum *total_overlap = NULL;
  MVMCKrylovStreamingComplexSum *total_hamiltonian = NULL;
  MVMCKrylovStreamingComplexSum *total_hamiltonian_reverse = NULL;
  MVMCKrylovStreamingComplexSum *total_squared = NULL;
  MVMCKrylovMatrixMeasurementAccumulator *block_storage = NULL;
  MVMCKrylovStreamingComplexSum *block_overlap = NULL;
  MVMCKrylovStreamingComplexSum *block_hamiltonian = NULL;
  MVMCKrylovStreamingComplexSum *block_hamiltonian_reverse = NULL;
  MVMCKrylovStreamingComplexSum *block_squared = NULL;
  unsigned int *masks = NULL;
  MVMCScaledComplex (*exact_values)[PROFILE_DEPTH_COUNT] = NULL;
  double complex (*raw_values)[PROFILE_DEPTH_COUNT] = NULL;
  uint64_t *configurations = NULL;
  size_t *configuration_lookup = NULL;
  double *slater = NULL;
  double complex *weights = NULL;
  double norms[PROFILE_DEPTH_COUNT];
  double log_basis_scale[PROFILE_DEPTH_COUNT];
  double eta = NAN;
  double projection_parameter = -0.27;
  size_t mask_count = 0;
  size_t sector_size = 0;
  size_t dimension = 0;
  size_t upper_count = 0;
  size_t lookup_count = 0;
  P4FMatrixTotals full;
  P4FMatrixTotals exact;
  P4FMatrixTotals blocks[MARKOV_OFFICIAL_BLOCKS];
  int return_code = 1;
  int all_match = 1;
  double maximum_denominator_error = 0.0;
  double maximum_entry_error = 0.0;
  double maximum_exact_error = 0.0;
  double legacy_antihermitian_residual = NAN;
  double source_antihermitian_residual = NAN;
  const char *source_antihermitian_convention = NULL;
  size_t sample;
  int order;

  memset(&trace, 0, sizeof(trace));
  memset(&fixture, 0, sizeof(fixture));
  memset(&full, 0, sizeof(full));
  memset(&exact, 0, sizeof(exact));
  memset(blocks, 0, sizeof(blocks));
  (void)&run_markov;
  (void)&parse_size_arg;
  if (world_size() != 1 || !p4f_parse_trace(stream, &trace)) {
    fprintf(stderr, "P4-F trace parse/preflight failed\n");
    goto cleanup;
  }
  mask_count = fixed_masks(trace.header.site_count, &masks);
  if (mask_count == 0 ||
      !checked_multiply_size(mask_count, mask_count, &sector_size) ||
      sector_size != trace.header.sector_size ||
      !initialize_model(&fixture, trace.header.site_count) ||
      !initialize_slater(trace.header.site_count, trace.header.qp_total,
                         &slater, &weights)) {
    fprintf(stderr, "P4-F fixture initialization failed\n");
    goto cleanup;
  }
  exact_values = (MVMCScaledComplex (*)[PROFILE_DEPTH_COUNT])calloc(
      sector_size, sizeof(*exact_values));
  raw_values = (double complex (*)[PROFILE_DEPTH_COUNT])calloc(
      sector_size, sizeof(*raw_values));
  configurations = (uint64_t *)calloc(sector_size, sizeof(*configurations));
  lookup_count = (size_t)1U << (2U * (unsigned int)trace.header.site_count);
  configuration_lookup =
      (size_t *)malloc(lookup_count * sizeof(*configuration_lookup));
  if (exact_values == NULL || raw_values == NULL || configurations == NULL ||
      configuration_lookup == NULL) {
    fprintf(stderr, "P4-F exact-sector allocation failed\n");
    goto cleanup;
  }
  for (sample = 0; sample < lookup_count; ++sample) {
    configuration_lookup[sample] = SIZE_MAX;
  }
  status = mvmc_classic_krylov_model_workspace_create_from_root(
      &fixture.raw, world_communicator(), &model_workspace);
  if (status != MVMC_KRYLOV_STATUS_OK ||
      (model = mvmc_classic_krylov_model(model_workspace)) == NULL ||
      !profile_limits(trace.header.cache_bytes, model, &limits) ||
      mvmc_bounded_krylov_plan_create(model, &limits, &plan) !=
          MVMC_KRYLOV_STATUS_OK ||
      mvmc_bounded_krylov_workspace_create(plan, &krylov_workspace) !=
          MVMC_KRYLOV_STATUS_OK ||
      mvmc_bounded_krylov_collective_workspace_create(
          world_communicator(), (size_t)trace.header.qp_total,
          PROFILE_DEPTH_COUNT, &collective_workspace) !=
          MVMC_KRYLOV_STATUS_OK ||
      !create_profile_amplitude(
          trace.header.site_count, trace.header.qp_total, slater, weights,
          projection_parameter, collective_workspace,
          &amplitude_workspace)) {
    fprintf(stderr, "P4-F exact-sector workspace setup failed\n");
    goto cleanup;
  }
  status = mvmc_bounded_krylov_collective_synchronize(
      collective_workspace, MVMC_KRYLOV_STATUS_OK,
      mvmc_bounded_krylov_plan_hash(plan), 1, &collective_result);
  if (status != MVMC_KRYLOV_STATUS_OK || !collective_result.valid ||
      !markov_evaluate_exact_sector(
          krylov_workspace, amplitude_workspace, masks, mask_count,
          sector_size, trace.header.site_count, exact_values, raw_values,
          configurations, norms) ||
      !markov_compute_scales_and_eta(
          sector_size,
          (const double complex(*)[PROFILE_DEPTH_COUNT])raw_values, norms,
          trace.header.rho, log_basis_scale, &eta) ||
      !markov_init_measurement_policy(eta, log_basis_scale,
                                      &measurement_policy) ||
      mvmc_krylov_matrix_measurement_dimension(
          MARKOV_ORDER, &dimension, &upper_count) !=
          MVMC_KRYLOV_STATUS_OK ||
      dimension != MARKOV_DIMENSION || upper_count != P4F_UPPER_COUNT) {
    fprintf(stderr, "P4-F exact-sector calibration failed\n");
    goto cleanup;
  }
  for (order = 0; order <= MARKOV_ORDER; ++order) {
    if (!p4f_close(norms[order], trace.norm[order], 8192.0) ||
        !p4f_close(log_basis_scale[order], trace.log_basis_scale[order],
                   8192.0)) {
      fprintf(stderr, "P4-F SCALE mismatch at order %d\n", order);
      goto cleanup;
    }
  }
  if (!p4f_close(eta, trace.header.eta, 8192.0)) {
    fprintf(stderr, "P4-F eta mismatch\n");
    goto cleanup;
  }
  for (sample = 0; sample < sector_size; ++sample) {
    const uint64_t configuration = configurations[sample];
    if (configuration >= lookup_count ||
        configuration_lookup[configuration] != SIZE_MAX) {
      fprintf(stderr, "P4-F configuration lookup construction failed\n");
      goto cleanup;
    }
    configuration_lookup[configuration] = sample;
  }

  total_overlap = (MVMCKrylovStreamingComplexSum *)calloc(
      upper_count, sizeof(*total_overlap));
  total_hamiltonian = (MVMCKrylovStreamingComplexSum *)calloc(
      upper_count, sizeof(*total_hamiltonian));
  total_hamiltonian_reverse = (MVMCKrylovStreamingComplexSum *)calloc(
      upper_count, sizeof(*total_hamiltonian_reverse));
  total_squared = (MVMCKrylovStreamingComplexSum *)calloc(
      upper_count, sizeof(*total_squared));
  if (total_overlap == NULL || total_hamiltonian == NULL ||
      total_hamiltonian_reverse == NULL || total_squared == NULL ||
      mvmc_krylov_matrix_measurement_accumulator_init_with_diagnostics(
          MARKOV_ORDER, total_overlap, total_hamiltonian,
          total_hamiltonian_reverse, total_squared, upper_count,
          &total_accumulator) != MVMC_KRYLOV_STATUS_OK ||
      !markov_init_block_accumulator(
          trace.header.official_block_count,
          trace.header.official_block_length, upper_count,
          &block_accumulator, &block_storage, &block_overlap,
          &block_hamiltonian, &block_hamiltonian_reverse, &block_squared)) {
    fprintf(stderr, "P4-F matrix accumulator initialization failed\n");
    goto cleanup;
  }

  for (sample = 0; sample < trace.header.sample_count; ++sample) {
    const uint64_t configuration = trace.configuration[sample];
    const size_t sector_index =
        configuration < lookup_count ? configuration_lookup[configuration]
                                     : SIZE_MAX;
    MVMCKrylovMatrixMeasurementSampleDiagnostics diagnostics;
    double error;
    if (sector_index == SIZE_MAX ||
        mvmc_krylov_matrix_measurement_accumulator_add_sample(
            &total_accumulator, &measurement_policy,
            exact_values[sector_index], PROFILE_DEPTH_COUNT,
            &diagnostics) != MVMC_KRYLOV_STATUS_OK ||
        !diagnostics.valid ||
        mvmc_krylov_matrix_measurement_block_accumulator_add_sample(
            &block_accumulator, &measurement_policy,
            exact_values[sector_index], PROFILE_DEPTH_COUNT, NULL) !=
            MVMC_KRYLOV_STATUS_OK) {
      fprintf(stderr, "P4-F sample reconstruction failed at %zu\n", sample);
      goto cleanup;
    }
    error = fabs(diagnostics.denominator - trace.denominator[sample]);
    maximum_denominator_error = fmax(maximum_denominator_error, error);
    if (!p4f_close(diagnostics.denominator, trace.denominator[sample],
                   8192.0)) {
      fprintf(stderr,
              "P4-F denominator mismatch at %zu: trace=%.17g reparse=%.17g\n",
              sample, trace.denominator[sample], diagnostics.denominator);
      goto cleanup;
    }
  }
  if (!p4f_extract_totals(&total_accumulator, &full)) {
    fprintf(stderr, "P4-F full matrix extraction failed\n");
    goto cleanup;
  }
  for (sample = 0; sample < trace.header.official_block_count; ++sample) {
    if (!p4f_extract_totals(&block_storage[sample], &blocks[sample])) {
      fprintf(stderr, "P4-F block matrix extraction failed at %zu\n",
              sample);
      goto cleanup;
    }
  }

  for (sample = 0; sample < sector_size; ++sample) {
    double complex overlap[P4F_UPPER_COUNT];
    double complex hamiltonian[P4F_UPPER_COUNT];
    double complex reverse[P4F_UPPER_COUNT];
    double complex squared[P4F_UPPER_COUNT];
    MVMCKrylovMatrixMeasurementSampleDiagnostics diagnostics;
    const MVMCKrylovStatus measurement_status =
        mvmc_krylov_matrix_measurement_sample_with_adjoint(
            &measurement_policy, exact_values[sample], PROFILE_DEPTH_COUNT,
            overlap, hamiltonian, reverse, squared, upper_count,
            &diagnostics);
    double guide;
    P4FMatrixFamily family;
    if (measurement_status != MVMC_KRYLOV_STATUS_OK ||
        !diagnostics.valid || !isfinite(diagnostics.log_guide) ||
        !isfinite(guide = exp(diagnostics.log_guide)) || guide <= 0.0) {
      fprintf(stderr, "P4-F exact matrix integration failed at %zu\n",
              sample);
      goto cleanup;
    }
    exact.denominator += diagnostics.denominator * guide;
    ++exact.sample_count;
    for (family = P4F_FAMILY_S; family <= P4F_FAMILY_B; ++family) {
      double complex *values = p4f_sample_family(
          family, overlap, hamiltonian, reverse, squared);
      size_t entry;
      for (entry = 0; entry < upper_count; ++entry) {
        exact.matrix[family][entry] += values[entry] * guide;
      }
    }
  }

  legacy_antihermitian_residual =
      p4f_legacy_antihermitian_residual(&full);
  if (mvmc_krylov_matrix_measurement_diagnostics_summary(
          &total_accumulator, 1.0e-12, 1.0,
          &diagnostics_summary) != MVMC_KRYLOV_STATUS_OK ||
      !diagnostics_summary.valid) {
    fprintf(stderr, "P4-F anti-Hermitian summary reconstruction failed\n");
    goto cleanup;
  }
  if (trace.header.schema == 3 &&
      strcmp(trace.header.fixture,
             "p4s9_long_direct_session_official") == 0) {
    source_antihermitian_convention = "legacy_upper_once";
    source_antihermitian_residual = legacy_antihermitian_residual;
  } else {
    source_antihermitian_convention = "corrected_full_frobenius";
    source_antihermitian_residual =
        diagnostics_summary.hamiltonian_antihermitian_residual;
  }
  if (!p4f_close(source_antihermitian_residual,
                 trace.original_antihermitian_residual, 8192.0)) {
    fprintf(stderr,
            "P4-F source anti-Hermitian summary mismatch: trace=%.17g "
            "source_reparse=%.17g legacy_reparse=%.17g "
            "corrected_reparse=%.17g delta=%.17g\n",
            trace.original_antihermitian_residual,
            source_antihermitian_residual,
            legacy_antihermitian_residual,
            diagnostics_summary.hamiltonian_antihermitian_residual,
            fabs(source_antihermitian_residual -
                 trace.original_antihermitian_residual));
    goto cleanup;
  }

  {
    MarkovMatrixFamily family;
    for (family = MARKOV_FAMILY_S; family <= MARKOV_FAMILY_B; ++family) {
      const P4FMatrixFamily p4f_family =
          family == MARKOV_FAMILY_S
              ? P4F_FAMILY_S
              : (family == MARKOV_FAMILY_K ? P4F_FAMILY_K : P4F_FAMILY_B);
      size_t entry;
      for (entry = 0; entry < upper_count; ++entry) {
        const double complex estimate =
            full.matrix[p4f_family][entry] / full.denominator;
        const double complex exact_estimate =
            exact.matrix[p4f_family][entry] / exact.denominator;
        maximum_entry_error =
            fmax(maximum_entry_error,
                 cabs(estimate - trace.original_estimate[family][entry]));
        maximum_exact_error =
            fmax(maximum_exact_error,
                 cabs(exact_estimate - trace.original_exact[family][entry]));
        if (!p4f_close_complex(estimate,
                               trace.original_estimate[family][entry],
                               8192.0) ||
            !p4f_close_complex(exact_estimate,
                               trace.original_exact[family][entry],
                               8192.0)) {
          all_match = 0;
        }
      }
    }
  }
  if (!all_match) {
    fprintf(stderr, "P4-F full ENTRY reconstruction mismatch\n");
    goto cleanup;
  }

  printf("REPARSE schema=%d source_schema=%d fixture=%s site_count=%d qp_total=%d",
         P4F_REPARSE_SCHEMA_VERSION, trace.header.schema,
         trace.header.fixture, trace.header.site_count,
         trace.header.qp_total);
  printf(" sample_count=%zu sector_size=%zu block_count=%zu",
         trace.header.sample_count, trace.header.sector_size,
         trace.header.official_block_count);
  printf(" block_length=%zu sequence_fnv1a64=%" PRIu64,
         trace.header.official_block_length, trace.sequence_fnv1a64);
  printf(" maximum_denominator_error=%.17g maximum_entry_error=%.17g",
         maximum_denominator_error, maximum_entry_error);
  printf(" maximum_exact_error=%.17g source_antihermitian_error=%.17g",
         maximum_exact_error,
         fabs(source_antihermitian_residual -
              trace.original_antihermitian_residual));
  printf(" source_antihermitian_convention=%s",
         source_antihermitian_convention);
  printf(" legacy_antihermitian_residual=%.17g",
         legacy_antihermitian_residual);
  printf(" corrected_antihermitian_residual=%.17g",
         diagnostics_summary.hamiltonian_antihermitian_residual);
  printf(" all_match=1\n");
  for (order = 0; order <= MARKOV_ORDER; ++order) {
    printf("REPARSE_SCALE order=%d norm=%.17g log_basis_scale=%.17g\n",
           order, norms[order], log_basis_scale[order]);
  }
  {
    P4FMatrixFamily family;
    for (family = P4F_FAMILY_S; family <= P4F_FAMILY_B; ++family) {
      size_t entry;
      for (entry = 0; entry < upper_count; ++entry) {
        size_t row = 0;
        size_t column = 0;
        (void)markov_entry_indices(upper_count, entry, &row, &column);
        printf("REPARSE_ENTRY family=%s row=%zu column=%zu",
               p4f_family_name(family), row, column);
        printf(" numerator_re=%.17g numerator_im=%.17g denominator=%.17g",
               creal(full.matrix[family][entry]),
               cimag(full.matrix[family][entry]), full.denominator);
        printf(" estimate_re=%.17g estimate_im=%.17g\n",
               creal(full.matrix[family][entry] / full.denominator),
               cimag(full.matrix[family][entry] / full.denominator));
      }
    }
    for (sample = 0; sample < trace.header.official_block_count; ++sample) {
      for (family = P4F_FAMILY_S; family <= P4F_FAMILY_B; ++family) {
        size_t entry;
        for (entry = 0; entry < upper_count; ++entry) {
          size_t row = 0;
          size_t column = 0;
          (void)markov_entry_indices(upper_count, entry, &row, &column);
          printf("REPARSE_BLOCK block=%zu family=%s row=%zu column=%zu",
                 sample, p4f_family_name(family), row, column);
          printf(" sample_count=%" PRIu64 " denominator=%.17g",
                 blocks[sample].sample_count, blocks[sample].denominator);
          printf(" numerator_re=%.17g numerator_im=%.17g\n",
                 creal(blocks[sample].matrix[family][entry]),
                 cimag(blocks[sample].matrix[family][entry]));
        }
      }
    }
  }
  if (!p4f_emit_solver_results(&full, &exact, blocks,
                               trace.header.official_block_count)) {
    fprintf(stderr, "P4-F GEVP solve failed\n");
    goto cleanup;
  }
  printf("REPARSE_DECISION decision=GO\n");
  return_code = 0;

cleanup:
  markov_free_block_storage(block_storage, block_overlap, block_hamiltonian,
                            block_hamiltonian_reverse, block_squared);
  free(total_squared);
  free(total_hamiltonian_reverse);
  free(total_hamiltonian);
  free(total_overlap);
  free(configuration_lookup);
  free(configurations);
  free(raw_values);
  free(exact_values);
  destroy_profile_amplitude(amplitude_workspace);
  mvmc_bounded_krylov_collective_workspace_destroy(collective_workspace);
  mvmc_bounded_krylov_workspace_destroy(krylov_workspace);
  mvmc_bounded_krylov_plan_destroy(plan);
  mvmc_classic_krylov_model_workspace_destroy(model_workspace);
  free(weights);
  free(slater);
  destroy_model(&fixture);
  free(masks);
  free(trace.denominator);
  free(trace.configuration);
  return return_code;
}

int main(int argc, char **argv) {
  FILE *stream = NULL;
  int result;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
#endif
  if (argc != 2) {
    if (world_rank() == 0) {
      fprintf(stderr, "usage: %s TRACE_LOG|-\n", argv[0]);
    }
    result = EXIT_FAILURE;
  } else if (world_rank() != 0) {
    result = EXIT_FAILURE;
  } else {
    stream = strcmp(argv[1], "-") == 0 ? stdin : fopen(argv[1], "r");
    if (stream == NULL) {
      fprintf(stderr, "cannot open P4-F trace: %s\n", argv[1]);
      result = EXIT_FAILURE;
    } else {
      result = p4f_reparse(stream) == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
      if (stream != stdin) fclose(stream);
    }
  }
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return result;
}
