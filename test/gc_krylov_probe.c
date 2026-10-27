/*
 * Enumeration probe for the GC first-step Lanczos oracle (Testing only).
 *
 * usage: gc_krylov_probe <input.txt> <output.tsv>
 *
 * The input is one fixed-order whitespace token stream written by
 * runtest_gc_lanczos.py:
 *   GC_LANCZOS_PROBE_V3
 *   L domain QP R active G S symmetry   spin support descriptor (symmetry
 *                                       0 general, 1 singlet, 2 triplet,
 *                                       3 zero; inactive is 0 1 0 0)
 *   R x (re im)                    raw translation weights
 *   QP x (re im)                   QP weights
 *   QP x (2L)^2 x (re im)          row-major QP Slater matrices
 *   T, T x (out_site out_spin in_site in_spin re im)
 *   I, I x (site value)            CoulombIntra
 *   E, E x (i j value)             CoulombInter
 *   H, H x (i j value)             Hund
 *   G J, G>0: L gutzwiller indices, J>0: L*L jastrow indices
 *   P, P x (re im)                 projection parameters
 *   END
 * Counts are validated, extra tokens and NaN/Inf values are rejected.  A
 * V2 stream (spin_triplet_s0 instead of the descriptor) is an explicit
 * error; it is never read as V3.
 * The output (rank 0) is one TSV: Q rows with the scaled q0..q2 of every
 * domain configuration, the guide-reweighted S/K/K^adj/T packed entries from
 * mvmc_krylov_matrix_measurement_sample_with_adjoint and the result of
 * mvmc_krylov_gevp_solve_complex_packed.  The probe runs the production GC
 * bridge, bounded engine, measurement and GEVP code; it never builds the
 * reference answer.
 */
#include "bounded_krylov_engine.h"
#include "krylov_gevp_solver.h"
#include "krylov_matrix_measurement.h"
#include "power_lanczos_gc_bridge.h"

#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

typedef struct {
  FILE *stream;
  int failed;
} Reader;

static int read_token(Reader *reader, char *buffer, size_t size) {
  if (reader->failed || fscanf(reader->stream, "%63s", buffer) != 1) {
    reader->failed = 1;
    return 0;
  }
  (void)size;
  return 1;
}

static long read_long(Reader *reader, long minimum, long maximum) {
  char token[64];
  char *end = NULL;
  long value;
  if (!read_token(reader, token, sizeof(token))) return 0;
  errno = 0;
  value = strtol(token, &end, 10);
  if (errno != 0 || end == token || *end != '\0' || value < minimum ||
      value > maximum) {
    reader->failed = 1;
    return 0;
  }
  return value;
}

static double read_double(Reader *reader) {
  char token[64];
  char *end = NULL;
  double value;
  if (!read_token(reader, token, sizeof(token))) return 0.0;
  errno = 0;
  value = strtod(token, &end);
  if (errno == ERANGE || end == token || *end != '\0' || !isfinite(value)) {
    reader->failed = 1;
    return 0.0;
  }
  return value;
}

static double complex read_complex(Reader *reader) {
  const double real = read_double(reader);
  const double imag = read_double(reader);
  return real + imag * I;
}

typedef struct {
  int site_count;
  int domain;
  int qp_total;
  int translation_count;
  MVMCGCSpinSupport spin_support;
  double complex *translation;
  double complex *qp_weights;
  double complex *slater;
  int transfer_count;
  int (*transfer_rows)[4];
  int **transfer_index;
  double complex *transfer_value;
  int intra_count;
  int *intra_index;
  double *intra_value;
  int inter_count;
  int (*inter_rows)[2];
  int **inter_index;
  double *inter_value;
  int hund_count;
  int (*hund_rows)[2];
  int **hund_index;
  double *hund_value;
  int ngutzwiller;
  int njastrow;
  int *gutzwiller_idx;
  int *jastrow_table;
  int **jastrow_idx;
  int nproj;
  double complex *projection;
} ProbeInput;

static void *checked_calloc(size_t count, size_t size) {
  return calloc(count == 0 ? 1 : count, size);
}

static int parse_pairs(Reader *reader, int count, int site_count,
                       int (**rows)[2], int ***index, double **value) {
  int k;
  *rows = (int(*)[2])checked_calloc((size_t)count, sizeof(**rows));
  *index = (int **)checked_calloc((size_t)count, sizeof(**index));
  *value = (double *)checked_calloc((size_t)count, sizeof(**value));
  if (*rows == NULL || *index == NULL || *value == NULL) return 0;
  for (k = 0; k < count; ++k) {
    (*rows)[k][0] = (int)read_long(reader, 0, site_count - 1);
    (*rows)[k][1] = (int)read_long(reader, 0, site_count - 1);
    (*value)[k] = read_double(reader);
    (*index)[k] = (*rows)[k];
  }
  return !reader->failed;
}

/* 1 parsed, 0 invalid, -1 obsolete V2 magic. */
static int parse_input(const char *path, ProbeInput *input) {
  Reader reader;
  char token[64];
  size_t orbital_count, matrix_count, k;
  int index;
  memset(input, 0, sizeof(*input));
  reader.stream = fopen(path, "r");
  reader.failed = reader.stream == NULL;
  if (reader.failed) return 0;
  if (!read_token(&reader, token, sizeof(token)) ||
      strcmp(token, "GC_LANCZOS_PROBE_V3") != 0) {
    const int obsolete = !reader.failed &&
                         strcmp(token, "GC_LANCZOS_PROBE_V2") == 0;
    fclose(reader.stream);
    return obsolete ? -1 : 0;
  }
  input->site_count = (int)read_long(&reader, 1, 16);
  input->domain = (int)read_long(&reader, 1, 2);
  input->qp_total = (int)read_long(&reader, 1, 4096);
  input->translation_count = (int)read_long(&reader, 1, 4096);
  /* Range only; the bridge validates the descriptor itself. */
  input->spin_support.active = (int)read_long(&reader, 0, 1);
  input->spin_support.gauss_count = (int)read_long(&reader, 1, 4096);
  input->spin_support.total_spin = (int)read_long(&reader, 0, 8);
  input->spin_support.pair_symmetry =
      (MVMCGCPairSymmetry)read_long(&reader, 0, 3);
  if (reader.failed) {
    fclose(reader.stream);
    return 0;
  }
  orbital_count = 2 * (size_t)input->site_count;
  matrix_count = orbital_count * orbital_count;
  input->translation = (double complex *)checked_calloc(
      (size_t)input->translation_count, sizeof(double complex));
  input->qp_weights = (double complex *)checked_calloc(
      (size_t)input->qp_total, sizeof(double complex));
  input->slater = (double complex *)checked_calloc(
      (size_t)input->qp_total * matrix_count, sizeof(double complex));
  if (input->translation == NULL || input->qp_weights == NULL ||
      input->slater == NULL) {
    fclose(reader.stream);
    return 0;
  }
  for (index = 0; index < input->translation_count; ++index) {
    input->translation[index] = read_complex(&reader);
  }
  for (index = 0; index < input->qp_total; ++index) {
    input->qp_weights[index] = read_complex(&reader);
  }
  for (k = 0; k < (size_t)input->qp_total * matrix_count; ++k) {
    input->slater[k] = read_complex(&reader);
  }
  input->transfer_count = (int)read_long(&reader, 0, 100000);
  input->transfer_rows = (int(*)[4])checked_calloc(
      (size_t)input->transfer_count, sizeof(*input->transfer_rows));
  input->transfer_index = (int **)checked_calloc(
      (size_t)input->transfer_count, sizeof(*input->transfer_index));
  input->transfer_value = (double complex *)checked_calloc(
      (size_t)input->transfer_count, sizeof(double complex));
  if (input->transfer_rows == NULL || input->transfer_index == NULL ||
      input->transfer_value == NULL) {
    reader.failed = 1;
  }
  for (index = 0; !reader.failed && index < input->transfer_count; ++index) {
    input->transfer_rows[index][0] =
        (int)read_long(&reader, 0, input->site_count - 1);
    input->transfer_rows[index][1] = (int)read_long(&reader, 0, 1);
    input->transfer_rows[index][2] =
        (int)read_long(&reader, 0, input->site_count - 1);
    input->transfer_rows[index][3] = (int)read_long(&reader, 0, 1);
    input->transfer_value[index] = read_complex(&reader);
    input->transfer_index[index] = input->transfer_rows[index];
  }
  input->intra_count = (int)read_long(&reader, 0, 100000);
  input->intra_index =
      (int *)checked_calloc((size_t)input->intra_count, sizeof(int));
  input->intra_value =
      (double *)checked_calloc((size_t)input->intra_count, sizeof(double));
  if (input->intra_index == NULL || input->intra_value == NULL) {
    reader.failed = 1;
  }
  for (index = 0; !reader.failed && index < input->intra_count; ++index) {
    input->intra_index[index] =
        (int)read_long(&reader, 0, input->site_count - 1);
    input->intra_value[index] = read_double(&reader);
  }
  input->inter_count = (int)read_long(&reader, 0, 100000);
  if (!reader.failed &&
      !parse_pairs(&reader, input->inter_count, input->site_count,
                   &input->inter_rows, &input->inter_index,
                   &input->inter_value)) {
    reader.failed = 1;
  }
  input->hund_count = (int)read_long(&reader, 0, 100000);
  if (!reader.failed &&
      !parse_pairs(&reader, input->hund_count, input->site_count,
                   &input->hund_rows, &input->hund_index,
                   &input->hund_value)) {
    reader.failed = 1;
  }
  input->ngutzwiller = (int)read_long(&reader, 0, 100000);
  input->njastrow = (int)read_long(&reader, 0, 100000);
  input->gutzwiller_idx = (int *)checked_calloc(
      (size_t)input->site_count, sizeof(int));
  input->jastrow_table = (int *)checked_calloc(
      (size_t)input->site_count * (size_t)input->site_count, sizeof(int));
  input->jastrow_idx =
      (int **)checked_calloc((size_t)input->site_count, sizeof(int *));
  if (input->gutzwiller_idx == NULL || input->jastrow_table == NULL ||
      input->jastrow_idx == NULL) {
    reader.failed = 1;
  }
  if (!reader.failed && input->ngutzwiller > 0) {
    for (index = 0; index < input->site_count; ++index) {
      input->gutzwiller_idx[index] =
          (int)read_long(&reader, 0, input->ngutzwiller - 1);
    }
  }
  if (!reader.failed && input->njastrow > 0) {
    for (k = 0; k < (size_t)input->site_count * (size_t)input->site_count;
         ++k) {
      input->jastrow_table[k] =
          (int)read_long(&reader, -1, input->njastrow - 1);
    }
  }
  for (index = 0; !reader.failed && index < input->site_count; ++index) {
    input->jastrow_idx[index] =
        input->jastrow_table + (size_t)index * (size_t)input->site_count;
  }
  input->nproj = (int)read_long(&reader, 0, 100000);
  input->projection = (double complex *)checked_calloc(
      (size_t)input->nproj, sizeof(double complex));
  if (input->projection == NULL) reader.failed = 1;
  for (index = 0; !reader.failed && index < input->nproj; ++index) {
    input->projection[index] = read_complex(&reader);
  }
  if (!reader.failed &&
      (!read_token(&reader, token, sizeof(token)) ||
       strcmp(token, "END") != 0 || fscanf(reader.stream, "%63s", token) == 1)) {
    reader.failed = 1;
  }
  fclose(reader.stream);
  return !reader.failed;
}

static int popcount64(uint64_t value) {
  int count = 0;
  while (value) {
    count += (int)(value & 1);
    value >>= 1;
  }
  return count;
}

static int allowed(uint64_t mask, int site_count, int domain) {
  const uint64_t low = (UINT64_C(1) << site_count) - 1;
  if (domain == (int)MVMC_KRYLOV_DOMAIN_GC_BALANCED) {
    return popcount64(mask & low) == popcount64(mask >> site_count);
  }
  return popcount64(mask) % 2 == 0;
}

static const char *state_name(MVMCScaledComplexState state) {
  switch (state) {
    case MVMC_SCALED_COMPLEX_FINITE_NONZERO:
      return "FINITE";
    case MVMC_SCALED_COMPLEX_EXACT_ZERO:
      return "EXACT_ZERO";
    case MVMC_SCALED_COMPLEX_NUMERIC_ZERO:
      return "NUMERIC_ZERO";
    case MVMC_SCALED_COMPLEX_NONFINITE:
      return "NONFINITE";
  }
  return "INVALID";
}

int main(int argc, char **argv) {
  ProbeInput input;
  MVMCPowerLanczosGCView view;
  MVMCPowerLanczosGCBridge *bridge = NULL;
  MVMCKrylovBoundedLimits limits;
  MVMCKrylovBoundedPlan *plan = NULL;
  MVMCKrylovBoundedWorkspace *engine = NULL;
  MVMCKrylovMatrixMeasurementPolicy policy;
  MVMCKrylovGEVPPolicy gevp_policy;
  MVMCKrylovGEVPResult gevp;
  MVMCKrylovGEVPStatus gevp_status;
  double complex overlap[3] = {0, 0, 0}, hamiltonian[3] = {0, 0, 0};
  double complex adjoint[3] = {0, 0, 0}, squared[3] = {0, 0, 0};
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  FILE *output = NULL;
  uint64_t mask;
  int rank = 0, size = 1, parsed, index, exit_code = 0;
  int measurement_ok = 1;
#ifdef _mpi_use
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  if (argc != 3) {
    if (rank == 0) fprintf(stderr, "usage: gc_krylov_probe in out\n");
#ifdef _mpi_use
    MPI_Finalize();
#endif
    return 2;
  }
  parsed = parse_input(argv[1], &input);
  if (parsed != 1) {
    if (rank == 0)
      fprintf(stderr, parsed < 0
                          ? "gc_krylov_probe: obsolete GC_LANCZOS_PROBE_V2 "
                            "input (regenerate as V3)\n"
                          : "gc_krylov_probe: invalid input\n");
#ifdef _mpi_use
    MPI_Finalize();
#endif
    return 2;
  }
  {
    /* Rank 0 owns the output; every rank learns whether it is writable. */
    int output_ok = 1;
    if (rank == 0) {
      output = fopen(argv[2], "w");
      output_ok = output != NULL;
    }
#ifdef _mpi_use
    MPI_Bcast(&output_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    if (!output_ok) {
      if (rank == 0) fprintf(stderr, "gc_krylov_probe: cannot write output\n");
#ifdef _mpi_use
      MPI_Finalize();
#endif
      return 2;
    }
  }
  memset(&view, 0, sizeof(view));
  view.data.site_count = input.site_count;
  view.data.arithmetic = MVMC_POWER_LANCZOS_CLASSIC_COMPLEX;
  view.data.transfer_count = input.transfer_count;
  view.data.transfer_indices = input.transfer_index;
  view.data.transfer_parameters = input.transfer_value;
  view.data.coulomb_intra_count = input.intra_count;
  view.data.coulomb_intra_indices = input.intra_index;
  view.data.coulomb_intra_parameters = input.intra_value;
  view.data.coulomb_inter_count = input.inter_count;
  view.data.coulomb_inter_indices = input.inter_index;
  view.data.coulomb_inter_parameters = input.inter_value;
  view.data.hund_count = input.hund_count;
  view.data.hund_indices = input.hund_index;
  view.data.hund_parameters = input.hund_value;
  view.data.qp_total = input.qp_total;
  view.data.qp_start = input.qp_total * rank / size;
  view.data.qp_end = input.qp_total * (rank + 1) / size;
  view.data.nproj = input.nproj;
  view.data.ngutzwiller_idx = input.ngutzwiller;
  view.data.njastrow_idx = input.njastrow;
  view.data.gutzwiller_idx = input.gutzwiller_idx;
  view.data.jastrow_idx = input.jastrow_idx;
  view.data.projection_parameters = input.projection;
  view.data.qp_weights = input.qp_weights;
  view.data.slater_complex = input.slater;
  view.domain = (MVMCKrylovConfigurationDomain)input.domain;
  view.translation_count = input.translation_count;
  view.translation_weights = input.translation;
  view.spin_support = input.spin_support;
#ifdef _mpi_use
  status = mvmc_power_lanczos_gc_bridge_create(&view, MPI_COMM_WORLD,
                                               MPI_COMM_WORLD, &bridge);
#else
  status = mvmc_power_lanczos_gc_bridge_create(&view, 0, 0, &bridge);
#endif
  if (status == MVMC_KRYLOV_STATUS_OK) {
    memset(&limits, 0, sizeof(limits));
    limits.amplitude_policy_hash =
        mvmc_power_lanczos_gc_bridge_generation_hash(bridge);
    limits.cache_bytes = (size_t)1 << 20;
    limits.max_row_transitions = 1u << 16;
    limits.max_workspace_bytes = (size_t)256 * 1024 * 1024;
    limits.max_node_expansions = UINT64_C(100000000);
    limits.max_terminal_amplitude_calls = UINT64_C(100000000);
    limits.max_total_row_transitions = UINT64_C(1000000000);
    limits.max_order = 2;
    status = mvmc_bounded_krylov_plan_create(
        mvmc_power_lanczos_gc_bridge_model(bridge), &limits, &plan);
  }
  if (status == MVMC_KRYLOV_STATUS_OK) {
    status = mvmc_bounded_krylov_workspace_create(plan, &engine);
  }
  memset(&policy, 0, sizeof(policy));
  policy.order = 2;
  policy.eta = 0x1p-40 * 3.0;
  for (index = 0; index <= 2; ++index) {
    policy.guide_lambda[index] = 1.0;
    policy.target_weight[index] = 1.0;
    policy.log_basis_scale[index] = 0.0;
  }
  if (output != NULL) {
    fprintf(output, "# GC_LANCZOS_PROBE_V3 L=%d domain=%d qp=%d active=%d "
            "G=%d S=%d symmetry=%d\n",
            input.site_count, input.domain, input.qp_total,
            input.spin_support.active, input.spin_support.gauss_count,
            input.spin_support.total_spin,
            (int)input.spin_support.pair_symmetry);
    fprintf(output, "BRIDGE\t%s\n", mvmc_krylov_status_string(status));
  }
  for (mask = 0; status == MVMC_KRYLOV_STATUS_OK &&
                 mask < (UINT64_C(1) << (2 * input.site_count));
       ++mask) {
    MVMCKrylovBoundedResult values;
    MVMCKrylovMatrixMeasurementSampleDiagnostics diagnostics;
    double complex sample_overlap[3], sample_hamiltonian[3];
    double complex sample_adjoint[3], sample_squared[3];
    MVMCKrylovStatus evaluated;
    if (!allowed(mask, input.site_count, input.domain)) continue;
    evaluated = mvmc_bounded_krylov_evaluate(
        engine, &mask, 1, mvmc_power_lanczos_gc_bridge_amplitude(bridge),
        mvmc_power_lanczos_gc_bridge_amplitude_context(bridge), &values);
    if (output != NULL) {
      fprintf(output, "Q\t%llu\t%d\t%s", (unsigned long long)mask,
              popcount64(mask), mvmc_krylov_status_string(evaluated));
      for (index = 0; evaluated == MVMC_KRYLOV_STATUS_OK && index <= 2;
           ++index) {
        fprintf(output, "\t%s\t%.17e\t%.17e\t%.17e",
                state_name(values.value[index].state),
                creal(values.value[index].phase),
                cimag(values.value[index].phase),
                values.value[index].log_abs);
      }
      fputc('\n', output);
    }
    if (evaluated != MVMC_KRYLOV_STATUS_OK) {
      status = evaluated;
      break;
    }
    if (measurement_ok &&
        mvmc_krylov_matrix_measurement_sample_with_adjoint(
            &policy, values.value, 3, sample_overlap, sample_hamiltonian,
            sample_adjoint, sample_squared, 3, &diagnostics) !=
            MVMC_KRYLOV_STATUS_OK) {
      /* Keep enumerating q values; only the S/K/T output is dropped. */
      measurement_ok = 0;
    }
    /* Reweight by the guide: sum_x g(x) * contribution(x) / g(x). */
    for (index = 0; measurement_ok && index < 3; ++index) {
      const double weight = exp(diagnostics.log_guide);
      overlap[index] += weight * sample_overlap[index];
      hamiltonian[index] += weight * sample_hamiltonian[index];
      adjoint[index] += weight * sample_adjoint[index];
      squared[index] += weight * sample_squared[index];
    }
  }
  if (output != NULL) {
    fprintf(output, "STATUS\t%s\n", mvmc_krylov_status_string(status));
    fprintf(output, "MEASUREMENT\t%s\n", measurement_ok ? "ok" : "failed");
    if (status == MVMC_KRYLOV_STATUS_OK && measurement_ok) {
      for (index = 0; index < 3; ++index) {
        fprintf(output, "M\t%d\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e"
                        "\t%.17e\t%.17e\n",
                index, creal(overlap[index]), cimag(overlap[index]),
                creal(hamiltonian[index]), cimag(hamiltonian[index]),
                creal(adjoint[index]), cimag(adjoint[index]),
                creal(squared[index]), cimag(squared[index]));
      }
      mvmc_krylov_gevp_default_policy(0x1p-40, &gevp_policy);
      gevp_status = mvmc_krylov_gevp_solve_complex_packed(
          &gevp_policy, 2, overlap, hamiltonian, adjoint, squared, 3, &gevp);
      fprintf(output,
              "GEVP\t%s\t%d\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e\t%.17e"
              "\t%.17e\n",
              mvmc_krylov_gevp_status_string(gevp_status),
              gevp.retained_rank, gevp.energy, gevp.variance,
              creal(gevp.coefficient[0]), cimag(gevp.coefficient[0]),
              creal(gevp.coefficient[1]), cimag(gevp.coefficient[1]),
              gevp.gevp_relative_residual);
    }
    fclose(output);
  }
  if (status != MVMC_KRYLOV_STATUS_OK) exit_code = 1;
  mvmc_bounded_krylov_workspace_destroy(engine);
  mvmc_bounded_krylov_plan_destroy(plan);
  mvmc_power_lanczos_gc_bridge_destroy(bridge);
#ifdef _mpi_use
  MPI_Finalize();
#endif
  return exit_code;
}
