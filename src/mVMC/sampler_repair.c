/* Selective post-accept repair for the real, non-block first-Lanczos sampler.
 * The growth rules and their ordering match the validated P1c/post/min scheme.
 * Included by vmcmain.h after matrix.c; no diagnostic proposal replay is needed.
 */
#ifndef MVMC_SAMPLER_REPAIR_C
#define MVMC_SAMPLER_REPAIR_C
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SAMPLER_REPAIR_GROWTH 1000.0
#define SAMPLER_REPAIR_EVENT_LIMIT_PER_KIND 4
static struct {
  int enabled, logging, bin, n, proposal_finite;
  double *reference, *floor, *old_pf, *new_pf, *k, *old_k, *new_k;
  double *buffer, *work;
  int *iwork, *selected;
  FILE *file;
  long samples, nonfinite, gt6, gt3, gt1, accepts, refreshed, components, total;
  double max_d;
  int out_step, in_step, move;
  long proposal_zero, ratio_underflow, ratio_nonfinite, current_nonfinite;
  long component_zero_pre, component_zero_post, measurement_zero, repair_failed;
  long measure_weight_zero, skip_weight, skip_energy, measure_factor_failed;
  double measure_log, measure_projection, measure_stored;
  int measure_issue;
  long recovery_full, recovery_restored, recovery_rollback, recovery_failed;
  long events_written, events_suppressed;
  int events_by_kind[4];
} Sr;

static void SamplerRepairFail(const char *message);
static void SamplerRepairFlush(void);

/* Bounded exceptional-event log. Scalars are already available to the solver.
 * Values a..f depend on event kind and are documented with the output format. */
static void SamplerRepairEvent(const char *kind, int qp, int mask,
    double a, double b, double c, double d, double e, double f,
    int info, const char *action) {
  const int category = strcmp(kind, "proposal") == 0 ? 0 :
      strcmp(kind, "component") == 0 ? 1 : strcmp(kind, "current") == 0 ? 2 : 3;
  if (!Sr.file) return;
  if (Sr.events_by_kind[category] >= SAMPLER_REPAIR_EVENT_LIMIT_PER_KIND) {
    Sr.events_suppressed++;
    return;
  }
  fprintf(Sr.file, "E %d %d %d %d %s %d %d %.17e %.17e %.17e %.17e %.17e %.17e %d %s\n",
          Sr.bin, Sr.out_step, Sr.in_step, Sr.move, kind, qp, mask,
          a, b, c, d, e, f, info, action);
  Sr.events_written++;
  Sr.events_by_kind[category]++;
  if (fflush(Sr.file) != 0 || ferror(Sr.file))
    SamplerRepairFail("cannot write exceptional-event detail");
}

static void SamplerRepairFail(const char *message) {
  fprintf(stderr, "Error: sampler repair: %s\n", message);
  if (Sr.file) fflush(Sr.file);
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  exit(EXIT_FAILURE);
}

/* Preserve the normalization and exceptional-value handling of the validated
 * candidate. No sample or estimator is clipped or removed by this operation. */
static void SamplerRepairKappa(const double *pf, double *k) {
  int q;
  double maximum = 0.0;
  for (q = 0; q < Sr.n; ++q)
    if (fabs(pf[q]) > maximum) maximum = fabs(pf[q]);
  for (q = 0; q < Sr.n; ++q) {
    double value = maximum > 0.0 && isfinite(maximum)
        ? fabs(pf[q]) / maximum : 0.0;
    k[q] = isfinite(value) ? value : 0.0;
  }
}

static int SamplerRepairSupported(void) {
#ifdef _pf_block_update
  return 0;
#else
  return NVMCCalMode == 1 && NLanczosMode == 1 && NLanczosStep == 1 &&
      NLanczosEstimatorMode == 0 && AllComplexFlag == 0 &&
      iFlgOrbitalGeneral == 0 && NBackFlowIdx == 0 && !FlagGrandCanonical &&
      !FlagRBM && NSplitSize == 1 && NExUpdatePath == 1;
#endif
}

static int SamplerRepairSwitch(const char *name, int default_value) {
  const char *value = getenv(name);
  if (value == NULL || *value == '\0') return default_value;
  if (strcmp(value, "0") == 0) return 0;
  if (strcmp(value, "1") == 0) return 1;
  SamplerRepairFail("MVMC_SAMPLER_REPAIR and MVMC_SAMPLER_DRIFT_LOG accept only 0 or 1");
  return 0;
}

static void *SamplerRepairAlloc(size_t count, size_t size) {
  void *result = calloc(count, size);
  if (result == NULL) SamplerRepairFail("workspace allocation failed");
  return result;
}

static void SamplerRepairInit(void) {
  int rank, settings[2];
  const int supported = SamplerRepairSupported();
  char path[64];
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  /* Read configuration once, from world rank zero, before any sampler work. */
  if (rank == 0) {
    settings[0] = SamplerRepairSwitch("MVMC_SAMPLER_REPAIR", -1);
    settings[1] = SamplerRepairSwitch("MVMC_SAMPLER_DRIFT_LOG", -1);
  }
  MPI_Bcast(settings, 2, MPI_INT, 0, MPI_COMM_WORLD);
  if (settings[0] == 1 && !supported)
    SamplerRepairFail("enabled outside real non-FSZ/no-BF/no-RBM canonical, non-block, split1, hop/exchange, unguided legacy first-Lanczos measurement");
  Sr.enabled = supported && settings[0] != 0;
  Sr.logging = settings[1] < 0 ? Sr.enabled : settings[1];
  if (Sr.logging && !supported)
    SamplerRepairFail("drift logging requested outside the supported measurement path");
  Sr.bin = Sr.out_step = Sr.in_step = -1;
  Sr.n = NQPFull;
  if (rank == 0)
    fprintf(stdout, "Sampler component repair: %s; growth=1000; drift_log=%d\n",
            Sr.enabled ? "on" : "off", Sr.logging);
  if (Sr.enabled) {
    if (NQPFull <= 0 || Nsize <= 0 || LapackLWork <= 0)
      SamplerRepairFail("invalid workspace dimensions");
    Sr.reference = SamplerRepairAlloc(Sr.n, sizeof(double));
    Sr.floor = SamplerRepairAlloc(Sr.n, sizeof(double));
    Sr.old_pf = SamplerRepairAlloc(Sr.n, sizeof(double));
    Sr.new_pf = SamplerRepairAlloc(Sr.n, sizeof(double));
    Sr.k = SamplerRepairAlloc(Sr.n, sizeof(double));
    Sr.old_k = SamplerRepairAlloc(Sr.n, sizeof(double));
    Sr.new_k = SamplerRepairAlloc(Sr.n, sizeof(double));
    Sr.selected = SamplerRepairAlloc(Sr.n, sizeof(int));
    Sr.buffer = SamplerRepairAlloc((size_t)Nsize, (size_t)Nsize * sizeof(double));
    Sr.work = SamplerRepairAlloc((size_t)LapackLWork, sizeof(double));
    Sr.iwork = SamplerRepairAlloc((size_t)Nsize, sizeof(int));
  }
  if (Sr.logging) {
    snprintf(path, sizeof(path), "sampler_drift_r%04d.dat", rank);
    Sr.file = fopen(path, "w");
    if (Sr.file == NULL) SamplerRepairFail("cannot open drift summary");
    fprintf(Sr.file, "# Z bin proposalZero ratioUnderflow ratioNonfinite currentNonfinite componentZeroPre componentZeroPost measurementZero repairFailed measureWeightZero skipWeight skipEnergy measureFactorFailed detailsTotal suppressedTotal\n");
    fprintf(Sr.file, "# E bin out in move kind qp mask a b c d e f info action; max 4 per kind/rank/run\n");
    fprintf(Sr.file, "# R bin fullRebuild recovered rollback failed\n");
    fprintf(Sr.file, "# Q bin samples nonfinite maxAbsDrift gt1e6 gt1e3 gt1 accepts refreshAccepts components\n");
  }
}

static void SamplerRepairFlush(void) {
  if (Sr.file && Sr.bin >= 0) {
    fprintf(Sr.file, "Q %d %ld %ld %.17e %ld %ld %ld %ld %ld %ld\n", Sr.bin,
            Sr.samples, Sr.nonfinite, Sr.max_d, Sr.gt6, Sr.gt3, Sr.gt1,
            Sr.accepts, Sr.refreshed, Sr.components);
    if (fflush(Sr.file) != 0 || ferror(Sr.file))
      SamplerRepairFail("cannot write drift summary");
  }
  if (Sr.file && Sr.bin >= 0) {
    fprintf(Sr.file, "Z %d %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld\n", Sr.bin,
            Sr.proposal_zero, Sr.ratio_underflow, Sr.ratio_nonfinite,
            Sr.current_nonfinite, Sr.component_zero_pre, Sr.component_zero_post,
            Sr.measurement_zero, Sr.repair_failed, Sr.measure_weight_zero, Sr.skip_weight,
            Sr.skip_energy, Sr.measure_factor_failed, Sr.events_written, Sr.events_suppressed);
    if (fflush(Sr.file) != 0 || ferror(Sr.file))
      SamplerRepairFail("cannot write zero-event summary");
  }
  if (Sr.file && Sr.bin >= 0) {
    fprintf(Sr.file, "R %d %ld %ld %ld %ld\n", Sr.bin, Sr.recovery_full,
            Sr.recovery_restored, Sr.recovery_rollback, Sr.recovery_failed);
    if (fflush(Sr.file) != 0 || ferror(Sr.file))
      SamplerRepairFail("cannot write recovery summary");
  }
  Sr.recovery_full = Sr.recovery_restored = Sr.recovery_rollback = Sr.recovery_failed = 0;
  Sr.proposal_zero = Sr.ratio_underflow = Sr.ratio_nonfinite = Sr.current_nonfinite = 0;
  Sr.component_zero_pre = Sr.component_zero_post = Sr.measurement_zero = Sr.repair_failed = 0;
  Sr.measure_weight_zero = Sr.skip_weight = Sr.skip_energy = Sr.measure_factor_failed = 0;
  /* Detail counts and budget are cumulative for the entire run, not per bin. */
  Sr.samples = Sr.nonfinite = Sr.gt6 = Sr.gt3 = Sr.gt1 = 0;
  Sr.accepts = Sr.refreshed = Sr.components = 0;
  Sr.max_d = 0.0;
}

static void SamplerRepairSetBin(int bin) {
  SamplerRepairFlush();
  Sr.bin = bin;
  Sr.out_step = Sr.in_step = -1; Sr.move = 0;
}

static void SamplerRepairRecompute(void) {
  if (!Sr.enabled) return;
  SamplerRepairKappa(PfM_real, Sr.reference);
  memcpy(Sr.floor, Sr.reference, (size_t)Sr.n * sizeof(double));
}

static void SamplerRepairProposal(int out_step, int in_step, int move,
    double old_log, double proposed_log, double projection_ratio, double ratio, int accepted) {
  Sr.out_step = out_step; Sr.in_step = in_step; Sr.move = move;
  if (!Sr.logging) return;
  Sr.proposal_zero += proposed_log == -INFINITY;
  Sr.ratio_underflow += ratio == 0.0 && isfinite(old_log) &&
      isfinite(proposed_log) && isfinite(projection_ratio);
  Sr.ratio_nonfinite += !isfinite(ratio);
  if (proposed_log == -INFINITY || ratio == 0.0 || !isfinite(ratio))
    SamplerRepairEvent("proposal", -1, 0, old_log, proposed_log, projection_ratio,
                      ratio, 0.0, 0.0, 0, accepted ? "accept-log" : "reject");
}

static void SamplerRepairCurrent(int out_step, int in_step, int move,
    double before, double after, const char *action) {
  Sr.out_step = out_step; Sr.in_step = in_step; Sr.move = move;
  if (!Sr.logging || (isfinite(before) && isfinite(after))) return;
  Sr.current_nonfinite += !isfinite(after);
  SamplerRepairEvent("current", -1, 0, before, after, 0.0, 0.0, 0.0, 0.0, 0, action);
}

/* Compare in log space: positive overflow is acceptance, not NaN rejection.
 * The caller draws one uniform variate even for deterministic decisions. */
static int SamplerRepairLogAccept(double old_log, double proposed_log,
                                  double projection_ratio, double uniform) {
  double half_log_ratio;
  if (!isfinite(old_log) || !isfinite(projection_ratio) ||
      !isfinite(proposed_log)) return 0;
  half_log_ratio = projection_ratio + (proposed_log - old_log);
  if (isnan(half_log_ratio)) return 0;
  if (half_log_ratio >= 0.0) return 1;
  /* Finite input logs represent positive weights even if their subtraction
   * overflows. genrand_real2 can return zero; the true zero was rejected above. */
  if (uniform == 0.0) return 1;
  return 0.5 * log(uniform) < half_log_ratio;
}

static void SamplerRepairRecoveryFailure(double before, double after, int info,
                                         const char *action) {
  int rank;
  Sr.recovery_failed++;
  SamplerRepairEvent("current", -1, 0, before, after, NAN, NAN, NAN, NAN, info, action);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  fprintf(stderr, "sampler recovery failure rank=%d bin=%d out=%d in=%d move=%d info=%d log_before=%.17e log_after=%.17e action=%s\n",
          rank, Sr.bin, Sr.out_step, Sr.in_step, Sr.move, info, before, after, action);
  SamplerRepairFlush();
  SamplerRepairFail("cannot establish a finite current amplitude");
}

/* Used after an already-completed full factorization: never retry it blindly. */
static void SamplerRepairRequireCurrent(double log_value, int info,
                                        const char *action) {
  if (Sr.enabled && (info != 0 || !isfinite(log_value)))
    SamplerRepairRecoveryFailure(log_value, log_value, info, action);
}

/* Only the exceptional accepted state is rebuilt. No snapshot of all inverse
 * matrices is needed: projection counts have not been committed yet, and the
 * old configuration is reconstructed by undoing the hop(s). Returns commit. */
static int SamplerRepairResolveAccepted(int *eleIdx, int *eleCfg, int *eleNum,
    int mi, int mj, int ri, int rj, int spin, int move,
    int qpStart, int qpEnd, MPI_Comm comm, double *new_log, double *old_log) {
  int info;
  double before;
  if (!Sr.enabled || isfinite(*new_log)) return 1;
  before = *new_log;
  Sr.current_nonfinite++;
  Sr.recovery_full++;
  info = CalculateMAll_real(eleIdx, qpStart, qpEnd);
  *new_log = info == 0 ? CalculateLogIP_real(PfM_real, qpStart, qpEnd, comm) : NAN;
  if (info != 0 || isnan(*new_log) || *new_log == INFINITY)
    SamplerRepairRecoveryFailure(before, *new_log, info, "abort-proposal-rebuild");
  if (isfinite(*new_log)) {
    SamplerRepairRecompute();
    Sr.recovery_restored++;
    SamplerRepairEvent("current", -1, 0, before, *new_log, NAN, NAN, NAN, NAN,
                      0, "recovered-accepted");
    return 1;
  }
  /* Directly evaluated zero has no target weight. Undo the provisional accept;
   * a failed factorization is NOT used as evidence that this weight is zero. */
  SamplerRepairEvent("current", -1, 0, before, *new_log, NAN, NAN, NAN, NAN,
                    0, "rollback-zero");
  if (move == 2)
    revertEleConfig(mj, rj, ri, 1-spin, eleIdx, eleCfg, eleNum);
  revertEleConfig(mi, ri, rj, spin, eleIdx, eleCfg, eleNum);
  before = *old_log;
  Sr.recovery_full++;
  info = CalculateMAll_real(eleIdx, qpStart, qpEnd);
  *old_log = info == 0 ? CalculateLogIP_real(PfM_real, qpStart, qpEnd, comm) : NAN;
  if (info != 0 || !isfinite(*old_log))
    SamplerRepairRecoveryFailure(before, *old_log, info, "abort-rollback");
  SamplerRepairRecompute();
  Sr.recovery_rollback++;
  SamplerRepairEvent("current", -1, 0, before, *old_log, NAN, NAN, NAN, NAN,
                    0, "restored-old-reject");
  return 0;
}

/* Called only on an accepted proposal, before the fast inverse update. */
static void SamplerRepairAcceptPre(const double *proposed) {
  if (!Sr.enabled) return;
  memcpy(Sr.old_pf, PfM_real, (size_t)Sr.n * sizeof(double));
  memcpy(Sr.new_pf, proposed, (size_t)Sr.n * sizeof(double));
}

/* Called after the accepted hop/exchange update. The caller updates logIpNew
 * after any repair; the accept decision and RNG stream have already been set. */
static int SamplerRepairAcceptPost(const int *eleIdx, int qpStart, int qpEnd) {
  int q, any = 0;
  if (!Sr.enabled) return 0;
  SamplerRepairKappa(PfM_real, Sr.k);
  Sr.proposal_finite = 1;
  for (q = 0; q < Sr.n; ++q)
    if (!isfinite(Sr.old_pf[q]) || !isfinite(Sr.new_pf[q])) Sr.proposal_finite = 0;
  SamplerRepairKappa(Sr.old_pf, Sr.old_k);
  SamplerRepairKappa(Sr.new_pf, Sr.new_k);
  if (Sr.proposal_finite)
    for (q = 0; q < Sr.n; ++q)
      if (Sr.old_k[q] < Sr.floor[q]) Sr.floor[q] = Sr.old_k[q];
  for (q = 0; q < Sr.n; ++q) {
    const int grown = Sr.reference[q] <= 0.0 ||
        Sr.k[q] > SAMPLER_REPAIR_GROWTH * Sr.reference[q];
    const int post = Sr.proposal_finite &&
        Sr.new_k[q] > SAMPLER_REPAIR_GROWTH * Sr.old_k[q];
    const int minimum = Sr.proposal_finite &&
        Sr.new_k[q] > SAMPLER_REPAIR_GROWTH * Sr.floor[q];
    Sr.selected[q] = grown || post || minimum;
    if (Sr.selected[q]) {
      const double before = PfM_real[q];
      const int mask = grown + 2 * post + 4 * minimum;
      int info = calculateMAll_child_real(eleIdx, qpStart, qpEnd, q, Sr.buffer,
                                         Sr.iwork, Sr.work, LapackLWork,
                                         PfM_real, InvM_real);
      if (Sr.logging) {
        Sr.component_zero_pre += before == 0.0;
        Sr.component_zero_post += info == 0 && PfM_real[q] == 0.0;
        Sr.repair_failed += info != 0;
        if (before == 0.0 || Sr.old_pf[q] == 0.0 || PfM_real[q] == 0.0 ||
            !isfinite(before) || !isfinite(PfM_real[q]) || info != 0)
          SamplerRepairEvent("component", q, mask, Sr.old_pf[q], Sr.new_pf[q],
                            before, PfM_real[q], Sr.reference[q], Sr.floor[q],
                            info, info == 0 ? "refresh" : "abort");
      }
      if (info != 0) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        fprintf(stderr, "sampler repair failure rank=%d bin=%d out=%d in=%d qp=%d info=%d pf_before=%.17e pf_after=%.17e\n",
                rank, Sr.bin, Sr.out_step, Sr.in_step, q, info, before, PfM_real[q]);
        SamplerRepairFlush();
        SamplerRepairFail("component factorization failed");
      }
      any = 1;
      Sr.components++;
      Sr.total++;
      /* The validated reference uses pre-repair kappa here. */
      Sr.reference[q] = Sr.k[q];
    }
  }
  Sr.refreshed += any;
  if (any) SamplerRepairKappa(PfM_real, Sr.k);
  for (q = 0; q < Sr.n; ++q)
    if (Sr.selected[q] || Sr.k[q] < Sr.floor[q]) Sr.floor[q] = Sr.k[q];
  return any;
}

static void SamplerRepairMeasure(double log_abs_ip, double projection,
                                 double stored_log_weight) {
  double d;
  if (!Sr.logging) return;
  d = fabs(2.0 * (log_abs_ip + projection) - stored_log_weight);
  Sr.samples++;
  Sr.measurement_zero += log_abs_ip == -INFINITY;
  Sr.measure_log = log_abs_ip;
  Sr.measure_projection = projection;
  Sr.measure_stored = stored_log_weight;
  Sr.measure_issue = !isfinite(log_abs_ip) || !isfinite(stored_log_weight);
  if (!isfinite(d)) Sr.nonfinite++;
  else {
    if (d > Sr.max_d) Sr.max_d = d;
    Sr.gt6 += d > 1e-6;
    Sr.gt3 += d > 1e-3;
    Sr.gt1 += d > 1.0;
  }
}

/* Record the actual existing measurement outcome, including zero contribution
 * and pre-existing numeric skips, rather than inferring recovery from a count. */
static void SamplerRepairMeasureFinish(int sample, double weight, double energy, int outcome) {
  if (!Sr.logging) return;
  Sr.measure_weight_zero += weight == 0.0;
  Sr.skip_weight += outcome == 1;
  Sr.skip_energy += outcome == 2;
  if (Sr.measure_issue || weight == 0.0 || outcome != 0) {
    Sr.out_step = Sr.in_step = -1; Sr.move = 0;
    SamplerRepairEvent("measure", -1, 0, Sr.measure_log, Sr.measure_projection,
                      Sr.measure_stored, weight, energy, (double)sample, outcome,
                      outcome == 1 ? "skip-weight" : outcome == 2 ? "skip-energy" :
                      weight == 0.0 ? "zero-contribution" : "accumulate");
  }
}

static void SamplerRepairMeasureFactorFailure(int sample, int info) {
  if (!Sr.logging) return;
  Sr.measure_factor_failed++;
  Sr.out_step = Sr.in_step = -1; Sr.move = 0;
  SamplerRepairEvent("measure", -1, 0, NAN, NAN, NAN, NAN, NAN,
                    (double)sample, info, "skip-factorization");
}

static void SamplerRepairFinalize(void) {
  SamplerRepairFlush();
  if (Sr.file) {
    int close_status;
    fprintf(Sr.file, "F refreshTotal %ld\n", Sr.total);
    close_status = fclose(Sr.file);
    Sr.file = NULL;
    if (close_status != 0) SamplerRepairFail("cannot close drift summary");
  }
  free(Sr.reference); free(Sr.floor); free(Sr.old_pf); free(Sr.new_pf);
  free(Sr.k); free(Sr.old_k); free(Sr.new_k); free(Sr.selected);
  free(Sr.buffer); free(Sr.work); free(Sr.iwork);
  memset(&Sr, 0, sizeof(Sr));
}
#endif
