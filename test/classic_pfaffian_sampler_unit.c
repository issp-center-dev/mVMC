#include "classic_pfaffian_sampler.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _mpi_use
#include <mpi.h>
#define TEST_COMMUNICATOR MPI_COMM_WORLD
#else
#define TEST_COMMUNICATOR 0
#endif

static int failure_count = 0;
static int world_rank = 0;
static int world_size = 1;

#define CHECK(condition, message)                                         \
  do {                                                                    \
    if (!(condition)) {                                                   \
      fprintf(stderr, "FAIL rank %d: %s (line %d)\n", world_rank,       \
              (message), __LINE__);                                       \
      ++failure_count;                                                    \
    }                                                                     \
  } while (0)

typedef struct {
  double value;
  int count;
} DrawContext;

typedef struct {
  MVMCClassicProposalDecision decisions[4];
  double uniforms[4];
  uint64_t generations[4];
  int draw_counts[4];
} ProposalTrace;

typedef struct {
  uint64_t proposal_id;
  MVMCClassicMoveKind move;
  int index_count;
  int indices[MVMC_CLASSIC_TRACE_INDEX_CAPACITY];
  int uniform_draw_count;
  double uniform;
  double complex current_total;
  double complex proposal_total;
  size_t regular_count;
  size_t near_singular_count;
  size_t singular_count;
  MVMCClassicProposalDecision decision;
  uint64_t accepted_generation;
  int ele_idx[2];
} P2DTrajectoryRecord;

typedef struct {
  P2DTrajectoryRecord records[2];
} P2DTrajectory;

typedef struct {
  int call_count;
  int qp_start;
  int inject_global_qp;
  int inject_inverse;
  int fail_observer;
} AuditObserverContext;

static double draw_once(void *context) {
  DrawContext *draw = (DrawContext *)context;
  ++draw->count;
  return draw->value;
}

static MVMCPfaffianStatus observe_real_candidate(
    void *context, const MVMCClassicPfaffianState *absolute_candidate,
    const double *absolute_inverse,
    MVMCClassicPfaffianRealObservation *observation) {
  AuditObserverContext *audit = (AuditObserverContext *)context;
  const size_t inverse_count = observation->component_count *
                               observation->matrix_element_count;
  int local_qp;

  ++audit->call_count;
  if (audit->fail_observer) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  if (observation->component_count != 0) {
    memcpy(observation->components, absolute_candidate->components,
           observation->component_count * sizeof(*observation->components));
  }
  if (inverse_count != 0) {
    memcpy(observation->inverse, absolute_inverse,
           inverse_count * sizeof(*observation->inverse));
  }
  local_qp = audit->inject_global_qp - audit->qp_start;
  if (local_qp >= 0 &&
      (size_t)local_qp < observation->component_count) {
    if (audit->inject_inverse) {
      observation->inverse[(size_t)local_qp *
                           observation->matrix_element_count] = 1.0;
    } else {
      if (observation->components[local_qp].state ==
          MVMC_PFAFFIAN_SINGULAR) {
        observation->components[local_qp].state = MVMC_PFAFFIAN_REGULAR;
        observation->components[local_qp].inverse_valid = 1;
        observation->components[local_qp].pfaffian = 1.0;
      } else {
        observation->components[local_qp].pfaffian += 1.0;
      }
    }
  }
  observation->valid = 1;
  return MVMC_PFAFFIAN_STATUS_OK;
}

static MVMCPfaffianStatus observe_complex_candidate(
    void *context, const MVMCClassicPfaffianState *absolute_candidate,
    const double complex *absolute_inverse,
    MVMCClassicPfaffianComplexObservation *observation) {
  AuditObserverContext *audit = (AuditObserverContext *)context;
  const size_t inverse_count = observation->component_count *
                               observation->matrix_element_count;

  ++audit->call_count;
  if (audit->fail_observer) return MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT;
  if (observation->component_count != 0) {
    memcpy(observation->components, absolute_candidate->components,
           observation->component_count * sizeof(*observation->components));
  }
  if (inverse_count != 0) {
    memcpy(observation->inverse, absolute_inverse,
           inverse_count * sizeof(*observation->inverse));
  }
  observation->valid = 1;
  return MVMC_PFAFFIAN_STATUS_OK;
}

static void save_p2d_trajectory(P2DTrajectory *trajectory, int index,
                                const MVMCClassicProposalTrace *trace,
                                const int *ele_idx) {
  P2DTrajectoryRecord *record = &trajectory->records[index];
  int item;

  memset(record, 0, sizeof(*record));
  record->proposal_id = trace->proposal_id;
  record->move = trace->move;
  record->index_count = trace->index_count;
  for (item = 0; item < trace->index_count; ++item) {
    record->indices[item] = trace->indices[item];
  }
  record->uniform_draw_count = trace->uniform_draw_count;
  record->uniform = trace->uniform;
  record->current_total = trace->current_total;
  record->proposal_total = trace->proposal_total;
  record->regular_count = trace->regular_count;
  record->near_singular_count = trace->near_singular_count;
  record->singular_count = trace->singular_count;
  record->decision = trace->decision;
  record->accepted_generation = trace->accepted_generation;
  record->ele_idx[0] = ele_idx[0];
  record->ele_idx[1] = ele_idx[1];
}

static void rank_partition(int qp_total, int *qp_start, int *qp_end) {
  *qp_start = (qp_total * world_rank) / world_size;
  *qp_end = (qp_total * (world_rank + 1)) / world_size;
}

static void set_real_pair(double *matrix, int row, int column, double value) {
  matrix[(size_t)row * 4 + (size_t)column] = value;
  matrix[(size_t)column * 4 + (size_t)row] = -value;
}

static void set_complex_pair(double complex *matrix, int row, int column,
                             double complex value) {
  matrix[(size_t)row * 4 + (size_t)column] = value;
  matrix[(size_t)column * 4 + (size_t)row] = -value;
}

static void fill_real_slater(double *slater, int qp_total) {
  int qp;
  memset(slater, 0, (size_t)qp_total * 16 * sizeof(*slater));
  for (qp = 0; qp < qp_total; ++qp) {
    double *matrix = slater + (size_t)qp * 16;
    const double base = (double)(qp + 1);
    set_real_pair(matrix, 0, 2, base);
    set_real_pair(matrix, 1, 2, 1.25 * base);
    set_real_pair(matrix, 0, 3, qp == 0 ? 0.0 : 1.5 * base);
    set_real_pair(matrix, 1, 3, 0.0);
  }
}

static void fill_complex_slater(double complex *slater, int qp_total) {
  int qp;
  memset(slater, 0, (size_t)qp_total * 16 * sizeof(*slater));
  for (qp = 0; qp < qp_total; ++qp) {
    double complex *matrix = slater + (size_t)qp * 16;
    const double complex base =
        (double)(qp + 1) * (1.0 + 0.2 * (double)qp * I);
    set_complex_pair(matrix, 0, 2, base);
    set_complex_pair(matrix, 1, 2, 0.5 * base);
    set_complex_pair(matrix, 0, 3, qp == 0 ? 0.0 : 1.25 * base);
    set_complex_pair(matrix, 1, 3, 0.0);
  }
}

static void save_trace(ProposalTrace *trace, int index,
                       const MVMCClassicProposalResult *result) {
  trace->decisions[index] = result->decision;
  trace->uniforms[index] = result->uniform;
  trace->generations[index] = result->accepted_generation;
  trace->draw_counts[index] = result->uniform_draw_count;
}

static void run_real_sequence(int qp_total, ProposalTrace *trace) {
  const int initial_idx[2] = {0, 0};
  const int hopping_idx[2] = {1, 0};
  const int spin_hopping_idx[2] = {1, 1};
  const int exchange_idx[2] = {0, 1};
  double complex weights[4] = {1.0, 0.75, -0.25, 0.5};
  double slater[64];
  double legacy[20];
  MVMCClassicPfaffianRealWorkspace *state_workspace = NULL;
  MVMCClassicPfaffianCollectiveWorkspace *collective_workspace = NULL;
  MVMCClassicSamplerState sampler_state;
  MVMCClassicProposalResult result;
  DrawContext draw;
  int qp_start, qp_end;
  int index;

  memset(trace, 0, sizeof(*trace));
  for (index = 0; index < 20; ++index) legacy[index] = 71.0;
  rank_partition(qp_total, &qp_start, &qp_end);
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 2, qp_total, qp_start, qp_end, &state_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create real sampler state workspace");
  CHECK(mvmc_classic_pfaffian_collective_workspace_create(
            TEST_COMMUNICATOR, &collective_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create real sampler collective workspace");
  fill_real_slater(slater, qp_total);
  CHECK(mvmc_classic_sampler_real_initialize(
            state_workspace, collective_workspace, slater, initial_idx, weights,
            0.0, 0.0, legacy, (size_t)qp_total * 5, &sampler_state) ==
            MVMC_PFAFFIAN_STATUS_OK && sampler_state.valid,
        "initialize real sampler from absolute state");

  draw.value = 0.75;
  draw.count = 0;
  CHECK(mvmc_classic_sampler_real_propose(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_PAIR_HOPPING, slater, spin_hopping_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, legacy,
            (size_t)qp_total * 5, &result) == MVMC_PFAFFIAN_STATUS_OK &&
            result.decision == MVMC_CLASSIC_PROPOSAL_REJECTED &&
            result.proposal_total == 0.0,
        "rank-two pair hopping from a doubly occupied site rejects zero total");
  CHECK(draw.count == 1 && result.uniform_draw_count == 1,
        "exact-zero pair hopping still consumes exactly one uniform");
  save_trace(trace, 0, &result);

  draw.value = 0.25;
  draw.count = 0;
  CHECK(mvmc_classic_sampler_real_propose(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_HOPPING, slater, hopping_idx, weights, 0.0, 0.0,
            0.0, 0.0, draw_once, &draw, legacy, (size_t)qp_total * 5,
            &result) == MVMC_PFAFFIAN_STATUS_OK && result.valid &&
            result.decision == MVMC_CLASSIC_PROPOSAL_ACCEPTED,
        "rank-one hopping accepts authoritative absolute candidate");
  CHECK(draw.count == 1 && result.uniform_draw_count == 1,
        "accepted hopping consumes exactly one uniform");
  save_trace(trace, 1, &result);

  draw.value = 0.5;
  draw.count = 0;
  CHECK(mvmc_classic_sampler_real_propose(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_SPIN_HOPPING, slater, spin_hopping_idx, weights,
            0.0, 0.0, -INFINITY, 0.0, draw_once, &draw, legacy,
            (size_t)qp_total * 5, &result) == MVMC_PFAFFIAN_STATUS_OK &&
            result.decision == MVMC_CLASSIC_PROPOSAL_REJECTED,
        "spin hopping rejection preserves accepted state");
  CHECK(result.accepted_generation == 2 && draw.count == 1,
        "rejected rank-one move preserves generation and one draw");
  save_trace(trace, 2, &result);

  draw.value = 0.0;
  draw.count = 0;
  CHECK(mvmc_classic_sampler_real_propose(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_EXCHANGE, slater, exchange_idx, weights, 0.0, 0.0,
            0.0, 0.0, draw_once, &draw, legacy, (size_t)qp_total * 5,
            &result) == MVMC_PFAFFIAN_STATUS_OK,
        "rank-two exchange evaluates absolute candidate");
  if (qp_total == 1) {
    CHECK(result.decision == MVMC_CLASSIC_PROPOSAL_REJECTED,
          "single-QP singular exchange has exact-zero total");
  } else {
    CHECK(result.decision == MVMC_CLASSIC_PROPOSAL_ACCEPTED &&
              sampler_state.accepted_aggregate.singular_count == 1,
          "multi-QP singular component with nonzero total is accepted");
  }
  CHECK(draw.count == 1,
        "exchange consumes one uniform even for exact-zero proposal");
  save_trace(trace, 3, &result);

  draw.count = 0;
  CHECK(mvmc_classic_sampler_real_propose(
            state_workspace, collective_workspace, &sampler_state,
            (MVMCClassicMoveKind)99, slater, initial_idx, weights, 0.0, 0.0,
            0.0, 0.0, draw_once, &draw, legacy, (size_t)qp_total * 5,
            &result) == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            draw.count == 0,
        "unsupported move fails before numerical evaluation and RNG draw");

  slater[6] = NAN;
  slater[9] = NAN;
  draw.count = 0;
  CHECK(mvmc_classic_sampler_real_propose(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_HOPPING, slater, hopping_idx, weights, 0.0, 0.0,
            0.0, 0.0, draw_once, &draw, legacy, (size_t)qp_total * 5,
            &result) == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            draw.count == 0,
        "invalid candidate fails collectively before RNG draw");

  mvmc_classic_pfaffian_collective_workspace_destroy(collective_workspace);
  mvmc_classic_pfaffian_real_workspace_destroy(state_workspace);
}

static void run_complex_rbm_sequence(int qp_total) {
  const int initial_idx[2] = {0, 0};
  const int hopping_idx[2] = {1, 0};
  const int pair_idx[2] = {1, 1};
  const int exchange_idx[2] = {0, 1};
  const double complex weights[4] = {
      1.0, 0.5 + 0.25 * I, -0.75 + 0.5 * I, 1.25 - 0.5 * I};
  double complex slater[64];
  double complex legacy[20];
  MVMCClassicPfaffianComplexWorkspace *state_workspace = NULL;
  MVMCClassicPfaffianCollectiveWorkspace *collective_workspace = NULL;
  MVMCClassicSamplerState sampler_state;
  MVMCClassicProposalResult result;
  DrawContext draw;
  int qp_start, qp_end;

  memset(legacy, 0, sizeof(legacy));
  rank_partition(qp_total, &qp_start, &qp_end);
  CHECK(mvmc_classic_pfaffian_complex_workspace_create(
            2, 2, qp_total, qp_start, qp_end, &state_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create complex sampler state workspace");
  CHECK(mvmc_classic_pfaffian_collective_workspace_create(
            TEST_COMMUNICATOR, &collective_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create complex sampler collective workspace");
  fill_complex_slater(slater, qp_total);
  CHECK(mvmc_classic_sampler_complex_initialize(
            state_workspace, collective_workspace, slater, initial_idx, weights,
            0.0, 0.0, legacy, (size_t)qp_total * 5, &sampler_state) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "initialize complex sampler");

  draw.value = 0.75;
  draw.count = 0;
  CHECK(mvmc_classic_sampler_complex_propose(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_PAIR_HOPPING, slater, pair_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, legacy,
            (size_t)qp_total * 5, &result) == MVMC_PFAFFIAN_STATUS_OK &&
            result.decision == MVMC_CLASSIC_PROPOSAL_REJECTED &&
            result.proposal_total == 0.0 && draw.count == 1,
        "complex pair hopping mutates both electron indices and rejects zero");

  draw.value = 0.5;
  draw.count = 0;
  CHECK(mvmc_classic_sampler_complex_propose(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_HOPPING, slater, hopping_idx, weights, 0.0, 0.0,
            0.0, 0.0, draw_once, &draw, legacy, (size_t)qp_total * 5,
            &result) == MVMC_PFAFFIAN_STATUS_OK &&
            result.decision == MVMC_CLASSIC_PROPOSAL_REJECTED,
        "complex proposal without RBM factor rejects");
  CHECK(draw.count == 1, "complex rejection consumes one uniform");

  draw.count = 0;
  CHECK(mvmc_classic_sampler_complex_propose(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_HOPPING, slater, hopping_idx, weights, 0.0, 0.0,
            0.0, log(3.0), draw_once, &draw, legacy,
            (size_t)qp_total * 5, &result) == MVMC_PFAFFIAN_STATUS_OK &&
            result.decision == MVMC_CLASSIC_PROPOSAL_ACCEPTED,
        "complex RBM log ratio participates in total-amplitude acceptance");
  CHECK(draw.count == 1 && result.accepted_generation == 2,
        "RBM smoke publishes the absolute candidate once");

  draw.value = 0.0;
  draw.count = 0;
  CHECK(mvmc_classic_sampler_complex_propose(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_SPIN_HOPPING, slater, pair_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, legacy,
            (size_t)qp_total * 5, &result) == MVMC_PFAFFIAN_STATUS_OK &&
            result.decision == MVMC_CLASSIC_PROPOSAL_REJECTED &&
            draw.count == 1,
        "complex spin hopping uses the same absolute transaction");

  draw.count = 0;
  CHECK(mvmc_classic_sampler_complex_propose(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_EXCHANGE, slater, exchange_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, legacy,
            (size_t)qp_total * 5, &result) == MVMC_PFAFFIAN_STATUS_OK &&
            (qp_total == 1
                 ? result.decision == MVMC_CLASSIC_PROPOSAL_REJECTED
                 : result.decision == MVMC_CLASSIC_PROPOSAL_ACCEPTED) &&
            draw.count == 1,
        "complex exchange handles exact-zero and partial singularity");

  mvmc_classic_pfaffian_collective_workspace_destroy(collective_workspace);
  mvmc_classic_pfaffian_complex_workspace_destroy(state_workspace);
}

static void run_p2d_real_trajectory(MVMCClassicPfaffianAuditMode mode,
                                    P2DTrajectory *trajectory) {
  const int initial_idx[2] = {0, 0};
  const int singular_idx[2] = {0, 1};
  const double complex weights[4] = {1.0, 0.75, -0.25, 0.5};
  double slater[64];
  double legacy[20] = {0};
  MVMCClassicPfaffianRealWorkspace *state_workspace = NULL;
  MVMCClassicPfaffianCollectiveWorkspace *collective_workspace = NULL;
  MVMCClassicPfaffianRealAuditWorkspace *audit_workspace = NULL;
  MVMCClassicSamplerState sampler_state;
  MVMCClassicProposalResult result;
  MVMCClassicProposalTrace trace;
  MVMCClassicPfaffianAuditReport rebuild_report;
  MVMCClassicProposalMetadata metadata;
  MVMCClassicRealProposalAudit audit;
  AuditObserverContext observer = {0, 0, -1, 0, 0};
  DrawContext draw = {0.0, 0};
  int qp_start, qp_end;

  memset(trajectory, 0, sizeof(*trajectory));
  rank_partition(4, &qp_start, &qp_end);
  observer.qp_start = qp_start;
  fill_real_slater(slater, 4);
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 2, 4, qp_start, qp_end, &state_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create P2-D real trajectory state workspace");
  CHECK(mvmc_classic_pfaffian_collective_workspace_create(
            TEST_COMMUNICATOR, &collective_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create P2-D real trajectory collective workspace");
  CHECK(mvmc_classic_pfaffian_real_audit_workspace_create(
            state_workspace, mode, 1.0e-12, 1.0e-11,
            &audit_workspace) == MVMC_PFAFFIAN_STATUS_OK,
        "create P2-D real trajectory audit workspace");
  CHECK(mvmc_classic_sampler_real_initialize(
            state_workspace, collective_workspace, slater, initial_idx,
            weights, 0.0, 0.0, legacy, 20, &sampler_state) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "initialize P2-D real trajectory");

  memset(&metadata, 0, sizeof(metadata));
  metadata.proposal_id = 100;
  metadata.index_count = 2;
  metadata.indices[0] = 0;
  metadata.indices[1] = 1;
  audit.metadata = &metadata;
  audit.workspace = audit_workspace;
  audit.observer = observe_real_candidate;
  audit.observer_context = &observer;
  audit.trace = &trace;
  CHECK(mvmc_classic_sampler_real_propose_with_audit(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_EXCHANGE, slater, singular_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, &audit, legacy, 20,
            &result) == MVMC_PFAFFIAN_STATUS_OK && result.valid &&
            result.decision == MVMC_CLASSIC_PROPOSAL_ACCEPTED && trace.valid &&
            sampler_state.accepted_aggregate.singular_count == 1,
        "deterministic R-to-S proposal is accepted and traced");
  CHECK(trace.proposal_id == 100 && trace.index_count == 2 &&
            trace.indices[0] == 0 && trace.indices[1] == 1,
        "deterministic trace records proposal ID and indices");
  save_p2d_trajectory(trajectory, 0, &trace, singular_idx);
  CHECK(mvmc_classic_sampler_real_rebuild_audit(
            state_workspace, collective_workspace, &sampler_state, slater,
            singular_idx, weights, 0.0, 0.0, audit_workspace,
            &rebuild_report) == MVMC_PFAFFIAN_STATUS_OK &&
            rebuild_report.valid && rebuild_report.executed &&
            !rebuild_report.mismatch,
        "periodic rebuild matches accepted singular state without healing");

  metadata.proposal_id = 101;
  metadata.indices[0] = 1;
  metadata.indices[1] = 0;
  draw.count = 0;
  CHECK(mvmc_classic_sampler_real_propose_with_audit(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_EXCHANGE, slater, initial_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, &audit, legacy, 20,
            &result) == MVMC_PFAFFIAN_STATUS_OK && result.valid &&
            result.decision == MVMC_CLASSIC_PROPOSAL_ACCEPTED && trace.valid &&
            sampler_state.accepted_aggregate.singular_count == 0,
        "deterministic S-to-R proposal is accepted and traced");
  save_p2d_trajectory(trajectory, 1, &trace, initial_idx);
  CHECK(mvmc_classic_sampler_real_rebuild_audit(
            state_workspace, collective_workspace, &sampler_state, slater,
            initial_idx, weights, 0.0, 0.0, audit_workspace,
            &rebuild_report) == MVMC_PFAFFIAN_STATUS_OK &&
            rebuild_report.valid && !rebuild_report.mismatch,
        "periodic rebuild matches recovered regular state");
  CHECK(observer.call_count ==
            (mode == MVMC_CLASSIC_AUDIT_OFF ? 0 : 2),
        "audit off performs zero callbacks and always mode one per proposal");

  mvmc_classic_pfaffian_real_audit_workspace_destroy(audit_workspace);
  mvmc_classic_pfaffian_collective_workspace_destroy(collective_workspace);
  mvmc_classic_pfaffian_real_workspace_destroy(state_workspace);
}

static void test_p2d_audit_trajectory_identity(void) {
  P2DTrajectory audit_off;
  P2DTrajectory audit_on;

  run_p2d_real_trajectory(MVMC_CLASSIC_AUDIT_OFF, &audit_off);
  run_p2d_real_trajectory(MVMC_CLASSIC_AUDIT_ALWAYS, &audit_on);
  CHECK(memcmp(&audit_off, &audit_on, sizeof(audit_off)) == 0,
        "audit off/on trajectory is byte-identical outside audit metadata");
}

static void test_p2d_guarded_fallback(void) {
  const int initial_idx[2] = {0, 0};
  const int singular_idx[2] = {0, 1};
  const double complex weights[4] = {1.0, 0.75, -0.25, 0.5};
  double slater[64];
  double legacy[20] = {0};
  MVMCClassicPfaffianRealWorkspace *state_workspace = NULL;
  MVMCClassicPfaffianCollectiveWorkspace *collective_workspace = NULL;
  MVMCClassicPfaffianRealAuditWorkspace *audit_workspace = NULL;
  MVMCClassicSamplerState sampler_state;
  MVMCClassicProposalResult result;
  MVMCClassicProposalMetadata metadata = {200, 0, {0, 0, 0, 0}};
  MVMCClassicRealProposalAudit audit;
  AuditObserverContext observer = {0, 0, -1, 0, 0};
  DrawContext draw = {0.0, 0};
  int qp_start, qp_end;

  rank_partition(4, &qp_start, &qp_end);
  observer.qp_start = qp_start;
  fill_real_slater(slater, 4);
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 2, 4, qp_start, qp_end, &state_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create guarded audit state workspace");
  CHECK(mvmc_classic_pfaffian_collective_workspace_create(
            TEST_COMMUNICATOR, &collective_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create guarded audit collective workspace");
  CHECK(mvmc_classic_pfaffian_real_audit_workspace_create(
            state_workspace, MVMC_CLASSIC_AUDIT_GUARDED,
            1.0e-12, 1.0e-11, &audit_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create guarded audit workspace");
  CHECK(mvmc_classic_sampler_real_initialize(
            state_workspace, collective_workspace, slater, initial_idx,
            weights, 0.0, 0.0, legacy, 20, &sampler_state) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "initialize guarded audit trajectory");
  audit.metadata = &metadata;
  audit.workspace = audit_workspace;
  audit.observer = observe_real_candidate;
  audit.observer_context = &observer;
  audit.trace = NULL;
  CHECK(mvmc_classic_sampler_real_propose_with_audit(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_EXCHANGE, slater, singular_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, &audit, legacy, 20,
            &result) == MVMC_PFAFFIAN_STATUS_OK &&
            result.audit.executed && !result.audit.fallback &&
            observer.call_count == 1,
        "guarded audit executes from an all-regular accepted state");
  metadata.proposal_id = 201;
  draw.count = 0;
  CHECK(mvmc_classic_sampler_real_propose_with_audit(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_EXCHANGE, slater, initial_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, &audit, legacy, 20,
            &result) == MVMC_PFAFFIAN_STATUS_OK &&
            result.audit.fallback && !result.audit.executed &&
            observer.call_count == 1 && draw.count == 1,
        "one singular accepted component selects all-rank absolute fallback");

  mvmc_classic_pfaffian_real_audit_workspace_destroy(audit_workspace);
  mvmc_classic_pfaffian_collective_workspace_destroy(collective_workspace);
  mvmc_classic_pfaffian_real_workspace_destroy(state_workspace);
}

static void test_p2d_audit_mismatch_and_rebuild_no_heal(void) {
  const int initial_idx[2] = {0, 0};
  const int singular_idx[2] = {0, 1};
  const int wrong_rebuild_idx[2] = {1, 0};
  const double complex weights[4] = {1.0, 0.75, -0.25, 0.5};
  double slater[64];
  double legacy[20] = {0};
  double legacy_before[20];
  double inverse_before[16];
  MVMCAbsolutePfaffianResult components_before[4];
  MVMCClassicPfaffianRealWorkspace *state_workspace = NULL;
  MVMCClassicPfaffianCollectiveWorkspace *collective_workspace = NULL;
  MVMCClassicPfaffianRealAuditWorkspace *audit_workspace = NULL;
  MVMCClassicSamplerState sampler_state;
  MVMCClassicProposalResult result;
  MVMCClassicPfaffianAuditReport rebuild_report;
  MVMCClassicPfaffianAuditReport aggregate_report;
  MVMCProjectedAmplitudeResult wrong_aggregate;
  MVMCClassicProposalMetadata metadata = {300, 2, {0, 1, 0, 0}};
  MVMCClassicRealProposalAudit audit;
  AuditObserverContext observer = {0, 0, 0, 0, 0};
  DrawContext draw = {0.0, 0};
  const MVMCClassicPfaffianState *accepted;
  const double *accepted_inverse;
  size_t local_count;
  int qp_start, qp_end;

  rank_partition(4, &qp_start, &qp_end);
  local_count = (size_t)(qp_end - qp_start);
  observer.qp_start = qp_start;
  fill_real_slater(slater, 4);
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 2, 4, qp_start, qp_end, &state_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create mismatch audit state workspace");
  CHECK(mvmc_classic_pfaffian_collective_workspace_create(
            TEST_COMMUNICATOR, &collective_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create mismatch audit collective workspace");
  CHECK(mvmc_classic_pfaffian_real_audit_workspace_create(
            state_workspace, MVMC_CLASSIC_AUDIT_ALWAYS,
            1.0e-12, 1.0e-11, &audit_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create mismatch audit workspace");
  CHECK(mvmc_classic_sampler_real_initialize(
            state_workspace, collective_workspace, slater, initial_idx,
            weights, 0.0, 0.0, legacy, 20, &sampler_state) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "initialize mismatch audit state");
  audit.metadata = &metadata;
  audit.workspace = audit_workspace;
  audit.observer = observe_real_candidate;
  audit.observer_context = &observer;
  audit.trace = NULL;
  observer.inject_global_qp = 1;
  CHECK(mvmc_classic_sampler_real_propose_with_audit(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_EXCHANGE, slater, singular_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, &audit, legacy, 20,
            &result) == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            result.audit.mismatch && result.audit.failure_global_qp == 1 &&
            draw.count == 0,
        "one-rank regular Pfaffian mismatch is collective before RNG");
  observer.inject_global_qp = 0;
  CHECK(mvmc_classic_sampler_real_propose_with_audit(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_EXCHANGE, slater, singular_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, &audit, legacy, 20,
            &result) == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            result.audit.mismatch && result.audit.failure_global_qp == 0 &&
            draw.count == 0 &&
            mvmc_classic_pfaffian_real_accepted(state_workspace)->generation ==
                1,
        "one-rank singular Pfaffian mismatch is collective before RNG");
  observer.inject_inverse = 1;
  draw.count = 0;
  CHECK(mvmc_classic_sampler_real_propose_with_audit(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_EXCHANGE, slater, singular_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, &audit, legacy, 20,
            &result) == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            result.audit.mismatch && result.audit.failure_global_qp == 0 &&
            draw.count == 0,
        "nonregular nonzero inverse observation is never skipped");

  observer.inject_global_qp = -1;
  observer.inject_inverse = 0;
  observer.fail_observer = world_rank == world_size - 1;
  draw.count = 0;
  CHECK(mvmc_classic_sampler_real_propose_with_audit(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_EXCHANGE, slater, singular_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, &audit, legacy, 20,
            &result) == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            result.audit.mismatch &&
            result.audit.failure_rank == world_size - 1 && draw.count == 0,
        "one-rank observer failure is collective before RNG");
  observer.fail_observer = 0;

  accepted = mvmc_classic_pfaffian_real_accepted(state_workspace);
  accepted_inverse = mvmc_classic_pfaffian_real_accepted_inverse(
      state_workspace);
  wrong_aggregate = sampler_state.accepted_aggregate;
  wrong_aggregate.total += 1.0;
  CHECK(mvmc_classic_pfaffian_real_audit_compare(
            audit_workspace, collective_workspace, accepted,
            accepted_inverse, &wrong_aggregate, accepted,
            accepted_inverse, &aggregate_report) ==
                MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            aggregate_report.mismatch &&
            aggregate_report.failure_global_qp == -1,
        "global aggregate mismatch is detected independently of components");
  memcpy(legacy_before, legacy, sizeof(legacy));
  if (local_count != 0) {
    memcpy(components_before, accepted->components,
           local_count * sizeof(*components_before));
    memcpy(inverse_before, accepted_inverse,
           local_count * 4 * sizeof(*inverse_before));
  }
  CHECK(mvmc_classic_sampler_real_rebuild_audit(
            state_workspace, collective_workspace, &sampler_state, slater,
            wrong_rebuild_idx, weights, 0.0, 0.0, audit_workspace,
            &rebuild_report) == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            rebuild_report.mismatch,
        "periodic rebuild mismatch is an error instead of silent heal");
  accepted = mvmc_classic_pfaffian_real_accepted(state_workspace);
  accepted_inverse = mvmc_classic_pfaffian_real_accepted_inverse(
      state_workspace);
  CHECK(accepted->generation == 1 &&
            memcmp(legacy_before, legacy, sizeof(legacy)) == 0 &&
            (local_count == 0 ||
             (memcmp(components_before, accepted->components,
                     local_count * sizeof(*components_before)) == 0 &&
              memcmp(inverse_before, accepted_inverse,
                     local_count * 4 * sizeof(*inverse_before)) == 0)),
        "rebuild mismatch preserves accepted bytes, mirror, and generation");

  metadata.index_count = MVMC_CLASSIC_TRACE_INDEX_CAPACITY + 1;
  draw.count = 0;
  CHECK(mvmc_classic_sampler_real_propose_with_audit(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_HOPPING, slater, wrong_rebuild_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, &audit, legacy, 20,
            &result) == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            draw.count == 0 && accepted->generation == 1,
        "invalid trace metadata fails collectively before prepare or RNG");

  mvmc_classic_pfaffian_real_audit_workspace_destroy(audit_workspace);
  mvmc_classic_pfaffian_collective_workspace_destroy(collective_workspace);
  mvmc_classic_pfaffian_real_workspace_destroy(state_workspace);
}

static void test_p2d_complex_audit_smoke(void) {
  const int initial_idx[2] = {0, 0};
  const int exchange_idx[2] = {0, 1};
  const double complex weights[4] = {
      1.0, 0.5 + 0.25 * I, -0.75 + 0.5 * I, 1.25 - 0.5 * I};
  double complex slater[64];
  double complex legacy[20] = {0};
  MVMCClassicPfaffianComplexWorkspace *state_workspace = NULL;
  MVMCClassicPfaffianCollectiveWorkspace *collective_workspace = NULL;
  MVMCClassicPfaffianComplexAuditWorkspace *audit_workspace = NULL;
  MVMCClassicSamplerState sampler_state;
  MVMCClassicProposalResult result;
  MVMCClassicPfaffianAuditReport rebuild_report;
  MVMCClassicProposalMetadata metadata = {400, 2, {0, 1, 0, 0}};
  MVMCClassicComplexProposalAudit audit;
  AuditObserverContext observer = {0, 0, -1, 0, 0};
  DrawContext draw = {0.0, 0};
  int qp_start, qp_end;

  rank_partition(1, &qp_start, &qp_end);
  observer.qp_start = qp_start;
  fill_complex_slater(slater, 1);
  CHECK(mvmc_classic_pfaffian_complex_workspace_create(
            2, 2, 1, qp_start, qp_end, &state_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create complex P2-D state workspace");
  CHECK(mvmc_classic_pfaffian_collective_workspace_create(
            TEST_COMMUNICATOR, &collective_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create complex P2-D collective workspace");
  CHECK(mvmc_classic_pfaffian_complex_audit_workspace_create(
            state_workspace, MVMC_CLASSIC_AUDIT_ALWAYS,
            1.0e-12, 1.0e-11, &audit_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create complex P2-D audit workspace");
  CHECK(mvmc_classic_sampler_complex_initialize(
            state_workspace, collective_workspace, slater, initial_idx,
            weights, 0.0, 0.0, legacy, 20, &sampler_state) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "initialize complex P2-D audit state");
  audit.metadata = &metadata;
  audit.workspace = audit_workspace;
  audit.observer = observe_complex_candidate;
  audit.observer_context = &observer;
  audit.trace = NULL;
  CHECK(mvmc_classic_sampler_complex_propose_with_audit(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_EXCHANGE, slater, exchange_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, &audit, legacy, 20,
            &result) == MVMC_PFAFFIAN_STATUS_OK && result.audit.executed &&
            !result.audit.mismatch && observer.call_count == 1,
        "complex always-audit compares every component before acceptance");
  CHECK(mvmc_classic_sampler_complex_rebuild_audit(
            state_workspace, collective_workspace, &sampler_state, slater,
            initial_idx, weights, 0.0, 0.0, audit_workspace,
            &rebuild_report) == MVMC_PFAFFIAN_STATUS_OK &&
            !rebuild_report.mismatch,
        "complex periodic rebuild matches accepted state with empty slices");

  mvmc_classic_pfaffian_complex_audit_workspace_destroy(audit_workspace);
  mvmc_classic_pfaffian_collective_workspace_destroy(collective_workspace);
  mvmc_classic_pfaffian_complex_workspace_destroy(state_workspace);
}

static void test_p2d_rank_divergent_audit_mode(void) {
  const int initial_idx[2] = {0, 0};
  const int hopping_idx[2] = {1, 0};
  const double complex weights[4] = {1.0, 0.75, -0.25, 0.5};
  double slater[64];
  double legacy[20] = {0};
  MVMCClassicPfaffianRealWorkspace *state_workspace = NULL;
  MVMCClassicPfaffianCollectiveWorkspace *collective_workspace = NULL;
  MVMCClassicPfaffianRealAuditWorkspace *audit_workspace = NULL;
  MVMCClassicSamplerState sampler_state;
  MVMCClassicProposalResult result;
  MVMCClassicProposalMetadata metadata = {500, 0, {0, 0, 0, 0}};
  MVMCClassicRealProposalAudit audit;
  AuditObserverContext observer = {0, 0, -1, 0, 0};
  DrawContext draw = {0.0, 0};
  MVMCClassicPfaffianAuditMode local_mode;
  int qp_start, qp_end;

  if (world_size == 1) return;
  rank_partition(4, &qp_start, &qp_end);
  observer.qp_start = qp_start;
  local_mode = world_rank == world_size - 1
                   ? MVMC_CLASSIC_AUDIT_ALWAYS
                   : MVMC_CLASSIC_AUDIT_OFF;
  fill_real_slater(slater, 4);
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 2, 4, qp_start, qp_end, &state_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create rank-divergent audit state workspace");
  CHECK(mvmc_classic_pfaffian_collective_workspace_create(
            TEST_COMMUNICATOR, &collective_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create rank-divergent audit collective workspace");
  CHECK(mvmc_classic_pfaffian_real_audit_workspace_create(
            state_workspace, local_mode, 1.0e-12, 1.0e-11,
            &audit_workspace) == MVMC_PFAFFIAN_STATUS_OK,
        "create rank-local audit mode workspace");
  CHECK(mvmc_classic_sampler_real_initialize(
            state_workspace, collective_workspace, slater, initial_idx,
            weights, 0.0, 0.0, legacy, 20, &sampler_state) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "initialize rank-divergent audit mode fixture");
  audit.metadata = &metadata;
  audit.workspace = audit_workspace;
  audit.observer = observe_real_candidate;
  audit.observer_context = &observer;
  audit.trace = NULL;
  CHECK(mvmc_classic_sampler_real_propose_with_audit(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_HOPPING, slater, hopping_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, &audit, legacy, 20,
            &result) == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            draw.count == 0 && observer.call_count == 0 &&
            mvmc_classic_pfaffian_real_accepted(state_workspace)->generation ==
                1,
        "rank-divergent audit mode fails before prepare, callback, or RNG");

  mvmc_classic_pfaffian_real_audit_workspace_destroy(audit_workspace);
  mvmc_classic_pfaffian_collective_workspace_destroy(collective_workspace);
  mvmc_classic_pfaffian_real_workspace_destroy(state_workspace);
}

static void test_collective_preflight_mismatch(void) {
  MVMCClassicPfaffianCollectiveWorkspace *workspace = NULL;
  MVMCPfaffianStatus status;
  int all_true = -1;

  CHECK(mvmc_classic_pfaffian_collective_workspace_create(
            TEST_COMMUNICATOR, &workspace) == MVMC_PFAFFIAN_STATUS_OK,
        "create preflight mismatch workspace");
  status = mvmc_classic_pfaffian_collective_preflight(
      workspace, 1,
      world_size > 1 && world_rank == world_size - 1
          ? MVMC_CLASSIC_MOVE_EXCHANGE
          : MVMC_CLASSIC_MOVE_HOPPING);
  CHECK(status == (world_size > 1 ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT
                                  : MVMC_PFAFFIAN_STATUS_OK),
        "rank-divergent move kind is rejected before proposal collectives");
  status = mvmc_classic_pfaffian_collective_preflight(
      workspace, world_size == 1 || world_rank != world_size - 1,
      MVMC_CLASSIC_MOVE_HOPPING);
  CHECK(status == (world_size > 1 ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT
                                  : MVMC_PFAFFIAN_STATUS_OK),
        "one-rank invalid precondition is propagated before proposal");
  status = mvmc_classic_pfaffian_collective_all_true(
      workspace, world_size == 1 || world_rank != world_size - 1, &all_true);
  CHECK(status == MVMC_PFAFFIAN_STATUS_OK &&
            all_true == (world_size == 1 ? 1 : 0),
        "optional guarded branch uses communicator-wide logical AND");
  status = mvmc_classic_pfaffian_collective_all_true(
      workspace, world_size > 1 && world_rank == world_size - 1 ? 2 : 1,
      &all_true);
  CHECK(status == (world_size > 1 ? MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT
                                  : MVMC_PFAFFIAN_STATUS_OK),
        "malformed all-true input is rejected after collective participation");
  mvmc_classic_pfaffian_collective_workspace_destroy(workspace);
}

static void test_log_domain_metropolis_edges(void) {
  MVMCClassicPfaffianCollectiveWorkspace *workspace = NULL;
  MVMCClassicPfaffianMetropolisResult result;

  CHECK(mvmc_classic_pfaffian_collective_workspace_create(
            TEST_COMMUNICATOR, &workspace) == MVMC_PFAFFIAN_STATUS_OK,
        "create log-domain Metropolis workspace");
  CHECK(mvmc_classic_pfaffian_collective_metropolis(
            workspace, 1.0, 1.0, INFINITY, 0.999, &result) ==
            MVMC_PFAFFIAN_STATUS_OK && result.valid && result.accepted &&
            isinf(result.log_acceptance_ratio) &&
            result.log_acceptance_ratio > 0.0,
        "positive log-ratio overflow accepts without exp");
  CHECK(mvmc_classic_pfaffian_collective_metropolis(
            workspace, 1.0, 1.0, -INFINITY, 0.5, &result) ==
            MVMC_PFAFFIAN_STATUS_OK && !result.accepted,
        "negative log-ratio overflow rejects without exp");
  CHECK(mvmc_classic_pfaffian_collective_metropolis(
            workspace, 1.0, 1.0, -INFINITY, 0.0, &result) ==
            MVMC_PFAFFIAN_STATUS_OK && !result.accepted,
        "u equals zero cannot accept a zero-probability log ratio");
  CHECK(mvmc_classic_pfaffian_collective_metropolis(
            workspace, 1.0, 1.0, 0.0, 0.0, &result) ==
            MVMC_PFAFFIAN_STATUS_OK && result.accepted,
        "u equals zero accepts a positive-probability nonzero proposal");
  CHECK(mvmc_classic_pfaffian_collective_metropolis(
            workspace, 1.0, 0.0, INFINITY, 0.0, &result) ==
            MVMC_PFAFFIAN_STATUS_OK && result.valid && !result.accepted,
        "exact-zero proposal rejects independently of draw and log factor");
  CHECK(mvmc_classic_pfaffian_collective_metropolis(
            workspace, 0.0, 1.0, 0.0, 0.5, &result) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "exact-zero current total is invalid");
  CHECK(mvmc_classic_pfaffian_collective_metropolis(
            workspace, 1.0, 1.0, 0.0, 1.0, &result) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT,
        "uniform outside [0,1) is invalid");
  mvmc_classic_pfaffian_collective_workspace_destroy(workspace);
}

static void test_initialize_propose_preflight_mismatch(void) {
  const int initial_idx[2] = {0, 0};
  const int hopping_idx[2] = {1, 0};
  const double complex weights[2] = {1.0, 1.0};
  double slater[32];
  double legacy[10] = {0};
  MVMCClassicPfaffianRealWorkspace *state_workspace = NULL;
  MVMCClassicPfaffianCollectiveWorkspace *collective_workspace = NULL;
  MVMCClassicSamplerState sampler_state;
  MVMCClassicProposalResult result;
  DrawContext draw = {0.25, 0};
  MVMCPfaffianStatus status;
  int qp_start, qp_end;

  if (world_size == 1) return;
  rank_partition(2, &qp_start, &qp_end);
  fill_real_slater(slater, 2);
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 2, 2, qp_start, qp_end, &state_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create initialize/propose mismatch state workspace");
  CHECK(mvmc_classic_pfaffian_collective_workspace_create(
            TEST_COMMUNICATOR, &collective_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create initialize/propose mismatch collective workspace");
  CHECK(mvmc_classic_sampler_real_initialize(
            state_workspace, collective_workspace, slater, initial_idx,
            weights, 0.0, 0.0, legacy, 10, &sampler_state) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "initialize state before branch namespace mismatch");

  if (world_rank == world_size - 1) {
    status = mvmc_classic_sampler_real_initialize(
        state_workspace, collective_workspace, slater, initial_idx, weights,
        0.0, 0.0, legacy, 10, &sampler_state);
  } else {
    status = mvmc_classic_sampler_real_propose(
        state_workspace, collective_workspace, &sampler_state,
        MVMC_CLASSIC_MOVE_HOPPING, slater, hopping_idx, weights,
        0.0, 0.0, 0.0, 0.0, draw_once, &draw, legacy, 10, &result);
  }
  CHECK(status == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT && draw.count == 0 &&
            mvmc_classic_pfaffian_real_accepted(state_workspace)->generation ==
                1,
        "initialize/propose HOPPING divergence fails before prepare or draw");

  mvmc_classic_pfaffian_collective_workspace_destroy(collective_workspace);
  mvmc_classic_pfaffian_real_workspace_destroy(state_workspace);
}

static void test_collective_mirror_preflight(void) {
  const int initial_idx[2] = {0, 0};
  const int hopping_idx[2] = {1, 0};
  const double complex weights[2] = {1.0, 1.0};
  double slater[32];
  double legacy[10];
  MVMCClassicPfaffianRealWorkspace *state_workspace = NULL;
  MVMCClassicPfaffianCollectiveWorkspace *collective_workspace = NULL;
  MVMCClassicSamplerState sampler_state;
  MVMCClassicProposalResult result;
  DrawContext draw = {0.0, 0};
  size_t local_count = 10;
  int qp_start, qp_end;

  rank_partition(2, &qp_start, &qp_end);
  fill_real_slater(slater, 2);
  memset(legacy, 0, sizeof(legacy));
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 2, 2, qp_start, qp_end, &state_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create mirror preflight state workspace");
  CHECK(mvmc_classic_pfaffian_collective_workspace_create(
            TEST_COMMUNICATOR, &collective_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create mirror preflight collective workspace");
  if (world_size == 1 || world_rank == world_size - 1) local_count = 9;
  CHECK(mvmc_classic_sampler_real_initialize(
            state_workspace, collective_workspace, slater, initial_idx,
            weights, 0.0, 0.0, legacy, local_count, &sampler_state) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT && !sampler_state.valid &&
            !mvmc_classic_pfaffian_real_accepted(state_workspace)->valid,
        "one-rank short initialization mirror fails before any publish");

  CHECK(mvmc_classic_sampler_real_initialize(
            state_workspace, collective_workspace, slater, initial_idx,
            weights, 0.0, 0.0, legacy, 10, &sampler_state) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "valid mirror initializes after collective preflight failure");
  draw.count = 0;
  CHECK(mvmc_classic_sampler_real_propose(
            state_workspace, collective_workspace, &sampler_state,
            MVMC_CLASSIC_MOVE_HOPPING, slater, hopping_idx, weights,
            0.0, 0.0, 0.0, 0.0, draw_once, &draw, legacy, local_count,
            &result) == MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT &&
            draw.count == 0 &&
            mvmc_classic_pfaffian_real_accepted(state_workspace)->generation ==
                1,
        "one-rank short proposal mirror fails before evaluation or draw");

  mvmc_classic_pfaffian_collective_workspace_destroy(collective_workspace);
  mvmc_classic_pfaffian_real_workspace_destroy(state_workspace);
}

static void test_zero_current_initialization(void) {
  const int ele_idx[2] = {0, 0};
  const double complex weights[1] = {1.0};
  double slater[16] = {0};
  double legacy[5] = {0};
  MVMCClassicPfaffianRealWorkspace *state_workspace = NULL;
  MVMCClassicPfaffianCollectiveWorkspace *collective_workspace = NULL;
  MVMCClassicSamplerState sampler_state;
  int qp_start, qp_end;

  rank_partition(1, &qp_start, &qp_end);
  CHECK(mvmc_classic_pfaffian_real_workspace_create(
            2, 2, 1, qp_start, qp_end, &state_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create zero-current workspace");
  CHECK(mvmc_classic_pfaffian_collective_workspace_create(
            TEST_COMMUNICATOR, &collective_workspace) ==
            MVMC_PFAFFIAN_STATUS_OK,
        "create zero-current collective workspace");
  CHECK(mvmc_classic_sampler_real_initialize(
            state_workspace, collective_workspace, slater, ele_idx, weights,
            0.0, 0.0, legacy, 5, &sampler_state) ==
            MVMC_PFAFFIAN_STATUS_INVALID_ARGUMENT && !sampler_state.valid &&
            !mvmc_classic_pfaffian_real_accepted(state_workspace)->valid,
        "exact-zero current total is rejected before sampler start");
  mvmc_classic_pfaffian_collective_workspace_destroy(collective_workspace);
  mvmc_classic_pfaffian_real_workspace_destroy(state_workspace);
}

int main(int argc, char **argv) {
  ProposalTrace trace_a;
  ProposalTrace trace_b;
#ifdef _mpi_use
  int global_failures = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#else
  (void)argc;
  (void)argv;
#endif

  run_real_sequence(1, &trace_a);
  run_real_sequence(4, &trace_a);
  run_real_sequence(4, &trace_b);
  CHECK(memcmp(&trace_a, &trace_b, sizeof(trace_a)) == 0,
        "fixed-input RNG A/B replay is byte-identical");
  run_complex_rbm_sequence(1);
  run_complex_rbm_sequence(4);
  test_p2d_audit_trajectory_identity();
  test_p2d_guarded_fallback();
  test_p2d_audit_mismatch_and_rebuild_no_heal();
  test_p2d_complex_audit_smoke();
  test_p2d_rank_divergent_audit_mode();
  test_zero_current_initialization();
  test_collective_preflight_mismatch();
  test_log_domain_metropolis_edges();
  test_initialize_propose_preflight_mismatch();
  test_collective_mirror_preflight();

#ifdef _mpi_use
  MPI_Allreduce(&failure_count, &global_failures, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  if (world_rank == 0 && global_failures == 0) {
    printf("classic Pfaffian sampler unit: PASS (%d ranks)\n", world_size);
  }
  MPI_Finalize();
  return global_failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
#else
  if (failure_count == 0) printf("classic Pfaffian sampler unit: PASS\n");
  return failure_count == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
#endif
}
