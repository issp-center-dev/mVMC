/*
 * Exact contract of the grand-canonical symmetric proposals: every draw
 * sequence is enumerated for L = 1, 2, 3 and both GC domains, the resulting
 * transition probabilities are exact rationals and are compared with an
 * independent combinatorial expectation (detailed balance Q[x,y] = Q[y,x],
 * unit row sums, uniformity, initial-N conditioning, connectivity of the
 * 15/16 local + 1/16 uniform mixture).  Also: L = 33 word boundary,
 * vacuum/full filling, callback errors, out-of-range draws, invalid input.
 * Finite-S start: the sector selector against an independent reference,
 * and every bounded draw sequence of the structured initial draw against
 * the frozen leaf counts/probabilities and the combinatorial formula, with
 * every call bound >= 2 (a draw(1) mutation is rejected by that check).
 *
 * usage: gc_krylov_proposal_unit [--dump path]   (dump for runtest_gc_lanczos)
 */
#include "gc_krylov_proposal.h"
#include "gc_lanczos_task5_fixture.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failure_count = 0;

#define CHECK(condition, message)                                    \
  do {                                                               \
    if (!(condition)) {                                              \
      fprintf(stderr, "FAIL: %s (line %d)\n", (message), __LINE__); \
      ++failure_count;                                               \
    }                                                                \
  } while (0)

#define MAX_DEPTH 256
#define MAX_STATES 64

typedef struct {
  uint64_t num;
  uint64_t den;
} Fraction;

static uint64_t gcd64(uint64_t a, uint64_t b) {
  while (b) {
    const uint64_t t = a % b;
    a = b;
    b = t;
  }
  return a;
}

static Fraction fraction(uint64_t num, uint64_t den) {
  Fraction value;
  const uint64_t g = gcd64(num, den);
  value.num = g ? num / g : 0;
  value.den = g ? den / g : 1;
  if (value.num == 0) value.den = 1;
  return value;
}

static Fraction fraction_add(Fraction a, Fraction b) {
  const uint64_t g = gcd64(a.den, b.den);
  return fraction(a.num * (b.den / g) + b.num * (a.den / g),
                  a.den / g * b.den);
}

static uint64_t last_leaves = 0;
static int last_small_bound = 0;

static Fraction fraction_mul(Fraction a, Fraction b) {
  const uint64_t g1 = gcd64(a.num, b.den);
  const uint64_t g2 = gcd64(b.num, a.den);
  if (a.num == 0 || b.num == 0) return fraction(0, 1);
  return fraction((a.num / g1) * (b.num / g2), (a.den / g2) * (b.den / g1));
}

static int fraction_equal(Fraction a, Fraction b) {
  return a.num == b.num && a.den == b.den;
}

typedef struct {
  size_t path[MAX_DEPTH];
  size_t bounds[MAX_DEPTH];
  size_t depth;
  size_t used;
  int inconsistent;
  int small_bound; /* a call with bound < 2 was seen */
  uint64_t leaves;
} Script;

static MVMCKrylovStatus script_draw(void *context, size_t bound,
                                    size_t *value) {
  Script *script = (Script *)context;
  if (script->used >= MAX_DEPTH || bound == 0) {
    script->inconsistent = 1;
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
  /* Checked on every call, also before the prefix is extended. */
  if (bound < 2) script->small_bound = 1;
  if (script->used < script->depth) {
    if (script->bounds[script->used] != bound) script->inconsistent = 1;
  } else {
    script->bounds[script->used] = bound;
    script->path[script->used] = 0;
    script->depth = script->used + 1;
  }
  *value = script->path[script->used++];
  return MVMC_KRYLOV_STATUS_OK;
}

typedef enum {
  KIND_UNIFORM,
  KIND_INITIAL,
  KIND_LOCAL,
  KIND_TWO_PAIRS,
  KIND_SPIN_INITIAL,
  KIND_SPIN_INITIAL_MUTATED
} Kind;

typedef struct {
  const MVMCKrylovFockModel *model;
  Kind kind;
  int initial_n;
  uint64_t current;
  MVMCGCSpinSupport support;
} Proposal;

/* Test-local copy of the structured draw with one extra draw(1) for a
 * single D value: the same leaves and probabilities, a different RNG use. */
static MVMCKrylovStatus mutated_spin_initial(
    const Proposal *proposal, MVMCKrylovFockProposalBoundedDraw draw,
    void *context, uint64_t *out) {
  const int L = (int)proposal->model->site_count;
  const int m = proposal->initial_n / 2, S = proposal->support.total_spin;
  const int low = 2 * m - L > 0 ? 2 * m - L : 0, high = m - S;
  int D = low, pass, site;
  size_t value;
  if (draw(context, high > low ? (size_t)(high - low + 1) : 1, &value) !=
      MVMC_KRYLOV_STATUS_OK)
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  if (high > low) D = low + (int)value;
  *out = 0;
  for (pass = 0; pass < 3; ++pass) {
    int k = pass == 0 ? D : m - D, remaining = 0;
    for (site = 0; site < L; ++site)
      remaining += !((*out >> site) & 1) && !((*out >> (site + L)) & 1);
    for (site = 0; site < L && k > 0; ++site) {
      int selected;
      if (((*out >> site) & 1) || ((*out >> (site + L)) & 1)) continue;
      if (k == remaining) {
        selected = 1;
      } else {
        if (draw(context, (size_t)remaining, &value) != MVMC_KRYLOV_STATUS_OK)
          return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
        selected = (int)value < k;
      }
      if (selected) {
        if (pass != 2) *out |= UINT64_C(1) << site;
        if (pass != 1) *out |= UINT64_C(1) << (site + L);
        --k;
      }
      --remaining;
    }
  }
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus run_proposal(const Proposal *proposal, Script *script,
                                     uint64_t *out) {
  switch (proposal->kind) {
    case KIND_UNIFORM:
      return mvmc_gc_krylov_draw_uniform(proposal->model, script_draw, script,
                                         out, 1);
    case KIND_INITIAL:
      return mvmc_gc_krylov_draw_initial(proposal->model,
                                         proposal->initial_n, script_draw,
                                         script, out, 1);
    case KIND_LOCAL:
      return mvmc_gc_krylov_draw_local(proposal->model, &proposal->current,
                                       script_draw, script, out, 1);
    case KIND_TWO_PAIRS:
      return mvmc_gc_krylov_draw_two_pairs(proposal->model, &proposal->current,
                                          script_draw, script, out, 1);
    case KIND_SPIN_INITIAL:
      return mvmc_gc_krylov_draw_spin_initial(
          proposal->model, &proposal->support, proposal->initial_n,
          script_draw, script, out, 1);
    case KIND_SPIN_INITIAL_MUTATED:
      return mutated_spin_initial(proposal, script_draw, script, out);
  }
  return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
}

/* Enumerate every draw sequence; probabilities into row[mask]. */
static int enumerate(const Proposal *proposal, Fraction *row,
                     size_t state_count) {
  Script script;
  size_t index;
  int ok = 1;
  memset(&script, 0, sizeof(script));
  for (index = 0; index < state_count; ++index) row[index] = fraction(0, 1);
  for (;;) {
    uint64_t out = UINT64_MAX;
    uint64_t denominator = 1;
    size_t k;
    MVMCKrylovStatus status;
    script.used = 0;
    status = run_proposal(proposal, &script, &out);
    if (status != MVMC_KRYLOV_STATUS_OK || script.inconsistent ||
        out >= state_count) {
      ok = 0;
      break;
    }
    script.depth = script.used;
    ++script.leaves;
    for (k = 0; k < script.used; ++k) denominator *= script.bounds[k];
    row[out] = fraction_add(row[out], fraction(1, denominator));
    for (k = script.used; k > 0; --k) {
      if (++script.path[k - 1] < script.bounds[k - 1]) {
        script.depth = k;
        break;
      }
    }
    if (k == 0) break;
  }
  last_leaves = script.leaves;
  last_small_bound = script.small_bound;
  return ok;
}

static int popcount64(uint64_t value) {
  int count = 0;
  while (value) {
    count += (int)(value & 1);
    value >>= 1;
  }
  return count;
}

static int in_domain(MVMCKrylovConfigurationDomain domain, uint64_t mask,
                     int nsite) {
  const uint64_t low = (UINT64_C(1) << nsite) - 1;
  if (domain == MVMC_KRYLOV_DOMAIN_GC_BALANCED) {
    return popcount64(mask & low) == popcount64(mask >> nsite);
  }
  return popcount64(mask) % 2 == 0;
}

static uint64_t binomial(int n, int k) {
  uint64_t value = 1;
  int i;
  if (k < 0 || k > n) return 0;
  for (i = 1; i <= k; ++i) value = value * (uint64_t)(n - k + i) / (uint64_t)i;
  return value;
}

static MVMCKrylovFockModel gc_model(MVMCKrylovConfigurationDomain domain,
                                    size_t nsite) {
  MVMCKrylovFockModel model;
  memset(&model, 0, sizeof(model));
  model.site_count = nsite;
  model.hermitian = 1;
  model.domain = domain;
  return model;
}

/* Independent combinatorial expectation of the local kernel. */
static Fraction expected_local(MVMCKrylovConfigurationDomain domain,
                               int nsite, uint64_t x, uint64_t y) {
  const uint64_t diff = x ^ y;
  if (domain == MVMC_KRYLOV_DOMAIN_GC_EVEN) {
    const uint64_t pairs = (uint64_t)(2 * nsite) * (uint64_t)(2 * nsite - 1);
    return popcount64(diff) == 2 ? fraction(2, pairs) : fraction(0, 1);
  }
  {
    const uint64_t low = (UINT64_C(1) << nsite) - 1;
    const uint64_t up_diff = diff & low;
    const uint64_t down_diff = diff >> nsite;
    if (x == y) {
      /* Self: everything that does not move. */
      Fraction moved = fraction(0, 1);
      uint64_t other;
      for (other = 0; other < (UINT64_C(1) << (2 * nsite)); ++other) {
        if (other != x && in_domain(domain, other, nsite)) {
          moved = fraction_add(moved, expected_local(domain, nsite, x, other));
        }
      }
      return fraction(moved.den - moved.num, moved.den);
    }
    if (popcount64(up_diff) == 1 && popcount64(down_diff) == 1) {
      /* Pair flip: both added or both removed. */
      const int up_was = (x & up_diff) != 0;
      const int down_was = ((x >> nsite) & down_diff) != 0;
      if (up_was == down_was) {
        return fraction(1, 2 * (uint64_t)nsite * (uint64_t)nsite);
      }
      return fraction(0, 1);
    }
    if (nsite >= 2 &&
        ((popcount64(up_diff) == 2 && down_diff == 0 &&
          popcount64(x & up_diff) == 1) ||
         (popcount64(down_diff) == 2 && up_diff == 0 &&
          popcount64((x >> nsite) & down_diff) == 1))) {
      return fraction(1, 2 * (uint64_t)nsite * (uint64_t)(nsite - 1));
    }
    return fraction(0, 1);
  }
}

static FILE *dump_stream = NULL;

static void dump_row(const char *kind, int domain, int nsite, int extra,
                     uint64_t x, const Fraction *row, size_t state_count) {
  size_t y;
  if (dump_stream == NULL) return;
  for (y = 0; y < state_count; ++y) {
    if (row[y].num != 0) {
      fprintf(dump_stream, "%s %d %d %d %llu %llu %llu %llu\n", kind, domain,
              nsite, extra, (unsigned long long)x, (unsigned long long)y,
              (unsigned long long)row[y].num, (unsigned long long)row[y].den);
    }
  }
}

static void test_exact_kernels(void) {
  const MVMCKrylovConfigurationDomain domains[2] = {
      MVMC_KRYLOV_DOMAIN_GC_BALANCED, MVMC_KRYLOV_DOMAIN_GC_EVEN};
  static Fraction local[MAX_STATES][MAX_STATES];
  int d, nsite;
  for (d = 0; d < 2; ++d) {
    for (nsite = 1; nsite <= 3; ++nsite) {
      const MVMCKrylovFockModel model = gc_model(domains[d], (size_t)nsite);
      const size_t state_count = (size_t)1 << (2 * nsite);
      Fraction row[MAX_STATES];
      Proposal proposal;
      uint64_t domain_size = 0, x, y;
      int n;
      for (x = 0; x < state_count; ++x) {
        domain_size += (uint64_t)in_domain(domains[d], x, nsite);
      }
      CHECK(domain_size ==
                (domains[d] == MVMC_KRYLOV_DOMAIN_GC_EVEN
                     ? (UINT64_C(1) << (2 * nsite - 1))
                     : binomial(2 * nsite, nsite)),
            "domain size");
      /* Uniform over the whole domain. */
      memset(&proposal, 0, sizeof(proposal));
      proposal.model = &model;
      proposal.kind = KIND_UNIFORM;
      CHECK(enumerate(&proposal, row, state_count), "uniform enumeration");
      dump_row("UNIFORM", (int)domains[d], nsite, 0, 0, row, state_count);
      for (y = 0; y < state_count; ++y) {
        CHECK(fraction_equal(row[y], in_domain(domains[d], y, nsite)
                                         ? fraction(1, domain_size)
                                         : fraction(0, 1)),
              "uniform probability is constant on the domain");
      }
      /* Uniform conditioned on N. */
      for (n = 0; n <= 2 * nsite; n += 2) {
        uint64_t sector = 0;
        proposal.kind = KIND_INITIAL;
        proposal.initial_n = n;
        CHECK(enumerate(&proposal, row, state_count), "initial enumeration");
        dump_row("INITIAL", (int)domains[d], nsite, n, 0, row, state_count);
        sector = domains[d] == MVMC_KRYLOV_DOMAIN_GC_EVEN
                     ? binomial(2 * nsite, n)
                     : binomial(nsite, n / 2) * binomial(nsite, n / 2);
        for (y = 0; y < state_count; ++y) {
          const int member = in_domain(domains[d], y, nsite) &&
                             popcount64(y) == n;
          CHECK(fraction_equal(row[y], member ? fraction(1, sector)
                                              : fraction(0, 1)),
                "initial-N probability");
        }
      }
      /* Local kernel. */
      for (x = 0; x < state_count; ++x) {
        Fraction sum = fraction(0, 1);
        if (!in_domain(domains[d], x, nsite)) continue;
        proposal.kind = KIND_LOCAL;
        proposal.current = x;
        CHECK(enumerate(&proposal, local[x], state_count),
              "local enumeration");
        dump_row("LOCAL", (int)domains[d], nsite, 0, x, local[x],
                 state_count);
        for (y = 0; y < state_count; ++y) {
          if (!in_domain(domains[d], y, nsite)) {
            CHECK(local[x][y].num == 0, "local stays in the domain");
            continue;
          }
          CHECK(fraction_equal(local[x][y],
                               expected_local(domains[d], nsite, x, y)),
                "local probability equals the combinatorial expectation");
          sum = fraction_add(sum, local[x][y]);
        }
        CHECK(fraction_equal(sum, fraction(1, 1)), "local row sums to one");
      }
      for (x = 0; x < state_count; ++x) {
        for (y = 0; y < state_count; ++y) {
          if (in_domain(domains[d], x, nsite) &&
              in_domain(domains[d], y, nsite)) {
            const Fraction uniform = fraction(1, 16 * domain_size);
            const Fraction forward = fraction_add(
                fraction_mul(fraction(15, 16), local[x][y]), uniform);
            const Fraction reverse = fraction_add(
                fraction_mul(fraction(15, 16), local[y][x]), uniform);
            CHECK(fraction_equal(local[x][y], local[y][x]),
                  "detailed balance Q[x,y] = Q[y,x]");
            /* 15/16 local + 1/16 uniform: symmetric and positive on D x D,
             * so every positive-amplitude subset is connected directly. */
            CHECK(fraction_equal(forward, reverse) && forward.num != 0,
                  "mixture kernel is symmetric and connects the domain");
          }
        }
      }
    }
  }
}

typedef struct {
  unsigned seed;
  size_t used;
} Sequence;

/* Deterministic in-range integers for the word-boundary checks. */
static MVMCKrylovStatus sequence_draw(void *context, size_t bound,
                                      size_t *value) {
  Sequence *sequence = (Sequence *)context;
  sequence->seed = sequence->seed * 1103515245u + 12345u;
  *value = (size_t)((sequence->seed >> 8) % bound);
  ++sequence->used;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus failing_draw(void *context, size_t bound,
                                     size_t *value) {
  (void)context;
  (void)bound;
  *value = 0;
  return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
}

static MVMCKrylovStatus out_of_range_draw(void *context, size_t bound,
                                          size_t *value) {
  (void)context;
  *value = bound;
  return MVMC_KRYLOV_STATUS_OK;
}

static void test_word_boundary(void) {
  const MVMCKrylovConfigurationDomain domains[2] = {
      MVMC_KRYLOV_DOMAIN_GC_BALANCED, MVMC_KRYLOV_DOMAIN_GC_EVEN};
  int d;
  for (d = 0; d < 2; ++d) {
    const MVMCKrylovFockModel model = gc_model(domains[d], 33);
    uint64_t words[2], next[2];
    Sequence sequence;
    Script zeros;
    int round;
    sequence.seed = 17u + (unsigned)d;
    sequence.used = 0;
    for (round = 0; round < 256; ++round) {
      CHECK(mvmc_gc_krylov_draw_uniform(&model, sequence_draw, &sequence,
                                        words, 2) == MVMC_KRYLOV_STATUS_OK &&
                mvmc_krylov_fock_validate(&model, words, 2) ==
                    MVMC_KRYLOV_STATUS_OK,
            "L=33 uniform draw is in the domain");
      CHECK(mvmc_gc_krylov_draw_local(&model, words, sequence_draw,
                                      &sequence, next, 2) ==
                    MVMC_KRYLOV_STATUS_OK &&
                mvmc_krylov_fock_validate(&model, next, 2) ==
                    MVMC_KRYLOV_STATUS_OK,
            "L=33 local move stays in the domain");
      CHECK(mvmc_gc_krylov_draw_initial(&model, 2 * (round % 34),
                                        sequence_draw, &sequence, words,
                                        2) == MVMC_KRYLOV_STATUS_OK &&
                popcount64(words[0]) + popcount64(words[1]) ==
                    2 * (round % 34) &&
                mvmc_krylov_fock_validate(&model, words, 2) ==
                    MVMC_KRYLOV_STATUS_OK,
            "L=33 initial-N draw");
    }
    memset(&zeros, 0, sizeof(zeros));
    CHECK(mvmc_gc_krylov_draw_initial(&model, 66, script_draw, &zeros, words,
                                      2) == MVMC_KRYLOV_STATUS_OK &&
              words[0] == UINT64_MAX && words[1] == UINT64_C(3),
          "L=33 full filling");
    memset(&zeros, 0, sizeof(zeros));
    CHECK(mvmc_gc_krylov_draw_local(&model, words, script_draw, &zeros, next,
                                    2) == MVMC_KRYLOV_STATUS_OK &&
              popcount64(next[0]) + popcount64(next[1]) == 64,
          "full-filling local move removes a pair");
    memset(&zeros, 0, sizeof(zeros));
    CHECK(mvmc_gc_krylov_draw_initial(&model, 0, script_draw, &zeros, words,
                                      2) == MVMC_KRYLOV_STATUS_OK &&
              words[0] == 0 && words[1] == 0,
          "L=33 vacuum");
    memset(&zeros, 0, sizeof(zeros));
    CHECK(mvmc_gc_krylov_draw_local(&model, words, script_draw, &zeros, next,
                                    2) == MVMC_KRYLOV_STATUS_OK &&
              popcount64(next[0]) + popcount64(next[1]) == 2,
          "vacuum local move adds a pair");
  }
}

static void test_errors(void) {
  const MVMCKrylovFockModel even = gc_model(MVMC_KRYLOV_DOMAIN_GC_EVEN, 2);
  const MVMCKrylovFockModel balanced =
      gc_model(MVMC_KRYLOV_DOMAIN_GC_BALANCED, 2);
  MVMCKrylovFockModel fixed = gc_model(MVMC_KRYLOV_DOMAIN_FIXED, 2);
  MVMCKrylovFockModel invalid = gc_model((MVMCKrylovConfigurationDomain)5, 2);
  Script script;
  uint64_t words[2] = {UINT64_C(0x5), 0};
  uint64_t out = 0;
  memset(&script, 0, sizeof(script));
  fixed.up_electron_count = 1;
  fixed.down_electron_count = 1;
  CHECK(mvmc_gc_krylov_draw_uniform(&even, failing_draw, NULL, &out, 1) ==
            MVMC_KRYLOV_STATUS_RESOURCE_LIMIT,
        "callback error is returned");
  CHECK(mvmc_gc_krylov_draw_local(&balanced, words, out_of_range_draw, NULL,
                                  &out, 1) ==
            MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE,
        "out-of-range draw is an invariant failure");
  CHECK(mvmc_gc_krylov_draw_initial(&even, 2, out_of_range_draw, NULL, &out,
                                    1) ==
            MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE,
        "out-of-range initial draw");
  CHECK(mvmc_gc_krylov_draw_uniform(&fixed, script_draw, &script, &out, 1) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "fixed domain rejected");
  CHECK(mvmc_gc_krylov_draw_uniform(&invalid, script_draw, &script, &out,
                                    1) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "invalid domain rejected");
  CHECK(mvmc_gc_krylov_draw_initial(&even, -2, script_draw, &script, &out,
                                    1) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "negative initial N rejected");
  CHECK(mvmc_gc_krylov_draw_initial(&even, 3, script_draw, &script, &out, 1) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "odd initial N rejected");
  CHECK(mvmc_gc_krylov_draw_initial(&even, 6, script_draw, &script, &out, 1) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "too large initial N rejected");
  CHECK(mvmc_gc_krylov_draw_local(&even, words, script_draw, &script, words,
                                  1) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "aliased current/out rejected");
  words[0] = UINT64_C(0x1);
  CHECK(mvmc_gc_krylov_draw_local(&even, words, script_draw, &script, &out,
                                  1) == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "odd current configuration rejected");
  CHECK(mvmc_gc_krylov_draw_uniform(&even, script_draw, &script, &out, 2) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "word count rejected");
  CHECK(mvmc_gc_krylov_draw_uniform(&even, NULL, &script, &out, 1) ==
            MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
        "missing draw callback rejected");
}

/* Enumerate every random path and compare with the independent closed form
 * Q(x,y)=1/binom(L,2)^2 for valid additions/removals. Detailed balance is
 * exact, including unavailable moves and vacuum/full filling. */
static void test_two_pair_kernel(void) {
  int sites;
  for (sites = 1; sites <= 4; ++sites) {
    MVMCKrylovFockModel model = gc_model(MVMC_KRYLOV_DOMAIN_GC_BALANCED,
                                         (size_t)sites);
    const size_t states = (size_t)1 << (2 * sites);
    Fraction *matrix = calloc(states * states, sizeof(*matrix));
    Proposal proposal;
    size_t x, y;
    memset(&proposal, 0, sizeof(proposal));
    proposal.model = &model;
    proposal.kind = KIND_TWO_PAIRS;
    CHECK(matrix != NULL, "two-pair matrix allocation");
    if (matrix == NULL) return;
    for (x = 0; x < states; ++x) {
      Fraction off_diagonal = fraction(0, 1);
      if (!in_domain(model.domain, x, sites)) continue;
      proposal.current = x;
      CHECK(enumerate(&proposal, matrix + x * states, states),
            "two-pair enumeration");
      for (y = 0; y < states; ++y) {
        const uint64_t difference = x ^ y;
        const uint64_t low = (UINT64_C(1) << sites) - 1;
        const int changed_up = popcount64(difference & low);
        const int changed_down = popcount64(difference >> sites);
        const int added = (x & difference) == 0;
        const int removed = (y & difference) == 0;
        Fraction expected = fraction(0, 1);
        if (x == y) continue;
        if (changed_up == 2 && changed_down == 2 && (added || removed)) {
          const uint64_t pairs = binomial(sites, 2);
          expected = fraction(1, pairs * pairs);
        }
        CHECK(fraction_equal(matrix[x * states + y], expected),
              "two-pair probability agrees with combinatorics");
        off_diagonal = fraction_add(off_diagonal, expected);
      }
      CHECK(fraction_equal(matrix[x * states + x],
                fraction(off_diagonal.den - off_diagonal.num, off_diagonal.den)),
            "two-pair self-loop completes unit row sum");
    }
    for (x = 0; x < states; ++x)
      for (y = 0; y < states; ++y)
        if (in_domain(model.domain, x, sites) &&
            in_domain(model.domain, y, sites))
          CHECK(fraction_equal(matrix[x * states + y], matrix[y * states + x]),
                "two-pair exact detailed balance");
    free(matrix);
  }
}

/* ---- Finite-S start (Task 5) -------------------------------------------- */

static MVMCGCSpinSupport spin_support(int spin, int symmetry) {
  MVMCGCSpinSupport support;
  support.active = 1;
  support.gauss_count = 4;
  support.total_spin = spin;
  support.pair_symmetry = (MVMCGCPairSymmetry)symmetry;
  return support;
}

static int reference_allowed(int L, int S, int triplet, int n) {
  const int m = n / 2, low = 2 * m - L > 0 ? 2 * m - L : 0;
  return n >= 0 && n <= 2 * L && n % 2 == 0 && low <= m - S &&
         (!triplet || (m - S) % 2 == 0);
}

static void test_initial_sectors(void) {
  int sectors[16], expected[16];
  size_t count = 0, k;
  int L, S, sym, requested, explicit_n;
  MVMCKrylovFockModel model = gc_model(MVMC_KRYLOV_DOMAIN_GC_BALANCED, 4);
  MVMCGCSpinSupport support = spin_support(1, MVMC_GC_PAIR_TRIPLET);
  CHECK(mvmc_gc_krylov_initial_sectors(&model, &support, 4, 0, sectors, 5,
                                       &count) == MVMC_KRYLOV_STATUS_OK &&
            count == 2 && sectors[0] == 2 && sectors[1] == 6,
        "L4/S1/triplet auto default 4 -> [2,6]");
  model = gc_model(MVMC_KRYLOV_DOMAIN_GC_BALANCED, 5);
  support = spin_support(2, MVMC_GC_PAIR_GENERAL);
  CHECK(mvmc_gc_krylov_initial_sectors(&model, &support, 4, 0, sectors, 6,
                                       &count) == MVMC_KRYLOV_STATUS_OK &&
            count == 2 && sectors[0] == 4 && sectors[1] == 6,
        "L5/S2/general auto default 4 -> [4,6]");
  model = gc_model(MVMC_KRYLOV_DOMAIN_GC_BALANCED, 4);
  support = spin_support(2, MVMC_GC_PAIR_GENERAL);
  CHECK(mvmc_gc_krylov_initial_sectors(&model, &support, 0, 0, sectors, 5,
                                       &count) == MVMC_KRYLOV_STATUS_OK &&
            count == 1 && sectors[0] == 4,
        "helper example L4/S2 requested 0 auto -> [4]");
  support = spin_support(1, MVMC_GC_PAIR_TRIPLET);
  sectors[0] = -7;
  count = 77;
  CHECK(mvmc_gc_krylov_initial_sectors(&model, &support, 4, 1, sectors, 5,
                                       &count) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            sectors[0] == -7 && count == 77,
        "explicit forbidden N4 rejected, outputs unchanged");
  CHECK(mvmc_gc_krylov_initial_sectors(&model, &support, 2, 1, sectors, 5,
                                       &count) == MVMC_KRYLOV_STATUS_OK &&
            count == 1 && sectors[0] == 2,
        "explicit allowed N2 kept alone");
  sectors[0] = -7;
  count = 77;
  CHECK(mvmc_gc_krylov_initial_sectors(&model, &support, 2, 1, sectors, 4,
                                       &count) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            mvmc_gc_krylov_initial_sectors(&model, &support, 3, 0, sectors, 5,
                                           &count) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            mvmc_gc_krylov_initial_sectors(&model, &support, 10, 0, sectors,
                                           5, &count) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            mvmc_gc_krylov_initial_sectors(&model, &support, -2, 0, sectors,
                                           5, &count) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            mvmc_gc_krylov_initial_sectors(&model, &support, 2, 2, sectors, 5,
                                           &count) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            mvmc_gc_krylov_initial_sectors(&model, &support, 2, 0, NULL, 5,
                                           &count) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
            sectors[0] == -7 && count == 77,
        "capacity < L+1, odd/out-of-range N, explicit flag, NULL rejected");
  {
    const MVMCGCSpinSupport bad[] = {
        {0, 1, 0, MVMC_GC_PAIR_GENERAL}, {1, 4, 0, MVMC_GC_PAIR_TRIPLET},
        {1, 4, 3, MVMC_GC_PAIR_GENERAL}, {1, 1, 1, MVMC_GC_PAIR_GENERAL}};
    for (k = 0; k < sizeof(bad) / sizeof(bad[0]); ++k) {
      CHECK(mvmc_gc_krylov_initial_sectors(&model, bad + k, 2, 0, sectors, 5,
                                           &count) ==
                MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
            "inactive, S=0 and invalid descriptors rejected");
    }
    support = spin_support(1, MVMC_GC_PAIR_SINGLET);
    CHECK(mvmc_gc_krylov_initial_sectors(&model, &support, 2, 0, sectors, 5,
                                         &count) ==
              MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE,
          "S>0 singlet stops with AMPLITUDE_FAILURE");
    support = spin_support(1, MVMC_GC_PAIR_ZERO);
    CHECK(mvmc_gc_krylov_initial_sectors(&model, &support, 2, 0, sectors, 5,
                                         &count) ==
              MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE,
          "S>0 zero stops with AMPLITUDE_FAILURE");
    model = gc_model(MVMC_KRYLOV_DOMAIN_GC_EVEN, 4);
    support = spin_support(1, MVMC_GC_PAIR_GENERAL);
    CHECK(mvmc_gc_krylov_initial_sectors(&model, &support, 2, 0, sectors, 5,
                                         &count) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "GC_EVEN rejected");
  }
  /* Exhaustive against the reference ordering (|N - r|, N). */
  for (L = 1; L <= 6; ++L) {
    model = gc_model(MVMC_KRYLOV_DOMAIN_GC_BALANCED, (size_t)L);
    for (S = 1; S <= L / 2; ++S) {
      for (sym = 0; sym < 2; ++sym) {
        const int triplet = sym == 1;
        support = spin_support(S, triplet ? MVMC_GC_PAIR_TRIPLET
                                          : MVMC_GC_PAIR_GENERAL);
        for (requested = 0; requested <= 2 * L; requested += 2) {
          for (explicit_n = 0; explicit_n < 2; ++explicit_n) {
            size_t n_expected = 0;
            MVMCKrylovStatus status;
            int distance;
            if (explicit_n) {
              if (reference_allowed(L, S, triplet, requested))
                expected[n_expected++] = requested;
            } else {
              for (distance = 0; distance <= 2 * L; ++distance) {
                if (reference_allowed(L, S, triplet, requested - distance))
                  expected[n_expected++] = requested - distance;
                if (distance > 0 &&
                    reference_allowed(L, S, triplet, requested + distance))
                  expected[n_expected++] = requested + distance;
              }
            }
            count = 99;
            status = mvmc_gc_krylov_initial_sectors(
                &model, &support, requested, explicit_n, sectors,
                (size_t)L + 1, &count);
            if (n_expected == 0) {
              CHECK(explicit_n && status == MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
                        count == 99,
                    "only an explicit forbidden N has no sector");
            } else {
              CHECK(status == MVMC_KRYLOV_STATUS_OK && count == n_expected &&
                        memcmp(sectors, expected,
                               n_expected * sizeof(*sectors)) == 0,
                    "sector list equals the independent ordering");
            }
          }
        }
      }
    }
  }
}

static uint64_t task5_choose(int n, int k) { return binomial(n, k); }

static void test_spin_initial_draw(void) {
  size_t item;
  int L, S, n;
  /* Frozen evidence: leaf counts and configuration probabilities. */
  for (item = 0; item < sizeof(task5_draws) / sizeof(task5_draws[0]); ++item) {
    const Task5Draw *draw = task5_draws + item;
    MVMCKrylovFockModel model =
        gc_model(MVMC_KRYLOV_DOMAIN_GC_BALANCED, (size_t)draw->site_count);
    const size_t states = (size_t)1 << (2 * draw->site_count);
    Fraction *row = calloc(states, sizeof(*row));
    Proposal proposal;
    size_t k, x;
    Fraction total = fraction(0, 1);
    int pass;
    CHECK(row != NULL, "row allocation");
    if (row == NULL) return;
    for (pass = 0; pass < 2; ++pass) {
      memset(&proposal, 0, sizeof(proposal));
      proposal.model = &model;
      proposal.kind = pass == 0 ? KIND_SPIN_INITIAL : KIND_SPIN_INITIAL_MUTATED;
      proposal.initial_n = draw->electron_count;
      proposal.support = spin_support(draw->total_spin, MVMC_GC_PAIR_GENERAL);
      CHECK(enumerate(&proposal, row, states), "spin initial enumeration");
      total = fraction(0, 1);
      for (k = 0; k < draw->count; ++k) {
        CHECK(fraction_equal(row[draw->probability[k].configuration],
                             fraction(draw->probability[k].numerator,
                                      draw->probability[k].denominator)),
              "probability equals the frozen evidence");
        total = fraction_add(total, row[draw->probability[k].configuration]);
      }
      for (x = 0; x < states; ++x) {
        int listed = 0;
        for (k = 0; k < draw->count; ++k)
          listed = listed || draw->probability[k].configuration == x;
        CHECK(listed || row[x].num == 0, "no unlisted configuration");
      }
      CHECK(fraction_equal(total, fraction(1, 1)), "probabilities sum to 1");
      CHECK(last_leaves == draw->leaves, "leaf count equals the evidence");
      if (pass == 0) {
        CHECK(!last_small_bound, "every draw call has bound >= 2");
      } else {
        /* The mutation keeps leaves and probabilities; only the bound
         * check sees it (it applies when D has a single value). */
        const int m = draw->electron_count / 2;
        const int low = 2 * m - draw->site_count > 0 ? 2 * m - draw->site_count
                                                     : 0;
        CHECK(last_small_bound == (low == m - draw->total_spin),
              "a draw(1) mutation is rejected by the bound check");
      }
    }
    free(row);
  }
  /* Independent combinatorial formula for every allowed N up to L = 5. */
  for (L = 1; L <= 5; ++L) {
    MVMCKrylovFockModel model = gc_model(MVMC_KRYLOV_DOMAIN_GC_BALANCED,
                                         (size_t)L);
    const size_t states = (size_t)1 << (2 * L);
    Fraction *row = calloc(states, sizeof(*row));
    CHECK(row != NULL, "row allocation");
    if (row == NULL) return;
    for (S = 1; S <= L / 2; ++S) {
      for (n = 0; n <= 2 * L; n += 2) {
        Proposal proposal;
        const int m = n / 2, low = 2 * m - L > 0 ? 2 * m - L : 0;
        const int high = m - S;
        size_t x;
        if (!reference_allowed(L, S, 0, n)) continue;
        memset(&proposal, 0, sizeof(proposal));
        proposal.model = &model;
        proposal.kind = KIND_SPIN_INITIAL;
        proposal.initial_n = n;
        proposal.support = spin_support(S, MVMC_GC_PAIR_GENERAL);
        CHECK(enumerate(&proposal, row, states) && !last_small_bound,
              "structured draw enumeration");
        for (x = 0; x < states; ++x) {
          const uint64_t low_mask = (UINT64_C(1) << L) - 1;
          const uint64_t up = x & low_mask, down = (x >> L) & low_mask;
          const int D = popcount64(up & down);
          const int ups = popcount64(up & ~down), downs = popcount64(down & ~up);
          Fraction expected = fraction(0, 1);
          if (D >= low && D <= high && ups == m - D && downs == m - D) {
            expected = fraction(
                1, (uint64_t)(high - low + 1) * task5_choose(L, D) *
                       task5_choose(L - D, m - D) *
                       task5_choose(L - m, m - D));
          }
          CHECK(fraction_equal(row[x], expected),
                "structured draw equals the combinatorial probability");
        }
      }
    }
    free(row);
  }
  /* Errors: no draw is consumed for a rejected request. */
  {
    MVMCKrylovFockModel model = gc_model(MVMC_KRYLOV_DOMAIN_GC_BALANCED, 4);
    MVMCGCSpinSupport support = spin_support(1, MVMC_GC_PAIR_TRIPLET);
    Script script;
    uint64_t out = 0;
    memset(&script, 0, sizeof(script));
    CHECK(mvmc_gc_krylov_draw_spin_initial(&model, &support, 4, script_draw,
                                           &script, &out, 1) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              script.used == 0,
          "triplet parity-forbidden N4 is not drawn");
    CHECK(mvmc_gc_krylov_draw_spin_initial(&model, &support, 0, script_draw,
                                           &script, &out, 1) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              mvmc_gc_krylov_draw_spin_initial(&model, &support, 3,
                                               script_draw, &script, &out,
                                               1) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              mvmc_gc_krylov_draw_spin_initial(&model, &support, 10,
                                               script_draw, &script, &out,
                                               1) ==
                  MVMC_KRYLOV_STATUS_INVALID_ARGUMENT &&
              script.used == 0,
          "out-of-range and odd N are not drawn");
    support = spin_support(1, MVMC_GC_PAIR_SINGLET);
    CHECK(mvmc_gc_krylov_draw_spin_initial(&model, &support, 2, script_draw,
                                           &script, &out, 1) ==
              MVMC_KRYLOV_STATUS_AMPLITUDE_FAILURE,
          "S>0 singlet is an amplitude failure");
    support = mvmc_gc_spin_support_inactive();
    CHECK(mvmc_gc_krylov_draw_spin_initial(&model, &support, 2, script_draw,
                                           &script, &out, 1) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "inactive support rejected");
    support = spin_support(1, MVMC_GC_PAIR_GENERAL);
    CHECK(mvmc_gc_krylov_draw_spin_initial(&model, &support, 2, failing_draw,
                                           NULL, &out, 1) ==
                  MVMC_KRYLOV_STATUS_RESOURCE_LIMIT &&
              mvmc_gc_krylov_draw_spin_initial(&model, &support, 2,
                                               out_of_range_draw, NULL, &out,
                                               1) ==
                  MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE,
          "callback error propagated, out-of-range draw is an invariant "
          "failure");
    CHECK(mvmc_gc_krylov_draw_spin_initial(&model, &support, 2, script_draw,
                                           &script, &out, 2) ==
              MVMC_KRYLOV_STATUS_INVALID_ARGUMENT,
          "word count checked");
  }
}

int main(int argc, char **argv) {
  if (argc == 3 && strcmp(argv[1], "--dump") == 0) {
    dump_stream = fopen(argv[2], "w");
    if (dump_stream == NULL) return 2;
  } else if (argc != 1) {
    fprintf(stderr, "usage: gc_krylov_proposal_unit [--dump path]\n");
    return 2;
  }
  test_exact_kernels();
  test_word_boundary();
  test_errors();
  test_two_pair_kernel();
  test_initial_sectors();
  test_spin_initial_draw();
  if (dump_stream != NULL) fclose(dump_stream);
  printf(failure_count == 0 ? "gc_krylov_proposal_unit: PASS\n"
                            : "gc_krylov_proposal_unit: FAIL\n");
  return failure_count == 0 ? 0 : 1;
}
