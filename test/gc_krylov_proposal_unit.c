/*
 * Exact contract of the grand-canonical symmetric proposals: every draw
 * sequence is enumerated for L = 1, 2, 3 and both GC domains, the resulting
 * transition probabilities are exact rationals and are compared with an
 * independent combinatorial expectation (detailed balance Q[x,y] = Q[y,x],
 * unit row sums, uniformity, initial-N conditioning, connectivity of the
 * 15/16 local + 1/16 uniform mixture).  Also: L = 33 word boundary,
 * vacuum/full filling, callback errors, out-of-range draws, invalid input.
 *
 * usage: gc_krylov_proposal_unit [--dump path]   (dump for runtest_gc_lanczos)
 */
#include "gc_krylov_proposal.h"

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
} Script;

static MVMCKrylovStatus script_draw(void *context, size_t bound,
                                    size_t *value) {
  Script *script = (Script *)context;
  if (script->used >= MAX_DEPTH || bound == 0) {
    script->inconsistent = 1;
    return MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
  }
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

typedef enum { KIND_UNIFORM, KIND_INITIAL, KIND_LOCAL, KIND_TWO_PAIRS } Kind;

typedef struct {
  const MVMCKrylovFockModel *model;
  Kind kind;
  int initial_n;
  uint64_t current;
} Proposal;

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
  if (dump_stream != NULL) fclose(dump_stream);
  printf(failure_count == 0 ? "gc_krylov_proposal_unit: PASS\n"
                            : "gc_krylov_proposal_unit: FAIL\n");
  return failure_count == 0 ? 0 : 1;
}
