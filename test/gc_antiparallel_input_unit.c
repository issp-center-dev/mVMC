#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "gc_antiparallel_input.h"

#define NSITE 2
#define NORB 4
#define FIDX 1
#define OPT_GUARD 3
#define OPT_SIZE (2 * (FIDX + NORB) + OPT_GUARD)
#define CANARY 0x5a5a5a5a
#define IDX_GUARD 2

static int failures = 0;

#define CHECK(condition, ...)                                                   \
  do {                                                                          \
    if (!(condition)) {                                                         \
      fprintf(stderr, "GC_AntiPair_Input_Unit FAIL: ");                        \
      fprintf(stderr, __VA_ARGS__);                                             \
      fprintf(stderr, "\n");                                                   \
      failures++;                                                               \
    }                                                                           \
  } while (0)

/* ---- stderr capture so each negative case is matched to its diagnostic. */
static int savedStderr = -1;
static FILE *captureFile = NULL;
static char captured[8192];

static void capture_begin(void) {
  fflush(stderr);
  captureFile = tmpfile();
  if (captureFile == NULL) {
    perror("tmpfile");
    exit(EXIT_FAILURE);
  }
  savedStderr = dup(fileno(stderr));
  dup2(fileno(captureFile), fileno(stderr));
}

static void capture_end(void) {
  size_t length;
  fflush(stderr);
  dup2(savedStderr, fileno(stderr));
  close(savedStderr);
  rewind(captureFile);
  length = fread(captured, 1, sizeof(captured) - 1, captureFile);
  captured[length] = '\0';
  fclose(captureFile);
  captureFile = NULL;
}

static FILE *file_from_text(const char *text) {
  FILE *fp = tmpfile();
  if (fp == NULL) {
    perror("tmpfile");
    exit(EXIT_FAILURE);
  }
  fputs(text, fp);
  rewind(fp);
  return fp;
}

/* Orbital tables live inside a guarded buffer so writes outside the
 * Nsite x Nsite rows are visible. */
typedef struct {
  int idxStorage[NSITE * (NSITE + 2 * IDX_GUARD)];
  int sgnStorage[NSITE * (NSITE + 2 * IDX_GUARD)];
  int *idx[NSITE];
  int *sgn[NSITE];
  int opt[OPT_SIZE];
  int optCount;
} Tables;

static void reset_tables(Tables *tables) {
  int row;
  int k;
  for (k = 0; k < NSITE * (NSITE + 2 * IDX_GUARD); k++) {
    tables->idxStorage[k] = CANARY;
    tables->sgnStorage[k] = CANARY;
  }
  for (row = 0; row < NSITE; row++) {
    tables->idx[row] =
        tables->idxStorage + row * (NSITE + 2 * IDX_GUARD) + IDX_GUARD;
    tables->sgn[row] =
        tables->sgnStorage + row * (NSITE + 2 * IDX_GUARD) + IDX_GUARD;
  }
  for (k = 0; k < OPT_SIZE; k++) tables->opt[k] = CANARY;
  tables->optCount = 7;
}

static int guards_intact(const Tables *tables) {
  int row;
  int k;
  for (row = 0; row < NSITE; row++) {
    const int base = row * (NSITE + 2 * IDX_GUARD);
    for (k = 0; k < IDX_GUARD; k++) {
      if (tables->idxStorage[base + k] != CANARY ||
          tables->sgnStorage[base + k] != CANARY ||
          tables->idxStorage[base + IDX_GUARD + NSITE + k] != CANARY ||
          tables->sgnStorage[base + IDX_GUARD + NSITE + k] != CANARY) {
        return 0;
      }
    }
  }
  for (k = 0; k < 2 * FIDX; k++) {
    if (tables->opt[k] != CANARY) return 0;
  }
  for (k = 2 * (FIDX + NORB); k < OPT_SIZE; k++) {
    if (tables->opt[k] != CANARY) return 0;
  }
  return 1;
}

static int run_body_tables(const char *text, const int ap, const int expectOk,
                           const char *expectedMessage, Tables *tables) {
  FILE *fp = file_from_text(text);
  int result;
  reset_tables(tables);
  capture_begin();
  result = GCAntiReadOrbitals(fp, tables->idx, tables->sgn, tables->opt,
                              &tables->optCount, FIDX, 1, ap, NSITE, NORB,
                              "orbitalidx.def");
  capture_end();
  CHECK(guards_intact(tables), "write outside the tables for input:\n%s",
        text);
  CHECK(fp != NULL && ftell(fp) >= 0, "reader closed the caller's FILE");
  fclose(fp);
  if (expectOk) {
    CHECK(result == 0, "expected success (ap=%d) but got: %s\ninput:\n%s", ap,
          captured, text);
    CHECK(captured[0] == '\0', "unexpected diagnostic: %s", captured);
  } else {
    CHECK(result != 0, "expected failure (ap=%d) for input:\n%s", ap, text);
    CHECK(strstr(captured, "GC anti-parallel") != NULL &&
              strstr(captured, "orbitalidx.def") != NULL,
          "diagnostic lacks prefix/file: %s", captured);
    if (expectedMessage != NULL) {
      CHECK(strstr(captured, expectedMessage) != NULL,
            "diagnostic '%s' missing; got: %s\ninput:\n%s", expectedMessage,
            captured, text);
    }
  }
  return result;
}

static int run_body(const char *text, const int ap, const int expectOk) {
  Tables tables;
  return run_body_tables(text, ap, expectOk, NULL, &tables);
}

static void run_body_error(const char *text, const int ap,
                           const char *expectedMessage) {
  Tables tables;
  (void)run_body_tables(text, ap, 0, expectedMessage, &tables);
}

static const char *GOOD_OPT = "0 1\n1 1\n2 1\n3 1\n";

static void positive_cases(void) {
  Tables tables;
  const char *mixed =
      "0 0 0 1\n0 1 1 -1\n1 0 2 1\n1 1 3 1\n3 0\n0 1\n2 1\n1 0\n\n";
  int i;
  int j;
  /* AP: signs are kept, OptFlag rows are matched by index. */
  CHECK(run_body_tables(mixed, 1, 1, NULL, &tables) == 0, "AP mixed body");
  CHECK(tables.idx[0][0] == 0 && tables.idx[0][1] == 1 &&
            tables.idx[1][0] == 2 && tables.idx[1][1] == 3,
        "AP parameter indices");
  CHECK(tables.sgn[0][1] == -1 && tables.sgn[0][0] == 1 &&
            tables.sgn[1][0] == 1 && tables.sgn[1][1] == 1,
        "AP signs");
  CHECK(tables.opt[2 * (FIDX + 3)] == 0 && tables.opt[2 * (FIDX + 3) + 1] == 0,
        "OptFlag index 3");
  CHECK(tables.opt[2 * (FIDX + 0)] == 1 && tables.opt[2 * (FIDX + 0) + 1] == 1,
        "OptFlag index 0");
  CHECK(tables.opt[2 * (FIDX + 2)] == 1 && tables.opt[2 * (FIDX + 1)] == 0,
        "OptFlag indices 1 and 2");
  CHECK(tables.optCount == 7 + NORB, "OptFlag count");
  /* PBC: the same four-column input reads, and every sign is +1. */
  CHECK(run_body_tables(mixed, 0, 1, NULL, &tables) == 0, "PBC mixed body");
  for (i = 0; i < NSITE; i++) {
    for (j = 0; j < NSITE; j++) {
      CHECK(tables.sgn[i][j] == 1, "PBC sign normalized at (%d,%d)", i, j);
    }
  }
  /* The eight accepted forms. */
  run_body("0 0 0\n0 1 1\n1 0 2\n1 1 3\n0 1\n1 1\n2 1\n3 1\n", 0, 1);
  run_body("0 0 0 1\n0 1 1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n", 0, 1);
  run_body("0 0 0 0\n0 1 1 7\n1 0 2 -9\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n", 0, 1);
  run_body("0 0 0 1\n0 1 1 -1\n1 0 2 1\n1 1 3 -1\n0 1\n1 1\n2 1\n3 1\n", 1, 1);
  run_body("1 1 3 1\n0 1 1 1\n1 0 2 1\n0 0 0 1\n0 1\n1 1\n2 1\n3 1\n", 1, 1);
  run_body("0 0 0 1\n0 1 1 1\n1 0 2 1\n1 1 3 1\n3 1\n1 0\n0 1\n2 0\n", 1, 1);
  run_body("0 0 0 1\n0 1 1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n\n  \n\n",
           1, 1);
  run_body("0 0 0 1\n0 1 1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1", 1, 1);
  /* Leading and trailing whitespace and CRLF are accepted. */
  run_body("  0 0 0 1 \r\n0 1 1 1\r\n1 0 2 1\r\n1 1 3 1\r\n0 1\r\n1 1\r\n"
           "2 1\r\n3 1\r\n",
           1, 1);
  /* A shared parameter index used by two pairs is valid. */
  run_body("0 0 0 1\n0 1 0 -1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n", 1, 1);
  /* Real orbitals leave the imaginary OptFlag cleared. */
  {
    FILE *fp = file_from_text(
        "0 0 0 1\n0 1 1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n");
    reset_tables(&tables);
    CHECK(GCAntiReadOrbitals(fp, tables.idx, tables.sgn, tables.opt,
                             &tables.optCount, FIDX, 0, 1, NSITE, NORB,
                             "orbitalidx.def") == 0,
          "real orbital body");
    CHECK(tables.opt[2 * FIDX] == 1 && tables.opt[2 * FIDX + 1] == 0,
          "real orbital imaginary flag");
    fclose(fp);
  }
}

static void pair_negative_cases(void) {
  /* Plan corpus: AP requires the fourth column. */
  run_body_error("0 0 0\n0 1 1\n1 0 2\n1 1 3\n0 1\n1 1\n2 1\n3 1\n", 1,
                 "4 integer columns");
  run_body_error("0 0 4 1\n", 1, "parameter index");
  run_body_error("0 0 0 0\n", 1, "sign must be +1 or -1");
  run_body_error("-1 0 0 1\n", 1, "site index");
  run_body_error("0 0 9223372036854775808 1\n", 1, "malformed");
  /* Missing / surplus / duplicate rows. */
  run_body_error("0 0 0 1\n0 1 1 1\n1 0 2 1\n0 1\n1 1\n2 1\n3 1\n", 1,
                 "columns");
  run_body_error("0 0 0 1\n0 1 1 1\n1 0 2 1\n", 1, "ended after 3 of 4 pair");
  run_body_error(
      "0 0 0 1\n0 1 1 1\n1 0 2 1\n1 1 3 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n", 1,
      "columns");
  run_body_error("0 0 0 1\n0 1 1 1\n0 1 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n", 1,
                 "duplicate pair");
  run_body_error("0 0 0 1\n0 2 1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n", 1,
                 "site index");
  run_body_error("0 0 0 1\n0 1 -1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n",
                 1, "parameter index");
  run_body_error("0 0\n0 1 1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n", 0,
                 "columns");
  run_body_error("0 0 0 1 1\n0 1 1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n",
                 0, "columns");
  run_body_error("0 0 0 2\n0 1 1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n", 1,
                 "sign must be +1 or -1");
  run_body_error("0 0 0 1\n\n0 1 1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n",
                 1, "blank line");
  run_body_error("# comment\n0 0 0 1\n0 1 1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n"
                 "2 1\n3 1\n",
                 1, "malformed");
  /* PBC ignores the fourth column only when it is a valid int. */
  run_body_error("0 0 0 1.0\n0 1 1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n",
                 0, "malformed");
  run_body_error(
      "0 0 0 2147483648\n0 1 1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n", 0,
      "int range");
  run_body_error("0 0 0x1 1\n0 1 1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n",
                 1, "malformed");
  {
    char longLine[6000];
    memset(longLine, ' ', sizeof(longLine));
    memcpy(longLine, "0 0 0 1", 7);
    longLine[sizeof(longLine) - 2] = '\n';
    longLine[sizeof(longLine) - 1] = '\0';
    run_body_error(longLine, 1, "too long");
  }
}

static void opt_negative_cases(void) {
  const char *pairs = "0 0 0 1\n0 1 1 1\n1 0 2 1\n1 1 3 1\n";
  char text[512];
  const char *cases[][2] = {
      {"0 1\n1 1\n2 1\n", "ended after 3 of 4 OptFlag"},
      {"0 1\n1 1\n2 1\n3 1\n0 1\n", "after the OptFlag rows"},
      {"0 1\n1 1\n1 1\n3 1\n", "duplicate OptFlag"},
      {"-1 1\n1 1\n2 1\n3 1\n", "OptFlag index"},
      {"4 1\n1 1\n2 1\n3 1\n", "OptFlag index"},
      {"0 -1\n1 1\n2 1\n3 1\n", "flag must be 0 or 1"},
      {"0 2\n1 1\n2 1\n3 1\n", "flag must be 0 or 1"},
      {"0 1.0\n1 1\n2 1\n3 1\n", "malformed"},
      {"0 1 1\n1 1\n2 1\n3 1\n", "2 integer columns"},
      {"0 1\n\n1 1\n2 1\n3 1\n", "blank line"},
      {"0\n1 1\n2 1\n3 1\n", "2 integer columns"},
  };
  size_t k;
  for (k = 0; k < sizeof(cases) / sizeof(cases[0]); k++) {
    snprintf(text, sizeof(text), "%s%s", pairs, cases[k][0]);
    run_body_error(text, 1, cases[k][1]);
    run_body_error(text, 0, cases[k][1]);
  }
}

static void argument_cases(void) {
  Tables tables;
  FILE *fp;
  const char *body =
      "0 0 0 1\n0 1 1 1\n1 0 2 1\n1 1 3 1\n0 1\n1 1\n2 1\n3 1\n";
  struct {
    int fidx;
    int optCount;
    int complexFlag;
    int nsite;
    int norb;
    const char *message;
  } cases[] = {
      {-1, 0, 1, NSITE, NORB, "OptFlag offset"},
      {INT_MAX / 2, 0, 1, NSITE, NORB, "OptFlag offset"},
      {FIDX, INT_MAX - 1, 1, NSITE, NORB, "OptFlag count"},
      {FIDX, -1, 1, NSITE, NORB, "OptFlag count"},
      {FIDX, 0, 2, NSITE, NORB, "ComplexType"},
      {FIDX, 0, 1, 0, NORB, "Nsite"},
      {FIDX, 0, 1, -1, NORB, "Nsite"},
      {FIDX, 0, 1, INT_MIN, NORB, "Nsite"},
      {FIDX, 0, 1, 23171, NORB, "Nsite"},
      {FIDX, 0, 1, 46341, NORB, "Nsite"},
      {FIDX, 0, 1, INT_MAX, NORB, "Nsite"},
      {FIDX, 0, 1, NSITE, 0, "orbital count"},
  };
  size_t k;
  for (k = 0; k < sizeof(cases) / sizeof(cases[0]); k++) {
    fp = file_from_text(body);
    reset_tables(&tables);
    tables.optCount = cases[k].optCount;
    capture_begin();
    CHECK(GCAntiReadOrbitals(fp, tables.idx, tables.sgn, tables.opt,
                             &tables.optCount, cases[k].fidx,
                             cases[k].complexFlag, 1, cases[k].nsite,
                             cases[k].norb, "orbitalidx.def") != 0,
          "argument case %zu accepted", k);
    capture_end();
    CHECK(strstr(captured, cases[k].message) != NULL,
          "argument case %zu diagnostic '%s' missing: %s", k,
          cases[k].message, captured);
    CHECK(guards_intact(&tables), "argument case %zu wrote tables", k);
    {
      int i;
      for (i = 2 * FIDX; i < 2 * (FIDX + NORB); i++) {
        CHECK(tables.opt[i] == CANARY, "argument case %zu wrote OptFlag", k);
      }
    }
    fclose(fp);
  }
  capture_begin();
  CHECK(GCAntiReadOrbitals(NULL, tables.idx, tables.sgn, tables.opt,
                           &tables.optCount, FIDX, 1, 1, NSITE, NORB,
                           "orbitalidx.def") != 0,
        "NULL FILE accepted");
  capture_end();
}

static int run_header(const char *text, const int nsite, const int expectOk,
                      const char *expectedMessage, int *norb,
                      int *complexFlag) {
  FILE *fp = file_from_text(text);
  int result;
  *norb = -5;
  *complexFlag = -5;
  capture_begin();
  result = GCAntiReadHeader(fp, nsite, norb, complexFlag, "orbitalidx.def");
  capture_end();
  fclose(fp);
  if (expectOk) {
    CHECK(result == 0, "header expected success: %s\n%s", captured, text);
  } else {
    CHECK(result != 0, "header expected failure:\n%s", text);
    CHECK(strstr(captured, "GC anti-parallel") != NULL &&
              strstr(captured, "orbitalidx.def") != NULL,
          "header diagnostic lacks prefix/file: %s", captured);
    if (expectedMessage != NULL) {
      CHECK(strstr(captured, expectedMessage) != NULL,
            "header diagnostic '%s' missing: %s", expectedMessage, captured);
    }
    CHECK(*norb == -5 && *complexFlag == -5,
          "header wrote outputs on failure");
  }
  return result;
}

static void header_cases(void) {
  int norb;
  int complexFlag;
  char text[512];
  const char *bar = "=============================================\n";
  struct {
    const char *line2;
    const char *line3;
    const char *message;
  } bad[] = {
      {"NOrbitalIdx 0\n", "ComplexType 1\n", "orbital count"},
      {"NOrbitalIdx -3\n", "ComplexType 1\n", "orbital count"},
      {"NOrbitalIdx 2147483648\n", "ComplexType 1\n", "orbital count"},
      {"NOrbitalIdx 9223372036854775808\n", "ComplexType 1\n", "malformed"},
      {"NOrbitalIdx 4.0\n", "ComplexType 1\n", "malformed"},
      {"NOrbitalIdx 1073741824\n", "ComplexType 1\n", "orbital count"},
      {"NOrbitalIdx 4\n", "ComplexType 2\n", "ComplexType"},
      {"NOrbitalIdx 4\n", "ComplexType -1\n", "ComplexType"},
      {"NOrbitalIdx 4 5\n", "ComplexType 1\n", "label and one integer"},
      {"4\n", "ComplexType 1\n", "label and one integer"},
      {"NOrbitalIdx 4\n", "ComplexType\n", "label and one integer"},
  };
  size_t k;
  snprintf(text, sizeof(text), "%sNOrbitalIdx 4\nComplexType 1\n%s%s", bar,
           bar, bar);
  CHECK(run_header(text, NSITE, 1, NULL, &norb, &complexFlag) == 0 &&
            norb == 4 && complexFlag == 1,
        "valid header values");
  snprintf(text, sizeof(text), "%s  Label   16  \nAnyName 0\n%s%s", bar, bar,
           bar);
  CHECK(run_header(text, 4, 1, NULL, &norb, &complexFlag) == 0 && norb == 16 &&
            complexFlag == 0,
        "header label names are not inspected");
  for (k = 0; k < sizeof(bad) / sizeof(bad[0]); k++) {
    snprintf(text, sizeof(text), "%s%s%s%s%s", bar, bad[k].line2, bad[k].line3,
             bar, bar);
    run_header(text, NSITE, 0, bad[k].message, &norb, &complexFlag);
  }
  snprintf(text, sizeof(text), "%sNOrbitalIdx 4\nComplexType 1\n%s", bar, bar);
  run_header(text, NSITE, 0, "header", &norb, &complexFlag);
  run_header("", NSITE, 0, "header", &norb, &complexFlag);
  snprintf(text, sizeof(text), "%sNOrbitalIdx 4\nComplexType 1\n%s%s", bar,
           bar, bar);
  run_header(text, 0, 0, "Nsite", &norb, &complexFlag);
  run_header(text, -2, 0, "Nsite", &norb, &complexFlag);
  run_header(text, INT_MAX / 2 + 1, 0, "Nsite", &norb, &complexFlag);
  /* (2*Nsite)^2 must fit in int for the Slater element stride. */
  run_header(text, 23171, 0, "Nsite", &norb, &complexFlag);
  CHECK(run_header(text, 23170, 1, NULL, &norb, &complexFlag) == 0,
        "largest supported Nsite");
}

int main(void) {
  positive_cases();
  pair_negative_cases();
  opt_negative_cases();
  argument_cases();
  header_cases();
  if (failures != 0) {
    fprintf(stderr, "GC_AntiPair_Input_Unit: %d failure(s)\n", failures);
    return EXIT_FAILURE;
  }
  printf("GC_AntiPair_Input_Unit passed\n");
  return EXIT_SUCCESS;
}
