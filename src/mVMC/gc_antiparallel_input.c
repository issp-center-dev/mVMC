#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "include/gc_antiparallel_input.h"

#define ANTI_LINE_CAPACITY 4096
#define ANTI_HEADER_LINES 5

enum {
  ANTI_LINE_OK = 0,
  ANTI_LINE_EOF = 1,
  ANTI_LINE_ERROR = 2,
  ANTI_LINE_TOO_LONG = 3
};

/* 0=success; nonzero for EOF (1), read error (2) or an overlong line (3).
 * A short final line without a newline is accepted. */
static int AntiReadLine(FILE *fp, char line[ANTI_LINE_CAPACITY]) {
  if (fgets(line, ANTI_LINE_CAPACITY, fp) == NULL) {
    return ferror(fp) ? ANTI_LINE_ERROR : ANTI_LINE_EOF;
  }
  if (strchr(line, '\n') == NULL && !feof(fp)) return ANTI_LINE_TOO_LONG;
  return ANTI_LINE_OK;
}

/* Return the number of whitespace-separated tokens, storing the first
 * `capacity` values, or -1 when a token is not a whole base-10 integer in
 * long long range.  Surplus tokens are counted so the caller can report the
 * column count.  A blank line has zero tokens. */
static int AntiParseIntegers(const char *line, long long *value,
                             const int capacity) {
  const char *p = line;
  int count = 0;
  while (*p != '\0') {
    char *end = NULL;
    long long v;
    while (isspace((unsigned char)*p)) ++p;
    if (*p == '\0') break;
    errno = 0;
    v = strtoll(p, &end, 10);
    if (errno == ERANGE || end == p ||
        (*end != '\0' && !isspace((unsigned char)*end))) {
      return -1;
    }
    if (count < capacity) value[count] = v;
    count++;
    p = end;
  }
  return count;
}

/* One label token followed by one integer: 0 success, 1 wrong token count,
 * 2 malformed or out-of-range integer. */
static int AntiHeaderValue(const char *line, long long *value) {
  const char *p = line;
  const char *token[3] = {NULL, NULL, NULL};
  size_t length[3] = {0, 0, 0};
  int count = 0;
  char buffer[ANTI_LINE_CAPACITY];
  char *end = NULL;
  long long v;
  while (*p != '\0') {
    while (isspace((unsigned char)*p)) ++p;
    if (*p == '\0') break;
    if (count < 3) token[count] = p;
    while (*p != '\0' && !isspace((unsigned char)*p)) ++p;
    if (count < 3) length[count] = (size_t)(p - token[count]);
    count++;
  }
  if (count != 2) return 1;
  memcpy(buffer, token[1], length[1]);
  buffer[length[1]] = '\0';
  errno = 0;
  v = strtoll(buffer, &end, 10);
  if (errno == ERANGE || end == buffer || *end != '\0') return 2;
  *value = v;
  return 0;
}

static int AntiNsiteSupported(const int nsite) {
  long long side;
  if (nsite <= 0 || nsite > INT_MAX / 2) return 0;
  side = 2LL * (long long)nsite;
  return side * side <= (long long)INT_MAX;
}

static int AntiLineError(const char *filename, const int lineNumber,
                         const int status) {
  if (status == ANTI_LINE_TOO_LONG) {
    fprintf(stderr,
            "Error: GC anti-parallel orbital file %s line %d is too long "
            "(limit %d bytes).\n",
            filename, lineNumber, ANTI_LINE_CAPACITY - 1);
  } else if (status == ANTI_LINE_ERROR) {
    fprintf(stderr,
            "Error: GC anti-parallel orbital file %s: read error at line %d.\n",
            filename, lineNumber);
  }
  return 1;
}

int GCAntiReadHeader(FILE *fp, const int nsite, int *norb, int *complexFlag,
                     const char *filename) {
  char line[ANTI_LINE_CAPACITY];
  long long norbValue = 0;
  long long complexValue = 0;
  int row;
  if (filename == NULL) filename = "(unnamed)";
  if (fp == NULL || norb == NULL || complexFlag == NULL) {
    fprintf(stderr, "Error: GC anti-parallel orbital file %s: invalid "
                    "header reader arguments.\n", filename);
    return 1;
  }
  if (!AntiNsiteSupported(nsite)) {
    fprintf(stderr,
            "Error: GC anti-parallel orbital file %s: Nsite=%d is "
            "unsupported (need 0 < Nsite and (2*Nsite)^2 <= %d).\n",
            filename, nsite, INT_MAX);
    return 1;
  }
  for (row = 1; row <= ANTI_HEADER_LINES; row++) {
    const int status = AntiReadLine(fp, line);
    if (status == ANTI_LINE_EOF) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s: the header ends at "
              "line %d (5 header lines are required).\n",
              filename, row);
      return 1;
    }
    if (status != ANTI_LINE_OK) return AntiLineError(filename, row, status);
    if (row == 2 || row == 3) {
      long long *target = row == 2 ? &norbValue : &complexValue;
      const int parsed = AntiHeaderValue(line, target);
      if (parsed == 1) {
        fprintf(stderr,
                "Error: GC anti-parallel orbital file %s line %d must contain "
                "a label and one integer: %s",
                filename, row, line);
        return 1;
      }
      if (parsed != 0) {
        fprintf(stderr,
                "Error: GC anti-parallel orbital file %s line %d has a "
                "malformed or out-of-range integer: %s",
                filename, row, line);
        return 1;
      }
    }
  }
  if (norbValue <= 0 || norbValue > INT_MAX / 2) {
    fprintf(stderr,
            "Error: GC anti-parallel orbital file %s line 2: the orbital count "
            "must be in [1, %d] (got %lld).\n",
            filename, INT_MAX / 2, norbValue);
    return 1;
  }
  if (complexValue != 0 && complexValue != 1) {
    fprintf(stderr,
            "Error: GC anti-parallel orbital file %s line 3: ComplexType must "
            "be 0 or 1 (got %lld).\n",
            filename, complexValue);
    return 1;
  }
  *norb = (int)norbValue;
  *complexFlag = (int)complexValue;
  return 0;
}

int GCAntiReadOrbitals(FILE *fp, int **idx, int **sgn, int *opt,
                       int *optCount, const int fidx, const int complexFlag,
                       const int apFlag, const int nsite, const int norb,
                       const char *filename) {
  char line[ANTI_LINE_CAPACITY];
  unsigned char *seenPair = NULL;
  unsigned char *seenOpt = NULL;
  int pairRows;
  int lineNumber = ANTI_HEADER_LINES;
  int row;
  int status;
  int info = 0;
  if (filename == NULL) filename = "(unnamed)";
  if (fp == NULL || idx == NULL || sgn == NULL || opt == NULL ||
      optCount == NULL) {
    fprintf(stderr, "Error: GC anti-parallel orbital file %s: invalid body "
                    "reader arguments.\n", filename);
    return 1;
  }
  if (!AntiNsiteSupported(nsite)) {
    fprintf(stderr,
            "Error: GC anti-parallel orbital file %s: Nsite=%d is "
            "unsupported.\n",
            filename, nsite);
    return 1;
  }
  pairRows = nsite * nsite;
  if (norb <= 0 || norb > INT_MAX / 2) {
    fprintf(stderr,
            "Error: GC anti-parallel orbital file %s: the orbital count must "
            "be in [1, %d] (got %d).\n",
            filename, INT_MAX / 2, norb);
    return 1;
  }
  /* 2*(fidx+index)+1 must stay inside int for every index < norb. */
  if (fidx < 0 || fidx > INT_MAX / 2 - norb) {
    fprintf(stderr,
            "Error: GC anti-parallel orbital file %s: OptFlag offset %d is "
            "out of range for %d orbitals.\n",
            filename, fidx, norb);
    return 1;
  }
  if (*optCount < 0 || *optCount > INT_MAX - norb) {
    fprintf(stderr,
            "Error: GC anti-parallel orbital file %s: OptFlag count %d would "
            "overflow.\n",
            filename, *optCount);
    return 1;
  }
  if (complexFlag != 0 && complexFlag != 1) {
    fprintf(stderr,
            "Error: GC anti-parallel orbital file %s: ComplexType must be 0 "
            "or 1 (got %d).\n",
            filename, complexFlag);
    return 1;
  }
  seenPair = (unsigned char *)calloc((size_t)pairRows, sizeof(*seenPair));
  seenOpt = (unsigned char *)calloc((size_t)norb, sizeof(*seenOpt));
  if (seenPair == NULL || seenOpt == NULL) {
    fprintf(stderr,
            "Error: GC anti-parallel orbital file %s: failed to allocate "
            "reader state.\n",
            filename);
    free(seenPair);
    free(seenOpt);
    return 1;
  }

  for (row = 0; row < pairRows && info == 0; row++) {
    long long v[4];
    int ncol;
    size_t key;
    lineNumber++;
    status = AntiReadLine(fp, line);
    if (status == ANTI_LINE_EOF) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s ended after %d of %d "
              "pair rows.\n",
              filename, row, pairRows);
      info = 1;
      break;
    }
    if (status != ANTI_LINE_OK) {
      info = AntiLineError(filename, lineNumber, status);
      break;
    }
    ncol = AntiParseIntegers(line, v, 4);
    if (ncol < 0) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s pair row %d (line %d) "
              "has a malformed or out-of-range integer token: %s",
              filename, row + 1, lineNumber, line);
      info = 1;
    } else if (ncol == 0) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s: blank line inside the "
              "pair block (line %d).\n",
              filename, lineNumber);
      info = 1;
    } else if ((apFlag && ncol != 4) || (!apFlag && ncol != 3 && ncol != 4)) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s pair row %d (line %d) "
              "has %d columns; %s: %s",
              filename, row + 1, lineNumber, ncol,
              apFlag ? "anti-periodic input requires 4 integer columns "
                       "(i j index sign)"
                     : "periodic input requires 3 or 4 integer columns",
              line);
      info = 1;
    } else if (v[0] < 0 || v[0] >= nsite || v[1] < 0 || v[1] >= nsite) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s pair row %d (line %d) "
              "has a site index out of range [0, %d): %s",
              filename, row + 1, lineNumber, nsite, line);
      info = 1;
    } else if (v[2] < 0 || v[2] >= norb) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s pair row %d (line %d) "
              "has a parameter index out of range [0, %d): %s",
              filename, row + 1, lineNumber, norb, line);
      info = 1;
    } else if (ncol == 4 && (v[3] < INT_MIN || v[3] > INT_MAX)) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s pair row %d (line %d) "
              "has a fourth column outside the int range: %s",
              filename, row + 1, lineNumber, line);
      info = 1;
    } else if (apFlag && v[3] != -1 && v[3] != 1) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s pair row %d (line %d): "
              "the anti-periodic sign must be +1 or -1 (got %lld).\n",
              filename, row + 1, lineNumber, v[3]);
      info = 1;
    }
    if (info != 0) break;
    key = (size_t)v[0] * (size_t)nsite + (size_t)v[1];
    if (seenPair[key]) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s: duplicate pair (%lld,"
              "%lld) at line %d.\n",
              filename, v[0], v[1], lineNumber);
      info = 1;
      break;
    }
    seenPair[key] = 1;
    idx[v[0]][v[1]] = (int)v[2];
    sgn[v[0]][v[1]] = apFlag ? (int)v[3] : 1;
  }

  for (row = 0; row < norb && info == 0; row++) {
    long long v[2];
    int ncol;
    lineNumber++;
    status = AntiReadLine(fp, line);
    if (status == ANTI_LINE_EOF) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s ended after %d of %d "
              "OptFlag rows.\n",
              filename, row, norb);
      info = 1;
      break;
    }
    if (status != ANTI_LINE_OK) {
      info = AntiLineError(filename, lineNumber, status);
      break;
    }
    ncol = AntiParseIntegers(line, v, 2);
    if (ncol < 0) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s OptFlag row %d (line "
              "%d) has a malformed or out-of-range integer token: %s",
              filename, row + 1, lineNumber, line);
      info = 1;
    } else if (ncol == 0) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s: blank line inside the "
              "OptFlag block (line %d).\n",
              filename, lineNumber);
      info = 1;
    } else if (ncol != 2) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s OptFlag row %d (line "
              "%d) has %d columns; exactly 2 integer columns (index flag) are "
              "required: %s",
              filename, row + 1, lineNumber, ncol, line);
      info = 1;
    } else if (v[0] < 0 || v[0] >= norb) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s OptFlag row %d (line "
              "%d): OptFlag index %lld is out of range [0, %d).\n",
              filename, row + 1, lineNumber, v[0], norb);
      info = 1;
    } else if (v[1] != 0 && v[1] != 1) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s OptFlag row %d (line "
              "%d): the flag must be 0 or 1 (got %lld).\n",
              filename, row + 1, lineNumber, v[1]);
      info = 1;
    } else if (seenOpt[v[0]]) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s: duplicate OptFlag "
              "index %lld at line %d.\n",
              filename, v[0], lineNumber);
      info = 1;
    }
    if (info != 0) break;
    seenOpt[v[0]] = 1;
    opt[2 * (fidx + (int)v[0])] = (int)v[1];
    opt[2 * (fidx + (int)v[0]) + 1] = complexFlag ? (int)v[1] : 0;
    (*optCount)++;
  }

  while (info == 0) {
    long long unused[1];
    lineNumber++;
    status = AntiReadLine(fp, line);
    if (status == ANTI_LINE_EOF) break;
    if (status != ANTI_LINE_OK) {
      info = AntiLineError(filename, lineNumber, status);
      break;
    }
    if (AntiParseIntegers(line, unused, 1) != 0) {
      fprintf(stderr,
              "Error: GC anti-parallel orbital file %s: unexpected content "
              "after the OptFlag rows (line %d): %s",
              filename, lineNumber, line);
      info = 1;
    }
  }
  free(seenPair);
  free(seenOpt);
  return info;
}
