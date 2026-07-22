#ifndef _LSLOCGRN_STATUS_H
#define _LSLOCGRN_STATUS_H

#include <errno.h>
#include <stdlib.h>

typedef enum {
  LSLANCZOS_OK = 0,
  LSLANCZOS_STATE_FAILURE = 1,
  LSLANCZOS_REAL_GATE_FAILURE = 2,
  LSLANCZOS_NUMERIC_REJECT = 3
} LSLanczosStatus;

/* Test-only fault injection.  A positive value rejects this many candidate
 * rebuilds per process without changing production behavior when unset. */
static inline int LSLanczosTestConsumeRebuildFailure(void) {
  static int initialized = 0;
  static long remaining = 0;
  if(!initialized) {
    const char *value = getenv("MVMC_LANCZOS_TEST_REBUILD_FAILURES");
    char *end = NULL;
    initialized = 1;
    if(value != NULL && value[0] != '\0') {
      long parsed;
      errno = 0;
      parsed = strtol(value, &end, 10);
      if(errno == 0 && end != value && *end == '\0' && parsed > 0) {
        remaining = parsed;
      }
    }
  }
  if(remaining <= 0) return 0;
  remaining--;
  return 1;
}

static inline int LSLanczosTestProjectionBranchAuditEnabled(void) {
  const char *value = getenv("MVMC_LANCZOS_TEST_PROJECTION_BRANCH_AUDIT");
  return value != NULL && value[0] != '\0' && value[0] != '0';
}

/* Optional parent-communicator rank selector for non-finite injection.
 * Returns -1 when unset and -2 for malformed input. */
static inline long LSLanczosTestNonfiniteParentRank(void) {
  const char *value = getenv("MVMC_LANCZOS_TEST_NONFINITE_PARENT_RANK");
  char *end = NULL;
  long rank;
  if(value == NULL || value[0] == '\0') return -1;
  errno = 0;
  rank = strtol(value, &end, 10);
  if(errno != 0 || end == value || *end != '\0' || rank < 0) return -2;
  return rank;
}

#endif
