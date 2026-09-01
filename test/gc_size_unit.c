#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gc_size.h"

static int check_boundaries(void) {
  void *zeroAllocation;

  if (GCCheckedMulSize(0, SIZE_MAX) != 0) return 1;
  if (GCCheckedMulSize(1, SIZE_MAX) != SIZE_MAX) return 1;
  if (GCCheckedAddSize(SIZE_MAX - 1, 1) != SIZE_MAX) return 1;
  if (GCNonnegativeSize(0, "zero") != 0) return 1;
  if (GCNonnegativeSize(LONG_MAX, "long max") != (size_t)LONG_MAX) {
    return 1;
  }
  if (GCCheckedSizeToInt((size_t)INT_MAX, "int max") != INT_MAX) return 1;

  zeroAllocation = GCCheckedMalloc(0, "logical zero allocation");
  if (zeroAllocation == NULL) return 1;
  free(zeroAllocation);
  return 0;
}

int main(int argc, char **argv) {
  if (argc == 1) return check_boundaries();
  if (strcmp(argv[1], "mul_overflow") == 0) {
    (void)GCCheckedMulSize(SIZE_MAX, 2);
  } else if (strcmp(argv[1], "add_overflow") == 0) {
    (void)GCCheckedAddSize(SIZE_MAX, 1);
  } else if (strcmp(argv[1], "negative") == 0) {
    (void)GCNonnegativeSize(-1, "negative test value");
  } else if (strcmp(argv[1], "int_overflow") == 0) {
    (void)GCCheckedSizeToInt((size_t)INT_MAX + 1, "int overflow");
  } else {
    fprintf(stderr, "unknown gc_size_unit mode: %s\n", argv[1]);
    return 2;
  }
  fprintf(stderr, "expected gc_size helper to terminate\n");
  return 1;
}
