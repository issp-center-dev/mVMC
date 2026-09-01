#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _mpi_use
#include <mpi.h>
#endif

#include "include/gc_size.h"

static void GCSizeAbort(void) {
#ifdef _mpi_use
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized != 0) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  exit(EXIT_FAILURE);
}

size_t GCCheckedMulSize(const size_t a, const size_t b) {
  if (a != 0 && b > SIZE_MAX / a) {
    fprintf(stderr, "Error: allocation size overflow (mul %zu x %zu).\n", a,
            b);
    GCSizeAbort();
  }
  return a * b;
}

size_t GCCheckedAddSize(const size_t a, const size_t b) {
  if (b > SIZE_MAX - a) {
    fprintf(stderr, "Error: allocation size overflow (add %zu + %zu).\n", a,
            b);
    GCSizeAbort();
  }
  return a + b;
}

size_t GCNonnegativeSize(const long value, const char *name) {
  if (value < 0) {
    fprintf(stderr, "Error: %s must be non-negative (got %ld).\n", name,
            value);
    GCSizeAbort();
  }
  return (size_t)value;
}

void *GCCheckedMalloc(size_t nbytes, const char *name) {
  void *allocation;
  if (nbytes == 0) nbytes = 1;
  allocation = malloc(nbytes);
  if (allocation == NULL) {
    fprintf(stderr, "Error: allocation failed for %s (%zu bytes).\n", name,
            nbytes);
    GCSizeAbort();
  }
  return allocation;
}

int GCCheckedSizeToInt(const size_t value, const char *name) {
  if (value > (size_t)INT_MAX) {
    fprintf(stderr, "Error: %s exceeds workspace API INT_MAX (%zu).\n", name,
            value);
    GCSizeAbort();
  }
  return (int)value;
}
