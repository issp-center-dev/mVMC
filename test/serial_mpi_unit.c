/* Exercise the serial communication contract in both build configurations. */
#undef _mpi_use
#include "mvmc_mpi.h"
#include "mvmc_mpi.h"
#include <stdio.h>
#include "../src/mVMC/safempi.c"
#include "../src/mVMC/safempi_fcmp.c"

#define CHECK(expr) do { if (!(expr)) { \
  fprintf(stderr, "serial MPI check failed at line %d: %s\n", __LINE__, #expr); \
  return 1; } } while (0)

int main(int argc, char **argv) {
  MPI_Comm comm, split;
  int rank = -1, size = -1;
  const int ints[3] = {2, -4, 7};
  int received[3] = {0};
  double values[3] = {0.25, -3.5, 8.0}, real_out[3] = {0};
  double complex complex_values[2] = {1.0 + 2.0 * I, -3.0 + 0.5 * I};
  double complex complex_out[2] = {0};
  int i;
  CHECK(MPI_Init(&argc, &argv) == MPI_SUCCESS);
  CHECK(MPI_Comm_dup(MPI_COMM_WORLD, &comm) == MPI_SUCCESS);
  CHECK(MPI_Comm_split(comm, 0, 0, &split) == MPI_SUCCESS);
  CHECK(MPI_Comm_rank(split, &rank) == MPI_SUCCESS && rank == 0);
  CHECK(MPI_Comm_size(split, &size) == MPI_SUCCESS && size == 1);
  CHECK(MPI_Barrier(split) == MPI_SUCCESS);
  CHECK(MPI_Wtime() > 0.0);
  CHECK(MPI_Bcast(values, 3, MPI_DOUBLE, 0, comm) == MPI_SUCCESS);
  CHECK(values[0] == 0.25 && values[1] == -3.5 && values[2] == 8.0);
  CHECK(MPI_Allreduce(ints, received, 3, MPI_INT, MPI_MAX, comm) == MPI_SUCCESS);
  for (i = 0; i < 3; ++i) CHECK(received[i] == ints[i]);
  CHECK(MPI_Allreduce(values, real_out, 3, MPI_DOUBLE, MPI_MIN, comm) == MPI_SUCCESS);
  for (i = 0; i < 3; ++i) CHECK(real_out[i] == values[i]);
  CHECK(MPI_Allreduce(complex_values, complex_out, 2, MPI_DOUBLE_COMPLEX,
                      MPI_SUM, MPI_COMM_SELF) == MPI_SUCCESS);
  for (i = 0; i < 2; ++i) CHECK(complex_out[i] == complex_values[i]);
  CHECK(MPI_Allreduce(MPI_IN_PLACE, complex_out, 2, MPI_DOUBLE_COMPLEX,
                      MPI_SUM, comm) == MPI_SUCCESS);
  for (i = 0; i < 2; ++i) CHECK(complex_out[i] == complex_values[i]);
  CHECK(MPI_Allreduce(NULL, NULL, 0, MPI_INT, MPI_SUM, comm) == MPI_SUCCESS);
  CHECK(MPI_Allreduce(ints, received, -1, MPI_INT, MPI_SUM, comm) != MPI_SUCCESS);
  CHECK(MPI_Allreduce(ints, received, 3, (MPI_Datatype)99, MPI_SUM, comm) != MPI_SUCCESS);
  CHECK(MPI_Bcast(values, 3, MPI_DOUBLE, 1, comm) != MPI_SUCCESS);

  /* SR-CG calls SafeMpiAllReduce even when the communicator has size 1. */
  memset(real_out, 0, sizeof(real_out));
  SafeMpiAllReduce(values, real_out, 3, comm);
  for (i = 0; i < 3; ++i) CHECK(real_out[i] == values[i]);
  memset(real_out, 0, sizeof(real_out));
  SafeMpiReduce(values, real_out, 3, comm);
  for (i = 0; i < 3; ++i) CHECK(real_out[i] == values[i]);
  memset(complex_out, 0, sizeof(complex_out));
  SafeMpiAllReduce_fcmp(complex_values, complex_out, 2, comm);
  for (i = 0; i < 2; ++i) CHECK(complex_out[i] == complex_values[i]);
  memset(complex_out, 0, sizeof(complex_out));
  SafeMpiReduce_fcmp(complex_values, complex_out, 2, comm);
  for (i = 0; i < 2; ++i) CHECK(complex_out[i] == complex_values[i]);
  SafeMpiAllReduce(real_out, real_out, 3, comm);
  SafeMpiAllReduce_fcmp(complex_out, complex_out, 2, comm);
  SafeMpiReduce(NULL, NULL, 0, comm);
  SafeMpiReduce_fcmp(NULL, NULL, 0, comm);
  for (i = 0; i < 3; ++i) CHECK(real_out[i] == values[i]);
  for (i = 0; i < 2; ++i) CHECK(complex_out[i] == complex_values[i]);
  CHECK(MPI_Comm_free(&split) == MPI_SUCCESS && split == MPI_COMM_NULL);
  CHECK(MPI_Comm_free(&comm) == MPI_SUCCESS && comm == MPI_COMM_NULL);
  CHECK(MPI_Finalize() == MPI_SUCCESS);
  puts("serial_mpi_unit: PASS");
  return 0;
}
