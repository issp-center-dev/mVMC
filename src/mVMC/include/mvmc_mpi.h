/*
mVMC - A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is developed based on the mVMC-mini program
(https://github.com/fiber-miniapp/mVMC-mini)
which follows "The BSD 3-Clause License".

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.
*/
/* MPI entry points used by both MPI and serial builds. */
#ifndef MVMC_MPI_H
#define MVMC_MPI_H

#ifdef _mpi_use
#include <mpi.h>
#else
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <sys/time.h>

typedef int MPI_Comm;
#define MPI_COMM_WORLD 0
#define MPI_COMM_SELF 1
#define MPI_COMM_NULL (-1)
#define MPI_SUCCESS 0
#define MPI_IN_PLACE ((void *)1)
typedef enum {
  MPI_INT, MPI_UNSIGNED, MPI_LONG_LONG, MPI_DOUBLE, MPI_DOUBLE_COMPLEX
} MPI_Datatype;
typedef enum { MPI_SUM, MPI_MIN, MPI_MAX } MPI_Op;

static inline int MPI_Init(int *argc, char ***argv) {
  (void)argc;
  (void)argv;
  return 0;
}
static inline int MPI_Finalize(void) { return 0; }
static inline double MPI_Wtime(void) {
  struct timeval now;
  gettimeofday(&now, NULL);
  return (double)now.tv_sec + 1.0e-6 * (double)now.tv_usec;
}
static inline int MPI_Abort(MPI_Comm comm, int errorcode) {
  (void)comm;
  exit(errorcode);
}
static inline int MPI_Barrier(MPI_Comm comm) {
  (void)comm;
  return 0;
}
static inline int MPI_Comm_size(MPI_Comm comm, int *size) {
  (void)comm;
  *size = 1;
  return 0;
}
static inline int MPI_Comm_rank(MPI_Comm comm, int *rank) {
  (void)comm;
  *rank = 0;
  return 0;
}
static inline int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm) {
  *newcomm = comm;
  return 0;
}
static inline int MPI_Comm_split(MPI_Comm comm, int color, int key,
                                 MPI_Comm *newcomm) {
  (void)color;
  (void)key;
  *newcomm = comm;
  return 0;
}
static inline int MPI_Comm_free(MPI_Comm *comm) {
  *comm = MPI_COMM_NULL;
  return 0;
}
/* Only the datatypes used by the serial production paths are supported. */
static inline size_t mvmc_mpi_type_size(MPI_Datatype datatype) {
  switch (datatype) {
    case MPI_INT: return sizeof(int);
    case MPI_UNSIGNED: return sizeof(unsigned int);
    case MPI_LONG_LONG: return sizeof(long long);
    case MPI_DOUBLE: return sizeof(double);
    case MPI_DOUBLE_COMPLEX: return sizeof(double complex);
  }
  return 0;
}
static inline int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype,
                           int root, MPI_Comm comm) {
  (void)buffer;
  (void)comm;
  return count >= 0 && root == 0 && mvmc_mpi_type_size(datatype) != 0 ? 0 : 1;
}
static inline int MPI_Allreduce(const void *send, void *recv, int count,
                               MPI_Datatype datatype, MPI_Op op,
                               MPI_Comm comm) {
  const size_t width = mvmc_mpi_type_size(datatype);
  (void)comm;
  if (count < 0 || width == 0 || (size_t)count > (size_t)-1 / width ||
      (op != MPI_SUM && op != MPI_MIN && op != MPI_MAX)) return 1;
  if (count > 0 && send != MPI_IN_PLACE && send != recv)
    memcpy(recv, send, (size_t)count * width);
  return MPI_SUCCESS;
}
#endif /* _mpi_use */
#endif /* MVMC_MPI_H */
