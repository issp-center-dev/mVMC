/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#include "power_lanczos_sampling.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "power_lanczos_sampling.c requires the power-Lanczos core"
#endif

#include <limits.h>
#include <math.h>
#include <string.h>

static int valid_status(MVMCKrylovStatus status) {
  return status >= MVMC_KRYLOV_STATUS_OK &&
         status <= MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE;
}

MVMCKrylovStatus mvmc_power_lanczos_sampling_synchronize(
    const MVMCPowerLanczosSamplingTopology *topology,
    MVMCKrylovStatus local_status) {
  MVMCKrylovStatus effective =
      valid_status(local_status) ? local_status
                                 : MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
#ifdef _mpi_use
  int encoded;
  int reduced = 0;
  if (topology == NULL ||
      topology->world_communicator == MPI_COMM_NULL) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  encoded = (int)effective;
  if (MPI_Allreduce(&encoded, &reduced, 1, MPI_INT, MPI_MAX,
                    topology->world_communicator) != MPI_SUCCESS ||
      reduced < (int)MVMC_KRYLOV_STATUS_OK ||
      reduced > (int)MVMC_KRYLOV_STATUS_INTERNAL_INVARIANT_FAILURE) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
  return (MVMCKrylovStatus)reduced;
#else
  (void)topology;
  return effective;
#endif
}

MVMCKrylovStatus mvmc_power_lanczos_sampling_topology_create(
    MVMCClassicPfaffianCommunicator world_communicator,
    MVMCClassicPfaffianCommunicator chain_communicator,
    size_t samples_per_chain,
    MVMCPowerLanczosSamplingTopology *topology) {
  uint64_t sample_count;
  uint64_t total_samples;
  if (topology == NULL || samples_per_chain == 0) {
    return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
#if SIZE_MAX > UINT64_MAX
  if (samples_per_chain > UINT64_MAX) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
#endif
  sample_count = (uint64_t)samples_per_chain;
  memset(topology, 0, sizeof(*topology));
  topology->world_communicator = world_communicator;
  topology->chain_communicator = chain_communicator;
#ifdef _mpi_use
  {
    uint64_t root_flag;
    uint64_t root_prefix = 0;
    uint64_t minimum_samples = 0;
    uint64_t maximum_samples = 0;
    if (world_communicator == MPI_COMM_NULL ||
        chain_communicator == MPI_COMM_NULL ||
        MPI_Comm_rank(world_communicator, &topology->world_rank) !=
            MPI_SUCCESS ||
        MPI_Comm_size(world_communicator, &topology->world_size) !=
            MPI_SUCCESS ||
        MPI_Comm_rank(chain_communicator, &topology->chain_rank) !=
            MPI_SUCCESS ||
        MPI_Comm_size(chain_communicator, &topology->chain_size) !=
            MPI_SUCCESS ||
        topology->world_size <= 0 || topology->chain_size <= 0) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    root_flag = topology->chain_rank == 0 ? UINT64_C(1) : UINT64_C(0);
    if (MPI_Allreduce(&sample_count, &minimum_samples, 1, MPI_UINT64_T,
                      MPI_MIN, world_communicator) != MPI_SUCCESS ||
        MPI_Allreduce(&sample_count, &maximum_samples, 1, MPI_UINT64_T,
                      MPI_MAX, world_communicator) != MPI_SUCCESS ||
        MPI_Allreduce(&root_flag, &topology->chain_count, 1, MPI_UINT64_T,
                      MPI_SUM, world_communicator) != MPI_SUCCESS ||
        MPI_Exscan(&root_flag, &root_prefix, 1, MPI_UINT64_T, MPI_SUM,
                   world_communicator) != MPI_SUCCESS) {
      return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
    if (topology->world_rank == 0) root_prefix = 0;
    if (minimum_samples != maximum_samples || minimum_samples == 0 ||
        topology->chain_count == 0) {
      return MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
    }
    topology->chain_index =
        topology->chain_rank == 0 ? root_prefix : UINT64_C(0);
    if (MPI_Bcast(&topology->chain_index, 1, MPI_UINT64_T, 0,
                  chain_communicator) != MPI_SUCCESS) {
      return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
    }
  }
#else
  (void)world_communicator;
  (void)chain_communicator;
  topology->world_rank = 0;
  topology->world_size = 1;
  topology->chain_rank = 0;
  topology->chain_size = 1;
  topology->chain_index = 0;
  topology->chain_count = 1;
#endif
  if (topology->chain_index >= topology->chain_count ||
      topology->chain_count > UINT64_MAX / sample_count) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
  total_samples = topology->chain_count * sample_count;
#if UINT64_MAX > SIZE_MAX
  if (total_samples > SIZE_MAX) {
    return MVMC_KRYLOV_STATUS_RESOURCE_LIMIT;
  }
#endif
  topology->samples_per_chain = sample_count;
  topology->total_samples = total_samples;
  return MVMC_KRYLOV_STATUS_OK;
}

static MVMCKrylovStatus validate_reduction(
    const MVMCPowerLanczosSamplingTopology *topology,
    const void *values, size_t count, size_t item_size) {
  MVMCKrylovStatus status = MVMC_KRYLOV_STATUS_OK;
  if (topology == NULL || (count != 0 && values == NULL) ||
      item_size == 0 || count > SIZE_MAX / item_size ||
      count > (size_t)INT_MAX) {
    status = MVMC_KRYLOV_STATUS_INVALID_ARGUMENT;
  }
  return mvmc_power_lanczos_sampling_synchronize(topology, status);
}

MVMCKrylovStatus mvmc_power_lanczos_sampling_sum_complex(
    const MVMCPowerLanczosSamplingTopology *topology,
    double complex *values, size_t count) {
  MVMCKrylovStatus status =
      validate_reduction(topology, values, count, sizeof(*values));
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (count == 0) return MVMC_KRYLOV_STATUS_OK;
#ifdef _mpi_use
  if (topology->chain_rank != 0) {
    memset(values, 0, count * sizeof(*values));
  }
  if (MPI_Allreduce(MPI_IN_PLACE, values, (int)count,
                    MPI_C_DOUBLE_COMPLEX, MPI_SUM,
                    topology->world_communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#endif
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_sampling_sum_double(
    const MVMCPowerLanczosSamplingTopology *topology,
    double *values, size_t count) {
  MVMCKrylovStatus status =
      validate_reduction(topology, values, count, sizeof(*values));
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (count == 0) return MVMC_KRYLOV_STATUS_OK;
#ifdef _mpi_use
  if (topology->chain_rank != 0) {
    memset(values, 0, count * sizeof(*values));
  }
  if (MPI_Allreduce(MPI_IN_PLACE, values, (int)count, MPI_DOUBLE, MPI_SUM,
                    topology->world_communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#endif
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_sampling_sum_u64(
    const MVMCPowerLanczosSamplingTopology *topology,
    uint64_t *values, size_t count) {
  MVMCKrylovStatus status =
      validate_reduction(topology, values, count, sizeof(*values));
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (count == 0) return MVMC_KRYLOV_STATUS_OK;
#ifdef _mpi_use
  if (topology->chain_rank != 0) {
    memset(values, 0, count * sizeof(*values));
  }
  if (MPI_Allreduce(MPI_IN_PLACE, values, (int)count, MPI_UINT64_T, MPI_SUM,
                    topology->world_communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#endif
  return MVMC_KRYLOV_STATUS_OK;
}

MVMCKrylovStatus mvmc_power_lanczos_sampling_max_double(
    const MVMCPowerLanczosSamplingTopology *topology,
    double *values, size_t count) {
  MVMCKrylovStatus status =
      validate_reduction(topology, values, count, sizeof(*values));
  size_t index;
  if (status != MVMC_KRYLOV_STATUS_OK) return status;
  if (count == 0) return MVMC_KRYLOV_STATUS_OK;
#ifdef _mpi_use
  if (topology->chain_rank != 0) {
    for (index = 0; index < count; ++index) values[index] = -INFINITY;
  }
  if (MPI_Allreduce(MPI_IN_PLACE, values, (int)count, MPI_DOUBLE, MPI_MAX,
                    topology->world_communicator) != MPI_SUCCESS) {
    return MVMC_KRYLOV_STATUS_COLLECTIVE_FAILURE;
  }
#else
  (void)index;
#endif
  return MVMC_KRYLOV_STATUS_OK;
}

size_t mvmc_power_lanczos_sampling_block_index(
    const MVMCPowerLanczosSamplingTopology *topology,
    size_t local_sample, size_t block_count) {
  uint64_t global_sample;
  uint64_t samples_per_block;
  if (topology == NULL || block_count == 0 ||
      local_sample >= topology->samples_per_chain ||
      topology->total_samples % (uint64_t)block_count != 0) {
    return SIZE_MAX;
  }
  global_sample = topology->chain_index * topology->samples_per_chain +
                  (uint64_t)local_sample;
  samples_per_block = topology->total_samples / (uint64_t)block_count;
  return (size_t)(global_sample / samples_per_block);
}
