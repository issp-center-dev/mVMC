/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_POWER_LANCZOS_SAMPLING_H
#define MVMC_POWER_LANCZOS_SAMPLING_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "power_lanczos_classic_bridge.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

typedef struct {
  MVMCClassicPfaffianCommunicator world_communicator;
  MVMCClassicPfaffianCommunicator chain_communicator;
  int world_rank;
  int world_size;
  int chain_rank;
  int chain_size;
  uint64_t chain_index;
  uint64_t chain_count;
  uint64_t samples_per_chain;
  uint64_t total_samples;
} MVMCPowerLanczosSamplingTopology;

/*
 * Describe independent Markov chains embedded in the world communicator.
 * Rank zero of every chain communicator contributes that chain's statistics;
 * all ranks receive the same reduced result.
 */
MVMCKrylovStatus mvmc_power_lanczos_sampling_topology_create(
    MVMCClassicPfaffianCommunicator world_communicator,
    MVMCClassicPfaffianCommunicator chain_communicator,
    size_t samples_per_chain,
    MVMCPowerLanczosSamplingTopology *topology);

MVMCKrylovStatus mvmc_power_lanczos_sampling_synchronize(
    const MVMCPowerLanczosSamplingTopology *topology,
    MVMCKrylovStatus local_status);

/* In-place world reductions. Non-root chain ranks contribute zero. */
MVMCKrylovStatus mvmc_power_lanczos_sampling_sum_complex(
    const MVMCPowerLanczosSamplingTopology *topology,
    double complex *values, size_t count);
MVMCKrylovStatus mvmc_power_lanczos_sampling_sum_double(
    const MVMCPowerLanczosSamplingTopology *topology,
    double *values, size_t count);
MVMCKrylovStatus mvmc_power_lanczos_sampling_sum_u64(
    const MVMCPowerLanczosSamplingTopology *topology,
    uint64_t *values, size_t count);
MVMCKrylovStatus mvmc_power_lanczos_sampling_max_double(
    const MVMCPowerLanczosSamplingTopology *topology,
    double *values, size_t count);

size_t mvmc_power_lanczos_sampling_block_index(
    const MVMCPowerLanczosSamplingTopology *topology,
    size_t local_sample, size_t block_count);

#endif /* MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE */

#endif /* MVMC_POWER_LANCZOS_SAMPLING_H */
