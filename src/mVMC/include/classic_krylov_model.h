/*
mVMC - A numerical solver package for a wide range of quantum lattice models
Copyright (C) 2016 The University of Tokyo, All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef MVMC_CLASSIC_KRYLOV_MODEL_H
#define MVMC_CLASSIC_KRYLOV_MODEL_H

#if defined(MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)

#include "classic_pfaffian_collective.h"
#include "krylov_fock_reference.h"

#include <complex.h>
#include <stddef.h>

typedef enum {
  MVMC_CLASSIC_KRYLOV_SOURCE_TRANSFER = 1,
  MVMC_CLASSIC_KRYLOV_SOURCE_COULOMB_INTRA,
  MVMC_CLASSIC_KRYLOV_SOURCE_COULOMB_INTER,
  MVMC_CLASSIC_KRYLOV_SOURCE_HUND,
  MVMC_CLASSIC_KRYLOV_SOURCE_EXCHANGE
} MVMCClassicKrylovSourceKind;

typedef struct {
  int output_site;
  int output_spin;
  int input_site;
  int input_spin;
  double complex coefficient;
} MVMCClassicKrylovTransfer;

typedef struct {
  int site;
  double coefficient;
} MVMCClassicKrylovSiteCoupling;

typedef struct {
  int first_site;
  int second_site;
  double coefficient;
} MVMCClassicKrylovPairCoupling;

/*
 * Counts remain signed so a root-global hook can reject malicious negative
 * production counts before any size_t conversion.  Raw pointers are read on
 * communicator rank 0 only; non-root callers may pass raw == NULL.
 */
typedef struct {
  int site_count;
  int up_electron_count;
  int down_electron_count;
  int pure_spin;
  int transfer_count;
  const MVMCClassicKrylovTransfer *transfers;
  int coulomb_intra_count;
  const MVMCClassicKrylovSiteCoupling *coulomb_intra;
  int coulomb_inter_count;
  const MVMCClassicKrylovPairCoupling *coulomb_inter;
  int hund_count;
  const MVMCClassicKrylovPairCoupling *hund;
  int exchange_count;
  const MVMCClassicKrylovPairCoupling *exchange;
  int pair_hopping_count;
  int inter_all_count;
  int nbody_inter_all_count;
} MVMCClassicKrylovRawModel;

typedef struct MVMCClassicKrylovModelWorkspace
    MVMCClassicKrylovModelWorkspace;

/*
 * Collective creation on communicator.  Rank 0 validates raw metadata and
 * every referenced element before broadcasting status and checked sizes.
 * The resulting immutable Fock model is byte-identical on every rank and
 * does not retain a pointer to the root raw data.
 */
MVMCKrylovStatus mvmc_classic_krylov_model_workspace_create_from_root(
    const MVMCClassicKrylovRawModel *raw,
    MVMCClassicPfaffianCommunicator communicator,
    MVMCClassicKrylovModelWorkspace **workspace);

void mvmc_classic_krylov_model_workspace_destroy(
    MVMCClassicKrylovModelWorkspace *workspace);

const MVMCKrylovFockModel *mvmc_classic_krylov_model(
    const MVMCClassicKrylovModelWorkspace *workspace);

size_t mvmc_classic_krylov_model_workspace_bytes(
    const MVMCClassicKrylovModelWorkspace *workspace);

#endif /* MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE */

#endif /* MVMC_CLASSIC_KRYLOV_MODEL_H */
